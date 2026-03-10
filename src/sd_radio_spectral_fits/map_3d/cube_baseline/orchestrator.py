# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map.cube_baseline.orchestrator

One-iteration orchestration: line-free -> ripple freqs -> fit baseline -> update session.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from sd_radio_spectral_fits.map_3d.baseline_subtraction import (
    LineFreeConfig,
    RippleConfig,
    BaselineConfig,
    estimate_linefree_mask_from_cube,
    estimate_ripple_frequencies_fft,
)

from sd_radio_spectral_fits.utils import parse_windows, in_any_windows  # type: ignore

from .session import BaselineSession
from .engine import fit_cube_baseline


def create_manual_signal_mask_1d(
    v_axis: np.ndarray,
    v_windows: Union[List[str], List[Tuple[float, float]]],
) -> np.ndarray:
    """Manual velocity windows -> 1D boolean mask (True=signal)."""
    if not v_windows:
        return np.zeros(len(v_axis), dtype=bool)
    windows = parse_windows(v_windows) if isinstance(v_windows[0], str) else v_windows
    return in_any_windows(v_axis, windows)


def _aggregate_spectrum(cube: np.ndarray, *, max_pix: int = 200000, seed: int = 0) -> np.ndarray:
    """Spatial median spectrum for ripple frequency estimation."""
    nchan, ny, nx = cube.shape
    flat = cube.reshape(nchan, ny * nx)
    npix = ny * nx
    if npix > int(max_pix):
        rng = np.random.default_rng(int(seed))
        idx = rng.choice(npix, size=int(max_pix), replace=False)
        flat = flat[:, idx]
    with np.errstate(all="ignore"):
        spec = np.nanmedian(flat, axis=1)
    return np.asarray(spec, dtype=np.float32)


def _cube_for_scope(cube: np.ndarray, target_mask_2d: np.ndarray, scope: str) -> np.ndarray:
    """Return a cube view used for line-free/ripple estimation."""
    if scope == "global":
        return cube
    if scope != "target":
        raise ValueError(f"Unknown scope: {scope}")
    tm = np.asarray(target_mask_2d, dtype=bool)
    if tm.ndim != 2:
        raise ValueError(f"target_mask_2d must be 2D, got shape={tm.shape}")
    npix = int(tm.sum())
    if npix <= 0:
        raise ValueError("target_mask_2d selects zero pixels")
    return cube[:, tm][:, np.newaxis, :]


def _resolve_linefree_mask(
    session: BaselineSession,
    *,
    cube_work: np.ndarray,
    target_mask_2d: np.ndarray,
    auto_linefree: bool,
    linefree_cfg: LineFreeConfig,
    manual_v_windows: Optional[Union[List[str], List[Tuple[float, float]]]],
    linefree_mode: str,
    linefree_scope: str,
) -> np.ndarray:
    """
    Resolve the line-free mask used for the current fit.

    linefree_mode:
      - 'auto'    : estimate from cube_work
      - 'prior'   : use session prior, fallback to auto
      - 'current' : use current session state, fallback to prior, then auto
      - 'or'      : union(auto, prior/current) i.e. keep channels line-free if either source supports them
    """
    cube_detect = _cube_for_scope(cube_work, target_mask_2d, linefree_scope)

    lf_ref: Optional[np.ndarray] = None
    if linefree_mode == "current":
        lf_ref = session.linefree_mask_1d_current
        if lf_ref is None:
            lf_ref = session.linefree_mask_1d_prior
    elif linefree_mode == "prior":
        lf_ref = session.linefree_mask_1d_prior
    elif linefree_mode == "or":
        lf_ref = session.linefree_mask_1d_current if session.linefree_mask_1d_current is not None else session.linefree_mask_1d_prior
    elif linefree_mode != "auto":
        raise ValueError(f"Unknown linefree_mode: {linefree_mode}")

    need_auto = (linefree_mode in ("auto", "or")) or ((linefree_mode in ("prior", "current")) and (lf_ref is None))
    lf_auto: Optional[np.ndarray] = None
    if need_auto:
        if not auto_linefree:
            raise ValueError(
                f"linefree_mode={linefree_mode!r} requires automatic line-free estimation, "
                "but auto_linefree is disabled and no usable prior/current mask is available."
            )
        lf_auto = estimate_linefree_mask_from_cube(cube_detect, cfg=linefree_cfg, agg="median")

    if linefree_mode == "or":
        if lf_ref is None and lf_auto is None:
            lf = np.ones(session.nchan, dtype=bool)
        elif lf_ref is None:
            lf = np.asarray(lf_auto, dtype=bool)
        elif lf_auto is None:
            lf = np.asarray(lf_ref, dtype=bool)
        else:
            lf = np.asarray(lf_ref, dtype=bool) | np.asarray(lf_auto, dtype=bool)
    elif lf_ref is not None:
        lf = np.asarray(lf_ref, dtype=bool)
    elif lf_auto is not None:
        lf = np.asarray(lf_auto, dtype=bool)
    else:
        lf = np.ones(session.nchan, dtype=bool)

    if manual_v_windows:
        manual_signal = create_manual_signal_mask_1d(session.v_axis, manual_v_windows)
        lf = lf & (~manual_signal)

    return np.asarray(lf, dtype=bool)


def _resolve_ripple_freqs(
    session: BaselineSession,
    *,
    cube_work: np.ndarray,
    target_mask_2d: np.ndarray,
    lf: np.ndarray,
    enable_ripple: bool,
    ripple_cfg: RippleConfig,
    ripple_freqs: Optional[Sequence[float]],
    ripple_mode: str,
    ripple_scope: str,
    poly_order: int,
) -> Optional[Sequence[float]]:
    """Resolve ripple frequencies for the current fit."""
    if not enable_ripple:
        return None
    if ripple_freqs is not None:
        return [float(f) for f in ripple_freqs]

    ref = None
    if ripple_mode == "current":
        ref = session.ripple_freqs_current
        if ref is None:
            ref = session.ripple_freqs_prior
    elif ripple_mode == "prior":
        ref = session.ripple_freqs_prior
    elif ripple_mode != "auto":
        raise ValueError(f"Unknown ripple_mode: {ripple_mode}")

    if ref is not None and len(ref) > 0:
        return list(np.asarray(ref, dtype=float))

    cube_detect = _cube_for_scope(cube_work, target_mask_2d, ripple_scope)
    spec = _aggregate_spectrum(cube_detect, max_pix=200000, seed=0)
    freqs = estimate_ripple_frequencies_fft(spec, lf, rcfg=ripple_cfg, poly_order_pre=poly_order)
    if not freqs:
        return None
    return freqs


def run_one_iteration(
    session: BaselineSession,
    *,
    target_mask_2d: Optional[np.ndarray] = None,
    # line-free
    auto_linefree: bool = True,
    linefree_cfg: LineFreeConfig = LineFreeConfig(),
    manual_v_windows: Optional[Union[List[str], List[Tuple[float, float]]]] = None,
    linefree_mode: str = "auto",
    linefree_scope: str = "global",
    # ripple
    enable_ripple: bool = True,
    ripple_cfg: RippleConfig = RippleConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    ripple_mode: str = "auto",
    ripple_scope: str = "global",
    # model + compute
    poly_order: int = 1,
    robust: bool = False,
    chunk_pix: int = 65536,
) -> None:
    """Execute one iteration and update `session` in-place."""
    if target_mask_2d is None:
        target_mask_2d = np.ones((session.ny, session.nx), dtype=bool)
    else:
        target_mask_2d = np.asarray(target_mask_2d, dtype=bool)
        if target_mask_2d.shape != (session.ny, session.nx):
            raise ValueError(f"target_mask_2d shape mismatch: {target_mask_2d.shape} vs {(session.ny, session.nx)}")

    cube_work = session.get_full_cube_work()

    lf = _resolve_linefree_mask(
        session,
        cube_work=cube_work,
        target_mask_2d=target_mask_2d,
        auto_linefree=auto_linefree,
        linefree_cfg=linefree_cfg,
        manual_v_windows=manual_v_windows,
        linefree_mode=linefree_mode,
        linefree_scope=linefree_scope,
    )

    freqs = _resolve_ripple_freqs(
        session,
        cube_work=cube_work,
        target_mask_2d=target_mask_2d,
        lf=lf,
        enable_ripple=enable_ripple,
        ripple_cfg=ripple_cfg,
        ripple_freqs=ripple_freqs,
        ripple_mode=ripple_mode,
        ripple_scope=ripple_scope,
        poly_order=poly_order,
    )

    bcfg = BaselineConfig(poly_order=poly_order, ripple=enable_ripple, robust=robust, chunk_pix=int(chunk_pix))
    baseline_cube, stats = fit_cube_baseline(
        session.cube_original,
        linefree_mask_1d=lf,
        baseline_cfg=bcfg,
        ripple_freqs=freqs,
        target_mask_2d=target_mask_2d,
    )

    params: Dict[str, Any] = dict(
        poly_order=int(poly_order),
        robust=bool(robust),
        chunk_pix=int(chunk_pix),
        auto_linefree=bool(auto_linefree),
        linefree_cfg=linefree_cfg.__dict__,
        manual_v_windows=manual_v_windows,
        linefree_mode=str(linefree_mode),
        linefree_scope=str(linefree_scope),
        enable_ripple=bool(enable_ripple),
        ripple_cfg=ripple_cfg.__dict__,
        ripple_mode=str(ripple_mode),
        ripple_scope=str(ripple_scope),
        ripple_freqs=list(freqs) if freqs else None,
    )
    session.update_baseline(
        baseline_cube,
        target_mask_2d=target_mask_2d,
        stats=stats,
        params=params,
        linefree_mask_1d=lf,
        ripple_freqs=freqs,
    )
