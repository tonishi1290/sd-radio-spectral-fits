# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map_3d.baseline_subtraction

Baseline subtraction utilities for spectral cubes and OTF bundles.

This module intentionally separates three layers that are easy to confuse in
user-facing APIs.

1) Base line-free mask selection
   - ``linefree_mask``               : exact manual mask supplied by the caller
   - ``linefree_mode='manual'``      : exact manual mask built from
     ``linefree_velocity_windows_kms``
   - ``linefree_mode='prior'``       : reuse pre-existing ``LINEFREE`` from the
     input bundle / FITS when available
   - ``linefree_mode='auto'``        : run automatic line-free detection from a
     ``LineFreeConfig``
   - ``linefree_mode='or'``          : union of ``prior`` and ``auto``
   - ``linefree_mode='infer'``       : resolve automatically from the available
     inputs without silently enabling auto detection when nothing was requested

2) Post-selection velocity-window modifiers
   - ``linefree_velocity_windows_kms`` :
       * in ``manual`` mode, this is the exact line-free definition
       * in ``auto`` / ``prior`` / ``or`` / ``infer`` after resolution, this is
         an OR-style include modifier applied after the base mask is built
   - ``exclude_v_windows`` : final exclusion windows removed from the fit mask

3) Ripple handling
   - ``BaselineConfig.ripple`` : whether the *final fitted baseline model*
     includes sinusoidal ripple terms
   - ``ripple_freqs`` / ``ripple_mode`` / ``ripple_cfg`` : where ripple
     frequencies come from
   - ``ripple_apply_stage`` / ``safe_velocity_windows_kms`` : whether ripple is
     used only for the detection cube, only for the final fit, or both

Array convention
----------------
All cube-like arrays are handled as ``(nchan, ny, nx)``.

Design goal
-----------
The orchestration functions in this module are intentionally conservative:
manual inputs must stay exact when requested, auto detection should only run
when explicitly enabled or inferable from explicit auto inputs, and QC products
should preserve the mask that was actually used.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import logging
import os
import platform
import sys
import time
import warnings
from dataclasses import dataclass, replace
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
import scipy
import scipy.ndimage as ndi
import astropy
from astropy.io import fits
import astropy.units as u
from astropy.table import Table

try:
    import psutil  # type: ignore
except Exception:  # pragma: no cover
    psutil = None  # type: ignore

try:
    from spectral_cube import SpectralCube
except Exception:  # pragma: no cover
    SpectralCube = None  # type: ignore

from .otf_bundle import OTFBundle
from .otf_bundle_io import read_otf_bundle, write_otf_bundle, validate_otf_bundle
from .mosaic import attach_mosaic_products_from_mask


def _profile_rss_gb() -> Optional[float]:
    if psutil is None:
        return None
    try:
        return float(psutil.Process(os.getpid()).memory_info().rss) / 1.0e9
    except Exception:
        return None


def _profile_record(store: List[dict[str, Any]], stage: str, t_prev: float, t0: float, **extra: Any) -> float:
    now = time.perf_counter()
    rec: dict[str, Any] = {
        "stage": str(stage),
        "dt_stage_sec": float(now - t_prev),
        "dt_total_sec": float(now - t0),
    }
    rss = _profile_rss_gb()
    if rss is not None:
        rec["rss_gb"] = float(rss)
    rec.update(extra)
    store.append(rec)
    return now


def _write_profile_sidecar(profile_path_prefix: str, records: List[dict[str, Any]]) -> tuple[str, str]:
    json_path = profile_path_prefix + ".baseline_profile.json"
    txt_path = profile_path_prefix + ".baseline_profile.txt"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump({"records": records}, f, ensure_ascii=False, indent=2)
    lines = ["# baseline profile"]
    for rec in records:
        rss = rec.get("rss_gb")
        rss_txt = f"  rss_gb={rss:.6f}" if isinstance(rss, (int, float)) else ""
        lines.append(f"{rec['stage']}: dt_stage_sec={rec['dt_stage_sec']:.6f}  dt_total_sec={rec['dt_total_sec']:.6f}{rss_txt}")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return json_path, txt_path

__all__ = [
    "LineFreeConfig",
    "RippleConfig",
    "BaselineConfig",
    "estimate_linefree_mask_1d",
    "estimate_linefree_mask_from_cube",
    "estimate_linefree_mask_from_cube_3d",
    "estimate_ripple_frequencies_fft",
    "subtract_baseline_cube",
    "subtract_baseline_from_bundle",
    "subtract_baseline_from_fits",
    "build_baseline_diagnostic_payload",
    "write_baseline_diagnostic_files",
]


# -----------------------------------------------------------------------------
# Configs
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class LineFreeConfig:
    """
    Automatic line-free detection configuration.

    Ownership
    ---------
    ``LineFreeConfig`` belongs to the *auto-detection* layer only.  It does not
    decide whether auto detection is used; that decision belongs to the main API
    parameters ``linefree_mode`` and ``linefree_mask``.

    This config is consulted only when the effective base line-free source is
    ``auto`` or the auto part of ``or``.

    Core role
    ---------
    ``LineFreeConfig`` controls how a provisional line-free mask is estimated
    from the cube itself.  In the global-1D path this is based on aggregated
    spectra.  In the voxel-wise 3D path it controls smoothing, seed generation,
    thresholds, hysteresis, dilation, and related parameters used before the
    final baseline fit.

    Parameters
    ----------
    smooth_width : int
        Median-filter width (channels) used to estimate a slow trend before
        sigma-clipping in the global-1D auto path.
    sigma : float
        Sigma threshold for clipping residual excursions.
    iters : int
        Number of sigma-clipping iterations.
    pad_chan : int
        Spectral dilation applied to detected line channels.
    min_linefree_frac : float
        Minimum acceptable line-free fraction; smaller fractions indicate that
        the resulting fit may be ill-conditioned.
    sampling_mode, sampling_seed : see field names
        Spatial subsampling controls for aggregated auto detection.
    interpolate_nans : bool
        Whether NaN gaps in aggregated spectra are interpolated before auto
        detection.
    conservative_quantile, conservative_both_tails :
        Optional conservative detection aids for sparse emission/absorption.

    3D / voxel-wise auto-mask parameters
    ------------------------------------
    The ``lf3d_*`` fields control the voxel-wise auto path.  Public aliases use
    arcsec / km/s where possible, while the internal solver still works in pixel
    / channel units.  These fields are ignored unless the chosen auto method
    actually enters the 3D detection path.

    ``lf3d_min_run_kms`` is a post-detection cleanup on the provisional 3D
    signal mask: for each spatial pixel, signal runs shorter than the requested
    contiguous velocity width are returned to the baseline side.  ``None`` or a
    non-positive value disables this cleanup.
    """
    smooth_width: int = 51
    sigma: float = 4.0
    iters: int = 6
    pad_chan: int = 3
    min_linefree_frac: float = 0.35
    sampling_mode: str = "random"
    sampling_seed: int = 0
    interpolate_nans: bool = True
    conservative_quantile: Optional[float] = None
    conservative_both_tails: bool = True
    mask_kind: str = "global_1d"
    auto_method: str = "agg_1d"
    lf3d_max_iter: int = 3
    lf3d_threshold: float = 3.0
    lf3d_dilation: int = 4
    lf3d_chunk_pix: int = 4096
    lf3d_positive_only: bool = True
    lf3d_reg: float = 1.0e-7
    lf3d_solver: str = "qr"
    lf3d_sigma_mode: str = "negative"
    lf3d_monotone_line: bool = True
    lf3d_hysteresis_sigma: Optional[float] = 1.5
    lf3d_weak_spatial_sigma: float = 1.25
    lf3d_seed_sigma: float = 4.0
    lf3d_seed_spatial_sigma: float = 1.5
    # Detection-cube smoothing applied before Stage-A local seed generation.
    # Internal units remain pixel / channel widths; public aliases expose arcsec / km/s.
    # The spectral width default of 3 means the public API applies a *mild*
    # detection-only smoothing along the spectral axis unless the user overrides it.
    lf3d_detect_spatial_sigma: float = 0.0
    lf3d_detect_spectral_width: int = 3
    lf3d_min_run_kms: Optional[float] = None
    lf3d_seed_baseline_order: int = 1
    lf3d_seed_dilation: int = 0
    lf3d_seed_hysteresis_sigma: Optional[float] = None
    lf3d_seed_baseline_mode: str = "median"
    lf3d_seed_median_width: int = 41
    lf3d_prov_poly_order: int = 1
    # public-facing aliases (preferred in docs / user configs)
    lf3d_seed_spatial_fwhm_arcsec: Optional[float] = None
    lf3d_seed_spatial_fwhm_pix: Optional[float] = None
    lf3d_seed_median_width_kms: Optional[float] = None
    lf3d_seed_median_width_chan: Optional[float] = None
    lf3d_seed_threshold_sigma: Optional[float] = None
    lf3d_seed_hysteresis_threshold_sigma: Optional[float] = None
    lf3d_seed_dilation_kms: Optional[float] = None
    lf3d_seed_dilation_chan: Optional[float] = None
    lf3d_detect_threshold_sigma: Optional[float] = None
    lf3d_hysteresis_threshold_sigma: Optional[float] = None
    lf3d_detect_spatial_fwhm_arcsec: Optional[float] = None
    lf3d_detect_spatial_fwhm_pix: Optional[float] = None
    lf3d_detect_spectral_width_kms: Optional[float] = None
    lf3d_detect_spectral_width_chan: Optional[float] = None
    lf3d_weak_spatial_fwhm_arcsec: Optional[float] = None
    lf3d_weak_spatial_fwhm_pix: Optional[float] = None
    # legacy / compatibility knobs kept for older configs; not used by r9 local-seed core
    lf3d_degree: int = 1
    lf3d_spatial_sigma: float = 1.25
    lf3d_seed_quantile: Optional[float] = None
    lf3d_seed_max_pix: int = 200000
    lf3d_asym_tau: Optional[float] = None
    lf3d_asym_iters: int = 1
    lf3d_asym_floor: float = 0.05
    lf3d_seed_weight_floor: float = 0.15
    auto_mask_profile: str = "default"




def _expand_linefree_profile(cfg: LineFreeConfig) -> LineFreeConfig:
    """Expand optional high-level presets for voxel-wise auto line-free masking.

    Supported profiles
    ------------------
    default
        Keep the configuration as-is.
    no_ripple
        Use a more aggressive Stage-A/B tuning intended for data where strong
        standing-wave ripple is absent or already mitigated upstream. This only
        affects the voxel-wise 3D auto-mask path; final Stage-C fitting and the
        global-1D path are unchanged.
    """
    profile = str(getattr(cfg, 'auto_mask_profile', 'default') or 'default').strip().lower()
    if profile in {'', 'default', 'standard', 'none'}:
        return cfg
    if profile == 'no_ripple':
        return replace(
            cfg,
            lf3d_seed_sigma=3.0,
            lf3d_seed_hysteresis_sigma=1.25,
            lf3d_seed_median_width=31,
            lf3d_seed_spatial_sigma=1.25,
            lf3d_threshold=2.5,
            lf3d_hysteresis_sigma=1.0,
            lf3d_weak_spatial_sigma=1.5,
        )
    raise ValueError(
        f"Unknown auto_mask_profile={profile!r}; expected 'default' or 'no_ripple'."
    )


@dataclass(frozen=True)
class RippleConfig:
    """
    Ripple-frequency search configuration.

    Ownership
    ---------
    ``RippleConfig`` belongs to the *frequency-search* layer.  It does not by
    itself enable ripple terms in the final model and it does not decide which
    stage uses ripple.  Those decisions belong to ``BaselineConfig.ripple`` and
    the main-API parameter ``ripple_apply_stage``.

    Core role
    ---------
    This config controls how ripple frequencies are estimated when they are not
    supplied explicitly via ``ripple_freqs`` and not taken from a prior
    ``RIPFREQ`` table.  Internal fitting uses cycles/channel, while the search
    range may be expressed in channels, Hz, or km/s.

    Parameters
    ----------
    nfreq : int
        Number of dominant ripple frequencies to keep.
    period_range_chan : (float, float)
        Search range in channel periods.
    period_range_hz, period_range_kms : tuple or None
        Optional user-facing search ranges.  If provided, they override
        ``period_range_chan`` for the search stage only.
    min_separation : float
        Minimum separation between selected FFT peaks in cycles/channel.
    window : {'hann', 'none'}
        Window function applied before FFT.
    interpolate_masked : bool
        Whether masked channels are interpolated before FFT instead of zero
        filling.
    stable_sort : bool
        Whether FFT peaks are ranked with a stable sort.
    sampling_mode, sampling_seed : see field names
        Aggregation-side sampling controls used by higher-level orchestration.
    """
    nfreq: int = 2
    period_range_chan: Tuple[float, float] = (20.0, 400.0)
    period_range_hz: Optional[Tuple[float, float]] = None
    period_range_kms: Optional[Tuple[float, float]] = None
    min_separation: float = 0.002
    window: str = "hann"
    interpolate_masked: bool = False
    stable_sort: bool = True
    sampling_mode: str = "random"
    sampling_seed: int = 0


@dataclass(frozen=True)
class BaselineConfig:
    """
    Final baseline-model configuration.

    Ownership
    ---------
    ``BaselineConfig`` belongs to the *final fitting model* layer.

    It answers questions such as:
    - What polynomial order is fitted?
    - Are sinusoidal ripple terms included in the final fitted model?
    - Is robust refitting enabled?
    - Which numerical solver / dtype / chunking strategy is used?

    It does **not** decide:
    - where the line-free mask comes from
    - where ripple frequencies come from
    - whether ripple is used in the detection-only stage

    Those choices belong to the main API parameters and to ``RippleConfig``.

    Important note about ``ripple``
    -------------------------------
    ``BaselineConfig.ripple`` is the switch for including sinusoidal ripple
    terms in the *final baseline model*.  In the current implementation it also
    gates some detection-stage helper paths when ``ripple_apply_stage`` is
    ``'mask_only'`` or ``'both'``.  This is more powerful than the field name
    may suggest, so callers should treat it as the master ripple enable switch.

    Parameters
    ----------
    poly_order : int
        Polynomial order (0=constant, 1=linear, ...).
    ripple : bool
        Whether ripple sin/cos terms are included in the final model when
        ripple frequencies are available.
    robust, robust_mode, robust_iters, robust_early_stop,
    robust_coef_rtol, robust_coef_atol,
    robust_selective_sigma, robust_selective_frac,
    robust_selective_max_pixels, robust_batch_pixels :
        Controls for iterative robust refitting.
    rcond : float or None
        Least-squares cut-off passed to the linear solver.
    chunk_pix : int
        Number of spatial pixels processed per chunk.
    reproducible_mode : bool
        Deterministic / version-stable numerical choices where practical.
    compute_dtype : {'float32', 'float64'}
        Internal compute dtype.
    normalize_x : bool
        Whether the polynomial coordinate is normalized to [-1, 1].  Ripple
        frequencies remain expressed in cycles/channel regardless of this flag.
    strict_failures : bool
        Raise immediately on fitting failures.
    fallback_to_pixelwise : bool
        Whether a failed vectorized solve may fall back to pixel-wise fitting.
    voxel_solver, voxel_qr_batch_pix, voxel_grouping, voxel_group_min_size :
        Numerical controls for voxel-wise 3D-mask fitting paths.
    """
    poly_order: int = 1
    ripple: bool = True
    robust: bool = False
    robust_mode: str = "full"
    robust_iters: int = 3
    robust_early_stop: bool = False
    robust_coef_rtol: float = 1.0e-4
    robust_coef_atol: float = 1.0e-6
    robust_selective_sigma: float = 4.5
    robust_selective_frac: float = 0.10
    robust_selective_max_pixels: int = 2048
    robust_batch_pixels: int = 512
    rcond: Optional[float] = None
    chunk_pix: int = 65536
    reproducible_mode: bool = False
    compute_dtype: str = "float32"
    normalize_x: bool = False
    strict_failures: bool = False
    fallback_to_pixelwise: bool = True
    voxel_solver: str = "qr"
    voxel_qr_batch_pix: int = 2048
    voxel_grouping: str = "auto"
    voxel_group_min_size: int = 12


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------


def _robust_std(x: np.ndarray) -> float:
    """Robust std via MAD * 1.4826."""
    x = np.asarray(x, float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    s = 1.4826 * mad
    if not np.isfinite(s) or s <= 0:
        s = float(np.nanstd(x))
    if not np.isfinite(s) or s <= 0:
        s = 0.0
    return float(s)


def _median_filter_1d(y: np.ndarray, width: int) -> np.ndarray:
    w = int(width)
    if w < 3:
        return y
    if w % 2 == 0:
        w += 1
    return ndi.median_filter(y, size=w, mode="nearest")


def _dilate_1d(mask: np.ndarray, pad: int) -> np.ndarray:
    p = int(pad)
    if p <= 0:
        return mask
    struct = np.ones((2 * p + 1,), dtype=bool)
    return ndi.binary_dilation(mask, structure=struct)


def _as_float32(x: Any) -> np.ndarray:
    return np.asarray(x, dtype=np.float32)


def _capture_stdout(fn: Any) -> str:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            fn()
        except Exception as exc:
            return f"<unavailable: {exc}>"
    return buf.getvalue().strip()


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_array(arr: Any) -> str:
    a = np.ascontiguousarray(np.asarray(arr))
    return _sha256_bytes(a.tobytes())


def _json_ready(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _json_ready(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_ready(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _json_ready(obj.tolist())
    if isinstance(obj, (np.floating,)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def build_baseline_diagnostic_payload(
    *,
    linefree_mask: np.ndarray,
    ripple_freqs: Optional[Sequence[float]],
    rms_map: Optional[np.ndarray],
    flag_map: Optional[np.ndarray],
    extra: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    lf = np.asarray(linefree_mask, dtype=bool)
    rip = np.asarray([] if ripple_freqs is None else list(ripple_freqs), dtype=np.float64)

    payload: dict[str, Any] = {
        "environment": {
            "python_version": sys.version.replace("\n", " "),
            "platform": platform.platform(),
            "numpy_version": np.__version__,
            "scipy_version": getattr(scipy, "__version__", None),
            "astropy_version": getattr(astropy, "__version__", None),
            "numpy_show_runtime": _capture_stdout(getattr(np, "show_runtime", lambda: None)),
            "numpy_show_config": _capture_stdout(getattr(np, "show_config", lambda: None)),
        },
        "linefree_used": {
            "shape": [int(v) for v in lf.shape],
            "ndim": int(lf.ndim),
            "nchan": int(lf.shape[0]) if lf.ndim >= 1 else int(lf.size),
            "n_true": int(lf.sum()),
            "fraction_true": float(lf.mean()) if lf.size > 0 else None,
            "sha256": _sha256_array(lf.astype(np.uint8, copy=False)),
        },
        "ripple_used": {
            "nfreq": int(rip.size),
            "freq_cyc_per_ch": [float(v) for v in rip.tolist()],
            "period_ch": [float(np.inf if v == 0.0 else 1.0 / v) for v in rip.tolist()],
            "sha256": _sha256_array(rip.astype(np.float64, copy=False)),
        },
    }

    if rms_map is not None:
        rms = np.asarray(rms_map, dtype=np.float32)
        finite = np.isfinite(rms)
        if finite.any():
            with np.errstate(all="ignore"):
                rms_min = float(np.nanmin(rms))
                rms_med = float(np.nanmedian(rms))
                rms_mean = float(np.nanmean(rms))
                rms_max = float(np.nanmax(rms))
        else:
            rms_min = rms_med = rms_mean = rms_max = None
        payload["base_rms"] = {
            "shape": [int(v) for v in rms.shape],
            "sha256": _sha256_array(rms),
            "min": rms_min,
            "median": rms_med,
            "mean": rms_mean,
            "max": rms_max,
        }
    else:
        payload["base_rms"] = None

    if flag_map is not None:
        flg = np.asarray(flag_map, dtype=np.uint8)
        vals, cnts = np.unique(flg, return_counts=True)
        payload["base_flag"] = {
            "shape": [int(v) for v in flg.shape],
            "sha256": _sha256_array(flg),
            "counts": {str(int(v)): int(c) for v, c in zip(vals.tolist(), cnts.tolist())},
        }
    else:
        payload["base_flag"] = None

    if extra:
        payload["extra"] = _json_ready(extra)
    return payload


def _format_baseline_diagnostic_text(payload: dict[str, Any]) -> str:
    env = payload.get("environment", {})
    lf = payload.get("linefree_used", {})
    rip = payload.get("ripple_used", {})
    rms = payload.get("base_rms", {}) or {}
    flg = payload.get("base_flag", {}) or {}
    lines = [
        "Baseline diagnostic summary",
        "===========================",
        f"Python : {env.get('python_version')}",
        f"Platform: {env.get('platform')}",
        f"NumPy  : {env.get('numpy_version')}",
        f"SciPy  : {env.get('scipy_version')}",
        f"Astropy: {env.get('astropy_version')}",
        "",
        "Line-free mask used",
        "-------------------",
        f"nchan        : {lf.get('nchan')}",
        f"n_true       : {lf.get('n_true')}",
        f"fraction_true: {lf.get('fraction_true')}",
        f"sha256       : {lf.get('sha256')}",
        "",
        "Ripple frequencies used",
        "-----------------------",
        f"nfreq        : {rip.get('nfreq')}",
        f"freq_cyc/ch  : {rip.get('freq_cyc_per_ch')}",
        f"period_ch    : {rip.get('period_ch')}",
        f"sha256       : {rip.get('sha256')}",
        "",
        "BASE_RMS summary",
        "----------------",
        f"min/median/mean/max: {rms.get('min')} / {rms.get('median')} / {rms.get('mean')} / {rms.get('max')}",
        f"sha256            : {rms.get('sha256')}",
        "",
        "BASE_FLG summary",
        "----------------",
        f"counts       : {flg.get('counts')}",
        f"sha256       : {flg.get('sha256')}",
    ]
    extra = payload.get("extra")
    if extra is not None:
        lines.extend(["", "Extra", "-----", json.dumps(extra, ensure_ascii=False, indent=2, sort_keys=True)])
    runtime = env.get("numpy_show_runtime")
    if runtime:
        lines.extend(["", "NumPy runtime", "-------------", str(runtime)])
    return "\n".join(lines) + "\n"


def write_baseline_diagnostic_files(prefix: str, payload: dict[str, Any]) -> Tuple[str, str]:
    base = str(prefix)
    if base.lower().endswith('.fits'):
        base = base[:-5]
    json_path = base + '.baseline_diag.json'
    txt_path = base + '.baseline_diag.txt'
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(_json_ready(payload), f, ensure_ascii=False, indent=2, sort_keys=True)
        f.write('\n')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write(_format_baseline_diagnostic_text(payload))
    return json_path, txt_path


def _validate_sampling_mode(mode: str) -> str:
    mode = str(mode).lower()
    if mode not in {"random", "stride"}:
        raise ValueError(f"Unknown sampling_mode={mode!r}; expected 'random' or 'stride'.")
    return mode


def _select_sample_indices(npix: int, max_pix: int, *, mode: str, seed: int) -> Optional[np.ndarray]:
    max_pix = int(max_pix)
    if max_pix <= 0:
        raise ValueError(f"max_pix must be positive, got {max_pix}")
    if npix <= max_pix:
        return None

    mode = _validate_sampling_mode(mode)
    if mode == "random":
        rng = np.random.default_rng(int(seed))
        idx = np.sort(rng.choice(npix, size=max_pix, replace=False).astype(np.int64, copy=False))
        return idx

    # Deterministic evenly spaced subsampling.
    step = float(npix) / float(max_pix)
    idx = np.floor(np.arange(max_pix, dtype=np.float64) * step).astype(np.int64)
    idx = np.clip(idx, 0, npix - 1)
    # floor() with max_pix < npix is monotonic non-decreasing; enforce uniqueness just in case.
    _, uniq = np.unique(idx, return_index=True)
    idx = idx[np.sort(uniq)]
    if idx.size < max_pix:
        # pad deterministically from the tail if uniqueness removal dropped elements.
        need = max_pix - idx.size
        tail = np.setdiff1d(np.arange(npix - need, npix, dtype=np.int64), idx, assume_unique=False)
        idx = np.concatenate([idx, tail[:need]])
    return np.sort(idx[:max_pix])


def _interp_nans_1d(y: np.ndarray, *, dtype: np.dtype = np.float64) -> np.ndarray:
    arr = np.asarray(y, dtype=dtype)
    if arr.ndim != 1:
        raise ValueError(f"Expected 1D input, got shape={arr.shape}")
    finite = np.isfinite(arr)
    if finite.all():
        return arr
    if not finite.any():
        return np.zeros_like(arr, dtype=dtype)
    x = np.arange(arr.size, dtype=np.float64)
    out = arr.copy()
    out[~finite] = np.interp(x[~finite], x[finite], out[finite])
    return out


def _fill_nan_with_median_along_spec(data: np.ndarray) -> np.ndarray:
    """Fill NaNs in (nchan, npix) or (nchan, ny, nx) with per-pixel spectral median."""
    arr = np.asarray(data, dtype=np.float32)
    if not np.isnan(arr).any():
        return arr
    if arr.ndim == 2:
        med = np.nanmedian(arr, axis=0, keepdims=True)
        med = np.nan_to_num(med, nan=0.0).astype(np.float32, copy=False)
        return np.where(np.isfinite(arr), arr, med)
    if arr.ndim == 3:
        med = np.nanmedian(arr, axis=0, keepdims=True)
        med = np.nan_to_num(med, nan=0.0).astype(np.float32, copy=False)
        return np.where(np.isfinite(arr), arr, med)
    raise ValueError("data must be 2D or 3D")


def _normalize_linefree_mask_shape(mask: np.ndarray, *, nchan: int, ny: int, nx: int, name: str = "linefree_mask") -> np.ndarray:
    arr = np.asarray(mask, dtype=bool)
    if arr.shape == (int(nchan),):
        return arr
    if arr.shape == (int(nchan), int(ny), int(nx)):
        return arr
    raise ValueError(
        f"{name} shape mismatch: {arr.shape} vs ({nchan},) or ({nchan},{ny},{nx})"
    )


def _broadcast_channel_mask_to_cube(mask_1d: np.ndarray, ny: int, nx: int) -> np.ndarray:
    m = np.asarray(mask_1d, dtype=bool).reshape(-1)
    return np.broadcast_to(m[:, None, None], (m.shape[0], int(ny), int(nx))).copy()


def _merge_linefree_masks_or(
    mask_a: Optional[np.ndarray],
    mask_b: Optional[np.ndarray],
    *,
    nchan: int,
    ny: int,
    nx: int,
) -> np.ndarray:
    if mask_a is None and mask_b is None:
        return np.ones((int(nchan),), dtype=bool)
    if mask_a is None:
        return _normalize_linefree_mask_shape(np.asarray(mask_b, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="linefree_mask")
    if mask_b is None:
        return _normalize_linefree_mask_shape(np.asarray(mask_a, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="linefree_mask")
    a = _normalize_linefree_mask_shape(mask_a, nchan=nchan, ny=ny, nx=nx, name="linefree_mask")
    b = _normalize_linefree_mask_shape(mask_b, nchan=nchan, ny=ny, nx=nx, name="linefree_mask")
    if a.ndim == 1 and b.ndim == 1:
        return a | b
    if a.ndim == 1:
        a = _broadcast_channel_mask_to_cube(a, ny, nx)
    if b.ndim == 1:
        b = _broadcast_channel_mask_to_cube(b, ny, nx)
    return np.asarray(a, dtype=bool) | np.asarray(b, dtype=bool)


def _collapse_linefree_mask_for_ripple(linefree_mask: np.ndarray, *, threshold_frac: float = 0.5) -> np.ndarray:
    lf = np.asarray(linefree_mask, dtype=bool)
    if lf.ndim == 1:
        return lf.reshape(-1)
    if lf.ndim != 3:
        raise ValueError(f"linefree_mask must be 1D or 3D, got shape={lf.shape}")
    frac = np.mean(lf.astype(np.float32), axis=(1, 2))
    out = frac >= float(threshold_frac)
    if not np.any(out):
        out = np.any(lf, axis=(1, 2))
    return np.asarray(out, dtype=bool)


def _apply_linefree_include_exclude(
    linefree_mask: np.ndarray,
    *,
    velocity_axis_kms: Optional[np.ndarray],
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]],
    exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]],
) -> np.ndarray:
    lf = np.asarray(linefree_mask, dtype=bool)
    if linefree_velocity_windows_kms:
        if velocity_axis_kms is None:
            raise ValueError("Velocity axis is required when linefree_velocity_windows_kms is provided.")
        include_1d = _linefree_mask_from_axis(np.asarray(velocity_axis_kms, dtype=float), linefree_velocity_windows_kms)
        if lf.ndim == 1:
            lf = lf | include_1d
        elif lf.ndim == 3:
            lf = lf | include_1d[:, None, None]
        else:
            raise ValueError(f"linefree_mask must be 1D or 3D, got shape={lf.shape}")
    if exclude_v_windows:
        if velocity_axis_kms is None:
            raise ValueError("Velocity axis is required when exclude_v_windows is provided.")
        exclude_1d = _manual_signal_mask_from_axis(np.asarray(velocity_axis_kms, dtype=float), exclude_v_windows)
        if lf.ndim == 1:
            lf = lf & (~exclude_1d)
        elif lf.ndim == 3:
            lf = lf & (~exclude_1d[:, None, None])
        else:
            raise ValueError(f"linefree_mask must be 1D or 3D, got shape={lf.shape}")
    return np.asarray(lf, dtype=bool)


def _gaussian_filter_spatial_nan_normalized(cube: np.ndarray, sigma: float) -> np.ndarray:
    arr = np.asarray(cube, dtype=np.float32)
    s = float(sigma)
    if s <= 0.0:
        return arr
    valid = np.isfinite(arr).astype(np.float32)
    safe = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32, copy=False)
    num = ndi.gaussian_filter(safe, sigma=(0.0, s, s), mode="nearest")
    den = ndi.gaussian_filter(valid, sigma=(0.0, s, s), mode="nearest")
    out = np.full(arr.shape, np.nan, dtype=np.float32)
    np.divide(num, den, out=out, where=den > 0.0)
    return out


def _uniform_filter_spectral_nan_normalized(cube3d: np.ndarray, width: int) -> np.ndarray:
    """Apply NaN-aware boxcar smoothing along the spectral axis only."""
    arr = np.asarray(cube3d, dtype=np.float32)
    if arr.ndim != 3:
        raise ValueError(f"cube3d must be 3D, got shape={arr.shape}")
    w = int(width)
    if w <= 1:
        return arr.astype(np.float32, copy=True)
    if w % 2 == 0:
        w += 1
    valid = np.isfinite(arr).astype(np.float32)
    safe = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32, copy=False)
    num = ndi.uniform_filter1d(safe, size=w, axis=0, mode="nearest")
    den = ndi.uniform_filter1d(valid, size=w, axis=0, mode="nearest")
    out = np.full(arr.shape, np.nan, dtype=np.float32)
    np.divide(num, den, out=out, where=den > 0.0)
    return out


def _make_lf3d_detection_cube(cube_data: np.ndarray, cfg: LineFreeConfig) -> np.ndarray:
    """Build the Stage-A detection cube used only for local seed detection.

    Detection smoothing is intentionally separated from the final fit: Stage-A may
    smooth in (x,y) and mildly along the spectral axis to stabilize seed finding,
    while Stage-B / Stage-C continue to fit the original unsmoothed cube.

    Notes
    -----
    - Spatial smoothing uses ``lf3d_detect_spatial_sigma`` in *pixel* units after
      public aliases such as ``lf3d_detect_spatial_fwhm_arcsec`` have been
      normalized by :func:`_normalize_linefree_cfg_public_units`.
    - Spectral smoothing uses ``lf3d_detect_spectral_width`` in *channel* units
      after public aliases such as ``lf3d_detect_spectral_width_kms`` have been
      normalized. The default width of 3 is intentionally mild and detection-only.
    """
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(f"cube_data must be 3D (nchan,ny,nx), got shape={data.shape}")
    sigma_spatial = float(getattr(cfg, 'lf3d_detect_spatial_sigma', 0.0) or 0.0)
    if sigma_spatial <= 0.0:
        sigma_spatial = float(getattr(cfg, 'lf3d_seed_spatial_sigma', 0.0) or 0.0)
    out = data.astype(np.float32, copy=True)
    if sigma_spatial > 0.0:
        out = _gaussian_filter_spatial_nan_normalized(out, sigma_spatial)
    spectral_width = int(getattr(cfg, 'lf3d_detect_spectral_width', 0) or 0)
    if spectral_width > 1:
        out = _uniform_filter_spectral_nan_normalized(out, spectral_width)
    return out


def _spectral_dilate_no_wrap(mask: np.ndarray, pad: int) -> np.ndarray:
    arr = np.asarray(mask, dtype=bool)
    p = int(pad)
    if p <= 0:
        return arr
    if arr.ndim == 2:
        struct = np.ones((2 * p + 1, 1), dtype=bool)
    elif arr.ndim == 3:
        struct = np.ones((2 * p + 1, 1, 1), dtype=bool)
    else:
        raise ValueError(f"mask must be 2D or 3D, got shape={arr.shape}")
    return ndi.binary_dilation(arr, structure=struct)


def _spectral_hysteresis_no_wrap(strong: np.ndarray, weak: np.ndarray) -> np.ndarray:
    s = np.asarray(strong, dtype=bool)
    w = np.asarray(weak, dtype=bool)
    if s.shape != w.shape:
        raise ValueError(f"strong/weak shape mismatch: {s.shape} vs {w.shape}")
    if s.ndim == 2:
        struct = np.ones((3, 1), dtype=bool)
    elif s.ndim == 3:
        struct = np.ones((3, 1, 1), dtype=bool)
    else:
        raise ValueError(f"strong/weak mask must be 2D or 3D, got shape={s.shape}")
    return ndi.binary_propagation(s, structure=struct, mask=w)



def _prune_short_spectral_runs(mask: np.ndarray, *, velocity_axis_kms: Optional[np.ndarray], min_run_kms: Optional[float]) -> np.ndarray:
    """Remove per-LOS signal runs shorter than the requested contiguous velocity width.

    Parameters
    ----------
    mask : ndarray
        Boolean signal mask, shape ``(nchan, npix)`` or ``(nchan, ny, nx)``.
    velocity_axis_kms : ndarray or None
        Spectral axis in km/s. Required when ``min_run_kms`` is positive.
    min_run_kms : float or None
        Minimum contiguous velocity width to preserve. ``None`` or ``<= 0`` disables
        the cleanup. A run of ``N`` channels has width ``N * abs(dv)`` where
        ``dv`` is the median positive channel step of ``velocity_axis_kms``.
    """
    arr = np.asarray(mask, dtype=bool)
    thr = None if min_run_kms is None else float(min_run_kms)
    if thr is None or thr <= 0.0:
        return arr
    if velocity_axis_kms is None:
        raise ValueError("lf3d_min_run_kms requires a velocity axis in km/s for the current cube.")
    dv = _median_positive_step(np.asarray(velocity_axis_kms, dtype=float), name='velocity_axis_kms')
    if not np.isfinite(dv) or dv <= 0.0:
        raise ValueError("lf3d_min_run_kms requires a strictly positive spectral channel width in km/s.")
    if arr.ndim not in (2, 3):
        raise ValueError(f"mask must be 2D or 3D, got shape={arr.shape}")
    flat = arr.reshape(arr.shape[0], -1)
    if flat.size == 0:
        return arr
    structure = np.zeros((3, 3), dtype=bool)
    structure[:, 1] = True
    labels, nlab = ndi.label(flat, structure=structure)
    if nlab <= 0:
        return arr
    counts = np.bincount(labels.ravel(), minlength=nlab + 1)
    short = np.flatnonzero((np.arange(nlab + 1) > 0) & ((counts.astype(float) * float(dv)) < thr))
    if short.size == 0:
        return arr
    kill = np.zeros(nlab + 1, dtype=bool)
    kill[short] = True
    pruned = flat.copy()
    pruned[kill[labels]] = False
    return pruned.reshape(arr.shape)


def _masked_residual_location_sigma(
    masked_res: np.ndarray,
    *,
    positive_only: bool,
    sigma_mode: str = "negative",
) -> tuple[np.ndarray, np.ndarray]:
    """Per-column robust location and sigma from masked residuals.

    For positive-line emission, sigma_mode='negative' estimates sigma from the
    lower half only, which is less inflated by positive lines than symmetric MAD.
    """
    arr = np.asarray(masked_res, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        med = np.nanmedian(arr, axis=0)
    med = np.nan_to_num(med, nan=0.0).astype(np.float32, copy=False)

    mode = str(sigma_mode or "negative").strip().lower()
    use_negative = bool(positive_only) and mode in {"negative", "neg", "lower", "lower_half"}
    if use_negative:
        neg_dev = np.where(arr <= med[None, :], med[None, :] - arr, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mad = np.nanmedian(neg_dev, axis=0)
        neg_count = np.sum(np.isfinite(neg_dev), axis=0)
        need_fallback = (~np.isfinite(mad)) | (mad <= 0.0) | (neg_count < 4)
        if np.any(need_fallback):
            abs_dev = np.abs(arr - med[None, :])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                mad_sym = np.nanmedian(abs_dev, axis=0)
            mad = np.where(need_fallback, mad_sym, mad)
    else:
        abs_dev = np.abs(arr - med[None, :])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mad = np.nanmedian(abs_dev, axis=0)
    mad = np.nan_to_num(mad, nan=0.0).astype(np.float32, copy=False)
    sigma = (1.4826 * mad).astype(np.float32, copy=False)
    sigma = np.where(np.isfinite(sigma) & (sigma > 0.0), sigma, np.float32(1.0e-10)).astype(np.float32, copy=False)
    return med, sigma




def _diff_mad_sigma_per_spectrum(
    y: np.ndarray,
    valid: np.ndarray,
    *,
    dtype: np.dtype = np.dtype(np.float32),
) -> np.ndarray:
    """Estimate per-spectrum noise from channel differences.

    sigma ~= 1.4826 * MAD(diff) / sqrt(2)
    This is intentionally conservative against broad lines and slow baselines because
    differencing suppresses low-frequency structure before the robust scale estimate.
    """
    yd = np.asarray(y, dtype=dtype)
    vd = np.asarray(valid, dtype=bool)
    if yd.ndim != 2 or vd.shape != yd.shape:
        raise ValueError(f"y/valid shape mismatch in _diff_mad_sigma_per_spectrum: {yd.shape} vs {vd.shape}")
    if yd.shape[0] < 2:
        return np.full(yd.shape[1], dtype.type(1.0e-6), dtype=dtype)
    diff = yd[1:, :] - yd[:-1, :]
    diff_valid = vd[1:, :] & vd[:-1, :]
    diff_masked = np.where(diff_valid, diff, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        med = np.nanmedian(diff_masked, axis=0)
        mad = np.nanmedian(np.abs(diff_masked - med[None, :]), axis=0)
    sigma = (1.4826 * mad / np.sqrt(2.0)).astype(dtype, copy=False)
    # fallback for columns with too few valid diffs or zero MAD
    bad = (~np.isfinite(sigma)) | (sigma <= 0.0)
    if np.any(bad):
        masked = np.where(vd, yd, np.nan)
        _, sig_fb = _masked_residual_location_sigma(masked, positive_only=False, sigma_mode="symmetric")
        sigma[bad] = np.asarray(sig_fb, dtype=dtype)[bad]
    sigma = np.where(np.isfinite(sigma) & (sigma > 0.0), sigma, dtype.type(1.0e-6)).astype(dtype, copy=False)
    return sigma


def _lf3d_asymmetric_positive_weights(
    residual: np.ndarray,
    fit_mask: np.ndarray,
    *,
    positive_only: bool,
    sigma_mode: str,
    tau: Optional[float],
    floor: float,
) -> np.ndarray:
    """Asymmetric downweighting used only during 3D mask-estimation provisional fits.

    Positive residual excursions above ``tau * sigma`` are downweighted as
    ``tau*sigma / residual`` (Huber-like), while negative residuals keep weight 1.
    This biases provisional baselines toward the spectral floor and reduces the
    tendency to absorb drifting emission wings into the baseline model.
    """
    r = np.asarray(residual, dtype=np.float32)
    m = np.asarray(fit_mask, dtype=bool)
    if r.ndim != 2 or m.shape != r.shape:
        raise ValueError(f"residual/fit_mask shape mismatch: {r.shape} vs {m.shape}")
    tau_eff = None if tau is None else float(tau)
    if (not positive_only) or tau_eff is None or tau_eff <= 0.0:
        w = np.ones(r.shape, dtype=np.float32)
        w[~m] = 0.0
        return w
    masked_res = np.where(m, r, np.nan)
    med, sigma = _masked_residual_location_sigma(
        masked_res,
        positive_only=positive_only,
        sigma_mode=sigma_mode,
    )
    pos = r - med[None, :]
    thr = tau_eff * sigma[None, :]
    w = np.ones(r.shape, dtype=np.float32)
    valid = m & np.isfinite(pos) & np.isfinite(thr) & (thr > 0.0)
    hit = valid & (pos > thr)
    if np.any(hit):
        ratio = np.ones(r.shape, dtype=np.float32)
        np.divide(
            thr,
            np.maximum(pos, 1.0e-12),
            out=ratio,
            where=valid,
        )
        w[hit] = np.maximum(np.float32(max(0.0, float(floor))), ratio[hit])
    w[~m] = 0.0
    return w



def _make_lf3d_local_seed_3d(
    cube_data: np.ndarray,
    cfg: LineFreeConfig,
    *,
    positive_only: bool,
    pre_smoothed: bool = False,
) -> np.ndarray:
    """Build a 3D local hard seed for voxel-wise mask estimation.

    r11 seed variant:
    - use spatially-smoothed data only for Stage-A seed generation
    - estimate local sigma from channel differences (diff-MAD)
    - estimate the Stage-A seed baseline either with a stiff local linear fit
      (legacy r10 mode) or with a broad spectral median filter
    - expand the seed with seed-specific strong/weak hysteresis
    """
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(f"cube_data must be 3D (nchan,ny,nx), got shape={data.shape}")
    nchan, ny, nx = data.shape
    sigma_spatial = float(getattr(cfg, 'lf3d_seed_spatial_sigma', getattr(cfg, 'lf3d_spatial_sigma', 0.0)) or 0.0)
    if pre_smoothed:
        work = data.astype(np.float32, copy=False)
    else:
        work = _gaussian_filter_spatial_nan_normalized(data, sigma_spatial)
    flat = work.reshape(nchan, ny * nx)
    valid = np.isfinite(flat)
    Y = np.where(valid, flat, 0.0).astype(np.float32, copy=False)

    seed_order = int(getattr(cfg, 'lf3d_seed_baseline_order', 1))
    seed_sigma = float(getattr(cfg, 'lf3d_seed_sigma', 3.0))
    seed_hyst = getattr(cfg, 'lf3d_seed_hysteresis_sigma', 1.5)
    seed_hyst = None if seed_hyst is None else float(seed_hyst)
    seed_dilation = int(getattr(cfg, 'lf3d_seed_dilation', 0))
    seed_mode = str(getattr(cfg, 'lf3d_seed_baseline_mode', 'median') or 'median').strip().lower()
    seed_median_width = int(getattr(cfg, 'lf3d_seed_median_width', 71) or 71)
    if seed_median_width < 1:
        seed_median_width = 1
    if seed_median_width % 2 == 0:
        seed_median_width += 1
    seed_open_width = int(getattr(cfg, 'lf3d_seed_open_width', seed_median_width) or seed_median_width)
    if seed_open_width < 1:
        seed_open_width = 1
    if seed_open_width % 2 == 0:
        seed_open_width += 1
    chunk_pix = max(1, int(getattr(cfg, 'lf3d_chunk_pix', 4096)))

    x = np.linspace(-1.0, 1.0, nchan, dtype=np.float32)
    X_seed = np.vander(x, max(seed_order, 0) + 1, increasing=True).astype(np.float32, copy=False)

    seed = np.zeros((nchan, ny * nx), dtype=bool)
    for p0 in range(0, ny * nx, chunk_pix):
        p1 = min(ny * nx, p0 + chunk_pix)
        Yc = Y[:, p0:p1]
        valid_c = valid[:, p0:p1]

        sigma_true = _diff_mad_sigma_per_spectrum(Yc, valid_c, dtype=np.dtype(np.float32))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            med0 = np.nanmedian(np.where(valid_c, Yc, np.nan), axis=0)
        med0 = np.nan_to_num(med0, nan=0.0).astype(np.float32, copy=False)

        if seed_mode == 'median':
            Yseed = np.where(valid_c, Yc, med0[None, :]).astype(np.float32, copy=False)
            baseline = ndi.median_filter(Yseed, size=(seed_median_width, 1), mode='nearest').astype(np.float32, copy=False)
        elif seed_mode == 'opening':
            Yseed = np.where(valid_c, Yc, med0[None, :]).astype(np.float32, copy=False)
            baseline = ndi.grey_opening(Yseed, size=(seed_open_width, 1), mode='nearest').astype(np.float32, copy=False)
        elif seed_mode == 'linear':
            if positive_only:
                fit_mask = valid_c & (Yc <= (med0[None, :] + sigma_true[None, :]))
            else:
                fit_mask = valid_c & (np.abs(Yc - med0[None, :]) <= sigma_true[None, :])
            baseline = np.zeros_like(Yc, dtype=np.float32)
            min_pts = max(int(X_seed.shape[1]) + 2, 6)
            try:
                coef = _solve_weighted_least_squares_batch(
                    X_seed,
                    Yc,
                    fit_mask.astype(np.float32, copy=False),
                    dtype=np.dtype(np.float32),
                    method='qr',
                    rcond=None,
                    reg_eps=0.0,
                )
                baseline = np.einsum('vk,kn->vn', X_seed, coef, optimize=True).astype(np.float32, copy=False)
                counts = np.sum(fit_mask, axis=0)
                bad = counts < min_pts
                if np.any(bad):
                    baseline[:, bad] = med0[None, bad]
            except Exception:
                baseline[:] = med0[None, :]
        else:
            raise ValueError(f"Unknown lf3d_seed_baseline_mode={seed_mode!r}; expected 'median', 'opening' or 'linear'")

        resid = Yc - baseline
        if positive_only:
            strong = valid_c & (resid > (seed_sigma * sigma_true[None, :]))
            if seed_hyst is not None and seed_hyst > 0.0 and seed_hyst < seed_sigma:
                weak = valid_c & (resid > (seed_hyst * sigma_true[None, :]))
                line = _spectral_hysteresis_no_wrap(strong, weak)
            else:
                line = strong
        else:
            metric = np.abs(resid)
            strong = valid_c & (metric > (seed_sigma * sigma_true[None, :]))
            if seed_hyst is not None and seed_hyst > 0.0 and seed_hyst < seed_sigma:
                weak = valid_c & (metric > (seed_hyst * sigma_true[None, :]))
                line = _spectral_hysteresis_no_wrap(strong, weak)
            else:
                line = strong
        if seed_dilation > 0:
            line = _spectral_dilate_no_wrap(line, seed_dilation)
        seed[:, p0:p1] = line
    return seed.reshape(nchan, ny, nx)

def _solve_masked_least_squares_batch(
    A: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    dtype: np.dtype,
    reg_eps: float = 0.0,
) -> np.ndarray:
    Ad = np.asarray(A, dtype=dtype)
    Yd = np.asarray(Y, dtype=dtype)
    Wd = np.asarray(W, dtype=dtype)
    if Ad.ndim != 2 or Yd.ndim != 2 or Wd.shape != Yd.shape:
        raise ValueError("A/Y/W shape mismatch in _solve_masked_least_squares_batch")
    xtwx = np.einsum("vk,vn,vl->nkl", Ad, Wd, Ad, optimize=True)
    if reg_eps > 0.0:
        xtwx = xtwx + np.eye(Ad.shape[1], dtype=dtype)[None, :, :] * dtype.type(reg_eps)
    xtwy = np.einsum("vk,vn,vn->nk", Ad, Wd, Yd, optimize=True)
    coef = np.linalg.solve(xtwx, xtwy[:, :, None])[..., 0]
    return np.asarray(coef.T, dtype=dtype)


def _solve_weighted_least_squares_batch(
    A: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    dtype: np.dtype,
    method: str = "qr",
    rcond: Optional[float] = None,
    reg_eps: float = 0.0,
    max_batch_pix: Optional[int] = None,
) -> np.ndarray:
    """Solve batched weighted least-squares with selectable core solver."""
    method_eff = str(method or "qr").strip().lower()
    if method_eff not in {"qr", "normal"}:
        raise ValueError(f"Unknown weighted least-squares method={method!r}; expected 'qr' or 'normal'")
    if method_eff == "qr":
        return _solve_weighted_qr_batch(A, Y, W, rcond=rcond, dtype=dtype, max_batch_pix=max_batch_pix)
    return _solve_masked_least_squares_batch(A, Y, W, dtype=dtype, reg_eps=reg_eps)


def _iter_identical_mask_groups(mask: np.ndarray) -> Iterable[np.ndarray]:
    """Yield groups of spectrum indices that share the same boolean fit mask."""
    m = np.asarray(mask, dtype=bool)
    if m.ndim != 2:
        raise ValueError('mask must be 2D in _iter_identical_mask_groups')
    n_spec = int(m.shape[1])
    if n_spec == 0:
        return
    packed = np.packbits(np.ascontiguousarray(m.astype(np.uint8)), axis=0, bitorder='little').T
    _, inverse = np.unique(packed, axis=0, return_inverse=True)
    for gid in range(int(inverse.max()) + 1):
        idx = np.where(inverse == gid)[0]
        if idx.size > 0:
            yield idx.astype(np.int64, copy=False)


def _solve_unweighted_group_qr(
    A_lf: np.ndarray,
    Y_lf: np.ndarray,
    *,
    rcond: Optional[float],
    dtype: np.dtype,
) -> np.ndarray:
    """Solve multiple spectra sharing one binary fit mask with one QR factorization."""
    A64 = np.asarray(A_lf, dtype=np.float64)
    Y64 = np.asarray(Y_lf, dtype=np.float64)
    try:
        Q, R = np.linalg.qr(A64, mode='reduced')
        coef = np.linalg.solve(R, Q.T @ Y64)
    except Exception:
        coef, *_ = np.linalg.lstsq(A64, Y64, rcond=rcond)
    return np.asarray(coef, dtype=dtype)


def _solve_voxelmask_grouped_initial(
    A_full: np.ndarray,
    Y: np.ndarray,
    fit_mask: np.ndarray,
    *,
    dtype: np.dtype,
    solver: str,
    rcond: Optional[float],
    reg_eps: float,
    min_points: int,
    qr_batch_pix: Optional[int] = None,
    grouping: str = "auto",
    group_min_size: int = 12,
) -> np.ndarray:
    """Initial voxel-mask solve grouped by identical binary fit masks.

    Grouping is only beneficial when many spectra share exactly the same boolean
    mask. In the common case of mostly unique masks, the generic batched QR path
    is faster. Therefore the default is adaptive.
    """
    Af = np.asarray(A_full, dtype=dtype)
    Yd = np.asarray(Y, dtype=dtype)
    M = np.asarray(fit_mask, dtype=bool)
    if Af.ndim != 2 or Yd.ndim != 2 or M.shape != Yd.shape:
        raise ValueError('A/Y/fit_mask shape mismatch in _solve_voxelmask_grouped_initial')
    n_coef = int(Af.shape[1])
    n_spec = int(Yd.shape[1])
    solver_eff = str(solver or 'qr').strip().lower()
    grouping_eff = str(grouping or 'auto').strip().lower()
    if grouping_eff not in {'auto', 'always', 'off'}:
        raise ValueError(f'Unknown voxel_grouping={grouping!r}; expected auto/always/off')
    if grouping_eff == 'off':
        W = np.asarray(M, dtype=dtype)
        return _solve_weighted_least_squares_batch(
            Af, Yd, W, dtype=dtype, method=solver_eff, rcond=rcond, reg_eps=reg_eps, max_batch_pix=qr_batch_pix
        )
    packed = np.packbits(np.ascontiguousarray(M.astype(np.uint8)), axis=0, bitorder='little').T
    _, inverse = np.unique(packed, axis=0, return_inverse=True)
    counts = np.bincount(inverse)
    unique_count = int(counts.size)
    max_group = int(np.max(counts)) if counts.size > 0 else 0
    mean_group = float(n_spec) / float(max(1, unique_count))
    use_grouping = True
    if grouping_eff == 'auto':
        use_grouping = (
            (solver_eff == 'qr')
            and (n_spec >= 256)
            and (
                (max_group >= max(int(group_min_size), int(0.35 * n_spec)))
                or (mean_group >= float(max(8, group_min_size)))
                or (unique_count <= max(32, n_spec // 12))
            )
        )
    if not use_grouping:
        W = np.asarray(M, dtype=dtype)
        return _solve_weighted_least_squares_batch(
            Af, Yd, W, dtype=dtype, method=solver_eff, rcond=rcond, reg_eps=reg_eps, max_batch_pix=qr_batch_pix
        )
    coef = np.empty((n_coef, n_spec), dtype=dtype)
    for gid in range(unique_count):
        idx = np.where(inverse == gid)[0]
        if idx.size == 0:
            continue
        mj = M[:, int(idx[0])]
        if int(np.sum(mj)) < int(min_points):
            raise ValueError('insufficient points reached grouped initial solve unexpectedly')
        A_lf = np.asarray(Af[mj, :], dtype=dtype)
        Y_lf = np.asarray(Yd[mj, :][:, idx], dtype=dtype)
        if solver_eff == 'qr':
            coef[:, idx] = _solve_unweighted_group_qr(A_lf, Y_lf, rcond=rcond, dtype=dtype)
        else:
            W = np.ones_like(Y_lf, dtype=dtype)
            coef[:, idx] = _solve_weighted_least_squares_batch(
                A_lf,
                Y_lf,
                W,
                dtype=dtype,
                method=solver_eff,
                rcond=rcond,
                reg_eps=reg_eps,
                max_batch_pix=qr_batch_pix,
            )
    return coef

def _fit_single_spectrum_masked(
    y_full: np.ndarray,
    A_full: np.ndarray,
    fit_mask: np.ndarray,
    *,
    robust: bool,
    robust_iters: int,
    rcond: Optional[float],
    dtype: np.dtype,
    robust_early_stop: bool = False,
    robust_coef_rtol: float = 1.0e-4,
    robust_coef_atol: float = 1.0e-6,
) -> np.ndarray:
    yj = np.asarray(y_full, dtype=dtype)
    mj = np.asarray(fit_mask, dtype=bool).reshape(-1)
    if mj.shape[0] != yj.shape[0]:
        raise ValueError(f"fit_mask length mismatch: {mj.shape[0]} vs nchan={yj.shape[0]}")
    A_lf = np.asarray(A_full[mj, :], dtype=dtype)
    ylf = yj[mj]
    cj, *_ = np.linalg.lstsq(A_lf, ylf, rcond=rcond)
    cj = np.asarray(cj, dtype=dtype)
    if robust:
        early = bool(robust_early_stop)
        for _ in range(max(0, int(robust_iters))):
            rj = ylf - (A_lf @ cj)
            w = _robust_reweight(rj, dtype=dtype)
            sw = np.sqrt(np.clip(w, dtype.type(0.0), None)).astype(dtype, copy=False)
            Aw = A_lf * sw[:, None]
            bw = ylf * sw
            cj_new, *_ = np.linalg.lstsq(Aw, bw, rcond=rcond)
            cj_new = np.asarray(cj_new, dtype=dtype)
            if early and bool(_coef_converged_mask(cj, cj_new, rtol=robust_coef_rtol, atol=robust_coef_atol)[0]):
                cj = cj_new
                break
            cj = cj_new
    base = A_full @ cj
    return np.asarray(yj - base, dtype=np.float32)


def _fit_selected_robust_voxelmask(
    out_chunk: np.ndarray,
    Yc_fill: np.ndarray,
    fit_cols: np.ndarray,
    selected_rel: np.ndarray,
    *,
    fit_mask_chunk: np.ndarray,
    A_full: np.ndarray,
    coef_init: Optional[np.ndarray],
    rcond: Optional[float],
    dtype: np.dtype,
    robust_iters: int,
    robust_early_stop: bool,
    robust_coef_rtol: float,
    robust_coef_atol: float,
    batch_pixels: int,
    solver: str,
    reg_eps: float,
    qr_batch_pix: Optional[int],
) -> None:
    sel = np.asarray(selected_rel, dtype=np.int64)
    if sel.size == 0:
        return
    batch = max(1, int(batch_pixels))
    for p0 in range(0, sel.size, batch):
        p1 = min(sel.size, p0 + batch)
        sub_rel = sel[p0:p1]
        sub_cols = np.asarray(fit_cols[sub_rel], dtype=np.int64)
        Ysub = np.asarray(Yc_fill[:, sub_cols], dtype=dtype)
        W0 = np.asarray(fit_mask_chunk[:, sub_cols], dtype=dtype)
        try:
            if coef_init is not None and coef_init.ndim == 2 and coef_init.shape[1] >= int(np.max(sub_rel)) + 1:
                coef_sub = np.asarray(coef_init[:, sub_rel], dtype=dtype, copy=True)
            else:
                coef_sub = _solve_voxelmask_grouped_initial(
                    A_full,
                    Ysub,
                    W0 > 0,
                    dtype=dtype,
                    solver=solver,
                    rcond=rcond,
                    reg_eps=reg_eps,
                    min_points=max(int(A_full.shape[1]) + 2, 8),
                    qr_batch_pix=qr_batch_pix,
                    grouping="auto",
                    group_min_size=12,
                )
            active = np.arange(sub_cols.size, dtype=np.int64)
            early = bool(robust_early_stop)
            for _ in range(max(0, int(robust_iters))):
                if active.size == 0:
                    break
                coef_old = np.asarray(coef_sub[:, active], dtype=dtype, copy=True)
                Yact = np.asarray(Ysub[:, active], dtype=dtype)
                Wact0 = np.asarray(W0[:, active], dtype=dtype)
                resid = Yact - (A_full @ coef_old)
                resid_masked = np.where(Wact0 > 0, resid, np.nan)
                Wrob = _robust_reweight_matrix(resid_masked, dtype=dtype)
                Wact = np.asarray(Wact0 * Wrob, dtype=dtype)
                coef_new = _solve_weighted_least_squares_batch(
                    A_full,
                    Yact,
                    Wact,
                    dtype=dtype,
                    method=solver,
                    rcond=rcond,
                    reg_eps=reg_eps,
                    max_batch_pix=qr_batch_pix,
                )
                coef_sub[:, active] = coef_new
                if early:
                    converged = _coef_converged_mask(
                        coef_old,
                        coef_new,
                        rtol=robust_coef_rtol,
                        atol=robust_coef_atol,
                    )
                    if np.all(converged):
                        break
                    active = active[~converged]
            base = np.asarray(A_full @ coef_sub, dtype=dtype)
            out_chunk[:, sub_cols] = np.asarray(Ysub - base, dtype=np.float32)
        except Exception:
            for j in sub_cols:
                out_chunk[:, j] = _fit_single_spectrum_masked(
                    Yc_fill[:, j],
                    A_full,
                    fit_mask_chunk[:, j],
                    robust=True,
                    robust_iters=robust_iters,
                    rcond=rcond,
                    dtype=dtype,
                    robust_early_stop=robust_early_stop,
                    robust_coef_rtol=robust_coef_rtol,
                    robust_coef_atol=robust_coef_atol,
                )

def _combine_linefree_masks_conservative(masks: Sequence[np.ndarray], nchan: int) -> np.ndarray:
    """Combine line-free masks conservatively (intersection of line-free channels)."""
    combined = np.ones(int(nchan), dtype=bool)
    used = 0
    for m in masks:
        arr = np.asarray(m, dtype=bool)
        if arr.shape != (int(nchan),):
            raise ValueError(f"Line-free mask shape mismatch: {arr.shape} vs ({nchan},)")
        combined &= arr
        used += 1
    if used == 0:
        return np.ones(int(nchan), dtype=bool)
    return combined


def _parse_manual_windows(windows: Sequence[Union[str, Tuple[float, float]]]) -> List[Tuple[float, float]]:
    parsed: List[Tuple[float, float]] = []
    for w in windows:
        if isinstance(w, str):
            s = w.strip()
            if ":" in s:
                parts = s.split(":")
            elif "," in s:
                parts = s.split(",")
            else:
                parts = s.split()
            if len(parts) != 2:
                raise ValueError(f"Could not parse window {w!r}; expected 'lo:hi' or 'lo,hi'.")
            lo, hi = float(parts[0]), float(parts[1])
        else:
            lo, hi = float(w[0]), float(w[1])
        if lo > hi:
            lo, hi = hi, lo
        parsed.append((lo, hi))
    return parsed


def _manual_signal_mask_from_axis(
    axis: np.ndarray,
    windows: Sequence[Union[str, Tuple[float, float]]],
) -> np.ndarray:
    arr = np.asarray(axis, dtype=float)
    mask = np.zeros(arr.shape, dtype=bool)
    for lo, hi in _parse_manual_windows(windows):
        mask |= (arr >= float(lo)) & (arr <= float(hi))
    return mask


def _resolve_compute_dtype(bcfg: BaselineConfig) -> np.dtype:
    if bool(getattr(bcfg, "reproducible_mode", False)):
        return np.dtype(np.float64)
    dt = np.dtype(getattr(bcfg, "compute_dtype", "float32"))
    if dt.kind != "f" or dt.itemsize not in (4, 8):
        raise ValueError(f"compute_dtype must be float32 or float64, got {dt}")
    return dt


def _resolve_rcond(rcond: Optional[float], *, n_rows: int, n_cols: int, dtype: np.dtype, reproducible_mode: bool) -> Optional[float]:
    if rcond is not None:
        return float(rcond)
    if not reproducible_mode:
        return None
    eps = np.finfo(dtype).eps
    return float(eps * max(int(n_rows), int(n_cols)))


def _poly_coordinate(nchan: int, *, dtype: np.dtype, normalize: bool) -> np.ndarray:
    if not normalize:
        return np.arange(nchan, dtype=dtype)
    if nchan <= 1:
        return np.zeros(nchan, dtype=dtype)
    return np.linspace(-1.0, 1.0, nchan, dtype=dtype)


# -----------------------------------------------------------------------------
# 1) Line-free channel detection
# -----------------------------------------------------------------------------


def estimate_linefree_mask_1d(
    spectrum: np.ndarray,
    cfg: LineFreeConfig = LineFreeConfig(),
) -> np.ndarray:
    """
    Estimate line-free channels from a 1D spectrum by iterative robust sigma-clipping.

    Strategy
    --------
    - Estimate a slow baseline with a 1D median filter (smooth_width).
    - Compute residual = spectrum - baseline_est
    - Robust std from MAD, clip channels with |residual| > sigma * std
    - Iterate, recomputing std on remaining "candidate line-free" channels.
    - Dilate detected line channels by pad_chan.

    Returns
    -------
    linefree : np.ndarray (bool), shape (nchan,)
        True means "line-free and usable for baseline fitting".
    """
    y = np.asarray(spectrum, dtype=np.float64)
    n = y.size
    if n < 8:
        return np.isfinite(y)

    if bool(getattr(cfg, "interpolate_nans", True)):
        y_work = _interp_nans_1d(y, dtype=np.float64)
    else:
        y_work = np.nan_to_num(y, nan=0.0)

    base = _median_filter_1d(y_work, cfg.smooth_width)
    resid = y_work - base

    good = np.isfinite(y).copy()
    line = np.zeros(n, dtype=bool)

    for _ in range(int(cfg.iters)):
        cand = good & (~line)
        if cand.sum() < max(10, int(0.05 * n)):
            break
        s = _robust_std(resid[cand])
        if s <= 0:
            break
        new_line = np.abs(resid) > (float(cfg.sigma) * s)
        new_line &= good
        if np.array_equal(new_line, line):
            break
        line = new_line

    line = _dilate_1d(line, cfg.pad_chan)
    linefree = good & (~line)

    frac = float(linefree.sum()) / float(max(1, good.sum()))
    if frac < float(cfg.min_linefree_frac):
        raise ValueError(
            f"Too few line-free channels: frac={frac:.3f} < {cfg.min_linefree_frac:.3f}. "
            "Consider providing linefree_windows or loosening LineFreeConfig."
        )
    return linefree


def _estimate_linefree_mask_from_cube_agg1d(
    cube_data: np.ndarray,
    cfg: LineFreeConfig = LineFreeConfig(),
    *,
    agg: str = "median",
    max_pix: int = 200000,
    seed: Optional[int] = None,
    sample_mode: Optional[str] = None,
) -> np.ndarray:
    """Estimate a global 1D line-free mask from a 3D cube by spatial aggregation."""
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(f"cube_data must be 3D (nchan,ny,nx), got shape={data.shape}")

    nchan, ny, nx = data.shape
    npix = ny * nx
    flat = data.reshape(nchan, npix)

    mode = sample_mode if sample_mode is not None else getattr(cfg, "sampling_mode", "random")
    seed_eff = getattr(cfg, "sampling_seed", 0) if seed is None else int(seed)
    idx = _select_sample_indices(npix, int(max_pix), mode=mode, seed=seed_eff)
    flat_s = flat if idx is None else flat[:, idx]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if agg == "median":
            spec = np.nanmedian(flat_s, axis=1)
        elif agg == "mean":
            spec = np.nanmean(flat_s, axis=1)
        else:
            raise ValueError(f"Unknown agg={agg!r}")

    masks: List[np.ndarray] = [estimate_linefree_mask_1d(spec, cfg=cfg)]

    q = getattr(cfg, "conservative_quantile", None)
    if q is not None:
        q = float(q)
        if not (0.5 < q < 1.0):
            raise ValueError(f"conservative_quantile must satisfy 0.5 < q < 1.0, got {q}")
        quantiles: List[Tuple[str, float]] = [("upper", q)]
        if bool(getattr(cfg, "conservative_both_tails", True)):
            quantiles.append(("lower", 1.0 - q))
        for qname, qval in quantiles:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                spec_q = np.nanquantile(flat_s, qval, axis=1)
            try:
                masks.append(estimate_linefree_mask_1d(spec_q, cfg=cfg))
            except ValueError as exc:
                logging.warning(
                    "Skipping conservative %s-quantile line-candidate spectrum (q=%.3f): %s",
                    qname, qval, exc,
                )

    return _combine_linefree_masks_conservative(masks, nchan)


def estimate_linefree_mask_from_cube(
    cube_data: np.ndarray,
    cfg: LineFreeConfig = LineFreeConfig(),
    *,
    agg: str = "median",
    max_pix: int = 200000,
    seed: Optional[int] = None,
    sample_mode: Optional[str] = None,
    velocity_axis_kms: Optional[np.ndarray] = None,
    cell_arcsec: Optional[float] = None,
) -> np.ndarray:
    """Estimate line-free regions from a cube.

    Returns a 1D global mask when ``cfg.mask_kind='global_1d'`` / ``cfg.auto_method='agg_1d'``
    and a 3D voxel mask when ``cfg.mask_kind='voxel_3d'`` / ``cfg.auto_method='local_3d_poly'``.
    """
    cfg = _expand_linefree_profile(cfg)
    cfg = _normalize_linefree_cfg_public_units(
        cfg,
        velocity_axis_kms=velocity_axis_kms,
        cell_arcsec=cell_arcsec,
    )
    mask_kind = str(getattr(cfg, "mask_kind", "global_1d") or "global_1d").strip().lower()
    auto_method = str(getattr(cfg, "auto_method", "agg_1d") or "agg_1d").strip().lower()
    if auto_method == "local_3d_poly" or mask_kind == "voxel_3d":
        return estimate_linefree_mask_from_cube_3d(
            cube_data,
            cfg=cfg,
            velocity_axis_kms=velocity_axis_kms,
            cell_arcsec=cell_arcsec,
        )
    return _estimate_linefree_mask_from_cube_agg1d(
        cube_data,
        cfg=cfg,
        agg=agg,
        max_pix=max_pix,
        seed=seed,
        sample_mode=sample_mode,
    )



# -----------------------------------------------------------------------------
# 2) Ripple frequency estimation
# -----------------------------------------------------------------------------


def _normalize_ripple_apply_stage(stage: str) -> str:
    s = str(stage or "final").strip().lower()
    if s not in {"final", "mask_only", "both"}:
        raise ValueError(f"Unknown ripple_apply_stage={stage!r}; expected 'final', 'mask_only' or 'both'.")
    return s


def _positive_period_pair(val: Optional[Tuple[float, float]], *, name: str) -> Optional[Tuple[float, float]]:
    if val is None:
        return None
    a, b = float(val[0]), float(val[1])
    if (a <= 0.0) or (b <= 0.0) or (a > b):
        raise ValueError(f"Invalid {name}={val!r}; expected positive (lo<=hi).")
    return (a, b)


def _median_positive_step(axis: np.ndarray, *, name: str) -> float:
    arr = np.asarray(axis, dtype=float).reshape(-1)
    d = np.diff(arr)
    d = np.abs(d[np.isfinite(d)])
    d = d[d > 0.0]
    if d.size == 0:
        raise ValueError(f"Could not derive positive step from {name}.")
    return float(np.nanmedian(d))


_GAUSSIAN_FWHM_TO_SIGMA = 1.0 / np.sqrt(8.0 * np.log(2.0))


def _fwhm_to_sigma_pix(fwhm_pix: float, *, name: str) -> float:
    val = float(fwhm_pix)
    if not np.isfinite(val) or val < 0.0:
        raise ValueError(f"{name} must be finite and >= 0, got {fwhm_pix!r}")
    return float(val * _GAUSSIAN_FWHM_TO_SIGMA)


def _odd_width_from_chan(width_chan: float, *, name: str) -> int:
    val = float(width_chan)
    if not np.isfinite(val) or val <= 0.0:
        raise ValueError(f"{name} must be finite and > 0, got {width_chan!r}")
    w = max(1, int(round(val)))
    if w % 2 == 0:
        w += 1
    return int(w)


def _dilation_from_chan(width_chan: float, *, name: str) -> int:
    val = float(width_chan)
    if not np.isfinite(val) or val < 0.0:
        raise ValueError(f"{name} must be finite and >= 0, got {width_chan!r}")
    return max(0, int(round(val)))


def _header_cell_arcsec(header: fits.Header) -> Optional[float]:
    vals: List[float] = []
    for key in ('CDELT1', 'CDELT2'):
        try:
            v = float(header.get(key))
        except Exception:
            continue
        if np.isfinite(v) and v != 0.0:
            vals.append(abs(v) * 3600.0)
    if len(vals) >= 1:
        return float(np.nanmedian(np.asarray(vals, dtype=float)))
    try:
        cd11 = float(header.get('CD1_1'))
        cd21 = float(header.get('CD2_1'))
        cd12 = float(header.get('CD1_2'))
        cd22 = float(header.get('CD2_2'))
        s1 = np.hypot(cd11, cd21) * 3600.0
        s2 = np.hypot(cd12, cd22) * 3600.0
        vals2 = [v for v in (s1, s2) if np.isfinite(v) and v > 0.0]
        if vals2:
            return float(np.nanmedian(np.asarray(vals2, dtype=float)))
    except Exception:
        pass
    return None


def _linefree_cfg_needs_velocity_axis(cfg: Optional[LineFreeConfig]) -> bool:
    if cfg is None:
        return False
    return any(getattr(cfg, name, None) is not None for name in (
        'lf3d_seed_median_width_kms',
        'lf3d_seed_dilation_kms',
        'lf3d_detect_spectral_width_kms',
        'lf3d_min_run_kms',
    ))


def _linefree_cfg_needs_cell_arcsec(cfg: Optional[LineFreeConfig]) -> bool:
    if cfg is None:
        return False
    return any(getattr(cfg, name, None) is not None for name in (
        'lf3d_seed_spatial_fwhm_arcsec',
        'lf3d_detect_spatial_fwhm_arcsec',
        'lf3d_weak_spatial_fwhm_arcsec',
    ))


def _normalize_linefree_cfg_public_units(
    cfg: LineFreeConfig,
    *,
    velocity_axis_kms: Optional[np.ndarray] = None,
    cell_arcsec: Optional[float] = None,
) -> LineFreeConfig:
    """Resolve public-facing unit aliases to the legacy internal config fields.

    Public / documented names use:
    - arcsec for spatial smoothing FWHM
    - km/s for spectral widths
    - *_threshold_sigma for noise-based thresholds

    Internal implementation remains on pixel / channel units and legacy field names.
    """
    updates: dict[str, Any] = {}

    # Threshold aliases (dimensionless, no axis context required).
    alias_map = {
        'lf3d_seed_threshold_sigma': 'lf3d_seed_sigma',
        'lf3d_seed_hysteresis_threshold_sigma': 'lf3d_seed_hysteresis_sigma',
        'lf3d_detect_threshold_sigma': 'lf3d_threshold',
        'lf3d_hysteresis_threshold_sigma': 'lf3d_hysteresis_sigma',
    }
    for public_name, internal_name in alias_map.items():
        val = getattr(cfg, public_name, None)
        if val is not None:
            updates[internal_name] = float(val)

    # Spatial widths. Public form is FWHM; internal form is gaussian sigma in pixels.
    val = getattr(cfg, 'lf3d_seed_spatial_fwhm_pix', None)
    if val is not None:
        updates['lf3d_seed_spatial_sigma'] = _fwhm_to_sigma_pix(val, name='lf3d_seed_spatial_fwhm_pix')
    val = getattr(cfg, 'lf3d_seed_spatial_fwhm_arcsec', None)
    if val is not None:
        if cell_arcsec is None:
            raise ValueError('lf3d_seed_spatial_fwhm_arcsec requires a spatial cell size (arcsec/pix).')
        updates['lf3d_seed_spatial_sigma'] = _fwhm_to_sigma_pix(float(val) / float(cell_arcsec), name='lf3d_seed_spatial_fwhm_arcsec')

    val = getattr(cfg, 'lf3d_detect_spatial_fwhm_pix', None)
    if val is not None:
        updates['lf3d_detect_spatial_sigma'] = _fwhm_to_sigma_pix(val, name='lf3d_detect_spatial_fwhm_pix')
    val = getattr(cfg, 'lf3d_detect_spatial_fwhm_arcsec', None)
    if val is not None:
        if cell_arcsec is None:
            raise ValueError('lf3d_detect_spatial_fwhm_arcsec requires a spatial cell size (arcsec/pix).')
        updates['lf3d_detect_spatial_sigma'] = _fwhm_to_sigma_pix(float(val) / float(cell_arcsec), name='lf3d_detect_spatial_fwhm_arcsec')

    val = getattr(cfg, 'lf3d_weak_spatial_fwhm_pix', None)
    if val is not None:
        updates['lf3d_weak_spatial_sigma'] = _fwhm_to_sigma_pix(val, name='lf3d_weak_spatial_fwhm_pix')
    val = getattr(cfg, 'lf3d_weak_spatial_fwhm_arcsec', None)
    if val is not None:
        if cell_arcsec is None:
            raise ValueError('lf3d_weak_spatial_fwhm_arcsec requires a spatial cell size (arcsec/pix).')
        updates['lf3d_weak_spatial_sigma'] = _fwhm_to_sigma_pix(float(val) / float(cell_arcsec), name='lf3d_weak_spatial_fwhm_arcsec')

    # Spectral widths. Public form is km/s; internal form remains channel counts.
    dv = None
    if velocity_axis_kms is not None:
        dv = _median_positive_step(np.asarray(velocity_axis_kms, dtype=float), name='velocity_axis_kms')

    val = getattr(cfg, 'lf3d_seed_median_width_chan', None)
    if val is not None:
        updates['lf3d_seed_median_width'] = _odd_width_from_chan(val, name='lf3d_seed_median_width_chan')
    val = getattr(cfg, 'lf3d_seed_median_width_kms', None)
    if val is not None:
        if dv is None:
            raise ValueError('lf3d_seed_median_width_kms requires a velocity axis in km/s for the current cube.')
        updates['lf3d_seed_median_width'] = _odd_width_from_chan(float(val) / float(dv), name='lf3d_seed_median_width_kms')

    val = getattr(cfg, 'lf3d_seed_dilation_chan', None)
    if val is not None:
        updates['lf3d_seed_dilation'] = _dilation_from_chan(val, name='lf3d_seed_dilation_chan')
    val = getattr(cfg, 'lf3d_seed_dilation_kms', None)
    if val is not None:
        if dv is None:
            raise ValueError('lf3d_seed_dilation_kms requires a velocity axis in km/s for the current cube.')
        updates['lf3d_seed_dilation'] = _dilation_from_chan(float(val) / float(dv), name='lf3d_seed_dilation_kms')

    val = getattr(cfg, 'lf3d_detect_spectral_width_chan', None)
    if val is not None:
        updates['lf3d_detect_spectral_width'] = _odd_width_from_chan(val, name='lf3d_detect_spectral_width_chan')
    val = getattr(cfg, 'lf3d_detect_spectral_width_kms', None)
    if val is not None:
        if dv is None:
            raise ValueError('lf3d_detect_spectral_width_kms requires a velocity axis in km/s for the current cube.')
        updates['lf3d_detect_spectral_width'] = _odd_width_from_chan(float(val) / float(dv), name='lf3d_detect_spectral_width_kms')

    return replace(cfg, **updates) if updates else cfg


def _resolve_ripple_period_range_chan(
    rcfg: RippleConfig,
    *,
    velocity_axis_kms: Optional[np.ndarray] = None,
    spectral_axis_hz: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    pr_kms = _positive_period_pair(getattr(rcfg, "period_range_kms", None), name="period_range_kms")
    if pr_kms is not None:
        if velocity_axis_kms is None:
            raise ValueError("period_range_kms requires a velocity axis in km/s for the current cube.")
        dv = _median_positive_step(np.asarray(velocity_axis_kms, dtype=float), name="velocity_axis_kms")
        return (pr_kms[0] / dv, pr_kms[1] / dv)
    pr_hz = _positive_period_pair(getattr(rcfg, "period_range_hz", None), name="period_range_hz")
    if pr_hz is not None:
        if spectral_axis_hz is None:
            raise ValueError("period_range_hz requires a spectral axis in Hz for the current cube.")
        dnu = _median_positive_step(np.asarray(spectral_axis_hz, dtype=float), name="spectral_axis_hz")
        return (pr_hz[0] / dnu, pr_hz[1] / dnu)
    per = _positive_period_pair(getattr(rcfg, "period_range_chan", None), name="period_range_chan")
    if per is None:
        raise ValueError("RippleConfig must define one of period_range_chan/period_range_hz/period_range_kms.")
    return per


def _spectral_axis_hz_from_linear_header(header: fits.Header, nchan: int) -> np.ndarray:
    ctype3 = str(header.get("CTYPE3", "") or "").upper()
    if "FREQ" not in ctype3:
        raise ValueError("Header spectral axis is not frequency-like; cannot build Hz axis from linear header.")
    crpix = float(header.get("CRPIX3", 1.0))
    crval = float(header.get("CRVAL3", 0.0))
    cdelt = float(header.get("CDELT3", 1.0))
    cunit = str(header.get("CUNIT3", "Hz") or "Hz")
    unit = u.Unit(cunit)
    idx = np.arange(int(nchan), dtype=float) + 1.0
    vals = crval + (idx - crpix) * cdelt
    return np.asarray((vals * unit).to_value(u.Hz), dtype=float)


def _spectral_axis_hz_from_fits_context(
    input_fits: str,
    header: fits.Header,
    *,
    cube_ext: Optional[Union[int, str]],
    hdu_index: int,
    nchan: int,
) -> np.ndarray:
    last_err: Optional[Exception] = None
    if SpectralCube is not None:
        try:
            sc = SpectralCube.read(input_fits, hdu=(cube_ext if cube_ext is not None else hdu_index))
            axis_hz = np.asarray(sc.spectral_axis.to_value(u.Hz), dtype=float)
            if axis_hz.shape == (nchan,):
                return axis_hz
        except Exception as exc:
            last_err = exc
    try:
        axis_hz = _spectral_axis_hz_from_linear_header(header, nchan)
        if axis_hz.shape == (nchan,):
            return axis_hz
    except Exception as exc:
        last_err = exc
    if last_err is not None:
        raise ValueError("Could not build spectral axis in Hz for ripple period_range_hz.") from last_err
    raise ValueError("Could not build spectral axis in Hz for ripple period_range_hz.")


def _fit_detection_cube_with_safe_windows(
    cube_data: np.ndarray,
    *,
    safe_mask_1d: np.ndarray,
    baseline_cfg: BaselineConfig,
    ripple_freqs: Sequence[float] | None,
) -> np.ndarray:
    data = np.asarray(cube_data, dtype=np.float32)
    safe = np.asarray(safe_mask_1d, dtype=bool).reshape(-1)
    if data.ndim != 3:
        raise ValueError("cube_data must be 3D in _fit_detection_cube_with_safe_windows")
    if safe.shape != (data.shape[0],):
        raise ValueError(f"safe_mask_1d shape mismatch: {safe.shape} vs ({data.shape[0]},)")
    if not np.any(safe):
        return data.copy()
    det_cfg = BaselineConfig(
        poly_order=int(baseline_cfg.poly_order),
        ripple=bool(baseline_cfg.ripple),
        robust=False,
        robust_mode=str(getattr(baseline_cfg, "robust_mode", "selective")),
        robust_iters=0,
        robust_early_stop=False,
        robust_coef_rtol=float(getattr(baseline_cfg, "robust_coef_rtol", 1.0e-4)),
        robust_coef_atol=float(getattr(baseline_cfg, "robust_coef_atol", 1.0e-6)),
        robust_selective_sigma=float(getattr(baseline_cfg, "robust_selective_sigma", 4.5)),
        robust_selective_frac=float(getattr(baseline_cfg, "robust_selective_frac", 0.10)),
        robust_selective_max_pixels=int(getattr(baseline_cfg, "robust_selective_max_pixels", 2048)),
        robust_batch_pixels=int(getattr(baseline_cfg, "robust_batch_pixels", 512)),
        rcond=getattr(baseline_cfg, "rcond", None),
        chunk_pix=int(getattr(baseline_cfg, "chunk_pix", 4096)),
        reproducible_mode=bool(getattr(baseline_cfg, "reproducible_mode", False)),
        compute_dtype=str(getattr(baseline_cfg, "compute_dtype", "float32")),
        normalize_x=bool(getattr(baseline_cfg, "normalize_x", False)),
        strict_failures=False,
        fallback_to_pixelwise=True,
        voxel_solver=str(getattr(baseline_cfg, "voxel_solver", "qr")),
        voxel_qr_batch_pix=int(getattr(baseline_cfg, "voxel_qr_batch_pix", 2048)),
        voxel_grouping=str(getattr(baseline_cfg, "voxel_grouping", "auto")),
        voxel_group_min_size=int(getattr(baseline_cfg, "voxel_group_min_size", 12)),
    )
    det_cube, _, _ = subtract_baseline_cube(
        data,
        linefree_mask=np.asarray(safe, dtype=bool),
        bcfg=det_cfg,
        ripple_freqs=list(ripple_freqs or []),
        return_qc=False,
    )
    return np.asarray(det_cube, dtype=np.float32)


def estimate_ripple_frequencies_fft(
    spectrum: np.ndarray,
    linefree_mask: np.ndarray,
    rcfg: RippleConfig = RippleConfig(),
    *,
    poly_order_pre: int = 1,
    velocity_axis_kms: Optional[np.ndarray] = None,
    spectral_axis_hz: Optional[np.ndarray] = None,
) -> List[float]:
    """
    Estimate dominant ripple frequencies (cycles/channel) using FFT.

    Steps
    -----
    - Fit and subtract a low-order polynomial on line-free channels (poly_order_pre)
    - Fill or zero-out non line-free channels
    - Apply window (hann) to reduce leakage
    - FFT power spectrum
    - Pick top-N peaks within period_range_chan, enforcing min_separation

    Returns
    -------
    freqs : list of float
        Frequencies in cycles/channel (0..0.5). Length <= rcfg.nfreq.
    """
    y = np.asarray(spectrum, dtype=np.float64)
    m = np.asarray(linefree_mask, dtype=bool)
    n = y.size
    if n < 16:
        return []
    if m.shape != (n,):
        raise ValueError(f"linefree_mask shape mismatch: {m.shape} vs ({n},)")
    if m.sum() < max(8, poly_order_pre + 2):
        return []

    # Use a normalized coordinate for the polynomial pre-fit to improve conditioning,
    # but keep ripple frequencies in cycles/channel.
    x_poly = _poly_coordinate(n, dtype=np.float64, normalize=True)
    if poly_order_pre >= 0 and m.sum() >= (poly_order_pre + 2):
        V = np.vstack([x_poly[m] ** k for k in range(poly_order_pre + 1)]).T
        coeff, *_ = np.linalg.lstsq(V, y[m], rcond=float(np.finfo(np.float64).eps * max(V.shape)))
        base = np.zeros_like(y)
        for k in range(poly_order_pre + 1):
            base += float(coeff[k]) * (x_poly ** k)
        r = y - base
    else:
        r = y.copy()

    if bool(getattr(rcfg, "interpolate_masked", False)):
        r2 = np.full_like(r, np.nan)
        r2[m] = r[m]
        r2 = _interp_nans_1d(r2, dtype=np.float64)
    else:
        r2 = np.zeros_like(r)
        r2[m] = r[m]

    if rcfg.window.lower() == "hann":
        w = np.hanning(n).astype(np.float64)
        r2 = r2 * w
    elif rcfg.window.lower() != "none":
        raise ValueError(f"Unknown ripple window={rcfg.window!r}")

    spec = np.fft.rfft(r2)
    pwr = (spec.real ** 2 + spec.imag ** 2).astype(np.float64, copy=False)
    freqs = np.fft.rfftfreq(n, d=1.0)

    per_lo, per_hi = _resolve_ripple_period_range_chan(
        rcfg,
        velocity_axis_kms=velocity_axis_kms,
        spectral_axis_hz=spectral_axis_hz,
    )
    fmin = 1.0 / float(per_hi)
    fmax = 1.0 / float(per_lo)
    sel = (freqs >= fmin) & (freqs <= fmax)
    if sel.size > 0:
        sel[0] = False
    if sel.sum() < 1:
        return []

    cand_idx = np.where(sel)[0]
    if bool(getattr(rcfg, "stable_sort", True)):
        order = np.argsort(-pwr[cand_idx], kind="stable")
    else:
        order = np.argsort(pwr[cand_idx])[::-1]
    cand_idx = cand_idx[order]

    chosen: List[float] = []
    for i in cand_idx:
        f = float(freqs[i])
        if not np.isfinite(f) or f <= 0:
            continue
        if any(abs(f - fc) < float(rcfg.min_separation) for fc in chosen):
            continue
        chosen.append(f)
        if len(chosen) >= int(rcfg.nfreq):
            break

    return chosen


# -----------------------------------------------------------------------------
# 3) Baseline subtraction core
# -----------------------------------------------------------------------------


def _design_matrix(
    nchan: int,
    *,
    poly_order: int,
    ripple_freqs: Optional[Sequence[float]] = None,
    dtype: np.dtype = np.float32,
    normalize_poly_x: bool = False,
) -> np.ndarray:
    """
    Build design matrix A(x) for model:
      baseline(x) = sum_{k=0..poly_order} c_k x^k  +  sum_j (a_j sin(2π f_j x) + b_j cos(2π f_j x))

    Polynomial columns may use a normalized coordinate in [-1,1] to improve conditioning.
    Ripple terms always use the native channel index so that frequencies remain cycles/channel.

    Returns A with shape (nchan, ncoef).
    """
    x_poly = _poly_coordinate(int(nchan), dtype=dtype, normalize=bool(normalize_poly_x))
    x_phase = np.arange(int(nchan), dtype=dtype)
    cols: List[np.ndarray] = []
    for k in range(int(poly_order) + 1):
        cols.append(x_poly ** k)

    if ripple_freqs:
        for f in ripple_freqs:
            w = dtype.type(2.0 * np.pi * float(f))
            cols.append(np.sin(w * x_phase, dtype=dtype))
            cols.append(np.cos(w * x_phase, dtype=dtype))

    return np.vstack(cols).T.astype(dtype, copy=False)


def _robust_reweight(
    resid: np.ndarray,
    *,
    c: float = 1.345,
    dtype: np.dtype = np.float32,
) -> np.ndarray:
    """Huber weights for residuals (vector)."""
    r = np.asarray(resid, dtype=dtype)
    s = _robust_std(r)
    if s <= 0:
        return np.ones_like(r, dtype=dtype)
    t = np.abs(r) / (dtype.type(c * s))
    w = np.ones_like(r, dtype=dtype)
    m = t > 1
    w[m] = 1.0 / t[m]
    return w


def _robust_reweight_matrix(
    resid: np.ndarray,
    *,
    c: float = 1.345,
    dtype: np.dtype = np.float32,
) -> np.ndarray:
    """Huber weights for residual matrix with residuals in columns.

    Parameters
    ----------
    resid : ndarray, shape (n_samples, n_spec)
        Residuals for multiple spectra.
    """
    r = np.asarray(resid, dtype=np.float64)
    if r.ndim != 2:
        raise ValueError('resid must be 2D in _robust_reweight_matrix')
    med = np.nanmedian(r, axis=0, keepdims=True)
    mad = np.nanmedian(np.abs(r - med), axis=0)
    scale = 1.4826 * mad
    bad = (~np.isfinite(scale)) | (scale <= 0)
    if np.any(bad):
        scale2 = np.nanstd(r, axis=0)
        scale = np.where(bad, scale2, scale)
    bad = (~np.isfinite(scale)) | (scale <= 0)
    scale = np.where(bad, 1.0, scale)
    t = np.abs(r) / (float(c) * scale[None, :])
    w = np.ones_like(r, dtype=np.float64)
    m = t > 1.0
    w[m] = 1.0 / t[m]
    if np.any(bad):
        w[:, bad] = 1.0
    return np.asarray(w, dtype=dtype)


def _coef_converged_mask(
    coef_old: np.ndarray,
    coef_new: np.ndarray,
    *,
    rtol: float,
    atol: float,
) -> np.ndarray:
    """Return per-spectrum convergence mask based on coefficient change."""
    a = np.asarray(coef_old, dtype=np.float64)
    b = np.asarray(coef_new, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError('coef_old / coef_new shape mismatch in _coef_converged_mask')
    if a.ndim == 1:
        a = a[:, None]
        b = b[:, None]
    delta = np.max(np.abs(b - a), axis=0)
    scale = float(atol) + float(rtol) * np.maximum(np.max(np.abs(a), axis=0), np.max(np.abs(b), axis=0))
    return np.asarray(delta <= scale, dtype=bool)


def _solve_weighted_qr_batch(
    A: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    rcond: Optional[float],
    dtype: np.dtype,
    max_batch_pix: Optional[int] = None,
) -> np.ndarray:
    """Solve weighted least-squares for a batch of spectra using QR."""
    A64 = np.asarray(A, dtype=np.float64)
    Y64 = np.asarray(Y, dtype=np.float64)
    W64 = np.asarray(W, dtype=np.float64)
    if A64.ndim != 2 or Y64.ndim != 2 or W64.shape != Y64.shape:
        raise ValueError('A/Y/W shape mismatch in _solve_weighted_qr_batch')
    n_spec = int(Y64.shape[1])
    n_coef = int(A64.shape[1])
    if n_spec == 0:
        return np.empty((n_coef, 0), dtype=dtype)
    batch = n_spec if (max_batch_pix is None or int(max_batch_pix) <= 0) else max(1, int(max_batch_pix))
    coef = np.empty((n_coef, n_spec), dtype=np.float64)
    for p0 in range(0, n_spec, batch):
        p1 = min(n_spec, p0 + batch)
        sqrtW = np.sqrt(np.clip(W64[:, p0:p1], 0.0, None))
        Aw = sqrtW.T[:, :, None] * A64[None, :, :]
        bw = sqrtW.T[:, :, None] * Y64[:, p0:p1].T[:, :, None]
        try:
            Q, R = np.linalg.qr(Aw, mode='reduced')
            Qtb = np.matmul(np.swapaxes(Q, -1, -2), bw)
            coef[:, p0:p1] = np.linalg.solve(R, Qtb)[..., 0].T
        except Exception:
            for ib in range(p1 - p0):
                Ai = Aw[ib]
                bi = bw[ib, :, 0]
                try:
                    Qi, Ri = np.linalg.qr(Ai, mode='reduced')
                    coef[:, p0 + ib] = np.linalg.solve(Ri, Qi.T @ bi)
                except Exception:
                    coef[:, p0 + ib], *_ = np.linalg.lstsq(Ai, bi, rcond=rcond)
    return np.asarray(coef, dtype=dtype)

def _fit_selected_robust_batched(
    out_chunk: np.ndarray,
    Yc_fill: np.ndarray,
    fit_cols: np.ndarray,
    selected_rel: np.ndarray,
    *,
    A_full: np.ndarray,
    A_lf: np.ndarray,
    linefree_mask: np.ndarray,
    coef_init: Optional[np.ndarray],
    rcond: Optional[float],
    dtype: np.dtype,
    robust_iters: int,
    batch_pixels: int,
    robust_early_stop: bool,
    robust_coef_rtol: float,
    robust_coef_atol: float,
) -> None:
    """Robust re-fit selected spectra using batched IRLS on weighted least-squares via QR."""
    sel = np.asarray(selected_rel, dtype=np.int64)
    if sel.size == 0:
        return
    Ylf_all = np.asarray(Yc_fill[linefree_mask, :], dtype=dtype)
    batch = max(1, int(batch_pixels))
    early = bool(robust_early_stop)
    for p0 in range(0, sel.size, batch):
        p1 = min(sel.size, p0 + batch)
        sub_rel = sel[p0:p1]
        sub_cols = np.asarray(fit_cols[sub_rel], dtype=np.int64)
        Ysub = np.asarray(Ylf_all[:, sub_cols], dtype=dtype)
        if coef_init is not None and coef_init.shape[1] >= np.max(sub_rel) + 1:
            coef_sub = np.asarray(coef_init[:, sub_rel], dtype=dtype, copy=True)
        else:
            coef_sub, *_ = np.linalg.lstsq(A_lf, Ysub, rcond=rcond)
            coef_sub = np.asarray(coef_sub, dtype=dtype)
        active = np.arange(sub_rel.size, dtype=np.int64)
        for _ in range(max(0, int(robust_iters))):
            if active.size == 0:
                break
            coef_old = np.asarray(coef_sub[:, active], dtype=dtype, copy=True)
            Yact = np.asarray(Ysub[:, active], dtype=dtype)
            resid = Yact - (A_lf @ coef_old)
            W = _robust_reweight_matrix(resid, dtype=dtype)
            coef_new = _solve_weighted_qr_batch(A_lf, Yact, W, rcond=rcond, dtype=dtype)
            coef_sub[:, active] = coef_new
            if early:
                converged = _coef_converged_mask(
                    coef_old,
                    coef_new,
                    rtol=robust_coef_rtol,
                    atol=robust_coef_atol,
                )
                if np.all(converged):
                    break
                active = active[~converged]
        base = np.asarray(A_full @ coef_sub, dtype=dtype)
        out_chunk[:, sub_cols] = np.asarray(Yc_fill[:, sub_cols] - base, dtype=np.float32)


def _robust_loc_scale(x: np.ndarray) -> tuple[float, float]:
    arr = np.asarray(x, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0, 0.0
    med = float(np.nanmedian(arr))
    mad = float(np.nanmedian(np.abs(arr - med)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanstd(arr))
    if not np.isfinite(scale) or scale <= 0:
        scale = 0.0
    return med, scale


def _selective_robust_candidates(
    resid_lf: np.ndarray,
    finite_lf: np.ndarray,
    *,
    sigma: float,
    max_frac: float,
    max_pixels: int,
) -> np.ndarray:
    """Pick spectra that merit a robust re-fit from an initial ordinary fit.

    Parameters
    ----------
    resid_lf : array, shape (n_lf, n_spec)
        Residuals on line-free channels from the initial ordinary fit.
    finite_lf : bool array, shape (n_lf, n_spec)
        Finite mask on line-free channels.
    sigma : float
        Threshold in robust-sigma units.
    max_frac : float
        Maximum fraction of spectra to refit.
    max_pixels : int
        Hard cap on the number of spectra to refit.
    """
    r = np.asarray(resid_lf, dtype=np.float64)
    f = np.asarray(finite_lf, dtype=bool)
    if r.ndim != 2 or f.shape != r.shape:
        raise ValueError('resid_lf / finite_lf shape mismatch in selective robust picker')
    if r.shape[1] == 0:
        return np.empty(0, dtype=np.int64)

    abs_r = np.abs(r)
    abs_r[~f] = np.nan
    cnt = np.sum(np.isfinite(abs_r), axis=0)
    good = cnt > 0
    if not np.any(good):
        return np.empty(0, dtype=np.int64)

    sq = np.where(np.isfinite(abs_r), abs_r * abs_r, 0.0)
    rms = np.full(abs_r.shape[1], np.nan, dtype=np.float64)
    rms[good] = np.sqrt(np.sum(sq[:, good], axis=0) / cnt[good])

    peak_work = np.where(np.isfinite(abs_r), abs_r, -np.inf)
    peak = np.max(peak_work, axis=0)
    peak[~good] = np.nan

    med_rms, scl_rms = _robust_loc_scale(rms[good])
    med_peak, scl_peak = _robust_loc_scale(peak[good])
    cand = np.zeros(abs_r.shape[1], dtype=bool)

    score = np.zeros(abs_r.shape[1], dtype=np.float64)
    if scl_rms > 0:
        z_rms = (rms - med_rms) / scl_rms
        cand |= np.isfinite(z_rms) & (z_rms > float(sigma))
        score = np.maximum(score, np.where(np.isfinite(z_rms), z_rms, 0.0))
    if scl_peak > 0:
        z_peak = (peak - med_peak) / scl_peak
        cand |= np.isfinite(z_peak) & (z_peak > float(sigma))
        score = np.maximum(score, np.where(np.isfinite(z_peak), z_peak, 0.0))

    idx = np.where(cand)[0]
    if idx.size == 0:
        return idx.astype(np.int64, copy=False)

    limit_frac = int(np.ceil(float(max_frac) * int(np.sum(good)))) if float(max_frac) > 0 else idx.size
    limit = idx.size
    if limit_frac > 0:
        limit = min(limit, limit_frac)
    if int(max_pixels) > 0:
        limit = min(limit, int(max_pixels))
    if limit <= 0:
        return np.empty(0, dtype=np.int64)
    if idx.size > limit:
        order = np.argsort(score[idx], kind='mergesort')[::-1]
        idx = idx[order[:limit]]
    return np.asarray(np.sort(idx), dtype=np.int64)


def _fit_selected_robust(
    out_chunk: np.ndarray,
    Yc_fill: np.ndarray,
    fit_cols: np.ndarray,
    selected_rel: np.ndarray,
    *,
    A_full: np.ndarray,
    A_lf: np.ndarray,
    linefree_mask: np.ndarray,
    coef_init: Optional[np.ndarray] = None,
    rcond: Optional[float],
    dtype: np.dtype,
    robust_iters: int,
    batch_pixels: int,
    robust_early_stop: bool,
    robust_coef_rtol: float,
    robust_coef_atol: float,
) -> None:
    _fit_selected_robust_batched(
        out_chunk,
        Yc_fill,
        fit_cols,
        selected_rel,
        A_full=A_full,
        A_lf=A_lf,
        linefree_mask=linefree_mask,
        coef_init=coef_init,
        rcond=rcond,
        dtype=dtype,
        robust_iters=robust_iters,
        batch_pixels=batch_pixels,
        robust_early_stop=robust_early_stop,
        robust_coef_rtol=robust_coef_rtol,
        robust_coef_atol=robust_coef_atol,
    )


def _fit_single_spectrum(
    y_full: np.ndarray,
    A_full: np.ndarray,
    A_lf: np.ndarray,
    linefree_mask: np.ndarray,
    *,
    robust: bool,
    robust_iters: int,
    rcond: Optional[float],
    dtype: np.dtype,
    robust_early_stop: bool = False,
    robust_coef_rtol: float = 1.0e-4,
    robust_coef_atol: float = 1.0e-6,
) -> np.ndarray:
    yj = np.asarray(y_full, dtype=dtype)
    ylf = yj[linefree_mask]
    cj, *_ = np.linalg.lstsq(A_lf, ylf, rcond=rcond)
    cj = np.asarray(cj, dtype=dtype)
    if robust:
        early = bool(robust_early_stop)
        for _ in range(max(0, int(robust_iters))):
            rj = ylf - (A_lf @ cj)
            w = _robust_reweight(rj, dtype=dtype)
            sw = np.sqrt(np.clip(w, dtype.type(0.0), None)).astype(dtype, copy=False)
            Aw = A_lf * sw[:, None]
            bw = ylf * sw
            cj_new, *_ = np.linalg.lstsq(Aw, bw, rcond=rcond)
            cj_new = np.asarray(cj_new, dtype=dtype)
            if early and bool(_coef_converged_mask(cj, cj_new, rtol=robust_coef_rtol, atol=robust_coef_atol)[0]):
                cj = cj_new
                break
            cj = cj_new
    base = A_full @ cj
    return np.asarray(yj - base, dtype=np.float32)


def _subtract_baseline_cube_voxelmask(
    cube_data: np.ndarray,
    *,
    linefree_mask: np.ndarray,
    bcfg: BaselineConfig = BaselineConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    return_qc: bool = True,
) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError("cube_data must be 3D (nchan,ny,nx)")
    nchan, ny, nx = data.shape
    m3 = np.asarray(linefree_mask, dtype=bool)
    if m3.shape != (nchan, ny, nx):
        raise ValueError(f"linefree_mask shape mismatch: {m3.shape} vs ({nchan},{ny},{nx})")

    dtype_compute = _resolve_compute_dtype(bcfg)
    normalize_x = bool(getattr(bcfg, "normalize_x", False) or getattr(bcfg, "reproducible_mode", False))
    robust_mode = str(getattr(bcfg, "robust_mode", "full") or "full").strip().lower()
    if robust_mode not in {"full", "batched_full", "selective"}:
        raise ValueError(f"Unknown robust_mode={robust_mode!r}; expected 'full', 'batched_full' or 'selective'")
    robust_iters = max(0, int(getattr(bcfg, "robust_iters", 3)))
    robust_early_stop = bool(getattr(bcfg, "robust_early_stop", False))
    robust_coef_rtol = float(getattr(bcfg, "robust_coef_rtol", 1.0e-4))
    robust_coef_atol = float(getattr(bcfg, "robust_coef_atol", 1.0e-6))
    robust_selective_sigma = float(getattr(bcfg, "robust_selective_sigma", 4.5))
    robust_selective_frac = float(getattr(bcfg, "robust_selective_frac", 0.10))
    robust_selective_max_pixels = int(getattr(bcfg, "robust_selective_max_pixels", 2048))

    finite_orig = np.isfinite(data)
    data2 = _fill_nan_with_median_along_spec(data)

    npix = ny * nx
    Y_filled = data2.reshape(nchan, npix)
    Y_orig = data.reshape(nchan, npix)
    finite_flat = finite_orig.reshape(nchan, npix)
    mask_flat = m3.reshape(nchan, npix)

    freqs = list(ripple_freqs) if (bcfg.ripple and ripple_freqs) else []
    A_full = _design_matrix(
        nchan,
        poly_order=bcfg.poly_order,
        ripple_freqs=freqs,
        dtype=dtype_compute,
        normalize_poly_x=normalize_x,
    )
    ncoef = int(A_full.shape[1])
    min_points = max(ncoef + 2, 8)

    rcond_eff = _resolve_rcond(
        bcfg.rcond,
        n_rows=nchan,
        n_cols=ncoef,
        dtype=dtype_compute,
        reproducible_mode=bool(getattr(bcfg, "reproducible_mode", False)),
    )

    out = np.empty_like(Y_orig, dtype=np.float32)
    resid_rms = np.full(npix, np.nan, dtype=np.float32) if return_qc else None
    flag = np.zeros(npix, dtype=np.uint8) if return_qc else None

    chunk = int(max(1, bcfg.chunk_pix))
    reg_eps = float(np.finfo(dtype_compute).eps * max(nchan, ncoef))
    voxel_solver = str(getattr(bcfg, "voxel_solver", "qr") or "qr").strip().lower()
    if voxel_solver not in {"qr", "normal"}:
        raise ValueError(f"Unknown voxel_solver={voxel_solver!r}; expected 'qr' or 'normal'")
    voxel_qr_batch_pix = int(getattr(bcfg, "voxel_qr_batch_pix", 2048))
    for p0 in range(0, npix, chunk):
        p1 = min(npix, p0 + chunk)
        Yc_fill = np.asarray(Y_filled[:, p0:p1], dtype=dtype_compute)
        Yc_orig = Y_orig[:, p0:p1]
        finite_c = finite_flat[:, p0:p1]
        fit_mask_c = mask_flat[:, p0:p1] & finite_c
        k = p1 - p0

        out_chunk = Yc_orig.copy()
        flag_chunk = np.zeros(k, dtype=np.uint8)

        finite_lf_counts = np.sum(fit_mask_c, axis=0)
        all_nan = np.sum(finite_c, axis=0) == 0
        insufficient = (~all_nan) & (finite_lf_counts < min_points)
        fit_cols = np.where((~all_nan) & (~insufficient))[0]

        flag_chunk[insufficient] = 1
        flag_chunk[all_nan] = 3

        if fit_cols.size > 0:
            try:
                Yfit = np.asarray(Yc_fill[:, fit_cols], dtype=dtype_compute)
                Wfit = np.asarray(fit_mask_c[:, fit_cols], dtype=dtype_compute)
                coef = _solve_voxelmask_grouped_initial(
                    A_full,
                    Yfit,
                    fit_mask_c[:, fit_cols],
                    dtype=dtype_compute,
                    solver=voxel_solver,
                    rcond=rcond_eff,
                    reg_eps=reg_eps,
                    min_points=min_points,
                    qr_batch_pix=voxel_qr_batch_pix,
                    grouping=str(getattr(bcfg, "voxel_grouping", "auto")),
                    group_min_size=int(getattr(bcfg, "voxel_group_min_size", 12)),
                )
                base = A_full @ coef
                out_fit = np.asarray(Yfit - base, dtype=np.float32)
                out_chunk[:, fit_cols] = out_fit

                if bool(bcfg.robust):
                    if robust_mode in {"full", "batched_full"}:
                        selected_rel = np.arange(fit_cols.size, dtype=np.int64)
                    else:
                        selected_rel = _selective_robust_candidates(
                            out_fit,
                            fit_mask_c[:, fit_cols],
                            sigma=robust_selective_sigma,
                            max_frac=robust_selective_frac,
                            max_pixels=robust_selective_max_pixels,
                        )
                    if selected_rel.size > 0:
                        _fit_selected_robust_voxelmask(
                            out_chunk,
                            Yc_fill,
                            fit_cols,
                            selected_rel,
                            fit_mask_chunk=fit_mask_c,
                            A_full=A_full,
                            coef_init=coef[:, selected_rel] if coef.ndim == 2 else None,
                            rcond=rcond_eff,
                            dtype=dtype_compute,
                            robust_iters=robust_iters,
                            robust_early_stop=robust_early_stop,
                            robust_coef_rtol=robust_coef_rtol,
                            robust_coef_atol=robust_coef_atol,
                            batch_pixels=int(getattr(bcfg, "robust_batch_pixels", 512)),
                            solver=voxel_solver,
                            reg_eps=reg_eps,
                            qr_batch_pix=voxel_qr_batch_pix,
                        )
            except Exception:
                if bool(getattr(bcfg, "strict_failures", False)):
                    raise
                if not bool(getattr(bcfg, "fallback_to_pixelwise", True)):
                    flag_chunk[fit_cols] = 2
                else:
                    for j in fit_cols:
                        try:
                            out_chunk[:, j] = _fit_single_spectrum_masked(
                                Yc_fill[:, j],
                                A_full,
                                fit_mask_c[:, j],
                                robust=bool(bcfg.robust),
                                robust_iters=robust_iters,
                                rcond=rcond_eff,
                                dtype=dtype_compute,
                                robust_early_stop=robust_early_stop,
                                robust_coef_rtol=robust_coef_rtol,
                                robust_coef_atol=robust_coef_atol,
                            )
                        except Exception:
                            if bool(getattr(bcfg, "strict_failures", False)):
                                raise
                            flag_chunk[j] = 2

        out_chunk[~finite_c] = np.nan
        out[:, p0:p1] = out_chunk

        if return_qc:
            ok = (flag_chunk == 0)
            if np.any(ok):
                lf_finite = fit_mask_c & ok[None, :]
                r = out_chunk ** 2
                r[~lf_finite] = 0.0
                cnt = np.sum(lf_finite, axis=0)
                good = cnt > 0
                rms_local = np.full(k, np.nan, dtype=np.float32)
                if np.any(good):
                    rms_local[good] = np.sqrt(np.sum(r[:, good], axis=0) / cnt[good]).astype(np.float32, copy=False)
                resid_rms[p0:p1] = rms_local
            flag[p0:p1] = flag_chunk

    out_cube = out.reshape(nchan, ny, nx)
    if return_qc:
        return out_cube, resid_rms.reshape(ny, nx), flag.reshape(ny, nx)
    return out_cube, None, None


def subtract_baseline_cube(
    cube_data: np.ndarray,
    *,
    linefree_mask: np.ndarray,
    bcfg: BaselineConfig = BaselineConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    return_qc: bool = True,
) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Subtract baseline from a cube.

    Parameters
    ----------
    cube_data : np.ndarray
        (nchan, ny, nx).
    linefree_mask : np.ndarray
        (nchan,) bool. True=use for fitting.
    bcfg : BaselineConfig
        Model configuration.
    ripple_freqs : sequence of float or None
        Ripple frequencies (cycles/channel). If bcfg.ripple=True and ripple_freqs is None, will fit polynomial only.
    return_qc : bool
        If True, return (resid_rms_map, fit_flag_map).

    Returns
    -------
    out_cube : np.ndarray
        Baseline-subtracted cube (float32).
    resid_rms : np.ndarray or None
        2D RMS map on line-free channels (float32).
    flag : np.ndarray or None
        2D uint8 flag map: 0=OK, 1=insufficient finite line-free points, 2=lstsq failed, 3=all-NaN spectrum.
    """
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError("cube_data must be 3D (nchan,ny,nx)")
    nchan, ny, nx = data.shape
    m = _normalize_linefree_mask_shape(np.asarray(linefree_mask, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="linefree_mask")
    if m.ndim == 3:
        return _subtract_baseline_cube_voxelmask(
            data,
            linefree_mask=m,
            bcfg=bcfg,
            ripple_freqs=ripple_freqs,
            return_qc=return_qc,
        )

    dtype_compute = _resolve_compute_dtype(bcfg)
    normalize_x = bool(getattr(bcfg, "normalize_x", False) or getattr(bcfg, "reproducible_mode", False))
    robust_mode = str(getattr(bcfg, "robust_mode", "full") or "full").strip().lower()
    if robust_mode not in {"full", "batched_full", "selective"}:
        raise ValueError(f"Unknown robust_mode={robust_mode!r}; expected 'full', 'batched_full' or 'selective'")
    robust_iters = max(0, int(getattr(bcfg, "robust_iters", 3)))
    robust_early_stop = bool(getattr(bcfg, "robust_early_stop", False))
    robust_coef_rtol = float(getattr(bcfg, "robust_coef_rtol", 1.0e-4))
    robust_coef_atol = float(getattr(bcfg, "robust_coef_atol", 1.0e-6))
    robust_selective_sigma = float(getattr(bcfg, "robust_selective_sigma", 4.5))
    robust_selective_frac = float(getattr(bcfg, "robust_selective_frac", 0.10))
    robust_selective_max_pixels = int(getattr(bcfg, "robust_selective_max_pixels", 2048))
    robust_batch_pixels = int(getattr(bcfg, "robust_batch_pixels", 512))

    finite_orig = np.isfinite(data)
    data2 = _fill_nan_with_median_along_spec(data)

    npix = ny * nx
    Y_filled = data2.reshape(nchan, npix)
    Y_orig = data.reshape(nchan, npix)
    finite_flat = finite_orig.reshape(nchan, npix)

    freqs = list(ripple_freqs) if (bcfg.ripple and ripple_freqs) else []
    A_full = _design_matrix(
        nchan,
        poly_order=bcfg.poly_order,
        ripple_freqs=freqs,
        dtype=dtype_compute,
        normalize_poly_x=normalize_x,
    )
    A_lf = A_full[m, :]
    nlf_global = int(m.sum())
    ncoef = int(A_full.shape[1])
    min_points = max(ncoef + 2, 8)

    if nlf_global < min_points:
        raise ValueError(f"Too few line-free points for model: nlf={nlf_global}, ncoef={ncoef}")

    rcond_eff = _resolve_rcond(
        bcfg.rcond,
        n_rows=A_lf.shape[0],
        n_cols=A_lf.shape[1],
        dtype=dtype_compute,
        reproducible_mode=bool(getattr(bcfg, "reproducible_mode", False)),
    )

    out = np.empty_like(Y_orig, dtype=np.float32)
    resid_rms = np.full(npix, np.nan, dtype=np.float32) if return_qc else None
    flag = np.zeros(npix, dtype=np.uint8) if return_qc else None

    chunk = int(max(1, bcfg.chunk_pix))
    for p0 in range(0, npix, chunk):
        p1 = min(npix, p0 + chunk)
        Yc_fill = Y_filled[:, p0:p1]
        Yc_orig = Y_orig[:, p0:p1]
        finite_c = finite_flat[:, p0:p1]
        k = p1 - p0

        out_chunk = Yc_orig.copy()
        flag_chunk = np.zeros(k, dtype=np.uint8)

        finite_lf_counts = np.sum(finite_c[m, :], axis=0)
        all_nan = np.sum(finite_c, axis=0) == 0
        insufficient = (~all_nan) & (finite_lf_counts < min_points)
        fit_mask = (~all_nan) & (~insufficient)

        flag_chunk[insufficient] = 1
        flag_chunk[all_nan] = 3

        fit_cols = np.where(fit_mask)[0]
        if fit_cols.size > 0:
            Yc_fit = np.asarray(Yc_fill[:, fit_cols], dtype=dtype_compute)
            Ylf = Yc_fit[m, :]
            try:
                coef, *_ = np.linalg.lstsq(A_lf, Ylf, rcond=rcond_eff)
                coef = np.asarray(coef, dtype=dtype_compute)
                base = A_full @ coef
                out_fit = np.asarray(Yc_fit - base, dtype=np.float32)
                out_chunk[:, fit_cols] = out_fit

                if bool(bcfg.robust):
                    if robust_mode == "full":
                        for j in fit_cols:
                            out_chunk[:, j] = _fit_single_spectrum(
                                Yc_fill[:, j],
                                A_full,
                                A_lf,
                                m,
                                robust=True,
                                robust_iters=robust_iters,
                                rcond=rcond_eff,
                                dtype=dtype_compute,
                                robust_early_stop=robust_early_stop,
                                robust_coef_rtol=robust_coef_rtol,
                                robust_coef_atol=robust_coef_atol,
                            )
                    elif robust_mode == "batched_full":
                        _fit_selected_robust(
                            out_chunk,
                            Yc_fill,
                            fit_cols,
                            np.arange(fit_cols.size, dtype=np.int64),
                            A_full=A_full,
                            A_lf=A_lf,
                            linefree_mask=m,
                            coef_init=coef,
                            rcond=rcond_eff,
                            dtype=dtype_compute,
                            robust_iters=robust_iters,
                            batch_pixels=robust_batch_pixels,
                            robust_early_stop=robust_early_stop,
                            robust_coef_rtol=robust_coef_rtol,
                            robust_coef_atol=robust_coef_atol,
                        )
                    else:
                        finite_lf = finite_c[m, :][:, fit_cols]
                        selected_rel = _selective_robust_candidates(
                            out_fit[m, :],
                            finite_lf,
                            sigma=robust_selective_sigma,
                            max_frac=robust_selective_frac,
                            max_pixels=robust_selective_max_pixels,
                        )
                        if selected_rel.size > 0:
                            _fit_selected_robust(
                                out_chunk,
                                Yc_fill,
                                fit_cols,
                                selected_rel,
                                A_full=A_full,
                                A_lf=A_lf,
                                linefree_mask=m,
                                coef_init=coef[:, selected_rel] if coef.ndim == 2 else None,
                                rcond=rcond_eff,
                                dtype=dtype_compute,
                                robust_iters=robust_iters,
                                batch_pixels=robust_batch_pixels,
                                robust_early_stop=robust_early_stop,
                                robust_coef_rtol=robust_coef_rtol,
                                robust_coef_atol=robust_coef_atol,
                            )
            except Exception:
                if bool(getattr(bcfg, "strict_failures", False)):
                    raise
                if not bool(getattr(bcfg, "fallback_to_pixelwise", True)):
                    flag_chunk[fit_cols] = 2
                else:
                    ok_rel: list[int] = []
                    for jrel, j in enumerate(fit_cols):
                        try:
                            out_chunk[:, j] = _fit_single_spectrum(
                                Yc_fill[:, j],
                                A_full,
                                A_lf,
                                m,
                                robust=bool(bcfg.robust and robust_mode == "full"),
                                robust_iters=robust_iters,
                                rcond=rcond_eff,
                                dtype=dtype_compute,
                                robust_early_stop=robust_early_stop,
                                robust_coef_rtol=robust_coef_rtol,
                                robust_coef_atol=robust_coef_atol,
                            )
                            ok_rel.append(int(jrel))
                        except Exception:
                            if bool(getattr(bcfg, "strict_failures", False)):
                                raise
                            flag_chunk[j] = 2
                    if bool(bcfg.robust) and len(ok_rel) > 0:
                        ok_rel_arr = np.asarray(ok_rel, dtype=np.int64)
                        fit_cols_ok = fit_cols[ok_rel_arr]
                        if robust_mode == "batched_full":
                            _fit_selected_robust(
                                out_chunk,
                                Yc_fill,
                                fit_cols_ok,
                                np.arange(fit_cols_ok.size, dtype=np.int64),
                                A_full=A_full,
                                A_lf=A_lf,
                                linefree_mask=m,
                                coef_init=None,
                                rcond=rcond_eff,
                                dtype=dtype_compute,
                                robust_iters=robust_iters,
                                batch_pixels=robust_batch_pixels,
                                robust_early_stop=robust_early_stop,
                                robust_coef_rtol=robust_coef_rtol,
                                robust_coef_atol=robust_coef_atol,
                            )
                        elif robust_mode == "selective":
                            finite_lf = finite_c[m, :][:, fit_cols_ok]
                            selected_rel = _selective_robust_candidates(
                                out_chunk[m, :][:, fit_cols_ok],
                                finite_lf,
                                sigma=robust_selective_sigma,
                                max_frac=robust_selective_frac,
                                max_pixels=robust_selective_max_pixels,
                            )
                            if selected_rel.size > 0:
                                _fit_selected_robust(
                                    out_chunk,
                                    Yc_fill,
                                    fit_cols_ok,
                                    selected_rel,
                                    A_full=A_full,
                                    A_lf=A_lf,
                                    linefree_mask=m,
                                    coef_init=None,
                                    rcond=rcond_eff,
                                    dtype=dtype_compute,
                                    robust_iters=robust_iters,
                                    batch_pixels=robust_batch_pixels,
                                    robust_early_stop=robust_early_stop,
                                    robust_coef_rtol=robust_coef_rtol,
                                    robust_coef_atol=robust_coef_atol,
                                )

        # Preserve original invalid channels rather than inventing data on NaN inputs.
        out_chunk[~finite_c] = np.nan
        out[:, p0:p1] = out_chunk

        if return_qc:
            ok = (flag_chunk == 0)
            if np.any(ok):
                lf_finite = finite_c[m, :] & ok[None, :]
                r = out_chunk[m, :] ** 2
                r[~lf_finite] = 0.0
                cnt = np.sum(lf_finite, axis=0)
                good = cnt > 0
                rms_local = np.full(k, np.nan, dtype=np.float32)
                if np.any(good):
                    rms_local[good] = np.sqrt(np.sum(r[:, good], axis=0) / cnt[good]).astype(np.float32, copy=False)
                resid_rms[p0:p1] = rms_local
            flag[p0:p1] = flag_chunk

    out_cube = out.reshape(nchan, ny, nx)
    if return_qc:
        return out_cube, resid_rms.reshape(ny, nx), flag.reshape(ny, nx)
    return out_cube, None, None


# -----------------------------------------------------------------------------
# 4) FITS wrapper (chunked write)
# -----------------------------------------------------------------------------

def _get_cube_hdu(hdul: fits.HDUList, cube_ext: Optional[Union[int, str]]) -> Tuple[int, fits.ImageHDU]:
    image_types = (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)
    excluded_names = {
        "MASK3D", "BASESUP3D", "LINECAND3D",
        "RMS", "BASE_RMS", "WEIGHT", "HIT", "MASK", "TSYS", "TINT", "TIME",
        "MOMENT0", "MOM0_BASESUP", "MOM0_LINECAND", "SIGNAL_MASK_USED",
        "DATA_SIGNAL_ONLY", "DATA_LINEFREE_ONLY",
    }
    excluded_btypes = {
        "SIGNALMASK", "BASELINESUPPORTMASK", "LINECANDIDATEMASK",
        "MOMENT0", "MOMENT0BASELINESUPPORT", "MOMENT0LINECANDIDATE", "SIGNALMASKUSED",
        "DATASIGNALONLY", "DATALINEFREEONLY",
        "WEIGHT", "HITCOUNT", "VALIDMASK", "SYSTEMTEMP", "INTEGRATIONTIME", "OBSERVATIONTIME",
        "BASELINERESIDRMS", "BASELINERMS",
    }

    def _is_3d_image(hdu: object) -> bool:
        return isinstance(hdu, image_types) and getattr(hdu, "data", None) is not None and np.ndim(hdu.data) == 3

    def _looks_like_analysis_hdu(hdu: object) -> bool:
        name = str(getattr(hdu, "name", "") or "").upper()
        hdr = getattr(hdu, "header", None)
        btype = str(hdr.get("BTYPE", "") if hdr is not None else "").upper()
        return (name in excluded_names) or (btype in excluded_btypes)

    def _has_spectral_axis(hdu: object) -> bool:
        hdr = getattr(hdu, "header", None)
        if hdr is None:
            return False
        return _find_spectral_fits_axis(hdr) is not None

    if cube_ext is None:
        if _is_3d_image(hdul[0]) and (not _looks_like_analysis_hdu(hdul[0])) and _has_spectral_axis(hdul[0]):
            return 0, hdul[0]  # type: ignore
        for i, hdu in enumerate(hdul):
            if _is_3d_image(hdu) and (not _looks_like_analysis_hdu(hdu)) and _has_spectral_axis(hdu):
                return i, hdu  # type: ignore
        raise ValueError(
            "No suitable 3D spectral cube found in FITS. Refusing to fall back to a 3D HDU without spectral CTYPE."
        )

    hdu = hdul[cube_ext]
    if not _is_3d_image(hdu):
        shape = getattr(getattr(hdu, "data", None), "shape", None)
        raise ValueError(f"Selected cube_ext={cube_ext!r} is not a 3D image HDU (shape={shape}).")
    if _looks_like_analysis_hdu(hdu):
        raise ValueError(f"Selected cube_ext={cube_ext!r} points to an analysis/product HDU, not an input data cube.")
    if not _has_spectral_axis(hdu):
        raise ValueError(
            f"Selected cube_ext={cube_ext!r} does not advertise a spectral axis in CTYPE1..3; refusing to guess."
        )
    for i, hh in enumerate(hdul):
        if hh is hdu:
            return i, hdu  # type: ignore
    raise RuntimeError("Could not resolve cube_ext to index.")


def _find_spectral_fits_axis(header: fits.Header) -> Optional[int]:
    """Return the FITS axis number (1-based) that appears to be spectral."""
    for ax in (1, 2, 3):
        ctype = str(header.get(f"CTYPE{ax}", "")).upper()
        if any(tok in ctype for tok in ("FREQ", "VRAD", "VELO", "VOPT", "WAVE", "AWAV")):
            return ax
    return None


def _standardize_cube_for_processing(data: np.ndarray, header: fits.Header) -> Tuple[np.ndarray, str]:
    """
    Standardize FITS cube data to (nchan, ny, nx) using the header spectral axis.
    Fail loudly if the spectral axis cannot be identified from CTYPE1..3.

    Returns
    -------
    data_std : np.ndarray
        Cube with spectral axis moved to axis 0.
    axis_order_in : str
        One of 'v_y_x', 'y_v_x', 'y_x_v'.
    """
    arr = np.asarray(data)
    if arr.ndim != 3:
        raise ValueError(f"Cube data must be 3D, got shape={arr.shape}")

    fits_ax = _find_spectral_fits_axis(header)
    if fits_ax is None:
        raise ValueError(
            "Could not identify the spectral FITS axis from CTYPE1..3. "
            "Refusing to guess axis order from shape alone."
        )

    np_ax = arr.ndim - int(fits_ax)
    if np_ax == 0:
        return arr, "v_y_x"
    if np_ax == 1:
        return np.transpose(arr, (1, 0, 2)), "y_v_x"
    if np_ax == 2:
        return np.transpose(arr, (2, 0, 1)), "y_x_v"
    raise ValueError(f"Could not map spectral FITS axis {fits_ax} to numpy axis for shape={arr.shape}")


def _restore_cube_axis_order(data_std: np.ndarray, axis_order_in: str) -> np.ndarray:
    """Restore a standardized (nchan, ny, nx) cube to its original axis order."""
    arr = np.asarray(data_std)
    if axis_order_in == "v_y_x":
        return arr
    if axis_order_in == "y_v_x":
        return np.transpose(arr, (1, 0, 2))
    if axis_order_in == "y_x_v":
        return np.transpose(arr, (1, 2, 0))
    raise ValueError(f"Unknown axis_order_in: {axis_order_in}")


def _replace_or_append_hdu(
    hdul: fits.HDUList,
    hdu: Union[fits.ImageHDU, fits.CompImageHDU, fits.BinTableHDU],
) -> None:
    """Replace HDU with the same EXTNAME if present, otherwise append."""
    name = str(getattr(hdu, "name", "")).upper()
    if not name:
        hdul.append(hdu)
        return
    while name in hdul:
        del hdul[name]
    hdul.append(hdu)


def _remove_named_hdus(hdul: fits.HDUList, names: Sequence[str], *, protect_hdu: object | None = None) -> None:
    """Remove all HDUs whose EXTNAME matches one of `names` (case-insensitive)."""
    names_up = {str(n).upper() for n in names}
    remove_idx = []
    for i, hdu in enumerate(hdul):
        if i == 0:
            continue
        if protect_hdu is not None and hdu is protect_hdu:
            continue
        hname = str(getattr(hdu, 'name', '') or '').upper()
        if hname in names_up:
            remove_idx.append(i)
    for i in reversed(remove_idx):
        del hdul[i]


def _strip_checksum_all_hdus(hdul: fits.HDUList) -> None:
    """Remove CHECKSUM/DATASUM cards from all HDUs in-place."""
    for hdu in hdul:
        hdr = getattr(hdu, "header", None)
        if hdr is None:
            continue
        for key in ("CHECKSUM", "DATASUM", "ZHECKSUM", "ZDATASUM"):
            while key in hdr:
                del hdr[key]



def _resolve_exclude_v_windows(
    exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]],
) -> Optional[Sequence[Union[str, Tuple[float, float]]]]:
    """Normalize explicit velocity windows excluded from baseline fitting."""
    if exclude_v_windows is None or len(exclude_v_windows) == 0:
        return None
    return exclude_v_windows


def _resolve_linefree_velocity_windows(
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]],
) -> Optional[Sequence[Union[str, Tuple[float, float]]]]:
    if linefree_velocity_windows_kms is None or len(linefree_velocity_windows_kms) == 0:
        return None
    return linefree_velocity_windows_kms


def _normalize_linefree_mode(linefree_mode: Optional[str]) -> str:
    mode = "infer" if linefree_mode is None else str(linefree_mode).strip().lower()
    allowed = {"infer", "manual", "prior", "auto", "or"}
    if mode not in allowed:
        raise ValueError(f"linefree_mode must be one of {sorted(allowed)}, got {linefree_mode!r}.")
    return mode


def _validate_linefree_selection_inputs(
    *,
    linefree_mask: Optional[np.ndarray],
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]],
    linefree_mode: Optional[str],
) -> str:
    mode = _normalize_linefree_mode(linefree_mode)
    if mode == "manual":
        if linefree_velocity_windows_kms is None:
            raise ValueError("linefree_mode='manual' requires linefree_velocity_windows_kms.")
        if linefree_mask is not None:
            raise ValueError(
                "linefree_mode='manual' uses linefree_velocity_windows_kms as the exact manual selection and "
                "cannot be combined with linefree_mask. Omit linefree_mode to infer manual selection from linefree_mask."
            )
    elif linefree_mask is not None and mode in ("prior", "auto", "or"):
        raise ValueError(
            "Explicit line-free mask (linefree_mask) "
            f"cannot be combined with linefree_mode={mode!r}."
        )
    return mode


def _resolve_effective_linefree_mode(
    *,
    linefree_mode: str,
    linefree_mask: Optional[np.ndarray],
    linefree_cfg: Optional[LineFreeConfig],
    has_prior: bool,
) -> str:
    if linefree_mode != "infer":
        return linefree_mode
    if linefree_mask is not None:
        return "manual_mask"
    if has_prior:
        return "prior"
    if linefree_cfg is not None:
        return "auto"
    raise ValueError(
        "Could not infer line-free selection source. Provide linefree_mask, set linefree_mode='manual' with "
        "linefree_velocity_windows_kms, provide linefree_cfg for auto detection, or supply prior LINEFREE."
    )


def _linefree_cfg_to_meta_dict(cfg: Optional[LineFreeConfig]) -> Optional[dict[str, Any]]:
    if cfg is None:
        return None
    return dict(cfg.__dict__)



def _linefree_mask_from_axis(
    axis: np.ndarray,
    windows: Sequence[Union[str, Tuple[float, float]]],
) -> np.ndarray:
    return _manual_signal_mask_from_axis(axis, windows)


def _validate_baseline_viewer_mode(mode: str) -> str:
    mode_norm = str(mode).strip().lower()
    allowed = {"signal", "linefree", "signal3d", "linefree3d"}
    if mode_norm not in allowed:
        raise ValueError(f"mode must be one of {sorted(allowed)}, got {mode!r}.")
    return mode_norm


def _true_runs_1d(mask: np.ndarray) -> list[tuple[int, int]]:
    arr = np.asarray(mask, dtype=bool).reshape(-1)
    if arr.size == 0:
        return []
    padded = np.concatenate(([False], arr, [False]))
    diff = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(diff == 1)
    ends = np.flatnonzero(diff == -1)
    return [(int(s), int(e)) for s, e in zip(starts, ends)]


def _signal_mask_from_linefree_internal_gaps(linefree_mask_used: np.ndarray) -> np.ndarray:
    lf = np.asarray(linefree_mask_used, dtype=bool)
    if lf.ndim == 1:
        flat = lf.reshape(lf.shape[0], 1)
        reshape_to = None
    elif lf.ndim == 3:
        flat = lf.reshape(lf.shape[0], -1)
        reshape_to = lf.shape
    else:
        raise ValueError(f"linefree_mask_used must be 1D or 3D, got shape={lf.shape}")
    if flat.shape[0] == 0:
        sig_flat = np.zeros_like(flat, dtype=bool)
    else:
        has_before = np.cumsum(flat.astype(np.int32), axis=0) > 0
        has_after = np.cumsum(flat[::-1].astype(np.int32), axis=0)[::-1] > 0
        sig_flat = (~flat) & has_before & has_after
    if reshape_to is None:
        return sig_flat[:, 0]
    return sig_flat.reshape(reshape_to)


def _moment0_from_signal_mask(data_cube: np.ndarray, signal_mask_used: np.ndarray, velocity_axis_kms: np.ndarray) -> np.ndarray:
    data = np.asarray(data_cube, dtype=np.float32)
    sig = np.asarray(signal_mask_used, dtype=bool)
    v = np.asarray(velocity_axis_kms, dtype=float).reshape(-1)
    if data.ndim != 3:
        raise ValueError(f"data_cube must be 3D (nchan, ny, nx), got shape={data.shape}")
    if v.shape[0] != data.shape[0]:
        raise ValueError(
            f"velocity_axis_kms length mismatch: {v.shape[0]} vs nchan={data.shape[0]}"
        )
    if sig.ndim == 1:
        if sig.shape[0] != data.shape[0]:
            raise ValueError(
                f"signal_mask_used length mismatch: {sig.shape[0]} vs nchan={data.shape[0]}"
            )
        sig3 = sig[:, None, None]
    elif sig.ndim == 3:
        if sig.shape != data.shape:
            raise ValueError(
                f"signal_mask_used shape mismatch: {sig.shape} vs data_cube={data.shape}"
            )
        sig3 = sig
    else:
        raise ValueError(f"signal_mask_used must be 1D or 3D, got shape={sig.shape}")
    if data.shape[0] <= 1:
        dv = 0.0
    else:
        dv = float(np.nanmedian(np.abs(np.diff(v))))
        if not np.isfinite(dv):
            dv = 0.0
    if np.count_nonzero(sig3) == 0 or dv <= 0.0:
        return np.zeros(data.shape[1:], dtype=np.float32)
    cube_sel = np.where(sig3, data, 0.0)
    with np.errstate(invalid="ignore"):
        mom0 = np.nansum(cube_sel, axis=0, dtype=np.float64) * dv
    return np.asarray(mom0, dtype=np.float32)


def _viewer_cube_from_mask(data_cube: np.ndarray, mask: np.ndarray) -> np.ndarray:
    data = np.asarray(data_cube, dtype=np.float32)
    m = np.asarray(mask, dtype=bool)
    if data.ndim != 3:
        raise ValueError(f"data_cube must be 3D (nchan, ny, nx), got shape={data.shape}")
    if m.ndim == 1:
        if m.shape[0] != data.shape[0]:
            raise ValueError(f"mask length mismatch: {m.shape[0]} vs nchan={data.shape[0]}")
        m3 = m[:, None, None]
    elif m.ndim == 3:
        if m.shape != data.shape:
            raise ValueError(f"mask shape mismatch: {m.shape} vs data_cube={data.shape}")
        m3 = m
    else:
        raise ValueError(f"mask must be 1D or 3D, got shape={m.shape}")
    out = np.zeros_like(data, dtype=np.float32)
    if np.any(m3):
        out[m3] = data[m3]
    return out


def _derive_signal_products(
    data_cube: np.ndarray,
    *,
    linefree_mask_used: np.ndarray,
    velocity_axis_kms: np.ndarray | None,
) -> dict[str, np.ndarray]:
    lf = np.asarray(linefree_mask_used, dtype=bool)
    sig = _signal_mask_from_linefree_internal_gaps(lf)
    if lf.ndim == 1:
        sig_name = "SIGNAL_MASK_USED"
    elif lf.ndim == 3:
        sig_name = "SIGNAL_MASK3D_USED"
    else:
        raise ValueError(f"linefree_mask_used must be 1D or 3D, got shape={lf.shape}")
    products: dict[str, np.ndarray] = {
        sig_name: sig.astype(bool),
    }
    if velocity_axis_kms is not None:
        try:
            products["MOMENT0"] = _moment0_from_signal_mask(data_cube, sig, velocity_axis_kms)
        except Exception as exc:
            logging.warning("Could not compute MOMENT0 from %s: %s", sig_name, exc)
    return products

def _subtract_baseline_from_fits_legacy(
    input_fits: str,
    output_fits: str,
    *,
    cube_ext: Optional[Union[int, str]] = None,
    # line-free
    linefree_cfg: Optional[LineFreeConfig] = None,
    linefree_mask: Optional[np.ndarray] = None,
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    linefree_mode: Optional[str] = "infer",
    load_prior_from_input: bool = True,
    # ripple
    ripple_cfg: RippleConfig = RippleConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    ripple_mode: str = "auto",
    # baseline model
    baseline_cfg: BaselineConfig = BaselineConfig(),
    reproducible_mode: Optional[bool] = None,
    # output
    add_qc_hdus: bool = True,
    overwrite: bool = True,
    write_diagnostics: bool = False,
    diagnostics_prefix: Optional[str] = None,
    write_profile: bool = False,
    profile_prefix: Optional[str] = None,
) -> None:
    """
    Read cube FITS, estimate/reuse line-free + ripple frequencies, subtract baseline, write new FITS.

    Outputs (if add_qc_hdus=True)
    ----------------------------
    - LINEFREE (ImageHDU): uint8 (nchan,)  1=line-free, 0=line
    - RIPFREQ  (BinTable): ripple frequencies and periods (if ripple used)
    - BASE_RMS (ImageHDU): 2D residual RMS on line-free channels
    - BASE_FLG (ImageHDU): 2D flag map (0 ok, 1 insufficient finite line-free points, 2 fit failed, 3 all-NaN spectrum)

    Notes
    -----
    - Baseline subtraction is done in channel space.
    - If load_prior_from_input=True, existing LINEFREE / RIPFREQ HDUs are used as priors
      according to linefree_mode / ripple_mode.
    """
    _profile_records: List[dict[str, Any]] = []
    _t0 = time.perf_counter()
    _tprev = _t0
    if write_profile:
        _tprev = _profile_record(_profile_records, "start", _tprev, _t0, input_fits=str(input_fits), output_fits=str(output_fits), cube_ext=(None if cube_ext is None else str(cube_ext)))
    def _read_linefree_prior(hdul_local: fits.HDUList) -> Optional[np.ndarray]:
        for name in ("LINEFREE3D_USED", "LINEFREE3D", "LINEFREE3D_PRIOR", "LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR"):
            if name in hdul_local:
                arr = np.squeeze(np.asarray(hdul_local[name].data))
                if arr.ndim in (1, 3):
                    return arr.astype(bool)
        return None

    def _read_ripple_prior(hdul_local: fits.HDUList) -> Optional[List[float]]:
        for name in ("RIPFREQ", "RIPFREQ_USED", "RIPFREQ_PRIOR"):
            if name not in hdul_local:
                continue
            tbl = hdul_local[name].data
            if tbl is None:
                continue
            names = getattr(tbl, "names", [])
            if "FREQ_CYC_PER_CH" in names:
                arr = np.asarray(tbl["FREQ_CYC_PER_CH"], dtype=float)
                if arr.size > 0:
                    return [float(v) for v in arr]
        return None

    exclude_v_windows = _resolve_exclude_v_windows(exclude_v_windows)
    linefree_velocity_windows_kms = _resolve_linefree_velocity_windows(linefree_velocity_windows_kms)
    safe_velocity_windows_kms = _resolve_linefree_velocity_windows(safe_velocity_windows_kms)
    ripple_apply_stage = _normalize_ripple_apply_stage(ripple_apply_stage)
    linefree_mode = _validate_linefree_selection_inputs(
        linefree_mask=linefree_mask,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        linefree_mode=linefree_mode,
    )

    with fits.open(input_fits, mode="readonly", memmap=True) as hdul_in:
        idx, hdu_cube = _get_cube_hdu(hdul_in, cube_ext)
        data_raw = np.asarray(hdu_cube.data, dtype=np.float32)
        data, axis_order_in = _standardize_cube_for_processing(data_raw, hdu_cube.header)
        nchan, ny, nx = data.shape
        if write_profile:
            _tprev = _profile_record(_profile_records, "read_and_standardize", _tprev, _t0, cube_shape=[int(nchan), int(ny), int(nx)], axis_order_in=str(axis_order_in))

        lf_prior = _read_linefree_prior(hdul_in) if load_prior_from_input else None
        freqs_prior = _read_ripple_prior(hdul_in) if load_prior_from_input else None

        if lf_prior is not None:
            try:
                lf_prior = _normalize_linefree_mask_shape(np.asarray(lf_prior, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="LINEFREE prior")
            except Exception:
                logging.warning(
                    "Ignoring prior LINEFREE because shape=%s does not match (%d,) or (%d,%d,%d).",
                    np.asarray(lf_prior).shape,
                    nchan,
                    nchan,
                    ny,
                    nx,
                )
                lf_prior = None

        v_axis = None
        need_velocity_axis = bool(
            linefree_velocity_windows_kms
            or exclude_v_windows
            or safe_velocity_windows_kms
            or getattr(ripple_cfg, "period_range_kms", None) is not None
            or _linefree_cfg_needs_velocity_axis(linefree_cfg)
        )
        if need_velocity_axis:
            v_axis = _velocity_axis_from_fits_context(
                input_fits,
                hdu_cube.header,
                cube_ext=cube_ext,
                hdu_index=idx,
                nchan=nchan,
            )
            if v_axis.shape != (nchan,):
                raise ValueError(f"velocity axis shape mismatch: {v_axis.shape} vs ({nchan},)")
        spectral_axis_hz = None
        if getattr(ripple_cfg, "period_range_hz", None) is not None:
            spectral_axis_hz = _spectral_axis_hz_from_fits_context(
                input_fits,
                hdu_cube.header,
                cube_ext=cube_ext,
                hdu_index=idx,
                nchan=nchan,
            )
            if spectral_axis_hz.shape != (nchan,):
                raise ValueError(f"spectral axis shape mismatch: {spectral_axis_hz.shape} vs ({nchan},)")

        cell_arcsec = _header_cell_arcsec(hdu_cube.header) if _linefree_cfg_needs_cell_arcsec(linefree_cfg) else None

        data_for_linefree = np.asarray(data, dtype=np.float32)
        detection_ripple_freqs: List[float] = []
        if linefree_detection_data is not None:
            arr_det = np.asarray(linefree_detection_data, dtype=np.float32)
            if arr_det.shape != data.shape:
                raise ValueError(f"linefree_detection_data shape mismatch: {arr_det.shape} vs {data.shape}")
            data_for_linefree = arr_det
        elif safe_velocity_windows_kms and ripple_apply_stage in {"mask_only", "both"} and bool(baseline_cfg.ripple):
            if v_axis is None:
                raise ValueError("safe_velocity_windows_kms requires a velocity axis.")
            safe_mask_1d = _linefree_mask_from_axis(np.asarray(v_axis, dtype=float), safe_velocity_windows_kms)
            if np.any(safe_mask_1d):
                if ripple_freqs is not None:
                    detection_ripple_freqs = [float(f) for f in ripple_freqs]
                elif ripple_mode == "prior" and freqs_prior:
                    detection_ripple_freqs = [float(f) for f in freqs_prior]
                else:
                    flat_det = data.reshape(nchan, ny * nx)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        spec_det = np.nanmedian(flat_det, axis=1)
                    detection_ripple_freqs = estimate_ripple_frequencies_fft(
                        spec_det,
                        safe_mask_1d,
                        rcfg=ripple_cfg,
                        poly_order_pre=baseline_cfg.poly_order,
                        velocity_axis_kms=v_axis,
                        spectral_axis_hz=spectral_axis_hz,
                    )
                data_for_linefree = _fit_detection_cube_with_safe_windows(
                    data,
                    safe_mask_1d=safe_mask_1d,
                    baseline_cfg=baseline_cfg,
                    ripple_freqs=detection_ripple_freqs,
                )

        effective_linefree_mode = _resolve_effective_linefree_mode(
            linefree_mode=linefree_mode,
            linefree_mask=linefree_mask,
            linefree_cfg=linefree_cfg,
            has_prior=(lf_prior is not None),
        )
        if effective_linefree_mode == "manual_mask":
            lf = _normalize_linefree_mask_shape(np.asarray(linefree_mask, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="linefree_mask")
        elif effective_linefree_mode == "manual":
            if v_axis is None:
                raise ValueError("linefree_mode='manual' requires a velocity axis.")
            lf = _linefree_mask_from_axis(np.asarray(v_axis, dtype=float), linefree_velocity_windows_kms)
        else:
            lf_auto: Optional[np.ndarray] = None
            if effective_linefree_mode in ("auto", "or"):
                if linefree_cfg is None:
                    raise ValueError(f"linefree_mode={effective_linefree_mode!r} requires linefree_cfg.")
                logging.info("Estimating line-free mask from cube.")
                lf_auto = estimate_linefree_mask_from_cube(
                    data_for_linefree,
                    cfg=linefree_cfg,
                    velocity_axis_kms=v_axis,
                    cell_arcsec=cell_arcsec,
                )
                lf_auto = _normalize_linefree_mask_shape(np.asarray(lf_auto, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="auto linefree_mask")
            if effective_linefree_mode == "prior":
                if lf_prior is None:
                    raise ValueError("linefree_mode='prior' requested, but no prior LINEFREE was found in the input.")
                lf = np.asarray(lf_prior, dtype=bool)
            elif effective_linefree_mode == "or":
                if lf_prior is None and lf_auto is None:
                    raise ValueError("linefree_mode='or' requested, but neither prior LINEFREE nor auto linefree_cfg is available.")
                lf = _merge_linefree_masks_or(lf_prior, lf_auto, nchan=nchan, ny=ny, nx=nx)
            elif effective_linefree_mode == "auto":
                if lf_auto is None:
                    raise ValueError("linefree_mode='auto' requested, but auto detection did not produce a mask.")
                lf = np.asarray(lf_auto, dtype=bool)
            else:
                raise ValueError(f"Unknown effective linefree_mode: {effective_linefree_mode}")
        lf = _apply_linefree_include_exclude(
            np.asarray(lf, dtype=bool),
            velocity_axis_kms=v_axis,
            linefree_velocity_windows_kms=(None if effective_linefree_mode == "manual" else linefree_velocity_windows_kms),
            exclude_v_windows=exclude_v_windows,
        )

        logging.info("Line-free fraction: %.3f", lf.mean())
        if write_profile:
            _tprev = _profile_record(_profile_records, "prepare_linefree", _tprev, _t0, linefree_fraction=float(np.mean(lf)))

        freqs: List[float] = []
        if reproducible_mode is not None:
            baseline_cfg = replace(
                baseline_cfg,
                reproducible_mode=bool(reproducible_mode),
                compute_dtype=("float64" if bool(reproducible_mode) else str(getattr(baseline_cfg, "compute_dtype", "float32"))),
                normalize_x=(bool(reproducible_mode) or bool(getattr(baseline_cfg, "normalize_x", False))),
            )

        if baseline_cfg.ripple and ripple_apply_stage in {"final", "both"}:
            if ripple_freqs is not None:
                freqs = [float(f) for f in ripple_freqs]
                logging.info("Using user-provided ripple frequencies: %s", freqs)
            elif ripple_mode == "prior" and freqs_prior:
                freqs = [float(f) for f in freqs_prior]
                logging.info("Using prior ripple frequencies from input FITS: %s", freqs)
            else:
                logging.info("Estimating ripple frequencies from aggregated spectrum (FFT).")
                flat = data.reshape(nchan, ny * nx)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    spec = np.nanmedian(flat, axis=1)
                ripple_mask_1d = _collapse_linefree_mask_for_ripple(np.asarray(lf, dtype=bool))
                if safe_velocity_windows_kms and v_axis is not None:
                    ripple_mask_1d = _linefree_mask_from_axis(np.asarray(v_axis, dtype=float), safe_velocity_windows_kms)
                freqs = estimate_ripple_frequencies_fft(
                    spec,
                    ripple_mask_1d,
                    rcfg=ripple_cfg,
                    poly_order_pre=baseline_cfg.poly_order,
                    velocity_axis_kms=v_axis,
                    spectral_axis_hz=spectral_axis_hz,
                )
                logging.info("Estimated ripple frequencies (cycles/chan): %s", freqs)
        if write_profile:
            _tprev = _profile_record(_profile_records, "prepare_ripple", _tprev, _t0, nfreq=int(len(freqs)))

        out_cube, rms_map, flag_map = subtract_baseline_cube(
            data,
            linefree_mask=lf,
            bcfg=baseline_cfg,
            ripple_freqs=freqs,
            return_qc=add_qc_hdus,
        )
        if write_profile:
            _tprev = _profile_record(_profile_records, "fit_baseline", _tprev, _t0)

        hdul_out = fits.HDUList([h.copy() for h in hdul_in])
        if write_profile:
            _tprev = _profile_record(_profile_records, "copy_output_hdus", _tprev, _t0, nhdu=int(len(hdul_out)))

    out_cube_write = _restore_cube_axis_order(out_cube, axis_order_in)
    if write_profile:
        _tprev = _profile_record(_profile_records, "restore_axis_order", _tprev, _t0)
    if idx == 0:
        hdul_out[0].data = out_cube_write
    else:
        hdul_out[idx].data = out_cube_write

    # Remove stale baseline/provenance and downstream analysis products from copied input FITS.
    target_cube_hdu = hdul_out[idx]

    _remove_named_hdus(
        hdul_out,
        [
            "LINEFREE",
            "LINEFREE_USED",
            "LINEFREE3D",
            "LINEFREE3D_USED",
            "LINEFREE_PRIOR",
            "LINEFREE3D_PRIOR",
            "RIPFREQ",
            "RIPFREQ_USED",
            "RIPFREQ_PRIOR",
            "BASE_RMS",
            "BASE_FLG",
            "BASEFLAG",
            "BASE_DIAG",
            "RMS",
            "MASK3D",
            "MOMENT0",
            "SIGNAL_MASK_USED",
            "SIGNAL_MASK3D_USED",
            "DATA_SIGNAL_ONLY",
            "DATA_LINEFREE_ONLY",
            "BASESUP3D",
            "LINECAND3D",
            "MOM0_BASESUP",
            "MOM0_LINECAND",
        ],
        protect_hdu=target_cube_hdu,
    )

    _strip_checksum_all_hdus(hdul_out)

    signal_mask_used = _signal_mask_from_linefree_internal_gaps(np.asarray(lf, dtype=bool))
    products = _derive_signal_products(
        np.asarray(out_cube, dtype=np.float32),
        linefree_mask_used=np.asarray(lf, dtype=bool),
        velocity_axis_kms=(v_axis if "v_axis" in locals() and v_axis is not None and np.asarray(v_axis).shape == (nchan,) else None),
    )

    sig_key = "SIGNAL_MASK3D_USED" if np.asarray(signal_mask_used).ndim == 3 else "SIGNAL_MASK_USED"
    hdr_sig = fits.Header()
    hdr_sig["BTYPE"] = "SignalMaskUsed"
    if sig_key == "SIGNAL_MASK3D_USED":
        hdr_sig["COMMENT"] = "1=voxel used for MOMENT0 from internal gaps of LINEFREE3D_USED along the spectral axis; 0=otherwise."
    else:
        hdr_sig["COMMENT"] = "1=channel used for MOMENT0 from internal gaps of LINEFREE_USED; 0=otherwise."
    _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=np.asarray(products[sig_key], dtype=np.uint8), header=hdr_sig, name=sig_key))

    if "MOMENT0" in products:
        hdr_m0 = fits.Header()
        hdr_m0["BTYPE"] = "Moment0"
        in_bunit = str(target_cube_hdu.header.get("BUNIT", "") or "").strip()
        hdr_m0["BUNIT"] = (f"{in_bunit} km/s" if in_bunit else "K km/s")
        if sig_key == "SIGNAL_MASK3D_USED":
            hdr_m0["COMMENT"] = "Integrated intensity over SIGNAL_MASK3D_USED."
        else:
            hdr_m0["COMMENT"] = "Integrated intensity over SIGNAL_MASK_USED (internal gaps of LINEFREE_USED only)."
        _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=_as_float32(products["MOMENT0"]), header=hdr_m0, name="MOMENT0"))

    if add_qc_hdus:
        hdr_lf = fits.Header()
        hdr_lf["BTYPE"] = "LineFreeMask"
        hdr_lf["COMMENT"] = "1=line-free region used for baseline fitting; 0=line/ignored."
        if np.asarray(lf).ndim == 3:
            lf_hdu = fits.ImageHDU(data=np.asarray(lf, dtype=np.uint8), header=hdr_lf, name="LINEFREE3D")
            _replace_or_append_hdu(hdul_out, lf_hdu)
            lf_used_hdu = fits.ImageHDU(data=np.asarray(lf, dtype=np.uint8), header=hdr_lf.copy(), name="LINEFREE3D_USED")
            _replace_or_append_hdu(hdul_out, lf_used_hdu)
        else:
            lf_hdu = fits.ImageHDU(data=np.asarray(lf, dtype=np.uint8), header=hdr_lf, name="LINEFREE")
            _replace_or_append_hdu(hdul_out, lf_hdu)
            lf_used_hdu = fits.ImageHDU(data=np.asarray(lf, dtype=np.uint8), header=hdr_lf.copy(), name="LINEFREE_USED")
            _replace_or_append_hdu(hdul_out, lf_used_hdu)

        if baseline_cfg.ripple and freqs:
            freq_arr = np.asarray(freqs, dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                period_arr = np.where(freq_arr != 0.0, 1.0 / freq_arr, np.inf)
            cols = [
                fits.Column(name="FREQ_CYC_PER_CH", format="D", array=freq_arr),
                fits.Column(name="PERIOD_CH", format="D", array=period_arr),
            ]
            tbhdu = fits.BinTableHDU.from_columns(cols, name="RIPFREQ")
            tbhdu.header["COMMENT"] = "Ripple frequencies used in baseline model."
            _replace_or_append_hdu(hdul_out, tbhdu)
            tbhdu_used = fits.BinTableHDU.from_columns(cols, name="RIPFREQ_USED")
            tbhdu_used.header["COMMENT"] = "Ripple frequencies used in baseline model."
            _replace_or_append_hdu(hdul_out, tbhdu_used)

        if rms_map is not None:
            hdr = fits.Header()
            hdr["BTYPE"] = "BaselineResidRMS"
            hdr["BUNIT"] = "K"
            _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=_as_float32(rms_map), header=hdr, name="BASE_RMS"))

        if flag_map is not None:
            hdr = fits.Header()
            hdr["BTYPE"] = "BaselineFitFlag"
            hdr["COMMENT"] = "0=OK, 1=insufficient finite line-free points, 2=lstsq failed (spectrum left unchanged), 3=all-NaN spectrum."
            _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=np.asarray(flag_map, dtype=np.uint8), header=hdr, name="BASE_FLG"))

    if write_profile:
        _tprev = _profile_record(_profile_records, "add_qc_hdus", _tprev, _t0)

    hdul_out.writeto(output_fits, overwrite=overwrite)
    if write_profile:
        _tprev = _profile_record(_profile_records, "write_fits", _tprev, _t0)

    if write_diagnostics:
        payload = build_baseline_diagnostic_payload(
            linefree_mask=lf,
            ripple_freqs=freqs,
            rms_map=rms_map,
            flag_map=flag_map,
            extra={
                "input_fits": input_fits,
                "output_fits": output_fits,
                "cube_ext": cube_ext,
                "linefree_mode": linefree_mode,
                "ripple_mode": ripple_mode,
                "linefree_velocity_windows_kms": list(linefree_velocity_windows_kms) if linefree_velocity_windows_kms is not None else None,
                "safe_velocity_windows_kms": list(safe_velocity_windows_kms) if safe_velocity_windows_kms is not None else None,
                "exclude_v_windows": list(exclude_v_windows) if exclude_v_windows is not None else None,
                "ripple_apply_stage": str(ripple_apply_stage),
                "baseline_cfg": baseline_cfg.__dict__,
                "linefree_cfg": _linefree_cfg_to_meta_dict(linefree_cfg),
                "ripple_cfg": ripple_cfg.__dict__,
                "reproducible_mode": bool(getattr(baseline_cfg, "reproducible_mode", False)),
            },
        )
        prefix = diagnostics_prefix if diagnostics_prefix else output_fits
        json_path, txt_path = write_baseline_diagnostic_files(prefix, payload)
        logging.info("Wrote baseline diagnostics: %s | %s", json_path, txt_path)

    if write_profile:
        prefix = profile_prefix if profile_prefix else output_fits
        json_path, txt_path = _write_profile_sidecar(prefix, _profile_records)
        logging.info("Wrote baseline profile: %s | %s", json_path, txt_path)

    logging.info("Wrote baselined cube FITS: %s", output_fits)



_BASELINE_REPLACE_IMAGE_EXT = {"LINEFREE", "LINEFREE_USED", "LINEFREE3D", "LINEFREE3D_USED", "SIGNAL_MASK_USED", "SIGNAL_MASK3D_USED", "BASE_RMS", "BASE_FLG"}
_BASELINE_REPLACE_TABLE_EXT = {"RIPFREQ", "RIPFREQ_USED"}
_BASELINE_DROP_IMAGE_EXT = {
    "MOMENT0", "MOMENT1", "MOMENT2", "BASESUP3D", "LINECAND3D", "MOM0_BASESUP", "MOM0_LINECAND",
    "LINEFREE3D", "LINEFREE3D_USED", "SIGNAL_MASK3D_USED",
    "DATA_SIGNAL_ONLY", "DATA_LINEFREE_ONLY",
    "MOSAIC_RMS_OBS", "MOSAIC_RMS", "MOSAIC_WEIGHT",
}
_BASELINE_DROP_TABLE_EXT = {"MOSAIC_INFO"}


def _velocity_axis_from_linear_header(header: fits.Header, nchan: int) -> np.ndarray:
    ctype3 = str(header.get("CTYPE3", "") or "").upper()
    if not ctype3:
        raise ValueError("Header is missing CTYPE3; cannot build spectral axis for exclude_v_windows.")
    crpix = float(header.get("CRPIX3", 1.0))
    crval = float(header.get("CRVAL3", 0.0))
    cdelt = float(header.get("CDELT3", 1.0))
    idx = np.arange(int(nchan), dtype=float) + 1.0
    return crval + (idx - crpix) * cdelt


def _velocity_axis_from_fits_context(
    input_fits: str,
    header: fits.Header,
    *,
    cube_ext: Optional[Union[int, str]],
    hdu_index: int,
    nchan: int,
) -> np.ndarray:
    last_err: Optional[Exception] = None
    if SpectralCube is not None:
        try:
            sc = SpectralCube.read(input_fits, hdu=(cube_ext if cube_ext is not None else hdu_index))
            try:
                v_axis = np.asarray(sc.with_spectral_unit(u.km / u.s, velocity_convention="radio").spectral_axis.value, dtype=float)
            except Exception:
                v_axis = np.asarray(sc.spectral_axis.to(u.km / u.s).value, dtype=float)
            if v_axis.shape == (nchan,):
                return v_axis
        except Exception as exc:
            last_err = exc
    try:
        v_axis = _velocity_axis_from_linear_header(header, nchan)
        if v_axis.shape == (nchan,):
            return v_axis
    except Exception as exc:
        last_err = exc
    if last_err is not None:
        raise ValueError(
            "Could not build spectral axis in km/s for linefree/exclude windows from FITS input."
        ) from last_err
    raise ValueError("Could not build spectral axis in km/s for linefree/exclude windows from FITS input.")


def _read_linefree_prior_from_bundle(bundle: OTFBundle) -> Optional[np.ndarray]:
    for name in ("LINEFREE3D_USED", "LINEFREE3D", "LINEFREE3D_PRIOR", "LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR"):
        arr = bundle.image_ext.get(name)
        if arr is None:
            continue
        arrn = np.squeeze(np.asarray(arr))
        if arrn.ndim in (1, 3):
            return arrn.astype(bool)
    return None


def _read_ripple_prior_from_bundle(bundle: OTFBundle) -> Optional[List[float]]:
    for name in ("RIPFREQ", "RIPFREQ_USED", "RIPFREQ_PRIOR"):
        tbl = bundle.table_ext.get(name)
        if tbl is None:
            continue
        try:
            colnames = set(tbl.colnames)
        except Exception:
            continue
        if "FREQ_CYC_PER_CH" in colnames:
            arr = np.asarray(tbl["FREQ_CYC_PER_CH"], dtype=float)
            arr = arr[np.isfinite(arr)]
            if arr.size > 0:
                return [float(v) for v in arr]
    return None


def _drop_bundle_exts_for_baseline(bundle: OTFBundle) -> None:
    for key in list(_BASELINE_DROP_IMAGE_EXT | _BASELINE_REPLACE_IMAGE_EXT):
        bundle.image_ext.pop(key, None)
    for key in list(_BASELINE_DROP_TABLE_EXT | _BASELINE_REPLACE_TABLE_EXT):
        bundle.table_ext.pop(key, None)


def _attach_baseline_qc_exts_bundle(
    bundle: OTFBundle,
    *,
    linefree_mask_used: np.ndarray,
    signal_mask_used: np.ndarray,
    ripple_freqs_used: Sequence[float] | None,
    rms_map: np.ndarray | None,
    flag_map: np.ndarray | None,
    add_qc_ext: bool,
) -> None:
    lf = np.asarray(linefree_mask_used, dtype=bool)
    sig = np.asarray(signal_mask_used, dtype=bool)
    if sig.ndim == 1:
        bundle.image_ext["SIGNAL_MASK_USED"] = sig.copy()
    elif sig.ndim == 3:
        bundle.image_ext["SIGNAL_MASK3D_USED"] = sig.copy()
    else:
        raise ValueError(f"signal_mask_used must be 1D or 3D, got shape={sig.shape}")
    if not add_qc_ext:
        return
    if lf.ndim == 1:
        bundle.image_ext["LINEFREE"] = lf.copy()
        bundle.image_ext["LINEFREE_USED"] = lf.copy()
    elif lf.ndim == 3:
        bundle.image_ext["LINEFREE3D"] = lf.copy()
        bundle.image_ext["LINEFREE3D_USED"] = lf.copy()
    else:
        raise ValueError(f"linefree_mask_used must be 1D or 3D, got shape={lf.shape}")
    if ripple_freqs_used:
        freq_arr = np.asarray(list(ripple_freqs_used), dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            period_arr = np.where(freq_arr != 0.0, 1.0 / freq_arr, np.inf)
        tab = Table({"FREQ_CYC_PER_CH": freq_arr, "PERIOD_CH": period_arr})
        bundle.table_ext["RIPFREQ"] = tab.copy(copy_data=True)
        bundle.table_ext["RIPFREQ_USED"] = tab.copy(copy_data=True)
    if rms_map is not None:
        bundle.image_ext["BASE_RMS"] = np.asarray(rms_map, dtype=np.float32)
    if flag_map is not None:
        bundle.image_ext["BASE_FLG"] = np.asarray(flag_map, dtype=np.uint8)


def _attach_signal_products_after_baseline(
    bundle: OTFBundle,
    *,
    data_out: np.ndarray,
    linefree_mask_used: np.ndarray,
) -> None:
    velocity_axis_kms = None
    try:
        velocity_axis_kms = _velocity_axis_from_linear_header(bundle.header, int(np.asarray(data_out).shape[0]))
    except Exception as exc:
        logging.warning("Could not build velocity axis for MOMENT0; SIGNAL_MASK_USED will still be written. (%s)", exc)
    products = _derive_signal_products(
        np.asarray(data_out, dtype=np.float32),
        linefree_mask_used=np.asarray(linefree_mask_used, dtype=bool),
        velocity_axis_kms=velocity_axis_kms,
    )
    for name, arr in products.items():
        bundle.image_ext[name] = np.asarray(arr)


def _refresh_mosaic_after_baseline(
    bundle: OTFBundle,
    *,
    linefree_mask_used: np.ndarray,
    update_mosaic_products: bool,
    gain_min: float,
) -> None:
    for key in ("MOSAIC_RMS_OBS", "MOSAIC_RMS", "MOSAIC_WEIGHT"):
        bundle.image_ext.pop(key, None)
    bundle.table_ext.pop("MOSAIC_INFO", None)
    lf = np.asarray(linefree_mask_used, dtype=bool)
    if update_mosaic_products and lf.ndim == 1:
        attach_mosaic_products_from_mask(
            bundle,
            linefree_mask_1d=np.asarray(lf, dtype=bool),
            gain_min=float(gain_min),
            in_place=True,
        )
    elif update_mosaic_products and lf.ndim != 1:
        logging.info("Skipping mosaic product refresh for 3D voxel line-free mask.")



def _update_bundle_after_baseline(
    bundle_in: OTFBundle,
    data_out: np.ndarray,
    *,
    linefree_mask_used: np.ndarray,
    ripple_freqs_used: Sequence[float] | None,
    base_rms_map: np.ndarray | None,
    base_flag_map: np.ndarray | None,
    add_qc_ext: bool,
    update_mosaic_products: bool,
    gain_min: float,
) -> OTFBundle:
    out = bundle_in.copy(deep=True)
    out.data = np.asarray(data_out, dtype=np.float32)
    _drop_bundle_exts_for_baseline(out)
    signal_mask_used = _signal_mask_from_linefree_internal_gaps(np.asarray(linefree_mask_used, dtype=bool))
    _attach_baseline_qc_exts_bundle(
        out,
        linefree_mask_used=np.asarray(linefree_mask_used, dtype=bool),
        signal_mask_used=signal_mask_used,
        ripple_freqs_used=list(ripple_freqs_used or []),
        rms_map=base_rms_map,
        flag_map=base_flag_map,
        add_qc_ext=bool(add_qc_ext),
    )
    _attach_signal_products_after_baseline(
        out,
        data_out=np.asarray(data_out, dtype=np.float32),
        linefree_mask_used=np.asarray(linefree_mask_used, dtype=bool),
    )
    _refresh_mosaic_after_baseline(
        out,
        linefree_mask_used=np.asarray(linefree_mask_used, dtype=bool),
        update_mosaic_products=bool(update_mosaic_products),
        gain_min=float(gain_min),
    )
    lf_arr = np.asarray(linefree_mask_used, dtype=bool)
    sig_arr = np.asarray(signal_mask_used, dtype=bool)
    out.meta["baseline_linefree_mask_kind"] = ("global_1d" if lf_arr.ndim == 1 else "voxel_3d")
    out.meta["baseline_signal_mask_kind"] = ("global_1d" if sig_arr.ndim == 1 else "voxel_3d")
    out.meta["baseline_linefree_shape"] = tuple(int(v) for v in lf_arr.shape)
    out.meta["baseline_signal_shape"] = tuple(int(v) for v in sig_arr.shape)
    out.meta["baseline_linefree_ntrue"] = int(np.count_nonzero(lf_arr))
    out.meta["baseline_signal_ntrue"] = int(np.count_nonzero(sig_arr))
    out.meta["baseline_ripple_nfreq"] = int(len(list(ripple_freqs_used or [])))
    out.meta["baseline_mosaic_products_updated"] = bool(update_mosaic_products and lf_arr.ndim == 1)
    out.meta.setdefault("variance_note", "baseline fitting uncertainty is not propagated into variance")
    return out




def make_baseline_viewer_bundle(
    bundle: OTFBundle,
    *,
    mode: str = "signal",
    fill_value: float = 0.0,
) -> OTFBundle:
    """Return a copy of the input bundle with only selected baseline-viewer voxels kept in data.

    Parameters
    ----------
    bundle : OTFBundle
        Baseline-subtracted bundle containing line-free / signal masks.
    mode : {"signal", "linefree", "signal3d", "linefree3d"}
        Which spectral mask product to use for the viewer cube.
    fill_value : float
        Value written outside the selected mask, typically 0.0 for 3D viewers.
    """
    validate_otf_bundle(bundle, require_variance=False)
    mode = _validate_baseline_viewer_mode(mode)
    out = bundle.copy(deep=True)
    data = np.asarray(bundle.data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(f"bundle.data must be 3D (nchan, ny, nx), got shape={data.shape}")
    key_map = {
        "signal": "SIGNAL_MASK_USED",
        "linefree": "LINEFREE_USED",
        "signal3d": "SIGNAL_MASK3D_USED",
        "linefree3d": "LINEFREE3D_USED",
    }
    key = key_map[mode]
    if key not in bundle.image_ext:
        raise KeyError(f"Required image_ext {key!r} not found in bundle.")
    mask = np.asarray(bundle.image_ext[key], dtype=bool)
    if mask.ndim == 1 and mask.shape[0] != data.shape[0]:
        raise ValueError(f"{key} length mismatch: {mask.shape[0]} vs nchan={data.shape[0]}")
    if mask.ndim == 3 and mask.shape != data.shape:
        raise ValueError(f"{key} shape mismatch: {mask.shape} vs data shape={data.shape}")
    if mask.ndim not in (1, 3):
        raise ValueError(f"{key} must be 1D or 3D, got shape={mask.shape}")
    out.data = np.full_like(data, fill_value, dtype=np.float32)
    if mask.ndim == 1:
        if np.any(mask):
            out.data[mask, :, :] = data[mask, :, :]
    else:
        out.data[mask] = data[mask]
    out.meta["baseline_viewer_mode"] = mode
    out.meta["baseline_viewer_fill_value"] = float(fill_value)
    out.meta["baseline_viewer_mask_key"] = key
    return out

def subtract_baseline_from_bundle(
    bundle: OTFBundle,
    *,
    linefree_cfg: Optional[LineFreeConfig] = None,
    linefree_mask: Optional[np.ndarray] = None,
    linefree_detection_data: Optional[np.ndarray] = None,
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    linefree_mode: Optional[str] = "infer",
    load_prior_from_input: bool = True,
    ripple_cfg: RippleConfig = RippleConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    ripple_mode: str = "auto",
    safe_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    ripple_apply_stage: str = "final",
    baseline_cfg: BaselineConfig = BaselineConfig(),
    reproducible_mode: Optional[bool] = None,
    add_qc_ext: bool = True,
    add_qc_hdus: Optional[bool] = None,
    update_mosaic_products: bool = True,
    gain_min: float = 0.5,
) -> OTFBundle:
    """
    Subtract a baseline from an in-memory :class:`OTFBundle`.

    This is the main high-level orchestration entry point for bundle-based
    baseline subtraction.  It combines four conceptual layers:

    1. **Base line-free source selection**
       ``linefree_mode``, ``linefree_mask``, ``linefree_cfg``, and prior
       ``LINEFREE`` determine the initial fit mask.
    2. **Post-selection mask modifiers**
       ``linefree_velocity_windows_kms`` and ``exclude_v_windows`` modify the
       base mask after it has been chosen.
    3. **Ripple orchestration**
       ``ripple_freqs``, ``ripple_mode``, ``ripple_cfg``,
       ``safe_velocity_windows_kms``, and ``ripple_apply_stage`` determine how
       ripple frequencies are obtained and whether ripple is used for the
       detection-only cube, the final fit, or both.
    4. **Final fit model**
       ``baseline_cfg`` controls the polynomial/ripple model and numerical
       fitting behaviour.

    Parameter ownership map
    -----------------------
    Main API parameters (belong to this function)
        ``bundle``
        ``linefree_mask``
        ``linefree_detection_data``
        ``linefree_velocity_windows_kms``
        ``exclude_v_windows``
        ``linefree_mode``
        ``load_prior_from_input``
        ``ripple_freqs``
        ``ripple_mode``
        ``safe_velocity_windows_kms``
        ``ripple_apply_stage``
        ``reproducible_mode``
        ``add_qc_ext`` / ``add_qc_hdus``
        ``update_mosaic_products``
        ``gain_min``

    Config parameters
        ``linefree_cfg``   -> :class:`LineFreeConfig`
        ``ripple_cfg``     -> :class:`RippleConfig`
        ``baseline_cfg``   -> :class:`BaselineConfig`

    Effective line-free rules
    -------------------------
    ``linefree_mode='manual'``
        ``linefree_velocity_windows_kms`` is the **exact** line-free definition.
        Auto detection is not run.
    ``linefree_mode='prior'``
        Use prior ``LINEFREE`` from the input.  No auto detection unless the
        caller explicitly requests another mode.
    ``linefree_mode='auto'``
        Run automatic line-free detection from ``linefree_cfg``.
    ``linefree_mode='or'``
        Use the union of prior ``LINEFREE`` and the auto-detected mask.
    ``linefree_mode='infer'``
        Resolve automatically in this order:
        1. explicit ``linefree_mask``
        2. prior ``LINEFREE`` from the input
        3. auto detection if ``linefree_cfg`` is provided
        4. otherwise raise an error

    ``linefree_velocity_windows_kms`` behaviour
    -------------------------------------------
    - In ``manual`` mode, these windows define the exact line-free mask.
    - In all other effective modes, these windows are OR-added after the base
      mask is chosen.

    Parameters
    ----------
    bundle : OTFBundle
        Input cube bundle.  ``bundle.data`` must be shaped ``(nchan, ny, nx)``.
    linefree_cfg : LineFreeConfig or None
        Auto line-free detection configuration.  Only used when the effective
        mode uses auto detection.
    linefree_mask : np.ndarray or None
        Exact manual line-free mask, 1D ``(nchan,)`` or 3D ``(nchan, ny, nx)``.
        When provided and ``linefree_mode`` is ``infer``, this becomes the base
        line-free source.
    linefree_detection_data : np.ndarray or None
        Optional detection-only cube used for auto line-free estimation.  Shape
        must match ``bundle.data``.
    linefree_velocity_windows_kms : sequence or None
        Velocity windows in km/s.  Exact manual definition in ``manual`` mode;
        OR-style include windows in all other effective modes.
    exclude_v_windows : sequence or None
        Velocity windows removed from the final line-free mask after all other
        selection steps.
    linefree_mode : {'infer', 'manual', 'prior', 'auto', 'or'} or None
        Base line-free source selection.
    load_prior_from_input : bool
        Whether to read prior ``LINEFREE`` and ``RIPFREQ`` products from the
        input bundle when available.
    ripple_cfg : RippleConfig
        Configuration for automatic ripple-frequency search.
    ripple_freqs : sequence of float or None
        Explicit ripple frequencies in cycles/channel.  Highest priority source.
    ripple_mode : {'auto', 'prior'}
        Source of ripple frequencies when ``ripple_freqs`` is not supplied.
    safe_velocity_windows_kms : sequence or None
        Velocity windows used as a safe region when ripple is estimated for the
        detection cube or for the final FFT-based frequency search.
    ripple_apply_stage : {'final', 'mask_only', 'both'}
        Which stage may use ripple assistance.
    baseline_cfg : BaselineConfig
        Final baseline-model and numerical-fitting configuration.
    reproducible_mode : bool or None
        Optional override for deterministic numerical settings.
    add_qc_ext, add_qc_hdus : bool
        Whether QC image/table extensions are written into the output bundle.
        ``add_qc_hdus`` is accepted as a compatibility alias.
    update_mosaic_products : bool
        Whether derived mosaic products are refreshed after subtraction.
    gain_min : float
        Minimum edge-gain threshold used when mosaic products are refreshed.

    Returns
    -------
    OTFBundle
        Deep-copied output bundle containing the baseline-subtracted data and,
        when requested, QC products such as ``LINEFREE_USED`` and ``RIPFREQ``.
    """
    if add_qc_hdus is not None:
        add_qc_ext = bool(add_qc_hdus)
    validate_otf_bundle(bundle, require_variance=False)
    data = np.asarray(bundle.data, dtype=np.float32)
    nchan, ny, nx = data.shape

    exclude_v_windows = _resolve_exclude_v_windows(exclude_v_windows)
    linefree_velocity_windows_kms = _resolve_linefree_velocity_windows(linefree_velocity_windows_kms)
    safe_velocity_windows_kms = _resolve_linefree_velocity_windows(safe_velocity_windows_kms)
    ripple_apply_stage = _normalize_ripple_apply_stage(ripple_apply_stage)
    linefree_mode = _validate_linefree_selection_inputs(
        linefree_mask=linefree_mask,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        linefree_mode=linefree_mode,
    )

    lf_prior = _read_linefree_prior_from_bundle(bundle) if load_prior_from_input else None
    freqs_prior = _read_ripple_prior_from_bundle(bundle) if load_prior_from_input else None
    if lf_prior is not None:
        try:
            lf_prior = _normalize_linefree_mask_shape(np.asarray(lf_prior, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="LINEFREE prior")
        except Exception:
            logging.warning(
                "Ignoring prior LINEFREE because shape=%s does not match (%d,) or (%d,%d,%d).",
                np.asarray(lf_prior).shape,
                nchan,
                nchan,
                ny,
                nx,
            )
            lf_prior = None

    v_axis: Optional[np.ndarray] = None
    need_velocity_axis = bool(
        linefree_velocity_windows_kms
        or exclude_v_windows
        or safe_velocity_windows_kms
        or getattr(ripple_cfg, "period_range_kms", None) is not None
        or _linefree_cfg_needs_velocity_axis(linefree_cfg)
    )
    if need_velocity_axis:
        v_axis = _velocity_axis_from_linear_header(bundle.header, nchan)
        if v_axis.shape != (nchan,):
            raise ValueError(f"velocity axis shape mismatch: {v_axis.shape} vs ({nchan},)")
    spectral_axis_hz: Optional[np.ndarray] = None
    if getattr(ripple_cfg, "period_range_hz", None) is not None:
        spectral_axis_hz = _spectral_axis_hz_from_linear_header(bundle.header, nchan)
        if spectral_axis_hz.shape != (nchan,):
            raise ValueError(f"spectral axis shape mismatch: {spectral_axis_hz.shape} vs ({nchan},)")

    cell_arcsec = _header_cell_arcsec(bundle.header) if _linefree_cfg_needs_cell_arcsec(linefree_cfg) else None

    data_for_linefree = np.asarray(data, dtype=np.float32)
    detection_ripple_freqs: List[float] = []
    if linefree_detection_data is not None:
        arr_det = np.asarray(linefree_detection_data, dtype=np.float32)
        if arr_det.shape != data.shape:
            raise ValueError(f"linefree_detection_data shape mismatch: {arr_det.shape} vs {data.shape}")
        data_for_linefree = arr_det
    elif safe_velocity_windows_kms and ripple_apply_stage in {"mask_only", "both"} and bool(baseline_cfg.ripple):
        if v_axis is None:
            raise ValueError("safe_velocity_windows_kms requires a velocity axis.")
        safe_mask_1d = _linefree_mask_from_axis(np.asarray(v_axis, dtype=float), safe_velocity_windows_kms)
        if np.any(safe_mask_1d):
            if ripple_freqs is not None:
                detection_ripple_freqs = [float(f) for f in ripple_freqs]
            elif ripple_mode == "prior" and freqs_prior:
                detection_ripple_freqs = [float(f) for f in freqs_prior]
            else:
                flat_det = data.reshape(nchan, ny * nx)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    spec_det = np.nanmedian(flat_det, axis=1)
                detection_ripple_freqs = estimate_ripple_frequencies_fft(
                    spec_det,
                    safe_mask_1d,
                    rcfg=ripple_cfg,
                    poly_order_pre=baseline_cfg.poly_order,
                    velocity_axis_kms=v_axis,
                    spectral_axis_hz=spectral_axis_hz,
                )
            data_for_linefree = _fit_detection_cube_with_safe_windows(
                data,
                safe_mask_1d=safe_mask_1d,
                baseline_cfg=baseline_cfg,
                ripple_freqs=detection_ripple_freqs,
            )

    effective_linefree_mode = _resolve_effective_linefree_mode(
        linefree_mode=linefree_mode,
        linefree_mask=linefree_mask,
        linefree_cfg=linefree_cfg,
        has_prior=(lf_prior is not None),
    )
    if effective_linefree_mode == "manual_mask":
        lf = _normalize_linefree_mask_shape(np.asarray(linefree_mask, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="linefree_mask")
    elif effective_linefree_mode == "manual":
        if v_axis is None:
            raise ValueError("linefree_mode='manual' requires a velocity axis.")
        lf = _linefree_mask_from_axis(np.asarray(v_axis, dtype=float), linefree_velocity_windows_kms)
    else:
        lf_auto = None
        if effective_linefree_mode in ("auto", "or"):
            if linefree_cfg is None:
                raise ValueError(f"linefree_mode={effective_linefree_mode!r} requires linefree_cfg.")
            lf_auto = estimate_linefree_mask_from_cube(
                data_for_linefree,
                cfg=linefree_cfg,
                velocity_axis_kms=v_axis,
                cell_arcsec=cell_arcsec,
            )
            lf_auto = _normalize_linefree_mask_shape(np.asarray(lf_auto, dtype=bool), nchan=nchan, ny=ny, nx=nx, name="auto linefree_mask")

        if effective_linefree_mode == "prior":
            if lf_prior is None:
                raise ValueError("linefree_mode='prior' requested, but no prior LINEFREE was found in the bundle.")
            lf = np.asarray(lf_prior, dtype=bool)
        elif effective_linefree_mode == "or":
            if lf_prior is None and lf_auto is None:
                raise ValueError("linefree_mode='or' requested, but neither prior LINEFREE nor auto linefree_cfg is available.")
            lf = _merge_linefree_masks_or(lf_prior, lf_auto, nchan=nchan, ny=ny, nx=nx)
        elif effective_linefree_mode == "auto":
            if lf_auto is None:
                raise ValueError("linefree_mode='auto' requested, but auto detection did not produce a mask.")
            lf = np.asarray(lf_auto, dtype=bool)
        else:
            raise ValueError(f"Unknown effective linefree_mode: {effective_linefree_mode}")

    lf = _apply_linefree_include_exclude(
        np.asarray(lf, dtype=bool),
        velocity_axis_kms=v_axis,
        linefree_velocity_windows_kms=(None if effective_linefree_mode == "manual" else linefree_velocity_windows_kms),
        exclude_v_windows=exclude_v_windows,
    )

    freqs: List[float] = []
    if reproducible_mode is not None:
        baseline_cfg = replace(
            baseline_cfg,
            reproducible_mode=bool(reproducible_mode),
            compute_dtype=("float64" if bool(reproducible_mode) else str(getattr(baseline_cfg, "compute_dtype", "float32"))),
            normalize_x=(bool(reproducible_mode) or bool(getattr(baseline_cfg, "normalize_x", False))),
        )

    if baseline_cfg.ripple and ripple_apply_stage in {"final", "both"}:
        lf_for_ripple = _collapse_linefree_mask_for_ripple(np.asarray(lf, dtype=bool))
        if ripple_freqs is not None:
            freqs = [float(f) for f in ripple_freqs]
        elif ripple_mode == "prior" and freqs_prior:
            freqs = [float(f) for f in freqs_prior]
        else:
            ripple_mask_1d = lf_for_ripple
            if safe_velocity_windows_kms and v_axis is not None:
                ripple_mask_1d = _linefree_mask_from_axis(np.asarray(v_axis, dtype=float), safe_velocity_windows_kms)
            flat = data.reshape(nchan, ny * nx)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                spec = np.nanmedian(flat, axis=1)
            freqs = estimate_ripple_frequencies_fft(
                spec,
                ripple_mask_1d,
                rcfg=ripple_cfg,
                poly_order_pre=baseline_cfg.poly_order,
                velocity_axis_kms=v_axis,
                spectral_axis_hz=spectral_axis_hz,
            )

    out_cube, rms_map, flag_map = subtract_baseline_cube(
        data,
        linefree_mask=np.asarray(lf, dtype=bool),
        bcfg=baseline_cfg,
        ripple_freqs=freqs,
        return_qc=add_qc_ext,
    )
    out = _update_bundle_after_baseline(
        bundle,
        out_cube,
        linefree_mask_used=np.asarray(lf, dtype=bool),
        ripple_freqs_used=freqs,
        base_rms_map=rms_map,
        base_flag_map=flag_map,
        add_qc_ext=bool(add_qc_ext),
        update_mosaic_products=bool(update_mosaic_products),
        gain_min=float(gain_min),
    )
    out.meta["baseline_linefree_auto_method"] = None if linefree_cfg is None else str(getattr(linefree_cfg, "auto_method", "agg_1d"))
    out.meta["baseline_ripple_apply_stage"] = str(ripple_apply_stage)
    out.meta["baseline_safe_velocity_windows_kms"] = None if safe_velocity_windows_kms is None else list(safe_velocity_windows_kms)
    return out


def subtract_baseline_from_fits(
    input_fits: str,
    output_fits: str,
    *,
    cube_ext: Optional[Union[int, str]] = None,
    linefree_cfg: Optional[LineFreeConfig] = None,
    linefree_mask: Optional[np.ndarray] = None,
    linefree_detection_data: Optional[np.ndarray] = None,
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    linefree_mode: Optional[str] = "infer",
    load_prior_from_input: bool = True,
    ripple_cfg: RippleConfig = RippleConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    ripple_mode: str = "auto",
    safe_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    ripple_apply_stage: str = "final",
    baseline_cfg: BaselineConfig = BaselineConfig(),
    reproducible_mode: Optional[bool] = None,
    add_qc_hdus: bool = True,
    overwrite: bool = True,
    write_diagnostics: bool = False,
    diagnostics_prefix: Optional[str] = None,
    write_profile: bool = False,
    profile_prefix: Optional[str] = None,
) -> None:
    """
    Subtract a baseline from a FITS cube or OTF-bundle FITS on disk.

    This function exposes the same conceptual API as
    :func:`subtract_baseline_from_bundle`, but starts from a filename and writes
    the result back to disk.  When possible it uses the bundle path; otherwise
    it falls back to the legacy FITS path.

    Parameter ownership
    -------------------
    The ownership rules are the same as for
    :func:`subtract_baseline_from_bundle`:

    - main-API line-free / ripple orchestration parameters live directly on
      this function
    - ``linefree_cfg`` belongs to :class:`LineFreeConfig`
    - ``ripple_cfg`` belongs to :class:`RippleConfig`
    - ``baseline_cfg`` belongs to :class:`BaselineConfig`

    Practical interpretation
    ------------------------
    - ``linefree_mode='manual'`` means ``linefree_velocity_windows_kms`` is the
      exact line-free definition.
    - ``linefree_mode='infer'`` avoids silently enabling auto detection unless
      explicit auto inputs are provided.
    - ``linefree_velocity_windows_kms`` is OR-added only when the effective mode
      is not ``manual``.

    Parameters
    ----------
    input_fits, output_fits : str
        Input and output FITS paths.
    cube_ext : int, str, or None
        Optional legacy cube extension selector.  Supplying this forces the
        legacy FITS path.
    linefree_cfg, linefree_mask, linefree_detection_data,
    linefree_velocity_windows_kms, exclude_v_windows, linefree_mode,
    load_prior_from_input, ripple_cfg, ripple_freqs, ripple_mode,
    safe_velocity_windows_kms, ripple_apply_stage, baseline_cfg,
    reproducible_mode :
        Same meaning as in :func:`subtract_baseline_from_bundle`.
    add_qc_hdus : bool
        Whether QC HDUs/extensions are written.
    overwrite : bool
        Whether ``output_fits`` may be overwritten.
    write_diagnostics : bool
        Whether sidecar diagnostic files are written.
    diagnostics_prefix : str or None
        Prefix for diagnostic sidecars.
    write_profile : bool
        Whether the legacy path writes timing/profile sidecars.
    profile_prefix : str or None
        Prefix for profile sidecars.
    """
    use_legacy = (cube_ext is not None) or bool(write_profile)
    exclude_v_windows = _resolve_exclude_v_windows(exclude_v_windows)
    linefree_velocity_windows_kms = _resolve_linefree_velocity_windows(linefree_velocity_windows_kms)
    linefree_mode = _validate_linefree_selection_inputs(
        linefree_mask=linefree_mask,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        linefree_mode=linefree_mode,
    )
    if not use_legacy:
        try:
            bundle_in = read_otf_bundle(input_fits)
        except Exception as exc:
            logging.warning("Input is not readable as OTFBundle FITS (%s); falling back to legacy FITS path.", exc)
        else:
            bundle_out = subtract_baseline_from_bundle(
                bundle_in,
                linefree_cfg=linefree_cfg,
                linefree_mask=linefree_mask,
                linefree_detection_data=linefree_detection_data,
                linefree_velocity_windows_kms=linefree_velocity_windows_kms,
                exclude_v_windows=exclude_v_windows,
                linefree_mode=linefree_mode,
                load_prior_from_input=load_prior_from_input,
                ripple_cfg=ripple_cfg,
                ripple_freqs=ripple_freqs,
                ripple_mode=ripple_mode,
                safe_velocity_windows_kms=safe_velocity_windows_kms,
                ripple_apply_stage=ripple_apply_stage,
                baseline_cfg=baseline_cfg,
                reproducible_mode=reproducible_mode,
                add_qc_ext=add_qc_hdus,
                update_mosaic_products=True,
                gain_min=0.5,
            )
            write_otf_bundle(bundle_out, output_fits, overwrite=overwrite)
            if write_diagnostics:
                lf_src = bundle_out.image_ext.get("LINEFREE3D_USED", bundle_out.image_ext.get("LINEFREE3D", bundle_out.image_ext.get("LINEFREE_USED", bundle_out.image_ext.get("LINEFREE", []))))
                lf = np.asarray(lf_src, dtype=bool)
                freqs = []
                tbl = bundle_out.table_ext.get("RIPFREQ_USED") or bundle_out.table_ext.get("RIPFREQ")
                if tbl is not None and "FREQ_CYC_PER_CH" in getattr(tbl, 'colnames', []):
                    freqs = [float(v) for v in np.asarray(tbl["FREQ_CYC_PER_CH"], dtype=float) if np.isfinite(v)]
                rms_map = bundle_out.image_ext.get("BASE_RMS")
                flag_map = bundle_out.image_ext.get("BASE_FLG")
                payload = build_baseline_diagnostic_payload(
                    linefree_mask=lf,
                    ripple_freqs=freqs,
                    rms_map=None if rms_map is None else np.asarray(rms_map, dtype=float),
                    flag_map=None if flag_map is None else np.asarray(flag_map),
                    extra={
                        "input_fits": input_fits,
                        "output_fits": output_fits,
                        "cube_ext": cube_ext,
                        "linefree_mode": linefree_mode,
                        "ripple_mode": ripple_mode,
                        "linefree_velocity_windows_kms": list(linefree_velocity_windows_kms) if linefree_velocity_windows_kms is not None else None,
                        "exclude_v_windows": list(exclude_v_windows) if exclude_v_windows is not None else None,
                        "baseline_cfg": baseline_cfg.__dict__,
                        "linefree_cfg": _linefree_cfg_to_meta_dict(linefree_cfg),
                        "ripple_cfg": ripple_cfg.__dict__,
                        "reproducible_mode": bool(getattr(baseline_cfg, "reproducible_mode", False)),
                    },
                )
                prefix = diagnostics_prefix if diagnostics_prefix else output_fits
                write_baseline_diagnostic_files(prefix, payload)
            logging.info("Wrote baselined cube FITS via bundle wrapper: %s", output_fits)
            return None
    return _subtract_baseline_from_fits_legacy(
        input_fits,
        output_fits,
        cube_ext=cube_ext,
        linefree_cfg=linefree_cfg,
        linefree_mask=linefree_mask,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        exclude_v_windows=exclude_v_windows,
        linefree_mode=linefree_mode,
        load_prior_from_input=load_prior_from_input,
        ripple_cfg=ripple_cfg,
        ripple_freqs=ripple_freqs,
        ripple_mode=ripple_mode,
        baseline_cfg=baseline_cfg,
        reproducible_mode=reproducible_mode,
        add_qc_hdus=add_qc_hdus,
        overwrite=overwrite,
        write_diagnostics=write_diagnostics,
        diagnostics_prefix=diagnostics_prefix,
        write_profile=write_profile,
        profile_prefix=profile_prefix,
    )


# --- r14 override: Stage-B weak update with spatial coherence -----------------
def estimate_linefree_mask_from_cube_3d(
    cube_data: np.ndarray,
    cfg: LineFreeConfig = LineFreeConfig(),
    *,
    velocity_axis_kms: Optional[np.ndarray] = None,
    cell_arcsec: Optional[float] = None,
) -> np.ndarray:
    """Estimate a 3D voxel-wise line-free mask from a cube.

    r14 update:
    - Keep Stage-A local hard seed from r13b.
    - Stage-B provisional fit remains polynomial-only on the original cube.
    - Strong detection uses the unsmoothed residual metric.
    - Weak detection uses a *spatially smoothed positive metric* to add spatial
      coherence for broad / drifting weak emission without reviving global-1D seed.
    """
    cfg = _expand_linefree_profile(cfg)
    cfg = _normalize_linefree_cfg_public_units(
        cfg,
        velocity_axis_kms=velocity_axis_kms,
        cell_arcsec=cell_arcsec,
    )
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(f"cube_data must be 3D (nchan,ny,nx), got shape={data.shape}")
    nchan, ny, nx = data.shape
    npix = ny * nx
    flat = data.reshape(nchan, npix)
    valid = np.isfinite(flat)
    Y = np.where(valid, flat, 0.0).astype(np.float32, copy=False)

    prov_order = int(getattr(cfg, 'lf3d_prov_poly_order', 1))
    max_iter = int(getattr(cfg, 'lf3d_max_iter', 3))
    threshold = float(getattr(cfg, 'lf3d_threshold', 3.0))
    dilation = int(getattr(cfg, 'lf3d_dilation', 4))
    positive_only = bool(getattr(cfg, 'lf3d_positive_only', True))
    reg_eps = float(getattr(cfg, 'lf3d_reg', 1.0e-7))
    solver = str(getattr(cfg, 'lf3d_solver', 'qr') or 'qr').strip().lower()
    if solver not in {'qr', 'normal'}:
        raise ValueError(f"Unknown lf3d_solver={solver!r}; expected 'qr' or 'normal'")
    sigma_mode = str(getattr(cfg, 'lf3d_sigma_mode', 'negative') or 'negative').strip().lower()
    monotone_line = bool(getattr(cfg, 'lf3d_monotone_line', True))
    hyst_sigma = getattr(cfg, 'lf3d_hysteresis_sigma', None)
    if hyst_sigma is None:
        hyst_sigma = 1.5
    else:
        hyst_sigma = float(hyst_sigma)
    weak_spatial_sigma = float(getattr(cfg, 'lf3d_weak_spatial_sigma', 1.25) or 0.0)
    min_run_kms = getattr(cfg, 'lf3d_min_run_kms', None)
    if min_run_kms is not None:
        min_run_kms = float(min_run_kms)
    chunk_pix = max(1, int(getattr(cfg, 'lf3d_chunk_pix', 4096)))

    x = np.linspace(-1.0, 1.0, nchan, dtype=np.float32)
    X_prov = np.vander(x, prov_order + 1, increasing=True).astype(np.float32, copy=False)

    det_cube = _make_lf3d_detection_cube(data, cfg)
    seed3d = _make_lf3d_local_seed_3d(det_cube, cfg, positive_only=positive_only, pre_smoothed=True)
    if min_run_kms is not None and min_run_kms > 0.0:
        seed3d = _prune_short_spectral_runs(
            seed3d,
            velocity_axis_kms=velocity_axis_kms,
            min_run_kms=min_run_kms,
        )
    line_accum = np.asarray(seed3d, dtype=bool).reshape(nchan, npix).copy()
    W = valid & (~line_accum)

    baseline = np.zeros_like(Y, dtype=np.float32)
    min_pts = max(prov_order + 2, 3)
    qr_cap = int(getattr(cfg, 'lf3d_qr_batch_pix', getattr(cfg, 'voxel_qr_batch_pix', 1024)) or 1024)
    for _iter in range(max(0, max_iter)):
        for p0 in range(0, npix, chunk_pix):
            p1 = min(npix, p0 + chunk_pix)
            Yc = Y[:, p0:p1]
            Wc = W[:, p0:p1]
            try:
                coef = _solve_weighted_least_squares_batch(
                    X_prov,
                    Yc,
                    Wc.astype(np.float32, copy=False),
                    dtype=np.dtype(np.float32),
                    method=solver,
                    rcond=None,
                    reg_eps=reg_eps,
                    max_batch_pix=qr_cap,
                )
                baseline[:, p0:p1] = np.einsum('vk,kn->vn', X_prov, coef, optimize=True)
            except Exception:
                bc = np.zeros_like(Yc, dtype=np.float32)
                for j in range(p1 - p0):
                    mj = np.asarray(Wc[:, j], dtype=bool)
                    if np.count_nonzero(mj) < min_pts:
                        continue
                    try:
                        cj, *_ = np.linalg.lstsq(
                            X_prov[mj, :].astype(np.float64, copy=False),
                            Yc[mj, j].astype(np.float64, copy=False),
                            rcond=None,
                        )
                        bc[:, j] = np.asarray(X_prov @ cj, dtype=np.float32)
                    except Exception:
                        bc[:, j] = 0.0
                baseline[:, p0:p1] = bc

        residual = Y - baseline
        masked_res = np.where(W, residual, np.nan)
        med, sigma = _masked_residual_location_sigma(
            masked_res,
            positive_only=positive_only,
            sigma_mode=sigma_mode,
        )
        if positive_only:
            metric = residual - med[None, :]
            weak_source = np.where(valid, np.maximum(metric, 0.0), np.nan).reshape(nchan, ny, nx)
        else:
            metric = np.abs(residual - med[None, :])
            weak_source = np.where(valid, metric, np.nan).reshape(nchan, ny, nx)

        strong_mask = (metric > (threshold * sigma[None, :])) & valid

        if hyst_sigma is not None and hyst_sigma > 0.0 and hyst_sigma < threshold:
            if weak_spatial_sigma > 0.0:
                weak_cube = _gaussian_filter_spatial_nan_normalized(weak_source, weak_spatial_sigma)
                weak_metric = weak_cube.reshape(nchan, npix)
            else:
                weak_metric = metric
            weak_mask = (weak_metric > (hyst_sigma * sigma[None, :])) & valid
            line_mask = _spectral_hysteresis_no_wrap(strong_mask, weak_mask)
        else:
            line_mask = strong_mask

        if dilation > 0:
            line_mask = _spectral_dilate_no_wrap(line_mask, dilation)
        if min_run_kms is not None and min_run_kms > 0.0:
            line_mask = _prune_short_spectral_runs(
                line_mask,
                velocity_axis_kms=velocity_axis_kms,
                min_run_kms=min_run_kms,
            )

        if monotone_line:
            line_accum |= line_mask
            W_new = valid & (~line_accum)
        else:
            W_new = valid & (~line_mask)

        if np.array_equal(W_new, W):
            W = W_new
            break
        W = W_new

    return W.reshape(nchan, ny, nx)
