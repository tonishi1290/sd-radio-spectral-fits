# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map.baseline_subtraction

Single-dish spectral cube baseline subtraction with:
- automatic line-free channel detection (from spatially-aggregated spectrum)
- polynomial baseline
- multi-period ripple model as a sum of fixed-frequency sinusoids
- fast batch least-squares across pixels (when line-free mask is global or 3D vectorized)

Key ideas
---------
1) Fit only on "line-free" channels to avoid subtracting real spectral lines.
   - line-free mask can be user-specified OR auto-estimated via robust sigma-clipping.
   - Now supports both Global 1D masks and Vectorized 3D masks.

2) Standing-wave baseline ripple is quasi-sinusoidal; multiple periods can exist.
   - We estimate dominant ripple frequencies from the spatially aggregated spectrum.
   - Fixing frequencies makes per-pixel fitting linear (A*sin + B*cos).

3) For huge cubes: avoid 3D connected-component labeling.
   - Baseline fitting is independent per spectrum, chunked in pixels to save memory.
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
from dataclasses import dataclass
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


@dataclass(frozen=True)
class RippleConfig:
    nfreq: int = 2
    period_range_chan: Tuple[float, float] = (20.0, 400.0)
    min_separation: float = 0.002
    window: str = "hann"
    interpolate_masked: bool = False
    stable_sort: bool = True
    sampling_mode: str = "random"
    sampling_seed: int = 0


@dataclass(frozen=True)
class BaselineConfig:
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


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

def _robust_std(x: np.ndarray) -> float:
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
    if w < 3: return y
    if w % 2 == 0: w += 1
    return ndi.median_filter(y, size=w, mode="nearest")

def _dilate_1d(mask: np.ndarray, pad: int) -> np.ndarray:
    p = int(pad)
    if p <= 0: return mask
    struct = np.ones((2 * p + 1,), dtype=bool)
    return ndi.binary_dilation(mask, structure=struct)

def _as_float32(x: Any) -> np.ndarray:
    return np.asarray(x, dtype=np.float32)

def _capture_stdout(fn: Any) -> str:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try: fn()
        except Exception as exc: return f"<unavailable: {exc}>"
    return buf.getvalue().strip()

def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def _sha256_array(arr: Any) -> str:
    a = np.ascontiguousarray(np.asarray(arr))
    return _sha256_bytes(a.tobytes())

def _json_ready(obj: Any) -> Any:
    if isinstance(obj, dict): return {str(k): _json_ready(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)): return [_json_ready(v) for v in obj]
    if isinstance(obj, np.ndarray): return _json_ready(obj.tolist())
    if isinstance(obj, (np.floating,)): return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.integer,)): return int(obj)
    if isinstance(obj, (np.bool_,)): return bool(obj)
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
            "shape": list(lf.shape),
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
        f"shape        : {lf.get('shape')}",
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
    if base.lower().endswith('.fits'): base = base[:-5]
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
    if mode not in {"random", "stride"}: raise ValueError(f"Unknown sampling_mode={mode!r}")
    return mode

def _select_sample_indices(npix: int, max_pix: int, *, mode: str, seed: int) -> Optional[np.ndarray]:
    max_pix = int(max_pix)
    if max_pix <= 0: raise ValueError(f"max_pix must be positive, got {max_pix}")
    if npix <= max_pix: return None
    mode = _validate_sampling_mode(mode)
    if mode == "random":
        rng = np.random.default_rng(int(seed))
        return np.sort(rng.choice(npix, size=max_pix, replace=False).astype(np.int64, copy=False))
    step = float(npix) / float(max_pix)
    idx = np.floor(np.arange(max_pix, dtype=np.float64) * step).astype(np.int64)
    idx = np.clip(idx, 0, npix - 1)
    _, uniq = np.unique(idx, return_index=True)
    idx = idx[np.sort(uniq)]
    if idx.size < max_pix:
        need = max_pix - idx.size
        tail = np.setdiff1d(np.arange(npix - need, npix, dtype=np.int64), idx, assume_unique=False)
        idx = np.concatenate([idx, tail[:need]])
    return np.sort(idx[:max_pix])

def _interp_nans_1d(y: np.ndarray, *, dtype: np.dtype = np.float64) -> np.ndarray:
    arr = np.asarray(y, dtype=dtype)
    finite = np.isfinite(arr)
    if finite.all(): return arr
    if not finite.any(): return np.zeros_like(arr, dtype=dtype)
    x = np.arange(arr.size, dtype=np.float64)
    out = arr.copy()
    out[~finite] = np.interp(x[~finite], x[finite], out[finite])
    return out

def _fill_nan_with_median_along_spec(data: np.ndarray) -> np.ndarray:
    arr = np.asarray(data, dtype=np.float32)
    if not np.isnan(arr).any(): return arr
    med = np.nanmedian(arr, axis=0, keepdims=True)
    med = np.nan_to_num(med, nan=0.0).astype(np.float32, copy=False)
    return np.where(np.isfinite(arr), arr, med)

def _combine_linefree_masks_conservative(masks: Sequence[np.ndarray], nchan: int) -> np.ndarray:
    combined = np.ones(int(nchan), dtype=bool)
    used = 0
    for m in masks:
        arr = np.asarray(m, dtype=bool)
        combined &= arr
        used += 1
    if used == 0: return np.ones(int(nchan), dtype=bool)
    return combined

def _parse_manual_windows(windows: Sequence[Union[str, Tuple[float, float]]]) -> List[Tuple[float, float]]:
    parsed: List[Tuple[float, float]] = []
    for w in windows:
        if isinstance(w, str):
            s = w.strip()
            parts = s.split(":") if ":" in s else (s.split(",") if "," in s else s.split())
            lo, hi = float(parts[0]), float(parts[1])
        else:
            lo, hi = float(w[0]), float(w[1])
        if lo > hi: lo, hi = hi, lo
        parsed.append((lo, hi))
    return parsed

def _manual_signal_mask_from_axis(axis: np.ndarray, windows: Sequence[Union[str, Tuple[float, float]]]) -> np.ndarray:
    arr = np.asarray(axis, dtype=float)
    mask = np.zeros(arr.shape, dtype=bool)
    for lo, hi in _parse_manual_windows(windows):
        mask |= (arr >= float(lo)) & (arr <= float(hi))
    return mask

def _resolve_compute_dtype(bcfg: BaselineConfig) -> np.dtype:
    if bool(getattr(bcfg, "reproducible_mode", False)): return np.dtype(np.float64)
    return np.dtype(getattr(bcfg, "compute_dtype", "float32"))

def _resolve_rcond(rcond: Optional[float], *, n_rows: int, n_cols: int, dtype: np.dtype, reproducible_mode: bool) -> Optional[float]:
    if rcond is not None: return float(rcond)
    if not reproducible_mode: return None
    return float(np.finfo(dtype).eps * max(int(n_rows), int(n_cols)))

def _poly_coordinate(nchan: int, *, dtype: np.dtype, normalize: bool) -> np.ndarray:
    if not normalize: return np.arange(nchan, dtype=dtype)
    if nchan <= 1: return np.zeros(nchan, dtype=dtype)
    return np.linspace(-1.0, 1.0, nchan, dtype=dtype)


# -----------------------------------------------------------------------------
# 1) Line-free channel detection
# -----------------------------------------------------------------------------

def estimate_linefree_mask_1d(spectrum: np.ndarray, cfg: LineFreeConfig = LineFreeConfig()) -> np.ndarray:
    y = np.asarray(spectrum, dtype=np.float64)
    n = y.size
    if n < 8: return np.isfinite(y)

    y_work = _interp_nans_1d(y, dtype=np.float64) if bool(getattr(cfg, "interpolate_nans", True)) else np.nan_to_num(y, nan=0.0)
    base = _median_filter_1d(y_work, cfg.smooth_width)
    resid = y_work - base

    good = np.isfinite(y).copy()
    line = np.zeros(n, dtype=bool)

    for _ in range(int(cfg.iters)):
        cand = good & (~line)
        if cand.sum() < max(10, int(0.05 * n)): break
        s = _robust_std(resid[cand])
        if s <= 0: break
        new_line = np.abs(resid) > (float(cfg.sigma) * s)
        new_line &= good
        if np.array_equal(new_line, line): break
        line = new_line

    line = _dilate_1d(line, cfg.pad_chan)
    linefree = good & (~line)

    frac = float(linefree.sum()) / float(max(1, good.sum()))
    if frac < float(cfg.min_linefree_frac):
        raise ValueError(f"Too few line-free channels: frac={frac:.3f} < {cfg.min_linefree_frac:.3f}")
    return linefree


def estimate_linefree_mask_from_cube(
    cube_data: np.ndarray,
    cfg: LineFreeConfig = LineFreeConfig(),
    *,
    agg: str = "median",
    max_pix: int = 200000,
    seed: Optional[int] = None,
    sample_mode: Optional[str] = None,
) -> np.ndarray:
    data = np.asarray(cube_data, dtype=np.float32)
    nchan, ny, nx = data.shape
    npix = ny * nx
    flat = data.reshape(nchan, npix)

    mode = sample_mode if sample_mode is not None else getattr(cfg, "sampling_mode", "random")
    seed_eff = getattr(cfg, "sampling_seed", 0) if seed is None else int(seed)
    idx = _select_sample_indices(npix, int(max_pix), mode=mode, seed=seed_eff)
    flat_s = flat if idx is None else flat[:, idx]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        spec = np.nanmedian(flat_s, axis=1) if agg == "median" else np.nanmean(flat_s, axis=1)

    masks: List[np.ndarray] = [estimate_linefree_mask_1d(spec, cfg=cfg)]

    q = getattr(cfg, "conservative_quantile", None)
    if q is not None:
        q = float(q)
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
                logging.warning("Skipping conservative %s-quantile line-candidate spectrum: %s", qname, exc)

    return _combine_linefree_masks_conservative(masks, nchan)


# -----------------------------------------------------------------------------
# 2) Ripple frequency estimation
# -----------------------------------------------------------------------------

def estimate_ripple_frequencies_fft(
    spectrum: np.ndarray,
    linefree_mask: np.ndarray,
    rcfg: RippleConfig = RippleConfig(),
    *,
    poly_order_pre: int = 1,
) -> List[float]:
    y = np.asarray(spectrum, dtype=np.float64)
    m = np.asarray(linefree_mask, dtype=bool)
    n = y.size
    if n < 16 or m.sum() < max(8, poly_order_pre + 2): return []

    x_poly = _poly_coordinate(n, dtype=np.float64, normalize=True)
    if poly_order_pre >= 0:
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

    if rcfg.window.lower() == "hann": r2 *= np.hanning(n).astype(np.float64)

    spec = np.fft.rfft(r2)
    pwr = (spec.real ** 2 + spec.imag ** 2).astype(np.float64, copy=False)
    freqs = np.fft.rfftfreq(n, d=1.0)

    per_lo, per_hi = rcfg.period_range_chan
    fmin, fmax = 1.0 / float(per_hi), 1.0 / float(per_lo)
    sel = (freqs >= fmin) & (freqs <= fmax)
    if sel.size > 0: sel[0] = False
    if sel.sum() < 1: return []

    cand_idx = np.where(sel)[0]
    order = np.argsort(-pwr[cand_idx], kind="stable") if bool(getattr(rcfg, "stable_sort", True)) else np.argsort(pwr[cand_idx])[::-1]
    cand_idx = cand_idx[order]

    chosen: List[float] = []
    for i in cand_idx:
        f = float(freqs[i])
        if not np.isfinite(f) or f <= 0: continue
        if any(abs(f - fc) < float(rcfg.min_separation) for fc in chosen): continue
        chosen.append(f)
        if len(chosen) >= int(rcfg.nfreq): break

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
    x_poly = _poly_coordinate(int(nchan), dtype=dtype, normalize=bool(normalize_poly_x))
    x_phase = np.arange(int(nchan), dtype=dtype)
    cols = [x_poly ** k for k in range(int(poly_order) + 1)]
    if ripple_freqs:
        for f in ripple_freqs:
            w = dtype.type(2.0 * np.pi * float(f))
            cols.extend([np.sin(w * x_phase, dtype=dtype), np.cos(w * x_phase, dtype=dtype)])
    return np.vstack(cols).T.astype(dtype, copy=False)


def _robust_reweight(resid: np.ndarray, *, c: float = 1.345, dtype: np.dtype = np.float32) -> np.ndarray:
    r = np.asarray(resid, dtype=dtype)
    s = _robust_std(r)
    if s <= 0: return np.ones_like(r, dtype=dtype)
    t = np.abs(r) / (dtype.type(c * s))
    w = np.ones_like(r, dtype=dtype)
    m = t > 1
    w[m] = 1.0 / t[m]
    return w


def _robust_reweight_matrix(resid: np.ndarray, *, c: float = 1.345, dtype: np.dtype = np.float32) -> np.ndarray:
    r = np.asarray(resid, dtype=np.float64)
    med = np.nanmedian(r, axis=0, keepdims=True)
    mad = np.nanmedian(np.abs(r - med), axis=0)
    scale = 1.4826 * mad
    bad = (~np.isfinite(scale)) | (scale <= 0)
    if np.any(bad):
        scale = np.where(bad, np.nanstd(r, axis=0), scale)
    bad = (~np.isfinite(scale)) | (scale <= 0)
    scale = np.where(bad, 1.0, scale)
    t = np.abs(r) / (float(c) * scale[None, :])
    w = np.ones_like(r, dtype=np.float64)
    m = t > 1.0
    w[m] = 1.0 / t[m]
    if np.any(bad): w[:, bad] = 1.0
    return np.asarray(w, dtype=dtype)


def _coef_converged_mask(coef_old: np.ndarray, coef_new: np.ndarray, *, rtol: float, atol: float) -> np.ndarray:
    a, b = np.asarray(coef_old, dtype=np.float64), np.asarray(coef_new, dtype=np.float64)
    if a.ndim == 1:
        a, b = a[:, None], b[:, None]
    delta = np.max(np.abs(b - a), axis=0)
    scale = float(atol) + float(rtol) * np.maximum(np.max(np.abs(a), axis=0), np.max(np.abs(b), axis=0))
    return np.asarray(delta <= scale, dtype=bool)


def _solve_weighted_qr_batch(A: np.ndarray, Y: np.ndarray, W: np.ndarray, *, rcond: Optional[float], dtype: np.dtype) -> np.ndarray:
    A64, Y64, W64 = np.asarray(A, dtype=np.float64), np.asarray(Y, dtype=np.float64), np.asarray(W, dtype=np.float64)
    sqrtW = np.sqrt(np.clip(W64, 0.0, None))
    Aw = sqrtW.T[:, :, None] * A64[None, :, :]
    bw = sqrtW.T[:, :, None] * Y64.T[:, :, None]
    try:
        Q, R = np.linalg.qr(Aw, mode='reduced')
        Qtb = np.matmul(np.swapaxes(Q, -1, -2), bw)
        return np.asarray(np.linalg.solve(R, Qtb)[..., 0].T, dtype=dtype)
    except Exception:
        coef = np.empty((A64.shape[1], Y64.shape[1]), dtype=np.float64)
        for ib in range(Y64.shape[1]):
            try:
                Qi, Ri = np.linalg.qr(Aw[ib], mode='reduced')
                coef[:, ib] = np.linalg.solve(Ri, Qi.T @ bw[ib, :, 0])
            except Exception:
                coef[:, ib], *_ = np.linalg.lstsq(Aw[ib], bw[ib, :, 0], rcond=rcond)
        return np.asarray(coef, dtype=dtype)


def _fit_selected_robust_batched_1d(
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
    sel = np.asarray(selected_rel, dtype=np.int64)
    if sel.size == 0: return
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
            if active.size == 0: break
            coef_old = np.asarray(coef_sub[:, active], dtype=dtype, copy=True)
            Yact = np.asarray(Ysub[:, active], dtype=dtype)
            resid = Yact - (A_lf @ coef_old)
            W = _robust_reweight_matrix(resid, dtype=dtype)
            coef_new = _solve_weighted_qr_batch(A_lf, Yact, W, rcond=rcond, dtype=dtype)
            coef_sub[:, active] = coef_new
            if early:
                converged = _coef_converged_mask(coef_old, coef_new, rtol=robust_coef_rtol, atol=robust_coef_atol)
                if np.all(converged): break
                active = active[~converged]
        base = np.asarray(A_full @ coef_sub, dtype=dtype)
        out_chunk[:, sub_cols] = np.asarray(Yc_fill[:, sub_cols] - base, dtype=np.float32)


def _fit_selected_robust_batched_3d(
    out_chunk: np.ndarray,
    Yc_fill: np.ndarray,
    fit_cols: np.ndarray,
    selected_rel: np.ndarray,
    *,
    A_full: np.ndarray,
    linefree_mask_chunk: np.ndarray,
    coef_init: Optional[np.ndarray],
    dtype: np.dtype,
    robust_iters: int,
    robust_early_stop: bool,
    robust_coef_rtol: float,
    robust_coef_atol: float,
) -> None:
    """3D Vectorized IRLS using Normal Equations and Einsum for speed over variable masks."""
    sel = np.asarray(selected_rel, dtype=np.int64)
    if sel.size == 0: return

    sub_cols = np.asarray(fit_cols[sel], dtype=np.int64)
    W_base = np.asarray(linefree_mask_chunk[:, sub_cols], dtype=dtype)
    Ysub = np.asarray(Yc_fill[:, sub_cols], dtype=dtype)

    reg = np.eye(A_full.shape[1], dtype=dtype) * 1e-7

    if coef_init is not None and coef_init.shape[1] >= np.max(sel) + 1:
        coef_sub = np.asarray(coef_init[:, sel], dtype=dtype, copy=True)
    else:
        XTWX = np.einsum('vi,vn,vj->nij', A_full, W_base, A_full) + reg
        XTWY = np.einsum('vi,vn,vn->ni', A_full, W_base, Ysub)
        try:
            coef_sub = np.linalg.solve(XTWX, XTWY).T
        except Exception:
            coef_sub = np.zeros((A_full.shape[1], sel.size), dtype=dtype)

    active = np.arange(sel.size, dtype=np.int64)
    early = bool(robust_early_stop)

    for _ in range(max(0, int(robust_iters))):
        if active.size == 0: break
        coef_old = coef_sub[:, active]
        Yact = Ysub[:, active]
        Wact_base = W_base[:, active]

        resid = Yact - (A_full @ coef_old)
        resid_lf = np.where(Wact_base > 0, resid, np.nan)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            med = np.nanmedian(resid_lf, axis=0)
            mad = np.nanmedian(np.abs(resid_lf - med), axis=0)

        scale = 1.4826 * mad
        bad = (~np.isfinite(scale)) | (scale <= 0)
        if np.any(bad):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                scale2 = np.nanstd(resid_lf, axis=0)
            scale = np.where(bad, scale2, scale)

        bad = (~np.isfinite(scale)) | (scale <= 0)
        scale = np.where(bad, 1.0, scale)

        t = np.abs(resid) / (1.345 * scale[None, :])
        w_huber = np.where(t > 1.0, 1.0 / t, 1.0)
        W_new = (Wact_base * w_huber).astype(dtype)

        XTWX = np.einsum('vi,vn,vj->nij', A_full, W_new, A_full) + reg
        XTWY = np.einsum('vi,vn,vn->ni', A_full, W_new, Yact)

        try:
            coef_new = np.linalg.solve(XTWX, XTWY).T
        except Exception:
            coef_new = coef_old.copy()

        coef_sub[:, active] = coef_new

        if early:
            converged = _coef_converged_mask(coef_old, coef_new, rtol=robust_coef_rtol, atol=robust_coef_atol)
            if np.all(converged): break
            active = active[~converged]

    base = np.asarray(A_full @ coef_sub, dtype=dtype)
    out_chunk[:, sub_cols] = np.asarray(Yc_fill[:, sub_cols] - base, dtype=np.float32)


def _fit_selected_robust(
    out_chunk: np.ndarray,
    Yc_fill: np.ndarray,
    fit_cols: np.ndarray,
    selected_rel: np.ndarray,
    *,
    A_full: np.ndarray,
    A_lf: Optional[np.ndarray],
    linefree_mask: np.ndarray,
    mask_is_1d: bool,
    coef_init: Optional[np.ndarray] = None,
    rcond: Optional[float],
    dtype: np.dtype,
    robust_iters: int,
    batch_pixels: int,
    robust_early_stop: bool,
    robust_coef_rtol: float,
    robust_coef_atol: float,
) -> None:
    if mask_is_1d:
        _fit_selected_robust_batched_1d(
            out_chunk, Yc_fill, fit_cols, selected_rel, A_full=A_full, A_lf=A_lf,
            linefree_mask=linefree_mask, coef_init=coef_init, rcond=rcond, dtype=dtype,
            robust_iters=robust_iters, batch_pixels=batch_pixels, robust_early_stop=robust_early_stop,
            robust_coef_rtol=robust_coef_rtol, robust_coef_atol=robust_coef_atol
        )
    else:
        _fit_selected_robust_batched_3d(
            out_chunk, Yc_fill, fit_cols, selected_rel, A_full=A_full,
            linefree_mask_chunk=linefree_mask, coef_init=coef_init, dtype=dtype,
            robust_iters=robust_iters, robust_early_stop=robust_early_stop,
            robust_coef_rtol=robust_coef_rtol, robust_coef_atol=robust_coef_atol
        )


def _robust_loc_scale(x: np.ndarray) -> tuple[float, float]:
    arr = np.asarray(x, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0: return 0.0, 0.0
    med = float(np.nanmedian(arr))
    mad = float(np.nanmedian(np.abs(arr - med)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 0: scale = float(np.nanstd(arr))
    if not np.isfinite(scale) or scale <= 0: scale = 0.0
    return med, scale


def _selective_robust_candidates(resid_lf: np.ndarray, finite_lf: np.ndarray, *, sigma: float, max_frac: float, max_pixels: int) -> np.ndarray:
    r, f = np.asarray(resid_lf, dtype=np.float64), np.asarray(finite_lf, dtype=bool)
    if r.shape[1] == 0: return np.empty(0, dtype=np.int64)

    abs_r = np.abs(r)
    abs_r[~f] = np.nan
    cnt = np.sum(np.isfinite(abs_r), axis=0)
    good = cnt > 0
    if not np.any(good): return np.empty(0, dtype=np.int64)

    sq = np.where(np.isfinite(abs_r), abs_r * abs_r, 0.0)
    rms = np.full(abs_r.shape[1], np.nan, dtype=np.float64)
    rms[good] = np.sqrt(np.sum(sq[:, good], axis=0) / cnt[good])

    peak = np.max(np.where(np.isfinite(abs_r), abs_r, -np.inf), axis=0)
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
    if idx.size == 0: return idx.astype(np.int64, copy=False)

    limit = min(idx.size, int(np.ceil(float(max_frac) * int(np.sum(good)))) if float(max_frac) > 0 else idx.size)
    if int(max_pixels) > 0: limit = min(limit, int(max_pixels))
    if idx.size > limit:
        idx = idx[np.argsort(score[idx], kind='mergesort')[::-1][:limit]]
    return np.asarray(np.sort(idx), dtype=np.int64)


def _fit_single_spectrum(
    y_full: np.ndarray,
    A_full: np.ndarray,
    linefree_mask_1d: np.ndarray,
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
    A_lf = A_full[linefree_mask_1d, :]
    ylf = yj[linefree_mask_1d]
    cj, *_ = np.linalg.lstsq(A_lf, ylf, rcond=rcond)
    cj = np.asarray(cj, dtype=dtype)

    if robust:
        early = bool(robust_early_stop)
        for _ in range(max(0, int(robust_iters))):
            rj = ylf - (A_lf @ cj)
            w = _robust_reweight(rj, dtype=dtype)
            cj_new, *_ = np.linalg.lstsq(A_lf * w[:, None], ylf * w, rcond=rcond)
            cj_new = np.asarray(cj_new, dtype=dtype)
            if early and bool(_coef_converged_mask(cj, cj_new, rtol=robust_coef_rtol, atol=robust_coef_atol)[0]):
                cj = cj_new
                break
            cj = cj_new
    return np.asarray(yj - (A_full @ cj), dtype=np.float32)


def subtract_baseline_cube(
    cube_data: np.ndarray,
    *,
    linefree_mask: np.ndarray,
    bcfg: BaselineConfig = BaselineConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    return_qc: bool = True,
) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3: raise ValueError("cube_data must be 3D")
    nchan, ny, nx = data.shape
    
    m_in = np.asarray(linefree_mask, dtype=bool)
    if m_in.ndim == 1:
        if m_in.shape != (nchan,): raise ValueError(f"1D linefree_mask shape mismatch: {m_in.shape}")
        mask_is_1d = True
        m_flat = m_in
    elif m_in.ndim == 3:
        if m_in.shape != (nchan, ny, nx): raise ValueError(f"3D linefree_mask shape mismatch: {m_in.shape}")
        mask_is_1d = False
        m_flat = m_in.reshape(nchan, ny * nx)
    else:
        raise ValueError("linefree_mask must be 1D or 3D")

    dtype_compute = _resolve_compute_dtype(bcfg)
    normalize_x = bool(getattr(bcfg, "normalize_x", False) or getattr(bcfg, "reproducible_mode", False))
    robust_mode = str(getattr(bcfg, "robust_mode", "full") or "full").strip().lower()
    if robust_mode not in {"full", "batched_full", "selective"}: raise ValueError(f"Unknown robust_mode={robust_mode!r}")
    
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
    A_full = _design_matrix(nchan, poly_order=bcfg.poly_order, ripple_freqs=freqs, dtype=dtype_compute, normalize_poly_x=normalize_x)
    ncoef = int(A_full.shape[1])
    min_points = max(ncoef + 2, 8)

    A_lf = A_full[m_flat, :] if mask_is_1d else None
    
    if mask_is_1d and int(m_flat.sum()) < min_points:
        raise ValueError(f"Too few line-free points globally: {int(m_flat.sum())}")

    rcond_eff = _resolve_rcond(bcfg.rcond, n_rows=nchan, n_cols=ncoef, dtype=dtype_compute, reproducible_mode=bool(getattr(bcfg, "reproducible_mode", False)))

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
        
        m_chunk = m_flat if mask_is_1d else m_flat[:, p0:p1]
        finite_lf_counts = np.sum(finite_c[m_flat, :], axis=0) if mask_is_1d else np.sum(finite_c & m_chunk, axis=0)

        all_nan = np.sum(finite_c, axis=0) == 0
        insufficient = (~all_nan) & (finite_lf_counts < min_points)
        fit_mask = (~all_nan) & (~insufficient)

        flag_chunk[insufficient] = 1
        flag_chunk[all_nan] = 3

        fit_cols = np.where(fit_mask)[0]
        if fit_cols.size > 0:
            Yc_fit = np.asarray(Yc_fill[:, fit_cols], dtype=dtype_compute)
            try:
                if mask_is_1d:
                    Ylf = Yc_fit[m_flat, :]
                    coef, *_ = np.linalg.lstsq(A_lf, Ylf, rcond=rcond_eff)
                    coef = np.asarray(coef, dtype=dtype_compute)
                else:
                    W_fit = (m_chunk[:, fit_cols] & finite_c[:, fit_cols]).astype(dtype_compute)
                    reg = np.eye(ncoef, dtype=dtype_compute) * 1e-7
                    XTWX = np.einsum('vi,vn,vj->nij', A_full, W_fit, A_full) + reg
                    XTWY = np.einsum('vi,vn,vn->ni', A_full, W_fit, Yc_fit)
                    try:
                        coef = np.linalg.solve(XTWX, XTWY).T
                    except np.linalg.LinAlgError:
                        coef = np.zeros((ncoef, fit_cols.size), dtype=dtype_compute)

                base = A_full @ coef
                out_fit = np.asarray(Yc_fit - base, dtype=np.float32)
                out_chunk[:, fit_cols] = out_fit

                if bool(bcfg.robust):
                    if robust_mode == "full":
                        for j in fit_cols:
                            m_single = m_flat if mask_is_1d else m_chunk[:, j]
                            out_chunk[:, j] = _fit_single_spectrum(
                                Yc_fill[:, j], A_full, m_single, robust=True, robust_iters=robust_iters,
                                rcond=rcond_eff, dtype=dtype_compute, robust_early_stop=robust_early_stop,
                                robust_coef_rtol=robust_coef_rtol, robust_coef_atol=robust_coef_atol,
                            )
                    elif robust_mode == "batched_full":
                        _fit_selected_robust(
                            out_chunk, Yc_fill, fit_cols, np.arange(fit_cols.size, dtype=np.int64),
                            A_full=A_full, A_lf=A_lf, linefree_mask=m_chunk, mask_is_1d=mask_is_1d,
                            coef_init=coef, rcond=rcond_eff, dtype=dtype_compute,
                            robust_iters=robust_iters, batch_pixels=robust_batch_pixels,
                            robust_early_stop=robust_early_stop, robust_coef_rtol=robust_coef_rtol,
                            robust_coef_atol=robust_coef_atol,
                        )
                    else:
                        if mask_is_1d:
                            finite_lf = finite_c[m_flat, :][:, fit_cols]
                            resid_lf_init = out_fit[m_flat, :]
                        else:
                            m_chunk_fit = m_chunk[:, fit_cols]
                            finite_lf = finite_c[:, fit_cols] & m_chunk_fit
                            resid_lf_init = np.where(m_chunk_fit, out_fit, np.nan)
                            
                        selected_rel = _selective_robust_candidates(
                            resid_lf_init, finite_lf, sigma=robust_selective_sigma,
                            max_frac=robust_selective_frac, max_pixels=robust_selective_max_pixels,
                        )
                        if selected_rel.size > 0:
                            _fit_selected_robust(
                                out_chunk, Yc_fill, fit_cols, selected_rel,
                                A_full=A_full, A_lf=A_lf, linefree_mask=m_chunk, mask_is_1d=mask_is_1d,
                                coef_init=coef[:, selected_rel] if coef.ndim == 2 else None, rcond=rcond_eff, dtype=dtype_compute,
                                robust_iters=robust_iters, batch_pixels=robust_batch_pixels,
                                robust_early_stop=robust_early_stop, robust_coef_rtol=robust_coef_rtol,
                                robust_coef_atol=robust_coef_atol,
                            )
            except Exception:
                if bool(getattr(bcfg, "strict_failures", False)): raise
                if not bool(getattr(bcfg, "fallback_to_pixelwise", True)):
                    flag_chunk[fit_cols] = 2
                else:
                    ok_rel: list[int] = []
                    for jrel, j in enumerate(fit_cols):
                        try:
                            m_single = m_flat if mask_is_1d else m_chunk[:, j]
                            out_chunk[:, j] = _fit_single_spectrum(
                                Yc_fill[:, j], A_full, m_single,
                                robust=bool(bcfg.robust and robust_mode == "full"),
                                robust_iters=robust_iters, rcond=rcond_eff, dtype=dtype_compute,
                                robust_early_stop=robust_early_stop, robust_coef_rtol=robust_coef_rtol,
                                robust_coef_atol=robust_coef_atol,
                            )
                            ok_rel.append(int(jrel))
                        except Exception:
                            if bool(getattr(bcfg, "strict_failures", False)): raise
                            flag_chunk[j] = 2
                    if bool(bcfg.robust) and len(ok_rel) > 0:
                        fit_cols_ok = fit_cols[np.asarray(ok_rel, dtype=np.int64)]
                        if robust_mode == "batched_full":
                            _fit_selected_robust(
                                out_chunk, Yc_fill, fit_cols_ok, np.arange(fit_cols_ok.size, dtype=np.int64),
                                A_full=A_full, A_lf=A_lf, linefree_mask=m_chunk, mask_is_1d=mask_is_1d,
                                coef_init=None, rcond=rcond_eff, dtype=dtype_compute,
                                robust_iters=robust_iters, batch_pixels=robust_batch_pixels,
                                robust_early_stop=robust_early_stop, robust_coef_rtol=robust_coef_rtol, robust_coef_atol=robust_coef_atol,
                            )
                        elif robust_mode == "selective":
                            if mask_is_1d:
                                finite_lf = finite_c[m_flat, :][:, fit_cols_ok]
                                resid_lf_init = out_chunk[m_flat, :][:, fit_cols_ok]
                            else:
                                m_chunk_ok = m_chunk[:, fit_cols_ok]
                                finite_lf = finite_c[:, fit_cols_ok] & m_chunk_ok
                                resid_lf_init = np.where(m_chunk_ok, out_chunk[:, fit_cols_ok], np.nan)
                                
                            selected_rel = _selective_robust_candidates(
                                resid_lf_init, finite_lf, sigma=robust_selective_sigma,
                                max_frac=robust_selective_frac, max_pixels=robust_selective_max_pixels,
                            )
                            if selected_rel.size > 0:
                                _fit_selected_robust(
                                    out_chunk, Yc_fill, fit_cols_ok, selected_rel,
                                    A_full=A_full, A_lf=A_lf, linefree_mask=m_chunk, mask_is_1d=mask_is_1d,
                                    coef_init=None, rcond=rcond_eff, dtype=dtype_compute,
                                    robust_iters=robust_iters, batch_pixels=robust_batch_pixels,
                                    robust_early_stop=robust_early_stop, robust_coef_rtol=robust_coef_rtol, robust_coef_atol=robust_coef_atol,
                                )

        out_chunk[~finite_c] = np.nan
        out[:, p0:p1] = out_chunk

        if return_qc:
            ok = (flag_chunk == 0)
            if np.any(ok):
                lf_finite = finite_c[m_flat, :] & ok[None, :] if mask_is_1d else finite_c & m_chunk & ok[None, :]
                r = out_chunk[m_flat, :] ** 2 if mask_is_1d else out_chunk ** 2
                r[~lf_finite] = 0.0
                cnt = np.sum(lf_finite, axis=0)
                good = cnt > 0
                rms_local = np.full(k, np.nan, dtype=np.float32)
                if np.any(good):
                    rms_local[good] = np.sqrt(np.sum(r[:, good], axis=0) / cnt[good]).astype(np.float32, copy=False)
                resid_rms[p0:p1] = rms_local
            flag[p0:p1] = flag_chunk

    out_cube = out.reshape(nchan, ny, nx)
    if return_qc: return out_cube, resid_rms.reshape(ny, nx), flag.reshape(ny, nx)
    return out_cube, None, None


# -----------------------------------------------------------------------------
# 4) FITS wrapper (chunked write)
# -----------------------------------------------------------------------------

def _get_cube_hdu(hdul: fits.HDUList, cube_ext: Optional[Union[int, str]]) -> Tuple[int, fits.ImageHDU]:
    image_types = (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)
    excluded_names = {
        "MASK3D", "BASESUP3D", "LINECAND3D", "RMS", "BASE_RMS", "WEIGHT", "HIT", "MASK", "TSYS", "TINT", "TIME",
        "MOMENT0", "MOM0_BASESUP", "MOM0_LINECAND", "SIGNAL_MASK_USED", "DATA_SIGNAL_ONLY", "DATA_LINEFREE_ONLY",
    }
    excluded_btypes = {
        "SIGNALMASK", "BASELINESUPPORTMASK", "LINECANDIDATEMASK", "MOMENT0", "MOMENT0BASELINESUPPORT",
        "MOMENT0LINECANDIDATE", "SIGNALMASKUSED", "DATASIGNALONLY", "DATALINEFREEONLY", "WEIGHT", "HITCOUNT",
        "VALIDMASK", "SYSTEMTEMP", "INTEGRATIONTIME", "OBSERVATIONTIME", "BASELINERESIDRMS", "BASELINERMS",
    }

    def _is_3d_image(hdu: object) -> bool: return isinstance(hdu, image_types) and getattr(hdu, "data", None) is not None and np.ndim(hdu.data) == 3
    def _looks_like_analysis_hdu(hdu: object) -> bool: return (str(getattr(hdu, "name", "") or "").upper() in excluded_names) or (str(getattr(hdu, "header", {}).get("BTYPE", "")).upper() in excluded_btypes)
    def _has_spectral_axis(hdu: object) -> bool: return getattr(hdu, "header", None) is not None and _find_spectral_fits_axis(hdu.header) is not None

    if cube_ext is None:
        if _is_3d_image(hdul[0]) and not _looks_like_analysis_hdu(hdul[0]) and _has_spectral_axis(hdul[0]): return 0, hdul[0]
        for i, hdu in enumerate(hdul):
            if _is_3d_image(hdu) and not _looks_like_analysis_hdu(hdu) and _has_spectral_axis(hdu): return i, hdu
        raise ValueError("No suitable 3D spectral cube found in FITS.")

    hdu = hdul[cube_ext]
    if not _is_3d_image(hdu): raise ValueError(f"cube_ext={cube_ext!r} is not a 3D image HDU.")
    if _looks_like_analysis_hdu(hdu): raise ValueError(f"cube_ext={cube_ext!r} points to an analysis HDU.")
    if not _has_spectral_axis(hdu): raise ValueError(f"cube_ext={cube_ext!r} does not advertise a spectral axis.")
    for i, hh in enumerate(hdul):
        if hh is hdu: return i, hdu
    raise RuntimeError("Could not resolve cube_ext.")

def _find_spectral_fits_axis(header: fits.Header) -> Optional[int]:
    for ax in (1, 2, 3):
        if any(tok in str(header.get(f"CTYPE{ax}", "")).upper() for tok in ("FREQ", "VRAD", "VELO", "VOPT", "WAVE", "AWAV")): return ax
    return None

def _standardize_cube_for_processing(data: np.ndarray, header: fits.Header) -> Tuple[np.ndarray, str]:
    arr = np.asarray(data)
    if arr.ndim != 3: raise ValueError("Cube data must be 3D")
    fits_ax = _find_spectral_fits_axis(header)
    if fits_ax is None: raise ValueError("Could not identify spectral axis")
    np_ax = arr.ndim - int(fits_ax)
    if np_ax == 0: return arr, "v_y_x"
    if np_ax == 1: return np.transpose(arr, (1, 0, 2)), "y_v_x"
    if np_ax == 2: return np.transpose(arr, (2, 0, 1)), "y_x_v"
    raise ValueError("Could not map spectral axis")

def _restore_cube_axis_order(data_std: np.ndarray, axis_order_in: str) -> np.ndarray:
    arr = np.asarray(data_std)
    if axis_order_in == "v_y_x": return arr
    if axis_order_in == "y_v_x": return np.transpose(arr, (1, 0, 2))
    if axis_order_in == "y_x_v": return np.transpose(arr, (1, 2, 0))
    raise ValueError("Unknown axis order")

def _replace_or_append_hdu(hdul: fits.HDUList, hdu: Union[fits.ImageHDU, fits.CompImageHDU, fits.BinTableHDU]) -> None:
    name = str(getattr(hdu, "name", "")).upper()
    if not name:
        hdul.append(hdu)
        return
    while name in hdul: del hdul[name]
    hdul.append(hdu)

def _remove_named_hdus(hdul: fits.HDUList, names: Sequence[str], *, protect_hdu: object | None = None) -> None:
    names_up = {str(n).upper() for n in names}
    for i in reversed([j for j, h in enumerate(hdul) if j > 0 and h is not protect_hdu and str(getattr(h, 'name', '')).upper() in names_up]):
        del hdul[i]

def _strip_checksum_all_hdus(hdul: fits.HDUList) -> None:
    for hdu in hdul:
        if getattr(hdu, "header", None):
            for key in ("CHECKSUM", "DATASUM", "ZHECKSUM", "ZDATASUM"):
                while key in hdu.header: del hdu.header[key]

def _resolve_exclude_v_windows(w: Optional[Sequence[Union[str, Tuple[float, float]]]]) -> Optional[Sequence[Union[str, Tuple[float, float]]]]:
    return w if w and len(w) > 0 else None

def _resolve_linefree_velocity_windows(w: Optional[Sequence[Union[str, Tuple[float, float]]]]) -> Optional[Sequence[Union[str, Tuple[float, float]]]]:
    return w if w and len(w) > 0 else None

def _validate_linefree_selection_inputs(*, linefree_mask: Optional[np.ndarray], linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]], linefree_mode: str) -> None:
    if int(linefree_mask is not None) + int(linefree_velocity_windows_kms is not None) > 1: raise ValueError("Specify only one of linefree_mask and linefree_velocity_windows_kms.")
    if (linefree_mask is not None or linefree_velocity_windows_kms is not None) and linefree_mode in ("prior", "or"): raise ValueError(f"Explicit line-free specification cannot be combined with {linefree_mode!r}.")

def _linefree_mask_from_axis(axis: np.ndarray, windows: Sequence[Union[str, Tuple[float, float]]]) -> np.ndarray:
    return _manual_signal_mask_from_axis(axis, windows)

def _validate_baseline_viewer_mode(mode: str) -> str:
    m = str(mode).strip().lower()
    if m not in {"signal", "linefree"}: raise ValueError(f"Invalid mode {mode}")
    return m

def _true_runs_1d(mask: np.ndarray) -> list[tuple[int, int]]:
    arr = np.asarray(mask, dtype=bool).reshape(-1)
    if arr.size == 0: return []
    padded = np.concatenate(([False], arr, [False]))
    diff = np.diff(padded.astype(np.int8))
    return list(zip(np.flatnonzero(diff == 1), np.flatnonzero(diff == -1)))


def _signal_mask_from_linefree_internal_gaps(linefree_mask_used: np.ndarray) -> np.ndarray:
    """Vectorized internal gap filling for 1D and 3D line-free masks."""
    lf = np.asarray(linefree_mask_used, dtype=bool)
    if lf.ndim == 1:
        runs = _true_runs_1d(lf)
        sig = np.zeros(lf.shape, dtype=bool)
        if len(runs) < 2: return sig
        for (_, e0), (s1, _) in zip(runs[:-1], runs[1:]):
            if s1 > e0: sig[e0:s1] = True
        return sig & (~lf)
        
    nchan = lf.shape[0]
    lf_flat = lf.reshape(nchan, -1)
    idx = np.arange(nchan)[:, None]
    
    first_true = np.min(np.where(lf_flat, idx, nchan), axis=0)
    last_true = np.max(np.where(lf_flat, idx, -1), axis=0)
    
    in_range = (idx >= first_true) & (idx <= last_true)
    sig_flat = in_range & (~lf_flat)
    return sig_flat.reshape(lf.shape)


def _moment0_from_signal_mask(data_cube: np.ndarray, signal_mask_used: np.ndarray, velocity_axis_kms: np.ndarray) -> np.ndarray:
    data = np.asarray(data_cube, dtype=np.float32)
    sig = np.asarray(signal_mask_used, dtype=bool)
    v = np.asarray(velocity_axis_kms, dtype=float).reshape(-1)
    if sig.shape[0] != data.shape[0] or v.shape[0] != data.shape[0]: raise ValueError("Shape mismatch")
    dv = float(np.nanmedian(np.abs(np.diff(v)))) if data.shape[0] > 1 else 0.0
    if np.count_nonzero(sig) == 0 or not np.isfinite(dv) or dv <= 0.0: return np.zeros(data.shape[1:], dtype=np.float32)
    cube_sel = np.where(sig if sig.ndim == 3 else sig[:, None, None], data, 0.0)
    with np.errstate(invalid="ignore"): return np.asarray(np.nansum(cube_sel, axis=0, dtype=np.float64) * dv, dtype=np.float32)

def _viewer_cube_from_mask(data_cube: np.ndarray, mask_1d: np.ndarray) -> np.ndarray:
    data, m = np.asarray(data_cube, dtype=np.float32), np.asarray(mask_1d, dtype=bool).reshape(-1)
    out = np.zeros_like(data, dtype=np.float32)
    if np.any(m): out[m, :, :] = data[m, :, :]
    return out

def _derive_signal_products(data_cube: np.ndarray, *, linefree_mask_used: np.ndarray, velocity_axis_kms: np.ndarray | None) -> dict[str, np.ndarray]:
    lf = np.asarray(linefree_mask_used, dtype=bool)
    sig = _signal_mask_from_linefree_internal_gaps(lf)
    products = {"SIGNAL_MASK_USED": sig}
    if velocity_axis_kms is not None:
        try: products["MOMENT0"] = _moment0_from_signal_mask(data_cube, sig, velocity_axis_kms)
        except Exception as exc: logging.warning("Could not compute MOMENT0: %s", exc)
    return products

def _subtract_baseline_from_fits_legacy(
    input_fits: str, output_fits: str, *, cube_ext: Optional[Union[int, str]] = None,
    linefree_cfg: LineFreeConfig = LineFreeConfig(), linefree_mask: Optional[np.ndarray] = None,
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    linefree_mode: str = "auto", load_prior_from_input: bool = True,
    ripple_cfg: RippleConfig = RippleConfig(), ripple_freqs: Optional[Sequence[float]] = None,
    ripple_mode: str = "auto", baseline_cfg: BaselineConfig = BaselineConfig(),
    reproducible_mode: Optional[bool] = None, add_qc_hdus: bool = True, overwrite: bool = True,
    write_diagnostics: bool = False, diagnostics_prefix: Optional[str] = None,
    write_profile: bool = False, profile_prefix: Optional[str] = None,
) -> None:
    _profile_records: List[dict[str, Any]] = []
    _t0 = _tprev = time.perf_counter()
    if write_profile: _tprev = _profile_record(_profile_records, "start", _tprev, _t0, input_fits=str(input_fits))

    def _read_linefree_prior(hdul_local: fits.HDUList) -> Optional[np.ndarray]:
        for name in ("LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR"):
            if name in hdul_local:
                arr = np.squeeze(np.asarray(hdul_local[name].data))
                if arr.ndim in (1, 3): return arr.astype(bool)
        return None

    def _read_ripple_prior(hdul_local: fits.HDUList) -> Optional[List[float]]:
        for name in ("RIPFREQ", "RIPFREQ_USED", "RIPFREQ_PRIOR"):
            if name in hdul_local and hdul_local[name].data is not None and "FREQ_CYC_PER_CH" in getattr(hdul_local[name].data, "names", []):
                arr = np.asarray(hdul_local[name].data["FREQ_CYC_PER_CH"], dtype=float)
                if arr.size > 0: return [float(v) for v in arr]
        return None

    exclude_v_windows = _resolve_exclude_v_windows(exclude_v_windows)
    linefree_velocity_windows_kms = _resolve_linefree_velocity_windows(linefree_velocity_windows_kms)
    _validate_linefree_selection_inputs(linefree_mask=linefree_mask, linefree_velocity_windows_kms=linefree_velocity_windows_kms, linefree_mode=linefree_mode)

    with fits.open(input_fits, mode="readonly", memmap=True) as hdul_in:
        idx, hdu_cube = _get_cube_hdu(hdul_in, cube_ext)
        data, axis_order_in = _standardize_cube_for_processing(np.asarray(hdu_cube.data, dtype=np.float32), hdu_cube.header)
        nchan, ny, nx = data.shape
        if write_profile: _tprev = _profile_record(_profile_records, "read", _tprev, _t0)

        lf_prior = _read_linefree_prior(hdul_in) if load_prior_from_input else None
        freqs_prior = _read_ripple_prior(hdul_in) if load_prior_from_input else None
        if lf_prior is not None and np.asarray(lf_prior).shape not in ((nchan,), (nchan, ny, nx)): lf_prior = None

        if linefree_mask is not None:
            lf = np.asarray(linefree_mask, dtype=bool)
            if lf.ndim not in (1, 3): raise ValueError("linefree_mask must be 1D or 3D")
        elif linefree_velocity_windows_kms is not None:
            v_axis = _velocity_axis_from_fits_context(input_fits, hdu_cube.header, cube_ext=cube_ext, hdu_index=idx, nchan=nchan)
            lf = _linefree_mask_from_axis(v_axis, linefree_velocity_windows_kms)
        else:
            lf_auto = estimate_linefree_mask_from_cube(data, cfg=linefree_cfg) if linefree_mode in ("auto", "or") or (linefree_mode == "prior" and lf_prior is None) else None
            if linefree_mode == "prior": lf = np.asarray(lf_prior if lf_prior is not None else lf_auto, dtype=bool)
            elif linefree_mode == "or": lf = np.ones(nchan, dtype=bool) if lf_prior is None and lf_auto is None else np.asarray(lf_auto if lf_prior is None else (lf_prior if lf_auto is None else lf_prior | lf_auto), dtype=bool)
            elif linefree_mode == "auto": lf = np.asarray(lf_auto if lf_auto is not None else np.ones(nchan, dtype=bool), dtype=bool)
            else: raise ValueError("Unknown linefree_mode")

        if exclude_v_windows:
            v_axis = _velocity_axis_from_fits_context(input_fits, hdu_cube.header, cube_ext=cube_ext, hdu_index=idx, nchan=nchan)
            lf = np.asarray(lf, dtype=bool) & (~_manual_signal_mask_from_axis(v_axis, exclude_v_windows))

        logging.info("Line-free fraction: %.3f", lf.mean())

        freqs = []
        if reproducible_mode is not None:
            baseline_cfg = BaselineConfig(**{**baseline_cfg.__dict__, "reproducible_mode": bool(reproducible_mode), "compute_dtype": "float64" if reproducible_mode else baseline_cfg.compute_dtype})

        if baseline_cfg.ripple:
            if ripple_freqs is not None: freqs = [float(f) for f in ripple_freqs]
            elif ripple_mode == "prior" and freqs_prior: freqs = [float(f) for f in freqs_prior]
            else:
                spec = np.nanmedian(data.reshape(nchan, -1), axis=1)
                lf_ripple = lf if lf.ndim == 1 else (np.count_nonzero(lf, axis=(1, 2)) >= (ny * nx) / 2)
                freqs = estimate_ripple_frequencies_fft(spec, lf_ripple, rcfg=ripple_cfg, poly_order_pre=baseline_cfg.poly_order)

        out_cube, rms_map, flag_map = subtract_baseline_cube(data, linefree_mask=lf, bcfg=baseline_cfg, ripple_freqs=freqs, return_qc=add_qc_hdus)
        hdul_out = fits.HDUList([h.copy() for h in hdul_in])

    out_cube_write = _restore_cube_axis_order(out_cube, axis_order_in)
    hdul_out[idx].data = out_cube_write

    _remove_named_hdus(hdul_out, ["LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR", "RIPFREQ", "RIPFREQ_USED", "RIPFREQ_PRIOR", "BASE_RMS", "BASE_FLG", "BASEFLAG", "BASE_DIAG", "RMS", "MASK3D", "MOMENT0", "SIGNAL_MASK_USED", "DATA_SIGNAL_ONLY", "DATA_LINEFREE_ONLY", "BASESUP3D", "LINECAND3D", "MOM0_BASESUP", "MOM0_LINECAND"], protect_hdu=hdul_out[idx])
    _strip_checksum_all_hdus(hdul_out)

    products = _derive_signal_products(np.asarray(out_cube, dtype=np.float32), linefree_mask_used=np.asarray(lf, dtype=bool), velocity_axis_kms=(v_axis if "v_axis" in locals() else None))

    hdr_sig = fits.Header({"BTYPE": "SignalMaskUsed", "COMMENT": "1=used for MOMENT0"})
    _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=np.asarray(products["SIGNAL_MASK_USED"], dtype=np.uint8), header=hdr_sig, name="SIGNAL_MASK_USED"))

    if "MOMENT0" in products:
        hdr_m0 = fits.Header({"BTYPE": "Moment0", "BUNIT": f"{str(hdul_out[idx].header.get('BUNIT', '')).strip()} km/s".strip()})
        _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=_as_float32(products["MOMENT0"]), header=hdr_m0, name="MOMENT0"))

    if add_qc_hdus:
        hdr_lf = fits.Header({"BTYPE": "LineFreeMask", "COMMENT": "1=line-free channel"})
        _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=lf.astype(np.uint8), header=hdr_lf, name="LINEFREE"))
        _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=lf.astype(np.uint8), header=hdr_lf.copy(), name="LINEFREE_USED"))
        if baseline_cfg.ripple and freqs:
            f_arr = np.asarray(freqs, dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"): p_arr = np.where(f_arr != 0, 1.0 / f_arr, np.inf)
            tb = fits.BinTableHDU.from_columns([fits.Column(name="FREQ_CYC_PER_CH", format="D", array=f_arr), fits.Column(name="PERIOD_CH", format="D", array=p_arr)], name="RIPFREQ")
            _replace_or_append_hdu(hdul_out, tb)
            _replace_or_append_hdu(hdul_out, fits.BinTableHDU.from_columns(tb.columns, name="RIPFREQ_USED"))
        if rms_map is not None: _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=_as_float32(rms_map), header=fits.Header({"BTYPE": "BaselineResidRMS"}), name="BASE_RMS"))
        if flag_map is not None: _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=np.asarray(flag_map, dtype=np.uint8), header=fits.Header({"BTYPE": "BaselineFitFlag"}), name="BASE_FLG"))

    hdul_out.writeto(output_fits, overwrite=overwrite)
    if write_diagnostics: write_baseline_diagnostic_files(diagnostics_prefix or output_fits, build_baseline_diagnostic_payload(linefree_mask=lf, ripple_freqs=freqs, rms_map=rms_map, flag_map=flag_map))
    if write_profile: _write_profile_sidecar(profile_prefix or output_fits, _profile_records)


_BASELINE_REPLACE_IMAGE_EXT = {"LINEFREE", "LINEFREE_USED", "SIGNAL_MASK_USED", "BASE_RMS", "BASE_FLG"}
_BASELINE_REPLACE_TABLE_EXT = {"RIPFREQ", "RIPFREQ_USED"}
_BASELINE_DROP_IMAGE_EXT = {"MOMENT0", "MOMENT1", "MOMENT2", "BASESUP3D", "LINECAND3D", "MOM0_BASESUP", "MOM0_LINECAND", "DATA_SIGNAL_ONLY", "DATA_LINEFREE_ONLY", "MOSAIC_RMS_OBS", "MOSAIC_RMS", "MOSAIC_WEIGHT"}
_BASELINE_DROP_TABLE_EXT = {"MOSAIC_INFO"}

def _velocity_axis_from_linear_header(header: fits.Header, nchan: int) -> np.ndarray:
    if not str(header.get("CTYPE3", "")).upper(): raise ValueError("Missing CTYPE3")
    return float(header.get("CRVAL3", 0)) + (np.arange(nchan, dtype=float) + 1.0 - float(header.get("CRPIX3", 1))) * float(header.get("CDELT3", 1))

def _velocity_axis_from_fits_context(input_fits: str, header: fits.Header, *, cube_ext: Optional[Union[int, str]], hdu_index: int, nchan: int) -> np.ndarray:
    try: return _velocity_axis_from_linear_header(header, nchan)
    except Exception: raise ValueError("Could not build spectral axis in km/s")

def _read_linefree_prior_from_bundle(bundle: OTFBundle) -> Optional[np.ndarray]:
    for name in ("LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR"):
        if (arr := bundle.image_ext.get(name)) is not None:
            arr1 = np.squeeze(np.asarray(arr))
            if arr1.ndim in (1, 3): return np.asarray(arr1, dtype=bool)
    return None

def _read_ripple_prior_from_bundle(bundle: OTFBundle) -> Optional[List[float]]:
    for name in ("RIPFREQ", "RIPFREQ_USED", "RIPFREQ_PRIOR"):
        if (tbl := bundle.table_ext.get(name)) is not None and "FREQ_CYC_PER_CH" in set(getattr(tbl, 'colnames', [])):
            arr = np.asarray(tbl["FREQ_CYC_PER_CH"], dtype=float)
            if (arr := arr[np.isfinite(arr)]).size > 0: return [float(v) for v in arr]
    return None

def _drop_bundle_exts_for_baseline(bundle: OTFBundle) -> None:
    for key in _BASELINE_DROP_IMAGE_EXT | _BASELINE_REPLACE_IMAGE_EXT: bundle.image_ext.pop(key, None)
    for key in _BASELINE_DROP_TABLE_EXT | _BASELINE_REPLACE_TABLE_EXT: bundle.table_ext.pop(key, None)

def _attach_baseline_qc_exts_bundle(bundle: OTFBundle, *, linefree_mask_used: np.ndarray, signal_mask_used: np.ndarray, ripple_freqs_used: Sequence[float] | None, rms_map: np.ndarray | None, flag_map: np.ndarray | None, add_qc_ext: bool) -> None:
    bundle.image_ext["SIGNAL_MASK_USED"] = np.asarray(signal_mask_used, dtype=bool)
    if not add_qc_ext: return
    bundle.image_ext["LINEFREE"] = bundle.image_ext["LINEFREE_USED"] = np.asarray(linefree_mask_used, dtype=bool)
    if ripple_freqs_used:
        f = np.asarray(list(ripple_freqs_used), dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"): p = np.where(f != 0, 1.0 / f, np.inf)
        bundle.table_ext["RIPFREQ"] = bundle.table_ext["RIPFREQ_USED"] = Table({"FREQ_CYC_PER_CH": f, "PERIOD_CH": p})
    if rms_map is not None: bundle.image_ext["BASE_RMS"] = np.asarray(rms_map, dtype=np.float32)
    if flag_map is not None: bundle.image_ext["BASE_FLG"] = np.asarray(flag_map, dtype=np.uint8)

def _attach_signal_products_after_baseline(bundle: OTFBundle, *, data_out: np.ndarray, linefree_mask_used: np.ndarray) -> None:
    v_axis = None
    try: v_axis = _velocity_axis_from_linear_header(bundle.header, int(data_out.shape[0]))
    except Exception: pass
    for name, arr in _derive_signal_products(data_out, linefree_mask_used=linefree_mask_used, velocity_axis_kms=v_axis).items():
        bundle.image_ext[name] = arr

def _refresh_mosaic_after_baseline(bundle: OTFBundle, *, linefree_mask_used: np.ndarray, update_mosaic_products: bool, gain_min: float) -> None:
    for key in ("MOSAIC_RMS_OBS", "MOSAIC_RMS", "MOSAIC_WEIGHT"): bundle.image_ext.pop(key, None)
    bundle.table_ext.pop("MOSAIC_INFO", None)
    if update_mosaic_products:
        # Pass 1D flattened representation if 3D, or standard 1D mask
        lf_1d = np.asarray(linefree_mask_used, dtype=bool)
        if lf_1d.ndim == 3: lf_1d = np.count_nonzero(lf_1d, axis=(1, 2)) >= (lf_1d.shape[1] * lf_1d.shape[2]) / 2
        attach_mosaic_products_from_mask(bundle, linefree_mask_1d=lf_1d, gain_min=float(gain_min), in_place=True)

def _update_bundle_after_baseline(bundle_in: OTFBundle, data_out: np.ndarray, *, linefree_mask_used: np.ndarray, ripple_freqs_used: Sequence[float] | None, base_rms_map: np.ndarray | None, base_flag_map: np.ndarray | None, add_qc_ext: bool, update_mosaic_products: bool, gain_min: float) -> OTFBundle:
    out = bundle_in.copy(deep=True)
    out.data = np.asarray(data_out, dtype=np.float32)
    _drop_bundle_exts_for_baseline(out)
    signal = _signal_mask_from_linefree_internal_gaps(linefree_mask_used)
    _attach_baseline_qc_exts_bundle(out, linefree_mask_used=linefree_mask_used, signal_mask_used=signal, ripple_freqs_used=ripple_freqs_used, rms_map=base_rms_map, flag_map=base_flag_map, add_qc_ext=add_qc_ext)
    _attach_signal_products_after_baseline(out, data_out=data_out, linefree_mask_used=linefree_mask_used)
    _refresh_mosaic_after_baseline(out, linefree_mask_used=linefree_mask_used, update_mosaic_products=update_mosaic_products, gain_min=gain_min)
    out.meta["baseline_linefree_nchan"] = int(np.count_nonzero(linefree_mask_used))
    out.meta["baseline_signal_nchan"] = int(np.count_nonzero(signal))
    return out

def make_baseline_viewer_bundle(bundle: OTFBundle, *, mode: str = "signal", fill_value: float = 0.0) -> OTFBundle:
    validate_otf_bundle(bundle, require_variance=False)
    out = bundle.copy(deep=True)
    key = "SIGNAL_MASK_USED" if _validate_baseline_viewer_mode(mode) == "signal" else "LINEFREE_USED"
    mask = np.asarray(bundle.image_ext[key], dtype=bool)
    if mask.ndim == 1:
        out.data = np.full_like(bundle.data, fill_value, dtype=np.float32)
        if np.any(mask): out.data[mask, :, :] = bundle.data[mask, :, :]
    else:
        out.data = np.where(mask, bundle.data, fill_value).astype(np.float32)
    return out

def subtract_baseline_from_bundle(
    bundle: OTFBundle, *, linefree_cfg: LineFreeConfig = LineFreeConfig(), linefree_mask: Optional[np.ndarray] = None,
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None, exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    linefree_mode: str = "auto", load_prior_from_input: bool = True, ripple_cfg: RippleConfig = RippleConfig(),
    ripple_freqs: Optional[Sequence[float]] = None, ripple_mode: str = "auto", baseline_cfg: BaselineConfig = BaselineConfig(),
    reproducible_mode: Optional[bool] = None, add_qc_ext: bool = True, add_qc_hdus: Optional[bool] = None, update_mosaic_products: bool = True, gain_min: float = 0.5,
) -> OTFBundle:
    if add_qc_hdus is not None: add_qc_ext = bool(add_qc_hdus)
    validate_otf_bundle(bundle, require_variance=False)
    data = np.asarray(bundle.data, dtype=np.float32)
    nchan, ny, nx = data.shape

    _validate_linefree_selection_inputs(linefree_mask=linefree_mask, linefree_velocity_windows_kms=linefree_velocity_windows_kms, linefree_mode=linefree_mode)
    lf_prior = _read_linefree_prior_from_bundle(bundle) if load_prior_from_input else None
    
    if linefree_mask is not None:
        lf = np.asarray(linefree_mask, dtype=bool)
    elif linefree_velocity_windows_kms is not None:
        lf = _linefree_mask_from_axis(_velocity_axis_from_linear_header(bundle.header, nchan), linefree_velocity_windows_kms)
    else:
        lf_auto = estimate_linefree_mask_from_cube(data, cfg=linefree_cfg) if linefree_mode in ("auto", "or") or (linefree_mode == "prior" and lf_prior is None) else None
        if linefree_mode == "prior": lf = np.asarray(lf_prior if lf_prior is not None else lf_auto, dtype=bool)
        elif linefree_mode == "or": lf = np.ones(nchan, dtype=bool) if lf_prior is None and lf_auto is None else np.asarray(lf_auto if lf_prior is None else (lf_prior if lf_auto is None else lf_prior | lf_auto), dtype=bool)
        elif linefree_mode == "auto": lf = np.asarray(lf_auto if lf_auto is not None else np.ones(nchan, dtype=bool), dtype=bool)
        else: raise ValueError("Unknown linefree_mode")

    if exclude_v_windows: lf = np.asarray(lf, dtype=bool) & (~_manual_signal_mask_from_axis(_velocity_axis_from_linear_header(bundle.header, nchan), exclude_v_windows))

    freqs = []
    if baseline_cfg.ripple:
        if ripple_freqs is not None: freqs = [float(f) for f in ripple_freqs]
        elif ripple_mode == "prior" and (fp := _read_ripple_prior_from_bundle(bundle) if load_prior_from_input else None): freqs = [float(f) for f in fp]
        else:
            lf_rip = lf if lf.ndim == 1 else (np.count_nonzero(lf, axis=(1, 2)) >= (ny * nx) / 2)
            freqs = estimate_ripple_frequencies_fft(np.nanmedian(data.reshape(nchan, -1), axis=1), lf_rip, rcfg=ripple_cfg, poly_order_pre=baseline_cfg.poly_order)

    out_cube, rms_map, flag_map = subtract_baseline_cube(data, linefree_mask=lf, bcfg=baseline_cfg, ripple_freqs=freqs, return_qc=add_qc_ext)
    return _update_bundle_after_baseline(bundle, out_cube, linefree_mask_used=lf, ripple_freqs_used=freqs, base_rms_map=rms_map, base_flag_map=flag_map, add_qc_ext=add_qc_ext, update_mosaic_products=update_mosaic_products, gain_min=gain_min)


def subtract_baseline_from_fits(
    input_fits: str, output_fits: str, *, cube_ext: Optional[Union[int, str]] = None,
    linefree_cfg: LineFreeConfig = LineFreeConfig(), linefree_mask: Optional[np.ndarray] = None,
    linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None, exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    linefree_mode: str = "auto", load_prior_from_input: bool = True, ripple_cfg: RippleConfig = RippleConfig(),
    ripple_freqs: Optional[Sequence[float]] = None, ripple_mode: str = "auto", baseline_cfg: BaselineConfig = BaselineConfig(),
    reproducible_mode: Optional[bool] = None, add_qc_hdus: bool = True, overwrite: bool = True,
    write_diagnostics: bool = False, diagnostics_prefix: Optional[str] = None, write_profile: bool = False, profile_prefix: Optional[str] = None,
) -> None:
    if cube_ext is not None or write_profile:
        return _subtract_baseline_from_fits_legacy(input_fits, output_fits, cube_ext=cube_ext, linefree_cfg=linefree_cfg, linefree_mask=linefree_mask, linefree_velocity_windows_kms=linefree_velocity_windows_kms, exclude_v_windows=exclude_v_windows, linefree_mode=linefree_mode, load_prior_from_input=load_prior_from_input, ripple_cfg=ripple_cfg, ripple_freqs=ripple_freqs, ripple_mode=ripple_mode, baseline_cfg=baseline_cfg, reproducible_mode=reproducible_mode, add_qc_hdus=add_qc_hdus, overwrite=overwrite, write_diagnostics=write_diagnostics, diagnostics_prefix=diagnostics_prefix, write_profile=write_profile, profile_prefix=profile_prefix)
    try: bundle_in = read_otf_bundle(input_fits)
    except Exception: return _subtract_baseline_from_fits_legacy(input_fits, output_fits, cube_ext=cube_ext, linefree_cfg=linefree_cfg, linefree_mask=linefree_mask, linefree_velocity_windows_kms=linefree_velocity_windows_kms, exclude_v_windows=exclude_v_windows, linefree_mode=linefree_mode, load_prior_from_input=load_prior_from_input, ripple_cfg=ripple_cfg, ripple_freqs=ripple_freqs, ripple_mode=ripple_mode, baseline_cfg=baseline_cfg, reproducible_mode=reproducible_mode, add_qc_hdus=add_qc_hdus, overwrite=overwrite, write_diagnostics=write_diagnostics, diagnostics_prefix=diagnostics_prefix, write_profile=write_profile, profile_prefix=profile_prefix)
    bundle_out = subtract_baseline_from_bundle(bundle_in, linefree_cfg=linefree_cfg, linefree_mask=linefree_mask, linefree_velocity_windows_kms=linefree_velocity_windows_kms, exclude_v_windows=exclude_v_windows, linefree_mode=linefree_mode, load_prior_from_input=load_prior_from_input, ripple_cfg=ripple_cfg, ripple_freqs=ripple_freqs, ripple_mode=ripple_mode, baseline_cfg=baseline_cfg, reproducible_mode=reproducible_mode, add_qc_ext=add_qc_hdus)
    write_otf_bundle(bundle_out, output_fits, overwrite=overwrite)
    if write_diagnostics:
        lf = np.asarray(bundle_out.image_ext.get("LINEFREE_USED", bundle_out.image_ext.get("LINEFREE", [])), dtype=bool)
        freqs = [float(v) for v in np.asarray(tbl["FREQ_CYC_PER_CH"], dtype=float) if np.isfinite(v)] if (tbl := bundle_out.table_ext.get("RIPFREQ_USED") or bundle_out.table_ext.get("RIPFREQ")) is not None and "FREQ_CYC_PER_CH" in getattr(tbl, 'colnames', []) else []
        write_baseline_diagnostic_files(diagnostics_prefix or output_fits, build_baseline_diagnostic_payload(linefree_mask=lf, ripple_freqs=freqs, rms_map=bundle_out.image_ext.get("BASE_RMS"), flag_map=bundle_out.image_ext.get("BASE_FLG")))
