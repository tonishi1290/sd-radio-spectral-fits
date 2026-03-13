# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map.baseline_subtraction

Single-dish spectral cube baseline subtraction with:
- automatic line-free channel detection (from spatially-aggregated spectrum)
- polynomial baseline
- multi-period ripple model as a sum of fixed-frequency sinusoids
- fast batch least-squares across pixels (when line-free mask is global)

Key ideas
---------
1) Fit only on "line-free" channels to avoid subtracting real spectral lines.
   - line-free mask can be user-specified OR auto-estimated via robust sigma-clipping.
   - GILDAS/CLASS uses the same practical idea: baseline fits exclude user-defined windows. (see CLASS manual)

2) Standing-wave baseline ripple is quasi-sinusoidal; multiple periods can exist.
   - We estimate dominant ripple frequencies from the spatially aggregated spectrum using FFT on line-free channels,
     then fit each pixel with those frequencies fixed.
   - Fixing frequencies makes per-pixel fitting linear (A*sin + B*cos), which is fast and stable.

3) For huge cubes: avoid 3D connected-component labeling.
   - Baseline fitting is independent per spectrum, so we can chunk in pixels and write output incrementally.

Notes
-----
- Axis convention: data arrays are handled in (nchan, ny, nx).
- This module is intentionally conservative and "pipeline-safe":
  it provides explicit QC outputs and lets the user override auto decisions.

Dependencies
------------
- numpy, scipy, astropy
- (optional) spectral-cube for reading and spectral-axis unit conversions

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

try:
    import psutil  # type: ignore
except Exception:  # pragma: no cover
    psutil = None  # type: ignore

try:
    from spectral_cube import SpectralCube
except Exception:  # pragma: no cover
    SpectralCube = None  # type: ignore




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
    Configuration for automatic line-free channel detection.

    Parameters
    ----------
    smooth_width : int
        Median-filter width (channels) used to estimate slow baseline before sigma-clipping.
        Larger => less sensitive to broad lines, but may smooth ripple.
    sigma : float
        Sigma threshold for clipping |residual| > sigma * robust_std.
    iters : int
        Iterations of sigma-clipping.
    pad_chan : int
        Dilate detected line channels along spectral axis by this many channels.
    min_linefree_frac : float
        If the resulting line-free fraction is below this, raise (or warn) because baseline fit becomes ill-posed.
    sampling_mode : {'random', 'stride'}
        Spatial subsampling mode when a cube is too large to aggregate fully.
    sampling_seed : int
        RNG seed used when sampling_mode='random'.
    interpolate_nans : bool
        If True, fill NaN gaps in the aggregated 1D spectrum by linear interpolation before line-free detection.
    conservative_quantile : float or None
        If set to a quantile between 0.5 and 1.0 (for example 0.9), estimate additional
        line-candidate masks from upper-quantile and optionally lower-quantile aggregated spectra,
        then combine them conservatively with the primary aggregated spectrum. This helps catch
        spatially sparse emission / absorption lines that a global median can miss.
    conservative_both_tails : bool
        If True and conservative_quantile is enabled, also inspect the symmetric lower quantile
        (1 - conservative_quantile) to catch absorption-like sparse features.
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


@dataclass(frozen=True)
class RippleConfig:
    """
    Configuration for ripple (standing wave) modeling and auto frequency estimation.

    Units
    -----
    Frequencies are expressed in cycles per channel (0..0.5, Nyquist).
    Period in channels is 1/frequency.

    Parameters
    ----------
    nfreq : int
        Number of ripple frequencies to estimate (top peaks).
    period_range_chan : (float, float)
        Period range (channels) to consider. Example (20, 400).
    min_separation : float
        Minimum separation in frequency (cycles/channel) between selected peaks.
    window : str
        Window applied before FFT: 'hann' or 'none'.
    interpolate_masked : bool
        If True, fill non-linefree channels by linear interpolation before FFT instead of zero-filling.
    stable_sort : bool
        If True, use a stable sort when ranking FFT peaks.
    sampling_mode : {'random', 'stride'}
        Spatial subsampling mode used by higher-level orchestration when an aggregated spectrum is built.
    sampling_seed : int
        RNG seed used when sampling_mode='random'.
    """
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
    """
    Baseline model configuration.

    Parameters
    ----------
    poly_order : int
        Polynomial order (0=constant, 1=linear, ...).
    ripple : bool
        If True, include sinusoid terms with frequencies from RippleConfig (or user-provided).
    robust : bool
        If True, enable robust refitting. The exact behaviour is controlled by ``robust_mode``.
    robust_mode : {'full', 'batched_full', 'selective'}
        'full' keeps the original behaviour and robust-refits every fitted pixel.
        'batched_full' robust-refits every fitted pixel too, but solves weighted least-squares
        in small sub-batches via batched QR where possible. It is usually faster than
        'full', but may not be bit-identical.
        'selective' first performs a fast ordinary least-squares solve for the whole chunk and
        then robust-refits only the worst residual spectra.
    robust_iters : int
        Number of reweighting iterations for robust refits.
    robust_early_stop : bool
        If True, stop robust reweighting early when coefficients have converged.
    robust_coef_rtol : float
        Relative tolerance for coefficient-based early stopping.
    robust_coef_atol : float
        Absolute tolerance for coefficient-based early stopping.
    robust_selective_sigma : float
        Threshold (in robust-sigma units) for selecting spectra to robust-refit in
        ``robust_mode='selective'``. The threshold is applied to both line-free RMS and
        line-free max|residual| statistics.
    robust_selective_frac : float
        At most this fraction of spectra in a chunk are robust-refit in selective mode.
    robust_selective_max_pixels : int
        Hard cap on the number of spectra robust-refit per chunk in selective mode.
    rcond : float or None
        rcond for numpy.linalg.lstsq. If None and reproducible_mode=True, an explicit dtype-based threshold is used.
    chunk_pix : int
        Pixel chunk size when processing large cubes in memory.
    reproducible_mode : bool
        If True, use deterministic / version-stable numerical settings where practical.
    compute_dtype : {'float32', 'float64'}
        Internal dtype for least-squares and matrix construction.
    normalize_x : bool
        If True, normalize the polynomial coordinate to [-1, 1] before building x^k columns.
        Ripple sin/cos terms always use the native channel index so that frequencies remain cycles/channel.
    strict_failures : bool
        If True, re-raise fitting failures instead of silently flagging affected pixels.
    fallback_to_pixelwise : bool
        If True, fall back to pixel-wise fitting when a vectorized chunk solve fails.
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
            "nchan": int(lf.size),
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


def estimate_linefree_mask_from_cube(
    cube_data: np.ndarray,
    cfg: LineFreeConfig = LineFreeConfig(),
    *,
    agg: str = "median",
    max_pix: int = 200000,
    seed: Optional[int] = None,
    sample_mode: Optional[str] = None,
) -> np.ndarray:
    """
    Estimate a global line-free mask from a 3D cube by aggregating spectra spatially.

    Parameters
    ----------
    cube_data : np.ndarray
        (nchan, ny, nx) cube.
    agg : {'median','mean'}
        Spatial aggregation method for the primary spectrum.
    max_pix : int
        If ny*nx is huge, sample up to max_pix pixels for aggregation.
    seed : int or None
        RNG seed for sampling when sample_mode='random'. If None, cfg.sampling_seed is used.
    sample_mode : {'random','stride'} or None
        Sampling mode for large cubes. If None, cfg.sampling_mode is used.

    Returns
    -------
    linefree : np.ndarray (bool), shape (nchan,)

    Notes
    -----
    If cfg.conservative_quantile is enabled, additional masks are estimated from upper
    (and optionally lower) quantile spectra and combined conservatively via intersection
    of line-free channels. This is intentionally biased toward flagging more candidate line
    channels so that spatially sparse features are less likely to leak into the baseline fit.
    """
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

    per_lo, per_hi = rcfg.period_range_chan
    if (per_lo <= 0) or (per_hi <= 0) or (per_lo > per_hi):
        raise ValueError(f"Invalid period_range_chan={rcfg.period_range_chan}; expected positive (lo<=hi).")
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
) -> np.ndarray:
    """Solve weighted least-squares for a batch of spectra using QR when possible.

    Parameters
    ----------
    A : ndarray, shape (n_samples, n_coef)
    Y : ndarray, shape (n_samples, n_spec)
    W : ndarray, shape (n_samples, n_spec)
    rcond : float or None
        rcond used only in per-spectrum fallback ``np.linalg.lstsq`` solves.
    """
    A64 = np.asarray(A, dtype=np.float64)
    Y64 = np.asarray(Y, dtype=np.float64)
    W64 = np.asarray(W, dtype=np.float64)
    if A64.ndim != 2 or Y64.ndim != 2 or W64.shape != Y64.shape:
        raise ValueError('A/Y/W shape mismatch in _solve_weighted_qr_batch')

    sqrtW = np.sqrt(np.clip(W64, 0.0, None))
    # Stack spectra on the leading axis for batched QR: (n_spec, n_samples, n_coef).
    Aw = sqrtW.T[:, :, None] * A64[None, :, :]
    bw = sqrtW.T[:, :, None] * Y64.T[:, :, None]

    try:
        Q, R = np.linalg.qr(Aw, mode='reduced')
        Qtb = np.matmul(np.swapaxes(Q, -1, -2), bw)
        coef = np.linalg.solve(R, Qtb)[..., 0].T
        return np.asarray(coef, dtype=dtype)
    except Exception:
        n_spec = Y64.shape[1]
        n_coef = A64.shape[1]
        coef = np.empty((n_coef, n_spec), dtype=np.float64)
        for ib in range(n_spec):
            Ai = Aw[ib]
            bi = bw[ib, :, 0]
            try:
                Qi, Ri = np.linalg.qr(Ai, mode='reduced')
                coef[:, ib] = np.linalg.solve(Ri, Qi.T @ bi)
            except Exception:
                coef[:, ib], *_ = np.linalg.lstsq(Ai, bi, rcond=rcond)
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
            Aw = A_lf * w[:, None]
            bw = ylf * w
            cj_new, *_ = np.linalg.lstsq(Aw, bw, rcond=rcond)
            cj_new = np.asarray(cj_new, dtype=dtype)
            if early and bool(_coef_converged_mask(cj, cj_new, rtol=robust_coef_rtol, atol=robust_coef_atol)[0]):
                cj = cj_new
                break
            cj = cj_new
    base = A_full @ cj
    return np.asarray(yj - base, dtype=np.float32)


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
    m = np.asarray(linefree_mask, dtype=bool)
    if m.shape != (nchan,):
        raise ValueError(f"linefree_mask shape mismatch: {m.shape} vs ({nchan},)")

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
        "MOMENT0", "MOM0_BASESUP", "MOM0_LINECAND",
    }
    excluded_btypes = {
        "SIGNALMASK", "BASELINESUPPORTMASK", "LINECANDIDATEMASK",
        "MOMENT0", "MOMENT0BASELINESUPPORT", "MOMENT0LINECANDIDATE",
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


def subtract_baseline_from_fits(
    input_fits: str,
    output_fits: str,
    *,
    cube_ext: Optional[Union[int, str]] = None,
    # line-free
    linefree_cfg: LineFreeConfig = LineFreeConfig(),
    linefree_mask: Optional[np.ndarray] = None,
    manual_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]] = None,
    linefree_mode: str = "auto",
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
        for name in ("LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR"):
            if name in hdul_local:
                arr = np.squeeze(np.asarray(hdul_local[name].data))
                if arr.ndim == 1:
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

    with fits.open(input_fits, mode="readonly", memmap=True) as hdul_in:
        idx, hdu_cube = _get_cube_hdu(hdul_in, cube_ext)
        data_raw = np.asarray(hdu_cube.data, dtype=np.float32)
        data, axis_order_in = _standardize_cube_for_processing(data_raw, hdu_cube.header)
        nchan, ny, nx = data.shape
        if write_profile:
            _tprev = _profile_record(_profile_records, "read_and_standardize", _tprev, _t0, cube_shape=[int(nchan), int(ny), int(nx)], axis_order_in=str(axis_order_in))

        lf_prior = _read_linefree_prior(hdul_in) if load_prior_from_input else None
        freqs_prior = _read_ripple_prior(hdul_in) if load_prior_from_input else None

        if lf_prior is not None and np.asarray(lf_prior).shape != (nchan,):
            logging.warning(
                "Ignoring prior LINEFREE because shape=%s does not match nchan=%d.",
                np.asarray(lf_prior).shape,
                nchan,
            )
            lf_prior = None

        if linefree_mask is not None:
            lf = np.asarray(linefree_mask, dtype=bool)
            if lf.shape != (nchan,):
                raise ValueError(f"linefree_mask shape mismatch: {lf.shape} vs ({nchan},)")
        else:
            lf_auto: Optional[np.ndarray] = None
            if linefree_mode in ("auto", "or") or (linefree_mode == "prior" and lf_prior is None):
                logging.info("Estimating global line-free mask from cube.")
                lf_auto = estimate_linefree_mask_from_cube(data, cfg=linefree_cfg)
            if linefree_mode == "prior":
                lf = np.asarray(lf_prior, dtype=bool) if lf_prior is not None else np.asarray(lf_auto, dtype=bool)
            elif linefree_mode == "or":
                if lf_prior is None and lf_auto is None:
                    lf = np.ones(nchan, dtype=bool)
                elif lf_prior is None:
                    lf = np.asarray(lf_auto, dtype=bool)
                elif lf_auto is None:
                    lf = np.asarray(lf_prior, dtype=bool)
                else:
                    lf = np.asarray(lf_prior, dtype=bool) | np.asarray(lf_auto, dtype=bool)
            elif linefree_mode == "auto":
                lf = np.asarray(lf_auto if lf_auto is not None else np.ones(nchan, dtype=bool), dtype=bool)
            else:
                raise ValueError(f"Unknown linefree_mode: {linefree_mode}")
        if manual_v_windows:
            if SpectralCube is None:
                raise ValueError("manual_v_windows requires spectral_cube to read the spectral axis in km/s.")
            sc = SpectralCube.read(input_fits, hdu=(cube_ext if cube_ext is not None else idx))
            try:
                v_axis = np.asarray(sc.with_spectral_unit(u.km / u.s, velocity_convention="radio").spectral_axis.value, dtype=float)
            except Exception:
                v_axis = np.asarray(sc.spectral_axis.to(u.km / u.s).value, dtype=float)
            if v_axis.shape != (nchan,):
                raise ValueError(f"manual_v_windows spectral axis shape mismatch: {v_axis.shape} vs ({nchan},)")
            manual_signal = _manual_signal_mask_from_axis(v_axis, manual_v_windows)
            lf = np.asarray(lf, dtype=bool) & (~manual_signal)

        logging.info("Line-free fraction: %.3f", lf.mean())
        if write_profile:
            _tprev = _profile_record(_profile_records, "prepare_linefree", _tprev, _t0, linefree_fraction=float(np.mean(lf)))

        freqs: List[float] = []
        if reproducible_mode is not None:
            baseline_cfg = BaselineConfig(
                poly_order=int(baseline_cfg.poly_order),
                ripple=bool(baseline_cfg.ripple),
                robust=bool(baseline_cfg.robust),
                rcond=baseline_cfg.rcond,
                chunk_pix=int(baseline_cfg.chunk_pix),
                reproducible_mode=bool(reproducible_mode),
                compute_dtype=("float64" if bool(reproducible_mode) else str(getattr(baseline_cfg, "compute_dtype", "float32"))),
                normalize_x=(bool(reproducible_mode) or bool(getattr(baseline_cfg, "normalize_x", False))),
                strict_failures=bool(getattr(baseline_cfg, "strict_failures", False)),
                fallback_to_pixelwise=bool(getattr(baseline_cfg, "fallback_to_pixelwise", True)),
            )

        if baseline_cfg.ripple:
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
                freqs = estimate_ripple_frequencies_fft(spec, lf, rcfg=ripple_cfg, poly_order_pre=baseline_cfg.poly_order)
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
            "LINEFREE_PRIOR",
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
            "BASESUP3D",
            "LINECAND3D",
            "MOM0_BASESUP",
            "MOM0_LINECAND",
        ],
        protect_hdu=target_cube_hdu,
    )

    _strip_checksum_all_hdus(hdul_out)

    if add_qc_hdus:
        hdr_lf = fits.Header()
        hdr_lf["BTYPE"] = "LineFreeMask"
        hdr_lf["COMMENT"] = "1=line-free channel used for baseline fitting; 0=line/ignored."
        lf_hdu = fits.ImageHDU(data=lf.astype(np.uint8), header=hdr_lf, name="LINEFREE")
        _replace_or_append_hdu(hdul_out, lf_hdu)
        lf_used_hdu = fits.ImageHDU(data=lf.astype(np.uint8), header=hdr_lf.copy(), name="LINEFREE_USED")
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
                "manual_v_windows": list(manual_v_windows) if manual_v_windows is not None else None,
                "baseline_cfg": baseline_cfg.__dict__,
                "linefree_cfg": linefree_cfg.__dict__,
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
