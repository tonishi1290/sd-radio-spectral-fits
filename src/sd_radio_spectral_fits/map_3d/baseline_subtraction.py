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

import logging
import warnings
from dataclasses import dataclass
from typing import Any, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.ndimage as ndi
from astropy.io import fits
import astropy.units as u

try:
    from spectral_cube import SpectralCube
except Exception:  # pragma: no cover
    SpectralCube = None  # type: ignore


__all__ = [
    "LineFreeConfig",
    "RippleConfig",
    "BaselineConfig",
    "estimate_linefree_mask_1d",
    "estimate_linefree_mask_from_cube",
    "estimate_ripple_frequencies_fft",
    "subtract_baseline_cube",
    "subtract_baseline_from_fits",
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
    """
    smooth_width: int = 51
    sigma: float = 4.0
    iters: int = 6
    pad_chan: int = 3
    min_linefree_frac: float = 0.35


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
    """
    nfreq: int = 2
    period_range_chan: Tuple[float, float] = (20.0, 400.0)
    min_separation: float = 0.002
    window: str = "hann"


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
        If True, apply a lightweight robust reweighting per pixel (slower). Default False for huge cubes.
    rcond : float or None
        rcond for numpy.linalg.lstsq.
    chunk_pix : int
        Pixel chunk size when processing large cubes in memory.
    """
    poly_order: int = 1
    ripple: bool = True
    robust: bool = False
    rcond: Optional[float] = None
    chunk_pix: int = 65536


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
    y = np.asarray(spectrum, dtype=np.float32)
    n = y.size
    if n < 8:
        return np.isfinite(y)

    base = _median_filter_1d(np.nan_to_num(y, nan=0.0), cfg.smooth_width)
    resid = y - base

    good = np.isfinite(resid).copy()
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
    seed: int = 0,
) -> np.ndarray:
    """
    Estimate a global line-free mask from a 3D cube by aggregating spectra spatially.

    Parameters
    ----------
    cube_data : np.ndarray
        (nchan, ny, nx) cube.
    agg : {'median','mean'}
        Spatial aggregation method.
    max_pix : int
        If ny*nx is huge, randomly sample up to max_pix pixels for aggregation.
    seed : int
        RNG seed for sampling.

    Returns
    -------
    linefree : np.ndarray (bool), shape (nchan,)
    """
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(f"cube_data must be 3D (nchan,ny,nx), got shape={data.shape}")

    nchan, ny, nx = data.shape
    npix = ny * nx
    flat = data.reshape(nchan, npix)

    if npix > int(max_pix):
        rng = np.random.default_rng(int(seed))
        idx = rng.choice(npix, size=int(max_pix), replace=False)
        flat_s = flat[:, idx]
    else:
        flat_s = flat

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if agg == "median":
            spec = np.nanmedian(flat_s, axis=1)
        elif agg == "mean":
            spec = np.nanmean(flat_s, axis=1)
        else:
            raise ValueError(f"Unknown agg={agg!r}")

    return estimate_linefree_mask_1d(spec, cfg=cfg)


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
    - Zero-out line channels (non line-free)
    - Apply window (hann) to reduce leakage
    - FFT power spectrum
    - Pick top-N peaks within period_range_chan, enforcing min_separation

    Returns
    -------
    freqs : list of float
        Frequencies in cycles/channel (0..0.5). Length <= rcfg.nfreq.
    """
    y = np.asarray(spectrum, dtype=np.float32)
    m = np.asarray(linefree_mask, dtype=bool)
    n = y.size
    if n < 16:
        return []

    x = np.arange(n, dtype=np.float32)
    # polynomial pre-fit on line-free
    if poly_order_pre >= 0 and m.sum() >= (poly_order_pre + 2):
        V = np.vstack([x[m] ** k for k in range(poly_order_pre + 1)]).T  # (nlf, p)
        coeff, *_ = np.linalg.lstsq(V, y[m], rcond=None)
        base = np.zeros_like(y)
        for k in range(poly_order_pre + 1):
            base += float(coeff[k]) * (x ** k)
        r = y - base
    else:
        r = y.copy()

    r2 = np.zeros_like(r)
    r2[m] = r[m]

    if rcfg.window.lower() == "hann":
        w = np.hanning(n).astype(np.float32)
        r2 = r2 * w

    # FFT
    spec = np.fft.rfft(r2)
    pwr = (spec.real ** 2 + spec.imag ** 2).astype(np.float64)
    freqs = np.fft.rfftfreq(n, d=1.0)  # cycles per channel

    per_lo, per_hi = rcfg.period_range_chan
    fmin = 1.0 / float(per_hi)
    fmax = 1.0 / float(per_lo)
    sel = (freqs >= fmin) & (freqs <= fmax)
    sel[0] = False  # ignore DC
    if sel.sum() < 1:
        return []

    cand_idx = np.where(sel)[0]
    # Sort by descending power
    cand_idx = cand_idx[np.argsort(pwr[cand_idx])[::-1]]

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
    x: np.ndarray,
    *,
    poly_order: int,
    ripple_freqs: Optional[Sequence[float]] = None,
) -> np.ndarray:
    """
    Build design matrix A(x) for model:
      baseline(x) = sum_{k=0..poly_order} c_k x^k  +  sum_j (a_j sin(2π f_j x) + b_j cos(2π f_j x))

    Returns A with shape (len(x), ncoef).
    """
    x = np.asarray(x, dtype=np.float32)
    cols: List[np.ndarray] = []
    for k in range(int(poly_order) + 1):
        cols.append(x ** k)

    if ripple_freqs:
        for f in ripple_freqs:
            w = 2.0 * np.pi * float(f)
            cols.append(np.sin(w * x))
            cols.append(np.cos(w * x))

    return np.vstack(cols).T.astype(np.float32, copy=False)


def _robust_reweight(
    resid: np.ndarray,
    *,
    c: float = 1.345,
) -> np.ndarray:
    """
    Huber weights for residuals (vector).
    """
    r = np.asarray(resid, dtype=np.float32)
    s = _robust_std(r)
    if s <= 0:
        return np.ones_like(r, dtype=np.float32)
    t = np.abs(r) / (c * s)
    w = np.ones_like(r, dtype=np.float32)
    m = t > 1
    w[m] = 1.0 / t[m]
    return w


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
        2D uint8 flag map: 0=OK, 1=insufficient points, 2=lstsq failed.
    """
    data = np.asarray(cube_data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError("cube_data must be 3D (nchan,ny,nx)")
    nchan, ny, nx = data.shape
    m = np.asarray(linefree_mask, dtype=bool)
    if m.shape != (nchan,):
        raise ValueError(f"linefree_mask shape mismatch: {m.shape} vs ({nchan},)")

    # Fill NaNs (keeps global mask approach viable)
    data2 = _fill_nan_with_median_along_spec(data)

    # Flatten spatial dims
    npix = ny * nx
    Y = data2.reshape(nchan, npix)

    # Prepare design matrices
    x = np.arange(nchan, dtype=np.float32)
    freqs = list(ripple_freqs) if (bcfg.ripple and ripple_freqs) else []
    A_full = _design_matrix(x, poly_order=bcfg.poly_order, ripple_freqs=freqs)          # (nchan, ncoef)
    A_lf = A_full[m, :]                                                                 # (nlf, ncoef)
    nlf = int(m.sum())
    ncoef = A_full.shape[1]

    if nlf < max(ncoef + 2, 8):
        raise ValueError(f"Too few line-free points for model: nlf={nlf}, ncoef={ncoef}")

    out = np.empty_like(Y, dtype=np.float32)
    resid_rms = np.full(npix, np.nan, dtype=np.float32) if return_qc else None
    flag = np.zeros(npix, dtype=np.uint8) if return_qc else None

    # Batch processing in chunks of pixels
    chunk = int(max(1, bcfg.chunk_pix))
    for p0 in range(0, npix, chunk):
        p1 = min(npix, p0 + chunk)
        Yc = Y[:, p0:p1]            # (nchan, k)
        Ylf = Yc[m, :]              # (nlf, k)

        try:
            if not bcfg.robust:
                # Multiple RHS least squares
                coef, *_ = np.linalg.lstsq(A_lf, Ylf, rcond=bcfg.rcond)  # (ncoef, k)
            else:
                # Robust: per-pixel IRLS (slower but safer when line mask is imperfect)
                k = p1 - p0
                coef = np.zeros((ncoef, k), dtype=np.float32)
                for j in range(k):
                    yj = Ylf[:, j]
                    # init
                    cj, *_ = np.linalg.lstsq(A_lf, yj, rcond=bcfg.rcond)
                    cj = cj.astype(np.float32, copy=False)
                    for _ in range(3):
                        rj = yj - (A_lf @ cj)
                        w = _robust_reweight(rj)
                        Aw = A_lf * w[:, None]
                        bw = yj * w
                        cj, *_ = np.linalg.lstsq(Aw, bw, rcond=bcfg.rcond)
                        cj = cj.astype(np.float32, copy=False)
                    coef[:, j] = cj

            base = A_full @ coef   # (nchan, k)
            out[:, p0:p1] = Yc - base

            if return_qc:
                r = out[m, p0:p1]
                resid_rms[p0:p1] = np.sqrt(np.mean(r * r, axis=0)).astype(np.float32, copy=False)
        except Exception:
            out[:, p0:p1] = Yc
            if return_qc:
                flag[p0:p1] = 2

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
    linefree_mode: str = "auto",
    load_prior_from_input: bool = True,
    # ripple
    ripple_cfg: RippleConfig = RippleConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    ripple_mode: str = "auto",
    # baseline model
    baseline_cfg: BaselineConfig = BaselineConfig(),
    # output
    add_qc_hdus: bool = True,
    overwrite: bool = True,
) -> None:
    """
    Read cube FITS, estimate/reuse line-free + ripple frequencies, subtract baseline, write new FITS.

    Outputs (if add_qc_hdus=True)
    ----------------------------
    - LINEFREE (ImageHDU): uint8 (nchan,)  1=line-free, 0=line
    - RIPFREQ  (BinTable): ripple frequencies and periods (if ripple used)
    - BASE_RMS (ImageHDU): 2D residual RMS on line-free channels
    - BASE_FLG (ImageHDU): 2D flag map (0 ok, 2 fit failed)

    Notes
    -----
    - Baseline subtraction is done in channel space.
    - If load_prior_from_input=True, existing LINEFREE / RIPFREQ HDUs are used as priors
      according to linefree_mode / ripple_mode.
    """
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
        logging.info("Line-free fraction: %.3f", lf.mean())

        freqs: List[float] = []
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

        out_cube, rms_map, flag_map = subtract_baseline_cube(
            data,
            linefree_mask=lf,
            bcfg=baseline_cfg,
            ripple_freqs=freqs,
            return_qc=add_qc_hdus,
        )

        hdul_out = fits.HDUList([h.copy() for h in hdul_in])

    out_cube_write = _restore_cube_axis_order(out_cube, axis_order_in)
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
        _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=lf.astype(np.uint8), header=hdr_lf, name="LINEFREE"))

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

        if rms_map is not None:
            hdr = fits.Header()
            hdr["BTYPE"] = "BaselineResidRMS"
            hdr["BUNIT"] = "K"
            _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=_as_float32(rms_map), header=hdr, name="BASE_RMS"))

        if flag_map is not None:
            hdr = fits.Header()
            hdr["BTYPE"] = "BaselineFitFlag"
            hdr["COMMENT"] = "0=OK, 2=lstsq failed in chunk (spectrum left unchanged)."
            _replace_or_append_hdu(hdul_out, fits.ImageHDU(data=np.asarray(flag_map, dtype=np.uint8), header=hdr, name="BASE_FLG"))

    hdul_out.writeto(output_fits, overwrite=overwrite)
    logging.info("Wrote baselined cube FITS: %s", output_fits)
