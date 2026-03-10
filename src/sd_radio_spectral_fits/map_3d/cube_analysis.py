# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map.cube_analysis

spectral-cube を用いた 3D マスク生成 + Moment0 生成 + FITS 追記ユーティリティ。

本モジュールは spectral-cube の「軸順序は常に (n_spectral, n_y, n_x)」に揃えます。
内部計算は numpy ndarray（unitless）に落としてから行い、Quantity/単位の混在による例外を回避します。

設計方針（重要）
- 速度軸変換・モーメント計算・WCS生成は spectral-cube に委譲する
- ただし速度変換は rest_value が必要になる場合があるため、ヘッダの RESTFRQ/RESTFREQ を拾うか明示指定を許可する
- km/s 変換に失敗した場合（require_kms=False）、以後の slab 解釈は「実際に残ったスペクトル軸単位」で行う
- FITS 追記時は NAXIS/BSCALE/BZERO 等の衝突を避けるため、ヘッダをクリーニングしてから新 HDU を作る
- 既存HDUsに CHECKSUM/DATASUM がある場合、追記後は不整合になり得るため、デフォルトで削除する（strip_checksum=True）

マスク生成 (method)
- 'smooth_mask'      : 高SNRコア + 低SNR拡張 → 3D label で「コアを含む成分」を選択（重いが厳密）
- 'smooth_mask_lite' : 3D label を避ける軽量版。binary_propagation による「コアから低SNRへ伝播」で復元（巨大キューブ向け）
- 'simple'           : rms_map を閾値に単純比較
- 'derivative'       : スペクトル2階微分を閾値に検出

依存
- spectral-cube
- astropy
- numpy
- scipy
"""

from __future__ import annotations

import logging
import warnings
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.ndimage as ndi
from astropy.io import fits
import astropy.units as u

from spectral_cube import SpectralCube


__all__ = [
    "estimate_robust_rms",
    "generate_cube_mask",
    "append_analysis_hdus_to_fits",
    "append_provisional_baseline_hdus_to_fits",
    "make_3d_mask_for_existing_fits",
]


# -----------------------------------------------------------------------------
# 0. Internal helpers
# -----------------------------------------------------------------------------

def _as_float_ndarray(x: Any, dtype: np.dtype = np.float32) -> np.ndarray:
    """Convert input to a float ndarray (copy avoided when possible)."""
    arr = np.asarray(x)
    if not np.issubdtype(arr.dtype, np.floating):
        arr = arr.astype(dtype)
    else:
        arr = arr.astype(dtype, copy=False)
    return arr


def _quantity_to_value(x: Any) -> Any:
    """If x is an astropy Quantity, return x.value, otherwise return x as-is."""
    return getattr(x, "value", x)


def _clean_header_for_image_like_hdu(header: fits.Header) -> fits.Header:
    """Remove keywords that can conflict with new HDU creation."""
    hdr = header.copy()

    for key in ("SIMPLE", "EXTEND", "BITPIX", "CHECKSUM", "DATASUM", "XTENSION", "PCOUNT", "GCOUNT", "EXTNAME"):
        if key in hdr:
            del hdr[key]

    if "NAXIS" in hdr:
        naxis = int(hdr.get("NAXIS", 0))
        del hdr["NAXIS"]
        for i in range(1, naxis + 1):
            k = f"NAXIS{i}"
            if k in hdr:
                del hdr[k]
    else:
        for i in range(1, 8):
            k = f"NAXIS{i}"
            if k in hdr:
                del hdr[k]

    for key in ("BSCALE", "BZERO", "BLANK", "DATAMIN", "DATAMAX"):
        if key in hdr:
            del hdr[key]

    # Remove tile-compression bookkeeping if the source header came from a CompImageHDU.
    # Keeping these cards on newly created image HDUs can leave inconsistent metadata.
    compression_keys = []
    for key in list(hdr.keys()):
        k = str(key).upper()
        if (
            k in {"ZIMAGE", "ZCMPTYPE", "ZBITPIX", "ZNAXIS", "ZHECKSUM", "ZDATASUM"}
            or k.startswith("ZNAXIS")
            or k.startswith("ZTILE")
            or k.startswith("ZNAME")
            or k.startswith("ZVAL")
            or k.startswith("ZQUANTIZ")
            or k.startswith("ZDITHER")
            or k.startswith("ZMASKCMP")
        ):
            compression_keys.append(key)
    for key in compression_keys:
        if key in hdr:
            del hdr[key]

    return hdr


def _strip_checksum_all_hdus(hdul: fits.HDUList) -> None:
    """Remove CHECKSUM/DATASUM from all HDUs (in-place)."""
    for hdu in hdul:
        hdr = getattr(hdu, "header", None)
        if hdr is None:
            continue
        for key in ("CHECKSUM", "DATASUM", "ZHECKSUM", "ZDATASUM"):
            if key in hdr:
                del hdr[key]


def _find_existing_hdu_name(hdul: fits.HDUList, candidates: Sequence[str]) -> Optional[str]:
    """Return the first existing HDU name from candidates."""
    existing = {str(getattr(hdu, "name", "")).upper() for hdu in hdul}
    for name in candidates:
        if str(name).upper() in existing:
            return str(name)
    return None


def _append_or_replace_hdu(hdul: fits.HDUList, hdu: Union[fits.ImageHDU, fits.CompImageHDU, fits.BinTableHDU]) -> None:
    """Replace HDU with the same EXTNAME if present, otherwise append."""
    name = str(getattr(hdu, "name", "")).upper()
    if not name:
        hdul.append(hdu)
        return
    while name in hdul:
        del hdul[name]
    hdul.append(hdu)


def _remove_named_hdus(hdul: fits.HDUList, names: Sequence[str]) -> None:
    """Remove all HDUs whose EXTNAME matches one of `names` (case-insensitive)."""
    names_up = {str(n).upper() for n in names}
    remove_idx = []
    for i, hdu in enumerate(hdul):
        if i == 0:
            continue
        hname = str(getattr(hdu, 'name', '') or '').upper()
        if hname in names_up:
            remove_idx.append(i)
    for i in reversed(remove_idx):
        del hdul[i]


def _read_linefree_mask_from_hdul(hdul: fits.HDUList) -> Optional[np.ndarray]:
    name = _find_existing_hdu_name(hdul, ("LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR"))
    if name is None:
        return None
    arr = np.asarray(hdul[name].data)
    arr = np.squeeze(arr)
    if arr.ndim != 1:
        return None
    return arr.astype(bool)


def _read_base_flag_map_from_hdul(hdul: fits.HDUList) -> Optional[np.ndarray]:
    name = _find_existing_hdu_name(hdul, ("BASE_FLG", "BASEFLAG"))
    if name is None:
        return None
    arr = np.asarray(hdul[name].data)
    arr = np.squeeze(arr)
    if arr.ndim != 2:
        return None
    return arr.astype(np.uint8)


def _broadcast_linefree_mask(linefree_mask_1d: np.ndarray, cube_shape: Tuple[int, int, int]) -> np.ndarray:
    lf = np.asarray(linefree_mask_1d, dtype=bool)
    nchan, ny, nx = cube_shape
    if lf.shape != (nchan,):
        raise ValueError(f"linefree_mask shape mismatch: {lf.shape} vs {(nchan,)}")
    return np.broadcast_to(lf[:, np.newaxis, np.newaxis], (nchan, ny, nx))


def _make_provisional_masks(
    cube_shape: Tuple[int, int, int],
    linefree_mask_1d: np.ndarray,
    base_flag_2d: Optional[np.ndarray] = None,
    *,
    good_flag_values: Sequence[int] = (0,),
) -> Tuple[np.ndarray, np.ndarray]:
    """Create baseline-support and line-candidate 3D masks from 1D line-free info."""
    basesup = _broadcast_linefree_mask(linefree_mask_1d, cube_shape)
    linecand = ~basesup

    if base_flag_2d is not None:
        good = np.isin(np.asarray(base_flag_2d), np.asarray(list(good_flag_values), dtype=np.uint8))
        if good.shape != cube_shape[1:]:
            raise ValueError(f"base_flag_2d shape mismatch: {good.shape} vs {cube_shape[1:]}")
        good3d = good[np.newaxis, :, :]
        basesup = basesup & good3d
        linecand = linecand & good3d

    return basesup.astype(np.uint8), linecand.astype(np.uint8)


def _moment0_from_mask(
    cube: SpectralCube,
    mask_3d: np.ndarray,
    *,
    how: str = "auto",
    spectral_slab: Optional[Tuple[Any, Any]] = None,
) -> Any:
    masked_cube = cube.with_mask(np.asarray(mask_3d, dtype=bool))
    if spectral_slab is not None:
        lo_q, hi_q = spectral_slab
        masked_cube = masked_cube.spectral_slab(lo_q, hi_q)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        try:
            return masked_cube.moment0(how=how)
        except TypeError:
            return masked_cube.moment0()




def _resolve_cube_ext_for_spectralcube(input_fits: str, cube_ext: Optional[Union[int, str]] = None) -> Optional[Union[int, str]]:
    """Resolve which FITS HDU should be read as the 3D cube for SpectralCube."""
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
        for ax in (1, 2, 3):
            ctype = str(hdr.get(f"CTYPE{ax}", "")).upper()
            if any(tok in ctype for tok in ("FREQ", "VRAD", "VELO", "VOPT", "WAVE", "AWAV")):
                return True
        return False

    with fits.open(input_fits, mode="readonly") as hdul:
        if cube_ext is not None:
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
            try:
                return int(cube_ext) if isinstance(cube_ext, int) else hdul.index_of(cube_ext)
            except Exception:
                return cube_ext

        if _is_3d_image(hdul[0]) and (not _looks_like_analysis_hdu(hdul[0])) and _has_spectral_axis(hdul[0]):
            return None
        for i, hdu in enumerate(hdul):
            if _is_3d_image(hdu) and (not _looks_like_analysis_hdu(hdu)) and _has_spectral_axis(hdu):
                return i

    raise ValueError(
        "No suitable 3D spectral cube found in FITS. Refusing to fall back to a 3D HDU without spectral CTYPE."
    )

def _read_spectral_cube(
    input_fits: str,
    cube_ext: Optional[Union[int, str]] = None,
    use_dask: Optional[bool] = None,
    file_format: Optional[str] = None,
) -> SpectralCube:
    """Read FITS into a SpectralCube in a signature-tolerant way."""
    read_kwargs: dict[str, Any] = {}
    if use_dask is not None:
        read_kwargs["use_dask"] = use_dask
    if file_format is not None:
        read_kwargs["format"] = file_format

    if cube_ext is None:
        return SpectralCube.read(input_fits, **read_kwargs)

    try:
        return SpectralCube.read(input_fits, hdu=cube_ext, **read_kwargs)
    except TypeError:
        return SpectralCube.read(input_fits, ext=cube_ext, **read_kwargs)


def _cube_to_numpy_filled(
    cube: SpectralCube,
    fill_value: float = np.nan,
    dtype: np.dtype = np.float32,
) -> Tuple[np.ndarray, Optional[u.Unit]]:
    """
    Return (data, unit) where data is a unitless ndarray with shape (nchan, ny, nx).

    注意: この関数は cube を numpy 配列に実体化します。
    超巨大データで Dask 遅延評価を維持したい場合は、マスク生成手法自体（smooth_mask系）を別設計にする必要があります。
    """
    cube_filled = cube.with_fill_value(fill_value)
    q = cube_filled.filled_data[:]
    unit = getattr(q, "unit", None) or getattr(cube, "unit", None)
    data = _as_float_ndarray(_quantity_to_value(q), dtype=dtype)
    return data, unit


def _header_rest_value_hz(header: fits.Header) -> Optional[u.Quantity]:
    """Try to get rest frequency from FITS header keys. Return Quantity[Hz] or None."""
    for key in ("RESTFRQ", "RESTFREQ", "RESTFREQU", "REST_FRQ", "REST_FREQ"):
        if key in header:
            try:
                return float(header[key]) * u.Hz
            except Exception:
                continue
    return None


def _maybe_convert_spectral_axis_to_kms(
    cube: SpectralCube,
    *,
    velocity_convention: str = "radio",
    rest_value: Optional[u.Quantity] = None,
    require_kms: bool = False,
) -> SpectralCube:
    """Attempt to convert/verify spectral axis to km/s using spectral-cube."""
    try:
        if rest_value is not None:
            cube2 = cube.with_spectral_unit(
                u.km / u.s,
                velocity_convention=velocity_convention,
                rest_value=rest_value,
            )
        else:
            cube2 = cube.with_spectral_unit(
                u.km / u.s,
                velocity_convention=velocity_convention,
            )
        logging.info("Spectral axis converted/verified as km/s (convention=%s).", velocity_convention)
        return cube2
    except Exception as e:
        msg = f"Could not convert spectral axis to km/s: {e}"
        if require_kms:
            raise RuntimeError(msg) from e
        logging.warning(msg)
        return cube


def _spectral_axis_unit(cube: SpectralCube) -> u.Unit:
    """Best-effort: return cube.spectral_axis.unit, fallback to u.one."""
    try:
        return cube.spectral_axis.unit
    except Exception:
        return u.one


def _spectral_axis_is_kms(cube: SpectralCube) -> bool:
    """Return True if cube.spectral_axis.unit is equivalent to km/s."""
    try:
        return cube.spectral_axis.unit.is_equivalent(u.km / u.s)
    except Exception:
        return False


def _parse_spectral_slab(
    slab: Optional[Tuple[Union[float, u.Quantity], Union[float, u.Quantity]]],
    *,
    assume_unit: u.Unit,
) -> Optional[Tuple[u.Quantity, u.Quantity]]:
    """Normalize spectral slab specification."""
    if slab is None:
        return None
    lo, hi = slab
    lo_q = lo if isinstance(lo, u.Quantity) else (float(lo) * assume_unit)
    hi_q = hi if isinstance(hi, u.Quantity) else (float(hi) * assume_unit)
    return lo_q, hi_q


# -----------------------------------------------------------------------------
# 1. Core analysis logic
# -----------------------------------------------------------------------------

def estimate_robust_rms(
    cube: SpectralCube,
    axis: int = 0,
    *,
    how: str = "auto",
    ignore_nan: Optional[bool] = None,
) -> np.ndarray:
    """
    spectral-cube の mad_std を用いた MAD ベースのロバスト RMS 推定。

    Returns
    -------
    rms_map : np.ndarray
        2D array (ny, nx), unitless (cube の強度単位に対応する値)。
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        try:
            if ignore_nan is None:
                rms_q = cube.mad_std(axis=axis, how=how)
            else:
                rms_q = cube.mad_std(axis=axis, how=how, ignore_nan=ignore_nan)
        except TypeError:
            rms_q = cube.mad_std(axis=axis)

    rms = _as_float_ndarray(_quantity_to_value(rms_q), dtype=np.float32)
    return rms


def generate_cube_mask(
    cube: SpectralCube,
    rms_map: Union[np.ndarray, Any],
    method: str = "smooth_mask",
    *,
    fill_value: float = np.nan,
    data_dtype: np.dtype = np.float32,
    **kwargs: Any,
) -> np.ndarray:
    """
    指定された手法で 3D シグナルマスク (0/1) を生成する。

    Returns
    -------
    mask_3d : np.ndarray
        uint8 array with shape (nchan, ny, nx). 1=signal, 0=non-signal.
    """
    data, _unit = _cube_to_numpy_filled(cube, fill_value=fill_value, dtype=data_dtype)
    rms2d = _as_float_ndarray(_quantity_to_value(rms_map), dtype=np.float32)

    if data.ndim != 3:
        raise ValueError(f"cube must be 3D, got {data.ndim}D with shape {data.shape}")
    if rms2d.ndim != 2:
        raise ValueError(f"rms_map must be 2D, got {rms2d.ndim}D with shape {rms2d.shape}")
    if data.shape[1:] != rms2d.shape:
        raise ValueError(f"Shape mismatch: cube spatial {data.shape[1:]} != rms_map {rms2d.shape}")

    rms3d = rms2d[np.newaxis, :, :]

    # Fill NaNs with per-pixel median along spectral axis (stable for derivative-based ops)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        med = np.nanmedian(data, axis=0, keepdims=True)
    med = np.nan_to_num(med, nan=0.0).astype(np.float32, copy=False)
    data_filled_med = np.where(np.isfinite(data), data, med)

    if method == "simple":
        sigma = float(kwargs.get("sigma", 3.0))
        with np.errstate(invalid="ignore"):
            mask_final = data_filled_med > (sigma * rms3d)

    elif method == "smooth_mask":
        # Heavy but strict: 3D smoothing + (high/low) + connected-component selection using 3D labels
        smooth_sigma = kwargs.get("smooth_sigma", (2.0, 1.0, 1.0))  # (spectral, y, x)
        high_snr = float(kwargs.get("high_snr", 3.0))
        low_snr = float(kwargs.get("low_snr", 1.5))
        min_vol = int(kwargs.get("min_vol", 27))
        eps = float(kwargs.get("den_eps", 1e-4))

        data0 = np.nan_to_num(data, nan=0.0).astype(np.float32, copy=False)
        w = np.isfinite(data).astype(np.float32)

        num = ndi.gaussian_filter(data0 * w, sigma=smooth_sigma)
        den = ndi.gaussian_filter(w, sigma=smooth_sigma)

        smoothed = np.full_like(data0, np.nan, dtype=np.float32)
        valid_den = den > eps
        smoothed[valid_den] = num[valid_den] / den[valid_den]

        # Robust RMS of smoothed cube (MAD)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            sm_med = np.nanmedian(smoothed, axis=0, keepdims=True)
            sm_mad = np.nanmedian(np.abs(smoothed - sm_med), axis=0)
        rms_sm = (sm_mad * 1.4826).astype(np.float32, copy=False)
        rms_sm3d = rms_sm[np.newaxis, :, :]

        valid_sm = np.isfinite(rms_sm3d) & (rms_sm3d > 0) & valid_den

        with np.errstate(invalid="ignore"):
            mask_core = (smoothed > (high_snr * rms_sm3d)) & valid_sm
            mask_ext = (smoothed > (low_snr * rms_sm3d)) & valid_sm

        labeled_ext, _ = ndi.label(mask_ext)
        core_labels = np.unique(labeled_ext[mask_core])
        core_labels = core_labels[core_labels > 0]
        mask_recon = np.isin(labeled_ext, core_labels)

        labeled_final, _ = ndi.label(mask_recon)
        bincounts = np.bincount(labeled_final.ravel())
        valid_labels = np.where(bincounts >= min_vol)[0]
        valid_labels = valid_labels[valid_labels > 0]
        mask_final = np.isin(labeled_final, valid_labels)

        del num, den, smoothed, labeled_ext, labeled_final, bincounts

    elif method == "smooth_mask_lite":
        # Lightweight variant:
        # - avoid 3D label arrays (int32) by using binary_propagation(core within ext)
        # - additionally avoid materializing "smoothed" cube by comparing num/den via cross-multiplication
        #
        # Trade-offs:
        # - No strict 3D per-component min_vol filtering (optional 2D area filter only)
        smooth_sigma = kwargs.get("smooth_sigma", (2.0, 1.0, 1.0))  # (spectral, y, x)
        high_snr = float(kwargs.get("high_snr", 4.0))
        low_snr = float(kwargs.get("low_snr", 2.0))
        eps = float(kwargs.get("den_eps", 1e-4))

        # Optional post-filters (cheap):
        min_area = int(kwargs.get("min_area", 0))       # spatial footprint min area (pixels), 0 disables
        min_nchan = int(kwargs.get("min_nchan", 0))     # per-pixel min channels, 0 disables

        # Propagation footprint
        # connectivity=1 (6-neighbor) is a good default; reach>1 expands propagation neighborhood
        conn = int(kwargs.get("prop_connectivity", 1))
        reach = int(kwargs.get("prop_reach", 1))
        base_struct = ndi.generate_binary_structure(rank=3, connectivity=conn)
        struct = ndi.iterate_structure(base_struct, reach) if reach > 1 else base_struct

        data0 = np.nan_to_num(data, nan=0.0).astype(np.float32, copy=False)
        w = np.isfinite(data).astype(np.float32)

        num = ndi.gaussian_filter(data0 * w, sigma=smooth_sigma)
        den = ndi.gaussian_filter(w, sigma=smooth_sigma)

        valid_den = den > eps

        # Compare smoothed = num/den > thr * rms2d  <=>  num > thr * rms2d * den
        thr_hi = high_snr * rms3d
        thr_lo = low_snr * rms3d

        with np.errstate(invalid="ignore"):
            mask_core = valid_den & (num > (thr_hi * den))
            mask_ext = valid_den & (num > (thr_lo * den))

        # Geodesic dilation: propagate core inside ext
        # Equivalent to selecting ext-components connected to core (but without integer labeling)
        mask_recon = ndi.binary_propagation(mask_core, structure=struct, mask=mask_ext)

        # Optional: spatial area filter (2D labeling is much cheaper than 3D)
        if min_area > 0:
            footprint2d = np.any(mask_recon, axis=0)
            lab2d, _ = ndi.label(footprint2d)
            bc = np.bincount(lab2d.ravel())
            keep = np.where(bc >= min_area)[0]
            keep = keep[keep > 0]
            keep2d = np.isin(lab2d, keep)
            mask_recon = mask_recon & keep2d[np.newaxis, :, :]

        # Optional: per-pixel min channel count
        if min_nchan > 0:
            nchan_map = np.sum(mask_recon, axis=0)
            keep_pix = nchan_map >= min_nchan
            mask_recon = mask_recon & keep_pix[np.newaxis, :, :]

        mask_final = mask_recon
        del num, den

    elif method == "derivative":
        sigma_v = float(kwargs.get("sigma_v", 2.0))
        deriv_snr = float(kwargs.get("deriv_snr", 3.0))
        dilation_iters = int(kwargs.get("dilation_iters", 2))

        d2 = -ndi.gaussian_filter1d(data_filled_med, sigma=sigma_v, axis=0, order=2)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            d2_med = np.nanmedian(d2, axis=0, keepdims=True)
            d2_mad = np.nanmedian(np.abs(d2 - d2_med), axis=0, keepdims=True)
        d2_rms = (d2_mad * 1.4826).astype(np.float32, copy=False)

        valid_d2 = np.isfinite(d2_rms) & (d2_rms > 0)

        with np.errstate(invalid="ignore"):
            mask_core = (d2 > (deriv_snr * d2_rms)) & valid_d2

        if dilation_iters > 0:
            struct = np.zeros((3, 3, 3), dtype=bool)
            struct[:, 1, 1] = True  # spectral axis only
            mask_final = ndi.binary_dilation(mask_core, structure=struct, iterations=dilation_iters)
        else:
            mask_final = mask_core

    else:
        raise ValueError(f"Unknown masking method: '{method}'")

    mask_final = np.asarray(mask_final, dtype=bool)
    mask_final[~np.isfinite(data)] = False
    return mask_final.astype(np.uint8)


# -----------------------------------------------------------------------------
# 2. FITS IO helpers
# -----------------------------------------------------------------------------

def append_analysis_hdus_to_fits(
    hdul: fits.HDUList,
    mask_3d: np.ndarray,
    mom0_hdu: Union[fits.ImageHDU, fits.PrimaryHDU],
    base_header: fits.Header,
    *,
    mask_compression: Optional[str] = "PLIO_1",
) -> None:
    """Append/replace final signal MASK3D and MOMENT0 HDUs."""
    mask_hdr = _clean_header_for_image_like_hdu(base_header)
    mask_hdr["BTYPE"] = "SignalMask"
    mask_hdr["BUNIT"] = ""

    mask_data = np.asarray(mask_3d, dtype=np.uint8)
    if mask_compression:
        hdu_mask = fits.CompImageHDU(
            data=mask_data,
            header=mask_hdr,
            name="MASK3D",
            compression_type=mask_compression,
        )
    else:
        hdu_mask = fits.ImageHDU(data=mask_data, header=mask_hdr, name="MASK3D")
    _append_or_replace_hdu(hdul, hdu_mask)

    mom0_hdr = _clean_header_for_image_like_hdu(mom0_hdu.header)
    mom0_hdr["BTYPE"] = "Moment0"
    mom0_data = np.asarray(mom0_hdu.data, dtype=np.float32)
    _append_or_replace_hdu(hdul, fits.ImageHDU(data=mom0_data, header=mom0_hdr, name="MOMENT0"))


def append_provisional_baseline_hdus_to_fits(
    hdul: fits.HDUList,
    *,
    base_header: fits.Header,
    basesup_3d: Optional[np.ndarray] = None,
    linecand_3d: Optional[np.ndarray] = None,
    mom0_basesup_hdu: Optional[Union[fits.ImageHDU, fits.PrimaryHDU]] = None,
    mom0_linecand_hdu: Optional[Union[fits.ImageHDU, fits.PrimaryHDU]] = None,
    mask_compression: Optional[str] = "PLIO_1",
) -> None:
    """Append/replace provisional baseline-derived masks and moments."""
    def _make_mask_hdu(name: str, data: np.ndarray, btype: str) -> Union[fits.ImageHDU, fits.CompImageHDU]:
        hdr = _clean_header_for_image_like_hdu(base_header)
        hdr["BTYPE"] = btype
        hdr["BUNIT"] = ""
        if mask_compression:
            return fits.CompImageHDU(data=np.asarray(data, dtype=np.uint8), header=hdr, name=name, compression_type=mask_compression)
        return fits.ImageHDU(data=np.asarray(data, dtype=np.uint8), header=hdr, name=name)

    def _make_mom_hdu(name: str, hdu_in: Union[fits.ImageHDU, fits.PrimaryHDU], btype: str) -> fits.ImageHDU:
        hdr = _clean_header_for_image_like_hdu(hdu_in.header)
        hdr["BTYPE"] = btype
        return fits.ImageHDU(data=np.asarray(hdu_in.data, dtype=np.float32), header=hdr, name=name)

    if basesup_3d is not None:
        _append_or_replace_hdu(hdul, _make_mask_hdu("BASESUP3D", basesup_3d, "BaselineSupportMask"))
    if linecand_3d is not None:
        _append_or_replace_hdu(hdul, _make_mask_hdu("LINECAND3D", linecand_3d, "LineCandidateMask"))
    if mom0_basesup_hdu is not None:
        _append_or_replace_hdu(hdul, _make_mom_hdu("MOM0_BASESUP", mom0_basesup_hdu, "Moment0BaselineSupport"))
    if mom0_linecand_hdu is not None:
        _append_or_replace_hdu(hdul, _make_mom_hdu("MOM0_LINECAND", mom0_linecand_hdu, "Moment0LineCandidate"))


# -----------------------------------------------------------------------------
# 3. Standalone wrapper
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

def make_3d_mask_for_existing_fits(
    input_fits: str,
    output_fits: Optional[str] = None,
    *,
    cube_ext: Optional[Union[int, str]] = None,
    rms_ext: str = "auto",
    method: str = "smooth_mask",
    mask_compression: Optional[str] = "PLIO_1",
    convert_to_kms: bool = True,
    velocity_convention: str = "radio",
    rest_value_hz: Optional[float] = None,
    require_kms: bool = False,
    file_format: Optional[str] = None,
    use_dask: Optional[bool] = None,
    fill_value: float = np.nan,
    overwrite: bool = True,
    strip_checksum: bool = True,
    # provisional baseline-derived masks
    write_provisional_masks: bool = False,
    provisional_good_flag_values: Sequence[int] = (0,),
    # RMS controls
    rms_how: str = "auto",
    rms_ignore_nan: Optional[bool] = None,
    rms_spectral_slab: Optional[Tuple[Union[float, u.Quantity], Union[float, u.Quantity]]] = None,
    # Moment0 controls
    moment_how: str = "auto",
    moment_spectral_slab: Optional[Tuple[Union[float, u.Quantity], Union[float, u.Quantity]]] = None,
    **kwargs: Any,
) -> None:
    """
    保存済みの 3D FITS を読み込み、必要に応じて provisional baseline masks、
    final signal MASK3D、および各種 moment を生成して追記し保存する。

    Parameters
    ----------
    rms_ext : str
        'auto' のときは BASE_RMS -> RMS の順で既存 RMS HDU を探す。
    write_provisional_masks : bool
        True のとき、LINEFREE/BASE_FLG から BASESUP3D / LINECAND3D / provisional moments を追加する。
    """
    if output_fits is None:
        output_fits = input_fits

    logging.info("Reading FITS with spectral-cube: %s", input_fits)
    cube_ext_resolved = _resolve_cube_ext_for_spectralcube(input_fits, cube_ext=cube_ext)
    cube = _read_spectral_cube(input_fits, cube_ext=cube_ext_resolved, use_dask=use_dask, file_format=file_format)

    rest_value: Optional[u.Quantity] = None
    if rest_value_hz is not None:
        rest_value = float(rest_value_hz) * u.Hz
    else:
        hdr = getattr(cube, "header", None)
        if hdr is not None:
            rest_value = _header_rest_value_hz(hdr)
        if rest_value is None:
            with fits.open(input_fits, mode="readonly") as hdul_tmp:
                if cube_ext_resolved is None:
                    base_hdr_tmp = hdul_tmp[0].header
                else:
                    try:
                        base_hdr_tmp = hdul_tmp[cube_ext_resolved].header
                    except Exception:
                        base_hdr_tmp = hdul_tmp[0].header
                rest_value = _header_rest_value_hz(base_hdr_tmp)

    if convert_to_kms:
        cube = _maybe_convert_spectral_axis_to_kms(
            cube,
            velocity_convention=velocity_convention,
            rest_value=rest_value,
            require_kms=require_kms,
        )

    spec_unit = _spectral_axis_unit(cube)
    spec_is_kms = _spectral_axis_is_kms(cube)
    if convert_to_kms and (not spec_is_kms):
        logging.warning(
            "Spectral axis is not km/s after conversion attempt (unit=%s). Slabs will use this unit.",
            spec_unit,
        )
    assume_unit = (u.km / u.s) if spec_is_kms else spec_unit

    cube_for_rms = cube
    rms_slab_q = _parse_spectral_slab(rms_spectral_slab, assume_unit=assume_unit)
    if rms_slab_q is not None:
        lo_q, hi_q = rms_slab_q
        try:
            cube_for_rms = cube.spectral_slab(lo_q, hi_q)
            logging.info("RMS estimated in spectral_slab: [%s, %s].", lo_q, hi_q)
        except Exception as e:
            logging.warning("Could not apply rms_spectral_slab (%s, %s): %s", lo_q, hi_q, e)

    with fits.open(input_fits, mode="readonly") as hdul_in:
        hdul_out = fits.HDUList([hdu.copy() for hdu in hdul_in])
        linefree_prior = _read_linefree_mask_from_hdul(hdul_in)
        base_flag_2d = _read_base_flag_map_from_hdul(hdul_in)

        rms_map: Optional[np.ndarray] = None
        rms_candidates: List[str]
        if str(rms_ext).lower() == "auto":
            # Prefer BASE_RMS for baseline-subtracted cubes. A copied/stale gridding RMS
            # can coexist with a fresh BASE_RMS after baseline subtraction.
            rms_candidates = ["BASE_RMS", "RMS"]
        else:
            rms_candidates = [str(rms_ext)]

        for name in rms_candidates:
            if name in hdul_in:
                logging.info("Using existing RMS map from FITS HDU '%s'.", name)
                rms_data = np.squeeze(hdul_in[name].data)
                try:
                    if rms_data.ndim == 2 and tuple(rms_data.shape) == tuple(cube.shape[1:]):
                        rms_map = _as_float_ndarray(rms_data, dtype=np.float32)
                        break
                    logging.warning(
                        "FITS RMS shape mismatch or invalid: rms=%s cube_spatial=%s. Recalculating...",
                        getattr(rms_data, "shape", None),
                        cube.shape[1:],
                    )
                except Exception:
                    logging.warning("FITS RMS invalid. Recalculating...")

    if rms_map is None:
        logging.info("Estimating robust RMS from cube via spectral-cube (MAD std).")
        rms_map = estimate_robust_rms(cube_for_rms, axis=0, how=rms_how, ignore_nan=rms_ignore_nan)

    logging.info("Generating 3D mask (method=%s).", method)
    mask_3d = generate_cube_mask(cube, rms_map, method=method, fill_value=fill_value, **kwargs)

    masked_cube = cube.with_mask(mask_3d.astype(bool))
    cube_for_mom = masked_cube
    mom_slab_q = _parse_spectral_slab(moment_spectral_slab, assume_unit=assume_unit)
    if mom_slab_q is not None:
        lo_q, hi_q = mom_slab_q
        try:
            cube_for_mom = masked_cube.spectral_slab(lo_q, hi_q)
            logging.info("Moment0 integrated in spectral_slab: [%s, %s].", lo_q, hi_q)
        except Exception as e:
            logging.warning("Could not apply moment_spectral_slab (%s, %s): %s", lo_q, hi_q, e)

    logging.info("Calculating Moment0 via spectral-cube (how=%s).", moment_how)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        try:
            mom0 = cube_for_mom.moment0(how=moment_how)
        except TypeError:
            mom0 = cube_for_mom.moment0()

    if strip_checksum:
        _strip_checksum_all_hdus(hdul_out)

    base_header = getattr(cube, "header", None)
    if base_header is None:
        with fits.open(input_fits, mode="readonly") as hdul_in:
            base_header = hdul_in[0].header.copy()
    else:
        base_header = base_header.copy()

    # Always clear previous provisional products first so reruns cannot keep stale results.
    _remove_named_hdus(hdul_out, ["BASESUP3D", "LINECAND3D", "MOM0_BASESUP", "MOM0_LINECAND"])
    if write_provisional_masks and (linefree_prior is not None):
        try:
            basesup_3d, linecand_3d = _make_provisional_masks(
                tuple(int(v) for v in cube.shape),
                linefree_prior,
                base_flag_2d=base_flag_2d,
                good_flag_values=provisional_good_flag_values,
            )
            # Keep provisional moments consistent with the final MOMENT0 integration slab/window.
            mom0_basesup = _moment0_from_mask(cube, basesup_3d, how=moment_how, spectral_slab=mom_slab_q)
            mom0_linecand = _moment0_from_mask(cube, linecand_3d, how=moment_how, spectral_slab=mom_slab_q)
            append_provisional_baseline_hdus_to_fits(
                hdul_out,
                base_header=base_header,
                basesup_3d=basesup_3d,
                linecand_3d=linecand_3d,
                mom0_basesup_hdu=mom0_basesup.hdu,
                mom0_linecand_hdu=mom0_linecand.hdu,
                mask_compression=mask_compression,
            )
        except Exception as e:
            logging.warning("Could not build provisional baseline-derived masks: %s", e)

    append_analysis_hdus_to_fits(
        hdul_out,
        mask_3d,
        mom0.hdu,
        base_header=base_header,
        mask_compression=mask_compression,
    )

    logging.info("Saving updated FITS to: %s", output_fits)
    hdul_out.writeto(output_fits, overwrite=overwrite)
