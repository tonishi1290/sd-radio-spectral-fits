from __future__ import annotations

import builtins
import warnings
from dataclasses import dataclass

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy import fft as spfft
from scipy.ndimage import distance_transform_edt
from scipy.signal.windows import tukey

from .otf_bundle import OTFBundle
from .otf_bundle_io import validate_otf_bundle
from ..ranges import parse_windows, window_to_mask


@dataclass
class FamilyNoiseSummary:
    rms_map: np.ndarray
    sigma_rep: float
    linefree_mask_1d: np.ndarray
    sigma_ch: np.ndarray | None = None
    noise_mask_2d: np.ndarray | None = None
    noise_mask_npix: int = 0


@dataclass
class PlaitSelectionDiagnostics:
    sigma_plait: float
    sigma_average: float
    selected_method: str


def _velocity_axis_from_header(header: fits.Header, nchan: int) -> np.ndarray:
    crpix = float(header.get("CRPIX3", 1.0))
    crval = float(header.get("CRVAL3", 0.0))
    cdelt = float(header.get("CDELT3", 1.0))
    idx = np.arange(nchan, dtype=float) + 1.0
    return crval + (idx - crpix) * cdelt


def _resolve_velocity_windows_spec(windows: list[str] | list[tuple[float, float]] | None) -> list[tuple[float, float]] | None:
    if windows is None:
        return None
    if len(windows) == 0:
        return []
    first = windows[0]
    if isinstance(first, str):
        return parse_windows(windows)
    out: list[tuple[float, float]] = []
    for lo, hi in windows:
        lo_f = float(lo)
        hi_f = float(hi)
        if hi_f < lo_f:
            lo_f, hi_f = hi_f, lo_f
        out.append((lo_f, hi_f))
    return out


def _make_linefree_mask_from_velocity_windows(vel_kms: np.ndarray, windows: list[str] | list[tuple[float, float]] | None) -> np.ndarray:
    if windows is None:
        raise ValueError("linefree_velocity_windows_kms must be provided for FFT/PLAIT v1")
    resolved = _resolve_velocity_windows_spec(windows)
    mask = window_to_mask(np.asarray(vel_kms, dtype=float), resolved)
    if int(np.count_nonzero(mask)) < 3:
        raise ValueError("Too few line-free channels after applying linefree_velocity_windows_kms")
    return mask


def _robust_rms(values: np.ndarray, *, floor: float | None = None) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan
    med = float(np.nanmedian(arr))
    mad = float(np.nanmedian(np.abs(arr - med)))
    sigma = 1.4826 * mad
    if floor is not None and np.isfinite(floor):
        sigma = builtins.max(float(floor), float(sigma))
    return float(sigma)


def _robust_location(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan
    return float(np.nanmedian(arr))


def _wcs_key(header: fits.Header) -> tuple:
    keys = (
        "CTYPE1", "CTYPE2", "CTYPE3",
        "CRPIX1", "CRPIX2", "CRPIX3",
        "CRVAL1", "CRVAL2", "CRVAL3",
        "CDELT1", "CDELT2", "CDELT3",
        "CUNIT1", "CUNIT2", "CUNIT3",
    )
    return tuple(header.get(k) for k in keys)


def _normalize_support_mask(bundle: OTFBundle, *, relative_weight_threshold: float = 0.10) -> np.ndarray:
    ny, nx = bundle.data.shape[1:]
    if bundle.support_mask is not None:
        support = np.asarray(bundle.support_mask, dtype=bool).reshape(ny, nx)
    elif "SUPPORT_MASK" in bundle.image_ext:
        support = np.asarray(bundle.image_ext["SUPPORT_MASK"], dtype=bool).reshape(ny, nx)
    elif "MASK" in bundle.image_ext:
        support = np.asarray(bundle.image_ext["MASK"], dtype=bool).reshape(ny, nx)
    else:
        support = np.any(np.isfinite(bundle.data), axis=0).reshape(ny, nx)

    w = None
    for key in ("WEIGHT_SUM", "WEIGHT", "WSUM"):
        if key in bundle.image_ext:
            w = np.asarray(bundle.image_ext[key], dtype=float)
            break
    if w is not None:
        finite = w[np.isfinite(w) & (w > 0)]
        if finite.size:
            thr = float(np.nanmax(finite)) * float(relative_weight_threshold)
            support = support & np.isfinite(w) & (w >= thr)
    return support


def _normalize_valid_mask(bundle: OTFBundle) -> np.ndarray:
    nchan, ny, nx = bundle.data.shape
    if bundle.valid_mask is None:
        support = _normalize_support_mask(bundle)
        return np.isfinite(bundle.data) & support[None, :, :]
    vm = np.asarray(bundle.valid_mask, dtype=bool)
    if vm.shape == (ny, nx):
        return np.broadcast_to(vm[None, :, :], (nchan, ny, nx)).copy()
    if vm.shape == (nchan, ny, nx):
        return vm.copy()
    raise ValueError(f"Unsupported valid_mask shape={vm.shape}")


def _resolve_taper_width_pix(shape_2d: tuple[int, int], width_pix: int | None = None) -> int:
    ny, nx = int(shape_2d[0]), int(shape_2d[1])
    if width_pix is not None:
        return builtins.max(1, int(width_pix))
    # Hybrid rule: at most 5 pix and at most 10% of the smaller map size.
    return builtins.max(1, builtins.min(5, int(np.floor(0.10 * builtins.min(ny, nx)))))


def _tukey_alpha_for_size(n: int, width_pix: int) -> float:
    if n <= 1:
        return 0.0
    alpha = (2.0 * float(width_pix)) / float(builtins.max(1, n - 1))
    return float(builtins.min(1.0, builtins.max(0.0, alpha)))


def _estimate_rms_map_from_arrays(
    data: np.ndarray,
    valid_mask: np.ndarray,
    support_mask_2d: np.ndarray,
    linefree_mask_1d: np.ndarray,
) -> np.ndarray:
    # This function must see pre-taper data only.
    use = valid_mask & support_mask_2d[None, :, :] & linefree_mask_1d[:, None, None]
    if int(np.count_nonzero(linefree_mask_1d)) < 3:
        raise ValueError("At least 3 line-free channels are required for empirical RMS estimation.")
    data_lf = np.where(use, np.asarray(data, dtype=float), np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        med = np.nanmedian(data_lf, axis=0)
        mad = np.nanmedian(np.abs(data_lf - med[None, :, :]), axis=0)
    sigma = 1.4826 * mad
    cnt = np.sum(np.isfinite(data_lf), axis=0)
    sigma[cnt < 3] = np.nan
    floor_base = _robust_rms(data_lf[np.isfinite(data_lf)])
    if np.isfinite(floor_base) and floor_base > 0:
        floor = builtins.max(1.0e-12, floor_base * 1.0e-6)
        sigma = np.where(np.isfinite(sigma), np.maximum(sigma, floor), np.nan)
    return np.asarray(sigma, dtype=float)


def _estimate_noise_spatial_mask(
    rms_map: np.ndarray,
    support_mask_2d: np.ndarray,
    *,
    min_pixels: int = 8,
) -> np.ndarray:
    support = np.asarray(support_mask_2d, dtype=bool)
    rms = np.asarray(rms_map, dtype=float)
    good = support & np.isfinite(rms) & (rms > 0)
    if int(np.count_nonzero(good)) <= int(min_pixels):
        return good
    vals = rms[good]
    med = float(np.nanmedian(vals))
    mad = float(np.nanmedian(np.abs(vals - med)))
    sigma = 1.4826 * mad
    if np.isfinite(sigma) and sigma > 0:
        upper = med + 3.0 * sigma
        mask = good & (rms <= upper)
        if int(np.count_nonzero(mask)) >= int(min_pixels):
            return mask
    return good


def _estimate_channel_noise_spectrum(
    data: np.ndarray,
    valid_mask: np.ndarray,
    support_mask_2d: np.ndarray,
    noise_mask_2d: np.ndarray,
    *,
    fallback_sigma: float,
) -> np.ndarray:
    nchan = int(data.shape[0])
    support = np.asarray(support_mask_2d, dtype=bool)
    noise_mask = np.asarray(noise_mask_2d, dtype=bool) & support
    out = np.full(nchan, np.nan, dtype=float)
    floor = builtins.max(1.0e-12, float(fallback_sigma) * 1.0e-6) if np.isfinite(fallback_sigma) and fallback_sigma > 0 else None
    for k in range(nchan):
        use = np.asarray(valid_mask[k], dtype=bool) & noise_mask
        vals = np.asarray(data[k], dtype=float)[use]
        sigma = _robust_rms(vals, floor=floor)
        if np.isfinite(sigma) and sigma > 0:
            out[k] = float(sigma)
    bad = ~np.isfinite(out) | (out <= 0)
    if np.any(bad):
        out[bad] = float(fallback_sigma)
    return out


def _estimate_bundle_noise(
    bundle: OTFBundle,
    *,
    linefree_mask_1d: np.ndarray,
    plait_noise_mode: str = "family_scalar",
) -> FamilyNoiseSummary:
    data = np.asarray(bundle.data, dtype=float)
    valid = _normalize_valid_mask(bundle)
    support = _normalize_support_mask(bundle)
    rms_map = _estimate_rms_map_from_arrays(data, valid, support, linefree_mask_1d)
    sigma_rep = _robust_location(rms_map[support & np.isfinite(rms_map) & (rms_map > 0)])
    if not np.isfinite(sigma_rep) or sigma_rep <= 0:
        sigma_rep = _robust_location(rms_map[np.isfinite(rms_map) & (rms_map > 0)])
    if not np.isfinite(sigma_rep) or sigma_rep <= 0:
        raise ValueError("Could not estimate representative empirical RMS for FFT/PLAIT.")

    sigma_ch = None
    noise_mask_2d = _estimate_noise_spatial_mask(rms_map, support)
    if str(plait_noise_mode).strip().lower() == "family_channel":
        sigma_ch = _estimate_channel_noise_spectrum(
            data,
            valid,
            support,
            noise_mask_2d,
            fallback_sigma=float(sigma_rep),
        )

    return FamilyNoiseSummary(
        rms_map=rms_map,
        sigma_rep=float(sigma_rep),
        linefree_mask_1d=np.asarray(linefree_mask_1d, dtype=bool),
        sigma_ch=None if sigma_ch is None else np.asarray(sigma_ch, dtype=float),
        noise_mask_2d=np.asarray(noise_mask_2d, dtype=bool),
        noise_mask_npix=int(np.count_nonzero(noise_mask_2d)),
    )


def _estimate_output_sigma_rep(data: np.ndarray, valid_mask: np.ndarray, support_mask: np.ndarray, linefree_mask_1d: np.ndarray) -> float:
    rms_map = _estimate_rms_map_from_arrays(data, valid_mask, support_mask, linefree_mask_1d)
    sigma = _robust_location(rms_map[support_mask & np.isfinite(rms_map) & (rms_map > 0)])
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = _robust_location(rms_map[np.isfinite(rms_map) & (rms_map > 0)])
    return float(sigma)


def _build_support_taper(support_mask_2d: np.ndarray, *, width_pix: int | None = None) -> np.ndarray:
    support = np.asarray(support_mask_2d, dtype=bool)
    dist = distance_transform_edt(support)
    width = _resolve_taper_width_pix(support.shape, width_pix=width_pix)
    maxdist = int(np.floor(float(np.nanmax(dist)))) if np.any(np.isfinite(dist)) else 0
    width = builtins.max(1, builtins.min(int(width), builtins.max(1, maxdist)))
    taper = np.zeros_like(dist, dtype=float)
    inside = support & (dist >= width)
    edge = support & (~inside)
    taper[inside] = 1.0
    taper[edge] = 0.5 * (1.0 - np.cos(np.pi * dist[edge] / float(width)))
    return taper


def _build_apodization_window_2d(shape: tuple[int, int], *, alpha: float | None = 0.1, width_pix: int | None = None) -> np.ndarray:
    ny, nx = int(shape[0]), int(shape[1])
    width_y = _resolve_taper_width_pix((ny, nx), width_pix=width_pix)
    width_x = width_y
    if alpha is not None:
        width_y = builtins.min(width_y, builtins.max(1, int(round(float(alpha) * builtins.max(1, ny - 1) / 2.0))))
        width_x = builtins.min(width_x, builtins.max(1, int(round(float(alpha) * builtins.max(1, nx - 1) / 2.0))))
    ay = _tukey_alpha_for_size(ny, width_y)
    ax = _tukey_alpha_for_size(nx, width_x)
    wy = tukey(ny, alpha=ay)
    wx = tukey(nx, alpha=ax)
    return np.asarray(wy[:, None] * wx[None, :], dtype=float)


def _apply_apodization(data: np.ndarray, *, alpha: float | None = 0.1, width_pix: int | None = None) -> np.ndarray:
    window2d = _build_apodization_window_2d(data.shape[1:], alpha=alpha, width_pix=width_pix)
    return data * window2d[None, :, :]


def _build_fft_input_window_2d(
    support_mask_2d: np.ndarray,
    *,
    apodize: bool,
    apodize_alpha: float | None,
    support_taper: bool,
    support_taper_width_pix: int | None,
) -> np.ndarray:
    support = np.asarray(support_mask_2d, dtype=bool)
    window = support.astype(float)
    if support_taper:
        window = _build_support_taper(support, width_pix=support_taper_width_pix)
    if apodize:
        window = window * _build_apodization_window_2d(support.shape, alpha=apodize_alpha, width_pix=support_taper_width_pix)
    window[~support] = 0.0
    return np.asarray(window, dtype=float)


def _estimate_bw_flatfield_gain_approx(
    x_support_mask_2d: np.ndarray,
    y_support_mask_2d: np.ndarray,
    *,
    wx: np.ndarray,
    wy: np.ndarray,
    crop: tuple[slice, slice],
    pad_frac: float,
    apodize: bool,
    apodize_alpha: float | None,
    support_taper: bool,
    support_taper_width_pix: int | None,
    fft_workers: int,
    support_out_mask_2d: np.ndarray,
) -> np.ndarray:
    x_window = _build_fft_input_window_2d(
        x_support_mask_2d,
        apodize=apodize,
        apodize_alpha=apodize_alpha,
        support_taper=support_taper,
        support_taper_width_pix=support_taper_width_pix,
    )
    y_window = _build_fft_input_window_2d(
        y_support_mask_2d,
        apodize=apodize,
        apodize_alpha=apodize_alpha,
        support_taper=support_taper,
        support_taper_width_pix=support_taper_width_pix,
    )
    x_pad, crop_check = _pad_spatial(x_window[None, :, :], pad_frac=pad_frac)
    y_pad, _ = _pad_spatial(y_window[None, :, :], pad_frac=pad_frac)
    if crop_check != crop:
        raise RuntimeError("Internal crop mismatch while estimating BW flat-field gain.")
    fx = spfft.fftn(x_pad, axes=(-2, -1), workers=fft_workers)
    fy = spfft.fftn(y_pad, axes=(-2, -1), workers=fft_workers)
    out_pad = spfft.ifftn(fx * wx[None, :, :] + fy * wy[None, :, :], axes=(-2, -1), workers=fft_workers).real
    gain = np.asarray(out_pad[0, crop[0], crop[1]], dtype=float)
    support_out = np.asarray(support_out_mask_2d, dtype=bool)
    finite = gain[support_out & np.isfinite(gain)]
    if finite.size > 0:
        peak = float(np.nanmax(finite))
        if np.isfinite(peak) and peak > 0:
            gain = gain / peak
    gain[~support_out] = 0.0
    return np.clip(np.asarray(gain, dtype=float), 0.0, 1.0)


def _pad_spatial(data: np.ndarray, *, pad_frac: float = 0.25) -> tuple[np.ndarray, tuple[slice, slice]]:
    nchan, ny, nx = data.shape
    pady = builtins.max(0, int(round(ny * float(pad_frac))))
    padx = builtins.max(0, int(round(nx * float(pad_frac))))
    ny_pad = spfft.next_fast_len(ny + 2 * pady)
    nx_pad = spfft.next_fast_len(nx + 2 * padx)
    y0 = (ny_pad - ny) // 2
    x0 = (nx_pad - nx) // 2
    out = np.zeros((nchan, ny_pad, nx_pad), dtype=data.dtype)
    out[:, y0:y0 + ny, x0:x0 + nx] = data
    return out, (slice(y0, y0 + ny), slice(x0, x0 + nx))


def _build_fourier_weights(ny: int, nx: int, *, sigma_x: float, sigma_y: float) -> tuple[np.ndarray, np.ndarray]:
    fy = np.fft.fftfreq(ny)
    fx = np.fft.fftfreq(nx)
    u, v = np.meshgrid(fx, fy)
    denom = (sigma_y * sigma_y) * (u * u) + (sigma_x * sigma_x) * (v * v)
    wx = np.zeros((ny, nx), dtype=float)
    wy = np.zeros((ny, nx), dtype=float)
    good = denom > 0
    wx[good] = (sigma_y * sigma_y) * (u[good] * u[good]) / denom[good]
    wy[good] = (sigma_x * sigma_x) * (v[good] * v[good]) / denom[good]
    nu_x = 1.0 / (sigma_x * sigma_x)
    nu_y = 1.0 / (sigma_y * sigma_y)
    wx[0, 0] = nu_x / (nu_x + nu_y)
    wy[0, 0] = nu_y / (nu_x + nu_y)
    return wx.astype(float), wy.astype(float)


def _build_fourier_weights_family_channel(ny: int, nx: int, *, sigma_x_ch: np.ndarray, sigma_y_ch: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    fy = np.fft.fftfreq(ny)
    fx = np.fft.fftfreq(nx)
    u, v = np.meshgrid(fx, fy)
    u2 = u * u
    v2 = v * v
    sigma_x_ch = np.asarray(sigma_x_ch, dtype=float)
    sigma_y_ch = np.asarray(sigma_y_ch, dtype=float)
    nchan = int(sigma_x_ch.size)
    wx_mean = np.zeros((ny, nx), dtype=float)
    wy_mean = np.zeros((ny, nx), dtype=float)
    mean_wx2 = np.zeros(nchan, dtype=float)
    mean_wy2 = np.zeros(nchan, dtype=float)
    for k in range(nchan):
        sx = float(sigma_x_ch[k])
        sy = float(sigma_y_ch[k])
        denom = (sy * sy) * u2 + (sx * sx) * v2
        wx = np.zeros((ny, nx), dtype=float)
        wy = np.zeros((ny, nx), dtype=float)
        good = denom > 0
        wx[good] = (sy * sy) * u2[good] / denom[good]
        wy[good] = (sx * sx) * v2[good] / denom[good]
        nu_x = 1.0 / (sx * sx)
        nu_y = 1.0 / (sy * sy)
        wx[0, 0] = nu_x / (nu_x + nu_y)
        wy[0, 0] = nu_y / (nu_x + nu_y)
        wx_mean += wx
        wy_mean += wy
        mean_wx2[k] = float(np.mean(np.abs(wx) ** 2))
        mean_wy2[k] = float(np.mean(np.abs(wy) ** 2))
    wx_mean /= float(builtins.max(1, nchan))
    wy_mean /= float(builtins.max(1, nchan))
    return wx_mean, wy_mean, mean_wx2, mean_wy2


def _plait_combine_family_channel(
    fx: np.ndarray,
    fy: np.ndarray,
    *,
    sigma_x_ch: np.ndarray,
    sigma_y_ch: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nchan, ny, nx = fx.shape
    fy_grid = np.fft.fftfreq(ny)
    fx_grid = np.fft.fftfreq(nx)
    u, v = np.meshgrid(fx_grid, fy_grid)
    u2 = u * u
    v2 = v * v
    sigma_x_ch = np.asarray(sigma_x_ch, dtype=float)
    sigma_y_ch = np.asarray(sigma_y_ch, dtype=float)
    fout = np.empty_like(fx)
    wx_mean = np.zeros((ny, nx), dtype=float)
    wy_mean = np.zeros((ny, nx), dtype=float)
    mean_wx2 = np.zeros(nchan, dtype=float)
    mean_wy2 = np.zeros(nchan, dtype=float)
    for k in range(nchan):
        sx = float(sigma_x_ch[k])
        sy = float(sigma_y_ch[k])
        denom = (sy * sy) * u2 + (sx * sx) * v2
        wx = np.zeros((ny, nx), dtype=float)
        wy = np.zeros((ny, nx), dtype=float)
        good = denom > 0
        wx[good] = (sy * sy) * u2[good] / denom[good]
        wy[good] = (sx * sx) * v2[good] / denom[good]
        nu_x = 1.0 / (sx * sx)
        nu_y = 1.0 / (sy * sy)
        wx[0, 0] = nu_x / (nu_x + nu_y)
        wy[0, 0] = nu_y / (nu_x + nu_y)
        fout[k] = fx[k] * wx + fy[k] * wy
        wx_mean += wx
        wy_mean += wy
        mean_wx2[k] = float(np.mean(np.abs(wx) ** 2))
        mean_wy2[k] = float(np.mean(np.abs(wy) ** 2))
    wx_mean /= float(builtins.max(1, nchan))
    wy_mean /= float(builtins.max(1, nchan))
    return fout, wx_mean, wy_mean, mean_wx2, mean_wy2


def _approx_output_variance_cube(shape: tuple[int, int, int], *, sigma_x: float, sigma_y: float, wx: np.ndarray, wy: np.ndarray) -> np.ndarray:
    mean_wx2 = float(np.mean(np.abs(wx) ** 2))
    mean_wy2 = float(np.mean(np.abs(wy) ** 2))
    sigma2 = (sigma_x * sigma_x) * mean_wx2 + (sigma_y * sigma_y) * mean_wy2
    nchan, ny, nx = shape
    return np.full((nchan, ny, nx), sigma2, dtype=float)


def _approx_output_variance_cube_family_channel(
    shape: tuple[int, int, int],
    *,
    sigma_x_ch: np.ndarray,
    sigma_y_ch: np.ndarray,
    mean_wx2_ch: np.ndarray,
    mean_wy2_ch: np.ndarray,
) -> np.ndarray:
    sigma_x_ch = np.asarray(sigma_x_ch, dtype=float)
    sigma_y_ch = np.asarray(sigma_y_ch, dtype=float)
    mean_wx2_ch = np.asarray(mean_wx2_ch, dtype=float)
    mean_wy2_ch = np.asarray(mean_wy2_ch, dtype=float)
    sigma2 = (sigma_x_ch * sigma_x_ch) * mean_wx2_ch + (sigma_y_ch * sigma_y_ch) * mean_wy2_ch
    sigma2 = np.where(np.isfinite(sigma2) & (sigma2 > 0), sigma2, np.nan)
    nchan, ny, nx = shape
    out = np.broadcast_to(sigma2[:, None, None], (nchan, ny, nx)).copy()
    return out


def _family_weighted_average_output(x_data: np.ndarray, y_data: np.ndarray, x_valid: np.ndarray, y_valid: np.ndarray, *, sigma_x: float, sigma_y: float) -> tuple[np.ndarray, np.ndarray]:
    # Conservative fallback / comparator: inverse-variance family average.
    w_x = np.where(x_valid, 1.0 / float(sigma_x * sigma_x), 0.0)
    w_y = np.where(y_valid, 1.0 / float(sigma_y * sigma_y), 0.0)
    denom = w_x + w_y
    numer = np.where(x_valid, x_data * w_x, 0.0) + np.where(y_valid, y_data * w_y, 0.0)
    out = np.full_like(x_data, np.nan, dtype=float)
    good = denom > 0
    out[good] = numer[good] / denom[good]
    var = np.full_like(x_data, np.nan, dtype=float)
    var[good] = 1.0 / denom[good]
    return out, var


def _family_weighted_average_output_family_channel(
    x_data: np.ndarray,
    y_data: np.ndarray,
    x_valid: np.ndarray,
    y_valid: np.ndarray,
    *,
    sigma_x_ch: np.ndarray,
    sigma_y_ch: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    nchan = int(x_data.shape[0])
    out = np.full_like(x_data, np.nan, dtype=float)
    var = np.full_like(x_data, np.nan, dtype=float)
    sigma_x_ch = np.asarray(sigma_x_ch, dtype=float)
    sigma_y_ch = np.asarray(sigma_y_ch, dtype=float)
    for k in range(nchan):
        sx = float(sigma_x_ch[k])
        sy = float(sigma_y_ch[k])
        wx = np.where(x_valid[k], 1.0 / (sx * sx), 0.0)
        wy = np.where(y_valid[k], 1.0 / (sy * sy), 0.0)
        denom = wx + wy
        numer = np.where(x_valid[k], x_data[k] * wx, 0.0) + np.where(y_valid[k], y_data[k] * wy, 0.0)
        good = denom > 0
        out[k, good] = numer[good] / denom[good]
        var[k, good] = 1.0 / denom[good]
    return out, var




def _collect_input_family_diagnostics(x_bundle: OTFBundle, y_bundle: OTFBundle) -> dict[str, np.ndarray]:
    mapping = (
        "WEIGHT",
        "WEIGHT_SUM",
        "HIT",
        "NSAMP",
        "WSUM",
        "WABS",
        "CANCEL",
        "WREL",
        "RMS",
        "TSYS",
        "TINT",
        "TIME",
        "NEFF",
        "XEFF",
        "YEFF",
        "BIAS_PIX",
    )
    out: dict[str, np.ndarray] = {}
    for src_bundle, suffix in ((x_bundle, "X"), (y_bundle, "Y")):
        for key in mapping:
            if key in src_bundle.image_ext:
                out[f"{key}_{suffix}_IN"] = np.asarray(src_bundle.image_ext[key])
    return out
def _make_channel_noise_tables(
    *,
    velocities_kms: np.ndarray,
    linefree_mask: np.ndarray,
    noise_x: FamilyNoiseSummary,
    noise_y: FamilyNoiseSummary,
) -> dict[str, Table]:
    if noise_x.sigma_ch is None or noise_y.sigma_ch is None:
        return {}
    chan = np.arange(int(velocities_kms.size), dtype=np.int32)
    linefree_u8 = np.asarray(linefree_mask, dtype=np.uint8)
    common = {
        "CHAN": chan,
        "VELOCITY_KMS": np.asarray(velocities_kms, dtype=np.float64),
        "LINEFREE": linefree_u8,
    }
    return {
        "CHANNEL_NOISE_X": Table({
            **common,
            "SIGMA": np.asarray(noise_x.sigma_ch, dtype=np.float64),
        }),
        "CHANNEL_NOISE_Y": Table({
            **common,
            "SIGMA": np.asarray(noise_y.sigma_ch, dtype=np.float64),
        }),
    }


def plait_fft_cubes(
    x_bundle: OTFBundle,
    y_bundle: OTFBundle,
    *,
    linefree_velocity_windows_kms,
    noise_mode: str = "empirical_rms",
    plait_noise_mode: str = "family_scalar",
    pad_frac: float = 0.25,
    apodize: bool = True,
    apodize_alpha: float = 0.1,
    support_taper: bool = True,
    support_taper_width_pix: int | None = None,
    science_mask_mode: str = "and",
    fft_workers: int | None = None,
    dtype: str | np.dtype = np.float64,
    diagnostics: bool = True,
    min_plait_size_pix: int = 32,
    small_map_policy: str = "fallback_average",
    quality_gate_mode: str = "none",
    min_improvement_frac: float = 0.0,
) -> OTFBundle:
    if str(noise_mode).strip().lower() != "empirical_rms":
        raise ValueError("FFT/PLAIT supports only noise_mode='empirical_rms'.")
    plait_noise_mode_norm = str(plait_noise_mode).strip().lower()
    if plait_noise_mode_norm not in {"family_scalar", "family_channel"}:
        raise ValueError("plait_noise_mode must be 'family_scalar' or 'family_channel'.")
    validate_otf_bundle(x_bundle, require_variance=False)
    validate_otf_bundle(y_bundle, require_variance=False)
    if x_bundle.data.shape != y_bundle.data.shape:
        raise ValueError(f"Shape mismatch: X={x_bundle.data.shape}, Y={y_bundle.data.shape}")
    if _wcs_key(x_bundle.header) != _wcs_key(y_bundle.header):
        raise ValueError("X and Y bundles must share the same WCS and spectral axis definition.")
    if x_bundle.unit != y_bundle.unit:
        raise ValueError(f"Unit mismatch: X={x_bundle.unit!r}, Y={y_bundle.unit!r}")

    nchan, ny, nx = x_bundle.data.shape
    vel = _velocity_axis_from_header(x_bundle.header, nchan)
    linefree_mask = _make_linefree_mask_from_velocity_windows(vel, linefree_velocity_windows_kms)

    noise_x = _estimate_bundle_noise(x_bundle, linefree_mask_1d=linefree_mask, plait_noise_mode=plait_noise_mode_norm)
    noise_y = _estimate_bundle_noise(y_bundle, linefree_mask_1d=linefree_mask, plait_noise_mode=plait_noise_mode_norm)

    x_valid = _normalize_valid_mask(x_bundle)
    y_valid = _normalize_valid_mask(y_bundle)
    x_support = _normalize_support_mask(x_bundle)
    y_support = _normalize_support_mask(y_bundle)

    if science_mask_mode.lower() == "and":
        valid_out = x_valid & y_valid
        support_out = x_support & y_support
    elif science_mask_mode.lower() == "union":
        valid_out = x_valid | y_valid
        support_out = x_support | y_support
    else:
        raise ValueError("science_mask_mode must be 'and' or 'union'.")

    x_data_raw = np.asarray(x_bundle.data, dtype=float)
    y_data_raw = np.asarray(y_bundle.data, dtype=float)

    # Average comparator / small-map fallback always available.
    if plait_noise_mode_norm == "family_channel":
        avg_out, var_avg = _family_weighted_average_output_family_channel(
            x_data_raw, y_data_raw, x_valid, y_valid,
            sigma_x_ch=noise_x.sigma_ch, sigma_y_ch=noise_y.sigma_ch,
        )
    else:
        avg_out, var_avg = _family_weighted_average_output(
            x_data_raw, y_data_raw, x_valid, y_valid,
            sigma_x=noise_x.sigma_rep, sigma_y=noise_y.sigma_rep,
        )
    avg_out[~valid_out] = np.nan
    var_avg[~valid_out] = np.nan

    if builtins.min(ny, nx) < int(min_plait_size_pix):
        if str(small_map_policy).lower() == "error":
            raise ValueError(f"FFT/PLAIT requires min(ny,nx) >= {min_plait_size_pix}; got {(ny, nx)}")
        header = x_bundle.header.copy()
        header["BWTYPE"] = ("AVG_FALLBACK", "small-map fallback instead of FFT/PLAIT")
        header["PLWGMODE"] = (plait_noise_mode_norm, "PLAIT Fourier weight mode")
        image_ext = {
            "RMS_MAP_X": noise_x.rms_map.astype(np.float32),
            "RMS_MAP_Y": noise_y.rms_map.astype(np.float32),
            "VALID_MASK_AND": (x_valid & y_valid).astype(np.uint8),
            "VALID_MASK_UNION": (x_valid | y_valid).astype(np.uint8),
        }
        image_ext.update(_collect_input_family_diagnostics(x_bundle, y_bundle))
        table_ext = {}
        if diagnostics:
            table_ext["CHANNEL_NOISE"] = Table({
                "PLWGMODE": [plait_noise_mode_norm],
                "SIGMA_X_REP": [float(noise_x.sigma_rep)],
                "SIGMA_Y_REP": [float(noise_y.sigma_rep)],
                "LINEFREE_NCHAN": [int(np.count_nonzero(linefree_mask))],
                "SELECTED": ["fallback_average_small_map"],
                "SIGMA_OUT_PLAIT": [np.nan],
                "SIGMA_OUT_AVG": [float(_estimate_output_sigma_rep(avg_out, valid_out, support_out, linefree_mask))],
                "NOISE_PIX_X": [int(noise_x.noise_mask_npix)],
                "NOISE_PIX_Y": [int(noise_y.noise_mask_npix)],
            })
            table_ext.update(_make_channel_noise_tables(velocities_kms=vel, linefree_mask=linefree_mask, noise_x=noise_x, noise_y=noise_y))
        out_bundle = OTFBundle(
            data=avg_out,
            header=header,
            variance=var_avg,
            valid_mask=valid_out,
            support_mask=support_out,
            unit=x_bundle.unit,
            family_label="PLAIT",
            image_ext=image_ext,
            table_ext=table_ext,
            meta={
                "basketweave_method": "fallback_average_small_map",
                "noise_mode": "empirical_rms",
                "plait_noise_mode": plait_noise_mode_norm,
                "linefree_velocity_windows_kms": list(_resolve_velocity_windows_spec(linefree_velocity_windows_kms) or []),
                "baseline_subtracted": bool(x_bundle.meta.get("baseline_subtracted", False) or y_bundle.meta.get("baseline_subtracted", False)),
                "sigma_x_rep": float(noise_x.sigma_rep),
                "sigma_y_rep": float(noise_y.sigma_rep),
                "small_map_fallback": True,
                "min_plait_size_pix": int(min_plait_size_pix),
                "input_family_diagnostic_ext": sorted(image_ext.keys()),
            },
        )
        from .mosaic import attach_mosaic_products
        attach_mosaic_products(
            out_bundle,
            linefree_velocity_windows_kms=linefree_velocity_windows_kms,
            gain_map=np.ones((ny, nx), dtype=float),
            trust_map=np.ones((ny, nx), dtype=float),
            gain_source="unity",
            trust_source="unity",
            gain_min=0.5,
            overwrite=True,
            in_place=True,
        )
        return out_bundle

    # FFT path uses taper/apodization only on data sent to FFT, not on RMS estimation.
    x_data = np.where(x_valid, np.asarray(x_bundle.data, dtype=dtype), 0.0)
    y_data = np.where(y_valid, np.asarray(y_bundle.data, dtype=dtype), 0.0)

    if support_taper:
        x_taper = _build_support_taper(x_support, width_pix=support_taper_width_pix)
        y_taper = _build_support_taper(y_support, width_pix=support_taper_width_pix)
        x_data = x_data * x_taper[None, :, :]
        y_data = y_data * y_taper[None, :, :]
    if apodize:
        x_data = _apply_apodization(x_data, alpha=apodize_alpha, width_pix=support_taper_width_pix)
        y_data = _apply_apodization(y_data, alpha=apodize_alpha, width_pix=support_taper_width_pix)

    x_pad, crop = _pad_spatial(x_data, pad_frac=pad_frac)
    y_pad, _ = _pad_spatial(y_data, pad_frac=pad_frac)

    fx = spfft.fftn(x_pad, axes=(-2, -1), workers=fft_workers)
    fy = spfft.fftn(y_pad, axes=(-2, -1), workers=fft_workers)
    if plait_noise_mode_norm == "family_channel":
        fout, wx, wy, mean_wx2_ch, mean_wy2_ch = _plait_combine_family_channel(
            fx,
            fy,
            sigma_x_ch=noise_x.sigma_ch,
            sigma_y_ch=noise_y.sigma_ch,
        )
    else:
        wx, wy = _build_fourier_weights(x_pad.shape[1], x_pad.shape[2], sigma_x=noise_x.sigma_rep, sigma_y=noise_y.sigma_rep)
        fout = fx * wx[None, :, :] + fy * wy[None, :, :]
        mean_wx2_ch = mean_wy2_ch = None
    out_pad = spfft.ifftn(fout, axes=(-2, -1), workers=fft_workers).real
    out = np.asarray(out_pad[:, crop[0], crop[1]], dtype=float)
    out[~valid_out] = np.nan
    if plait_noise_mode_norm == "family_channel":
        var_out = _approx_output_variance_cube_family_channel(
            out.shape,
            sigma_x_ch=noise_x.sigma_ch,
            sigma_y_ch=noise_y.sigma_ch,
            mean_wx2_ch=mean_wx2_ch,
            mean_wy2_ch=mean_wy2_ch,
        )
    else:
        var_out = _approx_output_variance_cube(out.shape, sigma_x=noise_x.sigma_rep, sigma_y=noise_y.sigma_rep, wx=wx, wy=wy)
    var_out[~valid_out] = np.nan

    # Output empirical RMS diagnostics are always computed on the final candidate cubes.
    sigma_plait = _estimate_output_sigma_rep(out, valid_out, support_out, linefree_mask)
    sigma_avg = _estimate_output_sigma_rep(avg_out, valid_out, support_out, linefree_mask)
    selection = PlaitSelectionDiagnostics(sigma_plait=float(sigma_plait), sigma_average=float(sigma_avg), selected_method="plait_fft")
    if str(quality_gate_mode).strip().lower() == "lower_empirical_rms":
        if np.isfinite(sigma_avg) and np.isfinite(sigma_plait):
            if sigma_plait >= sigma_avg * (1.0 - float(min_improvement_frac)):
                out = avg_out
                var_out = var_avg
                selection.selected_method = "average_quality_gate"

    header = x_bundle.header.copy()
    header["BWTYPE"] = ("PLAITFFT", "basketweave method")
    header["PLWGMODE"] = (plait_noise_mode_norm, "PLAIT Fourier weight mode")
    image_ext = {
        "RMS_MAP_X": noise_x.rms_map.astype(np.float32),
        "RMS_MAP_Y": noise_y.rms_map.astype(np.float32),
        "VALID_MASK_AND": (x_valid & y_valid).astype(np.uint8),
        "VALID_MASK_UNION": (x_valid | y_valid).astype(np.uint8),
        "WEIGHT_X_FFT": wx.astype(np.float32),
        "WEIGHT_Y_FFT": wy.astype(np.float32),
    }
    image_ext.update(_collect_input_family_diagnostics(x_bundle, y_bundle))
    table_ext = {}
    if diagnostics:
        table_ext["CHANNEL_NOISE"] = Table({
            "PLWGMODE": [plait_noise_mode_norm],
            "SIGMA_X_REP": [float(noise_x.sigma_rep)],
            "SIGMA_Y_REP": [float(noise_y.sigma_rep)],
            "LINEFREE_NCHAN": [int(np.count_nonzero(linefree_mask))],
            "SIGMA_OUT_PLAIT": [float(selection.sigma_plait)],
            "SIGMA_OUT_AVG": [float(selection.sigma_average)],
            "SELECTED": [selection.selected_method],
            "NOISE_PIX_X": [int(noise_x.noise_mask_npix)],
            "NOISE_PIX_Y": [int(noise_y.noise_mask_npix)],
        })
        table_ext.update(_make_channel_noise_tables(velocities_kms=vel, linefree_mask=linefree_mask, noise_x=noise_x, noise_y=noise_y))

    out_bundle = OTFBundle(
        data=out,
        header=header,
        variance=var_out,
        valid_mask=valid_out,
        support_mask=support_out,
        unit=x_bundle.unit,
        family_label="PLAIT",
        image_ext=image_ext,
        table_ext=table_ext,
        meta={
            "basketweave_method": "plait_fft",
            "noise_mode": "empirical_rms",
            "plait_noise_mode": plait_noise_mode_norm,
            "linefree_velocity_windows_kms": list(_resolve_velocity_windows_spec(linefree_velocity_windows_kms) or []),
            "baseline_subtracted": bool(x_bundle.meta.get("baseline_subtracted", False) or y_bundle.meta.get("baseline_subtracted", False)),
            "sigma_x_rep": float(noise_x.sigma_rep),
            "sigma_y_rep": float(noise_y.sigma_rep),
            "selection_method": selection.selected_method,
            "input_family_diagnostic_ext": sorted(image_ext.keys()),
        },
    )
    bw_gain_map = _estimate_bw_flatfield_gain_approx(
        x_support,
        y_support,
        wx=wx,
        wy=wy,
        crop=crop,
        pad_frac=pad_frac,
        apodize=bool(apodize),
        apodize_alpha=apodize_alpha,
        support_taper=bool(support_taper),
        support_taper_width_pix=support_taper_width_pix,
        fft_workers=fft_workers,
        support_out_mask_2d=support_out,
    )
    gain_map_for_output = bw_gain_map if selection.selected_method == "plait_fft" else np.ones((ny, nx), dtype=float)

    from .mosaic import attach_mosaic_products
    attach_mosaic_products(
        out_bundle,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        gain_map=gain_map_for_output,
        trust_map=np.ones((ny, nx), dtype=float),
        gain_source=("bw_flatfield_response" if selection.selected_method == "plait_fft" else "unity"),
        trust_source="unity",
        gain_min=0.5,
        overwrite=True,
        in_place=True,
    )
    return out_bundle
