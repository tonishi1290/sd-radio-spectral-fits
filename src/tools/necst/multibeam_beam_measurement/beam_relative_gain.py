from __future__ import annotations

import builtins
import warnings

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from astropy.table import Table

from .otf_bundle import OTFBundle

__all__ = ["estimate_relative_beam_intensities"]


# =========================
# Low-level helpers
# =========================

def _header_get_case_insensitive(header: Any, key: str, default: Any = None) -> Any:
    if header is None:
        return default
    try:
        if key in header:
            return header[key]
    except Exception:
        pass
    ukey = str(key).upper()
    try:
        for k in header.keys():
            if str(k).upper() == ukey:
                return header[k]
    except Exception:
        pass
    return default


def _meta_get_case_insensitive(meta: Any, key: str, default: Any = None) -> Any:
    if not isinstance(meta, Mapping):
        return default
    if key in meta:
        return meta[key]
    ukey = str(key).upper()
    lkey = str(key).lower()
    for k, v in meta.items():
        ks = str(k)
        if ks.upper() == ukey or ks.lower() == lkey:
            return v
    return default


def _canonical_axis_unit(unit: Any) -> str:
    u = str(unit or "").strip().lower().replace(" ", "")
    if u in {"deg", "degree", "degrees"}:
        return "deg"
    if u in {"arcmin", "arcminute", "arcminutes", "amin"}:
        return "arcmin"
    if u in {"arcsec", "arcsecond", "arcseconds", "asec", '"'}:
        return "arcsec"
    if u in {"rad", "radian", "radians"}:
        return "rad"
    if u in {"km/s", "kms-1", "kmsec-1", "km.s-1", "kms^-1", "km.s^-1"}:
        return "km/s"
    if u in {"m/s", "ms-1", "msec-1", "m.s-1", "ms^-1", "m.s^-1"}:
        return "m/s"
    return u


def _unit_to_arcsec_scale(unit: Any) -> float:
    u = _canonical_axis_unit(unit or "deg")
    if u == "deg":
        return 3600.0
    if u == "arcmin":
        return 60.0
    if u == "arcsec":
        return 1.0
    if u == "rad":
        return 180.0 * 3600.0 / np.pi
    raise ValueError(f"Unsupported spatial unit for CUNIT1/2: {unit!r}")


def _unit_to_kms_scale(unit: Any) -> float:
    u = _canonical_axis_unit(unit or "km/s")
    if u == "km/s":
        return 1.0
    if u == "m/s":
        return 1.0e-3
    raise ValueError(f"Unsupported spectral unit for CUNIT3: {unit!r}")


def _world_value(value: Any, unit: Any, *, target: str) -> float:
    fv = float(value)
    if target == "arcsec":
        return fv * _unit_to_arcsec_scale(unit)
    if target == "km/s":
        return fv * _unit_to_kms_scale(unit)
    raise ValueError(target)


def _float_close(a: Any, b: Any, *, rtol: float = 1.0e-8, atol: float = 1.0e-10) -> bool:
    try:
        af = float(a)
        bf = float(b)
    except Exception:
        return False
    return bool(np.isclose(af, bf, rtol=rtol, atol=atol, equal_nan=True))


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


def _finite_spatial_footprint(bundle: OTFBundle, *, relative_weight_threshold: float) -> np.ndarray:
    support = _normalize_support_mask(bundle, relative_weight_threshold=relative_weight_threshold)
    valid3 = _normalize_valid_mask(bundle)
    valid_any = np.any(valid3 & np.isfinite(bundle.data), axis=0)
    return support & valid_any


def _spatial_cell_arcsec(header: Any) -> tuple[float, float]:
    cunit1 = _header_get_case_insensitive(header, "CUNIT1", "deg")
    cunit2 = _header_get_case_insensitive(header, "CUNIT2", cunit1)
    sx = abs(_world_value(_header_get_case_insensitive(header, "CDELT1", np.nan), cunit1, target="arcsec"))
    sy = abs(_world_value(_header_get_case_insensitive(header, "CDELT2", np.nan), cunit2, target="arcsec"))
    if not np.isfinite(sx) or not np.isfinite(sy) or sx <= 0.0 or sy <= 0.0:
        raise ValueError("Header is missing valid spatial cell size information.")
    return float(sx), float(sy)


def _beam_fwhm_arcsec(bundle: OTFBundle) -> float | None:
    header = bundle.header
    candidates: list[float] = []
    for key in ("BMAJ", "BMIN"):
        val = _header_get_case_insensitive(header, key, None)
        if val is None:
            continue
        try:
            fval = float(val)
        except Exception:
            continue
        if np.isfinite(fval) and fval > 0.0:
            candidates.append(fval * 3600.0)
    meta = bundle.meta if isinstance(bundle.meta, Mapping) else {}
    for key in (
        "beam_fwhm_arcsec",
        "bmaj_empirical_arcsec",
        "bmin_empirical_arcsec",
        "bmaj_nominal_arcsec",
        "bmin_nominal_arcsec",
        "bmaj_eff_arcsec",
        "bmin_eff_arcsec",
    ):
        val = _meta_get_case_insensitive(meta, key, None)
        if val is None:
            continue
        try:
            fval = float(val)
        except Exception:
            continue
        if np.isfinite(fval) and fval > 0.0:
            candidates.append(fval)
    if not candidates:
        return None
    return float(np.nanmax(np.asarray(candidates, dtype=float)))


def _binary_erode_disk(mask: np.ndarray, radius_pix: int) -> np.ndarray:
    src = np.asarray(mask, dtype=bool)
    if radius_pix <= 0:
        return src.copy()
    ny, nx = src.shape
    out = np.ones_like(src, dtype=bool)
    r2 = float(radius_pix * radius_pix)
    for dy in range(-radius_pix, radius_pix + 1):
        for dx in range(-radius_pix, radius_pix + 1):
            if float(dy * dy + dx * dx) > r2:
                continue
            shifted = np.zeros_like(src, dtype=bool)
            if dy >= 0:
                ys_src = slice(0, ny - dy)
                ys_dst = slice(dy, ny)
            else:
                ys_src = slice(-dy, ny)
                ys_dst = slice(0, ny + dy)
            if dx >= 0:
                xs_src = slice(0, nx - dx)
                xs_dst = slice(dx, nx)
            else:
                xs_src = slice(-dx, nx)
                xs_dst = slice(0, nx + dx)
            shifted[ys_dst, xs_dst] = src[ys_src, xs_src]
            out &= shifted
    return out


def _edge_trim_mask(bundle: OTFBundle, *, edge_margin_beam: float, relative_weight_threshold: float) -> tuple[np.ndarray, list[str]]:
    qf: list[str] = []
    footprint = _finite_spatial_footprint(bundle, relative_weight_threshold=relative_weight_threshold)
    if edge_margin_beam <= 0.0:
        return footprint, qf
    beam = _beam_fwhm_arcsec(bundle)
    if beam is None:
        qf.append("beam_fwhm_missing:no_edge_trim")
        return footprint, qf
    sx, sy = _spatial_cell_arcsec(bundle.header)
    pix_scale = builtins.min(sx, sy)
    radius_pix = int(np.ceil(float(edge_margin_beam) * float(beam) / pix_scale))
    return _binary_erode_disk(footprint, radius_pix), qf


def _spectral_axis_kms(bundle: OTFBundle) -> np.ndarray:
    header = bundle.header
    nchan = int(bundle.data.shape[0])
    crpix = float(_header_get_case_insensitive(header, "CRPIX3", 1.0))
    crval = float(_header_get_case_insensitive(header, "CRVAL3", 0.0))
    cdelt = float(_header_get_case_insensitive(header, "CDELT3", 1.0))
    cunit = _header_get_case_insensitive(header, "CUNIT3", "km/s")
    scale = _unit_to_kms_scale(cunit)
    idx = np.arange(nchan, dtype=float) + 1.0
    return (crval + (idx - crpix) * cdelt) * scale


def _same_spatial_grid(a: OTFBundle, b: OTFBundle) -> bool:
    ha = a.header
    hb = b.header
    if a.data.shape[1:] != b.data.shape[1:]:
        return False
    for key in ("CTYPE1", "CTYPE2"):
        if str(_header_get_case_insensitive(ha, key, "")).upper() != str(_header_get_case_insensitive(hb, key, "")).upper():
            return False
    for key in ("CRPIX1", "CRPIX2"):
        if not _float_close(_header_get_case_insensitive(ha, key, np.nan), _header_get_case_insensitive(hb, key, np.nan)):
            return False
    for key, unit_key in (("CRVAL1", "CUNIT1"), ("CRVAL2", "CUNIT2"), ("CDELT1", "CUNIT1"), ("CDELT2", "CUNIT2")):
        av = _world_value(_header_get_case_insensitive(ha, key, np.nan), _header_get_case_insensitive(ha, unit_key, "deg"), target="arcsec")
        bv = _world_value(_header_get_case_insensitive(hb, key, np.nan), _header_get_case_insensitive(hb, unit_key, "deg"), target="arcsec")
        if not _float_close(av, bv, rtol=1.0e-8, atol=1.0e-6):
            return False
    return True


def _same_spectral_grid(a: OTFBundle, b: OTFBundle) -> bool:
    ha = a.header
    hb = b.header
    if a.data.shape[0] != b.data.shape[0]:
        return False
    if str(_header_get_case_insensitive(ha, "CTYPE3", "")).upper() != str(_header_get_case_insensitive(hb, "CTYPE3", "")).upper():
        return False
    if not _float_close(_header_get_case_insensitive(ha, "CRPIX3", np.nan), _header_get_case_insensitive(hb, "CRPIX3", np.nan)):
        return False
    for key in ("CRVAL3", "CDELT3"):
        av = _world_value(_header_get_case_insensitive(ha, key, np.nan), _header_get_case_insensitive(ha, "CUNIT3", "km/s"), target="km/s")
        bv = _world_value(_header_get_case_insensitive(hb, key, np.nan), _header_get_case_insensitive(hb, "CUNIT3", "km/s"), target="km/s")
        if not _float_close(av, bv, rtol=1.0e-8, atol=1.0e-8):
            return False
    return True


def _integrate_masked_moment0(data_cube: np.ndarray, mask: np.ndarray, velocity_axis_kms: np.ndarray) -> np.ndarray:
    data = np.asarray(data_cube, dtype=float)
    m = np.asarray(mask, dtype=bool)
    if data.ndim != 3:
        raise ValueError(f"data_cube must be 3D, got shape={data.shape}")
    if m.ndim == 1:
        if m.shape[0] != data.shape[0]:
            raise ValueError("1D mask length does not match nchan")
        m3 = np.broadcast_to(m[:, None, None], data.shape)
    elif m.ndim == 3:
        if m.shape != data.shape:
            raise ValueError(f"3D mask shape mismatch: {m.shape} vs {data.shape}")
        m3 = m
    else:
        raise ValueError(f"mask must be 1D or 3D, got shape={m.shape}")

    v = np.asarray(velocity_axis_kms, dtype=float)
    if v.shape != (data.shape[0],):
        raise ValueError("velocity_axis_kms shape mismatch")
    if v.size > 1:
        dv = float(np.abs(np.nanmedian(np.diff(v))))
    else:
        dv = 1.0
    if not np.isfinite(dv) or dv <= 0.0:
        raise ValueError("Could not determine positive channel width in km/s")

    valid_sel = m3 & np.isfinite(data)
    cube_sel = np.where(valid_sel, data, 0.0)
    with np.errstate(invalid="ignore"):
        mom0 = np.sum(cube_sel, axis=0, dtype=np.float64) * dv
    mom0 = np.asarray(mom0, dtype=float)
    nuse = np.count_nonzero(valid_sel, axis=0)
    mom0[nuse <= 0] = np.nan
    return mom0


def _integrated_sigma_map(bundle: OTFBundle, mask: np.ndarray, velocity_axis_kms: np.ndarray) -> np.ndarray | None:
    rms_map = bundle.image_ext.get("BASE_RMS")
    if rms_map is None:
        return None
    rms = np.asarray(rms_map, dtype=float)
    if rms.shape != bundle.data.shape[1:]:
        return None
    v = np.asarray(velocity_axis_kms, dtype=float)
    if v.size > 1:
        dv = float(np.abs(np.nanmedian(np.diff(v))))
    else:
        dv = 1.0
    if not np.isfinite(dv) or dv <= 0.0:
        return None
    m = np.asarray(mask, dtype=bool)
    data = np.asarray(bundle.data, dtype=float)
    if m.ndim == 1:
        if m.shape[0] != data.shape[0]:
            return None
        m3 = np.broadcast_to(m[:, None, None], data.shape)
    elif m.ndim == 3:
        if m.shape != data.shape:
            return None
        m3 = m
    else:
        return None
    nuse = np.count_nonzero(m3 & np.isfinite(data), axis=0).astype(float)
    sigma = np.asarray(rms, dtype=float) * dv * np.sqrt(np.maximum(nuse, 0.0))
    sigma[(nuse <= 0) | ~np.isfinite(rms)] = np.nan
    return sigma


def _resolve_identifiers(bundle: OTFBundle, index: int) -> dict[str, Any]:
    meta = bundle.meta if isinstance(bundle.meta, Mapping) else {}
    header = bundle.header
    out = {
        "fdnum": _meta_get_case_insensitive(meta, "fdnum", _header_get_case_insensitive(header, "FDNUM", index)),
        "ifnum": _meta_get_case_insensitive(meta, "ifnum", _header_get_case_insensitive(header, "IFNUM", 0)),
        "plnum": _meta_get_case_insensitive(meta, "plnum", _header_get_case_insensitive(header, "PLNUM", 0)),
        "stream_name": _meta_get_case_insensitive(meta, "stream_name", None),
    }
    if out["stream_name"] is None:
        for key in ("STREAM_NAME", "STREAM", "SPECTRAL_NAME", "spectral_name", "name"):
            val = _meta_get_case_insensitive(meta, key, None)
            if val is not None:
                out["stream_name"] = str(val)
                break
    if out["stream_name"] is None:
        for key in ("STREAM", "STREAMNM", "EXTNAME"):
            val = _header_get_case_insensitive(header, key, None)
            if val is not None:
                out["stream_name"] = str(val)
                break
    if out["stream_name"] is None:
        out["stream_name"] = f"bundle_{index}"
    return out


def _common_signal_mask_if_available(bundles: Sequence[OTFBundle]) -> np.ndarray | None:
    masks = []
    for bundle in bundles:
        if "SIGNAL_MASK3D_USED" in bundle.image_ext:
            m = np.asarray(bundle.image_ext["SIGNAL_MASK3D_USED"], dtype=bool)
            if m.shape != bundle.data.shape:
                raise ValueError(f"SIGNAL_MASK3D_USED shape mismatch: {m.shape} vs {bundle.data.shape}")
        elif "SIGNAL_MASK_USED" in bundle.image_ext:
            m1 = np.asarray(bundle.image_ext["SIGNAL_MASK_USED"], dtype=bool)
            if m1.ndim != 1 or m1.shape[0] != bundle.data.shape[0]:
                raise ValueError(f"SIGNAL_MASK_USED shape mismatch: {m1.shape} vs nchan={bundle.data.shape[0]}")
            m = np.broadcast_to(m1[:, None, None], bundle.data.shape).copy()
        else:
            return None
        masks.append(m)
    common = np.logical_and.reduce(masks)
    if not np.any(common):
        raise ValueError("Common intersection of SIGNAL_MASK_USED/SIGNAL_MASK3D_USED is empty.")
    return common


def _moment_maps_from_bundles(
    bundles: Sequence[OTFBundle],
    *,
    velocity_range_kms: tuple[float, float] | None,
) -> tuple[list[np.ndarray], list[np.ndarray | None], str, list[str]]:
    quality_flags: list[str] = []

    common_signal = _common_signal_mask_if_available(bundles)
    if common_signal is not None:
        ref = bundles[0]
        for other in bundles[1:]:
            if not _same_spectral_grid(ref, other):
                raise ValueError("All bundles must share the same spectral grid when using common SIGNAL_MASK products.")
        maps: list[np.ndarray] = []
        sigmas: list[np.ndarray | None] = []
        for bundle in bundles:
            v = _spectral_axis_kms(bundle)
            maps.append(_integrate_masked_moment0(bundle.data, common_signal, v))
            sigmas.append(_integrated_sigma_map(bundle, common_signal, v))
        quality_flags.append("moment_source:common_signal_mask")
        return maps, sigmas, "common_signal_mask", quality_flags

    if velocity_range_kms is not None:
        v0 = float(velocity_range_kms[0])
        v1 = float(velocity_range_kms[1])
        vmin = builtins.min(v0, v1)
        vmax = builtins.max(v0, v1)
        maps = []
        sigmas = []
        for i, bundle in enumerate(bundles):
            v = _spectral_axis_kms(bundle)
            mask1 = np.isfinite(v) & (v >= vmin) & (v <= vmax)
            if not np.any(mask1):
                raise ValueError(
                    f"velocity_range_kms={velocity_range_kms} selects no channels for stream={_resolve_identifiers(bundle, i)['stream_name']}"
                )
            maps.append(_integrate_masked_moment0(bundle.data, mask1, v))
            sigmas.append(_integrated_sigma_map(bundle, mask1, v))
        quality_flags.append("moment_source:velocity_range")
        return maps, sigmas, "velocity_range", quality_flags

    if all("MOMENT0" in b.image_ext for b in bundles):
        maps = []
        for bundle in bundles:
            mom0 = np.asarray(bundle.image_ext["MOMENT0"], dtype=float)
            if mom0.shape != bundle.data.shape[1:]:
                raise ValueError(f"MOMENT0 shape mismatch: {mom0.shape} vs {(bundle.data.shape[1], bundle.data.shape[2])}")
            maps.append(mom0)
        quality_flags.extend(["moment_source:existing_moment0", "weight:none_for_existing_moment0", "noncommon_signal_mask_risk"])
        return maps, [None] * len(bundles), "existing_moment0", quality_flags

    raise ValueError(
        "Could not determine comparison maps. Provide bundles with SIGNAL_MASK_USED/SIGNAL_MASK3D_USED, or specify velocity_range_kms, or provide MOMENT0 in every bundle."
    )


def _build_common_comparison_mask(
    moment_maps: Sequence[np.ndarray],
    sigma_maps: Sequence[np.ndarray | None],
    common_support: np.ndarray,
    *,
    positive_only: bool,
    min_snr: float,
    min_pixels: int,
    min_template_fraction: float,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    qf: list[str] = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        template = np.nanmedian(np.stack([np.asarray(m, dtype=float) for m in moment_maps], axis=0), axis=0)
    finite_all = np.asarray(common_support, dtype=bool)
    for m in moment_maps:
        finite_all &= np.isfinite(np.asarray(m, dtype=float))
    use = finite_all.copy()

    template_outlier_cap_factor = 5.0
    if positive_only:
        use &= np.isfinite(template) & (template > 0.0)
        tref = float(np.nanpercentile(template[use], 95.0)) if np.any(use) else np.nan
        if np.isfinite(tref) and tref > 0.0 and float(min_template_fraction) > 0.0:
            use &= template >= float(min_template_fraction) * tref
            qf.append(f"template_fraction:{float(min_template_fraction):g};template_ref:p95")
        if np.isfinite(tref) and tref > 0.0:
            cap = float(template_outlier_cap_factor) * tref
            before = int(np.count_nonzero(use))
            use &= template <= cap
            if int(np.count_nonzero(use)) < before:
                qf.append(f"template_upper_cap:{float(template_outlier_cap_factor):g}xp95")
    else:
        use &= np.isfinite(template) & (np.abs(template) > 0.0)
        at = np.abs(template)
        tref = float(np.nanpercentile(at[use], 95.0)) if np.any(use) else np.nan
        if np.isfinite(tref) and tref > 0.0 and float(min_template_fraction) > 0.0:
            use &= at >= float(min_template_fraction) * tref
            qf.append(f"template_fraction:{float(min_template_fraction):g};template_ref:p95")
        if np.isfinite(tref) and tref > 0.0:
            cap = float(template_outlier_cap_factor) * tref
            before = int(np.count_nonzero(use))
            use &= at <= cap
            if int(np.count_nonzero(use)) < before:
                qf.append(f"template_upper_cap:{float(template_outlier_cap_factor):g}xp95")

    if all(s is not None for s in sigma_maps):
        snr_stack = []
        for m, s in zip(moment_maps, sigma_maps):
            marr = np.asarray(m, dtype=float)
            sarr = np.asarray(s, dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                if positive_only:
                    snr = marr / sarr
                else:
                    snr = np.abs(marr) / sarr
            snr_stack.append(snr)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            ref_snr = np.nanmedian(np.stack(snr_stack, axis=0), axis=0)
        use &= np.isfinite(ref_snr) & (ref_snr >= float(min_snr))
        qf.append(f"snr_cut:{float(min_snr):g}")
    else:
        qf.append("snr_cut:none")

    nuse = int(np.count_nonzero(use))
    if nuse < int(min_pixels):
        raise ValueError(
            f"Comparison region is too small after common-region/SNR selection: {nuse} pixels < min_pixels={min_pixels}."
        )
    return use, template, qf


# =========================
# Sum-based strength
# =========================

def _strengths_from_sum(
    moment_maps: Sequence[np.ndarray],
    sigma_maps: Sequence[np.ndarray | None],
    compare_mask: np.ndarray,
    common_support: np.ndarray,
    *,
    min_background_pixels: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    strengths = []
    errs = []
    offsets = []
    qf: list[str] = ["estimator:sum"]
    use = np.asarray(compare_mask, dtype=bool)
    finite_all = np.asarray(common_support, dtype=bool)
    for m in moment_maps:
        finite_all &= np.isfinite(np.asarray(m, dtype=float))
    bg = finite_all & (~use)
    if int(np.count_nonzero(bg)) >= int(min_background_pixels):
        qf.append("sum_offset:background_median")
    else:
        qf.append("sum_offset:none")
    for m, s in zip(moment_maps, sigma_maps):
        marr = np.asarray(m, dtype=float)
        if int(np.count_nonzero(bg)) >= int(min_background_pixels):
            offset = float(np.nanmedian(marr[bg]))
        else:
            offset = 0.0
        offsets.append(offset)
        strengths.append(float(np.nansum(marr[use] - offset, dtype=np.float64)))
        if s is not None:
            sarr = np.asarray(s, dtype=float)
            errs.append(float(np.sqrt(np.nansum(np.square(sarr[use]), dtype=np.float64))))
        else:
            errs.append(np.nan)
    return np.asarray(strengths, dtype=float), np.asarray(errs, dtype=float), np.asarray(offsets, dtype=float), qf


# =========================
# Fit-based strength
# =========================

def _weighted_ratio_fit(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    denom = float(np.sum(w * x * x, dtype=np.float64))
    if not np.isfinite(denom) or denom <= 0.0:
        raise ValueError("Degenerate denominator in ratio fit.")
    a = float(np.sum(w * x * y, dtype=np.float64) / denom)
    cov = 1.0 / denom
    return a, float(np.sqrt(cov))


def _weighted_affine_fit(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[float, float, float, np.ndarray]:
    A = np.column_stack((x, np.ones_like(x, dtype=float)))
    sw = np.sqrt(np.asarray(w, dtype=float))
    Aw = A * sw[:, None]
    yw = y * sw
    coef, _, _, _ = np.linalg.lstsq(Aw, yw, rcond=None)
    a = float(coef[0])
    b = float(coef[1])
    yfit = a * x + b
    resid = y - yfit
    dof = builtins.max(1, x.size - 2)
    chi2 = float(np.sum(w * resid * resid, dtype=np.float64))
    sigma2 = chi2 / float(dof)
    xtwx = (A.T * w[None, :]) @ A
    cov = np.linalg.pinv(xtwx) * sigma2
    err_a = float(np.sqrt(cov[0, 0])) if np.isfinite(cov[0, 0]) and cov[0, 0] >= 0.0 else np.nan
    return a, b, err_a, resid


def _robust_fit_to_reference(
    ref_map: np.ndarray,
    tgt_map: np.ndarray,
    compare_mask: np.ndarray,
    *,
    sigma_ref: np.ndarray | None,
    sigma_tgt: np.ndarray | None,
    clip_sigma: float,
    max_clip_iter: int,
) -> tuple[float, float, float, str]:
    use = np.asarray(compare_mask, dtype=bool)
    x = np.asarray(ref_map, dtype=float)[use]
    y = np.asarray(tgt_map, dtype=float)[use]
    keep = np.isfinite(x) & np.isfinite(y)
    if sigma_ref is not None:
        sx = np.asarray(sigma_ref, dtype=float)[use]
        keep &= np.isfinite(sx) & (sx > 0.0)
    else:
        sx = None
    if sigma_tgt is not None:
        sy = np.asarray(sigma_tgt, dtype=float)[use]
        keep &= np.isfinite(sy) & (sy > 0.0)
    else:
        sy = None
    if np.count_nonzero(keep) < 3:
        raise ValueError("Not enough valid pixels for fit.")

    x = x[keep]
    y = y[keep]
    if sx is not None:
        sx = sx[keep]
    if sy is not None:
        sy = sy[keep]
    if sy is not None and sx is None:
        w = 1.0 / np.square(sy)
    elif sy is None and sx is None:
        w = np.ones_like(x, dtype=float)
    else:
        # effective-variance initialization; refined iteratively below using current slope
        a0_num = float(np.sum(x * y, dtype=np.float64))
        a0_den = float(np.sum(x * x, dtype=np.float64))
        a0 = a0_num / a0_den if np.isfinite(a0_den) and a0_den > 0.0 else 1.0
        var_eff = np.zeros_like(x, dtype=float)
        if sy is not None:
            var_eff += np.square(sy)
        if sx is not None:
            var_eff += (a0 * sx) ** 2
        bad = ~np.isfinite(var_eff) | (var_eff <= 0.0)
        var_eff[bad] = np.nan
        if np.all(~np.isfinite(var_eff)):
            w = np.ones_like(x, dtype=float)
        else:
            medv = float(np.nanmedian(var_eff[np.isfinite(var_eff)]))
            var_eff[~np.isfinite(var_eff)] = medv if np.isfinite(medv) and medv > 0.0 else 1.0
            w = 1.0 / var_eff

    span = float(np.nanpercentile(x, 95.0) - np.nanpercentile(x, 5.0)) if x.size >= 2 else 0.0
    if not np.isfinite(span) or span <= 0.0:
        a, err_a = _weighted_ratio_fit(x, y, w)
        return a, 0.0, err_a, "fit:ratio_fallback"

    current = np.ones_like(x, dtype=bool)
    a = np.nan
    b = np.nan
    err_a = np.nan
    fit_mode = "fit:affine"
    for _ in range(builtins.max(1, int(max_clip_iter))):
        xk = x[current]
        yk = y[current]
        wk = w[current]
        if xk.size < 3:
            break
        span_k = float(np.nanpercentile(xk, 95.0) - np.nanpercentile(xk, 5.0)) if xk.size >= 2 else 0.0
        if not np.isfinite(span_k) or span_k <= 0.0:
            a, err_a = _weighted_ratio_fit(xk, yk, wk)
            b = 0.0
            fit_mode = "fit:ratio_fallback"
            return a, b, err_a, fit_mode
        a, b, err_a, resid = _weighted_affine_fit(xk, yk, wk)
        if sy is not None or sx is not None:
            if sy is not None:
                syk = sy[current]
            else:
                syk = None
            if sx is not None:
                sxk = sx[current]
            else:
                sxk = None
            var_eff = np.zeros_like(xk, dtype=float)
            if syk is not None:
                var_eff += np.square(syk)
            if sxk is not None:
                var_eff += (a * sxk) ** 2
            bad = ~np.isfinite(var_eff) | (var_eff <= 0.0)
            if not np.all(bad):
                medv = float(np.nanmedian(var_eff[~bad]))
                var_eff[bad] = medv if np.isfinite(medv) and medv > 0.0 else 1.0
                wk = 1.0 / var_eff
                a, b, err_a, resid = _weighted_affine_fit(xk, yk, wk)
        med = float(np.nanmedian(resid))
        mad = float(1.4826 * np.nanmedian(np.abs(resid - med)))
        if not np.isfinite(mad) or mad <= 0.0:
            break
        new_local = np.abs(resid - med) <= float(clip_sigma) * mad
        if np.count_nonzero(new_local) < 3:
            break
        new_current = np.zeros_like(current, dtype=bool)
        idx = np.flatnonzero(current)
        new_current[idx[new_local]] = True
        if np.array_equal(new_current, current):
            break
        current = new_current

    return a, b, err_a, fit_mode


def _strengths_from_fit(
    moment_maps: Sequence[np.ndarray],
    sigma_maps: Sequence[np.ndarray | None],
    compare_mask: np.ndarray,
    sum_strengths: np.ndarray,
    sum_offsets: np.ndarray,
    *,
    positive_only: bool,
    clip_sigma: float,
    max_clip_iter: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, list[str]]:
    use = np.asarray(compare_mask, dtype=bool)
    ref_metric = np.full(len(moment_maps), np.nan, dtype=float)
    for i, m in enumerate(moment_maps):
        arr = np.asarray(m, dtype=float)
        vals = arr[use]
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            continue
        off = float(sum_offsets[i]) if np.isfinite(sum_offsets[i]) else 0.0
        centered = vals - off
        if positive_only:
            centered = centered[centered > 0.0]
            if centered.size == 0:
                continue
            ref_metric[i] = float(np.nanpercentile(centered, 90.0))
        else:
            ref_metric[i] = float(np.nanpercentile(np.abs(centered), 90.0))
    if not np.any(np.isfinite(ref_metric)):
        strength_metric = np.asarray(sum_strengths if positive_only else np.abs(sum_strengths), dtype=float)
        ref_idx = int(np.nanargmax(strength_metric))
    else:
        ref_idx = int(np.nanargmax(ref_metric))
    ref_map = np.asarray(moment_maps[ref_idx], dtype=float)
    slopes = np.full(len(moment_maps), np.nan, dtype=float)
    errs = np.full(len(moment_maps), np.nan, dtype=float)
    offsets = np.full(len(moment_maps), np.nan, dtype=float)
    qf = [f"reference_index:{ref_idx}"]
    slopes[ref_idx] = 1.0
    errs[ref_idx] = 0.0
    offsets[ref_idx] = 0.0
    mode_tags = []
    for i, (m, s) in enumerate(zip(moment_maps, sigma_maps)):
        if i == ref_idx:
            continue
        a, b, err_a, mode = _robust_fit_to_reference(
            ref_map,
            np.asarray(m, dtype=float),
            compare_mask,
            sigma_ref=sigma_maps[ref_idx],
            sigma_tgt=s,
            clip_sigma=clip_sigma,
            max_clip_iter=max_clip_iter,
        )
        slopes[i] = a
        errs[i] = err_a
        offsets[i] = b
        mode_tags.append(mode)
    qf.append("estimator:fit")
    if mode_tags:
        qf.extend(sorted(set(mode_tags)))
    return slopes, errs, offsets, ref_idx, qf


# =========================
# Public API
# =========================

def estimate_relative_beam_intensities(
    bundles: Sequence[OTFBundle],
    *,
    method: str = "sum",
    edge_margin_beam: float = 1.5,
    velocity_range_kms: tuple[float, float] | None = None,
    min_snr: float = 5.0,
    min_pixels: int = 16,
    positive_only: bool = True,
    min_template_fraction: float = 0.2,
    relative_weight_threshold: float = 0.10,
    clip_sigma: float = 4.0,
    max_clip_iter: int = 5,
) -> Table:
    """
    Estimate relative beam intensities from multiple OTFBundle maps.

    Definitions
    -----------
    Let ``M_b(y,x)`` be the comparison integrated-intensity map for beam ``b``.
    The comparison region ``Q`` is defined as the intersection of

    1. each bundle's trimmed valid spatial support,
    2. finite pixels in all ``M_b``, and
    3. when available, a common high-S/N region.

    Default estimator
    -----------------
    ``method='sum'`` computes

    ``A_b = sum_{(y,x) in Q} M_b(y,x)``

    and returns ``A_b / max_j A_j``.

    Optional fit estimator
    ----------------------
    ``method='fit'`` first chooses the strongest beam according to the sum-based
    estimator, then fits

    ``M_b = a_b * M_ref + c_b``

    on ``Q`` with iterative residual clipping.  If the reference-map dynamic
    range is too small, it falls back to ``M_b = a_b * M_ref``.

    Priority of available products
    ------------------------------
    1. common ``SIGNAL_MASK_USED`` / ``SIGNAL_MASK3D_USED`` across all bundles
    2. explicit ``velocity_range_kms``
    3. existing ``MOMENT0`` in all bundles

    No weak automatic signal inference is performed.
    """
    if not isinstance(bundles, Sequence) or len(bundles) == 0:
        raise ValueError("bundles must be a non-empty sequence of OTFBundle")
    if method not in {"sum", "fit"}:
        raise ValueError("method must be 'sum' or 'fit'")
    if int(min_pixels) < 1:
        raise ValueError("min_pixels must be >= 1")

    bundles = list(bundles)
    ref0 = bundles[0]
    for b in bundles[1:]:
        if not _same_spatial_grid(ref0, b):
            raise ValueError("All bundles must share the same spatial grid.")

    moment_maps, sigma_maps, moment_source, global_qf = _moment_maps_from_bundles(
        bundles,
        velocity_range_kms=velocity_range_kms,
    )

    trimmed_supports = []
    row_quality_flags: list[list[str]] = []
    for bundle in bundles:
        trimmed, qf = _edge_trim_mask(
            bundle,
            edge_margin_beam=float(edge_margin_beam),
            relative_weight_threshold=float(relative_weight_threshold),
        )
        trimmed_supports.append(trimmed)
        row_quality_flags.append(list(qf))
    common_support = np.logical_and.reduce(trimmed_supports)
    if not np.any(common_support):
        raise ValueError("Common trimmed support region is empty.")

    compare_mask, template_map, mask_qf = _build_common_comparison_mask(
        moment_maps,
        sigma_maps,
        common_support,
        positive_only=bool(positive_only),
        min_snr=float(min_snr),
        min_pixels=int(min_pixels),
        min_template_fraction=float(min_template_fraction),
    )
    global_qf.extend(mask_qf)

    sum_strengths, sum_errs, sum_offsets, sum_qf = _strengths_from_sum(
        moment_maps,
        sigma_maps,
        compare_mask,
        common_support,
        min_background_pixels=int(min_pixels),
    )
    global_qf.extend(sum_qf)

    if method == "sum":
        raw_strength = sum_strengths
        raw_err = sum_errs
        fit_offsets = np.full(len(bundles), np.nan, dtype=float)
        offset_correction = sum_offsets
        ref_idx = int(np.nanargmax(raw_strength))
    else:
        raw_strength, raw_err, fit_offsets, ref_idx, fit_qf = _strengths_from_fit(
            moment_maps,
            sigma_maps,
            compare_mask,
            sum_strengths,
            sum_offsets,
            positive_only=bool(positive_only),
            clip_sigma=float(clip_sigma),
            max_clip_iter=int(max_clip_iter),
        )
        global_qf.extend(fit_qf)
        offset_correction = np.full(len(bundles), np.nan, dtype=float)

    strength_metric = np.asarray(raw_strength if positive_only else np.abs(raw_strength), dtype=float)
    max_strength = float(np.nanmax(strength_metric))
    if not np.isfinite(max_strength) or max_strength == 0.0:
        raise ValueError("Could not normalize relative beam strengths because the maximum raw strength is zero or non-finite.")
    if bool(positive_only) and max_strength <= 0.0:
        raise ValueError(
            "All raw strengths are non-positive while positive_only=True. This usually means the selected comparison region or background correction is inconsistent with emission-dominated comparison."
        )
    relative = strength_metric / max_strength
    relative_err = np.full_like(raw_err, np.nan, dtype=float)
    strongest_idx = int(np.nanargmax(strength_metric))
    denom_err = float(raw_err[strongest_idx]) if np.isfinite(raw_err[strongest_idx]) else np.nan
    for i in range(len(bundles)):
        num = float(strength_metric[i])
        num_err = float(raw_err[i]) if np.isfinite(raw_err[i]) else np.nan
        if i == strongest_idx:
            relative_err[i] = 0.0
            continue
        if not (np.isfinite(num) and np.isfinite(num_err) and num > 0.0 and np.isfinite(max_strength) and max_strength > 0.0):
            relative_err[i] = np.nan
            continue
        frac2 = (num_err / num) ** 2
        if np.isfinite(denom_err) and denom_err > 0.0:
            frac2 += (denom_err / max_strength) ** 2
        relative_err[i] = float(relative[i] * np.sqrt(frac2))

    rows = []
    for i, bundle in enumerate(bundles):
        ident = _resolve_identifiers(bundle, i)
        qf = list(row_quality_flags[i])
        qf.extend(global_qf)
        rows.append(
            {
                "bundle_index": i,
                "fdnum": ident["fdnum"],
                "ifnum": ident["ifnum"],
                "plnum": ident["plnum"],
                "stream_name": ident["stream_name"],
                "method": method,
                "moment_source": moment_source,
                "raw_strength": float(raw_strength[i]),
                "raw_strength_err": float(raw_err[i]) if np.isfinite(raw_err[i]) else np.nan,
                "relative_strength": float(relative[i]),
                "relative_strength_err": float(relative_err[i]) if np.isfinite(relative_err[i]) else np.nan,
                "fit_offset": float(fit_offsets[i]) if np.isfinite(fit_offsets[i]) else np.nan,
                "offset_correction": float(offset_correction[i]) if np.isfinite(offset_correction[i]) else np.nan,
                "comparison_pixels": int(np.count_nonzero(compare_mask)),
                "support_pixels": int(np.count_nonzero(common_support)),
                "is_reference_for_fit": bool(i == ref_idx),
                "is_strongest": bool(i == strongest_idx),
                "quality_flag": ";".join(sorted(set(qf))),
            }
        )

    table = Table(rows=rows)
    # Put the strongest beam first for readability, keeping stable order among ties.
    order = np.argsort(np.asarray(table["relative_strength"], dtype=float))[::-1]
    table = table[order]
    return table
