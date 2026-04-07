from __future__ import annotations

import builtins
import warnings

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from astropy.table import Table

from .otf_bundle import OTFBundle

__all__ = ["BeamIntensityEstimateResult", "estimate_relative_beam_intensities"]


class BeamIntensityEstimateResult:
    """Container for relative beam intensity results."""

    def __init__(self, summary: Table, detail: Table, config: dict[str, Any]):
        self.summary = summary
        self.detail = detail
        self.config = config
        self.table = summary

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        lines = self.summary.pformat(max_width=-1, max_lines=-1)
        if self.detail.colnames:
            lines.append("")
            lines.append("詳細は result.detail、設定は result.config を参照してください。")
        return "\n".join(lines)

    def __len__(self) -> int:
        return len(self.summary)

    def __getitem__(self, item):
        return self.summary[item]

    @property
    def colnames(self):
        return self.summary.colnames

    @property
    def columns(self):
        return self.summary.columns

    def pprint(self, *args, **kwargs):
        return self.summary.pprint(*args, **kwargs)

    def pprint_detail(self, *args, **kwargs):
        return self.detail.pprint(*args, **kwargs)

    def as_table(self) -> Table:
        return self.summary


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



def _normalize_valid_mask(bundle: OTFBundle, *, relative_weight_threshold: float = 0.10) -> np.ndarray:
    nchan, ny, nx = bundle.data.shape
    support = _normalize_support_mask(bundle, relative_weight_threshold=relative_weight_threshold)
    if bundle.valid_mask is None:
        return np.isfinite(bundle.data) & support[None, :, :]
    vm = np.asarray(bundle.valid_mask, dtype=bool)
    if vm.shape == (ny, nx):
        return np.broadcast_to(vm[None, :, :], (nchan, ny, nx)).copy() & support[None, :, :]
    if vm.shape == (nchan, ny, nx):
        return vm.copy() & support[None, :, :]
    raise ValueError(f"Unsupported valid_mask shape={vm.shape}")



def _finite_spatial_footprint(bundle: OTFBundle, *, relative_weight_threshold: float) -> np.ndarray:
    support = _normalize_support_mask(bundle, relative_weight_threshold=relative_weight_threshold)
    valid3 = _normalize_valid_mask(bundle, relative_weight_threshold=relative_weight_threshold)
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
    if ny <= 0 or nx <= 0:
        return np.zeros_like(src, dtype=bool)
    if radius_pix >= builtins.max(ny, nx):
        return np.zeros_like(src, dtype=bool)

    out = np.ones_like(src, dtype=bool)
    r2 = float(radius_pix * radius_pix)
    for dy in range(-radius_pix, radius_pix + 1):
        for dx in range(-radius_pix, radius_pix + 1):
            if float(dy * dy + dx * dx) > r2:
                continue

            if abs(dy) >= ny or abs(dx) >= nx:
                return np.zeros_like(src, dtype=bool)

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

            if (ys_src.stop - ys_src.start) <= 0 or (xs_src.stop - xs_src.start) <= 0:
                return np.zeros_like(src, dtype=bool)

            shifted[ys_dst, xs_dst] = src[ys_src, xs_src]
            out &= shifted
            if not np.any(out):
                return out
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



def _base_rms_cube_if_available(bundle: OTFBundle) -> np.ndarray | None:
    rms_map = bundle.image_ext.get("BASE_RMS")
    if rms_map is None:
        return None
    rms = np.asarray(rms_map, dtype=float)
    if rms.shape != bundle.data.shape[1:]:
        return None
    return np.broadcast_to(rms[None, :, :], bundle.data.shape).astype(float, copy=False)



def _resolve_method(method: str) -> tuple[str, list[str]]:
    qf: list[str] = []
    if method == "sum":
        warnings.warn(
            "method='sum' is deprecated; it now maps to 'common_voxel_sum'.",
            DeprecationWarning,
            stacklevel=3,
        )
        qf.append("deprecated_method_alias:sum->common_voxel_sum")
        return "common_voxel_sum", qf
    if method == "common_voxel_sum":
        return method, qf
    if method == "fit":
        raise ValueError(
            "method='fit' has been removed. Use method='common_voxel_sum'."
        )
    raise ValueError("method must be 'common_voxel_sum' (or deprecated alias 'sum')")



def _resolve_template_fraction(
    *,
    template_fraction: float | None,
    min_template_fraction: float,
) -> float:
    if template_fraction is not None:
        return float(template_fraction)
    return float(min_template_fraction)



def _resolve_min_voxels(*, min_voxels: int | None, min_pixels: int | None) -> int:
    if min_voxels is not None:
        return int(min_voxels)
    if min_pixels is not None:
        warnings.warn(
            "min_pixels is deprecated; use min_voxels instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        return int(min_pixels)
    return 16



def _velocity_channel_mask(bundle: OTFBundle, velocity_range_kms: tuple[float, float] | None) -> np.ndarray:
    nchan = int(bundle.data.shape[0])
    if velocity_range_kms is None:
        return np.ones(nchan, dtype=bool)
    v0 = float(velocity_range_kms[0])
    v1 = float(velocity_range_kms[1])
    vmin = builtins.min(v0, v1)
    vmax = builtins.max(v0, v1)
    vaxis = _spectral_axis_kms(bundle)
    mask = np.isfinite(vaxis) & (vaxis >= vmin) & (vaxis <= vmax)
    if not np.any(mask):
        raise ValueError(
            f"velocity_range_kms={velocity_range_kms} selects no channels for stream={bundle.meta.get('stream_name', 'unknown') if isinstance(bundle.meta, Mapping) else 'unknown'}"
        )
    return mask



def _common_voxel_candidates(
    bundles: Sequence[OTFBundle],
    *,
    common_support: np.ndarray,
    velocity_range_kms: tuple[float, float] | None,
    relative_weight_threshold: float,
) -> tuple[np.ndarray, str, list[np.ndarray], list[np.ndarray | None], list[str]]:
    qf: list[str] = []
    valids: list[np.ndarray] = []
    rms_cubes: list[np.ndarray | None] = []
    channel_masks: list[np.ndarray] = []

    common_signal = _common_signal_mask_if_available(bundles)
    if common_signal is not None:
        qf.append("candidate_mask:common_signal")
        candidate_source = "common_signal"
    elif velocity_range_kms is not None:
        qf.append("candidate_mask:velocity_range")
        candidate_source = "velocity_range"
    else:
        qf.append("candidate_mask:all_channels")
        candidate_source = "all_channels"

    for i, bundle in enumerate(bundles):
        valid = _normalize_valid_mask(bundle, relative_weight_threshold=relative_weight_threshold)
        valid &= np.isfinite(np.asarray(bundle.data, dtype=float))
        valid &= common_support[None, :, :]
        ch_mask = _velocity_channel_mask(bundle, velocity_range_kms)
        channel_masks.append(ch_mask)
        valid &= ch_mask[:, None, None]
        if common_signal is not None:
            valid &= common_signal
        valids.append(valid)
        rms_cubes.append(_base_rms_cube_if_available(bundle))

    common_valid = np.logical_and.reduce(valids)
    if not np.any(common_valid):
        raise ValueError("Common valid 3D voxel region is empty.")
    return common_valid, candidate_source, valids, rms_cubes, qf



def _template_cube_from_common_valid(
    bundles: Sequence[OTFBundle],
    common_valid: np.ndarray,
) -> np.ndarray:
    stack = []
    for bundle in bundles:
        arr = np.asarray(bundle.data, dtype=float)
        masked = np.where(common_valid, arr, np.nan)
        stack.append(masked)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        template = np.nanmedian(np.stack(stack, axis=0), axis=0)
    template = np.asarray(template, dtype=float)
    template[~common_valid] = np.nan
    return template



def _template_sigma_cube(
    rms_cubes: Sequence[np.ndarray | None],
    common_valid: np.ndarray,
) -> np.ndarray | None:
    if not all(r is not None for r in rms_cubes):
        return None
    stack = np.stack([np.asarray(r, dtype=float) for r in rms_cubes if r is not None], axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sigma = np.nanmedian(stack, axis=0)
    sigma = np.asarray(sigma, dtype=float)
    sigma[~common_valid] = np.nan
    return sigma



def _build_common_voxel_mask(
    template_cube: np.ndarray,
    template_sigma_cube: np.ndarray | None,
    common_valid: np.ndarray,
    *,
    positive_only: bool,
    min_snr: float,
    min_voxels: int,
    template_fraction: float,
) -> tuple[np.ndarray, list[str]]:
    qf: list[str] = []
    amp = np.asarray(template_cube, dtype=float) if positive_only else np.abs(np.asarray(template_cube, dtype=float))
    use = np.asarray(common_valid, dtype=bool) & np.isfinite(amp)
    if positive_only:
        use &= template_cube > 0.0
    else:
        use &= amp > 0.0

    if not np.any(use):
        raise ValueError("Template-based comparison mask is empty before thresholding.")

    amax = float(np.nanmax(amp[use]))
    if not np.isfinite(amax) or amax <= 0.0:
        raise ValueError("Template maximum is not positive/finite; cannot build comparison voxel mask.")
    if float(template_fraction) > 0.0:
        use &= amp >= float(template_fraction) * amax
        qf.append(f"template_fraction:{float(template_fraction):g};template_ref:max")

    if template_sigma_cube is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            snr = amp / np.asarray(template_sigma_cube, dtype=float)
        use &= np.isfinite(snr) & (snr >= float(min_snr))
        qf.append(f"snr_cut:{float(min_snr):g}")
    else:
        qf.append("snr_cut:none")

    nuse = int(np.count_nonzero(use))
    if nuse < int(min_voxels):
        raise ValueError(
            f"Comparison voxel set is too small after template/SNR selection: {nuse} voxels < min_voxels={int(min_voxels)}."
        )
    return use, qf



def _strength_err_from_rms_cube(rms_cube: np.ndarray | None, compare_mask: np.ndarray) -> float:
    if rms_cube is None:
        return np.nan
    vals = np.asarray(rms_cube, dtype=float)[np.asarray(compare_mask, dtype=bool)]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    return float(np.sqrt(np.sum(vals * vals, dtype=np.float64)))



def _strengths_from_common_voxel_sum(
    bundles: Sequence[OTFBundle],
    compare_mask: np.ndarray,
    rms_cubes: Sequence[np.ndarray | None],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    qf = ["estimator:common_voxel_sum", "offset_subtraction:none", "normalization:fixed_reference"]
    strengths = np.full(len(bundles), np.nan, dtype=float)
    errs = np.full(len(bundles), np.nan, dtype=float)
    use = np.asarray(compare_mask, dtype=bool)
    for i, bundle in enumerate(bundles):
        arr = np.asarray(bundle.data, dtype=float)
        vals = arr[use]
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            continue
        strengths[i] = float(np.sum(vals, dtype=np.float64))
        errs[i] = _strength_err_from_rms_cube(rms_cubes[i], use)
    return strengths, errs, qf



def _beam_fwhm_consistency_flags(
    bundles: Sequence[OTFBundle],
    *,
    warn_fractional_mismatch: float = 0.05,
) -> tuple[list[str], list[float | None]]:
    vals: list[float | None] = [_beam_fwhm_arcsec(b) for b in bundles]
    finite = np.asarray([v for v in vals if v is not None and np.isfinite(v) and v > 0.0], dtype=float)
    qf: list[str] = []
    if finite.size != len(vals):
        qf.append("beam_fwhm_check:partial_or_missing")
        return qf, vals
    vmin = float(np.nanmin(finite))
    vmax = float(np.nanmax(finite))
    if vmin <= 0.0 or not np.isfinite(vmin) or not np.isfinite(vmax):
        qf.append("beam_fwhm_check:invalid")
        return qf, vals
    frac = vmax / vmin - 1.0
    qf.append(f"beam_fwhm_ratio_maxmin:{vmax / vmin:.6g}")
    if frac > float(warn_fractional_mismatch):
        qf.append(f"beam_fwhm_mismatch:>{100.0 * float(warn_fractional_mismatch):.1f}%")
    return qf, vals



def _relative_to_reference(
    raw_strength: np.ndarray,
    raw_err: np.ndarray,
    *,
    reference_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    ref = float(raw_strength[reference_index])
    ref_err = float(raw_err[reference_index]) if np.isfinite(raw_err[reference_index]) else np.nan
    if not np.isfinite(ref) or ref == 0.0:
        raise ValueError(
            f"Reference beam raw_strength is zero or non-finite at reference_index={reference_index}."
        )
    relative = np.asarray(raw_strength, dtype=float) / ref
    relative_err = np.full_like(relative, np.nan, dtype=float)
    relative_err[reference_index] = 0.0
    for i in range(len(relative)):
        if i == reference_index:
            continue
        num = float(raw_strength[i])
        num_err = float(raw_err[i]) if np.isfinite(raw_err[i]) else np.nan
        if not (np.isfinite(num) and np.isfinite(num_err) and num != 0.0):
            continue
        frac2 = (num_err / num) ** 2
        if np.isfinite(ref_err) and ref_err > 0.0:
            frac2 += (ref_err / ref) ** 2
        relative_err[i] = float(abs(relative[i]) * np.sqrt(frac2))
    return relative, relative_err



def estimate_relative_beam_intensities(
    bundles: Sequence[OTFBundle],
    *,
    input_labels: Sequence[str] | None = None,
    method: str = "common_voxel_sum",
    reference_index: int = 0,
    edge_margin_beam: float = 1.5,
    velocity_range_kms: tuple[float, float] | None = None,
    min_snr: float = 5.0,
    min_voxels: int | None = None,
    min_pixels: int | None = None,
    positive_only: bool = True,
    template_fraction: float | None = None,
    min_template_fraction: float = 0.2,
    relative_weight_threshold: float = 0.10,
) -> BeamIntensityEstimateResult:
    """
    Estimate relative beam intensities from multiple baseline-subtracted 3D cubes.

    The comparison is performed on a single common 3D voxel set. The mask is
    built from a template cube, while the actual strength ratio is always
    computed from the original input cubes.

    Notes
    -----
    - ``method='common_voxel_sum'`` is the primary and only supported method.
    - ``method='sum'`` is accepted as a deprecated alias to the same method.
    - ``method='fit'`` has been removed.
    - The reference beam is fixed by ``reference_index``; strongest-beam
      normalization is not used.
    """
    if not isinstance(bundles, Sequence) or len(bundles) == 0:
        raise ValueError("bundles must be a non-empty sequence of OTFBundle")

    bundles = list(bundles)
    method, method_qf = _resolve_method(method)
    min_voxels_resolved = _resolve_min_voxels(min_voxels=min_voxels, min_pixels=min_pixels)
    if min_voxels_resolved < 1:
        raise ValueError("min_voxels must be >= 1")
    template_fraction_resolved = _resolve_template_fraction(
        template_fraction=template_fraction,
        min_template_fraction=min_template_fraction,
    )
    if template_fraction_resolved < 0.0:
        raise ValueError("template_fraction must be >= 0")
    if not (0 <= int(reference_index) < len(bundles)):
        raise ValueError("reference_index is out of range")

    if input_labels is not None:
        input_labels = list(input_labels)
        if len(input_labels) != len(bundles):
            raise ValueError("input_labels must have the same length as bundles")
    else:
        input_labels = [None] * len(bundles)

    ref0 = bundles[0]
    for b in bundles[1:]:
        if not _same_spatial_grid(ref0, b):
            raise ValueError("All bundles must share the same spatial grid.")
        if not _same_spectral_grid(ref0, b):
            raise ValueError("All bundles must share the same spectral grid.")

    beam_qf, beam_fwhm_arcsec = _beam_fwhm_consistency_flags(bundles)

    trimmed_supports: list[np.ndarray] = []
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

    common_valid, candidate_source, valids, rms_cubes, global_qf = _common_voxel_candidates(
        bundles,
        common_support=common_support,
        velocity_range_kms=velocity_range_kms,
        relative_weight_threshold=float(relative_weight_threshold),
    )
    global_qf.extend(method_qf)
    global_qf.extend(beam_qf)

    template_cube = _template_cube_from_common_valid(bundles, common_valid)
    template_sigma = _template_sigma_cube(rms_cubes, common_valid)
    if template_sigma is None:
        global_qf.append("base_rms_usage:none")
    else:
        global_qf.append("base_rms_usage:mask_snr_and_error")

    compare_mask, mask_qf = _build_common_voxel_mask(
        template_cube,
        template_sigma,
        common_valid,
        positive_only=bool(positive_only),
        min_snr=float(min_snr),
        min_voxels=int(min_voxels_resolved),
        template_fraction=float(template_fraction_resolved),
    )
    global_qf.extend(mask_qf)

    raw_strength, raw_err, strength_qf = _strengths_from_common_voxel_sum(
        bundles,
        compare_mask,
        rms_cubes,
    )
    global_qf.extend(strength_qf)

    if bool(positive_only) and float(raw_strength[reference_index]) <= 0.0:
        raise ValueError(
            "Reference beam raw_strength is non-positive while positive_only=True. "
            "This usually means the comparison voxel mask is inconsistent with an emission-dominated region."
        )

    relative_strength, relative_strength_err = _relative_to_reference(
        raw_strength,
        raw_err,
        reference_index=int(reference_index),
    )

    def _warnings_from_flags(flags: list[str]) -> str:
        warns = []
        for tag in flags:
            if tag.startswith("beam_fwhm_missing:"):
                warns.append(tag)
            elif tag.startswith("beam_fwhm_mismatch:"):
                warns.append(tag)
        return ";".join(sorted(set(warns)))

    summary_rows = []
    detail_rows = []
    common_valid_voxels = int(np.count_nonzero(common_valid))
    comparison_voxels = int(np.count_nonzero(compare_mask))
    support_pixels = int(np.count_nonzero(common_support))
    for i, bundle in enumerate(bundles):
        ident = _resolve_identifiers(bundle, i)
        all_flags = list(row_quality_flags[i]) + list(global_qf)
        warning_flag = _warnings_from_flags(all_flags)
        user_label = input_labels[i]
        input_label = str(user_label) if user_label is not None else f"bundle_{i}"
        base = {
            "input_index": i,
            "input_label": input_label,
            "relative_strength": float(relative_strength[i]),
            "relative_strength_err": float(relative_strength_err[i]) if np.isfinite(relative_strength_err[i]) else np.nan,
            "raw_strength": float(raw_strength[i]),
            "raw_strength_err": float(raw_err[i]) if np.isfinite(raw_err[i]) else np.nan,
            "comparison_voxels": comparison_voxels,
            "common_valid_voxels": common_valid_voxels,
            "support_pixels": support_pixels,
            "is_reference": bool(i == int(reference_index)),
        }
        summary_rows.append(base | {"warning_flag": warning_flag})
        detail_rows.append(
            base
            | {
                "fdnum": ident["fdnum"],
                "ifnum": ident["ifnum"],
                "plnum": ident["plnum"],
                "stream_name": ident["stream_name"],
                "beam_fwhm_arcsec": float(beam_fwhm_arcsec[i]) if beam_fwhm_arcsec[i] is not None and np.isfinite(beam_fwhm_arcsec[i]) else np.nan,
                "method": method,
                "candidate_source": candidate_source,
                "base_rms_available": bool(rms_cubes[i] is not None),
                "warning_flag": warning_flag,
            }
        )

    summary = Table(rows=summary_rows)
    detail = Table(rows=detail_rows)

    config = {
        "method": method,
        "input_labels": [str(lbl) if lbl is not None else f"bundle_{i}" for i, lbl in enumerate(input_labels)],
        "reference_index": int(reference_index),
        "reference_label": str(input_labels[reference_index]) if input_labels[reference_index] is not None else f"bundle_{reference_index}",
        "edge_margin_beam": float(edge_margin_beam),
        "velocity_range_kms": tuple(map(float, velocity_range_kms)) if velocity_range_kms is not None else None,
        "min_snr": float(min_snr) if template_sigma is not None else None,
        "min_voxels": int(min_voxels_resolved),
        "positive_only": bool(positive_only),
        "template_fraction": float(template_fraction_resolved),
        "relative_weight_threshold": float(relative_weight_threshold),
        "candidate_source": candidate_source,
        "beam_fwhm_arcsec": [float(v) if v is not None and np.isfinite(v) else None for v in beam_fwhm_arcsec],
        "comparison_voxels": comparison_voxels,
        "common_valid_voxels": common_valid_voxels,
        "support_pixels": support_pixels,
        "quality_flags": sorted(set(global_qf)),
    }

    return BeamIntensityEstimateResult(summary=summary, detail=detail, config=config)
