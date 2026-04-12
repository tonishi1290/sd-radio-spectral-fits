from __future__ import annotations

import warnings
import time

import numpy as np
from astropy.io import fits
from astropy.table import Table

from .otf_bundle import OTFBundle
from .otf_bundle_io import read_otf_bundle, validate_otf_bundle, write_otf_bundle
from .provenance import append_bundle_provenance_step
from .plait_fft import (
    _estimate_rms_map_from_arrays,
    _make_linefree_mask_from_velocity_windows,
    _normalize_support_mask,
    _normalize_valid_mask,
    _resolve_velocity_windows_spec,
    _velocity_axis_from_header,
)


def _wcs_key(header: fits.Header) -> tuple:
    keys = (
        "CTYPE1", "CTYPE2", "CTYPE3",
        "CRPIX1", "CRPIX2", "CRPIX3",
        "CRVAL1", "CRVAL2", "CRVAL3",
        "CDELT1", "CDELT2", "CDELT3",
        "CUNIT1", "CUNIT2", "CUNIT3",
    )
    return tuple(header.get(k) for k in keys)


def _variance_to_rms_map(bundle: OTFBundle) -> np.ndarray | None:
    if bundle.variance is None:
        return None
    var = np.asarray(bundle.variance, dtype=float)
    if var.ndim == 1:
        var = var[:, None, None]
    valid = _normalize_valid_mask(bundle) & np.isfinite(bundle.data) & np.isfinite(var) & (var >= 0)
    var_use = np.where(valid, var, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        var2d = np.nanmedian(var_use, axis=0)
    rms2d = np.sqrt(var2d)
    count = np.sum(np.isfinite(var_use), axis=0)
    rms2d[count < 1] = np.nan
    return np.asarray(rms2d, dtype=float)


def _fallback_rms_map_from_ext(bundle: OTFBundle) -> tuple[np.ndarray | None, str | None]:
    for key in ("MOSAIC_RMS_OBS", "MOSAIC_RMS", "RMS_MAP_EMP", "RMS"):
        arr = bundle.image_ext.get(key)
        if arr is not None:
            return np.asarray(arr, dtype=float), f"image_ext:{key}"
    return None, None


def _coerce_map(arr: np.ndarray | None, shape: tuple[int, int], *, default: float, clip01: bool, name: str) -> np.ndarray:
    if arr is None:
        out = np.full(shape, float(default), dtype=float)
    else:
        out = np.asarray(arr, dtype=float)
        if out.shape != shape:
            raise ValueError(f"{name} must have shape {shape}; got {out.shape}")
        out = np.where(np.isfinite(out), out, np.nan)
    if clip01:
        out = np.clip(out, 0.0, 1.0)
    return np.asarray(out, dtype=float)


def estimate_mosaic_products(
    bundle: OTFBundle,
    *,
    linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None = None,
    gain_map: np.ndarray | None = None,
    trust_map: np.ndarray | None = None,
    gain_source: str | None = None,
    trust_source: str | None = None,
    gain_min: float = 0.5,
    reuse_existing_rms: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    validate_otf_bundle(bundle, require_variance=False)
    support = _normalize_support_mask(bundle)
    valid = _normalize_valid_mask(bundle) & np.isfinite(bundle.data)
    valid2d = np.any(valid, axis=0)
    linefree_resolved = _resolve_velocity_windows_spec(linefree_velocity_windows_kms)

    shape2d = support.shape

    if gain_map is None and "MOSAIC_GAIN" in bundle.image_ext:
        gain_map = np.asarray(bundle.image_ext["MOSAIC_GAIN"], dtype=float)
        if gain_source is None:
            gain_source = str(bundle.meta.get("mosaic_gain_source", "existing_ext"))
    if trust_map is None and "MOSAIC_TRUST" in bundle.image_ext:
        trust_map = np.asarray(bundle.image_ext["MOSAIC_TRUST"], dtype=float)
        if trust_source is None:
            trust_source = str(bundle.meta.get("mosaic_trust_source", "existing_ext"))

    gain = _coerce_map(gain_map, shape2d, default=1.0, clip01=True, name="gain_map")
    trust = _coerce_map(trust_map, shape2d, default=1.0, clip01=True, name="trust_map")
    gain_source = str(gain_source or ("unity" if gain_map is None else "provided"))
    trust_source = str(trust_source or ("unity" if trust_map is None else "provided"))

    rms_obs = None
    rms_source = None
    linefree_nchan = 0

    if (linefree_resolved is None) and reuse_existing_rms and ("MOSAIC_RMS_OBS" in bundle.image_ext):
        rms_obs = np.asarray(bundle.image_ext["MOSAIC_RMS_OBS"], dtype=float)
        rms_source = str(bundle.meta.get("mosaic_rms_source", "existing_ext:MOSAIC_RMS_OBS"))
    else:
        if linefree_resolved is not None:
            vel = _velocity_axis_from_header(bundle.header, bundle.nchan)
            linefree_mask = _make_linefree_mask_from_velocity_windows(vel, linefree_resolved)
            linefree_nchan = int(np.count_nonzero(linefree_mask))
            rms_obs = _estimate_rms_map_from_arrays(
                bundle.data,
                valid,
                support,
                linefree_mask,
            )
            rms_source = "linefree_final_cube"
        else:
            rms_obs = _variance_to_rms_map(bundle)
            if rms_obs is not None:
                rms_source = "variance_median"
            else:
                rms_obs, rms_source = _fallback_rms_map_from_ext(bundle)
                if rms_obs is None:
                    raise ValueError(
                        "Could not determine MOSAIC_RMS_OBS. Pass linefree_velocity_windows_kms, or ensure variance/RMS ext exists."
                    )

    rms_obs = np.asarray(rms_obs, dtype=float)
    good = (
        support
        & valid2d
        & np.isfinite(rms_obs)
        & (rms_obs > 0.0)
        & np.isfinite(gain)
        & (gain >= float(gain_min))
        & np.isfinite(trust)
        & (trust > 0.0)
    )
    weight = np.zeros_like(rms_obs, dtype=float)
    weight[good] = trust[good] * (gain[good] * gain[good]) / (rms_obs[good] * rms_obs[good])

    info = {
        "rms_source": str(rms_source),
        "linefree_nchan": int(linefree_nchan),
        "linefree_velocity_windows_kms": list(linefree_resolved or []),
        "gain_source": gain_source,
        "trust_source": trust_source,
        "gain_min_used": float(gain_min),
        "weight_formula": "trust * gain^2 / rms_obs^2",
        "numerator_formula": "trust * gain / rms_obs^2 * data",
    }
    return gain, rms_obs, trust, weight, info




def estimate_mosaic_products_from_mask(
    bundle: OTFBundle,
    *,
    linefree_mask_1d: np.ndarray,
    gain_map: np.ndarray | None = None,
    trust_map: np.ndarray | None = None,
    gain_source: str | None = None,
    trust_source: str | None = None,
    gain_min: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    validate_otf_bundle(bundle, require_variance=False)
    support = _normalize_support_mask(bundle)
    valid = _normalize_valid_mask(bundle) & np.isfinite(bundle.data)
    valid2d = np.any(valid, axis=0)
    shape2d = support.shape

    mask = np.asarray(linefree_mask_1d, dtype=bool).reshape(-1)
    if mask.shape != (bundle.nchan,):
        raise ValueError(f"linefree_mask_1d must have shape ({bundle.nchan},); got {mask.shape}")
    linefree_nchan = int(np.count_nonzero(mask))
    if linefree_nchan < 3:
        raise ValueError("At least 3 line-free channels are required to build MOSAIC_RMS_OBS.")

    if gain_map is None and "MOSAIC_GAIN" in bundle.image_ext:
        gain_map = np.asarray(bundle.image_ext["MOSAIC_GAIN"], dtype=float)
        if gain_source is None:
            gain_source = str(bundle.meta.get("mosaic_gain_source", "existing_ext"))
    if trust_map is None and "MOSAIC_TRUST" in bundle.image_ext:
        trust_map = np.asarray(bundle.image_ext["MOSAIC_TRUST"], dtype=float)
        if trust_source is None:
            trust_source = str(bundle.meta.get("mosaic_trust_source", "existing_ext"))

    gain = _coerce_map(gain_map, shape2d, default=1.0, clip01=True, name="gain_map")
    trust = _coerce_map(trust_map, shape2d, default=1.0, clip01=True, name="trust_map")
    gain_source = str(gain_source or ("unity" if gain_map is None else "provided"))
    trust_source = str(trust_source or ("unity" if trust_map is None else "provided"))

    rms_obs = _estimate_rms_map_from_arrays(bundle.data, valid, support, mask)
    rms_source = "linefree_mask_final_cube"

    good = (
        support
        & valid2d
        & np.isfinite(rms_obs)
        & (rms_obs > 0.0)
        & np.isfinite(gain)
        & (gain >= float(gain_min))
        & np.isfinite(trust)
        & (trust > 0.0)
    )
    weight = np.zeros_like(rms_obs, dtype=float)
    weight[good] = trust[good] * (gain[good] * gain[good]) / (rms_obs[good] * rms_obs[good])

    info = {
        "rms_source": str(rms_source),
        "linefree_nchan": int(linefree_nchan),
        "linefree_velocity_windows_kms": [],
        "gain_source": gain_source,
        "trust_source": trust_source,
        "gain_min_used": float(gain_min),
        "weight_formula": "trust * gain^2 / rms_obs^2",
        "numerator_formula": "trust * gain / rms_obs^2 * data",
    }
    return gain, rms_obs, trust, weight, info


def attach_mosaic_products_from_mask(
    bundle: OTFBundle,
    *,
    linefree_mask_1d: np.ndarray,
    gain_map: np.ndarray | None = None,
    trust_map: np.ndarray | None = None,
    gain_source: str | None = None,
    trust_source: str | None = None,
    gain_min: float = 0.5,
    in_place: bool = False,
    _record_provenance: bool = True,
) -> OTFBundle:
    _t0_prov = time.perf_counter()
    out = bundle if in_place else bundle.copy(deep=True)
    gain, rms_obs, trust, weight, info = estimate_mosaic_products_from_mask(
        out,
        linefree_mask_1d=linefree_mask_1d,
        gain_map=gain_map,
        trust_map=trust_map,
        gain_source=gain_source,
        trust_source=trust_source,
        gain_min=float(gain_min),
    )

    out.image_ext["MOSAIC_GAIN"] = np.asarray(gain, dtype=np.float32)
    out.image_ext["MOSAIC_RMS_OBS"] = np.asarray(rms_obs, dtype=np.float32)
    out.image_ext["MOSAIC_TRUST"] = np.asarray(trust, dtype=np.float32)
    out.image_ext["MOSAIC_WEIGHT"] = np.asarray(weight, dtype=np.float32)
    out.image_ext["MOSAIC_RMS"] = np.asarray(rms_obs, dtype=np.float32)
    out.table_ext["MOSAIC_INFO"] = Table({
        "RMS_SOURCE": [str(info["rms_source"])],
        "LINEFREE_NCHAN": [int(info["linefree_nchan"])],
        "WINDOWS": [str(info["linefree_velocity_windows_kms"])],
        "GAIN_SOURCE": [str(info["gain_source"])],
        "TRUST_SOURCE": [str(info["trust_source"])],
        "GAIN_MIN_USED": [float(info["gain_min_used"])],
        "WEIGHT_FORMULA": [str(info["weight_formula"])],
        "NUMERATOR_FORMULA": [str(info["numerator_formula"])],
    })
    out.meta["mosaic_rms_source"] = str(info["rms_source"])
    out.meta["mosaic_linefree_velocity_windows_kms"] = []
    out.meta["mosaic_gain_source"] = str(info["gain_source"])
    out.meta["mosaic_trust_source"] = str(info["trust_source"])
    out.meta["mosaic_gain_min_used"] = float(info["gain_min_used"])
    out.meta["mosaic_weight_formula"] = str(info["weight_formula"])
    out.meta["mosaic_numerator_formula"] = str(info["numerator_formula"])
    if _record_provenance:
        append_bundle_provenance_step(
            out,
            input_bundles=None,
            op_id="otf.mosaic.attach_mask.v1",
            module=__name__,
            function="attach_mosaic_products_from_mask",
            kind="aux",
            aux="mosaic_products",
            params_input={
                "linefree_mask_1d_nchan": int(np.asarray(linefree_mask_1d, dtype=bool).shape[0]),
                "gain_source": gain_source,
                "trust_source": trust_source,
                "gain_min": float(gain_min),
                "in_place": bool(in_place),
            },
            params_resolved=info,
            results_summary={
                "shape_2d": [int(v) for v in gain.shape],
                "linefree_nchan": int(info["linefree_nchan"]),
                "image_ext_written": ["MOSAIC_GAIN", "MOSAIC_RMS_OBS", "MOSAIC_TRUST", "MOSAIC_WEIGHT", "MOSAIC_RMS"],
                "table_ext_written": ["MOSAIC_INFO"],
            },
            duration_sec=float(time.perf_counter() - _t0_prov),
        )
    return out

def attach_mosaic_products(
    bundle: OTFBundle,
    *,
    linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None = None,
    gain_map: np.ndarray | None = None,
    trust_map: np.ndarray | None = None,
    gain_source: str | None = None,
    trust_source: str | None = None,
    gain_min: float = 0.5,
    overwrite: bool = False,
    in_place: bool = False,
    _record_provenance: bool = True,
) -> OTFBundle:
    _t0_prov = time.perf_counter()
    out = bundle if in_place else bundle.copy(deep=True)
    have_all = all(k in out.image_ext for k in ("MOSAIC_GAIN", "MOSAIC_RMS_OBS", "MOSAIC_TRUST", "MOSAIC_WEIGHT"))
    if (not overwrite) and have_all and (gain_map is None) and (trust_map is None) and (linefree_velocity_windows_kms is None):
        return out

    gain, rms_obs, trust, weight, info = estimate_mosaic_products(
        out,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        gain_map=gain_map,
        trust_map=trust_map,
        gain_source=gain_source,
        trust_source=trust_source,
        gain_min=float(gain_min),
        reuse_existing_rms=(not overwrite) and (linefree_velocity_windows_kms is None),
    )

    out.image_ext["MOSAIC_GAIN"] = np.asarray(gain, dtype=np.float32)
    out.image_ext["MOSAIC_RMS_OBS"] = np.asarray(rms_obs, dtype=np.float32)
    out.image_ext["MOSAIC_TRUST"] = np.asarray(trust, dtype=np.float32)
    out.image_ext["MOSAIC_WEIGHT"] = np.asarray(weight, dtype=np.float32)
    # Backward-compatible alias for readers that still look for MOSAIC_RMS.
    out.image_ext["MOSAIC_RMS"] = np.asarray(rms_obs, dtype=np.float32)

    out.table_ext["MOSAIC_INFO"] = Table({
        "RMS_SOURCE": [str(info["rms_source"])],
        "LINEFREE_NCHAN": [int(info["linefree_nchan"])],
        "WINDOWS": [str(info["linefree_velocity_windows_kms"])],
        "GAIN_SOURCE": [str(info["gain_source"])],
        "TRUST_SOURCE": [str(info["trust_source"])],
        "GAIN_MIN_USED": [float(info["gain_min_used"])],
        "WEIGHT_FORMULA": [str(info["weight_formula"])],
        "NUMERATOR_FORMULA": [str(info["numerator_formula"])],
    })
    out.meta["mosaic_rms_source"] = str(info["rms_source"])
    out.meta["mosaic_linefree_velocity_windows_kms"] = list(info["linefree_velocity_windows_kms"])
    out.meta["mosaic_gain_source"] = str(info["gain_source"])
    out.meta["mosaic_trust_source"] = str(info["trust_source"])
    out.meta["mosaic_gain_min_used"] = float(info["gain_min_used"])
    out.meta["mosaic_weight_formula"] = str(info["weight_formula"])
    out.meta["mosaic_numerator_formula"] = str(info["numerator_formula"])
    if _record_provenance:
        append_bundle_provenance_step(
            out,
            input_bundles=None,
            op_id="otf.mosaic.attach.v1",
            module=__name__,
            function="attach_mosaic_products",
            kind="aux",
            aux="mosaic_products",
            params_input={
                "linefree_velocity_windows_kms": list(_resolve_velocity_windows_spec(linefree_velocity_windows_kms) or []),
                "gain_source": gain_source,
                "trust_source": trust_source,
                "gain_min": float(gain_min),
                "overwrite": bool(overwrite),
                "in_place": bool(in_place),
            },
            params_resolved=info,
            results_summary={
                "shape_2d": [int(v) for v in gain.shape],
                "linefree_nchan": int(info["linefree_nchan"]),
                "image_ext_written": ["MOSAIC_GAIN", "MOSAIC_RMS_OBS", "MOSAIC_TRUST", "MOSAIC_WEIGHT", "MOSAIC_RMS"],
                "table_ext_written": ["MOSAIC_INFO"],
            },
            duration_sec=float(time.perf_counter() - _t0_prov),
        )
    return out


def mosaic_bundles(
    bundles,
    *,
    linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None = None,
    strict_shape: bool = True,
    strict_wcs: bool = True,
    strict_unit: bool = True,
    overwrite_mosaic_products: bool = False,
    gain_min: float = 0.5,
) -> OTFBundle:
    if not isinstance(bundles, (list, tuple)):
        bundles = [bundles]
    if len(bundles) == 0:
        raise ValueError("No bundles were supplied.")

    prepared: list[OTFBundle] = []
    ref = bundles[0]
    validate_otf_bundle(ref, require_variance=False)
    ref_shape = ref.data.shape
    ref_wcs = _wcs_key(ref.header)
    ref_unit = ref.unit

    for idx, b in enumerate(bundles):
        validate_otf_bundle(b, require_variance=False)
        if strict_shape and b.data.shape != ref_shape:
            raise ValueError(f"Bundle at index {idx} has shape={b.data.shape}, expected {ref_shape}")
        if strict_wcs and _wcs_key(b.header) != ref_wcs:
            raise ValueError(f"Bundle at index {idx} has a different WCS/axis definition.")
        if strict_unit and b.unit != ref_unit:
            raise ValueError(f"Bundle at index {idx} has unit={b.unit!r}, expected {ref_unit!r}")
        prepared.append(
            attach_mosaic_products(
                b,
                linefree_velocity_windows_kms=linefree_velocity_windows_kms,
                gain_min=float(gain_min),
                overwrite=True,
                in_place=False,
                _record_provenance=False,
            )
        )

    ref_p = prepared[0]
    shape = ref_p.data.shape
    numer = np.zeros(shape, dtype=float)
    denom = np.zeros(shape, dtype=float)
    valid_union = np.zeros(shape, dtype=bool)
    support_union = np.zeros(shape[1:], dtype=bool)
    contrib_count = np.zeros(shape[1:], dtype=float)
    weight_sum_2d = np.zeros(shape[1:], dtype=float)

    for b in prepared:
        support = _normalize_support_mask(b)
        valid = _normalize_valid_mask(b) & np.isfinite(b.data)
        gain = np.asarray(b.image_ext["MOSAIC_GAIN"], dtype=float)
        rms_obs = np.asarray(b.image_ext["MOSAIC_RMS_OBS"], dtype=float)
        trust = np.asarray(b.image_ext["MOSAIC_TRUST"], dtype=float)
        weight = np.asarray(b.image_ext["MOSAIC_WEIGHT"], dtype=float)
        good2d = (
            support
            & np.isfinite(gain)
            & (gain >= float(gain_min))
            & np.isfinite(rms_obs)
            & (rms_obs > 0.0)
            & np.isfinite(trust)
            & (trust > 0.0)
            & np.isfinite(weight)
            & (weight > 0.0)
        )
        numcoef = np.zeros_like(weight, dtype=float)
        numcoef[good2d] = trust[good2d] * gain[good2d] / (rms_obs[good2d] * rms_obs[good2d])
        use = valid & (good2d[None, :, :])
        numer += np.where(use, np.asarray(b.data, dtype=float) * numcoef[None, :, :], 0.0)
        denom += np.where(use, weight[None, :, :], 0.0)
        valid_union |= use
        support_union |= support
        contrib2d = np.any(use, axis=0)
        contrib_count += contrib2d.astype(float)
        weight_sum_2d += np.where(contrib2d, weight, 0.0)

    out = np.full(shape, np.nan, dtype=float)
    var_out = np.full(shape, np.nan, dtype=float)
    good = denom > 0.0
    out[good] = numer[good] / denom[good]
    var_out[good] = 1.0 / denom[good]

    image_ext = dict(ref_p.image_ext)
    image_ext["MOSAIC_INPUT_COUNT"] = np.asarray(contrib_count, dtype=np.float32)
    image_ext["MOSAIC_WEIGHT_SUM"] = np.asarray(weight_sum_2d, dtype=np.float32)
    image_ext["VALID_MASK_UNION"] = valid_union.astype(np.uint8)

    linefree_nchan_ref = 0
    if "MOSAIC_INFO" in ref_p.table_ext and len(ref_p.table_ext["MOSAIC_INFO"]) > 0:
        try:
            linefree_nchan_ref = int(ref_p.table_ext["MOSAIC_INFO"]["LINEFREE_NCHAN"][0])
        except Exception:
            linefree_nchan_ref = 0

    table_ext = dict(ref_p.table_ext)
    table_ext["MOSAIC_COMBINE_INFO"] = Table({
        "N_INPUT": [int(len(prepared))],
        "LINEFREE_NCHAN": [int(linefree_nchan_ref)],
        "GAIN_MIN_USED": [float(gain_min)],
        "WEIGHT_FORMULA": ["trust * gain^2 / rms_obs^2"],
        "NUMERATOR_FORMULA": ["trust * gain / rms_obs^2 * data"],
    })

    out_bundle = OTFBundle(
        data=out,
        header=ref_p.header.copy(),
        variance=var_out,
        valid_mask=valid_union,
        support_mask=support_union,
        unit=ref_p.unit,
        family_label="MOSAIC",
        image_ext=image_ext,
        table_ext=table_ext,
        meta=dict(ref_p.meta),
    )
    out_bundle.meta["mosaic_n_input"] = int(len(prepared))
    out_bundle.meta["mosaic_gain_min_used"] = float(gain_min)
    out_bundle.meta["mosaic_weight_formula"] = "trust * gain^2 / rms_obs^2"
    out_bundle.meta["mosaic_numerator_formula"] = "trust * gain / rms_obs^2 * data"
    out_bundle.meta["linefree_velocity_windows_kms"] = list(_resolve_velocity_windows_spec(linefree_velocity_windows_kms) or [])
    append_bundle_provenance_step(
        out_bundle,
        input_bundles=prepared,
        op_id="otf.mosaic.combine.v1",
        module=__name__,
        function="mosaic_bundles",
        kind="main",
        params_input={
            "linefree_velocity_windows_kms": list(_resolve_velocity_windows_spec(linefree_velocity_windows_kms) or []),
            "strict_shape": bool(strict_shape),
            "strict_wcs": bool(strict_wcs),
            "strict_unit": bool(strict_unit),
            "overwrite_mosaic_products": bool(overwrite_mosaic_products),
            "gain_min": float(gain_min),
        },
        params_resolved={
            "n_input": int(len(prepared)),
            "family_label": "MOSAIC",
            "gain_min_used": float(gain_min),
        },
        results_summary={
            "cube_shape": [int(v) for v in out_bundle.data.shape],
            "valid_voxels": int(np.count_nonzero(np.isfinite(out_bundle.data))),
            "support_npix": int(np.count_nonzero(out_bundle.support_mask)) if out_bundle.support_mask is not None else None,
        },
    )
    attach_mosaic_products(
        out_bundle,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        gain_map=np.ones(shape[1:], dtype=float),
        trust_map=np.ones(shape[1:], dtype=float),
        gain_source="unity",
        trust_source="unity",
        gain_min=float(gain_min),
        overwrite=True,
        in_place=True,
        _record_provenance=True,
    )
    return out_bundle


def mosaic_fits(
    paths,
    *,
    output_fits: str | None = None,
    overwrite: bool = False,
    linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None = None,
    strict_shape: bool = True,
    strict_wcs: bool = True,
    strict_unit: bool = True,
    overwrite_mosaic_products: bool = False,
    gain_min: float = 0.5,
) -> OTFBundle:
    if not isinstance(paths, (list, tuple)):
        paths = [paths]
    bundles = [read_otf_bundle(str(p)) for p in paths]
    out = mosaic_bundles(
        bundles,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        strict_shape=strict_shape,
        strict_wcs=strict_wcs,
        strict_unit=strict_unit,
        overwrite_mosaic_products=overwrite_mosaic_products,
        gain_min=float(gain_min),
    )
    if output_fits is not None:
        write_otf_bundle(out, output_fits, overwrite=overwrite)
    return out
