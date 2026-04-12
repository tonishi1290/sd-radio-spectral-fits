from __future__ import annotations

import numpy as np
from astropy.io import fits
from astropy.table import Table

from .otf_bundle import OTFBundle
from .otf_bundle_io import validate_otf_bundle
from .provenance import append_bundle_provenance_step
from .plait_fft import (
    _velocity_axis_from_header,
    _make_linefree_mask_from_velocity_windows,
    _normalize_support_mask,
    _normalize_valid_mask,
    _estimate_bundle_noise,
    _resolve_velocity_windows_spec,
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


def coadd_family_cubes(
    bundles,
    *,
    linefree_velocity_windows_kms,
    strict_shape: bool = True,
    strict_wcs: bool = True,
) -> OTFBundle:
    if not isinstance(bundles, (list, tuple)):
        bundles = [bundles]
    if len(bundles) == 0:
        raise ValueError("No bundles were supplied.")
    for b in bundles:
        validate_otf_bundle(b, require_variance=False)

    ref = bundles[0]
    shape = ref.data.shape
    unit = ref.unit
    fam = ref.family_label
    wcs_ref = _wcs_key(ref.header)
    for idx, b in enumerate(bundles[1:], start=1):
        if strict_shape and b.data.shape != shape:
            raise ValueError(f"Bundle at index {idx} has shape={b.data.shape}, expected {shape}")
        if strict_wcs and _wcs_key(b.header) != wcs_ref:
            raise ValueError(f"Bundle at index {idx} has a different WCS/axis definition.")
        if unit != b.unit:
            raise ValueError(f"Bundle at index {idx} has unit={b.unit!r}, expected {unit!r}")
        if fam != b.family_label:
            raise ValueError(f"Bundle at index {idx} has family={b.family_label!r}, expected {fam!r}")

    nchan = shape[0]
    vel = _velocity_axis_from_header(ref.header, nchan)
    linefree_mask = _make_linefree_mask_from_velocity_windows(vel, linefree_velocity_windows_kms)

    numer = np.zeros(shape, dtype=float)
    denom = np.zeros(shape, dtype=float)
    support_union = np.zeros(shape[1:], dtype=bool)
    valid_union = np.zeros(shape, dtype=bool)
    hit_count = np.zeros(shape[1:], dtype=float)
    rms_maps = []

    for b in bundles:
        noise = _estimate_bundle_noise(b, linefree_mask_1d=linefree_mask)
        rms_map = np.asarray(noise.rms_map, dtype=float)
        weight2d = np.zeros_like(rms_map, dtype=float)
        goodw = np.isfinite(rms_map) & (rms_map > 0)
        weight2d[goodw] = 1.0 / (rms_map[goodw] * rms_map[goodw])
        valid = _normalize_valid_mask(b) & np.isfinite(b.data) & (weight2d[None, :, :] > 0)
        numer += np.where(valid, b.data * weight2d[None, :, :], 0.0)
        denom += np.where(valid, weight2d[None, :, :], 0.0)
        valid_union |= valid
        support_union |= _normalize_support_mask(b)
        hit_count += np.any(valid, axis=0).astype(float)
        rms_maps.append(rms_map)

    out = np.full(shape, np.nan, dtype=float)
    good = denom > 0
    out[good] = numer[good] / denom[good]
    var_out = np.full(shape, np.nan, dtype=float)
    var_out[good] = 1.0 / denom[good]
    out[~valid_union] = np.nan
    var_out[~valid_union] = np.nan

    image_ext = dict(ref.image_ext)
    image_ext["HIT_COUNT_COADD"] = hit_count.astype(np.float32)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        with np.errstate(invalid="ignore"):
            rms_map_emp = np.sqrt(np.nanmedian(var_out, axis=0))
    image_ext["RMS_MAP_EMP"] = rms_map_emp.astype(np.float32)
    image_ext["VALID_MASK_UNION"] = valid_union.astype(np.uint8)
    if len(rms_maps) == 1:
        image_ext["RMS_MAP_INPUT_0"] = rms_maps[0].astype(np.float32)
    else:
        for i, arr in enumerate(rms_maps):
            image_ext[f"RMS_MAP_INPUT_{i}"] = arr.astype(np.float32)

    table_ext = dict(ref.table_ext)
    table_ext["COADD_INFO"] = Table({
        "N_INPUT": [len(bundles)],
        "FAMILY": ["" if fam is None else str(fam)],
        "LINEFREE_NCHAN": [int(np.count_nonzero(linefree_mask))],
    })

    out_bundle = OTFBundle(
        data=out,
        header=ref.header.copy(),
        variance=var_out,
        valid_mask=valid_union,
        support_mask=support_union,
        unit=ref.unit,
        family_label=ref.family_label,
        image_ext=image_ext,
        table_ext=table_ext,
        meta=dict(ref.meta),
    )
    out_bundle.meta["coadd_n_input"] = int(len(bundles))
    out_bundle.meta["linefree_velocity_windows_kms"] = list(_resolve_velocity_windows_spec(linefree_velocity_windows_kms) or [])
    out_bundle.meta["baseline_subtracted"] = bool(any(bool(b.meta.get("baseline_subtracted", False)) for b in bundles))
    out_bundle.meta["noise_mode"] = "empirical_rms_per_spectrum"
    append_bundle_provenance_step(
        out_bundle,
        input_bundles=list(bundles),
        op_id="otf.coadd.family.v1",
        module=__name__,
        function="coadd_family_cubes",
        kind="main",
        params_input={
            "linefree_velocity_windows_kms": list(_resolve_velocity_windows_spec(linefree_velocity_windows_kms) or []),
            "strict_shape": bool(strict_shape),
            "strict_wcs": bool(strict_wcs),
        },
        params_resolved={
            "n_input": int(len(bundles)),
            "family_label": None if fam is None else str(fam),
            "noise_mode": "empirical_rms_per_spectrum",
        },
        results_summary={
            "cube_shape": [int(v) for v in out_bundle.data.shape],
            "valid_voxels": int(np.count_nonzero(np.isfinite(out_bundle.data))),
            "support_npix": int(np.count_nonzero(out_bundle.support_mask)) if out_bundle.support_mask is not None else None,
        },
    )
    from .mosaic import attach_mosaic_products
    attach_mosaic_products(
        out_bundle,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        overwrite=True,
        in_place=True,
        _record_provenance=True,
    )
    return out_bundle
