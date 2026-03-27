from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.table import Table

from .config import GridResult
from .wcs_proj import build_spatial_wcs_dict
from .fits_io import _fill_header_metadata
from .otf_bundle import OTFBundle

_MASK_NAMES = {"VALID_MASK", "SUPPORT_MASK", "MASK", "VALID_MASK_AND", "VALID_MASK_UNION", "LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR", "LINEFREE3D", "LINEFREE3D_USED", "LINEFREE3D_PRIOR", "SIGNAL_MASK_USED", "SIGNAL_MASK3D_USED"}
_IMAGE_CORE_NAMES = {"VARIANCE", "VALID_MASK", "SUPPORT_MASK"}
_U8_IMAGE_NAMES = {"BASE_FLG"}


def _bool_to_u8(arr: np.ndarray) -> np.ndarray:
    return np.asarray(arr, dtype=np.uint8)


def _u8_to_bool_if_mask(name: str, arr: np.ndarray) -> np.ndarray:
    if str(name).upper() in _MASK_NAMES:
        return np.asarray(arr) != 0
    return np.asarray(arr)


def make_otf_header(
    grid_res: GridResult,
    *,
    v_tgt: np.ndarray,
    coord_sys: str,
    projection: str,
    lon0: float,
    lat0: float,
    config,
    out_scale: str = "TA*",
    rep_beameff: float = 1.0,
    family_label: str | None = None,
) -> fits.Header:
    ny, nx, _nchan = grid_res.cube.shape
    wcs_hdr = build_spatial_wcs_dict(
        coord_sys, projection, lon0, lat0, config.x0, config.y0, config.cell_arcsec, nx, ny
    )
    header = fits.Header()
    header.update(wcs_hdr)
    _fill_header_metadata(header, coord_sys, lon0, lat0, out_scale, rep_beameff, grid_res)
    header["CTYPE3"] = "VRAD"
    header["CUNIT3"] = "km/s"
    header["CRPIX3"] = 1.0
    header["CRVAL3"] = float(v_tgt[0]) if len(v_tgt) else 0.0
    header["CDELT3"] = float(v_tgt[1] - v_tgt[0]) if len(v_tgt) > 1 else 1.0
    if family_label is not None:
        header["FAMILY"] = (str(family_label), "OTF family label")
    if getattr(grid_res, "meta", None):
        meta = grid_res.meta
        if meta.get("RESTFREQ", 0) > 0:
            header["RESTFRQ"] = (meta["RESTFREQ"], "Rest Frequency (Hz)")
        if "SPECSYS" in meta:
            header["SPECSYS"] = (str(meta["SPECSYS"]), "Spectral reference frame")
        header["VELDEF"] = (str(meta.get("VELDEF", "RADI-LSR")), "Velocity definition")
    return header


def gridresult_to_otf_bundle(
    grid_res: GridResult,
    *,
    v_tgt: np.ndarray,
    coord_sys: str,
    projection: str,
    lon0: float,
    lat0: float,
    config,
    out_scale: str = "TA*",
    rep_beameff: float = 1.0,
    family_label: str | None = None,
    extra_meta: dict[str, Any] | None = None,
) -> OTFBundle:
    """Convert current GridResult into the simplified v1 OTFBundle."""
    header = make_otf_header(
        grid_res,
        v_tgt=v_tgt,
        coord_sys=coord_sys,
        projection=projection,
        lon0=lon0,
        lat0=lat0,
        config=config,
        out_scale=out_scale,
        rep_beameff=rep_beameff,
        family_label=family_label,
    )
    data = np.transpose(np.asarray(grid_res.cube, dtype=float), (2, 0, 1))
    support_mask = np.asarray(grid_res.mask_map, dtype=bool)
    valid_mask = np.isfinite(data) & support_mask[None, :, :]

    variance = None
    rms_map = getattr(grid_res, "rms_map", None)
    if rms_map is not None:
        rms2d = np.asarray(rms_map, dtype=float)
        var2d = np.full(rms2d.shape, np.nan, dtype=float)
        good = np.isfinite(rms2d) & (rms2d > 0)
        var2d[good] = rms2d[good] * rms2d[good]
        variance = np.broadcast_to(var2d[None, :, :], data.shape).copy()
        variance[~valid_mask] = np.nan
    else:
        weight_map = getattr(grid_res, "weight_map", None)
        if weight_map is not None:
            w = np.asarray(weight_map, dtype=float)
            var2d = np.full(w.shape, np.nan, dtype=float)
            good = np.isfinite(w) & (w > 0)
            var2d[good] = 1.0 / w[good]
            variance = np.broadcast_to(var2d[None, :, :], data.shape).copy()
            variance[~valid_mask] = np.nan

    image_ext: dict[str, np.ndarray] = {}
    mapping = {
        "WEIGHT": getattr(grid_res, "weight_map", None),
        "WEIGHT_SUM": getattr(grid_res, "weight_map", None),
        "HIT": getattr(grid_res, "hit_map", None),
        "NSAMP": getattr(grid_res, "nsamp_map", None),
        "WSUM": getattr(grid_res, "wsum_map", None),
        "WABS": getattr(grid_res, "wabs_map", None),
        "CANCEL": getattr(grid_res, "cancel_map", None),
        "WREL": getattr(grid_res, "weight_rel_map", None),
        "MASK": np.asarray(grid_res.mask_map, dtype=bool),
        "TSYS": getattr(grid_res, "tsys_map", None),
        "TINT": getattr(grid_res, "tint_map", None),
        "TIME": getattr(grid_res, "time_map", None),
        "RMS": getattr(grid_res, "rms_map", None),
        "NEFF": getattr(grid_res, "neff_map", None),
        "XEFF": getattr(grid_res, "xeff_map", None),
        "YEFF": getattr(grid_res, "yeff_map", None),
        "BIAS_PIX": getattr(grid_res, "dr_eff_map_pix", None),
    }
    for name, arr in mapping.items():
        if arr is not None:
            image_ext[name] = np.asarray(arr)

    meta = dict(getattr(grid_res, "meta", {}) or {})
    if extra_meta:
        meta.update(extra_meta)
    if family_label is not None:
        meta["family_label"] = str(family_label)
    if variance is not None:
        meta.setdefault("variance_source", "rms_map^2" if rms_map is not None else "1/weight_map (fallback)")

    return OTFBundle(
        data=data,
        header=header,
        variance=variance,
        valid_mask=valid_mask,
        support_mask=support_mask,
        unit=str(header.get("BUNIT", "")) or None,
        family_label=family_label,
        image_ext=image_ext,
        table_ext={},
        meta=meta,
    )


def validate_otf_bundle(bundle: OTFBundle, *, require_variance: bool = False) -> None:
    if not isinstance(bundle, OTFBundle):
        raise TypeError(f"Expected OTFBundle; got {type(bundle)!r}")
    if bundle.data.ndim != 3:
        raise ValueError(f"bundle.data must be 3D; got shape={bundle.data.shape}")
    nchan, ny, nx = bundle.data.shape
    if require_variance and bundle.variance is None:
        raise ValueError("This operation requires bundle.variance, but it is missing.")
    if bundle.variance is not None:
        arr = np.asarray(bundle.variance)
        if arr.ndim not in {1, 3}:
            raise ValueError(f"bundle.variance must be 1D or 3D; got shape={arr.shape}")
        if arr.ndim == 1 and arr.shape != (nchan,):
            raise ValueError(f"1D variance must have length nchan={nchan}; got shape={arr.shape}")
        if arr.ndim == 3 and arr.shape not in {(nchan, ny, nx), (nchan, 1, 1)}:
            raise ValueError(f"3D variance must be (nchan,ny,nx) or (nchan,1,1); got {arr.shape}")
    if bundle.support_mask is not None and np.asarray(bundle.support_mask).shape != (ny, nx):
        raise ValueError(f"support_mask must have shape (ny,nx)=({ny},{nx}); got {np.asarray(bundle.support_mask).shape}")
    if bundle.valid_mask is not None:
        vm = np.asarray(bundle.valid_mask)
        if vm.shape not in {(ny, nx), (nchan, ny, nx)}:
            raise ValueError(f"valid_mask must have shape (ny,nx) or (nchan,ny,nx); got {vm.shape}")


def _table_to_bintable_hdu(table: Table, name: str):
    try:
        return fits.BinTableHDU(table, name=name)
    except Exception:
        pass
    if hasattr(table, "as_array"):
        try:
            return fits.BinTableHDU(table.as_array(), name=name)
        except Exception:
            pass
    try:
        return fits.BinTableHDU(np.asarray(table), name=name)
    except Exception as exc:
        raise TypeError(f"Could not convert table extension {name!r} into BinTableHDU") from exc


def write_otf_bundle(bundle: OTFBundle, path: str, *, overwrite: bool = False) -> None:
    validate_otf_bundle(bundle, require_variance=False)
    header = bundle.header.copy()
    if bundle.unit is not None and "BUNIT" not in header:
        header["BUNIT"] = str(bundle.unit)
    if bundle.family_label is not None:
        header["FAMILY"] = str(bundle.family_label)
    primary = fits.PrimaryHDU(data=np.asarray(bundle.data, dtype=np.float32), header=header)
    hdus = [primary]
    if bundle.variance is not None:
        hdus.append(fits.ImageHDU(data=np.asarray(bundle.variance, dtype=np.float32), name="VARIANCE"))
    if bundle.support_mask is not None:
        hdus.append(fits.ImageHDU(data=_bool_to_u8(bundle.support_mask), name="SUPPORT_MASK"))
    if bundle.valid_mask is not None:
        hdus.append(fits.ImageHDU(data=_bool_to_u8(bundle.valid_mask), name="VALID_MASK"))

    reserved = set(_IMAGE_CORE_NAMES)
    for name, arr in bundle.image_ext.items():
        uname = str(name).upper()
        if uname in reserved:
            continue
        if uname in _MASK_NAMES:
            data = _bool_to_u8(arr)
        elif uname in _U8_IMAGE_NAMES:
            data = np.asarray(arr, dtype=np.uint8)
        else:
            data = np.asarray(arr, dtype=np.float32)
        hdus.append(fits.ImageHDU(data=data, name=uname))

    for name, table in bundle.table_ext.items():
        if not isinstance(table, Table):
            table = Table(table)
        hdus.append(_table_to_bintable_hdu(table, str(name).upper()))

    fits.HDUList(hdus).writeto(path, overwrite=overwrite)


def read_otf_bundle(path: str) -> OTFBundle:
    with fits.open(path) as hdul:
        primary = hdul[0]
        data = np.asarray(primary.data)
        if data.ndim != 3:
            raise ValueError(f"PRIMARY HDU must be 3D; got shape={data.shape}")
        header = primary.header.copy()
        variance = None
        valid_mask = None
        support_mask = None
        image_ext: dict[str, np.ndarray] = {}
        table_ext: dict[str, Table] = {}
        for hdu in hdul[1:]:
            name = str(hdu.name).upper()
            if isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
                arr = _u8_to_bool_if_mask(name, hdu.data)
                if name == "VARIANCE":
                    variance = np.asarray(arr)
                elif name == "VALID_MASK":
                    valid_mask = np.asarray(arr, dtype=bool)
                elif name == "SUPPORT_MASK":
                    support_mask = np.asarray(arr, dtype=bool)
                else:
                    image_ext[name] = np.asarray(arr)
            elif isinstance(hdu, fits.BinTableHDU):
                table_ext[name] = Table(hdu.data)
        bundle = OTFBundle(
            data=np.asarray(data),
            header=header,
            variance=variance,
            valid_mask=valid_mask,
            support_mask=support_mask,
            unit=str(header.get("BUNIT", "")) or None,
            family_label=header.get("FAMILY"),
            image_ext=image_ext,
            table_ext=table_ext,
            meta={"source_path": str(Path(path))},
        )
    return bundle
