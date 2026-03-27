#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WCS-aware plotting helpers for 2D/3D radio astronomy FITS data.

Intended placement:
    src/sd_radio_spectral_fits/map_3d/plotting.py

Main features
-------------
- 2D FITS / extension / (header, data) / HDU / Projection-like input
- 3D cube -> 2D map generation by channel or velocity selection
- provisional moment from LINEFREE / LINECAND3D / BASESUP3D
- final moment from MASK3D
- WCSAxes plotting with pseudo color + contour overlays
- linear / sqrt / log / asinh / power normalization
- Gaussian smoothing in arcsec, including target HPBW support
  (smooth_fwhm_arcsec and target_hpbw_arcsec are treated as alternatives)
- native celestial WCS support for both equatorial and galactic maps
- simple RGB composition

Notes
-----
- Contours are drawn either on the same WCS/shape or via the WCSAxes
  transform when the target axes can interpret the input WCS directly.
- Image overlays can optionally be reprojected to the base WCS when the
  optional ``reproject`` package is installed.
- For provisional moments, this implementation assumes:
    * LINECAND3D : True means line/signal candidate -> integrate there
    * BASESUP3D  : True means baseline support      -> integrate complement
    * LINEFREE   : True means line-free/baseline    -> integrate complement
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import sys
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union
import warnings

import numpy as np

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.stats import gaussian_fwhm_to_sigma
import astropy.units as u
from astropy.visualization import (
    AsinhStretch,
    ImageNormalize,
    LinearStretch,
    LogStretch,
    PowerStretch,
    SqrtStretch,
)
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from spectral_cube import SpectralCube

try:
    from .otf_bundle import OTFBundle
except Exception:  # pragma: no cover - optional when used as a standalone script
    try:
        from otf_bundle import OTFBundle  # type: ignore
    except Exception:  # pragma: no cover - environment dependent
        OTFBundle = None  # type: ignore

try:
    from spectral_cube.lower_dimensional_structures import Projection
except Exception:  # pragma: no cover - environment dependent
    Projection = None

SourceLike = Any


__all__ = [
    "Map2D",
    "RGBMap",
    "CubeInput",
    "resolve_map_input",
    "resolve_cube_input",
    "build_normalize",
    "estimate_rms_robust",
    "compute_contour_levels",
    "apply_gaussian_smoothing_2d",
    "crop_map2d",
    "describe_map2d",
    "describe_source",
    "plotting_help_text",
    "print_plotting_help",
    "make_2d_map",
    "make_provisional_moment",
    "make_final_moment",
    "add_contours",
    "plot_map",
    "make_rgb_map",
    "plot_rgb",
    "add_scalebar",
    "add_north_arrow",
    "plot_scene",
    "quicklook",
    "main",
]


@dataclass
class Map2D:
    """Container for a 2D map and its metadata."""

    data: np.ndarray
    header: fits.Header
    wcs: WCS
    unit: Optional[u.UnitBase] = None
    meta: Dict[str, Any] = field(default_factory=dict)

    @property
    def shape(self) -> Tuple[int, int]:
        return tuple(self.data.shape)


@dataclass
class RGBMap:
    """Container for RGB image data on a 2D celestial WCS."""

    rgb: np.ndarray
    header: fits.Header
    wcs: WCS
    meta: Dict[str, Any] = field(default_factory=dict)


@dataclass
class CubeInput:
    """Container for a cube-like input."""

    cube: SpectralCube
    header: fits.Header
    unit: Optional[u.UnitBase]
    meta: Dict[str, Any] = field(default_factory=dict)


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------


def _is_pathlike(obj: Any) -> bool:
    return isinstance(obj, (str, Path))



def _as_unit(unit_like: Any) -> Optional[u.UnitBase]:
    if unit_like is None:
        return None
    if isinstance(unit_like, u.UnitBase):
        return unit_like
    if hasattr(unit_like, "unit") and isinstance(getattr(unit_like, "unit"), u.UnitBase):
        return getattr(unit_like, "unit")
    try:
        return u.Unit(str(unit_like))
    except Exception:
        return None



def _copy_header(header: Optional[fits.Header]) -> fits.Header:
    return fits.Header() if header is None else fits.Header(header)



def _get_bunit_from_header(header: Optional[fits.Header]) -> Optional[u.UnitBase]:
    if header is None:
        return None
    value = header.get("BUNIT")
    if value in (None, ""):
        return None
    return _as_unit(value)



def _drop_nonspatial_wcs_keywords(header: fits.Header) -> fits.Header:
    """Create a 2D celestial header from an arbitrary FITS header."""
    hdr = _copy_header(header)
    wcs = WCS(hdr)
    chdr = wcs.celestial.to_header()
    for key in (
        "BUNIT",
        "BMAJ",
        "BMIN",
        "BPA",
        "OBJECT",
        "TELESCOP",
        "INSTRUME",
        "EQUINOX",
        "RADESYS",
        "LONPOLE",
        "LATPOLE",
    ):
        if key in hdr and key not in chdr:
            chdr[key] = hdr[key]
    return chdr



def _get_celestial_ctype_list(wcs: WCS) -> List[str]:
    """Return the first two celestial CTYPE strings robustly.

    Some Astropy/WCSLIB builds expose ``wcs.ctype`` as an object that accepts
    integer indexing but not slicing. Therefore, do not use ``[:2]`` here.
    """
    cwcs = wcs.celestial.wcs
    out: List[str] = []

    try:
        naxis = int(getattr(cwcs, "naxis", 2) or 2)
    except Exception:
        naxis = 2

    for i in range(max(0, min(2, naxis))):
        try:
            out.append(str(cwcs.ctype[i]))
        except Exception:
            break

    if len(out) < 2:
        try:
            chdr = wcs.celestial.to_header()
            for key in ("CTYPE1", "CTYPE2"):
                if key in chdr and len(out) < 2:
                    out.append(str(chdr[key]))
        except Exception:
            pass

    while len(out) < 2:
        out.append("")
    return out[:2]



def _infer_native_coord_system(wcs: WCS) -> str:
    cwcs = wcs.celestial.wcs

    lngtyp = str(getattr(cwcs, "lngtyp", "") or "").upper()
    lattyp = str(getattr(cwcs, "lattyp", "") or "").upper()
    if lngtyp == "GLON" or lattyp == "GLAT":
        return "galactic"
    if lngtyp == "RA" or lattyp == "DEC":
        return "equatorial"

    ctype = [ct.upper() for ct in _get_celestial_ctype_list(wcs)]
    if any("GLON" in ct for ct in ctype) or any("GLAT" in ct for ct in ctype):
        return "galactic"
    if any("RA" in ct for ct in ctype) or any("DEC" in ct for ct in ctype):
        return "equatorial"
    return "celestial"



def _default_axis_labels(wcs: WCS) -> Tuple[str, str]:
    mode = _infer_native_coord_system(wcs)
    if mode == "galactic":
        return "Galactic Longitude", "Galactic Latitude"
    if mode == "equatorial":
        return "Right Ascension", "Declination"
    ctype = _get_celestial_ctype_list(wcs)
    return ctype[0], ctype[1]



def _pixel_scale_arcsec(wcs: WCS) -> float:
    scales = proj_plane_pixel_scales(wcs.celestial)
    return float(np.nanmean(np.abs(scales)) * 3600.0)



def _infer_orig_hpbw_arcsec(header: Optional[fits.Header]) -> Tuple[Optional[float], Optional[str]]:
    """Infer original HPBW from common FITS header conventions.

    Returns
    -------
    (value_arcsec, source_name)
    """
    if header is None:
        return None, None

    # Standard FITS beam keywords: degrees.
    bmaj = header.get("BMAJ")
    bmin = header.get("BMIN")
    if bmaj is not None or bmin is not None:
        vals = [float(v) for v in (bmaj, bmin) if v is not None]
        if vals:
            return float(np.mean(vals) * 3600.0), "BMAJ/BMIN(deg)"

    # Less standard but unambiguous arcsec-like keywords.
    for key in ("HPBW_ARCSEC", "HPBW_ASEC", "HPBWASEC", "ORIGHPBW", "BEAM_ARCSEC"):
        if key in header:
            try:
                return float(header[key]), key
            except Exception:
                pass

    return None, None


def _validate_smoothing_arguments(
    *,
    smooth_fwhm_arcsec: Optional[float] = None,
    target_hpbw_arcsec: Optional[float] = None,
) -> None:
    """Validate mutually exclusive smoothing specifications."""
    if smooth_fwhm_arcsec is not None and target_hpbw_arcsec is not None:
        raise ValueError(
            "smooth_fwhm_arcsec and target_hpbw_arcsec are alternative ways to specify 2D smoothing; "
            "please provide only one of them."
        )


def _validate_fill_flags(
    *,
    zero_fill: bool = False,
    nan_fill: bool = True,
) -> None:
    """Prevent ambiguous mask-fill requests."""
    if zero_fill and nan_fill:
        raise ValueError("zero_fill and nan_fill cannot both be True at the same time.")


def _effective_hpbw_arcsec(
    *,
    smooth_fwhm_arcsec: Optional[float] = None,
    target_hpbw_arcsec: Optional[float] = None,
    orig_hpbw_arcsec: Optional[float] = None,
) -> Tuple[Optional[float], Optional[str]]:
    """Estimate the effective HPBW after requested smoothing.

    Returns
    -------
    (effective_hpbw_arcsec, description)

    Notes
    -----
    - If target_hpbw_arcsec is given together with orig_hpbw_arcsec, the target is
      the effective beam.
    - If only smooth_fwhm_arcsec is given and the original beam is known, the beams
      are combined in quadrature.
    - If the original beam is unknown, the effective beam is generally unknown. In
      that case None is returned even though smoothing may still be applied.
    """
    if target_hpbw_arcsec is not None and orig_hpbw_arcsec is not None:
        return float(target_hpbw_arcsec), "target_hpbw"
    if smooth_fwhm_arcsec is not None and orig_hpbw_arcsec is not None:
        eff = np.sqrt(float(orig_hpbw_arcsec) ** 2 + float(smooth_fwhm_arcsec) ** 2)
        return float(eff), "orig_plus_smooth_quadrature"
    if target_hpbw_arcsec is None and smooth_fwhm_arcsec is None and orig_hpbw_arcsec is not None:
        return float(orig_hpbw_arcsec), "original"
    return None, None



def _get_hdu_from_hdul(hdul: fits.HDUList, ext: Optional[Union[int, str]] = None) -> fits.ImageHDU:
    if ext is None:
        for hdu in hdul:
            if getattr(hdu, "data", None) is not None:
                return hdu
        raise ValueError("No image HDU with data found in the FITS file.")
    return hdul[ext]



def _has_ext(hdul: fits.HDUList, ext: Union[int, str]) -> bool:
    try:
        hdul[ext]
        return True
    except Exception:
        return False



def _open_hdul_if_needed(source: SourceLike) -> Tuple[Optional[fits.HDUList], bool]:
    if isinstance(source, fits.HDUList):
        return source, False
    if _is_pathlike(source):
        return fits.open(source), True
    return None, False


def _is_otf_bundle(source: Any) -> bool:
    return OTFBundle is not None and isinstance(source, OTFBundle)


def _bundle_ext_name(ext: Optional[Union[int, str]]) -> str:
    if ext is None:
        return "DATA"
    if isinstance(ext, int):
        if ext == 0:
            return "DATA"
        raise ValueError(f"OTFBundle does not support integer ext={ext!r} except 0 for primary/data.")
    return str(ext).strip()


def _bundle_find_image_ext_key(bundle: Any, ext: str) -> Optional[str]:
    target = str(ext).strip().upper()
    for key in bundle.image_ext.keys():
        if str(key).strip().upper() == target:
            return key
    return None


def _bundle_spectral_axis_unit(header: fits.Header) -> Optional[u.UnitBase]:
    for key in ("CUNIT3", "CUNIT1"):
        value = header.get(key)
        if value not in (None, ""):
            unit = _as_unit(value)
            if unit is not None:
                return unit
    return None


def _bundle_unit_for_ext(bundle: Any, ext_name: str, arr: np.ndarray) -> Optional[u.UnitBase]:
    ext_upper = str(ext_name).strip().upper()
    data_unit = _as_unit(getattr(bundle, "unit", None)) or _get_bunit_from_header(bundle.header)
    spec_unit = _bundle_spectral_axis_unit(bundle.header)

    if ext_upper in {"DATA", "PRIMARY", "CUBE"}:
        return data_unit
    if ext_upper in {"VARIANCE"}:
        return None if data_unit is None else data_unit * data_unit
    if ext_upper in {"RMS", "BASE_RMS", "MOSAIC_RMS_OBS", "MOSAIC_RMS", "RMS_MAP_EMP", "RMS_MAP_INPUT_0", "TSYS"}:
        return data_unit
    if ext_upper in {"MOMENT0", "MOM0_BASESUP", "MOM0_LINECAND", "MOM0_3D_MASK"}:
        return None if (data_unit is None or spec_unit is None) else data_unit * spec_unit
    if ext_upper in {"MOMENT1"}:
        return spec_unit
    if ext_upper in {"MOMENT2"}:
        return None if spec_unit is None else spec_unit * spec_unit
    if ext_upper in {
        "WEIGHT", "WEIGHT_SUM", "HIT", "NSAMP", "WSUM", "WABS", "CANCEL", "WREL",
        "MASK", "VALID_MASK", "SUPPORT_MASK", "VALID_MASK_AND", "VALID_MASK_UNION",
        "LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR", "SIGNAL_MASK_USED",
        "BASE_FLG", "MASK3D", "LINECAND3D", "BASESUP3D",
        "MOSAIC_GAIN", "MOSAIC_TRUST", "MOSAIC_WEIGHT", "MOSAIC_INFO",
    }:
        return None
    if arr.ndim == 3:
        return data_unit
    return None


def _bundle_header_for_array(bundle: Any, arr: np.ndarray, *, unit: Optional[u.UnitBase]) -> fits.Header:
    if arr.ndim == 2:
        hdr = _drop_nonspatial_wcs_keywords(bundle.header)
    else:
        hdr = _copy_header(bundle.header)
    if unit is not None:
        try:
            hdr["BUNIT"] = unit.to_string()
        except Exception:
            pass
    elif "BUNIT" in hdr:
        del hdr["BUNIT"]
    return hdr


def _bundle_lookup_array(
    bundle: Any,
    ext: Optional[Union[int, str]] = None,
) -> Tuple[np.ndarray, fits.Header, Optional[u.UnitBase], Dict[str, Any]]:
    ext_name = _bundle_ext_name(ext)
    ext_upper = ext_name.upper()

    if ext_upper in {"DATA", "PRIMARY", "CUBE"}:
        arr = np.asarray(bundle.data)
    elif ext_upper == "VARIANCE":
        if bundle.variance is None:
            raise ValueError("OTFBundle has no VARIANCE cube.")
        arr = np.asarray(bundle.variance)
    elif ext_upper == "VALID_MASK":
        if bundle.valid_mask is not None:
            arr = np.asarray(bundle.valid_mask)
        else:
            key = _bundle_find_image_ext_key(bundle, ext_upper)
            if key is None:
                raise ValueError("OTFBundle has no VALID_MASK array.")
            arr = np.asarray(bundle.image_ext[key])
    elif ext_upper == "SUPPORT_MASK":
        if bundle.support_mask is not None:
            arr = np.asarray(bundle.support_mask)
        else:
            key = _bundle_find_image_ext_key(bundle, ext_upper)
            if key is None:
                raise ValueError("OTFBundle has no SUPPORT_MASK array.")
            arr = np.asarray(bundle.image_ext[key])
    else:
        key = _bundle_find_image_ext_key(bundle, ext_upper)
        if key is None:
            raise ValueError(f"Extension '{ext_name}' was not found in the OTFBundle.")
        arr = np.asarray(bundle.image_ext[key])
        ext_name = str(key)

    unit = _bundle_unit_for_ext(bundle, ext_name, arr)
    header = _bundle_header_for_array(bundle, arr, unit=unit)
    meta = {
        "source_kind": "otf_bundle",
        "ext": ext_name,
        "family_label": getattr(bundle, "family_label", None),
    }
    if isinstance(getattr(bundle, "meta", None), dict):
        baseline_viewer_mode = bundle.meta.get("baseline_viewer_mode")
        if baseline_viewer_mode is not None:
            meta["baseline_viewer_mode"] = baseline_viewer_mode
        baseline_viewer_mask_key = bundle.meta.get("baseline_viewer_mask_key")
        if baseline_viewer_mask_key is not None:
            meta["baseline_viewer_mask_key"] = baseline_viewer_mask_key
    return arr, header, unit, meta


def _projection_to_map2d(obj: Any, meta: Optional[Dict[str, Any]] = None) -> Map2D:
    header = None
    if hasattr(obj, "header"):
        try:
            header = fits.Header(obj.header)
        except Exception:
            header = None
    if header is None and hasattr(obj, "hdu"):
        header = fits.Header(obj.hdu.header)

    data = np.asarray(getattr(obj, "value", obj), dtype=float)
    if header is None and hasattr(obj, "wcs"):
        header = obj.wcs.celestial.to_header()
    header = _copy_header(header)
    unit = _as_unit(getattr(obj, "unit", None)) or _get_bunit_from_header(header)
    wcs = WCS(header).celestial if len(header) > 0 else obj.wcs.celestial
    if unit is not None and "BUNIT" not in header:
        header["BUNIT"] = unit.to_string()
    return Map2D(data=data, header=header, wcs=wcs, unit=unit, meta=dict(meta or {}))


# -----------------------------------------------------------------------------
# Input resolution
# -----------------------------------------------------------------------------


def resolve_map_input(
    source: SourceLike = None,
    *,
    data: Optional[np.ndarray] = None,
    header: Optional[fits.Header] = None,
    ext: Optional[Union[int, str]] = None,
) -> Map2D:
    """Resolve 2D input from file, HDU, (header, data), Projection, Map2D, or OTFBundle."""
    if isinstance(source, Map2D):
        return source

    if _is_otf_bundle(source):
        arr, hdr, unit, meta = _bundle_lookup_array(source, ext=ext)
        arr = np.asarray(arr, dtype=float)
        if arr.ndim != 2:
            raise ValueError(
                f"Selected OTFBundle source is not 2D (ext={meta.get('ext')!r}, ndim={arr.ndim})."
            )
        hdr2d = _drop_nonspatial_wcs_keywords(hdr)
        return Map2D(data=arr, header=hdr2d, wcs=WCS(hdr2d).celestial, unit=unit, meta=meta)

    if data is not None:
        if header is None:
            raise ValueError("header must be given together with data.")
        arr = np.asarray(data, dtype=float)
        if arr.ndim != 2:
            raise ValueError(f"Expected 2D data, got ndim={arr.ndim}.")
        hdr = _drop_nonspatial_wcs_keywords(header)
        unit = _get_bunit_from_header(header)
        return Map2D(
            data=arr,
            header=hdr,
            wcs=WCS(hdr).celestial,
            unit=unit,
            meta={"source_kind": "header_data"},
        )

    if isinstance(source, tuple) and len(source) == 2:
        return resolve_map_input(data=np.asarray(source[1], dtype=float), header=source[0], ext=ext)

    if isinstance(source, (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)):
        arr = np.asarray(source.data, dtype=float)
        if arr.ndim != 2:
            raise ValueError(f"Expected a 2D HDU, got ndim={arr.ndim}.")
        hdr = _drop_nonspatial_wcs_keywords(source.header)
        unit = _get_bunit_from_header(source.header)
        return Map2D(data=arr, header=hdr, wcs=WCS(hdr).celestial, unit=unit, meta={"source_kind": "hdu"})

    if Projection is not None and isinstance(source, Projection):
        return _projection_to_map2d(source, meta={"source_kind": "projection"})

    if hasattr(source, "wcs") and hasattr(source, "shape") and getattr(source, "ndim", None) == 2:
        return _projection_to_map2d(source, meta={"source_kind": type(source).__name__})

    hdul, close_after = _open_hdul_if_needed(source)
    if hdul is not None:
        try:
            hdu = _get_hdu_from_hdul(hdul, ext=ext)
            arr = np.asarray(hdu.data, dtype=float)
            if arr.ndim != 2:
                raise ValueError(f"Selected FITS HDU is not 2D (ndim={arr.ndim}).")
            hdr = _drop_nonspatial_wcs_keywords(hdu.header)
            unit = _get_bunit_from_header(hdu.header)
            return Map2D(
                data=arr,
                header=hdr,
                wcs=WCS(hdr).celestial,
                unit=unit,
                meta={"source_kind": "fits", "source": str(source), "ext": ext},
            )
        finally:
            if close_after:
                hdul.close()

    raise TypeError(
        "Could not interpret source as a 2D map. Supported inputs are file path, HDUList, "
        "2D HDU, (header, data), Projection-like object, Map2D, or OTFBundle."
    )


def _try_resolve_map_input(
    source: SourceLike = None,
    *,
    data: Optional[np.ndarray] = None,
    header: Optional[fits.Header] = None,
    ext: Optional[Union[int, str]] = None,
) -> Optional[Map2D]:
    """Best-effort 2D resolver that returns None instead of raising.

    This is used to detect whether the selected input/HDU is already 2D, so that
    make_2d_map() can pass it through even when mode='moment0' etc. are given.
    """
    try:
        return resolve_map_input(source=source, data=data, header=header, ext=ext)
    except Exception:
        return None



def _make_cube_from_header_and_data(header: fits.Header, data: np.ndarray, *, unit: Optional[u.UnitBase]) -> SpectralCube:
    qdata = data * unit if unit is not None else data
    return SpectralCube(data=qdata, wcs=WCS(header), header=header, meta={})



def resolve_cube_input(
    source: SourceLike,
    *,
    ext: Optional[Union[int, str]] = None,
) -> CubeInput:
    """Resolve a 3D cube from file, HDU, (header, data), SpectralCube, or OTFBundle."""
    if isinstance(source, SpectralCube):
        header = getattr(source, "header", None)
        if header is None:
            header = source.wcs.to_header()
        unit = _as_unit(getattr(source, "unit", None)) or _get_bunit_from_header(header)
        return CubeInput(cube=source, header=_copy_header(header), unit=unit, meta={"source_kind": "spectral_cube"})

    if isinstance(source, tuple) and len(source) == 2:
        header, data = source
        arr = np.asarray(data, dtype=float)
        if arr.ndim != 3:
            raise ValueError(f"Expected 3D data, got ndim={arr.ndim}.")
        hdr = _copy_header(header)
        unit = _get_bunit_from_header(hdr)
        cube = _make_cube_from_header_and_data(hdr, arr, unit=unit)
        return CubeInput(cube=cube, header=hdr, unit=unit, meta={"source_kind": "header_data"})

    if _is_otf_bundle(source):
        arr, hdr, unit, meta = _bundle_lookup_array(source, ext=ext)
        arr = np.asarray(arr, dtype=float)
        if arr.ndim != 3:
            raise ValueError(
                f"Selected OTFBundle source is not 3D (ext={meta.get('ext')!r}, ndim={arr.ndim})."
            )
        cube = _make_cube_from_header_and_data(hdr, arr, unit=unit)
        return CubeInput(cube=cube, header=hdr, unit=unit, meta=meta)

    if isinstance(source, (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)):
        arr = np.asarray(source.data, dtype=float)
        if arr.ndim != 3:
            raise ValueError(f"Expected a 3D HDU, got ndim={arr.ndim}.")
        hdr = _copy_header(source.header)
        unit = _get_bunit_from_header(hdr)
        cube = _make_cube_from_header_and_data(hdr, arr, unit=unit)
        return CubeInput(cube=cube, header=hdr, unit=unit, meta={"source_kind": "hdu"})

    hdul, close_after = _open_hdul_if_needed(source)
    if hdul is not None:
        try:
            hdu = _get_hdu_from_hdul(hdul, ext=ext)
            arr = np.asarray(hdu.data, dtype=float)
            if arr.ndim != 3:
                raise ValueError(f"Selected FITS HDU is not 3D (ndim={arr.ndim}).")
            hdr = _copy_header(hdu.header)
            unit = _get_bunit_from_header(hdr)
            cube = _make_cube_from_header_and_data(hdr, arr, unit=unit)
            return CubeInput(
                cube=cube,
                header=hdr,
                unit=unit,
                meta={"source_kind": "fits", "source": str(source), "ext": ext},
            )
        finally:
            if close_after:
                hdul.close()

    raise TypeError(
        "Could not interpret source as a 3D cube. Supported inputs are file path, HDUList, "
        "3D HDU, (header, data), SpectralCube, or OTFBundle."
    )


# -----------------------------------------------------------------------------
# Normalization and smoothing
# -----------------------------------------------------------------------------


def build_normalize(
    data: np.ndarray,
    *,
    mode: str = "asinh",
    percentile: Optional[Tuple[float, float]] = None,
    cmin: Optional[float] = None,
    cmax: Optional[float] = None,
    stretch_a: float = 0.1,
    power_gamma: float = 1.0,
    invalid: float = np.nan,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> ImageNormalize:
    """Build an ImageNormalize instance from data and display options.

    Parameters
    ----------
    cmin, cmax : float or None
        Display-range lower/upper limits for the color scale.
        These are not velocity limits. Use ``vel_range`` / ``chan_range`` for
        spectral selection.

    vmin, vmax : float or None
        Deprecated aliases of ``cmin`` / ``cmax`` kept for backward compatibility.
    """
    if vmin is not None or vmax is not None:
        warnings.warn(
            "vmin/vmax are deprecated and mean display color-scale limits; use cmin/cmax instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if cmin is None and vmin is not None:
            cmin = vmin
        if cmax is None and vmax is not None:
            cmax = vmax

    arr = np.asarray(data, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        finite = np.array([0.0])

    if percentile is not None:
        p_lo, p_hi = percentile
        auto_cmin, auto_cmax = np.nanpercentile(finite, [p_lo, p_hi])
        if cmin is None:
            cmin = float(auto_cmin)
        if cmax is None:
            cmax = float(auto_cmax)

    if cmin is None:
        cmin = float(np.nanmin(finite))
    if cmax is None:
        cmax = float(np.nanmax(finite))

    if not np.isfinite(cmin):
        cmin = 0.0
    if not np.isfinite(cmax):
        cmax = 1.0
    if cmax == cmin:
        cmax = cmin + 1.0

    mode = str(mode).lower()
    if mode == "linear":
        stretch = LinearStretch()
    elif mode == "sqrt":
        stretch = SqrtStretch()
    elif mode == "log":
        positive = finite[finite > 0]
        if positive.size == 0:
            raise ValueError("log normalization requested, but data contain no positive finite values.")
        smallest_positive = float(np.nanmin(positive))
        if cmin <= 0:
            cmin = smallest_positive
        stretch = LogStretch()
    elif mode == "asinh":
        stretch = AsinhStretch(a=stretch_a)
    elif mode == "power":
        stretch = PowerStretch(power_gamma)
    else:
        raise ValueError(f"Unsupported normalization mode: {mode}")

    return ImageNormalize(vmin=cmin, vmax=cmax, stretch=stretch, invalid=invalid)



def apply_gaussian_smoothing_2d(
    data: np.ndarray,
    wcs: WCS,
    *,
    smooth_fwhm_arcsec: Optional[float] = None,
    target_hpbw_arcsec: Optional[float] = None,
    orig_hpbw_arcsec: Optional[float] = None,
    boundary: str = "extend",
) -> np.ndarray:
    """Apply 2D Gaussian smoothing specified in arcsec."""
    _validate_smoothing_arguments(
        smooth_fwhm_arcsec=smooth_fwhm_arcsec,
        target_hpbw_arcsec=target_hpbw_arcsec,
    )
    if smooth_fwhm_arcsec is None and target_hpbw_arcsec is None:
        return np.asarray(data, dtype=float)

    arr = np.asarray(data, dtype=float)
    pix_arcsec = _pixel_scale_arcsec(wcs)

    kernel_fwhm_arcsec = smooth_fwhm_arcsec
    if target_hpbw_arcsec is not None:
        if orig_hpbw_arcsec is None:
            warnings.warn(
                "target_hpbw_arcsec was requested, but orig_hpbw_arcsec could not be inferred. "
                "Treating target_hpbw_arcsec as the smoothing kernel FWHM.",
                RuntimeWarning,
            )
            kernel_fwhm_arcsec = target_hpbw_arcsec if kernel_fwhm_arcsec is None else max(
                float(kernel_fwhm_arcsec), float(target_hpbw_arcsec)
            )
        else:
            delta2 = float(target_hpbw_arcsec) ** 2 - float(orig_hpbw_arcsec) ** 2
            if delta2 <= 0:
                return arr.copy()
            kernel_from_target = float(np.sqrt(delta2))
            if kernel_fwhm_arcsec is None:
                kernel_fwhm_arcsec = kernel_from_target
            else:
                kernel_fwhm_arcsec = max(float(kernel_fwhm_arcsec), kernel_from_target)

    if kernel_fwhm_arcsec is None or kernel_fwhm_arcsec <= 0:
        return arr.copy()

    sigma_pix = (float(kernel_fwhm_arcsec) / pix_arcsec) * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(x_stddev=sigma_pix)
    return convolve(arr, kernel, boundary=boundary, preserve_nan=True, normalize_kernel=True)


# -----------------------------------------------------------------------------
# Description / cropping / convenience helpers
# -----------------------------------------------------------------------------


def _safe_header_get_float(header: Optional[fits.Header], key: str) -> Optional[float]:
    if header is None:
        return None
    try:
        value = header.get(key)
        return None if value is None else float(value)
    except Exception:
        return None


def _format_range_value(value: Any) -> str:
    if value is None:
        return "none"
    if isinstance(value, tuple) and len(value) == 2:
        return f"({value[0]}, {value[1]})"
    return str(value)


def describe_map2d(map2d: Map2D, *, name: Optional[str] = None) -> str:
    """Return a short human-readable description of a 2D map."""
    lines: List[str] = []
    if name:
        lines.append(f"name: {name}")
    lines.append(f"shape: {map2d.data.shape[1]} x {map2d.data.shape[0]} pixels")
    lines.append(f"coord_system: {_infer_native_coord_system(map2d.wcs)}")
    if map2d.unit is not None:
        lines.append(f"unit: {map2d.unit.to_string()}")
    mode = map2d.meta.get("mode")
    if mode is not None:
        lines.append(f"mode: {mode}")
    range_kind = map2d.meta.get("range_kind")
    range_value = map2d.meta.get("range_value")
    if range_kind is not None:
        lines.append(f"range: {range_kind} {_format_range_value(range_value)}")
    smooth_info = map2d.meta.get("smooth_info", {}) if isinstance(map2d.meta, dict) else {}
    eff_hpbw = smooth_info.get("effective_hpbw_arcsec")
    if eff_hpbw is None:
        bmaj = _safe_header_get_float(map2d.header, "BMAJ")
        bmin = _safe_header_get_float(map2d.header, "BMIN")
        if bmaj is not None:
            if bmin is None:
                bmin = bmaj
            lines.append(f"beam: {bmaj * 3600.0:.3g} x {bmin * 3600.0:.3g} arcsec")
    else:
        lines.append(f"effective_beam: {float(eff_hpbw):.3g} arcsec")

    finite = np.asarray(map2d.data, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size > 0:
        lines.append(
            "value_range: "
            f"min={float(np.nanmin(finite)):.6g}, max={float(np.nanmax(finite)):.6g}, median={float(np.nanmedian(finite)):.6g}"
        )
    return "\n".join(lines)


def describe_source(
    source: SourceLike = None,
    *,
    data: Optional[np.ndarray] = None,
    header: Optional[fits.Header] = None,
    ext: Optional[Union[int, str]] = None,
    mode: str = "moment0",
    chan_range: Optional[Tuple[int, int]] = None,
    vel_range: Optional[Tuple[float, float]] = None,
    spectral_unit: Union[str, u.UnitBase] = "km/s",
    mask: Optional[np.ndarray] = None,
    mask_mode: Optional[str] = None,
    zero_fill: bool = False,
    nan_fill: bool = True,
    smooth_fwhm_arcsec: Optional[float] = None,
    target_hpbw_arcsec: Optional[float] = None,
    orig_hpbw_arcsec: Optional[float] = None,
    linefree_ext: str = "LINEFREE",
    linecand_ext: str = "LINECAND3D",
    basesup_ext: str = "BASESUP3D",
    final_mask_ext: str = "MASK3D",
) -> str:
    """Resolve the source to 2D and return a short description.

    Either ``source`` or ``(header, data)`` may be supplied. When ``data`` is
    given, ``header`` must also be given. If both are supplied, ``(header, data)``
    takes precedence because it is the more explicit in-memory representation.
    """
    actual_source = source
    if data is not None:
        if header is None:
            raise ValueError("header must be given together with data.")
        if source is not None:
            warnings.warn(
                "describe_source received both source=... and data/header=...; using the in-memory (header, data) pair.",
                RuntimeWarning,
            )
        actual_source = (header, data)

    map2d = make_2d_map(
        source=actual_source,
        ext=ext,
        mode=mode,
        chan_range=chan_range,
        vel_range=vel_range,
        spectral_unit=spectral_unit,
        mask=mask,
        mask_mode=mask_mode,
        zero_fill=zero_fill,
        nan_fill=nan_fill,
        smooth_fwhm_arcsec=smooth_fwhm_arcsec,
        target_hpbw_arcsec=target_hpbw_arcsec,
        orig_hpbw_arcsec=orig_hpbw_arcsec,
        linefree_ext=linefree_ext,
        linecand_ext=linecand_ext,
        basesup_ext=basesup_ext,
        final_mask_ext=final_mask_ext,
    )
    return describe_map2d(map2d)
def plotting_help_text() -> str:
    """Return a compact programmatic help text for the plotting helpers."""
    return (
        "Minimal API\n"
        "-----------\n"
        "quicklook(source, ...)\n"
        "    Fast entry point for showing one 2D map or one reduced 3D cube with sensible defaults.\n\n"
        "plot_scene(base=..., overlays=[...], ...)\n"
        "    High-level API with one base layer and optional overlay layers.\n"
        "    Supported overlay kinds: image, contour, marker, catalog, text, beam, scalebar, north_arrow.\n"
        "    catalog layers can use coords/lon+lat/xy directly, or table-like input (pandas / astropy.table / dict of arrays).\n\n"
        "describe_source(source, ...)\n"
        "    Return a short text summary of what would be plotted.\n\n"
        "Typical quicklook examples\n"
        "--------------------------\n"
        "quicklook('cube.fits')\n"
        "quicklook('cube.fits', vel_range=(20, 40), title='12CO integrated intensity')\n"
        "quicklook((header, data), grid=False, save='map.pdf')\n"
        "quicklook('cube.fits', contours=[{'levels': 'auto', 'colors': 'white'}])\n\n"
        "Typical scene example\n"
        "---------------------\n"
        "plot_scene(\n"
        "    base={'kind': 'image', 'source': 'optical.fits'},\n"
        "    overlays=[\n"
        "        {'kind': 'contour', 'source': 'co.fits', 'mode': 'moment0', 'vel_range': (20, 40), 'colors': 'cyan'},\n"
        "        {'kind': 'beam', 'beam': 'auto'},\n"
        "    ],\n"
        "    title='Optical + CO contours',\n"
        ")\n"
    )


def print_plotting_help() -> None:
    """Print the compact programmatic help text."""
    print(plotting_help_text())


def crop_map2d(
    map2d: Map2D,
    *,
    center: Optional[Any] = None,
    size: Optional[Any] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    mode: str = "trim",
    fill_value: float = np.nan,
) -> Map2D:
    """Return a cropped cutout of a 2D map.

    Supported forms
    ---------------
    1. center + size
       - center: (xpix, ypix) or SkyCoord
       - size  : scalar / (ny, nx) in pixels, or astropy Quantity for world-size cutouts
    2. xlim + ylim
       - pixel-coordinate bounds; converted internally to center + size
    """
    if center is None and size is None and xlim is None and ylim is None:
        return map2d

    if (xlim is None) ^ (ylim is None):
        raise ValueError("xlim and ylim must be provided together.")

    position = center
    cutout_size = size
    if xlim is not None and ylim is not None:
        x0, x1 = float(xlim[0]), float(xlim[1])
        y0, y1 = float(ylim[0]), float(ylim[1])
        xc = 0.5 * (x0 + x1)
        yc = 0.5 * (y0 + y1)
        nx = max(1, int(round(abs(x1 - x0))))
        ny = max(1, int(round(abs(y1 - y0))))
        position = (xc, yc)
        cutout_size = (ny, nx)

    if position is None or cutout_size is None:
        raise ValueError("Either center+size or xlim+ylim must be provided for cropping.")

    cutout = Cutout2D(
        data=np.asarray(map2d.data, dtype=float),
        position=position,
        size=cutout_size,
        wcs=map2d.wcs,
        mode=mode,
        fill_value=fill_value,
    )
    hdr = fits.Header(map2d.header)
    hdr.update(cutout.wcs.to_header())
    meta = dict(map2d.meta)
    meta["crop"] = {
        "center": repr(center),
        "size": repr(size),
        "xlim": xlim,
        "ylim": ylim,
        "mode": mode,
    }
    return Map2D(
        data=np.asarray(cutout.data, dtype=float),
        header=hdr,
        wcs=cutout.wcs.celestial,
        unit=map2d.unit,
        meta=meta,
    )


# -----------------------------------------------------------------------------
# Spectral helpers
# -----------------------------------------------------------------------------


def _cube_data_array(cube: SpectralCube) -> np.ndarray:
    data = cube.unmasked_data[:]
    if hasattr(data, "value"):
        return np.asarray(data.value, dtype=float)
    return np.asarray(data, dtype=float)



def _cube_unit(cube: SpectralCube, fallback: Optional[u.UnitBase] = None) -> Optional[u.UnitBase]:
    unit = _as_unit(getattr(cube, "unit", None))
    return unit if unit is not None else fallback



def _native_spectral_axis_and_data(cube_input: CubeInput) -> Tuple[SpectralCube, np.ndarray, u.UnitBase, np.ndarray, Optional[u.UnitBase]]:
    cube = cube_input.cube
    spectral_axis = cube.spectral_axis
    data = _cube_data_array(cube)
    return cube, spectral_axis.value, spectral_axis.unit, data, _cube_unit(cube, fallback=cube_input.unit)



def _velocity_spectral_axis_and_data(
    cube_input: CubeInput,
    *,
    spectral_unit: Union[str, u.UnitBase] = "km/s",
) -> Tuple[SpectralCube, np.ndarray, u.UnitBase, np.ndarray, Optional[u.UnitBase]]:
    cube2 = cube_input.cube.with_spectral_unit(u.Unit(spectral_unit), velocity_convention="radio")
    spectral_axis = cube2.spectral_axis
    data = _cube_data_array(cube2)
    return cube2, spectral_axis.value, spectral_axis.unit, data, _cube_unit(cube2, fallback=cube_input.unit)



def _chan_indices_from_range(nchan: int, chan_range: Optional[Tuple[int, int]]) -> np.ndarray:
    if chan_range is None:
        return np.arange(nchan, dtype=int)
    lo, hi = chan_range
    lo = int(lo)
    hi = int(hi)
    if hi < lo:
        lo, hi = hi, lo
    lo = max(0, lo)
    hi = min(nchan - 1, hi)
    if hi < lo:
        raise ValueError("chan_range selects no channels.")
    return np.arange(lo, hi + 1, dtype=int)



def _vel_indices_from_range(spec_axis: np.ndarray, vel_range: Optional[Tuple[float, float]]) -> np.ndarray:
    if vel_range is None:
        return np.arange(spec_axis.size, dtype=int)
    lo, hi = float(vel_range[0]), float(vel_range[1])
    if hi < lo:
        lo, hi = hi, lo
    mask = (spec_axis >= lo) & (spec_axis <= hi)
    idx = np.where(mask)[0]
    if idx.size == 0:
        raise ValueError("vel_range selects no channels in the current spectral axis.")
    return idx.astype(int)



def _resolve_channel_indices(
    nchan: int,
    spec_axis: np.ndarray,
    *,
    chan_range: Optional[Tuple[int, int]] = None,
    vel_range: Optional[Tuple[float, float]] = None,
) -> np.ndarray:
    if chan_range is not None:
        return _chan_indices_from_range(nchan, chan_range)
    return _vel_indices_from_range(spec_axis, vel_range)



def _spectral_deltas(spec_axis: np.ndarray, spec_unit: u.UnitBase) -> u.Quantity:
    vals = np.asarray(spec_axis, dtype=float)
    if vals.size == 1:
        return np.array([1.0]) * spec_unit
    dvals = np.abs(np.gradient(vals))
    return dvals * spec_unit



def _moment0_from_subset(
    data3d: np.ndarray,
    spec_axis: np.ndarray,
    spec_unit: u.UnitBase,
    *,
    chan_idx: np.ndarray,
    data_unit: Optional[u.UnitBase],
    nan_fill: bool = True,
) -> Tuple[np.ndarray, Optional[u.UnitBase]]:
    sub = np.asarray(data3d[chan_idx, :, :], dtype=float)
    dv = _spectral_deltas(spec_axis[chan_idx], spec_unit).value[:, None, None]
    out = np.nansum(sub * dv, axis=0)
    if nan_fill:
        all_nan = np.all(~np.isfinite(sub), axis=0)
        out[all_nan] = np.nan
    out_unit = None if data_unit is None else data_unit * spec_unit
    return out, out_unit



def _linefree_to_channel_mask(linefree: np.ndarray, nchan: int) -> np.ndarray:
    arr = np.asarray(linefree)
    if arr.ndim == 1:
        if arr.shape[0] != nchan:
            raise ValueError(f"LINEFREE length ({arr.shape[0]}) does not match cube channels ({nchan}).")
        return arr.astype(bool)
    if arr.ndim == 2 and 1 in arr.shape:
        vec = arr.ravel()
        if vec.size != nchan:
            raise ValueError(f"LINEFREE size ({vec.size}) does not match cube channels ({nchan}).")
        return vec.astype(bool)
    if arr.ndim == 3:
        if arr.shape[0] == nchan and arr.shape[1] == 1 and arr.shape[2] == 1:
            return arr[:, 0, 0].astype(bool)
        if arr.shape[0] == nchan:
            return np.all(arr.astype(bool), axis=(1, 2))
    raise ValueError("Unsupported LINEFREE shape. Expected 1D, 2D singleton, or 3D compatible array.")



def _mask_from_associated_hdu(
    source: SourceLike,
    *,
    ext: Union[int, str],
) -> np.ndarray:
    if _is_otf_bundle(source):
        arr, _hdr, _unit, meta = _bundle_lookup_array(source, ext=ext)
        return np.asarray(arr)

    hdul, close_after = _open_hdul_if_needed(source)
    if hdul is None:
        raise ValueError(f"Associated mask HDU '{ext}' requires source to be a FITS path, HDUList, or OTFBundle.")
    try:
        if not _has_ext(hdul, ext):
            raise ValueError(f"Extension '{ext}' was not found in the FITS source.")
        return np.asarray(hdul[ext].data)
    finally:
        if close_after:
            hdul.close()



def _select_default_provisional_mask_mode(
    source: SourceLike,
    *,
    linefree_ext: str,
    linecand_ext: str,
    basesup_ext: str,
) -> str:
    """Choose a provisional-moment mask mode automatically.

    Search order is:
        1. linecand3d       -> use True voxels directly as signal candidates
        2. basesup3d        -> use the complement of baseline-support voxels
        3. linefree_complement -> use channels not marked as line-free

    Rationale:
        LINECAND3D is the most direct statement of where signal is expected.
        BASESUP3D and LINEFREE describe baseline / line-free regions, so their complement
        is used as a provisional signal region when LINECAND3D is unavailable.
    """
    if _is_otf_bundle(source):
        if _bundle_find_image_ext_key(source, linecand_ext) is not None:
            return "linecand3d"
        if _bundle_find_image_ext_key(source, basesup_ext) is not None:
            return "basesup3d"
        if _bundle_find_image_ext_key(source, linefree_ext) is not None:
            return "linefree_complement"
        raise ValueError(
            f"No provisional-mask candidate found in the OTFBundle ({linecand_ext} / {basesup_ext} / {linefree_ext})."
        )

    hdul, close_after = _open_hdul_if_needed(source)
    if hdul is None:
        raise ValueError(
            "Automatic provisional moment requires source to be a FITS path, HDUList, or OTFBundle so that "
            "LINECAND3D / BASESUP3D / LINEFREE can be searched."
        )
    try:
        if _has_ext(hdul, linecand_ext):
            return "linecand3d"
        if _has_ext(hdul, basesup_ext):
            return "basesup3d"
        if _has_ext(hdul, linefree_ext):
            return "linefree_complement"
        raise ValueError(
            f"No provisional-mask candidate found ({linecand_ext} / {basesup_ext} / {linefree_ext})."
        )
    finally:
        if close_after:
            hdul.close()



def _build_mask3d(
    source: SourceLike,
    cube_shape: Tuple[int, int, int],
    *,
    mask_mode: str,
    external_mask: Optional[np.ndarray] = None,
    linefree_ext: str = "LINEFREE",
    linecand_ext: str = "LINECAND3D",
    basesup_ext: str = "BASESUP3D",
    final_mask_ext: str = "MASK3D",
) -> np.ndarray:
    nchan, ny, nx = cube_shape
    mode = str(mask_mode).lower()

    if mode == "external":
        if external_mask is None:
            raise ValueError("mask_mode='external' requires mask=... to be provided.")
        mask = np.asarray(external_mask).astype(bool)
    elif mode == "linefree_complement":
        linefree = _mask_from_associated_hdu(source, ext=linefree_ext)
        chmask = _linefree_to_channel_mask(linefree, nchan)
        mask = (~chmask)[:, None, None] * np.ones((1, ny, nx), dtype=bool)
    elif mode == "linecand3d":
        mask = _mask_from_associated_hdu(source, ext=linecand_ext).astype(bool)
    elif mode == "basesup3d":
        basesup = _mask_from_associated_hdu(source, ext=basesup_ext).astype(bool)
        mask = ~basesup
    elif mode == "mask3d":
        mask = _mask_from_associated_hdu(source, ext=final_mask_ext).astype(bool)
    else:
        raise ValueError(f"Unsupported mask_mode: {mask_mode}")

    if mask.shape != cube_shape:
        if mask.ndim == 1 and mask.size == nchan:
            mask = mask[:, None, None] * np.ones((1, ny, nx), dtype=bool)
        elif mask.ndim == 3 and mask.shape[0] == nchan and mask.shape[1:] == (1, 1):
            mask = mask * np.ones((1, ny, nx), dtype=bool)
        else:
            raise ValueError(f"Mask shape {mask.shape} is incompatible with cube shape {cube_shape}.")
    return mask.astype(bool)



def _apply_mask_fill(
    data3d: np.ndarray,
    mask3d: Optional[np.ndarray],
    *,
    zero_fill: bool = False,
    nan_fill: bool = True,
    fill_value: float = np.nan,
) -> np.ndarray:
    arr = np.asarray(data3d, dtype=float).copy()
    if mask3d is None:
        return arr
    if zero_fill:
        arr[~mask3d] = 0.0
    elif nan_fill:
        arr[~mask3d] = np.nan
    else:
        arr[~mask3d] = fill_value
    return arr


# -----------------------------------------------------------------------------
# 2D generation
# -----------------------------------------------------------------------------


def make_2d_map(
    source: SourceLike,
    *,
    ext: Optional[Union[int, str]] = None,
    mode: str = "moment0",
    chan_range: Optional[Tuple[int, int]] = None,
    vel_range: Optional[Tuple[float, float]] = None,
    spectral_unit: Union[str, u.UnitBase] = "km/s",
    mask: Optional[np.ndarray] = None,
    mask_mode: Optional[str] = None,
    zero_fill: bool = False,
    nan_fill: bool = True,
    smooth_fwhm_arcsec: Optional[float] = None,
    target_hpbw_arcsec: Optional[float] = None,
    orig_hpbw_arcsec: Optional[float] = None,
    fill_value: float = np.nan,
    linefree_ext: str = "LINEFREE",
    linecand_ext: str = "LINECAND3D",
    basesup_ext: str = "BASESUP3D",
    final_mask_ext: str = "MASK3D",
) -> Map2D:
    """Create a 2D map from 2D input or from a 3D cube.

    If the selected input/HDU is already 2D (for example ext="MOMENT0"), it is
    passed through as a precomputed map. In that case, 3D-reduction arguments
    such as vel_range, chan_range, and mask_mode are ignored with a warning.
    """
    mode = str(mode).lower()
    _validate_smoothing_arguments(
        smooth_fwhm_arcsec=smooth_fwhm_arcsec,
        target_hpbw_arcsec=target_hpbw_arcsec,
    )
    _validate_fill_flags(zero_fill=zero_fill, nan_fill=nan_fill)

    # 2D passthrough path:
    # - explicit 2D/identity mode
    # - or the selected input/HDU itself is already 2D (e.g. ext='MOMENT0')
    candidate_2d = _try_resolve_map_input(source, ext=ext)
    if mode in {"identity", "map", "2d"} or candidate_2d is not None:
        map2d = candidate_2d if candidate_2d is not None else resolve_map_input(source, ext=ext)

        ignored_args = []
        ignored_reason = {}
        if mode not in {"identity", "map", "2d"}:
            ignored_args.append(f"mode='{mode}'")
            ignored_reason["mode"] = (
                "the selected input/HDU is already 2D, so no 3D reduction is performed"
            )
        if chan_range is not None:
            ignored_args.append("chan_range")
            ignored_reason["chan_range"] = (
                "channel selection is only used when collapsing a 3D cube into 2D"
            )
        if vel_range is not None:
            ignored_args.append("vel_range")
            ignored_reason["vel_range"] = (
                "velocity selection is only used when collapsing a 3D cube into 2D"
            )
        if mask is not None:
            ignored_args.append("mask")
            ignored_reason["mask"] = (
                "mask application here is defined for 3D cube reduction, not for precomputed 2D maps"
            )
        if mask_mode is not None:
            ignored_args.append("mask_mode")
            ignored_reason["mask_mode"] = (
                "mask_mode selects how a 3D cube is reduced; it is not applied to a precomputed 2D map"
            )
        if ignored_args:
            msg = (
                "Selected input is already 2D (for example, ext='MOMENT0'). "
                "It will be passed through as a precomputed 2D map, and 3D-reduction arguments are ignored. "
                "Ignored arguments: " + ", ".join(ignored_args) + "."
            )
            if vel_range is not None or chan_range is not None:
                msg += (
                    " If you want to remake a moment map over a specific velocity/channel range, "
                    "select the original 3D cube HDU instead of a precomputed 2D extension."
                )
            warnings.warn(msg, RuntimeWarning)

        orig = orig_hpbw_arcsec
        orig_source = "explicit"
        if orig is None:
            orig, orig_source = _infer_orig_hpbw_arcsec(map2d.header)
        if smooth_fwhm_arcsec is not None or target_hpbw_arcsec is not None:
            smoothed = apply_gaussian_smoothing_2d(
                map2d.data,
                map2d.wcs,
                smooth_fwhm_arcsec=smooth_fwhm_arcsec,
                target_hpbw_arcsec=target_hpbw_arcsec,
                orig_hpbw_arcsec=orig,
            )
            eff_hpbw, eff_source = _effective_hpbw_arcsec(
                smooth_fwhm_arcsec=smooth_fwhm_arcsec,
                target_hpbw_arcsec=target_hpbw_arcsec,
                orig_hpbw_arcsec=orig,
            )
            hdr = fits.Header(map2d.header)
            if eff_hpbw is not None:
                hdr["BMAJ"] = float(eff_hpbw) / 3600.0
                hdr["BMIN"] = float(eff_hpbw) / 3600.0
            meta = dict(map2d.meta)
            meta.update(
                {
                    "selected_input_dim": 2,
                    "source_is_precomputed_2d": True,
                    "ignored_make_2d_args": ignored_args,
                    "ignored_make_2d_reason": ignored_reason,
                    "smooth_info": {
                        "smooth_fwhm_arcsec": smooth_fwhm_arcsec,
                        "target_hpbw_arcsec": target_hpbw_arcsec,
                        "orig_hpbw_arcsec": orig,
                        "orig_hpbw_source": orig_source,
                        "effective_hpbw_arcsec": eff_hpbw,
                        "effective_hpbw_source": eff_source,
                    },
                }
            )
            return Map2D(data=smoothed, header=hdr, wcs=WCS(hdr).celestial, unit=map2d.unit, meta=meta)

        meta = dict(map2d.meta)
        meta.update(
            {
                "selected_input_dim": 2,
                "source_is_precomputed_2d": True,
                "ignored_make_2d_args": ignored_args,
                "ignored_make_2d_reason": ignored_reason,
            }
        )
        return Map2D(data=np.asarray(map2d.data, dtype=float), header=fits.Header(map2d.header), wcs=WCS(map2d.header).celestial, unit=map2d.unit, meta=meta)

    cube_input = resolve_cube_input(source, ext=ext)
    cube = cube_input.cube
    data3d_full = _cube_data_array(cube)

    use_velocity_axis = (vel_range is not None) or (mode == "velocity_slice") or (mode in {"moment0", "provisional_moment", "final_moment"})
    spectral_axis_mode = "native"
    spectral_axis_warning = None
    if use_velocity_axis:
        try:
            _, spec_axis, spec_unit, data3d, data_unit = _velocity_spectral_axis_and_data(
                cube_input,
                spectral_unit=spectral_unit,
            )
            spectral_axis_mode = "velocity"
        except Exception as exc:
            if vel_range is not None or mode == "velocity_slice":
                raise
            warnings.warn(
                "Velocity-axis conversion failed for this cube; falling back to the native spectral axis. "
                f"Original error: {exc}",
                RuntimeWarning,
            )
            _, spec_axis, spec_unit, data3d, data_unit = _native_spectral_axis_and_data(cube_input)
            spectral_axis_mode = "native_fallback"
            spectral_axis_warning = str(exc)
    else:
        _, spec_axis, spec_unit, data3d, data_unit = _native_spectral_axis_and_data(cube_input)

    mask3d = None
    moment_kind = "plain"
    if mode == "provisional_moment":
        moment_kind = "provisional"
        if mask_mode is None:
            mask_mode = _select_default_provisional_mask_mode(
                source,
                linefree_ext=linefree_ext,
                linecand_ext=linecand_ext,
                basesup_ext=basesup_ext,
            )
        mask3d = _build_mask3d(
            source,
            data3d.shape,
            mask_mode=mask_mode,
            external_mask=mask,
            linefree_ext=linefree_ext,
            linecand_ext=linecand_ext,
            basesup_ext=basesup_ext,
            final_mask_ext=final_mask_ext,
        )
        data3d = _apply_mask_fill(data3d, mask3d, zero_fill=zero_fill, nan_fill=nan_fill, fill_value=fill_value)
    elif mode == "final_moment":
        moment_kind = "final_signal"
        mask_mode = "mask3d" if mask_mode is None else mask_mode
        mask3d = _build_mask3d(
            source,
            data3d.shape,
            mask_mode=mask_mode,
            external_mask=mask,
            linefree_ext=linefree_ext,
            linecand_ext=linecand_ext,
            basesup_ext=basesup_ext,
            final_mask_ext=final_mask_ext,
        )
        data3d = _apply_mask_fill(data3d, mask3d, zero_fill=zero_fill, nan_fill=nan_fill, fill_value=fill_value)
    elif mask_mode is not None:
        mask3d = _build_mask3d(
            source,
            data3d.shape,
            mask_mode=mask_mode,
            external_mask=mask,
            linefree_ext=linefree_ext,
            linecand_ext=linecand_ext,
            basesup_ext=basesup_ext,
            final_mask_ext=final_mask_ext,
        )
        data3d = _apply_mask_fill(data3d, mask3d, zero_fill=zero_fill, nan_fill=nan_fill, fill_value=fill_value)

    nchan = data3d.shape[0]
    chan_idx = _resolve_channel_indices(nchan, spec_axis, chan_range=chan_range, vel_range=vel_range)

    if mode in {"moment0", "provisional_moment", "final_moment"}:
        map_data, out_unit = _moment0_from_subset(
            data3d,
            spec_axis,
            spec_unit,
            chan_idx=chan_idx,
            data_unit=data_unit,
            nan_fill=nan_fill,
        )
    elif mode == "channel_sum":
        sub = data3d[chan_idx, :, :]
        map_data = np.nansum(sub, axis=0)
        if nan_fill:
            all_nan = np.all(~np.isfinite(sub), axis=0)
            map_data[all_nan] = np.nan
        out_unit = data_unit
    elif mode == "channel_mean":
        sub = data3d[chan_idx, :, :]
        with np.errstate(invalid="ignore"):
            map_data = np.nanmean(sub, axis=0)
        if nan_fill:
            all_nan = np.all(~np.isfinite(sub), axis=0)
            map_data[all_nan] = np.nan
        out_unit = data_unit
    elif mode == "channel_slice":
        idx = int(chan_idx[0])
        map_data = np.asarray(data3d[idx, :, :], dtype=float)
        out_unit = data_unit
    elif mode == "velocity_slice":
        if vel_range is None:
            raise ValueError("velocity_slice requires vel_range to be given.")
        vel_target = float(np.mean(vel_range))
        idx = int(np.argmin(np.abs(spec_axis - vel_target)))
        map_data = np.asarray(data3d[idx, :, :], dtype=float)
        out_unit = data_unit
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    hdr2d = _drop_nonspatial_wcs_keywords(cube_input.header)
    wcs2d = WCS(hdr2d).celestial

    orig = orig_hpbw_arcsec
    orig_source = "explicit"
    if orig is None:
        orig, orig_source = _infer_orig_hpbw_arcsec(cube_input.header)
    if smooth_fwhm_arcsec is not None or target_hpbw_arcsec is not None:
        map_data = apply_gaussian_smoothing_2d(
            map_data,
            wcs2d,
            smooth_fwhm_arcsec=smooth_fwhm_arcsec,
            target_hpbw_arcsec=target_hpbw_arcsec,
            orig_hpbw_arcsec=orig,
        )

    eff_hpbw, eff_source = _effective_hpbw_arcsec(
        smooth_fwhm_arcsec=smooth_fwhm_arcsec,
        target_hpbw_arcsec=target_hpbw_arcsec,
        orig_hpbw_arcsec=orig,
    )
    if eff_hpbw is not None:
        hdr2d["BMAJ"] = float(eff_hpbw) / 3600.0
        hdr2d["BMIN"] = float(eff_hpbw) / 3600.0

    if out_unit is not None:
        hdr2d["BUNIT"] = out_unit.to_string()

    meta = dict(cube_input.meta)
    meta.update(
        {
            "mode": mode,
            "moment_kind": moment_kind,
            "range_kind": (
                "channel" if chan_range is not None else
                "velocity" if vel_range is not None else
                "all"
            ),
            "range_value": chan_range if chan_range is not None else vel_range,
            "mask_mode": mask_mode,
            "smooth_info": {
                "smooth_fwhm_arcsec": smooth_fwhm_arcsec,
                "target_hpbw_arcsec": target_hpbw_arcsec,
                "orig_hpbw_arcsec": orig,
                "orig_hpbw_source": orig_source,
                "effective_hpbw_arcsec": eff_hpbw,
                "effective_hpbw_source": eff_source,
            },
            "spectral_axis_unit": str(spec_unit),
            "spectral_axis_mode": spectral_axis_mode,
            "spectral_axis_warning": spectral_axis_warning,
        }
    )
    return Map2D(data=np.asarray(map_data, dtype=float), header=hdr2d, wcs=wcs2d, unit=out_unit, meta=meta)



def make_provisional_moment(
    source: SourceLike,
    *,
    ext: Optional[Union[int, str]] = None,
    linefree_ext: str = "LINEFREE",
    basesup_ext: str = "BASESUP3D",
    linecand_ext: str = "LINECAND3D",
    prefer: str = "auto",
    vel_range: Optional[Tuple[float, float]] = None,
    chan_range: Optional[Tuple[int, int]] = None,
    spectral_unit: Union[str, u.UnitBase] = "km/s",
    zero_fill: bool = False,
    nan_fill: bool = True,
    smooth_fwhm_arcsec: Optional[float] = None,
    target_hpbw_arcsec: Optional[float] = None,
    orig_hpbw_arcsec: Optional[float] = None,
) -> Map2D:
    """Create a provisional moment map from baseline-related information.

    Meaning of prefer='auto'
    ------------------------
    The function looks for associated baseline/signal helper extensions in this order:

    1. LINECAND3D
       A 3D boolean mask where True means "this voxel is a likely signal/line candidate".
       If present, integrate only where LINECAND3D is True.

    2. BASESUP3D
       A 3D boolean mask where True means "this voxel is suitable for baseline support"
       (i.e. likely baseline / line-free). If present and LINECAND3D is absent, integrate the
       complement of BASESUP3D.

    3. LINEFREE
       Usually a 1D channel mask where True means "this channel is line-free". If neither
       LINECAND3D nor BASESUP3D is present, integrate the complement of LINEFREE, namely the
       channels that are *not* marked line-free.

    In short, auto means: use the most direct signal mask if available; otherwise infer the
    signal region by taking the complement of a baseline mask.
    """
    prefer = str(prefer).lower()
    if prefer == "auto":
        mask_mode = None
    else:
        if prefer not in {"linecand3d", "basesup3d", "linefree_complement"}:
            raise ValueError("prefer must be one of: auto, linecand3d, basesup3d, linefree_complement")
        mask_mode = prefer

    return make_2d_map(
        source,
        ext=ext,
        mode="provisional_moment",
        chan_range=chan_range,
        vel_range=vel_range,
        spectral_unit=spectral_unit,
        mask_mode=mask_mode,
        zero_fill=zero_fill,
        nan_fill=nan_fill,
        smooth_fwhm_arcsec=smooth_fwhm_arcsec,
        target_hpbw_arcsec=target_hpbw_arcsec,
        orig_hpbw_arcsec=orig_hpbw_arcsec,
        linefree_ext=linefree_ext,
        linecand_ext=linecand_ext,
        basesup_ext=basesup_ext,
    )



def make_final_moment(
    source: SourceLike,
    *,
    ext: Optional[Union[int, str]] = None,
    final_mask_ext: str = "MASK3D",
    vel_range: Optional[Tuple[float, float]] = None,
    chan_range: Optional[Tuple[int, int]] = None,
    spectral_unit: Union[str, u.UnitBase] = "km/s",
    zero_fill: bool = False,
    nan_fill: bool = True,
    smooth_fwhm_arcsec: Optional[float] = None,
    target_hpbw_arcsec: Optional[float] = None,
    orig_hpbw_arcsec: Optional[float] = None,
) -> Map2D:
    """Create a final signal moment map using MASK3D."""
    return make_2d_map(
        source,
        ext=ext,
        mode="final_moment",
        chan_range=chan_range,
        vel_range=vel_range,
        spectral_unit=spectral_unit,
        mask_mode="mask3d",
        zero_fill=zero_fill,
        nan_fill=nan_fill,
        smooth_fwhm_arcsec=smooth_fwhm_arcsec,
        target_hpbw_arcsec=target_hpbw_arcsec,
        orig_hpbw_arcsec=orig_hpbw_arcsec,
        final_mask_ext=final_mask_ext,
    )


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------


def _looks_like_layer_mapping(obj: Any) -> bool:
    return isinstance(obj, Mapping)


def _layer_dict(obj: Any, *, default_kind: Optional[str] = None) -> Dict[str, Any]:
    if isinstance(obj, Mapping):
        out = dict(obj)
    else:
        out = {"source": obj}
    if default_kind is not None and "kind" not in out:
        out["kind"] = default_kind
    return out


def _layer_to_map2d(layer: Any, *, default_mode: str = "identity") -> Map2D:
    layer = _layer_dict(layer)
    contour_map_keys = {
        "mode",
        "chan_range",
        "vel_range",
        "spectral_unit",
        "mask",
        "mask_mode",
        "zero_fill",
        "nan_fill",
        "smooth_fwhm_arcsec",
        "target_hpbw_arcsec",
        "orig_hpbw_arcsec",
        "linefree_ext",
        "linecand_ext",
        "basesup_ext",
        "final_mask_ext",
    }
    mode = layer.get("mode", default_mode)
    if any(k in layer for k in contour_map_keys) or mode not in {None, "identity", "map", "2d"}:
        return make_2d_map(
            source=layer.get("source"),
            ext=layer.get("ext"),
            mode=mode or "identity",
            chan_range=layer.get("chan_range"),
            vel_range=layer.get("vel_range"),
            spectral_unit=layer.get("spectral_unit", "km/s"),
            mask=layer.get("mask"),
            mask_mode=layer.get("mask_mode"),
            zero_fill=layer.get("zero_fill", False),
            nan_fill=layer.get("nan_fill", True),
            smooth_fwhm_arcsec=layer.get("smooth_fwhm_arcsec"),
            target_hpbw_arcsec=layer.get("target_hpbw_arcsec"),
            orig_hpbw_arcsec=layer.get("orig_hpbw_arcsec"),
            linefree_ext=layer.get("linefree_ext", "LINEFREE"),
            linecand_ext=layer.get("linecand_ext", "LINECAND3D"),
            basesup_ext=layer.get("basesup_ext", "BASESUP3D"),
            final_mask_ext=layer.get("final_mask_ext", "MASK3D"),
        )
    return resolve_map_input(
        source=layer.get("source"),
        data=layer.get("data"),
        header=layer.get("header"),
        ext=layer.get("ext"),
    )


def _resolve_base_spec(base: Any) -> Dict[str, Any]:
    layer = _layer_dict(base, default_kind="image")
    kind = str(layer.get("kind", "image")).lower()
    if kind not in {"image", "rgb"}:
        raise ValueError("base kind must be 'image' or 'rgb'.")
    return layer


def _save_figure(
    fig: Any,
    save: Optional[Union[str, Path]],
    *,
    dpi: int = 250,
    bbox_inches: str = "tight",
    transparent: bool = False,
) -> None:
    if save is None:
        return
    fig.savefig(str(save), dpi=dpi, bbox_inches=bbox_inches, transparent=transparent)


def _format_readout_value(value: Any) -> str:
    try:
        v = float(value)
    except Exception:
        return str(value)
    if not np.isfinite(v):
        return "nan"
    av = abs(v)
    if av == 0:
        return "0"
    if av >= 1.0e6 or av < 1.0e-6:
        return f"{v:.5g}"
    return f"{v:.5g}"



def _estimate_readout_coord_precision(map2d: Map2D) -> Tuple[int, int]:
    try:
        scales = np.asarray(proj_plane_pixel_scales(map2d.wcs.celestial), dtype=float)
        finite = scales[np.isfinite(scales) & (scales > 0)]
        if finite.size == 0:
            raise ValueError
        deg_per_pix = float(np.nanmin(np.abs(finite)))
        arcsec_per_pix = deg_per_pix * 3600.0
    except Exception:
        arcsec_per_pix = 1.0
        deg_per_pix = 1.0 / 3600.0

    if arcsec_per_pix >= 1.0:
        sexa_precision = 0
    elif arcsec_per_pix >= 0.1:
        sexa_precision = 1
    else:
        sexa_precision = 2

    decimal_deg_precision = int(np.clip(np.ceil(-np.log10(max(deg_per_pix, 1.0e-12))) + 2, 4, 8))
    return sexa_precision, decimal_deg_precision



def _format_ra_deg_to_hms(ra_deg: float, precision: int = 0) -> str:
    ang = Angle(float(ra_deg), unit=u.deg).wrap_at(360.0 * u.deg)
    return ang.to_string(unit=u.hourangle, sep=":", pad=True, precision=precision)



def _format_dec_deg_to_dms(dec_deg: float, precision: int = 0) -> str:
    ang = Angle(float(dec_deg), unit=u.deg)
    return ang.to_string(unit=u.deg, sep=":", pad=True, alwayssign=True, precision=precision)



def _format_lonlat_readout(map2d: Map2D, lon_deg: float, lat_deg: float) -> Tuple[str, str]:
    native = _infer_native_coord_system(map2d.wcs)
    sexa_precision, decimal_deg_precision = _estimate_readout_coord_precision(map2d)

    if native == "equatorial":
        return (
            f"ra={_format_ra_deg_to_hms(lon_deg, precision=sexa_precision)}",
            f"dec={_format_dec_deg_to_dms(lat_deg, precision=sexa_precision)}",
        )

    lon_s = f"{float(lon_deg):.{decimal_deg_precision}f}"
    lat_s = f"{float(lat_deg):.{decimal_deg_precision}f}"
    return (f"lon={lon_s} deg", f"lat={lat_s} deg")



def _disable_artist_cursor_data(artist: Any) -> None:
    try:
        artist.get_cursor_data = lambda event: None
    except Exception:
        pass
    try:
        artist.format_cursor_data = lambda data: ""
    except Exception:
        pass



def _disable_axes_cursor_readout(ax: Any) -> None:
    try:
        ax.format_coord = lambda x, y: ""
    except Exception:
        pass
    try:
        for child in ax.get_children():
            _disable_artist_cursor_data(child)
    except Exception:
        pass



def _apply_readout_formatter(ax: Any, map2d: Map2D) -> None:
    def _fmt(x: float, y: float) -> str:
        try:
            ix = int(np.round(x))
            iy = int(np.round(y))
            if 0 <= ix < map2d.data.shape[1] and 0 <= iy < map2d.data.shape[0]:
                value = map2d.data[iy, ix]
            else:
                value = np.nan
            world = map2d.wcs.pixel_to_world_values(x, y)
            lon = float(world[0])
            lat = float(world[1])
            lon_s, lat_s = _format_lonlat_readout(map2d, lon, lat)
            value_s = _format_readout_value(value)
            return f"x={x:.2f}, y={y:.2f} | {lon_s}, {lat_s} | value={value_s}"
        except Exception:
            return f"x={x:.2f}, y={y:.2f}"

    ax.format_coord = _fmt


def _beam_config_from_map(map2d: Map2D, beam: Union[str, Mapping[str, Any]]) -> Dict[str, Any]:
    if isinstance(beam, str) and beam.lower() in {"header", "auto"}:
        smooth_info = map2d.meta.get("smooth_info", {}) if isinstance(map2d.meta, dict) else {}
        eff_hpbw = smooth_info.get("effective_hpbw_arcsec")
        if eff_hpbw is not None:
            return {"major": float(eff_hpbw), "minor": float(eff_hpbw), "angle": 0.0, "units": "arcsec"}
        bmaj = map2d.header.get("BMAJ")
        bmin = map2d.header.get("BMIN", bmaj)
        bpa = map2d.header.get("BPA", 0.0)
        if bmaj is not None:
            return {"major": float(bmaj), "minor": float(bmin), "angle": float(bpa), "units": "deg"}
        return {}
    return dict(beam)


def _add_beam_artist(ax: Any, map2d: Map2D, beam: Union[str, Mapping[str, Any]]) -> Optional[Any]:
    from matplotlib.patches import Ellipse

    beam_cfg = _beam_config_from_map(map2d, beam)
    if not beam_cfg:
        return None

    pix_arcsec = _pixel_scale_arcsec(map2d.wcs)
    major = beam_cfg.get("major")
    minor = beam_cfg.get("minor", major)
    angle = beam_cfg.get("angle", 0.0)
    units = str(beam_cfg.get("units", "arcsec")).lower()
    if major is None:
        return None

    if units in {"arcsec", "asec", "arcseconds"}:
        major_pix = float(major) / pix_arcsec
        minor_pix = float(minor) / pix_arcsec
    elif units in {"deg", "degree", "degrees"}:
        major_pix = float(major) * 3600.0 / pix_arcsec
        minor_pix = float(minor) * 3600.0 / pix_arcsec
    elif units in {"pix", "pixel", "pixels"}:
        major_pix = float(major)
        minor_pix = float(minor)
    else:
        raise ValueError(f"Unsupported beam units: {units}")

    ny, nx = map2d.data.shape
    x0 = beam_cfg.get("x", 0.12 * nx)
    y0 = beam_cfg.get("y", 0.12 * ny)
    artist = Ellipse(
        (x0, y0),
        width=major_pix,
        height=minor_pix,
        angle=angle,
        facecolor=beam_cfg.get("facecolor", "none"),
        edgecolor=beam_cfg.get("edgecolor", "white"),
        linewidth=beam_cfg.get("linewidth", 1.2),
        alpha=beam_cfg.get("alpha", 0.9),
    )
    ax.add_patch(artist)
    return artist


def _world_transform_for(ax: Any, wcs: WCS) -> Any:
    """Return the native world-coordinate transform for the current WCSAxes.

    World overlays handled by this module are first converted into the base
    celestial frame of the axes and then passed as native longitude/latitude
    values in degrees. In that situation, WCSAxes expects the generic
    ``'world'`` transform rather than a transform constructed from the base WCS
    object itself. Using ``ax.get_transform('world')`` avoids ambiguities
    between pixel/data and world inputs for marker/text/catalog overlays.
    """
    return ax.get_transform("world")


def _normalize_frame_name(frame: Optional[Any]) -> Optional[str]:
    if frame is None:
        return None
    text = str(frame).strip().lower()
    if text in {"", "none", "auto", "world", "native"}:
        return None
    aliases = {
        "equatorial": "icrs",
        "j2000": "icrs",
        "radec": "icrs",
        "ra/dec": "icrs",
        "fk5j2000": "fk5",
        "gal": "galactic",
        "gal_lonlat": "galactic",
        "glon/glat": "galactic",
        "l/b": "galactic",
    }
    return aliases.get(text, text)


def _target_frame_name_from_wcs(wcs: WCS) -> str:
    mode = _infer_native_coord_system(wcs)
    if mode == "galactic":
        return "galactic"
    cwcs = wcs.celestial.wcs
    radesys = _normalize_frame_name(getattr(cwcs, "radesys", None))
    if radesys in {"icrs", "fk5", "fk4"}:
        return str(radesys)
    return "icrs"


def _transform_skycoord_to_frame(sc: SkyCoord, frame_name: Optional[str]) -> SkyCoord:
    target = _normalize_frame_name(frame_name)
    if target is None:
        return sc
    if target == "icrs":
        return sc.icrs
    if target == "fk5":
        return sc.fk5
    if target == "fk4":
        return sc.fk4
    if target == "galactic":
        return sc.galactic
    try:
        return sc.transform_to(target)
    except Exception as exc:
        raise ValueError(f"Unsupported or unavailable coordinate frame: {frame_name}") from exc


def _skycoord_from_lonlat(
    lon: np.ndarray,
    lat: np.ndarray,
    *,
    frame_name: Optional[str],
    coord_unit: Union[str, u.UnitBase] = "deg",
) -> SkyCoord:
    frame = _normalize_frame_name(frame_name) or "icrs"
    unit = u.Unit(coord_unit)
    lon_q = np.asarray(lon, dtype=float) * unit
    lat_q = np.asarray(lat, dtype=float) * unit
    if frame == "galactic":
        return SkyCoord(l=lon_q, b=lat_q, frame="galactic")
    if frame in {"icrs", "fk5", "fk4"}:
        return SkyCoord(ra=lon_q, dec=lat_q, frame=frame)
    return SkyCoord(lon_q, lat_q, frame=frame)


def _skycoord_to_base_lonlat(sc: SkyCoord, base_wcs: WCS) -> np.ndarray:
    target = _target_frame_name_from_wcs(base_wcs)
    sc_t = _transform_skycoord_to_frame(sc, target)
    return np.column_stack([
        np.asarray(sc_t.spherical.lon.deg, dtype=float).ravel(),
        np.asarray(sc_t.spherical.lat.deg, dtype=float).ravel(),
    ])


def _world_array_to_base_lonlat(
    arr: np.ndarray,
    *,
    base_wcs: WCS,
    frame_name: Optional[str] = None,
    coord_unit: Union[str, u.UnitBase] = "deg",
) -> np.ndarray:
    target = _target_frame_name_from_wcs(base_wcs)
    source = _normalize_frame_name(frame_name) or target
    sc = _skycoord_from_lonlat(
        np.asarray(arr, dtype=float)[:, 0],
        np.asarray(arr, dtype=float)[:, 1],
        frame_name=source,
        coord_unit=coord_unit,
    )
    return _skycoord_to_base_lonlat(sc, base_wcs)


def _infer_frame_from_position_column_names(lon_name: Optional[str], lat_name: Optional[str]) -> Optional[str]:
    lon_lower = "" if lon_name is None else str(lon_name).strip().lower()
    lat_lower = "" if lat_name is None else str(lat_name).strip().lower()
    gal_lon = {"glon", "glon_deg", "l", "l_deg"}
    gal_lat = {"glat", "glat_deg", "b", "b_deg"}
    eq_lon = {"ra", "ra_deg"}
    eq_lat = {"dec", "dec_deg"}
    if lon_lower in gal_lon or lat_lower in gal_lat:
        return "galactic"
    if lon_lower in eq_lon or lat_lower in eq_lat:
        return "icrs"
    return None


def _layer_plot_positions(layer: Mapping[str, Any], base_wcs: WCS) -> Tuple[np.ndarray, bool]:
    positions, mode, frame_name = _extract_marker_positions(layer)
    if mode == "skycoord":
        return _skycoord_to_base_lonlat(positions, base_wcs), True
    if mode == "world":
        return _world_array_to_base_lonlat(
            np.asarray(positions, dtype=float),
            base_wcs=base_wcs,
            frame_name=frame_name,
            coord_unit=layer.get("coord_unit", "deg"),
        ), True
    return np.asarray(positions, dtype=float), False
def _table_like_column_names(obj: Any) -> List[str]:
    if obj is None:
        return []
    try:
        if hasattr(obj, "colnames"):
            return [str(v) for v in list(obj.colnames)]
    except Exception:
        pass
    try:
        if hasattr(obj, "columns"):
            cols = getattr(obj, "columns")
            if hasattr(cols, "keys"):
                return [str(v) for v in list(cols.keys())]
            return [str(v) for v in list(cols)]
    except Exception:
        pass
    try:
        if isinstance(obj, Mapping) or hasattr(obj, "keys"):
            return [str(v) for v in list(obj.keys())]
    except Exception:
        pass
    return []



def _table_like_get_column(obj: Any, name: str) -> Optional[np.ndarray]:
    try:
        col = obj[name]
    except Exception:
        return None
    try:
        if hasattr(col, "to_numpy"):
            col = col.to_numpy()
    except Exception:
        pass
    try:
        return np.asarray(col)
    except Exception:
        try:
            return np.asarray(list(col))
        except Exception:
            return None



def _find_table_like_column(obj: Any, candidates: Sequence[str]) -> Tuple[Optional[str], Optional[np.ndarray]]:
    names = _table_like_column_names(obj)
    if not names:
        return None, None

    name_map = {str(n).lower(): str(n) for n in names}
    for cand in candidates:
        key = str(cand)
        if key in names:
            return key, _table_like_get_column(obj, key)
        lowered = key.lower()
        if lowered in name_map:
            actual = name_map[lowered]
            return actual, _table_like_get_column(obj, actual)
    return None, None



def _extract_positions_from_table_like(layer: Mapping[str, Any]) -> Tuple[np.ndarray, Optional[str], Optional[str]]:
    table = layer.get("table", layer.get("catalog"))
    if table is None:
        raise ValueError("No table-like input was supplied.")

    mode = str(layer.get("coord_mode", "auto")).lower()

    lon_candidates = [v for v in (layer.get("lon_col"),) if v is not None] + [
        "lon", "glon", "l", "longitude", "ra", "ra_deg", "lon_deg", "glon_deg", "l_deg"
    ]
    lat_candidates = [v for v in (layer.get("lat_col"),) if v is not None] + [
        "lat", "glat", "b", "latitude", "dec", "dec_deg", "lat_deg", "glat_deg", "b_deg"
    ]
    x_candidates = [v for v in (layer.get("x_col"),) if v is not None] + [
        "x", "xpix", "x_pix", "xpixel", "col", "ix"
    ]
    y_candidates = [v for v in (layer.get("y_col"),) if v is not None] + [
        "y", "ypix", "y_pix", "ypixel", "row", "iy"
    ]

    explicit_frame = _normalize_frame_name(layer.get("frame"))

    if mode in {"auto", "world"}:
        lon_name, lon = _find_table_like_column(table, lon_candidates)
        lat_name, lat = _find_table_like_column(table, lat_candidates)
        if lon is not None and lat is not None:
            lon = np.asarray(lon, dtype=float).ravel()
            lat = np.asarray(lat, dtype=float).ravel()
            if lon.size != lat.size:
                raise ValueError("table-like lon/lat columns must have the same length.")
            inferred_frame = _infer_frame_from_position_column_names(lon_name, lat_name)
            return np.column_stack([lon, lat]), "world", explicit_frame or inferred_frame
        if mode == "world":
            raise ValueError("coord_mode='world' was requested, but lon/lat columns were not found.")

    _, x = _find_table_like_column(table, x_candidates)
    _, y = _find_table_like_column(table, y_candidates)
    if x is not None and y is not None:
        x = np.asarray(x, dtype=float).ravel()
        y = np.asarray(y, dtype=float).ravel()
        if x.size != y.size:
            raise ValueError("table-like x/y columns must have the same length.")
        return np.column_stack([x, y]), None, None

    available = ", ".join(_table_like_column_names(table))
    raise ValueError(
        "Could not infer positions from the table-like input. "
        f"Available columns: {available if available else '(none)'}"
    )



def _extract_catalog_labels(layer: Mapping[str, Any], npos: Optional[int] = None) -> Optional[List[str]]:
    labels = layer.get("labels")
    if labels is None:
        table = layer.get("table", layer.get("catalog"))
        if table is not None:
            label_candidates = [v for v in (layer.get("label_col"),) if v is not None] + [
                "label", "labels", "name", "names", "id", "ids", "source", "object", "obj"
            ]
            _, col = _find_table_like_column(table, label_candidates)
            if col is not None:
                labels = [str(v) for v in np.asarray(col).ravel().tolist()]
    if labels is None:
        return None
    out = [str(v) if v is not None else "" for v in list(labels)]
    if npos is not None and len(out) != int(npos):
        raise ValueError("catalog labels must have the same length as the positions.")
    return out


def _extract_marker_positions(layer: Mapping[str, Any]) -> Tuple[Any, Optional[str], Optional[str]]:
    if layer.get("table", layer.get("catalog")) is not None:
        return _extract_positions_from_table_like(layer)

    if layer.get("coords") is not None:
        arr = np.asarray(layer.get("coords"), dtype=float)
        if arr.ndim == 1:
            arr = arr[None, :]
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("marker coords must have shape (N, 2).")
        return arr, "world", _normalize_frame_name(layer.get("frame"))

    if layer.get("skycoord") is not None:
        sc = layer.get("skycoord")
        if not isinstance(sc, SkyCoord):
            raise TypeError("skycoord must be an astropy.coordinates.SkyCoord instance.")
        return sc, "skycoord", None

    if layer.get("xy") is not None:
        arr = np.asarray(layer.get("xy"), dtype=float)
        if arr.ndim == 1:
            arr = arr[None, :]
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("marker xy must have shape (N, 2).")
        return arr, None, None

    if layer.get("lon") is not None and layer.get("lat") is not None:
        lon = np.asarray(layer.get("lon"), dtype=float).ravel()
        lat = np.asarray(layer.get("lat"), dtype=float).ravel()
        if lon.size != lat.size:
            raise ValueError("lon and lat must have the same number of elements.")
        return np.column_stack([lon, lat]), "world", _normalize_frame_name(layer.get("frame"))

    raise ValueError("marker/catalog layer requires coords=..., skycoord=..., xy=..., or lon=... and lat=... .")
def _set_artist_label(artist: Any, label: Optional[str]) -> Any:
    if label in (None, ""):
        return artist
    try:
        artist.set_label(str(label))
        return artist
    except Exception:
        pass
    try:
        artist.collections[0].set_label(str(label))
    except Exception:
        pass
    return artist



def _annotation_map_from_base(base_obj: Any) -> Optional[Map2D]:
    if isinstance(base_obj, Map2D):
        return base_obj
    if isinstance(base_obj, RGBMap):
        return Map2D(
            data=np.asarray(np.nanmean(base_obj.rgb, axis=2), dtype=float),
            header=_copy_header(base_obj.header),
            wcs=base_obj.wcs,
            unit=None,
            meta=dict(base_obj.meta),
        )
    return None



def _nice_length(value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        return 1.0
    exp = float(np.floor(np.log10(value)))
    base = 10.0 ** exp
    best = base
    for mult in (1.0, 2.0, 5.0, 10.0):
        cand = mult * base
        if cand <= value:
            best = cand
        else:
            break
    return float(best)



def _format_angular_label_arcsec(value_arcsec: float) -> str:
    value_arcsec = float(value_arcsec)
    if value_arcsec >= 3600.0 and abs(value_arcsec / 3600.0 - round(value_arcsec / 3600.0)) < 1e-6:
        return f"{int(round(value_arcsec / 3600.0))}°"
    if value_arcsec >= 60.0 and abs(value_arcsec / 60.0 - round(value_arcsec / 60.0)) < 1e-6:
        return f"{int(round(value_arcsec / 60.0))}′"
    if value_arcsec >= 60.0:
        return f"{value_arcsec / 60.0:.2g}′"
    if abs(value_arcsec - round(value_arcsec)) < 1e-6:
        return f"{int(round(value_arcsec))}″"
    return f"{value_arcsec:.2g}″"



def _resolve_anchor_xy(shape: Tuple[int, int], *, width_pix: float, height_pix: float = 0.0, location: str = "lower right", margin_frac: float = 0.06) -> Tuple[float, float]:
    ny, nx = shape
    mx = margin_frac * nx
    my = margin_frac * ny
    loc = str(location).lower().replace("_", " ")
    if "left" in loc:
        x0 = mx
    else:
        x0 = nx - mx - width_pix
    if "upper" in loc:
        y0 = ny - my - max(height_pix, 1.0)
    else:
        y0 = my + max(height_pix, 1.0)
    return float(x0), float(y0)



def add_scalebar(ax: Any, map2d: Map2D, config: Optional[Union[bool, Mapping[str, Any]]] = None) -> Dict[str, Any]:
    cfg = {} if config in (None, True) else dict(config)
    pix_arcsec = _pixel_scale_arcsec(map2d.wcs)
    if not np.isfinite(pix_arcsec) or pix_arcsec <= 0:
        raise ValueError("Could not determine a finite positive pixel scale for the scalebar.")

    units = str(cfg.get("units", "arcsec")).lower()
    length = cfg.get("length")
    if length is None:
        target_arcsec = 0.18 * map2d.data.shape[1] * pix_arcsec
        length_arcsec = _nice_length(target_arcsec)
    else:
        if units in {"arcsec", "asec", "arcseconds"}:
            length_arcsec = float(length)
        elif units in {"arcmin", "amin", "arcminutes"}:
            length_arcsec = float(length) * 60.0
        elif units in {"deg", "degree", "degrees"}:
            length_arcsec = float(length) * 3600.0
        elif units in {"pix", "pixel", "pixels"}:
            length_arcsec = float(length) * pix_arcsec
        else:
            raise ValueError(f"Unsupported scalebar units: {units}")

    width_pix = float(length_arcsec / pix_arcsec)
    x0, y0 = _resolve_anchor_xy(
        map2d.data.shape,
        width_pix=width_pix,
        height_pix=0.10 * map2d.data.shape[0],
        location=cfg.get("location", "lower right"),
        margin_frac=float(cfg.get("margin_frac", 0.06)),
    )
    x1 = x0 + width_pix

    line_kwargs = dict(cfg.get("line_style", {}))
    line_kwargs.setdefault("color", cfg.get("color", "white"))
    line_kwargs.setdefault("linewidth", cfg.get("linewidth", 2.0))
    line_kwargs.setdefault("solid_capstyle", cfg.get("solid_capstyle", "butt"))
    line = ax.plot([x0, x1], [y0, y0], **line_kwargs)[0]

    label = cfg.get("label")
    if label in (None, "auto"):
        label = _format_angular_label_arcsec(length_arcsec)
    text_kwargs = dict(cfg.get("text_style", {}))
    text_kwargs.setdefault("color", cfg.get("color", "white"))
    text_kwargs.setdefault("fontsize", cfg.get("fontsize", 9))
    text_kwargs.setdefault("ha", "center")
    text_kwargs.setdefault("va", cfg.get("text_va", "bottom"))
    text_dy = float(cfg.get("text_dy", 0.015 * map2d.data.shape[0]))
    text = ax.text(0.5 * (x0 + x1), y0 + text_dy, str(label), **text_kwargs)
    return {"line": line, "text": text, "length_arcsec": float(length_arcsec)}



def add_north_arrow(ax: Any, map2d: Map2D, config: Optional[Union[bool, Mapping[str, Any]]] = None) -> Dict[str, Any]:
    cfg = {} if config in (None, True) else dict(config)
    pix_arcsec = _pixel_scale_arcsec(map2d.wcs)
    length = cfg.get("length")
    units = str(cfg.get("units", "arcsec")).lower()
    if length is None:
        length_arcsec = _nice_length(0.12 * min(map2d.data.shape) * pix_arcsec)
    else:
        if units in {"arcsec", "asec", "arcseconds"}:
            length_arcsec = float(length)
        elif units in {"arcmin", "amin", "arcminutes"}:
            length_arcsec = float(length) * 60.0
        elif units in {"deg", "degree", "degrees"}:
            length_arcsec = float(length) * 3600.0
        elif units in {"pix", "pixel", "pixels"}:
            length_arcsec = float(length) * pix_arcsec
        else:
            raise ValueError(f"Unsupported north_arrow units: {units}")

    width_pix = float(length_arcsec / pix_arcsec)
    x0, y0 = _resolve_anchor_xy(
        map2d.data.shape,
        width_pix=width_pix,
        height_pix=width_pix,
        location=cfg.get("location", "upper left"),
        margin_frac=float(cfg.get("margin_frac", 0.06)),
    )

    world0 = map2d.wcs.pixel_to_world(float(x0), float(y0))
    if not isinstance(world0, SkyCoord):
        raise TypeError("north_arrow requires a celestial WCS that converts pixels to SkyCoord.")

    sep = float(length_arcsec) * u.arcsec
    wn = world0.directional_offset_by(0.0 * u.deg, sep)
    we = world0.directional_offset_by(90.0 * u.deg, sep)
    xn, yn = map2d.wcs.world_to_pixel(wn)
    xe, ye = map2d.wcs.world_to_pixel(we)

    arrowprops = dict(cfg.get("arrowprops", {}))
    arrowprops.setdefault("color", cfg.get("color", "white"))
    arrowprops.setdefault("lw", cfg.get("linewidth", 1.5))
    arrowprops.setdefault("arrowstyle", cfg.get("arrowstyle", "-|>"))

    ann_n = ax.annotate("", xy=(xn, yn), xytext=(x0, y0), arrowprops=arrowprops)
    ann_e = ax.annotate("", xy=(xe, ye), xytext=(x0, y0), arrowprops=arrowprops)

    text_kwargs = dict(cfg.get("text_style", {}))
    text_kwargs.setdefault("color", cfg.get("color", "white"))
    text_kwargs.setdefault("fontsize", cfg.get("fontsize", 9))
    txt_n = ax.text(float(xn), float(yn), str(cfg.get("north_label", "N")), ha="center", va="bottom", **text_kwargs)
    txt_e = ax.text(float(xe), float(ye), str(cfg.get("east_label", "E")), ha="left", va="center", **text_kwargs)
    return {"north_arrow": ann_n, "east_arrow": ann_e, "north_text": txt_n, "east_text": txt_e, "length_arcsec": float(length_arcsec)}



def _scatter_points(ax: Any, base_wcs: WCS, layer: Mapping[str, Any]) -> Any:
    kwargs = dict(layer.get("style", {}))
    kwargs.setdefault("s", layer.get("s", 40))
    kwargs.setdefault("marker", layer.get("marker", "+"))
    kwargs.setdefault("color", layer.get("color", "white"))
    kwargs.setdefault("linewidths", layer.get("linewidths", 1.0))

    arr, is_world = _layer_plot_positions(layer, base_wcs)
    if is_world:
        artist = ax.scatter(arr[:, 0], arr[:, 1], transform=_world_transform_for(ax, base_wcs), **kwargs)
    else:
        artist = ax.scatter(arr[:, 0], arr[:, 1], **kwargs)
    return _set_artist_label(artist, layer.get("label"))



def _scatter_catalog(ax: Any, base_wcs: WCS, layer: Mapping[str, Any]) -> Dict[str, Any]:
    marker_artist = _scatter_points(ax, base_wcs, layer)
    arr, is_world = _layer_plot_positions(layer, base_wcs)
    labels = _extract_catalog_labels(layer, npos=arr.shape[0])
    text_artists: List[Any] = []
    if labels is not None:
        text_kwargs = dict(layer.get("text_style", {}))
        text_kwargs.setdefault("color", layer.get("text_color", layer.get("color", "white")))
        text_kwargs.setdefault("fontsize", layer.get("fontsize", 9))
        text_kwargs.setdefault("ha", layer.get("ha", "left"))
        text_kwargs.setdefault("va", layer.get("va", "bottom"))
        dx = float(layer.get("label_dx", 6.0))
        dy = float(layer.get("label_dy", 4.0))
        offset_mode = str(layer.get("label_offset_mode", "pixel")).lower()
        for (x, y), label in zip(arr, labels):
            if label in (None, ""):
                continue
            if offset_mode == "world":
                if is_world:
                    art = ax.text(x + dx, y + dy, str(label), transform=_world_transform_for(ax, base_wcs), **text_kwargs)
                else:
                    art = ax.text(x + dx, y + dy, str(label), **text_kwargs)
            else:
                ann_kwargs = dict(text_kwargs)
                if is_world:
                    art = ax.annotate(
                        str(label),
                        xy=(x, y),
                        xycoords=_world_transform_for(ax, base_wcs),
                        xytext=(dx, dy),
                        textcoords="offset pixels",
                        **ann_kwargs,
                    )
                else:
                    art = ax.annotate(
                        str(label),
                        xy=(x, y),
                        xytext=(dx, dy),
                        textcoords="offset pixels",
                        **ann_kwargs,
                    )
            text_artists.append(art)
    return {"marker": marker_artist, "texts": text_artists}



def _add_text_artist(ax: Any, base_wcs: WCS, layer: Mapping[str, Any]) -> Any:
    kwargs = dict(layer.get("style", {}))
    kwargs.setdefault("color", layer.get("color", "white"))
    kwargs.setdefault("fontsize", layer.get("fontsize", 10))
    kwargs.setdefault("ha", layer.get("ha", "left"))
    kwargs.setdefault("va", layer.get("va", "bottom"))
    text = str(layer.get("text", ""))

    if layer.get("coord") is not None:
        arr = np.asarray(layer.get("coord"), dtype=float).ravel()
        if arr.size != 2:
            raise ValueError("text coord must contain exactly two world-coordinate values.")
        xy = _world_array_to_base_lonlat(
            arr.reshape(1, 2),
            base_wcs=base_wcs,
            frame_name=_normalize_frame_name(layer.get("frame")),
            coord_unit=layer.get("coord_unit", "deg"),
        )[0]
        return ax.text(float(xy[0]), float(xy[1]), text, transform=_world_transform_for(ax, base_wcs), **kwargs)
    if layer.get("skycoord") is not None:
        sc = layer.get("skycoord")
        if not isinstance(sc, SkyCoord):
            raise TypeError("skycoord must be an astropy.coordinates.SkyCoord instance.")
        xy = _skycoord_to_base_lonlat(sc, base_wcs)[0]
        return ax.text(float(xy[0]), float(xy[1]), text, transform=_world_transform_for(ax, base_wcs), **kwargs)
    if layer.get("xy") is not None:
        x, y = layer.get("xy")
        return ax.text(x, y, text, **kwargs)
    raise ValueError("text layer requires coord=..., skycoord=..., or xy=... .")
def _reproject_map2d_to(overlay: Map2D, *, target_wcs: WCS, target_shape: Tuple[int, int]) -> Map2D:
    try:
        from reproject import reproject_interp  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "Image reprojection requested, but the optional 'reproject' package is not available."
        ) from exc

    array, _ = reproject_interp((np.asarray(overlay.data, dtype=float), overlay.wcs), target_wcs, shape_out=target_shape)
    hdr = _copy_header(overlay.header)
    hdr.update(target_wcs.to_header())
    return Map2D(data=np.asarray(array, dtype=float), header=hdr, wcs=target_wcs.celestial, unit=overlay.unit, meta=dict(overlay.meta, reprojected=True))


# Compare only core celestial header keys for a cheap same-grid check.
def _same_celestial_wcs_and_shape(a: Map2D, target_wcs: WCS, target_shape: Tuple[int, int]) -> bool:
    if tuple(a.data.shape) != tuple(target_shape):
        return False
    try:
        ah = a.wcs.celestial.to_header()
        bh = target_wcs.celestial.to_header()
        keys = ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "CUNIT1", "CUNIT2")
        return all(str(ah.get(k)) == str(bh.get(k)) for k in keys)
    except Exception:
        return False


def _add_image_overlay(ax: Any, layer: Mapping[str, Any], *, target_wcs: Optional[WCS] = None, target_shape: Optional[Tuple[int, int]] = None) -> Dict[str, Any]:
    overlay = _layer_to_map2d(layer, default_mode="identity")
    reproject_mode = str(layer.get("reproject", "never")).lower()
    if reproject_mode not in {"never", "auto", "required"}:
        raise ValueError("image overlay reproject must be one of: never, auto, required")

    if target_wcs is not None and target_shape is not None and reproject_mode != "never":
        if not _same_celestial_wcs_and_shape(overlay, target_wcs, target_shape):
            try:
                overlay = _reproject_map2d_to(overlay, target_wcs=target_wcs, target_shape=target_shape)
            except Exception as exc:
                if reproject_mode == "required":
                    raise
                warnings.warn(f"Image overlay reprojection was skipped: {exc}", RuntimeWarning)

    norm = layer.get("norm")
    if norm is None:
        norm = build_normalize(
            overlay.data,
            mode=layer.get("norm_mode", "asinh"),
            percentile=layer.get("norm_percentile", (1.0, 99.5)),
            cmin=layer.get("cmin"),
            cmax=layer.get("cmax"),
            stretch_a=layer.get("stretch_a", 0.1),
            power_gamma=layer.get("power_gamma", 1.0),
        )
    artist = ax.imshow(
        overlay.data,
        origin=layer.get("origin", "lower"),
        cmap=layer.get("cmap", "magma"),
        norm=norm,
        alpha=layer.get("alpha", 0.5),
        interpolation=layer.get("interpolation", "nearest"),
        transform=ax.get_transform(overlay.wcs),
        zorder=layer.get("zorder"),
    )
    _disable_artist_cursor_data(artist)
    _set_artist_label(artist, layer.get("label"))
    return {"artist": artist, "map2d": overlay}


def _prepare_contour_layer(layer: Mapping[str, Any]) -> Dict[str, Any]:
    contour = dict(layer)
    contour.pop("kind", None)
    return contour



# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------


def _auto_contour_levels(data: np.ndarray) -> np.ndarray:
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return np.array([0.0])
    vmax = float(np.nanmax(finite))
    if vmax > 0:
        return vmax * np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    vmin = float(np.nanmin(finite))
    return np.linspace(vmin, vmax, 5)


def estimate_rms_robust(data: np.ndarray, *, finite_only: bool = True) -> float:
    """Estimate a robust RMS using MAD.

    The estimate is intended for contour-level defaults, not for precision
    noise analysis. NaN/Inf values are ignored by default.
    """
    arr = np.asarray(data, dtype=float)
    if finite_only:
        arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    med = float(np.nanmedian(arr))
    mad = float(np.nanmedian(np.abs(arr - med)))
    if not np.isfinite(mad) or mad == 0.0:
        std = float(np.nanstd(arr))
        return std if np.isfinite(std) else float("nan")
    return 1.4826 * mad


def compute_contour_levels(
    data: np.ndarray,
    *,
    levels: Optional[Any] = None,
    level_mode: str = "auto",
    rms: Optional[float] = None,
    sigma_levels: Optional[Sequence[float]] = None,
    negative_sigma_levels: Optional[Sequence[float]] = None,
    fraction_levels: Optional[Sequence[float]] = None,
    level_scale: float = 1.0,
    sigma_scale: float = 1.0,
    symmetric: bool = False,
) -> np.ndarray:
    """Compute contour levels for auto, fraction, manual, or RMS-based modes.

    Parameters
    ----------
    levels
        Explicit contour factors or explicit contour values, depending on
        ``level_mode``. In ``'manual'``/``'explicit'`` mode, the returned levels
        are ``level_scale * levels``. In the other modes, a numeric ``levels``
        input is interpreted as already-final contour values and returned as-is
        for backward compatibility.
    level_mode
        One of ``'auto'``, ``'fraction'``, ``'manual'``, ``'explicit'``,
        ``'rms'``, or ``'sigma'``. ``'auto'`` intentionally follows the legacy
        behaviour and means legacy automatic fractional levels based on the
        positive peak. ``'sigma'`` is accepted as an alias of ``'rms'`` for
        backward compatibility.
    rms
        Noise estimate to use for RMS-based contour levels. If omitted in
        ``'rms'``/``'sigma'`` mode, a robust MAD-based estimate is derived from
        the 2D contour image as a convenience fallback. This is not intended to
        replace physically motivated error propagation from 3D cubes.
    sigma_levels
        Multipliers for RMS-based contour levels. Default is ``[3, 5, 7]``.
        Combined with ``sigma_scale`` so that, for example,
        ``sigma_levels=[1,2,3,4], sigma_scale=3`` yields ``3,6,9,12 × rms``.
    negative_sigma_levels
        Optional negative RMS multipliers. If omitted and ``symmetric=True``,
        mirrored negative levels are generated from the positive RMS levels.
    fraction_levels
        Fractions of the positive peak for fractional levels. Default is
        ``[0.1, 0.3, 0.5, 0.7, 0.9]``.
    level_scale
        Multiplicative scale used in ``'manual'``/``'explicit'`` mode.
        For example, ``levels=[1,2,3,4], level_scale=3`` yields ``[3,6,9,12]``.
    sigma_scale
        Extra multiplicative factor applied before ``sigma_levels`` in
        ``'rms'``/``'sigma'`` mode. For example,
        ``rms=0.8, sigma_levels=[1,2,3,4], sigma_scale=3`` yields
        ``[2.4, 4.8, 7.2, 9.6]``.
    symmetric
        If True, include negative mirrored levels in fraction or RMS-based
        modes.
    """
    level_mode = str(level_mode).lower()
    if level_mode == 'sigma':
        level_mode = 'rms'
    if level_mode == 'explicit':
        level_mode = 'manual'

    arr_levels = None
    if levels is not None and not (isinstance(levels, str) and levels.lower() == 'auto'):
        arr_levels = np.asarray(levels, dtype=float).ravel()
        arr_levels = arr_levels[np.isfinite(arr_levels)]
        if arr_levels.size == 0:
            raise ValueError('Explicit contour levels contain no finite values.')
        if level_mode == 'manual':
            arr = float(level_scale) * arr_levels
            return np.unique(np.sort(arr))
        # Backward compatibility: if explicit numeric levels are provided while
        # using a non-manual mode, interpret them as already-final contour
        # values and return them unchanged.
        return np.unique(np.sort(arr_levels))

    finite = np.asarray(data, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.array([0.0])

    if sigma_levels is None:
        sigma_levels = [3.0, 5.0, 7.0]
    if fraction_levels is None:
        fraction_levels = [0.1, 0.3, 0.5, 0.7, 0.9]

    def _fraction_based() -> np.ndarray:
        vmax = float(np.nanmax(finite))
        if vmax > 0:
            pos = vmax * np.asarray(fraction_levels, dtype=float)
            if symmetric:
                return np.unique(np.sort(np.concatenate([-pos[::-1], pos])))
            return np.unique(np.sort(pos))
        vmin = float(np.nanmin(finite))
        return np.linspace(vmin, vmax, 5)

    def _rms_based() -> np.ndarray:
        local_rms = float(rms) if rms is not None else estimate_rms_robust(finite)
        if not np.isfinite(local_rms) or local_rms <= 0:
            raise ValueError('Could not determine a positive finite RMS for RMS contour levels.')
        base = local_rms * float(sigma_scale)
        pos = base * np.asarray(sigma_levels, dtype=float)
        neg = np.array([], dtype=float)
        if negative_sigma_levels is not None:
            neg = -base * np.asarray(negative_sigma_levels, dtype=float)
        elif symmetric:
            neg = -pos[::-1]
        arr = np.concatenate([neg, pos]) if neg.size else pos
        arr = arr[np.isfinite(arr)]
        arr = arr[arr != 0]
        if arr.size == 0:
            raise ValueError('RMS contour levels collapsed to an empty set.')
        return np.unique(np.sort(arr))

    if level_mode == 'manual':
        raise ValueError("level_mode='manual' requires explicit numeric 'levels'.")
    if level_mode == 'fraction':
        return _fraction_based()
    if level_mode == 'rms':
        return _rms_based()
    if level_mode == 'auto':
        return _auto_contour_levels(np.asarray(data, dtype=float))
    raise ValueError(f'Unsupported contour level_mode: {level_mode}')



def _default_colorbar_label(map2d: Map2D) -> Optional[str]:
    if map2d.unit is None:
        return None
    if map2d.meta.get("mode") in {"moment0", "provisional_moment", "final_moment"}:
        return f"Integrated Intensity [{map2d.unit.to_string('latex')}]"
    return map2d.unit.to_string("latex")



def add_contours(
    ax: Any,
    contours: Sequence[Mapping[str, Any]],
    *,
    base_map: Optional[Map2D] = None,
) -> List[Any]:
    """Add one or more contour layers to an existing WCSAxes."""
    artists: List[Any] = []
    if not contours:
        return artists

    for layer in contours:
        layer = dict(layer)
        already_smoothed = False
        if base_map is not None and not any(k in layer for k in ("source", "data", "header")):
            cmap = base_map
        else:
            contour_map_keys = {
                "mode",
                "chan_range",
                "vel_range",
                "spectral_unit",
                "mask",
                "mask_mode",
                "zero_fill",
                "nan_fill",
                "linefree_ext",
                "linecand_ext",
                "basesup_ext",
                "final_mask_ext",
            }
            if any(k in layer for k in contour_map_keys):
                cmap = make_2d_map(
                    source=layer.get("source"),
                    ext=layer.get("ext"),
                    mode=layer.get("mode", "moment0"),
                    chan_range=layer.get("chan_range"),
                    vel_range=layer.get("vel_range"),
                    spectral_unit=layer.get("spectral_unit", "km/s"),
                    mask=layer.get("mask"),
                    mask_mode=layer.get("mask_mode"),
                    zero_fill=layer.get("zero_fill", False),
                    nan_fill=layer.get("nan_fill", True),
                    smooth_fwhm_arcsec=layer.get("smooth_fwhm_arcsec"),
                    target_hpbw_arcsec=layer.get("target_hpbw_arcsec"),
                    orig_hpbw_arcsec=layer.get("orig_hpbw_arcsec"),
                    linefree_ext=layer.get("linefree_ext", "LINEFREE"),
                    linecand_ext=layer.get("linecand_ext", "LINECAND3D"),
                    basesup_ext=layer.get("basesup_ext", "BASESUP3D"),
                    final_mask_ext=layer.get("final_mask_ext", "MASK3D"),
                )
                already_smoothed = (
                    layer.get("smooth_fwhm_arcsec") is not None
                    or layer.get("target_hpbw_arcsec") is not None
                )
            else:
                cmap = resolve_map_input(
                    source=layer.get("source"),
                    data=layer.get("data"),
                    header=layer.get("header"),
                    ext=layer.get("ext"),
                )

        if layer.get("crop"):
            cmap = crop_map2d(cmap, **dict(layer.get("crop")))

        cdata = np.asarray(cmap.data, dtype=float)
        if (not already_smoothed) and (
            layer.get("smooth_fwhm_arcsec") is not None or layer.get("target_hpbw_arcsec") is not None
        ):
            orig_hpbw, _ = _infer_orig_hpbw_arcsec(cmap.header)
            cdata = apply_gaussian_smoothing_2d(
                cdata,
                cmap.wcs,
                smooth_fwhm_arcsec=layer.get("smooth_fwhm_arcsec"),
                target_hpbw_arcsec=layer.get("target_hpbw_arcsec"),
                orig_hpbw_arcsec=layer.get("orig_hpbw_arcsec") or orig_hpbw,
            )

        levels = compute_contour_levels(
            cdata,
            levels=layer.get("levels", "auto"),
            level_mode=layer.get("level_mode", "auto"),
            rms=layer.get("rms"),
            sigma_levels=layer.get("sigma_levels"),
            negative_sigma_levels=layer.get("negative_sigma_levels"),
            fraction_levels=layer.get("fraction_levels"),
            level_scale=layer.get("level_scale", 1.0),
            sigma_scale=layer.get("sigma_scale", 1.0),
            symmetric=layer.get("symmetric", False),
        )

        artist = ax.contour(
            cdata,
            levels=levels,
            colors=layer.get("colors", "white"),
            linewidths=layer.get("linewidths", 1.0),
            linestyles=layer.get("linestyles", "solid"),
            alpha=layer.get("alpha", 0.7),
            transform=ax.get_transform(cmap.wcs),
        )
        label = layer.get("label")
        if label is not None:
            try:
                artist.collections[0].set_label(label)
            except Exception:
                pass
        artists.append(artist)
    return artists



def _add_legend(ax: Any, legend: Union[bool, Mapping[str, Any]]) -> Optional[Any]:
    if not legend:
        return None
    kwargs = dict(legend) if isinstance(legend, Mapping) else {}
    handles, labels = ax.get_legend_handles_labels()
    unique_h: List[Any] = []
    unique_l: List[str] = []
    seen = set()
    for handle, label in zip(handles, labels):
        if label in (None, "") or str(label).startswith("_"):
            continue
        label_str = str(label)
        if label_str in seen:
            continue
        seen.add(label_str)
        unique_h.append(handle)
        unique_l.append(label_str)
    if not unique_h:
        return None
    kwargs.setdefault("loc", "best")
    return ax.legend(unique_h, unique_l, **kwargs)



def plot_map(
    source: SourceLike = None,
    *,
    data: Optional[np.ndarray] = None,
    header: Optional[fits.Header] = None,
    ext: Optional[Union[int, str]] = None,
    ax: Optional[Any] = None,
    projection: Optional[WCS] = None,
    cmap: str = "viridis",
    norm: Optional[ImageNormalize] = None,
    norm_mode: str = "asinh",
    norm_percentile: Optional[Tuple[float, float]] = None,
    cmin: Optional[float] = None,
    cmax: Optional[float] = None,
    stretch_a: float = 0.1,
    power_gamma: float = 1.0,
    colorbar: bool = True,
    colorbar_label: Optional[str] = None,
    title: Optional[str] = None,
    origin: str = "lower",
    contours: Optional[Sequence[Mapping[str, Any]]] = None,
    grid: bool = True,
    beam: Optional[Union[str, Mapping[str, Any]]] = None,
    scalebar: Optional[Union[bool, Mapping[str, Any]]] = None,
    north_arrow: Optional[Union[bool, Mapping[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = False,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    crop: Optional[Mapping[str, Any]] = None,
    show_readout: bool = True,
    describe: bool = False,
    save: Optional[Union[str, Path]] = None,
    save_dpi: int = 250,
    save_transparent: bool = False,
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
    interpolation: str = "nearest",
    alpha: float = 1.0,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> Dict[str, Any]:
    """Plot a 2D map on WCSAxes with optional contours and convenience helpers.

    ``cmin`` / ``cmax`` set the display color-scale limits.
    Deprecated aliases ``vmin`` / ``vmax`` are still accepted for backward compatibility.
    """
    import matplotlib.pyplot as plt

    if vmin is not None or vmax is not None:
        warnings.warn(
            "plot_map(): vmin/vmax are deprecated and mean display color-scale limits; use cmin/cmax instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if cmin is None and vmin is not None:
            cmin = vmin
        if cmax is None and vmax is not None:
            cmax = vmax

    map2d = resolve_map_input(source=source, data=data, header=header, ext=ext)
    if crop:
        map2d = crop_map2d(map2d, **dict(crop))
    proj = projection or map2d.wcs

    created_fig = False
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=proj)
        created_fig = True
    else:
        fig = ax.figure

    normalize_info = {
        "mode": norm_mode,
        "percentile": norm_percentile,
        "cmin": cmin,
        "cmax": cmax,
        "stretch_a": stretch_a,
        "power_gamma": power_gamma,
    }
    if norm is None:
        norm = build_normalize(
            map2d.data,
            mode=norm_mode,
            percentile=norm_percentile,
            cmin=cmin,
            cmax=cmax,
            vmin=vmin,
            vmax=vmax,
            stretch_a=stretch_a,
            power_gamma=power_gamma,
        )
    else:
        normalize_info["mode"] = "external_norm"

    im = ax.imshow(
        map2d.data,
        origin=origin,
        cmap=cmap,
        norm=norm,
        interpolation=interpolation,
        alpha=alpha,
        transform=ax.get_transform(map2d.wcs),
    )
    _disable_artist_cursor_data(im)

    cbar = None
    if colorbar:
        cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        _disable_axes_cursor_readout(cbar.ax)
        if colorbar_label is None:
            colorbar_label = _default_colorbar_label(map2d)
        if colorbar_label:
            cbar.set_label(colorbar_label, rotation=270, va="bottom")

    if grid:
        ax.coords.grid(color="white", alpha=0.25, linestyle="solid")

    default_xlabel, default_ylabel = _default_axis_labels(map2d.wcs)
    ax.set_xlabel(xlabel or default_xlabel)
    ax.set_ylabel(ylabel or default_ylabel)

    if title is not None:
        ax.set_title(title)

    contour_artists = add_contours(ax, contours or [], base_map=map2d)

    beam_artist = None
    if beam is not None:
        beam_artist = _add_beam_artist(ax, map2d, beam)

    scalebar_artists = None
    if scalebar:
        scalebar_artists = add_scalebar(ax, map2d, scalebar)

    north_arrow_artists = None
    if north_arrow:
        north_arrow_artists = add_north_arrow(ax, map2d, north_arrow)

    if show_readout:
        _apply_readout_formatter(ax, map2d)

    if describe:
        print(describe_map2d(map2d))

    legend_artist = _add_legend(ax, legend)

    _save_figure(fig, save, dpi=save_dpi, transparent=save_transparent)

    if show and created_fig:
        plt.show()

    out_meta = dict(map2d.meta)
    out_meta["normalize_info"] = normalize_info

    return {
        "fig": fig,
        "ax": ax,
        "image": im,
        "colorbar": cbar,
        "contours": contour_artists,
        "beam": beam_artist,
        "scalebar": scalebar_artists,
        "north_arrow": north_arrow_artists,
        "legend": legend_artist,
        "description": describe_map2d(map2d),
        "map2d": Map2D(data=map2d.data, header=map2d.header, wcs=map2d.wcs, unit=map2d.unit, meta=out_meta),
    }


# -----------------------------------------------------------------------------
# RGB support
# -----------------------------------------------------------------------------


def _normalize_rgb_channel(
    data: np.ndarray,
    *,
    norm_mode: str = "asinh",
    percentile: Optional[Tuple[float, float]] = (1.0, 99.5),
    cmin: Optional[float] = None,
    cmax: Optional[float] = None,
    stretch_a: float = 0.1,
    power_gamma: float = 1.0,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> np.ndarray:
    norm = build_normalize(
        data,
        mode=norm_mode,
        percentile=percentile,
        cmin=cmin,
        cmax=cmax,
        vmin=vmin,
        vmax=vmax,
        stretch_a=stretch_a,
        power_gamma=power_gamma,
    )
    arr = np.asarray(norm(data), dtype=float)
    arr = np.where(np.isfinite(arr), arr, 0.0)
    return np.clip(arr, 0.0, 1.0)



def _resolve_rgb_component(component: Any) -> Map2D:
    if isinstance(component, Map2D):
        return component
    if isinstance(component, dict):
        if component.get("mode") in {None, "identity", "map", "2d"} and (
            component.get("data") is not None or component.get("source") is not None
        ):
            return resolve_map_input(
                source=component.get("source"),
                data=component.get("data"),
                header=component.get("header"),
                ext=component.get("ext"),
            )
        return make_2d_map(**component)
    return resolve_map_input(component)



def make_rgb_map(
    red: Any,
    green: Any,
    blue: Any,
    *,
    red_norm: Optional[Mapping[str, Any]] = None,
    green_norm: Optional[Mapping[str, Any]] = None,
    blue_norm: Optional[Mapping[str, Any]] = None,
) -> RGBMap:
    """Create an RGB image from three 2D maps or map specifications."""
    rmap = _resolve_rgb_component(red)
    gmap = _resolve_rgb_component(green)
    bmap = _resolve_rgb_component(blue)

    if rmap.data.shape != gmap.data.shape or rmap.data.shape != bmap.data.shape:
        raise ValueError("RGB inputs must have the same 2D shape in this first implementation.")
    if rmap.wcs.to_header_string() != gmap.wcs.to_header_string() or rmap.wcs.to_header_string() != bmap.wcs.to_header_string():
        warnings.warn(
            "RGB inputs do not have identical WCS headers. In this first implementation, arrays are stacked without reprojection.",
            RuntimeWarning,
        )

    r = _normalize_rgb_channel(rmap.data, **dict(red_norm or {}))
    g = _normalize_rgb_channel(gmap.data, **dict(green_norm or {}))
    b = _normalize_rgb_channel(bmap.data, **dict(blue_norm or {}))
    rgb = np.dstack([r, g, b])

    meta = {
        "red_meta": rmap.meta,
        "green_meta": gmap.meta,
        "blue_meta": bmap.meta,
    }
    return RGBMap(rgb=rgb, header=rmap.header, wcs=rmap.wcs, meta=meta)



def plot_rgb(
    rgb_source: Optional[RGBMap] = None,
    *,
    red: Any = None,
    green: Any = None,
    blue: Any = None,
    ax: Optional[Any] = None,
    title: Optional[str] = None,
    grid: bool = True,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    crop: Optional[Mapping[str, Any]] = None,
    show_readout: bool = False,
    scalebar: Optional[Union[bool, Mapping[str, Any]]] = None,
    north_arrow: Optional[Union[bool, Mapping[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = False,
    describe: bool = False,
    save: Optional[Union[str, Path]] = None,
    save_dpi: int = 250,
    save_transparent: bool = False,
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
    **rgb_kwargs: Any,
) -> Dict[str, Any]:
    """Plot an RGB map on WCSAxes."""
    import matplotlib.pyplot as plt

    if rgb_source is None:
        rgb_source = make_rgb_map(red, green, blue, **rgb_kwargs)

    if crop:
        # Apply an identical celestial cutout to each channel.
        chans = []
        cc: Optional[Map2D] = None
        for i in range(3):
            cm = Map2D(
                data=np.asarray(rgb_source.rgb[..., i], dtype=float),
                header=rgb_source.header,
                wcs=rgb_source.wcs,
                unit=None,
                meta=dict(rgb_source.meta),
            )
            cc = crop_map2d(cm, **dict(crop))
            chans.append(np.asarray(cc.data, dtype=float))
        if cc is None:
            raise RuntimeError("RGB crop requested, but no channels were processed.")
        rgb_source = RGBMap(rgb=np.dstack(chans), header=cc.header, wcs=cc.wcs, meta=dict(rgb_source.meta))

    created_fig = False
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=rgb_source.wcs)
        created_fig = True
    else:
        fig = ax.figure

    im = ax.imshow(rgb_source.rgb, origin="lower", transform=ax.get_transform(rgb_source.wcs))
    _disable_artist_cursor_data(im)
    if grid:
        ax.coords.grid(color="white", alpha=0.25, linestyle="solid")

    default_xlabel, default_ylabel = _default_axis_labels(rgb_source.wcs)
    ax.set_xlabel(xlabel or default_xlabel)
    ax.set_ylabel(ylabel or default_ylabel)
    if title is not None:
        ax.set_title(title)

    tmp = Map2D(data=np.asarray(np.nanmean(rgb_source.rgb, axis=2), dtype=float), header=rgb_source.header, wcs=rgb_source.wcs, unit=None, meta=dict(rgb_source.meta))

    if scalebar:
        add_scalebar(ax, tmp, scalebar)
    if north_arrow:
        add_north_arrow(ax, tmp, north_arrow)

    if show_readout:
        _apply_readout_formatter(ax, tmp)

    desc = f"RGB map\nshape: {rgb_source.rgb.shape[1]} x {rgb_source.rgb.shape[0]} pixels\ncoord_system: {_infer_native_coord_system(rgb_source.wcs)}"
    if describe:
        print(desc)

    legend_artist = _add_legend(ax, legend)

    _save_figure(fig, save, dpi=save_dpi, transparent=save_transparent)

    if show and created_fig:
        plt.show()

    return {"fig": fig, "ax": ax, "image": im, "rgb": rgb_source, "legend": legend_artist, "description": desc}


def plot_scene(
    base: Any,
    *,
    overlays: Optional[Sequence[Any]] = None,
    ax: Optional[Any] = None,
    title: Optional[str] = None,
    grid: bool = True,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    colorbar: bool = True,
    show_readout: bool = True,
    describe: bool = False,
    save: Optional[Union[str, Path]] = None,
    save_dpi: int = 250,
    save_transparent: bool = False,
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
    legend: Union[bool, Mapping[str, Any]] = False,
    crop: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """High-level plotting API with one base layer and optional overlays.

    Parameters
    ----------
    base : source or mapping
        Base layer specification. Supported base kinds are ``image`` and ``rgb``.
        Bare sources are interpreted as ``{'kind': 'image', 'source': source}``.
    overlays : sequence
        Optional layer specifications. Supported kinds in this implementation:
        ``image``, ``contour``, ``marker``, ``catalog``, ``text``, ``beam``,
        ``scalebar``, and ``north_arrow``.
    """
    overlays = list(overlays or [])
    scene_created_fig = ax is None
    base_layer = _resolve_base_spec(base)
    kind = str(base_layer.get("kind", "image")).lower()

    if kind == "rgb":
        result = plot_rgb(
            rgb_source=base_layer.get("rgb_source"),
            red=base_layer.get("red"),
            green=base_layer.get("green"),
            blue=base_layer.get("blue"),
            ax=ax,
            title=title or base_layer.get("title"),
            grid=grid if base_layer.get("grid") is None else base_layer.get("grid"),
            xlabel=xlabel or base_layer.get("xlabel"),
            ylabel=ylabel or base_layer.get("ylabel"),
            crop=crop or base_layer.get("crop"),
            show_readout=show_readout,
            describe=describe or base_layer.get("describe", False),
            save=None,
            show=False,
            figsize=figsize,
            scalebar=base_layer.get("scalebar"),
            north_arrow=base_layer.get("north_arrow"),
            legend=False,
            red_norm=base_layer.get("red_norm"),
            green_norm=base_layer.get("green_norm"),
            blue_norm=base_layer.get("blue_norm"),
        )
        fig = result["fig"]
        ax = result["ax"]
        base_wcs = result["rgb"].wcs
        base_description = result["description"]
        base_obj = result["rgb"]
    else:
        base_map = _layer_to_map2d(base_layer, default_mode=base_layer.get("mode", "identity") or "identity")
        result = plot_map(
            base_map,
            ax=ax,
            projection=base_layer.get("projection"),
            cmap=base_layer.get("cmap", "viridis"),
            norm=base_layer.get("norm"),
            norm_mode=base_layer.get("norm_mode", "asinh"),
            norm_percentile=base_layer.get("norm_percentile"),
            cmin=base_layer.get("cmin"),
            cmax=base_layer.get("cmax"),
            stretch_a=base_layer.get("stretch_a", 0.1),
            power_gamma=base_layer.get("power_gamma", 1.0),
            colorbar=colorbar if base_layer.get("colorbar") is None else base_layer.get("colorbar"),
            colorbar_label=base_layer.get("colorbar_label"),
            title=title or base_layer.get("title"),
            origin=base_layer.get("origin", "lower"),
            contours=None,
            grid=grid if base_layer.get("grid") is None else base_layer.get("grid"),
            beam=base_layer.get("beam"),
            xlabel=xlabel or base_layer.get("xlabel"),
            ylabel=ylabel or base_layer.get("ylabel"),
            crop=crop or base_layer.get("crop"),
            show_readout=show_readout,
            describe=describe or base_layer.get("describe", False),
            save=None,
            show=False,
            figsize=figsize,
            interpolation=base_layer.get("interpolation", "nearest"),
            alpha=base_layer.get("alpha", 1.0),
            scalebar=base_layer.get("scalebar"),
            north_arrow=base_layer.get("north_arrow"),
            legend=False,
        )
        fig = result["fig"]
        ax = result["ax"]
        base_wcs = result["map2d"].wcs
        base_description = result["description"]
        base_obj = result["map2d"]

    contour_specs: List[Dict[str, Any]] = []
    overlay_artists: List[Any] = []
    overlay_maps: List[Map2D] = []
    beam_artists: List[Any] = []
    base_shape: Optional[Tuple[int, int]] = None
    if isinstance(base_obj, Map2D):
        base_shape = tuple(base_obj.data.shape)
    elif isinstance(base_obj, RGBMap):
        base_shape = tuple(base_obj.rgb.shape[:2])

    for overlay in overlays:
        layer = _layer_dict(overlay)
        layer_kind = str(layer.get("kind", "contour")).lower()
        if layer_kind == "contour":
            contour_specs.append(_prepare_contour_layer(layer))
        elif layer_kind == "image":
            info = _add_image_overlay(ax, layer, target_wcs=base_wcs, target_shape=base_shape)
            overlay_artists.append(info["artist"])
            overlay_maps.append(info["map2d"])
        elif layer_kind == "marker":
            overlay_artists.append(_scatter_points(ax, base_wcs, layer))
        elif layer_kind == "catalog":
            info = _scatter_catalog(ax, base_wcs, layer)
            overlay_artists.append(info["marker"])
            overlay_artists.extend(info["texts"])
        elif layer_kind == "text":
            overlay_artists.append(_add_text_artist(ax, base_wcs, layer))
        elif layer_kind == "beam":
            beam_cfg = layer.get("beam", layer.get("config", "auto"))
            ann_map = _annotation_map_from_base(base_obj)
            if ann_map is not None:
                art = _add_beam_artist(ax, ann_map, beam_cfg)
                if art is not None:
                    _set_artist_label(art, layer.get("label"))
                    beam_artists.append(art)
            else:
                warnings.warn("beam overlay was requested, but the current base has no beam-compatible WCS image.", RuntimeWarning)
        elif layer_kind == "scalebar":
            ann_map = _annotation_map_from_base(base_obj)
            if ann_map is not None:
                overlay_artists.append(add_scalebar(ax, ann_map, layer.get("config", layer)))
        elif layer_kind in {"north_arrow", "compass"}:
            ann_map = _annotation_map_from_base(base_obj)
            if ann_map is not None:
                overlay_artists.append(add_north_arrow(ax, ann_map, layer.get("config", layer)))
        else:
            raise ValueError(f"Unsupported overlay kind: {layer_kind}")

    contour_artists: List[Any] = []
    if contour_specs:
        contour_artists = add_contours(ax, contour_specs, base_map=base_obj if isinstance(base_obj, Map2D) else None)

    legend_artist = _add_legend(ax, legend)

    _save_figure(fig, save, dpi=save_dpi, transparent=save_transparent)

    if show and scene_created_fig:
        import matplotlib.pyplot as plt
        plt.show()

    return {
        "fig": fig,
        "ax": ax,
        "base": base_obj,
        "base_description": base_description,
        "overlay_artists": overlay_artists,
        "overlay_maps": overlay_maps,
        "contours": contour_artists,
        "beam_artists": beam_artists,
        "legend": legend_artist,
    }


def quicklook(
    source: Any,
    *,
    ext: Optional[Union[int, str]] = None,
    mode: str = "moment0",
    chan_range: Optional[Tuple[int, int]] = None,
    vel_range: Optional[Tuple[float, float]] = None,
    spectral_unit: Union[str, u.UnitBase] = "km/s",
    smooth_fwhm_arcsec: Optional[float] = None,
    target_hpbw_arcsec: Optional[float] = None,
    orig_hpbw_arcsec: Optional[float] = None,
    cmap: str = "viridis",
    norm_mode: str = "asinh",
    norm_percentile: Optional[Tuple[float, float]] = (1.0, 99.5),
    cmin: Optional[float] = None,
    cmax: Optional[float] = None,
    contours: Optional[Sequence[Mapping[str, Any]]] = None,
    title: Optional[str] = None,
    grid: bool = True,
    beam: Optional[Union[str, Mapping[str, Any]]] = None,
    crop: Optional[Mapping[str, Any]] = None,
    describe: bool = False,
    help: bool = False,
    save: Optional[Union[str, Path]] = None,
    save_dpi: int = 250,
    save_transparent: bool = False,
    show_readout: bool = True,
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
) -> Dict[str, Any]:
    """Quick entry point for plotting one source with minimal required options.

    Parameters are intentionally biased toward the most common single-panel use case.
    ``contours`` is accepted as a convenience shortcut and is forwarded as contour
    overlays on top of the same base source unless the contour spec explicitly
    overrides ``source`` / reduction parameters.
    """
    if help:
        print_plotting_help()

    base = {
        "kind": "image",
        "source": source,
        "ext": ext,
        "mode": mode,
        "chan_range": chan_range,
        "vel_range": vel_range,
        "spectral_unit": spectral_unit,
        "smooth_fwhm_arcsec": smooth_fwhm_arcsec,
        "target_hpbw_arcsec": target_hpbw_arcsec,
        "orig_hpbw_arcsec": orig_hpbw_arcsec,
        "cmap": cmap,
        "norm_mode": norm_mode,
        "norm_percentile": norm_percentile,
        "cmin": cmin,
        "cmax": cmax,
        "beam": beam,
        "crop": crop,
    }

    overlay_list: Optional[List[Dict[str, Any]]] = None
    if contours:
        overlay_list = []
        base_for_contours = {
            "source": source,
            "ext": ext,
            "mode": mode,
            "chan_range": chan_range,
            "vel_range": vel_range,
            "spectral_unit": spectral_unit,
            "smooth_fwhm_arcsec": smooth_fwhm_arcsec,
            "target_hpbw_arcsec": target_hpbw_arcsec,
            "orig_hpbw_arcsec": orig_hpbw_arcsec,
            "crop": crop,
        }
        for spec in contours:
            ov = dict(base_for_contours)
            ov.update(dict(spec))
            ov["kind"] = "contour"
            overlay_list.append(ov)

    return plot_scene(
        base,
        overlays=overlay_list,
        title=title,
        grid=grid,
        describe=describe,
        save=save,
        save_dpi=save_dpi,
        save_transparent=save_transparent,
        show_readout=show_readout,
        show=show,
        figsize=figsize,
    )


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def _parse_comma_pair(text: Optional[str], cast=float) -> Optional[Tuple[Any, Any]]:
    if text is None:
        return None
    parts = [p.strip() for p in str(text).split(",")]
    if len(parts) != 2:
        raise ValueError(f"Expected two comma-separated values, got: {text}")
    return cast(parts[0]), cast(parts[1])



def main(argv: Optional[Sequence[str]] = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Plot 2D/3D radio astronomy FITS data on WCSAxes.")
    parser.add_argument("--api-help", action="store_true", help="Print compact Python API help and exit.")
    parser.add_argument("source", nargs="?", help="Input FITS file.")
    parser.add_argument("--ext", default=None, help="Primary HDU or extension name/index.")
    parser.add_argument(
        "--mode",
        default="moment0",
        choices=[
            "identity",
            "moment0",
            "channel_sum",
            "channel_mean",
            "channel_slice",
            "velocity_slice",
            "provisional_moment",
            "final_moment",
        ],
        help="How to make the 2D map.",
    )
    parser.add_argument("--chan-range", default=None, help="Channel range as lo,hi.")
    parser.add_argument("--vel-range", default=None, help="Velocity range as lo,hi in --spectral-unit.")
    parser.add_argument("--spectral-unit", default="km/s", help="Spectral unit for velocity-axis operations.")
    parser.add_argument("--smooth-fwhm-arcsec", type=float, default=None, help="Additional smoothing FWHM in arcsec.")
    parser.add_argument("--target-hpbw-arcsec", type=float, default=None, help="Target HPBW in arcsec.")
    parser.add_argument("--orig-hpbw-arcsec", type=float, default=None, help="Original HPBW in arcsec.")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap name.")
    parser.add_argument("--norm-mode", default="asinh", choices=["linear", "sqrt", "log", "asinh", "power"])
    parser.add_argument("--norm-percentile", default="1,99.5", help="Percentile clip as lo,hi.")
    parser.add_argument("--cmin", type=float, default=None, help="Display color-scale minimum.")
    parser.add_argument("--cmax", type=float, default=None, help="Display color-scale maximum.")
    parser.add_argument("--vmin", type=float, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--vmax", type=float, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--stretch-a", type=float, default=0.1)
    parser.add_argument("--power-gamma", type=float, default=1.0)
    parser.add_argument("--title", default=None)
    parser.add_argument("--describe", action="store_true", help="Print a short summary of the plotted map.")
    parser.add_argument("--beam", default=None, help="Beam setting: auto/header or omit.")
    parser.add_argument("--legend", action="store_true", help="Show a legend for labeled overlays/contours.")
    parser.add_argument("--scalebar", action="store_true", help="Draw a simple angular scale bar.")
    parser.add_argument("--scalebar-length", type=float, default=None, help="Scale-bar length in arcsec.")
    parser.add_argument("--north-arrow", action="store_true", help="Draw north/east orientation arrows.")
    parser.add_argument("--no-grid", dest="grid", action="store_false", default=True, help="Disable coordinate grid.")
    parser.add_argument("--output", default=None, help="Save figure to this path.")
    parser.add_argument("--dpi", type=int, default=250)
    parser.add_argument("--transparent", action="store_true", help="Save figure with transparent background.")
    parser.add_argument("--contour-levels", default=None, help="Contour levels as v1,v2,v3 or 'auto'. In manual mode these are treated as factors and multiplied by --contour-level-scale.")
    parser.add_argument("--contour-level-mode", default="auto", choices=["auto", "fraction", "manual", "explicit", "rms", "sigma"], help="How contour levels are generated: auto=fraction (legacy), fraction, manual/explicit, or rms/sigma.")
    parser.add_argument("--contour-level-scale", type=float, default=1.0, help="Scale factor for manual contour levels. Example: levels=1,2,3,4 with scale=3 gives 3,6,9,12.")
    parser.add_argument("--contour-fractions", default=None, help="Fraction levels as f1,f2,f3 for fraction/auto contour mode.")
    parser.add_argument("--contour-sigmas", default=None, help="RMS multipliers as s1,s2,s3 for rms/sigma contour mode.")
    parser.add_argument("--contour-sigma-scale", type=float, default=1.0, help="Extra scale factor for rms/sigma contour mode. Example: sigmas=1,2,3,4 with sigma-scale=3 gives 3,6,9,12 times rms.")
    parser.add_argument("--contour-negative", action="store_true", help="Add mirrored negative contour levels in fraction or rms contour modes.")
    parser.add_argument("--contour-rms", type=float, default=None, help="Explicit RMS to use for rms/sigma contour mode.")
    parser.add_argument("--contour-color", default="white")
    parser.add_argument("--contour-linewidth", type=float, default=0.8)
    parser.add_argument("--zero-fill", action="store_true", help="Use zero fill outside the mask.")
    parser.add_argument(
        "--no-nan-fill",
        dest="nan_fill",
        action="store_false",
        default=True,
        help="Disable NaN fill outside the mask.",
    )
    argv_list = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv_list)

    if args.api_help:
        print_plotting_help()
        if args.source is None:
            return

    if args.source is None:
        parser.error("source is required unless --api-help is used by itself.")

    ext = args.ext
    if isinstance(ext, str) and ext.isdigit():
        ext = int(ext)

    chan_range = _parse_comma_pair(args.chan_range, cast=int)
    vel_range = _parse_comma_pair(args.vel_range, cast=float)
    percentile = _parse_comma_pair(args.norm_percentile, cast=float) if args.norm_percentile else None

    map2d = make_2d_map(
        args.source,
        ext=ext,
        mode=args.mode,
        chan_range=chan_range,
        vel_range=vel_range,
        spectral_unit=args.spectral_unit,
        zero_fill=args.zero_fill,
        nan_fill=args.nan_fill,
        smooth_fwhm_arcsec=args.smooth_fwhm_arcsec,
        target_hpbw_arcsec=args.target_hpbw_arcsec,
        orig_hpbw_arcsec=args.orig_hpbw_arcsec,
    )

    contour_flag_names = (
        "--contour-levels",
        "--contour-level-mode",
        "--contour-level-scale",
        "--contour-fractions",
        "--contour-sigmas",
        "--contour-sigma-scale",
        "--contour-negative",
        "--contour-rms",
        "--contour-color",
        "--contour-linewidth",
    )
    wants_contours = any(
        str(token).startswith(contour_flag_names)
        for token in argv_list
    )

    contours = None
    if wants_contours:
        levels = None
        if args.contour_levels is not None:
            if args.contour_levels.strip().lower() == "auto":
                levels = "auto"
            else:
                levels = [float(v.strip()) for v in args.contour_levels.split(",") if v.strip()]
        sigma_levels = None
        if args.contour_sigmas:
            sigma_levels = [float(v.strip()) for v in args.contour_sigmas.split(",") if v.strip()]
        fraction_levels = None
        if args.contour_fractions:
            fraction_levels = [float(v.strip()) for v in args.contour_fractions.split(",") if v.strip()]
        contours = [
            {
                "levels": levels,
                "level_mode": args.contour_level_mode,
                "level_scale": args.contour_level_scale,
                "sigma_levels": sigma_levels,
                "sigma_scale": args.contour_sigma_scale,
                "fraction_levels": fraction_levels,
                "symmetric": args.contour_negative,
                "rms": args.contour_rms,
                "colors": args.contour_color,
                "linewidths": args.contour_linewidth,
            }
        ]

    if args.vmin is not None or args.vmax is not None:
        warnings.warn(
            "CLI options --vmin/--vmax are deprecated and mean display color-scale limits; use --cmin/--cmax instead.",
            DeprecationWarning,
            stacklevel=2,
        )
    cmin = args.cmin if args.cmin is not None else args.vmin
    cmax = args.cmax if args.cmax is not None else args.vmax

    result = plot_map(
        map2d,
        cmap=args.cmap,
        norm_mode=args.norm_mode,
        norm_percentile=percentile,
        cmin=cmin,
        cmax=cmax,
        stretch_a=args.stretch_a,
        power_gamma=args.power_gamma,
        contours=contours,
        title=args.title,
        grid=args.grid,
        beam=args.beam,
        describe=args.describe,
        save=args.output,
        save_dpi=args.dpi,
        save_transparent=args.transparent,
        show=args.output is None,
        scalebar={"length": args.scalebar_length, "units": "arcsec"} if args.scalebar else None,
        north_arrow=True if args.north_arrow else None,
        legend=args.legend,
    )



if __name__ == "__main__":
    main()