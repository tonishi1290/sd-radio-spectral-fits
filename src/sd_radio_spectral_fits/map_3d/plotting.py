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
- This first implementation intentionally does not perform reprojection.
- Contours are drawn either on the same WCS/shape or via the WCSAxes
  transform when the target axes can interpret the input WCS directly.
- For provisional moments, this implementation assumes:
    * LINECAND3D : True means line/signal candidate -> integrate there
    * BASESUP3D  : True means baseline support      -> integrate complement
    * LINEFREE   : True means line-free/baseline    -> integrate complement
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union
import warnings

import numpy as np

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
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
    "apply_gaussian_smoothing_2d",
    "make_2d_map",
    "make_provisional_moment",
    "make_final_moment",
    "add_contours",
    "plot_map",
    "make_rgb_map",
    "plot_rgb",
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
    """Resolve 2D input from file, HDU, (header, data), Projection, or Map2D."""
    if isinstance(source, Map2D):
        return source

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
        "2D HDU, (header, data), Projection-like object, or Map2D."
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
    """Resolve a 3D cube from file, HDU, (header, data), or SpectralCube."""
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
        "3D HDU, (header, data), or SpectralCube."
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
    hdul, close_after = _open_hdul_if_needed(source)
    if hdul is None:
        raise ValueError(f"Associated mask HDU '{ext}' requires source to be a FITS path or HDUList.")
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
    hdul, close_after = _open_hdul_if_needed(source)
    if hdul is None:
        raise ValueError(
            "Automatic provisional moment requires source to be a FITS path or HDUList so that "
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
            else:
                cmap = resolve_map_input(
                    source=layer.get("source"),
                    data=layer.get("data"),
                    header=layer.get("header"),
                    ext=layer.get("ext"),
                )

        cdata = np.asarray(cmap.data, dtype=float)
        if layer.get("smooth_fwhm_arcsec") is not None or layer.get("target_hpbw_arcsec") is not None:
            orig_hpbw, _ = _infer_orig_hpbw_arcsec(cmap.header)
            cdata = apply_gaussian_smoothing_2d(
                cdata,
                cmap.wcs,
                smooth_fwhm_arcsec=layer.get("smooth_fwhm_arcsec"),
                target_hpbw_arcsec=layer.get("target_hpbw_arcsec"),
                orig_hpbw_arcsec=layer.get("orig_hpbw_arcsec") or orig_hpbw,
            )

        levels = layer.get("levels", "auto")
        if isinstance(levels, str) and levels.lower() == "auto":
            levels = _auto_contour_levels(cdata)

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
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> Dict[str, Any]:
    """Plot a 2D map on WCSAxes with optional contours.

    ``cmin`` / ``cmax`` set the display color-scale limits.
    Deprecated aliases ``vmin`` / ``vmax`` are still accepted for backward compatibility.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

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

    im = ax.imshow(map2d.data, origin=origin, cmap=cmap, norm=norm, transform=ax.get_transform(map2d.wcs))

    cbar = None
    if colorbar:
        cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        if colorbar_label is None:
            colorbar_label = _default_colorbar_label(map2d)
        if colorbar_label:
            cbar.set_label(colorbar_label)

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
        beam_cfg: Dict[str, Any]
        if isinstance(beam, str) and beam.lower() in {"header", "auto"}:
            smooth_info = map2d.meta.get("smooth_info", {}) if isinstance(map2d.meta, dict) else {}
            eff_hpbw = smooth_info.get("effective_hpbw_arcsec")
            if eff_hpbw is not None:
                beam_cfg = {"major": float(eff_hpbw), "minor": float(eff_hpbw), "angle": 0.0, "units": "arcsec"}
            else:
                bmaj = map2d.header.get("BMAJ")
                bmin = map2d.header.get("BMIN", bmaj)
                bpa = map2d.header.get("BPA", 0.0)
                if bmaj is not None:
                    beam_cfg = {"major": float(bmaj), "minor": float(bmin), "angle": float(bpa), "units": "deg"}
                else:
                    beam_cfg = {}
        else:
            beam_cfg = dict(beam)

        if beam_cfg:
            pix_arcsec = _pixel_scale_arcsec(map2d.wcs)
            major = beam_cfg.get("major")
            minor = beam_cfg.get("minor", major)
            angle = beam_cfg.get("angle", 0.0)
            units = str(beam_cfg.get("units", "arcsec")).lower()
            if major is not None:
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
                beam_artist = Ellipse(
                    (x0, y0),
                    width=major_pix,
                    height=minor_pix,
                    angle=angle,
                    facecolor=beam_cfg.get("facecolor", "none"),
                    edgecolor=beam_cfg.get("edgecolor", "white"),
                    linewidth=beam_cfg.get("linewidth", 1.2),
                    alpha=beam_cfg.get("alpha", 0.9),
                )
                ax.add_patch(beam_artist)

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
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
    **rgb_kwargs: Any,
) -> Dict[str, Any]:
    """Plot an RGB map on WCSAxes."""
    import matplotlib.pyplot as plt

    if rgb_source is None:
        rgb_source = make_rgb_map(red, green, blue, **rgb_kwargs)

    created_fig = False
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=rgb_source.wcs)
        created_fig = True
    else:
        fig = ax.figure

    im = ax.imshow(rgb_source.rgb, origin="lower", transform=ax.get_transform(rgb_source.wcs))
    if grid:
        ax.coords.grid(color="white", alpha=0.25, linestyle="solid")

    default_xlabel, default_ylabel = _default_axis_labels(rgb_source.wcs)
    ax.set_xlabel(xlabel or default_xlabel)
    ax.set_ylabel(ylabel or default_ylabel)
    if title is not None:
        ax.set_title(title)

    if show and created_fig:
        plt.show()

    return {"fig": fig, "ax": ax, "image": im, "rgb": rgb_source}


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
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description="Plot 2D/3D radio astronomy FITS data on WCSAxes.")
    parser.add_argument("source", help="Input FITS file.")
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
    parser.add_argument("--output", default=None, help="Save figure to this path.")
    parser.add_argument("--dpi", type=int, default=250)
    parser.add_argument("--contour-levels", default=None, help="Contour levels as v1,v2,v3 or 'auto'.")
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
    args = parser.parse_args(argv)

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

    contours = None
    if args.contour_levels is not None:
        if args.contour_levels.strip().lower() == "auto":
            levels: Union[str, List[float]] = "auto"
        else:
            levels = [float(v.strip()) for v in args.contour_levels.split(",") if v.strip()]
        contours = [
            {
                "levels": levels,
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
        show=args.output is None,
    )

    if args.output is not None:
        result["fig"].savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()
