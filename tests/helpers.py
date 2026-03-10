from __future__ import annotations

import numpy as np
from astropy.io import fits


ANALYSIS_NAMES = {
    "MASK3D", "BASESUP3D", "LINECAND3D", "RMS", "BASE_RMS",
    "MOMENT0", "MOM0_BASESUP", "MOM0_LINECAND",
}


def make_cube(shape=(4, 3, 2), value=0.0, dtype=np.float32):
    arr = np.full(shape, value, dtype=dtype)
    return arr



def make_spectral_header(
    *,
    fits_spectral_axis: int = 3,
    ctype: str = "VRAD",
    cunit: str = "m/s",
    crval: float = 0.0,
    cdelt: float = 1000.0,
    crpix: float = 1.0,
    restfreq: float | None = None,
) -> fits.Header:
    """
    Build a minimal FITS header for a 3D cube.

    FITS axis numbering is 1-based; numpy shape is reversed relative to FITS axes.
    """
    hdr = fits.Header()
    spatial = ["RA---CAR", "DEC--CAR", "STOKES"]
    for ax in (1, 2, 3):
        hdr[f"CTYPE{ax}"] = spatial[ax - 1]
        hdr[f"CRVAL{ax}"] = 0.0
        hdr[f"CDELT{ax}"] = 1.0
        hdr[f"CRPIX{ax}"] = 1.0
        hdr[f"CUNIT{ax}"] = "deg"
    hdr[f"CTYPE{fits_spectral_axis}"] = ctype
    hdr[f"CRVAL{fits_spectral_axis}"] = crval
    hdr[f"CDELT{fits_spectral_axis}"] = cdelt
    hdr[f"CRPIX{fits_spectral_axis}"] = crpix
    hdr[f"CUNIT{fits_spectral_axis}"] = cunit
    if restfreq is not None:
        hdr["RESTFRQ"] = float(restfreq)
    return hdr



def make_image_hdu(
    data,
    *,
    name: str | None = None,
    header: fits.Header | None = None,
    btype: str | None = None,
    primary: bool = False,
):
    hdr = fits.Header() if header is None else header.copy()
    if btype is not None:
        hdr["BTYPE"] = btype
    if primary:
        return fits.PrimaryHDU(data=data, header=hdr)
    return fits.ImageHDU(data=data, header=hdr, name=name)
