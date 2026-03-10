from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from tests.helpers import make_cube, make_image_hdu, make_spectral_header
from sd_radio_spectral_fits.map_3d.cube_baseline.cli import (
    _get_cube_from_hdul,
    build_v_axis_kms_from_header,
)


def test_build_v_axis_velocity_blank_cunit_defaults_to_mps() -> None:
    hdr = make_spectral_header(fits_spectral_axis=3, ctype="VRAD", cunit="", crval=1000.0, cdelt=1000.0, crpix=1.0)
    v = build_v_axis_kms_from_header(hdr, 3, spectral_fits_axis=3)
    assert np.allclose(v, [1.0, 2.0, 3.0])



def test_build_v_axis_freq_with_restfreq_and_mhz() -> None:
    rest = 115271.2018e6
    hdr = make_spectral_header(
        fits_spectral_axis=3,
        ctype="FREQ",
        cunit="MHz",
        crval=115271.2018,
        cdelt=-1.0,
        crpix=1.0,
        restfreq=rest,
    )
    v = build_v_axis_kms_from_header(hdr, 2, spectral_fits_axis=3)
    expected_hz = np.array([115271.2018, 115270.2018]) * 1e6
    expected_v = (rest - expected_hz) / rest * 299792.458
    assert np.allclose(v, expected_v)



def test_build_v_axis_wavelength_is_rejected() -> None:
    hdr = make_spectral_header(fits_spectral_axis=3, ctype="WAVE", cunit="um")
    with pytest.raises(ValueError, match="wavelength-like"):
        build_v_axis_kms_from_header(hdr, 4, spectral_fits_axis=3)



def test_get_cube_from_hdul_skips_analysis_and_unusable_candidates() -> None:
    bad_mask_hdr = make_spectral_header(fits_spectral_axis=3, ctype="VRAD", cunit="m/s")
    bad_mask_hdu = make_image_hdu(make_cube(value=1.0), name="MASK3D", header=bad_mask_hdr, btype="SignalMask")

    bad_freq_hdr = make_spectral_header(fits_spectral_axis=3, ctype="FREQ", cunit="Hz", restfreq=None)
    bad_freq_hdu = make_image_hdu(make_cube(value=2.0), name="BADFREQ", header=bad_freq_hdr)

    good_hdr = make_spectral_header(fits_spectral_axis=3, ctype="VRAD", cunit="m/s")
    good_hdu = make_image_hdu(make_cube(value=3.0), name="SCI", header=good_hdr)

    hdul = fits.HDUList([fits.PrimaryHDU(), bad_mask_hdu, bad_freq_hdu, good_hdu])
    data, hdr, idx = _get_cube_from_hdul(hdul)
    assert idx == 3
    assert hdul[idx].name == "SCI"
    assert data.shape == (4, 3, 2)
    assert str(hdr["CTYPE3"]).upper() == "VRAD"



def test_get_cube_from_hdul_rejects_explicit_analysis_ext() -> None:
    hdr = make_spectral_header(fits_spectral_axis=3, ctype="VRAD", cunit="m/s")
    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        make_image_hdu(make_cube(), name="MASK3D", header=hdr, btype="SignalMask"),
    ])
    with pytest.raises(ValueError, match="analysis/product HDU"):
        _get_cube_from_hdul(hdul, cube_ext="MASK3D")



def test_get_cube_from_hdul_rejects_explicit_freq_without_restfreq() -> None:
    hdr = make_spectral_header(fits_spectral_axis=3, ctype="FREQ", cunit="MHz", restfreq=None)
    hdul = fits.HDUList([fits.PrimaryHDU(), make_image_hdu(make_cube(), name="SCI", header=hdr)])
    with pytest.raises(ValueError, match="RESTFRQ/RESTFREQ"):
        _get_cube_from_hdul(hdul, cube_ext="SCI")
