from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from tests.helpers import make_cube, make_image_hdu, make_spectral_header
from sd_radio_spectral_fits.map_3d import cube_analysis as ca



def test_clean_header_removes_conflicting_and_compression_keywords() -> None:
    hdr = fits.Header()
    hdr["SIMPLE"] = True
    hdr["EXTEND"] = True
    hdr["BITPIX"] = -32
    hdr["XTENSION"] = "IMAGE"
    hdr["PCOUNT"] = 0
    hdr["GCOUNT"] = 1
    hdr["EXTNAME"] = "OLD"
    hdr["NAXIS"] = 3
    hdr["NAXIS1"] = 2
    hdr["NAXIS2"] = 3
    hdr["NAXIS3"] = 4
    hdr["CHECKSUM"] = "AAAA"
    hdr["DATASUM"] = "BBBB"
    hdr["ZIMAGE"] = True
    hdr["ZCMPTYPE"] = "RICE_1"
    hdr["ZNAXIS"] = 3
    hdr["ZNAXIS1"] = 2
    hdr["ZHECKSUM"] = "CCCC"
    hdr["ZDATASUM"] = "DDDD"
    hdr["BUNIT"] = "K"

    out = ca._clean_header_for_image_like_hdu(hdr)

    for key in (
        "SIMPLE", "EXTEND", "BITPIX", "XTENSION", "PCOUNT", "GCOUNT", "EXTNAME",
        "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "CHECKSUM", "DATASUM",
        "ZIMAGE", "ZCMPTYPE", "ZNAXIS", "ZNAXIS1", "ZHECKSUM", "ZDATASUM",
    ):
        assert key not in out
    assert out["BUNIT"] == "K"



def test_make_provisional_masks_respects_linefree_and_good_flags() -> None:
    linefree = np.array([True, True, False, False], dtype=bool)
    base_flag = np.array([[0, 2], [0, 0]], dtype=np.uint8)
    basesup, linecand = ca._make_provisional_masks((4, 2, 2), linefree, base_flag_2d=base_flag, good_flag_values=(0,))

    assert basesup.shape == (4, 2, 2)
    assert linecand.shape == (4, 2, 2)

    # bad pixel (0,1) must be blank in both provisional masks
    assert np.all(basesup[:, 0, 1] == 0)
    assert np.all(linecand[:, 0, 1] == 0)

    # good pixel (1,1): first two channels support baseline, last two are line candidates
    assert basesup[:, 1, 1].tolist() == [1, 1, 0, 0]
    assert linecand[:, 1, 1].tolist() == [0, 0, 1, 1]



def test_resolve_cube_ext_for_spectralcube_skips_analysis_products(tmp_path) -> None:
    path = tmp_path / "cube.fits"
    hdr = make_spectral_header(fits_spectral_axis=3, ctype="VRAD", cunit="m/s")

    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        make_image_hdu(make_cube(value=1.0), name="MASK3D", header=hdr, btype="SignalMask"),
        make_image_hdu(make_cube(value=2.0), name="SCI", header=hdr),
    ])
    hdul.writeto(path)

    resolved = ca._resolve_cube_ext_for_spectralcube(str(path))
    assert resolved == 2



def test_resolve_cube_ext_for_spectralcube_rejects_explicit_analysis_ext(tmp_path) -> None:
    path = tmp_path / "cube2.fits"
    hdr = make_spectral_header(fits_spectral_axis=3, ctype="VRAD", cunit="m/s")
    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        make_image_hdu(make_cube(value=1.0), name="MASK3D", header=hdr, btype="SignalMask"),
    ])
    hdul.writeto(path)

    with pytest.raises(ValueError, match="analysis/product HDU"):
        ca._resolve_cube_ext_for_spectralcube(str(path), cube_ext="MASK3D")
