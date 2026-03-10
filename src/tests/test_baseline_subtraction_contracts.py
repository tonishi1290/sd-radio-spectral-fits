from __future__ import annotations

import numpy as np
from astropy.io import fits

from tests.helpers import make_cube, make_image_hdu, make_spectral_header
from sd_radio_spectral_fits.map_3d.baseline_subtraction import (
    BaselineConfig,
    _restore_cube_axis_order,
    _standardize_cube_for_processing,
    subtract_baseline_from_fits,
)



def test_standardize_and_restore_cube_axis_roundtrip_for_y_x_v() -> None:
    data_y_x_v = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
    hdr = make_spectral_header(fits_spectral_axis=1, ctype="VRAD", cunit="m/s")

    data_std, axis_order = _standardize_cube_for_processing(data_y_x_v, hdr)
    assert axis_order == "y_x_v"
    assert data_std.shape == (4, 2, 3)

    restored = _restore_cube_axis_order(data_std, axis_order)
    assert restored.shape == data_y_x_v.shape
    assert np.array_equal(restored, data_y_x_v)



def test_subtract_baseline_from_fits_replaces_stale_qc_hdus(tmp_path) -> None:
    input_fits = tmp_path / "in.fits"
    output_fits = tmp_path / "out.fits"

    cube = np.zeros((4, 2, 2), dtype=np.float32)
    hdr = make_spectral_header(fits_spectral_axis=3, ctype="VRAD", cunit="m/s")

    hdul = fits.HDUList([
        fits.PrimaryHDU(data=cube, header=hdr),
        make_image_hdu(np.array([1, 0, 1, 0], dtype=np.uint8), name="LINEFREE"),
        make_image_hdu(np.ones((2, 2), dtype=np.float32), name="RMS"),
        make_image_hdu(np.ones((2, 2), dtype=np.float32), name="BASE_RMS"),
        make_image_hdu(np.zeros((2, 2), dtype=np.uint8), name="BASE_FLG"),
    ])
    hdul.writeto(input_fits)

    subtract_baseline_from_fits(
        str(input_fits),
        str(output_fits),
        linefree_mask=np.ones(4, dtype=bool),
        baseline_cfg=BaselineConfig(poly_order=0, ripple=False, robust=False, chunk_pix=8),
        add_qc_hdus=True,
        overwrite=True,
    )

    with fits.open(output_fits) as hdul_out:
        names = [str(getattr(h, "name", "") or "") for h in hdul_out]
        assert names.count("LINEFREE") == 1
        assert names.count("BASE_RMS") == 1
        assert names.count("BASE_FLG") == 1
        assert "RMS" not in names
        assert hdul_out[0].data.shape == cube.shape
