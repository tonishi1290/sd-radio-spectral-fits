import numpy as np

from sdrsf_pipeline.axis import slice_channels, freq_axis_from_wcs

def test_slice_channels_updates_wcs():
    meta = dict(
        nchan=10,
        timesys="UTC",
        ctype1="FREQ", cunit1="Hz",
        crval1_hz=100.0,
        crpix1=1.0,
        cdelt1_hz=-2.0,
        rest_hz=100.0,
    )
    A = np.zeros((2,10))
    f0 = freq_axis_from_wcs(meta)
    meta2, A2 = slice_channels(meta, A, 3, 8)  # keep ch=3..7 (5 channels)
    f1 = freq_axis_from_wcs(meta2)
    assert A2.shape == (2,5)
    # f1[0] should equal original f0[3]
    assert abs(float(f1[0]) - float(f0[3])) < 1e-9
