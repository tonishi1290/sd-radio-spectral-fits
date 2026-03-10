from __future__ import annotations
import numpy as np
from sdrsf_pipeline.regrid_vlsrk import default_dv_kms_from_meta, vrange_from_meta_and_vcorr

def test_default_dv_positive():
    meta = dict(cdelt1_hz=-30000.0, rest_hz=115.2712018e9, crval1_hz=115.2712018e9, crpix1=1.0, nchan=10)
    dv = default_dv_kms_from_meta(meta)
    assert dv > 0

def test_vrange_from_meta_and_vcorr():
    meta = dict(cdelt1_hz=-1.0e6, rest_hz=115.0e9, crval1_hz=115.0e9, crpix1=1.0, nchan=4)
    vc = np.array([-5.0, 0.0, 7.0])
    vmin, vmax = vrange_from_meta_and_vcorr(meta, vc, nchan=4)
    assert vmax > vmin
