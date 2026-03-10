from __future__ import annotations

import numpy as np
from sdrsf_pipeline.regrid_vlsrk import make_vgrid, interp_to_vgrid

def test_interp_identity_increasing():
    v_src = np.arange(-10.0, 11.0, 1.0)
    y_src = v_src**2
    v_tgt = make_vgrid(-10, 10, 1.0).axis()
    y = interp_to_vgrid(v_src, y_src, v_tgt)
    assert np.allclose(y, y_src, atol=0, rtol=0)

def test_interp_identity_decreasing():
    v_src = np.arange(10.0, -11.0, -1.0)
    y_src = (v_src + 1.0)**2
    v_tgt = make_vgrid(-10, 10, 1.0).axis()
    y = interp_to_vgrid(v_src, y_src, v_tgt)
    # After reversing, should match exact analytic on integer grid
    assert np.allclose(y, (v_tgt + 1.0)**2, atol=0, rtol=0)
