from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from sdrsf_pipeline.rawspec import build_rawspec
from sdrsf_pipeline.calibrate import make_tastar_dumps

def _has_astropy():
    try:
        import astropy  # noqa
        return True
    except Exception:
        return False

@pytest.mark.skipif(not _has_astropy(), reason="astropy not installed")
def test_vcorr_autocompute_for_vlsrk_slice():
    nchan = 64
    t0 = pd.Timestamp("2026-02-11T00:00:00", tz="UTC")
    hot = pd.DataFrame(np.full((2, nchan), 1000.0), index=pd.date_range(t0, periods=2, freq="10S"))
    off = pd.DataFrame(np.full((3, nchan),  800.0), index=pd.date_range(t0, periods=3, freq="5S"))
    on  = pd.DataFrame(np.full((4, nchan),  820.0), index=pd.date_range(t0+pd.Timedelta(seconds=6), periods=4, freq="2S"))

    meta = dict(
        nchan=nchan, timesys="UTC",
        crval1_hz=115e9, crpix1=1.0, cdelt1_hz=-1e6, rest_hz=115e9,
        site_lat_deg=35.0, site_lon_deg=135.0, site_height_m=0.0,
    )
    mapping = pd.DataFrame(dict(timestamp=on.index, pos_id=0, ra_deg=83.6, dec_deg=-5.39)).set_index("timestamp")
    raw = build_rawspec(hot=hot, on=on, off=off, meta=meta, mapping=mapping)

    res = make_tastar_dumps(raw, vlsrk_range_kms=(-50, 50))
    assert res.data.shape[0] == 4
    assert res.data.shape[1] <= nchan
    assert "v_corr_kms" in res.table.columns
