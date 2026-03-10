from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from sdrsf_pipeline.rawspec import build_rawspec, save_rawspec, load_rawspec

def test_build_and_roundtrip(tmp_path):
    nchan = 64
    t0 = pd.Timestamp("2026-02-11T00:00:00", tz="UTC")
    hot = pd.DataFrame(np.ones((2, nchan)), index=pd.date_range(t0, periods=2, freq="1S"))
    on  = pd.DataFrame(np.ones((3, nchan))*2, index=pd.date_range(t0, periods=3, freq="1S"))
    off = pd.DataFrame(np.ones((4, nchan))*0.5, index=pd.date_range(t0, periods=4, freq="1S"))

    meta = dict(nchan=nchan, timesys="UTC", crval1_hz=1.0, crpix1=1.0, cdelt1_hz=1.0, rest_hz=1.0)
    mapping = pd.DataFrame(dict(timestamp=on.index, pos_id=[0,0,1])).set_index("timestamp")

    raw = build_rawspec(hot=hot, on=on, off=off, meta=meta, mapping=mapping)
    p = tmp_path / "raw.pkl"
    save_rawspec(raw, str(p))
    raw2 = load_rawspec(str(p))

    assert raw2["on"].shape == (3, nchan)
    assert raw2["meta"]["nchan"] == nchan

def test_nchan_mismatch_raises():
    nchan = 8
    t0 = pd.Timestamp("2026-02-11T00:00:00", tz="UTC")
    hot = pd.DataFrame(np.ones((2, nchan)), index=pd.date_range(t0, periods=2, freq="1S"))
    on  = pd.DataFrame(np.ones((3, nchan+1)), index=pd.date_range(t0, periods=3, freq="1S"))
    off = pd.DataFrame(np.ones((4, nchan)), index=pd.date_range(t0, periods=4, freq="1S"))
    meta = dict(nchan=nchan, timesys="UTC", crval1_hz=1.0, crpix1=1.0, cdelt1_hz=1.0, rest_hz=1.0)
    with pytest.raises(ValueError):
        build_rawspec(hot=hot, on=on, off=off, meta=meta, mapping=None)
