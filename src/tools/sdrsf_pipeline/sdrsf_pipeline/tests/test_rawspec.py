import numpy as np
import pandas as pd

from sdrsf_pipeline.rawspec import build_rawspec

def test_build_rawspec_shapes_and_keys():
    Nchan = 64
    t0 = pd.Timestamp("2026-02-11T00:00:00", tz="UTC")
    hot_t = pd.date_range(t0, periods=3, freq="5S")
    off_t = pd.date_range(t0, periods=4, freq="4S")
    on_t  = pd.date_range(t0, periods=10, freq="1S")

    hot = pd.DataFrame(np.ones((3, Nchan)), index=hot_t)
    off = pd.DataFrame(np.ones((4, Nchan))*2, index=off_t)
    on  = pd.DataFrame(np.ones((10, Nchan))*3, index=on_t)

    meta = dict(
        nchan=Nchan,
        timesys="UTC",
        ctype1="FREQ", cunit1="Hz",
        crval1_hz=1.0e11, crpix1=1.0, cdelt1_hz=-1.0e6,
        rest_hz=1.0e11,
        site_lat_deg=35.0, site_lon_deg=135.0, site_height_m=0.0,
        coord_frame="icrs",
    )
    mapping = pd.DataFrame(dict(timestamp=on_t, ra_deg=np.ones(len(on_t))*10, dec_deg=np.ones(len(on_t))*20)).set_index("timestamp")

    raw = build_rawspec(hot=hot, on=on, off=off, meta=meta, mapping=mapping)

    assert set(raw.keys()) >= {"meta","hot","on","off","mapping"}
    assert raw["hot"].shape == (3, Nchan)
    assert raw["on"].shape  == (10, Nchan)
    assert raw["off"].shape == (4, Nchan)
    assert raw["mapping"].shape[0] == 10
