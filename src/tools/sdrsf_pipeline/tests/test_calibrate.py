import numpy as np
import pandas as pd

from sdrsf_pipeline.rawspec import build_rawspec
from sdrsf_pipeline.calibrate import tastar_from_rawspec


def test_make_tastar_dumps_basic():
    # Make small synthetic rawspec with explicit HOT/OFF/ON time series.
    idx_hot = pd.date_range("2026-01-01", periods=2, freq="1s", tz="UTC")
    idx_off = pd.date_range("2026-01-01", periods=2, freq="1s", tz="UTC")
    idx_on = pd.date_range("2026-01-01", periods=3, freq="1s", tz="UTC")

    hot = pd.DataFrame(np.ones((2, 4)) * 10.0, index=idx_hot)
    off = pd.DataFrame(np.ones((2, 4)) * 2.0, index=idx_off)
    on  = pd.DataFrame(np.ones((3, 4)) * 5.0, index=idx_on)

    # mapping provided here is ON-only; build_rawspec expands it to full mapping_all.
    mapping_on = pd.DataFrame({"pos_id": [0, 1, 2], "v_corr_kms": [0.0, 1.0, 2.0]}, index=idx_on)

    meta = {"nchan": 4, "crval1_hz": 1.0, "cdelt1_hz": 1.0, "crpix1": 1.0, "restfrq_hz": 1.0}
    rawspec = build_rawspec(hot=hot, on=on, off=off, meta=meta, mapping=mapping_on)

    res = tastar_from_rawspec(rawspec, tamb_k=300.0)

    assert res.data.shape == (3, 4)
    assert len(res.table) == 3
    assert "pos_id" in res.table.columns
    assert res.table["pos_id"].tolist() == [0, 1, 2]
