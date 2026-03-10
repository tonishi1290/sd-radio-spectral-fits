#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import json
import numpy as np
import pandas as pd

from sdrsf_pipeline.rawspec import build_rawspec, save_rawspec

def main():
    Ndump_hot, Ndump_on, Ndump_off = 5, 30, 12
    Nchan = 1024

    t0 = pd.Timestamp("2026-02-11T00:00:00", tz="UTC")
    hot_time = pd.date_range(t0, periods=Ndump_hot, freq="5S")
    on_time  = pd.date_range(t0 + pd.Timedelta(seconds=20), periods=Ndump_on,  freq="1S")
    off_time = pd.date_range(t0 + pd.Timedelta(seconds=10), periods=Ndump_off, freq="2S")

    rng = np.random.default_rng(0)
    hot = 1000.0 + 3.0 * rng.normal(size=(Ndump_hot, Nchan))
    off =  800.0 + 3.0 * rng.normal(size=(Ndump_off, Nchan))
    on  =  820.0 + 3.0 * rng.normal(size=(Ndump_on,  Nchan))

    hot_df = pd.DataFrame(hot, index=hot_time)
    off_df = pd.DataFrame(off, index=off_time)
    on_df  = pd.DataFrame(on,  index=on_time)

    # Minimal meta with FITS-like WCS for frequency axis
    meta = dict(
        telescope="1p85m",
        backend="xffts",
        nchan=Nchan,
        timesys="UTC",
        ctype1="FREQ",
        cunit1="Hz",
        crval1_hz=115.2712018e9,     # frequency at reference pixel
        crpix1=1.0,                  # 1-based reference pixel
        cdelt1_hz=-61.03515625e3,    # can be negative (LSB)
        rest_hz=115.2712018e9,
        specsys="LSRK",
        # site information (recommended for Doppler/VLSRK correction)
        # site_name="...",
        site_lat_deg=35.0,
        site_lon_deg=135.0,
        site_height_m=0.0,
    )

    # Mapping for ON dumps
    mapping = pd.DataFrame(
        dict(
            timestamp=on_time,
            ra_deg=np.full(Ndump_on, 83.6331),
            dec_deg=np.full(Ndump_on, -5.3911),
            az_deg=np.linspace(120, 121, Ndump_on),
            el_deg=np.linspace(45, 46, Ndump_on),
            pos_id=np.zeros(Ndump_on, dtype=int),
            # Example constant velocity correction (km/s) to build VLSRK axis:
            v_corr_kms=np.full(Ndump_on, 0.0),
        )
    ).set_index("timestamp")

    raw = build_rawspec(hot=hot_df, on=on_df, off=off_df, meta=meta, mapping=mapping)
    save_rawspec(raw, "rawspec_example.pkl")
    print("saved: rawspec_example.pkl")

    # Also save meta as JSON for CLI demo
    with open("meta_example.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    print("saved: meta_example.json")

if __name__ == "__main__":
    main()
