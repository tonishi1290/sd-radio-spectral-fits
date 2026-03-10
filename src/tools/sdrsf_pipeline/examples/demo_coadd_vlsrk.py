#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import numpy as np
import pandas as pd

from sdrsf_pipeline.rawspec import build_rawspec
from sdrsf_pipeline.calibrate import make_tastar_dumps
from sdrsf_pipeline.fitsio import write_tastar_fits

def _gauss(x, mu, sig):
    return np.exp(-0.5*((x-mu)/sig)**2)

def main():
    # Build two rawspecs with different v_corr_kms in mapping, then calibrate and write two FITS.
    Ndump = 20
    Nchan = 2048
    t0 = pd.Timestamp("2026-02-11T00:00:00", tz="UTC")
    on_time  = pd.date_range(t0, periods=Ndump, freq="1S")
    off_time = pd.date_range(t0 - pd.Timedelta(seconds=10), periods=10, freq="2S")
    hot_time = pd.date_range(t0 - pd.Timedelta(seconds=20), periods=5,  freq="5S")

    # WCS for frequency axis (arbitrary example)
    meta = dict(
        telescope="1p85m",
        backend="xffts",
        nchan=Nchan,
        timesys="UTC",
        ctype1="FREQ",
        cunit1="Hz",
        crval1_hz=115.2712018e9,
        crpix1=1.0,
        cdelt1_hz=-30e3,
        rest_hz=115.2712018e9,
        specsys="LSRK",
    )

    # Construct synthetic spectra: (HOT, OFF constant), ON contains a line in radio-velocity space
    # We'll fake a line centered at v=10 km/s in *LSRK*; different v_corr will shift apparent frequency.
    # For demo, we do not model true Doppler; this is just to exercise the pipeline.

    rng = np.random.default_rng(0)
    hot = 1000.0 + rng.normal(scale=2.0, size=(len(hot_time), Nchan))
    off =  800.0 + rng.normal(scale=2.0, size=(len(off_time), Nchan))

    # Make ON as OFF + line (in channels)
    ch = np.arange(Nchan)
    on_base = 820.0 + rng.normal(scale=2.0, size=(Ndump, Nchan))
    # Put a line near channel 900
    on = on_base + 50.0 * _gauss(ch[None,:], 900.0, 10.0)

    hot_df = pd.DataFrame(hot, index=hot_time)
    off_df = pd.DataFrame(off, index=off_time)
    on_df  = pd.DataFrame(on,  index=on_time)

    # Two mappings with different v_corr_kms
    mapping1 = pd.DataFrame(dict(
        timestamp=on_time,
        ra_deg=83.6, dec_deg=-5.39,
        pos_id=0,
        v_corr_kms=np.full(Ndump, -5.0),
    )).set_index("timestamp")

    mapping2 = pd.DataFrame(dict(
        timestamp=on_time + pd.Timedelta(seconds=100),
        ra_deg=83.6, dec_deg=-5.39,
        pos_id=0,
        v_corr_kms=np.full(Ndump, +7.0),
    )).set_index("timestamp")

    # rawspec1/2
    raw1 = build_rawspec(hot=hot_df, on=on_df, off=off_df, meta=meta, mapping=mapping1)
    raw2 = build_rawspec(hot=hot_df, on=pd.DataFrame(on, index=mapping2.index), off=off_df, meta=meta, mapping=mapping2)

    res1 = make_tastar_dumps(raw1, tamb_k=300.0)
    res2 = make_tastar_dumps(raw2, tamb_k=300.0)

    write_tastar_fits("tastar_a.fits", meta=res1.meta, data=res1.data, table=res1.table, history={"demo":"a"})
    write_tastar_fits("tastar_b.fits", meta=res2.meta, data=res2.data, table=res2.table, history={"demo":"b"})
    print("saved: tastar_a.fits, tastar_b.fits")
    print("Now run:")
    print("  sdrsf-coadd-fits-vlsrk --out coadd_vlsrk.fits --vmin -100 --vmax 100 --dv 0.2 tastar_a.fits tastar_b.fits")

if __name__ == "__main__":
    main()
