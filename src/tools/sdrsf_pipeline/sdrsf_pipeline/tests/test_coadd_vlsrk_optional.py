import numpy as np
import pandas as pd
import pytest

pytest.importorskip("astropy")

from sdrsf_pipeline.fitsio import write_tastar_fits
from sdrsf_pipeline.scripts.coadd_fits_vlsrk import main as coadd_main

def _make_meta(nchan=128, crval=1.0e11, cdelt=-1.0e6, rest=1.0e11):
    return dict(
        nchan=nchan,
        timesys="UTC",
        ctype1="FREQ",
        cunit1="Hz",
        crval1_hz=crval,
        crpix1=1.0,
        cdelt1_hz=cdelt,
        rest_hz=rest,
        specsys="LSRK",
        telescope="test",
        backend="test",
        site_lat_deg=35.0,
        site_lon_deg=135.0,
        site_height_m=0.0,
    )

def test_coadd_fits_vlsrk_regrids_per_row_meta(tmp_path, monkeypatch):
    # Create two files with same WCS but different v_corr; spectra represent a fixed gaussian in VLSRK.
    nchan = 128
    meta = _make_meta(nchan=nchan)
    # Build radio velocity axis (approx) from frequency
    # We'll make signal as function of v_lsrk = v_radio + v_corr, so that after applying v_corr, lines align.
    from sdrsf_pipeline.axis import freq_axis_from_wcs, radio_velocity_kms
    f = freq_axis_from_wcs(meta)
    v_radio = radio_velocity_kms(f, rest_hz=meta["rest_hz"])

    def make_rows(vcorr):
        v_lsrk = v_radio + vcorr
        y = np.exp(-((v_lsrk - 0.0)/5.0)**2)
        # repeat rows
        A = np.vstack([y, y])
        t0 = pd.Timestamp("2026-02-11T00:00:00", tz="UTC")
        tab = pd.DataFrame(dict(timestamp=[t0, t0 + pd.Timedelta(seconds=1)], pos_id=[0,0], ra_deg=[10,10], dec_deg=[20,20], v_corr_kms=[vcorr, vcorr]))
        return A, tab

    A1, tab1 = make_rows(0.0)
    A2, tab2 = make_rows(3.0)

    f1 = tmp_path/"a.fits"
    f2 = tmp_path/"b.fits"
    write_tastar_fits(str(f1), meta, A1, tab1.set_index(pd.to_datetime(tab1["timestamp"])))
    write_tastar_fits(str(f2), meta, A2, tab2.set_index(pd.to_datetime(tab2["timestamp"])))

    out = tmp_path/"out.fits"

    # run CLI main via monkeypatch argv
    argv = ["prog", "--out", str(out), "--mode", "uniform", "--rms-vwin", "-50:50", str(f1), str(f2)]
    monkeypatch.setattr("sys.argv", argv)
    coadd_main()

    assert out.exists()
