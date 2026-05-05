import numpy as np
import pytest

pytest.importorskip("astropy")
from astropy.io import fits

from sd_radio_spectral_fits import (
    SDRadioSpectralSDFITSWriter,
    Site,
    DatasetInfo,
    SpectralAxisUniform,
)


def test_topocentric_multiline_freq_rows_emit_restfreq_and_veldef(tmp_path):
    """Topocentric FREQ rows with RESTFREQ must carry row-level line context.

    Snapshot-driven multi-line outputs can contain different rest frequencies
    and different CRVAL1/CDELT1 per row while still using CTYPE1='FREQ' and
    SPECSYS='TOPOCENT'.  The primary-header RESTFREQ default is not enough in
    that case; the SINGLE DISH table must include RESTFREQ/VELDEF columns.
    """

    site = Site(lat_deg=35.940874, lon_deg=138.472153, elev_m=1386.0)
    info = DatasetInfo(
        telescope="OMU1P85M",
        object_name="test",
        spectral_axis=SpectralAxisUniform(
            crval1_hz=230.0e9,
            cdelt1_hz=1.0e6,
            crpix1=1.0,
            restfreq_hz=230.538e9,
            specsys="TOPOCENT",
            veldef="RADIO",
            ctype1="FREQ",
            cunit1="Hz",
        ),
    )
    writer = SDRadioSpectralSDFITSWriter(
        n_chan=4,
        site=site,
        info=info,
        store_freq_column=False,
    )

    common = dict(
        time_mjd=60000.0,
        scanid=1,
        subscan=1,
        intgrp=1,
        obsmode="ON",
        data=np.arange(4, dtype=np.float32),
        exposure_s=0.1,
        polariza="XX",
        ra_deg=83.8,
        dec_deg=-5.4,
        crpix1=1.0,
        ctype1="FREQ",
        cunit1="Hz",
        specsys="TOPOCENT",
        veldef="RADIO",
    )
    writer.add_row(
        **common,
        fdnum=1,
        ifnum=0,
        plnum=0,
        restfreq_hz=230.538e9,
        crval1_hz=230.38744773705253e9,
        cdelt1_hz=76296.27368999299,
    )
    writer.add_row(
        **common,
        fdnum=4,
        ifnum=0,
        plnum=0,
        restfreq_hz=110.201e9,
        crval1_hz=110.24657277138584e9,
        cdelt1_hz=-76296.27368999299,
    )

    out = tmp_path / "multi_line_restfreq.fits"
    writer.write(str(out), overwrite=True)

    with fits.open(out) as hdul:
        data = hdul[1].data
        names = set(hdul[1].columns.names)
        assert "RESTFREQ" in names
        assert "VELDEF" in names
        assert np.allclose(data["RESTFREQ"], [230.538e9, 110.201e9])
        # The writer normalizes the user-friendly input veldef="RADIO"
        # with specsys="TOPOCENT" to the FITS/SDFITS-style row value
        # "RADI-OBS".  The important point for multi-line output is that
        # VELDEF is present row-by-row and is not omitted for topocentric
        # FREQ axes.
        assert [v.strip() for v in data["VELDEF"]] == ["RADI-OBS", "RADI-OBS"]
        assert np.allclose(data["CRVAL1"], [230.38744773705253e9, 110.24657277138584e9])
        assert np.allclose(data["CDELT1"], [76296.27368999299, -76296.27368999299])
