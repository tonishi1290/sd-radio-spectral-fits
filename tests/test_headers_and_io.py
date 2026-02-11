import numpy as np
from astropy.io import fits

from sd_radio_spectral_fits import (
    SWNAME,
    __version__,
    SDRadioSpectralSDFITSWriter,
    Site,
    DatasetInfo,
    SpectralAxisUniform,
    Efficiency,
)


def _make_writer(*, tmp_path, chunk_size=None, store_freq_column=False, specsys="TOPOCENT", ssysobs="TOPOCENT", veldef="RADIO"):
    site = Site(lat_deg=35.94, lon_deg=138.472, elev_m=1350.0)
    axis = SpectralAxisUniform(
        crval1_hz=115.2712018e9,
        cdelt1_hz=61e3,
        crpix1=1.0,
        restfreq_hz=115.2712018e9,
        specsys=specsys,
        ssysobs=ssysobs,
        veldef=veldef,
        refchan=1,
    )
    info = DatasetInfo(
        telescope="TESTSCOPE",
        observer="TESTER",
        project="TESTPROJ",
        object_name="TESTOBJ",
        radesys="ICRS",
        equinox=2000.0,
        src_radesys="ICRS",
        src_equinox=2000.0,
        eff=Efficiency(beameff=np.nan, apereff=np.nan, effstat="UNKNOWN"),
        refr_included_in_corr=True,
        doppler_tracking_applied=False,
        spectral_axis=axis,
        shared_meta={"RECEIVER": "RX", "BACKEND": "SPEC"},
    )
    w = SDRadioSpectralSDFITSWriter(
        n_chan=16,
        site=site,
        info=info,
        store_freq_column=store_freq_column,
        duplicate_data_columns=True,
        chunk_size=chunk_size,
        out_basename=str(tmp_path / "otf_test") if chunk_size is not None else "unused",
    )
    return w


def test_swname_swver_in_primary_and_extension_header(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    w = _make_writer(tmp_path=tmp_path, chunk_size=None, store_freq_column=False)

    mjd0 = 60200.0
    for i in range(3):
        spec = np.random.normal(0, 1, 16).astype(np.float32)
        w.add_row(
            time_mjd=mjd0 + i * (0.1 / 86400.0),
            scanid=1, subscan=0, intgrp=0,
            obsmode="ON",
            data=spec,
            exposure_s=0.1,
            calstat="RAW",
            tsys_k=np.nan,
            tau=np.nan,
            ra_deg=83.82208, dec_deg=-5.39111,
            srcframe="RADEC", src_long_deg=83.82208, src_lat_deg=-5.39111,
            scanframe="RADEC", scan_radesys="ICRS", scan_equinox=2000.0,
            scan_x_deg=0.0, scan_y_deg=0.0,
            az_enc_deg=120.0, el_enc_deg=45.0,
        )

    out = tmp_path / "small.fits"
    w.write(str(out), overwrite=True)

    with fits.open(out) as hdul:
        pri = hdul[0].header
        ext = hdul[1].header

        assert pri["SWNAME"] == SWNAME
        assert pri["SWVER"] == __version__

        assert ext["SWNAME"] == SWNAME
        assert ext["SWVER"] == __version__




def test_specsys_ssysobs_velref_written(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Doppler-corrected axis example: axis is LSRK, observation is TOPOCENT (common convention)
    w = _make_writer(
        tmp_path=tmp_path,
        chunk_size=None,
        store_freq_column=True,
        specsys="LSRK",
        ssysobs="TOPOCENT",
        veldef="RADIO",
    )

    mjd0 = 60200.0
    spec = np.random.normal(0, 1, 16).astype(np.float32)
    w.add_row(
        time_mjd=mjd0,
        scanid=1, subscan=0, intgrp=0,
        obsmode="ON",
        data=spec,
        exposure_s=0.1,
        calstat="RAW",
        ra_deg=83.0, dec_deg=-5.0,
        srcframe="RADEC", src_long_deg=83.0, src_lat_deg=-5.0,
        scanframe="RADEC", scan_radesys="ICRS", scan_equinox=2000.0,
        scan_x_deg=0.0, scan_y_deg=0.0,
        az_enc_deg=120.0, el_enc_deg=45.0,
    )

    out = tmp_path / "specsys.fits"
    w.write(str(out), overwrite=True)

    with fits.open(out) as hdul:
        ext = hdul[1].header
        assert ext["SPECSYS"] == "LSRK"
        assert ext["SSYSOBS"] == "TOPOCENT"
        assert ext["VELDEF"] == "RADIO"
        assert ext["VELREF"] == 257  # 1 (LSRK) + 256 (RADIO)

def test_parts_mode_writes_manifest_and_parts(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    w = _make_writer(tmp_path=tmp_path, chunk_size=2, store_freq_column=False)

    mjd0 = 60200.0
    for i in range(5):
        spec = np.random.normal(0, 1, 16).astype(np.float32)
        w.add_row(
            time_mjd=mjd0 + i * (0.1 / 86400.0),
            scanid=1, subscan=0, intgrp=0,
            obsmode="ON",
            data=spec,
            exposure_s=0.1,
            calstat="RAW",
            ra_deg=83.0, dec_deg=-5.0,
            srcframe="RADEC", src_long_deg=83.0, src_lat_deg=-5.0,
            scanframe="RADEC", scan_radesys="ICRS", scan_equinox=2000.0,
            scan_x_deg=0.0, scan_y_deg=0.0,
            az_enc_deg=120.0, el_enc_deg=45.0,
        )

    w.close()

    # parts and manifest
    assert (tmp_path / "otf_test_part0001.fits").exists()
    assert (tmp_path / "otf_test_part0002.fits").exists()
    assert (tmp_path / "otf_test_manifest.json").exists()

    # Check headers in a part
    with fits.open(tmp_path / "otf_test_part0001.fits") as hdul:
        assert hdul[0].header["SWNAME"] == SWNAME
        assert hdul[0].header["SWVER"] == __version__
