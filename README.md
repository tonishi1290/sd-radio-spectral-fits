# sd-radio-spectral-fits

Single-dish radio spectral **SDFITS-style** FITS writer.

- Distribution / repo name: `sd-radio-spectral-fits`
- Python import name: `sd_radio_spectral_fits`

## Install (editable)

```bash
pip install -e .
```

## Quickstart

```python
from sd_radio_spectral_fits import (
    Site, DatasetInfo, SpectralAxisUniform, Efficiency, SDRadioSpectralSDFITSWriter
)

site = Site(lat_deg=35.94, lon_deg=138.472, elev_m=1350.0)

axis = SpectralAxisUniform(
    crval1_hz=115.2712018e9,
    cdelt1_hz=61e3,
    crpix1=1.0,
    restfreq_hz=115.2712018e9,
    specsys="TOPOCENT",
    veldef="RADIO",
    refchan=1,
)

info = DatasetInfo(
    telescope="MyDish",
    observer="Observer",
    project="DEMO",
    object_name="Target",
    radesys="ICRS",
    equinox=2000.0,
    src_radesys="ICRS",
    src_equinox=2000.0,
    eff=Efficiency(effstat="UNKNOWN"),
    refr_included_in_corr=True,
    doppler_tracking_applied=False,
    spectral_axis=axis,
)

w = SDRadioSpectralSDFITSWriter(n_chan=1024, site=site, info=info)
# w.add_row(...)
# w.write("out.fits")
```

## Public API policy

The top-level package (`sd_radio_spectral_fits`) intentionally exports only the
main writer and a small set of dataclasses.

For advanced/low-level symbols (full schema constants, enums, helper types),
import from:

```python
from sd_radio_spectral_fits.sdfits_writer import SCHEMA
```


## Examples

Run examples after installing:

```bash
python -m examples.ex01_raw_on_off_hot
python -m examples.ex02_tastar_dump_multi_points
python -m examples.ex03_integrated_tastar_multi_points
python -m examples.ex04_integrated_doppler_corrected
python -m examples.ex05_otf_raw_hot_off_continuous_on
python -m examples.ex06_sun_scan_raw_radec_target_azel_offset
python -m examples.ex07_sun_scan_total_power_only
```


## Provenance (recommended)

This writer records the software name and version in FITS headers:

- `SWNAME = 'sd-radio-spectral-fits'`
- `SWVER  = '<package version>'`

These keywords also appear in console output when writing files.


## Header tracing (recommended)

To make future provenance tracking easier, the writer always records:

- `SWNAME = 'sd-radio-spectral-fits'`
- `SWVER  = '<package version>'`

in both the primary header and the binary-table extension header.

For spectral-frame handling, the extension header writes:

- `SPECSYS`  : spectral reference frame of the stored axis
- `SSYSOBS`  : observer frame (commonly `TOPOCENT`)
- `VELDEF`   : velocity definition (`RADIO`/`OPTI`/`RELA`)
- `VELREF`   : AIPS/casacore integer (adds `+256` when `VELDEF='RADIO'`)

