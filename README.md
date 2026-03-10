# sd-radio-spectral-fits

Single-dish radio spectral analysis package for SDFITS/Scantable workflows.

- Distribution / repository name: `sd-radio-spectral-fits`
- Python import name: `sd_radio_spectral_fits`

## Overview

`sd_radio_spectral_fits` is a practical package for single-dish radio spectroscopy.
It supports a consistent workflow from raw spectra to calibrated data products, including
SDFITS I/O, Scantable-based processing, interactive visualization, and 3D cube generation.

Main capabilities:

- SDFITS / Scantable read-write workflow
- Ta* calibration, baseline fitting, rest-frequency handling, and velocity coadd
- Interactive spectral viewers in `profile_view` (`viewer`, `montage`, `grid`)
- OTF / PS mapping and 3D FITS cube generation in `map_3d`
- NECST RawData tools for Sun-scan analysis and SDFITS conversion
- Single-beam and multi-beam beam-measurement workflows

## Installation

```bash
pip install -e .
```

## Quick start

```python
import sd_radio_spectral_fits.calibrate as cal
import sd_radio_spectral_fits.coadd as coadd
import sd_radio_spectral_fits.baseline as bsl
import sd_radio_spectral_fits.fitsio as fitsio
import sd_radio_spectral_fits.profile_view.viewer as pv

sc_raw = fitsio.read_scantable("example_raw.fits")

sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-200, 200),
    t_hot_k=300.0,
    vcorr_chunk_sec=5,
    dtype="float32",
)

sc_scan = coadd.run_velocity_coadd(
    inputs=sc_cal,
    group_mode="scan",
    mode="rms_weight",
)

sc_bsl = bsl.run_baseline_fit(
    input_data=sc_scan,
    poly_order=3,
    vwin=["-30:-10", "20:50"],
)

pv.view_spectra(sc_bsl, xrange=(-30, 50))
```

## 3D mapping example

```python
from sd_radio_spectral_fits.map_3d.config import GridConfig
from sd_radio_spectral_fits.map_3d.gridder import run_mapping_pipeline

cfg = GridConfig(
    x0=0.0,
    y0=0.0,
    nx=121,
    ny=121,
    cell_arcsec=7.5,
    beam_fwhm_arcsec=16.0,
    kernel="gjinc",
)

run_mapping_pipeline(
    scantable=scantable,
    config=cfg,
    output_fits="otf_map.fits",
    projection="SFL",
    out_scale="TA*",
    dv_kms=0.2,
)
```

## NECST-related command-line tools

The repository also includes CLI tools for NECST RawData workflows.

### Single-beam Sun-scan

```bash
sunscan_singlebeam RAWDATA \
  --spectral-name xffts-board1 \
  --outdir out_single
```

### NECST RawData to SDFITS

```bash
necst_v4_sdfits_converter RAWDATA \
  --spectrometer-config beams.toml \
  --out output.fits
```

### Multi-beam workflow

```bash
check_spectrometer_config beams.toml --out-csv config_check.csv

sunscan_extract_multibeam RAWDATA \
  --spectrometer-config beams.toml \
  --outdir out_extract

sunscan_fit_multibeam \
  out_extract/sunscan_multibeam_scan_summary_<TAG>.csv \
  --spectrometer-config beams.toml \
  --outdir out_fit \
  --center-beam-id B00 \
  --model both
```

## Documentation

Detailed manuals are available under `docs/`.
Most of the current documentation is written in Japanese.

Recommended entry points:

- `docs/sdrsf_package_manual_ja_v2.md`
- `docs/sdrsf_cookbook_ja_v2.md`
- `docs/sdfits_final_complete_manual_ja.md`
- `docs/profile_view_manual_ja_v2.md`
- `docs/package_manual_jp_map3d_complete_rev3_2026-03-10.md`
- `docs/package_cookbook_jp_map3d_rev3_2026-03-10.md`
- `docs/singlebeam_sunscan_converter_manual.md`
- `docs/multibeam_beam_measurement_manual_v2.md`

## Notes

- `profile_view` is intended for interactive spectral inspection.
- `map_3d` provides the practical pipeline for OTF / PS cube making.
- The Sun-scan and converter tools are primarily designed for NECST RawData.
- Multi-beam analysis is built on top of validated single-beam analysis and converter-compatible `beams.toml` settings.

## License

MIT License. See [LICENSE](LICENSE).
