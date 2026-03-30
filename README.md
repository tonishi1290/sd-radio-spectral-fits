# sd-radio-spectral-fits

Single-dish radio spectral analysis package for SDFITS / Scantable workflows.

- Distribution / repository name: `sd-radio-spectral-fits`
- Python import name: `sd_radio_spectral_fits`

## Overview

`sd_radio_spectral_fits` is a practical package for single-dish radio spectroscopy.
It supports a consistent workflow from raw spectra to calibrated data products, including
SDFITS I/O, Scantable-based processing, utility operations, interactive visualization,
and 3D cube generation.

Main capabilities:

- SDFITS / Scantable read-write workflow
- Ta* calibration, baseline fitting, rest-frequency handling, velocity regrid, and velocity coadd
- Scantable utilities for inspection, selection, coordinate offsets, beam-efficiency setup, and intensity scaling
- Interactive spectral viewers in `profile_view` (`viewer`, `montage`, `grid`)
- OTF / PS mapping and 3D FITS cube generation in `map_3d`
- FFT/PLAIT basketweave for OTF family cubes
- NECST RawData tools for Sun-scan analysis and SDFITS conversion
- Single-beam and multi-beam beam-measurement workflows

## Installation

```bash
pip install -e .
```

## Recommended processing order

For ordinary spectral-line analysis, the recommended order is:

1. Read SDFITS / Scantable
2. Run Ta* calibration
3. If needed, apply multi-beam relative scaling immediately after calibration
4. If needed, apply a later global scale factor
5. Regrid to a common velocity axis
6. Run baseline fitting
7. Coadd, inspect, or map

Important distinctions:

- `apply_relative_scale(...)` and `apply_global_scale(...)` multiply the actual spectral data.
- `set_beameff(...)` only sets physical `BEAMEFF` metadata and does **not** modify spectral data.
- Physical `TR*` conversion uses `BEAMEFF` as

$$
T_R^* = \frac{T_A^*}{\eta}
$$

so generic beam-to-beam intensity alignment factors should **not** be stored in `BEAMEFF`.

## Quick start: single-beam spectral analysis

```python
import sd_radio_spectral_fits as sd

sc_raw = sd.read_scantable("example_raw.fits")

sc_cal = sd.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-200, 200),
    t_hot_k=300.0,
    vcorr_chunk_sec=5,
    dtype="float32",
)

sc_vreg = sd.run_velocity_regrid(
    input_data=sc_cal,
    vmin_kms=-30.0,
    vmax_kms=55.0,
    dv_kms=0.2,
)

sc_bsl = sd.run_baseline_fit(
    input_data=sc_vreg,
    poly_order=3,
    vwin=["-30:-10", "20:50"],
)

sd.view_spectra(sc_bsl, xrange=(-30, 50))
```

## Quick start: multi-beam scaling after calibration

```python
import sd_radio_spectral_fits as sd

sc_raw = sd.read_scantable("example_multibeam_raw.fits")

sc_cal = sd.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-35, 60),
    vcorr_chunk_sec=10,
    dtype="float32",
)

# Relative intensity alignment between beams.
# These factors are multiplied into the spectral data.
relative_map = {
    (0, 0, 0): 1.00,
    (1, 0, 0): 0.97,
    (2, 0, 0): 1.03,
}

sc_rel = sd.apply_relative_scale(
    sc_cal,
    relative_map,
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
    strict=True,
)

# Optional later global factor.
sc_scaled = sd.apply_global_scale(sc_rel, 1.05)

# Physical BEAMEFF for later TR* handling.
beameff_map = {
    (0, 0, 0): 0.46,
    (1, 0, 0): 0.45,
    (2, 0, 0): 0.44,
}

sc_ready = sd.set_beameff(
    sc_scaled,
    beameff_map,
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
    strict=True,
)
```

## OTF workflow example

Single-beam OTF processing typically follows:

```python
import sd_radio_spectral_fits as sd
import sd_radio_spectral_fits.map_3d as m3d

sc_raw = sd.read_scantable("necst_otf_raw.fits")

sc_cal = sd.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-35, 60),
    vcorr_chunk_sec=10,
    dtype="float32",
)

sc_vreg = sd.run_velocity_regrid(
    input_data=sc_cal,
    vmin_kms=-30.0,
    vmax_kms=55.0,
    dv_kms=0.2,
)

sc_bsl = sd.run_baseline_fit(
    input_data=sc_vreg,
    poly_order=1,
    vwin=["-30:0", "20:55"],
)

cfg = m3d.MapConfig(
    x0=-3600.0,
    y0=-3600.0,
    nx=60,
    ny=60,
    cell_arcsec=120.0,
    beam_fwhm_arcsec=360.0,
    kernel="sf",
    convsupport=3,
    min_abs_weight_ratio=0.25,
    min_cancel_ratio=0.35,
    estimate_effective_beam=True,
)

bundle = m3d.grid_otf_family(
    sc_bsl,
    config=cfg,
    family_label="X",
    coord_sys="icrs",
    projection="SFL",
    ref_lon=83.809,
    ref_lat=-5.372639,
    dv_kms=0.2,
    vmin_kms=-30.0,
    vmax_kms=55.0,
    otf_input_state="with_turnarounds",
)
```

For multi-beam OTF analysis, apply `apply_relative_scale(...)` immediately after calibration and before velocity regrid / baseline / OTF gridding.

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

- `docs/sdrsf_package_manual_ja_scaling_fixed.md`
- `docs/sdrsf_cookbook_ja_scaling_rechecked.md`
- `docs/sdfits_final_complete_manual_ja_scaling_fixed.md`
- `docs/sdrsf_scantable_utilities_manual_ja_rechecked.md`
- `docs/profile_view_manual_ja_v2.md`
- `docs/OTF_gridding_and_plait_fft_manual_updated_2026-03-30_calibration_scaling_added.md`
- `docs/singlebeam_sunscan_converter_manual.md`
- `docs/multibeam_beam_measurement_manual_v2.md`

Suggested reading order:

1. `sdrsf_package_manual_ja_scaling_fixed.md`
2. `sdrsf_cookbook_ja_scaling_rechecked.md`
3. `sdfits_final_complete_manual_ja_scaling_fixed.md`
4. `sdrsf_scantable_utilities_manual_ja_rechecked.md`
5. `profile_view_manual_ja_v2.md`
6. `OTF_gridding_and_plait_fft_manual_updated_2026-03-30_calibration_scaling_added.md`

## Notes

- `profile_view` is intended for interactive spectral inspection.
- `map_3d` provides the practical pipeline for OTF / PS cube making.
- For OTF analysis, calibration is part of the practical workflow even though the cube-making package itself starts from calibrated scantables.
- For multi-beam analysis, relative scaling should normally be applied immediately after calibration.
- `set_beameff(...)` is for physical beam-efficiency metadata, not for generic beam-to-beam intensity alignment.
- The Sun-scan and converter tools are primarily designed for NECST RawData.
- Multi-beam analysis is built on top of validated single-beam analysis and converter-compatible `beams.toml` settings.

## License

MIT License. See [LICENSE](LICENSE).
