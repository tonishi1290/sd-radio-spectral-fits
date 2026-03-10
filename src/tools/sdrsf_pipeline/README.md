# sdrsf_pipeline (standalone tools)

This folder provides standalone scripts for:
- Building raw spectra container (HOT/ON/OFF) and saving it
- Calibrating dump-by-dump spectra to Ta* with 1-temperature chopper wheel
- Optional channel trimming by channel range or (if available) VLSRK range
- Coadding ON dumps per position with uniform or RMS-weighted averaging
- Baseline fitting (multi-window, polynomial) with machine-readable history
- FITS concatenation and optional velocity regridding for coaddition

Important:
- This does NOT modify the main package `sd_radio_spectral_fits`.
- Install this subpackage separately (recommended):

  ```bash
  pip install -e tools/sdrsf_pipeline
  ```

Quick start (synthetic demo):
  ```bash
  python tools/sdrsf_pipeline/examples/make_rawspec_example.py
  sdrsf-make-tastar --in rawspec_example.pkl --out tastar_dumps.fits
  sdrsf-coadd-points --in tastar_dumps.fits --out coadd_points.fits --mode uniform
  ```

Notes on spectral axis:
- Use FITS WCS-compatible keywords in raw meta: CTYPE1/CUNIT1/CRVAL1/CRPIX1/CDELT1.
- CDELT1 may be negative (e.g., LSB). That is supported.

## VLSRK regrid + coadd across FITS files


Use when different spectra have different Doppler corrections and you need to coadd on a common VLSRK grid.

Required DUMPS columns in inputs:
  - v_corr_kms (or set --v-corr-col to your column name)
  - pos_id (recommended) OR ra_deg/dec_deg + --pos-tol-arcsec

Example:
  ```bash
  sdrsf-coadd-fits-vlsrk --out coadd_vlsrk.fits \
    --vmin -100 --vmax 100 --dv 0.1 \
    --mode rms_weight \
    --rms-vwin "-100:-50" --rms-vwin "50:100" --rms-poly 1 \
    --baseline-vwin "-100:-50" --baseline-vwin "50:100" --baseline-poly 1 \
    tastar_*.fits
  ```


## Documentation
- Detailed English guide: `docs/USER_GUIDE.md`
- 詳細日本語ガイド: `docs/USER_GUIDE_ja.md`

## Notes about site (observatory) information
These tools prefer to read site information from the data itself (meta / FITS header).
The CLI does not require --site-lat/--site-lon/--site-height.
Store site keywords in:
- raw meta: site_lat_deg/site_lon_deg/site_height_m or site_name
- FITS header written by this pipeline: SITELAT/SITELON/SITEELEV/SITENAME


## Japanese docs
- Top-level Japanese README: `README_ja.md`
- Detailed Japanese guide: `docs/USER_GUIDE_ja.md`


### Explicit regrid then concat
- `sdrsf-concat-fits` : strict (WCS must be identical), no interpolation
- `sdrsf-regrid-and-concat` : explicit interpolation to the first file's frequency axis, then concat


## 0.1.2 updates

- Robust handling of missing/None site metadata (no accidental float(None) crashes).
- `sdrsf-make-rawspec` supports `--max-dumps` and channel slicing (`--ch-start/--ch-stop`) with consistent WCS updates.
- `sdrsf-regrid-and-concat` supports `--max-dumps` and pre-slicing channels.
- `sdrsf-baseline-fit-fits` supports per-row `v_corr_kms` (configurable via `--v-corr-col`) when selecting velocity windows.
