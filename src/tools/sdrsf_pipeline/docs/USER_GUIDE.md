# USER GUIDE (sdrsf_pipeline)

This guide is intentionally verbose and assumes the reader may be new to radio spectral-line calibration.

## 0. Scope / Non-goals
This toolset focuses on:
- Formatting HOT/ON/OFF raw spectra into a single container (`rawspec.pkl`)
- 1-temperature chopper-wheel calibration to Ta* (dump-by-dump)
- Channel trimming by channel range and (optionally) by VLSRK velocity range
- Baseline fitting (multiple windows, polynomial) with machine-readable history stored in FITS
- Coadding spectra at the same sky position with uniform or RMS-weighted averaging
- Coadding across multiple FITS with different Doppler corrections by regridding to a common VLSRK grid
- Concatenating FITS files when spectral WCS is identical

Non-goals (for now):
- A full observatory control/log ingestion system
- Automatic detection of OFF/HOT blocks from mixed streams (you provide HOT/ON/OFF already separated)

## 1. Data model

### 1.1 Raw spectra arrays
HOT, ON, OFF must each be:
- shape: (Ndump, Nchan)
- either a NumPy ndarray or a pandas DataFrame
- strongly recommended: DataFrame with UTC DatetimeIndex (one timestamp per dump)

### 1.2 Raw meta (critical)
The spectral axis is defined using FITS-WCS compatible keywords (frequency axis):
- nchan: int
- crval1_hz: float  (Hz at reference pixel)
- crpix1: float     (reference pixel, 1-based)
- cdelt1_hz: float  (Hz/channel, may be negative)
- rest_hz: float    (RESTFRQ, Hz)
- timesys: "UTC"

**Why WCS?**
Because when you slice channels, you must update the axis. This pipeline updates WCS safely:
- nchan_new = ch_stop - ch_start
- crval1_hz_new = crval1_hz + ch_start * cdelt1_hz
- crpix1 unchanged

This prevents the common mistake “data sliced but axis not updated”.

### 1.3 Site (observatory) information
For Doppler/VLSRK correction computed from information, you need a site:
- either meta['site_name'] (Astropy site registry)
- or meta['site_lat_deg'], meta['site_lon_deg'], meta['site_height_m']

The pipeline writes site keywords to FITS headers:
- SITELAT / SITELON / SITEELEV / SITENAME

### 1.4 Mapping table (per ON dump)
To calibrate and later coadd, you should provide a mapping table aligned to ON dumps.
Recommended columns:
- ra_deg, dec_deg  (ICRS/J2000 in degrees)
- pos_id           (integer position ID; recommended)
- optional: scan_id, az_deg, el_deg, etc.

If pos_id is absent, some tools can group by RA/Dec within a tolerance.

## 2. Stage A: Build rawspec.pkl

### 2.1 Example
See: `examples/make_rawspec_example.py`

### 2.2 CLI
`sdrsf-make-rawspec --hot hot.pkl --on on.pkl --off off.pkl --meta meta.json --mapping mapping.csv --out rawspec.pkl`

Input formats:
- pkl: pandas DataFrame
- csv: must include timestamp column or be indexable

## 3. Stage B: 1-temperature chopper wheel calibration (Ta*)

### 3.1 Calibration formula
Using HOT=ambient load (Tamb) and OFF=sky:
- gain = Tamb / (HOT - OFF)
- Ta*  = (ON - OFF) * gain

Tamb default: 300 K (change with --tamb)

### 3.2 Time interpolation for OFF/HOT
For each ON dump time:
- OFF(t) is linearly interpolated from surrounding OFF dumps
- If only one side exists, nearest OFF is used
Same for HOT(t).

This matches the requirement.

### 3.3 Channel trimming
You can trim to reduce data volume:

A) by channel index:
- --ch-start 0-based inclusive
- --ch-stop  0-based exclusive

B) by VLSRK velocity range (if RA/Dec + site are available):
- --vmin / --vmax in km/s
The pipeline computes per-dump VLSRK correction using a calc_vobs-style formula and selects a single channel slice covering the requested velocity window for all dumps.

## 4. Stage C: Baseline fitting (reusable)
Script: `sdrsf-baseline-fit`

- You provide one or multiple windows in velocity:
  --vwin "vmin:vmax"  (repeatable)
- Polynomial order: --poly (0 constant, 1 linear, ...)

History storage:
- saved as JSON in FITS extension `HISTORY`
- contains: windows, polynomial order, coefficients, created time
This is machine-readable and can be loaded later.

## 5. Stage D: Coadd per position
Script: `sdrsf-coadd-points` (for a single Ta* dump FITS)

Modes:
- uniform: simple mean
- rms_weight: weight=1/rms^2

RMS definition (per requirement):
- select RMS windows (velocity)
- fit polynomial (order configurable)
- compute residuals, then RMS = sqrt(mean((resid-mean(resid))^2))

Baseline subtraction before averaging is optional:
- --baseline-vwin (repeatable)
- --baseline-poly

## 6. Stage E: Coadd across multiple FITS with different Doppler correction
Script: `sdrsf-coadd-fits-vlsrk`

This is used when the Doppler correction differs between spectra. Steps:
1) For each spectrum row, compute its VLSRK axis:
   v_lsrk(ch) = v_radio(ch) + v_corr_kms
   v_corr_kms is computed from timestamp+ra/dec+site if not present
2) Interpolate each spectrum to a common VLSRK grid:
   --vmin/--vmax/--dv define the grid
3) Baseline (optional), RMS weighting (optional), then coadd per position

Output FITS has:
- axis: VRAD (km/s) with SPECSYS=LSRK
- POINTS table with group_id, ndump, coordinates, etc.
- HISTORY JSON for reproducibility

## 7. FITS concatenation
Script: `sdrsf-concat-fits`
Only works when spectral WCS is identical across inputs (safe vstack).
If WCS differs, use regridding/coadd instead.

## 8. “How to create raw HOT/ON/OFF correctly” (checklist)
1) Ensure each dump has the correct timestamp (UTC).
2) Ensure HOT/ON/OFF are separated correctly (you already label them).
3) Ensure Nchan is consistent across all blocks.
4) Store spectral WCS correctly (crval/crpix/cdelt, rest_hz).
5) Store site info in meta (site_name or lat/lon/height).
6) Store RA/Dec for ON dumps; prefer pos_id for mapping.
7) Validate by running:
   - examples + scripts
   - pytest


## 7b. Explicit regrid-and-concat (frequency axis)
Command: `sdrsf-regrid-and-concat`

Use this when you want to *concatenate* multiple FITS but their frequency WCS is not identical.
This command **regrids every non-reference input to the FIRST input's frequency axis** (linear interpolation),
then concatenates rows.

Important:
- This changes data by interpolation. It is NOT a “pure concat”.
- The command name is explicit to avoid accidental axis modification.
- RESTFRQ mismatch is an error by default (use --allow-rest-mismatch only if you understand the consequences).

Example:
```bash
sdrsf-regrid-and-concat --out concat_on_ref.fits a.fits b.fits c.fits
```


### Default VLSRK grid behavior (coadd)
For `sdrsf-coadd-fits-vlsrk`:
- If you provide `--vmin/--vmax/--dv`, those values are used.
- If any of them are omitted:
  - `dv` defaults to the value implied by the FIRST input's frequency WCS and RESTFRQ.
  - `vmin/vmax` default to the **overlap range** across all inputs (so every input contributes over the entire grid).
The chosen grid is recorded in HISTORY JSON.


---

## 11. Command reference

### 11.1 Common options (many commands)
- `--max-dumps N`  
  Process only the first N rows (dump/point). Useful for quick trials on large data.  
  **Important:** applied *before* Doppler/regridding (requirement).

- `--ch-start i --ch-stop j`  
  Slice channels `[i, j)` (0-based).  
  **Important:** WCS is updated so the frequency axis remains correct (requirement).

- `--on-fail {exit,skip,mask}`  
  What to do when a spectrum cannot be processed:
  - `exit`: stop with an error (recommended for reproducibility / safety)
  - `skip`: drop that spectrum and continue
  - `mask`: keep output row but fill with NaNs (escape hatch)

### 11.2 Multiple baseline windows
Baseline windows (`--vwin` / `--baseline-vwin` etc.) are **repeatable**:

```bash
--vwin "-100:-50" --vwin "50:100"
```

Line-mask windows are also repeatable:

```bash
--line-vwin "-10:10" --line-vwin "120:130"
```

Internally:
- effective baseline windows = (baseline windows) minus (line windows)

### 11.3 Iterative fitting (helper option)
If `--iter-max > 0` (or `--baseline-iter-max > 0`), the fitter can sigma-clip outliers and refit.
This is intended as a *helper*; explicit line masks are safer.

If baseline windows become empty (e.g., spectrum dominated by lines), the fit becomes impossible.
Default behavior is to stop (`exit`) unless you choose `--on-fail`.

### 11.4 --block-size guideline (large data)
Approx memory for float64 spectra:
- per row ≈ `8 * Nchan` bytes
- per block ≈ `8 * Nchan * block_size` bytes

Example: `Nchan=32768`, `block_size=512` → ~134 MB.

---

## 12. concat vs regrid vs coadd
- Safe “concat” requires identical WCS.
- `sdrsf-concat-fits --regrid-mismatch-to-first` will linearly interpolate mismatched inputs to the first file’s frequency axis (data changes).
- For summation into one spectrum per position, use `sdrsf-coadd-fits-vlsrk` (VLSRK-aware regridding).

---

## 13. sdrsf-coadd-fits-vlsrk behavior
- Builds per-row `v_lsrk` axis and interpolates to a common grid.
- Default grid = **common overlap range** of inputs.
- If `--vmin/--vmax/--dv` are given, they are used; by default the command errors if you go outside the overlap.
- `--allow-outside-overlap` permits extrapolated parts to become NaN for some inputs (not recommended).



## Notes on Ndump/Nchan trimming

Most scripts support **data reduction before heavy processing**:

- `--max-dumps N` trims the first N spectra (rows).
- `--ch-start S --ch-stop E` slices channels before further steps.
  The spectral WCS is updated consistently (CRVAL1 shift + NCHAN).

For example, to reduce the size of a rawspec early:

```bash
sdrsf-make-rawspec ... --max-dumps 2000 --ch-start 1024 --ch-stop 3072
```

The same concept is available in `sdrsf-regrid-and-concat` and `sdrsf-baseline-fit-fits`.


## Baseline fitting and VLSRK windows

`sdrsf-baseline-fit-fits` accepts baseline windows in **km/s** (radio definition).
When the FITS table includes a per-row velocity correction column (default `v_corr_kms`),
the script applies it as:

- `v_lsrk(ch,row) = v_radio(ch) + v_corr_kms[row]`

This keeps velocity windows stable even when Doppler correction varies with time.

If the column is missing, the script assumes `v_corr_kms=0` and still runs (explicitly recorded in history).
