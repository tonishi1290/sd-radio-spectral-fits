# Raw spectral data guide (HOT/ON/OFF)

This project assumes you start from **raw spectra per dump**:

- `HOT` : spectra while looking at the **ambient load** (Tamb)
- `ON`  : spectra while looking at the **sky/target**
- `OFF` : spectra while looking at an **OFF position** (blank sky)

Each is a 2D array-like object of shape **(Ndump, Nchan)**.
The recommended container is a pandas `DataFrame` with a `DatetimeIndex` in **UTC**.

## 1. Minimal required information

### 1.1 HOT/ON/OFF spectra tables

- Shape: `(Ndump, Nchan)`
- Index: dump timestamp in UTC (`DatetimeIndex` recommended)
- Values: raw spectrometer output (counts/power; linear units)
- Channel order: must be consistent across HOT/ON/OFF
  - If frequency decreases with increasing channel, that is OK as long as WCS uses the corresponding negative `CDELT1`.

### 1.2 Meta JSON (spectral WCS + observing basics)

The `meta.json` file supplies the spectral axis definition and some observing metadata.
At minimum:

- `rest_hz` : rest frequency (Hz)
- `crval1_hz` : reference value of frequency axis (Hz)
- `crpix1` : reference pixel (1-based, FITS convention)
- `cdelt1_hz` : frequency increment per channel (Hz/channel), may be negative
- `ctype1` : usually `"FREQ"`
- `cunit1` : usually `"Hz"`
- `timesys` : `"UTC"`
- `specsys` : `"LSRK"` (semantic label; the axis itself may still be topocentric)

Optional but strongly recommended (for Doppler correction reproducibility):

- `site_lat_deg`, `site_lon_deg`, `site_height_m` (or `site_name`)

## 2. Mapping table (per dump)

If you have sky coordinates and scanning information, provide a mapping table.
Recommended columns (one row per ON dump, aligned by timestamp):

- `timestamp` (UTC) or index as `DatetimeIndex`
- `ra_deg`, `dec_deg` (ICRS/J2000 degrees)
- `coord_frame` (string): `"icrs"`, `"fk5"`, `"galactic"` etc
- `pos_id` (optional): integer position/grouping label for coadds
- any additional engineering columns are allowed (az/el, scan_id, etc.)

The pipeline validates that `coord_frame` matches the expected frame.

## 3. Creating the rawspec container

Create a rawspec pickle (input for later steps):

```bash
sdrsf-make-rawspec --hot hot.pkl --on on.pkl --off off.pkl --meta meta.json --mapping mapping.csv --out rawspec.pkl
```

To reduce data volume **before** Doppler or calibration:

```bash
sdrsf-make-rawspec ... --max-dumps 2000 --ch-start 1024 --ch-stop 3072
```

This trims HOT/ON/OFF consistently and updates the WCS (CRVAL1 shift, NCHAN).

## 4. Common pitfalls

- **UTC timestamps**: mixing local time and UTC will break OFF interpolation and Doppler correction.
- **Mismatch in Nchan** across HOT/ON/OFF: must be identical after any slicing.
- **Wrong WCS sign**: if your channels run from high to low frequency, use a **negative** `cdelt1_hz`.
- **Missing site information**: VLSRK correction requires a site (lat/lon/height or a recognized `site_name`).
- **Frame mismatch**: if mapping says Galactic but you pass ICRS RA/Dec columns, the pipeline should error out.

