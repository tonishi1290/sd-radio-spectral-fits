# src/sd_radio_spectral_fits/doppler.py
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd

# NOTE: astropy is an optional dependency for this module.
# Doppler correction features require the package to be installed with the [fits] extra.

# 光速定数
from .constants import C_KMS, C_MS


def _require_astropy():
    try:
        import astropy.units as u
        from astropy.time import Time
        from astropy.coordinates import SkyCoord, EarthLocation
        import astropy.coordinates as coord
        return u, Time, SkyCoord, EarthLocation, coord
    except Exception as e:
        raise ImportError("Doppler/VLSRK correction requires astropy. Install: pip install 'sdrsf-pipeline[fits]'") from e


def earth_location_from_meta(meta: dict):
    """Build astropy.coordinates.EarthLocation from meta."""
    from astropy.coordinates import EarthLocation
    import astropy.units as u

    def _get_any(d, *keys):
        for k in keys:
            if k in d and d[k] not in (None, ""):
                return d[k]
        return None

    # 1) Geocentric XYZ
    x = _get_any(meta, "OBSGEO-X", "OBSGEO_X", "obsgeo_x_m", "site_x_m")
    y = _get_any(meta, "OBSGEO-Y", "OBSGEO_Y", "obsgeo_y_m", "site_y_m")
    z = _get_any(meta, "OBSGEO-Z", "OBSGEO_Z", "obsgeo_z_m", "site_z_m")
    if x is not None and y is not None and z is not None:
        return EarthLocation.from_geocentric(float(x) * u.m, float(y) * u.m, float(z) * u.m)

    # 2) SDFITS lat/lon/height keys
    lat = _get_any(meta, "SITELAT", "SITE_LAT", "OBS_LAT", "LAT", "site_lat_deg")
    lon = _get_any(meta, "SITELONG", "SITELON", "SITE_LON", "OBS_LON", "LON", "site_lon_deg")
    hgt = _get_any(meta, "SITEELEV", "SITEHGT", "SITE_ELEV", "OBS_HGT", "HEIGHT", "site_height_m")
    if lat is not None and lon is not None:
        if hgt is None:
            hgt = 0.0
        return EarthLocation.from_geodetic(float(lon) * u.deg, float(lat) * u.deg, float(hgt) * u.m)

    # 3) Named site
    site_name = _get_any(meta, "site_name", "SITENAME", "OBS_SITE", "SITE")
    if site_name is not None:
        return EarthLocation.of_site(str(site_name).strip())

    raise ValueError(
        "Site information missing: provide meta['site_name'] or GEO/SITE coordinates."
    )


def calc_vlsrk_correction_kms(
    *,
    ra_deg: float,
    dec_deg: float,
    unixtime: float,
    meta: dict,
    coord_frame: str = "icrs",
) -> float:
    """
    Compute VLSRK correction (km/s).
    Returns the projection of Observer velocity vector onto the Line-of-Sight.
    
    Sign Convention:
      Positive (+) : Observer is moving TOWARDS the target (Blue shift).
                     Observed Freq > Rest Freq.
                     Correction requires LOWERING the frequency (divide by k).
    """
    u, Time, SkyCoord, EarthLocation, coord = _require_astropy()
    loc = earth_location_from_meta(meta)
    
    # Time handling robust to float/int
    tobs = Time(float(unixtime), format="unix", scale="utc")

    target = SkyCoord(ra_deg*u.deg, dec_deg*u.deg, frame=coord_frame, obstime=tobs, location=loc)
    vobs_vec = SkyCoord(loc.get_gcrs(tobs)).transform_to(coord.LSRK()).velocity

    ra_rad = target.icrs.ra.to_value(u.rad)
    dec_rad = target.icrs.dec.to_value(u.rad)
    
    vx = vobs_vec.d_x.to_value(u.km/u.s)
    vy = vobs_vec.d_y.to_value(u.km/u.s)
    vz = vobs_vec.d_z.to_value(u.km/u.s)

    v_proj = vx*np.cos(dec_rad)*np.cos(ra_rad) + vy*np.cos(dec_rad)*np.sin(ra_rad) + vz*np.sin(dec_rad)
    
    return float(v_proj)


def compute_vcorr_series(
    *,
    times: pd.DatetimeIndex,
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    meta: dict,
    coord_frame: str = "icrs",
) -> np.ndarray:
    """Vectorized wrapper computing v_corr_kms per row."""
    if not isinstance(times, pd.DatetimeIndex):
        raise ValueError("times must be pandas.DatetimeIndex")

    if times.tz is None:
        t = times.tz_localize("UTC")
    else:
        t = times.tz_convert("UTC")

    t64 = t.tz_localize(None).to_numpy()
    epoch = np.datetime64("1970-01-01T00:00:00")
    unixt = (t64 - epoch) / np.timedelta64(1, "s")

    ra = np.asarray(ra_deg, float)
    dec = np.asarray(dec_deg, float)
    
    out = np.empty(unixt.shape[0], dtype=float)
    for i in range(out.size):
        out[i] = calc_vlsrk_correction_kms(
            ra_deg=float(ra[i]),
            dec_deg=float(dec[i]),
            unixtime=float(unixt[i]),
            meta=meta,
            coord_frame=coord_frame,
        )
    return out


# =========================================================
# Frequency/Velocity Frame Correction Logic
# =========================================================

def get_doppler_factor(v_kms: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    相対論的ドップラー係数 k を計算します（ベクトル化対応版）。
    k = sqrt((1 + beta) / (1 - beta))
    """
    v_kms = np.asarray(v_kms, dtype=float)
    beta = v_kms / C_KMS
    
    # 配列の中に光速を超えているものが1つでもあるかチェック
    if np.any(np.abs(beta) >= 1.0):
        # どの値が問題か特定してエラーメッセージを出す
        bad_val = v_kms[np.abs(beta) >= 1.0]
        raise ValueError(f"Velocity {bad_val} km/s exceeds speed of light!")
    
    # 計算自体はNumPyのユニバーサル関数なので、配列のまま一括で通ります
    k = np.sqrt((1.0 + beta) / (1.0 - beta))
    return k
    

def _normalize_unit(unit_str: str) -> str:
    """Normalize unit string to handle FITS variants."""
    u = str(unit_str).strip().lower()
    
    # Hz variants (e.g. Hz, KHZ, MHZ, GHZ) -> treat as frequency
    if "hz" in u:
        return "hz"
        
    # m/s variants
    if u in ("m/s", "m s-1", "ms-1", "meter/sec", "m/sec"):
        return "m/s"
        
    # km/s variants
    if u in ("km/s", "km s-1", "kms-1", "kilometer/sec", "km/sec"):
        return "km/s"
        
    # fallback
    return u


def _get_restfreq_or_raise(meta: dict) -> float:
    for k in ("RESTFRQ", "RESTFREQ", "rest_hz"):
        if k in meta and meta[k] not in (None, ""):
            return float(meta[k])
    raise ValueError(
        "Velocity axis detected but RESTFREQ/RESTFRQ missing in header. "
        "Cannot convert to Frequency for Doppler correction."
    )


def scale_frequency_wcs_by_velocity(meta: dict, v_corr_kms: float) -> dict:
    """
    Update WCS metadata to apply relativistic Doppler correction (Observer -> LSRK).
    
    Sign Conventions:
      - v_corr_kms: Internal definition (Positive = Approaching/Blue-shift).
      - FITS Velocity (VRAD/VOPT): Standard definition (Positive = Receding/Red-shift).
      - Correction Logic: f_lsrk = f_obs / k  (k_eff = 1/k)
      
    Args:
        meta: FITS header dictionary.
        v_corr_kms: Observer velocity [km/s] (positive = towards source).
        
    Returns:
        Updated meta dictionary (copy).
    """
    m = dict(meta)
    
    # 0. Check for negligible velocity (avoid float comparison k == 1.0)
    if abs(v_corr_kms) < 1e-9:
        return m
    
    k = get_doppler_factor(v_corr_kms)
    k_eff = 1.0 / k

    # 1. Parse Axis Type and Units
    raw_ctype = str(m.get("CTYPE1", "FREQ")).strip()
    ctype = raw_ctype.upper()
    
    raw_cunit = str(m.get("CUNIT1", "")) 
    unit_norm = _normalize_unit(raw_cunit)
    
    # Explicit Velocity Types (check prefix to handle -LSR, -HEL etc.)
    is_radio = ctype.startswith("VRAD")
    is_optical = ctype.startswith("VOPT")
    is_velo_rel = ctype.startswith("VELO")
    is_velo_any = is_radio or is_optical or is_velo_rel
    
    # --- [修正] 単位とCTYPEの矛盾チェック ---
    if unit_norm == "hz":
        if is_velo_any:
            # "Hz" unit but Velocity CTYPE -> Dangerous contradiction
            raise ValueError(
                f"Conflicting metadata: CUNIT1='{raw_cunit}' (Hz) but CTYPE1='{raw_ctype}' (Velocity). "
                "Cannot proceed safely."
            )
        if not ctype.startswith("FREQ"):
            # "Hz" unit but unknown CTYPE -> Warn but treat as FREQ
            warnings.warn(
                f"CUNIT1='{raw_cunit}' implies Frequency, but CTYPE1='{raw_ctype}' is not standard 'FREQ...'. "
                "Treating as Frequency axis."
            )
        is_freq = True
    else:
        # Not Hz: trust CTYPE
        is_freq = ctype.startswith("FREQ")
    # ----------------------------------------
    
    # Determine local speed of light based on CUNIT1
    if unit_norm == "m/s":
        c_local = C_MS
    elif unit_norm == "km/s":
        c_local = C_KMS
    else:
        c_local = C_KMS
    
    if is_velo_any:
        if unit_norm not in ("km/s", "m/s"):
            raise ValueError(
                f"Velocity WCS (CTYPE1='{raw_ctype}') requires explicit valid CUNIT1 ('km/s' or 'm/s'). "
                f"Found: '{raw_cunit}'"
            )
    
    # 2. Extract WCS
    try:
        crval = float(m["CRVAL1"])
        cdelt = float(m["CDELT1"])
    except KeyError:
        return m

    # 3. Apply Correction
    
    new_cdelt = cdelt # Placeholder
    
    if is_freq:
        m["CRVAL1"] = crval * k_eff
        m["CDELT1"] = cdelt * k_eff
        new_cdelt = m["CDELT1"]
        
    elif is_radio:
        # Radio Definition (VRAD): V = c * (1 - f / f0)
        # V_new = V_old * k_eff + c*(1 - k_eff)
        m["CRVAL1"] = crval * k_eff + c_local * (1.0 - k_eff)
        m["CDELT1"] = cdelt * k_eff
        new_cdelt = m["CDELT1"]

    elif is_optical:
        # Optical Definition (VOPT): V = c * (f0 / f - 1)
        # V_new = V_old / k_eff + c * (1/k_eff - 1)
        m["CRVAL1"] = crval / k_eff + c_local * (1.0 / k_eff - 1.0)
        m["CDELT1"] = cdelt / k_eff
        new_cdelt = m["CDELT1"]

    elif is_velo_rel:
        raise ValueError(
            f"CTYPE1='{raw_ctype}' (Relativistic) cannot be Doppler corrected while maintaining linear WCS. "
            "Please regrid the input data to 'VRAD' or 'FREQ' before processing."
        )
    else:
        # Fallback / Unsupported
        raise ValueError(f"Unsupported CTYPE1={raw_ctype} ({raw_cunit}) for rigorous Doppler correction.")

    # --- [修正] 軸方向の反転チェック ---
    # 通常のドップラー補正(k ~ 1)で符号が反転することは物理的に稀だが、
    # VOPTなど直感的でない軸での事故を防ぐためチェック
    if np.sign(new_cdelt) != np.sign(cdelt):
        warnings.warn(
            f"WCS CDELT1 sign flipped during Doppler correction (CTYPE1='{raw_ctype}'). "
            "Please verify the output axis direction."
        )

    return m
