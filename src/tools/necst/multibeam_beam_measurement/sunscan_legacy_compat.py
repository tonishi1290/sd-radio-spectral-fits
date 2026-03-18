#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sun scan reboot (NECST v4): selectable Az/El source

Why
---
You asked: "What if we interpolate altaz.lon/lat to t_spec? Encoder may have timestamp error."
This script lets you switch the Az/El source:

- encoder (recommended when timestamps are reliable): Az_true = Az_enc - dlon, El_true = El_enc - dlat
- altaz  (command stream): Az_true/El_true from altaz.lon/lat (optionally +/- dlon/dlat)

Common (always)
---------------
- t_spec from necstdb spectral table (direct read; no nercst dependency)
- tp1 is raw integrated total power (sum over channels)
- Sun center from astropy AltAz using config.toml in RawData
- Offsets:
    daz = wrap(Az_true - Az_sun) * cos(El_sun)
    d_el = El_true - El_sun

No binning. No smoothing by default. Notch ripple removal is enabled by default; disable with --ripple-no-remove.

Derivative
----------
Finite difference dy/dx in the given order.
Because daz/d_el can repeat (quantization/dwell), derivative skips duplicates by
looking for the nearest non-duplicate neighbors.

Trimming (recommended)
----------------------
To avoid dwell/pre-scan points, a motion-based trim can be enabled to keep the
largest contiguous "moving" block using vabs = |Δx/Δt|.

Outputs
-------
- AZ_dxdy_*.png / EL_dxdy_*.png                 (Sun offsets vs time; deltas; speeds)
- AZ_mount_azel_*.png / EL_mount_azel_*.png     (Az_true/El_true vs time; deltas; speeds)
- sun_scan_derivative_fits_*_pXXX.png              (all AZ/EL scan pairs; profile + derivative + fits)

Usage
-----
python sun_scan_v4_reboot_azelsource.py <RawDataDir> [--spectral xffts-board1] [--outdir .]
  --azel-source encoder|altaz
  --altaz-apply none|minus|plus
  --encoder-shift-sec <sec>
  --no-trim-scan
  --trim-vfrac 0.2 --trim-vmin 1e-4 --trim-gap 5 --trim-min-samples 100
  --no-strict-deriv
  --continue-on-error
"""

from __future__ import annotations

import argparse
import builtins
import os
import pathlib
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import math
import scipy.optimize
from scipy import signal

import necstdb

from astropy.time import Time
from astropy.coordinates import get_body, AltAz
from astropy import units as u
from astropy import constants as const

from neclib import config



def _decode_timestamp_text(v):
    if isinstance(v, (bytes, bytearray)):
        try:
            s = v.decode(errors="ignore")
        except Exception:
            s = str(v)
    else:
        s = str(v)
    return s.strip().strip("\x00")


def _find_first_nonempty_timestamp_text(values):
    arr = np.asarray(values, dtype=object)
    for v in arr:
        s = _decode_timestamp_text(v)
        if s and s.lower() != "nan":
            return s
    return None


def _extract_timestamp_suffix(text):
    s = _decode_timestamp_text(text)
    if (not s) or (s.lower() == "nan"):
        return None
    s_up = s.upper()
    if s_up.endswith("UTC"):
        return "UTC"
    if s_up.endswith("GPS"):
        return "GPS"
    if s_up.endswith("TAI"):
        return "TAI"
    if s_up.endswith("PC"):
        return "PC"
    return None


def _gps_timestamp_texts_to_unix(values):
    arr = np.asarray(values, dtype=object)
    epoch = datetime(1980, 1, 6, 0, 0, 0)
    gps_seconds = []
    for v in arr:
        s = _decode_timestamp_text(v)
        if (not s) or (s.lower() == "nan") or (len(s) < 3):
            raise ValueError(f"invalid GPS timestamp: {v!r}")
        body = s[:-3]
        dt = datetime.strptime(body, "%Y-%m-%dT%H:%M:%S.%f")
        gps_seconds.append((dt - epoch).total_seconds())
    return np.asarray(Time(np.asarray(gps_seconds, dtype=float), format="gps").utc.unix, dtype=float)


def _timestamp_texts_to_unix(values, suffix):
    arr = np.asarray(values, dtype=object)
    suf = str(suffix).strip().upper()
    texts = [_decode_timestamp_text(v) for v in arr]
    if suf == "UTC":
        bodies = [s[:-3] for s in texts]
        return np.asarray(Time(bodies, format="isot", scale="utc").unix, dtype=float)
    if suf == "TAI":
        bodies = [s[:-3] for s in texts]
        return np.asarray(Time(bodies, format="isot", scale="tai").utc.unix, dtype=float)
    if suf == "GPS":
        return _gps_timestamp_texts_to_unix(texts)
    raise ValueError(f"unsupported timestamp suffix: {suffix!r}")


def _select_spectral_time_from_structured(arr, *, numeric_fallback_candidates=("time", "t", "unix_time", "unixtime")):
    names = list(arr.dtype.names or [])
    timestamp_field = None
    for k in ("timestamp", "time_spectrometer"):
        if k in names:
            timestamp_field = k
            break

    fallback_field = None
    for k in numeric_fallback_candidates:
        if k in names:
            fallback_field = k
            break

    time_meta = {
        "applied": None,
        "timestamp_field": timestamp_field,
        "fallback_field": fallback_field,
        "suffix": None,
        "first_timestamp_text": None,
        "reason": None,
    }

    if timestamp_field is not None:
        first_text = _find_first_nonempty_timestamp_text(arr[timestamp_field])
        time_meta["first_timestamp_text"] = first_text
        suffix = _extract_timestamp_suffix(first_text)
        time_meta["suffix"] = suffix
        if suffix in ("UTC", "GPS", "TAI"):
            try:
                t_spec = _timestamp_texts_to_unix(arr[timestamp_field], suffix)
                if np.all(np.isfinite(t_spec)):
                    time_meta["applied"] = f"{timestamp_field}:{suffix}->unix"
                    time_meta["reason"] = "timestamp suffix accepted"
                    return np.asarray(t_spec, dtype=float), time_meta
                time_meta["reason"] = f"non-finite values after {suffix} conversion"
            except Exception as e:
                time_meta["reason"] = f"{suffix} conversion failed: {e}"
        elif suffix == "PC":
            time_meta["reason"] = "timestamp suffix is PC; falling back to numeric time"
        elif first_text is None:
            time_meta["reason"] = "timestamp field is empty/NaN; falling back to numeric time"
        else:
            time_meta["reason"] = f"unknown timestamp suffix in {first_text!r}; falling back to numeric time"

    if fallback_field is not None:
        t_spec = np.asarray(arr[fallback_field], dtype=float)
        time_meta["applied"] = fallback_field
        if time_meta["reason"] is None:
            time_meta["reason"] = "timestamp field unavailable; using numeric time"
        return t_spec, time_meta

    if timestamp_field is not None:
        raise RuntimeError(
            "timestamp-like field exists but no usable numeric fallback field is available; "
            f"timestamp_field={timestamp_field} available={names} reason={time_meta['reason']}"
        )

    raise RuntimeError("no usable spectral time field. available={}".format(names))


TELESCOPE = "OMU1P85M"
TEL_LOADDATA = "OMU1p85m"
DEFAULT_SPECTRAL_NAME = "xffts-board1"
DEFAULT_DB_NAMESPACE = "necst"
PLANET = "sun"




# -----------------------------
# Default parameters (CLI defaults)
# -----------------------------
# NOTE: These constants define the default values for argparse options.
#       Changing them changes script behavior, so keep them as the single source of truth.

DEFAULT_OUTDIR = "."

# Az/El source handling
DEFAULT_AZEL_SOURCE = "encoder"          # encoder|altaz
DEFAULT_ALTAZ_APPLY = "none"             # none|minus|plus (only for --azel-source altaz)
DEFAULT_ENCODER_SHIFT_SEC = 0.0
DEFAULT_ENCODER_VAVG_SEC = 0.0           # 0 disables smoothing (safe default)

# Chopper-wheel (1-load) Ta* calibration
DEFAULT_CHOPPER_WHEEL = True
DEFAULT_TAMB_FALLBACK_K = 300.0
DEFAULT_CHOPPER_WIN_SEC = 5.0            # deprecated; kept for backward compatibility
DEFAULT_CHOPPER_STAT = "median"          # deprecated; kept for backward compatibility

# Ripple (notch) removal
DEFAULT_RIPPLE_REMOVE = True
DEFAULT_RIPPLE_MODEL = "auto"            # auto|add|mul
DEFAULT_RIPPLE_TARGET_HZ = 1.2
DEFAULT_RIPPLE_SEARCH_HZ = 0.3
DEFAULT_RIPPLE_BW_HZ = 0.7
DEFAULT_RIPPLE_MAX_HARM = 8
DEFAULT_RIPPLE_ORDER = 2
DEFAULT_RIPPLE_NOTCH_PASS = 4
DEFAULT_RIPPLE_TREND_WIN_SEC = 5.0
DEFAULT_RIPPLE_RESAMPLE_DT_SEC = 0.0     # 0 means auto (median diff(t))

DEFAULT_RIPPLE_PRESET = "auto"            # auto|safe|normal|strong (preset + manual overrides)

# Ripple preset table:
# - Preset provides a *baseline* of notch settings.
# - Any explicit --ripple-* option overrides the preset value.
# - Nyquist-aware clipping is applied internally per scan, so overly large max_harm
#   will be safely reduced when dt is coarse or bw is wide.
RIPPLE_PRESETS = {
    # Signal-protective preset (good for fast scans): narrower notch, fewer passes/harmonics.
    "safe": {
        "bw_hz": 0.4,
        "max_harm": 6,
        "notch_passes": 2,
        "trend_win_sec": 4.0,
        "order": 2,
        "resample_dt_sec": 0.0,  # 0 => auto (median diff(t))
    },
    # Default preset (keeps DEFAULT_RIPPLE_* behavior).
    "normal": {
        "bw_hz": DEFAULT_RIPPLE_BW_HZ,
        "max_harm": DEFAULT_RIPPLE_MAX_HARM,
        "notch_passes": DEFAULT_RIPPLE_NOTCH_PASS,
        "trend_win_sec": DEFAULT_RIPPLE_TREND_WIN_SEC,
        "order": DEFAULT_RIPPLE_ORDER,
        "resample_dt_sec": DEFAULT_RIPPLE_RESAMPLE_DT_SEC,
    },
    # Aggressive preset (good for slow scans): wider notch and deeper attenuation.
    "strong": {
        "bw_hz": 0.8,
        "max_harm": 12,
        "notch_passes": 5,
        "trend_win_sec": 6.0,
        "order": 2,
        "resample_dt_sec": 0.0,  # 0 => auto (median diff(t))
    },
}


# Plot overlays
DEFAULT_OVERLAY_RAW_ALPHA = 0.25        # transparency for raw when overplotting raw vs corrected


# Profile x-range (offset axis)
DEFAULT_PROFILE_XLIM_DEG = 1.0        # plot/profile range [-xlim, +xlim] in degrees

# Debug plots (dxdy / mount_azel / summary_text)
DEFAULT_DEBUG_PLOT = False

DEFAULT_RIPPLE_EVAL_BAND_HZ = 0.0        # 0 means use DEFAULT_RIPPLE_BW_HZ

# Scan trimming
DEFAULT_TRIM_SCAN = True
DEFAULT_TRIM_VFRAC = 0.20
DEFAULT_TRIM_VMIN = 1e-4
DEFAULT_TRIM_GAP = 10
DEFAULT_TRIM_MIN_SAMPLES = 100

# Scan trimming (dominant-axis selection to exclude slews where both axes move)
DEFAULT_TRIM_DOMINANT_AXIS = True         # if True, prefer segments where main-axis speed dominates cross-axis speed
DEFAULT_TRIM_AXIS_RATIO_MIN = 3.0         # require |v_main| / (|v_cross|+eps) >= this
DEFAULT_TRIM_VPERCENTILE = 95.0           # percentile of |v_main| used to define "moving" threshold (robust against spikes)
# Scan trimming (steady scan selection near Sun using ON-only & dominant-axis speed ratio)
DEFAULT_TRIM_STEADY_SCAN = True             # try steady-scan selection first (recommended)
DEFAULT_TRIM_USE_ON_ONLY = True             # ignore HOT/OFF even if id matches (position=="ON" only)
DEFAULT_TRIM_SCAN_SPEED_MIN_ARCSEC_S = 20.0 # minimum |v_main| [arcsec/s] to be considered a scan (rejects near-stationary drift)
DEFAULT_TRIM_XWIN_FACTOR = 1.2              # consider points within |offset_main| <= profile_xlim_deg * factor as "near Sun"
DEFAULT_TRIM_CROSS_OFFSET_MAX_DEG = 0.5     # allow cross-axis offset deviation from its median within this range [deg]
DEFAULT_TRIM_STEADY_CV_MAX = 0.8            # max robust CV (sigma/|mean|) of v_main within selected scan block     # allow cross-axis offset deviation from its median within this range [deg]

# Derivative behavior
DEFAULT_STRICT_DERIV = True

# Error handling
DEFAULT_CONTINUE_ON_ERROR = False

# Edge (limb) fitting on derivative profiles (two Gaussians)
DEFAULT_EDGE_FIT = True
DEFAULT_EDGE_FIT_WIN_DEG = 0.15         # fit window half-width around each derivative peak [deg]
DEFAULT_EDGE_FIT_THRESHOLD = 0.20       # keep points with |dY/dX| >= threshold * |peak| within the window
DEFAULT_HPBW_INIT_ARCSEC = 324.0        # initial guess for HPBW [arcsec] used to set sigma_init
DEFAULT_EDGE_FIT_PLOT_MAX_SCANS = 3     # how many scan pairs to show per summary page

# -----------------------------
# Robust table reader
# -----------------------------
def _safe_read_table(db, table_name: str) -> Optional[pd.DataFrame]:
    """
    Robust table reader for necstdb.

    Preferred in your environment:
        db.open_table(name).read(astype="df")
    """
    try:
        if hasattr(db, "open_table"):
            tab = db.open_table(table_name)
            if tab is not None and hasattr(tab, "read"):
                df = tab.read(astype="df")
                if isinstance(df, pd.DataFrame):
                    return df
    except Exception:
        pass

    try:
        df = necstdb.read_table(db, table_name)
        if isinstance(df, pd.DataFrame):
            return df
    except Exception:
        pass

    return None


def _norm_name(x) -> str:
    try:
        return str(x).strip()
    except Exception:
        return ""


def _pick_field_name(names, preferred, candidates):
    if preferred:
        p = _norm_name(preferred)
        for n in names:
            if _norm_name(n) == p:
                return n
        for n in names:
            if _norm_name(n).lower() == p.lower():
                return n

    for c in candidates:
        cc = _norm_name(c)
        for n in names:
            if _norm_name(n) == cc:
                return n
        for n in names:
            if _norm_name(n).lower() == cc.lower():
                return n

    low = [(_norm_name(n), _norm_name(n).lower()) for n in names]
    for c in candidates:
        cl = _norm_name(c).lower()
        for (n0, nl) in low:
            if cl and cl in nl:
                return n0
    return None


def _get_field(arr, name, dtype=float):
    return np.asarray(arr[name], dtype=dtype)


def _normalize_time_units(t_ref, t_other, label):
    tr = np.asarray(t_ref, dtype=float)
    to = np.asarray(t_other, dtype=float)

    mr = float(np.nanmedian(tr)) if tr.size else np.nan
    mo = float(np.nanmedian(to)) if to.size else np.nan
    if (not np.isfinite(mr)) or (not np.isfinite(mo)) or mr <= 0 or mo <= 0:
        return to, 1.0

    ratio = mo / mr

    if 5e2 < ratio < 2e3:
        return to / 1e3, 1e-3
    if 5e5 < ratio < 2e6:
        return to / 1e6, 1e-6
    if 5e8 < ratio < 2e9:
        return to / 1e9, 1e-9

    if 5e2 < (1.0 / ratio) < 2e3:
        return to * 1e3, 1e3
    if 5e5 < (1.0 / ratio) < 2e6:
        return to * 1e6, 1e6
    if 5e8 < (1.0 / ratio) < 2e9:
        return to * 1e9, 1e9

    return to, 1.0


def _read_table_raw_bytes(table):
    for k in ("raw", "buffer"):
        try:
            b = table.read(astype=k)
            if isinstance(b, (bytes, bytearray)):
                return bytes(b)
        except Exception:
            pass
    raise RuntimeError("cannot read raw/buffer from table")


def _get_table_dtype(table):
    for attr in ("dtype", "_dtype", "data_dtype", "_data_dtype"):
        try:
            dt = getattr(table, attr)
            if dt is not None:
                return dt
        except Exception:
            pass
    try:
        arr = table.read(num=1, astype="array")
        return arr.dtype
    except Exception:
        return None


def _read_structured_array_tolerant(db, table_name: str):
    tbl = db.open_table(table_name)
    try:
        return tbl.read(astype="array")
    except Exception as e:
        msg = str(e)
        raw = _read_table_raw_bytes(tbl)
        dt = _get_table_dtype(tbl)
        if dt is None:
            raise RuntimeError(f"cannot determine dtype for table {table_name} (original error: {msg})")
        item = int(dt.itemsize)
        if item <= 0:
            raise RuntimeError(f"invalid dtype itemsize for table {table_name}: {item}")
        nrec = len(raw) // item
        if nrec <= 0:
            raise RuntimeError(f"raw buffer too small for table {table_name}: len={len(raw)} itemsize={item}")
        raw2 = raw[: nrec * item]
        try:
            return np.frombuffer(raw2, dtype=dt)
        except Exception as e2:
            raise RuntimeError(f"tolerant frombuffer failed for table {table_name}: {e2} (orig: {msg})")


def _extract_spectral_from_structured(arr):
    names = list(arr.dtype.names or [])
    s_field = None
    for k in ("data", "spectrum", "spec"):
        if k in names:
            s_field = k
            break
    if s_field is None:
        raise RuntimeError(f"no spectrum-like field in spectral table. available={names}")

    t_spec, time_meta = _select_spectral_time_from_structured(arr)
    spec = np.asarray(arr[s_field])
    if spec.ndim == 1 and spec.dtype == object:
        spec2d = np.stack(spec, axis=0).astype(np.float32, copy=False)
    else:
        spec2d = spec.astype(np.float32, copy=False)

    # 1-channel spectra may be stored as a scalar numeric field per row, which
    # appears here as shape (N,) instead of (N, 1). Treat that as a valid single-
    # channel spectrum table rather than rejecting it.
    if spec2d.ndim == 1:
        spec2d = spec2d.reshape(-1, 1)

    if spec2d.ndim != 2:
        raise RuntimeError(f"spectrum array must be 2D. got {spec2d.shape}")

    return np.asarray(t_spec, dtype=float), spec2d, time_meta, s_field


def _to_degC_guess(temp_arr):
    x = np.asarray(temp_arr, dtype=float)
    if x.size == 0:
        return np.full_like(x, np.nan, dtype=float)
    med = float(np.nanmedian(x))
    if np.isfinite(med) and med > 120.0:
        return x - 273.15
    return x


def _to_hpa_guess(press_arr):
    x = np.asarray(press_arr, dtype=float)
    if x.size == 0:
        return np.full_like(x, np.nan, dtype=float)
    med = float(np.nanmedian(x))
    if np.isfinite(med) and med > 2000.0:
        return x / 100.0
    return x


def _to_rh01_guess(humid_arr):
    x = np.asarray(humid_arr, dtype=float)
    if x.size == 0:
        return np.full_like(x, np.nan, dtype=float)
    med = float(np.nanmedian(x))
    if np.isfinite(med) and med <= 1.5:
        rh = x
    else:
        rh = x / 100.0
    return np.clip(rh, 0.0, 1.0)

def _normalize_azel_correction_apply(value: Optional[str], default: str = "none") -> str:
    if value is None:
        return str(default)
    s = str(value).strip().lower()
    if not s:
        return str(default)
    if s == "minus":
        return "subtract"
    if s == "plus":
        return "add"
    if s in {"subtract", "add", "none"}:
        return s
    return str(default)


def _sanitize_time_col(preferred_time_col: str, *, label: str = "time") -> str:
    pref = str(preferred_time_col).strip()
    if pref.lower() == "recorded_time":
        print(f"[warn] {label}=recorded_time is not supported; forcing to 'time'")
        pref = "time"
    return pref


def _resolve_table_name(*, db_namespace: str, telescope: str, full_table: Optional[str], suffix: Optional[str], default_suffix: str) -> str:
    if full_table is not None:
        s = str(full_table).strip()
        if s:
            return s
    suf = str(suffix).strip() if suffix is not None and str(suffix).strip() else str(default_suffix)
    return f"{str(db_namespace)}-{str(telescope)}-{suf}"


def _apply_azel_correction(lon: np.ndarray, lat: np.ndarray, dlon: np.ndarray, dlat: np.ndarray, mode: str) -> Tuple[np.ndarray, np.ndarray]:
    m = _normalize_azel_correction_apply(mode, default="none")
    if m == "add":
        return (np.asarray(lon, dtype=float) + np.asarray(dlon, dtype=float)) % 360.0, np.asarray(lat, dtype=float) + np.asarray(dlat, dtype=float)
    if m == "subtract":
        return (np.asarray(lon, dtype=float) - np.asarray(dlon, dtype=float)) % 360.0, np.asarray(lat, dtype=float) - np.asarray(dlat, dtype=float)
    return np.asarray(lon, dtype=float) % 360.0, np.asarray(lat, dtype=float)


def _sanitize_meteo_series(arr, t_spec: np.ndarray, *, default: float, vmin: Optional[float] = None, vmax: Optional[float] = None) -> Tuple[np.ndarray, bool]:
    out = np.full_like(np.asarray(t_spec, dtype=float), float(default), dtype=float)
    used_default = True
    if arr is not None:
        try:
            cand = np.asarray(arr, dtype=float)
            if cand.shape == out.shape:
                mask = np.isfinite(cand)
                if vmin is not None:
                    mask &= cand >= float(vmin)
                if vmax is not None:
                    mask &= cand <= float(vmax)
                if np.any(mask):
                    out[mask] = cand[mask]
                    used_default = bool(np.any(~mask))
                else:
                    used_default = True
            else:
                used_default = True
        except Exception:
            used_default = True
    return out, bool(used_default)


def _finalize_outside_weather_meteo(weather_meteo, t_spec: np.ndarray, *, default_temperature_c: float, default_pressure_hpa: float, default_humidity_pct: float,
                                    temperature_min_c: float, temperature_max_c: float, pressure_min_hpa: float, pressure_max_hpa: float,
                                    humidity_min_pct: float, humidity_max_pct: float):
    temp_c, used_t = _sanitize_meteo_series(None if weather_meteo is None else weather_meteo.get("temp_c"), t_spec, default=float(default_temperature_c), vmin=float(temperature_min_c), vmax=float(temperature_max_c))
    press_hpa, used_p = _sanitize_meteo_series(None if weather_meteo is None else weather_meteo.get("press_hpa"), t_spec, default=float(default_pressure_hpa), vmin=float(pressure_min_hpa), vmax=float(pressure_max_hpa))
    humid_pct, used_h = _sanitize_meteo_series(None if weather_meteo is None else weather_meteo.get("humid_pct"), t_spec, default=float(default_humidity_pct), vmin=float(humidity_min_pct), vmax=float(humidity_max_pct))
    source = "default"
    if weather_meteo is not None:
        source = "mixed" if (used_t or used_p or used_h) else "outside_weather"
    return {
        "table": None if weather_meteo is None else weather_meteo.get("table"),
        "time_col": None if weather_meteo is None else weather_meteo.get("time_col"),
        "temp_c": np.asarray(temp_c, dtype=float),
        "press_hpa": np.asarray(press_hpa, dtype=float),
        "humid_pct": np.asarray(humid_pct, dtype=float),
        "humid_rh": np.asarray(np.clip(np.asarray(humid_pct, dtype=float) / 100.0, 0.0, 1.0), dtype=float),
        "used_default_temperature": bool(used_t),
        "used_default_pressure": bool(used_p),
        "used_default_humidity": bool(used_h),
        "source": source,
    }


def _estimate_tamb_k_with_source(*, weather_meteo=None, structured_arr=None, fallback_k: float = 300.0, min_k: float = 250.0, max_k: float = 330.0) -> Tuple[float, str]:
    if weather_meteo is not None:
        try:
            temp_c = weather_meteo.get("temp_c")
            if temp_c is not None:
                med_k = float(np.nanmedian(np.asarray(temp_c, dtype=float))) + 273.15
                if np.isfinite(med_k) and (float(min_k) <= med_k <= float(max_k)):
                    return med_k, "inside_weather"
        except Exception:
            pass
    if structured_arr is not None:
        names = list(structured_arr.dtype.names or [])
        for field in ("Tamb", "temp_k", "temperature"):
            name = _pick_field_name(names, field, [field])
            if name is None:
                continue
            try:
                arr = np.asarray(structured_arr[name], dtype=float)
                med = float(np.nanmedian(arr))
                if field == "temperature" and np.isfinite(med) and med < 120.0:
                    med = med + 273.15
                if np.isfinite(med) and (float(min_k) <= med <= float(max_k)):
                    return med, "spectral"
            except Exception:
                continue
    return float(fallback_k), "default"


def _load_weather_meteo_for_tspec_from_table(
    db,
    table_name: str,
    t_spec: np.ndarray,
    *,
    preferred_time_col: str = "time",
):
    try:
        arr = _read_structured_array_tolerant(db, table_name)
    except Exception as e:
        print(f"[warn] cannot read weather table {table_name}: {e}")
        return None

    names = list(arr.dtype.names or [])
    pref = _sanitize_time_col(preferred_time_col, label="weather_time_col")

    t_name = _pick_field_name(names, pref, ["time", "timestamp"])
    if (t_name is None) or (str(t_name).strip().lower() == "recorded_time"):
        print(f"[warn] weather table has no usable time column (time/timestamp). recorded_time is ignored. available={names}")
        return None

    t_w = _get_field(arr, t_name, float)
    t_w, s_w = _normalize_time_units(t_spec, t_w, "weather")
    if s_w != 1.0:
        print(f"[info] weather time scaled by {s_w} to match unix seconds")

    out = {"table": table_name, "time_col": str(t_name)}

    if "pressure" in names:
        try:
            out["press_hpa"] = _to_hpa_guess(_interp_edge_hold(t_w, _get_field(arr, "pressure", float), t_spec))
        except Exception:
            out["press_hpa"] = None
    else:
        out["press_hpa"] = None

    if "temperature" in names:
        try:
            out["temp_c"] = _to_degC_guess(_interp_edge_hold(t_w, _get_field(arr, "temperature", float), t_spec))
        except Exception:
            out["temp_c"] = None
    else:
        out["temp_c"] = None

    if "humidity" in names:
        try:
            humid_raw = _interp_edge_hold(t_w, _get_field(arr, "humidity", float), t_spec)
            out["humid_pct"] = np.asarray(_to_rh01_guess(humid_raw) * 100.0, dtype=float)
            out["humid_rh"] = np.asarray(_to_rh01_guess(humid_raw), dtype=float)
        except Exception:
            out["humid_pct"] = None
            out["humid_rh"] = None
    else:
        out["humid_pct"] = None
        out["humid_rh"] = None

    ok_any = False
    for k in ("press_hpa", "temp_c", "humid_pct"):
        v = out.get(k)
        if v is not None and np.any(np.isfinite(np.asarray(v, dtype=float))):
            ok_any = True
            break
    if not ok_any:
        print("[warn] weather table fields exist but all are NaN after interpolation; ignoring weather meteo")
        return None
    return out


def _load_weather_meteo_for_tspec(
    db,
    telescope: str,
    t_spec: np.ndarray,
    *,
    db_namespace: str = DEFAULT_DB_NAMESPACE,
    preferred_time_col: str = "time",
):
    table_name = f"{str(db_namespace)}-{str(telescope)}-weather-ambient"
    return _load_weather_meteo_for_tspec_from_table(db, table_name, t_spec, preferred_time_col=preferred_time_col)


# -----------------------------
# Helpers
# -----------------------------
def _az_wrap_diff(a_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    """Az difference a-b in deg wrapped to [-180, +180)."""
    a = np.asarray(a_deg, dtype=float)
    b = np.asarray(b_deg, dtype=float)
    return (a - b + 180.0) % 360.0 - 180.0


def _interp_lin(t_src: np.ndarray, y_src: np.ndarray, t_dst: np.ndarray) -> np.ndarray:
    t_src = np.asarray(t_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)
    m = np.isfinite(t_src) & np.isfinite(y_src)
    if np.count_nonzero(m) < 2:
        raise RuntimeError("Not enough finite points for linear interpolation.")
    return np.interp(t_dst, t_src[m], y_src[m])


def _interp_az_deg(t_src: np.ndarray, az_src_deg: np.ndarray, t_dst: np.ndarray) -> np.ndarray:
    """Interpolate azimuth in deg robustly by unwrapping."""
    az = np.asarray(az_src_deg, dtype=float) % 360.0
    az_rad = np.deg2rad(az)
    az_unw = np.unwrap(az_rad)
    az_i = _interp_lin(t_src, az_unw, t_dst)
    return (np.rad2deg(az_i) % 360.0)


def _dedup_mean(df: pd.DataFrame, tcol: str, cols: List[str]) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Deduplicate by time (tcol) and mean over duplicated timestamps.
    Returns (t_unique, df_unique) where df_unique has the columns in 'cols' plus 't'.
    """
    d = df[[tcol] + cols].copy()
    d[tcol] = pd.to_numeric(d[tcol], errors="coerce")
    d = d[np.isfinite(d[tcol])]
    if len(d) == 0:
        raise RuntimeError("No finite timestamps after cleaning.")
    g = d.groupby(tcol, sort=True)
    out = g[cols].mean(numeric_only=True)
    out.insert(0, "t", out.index.values.astype(float))
    t_unique = out["t"].to_numpy(float)
    return t_unique, out.reset_index(drop=True)


def _decode_label(v) -> str:
    try:
        if isinstance(v, (bytes, bytearray)):
            return v.decode(errors="ignore").strip()
        return str(v).strip()
    except Exception:
        return ""


def _decode_position(v) -> str:
    """Decode a position/mode label strictly: bytes->str, upper+strip."""
    try:
        if isinstance(v, (bytes, bytearray)):
            s = v.decode(errors="ignore")
        else:
            s = str(v)
        return s.upper().strip()
    except Exception:
        return ""


def _estimate_tamb_k(
    *,
    weather_meteo=None,
    structured_arr=None,
    fallback_k: float = 300.0,
) -> float:
    """Estimate ambient load temperature Tamb [K] from necstdb-side information."""
    if weather_meteo is not None:
        try:
            temp_c = weather_meteo.get("temp_c")
            if temp_c is not None:
                med_c = float(np.nanmedian(np.asarray(temp_c, dtype=float)))
                med_k = med_c + 273.15
                if np.isfinite(med_k) and (50.0 <= med_k <= 400.0):
                    return med_k
        except Exception:
            pass

    if structured_arr is not None:
        names = list(structured_arr.dtype.names or [])
        for field in ("Tamb", "temp_k", "temperature"):
            name = _pick_field_name(names, field, [field])
            if name is None:
                continue
            try:
                arr = np.asarray(structured_arr[name], dtype=float)
                med = float(np.nanmedian(arr))
                if field == "temperature" and np.isfinite(med) and med < 120.0:
                    med = med + 273.15
                if np.isfinite(med) and (50.0 <= med <= 400.0):
                    return med
            except Exception:
                continue

    return float(fallback_k)


def _interp_edge_hold(t_known: np.ndarray, y_known: np.ndarray, t_dst: np.ndarray) -> np.ndarray:
    """np.interp with edge-hold, robust to duplicates and NaNs."""
    t_known = np.asarray(t_known, dtype=float)
    y_known = np.asarray(y_known, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)
    m = np.isfinite(t_known) & np.isfinite(y_known)
    if np.count_nonzero(m) < 2:
        raise RuntimeError("Need at least 2 finite points for interpolation.")
    d = pd.DataFrame({"t": t_known[m], "y": y_known[m]})
    d = d.groupby("t", sort=True)["y"].mean().reset_index()
    tt = d["t"].to_numpy(float)
    yy = d["y"].to_numpy(float)
    if tt.size < 2:
        raise RuntimeError("Need at least 2 unique time points for interpolation.")
    return np.interp(t_dst, tt, yy, left=float(yy[0]), right=float(yy[-1]))



# -----------------------------
# Encoder smoothing (safe local position smoothing; no integration)
# -----------------------------
def _smooth_encoder_interp_local_by_id(
    t_spec: np.ndarray,
    id_vals: np.ndarray,
    az_enc_raw_deg: np.ndarray,
    el_enc_raw_deg: np.ndarray,
    *,
    win_sec: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Safely smooth interpolated encoder positions without cumulative drift.

    Design goals
    ------------
    - Do *not* differentiate encoder positions to rates
    - Do *not* integrate rates back to positions
    - Keep smoothing local in time so errors do not accumulate downstream
    - Isolate contiguous scan/id blocks so one scan cannot contaminate another

    Method
    ------
    1) Interpolate raw encoder lon/lat to t_spec (done upstream)
    2) For each contiguous id block, locally smooth az/el themselves using a
       centered continuous-time boxcar mean
    3) Azimuth is unwrapped only within each local block, then wrapped back

    Notes
    -----
    - Valid non-empty id labels are required. If id is unavailable or mismatched,
      the caller should disable smoothing and keep the raw interpolated encoder
      positions.
    - This helper never integrates over time, so it cannot create cumulative drift.
    """
    t = np.asarray(t_spec, dtype=float)
    az = np.asarray(az_enc_raw_deg, dtype=float).copy()
    el = np.asarray(el_enc_raw_deg, dtype=float).copy()

    w = builtins.max(0.0, float(win_sec))
    if w <= 0.0 or t.size < 4:
        return az, el

    labels = np.asarray(id_vals, dtype=object)
    if labels.shape[0] != t.shape[0]:
        return az, el
    labels = np.array([_decode_label(v) for v in labels], dtype=object)
    if not np.any(labels != ""):
        return az, el

    n = t.size
    i = 0
    while i < n:
        j = i + 1
        while j < n and labels[j] == labels[i]:
            j += 1

        tt = t[i:j]
        azb = az[i:j]
        elb = el[i:j]

        good = np.isfinite(tt) & np.isfinite(azb) & np.isfinite(elb)
        if np.count_nonzero(good) >= 4:
            az_unw = np.rad2deg(np.unwrap(np.deg2rad(azb[good] % 360.0)))
            az_sm = continuous_boxcar_mean_fast(tt[good], az_unw, T=w)
            el_sm = continuous_boxcar_mean_fast(tt[good], elb[good], T=w)

            az_out = azb.copy()
            el_out = elb.copy()
            az_out[good] = np.mod(az_sm, 360.0)
            el_out[good] = el_sm

            az[i:j] = az_out
            el[i:j] = el_out

        i = j

    return az, el



# -----------------------------
# Ripple removal (non-uniform -> uniform -> FFT peak -> notch -> restore)
# -----------------------------
def _pandas_freq_from_dt(dt_sec: float) -> str:
    """Convert dt_sec to a pandas resample frequency string with millisecond resolution."""
    dt = float(dt_sec)
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt_sec must be positive and finite.")
    ms = int(builtins.max(1, round(dt * 1000.0)))
    return f"{ms}ms"


def _resample_uniform_time(
    t: np.ndarray,
    y: np.ndarray,
    *,
    dt_sec: Optional[float] = None,
    dt_min: float = 0.001,
    dt_max: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Resample non-uniform (t, y) to a uniform grid using pandas time interpolation.

    Practical policy (notch-friendly but safe):
      - If dt_sec is not specified, use ~median(diff(t)) (i.e., do NOT force oversampling).
        Oversampling by large factors can create filter artifacts because the additional
        samples contain no new information (they are interpolated).
      - A small amount of oversampling can still be requested explicitly via dt_sec.
      - Filling is time-interpolation within the resampled span with edge fill.

    Parameters
    ----------
    t : seconds (float), non-uniform
    y : values (float)
    dt_sec : uniform step in seconds. If None/<=0, infer from median(diff(t)).
    dt_min, dt_max : absolute bounds for inferred dt_sec (sec).

    Returns
    -------
    t_u : uniform seconds array
    y_u : resampled values (no NaN within span)
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if t.ndim != 1 or y.ndim != 1 or t.size != y.size:
        raise ValueError("t and y must be 1D and same length.")
    m = np.isfinite(t) & np.isfinite(y)
    if np.count_nonzero(m) < 8:
        raise RuntimeError("Not enough finite samples for resampling.")
    t = t[m]
    y = y[m]
    o = np.argsort(t)
    t = t[o]
    y = y[o]

    if dt_sec is None or (not np.isfinite(dt_sec)) or dt_sec <= 0:
        dts = np.diff(t)
        dts = dts[np.isfinite(dts) & (dts > 0)]
        if dts.size < 4:
            raise RuntimeError("Cannot infer dt_sec from timestamps.")
        dt_sec = float(np.median(dts))

        # Keep dt within a reasonable absolute range.
        dt_sec = builtins.max(float(dt_min), builtins.min(float(dt_max), dt_sec))
    else:
        dt_sec = float(dt_sec)

    freq = _pandas_freq_from_dt(dt_sec)
    idx = pd.to_datetime(t, unit="s")
    ser = pd.Series(y, index=idx).sort_index()
    ser = ser[~ser.index.duplicated(keep="first")]

    # Resample to uniform grid then interpolate to fill gaps.
    ser_u = ser.resample(freq).mean()
    ser_u = ser_u.interpolate(method="time", limit_direction="both").bfill().ffill()

    y_u = ser_u.to_numpy(float)
    if not np.all(np.isfinite(y_u)):
        # last resort: linear fill
        ser_u = ser_u.interpolate(method="linear", limit_direction="both").bfill().ffill()
        y_u = ser_u.to_numpy(float)

    t_u = ser_u.index.view("int64").astype(np.float64) / 1e9
    return t_u, y_u


def _find_peak_frequency_fft(
    y_u: np.ndarray,
    *,
    fs_hz: float,
    target_hz: float,
    search_hz: float,
) -> float:
    """Find a dominant peak frequency near target_hz within ±search_hz using FFT.

    Notes
    -----
    - Uses a Hann window to reduce leakage.
    - Refines the peak estimate with a 3-point quadratic interpolation (in log-magnitude)
      to reduce bin-quantization error.
    """
    y = np.asarray(y_u, dtype=float)
    if y.ndim != 1 or y.size < 32:
        return float(target_hz)

    y0 = y - float(np.mean(y))
    w = np.hanning(y0.size)
    yf = np.fft.rfft(y0 * w)
    f = np.fft.rfftfreq(y0.size, d=1.0 / float(fs_hz))

    f0 = float(target_hz)
    sw = float(search_hz)
    lo = builtins.max(0.0, f0 - sw)
    hi = f0 + sw

    band = (f >= lo) & (f <= hi) & (f > 0)
    if not np.any(band):
        return f0

    band_idx = np.where(band)[0]
    mag = np.abs(yf[band_idx])
    if mag.size == 0:
        return f0

    k_rel = int(np.argmax(mag))
    k0 = int(band_idx[k_rel])
    f_peak = float(f[k0])
    if not np.isfinite(f_peak) or f_peak <= 0:
        return f0

    # Quadratic interpolation around the peak bin (parabolic fit in log-magnitude).
    if 1 <= k0 < (f.size - 1):
        y1 = float(np.abs(yf[k0 - 1]))
        y2 = float(np.abs(yf[k0]))
        y3 = float(np.abs(yf[k0 + 1]))
        if y1 > 0 and y2 > 0 and y3 > 0:
            a = float(np.log(y1))
            b = float(np.log(y2))
            c = float(np.log(y3))
            den = (a - 2.0 * b + c)
            if np.isfinite(den) and abs(den) > 1e-12:
                delta = 0.5 * (a - c) / den
                delta = float(np.clip(delta, -1.0, 1.0))
                df = float(f[1] - f[0])
                f_peak = float(f_peak + delta * df)

    if not np.isfinite(f_peak) or f_peak <= 0:
        return f0
    return f_peak


def _notch_filter_harmonics(
    y_u: np.ndarray,
    *,
    fs_hz: float,
    f_peak_hz: float,
    bw_hz: float,
    max_harm: int,
    order: int,
    passes: int = 1,
) -> np.ndarray:
    """Apply Butterworth bandstop (notch) around f_peak and its harmonics using filtfilt.

    Nyquist-aware policy
    --------------------
    For each harmonic k, we notch [k*f_peak - bw/2, k*f_peak + bw/2].
    If the upper edge reaches Nyquist, the notch is invalid and is skipped.
    The effective harmonic limit is clipped by:
        k_max = floor((0.99*nyq - bw/2) / f_peak)

    Robustness notes
    ----------------
    - For each harmonic and each pass, try filtfilt with a safe padlen; if it fails, retry with padlen=0.
    - If both fail, skip that harmonic and continue.
    """
    y = np.asarray(y_u, dtype=float)
    if y.ndim != 1:
        raise ValueError("y_u must be 1D.")
    if y.size < 32:
        return y.copy()

    fs = float(fs_hz)
    nyq = 0.5 * fs
    bw = float(bw_hz)
    f1 = float(f_peak_hz)

    if (not np.isfinite(bw)) or bw <= 0:
        return y.copy()
    if (not np.isfinite(f1)) or f1 <= 0:
        return y.copy()
    if (not np.isfinite(nyq)) or nyq <= 0:
        return y.copy()

    # Clip harmonics so that the notch band does not exceed Nyquist.
    k_max = int(math.floor((0.99 * nyq - 0.5 * bw) / f1))
    if k_max < 1:
        return y.copy()

    n_h_eff = int(builtins.max(1, builtins.min(int(max_harm), k_max)))
    ord2 = int(builtins.max(1, order))
    n_pass = int(builtins.max(1, passes))

    y_f = y.copy()
    for k in range(1, n_h_eff + 1):
        f0 = f1 * float(k)
        if not np.isfinite(f0) or f0 <= 0:
            continue
        if f0 >= 0.99 * nyq:
            break

        low = f0 - 0.5 * bw
        high = f0 + 0.5 * bw

        # Validity checks
        if (not np.isfinite(low)) or (not np.isfinite(high)):
            continue
        if low <= 1e-6:
            continue
        if high >= 0.99 * nyq:
            continue
        if low >= high:
            continue

        wn = [low / nyq, high / nyq]
        try:
            b, a = signal.butter(ord2, wn, btype="bandstop")
        except Exception:
            continue

        padlen = 3 * (builtins.max(len(a), len(b)) - 1)
        padlen = int(builtins.min(padlen, builtins.max(0, y_f.size - 1)))

        for _ in range(n_pass):
            try:
                y_f = signal.filtfilt(b, a, y_f, padlen=padlen)
            except Exception:
                try:
                    y_f = signal.filtfilt(b, a, y_f, padlen=0)
                except Exception:
                    break

    return y_f


def _bandpower_around(
    y_u: np.ndarray,
    *,
    fs_hz: float,
    f_peak_hz: float,
    band_hz: float,
    max_harm: int,
) -> float:
    """Compute FFT band power around f_peak and its harmonics (sum of |FFT|^2 in bands).

    Nyquist-aware clipping:
      For each harmonic k, sum power within [k*f_peak - band/2, k*f_peak + band/2],
      but only for harmonics that keep the band within Nyquist.
    """
    y = np.asarray(y_u, dtype=float)
    if y.ndim != 1 or y.size < 32:
        return float("inf")

    fs = float(fs_hz)
    nyq = 0.5 * fs
    bw = float(band_hz)
    f1 = float(f_peak_hz)

    if (not np.isfinite(nyq)) or nyq <= 0:
        return float("inf")
    if (not np.isfinite(f1)) or f1 <= 0:
        return float("inf")
    if (not np.isfinite(bw)) or bw <= 0:
        bw = 0.5

    # Clip harmonics so that the band does not exceed Nyquist.
    k_max = int(math.floor((0.99 * nyq - 0.5 * bw) / f1))
    if k_max < 1:
        return float("inf")
    k_max = int(builtins.max(1, builtins.min(int(max_harm), k_max)))

    y0 = y - float(np.mean(y))
    w = np.hanning(y0.size)
    yf = np.fft.rfft(y0 * w)
    f = np.fft.rfftfreq(y0.size, d=1.0 / fs)
    p = np.abs(yf) ** 2

    total = 0.0
    for k in range(1, k_max + 1):
        f0 = f1 * float(k)
        if f0 <= 0 or f0 >= 0.99 * nyq:
            break
        lo = builtins.max(0.0, f0 - 0.5 * bw)
        hi = builtins.min(0.99 * nyq, f0 + 0.5 * bw)
        m = (f >= lo) & (f <= hi)
        if np.any(m):
            total += float(np.sum(p[m]))
    return float(total)




def continuous_boxcar_mean_fast(t: np.ndarray, y: np.ndarray, T: float) -> np.ndarray:
    """Fast continuous-time boxcar mean for non-uniform samples.

    Implementation (same idea as the legacy script that worked well for NECST sun scans):
      - Resample to a fine uniform grid (auto dt ~ median(dt)/2, clamped to 1..10 ms)
      - Apply centered rolling mean with window length ~T
      - Interpolate the smoothed series back to the original timestamps

    Parameters
    ----------
    t : np.ndarray
        Timestamps [s] (non-uniform)
    y : np.ndarray
        Values
    T : float
        Boxcar window [s]

    Returns
    -------
    y_smooth : np.ndarray
        Smoothed values aligned to the original t.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if t.ndim != 1 or y.ndim != 1 or t.size != y.size:
        raise ValueError("t and y must be 1D and same length.")
    if t.size < 2:
        return y.astype(float)

    # sort by time
    o = np.argsort(t)
    t = t[o]
    y = y[o]

    # choose resample dt (1..10 ms)
    dts = np.diff(t)
    dts = dts[np.isfinite(dts) & (dts > 0)]
    dt_med = float(np.nanmedian(dts)) if dts.size else 0.01
    target_dt = 0.5 * dt_med
    target_dt = builtins.max(0.001, builtins.min(0.01, float(target_dt)))
    resample_ms = int(round(target_dt * 1000.0))
    resample_ms = builtins.max(1, resample_ms)
    rule = f"{resample_ms}ms"
    dt_u = resample_ms / 1000.0

    idx = pd.to_datetime(t, unit="s")
    s = pd.Series(y.astype(float), index=idx)
    s = s[~s.index.duplicated(keep="first")]
    s_u = s.resample(rule).mean().interpolate(method="time").bfill().ffill()

    w = int(round(float(T) / dt_u)) if (np.isfinite(T) and float(T) > 0) else 1
    w = builtins.max(1, w)
    if w % 2 == 0:
        w += 1

    s_sm = s_u.rolling(window=w, center=True, min_periods=1).mean()

    df_u = pd.DataFrame({"smooth": s_sm})
    df_o = pd.DataFrame({"dummy": np.zeros(len(idx), dtype=float)}, index=idx)
    full_idx = df_u.index.union(df_o.index).sort_values()
    df_m = df_u.reindex(full_idx).interpolate(method="time").bfill().ffill()

    out = df_m.loc[df_o.index, "smooth"].to_numpy(float)

    # restore original (unsorted) order
    inv = np.empty_like(o)
    inv[o] = np.arange(o.size)
    return out[inv]


def remove_ripple_notch(
    t: np.ndarray,
    y: np.ndarray,
    *,
    target_hz: float = 1.2,
    search_hz: float = 0.3,
    bw_hz: float = 0.5,
    max_harm: int = 8,
    trend_win_sec: float = 5.0,
    notch_passes: int = 2,
    resample_dt_sec: Optional[float] = None,
    order: int = 2,
    eval_band_hz: Optional[float] = None,
    model: str = "auto",
    # --- safety/robustness knobs (kept internal by default)
    baseline_frac: float = 0.05,
    baseline_abs: float = 0.0,
    onsource_frac: float = 0.50,
    min_improve_frac: float = 0.05,
    penalty_lambda: float = 2.0,
    min_cycles_on_source: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Remove periodic ripple near target_hz (and harmonics) from non-uniform time series via notch filtering.

    Important practical note (fast scans):
      If the scan is very fast, the astrophysical/profile signal can have substantial
      power around ~1 Hz (edge/plateau changes occur on second timescales). In that case,
      a notch filter can remove real signal and worsen the profile/derivative. To avoid
      this, the 'auto' mode may decide to *not apply* the correction (returns raw).

    Returns
    -------
    y_corr : corrected y (same shape as input)
    trend  : extracted trend (baseline)
    ripple : estimated ripple component in y units (y - y_corr)
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if t.ndim != 1 or y.ndim != 1 or t.size != y.size:
        raise ValueError("t and y must be 1D and same length.")
    if t.size < 32:
        return y.copy(), y.copy(), np.zeros_like(y, dtype=float)

    # Sort and keep mapping to original order
    o = np.argsort(t)
    t_s = t[o]
    y_s = y[o]

    # Build a datetime index for time-interpolation back to original timestamps
    idx = pd.to_datetime(t_s, unit="s")

    # Baseline / on-source masks (for robust evaluation in auto mode)
    y_abs = np.abs(y_s)
    y95 = float(np.nanpercentile(y_abs[np.isfinite(y_abs)], 95)) if np.any(np.isfinite(y_abs)) else 0.0
    base_thr = builtins.max(float(baseline_abs), float(baseline_frac) * y95)
    on_thr = float(onsource_frac) * y95

    base_mask = np.isfinite(y_s) & (y_abs <= base_thr)
    on_mask = np.isfinite(y_s) & (y_abs >= on_thr)

    # If baseline region is too small (e.g., scan never leaves the Sun), fall back to full.
    if np.count_nonzero(base_mask) < 32:
        base_mask = np.isfinite(y_s)

    # 1) Rough trend and peak frequency estimation (use multiplicative residual where stable)
    f0 = float(target_hz)
    T_rough = builtins.max(float(trend_win_sec), 2.0 / builtins.max(f0, 1e-6))
    s_rough = continuous_boxcar_mean_fast(t_s, y_s, T=T_rough)
    s_rough = np.asarray(s_rough, dtype=float)

    eps = float(np.nanmedian(np.abs(s_rough))) * 1e-6 + 1e-12
    s_rough_safe = np.where(np.isfinite(s_rough), s_rough, eps)
    s_rough_safe = np.where(np.abs(s_rough_safe) >= eps, s_rough_safe, eps)

    # Stabilize multiplicative residual where |trend| is sufficiently large.
    s_abs0 = np.abs(s_rough_safe)
    s_max0 = float(np.nanmax(s_abs0)) if np.any(np.isfinite(s_abs0)) else 0.0
    s_floor0 = builtins.max(eps, 0.05 * s_max0)
    m0 = np.isfinite(s_abs0) & (s_abs0 >= s_floor0) & np.isfinite(y_s)

    r_rough = np.zeros_like(y_s, dtype=float)
    if np.count_nonzero(m0) >= 20:
        r_rough[m0] = y_s[m0] / s_rough_safe[m0] - 1.0
    else:
        r_rough = y_s / s_rough_safe - 1.0

    # Peak search: prefer baseline-only to avoid contaminating the peak pick by scan-shape power.
    r_for_peak = r_rough.copy()
    r_for_peak[~base_mask] = 0.0

    try:
        t_u0, r_u0 = _resample_uniform_time(t_s, r_for_peak, dt_sec=resample_dt_sec)
        dt0 = float(np.median(np.diff(t_u0)))
        fs0 = 1.0 / dt0
        f_peak = _find_peak_frequency_fft(r_u0, fs_hz=fs0, target_hz=f0, search_hz=float(search_hz))
    except Exception:
        f_peak = f0

    if not np.isfinite(f_peak) or f_peak <= 0:
        f_peak = f0

    # 2) Cycle-locked trend (integer number of ripple cycles)
    ripple_period = 1.0 / builtins.max(float(f_peak), 1e-6)
    cycles = int(builtins.max(2, int(round(float(trend_win_sec) / ripple_period))))
    Ttrend = float(cycles) * ripple_period

    s_trend_s = continuous_boxcar_mean_fast(t_s, y_s, T=Ttrend).astype(float)

    # safe trend for multiplicative model
    eps2 = float(np.nanmedian(np.abs(s_trend_s))) * 1e-6 + 1e-12
    s_safe = s_trend_s.copy()
    s_safe[~np.isfinite(s_safe)] = eps2
    s_safe[np.abs(s_safe) < eps2] = eps2

    # Decide where the multiplicative model is well-conditioned.
    s_abs = np.abs(s_trend_s)
    s_max = float(np.nanmax(s_abs)) if np.any(np.isfinite(s_abs)) else 0.0
    s_floor = builtins.max(eps2, 0.05 * s_max)
    mul_mask = np.isfinite(s_abs) & (s_abs >= s_floor) & np.isfinite(y_s)

    # residuals
    r_add = y_s - s_trend_s
    r_mul = np.zeros_like(y_s, dtype=float)
    r_mul[mul_mask] = y_s[mul_mask] / s_safe[mul_mask] - 1.0

    # helper: filter residual and bring back to original time grid
    def _filter_residual(r: np.ndarray) -> Tuple[np.ndarray, float, float, float]:
        # Resample r and baseline mask together to keep alignment.
        # Use the same dt policy as _resample_uniform_time (median diff(t) if dt not specified),
        # but do it in one pandas resampling so r and mask share the exact grid.
        idx_src = pd.to_datetime(t_s, unit="s")
        df_src = pd.DataFrame({"r": r, "m": base_mask.astype(float)}, index=idx_src).sort_index()

        if resample_dt_sec is None or (not np.isfinite(resample_dt_sec)) or float(resample_dt_sec) <= 0:
            dts = np.diff(t_s)
            dts = dts[np.isfinite(dts) & (dts > 0)]
            if dts.size < 4:
                raise RuntimeError("Cannot infer resample dt from timestamps.")
            dt_use = float(np.median(dts))
            dt_use = builtins.max(0.001, builtins.min(1.0, dt_use))
        else:
            dt_use = float(resample_dt_sec)

        freq = _pandas_freq_from_dt(dt_use)

        df_u = df_src.resample(freq).mean()
        df_u = df_u.interpolate(method="time", limit_direction="both").bfill().ffill()

        r_u = df_u["r"].to_numpy(float)
        m_u = df_u["m"].to_numpy(float)
        t_u = df_u.index.view("int64").astype(np.float64) / 1e9

        if r_u.size < 32:
            # Too short for a stable notch; return no-op.
            ser_u = pd.Series(r_u, index=df_u.index)
            ser_back = ser_u.reindex(ser_u.index.union(idx)).interpolate(method="time", limit_direction="both").reindex(idx)
            r_back = ser_back.to_numpy(float)
            p0 = 0.0
            p1 = 0.0
            fs = 1.0 / builtins.max(1e-6, dt_use)
            return r_back, p0, p1, float(fs)

        dt = float(np.median(np.diff(t_u)))
        fs = 1.0 / dt

        # apply notch (and harmonics)
        r_u_f = _notch_filter_harmonics(
            r_u,
            fs_hz=fs,
            f_peak_hz=float(f_peak),
            bw_hz=float(bw_hz),
            max_harm=int(max_harm),
            order=int(order),
            passes=int(notch_passes),
        )

        # back to original timestamps
        ser_u = pd.Series(r_u_f, index=df_u.index)
        ser_back = ser_u.reindex(ser_u.index.union(idx)).interpolate(method="time", limit_direction="both").reindex(idx)
        r_back = ser_back.to_numpy(float)

        # band power (pre/post) evaluated mainly on baseline (mask-weighted)
        bw_eval = float(eval_band_hz) if (eval_band_hz is not None and np.isfinite(eval_band_hz) and eval_band_hz > 0) else float(bw_hz)
        r_u_w = r_u * m_u
        r_u_f_w = r_u_f * m_u
        p_before = _bandpower_around(r_u_w, fs_hz=fs, f_peak_hz=float(f_peak), band_hz=bw_eval, max_harm=int(max_harm))
        p_after = _bandpower_around(r_u_f_w, fs_hz=fs, f_peak_hz=float(f_peak), band_hz=bw_eval, max_harm=int(max_harm))

        return r_back, float(p_before), float(p_after), float(fs)

    # 3) notch filter in both models
    r_add_nf, p_add0, p_add1, fs_add = _filter_residual(r_add)
    r_mul_nf, p_mul0, p_mul1, fs_mul = _filter_residual(r_mul)

    # ripple estimates
    a_hat = np.nan_to_num(r_add - r_add_nf)  # additive ripple (y units)
    m_hat = np.nan_to_num(r_mul - r_mul_nf)  # multiplicative ripple (fraction)
    m_hat = np.where(mul_mask, m_hat, 0.0)

    y_add = y_s - a_hat

    # Robust multiplicative reconstruction
    denom = 1.0 + m_hat
    denom_floor = 0.2
    valid_mul = mul_mask & np.isfinite(denom) & (denom > denom_floor)

    y_mul = np.full_like(y_s, np.nan, dtype=float)
    y_mul[valid_mul] = y_s[valid_mul] / denom[valid_mul]
    y_mul = np.where(np.isfinite(y_mul), y_mul, y_add)

    # helper: normalized on-source distortion
    def _norm_diff(y_new: np.ndarray) -> float:
        if np.count_nonzero(on_mask) < 10:
            return 0.0
        ref = y_s[on_mask]
        den = float(np.sqrt(np.nanmean(ref * ref))) + 1e-12
        dif = y_new[on_mask] - ref
        num = float(np.sqrt(np.nanmean(dif * dif)))
        return num / den

    # compute improvement ratios (avoid divide by zero)
    def _ratio(p0: float, p1: float) -> float:
        if not (np.isfinite(p0) and p0 > 0 and np.isfinite(p1)):
            return 1.0
        return float(p1 / p0)

    ratio_add = _ratio(p_add0, p_add1)
    ratio_mul = _ratio(p_mul0, p_mul1)

    # Fast scan guard: if the on-source span contains too few ripple cycles,
    # penalize changing the on-source shape (notch may remove real structure).
    penalty = float(penalty_lambda)
    if np.count_nonzero(on_mask) >= 2:
        t_on = float(t_s[on_mask][-1] - t_s[on_mask][0])
        cycles_on = t_on * float(f_peak)
        if np.isfinite(cycles_on) and cycles_on < float(min_cycles_on_source):
            penalty = builtins.max(penalty, 5.0)

    diff_add = _norm_diff(y_add)
    diff_mul = _norm_diff(y_mul)

    # auto selection: allow "none" if correction is not beneficial or is too distortive
    mdl = str(model).lower().strip()
    if mdl not in ["auto", "add", "mul"]:
        mdl = "auto"

    if mdl == "add":
        y_corr_s = y_add
    elif mdl == "mul":
        y_corr_s = y_mul
    else:
        cost_none = 1.0  # baseline: no correction
        cost_add = ratio_add + penalty * diff_add
        cost_mul = ratio_mul + penalty * diff_mul

        # If multiplicative is mostly invalid, de-prioritize it.
        n_mul = int(np.count_nonzero(mul_mask))
        n_ok = int(np.count_nonzero(valid_mul))
        if n_mul >= 50:
            ok_frac = n_ok / float(n_mul)
            if ok_frac < 0.8:
                cost_mul = float("inf")

        best = min([(cost_none, "none"), (cost_add, "add"), (cost_mul, "mul")], key=lambda x: x[0])[1]

        # Require a minimum improvement in band power for add/mul; otherwise keep raw.
        if best in ["add", "mul"]:
            best_ratio = ratio_add if best == "add" else ratio_mul
            if not (np.isfinite(best_ratio) and best_ratio <= (1.0 - float(min_improve_frac))):
                best = "none"

        if best == "add":
            y_corr_s = y_add
        elif best == "mul":
            y_corr_s = y_mul
        else:
            y_corr_s = y_s.copy()

    ripple_y_s = y_s - y_corr_s

    # restore original order
    inv = np.empty_like(o)
    inv[o] = np.arange(o.size)
    y_corr = y_corr_s[inv]
    s_trend = s_trend_s[inv]
    ripple_y = ripple_y_s[inv]

    return y_corr.astype(float), s_trend.astype(float), ripple_y.astype(float)


# -----------------------------
# Ripple preset resolution (preset + manual overrides, per scan)
# -----------------------------
def make_ripple_policy_from_args(args) -> dict:
    """Build a ripple policy dict from argparse args.

    - Preset gives a baseline of parameters.
    - Manual --ripple-* options (when provided) override the preset.
    """
    return dict(
        preset=str(getattr(args, "ripple_preset", DEFAULT_RIPPLE_PRESET)).lower().strip(),
        model=str(getattr(args, "ripple_model", DEFAULT_RIPPLE_MODEL)).lower().strip(),
        target_hz=float(getattr(args, "ripple_target_hz", DEFAULT_RIPPLE_TARGET_HZ)),
        search_hz=float(getattr(args, "ripple_search_hz", DEFAULT_RIPPLE_SEARCH_HZ)),
        # overrides (None means use preset)
        bw_hz=getattr(args, "ripple_bw_hz", None),
        max_harm=getattr(args, "ripple_max_harm", None),
        order=getattr(args, "ripple_order", None),
        notch_passes=getattr(args, "ripple_notch_pass", None),
        trend_win_sec=getattr(args, "ripple_trend_win_sec", None),
        resample_dt_sec=getattr(args, "ripple_resample_dt_sec", None),
        eval_band_hz=getattr(args, "ripple_eval_band_hz", None),
    )


def _estimate_scan_speed_arcsec_s(t: np.ndarray, x_deg: np.ndarray) -> float:
    """Estimate scan speed from (t, x) as a robust median of |dx/dt| [arcsec/s]."""
    t = np.asarray(t, dtype=float)
    x = np.asarray(x_deg, dtype=float)
    if t.ndim != 1 or x.ndim != 1 or t.size != x.size or t.size < 5:
        return float("nan")
    dt = np.diff(t)
    dx = np.diff(x)
    m = np.isfinite(dt) & (dt > 0) & np.isfinite(dx)
    if np.count_nonzero(m) < 3:
        return float("nan")
    v = np.abs(dx[m] / dt[m]) * 3600.0
    # drop extreme outliers
    lo = float(np.nanpercentile(v, 5))
    hi = float(np.nanpercentile(v, 95))
    vv = v[(v >= lo) & (v <= hi)]
    if vv.size >= 3:
        return float(np.nanmedian(vv))
    return float(np.nanmedian(v))


def resolve_ripple_cfg_for_scan(
    ripple_policy: Optional[dict],
    *,
    t: np.ndarray,
    x_main_deg: np.ndarray,
    profile_xlim_deg: float,
) -> dict:
    """Resolve the effective ripple-notch configuration for this scan segment.

    Preset selection
    ----------------
    preset == auto:
      choose a baseline based on estimated scan speed (arcsec/s)
    preset in {safe, normal, strong}:
      use RIPPLE_PRESETS[preset] as baseline.

    Manual overrides
    ----------------
    If an override key is not None, it replaces the preset value.

    Notes
    -----
    Nyquist-aware harmonic clipping is applied later inside _notch_filter_harmonics/_bandpower_around.
    """
    pol = ripple_policy or {}
    preset = str(pol.get("preset", DEFAULT_RIPPLE_PRESET)).lower().strip()
    model = str(pol.get("model", DEFAULT_RIPPLE_MODEL)).lower().strip()
    if model not in ["auto", "add", "mul"]:
        model = "auto"

    # Estimate scan speed near the Sun (|x| <= profile_xlim_deg). Fallback to full span if too few points.
    t = np.asarray(t, dtype=float)
    x = np.asarray(x_main_deg, dtype=float)
    xl = float(profile_xlim_deg) if (np.isfinite(profile_xlim_deg) and profile_xlim_deg > 0) else float("inf")
    near = np.isfinite(t) & np.isfinite(x) & (np.abs(x) <= xl)
    if np.count_nonzero(near) >= 10:
        v_scan = _estimate_scan_speed_arcsec_s(t[near], x[near])
    else:
        v_scan = _estimate_scan_speed_arcsec_s(t, x)

    # Baseline from preset
    if preset == "auto":
        if np.isfinite(v_scan):
            if v_scan < 150.0:
                base = dict(bw_hz=0.8, max_harm=10, notch_passes=5, trend_win_sec=6.0, order=2, resample_dt_sec=0.0)
            elif v_scan < 300.0:
                base = dict(bw_hz=0.7, max_harm=8, notch_passes=4, trend_win_sec=5.0, order=2, resample_dt_sec=0.0)
            elif v_scan < 500.0:
                base = dict(bw_hz=0.6, max_harm=6, notch_passes=3, trend_win_sec=4.0, order=2, resample_dt_sec=0.0)
            elif v_scan < 800.0:
                base = dict(bw_hz=0.4, max_harm=4, notch_passes=2, trend_win_sec=3.0, order=2, resample_dt_sec=0.0)
            else:
                # Very fast scans: signal-protective baseline. (auto model may still decide to do 'none'.)
                base = dict(bw_hz=0.3, max_harm=3, notch_passes=2, trend_win_sec=3.0, order=2, resample_dt_sec=0.0)
        else:
            base = dict(RIPPLE_PRESETS["normal"])
    else:
        if preset not in ["safe", "normal", "strong"]:
            preset = "normal"
        base = dict(RIPPLE_PRESETS[preset])

    # Manual overrides (if provided)
    for k in ["bw_hz", "max_harm", "order", "notch_passes", "trend_win_sec", "resample_dt_sec"]:
        if pol.get(k, None) is not None:
            base[k] = pol[k]

    # Convert resample dt: None or <=0 => None (auto in remove_ripple_notch)
    dtv = base.get("resample_dt_sec", 0.0)
    if dtv is None:
        resample_dt = None
    else:
        try:
            dtf = float(dtv)
            resample_dt = None if (not np.isfinite(dtf) or dtf <= 0) else dtf
        except Exception:
            resample_dt = None

    # eval band: None or <=0 => None (use bw)
    ev = pol.get("eval_band_hz", None)
    eval_band = None
    if ev is not None:
        try:
            evf = float(ev)
            eval_band = None if (not np.isfinite(evf) or evf <= 0) else evf
        except Exception:
            eval_band = None

    return dict(
        target_hz=float(pol.get("target_hz", DEFAULT_RIPPLE_TARGET_HZ)),
        search_hz=float(pol.get("search_hz", DEFAULT_RIPPLE_SEARCH_HZ)),
        bw_hz=float(base.get("bw_hz", DEFAULT_RIPPLE_BW_HZ)),
        max_harm=int(base.get("max_harm", DEFAULT_RIPPLE_MAX_HARM)),
        trend_win_sec=float(base.get("trend_win_sec", DEFAULT_RIPPLE_TREND_WIN_SEC)),
        resample_dt_sec=resample_dt,
        order=int(base.get("order", DEFAULT_RIPPLE_ORDER)),
        notch_passes=int(base.get("notch_passes", DEFAULT_RIPPLE_NOTCH_PASS)),
        eval_band_hz=eval_band,
        model=model,
    )


# -----------------------------
# Spectral integration & loading
# -----------------------------
def integ_all_channels(spec2d: np.ndarray) -> np.ndarray:
    spec2d = np.asarray(spec2d, dtype=float)
    if spec2d.ndim != 2:
        raise ValueError(f"spectrum array must be 2D, got shape={spec2d.shape}")
    return np.sum(spec2d, axis=1)


def load_spectral(
    db,
    spectral_name: str,
    *,
    telescope: str = TELESCOPE,
    db_namespace: str = DEFAULT_DB_NAMESPACE,
):
    spec_table_name = f"{str(db_namespace)}-{str(telescope)}-data-spectral-{spectral_name}"
    arr_spec = _read_structured_array_tolerant(db, spec_table_name)
    t_spec, spec2d, time_meta, _s_field = _extract_spectral_from_structured(arr_spec)
    tp1 = integ_all_channels(spec2d)
    o = np.argsort(t_spec)
    return arr_spec, t_spec[o], tp1[o], o, time_meta


def _get_position_tags_from_structured(
    *,
    arr_spec,
    order: np.ndarray,
) -> np.ndarray:
    """Get position tags (HOT/OFF/ON) directly from necstdb spectral rows."""
    names = list(arr_spec.dtype.names or [])
    pos_name = _pick_field_name(names, "position", ["position", "pos", "obsmode", "mode"])
    if pos_name is None:
        raise RuntimeError(f"Cannot find position tags in spectral table. available={names}")

    pos_vals = np.asarray(arr_spec[pos_name])
    if pos_vals.shape[0] != order.shape[0]:
        raise RuntimeError(f"position length mismatch: {pos_vals.shape[0]} vs {order.shape[0]}")
    pos_vals = pos_vals[order]
    pos_str = np.array([_decode_position(v) for v in pos_vals], dtype=object)
    return pos_str


def _get_position_tags_with_fallback(
    db,
    *,
    arr_spec,
    order: np.ndarray,
    t_spec: np.ndarray,
    telescope: str = TELESCOPE,
    db_namespace: str = DEFAULT_DB_NAMESPACE,
) -> np.ndarray:
    """Get HOT/OFF/ON tags from spectral rows, with DB-table fallback.

    Preferred source is the necstdb spectral table itself. If the spectral rows do
    not contain usable position labels, fall back to common obsmode tables in the
    same DB namespace and align them to spectral timestamps.
    """
    try:
        pos_str = _get_position_tags_from_structured(arr_spec=arr_spec, order=order)
        good = np.array([_normalize_obsmode(v) for v in pos_str], dtype=object)
        if np.any(good == "HOT") or np.any(good == "OFF") or np.any(good == "ON"):
            return good
        print("[warn] spectral position tags exist but do not contain HOT/OFF/ON; trying DB obsmode fallback.")
    except Exception as e:
        print(f"[warn] cannot read position tags from spectral rows: {e}; trying DB obsmode fallback.")

    pos_db = _infer_obsmode_from_db(
        db,
        t_spec,
        telescope=telescope,
        db_namespace=db_namespace,
    )
    pos_str = np.array([_normalize_obsmode(v) for v in pos_db], dtype=object)
    return pos_str

# -----------------------------
# Chopper-wheel calibration (1-temperature, Ta*)
# -----------------------------
def _nearest_by_time(t_src: np.ndarray, v_src: np.ndarray, t_dst: np.ndarray) -> np.ndarray:
    """Nearest-neighbor assignment of values v_src at times t_src onto t_dst."""
    t_src = np.asarray(t_src, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)
    if t_src.ndim != 1 or t_dst.ndim != 1:
        raise ValueError("t_src and t_dst must be 1D")
    if t_src.size == 0:
        raise RuntimeError("empty t_src")
    o = np.argsort(t_src)
    t = t_src[o]
    v = v_src[o]
    idx = np.searchsorted(t, t_dst, side="left")
    idx0 = np.clip(idx - 1, 0, t.size - 1)
    idx1 = np.clip(idx, 0, t.size - 1)
    choose1 = np.abs(t[idx1] - t_dst) < np.abs(t[idx0] - t_dst)
    out = np.where(choose1, v[idx1], v[idx0])
    return out


def _normalize_obsmode(x) -> str:
    """Normalize observation mode labels to a compact uppercase string."""
    s = _decode_label(x)
    s = s.strip().upper()
    # common aliases
    if s in ["H", "HOT", "AMBIENT", "AMB", "LOAD", "CAL", "CALHOT", "CHOP", "CHOPOFF", "CHOPOUT"]:
        return "HOT"
    if s in ["O", "OFF", "SKY", "BLANK", "REF", "REFERENCE"]:
        return "OFF"
    if s in ["ON", "SRC", "SOURCE", "TARGET"]:
        return "ON"
    # if it contains these tokens
    if "HOT" in s or "AMB" in s or "LOAD" in s:
        return "HOT"
    if "OFF" in s or "SKY" in s:
        return "OFF"
    if "ON" in s:
        return "ON"
    return s if s else "UNKNOWN"


def _infer_obsmode_from_db(
    db,
    t_spec: np.ndarray,
    *,
    telescope: str = TELESCOPE,
    db_namespace: str = DEFAULT_DB_NAMESPACE,
) -> np.ndarray:
    """Try to read observation mode (HOT/OFF/ON) from common necstdb tables.

    Returns an array aligned to t_spec (nearest-neighbor in time).
    """
    # candidate tables and columns (best-effort; varies by deployment)
    table_candidates = [
        f"{str(db_namespace)}-{str(telescope)}-ctrl-obsmode",
        f"{str(db_namespace)}-{str(telescope)}-ctrl-observation-mode",
        f"{str(db_namespace)}-{str(telescope)}-ctrl-antenna-obsmode",
        f"{str(db_namespace)}-{str(telescope)}-ctrl-antenna-observation-mode",
    ]
    time_cols = ["time", "t", "timestamp"]
    mode_cols = ["obsmode", "OBSMODE", "mode", "MODE", "obs_mode", "OBS_MODE", "observation_mode"]

    for tb in table_candidates:
        dfm = _safe_read_table(db, tb)
        if dfm is None:
            continue
        tcol = None
        for c in time_cols:
            if c in dfm.columns:
                tcol = c
                break
        if tcol is None:
            continue
        mcol = None
        for c in mode_cols:
            if c in dfm.columns:
                mcol = c
                break
        if mcol is None:
            # sometimes the mode is stored in a single unnamed column; skip
            continue

        tt = pd.to_numeric(dfm[tcol], errors="coerce").to_numpy(float)
        mm = dfm[mcol].to_numpy()
        good = np.isfinite(tt)
        tt = tt[good]
        mm = mm[good]
        if tt.size < 2:
            continue
        tt, s_mode = _normalize_time_units(t_spec, tt, "obsmode")
        if s_mode != 1.0:
            print(f"[info] obsmode time scaled by {s_mode} to match unix seconds")
        return _nearest_by_time(tt, mm, np.asarray(t_spec, dtype=float))

    raise RuntimeError(
        "Cannot infer observation mode (HOT/OFF/ON). "
        "Tried common necstdb tables but none found with time+mode columns."
    )


def chopper_wheel_tastar(
    t_spec: np.ndarray,
    tp: np.ndarray,
    position: np.ndarray,
    *,
    tamb_k: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """1-temperature chopper-wheel calibration (1-load method): tp -> Ta* [K].

    User requirements (summary):
      - Identify HOT/OFF by *exact* matches on decoded position labels.
      - Compute total power tp (already given).
      - Interpolate HOT/OFF references to all t_spec using np.interp with edge-hold.
      - Ta* = Tamb * (tp - tp_off) / (tp_hot - tp_off)
      - If denom is non-finite or <=0, treat denom as NaN (avoid zero division).

    Returns
    -------
    ta_star : np.ndarray
        Ta* [K] aligned to t_spec.
    tp_hot_ref : np.ndarray
        Interpolated HOT reference total power aligned to t_spec.
    tp_off_ref : np.ndarray
        Interpolated OFF reference total power aligned to t_spec.
    """
    t_spec = np.asarray(t_spec, dtype=float)
    tp = np.asarray(tp, dtype=float)
    if t_spec.ndim != 1 or tp.ndim != 1 or t_spec.size != tp.size:
        raise ValueError("t_spec and tp must be 1D and same length.")
    if t_spec.size < 8:
        raise RuntimeError("Not enough samples for chopper-wheel calibration.")

    pos = np.asarray(position, dtype=object)
    if pos.size != t_spec.size:
        raise ValueError("position must have the same length as t_spec.")

    is_hot = np.array([p == "HOT" for p in pos], dtype=bool)
    is_off = np.array([p == "OFF" for p in pos], dtype=bool)
    n_hot = int(np.sum(is_hot))
    n_off = int(np.sum(is_off))
    if n_hot < 2 or n_off < 2:
        raise RuntimeError(f"Chopper-wheel needs >=2 HOT and >=2 OFF samples (got HOT={n_hot}, OFF={n_off}).")

    t_hot = t_spec[is_hot]
    tp_hot = tp[is_hot]
    t_off = t_spec[is_off]
    tp_off = tp[is_off]

    tp_hot_ref = _interp_edge_hold(t_hot, tp_hot, t_spec)
    tp_off_ref = _interp_edge_hold(t_off, tp_off, t_spec)

    denom = tp_hot_ref - tp_off_ref
    denom_bad = (~np.isfinite(denom)) | (denom <= 0)
    denom = denom.astype(float)
    denom[denom_bad] = np.nan

    tamb = float(tamb_k)
    if not np.isfinite(tamb) or tamb <= 0:
        raise ValueError("Tamb must be positive and finite.")

    ta_star = tamb * (tp - tp_off_ref) / denom
    return ta_star.astype(float), tp_hot_ref.astype(float), tp_off_ref.astype(float)



# -----------------------------
# Scan segmentation by spectral id
# -----------------------------
def make_id_index(id_values: np.ndarray) -> List[List[List[int]]]:
    """Build scan index from spectral id labels.

    Notes
    -----
    The previous implementation relied on a trailing sentinel-like label to flush
    the final AZ/EL pair. When the table ended without such a sentinel, the last
    scan pair was silently dropped. This version always flushes the final pending
    axis/cross blocks at end-of-loop.
    """
    id_index: List[List[List[int]]] = []
    cross_blocks: List[List[int]] = []
    axis_block: List[int] = []
    prev_cross: Optional[int] = None
    prev_axis: Optional[int] = None

    def _flush_axis() -> None:
        nonlocal axis_block, cross_blocks
        if len(axis_block) > 0:
            cross_blocks.append(axis_block)
            axis_block = []

    def _flush_cross() -> None:
        nonlocal cross_blocks, id_index
        _flush_axis()
        if len(cross_blocks) > 0:
            id_index.append(cross_blocks)
            cross_blocks = []

    for i, _x in enumerate(id_values):
        x = _decode_label(_x)
        if x == "":
            continue

        # Sentinel / trailer label: stop parsing, but keep the accumulated final scan.
        if ("-" not in x) and ("_" not in x):
            break

        try:
            axis = int(x[-1])
            cross = int(x[:-2])
        except Exception:
            continue

        if prev_cross is None:
            prev_cross = cross
            prev_axis = axis
        else:
            if cross != prev_cross:
                _flush_cross()
                prev_cross = cross
                prev_axis = axis
            elif axis != prev_axis:
                _flush_axis()
                prev_axis = axis

        axis_block.append(i)

    _flush_cross()
    return id_index


def build_scans_from_id(df: pd.DataFrame, id_values: np.ndarray) -> Tuple[Dict[int, pd.DataFrame], Dict[int, pd.DataFrame]]:
    id_index = make_id_index(id_values)
    az_scans: Dict[int, pd.DataFrame] = {}
    el_scans: Dict[int, pd.DataFrame] = {}

    for cross_id, pair in enumerate(id_index):
        if len(pair) < 2:
            continue
        if len(pair[0]) >= 2:
            i0, i1 = pair[0][0], pair[0][-1]
            seg = df.iloc[i0:i1 + 1].copy()
            if len(seg) > 0:
                az_scans[cross_id] = seg
        if len(pair[1]) >= 2:
            j0, j1 = pair[1][0], pair[1][-1]
            seg = df.iloc[j0:j1 + 1].copy()
            if len(seg) > 0:
                el_scans[cross_id] = seg

    return az_scans, el_scans


# -----------------------------
# Config helper (best effort)
# -----------------------------
def ensure_obs_and_config_files(rawdata_path: pathlib.Path) -> None:
    impl = list(rawdata_path.glob("*.toml"))
    if not impl:
        return

    config_dst = rawdata_path / "config.toml"
    for p in impl:
        src = rawdata_path / p.name
        if p.name.endswith("_config.toml"):
            if not config_dst.exists():
                try:
                    config_dst.write_bytes(src.read_bytes())
                except Exception:
                    pass
        elif p.name not in ["config.toml", "device_setting.toml", "pointing_param.toml"]:
            obs_dst = rawdata_path / f"{p.stem}.obs"
            if not obs_dst.exists():
                try:
                    obs_dst.write_bytes(src.read_bytes())
                except Exception:
                    pass


# -----------------------------
# Sun ephemeris at t_spec (astropy)
# -----------------------------
def sun_altaz_deg(
    t_unix: np.ndarray,
    rawdata_path: pathlib.Path,
    *,
    planet: str = PLANET,
    weather_meteo=None,
) -> Tuple[np.ndarray, np.ndarray]:
    os.environ["NECST_ROOT"] = str(rawdata_path)
    try:
        config.reload()
    except Exception:
        pass

    t = Time(t_unix, format="unix")

    if weather_meteo is not None and weather_meteo.get("press_hpa") is not None:
        try:
            press = np.asarray(weather_meteo.get("press_hpa"), dtype=float) * u.hPa
        except Exception:
            press = np.zeros_like(t_unix) * u.hPa
    else:
        press = np.zeros_like(t_unix) * u.hPa

    if weather_meteo is not None and weather_meteo.get("temp_c") is not None:
        try:
            temp = np.asarray(weather_meteo.get("temp_c"), dtype=float) * u.deg_C
        except Exception:
            temp = np.zeros_like(t_unix) * u.deg_C
    else:
        temp = np.zeros_like(t_unix) * u.deg_C

    if weather_meteo is not None and weather_meteo.get("humid_rh") is not None:
        try:
            humid = np.asarray(weather_meteo.get("humid_rh"), dtype=float)
        except Exception:
            humid = np.zeros_like(t_unix)
    else:
        humid = np.zeros_like(t_unix)

    try:
        obswl = const.c / config.observation_frequency
    except Exception:
        obswl = const.c / (100.0 * u.GHz)

    ref = get_body(str(planet), location=config.location, time=t)
    altaz = AltAz(obstime=t, location=config.location, pressure=press, temperature=temp, relative_humidity=humid, obswl=obswl)
    r = ref.transform_to(altaz)
    return r.az.to_value(u.deg), r.alt.to_value(u.deg)


# -----------------------------
# Motion-based trimming
# -----------------------------
def trim_scan_by_motion(
    seg: pd.DataFrame,
    xcol: str,
    *,
    enabled: bool,
    vfrac: float,
    vmin: float,
    gap_fill: int,
    min_samples: int,
) -> pd.DataFrame:
    if not enabled:
        return seg
    if xcol not in seg.columns or "t_unix" not in seg.columns:
        raise RuntimeError(f"trim_scan_by_motion requires columns: {xcol}, t_unix")

    t = seg["t_unix"].to_numpy(float)
    x = seg[xcol].to_numpy(float)
    if len(t) < 5:
        return seg

    dt = np.diff(t)
    dx = np.diff(x)
    good_dt = np.isfinite(dt) & (dt > 0) & np.isfinite(dx)

    vabs = np.full_like(dx, np.nan, dtype=float)
    vabs[good_dt] = np.abs(dx[good_dt] / dt[good_dt])

    vv = vabs[np.isfinite(vabs)]
    if vv.size < 10:
        return seg

    p90 = float(np.nanpercentile(vv, 90))
    thr = builtins.max(float(vmin), float(vfrac) * p90)

    good = np.isfinite(vabs) & (vabs >= thr)

    n = len(x)
    pm = np.zeros(n, dtype=bool)
    pm[1:] |= good
    pm[:-1] |= good

    # fill small gaps
    gap = int(gap_fill)
    if gap > 0:
        i = 0
        while i < n:
            if pm[i]:
                i += 1
                continue
            j = i
            while j < n and (not pm[j]):
                j += 1
            if i > 0 and j < n and (j - i) <= gap and pm[i - 1] and pm[j]:
                pm[i:j] = True
            i = j

    # largest contiguous True block
    best_len = 0
    best = (0, n)
    i = 0
    while i < n:
        if not pm[i]:
            i += 1
            continue
        j = i
        while j < n and pm[j]:
            j += 1
        if (j - i) > best_len:
            best_len = j - i
            best = (i, j)
        i = j

    if best_len < int(min_samples):
        raise RuntimeError(
            f"Trimmed scan too short: {best_len} samples (<{min_samples}). "
            f"Try --no-trim-scan or lower trim thresholds."
        )

    return seg.iloc[best[0]:best[1]].copy()

def trim_scan_by_dominant_axis(
    seg: pd.DataFrame,
    main_xcol: str,
    cross_xcol: str,
    *,
    enabled: bool,
    vfrac: float,
    vmin: float,
    gap_fill: int,
    min_samples: int,
    ratio_min: float = DEFAULT_TRIM_AXIS_RATIO_MIN,
    vpercentile: float = DEFAULT_TRIM_VPERCENTILE,
) -> pd.DataFrame:
    """Trim to the most likely *scan* segment by requiring dominant motion in main axis.

    Motivation
    ----------
    The legacy trim (1D |Δx/Δt| threshold) can pick slews where both axes move
    (e.g., moving to scan start), which is not the intended constant-speed scan.
    This function keeps intervals where:
        - |v_main| is sufficiently large, and
        - |v_main| / (|v_cross| + eps) is sufficiently large.

    Selection
    ---------
    Among contiguous candidate blocks, choose the block with the largest *span*
    in main axis (robust when true scan is short but covers a larger range).
    """
    if not enabled:
        return seg
    if (main_xcol not in seg.columns) or (cross_xcol not in seg.columns) or ("t_unix" not in seg.columns):
        raise RuntimeError(f"trim_scan_by_dominant_axis requires columns: {main_xcol}, {cross_xcol}, t_unix")

    t = seg["t_unix"].to_numpy(float)
    xm = seg[main_xcol].to_numpy(float)
    xc = seg[cross_xcol].to_numpy(float)
    n = len(t)
    if n < 8:
        return seg

    dt = np.diff(t)
    dxm = np.diff(xm)
    dxc = np.diff(xc)
    good_dt = np.isfinite(dt) & (dt > 0) & np.isfinite(dxm) & np.isfinite(dxc)

    v_main = np.full_like(dxm, np.nan, dtype=float)
    v_cross = np.full_like(dxc, np.nan, dtype=float)
    v_main[good_dt] = np.abs(dxm[good_dt] / dt[good_dt])
    v_cross[good_dt] = np.abs(dxc[good_dt] / dt[good_dt])

    vv = v_main[np.isfinite(v_main)]
    if vv.size < 10:
        # fallback to legacy
        return trim_scan_by_motion(seg, main_xcol, enabled=True, vfrac=vfrac, vmin=vmin, gap_fill=gap_fill, min_samples=min_samples)

    vp = float(np.nanpercentile(vv, float(vpercentile)))
    thr_main = builtins.max(float(vmin), float(vfrac) * vp)

    eps = 1e-12
    ratio = v_main / (v_cross + eps)
    good = np.isfinite(v_main) & (v_main >= thr_main) & np.isfinite(ratio) & (ratio >= float(ratio_min))

    pm = np.zeros(n, dtype=bool)
    pm[1:] |= good
    pm[:-1] |= good

    # fill small gaps
    gap = int(gap_fill)
    if gap > 0:
        i = 0
        while i < n:
            if pm[i]:
                i += 1
                continue
            j = i
            while j < n and (not pm[j]):
                j += 1
            if i > 0 and j < n and (j - i) <= gap and pm[i - 1] and pm[j]:
                pm[i:j] = True
            i = j

    # choose contiguous block with max span in main axis (then length)
    best = None
    best_span = -1.0
    best_len = 0
    i = 0
    while i < n:
        if not pm[i]:
            i += 1
            continue
        j = i
        while j < n and pm[j]:
            j += 1
        xm_seg = xm[i:j]
        if xm_seg.size > 0 and np.any(np.isfinite(xm_seg)):
            span = float(np.nanmax(xm_seg) - np.nanmin(xm_seg))
        else:
            span = 0.0
        L = j - i
        if (span > best_span) or (span == best_span and L > best_len):
            best_span = span
            best_len = L
            best = (i, j)
        i = j

    if best is None or best_len < int(min_samples):
        # fallback to legacy (keeps behavior if dominant-axis selection fails)
        return trim_scan_by_motion(seg, main_xcol, enabled=True, vfrac=vfrac, vmin=vmin, gap_fill=gap_fill, min_samples=min_samples)

    return seg.iloc[best[0]:best[1]].copy()



def trim_scan_by_on_constant_speed_near_sun(
    seg: pd.DataFrame,
    main_xcol: str,
    cross_xcol: str,
    *,
    enabled: bool,
    ratio_min: float,
    speed_min_deg_s: float,
    profile_xlim_deg: float,
    xwin_factor: float,
    cross_offset_max_deg: float,
    gap_fill: int,
    min_samples: int,
    use_on_only: bool = True,
    steady_cv_max: float = DEFAULT_TRIM_STEADY_CV_MAX,
) -> pd.DataFrame:
    """Select the most likely constant-speed *ON-source* scan segment near the Sun.

    Key requirements (motivated by your data)
    ----------------------------------------
    - Ignore HOT/OFF samples even when spectral id matches (use position=="ON" only).
    - Prefer segments where main-axis motion dominates cross-axis motion (ratio criterion).
    - Prefer near-Sun part (|offset_main| within profile_xlim_deg * xwin_factor).
    - Reject near-stationary drift using speed_min_deg_s.
    - Choose the best contiguous block by (main span) and speed stability.

    Notes
    -----
    This is intended to avoid picking the pre-scan slew (move-to-start) that can satisfy
    simple |Δx/Δt| thresholds. When this selection fails, callers may fall back to
    dominant-axis or legacy trimming.
    """
    if not enabled:
        return seg
    if ("t_unix" not in seg.columns) or (main_xcol not in seg.columns) or (cross_xcol not in seg.columns):
        raise RuntimeError(f"trim_scan_by_on_constant_speed_near_sun requires columns: t_unix, {main_xcol}, {cross_xcol}")

    t = seg["t_unix"].to_numpy(float)
    xm = seg[main_xcol].to_numpy(float)
    xc = seg[cross_xcol].to_numpy(float)
    n = int(len(t))
    if n < 8:
        raise RuntimeError("segment too short for steady-scan selection")

    pos = None
    if use_on_only and ("position" in seg.columns):
        pos = np.asarray(seg["position"], dtype=object)
        on = np.array([(_decode_position(v) == "ON") for v in pos], dtype=bool)
    else:
        on = np.ones(n, dtype=bool)

    # near-Sun window on main axis
    xl = float(profile_xlim_deg) if (profile_xlim_deg is not None and np.isfinite(profile_xlim_deg) and profile_xlim_deg > 0) else np.inf
    xwin = (xl * float(xwin_factor)) if np.isfinite(xl) else np.inf
    near = np.isfinite(xm) & (np.abs(xm) <= xwin)

    # cross-axis stability: allow constant offset but exclude large excursions
    m_for_cross = on & near & np.isfinite(xc)
    if np.count_nonzero(m_for_cross) >= 5:
        c0 = float(np.nanmedian(xc[m_for_cross]))
    else:
        c0 = float(np.nanmedian(xc[np.isfinite(xc)])) if np.any(np.isfinite(xc)) else 0.0
    cross_ok = np.isfinite(xc) & (np.abs(xc - c0) <= float(cross_offset_max_deg))

    # interval speeds
    dt = np.diff(t)
    dxm = np.diff(xm)
    dxc = np.diff(xc)

    good_dt = np.isfinite(dt) & (dt > 0) & np.isfinite(dxm) & np.isfinite(dxc)
    v_main = np.full_like(dxm, np.nan, dtype=float)
    v_cross = np.full_like(dxc, np.nan, dtype=float)
    v_main[good_dt] = dxm[good_dt] / dt[good_dt]
    v_cross[good_dt] = dxc[good_dt] / dt[good_dt]

    abs_vm = np.abs(v_main)
    abs_vc = np.abs(v_cross)
    eps = 1e-12
    ratio = abs_vm / (abs_vc + eps)

    # interval-level ON/near/cross_ok mask (both endpoints must satisfy)
    on_int = on[:-1] & on[1:]
    near_int = near[:-1] & near[1:]
    cross_int = cross_ok[:-1] & cross_ok[1:]
    base_int = good_dt & on_int & near_int & cross_int & np.isfinite(ratio)

    good_int = base_int & (abs_vm >= float(speed_min_deg_s)) & (ratio >= float(ratio_min))

    # convert to point mask
    pm = np.zeros(n, dtype=bool)
    pm[:-1] |= good_int
    pm[1:] |= good_int

    # fill small gaps
    gap = int(gap_fill)
    if gap > 0:
        i = 0
        while i < n:
            if pm[i]:
                i += 1
                continue
            j = i
            while j < n and (not pm[j]):
                j += 1
            if i > 0 and j < n and (j - i) <= gap and pm[i - 1] and pm[j]:
                pm[i:j] = True
            i = j

    # choose best contiguous block: maximize main span, then length; require stable speed
    best = None
    best_span = -1.0
    best_len = 0

    i = 0
    while i < n:
        if not pm[i]:
            i += 1
            continue
        j = i
        while j < n and pm[j]:
            j += 1
        L = j - i
        if L >= int(min_samples):
            xm_seg = xm[i:j]
            # robust span (ignore edge outliers)
            if np.any(np.isfinite(xm_seg)):
                span = float(np.nanpercentile(xm_seg, 95) - np.nanpercentile(xm_seg, 5))
            else:
                span = 0.0

            # speed stability within this block
            vm_seg = v_main[i:j-1]  # interval speeds inside
            vm_seg = vm_seg[np.isfinite(vm_seg)]
            if vm_seg.size >= 5:
                v_med = float(np.nanmedian(vm_seg))
                v_abs = abs(v_med)
                if v_abs >= float(speed_min_deg_s):
                    # sign consistency
                    s0 = 1.0 if v_med >= 0 else -1.0
                    frac_same = float(np.mean(np.sign(vm_seg) == s0))
                    # robust CV
                    mad = float(np.nanmedian(np.abs(vm_seg - v_med)))
                    sigma = 1.4826 * mad
                    cv = float(sigma / (v_abs + 1e-12))
                    if (frac_same >= 0.8) and (cv <= float(steady_cv_max)):
                        if (span > best_span) or (span == best_span and L > best_len):
                            best_span = span
                            best_len = L
                            best = (i, j)
        i = j

    if best is None:
        raise RuntimeError("steady-scan selection failed (no valid ON constant-speed block found near Sun)")

    return seg.iloc[best[0]:best[1]].copy()


def trim_scan_by_on_xwindow_speed(
    seg: pd.DataFrame,
    main_xcol: str,
    *,
    enabled: bool,
    profile_xlim_deg: float,
    speed_min_deg_s: float,
    gap_fill: int,
    min_samples: int,
    use_on_only: bool = True,
) -> pd.DataFrame:
    """Simple scan-region selection: ON-only & |offset_main| <= xlim, then keep moving parts.

    Motivation
    ----------
    Simplify scan detection by defining the scan region as the neighborhood of the
    expected Sun position (offset ~ 0) within +/- profile_xlim_deg.

    Minimal safeguard
    -----------------
    Require |v_main| >= speed_min_deg_s so that long tracking/drift segments inside
    the window do not dominate the selection.
    """
    if not enabled:
        return seg
    if ("t_unix" not in seg.columns) or (main_xcol not in seg.columns):
        raise RuntimeError(f"trim_scan_by_on_xwindow_speed requires columns: t_unix, {main_xcol}")

    t = seg["t_unix"].to_numpy(float)
    xm = seg[main_xcol].to_numpy(float)
    n = int(len(t))
    if n < 8:
        raise RuntimeError("segment too short for x-window scan selection")

    # ON-only mask
    if use_on_only and ("position" in seg.columns):
        pos = np.asarray(seg["position"], dtype=object)
        on = np.array([(_decode_position(v) == "ON") for v in pos], dtype=bool)
    else:
        on = np.ones(n, dtype=bool)

    # x-window around offset==0
    xl = float(profile_xlim_deg)
    if (not np.isfinite(xl)) or (xl <= 0):
        xl = np.inf
    near = np.isfinite(xm) & (np.abs(xm) <= xl)

    # interval speeds
    dt = np.diff(t)
    dxm = np.diff(xm)
    good_dt = np.isfinite(dt) & (dt > 0) & np.isfinite(dxm)
    v_main = np.full_like(dxm, np.nan, dtype=float)
    v_main[good_dt] = dxm[good_dt] / dt[good_dt]
    abs_vm = np.abs(v_main)

    on_int = on[:-1] & on[1:]
    near_int = near[:-1] & near[1:]
    good_int = good_dt & on_int & near_int & (abs_vm >= float(speed_min_deg_s))

    # point mask
    pm = np.zeros(n, dtype=bool)
    pm[:-1] |= good_int
    pm[1:] |= good_int

    # fill small gaps
    gap = int(gap_fill)
    if gap > 0:
        i = 0
        while i < n:
            if pm[i]:
                i += 1
                continue
            j = i
            while j < n and (not pm[j]):
                j += 1
            if i > 0 and j < n and (j - i) <= gap and pm[i - 1] and pm[j]:
                pm[i:j] = True
            i = j

    # choose best block by span then length
    best = None
    best_span = -1.0
    best_len = 0

    i = 0
    while i < n:
        if not pm[i]:
            i += 1
            continue
        j = i
        while j < n and pm[j]:
            j += 1
        L = j - i
        if L >= int(min_samples):
            xm_seg = xm[i:j]
            if np.any(np.isfinite(xm_seg)):
                span = float(np.nanpercentile(xm_seg, 95) - np.nanpercentile(xm_seg, 5))
            else:
                span = 0.0
            if (span > best_span) or (span == best_span and L > best_len):
                best_span = span
                best_len = L
                best = (i, j)
        i = j

    if best is None:
        raise RuntimeError("x-window scan selection failed (no ON moving block found within +/-xlim)")

    return seg.iloc[best[0]:best[1]].copy()


def trim_scan_segment(
    seg: pd.DataFrame,
    main_xcol: str,
    cross_xcol: str,
    *,
    enabled: bool,
    trim_params: dict,
) -> pd.DataFrame:
    """Unified trimming entry point used throughout this script."""
    if not enabled:
        return seg
    vfrac = float(trim_params.get("vfrac", DEFAULT_TRIM_VFRAC))
    vmin = float(trim_params.get("vmin", DEFAULT_TRIM_VMIN))
    gap_fill = int(trim_params.get("gap_fill", DEFAULT_TRIM_GAP))
    min_samples = int(trim_params.get("min_samples", DEFAULT_TRIM_MIN_SAMPLES))
    ratio_min = float(trim_params.get("ratio_min", DEFAULT_TRIM_AXIS_RATIO_MIN))
    vpercentile = float(trim_params.get("vpercentile", DEFAULT_TRIM_VPERCENTILE))

    # --- NEW: prefer ON-only steady-scan selection near Sun (avoids slews and HOT/OFF) ---
    profile_xlim_deg = float(trim_params.get("profile_xlim_deg", DEFAULT_PROFILE_XLIM_DEG))
    xwin_factor = float(trim_params.get("xwin_factor", DEFAULT_TRIM_XWIN_FACTOR))
    cross_offset_max_deg = float(trim_params.get("cross_offset_max_deg", DEFAULT_TRIM_CROSS_OFFSET_MAX_DEG))
    speed_min_deg_s = float(trim_params.get("speed_min_deg_s", DEFAULT_TRIM_SCAN_SPEED_MIN_ARCSEC_S / 3600.0))
    use_on_only = bool(trim_params.get("use_on_only", DEFAULT_TRIM_USE_ON_ONLY))
    steady_scan = bool(trim_params.get("steady_scan", DEFAULT_TRIM_STEADY_SCAN))
    steady_cv_max = float(trim_params.get("steady_cv_max", DEFAULT_TRIM_STEADY_CV_MAX))

    if steady_scan:
        try:
            # --- simplest per your request: use +/- profile_xlim_deg around 0 (ON-only) ---
            return trim_scan_by_on_xwindow_speed(
                seg,
                main_xcol,
                enabled=True,
                profile_xlim_deg=profile_xlim_deg,
                speed_min_deg_s=speed_min_deg_s,
                gap_fill=gap_fill,
                min_samples=min_samples,
                use_on_only=use_on_only,
            )
        except Exception:
            # fall back to the more constrained steady-scan selection logic
            pass

        try:
            return trim_scan_by_on_constant_speed_near_sun(
                seg,
                main_xcol,
                cross_xcol,
                enabled=True,
                ratio_min=ratio_min,
                speed_min_deg_s=speed_min_deg_s,
                profile_xlim_deg=profile_xlim_deg,
                xwin_factor=xwin_factor,
                cross_offset_max_deg=cross_offset_max_deg,
                gap_fill=gap_fill,
                min_samples=min_samples,
                use_on_only=use_on_only,
                steady_cv_max=steady_cv_max,
            )
        except Exception:
            # fall back to the legacy trimming logic
            pass

    # --- fallback: dominant-axis trimming (can still reject slews) ---
    if bool(trim_params.get("dominant_axis", DEFAULT_TRIM_DOMINANT_AXIS)):
        return trim_scan_by_dominant_axis(
            seg,
            main_xcol,
            cross_xcol,
            enabled=True,
            vfrac=vfrac,
            vmin=vmin,
            gap_fill=gap_fill,
            min_samples=min_samples,
            ratio_min=ratio_min,
            vpercentile=vpercentile,
        )

    # --- last resort: 1D motion trimming (may pick slews in some datasets) ---
    return trim_scan_by_motion(seg, main_xcol, enabled=True, vfrac=vfrac, vmin=vmin, gap_fill=gap_fill, min_samples=min_samples)


# -----------------------------
# Finite difference derivative (skip duplicates)
# -----------------------------
def finite_difference_skip_duplicates(
    x: np.ndarray,
    y: np.ndarray,
    *,
    eps: float = None,
    jump_max_deg: float = 1.0,
    strict: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = int(x.size)
    if x.ndim != 1 or y.ndim != 1 or n != y.size:
        raise RuntimeError("derivative requires 1D x/y with same length.")
    if n < 3:
        raise RuntimeError("Need at least 3 points for derivative.")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise RuntimeError("Non-finite x/y encountered.")

    dx1 = np.diff(x)
    if np.any(~np.isfinite(dx1)):
        raise RuntimeError("Non-finite dx encountered.")
    if np.any(np.abs(dx1) > float(jump_max_deg)):
        raise RuntimeError(f"Too large x jump detected (>|{jump_max_deg}| deg). Possible wrap/glitch.")

    if eps is None:
        adx = np.abs(dx1)
        adx = adx[adx > 0]
        med = float(np.median(adx)) if adx.size else 0.0
        eps = builtins.max(1e-8, 1e-3 * med)
    eps = float(eps)

    dydx = np.full(n, np.nan, dtype=float)

    def _find_left(i: int) -> Optional[int]:
        xi = x[i]
        j = i - 1
        while j >= 0:
            if abs(xi - x[j]) >= eps:
                return j
            j -= 1
        return None

    def _find_right(i: int) -> Optional[int]:
        xi = x[i]
        k = i + 1
        while k < n:
            if abs(x[k] - xi) >= eps:
                return k
            k += 1
        return None

    for i in range(1, n - 1):
        j = _find_left(i)
        k = _find_right(i)
        if j is None or k is None:
            continue
        den = x[k] - x[j]
        if abs(den) < eps:
            continue
        dydx[i] = (y[k] - y[j]) / den

    k0 = _find_right(0)
    if k0 is not None and abs(x[k0] - x[0]) >= eps:
        dydx[0] = (y[k0] - y[0]) / (x[k0] - x[0])
    jn = _find_left(n - 1)
    if jn is not None and abs(x[n - 1] - x[jn]) >= eps:
        dydx[n - 1] = (y[n - 1] - y[jn]) / (x[n - 1] - x[jn])

    if strict and np.any(~np.isfinite(dydx)):
        n_nan = int(np.sum(~np.isfinite(dydx)))
        raise RuntimeError(f"Derivative undefined at {n_nan}/{n} points due to repeated x or insufficient span (dx≈0).")

    return x, dydx


# -----------------------------
# Data acquisition: RawData -> df
# -----------------------------
def _safe_int(value, default=-1) -> int:
    try:
        return int(value)
    except Exception:
        return int(default)


def _pick_tracking_point(seg: pd.DataFrame, main_xcol: str, *, trim_enabled: bool, trim_params: dict) -> dict:
    """Pick the point where the telescope most nearly thinks it is on the Sun.

    Operational definition:
      - trim to the usable scan segment (if enabled)
      - within that segment, choose the sample with minimum |main-axis offset|
        (AZ scan -> |daz| minimum, EL scan -> |d_el| minimum)

    This avoids the ambiguity of a simple midpoint and matches the requested
    point where the telescope believes it is looking at the Sun.
    """
    cross_xcol = "d_el" if str(main_xcol) == "daz" else "daz"
    seg2 = trim_scan_segment(seg, main_xcol, cross_xcol, enabled=trim_enabled, trim_params=trim_params) if trim_enabled else seg

    out = {
        "n_total": _safe_int(len(seg), 0),
        "n_used": _safe_int(len(seg2), 0),
        "t_unix": float("nan"),
        "az_true_deg": float("nan"),
        "el_true_deg": float("nan"),
        "az_sun_deg": float("nan"),
        "el_sun_deg": float("nan"),
        "main_offset_deg": float("nan"),
        "cross_offset_deg": float("nan"),
    }

    if len(seg2) == 0:
        return out

    x = seg2[main_xcol].to_numpy(float)
    good = np.isfinite(x)
    for col in ["t_unix", "az_true", "el_true", "az_sun", "el_sun"]:
        if col in seg2.columns:
            good &= np.isfinite(seg2[col].to_numpy(float))
    if cross_xcol in seg2.columns:
        cross = seg2[cross_xcol].to_numpy(float)
        good &= np.isfinite(cross)

    if not np.any(good):
        return out

    idx_good = np.where(good)[0]
    k = idx_good[int(np.argmin(np.abs(x[good])))]
    row = seg2.iloc[int(k)]

    out.update({
        "t_unix": float(row.get("t_unix", np.nan)),
        "az_true_deg": float(row.get("az_true", np.nan)),
        "el_true_deg": float(row.get("el_true", np.nan)),
        "az_sun_deg": float(row.get("az_sun", np.nan)),
        "el_sun_deg": float(row.get("el_sun", np.nan)),
        "main_offset_deg": float(row.get(main_xcol, np.nan)),
        "cross_offset_deg": float(row.get(cross_xcol, np.nan)),
    })
    return out


def _nanmean_pair(*vals: float) -> float:
    arr = np.asarray(vals, dtype=float)
    return float(np.nanmean(arr)) if np.any(np.isfinite(arr)) else float("nan")



def build_dataframe(
    rawdata_path: pathlib.Path,
    spectral_name: str,
    *,
    telescope: str = TELESCOPE,
    db_namespace: str = DEFAULT_DB_NAMESPACE,
    tel_loaddata: str = TEL_LOADDATA,
    planet: str = PLANET,
    azel_source: str,
    altaz_apply: Optional[str] = None,
    azel_correction_apply: Optional[str] = None,
    encoder_table: Optional[str] = None,
    encoder_table_suffix: str = "ctrl-antenna-encoder",
    altaz_table: Optional[str] = None,
    altaz_table_suffix: str = "ctrl-antenna-altaz",
    encoder_time_col: str = "time",
    altaz_time_col: str = "time",
    spectrometer_time_offset_sec: float = 0.0,
    encoder_shift_sec: float = 0.0,
    encoder_az_time_offset_sec: float = 0.0,
    encoder_el_time_offset_sec: float = 0.0,
    encoder_vavg_sec: float = 0.0,
    chopper_wheel: bool = True,
    tamb_k: Optional[float] = None,
    tamb_default_k: float = 300.0,
    tamb_min_k: float = 250.0,
    tamb_max_k: float = 330.0,
    weather_inside_table: Optional[str] = None,
    weather_inside_table_suffix: str = "weather-ambient",
    weather_inside_time_col: str = "time",
    weather_outside_table: Optional[str] = None,
    weather_outside_table_suffix: str = "weather-ambient",
    weather_outside_time_col: str = "time",
    outside_default_temperature_c: float = 0.0,
    outside_default_pressure_hpa: float = 760.0,
    outside_default_humidity_pct: float = 30.0,
    outside_temperature_min_c: float = -50.0,
    outside_temperature_max_c: float = 50.0,
    outside_pressure_min_hpa: float = 400.0,
    outside_pressure_max_hpa: float = 1100.0,
    outside_humidity_min_pct: float = 0.0,
    outside_humidity_max_pct: float = 100.0,
    chopper_win_sec: float = 5.0,
    chopper_stat: str = "median",
) -> Tuple[pd.DataFrame, Dict[int, pd.DataFrame], Dict[int, pd.DataFrame]]:
    ensure_obs_and_config_files(rawdata_path)

    db = necstdb.opendb(rawdata_path)
    arr_spec, t_spec, tp1, order, spec_time_meta = load_spectral(
        db,
        spectral_name,
        telescope=telescope,
        db_namespace=db_namespace,
    )
    t_spec = np.asarray(t_spec, dtype=float) + float(spectrometer_time_offset_sec)

    names_spec = list(arr_spec.dtype.names or [])
    id_valid = False
    id_name = _pick_field_name(names_spec, "id", ["id", "scan_id"])
    if id_name is not None:
        id_vals = np.asarray(arr_spec[id_name])
        if id_vals.shape[0] == order.shape[0]:
            id_vals = id_vals[order]
            id_valid = bool(np.any([_decode_label(v) != "" for v in id_vals]))
            if not id_valid:
                print("[warn] spectral-table id exists but all labels are empty; encoder smoothing will be disabled.")
        else:
            print(f"[warn] spectral-table id length mismatch: {id_vals.shape[0]} vs {order.shape[0]}; encoder smoothing will be disabled.")
            id_vals = np.array([b""] * len(t_spec), dtype=object)
    else:
        id_vals = np.array([b""] * len(t_spec), dtype=object)

    position_str = _get_position_tags_with_fallback(
        db,
        arr_spec=arr_spec,
        order=order,
        t_spec=t_spec,
        telescope=telescope,
        db_namespace=db_namespace,
    )
    src = str(azel_source).lower().strip()
    azel_apply_mode = _normalize_azel_correction_apply(
        azel_correction_apply if azel_correction_apply is not None else altaz_apply,
        default=("subtract" if src == "encoder" else "none"),
    )
    encoder_table_resolved = _resolve_table_name(db_namespace=db_namespace, telescope=telescope, full_table=encoder_table, suffix=encoder_table_suffix, default_suffix="ctrl-antenna-encoder")
    altaz_table_resolved = _resolve_table_name(db_namespace=db_namespace, telescope=telescope, full_table=altaz_table, suffix=altaz_table_suffix, default_suffix="ctrl-antenna-altaz")
    weather_inside_table_resolved = _resolve_table_name(db_namespace=db_namespace, telescope=telescope, full_table=weather_inside_table, suffix=weather_inside_table_suffix, default_suffix="weather-ambient")
    weather_outside_table_resolved = _resolve_table_name(db_namespace=db_namespace, telescope=telescope, full_table=weather_outside_table, suffix=weather_outside_table_suffix, default_suffix="weather-ambient")
    weather_inside_meteo = _load_weather_meteo_for_tspec_from_table(db, weather_inside_table_resolved, t_spec, preferred_time_col=weather_inside_time_col)
    if weather_outside_table_resolved == weather_inside_table_resolved and str(weather_outside_time_col).strip() == str(weather_inside_time_col).strip():
        weather_outside_raw = weather_inside_meteo
    else:
        weather_outside_raw = _load_weather_meteo_for_tspec_from_table(db, weather_outside_table_resolved, t_spec, preferred_time_col=weather_outside_time_col)
    weather_outside_meteo = _finalize_outside_weather_meteo(
        weather_outside_raw,
        t_spec,
        default_temperature_c=float(outside_default_temperature_c),
        default_pressure_hpa=float(outside_default_pressure_hpa),
        default_humidity_pct=float(outside_default_humidity_pct),
        temperature_min_c=float(outside_temperature_min_c),
        temperature_max_c=float(outside_temperature_max_c),
        pressure_min_hpa=float(outside_pressure_min_hpa),
        pressure_max_hpa=float(outside_pressure_max_hpa),
        humidity_min_pct=float(outside_humidity_min_pct),
        humidity_max_pct=float(outside_humidity_max_pct),
    )

    # Tamb: user provided OR estimate from inside weather/spectral necstdb fields else fallback.
    if tamb_k is None or (isinstance(tamb_k, float) and (not np.isfinite(tamb_k))):
        tamb_use, tamb_source = _estimate_tamb_k_with_source(
            weather_meteo=weather_inside_meteo,
            structured_arr=arr_spec,
            fallback_k=float(tamb_default_k),
            min_k=float(tamb_min_k),
            max_k=float(tamb_max_k),
        )
    else:
        tamb_use = float(tamb_k)
        tamb_source = "user"
    if (not np.isfinite(tamb_use)) or (tamb_use <= 0):
        tamb_use = float(tamb_default_k)
        tamb_source = "default"

    if bool(chopper_wheel):
        ta_star, tp_hot_ref, tp_off_ref = chopper_wheel_tastar(
            t_spec,
            tp1,
            position_str,
            tamb_k=float(tamb_use),
        )
    else:
        ta_star = tp1.copy().astype(float)
        tp_hot_ref = np.full_like(tp1, np.nan, dtype=float)
        tp_off_ref = np.full_like(tp1, np.nan, dtype=float)


    # altaz table (always needed for dlon/dlat and optional lon/lat)
    alt_name = altaz_table_resolved
    df_alt = _safe_read_table(db, alt_name)
    if df_alt is None:
        raise RuntimeError(f"cannot read altaz table: {alt_name}")
    alt_tcol = _sanitize_time_col(altaz_time_col, label="altaz_time_col")
    for k in [alt_tcol, "dlon", "dlat", "lon", "lat"]:
        if k not in df_alt.columns:
            raise RuntimeError(f"altaz table '{alt_name}' missing column '{k}'")
    t_alt, alt = _dedup_mean(df_alt, alt_tcol, ["lon", "lat", "dlon", "dlat"])
    t_alt, s_alt = _normalize_time_units(t_spec, t_alt, "altaz")
    if s_alt != 1.0:
        print(f"[info] altaz   time scaled by {s_alt} to match unix seconds")
    lon_alt = _interp_az_deg(t_alt, alt["lon"].to_numpy(float), t_spec)
    lat_alt = _interp_lin(t_alt, alt["lat"].to_numpy(float), t_spec)
    dlon = _interp_lin(t_alt, alt["dlon"].to_numpy(float), t_spec)
    dlat = _interp_lin(t_alt, alt["dlat"].to_numpy(float), t_spec)

    az_enc = np.full_like(t_spec, np.nan, dtype=float)
    el_enc = np.full_like(t_spec, np.nan, dtype=float)
    az_enc_raw = np.full_like(t_spec, np.nan, dtype=float)
    el_enc_raw = np.full_like(t_spec, np.nan, dtype=float)

    if src == "altaz":
        # command stream
        az_true, el_true = _apply_azel_correction(lon_alt, lat_alt, dlon, dlat, azel_apply_mode)
    else:
        # encoder stream
        enc_name = encoder_table_resolved
        df_enc = _safe_read_table(db, enc_name)
        if df_enc is None:
            raise RuntimeError(f"cannot read encoder table: {enc_name}")
        enc_tcol = _sanitize_time_col(encoder_time_col, label="encoder_time_col")
        for k in [enc_tcol, "lon", "lat"]:
            if k not in df_enc.columns:
                raise RuntimeError(f"encoder table '{enc_name}' missing column '{k}'")
        t_enc, enc = _dedup_mean(df_enc, enc_tcol, ["lon", "lat"])
        t_enc, s_enc = _normalize_time_units(t_spec, t_enc, "encoder")
        if s_enc != 1.0:
            print(f"[info] encoder time scaled by {s_enc} to match unix seconds")
        t_enc_base = t_enc + float(encoder_shift_sec)
        t_enc_az = t_enc_base + float(encoder_az_time_offset_sec)
        t_enc_el = t_enc_base + float(encoder_el_time_offset_sec)

        lon_enc = enc["lon"].to_numpy(float)
        lat_enc = enc["lat"].to_numpy(float)

        # raw interpolation (diagnostic)
        az_enc_raw = _interp_az_deg(t_enc_az, lon_enc, t_spec)
        el_enc_raw = _interp_lin(t_enc_el, lat_enc, t_spec)

        # safe local smoothing on interpolated positions (no cumulative integration)
        wsec = float(encoder_vavg_sec)
        if wsec > 0:
            if id_valid:
                az_enc, el_enc = _smooth_encoder_interp_local_by_id(
                    t_spec,
                    id_vals,
                    az_enc_raw,
                    el_enc_raw,
                    win_sec=wsec,
                )
            else:
                print("[warn] --encoder-vavg-sec > 0 requested, but valid spectral id labels are unavailable; smoothing disabled for safety.")
                az_enc = az_enc_raw
                el_enc = el_enc_raw
        else:
            az_enc = az_enc_raw
            el_enc = el_enc_raw

        az_true, el_true = _apply_azel_correction(az_enc, el_enc, dlon, dlat, azel_apply_mode)


    # Sun ephemeris at t_spec
    az_sun, el_sun = sun_altaz_deg(t_spec, rawdata_path, planet=planet, weather_meteo=weather_outside_meteo)

    # Sun-centered offsets
    off_az = _az_wrap_diff(az_true, az_sun)
    daz = off_az * np.cos(np.deg2rad(el_sun))
    d_el = el_true - el_sun

    idx = pd.to_datetime(t_spec, unit="s")
    df = pd.DataFrame(
        {
            "tp": tp1.astype(float),
            "tp1": tp1.astype(float),
            "ta_star": ta_star.astype(float),
            "position": position_str,
            "obsmode": position_str,
            "tp_hot": tp_hot_ref.astype(float),
            "tp_off": tp_off_ref.astype(float),
            "daz": daz.astype(float),
            "d_el": d_el.astype(float),
            "az_true": az_true.astype(float),
            "el_true": el_true.astype(float),
            "az_sun": az_sun.astype(float),
            "el_sun": el_sun.astype(float),
            "t_unix": t_spec.astype(float),
            "id": id_vals,
            # diagnostics
            "az_enc": az_enc.astype(float),
            "el_enc": el_enc.astype(float),
            "az_enc_raw": az_enc_raw.astype(float),
            "el_enc_raw": el_enc_raw.astype(float),
            "az_alt": lon_alt.astype(float),
            "el_alt": lat_alt.astype(float),
            "dlon": dlon.astype(float),
            "dlat": dlat.astype(float),
        },
        index=idx,
    )
    df.index.name = "timestamp"

    # store Tamb used (scalar) for downstream reference
    try:
        df.attrs["Tamb_k"] = float(tamb_use)
        df.attrs["tamb_source"] = str(tamb_source)
        df.attrs["encoder_table_resolved"] = str(encoder_table_resolved)
        df.attrs["altaz_table_resolved"] = str(altaz_table_resolved)
        df.attrs["weather_inside_table_resolved"] = str(weather_inside_table_resolved)
        df.attrs["weather_outside_table_resolved"] = str(weather_outside_table_resolved)
        df.attrs["encoder_time_col"] = str(encoder_time_col)
        df.attrs["altaz_time_col"] = str(altaz_time_col)
        df.attrs["weather_inside_time_col"] = str(weather_inside_time_col)
        df.attrs["weather_outside_time_col"] = str(weather_outside_time_col)
        df.attrs["outside_meteo_source"] = str(weather_outside_meteo.get("source"))
        df.attrs["outside_default_used_temperature"] = bool(weather_outside_meteo.get("used_default_temperature"))
        df.attrs["outside_default_used_pressure"] = bool(weather_outside_meteo.get("used_default_pressure"))
        df.attrs["outside_default_used_humidity"] = bool(weather_outside_meteo.get("used_default_humidity"))
        df.attrs["azel_correction_apply"] = str(azel_apply_mode)
    except Exception:
        pass
    try:
        df.attrs["spec_time_basis"] = str(spec_time_meta.get("applied"))
        df.attrs["spec_time_suffix"] = spec_time_meta.get("suffix")
        df.attrs["spec_time_fallback_field"] = spec_time_meta.get("fallback_field")
        df.attrs["spec_time_example"] = spec_time_meta.get("first_timestamp_text")
        df.attrs["spectrometer_time_offset_sec"] = float(spectrometer_time_offset_sec)
        df.attrs["encoder_shift_sec"] = float(encoder_shift_sec)
        df.attrs["encoder_az_time_offset_sec"] = float(encoder_az_time_offset_sec)
        df.attrs["encoder_el_time_offset_sec"] = float(encoder_el_time_offset_sec)
    except Exception:
        pass

    az_scans, el_scans = build_scans_from_id(df, id_vals)
    print(
        f"[info] built df: N={len(df)} scans: az={len(az_scans)} el={len(el_scans)}  "
        f"azel_source={src} azel_correction_apply={azel_apply_mode} encoder_vavg_sec={encoder_vavg_sec} "
        f"chopper_wheel={chopper_wheel} Tamb={float(tamb_use):.3g}K tamb_source={tamb_source}"
    )
    print(
        f"[info] spectral time basis: applied={spec_time_meta.get('applied')} "
        f"suffix={spec_time_meta.get('suffix')} fallback={spec_time_meta.get('fallback_field')} "
        f"example={spec_time_meta.get('first_timestamp_text')}"
    )
    print(
        f"[info] time offsets: spectrometer={float(spectrometer_time_offset_sec):+.6f}s "
        f"encoder_common={float(encoder_shift_sec):+.6f}s "
        f"encoder_az={float(encoder_az_time_offset_sec):+.6f}s "
        f"encoder_el={float(encoder_el_time_offset_sec):+.6f}s"
    )
    print(
        f"[info] tables: encoder={encoder_table_resolved} altaz={altaz_table_resolved} "
        f"weather_inside={weather_inside_table_resolved} weather_outside={weather_outside_table_resolved} "
        f"outside_source={weather_outside_meteo.get('source')}"
    )
    return df, az_scans, el_scans



# -----------------------------
# Limb fitting (two Gaussians on dY/dX)
# -----------------------------
def _gauss(x: np.ndarray, A: float, mu: float, sigma: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return A * np.exp(-(x - mu) ** 2 / (2.0 * sigma ** 2))


def fit_two_sides_edge(
    x: np.ndarray,
    dydx: np.ndarray,
    *,
    direct: float,
    fit_win_deg: float = DEFAULT_EDGE_FIT_WIN_DEG,
    fit_threshold: float = DEFAULT_EDGE_FIT_THRESHOLD,
    hpbw_init_arcsec: float = DEFAULT_HPBW_INIT_ARCSEC,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Fit two Gaussian peaks on a derivative profile (solar limb edges).

    Parameters
    ----------
    x, dydx : 1D arrays
        x coordinate (deg) and derivative dY/dX.
    direct : float
        Scan direction sign. Typically +1 if x increases with time, -1 otherwise.
        This is used only to decide which peak to fit first (positive vs negative);
        the returned pair is always sorted by mu (left then right).
    fit_win_deg : float
        Half-width (deg) of the local fitting window around each peak.
    fit_threshold : float
        Use only points with |dY/dX| >= fit_threshold * |peak| inside the window.
    hpbw_init_arcsec : float
        Initial guess of HPBW [arcsec] used to set sigma_init = HPBW/2.3548.

    Returns
    -------
    c_left, c_right : (A, mu, sigma) or None
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(dydx, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 20:
        return None, None

    # Trim ends (often include turn-around artifacts)
    trim = builtins.max(5, int(x.size * 0.05))
    if (x.size - 2 * trim) < 10:
        x_core, y_core = x, y
    else:
        x_core, y_core = x[trim:-trim], y[trim:-trim]

    # Light robust smoothing for peak detection only (median filter)
    try:
        y_smooth = signal.medfilt(y_core, kernel_size=5) if y_core.size >= 7 else y_core
    except Exception:
        y_smooth = y_core

    sigma_init = (float(hpbw_init_arcsec) / 3600.0) / 2.3548  # deg

    def _fit_peak(expected_sign: float) -> Optional[np.ndarray]:
        if y_smooth.size < 5:
            return None

        idx0 = int(np.argmax(y_smooth)) if expected_sign > 0 else int(np.argmin(y_smooth))
        peak_pos = float(x_core[idx0])
        peak_val_s = float(y_smooth[idx0])
        peak_val_r = float(y_core[idx0])

        win = float(fit_win_deg)
        mask_win = (x >= (peak_pos - win)) & (x <= (peak_pos + win))
        x_win = x[mask_win]
        y_win = y[mask_win]
        if x_win.size < 7:
            return None

        thr = float(fit_threshold)
        mask_thr = np.abs(y_win) >= (np.abs(peak_val_s) * thr)
        x_fit = x_win[mask_thr]
        y_fit = y_win[mask_thr]
        if x_fit.size < 7:
            # fallback: slightly looser threshold, then full window
            mask_thr = np.abs(y_win) >= (np.abs(peak_val_s) * 0.20)
            x_fit = x_win[mask_thr]
            y_fit = y_win[mask_thr]
            if x_fit.size < 7:
                x_fit, y_fit = x_win, y_win

        # Initial params
        p0 = [peak_val_r, peak_pos, sigma_init]

        # Bounds (keep it fairly loose but finite)
        if expected_sign > 0:
            A_min, A_max = peak_val_s * 0.2, peak_val_s * 3.0
        else:
            A_min, A_max = peak_val_s * 3.0, peak_val_s * 0.2
        A_min, A_max = (builtins.min(A_min, A_max), builtins.max(A_min, A_max))
        mu_min, mu_max = peak_pos - 0.10, peak_pos + 0.10
        sig_min, sig_max = sigma_init * 0.20, sigma_init * 4.00
        bounds = ([A_min, mu_min, sig_min], [A_max, mu_max, sig_max])

        try:
            # SciPy versions differ: curve_fit may not accept robust loss parameters.
            popt, _ = scipy.optimize.curve_fit(_gauss, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=20000)
            return popt
        except Exception:
            try:
                popt, _ = scipy.optimize.curve_fit(_gauss, x_fit, y_fit, p0=p0, maxfev=20000)
                return popt
            except Exception:
                return None

    c1 = _fit_peak(float(direct))
    c2 = _fit_peak(float(-direct))

    if (c1 is not None) and (c2 is not None) and (c1[1] > c2[1]):
        c1, c2 = c2, c1
    return c1, c2


def estimate_scan_speed_deg_s(seg: pd.DataFrame, xcol: str) -> float:
    """Estimate scan speed (deg/s) using the near-disk region around x≈0.

    This mimics the 'practical' method: fit a line x(t) in the region where
    the telescope is approximately in constant speed and crossing the solar disk.

    Notes
    -----
    - For daz (which already includes cos(el_ref)), this gives the *projected* sky speed.
    - If too few points are available in the default window, the window is widened.
    """
    if seg is None or len(seg) < 5 or (xcol not in seg.columns):
        return float("nan")

    if "t_unix" in seg.columns:
        t = pd.to_numeric(seg["t_unix"], errors="coerce").to_numpy(float)
    else:
        t = (seg.index.values.view("int64") / 1e9).astype(float)
    x = pd.to_numeric(seg[xcol], errors="coerce").to_numpy(float)

    m = np.isfinite(t) & np.isfinite(x)
    t, x = t[m], x[m]
    if t.size < 5:
        return float("nan")

    center_mask = (x >= -0.25) & (x <= 0.25)
    if np.count_nonzero(center_mask) < 10:
        center_mask = (x >= -0.40) & (x <= 0.40)
    if np.count_nonzero(center_mask) < 10:
        t_center, x_center = t, x
    else:
        t_center, x_center = t[center_mask], x[center_mask]

    try:
        p = np.polyfit(t_center - t_center[0], x_center, 1)
        return abs(float(p[0]))
    except Exception:
        return float("nan")


def save_text_summary(text: str, out_png: pathlib.Path) -> None:
    lines = text.split("\n")
    n_lines = len(lines)
    fig_h = builtins.max(4.0, n_lines * 0.25 + 1.0)
    fig_w = 12.0
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.text(
        0.03, 0.98, text,
        transform=ax.transAxes,
        fontsize=11,
        family="monospace",
        verticalalignment="top",
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.2, dpi=150)
    plt.close(fig)


def _plot_profile_and_derivative_panel(
    ax_prof,
    *,
    sid: int,
    axis_name: str,
    x_label: str,
    y_label: str,
    xlim_deg: Optional[float],
    raw_profile: Optional[Tuple[np.ndarray, np.ndarray]],
    corr_profile: Optional[Tuple[np.ndarray, np.ndarray]],
    raw_deriv: Optional[Tuple[np.ndarray, np.ndarray]],
    corr_deriv: Optional[Tuple[np.ndarray, np.ndarray]],
    fit_pair: Optional[Tuple[Optional[np.ndarray], Optional[np.ndarray]]],
    fit_labels: Tuple[str, str],
    result_row: Optional[dict],
) -> None:
    ax_der = ax_prof.twinx()

    handles = []
    labels = []

    if raw_profile is not None:
        xr, yr = raw_profile
        h, = ax_prof.plot(xr, yr, lw=0.8, alpha=DEFAULT_OVERLAY_RAW_ALPHA, label=f"{y_label} raw")
        handles.append(h); labels.append(h.get_label())
    if corr_profile is not None:
        xc, yc = corr_profile
        h, = ax_prof.plot(xc, yc, lw=1.2, label=f"{y_label} corr" if raw_profile is not None else y_label)
        handles.append(h); labels.append(h.get_label())

    if raw_deriv is not None:
        xr, dr = raw_deriv
        h, = ax_der.plot(xr, dr, lw=0.8, alpha=DEFAULT_OVERLAY_RAW_ALPHA, label="deriv raw")
        handles.append(h); labels.append(h.get_label())
    if corr_deriv is not None:
        xd, dd = corr_deriv
        h, = ax_der.plot(
            xd, dd,
            lw=1.4,
            color="deeppink",
            ls="-",
            label="deriv corr" if raw_deriv is not None else "derivative",
        )
        handles.append(h); labels.append(h.get_label())
        c1, c2 = fit_pair if fit_pair is not None else (None, None)
        if c1 is not None:
            h, = ax_der.plot(xd, _gauss(xd, *c1), ls="--", lw=1.0, label=fit_labels[0])
            handles.append(h); labels.append(h.get_label())
        if c2 is not None:
            h, = ax_der.plot(xd, _gauss(xd, *c2), ls="--", lw=1.0, label=fit_labels[1])
            handles.append(h); labels.append(h.get_label())

    ax_prof.grid(True)
    ax_prof.set_xlabel(x_label)
    ax_prof.set_ylabel(y_label)
    ax_der.set_ylabel(f"d({y_label})/d(x)")
    if xlim_deg is not None and float(xlim_deg) > 0:
        xl = float(xlim_deg)
        ax_prof.set_xlim(-xl, xl)
        ax_der.set_xlim(-xl, xl)

    title = f"{axis_name} scan {sid}"
    if result_row is not None:
        axis_key = axis_name.lower()
        title += (
            f"  repAz={result_row.get('rep_az_deg', float('nan')):.3f}"
            f" repEl={result_row.get('rep_el_deg', float('nan')):.3f}"
            f"  center={result_row.get('center_' + axis_key + '_deg', float('nan')) * 3600.0:.1f} arcsec"
            f"  HPBW={result_row.get('hpbw_' + axis_key + '_arcsec', float('nan')):.1f}"
            f"  v={result_row.get('speed_' + axis_key + '_arcsec_s', float('nan')):.0f}"
        )
    ax_prof.set_title(title)

    if handles:
        uniq_h = []
        uniq_l = []
        for h, l in zip(handles, labels):
            if l not in uniq_l:
                uniq_h.append(h)
                uniq_l.append(l)
        ax_prof.legend(uniq_h, uniq_l, fontsize=8, loc="best")


def plot_derivative_fits_paginated(
    outdir: pathlib.Path,
    *,
    tag: str,
    ycol: str,
    scan_ids: List[int],
    az_profile_raw: Dict[int, Tuple[np.ndarray, np.ndarray]],
    el_profile_raw: Dict[int, Tuple[np.ndarray, np.ndarray]],
    az_profile: Dict[int, Tuple[np.ndarray, np.ndarray]],
    el_profile: Dict[int, Tuple[np.ndarray, np.ndarray]],
    az_deriv_raw: Optional[Dict[int, Tuple[np.ndarray, np.ndarray]]] = None,
    el_deriv_raw: Optional[Dict[int, Tuple[np.ndarray, np.ndarray]]] = None,
    az_deriv: Optional[Dict[int, Tuple[np.ndarray, np.ndarray]]] = None,
    el_deriv: Optional[Dict[int, Tuple[np.ndarray, np.ndarray]]] = None,
    az_fit: Optional[Dict[int, Tuple[Optional[np.ndarray], Optional[np.ndarray]]]] = None,
    el_fit: Optional[Dict[int, Tuple[Optional[np.ndarray], Optional[np.ndarray]]]] = None,
    results: Optional[Dict[int, dict]] = None,
    xlim_deg: Optional[float] = None,
    scans_per_page: int = DEFAULT_EDGE_FIT_PLOT_MAX_SCANS,
) -> List[pathlib.Path]:
    scan_ids = sorted([int(s) for s in scan_ids])
    if len(scan_ids) == 0:
        return []

    y_label = _pretty_ylabel(ycol)
    scans_per_page = builtins.max(1, int(scans_per_page))
    out_paths: List[pathlib.Path] = []

    n_pages = int(math.ceil(len(scan_ids) / float(scans_per_page)))
    for ipage in range(n_pages):
        sub = scan_ids[ipage * scans_per_page:(ipage + 1) * scans_per_page]
        fig, axes = plt.subplots(len(sub), 2, figsize=(18, 4.8 * len(sub)), constrained_layout=True)
        if len(sub) == 1:
            axes = np.array([axes])
        if axes.ndim == 1:
            axes = axes.reshape(1, 2)

        for r, sid in enumerate(sub):
            rr = results.get(sid, {}) if isinstance(results, dict) else {}
            _plot_profile_and_derivative_panel(
                axes[r, 0],
                sid=sid,
                axis_name="AZ",
                x_label="daz [deg]",
                y_label=y_label,
                xlim_deg=xlim_deg,
                raw_profile=az_profile_raw.get(sid),
                corr_profile=az_profile.get(sid),
                raw_deriv=(az_deriv_raw or {}).get(sid),
                corr_deriv=(az_deriv or {}).get(sid),
                fit_pair=(az_fit or {}).get(sid),
                fit_labels=("fit left", "fit right"),
                result_row=rr,
            )
            _plot_profile_and_derivative_panel(
                axes[r, 1],
                sid=sid,
                axis_name="EL",
                x_label="d_el [deg]",
                y_label=y_label,
                xlim_deg=xlim_deg,
                raw_profile=el_profile_raw.get(sid),
                corr_profile=el_profile.get(sid),
                raw_deriv=(el_deriv_raw or {}).get(sid),
                corr_deriv=(el_deriv or {}).get(sid),
                fit_pair=(el_fit or {}).get(sid),
                fit_labels=("fit low", "fit high"),
                result_row=rr,
            )

        fig.suptitle(f"Sun scan profile + derivative summary: {tag}  page {ipage + 1}/{n_pages}", fontsize=14)
        out_png = outdir / f"sun_scan_derivative_fits_{tag}_p{ipage + 1:03d}.png"
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        out_paths.append(out_png)

    return out_paths


def write_scan_summary_csv(out_csv: pathlib.Path, scan_ids: List[int], results: Dict[int, dict]) -> None:
    rows = []
    for sid in sorted([int(s) for s in scan_ids]):
        rr = dict(results.get(sid, {}))
        rr = {k: v for k, v in rr.items() if not str(k).startswith("__")}
        rr["scan_id"] = int(sid)
        rows.append(rr)
    df_out = pd.DataFrame(rows)
    preferred = [
        "scan_id",
        "rep_az_deg", "rep_el_deg",
        "center_az_deg", "center_el_deg",
        "hpbw_az_arcsec", "hpbw_el_arcsec",
        "sun_az_deg", "sun_el_deg",
        "speed_az_arcsec_s", "speed_el_arcsec_s",
        "fit_ok_az", "fit_ok_el",
        "az_error", "el_error",
        "n_az", "n_az_used", "n_el", "n_el_used",
        "az_left_amp", "az_left_center_deg", "az_left_sigma_deg",
        "az_right_amp", "az_right_center_deg", "az_right_sigma_deg",
        "el_low_amp", "el_low_center_deg", "el_low_sigma_deg",
        "el_high_amp", "el_high_center_deg", "el_high_sigma_deg",
        "ripple_applied_az", "ripple_applied_el",
        "az_track_t_unix", "el_track_t_unix",
        "az_track_az_deg", "az_track_el_deg", "az_track_main_offset_deg", "az_track_cross_offset_deg",
        "el_track_az_deg", "el_track_el_deg", "el_track_main_offset_deg", "el_track_cross_offset_deg",
        "data_tag", "y_axis",
        "spec_time_basis", "spec_time_suffix", "spec_time_fallback_field", "spec_time_example",
        "azel_source", "azel_correction_apply", "altaz_apply",
        "spectrometer_time_offset_sec",
        "encoder_shift_sec", "encoder_az_time_offset_sec", "encoder_el_time_offset_sec", "encoder_vavg_sec",
        "chopper_wheel",
        "ripple_remove", "ripple_preset",
        "edge_fit_win_deg", "edge_fit_threshold", "hpbw_init_arcsec",
        "trim_scan", "profile_xlim_deg",
    ]
    cols = [c for c in preferred if c in df_out.columns] + [c for c in df_out.columns if c not in preferred]
    if len(cols) > 0:
        df_out = df_out.loc[:, cols]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_csv, index=False)


def run_edge_fits_and_summarize(
    *,
    az_scans: Dict[int, pd.DataFrame],
    el_scans: Dict[int, pd.DataFrame],
    ycol: str,
    trim_enabled: bool,
    trim_params: dict,
    strict_deriv: bool,
    ripple_enabled: bool,
    ripple_policy: Optional[dict],
    profile_xlim_deg: float,
    fit_win_deg: float,
    fit_threshold: float,
    hpbw_init_arcsec: float,
) -> Tuple[Dict[int, Tuple[np.ndarray, np.ndarray]], Dict[int, Tuple[np.ndarray, np.ndarray]],
           Dict[int, Tuple[Optional[np.ndarray], Optional[np.ndarray]]], Dict[int, Tuple[Optional[np.ndarray], Optional[np.ndarray]]],
           Dict[int, dict]]:
    """Compute derivative profiles, fit limb Gaussians, and summarize per scan id."""
    az_deriv: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    el_deriv: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    az_deriv_raw: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    el_deriv_raw: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    az_profile: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    el_profile: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    az_profile_raw: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    el_profile_raw: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    az_fit: Dict[int, Tuple[Optional[np.ndarray], Optional[np.ndarray]]] = {}
    el_fit: Dict[int, Tuple[Optional[np.ndarray], Optional[np.ndarray]]] = {}
    results: Dict[int, dict] = {}

    scan_ids = sorted(set(az_scans.keys()) | set(el_scans.keys()))
    for sid in scan_ids:
        rr = {"scan_id": int(sid)}

        # --- AZ ---
        if sid in az_scans:
            seg = az_scans[sid]
            rr["n_az"] = _safe_int(len(seg), 0)
            try:
                seg2 = trim_scan_segment(seg, "daz", "d_el", enabled=trim_enabled, trim_params=trim_params) if trim_enabled else seg
                t = seg2["t_unix"].to_numpy(float)
                x = seg2["daz"].to_numpy(float)
                y = seg2[ycol].to_numpy(float)
                rr["n_az_used"] = _safe_int(len(seg2), 0)
                track = _pick_tracking_point(seg, "daz", trim_enabled=trim_enabled, trim_params=trim_params)
                rr["az_track_t_unix"] = track["t_unix"]
                rr["az_track_az_deg"] = track["az_true_deg"]
                rr["az_track_el_deg"] = track["el_true_deg"]
                rr["az_track_sun_az_deg"] = track["az_sun_deg"]
                rr["az_track_sun_el_deg"] = track["el_sun_deg"]
                rr["az_track_main_offset_deg"] = track["main_offset_deg"]
                rr["az_track_cross_offset_deg"] = track["cross_offset_deg"]
                y_used = y
                applied = False
                if x.size > 0 and y.size > 0:
                    az_profile_raw[sid] = (x.copy(), y.copy())
                if ripple_enabled:
                    try:
                        cfg = resolve_ripple_cfg_for_scan(ripple_policy, t=t, x_main_deg=x, profile_xlim_deg=float(profile_xlim_deg))
                        y_corr, _trend, _ripple = remove_ripple_notch(t, y, **cfg)
                        if (y_corr is not None) and (np.asarray(y_corr).shape == np.asarray(y).shape) and np.all(np.isfinite(y_corr)):
                            delta = float(np.nanmax(np.abs(y_corr - y))) if np.any(np.isfinite(y_corr - y)) else 0.0
                            scale = float(np.nanmax(np.abs(y))) if np.any(np.isfinite(y)) else 1.0
                            applied = (delta > (1e-12 * scale + 1e-6))
                            y_used = y_corr if applied else y
                    except Exception:
                        y_used = y
                        applied = False
                if x.size > 0 and np.asarray(y_used).size > 0:
                    az_profile[sid] = (x.copy(), np.asarray(y_used, dtype=float).copy())
                direct = 1.0 if np.nanmedian(np.diff(x[np.isfinite(x)])) >= 0 else -1.0
                rr["ripple_applied_az"] = bool(applied)
                x_d_raw, d_raw = finite_difference_skip_duplicates(x, y, strict=bool(strict_deriv))
                az_deriv_raw[sid] = (x_d_raw, d_raw)
                x_d, d = finite_difference_skip_duplicates(x, y_used, strict=bool(strict_deriv))
                az_deriv[sid] = (x_d, d)
                c1, c2 = fit_two_sides_edge(
                    x_d, d, direct=direct,
                    fit_win_deg=fit_win_deg, fit_threshold=fit_threshold, hpbw_init_arcsec=hpbw_init_arcsec
                )
                az_fit[sid] = (c1, c2)
                rr["fit_ok_az"] = bool((c1 is not None) and (c2 is not None))
                if c1 is not None:
                    rr["az_left_amp"] = float(c1[0]); rr["az_left_center_deg"] = float(c1[1]); rr["az_left_sigma_deg"] = float(c1[2])
                if c2 is not None:
                    rr["az_right_amp"] = float(c2[0]); rr["az_right_center_deg"] = float(c2[1]); rr["az_right_sigma_deg"] = float(c2[2])
                if (c1 is not None) and (c2 is not None):
                    center = 0.5 * (float(c1[1]) + float(c2[1]))
                    sun = float(c2[1]) - float(c1[1])
                    sigma = 0.5 * (abs(float(c1[2])) + abs(float(c2[2])))
                    hpbw = sigma * 2.3548 * 3600.0
                    rr["center_az_deg"] = center
                    rr["sun_az_deg"] = sun
                    rr["sigma_az_deg"] = sigma
                    rr["hpbw_az_arcsec"] = hpbw
                rr["speed_az_arcsec_s"] = estimate_scan_speed_deg_s(seg2, "daz") * 3600.0
            except Exception as e:
                rr.setdefault("fit_ok_az", False)
                rr["az_error"] = str(e)

        # --- EL ---
        if sid in el_scans:
            seg = el_scans[sid]
            rr["n_el"] = _safe_int(len(seg), 0)
            try:
                seg2 = trim_scan_segment(seg, "d_el", "daz", enabled=trim_enabled, trim_params=trim_params) if trim_enabled else seg
                t = seg2["t_unix"].to_numpy(float)
                x = seg2["d_el"].to_numpy(float)
                y = seg2[ycol].to_numpy(float)
                rr["n_el_used"] = _safe_int(len(seg2), 0)
                track = _pick_tracking_point(seg, "d_el", trim_enabled=trim_enabled, trim_params=trim_params)
                rr["el_track_t_unix"] = track["t_unix"]
                rr["el_track_az_deg"] = track["az_true_deg"]
                rr["el_track_el_deg"] = track["el_true_deg"]
                rr["el_track_sun_az_deg"] = track["az_sun_deg"]
                rr["el_track_sun_el_deg"] = track["el_sun_deg"]
                rr["el_track_main_offset_deg"] = track["main_offset_deg"]
                rr["el_track_cross_offset_deg"] = track["cross_offset_deg"]
                y_used = y
                applied = False
                if x.size > 0 and y.size > 0:
                    el_profile_raw[sid] = (x.copy(), y.copy())
                if ripple_enabled:
                    try:
                        cfg = resolve_ripple_cfg_for_scan(ripple_policy, t=t, x_main_deg=x, profile_xlim_deg=float(profile_xlim_deg))
                        y_corr, _trend, _ripple = remove_ripple_notch(t, y, **cfg)
                        if (y_corr is not None) and (np.asarray(y_corr).shape == np.asarray(y).shape) and np.all(np.isfinite(y_corr)):
                            delta = float(np.nanmax(np.abs(y_corr - y))) if np.any(np.isfinite(y_corr - y)) else 0.0
                            scale = float(np.nanmax(np.abs(y))) if np.any(np.isfinite(y)) else 1.0
                            applied = (delta > (1e-12 * scale + 1e-6))
                            y_used = y_corr if applied else y
                    except Exception:
                        y_used = y
                        applied = False
                if x.size > 0 and np.asarray(y_used).size > 0:
                    el_profile[sid] = (x.copy(), np.asarray(y_used, dtype=float).copy())
                direct = 1.0 if np.nanmedian(np.diff(x[np.isfinite(x)])) >= 0 else -1.0
                rr["ripple_applied_el"] = bool(applied)
                x_d_raw, d_raw = finite_difference_skip_duplicates(x, y, strict=bool(strict_deriv))
                el_deriv_raw[sid] = (x_d_raw, d_raw)
                x_d, d = finite_difference_skip_duplicates(x, y_used, strict=bool(strict_deriv))
                el_deriv[sid] = (x_d, d)
                c1, c2 = fit_two_sides_edge(
                    x_d, d, direct=direct,
                    fit_win_deg=fit_win_deg, fit_threshold=fit_threshold, hpbw_init_arcsec=hpbw_init_arcsec
                )
                el_fit[sid] = (c1, c2)
                rr["fit_ok_el"] = bool((c1 is not None) and (c2 is not None))
                if c1 is not None:
                    rr["el_low_amp"] = float(c1[0]); rr["el_low_center_deg"] = float(c1[1]); rr["el_low_sigma_deg"] = float(c1[2])
                if c2 is not None:
                    rr["el_high_amp"] = float(c2[0]); rr["el_high_center_deg"] = float(c2[1]); rr["el_high_sigma_deg"] = float(c2[2])
                if (c1 is not None) and (c2 is not None):
                    center = 0.5 * (float(c1[1]) + float(c2[1]))
                    sun = float(c2[1]) - float(c1[1])
                    sigma = 0.5 * (abs(float(c1[2])) + abs(float(c2[2])))
                    hpbw = sigma * 2.3548 * 3600.0
                    rr["center_el_deg"] = center
                    rr["sun_el_deg"] = sun
                    rr["sigma_el_deg"] = sigma
                    rr["hpbw_el_arcsec"] = hpbw
                rr["speed_el_arcsec_s"] = estimate_scan_speed_deg_s(seg2, "d_el") * 3600.0
            except Exception as e:
                rr.setdefault("fit_ok_el", False)
                rr["el_error"] = str(e)

        rr["rep_az_deg"] = _nanmean_pair(rr.get("az_track_az_deg", np.nan), rr.get("el_track_az_deg", np.nan))
        rr["rep_el_deg"] = _nanmean_pair(rr.get("az_track_el_deg", np.nan), rr.get("el_track_el_deg", np.nan))
        results[sid] = rr

    results["__scan_ids__"] = scan_ids
    results["__az_deriv_raw__"] = az_deriv_raw
    results["__el_deriv_raw__"] = el_deriv_raw
    results["__az_profile_raw__"] = az_profile_raw
    results["__el_profile_raw__"] = el_profile_raw
    results["__az_profile__"] = az_profile
    results["__el_profile__"] = el_profile

    return az_deriv, el_deriv, az_fit, el_fit, results

# -----------------------------
# Plotting
# -----------------------------
def _pretty_ylabel(ycol: str) -> str:
    """Human-friendly y-axis label."""
    ycol = str(ycol)
    if ycol == "ta_star":
        return "Ta* [K]"
    if ycol == "tp1":
        return "tp1 [arb]"
    return ycol

def plot_dxdy_simple(
    seg: pd.DataFrame,
    title: str,
    out_png: pathlib.Path,
    *,
    trim_enabled: bool,
    trim_params: dict,
    which_xcol: str,
) -> None:
    # optionally add trimmed overlay
    t = seg["t_unix"].to_numpy(float)
    x = seg["daz"].to_numpy(float)
    y = seg["d_el"].to_numpy(float)

    dt = np.diff(t)
    dx = np.diff(x)
    dy = np.diff(y)
    t_mid = 0.5 * (t[1:] + t[:-1])

    bad_dt = (~np.isfinite(dt)) | (dt <= 0)
    vx = np.full_like(dx, np.nan, dtype=float)
    vy = np.full_like(dy, np.nan, dtype=float)
    ok_dt = ~bad_dt
    vx[ok_dt] = dx[ok_dt] / dt[ok_dt]
    vy[ok_dt] = dy[ok_dt] / dt[ok_dt]

    footer = (
        f"median|Δdaz|={np.median(np.abs(dx)):.3g} deg   median|Δd_el|={np.median(np.abs(dy)):.3g} deg   "
        f"bad Δt: {int(np.sum(bad_dt))}/{len(dt)}"
    )

    fig, ax = plt.subplots(3, 1, figsize=(12, 10), sharex=False)

    ax[0].plot(t, x, lw=1.0, label="daz [deg]")
    ax[0].plot(t, y, lw=1.0, label="d_el [deg]")

    if trim_enabled:
        try:
            segt = trim_scan_segment(seg, which_xcol, ("d_el" if which_xcol == "daz" else "daz"), enabled=True, trim_params=trim_params)
            tt = segt["t_unix"].to_numpy(float)
            xx = segt["daz"].to_numpy(float)
            yy = segt["d_el"].to_numpy(float)
            ax[0].plot(tt, xx, lw=2.0, ls="--", label="daz (trimmed)")
            ax[0].plot(tt, yy, lw=2.0, ls="--", label="d_el (trimmed)")
        except Exception:
            pass

    ax[0].set_title(title)
    ax[0].set_xlabel("t_unix [s]")
    ax[0].set_ylabel("Sun offset [deg]")
    ax[0].grid(True)
    ax[0].legend()

    ax[1].plot(t_mid, dx, lw=1.0, label="Δdaz per sample [deg]")
    ax[1].plot(t_mid, dy, lw=1.0, label="Δd_el per sample [deg]")
    ax[1].set_xlabel("t_mid [s]")
    ax[1].set_ylabel("Δoffset [deg]")
    ax[1].grid(True)
    ax[1].legend()

    ax[2].plot(t_mid, vx, lw=1.0, label="Δdaz/Δt [deg/s]")
    ax[2].plot(t_mid, vy, lw=1.0, label="Δd_el/Δt [deg/s]")
    ax[2].set_xlabel("t_mid [s]")
    ax[2].set_ylabel("speed [deg/s]")
    ax[2].grid(True)
    ax[2].legend()

    fig.text(0.01, 0.01, footer, fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_mount_azel_simple(seg: pd.DataFrame, title: str, out_png: pathlib.Path) -> None:
    t = seg["t_unix"].to_numpy(float)
    az = seg["az_true"].to_numpy(float)
    el = seg["el_true"].to_numpy(float)

    dt = np.diff(t)
    daz = np.diff(az)
    delv = np.diff(el)
    t_mid = 0.5 * (t[1:] + t[:-1])

    bad_dt = (~np.isfinite(dt)) | (dt <= 0)
    vaz = np.full_like(daz, np.nan, dtype=float)
    vel = np.full_like(delv, np.nan, dtype=float)
    ok = ~bad_dt
    vaz[ok] = daz[ok] / dt[ok]
    vel[ok] = delv[ok] / dt[ok]

    footer = (
        f"median|ΔAz_true|={np.median(np.abs(daz)):.3g} deg   median|ΔEl_true|={np.median(np.abs(delv)):.3g} deg   "
        f"bad Δt: {int(np.sum(bad_dt))}/{len(dt)}"
    )

    fig, ax = plt.subplots(3, 1, figsize=(12, 10), sharex=False)

    ax[0].plot(t, az, lw=1.0, label="Az_true [deg]")
    ax[0].plot(t, el, lw=1.0, label="El_true [deg]")
    ax[0].set_title(title)
    ax[0].set_xlabel("t_unix [s]")
    ax[0].set_ylabel("mount coord [deg]")
    ax[0].grid(True)
    ax[0].legend()

    ax[1].plot(t_mid, daz, lw=1.0, label="ΔAz_true per sample [deg]")
    ax[1].plot(t_mid, delv, lw=1.0, label="ΔEl_true per sample [deg]")
    ax[1].set_xlabel("t_mid [s]")
    ax[1].set_ylabel("Δmount [deg]")
    ax[1].grid(True)
    ax[1].legend()

    ax[2].plot(t_mid, vaz, lw=1.0, label="ΔAz_true/Δt [deg/s]")
    ax[2].plot(t_mid, vel, lw=1.0, label="ΔEl_true/Δt [deg/s]")
    ax[2].set_xlabel("t_mid [s]")
    ax[2].set_ylabel("mount speed [deg/s]")
    ax[2].grid(True)
    ax[2].legend()

    fig.text(0.01, 0.01, footer, fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)



def plot_profile_and_deriv(
    seg: pd.DataFrame,
    xcol: str,
    ycol: str,
    title: str,
    out_png: pathlib.Path,
    *,
    strict_deriv: bool,
    trim_enabled: bool,
    trim_params: dict,
    xlim_deg: Optional[float] = None,
    ripple_enabled: bool = False,
    ripple_policy: Optional[dict] = None,
) -> None:
    seg2 = trim_scan_segment(seg, xcol, ("d_el" if xcol == "daz" else "daz"), enabled=trim_enabled, trim_params=trim_params) if trim_enabled else seg

    x = seg2[xcol].to_numpy(float)
    y = seg2[ycol].to_numpy(float)
    y_used = y
    y_base_label = _pretty_ylabel(ycol)
    y_label = f"{y_base_label} (raw)"
    note = ""
    applied = False

    if ripple_enabled:
        try:
            if "t_unix" not in seg2.columns:
                raise RuntimeError("ripple removal requires 't_unix' in segment")
            t = seg2["t_unix"].to_numpy(float)

            cfg = resolve_ripple_cfg_for_scan(
                ripple_policy,
                t=t,
                x_main_deg=x,
                profile_xlim_deg=(float(xlim_deg) if (xlim_deg is not None and float(xlim_deg) > 0) else float(DEFAULT_PROFILE_XLIM_DEG)),
            )
            y_corr, s_tr, rip = remove_ripple_notch(t, y, **cfg)

            if np.all(np.isfinite(y_corr)) and y_corr.size == y.size:
                # sanity: reject obviously unstable results (e.g., multiplicative blow-up)
                max_raw = float(np.nanmax(np.abs(y))) if np.any(np.isfinite(y)) else 0.0
                max_corr = float(np.nanmax(np.abs(y_corr))) if np.any(np.isfinite(y_corr)) else np.inf
                if max_raw > 0 and (not np.isfinite(max_corr) or max_corr > 50.0 * max_raw):
                    raise RuntimeError(
                        f"unstable ripple result: max|corr|={max_corr:.3g} > 50*max|raw|={max_raw:.3g}"
                    )
                                # Apply only if correction actually changes the series (auto may decide 'none').
                delta = float(np.nanmax(np.abs(y_corr - y))) if np.any(np.isfinite(y_corr - y)) else 0.0
                scale = float(np.nanmax(np.abs(y))) if np.any(np.isfinite(y)) else 1.0
                applied = (delta > (1e-12 * scale + 1e-6))
                if applied:
                    y_used = y_corr
                    y_label = f"{y_base_label} (corr; model={cfg.get('model','auto')})"
                    note = " (ripple removed)"
        except Exception as e:
            # keep raw on failure
            note = f" (ripple removal failed: {e})"
    fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=False)

    if ripple_enabled and applied:
        # Profile overlay (raw vs corrected)
        ax[0].plot(x, y, lw=0.8, alpha=DEFAULT_OVERLAY_RAW_ALPHA, label=f"{y_base_label} (raw)")
        ax[0].plot(x, y_used, lw=1.2, label=y_label)
        ax[0].legend()

        # Derivative overlay (raw vs corrected)
        x_d_raw, dydx_raw = finite_difference_skip_duplicates(x, y, strict=bool(strict_deriv))
        x_d_corr, dydx_corr = finite_difference_skip_duplicates(x, y_used, strict=bool(strict_deriv))
        ax[1].plot(x_d_raw, dydx_raw, lw=0.8, alpha=DEFAULT_OVERLAY_RAW_ALPHA, label="raw")
        ax[1].plot(x_d_corr, dydx_corr, lw=1.2, label="corr")
        ax[1].legend()
    else:
        ax[0].plot(x, y_used, lw=1.0)
        x_d_corr, dydx_corr = finite_difference_skip_duplicates(x, y_used, strict=bool(strict_deriv))
        ax[1].plot(x_d_corr, dydx_corr, lw=1.0)

    ax[0].set_title(title + note)
    ax[0].set_xlabel(f"{xcol} [deg]")
    ax[0].set_ylabel(y_label)
    ax[0].grid(True)

    # x-range (offset axis)
    if xlim_deg is not None and float(xlim_deg) > 0:
        xl = float(xlim_deg)
        ax[0].set_xlim(-xl, xl)
        ax[1].set_xlim(-xl, xl)

    ax[1].set_xlabel(f"{xcol} [deg]")
    ax[1].set_ylabel(f"d({ycol})/d({xcol}) (finite diff; skip duplicates)")
    ax[1].grid(True)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("rawdata", help="RawData directory")
    ap.add_argument("--spectral", default=DEFAULT_SPECTRAL_NAME, help=f"Spectral name (default: {DEFAULT_SPECTRAL_NAME})")
    ap.add_argument("--db-namespace", default=DEFAULT_DB_NAMESPACE, help=f"Database namespace/prefix for NECST tables (default: {DEFAULT_DB_NAMESPACE})")
    ap.add_argument("--telescope", default=TELESCOPE, help=f"Telescope name used in NECST DB table names (default: {TELESCOPE})")
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR, help="Output directory for PNGs (default: current dir)")

    ap.add_argument("--debug-plot", action="store_true", default=DEFAULT_DEBUG_PLOT,
                    help="Enable debug plots (dxdy_*, mount_azel_*, summary_text_*).")
    ap.add_argument("--profile-xlim-deg", type=float, default=DEFAULT_PROFILE_XLIM_DEG,
                    help="Profile/derivative x-axis range in degrees: plot [-xlim, +xlim].")


    ap.add_argument("--azel-source", choices=["encoder", "altaz"], default=DEFAULT_AZEL_SOURCE,
                    help="Az/El source. encoder uses encoder lon/lat; altaz uses altaz lon/lat (commands).")
    ap.add_argument("--azel-correction-apply", choices=["none", "subtract", "add"], default=None,
                    help="Apply dlon/dlat to the selected Az/El source. Works for both --azel-source encoder and altaz.")
    ap.add_argument("--altaz-apply", choices=["none", "minus", "plus"], default=None,
                    help="Legacy alias for --azel-correction-apply.")
    ap.add_argument("--encoder-table", default=None, help="Full encoder table name. If omitted, resolve from db_namespace/telescope and suffix.")
    ap.add_argument("--encoder-table-suffix", default="ctrl-antenna-encoder", help="Encoder table suffix used with db_namespace/telescope when --encoder-table is omitted.")
    ap.add_argument("--altaz-table", default=None, help="Full altaz/cmd table name. If omitted, resolve from db_namespace/telescope and suffix.")
    ap.add_argument("--altaz-table-suffix", default="ctrl-antenna-altaz", help="Altaz/cmd table suffix used with db_namespace/telescope when --altaz-table is omitted.")
    ap.add_argument("--encoder-time-col", default="time", help="Encoder time column (recommended: time).")
    ap.add_argument("--altaz-time-col", default="time", help="Altaz/cmd time column (recommended: time).")
    ap.add_argument("--spectrometer-time-offset-sec", type=float, default=0.0,
                    help="Time offset (sec) applied to spectral timestamps before interpolation. Corrected time = recorded spectral time + offset.")
    ap.add_argument("--encoder-shift-sec", type=float, default=DEFAULT_ENCODER_SHIFT_SEC,
                    help="Common time shift (sec) added to encoder.time before interpolation.")
    ap.add_argument("--encoder-az-time-offset-sec", type=float, default=0.0,
                    help="Additional time offset (sec) applied only to encoder Az interpolation.")
    ap.add_argument("--encoder-el-time-offset-sec", type=float, default=0.0,
                    help="Additional time offset (sec) applied only to encoder El interpolation.")

    ap.add_argument("--tel-loaddata", default=TEL_LOADDATA,
                    help="Deprecated compatibility option from the old nercst path. It is now ignored because spectral data are read directly from necstdb.")

    ap.add_argument("--encoder-vavg-sec", type=float, default=DEFAULT_ENCODER_VAVG_SEC,
                    help="Only for --azel-source encoder: time-window (sec) for safe local smoothing of interpolated encoder az/el within each contiguous id block. No differentiation/integration. 0 disables.")

    # Chopper-wheel Ta* calibration (1-temperature; Tamb assumed)
    ap.add_argument("--no-chopper-wheel", dest="chopper_wheel", action="store_false",
                    help="Disable 1-temp chopper-wheel calibration (Ta*). If disabled, use raw tp1 as the vertical axis.")
    ap.add_argument(
        "--tamb-k",
        type=float,
        default=None,
        help=(
            "Ambient load temperature Tamb [K]. If omitted, estimate from necstdb weather/spectral data "
            "(temperature/temp_k/Tamb) and use it if within configured limits; otherwise fallback to --tamb-default-k."
        ),
    )
    ap.add_argument("--tamb-default-k", type=float, default=300.0, help="Fallback Tamb [K] used when inside weather/spectral temperature is missing or invalid.")
    ap.add_argument("--tamb-min-k", type=float, default=250.0, help="Minimum allowed Tamb [K] before fallback is used.")
    ap.add_argument("--tamb-max-k", type=float, default=330.0, help="Maximum allowed Tamb [K] before fallback is used.")
    ap.add_argument("--weather-table", default=None, help="Legacy weather table alias applied to both inside/outside weather when specific tables are not set.")
    ap.add_argument("--weather-time-col", default=None, help="Legacy weather time column alias applied to both inside/outside weather time columns.")
    ap.add_argument("--weather-inside-table", default=None, help="Full inside-weather table name for Tamb estimation.")
    ap.add_argument("--weather-inside-table-suffix", default="weather-ambient", help="Inside-weather table suffix used with db_namespace/telescope when --weather-inside-table is omitted.")
    ap.add_argument("--weather-inside-time-col", default=None, help="Inside-weather time column.")
    ap.add_argument("--weather-outside-table", default=None, help="Full outside-weather table name for refraction meteorology.")
    ap.add_argument("--weather-outside-table-suffix", default="weather-ambient", help="Outside-weather table suffix used with db_namespace/telescope when --weather-outside-table is omitted.")
    ap.add_argument("--weather-outside-time-col", default=None, help="Outside-weather time column.")
    ap.add_argument("--outside-default-temperature-c", type=float, default=0.0, help="Fallback outside temperature [C] for refraction.")
    ap.add_argument("--outside-default-pressure-hpa", type=float, default=760.0, help="Fallback outside pressure [hPa] for refraction.")
    ap.add_argument("--outside-default-humidity-pct", type=float, default=30.0, help="Fallback outside relative humidity [%] for refraction.")
    ap.add_argument("--outside-temperature-min-c", type=float, default=-50.0, help="Minimum allowed outside temperature [C] before fallback is used.")
    ap.add_argument("--outside-temperature-max-c", type=float, default=50.0, help="Maximum allowed outside temperature [C] before fallback is used.")
    ap.add_argument("--outside-pressure-min-hpa", type=float, default=400.0, help="Minimum allowed outside pressure [hPa] before fallback is used.")
    ap.add_argument("--outside-pressure-max-hpa", type=float, default=1100.0, help="Maximum allowed outside pressure [hPa] before fallback is used.")
    ap.add_argument("--outside-humidity-min-pct", type=float, default=0.0, help="Minimum allowed outside humidity [%] before fallback is used.")
    ap.add_argument("--outside-humidity-max-pct", type=float, default=100.0, help="Maximum allowed outside humidity [%] before fallback is used.")
    # The following two options are kept only for backward compatibility; the current
    # 1-load chopper-wheel implementation uses direct np.interp(edge-hold) without smoothing.
    ap.add_argument("--chopper-win-sec", type=float, default=DEFAULT_CHOPPER_WIN_SEC,
                    help="(deprecated; ignored) Time window (sec) to smooth HOT/OFF reference powers before interpolation.")
    ap.add_argument("--chopper-stat", choices=["median", "mean"], default=DEFAULT_CHOPPER_STAT,
                    help="(deprecated; ignored) Statistic used to smooth HOT/OFF reference powers.")
    ap.set_defaults(chopper_wheel=DEFAULT_CHOPPER_WHEEL)


    # Ripple removal (periodic noise near target frequency and harmonics)
    # Enabled by default (DEFAULT_RIPPLE_REMOVE=True). Use --ripple-no-remove to disable.
    ap.add_argument(
        "--ripple-no-remove",
        dest="ripple_remove",
        action="store_false",
        help="Disable notch-based ripple removal on the vertical-axis signal within each scan segment (Ta* if chopper-wheel enabled; otherwise tp1).",
    )
    ap.set_defaults(ripple_remove=DEFAULT_RIPPLE_REMOVE)

    # Preset controls the baseline notch behavior. Any explicit --ripple-* option overrides it.
    ap.add_argument(
        "--ripple-preset",
        choices=["auto", "safe", "normal", "strong"],
        default=DEFAULT_RIPPLE_PRESET,
        help="Ripple-notch tuning preset. 'auto' selects a baseline per scan from scan speed; manual --ripple-* options override preset values.",
    )

    ap.add_argument("--ripple-model", choices=["auto", "add", "mul"], default=DEFAULT_RIPPLE_MODEL,
                    help="Residual model for ripple removal: auto may choose add/mul/none; add uses y-trend; mul uses y/trend-1.")
    ap.add_argument("--ripple-target-hz", type=float, default=DEFAULT_RIPPLE_TARGET_HZ, help="Target ripple frequency (Hz).")
    ap.add_argument("--ripple-search-hz", type=float, default=DEFAULT_RIPPLE_SEARCH_HZ, help="Peak search half-width around target (Hz).")

    # The following are overrides (default None => use preset values)
    ap.add_argument("--ripple-bw-hz", type=float, default=None,
                    help="Override notch bandwidth (Hz) around f_peak (and harmonics). If omitted, uses preset.")
    ap.add_argument("--ripple-max-harm", type=int, default=None,
                    help="Override max harmonic multiplier to notch (k=1..N). If omitted, uses preset (Nyquist-clipped).")
    ap.add_argument("--ripple-order", type=int, default=None,
                    help="Override Butterworth bandstop filter order. If omitted, uses preset.")
    ap.add_argument("--ripple-notch-pass", type=int, default=None,
                    help="Override number of repeated notch applications per harmonic. If omitted, uses preset.")
    ap.add_argument("--ripple-trend-win-sec", type=float, default=None,
                    help="Override trend (baseline) rolling-mean window (sec). If omitted, uses preset.")
    ap.add_argument("--ripple-resample-dt-sec", type=float, default=None,
                    help="Override uniform resampling dt (sec). <=0 means auto (median diff(t)). If omitted, uses preset.")
    ap.add_argument("--ripple-eval-band-hz", type=float, default=None,
                    help="Override band width (Hz) for FFT band-power evaluation. <=0 means use bw. If omitted, uses bw.")


    # Limb fitting (two Gaussians on derivative profiles) and summary outputs
    ap.add_argument("--no-edge-fit", dest="edge_fit", action="store_false",
                    help="Disable limb edge Gaussian fitting and summary outputs.")
    ap.add_argument("--edge-fit-win-deg", type=float, default=DEFAULT_EDGE_FIT_WIN_DEG,
                    help="Half-width (deg) of local window around each derivative peak for Gaussian fitting.")
    ap.add_argument("--edge-fit-threshold", type=float, default=DEFAULT_EDGE_FIT_THRESHOLD,
                    help="Keep points with |dY/dX| >= threshold*|peak| for Gaussian fitting.")
    ap.add_argument("--hpbw-init-arcsec", type=float, default=DEFAULT_HPBW_INIT_ARCSEC,
                    help="Initial HPBW guess [arcsec] used to set sigma_init for Gaussian fitting.")
    ap.add_argument("--edge-fit-plot-max-scans", type=int, default=DEFAULT_EDGE_FIT_PLOT_MAX_SCANS,
                    help="How many scan pairs to include per derivative-summary page. All scans are written across multiple pages if needed.")
    ap.set_defaults(edge_fit=DEFAULT_EDGE_FIT)

    ap.add_argument("--no-trim-scan", dest="trim_scan", action="store_false",
                    help="Disable motion-based trimming of scan (default: enabled).")
    ap.add_argument("--trim-vfrac", type=float, default=DEFAULT_TRIM_VFRAC, help="Keep vabs >= vfrac * P90(vabs).")
    ap.add_argument("--trim-vmin", type=float, default=DEFAULT_TRIM_VMIN, help="Absolute floor for trimming threshold [deg/s].")
    ap.add_argument("--trim-gap", type=int, default=DEFAULT_TRIM_GAP, help="Fill gaps up to this many samples in trim mask.")
    ap.add_argument("--trim-min-samples", type=int, default=DEFAULT_TRIM_MIN_SAMPLES, help="Minimum samples after trimming.")
    
    ap.add_argument(
        "--trim-dominant-axis",
        dest="trim_dominant_axis",
        action="store_true",
        help="Prefer dominant-axis trimming to exclude slews where both axes move (recommended).",
    )
    ap.add_argument(
        "--trim-no-dominant-axis",
        dest="trim_dominant_axis",
        action="store_false",
        help="Disable dominant-axis trimming and use legacy 1D trimming (|Δx/Δt| only).",
    )
    ap.add_argument("--trim-axis-ratio-min", type=float, default=DEFAULT_TRIM_AXIS_RATIO_MIN,
                    help="For dominant-axis trimming: require |v_main|/(|v_cross|+eps) >= this (default: %(default)s).")
    ap.add_argument("--trim-vpercentile", type=float, default=DEFAULT_TRIM_VPERCENTILE,
                    help="For dominant-axis trimming: percentile of |v_main| used for threshold (default: %(default)s).")
    # Steady-scan selection near Sun (ON-only; avoids slews / HOT-OFF even if id matches)
    ap.add_argument("--trim-no-steady-scan", dest="trim_steady_scan", action="store_false",
                    help="Disable steady-scan selection near Sun and use only the legacy trimming logic.")
    ap.add_argument("--trim-include-hotoff", dest="trim_use_on_only", action="store_false",
                    help="Include HOT/OFF samples when selecting scan segment (not recommended).")
    ap.add_argument("--trim-scan-speed-min-arcsec", type=float, default=DEFAULT_TRIM_SCAN_SPEED_MIN_ARCSEC_S,
                    help="Minimum |v_main| [arcsec/s] to consider a segment a scan (rejects near-stationary drift).")
    ap.add_argument("--trim-xwin-factor", type=float, default=DEFAULT_TRIM_XWIN_FACTOR,
                    help="Consider points within |offset_main| <= profile_xlim_deg * factor as near-Sun candidates.")
    ap.add_argument("--trim-cross-offset-max-deg", type=float, default=DEFAULT_TRIM_CROSS_OFFSET_MAX_DEG,
                    help="Allow cross-axis offset deviation from its median within this range [deg] when selecting scan.")
    ap.add_argument("--trim-steady-cv-max", type=float, default=DEFAULT_TRIM_STEADY_CV_MAX,
                    help="Max robust CV (sigma/|mean|) of v_main within the selected scan block.")
    ap.set_defaults(trim_steady_scan=DEFAULT_TRIM_STEADY_SCAN)
    ap.set_defaults(trim_use_on_only=DEFAULT_TRIM_USE_ON_ONLY)

    ap.set_defaults(trim_scan=DEFAULT_TRIM_SCAN)
    ap.set_defaults(trim_dominant_axis=DEFAULT_TRIM_DOMINANT_AXIS)

    ap.add_argument("--strict-deriv", action="store_true", help="If set, undefined derivative points raise error (default).")
    ap.add_argument("--no-strict-deriv", dest="strict_deriv", action="store_false", help="If set, undefined points become NaN.")
    ap.set_defaults(strict_deriv=DEFAULT_STRICT_DERIV)

    ap.add_argument("--continue-on-error", action="store_true", help="Continue processing other scans even if one scan errors.")
    ap.set_defaults(continue_on_error=DEFAULT_CONTINUE_ON_ERROR)
    args = ap.parse_args()

    rawdata_path = pathlib.Path(args.rawdata).expanduser().resolve()
    if not rawdata_path.exists():
        raise SystemExit(f"RawData path not found: {rawdata_path}")

    outdir = pathlib.Path(args.outdir).expanduser().resolve()
    tag = rawdata_path.name
    outdir.mkdir(parents=True, exist_ok=True)

    legacy_weather_table = args.weather_table
    legacy_weather_time_col = args.weather_time_col or "time"
    weather_inside_table = args.weather_inside_table if args.weather_inside_table is not None else legacy_weather_table
    weather_outside_table = args.weather_outside_table if args.weather_outside_table is not None else legacy_weather_table
    weather_inside_time_col = args.weather_inside_time_col if args.weather_inside_time_col is not None else legacy_weather_time_col
    weather_outside_time_col = args.weather_outside_time_col if args.weather_outside_time_col is not None else legacy_weather_time_col

    df, az_scans, el_scans = build_dataframe(
        rawdata_path,
        args.spectral,
        db_namespace=args.db_namespace,
        telescope=args.telescope,
        azel_source=args.azel_source,
        azel_correction_apply=args.azel_correction_apply,
        altaz_apply=args.altaz_apply,
        encoder_table=args.encoder_table,
        encoder_table_suffix=args.encoder_table_suffix,
        altaz_table=args.altaz_table,
        altaz_table_suffix=args.altaz_table_suffix,
        encoder_time_col=args.encoder_time_col,
        altaz_time_col=args.altaz_time_col,
        spectrometer_time_offset_sec=args.spectrometer_time_offset_sec,
        encoder_shift_sec=args.encoder_shift_sec,
        encoder_az_time_offset_sec=args.encoder_az_time_offset_sec,
        encoder_el_time_offset_sec=args.encoder_el_time_offset_sec,
        encoder_vavg_sec=args.encoder_vavg_sec,
        chopper_wheel=bool(args.chopper_wheel),
        tamb_k=args.tamb_k,
        weather_inside_table=weather_inside_table,
        weather_inside_table_suffix=args.weather_inside_table_suffix,
        weather_inside_time_col=weather_inside_time_col,
        weather_outside_table=weather_outside_table,
        weather_outside_table_suffix=args.weather_outside_table_suffix,
        weather_outside_time_col=weather_outside_time_col,
        tamb_default_k=float(args.tamb_default_k),
        tamb_min_k=float(args.tamb_min_k),
        tamb_max_k=float(args.tamb_max_k),
        outside_default_temperature_c=float(args.outside_default_temperature_c),
        outside_default_pressure_hpa=float(args.outside_default_pressure_hpa),
        outside_default_humidity_pct=float(args.outside_default_humidity_pct),
        outside_temperature_min_c=float(args.outside_temperature_min_c),
        outside_temperature_max_c=float(args.outside_temperature_max_c),
        outside_pressure_min_hpa=float(args.outside_pressure_min_hpa),
        outside_pressure_max_hpa=float(args.outside_pressure_max_hpa),
        outside_humidity_min_pct=float(args.outside_humidity_min_pct),
        outside_humidity_max_pct=float(args.outside_humidity_max_pct),
        chopper_win_sec=float(args.chopper_win_sec),
        chopper_stat=str(args.chopper_stat),
    )

    trim_params = dict(
        vfrac=float(args.trim_vfrac),
        vmin=float(args.trim_vmin),
        gap_fill=int(args.trim_gap),
        min_samples=int(args.trim_min_samples),
        ratio_min=float(args.trim_axis_ratio_min),
        vpercentile=float(args.trim_vpercentile),
        dominant_axis=bool(args.trim_dominant_axis),
        # steady-scan selection
        steady_scan=bool(args.trim_steady_scan),
        use_on_only=bool(args.trim_use_on_only),
        profile_xlim_deg=float(args.profile_xlim_deg),
        xwin_factor=float(args.trim_xwin_factor),
        cross_offset_max_deg=float(args.trim_cross_offset_max_deg),
        speed_min_deg_s=float(args.trim_scan_speed_min_arcsec) / 3600.0,
        steady_cv_max=float(args.trim_steady_cv_max),
    )

    ripple_policy = make_ripple_policy_from_args(args)

    # Choose vertical-axis signal
    ycol = "ta_star" if bool(args.chopper_wheel) else "tp1"
    ytitle = "Ta*" if bool(args.chopper_wheel) else "tp1"

    for i, seg in az_scans.items():
        try:
            if bool(args.debug_plot):
                plot_dxdy_simple(seg, title=f"AZ scan {i}: daz/d_el and Δ per sample", out_png=outdir / f"AZ_dxdy_{i}_{tag}.png",
                                 trim_enabled=bool(args.trim_scan), trim_params=trim_params, which_xcol="daz")
                plot_mount_azel_simple(seg, title=f"AZ scan {i}: Az_true/El_true and Δ per sample", out_png=outdir / f"AZ_mount_azel_{i}_{tag}.png")
        except Exception as e:
            print(f"[error] AZ debug plot scan {i}: {e}")
            if not args.continue_on_error:
                raise

    for i, seg in el_scans.items():
        try:
            if bool(args.debug_plot):
                plot_dxdy_simple(seg, title=f"EL scan {i}: daz/d_el and Δ per sample", out_png=outdir / f"EL_dxdy_{i}_{tag}.png",
                                 trim_enabled=bool(args.trim_scan), trim_params=trim_params, which_xcol="d_el")
                plot_mount_azel_simple(seg, title=f"EL scan {i}: Az_true/El_true and Δ per sample", out_png=outdir / f"EL_mount_azel_{i}_{tag}.png")
        except Exception as e:
            print(f"[error] EL debug plot scan {i}: {e}")
            if not args.continue_on_error:
                raise


    # ------------------------------------------------------------
    # Limb edge fitting (two Gaussians on derivative) + final summary
    # ------------------------------------------------------------
    if bool(args.edge_fit):
        try:
            az_deriv, el_deriv, az_fit, el_fit, results = run_edge_fits_and_summarize(
                az_scans=az_scans,
                el_scans=el_scans,
                ycol=ycol,
                trim_enabled=bool(args.trim_scan),
                trim_params=trim_params,
                strict_deriv=bool(args.strict_deriv),
                ripple_enabled=bool(args.ripple_remove),
                ripple_policy=ripple_policy,
                profile_xlim_deg=float(args.profile_xlim_deg),
                fit_win_deg=float(args.edge_fit_win_deg),
                fit_threshold=float(args.edge_fit_threshold),
                hpbw_init_arcsec=float(args.hpbw_init_arcsec),
            )

            scan_ids = results.pop("__scan_ids__", sorted(set(az_scans.keys()) | set(el_scans.keys()))) if isinstance(results, dict) else []
            az_deriv_raw = results.pop("__az_deriv_raw__", {}) if isinstance(results, dict) else {}
            el_deriv_raw = results.pop("__el_deriv_raw__", {}) if isinstance(results, dict) else {}
            az_profile_raw = results.pop("__az_profile_raw__", {}) if isinstance(results, dict) else {}
            el_profile_raw = results.pop("__el_profile_raw__", {}) if isinstance(results, dict) else {}
            az_profile = results.pop("__az_profile__", {}) if isinstance(results, dict) else {}
            el_profile = results.pop("__el_profile__", {}) if isinstance(results, dict) else {}

            lines = []
            lines.append("=" * 112)
            lines.append(f"Sun scan limb-fit summary: {tag}")
            lines.append(f"  y-axis: {ytitle}  ripple_remove={bool(args.ripple_remove)}  ripple_preset={str(args.ripple_preset)}  chopper_wheel={bool(args.chopper_wheel)}")
            lines.append(f"  edge_fit_win={float(args.edge_fit_win_deg):.3f} deg  threshold={float(args.edge_fit_threshold):.2f}  HPBW_init={float(args.hpbw_init_arcsec):.1f} arcsec")
            lines.append("-" * 112)
            lines.append("per-scan results (rep Az/El = mean of the AZ/EL tracking points nearest zero main-axis offset):")
            lines.append("  scan   rep_Az   rep_El | center_az  HPBW_az  v_az | center_el  HPBW_el  v_el | fitAZ fitEL")
            lines.append("  ----  -------  ------- | ---------  -------  ---- | ---------  -------  ---- | ----- -----")
            for sid in sorted([int(s) for s in scan_ids]):
                rr = results.get(sid, {})
                caz = rr.get("center_az_deg", float("nan")) * 3600.0
                haz = rr.get("hpbw_az_arcsec", float("nan"))
                vaz = rr.get("speed_az_arcsec_s", float("nan"))
                cel = rr.get("center_el_deg", float("nan")) * 3600.0
                hel = rr.get("hpbw_el_arcsec", float("nan"))
                vel = rr.get("speed_el_arcsec_s", float("nan"))
                raz = rr.get("rep_az_deg", float("nan"))
                rel = rr.get("rep_el_deg", float("nan"))
                faz = "Y" if rr.get("fit_ok_az", False) else "N"
                fel = "Y" if rr.get("fit_ok_el", False) else "N"
                lines.append(f"  {sid:4d}  {raz:7.3f}  {rel:7.3f} | {caz:9.1f}  {haz:7.1f}  {vaz:4.0f} | {cel:9.1f}  {hel:7.1f}  {vel:4.0f} |   {faz}     {fel}")

            def _nanmedian_from_keys(key):
                vals = [results.get(s, {}).get(key, float("nan")) for s in scan_ids]
                vals = np.asarray(vals, dtype=float)
                return float(np.nanmedian(vals)) if np.any(np.isfinite(vals)) else float("nan")

            center_az_med = _nanmedian_from_keys("center_az_deg") * 3600.0
            center_el_med = _nanmedian_from_keys("center_el_deg") * 3600.0
            hpbw_az_med = _nanmedian_from_keys("hpbw_az_arcsec")
            hpbw_el_med = _nanmedian_from_keys("hpbw_el_arcsec")
            rep_az_med = _nanmedian_from_keys("rep_az_deg")
            rep_el_med = _nanmedian_from_keys("rep_el_deg")

            lines.append("-" * 112)
            lines.append("median over scans:")
            lines.append(f"  rep_Az = {rep_az_med:7.3f} deg   rep_El = {rep_el_med:7.3f} deg")
            lines.append(f"  center_az = {center_az_med:7.1f} arcsec   HPBW_az = {hpbw_az_med:6.1f} arcsec")
            lines.append(f"  center_el = {center_el_med:7.1f} arcsec   HPBW_el = {hpbw_el_med:6.1f} arcsec")
            lines.append("=" * 112)

            print("\n".join(lines))

            if bool(args.debug_plot):
                save_text_summary("\n".join(lines), outdir / f"summary_text_{tag}.png")

            csv_path = outdir / f"sun_scan_summary_{tag}.csv"
            write_scan_summary_csv(csv_path, scan_ids, results)

            out_pngs = plot_derivative_fits_paginated(
                outdir,
                tag=tag,
                ycol=ycol,
                scan_ids=scan_ids,
                az_profile_raw=az_profile_raw,
                el_profile_raw=el_profile_raw,
                az_profile=az_profile,
                el_profile=el_profile,
                az_deriv_raw=az_deriv_raw,
                el_deriv_raw=el_deriv_raw,
                az_deriv=az_deriv,
                el_deriv=el_deriv,
                az_fit=az_fit,
                el_fit=el_fit,
                results=results,
                xlim_deg=float(args.profile_xlim_deg),
                scans_per_page=int(args.edge_fit_plot_max_scans),
            )
            if out_pngs:
                print(f"[fit] wrote {len(out_pngs)} page(s) of sun_scan_derivative_fits_{tag}_pXXX.png")
            print(f"[fit] wrote {csv_path.name}")
            if bool(args.debug_plot):
                print(f"[fit] wrote summary_text_{tag}.png")
        except Exception as e:
            print(f"[warn] edge fitting failed: {e}")
    print("[done]")
    print(f"  outputs in: {outdir}")
    print(f"  df rows: {len(df)}  az_scans: {len(az_scans)}  el_scans: {len(el_scans)}")
    print(
        "  spectral time basis: "
        f"applied={df.attrs.get('spec_time_basis')} "
        f"suffix={df.attrs.get('spec_time_suffix')} "
        f"fallback={df.attrs.get('spec_time_fallback_field')} "
        f"example={df.attrs.get('spec_time_example')}"
    )


if __name__ == "__main__":
    main()
