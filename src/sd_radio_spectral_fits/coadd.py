# src/sd_radio_spectral_fits/coadd.py
from __future__ import annotations

import builtins
import datetime
import warnings
from typing import Dict, List, Optional, Tuple, Union, Sequence, Any

import numpy as np
import pandas as pd

from .fitsio import Scantable, read_scantable, write_scantable
from .axis import slice_channels
from .baseline import fit_polynomial_baseline
# [MODIFIED] Use new Standardizer
from .regrid_vlsrk import Standardizer, VGrid, make_vgrid
from .doppler import compute_vcorr_series
from .ranges import parse_windows
from .utils import FailPolicy, subtract_windows
# [ADDED] Row selector and Endian utilities
from .scantable_utils import _parse_row_selector, _df_to_native_endian, _resolve_table_timestamps

from .restfreq import apply_restfreq_override

from .tempscale import (
    normalize_tempscal,
    ensure_tempscal_column,
    ensure_beameff_column,
    tempscal_array,
    beameff_array,
    require_beameff,
    is_beameff_mixed,
    representative_beameff,
    append_scale_history,
)


BSL_METADATA_COLUMNS = [
    "BSL_DONE",
    "BSL_APPLIED",
    "BSL_STAGE",
    "BSL_POLY",
    "BSL_WINF",
    "BSL_RMS",
    "BSL_STAT",
    "BSL_NUSED",
    "BSL_COEF",
    "BSL_SCALE",
]


def _purge_bsl_metadata_in_row(row: dict) -> dict:
    for col in BSL_METADATA_COLUMNS:
        row.pop(col, None)
    return row

# -----------------------------------------------------------------------------
# QC Presets & Constants
# -----------------------------------------------------------------------------
QC_PRESETS = {
    "robust":   {"hard": 4.0, "soft_a": 1.0, "soft_p": 4.0, "clip_k": 3.0, "clip_iter": 3},
    "gentle":   {"hard": 0.0, "soft_a": 1.0, "soft_p": 3.0, "clip_k": 3.0, "clip_iter": 2},
    "hardonly": {"hard": 4.0, "soft_a": 0.0, "soft_p": 0.0, "clip_k": 0.0, "clip_iter": 0},
    "noclip":   {"hard": 4.0, "soft_a": 1.0, "soft_p": 4.0, "clip_k": 0.0, "clip_iter": 0},
}


def _parse_qc_spec(spec: str) -> dict:
    config = QC_PRESETS["robust"].copy()
    parts = [p.strip() for p in str(spec).split(";") if p.strip()]
    if not parts:
        return config

    if parts[0] in QC_PRESETS:
        config = QC_PRESETS[parts[0]].copy()
        parts = parts[1:]

    for p in parts:
        if "=" not in p: continue
        k, v = p.split("=", 1)
        k = k.strip().lower()
        v = v.strip()
        try:
            if k == "win": config["win"] = int(v)
            elif k == "winv": config["winv"] = float(v)
            elif k == "hard": config["hard"] = float(v)
            elif k == "soft":
                vals = [float(x) for x in v.split(",")]
                if len(vals) >= 2: config["soft_a"], config["soft_p"] = vals[0], vals[1]
            elif k == "clip":
                vals = [float(x) for x in v.split(",")]
                if len(vals) >= 2: config["clip_k"], config["clip_iter"] = float(vals[0]), int(vals[1])
        except ValueError:
            pass
    return config


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _get_coord_from_header(meta: dict, axis_type: str) -> float | None:
    keys = [axis_type, axis_type.lower(), "OBS" + axis_type, axis_type + "_DEG"]
    for k in keys:
        if k in meta and meta[k] not in (None, ""): return float(meta[k])
    for ax in (2, 3):
        ctype = str(meta.get(f"CTYPE{ax}", "")).upper()
        val = meta.get(f"CRVAL{ax}")
        if val is None: continue
        if axis_type == "RA":
            if ctype.startswith("RA"): return float(val)
        elif axis_type == "DEC":
            if ctype.startswith("DEC"): return float(val)
        elif axis_type == "GLON":
            if ctype.startswith("GLON"): return float(val)
        elif axis_type == "GLAT":
            if ctype.startswith("GLAT"): return float(val)
    return None


def _freq_unit_scale_to_hz(unit: Any) -> float:
    u = str(unit or "").strip().lower()
    if u in ("", "hz"):
        return 1.0
    if u == "khz":
        return 1.0e3
    if u == "mhz":
        return 1.0e6
    if u == "ghz":
        return 1.0e9
    if u == "thz":
        return 1.0e12
    return 1.0


def _normalize_frequency_wcs_units(table: pd.DataFrame, meta: dict) -> tuple[pd.DataFrame, dict]:
    """Normalize frequency-axis WCS units to Hz in both table and meta."""
    out = table.copy()
    m = dict(meta)

    ctype = str(m.get("CTYPE1", "")).strip().upper()
    if not ctype and "CTYPE1" in out.columns:
        vals = out["CTYPE1"].dropna().astype(str).str.strip().str.upper()
        if len(vals) > 0:
            ctype = vals.iloc[0]

    if not ctype.startswith("FREQ"):
        return out, m

    cunit_meta = str(m.get("CUNIT1", "Hz")).strip()
    scale_meta = _freq_unit_scale_to_hz(cunit_meta)
    if scale_meta != 1.0:
        for key in ("CRVAL1", "CDELT1"):
            if key in m and m[key] not in (None, ""):
                m[key] = float(m[key]) * scale_meta
        m["CUNIT1"] = "Hz"

    if "CUNIT1" in out.columns:
        units = out["CUNIT1"].astype(str).str.strip()
        for unit_name in pd.unique(units.dropna()):
            scale = _freq_unit_scale_to_hz(unit_name)
            if scale == 1.0:
                continue
            mask = units.str.lower() == str(unit_name).lower()
            for key in ("CRVAL1", "CDELT1"):
                if key in out.columns:
                    out.loc[mask, key] = pd.to_numeric(out.loc[mask, key], errors="coerce") * scale
        out["CUNIT1"] = "Hz"
    elif "CRVAL1" in out.columns or "CDELT1" in out.columns:
        out["CUNIT1"] = "Hz"

    return out, m


def _infer_coord_frame(table: pd.DataFrame, meta: dict, coord_frame: Optional[str]) -> str:
    if coord_frame not in (None, ""):
        return str(coord_frame).strip().lower()

    meta_frame = str(meta.get("coord_frame", meta.get("COORD_FRAME", ""))).strip().lower()
    if meta_frame in ("icrs", "galactic"):
        return meta_frame

    def _has_pair(lon: str, lat: str) -> bool:
        return (lon in table.columns and lat in table.columns and
                np.isfinite(pd.to_numeric(table[lon], errors="coerce")).any() and
                np.isfinite(pd.to_numeric(table[lat], errors="coerce")).any())

    has_icrs = _has_pair("RA", "DEC") or _has_pair("ra_deg", "dec_deg") or _has_pair("ra", "dec") or _has_pair("OBSRA", "OBSDEC")
    has_gal = _has_pair("GLON", "GLAT") or _has_pair("glon", "glat")

    if has_gal and not has_icrs:
        return "galactic"
    return "icrs"


def _ensure_timestamp_column(table: pd.DataFrame) -> pd.DataFrame:
    """Ensure table has a timezone-naive datetime64[ns] 'timestamp' column.

    Supports both the new standard-like convention
    (DATE-OBS + TIME[offset sec], plus optional MJD/TIMESTAMP)
    and legacy TIME=MJD-days files.
    """
    tab = _df_to_native_endian(table).copy()

    def _to_utc_naive(x):
        ts = pd.to_datetime(x, utc=True, errors="coerce")
        if isinstance(ts, pd.DatetimeIndex):
            return ts.tz_convert("UTC").tz_localize(None)
        return ts.dt.tz_localize(None)

    if "timestamp" in tab.columns:
        tab["timestamp"] = _to_utc_naive(tab["timestamp"])
        return tab

    ts = _resolve_table_timestamps(tab)
    if ts is not None and not ts.isna().all():
        if getattr(ts, "tz", None) is not None:
            ts = ts.tz_convert("UTC").tz_localize(None)
        tab["timestamp"] = ts.to_numpy()
        return tab

    tab["timestamp"] = pd.Timestamp("1970-01-01")
    tab["_time_is_dummy"] = True
    return tab


def _ensure_standard_columns_and_fill(tab: pd.DataFrame, meta: dict) -> pd.DataFrame:
    out = tab.copy()
    
    # ICRS aliases
    if "RA" not in out.columns:
        if "ra_deg" in out.columns: out["RA"] = out["ra_deg"]
        elif "ra" in out.columns: out["RA"] = pd.to_numeric(out["ra"], errors="coerce")
        elif "OBSRA" in out.columns: out["RA"] = pd.to_numeric(out["OBSRA"], errors="coerce")
        else: out["RA"] = np.nan
    if "DEC" not in out.columns:
        if "dec_deg" in out.columns: out["DEC"] = out["dec_deg"]
        elif "dec" in out.columns: out["DEC"] = pd.to_numeric(out["dec"], errors="coerce")
        elif "OBSDEC" in out.columns: out["DEC"] = pd.to_numeric(out["OBSDEC"], errors="coerce")
        else: out["DEC"] = np.nan
    if out["RA"].isna().all():
        ra_val = _get_coord_from_header(meta, "RA")
        if ra_val is not None: out["RA"] = float(ra_val)
    if out["DEC"].isna().all():
        dec_val = _get_coord_from_header(meta, "DEC")
        if dec_val is not None: out["DEC"] = float(dec_val)
        
    # Galactic aliases (for coord_frame='galactic' Doppler correction)
    if "GLON" not in out.columns:
        if "glon" in out.columns: out["GLON"] = pd.to_numeric(out["glon"], errors="coerce")
        else: out["GLON"] = np.nan
    if "GLAT" not in out.columns:
        if "glat" in out.columns: out["GLAT"] = pd.to_numeric(out["glat"], errors="coerce")
        else: out["GLAT"] = np.nan
    if out["GLON"].isna().all():
        glon_val = _get_coord_from_header(meta, "GLON")
        if glon_val is None and str(meta.get("coord_frame", "")).lower() == "galactic":
            glon_val = _get_coord_from_header(meta, "RA")  # backward compatibility
        if glon_val is not None:
            out["GLON"] = float(glon_val)
    if out["GLAT"].isna().all():
        glat_val = _get_coord_from_header(meta, "GLAT")
        if glat_val is None and str(meta.get("coord_frame", "")).lower() == "galactic":
            glat_val = _get_coord_from_header(meta, "DEC")  # backward compatibility
        if glat_val is not None:
            out["GLAT"] = float(glat_val)

    return out


def _get_rest_hz(meta: dict, table: Optional[pd.DataFrame] = None, override_hz: float | None = None) -> float:
    if override_hz is not None and override_hz > 0:
        return float(override_hz)
    for k in ("rest_hz", "RESTFRQ", "RESTFREQ", "restfrq_hz", "restfreq_hz"):
        if k in meta and meta[k] not in (None, ""):
            return float(meta[k])
    if table is not None:
        for k in ("RESTFRQ", "RESTFREQ"):
            if k in table.columns:
                vals = pd.to_numeric(table[k], errors="coerce")
                vals = vals[np.isfinite(vals)]
                if len(vals) > 0:
                    return float(vals.iloc[0])
    return 0.0

def _resolve_vcorr_lonlat_columns(table: pd.DataFrame, coord_frame: str) -> tuple[str | None, str | None]:
    frame = str(coord_frame or "icrs").strip().lower()
    if frame == "galactic":
        pairs = (("GLON", "GLAT"), ("glon", "glat"))
    else:
        pairs = (("RA", "DEC"), ("ra_deg", "dec_deg"), ("ra", "dec"), ("OBSRA", "OBSDEC"))

    for lon_col, lat_col in pairs:
        if lon_col in table.columns and lat_col in table.columns:
            return lon_col, lat_col
    return None, None


def _resolve_specsys(meta: dict, table: pd.DataFrame | None = None) -> str:
    if table is not None:
        for key in ("SPECSYS", "SSYSOBS"):
            if key in table.columns:
                ser = table[key]
                vals = ser.dropna().astype(str).str.strip()
                vals = vals[vals != ""]
                uniq = pd.unique(vals)
                if len(uniq) > 1:
                    raise ValueError(f"Input contains mixed {key} values: {list(uniq)}")
                if len(uniq) == 1:
                    return str(uniq[0]).upper()
    for key in ("SPECSYS", "SSYSOBS"):
        val = meta.get(key)
        if val not in (None, ""):
            return str(val).strip().upper()
    return ""


def _resolve_ssysobs(meta: dict, table: pd.DataFrame | None = None) -> str:
    if table is not None and "SSYSOBS" in table.columns:
        ser = table["SSYSOBS"]
        vals = ser.dropna().astype(str).str.strip()
        vals = vals[vals != ""]
        uniq = pd.unique(vals)
        if len(uniq) == 1:
            return str(uniq[0]).upper()
    val = meta.get("SSYSOBS")
    if val not in (None, ""):
        return str(val).strip().upper()
    if table is not None and "SPECSYS" in table.columns:
        ser = table["SPECSYS"]
        vals = ser.dropna().astype(str).str.strip()
        vals = vals[vals != ""]
        uniq = pd.unique(vals)
        if len(uniq) == 1:
            return str(uniq[0]).upper()
    val = meta.get("SPECSYS")
    if val not in (None, ""):
        return str(val).strip().upper()
    return ""


def _resolve_ssysobs_series(table: pd.DataFrame, meta: dict) -> pd.Series:
    n = len(table)
    if "SSYSOBS" in table.columns:
        out = table["SSYSOBS"].astype(object).copy()
    else:
        out = pd.Series([None] * n, index=table.index, dtype=object)

    if "SPECSYS" in table.columns:
        row_specsys = table["SPECSYS"].astype(object)
    else:
        row_specsys = pd.Series([None] * n, index=table.index, dtype=object)

    meta_ssysobs = meta.get("SSYSOBS", None)
    meta_specsys = meta.get("SPECSYS", None)

    def _is_missing(v: Any) -> bool:
        if v is None:
            return True
        if isinstance(v, float) and not np.isfinite(v):
            return True
        return str(v).strip() == ""

    for idx in out.index:
        cur = out.at[idx]
        if _is_missing(cur):
            if meta_ssysobs not in (None, ""):
                cur = meta_ssysobs
            else:
                row_val = row_specsys.at[idx]
                if not _is_missing(row_val):
                    cur = row_val
                elif meta_specsys not in (None, ""):
                    cur = meta_specsys
                else:
                    cur = ""
        out.at[idx] = str(cur).strip().upper() if cur not in (None, "") else ""
    return out


def _unique_nonempty_upper(values) -> str:
    ser = pd.Series(values, dtype=object).dropna().astype(str).str.strip()
    ser = ser[ser != ""]
    uniq = pd.unique(ser)
    if len(uniq) == 1:
        return str(uniq[0]).upper()
    return ""


def _fill_constant_column(out: pd.DataFrame, key: str, value: Any) -> pd.DataFrame:
    if value in (None, ""):
        return out
    if key not in out.columns:
        out[key] = value
        return out
    col = out[key]
    if pd.api.types.is_numeric_dtype(col):
        mask = col.isna()
    else:
        mask = col.isna() | (col.astype(str).str.strip() == "")
    if bool(mask.any()):
        out.loc[mask, key] = value
    return out


def _ensure_standard_spectral_columns(table: pd.DataFrame, meta: dict) -> pd.DataFrame:
    out = table.copy()

    rest = None
    for key in ("RESTFRQ", "RESTFREQ"):
        if meta.get(key) not in (None, ""):
            rest = float(meta[key])
            break
    if rest is None:
        for key in ("RESTFRQ", "RESTFREQ"):
            if key in out.columns:
                vals = pd.to_numeric(out[key], errors="coerce")
                vals = vals[np.isfinite(vals)]
                if len(vals) > 0:
                    rest = float(vals.iloc[0])
                    break

    cdelt1 = meta.get("CDELT1", meta.get("CD1_1", np.nan))
    defaults = {
        "CTYPE1": meta.get("CTYPE1", "FREQ"),
        "CUNIT1": meta.get("CUNIT1", "Hz"),
        "CRVAL1": meta.get("CRVAL1", np.nan),
        "CDELT1": cdelt1,
        "CRPIX1": meta.get("CRPIX1", 1.0),
        "SPECSYS": meta.get("SPECSYS", meta.get("SSYSOBS", "")),
        "SSYSOBS": meta.get("SSYSOBS", meta.get("SPECSYS", "")),
    }
    for key, value in defaults.items():
        out = _fill_constant_column(out, key, value)

    if rest is not None:
        out = _fill_constant_column(out, "RESTFRQ", rest)
        out = _fill_constant_column(out, "RESTFREQ", rest)

    if "RESTFRQ" not in out.columns and "RESTFREQ" in out.columns:
        out["RESTFRQ"] = pd.to_numeric(out["RESTFREQ"], errors="coerce")
    if "RESTFREQ" not in out.columns and "RESTFRQ" in out.columns:
        out["RESTFREQ"] = pd.to_numeric(out["RESTFRQ"], errors="coerce")
    if "SSYSOBS" not in out.columns and "SPECSYS" in out.columns:
        out["SSYSOBS"] = out["SPECSYS"]
    if "SPECSYS" not in out.columns and "SSYSOBS" in out.columns:
        out["SPECSYS"] = out["SSYSOBS"]

    return out


def _velocity_candidates_mps(table: pd.DataFrame, v_corr_col: str) -> list[tuple[str, float]]:
    cand: list[tuple[str, float]] = [("VELOSYS", 1.0)]
    ku = str(v_corr_col or "").strip().upper()
    if ku:
        scale = 1000.0 if (ku.endswith("_KMS") or ku == "V_CORR_KMS") else 1.0
        if ku != "VELOSYS":
            cand.append((v_corr_col, scale))
    if ku != "VFRAME":
        cand.append(("VFRAME", 1.0))
    if ku != "V_CORR_KMS":
        cand.append(("V_CORR_KMS", 1000.0))
    uniq = []
    seen = set()
    for key, scale in cand:
        if key in table.columns and key not in seen:
            uniq.append((key, scale))
            seen.add(key)
    return uniq


def _extract_existing_velocity_mps(table: pd.DataFrame, v_corr_col: str) -> np.ndarray | None:
    candidates = _velocity_candidates_mps(table, v_corr_col)
    arrays = []
    names = []
    for key, scale in candidates:
        vals = pd.to_numeric(table[key], errors="coerce").to_numpy(dtype=float) * scale
        if np.isfinite(vals).all():
            arrays.append(vals)
            names.append(key)
    if not arrays:
        return None
    ref = arrays[0]
    for key, arr in zip(names[1:], arrays[1:]):
        if not np.allclose(ref, arr, rtol=0.0, atol=1.0e-6, equal_nan=False):
            raise ValueError(f"Velocity columns disagree: {names[0]} vs {key}")
    return ref


def _validate_no_unapplied_velocity(table: pd.DataFrame, specsys: str, v_corr_col: str) -> None:
    if "TOPO" in specsys:
        return
    for key, scale in _velocity_candidates_mps(table, v_corr_col):
        vals = pd.to_numeric(table[key], errors="coerce").to_numpy(dtype=float) * scale
        finite = vals[np.isfinite(vals)]
        if finite.size > 0 and not np.allclose(finite, 0.0, rtol=0.0, atol=1.0e-6):
            raise ValueError(
                f"Input has SPECSYS={specsys} but contains non-zero unapplied velocity column '{key}'."
            )



def _group_positions(table: pd.DataFrame, *, pos_col: str, tol_arcsec: Optional[float] = None) -> Dict[int, np.ndarray]:
    if pos_col in table.columns:
        groups: Dict[int, np.ndarray] = {}
        for pid, idx in table.groupby(pos_col).groups.items():
            groups[int(pid)] = np.asarray(list(idx), int)
        return groups
        
    if tol_arcsec is None:
        raise ValueError(f"table lacks '{pos_col}' and tol_arcsec not given")
    if "RA" not in table.columns or "DEC" not in table.columns:
        raise ValueError("need RA/DEC to group without pos_id")
        
    # 許容距離（秒角）をラジアンに変換
    tol_rad = np.deg2rad(float(tol_arcsec) / 3600.0)

    ra = pd.to_numeric(table["RA"], errors="coerce").to_numpy(dtype=float)
    dec = pd.to_numeric(table["DEC"], errors="coerce").to_numpy(dtype=float)
    if (not np.all(np.isfinite(ra))) or (not np.all(np.isfinite(dec))):
        raise ValueError("RA/DEC contains non-finite values; cannot group positions.")

    # 座標をラジアンに変換
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)
    n = len(ra)
    if n > 1:
        print("\n--- [DEBUG] プログラムが計算した実際の座標間距離 ---")
        for i in range(min(5, n)):  # 最初の5点だけ確認
            d_ra = ra_rad - ra_rad[i]
            d_dec = dec_rad - dec_rad[i]
            a = np.sin(d_dec / 2.0)**2 + np.cos(dec_rad[i]) * np.cos(dec_rad) * np.sin(d_ra / 2.0)**2
            a = np.clip(a, 0.0, 1.0)
            dist_arcsec = np.rad2deg(2.0 * np.arcsin(np.sqrt(a))) * 3600.0
            
            # 自分自身（距離0）を除外して最小値を探す
            dist_arcsec[i] = np.inf
            nearest_idx = np.argmin(dist_arcsec)
            min_dist = dist_arcsec[nearest_idx]
            
            print(f"Point {i} (RA={ra[i]:.5f}, DEC={dec[i]:.5f}):")
            print(f"  -> 最も近い別のデータ(Index {nearest_idx})までの距離 = {min_dist:.2f} 秒角")
        print("----------------------------------------------------\n")
    # ====================================================================

    visited = np.zeros(n, dtype=bool)
    groups: Dict[int, np.ndarray] = {}
    gid = 0
    
    for i in range(n):
        if visited[i]:
            continue
            
        d_ra = ra_rad - ra_rad[i]
        d_dec = dec_rad - dec_rad[i]
        a = np.sin(d_dec / 2.0)**2 + np.cos(dec_rad[i]) * np.cos(dec_rad) * np.sin(d_ra / 2.0)**2
        a = np.clip(a, 0.0, 1.0)
        dist_rad = 2.0 * np.arcsin(np.sqrt(a))
        
        # 基準点から tol_rad 以内にあるデータをグループ化
        mask = (dist_rad <= tol_rad) & (~visited)
        
        groups[gid] = np.where(mask)[0]
        visited[mask] = True
        gid += 1
        
    return groups
    

def _group_indices(table: pd.DataFrame, mode: str, pos_col: str, pos_tol_arcsec: float | None) -> Dict[int, np.ndarray]:
    """
    mode: 'position' | 'scan' | 'intgrp'

    For mode='scan', spectra are grouped by the composite key
    (SCAN, OBSMODE, FDNUM, IFNUM, PLNUM).
    This prevents spectra from different spectrometers/IFs/polarization paths
    from being coadded together merely because they share the same SCAN number.

    Backward compatibility:
      - SCAN is still required.
      - If OBSMODE is absent, a single empty-string mode is assumed.
      - If FDNUM/IFNUM/PLNUM are absent, they are treated as 0.
    """
    mode = str(mode).lower().strip()
    
    if mode == "position":
        return _group_positions(table, pos_col=pos_col, tol_arcsec=pos_tol_arcsec)
    
    elif mode == "scan":
        # Required scan identifier
        scan_col = None
        for c in ["SCAN", "SCANID", "scan", "scanid", "scan_id"]:
            if c in table.columns:
                scan_col = c
                break
        if not scan_col:
            raise ValueError("Group mode 'scan' requested but no SCAN/SCANID column found in table.")

        # Optional grouping columns for instrument separation.
        obsmode_col = None
        for c in ["OBSMODE", "obsmode", "OBS_MODE", "obs_mode"]:
            if c in table.columns:
                obsmode_col = c
                break

        def _find_first(existing: list[str]) -> str | None:
            for c in existing:
                if c in table.columns:
                    return c
            return None

        fdnum_col = _find_first(["FDNUM", "fdnum"])
        ifnum_col = _find_first(["IFNUM", "ifnum"])
        plnum_col = _find_first(["PLNUM", "plnum"])

        # Build a compact grouping frame using positional row indices.
        work = pd.DataFrame({"_rowpos": np.arange(len(table), dtype=int)})
        work["_scan"] = pd.to_numeric(table[scan_col], errors="raise").astype(int)
        if obsmode_col is not None:
            work["_obsmode"] = table[obsmode_col].astype(str).str.strip().str.upper()
        else:
            work["_obsmode"] = ""

        input_id_col = _find_first(["_INPUT_ID", "INPUT_ID", "input_id"])

        def _optional_int_series(colname: str | None) -> pd.Series:
            if colname is None:
                return pd.Series(np.zeros(len(table), dtype=int), index=table.index)
            return pd.to_numeric(table[colname], errors="raise").fillna(0).astype(int)

        work["_fdnum"] = _optional_int_series(fdnum_col).to_numpy()
        work["_ifnum"] = _optional_int_series(ifnum_col).to_numpy()
        work["_plnum"] = _optional_int_series(plnum_col).to_numpy()
        if input_id_col is not None:
            work["_input_id"] = _optional_int_series(input_id_col).to_numpy()
        else:
            work["_input_id"] = 0

        groups: Dict[int, np.ndarray] = {}
        group_cols = ["_input_id", "_scan", "_obsmode", "_fdnum", "_ifnum", "_plnum"]
        for gid, (_key, sub) in enumerate(work.groupby(group_cols, sort=False, dropna=False)):
            groups[gid] = sub["_rowpos"].to_numpy(dtype=int)
        return groups
    
    elif mode == "intgrp":
        intgrp_col = None
        for c in ["INTGRP", "intgrp", "INTEGRATION_GROUP"]:
            if c in table.columns:
                intgrp_col = c
                break
        if not intgrp_col:
             raise ValueError("Group mode 'intgrp' requested but no INTGRP column found in table.")
        
        groups = {}
        for gid, idx in table.groupby(intgrp_col).groups.items():
            groups[int(gid)] = np.asarray(list(idx), int)
        return groups
    
    else:
        raise ValueError(f"Unknown group mode: {mode}")


# -----------------------------------------------------------------------------
# Fast rolling mean (NaN-safe) for 2D, segment-safe
# -----------------------------------------------------------------------------

def _mask_to_segments(mask: np.ndarray) -> list[tuple[int, int]]:
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0: return []
    m = mask.astype(np.int8)
    dm = np.diff(m)
    starts = list(np.where(dm == 1)[0] + 1)
    ends = list(np.where(dm == -1)[0] + 1)
    if mask[0]: starts = [0] + starts
    if mask[-1]: ends = ends + [mask.size]
    return [(int(s), int(e)) for s, e in zip(starts, ends) if e > s]


def _fast_nan_boxcar_2d(data: np.ndarray, width: int) -> np.ndarray:
    if width <= 1: return data
    if width % 2 == 0: width += 1
    half = width // 2
    data = np.asarray(data, dtype=float)
    n, c = data.shape
    valid = np.isfinite(data)
    filled = np.where(valid, data, 0.0)
    cum_val = np.concatenate([np.zeros((n, 1), dtype=float), np.cumsum(filled, axis=1)], axis=1)
    cum_cnt = np.concatenate([np.zeros((n, 1), dtype=float), np.cumsum(valid.astype(np.int64), axis=1)], axis=1)
    idx = np.arange(c)
    left = np.clip(idx - half, 0, c)
    right = np.clip(idx + half + 1, 0, c)
    sum_val = cum_val[:, right] - cum_val[:, left]
    sum_cnt = cum_cnt[:, right] - cum_cnt[:, left]
    out = np.full((n, c), np.nan, dtype=float)
    m = sum_cnt > 0
    out[m] = sum_val[m] / sum_cnt[m]
    return out


def _mad_along_axis(a: np.ndarray, axis: int) -> np.ndarray:
    med = np.nanmedian(a, axis=axis, keepdims=True)
    mad = np.nanmedian(np.abs(a - med), axis=axis)
    return mad


def _evaluate_rms_for_weight(
    v_axis: np.ndarray, spec: np.ndarray, windows: List[Tuple[float, float]],
    poly_order: int, bin_width: int = 1,
) -> Tuple[float, float]:
    """指定された次数でベースラインを差し引き、任意でbinningした後の残差からRMSと重みを計算する"""
    _, base, info = fit_polynomial_baseline(
        v_axis, spec, base_windows=windows, poly_order=int(poly_order), iter_max=0
    )
    mask = info.get("mask", np.zeros_like(spec, dtype=bool))
    if np.count_nonzero(mask) < builtins.max(2, int(poly_order) + 1):
        return 0.0, np.nan
        
    resid = spec - base
    resid_valid = resid[mask]
    
    # [MODIFIED] 残差に対して指定幅のBoxcarスムージングを適用
    if int(bin_width) > 1:
        resid_valid = _fast_nan_boxcar_2d(resid_valid[None, :], int(bin_width))[0]
        
    resid_fin = resid_valid[np.isfinite(resid_valid)]
    if resid_fin.size < 2:
        return 0.0, np.nan
        
    rms = float(np.std(resid_fin, ddof=1))
    if not np.isfinite(rms) or rms <= 0:
        return 0.0, np.nan
        
    return 1.0 / (rms * rms), rms


# -----------------------------------------------------------------------------
# Core QC
# -----------------------------------------------------------------------------

def _compute_qc_weights(
    spec_stack: np.ndarray,
    v_axis: np.ndarray,
    rms_windows: List[Tuple[float, float]],
    qc_config: dict,
    dv_val: float,
) -> Tuple[np.ndarray, dict]:
    n_dump, n_chan = spec_stack.shape
    win_ch = 51
    if "win" in qc_config: win_ch = int(qc_config["win"])
    elif "winv" in qc_config: win_ch = int(round(float(qc_config["winv"]) / float(dv_val)))
    else:
        total_width_v = sum(max(w) - min(w) for w in rms_windows)
        win_ch = int(round((total_width_v / 8.0) / float(dv_val)))
        win_ch = builtins.max(51, builtins.min(301, win_ch))
    if win_ch % 2 == 0: win_ch += 1
    if win_ch > n_chan: win_ch = builtins.max(3, (n_chan // 2) * 2 - 1)

    mask_rms = np.zeros(n_chan, dtype=bool)
    for (w0, w1) in rms_windows:
        lo, hi = (w0, w1) if w0 <= w1 else (w1, w0)
        mask_rms |= (v_axis >= lo) & (v_axis <= hi)

    if np.count_nonzero(mask_rms) < 2:
        return np.zeros(n_dump, dtype=float), {
            "n_dump": n_dump, "n_keep": 0, "win": win_ch,
            "q_median": np.nan, "q_mad": np.nan, "q_max": np.nan,
            "w_min": 0, "w_max": 0, "w_median": 0,
        }

    T_full = np.asarray(spec_stack, dtype=float).copy()
    T_full[:, ~mask_rms] = np.nan
    b_full = np.full_like(T_full, np.nan)
    segments = _mask_to_segments(mask_rms)
    for s, e in segments:
        b_full[:, s:e] = _fast_nan_boxcar_2d(T_full[:, s:e], win_ch)

    r_full = T_full - b_full
    diff_full = np.diff(r_full, axis=1)
    sigma_white = _mad_along_axis(diff_full, axis=1) / 1.41421356
    sigma_white = np.where(np.isfinite(sigma_white), sigma_white, np.inf)
    sigma_white = np.maximum(sigma_white, 1e-12)
    sigma_ripple = _mad_along_axis(b_full, axis=1)
    sigma_ripple = np.where(np.isfinite(sigma_ripple), sigma_ripple, np.inf)

    with np.errstate(divide="ignore", invalid="ignore"):
        q = sigma_ripple / sigma_white
    q[~np.isfinite(q)] = 1e9

    valid_q_mask = np.isfinite(q) & np.isfinite(sigma_white) & np.isfinite(sigma_ripple)
    valid_q_mask &= (sigma_white > 0) & (sigma_white < np.inf) & (sigma_ripple < np.inf)

    if np.count_nonzero(valid_q_mask) >= 3:
        q_for_stat = q[valid_q_mask]
        q_median = float(np.nanmedian(q_for_stat))
        q_mad = float(np.nanmedian(np.abs(q_for_stat - q_median)))
    else:
        q_median = float(np.nanmedian(q))
        q_mad = float(np.nanmedian(np.abs(q - q_median)))

    hard_thresh = float(qc_config.get("hard", 0))
    keep_mask = np.ones(n_dump, dtype=bool)
    if hard_thresh > 0 and q_mad > 1e-12:
        cutoff = q_median + hard_thresh * q_mad
        keep_mask &= (q <= cutoff)

    w_stat = np.zeros(n_dump, dtype=float)
    valid_sig = np.isfinite(sigma_white) & (sigma_white > 0) & (sigma_white < np.inf)
    w_stat[valid_sig] = 1.0 / (sigma_white[valid_sig] ** 2)

    if np.any(w_stat > 0):
        w_clip_val = float(np.nanpercentile(w_stat[w_stat > 0], 99.5))
        w_stat = np.clip(w_stat, 0, w_clip_val)

    soft_a = float(qc_config.get("soft_a", 0))
    soft_p = float(qc_config.get("soft_p", 0))
    if soft_a > 0:
        q0 = q_median * soft_a
        w_qual = np.zeros(n_dump, dtype=float)
        if q0 > 0:
            ratio = q / q0
            safe = (ratio >= 0) & (ratio < 100)
            w_qual[safe] = 1.0 / (1.0 + ratio[safe] ** soft_p)
            w_qual[ratio >= 100] = 0.0
        w_final = w_stat * w_qual
    else:
        w_final = w_stat

    w_final[~keep_mask] = 0.0
    w_final[~np.isfinite(w_final)] = 0.0

    if np.count_nonzero(w_final > 0) == 0 and n_dump > 0:
        best_idx = int(np.argmax(w_stat))
        if w_stat[best_idx] > 0: w_final[best_idx] = w_stat[best_idx]

    stats = {
        "n_dump": int(n_dump), "n_keep": int(np.count_nonzero(w_final > 0)),
        "q_median": q_median, "q_mad": q_mad,
        "q_max": float(np.nanmax(q)) if n_dump > 0 else 0.0, "win": int(win_ch),
        "w_min": float(np.nanmin(w_final[w_final > 0])) if np.any(w_final > 0) else 0.0,
        "w_max": float(np.nanmax(w_final)) if n_dump > 0 else 0.0,
        "w_median": float(np.nanmedian(w_final[w_final > 0])) if np.any(w_final > 0) else 0.0,
    }
    return w_final, stats


def _robust_channel_coadd(
    spectra: np.ndarray, weights: np.ndarray, clip_k: float, clip_iter: int,
) -> Tuple[np.ndarray, np.ndarray]:
    n_dump, n_chan = spectra.shape
    valid_data = np.isfinite(spectra)
    w_2d = weights[:, None] * valid_data
    sum_w = np.sum(w_2d, axis=0)
    sum_d = np.sum(np.where(valid_data, spectra, 0.0) * w_2d, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"): mean_spec = sum_d / sum_w
    if clip_iter <= 0 or clip_k <= 0: return mean_spec, sum_w

    current_mask = valid_data & (weights[:, None] > 0)
    for _ in range(int(clip_iter)):
        resid = spectra - mean_spec[None, :]
        resid_masked = resid.copy()
        resid_masked[~current_mask] = np.nan
        med_resid = np.nanmedian(resid_masked, axis=0)
        abs_dev = np.abs(resid_masked - med_resid[None, :])
        mad = np.nanmedian(abs_dev, axis=0)
        sigma = mad * 1.4826
        threshold = clip_k * sigma[None, :]
        z = threshold == 0
        if np.any(z):
            threshold = threshold.copy()
            threshold[z] = np.inf
        bad = np.isfinite(abs_dev) & (abs_dev > threshold)
        if np.count_nonzero(bad) == 0: break
        current_mask &= (~bad)
        w_eff = weights[:, None] * current_mask
        sum_w = np.sum(w_eff, axis=0)
        sum_d = np.sum(np.where(current_mask, spectra, 0.0) * w_eff, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"): mean_spec = sum_d / sum_w
    return mean_spec, sum_w


# -----------------------------------------------------------------------------
# Main Logic: Velocity Coadd (Standardizer Integration)
# -----------------------------------------------------------------------------

def run_velocity_coadd(
    inputs: Sequence[Union[str, Scantable]],
    output_path: Optional[str] = None,
    *,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    mode: str = "uniform",
    group_mode: str = "position",
    pos_col: str = "pos_id",
    pos_tol_arcsec: Optional[float] = None,
    v_corr_col: str = "VELOSYS",
    coord_frame: Optional[str] = None,
    vcorr_chunk_sec: Optional[float] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    dv: Optional[float] = None,
    allow_outside_overlap: bool = False,
    axis_type: str = "freq",
    rest_freq: Optional[float] = None,
    baseline_vwin: Union[str, List[str], None] = None,  # [MODIFIED]
    baseline_poly: int = 0,
    baseline_iter_max: int = 0,
    baseline_iter_sigma: float = 3.0,
    post_baseline_mode: Optional[str] = None,
    post_baseline_vwin: Union[str, List[str], None] = "inherit",
    post_baseline_poly: Union[str, int] = "inherit",
    post_baseline_iter_max: Union[str, int] = "inherit",
    post_baseline_iter_sigma: Union[str, float] = "inherit",
    rms_vwin: Union[str, List[str], None] = None,       # [MODIFIED]
    rms_poly: int = 1,   # [MODIFIED] 指定がない場合のデフォルトを 1 に変更
    rms_bin: int = 1,
    weight_zero_policy: str = "error",
    coadd_qc: Optional[str] = None,
    line_vwin: Union[str, List[str], None] = None,      # [MODIFIED]
    block_size: int = 0,
    max_dumps: int = 0,
    ch_start: Optional[int] = None,
    ch_stop: Optional[int] = None,
    on_fail: str = "exit",
    overwrite: bool = True,
    out_scale: str = "TA*",
    normalize_if_mixed: str = "auto",
    beameff_tol: float = 1e-6,
    sigma_scale: str = "TA*",
    verbose: bool = True,  # [ADDED] 実行時の情報出力を制御するフラグ
) -> Scantable:
    """
    Run velocity regridding and coaddition pipeline using Standardizer.
    
    group_mode: 'position' (stack repeats), 'scan' (reduce raw scans), 'intgrp'
    """
    if isinstance(inputs, (str, Scantable)):
        inputs = [inputs]

    policy = FailPolicy(on_fail)

    # [ADDED] 文字列が直接渡された場合は自動的にリスト化する（利便性向上）
    if isinstance(baseline_vwin, str): baseline_vwin = [baseline_vwin]
    post_baseline_mode_norm = None if post_baseline_mode is None else str(post_baseline_mode).strip().lower()
    if post_baseline_mode_norm in {"", "none"}:
        post_baseline_mode_norm = None
    if post_baseline_mode_norm not in {None, "inherit_all"}:
        raise ValueError("post_baseline_mode must be None or 'inherit_all'.")
    if isinstance(post_baseline_vwin, str):
        if str(post_baseline_vwin).strip().lower() != "inherit":
            post_baseline_vwin = [post_baseline_vwin]
        else:
            post_baseline_vwin = "inherit"
    if isinstance(rms_vwin, str): rms_vwin = [rms_vwin]
    if isinstance(line_vwin, str): line_vwin = [line_vwin]
    
    # [MODIFIED] Mutual Exclusion Check
    if baseline_vwin and rms_vwin:
        raise ValueError(
            "Cannot specify both 'baseline_vwin' and 'rms_vwin' simultaneously. "
            "Use 'baseline_vwin' for baseline subtraction and weighting, or "
            "'rms_vwin' for RMS weighting without subtraction."
        )

    weight_zero_policy = str(weight_zero_policy).strip().lower()
    if weight_zero_policy not in {"error", "drop", "impute_median"}:
        raise ValueError("weight_zero_policy must be 'error', 'drop', or 'impute_median'.")

    # [ADDED] Mode Validation
    valid_modes = ["uniform", "rms", "rms_weight"]
    mode_lower = str(mode).lower().strip()
    if mode_lower not in valid_modes:
        raise ValueError(f"Invalid mode: '{mode}'. Must be one of {valid_modes}.")
    is_rms_mode = mode_lower in ["rms", "rms_weight"]

    # [ADDED] QC vs Uniform Conflict Check
    if coadd_qc is not None and mode_lower == "uniform":
         raise ValueError(
             "Cannot use QC mode (--coadd-qc) when mode is 'uniform'. "
             "QC mode intrinsically applies statistical weighting. Please set mode to 'rms' or remove QC."
         )

    out_scale = normalize_tempscal(out_scale, default="TA*")
    sigma_scale = normalize_tempscal(sigma_scale, default="TA*")

    if sigma_scale != "TA*":
        raise NotImplementedError(
            "sigma_scale other than 'TA*' is not implemented. "
            "Provide TA*-scale RMS/weights, and request out_scale='TR*' to convert at output."
        )
    normalize_if_mixed = str(normalize_if_mixed).lower().strip()
    

    # 1. Read Inputs into a single Scantable (Raw / VLA)
    metas, data_list, tabs, hists = [], [], [], []
    rest_inputs_hz: list[float] = []
    
    for fi, inp in enumerate(inputs):
        if isinstance(inp, str):
            sc = read_scantable(inp)
        else:
            sc = inp
        
        meta_i = dict(sc.meta)
        data_i = sc.data
        table_i = _df_to_native_endian(sc.table).copy()
        hist_i = dict(sc.history)
        
        table_i = _ensure_timestamp_column(table_i)
        table_i = _ensure_standard_columns_and_fill(table_i, meta_i)

        if rest_freq is not None:
            apply_restfreq_override(meta_i, table_i, float(rest_freq), require_wcs_for_vrad=True)

        table_i, meta_i = _normalize_frequency_wcs_units(table_i, meta_i)
        table_i = _ensure_standard_spectral_columns(table_i, meta_i)
        specsys = _resolve_specsys(meta_i, table_i)
        ssysobs_i = _resolve_ssysobs(meta_i, table_i)
        if not specsys:
            raise ValueError("Each coadd input must define SPECSYS or SSYSOBS.")
        rest_i = _get_rest_hz(meta_i, table_i, override_hz=rest_freq)
        if not np.isfinite(rest_i) or rest_i <= 0:
            raise ValueError("Each coadd input must define RESTFRQ/RESTFREQ > 0.")
        rest_inputs_hz.append(float(rest_i))
        is_topo = "TOPO" in specsys

        existing_vel_mps = _extract_existing_velocity_mps(table_i, v_corr_col)

        if is_topo:
            if existing_vel_mps is not None:
                velosys_mps = np.asarray(existing_vel_mps, dtype=float)
            else:
                if ("_time_is_dummy" in table_i.columns) and bool(np.any(table_i["_time_is_dummy"].to_numpy(dtype=bool))):
                    raise ValueError("TOPOCENT input requires VELOSYS/VFRAME or valid timestamps for velocity correction.")

                times_utc = pd.DatetimeIndex(pd.to_datetime(table_i["timestamp"], utc=True, errors="coerce"))
                if bool(times_utc.isna().any()):
                    raise ValueError("TOPOCENT input contains invalid timestamps; cannot compute VELOSYS.")

                frame_used = _infer_coord_frame(table_i, meta_i, coord_frame)
                lon_col, lat_col = _resolve_vcorr_lonlat_columns(table_i, frame_used)
                if lon_col is None or lat_col is None:
                    raise ValueError("TOPOCENT input lacks usable coordinate columns for VELOSYS computation.")

                ra_vals = pd.to_numeric(table_i[lon_col], errors="coerce").to_numpy(dtype=float)
                dec_vals = pd.to_numeric(table_i[lat_col], errors="coerce").to_numpy(dtype=float)
                if (not np.all(np.isfinite(ra_vals))) or (not np.all(np.isfinite(dec_vals))):
                    raise ValueError("TOPOCENT input contains non-finite coordinates; cannot compute VELOSYS.")

                should_decimate = (vcorr_chunk_sec is not None and vcorr_chunk_sec > 0 and len(times_utc) > 2)
                if should_decimate:
                    t_vals = times_utc.values.astype(float) / 1e9
                    t_start = t_vals[0]
                    t_end = t_vals[-1]
                    duration = t_end - t_start
                    n_steps = int(np.ceil(duration / vcorr_chunk_sec))
                    if n_steps < 1:
                        n_steps = 1
                    grid_times = np.linspace(t_start, t_end, num=n_steps + 1)
                    idxs = np.searchsorted(t_vals, grid_times)
                    idxs = np.clip(idxs, 0, len(t_vals) - 1)
                    idxs = np.unique(idxs)
                    vc_sparse = compute_vcorr_series(
                        times=times_utc[idxs],
                        ra_deg=ra_vals[idxs],
                        dec_deg=dec_vals[idxs],
                        meta=meta_i,
                        coord_frame=str(frame_used),
                    )
                    vc_kms = np.interp(t_vals, t_vals[idxs], vc_sparse)
                else:
                    vc_kms = compute_vcorr_series(
                        times=times_utc,
                        ra_deg=ra_vals,
                        dec_deg=dec_vals,
                        meta=meta_i,
                        coord_frame=str(frame_used),
                    )
                velosys_mps = np.asarray(vc_kms, dtype=float) * 1000.0

            table_i = table_i.copy()
            table_i["VELOSYS"] = velosys_mps
            table_i["VFRAME"] = velosys_mps
            if v_corr_col not in ("VELOSYS", "VFRAME"):
                table_i[v_corr_col] = velosys_mps
        else:
            _validate_no_unapplied_velocity(table_i, specsys, v_corr_col)
            table_i = table_i.copy()
            table_i["VELOSYS"] = 0.0
            if "VFRAME" in table_i.columns:
                table_i["VFRAME"] = 0.0
            if v_corr_col not in ("VELOSYS", "VFRAME") and v_corr_col in table_i.columns:
                table_i[v_corr_col] = 0.0

        table_i["SPECSYS"] = specsys if specsys else table_i.get("SPECSYS", "")
        table_i["SSYSOBS"] = _resolve_ssysobs_series(table_i, meta_i).to_numpy(dtype=object)
        table_i["_local_idx"] = np.arange(len(table_i), dtype=int)
        
        if isinstance(data_i, list):
            if len(data_i) > 0:
                nchan_vec = np.fromiter((len(x) for x in data_i), dtype=int, count=len(data_i))
            else:
                nchan_vec = np.zeros(len(table_i), dtype=int)
            data_list.extend(data_i)
        else:
            if ch_start is not None or ch_stop is not None:
                meta_i, data_i = slice_channels(meta_i, data_i, ch_start, ch_stop)
            nchan_vec = np.full(len(table_i), int(data_i.shape[1]), dtype=int)
            data_list.extend(list(data_i))

        table_i = table_i.copy()
        table_i["_INPUT_ID"] = int(len(tabs))

        metas.append(meta_i)
        tabs.append(table_i)
        hists.append(hist_i)

    if rest_inputs_hz:
        rest_ref = rest_inputs_hz[0]
        if not np.allclose(np.asarray(rest_inputs_hz, dtype=float), rest_ref, rtol=0.0, atol=1.0e-6):
            raise ValueError("All coadd inputs must share the same RESTFRQ/RESTFREQ.")

    table_all = _df_to_native_endian(pd.concat(tabs, axis=0, ignore_index=True))
    meta_ref = metas[0]
    table_all = ensure_tempscal_column(table_all, default=str(meta_ref.get('TEMPSCAL', 'TA*')))
    table_all = ensure_beameff_column(table_all, default=meta_ref.get('BEAMEFF', float('nan')))
    beameff_vec = beameff_array(table_all, meta_ref, len(table_all), default=float('nan'))
    tempscal_vec = tempscal_array(table_all, meta_ref, len(table_all), default=str(meta_ref.get('TEMPSCAL', 'TA*')))
    
    if rows is not None and exclude_rows is not None:
        raise ValueError("Cannot specify both 'rows' and 'exclude_rows'.")

    total_len = len(table_all)
    final_idxs = None
    
    if rows is not None:
        final_idxs = _parse_row_selector(rows, total_len)
    elif exclude_rows is not None:
        exclude_idxs = _parse_row_selector(exclude_rows, total_len)
        final_idxs = np.setdiff1d(np.arange(total_len), exclude_idxs)
        
    if final_idxs is not None:
        final_idxs = np.unique(final_idxs)
        if len(final_idxs) == 0:
             raise ValueError("Selector resulted in empty coadd input.")

        table_all = table_all.iloc[final_idxs].reset_index(drop=True)
        data_list = [data_list[i] for i in final_idxs]
        beameff_vec = beameff_vec[final_idxs]
        tempscal_vec = tempscal_vec[final_idxs]

    sc_all = Scantable(meta=meta_ref, data=data_list, table=table_all)
    
    axis_type = str(axis_type).strip().lower()
    if axis_type not in ("freq", "vel"):
        raise ValueError("axis_type must be 'freq' or 'vel'.")

    std = Standardizer(sc_all, v_corr_col=v_corr_col)
    full_matrix, v_tgt = std.get_matrix(dv=dv, vmin=vmin, vmax=vmax)
    n_tgt = len(v_tgt)
    dv_val = float(v_tgt[1] - v_tgt[0]) if n_tgt > 1 else 0.0
    

    if verbose:
        print("v_tgt[0:3] =", v_tgt[:3])
        print("v_tgt[-3:] =", v_tgt[-3:])
        dv_median = np.nanmedian(np.diff(v_tgt)) if len(v_tgt) > 1 else np.nan
        print("dv(median) =", dv_median, " range =", (np.nanmin(v_tgt), np.nanmax(v_tgt)))
        if len(v_tgt) > 2:
            print("nonlinear check (max|d-dmed|) =", np.nanmax(np.abs(np.diff(v_tgt) - dv_median)))

    
    baseline_windows_parsed = parse_windows(baseline_vwin) if baseline_vwin else []
    rms_windows_parsed = parse_windows(rms_vwin) if rms_vwin else []
    line_windows_parsed = parse_windows(line_vwin) if line_vwin else []

    if post_baseline_mode_norm == "inherit_all":
        if not baseline_windows_parsed:
            raise ValueError(
                "post_baseline_mode='inherit_all' requires baseline_vwin to be set, "
                "because there is no pre-coadd baseline configuration to inherit."
            )
        # inherit_all provides defaults for post-baseline parameters, but explicit
        # post_baseline_vwin still overrides them. In particular, post_baseline_vwin=None
        # must disable post-baseline explicitly.
        if post_baseline_vwin is None:
            post_baseline_windows_parsed = []
        elif post_baseline_vwin == "inherit":
            post_baseline_windows_parsed = list(baseline_windows_parsed)
        else:
            post_baseline_windows_parsed = parse_windows(post_baseline_vwin) if post_baseline_vwin else []
    elif post_baseline_vwin == "inherit":
        post_baseline_windows_parsed = list(baseline_windows_parsed)
    else:
        post_baseline_windows_parsed = parse_windows(post_baseline_vwin) if post_baseline_vwin else []
    
    if line_windows_parsed:
        if baseline_windows_parsed:
            baseline_windows_parsed = subtract_windows(baseline_windows_parsed, line_windows_parsed)
        if rms_windows_parsed:
            rms_windows_parsed = subtract_windows(rms_windows_parsed, line_windows_parsed)
        if post_baseline_windows_parsed:
            post_baseline_windows_parsed = subtract_windows(post_baseline_windows_parsed, line_windows_parsed)

    def _inherit_or_cast_int(value, fallback: int) -> int:
        if isinstance(value, str) and value.strip().lower() == "inherit":
            return int(fallback)
        return int(value)

    def _inherit_or_cast_float(value, fallback: float) -> float:
        if isinstance(value, str) and value.strip().lower() == "inherit":
            return float(fallback)
        return float(value)

    if post_baseline_mode_norm == "inherit_all":
        resolved_post_poly = int(baseline_poly)
        resolved_post_iter_max = int(baseline_iter_max)
        resolved_post_iter_sigma = float(baseline_iter_sigma)
        if not (isinstance(post_baseline_poly, str) and post_baseline_poly.strip().lower() == "inherit"):
            resolved_post_poly = int(post_baseline_poly)
        if not (isinstance(post_baseline_iter_max, str) and post_baseline_iter_max.strip().lower() == "inherit"):
            resolved_post_iter_max = int(post_baseline_iter_max)
        if not (isinstance(post_baseline_iter_sigma, str) and post_baseline_iter_sigma.strip().lower() == "inherit"):
            resolved_post_iter_sigma = float(post_baseline_iter_sigma)
    else:
        resolved_post_poly = _inherit_or_cast_int(post_baseline_poly, int(baseline_poly))
        resolved_post_iter_max = _inherit_or_cast_int(post_baseline_iter_max, int(baseline_iter_max))
        resolved_post_iter_sigma = _inherit_or_cast_float(post_baseline_iter_sigma, float(baseline_iter_sigma))

    use_qc = coadd_qc is not None
    qc_config: dict = {}
    if use_qc:
        if not rms_windows_parsed:
            raise ValueError("QC Mode (--coadd-qc) requires --rms-vwin to be set.")
        qc_config = _parse_qc_spec(coadd_qc)
        if verbose:
            print(f"QC Mode Enabled: {qc_config}")

    if pos_col not in table_all.columns:
        # pos_tol_arcsec が指定されていない場合のみ、全体を1つのグループにする
        if pos_tol_arcsec is None:
            table_all["_dummy_group"] = 0
            pos_col = "_dummy_group"

    groups = _group_indices(table_all, group_mode, pos_col, pos_tol_arcsec)
    
    # [ADDED] 重み付け方針の判定と標準出力
    weight_strategy = "Uniform (weight=1.0)"
    if use_qc:
        weight_strategy = f"QC statistical weighting ({coadd_qc})"
    elif is_rms_mode:
        if baseline_windows_parsed:
            weight_strategy = "1/RMS^2 (calculated from baseline_vwin)"
        elif rms_windows_parsed:
            weight_strategy = "1/RMS^2 (calculated from rms_vwin)"
        else:
            weight_strategy = "1/RMS^2 (from header BSL_RMS)"

    if verbose:
        print(f"--- Coadd Configuration ---")
        print(f"Weighting Strategy : {weight_strategy}")
        if baseline_windows_parsed and not use_qc:
            print(f"Baseline Subtraction: Yes (poly_order={baseline_poly})")
        else:
            print(f"Baseline Subtraction: No")
        print(f"Total Groups       : {len(groups)}")
        print(f"---------------------------")
        
    out_rows: list[dict] = []
    out_specs: list[np.ndarray] = []
    
    mixed_groups_auto: List[dict] = []
    mixed_groups_never: List[dict] = []

    for gid, idxs in groups.items():
        idxs = np.asarray(idxs, int)
        rep = int(idxs[0])
        if "SSYSOBS" in table_all.columns:
            grp_vals = pd.Series(table_all.iloc[idxs]["SSYSOBS"], dtype=object).dropna().astype(str).str.strip()
            grp_vals = grp_vals[grp_vals != ""]
            if len(pd.unique(grp_vals)) > 1:
                raise ValueError(
                    f"Group {gid} mixes multiple SSYSOBS values {list(pd.unique(grp_vals))}; cannot assign a unique observed frame to one coadded row."
                )
        specs_arr = full_matrix[idxs]
        
        final_grp_rms = np.nan
        
        if use_qc:
            # [MODIFIED] QC用: 指定された rms_poly で事前にフィットして差し引く
            qc_specs_arr = specs_arr.copy()
            for i in range(len(qc_specs_arr)):
                _, base_qc, _ = fit_polynomial_baseline(
                    v_tgt, qc_specs_arr[i],
                    base_windows=rms_windows_parsed,
                    poly_order=int(rms_poly), iter_max=0
                )
                qc_specs_arr[i] -= base_qc

            weights, stats = _compute_qc_weights(
                qc_specs_arr, v_tgt, rms_windows_parsed, qc_config, dv_val
            )
            
            nonzero_w = weights[weights > 0]
            if len(nonzero_w) > 0:
                median_w = float(np.median(nonzero_w))
                has_data = np.any(np.isfinite(specs_arr), axis=1)
                rescue_mask = (weights <= 0) & has_data
                if np.any(rescue_mask):
                    n_rescued = int(np.count_nonzero(rescue_mask))
                    weights[rescue_mask] = median_w
                    print(f"Warning: {n_rescued} spectra in Group {gid} had valid data but 0 weight "
                          f"(likely outside QC RMS window). Imputed median weight ({median_w:.3g}) "
                          f"to prevent data loss.")
                    stats["n_keep"] = int(np.count_nonzero(weights > 0))

            if verbose and (len(groups) < 20 or stats['n_keep'] == 0):
                print(f"Grp {gid}: N_in={stats['n_dump']}, N_used={stats['n_keep']}, "
                      f"q(med)={stats['q_median']:.2f}, Strategy=QC")
                      
            clip_k = float(qc_config.get("clip_k", 0))
            clip_iter = int(qc_config.get("clip_iter", 0))
            
            out_spec, _ = _robust_channel_coadd(specs_arr, weights, clip_k, clip_iter)
            n_used = int(stats["n_keep"])
            mode_label = f"QC:{coadd_qc}"
            
        else:
            # === Normal / Legacy Logic ===
            n_rows_local = len(specs_arr)
            weights = np.ones(n_rows_local, dtype=float)
            processed_specs = []

            use_header_rms = False
            grp_rms = None
            if is_rms_mode and not baseline_windows_parsed and not rms_windows_parsed:
                if "BSL_RMS" not in table_all.columns:
                    raise ValueError("mode='rms' requested but no windows specified and 'BSL_RMS' column not found in input data.")
                grp_rms = table_all["BSL_RMS"].to_numpy()[idxs]
                if np.any(np.isnan(grp_rms)) or np.any(grp_rms <= 0):
                    raise ValueError(f"mode='rms' requested but Group {gid} has missing or invalid (<=0) 'BSL_RMS' values.")
                use_header_rms = True

            for i in range(n_rows_local):
                y_in = specs_arr[i]
                y_for_coadd = y_in.copy()
                w = 1.0

                if baseline_windows_parsed:
                    # 1. 個別フィット・差し引き
                    _, base, info = fit_polynomial_baseline(
                        v_tgt, y_in,
                        base_windows=baseline_windows_parsed,
                        poly_order=int(baseline_poly),
                        iter_max=int(baseline_iter_max),
                        iter_sigma=float(baseline_iter_sigma)
                    )
                    y_for_coadd = y_in - base
                    
                    if is_rms_mode:
                        # [MODIFIED] 既に差し引いた残差に対して、rms_bin を適用して重みを計算
                        resid_valid = y_for_coadd[info.get("mask", np.zeros_like(y_for_coadd, dtype=bool))]
                        if int(rms_bin) > 1:
                            resid_valid = _fast_nan_boxcar_2d(resid_valid[None, :], int(rms_bin))[0]
                        resid_fin = resid_valid[np.isfinite(resid_valid)]
                        
                        if resid_fin.size >= 2:
                            rms_val = float(np.std(resid_fin, ddof=1))
                            if np.isfinite(rms_val) and rms_val > 0:
                                w = 1.0 / (rms_val**2)
                            else:
                                w = 0.0
                        else:
                            w = 0.0
                
                elif rms_windows_parsed:
                    # 2. RMS評価用のフィット (データからは引かない)
                    if is_rms_mode:
                        w, _ = _evaluate_rms_for_weight(
                            v_tgt, y_in, rms_windows_parsed, int(rms_poly), int(rms_bin)
                        )
                
                elif use_header_rms:
                    # 3. ヘッダのRMSを利用
                    w = 1.0 / (grp_rms[i]**2)

                else:
                    # 4. Uniform
                    w = 1.0

                weights[i] = w
                processed_specs.append(y_for_coadd)

            if (rms_windows_parsed or baseline_windows_parsed) and is_rms_mode:
                invalid_weight_mask = np.array([
                    (weights[i] <= 0) and np.any(np.isfinite(processed_specs[i]))
                    for i in range(n_rows_local)
                ], dtype=bool)
                invalid_count = int(np.count_nonzero(invalid_weight_mask))
                if invalid_count > 0:
                    if weight_zero_policy == "error":
                        raise ValueError(
                            f"Group {gid} contains {invalid_count} spectra with valid data but non-positive RMS weight. "
                            "Use weight_zero_policy='drop' to discard them or 'impute_median' to force a replacement weight."
                        )
                    elif weight_zero_policy == "impute_median":
                        valid_weights = weights[weights > 0]
                        if len(valid_weights) == 0:
                            raise ValueError(
                                f"Group {gid} has no positive RMS-derived weights, so impute_median cannot be applied."
                            )
                        typical_weight = float(np.median(valid_weights))
                        for i in range(n_rows_local):
                            if invalid_weight_mask[i]:
                                weights[i] = typical_weight
                        print(f"Warning: {invalid_count} spectra in Group {gid} had valid data but 0 weight. "
                              f"Imputed median weight ({typical_weight:.3g}) as requested by weight_zero_policy='impute_median'.")
                    elif weight_zero_policy == "drop":
                        print(f"Warning: Dropping {invalid_count} spectra in Group {gid} because weight_zero_policy='drop'.")

            # Weighted Sum
            sum_spec = np.zeros(n_tgt, dtype=float)
            sum_w_chan = np.zeros(n_tgt, dtype=float)
            n_used = 0

            for i in range(n_rows_local):
                if weights[i] > 0:
                    y_in = processed_specs[i]
                    valid = np.isfinite(y_in)
                    if np.any(valid):
                        sum_spec[valid] += weights[i] * y_in[valid]
                        sum_w_chan[valid] += weights[i]
                        n_used += 1

            out_spec = np.full(n_tgt, np.nan)
            m = sum_w_chan > 0
            out_spec[m] = sum_spec[m] / sum_w_chan[m]
            mode_label = "rms" if is_rms_mode else "uniform"
            
            # [ADDED] 各グループの処理結果を出力
            if verbose and len(groups) < 20:
                print(f"Grp {gid}: N_in={n_rows_local}, N_used={n_used}, Strategy={weight_strategy}")

        # --------------------------------------------------------------
        group_beameff = beameff_vec[np.array(idxs, dtype=int)]
        group_mixed = is_beameff_mixed(group_beameff, tol=float(beameff_tol))
        n_rows_local = len(specs_arr)
        
        if group_mixed:
            stats = dict(
                gid=int(gid),
                n_spec=int(len(idxs)),
                beameff_min=float(np.nanmin(group_beameff)),
                beameff_max=float(np.nanmax(group_beameff)),
            )
            if normalize_if_mixed == "auto":
                mixed_groups_auto.append(stats)
            elif normalize_if_mixed == "never":
                mixed_groups_never.append(stats)
                warnings.warn(
                    "BEAMEFF is mixed within a coadd group but normalize_if_mixed='never' was requested. "
                    "Proceeding can lead to inconsistent scaling across spectra.",
                    RuntimeWarning,
                )

        do_tr_norm = (out_scale == "TR*") or (group_mixed and normalize_if_mixed == "auto")
        
        if do_tr_norm:
            b_final = 1.0
            scale_label = "TR*"
        else:
            b_final = float(representative_beameff(group_beameff))
            scale_label = "TA*"

        if do_tr_norm:
            sum_num = np.zeros(n_tgt, dtype=float)
            sum_den = np.zeros(n_tgt, dtype=float)
            
            valid_sum = False
            for i in range(n_rows_local):
                w_i = weights[i]
                if w_i > 0:
                    y_in = processed_specs[i] if not use_qc else specs_arr[i]
                    valid = np.isfinite(y_in)
                    if np.any(valid):
                        global_idx = idxs[i]
                        input_scale = tempscal_vec[global_idx]
                        bi = float(group_beameff[i])
                        
                        if input_scale == "TR*":
                            tr_val = y_in[valid]
                        else:
                            if bi > 0:
                                tr_val = y_in[valid] / bi
                            else:
                                tr_val = np.nan

                        if np.all(np.isfinite(tr_val)):
                            sum_num[valid] += w_i * tr_val
                            sum_den[valid] += w_i
                            valid_sum = True
            
            if valid_sum:
                tr_spec = np.full(n_tgt, np.nan)
                m2 = sum_den > 0
                tr_spec[m2] = sum_num[m2] / sum_den[m2]
                out_spec = tr_spec
                
            mode_label = f"{mode_label}+TRout"

        do_post_baseline = bool(post_baseline_windows_parsed) and (not use_qc)
        post_bsl_coeff = None
        post_bsl_nused = 0
        post_bsl_winf = ""
        if do_post_baseline:
            post_bsl_coeff, base_post, info_post = fit_polynomial_baseline(
                v_tgt, out_spec,
                base_windows=post_baseline_windows_parsed,
                poly_order=int(resolved_post_poly),
                iter_max=int(resolved_post_iter_max),
                iter_sigma=float(resolved_post_iter_sigma)
            )
            out_spec = out_spec - base_post
            final_grp_rms = info_post.get("std", np.nan)
            post_bsl_nused = int(np.count_nonzero(info_post.get("mask", np.zeros_like(out_spec, dtype=bool))))
            if isinstance(post_baseline_vwin, list):
                post_bsl_winf = ";".join(post_baseline_vwin)
            else:
                post_bsl_winf = ";".join(baseline_vwin) if baseline_vwin else ""

        out_specs.append(out_spec)

        row = table_all.iloc[rep].to_dict()
        row = _purge_bsl_metadata_in_row(row)
        if pos_col in row: row[pos_col] = int(row[pos_col])
        
        row["TEMPSCAL"] = scale_label
        row["BEAMEFF"] = float(b_final) if np.isfinite(b_final) else np.nan
        row["N_IN"] = int(len(idxs))
        row["N_USED"] = int(n_used)
        row["WEIGHT_MODE"] = mode_label

        row["WGT_MODE"] = "inverse_variance" if is_rms_mode else "uniform"
        if baseline_windows_parsed:
            row["WGT_SRC"] = "baseline_vwin" if is_rms_mode else "uniform"
            row["WGT_VWIN"] = ";".join(baseline_vwin) if baseline_vwin else ""
            row["WGT_POLY"] = int(baseline_poly)
        elif rms_windows_parsed:
            row["WGT_SRC"] = "rms_vwin"
            row["WGT_VWIN"] = ";".join(rms_vwin) if rms_vwin else ""
            row["WGT_POLY"] = int(rms_poly)
        elif use_header_rms:
            row["WGT_SRC"] = "input_bsl_rms"
            row["WGT_VWIN"] = ""
            row["WGT_POLY"] = 0
        else:
            row["WGT_SRC"] = "uniform"
            row["WGT_VWIN"] = ""
            row["WGT_POLY"] = 0
        row["WGT_BIN"] = int(rms_bin) if is_rms_mode else 1
        row["WGT_STAT"] = "std" if is_rms_mode else ""
        row["WGT_INPUT_SUB"] = bool(baseline_windows_parsed)
        row["WGT_ZERO_POLICY"] = str(weight_zero_policy) if is_rms_mode else ""
        
        if do_post_baseline:
            row["BSL_DONE"] = True
            row["BSL_APPLIED"] = True
            row["BSL_STAGE"] = "post_coadd"
            row["BSL_POLY"] = int(resolved_post_poly)
            row["BSL_WINF"] = post_bsl_winf
            row["BSL_RMS"] = float(final_grp_rms)
            row["BSL_STAT"] = "std"
            row["BSL_NUSED"] = int(post_bsl_nused)
            row["BSL_COEF"] = None if post_bsl_coeff is None else np.asarray(post_bsl_coeff, dtype=float)
            row["BSL_SCALE"] = scale_label
            
        out_rows.append(row)

    if not out_specs: raise ValueError("No output spectra produced.")

    out_data = np.vstack(out_specs)
    out_table = pd.DataFrame(out_rows)

    if "SSYSOBS" not in out_table.columns:
        out_table["SSYSOBS"] = _resolve_ssysobs_series(out_table, meta_ref).to_numpy(dtype=object)
    ssysobs_out = _unique_nonempty_upper(out_table["SSYSOBS"])

    sig_cols = [c for c in out_table.columns if str(c).upper().startswith("SIG_")]
    out_table = out_table.drop(columns=sig_cols + ["_local_idx", "_dummy_group", "_time_is_dummy", "_INPUT_ID", "INPUT_ID", "input_id"], errors="ignore")

    out_table = ensure_tempscal_column(out_table, default="TA*")
    out_table = ensure_beameff_column(out_table)
    
    if "VELOSYS" in out_table.columns:
        out_table["VELOSYS_OBS"] = pd.to_numeric(out_table["VELOSYS"], errors="coerce")
    out_table["VELOSYS"] = 0.0
    if "VFRAME" in out_table.columns:
        out_table["VFRAME_OBS"] = pd.to_numeric(out_table["VFRAME"], errors="coerce")
        out_table["VFRAME"] = 0.0
    elif str(v_corr_col).strip().upper() == "VFRAME":
        out_table["VFRAME"] = 0.0
    if v_corr_col in out_table.columns and v_corr_col not in ("VELOSYS", "VFRAME"):
        out_table[f"{v_corr_col}_OBS"] = pd.to_numeric(out_table[v_corr_col], errors="coerce")
        out_table[v_corr_col] = 0.0

    rest0 = _get_rest_hz(meta_ref, table_all, override_hz=rest_freq)
    if not np.isfinite(rest0) or rest0 <= 0:
        raise ValueError("RESTFRQ/RESTFREQ is required for coadd output.")

    meta_out = dict(meta_ref)

    meta_out["TEMPSCAL"] = str(out_table["TEMPSCAL"].iloc[0]) if not out_table.empty else "TA*"
    meta_out["BEAMEFF"] = float(out_table["BEAMEFF"].iloc[0]) if not out_table.empty else 1.0
    meta_out.pop("VELREF", None)

    if axis_type == "vel":
        crval1 = float(v_tgt[0] * 1000.0) if n_tgt > 0 else 0.0
        cdelt1 = float(dv_val * 1000.0)
        out_table["CTYPE1"] = "VRAD"
        out_table["CUNIT1"] = "m/s"
        out_table["CRVAL1"] = crval1
        out_table["CDELT1"] = cdelt1
        out_table["CRPIX1"] = 1.0
        out_table["SPECSYS"] = "LSRK"
        out_table["RESTFRQ"] = float(rest0)
        out_table["RESTFREQ"] = float(rest0)
        meta_out.update(dict(
            CTYPE1="VRAD", CUNIT1="m/s",
            CRVAL1=crval1, CDELT1=cdelt1, CRPIX1=float(1.0),
            SPECSYS="LSRK", VELOSYS=0.0,
            TIMESYS="UTC", NCHAN=int(n_tgt),
            RESTFRQ=float(rest0), RESTFREQ=float(rest0),
            VELDEF="RADIO", BUNIT="K",
        ))
        if ssysobs_out:
            meta_out["SSYSOBS"] = ssysobs_out
        else:
            meta_out.pop("SSYSOBS", None)
    else:
        c_kms = 299792.458
        freq_axis = rest0 * (1.0 - v_tgt / c_kms)
        crval1 = float(freq_axis[0]) if n_tgt > 0 else float(rest0)
        cdelt1 = float(freq_axis[1] - freq_axis[0]) if n_tgt > 1 else 0.0
        out_table["CTYPE1"] = "FREQ"
        out_table["CUNIT1"] = "Hz"
        out_table["CRVAL1"] = crval1
        out_table["CDELT1"] = cdelt1
        out_table["CRPIX1"] = 1.0
        out_table["SPECSYS"] = "LSRK"
        out_table["RESTFRQ"] = float(rest0)
        out_table["RESTFREQ"] = float(rest0)
        meta_out.update(dict(
            CTYPE1="FREQ", CUNIT1="Hz",
            CRVAL1=crval1,
            CDELT1=cdelt1,
            CRPIX1=float(1.0),
            SPECSYS="LSRK", VELOSYS=0.0,
            TIMESYS="UTC", NCHAN=int(n_tgt),
            RESTFRQ=float(rest0), RESTFREQ=float(rest0),
            VELDEF="RADIO", BUNIT="K",
        ))
        if ssysobs_out:
            meta_out["SSYSOBS"] = ssysobs_out
        else:
            meta_out.pop("SSYSOBS", None)
    if "NAXIS1" in meta_out: del meta_out["NAXIS1"]

    history = dict(
        stage="coadd_fits_vlsrk_qc", 
        out_scale=str(out_scale),
        normalize_if_mixed=str(normalize_if_mixed),
        beameff_tol=float(beameff_tol),
        sigma_scale=str(sigma_scale),
        created_at_utc=datetime.datetime.utcnow().isoformat() + "Z",
        mode=str(mode),
        group_mode=str(group_mode),
        rows=str(rows),
        exclude_rows=str(exclude_rows),
        vcorr_chunk_sec=str(vcorr_chunk_sec),
        qc_config=qc_config if use_qc else None,
        # [MODIFIED] rmsの評価に使った次数とビニング幅も追加
        rms=dict(windows=rms_windows_parsed, poly=int(rms_poly), bin=int(rms_bin)),
        baseline=dict(windows=baseline_windows_parsed, poly=int(baseline_poly)),
        post_baseline=dict(
            mode=post_baseline_mode_norm,
            windows=post_baseline_windows_parsed,
            poly=int(resolved_post_poly),
            iter_max=int(resolved_post_iter_max),
            iter_sigma=float(resolved_post_iter_sigma),
        ),
        weight_zero_policy=str(weight_zero_policy),
    )
    
    if mixed_groups_auto:
        history = append_scale_history(history, {
            "stage": "coadd",
            "event": "beameff_mixed_detected",
            "normalize_if_mixed": "auto",
            "groups": mixed_groups_auto[:8],
            "formula_tr": "Tr* = Ta* / BEAMEFF",
            "note": "Mixed BEAMEFF groups were normalized to TR* (output BEAMEFF=1).",
        })
    if mixed_groups_never:
        history = append_scale_history(history, {
            "stage": "coadd",
            "event": "beameff_mixed_detected",
            "normalize_if_mixed": "never",
            "groups": mixed_groups_never[:8],
            "note": "Mixed BEAMEFF was present but normalization was explicitly disabled (danger mode).",
        })

    res = Scantable(meta=meta_out, data=out_data, table=out_table, history=history)

    if output_path:
        write_scantable(output_path, res, overwrite=overwrite, tempscal=meta_out["TEMPSCAL"], data_scale=meta_out["TEMPSCAL"])

    return res
