#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
necst_v4_sdfits_converter_v1_40.py

NECST v4 (RawData folder: necstdb + nercst) の分光データを SDFITS に変換します。

v1_40:
- v1_38 の最終点検で見つかった OBSFREQ fallback / heterodyne context の不整合を修正
- refraction 用代表周波数に書き出し軸の中心周波数を使うよう調整
- そのほか v1_38 の機能を維持
- HISTORY の旧版表記を整理
- TOML ベースの multi-spectrometer 設定を導入
- stream ごとに beam / IF / pol / frontend / backend / sampler / LO chain / 周波数軸を定義可能
- 周波数軸は stream ごとに導出し、row ごとに WCS metadata を writer に渡す
- Az/El 系の意味を整理
  - AZIMUTH/ELEVATIO : row が最終的に表す beam-center Az/El
  - AZ_CMD/EL_CMD    : altaz / cmd を t_spec に内挿した指令値
  - CORR_AZ/CORR_EL  : dlon / dlat
- boresight (= encoder - correction) と beam-center を分離
- beam offset / beam rotation を converter 側で計算
- RA/DEC fallback は beam-center Az/El から計算
- `OBSFREQ` は周波数軸そのものではなく heterodyne metadata として扱う
  - 未指定時は `RESTFREQ` から補完可能
- 既定の出力配置は merged_time
  - 複数 stream の row を 1 つの SINGLE DISH テーブルへ集約し、時刻順に並べる
- converter 時点での channel slice をサポート
- global.db_namespace / global.telescope / spectrometers.db_stream_name / spectrometers.db_table_name を実際の DB 名解決に使用
- --db-namespace を追加し、CLI 明示指定 > TOML > 既定値 の優先順位で適用
  - `channel_slice='[chmin,chmax)'` 等で stream ごとに切り出し可能
  - 切り出し後の `NCHAN/CRVAL1/CRPIX1` は converter 側で更新

互換性
- `--spectrometer-config` 未指定時は legacy single-stream CLI でも動作
- writer 更新前でも動くよう、現 writer が受けられる範囲で metadata を書き出す
"""

import argparse
from dataclasses import dataclass, field, replace
import pathlib
import sys
import re
import math
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import tomllib  # Python 3.11+
except Exception:  # pragma: no cover
    tomllib = None

import numpy as np
import pandas as pd

import necstdb

# NOTE: nercst に依存せず、necstdb から直接読み出します（壊れた末尾レコードを自前で切り捨て可能）。

from astropy.time import Time
from astropy.coordinates import AltAz, SkyCoord
from astropy import units as u

from sd_radio_spectral_fits import (
    SDRadioSpectralSDFITSWriter,
    Site, DatasetInfo, SpectralAxisUniform, Efficiency,
)

# -----------------------------------------------------------------------------
# Utilities（小さく・依存を減らす）
# -----------------------------------------------------------------------------
def _decode_label(v):
    if isinstance(v, (bytes, bytearray)):
        try:
            s = v.decode(errors="ignore")
        except Exception:
            s = str(v)
    else:
        s = str(v)
    return s.strip().upper()


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


def _sanitize_for_filename(s, maxlen=80):
    """
    Make a safe filename fragment from an arbitrary string.
    """
    s = str(s) if s is not None else ""
    s = s.strip()
    if not s:
        return "empty"
    # replace path separators and spaces
    s = s.replace("/", "_").replace("\\", "_").replace(" ", "_")
    # keep a conservative charset
    s = re.sub(r"[^0-9A-Za-z._-]+", "_", s)
    s = re.sub(r"_+", "_", s)
    s = s.strip("._-")
    if not s:
        s = "empty"
    return s[:int(maxlen)]
def _norm_name(x):
    try:
        return str(x).strip()
    except Exception:
        return ""

def _pick_field_name(names, preferred, candidates):
    """
    names: list[str]
    preferred: optional string (user-specified)
    candidates: list[str] fallbacks
    """
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
    """
    t_ref: reference unix seconds (around 1e9)
    t_other: time array that may be in ms/us/ns.
    Returns (t_other_seconds, scale_applied)
      scale_applied: 1 (no change), 1e-3 (ms->s), 1e-6 (us->s), 1e-9 (ns->s), etc.
    """
    tr = np.asarray(t_ref, dtype=float)
    to = np.asarray(t_other, dtype=float)

    mr = float(np.nanmedian(tr)) if tr.size else np.nan
    mo = float(np.nanmedian(to)) if to.size else np.nan
    if (not np.isfinite(mr)) or (not np.isfinite(mo)) or mr <= 0 or mo <= 0:
        return to, 1.0

    ratio = mo / mr

    # other is larger -> likely ms/us/ns
    if 5e2 < ratio < 2e3:
        return to / 1e3, 1e-3
    if 5e5 < ratio < 2e6:
        return to / 1e6, 1e-6
    if 5e8 < ratio < 2e9:
        return to / 1e9, 1e-9

    # other is smaller -> likely ref is ms/us/ns (rare)
    if 5e2 < (1.0/ratio) < 2e3:
        return to * 1e3, 1e3
    if 5e5 < (1.0/ratio) < 2e6:
        return to * 1e6, 1e6
    if 5e8 < (1.0/ratio) < 2e9:
        return to * 1e9, 1e9

    return to, 1.0

def _safe_read_table(db, name):
    try:
        return db.open_table(name).read(astype="pandas")
    except Exception as e:
        print("[warn] cannot read table {}: {}".format(name, e))
        return None

def _read_table_raw_bytes(table):
    """
    necstdb の read(astype="raw"/"buffer") で生バッファを取得する。
    structured array への変換は行わないので、末尾が壊れていても読み出せることが多い。
    """
    for k in ("raw", "buffer"):
        try:
            b = table.read(astype=k)
            if isinstance(b, (bytes, bytearray)):
                return bytes(b)
        except Exception:
            pass
    raise RuntimeError("cannot read raw/buffer from table")

def _get_table_dtype(table):
    """
    necstdb table オブジェクトから dtype を取得（実装差に備えて複数候補）。
    """
    for attr in ("dtype", "_dtype", "data_dtype", "_data_dtype"):
        try:
            dt = getattr(table, attr)
            if dt is not None:
                return dt
        except Exception:
            pass
    # 最後の手段: 少数レコードで read(astype="array") が通る場合はそこから推定
    try:
        arr = table.read(num=1, astype="array")
        return arr.dtype
    except Exception:
        return None

def _read_structured_array_tolerant(db, table_name):
    """
    まず通常の read(astype="array") を試し、失敗したら raw/buffer から復旧する。
    復旧では「dtype.itemsize の倍数になるよう末尾バイトを切り捨て」る。

    """
    tbl = db.open_table(table_name)
    try:
        return tbl.read(astype="array")
    except Exception as e:
        msg = str(e)
        raw = _read_table_raw_bytes(tbl)
        dt = _get_table_dtype(tbl)
        if dt is None:
            raise RuntimeError("cannot determine dtype for table {} (original error: {})".format(table_name, msg))
        item = int(dt.itemsize)
        if item <= 0:
            raise RuntimeError("invalid dtype itemsize for table {}: {}".format(table_name, item))
        nrec = len(raw) // item
        if nrec <= 0:
            raise RuntimeError("raw buffer too small for table {}: len={} itemsize={}".format(table_name, len(raw), item))
        raw2 = raw[: nrec * item]  # drop incomplete tail
        try:
            return np.frombuffer(raw2, dtype=dt)
        except Exception as e2:
            raise RuntimeError("tolerant frombuffer failed for table {}: {} (orig: {})".format(table_name, e2, msg))

def _structured_to_dataframe(arr):
    """
    structured array -> pandas.DataFrame（配列フィールドは object になる可能性あり）
    """
    try:
        return pd.DataFrame.from_records(arr)
    except Exception:
        d = {name: arr[name] for name in (arr.dtype.names or [])}
        return pd.DataFrame(d)

def _extract_spectral_from_structured(arr, nchan):
    """
    spectral structured array から t_spec と spec2d を取り出す。
    - time basis is decided once from the first spectral timestamp row
    - accepted timestamp suffixes: UTC / GPS / TAI
    - PC / unknown / NaN timestamp falls back to numeric time field
    - spec field: prefer data > spectrum > spec
    """
    names = list(arr.dtype.names or [])
    s_field = None
    for k in ("data", "spectrum", "spec"):
        if k in names:
            s_field = k
            break
    if s_field is None:
        raise RuntimeError("no spectrum-like field in spectral table. available={}".format(names))

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
        raise RuntimeError("spectrum array must be 2D. got {}".format(spec2d.shape))
    if int(spec2d.shape[1]) != int(nchan):
        raise RuntimeError("NCHAN mismatch in spectral table: got {}, expected {}".format(spec2d.shape[1], nchan))

    return np.asarray(t_spec, dtype=float), spec2d, time_meta, s_field


def _interp_lin(t_src, y_src, t_dst, extrap="nan"):
    t_src = np.asarray(t_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)

    m = np.isfinite(t_src) & np.isfinite(y_src)
    if np.count_nonzero(m) < 2:
        return np.full_like(t_dst, np.nan, dtype=float)

    t = t_src[m]
    y = y_src[m]
    o = np.argsort(t)
    t = t[o]
    y = y[o]

    mode = str(extrap).lower().strip()
    if mode == "hold":
        return np.interp(t_dst, t, y, left=float(y[0]), right=float(y[-1]))
    return np.interp(t_dst, t, y, left=np.nan, right=np.nan)

def _interp_az_deg(t_src, az_deg, t_dst, extrap="nan"):
    """
    Circular interpolation for Az(deg) with unwrap.
    """
    t_src = np.asarray(t_src, dtype=float)
    az_deg = np.asarray(az_deg, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)

    m = np.isfinite(t_src) & np.isfinite(az_deg)
    if np.count_nonzero(m) < 2:
        return np.full_like(t_dst, np.nan, dtype=float)

    t = t_src[m]
    az = np.deg2rad(az_deg[m])
    az_unw = np.unwrap(az)

    o = np.argsort(t)
    t = t[o]
    az_unw = az_unw[o]

    mode = str(extrap).lower().strip()
    if mode == "hold":
        out = np.interp(t_dst, t, az_unw, left=float(az_unw[0]), right=float(az_unw[-1]))
    else:
        out = np.interp(t_dst, t, az_unw, left=np.nan, right=np.nan)

    return (np.rad2deg(out) % 360.0)


def _wrap180_deg(d_deg):
    """
    wrap angle difference (deg) into [-180, 180)
    """
    d = np.asarray(d_deg, dtype=float)
    return ((d + 180.0) % 360.0) - 180.0

def _estimate_exposure_s(t_unix_sorted, default=0.1):
    t = np.asarray(t_unix_sorted, dtype=float)
    if t.size < 2:
        return float(default)
    dt = np.diff(t)
    med = float(np.nanmedian(dt))
    if (not np.isfinite(med)) or (med <= 0):
        return float(default)
    return float(med)

def _to_degC_guess(temp_arr):
    """
    temp_arr が Kelvin なら degC に変換（互換優先の簡易判定）
    """
    x = np.asarray(temp_arr, dtype=float)
    if x.size == 0:
        return np.full_like(x, np.nan, dtype=float)
    med = float(np.nanmedian(x))
    if np.isfinite(med) and med > 120.0:
        return x - 273.15
    return x

def _to_hpa_guess(press_arr):
    """pressure が Pa の可能性を簡易判定し、hPa に揃える。"""
    x = np.asarray(press_arr, dtype=float)
    if x.size == 0:
        return np.full_like(x, np.nan, dtype=float)
    med = float(np.nanmedian(x))
    if np.isfinite(med) and med > 2000.0:
        return x / 100.0
    return x

def _to_rh01_guess(humid_pct_arr):
    """humidity[%] を 0..1 の relative humidity に揃える（範囲外は clip）。"""
    x = np.asarray(humid_pct_arr, dtype=float)
    if x.size == 0:
        return np.full_like(x, np.nan, dtype=float)
    rh = x / 100.0
    return np.clip(rh, 0.0, 1.0)

def _pressure_from_elev_hpa(elev_m):
    """高度[m]からの粗い気圧推定（標準大気の指数近似）。"""
    try:
        z = float(elev_m)
    except Exception:
        z = 0.0
    # scale height ~ 8434.5 m
    return 1013.25 * float(np.exp(-z / 8434.5))

def _try_radec_timeseries_from_wcs_table(db, wcs_table_name):
    """
    NECST の wcs テーブルが取れれば timestamp-index の ra/dec series を返す。
    形式が違えば None。
    """
    try:
        arr = db.open_table(wcs_table_name).read(astype="array")
        ts = pd.to_datetime(arr["timestamp"].astype(np.float64), unit="s")
        ra = np.array([r[1][1] for r in arr], dtype=np.float64)
        dec = np.array([r[1][2] for r in arr], dtype=np.float64)
        df = pd.DataFrame({"ra_deg": ra, "dec_deg": dec}, index=ts).sort_index()
        df.index.name = "timestamp"
        return df
    except Exception:
        return None

def _try_meteo_timeseries_from_weather_table(db, table_name, t_spec,
                                             preferred_time_col="time",
                                             extrap="hold"):
    """weather-ambient テーブルから pressure/temperature/humidity を t_spec に内挿して返す。

    注意:
      - recorded_time は使用しない（存在しても無視）
      - time 列が取得できない場合は None を返す
    戻り値:
      dict or None:
        {"press_hpa":..., "temp_c":..., "humid_pct":..., "table":..., "time_col":...}
    """
    try:
        arr = _read_structured_array_tolerant(db, str(table_name))
    except Exception as e:
        print("[warn] cannot read weather table {}: {}".format(table_name, e))
        return None

    names = list(arr.dtype.names or [])
    pref = str(preferred_time_col).strip()
    if pref.lower() == "recorded_time":
        print("[warn] weather_time_col=recorded_time is not supported; forcing to 'time'")
        pref = "time"

    t_name = _pick_field_name(names, pref, ["time", "timestamp"])
    if (t_name is None) or (str(t_name).strip().lower() == "recorded_time"):
        # recorded_time は使わない
        print("[warn] weather table has no usable time column (time/timestamp). recorded_time is ignored. available={}".format(names))
        return None

    # required fields
    if ("pressure" not in names) and ("temperature" not in names) and ("humidity" not in names):
        print("[warn] weather table has no pressure/temperature/humidity fields. available={}".format(names))
        return None

    t_w = _get_field(arr, t_name, float)
    t_w, s_w = _normalize_time_units(t_spec, t_w, "weather")
    if s_w != 1.0:
        print("[info] weather time scaled by {} to match unix seconds".format(s_w))

    press_hpa = None
    temp_c = None
    humid_pct = None

    if "pressure" in names:
        p = _get_field(arr, "pressure", float)
        p_t = _interp_lin(t_w, p, t_spec, extrap=extrap)
        press_hpa = _to_hpa_guess(p_t)

    if "temperature" in names:
        tc = _get_field(arr, "temperature", float)
        tc_t = _interp_lin(t_w, tc, t_spec, extrap=extrap)
        temp_c = _to_degC_guess(tc_t)

    if "humidity" in names:
        h = _get_field(arr, "humidity", float)
        h_t = _interp_lin(t_w, h, t_spec, extrap=extrap)
        med = float(np.nanmedian(h_t)) if h_t.size > 0 else np.nan
        humid_pct = h_t * 100.0 if (np.isfinite(med) and med <= 1.5) else h_t

    # If everything is NaN, treat as unusable
    ok_any = False
    for v in (press_hpa, temp_c, humid_pct):
        if v is not None and np.any(np.isfinite(np.asarray(v, dtype=float))):
            ok_any = True
            break
    if not ok_any:
        print("[warn] weather table fields exist but all are NaN after interpolation; ignoring weather meteo")
        return None

    return {
        "press_hpa": press_hpa,
        "temp_c": temp_c,
        "humid_pct": humid_pct,
        "table": str(table_name),
        "time_col": str(t_name),
    }

def _nearest_series(wcs_df, t_index, col, fallback):
    if wcs_df is None or wcs_df.empty:
        return np.full(len(t_index), float(fallback), dtype=float)

    left = pd.DataFrame({"timestamp": t_index})
    right = wcs_df.reset_index()[["timestamp", col]].sort_values("timestamp")

    try:
        merged = pd.merge_asof(left.sort_values("timestamp"), right, on="timestamp", direction="nearest")
        out = merged[col].to_numpy(dtype=float)
        if out.shape[0] != len(t_index):
            return np.full(len(t_index), float(fallback), dtype=float)
        return out
    except Exception:
        # slow fallback
        out = np.empty(len(t_index), dtype=float)
        for i, ti in enumerate(t_index):
            j = wcs_df.index.get_indexer([ti], method="nearest")[0]
            v = float(wcs_df.iloc[j][col])
            out[i] = v if np.isfinite(v) else float(fallback)
        return out

def _radec_from_azel(site, t_unix, az_deg, el_deg,
                     apply_refraction=True,
                     press_hpa=None,
                     temp_c=None,
                     rh01=None,
                     obswl_um=None):
    """Az/El(=AltAz) から ICRS RA/DEC を計算。

    apply_refraction=True の場合、press/temperature/RH/obswl を与えることで
    atmospheric refraction を考慮した変換を行う。
    """
    loc = site.to_earthlocation()
    t = Time(np.asarray(t_unix, dtype=float), format="unix", scale="utc")
    az = np.asarray(az_deg, dtype=float) * u.deg
    alt = np.asarray(el_deg, dtype=float) * u.deg

    if (not apply_refraction) or (press_hpa is None):
        frame = AltAz(obstime=t, location=loc, pressure=0 * u.hPa)
    else:
        p = np.asarray(press_hpa, dtype=float) * u.hPa
        tc = np.asarray(temp_c, dtype=float) * u.deg_C if temp_c is not None else 10.0 * u.deg_C
        rh = np.asarray(rh01, dtype=float) if rh01 is not None else 0.5
        # obswl は micron を想定（未指定なら 2600um 程度の mm 波を仮定）
        ow = (float(obswl_um) if (obswl_um is not None and np.isfinite(obswl_um)) else 2600.0) * u.micron
        frame = AltAz(obstime=t, location=loc, pressure=p, temperature=tc, relative_humidity=rh, obswl=ow)

    aa = SkyCoord(az=az, alt=alt, frame=frame)
    icrs = aa.icrs
    return np.asarray(icrs.ra.to_value(u.deg), dtype=float), np.asarray(icrs.dec.to_value(u.deg), dtype=float)

def _mode_from_pos(s):
    """Normalize position/OBSMODE label.
    ON/OFF/HOT/SKY はそのまま返す。空文字は UNKNOWN に落とす。
    その他の値は raw を保持して返す（移動中判定などに使える）。
    """
    s0 = str(s).strip().upper()
    if s0 in ("ON", "OFF", "HOT", "SKY"):
        return s0
    if s0 == "":
        return "UNKNOWN"
    return s0




# -----------------------------------------------------------------------------
# 1) Config / stream definition
# -----------------------------------------------------------------------------
@dataclass
class BeamConfig:
    beam_id: str = "B00"
    az_offset_arcsec: float = 0.0  # tangent-plane X = dAz*cosEl
    el_offset_arcsec: float = 0.0  # tangent-plane Y = dEl
    rotation_mode: str = "none"    # none | elevation
    reference_angle_deg: float = 0.0
    rotation_sign: float = 1.0
    dewar_angle_deg: float = 0.0
    beam_model_version: Optional[str] = None


@dataclass
class StreamWCS:
    nchan: int
    crval1_hz: float
    cdelt1_hz: float
    crpix1: float
    ctype1: str
    cunit1: str
    restfreq_hz: float
    specsys: str
    veldef: Optional[str]
    obsfreq_hz: Optional[float]
    imagfreq_hz: Optional[float]
    sideband: Optional[str]
    lo1_hz: Optional[float]
    lo2_hz: Optional[float]
    lo3_hz: Optional[float]
    sb1: Optional[str]
    sb2: Optional[str]
    sb3: Optional[str]
    store_freq_column: Any = "auto"


@dataclass
class StreamConfig:
    name: str
    fdnum: int
    ifnum: int
    plnum: int
    polariza: str
    beam: BeamConfig
    frontend: Optional[str] = None
    backend: Optional[str] = None
    sampler: Optional[str] = None
    db_stream_name: Optional[str] = None
    db_table_name: Optional[str] = None
    frequency_axis: Dict[str, Any] = field(default_factory=dict)
    local_oscillators: Dict[str, Any] = field(default_factory=dict)
    override: Dict[str, Any] = field(default_factory=dict)
    channel_slice_spec: Any = None
    channel_slice_bounds: Optional[Tuple[int, int]] = None
    channel_slice_label: Optional[str] = None
    stream_index: int = 0
    wcs_full: Optional[StreamWCS] = None
    wcs: Optional[StreamWCS] = None


@dataclass
class CommonInputs:
    db: Any
    arr_enc: np.ndarray
    arr_alt: np.ndarray
    wcs_df: Optional[pd.DataFrame]
    weather_inside_table: str
    weather_inside_time_col: str
    weather_outside_table: str
    weather_outside_time_col: str


@dataclass
class StreamData:
    stream: StreamConfig
    t_spec: np.ndarray
    spec2d: np.ndarray
    id_str: np.ndarray
    pos_str: np.ndarray
    temp_c_spectral: Optional[np.ndarray]
    press_hpa_spectral: Optional[np.ndarray]
    humid_pct_spectral: Optional[np.ndarray]
    spec_table_name: str
    spec_time_field: str
    spec_time_basis: str
    spec_time_suffix: Optional[str]
    spec_time_fallback_field: Optional[str]
    spec_time_example: Optional[str]
    spec_data_field: str



def _is_meaningful_str(v):
    return v is not None and str(v).strip() != ""


def _is_meaningful_scalar(v):
    try:
        return v is not None and np.isfinite(float(v))
    except Exception:
        return False


def _load_toml_file(path):
    if tomllib is None:
        raise RuntimeError("tomllib is unavailable. Python 3.11+ is required for --spectrometer-config.")
    with open(path, "rb") as f:
        return tomllib.load(f)


def _as_bool(v, default=False):
    if v is None:
        return bool(default)
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in ("1", "true", "yes", "on"):
        return True
    if s in ("0", "false", "no", "off"):
        return False
    return bool(default)




def _nonempty_str(value, default=None):
    if value is None:
        return default
    s = str(value).strip()
    return s if s else default


def _normalize_table_suffix(value, default):
    s = _nonempty_str(value, default=None)
    return s if s is not None else str(default)


def _resolve_prefixed_table_name(db_namespace, telescope, full_table, suffix, default_suffix):
    full = _nonempty_str(full_table, default=None)
    if full is not None:
        return full
    return _table_name(db_namespace, telescope, _normalize_table_suffix(suffix, default_suffix))

def _table_name(*parts):
    cleaned = []
    for p in parts:
        if p is None:
            continue
        s = str(p).strip()
        if s:
            cleaned.append(s)
    return "-".join(cleaned)
def _clean_upper(v):
    if v is None:
        return None
    s = str(v).strip().upper()
    return s if s else None


def _get_nested(dct, *keys, default=None):
    cur = dct
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def _validate_required(dct, keys, context):
    miss = [k for k in keys if k not in dct or dct[k] is None]
    if miss:
        raise ValueError("{} is missing required keys: {}".format(context, ", ".join(miss)))


def _normalize_rotation_mode(v):
    s = str(v if v is not None else "none").strip().lower()
    if s not in ("none", "elevation"):
        raise ValueError("rotation_mode must be one of ['none', 'elevation'], got {!r}".format(v))
    return s


def _normalize_definition_mode(v):
    s = str(v).strip().lower()
    alias = {
        "band_edges": "band_start_stop",
        "band_start_stop": "band_start_stop",
        "explicit_wcs": "explicit_wcs",
        "first_center_and_delta": "first_center_and_delta",
    }
    if s not in alias:
        raise ValueError("definition_mode must be one of explicit_wcs / first_center_and_delta / band_start_stop, got {!r}".format(v))
    return alias[s]


def _normalize_channel_origin(v):
    s = str(v).strip().lower()
    if s not in ("center", "edge"):
        raise ValueError("channel_origin must be 'center' or 'edge', got {!r}".format(v))
    return s


def _validate_sideband_label(name, value):
    if not _is_meaningful_str(value):
        return None
    s = str(value).strip().upper()
    if s not in ("USB", "LSB"):
        raise ValueError("{} must be USB/LSB, got {!r}".format(name, value))
    return s


def _derive_sideband_from_chain(sb1=None, sb2=None, sb3=None):
    # Final sky sideband can be derived only when the first conversion sideband
    # is known. Additional conversion stages flip the parity when they are LSB.
    s1 = _validate_sideband_label("sb1", sb1)
    s2 = _validate_sideband_label("sb2", sb2)
    s3 = _validate_sideband_label("sb3", sb3)
    if s1 is None:
        return None
    sbs = [s for s in (s1, s2, s3) if s is not None]
    n_lsb = sum(1 for s in sbs if s == "LSB")
    return "USB" if (n_lsb % 2 == 0) else "LSB"


def _derive_imagfreq_hz(obsfreq_hz, lo1_hz):
    if not (_is_meaningful_scalar(obsfreq_hz) and _is_meaningful_scalar(lo1_hz)):
        return None
    return float(2.0 * float(lo1_hz) - float(obsfreq_hz))



def _channel_centers_from_wcs(sw):
    pix = np.arange(int(sw.nchan), dtype=np.float64) + 1.0
    return float(sw.crval1_hz) + (pix - float(sw.crpix1)) * float(sw.cdelt1_hz)


def _parse_channel_slice_spec(value, nchan_full):
    if value is None:
        return None, None

    def _to_endpoint(x, name):
        if isinstance(x, bool):
            raise ValueError(f"channel_slice {name} must be an integer, got boolean {x!r}")
        if isinstance(x, str):
            if x.strip() == '':
                raise ValueError(f"channel_slice {name} must not be empty")
            try:
                return int(x)
            except Exception as e:
                raise ValueError(f"channel_slice {name} must be an integer, got {x!r}") from e
        if isinstance(x, (int, np.integer)):
            return int(x)
        if isinstance(x, (float, np.floating)):
            if not np.isfinite(float(x)) or not float(x).is_integer():
                raise ValueError(f"channel_slice {name} must be an integer, got {x!r}")
            return int(x)
        raise ValueError(f"channel_slice {name} must be an integer, got {x!r}")

    label = None
    if isinstance(value, str):
        s = value.strip()
        m = re.fullmatch(r'([\[(])\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*([\])])', s)
        if not m:
            raise ValueError(
                "channel_slice string must look like '[chmin,chmax]', '[chmin,chmax)', '(chmin,chmax]' or '(chmin,chmax)'"
            )
        left, a_s, b_s, right = m.groups()
        a = _to_endpoint(a_s, 'start')
        b = _to_endpoint(b_s, 'end')
        start = a if left == '[' else a + 1
        stop = b + 1 if right == ']' else b
        label = f"{left}{a},{b}{right}"
    elif isinstance(value, (list, tuple)):
        if len(value) != 2:
            raise ValueError(f"channel_slice sequence must have length 2, got {value!r}")
        a = _to_endpoint(value[0], 'start')
        b = _to_endpoint(value[1], 'end')
        start = a
        stop = b + 1
        label = f"[{a},{b}]"
    else:
        raise ValueError(f"unsupported channel_slice value type: {type(value).__name__}")

    nchan_full = int(nchan_full)
    if start < 0 or stop < 0:
        raise ValueError(f"channel_slice must be non-negative after normalization, got start={start}, stop={stop}")
    if start >= stop:
        raise ValueError(f"channel_slice selects no channels after normalization: start={start}, stop={stop}")
    if stop > nchan_full:
        raise ValueError(
            f"channel_slice stop={stop} exceeds available channel count {nchan_full}; label={label}"
        )
    return (int(start), int(stop)), label


def _slice_stream_wcs(sw, channel_slice_bounds):
    if channel_slice_bounds is None:
        return sw
    start, stop = map(int, channel_slice_bounds)
    if start == 0 and stop == int(sw.nchan):
        return sw
    centers_full = _channel_centers_from_wcs(sw)
    centers = np.asarray(centers_full[start:stop], dtype=np.float64)
    if centers.size == 0:
        raise ValueError(f"empty channel slice after applying bounds {channel_slice_bounds}")
    if centers.size > 2:
        d = np.diff(centers)
        if not np.allclose(d, d[0], rtol=0.0, atol=max(1e-6, abs(d[0]) * 1e-9)):
            raise ValueError(
                f"channel_slice {channel_slice_bounds} produced non-uniform channel centers; this converter supports only uniform sliced axes"
            )
    return replace(
        sw,
        nchan=int(centers.size),
        crval1_hz=float(centers[0]),
        cdelt1_hz=float(sw.cdelt1_hz),
        crpix1=1.0,
    )


def apply_channel_slice_config(streams, global_cfg, cli_channel_slice=None):
    global_slice = None if global_cfg is None else global_cfg.get('channel_slice', None)
    for stream in streams:
        raw_spec = cli_channel_slice if cli_channel_slice is not None else (
            stream.channel_slice_spec if stream.channel_slice_spec is not None else global_slice
        )
        stream.wcs_full = stream.wcs_full or stream.wcs
        if raw_spec is None:
            stream.channel_slice_bounds = None
            stream.channel_slice_label = None
            stream.wcs = stream.wcs_full
            continue
        bounds, label = _parse_channel_slice_spec(raw_spec, int(stream.wcs_full.nchan))
        stream.channel_slice_bounds = bounds
        stream.channel_slice_label = label
        stream.wcs = _slice_stream_wcs(stream.wcs_full, bounds)


def _apply_channel_slice_to_spec2d(spec2d, channel_slice_bounds, expected_nchan):
    arr = np.asarray(spec2d)
    if arr.ndim != 2:
        raise RuntimeError(f"spectrum array must be 2D before channel slicing. got {arr.shape}")
    if int(arr.shape[1]) != int(expected_nchan):
        raise RuntimeError(
            f"NCHAN mismatch in spectral table before channel slicing: got {arr.shape[1]}, expected {expected_nchan}"
        )
    if channel_slice_bounds is None:
        return arr
    start, stop = map(int, channel_slice_bounds)
    out = arr[:, start:stop]
    if out.shape[1] != int(stop - start):
        raise RuntimeError(
            f"channel_slice {channel_slice_bounds} produced inconsistent sliced shape {out.shape}"
        )
    return out


def _build_baseband_centers(freq_cfg, nchan):
    mode = _normalize_definition_mode(freq_cfg.get("definition_mode"))
    crpix1 = float(freq_cfg.get("crpix1", 1.0))
    if mode == "explicit_wcs":
        crval1 = float(freq_cfg["crval1_hz"])
        cdelt1 = float(freq_cfg["cdelt1_hz"])
        pix = np.arange(int(nchan), dtype=np.float64) + 1.0
        centers = crval1 + (pix - crpix1) * cdelt1
        return centers, crval1, cdelt1, crpix1

    if mode == "first_center_and_delta":
        first_center = float(freq_cfg["first_channel_center_hz"])
        spacing = float(freq_cfg["channel_spacing_hz"])
        pix = np.arange(int(nchan), dtype=np.float64)
        centers = first_center + pix * spacing
        return centers, float(centers[0]), float(spacing), 1.0

    # band_start_stop
    origin = _normalize_channel_origin(freq_cfg["channel_origin"])
    band_start = float(freq_cfg["band_start_hz"])
    band_stop = float(freq_cfg["band_stop_hz"])
    if int(nchan) <= 0:
        raise ValueError("nchan must be > 0")
    if origin == "center":
        if int(nchan) == 1:
            if not np.isclose(band_start, band_stop, rtol=0.0, atol=max(1e-6, max(abs(band_start), abs(band_stop)) * 1e-12)):
                raise ValueError("frequency_axis with nchan=1 and channel_origin='center' requires band_start_hz == band_stop_hz")
            centers = np.array([band_start], dtype=np.float64)
            spacing = 0.0
        else:
            spacing = (band_stop - band_start) / float(int(nchan) - 1)
            centers = band_start + spacing * np.arange(int(nchan), dtype=np.float64)
        return centers, float(centers[0]), float(spacing), 1.0

    spacing = (band_stop - band_start) / float(int(nchan))
    centers = band_start + (np.arange(int(nchan), dtype=np.float64) + 0.5) * spacing
    return centers, float(centers[0]), float(spacing), 1.0


def _apply_reverse_to_centers(centers, reverse):
    arr = np.asarray(centers, dtype=np.float64)
    if _as_bool(reverse, False):
        arr = arr[::-1].copy()
    return arr


def _apply_lo_chain_to_sky(baseband_centers_hz, lo_cfg):
    arr = np.asarray(baseband_centers_hz, dtype=np.float64)
    stages = []
    for i in (1, 2, 3):
        lo_key = "lo{}_hz".format(i)
        sb_key = "sb{}".format(i)
        lo = lo_cfg.get(lo_key)
        sb = lo_cfg.get(sb_key)
        if lo is None and sb is None:
            continue
        if lo is None or not _is_meaningful_str(sb):
            raise ValueError("LO chain requires both {} and {}".format(lo_key, sb_key))
        sb_u = str(sb).strip().upper()
        if sb_u not in ("USB", "LSB"):
            raise ValueError("{} must be USB/LSB, got {!r}".format(sb_key, sb))
        stages.append((float(lo), sb_u))
    out = arr.copy()
    for lo, sb in reversed(stages):
        if sb == "USB":
            out = lo + out
        else:
            out = lo - out
    return out


def derive_stream_wcs(stream):
    freq_cfg = dict(stream.frequency_axis or {})
    lo_cfg = dict(stream.local_oscillators or {})
    nchan = int(freq_cfg.get("nchan", 0))
    if nchan <= 0:
        raise ValueError("stream '{}' has invalid nchan={}".format(stream.name, nchan))

    ctype1 = str(freq_cfg.get("ctype1", "FREQ")).strip().upper()
    cunit1 = str(freq_cfg.get("cunit1", "Hz")).strip()
    specsys = str(freq_cfg.get("specsys", "TOPOCENT")).strip().upper()
    if ctype1 != "FREQ":
        raise ValueError("stream '{}' currently supports only CTYPE1='FREQ'; non-frequency output axes are not implemented in the converter".format(stream.name))
    if specsys not in ("TOPOCENT", "LSRK"):
        raise ValueError("stream '{}' currently supports only SPECSYS='TOPOCENT' or 'LSRK'".format(stream.name))
    veldef = freq_cfg.get("veldef", "RADIO")
    veldef = None if veldef is None else str(veldef).strip().upper()
    restfreq_hz = float(freq_cfg.get("restfreq_hz", np.nan))
    store_freq_column = freq_cfg.get("store_freq_column", "auto")

    mode = _normalize_definition_mode(freq_cfg.get("definition_mode"))
    if mode == "explicit_wcs":
        crval1 = float(freq_cfg["crval1_hz"])
        cdelt1 = float(freq_cfg["cdelt1_hz"])
        crpix1 = float(freq_cfg.get("crpix1", 1.0))
        axis_centers = crval1 + ((np.arange(nchan, dtype=np.float64) + 1.0) - crpix1) * cdelt1
    else:
        _validate_required(freq_cfg, ["definition_mode"], "frequency_axis")
        if mode == "first_center_and_delta":
            _validate_required(freq_cfg, ["first_channel_center_hz", "channel_spacing_hz"], "frequency_axis")
        elif mode == "band_start_stop":
            _validate_required(freq_cfg, ["band_start_hz", "band_stop_hz", "channel_origin"], "frequency_axis")
        baseband_centers, _, _, _ = _build_baseband_centers(freq_cfg, nchan)
        baseband_centers = _apply_reverse_to_centers(baseband_centers, freq_cfg.get("reverse", False))
        axis_centers = _apply_lo_chain_to_sky(baseband_centers, lo_cfg)
        crpix1 = 1.0
        crval1 = float(axis_centers[0])
        if nchan > 1:
            cdelt1 = float(axis_centers[1] - axis_centers[0])
        else:
            cdelt1 = 0.0
        if nchan > 2:
            d = np.diff(axis_centers)
            if not np.allclose(d, d[0], rtol=0.0, atol=max(1e-6, abs(d[0]) * 1e-9)):
                raise ValueError("stream '{}' produced non-uniform sky-frequency centers from LO chain".format(stream.name))

    sb1 = _validate_sideband_label("sb1", lo_cfg.get("sb1"))
    sb2 = _validate_sideband_label("sb2", lo_cfg.get("sb2"))
    sb3 = _validate_sideband_label("sb3", lo_cfg.get("sb3"))
    lo1_hz = lo_cfg.get("lo1_hz")
    lo2_hz = lo_cfg.get("lo2_hz")
    lo3_hz = lo_cfg.get("lo3_hz")
    obsfreq_explicit = _is_meaningful_scalar(lo_cfg.get("obsfreq_hz", None))
    imagfreq_explicit = _is_meaningful_scalar(lo_cfg.get("imagfreq_hz", None))
    sideband = _derive_sideband_from_chain(sb1, sb2, sb3)
    heterodyne_context_requested = any([
        obsfreq_explicit,
        imagfreq_explicit,
        _is_meaningful_scalar(lo1_hz),
        _is_meaningful_scalar(lo2_hz),
        _is_meaningful_scalar(lo3_hz),
        _is_meaningful_str(sb1),
        _is_meaningful_str(sb2),
        _is_meaningful_str(sb3),
    ])
    if heterodyne_context_requested and sideband is None:
        raise ValueError(
            "stream '{}' has LO/sideband metadata but final sideband cannot be derived; please specify a valid sb1 (and additional sb2/sb3 if needed)".format(stream.name)
        )

    obsfreq_hz = lo_cfg.get("obsfreq_hz", None)
    if obsfreq_explicit:
        obsfreq_hz = float(obsfreq_hz)
    elif heterodyne_context_requested and np.isfinite(restfreq_hz) and float(restfreq_hz) > 0:
        # Compatibility fallback for writer-side heterodyne metadata.  This is
        # not used to build the actual spectral axis, and is enabled only when
        # the row actually carries LO/sideband semantics.
        obsfreq_hz = float(restfreq_hz)
    else:
        obsfreq_hz = None
    if heterodyne_context_requested and not _is_meaningful_scalar(obsfreq_hz):
        raise ValueError(
            "stream '{}' has LO/sideband metadata but no OBSFREQ can be determined; specify local_oscillators.obsfreq_hz or a positive restfreq_hz".format(stream.name)
        )

    imagfreq_hz = lo_cfg.get("imagfreq_hz", None)
    if imagfreq_explicit:
        imagfreq_hz = float(imagfreq_hz)
    elif obsfreq_explicit:
        # Derive IMAGFREQ only from an explicit OBSFREQ.  When OBSFREQ is a
        # compatibility fallback from RESTFREQ, deriving IMAGFREQ would create
        # misleading tuning metadata.
        imagfreq_hz = _derive_imagfreq_hz(obsfreq_hz, lo1_hz)
    else:
        imagfreq_hz = None

    return StreamWCS(
        nchan=nchan,
        crval1_hz=float(crval1),
        cdelt1_hz=float(cdelt1),
        crpix1=float(crpix1),
        ctype1=ctype1,
        cunit1=cunit1,
        restfreq_hz=float(restfreq_hz) if np.isfinite(restfreq_hz) else np.nan,
        specsys=specsys,
        veldef=veldef,
        obsfreq_hz=obsfreq_hz,
        imagfreq_hz=imagfreq_hz,
        sideband=sideband,
        lo1_hz=(float(lo1_hz) if _is_meaningful_scalar(lo1_hz) else None),
        lo2_hz=(float(lo2_hz) if _is_meaningful_scalar(lo2_hz) else None),
        lo3_hz=(float(lo3_hz) if _is_meaningful_scalar(lo3_hz) else None),
        sb1=sb1,
        sb2=sb2,
        sb3=sb3,
        store_freq_column=store_freq_column,
    )


def _normalize_stream_block(block, stream_index):
    if "name" not in block:
        raise ValueError("Each [[spectrometers]] block must have 'name'")
    name = str(block["name"]).strip()
    if not name:
        raise ValueError("spectrometer.name must be non-empty")
    fdnum = int(block.get("fdnum", 0))
    ifnum = int(block.get("ifnum", 0))
    plnum = int(block.get("plnum", 0))
    polariza = _clean_upper(block.get("polariza"))
    if polariza is None:
        raise ValueError("stream '{}' requires polariza".format(name))
    beam_id = str(block.get("beam_id", "B{:02d}".format(fdnum))).strip() or "B{:02d}".format(fdnum)

    beam_block = dict(block.get("beam", {}) or {})
    beam = BeamConfig(
        beam_id=beam_id,
        az_offset_arcsec=float(beam_block.get("az_offset_arcsec", 0.0)),
        el_offset_arcsec=float(beam_block.get("el_offset_arcsec", 0.0)),
        rotation_mode=_normalize_rotation_mode(beam_block.get("rotation_mode", "none")),
        reference_angle_deg=float(beam_block.get("reference_angle_deg", beam_block.get("reference_el_deg", 0.0))),
        rotation_sign=float(beam_block.get("rotation_sign", 1.0)),
        dewar_angle_deg=float(beam_block.get("dewar_angle_deg", 0.0)),
        beam_model_version=(str(beam_block.get("beam_model_version")).strip() if beam_block.get("beam_model_version") is not None else None),
    )

    frequency_axis = dict(block.get("frequency_axis", {}) or {})
    local_oscillators = dict(block.get("local_oscillators", {}) or {})
    override = dict(block.get("override", {}) or {})
    channel_slice_spec = block.get("channel_slice", None)
    if channel_slice_spec is None and "channel_slice" in frequency_axis:
        channel_slice_spec = frequency_axis.pop("channel_slice")

    stream = StreamConfig(
        name=name,
        fdnum=fdnum,
        ifnum=ifnum,
        plnum=plnum,
        polariza=polariza,
        beam=beam,
        frontend=(str(block.get("frontend")).strip() if block.get("frontend") is not None else None),
        backend=(str(block.get("backend")).strip() if block.get("backend") is not None else None),
        sampler=(str(block.get("sampler")).strip() if block.get("sampler") is not None else None),
        db_stream_name=_nonempty_str(block.get("db_stream_name"), default=name),
        db_table_name=_nonempty_str(block.get("db_table_name"), default=None),
        frequency_axis=frequency_axis,
        local_oscillators=local_oscillators,
        override=override,
        channel_slice_spec=channel_slice_spec,
        stream_index=int(stream_index),
    )
    stream.wcs_full = derive_stream_wcs(stream)
    stream.wcs = stream.wcs_full
    return stream


def load_spectrometer_config(config_path):
    cfg = _load_toml_file(str(config_path))
    version_raw = cfg.get("schema_version", cfg.get("config_version", 1))
    schema_version = int(version_raw)
    if schema_version != 1:
        raise ValueError("unsupported schema_version/config_version={} (expected 1)".format(schema_version))
    blocks = list(cfg.get("spectrometers", []) or [])
    if not blocks:
        raise ValueError("TOML requires at least one [[spectrometers]] block")
    streams = [_normalize_stream_block(b, i) for i, b in enumerate(blocks)]

    seen_name = set()
    seen_key = set()
    for s in streams:
        if s.name in seen_name:
            raise ValueError("duplicate stream name '{}'".format(s.name))
        seen_name.add(s.name)
        key = (s.fdnum, s.ifnum, s.plnum, s.polariza, s.sampler)
        if key in seen_key:
            raise ValueError("duplicate (fdnum, ifnum, plnum, polariza, sampler) combination: {}".format(key))
        seen_key.add(key)

    global_cfg = dict(cfg.get("global", {}) or {})
    provenance = dict(cfg.get("provenance", {}) or {})
    return {
        "schema_version": schema_version,
        "config_name": cfg.get("config_name", None),
        "config_description": cfg.get("config_description", None),
        "global": global_cfg,
        "provenance": provenance,
        "streams": streams,
    }


def build_legacy_single_stream_config(args):
    beam = BeamConfig(
        beam_id="B{:02d}".format(0),
        az_offset_arcsec=0.0,
        el_offset_arcsec=0.0,
        rotation_mode="none",
    )
    stream = StreamConfig(
        name=str(args.spectral),
        fdnum=0,
        ifnum=0,
        plnum=0,
        polariza="XX",
        beam=beam,
        frontend=None,
        backend="XFFTS",
        sampler=str(args.spectral),
        frequency_axis={
            "nchan": int(args.nchan),
            "definition_mode": "first_center_and_delta",
            "first_channel_center_hz": float(args.if0_ghz) * 1e9,
            "channel_spacing_hz": ((float(args.if1_ghz) - float(args.if0_ghz)) * 1e9) / max(int(args.nchan) - 1, 1),
            "crpix1": 1.0,
            "restfreq_hz": float(args.restfreq_hz),
            "ctype1": "FREQ",
            "cunit1": "Hz",
            "specsys": "TOPOCENT",
            "veldef": "RADIO",
            "store_freq_column": False,
        },
        local_oscillators={
            "lo1_hz": float(args.lo1_ghz) * 1e9,
            "lo2_hz": float(args.lo2_ghz) * 1e9,
            "sb1": "USB",
            "sb2": "USB",
        },
        channel_slice_spec=None,
        stream_index=0,
    )
    stream.wcs_full = derive_stream_wcs(stream)
    stream.wcs = stream.wcs_full
    return {
        "schema_version": 1,
        "config_name": "legacy_single_stream",
        "config_description": "Auto-generated from legacy CLI arguments",
        "global": {"output_layout": "merged", "time_sort": True, "db_namespace": "necst", "telescope": str(args.telescope), "encoder_shift_sec": float(args.encoder_shift_sec)},
        "provenance": {},
        "streams": [stream],
    }


# -----------------------------------------------------------------------------
# 2) Extract
# -----------------------------------------------------------------------------
def extract_common_inputs(rawdata_path, db_namespace, telescope, wcs_table,
                          weather_inside_table, weather_inside_time_col,
                          weather_outside_table, weather_outside_time_col,
                          encoder_table=None, altaz_table=None):
    db = necstdb.opendb(str(rawdata_path))
    enc_table = _nonempty_str(encoder_table, default=_table_name(db_namespace, telescope, "ctrl", "antenna", "encoder"))
    alt_table = _nonempty_str(altaz_table, default=_table_name(db_namespace, telescope, "ctrl", "antenna", "altaz"))
    arr_enc = _read_structured_array_tolerant(db, enc_table)
    arr_alt = _read_structured_array_tolerant(db, alt_table)
    wcs_df = _try_radec_timeseries_from_wcs_table(db, str(wcs_table))
    return CommonInputs(
        db=db,
        arr_enc=arr_enc,
        arr_alt=arr_alt,
        wcs_df=wcs_df,
        weather_inside_table=str(weather_inside_table),
        weather_inside_time_col=str(weather_inside_time_col),
        weather_outside_table=str(weather_outside_table),
        weather_outside_time_col=str(weather_outside_time_col),
    )


def extract_spectral_stream(common, db_namespace, telescope, stream):
    if _is_meaningful_str(stream.db_table_name):
        spec_table_name = str(stream.db_table_name)
    else:
        spec_stream_name = _nonempty_str(stream.db_stream_name, default=stream.name)
        spec_table_name = _table_name(db_namespace, telescope, "data", "spectral", spec_stream_name)
    arr_spec = _read_structured_array_tolerant(common.db, spec_table_name)
    raw_nchan = int(stream.wcs_full.nchan if stream.wcs_full is not None else stream.wcs.nchan)
    t_spec, spec2d, time_meta, s_field = _extract_spectral_from_structured(arr_spec, nchan=raw_nchan)
    spec2d = _apply_channel_slice_to_spec2d(spec2d, stream.channel_slice_bounds, raw_nchan)
    if int(spec2d.shape[1]) != int(stream.wcs.nchan):
        raise RuntimeError(
            "Sliced NCHAN mismatch for stream '{}': got {}, expected {} after channel_slice={}".format(
                stream.name, spec2d.shape[1], stream.wcs.nchan, stream.channel_slice_label
            )
        )

    names = list(arr_spec.dtype.names or [])
    pos_field = None
    for cand in ("position", "obsmode", "obs_mode", "obsMode", "mode", "state", "status", "label"):
        if cand in names:
            pos_field = cand
            break

    if "id" in names:
        id_str = np.array([_decode_label(v) for v in np.asarray(arr_spec["id"], dtype=object)], dtype=object)
    else:
        id_str = np.array([""] * len(t_spec), dtype=object)

    if pos_field is not None:
        pos_str = np.array([_decode_label(v) for v in np.asarray(arr_spec[pos_field], dtype=object)], dtype=object)
    else:
        pos_str = np.array(["UNKNOWN"] * len(t_spec), dtype=object)
    if np.all([str(x).strip() == "" for x in pos_str]):
        pos_str = np.array(["UNKNOWN"] * len(t_spec), dtype=object)
    else:
        pos_str = np.array([("UNKNOWN" if str(x).strip() == "" else str(x).strip().upper()) for x in pos_str], dtype=object)

    temp_c = None
    press_hpa = None
    humid_pct = None
    if "temperature" in names:
        temp_c = _to_degC_guess(np.asarray(arr_spec["temperature"], dtype=float))
    if "pressure" in names:
        press_hpa = _to_hpa_guess(np.asarray(arr_spec["pressure"], dtype=float))
    if "humidity" in names:
        x = np.asarray(arr_spec["humidity"], dtype=float)
        med = float(np.nanmedian(x)) if x.size > 0 else np.nan
        humid_pct = x * 100.0 if (np.isfinite(med) and med <= 1.5) else x

    order = np.argsort(t_spec)
    return StreamData(
        stream=stream,
        t_spec=t_spec[order],
        spec2d=spec2d[order, :],
        id_str=id_str[order],
        pos_str=pos_str[order],
        temp_c_spectral=(temp_c[order] if temp_c is not None else None),
        press_hpa_spectral=(press_hpa[order] if press_hpa is not None else None),
        humid_pct_spectral=(humid_pct[order] if humid_pct is not None else None),
        spec_table_name=spec_table_name,
        spec_time_field=(time_meta.get("timestamp_field") if str(time_meta.get("applied", "")).startswith(("timestamp", "time_spectrometer")) else str(time_meta.get("applied"))),
        spec_time_basis=str(time_meta.get("applied")),
        spec_time_suffix=time_meta.get("suffix"),
        spec_time_fallback_field=time_meta.get("fallback_field"),
        spec_time_example=time_meta.get("first_timestamp_text"),
        spec_data_field=s_field,
    )


def _choose_meteo_array(met_source, weather_arr, spectral_arr):
    if met_source == "weather":
        if weather_arr is not None and np.any(np.isfinite(np.asarray(weather_arr, dtype=float))):
            return np.asarray(weather_arr, dtype=float)
        return spectral_arr
    if met_source == "spectral":
        if spectral_arr is not None and np.any(np.isfinite(np.asarray(spectral_arr, dtype=float))):
            return np.asarray(spectral_arr, dtype=float)
        return weather_arr
    if met_source == "fallback":
        return None
    # auto
    if weather_arr is not None and np.any(np.isfinite(np.asarray(weather_arr, dtype=float))):
        return np.asarray(weather_arr, dtype=float)
    if spectral_arr is not None and np.any(np.isfinite(np.asarray(spectral_arr, dtype=float))):
        return np.asarray(spectral_arr, dtype=float)
    return None


def resolve_stream_meteorology(common, stream_data, met_source, interp_extrap):
    inside_weather = None
    outside_weather = None
    if met_source in ("auto", "weather"):
        inside_weather = _try_meteo_timeseries_from_weather_table(
            common.db,
            common.weather_inside_table,
            stream_data.t_spec,
            preferred_time_col=common.weather_inside_time_col,
            extrap=interp_extrap,
        )
        if inside_weather is not None:
            print("[info] stream={} inside meteorology from weather table: {} (time_col={})".format(
                stream_data.stream.name, inside_weather.get("table"), inside_weather.get("time_col")
            ))
        if (str(common.weather_outside_table) == str(common.weather_inside_table)) and (str(common.weather_outside_time_col) == str(common.weather_inside_time_col)):
            outside_weather = inside_weather
        else:
            outside_weather = _try_meteo_timeseries_from_weather_table(
                common.db,
                common.weather_outside_table,
                stream_data.t_spec,
                preferred_time_col=common.weather_outside_time_col,
                extrap=interp_extrap,
            )
            if outside_weather is not None:
                print("[info] stream={} outside meteorology from weather table: {} (time_col={})".format(
                    stream_data.stream.name, outside_weather.get("table"), outside_weather.get("time_col")
                ))

    tamb_c = _choose_meteo_array(
        met_source,
        inside_weather.get("temp_c") if inside_weather is not None else None,
        stream_data.temp_c_spectral,
    )
    temp_c = _choose_meteo_array(
        met_source,
        outside_weather.get("temp_c") if outside_weather is not None else None,
        stream_data.temp_c_spectral,
    )
    press_hpa = _choose_meteo_array(
        met_source,
        outside_weather.get("press_hpa") if outside_weather is not None else None,
        stream_data.press_hpa_spectral,
    )
    humid_pct = _choose_meteo_array(
        met_source,
        outside_weather.get("humid_pct") if outside_weather is not None else None,
        stream_data.humid_pct_spectral,
    )
    return tamb_c, temp_c, press_hpa, humid_pct, inside_weather, outside_weather


# -----------------------------------------------------------------------------
# 3) Normalize
# -----------------------------------------------------------------------------
def normalize_pointing(t_spec, arr_enc, arr_alt, encoder_time_col, altaz_time_col, extrap, encoder_shift_sec=0.0):
    enc_names = list(arr_enc.dtype.names or [])
    alt_names = list(arr_alt.dtype.names or [])

    enc_t_name = _pick_field_name(enc_names, encoder_time_col, ["time", "timestamp"])
    alt_t_name = _pick_field_name(alt_names, altaz_time_col, ["time", "timestamp"])
    if enc_t_name is None:
        raise RuntimeError("encoder time column not found. preferred='{}' available={}".format(encoder_time_col, enc_names))
    if alt_t_name is None:
        raise RuntimeError("altaz time column not found. preferred='{}' available={}".format(altaz_time_col, alt_names))

    enc_lon_name = _pick_field_name(enc_names, None, ["lon", "az", "azimuth"])
    enc_lat_name = _pick_field_name(enc_names, None, ["lat", "el", "elevation", "alt"])
    if enc_lon_name is None or enc_lat_name is None:
        raise RuntimeError("encoder lon/lat not found. available={}".format(enc_names))

    alt_lon_name = _pick_field_name(alt_names, None, ["lon", "az", "azimuth"])
    alt_lat_name = _pick_field_name(alt_names, None, ["lat", "el", "elevation", "alt"])
    dlon_name = _pick_field_name(alt_names, None, ["dlon", "daz", "d_az"])
    dlat_name = _pick_field_name(alt_names, None, ["dlat", "del", "d_el"])
    if dlon_name is None or dlat_name is None:
        raise RuntimeError("altaz dlon/dlat not found. available={}".format(alt_names))

    t_enc = _get_field(arr_enc, enc_t_name, float)
    t_enc, s_enc = _normalize_time_units(t_spec, t_enc, "encoder")
    # Keep the sign convention identical to sunscan: shift encoder timestamps,
    # then interpolate encoder Az/El onto spectral timestamps t_spec.
    # Positive encoder_shift_sec means the encoder stream itself is moved later
    # in time (t_enc <- t_enc + shift).
    t_enc = t_enc + float(encoder_shift_sec)
    az_enc = _get_field(arr_enc, enc_lon_name, float)
    el_enc = _get_field(arr_enc, enc_lat_name, float)

    t_alt = _get_field(arr_alt, alt_t_name, float)
    t_alt, s_alt = _normalize_time_units(t_spec, t_alt, "altaz")
    dlon = _get_field(arr_alt, dlon_name, float)
    dlat = _get_field(arr_alt, dlat_name, float)

    if s_enc != 1.0:
        print("[info] encoder time scaled by {} to match unix seconds".format(s_enc))
    if s_alt != 1.0:
        print("[info] altaz   time scaled by {} to match unix seconds".format(s_alt))

    az_enc_t = _interp_az_deg(t_enc, az_enc, t_spec, extrap=extrap)
    el_enc_t = _interp_lin(t_enc, el_enc, t_spec, extrap=extrap)

    dlon_t = _interp_lin(t_alt, dlon, t_spec, extrap=extrap)
    dlat_t = _interp_lin(t_alt, dlat, t_spec, extrap=extrap)

    if (alt_lon_name is not None) and (alt_lat_name is not None):
        az_alt = _get_field(arr_alt, alt_lon_name, float)
        el_alt = _get_field(arr_alt, alt_lat_name, float)
        az_cmd_t = _interp_az_deg(t_alt, az_alt, t_spec, extrap=extrap)
        el_cmd_t = _interp_lin(t_alt, el_alt, t_spec, extrap=extrap)
    else:
        az_cmd_t = np.full_like(t_spec, np.nan, dtype=float)
        el_cmd_t = np.full_like(t_spec, np.nan, dtype=float)

    boresight_az = (az_enc_t - dlon_t) % 360.0
    boresight_el = el_enc_t - dlat_t

    d_az = _wrap180_deg(az_enc_t - az_cmd_t)
    d_el = (el_enc_t - el_cmd_t)
    cos_el_bore = np.cos(np.deg2rad(boresight_el))
    pe_x_arcsec = d_az * cos_el_bore * 3600.0
    pe_y_arcsec = d_el * 3600.0
    pe_r_arcsec = np.sqrt(pe_x_arcsec**2 + pe_y_arcsec**2)

    return {
        "az_cmd_t": az_cmd_t,
        "el_cmd_t": el_cmd_t,
        "az_enc_t": az_enc_t,
        "el_enc_t": el_enc_t,
        "dlon_t": dlon_t,
        "dlat_t": dlat_t,
        "boresight_az": boresight_az,
        "boresight_el": boresight_el,
        "pe_x_arcsec": pe_x_arcsec,
        "pe_y_arcsec": pe_y_arcsec,
        "pe_r_arcsec": pe_r_arcsec,
        "enc_fields": {"time": enc_t_name, "lon": enc_lon_name, "lat": enc_lat_name},
        "alt_fields": {"time": alt_t_name, "lon": alt_lon_name, "lat": alt_lat_name, "dlon": dlon_name, "dlat": dlat_name},
    }


def _beam_rotation_angle_deg(boresight_el_deg, beam):
    if beam.rotation_mode == "none":
        return float(beam.dewar_angle_deg)
    if beam.rotation_mode == "elevation":
        return float(beam.rotation_sign) * (np.asarray(boresight_el_deg, dtype=float) - float(beam.reference_angle_deg)) + float(beam.dewar_angle_deg)
    raise ValueError("unsupported rotation_mode={!r}".format(beam.rotation_mode))


def apply_beam_offset(boresight_az_deg, boresight_el_deg, beam):
    dx0 = float(beam.az_offset_arcsec)
    dy0 = float(beam.el_offset_arcsec)
    el = np.asarray(boresight_el_deg, dtype=float)
    az = np.asarray(boresight_az_deg, dtype=float)
    theta_deg = np.asarray(_beam_rotation_angle_deg(el, beam), dtype=float)
    if theta_deg.shape == ():
        theta_deg = np.full_like(el, float(theta_deg), dtype=float)
    theta = np.deg2rad(theta_deg)

    dx_rot = np.asarray(dx0 * np.cos(theta) - dy0 * np.sin(theta), dtype=float)
    dy_rot = np.asarray(dx0 * np.sin(theta) + dy0 * np.cos(theta), dtype=float)
    if dx_rot.shape == ():
        dx_rot = np.full_like(el, float(dx_rot), dtype=float)
    if dy_rot.shape == ():
        dy_rot = np.full_like(el, float(dy_rot), dtype=float)

    cos_el = np.cos(np.deg2rad(el))
    tiny = np.cos(np.deg2rad(89.9))
    safe = np.abs(cos_el) >= tiny
    d_az_deg = np.full_like(el, np.nan, dtype=float)
    d_az_deg[safe] = dx_rot[safe] / (3600.0 * cos_el[safe])
    d_el_deg = dy_rot / 3600.0

    beam_az = (az + d_az_deg) % 360.0
    beam_el = el + d_el_deg

    return {
        "beam_az_deg": beam_az,
        "beam_el_deg": beam_el,
        "beam_dx_arcsec": dx_rot,
        "beam_dy_arcsec": dy_rot,
        "beam_rot_deg": theta_deg,
    }


def select_radec_azel_source(pointing, beam, radec_azel_source):
    source = str(radec_azel_source if radec_azel_source is not None else "beam").strip().lower()

    if source == "beam":
        return {
            "source": "beam",
            "az_deg": np.asarray(pointing["beam_az_deg"], dtype=float),
            "el_deg": np.asarray(pointing["beam_el_deg"], dtype=float),
            "meaning": "beam-center Az/El from (encoder - correction) boresight + beam offset",
        }

    if source == "true":
        return {
            "source": "true",
            "az_deg": np.asarray(pointing["boresight_az"], dtype=float),
            "el_deg": np.asarray(pointing["boresight_el"], dtype=float),
            "meaning": "boresight Az/El from encoder - correction (no beam offset)",
        }

    if source == "encoder":
        applied = apply_beam_offset(pointing["az_enc_t"], pointing["el_enc_t"], beam)
        return {
            "source": "encoder",
            "az_deg": np.asarray(applied["beam_az_deg"], dtype=float),
            "el_deg": np.asarray(applied["beam_el_deg"], dtype=float),
            "meaning": "beam-offset Az/El from measured encoder lon/lat",
        }

    if source == "altaz":
        applied = apply_beam_offset(pointing["az_cmd_t"], pointing["el_cmd_t"], beam)
        return {
            "source": "altaz",
            "az_deg": np.asarray(applied["beam_az_deg"], dtype=float),
            "el_deg": np.asarray(applied["beam_el_deg"], dtype=float),
            "meaning": "beam-offset Az/El from raw altaz/cmd lon/lat",
        }

    if source == "cmd":
        cmd_boresight_az = (np.asarray(pointing["az_cmd_t"], dtype=float) - np.asarray(pointing["dlon_t"], dtype=float)) % 360.0
        cmd_boresight_el = np.asarray(pointing["el_cmd_t"], dtype=float) - np.asarray(pointing["dlat_t"], dtype=float)
        applied = apply_beam_offset(cmd_boresight_az, cmd_boresight_el, beam)
        return {
            "source": "cmd",
            "az_deg": np.asarray(applied["beam_az_deg"], dtype=float),
            "el_deg": np.asarray(applied["beam_el_deg"], dtype=float),
            "meaning": "beam-offset Az/El from (altaz/cmd - correction) boresight",
        }

    raise ValueError("unsupported radec_azel_source={!r}".format(radec_azel_source))


def normalize_radec(wcs_df, t_spec, beam_az_deg, beam_el_deg, site,
                    radec_method="wcs_first",
                    apply_refraction=True,
                    press_hpa=None,
                    temp_c=None,
                    humid_pct=None,
                    obswl_um=None):
    t_index = pd.to_datetime(t_spec, unit="s")
    use_wcs = (str(radec_method).strip().lower() == "wcs_first") and (wcs_df is not None)

    ra_deg = None
    dec_deg = None
    used_wcs = np.zeros(len(t_spec), dtype=bool)

    if use_wcs:
        try:
            ra0 = float(wcs_df["ra_deg"].iloc[-1])
            dec0 = float(wcs_df["dec_deg"].iloc[-1])
            ra_deg = _nearest_series(wcs_df, t_index, "ra_deg", ra0)
            dec_deg = _nearest_series(wcs_df, t_index, "dec_deg", dec0)
            used_wcs[:] = np.isfinite(ra_deg) & np.isfinite(dec_deg)
        except Exception:
            ra_deg = None
            dec_deg = None
            used_wcs[:] = False

    if ra_deg is None or dec_deg is None or np.any(~used_wcs):
        rh01 = _to_rh01_guess(humid_pct) if humid_pct is not None else None
        ra2, dec2 = _radec_from_azel(
            site, t_spec, beam_az_deg, beam_el_deg,
            apply_refraction=bool(apply_refraction),
            press_hpa=press_hpa,
            temp_c=temp_c,
            rh01=rh01,
            obswl_um=obswl_um,
        )
        if ra_deg is None or dec_deg is None:
            ra_deg, dec_deg = ra2, dec2
            used_wcs[:] = False
        else:
            ra_deg = np.where(used_wcs, ra_deg, ra2)
            dec_deg = np.where(used_wcs, dec_deg, dec2)

    if np.any(~np.isfinite(ra_deg)) or np.any(~np.isfinite(dec_deg)):
        nbad = int(np.count_nonzero(~np.isfinite(ra_deg)) + np.count_nonzero(~np.isfinite(dec_deg)))
        raise RuntimeError("RA/DEC still contains NaN after conversion ({} elements). Check Az/El validity and column mapping.".format(nbad))

    return {"ra_deg": ra_deg, "dec_deg": dec_deg, "used_wcs": used_wcs}


def _gal_from_radec(ra_deg, dec_deg):
    ra = np.asarray(ra_deg, dtype=float)
    dec = np.asarray(dec_deg, dtype=float)
    if ra.shape != dec.shape:
        raise ValueError("RA/DEC shape mismatch for Galactic conversion: {} vs {}".format(ra.shape, dec.shape))
    sc = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    gal = sc.galactic
    return {
        "glon_deg": np.asarray(gal.l.to_value(u.deg), dtype=float),
        "glat_deg": np.asarray(gal.b.to_value(u.deg), dtype=float),
    }


def _calc_refraction_deg(site, t_spec, az_deg, el_deg, press_hpa, temp_c, humid_pct, obswl_um):
    if press_hpa is None or temp_c is None or humid_pct is None:
        return np.full_like(np.asarray(t_spec, dtype=float), np.nan, dtype=float)
    t = Time(np.asarray(t_spec, dtype=float), format="unix", scale="utc")
    loc = site.to_earthlocation()
    az = np.asarray(az_deg, dtype=float) * u.deg
    alt = np.asarray(el_deg, dtype=float) * u.deg
    p = np.asarray(press_hpa, dtype=float) * u.hPa
    tc = np.asarray(temp_c, dtype=float) * u.deg_C
    rh = _to_rh01_guess(humid_pct)
    ow = (float(obswl_um) if (obswl_um is not None and np.isfinite(obswl_um)) else 2600.0) * u.micron
    try:
        aa_refr = SkyCoord(az=az, alt=alt, frame=AltAz(obstime=t, location=loc, pressure=p, temperature=tc, relative_humidity=rh, obswl=ow))
        aa_no = aa_refr.transform_to(AltAz(obstime=t, location=loc, pressure=0 * u.hPa))
        return np.asarray((aa_refr.alt - aa_no.alt).to_value(u.deg), dtype=float)
    except Exception:
        return np.full_like(np.asarray(t_spec, dtype=float), np.nan, dtype=float)


def _obswl_um_from_stream_or_args(stream_wcs, args):
    if getattr(args, "met_obswl_um", None) is not None and np.isfinite(float(args.met_obswl_um)):
        return float(args.met_obswl_um)
    c_ms = 299792458.0
    # Refraction should use the observing wavelength, not the analysis RESTFREQ.
    # Prefer explicit OBSFREQ, then a representative frequency near the written
    # axis center, and only finally RESTFREQ.
    axis_center_hz = np.nan
    try:
        nchan = int(getattr(stream_wcs, "nchan", 0))
        crval1 = float(getattr(stream_wcs, "crval1_hz", np.nan))
        cdelt1 = float(getattr(stream_wcs, "cdelt1_hz", np.nan))
        crpix1 = float(getattr(stream_wcs, "crpix1", np.nan))
        if nchan > 0 and np.isfinite(crval1) and np.isfinite(cdelt1) and np.isfinite(crpix1):
            mid_pix = 0.5 * (nchan + 1.0)
            axis_center_hz = crval1 + (mid_pix - crpix1) * cdelt1
    except Exception:
        axis_center_hz = np.nan
    candidates = [
        getattr(stream_wcs, "obsfreq_hz", np.nan),
        axis_center_hz,
        getattr(stream_wcs, "restfreq_hz", np.nan),
    ]
    for nu in candidates:
        try:
            nu_hz = float(nu)
        except Exception:
            nu_hz = np.nan
        if np.isfinite(nu_hz) and nu_hz > 0:
            return (c_ms / nu_hz) * 1e6
    return 2600.0


def _fill_with_fallback(arr, fallback):
    if arr is None:
        return np.asarray(fallback, dtype=float)
    x = np.asarray(arr, dtype=float)
    out = x.copy()
    if np.all(~np.isfinite(out)):
        return np.asarray(fallback, dtype=float)
    med = float(np.nanmedian(out))
    fb = med if np.isfinite(med) else np.nanmedian(np.asarray(fallback, dtype=float))
    out = np.where(np.isfinite(out), out, fb)
    return out


def finalize_inside_temperature_for_tamb(t_spec, tamb_c, args):
    inside_default = getattr(args, "inside_default_temperature_c", None)
    if inside_default is None:
        inside_default = getattr(args, "met_temperature_c", 10.0)
    temp_fb = float(inside_default)
    out = _fill_with_fallback(tamb_c, np.full_like(t_spec, temp_fb, dtype=float))
    tmin = float(getattr(args, "inside_temperature_min_c", -50.0))
    tmax = float(getattr(args, "inside_temperature_max_c", 50.0))
    bad = (~np.isfinite(out)) | (out < tmin) | (out > tmax)
    out = np.where(bad, temp_fb, out)
    return out


def finalize_meteorology_for_refraction(t_spec, temp_c, press_hpa, humid_pct, args):
    outside_press_default = getattr(args, "outside_default_pressure_hpa", None)
    if outside_press_default is None:
        outside_press_default = float(args.met_pressure_hpa) if getattr(args, "met_pressure_hpa", None) is not None else _pressure_from_elev_hpa(float(args.site_elev))
    outside_temp_default = getattr(args, "outside_default_temperature_c", None)
    if outside_temp_default is None:
        outside_temp_default = getattr(args, "met_temperature_c", 10.0)
    outside_humid_default = getattr(args, "outside_default_humidity_pct", None)
    if outside_humid_default is None:
        outside_humid_default = getattr(args, "met_humidity_pct", 50.0)
    press_fb = float(outside_press_default)
    temp_fb = float(outside_temp_default)
    humid_fb = float(outside_humid_default)
    press = _fill_with_fallback(press_hpa, np.full_like(t_spec, press_fb, dtype=float))
    temp = _fill_with_fallback(temp_c, np.full_like(t_spec, temp_fb, dtype=float))
    humid = _fill_with_fallback(humid_pct, np.full_like(t_spec, humid_fb, dtype=float))
    pmin = float(getattr(args, "outside_pressure_min_hpa", 400.0))
    pmax = float(getattr(args, "outside_pressure_max_hpa", 1100.0))
    tmin = float(getattr(args, "outside_temperature_min_c", -50.0))
    tmax = float(getattr(args, "outside_temperature_max_c", 50.0))
    hmin = float(getattr(args, "outside_humidity_min_pct", 0.0))
    hmax = float(getattr(args, "outside_humidity_max_pct", 100.0))
    press = np.where((~np.isfinite(press)) | (press < pmin) | (press > pmax), press_fb, press)
    temp = np.where((~np.isfinite(temp)) | (temp < tmin) | (temp > tmax), temp_fb, temp)
    humid = np.where((~np.isfinite(humid)) | (humid < hmin) | (humid > hmax), humid_fb, humid)
    return temp, press, humid


def _scan_key(i, id_str, mode_str, scan_key):
    if str(scan_key) == "id":
        return (str(id_str[i]).strip(),)
    return (str(id_str[i]).strip(), str(mode_str[i]).strip())


def iter_scan_blocks_indices(idx, id_str, mode_str, scan_key):
    idx = np.asarray(idx, dtype=int)
    if idx.size == 0:
        return
    k0 = 0
    key0 = _scan_key(idx[0], id_str, mode_str, scan_key)
    for k in range(1, idx.size):
        key = _scan_key(idx[k], id_str, mode_str, scan_key)
        if key != key0:
            yield (k0, k)
            k0 = k
            key0 = key
    yield (k0, idx.size)


def build_rows_for_stream(stream_data, pointing, radec, calc_refr_deg, temp_c, press_hpa, humid_pct,
                          collapse, scan_key, skip_nonstandard_position, use_modes_csv, include_unknown,
                          pe_rms_warn_arcsec, pe_max_warn_arcsec, verbose_scan_pe):
    t_spec = stream_data.t_spec
    spec2d = stream_data.spec2d
    id_str = stream_data.id_str
    pos_str = stream_data.pos_str
    mode_str = np.array([_mode_from_pos(x) for x in pos_str], dtype=object)
    used_modes = [m.strip().upper() for m in str(use_modes_csv).split(',') if m.strip()]
    used_set = set(used_modes)
    if include_unknown:
        used_set.add('UNKNOWN')

    uniq, cnt = np.unique(mode_str, return_counts=True)
    extra = []
    n_unknown = 0
    for u, c in zip(uniq, cnt):
        us = str(u)
        if us == "UNKNOWN":
            n_unknown = int(c)
        elif us not in used_set:
            extra.append((us, int(c)))
    if n_unknown > 0:
        print("[info] stream={} UNKNOWN position count: {}".format(stream_data.stream.name, int(n_unknown)))
    if extra:
        extra_s = ", ".join(["{}:{}".format(u, c) for (u, c) in extra[:30]])
        print("[info] stream={} nonstandard position labels found: {}".format(stream_data.stream.name, extra_s))

    m_used = np.array([str(x) in used_set for x in mode_str], dtype=bool)
    idx = np.where(m_used)[0]
    if idx.size == 0:
        raise RuntimeError("No samples left after mode selection for stream '{}'".format(stream_data.stream.name))

    if skip_nonstandard_position:
        strict = set(['ON', 'OFF', 'HOT', 'SKY'])
        m2 = np.array([str(x) in strict for x in mode_str[idx]], dtype=bool)
        idx = idx[m2]
        if idx.size == 0:
            raise RuntimeError("All samples were excluded by --skip-nonstandard-position for stream '{}'".format(stream_data.stream.name))

    rows = []
    scanid = -1
    for (k0, k1) in iter_scan_blocks_indices(idx, id_str, mode_str, scan_key):
        scanid += 1
        subscan = 0
        ii = idx[k0:k1]
        if ii.size == 0:
            continue

        ra_src = float(np.nanmedian(radec["ra_deg"][ii]))
        dec_src = float(np.nanmedian(radec["dec_deg"][ii]))
        pe_r = np.asarray(pointing["pe_r_arcsec"][ii], dtype=float)
        pe_x = np.asarray(pointing["pe_x_arcsec"][ii], dtype=float)
        pe_y = np.asarray(pointing["pe_y_arcsec"][ii], dtype=float)
        pe_mean = float(np.nanmean(pe_r)) if np.any(np.isfinite(pe_r)) else float("nan")
        pe_std = float(np.nanstd(pe_r)) if np.any(np.isfinite(pe_r)) else float("nan")
        pe_x_mean = float(np.nanmean(pe_x)) if np.any(np.isfinite(pe_x)) else float("nan")
        pe_y_mean = float(np.nanmean(pe_y)) if np.any(np.isfinite(pe_y)) else float("nan")
        id_label = str(id_str[int(ii[0])]).strip()
        mode_label = str(mode_str[int(ii[0])]).strip()
        t0 = float(t_spec[int(ii[0])])
        t1 = float(t_spec[int(ii[-1])])

        if verbose_scan_pe:
            pe2 = pe_x**2 + pe_y**2
            pe2 = pe2[np.isfinite(pe2)]
            pe_rms0 = float(np.sqrt(np.nanmean(pe2))) if pe2.size > 0 else float("nan")
            pe_abs = np.sqrt(pe2) if pe2.size > 0 else np.asarray([], dtype=float)
            pe_max0 = float(np.nanmax(pe_abs)) if pe_abs.size > 0 else float("nan")
            key_label = id_label if str(scan_key) == "id" else "{}|{}".format(id_label, mode_label)
            thr_rms = float(pe_rms_warn_arcsec)
            thr_max = float(pe_max_warn_arcsec)
            warn_rms = (np.isfinite(pe_rms0) and pe_rms0 > thr_rms)
            warn_max = (np.isfinite(pe_max0) and pe_max0 > thr_max)
            tag = "[WARN]" if (warn_rms or warn_max) else ""
            print("[scan {:04d}] stream={} key='{}' N={} PE_RMS0={:.2f} arcsec (thr {:.2f}) PE_MAX0={:.2f} arcsec (thr {:.2f}) {} t=[{:.3f},{:.3f}]".format(
                int(scanid), stream_data.stream.name, key_label, int(ii.size), pe_rms0, thr_rms, pe_max0, thr_max, tag, t0, t1
            ))
            print("          PE(mean±std)={:.2f}±{:.2f} arcsec  (dAz*cosEl mean={:.2f}, dEl mean={:.2f} arcsec)".format(
                pe_mean, pe_std, pe_x_mean, pe_y_mean
            ))

        exp0 = _estimate_exposure_s(t_spec[ii], default=0.1)

        if collapse:
            groups = []
            start = 0
            while start < ii.size:
                stop = start + 1
                mode0 = str(mode_str[int(ii[start])])
                while stop < ii.size and str(mode_str[int(ii[stop])]) == mode0:
                    stop += 1
                jj = ii[start:stop]
                groups.append((float(np.nanmin(t_spec[jj])), mode0, jj))
                start = stop
            for (_, mode, jj) in groups:
                rep = int(jj[np.nanargmin(t_spec[jj])])
                rows.append({
                    "stream": stream_data.stream,
                    "i": rep,
                    "time_unix": float(t_spec[rep]),
                    "scanid": int(scanid),
                    "subscan": int(subscan),
                    "intgrp": int(scanid),
                    "obsmode": str(mode),
                    "data": np.nanmean(spec2d[jj, :].astype(np.float64), axis=0).astype(np.float32),
                    "exposure_s": float(exp0) * float(len(jj)),
                    "ra_src": ra_src,
                    "dec_src": dec_src,
                    "ra_deg": float(radec["ra_deg"][rep]),
                    "dec_deg": float(radec["dec_deg"][rep]),
                    "glon_deg": float(radec["glon_deg"][rep]),
                    "glat_deg": float(radec["glat_deg"][rep]),
                    "beam_az_deg": float(pointing["beam_az_deg"][rep]),
                    "beam_el_deg": float(pointing["beam_el_deg"][rep]),
                    "boresight_az_deg": float(pointing["boresight_az"][rep]),
                    "boresight_el_deg": float(pointing["boresight_el"][rep]),
                    "beam_xoff_arcsec": float(pointing["beam_dx_arcsec"][rep]),
                    "beam_yoff_arcsec": float(pointing["beam_dy_arcsec"][rep]),
                    "beam_rot_deg": float(pointing["beam_rot_deg"][rep]),
                    "az_cmd_deg": float(pointing["az_cmd_t"][rep]),
                    "el_cmd_deg": float(pointing["el_cmd_t"][rep]),
                    "az_meas_deg": float(pointing["az_enc_t"][rep]),
                    "el_meas_deg": float(pointing["el_enc_t"][rep]),
                    "corr_az_deg": float(pointing["dlon_t"][rep]),
                    "corr_el_deg": float(pointing["dlat_t"][rep]),
                    "calc_refr_deg": float(calc_refr_deg[rep]) if calc_refr_deg is not None else np.nan,
                    "tamb_c": float(temp_c[rep]) if temp_c is not None else np.nan,
                    "pressure_hpa": float(press_hpa[rep]) if press_hpa is not None else np.nan,
                    "humidity_pct": float(humid_pct[rep]) if humid_pct is not None else np.nan,
                    "pe_x_arcsec": float(pointing["pe_x_arcsec"][rep]),
                    "pe_y_arcsec": float(pointing["pe_y_arcsec"][rep]),
                    "pe_r_arcsec": float(pointing["pe_r_arcsec"][rep]),
                })
                subscan += 1
        else:
            for i in ii:
                i = int(i)
                rows.append({
                    "stream": stream_data.stream,
                    "i": i,
                    "time_unix": float(t_spec[i]),
                    "scanid": int(scanid),
                    "subscan": int(subscan),
                    "intgrp": int(scanid),
                    "obsmode": str(mode_str[i]),
                    "data": np.asarray(spec2d[i, :], dtype=np.float32),
                    "exposure_s": float(exp0),
                    "ra_src": ra_src,
                    "dec_src": dec_src,
                    "ra_deg": float(radec["ra_deg"][i]),
                    "dec_deg": float(radec["dec_deg"][i]),
                    "glon_deg": float(radec["glon_deg"][i]),
                    "glat_deg": float(radec["glat_deg"][i]),
                    "beam_az_deg": float(pointing["beam_az_deg"][i]),
                    "beam_el_deg": float(pointing["beam_el_deg"][i]),
                    "boresight_az_deg": float(pointing["boresight_az"][i]),
                    "boresight_el_deg": float(pointing["boresight_el"][i]),
                    "beam_xoff_arcsec": float(pointing["beam_dx_arcsec"][i]),
                    "beam_yoff_arcsec": float(pointing["beam_dy_arcsec"][i]),
                    "beam_rot_deg": float(pointing["beam_rot_deg"][i]),
                    "az_cmd_deg": float(pointing["az_cmd_t"][i]),
                    "el_cmd_deg": float(pointing["el_cmd_t"][i]),
                    "az_meas_deg": float(pointing["az_enc_t"][i]),
                    "el_meas_deg": float(pointing["el_enc_t"][i]),
                    "corr_az_deg": float(pointing["dlon_t"][i]),
                    "corr_el_deg": float(pointing["dlat_t"][i]),
                    "calc_refr_deg": float(calc_refr_deg[i]) if calc_refr_deg is not None else np.nan,
                    "tamb_c": float(temp_c[i]) if temp_c is not None else np.nan,
                    "pressure_hpa": float(press_hpa[i]) if press_hpa is not None else np.nan,
                    "humidity_pct": float(humid_pct[i]) if humid_pct is not None else np.nan,
                    "pe_x_arcsec": float(pointing["pe_x_arcsec"][i]),
                    "pe_y_arcsec": float(pointing["pe_y_arcsec"][i]),
                    "pe_r_arcsec": float(pointing["pe_r_arcsec"][i]),
                })
                subscan += 1
    return rows


# -----------------------------------------------------------------------------
# 4) Writer bridge
# -----------------------------------------------------------------------------
def make_writer(site, telescope, observer, project, object_name, config_info, streams):
    nchans = sorted(set(int(s.wcs.nchan) for s in streams))
    max_nchan = int(max(nchans))
    cunit1s = sorted(set(str(s.wcs.cunit1).strip() for s in streams))
    if len(cunit1s) != 1:
        raise RuntimeError("All streams must share one CUNIT1 for a single output dataset; got {}".format(cunit1s))

    first = streams[0].wcs
    axis = SpectralAxisUniform(
        crval1_hz=float(first.crval1_hz),
        cdelt1_hz=float(first.cdelt1_hz),
        crpix1=float(first.crpix1),
        restfreq_hz=float(first.restfreq_hz) if np.isfinite(first.restfreq_hz) else 0.0,
        specsys=str(first.specsys),
        ssysobs="TOPOCENT",
        veldef=("RADIO" if first.veldef is None else str(first.veldef)),
        ctype1=str(first.ctype1),
        cunit1=str(first.cunit1),
        refchan=1,
    )

    backend_values = sorted(set(str(s.backend).strip() for s in streams if _is_meaningful_str(s.backend)))
    backend_meta = backend_values[0] if len(backend_values) == 1 else "MIXED"

    shared_meta = {
        "BACKEND": backend_meta,
        "NOTE": "Converted from NECST v4 RawData (multi-stream, beam-center Az/El, command/encoder separated)",
        "CONFIG": str(config_info.get("config_name") or ""),
        "N_STREAM": int(len(streams)),
        "NCHAN_MAX": int(max_nchan),
    }
    info = DatasetInfo(
        telescope=str(telescope),
        observer=str(observer),
        project=str(project),
        object_name=str(object_name),
        radesys="ICRS",
        equinox=2000.0,
        src_radesys="ICRS",
        src_equinox=2000.0,
        eff=Efficiency(effstat="UNKNOWN"),
        refr_included_in_corr=False,
        spectral_axis=axis,
        shared_meta=shared_meta,
    )

    need_freq_column = False
    for s in streams:
        sf = str(s.wcs.specsys).upper()
        store = s.wcs.store_freq_column
        if sf == "LSRK":
            need_freq_column = True
            break
        if isinstance(store, str):
            if store.strip().lower() == "true":
                need_freq_column = True
                break
        elif bool(store):
            need_freq_column = True
            break

    return SDRadioSpectralSDFITSWriter(
        n_chan=int(max_nchan),
        site=site,
        info=info,
        store_freq_column=need_freq_column,
    )


def write_sdfits(out_fits, writer, rows):
    last_key = None
    for row in rows:
        t_unix = float(row["time_unix"])
        mjd = float(Time(t_unix, format="unix", scale="utc").mjd)
        key = (t_unix, int(row.get("scanid_out", row["scanid"])), int(row["subscan"]), int(row["stream"].stream_index))
        if last_key is not None and key < last_key:
            raise RuntimeError("Row order violated: {} < {}".format(key, last_key))
        last_key = key

        sw = row["stream"].wcs
        writer.add_row(
            time_mjd=mjd,
            scanid=int(row.get("scanid_out", row["scanid"])),
            subscan=int(row["subscan"]),
            intgrp=int(row.get("intgrp_out", row["intgrp"])),
            obsmode=str(row["obsmode"]),
            data=np.asarray(row["data"], dtype=np.float32),
            exposure_s=float(row["exposure_s"]),
            polariza=str(row["stream"].polariza),
            fdnum=int(row["stream"].fdnum),
            ifnum=int(row["stream"].ifnum),
            plnum=int(row["stream"].plnum),
            backend=row["stream"].backend,
            sampler=row["stream"].sampler,
            frontend=row["stream"].frontend,
            obsfreq_hz=sw.obsfreq_hz,
            imagfreq_hz=sw.imagfreq_hz,
            lo1freq_hz=sw.lo1_hz,
            lo2freq_hz=sw.lo2_hz,
            lo3freq_hz=sw.lo3_hz,
            sideband=sw.sideband,
            sb1=sw.sb1,
            sb2=sw.sb2,
            sb3=sw.sb3,
            object_name=None,
            calstat="RAW",
            ra_deg=float(row["ra_deg"]),
            dec_deg=float(row["dec_deg"]),
            glon_deg=float(row["glon_deg"]),
            glat_deg=float(row["glat_deg"]),
            srcframe=None,
            src_radesys=None,
            src_equinox=None,
            src_long_deg=None,
            src_lat_deg=None,
            scanframe=None,
            scan_radesys=None,
            scan_equinox=None,
            scan_x_deg=None,
            scan_y_deg=None,
            az_center_deg=float(row["beam_az_deg"]),
            el_center_deg=float(row["beam_el_deg"]),
            az_enc_deg=float(row["az_meas_deg"]),
            el_enc_deg=float(row["el_meas_deg"]),
            boresight_az_deg=float(row["boresight_az_deg"]),
            boresight_el_deg=float(row["boresight_el_deg"]),
            az_cmd_deg=float(row["az_cmd_deg"]),
            el_cmd_deg=float(row["el_cmd_deg"]),
            corr_az_deg=float(row["corr_az_deg"]),
            corr_el_deg=float(row["corr_el_deg"]),
            calc_refr_deg=float(row["calc_refr_deg"]) if np.isfinite(row["calc_refr_deg"]) else np.nan,
            beam_xoff_arcsec=float(row["beam_xoff_arcsec"]),
            beam_yoff_arcsec=float(row["beam_yoff_arcsec"]),
            beam_rot_deg=float(row["beam_rot_deg"]),
            tamb_c=float(row["tamb_c"]) if np.isfinite(row["tamb_c"]) else np.nan,
            pressure_hpa=float(row["pressure_hpa"]) if np.isfinite(row["pressure_hpa"]) else np.nan,
            humidity_pct=float(row["humidity_pct"]) if np.isfinite(row["humidity_pct"]) else np.nan,
            restfreq_hz=(float(sw.restfreq_hz) if np.isfinite(sw.restfreq_hz) else None),
            crval1_hz=float(sw.crval1_hz),
            cdelt1_hz=float(sw.cdelt1_hz),
            crpix1=float(sw.crpix1),
            ctype1=str(sw.ctype1),
            cunit1=str(sw.cunit1),
            specsys=str(sw.specsys),
            veldef=(None if sw.veldef is None else str(sw.veldef)),
        )
    writer.write(out_fits, overwrite=True)


# -----------------------------------------------------------------------------
# 5) Diagnostics
# -----------------------------------------------------------------------------
def _global_pe_summary(rows):
    print("[summary] Global PE over all written samples in each mode (0-centered):")
    by_mode = {}
    for row in rows:
        by_mode.setdefault(str(row["obsmode"]).upper(), []).append(row)
    for mode in ("HOT", "ON", "OFF"):
        rr = by_mode.get(mode, [])
        if not rr:
            print("  {:>3s}: N=0".format(mode))
            continue
        px = np.array([x["pe_x_arcsec"] for x in rr], dtype=float)
        py = np.array([x["pe_y_arcsec"] for x in rr], dtype=float)
        pe2 = px**2 + py**2
        pe2 = pe2[np.isfinite(pe2)]
        if pe2.size == 0:
            print("  {:>3s}: N={}  PE_RMS0=nan  PE_MAX0=nan".format(mode, len(rr)))
            continue
        rms0 = float(np.sqrt(np.nanmean(pe2)))
        max0 = float(np.sqrt(np.nanmax(pe2)))
        print("  {:>3s}: N={}  PE_RMS0={:.2f} arcsec  PE_MAX0={:.2f} arcsec".format(mode, len(rr), rms0, max0))


# -----------------------------------------------------------------------------
# 6) CLI / main
# -----------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(
        prog=pathlib.Path(__file__).name,
        description="Convert NECST v4 RawData spectral DB to SDFITS (multi-stream capable).",
    )
    p.add_argument("rawdata", help="RawData folder path (necstdb + nercst)")
    p.add_argument("--spectrometer-config", default=None, help="TOML file describing spectrometers / beams / LO / WCS")
    p.add_argument("--strict-config", action="store_true", help="Treat stream mismatches as errors.")
    p.add_argument("--spectral", default="xffts-board1", help="Legacy single-stream spectral name")
    p.add_argument("--telescope", default="OMU1P85M", help="Telescope name used in table names")
    p.add_argument("--db-namespace", default=None, help="Database namespace/prefix used in table names (default: necst or TOML global.db_namespace)")
    p.add_argument("--tel-loaddata", default="OMU1p85m", help="legacy compatibility; unused")
    p.add_argument("--out", default=None, help="Output FITS name")
    p.add_argument("--object", default=None, help="OBJECT name (default: rawdata basename)")
    p.add_argument("--project", default="ProjectID", help="PROJID")
    p.add_argument("--observer", default="Unknown", help="OBSERVER")

    p.add_argument("--site-lat", type=float, default=35.940874, help="Site latitude [deg]")
    p.add_argument("--site-lon", type=float, default=138.472153, help="Site longitude [deg, east+]")
    p.add_argument("--site-elev", type=float, default=1386.0, help="Site elevation [m]")

    # legacy single-stream spectral axis defaults
    p.add_argument("--nchan", type=int, default=2**15, help="Legacy: number of channels")
    p.add_argument("--if0-ghz", type=float, default=0.0, help="Legacy: IF start [GHz]")
    p.add_argument("--if1-ghz", type=float, default=2.5, help="Legacy: IF end [GHz]")
    p.add_argument("--lo1-ghz", type=float, default=109.8, help="Legacy: LO1 [GHz]")
    p.add_argument("--lo2-ghz", type=float, default=4.0, help="Legacy: LO2 [GHz]")
    p.add_argument("--restfreq-hz", type=float, default=115.271e9, help="Legacy: RESTFREQ [Hz]")
    p.add_argument("--channel-slice", default=None, help="Optional channel slice applied at conversion time, e.g. '[1024,8192)' or '[1024,8191]'")

    p.add_argument("--encoder-time-col", default="time", help="encoder time column (recommended: time)")
    p.add_argument("--encoder-shift-sec", type=float, default=0.0, help="Shift applied to encoder timestamps before interpolation to spectral time: t_enc <- t_enc + shift [s]. Use the same sign convention as sunscan --encoder-shift-sec.")
    p.add_argument("--altaz-time-col", default="time", help="altaz time column (recommended: time)")
    p.add_argument("--interp-extrap", default="hold", choices=["nan", "hold"], help="Extrapolation policy for interpolation")
    p.add_argument("--use-modes", default="ON,HOT,OFF", help="Comma-separated OBSMODE list to include")
    p.add_argument("--include-unknown", action="store_true", help="Also include UNKNOWN samples")
    p.add_argument("--pe-rms-warn-arcsec", type=float, default=15.0)
    p.add_argument("--pe-max-warn-arcsec", type=float, default=30.0)
    p.add_argument("--scan-key", default="id_mode", choices=["id", "id_mode"])
    p.add_argument("--skip-nonstandard-position", action="store_true")
    p.add_argument("--wcs-table", default="necst-telescope-coordinate-wcs", help="WCS table name for RA/DEC time series")

    p.add_argument("--radec-method", default="wcs_first", choices=["wcs_first", "azel"])
    p.add_argument(
        "--radec-azel-source",
        default="beam",
        choices=["beam", "true", "encoder", "altaz", "cmd"],
        help=(
            "Az/El source used when RA/DEC is derived from Az/El (wcs_first fallback or --radec-method azel). "
            "beam=current beam-center; true=encoder-corrected boresight; encoder=measured encoder with beam offset; "
            "altaz=raw altaz/cmd with beam offset; cmd=(altaz/cmd - dlon/dlat) with beam offset."
        ),
    )
    p.add_argument("--refraction", default="on", choices=["on", "off"])
    p.add_argument("--met-pressure-hpa", type=float, default=None)
    p.add_argument("--met-temperature-c", type=float, default=10.0)
    p.add_argument("--met-humidity-pct", type=float, default=50.0)
    p.add_argument("--met-obswl-um", type=float, default=None)

    p.add_argument("--met-source", default="auto", choices=["auto", "weather", "spectral", "fallback"])
    p.add_argument("--encoder-table", default=None)
    p.add_argument("--encoder-table-suffix", default=None)
    p.add_argument("--altaz-table", default=None)
    p.add_argument("--altaz-table-suffix", default=None)
    p.add_argument("--weather-inside-table", default=None)
    p.add_argument("--weather-inside-table-suffix", default=None)
    p.add_argument("--weather-inside-time-col", default=None)
    p.add_argument("--weather-outside-table", default=None)
    p.add_argument("--weather-outside-table-suffix", default=None)
    p.add_argument("--weather-outside-time-col", default=None)
    p.add_argument("--weather-table", default="auto")
    p.add_argument("--weather-time-col", default="time")
    p.add_argument("--inside-default-temperature-c", type=float, default=None)
    p.add_argument("--inside-temperature-min-c", type=float, default=None)
    p.add_argument("--inside-temperature-max-c", type=float, default=None)
    p.add_argument("--outside-default-temperature-c", type=float, default=None)
    p.add_argument("--outside-default-pressure-hpa", type=float, default=None)
    p.add_argument("--outside-default-humidity-pct", type=float, default=None)
    p.add_argument("--outside-temperature-min-c", type=float, default=None)
    p.add_argument("--outside-temperature-max-c", type=float, default=None)
    p.add_argument("--outside-pressure-min-hpa", type=float, default=None)
    p.add_argument("--outside-pressure-max-hpa", type=float, default=None)
    p.add_argument("--outside-humidity-min-pct", type=float, default=None)
    p.add_argument("--outside-humidity-max-pct", type=float, default=None)

    p.add_argument("--collapse", action="store_true")
    p.add_argument("--pe-cor", action="store_true", help="legacy compatibility; unused")
    p.add_argument("--plot-pe", action="store_true", help="currently not supported for multi-stream; reserved")
    p.add_argument("--print-pe", action="store_true")
    p.add_argument("--plot-pe-outdir", default=".")
    p.add_argument("--plot-pe-max-scans", type=int, default=0)
    return p.parse_args(argv)


def _default_outfile(rawbase, config_dict, args):
    if args.out:
        return args.out
    if args.spectrometer_config:
        cfg_name = config_dict.get("config_name")
        if cfg_name:
            return "{}_{}_sdfits.fits".format(rawbase, _sanitize_for_filename(cfg_name))
        return "{}_sdfits.fits".format(rawbase)
    return "{}_{}_sdfits.fits".format(rawbase, str(args.spectral).replace("/", "_"))


def _argv_has_option(argv, *opts):
    sargv = [str(x) for x in (argv or [])]
    for a in sargv:
        for opt in opts:
            if a == opt or a.startswith(opt + "="):
                return True
    return False


def _resolve_runtime_naming(args, config_dict, argv=None):
    global_cfg = dict(config_dict.get("global", {}) or {})
    db_namespace = _nonempty_str(getattr(args, "db_namespace", None), default=None)
    if db_namespace is None:
        db_namespace = _nonempty_str(global_cfg.get("db_namespace"), default="necst")

    telescope_cfg = _nonempty_str(global_cfg.get("telescope"), default=None)
    telescope = telescope_cfg if telescope_cfg is not None else str(args.telescope)
    # Prefer explicit CLI telescope over TOML.  We detect whether the option was
    # present in argv so that even `--telescope OMU1P85M` overrides a TOML value.
    if (not args.spectrometer_config) or _argv_has_option(argv, "--telescope") or (telescope_cfg is None):
        telescope = str(args.telescope)

    encoder_table = _nonempty_str(getattr(args, "encoder_table", None), default=None)
    if encoder_table is None:
        encoder_table = _nonempty_str(global_cfg.get("encoder_table"), default=None)
    encoder_table_suffix = _nonempty_str(getattr(args, "encoder_table_suffix", None), default=_nonempty_str(global_cfg.get("encoder_table_suffix"), default="ctrl-antenna-encoder"))
    encoder_time_col_cfg = _nonempty_str(global_cfg.get("encoder_time_col"), default=None)
    if (not args.spectrometer_config) or _argv_has_option(argv, "--encoder-time-col") or (encoder_time_col_cfg is None):
        encoder_time_col = str(getattr(args, "encoder_time_col", "time"))
    else:
        encoder_time_col = str(encoder_time_col_cfg)
    altaz_table = _nonempty_str(getattr(args, "altaz_table", None), default=None)
    if altaz_table is None:
        altaz_table = _nonempty_str(global_cfg.get("altaz_table"), default=None)
    altaz_table_suffix = _nonempty_str(getattr(args, "altaz_table_suffix", None), default=_nonempty_str(global_cfg.get("altaz_table_suffix"), default="ctrl-antenna-altaz"))
    altaz_time_col_cfg = _nonempty_str(global_cfg.get("altaz_time_col"), default=None)
    if (not args.spectrometer_config) or _argv_has_option(argv, "--altaz-time-col") or (altaz_time_col_cfg is None):
        altaz_time_col = str(getattr(args, "altaz_time_col", "time"))
    else:
        altaz_time_col = str(altaz_time_col_cfg)
    weather_table_cfg = _nonempty_str(global_cfg.get("weather_table"), default=None)
    legacy_weather_time_col_cfg = _nonempty_str(global_cfg.get("weather_time_col"), default=None)
    weather_inside_table = _nonempty_str(getattr(args, "weather_inside_table", None), default=_nonempty_str(global_cfg.get("weather_inside_table"), default=None))
    weather_inside_table_suffix = _nonempty_str(getattr(args, "weather_inside_table_suffix", None), default=_nonempty_str(global_cfg.get("weather_inside_table_suffix"), default="weather-ambient"))
    weather_outside_table = _nonempty_str(getattr(args, "weather_outside_table", None), default=_nonempty_str(global_cfg.get("weather_outside_table"), default=None))
    weather_outside_table_suffix = _nonempty_str(getattr(args, "weather_outside_table_suffix", None), default=_nonempty_str(global_cfg.get("weather_outside_table_suffix"), default="weather-ambient"))

    weather_inside_time_col = _nonempty_str(getattr(args, "weather_inside_time_col", None), default=_nonempty_str(global_cfg.get("weather_inside_time_col"), default=legacy_weather_time_col_cfg))
    weather_outside_time_col = _nonempty_str(getattr(args, "weather_outside_time_col", None), default=_nonempty_str(global_cfg.get("weather_outside_time_col"), default=legacy_weather_time_col_cfg))

    def _choose_float_runtime(cli_attr, global_key, default):
        cli_val = getattr(args, cli_attr, None)
        if cli_val is not None:
            return float(cli_val)
        gval = global_cfg.get(global_key, None)
        if gval is not None:
            return float(gval)
        return float(default)

    inside_default_temperature_c = _choose_float_runtime("inside_default_temperature_c", "inside_default_temperature_c", getattr(args, "met_temperature_c", 10.0))
    inside_temperature_min_c = _choose_float_runtime("inside_temperature_min_c", "inside_temperature_min_c", -50.0)
    inside_temperature_max_c = _choose_float_runtime("inside_temperature_max_c", "inside_temperature_max_c", 50.0)
    outside_default_temperature_c = _choose_float_runtime("outside_default_temperature_c", "outside_default_temperature_c", getattr(args, "met_temperature_c", 10.0))
    outside_default_pressure_hpa = _choose_float_runtime("outside_default_pressure_hpa", "outside_default_pressure_hpa", getattr(args, "met_pressure_hpa", None) if getattr(args, "met_pressure_hpa", None) is not None else _pressure_from_elev_hpa(float(args.site_elev)))
    outside_default_humidity_pct = _choose_float_runtime("outside_default_humidity_pct", "outside_default_humidity_pct", getattr(args, "met_humidity_pct", 50.0))
    outside_temperature_min_c = _choose_float_runtime("outside_temperature_min_c", "outside_temperature_min_c", -50.0)
    outside_temperature_max_c = _choose_float_runtime("outside_temperature_max_c", "outside_temperature_max_c", 50.0)
    outside_pressure_min_hpa = _choose_float_runtime("outside_pressure_min_hpa", "outside_pressure_min_hpa", 400.0)
    outside_pressure_max_hpa = _choose_float_runtime("outside_pressure_max_hpa", "outside_pressure_max_hpa", 1100.0)
    outside_humidity_min_pct = _choose_float_runtime("outside_humidity_min_pct", "outside_humidity_min_pct", 0.0)
    outside_humidity_max_pct = _choose_float_runtime("outside_humidity_max_pct", "outside_humidity_max_pct", 100.0)

    encoder_shift_sec_cfg = global_cfg.get("encoder_shift_sec", None)
    if (not args.spectrometer_config) or _argv_has_option(argv, "--encoder-shift-sec") or (encoder_shift_sec_cfg is None):
        encoder_shift_sec = float(getattr(args, "encoder_shift_sec", 0.0))
    else:
        encoder_shift_sec = float(encoder_shift_sec_cfg)

    return {
        "db_namespace": db_namespace,
        "telescope": telescope,
        "encoder_table": encoder_table,
        "encoder_table_suffix": encoder_table_suffix,
        "encoder_time_col": encoder_time_col,
        "altaz_table": altaz_table,
        "altaz_table_suffix": altaz_table_suffix,
        "altaz_time_col": altaz_time_col,
        "weather_table_cfg": weather_table_cfg,
        "weather_inside_table": weather_inside_table,
        "weather_inside_table_suffix": weather_inside_table_suffix,
        "weather_outside_table": weather_outside_table,
        "weather_outside_table_suffix": weather_outside_table_suffix,
        "weather_inside_time_col": weather_inside_time_col,
        "weather_outside_time_col": weather_outside_time_col,
        "inside_default_temperature_c": inside_default_temperature_c,
        "inside_temperature_min_c": inside_temperature_min_c,
        "inside_temperature_max_c": inside_temperature_max_c,
        "outside_default_temperature_c": outside_default_temperature_c,
        "outside_default_pressure_hpa": outside_default_pressure_hpa,
        "outside_default_humidity_pct": outside_default_humidity_pct,
        "outside_temperature_min_c": outside_temperature_min_c,
        "outside_temperature_max_c": outside_temperature_max_c,
        "outside_pressure_min_hpa": outside_pressure_min_hpa,
        "outside_pressure_max_hpa": outside_pressure_max_hpa,
        "outside_humidity_min_pct": outside_humidity_min_pct,
        "outside_humidity_max_pct": outside_humidity_max_pct,
        "encoder_shift_sec": encoder_shift_sec,
    }


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = parse_args(argv)

    for attr in ("encoder_time_col", "altaz_time_col", "weather_time_col", "weather_inside_time_col", "weather_outside_time_col"):
        cur = getattr(args, attr, None)
        if cur is None:
            continue
        if str(cur).strip().lower() == "recorded_time":
            print("[warn] {}=recorded_time is not supported; forcing to 'time'".format(attr))
            setattr(args, attr, "time")

    radec_azel_source_now = str(getattr(args, "radec_azel_source", "beam")).strip().lower()
    if radec_azel_source_now not in ("beam", "true", "encoder", "altaz", "cmd"):
        raise ValueError("unsupported --radec-azel-source={!r}".format(radec_azel_source_now))
    if radec_azel_source_now != "beam":
        print("[info] using --radec-azel-source={} for Az/El-derived RA/DEC".format(radec_azel_source_now))
    if bool(getattr(args, "strict_config", False)):
        print("[warn] --strict-config is reserved in v1_40 and is not yet enforced beyond current validation.")

    rawdata_path = pathlib.Path(args.rawdata).expanduser().resolve()
    rawbase = rawdata_path.name
    object_name = args.object if args.object else rawbase

    if args.spectrometer_config:
        config_dict = load_spectrometer_config(args.spectrometer_config)
    else:
        config_dict = build_legacy_single_stream_config(args)

    global_cfg = dict(config_dict.get("global", {}) or {})
    runtime_naming = _resolve_runtime_naming(args, config_dict, argv=argv)
    db_namespace = runtime_naming["db_namespace"]
    telescope_name = runtime_naming["telescope"]
    args.encoder_shift_sec = float(runtime_naming.get("encoder_shift_sec", getattr(args, "encoder_shift_sec", 0.0)))
    args.encoder_time_col = str(runtime_naming.get("encoder_time_col", getattr(args, "encoder_time_col", "time")))
    args.altaz_time_col = str(runtime_naming.get("altaz_time_col", getattr(args, "altaz_time_col", "time")))
    for attr in ("encoder_time_col", "altaz_time_col", "weather_inside_time_col", "weather_outside_time_col"):
        cur = runtime_naming.get(attr, getattr(args, attr, None))
        if cur is None:
            continue
        if str(cur).strip().lower() == "recorded_time":
            label = attr if attr not in ("weather_inside_time_col", "weather_outside_time_col") else "weather_time_col"
            print("[warn] {}=recorded_time is not supported; forcing to 'time'".format(label))
            runtime_naming[attr] = "time"
    args.encoder_time_col = str(runtime_naming.get("encoder_time_col", getattr(args, "encoder_time_col", "time")))
    args.altaz_time_col = str(runtime_naming.get("altaz_time_col", getattr(args, "altaz_time_col", "time")))
    if float(args.encoder_shift_sec) != 0.0:
        print("[info] encoder_shift_sec = {} s (applied as t_enc <- t_enc + shift before interpolation to spectral time)".format(float(args.encoder_shift_sec)))
    apply_channel_slice_config(config_dict["streams"], global_cfg, cli_channel_slice=getattr(args, "channel_slice", None))
    layout = str(global_cfg.get("output_layout", "merged")).strip().lower()
    if layout not in ("merged", "merged_time"):
        raise RuntimeError("Only merged/merged_time output_layout is supported in v1_40")
    if layout in ("merged", "merged_time") and not _as_bool(global_cfg.get("time_sort", True), True):
        raise RuntimeError("merged/merged_time output_layout requires global.time_sort=true in v1_40")
    streams = list(config_dict["streams"])
    out_fits = _default_outfile(rawbase, config_dict, args)

    w_table_arg = str(getattr(args, "weather_table", "auto")).strip()
    encoder_table_resolved = _resolve_prefixed_table_name(db_namespace, telescope_name, runtime_naming.get("encoder_table"), runtime_naming.get("encoder_table_suffix"), "ctrl-antenna-encoder")
    altaz_table_resolved = _resolve_prefixed_table_name(db_namespace, telescope_name, runtime_naming.get("altaz_table"), runtime_naming.get("altaz_table_suffix"), "ctrl-antenna-altaz")
    weather_inside_table = _nonempty_str(runtime_naming.get("weather_inside_table"), default=None)
    weather_outside_table = _nonempty_str(runtime_naming.get("weather_outside_table"), default=None)
    legacy_weather_auto = _nonempty_str(runtime_naming.get("weather_table_cfg"), default=None)
    legacy_weather_table_cli = _argv_has_option(argv, "--weather-table") and (w_table_arg.lower() != "auto")
    weather_inside_table_cli = _argv_has_option(argv, "--weather-inside-table")
    weather_outside_table_cli = _argv_has_option(argv, "--weather-outside-table")
    if legacy_weather_table_cli:
        legacy_weather_auto = w_table_arg
        if not weather_inside_table_cli:
            weather_inside_table = w_table_arg
        if not weather_outside_table_cli:
            weather_outside_table = w_table_arg
    weather_inside_table_resolved = _resolve_prefixed_table_name(db_namespace, telescope_name, weather_inside_table or legacy_weather_auto, runtime_naming.get("weather_inside_table_suffix"), "weather-ambient")
    weather_outside_table_resolved = _resolve_prefixed_table_name(db_namespace, telescope_name, weather_outside_table or legacy_weather_auto, runtime_naming.get("weather_outside_table_suffix"), "weather-ambient")
    legacy_weather_time_cli = _argv_has_option(argv, "--weather-time-col")
    weather_inside_time_col_cli = _argv_has_option(argv, "--weather-inside-time-col")
    weather_outside_time_col_cli = _argv_has_option(argv, "--weather-outside-time-col")
    legacy_weather_time = _nonempty_str(getattr(args, "weather_time_col", None), default="time")
    weather_inside_time_col = _nonempty_str(runtime_naming.get("weather_inside_time_col"), default=legacy_weather_time)
    weather_outside_time_col = _nonempty_str(runtime_naming.get("weather_outside_time_col"), default=legacy_weather_time)
    if legacy_weather_time_cli:
        if not weather_inside_time_col_cli:
            weather_inside_time_col = legacy_weather_time
        if not weather_outside_time_col_cli:
            weather_outside_time_col = legacy_weather_time
    args.inside_default_temperature_c = float(runtime_naming.get("inside_default_temperature_c"))
    args.inside_temperature_min_c = float(runtime_naming.get("inside_temperature_min_c"))
    args.inside_temperature_max_c = float(runtime_naming.get("inside_temperature_max_c"))
    args.outside_default_temperature_c = float(runtime_naming.get("outside_default_temperature_c"))
    args.outside_default_pressure_hpa = float(runtime_naming.get("outside_default_pressure_hpa"))
    args.outside_default_humidity_pct = float(runtime_naming.get("outside_default_humidity_pct"))
    args.outside_temperature_min_c = float(runtime_naming.get("outside_temperature_min_c"))
    args.outside_temperature_max_c = float(runtime_naming.get("outside_temperature_max_c"))
    args.outside_pressure_min_hpa = float(runtime_naming.get("outside_pressure_min_hpa"))
    args.outside_pressure_max_hpa = float(runtime_naming.get("outside_pressure_max_hpa"))
    args.outside_humidity_min_pct = float(runtime_naming.get("outside_humidity_min_pct"))
    args.outside_humidity_max_pct = float(runtime_naming.get("outside_humidity_max_pct"))
    common = extract_common_inputs(
        rawdata_path=rawdata_path,
        db_namespace=str(db_namespace),
        telescope=str(telescope_name),
        wcs_table=str(args.wcs_table),
        weather_inside_table=weather_inside_table_resolved,
        weather_inside_time_col=str(weather_inside_time_col),
        weather_outside_table=weather_outside_table_resolved,
        weather_outside_time_col=str(weather_outside_time_col),
        encoder_table=encoder_table_resolved,
        altaz_table=altaz_table_resolved,
    )

    site = Site(lat_deg=float(args.site_lat), lon_deg=float(args.site_lon), elev_m=float(args.site_elev))

    writer = make_writer(
        site=site,
        telescope=str(telescope_name),
        observer=str(args.observer),
        project=str(args.project),
        object_name=str(object_name),
        config_info=config_dict,
        streams=streams,
    )
    writer.add_history("converter", pathlib.Path(__file__).name)
    writer.add_history("config_name", config_dict.get("config_name", ""))
    if args.spectrometer_config:
        writer.add_history("config_path", str(pathlib.Path(args.spectrometer_config).expanduser().resolve()))
        if global_cfg.get("telescope") is not None:
            writer.add_history("config_global_telescope", str(global_cfg.get("telescope")))
        if global_cfg.get("db_namespace") is not None:
            writer.add_history("config_global_db_namespace", str(global_cfg.get("db_namespace")))
    writer.add_history("output_layout", layout)
    writer.add_history("db_namespace", str(db_namespace))
    writer.add_history("telescope_tables", str(telescope_name))
    if getattr(args, "channel_slice", None) is not None:
        writer.add_history("channel_slice_cli", str(args.channel_slice))
    elif global_cfg.get("channel_slice", None) is not None:
        writer.add_history("channel_slice_global", str(global_cfg.get("channel_slice")))
    writer.add_history("azimuth_meaning", "AZIMUTH/ELEVATIO=beam-center Az/El")
    writer.add_history("az_enc_meaning", "AZ_ENC/EL_ENC=measured encoder Az/El interpolated to t_spec")
    writer.add_history("encoder_shift_sec", float(getattr(args, "encoder_shift_sec", 0.0)))
    writer.add_history("encoder_time_col", str(getattr(args, "encoder_time_col", "time")))
    writer.add_history("altaz_time_col", str(getattr(args, "altaz_time_col", "time")))
    writer.add_history("bore_meaning", "BORE_AZ/BORE_EL=boresight Az/El = encoder - correction")
    writer.add_history("az_cmd_meaning", "AZ_CMD/EL_CMD=commanded altaz interpolated to t_spec")
    writer.add_history("corr_meaning", "CORR_AZ/CORR_EL=dlon/dlat")
    writer.add_history("beam_model_columns", "BEAMXOFF/BEAMYOFF[arcsec],BEAMROT[deg]=applied beam model")
    writer.add_history("galactic_columns", "GLON/GLAT=galactic longitude/latitude derived from ICRS RA/DEC in converter")
    writer.add_history("encoder_table", str(encoder_table_resolved))
    writer.add_history("altaz_table", str(altaz_table_resolved))
    writer.add_history("weather_inside_table", str(weather_inside_table_resolved))
    writer.add_history("weather_outside_table", str(weather_outside_table_resolved))
    writer.add_history("weather_inside_time_col", str(weather_inside_time_col))
    writer.add_history("weather_outside_time_col", str(weather_outside_time_col))
    writer.add_history("source_block", "omitted unless true source coordinates become available")
    writer.add_history("scan_block", "omitted unless true scan offsets become available")

    all_rows = []
    spectral_time_summaries = []

    for stream in streams:
        print("[info] processing stream '{}'".format(stream.name))
        stream_data = extract_spectral_stream(common, str(db_namespace), str(telescope_name), stream)
        spectral_time_summaries.append({
            "stream": stream.name,
            "applied": stream_data.spec_time_basis,
            "suffix": stream_data.spec_time_suffix,
            "fallback": stream_data.spec_time_fallback_field,
            "example": stream_data.spec_time_example,
            "field": stream_data.spec_time_field,
            "table": stream_data.spec_table_name,
        })
        pointing = normalize_pointing(
            stream_data.t_spec,
            common.arr_enc,
            common.arr_alt,
            str(args.encoder_time_col),
            str(args.altaz_time_col),
            str(args.interp_extrap),
            float(args.encoder_shift_sec),
        )
        beam_applied = apply_beam_offset(pointing["boresight_az"], pointing["boresight_el"], stream.beam)
        pointing.update(beam_applied)
        if np.any(~np.isfinite(pointing["beam_az_deg"])) or np.any(~np.isfinite(pointing["beam_el_deg"])):
            nbad = int(np.count_nonzero(~np.isfinite(pointing["beam_az_deg"])) + np.count_nonzero(~np.isfinite(pointing["beam_el_deg"])))
            raise RuntimeError(
                "Beam-center Az/El contains NaN after applying beam offset/rotation for stream '{}'. Check beam offsets and elevations near the pole.".format(stream.name)
            )

        radec_azel = select_radec_azel_source(pointing, stream.beam, str(args.radec_azel_source))
        if np.any(~np.isfinite(radec_azel["az_deg"])) or np.any(~np.isfinite(radec_azel["el_deg"])):
            nbad = int(np.count_nonzero(~np.isfinite(radec_azel["az_deg"])) + np.count_nonzero(~np.isfinite(radec_azel["el_deg"])))
            raise RuntimeError(
                "Selected RA/DEC Az/El source '{}' contains NaN for stream '{}' ({} elements).".format(
                    radec_azel["source"], stream.name, nbad
                )
            )

        tamb_c, temp_c, press_hpa, humid_pct, inside_weather, outside_weather = resolve_stream_meteorology(
            common, stream_data, str(args.met_source).strip().lower(), str(args.interp_extrap)
        )
        tamb_c_use = finalize_inside_temperature_for_tamb(stream_data.t_spec, tamb_c, args)

        apply_refraction = (str(getattr(args, "refraction", "on")).strip().lower() == "on")
        if apply_refraction:
            temp_use, press_use, humid_use = finalize_meteorology_for_refraction(stream_data.t_spec, temp_c, press_hpa, humid_pct, args)
        else:
            temp_use, press_use, humid_use = temp_c, press_hpa, humid_pct

        obswl_um = _obswl_um_from_stream_or_args(stream.wcs, args)
        radec = normalize_radec(
            wcs_df=common.wcs_df if str(args.radec_method).strip().lower() == "wcs_first" else None,
            t_spec=stream_data.t_spec,
            beam_az_deg=radec_azel["az_deg"],
            beam_el_deg=radec_azel["el_deg"],
            site=site,
            radec_method=str(args.radec_method),
            apply_refraction=apply_refraction,
            press_hpa=press_use,
            temp_c=temp_use,
            humid_pct=humid_use,
            obswl_um=obswl_um,
        )
        radec.update(_gal_from_radec(radec["ra_deg"], radec["dec_deg"]))
        if np.any(~np.isfinite(radec["glon_deg"])) or np.any(~np.isfinite(radec["glat_deg"])):
            nbad = int(np.count_nonzero(~np.isfinite(radec["glon_deg"])) + np.count_nonzero(~np.isfinite(radec["glat_deg"])))
            raise RuntimeError("GLON/GLAT contains NaN after RA/DEC->Galactic conversion for stream '{}' ({} elements).".format(stream.name, nbad))

        if apply_refraction and press_use is not None and temp_use is not None and humid_use is not None:
            calc_refr_deg = _calc_refraction_deg(
                site, stream_data.t_spec, pointing["beam_az_deg"], pointing["beam_el_deg"],
                press_use, temp_use, humid_use, obswl_um
            )
        else:
            calc_refr_deg = np.full_like(stream_data.t_spec, np.nan, dtype=float)

        temp_out = tamb_c_use
        if apply_refraction:
            press_out, humid_out = press_use, humid_use
        else:
            press_out, humid_out = press_hpa, humid_pct

        rows = build_rows_for_stream(
            stream_data=stream_data,
            pointing=pointing,
            radec=radec,
            calc_refr_deg=calc_refr_deg,
            temp_c=temp_out,
            press_hpa=press_out,
            humid_pct=humid_out,
            collapse=bool(args.collapse),
            scan_key=str(args.scan_key),
            skip_nonstandard_position=bool(args.skip_nonstandard_position),
            use_modes_csv=str(args.use_modes),
            include_unknown=bool(args.include_unknown),
            pe_rms_warn_arcsec=float(args.pe_rms_warn_arcsec),
            pe_max_warn_arcsec=float(args.pe_max_warn_arcsec),
            verbose_scan_pe=bool(args.print_pe),
        )
        all_rows.extend(rows)

        writer.add_history(
            "stream_{}".format(stream.stream_index),
            "name={};fdnum={};ifnum={};plnum={};polariza={};beam_id={};rotation_mode={};beam_model_version={};channel_slice={};nchan_full={};nchan_written={};obsfreq_fallback={};spec_time_basis={};spec_time_suffix={};spec_time_fallback={};radec_azel_source={};radec_azel_meaning={}".format(
                stream.name,
                stream.fdnum,
                stream.ifnum,
                stream.plnum,
                stream.polariza,
                stream.beam.beam_id,
                stream.beam.rotation_mode,
                (stream.beam.beam_model_version or ""),
                (stream.channel_slice_label or "full"),
                (int(stream.wcs_full.nchan) if stream.wcs_full is not None else int(stream.wcs.nchan)),
                int(stream.wcs.nchan),
                ("RESTFREQ->OBSFREQ" if (_is_meaningful_scalar(stream.wcs.obsfreq_hz) and not _is_meaningful_scalar(stream.local_oscillators.get("obsfreq_hz", None))) else "explicit_or_none"),
                stream_data.spec_time_basis,
                (stream_data.spec_time_suffix or ""),
                (stream_data.spec_time_fallback_field or ""),
                radec_azel["source"],
                radec_azel["meaning"],
            ),
        )

    all_rows.sort(key=lambda r: (float(r["time_unix"]), int(r["scanid"]), int(r["subscan"]), int(r["stream"].stream_index)))
    scan_map = {}
    next_scanid = 0
    for row in all_rows:
        key = (int(row["stream"].stream_index), int(row["scanid"]))
        if key not in scan_map:
            scan_map[key] = next_scanid
            next_scanid += 1
        row["scanid_out"] = int(scan_map[key])
        row["intgrp_out"] = int(scan_map[key])

    if bool(args.plot_pe):
        print("[warn] --plot-pe is not implemented in v1_40 for multi-stream mode; skipping plot generation.")

    write_sdfits(str(out_fits), writer, all_rows)
    _global_pe_summary(all_rows)
    print("Saved: {}".format(out_fits))
    if spectral_time_summaries:
        print("[info] spectral time basis summary:")
        for ent in spectral_time_summaries:
            print(
                "  stream='{stream}' applied={applied} suffix={suffix} fallback={fallback} field={field} example={example}".format(
                    stream=ent.get("stream"),
                    applied=ent.get("applied"),
                    suffix=ent.get("suffix"),
                    fallback=ent.get("fallback"),
                    field=ent.get("field"),
                    example=ent.get("example"),
                )
            )


if __name__ == "__main__":
    main()
