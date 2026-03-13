# src/sd_radio_spectral_fits/calibrate.py
from __future__ import annotations

import builtins
from typing import Optional, Tuple, Union, Any, Sequence

import numpy as np
import pandas as pd

from .axis import wcs_slice_channels, channel_slice_from_vrange_union
from .doppler import compute_vcorr_series
from .rawspec import load_rawspec_auto

from .fitsio import Scantable, write_scantable
from .scantable_utils import _parse_row_selector, _df_to_native_endian, _resolve_table_timestamps
from .restfreq import apply_restfreq_override

from .utils import validate_mapping_frame

# --- 大気補正用モジュールのインポート ---
import warnings
    # 新しく作成した extract_meta_array もインポートするように修正
from .atmosphere import compute_t_cal_array, estimate_t_atm, extract_meta_value, extract_meta_array


def _maybe_validate_mapping_frame(meta: dict, mapping: pd.DataFrame) -> str | None:
    coord_cols = {str(c).strip().upper() for c in mapping.columns}
    has_any = bool(coord_cols.intersection({"RA", "DEC", "RA_DEG", "DEC_DEG", "GLON", "GLAT", "OBSRA", "OBSDEC"})) or ("coord_frame" in meta)
    if not has_any:
        return None
    return validate_mapping_frame(meta, mapping)

def _resolve_lonlat_columns(mapping: pd.DataFrame, coord_frame: str) -> tuple[str | None, str | None]:
    frame = str(coord_frame or "icrs").strip().lower()
    if frame == "galactic":
        pairs = (("GLON", "GLAT"), ("glon", "glat"))
    else:
        pairs = (("ra_deg", "dec_deg"), ("RA", "DEC"), ("ra", "dec"), ("OBSRA", "OBSDEC"))

    for lon_col, lat_col in pairs:
        if lon_col in mapping.columns and lat_col in mapping.columns:
            return lon_col, lat_col
    return None, None


def _choose_effective_coord_frame(
    mapping: pd.DataFrame,
    coord_frame: str | None,
    meta: dict,
) -> str:
    """Choose the coordinate frame for VELOSYS computation.

    Priority:
      1) explicit coord_frame argument if provided and matching columns exist
      2) meta['coord_frame'] if matching columns exist
      3) infer from available columns (GLON/GLAT -> galactic, else RA/DEC -> icrs)

    Notes
    -----
    The public API default is None (auto). This avoids accidentally forcing
    ICRS when both RA/DEC and GLON/GLAT are present and the metadata declares
    Galactic coordinates.
    """
    candidates = []
    if coord_frame not in (None, "", "auto"):
        candidates.append(str(coord_frame).strip().lower())
    meta_frame = meta.get("coord_frame")
    if meta_frame not in (None, ""):
        candidates.append(str(meta_frame).strip().lower())
    candidates.extend(["galactic", "icrs"])

    seen = set()
    for cand in candidates:
        if cand in seen:
            continue
        seen.add(cand)
        lon_col, lat_col = _resolve_lonlat_columns(mapping, cand)
        if lon_col is not None and lat_col is not None:
            return cand
    return str(meta.get("coord_frame") or coord_frame or "icrs").strip().lower()

def _is_missing_value(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    try:
        return bool(pd.isna(value))
    except Exception:
        return False


def _first_valid_value(meta: dict, table: pd.DataFrame | None, keys: Sequence[str]) -> Any:
    for key in keys:
        if key in meta and not _is_missing_value(meta[key]):
            return meta[key]
    if table is not None and not table.empty:
        for key in keys:
            if key not in table.columns:
                continue
            col = table[key]
            for value in col.tolist():
                if not _is_missing_value(value):
                    return value
    return None


def _normalize_specsys_name(value: Any, *, default: str = "TOPOCENT") -> str:
    if _is_missing_value(value):
        value = default
    s = str(value).strip().upper().replace(" ", "")
    aliases = {
        "TOPO": "TOPOCENT",
        "TOPOCENTER": "TOPOCENT",
        "LSR": "LSRK",
    }
    s = aliases.get(s, s)
    if s not in ("TOPOCENT", "LSRK"):
        raise ValueError(
            f"Unsupported SPECSYS/SSYSOBS={value!r}. "
            "run_tastar_calibration currently supports only frequency-axis TOPOCENT or LSRK input."
        )
    return s


def _frequency_unit_scale_to_hz(unit: Any) -> float:
    if _is_missing_value(unit):
        return 1.0
    s = str(unit).strip().lower().replace(" ", "")
    table = {
        "hz": 1.0,
        "khz": 1.0e3,
        "mhz": 1.0e6,
        "ghz": 1.0e9,
    }
    if s in table:
        return table[s]
    raise ValueError(
        f"Unsupported frequency-axis CUNIT1={unit!r}. "
        "run_tastar_calibration currently supports only frequency input convertible to Hz."
    )


def _normalize_frequency_axis_ctype(value: Any) -> str:
    if _is_missing_value(value):
        return "FREQ"
    s = str(value).strip().upper()
    if s.startswith("FREQ"):
        return "FREQ"
    raise ValueError(
        f"Unsupported CTYPE1={value!r}. "
        "run_tastar_calibration currently supports only frequency-axis input."
    )


def _canonicalize_frequency_axis_meta(meta: dict, table: pd.DataFrame | None = None) -> dict:
    m = dict(meta)

    ctype1 = _normalize_frequency_axis_ctype(_first_valid_value(m, table, ("CTYPE1", "CTYPE3")))
    cunit1_raw = _first_valid_value(m, table, ("CUNIT1", "CUNIT3"))
    cunit_scale = _frequency_unit_scale_to_hz(cunit1_raw)

    crval1 = _first_valid_value(m, table, ("CRVAL1", "CRVAL3"))
    cdelt1 = _first_valid_value(m, table, ("CDELT1", "CDELT3"))
    crpix1 = _first_valid_value(m, table, ("CRPIX1", "CRPIX3"))

    if not _is_missing_value(crval1):
        m["CRVAL1"] = float(crval1) * cunit_scale
    if not _is_missing_value(cdelt1):
        m["CDELT1"] = float(cdelt1) * cunit_scale
    if not _is_missing_value(crpix1):
        m["CRPIX1"] = float(crpix1)

    rest = _first_valid_value(m, table, ("RESTFRQ", "RESTFREQ", "rest_hz", "restfrq_hz", "restfreq_hz"))
    if _is_missing_value(rest):
        raise ValueError(
            "Frequency-axis calibration requires RESTFRQ/RESTFREQ (or rest_freq override) in metadata."
        )
    rest_hz = float(rest)
    if not np.isfinite(rest_hz) or rest_hz <= 0.0:
        raise ValueError(f"Invalid RESTFRQ/RESTFREQ value: {rest!r}")

    specsys = _normalize_specsys_name(_first_valid_value(m, table, ("SPECSYS", "SSYSOBS")), default="TOPOCENT")
    # For this calibrator, frequency-axis inputs are normalized such that the
    # observation-frame bookkeeping matches the current axis frame.
    ssysobs = specsys

    m["CTYPE1"] = ctype1
    m["CUNIT1"] = "Hz"
    m["SPECSYS"] = specsys
    m["SSYSOBS"] = ssysobs
    m["VELDEF"] = "RADIO"
    m["RESTFRQ"] = rest_hz
    m["RESTFREQ"] = rest_hz
    return m


def _coerce_velocity_array_to_ms(values: np.ndarray, *, source_name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"Velocity column {source_name!r} contains non-finite values.")

    name = str(source_name).strip().upper()

    # Standardized semantics in this pipeline:
    #   VELOSYS / VFRAME  -> m/s
    #   V_CORR_KMS/*_KMS  -> km/s
    # Avoid magnitude-based guessing because genuine small m/s corrections can be
    # mistaken for km/s near a velocity-crossing.
    if name in ("VELOSYS", "VFRAME"):
        return arr
    if name == "V_CORR_KMS" or name.endswith("_KMS"):
        return arr * 1000.0

    return arr


def _extract_existing_velocity_ms(
    mapping_on: pd.DataFrame,
    *,
    preferred_key: str = "VFRAME",
) -> tuple[np.ndarray | None, str | None]:
    preferred = str(preferred_key or "").strip().upper()
    candidate_keys: list[str] = []
    for key in (preferred, "VFRAME", "V_CORR_KMS", "VELOSYS"):
        if key and key not in candidate_keys:
            candidate_keys.append(key)

    found: list[tuple[str, np.ndarray]] = []
    for key in candidate_keys:
        if key not in mapping_on.columns:
            continue
        col = mapping_on[key]
        # If duplicated column names exist, prefer not to guess silently.
        if isinstance(col, pd.DataFrame):
            raise ValueError(f"Duplicate velocity columns found for {key!r}; please clean the input table.")
        arr = pd.to_numeric(col, errors="coerce").to_numpy(dtype=float)
        if np.all(np.isfinite(arr)):
            found.append((key, _coerce_velocity_array_to_ms(arr, source_name=key)))

    if not found:
        return None, None

    ref_key, ref_arr = found[0]
    for key, arr in found[1:]:
        if not np.allclose(arr, ref_arr, rtol=0.0, atol=1e-6):
            raise ValueError(
                "Inconsistent per-row velocity bookkeeping: "
                f"{ref_key} and {key} disagree. "
                "Expected VELOSYS/VFRAME/V_CORR_KMS to represent the same correction."
            )

    return ref_arr, ref_key



def _fast_to_datetime(series: pd.Series) -> pd.DatetimeIndex:
    if series.empty:
        return pd.DatetimeIndex([])
    dtype = series.dtype
    if pd.api.types.is_datetime64_any_dtype(dtype):
        return pd.DatetimeIndex(pd.to_datetime(series, utc=True, errors="coerce"))
    vals = series.values
    if pd.api.types.is_numeric_dtype(dtype):
        arr = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
        finite = arr[np.isfinite(arr)]
        if finite.size > 0 and float(np.nanmedian(finite)) < 1.0e6:
            # Likely MJD days (SDFITS TIME)
            return pd.DatetimeIndex(pd.to_datetime(arr, unit="D", origin=pd.Timestamp("1858-11-17"), utc=True, errors="coerce"))
        return pd.DatetimeIndex(pd.to_datetime(arr, unit="s", utc=True, errors="coerce"))
    return pd.DatetimeIndex(pd.to_datetime(vals, utc=True, errors="coerce"))


def _datetime_index_to_unix_seconds(index: pd.DatetimeIndex) -> np.ndarray:
    idx = pd.DatetimeIndex(pd.to_datetime(index, utc=True, errors="coerce"))
    if len(idx) == 0:
        return np.empty(0, dtype=float)
    epoch = pd.Timestamp("1970-01-01", tz="UTC")
    secs = (idx - epoch) / pd.Timedelta(seconds=1)
    return np.asarray(secs, dtype=float)


def _mean_datetime_index(index: pd.DatetimeIndex) -> pd.Timestamp:
    secs = _datetime_index_to_unix_seconds(index)
    finite = secs[np.isfinite(secs)]
    if finite.size == 0:
        return pd.Timestamp("NaT", tz="UTC")
    return pd.Timestamp(pd.to_datetime(float(np.mean(finite)), unit="s", utc=True))


def _interp_spectra_in_time(
    target_times: pd.DatetimeIndex,
    ref_times: pd.DatetimeIndex,
    ref_data: np.ndarray
) -> np.ndarray:
    t_tgt = _datetime_index_to_unix_seconds(target_times)
    t_ref = _datetime_index_to_unix_seconds(ref_times)
    
    if np.any(np.diff(t_ref) < 0):
        sort_idx = np.argsort(t_ref)
        t_ref = t_ref[sort_idx]
        ref_data = ref_data[sort_idx]
        
    n_tgt = len(t_tgt)
    n_chan = ref_data.shape[1]
    out_data = np.zeros((n_tgt, n_chan), dtype=ref_data.dtype)
    
    mask_before = t_tgt <= t_ref[0]
    if np.any(mask_before):
        out_data[mask_before] = ref_data[0]
        
    mask_after = t_tgt >= t_ref[-1]
    if np.any(mask_after):
        out_data[mask_after] = ref_data[-1]
        
    mask_mid = (~mask_before) & (~mask_after)
    if np.any(mask_mid):
        t_mid = t_tgt[mask_mid]
        idxs = np.searchsorted(t_ref, t_mid, side='left')
        
        t0 = t_ref[idxs - 1]
        t1 = t_ref[idxs]
        
        denom = t1 - t0
        denom[denom == 0] = 1.0 
        w = (t_mid - t0) / denom
        w = w[:, None].astype(ref_data.dtype)
        
        d0 = ref_data[idxs - 1]
        d1 = ref_data[idxs]
        
        out_data[mask_mid] = (1.0 - w) * d0 + w * d1
        
    return out_data

def _average_mode_by_scan(
    times: pd.DatetimeIndex, 
    data: np.ndarray, 
    mapping_subset: pd.DataFrame
) -> tuple[pd.DatetimeIndex, np.ndarray]:
    """スキャンごとにデータを時間平均し、S/N比を向上させるためのヘルパー関数"""
    scan_col = None
    for c in ["SCAN", "SCANID", "scan", "scanid", "scan_id"]:
        if c in mapping_subset.columns:
            scan_col = c
            break
            
    if not scan_col or len(times) == 0:
        return times, data
        
    t_out = []
    d_out = []
    for scan_id, group_idxs in mapping_subset.groupby(scan_col).groups.items():
        idx_array = np.array(list(group_idxs), dtype=int)
        if len(idx_array) == 0: continue
        
        # 時刻の平均値算出（datetime内部単位に依存しない秒基準）
        mean_time = _mean_datetime_index(times[idx_array])
        
        # データの積分（平均）
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_data = np.nanmean(data[idx_array], axis=0)
        
        t_out.append(mean_time)
        d_out.append(mean_data)
        
    if not t_out:
        return times, data
        
    # 時系列順を保証して返す
    sort_idx = np.argsort(_datetime_index_to_unix_seconds(pd.DatetimeIndex(t_out)))
    t_out_sorted = pd.DatetimeIndex([t_out[i] for i in sort_idx])
    d_out_sorted = np.stack([d_out[i] for i in sort_idx])
    
    return t_out_sorted, d_out_sorted


def make_tastar_dumps(
    raw,
    *,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    
    # --- キャリブレーションおよび大気モデル(ATM)用パラメータ ---
    t_hot_k: Optional[float] = None,              # 常温負荷（ホットロード）の物理温度 [K]
    tau_zenith: Union[float, str, None] = None,   # 天頂光学的厚み。数値か "auto"。未指定時は Basic 1-Temp モデル
    t_surface_k: Optional[float] = None,          # 屋外の地表気温 [K]。未指定時はヘッダ(WXTEMP等)から自動取得
    t_atm_model: str = "offset",                  # 空の実効温度 T_atm の推定モデル ("offset" or "ratio")
    t_atm_delta_k: float = 15.0,                  # offsetモデル使用時の差し引き温度差 [K]
    t_atm_eta: float = 0.95,                      # ratioモデル使用時の係数
    gain_mode: str = "hybrid",                    # "hybrid" (推奨) または "independent"
    verbose: bool = True,                         # ログ出力フラグ
    # -----------------------------------------------------
    
    ch_range: Optional[Tuple[int, int]] = None,
    vlsrk_range_kms: Optional[Tuple[float, float]] = None,
    v_corr_col: str = "VFRAME",
    coord_frame: str | None = None,
    vcorr_chunk_sec: Optional[float] = None,
    dtype: Optional[Union[type, str]] = None, 
    rest_freq: Optional[float] = None,
) -> Scantable:
    """
    Calibrate ON dumps to Ta* using chopper wheel calibration with optional Atmospheric (ATM) Model.
    """
    meta = dict(raw["meta"])
    on_arr = raw["on"]
    off_arr = raw["off"]
    hot_arr = raw["hot"]
    mapping_all: pd.DataFrame = _df_to_native_endian(raw.get("mapping", pd.DataFrame()).copy())

    # gain_mode の正規化と厳格チェック（typoで黙って hybrid に落ちる事故を防ぐ）
    gain_mode_norm = str(gain_mode).strip().lower()
    if gain_mode_norm not in ("hybrid", "independent"):
        raise ValueError(
            f"Invalid gain_mode: {gain_mode!r}. Must be 'hybrid' or 'independent'."
        )


    v_corr_col_norm = str(v_corr_col or "VFRAME").strip().upper()

    # 静止周波数の上書き処理
    if rest_freq is not None:
        apply_restfreq_override(meta, mapping_all, float(rest_freq), require_wcs_for_vrad=True)

    meta = _canonicalize_frequency_axis_meta(meta, mapping_all)

    # タイムスタンプの取得
    ts_all = _resolve_table_timestamps(mapping_all)
    if ts_all is None:
        raise ValueError("mapping table must contain TIMESTAMP, MJD, DATE-OBS/DATEOBS, or legacy TIME")

    if ts_all.isna().any():
        raise ValueError("mapping timestamp contains NaN; cannot calibrate reliably")

    # 行のフィルタリング処理
    if rows is not None and exclude_rows is not None:
        raise ValueError("Cannot specify both 'rows' and 'exclude_rows'.")

    total_len = len(mapping_all)
    final_idxs = None
    
    if rows is not None:
        final_idxs = _parse_row_selector(rows, total_len)
    elif exclude_rows is not None:
        exclude_idxs = _parse_row_selector(exclude_rows, total_len)
        final_idxs = np.setdiff1d(np.arange(total_len), exclude_idxs)
        
    if final_idxs is not None:
        final_idxs = np.unique(final_idxs)
        if len(final_idxs) == 0:
            raise ValueError("Selector resulted in empty input table.")

        keep_mask_global = np.zeros(total_len, dtype=bool)
        keep_mask_global[final_idxs] = True
        orig_obsmode = mapping_all["OBSMODE"].astype(str).str.strip().str.upper()

        on_arr = on_arr[keep_mask_global[orig_obsmode == 'ON']]
        off_arr = off_arr[keep_mask_global[orig_obsmode == 'OFF']]
        hot_arr = hot_arr[keep_mask_global[orig_obsmode == 'HOT']]

        mapping_all = mapping_all.iloc[final_idxs].reset_index(drop=True)
        ts_all = ts_all[final_idxs]
    
    # 観測モードごとに分割
    obsmode = mapping_all["OBSMODE"].astype(str).str.strip().str.upper()

    def _take_mode_times(mode: str, n: int) -> pd.DatetimeIndex:
        t = ts_all[obsmode == mode][:n]
        if len(t) != n:
            raise ValueError(f"OBSMODE={mode} rows ({len(t)}) do not match spectra count ({n})")
        return t

    t_hot = _take_mode_times("HOT", int(hot_arr.shape[0]))
    t_off = _take_mode_times("OFF", int(off_arr.shape[0]))
    t_on = _take_mode_times("ON", int(on_arr.shape[0]))

    if len(t_on) == 0: raise ValueError("No ON scans remaining.")
    if len(t_hot) == 0: raise ValueError("No HOT scans remaining.")
    if len(t_off) == 0: raise ValueError("No OFF scans remaining.")

    mapping_on = mapping_all.loc[obsmode == "ON"].iloc[: len(t_on)].copy()
    mapping_on.index = t_on

    # Canonicalize coordinate columns when available.
    # TOPOCENT inputs require coordinates to compute VELOSYS unless a valid row-wise
    # VELOSYS/VFRAME column is already present. LSRK inputs do not require coordinates.
    effective_coord_frame = _choose_effective_coord_frame(mapping_on, coord_frame, meta)
    lon_col0, lat_col0 = _resolve_lonlat_columns(mapping_on, effective_coord_frame)
    have_coords = (lon_col0 is not None and lat_col0 is not None)
    need_coords = str(meta.get("SPECSYS", "TOPOCENT")).strip().upper() == "TOPOCENT"

    if have_coords:
        if effective_coord_frame == "galactic":
            if "GLON" not in mapping_on.columns:
                mapping_on["GLON"] = pd.to_numeric(mapping_on[lon_col0], errors="coerce")
            if "GLAT" not in mapping_on.columns:
                mapping_on["GLAT"] = pd.to_numeric(mapping_on[lat_col0], errors="coerce")
        else:
            if "RA" not in mapping_on.columns:
                mapping_on["RA"] = pd.to_numeric(mapping_on[lon_col0], errors="coerce")
            if "DEC" not in mapping_on.columns:
                mapping_on["DEC"] = pd.to_numeric(mapping_on[lat_col0], errors="coerce")

        _maybe_validate_mapping_frame(meta, mapping_on)
    elif need_coords and _extract_existing_velocity_ms(mapping_on, preferred_key=v_corr_col_norm)[0] is None:
        if effective_coord_frame == "galactic":
            raise ValueError("coord_frame=galactic but mapping lacks GLON/GLAT (or aliases).")
        raise ValueError("coord_frame=icrs but mapping lacks RA/DEC (or aliases).")

    # -------------------------------------------------------------------------
    # Spectral-frame handling and per-row velocity bookkeeping
    # -------------------------------------------------------------------------
    specsys_in = str(meta["SPECSYS"]).strip().upper()
    ssysobs_in = specsys_in

    velosys_ms = None
    v_corr_kms = None

    if specsys_in == "TOPOCENT":
        # Prefer an existing per-row standard/legacy velocity column when present.
        # This avoids silently changing already-computed corrections and reduces
        # dependence on coordinate-frame inference for backward-compatible inputs.
        velosys_ms, vel_source = _extract_existing_velocity_ms(mapping_on, preferred_key=v_corr_col_norm)

        if velosys_ms is None:
            lon_col, lat_col = _resolve_lonlat_columns(mapping_on, effective_coord_frame)

            if lon_col and lat_col:
                ra_arr = pd.to_numeric(mapping_on[lon_col], errors="coerce").to_numpy(dtype=float)
                dec_arr = pd.to_numeric(mapping_on[lat_col], errors="coerce").to_numpy(dtype=float)

                if np.all(np.isfinite(ra_arr)) and np.all(np.isfinite(dec_arr)):
                    should_decimate = (vcorr_chunk_sec is not None and vcorr_chunk_sec > 0 and len(t_on) > 2)

                    if should_decimate:
                        t_vals = _datetime_index_to_unix_seconds(t_on)
                        t_start, t_end = t_vals[0], t_vals[-1]
                        duration = t_end - t_start

                        n_steps = builtins.max(1, int(np.ceil(duration / vcorr_chunk_sec)))
                        grid_times = np.linspace(t_start, t_end, num=n_steps + 1)
                        idxs = np.clip(np.unique(np.searchsorted(t_vals, grid_times)), 0, len(t_vals) - 1)

                        v_corr_sparse = compute_vcorr_series(
                            times=t_on[idxs], ra_deg=ra_arr[idxs], dec_deg=dec_arr[idxs],
                            meta=meta, coord_frame=effective_coord_frame,
                        )
                        v_corr_kms = np.interp(t_vals, t_vals[idxs], v_corr_sparse)
                    else:
                        v_corr_kms = compute_vcorr_series(
                            times=t_on, ra_deg=ra_arr, dec_deg=dec_arr,
                            meta=meta, coord_frame=effective_coord_frame,
                        )
                    velosys_ms = np.asarray(v_corr_kms, dtype=float) * 1000.0

        if velosys_ms is None:
            raise ValueError(
                "TOPOCENT frequency-axis input requires per-row VELOSYS/VFRAME, "
                "or enough coordinate/time/site metadata to compute it."
            )

        v_corr_kms = np.asarray(velosys_ms, dtype=float) / 1000.0
        mapping_on["VELOSYS"] = np.asarray(velosys_ms, dtype=float)
        mapping_on["VFRAME"] = np.asarray(velosys_ms, dtype=float)
    else:
        v_corr_kms = np.zeros(len(t_on), dtype=float)

    # -------------------------------------------------------------------------
    # Channel Slice
    # -------------------------------------------------------------------------
    ch_start, ch_stop = 0, int(meta["NAXIS1"])

    if vlsrk_range_kms is not None:
        a, b = vlsrk_range_kms
        s, e = channel_slice_from_vrange_union(meta, v_corr_kms, a, b)
        ch_start, ch_stop = s, e

    if ch_range is not None:
        s, e = int(ch_range[0]), int(ch_range[1])
        if e <= s: raise ValueError("invalid ch_range")
        ch_start, ch_stop = builtins.max(ch_start, s), builtins.min(ch_stop, e)

    if ch_stop <= ch_start: raise ValueError("channel slice became empty")

    meta2 = wcs_slice_channels(meta, ch_start, ch_stop)
    rest = None
    for k in ("RESTFRQ", "RESTFREQ"):
        if k in meta and meta[k] not in (None, ""):
            rest = float(meta[k])
            break
    if rest is not None:
        meta2["RESTFRQ"] = rest
        meta2["RESTFREQ"] = rest

    for k in ("SITELAT", "SITELONG", "SITEELEV", "OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z",
              "site_lat_deg", "site_lon_deg", "site_height_m", "obsgeo_x_m", "obsgeo_y_m", "obsgeo_z_m", "site_name"):
        if k in meta and meta[k] not in (None, ""):
            meta2[k] = meta[k]

    meta2["CTYPE1"] = "FREQ"
    meta2["CUNIT1"] = "Hz"
    meta2["SPECSYS"] = specsys_in
    meta2["SSYSOBS"] = ssysobs_in
    meta2["VELDEF"] = "RADIO"
    meta2["RESTFRQ"] = float(meta["RESTFRQ"])
    meta2["RESTFREQ"] = float(meta["RESTFREQ"])
    meta2.pop("VELOSYS", None)
    meta2.pop("VFRAME", None)

    # -------------------------------------------------------------------------
    # [NEW] Calculate Equivalent Calibration Temperature (T_cal)
    # -------------------------------------------------------------------------
    # 1. t_hot_k の自動取得 (引数優先 -> mapping_all -> meta -> エラー)
    # NOTE: t_hot は HOT の時刻配列として既に使用しているため、温度は別名にする
    t_hot_load = t_hot_k
    if t_hot_load is None:
        t_hot_load = extract_meta_value(mapping_all, meta, ("THOT", "TLOAD", "TAMB_K"))
        if t_hot_load is None:
            # [FIX] 勝手に 300.0 K にせず、厳格にエラーを返す
            raise ValueError(
                "Cannot determine hot load temperature. "
                "'t_hot_k' must be provided as an argument, or 'THOT'/'TLOAD' must exist in metadata."
            )
    else:
        t_hot_load = float(t_hot_load)

    # 【追加】T_HOT が摂氏 (例えば 20.0) で入力・記録されていた場合の保護
    if t_hot_load < 100.0:
        import warnings
        warnings.warn(f"T_HOT (or t_hot_k) is {t_hot_load}. Assuming Celsius and converting to Kelvin.")
        t_hot_load += 273.15
    calstat_str = "TASTAR_BASIC"
    used_tau = None
    
    if tau_zenith is not None:
        # --- ATM 補正モデル発動 ---
        
        # 1. tau の特定 (配列対応 ＆ mapping_allフォールバック)
        if isinstance(tau_zenith, str) and tau_zenith.strip().lower() in ("auto", "header"):
            # まず ON点 の時系列配列として探す
            used_tau = extract_meta_array(mapping_on, meta, ("TAU0", "TAU", "OPACITY"))
            # 見つからない、または全てNaNなら、全体(HOT等)のスカラー値を探す
            if used_tau is None or (isinstance(used_tau, np.ndarray) and np.all(np.isnan(used_tau))):
                used_tau = extract_meta_value(mapping_all, meta, ("TAU0", "TAU", "OPACITY"))
                
            if used_tau is None or (isinstance(used_tau, np.ndarray) and np.all(np.isnan(used_tau))):
                raise ValueError("tau_zenith='auto' requested, but 'TAU0' (or alias) not found in mapping or meta.")
        else:
            try:
                used_tau = float(tau_zenith)
            except (ValueError, TypeError):  # [FIX] TypeErrorも安全にキャッチ
                raise ValueError(f"Invalid tau_zenith value: {tau_zenith}. Must be float or 'auto'.")
                
        # 2. 外気温 (T_surface) の特定 (配列対応 ＆ mapping_allフォールバック)
        t_surf = t_surface_k
        if t_surf is None:
            t_surf = extract_meta_array(mapping_on, meta, ("WXTEMP", "TAMBIENT", "TAMB", "T_SURF"))
            if t_surf is None or (isinstance(t_surf, np.ndarray) and np.all(np.isnan(t_surf))):
                t_surf = extract_meta_value(mapping_all, meta, ("WXTEMP", "TAMBIENT", "TAMB", "T_SURF"))
                
            if t_surf is None or (isinstance(t_surf, np.ndarray) and np.all(np.isnan(t_surf))):
                 raise ValueError("t_surface_k must be provided or 'WXTEMP' / 'TAMBIENT' must exist in mapping for ATM model.")

        # 【安全装置】摂氏・ケルビン変換 (配列・スカラー両対応)
        if isinstance(t_surf, np.ndarray):
            t_surf = np.where(t_surf < 100.0, t_surf + 273.15, t_surf)
        else:
            if t_surf < 100.0:
                t_surf += 273.15
        
        # 3. 仰角 (Elevation) の取得
        el_col = next((c for c in ("ELEVATIO", "EL", "ELEVATION") if c in mapping_on.columns), None)
        if el_col is None:
            raise ValueError("Elevation column ('ELEVATIO', 'EL', etc.) is required for ATM model.")
            
        el_deg_arr = pd.to_numeric(mapping_on[el_col], errors="coerce").to_numpy(dtype=float)
        if np.any(np.isnan(el_deg_arr)):
             raise ValueError("Elevation data contains NaN.")

        # 4. T_atm の推定と T_cal 配列の計算
        t_atm = estimate_t_atm(
            t_surface_k=t_surf, 
            model=t_atm_model, 
            delta_t=t_atm_delta_k, 
            eta=t_atm_eta
        )
        t_cal = compute_t_cal_array(t_hot_load, t_atm, used_tau, el_deg_arr)
        calstat_str = "TASTAR_ATM"
        
        if verbose:
            mean_tsurf = np.mean(t_surf) if isinstance(t_surf, np.ndarray) else t_surf
            mean_tatm = np.mean(t_atm) if isinstance(t_atm, np.ndarray) else t_atm
            mean_tau = np.mean(used_tau) if isinstance(used_tau, np.ndarray) else used_tau
            print(f"[CALIBRATE] Applied ATM Model: tau={mean_tau:.3f}, T_surf={mean_tsurf:.1f}K, "
                  f"T_atm={mean_tatm:.1f}K (model='{t_atm_model}'), T_hot={t_hot_load}K.")
            print(f"            Mean T_cal = {np.mean(t_cal):.2f} K")
            
    else:
        # --- Basic 1-Temp モデル (デフォルト) ---
        if verbose:
            print(f"[CALIBRATE] Applied Basic 1-Temp Model (T_cal = {t_hot_load} K).")
        t_cal = np.full(len(mapping_on), float(t_hot_load))
        
    # -------------------------------------------------------------------------
    # Load and Calibrate Data
    # -------------------------------------------------------------------------
    store_dt = dtype if dtype is not None else np.float64
    calc_dt = np.float64 

    on2_arr = np.asarray(on_arr, dtype=store_dt)[:, ch_start:ch_stop]
    off2_arr = np.asarray(off_arr, dtype=store_dt)[:, ch_start:ch_stop]
    hot2_arr = np.asarray(hot_arr, dtype=store_dt)[:, ch_start:ch_stop]

    # [NEW] HOTとOFFをスキャンごとに時間平均(積分)してS/Nを向上させる
    mapping_hot = mapping_all.loc[obsmode == "HOT"].iloc[: len(t_hot)].reset_index(drop=True)
    mapping_off = mapping_all.loc[obsmode == "OFF"].iloc[: len(t_off)].reset_index(drop=True)
    
    t_hot_avg, hot2_avg = _average_mode_by_scan(t_hot, hot2_arr, mapping_hot)
    t_off_avg, off2_avg = _average_mode_by_scan(t_off, off2_arr, mapping_off)

    # 分子（空の引き算）用: ON点時刻におけるOFFスペクトルの内挿
    off_on_interp = _interp_spectra_in_time(t_on, t_off_avg, off2_avg)

    if gain_mode_norm == "independent":
        if verbose: print("[CALIBRATE] Using 'independent' gain mode.")
        # 従来方式: HOTとOFFをON時刻で別々に内挿してから引き算
        hot_on_interp = _interp_spectra_in_time(t_on, t_hot_avg, hot2_avg)
        denom = hot_on_interp.astype(calc_dt) - off_on_interp.astype(calc_dt)
    else:
        if verbose: print("[CALIBRATE] Using 'hybrid' gain mode.")
        # ハイブリッド方式: HOT取得時刻におけるOFFを内挿し、定在波を相殺した純粋な分母(ゲインアレイ)を作成
        off_at_hot = _interp_spectra_in_time(t_hot_avg, t_off_avg, off2_avg)
        denom_at_hot = hot2_avg.astype(calc_dt) - off_at_hot.astype(calc_dt)
        # そのクリーンな分母アレイをON点の時刻に内挿して適用
        denom = _interp_spectra_in_time(t_on, t_hot_avg, denom_at_hot)

    # ゲインと Ta* の計算 (T_cal 配列を 2D にブロードキャストして適用)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_cal_2d = t_cal[:, None].astype(calc_dt)
        gain = t_cal_2d / denom
        ta_calc = (on2_arr.astype(calc_dt) - off_on_interp.astype(calc_dt)) * gain
        
    ta_calc[~np.isfinite(ta_calc)] = np.nan
    ta = ta_calc.astype(store_dt)

    # -------------------------------------------------------------------------
    # Metadata Table Construction
    # -------------------------------------------------------------------------
    table = mapping_on.copy()
    if not table.index.equals(t_on):
        table = table.reindex(t_on)

    # メタデータへ履歴を記録 (Undo / Recalibrate用)
    table["THOT"] = float(t_hot_load)
    table["TAMB_K"] = float(t_hot_load) # 既存・外部ツールとの後方互換性のため維持
    table["TCAL"] = t_cal          # 適用された T_cal 配列 (各行スカラー)
    table["CALSTAT"] = calstat_str
    
    if used_tau is not None:
        # [FIX] 配列ならそのまま、スカラーならfloat変換
        if isinstance(used_tau, np.ndarray):
            table["TAU0"] = used_tau
        else:
            table["TAU0"] = float(used_tau)
    table["CH_START"] = int(ch_start)
    table["CH_STOP"] = int(ch_stop)
    table["NCHAN_SEL"] = int(ch_stop - ch_start)

    meta2["TEMPSCAL"] = "TA*"
    table["TEMPSCAL"] = "TA*"

    if rest_freq is not None:
        apply_restfreq_override(meta2, table, float(rest_freq), require_wcs_for_vrad=True)

    # Standardized spectral-axis metadata (public columns only)
    rest_hz = float(meta2["RESTFRQ"])

    # Drop internal-only columns from prior processing chains
    for col in list(table.columns):
        if str(col).upper().startswith("SIG_"):
            del table[col]
    for col in ("V_CORR_KMS", "v_corr_kms"):
        if col in table.columns:
            del table[col]

    table["CRVAL1"] = float(meta2.get("CRVAL1", table["CRVAL1"].iloc[0] if "CRVAL1" in table.columns else np.nan))
    table["CDELT1"] = float(meta2.get("CDELT1", table["CDELT1"].iloc[0] if "CDELT1" in table.columns else np.nan))
    table["CRPIX1"] = float(meta2.get("CRPIX1", table["CRPIX1"].iloc[0] if "CRPIX1" in table.columns else 1.0))
    table["CTYPE1"] = "FREQ"
    table["CUNIT1"] = "Hz"
    table["RESTFREQ"] = rest_hz
    table["RESTFRQ"] = rest_hz
    table["SPECSYS"] = specsys_in
    table["SSYSOBS"] = ssysobs_in
    table["VELDEF"] = "RADIO"

    if specsys_in == "TOPOCENT":
        table["VELOSYS"] = np.asarray(velosys_ms, dtype=float)
        table["VFRAME"] = np.asarray(velosys_ms, dtype=float)
    else:
        for col in ("VELOSYS", "VFRAME"):
            if col in table.columns:
                del table[col]

    return Scantable(meta=meta2, data=ta, table=table)


def tastar_from_rawspec(
    raw: dict[str, Any],
    t_hot_k: Optional[float] = None, 
    *,
    ch_range: tuple[int, int] | None = None,
    vlsrk_range_kms: tuple[float, float] | None = None,
    vcorr_chunk_sec: float | None = None,
    dtype: type | str | None = None,
    rest_freq: float | None = None,
    v_corr_col: str = "VFRAME",
    coord_frame: str | None = None,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    
    # 大気モデル(ATM)用パラメータのパススルー
    tau_zenith: Union[float, str, None] = None,
    t_surface_k: Optional[float] = None,
    t_atm_model: str = "offset",
    t_atm_delta_k: float = 15.0,
    t_atm_eta: float = 0.95,
    gain_mode: str = "hybrid",
    verbose: bool = True,
) -> Scantable:
    """Public helper used by CLI scripts."""
    return make_tastar_dumps(
        raw,
        t_hot_k=t_hot_k,
        ch_range=ch_range,
        vlsrk_range_kms=vlsrk_range_kms,
        vcorr_chunk_sec=vcorr_chunk_sec,
        dtype=dtype,
        rest_freq=rest_freq,
        v_corr_col=v_corr_col,
        coord_frame=(None if coord_frame in (None, "", "auto") else str(coord_frame)),
        rows=rows,
        exclude_rows=exclude_rows,
        tau_zenith=tau_zenith,
        t_surface_k=t_surface_k,
        t_atm_model=t_atm_model,
        t_atm_delta_k=t_atm_delta_k,
        t_atm_eta=t_atm_eta,
        gain_mode=gain_mode,
        verbose=verbose,
    )


def run_tastar_calibration(
    input_data: Union[str, Scantable, dict],
    output_path: Optional[str] = None,
    t_hot_k: Optional[float] = None,
    ch_range: Optional[Tuple[int, int]] = None,
    vlsrk_range_kms: Optional[Tuple[float, float]] = None,
    coord_frame: str | None = None,
    spectrum_column: str = "DATA",
    overwrite: bool = False,
    store_freq_column: bool = False,
    vcorr_chunk_sec: Optional[float] = None,
    dtype: Optional[Union[type, str]] = None,
    rest_freq: Optional[float] = None,
    v_corr_col: str = "VFRAME",
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    
    # 大気モデル(ATM)用パラメータのパススルー
    tau_zenith: Union[float, str, None] = None,
    t_surface_k: Optional[float] = None,
    t_atm_model: str = "offset",
    t_atm_delta_k: float = 15.0,
    t_atm_eta: float = 0.95,
    gain_mode: str = "hybrid",
    verbose: bool = True,
) -> Scantable:
    """
    High-level function: Load Raw data, calibrate to Ta*, and optionally save to SDFITS.
    """
    if isinstance(input_data, str):
        raw = load_rawspec_auto(input_data)
        input_name = input_data
    elif isinstance(input_data, Scantable):
        input_name = "scantable_object"
        table = input_data.table
        data = input_data.data
        if "OBSMODE" not in table.columns: raise ValueError("Must contain 'OBSMODE'")
        
        def _get_data_by_mode(mode_str):
            return data[table["OBSMODE"].astype(str).str.strip().str.upper() == mode_str]

        raw = {
            "meta": input_data.meta, "mapping": table,
            "on": _get_data_by_mode("ON"), "off": _get_data_by_mode("OFF"), "hot": _get_data_by_mode("HOT")
        }
    else:
        raw = input_data
        input_name = "memory_object"

    res = make_tastar_dumps(
        raw,
        t_hot_k=t_hot_k,
        ch_range=ch_range,
        vlsrk_range_kms=vlsrk_range_kms,
        coord_frame=(None if coord_frame in (None, "", "auto") else str(coord_frame)),
        vcorr_chunk_sec=vcorr_chunk_sec, 
        dtype=dtype,
        rest_freq=rest_freq, 
        v_corr_col=v_corr_col,
        rows=rows,
        exclude_rows=exclude_rows,
        tau_zenith=tau_zenith,
        t_surface_k=t_surface_k,
        t_atm_model=t_atm_model,
        t_atm_delta_k=t_atm_delta_k,
        t_atm_eta=t_atm_eta,
        gain_mode=gain_mode,
        verbose=verbose,
    )

    # 3. Create History
    # [FIX] 実際に適用された THOT を結果テーブルから取得して記録する
    actual_thot = float(res.table["THOT"].iloc[0])
    
    history = {
        'task': 'sdrsf-make-tastar', 'input': input_name,
        't_hot_k': actual_thot, 't_hot_k_input': str(t_hot_k), 'tau_zenith': str(tau_zenith),
        't_surface_k': str(t_surface_k), 't_atm_model': t_atm_model,
        't_atm_delta_k': float(t_atm_delta_k), 't_atm_eta': float(t_atm_eta),
        'gain_mode': gain_mode,
        'notes': 'Ta* calibration.',
        'ch_range': str(ch_range) if ch_range else "all",
        'vlsrk_range': str(vlsrk_range_kms) if vlsrk_range_kms else "none",
        'vcorr_chunk_sec': str(vcorr_chunk_sec), 'dtype': str(dtype),
        'rest_freq': str(rest_freq), 'v_corr_col': str(v_corr_col), 'rows': str(rows), 'exclude_rows': str(exclude_rows),
    }
    res.history = history

    if output_path:
        write_scantable(
            path=output_path, scantable=res,
            spectrum_column=str(spectrum_column).upper(), overwrite=overwrite, store_freq_column=store_freq_column 
        )
        if verbose: print(f"Saved: {output_path}")

    return res


def recalibrate_tastar(
    scantable: Scantable,
    new_tau: Union[float, None] = None,
    new_t_surface_k: Optional[float] = None,
    new_t_atm_model: str = "offset",
    new_t_atm_delta_k: float = 15.0,
    new_t_atm_eta: float = 0.95,
    verbose: bool = True
) -> Scantable:
    """
    すでにキャリブレーションされた Scantable の大気補正パラメータを、
    RAWデータを読み直すことなく瞬時に再適用（または解除）するユーティリティ。
    
    Args:
        scantable (Scantable): 対象の Scantable オブジェクト。
        new_tau (float | None): 新しい tau。None の場合は Basic 1-Temp モデルに Undo する。
        new_t_surface_k (float): 新しい外気温 [K]。未指定の場合はヘッダから探す。
        new_t_atm_model (str): T_atm の推定モデル。
        new_t_atm_delta_k (float): offsetモデルのパラメータ。
        new_t_atm_eta (float): ratioモデルのパラメータ。
        verbose (bool): ログ出力フラグ。
        
    Returns:
        Scantable: 再計算された新しいオブジェクト。
    """
    table = scantable.table.copy()
    data = scantable.data.copy()
    meta = scantable.meta.copy()
    
    # 過去データで THOT が無い場合は TAMB_K にフォールバックする安全策
    if "THOT" in table.columns:
        t_hot_k = float(table["THOT"].iloc[0])
    elif "TAMB_K" in table.columns:
        t_hot_k = float(table["TAMB_K"].iloc[0])
        if verbose: print("[WARNING] Using legacy 'TAMB_K' column as T_HOT.")
    else:
        raise ValueError("Cannot recalibrate: 'THOT' missing from mapping table.")

    # 【追加】レガシーデータで T_HOT が摂氏で記録されていた場合の保護
    if t_hot_k < 100.0:
        t_hot_k += 273.15

    if "TCAL" not in table.columns:
        raise ValueError("Cannot recalibrate: 'TCAL' missing from mapping table.")
        
    old_t_cal = table["TCAL"].to_numpy(dtype=float)
    
    if new_tau is None:
        # --- Undo (Basic モデルに戻す) ---
        if verbose: print(f"[RECALIBRATE] Reverting to Basic 1-Temp Model (T_cal = T_hot = {t_hot_k} K).")
        new_t_cal = np.full(len(table), t_hot_k)
        table["CALSTAT"] = "TASTAR_BASIC"
        if "TAU0" in table.columns: del table["TAU0"]
        
    else:
        # --- 新しいパラメータで ATM モデルを適用 ---
        t_surf = new_t_surface_k
        if t_surf is None:
            # [FIX] extract_meta_value ではなく extract_meta_array を使用し、時変データを引き継ぐ
            t_surf = extract_meta_array(table, meta, ("WXTEMP", "TAMBIENT", "TAMB", "T_SURF"))
            
            # [FIX] 全てが NaN の場合のエラーハンドリングを追加
            if t_surf is None or (isinstance(t_surf, np.ndarray) and np.all(np.isnan(t_surf))):
                 raise ValueError("new_t_surface_k must be provided or 'WXTEMP' / 'TAMBIENT' must exist in mapping.")
            
        # ケルビン変換の安全装置
        if isinstance(t_surf, np.ndarray):
            t_surf = np.where(t_surf < 100.0, t_surf + 273.15, t_surf)
        else:
            if t_surf < 100.0:
                t_surf += 273.15
                    

        el_col = next((c for c in ("ELEVATIO", "EL", "ELEVATION") if c in table.columns), None)
        if el_col is None: raise ValueError("Elevation column missing. Cannot apply ATM model.")
        
        el_deg = pd.to_numeric(table[el_col], errors="coerce").to_numpy(dtype=float)
        
        t_atm = estimate_t_atm(
            t_surface_k=t_surf, 
            model=new_t_atm_model, 
            delta_t=new_t_atm_delta_k, 
            eta=new_t_atm_eta
        )
        new_t_cal = compute_t_cal_array(t_hot_k, t_atm, new_tau, el_deg)
        
        # [修正後] 配列・スカラー両対応にする
        if isinstance(new_tau, np.ndarray):
            table["TAU0"] = new_tau
        else:
            table["TAU0"] = float(new_tau)
        
        table["CALSTAT"] = "TASTAR_ATM"
        
        if verbose:
            # ログ出力用：配列なら平均値、スカラーならそのまま
            mean_tau = np.mean(new_tau) if isinstance(new_tau, np.ndarray) else new_tau
            mean_tsurf = np.mean(t_surf) if isinstance(t_surf, np.ndarray) else t_surf
            mean_tatm = np.mean(t_atm) if isinstance(t_atm, np.ndarray) else t_atm
        
            # フォーマットには平均値(mean_tau)を使う
            print(f"[RECALIBRATE] Updated ATM Model: tau={mean_tau:.3f}, T_surf={mean_tsurf:.1f}K, "
              f"T_atm={mean_tatm:.1f}K (model='{new_t_atm_model}').")
            print(f"              Mean New T_cal = {np.mean(new_t_cal):.2f} K")
        
            
    # --- 高速再スケール処理 (比率のブロードキャスト) ---
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = new_t_cal / old_t_cal
    ratio[~np.isfinite(ratio)] = 0.0
    data = data * ratio[:, None].astype(data.dtype)
    table["TCAL"] = new_t_cal
    
    res = Scantable(meta=meta, data=data, table=table)
    res.history = scantable.history.copy()
    res.history["notes"] = res.history.get("notes", "") + f" | Recalibrated to tau={new_tau}"
    
    return res
