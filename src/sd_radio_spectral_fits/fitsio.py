# src/sd_radio_spectral_fits/fitsio.py
from __future__ import annotations

from typing import Optional, Dict, List, Tuple, Union, Any
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time

# Temperature scale helpers (package-local; keep Single Source of Truth)
from .tempscale import (
    normalize_tempscal,
    ensure_tempscal_column,
    ensure_beameff_column,
    tempscal_array,
    beameff_array,
    require_beameff,
    ta_to_tr,
    tr_to_ta,
    convert_rowwise_vla,
    append_scale_history,
)


# NOTE: Single Source of Truth for FITS BinTable writing
try:
    from .sdfits_bintable import build_single_dish_table_hdu, build_history_hdu, set_meta_keyword
except Exception:  # pragma: no cover
    # Allow running this module standalone
    from sdfits_bintable import build_single_dish_table_hdu, build_history_hdu, set_meta_keyword

# =========================================================
# 1. データコンテナ: Scantable
# =========================================================
@dataclass
class Scantable:
    """
    SDFITSデータのメモリ上の表現。
    
    ヘッダー情報(meta)、スペクトルデータ(data)、観測テーブル(table)、
    および処理履歴(history)を一括管理するコンテナです。

    【VLA (可変長配列) サポートについて】
    `data` 属性には、通常の2次元NumPy配列 (N_row x N_chan) だけでなく、
    各行でチャンネル数が異なる1次元配列のリスト `List[np.ndarray]` を格納できます。
    これをそのまま `write_scantable` に渡すと、自動的にFITSのVLA (TFORM='PE'等) として
    安全にディスクへ保存されます。
    """
    
    meta: dict
    # [MODIFIED] Allow ragged arrays (List of arrays) or standard 2D array
    data: Union[np.ndarray, List[np.ndarray]]            
    table: pd.DataFrame         # 付随情報 (N_row, columns...)
    history: dict = field(default_factory=dict)

    def __repr__(self):
        # Handle list vs ndarray in repr
        if isinstance(self.data, np.ndarray):
            shape_str = f"rows={self.data.shape[0]}, chans={self.data.shape[1]}"
        else:
            shape_str = f"rows={len(self.data)}, chans=VLA"
            
        return (f"<Scantable {shape_str}, "
                f"history_keys={list(self.history.keys())}>")

    def copy(self) -> Scantable:
        """浅いコピーを作成して返します"""
        # Handle data copy for both types
        if isinstance(self.data, np.ndarray):
            d = self.data.copy()
        else:
            d = [x.copy() for x in self.data]
            
        return Scantable(
            meta=self.meta.copy(),
            data=d,
            table=self.table.copy(),
            history=self.history.copy()
        )


# =========================================================
# 2. 読み書き関数 (I/O) - Scantable Interface
# =========================================================

def read_scantable(path: str, *, tr_input_policy: str = "preserve") -> Scantable:
    """
    FITSファイルを読み込み、Scantable オブジェクトとして返します。
    """
    # 既存の読み込みロジックを利用
    meta, data, table, hist = read_tastar_fits(path)

    # =================================================================
    # [ADDED] レガシーカラムの自動マイグレーション (Backward Compatibility)
    # 過去に "V_CORR_KMS" 等で保存された速度列を、
    # 公開意味論 VELOSYS/VFRAME [m/s] へ安全に正規化する。
    # =================================================================
    rename_map = {}
    migrated_msgs = []
    legacy_vcorr_cols = []
    for col in table.columns:
        up = str(col).upper()
        if up == "V_CORR_KMS":
            legacy_vcorr_cols.append(col)
        elif up == "TAMB_K":
            rename_map[col] = "THOT"
        elif up == "TAU":
            rename_map[col] = "TAU0"

    if rename_map:
        table = table.rename(columns=rename_map)

    if legacy_vcorr_cols:
        # Use the first finite legacy column; multiple legacy columns must agree.
        legacy_ms = None
        for col in legacy_vcorr_cols:
            series = pd.to_numeric(table[col], errors="coerce")
            arr_ms = series.to_numpy(dtype=float) * 1000.0
            if legacy_ms is None:
                legacy_ms = arr_ms
            elif not np.allclose(arr_ms, legacy_ms, rtol=0.0, atol=1e-6, equal_nan=True):
                raise ValueError("Multiple V_CORR_KMS columns disagree; cannot migrate safely.")
        if legacy_ms is not None:
            if "VFRAME" not in table.columns:
                table["VFRAME"] = legacy_ms
                migrated_msgs.append("Converted legacy V_CORR_KMS [km/s] to VFRAME [m/s]")
            if "VELOSYS" not in table.columns:
                table["VELOSYS"] = legacy_ms
                migrated_msgs.append("Converted legacy V_CORR_KMS [km/s] to VELOSYS [m/s]")
        table = table.drop(columns=legacy_vcorr_cols)

    if "VFRAME" in table.columns and "VELOSYS" not in table.columns:
        table["VELOSYS"] = pd.to_numeric(table["VFRAME"], errors="coerce")
        migrated_msgs.append("Mirrored VFRAME into VELOSYS")

    if migrated_msgs:
        hist["fitsio_migration"] = "; ".join(migrated_msgs)

    # =================================================================
    
    # -----------------------------------------------------------------
    # Temperature scale normalization (internal policy)
    # - New Policy (2026): Non-destructive read.
    # - If on-disk TEMPSCAL indicates TR*, keep it as TR*.
    # - Do NOT automatically convert to TA*.
    # -----------------------------------------------------------------
    try:
        # Ensure columns exist but respect original values
        table = ensure_tempscal_column(table, default=normalize_tempscal(meta.get("TEMPSCAL", "TA*"), default="TA*"))
        table = ensure_beameff_column(table)
        
        # [MODIFIED] Removed the logic that forced TR* -> TA* conversion.
        # We now simply respect whatever is in the file.
        # Users can use Viewer to toggle display, or tools to process accordingly.
        
    except Exception as e:
        # Non-fatal: keep legacy behavior if scale metadata is absent/broken.
        # Still, try to set a sensible TEMPSCAL for downstream tools.
        try:
            if isinstance(table, pd.DataFrame) and "TEMPSCAL" not in table.columns:
                table = table.copy()
                table["TEMPSCAL"] = normalize_tempscal(meta.get("TEMPSCAL", "TA*"), default="TA*")
        except Exception:
            pass

    if hist is None:
        hist = {}

    return Scantable(
        meta=meta,
        data=data,
        table=table,
        history=hist
    )


def write_scantable(
    path: str,
    scantable: Scantable,
    spectrum_column: str = "DATA",
    overwrite: bool = True,
    **kwargs
) -> None:
    """
    Scantable オブジェクトを FITS ファイルとして保存します。
    """
    # 元のオブジェクトを壊さないよう、テーブルのみコピーして作業する
    df = scantable.table.copy()
    meta = scantable.meta
    
    # データ列へ格上げ(Promote)すべき重要キーワード
    promote_keys = [
        "RESTFRQ", "RESTFREQ", "CRVAL1", "CDELT1", "CRPIX1", "CTYPE1",
        "CUNIT1", "SPECSYS", "SSYSOBS", "VELDEF", "VELOSYS", "VFRAME",
        "TEMPSCAL", "BEAMEFF"
    ]
    
    for key in promote_keys:
        # DataFrameに列がなく、Headerに値がある場合に追加
        if (key not in df.columns) and (key in meta):
            val = meta[key]
            # 値を全行に定数として埋める（型推論はPandasに任せる）
            df[key] = val
            
    # 修正したDataFrameを使って書き込み
    write_sdfits(
        out_path=path,
        meta=meta,
        data=scantable.data,
        table=df,
        history=scantable.history,
        spectrum_column=spectrum_column,
        overwrite=overwrite,
        **kwargs
    )


# =========================================================
# 3. 内部ロジック / 旧インターフェース
# =========================================================

def _get_any(d: dict, *keys):
    for k in keys:
        if k in d and d[k] not in (None, ""):
            return d[k]
    return None


def _dedup_case_insensitive_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate columns that would collide in FITS after uppercasing.
    """
    cols = list(df.columns)
    prefer_keep = {"TIMESTAMP": "TIMESTAMP"}

    drop = set()
    up_to_originals: Dict[str, List[str]] = {}
    for c in cols:
        up = str(c).upper()
        up_to_originals.setdefault(up, []).append(c)

    for up, keep_name in prefer_keep.items():
        if up in up_to_originals and keep_name in up_to_originals[up]:
            for c in up_to_originals[up]:
                if c != keep_name:
                    drop.add(c)

    ordered_keep: List[str] = []
    kept_upper = set()
    for c in cols:
        if c in drop:
            continue
        up = str(c).upper()
        if up in kept_upper:
            continue
        kept_upper.add(up)
        ordered_keep.append(c)

    return df.loc[:, ordered_keep].copy()


def _normalize_meta_wcs_and_rest(meta: dict) -> dict:
    """Normalize spectral WCS/rest keywords while keeping FITS-standard keys primary."""
    m = dict(meta)

    rest = None
    for k in ("RESTFRQ", "RESTFREQ"):
        if k in m and m[k] not in (None, ""):
            rest = float(m[k])
            break
    if rest is not None:
        m["RESTFRQ"] = rest
        m["RESTFREQ"] = rest
        m.setdefault("rest_hz", rest)

    if "CTYPE1" in m and m["CTYPE1"] not in (None, ""):
        ctype1 = str(m["CTYPE1"]).strip().upper()
        if ctype1.startswith("FREQ"):
            m["CTYPE1"] = "FREQ"

    if "CUNIT1" in m and m["CUNIT1"] not in (None, ""):
        cunit1 = str(m["CUNIT1"]).strip().lower().replace(" ", "")
        if cunit1 in ("hz", "khz", "mhz", "ghz"):
            # Reader normalization keeps the keyword spelling simple;
            # numeric axis conversion is handled in the calibrator where needed.
            if cunit1 == "hz":
                m["CUNIT1"] = "Hz"

    def _norm_specsys(value):
        s = str(value).strip().upper().replace(" ", "")
        aliases = {"TOPO": "TOPOCENT", "TOPOCENTER": "TOPOCENT", "LSR": "LSRK"}
        return aliases.get(s, s)

    if "SPECSYS" in m and m["SPECSYS"] not in (None, ""):
        m["SPECSYS"] = _norm_specsys(m["SPECSYS"])
    if "SSYSOBS" in m and m["SSYSOBS"] not in (None, ""):
        m["SSYSOBS"] = _norm_specsys(m["SSYSOBS"])
    elif "SPECSYS" in m and m["SPECSYS"] not in (None, ""):
        m["SSYSOBS"] = m["SPECSYS"]
    elif "SSYSOBS" in m and m["SSYSOBS"] not in (None, ""):
        m["SPECSYS"] = m["SSYSOBS"]

    return m
    


def _to_native_endian_array(arr: np.ndarray) -> np.ndarray:
    """Return a native-endian view/copy of a NumPy array when needed."""
    a = np.asarray(arr)
    dt = a.dtype
    if dt.kind == "O" or not hasattr(dt, "byteorder"):
        return a
    bo = dt.byteorder
    if bo in ("=", "|"):
        return a
    native = "<" if np.little_endian else ">"
    if bo == native:
        return a
    return a.byteswap().view(dt.newbyteorder("="))


def _df_to_native_endian(df: pd.DataFrame) -> pd.DataFrame:
    """Convert numeric big-endian columns to native-endian for pandas ops."""
    out = df.copy()
    for c in out.columns:
        v = out[c]
        try:
            arr = np.asarray(v)
        except Exception:
            continue
        if getattr(arr, "ndim", 1) != 1:
            continue
        if arr.dtype.kind in ("b", "i", "u", "f", "c", "m", "M"):
            out[c] = _to_native_endian_array(arr)
    return out


def _resolve_table_timestamps(table: pd.DataFrame) -> Optional[pd.DatetimeIndex]:
    """Resolve per-row UTC timestamps from standard/legacy time columns.

    Priority:
    1. TIMESTAMP
    2. MJD
    3. DATE-OBS / DATEOBS (+ TIME seconds offset if present)
    4. legacy TIME (MJD days or Unix seconds)
    """
    if table is None or len(table) == 0:
        return pd.DatetimeIndex([])

    if "TIMESTAMP" in table.columns:
        ts = pd.to_datetime(table["TIMESTAMP"].astype(str), utc=True, errors="coerce")
        if ts.notna().any():
            return pd.DatetimeIndex(ts)

    if "MJD" in table.columns:
        arr = pd.to_numeric(table["MJD"], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(arr).any():
            return pd.DatetimeIndex(pd.to_datetime(Time(arr, format="mjd", scale="utc").to_datetime(), utc=True))

    for dcol in ("DATE-OBS", "DATEOBS"):
        if dcol in table.columns:
            base = pd.to_datetime(table[dcol].astype(str), utc=True, errors="coerce")
            if base.notna().any():
                ts = pd.DatetimeIndex(base)
                if "TIME" in table.columns:
                    offs = pd.to_numeric(table["TIME"], errors="coerce").to_numpy(dtype=float)
                    if np.isfinite(offs).any():
                        # Standard SDFITS semantics: TIME is seconds from DATE-OBS.
                        ts = pd.DatetimeIndex(ts + pd.to_timedelta(np.nan_to_num(offs, nan=0.0), unit="s"))
                return ts

    if "TIME" in table.columns:
        arr = pd.to_numeric(table["TIME"], errors="coerce").to_numpy(dtype=float)
        finite = arr[np.isfinite(arr)]
        if finite.size:
            if float(np.nanmedian(finite)) < 1.0e6:
                ts = pd.to_datetime(Time(arr, format="mjd", scale="utc").to_datetime(), utc=True)
            else:
                ts = pd.to_datetime(arr, unit="s", utc=True, errors="coerce")
            return pd.DatetimeIndex(ts)

    return None

def _get_single_dish_bintable(hdul):
    for hdu in hdul:
        if isinstance(hdu, fits.BinTableHDU) and hdu.name.upper() == "SINGLE DISH":
            return hdu
    return None


def read_tastar_fits(path: str):
    """
    Read Ta* dumps FITS. (Supports VLA 'PE' format)
    """
    with fits.open(path) as hdul:
        meta = dict(hdul[0].header)

        # --- (B) SDFITS-like: SINGLE DISH BinTable ---
        if "SINGLE DISH" in hdul:
            hdu = hdul["SINGLE DISH"]
            
            header_dict = dict(hdu.header)
            for k in ["NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT"]:
                if k in header_dict:
                    del header_dict[k]
            meta.update(header_dict)
            
            meta = _normalize_meta_wcs_and_rest(meta)

            tab_hdu = _get_single_dish_bintable(hdul)
            tab = None if tab_hdu is None else tab_hdu.data
            
            if tab is None:
                raise ValueError("SINGLE DISH table is empty.")

            spec_col = None
            for cand in ("DATA", "SPECTRUM"):
                if cand in tab.dtype.names:
                    spec_col = cand
                    break
            if spec_col is None:
                raise KeyError("SINGLE DISH table has neither DATA nor SPECTRUM column.")

            # [MODIFIED] Read data: handle potential object array (VLA)
            raw_data = tab[spec_col]
            
            # If it's VLA, astropy returns a structured array of objects or similar
            # Try to convert to float array if possible (fixed length), otherwise keep as object
            try:
                # If fixed length, this works
                data = np.asarray(raw_data, dtype=np.float32)
            except (ValueError, TypeError):
                # VLA case: keep as object array or list of arrays
                data = raw_data
            
            # Basic NCHAN estimation from first row if possible
            if isinstance(data, np.ndarray) and data.ndim == 2:
                nchan_loaded = data.shape[-1]
            elif len(data) > 0:
                nchan_loaded = len(data[0])
            else:
                nchan_loaded = 0
                
            meta["NAXIS1"] = nchan_loaded
            meta["NCHAN"] = nchan_loaded
            meta["NCHANSEL"] = nchan_loaded

            exclude = {spec_col}
            for cand in ("FLAG", "FLAGS"):
                if cand in tab.dtype.names:
                    exclude.add(cand)

            cols = {}
            for name in tab.dtype.names:
                if name in exclude:
                    continue
                arr = np.asarray(tab[name])
                
                # --- 多次元カラムの読み込み対応 (BSL_COEFなど) ---
                if getattr(arr, "ndim", 1) > 1:
                    # 2次元以上の配列はリストに変換してDataFrameに格納する
                    cols[name] = list(arr)
                    continue
                # ----------------------------------------------------

                if arr.dtype.kind == "S":
                    arr = arr.astype(str)
                cols[name] = arr
            table = _df_to_native_endian(pd.DataFrame(cols))

            ts = _resolve_table_timestamps(table)
            if ts is not None and len(ts) == len(table):
                table = table.set_index(ts)
                table.index.name = "TIMESTAMP"

            history = {}
            if "HISTORY" in hdul:
                htab = hdul["HISTORY"].data
                if htab is not None and "KEY" in htab.dtype.names and "VALUE" in htab.dtype.names:
                    for k, v in zip(htab["KEY"], htab["VALUE"]):
                        kk = (k.decode() if isinstance(k, (bytes, bytearray)) else str(k)).strip()
                        vv = (v.decode() if isinstance(v, (bytes, bytearray)) else str(v)).strip()
                        history[kk] = vv

            return meta, data, table, history

        # --- (A) Internal format: DATA Image + DUMPS table ---
        if "DATA" in hdul:
            meta.update(dict(hdul["DATA"].header))
        if "DUMPS" in hdul:
            meta.update(dict(hdul["DUMPS"].header))

        meta = _normalize_meta_wcs_and_rest(meta)
        
        if "DATA" in hdul:
            data = np.asarray(hdul["DATA"].data)
        else:
            data_hdu = None
            for h in hdul:
                if type(h).__name__ == "ImageHDU":
                    data_hdu = h
                    break
            if data_hdu is None or data_hdu.data is None:
                raise KeyError("Missing DATA ImageHDU in Ta* FITS.")
            data = np.asarray(data_hdu.data)

        if "DUMPS" not in hdul:
            raise KeyError("Missing DUMPS BinTableHDU in Ta* FITS.")
        tab = hdul["DUMPS"].data
        cols = {}
        for name in tab.dtype.names:
            cols[name] = np.asarray(tab[name])
        table = _df_to_native_endian(pd.DataFrame(cols))

        ts = _resolve_table_timestamps(table)
        if ts is not None and len(ts) == len(table):
            table = table.set_index(ts)
            table.index.name = "TIMESTAMP"

        history = {}
        if "HISTORY" in hdul:
            htab = hdul["HISTORY"].data
            if htab is not None and "KEY" in htab.dtype.names and "VALUE" in htab.dtype.names:
                for k, v in zip(htab["KEY"], htab["VALUE"]):
                    kk = (k.decode() if isinstance(k, (bytes, bytearray)) else str(k)).strip()
                    vv = (v.decode() if isinstance(v, (bytes, bytearray)) else str(v)).strip()
                    history[kk] = vv

    return meta, data, table, history


def write_sdfits(
    out_path: str,
    meta: dict,
    data: Union[np.ndarray, List[np.ndarray]],
    table: pd.DataFrame,
    history: dict | None = None,
    *,
    spectrum_column: str = "DATA",
    include_flag: bool = True,
    overwrite: bool = True,
    string_widths: dict[str, int] | None = None,
    **kwargs,
) -> None:
    """
    Write CASA-friendly SDFITS-like output.

    This is the recommended write API going forward.

    Parameters (main ones)
    ----------------------
    out_path:
        Output FITS path.
    meta:
        Global header-like metadata (PRIMARY). Non-FITS keys are moved into HISTORY as `META:<KEY>`.
    data:
        Spectrum per row. Either a 2D ndarray (fixed length) or list of 1D arrays (ragged OK; VLA used).
        - チャンネル数が固定の場合: 2次元の `np.ndarray`
        - チャンネル数が変動する場合 (VLA): 1次元の `np.ndarray` を格納した `List`
          (FITS出力時に自動で TFORM='PE' などの可変長配列として処理されます)
    table:
        Per-row metadata (pandas DataFrame). Extra columns are preserved. Vector-in-cell columns
        (list/ndarray in each cell) are automatically written as fixed vectors or VLA.
    history:
        Processing history to write to HISTORY extension (key/value table).
    spectrum_column:
        FITS column name for spectrum (default 'DATA').
    include_flag:
        If True, write FLAG column (created as zeros if not provided elsewhere).
    overwrite:
        If True, overwrite existing file.
    string_widths:
        Optional fixed widths for specific string columns (e.g. {'OBJECT': 32}).
        Keys are case-insensitive.

    Notes
    -----
    - This function delegates the BinTable construction to `sdfits_bintable.build_single_dish_table_hdu()`
      so that both "analysis I/O" (fitsio) and "observing/convert writer" (sdfits_writer) share the same
      FITS column inference and VLA/fixed-length logic.
    """
    # Keep kwargs for backward compatibility (some callers may pass future options).
    _ = kwargs

    meta = _normalize_meta_wcs_and_rest(meta)

    # ---- ensure TIME/DATEOBS columns ----
    df = table.copy()
    # -----------------------------------------------------------------
    # TEMPSCAL / BEAMEFF handling (SDFITS/CASA compatibility)
    # - `tempscal` (or `out_scale`) declares the on-disk DATA scale ("TA*"|"TR*").
    # - `data_scale` declares the in-memory scale of `data` (default "TA*").
    # - If conversion is required, it is applied only to the written data (non-destructive),
    #   and recorded to HISTORY.
    # -----------------------------------------------------------------
    desired_scale = normalize_tempscal(kwargs.pop("tempscal", kwargs.pop("out_scale", meta.get("TEMPSCAL", "TA*"))), default="TA*")
    data_scale = normalize_tempscal(kwargs.pop("data_scale", "TA*"), default="TA*")

    df = ensure_tempscal_column(df, default=normalize_tempscal(meta.get("TEMPSCAL", "TA*"), default="TA*"))
    df = ensure_beameff_column(df)

    # Force on-disk declaration
    df["TEMPSCAL"] = desired_scale
    meta["TEMPSCAL"] = desired_scale

    # Convert data for writing if needed
    data_to_write = data
    if desired_scale != data_scale:
        n_rows = len(df)
        beameff = beameff_array(df, meta, n_rows)
        require_beameff(beameff, on_fail=str(kwargs.get("beameff_on_fail", "error")))
        if isinstance(data, list):
            # VLA list: rowwise conversion
            if desired_scale == "TR*" and data_scale == "TA*":
                data_to_write = convert_rowwise_vla(list(data), beameff, direction="ta_to_tr")
            elif desired_scale == "TA*" and data_scale == "TR*":
                data_to_write = convert_rowwise_vla(list(data), beameff, direction="tr_to_ta")
        else:
            arr = np.asarray(data, dtype=float)
            if desired_scale == "TR*" and data_scale == "TA*":
                data_to_write = ta_to_tr(arr, beameff).astype(np.float32, copy=False)
            elif desired_scale == "TA*" and data_scale == "TR*":
                data_to_write = tr_to_ta(arr, beameff).astype(np.float32, copy=False)

        if history is None:
            history = {}
        history = append_scale_history(history, {
            "stage": "write_sdfits",
            "action": "convert_on_write",
            "data_scale": data_scale,
            "tempscal": desired_scale,
            "formula": ("Tr* = Ta* / BEAMEFF" if desired_scale == "TR*" else "Ta* = Tr* * BEAMEFF"),
        })
    ts = None
    if isinstance(df.index, pd.DatetimeIndex):
        ts = pd.DatetimeIndex(df.index)
    elif "TIMESTAMP" in df.columns:
        try:
            ts = pd.DatetimeIndex(pd.to_datetime(df["TIMESTAMP"], utc=True, errors="coerce"))
        except Exception:
            ts = None
    elif "MJD" in df.columns:
        try:
            arr = pd.to_numeric(df["MJD"], errors="coerce").to_numpy(dtype=float)
            ts = pd.DatetimeIndex(pd.to_datetime(Time(arr, format="mjd", scale="utc").to_datetime(), utc=True))
        except Exception:
            ts = None
    elif "DATE-OBS" in df.columns or "DATEOBS" in df.columns or "TIME" in df.columns:
        ts = _resolve_table_timestamps(df)

    if ts is not None and len(ts) == len(df):
        dateobs_iso = np.array([t.isoformat(timespec="milliseconds") for t in ts.to_pydatetime()], dtype="U")
        mjd = Time(ts.to_pydatetime(), scale="utc").mjd.astype(np.float64)
        if "DATE-OBS" not in df.columns:
            df["DATE-OBS"] = dateobs_iso
        if "DATEOBS" not in df.columns:
            df["DATEOBS"] = dateobs_iso
        if "TIMESTAMP" not in df.columns:
            df["TIMESTAMP"] = dateobs_iso
        if "MJD" not in df.columns:
            df["MJD"] = mjd
        if "TIME" not in df.columns:
            # Standard semantics: seconds from DATE-OBS. Because DATE-OBS is stored
            # here as the full per-row UTC timestamp, TIME is zero for each row.
            df["TIME"] = np.zeros(len(df), dtype=np.float64)

    # ---- build SINGLE DISH table via core ----
    tbl, hdr_nchan = build_single_dish_table_hdu(
        table=df,
        data=data_to_write,
        spectrum_column=spectrum_column,
        include_flag=include_flag,
        flag=None,
        string_widths=string_widths,
        normalize_columns=True,
        bunit="K",
        extname="SINGLE DISH",
    )

    # ---- headers ----
    hdr = fits.Header()
    hist_extra: dict[str, str] = {}

    for k, v in meta.items():
        if v is None:
            continue
        ku = str(k).strip()
        kuu = ku.upper()
        if kuu in ("EXTNAME", "EXTVER", "XTENSION"):
            continue
        try:
            if isinstance(v, tuple) and len(v) >= 2:
                val, com = v[0], v[1]
            else:
                val, com = v, None
            set_meta_keyword(hdr, kuu, val, comment=(str(com) if com is not None else None))
        except Exception as e:
            # keep going; record into HISTORY as a diagnostic
            hist_extra[f"META:{kuu}"] = str(v)

    for k in ["NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT"]:
        if k in hdr:
            del hdr[k]

    hdr["NCHAN"] = hdr_nchan
    hdr["NCHANSEL"] = hdr_nchan

    pri = fits.PrimaryHDU(header=hdr)

    # Promote a subset of keys to SINGLE DISH header if missing there.
    for k in (
        "SITELAT", "SITELONG", "SITELON", "SITEELEV",
        "OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z",
        "CTYPE1", "CUNIT1", "CRVAL1", "CDELT1", "CRPIX1",
        "RESTFRQ", "RESTFREQ", "SPECSYS", "SSYSOBS", "VELDEF", "VELOSYS",
        "TIMESYS", "RADESYS", "EQUINOX",
        "NCHAN", "NCHANSEL",
    ):
        if k in pri.header and k not in tbl.header:
            tbl.header[k] = pri.header[k]

    if hist_extra:
        if history is None:
            history = {}
        if isinstance(history, dict):
            history.update(hist_extra)
        else:
            history = {"_history": str(history), **hist_extra}

    tbl.header["BUNIT"] = "K"

    # ---- HISTORY ----
    hdus: list[fits.hdu.base.ExtensionHDU] = [pri, tbl]

    hist_hdu = build_history_hdu(history)
    if hist_hdu is not None:
        hdus.append(hist_hdu)

    fits.HDUList(hdus).writeto(out_path, overwrite=overwrite)


