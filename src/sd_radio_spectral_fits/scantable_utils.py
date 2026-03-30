# src/sd_radio_spectral_fits/scantable_utils.py
from __future__ import annotations

"""
scantable_utils.py
==================

Scantableの検査、変更、結合、検索、座標計算を行うユーティリティモジュール。
大規模データ(10万行〜)を想定し、Pandas/Numpy/Astropyによるベクトル化処理で高速化されています。

主な機能:
- Inspection: show_scantable, describe_columns
- Coordinates: calc_mapping_offsets (Immutable: 元データを変更せず計算結果を返す)
- Merge: merge_scantables
- Modification: update_metadata (Safety checks included)
- Query: find_scans, filter_scantable

Dependencies:
    pandas, numpy, astropy
    .fitsio (Scantable, read_scantable)
"""

from typing import Any, Dict, List, Optional, Sequence, Union, Set
from collections.abc import Mapping
import warnings
import numpy as np
import pandas as pd

from astropy.coordinates import SkyCoord
import astropy.units as u

# fitsioからのインポート
# 循環参照を防ぐため、本来はTYPE_CHECKINGブロックを使うか遅延インポートが望ましいですが、
# 実行時にクラス定義が必要なためトップレベルでインポートします。
from .fitsio import Scantable, read_scantable


# =========================================================
# Configuration (設定)
# =========================================================

# show_scantableでデフォルト表示されるカラムと説明
DEFAULT_SHOW_COLS: Dict[str, str] = {
    "SCAN": "Scan ID",
    "TIMESTAMP": "Resolved UTC timestamp",
    "OBSMODE": "Observation Mode",
    "OBJECT": "Target Name",
    "OFS_LON": "Map Offset X [arcsec] (calc)",
    "OFS_LAT": "Map Offset Y [arcsec] (calc)",
    "EXPOSURE": "Integration Time [s]",
    "TEMPSCAL": "Temperature Scale",
}

# 変更時に警告/ブロックするカラム（座標系やデータ本体など）
DANGER_COLS: Set[str] = {
    "CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1",
    "RESTFREQ", "RESTFRQ", "SPECSYS", "DATA", "FLAG", "SPECTRUM"
}

# バリデーション用の許容値リスト
ENUM_VALS: Dict[str, Set[str]] = {
    "OBSMODE": {"ON", "OFF", "HOT", "SKY", "FS_SIG", "FS_REF", "CAL_ON", "CAL_OFF", "ZERO", "UNKNOWN"},
    "TEMPSCAL": {"TA*", "TR*", "TMB", "Tu"},
    "CALSTAT": {"RAW", "TASTAR", "TA", "TMB", "INTEGRATED"},
}

# カラム名のエイリアス（同義語）定義
# 読み込み時・書き込み時に相互に参照・同期するために使用
COLUMN_ALIASES: Dict[str, List[str]] = {
    "RESTFREQ": ["RESTFRQ"],
    "RESTFRQ": ["RESTFREQ"],
    # 必要に応じて他のエイリアスも追加可能 (例: "OBSRA": ["RA"] など)
}


# =========================================================
# Endianness Helpers (FITS由来のBig-endian対策)
# =========================================================

def _to_native_endian_array(a: np.ndarray) -> np.ndarray:
    """
    NumPy 2.x 互換: ndarray.newbyteorder() は廃止されたため dtype.newbyteorder() を使う。

    Big-endian 配列を native endian に変換して返す。
    (byteswap はコピーを伴う)
    """
    a = np.asarray(a)
    dt = a.dtype
    # '|' はエンディアン非依存 (bytes/str など)
    if dt.byteorder == "|" or dt.isnative:
        return a
    # byteswap して、native byteorder の dtype を view する
    return a.byteswap().view(dt.newbyteorder("="))


def _df_to_native_endian(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pandas が内部で take_nd を行う際に Big-endian バッファを扱えず
    `ValueError: Big-endian buffer not supported on little-endian compiler`
    が出るケースの対策。

    - Big-endian な数値カラム(i/u/f/c)だけを native endian に置き換えた DataFrame を返す。
    - 変換が不要なら元の df をそのまま返す（余計なコピーをしない）。

    注意: 元の DataFrame を変更しない（浅いコピーに対して列を置換）。
    """
    need = False
    for col in df.columns:
        dt = getattr(df[col], "dtype", None)
        if isinstance(dt, np.dtype) and dt.kind in ("i", "u", "f", "c") and (not dt.isnative):
            need = True
            break
    if not need:
        return df

    df2 = df.copy(deep=False)
    for col in df2.columns:
        s = df2[col]
        dt = getattr(s, "dtype", None)
        if not isinstance(dt, np.dtype):
            continue
        if dt.kind in ("i", "u", "f", "c") and (not dt.isnative):
            a = s.to_numpy(copy=False)
            if isinstance(a, np.ndarray):
                df2[col] = _to_native_endian_array(a)
    return df2


def _series_to_datetime_utc(series: pd.Series) -> pd.DatetimeIndex:
    """Best-effort conversion of a table column to UTC datetime.

    Supports:
    - datetime-like columns
    - ISO strings
    - legacy numeric TIME/MJD columns (MJD days or Unix seconds)
    """
    if series.empty:
        return pd.DatetimeIndex([])
    dtype = getattr(series, "dtype", None)
    if pd.api.types.is_datetime64_any_dtype(dtype):
        return pd.DatetimeIndex(pd.to_datetime(series, utc=True, errors="coerce"))
    if pd.api.types.is_numeric_dtype(dtype):
        arr = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
        finite = arr[np.isfinite(arr)]
        if finite.size > 0 and float(np.nanmedian(np.abs(finite))) < 1.0e6:
            return pd.to_datetime(arr, unit="D", origin=pd.Timestamp("1858-11-17"), utc=True, errors="coerce")
        return pd.to_datetime(arr, unit="s", utc=True, errors="coerce")
    return pd.DatetimeIndex(pd.to_datetime(series, utc=True, errors="coerce"))


def _resolve_table_timestamps(table: pd.DataFrame) -> Optional[pd.DatetimeIndex]:
    """Resolve per-row UTC timestamps from a Scantable-like table.

    Priority:
    1. TIMESTAMP
    2. MJD
    3. DATE-OBS / DATEOBS (+ TIME seconds if DATE column is midnight-only)
    4. legacy TIME (MJD days or Unix seconds)
    5. DatetimeIndex
    """
    tab = _df_to_native_endian(table)

    if "TIMESTAMP" in tab.columns:
        ts = _series_to_datetime_utc(tab["TIMESTAMP"])
        if not ts.isna().all():
            return ts

    if "MJD" in tab.columns:
        ts = _series_to_datetime_utc(tab["MJD"])
        if not ts.isna().all():
            return ts

    for dcol in ("DATE-OBS", "DATEOBS"):
        if dcol not in tab.columns:
            continue
        base = _series_to_datetime_utc(tab[dcol])
        if base.isna().all():
            continue
        if "TIME" in tab.columns:
            tnum = pd.to_numeric(tab["TIME"], errors="coerce").to_numpy(dtype=float)
            finite = np.isfinite(tnum)
            if finite.any():
                non_na = ~base.isna()
                midnight_only = bool(non_na.any() and (base[non_na] == base[non_na].normalize()).all())
                if midnight_only:
                    base = base + pd.to_timedelta(np.where(finite, tnum, 0.0), unit="s")
        return base

    if "TIME" in tab.columns:
        ts = _series_to_datetime_utc(tab["TIME"])
        if not ts.isna().all():
            return ts

    if isinstance(tab.index, pd.DatetimeIndex):
        ts = tab.index
        if ts.tz is None:
            ts = ts.tz_localize("UTC")
        else:
            ts = ts.tz_convert("UTC")
        return pd.DatetimeIndex(ts)

    return None


# =========================================================
# 1. Inspection (確認・表示)
# =========================================================

def describe_columns(sc: Scantable) -> None:
    """
    Scantableに含まれる全カラム名、データ型、および説明を表示する。
    また、Scantableのメタデータ（ヘッダー情報）に含まれるキーと値も一覧表示する。
    """
    df = _df_to_native_endian(sc.table)
    print(f"--- Columns in Scantable ({len(df)} rows) ---")
    print(f"{'Column Name':<15} | {'Dtype':<10} | Description")
    print("-" * 70)
    for col in df.columns:
        desc = DEFAULT_SHOW_COLS.get(col, "")
        dtype = str(df[col].dtype)
        print(f"{col:<15} | {dtype:<10} | {desc}")
    print("-" * 70)

    print(f"--- Metadata (Header) ---")
    if not sc.meta:
        print("(No metadata found)")
    else:
        print(f"{'Key':<15} | Value")
        print("-" * 70)
        for k, v in sc.meta.items():
            print(f"{k:<15} | {v}")
    print("-" * 70)


def show_scantable(
    inputs: Union[Scantable, str, Sequence[Union[Scantable, str]]],
    rows: Union[str, slice, List[int], None] = None,
    columns: Union[List[str], str] = "default",
    head: Optional[int] = 20,
    show_legend: bool = False,
    # --- 拡張データ用 ---
    extra_data: Optional[pd.DataFrame] = None,
    # --- その場で座標計算するためのショートカット ---
    ref_coord: Optional[Union[str, SkyCoord]] = None,
    frame: str = "ICRS",
    projection: str = "GLS",
    unit: str = "arcsec",
) -> None:
    """
    Scantableの中身を表形式で確認する。
    ref_coordを指定すると、その場でオフセットを計算して表示する（データは変更しない）。
    calc_mapping_offsetsで計算済みの外部データ(extra_data)を一時的に結合して表示することも可能。

    テーブルに含まれないカラム（例: RESTFREQ）が指定された場合、以下の順序で検索し、
    最初に見つかった値を全行に結合して表示する。

    1. テーブル内の同名カラム
    2. テーブル内のエイリアスカラム (例: RESTFREQ要求 -> RESTFRQ確認)
    3. メタデータ(ヘッダー)内の同名キー
    4. メタデータ内のエイリアスキー
    """
    # 入力の正規化（リスト化）
    if isinstance(inputs, (Scantable, str)):
        target_list: List[Union[Scantable, str]] = [inputs]
    else:
        target_list = list(inputs)

    for i, target in enumerate(target_list):
        # ロードまたはオブジェクト取得
        if isinstance(target, str):
            print(f"\n>>> Loading: {target} ...")
            sc = read_scantable(target)
            name = target
        else:
            sc = target
            name = f"Object #{i}"

        # DataFrameの準備（表示用: Big-endian を native にしておく）
        df_view = _df_to_native_endian(sc.table).copy()
        # Resolve a reliable UTC timestamp column from TIMESTAMP / MJD / DATE-OBS / TIME.
        if "TIMESTAMP" not in df_view.columns:
            ts_resolved = _resolve_table_timestamps(df_view)
            if ts_resolved is not None and not ts_resolved.isna().all():
                df_view.insert(0, "TIMESTAMP", ts_resolved)

        # 1. オンザフライ座標計算 (ref_coordがある場合)
        calc_df = None
        if ref_coord is not None:
            try:
                calc_df = calc_mapping_offsets(
                    sc, ref_coord=ref_coord, frame=frame,
                    projection=projection, unit=unit, verbose=False
                )
            except Exception as e:
                print(f"Warning: Coordinate calc failed: {e}")

        # 2. 外部データの結合（計算結果 or extra_data）
        extras = []
        if extra_data is not None:
            extras.append(extra_data)
        if calc_df is not None:
            extras.append(calc_df)

        if extras:
            current_len = len(df_view)
            to_concat = [df_view.reset_index(drop=True)]
            valid_extras = []
            for e in extras:
                if len(e) == current_len:
                    valid_extras.append(e.reset_index(drop=True))
                else:
                    print(f"Warning: Extra data length ({len(e)}) mismatch with table ({current_len}). Skipping.")
            if valid_extras:
                to_concat.extend(valid_extras)
                df_view = pd.concat(to_concat, axis=1)

        # 3. カラム選択ロジック
        available_cols = list(df_view.columns)

        # ユーザー指定カラムを大文字化して扱う（FITS系は通常大文字）
        desired_cols: List[str] = []
        if columns == "default":
            defaults = list(DEFAULT_SHOW_COLS.keys())
            if (ref_coord is not None) or (extra_data is not None):
                if "OFS_LON" not in defaults:
                    defaults.append("OFS_LON")
                if "OFS_LAT" not in defaults:
                    defaults.append("OFS_LAT")
            desired_cols = defaults
        elif columns == "all":
            desired_cols = available_cols
        else:
            if isinstance(columns, str):
                desired_cols = [c.strip().upper() for c in columns.split(",") if c.strip()]
            else:
                desired_cols = [str(c).strip().upper() for c in columns if str(c).strip()]

        # テーブルにないが、メタデータやエイリアスにあるカラムを検索して結合
        # (desired_cols は原則として大文字)
        for col in desired_cols:
            if col in available_cols:
                continue

            found_val = None
            found_source = None  # 'table_alias', 'meta', 'meta_alias'

            aliases = COLUMN_ALIASES.get(col, [])

            # 1. テーブル内のエイリアス検索
            for alias in aliases:
                if alias in available_cols:
                    found_val = df_view[alias]  # Series
                    found_source = "table_alias"
                    break

            if found_source is None:
                # 2. メタデータ(Header)検索
                if col in sc.meta:
                    found_val = sc.meta[col]
                    found_source = "meta"
                else:
                    # 3. メタデータ内のエイリアス検索
                    for alias in aliases:
                        if alias in sc.meta:
                            found_val = sc.meta[alias]
                            found_source = "meta_alias"
                            break

            if found_source is not None:
                df_view[col] = found_val
                available_cols.append(col)

        missing = [c for c in desired_cols if c not in available_cols]
        if missing and columns not in ("default", "all"):
            print(f"Warning: Columns not found in {name}: {missing}")

        use_cols = [c for c in desired_cols if c in available_cols]

        # 4. 行選択ロジック
        idxs = _parse_row_selector(rows, len(df_view))
        subset = df_view.iloc[idxs]

        print(f"\n=== {name} (Selected {len(subset)} / {len(df_view)} rows) ===")

        if show_legend:
            print("[Column Legend]")
            for c in use_cols:
                desc = DEFAULT_SHOW_COLS.get(c, "(No description)")
                print(f"  * {c:<10}: {desc}")
            print("-" * 40)

        if not use_cols:
            print("(No columns selected to display)")
        else:
            display_df = subset[use_cols].copy()

            is_already_row_numbers = False
            if pd.api.types.is_integer_dtype(display_df.index):
                if np.array_equal(display_df.index.values, idxs):
                    is_already_row_numbers = True

            if not is_already_row_numbers:
                idx_name = display_df.index.name
                if idx_name is None:
                    if isinstance(display_df.index, pd.DatetimeIndex):
                        idx_name = "TIMESTAMP"
                    else:
                        idx_name = "Index"

                if idx_name not in display_df.columns:
                    display_df.insert(0, idx_name, display_df.index)

                display_df.index = idxs
                display_df.index.name = "idx"

            if head is not None and len(display_df) > head:
                print(display_df.head(head).to_string(index=True))
                print(f"... (truncating to first {head} rows) ...")
            else:
                print(display_df.to_string(index=True))


# =========================================================
# 2. Coordinates (Immutable)
# =========================================================

def calc_mapping_offsets(
    sc: Scantable,
    ref_coord: Optional[Union[str, SkyCoord]] = None,
    frame: str = "ICRS",
    projection: str = "GLS",
    unit: str = "arcsec",
    cos_mode: str = "point",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    指定された座標系・投影法に基づいてマッピング用オフセット (OFS_LON, OFS_LAT) を計算する。

    【重要】Scantable自体は変更しません。計算結果を新しいDataFrameとして返します。

    Parameters
    ----------
    sc : Scantable
        入力データ。RA/DEC, GLON/GLAT, または AZ/EL 列が必要。
    ref_coord : str, SkyCoord, optional
        参照点（中心座標）。Noneの場合はOBSRA/DECまたはデータの平均値を使用。
    frame : str
        計算に使用する座標系 ("ICRS", "FK5", "Galactic", "AltAz" 等)。
    projection : str
        投影法。
        - "GLS" (Global Sinusoidal): dlon*cos(lat) を使用
        - "TAN" (Gnomonic): tangent plane（SkyCoord.spherical_offsets_to を使用）
        - "SIN": simple difference (dlon, dlat)
    unit : str
        出力単位 ("arcsec", "arcmin", "deg")。
    cos_mode : str
        GLS時の cos(lat) をどこで評価するか。
        - "point": 各点の lat (従来と同じ)
        - "ref"  : 参照点の lat0（小領域近似で一般的）
    verbose : bool
        詳細表示。

    Returns
    -------
    pd.DataFrame
        columns=["OFS_LON", "OFS_LAT"] を持つDataFrame。indexはsc.tableと同じ。
    """
    df = _df_to_native_endian(sc.table)

    # 1. 座標配列の作成 (ベクトル化)
    #    FITS 8文字制限等で AZ/EL 系のカラム名が揺れることがあるため、候補を広く見る。
    lon_col = lat_col = None
    src_frame = None

    if {"RA", "DEC"}.issubset(df.columns):
        lon_col, lat_col = "RA", "DEC"
        src_frame = "icrs"
    elif {"GLON", "GLAT"}.issubset(df.columns):
        lon_col, lat_col = "GLON", "GLAT"
        src_frame = "galactic"
    else:
        # AltAz / raw az-el
        for cand_lon, cand_lat in (("AZIMUTH", "ELEVATIO"), ("AZIMUTH", "ELEVATION"), ("AZ", "EL")):
            if {cand_lon, cand_lat}.issubset(df.columns):
                lon_col, lat_col = cand_lon, cand_lat
                src_frame = "altaz"
                break

    if lon_col is None or lat_col is None or src_frame is None:
        raise ValueError("Source coordinates (RA/DEC, GLON/GLAT, or AZ/EL) not found in table.")

    lon = df[lon_col].to_numpy(dtype=float, copy=False)
    lat = df[lat_col].to_numpy(dtype=float, copy=False)

    # 2. 座標変換 (Astropy内部で最適化済み)
    frame_map = {"j2000": "fk5", "b1950": "fk4", "lb": "galactic"}
    tgt_frame = frame_map.get(frame.lower(), frame.lower())

    target_coords: Optional[SkyCoord] = None
    ref_sky: Optional[SkyCoord] = None

    if src_frame in ("icrs", "galactic"):
        coords = SkyCoord(lon * u.deg, lat * u.deg, frame=src_frame)
        try:
            target_coords = coords.transform_to(tgt_frame)
        except Exception as e:
            raise ValueError(f"Coordinate transformation to '{tgt_frame}' failed: {e}")

        lon = target_coords.data.lon.deg
        lat = target_coords.data.lat.deg
    else:
        # altaz: ここでは差分計算のみを目的として、変換は行わない（site/timeが無いと不可）
        if tgt_frame not in ("altaz",):
            raise ValueError(
                "AltAz -> other frame transformation requires site/time metadata. "
                "Please provide RA/DEC or GLON/GLAT, or enrich metadata (EarthLocation + obstime)."
            )

    # 3. リファレンス設定
    if ref_coord is None:
        # Header (OBSRA/DEC) または データ平均
        if src_frame == "icrs" and ("OBSRA" in sc.meta and "OBSDEC" in sc.meta):
            try:
                ref_sky = SkyCoord(sc.meta["OBSRA"] * u.deg, sc.meta["OBSDEC"] * u.deg, frame="icrs").transform_to(tgt_frame)
                lon0 = float(ref_sky.data.lon.deg)
                lat0 = float(ref_sky.data.lat.deg)
            except Exception:
                ref_sky = None
                lon0, lat0 = _circular_mean_deg(lon), float(np.nanmean(lat))
        else:
            lon0, lat0 = _circular_mean_deg(lon), float(np.nanmean(lat))
    else:
        if isinstance(ref_coord, SkyCoord):
            ref_sky = ref_coord
        else:
            # 文字列はSkyCoordに委ねる（例: "05h35m14.5s -05d22m30s" 等）
            ref_sky = SkyCoord(ref_coord)

        try:
            ref_sky = ref_sky.transform_to(tgt_frame)
            lon0 = float(ref_sky.data.lon.deg)
            lat0 = float(ref_sky.data.lat.deg)
        except Exception as e:
            raise ValueError(f"ref_coord transformation to '{tgt_frame}' failed: {e}")

    # 4. オフセット計算
    # 周期境界処理 (-180 ~ 180)
    dlon = (lon - lon0 + 180.0) % 360.0 - 180.0
    dlat = lat - lat0

    cos_mode_n = str(cos_mode).strip().lower()
    if cos_mode_n in ("point", "local", "each"):
        cos_lat = lat
    elif cos_mode_n in ("ref", "reference", "center"):
        cos_lat = lat0
    else:
        raise ValueError("cos_mode must be 'point' or 'ref'.")

    proj = projection.upper()
    if proj in ("GLS", "SFL"):
        ofs_x = dlon * np.cos(np.deg2rad(cos_lat))
        ofs_y = dlat
    elif proj in ("CAR", "NONE"):
        ofs_x = dlon
        ofs_y = dlat
    else:
        raise ValueError(f"Unknown projection '{proj}'. Supported: GLS, SFL, CAR, SIN.")
        

    # 5. 単位変換
    unit_l = str(unit).lower()
    scale_map = {"arcsec": 3600.0, "arcmin": 60.0, "deg": 1.0}
    if unit_l not in scale_map:
        raise ValueError(f"Unknown unit '{unit}'. Choose from {list(scale_map)}.")
    scale = scale_map[unit_l]

    ofs_x = ofs_x * scale
    ofs_y = ofs_y * scale

    if verbose:
        if ref_sky is not None:
            ref_str = ref_sky.to_string("decimal")
        else:
            ref_str = f"lon0={lon0}, lat0={lat0}"
        print(f"Offsets calculated (Ref={ref_str}, Frame={tgt_frame}, Proj={projection}, cos_mode={cos_mode_n})")

    return pd.DataFrame({"OFS_LON": ofs_x, "OFS_LAT": ofs_y}, index=df.index)


def _circular_mean_deg(angle_deg: np.ndarray) -> float:
    """0/360 の境界を跨ぐ場合でも破綻しにくい circular mean (deg)."""
    ang = np.deg2rad(np.asarray(angle_deg, dtype=float))
    # NaN を除外
    m = np.isfinite(ang)
    if not np.any(m):
        return float("nan")
    ang = ang[m]
    s = np.nanmean(np.sin(ang))
    c = np.nanmean(np.cos(ang))
    return float(np.rad2deg(np.arctan2(s, c)) % 360.0)


# =========================================================
# 3. Merge (結合)
# =========================================================

def merge_scantables(
    inputs: Sequence[Union[Scantable, str]],
    sort_by_time: bool = False,
    shift_scan_id: bool = True  # <--- 追加: デフォルトでSCAN IDを自動シフトする
) -> Scantable:
    """
    複数のScantable（またはファイルパス）を高速に結合する。
    VLAデータと固定長データの混在も自動処理する。

    [FIX] FITS由来の Big-endian table が混じると iloc/sort で落ちるため、
    table は native endian に揃える。
    [NEW] shift_scan_id=True の場合、結合時に SCAN ID の重複を防ぐために自動でオフセットを加算する。
    """
    sc_list: List[Scantable] = []

    # ロード処理
    for inp in inputs:
        if isinstance(inp, str):
            sc_list.append(read_scantable(inp))
        else:
            sc_list.append(inp)

    if not sc_list:
        raise ValueError("No inputs to merge.")

    # ベースメタデータ
    base_meta = sc_list[0].meta.copy()

    # テーブル結合とSCAN IDの自動シフト処理
    tables = []
    scan_offset = 0
    
    for s in sc_list:
        # 元のテーブルを汚さないように浅いコピーを作成
        df = _df_to_native_endian(s.table).copy()
        
        if shift_scan_id and "SCAN" in df.columns:
            max_scan = df["SCAN"].max()
            if pd.notna(max_scan):
                df["SCAN"] = df["SCAN"] + scan_offset
                # 次のテーブルのためにオフセットを更新
                scan_offset += int(max_scan) + 1
                
        tables.append(df)

    merged_table = pd.concat(tables, axis=0, ignore_index=True)

    # データ結合
    all_data: List[np.ndarray] = []
    for s in sc_list:
        d = s.data
        if isinstance(d, list):
            all_data.extend(d)
        elif isinstance(d, np.ndarray):
            all_data.extend(list(d))
        else:
            raise TypeError(f"Unknown data type: {type(d)}")

    # ソート処理
    if sort_by_time:
        ts = _resolve_table_timestamps(merged_table)
        if ts is not None and not ts.isna().all():
            sort_vals = ts.view("i8")
            sort_idx = np.argsort(sort_vals)
            merged_table = merged_table.iloc[sort_idx].reset_index(drop=True)
            all_data = [all_data[i] for i in sort_idx]

    merged_hist = {"merged_count": len(inputs), "shifted_scan_id": shift_scan_id}

    return Scantable(
        meta=base_meta,
        data=all_data,
        table=merged_table,
        history=merged_hist
    )

# =========================================================
# 4. Modification (変更)
# =========================================================

def update_metadata(
    sc: Scantable,
    column: str,
    value: Any,
    rows: Union[str, List[int], None] = None,
    force: bool = False,
    verbose: bool = True
) -> None:
    """
    指定した行(ID)のカラム値を安全に変更する。
    ScantableをIn-placeで更新する。

    エイリアス対応:
    "RESTFREQ" を更新した場合、テーブル内に "RESTFRQ" が存在すればそれも同時に更新する（逆も同様）。
    """
    df = sc.table
    col_upper = column.upper()

    # 安全性チェック
    if not force:
        if col_upper in DANGER_COLS:
            if verbose:
                print(f"[Block] Modification of critical column '{column}' is protected. Use force=True.")
            return

        if col_upper in ENUM_VALS:
            valid = np.array(list(ENUM_VALS[col_upper]))
            vals = np.ravel(value)
            if not np.all(np.isin(vals, valid)):
                if verbose:
                    print(f"[Block] Invalid values for {column}. Allowed: {valid}")
                return

    # ロジックチェック (TEMPSCAL vs BEAMEFF)
    if col_upper == "TEMPSCAL":
        vals = np.ravel(value)
        if np.any(vals == "TR*"):
            has_beam = ("BEAMEFF" in df and df["BEAMEFF"].notna().any()) or \
                       ("BEAMEFF" in sc.meta and float(sc.meta.get("BEAMEFF", 0)) > 0)
            if not has_beam:
                if not force:
                    print("[Warning] Setting TEMPSCAL='TR*' but valid BEAMEFF not found. Blocked.")
                    return
                elif verbose:
                    print("[Warning] Forcing TEMPSCAL='TR*' despite missing BEAMEFF.")

    idxs = _parse_row_selector(rows, len(df))
    if len(idxs) == 0:
        if verbose:
            print("No rows selected.")
        return

    targets = [column]
    if col_upper in COLUMN_ALIASES:
        aliases = COLUMN_ALIASES[col_upper]
        for alias in aliases:
            if alias in df.columns:
                targets.append(alias)

    for target_col in targets:
        if target_col not in df.columns:
            if verbose:
                print(f"Creating new column '{target_col}'")
            df[target_col] = np.nan

        try:
            df.loc[df.index[idxs], target_col] = value
            if verbose:
                extra_msg = " (Alias sync)" if target_col != column else ""
                print(f"Updated '{target_col}'{extra_msg} for {len(idxs)} rows.")
        except Exception as e:
            print(f"Update failed for '{target_col}': {e}")



def set_beameff(
    sc: Scantable,
    efficiency: Union[float, np.ndarray, Sequence[float], Mapping[Any, float]],
    rows: Union[str, List[int], int, slice, None] = None,
    *,
    key_columns: Union[str, Sequence[str]] = ("FDNUM", "IFNUM", "PLNUM"),
    strict: bool = False,
    verbose: bool = True,
) -> None:
    """Thin wrapper around ``tempscale.set_beameff()``."""
    from .tempscale import set_beameff as _set_beameff_impl
    return _set_beameff_impl(
        sc,
        efficiency=efficiency,
        rows=rows,
        key_columns=key_columns,
        strict=strict,
        verbose=verbose,
    )


def apply_relative_scale(
    sc: Scantable,
    scale: Union[float, np.ndarray, Sequence[float], Mapping[Any, float]],
    rows: Union[str, List[int], int, slice, None] = None,
    *,
    key_columns: Union[str, Sequence[str]] = ("FDNUM", "IFNUM", "PLNUM"),
    strict: bool = False,
    invalidate_derived: bool = True,
    verbose: bool = True,
) -> None:
    """Thin wrapper around ``tempscale.apply_relative_scale()``."""
    from .tempscale import apply_relative_scale as _impl
    return _impl(
        sc,
        scale=scale,
        rows=rows,
        key_columns=key_columns,
        strict=strict,
        invalidate_derived=invalidate_derived,
        verbose=verbose,
    )


def apply_global_scale(
    sc: Scantable,
    factor: float,
    rows: Union[str, List[int], int, slice, None] = None,
    *,
    invalidate_derived: bool = True,
    verbose: bool = True,
) -> None:
    """Thin wrapper around ``tempscale.apply_global_scale()``."""
    from .tempscale import apply_global_scale as _impl
    return _impl(
        sc,
        factor=factor,
        rows=rows,
        invalidate_derived=invalidate_derived,
        verbose=verbose,
    )

# =========================================================
# 5. Filter / Query (検索・抽出)
# =========================================================

def find_scans(
    sc: Scantable,
    query: Optional[str] = None,
    extra_data: Optional[pd.DataFrame] = None,
    **kwargs
) -> np.ndarray:
    """
    条件に合う行ID (0-based positional index) を返す。
    extra_dataを指定すると、そのカラムも含めて検索可能。

    [FIX] Big-endian table だと query/isin の内部処理で落ちることがあるため、native endian に揃えてから検索する。
    """
    df = _df_to_native_endian(sc.table)

    if extra_data is not None:
        if len(extra_data) == len(df):
            search_df = pd.concat([df.reset_index(drop=True), extra_data.reset_index(drop=True)], axis=1)
        else:
            warnings.warn("extra_data size mismatch. Ignoring in query.")
            search_df = df
    else:
        search_df = df

    mask = np.ones(len(search_df), dtype=bool)

    if query:
        try:
            subset = search_df.query(query)
            mask &= search_df.index.isin(subset.index)
        except Exception as e:
            raise ValueError(f"Invalid query string '{query}': {e}")

    for k, v in kwargs.items():
        col = k if k in search_df else k.upper()
        if col in search_df:
            if isinstance(v, (list, tuple, np.ndarray)):
                mask &= search_df[col].isin(v)
            else:
                mask &= (search_df[col] == v)
        else:
            available = list(search_df.columns)
            raise KeyError(f"Column '{k}' not found. Available: {available}")

    return np.where(mask)[0]


def filter_scantable(
    sc: Scantable,
    query: Optional[str] = None,
    extra_data: Optional[pd.DataFrame] = None,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    **kwargs
) -> Scantable:
    """
    条件抽出した新しいScantableを返す。
    Scantable自体は不変、メタデータはコピーされる。

    [NEW] rows を追加:
        - rows=None のときは find_scans の結果をそのまま採用
        - rows が指定された場合は、その行集合と find_scans の結果の交差を取る
          (飛び飛びの行番号選択に便利)
    """
    idxs_found = find_scans(sc, query, extra_data=extra_data, **kwargs)
    if rows is None:
        idxs = idxs_found
    else:
        idxs_rows = _parse_row_selector(rows, len(sc.table))
        # 交差（順序は rows 側を優先）
        found_set = set(map(int, idxs_found))
        idxs = np.array([int(i) for i in idxs_rows if int(i) in found_set], dtype=int)

    if len(idxs) == 0:
        warnings.warn("Filter resulted in empty Scantable.")

    # テーブル抽出（endian を直してから iloc）
    tab = _df_to_native_endian(sc.table)
    new_table = tab.iloc[idxs].reset_index(drop=True)

    # データ抽出
    d = sc.data
    if isinstance(d, list):
        new_data = [d[int(i)] for i in idxs]
    elif isinstance(d, np.ndarray):
        new_data = d[idxs]
    else:
        new_data = []

    new_meta = sc.meta.copy()
    new_hist = sc.history.copy()
    new_hist["filter_query"] = str(query)
    if rows is not None:
        new_hist["filter_rows"] = str(rows)

    return Scantable(meta=new_meta, data=new_data, table=new_table, history=new_hist)


def select_rows(sc: Scantable, rows: Union[str, slice, Sequence[int], int]) -> Scantable:
    """
    行番号(0-based positional)で Scantable を抽出するショートカット。
    filter_scantable(rows=...) と同等だが、クエリ無しで使う用途向け。
    """
    idxs = _parse_row_selector(rows, len(sc.table))
    if len(idxs) == 0:
        warnings.warn("select_rows resulted in empty Scantable.")

    tab = _df_to_native_endian(sc.table)
    new_table = tab.iloc[idxs].reset_index(drop=True)

    d = sc.data
    if isinstance(d, list):
        new_data = [d[int(i)] for i in idxs]
    elif isinstance(d, np.ndarray):
        new_data = d[idxs]
    else:
        new_data = []

    new_meta = sc.meta.copy()
    new_hist = sc.history.copy()
    new_hist["select_rows"] = str(rows)

    return Scantable(meta=new_meta, data=new_data, table=new_table, history=new_hist)


# =========================================================
# Internal Helper
# =========================================================

def _parse_row_selector(rows: Union[str, slice, Sequence[int], int, None], max_len: int) -> np.ndarray:
    """行指定をパースしてインデックス配列(0-based positional)を返す"""
    if rows is None:
        return np.arange(max_len)

    if isinstance(rows, int):
        return np.array([rows]) if 0 <= rows < max_len else np.array([], dtype=int)

    if isinstance(rows, slice):
        return np.arange(max_len)[rows]

    if isinstance(rows, (list, tuple, np.ndarray)):
        arr = np.array(rows, dtype=int)
        valid = arr[(arr >= 0) & (arr < max_len)]
        # 順序を保ちたい場合があるので unique はしない（重複は許可）
        return valid

    if isinstance(rows, str):
        indices = []
        base_indices = np.arange(max_len)
        for p in rows.split(","):
            p = p.strip()
            if not p:
                continue
            if ":" in p:
                try:
                    slice_parts = [int(x) if x else None for x in p.split(":")]
                    s_obj = slice(*slice_parts)
                    indices.append(base_indices[s_obj])
                except ValueError:
                    print(f"Warning: Invalid slice syntax '{p}'")
            else:
                try:
                    idx = int(p)
                    if 0 <= idx < max_len:
                        indices.append(np.array([idx], dtype=int))
                except ValueError:
                    print(f"Warning: Invalid index '{p}'")

        if not indices:
            return np.array([], dtype=int)

        return np.concatenate(indices)

    raise TypeError(f"Invalid type for rows selector: {type(rows)}")
