# src/sd_radio_spectral_fits/map/gridder.py
import numpy as np
import pandas as pd
from ..regrid_vlsrk import Standardizer
from .core import grid_otf
from .wcs_proj import project_to_plane
from .config import MapConfig, GridInput
from .fits_io import save_map_fits
from ..tempscale import beameff_array, tempscal_array, ta_to_tr
from ..scantable_utils import calc_mapping_offsets, _df_to_native_endian


def _time_series_to_mjd_utc(series: pd.Series) -> np.ndarray:
    """Convert a timestamp-like Series to MJD UTC days (float)."""
    ts = pd.to_datetime(series, errors="coerce", utc=True)
    if len(ts) == 0:
        return np.array([], dtype=float)
    return ts.view("int64").astype(np.float64) / 86400e9 + 40587.0


def _resolve_time_array(table: pd.DataFrame) -> np.ndarray:
    """Resolve per-row time to MJD UTC days from preferred columns."""
    # Work on native-endian, pandas-safe table before extracting numeric arrays.
    if "MJD" in table.columns:
        mjd = pd.to_numeric(table["MJD"], errors="coerce").to_numpy(float)
        if np.isfinite(mjd).any():
            return mjd

    if "TIMESTAMP" in table.columns:
        mjd = _time_series_to_mjd_utc(table["TIMESTAMP"])
        if np.isfinite(mjd).any():
            return mjd

    date_col = None
    if "DATE-OBS" in table.columns:
        date_col = "DATE-OBS"
    elif "DATEOBS" in table.columns:
        date_col = "DATEOBS"

    if date_col is not None:
        base = _time_series_to_mjd_utc(table[date_col])
        if "TIME" in table.columns:
            sec = pd.to_numeric(table["TIME"], errors="coerce").to_numpy(float)
            if len(sec) == len(base):
                base = base + np.where(np.isfinite(sec), sec / 86400.0, np.nan)
        if np.isfinite(base).any():
            return base

    # Legacy fallback: TIME stored directly as MJD day.
    if "TIME" in table.columns:
        t = pd.to_numeric(table["TIME"], errors="coerce").to_numpy(float)
        if np.isfinite(t).any():
            finite = t[np.isfinite(t)]
            if finite.size and np.nanmedian(finite) > 1e4:
                return t

    return np.arange(len(table), dtype=float)

def run_mapping_pipeline(
    scantable,
    config: MapConfig,
    output_fits: str,
    coord_sys: str = "icrs",
    projection: str = "SFL",
    out_scale: str = "TA*",
    dv_kms: float = None,
    ref_lon: float = None,
    ref_lat: float = None,
):
    """
    Scantable から 3D FITS キューブを生成する汎用統合パイプライン。
    旧 run_otf_pipeline の全機能を包含し、PS観測にも対応可能。
    """

    # 0. SIG_* を除去（coadd後にVRADへ確定している場合は残骸になる）
    t = scantable.table
    if hasattr(t, "columns"):
        meta_ctype1 = str(getattr(scantable, "meta", {}).get("CTYPE1", "")).upper()
        if ("VEL" in meta_ctype1) or ("VRAD" in meta_ctype1):
            sig_cols = [c for c in t.columns if str(c).startswith("SIG_")]
            if sig_cols:
                if getattr(config, "verbose", False):
                    print(f"[info] drop SIG_* columns: {sig_cols}")
                scantable.table = t.drop(columns=sig_cols)

    # 1. 速度軸標準化 (Standardizer)
    print("1. Regridding velocity axis (Standardizer)...")
    std = Standardizer(scantable, v_corr_col="VFRAME")
    full_matrix, v_tgt = std.get_matrix(dv=dv_kms)

    import numpy as np
    print("v_tgt[0:3] =", v_tgt[:3])
    print("v_tgt[-3:] =", v_tgt[-3:])
    dv = np.nanmedian(np.diff(v_tgt))
    print("dv(median) =", dv, " range =", (np.nanmin(v_tgt), np.nanmax(v_tgt)))
    print("nonlinear check (max|d-dmed|) =", np.nanmax(np.abs(np.diff(v_tgt) - dv)))


    # 2. 座標投影 (deg -> arcsec)
    print("2. Projecting coordinates to local plane...")
    table = _df_to_native_endian(scantable.table).copy()
    is_galactic = coord_sys.lower() in ("galactic", "gal")
    lon_deg = pd.to_numeric(table["GLON" if is_galactic else "RA"], errors="coerce").to_numpy(float)
    lat_deg = pd.to_numeric(table["GLAT" if is_galactic else "DEC"], errors="coerce").to_numpy(float)

    valid_mask = np.isfinite(lon_deg) & np.isfinite(lat_deg)
    if not np.any(valid_mask):
        raise ValueError("No finite coordinates found.")

    lon0 = float(ref_lon) if ref_lon is not None else float(np.nanmedian(lon_deg[valid_mask]))
    lat0 = float(ref_lat) if ref_lat is not None else float(np.nanmedian(lat_deg[valid_mask]))
    
    x_arcsec, y_arcsec = project_to_plane(lon_deg, lat_deg, lon0, lat0, projection=projection)

    # 3. 温度スケールと BEAMEFF の処理
    print("3. Handling BEAMEFF and Temperature Scale...")
    out_scale_norm = str(out_scale).strip().upper()
    n_rows = len(table)
    
    beameff_vec = beameff_array(table, scantable.meta, n_rows)
    tempscal_vec = tempscal_array(table, scantable.meta, n_rows)

    rep_beameff = np.nan
    if out_scale_norm == "TR*":
        mask_ta = (tempscal_vec == "TA*")
        if np.any(mask_ta):
            full_matrix[mask_ta] = ta_to_tr(full_matrix[mask_ta], beameff_vec[mask_ta, None])
            print(f" -> Converted {np.count_nonzero(mask_ta)} spectra to TR*.")
        rep_beameff = 1.0 # TR*時は効率1.0として扱う
    else:
        valid_beameff = beameff_vec[np.isfinite(beameff_vec)]
        rep_beameff = float(np.median(valid_beameff)) if len(valid_beameff) > 0 else np.nan

    # 4. GridInput の詳細な組み立て (診断マップソースの抽出)
    grid_input = GridInput(
        x=x_arcsec,
        y=y_arcsec,
        spec=full_matrix,
        # ここに .to_numpy(dtype=bool) を追加して型安全にする
        flag=(table["OBSMODE"].astype(str).str.upper() == "ON").to_numpy(dtype=bool),
        time=_resolve_time_array(table),
        rms=_extract_rms(table, out_scale_norm, beameff_vec),
        tint=_extract_meta_col(table, ("EXPOSURE", "INTTIME", "DUR")),
        tsys=_extract_tsys_scalar(table),
        scan_id=pd.to_numeric(table["SCAN"], errors="coerce").to_numpy() if "SCAN" in table.columns else None,
        is_turnaround=np.asarray(table["IS_TURN"]).astype(bool) if "IS_TURN" in table.columns else None
    )
    
    # 5. グリッディング実行
    print("4. Executing Gridding Engine...")
    res = grid_otf(grid_input, config)

    # RESTFREQ と SPECSYS をメタデータに引き継ぐ
    
    if not res.meta:
        res.meta = {}
    res.meta["RESTFREQ"] = float(scantable.meta.get("RESTFREQ", scantable.meta.get("RESTFRQ", 0.0)))
    
    # StandardizerがVFRAMEを適用して速度軸を揃えているため、出力は必ずLSRKになる
    res.meta["SPECSYS"] = "LSRK"
    
    # 6. FITS 書き出し
    print("5. Writing Multi-Extension FITS...")
    save_map_fits(
        res, v_tgt, coord_sys, projection, lon0, lat0, config, 
        output_fits, out_scale=out_scale_norm, rep_beameff=rep_beameff
    )
    print(f"Done! Saved to {output_fits}")
    return res

# --- 内部ヘルパー (pipeline.py の抽出ロジックを整理) ---

def _extract_rms(table, out_scale, beameff_vec):
    if "BSL_RMS" not in table.columns:
        return np.ones(len(table))
    rms = pd.to_numeric(table["BSL_RMS"], errors="coerce").to_numpy(float, copy=True)
    rms[~np.isfinite(rms) | (rms <= 0)] = 1.0
    # TR* 指定時は RMS も BEAMEFF でスケールさせる
    if out_scale == "TR*":
        rms /= np.where(beameff_vec > 0, beameff_vec, 1.0)
    return rms

def _extract_tsys_scalar(table):
    if "TSYS" not in table.columns: return None
    return np.array([np.nanmean(x) if isinstance(x, (np.ndarray, list)) else x for x in table["TSYS"]])

def _extract_meta_col(table, names):
    for n in names:
        if n in table.columns:
            return pd.to_numeric(table[n], errors="coerce").to_numpy(float)
    return None

def create_grid_input(scantable, ref_coord=None, frame="ICRS", projection="SFL"):
    """
    ScantableからGridInputを安全に生成する。
    ※ 将来的に src/map/gridder.py に移設することを推奨します。
    """
    table = _df_to_native_endian(scantable.table).copy()
    
    # 投影座標(arcsec)の計算
    offsets_df = calc_mapping_offsets(
        scantable, 
        ref_coord=ref_coord, 
        frame=frame, 
        projection=projection, 
        unit="arcsec", 
        verbose=False
    )
    
    x_arcsec = offsets_df["OFS_LON"].to_numpy()
    y_arcsec = offsets_df["OFS_LAT"].to_numpy()

    obsmode = table["OBSMODE"].astype(str).str.upper().to_numpy() if "OBSMODE" in table.columns else np.full(len(table), "ON")
    time_arr = _resolve_time_array(table)
    scan_id = pd.to_numeric(table["SCAN"], errors="coerce").to_numpy(float) if "SCAN" in table.columns else None
    is_turn = np.asarray(table["IS_TURN"]).astype(bool) if "IS_TURN" in table.columns else None

    return GridInput(
        x=x_arcsec, 
        y=y_arcsec, 
        spec=np.asarray(scantable.data, dtype=float),
        flag=(obsmode == "ON"), 
        time=time_arr, 
        scan_id=scan_id, 
        is_turnaround=is_turn
    )
