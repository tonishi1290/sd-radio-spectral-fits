# src/sd_radio_spectral_fits/map/ps_gridder.py
from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Optional, Literal, Tuple
from astropy.io import fits

from ..regrid_vlsrk import Standardizer
from .wcs_proj import project_to_plane
from .config import GridInput, GridResult
from ..tempscale import beameff_array, tempscal_array, ta_to_tr
from ..scantable_utils import _df_to_native_endian
from .gridder import _extract_rms, _extract_tsys_scalar, _extract_meta_col, _resolve_time_array
from .fits_io import _fill_header_metadata, _add_diagnostic_hdus

@dataclass
class PSMapConfig:
    """Position Switch観測用のグリッディング設定クラス"""
    coord_sys: str = "icrs"
    projection: str = "SFL"
    ref_lon: Optional[float] = None
    ref_lat: Optional[float] = None
    
    # --- 1. 格子の定義 ---
    x_grid: float = 30.0  # X方向の格子間隔 [arcsec]
    y_grid: float = 30.0  # Y方向の格子間隔 [arcsec]
    grid_anchor_offsets: Tuple[float, float] = (0.0, 0.0)  # 格子の位相固定用アンカー [arcsec]
    
    # --- 2. 範囲の定義 ---
    grid_bounds_offsets: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None
    
    # --- 3. 割当・採用条件 ---
    grid_tol: float = 0.3          # 近傍格子点への許容率 (例: 0.3なら格子間隔の30%以内)
    invert_x: bool = True          # FITS表示規則に従い、左方向をX軸の正とするか
    combine: Literal["mean", "median", "sum"] = "mean"
    weight_mode: Literal["uniform", "rms", "tint"] = "rms"
    return_cell_members: bool = False
    
    # --- 4. ログ出力 ---
    verbose: bool = True           # 処理の詳細（重複点の情報など）を出力するか

def grid_ps(input_data: GridInput, config: PSMapConfig) -> GridResult:
    """PS観測用のグリッディング・コアエンジン"""
    if config.verbose:
        print("\n--- Starting PS Gridding ---")
        
    x = np.asarray(input_data.x, dtype=float)
    y = np.asarray(input_data.y, dtype=float)
    spec = np.asarray(input_data.spec, dtype=float)
    
    # フィルタリング (不正なデータや折り返しスキャンを除外)
    valid = np.asarray(input_data.flag) > 0 if input_data.flag is not None else np.ones(len(x), dtype=bool)
    if input_data.is_turnaround is not None:
        valid &= ~np.asarray(input_data.is_turnaround, dtype=bool)
    valid &= np.isfinite(x) & np.isfinite(y)
    valid &= np.any(np.isfinite(spec), axis=1)
    
    initial_count = len(x)
    valid_count = np.count_nonzero(valid)
    
    if config.verbose:
        print(f"Data filtering: {valid_count} / {initial_count} spectra remain valid.")
    
    if valid_count == 0:
        raise ValueError("No valid spectra remain after filtering.")
        
    x_v, y_v, spec_v = x[valid], y[valid], spec[valid]
    
    rms = np.asarray(input_data.rms)[valid] if input_data.rms is not None else np.ones(len(x_v))
    tint = np.asarray(input_data.tint)[valid] if input_data.tint is not None else np.ones(len(x_v))
    time_arr = np.asarray(input_data.time)[valid] if input_data.time is not None else np.zeros(len(x_v))
    tsys_arr = np.asarray(input_data.tsys)[valid] if input_data.tsys is not None else np.full(len(x_v), np.nan)
    
    # 重みの計算
    w_v = np.ones(len(x_v), dtype=float)
    if config.weight_mode == "rms":
        w_v = 1.0 / (rms**2 + 1e-12)
    elif config.weight_mode == "tint":
        w_v = tint
        
    # --- 範囲の決定 ---
    x_min, x_max = np.min(x_v), np.max(x_v)
    y_min, y_max = np.min(y_v), np.max(y_v)
    
    if config.grid_bounds_offsets is not None:
        (bx_min, by_min), (bx_max, by_max) = config.grid_bounds_offsets
        x_min, x_max = min(bx_min, bx_max), max(bx_min, bx_max)
        y_min, y_max = min(by_min, by_max), max(by_min, by_max)
        
    x_anchor, y_anchor = config.grid_anchor_offsets
    xg, yg = config.x_grid, config.y_grid
    
    # アンカーを基準としたインデックス範囲
    eps = 1e-4
    ix_min = int(np.floor((x_min - x_anchor) / xg + eps))
    ix_max = int(np.ceil((x_max - x_anchor) / xg - eps))
    iy_min = int(np.floor((y_min - y_anchor) / yg + eps))
    iy_max = int(np.ceil((y_max - y_anchor) / yg - eps))
    
    nx = ix_max - ix_min + 1
    ny = iy_max - iy_min + 1
    
    if config.verbose:
        print(f"Grid Size: {nx} (X) x {ny} (Y) pixels")
        print(f"Grid Anchor: X={x_anchor:.1f}, Y={y_anchor:.1f} arcsec")
        print(f"Grid Cell Size: {xg:.1f} x {yg:.1f} arcsec")
    
    x_centers_asc = x_anchor + np.arange(ix_min, ix_max + 1) * xg
    y_centers_asc = y_anchor + np.arange(iy_min, iy_max + 1) * yg
    
    # FITS規則への準拠 (invert_x=Trueなら、ピクセルix=0が左端(最大X)になるよう反転)
    if config.invert_x:
        x_centers = x_centers_asc[::-1]
        cdelt1_arcsec = -xg
    else:
        x_centers = x_centers_asc
        cdelt1_arcsec = xg
        
    y_centers = y_centers_asc
    cdelt2_arcsec = yg
    
    # --- セル割り当て ---
    idx_x_asc = np.round((x_v - x_anchor) / xg).astype(int) - ix_min
    idx_y_asc = np.round((y_v - y_anchor) / yg).astype(int) - iy_min
    
    dx = np.abs(x_v - (x_anchor + (idx_x_asc + ix_min) * xg))
    dy = np.abs(y_v - (y_anchor + (idx_y_asc + iy_min) * yg))
    
    # Tolチェックと境界チェック
    valid_tol = (dx <= xg * config.grid_tol) & (dy <= yg * config.grid_tol)
    valid_bounds = (idx_x_asc >= 0) & (idx_x_asc < nx) & (idx_y_asc >= 0) & (idx_y_asc < ny)
    valid_all = valid_tol & valid_bounds
    
    x_f, y_f, spec_f, w_f = x_v[valid_all], y_v[valid_all], spec_v[valid_all], w_v[valid_all]
    rms_f, tint_f, time_f, tsys_f = rms[valid_all], tint[valid_all], time_arr[valid_all], tsys_arr[valid_all]
    
    idx_x_f_asc, idx_y_f = idx_x_asc[valid_all], idx_y_asc[valid_all]
    idx_x_f = (nx - 1) - idx_x_f_asc if config.invert_x else idx_x_f_asc
        
    flat_idx = idx_y_f * nx + idx_x_f
    
    dropped_tol = valid_count - np.count_nonzero(valid_all)
    if config.verbose:
        print(f"Cell Assignment: {np.count_nonzero(valid_all)} spectra mapped.")
        if dropped_tol > 0:
            print(f"  -> Dropped {dropped_tol} spectra (exceeded grid_tol={config.grid_tol} or out of bounds).")
    
    # --- 合成 (Combine) ---
    nchan = spec.shape[1]
    cube_flat = np.full((ny * nx, nchan), np.nan, dtype=np.float32)
    weight_flat = np.zeros(ny * nx, dtype=np.float32)
    hit_flat = np.zeros(ny * nx, dtype=np.int32)
    mask_flat = np.zeros(ny * nx, dtype=bool)
    
    tint_flat = np.zeros(ny * nx, dtype=np.float32)
    time_flat = np.full(ny * nx, np.nan, dtype=np.float32)
    rms_flat = np.full(ny * nx, np.nan, dtype=np.float32)
    tsys_flat = np.full(ny * nx, np.nan, dtype=np.float32)
    
    # 高速グループ化
    order = np.argsort(flat_idx)
    unique_idx, start_indices, counts = np.unique(flat_idx[order], return_index=True, return_counts=True)
    
    if config.verbose:
        active_cells = len(unique_idx)
        total_cells = nx * ny
        fill_factor = (active_cells / total_cells) * 100 if total_cells > 0 else 0
        max_hits = np.max(counts) if len(counts) > 0 else 0
        avg_hits = np.mean(counts) if len(counts) > 0 else 0
        multi_hit_cells = np.count_nonzero(counts > 1)
        
        print("\n--- Cell Duplication & Combine Stats ---")
        print(f"Active cells : {active_cells} / {total_cells} ({fill_factor:.1f}% filled)")
        print(f"Combine mode : '{config.combine}' with weight_mode='{config.weight_mode}'")
        print(f"Max hits/cell: {max_hits}")
        print(f"Avg hits/cell: {avg_hits:.2f}")
        print(f"Cells with multiple spectra: {multi_hit_cells}")
    
    split_spec = np.split(spec_f[order], start_indices[1:])
    split_w = np.split(w_f[order], start_indices[1:])
    split_rms = np.split(rms_f[order], start_indices[1:])
    split_tint = np.split(tint_f[order], start_indices[1:])
    split_time = np.split(time_f[order], start_indices[1:])
    split_tsys = np.split(tsys_f[order], start_indices[1:])
    
    for i, f_idx in enumerate(unique_idx):
        sp = split_spec[i]
        wv = split_w[i]
        v_sp = np.isfinite(sp)
        wv_2d = wv[:, None] * v_sp
        W_sum_chan = np.sum(wv_2d, axis=0)
        W_sum = np.sum(wv)
        
        hit_flat[f_idx] = len(sp)
        mask_flat[f_idx] = True
        
        if config.combine == "mean":
            val = np.sum(np.where(v_sp, sp * wv_2d, 0), axis=0)
            m_chan = W_sum_chan > 0
            cube_view = cube_flat[f_idx]
            cube_view[m_chan] = val[m_chan] / W_sum_chan[m_chan]
            weight_flat[f_idx] = W_sum
        elif config.combine == "median":
            cube_flat[f_idx] = np.nanmedian(sp, axis=0)
            weight_flat[f_idx] = W_sum
        elif config.combine == "sum":
            cube_flat[f_idx] = np.nansum(sp, axis=0)
            weight_flat[f_idx] = W_sum
            
        time_flat[f_idx] = np.sum(split_time[i] * wv) / W_sum if W_sum > 0 else np.nan
        tint_flat[f_idx] = np.sum(split_tint[i])  # 積分時間は合算
        rms_flat[f_idx] = np.sqrt(np.sum((split_rms[i]**2) * wv) / W_sum) if W_sum > 0 else np.nan
        tsys_flat[f_idx] = np.nanmedian(split_tsys[i])
        
    res = GridResult(
        cube=cube_flat.reshape(ny, nx, nchan),
        weight_map=weight_flat.reshape(ny, nx),
        hit_map=hit_flat.reshape(ny, nx),
        mask_map=mask_flat.reshape(ny, nx),
        tint_map=tint_flat.reshape(ny, nx),
        time_map=time_flat.reshape(ny, nx),
        rms_map=rms_flat.reshape(ny, nx),
        tsys_map=tsys_flat.reshape(ny, nx),
    )
    
    # WCS用メタデータの保存
    res.meta = {
        "x0_arcsec": x_centers[0],
        "y0_arcsec": y_centers[0],
        "cdelt1_arcsec": cdelt1_arcsec,
        "cdelt2_arcsec": cdelt2_arcsec,
        "nx": nx,
        "ny": ny,
        "coord_sys": config.coord_sys,
        "projection": config.projection
    }
    
    if config.verbose:
        print("--- PS Gridding Complete ---\n")
        
    return res

def build_ps_wcs_dict(
    coord_sys: str, projection: str, lon0: float, lat0: float, 
    x0_arcsec: float, y0_arcsec: float, cdelt1_arcsec: float, cdelt2_arcsec: float, 
    nx: int, ny: int
) -> dict:
    """PS用の空間WCSヘッダ構築 (x, y個別のCDELTに対応)"""
    sys_upper = coord_sys.upper()
    proj_upper = projection.upper()
    
    proj_fits = "SFL" if proj_upper in ("GLS", "SFL", "SINE") else "CAR"
    if sys_upper in ("RADEC", "ICRS"):
        ctype1, ctype2 = f"RA---{proj_fits}", f"DEC--{proj_fits}"
    elif sys_upper in ("GALACTIC", "GAL"):
        ctype1, ctype2 = f"GLON-{proj_fits}", f"GLAT-{proj_fits}"
    else:
        ctype1, ctype2 = f"LON--{proj_fits}", f"LAT--{proj_fits}"

    cdelt1, cdelt2 = cdelt1_arcsec / 3600.0, cdelt2_arcsec / 3600.0
    crpix1 = 1.0 - (x0_arcsec / cdelt1_arcsec)
    crpix2 = 1.0 - (y0_arcsec / cdelt2_arcsec)

    if proj_fits == "SFL":
        crval1, crval2 = float(lon0), 0.0
        crpix2 = float(crpix2 - (lat0 / cdelt2))
    else:
        crval1, crval2 = float(lon0), float(lat0)

    return {
        "CTYPE1": ctype1, "CRVAL1": crval1, "CDELT1": cdelt1, "CRPIX1": crpix1, "CUNIT1": "deg",
        "CTYPE2": ctype2, "CRVAL2": crval2, "CDELT2": cdelt2, "CRPIX2": crpix2, "CUNIT2": "deg",
    }


def _add_ps_diagnostic_hdus(hdul, grid_res: GridResult, base_header: fits.Header):
    """PS観測専用の厳選された2D診断マップだけをFITSに追加する"""
    header_2d = base_header.copy()
    # 2Dマップには速度軸(3軸目)の情報は不要なので削除
    for k in ["CTYPE3", "CUNIT3", "CRPIX3", "CRVAL3", "CDELT3"]:
        if k in header_2d: 
            del header_2d[k]

    # PS観測でエンドユーザーが本当に必要とするマップのみに限定
    maps = [
        ("HIT", grid_res.hit_map, "HitCount", ""),                 # 観測点のカバー状況と積分回数の確認
        ("RMS", grid_res.rms_map, "BaselineRMS", "K"),             # ノイズレベルの確認
        ("TSYS", grid_res.tsys_map, "SystemTemp", "K"),            # 観測中のシステム雑音温度の確認
        ("TINT", grid_res.tint_map, "IntegrationTime", "s"),       # 実効的な総積分時間の確認
    ]

    for name, data, btype, bunit in maps:
        if data is not None:
            h = header_2d.copy()
            h["BTYPE"], h["BUNIT"] = btype, bunit
            hdul.append(fits.ImageHDU(data=data.astype(np.float32), header=h, name=name))
            


def save_ps_map_fits(
    grid_res: GridResult, v_tgt: np.ndarray, lon0: float, lat0: float,
    out_path: str, out_scale: str = "TA*", rep_beameff: float = 1.0,
):
    """PS専用のFITS出力関数"""
    ny, nx, nchan = grid_res.cube.shape
    m = grid_res.meta
    
    wcs_hdr = build_ps_wcs_dict(
        m["coord_sys"], m["projection"], lon0, lat0, 
        m["x0_arcsec"], m["y0_arcsec"], m["cdelt1_arcsec"], m["cdelt2_arcsec"], nx, ny
    )
    
    header = fits.Header()
    header.update(wcs_hdr)
    _fill_header_metadata(header, m["coord_sys"], lon0, lat0, out_scale, rep_beameff, grid_res)

    header["CTYPE3"] = "VRAD"
    header["CUNIT3"] = "km/s"
    header["CRPIX3"] = 1.0
    header["CRVAL3"] = float(v_tgt[0])
    header["CDELT3"] = float(v_tgt[1] - v_tgt[0]) if len(v_tgt) > 1 else 1.0
    header["SPECSYS"] = "LSRK"
    header["VELDEF"] = ("RADI-LSR", "Radio velocity definition")

    cube_fits = np.transpose(grid_res.cube, (2, 0, 1))
    hdu_primary = fits.PrimaryHDU(data=cube_fits.astype(np.float32), header=header)
    hdul = fits.HDUList([hdu_primary])

    # OTF用の汎用関数ではなく、PS専用の厳選関数を呼ぶ
    _add_ps_diagnostic_hdus(hdul, grid_res, header)
    hdul.writeto(out_path, overwrite=True)
    
    

def run_ps_mapping_pipeline(
    scantable, config: PSMapConfig, output_fits: str,
    out_scale: str = "TA*", dv_kms: float = None,
    vmin_kms: float = None,  # <--- 追加
    vmax_kms: float = None,  # <--- 追加
):
    """PS観測用のデータ抽出からFITS出力までを行う統合パイプライン"""
    if config.verbose:
        print("1. Regridding velocity axis (Standardizer)...")
        if vmin_kms is not None and vmax_kms is not None:
            print(f"   -> Requested Velocity Range: {vmin_kms} to {vmax_kms} km/s")

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


    std = Standardizer(scantable, v_corr_col="VFRAME")

    full_matrix, v_tgt = std.get_matrix(dv=dv_kms, vmin=vmin_kms, vmax=vmax_kms)

    if config.verbose:
        print("2. Projecting coordinates to local plane...")
        
    table = _df_to_native_endian(scantable.table).copy()
    is_galactic = config.coord_sys.lower() in ("galactic", "gal")
    lon_deg = pd.to_numeric(table["GLON" if is_galactic else "RA"], errors="coerce").to_numpy(float)
    lat_deg = pd.to_numeric(table["GLAT" if is_galactic else "DEC"], errors="coerce").to_numpy(float)

    lon0 = float(config.ref_lon) if config.ref_lon is not None else float(np.nanmedian(lon_deg))
    lat0 = float(config.ref_lat) if config.ref_lat is not None else float(np.nanmedian(lat_deg))
    x_arcsec, y_arcsec = project_to_plane(lon_deg, lat_deg, lon0, lat0, projection=config.projection)

    if config.verbose:
        print("3. Handling Temperature Scale...")
        
    out_scale_norm = str(out_scale).strip().upper()
    n_rows = len(table)
    beameff_vec = beameff_array(table, scantable.meta, n_rows)
    tempscal_vec = tempscal_array(table, scantable.meta, n_rows)

    rep_beameff = np.nan
    if out_scale_norm == "TR*":
        mask_ta = (tempscal_vec == "TA*")
        if np.any(mask_ta):
            full_matrix[mask_ta] = ta_to_tr(full_matrix[mask_ta], beameff_vec[mask_ta, None])
            if config.verbose:
                print(f" -> Converted {np.count_nonzero(mask_ta)} spectra to TR*.")
        rep_beameff = 1.0
    else:
        valid_beameff = beameff_vec[np.isfinite(beameff_vec)]
        rep_beameff = float(np.median(valid_beameff)) if len(valid_beameff) > 0 else np.nan

    grid_input = GridInput(
        x=x_arcsec, y=y_arcsec, spec=full_matrix,
        flag=(table["OBSMODE"].astype(str).str.upper() == "ON").to_numpy(dtype=bool) if "OBSMODE" in table.columns else np.ones(n_rows, dtype=bool),
        time=_resolve_time_array(table),
        rms=_extract_rms(table, out_scale_norm, beameff_vec),
        tint=_extract_meta_col(table, ("EXPOSURE", "INTTIME", "DUR")),
        tsys=_extract_tsys_scalar(table),
    )
    
    if config.verbose:
        print("4. Executing PS Gridding Engine...")
        
    res = grid_ps(grid_input, config)
    res.meta["RESTFREQ"] = float(scantable.meta.get("RESTFREQ", scantable.meta.get("RESTFRQ", 0.0)))

    if config.verbose:
        print("5. Writing PS Multi-Extension FITS...")
        
    save_ps_map_fits(res, v_tgt, lon0, lat0, output_fits, out_scale_norm, rep_beameff)
    
    if config.verbose:
        print(f"Done! Saved to {output_fits}")
        
    return res
