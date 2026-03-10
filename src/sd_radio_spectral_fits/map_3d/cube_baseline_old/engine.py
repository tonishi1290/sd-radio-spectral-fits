# src/sd_radio_spectral_fits/cube_baseline/engine.py
import numpy as np
import warnings
from typing import Optional, Tuple, Dict

def fit_cube_baseline(
    cube_original: np.ndarray,
    v_axis: np.ndarray,
    signal_mask_3d: Optional[np.ndarray] = None,
    manual_mask_1d: Optional[np.ndarray] = None,
    target_mask_2d: Optional[np.ndarray] = None,
    poly_order: int = 1
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    3Dキューブに対してベースラインフィッティングを行うコアエンジン。
    """
    if cube_original.ndim != 3:
        raise ValueError(f"cube_original must be 3D, got {cube_original.ndim}D.")
    
    ny, nx, nchan = cube_original.shape
    
    if len(v_axis) != nchan:
        raise ValueError(f"Length of v_axis ({len(v_axis)}) does not match nchan ({nchan}).")
    if signal_mask_3d is not None and signal_mask_3d.shape != cube_original.shape:
        raise ValueError(f"signal_mask_3d shape {signal_mask_3d.shape} does not match cube {cube_original.shape}.")
    if manual_mask_1d is not None and manual_mask_1d.shape != (nchan,):
        raise ValueError(f"manual_mask_1d shape {manual_mask_1d.shape} does not match nchan ({nchan}).")
    if target_mask_2d is not None and target_mask_2d.shape != (ny, nx):
        raise ValueError(f"target_mask_2d shape {target_mask_2d.shape} does not match spatial shape {(ny, nx)}.")

    # 速度軸の数値安定化 (中心化と規格化) - NaN対策を強化
    v0 = np.nanmedian(v_axis)
    dv = float(np.nanmax(np.abs(v_axis - v0)))
    if (not np.isfinite(dv)) or dv <= 0:
        dv = 1.0
    x_axis = (v_axis - v0) / dv
    
    cube_baseline = np.zeros(cube_original.shape, dtype=np.float32)
    
    rms_map = np.full((ny, nx), np.nan, dtype=np.float32)
    n_used_map = np.zeros((ny, nx), dtype=np.int32)
    fit_ok_map = np.zeros((ny, nx), dtype=bool)
    
    if target_mask_2d is None:
        target_mask_2d = np.ones((ny, nx), dtype=bool)
        
    y_indices, x_indices = np.where(target_mask_2d)
    
    for iy, ix in zip(y_indices, x_indices):
        spec = cube_original[iy, ix, :]
        if np.all(np.isnan(spec)):
            continue
            
        # 1. 自動マスク(3D)の反転
        if signal_mask_3d is not None:
            fit_mask = (~signal_mask_3d[iy, ix, :].astype(bool)) & np.isfinite(spec)
        else:
            fit_mask = np.isfinite(spec)
            
        # 2. 手動マスク(1D)の合成 (3Dで実体化させず、ここで落とす)
        if manual_mask_1d is not None:
            fit_mask &= (~manual_mask_1d)
            
        n_used = np.count_nonzero(fit_mask)
        if n_used <= poly_order + 1:
            continue
            
        # 規格化された x_axis を使ってフィッティング
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', np.RankWarning)
            coeffs = np.polyfit(x_axis[fit_mask], spec[fit_mask], poly_order)
            
        baseline = np.polyval(coeffs, x_axis).astype(np.float32)
        cube_baseline[iy, ix, :] = baseline
        
        # ロバストRMS (MAD) の計算
        residual = (spec - baseline)[fit_mask]
        if len(residual) > 0:
            med = np.nanmedian(residual)
            mad = np.nanmedian(np.abs(residual - med))
            rms_map[iy, ix] = mad * 1.4826
            
        n_used_map[iy, ix] = n_used
        fit_ok_map[iy, ix] = True
        
    stats = {
        'rms_map': rms_map,
        'n_used_map': n_used_map,
        'fit_ok_map': fit_ok_map
    }
    
    return cube_baseline, stats
