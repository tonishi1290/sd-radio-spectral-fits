# src/sd_radio_spectral_fits/otf/basketweave.py
from __future__ import annotations

import logging
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree

from .config import GridInput
from ..utils import parse_windows, in_any_windows

def solve_basket_weave_offsets(
    input_data: GridInput,
    search_radius_arcsec: float = 3.0,
    damp: float = 1e-2,
    v_axis: np.ndarray | None = None,
    v_windows_kms: list[str] | list[tuple[float, float]] | None = None,
    channel_mask: np.ndarray | None = None,
) -> np.ndarray:
    """
    Least-squares basket-weaving法を用いて、スキャンごとのベースラインオフセットを推定する。
    """
    x_all = np.asarray(input_data.x, dtype=float)
    y_all = np.asarray(input_data.y, dtype=float)
    scan_ids_all = input_data.scan_id

    if scan_ids_all is None:
        raise ValueError("scan_id must be provided in GridInput for basket-weaving.")
    
    scan_ids_all = np.asarray(scan_ids_all).astype(np.int64, copy=False)
    spec = np.asarray(input_data.spec)

    # --- 1. 速度ウィンドウからのチャネルマスク生成 (utils.py を活用) ---
    if v_windows_kms is not None:
        if v_axis is None:
            raise ValueError("v_axis must be provided when v_windows_kms is specified.")
        
        # 文字列のリスト ["-30:-5", "30:55"] の場合は utils.py でパース
        if len(v_windows_kms) > 0 and isinstance(v_windows_kms[0], str):
            windows = parse_windows(v_windows_kms)
        else:
            windows = v_windows_kms  # 既にタプルのリストの場合
            
        # utils.py の機能を使って複数ウィンドウの boolean マスクを一括生成
        channel_mask = in_any_windows(v_axis, windows)
        
        # 安全対策: 指定範囲にチャネルが1つも存在しない場合
        if not np.any(channel_mask):
            logging.error(f"Basket-weave aborted: No channels found in windows {v_windows_kms}.")
            logging.error(f"(Available velocity range is approx {np.nanmin(v_axis):.1f} to {np.nanmax(v_axis):.1f} km/s)")
            raise ValueError("Invalid v_windows_kms: No matching channels. Please check your velocity range.")

    # --- 2. スペクトル平均の計算 ---
    if spec.ndim == 2:
        if channel_mask is not None:
            # マスクされた Line-free 領域のみで平均を計算
            d_all = np.nanmean(spec[:, channel_mask], axis=1)
        else:
            # 範囲指定がない場合は警告を出して全チャネル平均
            logging.warning(
                "Basket-weaving is averaging across ALL channels. "
                "Strong emission lines might bias the baseline offset estimation. "
                "Consider specifying 'v_windows_kms' (Line-free range) for better accuracy."
            )
            d_all = np.nanmean(spec, axis=1)
    else:
        d_all = spec
        
    d_all = np.asarray(d_all, dtype=float)

    # --- 3. データ品質フラグのチェック ---
    good = np.ones_like(x_all, dtype=bool)
    if input_data.flag is not None:
        good &= np.asarray(input_data.flag) > 0
    if input_data.is_turnaround is not None:
        good &= ~np.asarray(input_data.is_turnaround, dtype=bool)
    good &= np.isfinite(x_all) & np.isfinite(y_all) & np.isfinite(d_all)

    if np.count_nonzero(good) < 2:
        return np.zeros(int(np.max(scan_ids_all)) + 1, dtype=float)

    x, y, d, scan_ids = x_all[good], y_all[good], d_all[good], scan_ids_all[good]
    uniq_scan_ids, inv_scan = np.unique(scan_ids, return_inverse=True)
    num_scans_compact = int(uniq_scan_ids.size)

    # --- 4. 交差点の探索と行列方程式の構築 ---
    tree = cKDTree(np.c_[x, y])
    try:
        pairs_arr = tree.query_pairs(r=float(search_radius_arcsec), output_type='ndarray')
    except TypeError:
        pairs_set = tree.query_pairs(r=float(search_radius_arcsec))
        if not pairs_set: 
            return np.zeros(int(np.max(scan_ids_all)) + 1, dtype=float)
        pairs_arr = np.array(list(pairs_set), dtype=np.int32)

    if len(pairs_arr) == 0:
        logging.warning("No cross points found for basket-weaving. Returning zero offsets.")
        return np.zeros(int(np.max(scan_ids_all)) + 1, dtype=float)
        
    i_idx, j_idx = pairs_arr[:, 0], pairs_arr[:, 1]
    scan_i, scan_j = inv_scan[i_idx], inv_scan[j_idx]
    valid_mask = (scan_i != scan_j)
    i_idx, j_idx, scan_i, scan_j = i_idx[valid_mask], j_idx[valid_mask], scan_i[valid_mask], scan_j[valid_mask]

    if len(i_idx) == 0:
        return np.zeros(int(np.max(scan_ids_all)) + 1, dtype=float)

    mean_diff = d[i_idx] - d[j_idx]
    eq_idx = len(i_idx)
    row_indices = np.repeat(np.arange(eq_idx), 2)
    col_indices = np.stack([scan_i, scan_j], axis=1).ravel()
    data_vals = np.zeros(eq_idx * 2, dtype=float)
    data_vals[0::2], data_vals[1::2] = 1.0, -1.0

    A = sp.coo_matrix((data_vals, (row_indices, col_indices)), shape=(eq_idx, num_scans_compact))

    # --- 5. 最小二乗法でオフセットを推定 ---
    res = lsqr(A, mean_diff, damp=damp)
    offsets_compact = res[0]
    offsets_compact -= np.mean(offsets_compact)

    num_scans_full = int(np.max(scan_ids_all)) + 1
    offsets_full = np.zeros(num_scans_full, dtype=float)
    offsets_full[uniq_scan_ids] = offsets_compact
    
    return offsets_full

def apply_basket_weave_correction(input_data: GridInput, offsets: np.ndarray) -> None:
    """
    推定されたオフセットを入力データ (GridInput.spec) から差し引く。
    """
    if input_data.scan_id is None:
        raise ValueError("scan_id is required in GridInput to apply corrections.")

    # 整数型のスキャンIDとして扱う
    scan_ids = np.asarray(input_data.scan_id).astype(np.int64, copy=False)
    
    # オフセット配列から各行に対応する値を取得
    valid_offset_mask = (scan_ids >= 0) & (scan_ids < len(offsets))
    
    # 補正値を適用（スペクトル全体から各行のオフセットを引く）
    correction_vector = np.zeros(len(scan_ids))
    correction_vector[valid_offset_mask] = offsets[scan_ids[valid_offset_mask]]
    
    # spec データの各チャンネルからオフセットを減算
    input_data.spec -= correction_vector[:, np.newaxis]
