# src/sd_radio_spectral_fits/otf/basketweave.py
from __future__ import annotations

import logging
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree

from .config import GridInput
from ..utils import parse_windows, in_any_windows


def _safe_bool_array(values, default: bool = False) -> np.ndarray:
    """Convert mixed boolean-like values into a strict bool ndarray."""
    a = np.asarray(values)
    if a.dtype.kind == "b":
        return a.astype(bool, copy=False)

    out = np.full(a.shape, bool(default), dtype=bool)
    if a.dtype.kind in "iuf":
        finite = np.isfinite(a)
        out[finite] = a[finite] != 0
        return out

    flat = a.astype(object, copy=False).ravel()
    flat_out = out.ravel()
    true_vals = {"true", "t", "1", "yes", "y", "on"}
    false_vals = {"false", "f", "0", "no", "n", "off", "", "nan", "none", "null", "<na>"}
    for i, v in enumerate(flat):
        if v is None:
            flat_out[i] = bool(default)
            continue
        s = str(v).strip().lower()
        if s in true_vals:
            flat_out[i] = True
        elif s in false_vals:
            flat_out[i] = False
        else:
            flat_out[i] = bool(default)
    return out


def _safe_scan_ids_int64(values) -> np.ndarray:
    """Convert scan ids to int64, mapping invalid entries to -1."""
    arr = np.asarray(values)
    if arr.size == 0:
        return np.array([], dtype=np.int64)

    num = np.empty(arr.shape, dtype=np.float64)
    num.fill(np.nan)
    try:
        num = arr.astype(np.float64)
    except (TypeError, ValueError):
        flat_src = arr.astype(object, copy=False).ravel()
        flat_dst = num.ravel()
        for i, v in enumerate(flat_src):
            try:
                flat_dst[i] = float(v)
            except (TypeError, ValueError):
                flat_dst[i] = np.nan

    out = np.full(num.shape, -1, dtype=np.int64)
    finite = np.isfinite(num)
    out[finite] = num[finite].astype(np.int64)
    return out


def _zero_offsets_from_scan_ids(scan_ids_all: np.ndarray) -> np.ndarray:
    max_scan = int(np.max(scan_ids_all)) if np.any(scan_ids_all >= 0) else -1
    return np.zeros(max_scan + 1, dtype=float)


def _sorted_pairs_array(pairs_arr) -> np.ndarray:
    arr = np.asarray(pairs_arr, dtype=np.int64)
    if arr.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    arr = arr.reshape(-1, 2)
    order = np.lexsort((arr[:, 1], arr[:, 0]))
    return arr[order]


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
    if search_radius_arcsec <= 0:
        raise ValueError("search_radius_arcsec must be positive")
    if damp < 0:
        raise ValueError("damp must be non-negative")

    x_all = np.asarray(input_data.x, dtype=float)
    y_all = np.asarray(input_data.y, dtype=float)
    scan_ids_all = input_data.scan_id

    if scan_ids_all is None:
        raise ValueError("scan_id must be provided in GridInput for basket-weaving.")

    scan_ids_all = _safe_scan_ids_int64(scan_ids_all)
    spec = np.asarray(input_data.spec)

    # --- 1. 速度ウィンドウからのチャネルマスク生成 ---
    if v_windows_kms is not None:
        if v_axis is None:
            raise ValueError("v_axis must be provided when v_windows_kms is specified.")

        if len(v_windows_kms) > 0 and isinstance(v_windows_kms[0], str):
            windows = parse_windows(v_windows_kms)
        else:
            windows = v_windows_kms

        channel_mask = in_any_windows(v_axis, windows)

        if not np.any(channel_mask):
            logging.error(f"Basket-weave aborted: No channels found in windows {v_windows_kms}.")
            logging.error(f"(Available velocity range is approx {np.nanmin(v_axis):.1f} to {np.nanmax(v_axis):.1f} km/s)")
            raise ValueError("Invalid v_windows_kms: No matching channels. Please check your velocity range.")

    if channel_mask is not None:
        channel_mask = np.asarray(channel_mask, dtype=bool)
        if spec.ndim != 2:
            raise ValueError("channel_mask can only be used when input_data.spec is 2D.")
        if channel_mask.ndim != 1 or channel_mask.shape[0] != spec.shape[1]:
            raise ValueError(
                f"channel_mask must have shape ({spec.shape[1]},), got {channel_mask.shape}"
            )
        if not np.any(channel_mask):
            raise ValueError("channel_mask selected zero channels.")

    # --- 2. スペクトル平均の計算 ---
    if spec.ndim == 2:
        if channel_mask is not None:
            d_all = np.nanmean(spec[:, channel_mask], axis=1)
        else:
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
        good &= ~_safe_bool_array(input_data.is_turnaround, default=False)
    good &= np.isfinite(x_all) & np.isfinite(y_all) & np.isfinite(d_all)
    good &= scan_ids_all >= 0

    if np.count_nonzero(good) < 2:
        return _zero_offsets_from_scan_ids(scan_ids_all)

    x, y, d, scan_ids = x_all[good], y_all[good], d_all[good], scan_ids_all[good]
    uniq_scan_ids, inv_scan = np.unique(scan_ids, return_inverse=True)
    num_scans_compact = int(uniq_scan_ids.size)

    # --- 4. 交差点の探索と行列方程式の構築 ---
    tree = cKDTree(np.c_[x, y])
    try:
        pairs_arr = tree.query_pairs(r=float(search_radius_arcsec), output_type="ndarray")
    except TypeError:
        pairs_set = tree.query_pairs(r=float(search_radius_arcsec))
        if not pairs_set:
            return _zero_offsets_from_scan_ids(scan_ids_all)
        pairs_arr = np.array(list(pairs_set), dtype=np.int64)

    pairs_arr = _sorted_pairs_array(pairs_arr)
    if len(pairs_arr) == 0:
        logging.warning("No cross points found for basket-weaving. Returning zero offsets.")
        return _zero_offsets_from_scan_ids(scan_ids_all)

    i_idx, j_idx = pairs_arr[:, 0], pairs_arr[:, 1]
    scan_i, scan_j = inv_scan[i_idx], inv_scan[j_idx]
    valid_mask = scan_i != scan_j
    i_idx, j_idx, scan_i, scan_j = i_idx[valid_mask], j_idx[valid_mask], scan_i[valid_mask], scan_j[valid_mask]

    if len(i_idx) == 0:
        return _zero_offsets_from_scan_ids(scan_ids_all)

    mean_diff = d[i_idx] - d[j_idx]
    eq_idx = len(i_idx)
    row_indices = np.repeat(np.arange(eq_idx), 2)
    col_indices = np.stack([scan_i, scan_j], axis=1).ravel()
    data_vals = np.zeros(eq_idx * 2, dtype=float)
    data_vals[0::2], data_vals[1::2] = 1.0, -1.0

    A = sp.coo_matrix((data_vals, (row_indices, col_indices)), shape=(eq_idx, num_scans_compact))

    # --- 5. 最小二乗法でオフセットを推定 ---
    res = lsqr(A, mean_diff, damp=damp)
    offsets_compact = np.asarray(res[0], dtype=float)
    offsets_compact -= np.mean(offsets_compact)

    num_scans_full = int(np.max(scan_ids)) + 1 if scan_ids.size else 0
    offsets_full = np.zeros(num_scans_full, dtype=float)
    offsets_full[uniq_scan_ids] = offsets_compact
    return offsets_full


def apply_basket_weave_correction(input_data: GridInput, offsets: np.ndarray) -> None:
    """
    推定されたオフセットを入力データ (GridInput.spec) から差し引く。
    """
    if input_data.scan_id is None:
        raise ValueError("scan_id is required in GridInput to apply corrections.")

    spec = np.asarray(input_data.spec)
    if spec.ndim != 2:
        raise ValueError(f"input_data.spec must be 2D, got shape={spec.shape}")
    if not np.issubdtype(spec.dtype, np.floating):
        spec = spec.astype(np.float64)
        input_data.spec = spec

    scan_ids = _safe_scan_ids_int64(input_data.scan_id)
    valid_offset_mask = (scan_ids >= 0) & (scan_ids < len(offsets))

    correction_vector = np.zeros(len(scan_ids), dtype=spec.dtype)
    correction_vector[valid_offset_mask] = offsets[scan_ids[valid_offset_mask]].astype(spec.dtype, copy=False)

    input_data.spec -= correction_vector[:, np.newaxis]
