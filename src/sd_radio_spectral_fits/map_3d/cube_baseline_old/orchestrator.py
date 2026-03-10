# src/sd_radio_spectral_fits/cube_baseline/orchestrator.py
import numpy as np
from typing import Optional, List, Union, Tuple
from .session import BaselineSession
from .engine import fit_cube_baseline
from ..map.cube_analysis import estimate_robust_rms, generate_cube_mask
from ..ranges import parse_windows, in_any_windows

def create_manual_mask_1d(v_axis: np.ndarray, v_windows: Union[List[str], List[Tuple[float, float]]]) -> np.ndarray:
    """手動指定の速度窓から1DのBooleanマスク(True=シグナル)を生成する"""
    if not v_windows:
        return np.zeros(len(v_axis), dtype=bool)
    windows = parse_windows(v_windows) if isinstance(v_windows[0], str) else v_windows
    return in_any_windows(v_axis, windows)

def run_one_iteration(
    session: BaselineSession,
    *,
    auto_mask_method: Optional[str] = "derivative",
    auto_mask_kwargs: Optional[dict] = None,
    manual_v_windows: Optional[Union[List[str], List[Tuple[float, float]]]] = None,
    poly_order: int = 1,
    target_mask_2d: Optional[np.ndarray] = None,
) -> None:
    """
    反復ベースライン処理の1サイクルを実行する。
    自動マスク生成 -> 手動マスク合成 -> フィッティング -> Sessionの更新(in-place)
    """
    if auto_mask_kwargs is None:
        auto_mask_kwargs = {}
        
    if target_mask_2d is None:
        ny, nx = session.cube_original.shape[:2]
        target_mask_2d = np.ones((ny, nx), dtype=bool)

    # 1. 自動マスクの生成 (現在の補正済みキューブに対して実行)
    if auto_mask_method:
        current_cube = session.get_full_cube_work()
        rms_map = estimate_robust_rms(current_cube)
        mask_3d = generate_cube_mask(current_cube, rms_map, method=auto_mask_method, **auto_mask_kwargs)
        
        session.update_mask_3d(mask_3d, params={"method": auto_mask_method, "kwargs": auto_mask_kwargs})
    else:
        mask_3d = session.mask_3d_current

    # 事故防止: 自動マスクも手動マスクも一切ない場合、全領域がフィットされてしまうのを防ぐ
    if mask_3d is None and not manual_v_windows:
        raise ValueError("No signal mask provided. Specify auto_mask_method or manual_v_windows to avoid fitting the entire spectrum including lines.")

    # 2. 手動マスク(1D)の準備
    manual_mask_1d = create_manual_mask_1d(session.v_axis, manual_v_windows) if manual_v_windows else None

    # 3. CubeBaselineEngine の実行 (常にオリジナルデータに対してフィット)
    new_baseline, stats = fit_cube_baseline(
        cube_original=session.cube_original,
        v_axis=session.v_axis,
        signal_mask_3d=mask_3d,
        manual_mask_1d=manual_mask_1d,
        target_mask_2d=target_mask_2d,
        poly_order=poly_order
    )

    # 4. Session の in-place 更新
    params = {
        "poly_order": poly_order,
        "manual_v_windows": manual_v_windows,
        "auto_mask_method": auto_mask_method
    }
    session.update_baseline(new_baseline, target_mask_2d, stats, params)
