# src/sd_radio_spectral_fits/cube_baseline/session.py
import numpy as np
from typing import Optional, Dict, Any, List

class BaselineSession:
    """ベースライン再処理の状態を一元管理するセッションクラス。"""
    
    def __init__(self, cube_original: np.ndarray, v_axis: np.ndarray):
        self.cube_original = cube_original
        # memmap や view の場合に書き込み不可フラグの変更で例外が出るのを防ぐ
        try:
            self.cube_original.flags.writeable = False
        except Exception:
            pass
            
        self.v_axis = v_axis
        
        self.baseline_cube = np.zeros(cube_original.shape, dtype=np.float32)
        self._cube_work_cache: Optional[np.ndarray] = None
        
        self.mask_3d_current: Optional[np.ndarray] = None
        self.target_mask_2d: Optional[np.ndarray] = None
        
        self.fit_stats: Dict[str, np.ndarray] = {}
        self.edit_history: List[Dict[str, Any]] = []
        self.grid_provenance = None

    def get_corrected_spectrum(self, iy: int, ix: int) -> np.ndarray:
        return self.cube_original[iy, ix, :] - self.baseline_cube[iy, ix, :]

    def get_corrected_subcube(self, y_slice: slice, x_slice: slice) -> np.ndarray:
        return self.cube_original[y_slice, x_slice, :] - self.baseline_cube[y_slice, x_slice, :]

    def get_full_cube_work(self) -> np.ndarray:
        if self._cube_work_cache is None:
            self._cube_work_cache = self.cube_original - self.baseline_cube
        return self._cube_work_cache

    def update_baseline(self, new_baseline: np.ndarray, target_mask_2d: np.ndarray, stats: Dict[str, np.ndarray], params: Dict[str, Any]):
        self.baseline_cube[target_mask_2d, :] = new_baseline[target_mask_2d, :]
        
        for key, new_stat_map in stats.items():
            if key not in self.fit_stats:
                dtype = new_stat_map.dtype
                fill_val = np.nan if np.issubdtype(dtype, np.floating) else 0
                self.fit_stats[key] = np.full(target_mask_2d.shape, fill_val, dtype=dtype)
            self.fit_stats[key][target_mask_2d] = new_stat_map[target_mask_2d]
            
        self._cube_work_cache = None
        self.edit_history.append({"action": "fit_baseline", "params": params})
        
    def update_mask_3d(self, new_mask: np.ndarray, params: Dict[str, Any]):
        self.mask_3d_current = new_mask
        self.edit_history.append({"action": "update_mask_3d", "params": params})
