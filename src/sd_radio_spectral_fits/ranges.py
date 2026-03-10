# src/sd_radio_spectral_fits/ranges.py
from __future__ import annotations

from typing import List, Tuple, Optional
import numpy as np

def parse_windows(specs: List[str]) -> List[Tuple[float, float]]:
    """
    'min:max' 形式の文字列リストを数値のタプルリストにパースする。
    
    Example:
        [' -100 : -50 ', '50:100'] -> [(-100.0, -50.0), (50.0, 100.0)]
    """
    if not specs:
        return []
        
    out: List[Tuple[float, float]] = []
    for s in specs:
        if not s or ":" not in s:
            continue
        try:
            # 空白を除去して分割
            parts = s.split(":")
            if len(parts) != 2:
                continue
            v1, v2 = float(parts[0]), float(parts[1])
            # 常に (小さい値, 大きい値) の順にする
            out.append((min(v1, v2), max(v1, v2)))
        except (ValueError, IndexError):
            continue
    return out

def window_to_mask(x_axis: np.ndarray, windows: List[Tuple[float, float]]) -> np.ndarray:
    """
    座標軸（x_axis）に対して、指定された範囲内（windows）をTrueとするマスクを作成する。
    
    Args:
        x_axis: 速度(km/s)や周波数の配列。
        windows: [(min, max), ...] 形式の数値タプルリスト。
    """
    x = np.asarray(x_axis, float)
    mask = np.zeros_like(x, dtype=bool)
    
    if not windows:
        return mask
        
    for a, b in windows:
        # 指定範囲内（境界含む）をTrueに設定
        mask |= (x >= a) & (x <= b)
        
    return mask

# --- 互換性のためのエイリアス ---
def windows_to_mask(x: np.ndarray, windows: List[Tuple[float, float]]) -> np.ndarray:
    """以前の複数形名での呼び出しにも対応する。"""
    return window_to_mask(x, windows)
