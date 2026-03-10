# src/sd_radio_spectral_fits/atmosphere.py
from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Optional, Union, Tuple

# 宇宙背景放射の温度 [K]
from .constants import T_BG_K

def get_airmass(elevation_deg: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    仰角 (Elevation) からエアマス (sec z) を計算する。
    
    Args:
        elevation_deg (float | ndarray): 望遠鏡の仰角 [度]。
        
    Returns:
        airmass (float | ndarray): エアマス X = sec(z)。
        ※低仰角 (EL < 5度) の場合、物理的に非現実的な値になるのを防ぐため、
        sin(5度) 相当のエアマス値にクリッピング（制限）されます。
    """
    el_rad = np.deg2rad(elevation_deg)
    # 0割り防止のため、最小仰角を5度に制限
    sin_el = np.clip(np.sin(el_rad), np.sin(np.deg2rad(5.0)), 1.0)
    airmass = 1.0 / sin_el
    return airmass

def estimate_t_atm(
    t_surface_k: float,
    model: str = "offset",
    delta_t: float = 15.0,
    eta: float = 0.95
) -> float:
    """
    地表の気温から、空の実効的な物理温度 (T_atm) を推定する経験則モデル。
    対流圏の上空の冷たい大気層からの放射も含むため、T_atm は地表気温より低くなります。
    
    Args:
        t_surface_k (float): 屋外の地表気温（外気温）[K]。
        model (str): 推定モデル。"offset" (引き算) または "ratio" (掛け算)。
        delta_t (float): offset モデル使用時の温度差 [K] (デフォルト: 15.0 K)。
        eta (float): ratio モデル使用時の係数 (デフォルト: 0.95)。
        
    Returns:
        t_atm (float): 推定された空の実効温度 [K]。
    """
    model = model.strip().lower()
    if model == "offset":
        return t_surface_k - delta_t
    elif model == "ratio":
        return t_surface_k * eta
    else:
        raise ValueError(f"Unknown T_atm model: '{model}'. Use 'offset' or 'ratio'.")

def compute_t_cal_array(
    t_hot_k: float,
    t_atm_k: float,
    tau_zenith: float,
    elevation_deg: np.ndarray,
    t_bg_k: float = T_BG_K
) -> np.ndarray:
    """
    Ulich & Haas (1976) の大気モデルに基づき、各データ行（仰角）ごとの
    厳密な等価キャリブレーション温度 (T_cal) を計算する。
    
    T_cal とは、「Chopper-Wheel 法で得られた (Hot - Sky) の強度差が、
    大気外 (天体と同じ位置) に置かれた場合、何度 [K] の信号に相当するか」を表す係数です。
    
    Args:
        t_hot_k (float): ホットロードの物理温度 [K]。
        t_atm_k (float): 空の実効物理温度 [K]。estimate_t_atm() で求めた値。
        tau_zenith (float): 天頂方向の光学的厚み (tau)。
        elevation_deg (ndarray): 各データ行の仰角配列 [度]。
        t_bg_k (float): 宇宙背景放射の温度 [K]。
        
    Returns:
        t_cal_array (ndarray): 各データ行に対応する T_cal の配列 [K]。
                               仰角の時間変化に伴い、行ごとに値が変動します。
    """
    # 1. 各行のエアマス X を計算
    x_arr = get_airmass(elevation_deg)
    
    # 2. 大気による減衰の逆数 (補正係数) e^(tau * X)
    exp_tau_x = np.exp(tau_zenith * x_arr)
    
    # 3. 厳密な T_cal の計算式
    t_cal_array = t_hot_k * exp_tau_x - t_atm_k * (exp_tau_x - 1.0) - t_bg_k
    
    return t_cal_array

def extract_meta_value(
    mapping_df: pd.DataFrame,
    meta_dict: dict,
    valid_keys: Tuple[str, ...]
) -> Optional[float]:
    """
    FITSのテーブル (mapping) または 全体メタデータ (meta) から、
    エイリアス（カラム名の揺らぎ）を許容して数値を安全に抽出するユーティリティ。
    
    Args:
        mapping_df (pd.DataFrame): 観測データのテーブル (ON点など)。
        meta_dict (dict): 全体のメタデータ辞書 (ヘッダ)。
        valid_keys (Tuple[str]): 検索を許容するキー名の優先順位タプル。
                                 例: ("TAU0", "TAU", "OPACITY")
        
    Returns:
        見つかった場合はその数値 (float)。見つからない場合は None。
    """
    # 1. テーブル (mapping) のカラムを優先して探す (時間変化パラメータなど)
    for key in valid_keys:
        if key in mapping_df.columns:
            # np.nan 等を排除し、最初の有効な値を取得
            val = pd.to_numeric(mapping_df[key], errors="coerce").dropna().values
            if len(val) > 0:
                return float(val[0])
                
    # 2. 全体メタデータ (meta) を探す (固定パラメータなど)
    for key in valid_keys:
        if key in meta_dict:
            try:
                return float(meta_dict[key])
            except (ValueError, TypeError):
                continue
                
    return None


def extract_meta_array(
    mapping_df: pd.DataFrame,
    meta_dict: dict,
    valid_keys: Tuple[str, ...]
) -> Union[np.ndarray, float, None]:
    """
    時間変化（行ごと）のパラメータを配列として取得するユーティリティ。
    
    Args:
        mapping_df (pd.DataFrame): 観測データのテーブル (ON点など)。
        meta_dict (dict): 全体のメタデータ辞書 (ヘッダ)。
        valid_keys (Tuple[str]): 検索を許容するキー名の優先順位タプル。
        
    Returns:
        テーブルにカラムがあればその配列(ndarray)を、無ければメタデータのスカラー値(float)を返す。
        どちらにも無ければ None を返す。
    """
    # 1. テーブルから時系列配列を探す
    for key in valid_keys:
        if key in mapping_df.columns:
            arr = pd.to_numeric(mapping_df[key], errors="coerce").to_numpy(dtype=float)
            # 全てNaNでなければ採用（一部の欠損値は前後の値で補間）
            if not np.all(np.isnan(arr)):
                s = pd.Series(arr).ffill().bfill()
                return s.to_numpy(dtype=float)
                
    # 2. 全体メタデータからスカラー値を探す
    for key in valid_keys:
        if key in meta_dict:
            try:
                return float(meta_dict[key])
            except (ValueError, TypeError):
                continue
                
    return None
