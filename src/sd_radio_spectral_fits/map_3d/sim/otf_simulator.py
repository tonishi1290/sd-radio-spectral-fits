# src/sd_radio_spectral_fits/otf/sim/otf_simulator.py
from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Tuple

from ...fitsio import Scantable

# ==========================================
# 1. 天体モデルの定義
# ==========================================

@dataclass
class GaussianSource:
    """空間的・速度的な広がりを持つ1つのガウシアン雲（コンポーネント）"""
    ra_offset_arcsec: float
    dec_offset_arcsec: float
    sig_major_arcsec: float
    sig_minor_arcsec: float
    pa_deg: float
    v_center_kms: float
    v_width_kms: float
    peak_t: float


class SkyModel:
    """複数のガウシアン雲を合成した天体モデル"""
    
    def __init__(self):
        self.sources: List[GaussianSource] = []
        
    def add_source(self, src: GaussianSource) -> None:
        self.sources.append(src)
        
    def evaluate(
        self, 
        ra_deg: np.ndarray, 
        dec_deg: np.ndarray, 
        v_axis: np.ndarray, 
        ref_ra: float, 
        ref_dec: float
    ) -> np.ndarray:
        """指定された座標と速度軸における、モデルの理論スペクトル(強度)を計算する"""
        spec = np.zeros((len(ra_deg), len(v_axis)), dtype=float)
        
        # 基準点からのオフセット [arcsec]
        dx_base = (ra_deg - ref_ra) * np.cos(np.deg2rad(ref_dec)) * 3600.0
        dy_base = (dec_deg - ref_dec) * 3600.0
        
        for src in self.sources:
            dx = dx_base - src.ra_offset_arcsec
            dy = dy_base - src.dec_offset_arcsec
            
            # 回転（位置角）
            theta = np.deg2rad(src.pa_deg)
            x_rot = dx * np.cos(theta) + dy * np.sin(theta)
            y_rot = -dx * np.sin(theta) + dy * np.cos(theta)
            
            # 空間プロファイル
            spatial = np.exp(-0.5 * ((x_rot / src.sig_major_arcsec)**2 + (y_rot / src.sig_minor_arcsec)**2))
            
            # 速度プロファイル
            spectral = np.exp(-0.5 * ((v_axis - src.v_center_kms) / src.v_width_kms)**2)
            
            # (N_pos, 1) * (1, N_vel) -> (N_pos, N_vel)
            spec += src.peak_t * (spatial[:, None] * spectral[None, :])
            
        return spec


def create_complex_sky_model() -> SkyModel:
    """
    大きく複雑な天体（巨大分子雲複合体のようなもの）のプリセット。
    空間構造だけでなく、速度勾配（成分ごとに速度が違う）も持つ。
    """
    model = SkyModel()
    
    # 1. Broad Halo (広がった淡いガス)
    model.add_source(GaussianSource(0, 0, 2000, 1500, 20, v_center_kms=10.0, v_width_kms=5.0, peak_t=1.5))
    
    # 2. Main Filament (中心を貫く高密度のフィラメント)
    model.add_source(GaussianSource(0, 0, 3000, 250, -45, v_center_kms=12.0, v_width_kms=2.0, peak_t=6.0))
    
    # 3. Dense Core (フィラメント内にある非常に明るいコア)
    model.add_source(GaussianSource(-500, 500, 200, 200, 0, v_center_kms=12.5, v_width_kms=1.0, peak_t=10.0))
    
    return model

# ==========================================
# 2. 観測シミュレータ本体
# ==========================================

def simulate_otf_observation(
    sky_model: SkyModel,
    ra_center: float = 200.0,
    dec_center: float = -30.0,
    map_width_arcsec: float = 1800.0,
    map_height_arcsec: float = 1800.0,
    scan_spacing_arcsec: float = 15.0,
    scan_speed_arcsec_s: float = 60.0,
    dump_interval_s: float = 0.25,
    noise_rms: float = 0.8,
    stripe_amplitude: float = 1.5,
    vframe_amplitude: float = 5.0,
    beam_offsets_arcsec: List[Tuple[float, float]] = None,
    beam_pa_deg: float = 0.0,
    pointing_jitter_arcsec: float = 0.0,
    scan_direction: str = "X",
) -> Scantable:
    """
    指定された天体モデルを元に、リアルなOTF観測（ラスター走査）をシミュレートする。
    """
    print(f"Simulating OTF {scan_direction.upper()}-scan trajectory (Boresight)...")
    
    if beam_offsets_arcsec is None:
        beam_offsets_arcsec = [(0.0, 0.0)]

    # --- 走査方向に応じて、高速軸・低速軸の長さを決定 ---
    if scan_direction.upper() == "X":
        fast_length = map_width_arcsec
        slow_length = map_height_arcsec
    else:
        fast_length = map_height_arcsec
        slow_length = map_width_arcsec

    # グリッドの分割数を計算
    dx_fast = scan_speed_arcsec_s * dump_interval_s
    n_fast = int(fast_length / dx_fast)
    n_slow = int(slow_length / scan_spacing_arcsec)
    
    fast_grid = np.linspace(-fast_length / 2, fast_length / 2, n_fast)
    slow_grid = np.linspace(-slow_length / 2, slow_length / 2, n_slow)
    
    ra_list, dec_list, scan_list, turn_list, time_list = [], [], [], [], []
    scan_id = 0
    current_time = 0.0
    cos_dec = np.cos(np.deg2rad(dec_center))
    
    # 1. Boresight (望遠鏡の中心座標) 軌跡の計算
    for i, slow_val in enumerate(slow_grid):
        # 行ごとに左右(上下)を反転させてジグザグ走査を再現
        row_fast = fast_grid if i % 2 == 0 else fast_grid[::-1]
        
        # ターンアラウンド（スキャン両端の加速・減速区間）のフラグ立て
        turn = np.zeros(n_fast, dtype=bool)
        margin = int(n_fast * 0.05)
        turn[:margin] = True
        turn[-margin:] = True
        
        # XスキャンとYスキャンで座標を割り当て
        if scan_direction.upper() == "X":
            x_pts, y_pts = row_fast, np.full(n_fast, slow_val)
        else:
            x_pts, y_pts = np.full(n_fast, slow_val), row_fast
            
        ra_list.append(ra_center + (x_pts / 3600.0) / cos_dec)
        dec_list.append(dec_center + (y_pts / 3600.0))
        scan_list.append(np.full(n_fast, scan_id))
        turn_list.append(turn)
        
        # 時刻の更新 (スキャン間に5秒のオーバーヘッドを仮定)
        t_arr = current_time + np.arange(n_fast) * dump_interval_s
        time_list.append(t_arr)
        current_time = t_arr[-1] + 5.0
        scan_id += 1
        
    ra_bore = np.concatenate(ra_list)
    dec_bore = np.concatenate(dec_list)
    scan_bore = np.concatenate(scan_list)
    turn_bore = np.concatenate(turn_list)
    time_bore = np.concatenate(time_list)
    
    # 2. マルチビーム展開と回転・誤差の付加
    print(f"Expanding for {len(beam_offsets_arcsec)} beams with PA={beam_pa_deg}deg...")
    ra_all, dec_all, scan_all, turn_all, time_all, beam_id_all = [], [], [], [], [], []
    
    # 回転角の準備
    theta = np.deg2rad(beam_pa_deg)
    cos_pa = np.cos(theta)
    sin_pa = np.sin(theta)
    
    for b_idx, (bx, by) in enumerate(beam_offsets_arcsec):
        # ビーム配置の回転 (Rotation Matrix)
        bx_rot = bx * cos_pa - by * sin_pa
        by_rot = bx * sin_pa + by * cos_pa
        
        # Boresight にオフセットを加算
        ra_beam = ra_bore + (bx_rot / 3600.0) / cos_dec
        dec_beam = dec_bore + (by_rot / 3600.0)
        
        # ポインティングジッター（揺らぎ）の付加
        if pointing_jitter_arcsec > 0:
            np.random.seed(42 + b_idx)
            ra_beam += np.random.normal(0, pointing_jitter_arcsec / 3600.0, len(ra_beam)) / cos_dec
            dec_beam += np.random.normal(0, pointing_jitter_arcsec / 3600.0, len(dec_beam))
            
        ra_all.append(ra_beam)
        dec_all.append(dec_beam)
        
        # ビームごとに一意のSCAN IDを付与（BoresightのSCAN数分だけ下駄を履かせる）
        scan_all.append(scan_bore + (b_idx * (scan_id + 1)))
        turn_all.append(turn_bore)
        time_all.append(time_bore)
        beam_id_all.append(np.full(len(ra_beam), b_idx))

    ra_arr = np.concatenate(ra_all)
    dec_arr = np.concatenate(dec_all)
    scan_arr = np.concatenate(scan_all)
    turn_arr = np.concatenate(turn_all)
    time_arr = np.concatenate(time_all)
    beam_arr = np.concatenate(beam_id_all)
    
    ndump = len(ra_arr)
    nchan = 200
    v_axis = np.linspace(-20, 40, nchan)
    
    # 3. 理論スペクトルの計算
    print(f"Evaluating sky model for {ndump} total dumps...")
    spec_data = sky_model.evaluate(ra_arr, dec_arr, v_axis, ra_center, dec_center)
    
    # 4. ノイズとストライプ（走査ごとのベースライン変動）の付加
    print("Adding noise, stripes, and Doppler shifts...")
    np.random.seed(42)
    spec_data += np.random.normal(0, noise_rms, (ndump, nchan))
    
    # スキャンごとに一定のオフセット（ストライプ）を乗せる
    stripes = np.random.normal(0, stripe_amplitude, int(np.max(scan_arr)) + 1)
    spec_data += stripes[scan_arr][:, None]
    
    # 5. 地球の自転・公転による視線速度変化（ドップラーシフト）の付加
    vframe_arr = vframe_amplitude * np.sin(time_arr / time_arr[-1] * 2 * np.pi)
    dv = v_axis[1] - v_axis[0]
    shift_bins = np.rint(vframe_arr / dv).astype(int)
    
    for i in range(ndump):
        spec_data[i] = np.roll(spec_data[i], -shift_bins[i])
        
    # 6. Scantable 形式に組み立てて出力
    print("Building Scantable...")
    table = pd.DataFrame({
        "RA": ra_arr, "DEC": dec_arr,
        "VFRAME": vframe_arr, "MJD": 60000.0 + time_arr / 86400.0,
        "SCAN": scan_arr, "IS_TURN": turn_arr, "BEAM": beam_arr,
        "OBSMODE": ["ON"] * ndump, "BSL_RMS": np.full(ndump, noise_rms)
    })
    
    # 周波数軸の設定 (115GHz CO(1-0) を想定)
    restfreq = 115271201800.0
    c_kms = 299792.458
    f_axis = restfreq * (1.0 - v_axis / c_kms)
    df_freq = f_axis[1] - f_axis[0]
    
    meta = {
        "RESTFRQ": restfreq, 
        "CRVAL1": float(f_axis[0]), 
        "CDELT1": float(df_freq),
        "CTYPE1": "FREQ", 
        "CUNIT1": "Hz", 
        "CRPIX1": 1.0, 
        "SPECSYS": "TOPO"
    }
    
    return Scantable(meta=meta, data=spec_data, table=table)
