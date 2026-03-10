# src/sd_radio_spectral_fits/simulation/grid_sim.py
from __future__ import annotations

import datetime
from typing import Tuple, Optional, List

import numpy as np
import pandas as pd

# 親ディレクトリのモジュールを読み込むため .. を使用
from ..fitsio import Scantable

def generate_simulated_grid_data(
    center_coord: Tuple[float, float] = (75.0, 30.0), # (RA, DEC) in deg
    grid_n: Tuple[int, int] = (5, 5),                 # (nx, ny)
    grid_step: Tuple[float, float] = (30.0, 30.0),    # (dx, dy) in arcsec
    projection: Optional[str] = "TAN",              # WCS projection: TAN/SIN/GNOMONIC, GLS/SFL/SINE, else CAR
    pos_jitter: float = 5.0,                          # 座標のランダムズレ (arcsec)
    axis_type: str = "velo",                          # "velo" or "freq" (今回はheader生成ロジックでFREQ優先に修正)
    rest_freq: float = 115.2712018e9,                 # 12CO(1-0) [Hz]
    vel_range: Tuple[float, float] = (-50.0, 60.0),   # 速度範囲 [km/s]
    vel_res: float = 0.5,                             # 速度分解能 [km/s]
    line_center: float = 10.0,                        # 輝線中心速度 [km/s]
    line_fwhm: float = 5.0,                           # 輝線幅 [km/s]
    peak_intensity: float = 10.0,                     # ピーク強度 [K]
    edge_factor: float = 0.3,                         # マップ端での強度比
    exposure_nominal: float = 30.0,                   # 標準積分時間 [sec]
    exposure_sigma: float = 5.0,                      # 積分時間のばらつき
    rms_at_nominal: float = 0.5,                      # RMS [K]
    baseline_windows: Optional[List[str]] = None,
    random_seed: Optional[int] = None
) -> Scantable:
    """
    Profile Map / MontageViewer テスト用のシミュレーションデータを生成する。
    Standardizerが誤動作しないよう、FITSヘッダー(WCS)を厳密に定義する。
    projection に応じて dRA の cos 補正を grid.py と同じ規則で切り替える。
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    # --- 1. 軸の生成 (Master Grid) ---
    c_kms = 299792.458
    n_chan = int(np.ceil((vel_range[1] - vel_range[0]) / vel_res))
    
    # 速度軸 (LSRK)
    v_axis = np.linspace(vel_range[0], vel_range[1], n_chan)
    
    # 周波数軸 (LSRK) の計算
    # f = f_rest * (1 - v / c)
    # df = -f_rest * dv / c
    f_axis = rest_freq * (1.0 - v_axis / c_kms)
    
    # チャンネル増分 (Hz) ※通常は負の値になる（速度が増えると周波数は減るため）
    df_hz = f_axis[1] - f_axis[0] 

    # --- 2. 空間グリッド計算 (省略: 元コードと同じ) ---
    nx, ny = grid_n
    dx_sec, dy_sec = grid_step
    ra0, dec0 = center_coord
    
    x_idx = np.arange(nx) - (nx - 1) / 2.0
    y_idx = np.arange(ny) - (ny - 1) / 2.0
    xx, yy = np.meshgrid(x_idx, y_idx)
    xx, yy = xx.flatten(), yy.flatten()

    # --- Projection handling (match grid.py logic) ---
    proj_mode = str(projection).upper() if projection is not None else "CAR"
    
    n_points = len(xx)
    ra_list, dec_list = [], []
    x_offsets, y_offsets = [], []

    for i in range(n_points):
        # 座標の揺らぎ
        off_x = xx[i] * dx_sec + np.random.normal(0, pos_jitter)
        off_y = yy[i] * dy_sec + np.random.normal(0, pos_jitter)

        # Dec/Lat offset is always linear in these small-angle approximations
        d_dec = off_y / 3600.0
        lat = dec0 + d_dec
        lat0 = dec0

        # cos-factor depends on projection
        if proj_mode in ("TAN", "SIN", "GNOMONIC"):
            # (2) Reference Point Dec: Orthogonal grid (Keeps square pixels)
            cosf = float(np.cos(np.deg2rad(lat0)))
        elif proj_mode in ("GLS", "SFL", "SINE"):
            # (1) Local Dec: Equal Area (GLS), standard for large single-dish maps
            cosf = float(np.cos(np.deg2rad(lat)))
        else:
            # CAR / NONE / raw
            cosf = 1.0

        # guard against division by (near-)zero at poles
        if abs(cosf) < 1e-12:
            cosf = 1e-12 if cosf >= 0 else -1e-12

        d_ra = -off_x / (3600.0 * cosf)
        
        ra_list.append(ra0 + d_ra)
        dec_list.append(dec0 + d_dec)
        x_offsets.append(off_x)
        y_offsets.append(off_y)

    # --- 3. 空間強度分布 (Gaussian) ---
    max_dist = ((max(nx, ny) - 1) / 2.0) * max(dx_sec, dy_sec)
    if max_dist > 0 and 0 < edge_factor < 1:
        spatial_sigma2 = -(max_dist**2) / (2.0 * np.log(edge_factor))
    else:
        spatial_sigma2 = np.inf

    # --- 4. データ生成ループ (中身は元コードとほぼ同じ) ---
    data_list = []
    rows = []
    line_sigma = line_fwhm / 2.35482

    def parse_win_str(ws):
        ranges = []
        for p in ws:
            try:
                s, e = map(float, p.split(":"))
                ranges.append((s, e))
            except: pass
        return ranges

    # --- 4a. Baseline fitting windows (BSL_WINF) ---
    # If baseline_windows is not provided, choose two velocity ranges that bracket the line,
    # avoiding the emission region. BSL_WINF is always stored in velocity [km/s].
    if baseline_windows is None:
        vmin, vmax = float(vel_range[0]), float(vel_range[1])
        span = vmax - vmin

        # Exclude region around the line: ~4*FWHM (at least 10 km/s)
        guard = max(4.0 * float(line_fwhm), 10.0)
        l0, l1 = vmin, float(line_center) - guard
        r0, r1 = float(line_center) + guard, vmax

        # Ensure each window has at least ~10 km/s (or >= 8 channels) width
        minw = max(10.0, 8.0 * float(vel_res))
        if (l1 - l0) < minw:
            l1 = vmin + minw
        if (r1 - r0) < minw:
            r0 = vmax - minw

        # Final guard against overlap / out-of-range
        l1 = min(l1, float(line_center) - max(2.0 * float(line_fwhm), 5.0))
        r0 = max(r0, float(line_center) + max(2.0 * float(line_fwhm), 5.0))
        l1 = max(l0, l1)
        r0 = min(r1, r0)

        # If still invalid, fall back to edge windows
        if not (l1 > l0 and r1 > r0 and l1 < r0):
            edge = max(minw, 0.25 * span)
            l0, l1 = vmin, vmin + edge
            r0, r1 = vmax - edge, vmax

        baseline_windows = [f"{l0:.1f}:{l1:.1f}", f"{r0:.1f}:{r1:.1f}"]

    base_ranges = parse_win_str(baseline_windows)
    for i in range(n_points):
        # A. 強度
        r2 = x_offsets[i]**2 + y_offsets[i]**2
        intensity = peak_intensity * np.exp(-r2 / (2 * spatial_sigma2)) if np.isfinite(spatial_sigma2) else peak_intensity
        
        exposure = max(1.0, np.random.normal(exposure_nominal, exposure_sigma))
        rms = rms_at_nominal * np.sqrt(exposure_nominal / exposure)
        
        local_v_center = line_center + np.random.uniform(-0.5, 0.5)
        signal = intensity * np.exp(-(v_axis - local_v_center)**2 / (2 * line_sigma**2))
        
        x_norm = np.linspace(-1, 1, n_chan)
        ripple_amp = rms * np.random.uniform(0.5, 2.0)
        ripple_freq = np.random.uniform(1.0, 3.0)
        ripple_phase = np.random.uniform(0, 2*np.pi)
        
        baseline_true = ripple_amp * np.sin(ripple_freq * np.pi * x_norm + ripple_phase) + \
                        np.random.uniform(-rms, rms) * x_norm + \
                        np.random.normal(0, rms/2)
        
        noise = np.random.normal(0, rms, n_chan)
        raw_spectrum = signal + baseline_true + noise
        
        # B. ベースライン除去 (Poly)
        mask_fit = np.zeros(n_chan, dtype=bool)
        current_wins_str = []
        for (s, e) in base_ranges:
            s_jit = s + np.random.uniform(-2.0, 2.0)
            e_jit = e + np.random.uniform(-2.0, 2.0)
            if s_jit > e_jit: s_jit, e_jit = e_jit, s_jit
            mask_fit |= (v_axis >= s_jit) & (v_axis <= e_jit)
            current_wins_str.append(f"{s_jit:.1f}:{e_jit:.1f}")
            
        poly_order = np.random.choice([1, 2, 3])
        if np.count_nonzero(mask_fit) > poly_order + 1:
            coeffs = np.polyfit(v_axis[mask_fit], raw_spectrum[mask_fit], deg=poly_order)
            fitted_baseline = np.polyval(coeffs, v_axis)
            bsl_rms_calc = np.std(raw_spectrum[mask_fit] - fitted_baseline[mask_fit])
        else:
            coeffs = np.zeros(poly_order + 1)
            fitted_baseline = np.zeros_like(v_axis)
            bsl_rms_calc = rms
            
        final_data = raw_spectrum - fitted_baseline
        data_list.append(final_data.astype(np.float32))
        
        obs_time = datetime.datetime.utcnow() + datetime.timedelta(seconds=i*exposure_nominal)
        
        rows.append({
            "SCAN": 1,
            "RA": ra_list[i],
            "DEC": dec_list[i],
            "timestamp": obs_time,
            "EXPOSURE": exposure,
            "TSYS": 150.0 + np.random.normal(0, 10.0),
            "BSL_DONE": True,
            "BSL_POLY": int(poly_order),
            "BSL_WINF": ";".join(current_wins_str),
            "BSL_COEF": coeffs.tolist() if hasattr(coeffs, "tolist") else list(coeffs),
            "BSL_RMS": float(bsl_rms_calc),
            "v_corr_kms": 0.0
        })

    data_arr = np.vstack(data_list)
    table_df = pd.DataFrame(rows)
    table_df["_local_idx"] = np.arange(n_points)

    # --- 5. メタデータ (Header) の修正 [重要] ---
    # Standardizerが誤読しないよう、FREQ基準で厳密に定義する
    header = {
        "TELESCOP": "SIMULATOR",
        "OBJECT":   "SimGrid_Corrected",
        "EQUINOX":  2000.0,
        "SPECSYS":  "LSRK",
        "RESTFREQ": float(rest_freq),
        "RESTFRQ":  float(rest_freq), # 互換性のため両方入れる
        "BUNIT":    "K",
        "NCHAN":    int(n_chan),
        
        # --- WCS 定義 (FREQ基準) ---
        "CTYPE1":   "FREQ-LSR",      # Standardizerがこれを認識すればHz計算を行う
        "CUNIT1":   "Hz",
        "CRPIX1":   1.0,             # 1-based index の先頭
        "CRVAL1":   float(f_axis[0]), # 1ピクセル目の周波数 (Hz)
        "CDELT1":   float(df_hz),     # 1チャンネルあたりの周波数増分 (Hz)
        
        # --- 補助情報 (速度変換用) ---
        # 多くのFITSリーダーはこれを見て速度軸を自動計算する
        "VELREF":   257,              # Radio LSRK
        "ALTRVAL":  float(v_axis[0] * 1000.0), # 参考: 始点速度 (m/s)
        "ALTRPIX":  1.0,
        "DELTAV":   float(vel_res * 1000.0)    # 参考: 速度分解能 (m/s)
    }
    
    # ユーザーが明示的に "velo" を指定しても、基本はFREQで記述し、
    # ソフトウェア側で変換させるのがFITSの流儀として安全
    
    history = {"task": "generate_simulated_grid_data_corrected"}
    return Scantable(meta=header, data=data_arr, table=table_df, history=history)
