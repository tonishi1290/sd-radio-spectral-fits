# src/sd_radio_spectral_fits/map/config.py
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Optional, Literal, Union


@dataclass
class MapConfig:
    """
    画像化(Gridding)およびマップ解析の統合設定クラス。
    旧 GridConfig の全機能を包含し、解析用パラメータを追加。
    """

    # --- 1. 空間グリッド定義 (Geometry) ---
    x0: float                # マップ中心からのオフセット [arcsec]
    y0: float
    nx: int                  # ピクセル数
    ny: int
    cell_arcsec: float       # ピクセルサイズ [arcsec]
    beam_fwhm_arcsec: float  # 望遠鏡のビームサイズ [arcsec]

    # --- 2. グリッディング・カーネル設定 ---
    kernel: Literal['gjinc', 'gauss'] = 'gjinc'
    gwidth_pix: Optional[float] = 2.10
    gwidth_beam: Optional[float] = None
    jwidth_pix: Optional[float] = None
    jwidth_beam: Optional[float] = None
    truncate: Union[Literal['first_null'], float] = 'first_null'
    support_radius_pix: Optional[float] = None
    chunk_ch: int = 256      # メモリ節約のためのチャンネル分割処理（デフォルト256）
    dtype: str = "float32"

    # --- 3. 重み付け・品質管理 ---
    alpha_rms: float = 0.5
    beta_tint: float = 0.0
    weight_clip_quantile: Optional[float] = 0.95
    weight_clip_max: Optional[float] = None
    exclude_turnaround: bool = True

    # --- 4. 推定器 / QC / 安全係数 ---
    estimator: Literal['avg', 'plane'] = 'avg'
    n_min_avg: int = 2
    n_min_plane: int = 6
    cond_max: float = 1e6
    dr_eff_warn_pix: float = 0.2
    eps_u0: float = 0.01
    eps_weight_sum: float = 1e-8

    # --- 5. 出力マップ制御 (core.py が要求するフラグ群) ---
    fill_nan_for_invalid: bool = True
    emit_diag_maps: bool = True
    emit_neff_map: bool = True
    emit_rms_map: bool = True
    emit_time_map: bool = True
    emit_tint_map: bool = True
    emit_tsys_map: bool = True

    # =========================================================
    # --- 6. 解析・出力設定 (Step 3以降で使用) ---
    # =========================================================
    generate_mask: bool = False      # 3DマスクとMoment0を自動生成するか
    mask_method: str = "smooth_mask" # 'simple', 'smooth_mask', 'derivative'

    # 共通 / simple用
    mask_sigma: float = 3.0          # マスク生成のベース閾値

    # smooth_mask用
    mask_high_snr: float = 3.0       # コア抽出の閾値
    mask_low_snr: float = 1.5        # 拡張マスクの閾値
    mask_min_vol: int = 27           # 許容する最小体積(ピクセル数)

    # derivative用
    mask_sigma_v: float = 2.0        # 速度方向のフィルタ幅(チャンネル数)
    mask_deriv_snr: float = 3.0      # 2次微分の検出閾値
    mask_dilation: int = 2           # 速度方向へのマスク拡張幅

    # 運用・出力用
    mask_compression: Optional[str] = 'PLIO_1' # FITS保存時の圧縮形式 (互換性重視ならNoneも可)
    dv_kms: Optional[float] = None   # 速度リグリッドの分解能 [km/s]
    output_prefix: str = "map"      # 出力ファイルの接頭辞

    # --- 7. 実行バックエンド ---
    backend: str = 'numpy'           # 'numpy' or 'numba' (将来用)

    # --- 8. 実行時オプション ---
    verbose: bool = False
    workers: int = -1                # cKDTree query_ball_point 用。旧挙動は -1。
    sort_neighbors: bool = False     # 旧挙動維持のため既定は False。
    reproducible_mode: bool = False  # True のとき float64 / workers=1 / sort固定。


# 以前のコードとの互換性のためにエイリアスを残す
GridConfig = MapConfig


@dataclass
class GridInput:
    """入力データ構造"""
    x: np.ndarray
    y: np.ndarray
    spec: np.ndarray
    flag: np.ndarray
    time: np.ndarray
    rms: Optional[np.ndarray] = None
    tint: Optional[np.ndarray] = None
    tsys: Optional[np.ndarray] = None
    scan_id: Optional[np.ndarray] = None
    subscan_id: Optional[np.ndarray] = None
    scan_dir: Optional[np.ndarray] = None
    is_turnaround: Optional[np.ndarray] = None


@dataclass
class GridResult:
    """出力データ構造"""
    cube: np.ndarray
    weight_map: np.ndarray
    hit_map: np.ndarray
    mask_map: np.ndarray
    # 診断用マップ群
    xeff_map: Optional[np.ndarray] = None
    yeff_map: Optional[np.ndarray] = None
    dr_eff_map_pix: Optional[np.ndarray] = None
    neff_map: Optional[np.ndarray] = None
    rms_map: Optional[np.ndarray] = None
    time_map: Optional[np.ndarray] = None
    tint_map: Optional[np.ndarray] = None
    tsys_map: Optional[np.ndarray] = None
    meta: Optional[dict] = None

    # ==========================================
    # --- Step 3 解析結果格納用 ---
    # ==========================================
    mask_3d: Optional[np.ndarray] = None
    mom0_map: Optional[np.ndarray] = None
