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

    Notes
    -----
    公開 kernel preset は pixel-based で解釈する。
    既定 kernel は ``kernel='sf'``、既定 ``convsupport=3`` とする。

    - ``cell_arcsec`` は明示指定を推奨する。
      ``cell_arcsec=None`` のときは ``beam_fwhm_arcsec / 3`` を自動採用する。
    - ``kernel='gjinc'`` / ``kernel='gjinc_beam'`` では ``kernel_preset='mangum'`` または
      ``kernel_preset='casa'`` を使う。互換 alias として
      ``mangum2007 -> mangum``、``legacy -> casa`` を受け付ける。
    - ``kernel='sf'`` では support は ``convsupport`` だけで制御する。
      ``support_radius_*`` や ``truncate`` は使わない。
    - 明示的な ``*_pix`` / ``*_beam`` / ``support_radius_*`` は、
      後方互換のため残す。
    """

    # --- 1. 空間グリッド定義 (Geometry) ---
    x0: float                # マップ中心からのオフセット [arcsec]
    y0: float
    nx: int                  # ピクセル数
    ny: int
    cell_arcsec: Optional[float] = None       # ピクセルサイズ [arcsec]; None のとき beam/3
    beam_fwhm_arcsec: Optional[float] = None  # 望遠鏡のビームサイズ [arcsec]

    # --- 2. グリッディング・カーネル設定 ---
    kernel: Literal['sf', 'sf_beam', 'gjinc', 'gjinc_beam', 'gauss'] = 'sf'
    kernel_preset: Optional[str] = None
    kernel_sign: Literal['auto', 'signed', 'positive_only'] = 'auto'
    convsupport: int = 3

    # sf_beam は連続 beam-aware spheroidal kernel の最小実装。
    # 基本は ``sf_beam_match_convsupport`` により design point (cell=beam/3) で
    # official ``SF(n)`` に合わせた preset を使う。
    # 個別に指定したいときは ``sf_beam_support_beam`` と
    # ``sf_beam_shape_c`` を明示する。
    sf_beam_match_convsupport: Optional[int] = 3
    sf_beam_support_beam: Optional[float] = None
    sf_beam_shape_c: Optional[float] = None

    # gjinc_beam は連続 beam-aware GJINC kernel。
    # 基本は ``kernel_preset='mangum'`` または ``'casa'`` により
    # design point (cell=beam/3) の public pixel-based gjinc default を
    # beam 単位へ写した preset を使う。
    # 必要なら beam 単位の幅・support を個別指定できる。
    gjinc_beam_gwidth_beam: Optional[float] = None
    gjinc_beam_jwidth_beam: Optional[float] = None
    gjinc_beam_support_beam: Optional[float] = None

    # kernel_sign は GJINC/GJINC_BEAM/GAUSS の負重みをどう扱うかを表す。
    #   auto          : mangum -> signed, casa -> positive_only, sf/sf_beam/gauss -> positive_only
    #   signed        : 有限な符号付き重みをそのまま使う
    #   positive_only : 正の重みだけを使う
    # 明示指定がある場合は *_pix > *_beam > preset default の順で解決する
    gwidth_pix: Optional[float] = None
    gwidth_beam: Optional[float] = None
    jwidth_pix: Optional[float] = None
    jwidth_beam: Optional[float] = None

    # support は kernel をどの半径で打ち切るかを表す cutoff radius。
    # sf では support = convsupport * cell。
    # gjinc/gauss では pix / beam / truncate の順で解決する。
    # gjinc_beam では dedicated な *_beam パラメータか preset を使う。
    # truncate=None のときは preset default を用いる。
    #   sf               : support = convsupport * cell
    #   mangum + gjinc   : support = 3 * cell
    #   casa   + gjinc   : support = first null of resolved jinc width
    #   gauss            : support = 3 * gwidth (3 * HWHM)
    truncate: Optional[Union[Literal['first_null'], float]] = None
    support_radius_pix: Optional[float] = None
    support_radius_beam: Optional[float] = None

    chunk_ch: int = 256      # メモリ節約のためのチャンネル分割処理（デフォルト256）
    dtype: str = "float32"

    # --- 3. 重み付け・品質管理 ---
    # weight_mode は dump ごとの追加重み q_i の作り方を決める。
    #   uniform : q_i = 1
    #   rms     : q_i *= (1 / rms_i^2)^alpha_rms
    #             ここで rms_i は BSL_RMS 列から与える。
    #             weight_mode='rms' のとき BSL_RMS 欠落・NaN・非正値は error。
    # beta_tint は RMS に積分時間の効果がすでに入っている前提で既定 0.0 とする。
    # weight_clip_quantile は既定では無効化し、必要時のみ明示指定で使う。
    weight_mode: Literal['uniform', 'rms'] = 'uniform'
    alpha_rms: float = 1.0
    beta_tint: float = 0.0
    weight_clip_quantile: Optional[float] = None
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

    # 追加の post-gridding 妥当性判定
    # min_abs_weight_ratio > 0 のとき、|WEIGHT| / median(|WEIGHT|) が
    # この値未満の画素を無効化する。CASA minweight と同趣旨。
    min_abs_weight_ratio: float = 0.0
    # min_cancel_ratio > 0 のとき、|WSUM| / WABS がこの値未満の
    # 画素を無効化する。signed kernel の強い相殺を検出する。
    min_cancel_ratio: float = 0.0

    # beam / sampling に関する補助
    warn_if_cell_coarse: bool = True       # cell > beam/3 のとき warning
    estimate_effective_beam: bool = True   # nominal / empirical beam を推定

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
    vmin_kms: Optional[float] = None # 出力速度軸の下限 [km/s]
    vmax_kms: Optional[float] = None # 出力速度軸の上限 [km/s]
    output_prefix: str = "map"      # 出力ファイルの接頭辞

    # --- 7. 実行バックエンド ---
    backend: str = 'numpy'           # 'numpy' or 'numba' (将来用)

    # --- 8. 実行時オプション ---
    verbose: bool = False
    workers: int = -1                # cKDTree query_ball_point 用。旧挙動は -1。
    sort_neighbors: bool = False     # 旧挙動維持のため既定は False。
    reproducible_mode: bool = False  # True のとき float64 / workers=1 / sort固定。
    write_diagnostics: bool = False  # True のとき OTF 診断 sidecar を出力。
    diagnostics_prefix: Optional[str] = None  # 診断ファイルの接頭辞。

    def __post_init__(self) -> None:
        if isinstance(self.beam_fwhm_arcsec, (bool, np.bool_)) or self.beam_fwhm_arcsec is None:
            raise ValueError("beam_fwhm_arcsec must be specified and positive")
        beam = float(self.beam_fwhm_arcsec)
        if not np.isfinite(beam) or beam <= 0:
            raise ValueError("beam_fwhm_arcsec must be positive")
        self.beam_fwhm_arcsec = beam

        if isinstance(self.cell_arcsec, (bool, np.bool_)):
            raise ValueError("cell_arcsec must be a positive float or None")
        if self.cell_arcsec is None:
            self.cell_arcsec = beam / 3.0
        else:
            cell = float(self.cell_arcsec)
            if not np.isfinite(cell) or cell <= 0:
                raise ValueError("cell_arcsec must be positive")
            self.cell_arcsec = cell


# 以前のコードとの互換性のためにエイリアスを残す
GridConfig = MapConfig


def normalize_row_flag_mask(values, ndump: Optional[int] = None, *, allow_none: bool = False, name: str = "flag") -> np.ndarray:
    """
    Normalize a row-level dump selection mask to a strict 1D bool array.

    Public semantics:
      - shape = (ndump,)
      - True  = use this dump row
      - False = drop this dump row

    Compatibility:
      - bool arrays are accepted directly.
      - integer / float arrays containing only 0 and 1 are accepted and converted
        to bool. Other numeric values are rejected to avoid the old ambiguous
        ``flag > 0`` semantics.
    """
    if values is None:
        if allow_none:
            if ndump is None:
                raise ValueError(f"{name} requires ndump when allow_none=True")
            return np.ones(int(ndump), dtype=bool)
        raise ValueError(f"{name} must be provided as a 1D bool row mask")

    arr = np.asarray(values)
    if arr.ndim != 1:
        raise ValueError(f"{name} must have shape ({ndump},) and be a 1D row mask; got {arr.shape}")
    if ndump is not None and arr.shape[0] != int(ndump):
        raise ValueError(f"{name} must have shape ({int(ndump)},), got {arr.shape}")

    if arr.dtype.kind == 'b':
        return arr.astype(bool, copy=False)

    if arr.dtype.kind in 'iu':
        bad = (arr != 0) & (arr != 1)
        if np.any(bad):
            unique_bad = np.unique(arr[bad])
            raise ValueError(
                f"{name} must be a 1D bool row mask. Integer compatibility accepts only 0/1; "
                f"found invalid values {unique_bad[:8]!r}."
            )
        return arr.astype(bool, copy=False)

    if arr.dtype.kind == 'f':
        finite = np.isfinite(arr)
        bad = (~finite) | ((arr != 0.0) & (arr != 1.0))
        if np.any(bad):
            bad_vals = arr[bad]
            preview = bad_vals[:8]
            raise ValueError(
                f"{name} must be a 1D bool row mask. Float compatibility accepts only 0.0/1.0; "
                f"found invalid values {preview!r}."
            )
        return arr.astype(bool)

    # Avoid silently treating arbitrary objects / strings as bools.
    raise ValueError(
        f"{name} must be a 1D bool row mask with shape ({int(ndump) if ndump is not None else 'ndump'},). "
        f"Got dtype={arr.dtype!r}."
    )


@dataclass
class GridInput:
    """入力データ構造

    Notes
    -----
    ``flag`` is a row-level dump mask with shape ``(ndump,)``.

    - ``True``  : use this dump row
    - ``False`` : drop this dump row

    Channel-level invalid data should be represented by ``NaN`` in ``spec``
    rather than by a 2D flag array.
    """
    x: np.ndarray
    y: np.ndarray
    spec: np.ndarray
    flag: Optional[np.ndarray]
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
    hit_map: np.ndarray   # 実際に gridding 分母に使ったサンプル数（sign filter 後）
    mask_map: np.ndarray  # 空間的な分母有効性（channel ごとの NaN マスク前）
    nsamp_map: Optional[np.ndarray] = None  # support 内の有限サンプル数（sign filter 前）
    wsum_map: Optional[np.ndarray] = None
    wabs_map: Optional[np.ndarray] = None
    cancel_map: Optional[np.ndarray] = None  # |WSUM| / WABS
    weight_rel_map: Optional[np.ndarray] = None  # |WEIGHT| / median(|WEIGHT|)
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
