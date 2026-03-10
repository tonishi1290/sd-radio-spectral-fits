# src/sd_radio_spectral_fits/otf/core.py
from __future__ import annotations

import builtins
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import j1

from .config import MapConfig, GridInput, GridResult

# ==========================================
# 1. Validation / Normalization Helpers
# ==========================================

def _ensure_1d(name: str, arr: np.ndarray, ndump: int) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim != 1 or a.shape[0] != ndump:
        raise ValueError(f"{name} must have shape ({ndump},), got {a.shape}")
    return a


def _validate_input(input_data: GridInput) -> None:
    x = np.asarray(input_data.x)
    y = np.asarray(input_data.y)
    spec = np.asarray(input_data.spec)
    flag = np.asarray(input_data.flag)
    time = np.asarray(input_data.time)

    if spec.ndim != 2:
        raise ValueError(f"spec must be 2D (ndump, nchan), got shape={spec.shape}")
    ndump = spec.shape[0]

    _ensure_1d("x", x, ndump)
    _ensure_1d("y", y, ndump)
    _ensure_1d("flag", flag, ndump)
    _ensure_1d("time", time, ndump)

    if input_data.rms is not None:
        _ensure_1d("rms", np.asarray(input_data.rms), ndump)
    if input_data.tint is not None:
        _ensure_1d("tint", np.asarray(input_data.tint), ndump)
    if input_data.tsys is not None:
        _ensure_1d("tsys", np.asarray(input_data.tsys), ndump)
    if input_data.is_turnaround is not None:
        _ensure_1d("is_turnaround", np.asarray(input_data.is_turnaround), ndump)


def _normalize_backend_and_kernel(config: GridConfig) -> tuple[str, str]:
    backend = str(config.backend)
    kernel = str(config.kernel).lower()

    # バックエンドとカーネルのエイリアス解決
    if backend == "numpy_gjinc":
        backend = "numpy"
        if kernel != "gjinc":
            raise ValueError("backend='numpy_gjinc' requires kernel='gjinc'")
    elif backend == "numpy_gauss":
        backend = "numpy"
        if kernel != "gauss":
            raise ValueError("backend='numpy_gauss' requires kernel='gauss'")
    elif backend == "cygrid_gauss":
        backend = "cygrid"
        if kernel != "gauss":
            raise ValueError("backend='cygrid_gauss' requires kernel='gauss'")

    if kernel not in ("gjinc", "gauss"):
        raise ValueError(f"Unknown kernel: {config.kernel!r}. Use 'gjinc' or 'gauss'.")

    return backend, kernel


def _validate_and_resolve_config(config: GridConfig) -> dict[str, float | str]:
    backend_impl, kernel_name = _normalize_backend_and_kernel(config)

    if config.estimator == 'plane':
        raise NotImplementedError("estimator='plane' is not implemented yet. Please use 'avg'.")

    # 基本的な正の数チェック
    if config.nx <= 0 or config.ny <= 0:
        raise ValueError("nx and ny must be positive")
    if config.cell_arcsec <= 0 or config.beam_fwhm_arcsec <= 0:
        raise ValueError("cell_arcsec and beam_fwhm_arcsec must be positive")
    if config.chunk_ch <= 0:
        raise ValueError("chunk_ch must be positive")
    
    # 空間カーネルのパラメータ解決
    if config.gwidth_pix is not None:
        gwidth_pix = float(config.gwidth_pix)
    elif config.gwidth_beam is not None:
        gwidth_pix = float(config.gwidth_beam) * (config.beam_fwhm_arcsec / config.cell_arcsec)
    else:
        gwidth_pix = 2.10
        
    if gwidth_pix <= 0:
        raise ValueError("gwidth must be positive")

    # GJINC特有のパラメータ解決
    if kernel_name == "gjinc":
        if config.jwidth_pix is not None:
            jwidth_pix = float(config.jwidth_pix)
        elif config.jwidth_beam is not None:
            jwidth_pix = float(config.jwidth_beam) * (config.beam_fwhm_arcsec / config.cell_arcsec)
        else:
            jwidth_pix = 1.55
    else:
        jwidth_pix = np.nan

    # 打ち切り半径 (support radius) の解決
    if config.support_radius_pix is not None:
        support_radius_pix = float(config.support_radius_pix)
    else:
        if kernel_name == "gjinc" and config.truncate == "first_null":
            support_radius_pix = 1.21967 * jwidth_pix
        elif isinstance(config.truncate, (int, float)):
            support_radius_pix = float(config.truncate)
        else:
            raise ValueError("Invalid truncate mode.")

    return {
        "backend_impl": backend_impl,
        "kernel_name": kernel_name,
        "gwidth_pix": gwidth_pix,
        "jwidth_pix": jwidth_pix,
        "support_radius_pix": support_radius_pix,
    }

# ==========================================
# 2. Kernels
# ==========================================

def _kernel_gauss(r_arcsec: np.ndarray, gwidth_arcsec: float) -> np.ndarray:
    """Gaussian kernel with HWHM parameter gwidth_arcsec."""
    return np.exp(-np.log(2.0) * (r_arcsec / gwidth_arcsec) ** 2)


def _kernel_gjinc(r_arcsec: np.ndarray, gwidth_arcsec: float, jwidth_arcsec: float, eps_u0: float) -> np.ndarray:
    """GJINC kernel (CASA-like form) with safe handling at r=0."""
    u = np.pi * r_arcsec / jwidth_arcsec
    jinc = np.empty_like(u)
    
    # r=0 付近の発散を防ぐ安全処理
    m0 = np.abs(u) < eps_u0
    jinc[m0] = 0.5
    
    mu = ~m0
    if np.any(mu):
        jinc[mu] = j1(u[mu]) / u[mu]
        
    gauss = np.exp(-np.log(2.0) * (r_arcsec / gwidth_arcsec) ** 2)
    return jinc * gauss

# ==========================================
# 3. Main API & Implementation
# ==========================================

def grid_otf(input_data: GridInput, config: GridConfig) -> GridResult:
    """
    OTF Gridding コアエンジン。指定された空間カーネルと設定に基づいてグリッディングを実行する。
    """
    _validate_input(input_data)
    kres = _validate_and_resolve_config(config)

    return _backend_numpy_avg(input_data, config, kres)


def _backend_numpy_avg(input_data: GridInput, config: GridConfig, kernel_resolved: dict[str, float | str]) -> GridResult:
    """NumPyベースの加重平均グリッディング実装"""
    dtype_np = np.float64 if config.dtype == "float64" else np.float32

    x_all = np.asarray(input_data.x)
    y_all = np.asarray(input_data.y)
    spec_all = np.asarray(input_data.spec)
    flag_all = np.asarray(input_data.flag)
    time_all = np.asarray(input_data.time)
    ndump, nchan = spec_all.shape

    # --- 1. データの前処理・フィルタリング ---
    valid = np.asarray(flag_all > 0, dtype=bool)
    if config.exclude_turnaround and input_data.is_turnaround is not None:
        valid &= ~np.asarray(input_data.is_turnaround, dtype=bool)

    # 空間・スペクトルの有効性チェック
    valid &= np.isfinite(x_all) & np.isfinite(y_all)
    valid &= np.any(np.isfinite(spec_all), axis=1)

    rms_all, tint_all, tsys_all = None, None, None
    if input_data.rms is not None:
        rms_all = np.asarray(input_data.rms)
        valid &= np.isfinite(rms_all) & (rms_all > 0)
    if input_data.tint is not None:
        tint_all = np.asarray(input_data.tint)
        valid &= np.isfinite(tint_all) & (tint_all >= 0)
    if input_data.tsys is not None:
        tsys_all = np.asarray(input_data.tsys) # TSYSは一部欠損でも全データは落とさない

    if not np.any(valid):
        raise ValueError("No valid dumps remain after filtering.")

    # フィルタリングの適用
    x = x_all[valid].astype(dtype_np, copy=False)
    y = y_all[valid].astype(dtype_np, copy=False)
    spec = spec_all[valid].astype(dtype_np, copy=False)
    time_vals = time_all[valid].astype(dtype_np, copy=False)
    rms = rms_all[valid].astype(dtype_np, copy=False) if rms_all is not None else None
    tint = tint_all[valid].astype(dtype_np, copy=False) if tint_all is not None else None
    tsys = tsys_all[valid].astype(dtype_np, copy=False) if tsys_all is not None else None

    # 重み (q) の計算
    q = np.ones(len(x), dtype=dtype_np)
    if rms is not None:
        q *= (1.0 / (rms * rms + dtype_np(1e-12))) ** dtype_np(config.alpha_rms)
    if tint is not None and config.beta_tint > 0:
        q *= tint ** dtype_np(config.beta_tint)

    # 重みのクリップ（異常に重い1点が結果を支配するのを防ぐ）
    if config.weight_clip_quantile is not None and q.size > 0:
        q = np.minimum(q, dtype_np(np.quantile(q, config.weight_clip_quantile)))
    if config.weight_clip_max is not None:
        q = np.minimum(q, dtype_np(config.weight_clip_max))

    # --- 2. 幾何学パラメータとグリッドの準備 ---
    kernel_name = str(kernel_resolved["kernel_name"])
    backend_impl = str(kernel_resolved["backend_impl"])
    gwidth_arcsec = dtype_np(float(kernel_resolved["gwidth_pix"]) * config.cell_arcsec)
    jwidth_arcsec = dtype_np(float(kernel_resolved["jwidth_pix"]) * config.cell_arcsec) if kernel_name == "gjinc" else dtype_np(np.nan)
    r_trunc_arcsec = dtype_np(float(kernel_resolved["support_radius_pix"]) * config.cell_arcsec)

    # 出力グリッド座標 [arcsec] の生成
    xg_1d = (config.x0 + np.arange(config.nx, dtype=dtype_np) * dtype_np(config.cell_arcsec)).astype(dtype_np)
    yg_1d = (config.y0 + np.arange(config.ny, dtype=dtype_np) * dtype_np(config.cell_arcsec)).astype(dtype_np)
    xg2d, yg2d = np.meshgrid(xg_1d, yg_1d)
    grid_pts = np.column_stack((xg2d.ravel(), yg2d.ravel()))
    ngrid = grid_pts.shape[0]

    # KDTree で近傍探索
    tree = cKDTree(np.column_stack((x, y)))
    try:
        neighbors_list = tree.query_ball_point(grid_pts, r=float(r_trunc_arcsec), workers=-1)
    except TypeError:
        neighbors_list = tree.query_ball_point(grid_pts, r=float(r_trunc_arcsec))

    # --- 3. 出力配列の初期化 ---
    if config.fill_nan_for_invalid:
        cube_flat = np.full((ngrid, nchan), np.nan, dtype=dtype_np)
    else:
        cube_flat = np.zeros((ngrid, nchan), dtype=dtype_np)

    weight_flat = np.zeros(ngrid, dtype=dtype_np)
    hit_flat = np.zeros(ngrid, dtype=np.int32)
    mask_flat = np.zeros(ngrid, dtype=bool)

    # 各種診断マップ
    xeff_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_diag_maps else None
    yeff_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_diag_maps else None
    dr_eff_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_diag_maps else None
    neff_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_neff_map else None
    rms_flat = np.full(ngrid, np.nan, dtype=dtype_np) if (config.emit_rms_map and rms is not None) else None
    
    # 新規追加の診断マップ (TIME, TINT, TSYS)
    time_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_time_map else None
    tint_flat = np.full(ngrid, np.nan, dtype=dtype_np) if (config.emit_tint_map and tint is not None) else None
    tsys_flat = np.full(ngrid, np.nan, dtype=dtype_np) if (config.emit_tsys_map and tsys is not None) else None

    # --- 4. メインループ (グリッド点ごと) ---
    for g_idx, nbrs_list in enumerate(neighbors_list):
        if len(nbrs_list) < config.n_min_avg:
            continue
            
        nbrs = np.asarray(nbrs_list, dtype=np.int64)

        dx = x[nbrs] - grid_pts[g_idx, 0]
        dy = y[nbrs] - grid_pts[g_idx, 1]
        r = np.sqrt(dx * dx + dy * dy)

        # 厳密なサポート半径の判定
        m_keep = (r <= r_trunc_arcsec)
        if not np.any(m_keep):
            continue
        if not np.all(m_keep):
            nbrs = nbrs[m_keep]
            r = r[m_keep]
            
        if nbrs.size < config.n_min_avg:
            continue

        # 空間カーネルの計算
        if kernel_name == "gjinc":
            k = _kernel_gjinc(r, float(gwidth_arcsec), float(jwidth_arcsec), float(config.eps_u0)).astype(dtype_np, copy=False)
        else:
            k = _kernel_gauss(r, float(gwidth_arcsec)).astype(dtype_np, copy=False)

        # 重みの総合 (カーネル x 品質重み)
        w = k * q[nbrs]
        m_w = np.isfinite(w) & (w > 0)
        
        if not np.any(m_w):
            continue
        if not np.all(m_w):
            nbrs = nbrs[m_w]
            w = w[m_w]
            
        if nbrs.size < config.n_min_avg:
            continue

        W_sum = np.sum(w, dtype=dtype_np)
        if (not np.isfinite(W_sum)) or (W_sum <= dtype_np(config.eps_weight_sum)):
            continue

        # 基本マップへの記録
        weight_flat[g_idx] = W_sum
        hit_flat[g_idx] = int(np.count_nonzero(w > 0))
        mask_flat[g_idx] = bool((hit_flat[g_idx] >= config.n_min_avg) and (W_sum > dtype_np(config.eps_weight_sum)))
        
        if not mask_flat[g_idx]:
            continue

        # 診断マップの計算
        if config.emit_diag_maps:
            xeff = np.sum(w * x[nbrs], dtype=dtype_np) / W_sum
            yeff = np.sum(w * y[nbrs], dtype=dtype_np) / W_sum
            xeff_flat[g_idx] = xeff
            yeff_flat[g_idx] = yeff
            dr_eff_flat[g_idx] = np.sqrt((xeff - grid_pts[g_idx, 0]) ** 2 + (yeff - grid_pts[g_idx, 1]) ** 2) / dtype_np(config.cell_arcsec)

        if config.emit_neff_map:
            w2_sum = np.sum(w * w, dtype=dtype_np)
            if w2_sum > 0:
                neff_flat[g_idx] = (W_sum * W_sum) / w2_sum

        if rms_flat is not None:
            rms_flat[g_idx] = np.sqrt(np.sum((w * w) * (rms[nbrs] * rms[nbrs]), dtype=dtype_np) / (W_sum * W_sum))
            
        # --- 新規マップの計算 (TIME, TINT, TSYS) ---
        if time_flat is not None:
            time_flat[g_idx] = np.sum(w * time_vals[nbrs], dtype=dtype_np) / W_sum
            
        # TINTの計算条件から neff_flat 依存を外す
        if tint_flat is not None:
            if neff_flat is not None:
                neff_val = neff_flat[g_idx]
            else:
                w2_sum = np.sum(w * w, dtype=dtype_np)
                neff_val = (W_sum * W_sum) / w2_sum if w2_sum > 0 else 1.0
            tint_flat[g_idx] = (np.sum(w * tint[nbrs], dtype=dtype_np) / W_sum) * neff_val
            
        if tsys_flat is not None:
            w_tsys = w * tsys[nbrs]
            m_tsys = np.isfinite(w_tsys)
            if np.any(m_tsys):
                # NaNを落としてTSYSが存在するデータだけで平均をとる
                tsys_flat[g_idx] = np.sum(w_tsys[m_tsys], dtype=dtype_np) / np.sum(w[m_tsys], dtype=dtype_np)

        # スペクトル合成 (チャネルごとのNaN考慮とチャンク処理)
        spec_nbrs = spec[nbrs, :]
        for ch0 in range(0, nchan, config.chunk_ch):
            ch1 = builtins.min(ch0 + config.chunk_ch, nchan)
            
            spec_chunk = spec_nbrs[:, ch0:ch1]
            valid_data = np.isfinite(spec_chunk)
            
            # (N_nbrs, N_chunk_channels) の重み行列
            w_2d = w[:, None] * valid_data
            
            W_sum_chan = np.sum(w_2d, axis=0)
            spec_sum_chan = np.sum(w_2d * np.where(valid_data, spec_chunk, 0.0), axis=0)
            
            # 重みの合計が閾値を超えているチャネルのみ平均値を記録
            m_chan = W_sum_chan > dtype_np(config.eps_weight_sum)
            cube_view = cube_flat[g_idx, ch0:ch1]
            cube_view[m_chan] = spec_sum_chan[m_chan] / W_sum_chan[m_chan]
            cube_flat[g_idx, ch0:ch1] = cube_view

    # --- 5. 2次元/3次元配列への変形と返却 ---
    ny, nx = config.ny, config.nx

    n_warn = None
    if dr_eff_flat is not None:
        n_warn = int(np.sum(np.isfinite(dr_eff_flat) & (dr_eff_flat > dtype_np(config.dr_eff_warn_pix))))

    # 実効ビームサイズ (FWHM) の理論計算: sqrt(Beam^2 + Kernel_FWHM^2)
    # ガウシアン/GJINCともに、gwidth は HWHM を意味するため、Kernel_FWHM = 2.0 * gwidth
    bmaj_eff_arcsec = float(np.sqrt(config.beam_fwhm_arcsec**2 + (2.0 * float(gwidth_arcsec))**2))

    return GridResult(
        cube=cube_flat.reshape(ny, nx, nchan),
        weight_map=weight_flat.reshape(ny, nx),
        hit_map=hit_flat.reshape(ny, nx),
        mask_map=mask_flat.reshape(ny, nx),
        xeff_map=xeff_flat.reshape(ny, nx) if xeff_flat is not None else None,
        yeff_map=yeff_flat.reshape(ny, nx) if yeff_flat is not None else None,
        dr_eff_map_pix=dr_eff_flat.reshape(ny, nx) if dr_eff_flat is not None else None,
        neff_map=neff_flat.reshape(ny, nx) if neff_flat is not None else None,
        rms_map=rms_flat.reshape(ny, nx) if rms_flat is not None else None,
        time_map=time_flat.reshape(ny, nx) if time_flat is not None else None,
        tint_map=tint_flat.reshape(ny, nx) if tint_flat is not None else None,
        tsys_map=tsys_flat.reshape(ny, nx) if tsys_flat is not None else None,
        meta={
            "kernel": kernel_name, 
            "backend": backend_impl,
            "cell_arcsec": float(config.cell_arcsec),
            "beam_fwhm_arcsec": float(config.beam_fwhm_arcsec),
            "bmaj_eff_arcsec": bmaj_eff_arcsec,
            "alpha_rms": float(config.alpha_rms),
            "beta_tint": float(config.beta_tint),
            "support_radius_pix": float(kernel_resolved["support_radius_pix"]),
            "dr_eff_warn_count": n_warn,
        },
    )
