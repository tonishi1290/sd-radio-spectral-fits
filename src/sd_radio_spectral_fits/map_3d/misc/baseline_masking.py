import warnings
import numpy as np
import scipy.ndimage as nd
import concurrent.futures
import multiprocessing

# =====================================================================
# 1. 超高速ベクトル化コアエンジン (内部関数)
# =====================================================================
def _core_baseline_engine(data_chunk, degree, max_iter, threshold, dilation):
    """
    1Dに平坦化されたチャンク (V, N) のベースラインをベクトル演算で特定するエンジン
    """
    V, N_chunk = data_chunk.shape
    
    # 計算精度と省メモリのバランスをとるため float32 で計算
    # 欠損値(NaN)は一旦0.0で埋める
    Y = np.nan_to_num(data_chunk, nan=0.0).astype(np.float32)
    original_valid_mask = np.isfinite(data_chunk)
    W = original_valid_mask.astype(np.float32)
    
    # ヴァンデルモンド行列による設計行列Xの生成（丸め誤差防止のため-1から1に規格化）
    v_idx = np.linspace(-1, 1, V, dtype=np.float32)
    X = np.vander(v_idx, degree + 1, increasing=True)
    reg = np.eye(degree + 1, dtype=np.float32) * 1e-7 
    
    for _ in range(max_iter):
        # np.einsumを用いたループなしの正規方程式の一括生成
        XTWX = np.einsum('vi,vn,vj->nij', X, W, X) + reg
        XTWY = np.einsum('vi,vn,vn->ni', X, W, Y)
        
        try:
            # 全空間ピクセルの係数を一気に解く
            beta = np.linalg.solve(XTWX, XTWY)
        except np.linalg.LinAlgError:
            # 特異行列エラーの回避（全データ欠損ピクセルなど）
            beta = np.zeros((N_chunk, degree + 1), dtype=np.float32)
            
        # ベースラインを復元し残差を計算
        baseline = np.einsum('vi,ni->vn', X, beta)
        residuals = Y - baseline
        
        # 有効領域(W>0)のみを残差として抽出
        masked_res = np.where(W > 0, residuals, np.nan)
        
        # 全NaNピクセルに対する無害な RuntimeWarning をミュート
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            median_res = np.nanmedian(masked_res, axis=0)
            mad = np.nanmedian(np.abs(masked_res - median_res), axis=0)
            
        # 欠損ピクセルでの演算結果(NaN)を0.0に戻し、後段の比較演算を安全にする
        median_res = np.nan_to_num(median_res, nan=0.0)
        mad = np.nan_to_num(mad, nan=0.0)
        
        # ノイズレベル推定と正の方向（輝線）に対するシグマクリッピング
        sigma = 1.4826 * np.where(mad == 0, 1e-10, mad)
        line_mask = (residuals - median_res) > (threshold * sigma)
        
        # np.rollを用いた超高速ビットシフトDilation（裾野の保護）
        if dilation > 0:
            dilated = line_mask.copy()
            for d in range(1, dilation + 1):
                dilated |= np.roll(line_mask, d, axis=0)
                dilated |= np.roll(line_mask, -d, axis=0)
            line_mask = dilated
            
        # 重みの更新 (輝線を0、ベースラインを1にする)
        W = np.where(line_mask, 0.0, 1.0).astype(np.float32)
        # もともとデータが無かった場所は強制的に無効(0)に戻す
        W *= original_valid_mask.astype(np.float32)
        
    # 最終的なline-free領域を True としたBoolean配列を返す
    return (W == 1.0)

# =====================================================================
# 2. 並列処理用ワーカー関数 (内部関数)
# =====================================================================
def _process_single_chunk_worker(sub_cube, y_crop_start, y_crop_end, x_crop_start, x_crop_end,
                                 degree, max_iter, threshold, dilation, spatial_sigma,
                                 y_dest_start, y_dest_end, x_dest_start, x_dest_end):
    """
    のりしろ込みのキューブを受け取り、スムージング、エンジン実行、クロップを行って返すワーカー
    """
    # NaN汚染を防ぐ安全な空間スムージング処理
    if spatial_sigma > 0:
        # NaNを0に置換してフィルター実行中のドミノ汚染を防止
        safe_sub_cube = np.nan_to_num(sub_cube, nan=0.0)
        sub_cube_proc = nd.gaussian_filter(safe_sub_cube, sigma=(0, spatial_sigma, spatial_sigma))
        
        # 元々NaNだった位置にNaNを復元し、無効領域として正しく認識させる
        sub_cube_proc = np.where(np.isfinite(sub_cube), sub_cube_proc, np.nan)
    else:
        sub_cube_proc = sub_cube
        
    V_sub, Y_sub, X_sub = sub_cube_proc.shape
    chunk_flat = sub_cube_proc.reshape(V_sub, Y_sub * X_sub)
    
    # コアエンジンの実行
    mask_flat = _core_baseline_engine(chunk_flat, degree, max_iter, threshold, dilation)
    mask_sub_cube = mask_flat.reshape(V_sub, Y_sub, X_sub)
    
    # のりしろ部分を切り捨てる (Crop)
    final_mask_chunk = mask_sub_cube[:, y_crop_start:y_crop_end, x_crop_start:x_crop_end]
    
    return final_mask_chunk, y_dest_start, y_dest_end, x_dest_start, x_dest_end

# =====================================================================
# 3. メイン関数 (ユーザーAPI)
# =====================================================================
def get_baseline_mask_3d(data_cube, degree=1, max_iter=3, threshold=3.0, dilation=2, 
                         spatial_sigma=1.0, chunk_size=50, max_workers=None):
    """
    3Dデータキューブから、ベースライン(line-free)領域の3D bit配列(Boolean)を生成する関数。
    IPython/Jupyter環境から安全に呼び出せるマルチスレッド対応版。
    
    Parameters:
    -----------
    data_cube : np.ndarray
        入力の3Dデータ配列。形状は (V, Y, X) を想定。
    degree : int
        ベースラインの多項式フィッティングの次数（デフォルト: 1次）
    max_iter : int
        シグマ・クリッピングの反復回数（デフォルト: 3回）
    threshold : float
        輝線を判定するためのMADの閾値（デフォルト: 3.0シグマ）
    dilation : int
        輝線の裾野として拡張（マスク）する前後チャンネル数（デフォルト: 2）
    spatial_sigma : float
        空間方向のスムージングの強さ（ピクセル単位）。0でスムージングなし。
    chunk_size : int
        一度に処理する空間の一辺のピクセル数。メモリ使用量に影響（デフォルト: 50）
    max_workers : int or None
        並列実行するスレッド数。Noneの場合は搭載CPUコア数を自動適用。
        
    Returns:
    --------
    mask_3d : np.ndarray (dtype=bool)
        入力と同じ (V, Y, X) 形状のBoolean配列。
        True: ベースライン領域, False: 輝線・無効領域
    """
    V, Y_dim, X_dim = data_cube.shape
    mask_3d = np.zeros((V, Y_dim, X_dim), dtype=bool)
    
    # ガウシアンフィルタの影響範囲(のりしろ)の計算
    margin = int(np.ceil(3 * spatial_sigma)) if spatial_sigma > 0 else 0
    
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
        
    futures = []
    
    # ThreadPoolExecutorにより、IPython上でも安全に並列処理を実行
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        for y0 in range(0, Y_dim, chunk_size):
            for x0 in range(0, X_dim, chunk_size):
                
                # のりしろ込みの切り出しインデックス
                y_start = max(0, y0 - margin)
                y_end   = min(Y_dim, y0 + chunk_size + margin)
                x_start = max(0, x0 - margin)
                x_end   = min(X_dim, x0 + chunk_size + margin)
                
                # データの切り出し（View）
                sub_cube = data_cube[:, y_start:y_end, x_start:x_end]
                
                # クロッピング用ローカルインデックス
                y_crop_start = y0 - y_start
                y_crop_end   = y_crop_start + min(chunk_size, Y_dim - y0)
                x_crop_start = x0 - x_start
                x_crop_end   = x_crop_start + min(chunk_size, X_dim - x0)
                
                # 全体配列への書き込みインデックス
                y_dest_start = y0
                y_dest_end   = min(y0 + chunk_size, Y_dim)
                x_dest_start = x0
                x_dest_end   = min(x0 + chunk_size, X_dim)
                
                # スレッドプールへタスクを投げる
                future = executor.submit(
                    _process_single_chunk_worker,
                    sub_cube, y_crop_start, y_crop_end, x_crop_start, x_crop_end,
                    degree, max_iter, threshold, dilation, spatial_sigma,
                    y_dest_start, y_dest_end, x_dest_start, x_dest_end
                )
                futures.append(future)
                
        # 完了したタスクから順次結果を全体配列にマージ
        for future in concurrent.futures.as_completed(futures):
            res_mask, yd_s, yd_e, xd_s, xd_e = future.result()
            mask_3d[:, yd_s:yd_e, xd_s:xd_e] = res_mask
            
    return mask_3d