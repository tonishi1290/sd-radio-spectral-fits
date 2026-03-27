# sd_radio_spectral_fits.map_3d 完全説明書（日本語・詳細版）

## 0. この説明書の目的

この文書は、`sd_radio_spectral_fits.map_3d` で 3D 電波スペクトルキューブを扱う際に、**どの関数を何のために使うか、各パラメーターが何を意味するか、どのような出力が得られるか**を、できるだけ省略せずにまとめた完全版の説明書です。

本説明書は、今回の更新方針も含めた **map_3d 系の現行設計** を前提にしています。特に次の方針を重視しています。

- `MASK3D` は **final signal mask**（最終的な信号マスク）として扱う
- baseline 由来の情報と signal mask は **役割を分離**して考える
- baseline 再解析では `BASE_RMS` を優先し、必要に応じて `RMS` を使う
- 空間 WCS とスペクトル軸は、可能な限り **header を正**として扱う
- stale な HDU（古い `MASK3D`、`MOMENT0`、`BASE_RMS` など）を残して解析を混乱させない
- 温度スケールは `TA*` / `TR*` の概念を保ちつつ、単位系はできるだけ明確に分ける

この説明書では、主として **ユーザーが直接使う公開関数・設定クラス** を対象にします。先頭に `_` が付いた非公開 helper は原則として説明対象外ですが、公開関数の挙動を理解するのに重要な内部仕様は必要に応じて補足します。

---

## 1. パッケージ全体の見取り図

`map_3d` での処理は、大きく次の 5 段階に分かれます。

1. **scantable から速度軸をそろえたスペクトル行列を作る**
   - `run_mapping_pipeline()`
   - `run_ps_mapping_pipeline()`
   - `run_otf_full_pipeline()`

2. **空間的にグリッディングして 3D FITS を作る**
   - OTF: `grid_otf()`（内部）→ `save_map_fits()`
   - PS : `grid_ps()` → `save_ps_map_fits()`

3. **baseline を推定・除去する**
   - `estimate_linefree_mask_1d()`
   - `estimate_linefree_mask_from_cube()`
   - `estimate_ripple_frequencies_fft()`
   - `subtract_baseline_cube()`
   - `subtract_baseline_from_fits()`

4. **3D signal mask と moment map を作る**
   - `estimate_robust_rms()`
   - `generate_cube_mask()`
   - `make_3d_mask_for_existing_fits()`

5. **OTF 特有の補正を行う**
   - `solve_basket_weave_offsets()`
   - `apply_basket_weave_correction()`
   - `run_otf_full_pipeline()`

### 1.1 典型的な使い方

#### OTF の基本フロー

1. `run_mapping_pipeline()` あるいは `run_otf_full_pipeline()` で 3D FITS を作る
2. 必要なら `subtract_baseline_from_fits()` で baseline を引く
3. `make_3d_mask_for_existing_fits()` で `MASK3D` と `MOMENT0` を作る

#### PS の基本フロー

1. `run_ps_mapping_pipeline()` で 3D FITS を作る
2. 必要なら `subtract_baseline_from_fits()` で baseline を引く
3. `make_3d_mask_for_existing_fits()` で `MASK3D` と `MOMENT0` を作る

#### baseline を重点的に見たい場合

1. 既存 FITS に対して `subtract_baseline_from_fits()` を使う
2. `LINEFREE`、`RIPFREQ`、`BASE_RMS`、`BASE_FLG` を確認する
3. signal mask を作りたい場合だけ `make_3d_mask_for_existing_fits()` を使う

---

## 2. 入力データ・出力データの基本概念

### 2.1 scantable とは何か

このパッケージは、`scantable` オブジェクトを入力として受け取る関数が多くあります。想定している `scantable` は、少なくとも次を持つオブジェクトです。

- `scantable.table` : pandas DataFrame 相当の表
- `scantable.data` : スペクトル本体
- `scantable.meta` : メタデータ辞書

よく参照される列・メタデータは次です。

- 座標: `RA`, `DEC`, `GLON`, `GLAT`
- 観測状態: `OBSMODE`, `SCAN`, `IS_TURN`
- 時間: `MJD`, `TIMESTAMP`, `DATE-OBS`, `TIME`
- ノイズ: `BSL_RMS`, `TSYS`, `EXPOSURE`, `INTTIME`, `DUR`
- スペクトル軸: `RESTFREQ` / `RESTFRQ`, `CTYPE1`, `VFRAME`

### 2.2 `GridInput` とは何か

`GridInput` は、グリッディングエンジンに渡す入力データの標準形です。

- `x`, `y` : 各スペクトルの空間位置 [arcsec]
- `spec` : 各スペクトル本体
- `flag` : 採用するかどうかのフラグ
- `time` : 各スペクトルの時刻
- `rms`, `tint`, `tsys` : 重み付けや診断用の補助量
- `scan_id`, `subscan_id`, `scan_dir`, `is_turnaround` : OTF 特有の補助情報

OTF では `spec` は **(ndump, nchan)** を想定します。PS でも同様です。

### 2.3 `GridResult` とは何か

`GridResult` はグリッディング結果の標準形です。

- `cube` : 出力キューブ。内部では通常 **(ny, nx, nchan)**
- `weight_map` : 総重み
- `hit_map` : 何本のスペクトルがそのピクセルに寄与したか
- `mask_map` : 有効ピクセルかどうか
- `rms_map`, `time_map`, `tint_map`, `tsys_map` : 診断用 2D マップ
- `xeff_map`, `yeff_map`, `dr_eff_map_pix`, `neff_map` : OTF グリッディング品質評価用
- `meta` : WCS や `RESTFREQ` などを保持する辞書
- `mask_3d`, `mom0_map` : 解析結果を追加で格納するための領域

### 2.4 FITS に保存されたときの基本 HDU

典型的には次の HDU が現れます。

- Primary HDU: 3D cube
- `WEIGHT`, `HIT`, `MASK`, `TSYS`, `TINT`, `TIME`, `RMS`, `NEFF`, `XEFF`, `YEFF`, `BIAS_PIX`
- baseline 後: `LINEFREE`, `RIPFREQ`, `BASE_RMS`, `BASE_FLG`
- signal mask 後: `MASK3D`, `MOMENT0`

### 2.5 `RMS` と `BASE_RMS` の違い

- `RMS`
  - 主にグリッディング結果としての 2D ノイズ指標
  - OTF / PS の 2D 診断マップとして使う
- `BASE_RMS`
  - baseline 除去後の line-free channels 上の residual RMS
  - baseline の成否や、baseline 後マスク生成の出発点として使う

運用上は、baseline を引いた後の解析では **`BASE_RMS` を優先**するのが自然です。

---

## 3. 設定クラス（dataclass）完全説明

## 3.1 `MapConfig`

`MapConfig` は、主に OTF / 汎用マッピングで使う総合設定クラスです。

### 必須パラメーター

#### `x0: float`
出力グリッドの X 方向の原点オフセット [arcsec] です。`build_spatial_wcs_dict()` と組み合わせて、ローカル平面座標の 0 点がどこに来るかを決めます。

- 単位: arcsec
- 正の意味: `invert_x=True` の運用では見かけ上右左が反転するので、プロットと WCS の見え方に注意が必要です

#### `y0: float`
出力グリッドの Y 方向の原点オフセット [arcsec] です。

#### `nx: int`
X 方向の出力ピクセル数です。

- 1 以上の整数が必要です
- 小さすぎると地図が切れます
- 大きすぎると空白が増え、計算量も増えます

#### `ny: int`
Y 方向の出力ピクセル数です。

#### `cell_arcsec: float`
出力ピクセルサイズ [arcsec] です。

- 小さすぎると oversampling が強くなり、hit が薄くなります
- 大きすぎるとビームより粗くなり、分解能を失います

#### `beam_fwhm_arcsec: float`
望遠鏡ビーム FWHM [arcsec] です。

- kernel 幅をビーム基準で決めるときに使います
- `gwidth_beam`, `jwidth_beam` を使う場合に重要です

### カーネル関連パラメーター

#### `kernel: Literal['gjinc', 'gauss'] = 'gjinc'`
空間グリッディングカーネルの種類です。

- `gjinc`: 一般に OTF グリッディングでよく使うカーネル
- `gauss`: 単純な Gaussian カーネル

#### `gwidth_pix: Optional[float] = 2.10`
Gaussian 部分の幅をピクセル単位で直接指定します。

- `gwidth_beam` を使わない場合はこちらが使われます
- 0 以下は不可です

#### `gwidth_beam: Optional[float] = None`
Gaussian 部分の幅を「ビーム FWHM の何倍か」で指定します。

- 内部では `beam_fwhm_arcsec / cell_arcsec` を掛けてピクセル幅へ変換します
- `gwidth_pix` が指定されていると、通常そちらが優先されます

#### `jwidth_pix: Optional[float] = None`
`gjinc` の Jinc 部分の幅をピクセル単位で指定します。

- `kernel='gjinc'` のときにのみ意味があります
- 未指定時は既定値 1.55 が使われます

#### `jwidth_beam: Optional[float] = None`
Jinc 部分の幅をビーム単位で指定します。

#### `truncate: Union[Literal['first_null'], float] = 'first_null'`
カーネルの打ち切り半径の指定です。

- `'first_null'`: GJINC の first null で切る
- 数値: support radius をピクセル単位で直接指定する

#### `support_radius_pix: Optional[float] = None`
打ち切り半径を明示的にピクセル単位で指定します。

- 指定した場合は `truncate` より優先されます
- support を大きくしすぎると計算量が増えます

#### `chunk_ch: int = 256`
スペクトル軸を何チャネルずつ分割して処理するかの目安です。

- メモリ節約用です
- 大きすぎるとメモリ消費が増えます
- 小さすぎると Python 側 overhead が増えます

#### `dtype: str = 'float32'`
内部計算や出力で使う浮動小数点型です。

- 典型値: `'float32'`, `'float64'`
- `float64` は安全ですがメモリを多く使います

### 重み付けと品質管理

#### `alpha_rms: float = 0.5`
RMS に対する重み指数です。実装では概ね

`q *= (1 / rms^2) ** alpha_rms`

の形で使われます。

- 0 に近いと RMS の違いをあまり重みに反映しません
- 1 に近いと RMS の小さい点を強く重視します

#### `beta_tint: float = 0.0`
積分時間 `tint` の重み指数です。概ね

`q *= tint ** beta_tint`

の形で使われます。

- 0 なら tint を無視します
- 1 に近づくほど長時間積分を重視します

#### `weight_clip_quantile: Optional[float] = 0.95`
重み `q` の上位分位点でクリップします。

- 一部の極端な高重み点が全体を支配するのを防ぎます
- `None` なら使いません

#### `weight_clip_max: Optional[float] = None`
重みの絶対上限を指定します。

#### `exclude_turnaround: bool = True`
`is_turnaround=True` の dump を除外するかどうかです。

- OTF の折り返し区間は通常品質が悪いので `True` が推奨です

### 推定器 / QC / 数値安全係数

#### `estimator: Literal['avg', 'plane'] = 'avg'`
グリッディング時の推定器です。

- 現在の実装では `plane` は **未実装** で、指定すると `NotImplementedError` です
- 実質 `avg` を使います

#### `n_min_avg: int = 2`
加重平均を成立させるために必要な最低近傍点数です。

- これ未満だと有効ピクセルと見なされません

#### `n_min_plane: int = 6`
平面近似用の最低点数です。

- 現在の `plane` 未実装のため、将来用の意味合いが強いです

#### `cond_max: float = 1e6`
行列条件数の上限です。

- 将来的な平面フィット等の数値安定性評価用です

#### `dr_eff_warn_pix: float = 0.2`
有効位置ずれ `dr_eff_map_pix` に対する警告閾値 [pixel] です。

#### `eps_u0: float = 0.01`
GJINC カーネルの `u=0` 近傍での安全処理に使うしきい値です。

#### `eps_weight_sum: float = 1e-8`
総重みがほぼ 0 かどうかを判定する数値下限です。

### 出力制御

#### `fill_nan_for_invalid: bool = True`
無効ピクセルを `NaN` で埋めるかどうかです。

- `True` なら無効画素が明示的に `NaN` になります
- `False` なら 0 埋めになるため、後段解析で誤解を招くことがあります

#### `emit_diag_maps: bool = True`
`XEFF`, `YEFF`, `BIAS_PIX` などの診断マップを作るかどうかです。

#### `emit_neff_map: bool = True`
`NEFF` を出すかどうかです。

#### `emit_rms_map: bool = True`
`RMS` を出すかどうかです。

#### `emit_time_map: bool = True`
`TIME` を出すかどうかです。

#### `emit_tint_map: bool = True`
`TINT` を出すかどうかです。

#### `emit_tsys_map: bool = True`
`TSYS` を出すかどうかです。

### 解析・出力設定

#### `generate_mask: bool = False`
グリッディング後にマスクと moment を自動生成する構想用フラグです。

- 現状の公開フローでは、通常 `make_3d_mask_for_existing_fits()` を別段で呼ぶ方が明確です

#### `mask_method: str = 'smooth_mask'`
自動マスク生成を使う場合の方法名です。

#### `mask_sigma: float = 3.0`
`simple` 法の閾値です。

#### `mask_high_snr: float = 3.0`
`smooth_mask` 系の high threshold です。

#### `mask_low_snr: float = 1.5`
`smooth_mask` 系の low threshold です。

#### `mask_min_vol: int = 27`
`smooth_mask` の 3D 最小体積閾値です。

#### `mask_sigma_v: float = 2.0`
`derivative` 法の速度方向 Gaussian 幅です。

#### `mask_deriv_snr: float = 3.0`
`derivative` 法の検出閾値です。

#### `mask_dilation: int = 2`
`derivative` 法の速度方向 dilation 回数です。

#### `mask_compression: Optional[str] = 'PLIO_1'`
マスクを FITS へ書く際の圧縮形式です。

- 互換性優先なら `None`
- 圧縮したいなら `PLIO_1` など

#### `dv_kms: Optional[float] = None`
標準化後の速度分解能 [km/s] です。

- `Standardizer.get_matrix(dv=...)` に渡すための値です
- `None` なら入力依存になります

#### `output_prefix: str = 'map'`
出力ファイル名を組み立てる際の接頭辞です。

#### `backend: str = 'numpy'`
内部バックエンドです。

- 現状の主実装は `numpy`
- 将来 `numba` 等に拡張する余地があります

---

## 3.2 `GridInput`

`GridInput` の各フィールドは次の通りです。

#### `x: np.ndarray`
各 dump の X 座標 [arcsec]。shape は `(ndump,)`。

#### `y: np.ndarray`
各 dump の Y 座標 [arcsec]。shape は `(ndump,)`。

#### `spec: np.ndarray`
各 dump のスペクトル本体。shape は通常 `(ndump, nchan)`。

#### `flag: np.ndarray`
各 dump を採用するかどうか。`True` / 非 0 が有効。

#### `time: np.ndarray`
各 dump の時刻。MJD UTC day を推奨。

#### `rms: Optional[np.ndarray]`
各 dump の基礎 RMS。重み付けに使います。

#### `tint: Optional[np.ndarray]`
各 dump の積分時間 [s]。

#### `tsys: Optional[np.ndarray]`
各 dump のシステム温度 [K]。

#### `scan_id: Optional[np.ndarray]`
OTF で basket-weave 補正を行う際に重要です。

#### `subscan_id: Optional[np.ndarray]`
サブスキャン単位の識別子。現状の主フローでは必須ではありません。

#### `scan_dir: Optional[np.ndarray]`
走査方向の識別子。

#### `is_turnaround: Optional[np.ndarray]`
折り返し dump を示します。`exclude_turnaround=True` と組み合わせます。

---

## 3.3 `GridResult`

#### `cube: np.ndarray`
グリッディング後の 3D キューブ。内部形状は通常 `(ny, nx, nchan)`。

#### `weight_map: np.ndarray`
各ピクセルの総重み。

#### `hit_map: np.ndarray`
各ピクセルに何本の dump が寄与したか。

#### `mask_map: np.ndarray`
そのピクセルが有効かどうか。

#### `xeff_map`, `yeff_map`
有効中心位置のずれに関する情報。

#### `dr_eff_map_pix`
実効位置ずれの大きさ [pixel]。

#### `neff_map`
有効サンプル数。

#### `rms_map`
2D のノイズマップ。

#### `time_map`
平均観測時刻。

#### `tint_map`
実効積分時間。

#### `tsys_map`
システム温度。

#### `meta: Optional[dict]`
FITS 書き出し時に使うメタ情報。典型的には

- `RESTFREQ`
- `SPECSYS`
- `coord_sys`
- `projection`
- `x0_arcsec`, `y0_arcsec`
- `cdelt1_arcsec`, `cdelt2_arcsec`
- `bmaj_eff_arcsec`

などが入ります。

#### `mask_3d`, `mom0_map`
後段解析でマスクや moment0 を一時保持するための領域です。

---

## 3.4 `PSMapConfig`

PS マッピング専用の設定です。

#### `coord_sys: str = 'icrs'`
入力座標系です。`'icrs'` や `'galactic'` などを想定します。

#### `projection: str = 'SFL'`
局所平面への投影法です。

- `SFL` / `GLS` / `SINE` 系
- `CAR`

#### `ref_lon: Optional[float] = None`
基準経度 [deg]。

- `None` のときは入力座標の中央値から自動決定します

#### `ref_lat: Optional[float] = None`
基準緯度 [deg]。

#### `x_grid: float = 30.0`
X 方向の格子間隔 [arcsec]。

#### `y_grid: float = 30.0`
Y 方向の格子間隔 [arcsec]。

#### `grid_anchor_offsets: Tuple[float, float] = (0.0, 0.0)`
格子位相を固定するためのアンカー位置 [arcsec]。

- 観測ごとにグリッド位相がずれるのを避けたいときに重要です

#### `grid_bounds_offsets: Optional[((x1,y1),(x2,y2))] = None`
格子の採用範囲を明示指定します。

- 指定しなければ入力点の min/max から決まります

#### `grid_tol: float = 0.3`
近傍格子点へ割り当てる許容率です。

- 例えば `0.3` なら、格子間隔の 30% 以内の点だけ採用します

#### `invert_x: bool = True`
FITS 表示規則に合わせて X 軸の正方向を左にするかどうかです。

#### `combine: Literal['mean','median','sum'] = 'mean'`
1 セルに複数スペクトルが入ったときの結合方法です。

- `mean`: 加重平均
- `median`: 中央値
- `sum`: 和

#### `weight_mode: Literal['uniform','rms','tint'] = 'rms'`
PS セル内結合時の重みの付け方です。

- `uniform`: 全点同じ重み
- `rms`: `1/rms^2`
- `tint`: 積分時間重み

#### `return_cell_members: bool = False`
セルごとの構成メンバーを返す構想用フラグです。現状主フローではあまり使いません。

#### `verbose: bool = True`
ログ出力を出すかどうかです。

---

## 3.5 `LineFreeConfig`

baseline 用の line-free channel 自動推定設定です。

#### `smooth_width: int = 51`
中央値フィルターの幅 [channel]。

- 大きいほどゆっくりしたベースラインの推定に向きます
- 広い線を baseline と誤解しにくくなります
- ただし ripple を平滑化しすぎる危険もあります

#### `sigma: float = 4.0`
外れ値判定のしきい値です。

#### `iters: int = 6`
反復回数です。

#### `pad_chan: int = 3`
検出した line channel を何チャネル分拡張するかです。

#### `min_linefree_frac: float = 0.35`
最終的な line-free fraction の下限です。

- これを下回ると baseline fit が不安定になるため、実装では `ValueError` を投げます

---

## 3.6 `RippleConfig`

standing wave / ripple 推定設定です。

#### `nfreq: int = 2`
推定する ripple 周波数の本数です。

#### `period_range_chan: Tuple[float, float] = (20.0, 400.0)`
探索する周期範囲 [channel] です。

- 周波数ではなく周期で指定する点に注意してください
- 例えば `(20, 400)` は比較的遅い ripple を想定します

#### `min_separation: float = 0.002`
選ぶピーク間の最小周波数分離 [cycles/channel] です。

#### `window: str = 'hann'`
FFT 前に掛ける窓関数です。

- `'hann'`
- `'none'`

---

## 3.7 `BaselineConfig`

baseline モデルの設定です。

#### `poly_order: int = 1`
多項式 baseline の次数です。

- 0: 定数
- 1: 直線
- 2 以上: 曲率を持つ基線

次数を上げすぎると、実線成分まで baseline が吸収する危険があります。

#### `ripple: bool = True`
ripple 項をモデルに含めるかどうかです。

#### `robust: bool = False`
IRLS 的な軽量 robust reweighting を行うかどうかです。

- `True` の方が line-free mask の誤差に強いことがあります
- ただし大きい cube では遅くなります

#### `rcond: Optional[float] = None`
`numpy.linalg.lstsq()` に渡す `rcond` です。

#### `chunk_pix: int = 65536`
1 回に処理する空間ピクセル数です。

- 大きいほど速いがメモリを使います
- 小さいほど安全だが遅くなります

---

## 4. OTF / 汎用マッピング API

## 4.1 `run_mapping_pipeline()`

### 役割

scantable から速度軸を標準化し、座標投影し、温度スケールを処理して、OTF グリッディングを実行し、最終的に 3D FITS を保存する統合関数です。

### シグネチャ

```python
run_mapping_pipeline(
    scantable,
    config: MapConfig,
    output_fits: str,
    coord_sys: str = 'icrs',
    projection: str = 'SFL',
    out_scale: str = 'TA*',
    dv_kms: float = None,
    ref_lon: float = None,
    ref_lat: float = None,
)
```

### パラメーター

#### `scantable`
入力 scantable です。

#### `config: MapConfig`
グリッディング設定一式です。

#### `output_fits: str`
出力 FITS ファイル名です。

#### `coord_sys: str = 'icrs'`
座標系。典型的には `'icrs'` または `'galactic'`。

#### `projection: str = 'SFL'`
空間投影法。`project_to_plane()` と `build_spatial_wcs_dict()` に渡されます。

#### `out_scale: str = 'TA*'`
出力温度スケール。

- `'TA*'`: 入力を主に TA* として扱う
- `'TR*'`: `TA* -> TR*` 変換を行う

#### `dv_kms: float = None`
速度軸標準化後のチャネル幅 [km/s]。

- `None` なら `Standardizer` に依存します
- 明示した方が再現性が高くなります

#### `ref_lon: float = None`
投影中心経度 [deg]。未指定なら入力座標の中央値。

#### `ref_lat: float = None`
投影中心緯度 [deg]。未指定なら入力座標の中央値。

### 戻り値

`GridResult` を返します。

### 処理内容

1. 必要なら `SIG_*` 列を落とす
2. `Standardizer` で速度軸を標準化する
3. `RA/DEC` あるいは `GLON/GLAT` を局所平面へ投影する
4. `BEAMEFF` と温度スケールを処理する
5. `GridInput` を組み立てる
6. `grid_otf()` を実行する
7. `save_map_fits()` で Multi-Extension FITS を保存する

### 注意点

- `coord_sys` と `projection` は WCS の見え方に直結します
- `out_scale='TR*'` では `ta_to_tr()` が適用されます
- `BSL_RMS` があると `_extract_rms()` で重みへ反映されます
- 実装中に速度軸の debug print が残っている版もあります

---

## 4.2 `create_grid_input()`

### 役割

scantable から `GridInput` を安全に組み立てる helper です。basket-weave や独自グリッディング前処理で便利です。

### シグネチャ

```python
create_grid_input(scantable, ref_coord=None, frame='ICRS', projection='SFL')
```

### パラメーター

#### `scantable`
入力 scantable。

#### `ref_coord = None`
`calc_mapping_offsets()` に渡す基準座標。

#### `frame = 'ICRS'`
座標系フレーム。

#### `projection = 'SFL'`
投影法。

### 戻り値

`GridInput`

### 内部で行うこと

- `calc_mapping_offsets()` を用いて `OFS_LON`, `OFS_LAT` を計算
- `OBSMODE` から `flag` を作成
- `_resolve_time_array()` で時刻配列を作成
- `SCAN`, `IS_TURN` があれば補助情報として格納

---

## 4.3 `run_otf_full_pipeline()`

### 役割

OTF 専用の上位関数で、必要なら basket-weave 補正を行ってから `run_mapping_pipeline()` に処理を委譲します。

### シグネチャ

```python
run_otf_full_pipeline(
    scantable,
    config,
    output_fits: str,
    do_basket_weave: bool = True,
    **kwargs,
)
```

### パラメーター

#### `scantable`
入力 scantable。

#### `config`
通常 `MapConfig` を想定します。

#### `output_fits: str`
出力 FITS ファイル名。

#### `do_basket_weave: bool = True`
`True` のとき basket-weave 補正を実行します。

#### `**kwargs`
そのまま `run_mapping_pipeline()` に渡します。代表例:

- `coord_sys`
- `projection`
- `out_scale`
- `dv_kms`
- `ref_lon`
- `ref_lat`

### 処理内容

1. `create_grid_input()` で `GridInput` を作る
2. `solve_basket_weave_offsets()` で scan ごとの定数オフセットを推定
3. scantable のデータ本体へオフセットを直接引く
4. `run_mapping_pipeline()` を呼ぶ

### 注意点

- basket-weave は **scantable.data を直接変更**します
- scan ID が正しくないとオフセットが崩れます

---

## 5. PS マッピング API

## 5.1 `grid_ps()`

### 役割

PS 観測用のグリッディングコアです。近傍格子へ割り当ててセルごとに結合します。

### シグネチャ

```python
grid_ps(input_data: GridInput, config: PSMapConfig) -> GridResult
```

### パラメーター

#### `input_data: GridInput`
PS 用に整形された入力データ。

#### `config: PSMapConfig`
PS 専用設定。

### 主な処理

- `flag`, `is_turnaround`, NaN の有無で有効データを絞る
- `weight_mode` に応じて重みを決める
- `grid_bounds_offsets` または入力範囲から格子境界を決める
- 各スペクトルを最寄りの格子セルへ割り当てる
- `combine` に応じてセル内スペクトルを統合する
- `RMS`, `TSYS`, `TINT`, `TIME` を 2D マップにする

### `combine` の意味

- `mean`: 重み付き平均。通常の推奨値
- `median`: 外れ値に強いが、重みの扱いが単純ではない
- `sum`: 合算したいとき用

### `weight_mode` の意味

- `uniform`: 全点同重み
- `rms`: `1 / rms^2`
- `tint`: 積分時間比例

### 注意点

- `grid_tol` が小さすぎると点がどのセルにも割り当たらず、空白が増えます
- `invert_x=True` は FITS 表示の慣例に合わせるためのものです

---

## 5.2 `build_ps_wcs_dict()`

### 役割

PS 用の空間 WCS ヘッダ辞書を組み立てます。

### シグネチャ

```python
build_ps_wcs_dict(
    coord_sys: str,
    projection: str,
    lon0: float,
    lat0: float,
    x0_arcsec: float,
    y0_arcsec: float,
    cdelt1_arcsec: float,
    cdelt2_arcsec: float,
    nx: int,
    ny: int,
) -> dict
```

### パラメーター

#### `coord_sys`
`ICRS`, `RADEC`, `GALACTIC` など。

#### `projection`
`SFL`, `GLS`, `CAR` など。

#### `lon0`, `lat0`
基準天球座標 [deg]。

#### `x0_arcsec`, `y0_arcsec`
局所平面での基準点 [arcsec]。

#### `cdelt1_arcsec`, `cdelt2_arcsec`
ピクセルスケール [arcsec/pix]。

#### `nx`, `ny`
画像サイズ。現実装では主に整合確認用です。

### 注意点

- `SFL` のとき `CRVAL2=0` とし、`CRPIX2` をずらして実座標へ合わせる流儀を取っています
- `cdelt1_arcsec` は `invert_x` の結果に応じて負になることがあります

---

## 5.3 `save_ps_map_fits()`

### 役割

PS グリッディング結果を FITS に保存します。

### シグネチャ

```python
save_ps_map_fits(
    grid_res: GridResult,
    v_tgt: np.ndarray,
    lon0: float,
    lat0: float,
    out_path: str,
    out_scale: str = 'TA*',
    rep_beameff: float = 1.0,
)
```

### パラメーター

#### `grid_res: GridResult`
PS グリッディング結果。

#### `v_tgt: np.ndarray`
出力速度軸 [km/s]。

#### `lon0`, `lat0`
WCS 用基準座標 [deg]。

#### `out_path: str`
出力ファイル名。

#### `out_scale: str = 'TA*'`
出力温度スケールラベル。

#### `rep_beameff: float = 1.0`
代表的な主ビーム効率。

### 出力される主な HDU

- Primary: cube
- `HIT`
- `RMS`
- `TSYS`
- `TINT`

PS では OTF より診断マップを絞っています。

---

## 5.4 `run_ps_mapping_pipeline()`

### 役割

PS 観測用の統合パイプラインです。速度軸の標準化から FITS 出力まで行います。

### シグネチャ

```python
run_ps_mapping_pipeline(
    scantable,
    config: PSMapConfig,
    output_fits: str,
    out_scale: str = 'TA*',
    dv_kms: float = None,
    vmin_kms: float = None,
    vmax_kms: float = None,
)
```

### パラメーター

#### `scantable`
入力 scantable。

#### `config: PSMapConfig`
PS グリッディング設定。

#### `output_fits: str`
出力 FITS 名。

#### `out_scale: str = 'TA*'`
出力温度スケール。

#### `dv_kms: float = None`
速度軸標準化後の分解能 [km/s]。

#### `vmin_kms: float = None`
速度範囲下限 [km/s]。

#### `vmax_kms: float = None`
速度範囲上限 [km/s]。

### 処理内容

- 必要に応じて `SIG_*` 列を除去
- `Standardizer` で速度軸をそろえる
- 座標投影を行う
- `BEAMEFF` と温度スケールを処理する
- `GridInput` を構築する
- `grid_ps()` を実行する
- `save_ps_map_fits()` で保存する

---



---

## 5.4.1 baseline 済み PS scantable `sc_all_coadded_12co` を解析するときの実践例

ここでは、Position Switching 観測で baseline fitting 済みの scantable `sc_all_coadded_12co` を想定します。

### どこから入るべきか

このケースでは、最初の入口は `run_ps_mapping_pipeline()` です。

理由は次の通りです。

- `sc_all_coadded_12co` は scantable である
- baseline fitting はすでに終わっている
- `subtract_baseline_from_fits()` は FITS cube に対する baseline 関数である

したがって、まずやるべきことは **cube 化** です。

### 標準的な流れ

1. `run_ps_mapping_pipeline()` で 3D FITS を作る
2. 出力された `HIT`, `RMS`, `TSYS`, `TINT` を確認する
3. 必要なら `make_3d_mask_for_existing_fits()` で `MASK3D` と `MOMENT0` を追加する

### 実行例

```python
import numpy as np
from sd_radio_spectral_fits.map_3d.ps_gridder import PSMapConfig, run_ps_mapping_pipeline
from sd_radio_spectral_fits.map_3d.cube_analysis import make_3d_mask_for_existing_fits

x_grid = np.arange(-360.0, 360.0 + 36.0, 36.0)
y_grid = np.arange(-360.0, 360.0 + 36.0, 36.0)

ps_cfg = PSMapConfig(
    coord_sys='icrs',
    projection='SFL',
    ref_lon=83.809,
    ref_lat=-5.372639,
    x_grid=x_grid,
    y_grid=y_grid,
    grid_tol=12.0,
    combine='mean',
    weight_mode='rms',
    invert_x=True,
    verbose=True,
)

run_ps_mapping_pipeline(
    scantable=sc_all_coadded_12co,
    config=ps_cfg,
    output_fits='sc_all_coadded_12co_cube.fits',
    out_scale='TA*',
    dv_kms=0.2,
    vmin_kms=-30.0,
    vmax_kms=55.0,
)

make_3d_mask_for_existing_fits(
    input_fits='sc_all_coadded_12co_cube.fits',
    output_fits='sc_all_coadded_12co_cube_masked.fits',
    rms_ext='RMS',
    method='smooth_mask_lite',
    moment_spectral_slab=(-10.0, 30.0),
    overwrite=True,
)
```

### ここでの判断基準

- baseline 済み scantable から始めるなら、最初から cube に baseline を掛け直さない
- まず `RMS` と `HIT` を見て、PS グリッディング結果の品質を確認する
- signal mask は cube のあとで作る
- もし cube 化後に広域の基線ずれが気になるなら、元の cube を残した上で `subtract_baseline_from_fits()` を試す

### `RMS` と `BASE_RMS` の使い分け

この例で最初に見るべきは `RMS` です。

- `RMS`
  - PS pipeline が出す画素単位の雑音指標
- `BASE_RMS`
  - 後段で cube baseline を掛けたときの line-free 残差指標

つまり、`sc_all_coadded_12co` から素直に cube を作った段階では、まず `RMS` を参照します。

## 6. baseline 解析 API

## 6.1 `estimate_linefree_mask_1d()`

### 役割

1 本の 1D スペクトルから line-free channel を推定します。

### シグネチャ

```python
estimate_linefree_mask_1d(
    spectrum: np.ndarray,
    cfg: LineFreeConfig = LineFreeConfig(),
) -> np.ndarray
```

### パラメーター

#### `spectrum: np.ndarray`
1 次元スペクトルです。shape は `(nchan,)`。

#### `cfg: LineFreeConfig`
line-free 推定設定。

### 戻り値

`linefree: np.ndarray[bool]`

- `True`: baseline fit に使ってよいチャネル
- `False`: line の可能性があるため除外すべきチャネル

### 処理の流れ

1. 1D median filter で遅い baseline を推定
2. residual を作る
3. MAD で robust std を求める
4. `sigma` を超えるチャネルを line とみなす
5. それを反復し、最後に `pad_chan` 分だけ拡張する

### 注意点

- line が非常に広いと median filter 幅の選び方に依存します
- `min_linefree_frac` を下回ると止まります

---

## 6.2 `estimate_linefree_mask_from_cube()`

### 役割

3D cube 全体から、空間的に集約した 1 本の代表スペクトルを作り、そこから **global な line-free mask** を推定します。

### シグネチャ

```python
estimate_linefree_mask_from_cube(
    cube_data: np.ndarray,
    cfg: LineFreeConfig = LineFreeConfig(),
    *,
    agg: str = 'median',
    max_pix: int = 200000,
    seed: int = 0,
) -> np.ndarray
```

### パラメーター

#### `cube_data: np.ndarray`
shape `(nchan, ny, nx)` の cube。

#### `cfg: LineFreeConfig`
line-free 推定設定。

#### `agg: str = 'median'`
空間集約方法。

- `'median'`: 外れ値や強い局所線に強い
- `'mean'`: 平均。広がった弱線の感度はよいが外れ値に弱い

#### `max_pix: int = 200000`
空間画素数が大きいときのサンプリング上限。

- 巨大 cube で計算を軽くするため、最大この数だけランダムサンプリングします

#### `seed: int = 0`
サンプリングの乱数 seed。

### 戻り値

global な `(nchan,)` の line-free mask。

### 注意点

- ここで得られるのは **global mask** です
- 位置ごとに異なる line-free を作るわけではありません

---

## 6.3 `estimate_ripple_frequencies_fft()`

### 役割

line-free channel 上で FFT を行い、standing wave / ripple の代表周波数を推定します。

### シグネチャ

```python
estimate_ripple_frequencies_fft(
    spectrum: np.ndarray,
    linefree_mask: np.ndarray,
    rcfg: RippleConfig = RippleConfig(),
    *,
    poly_order_pre: int = 1,
) -> List[float]
```

### パラメーター

#### `spectrum: np.ndarray`
1D スペクトル `(nchan,)`。

#### `linefree_mask: np.ndarray`
`True` が line-free のマスク `(nchan,)`。

#### `rcfg: RippleConfig`
ripple 推定設定。

#### `poly_order_pre: int = 1`
FFT 前に差し引く低次多項式の次数です。

- slow baseline を除いた上で ripple だけを見やすくします

### 戻り値

`List[float]`

- 単位は `cycles/channel`
- 周期 [channel] は `1 / freq`

### 注意点

- linefree_mask が悪いと実線を ripple と誤認します
- `period_range_chan` は探索結果に強く効きます

---

## 6.4 `subtract_baseline_cube()`

### 役割

既に line-free mask と必要なら ripple 周波数が与えられた cube に対して、各スペクトルの baseline を引きます。

### シグネチャ

```python
subtract_baseline_cube(
    cube_data: np.ndarray,
    *,
    linefree_mask: np.ndarray,
    bcfg: BaselineConfig = BaselineConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    return_qc: bool = True,
) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]
```

### パラメーター

#### `cube_data: np.ndarray`
入力 cube。shape は **(nchan, ny, nx)**。

#### `linefree_mask: np.ndarray`
shape `(nchan,)` の bool。`True` のチャネルだけ baseline fit に使います。

#### `bcfg: BaselineConfig`
baseline モデル設定。

#### `ripple_freqs: Optional[Sequence[float]] = None`
ripple 周波数リスト [cycles/channel]。

- `bcfg.ripple=True` でも `ripple_freqs=None` なら、実質的には多項式 baseline のみになります

#### `return_qc: bool = True`
`True` のとき residual RMS map と flag map も返します。

### 戻り値

1. `out_cube`
   - baseline 除去後 cube `(nchan, ny, nx)`
2. `resid_rms`
   - line-free 上での residual RMS map `(ny, nx)`
3. `flag`
   - fit 状態の flag map `(ny, nx)`
   - 典型的には `0=OK`, `2=lstsq failed`

### 重要な制約

- line-free point 数 `nlf` は最低限必要で、実装では
  `nlf < max(ncoef + 2, 8)`
  なら止まります
- つまり、定数 baseline でも最低 8 チャネル以上の line-free が必要です

### 注意点

- この関数は **チャネル空間**で baseline を引きます
- velocity window を使いたい場合は、事前にチャネルへ変換するか wrapper を使います

---

## 6.5 `subtract_baseline_from_fits()`

### 役割

既存の 3D FITS を読み、line-free mask と ripple を推定または読み込み、baseline を引き、QC HDU を付けて新しい FITS を保存する wrapper です。

### シグネチャ

```python
subtract_baseline_from_fits(
    input_fits: str,
    output_fits: str,
    *,
    cube_ext: Optional[Union[int, str]] = None,
    linefree_cfg: LineFreeConfig = LineFreeConfig(),
    linefree_mask: Optional[np.ndarray] = None,
    ripple_cfg: RippleConfig = RippleConfig(),
    ripple_freqs: Optional[Sequence[float]] = None,
    baseline_cfg: BaselineConfig = BaselineConfig(),
    add_qc_hdus: bool = True,
    overwrite: bool = True,
) -> None
```

### パラメーター

#### `input_fits: str`
入力 FITS パス。

#### `output_fits: str`
出力 FITS パス。

#### `cube_ext: Optional[Union[int, str]] = None`
どの HDU を cube とみなすかです。

- `None`: Primary を優先し、必要なら最初の image HDU を探す
- `int`: HDU index
- `str`: HDU name

#### `linefree_cfg: LineFreeConfig`
`linefree_mask` を自動推定する場合の設定です。

#### `linefree_mask: Optional[np.ndarray] = None`
手動の line-free mask。

- 与えた場合は `linefree_cfg` を使いません
- shape は `(nchan,)` が必要です

#### `ripple_cfg: RippleConfig`
`ripple_freqs` を自動推定する場合の設定。

#### `ripple_freqs: Optional[Sequence[float]] = None`
手動で指定する ripple 周波数 [cycles/channel]。

#### `baseline_cfg: BaselineConfig`
baseline fit 設定。

#### `add_qc_hdus: bool = True`
`LINEFREE`, `RIPFREQ`, `BASE_RMS`, `BASE_FLG` を追加するかどうか。

#### `overwrite: bool = True`
既存ファイルを上書きするかどうか。

### 出力内容

- cube 本体は baseline 除去後 data で置き換えられます
- `LINEFREE`: 1D uint8
- `RIPFREQ`: BinTable
- `BASE_RMS`: 2D image
- `BASE_FLG`: 2D image

### 注意点

- baseline fit は channel space で行います
- velocity slab から line-free を定めたい場合は、別途チャネルマスクに落としてから使います
- 再処理時には、古い QC HDU を残さない運用が重要です

---

## 7. 3D マスク・Moment 解析 API

## 7.1 `estimate_robust_rms()`

### 役割

`spectral-cube` の `mad_std` を使って、スペクトル方向に沿った robust RMS を 2D map として求めます。

### シグネチャ

```python
estimate_robust_rms(
    cube: SpectralCube,
    axis: int = 0,
    *,
    how: str = 'auto',
    ignore_nan: Optional[bool] = None,
) -> np.ndarray
```

### パラメーター

#### `cube: SpectralCube`
入力 cube。

#### `axis: int = 0`
どの軸に沿って RMS を取るか。通常 `0` はスペクトル軸です。

#### `how: str = 'auto'`
`spectral-cube` 側の計算方法指定。

#### `ignore_nan: Optional[bool] = None`
NaN を無視するかどうか。`spectral-cube` の版により使えたり使えなかったりします。

### 戻り値

2D `rms_map`。

---

## 7.2 `generate_cube_mask()`

### 役割

指定した方法で 3D signal mask を作ります。

### シグネチャ

```python
generate_cube_mask(
    cube: SpectralCube,
    rms_map,
    method: str = 'smooth_mask',
    *,
    fill_value: float = np.nan,
    data_dtype: np.dtype = np.float32,
    **kwargs,
) -> np.ndarray
```

### 共通パラメーター

#### `cube: SpectralCube`
入力 cube。

#### `rms_map`
2D RMS map。

#### `method: str = 'smooth_mask'`
方法名。

- `'simple'`
- `'smooth_mask'`
- `'smooth_mask_lite'`
- `'derivative'`

#### `fill_value: float = np.nan`
欠損値の埋め値。

#### `data_dtype: np.dtype = np.float32`
内部計算型。

### 戻り値

`mask_3d: np.ndarray`。shape `(nchan, ny, nx)`、dtype `uint8`。

### 方法ごとの `kwargs`

#### `method='simple'`

- `sigma: float = 3.0`
  - 単純に `data > sigma * rms` でマスクを作ります

#### `method='smooth_mask'`

- `smooth_sigma = (2.0, 1.0, 1.0)`
  - Gaussian smoothing の幅 `(spec, y, x)`
- `high_snr: float = 3.0`
  - コア領域抽出閾値
- `low_snr: float = 1.5`
  - 拡張領域抽出閾値
- `min_vol: int = 27`
  - 3D connected component の最小体積
- `den_eps: float = 1e-4`
  - smoothing 時の分母ゼロ防止閾値

特徴:
- 厳密だが重い
- 3D labeling を使うため大きい cube ではメモリを使う

#### `method='smooth_mask_lite'`

- `smooth_sigma = (2.0, 1.0, 1.0)`
- `high_snr: float = 4.0`
- `low_snr: float = 2.0`
- `den_eps: float = 1e-4`
- `min_area: int = 0`
  - 2D footprint の最小面積
- `min_nchan: int = 0`
  - 各空間画素で最低必要なチャネル数
- `prop_connectivity: int = 1`
  - binary propagation の連結性
- `prop_reach: int = 1`
  - 伝播近傍の広がり

特徴:
- 3D labeling を避ける軽量版
- 大きい cube で使いやすい

#### `method='derivative'`

- `sigma_v: float = 2.0`
  - 速度方向 smoothing 幅
- `deriv_snr: float = 3.0`
  - 2 次微分検出閾値
- `dilation_iters: int = 2`
  - スペクトル軸方向 dilation 回数

特徴:
- 線幅やピーク構造に敏感
- baseline ripple やノイズ構造にも反応し得る

### 注意点

- `MASK3D` は **最終 signal mask** として扱うべきです
- provisional mask を作りたい場合は、別系統の HDU 名を使う方が混乱が少ないです

---

## 7.3 `append_analysis_hdus_to_fits()`

### 役割

既存の HDUList に `MASK3D` と `MOMENT0` を追加または置換します。

### シグネチャ

```python
append_analysis_hdus_to_fits(
    hdul: fits.HDUList,
    mask_3d: np.ndarray,
    mom0_hdu,
    base_header: fits.Header,
    *,
    mask_compression: Optional[str] = 'PLIO_1',
) -> None
```

### パラメーター

#### `hdul: fits.HDUList`
追記対象。

#### `mask_3d: np.ndarray`
3D signal mask。

#### `mom0_hdu`
`spectral-cube` の `moment0()` が返した HDU 相当。

#### `base_header: fits.Header`
元 cube のヘッダ。`MASK3D` 用の WCS を引き継ぐために使います。

#### `mask_compression: Optional[str] = 'PLIO_1'`
マスク圧縮形式。

### 挙動

- 既存 `MASK3D`, `MOMENT0` は削除してから追加します
- `base_header` は `_clean_header_for_image_like_hdu()` でクリーニングして使います

---

## 7.4 `make_3d_mask_for_existing_fits()`

### 役割

保存済み 3D FITS を読み、必要に応じて km/s 変換し、RMS map を決め、3D mask を作り、`MOMENT0` を積分して FITS に追記保存します。

### シグネチャ

```python
make_3d_mask_for_existing_fits(
    input_fits: str,
    output_fits: Optional[str] = None,
    *,
    cube_ext: Optional[Union[int, str]] = None,
    rms_ext: str = 'RMS',
    method: str = 'smooth_mask',
    mask_compression: Optional[str] = 'PLIO_1',
    convert_to_kms: bool = True,
    velocity_convention: str = 'radio',
    rest_value_hz: Optional[float] = None,
    require_kms: bool = False,
    file_format: Optional[str] = None,
    use_dask: Optional[bool] = None,
    fill_value: float = np.nan,
    overwrite: bool = True,
    strip_checksum: bool = True,
    rms_how: str = 'auto',
    rms_ignore_nan: Optional[bool] = None,
    rms_spectral_slab = None,
    moment_how: str = 'auto',
    moment_spectral_slab = None,
    **kwargs,
) -> None
```

### 出力先関連

#### `input_fits: str`
入力 FITS。

#### `output_fits: Optional[str] = None`
出力先。`None` なら入力ファイルを上書きします。

#### `cube_ext: Optional[Union[int, str]] = None`
どの HDU を cube として読むか。

### RMS 関連

#### `rms_ext: str = 'RMS'`
既存 RMS HDU 名です。

- 仕様上は baseline 後解析では `BASE_RMS` を優先したいことが多いです
- 実装ブランチによっては `'auto'` 相当の扱いに拡張されている場合があります

#### `rms_how: str = 'auto'`
`spectral-cube` で RMS を推定するときの `how`。

#### `rms_ignore_nan: Optional[bool] = None`
NaN の扱い。

#### `rms_spectral_slab`
RMS 推定に使う速度範囲です。

- `(-20, 0)` のような float 2 つ、または `Quantity` で指定します
- float の単位は、その時点でのスペクトル軸単位に依存します

### マスク関連

#### `method: str = 'smooth_mask'`
マスク方法。

#### `mask_compression: Optional[str] = 'PLIO_1'`
`MASK3D` の圧縮形式。

#### `**kwargs`
`generate_cube_mask()` に渡される method-specific パラメーターです。

### スペクトル軸関連

#### `convert_to_kms: bool = True`
スペクトル軸を km/s に変換しようとするかどうか。

#### `velocity_convention: str = 'radio'`
速度定義です。通常 `radio` を使います。

#### `rest_value_hz: Optional[float] = None`
rest frequency を手動指定する場合の値 [Hz]。

#### `require_kms: bool = False`
`True` なら km/s 変換に失敗した時点で停止します。

### I/O と安定性

#### `file_format: Optional[str] = None`
`spectral-cube` 読み込み形式を明示したい場合に使います。

#### `use_dask: Optional[bool] = None`
Dask 読み込みを使うかどうか。

#### `fill_value: float = np.nan`
NumPy 展開時の fill 値。

#### `overwrite: bool = True`
上書き保存するかどうか。

#### `strip_checksum: bool = True`
既存 HDU の `CHECKSUM/DATASUM` を削除してから保存するかどうか。

### Moment 関連

#### `moment_how: str = 'auto'`
`spectral-cube.moment0()` の `how`。

#### `moment_spectral_slab`
moment0 を積分する速度範囲。

### 注意点

- `convert_to_kms=True` でも、rest frequency が足りないと Hz のまま残る場合があります
- その場合、`rms_spectral_slab` や `moment_spectral_slab` の float は **残っている単位** で解釈されます
- km/s で確実に指定したい場合は `Quantity` を使う方が安全です

---

## 8. FITS 出力 API

## 8.1 `save_map_fits()`

### 役割

OTF / 汎用マッピング結果 `GridResult` を 3D Multi-Extension FITS として保存します。

### シグネチャ

```python
save_map_fits(
    grid_res,
    v_tgt,
    coord_sys,
    projection,
    lon0,
    lat0,
    config,
    out_path,
    out_scale='TA*',
    rep_beameff=1.0,
)
```

### パラメーター

#### `grid_res`
`GridResult`。

#### `v_tgt`
出力速度軸。

#### `coord_sys`
座標系。

#### `projection`
投影法。

#### `lon0`, `lat0`
基準座標。

#### `config`
通常は `MapConfig`。空間 WCS の `x0`, `y0`, `cell_arcsec` に使います。

#### `out_path`
出力パス。

#### `out_scale='TA*'`
温度スケールラベル。

#### `rep_beameff=1.0`
代表的主ビーム効率。

### 出力内容

- Primary HDU: cube（`(nchan, ny, nx)` に転置して保存）
- `WEIGHT`, `HIT`, `MASK`, `TSYS`, `TINT`, `TIME`, `RMS`, `NEFF`, `XEFF`, `YEFF`, `BIAS_PIX`

### 現行設計上の注意

- 元コードでは `BUNIT='K (TA*)'` 形式が残っている版があります
- 更新仕様では、単位そのものは `K` とし、温度スケールは別 keyword で明示する方が望ましいです

---

## 9. basket-weave API

## 9.1 `solve_basket_weave_offsets()`

### 役割

交差点での差を使い、scan ごとの定数オフセットを最小二乗で推定します。

### シグネチャ

```python
solve_basket_weave_offsets(
    input_data: GridInput,
    search_radius_arcsec: float = 3.0,
    damp: float = 1e-2,
    v_axis: np.ndarray | None = None,
    v_windows_kms: list[str] | list[tuple[float, float]] | None = None,
    channel_mask: np.ndarray | None = None,
) -> np.ndarray
```

### パラメーター

#### `input_data: GridInput`
`scan_id` を持つ `GridInput`。

#### `search_radius_arcsec: float = 3.0`
交差点探索半径 [arcsec]。

#### `damp: float = 1e-2`
LSQR の damping 係数。

#### `v_axis: np.ndarray | None = None`
速度軸 [km/s]。`v_windows_kms` を使うときに必要です。

#### `v_windows_kms`
line-free velocity windows。例:

- `['-30:-5', '30:55']`
- `[(-30, -5), (30, 55)]`

#### `channel_mask: np.ndarray | None = None`
すでにチャネルマスクがある場合はそれを使います。

### 戻り値

各 `scan_id` に対応するオフセット配列。

### 注意点

- `scan_id` が無いと実行できません
- `v_windows_kms` を使わない場合、全チャネル平均を使うため、強線があると推定が偏る可能性があります

---

## 9.2 `apply_basket_weave_correction()`

### 役割

推定済みオフセットを `GridInput.spec` に直接適用します。

### シグネチャ

```python
apply_basket_weave_correction(input_data: GridInput, offsets: np.ndarray) -> None
```

### パラメーター

#### `input_data: GridInput`
`scan_id` を含む `GridInput`。

#### `offsets: np.ndarray`
scan ごとの補正量。

### 挙動

各スペクトルの全チャネルから、その scan に対応する定数オフセットを引きます。

---

## 10. WCS / 座標変換 API

## 10.1 `project_to_plane()`

### 役割

天球座標を局所平面 `(x, y)` [arcsec] に投影します。

### シグネチャ

```python
project_to_plane(
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    lon0: float,
    lat0: float,
    projection: str = 'SFL',
    invert_x: bool = True,
) -> tuple[np.ndarray, np.ndarray]
```

### パラメーター

#### `lon_deg`, `lat_deg`
入力座標 [deg]。

#### `lon0`, `lat0`
投影中心 [deg]。

#### `projection: str = 'SFL'`
`CAR` か `SFL` 系。

#### `invert_x: bool = True`
X 軸の正方向を FITS 表示に合わせて左向きにするかどうか。

### 戻り値

`(x_arcsec, y_arcsec)`

### 注意点

- `SFL` では `cos(lat)` を掛けます
- `CAR` では掛けません

---

## 10.2 `build_spatial_wcs_dict()`

### 役割

空間 WCS ヘッダ辞書を構築します。

### シグネチャ

```python
build_spatial_wcs_dict(
    coord_sys: str,
    projection: str,
    lon0: float,
    lat0: float,
    config_x0: float,
    config_y0: float,
    cell_arcsec: float,
    nx: int,
    ny: int,
) -> dict
```

### パラメーター

#### `coord_sys`
座標系。

#### `projection`
投影法。`CAR` か `SFL`。

#### `lon0`, `lat0`
基準座標 [deg]。

#### `config_x0`, `config_y0`
局所平面での基準オフセット [arcsec]。

#### `cell_arcsec`
ピクセルスケール [arcsec/pix]。

#### `nx`, `ny`
画像サイズ。

### 注意点

- `SFL` では `CRVAL2=0` に固定し、`CRPIX2` をずらします
- 不明な projection は無言フォールバックせず `ValueError` です

---

## 11. 推奨ワークフロー

## 11.1 OTF でまず地図を作る

1. `MapConfig` を作る
2. `run_mapping_pipeline()` または `run_otf_full_pipeline()`
3. 出力 FITS の `RMS`, `HIT`, `WEIGHT` を確認

## 11.2 baseline を引く

1. `subtract_baseline_from_fits()` を使う
2. まずは自動 line-free / 自動 ripple で試す
3. 問題があれば `linefree_mask`, `ripple_freqs` を手動指定する
4. `BASE_RMS`, `BASE_FLG`, `LINEFREE`, `RIPFREQ` を確認する

## 11.3 signal mask を作る

1. baseline 後の FITS に対して `make_3d_mask_for_existing_fits()` を使う
2. `rms_ext` は baseline 後なら `BASE_RMS` 優先が自然
3. `method='simple'` で荒く確認し、必要に応じて `smooth_mask` 系へ進む

## 11.4 PS で規則的な格子化をする

1. `PSMapConfig` を用意する
2. `run_ps_mapping_pipeline()`
3. `grid_anchor_offsets` と `grid_tol` を意識して、セル割当が妥当か確認する

---

## 12. 重要な設計上の注意

### 12.1 `MASK3D` と baseline 情報は同じではない

`MASK3D` は、あくまで **最終的に信号があると判断した領域** です。`LINEFREE` の補集合や baseline support をそのまま `MASK3D` と見なすべきではありません。

### 12.2 baseline の line-free は global であることが多い

`estimate_linefree_mask_from_cube()` は cube 全体から 1 本の global mask を作ります。空間位置ごとの差を持たせたい場合は、別設計が必要です。

### 12.3 velocity slab の単位に注意

`spectral-cube` の km/s 変換に失敗すると、slab の float 指定は km/s ではなく元の単位（例えば Hz）のまま解釈されます。安全を重視するなら `Quantity` を使ってください。

### 12.4 `BUNIT` と温度スケールは分けて考える

`TA*`, `TR*` は物理量の種類です。一方 `K` は単位です。実務上は

- `BUNIT='K'`
- `TEMPSCAL='TA*'` または `'TR*'`

のように分けた方が誤解が減ります。

### 12.5 stale HDU を残さない

再処理するときは、古い `MASK3D`, `MOMENT0`, `BASE_RMS`, `BASE_FLG`, `LINEFREE`, `RIPFREQ` などが残ると混乱の元になります。常に「どの HDU が最新か」をはっきりさせてください。

---

## 13. よくある質問

### Q1. `RMS` と `BASE_RMS` のどちらを使えばよいですか。

- baseline 前の地図品質を見るなら `RMS`
- baseline 後の signal mask や baseline 成否を見るなら `BASE_RMS`

### Q2. `simple` と `smooth_mask` のどちらがよいですか。

- まず挙動確認なら `simple`
- 本格解析なら `smooth_mask` か `smooth_mask_lite`
- 巨大 cube なら `smooth_mask_lite`

### Q3. `TR*` にすべきですか。

- 物理解釈でビーム効率補正込みの温度を使いたいなら `TR*`
- 生に近い扱いなら `TA*`

### Q4. `cube_ext` は明示した方がよいですか。

はい。特に Primary 以外に cube がある FITS では明示した方が安全です。

### Q5. basket-weave は必ず必要ですか。

走査方向ごとに定数オフセットが出る場合は有効ですが、常に必要ではありません。交差点が少ない場合や line-free 範囲が不適切な場合は不安定になります。

---

## 14. 最後に

このパッケージは、

- グリッディング
- baseline 除去
- signal mask 生成
- moment 解析

を 1 つの流れで扱えるように設計されています。ただし、それぞれの段階は**役割が違う**ので、

- グリッディングの診断量
- baseline の QC
- final signal mask

を混ぜずに扱うことが重要です。

本説明書は、`map_3d` 系の公開 API と主要パラメーターを、省略を減らして整理したものです。実際の運用では、この説明書と cookbook を併用すると理解しやすくなります。


---

## 15. OTF 観測の実践例: Orion KL の `^{12}CO` キューブを作って baseline を引く

この節では、実際に使われている OTF 解析スクリプトを、`sd_radio_spectral_fits.map_3d` の考え方に沿って整理し直して説明します。

### 15.1 この例でやりたいこと

生データ `sc_raw` から始めて、次の順に処理します。

1. chopper wheel calibration で `Ta*` の scantable を作る
2. basket-weave 補正で scan ごとのオフセットを補正する
3. OTF グリッディングで 3D FITS cube を作る
4. cube に対して baseline fitting を行う
5. 必要なら `MASK3D` と `MOMENT0` を作る

### 15.2 元の処理の意味

ユーザー例の処理は、概念的には次の 4 段階に分かれます。

- **較正段階**: `run_tastar_calibration()`
- **OTF 特有の補正**: `create_grid_input()` → `solve_basket_weave_offsets()` → `apply_basket_weave_correction()`
- **キューブ化**: `run_mapping_pipeline()`（環境によっては `run_otf_pipeline()` というラッパー名を使っている場合があります）
- **cube baseline 補正**: `subtract_baseline_from_fits()`

### 15.3 整理した実用コード例

```python
import numpy as np
from astropy.io import fits

from sd_radio_spectral_fits import calibration as cal
from sd_radio_spectral_fits.map_3d.gridder import create_grid_input, run_mapping_pipeline
from sd_radio_spectral_fits.map_3d.basketweave import (
    solve_basket_weave_offsets,
    apply_basket_weave_correction,
)
from sd_radio_spectral_fits.map_3d.config import GridConfig
from sd_radio_spectral_fits.map_3d.baseline_subtraction import (
    subtract_baseline_from_fits,
    BaselineConfig,
    RippleConfig,
)

# ------------------------------------------------------------
# 1. Chopper wheel calibration
# ------------------------------------------------------------
sc_cal_12co = cal.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-30.0, 55.0),
    t_hot_k=300.0,
    vcorr_chunk_sec=10.0,
    dtype="float32",
)

# ------------------------------------------------------------
# 2. Basket-weave correction
# ------------------------------------------------------------
# projection は後段のマッピングとそろえるのが分かりやすいです
# ここでは例として CAR を使います。

tmp_in = create_grid_input(sc_cal_12co, projection="CAR")
offsets = solve_basket_weave_offsets(
    tmp_in,
    search_radius_arcsec=35.0,
)
apply_basket_weave_correction(tmp_in, offsets)

# 補正後のスペクトル行列を scantable 側へ戻す
sc_cal_12co.data = tmp_in.spec

# ------------------------------------------------------------
# 3. OTF グリッディングして 3D FITS を作る
# ------------------------------------------------------------
config = GridConfig(
    x0=-3600.0,
    y0=-3600.0,
    nx=201,
    ny=201,
    cell_arcsec=36.0,
    beam_fwhm_arcsec=350.0,
    kernel="gjinc",
    truncate="first_null",
)

cube_fits = "ori-kl-12co.fits"

res = run_mapping_pipeline(
    scantable=sc_cal_12co,
    config=config,
    output_fits=cube_fits,
    projection="CAR",
    ref_lon=83.809,
    ref_lat=-5.372639,
    dv_kms=0.2,
    out_scale="TA*",
)

# ------------------------------------------------------------
# 4. cube baseline fitting
# ------------------------------------------------------------
with fits.open(cube_fits) as hdul:
    nchan = hdul[0].data.shape[0]  # (nchan, ny, nx) を想定

lf_mask = np.ones(nchan, dtype=bool)
lf_mask[169:215] = False

r_cfg = RippleConfig(
    nfreq=6,
    period_range_chan=(10.0, 1000.0),
)

b_cfg = BaselineConfig(
    poly_order=1,
    ripple=True,
    robust=True,
)

cube_bsub_fits = "ori-kl-12co_bsub.fits"

subtract_baseline_from_fits(
    input_fits=cube_fits,
    output_fits=cube_bsub_fits,
    linefree_mask=lf_mask,
    baseline_cfg=b_cfg,
    ripple_cfg=r_cfg,
    add_qc_hdus=True,
    overwrite=True,
)
```

### 15.4 各段階で何をしているか

#### 15.4.1 `run_tastar_calibration()`

この段階では、生のスペクトル列から `Ta*` スケールの scantable を作ります。

- `vlsrk_range_kms=(-30, 55)`
  - 速度補正や出力軸の想定範囲です。
  - Orion KL の `^{12}CO` で、信号と基線の両方を十分含む範囲として妥当です。
- `t_hot_k=300`
  - HOT load の温度です。
- `vcorr_chunk_sec=10`
  - 速度補正を何秒ごとにまとめて評価するかです。
- `dtype='float32'`
  - メモリ節約と通常の解析精度のバランスがよい設定です。

#### 15.4.2 `create_grid_input()` と basket-weave 補正

ここでは scantable を、グリッディングエンジンが扱いやすい `GridInput` に直し、その上で scan ごとのオフセットを推定します。

- `projection='CAR'`
  - 座標の局所平面投影法です。
  - この例では最終的なマッピングも `CAR` なので、途中もそろえておくと理解しやすくなります。
- `search_radius_arcsec=35.0`
  - 交差点を探す半径です。
  - 小さすぎると交差点が減り、大きすぎると無関係な点を結びやすくなります。

`apply_basket_weave_correction()` は `tmp_in.spec` をその場で書き換えます。そのため、最後に `sc_cal_12co.data = tmp_in.spec` として scantable 本体へ反映しています。

#### 15.4.3 `GridConfig`

この設定は、「どの範囲を、どのピクセルサイズで、どのカーネルで画像化するか」を決めます。

- `x0=-3600.0`, `y0=-3600.0`
  - マップ左下に相当する開始オフセットです。
- `nx=201`, `ny=201`
  - ピクセル数です。
  - この例では 201 × 201 ピクセルの地図になります。
- `cell_arcsec=36.0`
  - 1 ピクセルの大きさです。
- `beam_fwhm_arcsec=350.0`
  - 望遠鏡ビームの FWHM です。
- `kernel='gjinc'`
  - 典型的な single-dish OTF に向くカーネルです。
- `truncate='first_null'`
  - GJINC カーネルの打ち切り半径です。

#### 15.4.4 `run_mapping_pipeline()`

この関数は、速度軸再サンプル、座標投影、温度スケール処理、OTF グリッディング、FITS 書き出しを一括して行います。

- `projection='CAR'`
  - 局所平面投影法です。
- `ref_lon=83.809`, `ref_lat=-5.372639`
  - 地図の基準位置です。
  - Orion KL の中心座標を手で与えています。
- `dv_kms=0.2`
  - 出力 cube の速度分解能です。
- `out_scale='TA*'`
  - 出力温度スケールです。

この段階の出力 `ori-kl-12co.fits` には、通常 Primary HDU の 3D cube に加えて、`HIT`、`RMS`、`TSYS`、`TINT` などの 2D HDU が付きます。

#### 15.4.5 `subtract_baseline_from_fits()`

この段階では、**scantable ではなく、すでに作った 3D FITS cube** に baseline fitting を行います。

今回の例では、line-free 領域を手動で与えています。

- `lf_mask = np.ones(nchan, dtype=bool)`
  - まずは全チャネルを line-free 候補にする
- `lf_mask[169:215] = False`
  - 本物の輝線がある範囲だけ baseline fitting から除外する

`RippleConfig` では、定在波成分の探索方法を決めます。

- `nfreq=6`
  - 抽出する ripple 周波数の数です。
- `period_range_chan=(10, 1000)`
  - 探索する周期範囲です。
  - 複雑な ripple を疑うなら広めに取ります。

`BaselineConfig` では、baseline モデルそのものを決めます。

- `poly_order=1`
  - 1 次多項式です。
- `ripple=True`
  - 多項式に加えて sinusoid 成分も使います。
- `robust=True`
  - 外れ値に引きずられにくい fitting にします。

この処理の出力 `ori-kl-12co_bsub.fits` には、cube 本体が baseline 補正後の値で保存され、必要に応じて `LINEFREE`、`RIPFREQ`、`BASE_RMS`、`BASE_FLG` などの QC HDU が付きます。

### 15.5 この例のあとに何をすればよいか

次の標準的な流れは、baseline 後 cube に対して `MASK3D` と `MOMENT0` を作ることです。

```python
from sd_radio_spectral_fits.map_3d.cube_analysis import make_3d_mask_for_existing_fits

make_3d_mask_for_existing_fits(
    input_fits="ori-kl-12co_bsub.fits",
    output_fits="ori-kl-12co_bsub_masked.fits",
    rms_ext="auto",
    method="smooth_mask_lite",
    moment_spectral_slab=(-20.0, 40.0),
    overwrite=True,
)
```

ここで `rms_ext='auto'` にしておくと、baseline 後 cube では `BASE_RMS` を優先して使い、無ければ `RMS`、さらに無ければ cube から robust RMS を再推定します。

### 15.6 この例で特に重要な注意

1. **`run_otf_pipeline` という名前を使っている環境もある**
   - 提供されたソースでは、高水準関数は主に `run_mapping_pipeline()` と `run_otf_full_pipeline()` です。
   - お手元で `run_otf_pipeline()` が生きているなら、そのラッパーを使って構いません。
   - ただし概念的には、ここで示した流れと同じです。

2. **baseline は scantable に対してではなく cube に対して行っている**
   - この例では、まず OTF cube を作ってから `subtract_baseline_from_fits()` をかけています。
   - したがって baseline の line-free mask は「チャネル番号ベース」で与えています。

3. **line-free の手動指定は最初の 1 回としては非常に有効**
   - 強い線・広い翼・複雑な ripple がある領域では、自動推定よりも手動の方が安定することがあります。
   - ただし、深い積分で弱い線が後から見えてくることもあるので、必要なら後日 `LINEFREE` を見直します。

4. **basket-weave を行うなら、その後の cube baseline と役割を混同しない**
   - basket-weave は scan ごとのオフセット補正です。
   - cube baseline はチャネル方向の基線補正です。
   - 目的が違うので、両方とも有効な場合があります。
