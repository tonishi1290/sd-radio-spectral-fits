# baseline 自動検出詳細ノート（3D auto line-free）

## 0. この文書の位置づけ

この文書は、`baseline_subtraction_detailed_manual_release_candidate_2026-03-27_rev2.md` の補助資料である。  
主マニュアルでは運用・入出力・代表的な使い方を中心に説明し、本書では **3D auto line-free の設計思想・内部フロー・主要パラメータの意味・調整指針・非保証条件** を詳しく説明する。

本書の対象は主に以下の機能である。

- `mask_kind='voxel_3d'`
- `auto_method='local_3d_poly'`
- `auto_mask_profile='default'|'no_ripple'`
- `safe_velocity_windows_kms`
- `ripple_apply_stage='final'|'mask_only'|'both'`
- `linefree_detection_data`

本書では、公開向けパラメータ名（`arcsec` / `km/s` / `*_threshold_sigma`）を優先して説明する。内部実装では最終的に pixel / channel 単位へ正規化して処理する。

---

## 1. 目的と基本方針

### 1.1 目的

3D auto line-free の目的は、入力 cube

\[
D(v,y,x)
\]

に対して、**各 spatial pixel ごとに異なる line-free 領域** を推定し、最終 baseline fit に使うための 3D boolean mask

\[
L(v,y,x)
\]

を作ることである。ここで

- `True`: baseline fit に使ってよい voxel
- `False`: line 候補、または baseline fit から除外すべき voxel

を表す。

### 1.2 役割分担

本実装では、責務を以下の 3 段に明確に分ける。

#### Stage A: detection cube と 3D local hard seed 生成

- まず、auto 判定専用の **detection cube** を作る
- detection cube では、必要に応じて **x,y 方向 smoothing** と **v 方向 smoothing** を使う
- その detection cube 上で、強い line の芯や明瞭な wing を **hard mask** として保護する
- ここで作るのは **3D local seed** であり、global 1D seed は使わない

#### Stage B: low-capacity provisional fit による mask 更新

- Stage A の seed を初期条件として、低自由度モデルで provisional baseline を作る
- その残差から line 候補を更新していく
- ここでは **mask を作ること** が目的であり、最終 science 用 baseline を決める段階ではない

#### Stage C: 最終 baseline subtraction

- Stage B で得た `L(v,y,x)` を使い、元の cube に対して最終 fit を行う
- ripple 項や高次モデルを使うのはこの段階

### 1.3 なぜこの分離が必要か

この分離は、以下の 2 つの失敗モードを避けるために重要である。

1. **強い line に baseline が持ち上げられる問題**  
   強 line を soft weight だけで扱うと、L2 最小二乗が line に引っ張られて baseline が上に歪みやすい。

2. **広い欠損 + 高自由度モデルの発散**  
   広い hard mask を掛けた状態で ripple 項を含む高表現力モデルを動かすと、masked 領域で発散しやすい。

したがって、

- Stage A/B では **低自由度・剛性重視**
- Stage C では **最終品質重視**

とする。

---

## 2. 数学的定義

### 2.1 入力と出力

- 入力 cube:  \( D(v,y,x) \)
- detection cube: \( D_{\rm det}(v,y,x) \)
- seed line mask: \( S(v,y,x) \)
- provisional baseline: \( B_{\rm prov}(v,y,x) \)
- 最終 baseline: \( B(v,y,x) \)
- residual: \( R(v,y,x)=D-B \)
- line-free mask: \( L(v,y,x) \)

### 2.2 line-free / signal の意味

`LINEFREE3D_USED` は 3D voxel mask であり、`True` が baseline fit に使われた voxel を表す。

`SIGNAL_MASK3D_USED` は各 `(y,x)` ごとに channel 軸方向の **internal gap** を取り、signal 用 voxel を表す。

概念的には、各 pixel で

- `LINEFREE3D_USED[:,y,x]` = baseline fit に使用した channel
- その complement 全体ではなく、**line-free の内部 gap** を `SIGNAL_MASK3D_USED[:,y,x]` とする

という整理である。

`MOMENT0` は、その run における **正準な integrated intensity** とし、別名にはしない。

---

## 3. detection cube とは何か

### 3.1 detection cube の考え方

本実装では、最終 fit に使う元の cube \(D\) とは別に、auto 判定専用の cube

\[
D_{\rm det}(v,y,x)
\]

を作る。

この detection cube の目的は、

- 近傍 pixel の情報を借りて S/N を上げる
- 高周波ノイズを少し抑える
- weak line や velocity gradient を検出しやすくする
- ripple が極端に強くない限り、line 判定を安定化する

ことである。

### 3.2 基本形

概念的には

\[
D_{\rm det} = S_v\bigl(S_{xy}(D)\bigr)
\]

である。

ここで

- \(S_{xy}\): x,y 方向の smoothing
- \(S_v\): v 方向の smoothing

である。

### 3.3 重要な原則

- detection cube は **auto 判定専用** である
- Stage B / Stage C の fit は、原則として **元の unsmoothed cube** \(D\) に対して行う
- したがって、smoothing は auto 判定を安定化するための補助手段であり、最終 science cube を改変するものではない

---

## 4. Stage A: detection cube と 3D local hard seed

### 4.1 Stage A の目的

Stage A の目的は、Stage B の初手 provisional fit を強 line に汚染させないことである。  
ここでの出力は 3D boolean mask

\[
S(v,y,x)
\]

であり、`True` が **初手から line と見なして fit から除外する voxel** を表す。

### 4.2 Stage A の入力データ

Stage A では detection cube \(D_{\rm det}\) を使う。  
現実の実装では、以下の 2 種類の smoothing を使える。

#### (A) x,y 方向 smoothing

目的:
- 近傍 pixel の情報を借りて局所 S/N を上げる
- velocity gradient や weak extended emission の検出を安定化する

公開パラメータ:
- `lf3d_detect_spatial_fwhm_arcsec`
- `lf3d_detect_spatial_fwhm_pix`

内部では Gaussian kernel に変換して使う。  
無効値を含む場合でも bias を出しにくいように、**valid-weight 正規化付き Gaussian smoothing** を行う。

#### (B) v 方向 smoothing

目的:
- 高周波ノイズを少し抑える
- narrow spike に引っ張られにくくする
- seed 生成前の detection を安定化する

公開パラメータ:
- `lf3d_detect_spectral_width_kms`
- `lf3d_detect_spectral_width_chan`

内部では **boxcar smoothing** を用いる。  
これは Stage A の検出補助であり、最終 fit には使わない。

### 4.3 seed baseline

Stage A の baseline は、最終 baseline ではない。  
目的は **line を粗く見つけるための「底」** を作ることである。

現行の公開候補では、Stage A の baseline 候補として主に以下を使う。

- `lf3d_seed_baseline_mode='median'`
- （参考実装）`'opening'`

現時点の推奨は `median` である。  
`opening` は一部ケースには効くが、Case1/4/5 を壊しやすかったため、既定値にはしない。

### 4.4 seed 閾値

以前の `lf3d_seed_quantile` は global 1D seed 由来であり、local 3D seed では物理的意味が薄い。  
現在は **局所ノイズに対する直接閾値** として

- `lf3d_seed_threshold_sigma`

を用いる。

seed line 判定の概念は

\[
S_0(v,y,x) =
\Bigl(D_{\rm det}(v,y,x) - B_{\rm seed}(v,y,x)
> k_{\rm seed}\,\sigma_{\rm loc}(y,x)\Bigr)
\]

である。

### 4.5 seed hysteresis / dilation

Stage A では、狭い強ピークだけでなく、その wing を取りこぼさないことが重要である。  
そのために

- `lf3d_seed_hysteresis_threshold_sigma`
- `lf3d_seed_dilation_kms`
- `lf3d_seed_dilation_chan`

を用いる。

ただし、固定幅 dilation は

- 狭線には過剰
- 広い wing には不足

になりやすいので、必要なら seed hysteresis を優先する。

### 4.6 Stage A の主要パラメータと意味

#### `lf3d_detect_spatial_fwhm_arcsec` / `lf3d_detect_spatial_fwhm_pix`

- 役割: detection cube を作る前段の **x,y 方向 smoothing 幅**
- 大きくすると:
  - 近傍 pixel の line 情報を借りやすくなる
  - velocity gradient や弱い extended emission を検出しやすい
  - 一方で spatial に孤立した構造や sharp な境界はぼける
- 小さくすると:
  - local 性が上がる
  - 強い孤立 line には有利
  - ただし低 S/N や velocity gradient では取りこぼしやすい

**実務的目安**
- `0` または `None`: x,y smoothing なし
- `0.5 beam 前後`: 軽い smoothing
- `1 beam 前後`: 強めの smoothing

cell size が細かい場合は `arcsec` 指定が推奨である。

#### `lf3d_detect_spectral_width_kms` / `lf3d_detect_spectral_width_chan`

- 役割: detection cube を作る前段の **v 方向 smoothing 幅**
- これは final fit の spectral smoothing ではなく、**Stage A の line 検出を安定化するためだけの弱い smoothing** である
- 大きくすると:
  - 高周波ノイズや narrow spike を抑えやすい
  - weak line の芯が安定する
  - ただし近接線をつなげたり、狭線を太らせる危険がある
- 小さくすると:
  - local な spectral structure を保ちやすい
  - ただし noisy な detection になりやすい

**実務的目安**
- `0` または `1 chan`: 実質的に smoothing なし
- `数 channel` または `~0.2–0.5 km/s`: 軽い smoothing
- あまり大きくしすぎない方が安全

#### `lf3d_seed_spatial_fwhm_arcsec` / `lf3d_seed_spatial_fwhm_pix`

- 役割: Stage A の seed 生成そのものに使う **x,y 方向 smoothing 幅**
- `lf3d_detect_*` と役割が近いが、実装や将来の整理では detection cube 生成と seed 生成を分けて調整するために独立に持てるようにしている
- 現状の主眼は、seed 生成前の spatial 情報の混合量をどれだけ与えるかである

#### `lf3d_seed_median_width_kms` / `lf3d_seed_median_width_chan`

- 役割: Stage A の `median` baseline で使う **v 方向の broad trend 幅**
- これは「速度方向 smoothing」に最も近いパラメータの 1 つである
- 大きくすると:
  - broad ripple や長周期の揺れを baseline 側へ落としやすい
  - ただし broad line も baseline に食わせやすい
- 小さくすると:
  - narrow line を seed しやすい
  - ただし ripple / broad undulation を line と誤認しやすい

**実務的目安**
- 狭線・攻めた検出: やや小さめ
- ripple や broad baseline を恐れる: やや大きめ

#### `lf3d_seed_threshold_sigma`

- 役割: seed を hard mask に入れる強さの閾値
- 単位: 局所ノイズ \(\sigma\) の倍数
- 大きくすると:
  - 保守的になる
  - ripple や弱い揺らぎを line 扱いしにくい
  - ただし weak line / wing を取りこぼしやすい
- 小さくすると:
  - line は拾いやすい
  - ただし false seed が増える

#### `lf3d_seed_hysteresis_threshold_sigma`

- 役割: strong seed の周囲を weak レベルまで追跡して hard seed を広げる
- これは **v 方向で line の裾野をどこまで hard seed に含めるか** を決める
- `None` にすると Stage A では hysteresis を使わない
- 小さくすると seed が広がりやすい

#### `lf3d_seed_dilation_kms` / `lf3d_seed_dilation_chan`

- 役割: no-wrap の固定 spectral dilation 幅
- `lf3d_seed_hysteresis_threshold_sigma` を使わない場合の簡単な widen 手段
- ただし固定幅なので、広い wing と狭い line を同時に最適には扱えない

### 4.7 Stage A の調整指針

#### velocity gradient を取りたい
- `lf3d_detect_spatial_fwhm_arcsec` と `lf3d_seed_spatial_fwhm_arcsec` をやや大きめ
- `lf3d_detect_spectral_width_kms` は軽く入れる
- `lf3d_seed_threshold_sigma` を少し下げる
- `lf3d_seed_hysteresis_threshold_sigma` を有効にする

#### ripple が気になる
- `lf3d_seed_median_width_kms` をやや広め
- `lf3d_seed_threshold_sigma` を高め
- `lf3d_seed_hysteresis_threshold_sigma` は切るか高めにする
- 必要なら `safe_velocity_windows_kms` と detection-only cube を使う

#### simple な強 line を守りたい
- `lf3d_seed_threshold_sigma` を必要以上に高くしない
- `lf3d_seed_hysteresis_threshold_sigma` を弱めに使う
- `lf3d_seed_median_width_kms` を広げすぎない

---

## 5. Stage B: provisional fit による 3D mask 更新

### 5.1 Stage B の目的

Stage B の目的は、Stage A の hard seed を出発点として、より完全な 3D line-free mask

\[
L(v,y,x)
\]

を得ることである。

ここでは **mask を作ること** が目的であり、science 用の最終 baseline を決める段階ではない。

### 5.2 provisional model は低自由度に限定する

Stage B では ripple 項を使わない。  
provisional fit は **polynomial only** に固定し、次数も最終 fit 用とは分離する。

現行では

- `lf3d_prov_poly_order`

を使う。推奨は `0` または `1` である。

`bcfg.poly_order` をそのまま Stage B に流し込まないことが重要である。

### 5.3 反復の基本式

初期 line mask を `line_accum = S` とする。  
各反復で、fit に使う重みは概念的に

\[
W = \mathrm{finite}(D) \cap \neg \mathrm{line\_accum}
\]

である。

この `W` を使って provisional baseline

\[
B_{\rm prov}
\]

を作り、残差

\[
R_{\rm prov}=D-B_{\rm prov}
\]

から strong / weak 判定を行う。

### 5.4 strong / weak 判定

現行の改修では

- **strong 判定**: unsmoothed residual
- **weak 判定**: spatially smoothed positive residual

を使う方向が最もバランスが良かった。

これに対応する主なパラメータは

- `lf3d_detect_threshold_sigma`
- `lf3d_hysteresis_threshold_sigma`
- `lf3d_weak_spatial_fwhm_arcsec`
- `lf3d_weak_spatial_fwhm_pix`

である。

この設計は、

- strong line の芯は local に確実に捉える
- weak extended や velocity gradient に沿う成分は spatial coherence を借りて拾う

ことを狙っている。

### 5.5 monotone update

Stage B では、line 領域は基本的に OR で累積し、単調に増やす。

すなわち概念的には

\[
\mathrm{line\_accum}_{n+1}
=
\mathrm{line\_accum}_n \,\cup\, \mathrm{new\_line}_n
\]

である。

これは、Stage B の model を低自由度・剛性重視にしているため、安全側に倒しやすいからである。

### 5.6 Stage B の主要パラメータと意味

#### `lf3d_prov_poly_order`

- 役割: provisional fit の硬さを決める
- `0`: 最も stiff
- `1`: 傾きまで吸収できる
- 高くしすぎると、mask 推定段階で line を baseline 側へ食いやすい

したがって Stage B では通常 `0` または `1` を推奨する。

#### `lf3d_detect_threshold_sigma`

- 役割: strong line 判定の基本閾値
- 単位: 局所ノイズ \(\sigma\) の倍数
- 小さくすると line 回収は強くなるが、false line も増える
- 大きくすると conservative になる

#### `lf3d_hysteresis_threshold_sigma`

- 役割: strong 判定で見つかった line の周囲を、weak レベルまで追跡して line_accum へ入れる
- 速度方向の「裾野の広がり」を制御する重要パラメータ
- `None` の場合、Stage B の hysteresis を使わない

#### `lf3d_weak_spatial_fwhm_arcsec` / `lf3d_weak_spatial_fwhm_pix`

- 役割: weak 判定のためだけに positive residual を x,y 方向に smoothing する幅
- 大きくすると:
  - weak extended や velocity gradient に強くなる
  - ただし ripple や広い偽構造も weak 扱いしやすくなる
- 小さくすると:
  - local な判定に寄る
  - ripple には安全だが、extended emission を拾いにくい

### 5.7 Stage B の調整指針

#### weak extended を拾いたい
- `lf3d_weak_spatial_fwhm_arcsec` をやや大きめ
- `lf3d_hysteresis_threshold_sigma` を有効化
- `lf3d_detect_threshold_sigma` は少し下げる

#### ripple を恐れる
- `lf3d_weak_spatial_fwhm_arcsec` を控えめにする
- `lf3d_hysteresis_threshold_sigma` を切るか高めにする
- Stage A の seed 側で broad median を強める

#### 速度勾配を追いたい
- `lf3d_weak_spatial_fwhm_arcsec` を `1.0–1.5 beam` 程度で試す
- Stage A の spatial smoothing とセットで見る

---

## 6. Stage C: 最終 baseline subtraction

Stage C では、Stage B で得られた `L(v,y,x)` を用いて、元の cube に対して最終 fit を行う。

ここで初めて

- 多項式
- ripple 項
- robust refit

を使う。

### 6.1 solver

3D voxel path の既定 solver は QR ベースである。

- `voxel_solver='qr'`

を既定とし、`normal` は optional に留める。

### 6.2 grouping

identical な line-free mask pattern が多い場合は grouping によって高速化できる。  
ただし、実データで本当に効くかは mask の繰り返し度合いに依存する。

そのため

- `voxel_grouping='auto'|'always'|'off'`

を持つ。

### 6.3 robust

3D path の robust は、全 pixel へ一律に掛けるのではなく、**selective robust** を基本とする。

---

## 7. ripple との責務分離

### 7.1 auto 3D line-free は ripple 専用解決器ではない

auto 3D line-free の主責務は、**正の line を fit から外すこと**である。  
強い ripple の存在下では、Stage A/B が ripple crest を line と誤認しやすい。

したがって本実装では

- auto line-free は ripple 専用モデルを持たない
- 強い ripple の下での auto 精度は保証しない

という立場を取る。

### 7.2 それでも ripple がある場合

強い ripple があり、かつ auto line-free を使いたい場合は、以下を推奨する。

1. **安全な速度範囲** `safe_velocity_windows_kms` を与える
2. その範囲だけで既存 ripple model を fit する
3. その結果を使って **mask 同定専用 cube** を作る
4. auto line-free はその cube で実行する
5. 最終 fit は元の cube に対して行う

つまり、前処理済み cube は **mask 同定専用**であり、science 用最終 cube ではない。

### 7.3 detection-only cube

必要に応じて、以下の入力を与えられるようにする。

- `linefree_detection_data`

意味は、

- auto mask 生成だけに使う cube
- 最終 baseline fit 対象は別に元の cube

である。

---

## 8. ripple parameter の単位系

ユーザー視点では ripple 周期を km/s で読むのが分かりやすい場合がある。  
一方、物理的には ripple は周波数周期である。

そこで本実装では、ripple の探索範囲として

- `period_range_chan`
- `period_range_hz`
- `period_range_kms`

のいずれかを受けられるようにする。

### 8.1 方針

- ユーザー入力: `chan / hz / kms` を許す
- 内部正規化: 必要に応じて Hz へ寄せる
- 実際の FFT / fit: channel 単位へ変換して使う

### 8.2 注意

- `period_range_chan` は sampling 依存であり、物理量そのものではない
- `period_range_kms` は見やすいが、基準周波数に依存する
- `period_range_hz` は物理的には最も自然

したがって、ユーザーには km/s も許すが、文書上は「内部では周波数・channel に正規化して使う」と明示する。

---

## 9. `auto_mask_profile`

現行の公開候補では、mask 側のパラメータ組として

- `default`
- `no_ripple`

を持つ。

### 9.1 `default`

一般用途向け。  
ripple が完全に無いとは仮定しないため、やや保守的。

### 9.2 `no_ripple`

ripple が十分弱い、または事前軽減済みであることを前提に、より攻めた seed / weak 判定を使う。

ただし synthetic では

- Case2 velocity gradient は改善
- Case1 simple, Case4 weak_extended, Case5 ripple は悪化

となりやすかったため、現時点では **experimental** とする。  
既定値は `default` である。

---

## 10. 主要パラメータ一覧（3D auto line-free）

### Stage A / detection cube

- `lf3d_detect_spatial_fwhm_arcsec` / `lf3d_detect_spatial_fwhm_pix`  
  detection cube 作成時の x,y smoothing 幅。
- `lf3d_detect_spectral_width_kms` / `lf3d_detect_spectral_width_chan`  
  detection cube 作成時の v smoothing 幅。
- `lf3d_seed_baseline_mode`  
  `median` / `opening`。既定は `median`。
- `lf3d_seed_threshold_sigma`  
  seed hard mask の強さ。
- `lf3d_seed_spatial_fwhm_arcsec` / `lf3d_seed_spatial_fwhm_pix`  
  Stage A seed 生成側の x,y smoothing 幅。
- `lf3d_seed_median_width_kms` / `lf3d_seed_median_width_chan`  
  Stage A の spectral median baseline 幅。
- `lf3d_seed_baseline_order`  
  代替 baseline mode 用の次数。通常は低次のみ。
- `lf3d_seed_hysteresis_threshold_sigma`  
  Stage A の seed hysteresis 閾値。
- `lf3d_seed_dilation_kms` / `lf3d_seed_dilation_chan`  
  Stage A の no-wrap fixed spectral dilation 幅。

### Stage B

- `lf3d_prov_poly_order`  
  provisional fit の polynomial order。
- `lf3d_detect_threshold_sigma`  
  strong line 判定閾値。
- `lf3d_hysteresis_threshold_sigma`  
  Stage B の weak expansion 閾値。
- `lf3d_weak_spatial_fwhm_arcsec` / `lf3d_weak_spatial_fwhm_pix`  
  weak 判定用 residual の x,y smoothing 幅。
- `lf3d_max_iter`  
  Stage B の反復回数。

### Stage C / fitting

- `voxel_solver`  
  `qr` を既定とする。
- `voxel_grouping`  
  repeated mask pattern を group して高速化。
- `voxel_qr_batch_pix`  
  QR を何 pixel ずつ batched に処理するか。
- `robust`, `robust_mode`  
  final fit の robust refit 設定。

### ripple / detection-only

- `safe_velocity_windows_kms`  
  ripple fit や detection-only cube 生成の安全窓。
- `ripple_apply_stage`  
  `final` / `mask_only` / `both`。
- `period_range_chan`, `period_range_hz`, `period_range_kms`  
  ripple の探索範囲。
- `linefree_detection_data`  
  mask 同定専用 cube。

---

## 11. 非保証条件

以下の条件では、auto 3D line-free の結果をそのまま強く信頼しない方がよい。

1. 強い standing-wave ripple があり、事前軽減や safe window が無い
2. broad weak emission が「安全窓」に入り込んでいる
3. 極端に広い line で、line-free channel がほとんど残らない
4. provisional fit に高次 polynomial や ripple 項を入れてしまった場合

---

## 12. 推奨の使い分け

### 12.1 通常の推奨

- `mask_kind='voxel_3d'`
- `auto_method='local_3d_poly'`
- `auto_mask_profile='default'`
- `voxel_solver='qr'`

### 12.2 ripple が強い場合

- auto line-free を過信しない
- `safe_velocity_windows_kms` を与える
- 必要なら `linefree_detection_data` に user-side 前処理 cube を与える
- final fit は元の cube で行う

### 12.3 ripple が十分弱い場合

- `auto_mask_profile='no_ripple'` を試せる
- ただし experimental であり、実データ比較を推奨する

---

## 13. 今後の改善候補

現時点で主に改善余地があるのは、**mask 同定側**である。

候補としては

- Stage A detection cube の smoothing の最適化
- Stage A seed baseline の改良
- Stage B weak update の spatial coherence 改良
- detection-only cube を使った強 ripple ケースの実データ検証

がある。

fitting 側は、QR solver と grouping を含め、かなり整理されている。

---

## 14. 主マニュアルとの関係

運用・入出力・実行例・bundle/FITS への書き戻し・viewer の使い方は、主に

- `baseline_subtraction_detailed_manual_release_candidate_2026-03-27_rev2.md`

を参照すること。
