# baseline_subtraction 詳細マニュアル

更新日: 2026-03-27
対象: `sd_radio_spectral_fits.map_3d.baseline_subtraction`
対象実装: `linefree_mode='manual'` と `linefree_mode='infer'` を含む 2026-03-27 改訂版

---

# 0. この版で特に重要な変更

今回の改訂で、`line-free` の決め方を次のように整理した。

1. `linefree_cfg` の既定は `None`
2. `linefree_mode` の既定は `'infer'`
3. `linefree_mode='manual'` を追加
4. `linefree_velocity_windows_kms` は
   - `linefree_mode='manual'` のときは **その窓だけで厳密に line-free を決める**
   - 自動探索系 (`auto`, `or`, `prior`, `infer` でそれらに解決された場合) では **最後に OR で追加する include 窓** として使う
5. `exclude_v_windows` は、どのモードでも **最後に除外** する
6. `LineFreeConfig` を与えていないのに、暗黙に auto detection が走ることはない
7. `ripple` 関係の設定について、**最終 fit モデル / 周波数の決め方 / どの段階で使うか** を分離して説明した
8. 本体 API の直引数と `BaselineConfig` / `LineFreeConfig` / `RippleConfig` の所属関係を総覧表として追加した

この文書では、特に `linefree_mask` / `linefree_velocity_windows_kms` / `linefree_cfg` / `linefree_mode` の関係を曖昧さなく定義し、加えて `BaselineConfig.ripple` / `ripple_mode` / `ripple_apply_stage` / `safe_velocity_windows_kms` の責務も分離して定義する。

---

# 1. 先に結論: 何をどう書けばよいか

## 1.1 速度窓を厳密に指定して、その窓だけで baseline fit したい

```python
import sd_radio_spectral_fits.map_3d.baseline_subtraction as m3d_bsl

b_cfg = m3d_bsl.BaselineConfig(
    poly_order=1,
    ripple=False,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=b_cfg,
    add_qc_hdus=True,
    update_mosaic_products=True,
    gain_min=0.5,
)
```

意味:
- `linefree_mode='manual'`
  - 自動探索も prior も使わない
- `linefree_velocity_windows_kms=['-30:0', '20:55']`
  - **この窓だけを line-free として使う**

このケースでは、`linefree_velocity_windows_kms` は include 補助ではなく、**厳密な manual 指定**である。

---

## 1.2 自動探索を使うが、ここは line-free に足したい

```python
lf_cfg = m3d_bsl.LineFreeConfig(
    smooth_width=51,
    sigma=4.0,
    iters=6,
    pad_chan=3,
    min_linefree_frac=0.35,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_cfg=lf_cfg,
    linefree_mode='auto',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

意味:
- auto でまず line-free を推定する
- その後 `['-30:0', '20:55']` を **OR で追加する**

このケースでは、`linefree_velocity_windows_kms` は **manual exact 指定ではない**。

---

## 1.3 prior LINEFREE を使いたい

```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='prior',
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

意味:
- 入力 bundle/FITS に既にある `LINEFREE` または `LINEFREE3D` 系を使う
- prior が存在しなければ error

---

## 1.4 prior と auto の和集合を使いたい

```python
lf_cfg = m3d_bsl.LineFreeConfig(
    smooth_width=51,
    sigma=4.0,
    iters=6,
    pad_chan=3,
    min_linefree_frac=0.35,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_cfg=lf_cfg,
    linefree_mode='or',
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

意味:
- base mask を `prior ∪ auto` とする
- その後、`linefree_velocity_windows_kms` が与えられていればさらに OR で追加
- `exclude_v_windows` が与えられていれば最後に除外

---

## 1.5 `linefree_mask` を直接与えたい

```python
lf = ...  # shape (nchan,) または (nchan, ny, nx) の bool mask

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mask=lf,
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

意味:
- `linefree_mode` を省略すると `'infer'` として解釈される
- `linefree_mask` があるので、内部的には **manual_mask** 扱いになる
- すなわち base mask は `linefree_mask`
- ただし、`linefree_velocity_windows_kms` を同時に与えた場合は **OR で追加**、`exclude_v_windows` は最後に除外される

**本当に `linefree_mask` だけをそのまま使いたいなら、`linefree_velocity_windows_kms` を一緒に渡さないこと。**

---

## 1.6 `infer` に任せたい

```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_cfg=lf_cfg,
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

意味:
- `linefree_mode` を省略すると `'infer'`
- `infer` は次の優先順位で解決する

1. `linefree_mask` があれば manual_mask
2. prior があれば prior
3. `linefree_cfg` があれば auto
4. どれも無ければ error

**重要**:
- `linefree_velocity_windows_kms` の有無だけでは `manual` にはならない
- `linefree_cfg=None` で prior も無ければ、`infer` は error になる

---

# 2. まず定義を固定する

この章では、曖昧さを避けるために、用語を固定する。

## 2.1 data の軸
内部では cube を

\[
(n_{chan}, n_y, n_x)
\]

として扱う。

- 第 0 軸: スペクトル軸
- 第 1 軸: y
- 第 2 軸: x

## 2.2 line-free
baseline fit に使ってよい channel / voxel の集合。

- 1D 経路では `(nchan,)` の bool mask
- 3D 経路では `(nchan, ny, nx)` の bool mask

## 2.3 include 窓
`linefree_velocity_windows_kms` によって与える速度窓。

ただしこの引数の意味は `linefree_mode` に依存する。

- `manual` では **厳密な manual line-free**
- それ以外では **最後に OR で追加する include 窓**

## 2.4 exclude 窓
`exclude_v_windows` により与える速度窓。

これは常に

- **fit に使わない**
- **最後に除外する**

という意味である。

## 2.5 prior
入力 bundle / FITS に既に存在する `LINEFREE` / `LINEFREE3D` 系の mask。

## 2.6 auto
`LineFreeConfig` に基づいて自動推定した line-free mask。

## 2.7 base mask と final mask
この文書では次の2段階を分ける。

- **base mask**: `manual` / `linefree_mask` / `prior` / `auto` / `prior ∪ auto` で最初に決まる mask
- **final mask**: その後に include / exclude を適用して最終的に fit に使う mask

---

# 3. 公開 API の line-free 関係引数

`subtract_baseline_from_bundle(...)` および `subtract_baseline_from_fits(...)` で、line-free の意味を持つ主要引数は次の通りである。

```python
linefree_cfg: Optional[LineFreeConfig] = None
linefree_mask: Optional[np.ndarray] = None
linefree_detection_data: Optional[np.ndarray] = None
linefree_velocity_windows_kms: Optional[Sequence[Union[str, Tuple[float, float]]]] = None
exclude_v_windows: Optional[Sequence[Union[str, Tuple[float, float]]]] = None
linefree_mode: Optional[str] = 'infer'
load_prior_from_input: bool = True
```


## 3.0 まず地図: どのパラメータがどこに属しているか

この節は、今回の混乱を避けるために追加した。

`subtract_baseline_from_bundle(...)` / `subtract_baseline_from_fits(...)` には多数の引数があるが、
それらは全部が同じ層の設定ではない。大きく分けると次の4層で読むと分かりやすい。

1. **本体 API の直引数**
   - 今回の実行で何を入力し、どの mode を使い、どの補助窓を適用するかを決める
2. **`BaselineConfig`**
   - 最終 baseline fit のモデルと数値解法を決める
3. **`LineFreeConfig`**
   - auto line-free detection の方法を決める
4. **`RippleConfig`**
   - ripple 周波数の auto 探索方法を決める

要するに、**本体 API は「何をするか」を束ねる窓口**であり、
`BaselineConfig` / `LineFreeConfig` / `RippleConfig` はその中で使う**下位設定オブジェクト**である。

### 3.0.1 本体 API の直引数

以下は `subtract_baseline_from_bundle(...)` / `subtract_baseline_from_fits(...)` の**本体側の引数**である。
ここで「本体側」というのは、dataclass の中ではなく、関数に直接書く引数という意味である。

#### A. 入出力・対象データ

- `bundle` / `input_fits`, `output_fits`
  - 何に対して baseline subtraction を行うか
- `cube_ext`
  - legacy FITS path でどの extension を cube として読むか

#### B. line-free base source を決める引数

- `linefree_mode`
  - `manual / prior / auto / or / infer` を選ぶ本体スイッチ
- `linefree_mask`
  - bool mask を直接与える manual source
- `linefree_cfg`
  - auto detection を使う場合の設定オブジェクト
- `load_prior_from_input`
  - prior `LINEFREE` / `LINEFREE3D` を入力から読むかどうか
- `linefree_detection_data`
  - line-free 判定だけ別データで行いたいときの detection-only cube

#### C. line-free base source の後で適用する修飾引数

- `linefree_velocity_windows_kms`
  - `manual` では exact manual windows
  - それ以外では final mask に OR 追加する include 窓
- `exclude_v_windows`
  - 最後に除外する窓

#### D. ripple 周波数の供給元を決める引数

- `ripple_freqs`
  - ripple 周波数を明示的に与える
- `ripple_mode`
  - `'auto'` / `'prior'` など、周波数をどこから持ってくるかを決める
- `ripple_cfg`
  - auto 探索時の探索設定

#### E. ripple をどの段階で使うかを決める引数

- `safe_velocity_windows_kms`
  - ripple 周波数推定や detection 用前処理に使う safe 窓
- `ripple_apply_stage`
  - detection 側 / final 側 / 両方のどこで ripple を使うか

#### F. 最終 fit そのものの設定

- `baseline_cfg`
  - baseline model 本体の設定オブジェクト
- `reproducible_mode`
  - 必要に応じて `baseline_cfg.reproducible_mode` を上書きする本体側フラグ

#### G. 出力・QC・bundle 更新

- `add_qc_ext`, `add_qc_hdus`
  - QC HDU を出力するか
- `update_mosaic_products`
  - baseline 後に mosaic 系の派生 products を更新するか
- `gain_min`
  - mosaic 系更新時の edge/gain 関係しきい値
- `overwrite`, `write_diagnostics`, `diagnostics_prefix`, `write_profile`, `profile_prefix`
  - FITS path / diagnostics path の出力制御

### 3.0.2 どの config に属するか

ここでは、「どの引数がどの config の中身なのか」を一覧で固定する。

#### `BaselineConfig` に属するもの

`BaselineConfig` は、**最終 baseline fit の model と solver**の設定である。

主な所属パラメータ:

- model 形状
  - `poly_order`
  - `ripple`
- robust fit
  - `robust`
  - `robust_mode`
  - `robust_iters`
  - `robust_early_stop`
  - `robust_coef_rtol`
  - `robust_coef_atol`
  - `robust_selective_sigma`
  - `robust_selective_frac`
  - `robust_selective_max_pixels`
  - `robust_batch_pixels`
- 数値解法・再現性
  - `rcond`
  - `chunk_pix`
  - `reproducible_mode`
  - `compute_dtype`
  - `normalize_x`
  - `strict_failures`
  - `fallback_to_pixelwise`
- voxel 3D fit 系
  - `voxel_solver`
  - `voxel_qr_batch_pix`
  - `voxel_grouping`
  - `voxel_group_min_size`

**重要**:
`ripple_mode` や `ripple_cfg` は `BaselineConfig` には属さない。
それらは「ripple model をどう解くか」ではなく、**ripple 周波数をどこから持ってくるか**の本体側設定である。

#### `LineFreeConfig` に属するもの

`LineFreeConfig` は、**auto line-free detection** の設定である。

典型的な所属パラメータ:

- 1D global auto 系
  - `smooth_width`
  - `sigma`
  - `iters`
  - `pad_chan`
  - `min_linefree_frac`
- 3D voxel auto 系
  - `use_3d`
  - `spatial_sigma_pix`
  - `velocity_sigma_chan`
  - `seed_sigma`
  - `grow_sigma`
  - `min_linefree_frac_global`
  - `seed_min_size`
  - `grow_max_iters`
  - `grow_connectivity`
  - `velocity_windows_kms`
  - `erode_chan`
  - `dilate_chan`
  - `min_component_voxels`
  - `fallback_to_1d`

**重要**:
`linefree_velocity_windows_kms` は `LineFreeConfig` には属さない。
これは auto detection の内部パラメータではなく、**本体 API 側の include/exact window 指定**である。

#### `RippleConfig` に属するもの

`RippleConfig` は、**ripple 周波数の auto 探索**の設定である。

典型的な所属パラメータ:

- 探索する本数
  - `nfreq`
- 探索範囲
  - `freq_min_cyc_per_ch`
  - `freq_max_cyc_per_ch`
  - `period_range_chan`
  - `period_range_kms`
  - `period_range_hz`
- 探索時の前処理
  - `window`
  - `detrend_poly_order`
  - `peak_prominence`
  - `peak_height`
  - `peak_distance`
  - `peak_width`

**重要**:
`safe_velocity_windows_kms` や `ripple_apply_stage` は `RippleConfig` には属さない。
それらは「探索設定」ではなく、**本体 API が ripple をどの段階で使うか**を決める引数である。

### 3.0.3 実務上、まず触るべき「本体パラメータ」はどれか

全部を毎回意識する必要はない。多くの場合、最初に考えるべき本体パラメータは次の順である。

1. `linefree_mode`
   - manual / prior / auto / or / infer のどれにしたいか
2. `linefree_velocity_windows_kms` と `exclude_v_windows`
   - exact manual windows なのか、include/exclude 修飾なのか
3. `baseline_cfg`
   - `poly_order`, `ripple`, 必要なら `robust` など
4. `linefree_cfg`
   - auto を使うときだけ調整する
5. `ripple_mode`, `ripple_cfg`, `safe_velocity_windows_kms`, `ripple_apply_stage`
   - ripple を使うときだけ調整する

つまり、通常の読解順は


a) **本体の mode と窓** を決める  
b) **baseline model** を決める  
c) 必要なら **auto line-free** を調整する  
d) 必要なら **ripple 周波数探索** を調整する

である。

### 3.0.4 一番よくある誤読

#### 誤読1: `linefree_velocity_windows_kms` は `LineFreeConfig` の一部である
違う。これは本体 API 側の引数であり、`manual` では exact、他 mode では include である。

#### 誤読2: `ripple_mode` は `BaselineConfig` の一部である
違う。`BaselineConfig.ripple` は model の on/off であり、`ripple_mode` は周波数 source の選択である。

#### 誤読3: `safe_velocity_windows_kms` は signal 窓である
違う。これは ripple 周波数推定や detection 前処理で安全に使う窓であり、`linefree_velocity_windows_kms` と役割が異なる。

### 3.0.5 まずこの表だけ見ればよい

| 何を決めるか | 主に触る場所 | 代表パラメータ |
|---|---|---|
| line-free の base source | 本体 API | `linefree_mode`, `linefree_mask`, `linefree_cfg`, `load_prior_from_input` |
| line-free の include / exclude | 本体 API | `linefree_velocity_windows_kms`, `exclude_v_windows` |
| 最終 baseline model | `BaselineConfig` | `poly_order`, `ripple`, `robust` |
| auto line-free の方法 | `LineFreeConfig` | `smooth_width`, `sigma`, `pad_chan`, `use_3d`, `seed_sigma` |
| ripple 周波数の source | 本体 API | `ripple_freqs`, `ripple_mode` |
| ripple auto 探索の細部 | `RippleConfig` | `nfreq`, `period_range_*`, `peak_*` |
| ripple をどの段階で使うか | 本体 API | `safe_velocity_windows_kms`, `ripple_apply_stage` |
| 出力/QC | 本体 API | `add_qc_ext`, `update_mosaic_products`, `gain_min` |


## 3.1 `linefree_cfg`
自動探索設定。

- 型: `LineFreeConfig | None`
- 既定値: `None`
- `None` のとき、auto detection はその引数だけでは起動しない

## 3.2 `linefree_mask`
channel / voxel ごとの bool mask を直接与える。

- 1D: `(nchan,)`
- 3D: `(nchan, ny, nx)`

## 3.3 `linefree_detection_data`
mask 同定専用に使う detection-only cube。

- 最終 baseline fit の対象データとは別に、line-free 判定だけ別データで行いたいときに使う
- strong ripple や補助 preprocessing を user-side で行った場合に有用

## 3.4 `linefree_velocity_windows_kms`
速度窓指定。

形式例:

```python
['-30:0', '20:55']
[( -30.0, 0.0 ), (20.0, 55.0)]
```

ただし意味は `linefree_mode` に依存する。

## 3.5 `exclude_v_windows`
明示的に baseline fit から外す速度窓。

## 3.6 `linefree_mode`
有効値:

- `'infer'`
- `'manual'`
- `'prior'`
- `'auto'`
- `'or'`

既定値は `'infer'`。

## 3.7 `load_prior_from_input`
`True` なら、入力から prior `LINEFREE` 系を読む。

- `False` なら prior を見ない
- `infer` の解決にも影響する

---

# 4. `linefree_mode` の意味

## 4.1 `'manual'`
**今回追加された新モード。**

このモードでは、`linefree_velocity_windows_kms` をそのまま厳密な line-free として使う。

### base mask
\[
L_{base} = W_{manual}
\]

ここで
- \(W_{manual}\): `linefree_velocity_windows_kms` から作られる 1D mask

### その後
- `linefree_velocity_windows_kms` の再 OR 追加はしない
- `exclude_v_windows` は最後に除外

したがって、最終的には

\[
L_{final} = W_{manual} \cap \neg E
\]

である。

### 注意
- `linefree_mode='manual'` には `linefree_velocity_windows_kms` が必須
- `linefree_mask` との併用は不可

---

## 4.2 `'prior'`
入力から読んだ prior line-free mask を base mask に使う。

### base mask
\[
L_{base} = L_{prior}
\]

### その後
- `linefree_velocity_windows_kms` があれば OR 追加
- `exclude_v_windows` があれば最後に除外

したがって

\[
L_{final} = (L_{prior} \cup W_{include}) \cap \neg E
\]

である。

### 注意
- prior が無ければ error
- `load_prior_from_input=False` なら prior は読まれない

---

## 4.3 `'auto'`
`LineFreeConfig` に基づいて自動探索した mask を base mask に使う。

### base mask
\[
L_{base} = L_{auto}
\]

### その後
- `linefree_velocity_windows_kms` があれば OR 追加
- `exclude_v_windows` があれば最後に除外

したがって

\[
L_{final} = (L_{auto} \cup W_{include}) \cap \neg E
\]

である。

### 注意
- `linefree_cfg` が無ければ error
- 1D global auto か 3D voxel auto かは `LineFreeConfig.mask_kind` と `auto_method` に依存する

---

## 4.4 `'or'`
`prior` と `auto` の和集合を base mask に使う。

### base mask
\[
L_{base} = L_{prior} \cup L_{auto}
\]

ただし、片方しか存在しない場合はその片方だけでもよい。

### その後
- `linefree_velocity_windows_kms` があれば OR 追加
- `exclude_v_windows` があれば最後に除外

したがって

\[
L_{final} = (L_{prior} \cup L_{auto} \cup W_{include}) \cap \neg E
\]

である。

### 注意
- prior も auto も両方無ければ error
- `linefree_cfg is None` なら auto は作れない
- `load_prior_from_input=False` なら prior は読まれない

---

## 4.5 `'infer'`
入力状況から自動的に base source を決める既定モード。

解決規則は次の通り。

1. `linefree_mask` があれば **manual_mask**
2. そうでなく prior があれば **prior**
3. そうでなく `linefree_cfg` があれば **auto**
4. どれも無ければ error

### 重要
- `linefree_velocity_windows_kms` の有無だけでは `manual` にはならない
- `infer` は「勝手に auto」ではなく、「使える情報から base source を決める」ためのモードである

### internal manual_mask
`linefree_mask` が与えられた場合、内部的には public な `manual` とは別に `manual_mask` として処理される。

このとき

\[
L_{base} = L_{mask}
\]

であり、その後

\[
L_{final} = (L_{mask} \cup W_{include}) \cap \neg E
\]

となる。

つまり、`linefree_mask` を与えた場合でも `linefree_velocity_windows_kms` を同時に与えると OR 追加される。

---

# 5. 真理値表ではなく、運用表で見る

## 5.1 base source の決まり方

| public 指定 | base source | auto は走るか | prior は読むか |
|---|---|---:|---:|
| `linefree_mode='manual'` | `linefree_velocity_windows_kms` | no | no |
| `linefree_mode='prior'` | prior | no | yes |
| `linefree_mode='auto'` | auto | yes | optional |
| `linefree_mode='or'` | prior ∪ auto | yes if cfg given | yes |
| `linefree_mode='infer'` + `linefree_mask` | `linefree_mask` | no | irrelevant |
| `linefree_mode='infer'` + prior only | prior | no | yes |
| `linefree_mode='infer'` + cfg only | auto | yes | no prior available |
| `linefree_mode='infer'` + prior + cfg | prior | no | yes |

最後の行は重要である。
`infer` で prior と cfg の両方がある場合、既定では **prior** が優先される。`prior ∪ auto` を使いたいなら `linefree_mode='or'` を明示する。

## 5.2 include / exclude の適用

| base source | `linefree_velocity_windows_kms` | `exclude_v_windows` |
|---|---|---|
| manual | exact selection として使用 | 最後に除外 |
| prior | OR 追加 | 最後に除外 |
| auto | OR 追加 | 最後に除外 |
| or | OR 追加 | 最後に除外 |
| manual_mask | OR 追加 | 最後に除外 |

---

# 6. `linefree_velocity_windows_kms` の意味をもう一度固定する

この引数は今回もっとも誤解されやすかった。

## 6.1 manual のとき

```python
linefree_mode='manual'
linefree_velocity_windows_kms=['-30:0', '20:55']
```

このときの意味は:

- `-30:0` と `20:55` **だけ**を fit に使う
- auto detection は走らない
- prior も使わない

## 6.2 manual 以外のとき

```python
linefree_mode='auto'
linefree_velocity_windows_kms=['-30:0', '20:55']
```

このときの意味は:

- base mask は auto で作る
- その後 `-30:0`, `20:55` を **include 窓として OR 追加**する

したがって、manual 以外でこの引数を与えても、「この窓だけを使う」意味にはならない。

---

# 7. `exclude_v_windows` の意味

`exclude_v_windows` は、どのモードでも **最後に除外**する。

例:

```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:20'],
    exclude_v_windows=['0:15'],
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

この場合、最終 line-free は

\[
[-30, 20] \setminus [0, 15]
\]

となる。

したがって、実質的には

- `[-30,0)`
- `(15,20]`

だけが fit に使われる。

---

# 8. 速度軸の意味と channel 範囲

速度軸は線形 header から

\[
v(i)=CRVAL3 + (i-CRPIX3) \times CDELT3
\]

で作る。

ここで
- `i` は FITS 1-origin channel index
- `CRVAL3`, `CRPIX3`, `CDELT3` は header から読む

例えば

- `CTYPE3 = VRAD`
- `CUNIT3 = km/s`
- `CRVAL3 = -30.0`
- `CDELT3 = 0.2`
- `CRPIX3 = 1.0`

なら、

- channel 1 は `-30.0 km/s`
- channel 151 は `0.0 km/s`
- channel 241 は `18.0 km/s`
- channel 426 は `55.0 km/s`

である。

したがって

```python
linefree_mode='manual'
linefree_velocity_windows_kms=['-30:0', '18:55']
```

は理想的には

- ch 1–151
- ch 241–426

を line-free にする。

このような channel 範囲の確認は、manual 指定が意図通り効いているかを見る上で有用である。

---

# 9. `LineFreeConfig` は何のためにあるか

`LineFreeConfig` は auto detection の設定である。

したがって:

- `linefree_mode='manual'` のときは base mask 決定には使わない
- `linefree_mode='auto'` のときは base mask を作るために使う
- `linefree_mode='or'` のときは auto 側の base mask を作るために使う
- `linefree_mode='infer'` のときは、prior が無く `linefree_mask` も無い場合に auto 側へ解決されるために使う

## 9.1 1D global auto
従来型。最終 mask は 1D `(nchan,)` である。

代表 spectrum を作り、

1. ゆっくりした baseline を推定
2. residual に対する robust sigma clip
3. `pad_chan` による拡張

などで global line-free mask を作る。

## 9.2 3D voxel auto
より新しい経路。最終 mask は `(nchan, ny, nx)` である。

- detection cube を作る
- `x,y` smoothing と `v` smoothing を使って Stage A local seed を安定化
- Stage B provisional fit で更新
- Stage C で final fit

最終 fit 自体は元の unsmoothed cube に対して行う。

---

# 10. `infer` をどう使うべきか

`infer` は convenience 用の既定モードだが、曖昧に使わない方がよい。

## 10.1 安全な使い方

### 既に prior が入っている bundle を再利用したい
```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='prior',
    baseline_cfg=..., 
)
```

### auto detection を使うと明示したい
```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_cfg=lf_cfg,
    linefree_mode='auto',
    baseline_cfg=..., 
)
```

### exact manual windows を使いたい
```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=..., 
)
```

### bool mask をそのまま使いたい
```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mask=lf,
    baseline_cfg=..., 
)
```

## 10.2 `infer` を省略形として使う場合
省略は次のときだけ分かりやすい。

- `linefree_mask` を与えるとき
- prior が確実に入っていると分かっているとき
- `linefree_cfg` を与えて auto にしたいとき

## 10.3 `infer` で避けるべき書き方

```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=..., 
)
```

これは今回の改訂版では **不適切**である。

理由:
- `linefree_velocity_windows_kms` の有無だけでは base source は決まらない
- prior も `linefree_cfg` も無ければ `infer` は error になる
- exact manual windows にしたいなら `linefree_mode='manual'` を明示すべき

---

# 11. baseline 後に何が更新されるか

baseline 後の bundle / FITS では、少なくとも line-free / signal / QC 関係の ext が更新される。

典型的には:

- `LINEFREE_USED`
- `SIGNAL_MASK_USED`
- `BASE_RMS`
- `BASE_FLG`
- `RIPFREQ_USED`
- `MOMENT0`
- `MOSAIC_*` 系（`update_mosaic_products=True` の場合）

3D voxel path を使った場合は:

- `LINEFREE3D_USED`
- `SIGNAL_MASK3D_USED`

も更新される。

---

# 12. `SIGNAL_MASK_USED` は何を意味するか

これは line-free の補集合全体ではない。

`LINEFREE_USED` の **内部 gap** のみを signal とする。

1D mask なら、概念的には

- 先頭側の外側 gap は signal にしない
- 末尾側の外側 gap も signal にしない
- line-free に挟まれた内部 gap だけを signal にする

という規則である。

したがって、

```python
bundle_sig = m3d_bsl.make_baseline_viewer_bundle(bundle_bl, mode='signal')
```

で見えるものは、**fit に使った channel そのものではない**。

fit に使った channel を見たいときは

```python
bundle_lf = m3d_bsl.make_baseline_viewer_bundle(bundle_bl, mode='linefree')
```

を見る。

3D voxel path なら

- `mode='signal3d'`
- `mode='linefree3d'`

を使う。

---

# 13. viewer bundle の使い方

## 13.1 line-free 側だけ見たい

```python
bundle_lf = m3d_bsl.make_baseline_viewer_bundle(
    bundle_bl,
    mode='linefree',
)
```

## 13.2 signal 側だけ見たい

```python
bundle_sig = m3d_bsl.make_baseline_viewer_bundle(
    bundle_bl,
    mode='signal',
)
```

## 13.3 3D voxel line-free 側だけ見たい

```python
bundle_lf3d = m3d_bsl.make_baseline_viewer_bundle(
    bundle_bl,
    mode='linefree3d',
)
```

## 13.4 3D voxel signal 側だけ見たい

```python
bundle_sig3d = m3d_bsl.make_baseline_viewer_bundle(
    bundle_bl,
    mode='signal3d',
)
```

保存:

```python
import sd_radio_spectral_fits.map_3d as m3d

m3d.write_otf_bundle(bundle_lf,   'linefree_only.fits',   overwrite=True)
m3d.write_otf_bundle(bundle_sig,  'signal_only.fits',     overwrite=True)
m3d.write_otf_bundle(bundle_lf3d, 'linefree3d_only.fits', overwrite=True)
m3d.write_otf_bundle(bundle_sig3d,'signal3d_only.fits',   overwrite=True)
```

---

# 14. Cookbook

## 14.1 exact manual windows + poly only

```python
import sd_radio_spectral_fits.map_3d.baseline_subtraction as m3d_bsl

b_cfg = m3d_bsl.BaselineConfig(
    poly_order=1,
    ripple=False,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=b_cfg,
    add_qc_hdus=True,
    update_mosaic_products=True,
    gain_min=0.5,
)
```

## 14.2 exact manual windows + poly only + exclude

```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:20'],
    exclude_v_windows=['0:15'],
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

## 14.3 exact manual windows + poly+ripple

```python
r_cfg = m3d_bsl.RippleConfig(
    nfreq=2,
    period_range_chan=(20.0, 400.0),
    min_separation=0.002,
    window='hann',
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    ripple_cfg=r_cfg,
    ripple_mode='auto',
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=True),
)
```

## 14.4 auto + include 窓 + exclude 窓

```python
lf_cfg = m3d_bsl.LineFreeConfig(
    smooth_width=51,
    sigma=4.0,
    iters=6,
    pad_chan=3,
    min_linefree_frac=0.35,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_cfg=lf_cfg,
    linefree_mode='auto',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    exclude_v_windows=['5:15'],
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

意味:
- auto でまず推定
- `-30:0` と `20:55` は include で追加
- `5:15` は最後に除外

## 14.5 prior + include

```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='prior',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

意味:
- prior を base mask にする
- 指定窓は OR で追加

## 14.6 prior ∪ auto

```python
lf_cfg = m3d_bsl.LineFreeConfig(
    smooth_width=51,
    sigma=4.0,
    iters=6,
    pad_chan=3,
    min_linefree_frac=0.35,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_cfg=lf_cfg,
    linefree_mode='or',
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

## 14.7 linefree_mask を直接使う

```python
lf = ...

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mask=lf,
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
)
```

## 14.8 3D voxel auto line-free

```python
lf_cfg = m3d_bsl.LineFreeConfig(
    mask_kind='voxel_3d',
    auto_method='local_3d_poly',
    lf3d_detect_spatial_fwhm_arcsec=18.0,
    lf3d_detect_spectral_width_kms=0.6,
    lf3d_seed_spatial_fwhm_arcsec=18.0,
    lf3d_seed_median_width_kms=4.0,
    lf3d_seed_threshold_sigma=4.0,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_cfg=lf_cfg,
    linefree_mode='auto',
    baseline_cfg=m3d_bsl.BaselineConfig(
        poly_order=1,
        ripple=False,
        voxel_solver='qr',
    ),
)
```

## 14.9 FITS を直接処理する

```python
m3d_bsl.subtract_baseline_from_fits(
    'in.fits',
    'out.fits',
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=m3d_bsl.BaselineConfig(poly_order=1, ripple=False),
    add_qc_hdus=True,
    overwrite=True,
)
```

---

# 15. よくある誤解

## 15.1 `linefree_velocity_windows_kms` を渡したのに、勝手に auto された

今回の改訂版では、これは public 仕様として起こらないように整理した。

- exact manual にしたいなら `linefree_mode='manual'`
- auto に include したいなら `linefree_mode='auto'` などを明示
- `infer` は base source を推定するだけで、`linefree_velocity_windows_kms` の有無だけでは manual にしない

## 15.2 `mode='signal'` を見たら line が見えた

正常であることがある。

`mode='signal'` は `SIGNAL_MASK_USED` を表示しており、これは `LINEFREE_USED` の内部 gap だからである。fit に使った側を見たいなら `mode='linefree'` を見る。

## 15.3 `infer` で `linefree_velocity_windows_kms` だけを渡した

不適切である。exact manual のつもりなら `linefree_mode='manual'` を付ける。

## 15.4 `linefree_mask` と `linefree_velocity_windows_kms` を一緒に渡した

可能だが、その場合は

\[
L_{final} = (L_{mask} \cup W_{include}) \cap \neg E
\]

になる。exact mask だけを使いたいなら `linefree_velocity_windows_kms` を渡さない。

---

# 16. 最後の実務上の推奨

1. **exact manual windows にしたいときは必ず `linefree_mode='manual'` を書く。**
2. **auto を使うなら `linefree_cfg` を与え、`linefree_mode='auto'` か `'or'` を明示する。**
3. **prior を使いたいときは `linefree_mode='prior'` を明示する。**
4. **`infer` は convenience 用と割り切る。** 共有コード・長期運用コードでは明示モードの方が安全。
5. **fit に使った側を確認したいときは `LINEFREE_USED` または `mode='linefree'` を見る。**
6. **`mode='signal'` は fit 側の確認用ではない。**
7. **broad line / strong line / strong ripple がある場合は、auto だけに頼らず `manual` や `exclude_v_windows` を併用する。**

---

# 17. 最小確認コード

## 17.1 実際に使われた line-free 区間を速度で表示する

```python
import numpy as np

hdr = bundle_bl.header
nchan = bundle_bl.data.shape[0]
v = hdr['CRVAL3'] + (np.arange(nchan) + 1.0 - hdr['CRPIX3']) * hdr['CDELT3']

lf = np.asarray(bundle_bl.image_ext['LINEFREE_USED'], dtype=bool)

runs = []
in_run = False
s = None
for i, x in enumerate(lf):
    if x and not in_run:
        s = i
        in_run = True
    if in_run and (i == len(lf) - 1 or not lf[i + 1]):
        e = i
        runs.append((s, e))
        in_run = False

for s, e in runs:
    print(f'linefree: ch {s+1:3d}-{e+1:3d}   v = {v[s]:7.3f} .. {v[e]:7.3f} km/s')
```

## 17.2 fit 側 / signal 側 viewer bundle を別々に書く

```python
bundle_lf = m3d_bsl.make_baseline_viewer_bundle(bundle_bl, mode='linefree')
bundle_sig = m3d_bsl.make_baseline_viewer_bundle(bundle_bl, mode='signal')

import sd_radio_spectral_fits.map_3d as m3d
m3d.write_otf_bundle(bundle_lf, 'linefree_only.fits', overwrite=True)
m3d.write_otf_bundle(bundle_sig, 'signal_only.fits', overwrite=True)
```

これで、fit に使った側と signal として積分される側を混同しにくくなる。

---

以上。

---

# 18. ripple 設定の責務分離

この章は、`BaselineConfig.ripple`、`ripple_mode`、`ripple_freqs`、`ripple_cfg`、`safe_velocity_windows_kms`、`ripple_apply_stage` が何を決めているのかを、**責務ごとに分けて**整理するための章である。  
名前が似ているため混乱しやすいが、少なくとも論理上は次の 3 層に分かれている。

1. **最終 baseline model に ripple 項を入れるか**
2. **ripple 周波数をどこから持ってくるか**
3. **その ripple 情報をどの段階で使うか**

この 3 つを混同しないことが重要である。

## 18.1 まず最小定義を固定する

### A. final model
最終的に各 spectrum に対して fit する baseline model の形である。

- 多項式だけか
- 多項式 + ripple(sin/cos) か

を決める。

### B. ripple frequency source
ripple 項を入れるなら、その周波数

\[
f_k \quad [\text{cycles/channel}]
\]

をどこから得るかを決める。

### C. ripple usage stage
求めた ripple 情報を

- line-free 自動判定の前処理に使うのか
- 最終 baseline fit に使うのか
- 両方に使うのか

を決める。

この章では、以後この A/B/C の役割名を使う。

## 18.2 各引数が属している責務

### 18.2.1 `BaselineConfig.ripple`
これは **A. final model** 側の設定である。

- `False` : 最終 baseline model は多項式のみ
- `True` : 最終 baseline model に ripple の sin/cos 項を入れる

したがって、`BaselineConfig(poly_order=1, ripple=True)` は

\[
T(i) \approx \sum_{m=0}^{1} a_m x(i)^m + \sum_k \left[b_k \sin(2\pi f_k i) + c_k \cos(2\pi f_k i)\right]
\]

のような形を意味する。ここで

- \(i\) : channel index
- \(x(i)\) : 多項式用の内部座標
- \(f_k\) : cycles/channel 単位の ripple 周波数

である。

**重要**: `BaselineConfig.ripple` は「周波数をどう推定するか」を決める引数ではない。  
あくまで **最終 fit model に ripple 項を入れるかどうか** が本来の責務である。

### 18.2.2 `ripple_freqs`
これは **B. ripple frequency source** 側の設定である。

- `None` でなければ、ユーザーが明示した周波数列を使う
- 単位は **cycles/channel**

最優先で使われるので、`ripple_freqs` を与えた場合は `ripple_mode` は実質的に補助的になる。

### 18.2.3 `ripple_mode`
これも **B. ripple frequency source** 側である。  
`ripple_freqs` を明示しないときに、周波数をどこから持つかの方針を決める。

現行コードでは主に次の意味になる。

- `'prior'` : 入力 bundle / FITS にある prior ripple 周波数を優先
- `'auto'` : prior があっても使わず、FFT による自動推定を使う

したがって `ripple_mode` は、**ripple 項を fit に入れるかどうか** ではなく、**周波数のソース**を決める設定である。

### 18.2.4 `ripple_cfg`
これも **B. ripple frequency source** 側である。  
FFT 探索時の設定であり、例えば

- 探索する周期範囲
- 何本まで採用するか (`nfreq`)
- 事前に落とす多項式次数

などを決める。

つまり `ripple_cfg` は、`ripple_mode='auto'` または prior/explicit が無い場合の **自動推定器の設定**である。

### 18.2.5 `safe_velocity_windows_kms`
これは A でも B でもなく、**B を安定に実行するための補助マスク**である。

役割は、ripple 周波数推定のときに

- line が強い領域
- wing が広い領域
- 明らかに baseline ではない領域

を避けて、比較的安全な速度窓だけを使うことである。

したがって `safe_velocity_windows_kms` は

- final line-free mask を直接決める引数ではない
- final signal mask を直接決める引数でもない
- **ripple 周波数推定用の安全窓**である

と理解するのが正しい。

### 18.2.6 `ripple_apply_stage`
これは **C. ripple usage stage** 側の設定である。

- `'final'` : 最終 baseline fit にだけ使う
- `'mask_only'` : line-free 自動判定のための detection 側だけで使う
- `'both'` : detection 側と final fit 側の両方で使う

ここでいう detection 側とは、auto line-free を安定化するために、必要なら一度 ripple を見込んだ前処理をした data を作ってから mask 推定に使う流れを指す。

## 18.3 現行コードでの実際の流れ

`subtract_baseline_from_bundle()` の ripple 関係の流れは、概略として次の順である。

### 18.3.1 detection 用 ripple 前処理
次の条件が揃うと、line-free 自動判定用に detection cube を作る前段で ripple を使う。

- `safe_velocity_windows_kms` が与えられている
- `ripple_apply_stage in {'mask_only', 'both'}`
- `BaselineConfig.ripple == True`
- `linefree_detection_data` を直接与えていない

このとき safe window 上で ripple 周波数を

1. `ripple_freqs`
2. `ripple_mode='prior'` の prior
3. それ以外なら FFT 自動推定

の順で決め、その周波数を使って detection 用 data を作る。

### 18.3.2 line-free mask 決定
その後、`linefree_mode` と `linefree_cfg` に基づいて base mask を決め、さらに

- `linefree_velocity_windows_kms` の include
- `exclude_v_windows` の exclude

を適用して最終 `LINEFREE_USED` を決める。

ここで大事なのは、**`safe_velocity_windows_kms` は line-free mask 本体を決める引数ではない**ことである。

### 18.3.3 final ripple fit
最後に、`BaselineConfig.ripple == True` かつ `ripple_apply_stage in {'final', 'both'}` のとき、最終 baseline fit 用の ripple 周波数を決める。

このときの周波数決定優先順も同じである。

1. `ripple_freqs`
2. `ripple_mode='prior'` の prior ripple 周波数
3. それ以外なら FFT 自動推定

そして FFT 自動推定のときは、通常は `LINEFREE_USED` から作った 1D mask を使うが、`safe_velocity_windows_kms` が与えられていればそちらを優先する。

## 18.4 重要: `BaselineConfig.ripple` は見た目より強いスイッチである

ここは特に誤解しやすい。

現行コードでは、`BaselineConfig.ripple=False` にすると

- 最終 baseline model に ripple 項が入らない
- だけでなく
- `ripple_apply_stage='mask_only'` や `'both'` の detection 側 ripple 補助も走らない

という形になっている。

つまり実装上は、`BaselineConfig.ripple` が

- final model on/off
- detection 側 ripple 前処理の許可/不許可

の両方を兼ねている。

これは名前だけを見ると分かりにくい。  
論理的には `BaselineConfig.ripple` は final model 側の設定と読むのが自然だが、**現行実装ではそれより少し強い意味を持っている**ことに注意する。

## 18.5 典型的な混乱と正しい読み方

### 18.5.1 `BaselineConfig.ripple=True` にしたのに `ripple_mode` も必要なのか
必要である。  
`BaselineConfig.ripple=True` は「ripple 項を使う」という宣言であり、`ripple_mode` は「その周波数をどこから得るか」を決める。

- `ripple=True` だけでは、どの周波数を使うかは未定
- `ripple_mode` / `ripple_freqs` / `ripple_cfg` が周波数決定を補う

と読む。

### 18.5.2 `ripple_mode='prior'` にしたのに `ripple=False` ならどうなるか
この場合、最終 baseline fit には ripple 項が入らない。  
したがって prior ripple 周波数が存在しても、final fit 側では使われない。

また現行コードでは、detection 側 ripple 補助も `BaselineConfig.ripple` でガードされるため、`ripple_apply_stage='mask_only'` でも `ripple=False` なら何も起きない。

### 18.5.3 `safe_velocity_windows_kms` は signal 窓なのか
違う。  
これは **ripple 周波数推定用の安全窓**である。  
line-free 本体や signal 本体を決める窓ではない。

### 18.5.4 `linefree_velocity_windows_kms` と `safe_velocity_windows_kms` は何が違うか
役割が違う。

- `linefree_velocity_windows_kms` : line-free mask に対する include 修飾子（manual では exact）
- `safe_velocity_windows_kms` : ripple 周波数推定を安全にするための窓

両者は速度窓という点では似ているが、**作用先が違う**。

## 18.6 実務上の推奨

### 18.6.1 ripple を使わない場合
```python
b_cfg = m3d_bsl.BaselineConfig(
    poly_order=1,
    ripple=False,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=b_cfg,
)
```

この場合、ripple 関係の設定はほぼ無関係である。

### 18.6.2 final fit だけで ripple を使いたい場合
```python
b_cfg = m3d_bsl.BaselineConfig(
    poly_order=1,
    ripple=True,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:0', '20:55'],
    baseline_cfg=b_cfg,
    ripple_mode='auto',
    ripple_apply_stage='final',
)
```

意味は

- line-free は manual で固定
- ripple 周波数は auto で決める
- ripple 項は final fit にだけ使う

である。

### 18.6.3 auto line-free を安定化するため detection 側にも ripple を使いたい場合
```python
lf_cfg = m3d_bsl.LineFreeConfig(
    mode='global_1d',
)

r_cfg = m3d_bsl.RippleConfig(
    nfreq=1,
)

b_cfg = m3d_bsl.BaselineConfig(
    poly_order=1,
    ripple=True,
)

bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='auto',
    linefree_cfg=lf_cfg,
    baseline_cfg=b_cfg,
    ripple_cfg=r_cfg,
    ripple_mode='auto',
    safe_velocity_windows_kms=['-30:0', '20:55'],
    ripple_apply_stage='both',
)
```

意味は

- line-free base mask は auto
- ripple 周波数は safe window 上で auto 推定
- その ripple 情報を detection 側にも final fit 側にも使う

である。

### 18.6.4 detection 側は user-side で別途処理済み data を使いたい場合
```python
bundle_bl = m3d_bsl.subtract_baseline_from_bundle(
    bundle_in,
    linefree_mode='auto',
    linefree_cfg=lf_cfg,
    linefree_detection_data=data_det,
    baseline_cfg=b_cfg,
    ripple_mode='auto',
    ripple_apply_stage='final',
)
```

この場合、auto line-free のための detection data は `linefree_detection_data` が優先される。  
したがって `mask_only` 側の内部 ripple 前処理には依存しない。

## 18.7 最後に一文でまとめる

- `BaselineConfig.ripple` は **最終 fit model に ripple 項を入れるか**
- `ripple_freqs` / `ripple_mode` / `ripple_cfg` は **その周波数をどこから得るか**
- `ripple_apply_stage` は **その ripple 情報を detection / final のどこで使うか**
- `safe_velocity_windows_kms` は **ripple 周波数推定の安全窓**

である。

この 4 行を混同しなければ、ripple 関係の設定は整理して理解できる。



---

# 19. voxel_3d の detection smoothing と `lf3d_min_run_kms`

この章では、`mask_kind='voxel_3d'` のときに特に混乱しやすい次の 3 つの引数を、**単位・内部表現・適用段階・境界条件まで含めて**整理する。

- `lf3d_detect_spatial_fwhm_arcsec`
- `lf3d_detect_spectral_width_chan`
- `lf3d_min_run_kms`

ここで最初に重要なことを明示しておく。

1. これら 3 つはすべて **`LineFreeConfig` に属する**。
2. これら 3 つはすべて **voxel_3d の auto line-free 推定のための設定**であり、`BaselineConfig` の設定ではない。
3. これら 3 つは **最終 baseline fit の model そのもの**を直接変える引数ではない。
4. ただし、auto line-free 推定結果を通じて、最終的にどの channel / voxel を baseline 側に残すかには影響する。

以下、変数を明示しておく。

- 入力 cube: `D(v, y, x)`
- detection 用 cube: `D_det(v, y, x)`
- provisional な 3D signal mask: `S(v, y, x)`
- 速度軸: `v_i` [km/s]
- channel 幅の代表値: 

\[
\Delta v = \mathrm{median}\left(|v_{i+1}-v_i|\right)
\]

- 空間 pixel size: `cell_arcsec` [arcsec/pix]

この章では、**実装上そうなっていること**と、**論理的にどう解釈すべきか**を分けて記す。

## 19.1 まず 3 つの引数の所属先と責務を一行でまとめる

### 19.1.1 `lf3d_detect_spatial_fwhm_arcsec`
- 所属先: `LineFreeConfig`
- 単位: arcsec
- 役割: detection cube を作るときに、**xy 平面に掛ける Gaussian smoothing kernel の FWHM** を指定する
- 注意: **target beam の FWHM ではなく、追加で掛ける smoothing kernel の FWHM** である

### 19.1.2 `lf3d_detect_spectral_width_chan`
- 所属先: `LineFreeConfig`
- 単位: channel
- 役割: detection cube を作るときに、**v 方向に掛ける boxcar smoothing の幅** を指定する
- 注意: **binning ではない**。channel 数は減らない

### 19.1.3 `lf3d_min_run_kms`
- 所属先: `LineFreeConfig`
- 単位: km/s
- 役割: provisional な 3D signal mask に対して、各 LOS ごとの **短すぎる連続速度 run を baseline 側へ戻す post-prune** を行う
- 注意: detection smoothing そのものではない。**smoothing の後、mask ができた後の cleanup** である

## 19.2 これら 3 つはどの段階で使われるか

`voxel_3d` の auto line-free は概念的に次の順で進む。

1. 元の cube `D(v,y,x)` から **detection 用 cube** `D_det(v,y,x)` を作る
2. `D_det` から seed / provisional な 3D signal mask `S` を作る
3. `S` を使って provisional baseline を作る
4. 反復しながら signal / baseline の分離を更新する
5. 最後に unsmoothed な元の cube に対して final baseline fit を行う

このとき、今回の 3 引数の適用位置は次の通りである。

- `lf3d_detect_spatial_fwhm_arcsec`
  - **1. detection cube 作成段階**で使う
- `lf3d_detect_spectral_width_chan`
  - **1. detection cube 作成段階**で使う
- `lf3d_min_run_kms`
  - **2. seed / provisional mask ができた後**と、**反復中に更新された line mask の直後**で使う

したがって、

- `lf3d_detect_*` は **data smoothing**
- `lf3d_min_run_kms` は **mask cleanup**

という役割分担である。

## 19.3 `lf3d_detect_spatial_fwhm_arcsec` の意味

### 19.3.1 論理的な意味

この引数は、detection cube を作るときに xy 平面へ掛ける Gaussian smoothing kernel の FWHM を arcsec 単位で与えるための引数である。

つまり、

\[
D_{det}(v,y,x)
=
G_{xy}(\mathrm{FWHM}=F_{det}) * D(v,y,x)
\]

の `F_det` に相当する。

ここで重要なのは、これは

- 「最終 beam を 350 arcsec にする」
- 「元の beam が何 arcsec なので追加 kernel は何 arcsec」

を自動で逆算する引数ではない、ということである。

**実装上は単純に “その FWHM の Gaussian kernel を掛ける”**。

### 19.3.2 実装上の処理

内部 solver は pixel 単位で動くので、public な arcsec 指定は内部で pixel 単位の `sigma` に変換される。

まず

\[
\mathrm{FWHM}_{pix} = \frac{\mathrm{lf3d\_detect\_spatial\_fwhm\_arcsec}}{\mathrm{cell\_arcsec}}
\]

を作り、次に Gaussian の関係式

\[
\mathrm{FWHM} = \sqrt{8\ln 2}\,\sigma
\]

から

\[
\sigma_{pix} = \frac{\mathrm{FWHM}_{pix}}{\sqrt{8\ln 2}}
\]

へ変換して、xy 方向の `gaussian_filter` に渡す。

したがって、たとえば

- `cell_arcsec = 36 arcsec/pix`
- `lf3d_detect_spatial_fwhm_arcsec = 350 arcsec`

なら

\[
\mathrm{FWHM}_{pix} \approx 350/36 \approx 9.72\ \mathrm{pix}
\]

\[
\sigma_{pix} \approx 9.72 / 2.3548 \approx 4.13\ \mathrm{pix}
\]

となる。

### 19.3.3 何をしていないか

この引数は **target beam 化** をしていない。

つまり、元の実効 beam が既に `B_in` だとしても、

\[
B_{out} = \sqrt{B_{in}^2 + F_{det}^2}
\]

に近い方向へ進む。ここで `F_det` はこの引数で与えた Gaussian kernel の FWHM である。

したがって、

- 「検出用に 350 arcsec の Gaussian でさらに滑らかにしたい」

なら意図通りであるが、

- 「最終的な detection 分解能を 350 arcsec にしたい」

という意味ではない。

### 19.3.4 `lf3d_detect_spatial_fwhm_pix` との関係

内部には pixel 単位の alias もある。

- `lf3d_detect_spatial_fwhm_pix`
- `lf3d_detect_spatial_fwhm_arcsec`

の両方を持つが、public には arcsec 版を使う方が自然である。

実装上は

1. まず `*_pix`
2. 次に `*_arcsec`

の順に処理されるため、**両方を同時に与えた場合は `*_arcsec` が後勝ち**である。

### 19.3.5 無効化条件

- `None` のときは alias としては何もしない
- 内部値 `lf3d_detect_spatial_sigma=0.0` なら smoothing は実質無効

## 19.4 `lf3d_detect_spectral_width_chan` の意味

### 19.4.1 論理的な意味

この引数は、detection cube を作るときに v 方向へ掛ける **boxcar smoothing 幅** を channel 数で与える引数である。

つまり、各 LOS に対して概念的には

\[
D_{det}(v_i,y,x)
=
\frac{1}{W}
\sum_{j \in \mathcal{N}_W(i)} D(v_j,y,x)
\]

の `W` を決める引数である。

ここで `\mathcal{N}_W(i)` は中心 `i` のまわりの `W` channel の近傍である。

### 19.4.2 これは binning ではない

この引数は **binning** ではない。

- channel 数は減らない
- downsample しない
- 5ch 指定でも出力 channel 数は元のまま

やっていることは **5ch moving-average smoothing** である。

したがって、コメントは

```python
lf3d_detect_spectral_width_chan=5  # detection 用に v 方向を 5ch boxcar smoothing
```

のように書くのが正しい。

`# 5ch で binning` と書くのは実装の意味と一致しない。

### 19.4.3 実装上の処理

内部では `uniform_filter1d(..., size=W, axis=0)` 相当の処理をしている。

したがって、幅 `W=5` なら 5 channel の boxcar smoothing である。

### 19.4.4 偶数を入れた場合どうなるか

中心対称の smoothing にしたいため、内部では **奇数幅へ正規化** される。

概念的には

\[
W = \mathrm{round}(W_{in})
\]

\[
W \leftarrow \max(1, W)
\]

\[
W\ \text{が偶数なら}\ W \leftarrow W+1
\]

である。

したがって

- `4` → 実際には `5`
- `5` → 実際には `5`
- `6` → 実際には `7`

になる。

**偶数入力は許されるが、そのまま偶数幅では動かない。**

### 19.4.5 float を入れた場合どうなるか

public alias は `float` を受け取れるが、内部で最終的には channel 幅へ変換され、奇数整数へ丸められる。

したがって

- `5.0` は問題ない
- `4.6` は丸め後に奇数化される

が、曖昧さを避けるためには **最初から奇数整数相当の値を入れる** のが望ましい。

### 19.4.6 `lf3d_detect_spectral_width_kms` との関係

同じ検出用 spectral smoothing には km/s 版 alias もある。

- `lf3d_detect_spectral_width_chan`
- `lf3d_detect_spectral_width_kms`

`*_kms` を使うと、まず

\[
W_{raw} = \frac{\mathrm{lf3d\_detect\_spectral\_width\_kms}}{\Delta v}
\]

を作り、その後に上と同じ奇数化を行う。

実装上は

1. まず `*_chan`
2. 次に `*_kms`

の順で処理されるため、**両方を同時に与えた場合は `*_kms` が後勝ち**である。

### 19.4.7 無効化条件

- `None` なら alias としては何もしない
- 内部値 `lf3d_detect_spectral_width=1` は実質 smoothing なしに近い
- `<=0` のような不正値は奇数化前に 1 以上へ正規化される設計である

## 19.5 `lf3d_min_run_kms` の意味

### 19.5.1 論理的な意味

この引数は、provisional な 3D signal mask `S(v,y,x)` ができた後に、**各 LOS ごとの連続速度 run が短すぎるものを baseline 側へ戻す**ための引数である。

ここで 1 本の LOS `S[:,y,x]` を見て、連続 True 区間を 1 つ取り、その channel 数を `N_run` とする。代表 channel 幅を `\Delta v` [km/s/ch] とすると、その run の速度幅を

\[
W_{run} = N_{run}\,\Delta v
\]

と定義する。

このとき

\[
W_{run} < \mathrm{lf3d\_min\_run\_kms}
\]

なら、その run 全体を baseline 側へ戻す。

### 19.5.2 これは体積判定ではない

この引数は 3D connected-component の体積判定ではない。

- xy 面積を見ない
- xyv 体積を見ない
- **各 LOS ごとの連続速度幅だけ**を見る

したがって、「非常に短い速度域の点々を落としたい」という目的に対して、最小パラメータで効く post-prune である。

### 19.5.3 実装上の処理

mask は内部で `(nchan, npix)` に reshape され、**速度方向にだけ連結**を取る。

重要なのは、これは

- 別の spatial pixel の run と連結しない
- xy 方向へは連結を見ない

ということである。

よって、たとえば `(x_1,y_1)` の短い run が `(x_2,y_2)` の run と誤って結合されて残ることはない。

### 19.5.4 `None` や `0` のとき

- `None` → 無効
- `<= 0` → 無効

である。

したがって、

```python
lf3d_min_run_kms=None
```

または

```python
lf3d_min_run_kms=0.0
```

なら、この post-prune は行われない。

### 19.5.5 `float` でよいか

はい。単位が km/s なので `float` でよい。

むしろ、

- `dv = 0.2 km/s/ch`
- `min_run_kms = 0.6`

のように、channel 幅の整数倍とは限らない指定を自然に書けるので `float` が適切である。

### 19.5.6 `<` か `<=` か

実装上の判定は

\[
W_{run} < \mathrm{lf3d\_min\_run\_kms}
\]

である。`<=` ではない。

したがって、たとえば

- `dv = 0.2 km/s/ch`
- `lf3d_min_run_kms = 0.6 km/s`

なら

- 1ch run → `0.2` → 除去
- 2ch run → `0.4` → 除去
- 3ch run → `0.6` → **残す**

である。

### 19.5.7 負の `CDELT3` でも大丈夫か

大丈夫である。

使うのは速度軸の増減方向ではなく、

\[
\Delta v = \mathrm{median}(|v_{i+1}-v_i|)
\]

という **絶対幅** だからである。

### 19.5.8 どの段階に効くか

この引数は最終出力の mask だけに効くのではない。

実装上は

1. Stage-A で作った `seed3d`
2. 各 iteration で更新される provisional `line_mask`

の両方に適用される。

したがって、**反復中の auto 判定そのもの**を少し保守的にする効果がある。

### 19.5.9 速度軸が不等間隔のとき

実装上は `\Delta v` に median step を使っているので、厳密には等間隔軸を想定した近似である。

ただし `map_3d` の通常の regridded cube では速度軸は線形であることを前提にしてよいので、通常運用では問題になりにくい。

## 19.6 3 つの引数を一緒に使うとどうなるか

たとえば

```python
lf_cfg = m3d_bsl.LineFreeConfig(
    mask_kind='voxel_3d',
    auto_method='local_3d_poly',
    lf3d_detect_spatial_fwhm_arcsec=350.0,
    lf3d_detect_spectral_width_chan=5,
    lf3d_min_run_kms=0.6,
)
```

と書いたとする。

このときの概念的な処理は次である。

1. 元の cube `D(v,y,x)` を用意する
2. xy 平面に **350 arcsec の Gaussian kernel** を掛ける
3. v 方向に **5ch boxcar smoothing** を掛ける
4. その detection cube から provisional signal mask を作る
5. 各 LOS で、連続速度幅が `0.6 km/s` 未満の短い run を落とす
6. その後、反復と final baseline fit へ進む

ここで、

- 350 arcsec は **kernel 幅**
- 5ch は **smoothing 幅**
- 0.6 km/s は **post-prune の最小 run 幅**

であり、意味が全く異なることに注意する。

## 19.7 実務上の推奨

### 19.7.1 最初の一歩

最初は次のような考え方がよい。

- 空間側: `lf3d_detect_spatial_fwhm_arcsec` は **beam と同程度か、やや小さめ**から試す
- 速度側: `lf3d_detect_spectral_width_chan` は **3 か 5** から試す
- 短い run の除去: `lf3d_min_run_kms` は **2〜3 channel 相当**から試す

たとえば `dv = 0.2 km/s/ch` なら

- `lf3d_detect_spectral_width_chan = 5`
- `lf3d_min_run_kms = 0.4` または `0.6`

が自然な出発点になる。

### 19.7.2 書き方の推奨

実務的には次のように書くと誤読しにくい。

```python
lf_cfg = m3d_bsl.LineFreeConfig(
    mask_kind='voxel_3d',
    auto_method='local_3d_poly',
    lf3d_detect_spatial_fwhm_arcsec=350.0,   # detection cube の xy Gaussian kernel FWHM
    lf3d_detect_spectral_width_chan=5,       # detection cube の v boxcar smoothing 幅
    lf3d_min_run_kms=0.6,                    # 各 LOS でこれ未満の短い run を baseline 側へ戻す
)
```

### 19.7.3 避けた方がよい表現

次のような表現は避けた方がよい。

- `lf3d_detect_spatial_fwhm_arcsec=350  # 350 arcsec beam にする`
- `lf3d_detect_spectral_width_chan=5    # 5ch で binning`

これらは実装の意味を狭義に誤解させやすい。

より正確には

- `350 arcsec の Gaussian kernel で detection 用に空間 smoothing`
- `5ch boxcar で detection 用に速度 smoothing`

である。

## 19.8 最後に一文でまとめる

- `lf3d_detect_spatial_fwhm_arcsec` は **detection cube の xy Gaussian kernel FWHM**
- `lf3d_detect_spectral_width_chan` は **detection cube の v boxcar smoothing 幅**
- `lf3d_min_run_kms` は **各 LOS における短い spectral run を落とす post-prune の閾値**

である。

この 3 つを混同しなければ、`voxel_3d` の detection 側の挙動はかなり追いやすくなる。
