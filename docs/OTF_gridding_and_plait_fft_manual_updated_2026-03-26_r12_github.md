# OTF gridding + FFT/PLAIT basketweave 実装マニュアル（現在の実装版 + 2026-03-26 増補改訂）


## 1. この文書の対象

この文書は、**現在の実装**に基づいて、OTF の gridding から FFT/PLAIT 型 basketweave までを一貫して説明する実務向けマニュアルです。

対象となる主な公開 API は次です。

- `grid_otf_family(...)`
- `coadd_family_cubes(...)`
- `basketweave_cubes(...)`
- `basketweave_fits(...)`
- `run_otf_plait_pipeline(...)`
- `read_otf_bundle(...)`
- `write_otf_bundle(...)`

この文書では、**現在の推奨経路**だけを扱います。過去の設計経緯や旧方式の説明は入れていません。

この版では、kernel に関する記述を **現在の code / spec / kernel 詳細説明書** に合わせて再点検し、特に `sf` の formal default を `convsupport=3`、推奨の 2 本立てを `sf + convsupport=3` と `gjinc + mangum + signed` にそろえています。

---

## 2. 全体像

現在の実装は、次の 3 段階で考えると分かりやすいです。

### 2.1 family ごとの OTF gridding

X 方向群、Y 方向群のように、**走査方向ごとに分けた入力**をそれぞれ gridding し、family cube を作ります。

```text
raw/scantable -> grid_otf_family(...) -> family cube
```

ここでの重みは、各 dump の `BSL_RMS` に基づく inverse-variance weighting が基本です。

### 2.2 family cube の coadd

同じ family に属する既存 cube や追加観測分の cube を、**各スペクトルごとの empirical RMS map** を使って coadd します。

$$
T_{\rm out}(k,y,x)=
\frac{\sum_n w_n(y,x)\,T_n(k,y,x)}{\sum_n w_n(y,x)}
$$

### 2.3 X/Y family から FFT/PLAIT

X family cube と Y family cube を受けて、FFT/PLAIT によって stripe を抑えた最終 cube を作ります。

`X family cube`, `Y family cube` → `basketweave_cubes(...)` → `PLAIT cube`

現在の実装では、FFT 重みモードは `plait_noise_mode` で切り替えます。

- `plait_noise_mode="family_scalar"`
  - X/Y それぞれ 1 個の代表 RMS を使う旧来互換経路
- `plait_noise_mode="family_channel"`
  - X/Y それぞれの channel 別 RMS スペクトル `sigma_ch(k)` を使って channel ごとに Fourier 重みを変える経路

**2026-03-26 増補**: 実装上は `family_scalar` と `family_channel` の両方が利用可能です。`family_scalar` は後方互換・比較用として残し、channel 依存の雲・帯域依存ノイズ・line 近傍の noise mismatch を扱いたいときは `family_channel` を使います。

---

## 3. 原理式

## 3.1 OTF gridding

各 dump を $i$、各 output pixel を $p$、各 channel を $k$ とします。

- dump スペクトル: $S_i(k)$
- gridding kernel: $K_{ip}$
- dump RMS: $\sigma_i$
- dump 重み: $q_i = 1 / \sigma_i^2$

現在の実装での基本形は

$$
T(k,p)=
\frac{\sum_i K_{ip}\,q_i\,S_i(k)}{\sum_i K_{ip}\,q_i}
$$

です。

ここで $\sigma_i$ は通常 `BSL_RMS` です。

現在の実装では、`estimator='avg'` のとき、この平均を**何点以上の dump で支えるか**を `n_min_avg` で制御します。概念的には、各 output pixel $p$ について、実際にその pixel に寄与した有効サンプル数を $N_{\rm eff}(p)$ とすると、

$$
N_{\rm eff}(p) < n_{\rm min,avg}
$$

の pixel は、平均値を採用せず invalid とみなします。したがって、

- `n_min_avg=1`
  - 1 dump でも値を出す
- `n_min_avg=2`
  - 少なくとも 2 dump 必要
- `n_min_avg` を大きくする
  - edge や sparse coverage の pixel は落ちやすくなるが、信頼性は上がる

という意味になります。現在の既定は `n_min_avg=2` です。したがって、synthetic な 1 点だけの入力では、kernel が正常でも出力 pixel が invalid になって NaN になることがあります。これは kernel バグではなく、QC 条件による動作です。

### 3.1.1 row flag と channel 欠損の扱い

現在の実装では、gridding 前の入力 `GridInput` は少なくとも次を持つと考えてください。

- `spec`: shape `(ndump, nchan)` の 2D 配列
- `x`, `y`, `time`, `flag`: shape `(ndump,)` の 1D 配列

ここで `flag` は **dump 行全体を使うかどうか**を表す row mask です。公開仕様としては

- `True` = その dump 行を gridding に使う
- `False` = その dump 行全体を除外する

です。現在の実装は後方互換のため `0/1` や `0.0/1.0` も受け付けますが、意味は同じです。`-1`, `2`, `0.5` のような曖昧な値や、`(ndump, nchan)` の 2D flag は想定していません。

重要なのは、`flag` は **channel ごとの mask ではない**という点です。ある dump の一部 channel だけを無効にしたい場合は、`flag` ではなく `spec[i, k] = NaN` を使います。実務上の役割分担は次です。

- 行全体を落としたい → `flag[i] = False`
- 一部 channel だけ落としたい → `spec[i, bad_channels] = NaN`

この仕様は `core.grid_otf(...)` だけでなく、OTF 入力の merge、basketweave 前段の geometry mask、`ps_gridder` でも同じです。したがって、行単位の採否と channel 欠損を分けて扱うことが、現在の実装を安全に使う上で重要です。

## 3.2 cube coadd

coadd 前に、各 cube の各 spatial pixel $(y,x)$ にある 1 本のスペクトル
$T(:,y,x)$ から、line-free channel 集合 $L$ を使って empirical RMS を計算します。

$$
\sigma_n(y,x)=
{m RMS}_{m robust}
\left(
\{T_n(k,y,x)\mid k\in L,\ \mathrm{finite}\}
\right)
$$

これから 2D 重み map

$$
w_n(y,x)=1/\sigma_n(y,x)^2
$$

を作り、channel ごとに

$$
T_{\rm out}(k,y,x)=
\frac{\sum_n w_n(y,x)\,T_n(k,y,x)}{\sum_n w_n(y,x)}
$$

で coadd します。

## 3.3 FFT/PLAIT

X family, Y family それぞれについて、まず各スペクトルごとの RMS map を作ります。

$$
\sigma_X(y,x),\qquad \sigma_Y(y,x)
$$

### 3.3.1 `family_scalar`

`plait_noise_mode="family_scalar"` では、FFT 重みを簡単化するため、それぞれの map から代表値

$$
\bar\sigma_X = {m median}_{m robust}\{\sigma_X(y,x)\}
$$

$$
\bar\sigma_Y = {m median}_{m robust}\{\sigma_Y(y,x)\}
$$

を作ります。

Fourier 空間 $(u,v)$ での重みは

$$
W_X(u,v)=
\frac{\bar\sigma_Y^2\,u^2}{\bar\sigma_Y^2\,u^2+\bar\sigma_X^2\,v^2}
$$

$$
W_Y(u,v)=
\frac{\bar\sigma_X^2\,v^2}{\bar\sigma_Y^2\,u^2+\bar\sigma_X^2\,v^2}
$$

です。

DC 原点 $(u,v)=(0,0)$ だけは inverse-variance 平均にします。

$$
W_X(0,0)=
\frac{1/\bar\sigma_X^2}{1/\bar\sigma_X^2+1/\bar\sigma_Y^2}
$$

$$
W_Y(0,0)=
\frac{1/\bar\sigma_Y^2}{1/\bar\sigma_X^2+1/\bar\sigma_Y^2}
$$

### 3.3.2 `family_channel`

`plait_noise_mode="family_channel"` では、各 family について line-free channel から channel 別 noise spectrum

$$
\sigma_{X,\mathrm{ch}}(k),\qquad \sigma_{Y,\mathrm{ch}}(k)
$$

を作り、channel ごとに Fourier 重みを作ります。概念的には

$$
W_X(k;u,v)=
\frac{\sigma_{Y,\mathrm{ch}}(k)^2\,u^2}{\sigma_{Y,\mathrm{ch}}(k)^2\,u^2+\sigma_{X,\mathrm{ch}}(k)^2\,v^2}
$$

$$
W_Y(k;u,v)=
\frac{\sigma_{X,\mathrm{ch}}(k)^2\,v^2}{\sigma_{Y,\mathrm{ch}}(k)^2\,u^2+\sigma_{X,\mathrm{ch}}(k)^2\,v^2}
$$

です。DC 原点の扱いは同様に channel ごとの inverse-variance 平均です。

実装上の注意:

- `family_channel` は channel 依存 noise / 雲 / line 近傍悪化を扱うための v2 相当経路です。
- `family_scalar` は後方互換・比較用として残しています。
- diagnostics では `CHANNEL_NOISE`, `CHANNEL_NOISE_X`, `CHANNEL_NOISE_Y`, header の `PLWGMODE` によってどちらの経路を使ったか追跡できます。

---

## 4. 現在の実装での標準オブジェクト: `OTFBundle`

`OTFBundle` は、Python 上で family cube や PLAIT 出力を扱う標準オブジェクトです。

最低限の主要フィールドは次です。

- `data`: science cube, shape `(nchan, ny, nx)`
- `header`: FITS header
- `variance`: variance cube またはそれに準ずる 3D 配列
- `valid_mask`: `(nchan, ny, nx)` か `(ny, nx)` の bool mask
- `support_mask`: `(ny, nx)` の bool mask
- `unit`: 単位
- `family_label`: `'X'`, `'Y'`, `'PLAIT'` など
- `image_ext`: 追加 image ext 相当
- `table_ext`: 追加 table ext 相当
- `meta`: 実行時 metadata

### 4.1 現在の実装でよく使う ext

`grid_otf_family(...)` で作られる bundle には、少なくとも次の ext 系情報が入り得ます。

- `WEIGHT`
- `WEIGHT_SUM`
- `HIT`
- `NSAMP`
- `WSUM`
- `WABS`
- `CANCEL`
- `WREL`
- `MASK`
- `TSYS`
- `TINT`
- `TIME`
- `RMS`
- `NEFF`
- `XEFF`
- `YEFF`
- `BIAS_PIX`

`basketweave_cubes(...)` の出力では、加えて例えば次が入ります。

- `RMS_MAP_X`
- `RMS_MAP_Y`
- `VALID_MASK_AND`
- `VALID_MASK_UNION`
- `WEIGHT_X_FFT`
- `WEIGHT_Y_FFT`
- table ext `CHANNEL_NOISE`

**2026-03-26 増補: BW 出力 / mosaic 対応で重要な ext**

1. **BW 後 provenance ext**
   - 入力 family bundle に存在していた diagnostics は、存在するものだけ suffix 付きで引き継がれます。
   - 代表例: `WEIGHT_X_IN`, `WEIGHT_Y_IN`, `HIT_X_IN`, `HIT_Y_IN`, `NSAMP_X_IN`, `NSAMP_Y_IN`, `RMS_X_IN`, `RMS_Y_IN` など
   - 実装上は `WEIGHT`, `WEIGHT_SUM`, `HIT`, `NSAMP`, `WSUM`, `WABS`, `CANCEL`, `WREL`, `RMS`, `TSYS`, `TINT`, `TIME`, `NEFF`, `XEFF`, `YEFF`, `BIAS_PIX` が入力にあれば `..._X_IN`, `..._Y_IN` として保持されます。

2. **channel noise tables**
   - `CHANNEL_NOISE`: 代表 RMS、line-free channel 数、選択経路、`PLWGMODE` など
   - `CHANNEL_NOISE_X`, `CHANNEL_NOISE_Y`: `family_channel` のときの channel 別 `SIGMA`

3. **mosaic 用共通 ext**
   - `MOSAIC_GAIN`
   - `MOSAIC_RMS_OBS`
   - `MOSAIC_TRUST`
   - `MOSAIC_WEIGHT`
   - table ext `MOSAIC_INFO`

4. **mosaic 後追加 ext**
   - `MOSAIC_INPUT_COUNT`
   - `MOSAIC_WEIGHT_SUM`
   - table ext `MOSAIC_COMBINE_INFO`

重要:

- `WEIGHT_X_FFT`, `WEIGHT_Y_FFT` は **tile 間モザイク重みではなく**、BW の Fourier 空間で X/Y をどう混ぜたかを示す diagnostics です。
- `HIT_X_IN`, `HIT_Y_IN` などは provenance です。最終モザイク重みそのものではありません。

---

## 5. 使うべき時 / 使うべきではない時

## 5.1 FFT/PLAIT を使うべき時

現在の実装とシミュレーション結果から、次のような場合には使う価値があります。

- X / Y のような **異なる走査方向の family** がある
- family ごとに gridding した後に、**方向性のある stripe** が残る
- baseline drift や scanning noise に由来する **走査方向依存の縞** が見える
- coverage が完全一致しなくても、X/Y の両方が同じ target / 同じ grid / 同じ速度軸を共有している

実装検証では、`moderate`, `strong`, `partial_coverage`, `noise_mismatch`, `edge_source`, `ripple_only`, `channel_missing`, `coadd_then_plait` のようなケースで、単純平均に対して有意な改善が見られました。

## 5.2 FFT/PLAIT を使うべきではない時

現在の実装では、少なくとも次の条件では保守的に考えるべきです。

- **きれいな noise-only / stripe がほぼ無い**データ
  - `clean_null` 型では悪化し得ます
- **小さい map**
  - 現在の実装は `min(ny,nx) < min_plait_size_pix` で fallback します
- X/Y の WCS, 速度軸, unit が一致していない
- line-free 範囲をまともに取れない
- そもそも X/Y で別々に gridding していない

### 5.2.1 clean field で悪化する理由

PLAIT は stripe を抑える処理であり、stripe が無い場合でも方向依存の Fourier 重みを掛けるため、**白色雑音だけのケースでは単純平均より不利**になることがあります。これは現在の実装のバグというより、PLAIT の性質です。

### 5.2.2 小さい map

現在の実装では、`min(ny,nx) < min_plait_size_pix` のとき FFT/PLAIT は行わず、**family weighted average** に自動 fallback します。

既定値:

- `min_plait_size_pix = 32`
- `small_map_policy = 'fallback_average'`

---

## 6. line-free 範囲の考え方

現在の実装では、line-free 範囲は次の処理で使います。

- `coadd_family_cubes(...)` の empirical RMS map
- `basketweave_cubes(...)` の X/Y RMS map
- `basketweave_fits(...)` の内部 basketweave
- `run_otf_plait_pipeline(...)` の内部 basketweave

### 6.1 現在の実装での基本方針

- **coadd 前に、その入力 cube から RMS を計算する**
- **basketweave 前に、その X/Y cube から RMS を計算する**
- 過度に複雑な provenance 管理はしない

### 6.2 line-free 範囲の指定法

現在の実装では `linefree_velocity_windows_kms` は、`sd_radio_spectral_fits.ranges.parse_windows()` と同じ入力規約にそろえています。

受け付ける代表例:

```python
linefree_velocity_windows_kms = ['-30:-10', '20:50']
```

または

```python
linefree_velocity_windows_kms = [(-30.0, -10.0), (20.0, 50.0)]
```

文字列指定では、最終的に `parse_windows()` により数値タプル列へ正規化されます。

```python
parse_windows([' -30 : -10 ', '50:20'])
# -> [(-30.0, -10.0), (20.0, 50.0)]
```

現在の実装上の重要点:

- `'min:max'` 形式の文字列リストを受け付けます
- 数値タプル列 `[(min, max), ...]` も受け付けます
- 文字列入力では、`parse_windows()` が **常に `(小さい値, 大きい値)` に並べ替えます**
- 境界は `window_to_mask()` により **両端を含む** として扱われます
- 空文字や `':'` を含まない文字列、数値化できない文字列は `parse_windows()` 側で無視されます

注意:

- 少なくとも 3 channel 以上の line-free channel が必要です
- 線が line-free に混入すると RMS が過大評価されます
- baseline subtraction 後の cube に対して RMS を測る場合も、line-free の選び方は重要です
- 文字列指定に無効な要素が混ざると、その要素だけ無視されることがあります。厳密に管理したい場合は、数値タプルか、整形式の `'min:max'` 文字列だけを使うのが安全です

### 6.3 同系統の窓入力

現在の実装では、`linefree_velocity_windows_kms` のような窓入力は、可能なものから順次 `ranges.py` の共通関数へそろえる方針です。

実務上は、同種の入力では次を原則にしてください。

- 文字列で与えるなら `'min:max'` 形式
- Python 側で明示するなら `[(min, max), ...]` 形式
- 逆順に書いてもよいが、内部では `(小さい値, 大きい値)` に正規化される

---

## 7. support / valid / taper の現在実装

## 7.1 support mask

support は「その pixel を science 領域として使うか」の 2D mask です。

現在の実装では、`WEIGHT`, `WEIGHT_SUM`, `WSUM` のような weight-like map がある場合、**最大値の 10% 未満**の pixel は support から外します。

つまり、単に `weight > 0` ではなく、

- coverage が薄い edge
- ほとんど重みが無い領域

を落とす安全装置が入っています。

## 7.2 valid mask

valid は channel ごとの 3D mask です。

現在の実装では、

- `np.isfinite(data)`
- support
- 必要に応じて weight 条件

を組み合わせて作られます。

## 7.3 taper

現在の実装では taper は **FFT 前処理**であり、RMS 推定には使いません。

重要:

- **RMS 推定は必ず taper 前の data / valid / support で行う**
- taper は FFT に渡す cube にだけ掛ける

### 7.3.1 taper 幅

現在の実装では taper 幅は固定割合ではなく、動的です。

- 最大 5 pixel
- かつ map サイズの 10% 以下

---

## 7.4 source-filling map / edge の注意

このチャットで確認した重要点として、**FFT/PLAIT の Fourier 重み自体が一様な広がった信号を消しているわけではない**一方、実装全体では `support_taper=True` と `apodize=True` が edge で値を減衰させます。

したがって:

- マップ全面に近い広がった成分では、中心はほぼ保たれても edge は弱くなり得ます。
- これは `family_channel` 固有の問題ではなく、`family_scalar` でも同様です。
- source-filling に近い場合は、少なくとも一度は `support_taper=False`, `apodize=False` との比較を行ってください。

実務上の整理:

- stripe 除去が主目的で、edge を mosaic overlap で捨てられるなら既定設定のままでよい
- edge まで flux 保存を見たい場合は、既定設定と no-taper/no-apodize を並べて比較する

## 7.5 `GridInput.flag` と `spec` の NaN

この節は、OTF gridding の入力段で最も間違えやすい点をまとめたものです。現在の実装では、入力 `GridInput` の `flag` と `spec` の NaN は別の役割を持ちます。

### 7.5.1 `flag` は 1D の row mask

`GridInput.flag` は shape `(ndump,)` の 1D 配列です。意味は

- `True` = その dump 行を使う
- `False` = その dump 行を除外する

です。`basketweave` の geometry mask、複数 `GridInput` の merge、`ps_gridder` でも同じ意味で扱われます。

後方互換として `0/1` や `0.0/1.0` も受け付けますが、これは bool row mask の別表現と考えてください。公開仕様としては bool を使うのが安全です。

### 7.5.2 `flag` は channel mask ではない

`flag` で指定できるのは **行全体の採否**だけです。したがって、例えば RFI で一部 channel だけ無効にしたい場合に、`flag` を 2D `(ndump, nchan)` で与える使い方は現在の実装では想定していません。

一部 channel を落としたい場合は、`spec[i, k] = NaN` を使います。現在の実装では finite な channel だけが gridding に寄与するので、row 全体を落とさずに channel 欠損だけを表現できます。

### 7.5.3 実務上の使い分け

- 行全体を捨てたい
  - `flag[i] = False`
- 一部 channel だけ捨てたい
  - `spec[i, bad_channels] = NaN`
- 全行を使う
  - `flag = np.ones(ndump, dtype=bool)`

この役割分担を守ると、`core.grid_otf(...)`, `basketweave`, `ps_gridder` の全てで一貫した挙動になります。

### 7.5.4 `n_min_avg` との関係

`flag=False` の行や、全 channel が NaN の行は、その output pixel に寄与しません。したがって sparse な観測や synthetic test では、kernel の形だけでなく `flag` と `spec` の NaN によって有効サンプル数が減り、`n_min_avg` に引っかかって pixel が invalid になることがあります。

つまり、pixel が NaN になったときに疑う順序は

1. 行全体が `flag=False` で落ちていないか
2. `spec` の NaN が多すぎないか
3. その結果 `n_min_avg` に届いているか

です。kernel を疑う前に、この 3 点を確認してください。

## 8. 現在の推奨ワークフロー

## 8.1 raw/scantable から family cube を作って PLAIT

```python
from sd_radio_spectral_fits.map_3d import (
    grid_otf_family,
    basketweave_cubes,
)

x_bundle = grid_otf_family(
    x_scantables,
    config=cfg,
    family_label='X',
    coord_sys='icrs',
    projection='SFL',
    ref_lon=83.809,
    ref_lat=-5.372639,
    dv_kms=0.2,
    vmin_kms=-30.0,
    vmax_kms=55.0,
    otf_input_state='with_turnarounds',
)

y_bundle = grid_otf_family(
    y_scantables,
    config=cfg,
    family_label='Y',
    coord_sys='icrs',
    projection='SFL',
    ref_lon=83.809,
    ref_lat=-5.372639,
    dv_kms=0.2,
    vmin_kms=-30.0,
    vmax_kms=55.0,
    otf_input_state='with_turnarounds',
)

bw_bundle = basketweave_cubes(
    x_bundle,
    y_bundle,
    linefree_velocity_windows_kms=['-30:-10', '20:50'],
)
```

## 8.2 family cube を FITS に保存してから basketweave

```python
from sd_radio_spectral_fits.map_3d import (
    grid_otf_family,
    basketweave_fits,
)

grid_otf_family(
    x_scantables,
    config=cfg,
    family_label='X',
    output_fits='x_family.fits',
    overwrite=True,
)

grid_otf_family(
    y_scantables,
    config=cfg,
    family_label='Y',
    output_fits='y_family.fits',
    overwrite=True,
)

basketweave_fits(
    'x_family.fits',
    'y_family.fits',
    'bw_output.fits',
    linefree_velocity_windows_kms=['-30:-10', '20:50'],
    overwrite=True,
)
```

## 8.3 既存 cube を coadd してから PLAIT

```python
from sd_radio_spectral_fits.map_3d import (
    read_otf_bundle,
    coadd_family_cubes,
    basketweave_cubes,
)

x1 = read_otf_bundle('x_day1.fits')
x2 = read_otf_bundle('x_day2.fits')
y1 = read_otf_bundle('y_day1.fits')
y2 = read_otf_bundle('y_day2.fits')

x_master = coadd_family_cubes(
    [x1, x2],
    linefree_velocity_windows_kms=['-30:-10', '20:50'],
)

y_master = coadd_family_cubes(
    [y1, y2],
    linefree_velocity_windows_kms=['-30:-10', '20:50'],
)

bw = basketweave_cubes(
    x_master,
    y_master,
    linefree_velocity_windows_kms=['-30:-10', '20:50'],
)
```

## 8.4 1 本で流す

```python
from sd_radio_spectral_fits.map_3d import run_otf_plait_pipeline

bw_bundle = run_otf_plait_pipeline(
    x_inputs=x_scantables,
    y_inputs=y_scantables,
    config=cfg,
    linefree_velocity_windows_kms=['-30:-10', '20:50'],
    coord_sys='icrs',
    projection='SFL',
    ref_lon=83.809,
    ref_lat=-5.372639,
    dv_kms=0.2,
    vmin_kms=-30.0,
    vmax_kms=55.0,
    output_fits='bw_output.fits',
    overwrite=True,
)
```

---

## 8.5 tile ごとに BW してから mosaic する

現在の実装とこのチャットでの最終方針では、広域観測は **tile ごとに OTF gridding / 必要なら coadd / 必要なら BW** を行い、その後で mosaic する使い方を推奨します。

理由:

- 観測日・気象・stripe 条件が tile ごとに違うと、巨大 1 枚の一括 BW より tile ごとの BW の方が素直
- edge の違和感は overlap を持った mosaic で吸収する方が現実的
- メモリ・FFT サイズの面でも tile 分割が有利

標準的な流れ:

1. `grid_otf_family(...)` で X/Y family cube を作る
2. 必要なら `coadd_family_cubes(...)` で同一 family をまとめる
3. `basketweave_cubes(...)` または `basketweave_fits(...)` で tile ごとに BW を行う
4. 各 tile の `MOSAIC_*` ext を用いて `mosaic_bundles(...)` / `mosaic_fits(...)` で後段合成する

BW が不要な tile は、そのまま `grid_otf_family(...)` 出力を mosaic 入力に使います。このとき `MOSAIC_GAIN=1` と扱われるので、BW あり/なしで同一の mosaic API を使えます。

## 9. 関数ごとの詳細説明

## 9.1 `grid_otf_family(...)`

### 役割

- 1 つ以上の scantable を受ける
- それらを 1 family としてまとめる
- velocity axis を揃える
- plane projection に落とす
- `grid_otf(...)` を呼ぶ
- `GridResult` を `OTFBundle` に変換する
- 必要なら FITS に書く

### シグネチャ

```python
grid_otf_family(
    scantables,
    config: MapConfig,
    *,
    family_label: str,
    coord_sys: str = 'icrs',
    projection: str = 'SFL',
    out_scale: str = 'TA*',
    dv_kms: float | None = None,
    vmin_kms: float | None = None,
    vmax_kms: float | None = None,
    linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None = None,
    ref_coord=None,
    ref_lon: float | None = None,
    ref_lat: float | None = None,
    reproducible_mode: bool | None = None,
    workers: int | None = None,
    sort_neighbors: bool | None = None,
    verbose: bool | None = None,
    otf_input_state=None,
    otf_scan_region=None,
    otf_scan_png=None,
    existing_turn_labels: str | None = None,
    otf_scan_existing_is_turn: str | None = None,
    output_fits: str | None = None,
    overwrite: bool = False,
)
```

### 引数説明

#### `scantables`
- 1 個または複数個の scantable
- list/tuple でもよい
- 実装内で 1 family として結合される

#### `config`
- `MapConfig` インスタンス
- spatial grid, kernel, weighting などを決める

#### `family_label`
- `'X'`, `'Y'` などの family 名
- header と bundle metadata に入る

#### `coord_sys`
- 投影前の座標系
- 既定 `'icrs'`

#### `projection`
- 2D plane projection
- 既定 `'SFL'`

#### `out_scale`
- 出力温度尺度
- 現在は `'TA*'` または `'TR*'` を想定

#### `dv_kms`, `vmin_kms`, `vmax_kms`
- velocity regrid の設定
- `Standardizer.get_matrix(...)` に渡される
- `None` の場合は `config` 側の `dv_kms`, `vmin_kms`, `vmax_kms` を優先

#### `linefree_velocity_windows_kms`
- **2026-03-26 増補**
- `MOSAIC_RMS_OBS` / `MOSAIC_WEIGHT` を作るときに使う line-free 窓
- `'min:max'` 形式の文字列リスト、または数値タプル列 `[(min, max), ...]`
- `None` なら variance / RMS ext fallback を使う

#### `ref_coord`, `ref_lon`, `ref_lat`
- 参照座標
- `ref_coord` を与えればそれを使う
- そうでなければ `ref_lon`, `ref_lat`
- さらに無ければ scantable 群から解決

#### `reproducible_mode`, `workers`, `sort_neighbors`, `verbose`
- `MapConfig` の runtime override
- `None` なら `config` 側の値を使う

#### `otf_input_state`
OTF scan 領域の扱いを決める。

代表的には次を使う。

- `'with_turnarounds'`
  - turnaround を含める
- `'scan_only'`
  - turnaround を除く
- `'use_existing_labels'`
  - 既存の `SCAN` / `IS_TURN` を使う

#### `otf_scan_region`
- scan 領域を明示したいときの追加指定

#### `otf_scan_png`
- scan 領域診断の PNG 出力
- `True` なら自動命名
- `str` ならそのパス

#### `existing_turn_labels`, `otf_scan_existing_is_turn`
- 既存ラベル利用時の挙動調整

#### `output_fits`
- 与えると bundle をそのまま FITS 出力する

#### `overwrite`
- FITS 上書き可否

### 出力
- `OTFBundle`
- **2026-03-26 増補**: 戻り値 bundle には `attach_mosaic_products(...)` が自動で適用され、少なくとも `MOSAIC_GAIN=1`, `MOSAIC_TRUST=1`、および line-free/variance 由来の `MOSAIC_RMS_OBS`, `MOSAIC_WEIGHT`, `MOSAIC_INFO` が追加されます。

### 現在の実装で重要な点
- family 内の複数 scantable は最初にマージされる
- `BSL_RMS` 重みは upstream `grid_otf(...)` の責務
- `linefree_velocity_windows_kms` はこの関数ではまだ使わない
- `baseline_subtracted` は `False` で bundle meta に入る
- `GridInput.flag` は 1D row mask として扱われる
  - `True` / `1` / `1.0` = 採用
  - `False` / `0` / `0.0` = 除外
- dump の一部 channel を落とす場合は `flag` ではなく `spec=NaN` を使う

---

## 9.2 `coadd_family_cubes(...)`

### 役割

- 同一 family の bundle 群を coadd する
- 各 bundle から **各スペクトルごとの 2D empirical RMS map** を計算する
- その RMS map に基づく 2D 重みで channel ごとに平均する

### シグネチャ

```python
coadd_family_cubes(
    bundles,
    *,
    linefree_velocity_windows_kms,
    strict_shape: bool = True,
    strict_wcs: bool = True,
) -> OTFBundle
```

### 引数説明

#### `bundles`
- `OTFBundle` または list/tuple of `OTFBundle`
- 同じ family の cube である必要がある

#### `linefree_velocity_windows_kms`
- empirical RMS map 計算に使う line-free 窓
- **必須**
- `'min:max'` 形式の文字列リスト、または数値タプル列 `[(min, max), ...]` を受け付ける
- 文字列入力は `sd_radio_spectral_fits.ranges.parse_windows()` で正規化される

#### `strict_shape`
- `True` なら shape mismatch を error

#### `strict_wcs`
- `True` なら WCS mismatch を error

### 出力
- coadd 後の `OTFBundle`

### 現在の実装での重み

各入力 bundle $n$ に対して、

$$
\sigma_n(y,x)={m RMS}_{m robust}\{T_n(k,y,x)\mid k\in L\}
$$

から

$$
w_n(y,x)=1/\sigma_n(y,x)^2
$$

を作って coadd する。

### 現在の ext / meta

- image ext:
  - `HIT_COUNT_COADD`
  - `RMS_MAP_EMP`
  - `VALID_MASK_UNION`
  - `RMS_MAP_INPUT_i`
  - **2026-03-26 増補**: `MOSAIC_GAIN`, `MOSAIC_RMS_OBS`, `MOSAIC_TRUST`, `MOSAIC_WEIGHT`
- table ext:
  - `COADD_INFO`
  - **2026-03-26 増補**: `MOSAIC_INFO`
- meta:
  - `coadd_n_input`
  - `linefree_velocity_windows_kms`
  - `baseline_subtracted`
  - `noise_mode='empirical_rms_per_spectrum'`
  - **2026-03-26 増補**: `mosaic_*` 系 metadata

### 実務上の注意
- baseline 後 cube を coadd してもよい
- その場合も、**coadd 前にその cube から line-free RMS を再計算**する
- つまり、古い RMS を引き継ぐより、その時点の入力から測り直す方針

---

## 9.3 `basketweave_cubes(...)`

### 役割

- X bundle と Y bundle を受ける
- line-free 範囲から X/Y の RMS map を計算する
- `plait_noise_mode` に応じて
  - `family_scalar`: X/Y の代表 RMS
  - `family_channel`: channel 別 RMS スペクトル
  を使う
- taper / padding / FFT / 重み付け / inverse FFT を行う
- `OTFBundle` として返す
- **2026-03-26 増補**: BW 出力に入力 family provenance ext と `MOSAIC_*` ext を付与する

### シグネチャ

```python
basketweave_cubes(
    x_bundle: OTFBundle,
    y_bundle: OTFBundle,
    *,
    linefree_velocity_windows_kms,
    output_fits: str | None = None,
    overwrite: bool = False,
    **kwargs,
) -> OTFBundle
```

内部では `plait_fft_cubes(...)` を呼ぶ。

### `plait_fft_cubes(...)` の主要パラメータ

```python
plait_fft_cubes(
    x_bundle,
    y_bundle,
    *,
    linefree_velocity_windows_kms,
    noise_mode: str = 'empirical_rms',
    plait_noise_mode: str = 'family_scalar',
    pad_frac: float = 0.25,
    apodize: bool = True,
    apodize_alpha: float = 0.1,
    support_taper: bool = True,
    support_taper_width_pix: int | None = None,
    science_mask_mode: str = 'and',
    fft_workers: int | None = None,
    dtype: str | np.dtype = np.float64,
    diagnostics: bool = True,
    min_plait_size_pix: int = 32,
    small_map_policy: str = 'fallback_average',
    quality_gate_mode: str = 'none',
    min_improvement_frac: float = 0.0,
)
```

### 引数説明

#### `linefree_velocity_windows_kms`
- X/Y の RMS map 計算に使う line-free 窓
- **必須**
- `'min:max'` 形式の文字列リスト、または数値タプル列 `[(min, max), ...]` を受け付ける
- 文字列入力は `sd_radio_spectral_fits.ranges.parse_windows()` で正規化される

#### `noise_mode`
- 現在の実装は `'empirical_rms'` のみ対応

#### `plait_noise_mode`
- **2026-03-26 増補**
- `'family_scalar'` または `'family_channel'`
- `family_scalar`:
  - X/Y それぞれ 1 個の代表 RMS を使う
  - 旧来互換・比較用
- `family_channel`:
  - channel 別 `sigma_ch(k)` を作って channel ごとに Fourier 重みを変える
  - 雲・帯域依存ノイズ・line 近傍悪化に対して v2 相当の改善を狙う

#### `pad_frac`
- spatial padding 量
- 既定 `0.25`

#### `apodize`
- Tukey 型 apodization を掛けるか

#### `apodize_alpha`
- apodization の強さの基準値
- 実際の幅は動的制御される

#### `support_taper`
- support 境界に沿った taper を掛けるか

#### `support_taper_width_pix`
- support taper の幅
- `None` のとき動的決定

#### `science_mask_mode`
- `'and'` または `'union'`
- 出力 valid/support の作り方

`'and'`:
- X/Y の両方が有効な pixel のみを science 領域とする

`'union'`:
- どちらか一方が有効なら science 領域とする

#### `fft_workers`
- SciPy FFT の worker 数

#### `dtype`
- FFT に使う型
- 既定 `float64`

#### `diagnostics`
- diagnostics ext / table を出すか

#### `min_plait_size_pix`
- FFT/PLAIT を許す最小 map サイズ
- 既定 `32`

#### `small_map_policy`
- 小マップ時の挙動
- 現在主に `'fallback_average'` を想定

#### `quality_gate_mode`
- 既定 `'none'`
- optional quality gate
- 実装はあるが既定では無効

#### `min_improvement_frac`
- quality gate の閾値

### 出力
- `OTFBundle`
- `family_label = 'PLAIT'`

### 現在の ext / meta

#### image ext
- `RMS_MAP_X`
- `RMS_MAP_Y`
- `VALID_MASK_AND`
- `VALID_MASK_UNION`
- `WEIGHT_X_FFT`
- `WEIGHT_Y_FFT`
- 入力 family provenance ext (`WEIGHT_X_IN`, `WEIGHT_Y_IN`, `HIT_X_IN`, `HIT_Y_IN`, `NSAMP_X_IN`, `NSAMP_Y_IN`, `RMS_X_IN`, `RMS_Y_IN` など。入力に存在するものだけ)
- **2026-03-26 増補**: `MOSAIC_GAIN`, `MOSAIC_RMS_OBS`, `MOSAIC_TRUST`, `MOSAIC_WEIGHT`

#### table ext
- `CHANNEL_NOISE`
  - `PLWGMODE`
  - `SIGMA_X_REP`
  - `SIGMA_Y_REP`
  - `LINEFREE_NCHAN`
  - `SIGMA_OUT_PLAIT`
  - `SIGMA_OUT_AVG`
  - `SELECTED`
  - `NOISE_PIX_X`, `NOISE_PIX_Y`
- `CHANNEL_NOISE_X`, `CHANNEL_NOISE_Y`
  - `CHAN`, `VELOCITY_KMS`, `LINEFREE`, `SIGMA`
- **2026-03-26 増補**: `MOSAIC_INFO`

#### meta
- `basketweave_method='plait_fft'` または fallback 名
- `noise_mode='empirical_rms'`
- `plait_noise_mode`
- `linefree_velocity_windows_kms`
- `baseline_subtracted`
- `sigma_x_rep`
- `sigma_y_rep`
- `selection_method`
- `small_map_fallback`
- `mosaic_*` 系 metadata

### 実務上の注意
- RMS 推定は **taper 前**
- taper は FFT 用だけ
- support は weight threshold で切られる
- `clean_null` では悪化し得る
- small map は fallback する

---

## 9.4 `basketweave_fits(...)`

### 役割

- X family FITS と Y family FITS を読む
- `basketweave_cubes(...)` を呼ぶ
- 出力 FITS を書く

### シグネチャ

```python
basketweave_fits(
    x_fits: str,
    y_fits: str,
    output_fits: str,
    *,
    linefree_velocity_windows_kms,
    overwrite: bool = False,
    **kwargs,
) -> OTFBundle
```

### 実務上の注意
- `x_fits`, `y_fits` は **すでに同じ grid / 同じ速度軸 / 同じ unit** である必要がある
- そうでない場合は reader 後に `basketweave_cubes(...)` 側で error になる

---

## 9.5 `run_otf_plait_pipeline(...)`

### 役割

- X family の gridding
- Y family の gridding
- そのまま FFT/PLAIT

を 1 本で行う高位関数。

### シグネチャ

```python
run_otf_plait_pipeline(
    x_inputs,
    y_inputs,
    config,
    *,
    linefree_velocity_windows_kms,
    output_fits: str | None = None,
    x_family_label: str = 'X',
    y_family_label: str = 'Y',
    x_output_fits: str | None = None,
    y_output_fits: str | None = None,
    overwrite: bool = False,
    grid_kwargs: dict | None = None,
    plait_kwargs: dict | None = None,
    **legacy_kwargs,
) -> OTFBundle
```

### 引数説明

#### `x_inputs`, `y_inputs`
- X family, Y family の raw/scantable 群

#### `config`
- `MapConfig`

#### `linefree_velocity_windows_kms`
- FFT/PLAIT の RMS 推定用 line-free 窓
- `'min:max'` 形式の文字列リスト、または数値タプル列 `[(min, max), ...]` を受け付ける
- 文字列入力は `sd_radio_spectral_fits.ranges.parse_windows()` で正規化される

#### `output_fits`
- 最終 PLAIT cube の出力先

#### `x_family_label`, `y_family_label`
- family 名

#### `x_output_fits`, `y_output_fits`
- 中間 family cube を書きたいときの出力先

#### `overwrite`
- FITS 上書き可否

#### `grid_kwargs`
- `grid_otf_family(...)` に渡す dict

#### `plait_kwargs`
- `basketweave_cubes(...)` / `plait_fft_cubes(...)` に渡す dict

#### `**legacy_kwargs`
- 実装内部で key によって `grid_kwargs` と `plait_kwargs` に振り分けるための後方互換的入口
- 新規コードでは、可能なら `grid_kwargs`, `plait_kwargs` を明示した方が分かりやすい

---

## 9.6 `read_otf_bundle(...)` / `write_otf_bundle(...)`

## `write_otf_bundle(...)`

### 役割
- `OTFBundle` を FITS に書く

### シグネチャ

```python
write_otf_bundle(bundle: OTFBundle, path: str, *, overwrite: bool = False) -> None
```

### 現在の書き方
- PRIMARY = `data`
- `VARIANCE`
- `SUPPORT_MASK`
- `VALID_MASK`
- その他 `image_ext`
- `table_ext`

## `read_otf_bundle(...)`

### 役割
- FITS から `OTFBundle` を復元する

### シグネチャ

```python
read_otf_bundle(path: str) -> OTFBundle
```

### 現在の読み方
- PRIMARY を `data`
- `VARIANCE`, `VALID_MASK`, `SUPPORT_MASK` は専用フィールドへ
- その他 image ext は `image_ext`
- bin table は `table_ext`

---

## 9.7 `mosaic.py` / `attach_mosaic_products(...)` / `mosaic_bundles(...)` / `mosaic_fits(...)`

### 役割

- `attach_mosaic_products(...)`: 1 つの `OTFBundle` から `MOSAIC_GAIN`, `MOSAIC_RMS_OBS`, `MOSAIC_TRUST`, `MOSAIC_WEIGHT`, `MOSAIC_INFO` を作る
- `mosaic_bundles(...)`: 複数 `OTFBundle` を gain-aware に合成する
- `mosaic_fits(...)`: FITS を読んで `mosaic_bundles(...)` を呼び、必要なら FITS に書く wrapper

### 現在の標準式

各 tile `m` について保存されている値を `D_m(k,y,x)`、observed RMS を `sigma_obs,m(y,x)`、mosaic gain を `g_m(y,x)`、trust を `t_m(y,x)` とします。現在のモザイクでは

$$
N_m(k,y,x)=t_m(y,x)\,\frac{g_m(y,x)}{\sigma_{m,\mathrm{obs}}(y,x)^2}\,D_m(k,y,x)
$$

$$
W_m(y,x)=t_m(y,x)\,\frac{g_m(y,x)^2}{\sigma_{m,\mathrm{obs}}(y,x)^2}
$$

を作り、

$$
S_{\mathrm{mosaic}}(k,y,x)=\frac{\sum_m N_m(k,y,x)}{\sum_m W_m(y,x)}
$$

で合成します。

ここで重要なのは、`MOSAIC_GAIN` は**重みではなく保存値に掛かっている multiplicative response** だという点です。したがって `MOSAIC_WEIGHT=1/rms^2` だけでは不十分です。現在の実装で保存している `MOSAIC_WEIGHT` は、記号で書けば

$$
W_{\rm mosaic}=t_{\rm mosaic}\times\frac{g_{\rm mosaic}^2}{\sigma_{\rm mosaic}^2}
$$

に対応します。ここで

- $t_{\rm mosaic}$ は `MOSAIC_TRUST`
- $g_{\rm mosaic}$ は `MOSAIC_GAIN`
- $\sigma_{\rm mosaic}$ は `MOSAIC_RMS_OBS`

です。

### `MOSAIC_*` ext の意味

- `MOSAIC_GAIN`
  - OTF gridding / coadd だけの tile では 1
  - BW tile では BW edge attenuation の 2D response
- `MOSAIC_RMS_OBS`
  - **その時点の最終 cube** から line-free で再計算した observed RMS
- `MOSAIC_TRUST`
  - 位置ずれや既知の怪しい領域を将来 heuristic に落とすための係数
  - 現在の標準では 1
- `MOSAIC_WEIGHT`
  - `trust * gain^2 / rms_obs^2`
- `MOSAIC_INFO`
  - `RMS_SOURCE`, `LINEFREE_NCHAN`, `WINDOWS`, `GAIN_SOURCE`, `TRUST_SOURCE`, `GAIN_MIN_USED`, `WEIGHT_FORMULA`, `NUMERATOR_FORMULA`

### `gain_min`

- 現在の既定は `gain_min=0.5`
- `gain < 0.5` の pixel は mosaic で使いません
- 理由: `D/g` 補正後の RMS が 2 倍を超える edge は無理に救わず、overlap 領域の隣接 tile に任せる方が安全だからです

### line-free 窓を変えたときの規約

重要:

- 後段 `mosaic_bundles(...)` / `mosaic_fits(...)` で **新しい** `linefree_velocity_windows_kms` を与えたときは、既存の `MOSAIC_WEIGHT` をそのまま使い回してはいけません。
- 現在の実装では、その都度 `attach_mosaic_products(..., overwrite=True)` を通し、**保存済み `MOSAIC_GAIN` / `MOSAIC_TRUST` は再利用しつつ、`MOSAIC_RMS_OBS` / `MOSAIC_WEIGHT` は current line-free 窓で再計算**します。

### 現在の適用範囲

- 同一 shape
- 同一 WCS
- 同一 unit
- reprojection / union-WCS はまだ未対応

### 出力 mosaic bundle

- `family_label="MOSAIC"`
- `MOSAIC_INPUT_COUNT`, `MOSAIC_WEIGHT_SUM`, `MOSAIC_COMBINE_INFO` を追加
- 出力 mosaic 自身には `attach_mosaic_products(...)` を再適用し、`MOSAIC_GAIN=1`, `MOSAIC_TRUST=1` の最終 product として整えます

### 実務上の注意

- `WEIGHT`, `WEIGHT_SUM`, `HIT`, `WEIGHT_X_FFT`, `WEIGHT_Y_FFT` を tile 間 mosaic 重みに直接使わないでください。これらは provenance / diagnostics です。
- `MOSAIC_WEIGHT` にさらに `1/rms^2` を掛けないでください。
- 分母和 `W_sum` が 0 の pixel は `NaN` とし、NumPy 実装では割る前に `valid = (W_sum > 0)` を切ってください。

## 9.8 verbose beam summary のローカル拡張（別パッチ）

このチャットでは、別パッチとして次の表示改善も提案しました。これは本体の数式には影響しませんが、適用している場合は便利です。

- `grid_otf_family(..., verbose=True)` で kernel summary, nominal effective beam, empirical center beam を表示
- `basketweave_cubes(..., verbose=True)` で最終出力 bundle に継承された beam 情報を 1 回だけ表示

注意:

- BW 後の beam 表示は **FFT/PLAIT 後に新しく empirical beam を再推定したものではなく**、family 側の beam 情報を継承して表示するだけです。
- 真に BW 後の中心 beam を再推定したいなら、flat-field / point response を別途流す実装が必要です。

## 10. `MapConfig` の説明（このワークフローで重要なもの）

`grid_otf_family(...)` では既存の `MapConfig` をそのまま使います。ここでは OTF gridding + FFT/PLAIT ワークフローで関係が深い順に説明します。

## 10.1 空間 grid 定義

- `x0`, `y0`
  - map 中心オフセット [arcsec]
- `nx`, `ny`
  - pixel 数
- `cell_arcsec`
  - pixel size [arcsec]
- `beam_fwhm_arcsec`
  - telescope beam FWHM [arcsec]

これらは X family / Y family で**完全に同じ**にする必要があります。

実務上は、`beam_fwhm_arcsec` を先に決め、そのうえで **`cell_arcsec ≈ beam_fwhm_arcsec / 3`** を最初の標準値にしてください。kernel の public preset は pixel-based なので、`cell_arcsec` は gridding kernel の物理幅に直接効きます。

## 10.2 kernel 関連

現在の実装では、OTF gridding に使う public kernel は

- `kernel='sf'`
- `kernel='gjinc'`
- `kernel='gauss'`

です。

このうち、実務上まず考えるべきなのは **`sf` と `gjinc`** です。`gauss` は比較用・単純参照用としては有用ですが、通常の標準推奨 kernel として最初に選ぶものではありません。

### 10.2.1 まず最初に知っておくべき推奨

OTF gridding の kernel は、最終 map の

- 空間分解能
- 点源や細い構造の peak response
- ringing（点源周囲の波紋状アーティファクト）
- aliasing 抑制（高空間周波数の折り返し偽構造の抑制）
- ノイズの平均化のされ方
- 欠損・edge・粗い sampling に対する頑健さ

に直接効きます。

このパッケージで、まずユーザーに推奨する設定は次の 2 つです。

#### A. 標準推奨（まず安全に全体像を見る）

```python
cfg = MapConfig(
    ...,
    beam_fwhm_arcsec=B,
    cell_arcsec=B / 3,
    kernel='sf',
    convsupport=3,
)
```

意味:

- `cell_arcsec` は **入力ビーム FWHM の 1/3 を目安**にする
- `sf` は positive-only の spheroidal kernel
- `convsupport=3` は現在の formal default
- ringing が小さく、aliasing 抑制と欠損・edge・粗い sampling に対する頑健さを重視する設定である
- `gjinc + mangum` より broad だが、まず破綻しにくい

#### B. 分解能重視の再イメージング

```python
cfg = MapConfig(
    ...,
    beam_fwhm_arcsec=B,
    cell_arcsec=B / 3,
    kernel='gjinc',
    kernel_preset='mangum',
    kernel_sign='signed',
)
```

意味:

- `gjinc` のうち、Mangum 系の signed kernel を使う
- effective beam の太りと point-source peak の低下をできるだけ抑え、細い構造を sharp に出したいときに向く
- 一方で、signed kernel なので noise amplification、cancel、ringing の評価が必要

現在の内部比較では、代表条件 `B=350 arcsec`, `cell=B/3` に対して nominal output PSF FWHM は概ね

- `gjinc + mangum` : 約 **383 arcsec**（入力ビーム 350 arcsec に対して **+9 % 程度**）
- `gjinc + casa`   : 約 **401 arcsec**（**+14 % 程度**）
- `sf + convsupport=3` : 約 **426 arcsec**（**+22 % 程度**）
- `sf + convsupport=6` : 約 **563 arcsec**（**+61 % 程度**）

でした。したがって、

- **まずは安全に作る**なら `sf + convsupport=3`
- **少しでも sharpness を詰めたい**なら `gjinc + mangum + signed`

という 2 段階運用が実務上分かりやすいです。

ここでの数値は、本パッケージ内部での nominal comparison です。kernel 理論、`SF` / `GJINC` の数式、`convsupport` / `truncate` / `first null` の意味、CASA / casacore との対応、未実装項目、内部比較表の詳細は、**別文書の kernel 詳細説明書**を参照してください。この OTF + FFT/PLAIT マニュアルでは、解析ユーザーが workflow 上すぐ使える要点だけに絞って書きます。

### 10.2.2 用語の整理

kernel 関連では、次の語を混同しないでください。

- **support radius（有効半径）**
  - 一般的な意味で、kernel を実際に計算に使う最大距離
- **truncate**
  - 主に `gauss` / `gjinc` で使う support radius の名前
- **convsupport**
  - 主に `sf` で使う support radius の名前で、**pixel 単位**
- **first null**
  - `gjinc` の主ローブが最初に 0 になる半径。`gjinc + casa` では自然な cutoff 候補

この manual では、一般論としては **support radius（有効半径）**、`sf` の public パラメータとしては **`convsupport`**、`gjinc` / `gauss` の cutoff 指定としては **`truncate`** という語を使います。

また、次の比較量を意識すると kernel の役割が分かりやすくなります。

- **effective beam FWHM**
  - 出力マップで実際にどれだけビームが太るか
- **peak response**
  - 点源や細い構造のピークがどれだけ残るか
- **ringing**
  - 強い点源や急な境界の周囲に出る波紋状アーティファクト
- **aliasing 抑制**
  - 表現できない高空間周波数が低周波側へ折り返して偽構造になるのをどれだけ防げるか

### 10.2.3 現在の public パラメータ

- `kernel`
  - `'sf'`, `'gjinc'`, `'gauss'`
- `kernel_preset`
  - `gjinc` に対して主に `'mangum'` または `'casa'`
  - 後方互換 alias として `'mangum2007' -> 'mangum'`, `'legacy' -> 'casa'`
- `kernel_sign`
  - `'auto'`, `'signed'`, `'positive_only'`
- `convsupport`
  - `sf` の public support radius パラメータ
  - 現在の formal default は `3`
  - 現在実装では integer 的な値を想定し、`convsupport < 3` は warning 対象
- `gwidth_pix`, `gwidth_beam`
- `jwidth_pix`, `jwidth_beam`
- `truncate`
- `support_radius_pix`, `support_radius_beam`

重要な現実装上の注意:

- `sf` では、public には `convsupport` を使います
- `sf` に対して `support_radius_pix`, `support_radius_beam`, `truncate` を与える使い方は、現在実装では禁止です
- `gjinc` では `kernel_preset`, `kernel_sign`, `jwidth_*`, `gwidth_*`, `truncate`, `support_radius_*` が意味を持ちます
- `gauss` と `sf` は正の kernel なので、`kernel_sign='signed'` を明示しても、現在実装では `positive_only` と同等に扱います

通常は、explicit width を細かくいじる前に、まず

- `sf, convsupport=3`
- `gjinc, kernel_preset='mangum', kernel_sign='signed'`

のどちらかから始めるのが安全です。

### 10.2.4 `cell_arcsec` と `beam_fwhm_arcsec`

kernel を理解するときに最も間違えやすいのは、`beam_fwhm_arcsec` と `cell_arcsec` の役割です。

- `beam_fwhm_arcsec`
  - 入力望遠鏡ビームの FWHM [arcsec]
- `cell_arcsec`
  - 出力 map の pixel size [arcsec/pixel]

現在の kernel 実装は public preset としては **pixel-based** です。したがって kernel の物理幅は `cell_arcsec` を通して決まります。一方、最終的な nominal effective beam は

- 入力ビーム
- gridding kernel の物理幅

の両方で決まります。

このため、`cell_arcsec` は小さすぎても大きすぎてもよくありません。実務上の第一近似として、**`cell_arcsec ≈ beam_fwhm_arcsec / 3`** を基本にしてください。現在実装では `cell_arcsec=None` のとき、この目安が既定として解決されます。

### 10.2.5 `sf` と `gjinc` と `gauss` の簡単な使い分け

#### `kernel='sf', convsupport=3`

- positive-only
- ringing が小さい
- aliasing 抑制と欠損・edge・粗い sampling に対する頑健さを重視する
- signed kernel より安全で、まず全体像を作りやすい
- ただし `gjinc + mangum` よりは broad

#### `kernel='gjinc', kernel_preset='mangum', kernel_sign='signed'`

- Gaussian-tapered jinc
- signed kernel を使い、高空間周波数保持と peak response の維持を狙う
- このパッケージ内比較では最も sharp
- ただしノイズ増幅や負ローブ由来の ringing、cancel に注意

#### `kernel='gjinc', kernel_preset='casa'`

- first-null cutoff の safe な GJINC
- `mangum` より broad
- `legacy` は現在この `casa` の互換 alias として読む

#### `kernel='gauss'`

- 単調な positive smoothing kernel
- ringing は小さく、強い点源周囲の負アーティファクトを避けやすい
- その代わり、高空間周波数を素直に落としやすく、分解能保持の意味では `gjinc` より弱い
- `sf(convsupport=3)` と比べても分解能の太りだけでは大差ない場合があるが、default としては `sf` の方が support 制御と役割の説明が明確

#### `kernel='sf', convsupport=6`

- 背景的には重要な値だが、現在の内部比較では broadening がかなり大きい
- formal default ではなく、明示的に broadening を許容する場合に限って使う

### 10.2.6 FFT/PLAIT との関係

FFT/PLAIT 側は **family cube がすでに同じ kernel で gridding されていること**を前提にします。したがって X family / Y family の両方で、

- `beam_fwhm_arcsec`
- `cell_arcsec`
- `kernel`
- `kernel_preset` / `kernel_sign` / `convsupport`

は揃えてください。

また、`grid_otf_family(..., verbose=True)` では kernel summary, nominal effective beam, empirical center beam が表示されます。kernel を切り替えたときは、この表示も確認してください。


## 10.3 chunk / dtype

- `chunk_ch`
  - channel chunk size
- `dtype`
  - gridding の内部型

## 10.4 重み付け

- `weight_mode`
  - `'uniform'` or `'rms'`
- `alpha_rms`
- `beta_tint`
- `weight_clip_quantile`
- `weight_clip_max`
- `exclude_turnaround`

通常は **`weight_mode='rms'`** が推奨です。

### 10.4.1 `GridInput.flag` の前提

`weight_mode` や `exclude_turnaround` を評価する前段で、入力 row 自体を `flag` で落とすことができます。現在の実装では `flag` は 1D bool row mask と考えるのが基本です。

- `True` = 使う
- `False` = 落とす

したがって、turnaround を除くかどうかは `exclude_turnaround`、観測条件や upstream 判定で row 自体を使うかどうかは `flag`、という役割分担になります。

## 10.5 estimator / QC

- `estimator`
  - `'avg'` or `'plane'`
- `n_min_avg`
- `n_min_plane`
- `cond_max`
- `dr_eff_warn_pix`
- `eps_u0`
- `eps_weight_sum`
- `min_abs_weight_ratio`
- `min_cancel_ratio`

### 10.5.1 `n_min_avg` とは何か

`n_min_avg` は、`estimator='avg'` のときに、ある output pixel を**平均値として採用するために必要な最小寄与サンプル数**です。kernel がその pixel に届いていても、実際に寄与した dump 数が少なすぎる場合は、その pixel を invalid とみなします。

実務上は次のように理解してください。

- `n_min_avg=1`
  - 1 dump でも値を出す
  - coverage を広く残したい確認用には便利
  - ただし局所的に不安定になりやすい
- `n_min_avg=2`
  - 少なくとも 2 dump 必要
  - 現在の既定であり、通常はこちらが自然
- `n_min_avg` をさらに大きくする
  - edge や sparse coverage の pixel は落ちやすくなる
  - その代わり、採用された pixel の信頼性は上がる

重要なのは、`n_min_avg` は kernel の種類そのものを変えるパラメータではなく、**平均を許可する条件**だという点です。したがって、1 点だけの synthetic 入力で全 pixel が NaN になる場合でも、まず `n_min_avg` を確認してください。現在の既定 `n_min_avg=2` では、1 dump だけの入力は kernel が正常でも invalid になります。

### 10.5.2 `n_min_plane`

`estimator='plane'` のときに、局所平面フィットを行うために必要な最小サンプル数です。`avg` より強い局所モデルを使うぶん、通常は `n_min_avg` より厳しい条件が必要になります。

## 10.6 beam / sampling 補助

- `warn_if_cell_coarse`
- `estimate_effective_beam`

## 10.7 diagnostics 出力

- `fill_nan_for_invalid`
- `emit_diag_maps`
- `emit_neff_map`
- `emit_rms_map`
- `emit_time_map`
- `emit_tint_map`
- `emit_tsys_map`

`grid_otf_family(...)` で後段の bundle ext を豊かにしたいなら、これらは有用です。

## 10.8 mask / analysis 系

- `generate_mask`
- `mask_method`
- `mask_sigma`
- `mask_high_snr`
- `mask_low_snr`
- `mask_min_vol`
- `mask_sigma_v`
- `mask_deriv_snr`
- `mask_dilation`
- `mask_compression`

これらは主に後段解析向けで、FFT/PLAIT の必須要素ではありません。

## 10.9 velocity axis

- `dv_kms`
- `vmin_kms`
- `vmax_kms`

X family / Y family の両方で共通にする必要があります。

## 10.10 実行 backend / runtime

- `backend`
- `verbose`
- `workers`
- `sort_neighbors`
- `reproducible_mode`
- `write_diagnostics`
- `diagnostics_prefix`

---

## 11. 実用 cookbook

## 11.1 まずは推奨 kernel 設定で X/Y を別 gridding

### 11.1.1 標準推奨: `sf + convsupport=3`

```python
cfg = MapConfig(
    x0=-3600.0,
    y0=-3600.0,
    nx=201,
    ny=201,
    beam_fwhm_arcsec=350.0,
    cell_arcsec=350.0 / 3.0,   # まずは beam の 1/3 を目安
    kernel='sf',
    convsupport=3,
    weight_mode='rms',
    reproducible_mode=True,
    verbose=True,
)
```

この設定は、

- まず安全に全体像を見たい
- signed kernel による ringing を避けたい
- aliasing や local instability に比較的強い方から始めたい

という場合の標準推奨です。

### 11.1.2 分解能重視: `gjinc + mangum + signed`

```python
cfg = MapConfig(
    x0=-3600.0,
    y0=-3600.0,
    nx=201,
    ny=201,
    beam_fwhm_arcsec=350.0,
    cell_arcsec=350.0 / 3.0,
    kernel='gjinc',
    kernel_preset='mangum',
    kernel_sign='signed',
    weight_mode='rms',
    reproducible_mode=True,
    verbose=True,
)
```

この設定は、

- 細い構造をできるだけ sharp に出したい
- final map 候補として分解能重視の再イメージングをしたい
- その代わり noise amplification や ringing を診断 map と合わせて確認する

という場合に向きます。

### 11.1.3 まず何を見るべきか

この段階でまず見るべきものは、

- `verbose=True` の kernel summary
- nominal effective beam
- empirical center beam
- `WEIGHT`, `WABS`, `CANCEL`, `WREL`
- sparse coverage なら `n_min_avg` によって pixel が落ちていないか

です。特に sparse な test や synthetic 入力では、kernel の違いより先に `n_min_avg` が効いて NaN が増えることがあります。


どちらの設定でも、まずは `grid_otf_family(..., verbose=True)` で出る

- kernel summary
- nominal effective beam
- empirical center beam

を確認してください。また、sharpness 重視設定では、少なくとも

- `WEIGHT`
- `CANCEL`
- family cube そのものの見た目

を合わせて確認するのが安全です。

このあと `grid_otf_family(...)` を X/Y それぞれに使う。

## 11.2 family cube を見てから PLAIT する

推奨順序:

1. X family cube を作る
2. Y family cube を作る
3. まず X/Y を単純に見て stripe の向きを確認する
4. line-free 範囲を決める
5. `basketweave_cubes(...)` または `basketweave_fits(...)`

## 11.3 baseline 後 cube を扱う

- すでに gridding した後で baseline subtraction した cube でも扱ってよい
- ただしその場合は、**coadd や basketweave の直前に、その cube から line-free RMS を再計算**する
- 現在の実装はその方針

## 11.4 既存 FITS と新しい観測分を混ぜる

1. 既存 FITS を `read_otf_bundle(...)`
2. 新しい family を `grid_otf_family(...)`
3. 両者を `coadd_family_cubes(...)`
4. X/Y が揃ったら `basketweave_cubes(...)`

## 11.5 stripe が弱い / clean field のとき

現在の実装と検証では、stripe が弱い場合は PLAIT によって悪化し得ます。

このときは

- family average をそのまま採用する
- `quality_gate_mode` を試す
- まず family cube を見て、stripe が本当にあるか確認する

の順が安全です。

## 11.6 小さい map のとき

現在の実装では fallback されるので、無理に FFT/PLAIT を通そうとしない方がよいです。

目安:

- `min(ny, nx) < 32` なら、まず平均で十分かを確認

---

## 11.7 tile ごとに BW してから mosaic する例

```python
# tile A
xA = grid_otf_family(..., family_label="X", linefree_velocity_windows_kms=["-30:-10", "20:50"])
yA = grid_otf_family(..., family_label="Y", linefree_velocity_windows_kms=["-30:-10", "20:50"])
bwA = basketweave_cubes(
    xA, yA,
    linefree_velocity_windows_kms=["-30:-10", "20:50"],
    plait_noise_mode="family_channel",
)

# tile B も同様に作る

out = mosaic_bundles(
    [bwA, bwB],
    linefree_velocity_windows_kms=["-30:-10", "20:50"],
    gain_min=0.5,
)
```

## 11.8 BW なし tile をそのまま mosaic する例

```python
tileA = grid_otf_family(..., family_label="X", linefree_velocity_windows_kms=["-30:-10", "20:50"])
tileB = grid_otf_family(..., family_label="X", linefree_velocity_windows_kms=["-30:-10", "20:50"])

out = mosaic_bundles(
    [tileA, tileB],
    linefree_velocity_windows_kms=["-30:-10", "20:50"],
    gain_min=0.5,
)
```

この場合 `MOSAIC_GAIN=1` なので、通常の inverse-variance mosaic と同じになります。

## 11.9 line-free 窓を後から変えるとき

```python
out = mosaic_fits(
    ["tileA_bw.fits", "tileB_bw.fits"],
    linefree_velocity_windows_kms=["-40:-20", "30:60"],
    gain_min=0.5,
)
```

このとき現在の実装は、保存済み `MOSAIC_GAIN` / `MOSAIC_TRUST` を再利用しつつ、`MOSAIC_RMS_OBS` / `MOSAIC_WEIGHT` を新しい窓で作り直します。

## 12. 現在の実装の制限

### 12.1 noise model

現在の実装では、noise の扱いは段階ごとに次です。

- OTF gridding: dump-level `BSL_RMS` に基づく inverse-variance weighting
- family coadd: 各スペクトルごとの 2D empirical RMS map
- BW: `family_scalar` または `family_channel`

したがって、厳密な full covariance model を入れているわけではありません。

### 12.2 quality gate

`quality_gate_mode` は実装されていますが、既定では無効です。理由は、単純な line-free RMS gate は、実際に stripe があるケースでも有効な PLAIT を弾くことがあるためです。

### 12.3 FITS/WCS roundtrip

コード上は `read_otf_bundle(...)` / `write_otf_bundle(...)` が実装されています。実運用では、実観測 FITS を用いた full end-to-end の roundtrip を本番環境で確認してください。

### 12.4 mosaic の現在制限

- `mosaic_bundles(...)` / `mosaic_fits(...)` は現在 **同一 shape / 同一 WCS / 同一 unit** 前提です。
- reprojection, union-WCS, partial overlap を含む一般モザイクの自動 resampling はまだ未対応です。
- `MOSAIC_GAIN` は BW に対する 2D 近似であり、完全な channel 依存 gain ではありません。
- `MOSAIC_TRUST` は現在 unity が標準で、位置ずれ等の heuristic downweight は今後の拡張余地です。

### 12.5 source-filling / edge flux

- 一様に広がった信号を Fourier 重みが直接消しているわけではありません。
- しかし `support_taper` / `apodize` により edge flux は弱くなり得ます。
- source-filling map では、既定設定と no-taper/no-apodize の比較確認を推奨します。

### 12.6 synthetic validation で確認したこと

このチャット内では synthetic test により少なくとも次を確認しました。

- `family_scalar` 経路は旧実装と数値一致すること
- `family_channel` は channel 依存 cloud / stripe / noise mismatch で改善が出ること
- broadband に一様な雲では `family_channel` の利点は小さいが、悪化は目立たないこと
- BW edge に対する gain-aware mosaic は bias を避けられること
- `gain_min=0.5` で危険な低 gain edge を適切に除外できること
- line-free 窓を変えたとき、`MOSAIC_RMS_OBS` / `MOSAIC_WEIGHT` を作り直す必要があること

---

## 13. 現在の実装での推奨設定

まずは次が無難です。

### gridding

- `weight_mode='rms'`
- X/Y で完全に同じ `MapConfig`
- 共通 `ref_coord` / `ref_lon`, `ref_lat`
- 同じ `dv_kms`, `vmin_kms`, `vmax_kms`
- mosaic を見据えるなら、tile 作成時から `linefree_velocity_windows_kms` を明示

### basketweave

- channel 依存 noise / 雲 / line 近傍悪化が疑われるなら `plait_noise_mode="family_channel"` をまず検討
- 旧比較や単純場では `plait_noise_mode="family_scalar"` も有効
- `linefree_velocity_windows_kms` を明示
- `science_mask_mode='and'`
- `pad_frac=0.25`
- `apodize=True`
- `support_taper=True`
- `diagnostics=True`
- `small_map_policy='fallback_average'`
- `quality_gate_mode='none'`
- source-filling のときは `support_taper` / `apodize` の影響を必ず点検

### mosaic

- 既定 `gain_min=0.5`
- 後段で line-free 窓を変えたら `MOSAIC_RMS_OBS` / `MOSAIC_WEIGHT` を再計算する
- tile 間の標準重みは `MOSAIC_WEIGHT` をそのまま使い、`1/rms^2` を追加で掛けない
- `WEIGHT`, `WEIGHT_X_FFT`, `HIT` などを mosaic 重みに流用しない

---

## 14. 最後の実務的まとめ

### 14.1 この改訂版で特に重要な追加点

- `plait_noise_mode="family_scalar" / "family_channel"` の 2 経路を区別すること
- BW 出力には入力 family provenance ext と `CHANNEL_NOISE_X/Y` が追加されたこと
- `grid_otf_family`, `coadd_family_cubes`, `basketweave_cubes`, `mosaic_bundles` が `MOSAIC_*` ext を共有するようにしたこと
- `MOSAIC_GAIN` は重みではなく multiplicative response であり、現在の mosaic は `trust * gain^2 / rms_obs^2` を用いること
- `gain_min=0.5` を既定とし、低 gain edge を無理に救わないこと

### 14.2 実務的な最短ルート

最も安全な使い方は、次です。

1. X/Y を別々に gridding する
2. family cube を実際に見る
3. line-free 範囲を決める
4. 必要なら `coadd_family_cubes(...)`
5. 必要なら `basketweave_cubes(...)` か `basketweave_fits(...)`
6. 各 tile の `MOSAIC_*` を確認して `mosaic_bundles(...)` / `mosaic_fits(...)`
7. mosaic 後の `MOSAIC_COMBINE_INFO`, `MOSAIC_WEIGHT_SUM`, `VALID_MASK_UNION` を確認する

### 14.3 一言で言うと

現在の実装を一言で言うと、

- **OTF gridding は従来の inverse-variance gridding を維持**
- **family coadd は各スペクトルごとの empirical RMS map**
- **FFT/PLAIT は `family_scalar` / `family_channel` を切り替え可能**
- **mosaic は gain-aware だが、現在は same-WCS / same-shape 前提**

です。

この順に進めれば、現在の実装の長所を活かしつつ、BW と後段 mosaic を無理なく運用できます。
