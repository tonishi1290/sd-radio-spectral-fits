# 3D電波観測データ解析パッケージ Cookbook

## 1. はじめに

本書は、この 3D 電波観測データ解析パッケージを **実際に使うための実践的なレシピ集** です。  
単なる関数一覧ではなく、

- 何をしたいときに
- どの関数を使い
- どんな入力を与え
- どんな出力が得られ
- 次に何をすればよいか

が分かるように、日本語でできるだけ平易にまとめています。

本 cookbook は、次の 4 つの処理層を行き来しながら使うことを想定しています。

1. **scantable から 3D FITS cube を作る**  
   OTF 観測や PS 観測のデータをグリッディングして FITS に保存する。
2. **cube の baseline を引く**  
   line-free 領域や ripple を推定し、スペクトルの基線を補正する。
3. **3D mask と moment を作る**  
   信号検出用の `MASK3D`、および `MOMENT0` を作る。
4. **baseline 情報を再利用する**  
   前回の line-free や ripple の情報を利用して再解析する。

---

## 2. 最初に押さえるべき考え方

### 2.1 このパッケージでよく出てくる 3 種類の「マスク」

混乱しやすいので、最初に整理します。

#### A. baseline 用の line-free 情報

これは「**線が無いと仮定して baseline fitting に使ったチャネル**」です。

- 代表的な保存先: `LINEFREE`
- 意味: `1 = line-free`, `0 = 線の可能性あり`

#### B. baseline 由来の provisional mask

これは「line-free の補集合」などから作る **暫定的な線候補領域** です。

- 代表的な保存先: `LINECAND3D`, `BASESUP3D`, `MOM0_LINECAND`, `MOM0_BASESUP`
- 意味: baseline 側の情報だけで作った暫定結果

#### C. final signal mask

これは最終的な信号検出結果です。

- 代表的な保存先: `MASK3D`
- 意味: 3D 信号マスク

**重要**  
`MASK3D` は final signal 専用です。  
baseline 情報から作る暫定マスクとは分けて扱います。

---

### 2.2 `RMS` と `BASE_RMS` の違い

#### `RMS`

通常は、マッピング後の cube に対して求めた RMS です。

#### `BASE_RMS`

baseline 補正後の residual RMS です。

baseline 補正後の解析では、一般に **`BASE_RMS` を優先**するのが自然です。  
なぜなら、baseline を引き直した後の最新のノイズ状態を反映しているからです。

---

### 2.3 `BUNIT` と `TEMPSCAL`

外部ツール互換性のため、主 cube の単位は次のように考えるのが安全です。

- `BUNIT = 'K'`
- `TEMPSCAL = 'TA*'` または `TR*`

つまり、`TA*` / `TR*` は `BUNIT` に書き込まず、別 keyword で表します。

---

### 2.4 spectral axis は shape ではなく FITS header で考える

3D cube の軸順序は、

- `(nchan, ny, nx)`
- `(ny, nx, nchan)`
- `(ny, nchan, nx)`

などが混在し得ます。  
このため、**配列 shape だけで決め打ちせず、FITS header の `CTYPE1..3` を見て spectral axis を決める**のが安全です。

---

## 3. 代表的な入口関数

この cookbook では、主に次の高水準関数を使います。

### 3.1 マッピング系

- `run_mapping_pipeline()`  
  OTF を含む一般的な scantable → 3D cube パイプライン
- `run_otf_full_pipeline()`  
  basket-weave 補正付き OTF パイプライン
- `run_ps_mapping_pipeline()`  
  Position Switch 観測用の 3D cube 作成

### 3.2 baseline 系

- `subtract_baseline_from_fits()`  
  既存 FITS cube に対する baseline 引き
- `run_cli_pipeline()`  
  FITS 入出力込みの反復 baseline パイプライン
- `run_one_iteration()`  
  `BaselineSession` を使った 1 反復の制御

### 3.3 解析系

- `make_3d_mask_for_existing_fits()`  
  保存済み FITS に対する `MASK3D` / `MOMENT0` 生成

---

## 4. 共通準備

以下のような import を共通で使うことを想定します。

```python
import numpy as np
from astropy.io import fits

from sd_radio_spectral_fits.map_3d.config import MapConfig
from sd_radio_spectral_fits.map_3d.gridder import run_mapping_pipeline
from sd_radio_spectral_fits.map_3d.otf_pipeline import run_otf_full_pipeline
from sd_radio_spectral_fits.map_3d.ps_gridder import PSMapConfig, run_ps_mapping_pipeline
from sd_radio_spectral_fits.map_3d.baseline_subtraction import (
    LineFreeConfig,
    RippleConfig,
    BaselineConfig,
    subtract_baseline_from_fits,
)
from sd_radio_spectral_fits.map_3d.cube_baseline.cli import run_cli_pipeline
from sd_radio_spectral_fits.map_3d.cube_analysis import make_3d_mask_for_existing_fits
```

実際の import path は、リポジトリ構成やインストール方法に応じて微調整してください。

---

# 5. Cookbook 本編

## Recipe 1. まずは OTF 観測の scantable から 3D cube を作る

### 目的

OTF 観測データを、標準的な 3D FITS cube にしたい。

### 使う関数

- `MapConfig`
- `run_mapping_pipeline()`

### 例

```python
from sd_radio_spectral_fits.map_3d.config import MapConfig
from sd_radio_spectral_fits.map_3d.gridder import run_mapping_pipeline

cfg = MapConfig(
    x0=0.0,
    y0=0.0,
    nx=121,
    ny=121,
    cell_arcsec=7.5,
    beam_fwhm_arcsec=16.0,
    kernel='gjinc',
    chunk_ch=256,
    dtype='float32',
    generate_mask=False,
)

res = run_mapping_pipeline(
    scantable=scantable,
    config=cfg,
    output_fits='otf_map.fits',
    coord_sys='icrs',
    projection='SFL',
    out_scale='TA*',
    dv_kms=0.2,
)
```

### 何が起こるか

この処理は大まかに次を行います。

1. 速度軸をそろえる
2. 天球座標を局所平面へ投影する
3. 温度スケールとビーム効率を整理する
4. グリッディングする
5. 3D FITS と各種診断 HDU を保存する

### 出力として期待するもの

- 主 cube
- `WEIGHT`
- `HIT`
- `RMS`
- `TIME`
- `TINT`
- `TSYS`
- 場合によっては `NEFF`, `XEFF`, `YEFF`, `BIAS_PIX`

### こんなときに向く

- まず cube を作りたい
- OTF でも PS でも、共通の入り口から扱いたい
- 後段の baseline や mask 解析へ渡す前段階として使いたい

---

## Recipe 2. basket-weave 付きで OTF を処理する

### 目的

交差走査のゼロレベル差を補正してから cube を作りたい。

### 使う関数

- `run_otf_full_pipeline()`

### 例

```python
from sd_radio_spectral_fits.map_3d.otf_pipeline import run_otf_full_pipeline

res = run_otf_full_pipeline(
    scantable=scantable,
    config=cfg,
    output_fits='otf_map_bw.fits',
    do_basket_weave=True,
    coord_sys='galactic',
    projection='SFL',
    out_scale='TA*',
    dv_kms=0.2,
)
```

### ポイント

この関数は、

1. basket-weave オフセット推定
2. scantable のデータ実体への補正
3. `run_mapping_pipeline()` への委譲

という流れで動きます。

### 注意

basket-weave は、scan ごとのオフセット差を抑えるのに有効ですが、

- そもそも ON/OFF の取り方が不安定
- 大きな baseline ripple が残っている
- scan の向きや coverage に偏りがある

ときには、万能ではありません。  
cube を見てから baseline や RMS の再確認も併用してください。

---

## Recipe 3. PS 観測から cube を作る

### 目的

Position Switch 観測を、PS 専用のグリッドへ載せたい。

### 使う関数

- `PSMapConfig`
- `run_ps_mapping_pipeline()`

### 例

```python
from sd_radio_spectral_fits.map_3d.ps_gridder import PSMapConfig, run_ps_mapping_pipeline
import numpy as np

x_grid = np.arange(-300, 301, 30)
y_grid = np.arange(-300, 301, 30)

ps_cfg = PSMapConfig(
    coord_sys='icrs',
    projection='SFL',
    ref_lon=83.82208,
    ref_lat=-5.39111,
    x_grid=x_grid,
    y_grid=y_grid,
    grid_tol=10.0,
    combine='mean',
    weight_mode='rms',
    verbose=True,
)

res = run_ps_mapping_pipeline(
    scantable=scantable,
    config=ps_cfg,
    output_fits='ps_map.fits',
    out_scale='TA*',
    dv_kms=0.2,
    vmin_kms=-20.0,
    vmax_kms=40.0,
)
```

### 何が得られるか

- 主 cube
- `HIT`
- `RMS`
- `TSYS`
- `TINT`

OTF よりシンプルな構成になりやすく、等間隔サンプリングに近いデータの整理に向いています。

### 向いている場面

- 規則的な格子点観測
- 特定位置の比較
- 簡潔なマップを早く作りたいとき

---

## Recipe 4. 既存 3D FITS cube から baseline を 1 回だけ引く

### 目的

まずは簡単に、FITS cube に baseline 補正を掛けたい。

### 使う関数

- `LineFreeConfig`
- `RippleConfig`
- `BaselineConfig`
- `subtract_baseline_from_fits()`

### 例

```python
from sd_radio_spectral_fits.map_3d.baseline_subtraction import (
    LineFreeConfig,
    RippleConfig,
    BaselineConfig,
    subtract_baseline_from_fits,
)

lf_cfg = LineFreeConfig(
    smooth_width=7,
    snr_high=4.0,
    snr_low=2.0,
    min_width=3,
)

rp_cfg = RippleConfig(
    nfreq=2,
    period_range_chan=(20, 400),
)

bl_cfg = BaselineConfig(
    poly_order=1,
    ripple=True,
    robust=True,
    chunk_pix=65536,
)

subtract_baseline_from_fits(
    input_fits='otf_map.fits',
    output_fits='otf_map_bsub.fits',
    cube_ext=None,
    linefree_cfg=lf_cfg,
    ripple_cfg=rp_cfg,
    baseline_cfg=bl_cfg,
    add_qc_hdus=True,
    overwrite=True,
)
```

### 何が保存されるか

baseline 後の cube に加えて、典型的には

- `LINEFREE`
- `RIPFREQ`
- `BASE_RMS`
- `BASE_FLG`

が保存されます。

### まず見るべきもの

- 主 cube のスペクトルが平らになっているか
- `LINEFREE` が妥当か
- `BASE_RMS` が極端に悪い場所がないか
- `BASE_FLG` に失敗画素が集中していないか

---

## Recipe 5. 反復 baseline fitting を行う

### 目的

1 回の baseline では不十分なので、line-free と ripple を更新しながら何回か回したい。

### 使う関数

- `run_cli_pipeline()`

### 例

```python
from sd_radio_spectral_fits.map_3d.cube_baseline.cli import run_cli_pipeline
from sd_radio_spectral_fits.map_3d.baseline_subtraction import LineFreeConfig, RippleConfig

run_cli_pipeline(
    input_fits='otf_map.fits',
    output_fits='otf_map_iterbsub.fits',
    cube_ext=None,
    iterations=3,
    poly_order=1,
    auto_linefree=True,
    linefree_cfg=LineFreeConfig(),
    manual_v_windows=['-5:15'],
    enable_ripple=True,
    ripple_cfg=RippleConfig(nfreq=2, period_range_chan=(20, 400)),
    robust=True,
    chunk_pix=65536,
    run_cube_analysis=False,
)
```

### 何が便利か

この高水準関数は、

- FITS を読む
- velocity axis を作る
- `BaselineSession` を初期化する
- `run_one_iteration()` を繰り返す
- baseline 後の FITS を書く

までを一括でやってくれます。

### 使いどころ

- とりあえず素早く反復 baseline をかけたい
- 後で `cube_analysis` は別途行いたい
- CLI 的にも Python API 的にも使いたい

---

## Recipe 6. baseline 時に「ここは線がありそう」と手で教える

### 目的

自動推定だけでは危ないので、手動で線領域を指定したい。

### 使う関数

- `run_cli_pipeline()` の `manual_v_windows`
- `create_manual_signal_mask_1d()`

### 例

```python
run_cli_pipeline(
    input_fits='otf_map.fits',
    output_fits='otf_map_manualwin.fits',
    iterations=2,
    auto_linefree=True,
    manual_v_windows=['-10:20', '45:60'],
    enable_ripple=True,
)
```

### どう解釈されるか

- `manual_v_windows` で与えた速度範囲は **signal 候補**として扱われます
- line-free mask は、それらを除いた領域になります

### こんなときに重要

- 非常に強い線がある
- 自動 line-free 推定が裾を見落としやすい
- 弱い wing を守りたい

---

## Recipe 7. baseline 結果を見ながら、さらに手で refinement したい

### 目的

自動推定 → 目視確認 → 再 baseline の流れを回したい。

### 基本的な考え方

1. まず baseline を 1 回実行
2. `LINEFREE`, `BASE_RMS`, スペクトルを確認
3. 必要なら manual velocity windows を追加
4. もう一度 `run_cli_pipeline()` または `subtract_baseline_from_fits()` を実行

### 例

```python
# 1回目
run_cli_pipeline(
    input_fits='otf_map.fits',
    output_fits='tmp_iter1.fits',
    iterations=1,
    auto_linefree=True,
    enable_ripple=True,
)

# LINEFREE やスペクトルを確認後、2回目
run_cli_pipeline(
    input_fits='otf_map.fits',
    output_fits='tmp_iter2.fits',
    iterations=2,
    auto_linefree=True,
    manual_v_windows=['-8:22'],
    enable_ripple=True,
)
```

### 運用のコツ

本当に重要なのは、

- 1 回で終わらせようとしない
- line-free と線候補を明確に分ける
- final signal mask と baseline 由来 provisional mask を混同しない

ことです。

---

## Recipe 8. baseline 後に `MASK3D` と `MOMENT0` を作る

### 目的

baseline 後の cube から final signal mask と moment0 を作りたい。

### 使う関数

- `make_3d_mask_for_existing_fits()`

### 例

```python
from sd_radio_spectral_fits.map_3d.cube_analysis import make_3d_mask_for_existing_fits

make_3d_mask_for_existing_fits(
    input_fits='otf_map_bsub.fits',
    output_fits='otf_map_bsub_masked.fits',
    cube_ext=None,
    rms_ext='BASE_RMS',
    method='smooth_mask_lite',
    convert_to_kms=True,
    require_kms=False,
    rms_spectral_slab=(-100.0, 100.0),
    moment_spectral_slab=(-20.0, 40.0),
)
```

### 重要ポイント

- baseline 後なら `rms_ext='BASE_RMS'` を使うのが自然
- `MASK3D` は final signal mask
- `MOMENT0` は `MASK3D` を掛けた後の積分結果

### 典型的な出力

- `MASK3D`
- `MOMENT0`

場合によっては、baseline 由来の provisional 系 HDU も併せて扱います。

---

## Recipe 9. provisional mask も一緒に使いたい

### 目的

baseline 情報から得た暫定的な線候補領域も、解析に活かしたい。

### 基本方針

次のように分けて考えると混乱しにくいです。

- `LINEFREE`  
  baseline に使った「線が無い」と考えた領域
- `LINECAND3D`  
  baseline 的には「線の可能性がある」暫定領域
- `BASESUP3D`  
  baseline 側の支持情報から作る補助マスク
- `MASK3D`  
  final signal mask

### 推奨運用

#### baseline の成否を見るとき

- `LINEFREE`
- `BASE_RMS`
- `BASE_FLG`
- `LINECAND3D`

をまず見る

#### 科学解析で moment を作るとき

- 最終的には `MASK3D` を使う
- `MOM0_LINECAND` は provisional として参考にする

### 重要

深積分すると、初期には baseline 側で「線が無い」と思っていたところに、弱い広がった信号が見えることがあります。  
そのため、**provisional mask と final signal mask を同一視しない**ことが大切です。

---

## Recipe 10. baseline 後にすぐ解析まで一気にやりたい

### 目的

1 本の高水準関数で、baseline → mask/moment まで一気に進みたい。

### 例

```python
run_cli_pipeline(
    input_fits='otf_map.fits',
    output_fits='otf_map_allinone.fits',
    iterations=2,
    poly_order=1,
    auto_linefree=True,
    enable_ripple=True,
    robust=True,
    run_cube_analysis=True,
    cube_analysis_method='smooth_mask_lite',
)
```

### 向いているケース

- まず全体像を見たい
- クイックルックが欲しい
- 反復 baseline 後の `MASK3D` と `MOMENT0` まで一度に欲しい

### 注意

最終的な論文品質の解析では、後から `rms_ext` や slab 範囲を明示して再解析する方が安心です。

---

## Recipe 11. `simple`, `smooth_mask_lite`, `smooth_mask`, `derivative` を使い分ける

### 目的

3D mask の作り方を、対象に応じて選びたい。

### 1. `simple`

#### 向く場面

- 強い線
- まずざっくり見たい
- 高 S/N データ

#### 例

```python
make_3d_mask_for_existing_fits(
    'cube_bsub.fits',
    output_fits='cube_simplemask.fits',
    method='simple',
    rms_ext='BASE_RMS',
    mask_sigma=3.0,
)
```

### 2. `smooth_mask_lite`

#### 向く場面

- 初期の標準選択肢
- 弱い信号も少し拾いたい
- ただし過剰に複雑にはしたくない

#### 例

```python
make_3d_mask_for_existing_fits(
    'cube_bsub.fits',
    output_fits='cube_smoothelite.fits',
    method='smooth_mask_lite',
    rms_ext='BASE_RMS',
    mask_high_snr=3.0,
    mask_low_snr=1.5,
    mask_min_vol=27,
)
```

### 3. `smooth_mask`

#### 向く場面

- 空間的にも速度的にも広がる信号
- 弱い構造を丁寧に取りたい

#### 注意

ノイズのつながりを拾いすぎる場合があるので、`high_snr`, `low_snr`, `min_vol` の調整が重要です。

### 4. `derivative`

#### 向く場面

- エッジや速度構造がはっきりしている
- 単純なしきい値では拾いにくい

#### 例

```python
make_3d_mask_for_existing_fits(
    'cube_bsub.fits',
    output_fits='cube_derivmask.fits',
    method='derivative',
    rms_ext='BASE_RMS',
    mask_sigma_v=2.0,
    mask_deriv_snr=3.0,
    mask_dilation=2,
)
```

---

## Recipe 12. RMS 計算と moment の積分範囲を分ける

### 目的

RMS は広い無信号域で見積もりたいが、moment は狭い速度範囲だけ積分したい。

### 例

```python
make_3d_mask_for_existing_fits(
    input_fits='cube_bsub.fits',
    output_fits='cube_bsub_mom.fits',
    method='smooth_mask_lite',
    rms_ext='BASE_RMS',
    rms_spectral_slab=(-100.0, 100.0),
    moment_spectral_slab=(-5.0, 20.0),
)
```

### どう使い分けるか

- `rms_spectral_slab`  
  ノイズ評価のために使う範囲
- `moment_spectral_slab`  
  実際に積分したい速度範囲

### 例の意味

- RMS は広い範囲で安定に見る
- Moment0 は目的の線だけに絞って作る

これは日常的によく使う運用です。

---

## Recipe 13. baseline の prior を別の解析に活かす

### 目的

同じマッピングを積み増ししたり、同一対象の別 cube に近い初期値を与えたい。

### 基本ルール

#### 再利用してよいもの

- `LINEFREE`
- `RIPFREQ`

これらは **prior** として使えます。

#### 再利用せず再計算すべきもの

- `BASE_RMS`
- `BASE_FLG`

これらは、その cube 自身のフィット結果だからです。

### 推奨の考え方

- 前回の `LINEFREE` を初期値にする
- 前回の `RIPFREQ` を初期値にする
- しかし最終的な baseline は今回の cube に対してもう一度求める

### 典型例

- 1 回目の浅い積分から provisional baseline を作る
- 深積分後に、その prior を使って baseline を引き直す
- 最終的な `MASK3D` と `MOMENT0` は深積分版で再作成する

---

## Recipe 14. `BaselineSession` を使って細かく制御する

### 目的

高水準 CLI ではなく、自分で 1 反復ずつ制御したい。

### 使う関数

- `BaselineSession`
- `run_one_iteration()`

### 例

```python
import numpy as np
from astropy.io import fits
from sd_radio_spectral_fits.map_3d.cube_baseline.session import BaselineSession
from sd_radio_spectral_fits.map_3d.cube_baseline.orchestrator import run_one_iteration
from sd_radio_spectral_fits.map_3d.cube_baseline.cli import build_v_axis_kms_from_header

with fits.open('cube.fits') as hdul:
    cube = hdul[0].data
    header = hdul[0].header

v_axis = build_v_axis_kms_from_header(header, cube.shape[0])
session = BaselineSession(cube, v_axis)

# 1回目
run_one_iteration(
    session,
    auto_linefree=True,
    manual_v_windows=['-5:15'],
    enable_ripple=True,
    poly_order=1,
    robust=True,
)

# target_mask_2d を使って局所的にやり直す例
ny, nx = session.ny, session.nx
target = np.zeros((ny, nx), dtype=bool)
target[30:60, 40:80] = True

run_one_iteration(
    session,
    target_mask_2d=target,
    auto_linefree=True,
    enable_ripple=True,
    poly_order=1,
    robust=True,
)

cube_corrected = session.get_full_cube_work()
```

### 向いている場面

- 一部の領域だけ再 baseline したい
- 自分で毎回 intermediate state を見たい
- GUI や notebook と組み合わせたい

---

## Recipe 15. baseline だけ先に見て、moment を provisional に作る

### 目的

final signal mask を作る前に、まず baseline 情報だけから暫定的に積分したい。

### 推奨方針

- provisional な `MOM0_LINECAND`
- provisional な `MOM0_BASESUP`

を作り、初期評価に使う。

### 使いどころ

- 積分前に線候補位置をざっと見たい
- signal mask 生成前の QC をしたい
- baseline の失敗・過剰除去をチェックしたい

### 注意

これは **最終結果ではありません**。  
論文や最終定量には `MASK3D` ベースの `MOMENT0` を使うのが原則です。

---

## Recipe 16. baseline 済み cube をもう一度 baseline し直す

### 目的

前回の baseline が不十分だったので、再 baseline したい。

### まず確認すること

- その FITS に `LINEFREE` が入っているか
- `RIPFREQ` が入っているか
- 古い `BASE_RMS` / `BASE_FLG` が残っていないか

### 推奨手順

1. 既存の baseline provenance を確認
2. prior として再利用するものを決める
3. baseline は今回の cube に対して再計算
4. 古い解析 HDU は混ぜない
5. 必要なら `MASK3D` / `MOMENT0` も作り直す

### 実務的な判断

- `LINEFREE` は prior として有用
- `BASE_RMS` は引き継がない
- final `MASK3D` も、通常は作り直す

---

## Recipe 17. PS と OTF を同じ流れで比較する

### 目的

同一領域の PS 観測と OTF 観測を、できるだけ同じ後処理で比較したい。

### 推奨フロー

#### OTF 側

```python
run_otf_full_pipeline(...)
subtract_baseline_from_fits(...)
make_3d_mask_for_existing_fits(...)
```

#### PS 側

```python
run_ps_mapping_pipeline(...)
subtract_baseline_from_fits(...)
make_3d_mask_for_existing_fits(...)
```

### 比較時の注意

- `out_scale` をそろえる
- velocity grid をそろえる
- mask method をそろえる
- moment の積分範囲をそろえる
- `BASE_RMS` / `RMS` のどちらを使ったかを明記する

---

## Recipe 18. `RMS` と `BASE_RMS` のどちらを使うか迷ったら

### 原則

#### baseline 前・通常のマップ品質評価

- `RMS`

#### baseline 後の signal mask / moment 解析

- `BASE_RMS`

### なぜか

baseline を引き直した後は、ノイズの見え方も変わるからです。  
古い `RMS` をそのまま使うと、現状の residual ノイズを反映しない場合があります。

### 実務メモ

迷ったら、

- baseline 前の比較には `RMS`
- baseline 後の解析には `BASE_RMS`

で始めるのが無難です。

---

## Recipe 19. FITS の中身を確認してから解析する

### 目的

何が入っているか分からない FITS に対して、安全に次の処理を決めたい。

### 例

```python
from astropy.io import fits

with fits.open('cube_bsub.fits') as hdul:
    hdul.info()
    for h in hdul:
        print(h.name, h.data.shape if h.data is not None else None)
```

### 最低限見たいもの

- 主 cube はどの HDU か
- `LINEFREE` があるか
- `RIPFREQ` があるか
- `BASE_RMS` があるか
- `MASK3D` / `MOMENT0` があるか
- 古い解析結果が残っていないか

### 重要

再解析するときは、**入力 cube そのもの**と **解析生成物 HDU** を取り違えないことが非常に重要です。

---

## Recipe 20. 何から始めればよいか分からない人向けの最短コース

### まずはこの 3 段階で十分

#### Step 1. cube を作る

```python
run_mapping_pipeline(...)
```

または

```python
run_ps_mapping_pipeline(...)
```

#### Step 2. baseline を引く

```python
subtract_baseline_from_fits(...)
```

#### Step 3. mask と moment を作る

```python
make_3d_mask_for_existing_fits(...)
```

### 最初のおすすめ設定

- baseline: `poly_order=1`, `ripple=True`, `robust=True`
- mask: `method='smooth_mask_lite'`
- RMS: `BASE_RMS`
- moment slab: 目的線に近い範囲

### 最後に確認するもの

- representative spectrum
- `BASE_RMS`
- `MASK3D`
- `MOMENT0`

---

# 6. 代表的な解析手段の整理表

| やりたいこと | 使う主関数 | 主な出力 |
|---|---|---|
| OTF から cube を作る | `run_mapping_pipeline` | cube, RMS, HIT, WEIGHT |
| basket-weave 付き OTF | `run_otf_full_pipeline` | 補正済み cube |
| PS 観測から cube | `run_ps_mapping_pipeline` | cube, RMS, HIT, TINT |
| 既存 cube の baseline 引き | `subtract_baseline_from_fits` | baseline 後 cube, LINEFREE, BASE_RMS |
| 反復 baseline | `run_cli_pipeline` | 反復 baseline 後 cube |
| 3D signal mask | `make_3d_mask_for_existing_fits` | MASK3D |
| Moment0 | `make_3d_mask_for_existing_fits` | MOMENT0 |
| provisional line candidate | baseline provenance + 解析 | LINECAND3D, MOM0_LINECAND |
| prior 再利用 | `LINEFREE`, `RIPFREQ` を読み再 baseline | 改良版 baseline |

---

# 7. よくある失敗と対処

## 7.1 `MASK3D` と `LINEFREE` を混同する

### 失敗例

- `LINEFREE` を final signal mask だと思う
- `MASK3D` を baseline 用の窓に再利用する

### 対処

- `LINEFREE` は baseline 用
- `MASK3D` は final signal 用
- provisional 系は別に扱う

---

## 7.2 古い `RMS` を使ってしまう

### 失敗例

baseline を引き直したのに、昔の `RMS` をそのまま使う。

### 対処

baseline 後の解析では `BASE_RMS` を優先する。

---

## 7.3 deep integration で弱い線が出てきて混乱する

### 失敗例

浅い積分では line-free だった領域に、深積分で信号が出てくる。

### 対処

- provisional と final を分ける
- `LINEFREE` を絶対視しない
- 深積分後に baseline も `MASK3D` も再検討する

---

## 7.4 `manual_v_windows` を入れずに強線で baseline が壊れる

### 対処

明らかに強い線や broad wing がある場合は、早い段階から `manual_v_windows` を使う。

---

## 7.5 解析生成物 HDU を入力 cube と取り違える

### 対処

FITS を開いて、どれが主 cube か先に確認する。  
`MASK3D`, `BASESUP3D`, `LINECAND3D`, `MOMENT0`, `RMS`, `BASE_RMS` は解析生成物です。

---

# 8. 初学者向けのおすすめワークフロー

## パターン A. 一番基本

```python
run_mapping_pipeline(...)
subtract_baseline_from_fits(...)
make_3d_mask_for_existing_fits(...)
```

### こんな人向け

- まずは一通り動かしたい
- OTF / PS を問わず基本を押さえたい

---

## パターン B. 慎重派

```python
run_mapping_pipeline(...)
subtract_baseline_from_fits(...)
# LINEFREE, BASE_RMS, representative spectra を確認
run_cli_pipeline(..., iterations=2, manual_v_windows=[...])
make_3d_mask_for_existing_fits(...)
```

### こんな人向け

- 弱線や broad wing が重要
- baseline の信頼性を重視したい

---

## パターン C. 積分を重ねながら洗練する

```python
# 浅い積分で provisional baseline
run_cli_pipeline(...)
# prior を参考に深積分版を再 baseline
run_cli_pipeline(...)
# 最終 mask / moment
make_3d_mask_for_existing_fits(...)
```

### こんな人向け

- 同一領域を繰り返し観測する
- 前回結果を次に活かしたい

---

# 9. まとめ

このパッケージは、

- scantable → cube
- cube → baseline 補正
- baseline 済み cube → signal mask / moment
- baseline prior の再利用

を一貫して扱えるように設計されています。

特に大切なのは、次の 5 点です。

1. `MASK3D` は final signal 専用
2. baseline 由来情報は provisional 系として分離
3. baseline 後の解析では `BASE_RMS` を優先
4. `BUNIT='K'` と `TEMPSCAL` を分ける
5. spectral axis は header ベースで解釈する

最初は、

- cube を作る
- baseline を引く
- `MASK3D` と `MOMENT0` を作る

の 3 段階から始めれば十分です。  
その後、必要に応じて

- manual velocity windows
- prior 再利用
- provisional mask
- method の切り替え

を加えていくのがおすすめです。

---

# 10. 付録: 最小サンプル集

## A. OTF → baseline → moment

```python
cfg = MapConfig(
    x0=0.0, y0=0.0,
    nx=121, ny=121,
    cell_arcsec=7.5,
    beam_fwhm_arcsec=16.0,
)

run_mapping_pipeline(scantable, cfg, 'map.fits', dv_kms=0.2)

subtract_baseline_from_fits(
    'map.fits', 'map_bsub.fits',
    baseline_cfg=BaselineConfig(poly_order=1, ripple=True, robust=True),
)

make_3d_mask_for_existing_fits(
    'map_bsub.fits',
    output_fits='map_bsub_masked.fits',
    rms_ext='BASE_RMS',
    method='smooth_mask_lite',
    moment_spectral_slab=(-5.0, 20.0),
)
```

## B. PS → baseline → quicklook

```python
ps_cfg = PSMapConfig(
    coord_sys='icrs',
    projection='SFL',
    ref_lon=83.82208,
    ref_lat=-5.39111,
    x_grid=np.arange(-300, 301, 30),
    y_grid=np.arange(-300, 301, 30),
)

run_ps_mapping_pipeline(scantable, ps_cfg, 'ps_map.fits', dv_kms=0.2)

run_cli_pipeline(
    input_fits='ps_map.fits',
    output_fits='ps_map_bsub.fits',
    iterations=2,
    manual_v_windows=['-5:15'],
    run_cube_analysis=True,
    cube_analysis_method='smooth_mask_lite',
)
```

## C. まず provisional に線候補を見る

```python
subtract_baseline_from_fits('cube.fits', 'cube_bsub.fits')

# その後、LINEFREE の補集合や provisional HDU を確認して
# 線候補がどこにあるかをざっと把握する
# final な科学解析は、MASK3D / MOMENT0 を作成してから行う
```



---

## 16. 実践例: Orion KL の OTF `^{12}CO` データを 1 本の流れで解析する

この節では、実際に使われている OTF 解析スクリプトを、**初めて見る人でも流れが分かるように**整理して紹介します。

### 16.1 何をするレシピか

このレシピでは、次の 5 段階を一気に行います。

1. 生の OTF scantable を `Ta*` に較正する
2. basket-weave で scan ごとの縞状オフセットを補正する
3. OTF グリッディングで 3D FITS cube を作る
4. cube に baseline fitting をかける
5. 必要なら `MASK3D` と `MOMENT0` を作る

「OTF の raw から、見やすい baseline 済み cube まで持っていきたい」ときの標準レシピです。

### 16.2 まず全体コードを見る

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
from sd_radio_spectral_fits.map_3d.cube_analysis import make_3d_mask_for_existing_fits

# ============================================================
# Step 1. Chopper wheel calibration
# ============================================================
sc_cal_12co = cal.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-30.0, 55.0),
    t_hot_k=300.0,
    vcorr_chunk_sec=10.0,
    dtype="float32",
)

# ============================================================
# Step 2. Basket-weave correction
# ============================================================
# OTF の scan ごとのオフセットを求めて補正する

grid_input = create_grid_input(sc_cal_12co, projection="CAR")
offsets = solve_basket_weave_offsets(
    grid_input,
    search_radius_arcsec=35.0,
)
apply_basket_weave_correction(grid_input, offsets)

# 補正結果を scantable に戻す
sc_cal_12co.data = grid_input.spec

# ============================================================
# Step 3. OTF cube を作る
# ============================================================
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

# ============================================================
# Step 4. baseline fitting
# ============================================================
with fits.open(cube_fits) as hdul:
    nchan = hdul[0].data.shape[0]

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

# ============================================================
# Step 5. 3D signal mask と moment0
# ============================================================
make_3d_mask_for_existing_fits(
    input_fits=cube_bsub_fits,
    output_fits="ori-kl-12co_bsub_masked.fits",
    rms_ext="auto",
    method="smooth_mask_lite",
    moment_spectral_slab=(-20.0, 40.0),
    overwrite=True,
)
```

### 16.3 各 step をやさしく説明する

#### Step 1. `run_tastar_calibration()`

ここでは、生の観測データを `Ta*` のスペクトルに直します。

- `input_data=sc_raw`
  - 生の scantable です。
- `vlsrk_range_kms=(-30, 55)`
  - 速度補正を考える範囲です。
- `t_hot_k=300`
  - HOT load 温度です。
- `vcorr_chunk_sec=10`
  - 速度補正を何秒ごとにまとめるかです。
- `dtype='float32'`
  - 大きいデータを扱いやすくするためです。

この結果 `sc_cal_12co` は、「Ta* に較正された OTF scantable」になります。

#### Step 2. `create_grid_input()` と basket-weave 補正

OTF の縞模様は、scan ごとの基線オフセットが原因で出ることがあります。そこでまず、scan 同士の交差点を利用してオフセットを見積もります。

- `create_grid_input(sc_cal_12co, projection='CAR')`
  - scantable を、マッピングエンジン向けの `GridInput` に変換します。
- `solve_basket_weave_offsets(..., search_radius_arcsec=35.0)`
  - 交差点とみなす半径を 35 arcsec に設定しています。
- `apply_basket_weave_correction()`
  - 推定したオフセットを各スペクトルから差し引きます。

最後に `sc_cal_12co.data = grid_input.spec` として、補正済みスペクトルを scantable 本体へ戻します。

#### Step 3. `GridConfig` を作って cube 化

ここでは「どの大きさの地図を、どの pixel サイズで作るか」を決めます。

- `x0=-3600`, `y0=-3600`
  - 地図の左下に相当する基準オフセットです。
- `nx=201`, `ny=201`
  - 201 × 201 pixel の地図にします。
- `cell_arcsec=36`
  - 1 pixel の大きさです。
- `beam_fwhm_arcsec=350`
  - 望遠鏡ビームです。
- `kernel='gjinc'`
  - OTF でよく使われる gridding kernel です。
- `truncate='first_null'`
  - GJINC kernel の打ち切り方です。

`run_mapping_pipeline()` は、

- 速度軸の再サンプリング
- 座標投影
- 温度スケール処理
- グリッディング
- FITS 保存

をまとめて行います。

#### Step 4. baseline fitting

ここが、この例の一番大事な段階です。

まず `nchan` を読み、長さ `nchan` の boolean 配列 `lf_mask` を作ります。

- `True` は「baseline fitting に使ってよい」
- `False` は「本物の線がありそうなので使わない」

という意味です。

この例では

```python
lf_mask[169:215] = False
```

として、線が入っているチャネルだけを外しています。

次に `RippleConfig` と `BaselineConfig` を与えます。

- `RippleConfig(nfreq=6, period_range_chan=(10, 1000))`
  - ripple の候補を 6 個まで探す
  - 周期は 10〜1000 channel を見る
- `BaselineConfig(poly_order=1, ripple=True, robust=True)`
  - 1 次 baseline
  - ripple も含める
  - 外れ値に強い fitting を使う

この結果、`ori-kl-12co_bsub.fits` には

- baseline 補正後の cube
- `LINEFREE`
- `RIPFREQ`
- `BASE_RMS`
- `BASE_FLG`

が保存されます。

#### Step 5. `MASK3D` と `MOMENT0`

最後に、baseline 後の cube から 3D の signal mask を作ります。

- `rms_ext='auto'`
  - `BASE_RMS` があればそれを優先
  - なければ `RMS`
  - さらに無ければ cube から推定
- `method='smooth_mask_lite'`
  - 軽量で使いやすい signal mask 生成法です。
- `moment_spectral_slab=(-20, 40)`
  - その速度範囲だけで moment0 を作ります。

### 16.4 このレシピを自分のデータに移すときの考え方

#### 16.4.1 `vlsrk_range_kms`

まず、線が十分に含まれ、かつ前後に baseline 領域が残るようにします。線の幅だけでなく、余白も必要です。

#### 16.4.2 `search_radius_arcsec`

scan 間隔や pointing 精度に合わせます。

- 小さすぎる → 交差点が見つからない
- 大きすぎる → 無関係な点を結びやすい

最初は beam の 0.1〜0.2 倍程度から試すと分かりやすいです。

#### 16.4.3 `cell_arcsec`

一般には beam の 1/3〜1/4 程度がよく使われます。今回の例は 350 arcsec に対して 36 arcsec なので、かなり細かめです。観測のサンプリング密度に対して無理がないかは `HIT` map で確認します。

#### 16.4.4 `lf_mask`

最初の 1 回は、**自動推定より手で決めた方が安全**なことがよくあります。特に

- 強い線
- 広い翼
- 立派な ripple
- line forest

があるときは、まず手動 line-free から入る方が安定します。

### 16.5 この例のあとに見るべきもの

このレシピのあとに、最低限次を確認します。

1. `HIT`
   - 地図のどこに十分データが入っているか
2. `RMS` または `BASE_RMS`
   - 雑音がどこで悪いか
3. `BASE_FLG`
   - baseline fitting に失敗した場所がないか
4. `MASK3D`
   - 信号マスクが広すぎないか、狭すぎないか
5. `MOMENT0`
   - 期待した構造が出ているか

### 16.6 よくあるつまずき

#### 16.6.1 basket-weave をかけたのに stripe が残る

原因候補:

- `search_radius_arcsec` が不適切
- scan ID の扱いが観測実態と合っていない
- stripe の主因が scan オフセットではなく、スペクトル baseline にある

#### 16.6.2 baseline 後に線まで削れてしまう

原因候補:

- `lf_mask` で線を十分に除外できていない
- `poly_order` が高すぎる
- ripple モデルが強すぎる

#### 16.6.3 `MOMENT0` が汚い

原因候補:

- baseline がまだ残っている
- `moment_spectral_slab` が広すぎる
- `MASK3D` が緩すぎる

### 16.7 このレシピの要点を一言で言うと

- **OTF の縞は basket-weave で抑える**
- **チャネル方向の基線は cube baseline で抑える**
- **最初の 1 回は line-free を手で決めると安定しやすい**
- **最後に `MASK3D` と `MOMENT0` を作って、見た目と QC を確認する**

この流れを覚えておけば、OTF の raw から「見られる 3D cube」までをかなり整理して進められます。


---

## Recipe 21. baseline 済み PS scantable `sc_all_coadded_12co` から cube と mask を作る

### 目的

Position Switching 観測で、すでに baseline fitting 済みの scantable `sc_all_coadded_12co` から、

1. まず 3D FITS cube を作る
2. 観測品質を `HIT` / `RMS` / `TSYS` / `TINT` で確認する
3. 必要なら `MASK3D` と `MOMENT0` を追加する

という一連の流れを、なるべく迷わず実行する。

### 使う関数

- `PSMapConfig`
- `run_ps_mapping_pipeline()`
- `make_3d_mask_for_existing_fits()`

### このレシピが向いている状況

- baseline fitting は scantable の段階ですでに終わっている
- `sc_all_coadded_12co` のように、複数点の PS 観測を 1 つの scantable にまとめている
- まずは素直に cube を作り、signal mask はその後で考えたい
- baseline をもう一度 cube に掛けたいわけではない

### まず押さえること

このケースでは、**最初に使うべき入口は `run_ps_mapping_pipeline()`** です。

理由は単純で、`subtract_baseline_from_fits()` は **3D FITS cube に対して baseline を掛ける関数**だからです。
今回の `sc_all_coadded_12co` はすでに baseline 済みの scantable なので、最初の仕事は baseline ではなく **cube 化** です。

### 実行例

```python
import numpy as np
from sd_radio_spectral_fits.map_3d.ps_gridder import PSMapConfig, run_ps_mapping_pipeline
from sd_radio_spectral_fits.map_3d.cube_analysis import make_3d_mask_for_existing_fits

# 例: Orion KL 周辺を 36 arcsec 格子で並べる場合
# x_grid, y_grid は「基準位置からのオフセット [arcsec]」
# 観測点に合わせて適宜変えてください。
x_grid = np.arange(-360.0, 360.0 + 36.0, 36.0)
y_grid = np.arange(-360.0, 360.0 + 36.0, 36.0)

ps_cfg = PSMapConfig(
    coord_sys='icrs',
    projection='SFL',
    ref_lon=83.809,          # 例: 基準 RA [deg]
    ref_lat=-5.372639,       # 例: 基準 Dec [deg]
    x_grid=x_grid,
    y_grid=y_grid,
    grid_tol=12.0,           # 格子点から何 arcsec ずれていても同一セルに入れるか
    combine='mean',
    weight_mode='rms',       # scantable の BSL_RMS / RMS を使って重み付け
    invert_x=True,
    verbose=True,
)

res = run_ps_mapping_pipeline(
    scantable=sc_all_coadded_12co,
    config=ps_cfg,
    output_fits='sc_all_coadded_12co_cube.fits',
    out_scale='TA*',
    dv_kms=0.2,
    vmin_kms=-30.0,
    vmax_kms=55.0,
)

# 必要なら signal mask と moment を追加
make_3d_mask_for_existing_fits(
    input_fits='sc_all_coadded_12co_cube.fits',
    output_fits='sc_all_coadded_12co_cube_masked.fits',
    cube_ext=None,
    rms_ext='RMS',
    method='smooth_mask_lite',
    moment_spectral_slab=(-10.0, 30.0),
    overwrite=True,
)
```

### 何をしているのか

#### 1. `PSMapConfig` を作る

ここでは「どの格子へ並べるか」を決めています。

- `coord_sys='icrs'`
  - 座標を RA/Dec 系として扱います。
  - すでに ICRS/J2000 系で解析したいときの標準的な指定です。
- `projection='SFL'`
  - 天球上の位置を局所平面へ写す投影法です。
  - 小領域マップではまずこれで十分です。
- `ref_lon`, `ref_lat`
  - マップの基準中心です。
  - `x_grid`, `y_grid` はこの基準からのオフセットになります。
- `x_grid`, `y_grid`
  - 格子点そのものです。
  - ここでは 36 arcsec 間隔の等間隔格子を明示しています。
- `grid_tol`
  - 観測点が格子点からどの程度ずれていても、そのセルへ入れてよいかの許容値です。
  - データ取得時のずれがあるなら、0 ではなく少し余裕を持たせます。
- `combine='mean'`
  - 同じセルに複数スペクトルが入ったとき、平均でまとめます。
- `weight_mode='rms'`
  - baseline 後の品質を使って重み付けします。
  - `BSL_RMS` があればそれが使われやすく、baseline 済み scantable にはまずこれが自然です。

#### 2. `run_ps_mapping_pipeline()` を呼ぶ

この 1 回で、次がまとめて行われます。

- 速度軸の標準化
- 必要なら `vmin_kms` / `vmax_kms` による速度範囲の制限
- 座標投影
- 温度スケールの整理
- PS グリッディング
- 3D FITS と診断 HDU の保存

#### 3. `make_3d_mask_for_existing_fits()` を呼ぶ

cube ができたあとで、signal mask と moment0 を追加します。

- `rms_ext='RMS'`
  - `run_ps_mapping_pipeline()` が出した `RMS` HDU を使います。
- `method='smooth_mask_lite'`
  - まずは軽量で使いやすい方法です。
- `moment_spectral_slab=(-10, 30)`
  - moment0 をその速度範囲でだけ積分します。

### 出力される FITS の見方

#### `sc_all_coadded_12co_cube.fits`

- Primary HDU
  - 3D cube 本体
- `HIT`
  - その画素に何本のスペクトルが入ったか
- `RMS`
  - 画素ごとの雑音の目安
- `TSYS`
  - 系雑音温度の代表値
- `TINT`
  - 実効積分時間

#### `sc_all_coadded_12co_cube_masked.fits`

上に加えて、通常は次も入ります。

- `MASK3D`
  - final signal mask
- `MOMENT0`
  - 積分強度マップ
- 条件によっては provisional 系 HDU
  - `BASESUP3D`
  - `LINECAND3D`
  - `MOM0_BASESUP`
  - `MOM0_LINECAND`

### 最初にどこを見るべきか

1. `HIT`
   - 被覆が極端に薄い場所がないか
2. `RMS`
   - 雑音が大きく悪化している領域がないか
3. `TSYS`
   - 観測そのものの品質が大きく変動していないか
4. cube の代表スペクトル
   - line が期待した速度範囲に出ているか
5. `MOMENT0`
   - mask が厳しすぎないか、甘すぎないか

### よくある判断

#### baseline 済みなのに、もう一度 `subtract_baseline_from_fits()` を掛けるべきか

まずは **掛けない** のが基本です。

理由は、scantable 段階ですでに baseline fitting しているなら、二重に基線補正を掛けることになるからです。

ただし、cube 化したあとに

- 大域的にわずかな傾きが残る
- basket-weave や gridding 後に広い ripple が見える
- spatial average で見たときに baseline が気になる

なら、**検証目的で 1 回だけ** `subtract_baseline_from_fits()` を試す価値はあります。
その場合も、元の cube は残してください。

#### `RMS` と `BASE_RMS` はどちらを見るべきか

このレシピでは、PS pipeline の直後なので、まず見るのは `RMS` です。
`BASE_RMS` は baseline 解析を cube に掛けたあとに現れる量で、意味が違います。

- `RMS`
  - グリッディング後の画素雑音の目安
- `BASE_RMS`
  - baseline fit の line-free 残差に基づく目安

### 変形例 1. `TR*` で出したい

```python
res = run_ps_mapping_pipeline(
    scantable=sc_all_coadded_12co,
    config=ps_cfg,
    output_fits='sc_all_coadded_12co_tr_cube.fits',
    out_scale='TR*',
    dv_kms=0.2,
    vmin_kms=-30.0,
    vmax_kms=55.0,
)
```

このとき、主 cube の `BUNIT` は `K` のままですが、header の `TEMPSCAL` が `TR*` になります。

### 変形例 2. とりあえず cube だけ作って mask は後日考える

```python
res = run_ps_mapping_pipeline(
    scantable=sc_all_coadded_12co,
    config=ps_cfg,
    output_fits='sc_all_coadded_12co_cube.fits',
    out_scale='TA*',
)
```

このように、まず cube だけを作って、`RMS` や代表スペクトルを見ながら後で `make_3d_mask_for_existing_fits()` を掛けても問題ありません。

### 変形例 3. 速度範囲を切らずに全部残す

```python
res = run_ps_mapping_pipeline(
    scantable=sc_all_coadded_12co,
    config=ps_cfg,
    output_fits='sc_all_coadded_12co_fullv_cube.fits',
    out_scale='TA*',
    dv_kms=0.2,
    vmin_kms=None,
    vmax_kms=None,
)
```

探索観測や、line の速度位置をまだ決め打ちしたくないときに向いています。

### このレシピの要点

- baseline 済み scantable から始めるなら、最初の入口は `run_ps_mapping_pipeline()`
- cube 化のあとで `RMS/HIT/TSYS/TINT` を確認する
- signal mask は `make_3d_mask_for_existing_fits()` で後段で作る
- baseline の再実行は、必要性を確認してから限定的に行う
