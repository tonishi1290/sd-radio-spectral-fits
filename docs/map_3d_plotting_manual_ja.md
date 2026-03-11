# `sd_radio_spectral_fits.map_3d.plotting` 詳細説明書（日本語・preliminary）

> **重要な注記**  
> 本文書および現在の `plotting.py` 実装は **preliminary** です。今後、API 名、引数名、返り値、内部メタデータ、CLI、拡張 HDU 名の扱い、RGB 合成、reprojection の扱いなどが**大きく変更される可能性があります**。  
> 特に、表示スケールの引数名、WCS 不一致時の扱い、mask 解釈、`map_3d` パイプライン出力との結び付きは、今後の利用経験に応じて再設計されうる部分です。

---

## 1. このモジュールの目的

`sd_radio_spectral_fits.map_3d.plotting` は、電波天文の 2D / 3D FITS データ、および Python 上の `(header, data)` を、**天球座標付き 2D 図**として統一的に扱うための preliminary 可視化ユーティリティです。

主な目的は次の通りです。

- 2D FITS をそのまま WCSAxes で描く
- 3D cube から channel / velocity 範囲を指定して 2D map を生成する
- `PrimaryHDU` と extension HDU の両方を扱う
- `MOMENT0`, `RMS`, `BASE_RMS`, `HIT`, `TSYS`, `TINT`, `MASK3D` などを描く
- pseudo color と contour を同じ天球座標系で重ねる
- 複数データセットから複数 contour を重ねる
- RGB 合成を行う
- baseline 由来の **provisional moment** と `MASK3D` 由来の **final moment** を区別して扱う
- smoothing, normalize, beam 表示, colorbar, grid をまとめて扱う

このモジュールは、一般的な FITS viewer を目指すものではなく、**`map_3d` 系の出力をそのまま扱う可視化補助**として設計されています。

---

## 2. 公開 API の全体像

現行実装で、外部から主に使う関数・データ構造は次の通りです。

### 2.1 データ構造

#### `Map2D`
2D map を保持するコンテナです。

- `data`: 2D `numpy.ndarray`
- `header`: `astropy.io.fits.Header`
- `wcs`: `astropy.wcs.WCS`
- `unit`: `astropy.units.Unit` または `None`
- `meta`: 辞書。生成モード、入力元、smoothing 情報、normalize 情報などを格納

#### `CubeInput`
3D cube を内部で扱うためのコンテナです。

- `cube`: `spectral_cube.SpectralCube`
- `header`: FITS header
- `unit`: データ unit
- `meta`: 補助情報

#### `RGBMap`
RGB 合成画像を保持するコンテナです。

- `rgb`: shape `(ny, nx, 3)` の 3 次元配列
- `header`: 基準とする header
- `wcs`: 基準とする 2D WCS
- `meta`: 各色チャネルの生成情報など

### 2.2 入力解決

#### `resolve_map_input(...)`
2D 入力を `Map2D` に解決します。

#### `resolve_cube_input(...)`
3D 入力を `CubeInput` に解決します。

### 2.3 2D map 生成

#### `make_2d_map(...)`
2D 入力をそのまま `Map2D` にするか、3D cube から 2D map を新規生成します。

#### `make_provisional_moment(...)`
baseline 由来情報から provisional moment を生成します。

#### `make_final_moment(...)`
`MASK3D` などの final mask から final moment を生成します。

### 2.4 描画

#### `build_normalize(...)`
表示用の `ImageNormalize` を構築します。

#### `plot_map(...)`
`Map2D` を WCSAxes 上に pseudo color で描画し、必要なら contour, beam, colorbar を追加します。

#### `add_contours(...)`
既存の `ax` に contour を追加します。

#### `make_rgb_map(...)`
3 枚の 2D map から RGB 合成用 `RGBMap` を作ります。

#### `plot_rgb(...)`
RGB 合成画像を WCSAxes 上に描画します。

### 2.5 CLI

`python -m sd_radio_spectral_fits.map_3d.plotting ...` あるいは相当する直接実行で、簡易 CLI が使えます。

---

## 3. 何を入力できるか

このモジュールでは、**2D と 3D の入力パターンをかなり幅広く受ける**ようにしています。

### 3.1 2D 入力として受けられるもの

`resolve_map_input(...)` や `plot_map(...)` の `source` に、概ね次を与えられます。

1. **FITS ファイルパス**  
   例: `"cube.fits"`, `"map2d.fits"`

2. **FITS HDUList**  
   `astropy.io.fits.open(...)` で開いたもの

3. **2D HDU**  
   `PrimaryHDU`, `ImageHDU` など

4. **`(header, data)` のタプル**  
   Python 上で持っている 2D 配列と header

5. **`Map2D` 自体**  
   すでに本モジュールで生成済みのもの

6. **2D の `Projection` 相当オブジェクト**  
   `spectral_cube.Projection` など、2D WCS と配列を持つもの

### 3.2 3D 入力として受けられるもの

`resolve_cube_input(...)`, `make_2d_map(...)`, `make_provisional_moment(...)`, `make_final_moment(...)` では、概ね次を与えられます。

1. **3D FITS ファイルパス**
2. **3D HDUList / 3D HDU**
3. **`(header, data)` のタプル（3D）**
4. **`SpectralCube`**

### 3.3 `ext` の役割

`ext` は、FITS のどの HDU を使うかを指定します。

- `None`: 既定の HDU を使う
- 整数: HDU index
- 文字列: extension name

例:

```python
ext=0
ext=1
ext="MOMENT0"
ext="RMS"
ext="MASK3D"
```

### 3.4 重要な注意: `ext="MOMENT0"` は 2D である

`ext="MOMENT0"` のように、既に 2D の拡張を指定した場合、それは **3D cube から新たに moment0 を作る**のではなく、**既に出来上がっている 2D map を読む**動作になります。

したがって、

- `vel_range`
- `chan_range`

は、その 2D map を読むだけの場合には使われません。これらは **3D cube を 2D に落とすときだけ有効**です。

---

## 4. `map_3d` パッケージの出力をどう入力するか

本モジュールの実運用では、`map_3d` 系パイプラインの出力 FITS をそのまま扱うことが重要です。

### 4.1 想定している典型的な出力

実装およびハンドオーバー上で強く想定しているのは、次のような出力です。

- 3D spectral cube 本体
- 2D extension
  - `MOMENT0`
  - `RMS`
  - `BASE_RMS`
  - `HIT`
  - `TSYS`
  - `TINT`
- 3D / 1D mask 系 extension
  - `MASK3D`
  - `LINEFREE`
  - `BASESUP3D`
  - `LINECAND3D`

ただし、**これらの extension 名や HDU 構成は preliminary であり、将来変更される可能性があります**。

### 4.2 もっとも典型的な読み方

#### 4.2.1 既に出来ている 2D map を読む

```python
map2d = make_2d_map("ps_12co_cube_masked.fits", ext="MOMENT0", mode="map")
```

この場合、`MOMENT0` 拡張をそのまま 2D map として読みます。

同様に、

```python
rms_map = make_2d_map("ps_12co_cube_masked.fits", ext="RMS", mode="map")
hit_map = make_2d_map("ps_12co_cube_masked.fits", ext="HIT", mode="map")
```

のように使えます。

#### 4.2.2 3D cube 本体から新たに moment0 を作る

```python
map2d = make_2d_map(
    "ps_12co_cube_masked.fits",
    ext=0,
    mode="moment0",
    vel_range=(0.0, 20.0),
    spectral_unit="km/s",
)
```

ここで `ext=0` は、3D cube が入っている HDU を表す例です。実際のファイル構造に応じて変更してください。

#### 4.2.3 provisional moment を作る

```python
pmap = make_provisional_moment("ps_12co_cube_masked.fits", ext=0)
```

#### 4.2.4 final moment を作る

```python
fmap = make_final_moment("ps_12co_cube_masked.fits", ext=0)
```

### 4.3 `LINECAND3D -> BASESUP3D -> LINEFREE 補集合` の意味

これは provisional moment を作るときの **自動探索の優先順位**です。

#### 第一候補: `LINECAND3D`
もし `LINECAND3D` があるなら、これは「**ここが line candidate である**」という 3D マスクを表すとみなし、これを最優先で使います。

#### 第二候補: `BASESUP3D`
`LINECAND3D` が無いが `BASESUP3D` があるなら、これは「**baseline 推定に使ってよい領域**」とみなします。したがって、その**補集合**、つまり baseline support ではない側を provisional な signal 候補として使います。

#### 第三候補: `LINEFREE`
`LINECAND3D` も `BASESUP3D` も無いが `LINEFREE` があるなら、これは「**line-free なチャネル**」を表します。したがって、その**補集合**、つまり line-free ではないチャネルを provisional な signal 候補として使います。

要するに、

1. signal を直接指定するマスクがあればそれを使う
2. 無ければ baseline 側の情報から signal 側を補集合で作る
3. さらに無ければ line-free 情報の補集合を使う

という意味です。

### 4.4 `MASK3D` の意味

`MASK3D` は final signal mask を想定しています。よって `make_final_moment(...)` は provisional な試行ではなく、**最終的な signal とみなした領域**で積分するための関数です。

### 4.5 2D extension と 3D 再計算を混同しない

とても重要です。

- `ext="MOMENT0"` は、**すでに出来ている 2D map を読む**
- `ext=0, mode="moment0", vel_range=(...)` は、**3D cube から新しく moment0 を作る**

この 2 つは同じではありません。

---

## 5. 2D map 生成: `make_2d_map(...)`

### 5.1 役割

`make_2d_map(...)` は、本モジュールの中核関数です。

- 入力が 2D なら、そのまま `Map2D` にする
- 入力が 3D なら、指定モードで 2D に落とす
- 必要に応じて mask, smoothing を適用する
- meta に生成情報を残す

### 5.2 シグネチャ

現行実装の概念的シグネチャは次です。

```python
make_2d_map(
    source,
    *,
    ext=None,
    mode="moment0",
    chan_range=None,
    vel_range=None,
    spectral_unit="km/s",
    mask=None,
    mask_mode=None,
    zero_fill=False,
    nan_fill=True,
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
    fill_value=0.0,
    linefree_ext="LINEFREE",
    linecand_ext="LINECAND3D",
    basesup_ext="BASESUP3D",
    final_mask_ext="MASK3D",
)
```

### 5.3 パラメーター全解説

#### `source`
入力データです。ファイルパス、HDU, `(header, data)`, `SpectralCube`, `Map2D` などを受けます。

#### `ext`
使う HDU を指定します。`None`, 整数 index, extension name を受けます。

#### `mode`
3D から 2D をどう作るかを指定します。現行実装では概ね次を想定します。

- `"map"`, `"2d"`, `"identity"`  
  2D 入力をそのまま使う
- `"channel"`, `"slice"`, `"channel_slice"` 系  
  単一チャネルを取り出す
- `"channel_sum"`  
  チャネル範囲和
- `"channel_mean"`  
  チャネル範囲平均
- `"moment0"`  
  速度軸あるいは spectral axis に沿って積分
- `"provisional_moment"`  
  baseline 由来 mask で積分
- `"final_moment"`  
  `MASK3D` など final mask で積分

注意として、**入力が既に 2D なら、`mode="moment0"` であっても再積分はせず、その 2D map を読むだけ**です。

#### `chan_range`
チャネル index の範囲です。通常 `(lo, hi)` を与えます。両端を含む設計です。`vel_range` より優先されます。

用途:
- `channel_sum`
- `channel_mean`
- `channel_slice`
- `moment0` のチャネル直接指定

#### `vel_range`
速度範囲です。通常 `(v_lo, v_hi)` を与えます。`spectral_unit` に基づいて解釈されます。

用途:
- `moment0`
- provisional / final moment の積分範囲制限

注意:
- 入力が 2D の場合は無効です
- cube の spectral axis が速度変換できない場合は、実装により native axis へフォールバックすることがあります

#### `spectral_unit`
`vel_range` をどの単位で与えるかを指定します。典型例は `"km/s"` です。

#### `mask`
追加で適用する mask です。配列または解釈可能なオブジェクトを想定します。用途は preliminary で、将来変更の可能性があります。

#### `mask_mode`
`mask` の意味づけを指定するための補助引数です。現段階では将来拡張用の性格が強いです。

#### `zero_fill`
mask 外や除外領域を 0 で埋めるかを指定します。

- `True`: 0 埋め
- `False`: 0 埋めしない

積分や平均に対する影響が大きいので注意が必要です。

#### `nan_fill`
mask 外や除外領域を `NaN` で埋めるかを指定します。

- `True`: `NaN` 埋め
- `False`: `NaN` 埋めしない

通常は `nan_fill=True` の方が「無効領域を見えやすく保持する」点で安全です。

#### `smooth_fwhm_arcsec`
追加 smoothing kernel の FWHM を秒角で指定します。これは「追加でこれだけ平滑化したい」ときに使います。

#### `target_hpbw_arcsec`
最終的に目指す HPBW を秒角で指定します。`orig_hpbw_arcsec` または header の `BMAJ/BMIN` が分かる場合、

```text
kernel^2 = target^2 - original^2
```

で追加 kernel を計算します。

#### `orig_hpbw_arcsec`
元の HPBW が header から推定できない場合に、明示的に秒角で指定します。

#### `fill_value`
必要に応じて使われる埋め値です。実装の将来変更余地があります。

#### `linefree_ext`
`LINEFREE` extension 名です。既定は `"LINEFREE"` です。

#### `linecand_ext`
`LINECAND3D` extension 名です。既定は `"LINECAND3D"` です。

#### `basesup_ext`
`BASESUP3D` extension 名です。既定は `"BASESUP3D"` です。

#### `final_mask_ext`
`MASK3D` extension 名です。既定は `"MASK3D"` です。

### 5.4 smoothing の重要な注意

#### `smooth_fwhm_arcsec` と `target_hpbw_arcsec` は概念が違う

- `smooth_fwhm_arcsec`: 追加 kernel の大きさを直接指定
- `target_hpbw_arcsec`: 最終 HPBW 目標値を指定

#### header に beam が無い場合

header の `BMAJ/BMIN` が無い場合、`target_hpbw_arcsec` 単独では追加 kernel を厳密に決められません。この場合は、

- `orig_hpbw_arcsec` を明示する
- それが無ければ warning の上で `target_hpbw_arcsec` を追加 kernel とみなす

という挙動になります。

#### 目標 HPBW が元より小さい場合

追加 smoothing は不要なので、実装上は元データのまま返します。

### 5.5 返り値

返り値は `Map2D` です。

`meta` には概ね次のような情報が入ります。

- 入力元
- `ext`
- `mode`
- `chan_range`
- `vel_range`
- `spectral_unit`
- smoothing 情報
- provisional/final の種別

---

## 6. provisional moment: `make_provisional_moment(...)`

### 6.1 役割

baseline 情報から provisional な signal map を作ります。

### 6.2 シグネチャ

概念的には次です。

```python
make_provisional_moment(
    source,
    *,
    ext=None,
    linefree_ext="LINEFREE",
    basesup_ext="BASESUP3D",
    linecand_ext="LINECAND3D",
    prefer="auto",
    vel_range=None,
    chan_range=None,
    spectral_unit="km/s",
    zero_fill=False,
    nan_fill=True,
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
)
```

### 6.3 パラメーター全解説

#### `source`, `ext`
`make_2d_map(...)` と同様です。通常は 3D cube 本体を指します。

#### `linefree_ext`, `basesup_ext`, `linecand_ext`
provisional mask を探す拡張名です。

#### `prefer`
どの情報を優先して provisional mask を作るかを指定します。

- `"auto"`: `LINECAND3D -> BASESUP3D -> LINEFREE 補集合` の順に探索
- `"linecand3d"` 相当: `LINECAND3D` を優先
- `"basesup3d"` 相当: `BASESUP3D` を優先
- `"linefree"` 相当: `LINEFREE` を優先

`auto` の意味は、前節の説明が最重要です。

#### `vel_range`, `chan_range`, `spectral_unit`
積分範囲の制限です。

#### `zero_fill`, `nan_fill`
除外領域の埋め方です。

- `zero_fill=True`: 除外領域を 0 として積分に参加させる
- `nan_fill=True`: 除外領域を無効値として扱う

通常は `nan_fill=True` の方が安全ですが、見せ方や後段処理によって選びます。

#### `smooth_fwhm_arcsec`, `target_hpbw_arcsec`, `orig_hpbw_arcsec`
生成後の 2D map に対する smoothing 条件です。

### 6.4 返り値と meta

返り値は `Map2D` です。`meta` には少なくとも「これは provisional である」こと、どの extension が使われたか、どの `prefer` で決まったかが入るのが望ましい設計です。

---

## 7. final moment: `make_final_moment(...)`

### 7.1 役割

`MASK3D` など final signal mask を使って final moment を作ります。

### 7.2 シグネチャ

概念的には次です。

```python
make_final_moment(
    source,
    *,
    ext=None,
    final_mask_ext="MASK3D",
    vel_range=None,
    chan_range=None,
    spectral_unit="km/s",
    zero_fill=False,
    nan_fill=True,
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
)
```

### 7.3 パラメーター全解説

#### `final_mask_ext`
final signal mask の extension 名です。通常は `"MASK3D"` を想定します。

#### `vel_range`, `chan_range`, `spectral_unit`
積分範囲の制限です。`moment0` と同様に使います。

#### `zero_fill`, `nan_fill`
無効領域の埋め方です。

#### smoothing 群
2D 化後に適用します。

### 7.4 provisional との違い

- provisional: baseline 側情報や line-free 情報から「仮の signal 候補」を作る
- final: final signal mask を使って「確定版 signal」として積分する

この違いは解析上非常に重要です。

---

## 8. normalize: `build_normalize(...)`

### 8.1 役割

`ImageNormalize` を構築し、pseudo color の表示スケールを決めます。

### 8.2 現行シグネチャ

```python
build_normalize(
    data,
    *,
    mode="asinh",
    percentile=None,
    vmin=None,
    vmax=None,
    stretch_a=0.1,
    power_gamma=1.0,
    invalid=np.nan,
)
```

### 8.3 重要な注意

現行 preliminary 実装では、表示スケール下限・上限の引数名が **`vmin`, `vmax`** です。しかし、ここでの `v` は **velocity ではありません**。意味は「color scale の下限・上限」です。

したがって、

- `vel_range`: 速度範囲
- `vmin`, `vmax`: 表示スケール下限・上限

です。今後、より分かりやすい `cmin/cmax` 等へ改名される可能性があります。

### 8.4 パラメーター全解説

#### `data`
表示対象配列です。有限値部分から自動範囲を計算します。

#### `mode`
normalize / stretch の手法です。

- `"linear"`
- `"sqrt"`
- `"log"`
- `"asinh"`
- `"power"`

#### `percentile`
自動 clip 用の percentile 範囲です。例: `(1.0, 99.5)`。

#### `vmin`, `vmax`
表示スケールの下限・上限です。`percentile` が与えられても、これらを明示した側はその値が優先されます。

#### `stretch_a`
`asinh` 用パラメーターです。小さいほど弱い構造を強調しやすい傾向があります。

#### `power_gamma`
`power` stretch 用パラメーターです。

#### `invalid`
無効値の扱いです。

### 8.5 挙動の優先順位

概ね次の順です。

1. `percentile` が与えられれば自動 `vmin/vmax` 候補を作る
2. ただし、明示 `vmin/vmax` があればその側を優先
3. どちらも無ければ有限値範囲から決める

### 8.6 実務的な使い分け

- `MOMENT0`: `asinh` + percentile clip が使いやすい
- `RMS`, `BASE_RMS`: `linear` または `log`
- `HIT`, `TINT`: `linear` が分かりやすい

---

## 9. 描画本体: `plot_map(...)`

### 9.1 役割

`Map2D` または 2D 入力を pseudo color で描き、必要なら colorbar, contour, beam, grid を追加します。

### 9.2 現行シグネチャ

```python
plot_map(
    source=None,
    *,
    data=None,
    header=None,
    ext=None,
    ax=None,
    projection=None,
    cmap="viridis",
    norm=None,
    norm_mode="asinh",
    norm_percentile=None,
    vmin=None,
    vmax=None,
    stretch_a=0.1,
    power_gamma=1.0,
    colorbar=True,
    colorbar_label=None,
    title=None,
    origin="lower",
    contours=None,
    grid=True,
    beam=None,
    xlabel=None,
    ylabel=None,
    show=True,
    figsize=(10, 8),
)
```

### 9.3 パラメーター全解説

#### `source`
描画対象です。`Map2D`, ファイルパス, HDU, `(header, data)` などを受けます。

#### `data`, `header`
`source` を使わずに生配列を渡したい場合に使います。`header` は 2D WCS を持つ必要があります。

#### `ext`
`source` が FITS のとき、どの HDU を使うかを指定します。

#### `ax`
既存の Matplotlib/Astropy WCSAxes を使う場合に指定します。`None` なら新しく Figure / Axes を作ります。

#### `projection`
描画に使う WCS を明示したいときに指定します。通常は `map2d.wcs` が使われます。

#### `cmap`
Matplotlib の colormap 名です。例: `"viridis"`, `"turbo"`, `"magma"`。

#### `norm`
外部で作った `ImageNormalize` をそのまま使う場合に指定します。これを指定すると `norm_mode`, `norm_percentile`, `vmin`, `vmax`, `stretch_a`, `power_gamma` は内部 normalize 生成には使われません。

#### `norm_mode`
内部で normalize を作るときのモードです。`build_normalize(...)` の `mode` に対応します。

#### `norm_percentile`
内部 normalize 用の percentile clip です。

#### `vmin`, `vmax`
表示スケール下限・上限です。速度ではありません。

#### `stretch_a`
`asinh` 用パラメーターです。

#### `power_gamma`
`power` 用パラメーターです。

#### `colorbar`
`True` なら colorbar を描きます。

#### `colorbar_label`
colorbar のラベルを明示します。`None` なら unit などから既定ラベルを推定します。

#### `title`
図タイトルです。

#### `origin`
`imshow` の `origin` です。通常は `"lower"`。

#### `contours`
contour 層のリストです。各要素は辞書で指定します。詳細は次節を参照してください。

#### `grid`
座標 grid を表示するかどうかです。

#### `beam`
beam 表示設定です。

- `None`: 描かない
- `"header"`, `"auto"`: `smooth_info` または header の `BMAJ/BMIN/BPA` から描く
- 辞書: beam 形状や位置を明示指定

#### `xlabel`, `ylabel`
軸ラベルの上書きです。`None` なら WCS から既定ラベルを決めます。

#### `show`
`True` なら、新規作成 Figure のとき `plt.show()` を呼びます。

#### `figsize`
Figure サイズです。

### 9.4 返り値

辞書を返します。主なキー:

- `fig`
- `ax`
- `image`
- `colorbar`
- `contours`
- `beam`
- `map2d`

返された `map2d.meta["normalize_info"]` に normalize 情報が追加されます。

---

## 10. contour: `add_contours(...)` と `contours` 辞書

### 10.1 役割

`plot_map(..., contours=[...])` で複数 contour を重ねられます。

### 10.2 contour 辞書の全パラメーター

各 contour layer は辞書で指定します。現行実装で意味を持つ主なキーは次です。

#### 入力指定

- `source`: 入力元
- `data`: 2D 生配列
- `header`: 2D header
- `ext`: FITS extension

#### 3D から contour map を作るための指定

- `mode`
- `chan_range`
- `vel_range`
- `spectral_unit`
- `mask`
- `mask_mode`
- `zero_fill`
- `nan_fill`
- `smooth_fwhm_arcsec`
- `target_hpbw_arcsec`
- `orig_hpbw_arcsec`
- `linefree_ext`
- `linecand_ext`
- `basesup_ext`
- `final_mask_ext`

これらが一つでも入っていると、その contour layer は `make_2d_map(...)` を通して 2D 化されます。

#### 線の見た目

- `levels`: contour level の配列、または `"auto"`
- `colors`: 線色
- `linewidths`: 線幅
- `linestyles`: 線種
- `alpha`: 透明度
- `label`: 凡例用ラベル

### 10.3 `levels="auto"`
データの最大値などから自動 level を作ります。厳密な level の選び方は将来変更の可能性があります。

### 10.4 同じ map を contour にする場合

`plot_map(...)` のベース map と同じものに contour を重ねるだけなら、最小限の指定で済みます。

```python
contours=[{"levels": "auto", "colors": "white"}]
```

### 10.5 別データセットから contour を重ねる場合

```python
contours=[
    {
        "source": "other_cube.fits",
        "ext": 0,
        "mode": "moment0",
        "vel_range": (0, 20),
        "spectral_unit": "km/s",
        "levels": [10, 20, 40, 80],
        "colors": "cyan",
        "linewidths": 0.8,
    }
]
```

---

## 11. RGB 合成: `make_rgb_map(...)` と `plot_rgb(...)`

### 11.1 `make_rgb_map(...)` の役割

3 枚の 2D map を正規化して RGB に割り当てます。

### 11.2 現行シグネチャ

```python
make_rgb_map(
    red,
    green,
    blue,
    *,
    red_norm=None,
    green_norm=None,
    blue_norm=None,
)
```

### 11.3 `red`, `green`, `blue` に与えられるもの

各色チャネルには、次のいずれかを与えられます。

1. `Map2D`
2. 2D 入力そのもの
3. `make_2d_map(...)` 用の辞書

たとえば、

```python
red={"source": "cube.fits", "ext": 0, "mode": "moment0", "vel_range": (0, 5)}
```

のように書けます。

### 11.4 `red_norm`, `green_norm`, `blue_norm`
各色チャネルの正規化辞書です。現行実装で使う主なキー:

- `norm_mode`
- `percentile`
- `vmin`
- `vmax`
- `stretch_a`
- `power_gamma`

### 11.5 重要な注意: reprojection はまだ行わない

現行 preliminary 実装では、RGB 各チャネルの **shape が一致すること**を要求します。WCS が厳密一致しない場合も、まずは warning を出した上でそのまま stack します。**reproject はまだ実装していません**。

### 11.6 `plot_rgb(...)` のパラメーター

```python
plot_rgb(
    rgb_source=None,
    *,
    red=None,
    green=None,
    blue=None,
    ax=None,
    title=None,
    grid=True,
    xlabel=None,
    ylabel=None,
    show=True,
    figsize=(10, 8),
    **rgb_kwargs,
)
```

- `rgb_source`: 既に作った `RGBMap`
- `red`, `green`, `blue`: `rgb_source` を作るための入力
- `ax`, `title`, `grid`, `xlabel`, `ylabel`, `show`, `figsize`: `plot_map(...)` と同様
- `**rgb_kwargs`: `make_rgb_map(...)` にそのまま渡されます

---

## 12. beam 表示の詳細

### 12.1 `beam="header"`
優先順位は概ね次です。

1. `map2d.meta["smooth_info"]["effective_hpbw_arcsec"]` があればそれを使う
2. 無ければ header の `BMAJ/BMIN/BPA` を使う
3. どちらも無ければ描かない

### 12.2 beam 辞書の全パラメーター

`beam` に辞書を与える場合の主なキー:

- `major`: 長軸
- `minor`: 短軸。省略時は `major`
- `angle`: 角度
- `units`: `"arcsec"`, `"deg"`, `"pix"`
- `x`, `y`: 描画位置（ピクセル）
- `facecolor`
- `edgecolor`
- `linewidth`
- `alpha`

### 12.3 HPBW が header に無い場合

観測データによっては `BMAJ/BMIN` が入っていないことがあります。この場合、beam 表示や `target_hpbw_arcsec` の厳密計算には `orig_hpbw_arcsec` の明示が有効です。

---

## 13. 銀河座標と赤道座標

現行実装では、WCS の celestial 軸情報から、概ね次を判定します。

- `GLON/GLAT` 系: 銀河座標
- `RA/DEC` 系: 赤道座標

これに応じて、既定の軸ラベルを自動選択します。

---

## 14. CLI 詳細

### 14.1 役割

簡単な一枚図をコマンドラインから作るための preliminary CLI です。

### 14.2 主な引数

現行実装で確認できる主要引数は次です。

- 位置引数
  - `source`: 入力 FITS ファイル

- 入力・2D 化
  - `--ext`
  - `--mode`
  - `--chan-range`
  - `--vel-range`
  - `--spectral-unit`

- smoothing
  - `--smooth-fwhm-arcsec`
  - `--target-hpbw-arcsec`
  - `--orig-hpbw-arcsec`

- 表示
  - `--cmap`
  - `--norm-mode`
  - `--norm-percentile`
  - `--vmin`
  - `--vmax`
  - `--stretch-a`
  - `--power-gamma`
  - `--title`
  - `--output`
  - `--dpi`

- contour
  - `--contour-levels`
  - `--contour-color`
  - `--contour-linewidth`

- mask / fill
  - `--zero-fill`
  - `--nan-fill` 相当の引数群

### 14.3 典型例

```bash
python -m sd_radio_spectral_fits.map_3d.plotting ps_12co_cube.fits \
  --ext 0 \
  --mode moment0 \
  --vel-range 0,20 \
  --spectral-unit km/s \
  --target-hpbw-arcsec 400 \
  --orig-hpbw-arcsec 350 \
  --cmap turbo \
  --norm-mode asinh \
  --norm-percentile 1,99.5 \
  --title "12CO moment0" \
  --output moment0.png
```

---

## 15. `map_3d` 出力の受け渡しパターンまとめ

### 15.1 既存 2D extension を描く

対象: `MOMENT0`, `RMS`, `BASE_RMS`, `HIT`, `TSYS`, `TINT` など

```python
m = make_2d_map("cube_masked.fits", ext="MOMENT0", mode="map")
plot_map(m, cmap="turbo", norm_mode="asinh", norm_percentile=(1, 99.5))
```

### 15.2 cube から新規再計算する

対象: 任意速度範囲の moment0 など

```python
m = make_2d_map(
    "cube.fits",
    ext=0,
    mode="moment0",
    vel_range=(0, 20),
    spectral_unit="km/s",
)
```

### 15.3 provisional moment を作る

```python
m = make_provisional_moment("cube_masked.fits", ext=0, prefer="auto")
```

### 15.4 final moment を作る

```python
m = make_final_moment("cube_masked.fits", ext=0, final_mask_ext="MASK3D")
```

### 15.5 Python で作った `(header, data)` を描く

```python
m = plot_map(data=my2d, header=my_header, cmap="viridis")
```

### 15.6 cube から RGB を作る

```python
rgb = make_rgb_map(
    red={"source": "cube.fits", "ext": 0, "mode": "moment0", "vel_range": (-10, 0)},
    green={"source": "cube.fits", "ext": 0, "mode": "moment0", "vel_range": (0, 10)},
    blue={"source": "cube.fits", "ext": 0, "mode": "moment0", "vel_range": (10, 20)},
)
plot_rgb(rgb, title="RGB moment0")
```

---

## 16. よくある混乱点

### 16.1 `ext="MOMENT0"` に `vel_range` を付けても意味があるか

通常はありません。すでに 2D map を読んでいるだけだからです。

### 16.2 `vmin/vmax` は速度か

違います。表示スケールの下限・上限です。

### 16.3 `target_hpbw_arcsec` だけでよいか

元ビームが分かっているならよいですが、header に beam が無い場合は `orig_hpbw_arcsec` を与えた方が安全です。

### 16.4 RGB の WCS が少し違ってもよいか

現行実装では warning の上でそのまま stack します。厳密な reproject はまだありません。

### 16.5 provisional と final は何が違うか

provisional は baseline 側情報から作る仮の signal map、final は final signal mask に基づく確定版 map です。

---

## 17. 将来の改良候補

本モジュールは preliminary なので、今後特に次の改良がありえます。

- `vmin/vmax` を `cmin/cmax` などへ改名
- reprojection の導入
- contour / RGB での WCS 整合性チェック強化
- background-noise ベースの normalize
- beam の自動描画強化
- region overlay
- interactive viewer
- `map_3d` 出力とのより厳密な統合
- extension 名の標準化

---

## 18. 最後に

このモジュールは、**`map_3d` の出力を、天球座標上で比較的少ない手間で確かめる**ことを第一目的にしています。厳密な FITS viewer というより、解析・開発・QC のための作業用描画ツールとして使うのが自然です。

一方で、provisional/final moment の区別、beam/HPBW、normalize、mask 解釈など、科学的解釈に直接関わる部分を含むため、preliminary 実装であることを常に意識して使う必要があります。
