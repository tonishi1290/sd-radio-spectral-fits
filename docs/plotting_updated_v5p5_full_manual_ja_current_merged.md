# plotting_updated_v5p5 詳細説明書（完全版・日本語）

> 対象コード: **`plotting_updated_v5p5.py`**  
> 本書は `plotting_updated_v5p5.py` を基準に、公開 API、主要データ構造、layer 辞書、CLI オプション、設計思想、注意点までを、できる限り省略せずに整理した詳細説明書です。  
> 既存の詳細説明書の書きぶりと説明の粒度をできるだけ維持しつつ、内容は **現行 v5p5 仕様** に統一しています。特に、**`quicklook(contours=...)`**、**contour level の `auto / fraction / manual / rms` 整理**、**CLI contour 条件の改善**、**`plot_scene()` / `quicklook()` を新規コードで推奨する方針**を反映しています。

---

## 0. 現在の推奨 API と役割分担

- **最初の 1 枚を出す** → `quicklook()`
- **複数 layer や拡張性を重視する** → `plot_scene()`
- **既存コードの移行・互換性維持** → `plot_map()` / `plot_rgb()`

### 0.1 `plot_map()` の位置づけ

`plot_map()` は現時点でも十分に使えますが、設計上は **過去との互換性を維持するために重要な API** という位置づけです。したがって、**新しいコードを今後書くときは、まず `quicklook()` と `plot_scene()` を優先して考える**ことを推奨します。

`plot_map()` を残している理由:

- 既存の解析コードの移行コストを下げるため
- 一枚絵 + contour の従来スタイルを壊さないため
- `make_2d_map() + plot_map()` という既存ワークフローを尊重するため

ただし、今後の新機能は基本的に `plot_scene()` 側に集約していく方針です。

### 0.2 現在の contour level 仕様

- `levels="auto"` または `level_mode="auto"` は、**正のピークに対する 10, 30, 50, 70, 90%** を意味します。
- `level_mode="fraction"` は `fraction_levels=[...]` をピーク値に掛けます。
- `level_mode="manual"` は `levels=[...]` を factor として扱い、必要なら `level_scale` を掛けます。  
  例: `levels=[1,2,3,4], level_scale=3` → `3,6,9,12`
- `level_mode="rms"` は `rms * sigma_scale * sigma_levels` で level を作ります。
- `level_mode="sigma"` は後方互換 alias として受理されます。
- `quicklook()` は **`contours=` に正式対応**しています。

### 0.3 本書の読み方

- **まず最短で 1 枚出したい** → 第 3 章
- **入力形式や 2D / 3D の扱いを知りたい** → 第 4～6 章
- **各関数の全パラメータを引きたい** → 第 7 章
- **CLI を使いたい** → 第 8 章後半と第 10 章
- **方針や注意点を確認したい** → 第 11 章以降

---

## 1. このモジュールは何をするものか

`plotting_updated_v5p5.py` は、電波天文学の 2D / 3D FITS データを **天球座標つきで描画するための WCS-aware plotting helper** です。主な用途は次のとおりです。

- 2D FITS の単純表示
- 3D cube からの 2D map 作成
- 速度積分図（moment0）や channel slice の作成
- provisional / final moment の作成
- pseudo-color 表示
- contour 重ね描き
- RGB 合成
- beam, scalebar, north/east arrow の描画
- marker / catalog / text の重ね描き
- optical / IR 背景画像の上に radio contour を重ねる図の作成
- Python API と CLI の両対応

設計上の重要な思想は次の 4 点です。

1. **1枚を簡単に出す**ための入口として `quicklook()` を用意する。  
2. **高機能な重ね描き**は `plot_scene()` にまとめる。  
3. 2D / 3D、ファイル / メモリ上データを **できるだけ同じ流儀で扱う**。  
4. contour と image overlay を分け、**必要な場面でだけ reprojection を使う**。

---

## 2. 依存関係

必須に近いもの:

- `numpy`
- `astropy`
  - `astropy.io.fits`
  - `astropy.wcs`
  - `astropy.visualization`
  - `astropy.coordinates`
  - `astropy.nddata.Cutout2D`
  - `astropy.convolution`
- `matplotlib`
- `spectral_cube`

任意:

- `reproject`
  - image overlay を base WCS に再投影したいときのみ必要

補足:

- contour overlay 自体は `reproject` なしでも使えます。contour は WCSAxes の transform で重ねる設計です。
- RGB は **同じ 2D shape 前提**です。v5 系では RGB の本格 reprojection はまだ入っていません。

---

## 3. まず最短で使う方法

### 3.1 2D FITS をすぐ見る

```python
import sd_radio_spectral_fits.map_3d.plotting as pl

pl.quicklook('moment0.fits')
```

### 3.2 3D cube から moment0 を作って見る

```python
pl.quicklook('cube.fits', mode='moment0', vel_range=(20, 40))
```

### 3.3 Python 上の `(header, data)` を直接描く

```python
pl.quicklook((header, data))
```

### 3.4 背景画像 + contour

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits'},
    overlays=[
        {'kind': 'contour', 'source': 'co_cube.fits', 'mode': 'moment0', 'vel_range': (20, 40), 'colors': 'cyan'},
        {'kind': 'beam', 'beam': 'auto'},
    ],
    title='Optical + CO contour',
)
```

---

## 4. 入力として何を渡せるか

このモジュールは、ユーザーが「ファイルの場合もあるし、Python 上の header/data の場合もあるし、2D/3D もある」という状況を前提に作られています。入力は大きく **2D 系** と **3D 系** に分かれます。

### 4.1 2D 入力として受けられるもの

`resolve_map_input()` が受ける主な形式:

- `Map2D`
- FITS file path (`str` / `Path`)
- `fits.HDUList`
- 2D `PrimaryHDU` / `ImageHDU` / `CompImageHDU`
- `(header, data)` タプル
- `spectral_cube` の `Projection`
- `.wcs`, `.shape`, `.ndim == 2` を持つ projection-like object

### 4.2 3D 入力として受けられるもの

`resolve_cube_input()` が受ける主な形式:

- `SpectralCube`
- FITS file path (`str` / `Path`)
- `fits.HDUList`
- 3D `PrimaryHDU` / `ImageHDU` / `CompImageHDU`
- `(header, data)` タプル

### 4.3 `ext=None` の意味

FITS file や `HDUList` を渡し、`ext=None` にすると、**最初に data を持つ HDU** が選ばれます。

### 4.4 2D passthrough の挙動

`make_2d_map()` に対して `mode='moment0'` を指定しても、選ばれた HDU がすでに 2D であれば、その 2D map をそのまま通します。たとえば `ext='MOMENT0'` を選んだ場合です。

このとき:

- `chan_range`
- `vel_range`
- `mask`
- `mask_mode`

などの **3D reduction 用パラメータは無視**されます。警告も出ます。

---

## 5. 主要データ構造

### 5.1 `Map2D`

```python
@dataclass
class Map2D:
    data: np.ndarray
    header: fits.Header
    wcs: WCS
    unit: Optional[u.UnitBase] = None
    meta: Dict[str, Any] = field(default_factory=dict)
```

意味:

- `data` : 2D array
- `header` : 2D celestial FITS header
- `wcs` : celestial WCS
- `unit` : 物理単位
- `meta` : reduction 情報や smoothing 情報など

### 5.2 `RGBMap`

```python
@dataclass
class RGBMap:
    rgb: np.ndarray
    header: fits.Header
    wcs: WCS
    meta: Dict[str, Any] = field(default_factory=dict)
```

- `rgb.shape == (ny, nx, 3)`
- `header`, `wcs` は RGB 全体の座標基準

### 5.3 `CubeInput`

```python
@dataclass
class CubeInput:
    cube: SpectralCube
    header: fits.Header
    unit: Optional[u.UnitBase]
    meta: Dict[str, Any] = field(default_factory=dict)
```

---

## 6. 公開 API の全体像

公開 API (`__all__`) に含まれている主要関数は次です。

補足として、現時点では **新規コードでは `quicklook()` と `plot_scene()` を優先し、`plot_map()` は互換 API として扱う**のが推奨です。

- `resolve_map_input`
- `resolve_cube_input`
- `build_normalize`
- `estimate_rms_robust`
- `compute_contour_levels`
- `apply_gaussian_smoothing_2d`
- `crop_map2d`
- `describe_map2d`
- `describe_source`
- `plotting_help_text`
- `print_plotting_help`
- `make_2d_map`
- `make_provisional_moment`
- `make_final_moment`
- `add_contours`
- `plot_map`
- `make_rgb_map`
- `plot_rgb`
- `add_scalebar`
- `add_north_arrow`
- `plot_scene`
- `quicklook`
- `main`

以下、すべて詳しく説明します。

---

# 7. 関数リファレンス（完全版）

## 7.1 `resolve_map_input()`

```python
resolve_map_input(source=None, *, data=None, header=None, ext=None) -> Map2D
```

2D 入力を `Map2D` に正規化します。

### パラメータ

#### `source`
2D source。本関数が解釈できる型は次です。

- `Map2D`
- FITS path
- `fits.HDUList`
- 2D HDU
- `(header, data)` タプル
- `Projection`
- projection-like object

#### `data`
NumPy 配列。`header` とセットで使います。`data` を渡した場合は **2D でなければなりません**。

#### `header`
`data` と組み合わせる FITS header。2D celestial header に整形されます。

#### `ext`
FITS extension name または index。`source` が file / `HDUList` のときに使います。

### 注意点

- `data` を渡した場合、`header` は必須です。
- `(header, data)` でも同じことができます。
- 非 spatial WCS keyword は落とし、celestial WCS だけを保持する方向に整形されます。

---

## 7.2 `resolve_cube_input()`

```python
resolve_cube_input(source, *, ext=None) -> CubeInput
```

3D cube 入力を `CubeInput` に正規化します。

### パラメータ

#### `source`
3D cube source。受け取れる型:

- `SpectralCube`
- FITS path
- `fits.HDUList`
- 3D HDU
- `(header, data)`

#### `ext`
FITS extension name/index。

### 注意点

- 3D 以外を渡すとエラーです。
- `BUNIT` があれば unit を拾います。

---

## 7.3 `build_normalize()`

```python
build_normalize(
    data,
    *,
    mode='asinh',
    percentile=None,
    cmin=None,
    cmax=None,
    stretch_a=0.1,
    power_gamma=1.0,
    invalid=np.nan,
    vmin=None,
    vmax=None,
) -> ImageNormalize
```

表示レンジと stretch を指定して `ImageNormalize` を返します。

### パラメータ

#### `data`
正規化対象の配列。

#### `mode`
正規化方式。

- `'linear'`
- `'sqrt'`
- `'log'`
- `'asinh'`
- `'power'`

#### `percentile`
`(low, high)` percentile。例えば `(1, 99.5)`。

- `cmin/cmax` を与えないときの自動レンジ決定に使います。
- `None` の場合は finite 値の min/max を使います。

#### `cmin`, `cmax`
表示 color scale の下限・上限。

重要:

- これは **速度範囲ではありません**。
- `vel_range` / `chan_range` とは全く別です。

#### `stretch_a`
`mode='asinh'` の `AsinhStretch(a=...)` の `a`。

#### `power_gamma`
`mode='power'` の指数。

#### `invalid`
`ImageNormalize` の `invalid` に渡す値。

#### `vmin`, `vmax`
**非推奨**。後方互換用の alias です。意味は `cmin/cmax` と同じです。

### 注意点

- `log` は正の有限値がないとエラーになります。
- `cmin == cmax` になった場合は安全のため `cmax = cmin + 1` にされます。

---

## 7.4 `apply_gaussian_smoothing_2d()`

```python
apply_gaussian_smoothing_2d(
    data,
    wcs,
    *,
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
    boundary='extend',
) -> np.ndarray
```

2D Gaussian smoothing を arcsec 単位で行います。

### パラメータ

#### `data`
2D array。

#### `wcs`
pixel scale を求めるための celestial WCS。

#### `smooth_fwhm_arcsec`
追加でかける smoothing kernel の FWHM [arcsec]。

#### `target_hpbw_arcsec`
最終的にそろえたい beam size [arcsec]。

#### `orig_hpbw_arcsec`
元の beam size [arcsec]。`target_hpbw_arcsec` と組み合わせて kernel を逆算するときに必要です。

#### `boundary`
`astropy.convolution.convolve()` に渡す boundary option。既定は `'extend'`。

### 重要なルール

- `smooth_fwhm_arcsec` と `target_hpbw_arcsec` は **同時指定不可**です。
- `target_hpbw_arcsec` を使う場合、`orig_hpbw_arcsec` が分からないと完全には意味が定まりません。
- `orig_hpbw_arcsec` が不明のときは警告を出し、`target_hpbw_arcsec` を実質 kernel FWHM とみなす実装になっています。

### beam の実効値

内部では `_effective_hpbw_arcsec()` により、次の考え方で effective beam を推定しています。

- `target_hpbw_arcsec` と `orig_hpbw_arcsec` がある → effective beam は target
- `smooth_fwhm_arcsec` と `orig_hpbw_arcsec` がある → beam を二乗和平方根で更新
- 元の beam が不明 → effective beam は不明のまま

---

## 7.5 `estimate_rms_robust()`

```python
estimate_rms_robust(data, *, finite_only=True) -> float
```

MAD ベースの robust RMS 推定です。主用途は contour level の自動決定です。

### パラメータ

#### `data`
入力配列。

#### `finite_only`
`True` のとき NaN / Inf を除外します。

### 注意点

- 精密なノイズ解析用ではなく、**描画用の実用的な RMS 推定**です。
- MAD が 0 や非有限なら標準偏差に fallback します。

---

## 7.6 `compute_contour_levels()`

```python
compute_contour_levels(
    data,
    *,
    levels=None,
    level_mode='auto',
    rms=None,
    sigma_levels=None,
    negative_sigma_levels=None,
    fraction_levels=None,
    level_scale=1.0,
    sigma_scale=1.0,
    symmetric=False,
) -> np.ndarray
```

contour level を **legacy auto / fraction / manual / rms** で作る関数です。v5 系の途中では `auto` が sigma 優先になっていた版もありましたが、**v5p5 では後方互換性を優先し、`auto` は legacy fraction に戻しています。**

### パラメータ

#### `data`
contour を引く対象の **2D 配列**です。

- `Map2D` そのものではなく、実際には `np.ndarray` を想定します。
- `NaN` を含んでいても構いません。内部では finite 値だけを評価します。
- contour の level 計算は **この data の値**に対して行われます。したがって、事前 smoothing や crop が level 計算に影響します。

#### `levels`
level の与え方です。意味は `level_mode` によって変わります。

1. **`levels='auto'`**
   - `level_mode='auto'` と組み合わせるのが最も自然です。
   - 現行版では **legacy auto = positive peak に対する fraction** を意味します。

2. **数値配列を渡す場合**
   - `level_mode='manual'` または `'explicit'` のときは、**factor 配列**として解釈します。
   - 実際の contour level は `level_scale * levels` です。
   - 例: `levels=[1,2,3,4], level_scale=3` → 実際の level は `[3,6,9,12]`

3. **非 manual モードで数値配列を渡す場合**
   - 後方互換性のため、**すでに最終値が入っている**と解釈してそのまま返します。
   - つまり `level_mode='fraction'` でも `levels=[2,5,10]` を渡した場合、それは fraction としてではなく、**最終 contour 値**として扱われます。
   - この仕様は便利な一方で、少し分かりにくいので、新規コードでは `fraction_levels` / `sigma_levels` を明示する方が安全です。

#### `level_mode`
次のいずれかです。

- `'auto'`
- `'fraction'`
- `'manual'`
- `'explicit'`
- `'rms'`
- `'sigma'`

意味は次のとおりです。

- **`'auto'`**  
  現行版では **legacy fraction mode** を意味します。正の peak を基準に、既定では `0.1, 0.3, 0.5, 0.7, 0.9` 倍の contour を作ります。

- **`'fraction'`**  
  `fraction_levels` を使って、**正の peak に対する比率**で contour を作ります。

- **`'manual'`**  
  `levels` を factor 配列とみなし、`level_scale * levels` を最終 level とします。  
  これは `3*[1,2,3,4]` 的な指定を、誤解のない形で表現するためのモードです。

- **`'explicit'`**  
  後方互換 alias です。内部では `manual` と同じ扱いになります。

- **`'rms'`**  
  `rms`, `sigma_levels`, `sigma_scale` を使って、`rms * sigma_scale * sigma_levels` で level を作ります。

- **`'sigma'`**  
  後方互換 alias です。内部では `rms` と同じ扱いになります。

#### `rms`
`level_mode='rms'` / `'sigma'` で使う RMS 値です。

- 明示値があればそれを使います。
- 省略時は `estimate_rms_robust()` による **MAD ベースの簡易推定**を使います。
- ただしこの自動推定は、**2D map 全体の値から決める便宜的 fallback** です。特に moment0 などでは、輝線領域を含むと大きめになることがあり、物理的に最適とは限りません。
- したがって、**本当に意味のある sigma contour を描きたい場合は `rms=` を明示する方がよい**です。

#### `sigma_levels`
RMS contour の倍率配列です。

- 既定は `[3, 5, 7]`
- `sigma_scale` と組み合わせて使います。
- 例: `sigma_levels=[1,2,3,4], sigma_scale=3` → `3,6,9,12 × rms`

#### `negative_sigma_levels`
負側 contour の倍率配列です。

- 明示すれば、その倍率で負側 contour を作ります。
- 省略かつ `symmetric=True` の場合は、正側の `sigma_levels` を反転して負側へ使います。
- `rms` モードでのみ意味があります。

#### `fraction_levels`
fraction contour で使う peak に対する比率です。

- 既定は `[0.1, 0.3, 0.5, 0.7, 0.9]`
- `auto` でも内部的にはこの既定値が使われます。
- 新規コードで fraction contour を明示したいときは、`levels` ではなく **こちらを使う**のが推奨です。

#### `level_scale`
manual / explicit contour の倍率です。

- `level_mode='manual'` のときにだけ意味があります。
- 例: `levels=[1,2,3,4], level_scale=3` → `[3,6,9,12]`
- 手元で基準だけ変えたいときに便利です。
- `fraction` には通常不要です。

#### `sigma_scale`
RMS contour の追加倍率です。

- `level_mode='rms'` のときに使います。
- 例: `rms=0.8, sigma_levels=[1,2,3,4], sigma_scale=3` → `[2.4, 4.8, 7.2, 9.6]`
- `sigma_levels` だけでなく、共通係数を切り替えたいときに便利です。

#### `symmetric`
`True` のとき、fraction または RMS contour に対して **負側も対称に追加**します。

- fraction: `[-0.9, -0.7, ..., +0.7, +0.9] × peak`
- rms: `[-7σ, -5σ, -3σ, +3σ, +5σ, +7σ]` のような形
- ただし `negative_sigma_levels` を明示した場合は、その指定を優先します。

### contour level の整理

現在の contour level の考え方は次の 4 系統です。

1. **`auto`**  
   旧版互換。既定の fraction。初心者向け、あるいは過去コードとの互換用途。

2. **`fraction`**  
   明示的な fraction contour。radio map の見やすい quicklook に向いています。

3. **`manual`**  
   `levels=[1,2,3,4], level_scale=3` のように、規則的な level 列を安全に書けます。

4. **`rms`**  
   `rms` が信頼できるときにだけ使う。moment0 の物理的な不確かさまで厳密に扱うものではなく、あくまで contour 指定法の一つです。

### 注意点

- **`auto` は v5p5 で fraction に戻っています。** 以前の試作版のような「まず sigma を試す」挙動ではありません。
- `rms` の自動推定は、moment0 の物理誤差を厳密に表すものではありません。
- 新規コードでは、
  - fraction → `fraction_levels=[...]`
  - manual → `levels=[...], level_scale=...`
  - rms → `rms=..., sigma_levels=[...], sigma_scale=...`
  のように **意味の違う引数を分ける**方が分かりやすいです。
- 非 manual モードで `levels=[...]` を入れると、後方互換のため **最終値扱い**になります。ここは仕様として残っていますが、分かりにくいと感じる場合は使わない方が安全です。

---
## 7.7 `crop_map2d()`

```python
crop_map2d(
    map2d,
    *,
    center=None,
    size=None,
    xlim=None,
    ylim=None,
    mode='trim',
    fill_value=np.nan,
) -> Map2D
```

2D map を切り出します。

### パラメータ

#### `map2d`
切り出し対象の `Map2D`。

#### `center`
cutout 中心。

- `(xpix, ypix)`
- `SkyCoord`

#### `size`
cutout size。

- scalar
- `(ny, nx)`
- `Quantity`（world size）

#### `xlim`, `ylim`
pixel 範囲。`(xmin, xmax)`, `(ymin, ymax)` を与えると、内部で center+size に変換します。

#### `mode`
`Cutout2D` の mode。通常は `'trim'` を使います。

#### `fill_value`
はみ出し領域の埋め値。

### ルール

- `xlim` と `ylim` は **両方必要**です。
- `center+size` 方式と `xlim+ylim` 方式のどちらかを使います。

---

## 7.8 `describe_map2d()`

```python
describe_map2d(map2d, *, name=None) -> str
```

2D map の簡単な説明文字列を返します。

出力には主に次が含まれます。

- 名前（`name` を与えた場合）
- shape
- coordinate system
- unit
- mode
- range kind / range value
- beam or effective beam
- value range（min / max / median）

### パラメータ

#### `map2d`
説明対象の `Map2D`。

#### `name`
説明に付ける名前。

---

## 7.9 `describe_source()`

```python
describe_source(
    source=None,
    *,
    data=None,
    header=None,
    ext=None,
    mode='moment0',
    chan_range=None,
    vel_range=None,
    spectral_unit='km/s',
    mask=None,
    mask_mode=None,
    zero_fill=False,
    nan_fill=True,
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
    linefree_ext='LINEFREE',
    linecand_ext='LINECAND3D',
    basesup_ext='BASESUP3D',
    final_mask_ext='MASK3D',
) -> str
```

「この source をこういう設定で 2D 化したら何が描かれるか」を短く説明します。

### パラメータ

ほぼ `make_2d_map()` と同じです。違いは **実際に図を描かず、説明テキストを返す**ことです。

#### `source`
対象 source。

#### `data`, `header`
メモリ上 2D/3D データを使う場合の組。

- `data` を使うとき `header` は必須
- `source` と同時に渡したときは `(header, data)` 側が優先されます

残りのパラメータは、後述の `make_2d_map()` と同じ意味です。

---

## 7.10 `plotting_help_text()` / `print_plotting_help()`

```python
plotting_help_text() -> str
print_plotting_help() -> None
```

コンパクトな API help を返す / 表示する関数です。

用途:

- notebook や対話環境で「簡単な使い方」を確認したいとき
- `quicklook(help=True)` の補助
- CLI の `--api-help`

---

## 7.11 `make_2d_map()`

```python
make_2d_map(
    source,
    *,
    ext=None,
    mode='moment0',
    chan_range=None,
    vel_range=None,
    spectral_unit='km/s',
    mask=None,
    mask_mode=None,
    zero_fill=False,
    nan_fill=True,
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
    fill_value=np.nan,
    linefree_ext='LINEFREE',
    linecand_ext='LINECAND3D',
    basesup_ext='BASESUP3D',
    final_mask_ext='MASK3D',
) -> Map2D
```

このモジュールの中核関数です。2D 入力なら passthrough、3D 入力なら 2D map を作ります。

### `mode`
サポートされる mode:

- `'identity'`
- `'map'`
- `'2d'`
- `'moment0'`
- `'channel_sum'`
- `'channel_mean'`
- `'channel_slice'`
- `'velocity_slice'`
- `'provisional_moment'`
- `'final_moment'`

### パラメータ

#### `source`
2D/3D source。

#### `ext`
対象 HDU。

#### `mode`
2D 化の方法。

意味:

- `identity` / `map` / `2d`  : 2D としてそのまま通す
- `moment0`                  : スペクトル軸方向に積分
- `channel_sum`              : channel 値を単純加算
- `channel_mean`             : channel 平均
- `channel_slice`            : 1 channel だけ抜き出す
- `velocity_slice`           : 指定速度に最も近い 1 面を抜き出す
- `provisional_moment`       : provisional mask を用いる moment
- `final_moment`             : `MASK3D` を使う final signal moment

#### `chan_range`
`(lo, hi)` のチャネル範囲。両端を含みます。

- `chan_range=(10, 20)` なら 10 〜 20 ch を使います
- `hi < lo` でも入れ替えて処理します
- channel index は 0-based です

#### `vel_range`
`(vlo, vhi)` の速度範囲。単位は `spectral_unit`。

- `hi < lo` でも入れ替えます
- `moment0` や `channel_sum` などで使用できます
- `velocity_slice` では `(v1, v2)` の平均値に最も近い 1 channel を取ります

#### `spectral_unit`
速度軸演算に使う単位。

- 既定は `'km/s'`
- 内部では `with_spectral_unit(..., velocity_convention='radio')` を試します
- 変換できない場合は warning 相当の情報を meta に残し、native axis fallback する設計です

#### `mask`
外部マスク配列。`mask_mode='external'` 相当で使います。

#### `mask_mode`
3D reduction 時のマスク指定。

サポート値:

- `'external'`
- `'linefree_complement'`
- `'linecand3d'`
- `'basesup3d'`
- `'mask3d'`
- `None`

意味:

- `external`             : `mask=` を使う
- `linefree_complement`  : `LINEFREE` の補集合
- `linecand3d`           : `LINECAND3D == True`
- `basesup3d`            : `BASESUP3D == False`
- `mask3d`               : `MASK3D == True`

#### `zero_fill`
mask 外を 0 で埋める。

#### `nan_fill`
mask 外を NaN で埋める。既定はこちらです。

#### `smooth_fwhm_arcsec`
追加 smoothing の FWHM。

#### `target_hpbw_arcsec`
目標 beam size。

#### `orig_hpbw_arcsec`
元 beam size。`target_hpbw_arcsec` と併用したいときに重要です。

#### `fill_value`
`zero_fill=False` かつ `nan_fill=False` のときの mask 外埋め値。

#### `linefree_ext`, `linecand_ext`, `basesup_ext`, `final_mask_ext`
associated mask extension 名。既定はそれぞれ:

- `LINEFREE`
- `LINECAND3D`
- `BASESUP3D`
- `MASK3D`

### 補足

- `mode='moment0'` では channel 幅（速度幅）を掛けて積分します。
- 出力 unit は概ね `data_unit * spectral_unit` になります。
- 2D passthrough のときは、3D reduction 引数は無視されます。

---

## 7.12 `make_provisional_moment()`

```python
make_provisional_moment(
    source,
    *,
    ext=None,
    linefree_ext='LINEFREE',
    basesup_ext='BASESUP3D',
    linecand_ext='LINECAND3D',
    prefer='auto',
    vel_range=None,
    chan_range=None,
    spectral_unit='km/s',
    zero_fill=False,
    nan_fill=True,
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
) -> Map2D
```

baseline / signal 補助 extension を使って provisional moment を作ります。

### `prefer`

- `'auto'`
- `'linecand3d'`
- `'basesup3d'`
- `'linefree_complement'`

`'auto'` の探索順:

1. `LINECAND3D`
2. `BASESUP3D`
3. `LINEFREE`

---

## 7.13 `make_final_moment()`

```python
make_final_moment(
    source,
    *,
    ext=None,
    final_mask_ext='MASK3D',
    vel_range=None,
    chan_range=None,
    spectral_unit='km/s',
    zero_fill=False,
    nan_fill=True,
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
) -> Map2D
```

`MASK3D` を使う final signal moment 作成関数です。

ほぼ `make_2d_map(..., mode='final_moment', mask_mode='mask3d')` のラッパーです。

---

## 7.14 `add_contours()`

```python
add_contours(ax, contours, *, base_map=None) -> List[Any]
```

既存の WCSAxes に contour layer を追加します。`plot_map()` / `plot_scene()` / `quicklook(contours=...)` の contour 描画は、最終的にこの関数を通ります。

### パラメータ

#### `ax`
WCSAxes。

#### `contours`
contour layer spec の列です。各要素は dict を想定します。

#### `base_map`
contour layer に `source` / `data` / `header` が一切ない場合の fallback `Map2D` です。

- `plot_map(map2d, contours=[...])` のようなケースでは、ここに base の `Map2D` が渡されます。
- `quicklook(contours=[...])` では、基本的に quicklook の base と同じ source / mode / vel_range などを引き継いだ contour layer が内部生成されます。

### contour layer で使える主なキー

#### 2D/3D map 生成系

- `source`
- `data`
- `header`
- `ext`
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
- `crop`

補足:

- これらのキーが含まれると、内部で `make_2d_map()` を呼んで contour 用 `Map2D` を作ります。
- `smooth_fwhm_arcsec` / `target_hpbw_arcsec` を `make_2d_map()` 側で既に使った場合、現在は **二重 smoothing を避ける**よう整理されています。
- `crop` は contour 側にも正しく反映されます。以前の版では、base には crop が効いても contour 側には level 計算時の crop が効かないことがありました。

#### contour style 系

- `levels`
- `level_mode`
- `rms`
- `sigma_levels`
- `negative_sigma_levels`
- `fraction_levels`
- `level_scale`
- `sigma_scale`
- `symmetric`
- `colors`
- `linewidths`
- `linestyles`
- `alpha`
- `label`

### contour level の考え方

現在は、以下の使い分けを推奨します。

- **旧版互換・迷ったらまずこれ**  
  `{"levels": "auto"}`

- **fraction を明示したい**  
  `{"level_mode": "fraction", "fraction_levels": [0.1,0.3,0.5,0.7,0.9]}`

- **manual による規則列**  
  `{"level_mode": "manual", "levels": [1,2,3,4], "level_scale": 3}`

- **rms を明示して σ contour**  
  `{"level_mode": "rms", "rms": 0.8, "sigma_levels": [1,2,3,4], "sigma_scale": 3}`

### 現在の仕様で特に重要な点

1. **`levels='auto'` は legacy fraction** です。  
   旧版と同じ見た目をできるだけ保つためです。

2. **`quicklook(contours=...)` が正式対応**しました。  
   contour layer に `source` や `mode` を省略した場合は、quicklook の base 設定を継承します。

3. **同じ smoothing / target HPBW が二重適用されない**ように整理しました。

4. **contour 側の `crop` も level 計算に反映**されます。

### 注意点

- contour は `ax.contour(..., transform=ax.get_transform(cmap.wcs))` で描きます。したがって contour 用 map の WCS が正しいことが重要です。
- `rms` モードは使えますが、moment0 では物理的な誤差伝播を厳密に行っているわけではありません。
- 本当に 3σ contour を物理的に言いたい場合は、3D から導いた RMS または error map を別途考える方が安全です。

---
## 7.15 `add_scalebar()`

```python
add_scalebar(ax, map2d, config=None) -> Dict[str, Any]
```

angular scalebar を描きます。

### `config` で使えるキー

#### `length`
長さ。

#### `units`
長さの単位。

- `'arcsec'`
- `'arcmin'`
- `'deg'`
- `'pix'`

#### `location`
配置位置。例:

- `'lower right'`
- `'lower left'`
- `'upper right'`
- `'upper left'`

#### `margin_frac`
余白の相対量。

#### `line_style`
`ax.plot()` に渡す style dict。

#### `color`
線と文字の基本色。

#### `linewidth`
線幅。

#### `solid_capstyle`
scalebar 線の端のスタイル。

#### `label`
表示文字列。`None` または `'auto'` なら自動表記。

#### `text_style`
文字 style dict。

#### `fontsize`
文字サイズ。

#### `text_va`
文字の縦位置。

#### `text_dy`
文字の上下オフセット（pixel 単位相当）。

### 戻り値

- `line`
- `text`
- `length_arcsec`

---

## 7.16 `add_north_arrow()`

```python
add_north_arrow(ax, map2d, config=None) -> Dict[str, Any]
```

north / east の向きを矢印で描きます。

### `config` で使えるキー

#### `length`
矢印長。

#### `units`
- `'arcsec'`
- `'arcmin'`
- `'deg'`
- `'pix'`

#### `location`
配置位置。既定は `'upper left'`。

#### `margin_frac`
余白比率。

#### `arrowprops`
`annotate()` の `arrowprops` dict。

#### `color`
矢印と文字の色。

#### `linewidth`
矢印線幅。

#### `arrowstyle`
矢印スタイル。既定 `'-|>'`。

#### `text_style`
文字 style dict。

#### `fontsize`
文字サイズ。

#### `north_label`
既定 `'N'`。

#### `east_label`
既定 `'E'`。

### 戻り値

- `north_arrow`
- `east_arrow`
- `north_text`
- `east_text`
- `length_arcsec`

---

## 7.17 `plot_map()`

```python
plot_map(
    source=None,
    *,
    data=None,
    header=None,
    ext=None,
    ax=None,
    projection=None,
    cmap='viridis',
    norm=None,
    norm_mode='asinh',
    norm_percentile=None,
    cmin=None,
    cmax=None,
    stretch_a=0.1,
    power_gamma=1.0,
    colorbar=True,
    colorbar_label=None,
    title=None,
    origin='lower',
    contours=None,
    grid=True,
    beam=None,
    scalebar=None,
    north_arrow=None,
    legend=False,
    xlabel=None,
    ylabel=None,
    crop=None,
    show_readout=True,
    describe=False,
    save=None,
    save_dpi=250,
    save_transparent=False,
    show=True,
    figsize=(10, 8),
    interpolation='nearest',
    alpha=1.0,
    vmin=None,
    vmax=None,
) -> Dict[str, Any]
```

2D map 描画の基本関数です。**ただし v5p5 では、新規コードの主 API としては `plot_scene()` または `quicklook()` を推奨**します。`plot_map()` は非常に重要な関数であり、今後もしばらく残すべきですが、位置づけとしては **互換 API** です。

### なぜ `plot_map()` を残すのか

- 過去コードとの互換性が高い
- `make_2d_map()` + `plot_map()` という流れが分かりやすい
- 単一 image + optional contour では依然として素直に書ける

### なぜ新規コードでは第一候補にしないのか

- 将来の layer 拡張は `plot_scene()` に寄せたい
- `quicklook()` で十分な用途が多い
- catalog / marker / image overlay / annotation を増やすなら base + overlays の方が自然

### 入力パラメータ

#### `source`, `data`, `header`, `ext`
`resolve_map_input()` と同じ意味です。

#### `ax`
既存 axes。省略時は新規 figure/axes を作ります。

#### `projection`
axes 作成時に使う投影 WCS。省略時は `map2d.wcs`。

### 表示パラメータ

#### `cmap`
matplotlib colormap 名。

#### `norm`
外部で作った `ImageNormalize`。これを与えると `norm_mode` などは実質使われません。

#### `norm_mode`
`build_normalize()` に渡す mode。

#### `norm_percentile`
percentile clip。

#### `cmin`, `cmax`
表示レンジ。

#### `stretch_a`, `power_gamma`
`asinh` / `power` stretch 用。

#### `colorbar`
colorbar を描くか。

#### `colorbar_label`
colorbar ラベル。`None` なら unit や mode から自動決定。

#### `title`
タイトル。

#### `origin`
`imshow()` の origin。通常 `'lower'`。

#### `interpolation`
`imshow()` の補間。既定 `'nearest'`。

#### `alpha`
画像の透明度。

#### `vmin`, `vmax`
**非推奨**。`cmin/cmax` の alias。新規コードでは使わない方がよいです。

### 重ね描き・annotation

#### `contours`
contour layer spec の列。

補足:

- `plot_map()` に `contours=` を与えると、内部で `add_contours()` を呼びます。
- 単一 image に対する contour 追加には便利ですが、今後 layer を増やすなら `plot_scene()` の方が一貫しています。

#### `grid`
座標グリッド ON/OFF。

#### `beam`
beam 表示設定。

- `'auto'` / `'header'` : header または smooth 情報から推定
- dict : beam を明示指定

beam dict の主なキー:

- `major`
- `minor`
- `angle`
- `units` (`arcsec`, `deg`, `pix`)
- `x`, `y`
- `facecolor`
- `edgecolor`
- `linewidth`
- `alpha`
- `label`

#### `scalebar`
scalebar 設定。`True` でも可。

#### `north_arrow`
north/east 矢印設定。`True` でも可。

#### `legend`
legend ON/OFF または `ax.legend()` 用 kwargs dict。

#### `xlabel`, `ylabel`
軸ラベル。省略時は WCS から自動。

### その他

#### `crop`
`crop_map2d()` に渡す dict。

#### `show_readout`
カーソル位置に応じた表示を `ax.format_coord` に設定します。

#### `describe`
`describe_map2d()` による短い説明を print します。

#### `save`
保存 path。

#### `save_dpi`
保存 DPI。

#### `save_transparent`
透明背景保存。

#### `show`
`True` なら figure を表示します。`ax` を外から渡した場合には、現在は **勝手に `plt.show()` しない**ように整理されています。

#### `figsize`
新規 figure 作成時のサイズ。

### 戻り値

- `fig`
- `ax`
- `image_artist`
- `colorbar`
- `map2d`
- `beam_artists`
- `contours`
- `scalebar`
- `north_arrow`
- `legend`
- `description`

### 実務上の推奨

- 既存コードを崩したくない場合は `plot_map()` をそのまま使ってよいです。
- 新しく書く notebook / script では、
  - 単純表示: `quicklook()`
  - 多層重ね描き: `plot_scene()`
  を第一候補にする方が、今後の保守に向いています。

---
## 7.18 `make_rgb_map()`

```python
make_rgb_map(red, green, blue, *, red_norm=None, green_norm=None, blue_norm=None) -> RGBMap
```

3 枚の 2D map を RGB 合成します。

### パラメータ

#### `red`, `green`, `blue`
各チャネル source。

受けられるもの:

- `Map2D`
- bare source
- map specification dict

#### `red_norm`, `green_norm`, `blue_norm`
各 channel の normalization dict。

使える主なキー:

- `norm_mode`
- `percentile`
- `cmin`
- `cmax`
- `stretch_a`
- `power_gamma`
- `vmin`
- `vmax`

### 制約

- 3 channel は **同じ 2D shape** が必要です。
- WCS header が一致しない場合、warning は出ますが **そのまま stack** します。
- v5 系では RGB の本格 reprojection は未対応です。

---

## 7.19 `plot_rgb()`

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
    crop=None,
    show_readout=False,
    scalebar=None,
    north_arrow=None,
    legend=False,
    describe=False,
    save=None,
    save_dpi=250,
    save_transparent=False,
    show=True,
    figsize=(10, 8),
    **rgb_kwargs,
) -> Dict[str, Any]
```

RGB 表示関数です。

### パラメータ

#### `rgb_source`
既に作成済みの `RGBMap`。

#### `red`, `green`, `blue`
`rgb_source` を与えない場合に使う 3 チャネル。

#### `crop`
各チャネルに同じ cutout を適用します。

#### `show_readout`
平均 RGB 強度を使った擬似 `Map2D` を介して readout を設定します。

#### `scalebar`, `north_arrow`, `legend`, `describe`, `save`, `show`, `figsize`
`plot_map()` と同様の意味です。

#### `**rgb_kwargs`
`make_rgb_map()` にそのまま渡されます。

### 戻り値

- `fig`
- `ax`
- `image`
- `rgb`
- `legend`
- `description`

---

## 7.20 `plot_scene()`

```python
plot_scene(
    base,
    *,
    overlays=None,
    ax=None,
    title=None,
    grid=True,
    xlabel=None,
    ylabel=None,
    colorbar=True,
    show_readout=True,
    describe=False,
    save=None,
    save_dpi=250,
    save_transparent=False,
    show=True,
    figsize=(10, 8),
    legend=False,
    crop=None,
) -> Dict[str, Any]
```

高水準 API です。**v5p5 では、新規コードの主 API は基本的にこの `plot_scene()`** だと考えてよいです。発想は **1 つの base layer + 複数 overlay** です。

### `base`
base は bare source でも dict でもよいです。

- bare source は `{'kind': 'image', 'source': source}` と解釈
- サポートされる `base kind` は:
  - `'image'`
  - `'rgb'`

### `overlays`
overlay spec の列。サポートされる `kind` は:

- `'image'`
- `'contour'`
- `'marker'`
- `'catalog'`
- `'text'`
- `'beam'`
- `'scalebar'`
- `'north_arrow'`
- `'compass'`（`north_arrow` の alias）

### top-level パラメータ

#### `ax`
既存 axes。

- `None` なら新規 figure を作成します。
- `ax` を外から渡した場合、現在は **`show=True` でも勝手に `plt.show()` しません。** これは `plot_map()` / `plot_rgb()` と揃えた挙動です。

#### `title`, `grid`, `xlabel`, `ylabel`
base 表示の共通制御。

#### `colorbar`
base が image のときの colorbar ON/OFF。

#### `show_readout`
base map の readout 設定。

#### `describe`
base の説明を表示。

#### `save`, `save_dpi`, `save_transparent`, `show`, `figsize`, `legend`, `crop`
図全体の保存と表示設定。

### 戻り値

- `fig`
- `ax`
- `base`
- `base_description`
- `overlay_artists`
- `overlay_maps`
- `contours`
- `beam_artists`
- `legend`

### `plot_scene()` を使うべき場面

- optical / IR 画像の上に radio contour を重ねたい
- 複数 line の contour を重ねたい
- marker / catalog / text / beam / scalebar / north_arrow を同時に使いたい
- 今後も拡張しやすい書き方にしたい

### `plot_map()` との関係

- `plot_map()` は重要な互換 API です
- しかし、新機能の置き場は基本的に `plot_scene()` に集約する方針です
- 実際、`quicklook()` も内部では `plot_scene()` を呼びます

このため、**将来を見据えたコードは `plot_scene()` 中心に考える**のがよいです。

---
## 7.21 `quicklook()`

```python
quicklook(
    source,
    *,
    ext=None,
    mode='moment0',
    chan_range=None,
    vel_range=None,
    spectral_unit='km/s',
    smooth_fwhm_arcsec=None,
    target_hpbw_arcsec=None,
    orig_hpbw_arcsec=None,
    cmap='viridis',
    norm_mode='asinh',
    norm_percentile=(1.0, 99.5),
    cmin=None,
    cmax=None,
    contours=None,
    title=None,
    grid=True,
    beam=None,
    crop=None,
    describe=False,
    help=False,
    save=None,
    save_dpi=250,
    save_transparent=False,
    show_readout=True,
    show=True,
    figsize=(10, 8),
) -> Dict[str, Any]
```

最小指定で 1 枚出すための入口です。**v5p5 では `contours=` に正式対応**しており、単純な image + contour を非常に短く書けるようになっています。

### パラメータ

#### `source`
描画対象 source。

#### `ext`, `mode`, `chan_range`, `vel_range`, `spectral_unit`
`make_2d_map()` の基本指定。

#### `smooth_fwhm_arcsec`, `target_hpbw_arcsec`, `orig_hpbw_arcsec`
beam / smoothing 関係。

#### `cmap`, `norm_mode`, `norm_percentile`, `cmin`, `cmax`
見た目。

#### `contours`
重要な引数です。contour layer spec の列を渡します。

基本挙動:

- `contours=None` のときは単純 image 表示
- `contours=[{...}]` を与えると、その内容を overlay contour に変換して `plot_scene()` に渡します
- contour spec 側で `source`, `mode`, `vel_range`, `chan_range`, `spectral_unit`, `smooth_fwhm_arcsec`, `target_hpbw_arcsec`, `orig_hpbw_arcsec`, `crop` を省略した場合は、**quicklook の base と同じ設定を自動継承**します
- contour spec 側でそれらを明示すれば、**別ソース・別積分範囲・別 smoothing** の contour も可能です

簡単な例:

```python
quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(0, 20),
    contours=[
        {
            'levels': 'auto',
            'colors': 'white',
            'linewidths': 0.8,
            'alpha': 0.6,
        }
    ],
)
```

#### `title`, `grid`, `beam`, `crop`
図の補助設定。

#### `describe`
説明を print。

#### `help`
`print_plotting_help()` を実行してから続行します。

#### `save`, `save_dpi`, `save_transparent`, `show_readout`, `show`, `figsize`
出力設定。

### 位置づけ

- 初心者向け入口
- notebook でまず 1 枚見る入口
- `plot_scene()` に進む前の quick inspection
- **新規コードで最初に試す第一候補**

### `plot_map()` との違い

- `quicklook()` は **簡便さ重視**
- `plot_map()` は **互換性重視**

新規コードでは、

- まず `quicklook()`
- 複雑になったら `plot_scene()`
- 過去コード互換が必要なら `plot_map()`

という使い分けが自然です。

---
## 7.22 `main()` と CLI

```bash
python plotting_updated_v5p5.py [options] source
```

CLI は、`plot_map()` を使った単純表示を command line からすぐ実行するための入口です。現状の CLI は **単一 image + optional contour** に強く、`plot_scene()` の全機能をそのまま expose しているわけではありません。

### CLI オプション完全一覧

#### `--api-help`
compact Python API help を表示して終了します。`source` が与えられていれば、そのまま処理続行も可能です。

#### `source`
入力 FITS file。

#### `--ext`
HDU index または extension name。

#### `--mode`
2D 化モード。選択肢:

- `identity`
- `moment0`
- `channel_sum`
- `channel_mean`
- `channel_slice`
- `velocity_slice`
- `provisional_moment`
- `final_moment`

#### `--chan-range`
`lo,hi`

#### `--vel-range`
`lo,hi`（`--spectral-unit` に従う）

#### `--spectral-unit`
既定 `km/s`

#### `--smooth-fwhm-arcsec`
追加 smoothing FWHM

#### `--target-hpbw-arcsec`
目標 HPBW

#### `--orig-hpbw-arcsec`
元 HPBW

#### `--cmap`
colormap 名

#### `--norm-mode`
`linear|sqrt|log|asinh|power`

#### `--norm-percentile`
`lo,hi`。既定は `1,99.5`

#### `--cmin`, `--cmax`
表示レンジ

#### `--vmin`, `--vmax`
非推奨 alias

#### `--stretch-a`
`asinh` 用

#### `--power-gamma`
`power` 用

#### `--title`
タイトル

#### `--describe`
短い説明を print

#### `--beam`
beam 設定。現状の help 文は `auto/header or omit`。CLI では主に `auto` / `header` を想定します。

#### `--legend`
legend 表示

#### `--scalebar`
scalebar 表示

#### `--scalebar-length`
scalebar 長さ [arcsec]

#### `--north-arrow`
north/east 矢印表示

#### `--no-grid`
座標グリッドを切る

#### `--output`
保存 path

#### `--dpi`
保存 DPI

#### `--transparent`
透明背景保存

#### `--contour-levels`
`v1,v2,v3` または `'auto'`

重要:

- `manual` モードでは、ここで与えた値は **最終値そのものではなく factor** として使われ、`--contour-level-scale` が掛かります
- 例: `--contour-level-mode manual --contour-levels 1,2,3,4 --contour-level-scale 3`
  → 実際の contour は `3,6,9,12`

#### `--contour-level-mode`
`auto|fraction|manual|explicit|rms|sigma`

現在の意味:

- `auto` : legacy fraction
- `fraction` : `--contour-fractions` を使う
- `manual` / `explicit` : `--contour-levels` と `--contour-level-scale` を使う
- `rms` / `sigma` : `--contour-rms`, `--contour-sigmas`, `--contour-sigma-scale` を使う

#### `--contour-level-scale`
manual contour 用の倍率

例:

```bash
--contour-level-mode manual --contour-levels 1,2,3,4 --contour-level-scale 3
```

→ 実際の level は `3,6,9,12`

#### `--contour-fractions`
`f1,f2,f3`。fraction contour 用。

例:

```bash
--contour-level-mode fraction --contour-fractions 0.05,0.1,0.2,0.4,0.8
```

#### `--contour-sigmas`
`s1,s2,s3`。rms contour の倍率。

#### `--contour-sigma-scale`
rms contour の追加倍率。

例:

```bash
--contour-level-mode rms --contour-rms 0.8 --contour-sigmas 1,2,3,4 --contour-sigma-scale 3
```

→ 実際の level は `3,6,9,12 × 0.8`

#### `--contour-negative`
負側 contour も追加

- fraction モードでは `symmetric=True` に相当
- rms モードでも同様

#### `--contour-rms`
明示 RMS

#### `--contour-color`
contour 色

#### `--contour-linewidth`
contour 線幅

#### `--zero-fill`
mask 外を 0 に

#### `--no-nan-fill`
mask 外 NaN fill を無効化

### CLI contour 仕様の注意

以前の版では、`--contour-levels` を与えないと contour 自体が生成されないことがありました。現在は **contour 関連の任意オプションが 1 つでも指定されれば contour を作る**ようになっています。

つまり、たとえば次のような指定も有効です。

```bash
python plotting_updated_v5p5.py cube.fits   --mode moment0   --vel-range 0,20   --contour-level-mode fraction   --contour-fractions 0.1,0.3,0.5,0.7,0.9
```

### CLI の位置づけ

- 単純な quick inspection には有用
- ただし、本格的な multi-layer scene を全て CLI で完結させるより、Python API の `quicklook()` / `plot_scene()` を使う方が柔軟です
- したがって CLI は、**軽量な入口**として考えるのが自然です

---
# 8. `plot_scene()` の layer specification 完全リファレンス

ここが最も重要です。`plot_scene()` では、base と overlay を辞書で指定できます。

## 8.1 共通原則

- bare object は `{ 'source': obj }` 風に解釈されることがあります
- `kind` が重要です
- 2D/3D の reduction パラメータを含む layer は `_layer_to_map2d()` によって `make_2d_map()` 経由で 2D 化されます

---

## 8.2 base layer: `kind='image'`

### 代表例

```python
base = {
    'kind': 'image',
    'source': 'cube.fits',
    'mode': 'moment0',
    'vel_range': (20, 40),
    'cmap': 'inferno',
}
```

### 使える主なキー

#### source 解決系

- `source`
- `data`
- `header`
- `ext`

#### 2D/3D reduction 系

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

#### 表示系

- `projection`
- `cmap`
- `norm`
- `norm_mode`
- `norm_percentile`
- `cmin`
- `cmax`
- `stretch_a`
- `power_gamma`
- `colorbar`
- `colorbar_label`
- `title`
- `origin`
- `grid`
- `beam`
- `xlabel`
- `ylabel`
- `crop`
- `describe`
- `interpolation`
- `alpha`
- `scalebar`
- `north_arrow`

---

## 8.3 base layer: `kind='rgb'`

### 代表例

```python
base = {
    'kind': 'rgb',
    'red': {'source': '12co.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
    'green': {'source': '13co.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
    'blue': {'source': 'c18o.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
}
```

### 使える主なキー

- `rgb_source`
- `red`
- `green`
- `blue`
- `red_norm`
- `green_norm`
- `blue_norm`
- `title`
- `grid`
- `xlabel`
- `ylabel`
- `crop`
- `describe`
- `scalebar`
- `north_arrow`

---

## 8.4 overlay: `kind='contour'`

### 代表例

```python
overlay = {
    'kind': 'contour',
    'source': '13co.fits',
    'mode': 'moment0',
    'vel_range': (20, 40),
    'level_mode': 'sigma',
    'sigma_levels': [3, 5, 7],
    'colors': 'cyan',
    'linewidths': 1.0,
    'label': '13CO',
}
```

### 使えるキー

#### map 生成系
`make_2d_map()` 相当のキーをほぼすべて使えます。

#### contour style 系

- `levels`
- `level_mode`
- `rms`
- `sigma_levels`
- `negative_sigma_levels`
- `fraction_levels`
- `symmetric`
- `colors`
- `linewidths`
- `linestyles`
- `alpha`
- `label`

### 注意点

- contour は **WCSAxes transform** で描きます。
- contour 自体は image のような reprojection を標準では行いません。

---

## 8.5 overlay: `kind='image'`

### 代表例

```python
overlay = {
    'kind': 'image',
    'source': 'radio_moment0.fits',
    'cmap': 'magma',
    'alpha': 0.45,
    'reproject': 'auto',
    'label': 'radio image',
}
```

### 使えるキー

- image base とほぼ同様の map 生成・表示キー
- さらに:
  - `reproject` : `'never' | 'auto' | 'required'`
  - `label`
  - `zorder`

### `reproject` の意味

- `never`    : 再投影しない
- `auto`     : 必要なら再投影を試す。失敗したら warning でスキップ
- `required` : 必要なら再投影。失敗したらエラー

---

## 8.6 overlay: `kind='marker'`

点群を scatter します。

### 座標の与え方

次のいずれか 1 つを使います。

- `coords=[[lon, lat], ...]` + `frame`
- `skycoord=SkyCoord(...)`
- `xy=[[x, y], ...]`
- `lon=[...], lat=[...]` + `frame`
- `table=...` / `catalog=...`（ただし label なしの点群用途なら marker 的にも使えます）

### style 系

- `style` : `ax.scatter()` に渡す dict
- `s`
- `marker`
- `color`
- `linewidths`
- `label`

### 座標系キー

- `frame`
  - `'icrs'`, `'fk5'`, `'fk4'`, `'galactic'`, `'equatorial'`, `'j2000'`, `'radec'`, `'l/b'` など
- `coord_unit`
  - world 値の単位。既定 `'deg'`

---

## 8.7 overlay: `kind='catalog'`

marker と似ていますが、ラベル付与に向いています。

### 位置の与え方

`marker` と同様です。

### table-like 入力

受けられる table-like の代表:

- `pandas.DataFrame`
- `astropy.table.Table`
- `dict of arrays`

位置列の自動探索候補:

#### world 座標列

lon 候補:

- `lon`
- `glon`
- `l`
- `longitude`
- `ra`
- `ra_deg`
- `lon_deg`
- `glon_deg`
- `l_deg`

lat 候補:

- `lat`
- `glat`
- `b`
- `latitude`
- `dec`
- `dec_deg`
- `lat_deg`
- `glat_deg`
- `b_deg`

#### pixel 座標列

x 候補:

- `x`
- `xpix`
- `x_pix`
- `xpixel`
- `col`
- `ix`

y 候補:

- `y`
- `ypix`
- `y_pix`
- `ypixel`
- `row`
- `iy`

### table 用の追加キー

- `table`
- `catalog`
- `coord_mode` : `'auto'|'world'` など
- `lon_col`, `lat_col`
- `x_col`, `y_col`
- `label_col`

### label 系

- `labels` : 明示 label list
- `label_col` : table から読む列名
- `text_style`
- `text_color`
- `fontsize`
- `ha`, `va`
- `label_dx`, `label_dy`
- `label_offset_mode` : `'pixel'` / `'world'`

### label の自動探索候補

`labels` が無い場合、table から次の列を探します。

- `label`
- `labels`
- `name`
- `names`
- `id`
- `ids`
- `source`
- `object`
- `obj`

### 注意

- `label_offset_mode='pixel'` が既定です。通常はこちらが安全です。
- `world` offset は座標値に直接足すので、慎重に使ってください。

---

## 8.8 overlay: `kind='text'`

任意位置に文字を置きます。

### 位置の与え方

- `coord=(lon, lat)` + `frame`
- `skycoord=SkyCoord(...)`
- `xy=(x, y)`

### style 系

- `text`
- `style`
- `color`
- `fontsize`
- `ha`
- `va`
- `frame`
- `coord_unit`

---

## 8.9 overlay: `kind='beam'`

### 使い方

```python
{'kind': 'beam', 'beam': 'auto'}
```

または

```python
{'kind': 'beam', 'beam': {'major': 20, 'minor': 15, 'angle': 30, 'units': 'arcsec'}}
```

### キー

- `beam`
- `config`
- `label`

`beam` / `config` の中身は `plot_map(beam=...)` と同じです。

---

## 8.10 overlay: `kind='scalebar'`

```python
{'kind': 'scalebar', 'config': {'length': 30, 'units': 'arcsec'}}
```

`config` を省くと layer 自体の dict をそのまま設定とみなします。

---

## 8.11 overlay: `kind='north_arrow'` / `'compass'`

```python
{'kind': 'north_arrow', 'config': {'location': 'upper left'}}
```

---

# 9. 座標系の扱い

## 9.1 base WCS と overlay WCS

- base image は自分の WCS を持つ
- contour は `ax.get_transform(cmap.wcs)` により、その WCS から重ねる
- marker / catalog / text は world 値を **base の座標系へ変換した後**、`ax.get_transform('world')` で描く

## 9.2 frame の正規化

受け付ける alias の例:

- `equatorial` → `icrs`
- `j2000` → `icrs`
- `radec`, `ra/dec` → `icrs`
- `gal`, `glon/glat`, `l/b` → `galactic`

## 9.3 galactic / equatorial 混在

catalog や marker で galactic を与え、base が RA/Dec でも構いません。内部で frame 変換します。

---

# 10. 見た目のベストプラクティス

## 10.1 まずは `quicklook()`

最初の 1 枚は `quicklook()` が最も楽です。

## 10.2 電波同士の比較

- まず contour 重ね
- 必要なら `target_hpbw_arcsec` で beam をそろえる

## 10.3 光学 + 電波

- base = 光学画像
- overlay = radio contour
- radio pseudo-color overlay は必要時のみ

## 10.4 `cmin/cmax` と `vel_range` を混同しない

- `vel_range` はデータの積分範囲
- `cmin/cmax` は表示レンジ

---

# 11. 制限と注意点

1. **RGB は shape 一致前提** です。  
2. contour は標準では reprojection しません。  
3. `reproject` は image overlay 専用の optional 機能です。  
4. `target_hpbw_arcsec` を厳密に使いたいなら、`orig_hpbw_arcsec` が分かっている方が良いです。  
5. `log` normalization は正値が必要です。  
6. `show_readout` は generic な `lon/lat` 表示です。RA/Dec の sexagesimal 表示まではしていません。  
7. 2D extension を `make_2d_map()` に渡した場合は 2D passthrough になるため、`vel_range` などは効きません。  

---

# 12. どの関数から使い始めるべきか

- **初心者**: `quicklook()`
- **2D map を丁寧に描く**: `plot_map()`
- **RGB**: `make_rgb_map()` + `plot_rgb()`
- **複数レイヤー**: `plot_scene()`
- **スクリプト前に仕様確認**: `describe_source()` と `plotting_help_text()`

---



# 12A. 現在の仕様で特に重要な点まとめ

この節は、過去版から移行する人が見落としやすい点をまとめたものです。

## 12A.1 `quicklook(contours=...)` の追加

v5p5 では、`quicklook()` は contour を受けられます。これにより、

- 画像 1 枚だけをまず見る
- 同じ source / mode / vel_range で contour を軽く重ねる

という用途では、`plot_map()` を使わなくてもかなり短く書けます。

## 12A.2 `auto` contour の意味

`levels='auto'` は **legacy fraction** です。既定では

- 0.1
- 0.3
- 0.5
- 0.7
- 0.9

の正の peak 比で contour を作ります。

## 12A.3 `manual` contour の考え方

`manual` では、`levels` は **factor 配列**です。これに `level_scale` を掛けて最終 level を作ります。

```python
levels=[1,2,3,4], level_scale=3
```

→ `3,6,9,12`

これは、`3*[1,2,3,4]` のような意図を、誤解のない形で実現するための設計です。

## 12A.4 `rms` contour の考え方

`rms` では、`sigma_levels` に `sigma_scale` と `rms` を掛けて level を作ります。

```python
rms=0.8, sigma_levels=[1,2,3,4], sigma_scale=3
```

→ `2.4, 4.8, 7.2, 9.6`

## 12A.5 `plot_map()` の位置づけ

`plot_map()` は重要な API ですが、将来の中心は `plot_scene()` です。したがって:

- 過去コードの維持 → `plot_map()`
- 新しく簡単に始める → `quicklook()`
- 多層表示・将来拡張 → `plot_scene()`

という考え方が自然です。

---

# 12B. 代表的な移行パターン

## 12B.1 旧来の `make_2d_map() + plot_map()` をそのまま維持する

既存 notebook や script を大きく変えたくない場合は、この形で構いません。

```python
map2d = make_2d_map(
    out_fits,
    mode='moment0',
    vel_range=(0.0, 20.0),
    spectral_unit='km/s',
    target_hpbw_arcsec=450.0,
)

result = plot_map(
    map2d,
    cmap='turbo',
    norm_mode='asinh',
    norm_percentile=(1.0, 99.5),
    stretch_a=0.10,
    contours=[
        {
            'levels': 'auto',
            'colors': 'white',
            'linewidths': 0.8,
            'alpha': 0.6,
        }
    ],
    beam='header',
    title='Moment 0 with target HPBW',
    show=False,
)
```

## 12B.2 同じことを `quicklook()` でより短く書く

```python
result = quicklook(
    out_fits,
    mode='moment0',
    vel_range=(0.0, 20.0),
    spectral_unit='km/s',
    target_hpbw_arcsec=450.0,
    cmap='turbo',
    norm_mode='asinh',
    norm_percentile=(1.0, 99.5),
    contours=[
        {
            'levels': 'auto',
            'colors': 'white',
            'linewidths': 0.8,
            'alpha': 0.6,
        }
    ],
    beam='header',
    title='Moment 0 with target HPBW',
    show=False,
)
```

## 12B.3 今後を見据えて `plot_scene()` に寄せる

```python
result = plot_scene(
    base={
        'kind': 'image',
        'source': out_fits,
        'mode': 'moment0',
        'vel_range': (0.0, 20.0),
        'spectral_unit': 'km/s',
        'target_hpbw_arcsec': 450.0,
        'cmap': 'turbo',
        'norm_mode': 'asinh',
        'norm_percentile': (1.0, 99.5),
        'beam': 'header',
    },
    overlays=[
        {
            'kind': 'contour',
            'source': out_fits,
            'mode': 'moment0',
            'vel_range': (0.0, 20.0),
            'spectral_unit': 'km/s',
            'target_hpbw_arcsec': 450.0,
            'levels': 'auto',
            'colors': 'white',
            'linewidths': 0.8,
            'alpha': 0.6,
        }
    ],
    title='Moment 0 with target HPBW',
    show=False,
)
```

---

# 13. 付録: 最小レシピ

### 2D FITS

```python
pl.quicklook('moment0.fits')
```

### 3D cube を速度積分

```python
pl.quicklook('cube.fits', vel_range=(20, 40), title='Integrated intensity')
```

### optical + radio contour

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits'},
    overlays=[
        {'kind': 'contour', 'source': '12co.fits', 'mode': 'moment0', 'vel_range': (20, 40), 'colors': 'cyan'},
        {'kind': 'beam', 'beam': 'auto'},
        {'kind': 'scalebar', 'config': {'length': 30, 'units': 'arcsec'}},
    ],
)
```

### RGB

```python
rgb = pl.make_rgb_map(
    {'source': '12co.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
    {'source': '13co.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
    {'source': 'c18o.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
)
pl.plot_rgb(rgb, title='12CO/13CO/C18O RGB')
```
