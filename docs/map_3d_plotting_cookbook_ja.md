# `sd_radio_spectral_fits.map_3d.plotting` Cookbook（日本語・preliminary）

> **重要な注記**  
> 本 cookbook の内容は **preliminary** な `plotting.py` 実装を前提としています。将来、API や推奨書き方が大きく変わる可能性があります。  
> 特に、`vmin/vmax` の名称、WCS 不一致時の扱い、extension 名、mask の意味づけは変更候補です。

---

## 1. まず最初に: 典型的な考え方

本モジュールの使い方は、まず次の 3 種に分けて考えると分かりやすいです。

1. **既にある 2D map を読む**  
   例: `MOMENT0`, `RMS`, `HIT`

2. **3D cube から新しく 2D map を作る**  
   例: `moment0`, `channel_sum`, `channel_mean`

3. **provisional / final moment を作る**  
   例: `LINEFREE`, `LINECAND3D`, `BASESUP3D`, `MASK3D` を使う

---

## 2. 最小例: 3D cube から moment0 を作って描く

```python
from sd_radio_spectral_fits.map_3d.plotting import make_2d_map, plot_map

map2d = make_2d_map(
    "ps_12co_cube.fits",
    ext=0,
    mode="moment0",
    vel_range=(0.0, 20.0),
    spectral_unit="km/s",
    target_hpbw_arcsec=400.0,
    orig_hpbw_arcsec=350.0,
)

result = plot_map(
    map2d,
    cmap="turbo",
    norm_mode="asinh",
    norm_percentile=(1.0, 99.5),
    stretch_a=0.10,
    title="12CO moment0",
    show=False,
)
```

### ポイント

- `ext=0` は 3D cube HDU の例です
- `vel_range` は 3D を 2D に積分するときに有効です
- `target_hpbw_arcsec` を使うなら、元ビームが不明な場合は `orig_hpbw_arcsec` を与えるのが安全です

---

## 3. 既存 `MOMENT0` extension をそのまま描く

```python
map2d = make_2d_map(
    "ps_12co_cube_masked.fits",
    ext="MOMENT0",
    mode="map",
    target_hpbw_arcsec=400.0,
    orig_hpbw_arcsec=350.0,
)

result = plot_map(
    map2d,
    cmap="turbo",
    norm_mode="asinh",
    norm_percentile=(1.0, 99.5),
    stretch_a=0.10,
    title="Existing MOMENT0 extension",
    show=False,
)
```

### ポイント

- `ext="MOMENT0"` はすでに 2D です
- この場合 `vel_range` は使われません
- もし `vel_range` を付けても、「3D ではないので無視される」というのが自然な解釈です

---

## 4. `RMS` や `HIT` をそのまま描く

### 4.1 `RMS`

```python
rms_map = make_2d_map("ps_12co_cube_masked.fits", ext="RMS", mode="map")

plot_map(
    rms_map,
    cmap="magma",
    norm_mode="log",
    norm_percentile=(1.0, 99.5),
    title="RMS map",
    show=False,
)
```

### 4.2 `HIT`

```python
hit_map = make_2d_map("ps_12co_cube_masked.fits", ext="HIT", mode="map")

plot_map(
    hit_map,
    cmap="viridis",
    norm_mode="linear",
    title="Hit count map",
    show=False,
)
```

### コメント

- `RMS` は `linear` でも `log` でもよいですが、ダイナミックレンジが大きいときは `log` が見やすいです
- `HIT` は通常 `linear` が分かりやすいです

---

## 5. 2D extension と 3D 再計算の違いを並べて確認する

```python
m_existing = make_2d_map("cube_masked.fits", ext="MOMENT0", mode="map")

m_recalc = make_2d_map(
    "cube_masked.fits",
    ext=0,
    mode="moment0",
    vel_range=(0, 20),
    spectral_unit="km/s",
)

plot_map(m_existing, title="Existing MOMENT0", show=False)
plot_map(m_recalc, title="Recalculated moment0 (0-20 km/s)", show=False)
```

### 何を見るか

- 既存 `MOMENT0` がどの速度範囲・mask 定義で作られたか
- 再計算版と一致するか
- 既存 `MOMENT0` が smoothing 済みかどうか

---

## 6. 単一チャネルを描く

```python
ch10 = make_2d_map(
    "ps_12co_cube.fits",
    ext=0,
    mode="channel_slice",
    chan_range=(10, 10),
)

plot_map(
    ch10,
    cmap="turbo",
    norm_mode="linear",
    title="Channel 10",
    show=False,
)
```

### コメント

実装によっては `mode="channel"` や `mode="slice"` も使える設計になっていることがあります。現在の preliminary 実装では、どの別名が許されるかは将来整理される可能性があります。

---

## 7. チャネル和とチャネル平均

### 7.1 channel sum

```python
csum = make_2d_map(
    "ps_12co_cube.fits",
    ext=0,
    mode="channel_sum",
    chan_range=(30, 50),
)

plot_map(csum, cmap="turbo", norm_mode="asinh", title="Channel sum 30-50", show=False)
```

### 7.2 channel mean

```python
cmean = make_2d_map(
    "ps_12co_cube.fits",
    ext=0,
    mode="channel_mean",
    chan_range=(30, 50),
)

plot_map(cmean, cmap="turbo", norm_mode="asinh", title="Channel mean 30-50", show=False)
```

---

## 8. provisional moment を作る

### 8.1 自動探索 `prefer="auto"`

```python
pmap = make_provisional_moment(
    "ps_12co_cube_masked.fits",
    ext=0,
    prefer="auto",
    vel_range=(0, 20),
    spectral_unit="km/s",
    nan_fill=True,
)

plot_map(
    pmap,
    cmap="turbo",
    norm_mode="asinh",
    norm_percentile=(1, 99.5),
    title="Provisional moment (auto)",
    show=False,
)
```

### `auto` の意味をもう一度

- `LINECAND3D` があればそれを使う
- 無ければ `BASESUP3D` の補集合を使う
- それも無ければ `LINEFREE` の補集合を使う

つまり、「signal 候補を直接示す情報があるならそれを優先し、無ければ baseline 側の情報から signal 側を組み立てる」という意味です。

### 8.2 `LINEFREE` を明示優先する

```python
pmap_linefree = make_provisional_moment(
    "ps_12co_cube_masked.fits",
    ext=0,
    prefer="linefree",
)
```

### 8.3 `BASESUP3D` を明示優先する

```python
pmap_basesup = make_provisional_moment(
    "ps_12co_cube_masked.fits",
    ext=0,
    prefer="basesup3d",
)
```

### 8.4 `LINECAND3D` を明示優先する

```python
pmap_linecand = make_provisional_moment(
    "ps_12co_cube_masked.fits",
    ext=0,
    prefer="linecand3d",
)
```

---

## 9. final moment を作る

```python
fmap = make_final_moment(
    "ps_12co_cube_masked.fits",
    ext=0,
    final_mask_ext="MASK3D",
    vel_range=(0, 20),
    spectral_unit="km/s",
    nan_fill=True,
)

plot_map(
    fmap,
    cmap="turbo",
    norm_mode="asinh",
    norm_percentile=(1, 99.5),
    title="Final moment (MASK3D)",
    show=False,
)
```

---

## 10. provisional と final を並べる

```python
pmap = make_provisional_moment("cube_masked.fits", ext=0, prefer="auto")
fmap = make_final_moment("cube_masked.fits", ext=0)

plot_map(pmap, title="Provisional moment", show=False)
plot_map(fmap, title="Final moment", show=False)
```

### 見るべき点

- provisional の方が広めに出るか
- final の方が保守的になるか
- baseline 起源の偽構造が provisional に残るか

---

## 11. 0 埋めと NaN 埋めを比べる

### 11.1 NaN 埋め

```python
pmap_nan = make_provisional_moment(
    "cube_masked.fits",
    ext=0,
    nan_fill=True,
    zero_fill=False,
)
```

### 11.2 0 埋め

```python
pmap_zero = make_provisional_moment(
    "cube_masked.fits",
    ext=0,
    nan_fill=False,
    zero_fill=True,
)
```

### コメント

- `NaN` 埋めは「ここは無効領域」と分かりやすいです
- 0 埋めは積分値や平均値の解釈に影響することがあります
- 科学的解釈が絡む場合、どちらで作ったかを必ず記録してください

---

## 12. smoothing の使い方

### 12.1 追加 kernel を直接指定する

```python
m = make_2d_map(
    "cube.fits",
    ext=0,
    mode="moment0",
    vel_range=(0, 20),
    smooth_fwhm_arcsec=120.0,
)
```

### 12.2 最終 HPBW を指定する

```python
m = make_2d_map(
    "cube.fits",
    ext=0,
    mode="moment0",
    vel_range=(0, 20),
    target_hpbw_arcsec=400.0,
    orig_hpbw_arcsec=350.0,
)
```

### 12.3 header に beam が無い場合

```python
m = make_2d_map(
    "cube.fits",
    ext=0,
    mode="moment0",
    vel_range=(0, 20),
    target_hpbw_arcsec=400.0,
    orig_hpbw_arcsec=350.0,
)
```

### コメント

header に `BMAJ/BMIN` が無いと、`target_hpbw_arcsec` だけでは厳密な追加 kernel が決まりません。`orig_hpbw_arcsec` を明示してください。

---

## 13. beam を図中に描く

### 13.1 header / smoothing 情報から自動表示

```python
plot_map(
    map2d,
    beam="header",
    title="Beam from header or smooth_info",
    show=False,
)
```

### 13.2 beam を手動指定

```python
plot_map(
    map2d,
    beam={
        "major": 400.0,
        "minor": 350.0,
        "angle": 30.0,
        "units": "arcsec",
        "x": 40,
        "y": 40,
        "edgecolor": "white",
        "linewidth": 1.2,
    },
    show=False,
)
```

---

## 14. normalize の典型例

### 14.1 moment0 で asinh

```python
plot_map(
    map2d,
    norm_mode="asinh",
    norm_percentile=(1.0, 99.5),
    stretch_a=0.10,
    show=False,
)
```

### 14.2 RMS で log

```python
plot_map(
    rms_map,
    norm_mode="log",
    norm_percentile=(1.0, 99.5),
    show=False,
)
```

### 14.3 手動 `vmin/vmax`

```python
plot_map(
    map2d,
    norm_mode="linear",
    vmin=0.0,
    vmax=200.0,
    show=False,
)
```

### 注意

現行実装の `vmin/vmax` は **速度ではなく color scale 範囲**です。

---

## 15. 同じ map に contour を重ねる

```python
plot_map(
    map2d,
    cmap="turbo",
    norm_mode="asinh",
    norm_percentile=(1.0, 99.5),
    contours=[
        {
            "levels": "auto",
            "colors": "white",
            "linewidths": 0.8,
            "alpha": 0.6,
        }
    ],
    title="Moment0 with contours",
    show=False,
)
```

---

## 16. 別データセットを contour にする

```python
plot_map(
    map2d,
    cmap="turbo",
    norm_mode="asinh",
    norm_percentile=(1.0, 99.5),
    contours=[
        {
            "source": "other_cube.fits",
            "ext": 0,
            "mode": "moment0",
            "vel_range": (0, 20),
            "spectral_unit": "km/s",
            "levels": [20, 40, 80, 160],
            "colors": "cyan",
            "linewidths": 0.8,
            "alpha": 0.7,
        }
    ],
    show=False,
)
```

---

## 17. 複数 contour を重ねる

```python
plot_map(
    map2d,
    cmap="gray",
    norm_mode="linear",
    contours=[
        {
            "levels": [20, 40, 80],
            "colors": "white",
            "linewidths": 0.8,
            "label": "self",
        },
        {
            "source": "13co_cube.fits",
            "ext": 0,
            "mode": "moment0",
            "vel_range": (0, 20),
            "spectral_unit": "km/s",
            "levels": [10, 20, 40],
            "colors": "cyan",
            "linewidths": 0.8,
            "label": "13CO",
        },
        {
            "source": "c18o_cube.fits",
            "ext": 0,
            "mode": "moment0",
            "vel_range": (0, 20),
            "spectral_unit": "km/s",
            "levels": [5, 10, 20],
            "colors": "magenta",
            "linewidths": 0.8,
            "label": "C18O",
        },
    ],
    show=False,
)
```

---

## 18. `header, data` から描く

```python
from astropy.io import fits

hdu = fits.open("map2d.fits")["MOMENT0"]
header = hdu.header
data = hdu.data

plot_map(
    data=data,
    header=header,
    cmap="turbo",
    norm_mode="asinh",
    norm_percentile=(1, 99.5),
    title="From (header, data)",
    show=False,
)
```

### 使いどころ

- Python 内で中間 map を作った
- しかし FITS に保存する前に描きたい
- 既存 pipeline の戻り値をそのまま可視化したい

---

## 19. Python で作った 2D 配列に contour を重ねる

```python
plot_map(
    data=my_data,
    header=my_header,
    cmap="viridis",
    contours=[
        {
            "data": my_contour_data,
            "header": my_contour_header,
            "levels": [1, 2, 4, 8],
            "colors": "white",
        }
    ],
    show=False,
)
```

---

## 20. 銀河座標の map を描く

```python
plot_map(
    "galactic_map.fits",
    ext=0,
    cmap="turbo",
    norm_mode="asinh",
    title="Galactic coordinates",
    show=False,
)
```

### コメント

WCS が `GLON/GLAT` なら、既定の軸ラベルは銀河座標系になります。

---

## 21. 赤道座標の map を描く

```python
plot_map(
    "equatorial_map.fits",
    ext=0,
    cmap="turbo",
    norm_mode="asinh",
    title="Equatorial coordinates",
    show=False,
)
```

---

## 22. RGB 合成の最小例

```python
from sd_radio_spectral_fits.map_3d.plotting import make_rgb_map, plot_rgb

rgb = make_rgb_map(
    red=make_2d_map("r_map.fits", ext=0, mode="map"),
    green=make_2d_map("g_map.fits", ext=0, mode="map"),
    blue=make_2d_map("b_map.fits", ext=0, mode="map"),
)

plot_rgb(rgb, title="RGB composite", show=False)
```

---

## 23. 3D cube から RGB 合成

```python
rgb = make_rgb_map(
    red={
        "source": "cube.fits",
        "ext": 0,
        "mode": "moment0",
        "vel_range": (-10, 0),
        "spectral_unit": "km/s",
    },
    green={
        "source": "cube.fits",
        "ext": 0,
        "mode": "moment0",
        "vel_range": (0, 10),
        "spectral_unit": "km/s",
    },
    blue={
        "source": "cube.fits",
        "ext": 0,
        "mode": "moment0",
        "vel_range": (10, 20),
        "spectral_unit": "km/s",
    },
    red_norm={"norm_mode": "asinh", "percentile": (1, 99.5)},
    green_norm={"norm_mode": "asinh", "percentile": (1, 99.5)},
    blue_norm={"norm_mode": "asinh", "percentile": (1, 99.5)},
)

plot_rgb(rgb, title="Velocity-split RGB", show=False)
```

---

## 24. extension ごとに RGB を作る

```python
rgb = make_rgb_map(
    red={"source": "cube_masked.fits", "ext": "MOMENT0", "mode": "map"},
    green={"source": "cube_masked.fits", "ext": "RMS", "mode": "map"},
    blue={"source": "cube_masked.fits", "ext": "HIT", "mode": "map"},
    red_norm={"norm_mode": "asinh", "percentile": (1, 99.5)},
    green_norm={"norm_mode": "log", "percentile": (1, 99.5)},
    blue_norm={"norm_mode": "linear", "percentile": (1, 99.5)},
)

plot_rgb(rgb, title="RGB from extensions", show=False)
```

### 注意

科学的な意味が直感的ではないので、RGB に何を割り当てたかを必ず図注に書くことを推奨します。

---

## 25. `Map2D` を受け渡して後から描く

```python
m = make_2d_map("cube.fits", ext=0, mode="moment0", vel_range=(0, 20))

# いったん別セル・別関数へ渡す
result = plot_map(m, cmap="turbo", norm_mode="asinh", show=False)
```

### 利点

- 2D 化と描画を分離できる
- 同じ `Map2D` を別設定で何度も描ける
- contour や RGB の元として再利用できる

---

## 26. 同じ `Map2D` を別 normalize で見比べる

```python
m = make_2d_map("cube.fits", ext=0, mode="moment0", vel_range=(0, 20))

plot_map(m, norm_mode="linear", title="linear", show=False)
plot_map(m, norm_mode="sqrt", title="sqrt", show=False)
plot_map(m, norm_mode="log", title="log", show=False)
plot_map(m, norm_mode="asinh", title="asinh", show=False)
```

### 何を見るか

- 強ピークと弱い拡散成分のバランス
- ノイズ床の見え方
- contour と pseudo color の整合

---

## 27. `norm` を外で作って渡す

```python
from sd_radio_spectral_fits.map_3d.plotting import build_normalize

norm = build_normalize(
    map2d.data,
    mode="asinh",
    percentile=(1, 99.5),
    stretch_a=0.08,
)

plot_map(
    map2d,
    cmap="turbo",
    norm=norm,
    title="External normalize",
    show=False,
)
```

---

## 28. 既存 WCSAxes に描く

```python
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=map2d.wcs)

plot_map(map2d, ax=ax, show=False)
```

### 使いどころ

- 複数 panel の図を作る
- subplot を組む
- 既存 Figure に追加する

---

## 29. Multi-panel を手で組む

```python
import matplotlib.pyplot as plt

m0 = make_2d_map("cube_masked.fits", ext="MOMENT0", mode="map")
rms = make_2d_map("cube_masked.fits", ext="RMS", mode="map")
hit = make_2d_map("cube_masked.fits", ext="HIT", mode="map")

fig = plt.figure(figsize=(15, 4))
ax1 = fig.add_subplot(131, projection=m0.wcs)
ax2 = fig.add_subplot(132, projection=rms.wcs)
ax3 = fig.add_subplot(133, projection=hit.wcs)

plot_map(m0, ax=ax1, cmap="turbo", norm_mode="asinh", show=False, title="MOMENT0")
plot_map(rms, ax=ax2, cmap="magma", norm_mode="log", show=False, title="RMS")
plot_map(hit, ax=ax3, cmap="viridis", norm_mode="linear", show=False, title="HIT")
```

---

## 30. CLI を使う

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

## 31. よくある落とし穴

### 31.1 `ext="MOMENT0"` に `vel_range` を付けた

これは通常、意味がありません。`MOMENT0` はすでに 2D だからです。

### 31.2 `vmin/vmax` を速度と勘違いした

違います。表示スケール範囲です。

### 31.3 `target_hpbw_arcsec` だけ指定して header に beam が無い

`orig_hpbw_arcsec` を明示するのが安全です。

### 31.4 RGB で shape が違う

現行実装では基本的に許されません。

### 31.5 RGB で WCS が少し違う

warning は出ますが、reproject はしません。科学的比較には注意が必要です。

### 31.6 provisional と final の意味を混同した

解析上かなり重要な違いです。図タイトルや保存ファイル名でも区別することを推奨します。

---

## 32. 実務的な命名の勧め

図や保存ファイル名では、少なくとも次を明記するのが安全です。

- `moment0` か `provisional` か `final` か
- 速度範囲またはチャネル範囲
- smoothing 条件
- normalize 条件
- 使用した extension / mask 名

例:

- `co10_moment0_v0to20_hpbw400_asinh.png`
- `co10_provisional_auto_linecand_or_basesup.png`
- `co10_final_mask3d_v0to20.png`

---

## 33. 最後に

この cookbook は、**map_3d 出力を現場で素早く確認する**ことを第一目的にしています。provisional/final moment や beam の扱いなど、解析解釈に関わる重要点を含むため、必ず「何を入力し、どの mask / extension を使い、どの速度範囲で 2D 化したか」を明示して使ってください。
