# plotting_updated_v5p5 超詳細 cookbook（日本語・拡充版）

> 対象コード: **`plotting_updated_v5p5.py`**  
> 本書は、初心者が最初の 1 枚を出すところから、研究用途での複雑な overlay・catalog・RGB・reprojection・CLI 利用まで、できるだけ多くの実例を集めた cookbook です。  
> 仕様の厳密な定義・全パラメータ一覧は `plotting_updated_v5p5_full_manual_ja_current_merged.md` を参照してください。本書は **使い方中心** です。  
> ただし、現行 v5p5 で重要になった **`quicklook(contours=...)`**、**contour level の `auto / fraction / manual / rms` 整理**、**CLI contour の新しい挙動**、**`plot_map()` は互換 API として残し、新規コードでは `quicklook()` / `plot_scene()` を優先する方針** は、本書でも実例を増やして明示しています。

---

# 0. この cookbook の考え方

このモジュールには大きく分けて 4 つの入口があります。

- **`quicklook()`**  
  とにかく 1 枚すぐ見たい時の入口。**現行版では最初に考えるべき入口**です。
- **`plot_map()`**  
  2D map を 1 枚描く低〜中水準 API。**既存コードとの互換性維持のために重要**です。
- **`plot_rgb()`**  
  RGB 合成を描く API。これも互換 API として重要です。
- **`plot_scene()`**  
  `base + overlays` で複数レイヤーを統一的に扱う高水準 API。**今後の新規コードの中心**です。

また、2D map を事前に作りたい時には次を使います。

- **`make_2d_map()`**
- **`make_provisional_moment()`**
- **`make_final_moment()`**

現時点でのおすすめの考え方は次のとおりです。

- **まず 1 枚出したい** → `quicklook()`
- **既存の `make_2d_map() + plot_map()` 系コードを保ちたい** → `plot_map()`
- **光学画像 + 電波 contour、catalog、annotation、将来の拡張まで見据える** → `plot_scene()`

`plot_map()` は今でも有用ですが、設計上は **過去との互換性を維持するために残している API** という位置づけです。したがって、今後新しくコードを書くなら、まず **`quicklook()` と `plot_scene()` を優先して考える**のが推奨です。

本書では、まず簡単なものから始め、後半に高度な例を増やしていきます。その際、旧来の `plot_map()` の例も残しますが、必要に応じて **「同じことを現行推奨 API で書くならどうするか」** も併記します。

---
# 1. まず最初に import

```python
from sd_radio_spectral_fits.map_3d.plotting as pl
```

Jupyter で使うなら、最初に以下もよく使います。

```python
%matplotlib inline
# あるいは
# %matplotlib widget
```

---

# 2. 最短で 1 枚出す

## 2.1 2D FITS をそのまま表示

```python
pl.quicklook('moment0.fits')
```

何が起こるか:

- 入力が 2D ならそのまま表示
- `asinh` 正規化
- percentile clip は既定 `(1.0, 99.5)`
- 座標グリッド ON
- readout ON

---

## 2.2 タイトル付きで表示

```python
pl.quicklook('moment0.fits', title='12CO moment0')
```

---

## 2.3 画面に出さず PDF 保存

```python
pl.quicklook('moment0.fits', save='moment0.pdf', show=False)
```

---

## 2.4 PNG 保存 + 透明背景

```python
pl.quicklook(
    'moment0.fits',
    save='moment0.png',
    save_transparent=True,
    show=False,
)
```

---

## 2.5 beam を自動表示

```python
pl.quicklook('moment0.fits', beam='auto')
```

beam はヘッダから拾える場合に描画されます。

---

## 2.6 すぐ説明を出す

```python
pl.quicklook('moment0.fits', describe=True)
```

これは描画と同時に、表示した map の簡単な説明も標準出力に出します。

---

## 2.7 Python API の簡単な help を出す

```python
pl.quicklook('moment0.fits', help=True)
```

または

```python
print(pl.plotting_help_text())
```

---

## 2.8 `quicklook()` で contour を重ねる（現行版の重要変更）

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(0, 20),
    contours=[
        {
            'levels': 'auto',
            'colors': 'white',
            'linewidths': 0.8,
            'alpha': 0.7,
        }
    ],
)
```

`quicklook()` は現行版で **`contours=` に正式対応**しています。
base 側と同じ source, mode, vel_range, spectral_unit, smoothing 設定を contour 側に継承するので、最短の重ね描きにはこれが便利です。

---

## 2.9 `quicklook()` で contour の source だけ変える

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(0, 20),
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (0, 20),
            'colors': 'cyan',
            'levels': 'auto',
        }
    ],
)
```

12CO を base に、13CO を contour にしたい時の最短例です。

---

## 2.10 `quicklook()` で contour の積分範囲だけ変える

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(0, 20),
    contours=[
        {
            'mode': 'moment0',
            'vel_range': (5, 15),
            'colors': 'yellow',
            'levels': 'auto',
        }
    ],
)
```

contour 側で `source` を省略すると base と同じ source を使います。

---

## 2.11 `quicklook()` で beam と contour を同時に使う

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(0, 20),
    beam='header',
    contours=[{'levels': 'auto', 'colors': 'white'}],
    title='quicklook with contour and beam',
)
```

---

## 2.12 `quicklook()` を使うべきか、`plot_map()` を使うべきか

- **最初の 1 枚**、あるいは **base + contour 程度まで**なら `quicklook()` が最も簡単です。
- **既存コード資産をそのまま活かしたい**なら `plot_map()` が自然です。
- **複数 contour、marker、catalog、annotation、image overlay まで育てたい**なら `plot_scene()` を最初から使う方が後で楽です。

---

# 3. 入力形式のバリエーション

このモジュールは、ファイルだけでなく Python 上のデータも扱えます。

## 3.1 path 文字列

```python
pl.quicklook('cube.fits')
```

---

## 3.2 `pathlib.Path`

```python
from pathlib import Path
pl.quicklook(Path('cube.fits'))
```

---

## 3.3 `(header, data)`

```python
pl.quicklook((header, data))
```

---

## 3.4 `plot_map(data=..., header=...)`

```python
pl.plot_map(data=data2d, header=header2d, title='in-memory 2D map')
```

---

## 3.5 HDU を使う

```python
from astropy.io import fits

hdu = fits.open('cube.fits')['MOMENT0']
pl.quicklook(hdu)
```

---

## 3.6 `Map2D` を使う

```python
m = pl.make_2d_map('cube.fits', mode='moment0', vel_range=(20, 40))
pl.plot_map(m)
```

---

## 3.7 `RGBMap` を使う

```python
rgb = pl.make_rgb_map('r.fits', 'g.fits', 'b.fits')
pl.plot_rgb(rgb)
```

---

## 3.8 `OTFBundle` を使う

```python
bundle = pl.read_otf_bundle('otf_bundle.fits')
pl.quicklook(bundle, ext='MOMENT0')
```

ここでのポイント:

- Python 上の `OTFBundle` を **そのまま** source として渡せます
- 2D product をそのまま見るときは `ext='MOMENT0'` のように指定します
- `ext` は `MOMENT0`, `RMS`, `HIT`, `MOSAIC_WEIGHT` など bundle 内の image ext 名を使えます

---

## 3.9 baseline viewer bundle を使う

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
viewer = pl.make_baseline_viewer_bundle(bundle, mode='signal', fill_value=0.0)
pl.quicklook(viewer, mode='moment0', vel_range=(20, 40))
```

重要:

- viewer bundle は特別な型ではなく、**通常の `OTFBundle`** です
- ただし `viewer.data` は 3D cube なので、2D 表示したいときは `mode='moment0'` などで 2D 化してから描きます
- `plot_map(viewer)` のように 3D bundle をそのまま 2D 扱いするのではありません

---

## 3.10 `HDUList` を使う

```python
from astropy.io import fits

hdul = fits.open('products.fits')
pl.quicklook(hdul, ext='MOMENT0')
```

`HDU` だけでなく `HDUList` もそのまま渡せます。

---

## 3.11 OTF bundle FITS を path のまま使う

```python
pl.quicklook('otf_bundle.fits', ext='MOMENT0')
```

ディスク上の OTF bundle FITS を使うだけなら、必ずしも `read_otf_bundle()` で Python object に直さなくても構いません。

- path / `HDUList` / `OTFBundle` のどれを使うかは、手元の処理フローに合わせて選べます
- ただし **viewer bundle を一度 Python 上で作ってから描きたい**ときは、`OTFBundle` object を直接渡せる方が自然です

---

# 4. 3D cube から 2D を作る基本例

## 4.1 全 channel の moment0

```python
pl.quicklook('cube.fits', mode='moment0')
```

---

## 4.2 速度範囲を指定した moment0

```python
pl.quicklook('cube.fits', mode='moment0', vel_range=(20, 40))
```

重要:

- `vel_range` は **積分する速度範囲**
- 表示色の範囲ではありません

---

## 4.3 順序が逆でもよい

```python
pl.quicklook('cube.fits', mode='moment0', vel_range=(40, 20))
```

内部で並べ替えられます。

---

## 4.4 channel 範囲を使う

```python
pl.quicklook('cube.fits', mode='moment0', chan_range=(100, 180))
```

---

## 4.5 channel mean

```python
pl.quicklook('cube.fits', mode='channel_mean', vel_range=(20, 40))
```

---

## 4.6 channel sum

```python
pl.quicklook('cube.fits', mode='channel_sum', vel_range=(20, 40))
```

---

## 4.7 1 channel だけ見る

```python
pl.quicklook('cube.fits', mode='channel_slice', chan_range=(120, 120))
```

---

## 4.8 1 velocity 面だけ見る

```python
pl.quicklook('cube.fits', mode='velocity_slice', vel_range=(25.0, 25.0))
```

---

## 4.9 幅を持たせた velocity_slice

```python
pl.quicklook('cube.fits', mode='velocity_slice', vel_range=(24.5, 25.5))
```

内部的には中心付近の 1 面が選ばれます。

---

## 4.10 `make_2d_map()` を使ってから描く

```python
m = pl.make_2d_map(
    'cube.fits',
    mode='moment0',
    vel_range=(20, 40),
)
print(pl.describe_map2d(m))
pl.plot_map(m)
```

---

# 5. precomputed 2D products を使う

## 5.1 `MOMENT0` extension をそのまま描く

```python
pl.quicklook('products.fits', ext='MOMENT0')
```

この場合は 2D passthrough です。

---

## 5.2 `LINEFREE` など別 extension を確認

```python
pl.quicklook('products.fits', ext='MOMENT1')
```

2D であればそのまま描けます。

---

## 5.3 `identity` を明示する

```python
pl.quicklook('products.fits', ext='MOMENT0', mode='identity')
```

---

## 5.4 OTF bundle の `MOMENT0` extension をそのまま描く

```python
bundle = pl.read_otf_bundle('otf_bundle.fits')
pl.quicklook(bundle, ext='MOMENT0')
```

これは bundle の 2D extension を **2D passthrough** でそのまま表示する例です。

---

## 5.5 OTF bundle の `RMS` や `HIT` を確認する

```python
bundle = pl.read_otf_bundle('otf_bundle.fits')
pl.quicklook(bundle, ext='RMS', cmap='magma')
pl.quicklook(bundle, ext='HIT', cmap='viridis')
```

観測の被覆やノイズの確認に便利です。

---

## 5.6 `LINEFREE_USED` や `SIGNAL_MASK_USED` は 2D map ではない

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
print(bundle.image_ext['LINEFREE_USED'].shape)
print(bundle.image_ext['SIGNAL_MASK_USED'].shape)
```

ここで重要なのは、現行仕様では

- `LINEFREE`, `LINEFREE_USED`
- `SIGNAL_MASK_USED`

は **global な 1D channel mask** であり、2D map でも 3D voxel mask でもないという点です。

したがって、次のような使い方はしません。

```python
# これは意図に合わない
# pl.quicklook(bundle, ext='LINEFREE_USED')
```

代わりに、viewer bundle を作ってから 3D cube を 2D 化して確認します。

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
viewer = pl.make_baseline_viewer_bundle(bundle, mode='linefree', fill_value=0.0)
pl.quicklook(viewer, mode='moment0', vel_range=(20, 40))
```

---

## 5.7 `identity` を OTF bundle の 2D ext に対して明示する

```python
bundle = pl.read_otf_bundle('otf_bundle.fits')
pl.quicklook(bundle, ext='MOMENT0', mode='identity')
```

2D ext をそのまま描いていることを、呼び出し側で明示したいときの書き方です。

---

## 5.8 `HDUList` + `ext` で OTF bundle 風 products を扱う

```python
from astropy.io import fits

hdul = fits.open('otf_bundle.fits')
pl.quicklook(hdul, ext='MOMENT0')
```

FITS 化された products を読むだけなら、この書き方でも十分です。

---

# 6. provisional / final moment の例

## 6.1 provisional moment を自動選択で作る

```python
m = pl.make_provisional_moment('cube_with_masks.fits')
pl.plot_map(m, title='provisional moment')
```

`prefer='auto'` では、概ね次の順で使われます。

1. `LINECAND3D`
2. `BASESUP3D` の補集合
3. `LINEFREE` の補集合

---

## 6.2 provisional moment + 速度範囲制限

```python
m = pl.make_provisional_moment(
    'cube_with_masks.fits',
    vel_range=(20, 40),
)
pl.plot_map(m)
```

---

## 6.3 final moment を `MASK3D` から作る

```python
m = pl.make_final_moment('cube_with_masks.fits')
pl.plot_map(m, title='final moment from MASK3D')
```

---

## 6.4 final moment + smoothing

```python
m = pl.make_final_moment(
    'cube_with_masks.fits',
    target_hpbw_arcsec=30.0,
)
pl.plot_map(m)
```

---

## 6.5 provisional moment を OTF bundle から作る

```python
bundle = pl.read_otf_bundle('cube_with_masks_bundle.fits')
m = pl.make_provisional_moment(bundle)
pl.plot_map(m, title='provisional moment from bundle')
```

`prefer='auto'` の探索順は、FITS path のときと同じです。

---

## 6.6 final moment を OTF bundle から作る

```python
bundle = pl.read_otf_bundle('cube_with_masks_bundle.fits')
m = pl.make_final_moment(bundle)
pl.plot_map(m, title='final moment from bundle')
```

`MASK3D` が bundle 内にあれば、そのまま使えます。

---

## 6.7 baseline viewer bundle を provisional / final moment と混同しない

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
viewer = pl.make_baseline_viewer_bundle(bundle, mode='signal')

# viewer は 3D cube source として使う
pl.quicklook(viewer, mode='moment0', vel_range=(20, 40))
```

viewer bundle は

- `SIGNAL_MASK_USED` や `LINEFREE_USED` を使って `data` を選別した 3D bundle
- `MASK3D` を直接持つ final-moment source そのものではない

という点を区別してください。

---

# 7. `make_2d_map()` の mask の使い方

## 7.1 外部 mask を使う

```python
m = pl.make_2d_map(
    'cube.fits',
    mode='moment0',
    mask=my_mask_3d,
    mask_mode='external',
)
pl.plot_map(m)
```

`my_mask_3d` は cube と整合する shape が必要です。

---

## 7.2 `LINEFREE` の補集合を使う

```python
m = pl.make_2d_map(
    'cube_with_helper.fits',
    mode='moment0',
    mask_mode='linefree_complement',
)
pl.plot_map(m)
```

---

## 7.3 `LINECAND3D` を直接使う

```python
m = pl.make_2d_map(
    'cube_with_helper.fits',
    mode='moment0',
    mask_mode='linecand3d',
)
pl.plot_map(m)
```

---

## 7.4 `BASESUP3D` の補集合を使う

```python
m = pl.make_2d_map(
    'cube_with_helper.fits',
    mode='moment0',
    mask_mode='basesup3d',
)
pl.plot_map(m)
```

---

## 7.5 `MASK3D` を使う

```python
m = pl.make_2d_map(
    'cube_with_helper.fits',
    mode='moment0',
    mask_mode='mask3d',
)
pl.plot_map(m)
```

---

## 7.6 マスク外を 0 埋め

```python
m = pl.make_2d_map(
    'cube_with_helper.fits',
    mode='moment0',
    mask_mode='mask3d',
    zero_fill=True,
    nan_fill=False,
)
pl.plot_map(m)
```

---

## 7.7 マスク外を NaN 埋め（既定）

```python
m = pl.make_2d_map(
    'cube_with_helper.fits',
    mode='moment0',
    mask_mode='mask3d',
    nan_fill=True,
)
pl.plot_map(m)
```

---

## 7.8 独自 fill 値を使う

```python
m = pl.make_2d_map(
    'cube_with_helper.fits',
    mode='moment0',
    mask_mode='mask3d',
    zero_fill=False,
    nan_fill=False,
    fill_value=-999.0,
)
pl.plot_map(m)
```

---

## 7.9 OTF bundle の `LINEFREE` 補集合を使う

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
m = pl.make_2d_map(
    bundle,
    mode='moment0',
    mask_mode='linefree_complement',
)
pl.plot_map(m)
```

---

## 7.10 OTF bundle の `LINECAND3D` を直接使う

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
m = pl.make_2d_map(
    bundle,
    mode='moment0',
    mask_mode='linecand3d',
)
pl.plot_map(m)
```

---

## 7.11 OTF bundle の `BASESUP3D` 補集合を使う

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
m = pl.make_2d_map(
    bundle,
    mode='moment0',
    mask_mode='basesup3d',
)
pl.plot_map(m)
```

---

## 7.12 `LINEFREE*` は現行仕様では global 1D mask

ここは特に重要です。

現行仕様では、`LINEFREE`, `LINEFREE_USED`, `SIGNAL_MASK_USED` は

- shape = `(nchan,)`
- cube 全体に共通な channel mask

です。つまり、pixel ごとに異なる `L(c, y, x)` を持つわけではありません。

したがって、

- `mask_mode='linefree_complement'` は **各 pixel で別々の 3D mask** を使う処理ではない
- channel 軸だけの 1D mask を cube 全体へ適用する処理

と理解してください。

---

## 7.13 pixel ごとに違う mask を使いたいなら外部 3D mask を渡す

```python
mask3d = my_per_voxel_mask.astype(bool)

m = pl.make_2d_map(
    bundle,
    mode='moment0',
    mask=mask3d,
    mask_mode='external',
)
pl.plot_map(m)
```

将来的な per-spectrum / per-voxel mask を今すぐ試したいなら、この形が最も明示的です。

---

# 8. smoothing / beam matching

## 8.1 追加 smoothing

```python
pl.quicklook('cube.fits', mode='moment0', smooth_fwhm_arcsec=10.0)
```

---

## 8.2 target HPBW を指定

```python
pl.quicklook('cube.fits', mode='moment0', target_hpbw_arcsec=30.0)
```

---

## 8.3 元の HPBW を自分で与える

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    orig_hpbw_arcsec=17.0,
    target_hpbw_arcsec=30.0,
)
```

---

## 8.4 smoothing した map を別変数に保存

```python
m_smooth = pl.make_2d_map(
    'cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    target_hpbw_arcsec=45.0,
)
pl.plot_map(m_smooth, title='45 arcsec matched')
```

---

# 9. 表示レンジと normalization

## 9.1 percentile clip を狭める

```python
pl.quicklook('moment0.fits', norm_percentile=(5, 99.5))
```

---

## 9.2 percentile clip を広げる

```python
pl.quicklook('moment0.fits', norm_percentile=(0.1, 99.9))
```

---

## 9.3 `linear`

```python
pl.quicklook('moment0.fits', norm_mode='linear')
```

---

## 9.4 `sqrt`

```python
pl.quicklook('moment0.fits', norm_mode='sqrt')
```

---

## 9.5 `log`

```python
pl.quicklook('moment0.fits', norm_mode='log', cmin=0.1)
```

---

## 9.6 `asinh`

```python
pl.quicklook('moment0.fits', norm_mode='asinh', stretch_a=0.05)
```

---

## 9.7 `power`

```python
pl.quicklook('moment0.fits', norm_mode='power', power_gamma=0.7)
```

---

## 9.8 表示色の最小値と最大値を固定

```python
pl.quicklook('moment0.fits', cmin=0.0, cmax=150.0)
```

---

## 9.9 `plot_map()` で `norm` を自分で渡す

```python
from astropy.visualization import ImageNormalize, AsinhStretch

m = pl.make_2d_map('cube.fits', mode='moment0')
custom_norm = ImageNormalize(m.data, stretch=AsinhStretch(a=0.03))
pl.plot_map(m, norm=custom_norm, colorbar=True)
```

---

## 9.10 colormap を変える

```python
pl.quicklook('moment0.fits', cmap='magma')
```

---

# 10. crop / cutout

## 10.1 pixel 中心と size で cutout

```python
pl.quicklook(
    'moment0.fits',
    crop={'center': (150, 160), 'size': (120, 100), 'mode': 'pixel'},
)
```

---

## 10.2 world 座標中心で cutout

```python
pl.quicklook(
    'moment0.fits',
    crop={'center': (83.8221, -5.3911), 'size': (0.20, 0.15), 'mode': 'world'},
)
```

`center` と `size` は度単位で与えます。

---

## 10.3 `crop_map2d()` を先に使う

```python
m = pl.make_2d_map('cube.fits', mode='moment0')
mc = pl.crop_map2d(m, center=(83.8221, -5.3911), size=(0.2, 0.2), mode='world')
pl.plot_map(mc)
```

---

## 10.4 crop 後に contour を重ねる

```python
base = pl.make_2d_map('12co_cube.fits', mode='moment0', vel_range=(20, 40))
base = pl.crop_map2d(base, center=(83.8221, -5.3911), size=(0.2, 0.2), mode='world')

ov = pl.make_2d_map('13co_cube.fits', mode='moment0', vel_range=(20, 40))

pl.plot_map(base, contours=[{'source': ov, 'colors': 'cyan'}])
```

---

# 11. `describe_source()` / `describe_map2d()` の使い方

## 11.1 source の説明を表示

```python
print(pl.describe_source('cube.fits', mode='moment0', vel_range=(20, 40)))
```

---

## 11.2 in-memory データの説明

```python
print(pl.describe_source(data=data3d, header=header3d, mode='moment0', vel_range=(20, 40)))
```

---

## 11.3 2D map の説明

```python
m = pl.make_2d_map('cube.fits', mode='moment0', vel_range=(20, 40))
print(pl.describe_map2d(m))
```

---

## 11.4 ログ用に説明を保存

```python
m = pl.make_2d_map('cube.fits', mode='moment0', vel_range=(20, 40))
with open('map_description.txt', 'w') as f:
    f.write(pl.describe_map2d(m))
```

---

# 12. `plot_map()` を直接使う例

## 12.1 最小例

```python
m = pl.make_2d_map('cube.fits', mode='moment0')
pl.plot_map(m)
```

---

## 12.2 colorbar を消す

```python
pl.plot_map(m, colorbar=False)
```

---

## 12.3 colorbar label を自分で指定

```python
pl.plot_map(m, colorbar_label='Integrated intensity [K km/s]')
```

---

## 12.4 grid を消す

```python
pl.plot_map(m, grid=False)
```

---

## 12.5 軸ラベルを自分で指定

```python
pl.plot_map(m, xlabel='RA (J2000)', ylabel='Dec (J2000)')
```

---

## 12.6 alpha を変える

```python
pl.plot_map(m, alpha=0.8)
```

---

## 12.7 `interpolation='bilinear'`

```python
pl.plot_map(m, interpolation='bilinear')
```

これは **表示上の補間** です。物理的な分解能は増えません。

---

## 12.8 `interpolation='bicubic'`

```python
pl.plot_map(m, interpolation='bicubic')
```

---

## 12.9 既存 Axes に描く

```python
import matplotlib.pyplot as plt

m = pl.make_2d_map('cube.fits', mode='moment0')
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=m.wcs)
pl.plot_map(m, ax=ax, show=False)
plt.show()
```

---

## 12.10 戻り値を利用する

```python
out = pl.plot_map(m, show=False)
print(out.keys())
# fig, ax, image, colorbar, contours, beam, scalebar, north_arrow, legend, description, map2d
```

---

## 12.11 `plot_map()` に OTF bundle の 2D ext を渡す

```python
bundle = pl.read_otf_bundle('otf_bundle.fits')
pl.plot_map(bundle, ext='MOMENT0', title='bundle MOMENT0')
```

`plot_map()` でも bundle を直接受けられます。

---

## 12.12 `plot_map()` に 3D bundle をそのまま渡さない

```python
bundle = pl.read_otf_bundle('otf_bundle.fits')

# 2D ext を選ぶ
pl.plot_map(bundle, ext='MOMENT0')

# あるいは 3D -> 2D 化してから渡す
m = pl.make_2d_map(bundle, mode='moment0', vel_range=(20, 40))
pl.plot_map(m)
```

`plot_map(bundle)` のように `ext` なしで 3D bundle を 2D map として描こうとするのは避けてください。

---

## 12.13 baseline viewer bundle を `plot_map()` で使う

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
viewer = pl.make_baseline_viewer_bundle(bundle, mode='signal')
m = pl.make_2d_map(viewer, mode='moment0', vel_range=(20, 40))
pl.plot_map(m, title='viewer-based moment0')
```

viewer bundle も 3D source なので、まず 2D 化してから `plot_map()` に渡します。

---

# 13. contour の基本

## 13.1 1 本だけ重ねる

```python
base = pl.make_2d_map('12co.fits', mode='moment0')
ov = pl.make_2d_map('13co.fits', mode='moment0')

pl.plot_map(
    base,
    contours=[{'source': ov, 'colors': 'cyan', 'linewidths': 1.0}],
)
```

---

## 13.2 複数 contour を重ねる

```python
pl.plot_map(
    base,
    contours=[
        {'source': '13co.fits', 'mode': 'moment0', 'colors': 'cyan', 'linewidths': 1.0},
        {'source': 'c18o.fits', 'mode': 'moment0', 'colors': 'magenta', 'linewidths': 1.0},
    ],
)
```

---

## 13.3 速度範囲を変えて contour

```python
pl.plot_map(
    base,
    contours=[
        {'source': '13co_cube.fits', 'mode': 'moment0', 'vel_range': (20, 30), 'colors': 'cyan'},
        {'source': '13co_cube.fits', 'mode': 'moment0', 'vel_range': (30, 40), 'colors': 'yellow'},
    ],
)
```

---

## 13.4 contour に label を付けて legend

```python
pl.plot_map(
    base,
    contours=[
        {'source': '13co.fits', 'mode': 'moment0', 'colors': 'cyan', 'label': '13CO'},
        {'source': 'c18o.fits', 'mode': 'moment0', 'colors': 'magenta', 'label': 'C18O'},
    ],
    legend=True,
)
```

---

## 13.5 contour source に OTF bundle を使う

```python
base = pl.make_2d_map('12co_cube.fits', mode='moment0', vel_range=(20, 40))
cont_bundle = pl.read_otf_bundle('13co_bundle.fits')

pl.plot_map(
    base,
    contours=[
        {
            'source': cont_bundle,
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'cyan',
            'levels': 'auto',
        }
    ],
)
```

---

## 13.6 color は bundle の 2D ext、contour は別 bundle の 3D cube

```python
base_bundle = pl.read_otf_bundle('12co_bundle.fits')
cont_bundle = pl.read_otf_bundle('13co_bundle.fits')

pl.plot_map(
    base_bundle,
    ext='MOMENT0',
    contours=[
        {
            'source': cont_bundle,
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'white',
            'levels': 'auto',
        }
    ],
)
```

このように、color と contour で source 型が違っていても構いません。

---

## 13.7 color は FITS、contour は bundle

```python
cont_bundle = pl.read_otf_bundle('13co_bundle.fits')

pl.plot_map(
    'optical_or_radio_moment0.fits',
    contours=[
        {
            'source': cont_bundle,
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'cyan',
        }
    ],
)
```

---

# 14. contour level の与え方

現行 v5p5 では contour level の考え方を次のように整理しています。

- **`levels='auto'` または `level_mode='auto'`**  
  旧版互換。**正のピーク値の 10, 30, 50, 70, 90%** を使います。
- **`level_mode='fraction'`**  
  `fraction_levels=[...]` をピーク値に掛けます。
- **`level_mode='manual'`**  
  `levels=[...]` を factor として扱い、必要なら `level_scale` を掛けます。  
  例: `levels=[1,2,3,4], level_scale=3` → `3,6,9,12`
- **`level_mode='rms'`**  
  `rms * sigma_scale * sigma_levels` で level を作ります。
- **`level_mode='explicit'` / `level_mode='sigma'`**  
  後方互換 alias として受理されます。内部的にはそれぞれ `manual` / `rms` に正規化されます。

以下では、現行仕様に沿って例を並べます。

## 14.1 旧版互換の `levels='auto'`

```python
pl.plot_map(
    base,
    contours=[{'source': ov, 'levels': 'auto', 'colors': 'white'}],
)
```

`auto` は現行版では **legacy fraction** です。従来どおり「ピークの 10, 30, 50, 70, 90%」を意味します。

---

## 14.2 `level_mode='auto'` を明示する

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'auto',
        'colors': 'white',
    }],
)
```

`levels='auto'` とほぼ同じですが、layer 辞書側の意図を明示したい時はこちらでも構いません。

---

## 14.3 fraction mode

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'fraction',
        'fraction_levels': [0.2, 0.4, 0.6, 0.8],
        'colors': 'yellow',
    }],
)
```

これは「ピークの 20, 40, 60, 80%」です。  
**fraction mode では `levels=[...]` ではなく `fraction_levels=[...]` を使う**、というのが現行版の整理です。

---

## 14.4 fraction mode を `quicklook()` で使う

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(0, 20),
    contours=[{
        'level_mode': 'fraction',
        'fraction_levels': [0.1, 0.25, 0.5, 0.75],
        'colors': 'white',
    }],
)
```

---

## 14.5 manual mode（そのまま値を与える）

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'manual',
        'levels': [5, 10, 20, 40],
        'colors': 'white',
    }],
)
```

この場合、level は `[5, 10, 20, 40]` がそのまま使われます。

---

## 14.6 manual mode + `level_scale`

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'manual',
        'levels': [1, 2, 3, 4],
        'level_scale': 3,
        'colors': 'white',
    }],
)
```

これは **`3*[1,2,3,4]` 的な発想を安全に表現する方法**です。  
実際の contour level は `3, 6, 9, 12` になります。

---

## 14.7 manual mode を `quicklook()` で使う

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(0, 20),
    contours=[{
        'level_mode': 'manual',
        'levels': [1, 2, 3, 4],
        'level_scale': 5,
        'colors': 'cyan',
    }],
)
```

この場合、実際の level は `5, 10, 15, 20` です。

---

## 14.8 `level_mode='explicit'`（後方互換 alias）

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'explicit',
        'levels': [3, 6, 9],
        'colors': 'white',
    }],
)
```

`explicit` は現行版では **`manual` の alias** と考えて構いません。

---

## 14.9 RMS mode（RMS 明示）

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'rms',
        'rms': 0.35,
        'sigma_levels': [3, 5, 7],
        'colors': 'cyan',
    }],
)
```

この場合、実際の level は `0.35 × [3, 5, 7]` です。

---

## 14.10 RMS mode + `sigma_scale`

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'rms',
        'rms': 0.8,
        'sigma_levels': [1, 2, 3, 4],
        'sigma_scale': 3,
        'colors': 'cyan',
    }],
)
```

この場合、実際の level は `0.8 × 3 × [1, 2, 3, 4]` なので `2.4, 4.8, 7.2, 9.6` です。  
**`sigma_scale` は「3σ, 6σ, 9σ, 12σ」などを簡潔に書きたい時」に便利**です。

---

## 14.11 `level_mode='sigma'`（後方互換 alias）

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'sigma',
        'rms': 0.25,
        'sigma_levels': [3, 5, 7],
        'colors': 'white',
    }],
)
```

`level_mode='sigma'` は現行版では **`rms` の alias** です。新しいコードでは `level_mode='rms'` と書く方が分かりやすいです。

---

## 14.12 RMS を自動推定に頼る例（補助機能）

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'rms',
        'sigma_levels': [3, 5, 7],
        'colors': 'cyan',
    }],
)
```

RMS を与えなければ内部推定が使われますが、moment0 や signal の広い地図では **物理的に最適とは限らない** ため、研究用途では `rms` を自分で与える、または fraction / manual を使う方が安全です。

---

## 14.13 symmetric=True で負側も描く

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'rms',
        'rms': 0.25,
        'sigma_levels': [3, 5, 7],
        'symmetric': True,
        'colors': 'white',
        'linestyles': 'solid',
    }],
)
```

---

## 14.14 negative_sigma_levels を自分で指定

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'rms',
        'rms': 0.25,
        'sigma_levels': [3, 5, 7],
        'negative_sigma_levels': [3],
        'colors': 'white',
    }],
)
```

---

## 14.15 `levels=[...]` を非 manual mode で渡すときの注意

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'fraction',
        'levels': [5, 10, 20],
        'colors': 'white',
    }],
)
```

これは**後方互換のために受理されます**が、現行版の整理としては推奨しません。  
`fraction` なら `fraction_levels=[...]`、`manual` なら `levels=[...]`、`rms` なら `sigma_levels=[...]` を使う方が意図が明確です。

---

## 14.16 同じ base に fraction と manual の 2 系統を重ねる

```python
pl.plot_map(
    base,
    contours=[
        {
            'source': ov,
            'level_mode': 'fraction',
            'fraction_levels': [0.2, 0.4, 0.6],
            'colors': 'yellow',
            'label': 'fraction',
        },
        {
            'source': ov,
            'level_mode': 'manual',
            'levels': [1, 2, 3, 4],
            'level_scale': 5,
            'colors': 'cyan',
            'label': 'manual x5',
        },
    ],
    legend=True,
)
```

---

## 14.17 `quicklook()` で旧版互換 auto contour を最短で使う

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(0, 20),
    contours=[{'levels': 'auto', 'colors': 'white'}],
)
```

これは現行版でも、**従来の「ピークの 10, 30, 50, 70, 90%」** と同じ意味です。

---
# 15. contour の見た目を変える

## 15.1 線幅

```python
pl.plot_map(
    base,
    contours=[{'source': ov, 'colors': 'white', 'linewidths': 1.8}],
)
```

---

## 15.2 線種

```python
pl.plot_map(
    base,
    contours=[{'source': ov, 'colors': 'white', 'linestyles': 'dashed'}],
)
```

---

## 15.3 アルファ

```python
pl.plot_map(
    base,
    contours=[{'source': ov, 'colors': 'white', 'alpha': 0.7}],
)
```

---

## 15.4 contour だけ別 colormap 的に色分けしたい時

```python
colors = ['cyan', 'lime', 'yellow', 'orange']
levels = [5, 10, 20, 40]

for c, lv in zip(colors, levels):
    pl.plot_map(
        base,
        contours=[{'source': ov, 'levels': [lv], 'colors': c}],
        show=False,
    )
```

通常は 1 回の `contours` にまとめる方がよいですが、単一 level ごとに色を変えたい時の考え方として使えます。

---

# 16. 光学画像の上に電波 contour

## 16.1 最も基本的な例

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits', 'cmap': 'gray'},
    overlays=[
        {'kind': 'contour', 'source': '12co_moment0.fits', 'colors': 'cyan', 'label': '12CO'},
    ],
    legend=True,
)
```

この使い方が、一般には最も安全です。

---

## 16.2 13CO と C18O も追加

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits', 'cmap': 'gray'},
    overlays=[
        {'kind': 'contour', 'source': '12co.fits', 'colors': 'cyan', 'label': '12CO'},
        {'kind': 'contour', 'source': '13co.fits', 'colors': 'yellow', 'label': '13CO'},
        {'kind': 'contour', 'source': 'c18o.fits', 'colors': 'magenta', 'label': 'C18O'},
    ],
    legend=True,
)
```

---

## 16.3 光学は薄い grayscale にする

```python
pl.plot_scene(
    base={
        'kind': 'image',
        'source': 'optical.fits',
        'cmap': 'gray',
        'norm_percentile': (5, 99.5),
    },
    overlays=[
        {'kind': 'contour', 'source': '12co.fits', 'colors': 'red'},
    ],
)
```

---

## 16.4 光学画像の上に OTF bundle contour

```python
co_bundle = pl.read_otf_bundle('12co_bundle.fits')

pl.plot_scene(
    base={
        'kind': 'image',
        'source': 'optical.fits',
        'cmap': 'gray',
    },
    overlays=[
        {
            'kind': 'contour',
            'source': co_bundle,
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'cyan',
        }
    ],
)
```

OTF bundle を contour source としてそのまま使えます。

---

## 16.5 光学画像の上に baseline viewer bundle 由来の contour

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
viewer = pl.make_baseline_viewer_bundle(bundle, mode='signal')

pl.plot_scene(
    base={
        'kind': 'image',
        'source': 'optical.fits',
        'cmap': 'gray',
    },
    overlays=[
        {
            'kind': 'contour',
            'source': viewer,
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'yellow',
        }
    ],
)
```

viewer bundle も contour source にできます。

---

# 17. image overlay と reprojection

## 17.1 透過 image overlay（同一 WCS 近傍）

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits', 'cmap': 'gray'},
    overlays=[
        {
            'kind': 'image',
            'source': '12co_moment0.fits',
            'cmap': 'inferno',
            'alpha': 0.45,
            'reproject': 'auto',
        }
    ],
)
```

---

## 17.2 reprojection を必須にする

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits', 'cmap': 'gray'},
    overlays=[
        {
            'kind': 'image',
            'source': '12co_moment0.fits',
            'alpha': 0.5,
            'reproject': 'required',
        }
    ],
)
```

`reproject` パッケージが無ければ error になります。

---

## 17.3 reprojection を禁止

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits', 'cmap': 'gray'},
    overlays=[
        {
            'kind': 'image',
            'source': '12co_moment0.fits',
            'alpha': 0.5,
            'reproject': 'never',
        }
    ],
)
```

---

## 17.4 image overlay の色レンジを個別指定

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits', 'cmap': 'gray'},
    overlays=[
        {
            'kind': 'image',
            'source': '12co_moment0.fits',
            'cmap': 'magma',
            'alpha': 0.6,
            'cmin': 0.0,
            'cmax': 80.0,
            'reproject': 'auto',
        }
    ],
)
```

---

# 18. `plot_scene()` の基本形

## 18.1 裸の source を base にする

```python
pl.plot_scene('moment0.fits')
```

これは内部で

```python
{'kind': 'image', 'source': 'moment0.fits'}
```

と解釈されます。

---

## 18.2 base を辞書で明示する

```python
pl.plot_scene({
    'kind': 'image',
    'source': 'moment0.fits',
    'cmap': 'viridis',
    'beam': 'auto',
})
```

---

## 18.3 base だけ、overlay なし

```python
pl.plot_scene(
    {'kind': 'image', 'source': 'moment0.fits', 'title': 'base only'},
)
```

---

## 18.4 overlays を後から足す

```python
scene = {
    'kind': 'image',
    'source': 'moment0.fits',
    'cmap': 'gray',
}

ov = [
    {'kind': 'beam', 'label': 'beam'},
    {'kind': 'scalebar'},
    {'kind': 'north_arrow'},
]

pl.plot_scene(scene, overlays=ov)
```

---

## 18.5 base に OTF bundle の 2D ext を使う

```python
bundle = pl.read_otf_bundle('otf_bundle.fits')

pl.plot_scene(
    base={
        'kind': 'image',
        'source': bundle,
        'ext': 'MOMENT0',
        'cmap': 'inferno',
        'beam': 'auto',
    },
    title='bundle MOMENT0 as base',
)
```

---

## 18.6 overlay contour に OTF bundle を使う

```python
base_bundle = pl.read_otf_bundle('12co_bundle.fits')
cont_bundle = pl.read_otf_bundle('13co_bundle.fits')

pl.plot_scene(
    base={
        'kind': 'image',
        'source': base_bundle,
        'ext': 'MOMENT0',
        'cmap': 'magma',
    },
    overlays=[
        {
            'kind': 'contour',
            'source': cont_bundle,
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'white',
        }
    ],
)
```

---

## 18.7 base に baseline viewer bundle を使う

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
viewer = pl.make_baseline_viewer_bundle(bundle, mode='signal')

pl.plot_scene(
    base={
        'kind': 'image',
        'source': viewer,
        'mode': 'moment0',
        'vel_range': (20, 40),
        'cmap': 'inferno',
    },
    title='viewer bundle as base source',
)
```

`plot_scene()` でも viewer bundle は通常の 3D source として扱います。

---

# 19. marker overlay

## 19.1 pixel 座標でマーカー

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[
        {'kind': 'marker', 'xy': [(100, 120), (140, 180)], 'color': 'red', 'marker': 'x', 's': 60},
    ],
)
```

---

## 19.2 world 座標（度）でマーカー

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[
        {
            'kind': 'marker',
            'coords': [(83.8221, -5.3911), (83.8000, -5.3700)],
            'frame': 'icrs',
            'color': 'cyan',
            'marker': '+',
        }
    ],
)
```

---

## 19.3 `lon=` と `lat=` を使う

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[
        {
            'kind': 'marker',
            'lon': [83.8221, 83.8000],
            'lat': [-5.3911, -5.3700],
            'frame': 'icrs',
            'color': 'yellow',
        }
    ],
)
```

---

## 19.4 `SkyCoord` を使う

```python
from astropy.coordinates import SkyCoord
import astropy.units as u

sc = SkyCoord([83.8221, 83.8000]*u.deg, [-5.3911, -5.3700]*u.deg, frame='icrs')

pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'marker', 'skycoord': sc, 'color': 'lime', 'marker': 'o', 's': 30}],
)
```

---

## 19.5 galactic 座標の marker

```python
pl.plot_scene(
    'galactic_map.fits',
    overlays=[
        {
            'kind': 'marker',
            'coords': [(210.1, -1.3), (210.4, -1.1)],
            'frame': 'galactic',
            'color': 'white',
        }
    ],
)
```

---

# 20. catalog overlay の基本

## 20.1 dict of arrays を使う

```python
catalog = {
    'ra': [83.8221, 83.8000],
    'dec': [-5.3911, -5.3700],
    'name': ['Src-A', 'Src-B'],
}

pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'catalog', 'table': catalog, 'color': 'cyan', 'marker': 'o'}],
)
```

---

## 20.2 `labels` を明示する

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[
        {
            'kind': 'catalog',
            'coords': [(83.8221, -5.3911), (83.8000, -5.3700)],
            'frame': 'icrs',
            'labels': ['A', 'B'],
            'color': 'yellow',
        }
    ],
)
```

---

## 20.3 pixel 座標の catalog

```python
catalog = {
    'x': [100, 160, 210],
    'y': [120, 175, 240],
    'label': ['P1', 'P2', 'P3'],
}

pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'catalog', 'table': catalog, 'coord_mode': 'pixel', 'color': 'white'}],
)
```

---

## 20.4 `pandas.DataFrame`

```python
import pandas as pd

df = pd.DataFrame({
    'ra': [83.8221, 83.8000],
    'dec': [-5.3911, -5.3700],
    'name': ['core1', 'core2'],
})

pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'catalog', 'table': df, 'label_col': 'name', 'color': 'red'}],
)
```

---

## 20.5 `astropy.table.Table`

```python
from astropy.table import Table

tab = Table({
    'glon': [210.1, 210.3],
    'glat': [-1.3, -1.2],
    'source': ['G1', 'G2'],
})

pl.plot_scene(
    'galactic_map.fits',
    overlays=[{'kind': 'catalog', 'table': tab, 'color': 'orange'}],
)
```

---

## 20.6 カラム名が標準でない時

```python
df = pd.DataFrame({
    'ra_deg_custom': [83.8221, 83.8000],
    'dec_deg_custom': [-5.3911, -5.3700],
    'srcname': ['A', 'B'],
})

pl.plot_scene(
    'moment0.fits',
    overlays=[{
        'kind': 'catalog',
        'table': df,
        'lon_col': 'ra_deg_custom',
        'lat_col': 'dec_deg_custom',
        'label_col': 'srcname',
        'frame': 'icrs',
    }],
)
```

---

## 20.7 catalog ラベルの位置を調整

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{
        'kind': 'catalog',
        'table': df,
        'label_col': 'name',
        'label_dx': 10,
        'label_dy': -8,
        'label_offset_mode': 'pixel',
        'text_color': 'white',
        'fontsize': 8,
    }],
)
```

---

## 20.8 world offset でラベルをずらす

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{
        'kind': 'catalog',
        'table': df,
        'label_col': 'name',
        'label_dx': 0.002,
        'label_dy': 0.001,
        'label_offset_mode': 'world',
        'frame': 'icrs',
    }],
)
```

world offset は度単位相当で効くため、通常は `pixel` の方が扱いやすいです。

---

# 21. text overlay

## 21.1 pixel 座標に文字を書く

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'text', 'xy': (120, 150), 'text': 'Core A', 'color': 'white'}],
)
```

---

## 21.2 world 座標に文字を書く

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{
        'kind': 'text',
        'coord': (83.8221, -5.3911),
        'frame': 'icrs',
        'text': 'IRc2',
        'color': 'yellow',
        'fontsize': 12,
    }],
)
```

---

## 21.3 `SkyCoord` で文字

```python
from astropy.coordinates import SkyCoord
import astropy.units as u

sc = SkyCoord(83.8221*u.deg, -5.3911*u.deg, frame='icrs')

pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'text', 'skycoord': sc, 'text': 'Peak', 'color': 'cyan'}],
)
```

---

# 22. beam / scalebar / north_arrow

## 22.1 base に beam を付ける

```python
pl.plot_scene({
    'kind': 'image',
    'source': 'moment0.fits',
    'beam': 'auto',
})
```

---

## 22.2 beam を overlay として追加

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'beam', 'beam': 'auto', 'label': 'Beam'}],
    legend=True,
)
```

---

## 22.3 scalebar を既定設定で付ける

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'scalebar'}],
)
```

---

## 22.4 scalebar の長さを指定

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'scalebar', 'config': {'length': 30, 'units': 'arcsec'}}],
)
```

---

## 22.5 arcmin 単位の scalebar

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'scalebar', 'config': {'length': 1.0, 'units': 'arcmin'}}],
)
```

---

## 22.6 pixel 単位の scalebar

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'scalebar', 'config': {'length': 40, 'units': 'pix'}}],
)
```

---

## 22.7 scalebar の位置と色を変える

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{
        'kind': 'scalebar',
        'config': {
            'length': 30,
            'units': 'arcsec',
            'location': 'lower left',
            'color': 'yellow',
            'fontsize': 10,
        }
    }],
)
```

---

## 22.8 north_arrow を付ける

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{'kind': 'north_arrow'}],
)
```

---

## 22.9 north_arrow の長さと位置を変える

```python
pl.plot_scene(
    'moment0.fits',
    overlays=[{
        'kind': 'north_arrow',
        'config': {
            'length': 20,
            'units': 'arcsec',
            'location': 'upper right',
            'color': 'white',
        }
    }],
)
```

---

## 22.10 publication 用の定番 annotation 一式

```python
pl.plot_scene(
    {'kind': 'image', 'source': 'moment0.fits', 'beam': 'auto'},
    overlays=[
        {'kind': 'scalebar', 'config': {'length': 30, 'units': 'arcsec'}},
        {'kind': 'north_arrow'},
    ],
    title='12CO integrated intensity',
)
```

---

# 23. RGB 合成の基本

## 23.1 最小例

```python
pl.plot_rgb(red='r.fits', green='g.fits', blue='b.fits')
```

---

## 23.2 3 本の moment0 を RGB にする

```python
pl.plot_rgb(
    red={'source': '12co_cube.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
    green={'source': '13co_cube.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
    blue={'source': 'c18o_cube.fits', 'mode': 'moment0', 'vel_range': (20, 40)},
)
```

---

## 23.3 各 channel の normalize を変える

```python
pl.plot_rgb(
    red='12co.fits',
    green='13co.fits',
    blue='c18o.fits',
    red_norm={'norm_mode': 'asinh', 'percentile': (1, 99.7)},
    green_norm={'norm_mode': 'sqrt', 'percentile': (1, 99.5)},
    blue_norm={'norm_mode': 'linear', 'percentile': (5, 99.0)},
)
```

---

## 23.4 RGB を crop する

```python
pl.plot_rgb(
    red='12co.fits',
    green='13co.fits',
    blue='c18o.fits',
    crop={'center': (83.8221, -5.3911), 'size': (0.2, 0.2), 'mode': 'world'},
)
```

---

## 23.5 title, grid, save

```python
pl.plot_rgb(
    red='12co.fits',
    green='13co.fits',
    blue='c18o.fits',
    title='RGB molecular-line composite',
    save='rgb_lines.pdf',
    show=False,
)
```

---

## 23.6 scalebar と north_arrow を付ける

```python
pl.plot_rgb(
    red='12co.fits',
    green='13co.fits',
    blue='c18o.fits',
    scalebar={'length': 30, 'units': 'arcsec'},
    north_arrow=True,
)
```

---

# 24. RGB を base にして plot_scene を使う

## 24.1 RGB base + contour

```python
pl.plot_scene(
    base={
        'kind': 'rgb',
        'red': '12co.fits',
        'green': '13co.fits',
        'blue': 'c18o.fits',
    },
    overlays=[
        {'kind': 'contour', 'source': 'dust_continuum.fits', 'colors': 'white', 'label': 'Dust'},
    ],
    legend=True,
)
```

---

## 24.2 RGB base + catalog + beam

```python
pl.plot_scene(
    base={
        'kind': 'rgb',
        'red': '12co.fits',
        'green': '13co.fits',
        'blue': 'c18o.fits',
        'scalebar': {'length': 30, 'units': 'arcsec'},
        'north_arrow': True,
    },
    overlays=[
        {'kind': 'catalog', 'table': df, 'label_col': 'name', 'color': 'white'},
        {'kind': 'beam', 'beam': 'auto'},
    ],
)
```

---

# 25. base + overlays の典型レシピ

## 25.1 12CO base + 13CO/C18O contour + catalog

```python
pl.plot_scene(
    base={
        'kind': 'image',
        'source': '12co_cube.fits',
        'mode': 'moment0',
        'vel_range': (20, 40),
        'cmap': 'inferno',
        'beam': 'auto',
        'colorbar_label': '12CO [K km/s]',
    },
    overlays=[
        {
            'kind': 'contour',
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'cyan',
            'level_mode': 'sigma',
            'sigma_levels': [3, 5, 7],
            'label': '13CO',
        },
        {
            'kind': 'contour',
            'source': 'c18o_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'lime',
            'level_mode': 'sigma',
            'sigma_levels': [3, 5],
            'label': 'C18O',
        },
        {
            'kind': 'catalog',
            'table': df,
            'label_col': 'name',
            'color': 'white',
            'marker': 'x',
        },
        {'kind': 'scalebar', 'config': {'length': 30, 'units': 'arcsec'}},
        {'kind': 'north_arrow'},
    ],
    legend=True,
    title='Line comparison',
)
```

---

## 25.2 光学 base + 電波 pseudo-color + 電波 contour

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'optical.fits', 'cmap': 'gray'},
    overlays=[
        {
            'kind': 'image',
            'source': '12co.fits',
            'cmap': 'inferno',
            'alpha': 0.35,
            'reproject': 'auto',
            'label': '12CO map',
        },
        {
            'kind': 'contour',
            'source': '13co.fits',
            'colors': 'cyan',
            'linewidths': 1.0,
            'label': '13CO contour',
        },
    ],
    legend=True,
)
```

---

## 25.3 continuum base + line contours

```python
pl.plot_scene(
    base={'kind': 'image', 'source': 'continuum.fits', 'cmap': 'magma', 'beam': 'auto'},
    overlays=[
        {'kind': 'contour', 'source': '12co.fits', 'colors': 'white', 'label': '12CO'},
        {'kind': 'contour', 'source': '13co.fits', 'colors': 'cyan', 'label': '13CO'},
    ],
    legend=True,
)
```

---

# 26. 座標系混在の例

## 26.1 ICRS base に galactic catalog を重ねる

```python
cat = {
    'glon': [209.9, 210.1],
    'glat': [-1.2, -1.3],
    'name': ['G1', 'G2'],
}

pl.plot_scene(
    base={'kind': 'image', 'source': 'icrs_map.fits'},
    overlays=[{'kind': 'catalog', 'table': cat, 'color': 'yellow'}],
)
```

`frame` はカラム名から自動推定されます。

---

## 26.2 明示的に frame を与える

```python
pl.plot_scene(
    'icrs_map.fits',
    overlays=[{
        'kind': 'catalog',
        'coords': [(209.9, -1.2), (210.1, -1.3)],
        'frame': 'galactic',
        'labels': ['G1', 'G2'],
        'color': 'cyan',
    }],
)
```

---

# 27. 返り値を活かす例

## 27.1 `plot_scene()` の返り値を確認

```python
out = pl.plot_scene('moment0.fits', show=False)
print(out.keys())
# fig, ax, base, base_description, overlay_artists, overlay_maps, contours, beam_artists, legend
```

---

## 27.2 返り値から追加の matplotlib 調整

```python
out = pl.plot_scene('moment0.fits', show=False)
out['ax'].set_title('Edited after plotting')
out['fig'].savefig('edited_after_plotting.pdf', dpi=300)
```

---

## 27.3 `base_description` をログに残す

```python
out = pl.plot_scene('moment0.fits', show=False)
print(out['base_description'])
```

---

# 28. サブルーチンとしての再利用

## 28.1 自分用の quick wrapper を作る

```python
def show_line_map(path, v0, v1, title=None):
    return pl.quicklook(
        path,
        mode='moment0',
        vel_range=(v0, v1),
        cmap='inferno',
        beam='auto',
        title=title,
    )

show_line_map('12co_cube.fits', 20, 40, title='12CO 20-40 km/s')
```

---

## 28.2 同じ style を複数 source に適用

```python
def standard_plot(m):
    return pl.plot_map(
        m,
        cmap='magma',
        norm_mode='asinh',
        norm_percentile=(1, 99.7),
        beam='auto',
        show=False,
    )

for path in ['map1.fits', 'map2.fits', 'map3.fits']:
    m = pl.make_2d_map(path, mode='moment0')
    out = standard_plot(m)
    out['fig'].savefig(path.replace('.fits', '.pdf'), dpi=300)
```

---

# 29. ループ処理の例

## 29.1 複数 velocity bin を順番に保存

```python
bins = [(10, 20), (20, 30), (30, 40), (40, 50)]
for v0, v1 in bins:
    pl.quicklook(
        'cube.fits',
        mode='moment0',
        vel_range=(v0, v1),
        title=f'{v0}-{v1} km/s',
        save=f'moment0_{v0}_{v1}.pdf',
        show=False,
    )
```

---

## 29.2 3 本の線を同じ velocity bin で比較

```python
bins = [(20, 30), (30, 40)]
for v0, v1 in bins:
    pl.plot_scene(
        base={'kind': 'image', 'source': '12co_cube.fits', 'mode': 'moment0', 'vel_range': (v0, v1)},
        overlays=[
            {'kind': 'contour', 'source': '13co_cube.fits', 'mode': 'moment0', 'vel_range': (v0, v1), 'colors': 'cyan'},
            {'kind': 'contour', 'source': 'c18o_cube.fits', 'mode': 'moment0', 'vel_range': (v0, v1), 'colors': 'lime'},
        ],
        title=f'{v0}-{v1} km/s',
        save=f'compare_{v0}_{v1}.png',
        show=False,
    )
```

---

# 30. publication 用の保存レシピ

## 30.1 PDF

```python
pl.plot_scene(
    'moment0.fits',
    save='figure.pdf',
    show=False,
)
```

---

## 30.2 高 DPI PNG

```python
pl.plot_scene(
    'moment0.fits',
    save='figure.png',
    save_dpi=400,
    show=False,
)
```

---

## 30.3 透明背景

```python
pl.plot_scene(
    'moment0.fits',
    save='figure_transparent.png',
    save_transparent=True,
    show=False,
)
```

---

## 30.4 最終論文図の定番

```python
pl.plot_scene(
    base={
        'kind': 'image',
        'source': 'moment0.fits',
        'cmap': 'inferno',
        'beam': 'auto',
    },
    overlays=[
        {'kind': 'scalebar', 'config': {'length': 30, 'units': 'arcsec'}},
        {'kind': 'north_arrow'},
    ],
    title='Integrated intensity',
    save='paper_figure.pdf',
    show=False,
)
```

---

# 31. CLI cookbook

以下では、ファイル名は適宜読み替えてください。

## 31.1 最小例

```bash
python plotting_updated_v5p5.py cube.fits
```

---

## 31.2 extension を指定

```bash
python plotting_updated_v5p5.py products.fits --ext MOMENT0
```

---

## 31.3 velocity 範囲付き moment0

```bash
python plotting_updated_v5p5.py cube.fits --mode moment0 --vel-range 20,40
```

---

## 31.4 channel 範囲付き moment0

```bash
python plotting_updated_v5p5.py cube.fits --mode moment0 --chan-range 100,180
```

---

## 31.5 smoothing

```bash
python plotting_updated_v5p5.py cube.fits --mode moment0 --vel-range 20,40 --target-hpbw-arcsec 30
```

---

## 31.6 colormap と norm

```bash
python plotting_updated_v5p5.py cube.fits --cmap magma --norm-mode asinh --norm-percentile 1,99.7
```

---

## 31.7 cmin/cmax を固定

```bash
python plotting_updated_v5p5.py cube.fits --cmin 0 --cmax 100
```

---

## 31.8 grid を消す

```bash
python plotting_updated_v5p5.py cube.fits --no-grid
```

---

## 31.9 contour を最短で重ねる（現行版では `--contour-levels` 必須ではない）

```bash
python plotting_updated_v5p5.py cube.fits --contour-color white
```

現行 v5p5 では、CLI は **contour 系オプションが 1 つでも与えられれば contour を作る**ようになっています。  
したがって `--contour-levels auto` を必ず付ける必要はありません。

---

## 31.10 旧版互換 auto contour

```bash
python plotting_updated_v5p5.py cube.fits --contour-levels auto
```

これは **ピークの 10, 30, 50, 70, 90%** です。

---

## 31.11 fraction contour

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode fraction --contour-fractions 0.2,0.4,0.6,0.8
```

---

## 31.12 manual contour

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode manual --contour-levels 5,10,20,40
```

---

## 31.13 manual contour + scale

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode manual --contour-levels 1,2,3,4 --contour-level-scale 3
```

これは `3, 6, 9, 12` を意味します。

---

## 31.14 explicit contour（後方互換 alias）

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode explicit --contour-levels 3,6,9
```

---

## 31.15 RMS contour

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode rms --contour-rms 0.25 --contour-sigmas 3,5,7
```

---

## 31.16 sigma contour（後方互換 alias）

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode sigma --contour-rms 0.25 --contour-sigmas 3,5,7
```

新しいコードやメモでは `rms` と書く方が分かりやすいです。

---

## 31.17 RMS contour + sigma_scale

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode rms --contour-rms 0.8 --contour-sigmas 1,2,3,4 --contour-sigma-scale 3
```

これは `2.4, 4.8, 7.2, 9.6` を意味します。

---

## 31.18 負側 contour を加える

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode rms --contour-rms 0.25 --contour-sigmas 3,5,7 --contour-negative
```

---

## 31.19 contour 色と線幅も同時指定

```bash
python plotting_updated_v5p5.py cube.fits --contour-level-mode fraction --contour-fractions 0.2,0.4,0.6 --contour-color cyan --contour-linewidth 1.2
```

---

## 31.20 beam / scalebar / north_arrow

```bash
python plotting_updated_v5p5.py cube.fits --beam auto --scalebar --scalebar-length 30 --north-arrow
```

---

## 31.21 説明を出力

```bash
python plotting_updated_v5p5.py cube.fits --describe
```

---

## 31.22 API help だけ見る

```bash
python plotting_updated_v5p5.py --api-help
```

---

## 31.23 publication 用に PDF 保存

```bash
python plotting_updated_v5p5.py cube.fits --mode moment0 --vel-range 20,40 --beam auto --scalebar --output moment0.pdf --no-grid
```

---

## 31.24 contour を使うなら Python API と CLI のどちらがよいか

- **試しに 1 枚出したい**なら CLI で十分です。
- **複数 contour source、catalog、annotation まで進めたい**なら Python API の `plot_scene()` が自然です。
- **base + contour を最短で試したい**なら Python API の `quicklook(contours=...)` が最も簡単です。


---

## 31.25 OTF bundle FITS を CLI で使う

```bash
python -m sd_radio_spectral_fits.map_3d.plotting otf_bundle.fits --ext MOMENT0
```

ディスク上の bundle FITS を products として読むだけなら、CLI でも自然に扱えます。

---

## 31.26 Python 上の `OTFBundle` object は CLI には直接渡せない

CLI は path ベースなので、Python 上で作った

- `OTFBundle`
- baseline viewer bundle

をそのまま引数には渡せません。

その場合は Python API を使います。

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
viewer = pl.make_baseline_viewer_bundle(bundle, mode='signal')
pl.quicklook(viewer, mode='moment0', vel_range=(20, 40))
```

---

# 32. `quicklook(contours=...)` の集中的な実例集

この章は、**現行 v5p5 での重要追加機能**である `quicklook(contours=...)` を中心に、
「最短で 1 枚出したいが contour も欲しい」状況に特化して例を増やした章です。

基本方針:

- まずは `quicklook()` を使う
- それで足りなくなったら `plot_scene()` へ進む
- `plot_map()` は互換 API として残っているが、新規コードでは主役にしない

---

## 32.1 base と同じ source を contour に使う最短例

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    contours=[{'levels': 'auto', 'colors': 'white'}],
)
```

この書き方では、contour 側の `source`, `mode`, `vel_range`, `spectral_unit`, `ext` などは、
**base の `quicklook()` 側の設定を自動継承**します。

---

## 32.2 base と同じ source だが contour だけ色を変える

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    cmap='inferno',
    contours=[
        {'levels': 'auto', 'colors': 'cyan', 'linewidths': 1.2},
    ],
)
```

---

## 32.3 base と同じ source だが contour だけ別の mode を使う

```python
pl.quicklook(
    'cube.fits',
    mode='channel_mean',
    vel_range=(20, 40),
    contours=[
        {
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'white',
        }
    ],
)
```

このように、base と contour で mode を変えることもできます。

---

## 32.4 base と contour で source を分ける

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    cmap='inferno',
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'cyan',
            'levels': 'auto',
        }
    ],
)
```

---

## 32.5 3 本の線を `quicklook()` だけで比較する

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    cmap='magma',
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'cyan',
            'levels': 'auto',
            'label': '13CO',
        },
        {
            'source': 'c18o_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'lime',
            'levels': 'auto',
            'label': 'C18O',
        },
    ],
    legend=True,
)
```

`quicklook()` でも、ここまでなら十分実用的です。

---

## 32.6 crop と contour を一緒に使う

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    crop={'center': (83.8221, -5.3911), 'size': (0.20, 0.15), 'mode': 'world'},
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'cyan',
            'levels': 'auto',
        }
    ],
)
```

現行版では、`quicklook(contours=...)` でも **crop が contour 側に反映**されます。

---

## 32.7 `target_hpbw_arcsec` と contour を一緒に使う

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    target_hpbw_arcsec=45.0,
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'target_hpbw_arcsec': 45.0,
            'colors': 'white',
            'levels': 'auto',
        }
    ],
)
```

重要:

- 現行版では、`quicklook()` 内で contour 側の smoothing が **二重適用されない** よう調整されています。
- したがって、`target_hpbw_arcsec` は base と contour に対して **それぞれ 1 回だけ**効きます。

---

## 32.8 `orig_hpbw_arcsec` を明示した beam matching

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    orig_hpbw_arcsec=17.0,
    target_hpbw_arcsec=30.0,
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'orig_hpbw_arcsec': 18.0,
            'target_hpbw_arcsec': 30.0,
            'colors': 'cyan',
            'levels': 'auto',
        }
    ],
)
```

複数線で beam をそろえたいときの基本形です。

---

## 32.9 2D moment0 を base にして、3D cube から contour を作る

```python
pl.quicklook(
    '12co_moment0.fits',
    cmap='inferno',
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'colors': 'white',
            'levels': 'auto',
        }
    ],
)
```

---

## 32.10 `quicklook()` で title / describe / help を同時に使う

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    contours=[{'levels': 'auto', 'colors': 'white'}],
    title='12CO with contour',
    describe=True,
)
```

`help=True` は通常それ単独で使う方が分かりやすいですが、
試行錯誤の途中で `describe=True` と併用すると便利です。

---

## 32.11 `quicklook()` と `plot_scene()` の境界

`quicklook()` で快適に扱えるのは、概ね次の範囲です。

- base image 1 枚
- contour 数本
- beam
- crop
- 保存
- 説明表示

一方、次が欲しくなったら `plot_scene()` に移る方が自然です。

- marker / catalog / text
- optical base + radio image overlay + radio contour
- scalebar / north_arrow / legend をまとめて管理
- 複数の image overlay
- base を RGB にする

---

## 32.12 base を OTF bundle の 2D ext にする

```python
base_bundle = pl.read_otf_bundle('12co_bundle.fits')

pl.quicklook(
    base_bundle,
    ext='MOMENT0',
    cmap='inferno',
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'white',
        }
    ],
)
```

---

## 32.13 base と contour の両方に bundle を使う

```python
base_bundle = pl.read_otf_bundle('12co_bundle.fits')
cont_bundle = pl.read_otf_bundle('13co_bundle.fits')

pl.quicklook(
    base_bundle,
    ext='MOMENT0',
    cmap='magma',
    contours=[
        {
            'source': cont_bundle,
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'cyan',
        }
    ],
)
```

---

## 32.14 baseline viewer bundle を base にする

```python
bundle = pl.read_otf_bundle('baseline_bundle.fits')
viewer = pl.make_baseline_viewer_bundle(bundle, mode='signal')

pl.quicklook(
    viewer,
    mode='moment0',
    vel_range=(20, 40),
    cmap='inferno',
    contours=[
        {
            'source': bundle,
            'ext': 'MOMENT0',
            'levels': 'auto',
            'colors': 'white',
        }
    ],
)
```

viewer bundle と元 bundle を比較したいときの基本形です。

---

## 32.15 contour に bundle の 2D ext を直接使う

```python
base_bundle = pl.read_otf_bundle('12co_bundle.fits')
cont_bundle = pl.read_otf_bundle('13co_bundle.fits')

pl.quicklook(
    base_bundle,
    ext='MOMENT0',
    cmap='inferno',
    contours=[
        {
            'source': cont_bundle,
            'ext': 'MOMENT0',
            'levels': 'auto',
            'colors': 'cyan',
        }
    ],
)
```

contour 側も 2D passthrough にできます。

---

# 33. contour level 設計の現行版まとめと詳細例

この章は、v5p5 の contour level 周りを、**なるべく誤解なく**使うための整理章です。

現行版の基本方針:

- 旧版互換を重視して `auto` は fraction legacy とする
- `fraction` は相対レベル
- `manual` は絶対値レベル
- `manual` には `level_scale` が非常に有用
- `rms` は補助的な高度機能として残す

---

## 33.1 `auto` は旧版互換

```python
pl.plot_map(
    base,
    contours=[{'source': ov, 'levels': 'auto', 'colors': 'white'}],
)
```

現行版では、`levels='auto'` は **旧版互換**として、
概ね正のピーク値の `10, 30, 50, 70, 90%` を使います。

これは、以前の利用コードを壊さないための方針です。

---

## 33.2 `level_mode='auto'` も同じ意味で使える

```python
pl.plot_map(
    base,
    contours=[{'source': ov, 'level_mode': 'auto', 'colors': 'white'}],
)
```

ただし、**旧版コードとの見た目互換を意識するなら `levels='auto'` の方が分かりやすい**です。

---

## 33.3 `fraction` はピークに対する割合

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'fraction',
        'fraction_levels': [0.1, 0.2, 0.4, 0.8],
        'colors': 'white',
    }],
)
```

これは「ピークの 10%, 20%, 40%, 80%」です。

---

## 33.4 `fraction` と `auto` の違い

```python
# auto
pl.plot_map(base, contours=[{'source': ov, 'levels': 'auto', 'colors': 'white'}])

# fraction を明示
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'fraction',
        'fraction_levels': [0.1, 0.3, 0.5, 0.7, 0.9],
        'colors': 'white',
    }],
)
```

この 2 つは、現行版ではほぼ同じ結果になります。
違いは、後者の方が**意図が明示的**である点です。

---

## 33.5 `manual` は絶対レベル値をそのまま使う

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'manual',
        'levels': [5, 10, 20, 40],
        'colors': 'cyan',
    }],
)
```

---

## 33.6 `manual + level_scale` は非常に有用

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'manual',
        'levels': [1, 2, 3, 4],
        'level_scale': 3.0,
        'colors': 'cyan',
    }],
)
```

この場合、実際の contour level は

- 3
- 6
- 9
- 12

になります。

この使い方は、
**「形」は `[1,2,3,4]` のまま保ち、scale だけ変えて複数試す**
場面で非常に便利です。

---

## 33.7 `manual + level_scale` をループで試す

```python
for s in [1.0, 2.0, 3.0, 5.0]:
    pl.plot_map(
        base,
        contours=[{
            'source': ov,
            'level_mode': 'manual',
            'levels': [1, 2, 3, 4],
            'level_scale': s,
            'colors': 'white',
        }],
        title=f'manual levels, scale={s}',
    )
```

---

## 33.8 `rms` は明示的に使うのが基本

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'rms',
        'rms': 0.35,
        'sigma_levels': [3, 5, 7],
        'colors': 'white',
    }],
)
```

この場合、実際の contour level は

- `3 * 0.35`
- `5 * 0.35`
- `7 * 0.35`

です。

---

## 33.9 `rms + sigma_scale` の意味

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'rms',
        'rms': 0.35,
        'sigma_levels': [1, 2, 3, 4],
        'sigma_scale': 3.0,
        'colors': 'white',
    }],
)
```

このとき、実際の contour level は

- `0.35 * 3 * 1`
- `0.35 * 3 * 2`
- `0.35 * 3 * 3`
- `0.35 * 3 * 4`

つまり、`3σ, 6σ, 9σ, 12σ` 相当です。

---

## 33.10 `3*[1,2,3,4]` 的な使い方を安全に表現する

Python そのものの `3*[1,2,3,4]` は、
`[1,2,3,4]` を 3 回繰り返す意味になってしまうので危険です。

現行版でおすすめなのは次の 2 通りです。

### manual の場合

```python
{
    'level_mode': 'manual',
    'levels': [1, 2, 3, 4],
    'level_scale': 3.0,
}
```

### rms の場合

```python
{
    'level_mode': 'rms',
    'rms': 0.35,
    'sigma_levels': [1, 2, 3, 4],
    'sigma_scale': 3.0,
}
```

---

## 33.11 `sigma` は `rms` の後方互換 alias

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'sigma',
        'rms': 0.35,
        'sigma_levels': [3, 5, 7],
        'colors': 'white',
    }],
)
```

現行版では、`sigma` は **後方互換 alias** として残っています。
ただし、新規コードでは `level_mode='rms'` の方を推奨します。

---

## 33.12 RMS 自動推定は補助機能と考える

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'rms',
        'sigma_levels': [3, 5, 7],
        'colors': 'white',
    }],
)
```

この書き方は可能ですが、RMS が明示されていないため、
内部では 2D データからの簡易推定に頼ります。

実務上の考え方:

- **まずは fraction / manual を優先**
- RMS は値を明示できるときに使う
- 自動推定は補助機能とみなす

---

## 33.13 3D 由来の moment0 では RMS contour の解釈に注意

理想的には、3D cube の baseline RMS と積分チャネル数から、
各画素ごとの moment0 RMS を推定するのが物理的に自然です。

ただし現行版の contour API は、まだ **一般的な scalar RMS** を中心に設計されています。
したがって、3D 由来の moment0 で厳密な S/N contour を描きたい場合には、
将来的な別実装や、あらかじめ作成した RMS map / S/N map を使う設計を検討してください。

---

## 33.14 `levels=[...]` を `fraction` や `rms` でも渡せてしまうことの注意

後方互換性のため、内部では `levels=[...]` が存在すると、それを最終値としてそのまま使う経路があります。

したがって、たとえば次のような書き方は可能ですが、**意図が分かりにくい**です。

```python
pl.plot_map(
    base,
    contours=[{
        'source': ov,
        'level_mode': 'fraction',
        'levels': [5, 10, 20],
        'colors': 'white',
    }],
)
```

新規コードでは、次の使い分けを推奨します。

- `fraction` → `fraction_levels=[...]`
- `manual` → `levels=[...]`
- `rms` → `sigma_levels=[...]` と `rms=...`

---

## 33.15 `quicklook()` で `manual + level_scale`

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    contours=[{
        'source': '13co_cube.fits',
        'mode': 'moment0',
        'vel_range': (20, 40),
        'level_mode': 'manual',
        'levels': [1, 2, 3, 4],
        'level_scale': 2.5,
        'colors': 'white',
    }],
)
```

---

## 33.16 `quicklook()` で `rms + sigma_scale`

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    contours=[{
        'source': '13co_cube.fits',
        'mode': 'moment0',
        'vel_range': (20, 40),
        'level_mode': 'rms',
        'rms': 0.25,
        'sigma_levels': [1, 2, 3, 4],
        'sigma_scale': 3.0,
        'colors': 'cyan',
    }],
)
```

---

## 33.17 `fraction` と `manual` を重ねる比較図

```python
pl.plot_map(
    base,
    contours=[
        {
            'source': ov,
            'level_mode': 'fraction',
            'fraction_levels': [0.2, 0.4, 0.6, 0.8],
            'colors': 'yellow',
            'label': 'fraction',
        },
        {
            'source': ov,
            'level_mode': 'manual',
            'levels': [5, 10, 20, 40],
            'colors': 'cyan',
            'label': 'manual',
        },
    ],
    legend=True,
)
```

---

## 33.18 どう使い分けるべきか

実務上のおすすめは次です。

- **とりあえず図を見たい** → `levels='auto'`
- **ピーク比で比較したい** → `fraction`
- **物理単位で contour を制御したい** → `manual`
- **ノイズ基準で contour を切りたい** → `rms`

特に、多くの図を素早く作る段階では
`auto` と `manual` の 2 系統が最も実用的です。

---

# 34. `plot_map()` から `quicklook()` / `plot_scene()` への移行例

この章は、過去コードを持っている人向けの移行メモです。

現行方針:

- `plot_map()` は **互換 API** として残す
- ただし、新規コードでは **`quicklook()` / `plot_scene()` を優先**する
- 単純な 1 枚 → `quicklook()`
- 複数レイヤー → `plot_scene()`

---

## 34.1 旧来の `make_2d_map() + plot_map()`

```python
m = pl.make_2d_map('cube.fits', mode='moment0', vel_range=(20, 40))
pl.plot_map(m, cmap='inferno', title='old style')
```

これは今でも有効です。

---

## 34.2 同じことを `quicklook()` で書く

```python
pl.quicklook(
    'cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    cmap='inferno',
    title='new style',
)
```

---

## 34.3 旧来の contour 付き `plot_map()`

```python
m = pl.make_2d_map('12co_cube.fits', mode='moment0', vel_range=(20, 40))
pl.plot_map(
    m,
    contours=[{
        'source': '13co_cube.fits',
        'mode': 'moment0',
        'vel_range': (20, 40),
        'levels': 'auto',
        'colors': 'white',
    }],
)
```

---

## 34.4 同じことを `quicklook(contours=...)` で書く

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    contours=[{
        'source': '13co_cube.fits',
        'mode': 'moment0',
        'vel_range': (20, 40),
        'levels': 'auto',
        'colors': 'white',
    }],
)
```

---

## 34.5 marker / catalog が要るなら `plot_scene()` へ進む

```python
pl.plot_scene(
    base={
        'kind': 'image',
        'source': '12co_cube.fits',
        'mode': 'moment0',
        'vel_range': (20, 40),
        'cmap': 'inferno',
    },
    overlays=[
        {
            'kind': 'contour',
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'white',
        },
        {'kind': 'catalog', 'table': df, 'label_col': 'name', 'color': 'cyan'},
        {'kind': 'beam', 'beam': 'auto'},
        {'kind': 'scalebar', 'config': {'length': 30, 'units': 'arcsec'}},
    ],
)
```

---

## 34.6 `plot_map()` を今後も使ってよい場面

次のような場面では、今後も `plot_map()` を使って問題ありません。

- 既存コードの保守
- すでに `Map2D` を明示的に作っているパイプライン
- 1 枚の 2D map を丁寧に制御したい
- `plot_scene()` ほどの複雑さは不要

一方、新規に書き始めるなら、まずは `quicklook()` / `plot_scene()` を選ぶ方が自然です。

---

# 35. 失敗しやすい例と対策

この章は、実際によく起こる混乱をまとめた実務的なメモです。

---

## 35.1 `vel_range` と `cmin/cmax` を混同する

誤り:

```python
pl.quicklook('cube.fits', vel_range=(0, 100), cmin=(20, 40))
```

`vel_range` は **積分する速度範囲**、
`cmin/cmax` は **表示色の範囲**です。

正しくは:

```python
pl.quicklook('cube.fits', mode='moment0', vel_range=(20, 40), cmin=0, cmax=100)
```

---

## 35.2 `plot_scene()` に行くべき場面で無理に `quicklook()` を使う

`quicklook()` は便利ですが、

- catalog
- marker
- text
- scalebar
- north_arrow
- 複数 image overlay

が要ると、`plot_scene()` の方が構成が自然です。

---

## 35.3 `fraction` なのに `levels=[...]` を使ってしまう

可能ではありますが、意図が伝わりにくいです。

推奨:

```python
{
    'level_mode': 'fraction',
    'fraction_levels': [0.2, 0.4, 0.6],
}
```

---

## 35.4 `manual` なのに `fraction_levels` を使ってしまう

これは意味の取り違えです。

- `manual` → `levels=[...]`
- `fraction` → `fraction_levels=[...]`
- `rms` → `sigma_levels=[...]`

と覚えるのが安全です。

---

## 35.5 `rms` に自動推定を期待しすぎる

2D map しかない状況では、RMS 自動推定は近似です。

- とりあえずの図 → `auto` / `fraction`
- 物理単位で見せたい → `manual`
- 信頼できるノイズ値がある → `rms`

が現実的です。

---

## 35.6 RGB の 3 channel が形状不一致

RGB は基本的に、3 枚が同じ shape であることを期待しています。
本格的な再投影は、将来拡張の対象です。

---

## 35.7 optical + radio では contour を第一候補にする

光学画像の上に電波を載せるとき、
まずは `contour` が最も安全です。

pseudo-color overlay は便利ですが、
WCS・beam・dynamic range の違いを意識する必要があります。

---

## 35.8 beam matching と display interpolation を混同しない

- `target_hpbw_arcsec` は **物理的に分解能をそろえる**操作
- `interpolation='bilinear'` は **表示上の補間**

この 2 つは全く別です。

---

## 35.9 crop の world / pixel を混同しない

```python
crop={'center': (150, 160), 'size': (100, 100), 'mode': 'pixel'}
```

と

```python
crop={'center': (83.82, -5.39), 'size': (0.2, 0.2), 'mode': 'world'}
```

は意味が全く違います。

---

## 35.10 `plot_map()` は古いから消すべき、ではない

現行方針は、

- `plot_map()` は互換 API として残す
- ただし、新規コードでは主役にしない

です。

つまり、既存コードの保守対象としては重要です。

---

## 35.11 `LINEFREE_USED` を per-spectrum 3D mask だと思い込まない

現行仕様では、`LINEFREE_USED` と `SIGNAL_MASK_USED` は **cube 全体で共通の 1D channel mask** です。

したがって、

- 「各 pixel のスペクトルごとに別の line-free が入っている」
- 「`LINEFREE_USED` をそのまま 3D voxel mask として contour できる」

とは考えないでください。

viewer が欲しいときは、

```python
viewer = pl.make_baseline_viewer_bundle(bundle, mode='signal')
```

のように、1D channel mask を使って 3D cube を作り直してから扱います。

---

## 35.12 `plot_map(bundle)` と書いて 3D bundle を 2D 扱いしない

安全な書き方は次の 2 つです。

```python
pl.plot_map(bundle, ext='MOMENT0')
```

または

```python
m = pl.make_2d_map(bundle, mode='moment0', vel_range=(20, 40))
pl.plot_map(m)
```

---

## 35.13 bundle の 1D ext を 2D 表示 API に直接渡さない

```python
# 避ける
# pl.quicklook(bundle, ext='LINEFREE_USED')
```

`LINEFREE_USED` などは 1D spectral mask なので、2D image として扱う API には向きません。

---

## 35.14 baseline viewer bundle は別の特別クラスではない

viewer bundle は通常の `OTFBundle` として plotting に渡します。

- 2D ext を持つとは限らない
- `data` は 3D cube
- したがって `mode='moment0'` などで 2D 化して使う

と覚えておくと混乱しにくいです。

---

# 36. 追加の publication / presentation レシピ

この章は、論文図・スライド図・比較図を作る実務例をさらに増やした章です。

---

## 36.1 1 枚の論文図を `quicklook()` だけで作る

```python
pl.quicklook(
    '12co_cube.fits',
    mode='moment0',
    vel_range=(20, 40),
    cmap='inferno',
    contours=[
        {
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'white',
            'linewidths': 0.9,
        }
    ],
    beam='auto',
    title='Integrated intensity',
    save='paper_like_quicklook.pdf',
    show=False,
)
```

---

## 36.2 `plot_scene()` で論文図を組む

```python
pl.plot_scene(
    base={
        'kind': 'image',
        'source': '12co_cube.fits',
        'mode': 'moment0',
        'vel_range': (20, 40),
        'cmap': 'inferno',
        'beam': 'auto',
        'colorbar_label': '12CO [K km/s]',
    },
    overlays=[
        {
            'kind': 'contour',
            'source': '13co_cube.fits',
            'mode': 'moment0',
            'vel_range': (20, 40),
            'levels': 'auto',
            'colors': 'white',
            'label': '13CO',
        },
        {'kind': 'scalebar', 'config': {'length': 30, 'units': 'arcsec'}},
        {'kind': 'north_arrow'},
    ],
    legend=True,
    title='Publication-ready example',
    save='paper_scene.pdf',
    show=False,
)
```

---

## 36.3 スライド用に high-contrast にする

```python
pl.quicklook(
    'moment0.fits',
    cmap='turbo',
    norm_mode='asinh',
    norm_percentile=(1, 99.8),
    contours=[{'levels': 'auto', 'colors': 'black', 'linewidths': 1.2}],
    title='Slide figure',
)
```

---

## 36.4 3 枚の図を同じ contour 定義で量産する

```python
contour_cfg = {
    'source': '13co_cube.fits',
    'mode': 'moment0',
    'vel_range': (20, 40),
    'level_mode': 'manual',
    'levels': [1, 2, 3, 4],
    'level_scale': 3.0,
    'colors': 'white',
}

for src in ['12co_cube.fits', '12co_cube_alt1.fits', '12co_cube_alt2.fits']:
    pl.quicklook(
        src,
        mode='moment0',
        vel_range=(20, 40),
        contours=[dict(contour_cfg)],
        save=f'{src}_compare.pdf',
        show=False,
    )
```

---

## 36.5 colorbar label を明示して比較図を揃える

```python
for src, out in [
    ('12co_cube.fits', '12co.pdf'),
    ('13co_cube.fits', '13co.pdf'),
    ('c18o_cube.fits', 'c18o.pdf'),
]:
    pl.quicklook(
        src,
        mode='moment0',
        vel_range=(20, 40),
        colorbar_label='Integrated intensity [K km/s]',
        save=out,
        show=False,
    )
```

---

# 37. CLI の追加例（現行 contour 仕様中心）

この章は、既存の CLI cookbook を補うための追加例です。

---

## 37.1 `--contour-levels` を与えずに最短 contour

```bash
python plotting_updated_v5p5.py cube.fits \
  --mode moment0 --vel-range 20,40 \
  --contour-color white
```

現行版では、contour 系オプションが 1 つでもあれば contour を作れます。
`--contour-levels` は必須ではありません。

---

## 37.2 旧版互換 auto contour を CLI で明示

```bash
python plotting_updated_v5p5.py cube.fits \
  --mode moment0 --vel-range 20,40 \
  --contour-level-mode auto \
  --contour-color white
```

---

## 37.3 fraction contour をより細かく指定

```bash
python plotting_updated_v5p5.py cube.fits \
  --mode moment0 --vel-range 20,40 \
  --contour-level-mode fraction \
  --contour-fractions 0.05,0.10,0.20,0.40,0.80 \
  --contour-color cyan
```

---

## 37.4 manual contour を絶対値で指定

```bash
python plotting_updated_v5p5.py cube.fits \
  --mode moment0 --vel-range 20,40 \
  --contour-level-mode manual \
  --contour-levels 5,10,20,40 \
  --contour-color white
```

---

## 37.5 manual contour + level_scale

```bash
python plotting_updated_v5p5.py cube.fits \
  --mode moment0 --vel-range 20,40 \
  --contour-level-mode manual \
  --contour-levels 1,2,3,4 \
  --contour-level-scale 3 \
  --contour-color white
```

---

## 37.6 RMS contour + sigma_scale

```bash
python plotting_updated_v5p5.py cube.fits \
  --mode moment0 --vel-range 20,40 \
  --contour-level-mode rms \
  --contour-rms 0.25 \
  --contour-sigmas 1,2,3,4 \
  --contour-sigma-scale 3 \
  --contour-color cyan
```

---

## 37.7 contour 線幅・アルファ・負側 contour を同時に指定

```bash
python plotting_updated_v5p5.py cube.fits \
  --mode moment0 --vel-range 20,40 \
  --contour-level-mode rms \
  --contour-rms 0.25 \
  --contour-sigmas 3,5,7 \
  --contour-negative \
  --contour-linewidth 1.3 \
  --contour-alpha 0.7 \
  --contour-color white
```

---

## 37.8 publication 用の CLI 例

```bash
python plotting_updated_v5p5.py 12co_cube.fits \
  --mode moment0 --vel-range 20,40 \
  --cmap inferno --norm-mode asinh --norm-percentile 1,99.7 \
  --beam --scalebar --north-arrow \
  --contour-level-mode auto --contour-color white \
  --title 'Integrated intensity' \
  --output figure.pdf
```

---

## 37.9 OTF bundle を使った比較図を量産する

```python
base_bundle = pl.read_otf_bundle('12co_bundle.fits')
cont_bundle = pl.read_otf_bundle('13co_bundle.fits')

for vr in [(10, 20), (20, 30), (30, 40)]:
    pl.quicklook(
        base_bundle,
        mode='moment0',
        vel_range=vr,
        cmap='inferno',
        contours=[
            {
                'source': cont_bundle,
                'mode': 'moment0',
                'vel_range': vr,
                'levels': 'auto',
                'colors': 'white',
            }
        ],
        save=f'compare_{vr[0]}_{vr[1]}.pdf',
        show=False,
    )
```

Python 上で bundle を保持したまま図を量産したいときに便利です。

---

# 38. この cookbook の現行版における使い分けの結論

最後に、使い分けを短くまとめます。

## 38.1 初心者

- まず `quicklook()`
- contour が欲しければ `quicklook(contours=...)`
- 何が起きているか知りたければ `describe=True`

## 38.2 中級者

- 2D map を明示的に作るなら `make_2d_map()` + `plot_map()`
- 旧コードを保守するときも `plot_map()` でよい
- ただし新規コードは `quicklook()` 寄りが自然

## 38.3 上級者

- 複数 line, optical, catalog, annotation を統合するなら `plot_scene()`
- reprojection や frame 混在を意識的に扱う

## 38.4 contour level

- まず `auto`
- ピーク比を明示したいなら `fraction`
- 絶対値で管理したいなら `manual`
- ノイズ基準なら `rms`

## 38.5 今後も残る方針

- `plot_map()` / `plot_rgb()` は互換 API として残す
- ただし中心 API は今後 `quicklook()` / `plot_scene()` に寄る
- cookbook もこの方針で今後さらに増補していく

## 38.6 OTF bundle / baseline viewer bundle

- OTF products を Python 上でそのまま扱うなら `OTFBundle` を直接渡す
- bundle の 2D ext をそのまま見たいなら `ext='MOMENT0'` などを使う
- bundle の 3D cube から 2D を作りたいなら `mode='moment0'` などを使う
- baseline viewer bundle も通常の `OTFBundle` として扱い、必要に応じて 2D 化する

## 38.7 1D spectral mask の理解

- 現行の `LINEFREE`, `LINEFREE_USED`, `SIGNAL_MASK_USED` は global 1D channel mask
- 2D image や per-spectrum 3D mask と混同しない
- viewer が欲しければ `make_baseline_viewer_bundle()` で 3D source に変換してから描く

## 38.8 path / HDU / HDUList / bundle の選び方

- 単にファイルを読んで表示するだけなら path / `HDUList` で十分
- Python 上で途中生成した products をそのまま渡したいなら `OTFBundle` が便利
- cookbook ではまず `quicklook()` を第一候補とし、複雑な overlay は `plot_scene()` に進む
