
# profile_view 詳細説明書（日本語）

## 0. この説明書の位置づけ

本書は、`sd_radio_spectral_fits.profile_view`（旧 `plotting` 系）に含まれる可視化機能の詳細説明書です。  
単なる使用例集ではなく、

- どのモジュールが何を担当するか
- どの関数・クラスが公開 API か
- 各パラメータが何を意味するか
- 速度軸・周波数軸・基線窓 `BSL_WINF`・`rest_freq` 上書きがどう相互作用するか
- `grid` の gridding がどのように座標をセルへ割り当てるか
- PDF 出力・固定紙面サイズ・余白指定がどう効くか

を、実装準拠で整理したものです。

本パッケージは、**Scantable を 1 本ずつ見る viewer、ページ単位に並べる montage、座標格子へ並べる grid** の 3 系統を持ちます。  
いずれも現在の実装は、**row ごとの軸情報（`CRVAL1/CDELT1/CRPIX1/CTYPE1/CUNIT1/RESTFREQ/RESTFRQ/SPECSYS/SSYSOBS/VELOSYS/VFRAME` など）を尊重して、表示のたびに物理軸を組み立てる**方針です。

---

## 1. モジュール構成

`profile_view` 直下の主要ファイルは以下です。

- `__init__.py`
  - 上位 API の窓口
  - `plot_profile_map()` を提供
  - `view_spectra`, `SpectralViewer`, `ProfileMapMontageViewer`, `ProfileMapGridViewer` を re-export
- `viewer.py`
  - 1 本のスペクトルを順送りで見る対話 viewer
- `montage.py`
  - 複数スペクトルをページ単位で並べる montage viewer
- `grid.py`
  - 行ごとの天球座標をもとに格子へ割り当てて profile map を作る grid viewer
- `windowmask.py`
  - `BSL_WINF` と `rest_freq` 上書きの整合を取るための窓変換・RMS 計算ユーティリティ
- `utils.py`
  - 軸変換、補間、平滑化、基線再構成、PDF 出力、紙面レイアウト補助

---

## 2. まず押さえるべき共通概念

## 2.1 入力の受け方

多くの API は次のいずれかを受け取れます。

- `Scantable` オブジェクト
- SDFITS などのファイルパス文字列
- それらのリスト

複数入力が与えられると、内部で `_merge_scantables()` が呼ばれ、**複数の Scantable を 1 つへ縦結合**します。  
このとき data は **VLA 互換の list-of-arrays** として統合されます。したがって、各 row のチャンネル数が異なっていても扱えます。

## 2.2 行選択

`rows` / `exclude_rows` は viewer, montage, grid で共通です。

- `rows`
  - 表示対象 row を指定
- `exclude_rows`
  - 表示対象から除外する row を指定
- 両者は同時指定不可

指定形式は、内部の `_parse_row_selector()` に委ねられていますが、典型的には以下が使えます。

- `0`
- `"0,2,5"`
- `"10:20"`
- `[0, 1, 5]`

内部では

- 重複除去
- 昇順化
- 範囲外チェック
- `ORIG_ROW_ID` の保存

が行われます。

## 2.3 x 軸の 3 モード

`xaxis` は共通して次の 3 種です。

- `"vel"` : 速度軸
- `"freq"` : 周波数軸
- `"chan"` : チャンネル番号

`axis_type` は旧 API 互換の別名で、`xaxis` が未指定または `"vel"` のときだけ `xaxis` の代わりに使われます。

## 2.4 速度軸の考え方

現在実装の要点は次の通りです。

1. row ごとの WCS と `SPECSYS/SSYSOBS`、`VELOSYS/VFRAME` を見て、まず **native な速度軸**を作る
2. `rest_freq` がユーザー指定され、row 側の `RESTFREQ/RESTFRQ` と異なる場合、`recalculate_velocity_axis()` で **display 用の速度軸**へ再計算する
3. `BSL_WINF` は database / row では native velocity 窓とみなし、表示時には **display velocity へ移してから x 軸へ変換**する

つまり、本パッケージは

- native velocity
- display velocity
- current xaxis

の 3 段階を意識して動いています。

## 2.5 `rest_freq` 上書きの意味

`rest_freq` は「WCS の `RESTFREQ/RESTFRQ` を書き換える」ものではなく、**表示上の速度解釈を別の静止周波数に置き換える**ために使われます。  
したがって、

- 速度軸
- fit window の位置
- `freq` 表示時の周波数軸
- RMS(win) の評価範囲

が連動して変化します。

## 2.6 `BSL_WINF` と `BSL_COEF`

本パッケージでは baseline 関連情報として主に次を参照します。

- `BSL_WINF`
  - 基線フィットに使った速度窓
- `BSL_COEF`
  - 基線モデル
  - ベクトルそのもの、または多項式係数のいずれか
- `BSL_SCALE`
  - baseline モデルがどの温度スケールで定義されているか
- `BSL_RMS`
  - baseline fit の RMS

`BSL_COEF` は実装上、

- 長さ = チャンネル数 → baseline ベクトルそのもの
- それ以外 → `np.polyval()` 用の多項式係数

として解釈されます。

## 2.7 `TEMPSCAL` / `BEAMEFF`

表示スケール切替 `t` は概ね

- `TA*`
- `TR*`

の切替です。  
row ごとに `TEMPSCAL` と `BEAMEFF` を見て、表示直前に変換します。

重要な点:

- viewer では `TR*` へ切り替える前に `BEAMEFF` の有無をチェックする
- montage / grid でも row ごとの `BEAMEFF` / `TEMPSCAL` を見て変換する
- baseline も `BSL_SCALE` があればそのスケールを尊重してから変換する

## 2.8 VLA（可変長スペクトル）

このパッケージは **可変長 `DATA`** を前提として動けるように書かれています。

- `viewer` は row ごとにその長さで軸を生成
- `montage` も row ごとに軸を作る
- `grid` はセル coadd 時、参照 row の velocity 軸へ補間して平均する

したがって、各 row のチャンネル数や WCS が完全一致しない場合でも、少なくとも表示は可能です。  
ただし `grid` の coadd は **参照 row の速度軸へ線形補間して `np.nanmean()` する**実装なので、厳密な再サンプリング戦略を自分で選びたい場合は別途 upstream で揃えておく方が安全です。

---

## 3. 公開 API の全体像

実際に外から使う主要 API は次の 5 つです。

- `plot_profile_map(...)`
- `view_spectra(...)`
- `SpectralViewer(...)`
- `ProfileMapMontageViewer(...)`
- `ProfileMapGridViewer(...)`

以下で順に説明します。

---

## 4. `plot_profile_map`

### シグネチャ

```python
plot_profile_map(input_data: Union[str, Scantable], **kwargs)
```

### 役割

profile map 表示の総合入口です。  
`mode` を見て、

- `mode="montage"` または `"mosaic"` → `ProfileMapMontageViewer`
- `mode="grid"` または `"gridding"` → `ProfileMapGridViewer`

を起動します。

### 主要引数

- `input_data`
  - `Scantable` またはファイルパス
- `mode`
  - `"montage"` が既定
  - `"grid"` を指定すると grid viewer
- `axis_type`
  - `xaxis` の別名
- `xranges`
  - `xrange` の別名
- `yranges`
  - `yrange` の別名

### 注意

この関数は薄いラッパーです。  
実際のパラメータの意味は `ProfileMapMontageViewer` または `ProfileMapGridViewer` の説明に従います。

### 使用例

```python
import sd_radio_spectral_fits.profile_view as pv

pv.plot_profile_map(sc, mode="montage", xaxis="vel", nrows=4, ncols=4)
pv.plot_profile_map(sc, mode="grid", coord="radec", ref_point=(83.82208, -5.39111), x_grid=30, y_grid=30)
```

---

## 5. `view_spectra` と `SpectralViewer`

## 5.1 `view_spectra`

### シグネチャ

```python
view_spectra(input_data: Union[Scantable, str, Sequence[Union[Scantable, str]]], **kwargs)
```

### 役割

`SpectralViewer(...)` の薄いラッパーです。  
戻り値は `SpectralViewer` インスタンスです。

---

## 5.2 `SpectralViewer`

### シグネチャ

```python
SpectralViewer(
    scantable,
    rows=None,
    exclude_rows=None,
    xaxis="vel",
    axis_type=None,
    rest_freq=None,
    xrange=None,
    yrange=None,
    autoscale_y=True,
    show_fitwin_rms=True,
    show_fit_windows=True,
    smooth_mode="none",
    smooth_width=1,
    box_downsample=False,
    box_policy="trim",
    rms_on="raw",
    show_top_axis=True,
    save_dir=".",
    save_prefix="spectrum",
    save_pdf=None,
    max_pdf_pages=100,
    show=True,
    figsize=None,
    content_aspect=None,
    paper_margins=None,
)
```

### 役割

1 本ずつスペクトルを表示し、キーボードで前後移動しながら確認する viewer です。  
個々の row を精査する用途に最も向いています。

### パラメータ詳細

#### 入力と row 選択

- `scantable`
  - `Scantable`、ファイルパス、またはそれらのリスト
- `rows`
  - 表示対象 row のみ残す
- `exclude_rows`
  - 表示対象から除外
  - `rows` と同時指定不可

#### 軸関係

- `xaxis`
  - `"vel"`, `"freq"`, `"chan"`
- `axis_type`
  - 旧名。`xaxis` の別名
- `rest_freq`
  - 表示用静止周波数 [Hz]
  - row 側の `RESTFREQ/RESTFRQ` と異なる場合は、速度軸と fit window を display 側へ再計算
- `xrange`
  - 表示 x 範囲
  - 入力は「現在の初期 xaxis の単位」で解釈され、内部で速度範囲へ変換して保持
  - その後 `xaxis` を切り替えても、同じ物理速度範囲が保たれる
- `yrange`
  - y 範囲固定

#### y 軸と RMS

- `autoscale_y`
  - `yrange` 未指定時に表示範囲から自動スケール
- `show_fitwin_rms`
  - タイトル行に `RMS(win)` を表示するか
- `show_fit_windows`
  - `BSL_WINF` を緑の帯で描くか
- `rms_on`
  - `"raw"` または `"display"`
  - RMS 計算を raw グリッドで行うか、平滑化・ダウンサンプル後の表示グリッドで行うか

#### 平滑化

- `smooth_mode`
  - `"none"` または `"boxcar"` が実用上の主対象
  - キー操作では `"none" <-> "boxcar"` を循環
- `smooth_width`
  - 平滑幅
- `box_downsample`
  - `boxcar` のとき、平滑のみではなく箱平均ダウンサンプルも行うか
- `box_policy`
  - `box_downsample=True` のときの端数処理
  - `"trim"`: 端数を切る
  - `"pad"`: NaN で埋めてから平均

#### 表示要素

- `show_top_axis`
  - 上側補助軸を出すか
  - `vel` 表示なら上に channel
  - `freq` 表示なら上に velocity
  - `chan` 表示なら上に velocity または frequency 相当

#### 保存

- `save_dir`
  - 手動保存先ディレクトリ
- `save_prefix`
  - 手動保存時ファイル名接頭辞
- `save_pdf`
  - 指定すると、全 row を順に PDF へ書き出す
- `max_pdf_pages`
  - PDF 出力の上限ページ数

#### 図サイズ・紙面

- `show`
  - `True` なら対話表示
- `figsize`
  - 通常の `(width, height)` でもよい
  - `"A4L"` などの紙サイズ指定も `parse_figsize()` で受理
- `content_aspect`
  - プロット本体の width/height を固定したいときに使う
  - 例: `"4:3"`, `"1:1"`, `1.3333`
- `paper_margins`
  - 固定紙面モードの余白
  - `(left, right, bottom, top)` を mm または inch 指定文字列で与える

### viewer の内部動作

1. 入力を merge
2. row 選択
3. 現在 row の `TEMPSCAL` と `BEAMEFF` を見て表示スケールを決定
4. row ごとに `_row_vel_lsrk_and_oriented()` で物理軸を生成
5. 必要なら `rest_freq` 上書きで velocity/frequency を再計算
6. `BSL_COEF` があれば baseline を再構成
7. `show_fit_windows=True` なら `BSL_WINF` を current xaxis へ変換して描画
8. `compute_rms_win()` で RMS(win) を算出

### viewer のキーボード操作

- `n / p` : 次 / 前
- `b` : baseline view 切替（Original <-> Baseline Added）
- `x` : x 軸切替
- `t` : `TA* <-> TR*`
- `m` : smoothing 切替
- `[` / `]` : smoothing 幅
- `d` : downsample 切替
- `r / R` : `rms_on` 切替
- `T` : 上軸切替
- `s` : 現在スペクトルを PDF 保存
- `q` : 終了

### viewer の使用例

#### 最も基本

```python
from sd_radio_spectral_fits.profile_view import view_spectra
view_spectra(sc)
```

#### 速度範囲を指定して見る

```python
view_spectra(
    sc,
    xaxis="vel",
    xrange=(-20, 40),
    show_fit_windows=True,
    smooth_mode="boxcar",
    smooth_width=3,
)
```

#### `rest_freq` を上書きして表示解釈を変える

```python
view_spectra(
    sc,
    xaxis="vel",
    rest_freq=115.2712018e9,
    xrange=(-10, 20),
)
```

#### A4 横で PDF 化

```python
view_spectra(
    sc,
    show=False,
    save_pdf="viewer_all.pdf",
    figsize="A4L",
    paper_margins="15mm,10mm,15mm,20mm",
    content_aspect="4:3",
)
```

---

## 6. `ProfileMapMontageViewer`

### シグネチャ

```python
ProfileMapMontageViewer(
    scantable,
    rows=None,
    exclude_rows=None,
    xaxis="vel",
    axis_type=None,
    rest_freq=None,
    imin=0,
    imax=None,
    nrows=5,
    ncols=4,
    xrange=None,
    yrange=None,
    autoscale_y=True,
    annotate_rms=True,
    show_fit_windows=True,
    smooth_mode="none",
    smooth_width=1,
    box_downsample=False,
    box_policy="trim",
    rms_on="raw",
    save_dir=".",
    save_prefix="montage",
    save_pdf=None,
    max_pdf_pages=100,
    show=True,
    figsize=None,
    content_aspect=None,
    paper_margins=None,
)
```

### 役割

複数スペクトルを `nrows x ncols` のページに並べて、ページ送りで確認する viewer です。  
「行番号順に順番に眺める」用途に向いています。

### パラメータ詳細

#### 入力・row 範囲

- `scantable`, `rows`, `exclude_rows`
  - viewer と同じ
- `imin`, `imax`
  - 表示対象 row の開始・終了インデックス
  - Python 的な半開区間
  - `rows/exclude_rows` の後にさらに適用される「表示スライス」と考えるとわかりやすい

#### ページ構成

- `nrows`, `ncols`
  - 1 ページあたりの行数・列数
  - `page_size = nrows * ncols`

#### 軸・平滑化・RMS

- `xaxis`, `axis_type`, `rest_freq`, `xrange`, `yrange`
  - viewer と同様
- `autoscale_y`
  - ページ中の表示範囲に合わせて共有 y を調整
- `annotate_rms`
  - 各パネル左上に RMS を書くか
- `show_fit_windows`
  - `BSL_WINF` を表示するか
- `smooth_mode`, `smooth_width`, `box_downsample`, `box_policy`, `rms_on`
  - viewer と同様

#### 保存・レイアウト

- `save_dir`, `save_prefix`, `save_pdf`, `max_pdf_pages`, `show`, `figsize`, `content_aspect`, `paper_margins`
  - viewer とほぼ同じ

### montage の特徴

- `Standardizer` を使わず、**常に row ごとの独立処理ルート**を使用
- VLA をそのまま扱える
- 各パネルの軸は row ごとに再計算
- `axes_flat = axes_arr[::-1, :].ravel()` としているため、ページ上の配置は「見た目の左上から右下」ではなく、実装上のフラット順序に注意が必要

### montage のキーボード操作

- `n / p` : 次 / 前ページ
- `b` : baseline 表示切替
- `w` : fit windows 表示切替
- `x` : x 軸切替
- `t` : `TA* <-> TR*`
- `m` : smoothing 切替
- `[` / `]` : smoothing 幅
- `d` : downsample
- `r / R` : `rms_on`
- `s` : 現ページ PDF 保存
- `q` : 終了

### montage の使用例

#### 4x4 のページで確認

```python
from sd_radio_spectral_fits.profile_view import ProfileMapMontageViewer

ProfileMapMontageViewer(
    sc,
    nrows=4,
    ncols=4,
    xaxis="vel",
    xrange=(-20, 40),
)
```

#### row 100〜199 だけを A4 横で保存

```python
ProfileMapMontageViewer(
    sc,
    imin=100,
    imax=200,
    nrows=5,
    ncols=4,
    show=False,
    save_pdf="montage_subset.pdf",
    figsize="A4L",
    paper_margins="12mm,10mm,12mm,18mm",
)
```

#### `rest_freq` を切り替えて fit window も追随させる

```python
ProfileMapMontageViewer(
    sc,
    xaxis="vel",
    rest_freq=110.201354e9,
    show_fit_windows=True,
    annotate_rms=True,
)
```

---

## 7. `ProfileMapGridViewer`

### シグネチャ

```python
ProfileMapGridViewer(
    st,
    rows=None,
    exclude_rows=None,
    mode="all",
    combine="mean",
    coord="radec",
    projection="SFL",
    ref_point=None,
    x_grid=30.0,
    y_grid=30.0,
    corner_offsets=None,
    grid_bounds_offsets=None,
    grid_anchor_offsets=None,
    grid_tol=0.5,
    invert_x=True,
    show_rms=True,
    show_coord_text=True,
    nrows=None,
    ncols=None,
    xaxis="vel",
    axis_type=None,
    rest_freq=None,
    imin=0,
    imax=None,
    xrange=None,
    yrange=None,
    annotate_rms=True,
    smooth_mode="none",
    smooth_width=1,
    box_downsample=False,
    box_policy="trim",
    rms_on="raw",
    show=True,
    square_aspect=True,
    box_padding=0.0,
    show_baseline=False,
    show_id=True,
    show_fit_windows=True,
    show_grid="auto",
    offset_unit="arcsec",
    save_dir=".",
    save_prefix="grid",
    save_pdf=None,
    max_pdf_pages=100,
    figsize=None,
    content_aspect=None,
    paper_margins=None,
)
```

### 役割

天球座標に基づいて各 row を格子セルへ割り当て、そのセルに属するスペクトルを表示する viewer です。  
このモジュールでは **gridding の考え方**が最も重要です。

---

## 8. grid の gridding は何をしているか

grid の内部処理は、概ね次の順です。

1. `coord` に応じて row から
   - `RA/DEC`
   - または `GLON/GLAT`
   を取り出す
2. `ref_point=(lon0, lat0)` を原点とする
3. `projection` に応じて x オフセットを計算する
4. そのオフセットを `x_grid`, `y_grid` の格子間隔で量子化する
5. `grid_tol` を使って「格子点から離れすぎた row」を弾く
6. できた `(ix, iy)` に row を投入する
7. 表示上は `nrows`, `ncols` ごとにページ分割して描く

ここで重要なのは、**この grid は一般的な意味での full WCS reprojection ではない**ことです。  
ただし `SFL/GLS` 分岐は単なる近似ではなく、`x = (lon-lon0) cos(lat)`, `y = lat-lat0` という **Sanson-Flamsteed / Sinusoidal 型の投影式** をそのまま用いています。  
一方で `TAN/SIN/GNOMONIC` 分岐は、完全な投影式ではなく `cos(lat0)` を使う局所直交近似です。

---

## 9. grid の gridding パラメータ詳細（最重要）

## 9.1 `coord`

- `"radec"`
  - `RA`, `DEC` 系列を使う
  - 候補列は大小文字を無視して
    - `RA`, `RA_DEG`, `RA2000`, `RAJ2000`
    - `DEC`, `DEC_DEG`, `DEC2000`, `DECJ2000`
- `"gal"`
  - `GLON`, `GLAT` 系列を使う
  - 候補列
    - `GLON`, `L`, `GAL_L`, `LON`
    - `GLAT`, `B`, `GAL_B`, `LAT`
- `"auto"`
  - **現実装では無効**
  - 指定すると例外

### 実務上の意味

- 赤経赤緯ベースの通常の観測なら `coord="radec"`
- 銀河面や銀経銀緯で見たいなら `coord="gal"`

---

## 9.2 `projection`

内部の投影モード指定です。  
ただし **すべてが同じ精度・意味で「完全な投影」なのではない** 点に注意が必要です。

### 実装されている分岐

- `"TAN"`, `"SIN"`, `"GNOMONIC"`
  - `cos(lat0)` を使う
  - 基準点緯度 `lat0` で固定した局所直交近似
  - 名前は投影名ですが、実装は full TAN / full SIN / full gnomonic ではありません
- `"GLS"`, `"SFL"`, `"SINE"`
  - `cos(lat)` を各 row ごとに使う
  - `x = (lon-lon0) cos(lat)`, `y = lat-lat0` 型
  - **Sanson-Flamsteed / Sinusoidal 型の投影式そのもの**として扱われます
  - ただし実装上は、x 側では reference latitude を明示的には使わず、事実上 `lat_ref = 0` の形です
- それ以外
  - cos 因子なし（`1.0`）
  - 生の差分座標に近い扱い

### 数式イメージ

`dlon = wrap(lon - lon0)` とすると、x オフセットは

- TAN/SIN/GNOMONIC:
  - `x = ± dlon * cos(lat0) * 3600`
- GLS/SFL/SINE:
  - `x = ± dlon * cos(lat) * 3600`

y は

- `y = (lat - lat0) * 3600`

です。

このため `SFL` は「reference latitude まわりの局所近似」ではなく、**各 row の緯度をそのまま使う投影**です。  
一方で y については `lat0` を差し引いているため、表示上の原点は `ref_point` に置かれます。

### 実務上の選び方

- 小さめの RA/Dec マップで、縦横の直交性を重視する
  - `projection="TAN"` または `"SIN"`
- 単一鏡の広めのマップで、SFL/GLS 的な球面投影をそのまま使いたい
  - `projection="SFL"` または `"GLS"`

### 注意

- `projection="SFL"` / `"GLS"` は、ここでは単なる cos 補正ではなく、投影式として解釈してよいです
- ただし外部 FITS 画像の WCS 全般と厳密一致する保証まではなく、回転や高次の WCS 要素までは扱っていません
- `projection="SIN"` は名前の印象ほど「完全な正射影」ではありません

---

## 9.3 `ref_point`

原点 `(lon0, lat0)` です。  
指定しない場合は、選ばれた座標列の有限値に対する **中央値**が使われます。

### `coord="radec"` の特殊処理

RA 側については、データ列の絶対値が 24.5 以下なら hour angle 風とみなし、`_maybe_hours_to_deg()` により **hours -> degrees 変換**される可能性があります。  
同様に `ref_point[0]` も 24.5 以下で、データ側も hours らしい場合は `* 15` されます。

### 実務上の推奨

RA/Dec のときは、誤解を避けるため **degree で与える**のが最も安全です。

```python
ref_point=(83.82208, -5.39111)
```

---

## 9.4 `x_grid`, `y_grid`

格子間隔 [arcsec] です。  
セルサイズそのものと考えてよいです。

- `x_grid`
  - x 方向格子間隔
- `y_grid`
  - y 方向格子間隔

### 実務上の意味

- 観測点が 30 秒角刻みなら `x_grid=30`, `y_grid=30`
- 5-beam 相対配置込みで、最終的に 15 秒角格子へ寄せたいなら `15, 15`
- ただし、この grid viewer は「見た目の配置用」であり、高度な再グリッド補間器ではありません

---

## 9.5 `corner_offsets` / `grid_bounds_offsets`

表示グリッドの物理範囲を、原点からのオフセットで指定します。  
単位は [arcsec] です。

### 形

```python
((x1, y1), (x2, y2))
```

内部ではこれを

- `xmin_req, xmax_req`
- `ymin_req, ymax_req`

へソートして使います。

### `corner_offsets` と `grid_bounds_offsets` の関係

現実装では `grid_bounds_offsets` が新しい名前で、`corner_offsets` は後方互換 alias です。

- 両方指定して値が異なる → 例外
- 片方だけ指定 → それを使う
- 未指定 → データの min/max から自動決定

### 指定した場合と未指定の場合の違い

#### 指定した場合
- その矩形範囲に収まる格子を使う
- page 構成が安定しやすい
- 欠測セルも含めて「同じ地図サイズ」で比較しやすい

#### 未指定の場合
- データが存在する範囲から自動推定
- 最小限のグリッドになる
- 別データ間で見た目の枠が揺れやすい

---

## 9.6 `grid_anchor_offsets`

格子原点を、オフセット座標系のどこに置くかを明示するためのパラメータです。

### 形

```python
(x_anchor, y_anchor)   # arcsec
```

### 意味

これを指定すると、「格子線は常にこの anchor を通る」ようにセルが決まります。  
観測データの最小値から勝手に始めるのではなく、**0,0 基準の揃った格子**にしたいときに重要です。

### 指定時の挙動

データ点 `(x_arcsec, y_arcsec)` に対し、内部では概ね

```python
ix_abs = round((x_arcsec - xa) / dx)
iy_abs = round((y_arcsec - ya) / dy)
```

のようにセル番号を決めます。

### 典型例

原点対称のマップを必ず 0,0 基準の 30" 格子へ載せたい:

```python
grid_anchor_offsets=(0.0, 0.0)
```

### 指定すべき場面

- 日ごとに別観測を比べたい
- 同じ target を複数回見て、全く同じ格子に載せたい
- center からの相対オフセットを厳密に保ちたい

---

## 9.7 `grid_tol`

row が「どれだけ格子点に近ければ、そのセルに属するとみなすか」の許容係数です。

### 実装

```python
tolx = abs(grid_tol) * x_grid
toly = abs(grid_tol) * y_grid
```

が使われ、row の実オフセットと最近傍格子点との差が

- `abs(dx_resid) <= tolx`
- `abs(dy_resid) <= toly`

のときだけ採用されます。

### デフォルト `0.5` の意味

セル中心から半セル以内なら採用、という感覚です。

### 実務上の選び方

- 観測点がほぼ規則格子上にある
  - `0.3 ~ 0.5`
- 少しずれていても拾いたい
  - `0.6 ~ 0.8`
- 誤割当てを厳しく避けたい
  - `0.2 ~ 0.3`

---

## 9.8 `invert_x`

x 軸の符号を反転するかです。  
内部では `sgn = -1 if invert_x else +1` が掛かります。

### 実務上の意味

RA の増加方向と表示の左右向きの慣習に合わせるために使います。

- `True`
  - RA increasing を左向きにしたい電波/天文的慣習に近い
- `False`
  - 数学的な右増加へ寄せる

### 迷ったら

RA/Dec マップならまず `True` から試すのが自然です。

---

## 9.9 `nrows`, `ncols`

ここは **表示ページの行列数** です。  
物理グリッドの `nx, ny` そのものではありません。

- 未指定なら
  - `nrows = sol.ny`
  - `ncols = sol.nx`
  - つまり全グリッドを 1 ページに出そうとする
- 指定した場合
  - そのサイズでページに切り、複数ページになる

### 重要

`nrows/ncols` は **格子解そのもの**を変えません。  
変えるのはページ分割だけです。

---

## 9.10 `mode`

セル内に複数 row が入ったときの表示方法です。

- `"all"`
  - 各 row をスパゲッティ表示
- `"coadd"`
  - セル内の row を 1 本へ平均して表示

### 実装詳細

`"coadd"` のとき `_cell_spectrum()` が使われ、

1. 参照 row の速度軸を採る
2. 他 row を必要ならその軸へ線形補間
3. `np.nanmean()` で平均

します。

---

## 9.11 `combine`

現在の実装では `self.combine = ...` と保存されますが、**coadd の計算分岐には実質使われていません**。  
現状の `_cell_spectrum()` は `np.nanmean()` 固定です。

したがって説明としては、

- API 上は「将来的な結合方法の切替」を意図したパラメータ
- **現バージョンでは実効値は mean 相当のみ**

と理解するのが正確です。

---

## 9.12 `show_rms`

この引数は受理されますが、**現コードでは実質使われていません**。  
RMS 表示の実制御は主に `annotate_rms` 側で行われています。  
したがって、現時点では「互換的なダミー引数に近い」と見てよいです。

---

## 9.13 `show_coord_text`

各セルに座標文字列を書くかどうかです。  
RA/Dec または GLON/GLAT をセル左上に入れるときに効きます。

---

## 9.14 `offset_unit`

セル注釈などでオフセットを出すときの単位です。

- `"arcsec"` → 秒角
- `"arcmin"` → 分角
- `"deg"` → 度

内部では単なる表示スケール変換で、格子割り当て自体は常に arcsec 基準です。

---

## 9.15 `show_grid`

セル内プロットに点線グリッドを出すかです。

- `"auto"`
  - 現状では有効扱い
- `True`
  - 有効
- `False`
  - 無効
- `"none"`, `"false"`, `"0"`
  - 無効

---

## 9.16 `show_baseline`

初期状態で baseline-added view を有効にするかです。  
`BSL_COEF` がない row では baseline 表示にはなりません。

---

## 9.17 `show_id`

セル内に row ID 表示を出すかどうかです。  
複数 row が混ざるときの識別補助に使います。

---

## 9.18 `square_aspect` と `box_padding`

ページレイアウト用です。

- `square_aspect=True`
  - 各パネルを正方形に近づける
- `box_padding`
  - パネル間ギャップ係数
  - `0.1` なら「パネル幅の 10% 分の隙間」相当

---

## 10. grid の内部関数 `_solve_grid_manual_from_lonlat`

### シグネチャ

```python
_solve_grid_manual_from_lonlat(
    lon_deg,
    lat_deg,
    coord,
    ref_point,
    dx_arcsec,
    dy_arcsec,
    corner_offsets_arcsec=None,
    grid_anchor_offsets_arcsec=None,
    tol_factor=0.5,
    invert_x=True,
    projection="SFL",
) -> Optional[GridSolution]
```

### 役割

grid viewer の heart 部分です。  
与えられた天球座標列を、**整数格子インデックス `(ix, iy)`** へ落とします。

### 出力 `GridSolution`

- `coord`
  - `"radec"` or `"gal"`
- `projection`
  - 実際に使った projection 名
- `x_label`, `y_label`
  - 表示用ラベル名
- `ix`, `iy`
  - 各 row が属するセル番号
- `nx`, `ny`
  - 物理グリッド総数
- `dx`, `dy`
  - 格子間隔 [arcsec]
- `x0`, `y0`
  - ref_point
- `x_deg`, `y_deg`
  - 入力座標の保持
- `score`
  - 占有率と重複率から作った簡易スコア

### `score` の意味

内部ではおおよそ

- occupied cell fraction
- duplicate rate

から

```python
score = occ - 0.35 * dup_rate
```

で作られています。  
現在の `ProfileMapGridViewer` では複数候補比較には使っていませんが、grid 解の妥当性評価の名残です。

---

## 11. `ProfileMapGridViewer` のその他パラメータ

gridding 以外のパラメータは概ね montage に近いです。

- `xaxis`, `axis_type`, `rest_freq`
- `imin`, `imax`
- `xrange`, `yrange`
- `annotate_rms`
- `smooth_mode`, `smooth_width`
- `box_downsample`, `box_policy`
- `rms_on`
- `save_dir`, `save_prefix`, `save_pdf`, `max_pdf_pages`
- `figsize`, `content_aspect`, `paper_margins`

ただし grid では、coadd のときだけ**複数 row を 1 本へ束ねる**点が montage と異なります。

### grid のキーボード操作

- `n / p` : 次 / 前ページ
- `x` : x 軸切替
- `m` : smoothing 切替
- `[` / `]` : smoothing 幅
- `d` : downsample
- `t` : `TA* <-> TR*`
- `r / R` : `rms_on`
- `b` : baseline view 切替
- `w` : fit windows 切替
- `g` : グリッド線切替
- `a` : square aspect 切替
- `s` : PDF 保存
- `h` : help
- `q` : 終了

---

## 12. grid の実務的な使い方

## 12.1 最小例

```python
from sd_radio_spectral_fits.profile_view import ProfileMapGridViewer

ProfileMapGridViewer(
    sc,
    coord="radec",
    ref_point=(83.82208, -5.39111),
    x_grid=30.0,
    y_grid=30.0,
)
```

これは

- RA/Dec を使う
- `(83.82208, -5.39111)` を原点にする
- 30" x 30" 格子へ割り当てる
- 格子全体を 1 ページに出す

という意味です。

## 12.2 原点 0,0 に揃えた格子へ必ず載せたい

```python
ProfileMapGridViewer(
    sc,
    coord="radec",
    ref_point=(83.82208, -5.39111),
    x_grid=30.0,
    y_grid=30.0,
    grid_anchor_offsets=(0.0, 0.0),
)
```

これにより、別日データでも同じ 0,0 基準格子に載りやすくなります。

## 12.3 物理範囲を固定して比較したい

```python
ProfileMapGridViewer(
    sc,
    coord="radec",
    ref_point=(83.82208, -5.39111),
    x_grid=30.0,
    y_grid=30.0,
    grid_bounds_offsets=((-120.0, -90.0), (120.0, 90.0)),
)
```

これで、データが欠けていても常に同じ 240" x 180" の枠が使われます。

## 12.4 セル内は重ね描きではなく平均したい

```python
ProfileMapGridViewer(
    sc,
    mode="coadd",
    coord="radec",
    ref_point=(83.82208, -5.39111),
    x_grid=15.0,
    y_grid=15.0,
)
```

### 注意
この coadd は現在 `nanmean` 固定です。`combine="median"` のような指定は現実装では効きません。

## 12.5 銀河座標で見る

```python
ProfileMapGridViewer(
    sc,
    coord="gal",
    ref_point=(0.0, 0.0),
    projection="SFL",
    x_grid=60.0,
    y_grid=60.0,
)
```

## 12.6 ページ分割したい

```python
ProfileMapGridViewer(
    sc,
    coord="radec",
    ref_point=(83.82208, -5.39111),
    x_grid=30.0,
    y_grid=30.0,
    nrows=4,
    ncols=5,
)
```

ここで `nrows/ncols` は grid 解を変えるのではなく、**表示を 4x5 ずつに区切るだけ**です。

---

## 13. `windowmask.py` の詳細

`windowmask.py` は、`BSL_WINF` と `rest_freq` 上書きの整合を取るための共通 API です。

## 13.1 `fit_windows_disp_vel`

```python
fit_windows_disp_vel(
    win_str,
    rest_hz_meta,
    rest_hz_user,
    *,
    shift_threshold_hz=1.0,
) -> List[Tuple[float, float]]
```

### 役割
`BSL_WINF`（native velocity 窓）を **display velocity 窓**へ変換します。

### 引数
- `win_str`
  - 例: `"(-10, -5); (20, 40)"` 相当の window 文字列
- `rest_hz_meta`
  - row / header 側 rest frequency
- `rest_hz_user`
  - ユーザー指定 `rest_freq`
- `shift_threshold_hz`
  - これ未満の差なら restfreq 差は無視する

### 使いどころ
- 速度窓をまず display velocity へ移したいとき

## 13.2 `vel_windows_to_xaxis`

```python
vel_windows_to_xaxis(wins_vel, xaxis, vel_disp, freq_ghz, chan)
```

### 役割
display velocity の窓を current xaxis の窓へ変換します。

- `xaxis="vel"` → そのまま
- `xaxis="freq"` → 周波数へ
- `xaxis="chan"` → channel へ

## 13.3 `fit_windows_xaxis`

```python
fit_windows_xaxis(
    win_str,
    rest_hz_meta,
    rest_hz_user,
    xaxis,
    vel_disp,
    freq_ghz,
    chan,
    *,
    shift_threshold_hz=1.0,
)
```

### 役割
最もよく使う総合関数です。

`BSL_WINF`
→ native velocity
→ display velocity
→ current xaxis

を一発で処理します。

## 13.4 `vel_range_to_xbounds`

```python
vel_range_to_xbounds(xaxis, vrange, vel_disp, freq_ghz, chan)
```

### 役割
display velocity 範囲 `(vmin, vmax)` を current xaxis 範囲へ変換します。  
`xrange` 管理で使われます。

## 13.5 `mask_from_xbounds`, `mask_from_windows`

x から直接 mask を作る補助関数です。  
単純ですが、RMS 評価で重要です。

## 13.6 `prepare_rms_xy`

```python
prepare_rms_xy(
    x_raw, y_raw, rms_on,
    smooth_mode, smooth_width,
    box_downsample, box_policy
)
```

### 役割
RMS を raw グリッドで取るか、表示グリッドで取るかを統一的に処理します。

- `rms_on="raw"` → 元の x, y
- `rms_on="display"` → `_process_spectrum()` 後の x, y

## 13.7 `compute_rms_win`

```python
compute_rms_win(
    *,
    x_raw, y_resid,
    rms_on, smooth_mode, smooth_width, box_downsample, box_policy,
    xaxis, vel_disp, freq_ghz, chan,
    xr_vel, win_str, rest_hz_meta, rest_hz_user,
    has_bsl,
    use_fit_windows=True,
) -> Tuple[rms, n, mask]
```

### 役割
最終的な RMS(win) を計算します。

### 振る舞い
- `xr_vel` があればまずその範囲に制限
- `has_bsl` かつ fit window が使えるなら
  - `(xr_vel) AND (fit windows)` 上で RMS
- fit window がなければ
  - `xr_vel` のみで RMS

### 返り値
- `rms`
- `n`
  - 使用サンプル数
- `mask`
  - 最終 mask

---

## 14. `utils.py` の詳細

`utils.py` は多機能なので、用途別に整理します。

## 14.1 軸・数値補助

### `AxisBundle.build(meta, nchan)`
固定長スペクトル向けに、WCS から channel / frequency 軸束を作ります。  
ただし現 viewer/montage/grid は row ごと軸生成へ寄っており、補助的な位置づけです。

### `_norm_range(r)`
`(min, max)` を正規化。順序が逆なら入れ替える。

### `_interp_extrap(x, y, xq)`
NaN を無視しつつ線形補間・外挿を行う。

### `_vel_from_freq_ghz`, `_freq_ghz_from_vel`
radio definition の周波数・速度変換。

### `recalculate_velocity_axis(vel_old, rest_old, rest_new)`
`rest_freq` 上書き時の速度再解釈に使う最重要関数。

## 14.2 平滑化・ダウンサンプル

### `_running_mean(y, w)`
NaN 対応の moving average。

### `_box_average_downsample(x, y, w, policy="trim")`
箱平均しつつ点数も減らす。

- `policy="trim"`
  - 端数を落とす
- `policy="pad"`
  - NaN パディング後に平均

### `_process_spectrum(x, y, smooth_mode, smooth_width, box_downsample, box_policy)`
viewer/montage/grid 共通の平滑化入口。

- `"none"`
- `"running"`
- `"boxcar"`

を解釈します。

## 14.3 baseline 関連

### `_parse_list_like(obj)`
`BSL_COEF` を ndarray / list / 文字列 / scalar から配列へ解釈。

### `_fit_baseline_poly(...)`
row に入った baseline 多項式係数を velocity 上で評価する補助。

### `_baseline_extrapolated_line(...)`
baseline の表示用再構成を行う補助関数。

### `_build_fitwin_mask_on_vel(vel, win_str)`
velocity 上に fit window mask を作る。

### `in_any_windows`, `subtract_windows`
window 幾何処理の補助。  
`subtract_windows()` はコメントにもある通り簡略実装寄りです。

## 14.4 x 軸範囲変換

### `_convert_vel_range_to_xrange(...)`
display velocity 範囲を current xaxis 範囲へ変換。

### `_convert_user_xrange_to_vel_range(...)`
ユーザーが与えた `xrange` を velocity 範囲へ戻す。  
これにより x 軸を後から切り替えても「同じ物理範囲」が保たれます。

## 14.5 PDF と紙面レイアウト

### `parse_figsize(size_spec, default=None)`
図サイズ指定を解釈。

受理例:

- `(10, 6)`
- `"A4"`
- `"A4L"`
- `"A3P"`
- `"LETTER"`

戻り値は `(fig_size, bbox_mode)` で、固定紙面モードか通常モードかも決まります。

### `drive_pdf_generation(...)`
ページ送りを行いながら PDF を順次保存する共通ルーチン。

### `parse_paper_margins(spec, default="20mm,15mm,20mm,30mm")`
紙面余白を inch へ変換。  
`(left, right, bottom, top)` 順です。

### `parse_content_aspect(spec)`
内容領域の縦横比 `width/height` を解釈。

受理例:

- `"4:3"`
- `"1:1"`
- `1.333333`
- `(4, 3)`

### `paper_inner_rect_frac(...)`
紙サイズと余白から、内容領域の figure fraction を返す。

### `place_single_axes_in_rect(...)`
単一 Axes を与えられた矩形へ最大化して配置。

### `place_grid_axes_in_rect(...)`
複数 Axes を grid として矩形へ配置。  
`montage` の固定紙面モードで重要です。

### `axes_grid_bbox_frac(...)`
Axes 群の bounding box を figure fraction で返す。

### `adjacent_label_positions_from_bbox(...)`
内容 bbox に隣接するタイトル・xlabel・ylabel の位置を計算する。  
`grid` の固定紙面出力で重要です。

---

## 15. `grid`, `montage`, `viewer` で共通に使われる内部 helper

これらは公開 API ではありませんが、挙動理解に重要です。

## 15.1 `_filter_scantable_by_rows`
row 選択、`ORIG_ROW_ID` 保存、history 追記。

## 15.2 `_merge_scantables`
複数 Scantable を縦結合し、data は VLA list に統一。

## 15.3 `_row_specsys`
row または meta から実効的な `SPECSYS` / `SSYSOBS` を読む。

## 15.4 `_row_meta_for_axis`
row にある WCS / spectral metadata を header meta へ上書きして、1 row 用メタデータを作る。

## 15.5 `_row_rest_hz_meta`
row / meta から `RESTFRQ` または `RESTFREQ` を有効値として取得。

## 15.6 `_row_vcorr_kms`
row / meta から LSRK 補正速度を解決。

優先順:

1. `VELOSYS`, `VFRAME` [m/s]
2. `V_CORR_KMS`, `v_corr_kms` [km/s]

ただし `SPECSYS/SSYSOBS` が明示的に TOPO 以外なら追加補正しません。

---

## 16. 実装上の注意点・癖

## 16.1 `coord="auto"` は使えない
grid では `auto` を指定すると例外です。

## 16.2 `combine` は将来用に近い
現状は `mean` 相当固定です。

## 16.3 `show_rms` は未使用
grid の `show_rms` は API に残っていますが、現実装での効きはありません。

## 16.4 `projection` の解釈上の注意
x オフセットに掛ける cos 因子の選択が主用途です。

## 16.5 RA の unit 推定に自動変換がある
RA 値が 24.5 以下なら hour 単位と推定される可能性があります。  
誤解を避けるには degree で統一するのが安全です。

## 16.6 `grid` の coadd は参照 row 軸へ補間
セル内 coadd の再サンプリング先は「最初の row の velocity 軸」です。

## 16.7 `BSL_COEF` の解釈は 2 通り
ベクトル baseline と多項式係数の両対応です。

---

## 17. cookbook（実務的な使い方）

## 17.1 まず 1 本ずつ確認して怪しい row を探す

```python
from sd_radio_spectral_fits.profile_view import view_spectra

view_spectra(
    sc,
    xaxis="vel",
    xrange=(-30, 50),
    show_fit_windows=True,
    show_fitwin_rms=True,
)
```

用途:
- baseline の入り方を確認
- `BSL_WINF` が妥当かを見る
- `rest_freq` のズレを疑う前の生確認

## 17.2 `rest_freq` を変えて線同定を確かめる

```python
view_spectra(
    sc,
    xaxis="vel",
    rest_freq=110.201354e9,
    xrange=(-20, 20),
    show_fit_windows=True,
)
```

用途:
- header の rest frequency と別の遷移で見たい
- 速度中心がどれだけずれるか確認したい

## 17.3 多数 row を一気に眺める

```python
from sd_radio_spectral_fits.profile_view import ProfileMapMontageViewer

ProfileMapMontageViewer(
    sc,
    imin=0,
    imax=80,
    nrows=4,
    ncols=5,
    xaxis="vel",
    xrange=(-20, 40),
    annotate_rms=True,
)
```

用途:
- mapping データ全体の quality check
- 突然おかしい baseline や雑音 row を探す

## 17.4 A4 横固定で montge PDF を作る

```python
ProfileMapMontageViewer(
    sc,
    nrows=4,
    ncols=5,
    show=False,
    save_pdf="montage_A4L.pdf",
    figsize="A4L",
    paper_margins="12mm,10mm,15mm,18mm",
    content_aspect="4:3",
)
```

用途:
- 配布資料
- 研究室内の確認用 PDF

## 17.5 RA/Dec の profile grid を作る

```python
from sd_radio_spectral_fits.profile_view import ProfileMapGridViewer

ProfileMapGridViewer(
    sc,
    coord="radec",
    projection="SFL",
    ref_point=(83.82208, -5.39111),
    x_grid=30.0,
    y_grid=30.0,
    mode="all",
    xaxis="vel",
    xrange=(-20, 40),
    show_coord_text=True,
    show_fit_windows=True,
)
```

用途:
- profile map を天球配置で直感的に見る
- セルごとのスペクトル分布を把握する

## 17.6 観測設計通りの 0,0 基準格子に固定する

```python
ProfileMapGridViewer(
    sc,
    coord="radec",
    ref_point=(83.82208, -5.39111),
    x_grid=15.0,
    y_grid=15.0,
    grid_anchor_offsets=(0.0, 0.0),
    grid_bounds_offsets=((-90.0, -90.0), (90.0, 90.0)),
    mode="coadd",
)
```

用途:
- 複数実験日で同じ格子に載せて比較
- ビーム配置や pointing offset を反映した profile map の見比べ

## 17.7 銀河座標で見る

```python
ProfileMapGridViewer(
    sc,
    coord="gal",
    ref_point=(0.0, 0.0),
    projection="SFL",
    x_grid=60.0,
    y_grid=60.0,
    offset_unit="arcmin",
)
```

用途:
- 銀経銀緯で見た方が自然な広域分布

## 17.8 限られた row だけ再確認する

```python
view_spectra(sc, rows="10, 12, 15")
ProfileMapMontageViewer(sc, rows="100:140", nrows=4, ncols=5)
ProfileMapGridViewer(sc, exclude_rows="0:5", coord="radec", ref_point=(83.82208, -5.39111))
```

---

## 18. どの viewer を選ぶべきか

- **個々の row の基線・速度軸・RMS を細かく見る**
  - `SpectralViewer`
- **多くの row を順番にざっと見る**
  - `ProfileMapMontageViewer`
- **空間配置とスペクトル形状を同時に見る**
  - `ProfileMapGridViewer`

実務上は、

1. `view_spectra()` で問題 row を掴む
2. `ProfileMapMontageViewer` で全体傾向を見る
3. `ProfileMapGridViewer` で空間配置を確認する

の順が使いやすいです。

---

## 19. まとめ

`profile_view` は、単なる表示ツールではなく、

- row ごとの spectral WCS
- `rest_freq` 上書き
- `BSL_WINF` の表示変換
- `TEMPSCAL/BEAMEFF` によるスケール変換
- VLA 対応
- grid への座標量子化
- 固定紙面 PDF 出力

を一体化した実務ツール群です。

特に `grid` は、一般的な画像 regridding とは違い、

- ref_point を決め
- 座標差を arcsec 化し
- `x_grid`, `y_grid` で丸め
- `grid_anchor_offsets` と `grid_bounds_offsets` で格子を固定し
- `grid_tol` で採否を決める

という **明快だが実装依存の manual gridding** を採っています。  
この性格を理解して使うと、観測データの profile map 確認には非常に有用です。
