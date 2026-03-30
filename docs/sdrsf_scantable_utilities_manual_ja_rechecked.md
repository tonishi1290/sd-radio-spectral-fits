# Scantable Utilities マニュアル
## `scantable_utils.py` / `tempscale.py` 実装準拠ガイド（最新版）

## 0. この文書について

本書は、最新版の `scantable_utils.py` と `tempscale.py` を実装ベースで説明する補助説明書です。  
対象は、Scantable の検査・部分抽出・オフセット計算・結合・補助更新、`BEAMEFF` 付与、および data 本体への相対/一律スケール補正です。

特に今回は、次の2点を重点的に明文化しています。

1. `calc_mapping_offsets()` が**何を入力として、どの式で、どの座標差を返すか**
2. `set_beameff()` と `apply_relative_scale()` / `apply_global_scale()` を**どう使い分けるべきか**

---

## 1. このモジュール群の役割

### 1.1 `scantable_utils.py`

`scantable_utils.py` は、Scantable の

- 列確認
- 行抽出
- 条件検索
- マッピング用オフセット計算
- 複数ファイル結合
- metadata / table 列の補助更新

を行う実務ユーティリティです。

最新版コードの主要機能は次です。

```python
describe_columns(sc)
show_scantable(inputs, rows=None, columns="default", head=20, show_legend=False,
               extra_data=None, ref_coord=None, frame="ICRS", projection="GLS", unit="arcsec")
calc_mapping_offsets(sc, ref_coord=None, frame="ICRS", projection="GLS", unit="arcsec",
                     cos_mode="point", verbose=True)
merge_scantables(inputs, sort_by_time=False, shift_scan_id=True)
update_metadata(sc, column, value, rows=None, force=False, verbose=True)
set_beameff(sc, efficiency, rows=None, key_columns=("FDNUM", "IFNUM", "PLNUM"), strict=False, verbose=True)
apply_relative_scale(sc, scale, rows=None, key_columns=("FDNUM", "IFNUM", "PLNUM"), strict=False, invalidate_derived=True, verbose=True)
apply_global_scale(sc, factor, rows=None, invalidate_derived=True, verbose=True)
find_scans(sc, query=None, extra_data=None, **kwargs)
filter_scantable(sc, query=None, extra_data=None, rows=None, **kwargs)
select_rows(sc, rows)
```

### 1.2 `tempscale.py`

`tempscale.py` は、温度スケールと `BEAMEFF` の意味論を担う補助モジュールです。  
主な役割は

- `TEMPSCAL` の正規化
- `BEAMEFF` / `TEMPSCAL` の row-wise 解決
- `TA*` と `TR*` の変換
- ragged / VLA データに対する row-wise 変換
- scale history の追記

です。

最新版では `tempscale.py` に `set_beameff()`、`apply_relative_scale()`、`apply_global_scale()` が定義され、`scantable_utils.py` はそれらの薄い wrapper を提供します。

---

## 2. インポート

### 2.1 もっとも確実な形

```python
import sd_radio_spectral_fits.scantable_utils as su
import sd_radio_spectral_fits.tempscale as ts
```

### 2.2 package top-level について

最新版の `__init__.py` では、次は top-level から再公開されています。

- `show_scantable`
- `describe_columns`
- `update_metadata`
- `merge_scantables`
- `filter_scantable`
- `find_scans`
- `calc_mapping_offsets`

最新版の `__init__.py` では、`set_beameff`、`apply_relative_scale`、`apply_global_scale` も top-level から再公開されます。したがって、

```python
import sd_radio_spectral_fits as sd
sd.set_beameff(...)
sd.apply_relative_scale(...)
sd.apply_global_scale(...)
```

でも使えます。もちろん、実務では `scantable_utils` 経由でも構いません。

### 2.3 どの入口を使うべきか

実装本体は `tempscale.py` にあり、`scantable_utils.py` は utility 入口としてそれを再公開します。通常は

- package top-level を使うなら `sd.set_beameff(...)`, `sd.apply_relative_scale(...)`, `sd.apply_global_scale(...)`
- utility モジュールをまとめて使うなら `su.set_beameff(...)`, `su.apply_relative_scale(...)`, `su.apply_global_scale(...)`

のどちらでもよいです。

---

## 3. まず押さえるべき設計思想

### 3.1 `table`, `data`, `meta` の役割分担

Scantable は大きく次の3つを持ちます。

- `sc.table`
  - 行ごとの表データ
- `sc.data`
  - 実スペクトル配列
- `sc.meta`
  - グローバル header 相当

`scantable_utils.py` は主に `table` と `meta` を扱います。  
ただし `filter_scantable()`, `select_rows()`, `merge_scantables()` は `data` も対応行だけ追随させます。

### 3.2 非破壊関数と破壊的関数

#### 非破壊
- `describe_columns()`
- `show_scantable()`
- `calc_mapping_offsets()`
- `find_scans()`
- `filter_scantable()`
- `select_rows()`
- `merge_scantables()`

#### 破壊的
- `update_metadata()`
- `set_beameff()`

### 3.3 `rows` は 0-based positional index

多くの関数の `rows` は、表示 index ではなく **0-based positional index** です。

受理される典型例:

- `None`
- `3`
- `slice(0, 10)`
- `[0, 5, 10]`
- `"0:10"`
- `"0:10,50,60:80"`

`"0:10"` は Python slice と同じで stop を含みません。

### 3.4 Big-endian FITS 対策

最新版には `_df_to_native_endian()` があり、FITS 由来の big-endian DataFrame が Pandas の `iloc`, `query`, `isin`, sort などで落ちる問題に対処しています。  
この保護は少なくとも

- `show_scantable()`
- `find_scans()`
- `filter_scantable()`
- `merge_scantables()`

で効いています。

### 3.5 `TIMESTAMP` の自動解決

`show_scantable()` や `merge_scantables(sort_by_time=True)` では、必要に応じて時刻列を内部解決します。優先順位は

1. `TIMESTAMP`
2. `MJD`
3. `DATE-OBS` / `DATEOBS` と `TIME`
4. 旧 `TIME`
5. DatetimeIndex

です。

---

## 4. `show_scantable()`

## 4.1 何をするか

Scantable の中身を表形式で確認します。

### 主な特徴
- 入力が `Scantable` でもファイルパスでもよい
- `rows` で位置インデックス指定できる
- table に無い列でも、alias や meta から拾えることがある
- `TIMESTAMP` が無ければ内部的に補う
- `extra_data` や `calc_mapping_offsets()` の結果をその場で一時結合できる

## 4.2 alias / meta 検索

たとえば `RESTFREQ` を要求したが table に無い場合、

1. table の `RESTFREQ`
2. table の alias `RESTFRQ`
3. meta の `RESTFREQ`
4. meta の alias `RESTFRQ`

の順で探します。

したがって、列揺れ確認に便利です。

## 4.3 実用例

```python
su.show_scantable(
    inputs="tastar_dump.fits",
    rows="0:10",
    columns="SCAN, TIMESTAMP, OBSMODE, RESTFRQ, RA, DEC, VELOSYS, TEMPSCAL, BEAMEFF",
    head=100,
)
```

全列をざっと見る場合:

```python
su.show_scantable(sc, rows="0:5", columns="all", head=5)
```

---

## 5. `calc_mapping_offsets()`

## 5.1 何をする関数か

`calc_mapping_offsets()` は、Scantable の各 row について、**ある参照点から見た相対オフセット**を計算し、

- `OFS_LON`
- `OFS_LAT`

の 2 列を持つ `DataFrame` を返す関数です。

ここで重要なのは、これは

- スペクトル data を変換する関数ではない
- `Scantable` 自体に列を追加する関数でもない
- WCS 再投影器でもない

という点です。

あくまで、

> 「各 row の座標が、参照点から見てどれだけずれているか」

を、指定した frame / projection / unit に従って数値化して返す関数です。

したがって用途は主に

- mapping 観測の中心からの相対位置を確認する
- `show_scantable()` に一時的に表示する
- `find_scans()` / `filter_scantable()` の `extra_data` と組み合わせて、中心付近だけ抽出する
- scan pattern の sanity check をする

です。

## 5.2 入力座標は何を使うか

最新版コードでは、row から次を探します。

### 優先順
1. `RA`, `DEC` がある
   - `src_frame = "icrs"`
2. `GLON`, `GLAT` がある
   - `src_frame = "galactic"`
3. 次のいずれかがある
   - `AZIMUTH`, `ELEVATIO`
   - `AZIMUTH`, `ELEVATION`
   - `AZ`, `EL`
   - `src_frame = "altaz"`

つまり、**RA/DEC があればそれを最優先**します。

## 5.3 参照点 `ref_coord`

### `ref_coord` を与えた場合
- `SkyCoord` ならそのまま使う
- 文字列なら `SkyCoord(ref_coord)` で解釈する

### `ref_coord=None` の場合
- `src_frame == "icrs"` かつ `meta` に `OBSRA`, `OBSDEC` があればそれを使おうとする
- 失敗したら、row の座標から
  - 経度側は circular mean
  - 緯度側は単純平均

を使います。

### 重要
`ref_coord="Orion KL"` のような名前解決は、最新版コードでは保証されません。  
確実なのは

- 座標文字列
- `SkyCoord`

です。

例:

```python
ref_coord="05h35m14.5s -05d22m30s"
```

## 5.4 frame 変換

### RA/DEC または GLON/GLAT の場合
`SkyCoord(...).transform_to(tgt_frame)` を使って target frame へ変換します。

`frame` 文字列には簡単な別名マップがあります。

- `"j2000"` -> `"fk5"`
- `"b1950"` -> `"fk4"`
- `"lb"` -> `"galactic"`

### AltAz の場合
AltAz は、site / time が無いと他 frame へ変換できません。  
最新版コードでは、AltAz 入力のとき

- `frame="AltAz"` 相当なら差分計算のみ許可
- それ以外の frame 変換は `ValueError`

です。

したがって、AltAz row しか無い Scantable に対して `frame="ICRS"` を指定して mapping offset を出す、という用途にはそのまま使えません。

## 5.5 実際に計算している量

変数を明示します。

- row 座標: $(\ell, b)$
- 参照点: $(\ell_0, b_0)$

まず経度差は 0/360 境界を跨いでも破綻しにくいように

$$
\Delta \ell = ((\ell - \ell_0 + 180) \bmod 360) - 180
$$

とします。

緯度差は

$$
\Delta b = b - b_0
$$

です。

### `projection="GLS"` または `"SFL"`

このとき最新版コードは

$$
\mathrm{OFS\_LON} = \Delta \ell \cos b_\ast
$$

$$
\mathrm{OFS\_LAT} = \Delta b
$$

を使います。

ここで $b_\ast$ は `cos_mode` で決まります。

- `cos_mode="point"`
  - $b_\ast = b$
- `cos_mode="ref"`
  - $b_\ast = b_0$

つまり、`GLS/SFL` の違いは現コード上はありません。  
どちらも同じ分岐です。

### `projection="CAR"` または `"NONE"`

このときは

$$
\mathrm{OFS\_LON} = \Delta \ell
$$

$$
\mathrm{OFS\_LAT} = \Delta b
$$

です。  
経度方向の $\cos b$ 補正をしません。

## 5.6 単位変換

最後に degree から

- `arcsec`
- `arcmin`
- `deg`

へ換算します。

たとえば `unit="arcsec"` なら scale は $3600$ です。

したがって最終出力は、たとえば `GLS` なら

$$
\mathrm{OFS\_LON[arcsec]} = \Delta \ell \cos b_\ast \times 3600
$$

$$
\mathrm{OFS\_LAT[arcsec]} = \Delta b \times 3600
$$

です。

## 5.7 何をしていないか

最新版コードに即して、明確に書くと次は**していません**。

- `TAN` 投影
- `SIN` 投影
- full WCS reprojection
- 画像の回転や PC/CD 行列処理
- `Scantable.table` への永続的列追加
- 名前解決による天体辞書参照
- AltAz -> ICRS の site/time 無し変換

docstring 中には `TAN` や `SIN` の説明が残っていますが、**現実装は `GLS/SFL` と `CAR/NONE` だけ**です。ここは manual 側で実装優先で読むべきです。

## 5.8 具体例

### 参照点まわりの offset を DataFrame で得る

```python
ofs = su.calc_mapping_offsets(
    sc=sc,
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
    cos_mode="ref",
)
print(ofs.head())
```

### 中心 60" 以内だけ抜く

```python
ofs = su.calc_mapping_offsets(
    sc,
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
)

sc_center = su.filter_scantable(
    sc,
    query="OFS_LON**2 + OFS_LAT**2 < 60**2",
    extra_data=ofs,
)
```

### `show_scantable()` に直接表示

```python
su.show_scantable(
    inputs=sc,
    rows="0:20",
    columns="TIMESTAMP, RA, DEC, OFS_LON, OFS_LAT",
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
)
```

## 5.9 `cos_mode="point"` と `"ref"` の使い分け

小さなマップでは差は小さいですが、考え方は違います。

### `point`
各 row の緯度 $b$ で $\cos b$ を掛けます。

$$
x_i = \Delta \ell_i \cos b_i
$$

### `ref`
参照点の緯度 $b_0$ で統一します。

$$
x_i = \Delta \ell_i \cos b_0
$$

小領域近似としては `ref` の方が「一貫した 1 枚の平面座標」として解釈しやすいです。  
一方、現行コードの既定値は `point` です。

実務上は、

- 小さい map で従来互換重視 -> `point`
- 中心基準の平面座標として使いたい -> `ref`

と考えるとわかりやすいです。

---

## 6. `merge_scantables()`

## 6.1 何をするか

複数の Scantable またはファイルパスを結合し、新しい 1 つの Scantable を返します。

### 特徴
- fixed-length / VLA 混在に対応
- table は native endian 化してから結合
- `shift_scan_id=True` が既定
- `sort_by_time=True` なら内部解決時刻でソート

## 6.2 `shift_scan_id=True` の意味

複数ファイルの `SCAN` が重複すると、後段で `group_mode="scan"` を使う際に危険です。  
そこで最新版では、既定で各ファイルの `SCAN` にオフセットを足して重複を防ぎます。

これはマッピングや複数日結合ではかなり重要です。

## 6.3 実用例

```python
sc_merged = su.merge_scantables(
    inputs=["day1.fits", "day2.fits", "day3.fits"],
    sort_by_time=True,
    shift_scan_id=True,
)
```

---

## 7. `update_metadata()`

## 7.1 何をするか

指定列または alias 列を破壊的に更新します。  
主用途は

- `OBJECT`
- 補助列
- 限定的な `RESTFRQ` 修正

です。

## 7.2 危険列

最新版には `DANGER_COLS` があり、少なくとも

- `CRVAL1`, `CDELT1`, `CRPIX1`
- `CTYPE1`, `CUNIT1`
- `RESTFREQ`, `RESTFRQ`
- `SPECSYS`
- `DATA`, `FLAG`, `SPECTRUM`

は危険列として扱われます。

これらを触るときは、何を変えるかを明確にすべきです。

---

## 8. `find_scans()`, `filter_scantable()`, `select_rows()`

## 8.1 `find_scans()`

条件に合う row の **0-based positional index** を返します。

### 条件指定
- `query=...`
- `kwargs`
- `extra_data`

を組み合わせられます。

```python
idx = su.find_scans(sc, OBSMODE="ON")
```

```python
idx = su.find_scans(sc, query="SCAN > 10 and OBSMODE == 'ON'")
```

```python
ofs = su.calc_mapping_offsets(sc, ref_coord="05h35m14.5s -05d22m30s")
idx = su.find_scans(
    sc,
    query="OFS_LON**2 + OFS_LAT**2 < 120**2",
    extra_data=ofs,
)
```

## 8.2 `filter_scantable()`

`find_scans()` 相当の条件で、新しい Scantable を返します。  
`table` だけでなく `data` も対応行だけに切り詰めます。

```python
sc_on = su.filter_scantable(sc, OBSMODE="ON")
```

```python
sc_scan = su.filter_scantable(sc, query="SCAN >= 10 and SCAN <= 20")
```

## 8.3 `select_rows()`

単純に位置インデックスで切り出すショートカットです。

```python
sc_head100 = su.select_rows(sc, "0:100")
```

---

## 9. `set_beameff()` / `apply_relative_scale()` / `apply_global_scale()`

## 9.1 3関数の役割分担

### `set_beameff()`
`BEAMEFF` を row ごとに設定する **metadata setter** です。data 本体は変えず、`TEMPSCAL` も変えません。

受け付ける `efficiency` の形は次の 3 通りです。

1. scalar
2. 1 次元 sequence / ndarray
3. mapping

mapping を使うと、たとえば `FDNUM / IFNUM / PLNUM` ごとに `BEAMEFF` を一括設定できます。

なお `BEAMEFF` は物理量なので、実装では

$$
0 < \eta \le 1
$$

を満たさない値は error になります。

### `apply_relative_scale()`
相対スケール係数 $g_i$ を row ごとに data へ直接**掛ける**関数です。

$$
T_i^{\mathrm{aligned}}(k) = g_i\,T_i(k)
$$

マルチビームで、まず全ビームの強度スケールを一致させたいときに使います。

### `apply_global_scale()`
全 row へ一律係数 $f$ を data へ直接**掛ける**関数です。

$$
T_i^{\mathrm{final}}(k) = f\,T_i^{\mathrm{aligned}}(k)
$$

後日、絶対スケール係数が決まったときに全体へ一括適用する用途です。

`apply_global_scale()` は scalar 専用です。row ごとに異なる係数を与えたい場合は `apply_relative_scale()` を使います。

## 9.2 なぜ `set_beameff()` と相対スケール補正を分けるのか

`BEAMEFF` は物理的には

$$
T_R^* = \frac{T_A^*}{\eta}
$$

の $\eta$ として使われます。つまり、`apply_relative_scale()` / `apply_global_scale()` が **掛ける**係数であるのに対し、`BEAMEFF` は `TA* -> TR*` 変換では **割る**側に入ります。

したがって、beam 間相対補正のための係数を `BEAMEFF` に流用すると、viewer / coadd / write が自動で `TR*` 変換に使ってしまい、意味が混線します。

そのため最新版では、

- 物理 `BEAMEFF` を設定する `set_beameff()`
- data 本体へ相対係数を掛ける `apply_relative_scale()`
- data 本体へ一律係数を掛ける `apply_global_scale()`

を分離します。

## 9.3 実装場所と公開場所

実装本体は `tempscale.py` にあります。`scantable_utils.py` 側の 3 関数は薄い wrapper です。

- `tempscale.set_beameff()`
- `tempscale.apply_relative_scale()`
- `tempscale.apply_global_scale()`

`scantable_utils.py` では同名 wrapper があり、さらに最新版の `__init__.py` では top-level からも再公開されます。したがって、次のどちらでも使えます。

```python
import sd_radio_spectral_fits as sd
sd.set_beameff(...)
sd.apply_relative_scale(...)
sd.apply_global_scale(...)
```

```python
import sd_radio_spectral_fits.scantable_utils as su
su.set_beameff(...)
su.apply_relative_scale(...)
su.apply_global_scale(...)
```

## 9.4 `apply_relative_scale()` / `apply_global_scale()` の重要仕様

- data 本体を **in-place** に更新する
- `TEMPSCAL` は変更しない
- 累積係数を `INTSCALE` 列へ記録する
- 既知のスケール依存派生列を既定で無効化する

既知の無効化対象は少なくとも次です。

- `BSL_RMS`
- `BSL_COEF`
- `BSL_SCALE`
- `MOMENT0`
- `MOM0`
- `INTEGRATED`

重要なのは、`apply_global_scale()` に物理 `TR*` 変換係数 $1/\eta$ をそのまま入れても、**`TEMPSCAL` は自動では `TR*` に変わらない**という点です。物理的な `TA* -> TR*` を後段で使いたい場合は、`set_beameff()` で `BEAMEFF=\eta` を設定し、viewer / coadd / write の `TR*` 変換へ渡す方が意味論として安全です。

## 9.5 係数入力の厳密な意味

### `set_beameff()`
`efficiency` は次のいずれかです。

- scalar
- 1 次元 sequence / ndarray
- mapping

sequence / ndarray の場合、長さは **選択された row 数と完全一致**しなければなりません。

mapping の場合、`key_columns` で決まる key ごとに値を引きます。たとえば

```python
key_columns=("FDNUM", "IFNUM", "PLNUM")
```

なら key は

$$
(FDNUM, IFNUM, PLNUM)
$$

です。単一列 key の場合は、mapping key は tuple でなくても構いません。

`strict=False` では、mapping 側に未使用 key があっても warning です。

`strict=True` では、未使用 key があると error になります。

一方、**選択された row に対応する key が mapping に存在しない**場合は、`strict` に関係なく error です。

### `apply_relative_scale()`
`scale` も同様に

- scalar
- 1 次元 sequence / ndarray
- mapping

を受けます。意味は `set_beameff()` と同じですが、こちらは data 本体へ直接掛けます。

### `apply_global_scale()`
`factor` は scalar だけです。

## 9.6 いつ適用すべきか

一番安全なのは、**calibration 直後**です。

1. `run_tastar_calibration()` で `Ta*` を作る
2. 必要なら `apply_relative_scale()` で beam 間相対補正
3. generic な一律係数が必要なら `apply_global_scale()` を掛ける
4. 物理 `TR*` 変換に必要な $\eta$ が分かっているなら、`set_beameff()` を設定する
5. その後に baseline / RMS / moment / coadd

この順なら、後段の派生量はすべて現在のスケールで計算されます。

実務上は、`set_beameff()` は data を変えないので、step 4 は calibration 直後でも coadd 前でも構いません。ただし `TR*` を viewer / coadd / write で使いたい時点までには設定しておく必要があります。

## 9.7 `set_beameff()` の使い方

### 全行へ同じ効率

```python
sd.set_beameff(sc, 0.42)
```

### 行ごと配列で設定

```python
sd.set_beameff(sc, [0.41, 0.42, 0.43], rows=[0, 1, 2])
```

### `FDNUM / IFNUM / PLNUM` ごとの mapping 指定

```python
beam_eff_map = {
    (0, 0, 0): 0.46,
    (0, 0, 1): 0.44,
    (1, 0, 0): 0.42,
    (1, 0, 1): 0.41,
}

sd.set_beameff(
    sc,
    beam_eff_map,
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
)
```

### 未使用 key も error にしたい

```python
sd.set_beameff(
    sc,
    beam_eff_map,
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
    strict=True,
)
```

## 9.8 `apply_relative_scale()` の使い方

```python
rel_map = {
    (0, 0, 0): 1.00,
    (0, 0, 1): 0.97,
    (1, 0, 0): 1.03,
    (1, 0, 1): 1.01,
}

sd.apply_relative_scale(
    sc,
    rel_map,
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
)
```

### 一部 row だけに sequence を掛ける

```python
sd.apply_relative_scale(
    sc,
    [1.00, 0.98, 1.02],
    rows=[10, 11, 12],
)
```

## 9.9 `apply_global_scale()` の使い方

```python
sd.apply_global_scale(sc, 0.91)
```

### 派生量 invalidation を止めたい場合

```python
sd.apply_global_scale(sc, 0.91, invalidate_derived=False)
```

ただし、通常は `invalidate_derived=True` の既定のままが安全です。

## 9.10 高速化上の注意

- 2D `numpy.ndarray` の場合は選択 row に対してベクトル化して掛ける
- VLA / list-of-arrays の場合は row ループになるが、各 row の掛け算自体は NumPy 配列上で行う
- mapping 指定の係数解決は、単一 key なら `map`、複数 key なら `merge` ベースで処理する
- したがって、beam 間相対補正は **baseline や moment 計算より前**にまとめて一度だけ行うのが最も効率的です

VLA では row ごとに長さが違うため、2D ndarray のような全面 broadcasting は本質的に使いにくいです。そのため VLA は 2D ndarray より遅くなりやすいですが、最新版実装では「係数 lookup の高速化」を先に行っているので、実用上はかなり軽くなっています。

## 10. `tempscale.py` と `set_beameff()` の関係

`tempscale.py` の本質は、`BEAMEFF` を使って

$$
T_R^* = \frac{T_A^*}{\eta}
$$

$$
T_A^* = T_R^* \eta
$$

を行うことです。

重要なのは、`set_beameff()` は**この変換を実行する関数ではない**ということです。  
変換自体は

- viewer
- coadd
- write
- 明示的変換 helper

が担います。

`set_beameff()` は、**変換に必要な $\eta$ を row ごとに準備するだけ**です。  
この責務分離は正しいです。

---

## 11. 現場でのおすすめ順序

実務では次の流れが安定です。

1. `describe_columns()` / `show_scantable()` で列確認
2. `find_scans()` / `filter_scantable()` / `select_rows()` で小さい試験セットを作る
3. 必要なら `calc_mapping_offsets()` で位置 sanity check
4. 複数日なら `merge_scantables(shift_scan_id=True)`
5. calibration 後、必要なら `apply_relative_scale()`
6. 必要なら `apply_global_scale()`
7. 物理 `TR*` を使うなら `set_beameff()`
8. baseline / RMS / moment / coadd / viewer / write

`set_beameff()` は metadata setter なので、step 7 は少し前後しても構いません。ただし `TR*` を downstream で解釈する前には設定しておく必要があります。


## 12. まとめ

### `calc_mapping_offsets()`
- 参照点からの相対 offset を返す
- immutable
- `GLS/SFL` と `CAR/NONE` が実装本体
- `TAN/SIN` は docstring に残っていても現実装では使えない
- mapping 観測の sanity check と位置抽出に非常に有用

### `set_beameff()`
- row ごとの `BEAMEFF` を設定するだけ
- data も `TEMPSCAL` も変えない
- downstream の `TR*` 利用準備として重要
- マルチビーム解析でも十分便利
- ただし大量の beam 群を扱うには convenience helper があるとさらに良い
