# 移動天体 OTF mapping 使用説明書: Moon/Sun offset, Moon disk, SPICE selenographic

この文書は、`sd_radio_spectral_fits.map_3d` に追加した移動天体用の OTF mapping 座標系をまとめたものです。月を `dAz, dEl` で観測した場合、太陽を `dAz, dEl` で観測した場合、月を見かけ円盤として回転補正して表示する場合、SPICE で月面固定座標へ投影する場合をカバーします。

## 1. 目的

通常の RA/Dec または Galactic mapping では、観測対象が固定天球上にあることを仮定します。一方、月や太陽は観測中に RA/Dec 上を移動します。そのため、月や太陽に対する `dAz, dEl` OTF 観測を通常の固定 RA/Dec map として扱うと、対象の移動や視野回転により、表示・合成・basket weave の解釈が不安定になります。

今回の修正では、次の座標系を追加しました。

| `coord_sys` | 対象 | 目的 | 軸 | 主な用途 |
|---|---:|---|---|---|
| `moon_azel_offset` | Moon | 月中心からの topocentric Az/El offset | `dAz cos El`, `dEl` | 生の観測パターン確認 |
| `sun_azel_offset` | Sun | 太陽中心からの topocentric Az/El offset | `dAz cos El`, `dEl` | 太陽 scan の観測確認 |
| `moon_sky_offset` | Moon | 月中心の天球接平面 offset | `dRA cos Dec`, `dDec` | AltAz 視野回転を除いた月円盤表示 |
| `sun_sky_offset` | Sun | 太陽中心の天球接平面 offset | `dRA cos Dec`, `dDec` | 太陽円盤を天球北基準で表示 |
| `moon_disk_offset` | Moon | 月円盤の向きを月北方向に固定 | `Moon disk X`, `Moon disk Y` | 欠けた月・月面模様が回らない表示 |
| `moon_selenographic` | Moon | SPICE による月面固定座標 | `SELON`, `SELAT` | 月面経度・緯度での合成 |

通常の月 mapping 表示では、まず `moon_disk_offset` を推奨します。これは月を約 30 分角の見かけ円盤のまま表示し、かつ月の北方向が表示上で固定されるためです。

## 2. 入力データの前提

各 row、各 beam について、次の情報が必要です。

- 実際の beam 方向を表す `RA`, `DEC` または `GLON`, `GLAT`
- 各 dump の時刻: `TIMESTAMP`, `MJD`, `DATE-OBS` など
- 観測地: `OBSGEO-X/Y/Z` または `SITELAT`, `SITELONG`, `SITEELEV`
- スペクトルデータ本体

重要な点として、マルチビーム受信機では `BORE_AZ`, `BORE_EL` は座標計算に使いません。各 beam の `RA/DEC` または `GLON/GLAT` が正しく計算されている前提で、その beam 方向を各時刻・観測地で評価してから Moon/Sun 中心との差を取ります。

## 3. RA/DEC の解釈

移動天体では月の視差が大きいため、`RA`, `DEC` が何を表すかが重要です。

既定では、`RA`, `DEC` は **topocentric apparent direction** として扱います。これは、実際の encoder/corrected pointing からその時刻の視線方向を戻した場合に対応します。

必要なら scantable metadata で次を指定できます。

```python
scantable.meta["MOVING_BODY_RADEC_FRAME"] = "apparent"  # default
```

または、もし far-field の ICRS 座標として保存している場合は、明示的に次を指定します。

```python
scantable.meta["MOVING_BODY_RADEC_FRAME"] = "ICRS"
```

対応値は以下です。

- `apparent` または `topocentric_apparent`
- `ICRS`
- `FK5`
- `FK4`

月や太陽に対する実観測 beam 方向であれば、通常は `apparent` が安全です。

## 4. 各座標系の定義

### 4.1 `moon_azel_offset` / `sun_azel_offset`

各 dump の時刻 $t_i$ で、beam 方向と対象中心の Az/El を計算します。

- beam 方向: $Az_i$, $El_i$
- 対象中心: $Az_b(t_i)$, $El_b(t_i)$

`projection="SFL"` の場合、gridding に使う座標は

$$
x_i = \left[Az_i - Az_b(t_i)\right]_{\rm wrap}\cos El_b(t_i)
$$

$$
y_i = El_i - El_b(t_i)
$$

です。単位は内部では arcsec、FITS WCS では deg です。

この座標系は、観測時の `dAz, dEl` scan pattern をそのまま確認するのに向いています。ただし AltAz 座標なので、長時間観測では月面模様や欠け方が回って見えることがあります。

### 4.2 `moon_sky_offset` / `sun_sky_offset`

各 dump の時刻 $t_i$ で、対象中心を原点にした天球接平面 offset を作ります。

小角近似で書けば、

$$
x_i \simeq \left[\alpha_i - \alpha_b(t_i)\right]\cos\delta_b(t_i)
$$

$$
y_i \simeq \delta_i - \delta_b(t_i)
$$

です。実装では Astropy の spherical offset 相当を使い、小角近似より安全に計算します。

- `x`: RA 増加方向、つまり天球東方向
- `y`: Dec 増加方向、つまり天球北方向

この座標系は、AltAz の視野回転を除去できます。ただし軸は天球北・天球東に固定されるため、月そのものの北方向に固定されるわけではありません。

### 4.3 `moon_disk_offset`

`moon_disk_offset` は今回の月観測で最も実用的な表示用座標系です。月を約 30 分角の見かけ円盤として表示しつつ、月の北方向が常に同じ向きになるように回転します。

処理は次の順序です。

1. 月中心の `moon_sky_offset` を計算する。
2. 必要なら pointing 補正 `dx_arcsec`, `dy_arcsec` を sky offset 座標で加える。
3. 月の見かけ北方向 position angle $P$ を計算する。
4. sky offset を $P$ で回転し、moon disk offset を得る。

sky offset を

$$
(x_{\rm sky}, y_{\rm sky})
$$

とし、$P$ を「天球北から東回りに測った月北極方向の position angle」とすると、

$$
x_{\rm disk} = x_{\rm sky}\cos P - y_{\rm sky}\sin P
$$

$$
y_{\rm disk} = x_{\rm sky}\sin P + y_{\rm sky}\cos P
$$

です。

軸の意味は次の通りです。

- `Moon disk Y`: 見かけの月北方向。月中心から見て月の北極が投影される方向。
- `Moon disk X`: 月北方向から position angle で +90 度回した方向。

これは selenographic longitude/latitude ではありません。月面球への展開ではなく、見かけの 30 分角円盤の中で軸だけを月に固定する座標です。

月が欠けている場合、`dAz/dEl` 表示では欠けている方向が日周運動で回って見えます。`moon_disk_offset` では、この回転を抑え、月の向きを揃えたまま表示できます。

### 4.4 `moon_selenographic`

`moon_selenographic` は SPICE を用いて、各 beam の視線と月の基準楕円体との交点を計算し、その点を selenographic longitude/latitude へ変換します。

- `SELON`: east-positive selenographic longitude
- `SELAT`: planetocentric selenographic latitude

概念的には次です。

```text
observer + line of sight
  -> Moon body-fixed frame
  -> intersection with Moon ellipsoid
  -> selenographic lon/lat
```

これは月面上の同じ経度・緯度へ複数日の観測を重ねたい場合に使います。ただし、ポインティングずれに敏感で、特に月縁付近では小さな sky offset が大きな lon/lat 誤差になります。

## 5. 推奨される使い分け

### 月を約 30 分角の円盤として表示し、月の向きを固定したい

```python
coord_sys = "moon_disk_offset"
projection = "SFL"
```

これが通常の推奨です。

### 月の生の dAz/dEl 観測パターンを確認したい

```python
coord_sys = "moon_azel_offset"
projection = "SFL"
```

これは scan pattern、turnaround、coverage の確認に向いています。

### 月中心の天球北基準 offset を見たい

```python
coord_sys = "moon_sky_offset"
projection = "SFL"
```

AltAz 視野回転は消えますが、月北ではなく天球北に固定されます。

### 太陽 map を移動天体中心で表示したい

観測確認なら、

```python
coord_sys = "sun_azel_offset"
projection = "SFL"
```

天球北基準の見かけ円盤なら、

```python
coord_sys = "sun_sky_offset"
projection = "SFL"
```

です。今回の実装では、太陽の自転軸方向へ回転する `sun_disk_offset` は未実装です。

### 月面経度・緯度に載せたい

```python
coord_sys = "moon_selenographic"
projection = "CAR"
```

SPICE kernel と `spiceypy` が必要です。

## 6. 基本的な使用例

### 6.1 単一 scantable の月円盤 map

```python
run_otf_full_pipeline(
    scantable,
    config,
    output_fits="moon_disk_offset.fits",
    coord_sys="moon_disk_offset",
    projection="SFL",
)
```

### 6.2 複数 scantable を合成する

複数 map は、個別 FITS を後で足すより、scantable 段階でまとめて処理することを推奨します。

```python
run_otf_full_pipeline_multi(
    [scantable_0, scantable_1, scantable_2],
    config,
    output_fits="moon_disk_offset_merged.fits",
    coord_sys="moon_disk_offset",
    projection="SFL",
    do_basket_weave=True,
)
```

理由は、各 dump の時刻で Moon/Sun 中心を引き、同じ座標系に載せたうえで gridding と basket weave を行えるためです。

### 6.3 太陽の dAz/dEl map

```python
run_otf_full_pipeline(
    scantable,
    config,
    output_fits="sun_azel_offset.fits",
    coord_sys="sun_azel_offset",
    projection="SFL",
)
```

### 6.4 太陽の天球接平面 offset map

```python
run_otf_full_pipeline(
    scantable,
    config,
    output_fits="sun_sky_offset.fits",
    coord_sys="sun_sky_offset",
    projection="SFL",
)
```

## 7. pointing 補正

解析時に `pointing_dx_arcsec`, `pointing_dy_arcsec` を与えることができます。ファイルを別に用意する必要はありません。

### 7.1 符号規約

補正は次の意味です。

$$
(x_{\rm corrected}, y_{\rm corrected})
=
(x_{\rm original}, y_{\rm original}) + (dx, dy)
$$

例えば、`moon_sky_offset` や `moon_disk_offset` の map で月中心が

```text
x_fit = +12 arcsec
y_fit = -8 arcsec
```

に見えていた場合、中心を $(0,0)$ に戻す補正は

```python
pointing_dx_arcsec = -12.0
pointing_dy_arcsec = +8.0
```

です。

### 7.2 単一 scantable の補正

```python
run_otf_full_pipeline(
    scantable,
    config,
    output_fits="moon_disk_offset_corrected.fits",
    coord_sys="moon_disk_offset",
    projection="SFL",
    pointing_dx_arcsec=-12.0,
    pointing_dy_arcsec=+8.0,
)
```

### 7.3 複数 scantable で補正量が異なる場合

list で渡せます。

```python
run_otf_full_pipeline_multi(
    [scantable_0, scantable_1, scantable_2],
    config,
    output_fits="moon_disk_offset_merged_corrected.fits",
    coord_sys="moon_disk_offset",
    projection="SFL",
    pointing_dx_arcsec=[-12.0, -4.0, +10.0],
    pointing_dy_arcsec=[ +8.0, +3.0,  -6.0],
    do_basket_weave=True,
)
```

または dict/list 形式の `pointing_correction` も使えます。

```python
pointing_correction = [
    {"dx_arcsec": -12.0, "dy_arcsec": +8.0},
    {"dx_arcsec":  -4.0, "dy_arcsec": +3.0},
    {"dx_arcsec": +10.0, "dy_arcsec": -6.0},
]

run_otf_full_pipeline_multi(
    [scantable_0, scantable_1, scantable_2],
    config,
    output_fits="moon_disk_offset_merged_corrected.fits",
    coord_sys="moon_disk_offset",
    projection="SFL",
    pointing_correction=pointing_correction,
    do_basket_weave=True,
)
```

空の scantable が途中にあって skip される場合でも、補正量は元の入力 index に対応するように処理されます。

### 7.4 補正が適用される座標系

補正は、次の座標系で意味を持ちます。

- `moon_sky_offset`
- `sun_sky_offset`
- `moon_disk_offset`
- `moon_selenographic`

`moon_disk_offset` では、補正はまず `moon_sky_offset` 軸で加え、その後で月北方向へ回転します。したがって、補正値は sky offset frame の `dRA cos Dec`, `dDec` として解釈します。

`moon_selenographic` では、補正後の視線方向で SPICE の月面交点を計算します。gridding 後の画像を平行移動するのではありません。

`moon_azel_offset` / `sun_azel_offset` では、pointing 補正は記録されますが、座標には適用されません。これは補正量が sky offset frame で定義されており、Az/El offset へ単純に足すと意味が変わるためです。

## 8. basket weave との関係

basket weave は最終的な `x, y` 座標と scan ID を用いて baseline offset を推定します。そのため、`moon_disk_offset` や `sun_sky_offset` でも原理的には使用できます。

推奨順序は次です。

1. 各 row の座標を Moon/Sun offset frame へ変換する。
2. pointing 補正を適用する。
3. 補正後座標で basket weave の交差・近接関係を評価する。
4. baseline 補正後に gridding する。

basket weave は baseline を補正する処理であり、pointing error 自体は補正しません。複数 scan で pointing error が異なる場合は、basket weave の前に scantable ごとの pointing 補正を入れてください。

## 9. FITS 出力の方針

今回の移動天体 offset map は、固定天球上の RA/Dec image ではありません。そのため、`RA---TAN` / `DEC--TAN` に偽装しません。

FITS としては、以下のような **linear WCS** として保存します。

```text
CTYPE1  = 'OFFSETX'
CTYPE2  = 'OFFSETY'
CUNIT1  = 'deg'
CUNIT2  = 'deg'
```

各座標系の意味は追加 header keyword で明示します。

### 9.1 `moon_azel_offset` / `sun_azel_offset`

```text
COORDSYS = 'AZEL_OFFSET'
OFFSYS   = 'AZEL'
OFFBODY  = 'MOON' or 'SUN'
OFFXDEF  = 'dAz*cos(El_body)'
OFFYDEF  = 'dEl'
```

### 9.2 `moon_sky_offset` / `sun_sky_offset`

```text
COORDSYS = 'SKY_OFFSET'
OFFSYS   = 'RADEC'
OFFBODY  = 'MOON' or 'SUN'
OFFXDEF  = 'dRA*cos(Dec_body)'
OFFYDEF  = 'dDec'
```

### 9.3 `moon_disk_offset`

```text
COORDSYS = 'MOON_DISK_OFFSET'
OFFSYS   = 'MOON_DISK'
OFFBODY  = 'MOON'
DISKREF  = 'APPARENT'
DISKYPOS = 'LUNAR_NORTH'
DISKXPOS = 'PA_NORTH_PLUS_90'
```

この FITS は、他の FITS viewer で通常の image/cube として表示できます。ただし celestial WCS ではないため、RA/Dec catalog overlay や celestial north arrow は使いません。

### 9.4 `moon_selenographic`

`projection="CAR"` では、

```text
CTYPE1  = 'SELON'
CTYPE2  = 'SELAT'
COORDSYS= 'SELENOGRAPHIC'
BODY    = 'MOON'
SELONDIR= 'EAST'
SELLAT  = 'PLANETOCENTRIC'
SURFACE = 'ELLIPSOID'
```

です。

`projection="SFL"` では、surface longitude/latitude を線形 offset として gridding するため、

```text
CTYPE1  = 'SELOFFX'
CTYPE2  = 'SELOFFY'
OFFXDEF = '(SELON-SELON0)*cos(SELAT)'
OFFYDEF = 'SELAT-SELAT0'
```

になります。

## 10. plotting.py の表示

`plotting.py` は今回の非天球 WCS を認識します。

軸ラベルは概ね以下になります。

| 座標系 | X 軸 | Y 軸 |
|---|---|---|
| `moon_azel_offset`, `sun_azel_offset` | `Az/El offset X = dAz cos(El_body) [deg]` | `Az/El offset Y = dEl [deg]` |
| `moon_sky_offset`, `sun_sky_offset` | `Moon/Sun sky offset X = dRA cos(Dec_body) [deg]` | `Moon/Sun sky offset Y = dDec [deg]` |
| `moon_disk_offset` | `Moon disk X: +90 deg eastward from lunar north [deg]` | `Moon disk Y: apparent lunar north [deg]` |
| `moon_selenographic` | `Selenographic Longitude [deg]` | `Selenographic Latitude [deg]` |

`moon_disk_offset` は非天球 WCS なので、`north_arrow` と `SkyCoord` overlay は自動的に抑制またはエラーになります。点や注釈を重ねる場合は、pixel/xy 座標を使ってください。

`spectral-cube` は celestial WCS を要求する場合があるため、`OFFSETX/OFFSETY` や `SELON/SELAT` の cube を直接読めないことがあります。`plotting.py` では plain FITS fallback で 2D 表示・moment 表示できるようにしています。

## 11. 依存関係

### 11.1 Astropy

以下の座標系では Astropy が必要です。

- `moon_azel_offset`
- `sun_azel_offset`
- `moon_sky_offset`
- `sun_sky_offset`
- `moon_disk_offset`

Astropy は、各 dump 時刻・観測地に対する Moon/Sun の apparent position、AltAz 変換、sky offset 計算に使います。

### 11.2 SpiceyPy

`spiceypy` は `moon_selenographic` を使う時だけ import します。通常の Moon/Sun offset map や `moon_disk_offset` では import しません。

`moon_selenographic` には SPICE kernel が必要です。代表例は次です。

```text
naif0012.tls
pck00011.tpc
de440s.bsp or de440.bsp
moon_pa*.bpc
moon_de*.tf
```

実際に使う kernel は、必要精度と運用環境に合わせて指定してください。

## 12. 推奨ワークフロー

### 12.1 月 map の通常運用

1. `moon_azel_offset` で scan pattern と coverage を確認する。
2. `moon_disk_offset` で月円盤を表示する。
3. 月縁または円盤中心を fit して pointing offset を推定する。
4. `pointing_dx_arcsec`, `pointing_dy_arcsec` を指定して再解析する。
5. 複数 scantable を合成する場合は、scantable ごとの補正量を list で与える。
6. 月面経度・緯度で比較したい場合だけ `moon_selenographic` を使う。

### 12.2 太陽 map の通常運用

1. `sun_azel_offset` で scan pattern と coverage を確認する。
2. `sun_sky_offset` で太陽中心の天球接平面 map を作る。
3. 必要なら `pointing_dx_arcsec`, `pointing_dy_arcsec` を指定して補正する。
4. 複数 scantable は `run_otf_full_pipeline_multi()` でまとめて gridding する。

## 13. 注意事項

- `moon_disk_offset` は月面経度・緯度ではありません。
- `moon_selenographic` は pointing error に敏感です。
- 月縁付近では `moon_selenographic` の lon/lat が急激に変化します。
- `BORE_AZ/BORE_EL` はマルチビームで危険なため座標計算には使いません。
- 移動天体 offset map は celestial WCS ではないため、RA/Dec overlay は使わないでください。
- 個別 FITS を後で足すより、scantable 段階で複数入力をまとめる方が安全です。

## 14. 最小確認例

### 月円盤を月北固定で表示

```python
run_otf_full_pipeline_multi(
    scantables,
    config,
    output_fits="moon_disk_merged.fits",
    coord_sys="moon_disk_offset",
    projection="SFL",
    do_basket_weave=True,
)
```

### scantable ごとに pointing 補正して表示

```python
run_otf_full_pipeline_multi(
    scantables,
    config,
    output_fits="moon_disk_merged_corrected.fits",
    coord_sys="moon_disk_offset",
    projection="SFL",
    pointing_dx_arcsec=[-12.0, -4.0, +10.0],
    pointing_dy_arcsec=[ +8.0, +3.0,  -6.0],
    do_basket_weave=True,
)
```

### 太陽中心 sky offset map

```python
run_otf_full_pipeline(
    scantable,
    config,
    output_fits="sun_sky_offset.fits",
    coord_sys="sun_sky_offset",
    projection="SFL",
)
```

### SPICE 月面固定 map

```python
run_otf_full_pipeline(
    scantable,
    config,
    output_fits="moon_selenographic.fits",
    coord_sys="moon_selenographic",
    projection="CAR",
    pointing_dx_arcsec=-12.0,
    pointing_dy_arcsec=+8.0,
)
```

`moon_selenographic` を本格運用する前に、`moon_disk_offset` で pointing offset を確認することを推奨します。
