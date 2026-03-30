# 単一鏡電波分光解析パッケージ Cookbook
## 実用例・推奨レシピ集（日本語, 2026-03-29 アーカイブ対応版）

> 方針  
> 旧版 cookbook の情報を落とさずに、最新版アーカイブの公開 API と内部意味論に合わせて更新した版です。  
> 例は基本的に **単発で意味が読める形** を保ち、必要なところだけ補足を増やしています。  
> 本書の対象は主に `fitsio.py`, `rawspec.py`, `calibrate.py`, `baseline.py`, `regrid_vlsrk.py`, `coadd.py`, `restfreq.py`, `tempscale.py`, `scantable_utils.py`, `axis.py`, `doppler.py`, `atmosphere.py`, `sdfits_writer.py` です。  
> `profile_view` 系の詳細仕様は、最新版の `profile_view_manual_ja_v2.md` を正本として参照してください。本書では重複を避け、解析後確認に使う簡単な例だけを残します。

---

## 0. まず最初に押さえる考え方

このパッケージの解析は、概ね次の流れです。

1. raw を読む  
2. `Ta*` キャリブレーションを行う  
3. 必要なら `rest_freq` を解釈上だけ切り替える  
4. baseline を引く、または baseline 窓で RMS を測る  
5. `LSRK` に揃えた共通速度グリッドへ再配置する  
6. scan 単位、position 単位、`INTGRP` 単位で coadd する  
7. 必要なら保存時に `TA*` / `TR*` を変換して書く

特に重要なのは、次の量の意味を混同しないことです。

- `DATA`: 実際のスペクトル配列
- `TEMPSCAL`: その配列が `TA*` か `TR*` かというラベル
- `BEAMEFF`: `TA*` と `TR*` を往復するときの効率
- `SPECSYS`: 今その軸が属している基準系
- `SSYSOBS`: 観測時に固定だった基準系。`regrid` / `coadd` / `baseline` を通しても、通常はこの意味で保持される
- `VELOSYS` / `VFRAME`: 未適用の行ごとの速度補正列
- `RESTFRQ` / `RESTFREQ`: 速度解釈に使う静止周波数

coadd の重み付き平均は、概念的には

$$
\bar T_k = \frac{\sum_i w_i T_{i,k}}{\sum_i w_i}
$$

であり、RMS 重みでは通常

$$
w_i = \frac{1}{\sigma_i^2}
$$

です。

---

## 1. 事前準備

多くの例では、まず次を入れておくと扱いやすいです。

```python
import astropy.utils.iers
astropy.utils.iers.conf.auto_download = False
astropy.utils.iers.conf.auto_max_age = None
astropy.utils.iers.conf.iers_degraded_accuracy = "ignore"

import sd_radio_spectral_fits.fitsio as fitsio
import sd_radio_spectral_fits.rawspec as rawspec
import sd_radio_spectral_fits.calibrate as cal
import sd_radio_spectral_fits.baseline as bsl
import sd_radio_spectral_fits.coadd as coadd
import sd_radio_spectral_fits.regrid_vlsrk as rv
import sd_radio_spectral_fits.restfreq as rf
import sd_radio_spectral_fits.tempscale as ts
import sd_radio_spectral_fits.scantable_utils as su
import sd_radio_spectral_fits.axis as axis
import sd_radio_spectral_fits.doppler as dop
import sd_radio_spectral_fits.atmosphere as atm
```

旧版 cookbook にあった表示系も、環境にあれば次のように読みます。

```python
# profile_view の詳細は別冊 profile_view_manual_ja_v2.md を参照
import sd_radio_spectral_fits.profile_view as pv
```

---

## 2. この cookbook で使う最小限の用語

### 2.1 Scantable

`fitsio.Scantable` は次の 4 つをまとめたコンテナです。

- `meta`: ファイル全体に近いメタデータ
- `data`: スペクトル配列。2 次元配列でも、行ごと長さの違う `list[np.ndarray]` でもよい
- `table`: 行ごとの補助表
- `history`: 処理履歴

### 2.2 RawSpec

`rawspec.RawSpec` は raw の `HOT`, `OFF`, `ON` を持つコンテナです。

- `hot`
- `off`
- `on`
- `mapping`
- `meta`

### 2.3 行ごとの WCS と file-level WCS

このパッケージでは、`meta` にある `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`, `RESTFRQ` などを、必要に応じて `table` の列にも昇格させて扱います。  
したがって、解析の途中では「file-level にしか無い」つもりの値が table 側にも見えることがありますが、これは仕様です。

### 2.4 `TA*` と `TR*`

定義は

$$
T_R^* = \frac{T_A^*}{\eta_{\rm B}}
$$

$$
T_A^* = T_R^* \eta_{\rm B}
$$

です。ここで $\eta_{\rm B}$ が `BEAMEFF` です。

最新版では **read 時に勝手に `TR*` を `TA*` に戻しません**。  
`TEMPSCAL` を尊重し、必要なら coadd 出力時や保存時に変換します。

---

## 3. 最小の 12CO 解析

### 目的
1 本の raw FITS から `Ta*` を作り、scan coadd、baseline、最終 coadd まで行う。

```python
raw_fits = "orikl_raw_12co_20260222_111212.fits"

vwin_baseline = ["-35:-5", "30:55"]
vwin_rms = "-200:0"
npoly = 7

# 1. 読み込み
sc_raw = fitsio.read_scantable(raw_fits)

# 2. Ta* キャリブレーション
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-200, 200),
    t_hot_k=300.0,
    vcorr_chunk_sec=5,
    dtype="float32",
    gain_mode="hybrid",
)

# 3. scan 単位で coadd
sc_scan = coadd.run_velocity_coadd(
    inputs=sc_cal,
    group_mode="scan",
    mode="rms_weight",
    rms_vwin=vwin_rms,
    rms_bin=21,
    vmin=-150,
    vmax=150,
    pos_tol_arcsec=1.0,
    overwrite=True,
)

# 4. baseline
sc_scan_bsl = bsl.run_baseline_fit(
    input_data=sc_scan,
    poly_order=npoly,
    vwin=vwin_baseline,
    overwrite=True,
)

# 5. 同一位置を最終 coadd
sc_final = coadd.run_velocity_coadd(
    inputs=sc_scan_bsl,
    mode="rms_weight",
    pos_tol_arcsec=1.0,
    overwrite=True,
)

# 6. 必要なら final baseline
sc_final_bsl = bsl.run_baseline_fit(
    input_data=sc_final,
    poly_order=npoly,
    vwin=vwin_baseline,
    overwrite=True,
)
```

### 重要: マルチビームの相対スケール合わせはここで行う
beam 間の相対スケール係数があるなら、**calibration 直後**に `apply_relative_scale()` を掛けてから baseline / coadd へ進むのが安全です。

```python
rel_map = {
    (0, 0, 0): 1.00,
    (0, 0, 1): 0.97,
    (1, 0, 0): 1.03,
    (1, 0, 1): 1.01,
}

su.apply_relative_scale(
    sc_cal,
    rel_map,
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
)
```

後で全体の絶対係数が決まったら、さらに

```python
su.apply_global_scale(sc_cal, 0.91)
```

のように一律 factor を掛けます。これらは data 本体へ直接掛かり、`INTSCALE` に累積係数が記録され、既知の scale 依存派生列は既定で無効化されます。

### ポイント
- `group_mode="scan"` により、まず scan 内の dump をまとめる
- そのあと baseline
- 最後に同位置の repeat をまとめる
- `gain_mode="hybrid"` が現行の推奨値

### 3.1 multi-beam で calibration 直後に相対補正を入れる完全版

相対スケール係数 $g_i$ と後日の一律係数 $f$ は、どちらも **data へ掛ける係数**です。

$$
T_i^{\mathrm{aligned}}(k) = g_i\,T_i(k)
$$

$$
T_i^{\mathrm{final}}(k) = f\,T_i^{\mathrm{aligned}}(k)
$$

一方、物理 `BEAMEFF` は

$$
T_R^* = \frac{T_A^*}{\eta}
$$

の $\eta$ なので、beam 間相対補正の係数に流用しないでください。

```python
raw_fits = "multibeam_raw.fits"

rel_map = {
    (0, 0, 0): 1.00,
    (0, 0, 1): 0.97,
    (1, 0, 0): 1.03,
    (1, 0, 1): 1.01,
}

sc_raw = fitsio.read_scantable(raw_fits)
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-200, 200),
    t_hot_k=300.0,
    vcorr_chunk_sec=5,
    dtype="float32",
)

# 1. まず calibration 直後に beam 間相対補正を掛ける
su.apply_relative_scale(
    sc_cal,
    rel_map,
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
)

# 2a. generic な一律係数が後で決まったなら data へ掛ける
# su.apply_global_scale(sc_cal, 0.91)

# 2b. 後で決まった係数が物理 TR* 変換の eta なら、こちらでもよい
# su.set_beameff(sc_cal, 0.42)
# ただし set_beameff() 自体は data を掛け直さず、後段の TR* 変換に使われる metadata を入れるだけ

sc_scan = coadd.run_velocity_coadd(
    inputs=sc_cal,
    group_mode="scan",
    mode="rms_weight",
    rms_vwin="-200:0",
    pos_tol_arcsec=1.0,
)

sc_scan_bsl = bsl.run_baseline_fit(
    input_data=sc_scan,
    poly_order=3,
    vwin=["-35:-5", "30:55"],
)

sc_final = coadd.run_velocity_coadd(
    inputs=sc_scan_bsl,
    mode="rms_weight",
    pos_tol_arcsec=1.0,
)
```

### この例のポイント
- 相対補正は **calibration 直後** に 1 回だけ掛ける
- そうすると `BSL_RMS`, moment0, coadd が最初から揃った強度スケールで計算される
- `apply_global_scale()` は generic な一律係数用
- `set_beameff()` は物理 `TR*` 変換用 metadata 用

---

## 4. Ta* キャリブレーションの式を明示して理解する

最新版の `run_tastar_calibration()` / `make_tastar_dumps()` では、概念的にはまず

$$
T_A^*(t) = \left(P_{\rm ON}(t) - \widetilde{P}_{\rm OFF}(t)\right)
\frac{T_{\rm cal}(t)}{D(t)}
$$

を作ります。

ここで

- $P_{\rm ON}(t)$ は ON スペクトル
- $\widetilde{P}_{\rm OFF}(t)$ は ON 時刻へ内挿した OFF
- $T_{\rm cal}(t)$ は chopper-wheel の等価温度
- $D(t)$ は分母のゲイン配列

です。

### 4.1 `independent`

`gain_mode="independent"` では

$$
D(t) = \widetilde{P}_{\rm HOT}(t) - \widetilde{P}_{\rm OFF}(t)
$$

です。

### 4.2 `hybrid`

`gain_mode="hybrid"` では、まず HOT 時刻で OFF を引いたクリーンな分母を作ります。

$$
D_{\rm hot}(t_{\rm hot}) = P_{\rm HOT}(t_{\rm hot}) - \widetilde{P}_{\rm OFF}(t_{\rm hot})
$$

その後、それを ON 時刻へ内挿して

$$
D(t_{\rm on}) = \widetilde{D}_{\rm hot}(t_{\rm on})
$$

を使います。  
実装上は HOT/OFF も scan ごとに平均してから時間補間するので、旧版より分母が安定です。

---

## 5. 複数ファイルの 12CO を統合する

### 目的
複数日にまたがる、または複数ファイルに分かれた 12CO を最終統合する。

```python
file_list = [
    "orikl_raw_12co_20260222_111212.fits",
    "orikl_raw_12co_20260222_112005.fits",
    "orikl_raw_12co_20260222_113301.fits",
]

vwin_baseline = ["-35:-5", "30:55"]
vwin_rms = "-200:0"
npoly = 7

sc1 = cal.run_tastar_calibration(file_list[0], vlsrk_range_kms=(-200, 200), t_hot_k=300, vcorr_chunk_sec=5, dtype="float32")
sc1 = coadd.run_velocity_coadd(sc1, group_mode="scan", mode="rms_weight", rms_vwin=vwin_rms, rms_bin=21, vmin=-150, vmax=150, pos_tol_arcsec=1.0)
sc1 = bsl.run_baseline_fit(sc1, poly_order=npoly, vwin=vwin_baseline)
sc1 = coadd.run_velocity_coadd(sc1, mode="rms_weight", pos_tol_arcsec=1.0)
sc1 = bsl.run_baseline_fit(sc1, poly_order=npoly, vwin=vwin_baseline)

sc2 = cal.run_tastar_calibration(file_list[1], vlsrk_range_kms=(-200, 200), t_hot_k=300, vcorr_chunk_sec=5, dtype="float32")
sc2 = coadd.run_velocity_coadd(sc2, group_mode="scan", mode="rms_weight", rms_vwin=vwin_rms, rms_bin=21, vmin=-150, vmax=150, pos_tol_arcsec=1.0)
sc2 = bsl.run_baseline_fit(sc2, poly_order=npoly, vwin=vwin_baseline)
sc2 = coadd.run_velocity_coadd(sc2, mode="rms_weight", pos_tol_arcsec=1.0)
sc2 = bsl.run_baseline_fit(sc2, poly_order=npoly, vwin=vwin_baseline)

sc3 = cal.run_tastar_calibration(file_list[2], vlsrk_range_kms=(-200, 200), t_hot_k=300, vcorr_chunk_sec=5, dtype="float32")
sc3 = coadd.run_velocity_coadd(sc3, group_mode="scan", mode="rms_weight", rms_vwin=vwin_rms, rms_bin=21, vmin=-150, vmax=150, pos_tol_arcsec=1.0)
sc3 = bsl.run_baseline_fit(sc3, poly_order=npoly, vwin=vwin_baseline)
sc3 = coadd.run_velocity_coadd(sc3, mode="rms_weight", pos_tol_arcsec=1.0)
sc3 = bsl.run_baseline_fit(sc3, poly_order=npoly, vwin=vwin_baseline)

sc_total = coadd.run_velocity_coadd(
    inputs=[sc1, sc2, sc3],
    mode="rms_weight",
    pos_tol_arcsec=1.0,
    overwrite=True,
)
```

### ポイント
- cookbook では loop を外して書いた
- 実運用では list や loop に戻してよい
- `inputs=[sc1, sc2, sc3]` としてそのまま最終統合できる

---

## 6. 13CO / C18O を同じ raw から分けて解析する

### 目的
同一 raw データを、異なる `rest_freq` で解釈して別線として処理する。

```python
raw_fits = "orikl_raw_13co_20260222_111212.fits"

f_13co = 110.2014e9
f_c18o = 109.7821734e9

vwin_baseline_13 = ["-10:5", "15:30"]
vwin_baseline_c18o = ["-10:2.5", "15:27"]
vwin_rms = "-200:0"

sc_raw = fitsio.read_scantable(raw_fits)

sc_13 = cal.run_tastar_calibration(
    input_data=sc_raw,
    rest_freq=f_13co,
    vlsrk_range_kms=(-200, 200),
    t_hot_k=300.0,
    vcorr_chunk_sec=5,
    dtype="float32",
)
sc_13 = coadd.run_velocity_coadd(
    inputs=sc_13,
    group_mode="scan",
    mode="rms_weight",
    rms_vwin=vwin_rms,
    rms_bin=21,
    vmin=-150,
    vmax=150,
    pos_tol_arcsec=1.0,
)
sc_13 = bsl.run_baseline_fit(sc_13, poly_order=7, vwin=vwin_baseline_13)
sc_13 = coadd.run_velocity_coadd(sc_13, mode="rms_weight", pos_tol_arcsec=1.0)
sc_13 = bsl.run_baseline_fit(sc_13, poly_order=7, vwin=vwin_baseline_13)

sc_c18o = cal.run_tastar_calibration(
    input_data=sc_raw,
    rest_freq=f_c18o,
    vlsrk_range_kms=(-200, 200),
    t_hot_k=300.0,
    vcorr_chunk_sec=5,
    dtype="float32",
)
sc_c18o = coadd.run_velocity_coadd(
    inputs=sc_c18o,
    group_mode="scan",
    mode="rms_weight",
    rms_vwin=vwin_rms,
    rms_bin=41,
    vmin=-150,
    vmax=150,
    pos_tol_arcsec=1.0,
)
sc_c18o = bsl.run_baseline_fit(sc_c18o, poly_order=7, vwin=vwin_baseline_c18o)
sc_c18o = coadd.run_velocity_coadd(sc_c18o, mode="rms_weight", pos_tol_arcsec=1.0)
sc_c18o = bsl.run_baseline_fit(sc_c18o, poly_order=7, vwin=vwin_baseline_c18o)
```

### ポイント
- `CTYPE1="FREQ"` のときは、`rest_freq` を変えても元の周波数 WCS をむやみに壊さない
- 13CO と C18O で baseline 窓や `rms_bin` を変えられる
- `rest_freq` は calibration と baseline / coadd のどちらにも渡せるが、流れの最初で明示しておく方が混乱しにくい

---

## 7. ATM モデルを使った Ta* キャリブレーション

### 目的
`tau_zenith` と外気温から `T_cal` を厳密化する。

大気モデルではまずエアマス

$$
X = \sec z = \frac{1}{\sin {\rm EL}}
$$

を使います。実装では低仰角暴走を防ぐため、`EL < 5 度` では `5 度` 相当にクリップします。

`T_atm` のモデルは 2 つあります。

$$
T_{\rm atm} = T_{\rm surf} - \Delta T
$$

または

$$
T_{\rm atm} = \eta T_{\rm surf}
$$

です。

`T_cal` は

$$
T_{\rm cal} = T_{\rm hot} e^{\tau_0 X} - T_{\rm atm}(e^{\tau_0 X}-1) - T_{\rm bg}
$$

で計算されます。

```python
sc_atm = cal.run_tastar_calibration(
    input_data="orikl_raw_12co_20260222_111212.fits",
    t_hot_k=300.0,
    tau_zenith="auto",      # header / table の TAU0 を使う
    t_atm_model="offset",
    t_atm_delta_k=15.0,
    vlsrk_range_kms=(-200, 200),
    vcorr_chunk_sec=5,
    dtype="float32",
    verbose=True,
)
```

### `ratio` モデルを使う

```python
sc_atm_ratio = cal.run_tastar_calibration(
    input_data="orikl_raw_12co_20260222_111212.fits",
    t_hot_k=300.0,
    tau_zenith=0.12,
    t_surface_k=278.0,
    t_atm_model="ratio",
    t_atm_eta=0.95,
    vlsrk_range_kms=(-200, 200),
    vcorr_chunk_sec=5,
    dtype="float32",
)
```

### ポイント
- `tau_zenith=None` なら Basic 1-Temp
- `tau_zenith="auto"` は `TAU0` がないと止まる
- `t_surface_k` 未指定なら `WXTEMP`, `TAMBIENT`, `TAMB`, `T_SURF` を探す
- `TCAL`, `THOT`, `TAU0`, `CALSTAT` は table 側に残る

---

## 8. 既存 Ta* を再キャリブレーションする

### 目的
RAW を読み直さずに ATM パラメータだけ更新する。

```python
sc_recal = cal.recalibrate_tastar(
    scantable=sc_atm,
    new_tau=0.10,
    new_t_surface_k=275.0,
    new_t_atm_model="offset",
    new_t_atm_delta_k=12.0,
    verbose=True,
)
```

### Basic 1-Temp に戻す

```python
sc_basic = cal.recalibrate_tastar(
    scantable=sc_atm,
    new_tau=None,
    verbose=True,
)
```

### ポイント
- `THOT`, `TCAL` が入っていれば高速に巻き戻せる
- 旧データで `THOT` が無い場合は `TAMB_K` にフォールバックする安全策がある
- 過去データで摂氏が混じる可能性にも保護が入っている

---

## 9. baseline のみ単独でかける

### 目的
すでに coadd 済みのデータへ baseline だけかける。

```python
sc_bsl = bsl.run_baseline_fit(
    input_data="coadded_before_baseline.fits",
    vwin=["-35:-5", "30:55"],
    poly_order=7,
    line_vwin=["5:20"],
    iter_max=3,
    iter_sigma=3.0,
    overwrite=True,
)
```

### ポイント
- `line_vwin` は baseline 窓から差し引かれる
- `iter_max > 0` にすると sigma clipping を回せる
- baseline fit は **入力が `TA*` でも `TR*` でも、そのままのスケールで行う**

---

## 10. baseline を評価だけして、データは変えない

### 目的
`BSL_RMS` や係数だけ欲しいときに、スペクトルは変えずに走らせる。

```python
sc_bsl_eval = bsl.run_baseline_fit(
    input_data="coadded_before_baseline.fits",
    vwin=["-35:-5", "30:55"],
    poly_order=5,
    apply=False,
    overwrite=True,
)
```

### 出力される代表列

- `BSL_DONE`
- `BSL_APPLIED`
- `BSL_STAGE`
- `BSL_POLY`
- `BSL_WINF`
- `BSL_RMS`
- `BSL_STAT`
- `BSL_NUSED`
- `BSL_COEF`
- `BSL_SCALE`

### ポイント
- `apply=False` でも `BSL_*` は保存される
- 後段の `coadd(mode="rms")` が `BSL_RMS` をそのまま使える

---

## 11. baseline で row 抽出、channel slice、失敗時方針を使う

```python
sc_sub_bsl = bsl.run_baseline_fit(
    input_data="coadded_before_baseline.fits",
    rows="0:20,50,60:80",
    vwin=["-35:-5", "30:55"],
    poly_order=5,
    ch_start=100,
    ch_stop=3500,
    on_fail="warn",
    overwrite=True,
)
```

### `exclude_rows` を使う

```python
sc_sub_bsl = bsl.run_baseline_fit(
    input_data="coadded_before_baseline.fits",
    exclude_rows="0:100",
    vwin=["-35:-5", "30:55"],
    poly_order=5,
    on_fail="skip",
)
```

### `bsl_overwrite="error"`

```python
sc_sub_bsl = bsl.run_baseline_fit(
    input_data="already_baselined.fits",
    vwin=["-35:-5", "30:55"],
    poly_order=3,
    bsl_overwrite="error",
)
```

### ポイント
- `rows` と `exclude_rows` は同時指定不可
- `on_fail` は `exit`, `warn`, `skip` の考え方で使う
- `bsl_overwrite="error"` は既存 `BSL_*` の accidental overwrite 防止に有効

---

## 12. coadd を「scan 単位」だけで使う

### 目的
まず raw dump を短時間平均したい。

```python
sc_scan = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="scan",
    mode="rms_weight",
    rms_vwin="-200:0",
    rms_poly=1,
    rms_bin=21,
    vmin=-150,
    vmax=150,
    axis_type="freq",
)
```

### ポイント
- scan ごとに 1 本のスペクトルへまとめる
- 後続の baseline が安定しやすい
- `axis_type="freq"` なら出力 WCS は `FREQ/Hz`

---

## 13. `baseline_vwin` を使う重み付き coadd

### 目的
個々のスペクトルで baseline を引き、その残差 RMS で重み付けしたい。

```python
sc = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="position",
    mode="rms_weight",
    baseline_vwin=["-35:-5", "30:55"],
    baseline_poly=3,
    baseline_iter_max=2,
    baseline_iter_sigma=3.0,
    line_vwin=["5:20"],
    axis_type="vel",
)
```

### ポイント
- `baseline_vwin` は、重み評価前に baseline subtraction を行う経路
- `line_vwin` は baseline 窓からさらに差し引かれる
- `axis_type="vel"` なら出力 WCS は `VRAD/m/s`

---

## 14. `rms_vwin` を使う重み付き coadd

### 目的
データ本体の baseline を引かず、ノイズ評価だけしたい。

```python
sc = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="position",
    mode="rms_weight",
    rms_vwin=["-200:0"],
    rms_poly=1,
    rms_bin=21,
    axis_type="vel",
)
```

### ポイント
- baseline をデータへ適用しない
- まず重みだけ決めたい場合に向く
- `baseline_vwin` と `rms_vwin` は同時指定不可

---

## 15. `BSL_RMS` を使って coadd する

### 目的
入力がすでに baseline 済みで、各行の `BSL_RMS` をそのまま使いたい。

```python
sc = coadd.run_velocity_coadd(
    inputs="baseline_done.fits",
    mode="rms_weight",
    group_mode="position",
    axis_type="vel",
)
```

### 前提
- `baseline_vwin` も `rms_vwin` も与えない
- 入力 table に `BSL_RMS` が入っている

### ポイント
- この場合、重みの情報源は input table です
- 出力 row には `WGT_SRC="input_bsl_rms"` などの情報が残る

---

## 16. QC モードを使う

### 目的
怪しいスペクトルを robust に落としながら coadd する。

```python
sc_qc = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="position",
    mode="rms_weight",
    rms_vwin=["-200:0"],
    coadd_qc="robust",
    axis_type="vel",
    verbose=True,
)
```

### パラメータ付き指定

```python
sc_qc = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="position",
    mode="rms_weight",
    rms_vwin=["-200:0"],
    coadd_qc="robust;clip=3,2;hard=4;soft=1,4",
)
```

### 代表 preset

- `robust`
- `gentle`
- `hardonly`
- `noclip`

### ポイント
- `coadd_qc` を使うと `mode="uniform"` は不可
- QC は統計的重み付けを内在するので、`rms` 系モードと組み合わせる

---

## 17. `weight_zero_policy` を使い分ける

### 目的
重み計算で RMS が壊れた行があるときの扱いを制御する。

```python
sc = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="position",
    mode="rms_weight",
    rms_vwin=["-200:0"],
    weight_zero_policy="drop",
)
```

### 許される値

- `error`
- `drop`
- `impute_median`

### ポイント
- `error` は最も保守的
- `drop` は壊れた行だけ捨てる
- `impute_median` は典型重みで置換する

---

## 18. post-baseline を coadd 後に入れる

### 目的
coadd の後にもう一度 baseline を当てたい。

```python
sc = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="position",
    mode="rms_weight",
    baseline_vwin=["-35:-5", "30:55"],
    baseline_poly=3,
    post_baseline_mode="inherit_all",
)
```

ここで `inherit_all` は

- `post_baseline_vwin`
- `post_baseline_poly`
- `post_baseline_iter_max`
- `post_baseline_iter_sigma`

を一括で継承します。したがって、同時に各項目へ `"inherit"` を書く必要はありません。API の既定値として内部的には `inherit` が見えても、利用者側でそれを毎回書く必要はありません。例では省略した方が誤解がありません。

### 一部だけ上書きする

```python
sc = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="position",
    mode="rms_weight",
    baseline_vwin=["-35:-5", "30:55"],
    baseline_poly=3,
    post_baseline_mode="inherit_all",
    post_baseline_vwin=["-40:-8", "32:58"],
    post_baseline_poly=5,
)
```

この例のように、`inherit_all` を有効にした上で、変えたい項目だけを明示値で上書きするのが一番誤解が少ない書き方です。`post_baseline_vwin="inherit"` などを併記する必要はありません。

### ポイント
- `post_baseline_mode` は現状 `None` または `inherit_all`
- `inherit_all` を使うとき、`post_baseline_vwin="inherit"` などの個別指定は冗長です
- 個別に変えたい項目があるときだけ、その項目に明示値を与えます
- `line_vwin` を与えると post-baseline 側にも差し引きが反映される
- 出力 row には `BSL_STAGE="post_coadd"` などが残る

---

## 19. `out_scale="TR*"` で出力する

### 目的
最終 coadd を `Tr*` として出力する。

```python
sc_tr = coadd.run_velocity_coadd(
    inputs="baseline_done.fits",
    mode="rms_weight",
    group_mode="position",
    rms_vwin=["-200:0"],
    out_scale="TR*",
    normalize_if_mixed="auto",
)
```

### ポイント
- `BEAMEFF` が無いと止まる
- `BEAMEFF` 混在グループは `normalize_if_mixed="auto"` が安全
- 実装上、`sigma_scale` は現状 `TA*` 以外は未実装

---

## 20. BEAMEFF 混在グループをどう扱うか

### 自動正規化

```python
sc_auto = coadd.run_velocity_coadd(
    inputs="baseline_done.fits",
    mode="rms_weight",
    group_mode="position",
    rms_vwin=["-200:0"],
    out_scale="TR*",
    normalize_if_mixed="auto",
)
```

### 正規化しない

```python
sc_never = coadd.run_velocity_coadd(
    inputs="baseline_done.fits",
    mode="rms_weight",
    group_mode="position",
    rms_vwin=["-200:0"],
    out_scale="TA*",
    normalize_if_mixed="never",
)
```

### ポイント
- `auto` では mixed group を `TR*` 正規化してから coadd し、出力 `BEAMEFF=1` 相当に整理する
- `never` は理論的に危険な場合がある

---

## 21. `set_beameff()` と `apply_relative_scale()` / `apply_global_scale()`

### 21.1 役割分担

- `set_beameff()` は `BEAMEFF` metadata を設定するだけ
- `apply_relative_scale()` は beam ごとの相対係数 $g_i$ を data へ直接**掛ける**
- `apply_global_scale()` は後で決まった一律 factor $f$ を data へ直接**掛ける**

$$
T_i^{\mathrm{aligned}}(k) = g_i\,T_i(k)
$$

$$
T_i^{\mathrm{final}}(k) = f\,T_i^{\mathrm{aligned}}(k)
$$

一方、物理 `BEAMEFF` は

$$
T_R^* = \frac{T_A^*}{\eta}
$$

の $\eta$ です。したがって、相対スケール合わせをしたいだけなら、`BEAMEFF` に流用せず `apply_relative_scale()` を使う方が安全です。

### 21.2 calibration 直後に適用するのが推奨

```python
sc_cal = cal.run_tastar_calibration(...)

su.apply_relative_scale(
    sc_cal,
    {
        (0, 0, 0): 1.00,
        (0, 0, 1): 0.97,
        (1, 0, 0): 1.03,
        (1, 0, 1): 1.01,
    },
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
)
```

この段階で掛けておくと、その後の baseline, `BSL_RMS`, moment0, coadd はすべて現在のスケールで再計算されます。

### 21.3 後日、一律 factor を掛ける

```python
su.apply_global_scale(sc_cal, 0.91)
```

この係数は generic な一律補正用です。後で決まった量が物理 `TR*` 変換の $\eta$ である場合は、

```python
su.set_beameff(sc_cal, 0.42)
```

のように `BEAMEFF` として与え、viewer / coadd / write の `TR*` 変換へ渡す方が意味論として自然です。`set_beameff()` 自体は data 本体を掛け直しません。

### 21.4 既に派生量があるとき

`apply_relative_scale()` / `apply_global_scale()` は既定で次の列を無効化します。

- `BSL_RMS`
- `BSL_COEF`
- `BSL_SCALE`
- `MOMENT0`
- `MOM0`
- `INTEGRATED`

したがって、scale を掛けた後は必要なら baseline や moment0 を再計算してください。

### 21.5 `set_beameff()` を使う場面

`TR*` 表示や `out_scale="TR*"` を使うために、物理 `BEAMEFF` を設定したいときです。

```python
su.set_beameff(sc_tr, 0.42)
```

### 21.6 multi-beam の mapping 指定

```python
su.set_beameff(
    sc_tr,
    {
        (0, 0, 0): 0.46,
        (0, 0, 1): 0.44,
        (1, 0, 0): 0.42,
        (1, 0, 1): 0.41,
    },
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
)
```

```python
su.apply_relative_scale(
    sc_tr,
    {
        (0, 0, 0): 1.00,
        (0, 0, 1): 0.97,
        (1, 0, 0): 1.03,
        (1, 0, 1): 1.01,
    },
    key_columns=("FDNUM", "IFNUM", "PLNUM"),
)
```

## 23. row 抽出を使う

### 目的
一部の scan だけ試験解析する。

```python
sc_sub = bsl.run_baseline_fit(
    input_data="coadded_before_baseline.fits",
    rows="0:20,50,60:80",
    vwin=["-35:-5", "30:55"],
    poly_order=5,
)
```

あるいは

```python
sc_sub = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    exclude_rows="0:100",
    mode="rms_weight",
    rms_vwin=["-200:0"],
)
```

### `select_rows()` を直接使う

```python
sc_head100 = su.select_rows(sc_tr, "0:100")
```

### 条件抽出で新しい Scantable を作る

```python
sc_on = su.filter_scantable(sc_tr, OBSMODE="ON")
```

```python
sc_scan = su.filter_scantable(sc_tr, query="SCAN >= 10 and SCAN <= 20")
```

### ポイント
- 文字列 selector は `0:20,50,60:80` のように書ける
- `rows` と `query` を役割で使い分けると混乱しにくい

---

## 24. `show_scantable()` で table を確認する

```python
su.show_scantable(
    inputs="tastar_dump.fits",
    rows="0:10",
    columns="SCAN, TIMESTAMP, OBSMODE, RESTFRQ, RA, DEC, VELOSYS, TEMPSCAL, BEAMEFF",
    head=100,
)
```

### ポイント
- 列名が table に無くても meta や alias から拾えることがある
- `RESTFRQ` / `RESTFREQ` の揺れ確認に便利
- `TIME=0` でも `TIMESTAMP` が補われることがある

### 全列をざっと見る

```python
su.show_scantable(sc_tr, rows="0:5", columns="all", head=5)
```

### 列の存在確認を先にする

```python
su.describe_columns(sc_tr)
```

---

## 25. 参照点からのオフセットをその場で見る

### `show_scantable()` に直接渡す

```python
su.show_scantable(
    inputs="tastar_dump.fits",
    rows="0:20",
    columns="TIMESTAMP, RA, DEC, OFS_LON, OFS_LAT",
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
)
```

### DataFrame として欲しい場合

```python
ofs = su.calc_mapping_offsets(
    sc=sc_tr,
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
    cos_mode="ref",
)
print(ofs.head())
```

### 位置で絞り込む

```python
ofs = su.calc_mapping_offsets(
    sc_tr,
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
)

sc_center = su.filter_scantable(
    sc_tr,
    query="OFS_LON**2 + OFS_LAT**2 < 60**2",
    extra_data=ofs,
)
```

### 重要
- `ref_coord="Orion KL"` のような名前解決は保証されない
- 確実なのは座標文字列か `SkyCoord`
- 現コードで確実に使える投影は `GLS/SFL` と `CAR/NONE`

---

## 26. 複数 Scantable を安全に結合する

```python
sc_merged = su.merge_scantables(
    inputs=["day1.fits", "day2.fits", "day3.fits"],
    sort_by_time=True,
    shift_scan_id=True,
)
```

### ポイント
- `shift_scan_id=True` が重要
- 後段で `group_mode="scan"` を使うなら scan 番号重複を避けられる
- `meta` は先頭ファイル基準なので、結合後に確認する

### 結合後の確認

```python
su.show_scantable(
    sc_merged,
    rows="0:10",
    columns="SCAN, TIMESTAMP, OBSMODE, OBJECT"
)
```

### metadata を限定的に修正する

```python
su.update_metadata(
    sc=sc_merged,
    column="OBJECT",
    value="Orion-KL",
    rows="0:100",
    verbose=True,
)
```

### `RESTFRQ` を強制修正する場合

```python
su.update_metadata(
    sc=sc_merged,
    column="RESTFRQ",
    value=115.2712018e9,
    force=True,
)
```

---

## 27. Standardizer を直接使う

### 目的
coadd の前に、再グリッド後の行列だけ欲しい。

```python
from sd_radio_spectral_fits.regrid_vlsrk import Standardizer

std = Standardizer(sc_tr, v_corr_col="VELOSYS")
mat, v = std.get_matrix()

print(mat.shape)
print(v[:5], v[-5:])
```

### グリッドを自分で指定

```python
mat, v = std.get_matrix(dv=0.2, vmin=-50, vmax=80)
```

### 一度決めた target grid を再利用する

```python
target = rv.make_vgrid(vmin_kms=-50, vmax_kms=80, dv_kms=0.2)
mat, v = std.get_matrix(target_grid=target)
```

### ポイント
- 異なる WCS を持つ入力でも共通 velocity matrix を作れる
- 独自統計や PCA などへつなげやすい
- `clear_caches()` でキャッシュを明示的に捨てられる

---

## 28. `run_velocity_regrid()` を単独で使う

### 目的
row を保ったまま、各スペクトルを共通 `LSRK` 速度グリッドへ再配置する。

```python
sc_rg = rv.run_velocity_regrid(
    input_data="tastar_dump.fits",
    vmin_kms=-50,
    vmax_kms=80,
    dv_kms=0.2,
    v_corr_col="VFRAME",
    overwrite=True,
)
```

### 一部 row だけ

```python
sc_rg = rv.run_velocity_regrid(
    input_data="tastar_dump.fits",
    rows="0:200",
    vmin_kms=-50,
    vmax_kms=80,
    dv_kms=0.2,
    ch_start=200,
    ch_stop=3800,
    keep_row_order=True,
    drop_allnan_rows=False,
)
```

### 出力の WCS

`axis_type="vel"` の coadd と同様に、regrid 後は概念的に

$$
{\rm CTYPE1} = {\rm VRAD}
$$

$$
{\rm CUNIT1} = {\rm m/s}
$$

$$
{\rm CRVAL1} = 1000 \times v_{\min}
$$

$$
{\rm CDELT1} = 1000 \times \Delta v
$$

です。

### ポイント
- 出力は row preserving
- `SPECSYS` は `LSRK` へ更新される
- `SSYSOBS` は観測時 frame を保持する。たとえば TOPOCENT 観測を regrid しても、通常 `SSYSOBS` は `TOPOCENT` のまま
- 未適用速度列は落とされる

---

## 29. `axis_type="freq"` と `axis_type="vel"` の違い

coadd では共通速度グリッドを内部で作りますが、最終出力の WCS を周波数で持つか速度で持つかを選べます。

### `axis_type="vel"`

出力は `VRAD/m/s` です。

### `axis_type="freq"`

出力は `FREQ/Hz` です。  
このとき、目標速度グリッド $v$ から周波数は

$$
\nu(v) = \nu_0 \left(1 - \frac{v}{c}\right)
$$

で計算されます。ここで $\nu_0$ は `RESTFRQ` です。

### 例

```python
sc_freq = coadd.run_velocity_coadd(
    inputs="tastar_dump.fits",
    group_mode="position",
    mode="rms_weight",
    rms_vwin=["-200:0"],
    axis_type="freq",
    vmin=-50,
    vmax=80,
    dv=0.2,
)
```

---

## 30. `restfreq.apply_restfreq_override()` を直接使う

### 目的
Scantable の `meta` と `table` に対して、静止周波数の上書きを自分で行う。

```python
sc = fitsio.read_scantable("some_file.fits")
info = rf.apply_restfreq_override(
    meta=sc.meta,
    table=sc.table,
    rest_freq_hz=109.7821734e9,
    require_wcs_for_vrad=True,
)
print(info)
```

### `FREQ` 軸と `VRAD` 軸の違い

- `FREQ` 軸では、周波数 WCS はそのままにして `RESTFRQ/RESTFREQ` だけ更新
- `VRAD` 軸では、速度軸 WCS も更新

`VRAD` 軸では

$$
v_2 = a v_1 + b
$$

$$
a = \frac{\nu_1}{\nu_2}
$$

$$
b = c \left(1 - \frac{\nu_1}{\nu_2}\right)
$$

です。

### ポイント
- `meta` だけでなく `table` に `CRVAL1`, `CDELT1`, `RESTFRQ` があればそちらも更新される
- `VRAD` 軸では既存 `RESTFRQ` が無いと止まる

---

## 31. RawSpec を自作する

### 目的
自分で `HOT / OFF / ON` DataFrame を持っている場合に RawSpec を作る。

```python
from sd_radio_spectral_fits.rawspec import build_rawspec

raw = build_rawspec(
    hot=hot_df,
    off=off_df,
    on=on_df,
    meta={
        "RESTFRQ": 115.2712018e9,
        "CRVAL1": 115.3e9,
        "CDELT1": -1.0e5,
        "CRPIX1": 1.0,
        "SPECSYS": "TOPOCENT",
        "SITELAT": 34.0,
        "SITELONG": 135.0,
        "SITEELEV": 100.0,
    },
    mapping=on_mapping_df,
)
```

### そのまま Ta* 化

```python
sc = cal.run_tastar_calibration(
    input_data=raw,
    t_hot_k=300.0,
    vlsrk_range_kms=(-200, 200),
)
```

### ポイント
- `hot`, `off`, `on` は `DatetimeIndex` を持つ DataFrame が必要
- `mapping` に `OBSMODE` が無ければ ON-only mapping と解釈される
- `mapping_all` には `HOT/OFF/ON` が統合される

---

## 32. RawSpec を保存・再読込する

### pickle で保存する

```python
rawspec.save_rawspec(raw, "rawspec.pkl")
```

### 自動読込

```python
raw2 = rawspec.load_rawspec_auto("rawspec.pkl")
raw3 = rawspec.load_rawspec_auto("raw_input.fits")
```

### ポイント
- `load_rawspec_auto()` は拡張子や FITS マジックから自動判定する
- FITS SDFITS-like も、pickle 形式も扱える

---

## 33. `read_scantable()` の legacy マイグレーションを知っておく

最新版の `read_scantable()` では、旧列をかなり自動整理します。

### 代表例

- `V_CORR_KMS` を `VFRAME` と `VELOSYS` の `m/s` 列へ移す
- `TAMB_K` を `THOT` へ寄せる
- `TAU` を `TAU0` へ寄せる

### 例

```python
sc = fitsio.read_scantable("legacy_file.fits")
print(sc.history.get("fitsio_migration"))
```

### ポイント
- 旧データを読みやすくするための救済であり、以後の公開意味論は `VELOSYS` / `VFRAME` / `THOT` / `TAU0`
- ただし列が本質的に矛盾している場合は止まる

---

## 34. SDFITS として保存する

```python
fitsio.write_scantable(
    path="final_12co.fits",
    scantable=sc_final_bsl,
    overwrite=True,
)
```

### 書き出し時だけ `TR*` にしたい

```python
fitsio.write_scantable(
    path="final_12co_tr.fits",
    scantable=sc_final_bsl,
    overwrite=True,
    tempscal="TR*",
    data_scale="TA*",
)
```

### 数式で書くと

保存時変換では

$$
T_R^* = \frac{T_A^*}{\eta_{\rm B}}
$$

または

$$
T_A^* = T_R^* \eta_{\rm B}
$$

を **書き出し時だけ** 適用します。メモリ上の `scantable.data` は壊しません。

### ポイント
- `tempscal` または `out_scale` が出力ラベル
- `data_scale` がメモリ上の実体のスケール
- 変換が実際に走った場合は `history` に追記される

---

## 35. `write_sdfits()` を直接使う

### 目的
`meta`, `data`, `table`, `history` を直接与えて SDFITS-like FITS を書く。

```python
fitsio.write_sdfits(
    out_path="manual_write.fits",
    meta=sc_final_bsl.meta,
    data=sc_final_bsl.data,
    table=sc_final_bsl.table,
    history=sc_final_bsl.history,
    spectrum_column="DATA",
    include_flag=True,
    overwrite=True,
)
```

### 文字列幅を明示する

```python
fitsio.write_sdfits(
    out_path="manual_write.fits",
    meta=sc_final_bsl.meta,
    data=sc_final_bsl.data,
    table=sc_final_bsl.table,
    history=sc_final_bsl.history,
    string_widths={"OBJECT": 32, "OBSMODE": 24},
    overwrite=True,
)
```

### ポイント
- `data` は 2 次元配列でも `list[np.ndarray]` でもよい
- 後者なら VLA として書かれる
- table の vector-in-cell も fixed vector か VLA として保存される

---

## 36. low-level writer を使って大きな raw/OTF を parts 書き出しする

現アーカイブには low-level の `SDRadioSpectralSDFITSWriter` も含まれています。  
これは converter 系や巨大 OTF 出力向けの writer です。

### 最小例

```python
from sd_radio_spectral_fits.sdfits_writer import (
    SDRadioSpectralSDFITSWriter,
    Site,
    DatasetInfo,
    SpectralAxisUniform,
)

site = Site(lat_deg=34.0, lon_deg=135.0, elev_m=100.0)
ax = SpectralAxisUniform(
    crval1_hz=115.3e9,
    cdelt1_hz=-1.0e5,
    crpix1=1.0,
    restfreq_hz=115.2712018e9,
    specsys="TOPOCENT",
    ssysobs="TOPOCENT",
    veldef="RADIO",
)
info = DatasetInfo(
    telescope="ExampleScope",
    observer="User",
    project="TEST",
    object_name="Orion-KL",
    spectral_axis=ax,
)

writer = SDRadioSpectralSDFITSWriter(
    n_chan=4096,
    site=site,
    info=info,
    store_freq_column=False,
)
```

### parts mode

```python
writer = SDRadioSpectralSDFITSWriter(
    n_chan=4096,
    site=site,
    info=info,
    store_freq_column=False,
    chunk_size=10000,
    out_basename="otf_run1",
)
```

### ポイント
- `chunk_size` を指定すると `*_part0000.fits` 形式で分割される
- 最後に manifest JSON が書かれる
- `store_freq_column=False` なら uniform WCS が前提
- `SPECSYS="LSRK"` の uniform axis を直接書くなら `store_freq_column=True` が必要

---

## 37. viewer 系を簡単に使う

> 注記  
> `profile_view` の詳細な API、キーボード操作、`viewer` / `montage` / `grid` の内部仕様は `profile_view_manual_ja_v2.md` を参照してください。ここでは解析後確認にすぐ使える簡単な例だけを載せます。

### 1 本ずつ見る

```python
pv.view_spectra(
    sc_final_bsl,
    xaxis="vel",
    xrange=(-30, 55),
)
```

### montage を開く

```python
pv.plot_profile_map(
    sc_final_bsl,
    mode="montage",
    xaxis="vel",
    nrows=5,
    ncols=4,
    xrange=(-10, 25),
)
```

### grid を開く

```python
pv.plot_profile_map(
    sc_final_bsl,
    mode="grid",
    coord="radec",
    projection="SFL",
    ref_point=(83.822, -5.391),
    x_grid=30.0,
    y_grid=30.0,
    xaxis="vel",
    xrange=(-10, 25),
)
```

### ポイント

- 詳細仕様の正本は `profile_view_manual_ja_v2.md`
- cookbook 側では最上位 API の簡単な呼び出しだけ覚えればよい
- row ごとの WCS を尊重して表示する点が重要

## 38. よくある失敗例と回避法

### 38.1 `baseline_vwin` と `rms_vwin` を同時指定する
不可です。どちらか一方にします。

### 38.2 `coadd_qc` を使うのに `mode="uniform"` にしている
不可です。`mode="rms"` か `mode="rms_weight"` にします。

### 38.3 `TOPOCENT` なのに座標も `VELOSYS` も無い
`LSRK` 補正を作れません。  
元データに `VELOSYS` / `VFRAME` があるか、RA/DEC と site 情報を補います。

### 38.4 `TR*` を出したいのに `BEAMEFF` が無い
先に `set_beameff()` してください。

### 38.5 同じ raw から 13CO/C18O を解析するときに、同じ rest frequency のまま処理してしまう
`rest_freq=` を明示してください。

### 38.6 scan 単位でまとめずにいきなり最終統合して baseline が不安定
まず `group_mode="scan"` を試してください。

### 38.7 `sigma_scale="TR*"` を指定したい
現実装では未実装です。重み評価は `TA*` を前提に考えます。

### 38.8 `read_scantable()` が `TR*` を勝手に `TA*` に戻すと思い込む
最新版では戻しません。`TEMPSCAL` を尊重します。

---

## 39. 研究メモとしての推奨レシピ

### 12CO の無難な初期値

```python
vlsrk_range_kms = (-200, 200)
vwin_rms = "-200:0"
vwin_baseline = ["-35:-5", "30:55"]
poly_order = 7
```

### 13CO / C18O の初期値例

```python
f_13co = 110.2014e9
f_c18o = 109.7821734e9

vwin_baseline_13 = ["-10:5", "15:30"]
vwin_baseline_c18o = ["-10:2.5", "15:27"]
```

### なぜこのような二段階か

- raw dump はノイジーで baseline 不安定
- scan coadd するとベースライン形状が見えやすい
- final coadd 前後で baseline を 1 回ずつ入れると、全体として安定しやすい

---

## 40. 最後に

この cookbook の本質は、単にコマンド例を並べることではなく、

- どの段階で何を固定し
- どの段階で何を再解釈し
- どの段階で何を平均するか

を明確にすることにあります。

特に、毎回意識するとよいのは次の 7 点です。

- `rest_freq`
- `SPECSYS`
- `VELOSYS` / `VFRAME`
- `baseline_vwin` / `rms_vwin`
- `TEMPSCAL` / `BEAMEFF`
- `axis_type`
- `out_scale`

この 7 点を固定して考えると、かなり複雑なデータでも混乱が減ります。
