# 単一鏡電波分光解析パッケージ Cookbook
## 実用例・推奨レシピ集（日本語）

> 方針  
> `obs_method.txt` にあったループ主体の書き方は、ここでは **意味が読みやすい単発レシピ** に直しています。  
> 複数ファイル処理をしたい場合は、まず単発例が理解できることを優先してください。

---

## 0. 事前準備

多くの例では、まず次を入れておくと扱いやすいです。

```python
import astropy.utils.iers
astropy.utils.iers.conf.auto_download = False
astropy.utils.iers.conf.auto_max_age = None
astropy.utils.iers.conf.iers_degraded_accuracy = "ignore"

import sd_radio_spectral_fits.calibrate as cal
import sd_radio_spectral_fits.coadd as coadd
import sd_radio_spectral_fits.baseline as bsl
import sd_radio_spectral_fits.fitsio as fitsio
import sd_radio_spectral_fits.scantable_utils as su

# profile_view（旧 plotting）
import sd_radio_spectral_fits.profile_view.viewer as pv
import sd_radio_spectral_fits.profile_view.montage as mnt
import sd_radio_spectral_fits.profile_view.grid as grd
```

---

## 1. 最小の 12CO 解析

### 目的
1 本の raw FITS から Ta* を作り、scan coadd、baseline、最終 coadd まで行う。

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

# 7. 表示
pv.view_spectra(sc_final_bsl, xrange=(-30, 55))
```

### ポイント
- `group_mode="scan"` により、まず scan 内の dump をまとめる
- そのあと baseline
- 最後に同位置の repeat をまとめる

---

## 2. 複数ファイルの 12CO を統合する

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

# 各ファイルを個別に「scan coadd → baseline → file内最終coadd」まで済ませる
# ここでは明示的に3本書く
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

# 全ファイル統合
sc_total = coadd.run_velocity_coadd(
    inputs=[sc1, sc2, sc3],
    baseline_poly=npoly,
    baseline_vwin=vwin_baseline,
    mode="rms_weight",
    pos_tol_arcsec=1.0,
    overwrite=True,
)

pv.view_spectra(sc_total, xrange=(-30, 55))
```

### ポイント
- cookbook では loop を外して書いた
- 実運用では list やループへ戻してよい
- `inputs=[sc1, sc2, sc3]` としてそのまま最終統合できる

---

## 3. 13CO / C18O を同じ raw から分けて解析する

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

# 13CO
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

# C18O
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

pv.view_spectra(sc_13, xrange=(-5, 25))
pv.view_spectra(sc_c18o, xrange=(-5, 25))
```

### ポイント
- FREQ 軸なら、rest frequency を変えても元の周波数 WCS を乱暴に壊さない
- 13CO と C18O で baseline 窓や `rms_bin` を変えられる

---

## 4. ATM モデルを使った Ta* キャリブレーション

### 目的
`tau_zenith` と外気温から `T_cal` をより厳密に計算する。

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

### もし header に TAU0 が無いなら
```python
sc_atm = cal.run_tastar_calibration(
    input_data="orikl_raw_12co_20260222_111212.fits",
    t_hot_k=300.0,
    tau_zenith=0.12,
    t_surface_k=278.0,
    t_atm_model="ratio",
    t_atm_eta=0.95,
)
```

### ポイント
- `tau_zenith=None` なら Basic 1-Temp
- `tau_zenith="auto"` は `TAU0` がないと止まる
- `t_surface_k` 未指定なら `WXTEMP`, `TAMBIENT`, `TAMB`, `T_SURF` を探す

---

## 5. 既存 Ta* を再キャリブレーションする

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
- 解析条件の比較に便利

---

## 6. baseline のみ単独でかける

### 目的
すでに coadd 済みのデータへ baseline だけかける。

```python
sc_bsl = bsl.run_baseline_fit(
    input_data="coadded_before_baseline.fits",
    vwin=["-35:-5", "30:55"],
    poly_order=7,
    line_vwin=["5:20"],  # 線を除外したい場合
    iter_max=3,
    iter_sigma=3.0,
    overwrite=True,
)
```

### ポイント
- `line_vwin` は baseline 窓から差し引かれる
- `iter_max > 0` にすると sigma clipping を回せる

---

## 7. coadd を「scan 単位」だけで使う

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

---

## 8. `baseline_vwin` を使う重み付き coadd

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
- baseline subtraction と RMS 評価を同時にしたいときはこちら
- `rms_vwin` と同時には使えない

---

## 9. `rms_vwin` を使う重み付き coadd

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

---

## 10. `BSL_RMS` を使って coadd する

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

---

## 11. QC モードを使う

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

### ポイント
- `coadd_qc` を使うと `mode="uniform"` は不可
- `rms_vwin` が必須

---

## 12. `out_scale="TR*"` で出力する

### 目的
最終 coadd を Tr* として出力する。

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
- `BEAMEFF` 混在グループは `auto` が安全

---

## 13. BEAMEFF をあとから設定する

`set_beameff()` は data 配列を変換せず、行ごとの `BEAMEFF` を table に与える関数です。  
したがって、まず効率を入れておき、その後 viewer や `out_scale="TR*"` で利用する、という流れになります。

### 基本形: 全行へ同じ効率を付与
```python
import sd_radio_spectral_fits.scantable_utils as su

su.set_beameff(sc_tr, efficiency=0.42)
```

### 一部行だけ
```python
su.set_beameff(sc_tr, efficiency=0.42, rows="0:20")
```

### ポイント
- data 実体は変わらない
- `TEMPSCAL` も変わらない
- `BEAMEFF` 列が無ければ自動作成される
- 0 より小さい値や 1 より大きい値には warning が出る

### 重要
この添付ソースの `set_beameff()` に **`target_scale` 引数はありません。**  
したがって
```python
# この形は現添付コードの API ではない
sd.set_beameff(sc, efficiency=0.45, target_scale="TR*")
```
のような呼び方は、そのままでは一致しません。  
現在の実装では、`BEAMEFF` を付与した後に viewer 側の TR* 表示、または coadd / write 側の scale 指定を使います。

---

## 14. FDNUM / IFNUM / PLNUM ごとに異なる BEAMEFF を設定する

多ビーム・多 IF・多偏波データでは、この使い方が最も重要です。

### 1組だけ設定
```python
idx = su.find_scans(sc_tr, FDNUM=0, IFNUM=1, PLNUM=0)
su.set_beameff(sc_tr, efficiency=0.45, rows=idx)
```

### 複数組に別々の値を設定
```python
beam_eff_map = {
    (0, 0, 0): 0.46,
    (0, 0, 1): 0.44,
    (1, 0, 0): 0.42,
    (1, 0, 1): 0.41,
}

for (fdnum, ifnum, plnum), eta in beam_eff_map.items():
    idx = su.find_scans(sc_tr, FDNUM=fdnum, IFNUM=ifnum, PLNUM=plnum)
    su.set_beameff(sc_tr, efficiency=eta, rows=idx, verbose=True)
```

### query でまとめて選ぶ
```python
idx = su.find_scans(
    sc_tr,
    query="FDNUM == 0 and IFNUM == 1 and PLNUM in [0, 1]"
)
su.set_beameff(sc_tr, efficiency=0.43, rows=idx)
```

### 設定結果を確認する
```python
su.show_scantable(
    sc_tr,
    rows="0:20",
    columns="FDNUM, IFNUM, PLNUM, TEMPSCAL, BEAMEFF"
)
```

---

## 15. row 抽出を使う

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

---

## 16. `show_scantable()` で table を確認する

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

## 17. 参照点からのオフセットをその場で見る

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
- `ref_coord="Orion KL"` のような名前解決は、現コードの `calc_mapping_offsets()` では保証されない
- 確実なのは座標文字列か `SkyCoord`
- 現コードで確実に使える投影は `GLS/SFL` と `CAR/NONE`

---

## 18. 複数 Scantable を安全に結合する

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

## 19. Standardizer を直接使う

### 目的
coadd の前に、再グリッド後の行列だけ欲しい。

```python
from sd_radio_spectral_fits.regrid_vlsrk import Standardizer

std = Standardizer(sc_tr, v_corr_col="VELOSYS")


### 目的
coadd の前に、再グリッド後の行列だけ欲しい。

```python
from sd_radio_spectral_fits.regrid_vlsrk import Standardizer

std = Standardizer(sc_tr, v_corr_col="VELOSYS")

# 自動グリッド
mat, v = std.get_matrix()

print(mat.shape)
print(v[:5], v[-5:])
```

### グリッドを自分で指定
```python
mat, v = std.get_matrix(dv=0.2, vmin=-50, vmax=80)
```

### ポイント
- 異なる WCS を持つ入力でも共通 velocity matrix を作れる
- 独自統計や PCA などへつなげやすい

---

## 20. 生データから RawSpec を自作する

### 目的
自分で HOT / OFF / ON DataFrame を持っている場合に RawSpec を作る。

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

---

## 21. SDFITS として保存する

```python
fitsio.write_scantable(
    path="final_12co.fits",
    scantable=sc_final_bsl,
    overwrite=True,
)
```

### 書き出し時だけ TR* にしたい
```python
fitsio.write_scantable(
    path="final_12co_tr.fits",
    scantable=sc_final_bsl,
    overwrite=True,
    tempscal="TR*",
    data_scale="TA*",
)
```

---

## 22. viewer 系を簡単に使う

### 1 本ずつ見る
```python
pv.view_spectra(
    sc_final_bsl,
    xrange=(-30, 55),
)
```

### montage
```python
viewer = mnt.ProfileMapMontageViewer(
    sc_final_bsl,
    xaxis="vel",
    nrows=5,
    ncols=4,
    xrange=(-10, 25),
)
```

### grid
```python
grid = grd.ProfileMapGridViewer(
    sc_final_bsl,
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
- profile_view は旧 plotting 相当
- 表示系の完全な説明は別冊向き
- ただし解析後確認として非常に重要

---

## 23. よくある失敗例と回避法

### 23.1 `baseline_vwin` と `rms_vwin` を同時指定する
不可です。どちらか一方にします。

### 23.2 TOPOCENT なのに RA/DEC がない
`VELOSYS` を再計算できません。  
元データに `VELOSYS` があるなら保持し、ないなら座標列と site 情報を補います。

### 23.3 TR* を出したいのに `BEAMEFF` が無い
先に `set_beameff()` してください。

### 23.4 同じ raw から 13CO/C18O を解析するときに、同じ rest frequency のまま処理してしまう
`rest_freq=` を明示してください。

### 23.5 scan 単位でまとめずにいきなり最終統合して baseline が不安定
まず `group_mode="scan"` を試してください。

---

## 24. 研究メモとしての推奨レシピ

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

## 25. 最後に

この cookbook の本質は、単にコマンド例を並べることではなく、

- どの段階で何を固定し
- どの段階で何を再解釈し
- どの段階で何を平均するか

を明確にすることにあります。

特に

- `rest_freq`
- `SPECSYS`
- `VELOSYS`
- `baseline_vwin` / `rms_vwin`
- `TEMPSCAL` / `BEAMEFF`

の 5 点は、毎回意識して使うことを強く勧めます。
