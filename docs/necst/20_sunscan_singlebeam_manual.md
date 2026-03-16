# sunscan singlebeam 詳細マニュアル

対象スクリプト:

- `multibeam_beam_measurement/sunscan_singlebeam.py`

---

## 1. singlebeam の目的

singlebeam は、1 本の spectrometer stream について太陽スキャンを解析し、主に次を求めます。

- スキャン中心位置
- HPBW
- スキャン速度
- limb fit の成功 / 失敗
- ripple 除去の影響
- trim の影響
- 解析に使った time offset や Az/El source

singlebeam の出力は、

- 1 stream の品質確認
- multi-beam の前段確認
- pseudo multibeam の元データ

として使えます。

---

## 2. 入力から出力までの流れ

概念的な流れは次です。

1. CLI を読む
2. 必要なら spectrometer config TOML を読み、`stream-name` で 1 stream を選ぶ
3. `SunScanAnalysisConfig` を構築する
4. NECST DB から spectrum / encoder / altaz / mode / weather を読む
5. spectrometer / encoder の時刻補正を適用する
6. 太陽中心に対する main-axis / cross-axis offset 系列を作る
7. 必要なら chopper wheel を適用して `Ta*` 相当を作る
8. ripple 除去を行う
9. scan trim を行う
10. limb fit / derivative fit を行う
11. summary CSV と derivative plot を出力する

---

## 3. `SunScanAnalysisConfig` の構造

singlebeam / multibeam extract は、共通の解析設定 dataclass を使います。

### 3.1 InputConfig

入力と座標・時刻に関する設定です。

- `rawdata_path`
- `db_namespace`
- `telescope`
- `tel_loaddata`
- `planet`
- `spectral_name`
- `azel_source`
- `altaz_apply`
- `spectrometer_time_offset_sec`
- `encoder_shift_sec`
- `encoder_az_time_offset_sec`
- `encoder_el_time_offset_sec`
- `encoder_vavg_sec`

### 3.2 CalibrationConfig

- `chopper_wheel`
- `tamb_k`
- `chopper_win_sec`
- `chopper_stat`

### 3.3 RippleConfig

- `enabled`
- `preset`
- `model`
- `target_hz`
- `search_hz`
- `bw_hz`
- `max_harm`
- `order`
- `notch_passes`
- `trend_win_sec`
- `resample_dt_sec`
- `eval_band_hz`

### 3.4 ProfileConfig

- `profile_xlim_deg`

### 3.5 TrimConfig

- `enabled`
- `vfrac`
- `vmin`
- `gap_fill`
- `min_samples`
- `dominant_axis`
- `ratio_min`
- `vpercentile`
- `steady_scan`
- `use_on_only`
- `xwin_factor`
- `cross_offset_max_deg`
- `speed_min_deg_s`
- `steady_cv_max`

### 3.6 EdgeFitConfig

- `enabled`
- `strict_deriv`
- `fit_win_deg`
- `fit_threshold`
- `hpbw_init_arcsec`

### 3.7 ReportConfig

- `outdir`
- `debug_plot`
- `edge_fit_plot_max_scans`
- `tag`

### 3.8 RuntimeConfig

- `continue_on_error`

### 3.9 BeamOverride

spectrometer config から取った stream metadata を解析結果へ移すための薄い構造です。

- `stream_name`
- `beam_id`
- `restfreq_hz`
- `hpbw_init_arcsec`
- `polariza`
- `fdnum`
- `ifnum`
- `plnum`
- `sampler`

---

## 4. CLI パラメータ完全説明

### 4.1 入力 / 出力 / stream 選択

#### `rawdata`

RawData ディレクトリ。

#### `--outdir`

出力ディレクトリ。

#### `--spectrometer-config`

converter-compatible spectrometer TOML。
与えた場合、singlebeam でも converter と同じ stream 定義を共有できます。

#### `--stream-name`

`--spectrometer-config` を使うときに、解析する 1 stream を選びます。

#### `--db-namespace`

DB table 名の namespace。
TOML の `global.db_namespace` と共有できます。

#### `--telescope`

DB table 名の telescope 名。

#### `--tel-loaddata`

legacy loader 互換の telescope 名。

#### `--planet`

太陽以外の天体へ拡張したい場合の body 名ですが、通常は `sun` を使います。

#### `--spectral-name`

spectrometer config を使わないときの stream 名。

---

### 4.2 Az/El / time offset

#### `--azel-source`

解析に使う Az/El source。

- `encoder`
- `altaz`

現行 singlebeam では、この 2 択が重要です。

#### `--altaz-apply`

altaz source をどう適用するか。

- `none`
- `minus`
- `plus`

実データに応じた補正の符号合わせに使います。

#### `--spectrometer-time-offset-sec`

spectrometer の時刻補正 [s]。

#### `--encoder-shift-sec`

encoder 共通の時刻シフト [s]。

#### `--encoder-az-time-offset-sec`

encoder Az 系列の追加オフセット [s]。

#### `--encoder-el-time-offset-sec`

encoder El 系列の追加オフセット [s]。

#### `--encoder-vavg-sec`

encoder 系列の時間平均窓 [s]。

---

### 4.3 chopper wheel / calibration

#### `--no-chopper-wheel`

chopper wheel を無効にします。
デフォルトは有効です。

#### `--tamb-k`

Tamb [K] を明示します。
省略時は DB / 推定値側を使います。

#### `--chopper-win-sec`

HOT / OFF / ON の周辺から統計量を取る時間窓 [s]。

#### `--chopper-stat`

chopper wheel の代表値の取り方。

- `median`
- `mean`

ノイズや外れ値に強いのは通常 `median` です。

---

### 4.4 profile / ripple

#### `--profile-xlim-deg`

profile plot / fit の横軸表示範囲 [deg]。

#### `--ripple-no-remove`

ripple 除去を無効にします。

#### `--ripple-preset`

preset を使います。

- `auto`
- `safe`
- `normal`
- `strong`

#### `--ripple-model`

ripple モデル。

- `auto`
- `add`
- `mul`

#### `--ripple-target-hz`

除去したい主 ripple 周波数 [Hz]。

#### `--ripple-search-hz`

探索帯域幅 [Hz]。

#### `--ripple-bw-hz`

notch の帯域幅 [Hz]。

#### `--ripple-max-harm`

扱う高調波数。

#### `--ripple-order`

notch / model の次数。

#### `--ripple-notch-pass`

notch pass 回数。

#### `--ripple-trend-win-sec`

trend 推定の窓 [s]。

#### `--ripple-resample-dt-sec`

resample の時間刻み [s]。

#### `--ripple-eval-band-hz`

評価に使う帯域 [Hz]。

---

### 4.5 edge fit

#### `--no-edge-fit`

edge fit を無効にします。

#### `--edge-fit-win-deg`

limb fit に使うウィンドウ幅 [deg]。

#### `--edge-fit-threshold`

fit 対象抽出の閾値。

#### `--hpbw-init-arcsec`

HPBW の初期値 [arcsec]。
最適化の初期推定として使われます。

#### `--edge-fit-plot-max-scans`

plot に載せる scan 数上限。

---

### 4.6 trim

#### `--no-trim-scan`

scan trim を無効にします。

#### `--trim-vfrac`

速度や微分に対する閾値の相対係数。

#### `--trim-vmin`

trim に使う最低速度閾値。

#### `--trim-gap`

gap の許容サンプル数。

#### `--trim-min-samples`

最低サンプル数。

#### `--trim-dominant-axis` / `--trim-no-dominant-axis`

支配軸ベースの trim を有効 / 無効にします。

#### `--trim-axis-ratio-min`

支配軸と副軸の速度比閾値。

#### `--trim-vpercentile`

速度 percentiles を使うときの百分位。

#### `--trim-no-steady-scan`

steady-scan 条件を無効にします。

#### `--trim-include-hotoff`

trim の対象に HOT/OFF も含めます。
既定では ON 優先です。

#### `--trim-scan-speed-min-arcsec`

最低スキャン速度 [arcsec/s]。

#### `--trim-xwin-factor`

main-axis window の余裕係数。

#### `--trim-cross-offset-max-deg`

cross-axis 許容幅 [deg]。

#### `--trim-steady-cv-max`

steady-scan 判定の変動係数上限。

---

### 4.7 strictness / runtime / optics

#### `--strict-deriv` / `--no-strict-deriv`

derivative fit の厳しさを切り替えます。

#### `--continue-on-error`

一部 scan に失敗しても継続します。

#### `--debug-plot`

debug plot を追加で出力します。

#### `--dish-diameter-m`

望遠鏡口径 [m]。

#### `--hpbw-factor`

口径から HPBW を見積もるときの比例係数です。

---

## 5. spectrometer config と singlebeam

singlebeam は spectrometer config を使うと、converter と同じ `[[spectrometers]]` 群から 1 stream を選べます。

### 5.1 `global` から読むもの

現行 singlebeam は、config が与えられると次を `global` から読めます。

- `db_namespace`
- `telescope`
- `tel_loaddata`
- `planet`
- `spectrometer_time_offset_sec`
- `encoder_shift_sec`
- `encoder_az_time_offset_sec`
- `encoder_el_time_offset_sec`

CLI で同名オプションを明示した場合は、そちらが優先です。

### 5.2 stream から読むもの

`stream-name` で選んだ stream から、少なくとも次が解析結果に引き継がれます。

- `stream_name`
- `beam_id`
- `restfreq_hz`
- `polariza`
- `fdnum`
- `ifnum`
- `plnum`
- `sampler`

---

## 6. singlebeam の出力物

主な出力は次です。

### `sun_scan_summary_<tag>.csv`

もっとも重要な出力です。

代表的な列:

- `scan_id`
- `rep_az_deg`, `rep_el_deg`
- `center_az_deg`, `center_el_deg`
- `hpbw_az_arcsec`, `hpbw_el_arcsec`
- `fit_ok_az`, `fit_ok_el`
- `az_error`, `el_error`
- `speed_az_arcsec_s`, `speed_el_arcsec_s`
- `spec_time_basis`, `spec_time_suffix`, `spec_time_fallback_field`, `spec_time_example`
- `azel_source`, `altaz_apply`
- `spectrometer_time_offset_sec`
- `encoder_shift_sec`, `encoder_az_time_offset_sec`, `encoder_el_time_offset_sec`, `encoder_vavg_sec`
- `db_namespace`, `telescope`, `tel_loaddata`, `planet`

### derivative fit PNG 群

edge fit / derivative fit の確認図です。
scan ごとの profile, derivative, fit が出ます。

### `summary_text_<tag>.png`

`debug_plot` 有効時に summary 文字列を画像化したものです。

---

## 7. singlebeam 解析の解釈

### 7.1 `rep_az_deg`, `rep_el_deg`

各 scan の代表的 Az/El です。
通常は zero main-axis offset 近傍の tracking 点から代表値を取ります。

### 7.2 `center_az_deg`, `center_el_deg`

太陽中心に対する fitted center です。
0 に近いほど中心合わせが良いことを示します。

### 7.3 `hpbw_*_arcsec`

各軸の半値全幅です。

### 7.4 `fit_ok_az`, `fit_ok_el`

各軸の fit 成功 / 失敗です。

### 7.5 time offset 列

summary CSV には、解析時に実際に使った spectrometer / encoder offset が書かれます。
これは、あとで multi-beam 結果と照合する際に非常に重要です。

---

## 8. singlebeam 特有の注意

### 8.1 `spectrometer-config` を使うとき

stream 名を必ず明示してください。
1 本だけしか無い config なら省略でもよいことがありますが、複数 stream 定義がある場合は混乱の元です。

### 8.2 `azel_source=altaz`

指令値側を使うため、encoder 系とは違う見え方になります。
符号合わせのために `altaz_apply` を確認してください。

### 8.3 ripple を強くしすぎない

ripple 除去は有効ですが、強くしすぎると profile 自体を崩します。
まずは `preset=auto` または `safe` を推奨します。

### 8.4 trim が強すぎる場合

scan の端が削られすぎると HPBW や中心が不安定になります。
`trim-vfrac`, `trim-min-samples`, `trim-cross-offset-max-deg` を確認してください。
