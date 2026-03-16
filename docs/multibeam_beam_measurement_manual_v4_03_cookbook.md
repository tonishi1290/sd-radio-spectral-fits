# Multi-beam 用 beam 測定パッケージ 説明書 v4-03  
## 詳細 cookbook・TOML 例・運用例集

## 1. この分冊の目的

本書は、現行 multi-beam package を**実際に動かすための例集**です。  
方針は次の通りです。

- 実運用でそのまま写しやすい
- single / extract / fit / pseudo / config 検査の流れが分かる
- converter-compatible TOML を使う
- usage flags を活用した stream 管理を例示する
- converter と sunscan の責務は分けつつ、**shared config** は分かるようにする

---

## 2. まず最初のおすすめ手順

1. `check_spectrometer_config` で config を検査  
2. `sunscan_singlebeam` で 1 本だけ成功させる  
3. `sunscan_extract_multibeam` を回す  
4. `sunscan_fit_multibeam --model both` を回す  
5. 必要なら `beam_model_*.toml` を converter へ反映  
6. 実 multi-beam データが無い段階では `sunscan_make_pseudo_multibeam` で dry-run

---

## 3. 最小 single-beam 用 TOML 例

ファイル名例: `single_beam_example.toml`

```toml
[global]
db_namespace = "necst"
telescope = "OMU1P85M"
tel_loaddata = "OMU1p85m"
planet = "sun"

spectrometer_time_offset_sec = -0.070
encoder_shift_sec = 0.0
encoder_az_time_offset_sec = 0.000
encoder_el_time_offset_sec = 0.040

[[spectrometers]]
name = "single115_xx"
fdnum = 0
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B00"
db_stream_name = "xffts-board1"

enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "none"
reference_angle_deg = 0.0
rotation_sign = 1.0
dewar_angle_deg = 0.0
```

### この例の意図
- 1 stream だけ定義
- global で time offset も定義
- singlebeam / extract / converter の shared config として使いやすい形

---

## 4. single-beam を 1 本だけ試す

### 4.1 config TOML を使う基本形

```bash
sunscan_singlebeam RAWDATA \
  --spectrometer-config single_beam_example.toml \
  --stream-name single115_xx \
  --outdir out_single
```

### 4.2 module 実行の形

```bash
python -m multibeam_beam_measurement.sunscan_singlebeam RAWDATA \
  --spectrometer-config single_beam_example.toml \
  --stream-name single115_xx \
  --outdir out_single
```

### 4.3 CLI で time offset を一時上書き

```bash
sunscan_singlebeam RAWDATA \
  --spectrometer-config single_beam_example.toml \
  --stream-name single115_xx \
  --spectrometer-time-offset-sec -0.050 \
  --encoder-el-time-offset-sec 0.020 \
  --outdir out_single_override
```

### 4.4 どう確認するか
- `sun_scan_summary_<tag>.csv`
- derivative fit PNG
- debug plot（必要なら）
- 標準出力 summary

---

## 5. multi-beam 用 TOML 例
以下は

- 中心 beam `B00` が 230 GHz
- 周辺 4 beam `B01`〜`B04` が 115 GHz
- 各 beam に XX / YY
- fit の primary は XX
- converter-compatible で共有しやすい

という想定です。

ファイル名例: `multibeam_5beam_example.toml`

```toml
[global]
db_namespace = "necst"
telescope = "OMU1P85M"
tel_loaddata = "OMU1p85m"
planet = "sun"

spectrometer_time_offset_sec = -0.070
encoder_shift_sec = 0.0
encoder_az_time_offset_sec = 0.000
encoder_el_time_offset_sec = 0.040

[[spectrometers]]
name = "center230_xx"
fdnum = 0
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B00"
db_stream_name = "xffts-board0"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 230.538e9

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "center230_yy"
fdnum = 0
ifnum = 0
plnum = 1
polariza = "YY"
beam_id = "B00"
db_stream_name = "xffts-board1"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = false
beam_fit_use = false

[spectrometers.frequency_axis]
restfreq_hz = 230.538e9

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "east115_xx"
fdnum = 1
ifnum = 1
plnum = 0
polariza = "XX"
beam_id = "B01"
db_stream_name = "xffts-board2"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 300.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "east115_yy"
fdnum = 1
ifnum = 1
plnum = 1
polariza = "YY"
beam_id = "B01"
db_stream_name = "xffts-board3"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = false
beam_fit_use = false

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 300.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "north115_xx"
fdnum = 2
ifnum = 1
plnum = 0
polariza = "XX"
beam_id = "B02"
db_stream_name = "xffts-board4"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 300.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "north115_yy"
fdnum = 2
ifnum = 1
plnum = 1
polariza = "YY"
beam_id = "B02"
db_stream_name = "xffts-board5"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = false
beam_fit_use = false

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 300.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "west115_xx"
fdnum = 3
ifnum = 1
plnum = 0
polariza = "XX"
beam_id = "B03"
db_stream_name = "xffts-board6"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = -300.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "west115_yy"
fdnum = 3
ifnum = 1
plnum = 1
polariza = "YY"
beam_id = "B03"
db_stream_name = "xffts-board7"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = false
beam_fit_use = false

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = -300.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "south115_xx"
fdnum = 4
ifnum = 1
plnum = 0
polariza = "XX"
beam_id = "B04"
db_stream_name = "xffts-board8"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = -300.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "south115_yy"
fdnum = 4
ifnum = 1
plnum = 1
polariza = "YY"
beam_id = "B04"
db_stream_name = "xffts-board9"
enabled = true
use_for_convert = true
use_for_sunscan = true
use_for_fit = false
beam_fit_use = false

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = -300.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0
```

---

## 6. config をまず検査する

```bash
check_spectrometer_config multibeam_5beam_example.toml
```

CSV も欲しければ:

```bash
check_spectrometer_config multibeam_5beam_example.toml \
  --out-csv config_check.csv
```

### ここで確認すること
- stream 名が期待通りか
- `beam_id` が XX/YY で正しく共有されているか
- `restfreq_hz` が有限値か
- `beam_fit_use` の primary 選択が意図通りか
- warning が出ていないか

---

## 7. extract を実行する

### 7.1 全 stream を解析

```bash
sunscan_extract_multibeam RAWDATA \
  --spectrometer-config multibeam_5beam_example.toml \
  --outdir out_extract
```

### 7.2 XX だけを解析

```bash
sunscan_extract_multibeam RAWDATA \
  --spectrometer-config multibeam_5beam_example.toml \
  --stream-name center230_xx \
  --stream-name east115_xx \
  --stream-name north115_xx \
  --stream-name west115_xx \
  --stream-name south115_xx \
  --outdir out_extract_xx
```

### 7.3 CLI で telescope / planet / offset を一時上書き

```bash
sunscan_extract_multibeam RAWDATA \
  --spectrometer-config multibeam_5beam_example.toml \
  --telescope OMU1P85M \
  --tel-loaddata OMU1p85m \
  --planet sun \
  --spectrometer-time-offset-sec -0.050 \
  --encoder-el-time-offset-sec 0.020 \
  --outdir out_extract_override
```

### 7.4 extract 後に確認するもの
- `sunscan_multibeam_scan_summary_<tag>.csv`
- `sunscan_multibeam_manifest_<tag>.csv`
- `spectrometer_stream_table_<tag>.csv`
- `analysis_config_snapshot_<tag>.json`
- `per_stream/<stream_name>/...`

---

## 8. disabled stream を残したまま一時的に外す

たとえば `north115_yy` を一時的に外すなら、削除せずこうします。

```toml
[[spectrometers]]
name = "north115_yy"
...
enabled = true
use_for_convert = true
use_for_sunscan = false
use_for_fit = false
beam_fit_use = false
```

### 意味
- converter では残す
- sunscan extract では通常は自動スキップ
- fit からも外す

### extract の current behavior
- 通常実行では manifest に `status="skipped"` として残る
- 明示 `--stream-name north115_yy` を付ければ one-off override 可能

---

## 9. YY は fit には使わず、extract には残したい

```toml
use_for_convert = true
use_for_sunscan = true
use_for_fit = false
beam_fit_use = false
```

こうしておくと

- extract では XX/YY とも individual summary が残る
- fit は XX のみを対象にできる

最初の実データでは、この運用が安全です。

---

## 10. fit を実行する

### 10.1 最初は `both` を推奨

```bash
sunscan_fit_multibeam \
  out_extract/sunscan_multibeam_scan_summary_RUN.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --center-beam-id B00 \
  --model both \
  --outdir out_fit
```

### 10.2 virtual center だけで試す

```bash
sunscan_fit_multibeam \
  out_extract/sunscan_multibeam_scan_summary_RUN.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --model virtual_center \
  --outdir out_fit_virtual
```

### 10.3 fit 対象 stream を明示指定する

```bash
sunscan_fit_multibeam \
  out_extract/sunscan_multibeam_scan_summary_RUN.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --center-beam-id B00 \
  --fit-stream-name center230_xx \
  --fit-stream-name east115_xx \
  --fit-stream-name north115_xx \
  --fit-stream-name west115_xx \
  --fit-stream-name south115_xx \
  --outdir out_fit_xx
```

### 10.4 fit 後に確認するもの
- `fit_summary.txt`
- `beam_fit_results_<model>.csv`
- `beam_fit_residuals_<model>.csv`
- `beam_fit_rejected_<model>.csv`
- `beam_model_<model>.toml`

---

## 11. pseudo multi-beam を作る

### 11.1 単純な dry-run

```bash
sunscan_make_pseudo_multibeam \
  out_single/sun_scan_summary_single115_xx.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --outdir out_pseudo
```

### 11.2 ノイズを加える

```bash
sunscan_make_pseudo_multibeam \
  out_single/sun_scan_summary_single115_xx.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --noise-arcsec 5.0 \
  --seed 1234 \
  --outdir out_pseudo_noise
```

### 11.3 representative EL を人工的に複数入れる

```bash
sunscan_make_pseudo_multibeam \
  out_single/sun_scan_summary_single115_xx.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --rep-el-deg 25 \
  --rep-el-deg 45 \
  --rep-el-deg 65 \
  --outdir out_pseudo_multi_el
```

### 11.4 注意
`sunscan_multibeam pseudo` では current 実装上、`--rep-el-deg` の forward が完全ではありません。  
この用途では **standalone の `sunscan_make_pseudo_multibeam` を直接使う**方が安全です。

---

## 12. pseudo から fit まで通す

```bash
sunscan_make_pseudo_multibeam \
  out_single/sun_scan_summary_single115_xx.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --rep-el-deg 25 \
  --rep-el-deg 45 \
  --rep-el-deg 65 \
  --outdir out_pseudo

sunscan_fit_multibeam \
  out_pseudo/pseudo_multibeam_summary_pseudo_multibeam.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --center-beam-id B00 \
  --model both \
  --outdir out_pseudo_fit
```

これで
- config の beam 配置
- fit の流れ
- output file 群
を、実 RawData が無い段階で検証できます。

---

## 13. `sunscan_multibeam` wrapper を使う例

### 13.1 extract
```bash
sunscan_multibeam extract RAWDATA \
  --spectrometer-config multibeam_5beam_example.toml \
  --outdir out_extract
```

### 13.2 fit
```bash
sunscan_multibeam fit \
  out_extract/sunscan_multibeam_scan_summary_RUN.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --center-beam-id B00 \
  --model both \
  --outdir out_fit
```

### 13.3 check-config
```bash
sunscan_multibeam check-config multibeam_5beam_example.toml
```

### 13.4 pseudo
```bash
sunscan_multibeam pseudo \
  out_single/sun_scan_summary_single115_xx.csv \
  --spectrometer-config multibeam_5beam_example.toml \
  --outdir out_pseudo
```

### 13.5 wrapper と standalone の使い分け
- 普段は wrapper でもよい
- 細かい挙動を制御したいときは standalone が明確

---

## 14. converter と共用するときの考え方

## 14.1 考え方
converter と sunscan は別機能ですが、同じ spectrometer config TOML を共有できます。  
共有の中核は

- `name`
- `beam_id`
- `db_stream_name`
- `db_table_name`
- `fdnum / ifnum / plnum / polariza`
- `restfreq_hz`
- beam offset / rotation
- usage flags

です。

## 14.2 共用時の実務上のおすすめ
- stream 定義は 1 つの TOML で一元管理
- sunscan 固有の解析条件は CLI で与える
- beam geometry fit 後は `beam_model_*.toml` を候補として converter 側へ反映
- 一時的に使えない stream は削除せず usage flags で管理

---

## 15. よくあるトラブルと対策

## 15.1 `stream_name not found`
原因:
- `--stream-name` が `name` と合っていない

対策:
- `check_spectrometer_config` で stream table を確認

## 15.2 summary にはあるのに fit で使われない
原因候補:
- `use_for_fit = false`
- primary stream 解決で他 stream が代表になっている
- `--fit-stream-name` で絞っている

対策:
- `check_spectrometer_config`
- `fit_summary.txt`
- `selected_rows_for_fit.csv`
を確認

## 15.3 expected beam offsets が全部 0
原因:
- nominal beam offset が未設定

対策:
- pseudo dry-run 用でもよいので nominal offsets を入れる

## 15.4 extract で stream が消えたように見える
原因:
- `use_for_sunscan = false`

対策:
- manifest を見る
- `status="skipped"` と `skip_reason` を確認

## 15.5 time offset の解釈が分からなくなる
対策:
- 新方式では `encoder_shift_sec = 0.0` を基本にする
- summary / snapshot に残る metadata を確認する

---

## 16. 初回運用のおすすめセット

### 16.1 最初の最小流れ
1. `check_spectrometer_config`
2. `sunscan_singlebeam`
3. `sunscan_extract_multibeam`
4. `sunscan_fit_multibeam --model both`

### 16.2 初回設定の考え方
- XX を primary にする
- YY は extract だけ残してもよい
- `encoder_shift_sec = 0.0`
- Az / El 個別補正を分ける
- `rotation_mode = "elevation"` は dry-run 検証に有用
- まず single でうまくいく条件を使う

---

## 17. まとめ

この cookbook で重要なのは次です。

- single で成功してから multi に進む
- usage flags で stream を削除せず管理する
- fit は最初 `both` で比較する
- pseudo は実 multi-beam RawData が無い段階でも非常に有用
- converter と sunscan は別機能だが、stream / beam / restfreq / beam model は shared config として整理すると運用しやすい
