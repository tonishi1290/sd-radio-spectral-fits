# cookbook

この cookbook は、現行 converter / singlebeam / multibeam を実務で使うための具体例集です。

---

## 1. 変数の置き方

```bash
RAW="/path/to/RawData_20260316_120000"
RUN="$(basename "$RAW")"
CFG_SINGLE="./omu_singlebeam_example.toml"
CFG_MULTI="./omu_multibeam_example.toml"
```

---

## 2. converter: 単一 stream を config なしで変換

```bash
python necst_v4_sdfits_converter.py \
  "$RAW" \
  --spectral xffts-board1 \
  --object ORION_KL \
  --project TEST \
  --observer OHNISHI \
  --out "./out/${RUN}_single_legacy.fits"
```

使いどころ:

- まず実 table 名が合っているか確認したい
- TOML を使う前に単一 stream の生存確認をしたい

---

## 3. converter: spectrometer config を使って変換

```bash
python necst_v4_sdfits_converter.py \
  "$RAW" \
  --spectrometer-config "$CFG_MULTI" \
  --observer "Your Name" \
  --project "MULTIBEAM_TEST" \
  --object "SUN" \
  --out "./out/${RUN}_multistream.fits"
```

---

## 4. converter: table 名が規則と合わない場合

TOML 側で `db_table_name` を直接指定します。

```toml
[[spectrometers]]
name = "single115_xx"
db_table_name = "necst-OMU1P85M-data-spectral-xffts-board1"
```

こうすると `db_stream_name` に依存せず、実 table 名へ直接アクセスします。

---

## 5. converter: position-switch 的に `cmd` を使って Az/El 由来 RA/DEC を作る

```bash
python necst_v4_sdfits_converter.py \
  "$RAW" \
  --spectrometer-config "$CFG_SINGLE" \
  --radec-method azel \
  --radec-azel-source cmd \
  --out "./out/${RUN}_cmd_radec.fits"
```

解釈:

- WCS ではなく Az/El 由来の RA/DEC を使う
- Az/El source は `cmd`
- つまり、altaz/cmd 側から補正量を引いた boresight に beam offset を適用したものを使う

---

## 6. singlebeam: config を使って 1 stream を解析

```bash
python -m multibeam_beam_measurement.sunscan_singlebeam \
  "$RAW" \
  --spectrometer-config "$CFG_SINGLE" \
  --stream-name single115_xx \
  --outdir "./out/singlebeam"
```

---

## 7. singlebeam: 時刻補正を CLI で上書き

```bash
python -m multibeam_beam_measurement.sunscan_singlebeam \
  "$RAW" \
  --spectrometer-config "$CFG_SINGLE" \
  --stream-name single115_xx \
  --spectrometer-time-offset-sec -0.050 \
  --encoder-el-time-offset-sec 0.020 \
  --outdir "./out/singlebeam_override"
```

使いどころ:

- config は固定しつつ、時刻補正だけ探索したい

---

## 8. singlebeam: altaz を使って見る

```bash
python -m multibeam_beam_measurement.sunscan_singlebeam \
  "$RAW" \
  --spectrometer-config "$CFG_SINGLE" \
  --stream-name single115_xx \
  --azel-source altaz \
  --altaz-apply minus \
  --outdir "./out/singlebeam_altaz"
```

使いどころ:

- encoder と altaz で結果がどう違うか見たい
- 符号の整合を確認したい

---

## 9. multibeam extract: 全 stream をまとめて抽出

```bash
python -m multibeam_beam_measurement.sunscan_extract_multibeam \
  "$RAW" \
  --spectrometer-config "$CFG_MULTI" \
  --outdir "./out/multi_extract_all"
```

ここで確認すべき出力:

- `sunscan_multibeam_scan_summary_*.csv`
- `sunscan_multibeam_manifest_*.csv`
- `spectrometer_stream_table_*.csv`
- `analysis_config_snapshot_*.json`

---

## 10. multibeam extract: XX だけを抽出

```bash
python -m multibeam_beam_measurement.sunscan_extract_multibeam \
  "$RAW" \
  --spectrometer-config "$CFG_MULTI" \
  --stream-name center230_xx \
  --stream-name east115_xx \
  --stream-name north115_xx \
  --stream-name west115_xx \
  --stream-name south115_xx \
  --outdir "./out/multi_extract_xx_only"
```

使いどころ:

- YY を切った dry-run
- まず 1 偏波だけで幾何を見たい

---

## 11. multibeam fit: `both` でモデル比較

```bash
python -m multibeam_beam_measurement.sunscan_fit_multibeam \
  "./out/multi_extract_all/sunscan_multibeam_scan_summary_${RUN}.csv" \
  --spectrometer-config "$CFG_MULTI" \
  --center-beam-id B00 \
  --model both \
  --outdir "./out/multi_fit_both"
```

確認すべき出力:

- `beam_fit_results_center_beam.csv`
- `beam_fit_results_virtual_center.csv`
- `fit_summary.txt`
- `beam_model_center_beam.toml`
- `beam_model_virtual_center.toml`

---

## 12. multibeam fit: YY だけを明示選択

```bash
python -m multibeam_beam_measurement.sunscan_fit_multibeam \
  "./out/multi_extract_all/sunscan_multibeam_scan_summary_${RUN}.csv" \
  --spectrometer-config "$CFG_MULTI" \
  --center-beam-id B00 \
  --model both \
  --fit-stream-name center230_yy \
  --fit-stream-name east115_yy \
  --fit-stream-name north115_yy \
  --fit-stream-name west115_yy \
  --fit-stream-name south115_yy \
  --outdir "./out/multi_fit_yy"
```

---

## 13. multibeam fit: outlier に厳しめ

```bash
python -m multibeam_beam_measurement.sunscan_fit_multibeam \
  "./out/multi_extract_all/sunscan_multibeam_scan_summary_${RUN}.csv" \
  --spectrometer-config "$CFG_MULTI" \
  --center-beam-id B00 \
  --model both \
  --sigma-clip 3.5 \
  --clip-iters 3 \
  --min-points-per-beam 3 \
  --min-scans-per-beam 2 \
  --outdir "./out/multi_fit_strict"
```

---

## 14. pseudo multibeam: singlebeam summary から疑似 multi を作る

```bash
python -m multibeam_beam_measurement.synthetic_multibeam \
  ./out/singlebeam/sun_scan_summary_single115_xx.csv \
  --spectrometer-config "$CFG_MULTI" \
  --outdir "./out/pseudo_multi"
```

使いどころ:

- 実 multi-beam データがまだ無い
- config と fit の流れだけ確認したい

---

## 15. pseudo multibeam: jitter を加える

```bash
python -m multibeam_beam_measurement.synthetic_multibeam \
  ./out/singlebeam/sun_scan_summary_single115_xx.csv \
  --spectrometer-config "$CFG_MULTI" \
  --noise-arcsec 5.0 \
  --seed 42 \
  --outdir "./out/pseudo_multi_noisy"
```

---

## 16. check-config: primary stream 解決を確認

```bash
python -m multibeam_beam_measurement.check_spectrometer_config \
  "$CFG_MULTI" \
  --out-csv ./out/stream_table.csv
```

この段階で確認すること:

- beam_id が重複していても意味的に正しいか
- `beam_fit_use` が 1 beam に 1 本だけ立っているか
- nominal offsets が 0 ばかりではないか

---

## 17. wrapper を使う

### extract

```bash
python -m multibeam_beam_measurement.sunscan_multibeam \
  extract \
  "$RAW" \
  --spectrometer-config "$CFG_MULTI" \
  --outdir ./out/wrapper_extract
```

### fit

```bash
python -m multibeam_beam_measurement.sunscan_multibeam \
  fit \
  ./out/wrapper_extract/sunscan_multibeam_scan_summary_${RUN}.csv \
  --spectrometer-config "$CFG_MULTI" \
  --center-beam-id B00 \
  --model both \
  --outdir ./out/wrapper_fit
```

### pseudo

```bash
python -m multibeam_beam_measurement.sunscan_multibeam \
  pseudo \
  ./out/singlebeam/sun_scan_summary_single115_xx.csv \
  --spectrometer-config "$CFG_MULTI" \
  --outdir ./out/wrapper_pseudo
```

### check-config

```bash
python -m multibeam_beam_measurement.sunscan_multibeam \
  check-config \
  "$CFG_MULTI" \
  --out-csv ./out/wrapper_stream_table.csv
```

---

## 18. 問題のある stream を一時的に無効化する

例えば `north115_yy` をしばらく外したいなら、削除せずに次のようにします。

```toml
[[spectrometers]]
name = "north115_yy"
beam_id = "B02"
polariza = "YY"
enabled = true
use_for_convert = false
use_for_sunscan = false
use_for_fit = false
beam_fit_use = false
```

これで定義は残したまま、通常の変換・抽出・fit から外せます。

---

## 19. 12CO singlebeam 用 config 例

```toml
[global]
db_namespace = "necst"
telescope = "OMU1P85M"
tel_loaddata = "OMU1p85m"
planet = "sun"

spectrometer_time_offset_sec = -0.070
encoder_shift_sec = 0.000
encoder_az_time_offset_sec = 0.000
encoder_el_time_offset_sec = 0.040

[[spectrometers]]
name = "single115_xx"
fdnum = 0
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B00"
frontend = "RX115"
backend = "XFFTS"
sampler = "XFFTS0"
db_stream_name = "xffts-board1"
enabled = true
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
nchan = 32768
definition_mode = "band_start_stop"
band_start_hz = 0.0
band_stop_hz = 2.5e9
channel_origin = "center"
reverse = false
ctype1 = "FREQ"
cunit1 = "Hz"
specsys = "TOPOCENT"
veldef = "RADIO"
store_freq_column = "auto"
restfreq_hz = 115.2712018e9

[spectrometers.local_oscillators]
lo1_hz = 109.8e9
lo2_hz = 4.0e9
sb1 = "USB"
sb2 = "USB"
obsfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "none"
reference_angle_deg = 0.0
rotation_sign = 1.0
dewar_angle_deg = 0.0
```

---

## 20. center 230 GHz + 周辺 115 GHz の multi-beam 例

```toml
[global]
db_namespace = "necst"
telescope = "OMU1P85M"
tel_loaddata = "OMU1p85m"
planet = "sun"

spectrometer_time_offset_sec = -0.070
encoder_shift_sec = 0.000
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
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 230.538000e9
nchan = 32768
definition_mode = "band_start_stop"
band_start_hz = 0.0
band_stop_hz = 2.5e9
channel_origin = "center"
reverse = false
ctype1 = "FREQ"
cunit1 = "Hz"
specsys = "TOPOCENT"
veldef = "RADIO"
store_freq_column = "auto"

[spectrometers.local_oscillators]
lo1_hz = 225.0e9
lo2_hz = 4.0e9
sb1 = "USB"
sb2 = "USB"
obsfreq_hz = 230.538000e9

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
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 115.2712018e9
nchan = 32768
definition_mode = "band_start_stop"
band_start_hz = 0.0
band_stop_hz = 2.5e9
channel_origin = "center"
reverse = false
ctype1 = "FREQ"
cunit1 = "Hz"
specsys = "TOPOCENT"
veldef = "RADIO"
store_freq_column = "auto"

[spectrometers.local_oscillators]
lo1_hz = 109.8e9
lo2_hz = 4.0e9
sb1 = "USB"
sb2 = "USB"
obsfreq_hz = 115.2712018e9

[spectrometers.beam]
az_offset_arcsec = 300.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 45.0
rotation_sign = 1.0
dewar_angle_deg = 0.0
```

実運用ではこれを north / west / south と XX / YY 分だけ展開してください。
