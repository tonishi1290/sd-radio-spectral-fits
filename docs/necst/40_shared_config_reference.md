# shared spectrometer config TOML 完全リファレンス

このファイルは、converter と sunscan が共有できる spectrometer config TOML の説明書です。

重要な注意:

- **shared** ではあるが、全キーが全コマンドで同じ意味で効くわけではありません
- 現行コードで実際に使われる範囲を明記します

---

## 1. 基本構造

```toml
[global]
...

[[spectrometers]]
name = "..."
fdnum = 0
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B00"
...

[spectrometers.frequency_axis]
...

[spectrometers.local_oscillators]
...

[spectrometers.override]
...

[spectrometers.beam]
...
```

---

## 2. `[global]`

### `db_namespace`

DB table 名の namespace。

- converter: 使用
- singlebeam: 使用
- multibeam extract: 使用
- fit / pseudo / check-config: 直接は不要

### `telescope`

DB table 名の telescope 名。

- converter: 使用
- singlebeam: 使用
- multibeam extract: 使用

### `tel_loaddata`

legacy loader 互換の telescope 名。

- converter: 実質未使用
- singlebeam: 使用
- multibeam extract: 使用

### `planet`

対象天体名。

- converter: 直接は使わない
- singlebeam: 使用
- multibeam extract: 使用

### `spectrometer_time_offset_sec`

spectrometer 側時刻補正 [s]。

- converter: 現行 global semantics としては中心ではない。少なくとも singlebeam/extract のような対称な global 処理としては扱われない
- singlebeam: 使用
- multibeam extract: 使用

### `encoder_shift_sec`

encoder 共通時刻シフト [s]。

- converter: 使用
- singlebeam: 使用
- multibeam extract: 使用

### `encoder_az_time_offset_sec`

encoder Az 系列の追加 offset [s]。

- converter: 現行 global semantics として singlebeam/extract ほど対称ではない
- singlebeam: 使用
- multibeam extract: 使用

### `encoder_el_time_offset_sec`

encoder El 系列の追加 offset [s]。

- converter: 現行 global semantics として singlebeam/extract ほど対称ではない
- singlebeam: 使用
- multibeam extract: 使用

### `encoder_table`

encoder table 名を固定します。

- converter: 使用
- singlebeam / extract: 通常は独自読み出し系で table 名を構成

### `altaz_table`

altaz table 名を固定します。

- converter: 使用

### `weather_table`

weather table 名を固定します。

- converter: 使用

### `output_layout`

converter の出力 layout。

- converter: 使用
- sunscan: 未使用

現行 converter では `merged` / `merged_time` のみが前提です。

### `time_sort`

merged output で時間ソートを要求する flag。

- converter: 使用
- sunscan: 未使用

### `channel_slice`

global な channel 範囲指定。

- converter: 使用
- singlebeam / extract: 直接は使わない

---

## 3. `[[spectrometers]]` scalar fields

### `name`

stream の内部名。必須。

- converter: 使用
- singlebeam: `--stream-name` で参照
- extract: `--stream-name` で参照
- fit: primary stream 解決で参照
- pseudo / check-config: 使用

### `fdnum`

FDNUM。

### `ifnum`

IFNUM。

### `plnum`

PLNUM。

### `polariza`

偏波。
許容される主な値は `RR LL RL LR XX YY XY YX` です。

### `beam_id`

beam 識別子。
同一 beam の XX/YY をまとめる鍵です。

### `frontend`

frontend 名。

### `backend`

backend 名。

### `sampler`

sampler 名。

### `db_stream_name`

自動 table 名構築に使う stream 名。

### `db_table_name`

spectral table 名の完全指定。
これがあれば converter の table 解決で最優先です。

### `enabled`

stream 総合有効 / 無効フラグ。

- converter: 現行では stream usage semantics としては未使用寄り
- singlebeam: `config_io` を通して stream 選択補助に反映されうる
- extract: 使用
- fit: 使用
- pseudo / check-config: 使用

### `use_for_convert`

converter で使うかの意図を表すフラグ。

- converter: 現行では自動スキップ用 semantics としては未統合
- extract: stream table / config validation では保持される
- fit: 直接未使用

### `use_for_sunscan`

extract の通常対象かどうか。

- extract: 使用
- fit: 直接未使用

### `use_for_fit`

fit の既定対象かどうか。

- fit: 使用
- check-config: 使用

### `beam_fit_use`

同一 beam の複数 stream から primary stream を 1 本選ぶための印です。

- fit: 使用
- check-config: 使用

---

## 4. `[spectrometers.frequency_axis]`

converter が周波数軸 WCS を作るための主表です。

### `nchan`

チャンネル数。

### `definition_mode`

周波数軸の定義モード。
現行例では `band_start_stop` をよく使います。

### `band_start_hz`

帯域開始 [Hz]。

### `band_stop_hz`

帯域終了 [Hz]。

### `channel_origin`

基準チャンネルの扱い。
例: `center`

### `reverse`

チャンネル軸を逆順にするか。

### `ctype1`

FITS `CTYPE1`。
通常 `FREQ`。

### `cunit1`

FITS `CUNIT1`。
通常 `Hz`。

### `specsys`

スペクトル系。

### `veldef`

速度定義。

### `store_freq_column`

周波数列を保持するかの方針。

### `restfreq_hz`

RESTFREQ [Hz]。
converter にも multi-beam 側にも重要です。

### `channel_slice`

stream 個別の channel slice。
converter 側で使われます。

---

## 5. `[spectrometers.local_oscillators]`

### `lo1_hz`, `lo2_hz`, `lo3_hz`

LO 周波数 [Hz]。

### `sb1`, `sb2`, `sb3`

各段の sideband。

### `obsfreq_hz`

観測周波数 [Hz]。

### `imagfreq_hz`

イメージ周波数 [Hz]。

### `sideband`

総合 sideband 記録用。

---

## 6. `[spectrometers.override]`

converter / sunscan の将来拡張や stream 個別 override 用の領域です。
現行コードでは強い一般仕様はなく、必要なメタデータをぶら下げるための保留領域と考えるのが安全です。

---

## 7. `[spectrometers.beam]`

### `az_offset_arcsec`

beam の nominal Az 方向 offset [arcsec]。

### `el_offset_arcsec`

beam の nominal El 方向 offset [arcsec]。

### `rotation_mode`

beam rotation モード。

- `none`
- `elevation`

### `reference_angle_deg`

回転モデルの基準角 [deg]。

### `rotation_sign`

回転の符号。通常 `-1`, `0`, `+1` のいずれか。

### `dewar_angle_deg`

dewar 固有の位相ずれ [deg]。

### `beam_model_version`

beam model の来歴用文字列。
multi-beam fit が書き戻す TOML では、ここに `sunscan_multibeam_<model>` のような値が入ります。

---

## 8. support の考え方

### converter を主に見るとき

重要なのは

- DB 名解決
- WCS
- LO
- beam model

です。

### singlebeam / extract を主に見るとき

重要なのは

- `global` の入力・時刻補正
- stream metadata
- `restfreq_hz`
- `beam_id`

です。

### fit を主に見るとき

重要なのは

- `beam_id`
- `use_for_fit`
- `beam_fit_use`
- nominal beam offsets / rotation_mode

です。

---

## 9. 典型的な記述パターン

### 9.1 single-beam 用

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
db_stream_name = "xffts-board1"

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
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "none"
reference_angle_deg = 0.0
rotation_sign = 1.0
dewar_angle_deg = 0.0
```

### 9.2 multi-beam 用

同一 `beam_id` の XX / YY を持つ場合は、fit の代表 stream を明示してください。

```toml
[[spectrometers]]
name = "east115_xx"
beam_id = "B01"
polariza = "XX"
enabled = true
use_for_sunscan = true
use_for_fit = true
beam_fit_use = true

[[spectrometers]]
name = "east115_yy"
beam_id = "B01"
polariza = "YY"
enabled = true
use_for_sunscan = true
use_for_fit = false
beam_fit_use = false
```

---

## 10. 安全な運用指針

1. 実 table 名が怪しいときは `db_table_name` を使う
2. XX / YY があるなら `beam_fit_use` を明示する
3. pseudo dry-run を使うなら nominal beam offsets を 0 のままにしない
4. rotation を評価したいなら `rotation_mode='elevation'` を設定する
