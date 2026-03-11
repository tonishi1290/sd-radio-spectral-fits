# Single-beam 解析・SDFITS 変換 説明書

## 1. 本書の目的

本書は、single-beam の太陽スキャン解析と SDFITS 変換を、次の 2 つのコマンドで実行するための利用者向け説明書です。

- `sunscan_singlebeam`
- `necst_v4_sdfits_converter`

本書では、これら 2 つが **`pyproject.toml` の console_scripts に登録されている**前提で説明します。したがって、実行例はすべて

```bash
sunscan_singlebeam ...
necst_v4_sdfits_converter ...
```

の形式で記載します。

`sun_scan_v4_v5.py` は **legacy / 検証用** と位置づけ、通常運用では `sunscan_singlebeam` を使用します。

---

## 2. 全体の位置づけ

single-beam の標準的な作業順は次のとおりです。

1. `sunscan_singlebeam` で太陽スキャンを解析する
2. 解析結果を見て、観測が妥当かを確認する
3. 必要に応じて `necst_v4_sdfits_converter` で RawData を SDFITS へ変換する

役割分担は明確です。

- `sunscan_singlebeam`
  - 太陽 limb fitting による pointing / beam 評価
  - summary CSV / derivative-fit PNG / 端末表示の生成
- `necst_v4_sdfits_converter`
  - NECST RawData の分光 DB を SDFITS に変換
  - spectrometer/beam/LO/WCS を TOML で与える

---

## 3. コマンドの前提

### 3.1 `sunscan_singlebeam`

- RawData ディレクトリを入力とする
- single-beam 用の標準解析コマンド
- `sun_scan_v4_v5.py` 相当の解析を package 側実装で実行する
- 出力は single-beam 互換の summary CSV / derivative-fit PNG / 端末 summary

### 3.2 `necst_v4_sdfits_converter`

- RawData ディレクトリを入力とする
- spectrometer 設定 TOML を与えると、stream / beam / 周波数軸 / LO / 偏波を明示的に定義できる
- TOML を与えない場合は legacy single-stream モードでも動作する

---

## 4. `sunscan_singlebeam` の使い方

## 4.1 基本形

```bash
sunscan_singlebeam RAWDATA \
  --telescope OMU1P85M \
  --tel-loaddata OMU1p85m \
  --planet sun \
  --spectral-name xffts-board1 \
  --outdir out_single
```

### 主な入出力

- 入力
  - `RAWDATA`: NECST RawData ディレクトリ
  - `--spectral-name`: 解析対象の分光 stream 名
- 出力
  - `sun_scan_summary_<tag>.csv`
  - `sun_scan_derivative_fits_<tag>_pXXX.png`
  - 端末への summary 表示

`<tag>` は通常、RawData ディレクトリ名をもとに決まります。

### DB 名解決に必要な 3 つのパラメータ

`sunscan_singlebeam` では、RawData 内の DB 名解決と天体位置計算のために、次の 3 つの値を内部で使います。

- `--telescope`
  - DB table 名の解決に使う望遠鏡名
  - 例: `OMU1P85M`
- `--tel-loaddata`
  - `loaddb(...)` に渡す望遠鏡名
  - 例: `OMU1p85m`
- `--planet`
  - `astropy` の `get_body(...)` に渡す天体名
  - 例: `sun`

現在の既定値は

- `--telescope OMU1P85M`
- `--tel-loaddata OMU1p85m`
- `--planet sun`

です。OMU 1.85 m の太陽スキャンでは通常この既定値のままで構いませんが、**望遠鏡名や対象天体が異なる場合はここを必ず合わせてください。**

---

## 4.2 実運用でよく使う例

```bash
sunscan_singlebeam "$RAW" \
  --telescope OMU1P85M \
  --tel-loaddata OMU1p85m \
  --planet sun \
  --spectral-name xffts-board1 \
  --outdir out_single \
  --azel-source encoder \
  --altaz-apply none \
  --profile-xlim-deg 1.0 \
  --ripple-preset auto \
  --edge-fit-win-deg 0.15 \
  --edge-fit-threshold 0.20 \
  --hpbw-init-arcsec 324.0
```

---

## 4.3 よく使うオプション

### 座標・時系列まわり

- `--telescope`
  - DB table 名の解決に使う望遠鏡名
  - 例: `OMU1P85M`
- `--tel-loaddata`
  - `loaddb(...)` に渡す望遠鏡名
  - 例: `OMU1p85m`
- `--planet`
  - 天体位置計算に使う対象名
  - 例: `sun`
- `--spectral-name`
  - 分光 stream 名
- `--azel-source {encoder,altaz}`
  - Az/El の基準系列
- `--altaz-apply {none,minus,plus}`
  - altaz 補正の適用法
- `--encoder-shift-sec`
  - encoder 時系列の時間シフト
- `--encoder-vavg-sec`
  - encoder 速度の平滑化窓

### 校正まわり

- `--no-chopper-wheel`
  - chopper-wheel 校正を無効化
- `--tamb-k`
  - 外気温 [K] を固定指定
- `--chopper-win-sec`
  - HOT/OFF 推定窓
- `--chopper-stat {median,mean}`
  - HOT/OFF 推定の代表値

### ripple 除去

- `--ripple-no-remove`
- `--ripple-preset {auto,safe,normal,strong}`
- `--ripple-model {auto,add,mul}`
- `--ripple-target-hz`
- `--ripple-search-hz`
- `--ripple-bw-hz`
- `--ripple-max-harm`
- `--ripple-order`
- `--ripple-notch-pass`
- `--ripple-trend-win-sec`
- `--ripple-resample-dt-sec`
- `--ripple-eval-band-hz`

通常は `--ripple-preset auto` を起点にし、必要になったときだけ個別パラメータを触るのが安全です。

### edge fit

- `--no-edge-fit`
- `--edge-fit-win-deg`
- `--edge-fit-threshold`
- `--hpbw-init-arcsec`
- `--edge-fit-plot-max-scans`

### trim / scan 選別

- `--no-trim-scan`
- `--trim-vfrac`
- `--trim-vmin`
- `--trim-gap`
- `--trim-min-samples`
- `--trim-dominant-axis` / `--trim-no-dominant-axis`
- `--trim-axis-ratio-min`
- `--trim-vpercentile`
- `--trim-no-steady-scan`

---

## 4.4 出力の見方

端末 summary では主に次を見ます。

- `rep_Az`, `rep_El`
  - その scan の代表 Az/El
- `center_az`, `center_el`
  - 太陽中心に対する推定 offset [arcsec]
- `HPBW_az`, `HPBW_el`
  - limb fit から得た見かけの HPBW [arcsec]
- `fitAZ`, `fitEL`
  - AZ / EL fit の成否

CSV では、少なくとも次の列を確認してください。

- `scan_id`
- `rep_az_deg`
- `rep_el_deg`
- `center_az_deg`
- `center_el_deg`
- `hpbw_az_arcsec`
- `hpbw_el_arcsec`
- `fit_ok_az`
- `fit_ok_el`

---

## 4.5 `sun_scan_v4_v5.py` との関係

- `sun_scan_v4_v5.py` は legacy / 確認用です
- 通常運用では `sunscan_singlebeam` を使います
- 必要なときだけ、`sun_scan_v4_v5.py` を用いて比較確認します

すなわち、**single-beam 解析の実質的な標準コマンドは `sunscan_singlebeam`** と考えてください。

---

## 5. `necst_v4_sdfits_converter` の使い方

## 5.1 基本形

```bash
necst_v4_sdfits_converter RAWDATA \
  --spectrometer-config beams.toml \
  --out output.fits
```

RawData から SDFITS を作成します。`--spectrometer-config` を与えると、stream ごとに beam / 偏波 / 周波数軸 / LO / DB 名を明示的に定義できます。

---

## 5.2 最小の legacy single-stream 例

TOML を使わずに legacy モードで動かす場合です。

```bash
necst_v4_sdfits_converter RAWDATA \
  --spectral xffts-board1 \
  --out output.fits \
  --restfreq-hz 115271201800.0
```

ただし、今後の運用では **TOML を与える方法を標準**とすることを勧めます。

---

## 5.3 実運用でよく使う例

```bash
necst_v4_sdfits_converter "$RAW" \
  --spectrometer-config beams.toml \
  --out out_singlebeam.fits \
  --object OrionKL \
  --project ProjectID \
  --observer YourName
```

---

## 5.4 主な CLI オプション

### 入力・設定

- `--spectrometer-config`
  - spectrometer / beam / WCS / LO を記述した TOML
- `--strict-config`
  - stream 不整合をエラー扱いにする
- `--spectral`
  - legacy single-stream 名
- `--telescope`
  - テーブル名解決に使う望遠鏡名
- `--db-namespace`
  - DB namespace/prefix

### 出力情報

- `--out`
  - 出力 FITS 名
- `--object`
  - FITS の `OBJECT`
- `--project`
  - FITS の `PROJID`
- `--observer`
  - FITS の `OBSERVER`

### 観測地

- `--site-lat`
- `--site-lon`
- `--site-elev`

### legacy single-stream 用周波数軸

- `--nchan`
- `--if0-ghz`
- `--if1-ghz`
- `--lo1-ghz`
- `--lo2-ghz`
- `--restfreq-hz`

### そのほか

- `--channel-slice`
- `--encoder-time-col`
- `--altaz-time-col`
- `--interp-extrap {nan,hold}`
- `--use-modes`
- `--include-unknown`
- `--wcs-table`
- `--radec-method {wcs_first,azel}`
- `--radec-azel-source {beam,true,encoder,altaz,cmd}`
- `--refraction {on,off}`
- `--met-source {auto,weather,spectral,fallback}`
- `--weather-table`
- `--collapse`

---

## 6. converter に渡す `beams.toml` の丁寧な説明

## 6.1 基本方針

single-beam でも、将来の拡張や再現性のため、**TOML で stream / beam / WCS / LO を明示**することを推奨します。

single-beam なら `[[spectrometers]]` は 1 ブロックで十分です。

---

## 6.2 最小例（single-beam）

```toml
schema_version = 1
config_name = "singlebeam_example"
config_description = "Single-beam example for necst_v4_sdfits_converter"

[global]
db_namespace = "necst"
telescope = "OMU1P85M"
output_layout = "merged"
time_sort = true

[[spectrometers]]
name = "xffts-board1"
fdnum = 0
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B00"
db_stream_name = "xffts-board1"
frontend = "RX"
backend = "XFFTS"
sampler = "xffts-board1"

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "none"
reference_angle_deg = 0.0
rotation_sign = 1.0
dewar_angle_deg = 0.0
beam_model_version = "singlebeam_v1"

[spectrometers.frequency_axis]
nchan = 32768
definition_mode = "first_center_and_delta"
first_channel_center_hz = 0.0
channel_spacing_hz = 76296.27368999299
crpix1 = 1.0
restfreq_hz = 115271201800.0
ctype1 = "FREQ"
cunit1 = "Hz"
specsys = "TOPOCENT"
veldef = "RADIO"
store_freq_column = false

[spectrometers.local_oscillators]
lo1_hz = 109800000000.0
lo2_hz = 4000000000.0
sb1 = "USB"
sb2 = "USB"
```

---

## 6.3 TOML 全体の構造

### トップレベル

- `schema_version`
  - 現在は `1`
- `config_name`
  - 任意の設定名
- `config_description`
  - 任意の説明文
- `[global]`
  - 全 stream 共通の設定
- `[provenance]`
  - 任意の来歴情報
- `[[spectrometers]]`
  - stream ごとの設定ブロック

---

## 6.4 `[global]` の主な項目

- `db_namespace`
  - DB テーブル名の namespace/prefix
  - 例: `necst`
- `telescope`
  - テーブル名解決に使う望遠鏡名
  - 例: `OMU1P85M`
- `output_layout`
  - 通常は `merged`
- `time_sort`
  - 通常は `true`

single-beam ではこの 2 つだけでも十分です。

```toml
[global]
db_namespace = "necst"
telescope = "OMU1P85M"
```

---

## 6.5 `[[spectrometers]]` の主な項目

### 必須に近い項目

- `name`
  - stream の論理名
- `fdnum`
  - feed 番号
- `ifnum`
  - IF 番号
- `plnum`
  - 偏波番号
- `polariza`
  - 偏波コード
- `beam_id`
  - beam の識別子
- `db_stream_name`
  - RawData 側の spectral stream 名

### あるとよい項目

- `frontend`
- `backend`
- `sampler`
- `db_table_name`
  - stream 名ではなく、DB テーブル名を直接指定したい場合
- `channel_slice`
  - converter 時点で channel を切る場合

### `polariza` の注意

`polariza` は空文字にせず、正式コードを使ってください。たとえば

- `XX`, `YY`, `XY`, `YX`
- `RR`, `LL`, `RL`, `LR`

です。

---

## 6.6 `[spectrometers.beam]` の説明

beam 幾何と EL 依存回転を表します。

- `az_offset_arcsec`
  - tangent-plane X = `dAz * cos(El)` [arcsec]
- `el_offset_arcsec`
  - tangent-plane Y = `dEl` [arcsec]
- `rotation_mode`
  - `none` または `elevation`
- `reference_angle_deg`
  - 回転の基準角
- `rotation_sign`
  - 回転の符号
- `dewar_angle_deg`
  - 固定の位相項
- `beam_model_version`
  - 任意の版情報

### single-beam の基本

single-beam では通常、beam offset は 0 でよく、回転も不要です。

```toml
[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "none"
reference_angle_deg = 0.0
rotation_sign = 1.0
dewar_angle_deg = 0.0
```

### 補足

`rotation_mode = "elevation"` を使うのは、主として multi-beam で beam 配置が EL に応じて回転する場合です。single-beam では通常不要です。

---

## 6.7 `[spectrometers.frequency_axis]` の説明

converter が spectral axis を構成するための情報です。

### 最低限必要な考え方

- `nchan`
  - チャンネル数
- `definition_mode`
  - 周波数軸の定義方法
- `restfreq_hz`
  - RESTFREQ
- `ctype1`, `cunit1`, `specsys`, `veldef`
  - 軸の種別と座標系

### よく使う定義法

`definition_mode = "first_center_and_delta"`

この場合、次を指定します。

- `first_channel_center_hz`
- `channel_spacing_hz`
- `crpix1`

### 例

```toml
[spectrometers.frequency_axis]
nchan = 32768
definition_mode = "first_center_and_delta"
first_channel_center_hz = 0.0
channel_spacing_hz = 76296.27368999299
crpix1 = 1.0
restfreq_hz = 115271201800.0
ctype1 = "FREQ"
cunit1 = "Hz"
specsys = "TOPOCENT"
veldef = "RADIO"
store_freq_column = false
```

### 注意

- `restfreq_hz` はできるだけ明示してください
- `ctype1` は現在 `FREQ` を前提とします
- `specsys` は通常 `TOPOCENT` または `LSRK`

---

## 6.8 `[spectrometers.local_oscillators]` の説明

LO チェーンと sideband を記述します。

- `lo1_hz`, `lo2_hz`, `lo3_hz`
- `sb1`, `sb2`, `sb3`
- 必要なら `obsfreq_hz`, `imagfreq_hz`

### 例

```toml
[spectrometers.local_oscillators]
lo1_hz = 109800000000.0
lo2_hz = 4000000000.0
sb1 = "USB"
sb2 = "USB"
```

### 注意

- LO/sideband 情報を入れるなら、整合した sideband 情報を入れてください
- `obsfreq_hz` を明示するか、少なくとも `restfreq_hz` を正にしてください

---

## 6.9 `[spectrometers.override]` の説明

必要に応じて row metadata を上書きするための領域です。single-beam の初期運用では通常不要です。

---

## 6.10 よくある注意点

### `beam_id`

single-beam でも、明示的に `beam_id = "B00"` と書くことを勧めます。

### `db_stream_name`

`db_stream_name` は **RawData 内の実際の spectral stream 名に一致**させてください。ここがずれると converter も extract 系も失敗します。

### `polariza`

空文字にしないでください。`XX` / `YY` などの正式コードを使ってください。

### beam 回転

single-beam では普通、beam 回転は不要です。迷ったら `rotation_mode = "none"` で始めてください。

---

## 7. single-beam の実務フロー

## 7.1 まず太陽スキャン解析

```bash
sunscan_singlebeam "$RAW" \
  --telescope OMU1P85M \
  --tel-loaddata OMU1p85m \
  --planet sun \
  --spectral-name xffts-board1 \
  --outdir out_single \
  --azel-source encoder \
  --altaz-apply none \
  --profile-xlim-deg 1.0 \
  --ripple-preset auto \
  --edge-fit-win-deg 0.15 \
  --edge-fit-threshold 0.20 \
  --hpbw-init-arcsec 324.0
```

確認するもの:

- summary の `center_az`, `center_el`
- `HPBW_az`, `HPBW_el`
- `fitAZ`, `fitEL`
- `sun_scan_summary_*.csv`
- derivative-fit PNG

## 7.2 次に SDFITS 変換

```bash
necst_v4_sdfits_converter "$RAW" \
  --spectrometer-config beams.toml \
  --out out_singlebeam.fits \
  --object OrionKL \
  --project ProjectID \
  --observer YourName
```

確認するもの:

- 出力 FITS が作成されるか
- `OBJECT`, `OBSERVER`, `PROJID` が期待通りか
- 周波数軸 (`RESTFREQ`, `CRVAL1`, `CDELT1`, `CTYPE1`, `CUNIT1`) が妥当か
- beam offset を入れていないなら、single-beam として中心で扱われているか

---

## 8. 推奨事項

- single-beam 解析の標準コマンドは `sunscan_singlebeam` に統一する
- `sun_scan_v4_v5.py` は legacy / 比較確認用として残す
- converter は、できるだけ `--spectrometer-config` を使う
- single-beam でも `beams.toml` を作っておくと、後の multi-beam 展開が容易になる

---

## 9. 次の段階

single-beam の運用が安定したら、次は multi-beam 用説明書へ進みます。multi-beam 側では、同じく `pyproject.toml` の console_scripts を前提に、次のコマンド群で説明する予定です。

- `sunscan_extract_multibeam`
- `sunscan_fit_multibeam`
- `sunscan_multibeam`
- `check_spectrometer_config`

