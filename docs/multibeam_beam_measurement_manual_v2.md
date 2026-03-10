# Multi-beam 用 beam 測定パッケージ 説明書

## 1. 本書の目的

本書は、multi-beam 受信機の太陽スキャン解析を、次のコマンド群で実行するための利用者向け説明書です。

- `check_spectrometer_config`
- `sunscan_extract_multibeam`
- `sunscan_fit_multibeam`
- `sunscan_make_pseudo_multibeam`
- `sunscan_multibeam`

本書では、これらが **`pyproject.toml` の `console_scripts` に登録されている**前提で説明します。したがって実行例はすべて

```bash
check_spectrometer_config ...
sunscan_extract_multibeam ...
sunscan_fit_multibeam ...
sunscan_make_pseudo_multibeam ...
```

の形式で記載します。

`sun_scan_v4_v5.py` は **legacy / 比較確認用** です。通常運用では

- single-beam: `sunscan_singlebeam`
- multi-beam: 本書のコマンド群

を使ってください。

---

## 2. まず何をしているのか

multi-beam 解析は、single-beam の解析を複数 stream に対して実行し、その結果から beam 配置を推定する 2 段階処理です。

### 第1段階: 各 stream の single-beam 解析

`sunscan_extract_multibeam` は、`beams.toml` に書かれた全 stream について `sunscan_singlebeam` 相当の解析を実行し、各 stream ごとに

- 太陽中心に対する offset
- HPBW
- fit の成否
- representative EL

を求めます。

### 第2段階: beam 幾何 fit

`sunscan_fit_multibeam` は、第1段階で得た結果をまとめて読み、

- beam がどのように並んでいるか
- EL に応じて offset が回転するかどうか
- 中心 beam を基準にした方が良いか、仮想中心で扱う方が良いか

を求めます。

したがって、**RawData からいきなり beam 配置が出るわけではありません。**
まず各 stream の太陽中心を測り、その後に beam 幾何を fit します。

---

## 3. 一般的な流れ

通常の流れは次のとおりです。

1. `beams.toml` を用意する
2. `check_spectrometer_config` で設定ファイルを検査する
3. `sunscan_extract_multibeam` で各 stream を解析する
4. `sunscan_fit_multibeam` で beam 幾何を求める
5. `beam_model_*.toml` を確認し、必要なら `necst_v4_sdfits_converter` に反映する

最初の実データでは、`--model both` で **center_beam** と **virtual_center** の両方を出して比較するのが安全です。

---

## 4. single-beam と multi-beam の違い

### single-beam

入力:
- RawData
- 1 本の spectrometer stream

コマンド:
- `sunscan_singlebeam`

主な出力:
- `sun_scan_summary_*.csv`
- derivative fit PNG
- 標準出力 summary

### multi-beam

入力:
- RawData
- 複数 stream の対応表 `beams.toml`

コマンド:
- `sunscan_extract_multibeam`
- `sunscan_fit_multibeam`

主な出力:
- 各 stream ごとの single-beam 相当結果
- `sunscan_multibeam_scan_summary_*.csv`
- `fit_summary.txt`
- `beam_fit_results_*.csv`
- `beam_fit_residuals_*.csv`
- `beam_model_*.toml`

つまり multi-beam は、**single-beam 解析を多数回実行して、その結果をまとめて beam 配置を求める仕組み**です。

---

## 5. 何を事前に用意すればよいか

最低限必要なのは次の 3 つです。

1. RawData ディレクトリ
2. `beams.toml`
3. single-beam で実績のある解析パラメータ

### 5.1 RawData

- NECST DB が入った観測ディレクトリ
- できれば全 beam / 全偏波が同時に記録されていること
- 最初は太陽スキャンが明瞭で、single-beam でもよく解析できるデータを使うのが安全です

### 5.2 解析パラメータ

初回の multi-beam 実行では、まず `sunscan_singlebeam` で成功している条件をそのまま使ってください。例:

```bash
--azel-source encoder \
--altaz-apply none \
--profile-xlim-deg 1.0 \
--ripple-preset auto \
--edge-fit-win-deg 0.15 \
--edge-fit-threshold 0.20 \
--hpbw-init-arcsec 324.0
```

---

## 6. `beams.toml` とは何か

### 6.1 一言でいうと

`beams.toml` は、**どの stream がどの beam か、周波数と偏波は何か、converter にどう渡すか**をまとめた設定ファイルです。

### 6.2 converter と何が同じか

`necst_v4_sdfits_converter` に渡す設定と、基本構造は同じです。共通する主な項目は次です。

- `[[spectrometers]]`
- `name`
- `fdnum`
- `ifnum`
- `plnum`
- `polariza`
- `beam_id`
- `db_stream_name`
- `[spectrometers.frequency_axis]`
  - `restfreq_hz`
- `[spectrometers.beam]`
  - `az_offset_arcsec`
  - `el_offset_arcsec`
  - `rotation_mode`
  - `reference_angle_deg`
  - `rotation_sign`
  - `dewar_angle_deg`

つまり、**beam 測定パッケージは converter 用設定をそのまま再利用できるように設計されています。**

### 6.3 converter と何が違うか

multi-beam 解析では、converter には無くても重要な情報があります。

代表例:

- `beam_fit_use = true`
  - 同じ `beam_id` に複数 stream（例: XX/YY）がある場合、どれを幾何 fit の primary stream に使うかを決めます。

また、multi-beam 解析では `[spectrometers.beam]` の値に 2 つの意味があります。

1. **観測前の名目値 / 初期値**
   - nominal 配置
   - pseudo multi-beam dry-run
2. **fit 後の最終値**
   - `beam_model_*.toml` に出力される
   - converter に反映する候補になる

つまり、最初に書く `beams.toml` は **仮設定** です。実データ解析後に、fit 結果に基づいて更新されることがあります。

---

## 7. `beams.toml` の最小例

最小限必要なのは、各 `[[spectrometers]]` について次です。

- `name`
- `fdnum`
- `ifnum`
- `plnum`
- `polariza`
- `beam_id`
- `db_stream_name`
- `restfreq_hz`
- `[spectrometers.beam]`

偏波コード `polariza` は、`XX`, `YY`, `XY`, `YX`, `RR`, `LL`, `RL`, `LR` の正式コードを使ってください。

---

## 8. 具体例: 中心 230 GHz 両偏波、外側 115 GHz 4 beam 両偏波

ここでは、実際に使いやすい例を **説明書の中に埋め込んで**示します。

想定する受信機構成:

- 物理 beam は 5 本
- 中心 beam `B00` は 230 GHz
- 外側 4 beam `B01`〜`B04` は 115 GHz
- 各 beam に XX/YY の 2 偏波がある
- 幾何 fit には各 beam の XX を primary stream として使う
- nominal offset は
  - B01: +Az 方向
  - B02: +El 方向
  - B03: -Az 方向
  - B04: -El 方向
- EL 依存回転を仮定する dry-run 用に `rotation_mode="elevation"`, `rotation_sign=1` を入れている

```toml
schema_version = 1

[[spectrometers]]
name = "b00_230_xx"
fdnum = 0
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B00"
db_stream_name = "b00_230_xx"
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 2.30538e11

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b00_230_yy"
fdnum = 1
ifnum = 0
plnum = 1
polariza = "YY"
beam_id = "B00"
db_stream_name = "b00_230_yy"

[spectrometers.frequency_axis]
restfreq_hz = 2.30538e11

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b01_115_xx"
fdnum = 2
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B01"
db_stream_name = "b01_115_xx"
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 1.152712018e11

[spectrometers.beam]
az_offset_arcsec = 600.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b01_115_yy"
fdnum = 3
ifnum = 0
plnum = 1
polariza = "YY"
beam_id = "B01"
db_stream_name = "b01_115_yy"

[spectrometers.frequency_axis]
restfreq_hz = 1.152712018e11

[spectrometers.beam]
az_offset_arcsec = 600.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b02_115_xx"
fdnum = 4
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B02"
db_stream_name = "b02_115_xx"
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 1.152712018e11

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 600.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b02_115_yy"
fdnum = 5
ifnum = 0
plnum = 1
polariza = "YY"
beam_id = "B02"
db_stream_name = "b02_115_yy"

[spectrometers.frequency_axis]
restfreq_hz = 1.152712018e11

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = 600.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b03_115_xx"
fdnum = 6
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B03"
db_stream_name = "b03_115_xx"
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 1.152712018e11

[spectrometers.beam]
az_offset_arcsec = -600.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b03_115_yy"
fdnum = 7
ifnum = 0
plnum = 1
polariza = "YY"
beam_id = "B03"
db_stream_name = "b03_115_yy"

[spectrometers.frequency_axis]
restfreq_hz = 1.152712018e11

[spectrometers.beam]
az_offset_arcsec = -600.0
el_offset_arcsec = 0.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b04_115_xx"
fdnum = 8
ifnum = 0
plnum = 0
polariza = "XX"
beam_id = "B04"
db_stream_name = "b04_115_xx"
beam_fit_use = true

[spectrometers.frequency_axis]
restfreq_hz = 1.152712018e11

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = -600.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0

[[spectrometers]]
name = "b04_115_yy"
fdnum = 9
ifnum = 0
plnum = 1
polariza = "YY"
beam_id = "B04"
db_stream_name = "b04_115_yy"

[spectrometers.frequency_axis]
restfreq_hz = 1.152712018e11

[spectrometers.beam]
az_offset_arcsec = 0.0
el_offset_arcsec = -600.0
rotation_mode = "elevation"
reference_angle_deg = 35.0
rotation_sign = 1.0
dewar_angle_deg = 0.0
```

### この例のポイント

- `beam_id` は物理 beam を表します
- `name` / `db_stream_name` は stream を表します
- `beam_fit_use = true` は、各 beam の XX を primary stream にしているという意味です
- `restfreq_hz` は 230 GHz / 115 GHz の違いを表します
- nominal offset は dry-run や初期値として使えます

---

## 9. どのモデルがあるのか

`sunscan_fit_multibeam` には、beam 幾何の扱いとして 2 つのモデルがあります。

### 9.1 `center_beam` モデル

**本当に中心 beam がある**ときのモデルです。

考え方:
- 中心 beam（例: `B00`）の位置を基準にする
- 各 outer beam の相対 offset を fit する
- 中心 beam が信頼できるなら、最も素直で解釈しやすい

向いているケース:
- 中心 beam が実在する
- 中心 beam の fit が安定している
- ポインティングの中心をその beam に取りたい

### 9.2 `virtual_center` モデル

**中心 beam が無い、または中心 beam を基準にしにくい**ときのモデルです。

考え方:
- 全 beam の配置から仮想的な回転中心を求める
- その仮想中心に対する各 beam の配置を fit する

向いているケース:
- 中心 beam が存在しない
- 中心 beam があっても信頼性が低い
- 配列全体の幾何をまず安定に決めたい

### 9.3 `--model both` とは何か

`--model both` は、**両方のモデルを走らせて比較する**という意味です。

出力:
- `beam_fit_results_center_beam.csv`
- `beam_fit_results_virtual_center.csv`
- `beam_model_center_beam.toml`
- `beam_model_virtual_center.toml`
- `fit_summary.txt`

`fit_summary.txt` では、両モデルの residual を比較して `recommended_model` が示されます。

最初の実データでは、**`--model both` を推奨**します。

---

## 10. 中心 beam がある場合と無い場合、結果はどう違うか

### 中心 beam がある場合

`--center-beam-id B00` を指定し、`--model both` で走らせると、通常は

- `center_beam` 結果
- `virtual_center` 結果

の両方が出ます。

このとき、**ポインティング補正や converter 反映に使う第一候補**は通常 `center_beam` 側です。
理由は、実在する中心 beam に対する相対配置として解釈しやすいからです。

ただし、
- 中心 beam の fit が悪い
- `center_beam` の residual が大きい
- `virtual_center` の方が明らかに安定

なら、`virtual_center` を採用した方がよい場合があります。

### 中心 beam がない場合

`center_beam` モデルは物理的な意味を持ちません。したがって、基本は

- `--model virtual_center`

または

- `--model both` を走らせても、最終採用は `virtual_center`

になります。

この場合、**配列全体の相対配置は求まりますが、絶対的なポインティング中心の扱いは別途考える必要があります。**

---

## 11. どの出力を何に使うべきか

### 11.1 まず見るべきもの

- `fit_summary.txt`
  - 全体の結果要約
  - `recommended_model`
  - residual の大きさ
  - `rotation_sign`
- `beam_fit_residuals_*.csv`
  - どの beam / どの scan が悪いか

### 11.2 ポインティングパラメータ解析に使うもの

通常、最終的に使う候補は次です。

- 中心 beam がある: `beam_model_center_beam.toml`
- 中心 beam がない、または不安定: `beam_model_virtual_center.toml`

ただし、これは **必ず `fit_summary.txt` と residual を確認したうえで** 決めてください。

### 11.3 converter に反映するもの

converter に渡す候補は、最終的に採用したモデルの

- `beam_model_center_beam.toml`
- または `beam_model_virtual_center.toml`

です。

この TOML は converter 互換形式です。必要なら既存の `beams.toml` の `[spectrometers.beam]` を、この結果で更新して使います。

---

## 12. オフセットが回転する場合と回転しない場合

ここは非常に重要です。

### 12.1 回転しない場合

beam offset が EL に依存せず固定なら、

```toml
rotation_mode = "none"
```

です。

このとき `az_offset_arcsec`, `el_offset_arcsec` は固定オフセットとして使われます。

### 12.2 回転する場合

beam offset が EL に応じて回転するなら、

```toml
rotation_mode = "elevation"
```

です。

converter 互換の回転式は、概念的には

\[
\theta(EL) = rotation\_sign \times (EL - reference\_angle\_deg) + dewar\_angle\_deg
\]

です。

ここで
- `reference_angle_deg`: 基準 EL
- `rotation_sign`: 回転方向（通常は -1, 0, +1）
- `dewar_angle_deg`: 基準姿勢での向き

です。

fit 結果 CSV には、診断用に

- `rotation_slope_deg_per_deg`

も出ます。これは **EL 1 度に対して回転角が何度変わるか** を見るための量です。converter の直接入力ではなく、fit の健全性確認用と考えてください。

### 12.3 実務上どう判断するか

- 実データで `rotation_sign = 0` に近い
- `rotation_slope_deg_per_deg` もほぼ 0

なら、回転しないモデルでよい可能性があります。

逆に
- 複数 EL を入れると outer beam の relative position が系統的に回る
- `rotation_sign` が安定して ±1 に出る

なら、`rotation_mode="elevation"` を採用すべきです。

---

## 13. 実行例: 一般的な流れ

### 13.1 設定ファイル検査

```bash
check_spectrometer_config beams.toml --out-csv config_check.csv
```

### 13.2 各 stream の抽出解析

```bash
sunscan_extract_multibeam RAWDATA \
  --spectrometer-config beams.toml \
  --outdir out_extract \
  --continue-on-error \
  --azel-source encoder \
  --altaz-apply none \
  --profile-xlim-deg 1.0 \
  --ripple-preset auto \
  --edge-fit-win-deg 0.15 \
  --edge-fit-threshold 0.20 \
  --hpbw-init-arcsec 324.0
```

### 13.3 beam 幾何 fit

```bash
sunscan_fit_multibeam \
  out_extract/sunscan_multibeam_scan_summary_<TAG>.csv \
  --spectrometer-config beams.toml \
  --outdir out_fit \
  --center-beam-id B00 \
  --model both
```

---

## 14. 実行後に最低限確認すること

### 抽出段階

- `sunscan_multibeam_manifest_*.csv`
- `spectrometer_stream_table_*.csv`
- `sunscan_multibeam_scan_summary_*.csv`

確認点:
- 失敗した stream がないか
- `beam_id`, `stream_name`, `restfreq_hz`, `polariza` が正しいか
- `center_az_deg`, `center_el_deg`, `hpbw_*`, `fit_ok_*` が妥当か

### fit 段階

- `fit_summary.txt`
- `beam_fit_results_*.csv`
- `beam_fit_residuals_*.csv`
- `beam_model_*.toml`

確認点:
- `selected_rows`, `rejected_rows` が妥当か
- `rotation_sign` が物理的におかしくないか
- `recommended_model` はどちらか
- residual が特定 beam / 特定 EL だけで大きくないか

---

## 15. 初心者向けの実務上の助言

1. 最初から高度な設定をしないでください。まずは single-beam で成功している条件をそのまま使います。
2. 最初の実データでは、必ず `--model both` を使ってください。
3. `beam_model_*.toml` をすぐ converter に入れる前に、`fit_summary.txt` と residual を確認してください。
4. 中心 beam があるなら、まず `center_beam` を第一候補にします。
5. 中心 beam が無い、または不安定なら、`virtual_center` を採用します。
6. 実データ初回では、XX/YY を無理に同時統合しないで、まず primary stream だけで幾何を確かめるのが安全です。

---

## 16. まとめ

multi-beam 解析で一番大事なのは、

- まず各 stream の太陽中心を正しく測ること
- その後に beam 幾何を fit すること
- 中心 beam がある場合と無い場合を区別すること
- 回転するモデルか、回転しないモデルかを意識すること
- 最終的に使う TOML を、`fit_summary.txt` と residual を見て選ぶこと

です。

最初は複雑に見えますが、流れとしては

1. `beams.toml` を用意する
2. `check_spectrometer_config`
3. `sunscan_extract_multibeam`
4. `sunscan_fit_multibeam --model both`
5. `beam_model_*.toml` を確認する

の 5 段階です。

この順で進めれば、観測経験がまだ多くなくても、どこで何を見ればよいかが分かるはずです。
