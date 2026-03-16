# NECST converter / sunscan ドキュメントセット概要

## 1. このドキュメントセットの目的

このドキュメントセットは、現行の以下のスクリプト群を対象にしています。

- `necst_v4_sdfits_converter.py`
- `multibeam_beam_measurement/sunscan_singlebeam.py`
- `multibeam_beam_measurement/sunscan_extract_multibeam.py`
- `multibeam_beam_measurement/sunscan_fit_multibeam.py`
- `multibeam_beam_measurement/sunscan_multibeam.py`
- `multibeam_beam_measurement/synthetic_multibeam.py`
- `multibeam_beam_measurement/check_spectrometer_config.py`

ここで重要なのは、**converter** と **sunscan** は役割が異なる、という点です。

### converter の役割

converter は、NECST RawData / necstdb から **SDFITS を作る**ためのツールです。

主な責務は次です。

- spectrometer stream の読み出し
- encoder / altaz / weather / WCS の統合
- beam model を適用した Az/El の算出
- 必要に応じた Az/El 由来の RA/DEC の算出
- SDFITS 行の生成
- multi-stream 構成の1本の FITS 出力

converter は、**観測データを FITS として整形・保存するツール**です。

### sunscan の役割

sunscan は、主に太陽スキャン観測から

- スキャン中心
- HPBW
- スキャン速度
- limb fit の品質
- multi-beam の相対配置・回転

を求めるための**解析ツール**です。

sunscan はさらに 2 層に分かれます。

- **singlebeam**: 1 本の stream を解析する
- **multibeam**: 複数 stream をまとめて扱い、extract / fit / pseudo / config-check を行う

### shared config の役割

converter と sunscan は、同じ TOML の `[[spectrometers]]` 定義を共有できます。

ただし、重要な注意があります。

- 同じキーでも、**converter と sunscan で意味が完全に同じとは限りません**
- あるキーは converter では有効、sunscan では未使用、あるいはその逆、ということがあります
- したがって、設定ファイルを 1 本で共用するときは、**shared config reference** を必ず参照してください

この点を明確にするため、本ドキュメントは以下のように分冊しています。

---

## 2. 推奨される読み方

### A. converter を使いたいとき

まず以下を読んでください。

1. `10_converter_manual.md`
2. 必要なら `40_shared_config_reference.md`
3. 実行例は `50_cookbook.md`

### B. single-beam の太陽スキャン解析をしたいとき

1. `20_sunscan_singlebeam_manual.md`
2. 必要なら `40_shared_config_reference.md`
3. 実行例は `50_cookbook.md`

### C. multi-beam の抽出・beam fitting をしたいとき

1. `30_sunscan_multibeam_manual.md`
2. `40_shared_config_reference.md`
3. 実行例は `50_cookbook.md`

---

## 3. 各ファイルの内容

### `10_converter_manual.md`

converter 専用の詳細マニュアルです。

- 役割
- RawData から SDFITS までの流れ
- converter の全 CLI パラメータ
- RA/DEC の決め方
- `--radec-azel-source`
- spectrometer config のうち converter が使う項目
- converter 特有の注意点
- troubleshooting

### `20_sunscan_singlebeam_manual.md`

single-beam 解析専用の詳細マニュアルです。

- 1 stream 解析の目的
- `SunScanAnalysisConfig` の構造
- singlebeam の全 CLI パラメータ
- scan trim / ripple / edge fit の意味
- 出力 CSV / PNG
- singlebeam 特有の注意点

### `30_sunscan_multibeam_manual.md`

multibeam 解析専用の詳細マニュアルです。

- extract / fit / pseudo / check-config / wrapper の役割
- stream usage flags
- primary stream 解決
- beam fitting の数学的意味
- `center_beam` と `virtual_center` の違い
- sigma clipping と minimum-count 条件
- multibeam の全 CLI パラメータ
- 出力 CSV / TOML / manifest / snapshot

### `40_shared_config_reference.md`

shared spectrometer config TOML の完全リファレンスです。

- `[global]`
- `[[spectrometers]]`
- `[spectrometers.frequency_axis]`
- `[spectrometers.local_oscillators]`
- `[spectrometers.override]`
- `[spectrometers.beam]`
- `enabled / use_for_convert / use_for_sunscan / use_for_fit / beam_fit_use`
- converter / singlebeam / multibeam のサポート状況

### `50_cookbook.md`

実際の運用例をまとめた cookbook です。

- converter 単独
- singlebeam
- multibeam extract / fit
- pseudo multibeam
- 問題のある stream の一時無効化
- `db_table_name` による table 名の固定
- `--radec-azel-source cmd` を使った position-switch 的な運用例

---

## 4. converter と sunscan の境界を一言で言うと

- converter は **保存系**
- sunscan は **解析系**

です。

### converter でやること

- spectrum を読む
- pointing / WCS / weather を合わせる
- FITS を作る

### sunscan でやること

- 太陽 limb fit をする
- scan 中心や HPBW を出す
- 複数ビームの相対配置を復元する
- 将来使う beam model TOML を更新する

---

## 5. 用語整理

### stream

ここでの stream は、1 本の spectrometer 入力系列です。典型的には

- board ごとの stream
- 偏波ごとの stream
- IF / FDNUM / PLNUM で区別される stream

を意味します。

### beam

beam は、物理的あるいは論理的な受信ビームです。

- 1 つの beam に XX / YY の 2 stream が対応することがある
- multi-beam fit では、同一 beam の複数 stream のうち、通常は 1 本だけを primary stream として使う

### primary stream

同一 `beam_id` に複数 stream が存在する場合、multi-beam fit に使う代表 stream です。

通常は

- `use_for_fit = true`
- かつ `beam_fit_use = true`

の stream が primary stream になります。

### boresight

converter における boresight は、現在の実装では **encoder から補正量 `dlon / dlat` を引いた中心軸**です。

### beam-center

boresight に stream 固有の beam offset / rotation を適用したものです。

### `cmd` Az/El source

converter における `--radec-azel-source cmd` は、**altaz/cmd 系列から `dlon / dlat` を引いた boresight を作り、その上に beam offset を適用した Az/El** を Az/El 由来の RA/DEC に使うモードです。

---

## 6. 現在の推奨運用

### converter

- spectrometer config TOML を使う
- stream 名と table 名が一致しない場合は `db_table_name` を使って固定する
- `--radec-azel-source` を使うときは、何を RA/DEC 化しているかを意識する

### singlebeam

- まず 1 stream を単独で十分に確認する
- ripple / trim / edge fit を無理に強くしすぎない
- summary CSV を保存し、後段で pseudo / fit に再利用する

### multibeam

- まず `check-config` で TOML を検証する
- extract の manifest / stream_table / config_snapshot を必ず確認する
- fit は primary stream の選び方を明示する
- 実データが少ないときは pseudo multibeam で dry-run を行う

---

## 7. 注意

このドキュメントセットは、**現行コード**を基準にしています。

- 過去に提案されたが未実装の案
- 将来実装したい改善案
- 現在のコードに存在しない CLI / TOML semantics

は、原則として書いていません。

今のスクリプトが実際にどう動くか、を最優先に記述しています。
