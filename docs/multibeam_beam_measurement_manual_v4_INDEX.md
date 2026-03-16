# Multi-beam 用 beam 測定パッケージ 説明書 v4  
## 分冊版インデックス

この版では、旧 `multibeam_beam_measurement_manual_v3.md` を、現行実装に合わせて**役割分離を明示する構成**へ組み替えています。

## 構成

### 1. `multibeam_beam_measurement_manual_v4_01_architecture_and_shared_config.md`
目的:
- converter と sunscan の役割分担を明示する
- spectrometer config TOML のうち、何が共有情報で、何が sunscan 固有で、何が converter 側の情報かを整理する
- multi-beam package の全体像と出力物を説明する

読むタイミング:
- 最初に読む

### 2. `multibeam_beam_measurement_manual_v4_02_sunscan_command_reference.md`
目的:
- multi-beam package 側の CLI 引数と各パラメータを、現行実装に基づいて詳細に説明する
- `sunscan_singlebeam`, `sunscan_extract_multibeam`, `sunscan_fit_multibeam`, `sunscan_make_pseudo_multibeam`, `sunscan_multibeam`, `check_spectrometer_config` を扱う

読むタイミング:
- 実行時に参照する

### 3. `multibeam_beam_measurement_manual_v4_03_cookbook.md`
目的:
- single / multi / pseudo / fit / config 検証 / stream の一時無効化など、具体的な運用例をまとめる
- converter-compatible TOML の具体例もここに載せる

読むタイミング:
- 実際に動かすとき

## この構成にした理由

旧版では「multi-beam 解析の流れ」の説明は非常に分かりやすかった一方で、現行コードでは次の整理が重要になっています。

- converter と sunscan は**別機能**
- ただし spectrometer config TOML の一部は**共有情報**
- `enabled / use_for_convert / use_for_sunscan / use_for_fit / beam_fit_use` が入ったため、**config semantics を先に説明した方が安全**
- `sunscan_singlebeam` が `--spectrometer-config` を直接扱えるようになったので、single と multi の関係を再整理した方がよい
- pseudo / wrapper / check-config まで含めると、1 本の説明書に全部を平坦に並べるより、分冊の方が参照しやすい

## 推奨読書順

1. 01 architecture / shared config
2. 02 command reference
3. 03 cookbook

## 使い方の最短ルート

- はじめて: 01 → 03
- 実行中の確認: 02
- 設定ファイルを作る: 01 → 03
