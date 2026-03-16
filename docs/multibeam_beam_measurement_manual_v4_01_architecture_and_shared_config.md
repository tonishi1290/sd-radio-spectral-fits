# Multi-beam 用 beam 測定パッケージ 説明書 v4-01  
## architecture・converter との境界・shared spectrometer config・出力仕様

## 1. 本書の位置づけ

本書は、multi-beam 用 beam 測定パッケージの**全体設計**と、converter-compatible spectrometer config TOML の**共有部分 / 非共有部分**を整理するための文書です。

この分冊で最も重要なのは、次の整理です。

- **converter** は converter の役割を持つ
- **sunscan package** は sunscan package の役割を持つ
- ただし **spectrometer config TOML の一部は共有情報**である
- したがって、説明書も
  - 何が converter の責務か
  - 何が sunscan の責務か
  - 何が shared information か
  を先に分けて書く方が分かりやすい

---

## 2. 役割分担

## 2.1 converter の役割

converter の主目的は、**RawData から SDFITS を作ること**です。  
主に扱うのは

- spectrometer stream の読込
- pointing / WCS / RA/DEC / meteorology の解決
- FITS への書き出し
- beam モデルの適用
- converter 独自の出力形式制御

です。

したがって converter では、

- DB table 名
- WCS
- LO / IF
- output layout
- weather / refraction
- RA/DEC 変換方法

などが重要です。

## 2.2 sunscan package の役割

sunscan package の主目的は、**太陽スキャンの single-beam 解析と multi-beam 幾何 fit**です。  
主に扱うのは

- Ta* または power プロファイルからの limb fit
- scan trimming
- ripple 除去
- representative Az / El
- beam center / HPBW の推定
- 複数 beam の相対幾何 fit

です。

したがって sunscan package では、

- scan の切り出し
- trim 条件
- ripple 条件
- edge fit 条件
- summary CSV / manifest / fit summary

などが重要です。

## 2.3 shared information の役割

converter と sunscan が共通で持つと便利なのは、**各 stream / beam の定義情報**です。

代表例:

- `name`
- `beam_id`
- `db_stream_name`
- `db_table_name`
- `fdnum`
- `ifnum`
- `plnum`
- `polariza`
- `restfreq_hz`
- beam offset / rotation model
- stream usage flags

これらは spectrometer config TOML に置くのが自然です。

---

## 3. 何が shared で、何が shared でないか

## 3.1 基本方針

spectrometer config TOML には、次の 3 種類の情報が混ざります。

### A. shared information
converter と sunscan の両方に意味がある情報。

### B. sunscan-centric information
sunscan package 側で主に意味を持つ情報。

### C. converter-centric information
converter 側で主に意味を持つ情報。

説明書では、この 3 つを混ぜない方が理解しやすいです。

---

## 4. shared information 一覧

## 4.1 `[global]` の shared information

### `db_namespace`
- DB table 名の接頭辞
- converter / singlebeam / extract で意味がある

### `telescope`
- telescope 名
- converter / singlebeam / extract で意味がある

### `tel_loaddata`
- observatory 名や legacy 名と整合を取るための telescope 名
- singlebeam / extract で記録・利用される
- converter-compatible config を human-readable に保つ意味でも有用

### `planet`
- sunscan package 側で使う
- converter では通常 object 名や source 情報と別概念
- ただし shared config に置くと single / multi で統一しやすい

### `spectrometer_time_offset_sec`
- spectrometer 側の時刻補正
- sunscan で解析に直接使う
- converter-compatible config に置いておくと意味が明確

### `encoder_shift_sec`
- Az/El 共通の encoder 時刻補正
- 互換用に残しておくと便利

### `encoder_az_time_offset_sec`
### `encoder_el_time_offset_sec`
- Az / El 個別時刻補正
- 今回の改修で重要になった shared information

## 4.2 `[[spectrometers]]` の shared information

### stream identity
- `name`
- `db_stream_name`
- `db_table_name`

### hardware / index identity
- `fdnum`
- `ifnum`
- `plnum`
- `sampler`
- `frontend`
- `backend`

### polarization / beam grouping
- `polariza`
- `beam_id`

### frequency identity
- `[spectrometers.frequency_axis].restfreq_hz`

### beam model
- `[spectrometers.beam]`
  - `az_offset_arcsec`
  - `el_offset_arcsec`
  - `rotation_mode`
  - `reference_angle_deg`
  - `rotation_sign`
  - `dewar_angle_deg`

### stream usage policy
- `enabled`
- `use_for_convert`
- `use_for_sunscan`
- `use_for_fit`
- `beam_fit_use`

---

## 5. shared でないが、TOML に同居していてよい情報

## 5.1 converter-centric

以下は主として converter 側で意味が大きい項目です。

- `channel_slice`
- `output_layout`
- `time_sort`
- weather table の制御
- WCS 詳細
- LO / IF 詳細
- converter 側の output control

これらは spectrometer config TOML に入っていてよいですが、sunscan package が本質的に使うわけではありません。

## 5.2 sunscan-centric

以下は主として sunscan CLI 側で意味が大きい項目です。

- ripple 除去条件
- trim 条件
- edge fit 条件
- profile 表示範囲
- debug plot
- `continue_on_error`

これらは現在、主に CLI / `SunScanAnalysisConfig` 側の責務であり、stream 定義そのものとは分けて考える方が分かりやすいです。

---

## 6. current config semantics

## 6.1 singlebeam / extract が `[global]` から読むもの

現行実装で singlebeam / extract が `[global]` から積極的に読むのは次です。

- `db_namespace`
- `telescope`
- `tel_loaddata`
- `planet`
- `spectrometer_time_offset_sec`
- `encoder_shift_sec`
- `encoder_az_time_offset_sec`
- `encoder_el_time_offset_sec`

CLI で同名 option を明示した場合は CLI が優先です。

## 6.2 usage flags の current semantics

### `enabled`
- 総合スイッチ
- `false` なら他の `use_for_*` が `true` でも最終的に無効とみなされる

### `use_for_convert`
- converter 側での使用可否を表す
- multi-beam package 側では、主に stream table / TOML 生成 / metadata で保持される

### `use_for_sunscan`
- singlebeam / extract での解析対象可否
- `sunscan_extract_multibeam` の通常実行で効く

### `use_for_fit`
- fit 候補集合に含めるか
- `beam_fit_use` とは別概念

### `beam_fit_use`
- 同じ `beam_id` の中で代表 stream を選ぶための印

重要:
- `use_for_fit=false` は「fit 候補から除外」
- `beam_fit_use=true` は「その beam の代表にする」

## 6.3 explicit stream selection の優先

現行実装では、`--stream-name` のような explicit selection は強い意味を持ちます。

- `sunscan_singlebeam`  
  明示 `--stream-name` があれば disabled stream でも one-off override として選択可能
- `sunscan_extract_multibeam`  
  明示 `--stream-name` があれば usage flags による自動除外をバイパスして選択可能
- `sunscan_fit_multibeam`  
  明示 `--fit-stream-name` があれば primary stream 解決を明示指定で上書き可能

---

## 7. time offset semantics

## 7.1 基本式

現行 sunscan package では概念的に

- `t_spec_use = t_spec_raw + spectrometer_time_offset_sec`
- `t_az_use = t_enc_raw + encoder_shift_sec + encoder_az_time_offset_sec`
- `t_el_use = t_enc_raw + encoder_shift_sec + encoder_el_time_offset_sec`

です。

定義はすべて

> 補正後時刻 = 記録時刻 + offset

です。

## 7.2 推奨初期値

最初の運用では、共通 offset を使いすぎないために

```toml
spectrometer_time_offset_sec = -0.07
encoder_shift_sec = 0.0
encoder_az_time_offset_sec = 0.00
encoder_el_time_offset_sec = +0.04
```

のようにする方が解釈しやすいです。

## 7.3 summary / snapshot に残るもの

現行実装では、少なくとも以下が summary / snapshot / manifest に残ります。

- `spec_time_basis`
- `spec_time_suffix`
- `spec_time_fallback_field`
- `spec_time_example`
- `spectrometer_time_offset_sec`
- `encoder_shift_sec`
- `encoder_az_time_offset_sec`
- `encoder_el_time_offset_sec`
- `encoder_vavg_sec`

これはあとから「何を適用していたか」を追跡するために重要です。

---

## 8. `beam_id` と `beam_fit_use` の考え方

## 8.1 典型形

たとえば 1 物理 beam あたり XX / YY の 2 stream がある場合:

- `center230_xx` → `beam_id="B00"`
- `center230_yy` → `beam_id="B00"`

のようになります。

このとき fit では代表 1 本が必要です。  
通常は XX を代表にして

- `center230_xx`: `beam_fit_use = true`
- `center230_yy`: `beam_fit_use = false` または未指定

とします。

## 8.2 `use_for_fit` との違い

たとえば YY は fit に使わず、convert / extract では残したいなら

```toml
use_for_convert = true
use_for_sunscan = true
use_for_fit = false
beam_fit_use = false
```

とします。

逆に「XX/YY とも fit 候補には残すが、代表 stream だけ XX にしたい」なら

- XX: `use_for_fit = true`, `beam_fit_use = true`
- YY: `use_for_fit = true`, `beam_fit_use = false`

です。

---

## 9. multi-beam package の出力物

## 9.1 `sunscan_singlebeam` の出力

通常は次が中心です。

- `sun_scan_summary_<tag>.csv`
- derivative fit PNG 群
- debug plot PNG 群（`--debug-plot` 時）
- 標準出力 summary

## 9.2 `sunscan_extract_multibeam` の出力

extract は次の 4 つが特に重要です。

### `sunscan_multibeam_scan_summary_<tag>.csv`
全 stream の single-beam 結果を縦にまとめた summary。

### `sunscan_multibeam_manifest_<tag>.csv`
各 stream の実行状態一覧。

主な列:
- `stream_name`
- `beam_id`
- `spectral_name`
- `summary_csv`
- `spec_time_basis`
- `spec_time_suffix`
- `spec_time_fallback_field`
- `status`
- `skip_reason`
- `enabled`
- `use_for_convert`
- `use_for_sunscan`
- `use_for_fit`
- `beam_fit_use`
- `error`

### `spectrometer_stream_table_<tag>.csv`
config から解釈した stream table。

### `analysis_config_snapshot_<tag>.json`
解析時に使った base_config と validation 情報の snapshot。

## 9.3 `sunscan_fit_multibeam` の出力

主に次です。

- `selected_rows_for_fit.csv`
- `beam_fit_results_<model>.csv`
- `run_shift_results_<model>.csv`
- `beam_fit_residuals_<model>.csv`
- `beam_fit_rejected_<model>.csv`
- `beam_model_<model>.toml`
- `fit_summary.txt`

## 9.4 `sunscan_make_pseudo_multibeam` の出力

- `pseudo_multibeam_summary_<tag>.csv`
- `pseudo_multibeam_manifest_<tag>.csv`
- `pseudo_multibeam_stream_table_<tag>.csv`

---

## 10. `beam_model_*.toml` の意味

fit 後に出る `beam_model_*.toml` は、元の spectrometer config TOML を土台にして、  
各 beam の `[spectrometers.beam]` を fit 結果で更新したものです。

この TOML には

- `enabled`
- `use_for_convert`
- `use_for_sunscan`
- `use_for_fit`
- `beam_fit_use`

も保持されます。

したがって、

- もとの stream 定義
- nominal / final beam offset
- converter で使うための beam model

を 1 ファイル系統で管理しやすくなります。

---

## 11. `check_spectrometer_config` が見るもの

現行実装の validator は主に次を見ます。

- stream table が空でないか
- `polariza` が許可値か
- `restfreq_hz` が有限か
- `rotation_sign` が数値か
- `rotation_sign` が converter-compatible の `-1, 0, +1` から外れていないか
- primary stream が解決できるか
- primary stream の nominal offset が全て 0 ではないか
- primary stream の `rotation_mode` が全て `none` ではないか

最後の 2 つは「完全にエラー」ではなく warning です。  
特に pseudo dry-run では nominal beam offset が全部 0 だと何も検証できないため重要です。

---

## 12. `center_beam` と `virtual_center`

## 12.1 `center_beam`
本当に中心 beam がある場合のモデル。

向いているケース:
- 中心 beam が物理的に存在
- 中心 beam の fit が安定
- 配列全体をその beam 基準で解釈したい

## 12.2 `virtual_center`
仮想的な回転中心を置くモデル。

向いているケース:
- 実中心 beam が無い
- 中心 beam が不安定
- 配列全体の幾何だけをまず安定に求めたい

## 12.3 `both`
最初の実データでは `--model both` が推奨です。  
両方走らせて、`fit_summary.txt` の residual を見て採用モデルを決めるのが安全です。

---

## 13. `sunscan_multibeam` ラッパの位置づけ

`sunscan_multibeam` は便利ですが、**常に単機能コマンドと完全に対称とは限りません。**

現行実装上の実務的な理解としては、

- 通常運用では wrapper を使ってもよい
- 細かな挙動確認や特殊オプションが必要なときは単機能コマンドの方が分かりやすい

という位置づけが安全です。

特に pseudo では、**スタンドアロンの `sunscan_make_pseudo_multibeam` を直接使う方が意図が明確**です。

---

## 14. 現場運用での基本戦略

1. まず `check_spectrometer_config` で config を検査する  
2. `sunscan_singlebeam` で 1 本だけ成功させる  
3. その条件で `sunscan_extract_multibeam` を回す  
4. `sunscan_fit_multibeam --model both` で比較する  
5. `beam_model_*.toml` を見て、converter に反映するか決める

この順番が最も安全です。

---

## 15. この分冊のまとめ

この package を理解するときに最も大事なのは、

- converter と sunscan は別機能である
- ただし stream / beam / restfreq / beam model / usage flags は shared information と考えるとよい
- `enabled / use_for_* / beam_fit_use` を理解すると config 管理が非常に楽になる
- time offset は 3 分解 + 共通 offset として整理するとよい
- single / extract / fit / pseudo は責務が違う
- まず single で確かめてから multi に進むのが安全

という点です。

次の分冊 `v4-02` では、各コマンドの CLI 引数と全パラメータを詳細に説明します。
