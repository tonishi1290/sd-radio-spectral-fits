# sun_scan機能温存を最優先にした共通化仕様（再確認版・確定）

## 1. 目的
本仕様は、`sun_scan_v4_v5.py` 相当の single-beam 太陽スキャン解析機能を、将来の multi-beam 解析ソフトでも再利用できるように共通化するための設計方針を定める。

最優先要件は **既存 `sun_scan_v4_v5.py` の機能・CLI・出力・既定挙動を損ねないこと** である。

本仕様は、将来の改修を **1か所の共通コアへ反映すれば single-beam / multi-beam の両方に自然に波及する** 状態を目標とする。

---

## 2. 互換性保持対象（v1で変更禁止）
以下は v1 リファクタリングで変更してはならない。

1. CLI 引数名
2. CLI の既定値
3. 既存出力ファイル名と命名規則
4. summary CSV の列名と列順
5. summary text の基本内容と並び
6. derivative fit PNG の基本命名規則
7. `daz = wrap(Az_true - Az_sun) * cos(El_sun)`, `d_el = El_true - El_sun` という座標定義
8. `--debug-plot` の on/off 条件
9. `continue_on_error` の現在の有効範囲

### 2.1 `continue_on_error` の現在の有効範囲
実装上、`continue_on_error` が直接効いているのは以下に限る。

- AZ/EL の debug plot 生成ループ
- `build_dataframe()` 以前ではなく、その後の debug plot 生成の継続可否

一方、`run_edge_fits_and_summarize()` の **scan単位例外**は内部で `az_error` / `el_error` に吸収されるため、ここは `continue_on_error` とは独立である。

したがって、v1 では **`continue_on_error` の意味を拡張してはいけない**。既存の実装範囲だけを保つ。

---

## 3. 現行 `sun_scan_v4_v5.py` の事実
- デフォルト値はファイル冒頭の `DEFAULT_*` 定数群に集約されており、事実上の single source of truth である。
- 解析本体は limb fit だけではなく、Az/El source 選択、Ta* 化、ripple 除去、trim、steady-scan 選別、導関数生成、2ピーク Gaussian fit、summary/CSV/PNG 出力までを含む。
- `main()` が argparse と解析実行の両方を担っている。
- `trim_params` と `ripple_policy` は `main()` で辞書化されて後段に渡されているため、ここが将来の `Config` 化の自然な境界である。
- `write_scan_summary_csv()` の `preferred` 列順が、現在の canonical な CSV 順序である。

---

## 4. 共通化の基本原則

### 4.1 守るべき原則
- `sun_scan_v4_v5.py` は将来的に **薄い CLI ラッパー**にする。
- 数値処理とレポート生成は別モジュールへ移す。
- multi-beam 側は `sun_scan` を再実装せず、同じ共通コアを beam ごとに繰り返し呼ぶ。
- single-beam と multi-beam の差は **入力の束ね方**と**追加出力**だけに留める。

### 4.2 やってはいけないこと
- multi-beam 側で trim / ripple / derivative / edge fit を別実装しない。
- 既存 `sun_scan` の CLI を変更しない。
- summary CSV の列順を変えない。
- 先に multi-beam 向け都合で single-beam の既定挙動を変えない。
- core モジュールに argparse 依存を残したまま固定しない。

---

## 5. モジュール責務分界（確定）

### 5.1 `sunscan_core.py`
**責務**
- RawData → DataFrame / scan 辞書への変換
- Az/El source 切替
- chopper-wheel Ta* 化
- ripple 除去
- trim / steady-scan 選別
- derivative 計算
- limb edge fit
- per-scan の結果辞書生成

**ここへ置く関数群**
- `build_dataframe`
- `trim_scan_segment`
- `finite_difference_skip_duplicates`
- `fit_two_sides_edge`
- `resolve_ripple_cfg_for_scan`
- `run_edge_fits_and_summarize`
- `_pick_tracking_point`
- `estimate_scan_speed_deg_s`
- ripple / smoothing / tracking point 選択に関わる補助関数

**備考**
- `make_ripple_policy_from_args()` は end-state では core に置かない。
- 代わりに core には argparse 非依存の純関数 `make_ripple_policy(...)` または `make_ripple_policy_from_config(...)` を置く。
- `make_ripple_policy_from_args()` は Phase 1 では互換性確保のため一時的に `sun_scan_v4_v5.py` 側に残してよい。

**置かないもの**
- argparse
- 既存 CLI 固有の標準出力文言
- 既存ファイル命名の最終決定

### 5.2 `sunscan_report.py`
**責務**
- summary text 生成
- summary text PNG 保存
- summary CSV 出力
- derivative fit summary PNG 出力
- single-beam / multi-beam 共通のレポート出力

**ここへ置く関数群**
- `write_scan_summary_csv`
- `plot_derivative_fits_paginated`
- `save_text_summary`
- `_plot_profile_and_derivative_panel`
- `_pretty_ylabel`

**重要条件**
- 既存 `sun_scan` と同じ列名・列順・ファイル名規則を維持する。
- single-beam 互換出力は multi-beam でも再利用可能にする。

### 5.3 `sun_scan_v4_v5.py`
**責務**
- 従来どおりの CLI 入口
- argparse での引数受理
- CLI 引数 → `SunScanAnalysisConfig` への変換
- 既存どおりの出力フロー制御
- 互換性維持のための最終ラッパー

**禁止事項**
- コア数値処理をここへ残し続けない。
- 将来の改修をここへ直書きしない。

### 5.4 `sunscan_multibeam.py`
**責務**
- 1つの RawData / DB フォルダから複数 beam を順に解析する上位ラッパー
- beam ごとの `spectral_name`, `beam_id`, `hpbw_init_arcsec` などを解決する
- 各 beam に対して `sunscan_core` を呼ぶ
- `sunscan_report` を呼び、single-beam 互換の出力を beam ごとに作る
- multi-beam 独自の集約出力（all-beam CSV, fit input CSV など）を追加する

**重要条件**
- single-beam 互換出力は「追加」してよいが、「置換」してはいけない。

---

## 6. `SunScanAnalysisConfig` 正式項目一覧（確定）

`SunScanAnalysisConfig` は以下の下位設定を束ねる親設定とする。

### 6.1 `InputConfig`
- `rawdata_path: Path`
- `spectral_name: str`
- `azel_source: Literal["encoder", "altaz"]`
- `altaz_apply: Literal["none", "minus", "plus"]`
- `encoder_shift_sec: float`
- `encoder_vavg_sec: float`

### 6.2 `CalibrationConfig`
- `chopper_wheel: bool`
- `tamb_k: Optional[float]`
- `chopper_win_sec: float`  
  ※ backward compatibility 用。現行実装では実質 deprecated。
- `chopper_stat: Literal["median", "mean"]`  
  ※ backward compatibility 用。現行実装では実質 deprecated。

### 6.3 `RippleConfig`
- `enabled: bool`
- `preset: Literal["auto", "safe", "normal", "strong"]`
- `model: Literal["auto", "add", "mul"]`
- `target_hz: float`
- `search_hz: float`
- `bw_hz: Optional[float]`
- `max_harm: Optional[int]`
- `order: Optional[int]`
- `notch_passes: Optional[int]`
- `trend_win_sec: Optional[float]`
- `resample_dt_sec: Optional[float]`
- `eval_band_hz: Optional[float]`

### 6.4 `ProfileConfig`
- `profile_xlim_deg: float`

**注記**
- `profile_xlim_deg` は InputConfig と TrimConfig の両方に重複して置かない。
- 現行実装では、`profile_xlim_deg` は
  - trim の near-Sun 候補選定
  - ripple auto の速度推定範囲
  - derivative/profile 図の x 範囲
  に共通で使われている。
- よって **単一の `ProfileConfig.profile_xlim_deg` を single source of truth とする。**

### 6.5 `TrimConfig`
- `enabled: bool`
- `vfrac: float`
- `vmin: float`
- `gap_fill: int`
- `min_samples: int`
- `dominant_axis: bool`
- `ratio_min: float`
- `vpercentile: float`
- `steady_scan: bool`
- `use_on_only: bool`
- `xwin_factor: float`
- `cross_offset_max_deg: float`
- `speed_min_deg_s: float`
- `steady_cv_max: float`

### 6.6 `EdgeFitConfig`
- `enabled: bool`
- `strict_deriv: bool`
- `fit_win_deg: float`
- `fit_threshold: float`
- `hpbw_init_arcsec: float`

### 6.7 `ReportConfig`
- `outdir: Path`
- `debug_plot: bool`
- `edge_fit_plot_max_scans: int`
- `tag: Optional[str]`

**注記**
- `tag` の既定挙動は現行互換とし、未指定時は `rawdata_path.name` を使う。
- `edge_fit_plot_max_scans` は edge fit そのものではなく、出力ページング条件なので `ReportConfig` に属する。

### 6.8 `RuntimeConfig`
- `continue_on_error: bool`

---

## 7. beam ごとの上書き規則
multi-beam 側では、共通 `SunScanAnalysisConfig` に対して beam ごとの override を許可する。

### 7.1 v1 で beam ごとに上書きを許可する項目
- `InputConfig.spectral_name`
- `EdgeFitConfig.hpbw_init_arcsec`
- 将来的な `restfreq_hz` 由来の HPBW 初期値
- `ReportConfig.tag`
- `beam_id`, `stream_name` などの識別情報

### 7.2 v1 で beam ごとに上書きを原則禁止する項目
- trim
- ripple
- edge-fit threshold / window
- azel source
- encoder shift / smoothing
- chopper-wheel on/off
- profile_xlim_deg

**理由**
同一観測中に beam 間で前処理が異なると、beam 位置差にソフト差分が混入するため。

---

## 8. 出力ファイル仕様（single-beam 互換）

### 8.1 既存互換として維持するファイル名
- `sun_scan_summary_{tag}.csv`
- `sun_scan_derivative_fits_{tag}_pXXX.png`
- `summary_text_{tag}.png`  ※ `debug_plot=True` のときのみ
- `AZ_dxdy_{scan_id}_{tag}.png`  ※ `debug_plot=True`
- `EL_dxdy_{scan_id}_{tag}.png`  ※ `debug_plot=True`
- `AZ_mount_azel_{scan_id}_{tag}.png`  ※ `debug_plot=True`
- `EL_mount_azel_{scan_id}_{tag}.png`  ※ `debug_plot=True`

### 8.2 summary CSV の canonical 列順（完全版）
以下を canonical な preferred 列順とする。

1. `scan_id`
2. `rep_az_deg`
3. `rep_el_deg`
4. `center_az_deg`
5. `center_el_deg`
6. `hpbw_az_arcsec`
7. `hpbw_el_arcsec`
8. `sun_az_deg`
9. `sun_el_deg`
10. `speed_az_arcsec_s`
11. `speed_el_arcsec_s`
12. `fit_ok_az`
13. `fit_ok_el`
14. `az_error`
15. `el_error`
16. `n_az`
17. `n_az_used`
18. `n_el`
19. `n_el_used`
20. `az_left_amp`
21. `az_left_center_deg`
22. `az_left_sigma_deg`
23. `az_right_amp`
24. `az_right_center_deg`
25. `az_right_sigma_deg`
26. `el_low_amp`
27. `el_low_center_deg`
28. `el_low_sigma_deg`
29. `el_high_amp`
30. `el_high_center_deg`
31. `el_high_sigma_deg`
32. `ripple_applied_az`
33. `ripple_applied_el`
34. `az_track_t_unix`
35. `el_track_t_unix`
36. `az_track_az_deg`
37. `az_track_el_deg`
38. `az_track_main_offset_deg`
39. `az_track_cross_offset_deg`
40. `el_track_az_deg`
41. `el_track_el_deg`
42. `el_track_main_offset_deg`
43. `el_track_cross_offset_deg`
44. `data_tag`
45. `y_axis`
46. `azel_source`
47. `altaz_apply`
48. `encoder_shift_sec`
49. `encoder_vavg_sec`
50. `chopper_wheel`
51. `ripple_remove`
52. `ripple_preset`
53. `edge_fit_win_deg`
54. `edge_fit_threshold`
55. `hpbw_init_arcsec`
56. `trim_scan`
57. `profile_xlim_deg`

### 8.3 multi-beam 側での追加列
multi-beam では canonical 列集合の前または後ろに以下を追加してよい。

- `beam_id`
- `stream_name`
- `run_id`
- `rawdata_path`

ただし、**single-beam 互換 CSV も別途出力すること**。

---

## 9. 検証要件（リファクタリング受け入れ条件）
同じ RawData と同じ CLI 引数に対して、少なくとも以下が旧版と一致または実質一致すること。

1. `sun_scan_summary_{tag}.csv` の列名と列順
2. `center_az_deg`, `center_el_deg`, `hpbw_az_arcsec`, `hpbw_el_arcsec`, `rep_az_deg`, `rep_el_deg`
3. `fit_ok_az`, `fit_ok_el`
4. derivative fit PNG のページ数とファイル命名規則
5. `summary_text_{tag}.png` の生成条件
6. `continue_on_error=False` 時の debug plot 停止挙動
7. `continue_on_error=True` 時の debug plot 継続挙動
8. scan 単位例外が `az_error` / `el_error` に記録される挙動

---

## 10. 実装順序

### Phase 1: 純粋リファクタリング
- `sunscan_core.py` と `sunscan_report.py` を新設
- 現行 `sun_scan_v4_v5.py` から関数を機械的に移す
- `sun_scan_v4_v5.py` は import して呼ぶだけにする
- 挙動差分ゼロを確認する

### Phase 2: Config 導入
- `SunScanAnalysisConfig` を導入
- `main()` では argparse → config 変換のみ行う
- core は config を受け取る API を追加する
- `make_ripple_policy_from_args()` を wrapper 側へ寄せ、core は argparse 非依存にする

### Phase 3: multi-beam ラッパー追加
- `sunscan_multibeam.py` を追加
- beam ごとの override を適用する
- single-beam 互換出力 + multi-beam 集約出力を作る

---

## 11. 本件に関する最終結論
- 共通化は **fit だけでは不十分**。
- **ripple / trim / steady-scan selection / derivative / fit / report** までをまとめて共通コア化すべきである。
- ただし `sun_scan_v4_v5.py` の既存 CLI と既存出力は、少なくとも v1 では互換性保持対象として凍結する。
- multi-beam は `sun_scan` の上位ラッパーであり、別解析器ではない。
- `profile_xlim_deg` は1か所にだけ定義する。
- `continue_on_error` の意味は現行スコープから拡張しない。
- argparse 依存を core に固定しない。

