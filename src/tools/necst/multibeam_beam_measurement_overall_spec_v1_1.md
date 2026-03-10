# Multi-beam 用 beam 測定パッケージ 全体仕様書
## 統合版 v1.1

---

## 0. 目的と適用範囲

本仕様書は、太陽スキャン観測を用いて multi-beam 受信機の各 beam の太陽に対するオフセットを測定し、その結果を beam 幾何モデルとして整理し、最終的に SDFITS 書き出し系で利用可能な形式へ落とし込むための**全体設計**を定める。

本仕様の対象は次の 3 系統である。

1. 既存 `sun_scan_v4_v5.py` の single-beam 太陽スキャン解析系
2. 既存 `necst_v4_sdfits_converter.py` の SDFITS 変換系
3. 新規に追加する multi-beam 用 beam 測定パッケージ

本仕様の最優先原則は次のとおりである。

- **既存 `sun_scan_v4_v5.py` を壊さない**
- **既存 `necst_v4_sdfits_converter.py` を初期段階では改変しない**
- multi-beam 側で `sun_scan` 相当の解析ロジックを再実装しない
- 将来 `sun_scan` 側を改善した際、その改善を multi-beam 側へ最小作業で反映できる構造にする
- beam 幾何 fit の結果は converter 互換の出力形式で保存できるようにする

本仕様書は、実装前の最終設計書であり、**まず既存 converter と既存 sun_scan が従来どおり動作することを確認し、その後に multi-beam 用 beam 測定パッケージの検証へ進む**前提で記述する。

---

## 1. システム全体像

システム全体は、論理的に次の 4 段階からなる。

### Phase 0: 既存系のベースライン確認

- `sun_scan_v4_v5.py` が従来どおり動作することを確認する
- `necst_v4_sdfits_converter.py` が従来どおり動作することを確認する
- 比較用の基準 RawData と基準出力を保存する

### Phase 1: `sun_scan` の共通コア化

- `sun_scan_v4_v5.py` の解析ロジックを `sunscan_*` モジュールへ切り出す
- ただし CLI・既定値・出力形式・数値結果は維持する
- この段階では新機能を足さない

### Phase 2: multi-beam 用 beam 測定パッケージの追加

- 1 つの RawData / DB フォルダに含まれる複数 beam を解析する
- 各解析対象 stream について、`sun_scan` 相当の解析を同一コアで実行する
- beam ごとの太陽に対するオフセットを取得する
- stream ごとの single-beam 互換出力と、all-stream / all-beam 集約出力を生成する

### Phase 3: beam 幾何 fit と converter 互換出力

- all-stream 集約出力を入力として beam 幾何 fit を行う
- center-beam model と virtual-center model を評価する
- converter 互換の beam model TOML を出力する
- 必要なら virtual center 補助ファイルも出力する

**重要**: 初期リリース範囲では、converter 本体の内部実装変更は必須条件ではない。converter 互換の出力を行えればよく、converter への共通部移植は後続段階の検討項目とする。

---

## 2. 既存系についての基本方針

### 2.1 `sun_scan_v4_v5.py` について

`sun_scan_v4_v5.py` は single-beam 太陽スキャン解析の**標準フロントエンド**として扱う。既存利用者にとって見えるインターフェースは維持する。

変更禁止対象は次のとおりである。

1. CLI 引数名
2. CLI の既定値
3. `tag` 未指定時の既定規則（実質 `rawdata_path.name`）
4. summary CSV の列名・列順
5. derivative-fit PNG の基本命名規則
6. standard output の基本構成
7. `continue_on_error` の**現行実装どおりの挙動**
8. オフセット定義
   - `daz = wrap(Az_true - Az_sun) * cos(El_sun)`
   - `d_el = El_true - El_sun`

補足:

- `--chopper-win-sec` と `--chopper-stat` は現行でも後方互換用に残っているが、現在のアルゴリズムでは実質無効である。共通化後もこの性質を維持する。
- `--debug-plot` は `AZ_dxdy_*`, `EL_dxdy_*`, `AZ_mount_azel_*`, `EL_mount_azel_*`, `summary_text_*` の出力を制御する。`sun_scan_derivative_fits_*_pXXX.png` は `edge_fit` が有効であれば `debug_plot` とは独立に出力される。
- 現行 `continue_on_error` は「すべての例外を黙って継続する」意味ではない。少なくとも debug plot ループでは効くが、per-scan fit 失敗は内部で `az_error` / `el_error` に落とされ、top-level edge-fit block は警告を出して継続する。この**現状の意味をそのまま維持**する。

### 2.2 `necst_v4_sdfits_converter.py` について

converter は初期段階では**現状維持**とする。multi-beam 用 beam 測定パッケージは converter 本体を変更しなくても利用できる構造にする。

本仕様における converter との関係は次の 2 点に限る。

1. beam 幾何 fit の結果を converter 互換 TOML として出力できること
2. spectrometer / beam の識別情報を converter と矛盾しない形で扱うこと

したがって、初期段階では

- converter の CLI は変更しない
- converter の既存出力は変更しない
- converter 内部へ新しい fit ロジックを入れない

---

## 3. multi-beam beam 測定パッケージの設計原則

### 3.1 single source of truth

太陽スキャン解析の single source of truth は `sun_scan` 系共通コアである。multi-beam 側はこのコアを解析対象 stream ごとに呼ぶ上位ラッパーであり、独自の trim / ripple / derivative / edge-fit ロジックを持たない。

### 3.2 パイプライン全体を共通化する

共通化対象は DB 読み込みだけでは不十分である。以下を共有サブシステムに含める。

- RawData / DB 読み込み
- Az/El source 選択
- Ta* 化
- trim / steady-scan 選別
- ripple 設定解決と ripple 除去
- derivative 計算
- edge fit
- per-scan summary 生成
- summary CSV / PNG / text 出力

ただし、これらすべてを `sunscan_core.py` 一枚に入れるという意味ではない。I/O、数値処理、report は別モジュールに分けてよいが、**解析ロジックとしては single source of truth を維持**する。

### 3.3 multi-beam 側は stream 列挙と beam 幾何差分のみを担当する

multi-beam 側が追加で持つ責務は次に限定する。

- 同一 RawData 内の複数解析対象 stream の列挙
- stream ごとの最小限の override 解決
- stream ごとの解析実行
- all-stream 集約表の生成
- beam 幾何 fit の実行
- converter 互換出力の保存

### 3.4 mixed-frequency multi-beam を正規ケースとする

multi-beam 受信機では、中心 beam と周辺 beam で周波数が異なる可能性を標準ケースとして扱う。従って beam ごとに HPBW 初期値を持てる設計にする。

### 3.5 解析の基本単位と物理 beam の区別

v1 では次を明確に区別する。

- **解析基本単位**: spectrometer 設定の `name` で一意に識別される解析対象 stream
- **物理 beam**: `beam_id` で識別される受信ビーム

同じ `beam_id` に複数の stream（例: 異なる偏波）が紐づくことはあり得る。抽出段階では stream 単位で解析する。beam 幾何 fit は beam_id 単位の物理 beam を対象とする。

---

## 4. モジュール構成

### 4.1 `sunscan_config.py`

役割:

- 解析設定の dataclass 群を定義する
- CLI 非依存であること
- single-beam と multi-beam の両方が共通に使う設定表現を与えること

### 4.2 `sunscan_io.py`

役割:

- RawData / DB から必要テーブルを読み出す
- 解析前の**生 bundle** を組み立てる
- ここで返すのは raw table 群、補間前 series、metadata、stream 識別情報などであり、`daz` / `d_el` を含む最終 `DataFrame` はまだ作らない
- fit や report は行わない

補足: `build_dataframe()` 相当の最終整形は `sunscan_core.py` に置く。これにより I/O と数値処理の境界を明確にする。

### 4.3 `sunscan_core.py`

役割:

- `sun_scan` 相当の数値処理本体を持つ
- raw bundle と解析設定から `DataFrame`, `az_scans`, `el_scans` を構築する
- trim、ripple、derivative、edge fit、per-scan summary を実行する
- single-beam / multi-beam の両方が直接使用する

### 4.4 `sunscan_report.py`

役割:

- summary text を整形する
- summary CSV を出力する
- derivative-fit PNG を出力する
- multi-beam 用の追加集約出力も補助できる

### 4.5 `sun_scan_v4_v5.py`

役割:

- 従来どおりの CLI を提供する thin wrapper
- `argparse` から `SunScanAnalysisConfig` を生成する
- `sunscan_io.py` / `sunscan_core.py` / `sunscan_report.py` を呼ぶ
- 既存互換のファイル名・標準出力・例外処理フローを維持する

### 4.6 `sunscan_multibeam.py`

役割:

- multi-beam 用の最上位 wrapper
- サブコマンドまたは内部呼び出しにより、extract / fit を統括する
- `sun_scan` 共通コアを解析対象 stream ごとに反復適用する

### 4.7 `sunscan_extract_multibeam.py`

役割:

- 1 つの RawData / DB フォルダに対し、解析対象 stream ごとに `sun_scan` 相当解析を実行する
- stream ごとの single-beam 互換出力を生成する
- all-stream 集約 scan summary を生成する

備考: `sunscan_multibeam.py extract` の thin wrapper として実装してよい。

### 4.8 `sunscan_fit_multibeam.py`

役割:

- all-stream 集約 scan summary を読み込み、beam 幾何 fit を行う
- center-beam model / virtual-center model を評価する
- beam fit 結果表・nuisance shift 表・converter 互換 TOML を出力する

備考: `sunscan_multibeam.py fit` の thin wrapper として実装してよい。

### 4.9 `beam_geometry.py`

役割:

- beam 幾何 fit で使う純粋な幾何計算関数を持つ
- 回転表現、offset 回転、converter 互換 export で共通に用いる
- 初期段階では converter 本体から直接利用させなくてもよい

---

## 5. 設定仕様

### 5.1 `SunScanAnalysisConfig` の正式構造

`SunScanAnalysisConfig` は以下の下位設定を束ねる。

#### 5.1.1 `InputConfig`

- `rawdata_path: Path`
- `spectral_name: str`
- `azel_source: Literal["encoder", "altaz"]`
- `altaz_apply: Literal["none", "minus", "plus"]`
- `encoder_shift_sec: float`
- `encoder_vavg_sec: float`

#### 5.1.2 `CalibrationConfig`

- `chopper_wheel: bool`
- `tamb_k: Optional[float]`
- `chopper_win_sec: float`  # 後方互換のため保持するが現アルゴリズムでは無効
- `chopper_stat: Literal["median", "mean"]`  # 後方互換のため保持するが現アルゴリズムでは無効

#### 5.1.3 `RippleConfig`

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

#### 5.1.4 `ProfileConfig`

- `profile_xlim_deg: float`

補足: `profile_xlim_deg` は表示範囲だけでなく、steady-scan 選別や ripple 自動判定の near-Sun 判定にも効くため、単なる plot option ではない。

#### 5.1.5 `TrimConfig`

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

#### 5.1.6 `EdgeFitConfig`

- `enabled: bool`
- `strict_deriv: bool`
- `fit_win_deg: float`
- `fit_threshold: float`
- `hpbw_init_arcsec: float`

#### 5.1.7 `ReportConfig`

- `outdir: Path`
- `debug_plot: bool`
- `edge_fit_plot_max_scans: int`
- `tag: Optional[str]`

補足: 現行 `sun_scan` には `tag` の CLI はない。`tag` は内部の report 命名用フィールドであり、single-beam wrapper では既定で `rawdata_path.name` に解決する。

#### 5.1.8 `RuntimeConfig`

- `continue_on_error: bool`

### 5.2 stream ごとの override 規則

#### 5.2.1 v1 で stream ごとに上書きを許可する項目

- `InputConfig.spectral_name`
- `EdgeFitConfig.hpbw_init_arcsec`
- `ReportConfig.tag`（multi-beam 内部でのみ使用）
- `beam_id`
- `stream_name`
- `restfreq_hz`
- 将来の HPBW 自動推定に必要な stream 固有補助情報

#### 5.2.2 v1 で stream ごとの上書きを原則禁止する項目

- trim 関連しきい値
- ripple preset / しきい値
- edge fit の threshold / window
- azel source
- encoder shift / smoothing

理由: beam 幾何差と解析差を混ぜないため。

### 5.3 mixed-frequency beam の HPBW 初期値

HPBW 初期値の決定優先順位は次のとおりとする。

1. stream 個別指定 `hpbw_init_arcsec`
2. `restfreq_hz` と望遠鏡口径 `dish_diameter_m` からの自動推定
3. 共通デフォルト値

自動推定式は v1 では次を標準とする。

`HPBW_init_arcsec = hpbw_factor * (c / restfreq_hz) / dish_diameter_m * 206265`

ここで `hpbw_factor` は設定値とする。

---

## 6. 識別子と設定ファイル

### 6.1 解析対象 stream と beam 識別子

- 解析対象 stream の主識別子は spectrometer 設定の `name` とする。これは v1 で一意でなければならない。
- 物理 beam の主識別子は `beam_id` とする。
- `fdnum` は beam_id 未指定時の既定割当規則にのみ用いる。

### 6.2 spectrometer 設定

multi-beam 解析では、stream ごとの設定を TOML で与えられるようにする。最小限、次の情報を持てること。

- `name`
- `db_stream_name`
- `beam_id`
- `polariza`
- `restfreq_hz`
- beam model export に必要な格納先

### 6.3 `polariza` の正式表記

仕様書上、`polariza` は次の正式コードのみを許容する。

- `RR`, `LL`, `RL`, `LR`
- `XX`, `YY`, `XY`, `YX`

`LCP`, `RCP` のような表現は仕様上の正式コードとしない。

### 6.4 同一 `beam_id` に複数 stream がある場合の v1 規則

同一 `beam_id` に複数の stream（例: 異なる偏波）が存在し得る。v1 では次の規則とする。

1. **抽出段階**では各 stream を独立に解析してよい。
2. **beam 幾何 fit 段階**では、1 つの `beam_id` に対して同時に複数 stream を無条件に混ぜない。
3. fit に用いる解析対象は、原則として `beam_id` ごとに高々 1 つの primary stream とする。
4. 同一 `beam_id` に複数 stream があり、primary の指定が無い場合は fit をエラーにする。

primary stream の指定方法は、v1 では次のいずれかを許容する。

- spectrometer-config に `beam_fit_use = true/false` を追加する
- fit CLI で `--fit-stream-names` または同等の selector を与える

---

## 7. 抽出段階の仕様

### 7.1 解析の基本単位

1 つの RawData / DB フォルダ中に複数 beam の spectral stream が含まれることを前提とする。抽出段階では、各**解析対象 stream** ごとに `sun_scan` 相当解析を実行する。

### 7.2 1 stream の解析で出すもの

各 stream ごとに次を出力可能とする。

- single-beam 互換 summary CSV
- single-beam 互換 derivative-fit PNG
- single-beam 互換 summary text
- sample-level 補助表（任意）

### 7.3 single-beam 互換出力の衝突回避

multi-beam では複数 stream を同一 run で処理するため、single-beam 互換出力を同一 directory に直接置くとファイル名が衝突する。v1 では次のどちらかを必須とする。

#### 推奨方式 A: stream ごとの subdirectory

- `outdir/per_stream/<stream_name>/sun_scan_summary_<base_tag>.csv`
- `outdir/per_stream/<stream_name>/sun_scan_derivative_fits_<base_tag>_pXXX.png`
- `outdir/per_stream/<stream_name>/AZ_dxdy_*_<base_tag>.png`
- `outdir/per_stream/<stream_name>/EL_dxdy_*_<base_tag>.png`
- `outdir/per_stream/<stream_name>/summary_text_<base_tag>.png`

この方式では、各 subdirectory 内のファイル名は single-beam と同一規則を保てる。

#### 許容方式 B: tag に stream 識別子を付加

- `base_tag = rawdata_path.name`
- `resolved_tag = <base_tag>__<stream_name>` など

ただし、single-beam 互換性の観点では方式 A を推奨する。

### 7.4 all-stream 集約 scan summary

multi-beam 用の主たる抽出出力は、1 行 = 1 run × 1 stream × 1 scan の scan summary table とする。

必須列は少なくとも次を含む。

#### 識別列
- `run_id`
- `tag`
- `stream_name`
- `beam_id`
- `scan_id`

#### 幾何列
- `rep_az_deg`
- `rep_el_deg`
- `center_az_deg`
- `center_el_deg`
- `x_arcsec = center_az_deg * 3600`
- `y_arcsec = center_el_deg * 3600`

#### 品質列
- `hpbw_az_arcsec`
- `hpbw_el_arcsec`
- `fit_ok_az`
- `fit_ok_el`

#### provenance 列
- `source_db_path`
- `spectral_name`
- `restfreq_hz`
- `hpbw_init_arcsec`
- `analysis_config_digest`
- `polariza`
- `fdnum`（利用可能なら）
- `ifnum`（利用可能なら）
- `plnum`（利用可能なら）
- `sampler`（利用可能なら）

### 7.5 canonical な列順

single-beam summary CSV の canonical 列順は維持する。all-stream 集約表では、その前に識別列を追加し、その後に single-beam canonical 列を続ける。

### 7.6 `run_id` の既定

`run_id` は、既定では RawData ディレクトリ名を用いる。複数回観測を後でまとめる場合は、明示指定で上書きできる。

---

## 8. beam 幾何 fit の仕様

### 8.1 fit の目的

beam 幾何 fit は、抽出段階で得られた各 beam の太陽に対するオフセットから、beam 配置と EL 依存回転を推定するために行う。

### 8.2 fit に用いる観測量

各 beam `b`、scan `s` について

- `x_obs(b,s) = 3600 * center_az_deg`
- `y_obs(b,s) = 3600 * center_el_deg`
- `E_s = rep_el_deg`

を用いる。

### 8.3 mixed-frequency と fit

mixed-frequency 系では beam ごとに HPBW が異なってもよいが、beam 幾何 fit に用いる位置量の定義は共通である。周波数差は主に edge fit の初期値と品質評価に効く。

### 8.4 center-beam model

中心 beam がある場合は、中心 beam との差を用いて相対配置を求める。

- `x_rel_obs(b,s) = x_obs(b,s) - x_obs(c,s)`
- `y_rel_obs(b,s) = y_obs(b,s) - y_obs(c,s)`

共通 pointing shift は差分で消える。

### 8.5 virtual-center model

中心 beam がない場合、または中心 beam があっても比較のため、仮想中心モデルを評価する。各 scan または各 run ごとの nuisance shift を導入して beam 間相対配置を求める。

### 8.6 光学軸ずれがある場合

中心 beam がなく、かつ光学軸の pointing ずれがある場合、Sun scan 単独では絶対中心は一意に決まらない。したがって v1 では

- beam 間相対配置
- scan / run ごとの nuisance shift

のみを推定対象とし、絶対中心補正は別較正とする。

### 8.7 `reference_angle_deg` と `dewar_angle_deg`

現行表現では

`theta(EL) = rotation_sign * (EL - reference_angle_deg) + dewar_angle_deg`

であるが、`reference_angle_deg` と `dewar_angle_deg` は独立 fit パラメータとして扱わない。さらに、共通位相と base offset ベクトルの向きも縮退するため、v1 の内部 fit では独立な `rotation_phase_deg` または固定ゲージ表現を用いる。

export 時にのみ converter 互換表現へ変換する。

### 8.8 `rotation_sign`

`rotation_sign` は v1 では離散候補 `+1` または `-1` を評価し、より良い方を採用する。

### 8.9 `rotation_mode`

v1 では `rotation_mode = "elevation"` を標準とし、`"none"` を許容する。高次モデルや一般角依存は v1 の範囲外とする。

### 8.10 fit の入力単位

beam 幾何 fit の primary 入力単位は `beam_id` とする。ただし、抽出表は stream 単位なので、fit 前に `beam_id` ごとの primary stream 選別を済ませる必要がある。

### 8.11 fit の出力

beam 幾何 fit の出力は次を含む。

- beam fit result table
- nuisance shift result table
- converter 互換 beam model TOML
- virtual center 補助 TOML（必要時）
- human-readable summary

---

## 9. converter 互換出力仕様

### 9.1 目的

beam 幾何 fit の結果を、converter が理解できる beam model 形式へ落とし込む。

### 9.2 export の原則

- export は converter 互換であること
- 初期段階では converter 内部変更を前提にしないこと
- fit 内部の gauge と export gauge を混同しないこと

### 9.3 export に含む代表キー

- `az_offset_arcsec`
- `el_offset_arcsec`
- `rotation_mode`
- `reference_angle_deg`
- `rotation_sign`
- `dewar_angle_deg`

ただし、`reference_angle_deg` / `dewar_angle_deg` は fit 内部の独立推定量ではなく、export 用表現である。

### 9.4 virtual center 補助ファイル

中心 beam がない場合や絶対中心補正を別扱いする場合は、virtual center 補助ファイルを別出力とする。これは converter 本体へ直接入力することを必須としない。

---

## 10. CLI 構成

### 10.1 `sun_scan_v4_v5.py`

- 従来どおりの CLI を維持する
- 引数・既定値・出力形式を変えない

### 10.2 `sunscan_multibeam.py`

以下のいずれかの設計を許容する。

#### 案 A: サブコマンド方式
- `sunscan_multibeam.py extract`
- `sunscan_multibeam.py fit`

#### 案 B: 個別スクリプト thin wrapper 方式
- `sunscan_extract_multibeam.py`
- `sunscan_fit_multibeam.py`

v1 では、どちらで実装してもよい。ただし内部ロジックの重複は禁止する。

### 10.3 抽出 CLI の必須引数

- `rawdata_path`
- `--spectrometer-config`

### 10.4 fit CLI の必須引数

- `input_paths...`
- `--spectrometer-config`

### 10.5 fit CLI の追加必須機能

同一 `beam_id` に複数 stream が存在する可能性に対応するため、fit CLI には少なくとも次のどちらかを持たせる。

- `--fit-stream-names` のような stream selector
- config 側の `beam_fit_use` を尊重する機能

---

## 11. 実装順序

### Step 1: 既存系のベースライン固定

- 既存 `sun_scan_v4_v5.py` の基準出力を保存
- 既存 `necst_v4_sdfits_converter.py` の基準出力を保存

### Step 2: `sun_scan` 共通コア化

- `sunscan_config.py`
- `sunscan_io.py`
- `sunscan_core.py`
- `sunscan_report.py`

を作成する。

この段階では新機能を追加しない。

### Step 3: `sun_scan_v4_v5.py` を thin wrapper 化

- 外見を維持したまま内部呼び出し先だけ差し替える
- 同一 RawData で旧版と数値一致を確認する

### Step 4: multi-beam 抽出段階を追加

- `sunscan_multibeam.py`
- `sunscan_extract_multibeam.py`

を追加し、stream ごとに single-beam 共通コアを反復適用する。

### Step 5: beam 幾何 fit を追加

- `beam_geometry.py`
- `sunscan_fit_multibeam.py`

を追加する。

### Step 6: converter 互換出力を確認

- 出力 TOML が converter で読めることを確認する
- この段階でも converter 本体の改変は必須としない

### Step 7: 必要なら後続で converter 共通部を検討

converter と beam 測定パッケージの間で純粋関数を共有したい場合は、**Phase 0〜6 の等価性確認後**に別タスクとして検討する。

---

## 12. 回帰試験

### 12.1 `sun_scan` の回帰試験

以下を旧版と新ラッパー版で比較する。

1. summary CSV の数値一致
2. `az_scans` / `el_scans` の scan 分割一致
3. `center_az_deg`, `center_el_deg`, `hpbw_*`, `fit_ok_*` の一致
4. derivative-fit PNG のページングと命名規則
5. `continue_on_error` の現行挙動一致
6. `debug_plot` の on/off による生成ファイル集合の一致

### 12.2 multi-beam 抽出の回帰試験

同じ stream / spectral を個別に `sun_scan` 解析した結果と、multi-beam 経由でその stream を解析した結果が一致すること。

### 12.3 beam 幾何 fit の検証

- center-beam model と virtual-center model の両方を実行できること
- mixed-frequency 系で beam ごとの HPBW 初期値が正しく使われること
- `beam_id` 主識別が正しく働くこと
- 同一 `beam_id` に複数 stream がある場合に、primary stream の選別が強制されること

### 12.4 converter 互換検証

- export TOML が converter 互換形式であること
- spectrometer / beam 識別情報が矛盾しないこと

---

## 13. v1 の範囲内 / 範囲外

### 13.1 v1 の範囲内

- `sun_scan` 共通コア化
- multi-beam 抽出
- beam 幾何 fit
- converter 互換 export
- mixed-frequency beam の HPBW 初期値対応
- center-beam / virtual-center の両モデル評価
- 同一 `beam_id` に複数 stream がある場合の primary stream 選別

### 13.2 v1 の範囲外

- converter 本体の本格的な内部共通化
- 高次の回転モデル
- beam ごとに異なる trim / ripple 戦略
- Sun scan 単独からの絶対光学軸推定
- 複数偏波 stream の自動統合平均を一般解として実装すること
- GUI

---

## 14. 最終的な要約

本件の本質は、**multi-beam 用 beam 測定パッケージを新規に別実装することではなく、`sun_scan` の解析本体を single source of truth として切り出し、それを single-beam / multi-beam の両方が共有する構造へ整理すること**にある。

そのうえで、stream ごとの解析結果を all-stream 集約し、beam 幾何 fit を別段階として実行し、converter 互換形式で保存する。

実装の最初の関門は新機能ではなく、**既存 `sun_scan` と既存 converter が従来どおり問題なく動くことをベースラインとして固定すること**である。

この方針により、

- `sun_scan` を壊さない
- converter を初期段階では壊さない
- multi-beam 測定パッケージを追加できる
- 将来の `sun_scan` 改善を multi-beam 側へ自然に反映できる

という 4 条件を同時に満たす。
