# Multi-beam 用 beam 測定パッケージ 説明書 v4-02  
## sunscan command reference・全 CLI パラメータ説明

## 1. 本書の役割

本書は multi-beam package 側の CLI 引数を、**現行実装ベースで**詳細に説明するための分冊です。

対象コマンド:

- `check_spectrometer_config`
- `sunscan_singlebeam`
- `sunscan_extract_multibeam`
- `sunscan_fit_multibeam`
- `sunscan_make_pseudo_multibeam`
- `sunscan_multibeam`

converter そのものの全 CLI は本書の対象外です。  
ただし、converter-compatible TOML との関係が深い部分は必要に応じて触れます。

---

## 2. 実行形式について

環境によっては console script が入っている場合と、module 実行が必要な場合があります。  
以下では両方の形を意識してください。

例:

```bash
sunscan_extract_multibeam ...
```

または

```bash
python -m multibeam_beam_measurement.sunscan_extract_multibeam ...
```

---

## 3. `check_spectrometer_config`

## 3.1 目的

- spectrometer config TOML を検査する
- stream table を表示する
- primary stream 解決を確認する
- warning を表示する

## 3.2 引数

### 位置引数 `spectrometer_config`
意味:
- spectrometer config TOML のパス

### `--stream-name`
型:
- repeatable string

意味:
- primary stream の explicit selection
- 同じ beam に複数 stream があるとき、fit 代表候補の解決確認を明示的に行いたい場合に使う

### `--out-csv`
型:
- string

意味:
- stream table CSV の出力先

## 3.3 出力

- 標準出力に stream table
- `primary_streams`
- `duplicate_beam_ids`
- warnings
- 必要なら `--out-csv` で CSV

---

## 4. `sunscan_singlebeam`

## 4.1 目的

- 1 本の stream を single-beam として解析する
- multi 用 spectrometer config TOML をそのまま流用できる

## 4.2 基本的な使い方

### config を使わない legacy 的な形
```bash
sunscan_singlebeam RAWDATA --spectral-name xffts-board1
```

### spectrometer config TOML から 1 本選ぶ
```bash
sunscan_singlebeam RAWDATA \
  --spectrometer-config beams.toml \
  --stream-name center230_xx
```

---

## 4.3 引数の完全説明

### 位置引数 `rawdata`
意味:
- NECST DB を含む RawData directory

### `--outdir`
既定:
- `.`

意味:
- 出力先 directory

### `--spectrometer-config`
意味:
- converter-compatible spectrometer config TOML
- 指定すると `[global]` の shared setting を読み、stream も選択できる

### `--stream-name`
意味:
- `--spectrometer-config` 内の logical stream 名で明示選択する

重要:
- 明示指定すると、`use_for_sunscan=false` の stream でも one-off override として選択できる

### `--db-namespace`
既定:
- `necst`

意味:
- DB namespace

優先:
- CLI > `[global].db_namespace` > default

### `--telescope`
既定:
- `OMU1P85M`

意味:
- telescope 名

優先:
- CLI > `[global].telescope` > default

### `--tel-loaddata`
既定:
- `OMU1p85m`

意味:
- telescope legacy/load name

優先:
- CLI > `[global].tel_loaddata` > default

### `--planet`
既定:
- `sun`

意味:
- 対象天体名

優先:
- CLI > `[global].planet` > default

### `--spectral-name`
既定:
- `xffts-board1`

意味:
- config を使わないときの spectral stream 名
- config を使う場合も、`--stream-name` を指定しないときの補助選択キーになる

### `--azel-source`
選択肢:
- `encoder`
- `altaz`

意味:
- Az/El の主たる元データをどちらにするか

### `--altaz-apply`
選択肢:
- `none`
- `minus`
- `plus`

意味:
- altaz correction の適用方法
- 詳細な物理解釈は観測系依存だが、single / extract で共通に扱う

### `--spectrometer-time-offset-sec`
既定:
- `0.0`

意味:
- spectrometer timestamp への補正秒

優先:
- CLI > `[global].spectrometer_time_offset_sec`

### `--encoder-shift-sec`
既定:
- `0.0`

意味:
- encoder Az / El 共通補正秒

優先:
- CLI > `[global].encoder_shift_sec`

### `--encoder-az-time-offset-sec`
既定:
- `0.0`

意味:
- encoder Az 追加補正秒

優先:
- CLI > `[global].encoder_az_time_offset_sec`

### `--encoder-el-time-offset-sec`
既定:
- `0.0`

意味:
- encoder El 追加補正秒

優先:
- CLI > `[global].encoder_el_time_offset_sec`

### `--encoder-vavg-sec`
既定:
- `0.0`

意味:
- encoder 平滑化 / 平均化に関する時間窓

### `--no-chopper-wheel`
既定:
- chopper wheel 使用 (`true`)

意味:
- 指定すると chopper wheel 補正を使わない
- その場合、`Ta*` ではなく raw power 系を使う解析になる

### `--tamb-k`
意味:
- 周囲温度の手動指定
- chopper wheel 系の fallback として使う

### `--chopper-win-sec`
既定:
- `5.0`

意味:
- HOT / OFF / ON の統計窓長

### `--chopper-stat`
選択肢:
- `median`
- `mean`

意味:
- chopper wheel 統計量

### `--profile-xlim-deg`
既定:
- `1.0`

意味:
- profile plot や fit 用に使う典型的な横軸範囲

### `--ripple-no-remove`
既定:
- ripple remove 有効

意味:
- ripple 除去を無効化

### `--ripple-preset`
選択肢:
- `auto`
- `safe`
- `normal`
- `strong`

意味:
- ripple 除去のプリセット

### `--ripple-model`
選択肢:
- `auto`
- `add`
- `mul`

意味:
- ripple モデルの形式

### `--ripple-target-hz`
既定:
- `1.2`

意味:
- 主対象 ripple 周波数

### `--ripple-search-hz`
既定:
- `0.3`

意味:
- ripple 周波数探索幅

### `--ripple-bw-hz`
意味:
- notch / band の幅

### `--ripple-max-harm`
意味:
- 高調波の最大次数

### `--ripple-order`
意味:
- ripple モデルの次数

### `--ripple-notch-pass`
意味:
- notch pass 回数

### `--ripple-trend-win-sec`
意味:
- trend 推定窓長

### `--ripple-resample-dt-sec`
意味:
- ripple 解析用の resample 間隔

### `--ripple-eval-band-hz`
意味:
- ripple 評価帯域

### `--no-edge-fit`
既定:
- edge fit 有効

意味:
- limb edge fit を無効化

### `--edge-fit-win-deg`
既定:
- `0.15`

意味:
- edge fit に使う窓幅

### `--edge-fit-threshold`
既定:
- `0.20`

意味:
- edge 選択 threshold

### `--hpbw-init-arcsec`
既定:
- `324.0`

意味:
- HPBW 初期値

重要:
- config から stream を選んだ場合、CLI で明示しなければ restfreq と dish diameter から自動推定で上書きされる

### `--edge-fit-plot-max-scans`
既定:
- `3`

意味:
- derivative fit plot の 1 ページあたり最大 scan 数

### `--no-trim-scan`
既定:
- trim 有効

意味:
- scan trimming を無効化

### `--trim-vfrac`
既定:
- `0.20`

意味:
- trimming 関連の速度閾値比

### `--trim-vmin`
既定:
- `1e-4`

意味:
- trimming 関連の最小速度閾値

### `--trim-gap`
既定:
- `10`

意味:
- gap fill 長

### `--trim-min-samples`
既定:
- `100`

意味:
- trimming に必要な最少サンプル数

### `--trim-dominant-axis`
既定:
- 有効

意味:
- dominant axis を使って trimming する

### `--trim-no-dominant-axis`
意味:
- dominant axis trimming を使わない

### `--trim-axis-ratio-min`
既定:
- `3.0`

意味:
- dominant axis とみなすための比率

### `--trim-vpercentile`
既定:
- `95.0`

意味:
- 速度 percentilե 判定閾値

### `--trim-no-steady-scan`
既定:
- steady scan 条件を使う

意味:
- steady scan 条件を無効化

### `--trim-include-hotoff`
既定:
- ON のみ使用 (`trim_use_on_only = true`)

意味:
- HOT / OFF も trimming 判定に含める

### `--trim-scan-speed-min-arcsec`
既定:
- `20.0`

意味:
- scan speed の下限 [arcsec/s]

### `--trim-xwin-factor`
既定:
- `1.2`

意味:
- trimming 横窓の factor

### `--trim-cross-offset-max-deg`
既定:
- `0.5`

意味:
- 許容する cross-offset 最大値 [deg]

### `--trim-steady-cv-max`
既定:
- `0.8`

意味:
- steady scan 判定の CV 上限

### `--strict-deriv`
既定:
- 有効

意味:
- derivative fit を厳密条件で行う

### `--no-strict-deriv`
意味:
- derivative fit の strict 条件を緩める

### `--continue-on-error`
意味:
- 一部の解析過程で例外が出ても継続を試みる

### `--debug-plot`
意味:
- debug plot を保存する

### `--dish-diameter-m`
既定:
- `1.85`

意味:
- dish diameter [m]
- HPBW 初期値自動推定にも使う

### `--hpbw-factor`
既定:
- `1.2`

意味:
- HPBW 推定係数

---

## 5. `sunscan_extract_multibeam`

## 5.1 目的

- config 内の複数 stream に対して single-beam 解析を実行
- 全 stream summary / manifest / snapshot を出力

## 5.2 基本例

```bash
sunscan_extract_multibeam RAWDATA \
  --spectrometer-config beams.toml \
  --outdir out_extract
```

## 5.3 引数

singlebeam と共通の解析パラメータが多数あります。  
相違点と multi 特有の意味を中心に説明します。

### 位置引数 `rawdata`
- RawData directory

### `--spectrometer-config`
- 必須
- converter-compatible spectrometer config TOML

### `--outdir`
- 出力先

### `--run-id`
意味:
- summary CSV に書く `run_id` を手動上書き

### `--db-namespace`
### `--telescope`
### `--tel-loaddata`
### `--planet`
意味:
- singlebeam と同じ
- `[global]` からも読める

### `--stream-name`
型:
- repeatable

意味:
- 明示的に対象 stream を絞る

重要:
- 明示指定した場合、usage flags による自動除外よりも explicit selection が優先される

### `--azel-source` ～ `--hpbw-factor`
意味:
- singlebeam と同じ
- extract は内部で各 stream ごとに singlebeam 相当解析を呼ぶため

## 5.4 extract 固有の current behavior

### usage flags による自動除外
明示 `--stream-name` が無い通常実行では、`use_for_sunscan=false` の stream は自動除外されます。

### skipped stream の記録
自動除外された stream は manifest に `status="skipped"` として残ります。

### per-stream 出力
`outdir/per_stream/<stream_name>/` に各 stream 個別出力を作ります。

### summary の先頭列
all-stream summary では主に次が先頭寄りに並びます。

- `run_id`
- `tag`
- `stream_name`
- `beam_id`
- `scan_id`
- `x_arcsec`
- `y_arcsec`
- `source_db_path`
- `spectral_name`
- `db_stream_name`
- `restfreq_hz`
- `analysis_config_digest`
- `polariza`
- `fdnum`
- `ifnum`
- `plnum`
- `sampler`

---

## 6. `sunscan_fit_multibeam`

## 6.1 目的

- extract summary CSV から beam geometry を fit

## 6.2 基本例

```bash
sunscan_fit_multibeam \
  out_extract/sunscan_multibeam_scan_summary_RUN.csv \
  --spectrometer-config beams.toml \
  --center-beam-id B00 \
  --model both \
  --outdir out_fit
```

## 6.3 引数

### 位置引数 `summary_csv`
型:
- 1 個以上

意味:
- all-stream summary CSV
- 複数ファイルを同時に与えると concat して fit する

### `--spectrometer-config`
- 必須
- primary stream 解決と beam model TOML 出力に必要

### `--outdir`
- 出力先

### `--center-beam-id`
意味:
- center-beam model の中心 beam ID

注意:
- `model=center_beam` では必須
- `model=both` では指定があれば `center_beam + virtual_center` の両方
- 指定が無ければ `both` でも実質 `virtual_center` だけになる

### `--fit-stream-name`
型:
- repeatable

意味:
- fit に使う stream を明示選択
- primary stream 解決を explicit override したいときに使う

### `--model`
選択肢:
- `both`
- `center_beam`
- `virtual_center`

意味:
- fit model の選択

### `--reference-angle-deg`
意味:
- fit 時の reference elevation angle を手動上書き

### `--sigma-clip`
既定:
- `4.5`

意味:
- residual clipping の sigma threshold
- `<=0` なら clipping 無効

### `--clip-iters`
既定:
- `2`

意味:
- sigma clipping 反復回数

### `--min-points-per-beam`
既定:
- `2`

意味:
- 1 beam あたり最低点数

現行挙動:
- 実データ点数がこれより少ない場合は、コードが自動的に可能な範囲まで緩和する

### `--min-scans-per-beam`
既定:
- `2`

意味:
- 1 beam あたり最低 scan 数

現行挙動:
- 実際の最大 scan 数に合わせて自動緩和されることがある

## 6.4 fit の current behavior

### primary stream 解決
通常は config から primary stream を解決します。  
明示 `--fit-stream-name` があればそれを優先します。

### single-beam compatibility mode
fit 対象に 2 つ未満の beam ID しか残らない場合は、本格幾何 fit をせず、compatibility mode になります。  
このとき

- `selected_rows_singlebeam_compat.csv`
- `fit_summary.txt`

は出ますが、beam geometry fit 自体は行いません。

### `recommended_model`
`fit_summary.txt` には residual に基づく `recommended_model` が記録されます。

---

## 7. `sunscan_make_pseudo_multibeam`

## 7.1 目的

- single-beam summary CSV から疑似 multi-beam summary を作る

## 7.2 基本例

```bash
sunscan_make_pseudo_multibeam \
  out_single/sun_scan_summary_tag.csv \
  --spectrometer-config beams.toml \
  --outdir out_pseudo
```

## 7.3 引数

### 位置引数 `singlebeam_summary_csv`
型:
- 1 個以上

意味:
- single-beam summary CSV

### `--spectrometer-config`
- 必須
- 擬似 multi-beam の target stream 配置を決める

### `--outdir`
- 出力先

### `--stream-name`
- 特定 stream だけ擬似生成したいときに使う

### `--noise-arcsec`
既定:
- `0.0`

意味:
- 擬似 jitter [arcsec]

### `--seed`
既定:
- `0`

意味:
- random seed

### `--tag`
意味:
- 出力ファイル名の base tag 上書き

### `--rep-el-deg`
型:
- repeatable float

意味:
- source summary の representative EL を人工的に差し替えて、複数 EL 条件の擬似データを増やす

重要:
- これは **standalone の `sunscan_make_pseudo_multibeam`** では使える
- しかし current `sunscan_multibeam pseudo` ラッパはこの値を実行側へ渡していない
- したがって `--rep-el-deg` を使いたいときは **standalone コマンドを直接使う**方が安全

---

## 8. `sunscan_multibeam`

## 8.1 目的

- extract / fit / pseudo / check-config を一つにまとめた wrapper

## 8.2 subcommand 一覧

### `sunscan_multibeam extract`
- `sunscan_extract_multibeam` の wrapper

### `sunscan_multibeam fit`
- `sunscan_fit_multibeam` の wrapper

### `sunscan_multibeam pseudo`
- `sunscan_make_pseudo_multibeam` の wrapper
- current 実装では `--rep-el-deg` の扱いに注意

### `sunscan_multibeam check-config`
- validator の wrapper

## 8.3 wrapper を使うべき場面

- 通常の運用でコマンド体系を統一したいとき

## 8.4 standalone を使うべき場面

- 個別コマンドの挙動を厳密に追いたいとき
- pseudo の `--rep-el-deg` のような細かい制御が必要なとき
- デバッグ時

---

## 9. config 側のパラメータの current interpretation

## 9.1 `beam_fit_use`
同じ beam 内で代表 stream を選ぶための印です。

## 9.2 `use_for_fit`
fit 候補集合への参加可否です。

## 9.3 `use_for_sunscan`
extract / singlebeam の自動選択に効きます。

## 9.4 `use_for_convert`
multi-beam package 内では主に metadata として保持されます。

## 9.5 `enabled`
総合無効化です。

---

## 10. よくある使い分け

### XX を fit 代表にし、YY は extract だけ残す
- XX: `use_for_fit=true`, `beam_fit_use=true`
- YY: `use_for_fit=false`, `beam_fit_use=false`
  または
- YY: `use_for_fit=true`, `beam_fit_use=false`
  でもよいが、意図を明確にしたいなら前者が分かりやすい

### ある分光計を一時的に解析から外す
```toml
enabled = true
use_for_convert = true
use_for_sunscan = false
use_for_fit = false
beam_fit_use = false
```

### 完全に無効化する
```toml
enabled = false
```

---

## 11. 解析結果に残る主要 metadata

single / extract の summary には、少なくとも次が残ります。

- `spec_time_basis`
- `spec_time_suffix`
- `spec_time_fallback_field`
- `spec_time_example`
- `azel_source`
- `altaz_apply`
- `spectrometer_time_offset_sec`
- `encoder_shift_sec`
- `encoder_az_time_offset_sec`
- `encoder_el_time_offset_sec`
- `encoder_vavg_sec`
- `chopper_wheel`
- `ripple_remove`
- `ripple_preset`
- `edge_fit_win_deg`
- `edge_fit_threshold`
- `hpbw_init_arcsec`
- `trim_scan`
- `profile_xlim_deg`
- `db_namespace`
- `telescope`
- `tel_loaddata`
- `planet`

これにより、summary CSV 単体からでも「どういう条件で解析したか」をかなり追跡できます。

---

## 12. 実務上の注意

1. 明示的に stream を指定しない限り、usage flags が効く  
2. `encoder_shift_sec` と Az/El 個別 offset を同時に入れるときは意味を整理しておく  
3. `beam_fit_use` と `use_for_fit` を混同しない  
4. single で成功してから extract に進む  
5. pseudo の高度な制御は standalone を優先する

---

## 13. この分冊のまとめ

本分冊では、multi-beam package 側の CLI 全体を説明しました。  
次の `v4-03 cookbook` では、実際にコピペしやすい TOML とコマンド例をまとめます。
