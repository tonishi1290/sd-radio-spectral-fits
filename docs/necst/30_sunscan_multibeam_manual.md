# sunscan multibeam 詳細マニュアル

対象スクリプト:

- `multibeam_beam_measurement/sunscan_extract_multibeam.py`
- `multibeam_beam_measurement/sunscan_fit_multibeam.py`
- `multibeam_beam_measurement/sunscan_multibeam.py`
- `multibeam_beam_measurement/synthetic_multibeam.py`
- `multibeam_beam_measurement/check_spectrometer_config.py`

---

## 1. multibeam の全体像

multibeam は 1 本のコマンドではなく、次の機能群からなります。

1. **extract**
   - 各 stream について singlebeam 解析を繰り返し
   - それらを 1 本の aggregated summary CSV にまとめる
2. **fit**
   - aggregated summary CSV から beam geometry を推定する
3. **pseudo**
   - singlebeam summary から擬似 multi-beam データを作る
4. **check-config**
   - spectrometer TOML の妥当性を検査する
5. **wrapper (`sunscan_multibeam.py`)**
   - `extract / fit / pseudo / check-config` を 1 つの入口にまとめる

---

## 2. extract の目的

extract は multi-beam 解析の前段です。

この段階ではまだ beam geometry fitting はしません。
やることは次です。

- spectrometer config から対象 stream 群を決める
- 各 stream ごとに singlebeam 解析を走らせる
- per-stream の `sun_scan_summary_*.csv` を保存する
- それらを連結して `sunscan_multibeam_scan_summary_*.csv` を作る
- どの stream を解析し、どの stream を skip したかを manifest に残す
- config snapshot と stream table を保存する

つまり extract は、**multi-beam fit に使える形へデータを正規化・集約する段階**です。

---

## 3. extract の CLI パラメータ完全説明

extract の CLI は singlebeam にかなり近いですが、stream 群を扱うための追加要素があります。

### 3.1 入出力 / stream selection

#### `rawdata`

RawData ディレクトリ。

#### `--spectrometer-config`

必須です。stream 群の定義を含む TOML を渡します。

#### `--outdir`

出力ディレクトリ。

#### `--run-id`

run ID を明示上書きします。
省略時は rawdata ディレクトリ名から作られます。

#### `--stream-name`

特定 stream 名だけを選ぶ repeatable option です。
明示した場合、通常の `use_for_sunscan` フィルタよりこちらが優先されます。

### 3.2 `global` から共有される入力系設定

extract は config が与えられると、以下を `global` から読めます。

- `db_namespace`
- `telescope`
- `tel_loaddata`
- `planet`
- `spectrometer_time_offset_sec`
- `encoder_shift_sec`
- `encoder_az_time_offset_sec`
- `encoder_el_time_offset_sec`

CLI 明示が優先です。

### 3.3 解析パラメータ

singlebeam と同じです。

- `--azel-source`
- `--altaz-apply`
- `--spectrometer-time-offset-sec`
- `--encoder-shift-sec`
- `--encoder-az-time-offset-sec`
- `--encoder-el-time-offset-sec`
- `--encoder-vavg-sec`
- `--no-chopper-wheel`
- `--tamb-k`
- `--chopper-win-sec`
- `--chopper-stat`
- `--profile-xlim-deg`
- ripple 関連一式
- edge fit 関連一式
- trim 関連一式
- `--continue-on-error`
- `--debug-plot`
- `--dish-diameter-m`
- `--hpbw-factor`

意味は singlebeam マニュアルと同じです。

---

## 4. stream usage flags

multi-beam では、各 `[[spectrometers]]` について使用方針を分けられます。

- `enabled`
- `use_for_convert`
- `use_for_sunscan`
- `use_for_fit`
- `beam_fit_use`

### 4.1 `enabled`

総合スイッチです。
`false` なら他の `use_for_*` が `true` でも最終的に無効寄りになります。

### 4.2 `use_for_sunscan`

extract の通常実行で、その stream を解析対象にするかを決めます。

### 4.3 `use_for_fit`

fit の既定 stream 選択に使います。

### 4.4 `beam_fit_use`

同じ `beam_id` に複数 fit-enabled stream がある場合、**どれを primary stream にするか**を示します。

典型例:

- `B01` に `east115_xx`, `east115_yy` がある
- fit に使うのは XX だけにしたい

その場合

- `east115_xx: use_for_fit=true, beam_fit_use=true`
- `east115_yy: use_for_fit=false, beam_fit_use=false`

とするのが明快です。

### 4.5 explicit `--stream-name` / `--fit-stream-name`

CLI 明示指定は最優先です。
そのため、一時的に `use_for_sunscan=false` の stream を試験的に選ぶこともできます。

---

## 5. extract の出力物

### 5.1 per-stream 出力

各 stream ごとに singlebeam 出力相当のファイルがサブディレクトリへ書かれます。

### 5.2 aggregated summary CSV

#### `sunscan_multibeam_scan_summary_<tag>.csv`

extract が最終的に fit 用へ渡す集約 CSV です。
各行は基本的に 1 scan × 1 stream です。

主な列:

- `run_id`
- `tag`
- `stream_name`
- `beam_id`
- `restfreq_hz`
- `polariza`
- `fdnum`, `ifnum`, `plnum`, `sampler`
- `x_arcsec`, `y_arcsec`
- `scan_id`
- `rep_el_deg`
- singlebeam summary の主要列

### 5.3 manifest

#### `sunscan_multibeam_manifest_<tag>.csv`

何を実行し、何を skip / error にしたかを記録する管理表です。

特に重要な status:

- `ok`
- `skipped`
- `error`

skip された stream も manifest に残るため、「設定にはあるが今回使わなかった stream」を後から確認できます。

### 5.4 stream table

#### `spectrometer_stream_table_<tag>.csv`

config から作られた stream 一覧です。

- beam_id
- restfreq_hz
- fdnum / ifnum / plnum / polariza
- `enabled / use_for_convert / use_for_sunscan / use_for_fit / beam_fit_use`
- nominal beam offsets / rotation

を確認できます。

### 5.5 config snapshot

#### `analysis_config_snapshot_<tag>.json`

extract 実行時に使った実際の解析設定を JSON へ保存します。

これは再現性のために非常に重要です。

---

## 6. fit の目的

fit は aggregated summary CSV から、複数 beam の相対配置と回転モデルを推定します。

fit が最終的に出したいものは、概ね次です。

- 各 beam の nominal offset `(az_offset_arcsec, el_offset_arcsec)`
- 回転モデルの基準角 `reference_angle_deg`
- 回転符号 `rotation_sign`
- 回転 slope `rotation_slope_deg_per_deg`
- dewar angle `dewar_angle_deg`

これらを beam model TOML として書き戻すことで、将来の converter / 解析で使える形にします。

---

## 7. fit の入力データが表すもの

fit が期待する aggregated summary CSV には最低限次が必要です。

- `beam_id`
- `stream_name`
- `x_arcsec`
- `y_arcsec`
- `rep_el_deg`
- `run_id`
- `scan_id`

また、`fit_ok_az`, `fit_ok_el` がある場合は、false な行は事前に落とされます。

つまり fit の入力は、「各 scan で各 beam がどこに見えたか」を並べた点群です。

---

## 8. primary stream 解決

fit は、spectrometer config と aggregated summary から **primary stream** を決めます。

### 8.1 既定動作

1. まず `use_for_fit=true` の stream だけを残す
2. 同じ `beam_id` ごとに group を作る
3. `beam_fit_use=true` が 1 本だけあればそれを primary にする
4. group に 1 本しかなければそれを使う
5. group に複数あるのに primary が無ければエラー

### 8.2 明示指定

`--fit-stream-name` を使えば、その指定がそのまま primary 扱いになります。

これは XX / YY を切り替えたいときに便利です。

---

## 9. fit モデルの意味

現在の fit モデル選択肢は次です。

- `center_beam`
- `virtual_center`
- `both`

### 9.1 `center_beam`

各 scan について、指定した `center_beam_id` の位置をその scan の基準 `(x0, y0)` とみなします。

つまり相対座標は

`x_rel = x_beam - x_center_beam`

`y_rel = y_beam - y_center_beam`

です。

このモデルは、中心 beam の位置が信頼できるときに自然です。

### 9.2 `virtual_center`

各 scan について、全 beam の平均位置をその scan の基準 `(x0, y0)` とみなします。

つまり

`x0 = mean(x_beams)`

`y0 = mean(y_beams)`

です。

中心 beam を固定したくない場合、または対称配列で幾何中心を使いたい場合に有効です。

### 9.3 `both`

可能なら両方を試し、平均 residual が小さい方を推奨モデルとして出します。

ただし `center_beam` を試すには `--center-beam-id` が必要です。

---

## 10. fit が実際に何をしているか

ここが最も重要です。

### 10.1 relative points の構築

まず各 scan について relative point cloud を作ります。

- `center_beam`: 中心 beam を原点化
- `virtual_center`: 全 beam の平均を原点化

### 10.2 template pattern の推定

各 beam の relative position を complex 数

`z = x + i y`

として扱い、scan ごとの回転を戻しながら平均して template pattern を作ります。

### 10.3 scan ごとの回転角推定

各 scan の beam 配置が template に対して何度回転しているかを推定します。

### 10.4 `theta(rep_el)` の線形近似

各 scan に representative elevation `rep_el_deg` があるので、

`theta_deg = slope * (rep_el_deg - reference_angle_deg) + intercept`

という線形モデルを当てます。

### 10.5 `rotation_sign` と `dewar_angle_deg`

slope の符号を `rotation_sign = ±1` に量子化し、残りを `dewar_angle_deg` としてまとめます。

### 10.6 residual の計算

template pattern を scan ごとの回転角で再回転し、観測点との差

- `dx`
- `dy`
- `resid_arcsec = sqrt(dx^2 + dy^2)`

を計算します。

### 10.7 sigma clipping

`resid_arcsec` に対し MAD ベースの sigma clipping を行い、外れ値を落として再推定します。

### 10.8 minimum-count filtering

各 beam ごとに

- `min_points_per_beam`
- `min_scans_per_beam`

を満たさないものは除外されます。

入力が少なすぎる場合は threshold が自動的に緩和されることがあります。

---

## 11. single-beam compatibility mode

fit に渡した summary CSV に beam_id が 1 種類しか無い場合、幾何 fit はできません。
その場合は single-beam compatibility mode になり、

- `selected_rows_singlebeam_compat.csv`
- `fit_summary.txt`

だけを出して終了します。

これは異常ではなく、singlebeam データで package 全体の流れを確認するときの想定動作です。

---

## 12. fit の CLI パラメータ完全説明

### 12.1 入力 / 選択

#### `summary_csv`

1 本以上の aggregated summary CSV。
複数 run をまとめて fit できます。

#### `--spectrometer-config`

必須。stream / beam 定義と primary stream 解決に使います。

#### `--outdir`

出力ディレクトリ。

#### `--center-beam-id`

`center_beam` model で中心とみなす beam ID。

#### `--fit-stream-name`

fit に使う stream を明示選択します。repeatable です。

### 12.2 model / geometry

#### `--model`

- `both`
- `center_beam`
- `virtual_center`

#### `--reference-angle-deg`

reference elevation を固定上書きします。
省略時は `rep_el_deg` の median を使います。

### 12.3 outlier rejection / minimum-count

#### `--sigma-clip`

residual clipping の閾値。
`<= 0` なら無効です。

#### `--clip-iters`

sigma clipping の反復回数。

#### `--min-points-per-beam`

各 beam に必要な最少点数。

#### `--min-scans-per-beam`

各 beam に必要な最少 scan 数。

---

## 13. fit の出力物

モデルごとに次が出ます。

### `beam_fit_results_<model>.csv`

各 beam の最終 offset と residual 統計。

主な列:

- `beam_id`
- `az_offset_arcsec`
- `el_offset_arcsec`
- `radius_arcsec`
- `phase_deg`
- `reference_angle_deg`
- `rotation_sign`
- `rotation_slope_deg_per_deg`
- `dewar_angle_deg`
- `rms_resid_arcsec`
- `max_resid_arcsec`

### `run_shift_results_<model>.csv`

各 scan の shift と推定回転角。

### `beam_fit_residuals_<model>.csv`

残した点の residual table。

### `beam_fit_rejected_<model>.csv`

sigma clipping で落ちた点。

### `beam_model_<model>.toml`

元の spectrometer config TOML をベースに、beam model を更新した TOML。
将来の converter / sunscan で使うための重要成果物です。

### `fit_summary.txt`

どのモデルが推奨か、平均 residual はどれかをまとめたテキストです。

---

## 14. pseudo multibeam

`synthetic_multibeam.py` は、singlebeam summary を複製・平行移動して疑似 multi-beam summary を作るツールです。

### 14.1 目的

- 実 multi-beam 生データが無いときに fit の流れを確認する
- beam offsets / rotation_mode の sanity check を行う
- config に入れた幾何が dry-run で破綻しないか確認する

### 14.2 仕組み

各 singlebeam summary の `(x_arcsec, y_arcsec)` に対し、config に書かれた nominal beam offset を適用し、必要なら乱数 jitter を加えて複数 stream の summary を人工的に生成します。

### 14.3 パラメータ

#### `singlebeam_summary_csv`

元となる singlebeam summary CSV。

#### `--spectrometer-config`

必須。pseudo で使う beam 配置を与えます。

#### `--outdir`

出力先。

#### `--stream-name`

pseudo 生成する stream を限定します。

#### `--noise-arcsec`

人工 jitter の標準偏差 [arcsec]。

#### `--seed`

乱数 seed。

#### `--tag`

出力タグを上書き。

#### `--rep-el-deg`

representative elevation を override して複数 scan を人工生成したいときに使います。repeatable です。

---

## 15. check-config

`check_spectrometer_config.py` は、config の静的検査ツールです。

### 15.1 やること

- stream table を表示する
- primary stream を解決する
- duplicate beam ID を確認する
- restfreq, rotation_sign, zero-offset dry-run 上の注意を warning として出す

### 15.2 パラメータ

#### `spectrometer_config`

TOML パス。

#### `--stream-name`

primary stream を明示的に与えたいときに使います。

#### `--out-csv`

stream table を CSV 保存します。

---

## 16. wrapper `sunscan_multibeam.py`

これは subcommand wrapper です。

- `extract`
- `fit`
- `pseudo`
- `check-config`

を 1 つの入口にまとめています。

中身は薄く、各サブコマンドの実体はそれぞれの専用スクリプトです。

---

## 17. multibeam 特有の注意

### 17.1 同じ beam に複数 stream がある場合

XX / YY の両方を `use_for_fit=true` にしたまま `beam_fit_use` を指定しないと、primary stream 解決でエラーになります。

### 17.2 nominal beam offsets が全部 0 の場合

pseudo dry-run では回転を全く拘束できません。`check-config` でも warning が出ます。

### 17.3 rotation_mode が全部 `none` の場合

elevation 依存の回転を検証できません。

### 17.4 `center_beam` と `virtual_center`

どちらが良いかは観測配置に依存します。

- 中心 beam を明確に信じるなら `center_beam`
- 全体の平均中心で安定化したいなら `virtual_center`

です。

`both` を使って比較するのが安全です。
