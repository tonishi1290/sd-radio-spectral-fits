# converter 詳細マニュアル

対象スクリプト:

- `necst_v4_sdfits_converter.py`

---

## 1. converter の目的

converter は NECST RawData を読み、1 本の SDFITS ファイルへ変換するツールです。

このスクリプトが扱うものは大きく次の 5 種類です。

1. spectrometer table
2. encoder table
3. altaz table
4. WCS table
5. weather table

そして最終的に、各 spectrometer stream について

- spectrum
- mode
- pointing
- weather
- beam model
- WCS / frequency axis

を統合し、writer に渡して FITS 行を生成します。

---

## 2. 全体フロー

現行実装の概念的な流れは次です。

1. CLI を読む
2. spectrometer config TOML を読む（または legacy single-stream config を内部生成する）
3. `global` 設定から table 名・namespace・telescope 名などを解決する
4. encoder / altaz / WCS / weather を開く
5. stream ごとに spectrometer table を読む
6. spectrometer 時刻列を選ぶ
7. encoder / altaz を spectrometer 時刻へ内挿する
8. boresight を求める
9. beam offset / rotation を適用して beam-center Az/El を求める
10. `--radec-method` と `--radec-azel-source` に応じて RA/DEC を求める
11. meteorology / refraction を処理する
12. mode ごとに FITS 行を構築する
13. 最終的に 1 本の FITS として保存する

---

## 3. converter が使う RawData / DB について

現行実装では、spectral / encoder / altaz / weather の主経路は **necstdb 直読みに寄っています**。

したがって converter は、RawData ディレクトリ中にある NECST DB table 名に強く依存します。

特に重要なのは、spectrometer config 中の以下です。

- `global.db_namespace`
- `global.telescope`
- `spectrometers.db_stream_name`
- `spectrometers.db_table_name`

`db_table_name` が明示されていれば、それが最優先です。
`db_table_name` が無ければ `db_stream_name` と `db_namespace`, `telescope` から table 名が組み立てられます。

そのため、**実際の RawData の table 名と config が一致していないと、最初に失敗しやすい**のが converter です。

---

## 4. CLI パラメータ完全説明

### 4.1 位置引数

#### `rawdata`

RawData ディレクトリのパスです。

ここで指定するのは、**観測 run 全体の RawData ディレクトリ**です。
分光計名や stream 名を裸で後ろに追加してはいけません。

---

### 4.2 spectrometer config / naming / output

#### `--spectrometer-config`

converter-compatible な spectrometer TOML を指定します。

これを与えると、legacy 単一 stream モードではなく、TOML の `[[spectrometers]]` に従って stream 群を処理します。

#### `--strict-config`

現在の実装では **予約パラメータ**です。
警告は出ますが、厳密チェックが大幅に強化されるわけではありません。

#### `--spectral`

legacy 単一 stream 用の spectral 名です。
`spectrometer-config` を使わずに 1 stream だけ変換するときに使います。

#### `--telescope`

table 名構築に使う telescope 名です。
例: `OMU1P85M`

#### `--db-namespace`

table 名構築に使う namespace です。
例: `necst`

#### `--tel-loaddata`

現行 converter では **legacy compatibility 用で、実質未使用**です。
help に残っていますが、主経路の DB 読み出しはこれに依存しません。

#### `--out`

出力 FITS 名を明示します。
指定しないと既定のファイル名が作られます。

#### `--object`

FITS の `OBJECT` に入れる名前です。
省略時は `rawdata` の basename が使われます。

#### `--project`

FITS の project ID 用文字列です。

#### `--observer`

FITS の observer 文字列です。

---

### 4.3 site / telescope constants

#### `--site-lat`

サイト緯度 [deg]。
Az/El から RA/DEC を計算するときに使います。

#### `--site-lon`

サイト経度 [deg, east+]。

#### `--site-elev`

サイト標高 [m]。

---

### 4.4 legacy frequency-axis parameters

これらは主に **legacy single-stream mode** 用の周波数軸パラメータです。
TOML で `frequency_axis` / `local_oscillators` を正しく与えるなら、通常はこちらより TOML を優先して運用します。

#### `--nchan`

legacy 単一 stream のチャンネル数。

#### `--if0-ghz`

legacy IF 開始周波数 [GHz]。

#### `--if1-ghz`

legacy IF 終了周波数 [GHz]。

#### `--lo1-ghz`

legacy LO1 [GHz]。

#### `--lo2-ghz`

legacy LO2 [GHz]。

#### `--restfreq-hz`

legacy RESTFREQ [Hz]。

---

### 4.5 channel selection

#### `--channel-slice`

変換時に channel 範囲を切り出します。

書式例:

- `[1024,8192)`
- `[1024,8191]`

CLI で指定した場合は、TOML の global / stream ごとの channel slice より強い扱いになります。

---

### 4.6 encoder / altaz / interpolation

#### `--encoder-time-col`

encoder table の時刻列名です。
推奨は `time` です。

`recorded_time` は現行実装では未対応で、指定しても `time` に強制されます。

#### `--encoder-shift-sec`

encoder timestamp に加える一様シフト [s] です。
内部的には

`t_enc <- t_enc + encoder_shift_sec`

としてから spectrometer 時刻へ内挿します。

#### `--altaz-time-col`

altaz table の時刻列名です。
encoder と同様に `recorded_time` は強制的に `time` へ落とされます。

#### `--interp-extrap`

内挿の外挿ポリシーです。

- `hold`: 端では値を保持
- `nan`: 範囲外は NaN

beam-center Az/El や RA/DEC を最後まで通したいなら、通常は `hold` の方が安全です。

---

### 4.7 mode / row selection

#### `--use-modes`

変換対象とする `OBSMODE` の一覧です。
既定は `ON,HOT,OFF` です。

#### `--include-unknown`

`UNKNOWN` サンプルも含めます。

#### `--scan-key`

scan grouping に使うキーです。

- `id`
- `id_mode`

#### `--skip-nonstandard-position`

非標準 position を落とすためのオプションです。
特殊な運用をしていない限り、省略で問題ありません。

---

### 4.8 pointing-error diagnostics

#### `--pe-rms-warn-arcsec`

pointing error RMS の warning 閾値 [arcsec]。

#### `--pe-max-warn-arcsec`

pointing error 最大値の warning 閾値 [arcsec]。

#### `--collapse`

writer 側の出力行 collapsing に関係する互換オプションです。
運用では既定のままでよいことが多いです。

#### `--pe-cor`

legacy compatibility 用で、現行実装では未使用です。

#### `--plot-pe`

現行 multi-stream 実装では未サポート寄りです。reserved と考えてください。

#### `--print-pe`

pointing error 情報を標準出力へ出します。

#### `--plot-pe-outdir`

PE plot の出力先ディレクトリ。

#### `--plot-pe-max-scans`

PE plot の最大 scan 数。

---

### 4.9 WCS / RADEC / refraction

#### `--wcs-table`

WCS table 名です。
既定は `necst-telescope-coordinate-wcs` です。

#### `--radec-method`

RA/DEC の求め方です。

- `wcs_first`: まず WCS table を使い、足りないときだけ Az/El 由来に fallback
- `azel`: WCS を使わず、Az/El 由来で計算

#### `--radec-azel-source`

Az/El 由来の RA/DEC を作るとき、どの Az/El を使うかを選びます。

- `beam`
- `true`
- `encoder`
- `altaz`
- `cmd`

詳細は後述します。

#### `--refraction`

Az/El→RA/DEC 計算時に refraction を使うかどうかです。

- `on`
- `off`

#### `--met-pressure-hpa`

refraction 用の pressure override [hPa]。

#### `--met-temperature-c`

refraction 用の temperature override [degC]。

#### `--met-humidity-pct`

refraction 用の humidity override [%]。

#### `--met-obswl-um`

refraction 計算に渡す観測波長 [um]。

#### `--met-source`

meteorology の採り方です。

- `auto`
- `weather`
- `spectral`
- `fallback`

#### `--weather-table`

weather table 名です。
`auto` なら `global.weather_table` または既定 table 名を使います。

#### `--weather-time-col`

weather table の時刻列名です。

---

## 5. spectrometer config TOML のうち converter が使うもの

### 5.1 `[global]`

現行 converter で実際に意味を持つ主なキーは次です。

- `db_namespace`
- `telescope`
- `encoder_table`
- `altaz_table`
- `weather_table`
- `encoder_shift_sec`
- `output_layout`
- `time_sort`
- `channel_slice`

このうち重要なのは次です。

#### `global.db_namespace`

table 名構築に使います。

#### `global.telescope`

table 名構築に使います。

#### `global.encoder_table`

encoder table 名を直接固定します。
既定の自動命名と違う場合に有効です。

#### `global.altaz_table`

altaz table 名を直接固定します。

#### `global.weather_table`

weather table 名を直接固定します。

#### `global.encoder_shift_sec`

converter で実際に使われる global 側の時刻補正です。
現行 converter の global 側で直接サポートされている time offset は、ここでは **encoder_shift_sec** が中心です。

#### `global.output_layout`

現行実装では `merged` / `merged_time` だけが許容されます。

#### `global.time_sort`

`merged` / `merged_time` では `true` が必要です。

#### `global.channel_slice`

global な channel 切り出しです。
stream ごとの channel slice や CLI `--channel-slice` と組み合わせて使います。

---

### 5.2 `[[spectrometers]]`

converter で重要なキーです。

- `name`
- `fdnum`
- `ifnum`
- `plnum`
- `polariza`
- `beam_id`
- `frontend`
- `backend`
- `sampler`
- `db_stream_name`
- `db_table_name`
- `channel_slice`

#### `name`

内部 stream 名です。

#### `fdnum`, `ifnum`, `plnum`, `polariza`

出力行や stream の識別に使います。

#### `beam_id`

beam model と関連づく beam 識別子です。

#### `db_stream_name`

自動 table 名構築に使う stream 名です。

#### `db_table_name`

spectral table 名を完全指定します。
実データの table 名が規則と合わない場合に最も安全です。

#### `frontend`, `backend`, `sampler`

メタデータとして保持されます。

---

### 5.3 `[spectrometers.frequency_axis]`

converter が WCS を作るための主表です。

主なキー:

- `nchan`
- `definition_mode`
- `band_start_hz`
- `band_stop_hz`
- `channel_origin`
- `reverse`
- `ctype1`
- `cunit1`
- `specsys`
- `veldef`
- `store_freq_column`
- `restfreq_hz`

特に重要なのは `restfreq_hz` です。

---

### 5.4 `[spectrometers.local_oscillators]`

LO / sideband 情報を保持します。

主なキー:

- `lo1_hz`
- `lo2_hz`
- `lo3_hz`
- `sb1`
- `sb2`
- `sb3`
- `obsfreq_hz`
- `imagfreq_hz`
- `sideband`

---

### 5.5 `[spectrometers.beam]`

beam model の中核です。

- `az_offset_arcsec`
- `el_offset_arcsec`
- `rotation_mode`
- `reference_angle_deg`
- `rotation_sign`
- `dewar_angle_deg`
- `beam_model_version`

converter はこれを使って boresight から beam-center を作ります。

---

## 6. boresight / beam-center / RADEC の意味

converter を正しく使う上で、ここが最も重要です。

### 6.1 encoder / altaz / correction

現行実装では、pointing 正規化の中で大まかに以下を作っています。

- encoder Az/El を spectrometer 時刻へ内挿した系列
- altaz/cmd 系列を spectrometer 時刻へ内挿した系列
- correction `dlon / dlat`

### 6.2 boresight

現在の boresight は、概念的には

`boresight = encoder - correction`

です。

### 6.3 beam-center

そこに stream 固有の beam offset / rotation を適用して

`beam-center = boresight + beam offset`

を作ります。

### 6.4 `--radec-azel-source`

#### `beam`

現在の標準です。

- boresight = `encoder - correction`
- そこに beam offset / rotation を適用
- その結果を RA/DEC fallback に使う

#### `true`

beam offset をかける前の boresight を使います。

- `true = encoder - correction`

中心ビームでは `beam` とほぼ同じになりますが、周辺ビームでは異なります。

#### `encoder`

encoder 系そのものに beam offset を適用したものを使うモードです。
補正量を引かないため、`beam` や `true` とは意味が違います。

#### `altaz`

raw の altaz/cmd 系列に beam offset を適用したものを使います。

#### `cmd`

ごく重要なモードです。
現行実装では、`cmd` は

- `altaz/cmd - correction` を boresight とし
- その後に beam offset / rotation を適用した Az/El

を Az/El 由来の RA/DEC に使います。

position-switch 的に「指令値側を基準に RA/DEC を作りたい」場合は、この `cmd` が最も近い意味になります。

---

## 7. 実運用上の要点

### 7.1 table 名が合わないとき

次のようなエラーは非常によく起きます。

- `Table 'necst-...-data-spectral-xffts-board0' does not exist.`

この場合は、config の `db_stream_name` が RawData 内の実 table 名と合っていません。

最も安全な対処は `db_table_name` を使うことです。

### 7.2 WCS より Az/El を優先したいとき

`--radec-method azel` を使います。
さらに source を切り替えたいなら `--radec-azel-source` を指定します。

### 7.3 position switch 的に `cmd` を使いたいとき

```bash
python necst_v4_sdfits_converter.py \
  --spectrometer-config 12co.conf \
  --radec-method azel \
  --radec-azel-source cmd \
  rawdata_dir
```

とします。

### 7.4 global と CLI の関係

現在の converter では、naming と table 解決に関しては `global` と CLI の両方が関与します。

- CLI を明示したもの
- `global` に書かれたもの
- 既定値

の順で解決される項目があります。

---

## 8. converter 特有の troubleshooting

### 8.1 `unrecognized arguments`

この converter の位置引数は `rawdata` 1 個です。
stream 名は裸で追加せず、必要なら `--spectral` を使ってください。

### 8.2 `Table ... does not exist`

- `db_stream_name` が違う
- `db_table_name` を書くべき
- `global.telescope` や `global.db_namespace` が違う

のいずれかが主因です。

### 8.3 `recorded_time is not supported`

現行版では `encoder_time_col`, `altaz_time_col`, `weather_time_col` に `recorded_time` は使えません。`time` に強制されます。

### 8.4 `Selected RA/DEC Az/El source ... contains NaN`

選んだ source の Az/El が NaN です。
`interp-extrap=hold`、table coverage、beam offset、極端な elevation を確認してください。

### 8.5 `Beam-center Az/El contains NaN`

beam offset / rotation 適用後に beam-center が壊れています。beam 定義または入力 elevation を確認してください。
