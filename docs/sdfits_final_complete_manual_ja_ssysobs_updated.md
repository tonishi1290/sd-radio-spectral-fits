# sd_radio_spectral_fits SDFITS 最終完全版説明書
## 0. この文書について
本書は、`sd_radio_spectral_fits` が前提とする SDFITS の仕様を、**writer / 外部 producer / reader / standardizer / coadd / baseline / calibration** の全体像が一冊で分かるように統合した最終版です。

この改訂版は、今回添付された最新版解析パッケージのコア実装に合わせて更新しています。現アーカイブのコア版は概ね `1.2.1` 系で、主対象は次のモジュールです。

- `sdfits_writer.py`
- `sdfits_bintable.py`
- `fitsio.py`
- `rawspec.py`
- `scantable_utils.py`
- `restfreq.py`
- `regrid_vlsrk.py`
- `coadd.py`
- `calibrate.py`
- `baseline.py`
- `tempscale.py`
- `doppler.py`
- `axis.py`
- `atmosphere.py`

本書では、次の3本の文書を統合しています。

1. writer 本体仕様
2. 列一覧の完全表形式版
3. reader / standardizer / coadd / baseline / calibration 側の前提契約版

重複する内容は完全には削らず、**読み手がどの部から読んでも必要情報に到達できる**ことを優先しています。ただし、構造は次のように整理しました。

- **第I部**: writer と外部 producer が満たすべき SDFITS 本体仕様
- **第II部**: 列・キーワード・履歴の完全表
- **第III部**: reader / standardizer / coadd / baseline / calibration が期待する契約

ここでいう「外部 producer」とは、converter、観測系 writer、あるいは別プロセスが生成する SDFITS 出力を指します。今回のアーカイブには特定 converter 本体は含まれていませんが、**外部 producer がこのコア実装へ渡すべき契約**は引き続き重要なので、本書では削らずに明示します。

特に重要な前提は以下です。

- `DATA` は **fixed-length でも VLA でもよい**
- `FLAG` と `FREQ` は `DATA` と row ごとに整合していなければならない
- `POLARIZA` は物理偏波ラベルとして必須
- `FDNUM` / `IFNUM` / `PLNUM` は row 識別用の中核列
- 時刻解決は `TIMESTAMP -> MJD -> DATE-OBS + TIME -> 旧 TIME` を採用する
- `RESTFREQ` / `VELDEF` は常設ではなく、文脈依存で必須
- 読込時の温度スケールは **非破壊** で、`TR*` を勝手に `TA*` へ戻さない
- `VELOSYS` / `VFRAME` の公開意味論は **row ごとの未適用速度補正 [m/s]** である

---

## 第I部: writer が前提とする SDFITS 本体仕様

## sd_radio_spectral_fits が前提とする SDFITS 仕様書

### 1. 目的

本書は、`sd_radio_spectral_fits` パッケージのうち、主に次のファイル群に基づいて、このパッケージが前提とする SDFITS の構造と、`sdfits_writer.py`、`fitsio.py`、`regrid_vlsrk.py`、`coadd.py`、`baseline.py` が想定している列・ヘッダ・履歴・時刻・偏波・beam・WCS の仕様を、実装ベースで詳細に説明するものです。

- `sdfits_writer.py`
- `sdfits_bintable.py`
- `fitsio.py`
- `restfreq.py`
- `regrid_vlsrk.py`
- `coadd.py`
- `calibrate.py`
- `baseline.py`
- `tempscale.py`
- `doppler.py`
- `axis.py`

ここでいう SDFITS は、古典的な標準をそのまま再現するというより、**単一鏡スペクトル解析に必要な provenance と row metadata を十分に保持するための SDFITS-like 実装**です。したがって、本書では「一般論としての SDFITS」ではなく、**このパッケージが実際に書く FITS の仕様**を定義します。

また、外部 converter / producer が存在する場合も、最終的にこのコア実装へ引き渡される以上、**producer 側が守るべき handoff 契約**も明示します。

---

### 2. 基本方針

#### 2.1 本実装は `SINGLE DISH` BinTable を核とする

出力 FITS は基本的に次の3層で構成されます。

1. `PRIMARY` HDU
2. `SINGLE DISH` BinTableHDU
3. 必要に応じて `HISTORY` BinTableHDU

`PRIMARY` にはファイル全体に共通する情報、`SINGLE DISH` には各スペクトル row の metadata とスペクトル本体、`HISTORY` には provenance を保存します。

#### 2.2 row は「1本のスペクトル + その row metadata」

1 row は単なる dump 番号ではなく、**その row を単独で解釈するのに必要な metadata を伴った1本のスペクトル**です。

したがって、同一時刻でも

- beam が違えば別 row
- IF が違えば別 row
- 偏波が違えば別 row
- sampler が違えば別 row

となります。

#### 2.3 `sdfits_bintable.py` が BinTable 生成の Single Source of Truth

`SINGLE DISH` テーブルの最終的な FITS Column 化は `sdfits_bintable.py` の `build_single_dish_table_hdu()` が担います。`sdfits_writer.py` は、

- どの列を出すか
- 各列にどんな意味を与えるか
- row ごとの値をどう蓄積するか
- どの header keyword を書くか

を担当し、最後の BinTable 生成を `sdfits_bintable.py` に委譲します。

---

### 3. この writer が書く FITS 全体構造

### 3.1 PRIMARY HDU

`PRIMARY` には主としてファイル共通情報が入ります。

代表的な keyword:

- `TELESCOP`
- `OBSERVER`
- `PROJID`
- `OBJECT`
- `DATE`
- `TIMESYS`
- `MJDSTART`, `MJDEND`
- `SWNAME`, `SWVER`, `ORIGIN`
- `SITELAT`, `SITELON`, `SITELONG`, `SITEELEV`
- `OBSGEO-X`, `OBSGEO-Y`, `OBSGEO-Z`
- 既定 WCS としての `CTYPE1`, `CUNIT1`, `CRVAL1`, `CDELT1`, `CRPIX1`, `SPECSYS`, `SSYSOBS`, `REFCHAN`
- 条件付きで `RESTFREQ`, `RESTFRQ`, `VELREF`, `VELDEF`

### 3.2 `SINGLE DISH` テーブル

`SINGLE DISH` には次が row ごとに入ります。

- スペクトル本体 `DATA`
- channel flag `FLAG`
- beam / IF / 偏波 / backend 情報
- 天球位置
- pointing / weather / source / scan 情報
- calibration 情報
- スペクトル軸 metadata (`CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`, `SPECSYS`, 必要に応じて `RESTFREQ`, `VELDEF`)

### 3.3 `HISTORY` テーブル

本実装では file-level history を FITS header card の `HISTORY` ではなく、**専用 BinTable (`EXTNAME='HISTORY'`)** として持ちます。

列は原則として

- `KEY`
- `VALUE`

の2列です。

この設計により、長い provenance や converter 設定を比較的安全に保存できます。

---

### 4. SDFITS-like であることの意味

本 writer は古典的標準にない convenience column も持ちます。たとえば

- `MJD`
- `DATEOBS`
- `TIMESTAMP`
- `FDNUM`, `IFNUM`, `PLNUM`, `POLARIZA`
- `FRONTEND`, `BACKEND`, `SAMPLER`
- `BORE_AZ`, `BORE_EL`, `BEAMXOFF`, `BEAMYOFF`, `BEAMROT`
- `AZ_CMD`, `EL_CMD`, `CORR_AZ`, `CORR_EL`, `CALC_REFR`

などです。

したがって、本パッケージでいう SDFITS は「SINGLE DISH テーブルを核にした実務向け単一鏡スペクトル FITS 実装」と理解してください。

---

### 5. 最重要事項: `DATA` は固定長でも可変長でもよい

これは本 writer の非常に重要な性質です。

#### 5.1 `DATA` は VLA を許す

`_ColumnBuffer.append_spectrum()` は、各 row の `data` を単に 1 次元 float32 配列として保持します。row ごとに長さが異なっても構いません。

その結果、この writer は次の両方に対応します。

- 全 row が同じチャンネル数の通常の fixed-length データ
- row ごとにチャンネル数が異なる可変長データ

#### 5.2 固定長なら fixed-length vector、異なれば VLA

`build_single_dish_table_hdu()` / `_spectrum_and_flag_columns()` は、row 長を見て自動的に判定します。

##### 固定長の場合

- `DATA` -> `nE`
- `FLAG` -> `nL`
- `FREQ` -> `nD`（有効時）

さらに `TDIMn=(nchan)` を書きます。

##### 可変長の場合

- `DATA` -> `PE`
- `FLAG` -> `PL`
- `FREQ` -> `PD`（有効時）

ここで

- `PE`: 可変長 float32 ベクトル
- `PD`: 可変長 float64 ベクトル
- `PL`: 可変長 logical ベクトル

です。

#### 5.3 `FLAG` は `DATA` と同長でなければならない

row ごとに `flag` を与える場合、長さは `data` と一致していなければなりません。与えない場合は自動で all-False が生成されます。

#### 5.4 `FREQ` も `DATA` と同長でなければならない

`store_freq_column=True` のとき `FREQ` ベクトル列を持てます。これも row ごとに `DATA` と同じ長さでなければなりません。

#### 5.5 `n_chan` は固定長強制値ではない

writer 初期化時の `n_chan` は、現在の VLA 対応実装では「既定または最大チャンネル数の目安」です。実際の row 長は `data` に従います。

#### 5.6 `NCHAN` / `NCHANSEL`

テーブル header に入る `NCHAN`, `NCHANSEL` は、fixed-length ならそのチャンネル数、VLA なら最大長のヒント値です。VLA では row ごとの実長は `TFORM='P*'` の配列を見て判断してください。

---

### 6. 列ポリシー

本 writer は列を4群に分けて扱います。

### 6.1 always

常に出す列です。

- `TIME`, `MJD`, `DATE-OBS`, `DATEOBS`, `TIMESTAMP`
- `SCAN`, `SUBSCAN`, `INTGRP`
- `OBJECT`, `OBSMODE`
- `EXPOSURE`, `CALSTAT`, `TEMPSCAL`, `FLAGROW`
- `RA`, `DEC`, `GLON`, `GLAT`
- `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`, `SPECSYS`
- `FDNUM`, `IFNUM`, `PLNUM`, `POLARIZA`
- `DATA`, `FLAG`

### 6.2 context-required

文脈依存で必須になる列です。

- `RESTFREQ`
- `VELDEF`
- heterodyne 文脈での `OBSFREQ`
- heterodyne 文脈での `SIDEBAND`

### 6.3 optional

意味のある値が1つでもあれば列全体を出します。

- `TCAL`, `THOT`, `TSYS`, `TAU0`
- `AZIMUTH`, `ELEVATIO`
- `VFRAME`, `FOFFSET`
- `BORE_AZ`, `BORE_EL`, `BEAMXOFF`, `BEAMYOFF`, `BEAMROT`
- `FRONTEND`, `BACKEND`, `SAMPLER`
- `IMAGFREQ`, `LO1FREQ`, `LO2FREQ`, `LO3FREQ`
- `SB1`, `SB2`, `SB3`
- `FREQ`

### 6.4 block-optional

関連列群のどれか1つが meaningful ならブロック全体を出します。

- source block: `SRCFRAME`, `SRCRDSYS`, `SRCEQNX`, `SRC_LONG`, `SRC_LAT`
- scan block: `SCANFRAM`, `SCANRDSYS`, `SCANEQNX`, `SCANX`, `SCANY`
- pointing block: `AZ_CMD`, `EL_CMD`, `CORR_AZ`, `CORR_EL`, `CALC_REFR`
- weather block: `TAMBIENT`, `PRESSURE`, `HUMIDITY`, `WINDSPD`, `WINDDIR`

### 6.5 未設定値の原則

optional / block-optional の列について、意味のない見せかけの既定値を入れるのは危険です。本実装では未設定をなるべく

- 数値列: `np.nan`
- 文字列列: `None`

で表し、meaningful な値があるときだけ列を出す方針です。

---

### 7. Time 仕様

これは重要なので、**採用仕様**と**現在の書き出し実装**を分けて説明します。

### 7.1 採用仕様

#### `DATE-OBS`

各 row の絶対 UTC 時刻を ISO 文字列で保存します。

#### `TIME`

`DATE-OBS` からの秒オフセットとして扱います。補助列です。

#### `MJD`

同じ時刻を MJD UTC days でも保存します。

#### `TIMESTAMP`

読込側で使う内部統一時刻です。読込側は原則として次の順で解決します。

1. `TIMESTAMP`
2. `MJD`
3. `DATE-OBS + TIME`
4. 旧 `TIME`

#### 旧互換

旧形式で `TIME=MJD days` を使っていたファイルも互換的に読めるようにします。

#### 実務上の推奨

絶対時刻は

- `TIMESTAMP`
- または `MJD`

を主として使い、`TIME` は補助列と考える、というのが本仕様です。

### 7.2 現在の writer 実装

現行 `sdfits_writer.py` は row 追加時に

- `MJD`: 実際の MJD UTC
- `DATEOBS`: ISO UTC 文字列
- `DATE-OBS`, `TIMESTAMP`: `DATEOBS` の複製

を保存しています。

一方で `TIME` は現実装では各 row に `0.0` を入れています。したがって、**書き出し側の現状では `TIME` はまだ秒オフセットを十分には表現していません**。

本書では、ユーザーが整理した仕様を正式な運用仕様として採用しつつ、現 writer 実装の現在地として

- `TIMESTAMP` / `MJD` / `DATE-OBS` は有効
- `TIME` は現状では将来仕様の受け皿であり補助列

と位置づけます。

---

### 8. beam / IF / polarization / spectrometer 情報の正式仕様

本仕様の中心は以下の列です。

- `FDNUM`
- `IFNUM`
- `PLNUM`
- `POLARIZA`
- `BACKEND`
- `SAMPLER`

必要に応じて `FRONTEND` も併用します。

### 8.1 `FDNUM`

- 意味: beam / feed 番号
- 型: int32 相当
- 推奨: 0 始まり
- 用途: beam 識別専用

### 8.2 `IFNUM`

- 意味: IF / spectral window 番号
- 型: int32 相当
- 推奨: 0 始まり
- 用途: IF 識別専用

### 8.3 `PLNUM`

- 意味: 装置内部の偏波チャネル番号
- 型: int32 相当
- 用途: 内部番号
- 重要: 物理偏波の意味は持たせない

### 8.4 `POLARIZA`

- 意味: 物理偏波ラベル
- 型: 文字列
- 必須: はい
- 許容値:
  - `XX`, `YY`
  - `RR`, `LL`
  - `XY`, `YX`
  - `RL`, `LR`
  - `I`, `Q`, `U`, `V`

`POLARIZA` は空文字不可、必ず確定した値を持たせます。

### 8.5 `BACKEND`

- 意味: 分光計装置名
- 例: `XFFTS`, `ROACH2`, `ACS`
- 型: 文字列
- 出力: optional

### 8.6 `SAMPLER`

- 意味: 個々の入力系列名
- 例: `A01`, `A02`, `IF0P0`
- 型: 文字列
- 出力: optional

### 8.7 `FRONTEND`

- 意味: 受信機・frontend 名
- 型: 文字列
- 出力: optional

### 8.8 旧 `BEAMNO`, `POLNO` を使わない理由

本仕様では `BEAMNO` と `POLNO` を正規列としません。今後の正規参照先は

- `FDNUM`
- `IFNUM`
- `PLNUM`
- `POLARIZA`
- `BACKEND`
- `SAMPLER`

です。

---

### 9. 偏波ラベル運用規約

#### 線偏波受信機

- `PLNUM=0` -> `POLARIZA='XX'`
- `PLNUM=1` -> `POLARIZA='YY'`

#### 円偏波受信機

- `PLNUM=0` -> `POLARIZA='RR'`
- `PLNUM=1` -> `POLARIZA='LL'`

#### 相関積

- 線偏波基底: `XY`, `YX`
- 円偏波基底: `RL`, `LR`

#### Stokes 量

変換済み row のみ

- `I`, `Q`, `U`, `V`

を使用します。

---

### 10. `Site`, `Efficiency`, `SpectralAxisUniform`, `DatasetInfo`

### 10.1 `Site`

- `lat_deg`: 緯度 [deg]
- `lon_deg`: 経度 [deg, east+]
- `elev_m`: 標高 [m]

用途:

- `EarthLocation` 生成
- `OBSGEO-X/Y/Z` 生成
- LSRK 補正の観測地点

### 10.2 `Efficiency`

- `beameff`
- `apereff`
- `mooneff`
- `suneff`
- `effstat`

`effstat` 許容値:

- `ASSUMED`
- `MEASURED`
- `UNKNOWN`

header には finite のものだけ書かれます。

### 10.3 `SpectralAxisUniform`

writer 全体の既定スペクトル軸です。

- `crval1_hz`
- `cdelt1_hz`
- `crpix1`
- `restfreq_hz`
- `specsys`
- `ssysobs`
- `veldef`
- `ctype1`
- `cunit1`
- `refchan`

現実装の重要制限:

- `ctype1` は `FREQ` のみを受け付ける
- `cunit1` は `Hz` のみ

したがって current writer / converter の主系列は **周波数軸保存** です。

### 10.4 `DatasetInfo`

ファイル全体既定 metadata:

- `telescope`, `observer`, `project`, `object_name`
- `radesys`, `equinox`
- `src_radesys`, `src_equinox`
- `bmaj_deg`, `bmin_deg`, `bpa_deg`
- `eff`
- `refr_included_in_corr`
- `doppler_tracking_applied`
- `spectral_axis`
- `shared_meta`

`shared_meta` はテーブル header に `HIERARCH` keyword として書かれます。

---

### 11. `SDRadioSpectralSDFITSWriter.__init__()` の引数

#### `n_chan`

- 型: int
- 必須
- 意味: 既定または最大チャンネル数の目安

#### `site`

- 型: `Site`
- 必須

#### `info`

- 型: `DatasetInfo`
- 必須

#### `store_freq_column`

- 型: bool
- 既定: `False`
- 意味: row ごとに `FREQ` ベクトル列を持つかどうか

#### `chunk_size`

- 型: `Optional[int]`
- 既定: `None`
- 意味: parts 分割モードの row 数

#### `out_basename`

- 型: `Optional[str]`
- 意味: parts 出力ファイル名の基底
- 条件: `chunk_size` を使うなら必須

#### `history`

- 型: dict / list / tuple / str / その他
- 意味: file-level history の初期値

---

### 12. `add_history()` の仕様

`add_history(key, value)` は file-level history を追加します。

内部的には次のように扱います。

- `_history is None` -> dict を新規作成
- dict -> `key -> value`
- list -> `{key:value}` を append
- その他 -> 文字列化して dict へ退避

最終的には `build_history_hdu()` で `KEY`, `VALUE` 2列の `HISTORY` BinTable になります。

---

### 13. `add_row()` の全引数と意味

以下、`add_row()` の引数を群ごとに説明します。

### 13.1 最小必須

#### `time_mjd`
- 型: float
- 意味: 観測時刻 [MJD UTC]
- 条件: finite

#### `scanid`
- 型: int
- 意味: scan 番号
- 条件: 0 以上

#### `subscan`
- 型: int
- 意味: subscan 番号
- 条件: 0 以上

#### `intgrp`
- 型: int
- 意味: integration group 番号
- 条件: 0 以上

#### `obsmode`
- 型: str
- 意味: 観測モード
- 実装: 非空文字列であること

#### `data`
- 型: `np.ndarray` など
- 意味: スペクトル配列
- 実装: 1 次元 float32 へ reshape

#### `exposure_s`
- 型: float
- 意味: 積分時間 [s]
- 条件: `>0`

#### `polariza`
- 型: str
- 意味: 物理偏波ラベル
- 条件: 非空、許容値のいずれか

### 13.2 row 補助

#### `object_name`
- 型: `Optional[str]`
- 意味: row の `OBJECT`
- `None` なら `DatasetInfo.object_name`

#### `flagrow`
- 型: int
- 意味: row 全体フラグ

#### `flag`
- 型: `Optional[np.ndarray]`
- 意味: channel-wise flag
- 条件: `data` と同長

### 13.3 beam / IF / backend / frontend

#### `fdnum`
- feed / beam 番号
- int
- 0 以上

#### `ifnum`
- IF / spectral window 番号
- int
- 0 以上

#### `plnum`
- 内部偏波番号
- int
- 0 以上

#### `backend`
- 分光計装置名
- optional string

#### `sampler`
- sampler / input 系列名
- optional string

#### `frontend`
- 受信機名
- optional string

### 13.4 heterodyne metadata

#### `obsfreq_hz`
- 信号側 sky/reference 周波数 [Hz]
- heterodyne 文脈では必須

#### `imagfreq_hz`
- image 側 sky/reference 周波数 [Hz]

#### `lo1freq_hz`, `lo2freq_hz`, `lo3freq_hz`
- 各 LO 周波数 [Hz]

#### `sideband`
- `DATA` が最終的に表す sideband
- 許容値: `USB`, `LSB`
- heterodyne 文脈では必須

#### `sb1`, `sb2`, `sb3`
- 各変換段の sideband
- 許容値: `USB`, `LSB`

#### heterodyne 文脈の判定

次のどれかがあると heterodyne 文脈が active とみなされます。

- `obsfreq_hz`
- `imagfreq_hz`
- `lo1freq_hz`, `lo2freq_hz`, `lo3freq_hz`
- `sideband`
- `sb1`, `sb2`, `sb3`

このとき少なくとも

- `OBSFREQ`
- `SIDEBAND`

が必要です。

### 13.5 calibration / temperature scale

#### `calstat`
- 行の校正状態
- 許容値: `RAW`, `TASTAR`, `TA`, `TMB`, `INTEGRATED`

#### `tempscal`
- temperature scale 文字列
- 明示すれば `normalize_tempscal()` によって正規化
- 未指定時は `calstat` を見て自動決定

#### 自動決定規則

- `calstat` に `TMB` または `TR` を含む -> `TR*`
- それ以外 -> `TA*`

#### `t_cal_k`
- `TCAL` [K]

#### `t_hot_k`
- `THOT` [K]

#### `tsys_k`
- `TSYS` [K]

#### `tau0`
- 天頂光学的厚み
- 出力列名は `TAU0`

### 13.6 sky coordinates

#### `ra_deg`, `dec_deg`
- ICRS/FK5/FK4 系の経緯ではなく天球座標 [deg]

#### `glon_deg`, `glat_deg`
- 銀河座標 [deg]

`add_row()` は

- RA/DEC があって GLON/GLAT がない
- GLON/GLAT があって RA/DEC がない

のどちらかを許し、欠けている側を補完します。両方欠けるとエラーです。

### 13.7 source block

#### `srcframe`
- source 座標ブロックの座標種別
- 許容値: `RADEC`, `GALACTIC`, `AZEL`, `AZEL_GEO`

#### `src_radesys`
- `srcframe='RADEC'` のときのみ有効
- 許容値: `ICRS`, `FK5`, `FK4`, `FK4-NO-E`, `GALACTIC`

#### `src_equinox`
- `srcframe='RADEC'` のときのみ有効

#### `src_long_deg`, `src_lat_deg`
- source 座標値 [deg]

#### 重要原則

source block が active なら `srcframe` は必須です。暗黙に `RADEC` を仮定することは許しません。

### 13.8 scan block

#### `scanframe`
- scan offset の座標種別
- 許容値: `RADEC`, `GALACTIC`, `AZEL`, `AZEL_GEO`

#### `scan_radesys`
- `scanframe='RADEC'` のときのみ有効

#### `scan_equinox`
- `scanframe='RADEC'` のときのみ有効

#### `scan_x_deg`, `scan_y_deg`
- scan offset [deg]

#### 原則

scan block が active なら `scanframe` を必須とし、暗黙既定値は認めません。

### 13.9 pointing / beam geometry

#### `az_center_deg`, `el_center_deg`
- row が最終的に表す beam-center Az/El [deg]

#### `az_enc_deg`, `el_enc_deg`
- 実測 encoder Az/El [deg]

#### `boresight_az_deg`, `boresight_el_deg`
- boresight Az/El [deg]

#### `beam_xoff_arcsec`, `beam_yoff_arcsec`
- 適用した beam offset [arcsec]

#### `beam_rot_deg`
- beam rotation [deg]

#### `az_cmd_deg`, `el_cmd_deg`
- 指令値 Az/El [deg]

#### `corr_az_deg`, `corr_el_deg`
- pointing correction [deg]

#### `calc_refr_deg`
- 計算した refraction 量 [deg]

### 13.10 Doppler / weather

#### `v_frame_mps`
- `VFRAME` [m/s]

#### `f_offset_hz`
- `FOFFSET` [Hz]

#### `tamb_c`
- 入力は摂氏 [degC]
- 出力列 `TAMBIENT` には K で保存

#### `pressure_hpa`
- 入力は hPa
- 出力列 `PRESSURE` には mmHg で保存

#### `humidity_pct`
- 入力は %
- 出力列 `HUMIDITY` には 0..1 の fraction で保存

#### `wind_spd_mps`
- `WINDSPD` [m/s]

#### `wind_dir_deg`
- `WINDDIR` [deg]

### 13.11 スペクトル WCS

#### `freq_hz`
- `FREQ` 列を明示的に与えるときの周波数ベクトル
- `store_freq_column=True` のときのみ使用

#### `restfreq_hz`
- 行の `RESTFREQ`
- 未指定時は `DatasetInfo.spectral_axis.restfreq_hz`

#### `crval1_hz`, `cdelt1_hz`, `crpix1`
- 行ごとの WCS scalar

#### `ctype1`
- 現系統では通常 `FREQ`

#### `cunit1`
- 通常 `Hz`

#### `specsys`
- 許容値: `TOPOCENT`, `LSRK`, `BARYCENT`, `GEOCENTR`, `HELIOCEN`

#### `veldef`
- 入力としては `RADIO`, `OPTICAL`, `RELATIVISTIC`, `RADI`, `OPTI`, `RELA` など
- velocity context では `RADI-OBS`, `OPTI-LSR` のような 8 文字形式へ正規化される

#### velocity context の判定

次のどれかで velocity context とみなします。

- `CTYPE1 != 'FREQ'`
- `SPECSYS == 'LSRK'`
- `VFRAME` が finite

このとき

- `RESTFREQ > 0`
- `VELDEF` がある

が必須です。

#### 現系統の注意

`SpectralAxisUniform` と converter は実質 `CTYPE1='FREQ'` を前提にしているので、現 writer の主要運用は「周波数軸保存 + 必要に応じて速度解釈 metadata を保持する」方式です。

---

### 14. `sdfits_bintable.py` の役割

### 14.1 列名の正規化

`normalize_columns=True` のとき、列名は大文字化され、FITS の case-insensitive な衝突を避けるために dedup されます。

#### 重要な挙動

- FITS 上は列名は case-insensitive
- `timestamp` と `TIMESTAMP` は衝突し得る
- dedup 時には `TIMESTAMP` を優先するロジックがある

### 14.2 DataFrame scalar 列の型推定

scalar 列は dtype に応じて FITS format へ変換されます。

- int16 相当 -> `I`
- int32 相当 -> `J`
- int64 相当 -> `K`
- float32 -> `E`
- float64 -> `D`
- bool scalar -> `B` として書く実装
- 文字列 -> `nA`

### 14.3 vector-in-cell 列の推定

DataFrame のセルに list / tuple / ndarray が入っている列は vector column とみなされます。

- 全 row 同長、かつ null なし -> fixed-length vector
- それ以外 -> VLA (`P*`) へ自動切替

### 14.4 all-null 列は落ちる

全値が `None` / `NaN` で型推定できない列は drop されます。これは optional 列の未設定を自然に列ごと省くための挙動です。

### 14.5 追加 vector 列

`extra_vector_columns` を使うと、DataFrame にない vector 列を追加できます。writer では主に `FREQ` に使われます。

### 14.6 `build_history_hdu()`

`build_history_hdu(history)` は次の形式を受け取れます。

- dict -> そのまま `KEY`, `VALUE`
- list / tuple -> dict は `i:key` に flatten、スカラーは `i` で保存
- その他 -> `0` を key として文字列化

出力は `EXTNAME='HISTORY'` の BinTableHDU です。

---

### 15. 時刻・座標・WCS 列の具体的意味

### 15.1 time columns

- `DATE-OBS`: row の絶対 UTC 時刻 (ISO)
- `DATEOBS`: 同じ情報の alias
- `TIMESTAMP`: 内部統一時刻用列
- `MJD`: MJD UTC days
- `TIME`: `DATE-OBS` からの秒オフセットとして扱う補助列

### 15.2 sky coordinates

- `RA`, `DEC`: 天球座標 [deg]
- `GLON`, `GLAT`: 銀河座標 [deg]

### 15.3 beam / pointing

- `AZIMUTH`, `ELEVATIO`: row が表す beam-center Az/El
- `BORE_AZ`, `BORE_EL`: boresight Az/El
- `BEAMXOFF`, `BEAMYOFF`: beam offset [arcsec]
- `BEAMROT`: beam rotation [deg]
- `AZ_CMD`, `EL_CMD`: 指令値
- `CORR_AZ`, `CORR_EL`: 補正量
- `CALC_REFR`: 計算 refraction

### 15.4 weather block

保存単位は入力と異なる点に注意してください。

- `TAMBIENT`: K
- `PRESSURE`: mmHg
- `HUMIDITY`: fraction 0..1
- `WINDSPD`: m/s
- `WINDDIR`: deg

### 15.5 calibration

- `CALSTAT`: `RAW`, `TASTAR`, `TA`, `TMB`, `INTEGRATED`
- `TEMPSCAL`: `TA*` / `TR*` など
- `TCAL`, `THOT`, `TSYS`, `TAU0`

---

### 16. `RESTFREQ`, `VELDEF`, `SPECSYS`, `VFRAME` の意味論

### 16.1 本 writer / converter の主系統は `CTYPE1='FREQ'`

現行 `SpectralAxisUniform` と `necst_v4_sdfits_converter.py` は、実質的に周波数軸保存を前提としています。したがって、多くの出力は

- `CTYPE1 = FREQ`
- `CUNIT1 = Hz`

です。

### 16.2 それでも `RESTFREQ` と `VELDEF` は重要

周波数軸保存でも、後段で速度解釈や LSRK 変換を行うなら、

- `RESTFREQ`
- `VELDEF`
- `SPECSYS`

は意味を持ちます。

### 16.3 velocity context での原則

この writer は、少なくとも row が velocity context を持つなら

- `RESTFREQ > 0`
- `VELDEF` あり

を要求します。

### 16.4 `VELDEF` の正規化

入力として `RADIO` を与えても、row 文脈では `RADI-OBS`, `RADI-LSR` のような 8 文字形式へ正規化されます。

ただし `PRIMARY` の既定 header 側では `SpectralAxisUniform.veldef` をそのまま書くので、row 列と header 既定値の表現が完全一致しない可能性があります。これは現実装上の注意点です。

### 16.5 `SPECSYS='LSRK'` と `FREQ` 列

`DatasetInfo.validate()` は、`spectral_axis.specsys == 'LSRK'` なら `store_freq_column=True` を要求します。これは、LSRK 周波数軸が時刻と方向に依存し、単一の固定軸だけでは表現しにくいためです。

### 16.6 `VFRAME`

`VFRAME` は optional ですが、finite であれば velocity context を活性化します。後段 coadd / standardize では、`SPECSYS='TOPOCENT'` のとき `VELOSYS` がなければ `VFRAME` を互換的に使う、という運用ポリシーが有効です。単位は m/s を前提にしてください。

---

### 17. 外部 converter / producer が writer に与えるもの

今回のアーカイブには特定の converter 本体は含まれていません。しかし、観測系や別プロセスが `SDRadioSpectralSDFITSWriter` または `fitsio.write_sdfits()` に渡す前段は、実務上きわめて重要です。ここでは、**外部 producer がこの writer に引き渡すべき契約**を整理します。

### 17.1 stream ごとの設定

multi-stream 系の producer は、各 spectrometer stream ごとに少なくとも次を明示できるのが望ましいです。

- `name`
- `fdnum`
- `ifnum`
- `plnum`
- `polariza`
- `frontend`
- `backend`
- `sampler`
- beam 情報
- spectral WCS 情報
- LO / sideband 情報
- 必要なら `channel_slice`

### 17.2 `polariza` は必須

writer コアは `POLARIZA` を first-class の必須属性として扱います。したがって外部 producer 側でも、`PLNUM` だけに依存せず、**物理偏波ラベル** を明示する必要があります。

### 17.3 stream ごとの WCS 導出

外部 producer は stream ごとに、少なくとも次を確定させてから writer に渡す必要があります。

- `CRVAL1`
- `CDELT1`
- `CRPIX1`
- `CTYPE1`
- `CUNIT1`
- `SPECSYS`
- 必要なら `RESTFREQ`, `VELDEF`

現行コア実装の主系統は `CTYPE1='FREQ'` です。ただし reader / baseline / regrid / coadd 側は `VRAD` などの速度軸も扱うため、外部 producer が速度軸で渡す場合は `RESTFREQ` と `VELDEF` を省略してはいけません。

### 17.4 heterodyne 文脈

LO / sideband 情報を持つ場合、writer コアの契約は次です。

- LO / image / sideband のどれかを出すなら heterodyne 文脈が active
- heterodyne 文脈が active なら `obsfreq_hz` は必須
- heterodyne 文脈が active なら最終 sky sideband `sideband` は必須

つまり、中途半端に LO だけ書いて `OBSFREQ` や `SIDEBAND` を欠くことは許されません。

### 17.5 writer への橋渡し

外部 producer は最終的に row ごとに少なくとも次を writer へ渡せるのが望ましいです。

- `polariza`, `fdnum`, `ifnum`, `plnum`
- `backend`, `sampler`, `frontend`
- `obsfreq_hz`, `imagfreq_hz`, `lo1/2/3freq_hz`, `sideband`, `sb1/2/3`
- `ra/dec/glon/glat`
- `az_center`, `az_enc`, `boresight`, `beam offset`, `command`, `correction`, `calc_refr`
- weather
- `restfreq_hz`, `crval1_hz`, `cdelt1_hz`, `crpix1`, `ctype1`, `cunit1`, `specsys`, `veldef`

### 17.6 `make_writer()` 相当で決めるべきこと

外部 producer が writer インスタンスを組み立てるときは、概ね次を事前に確定させる必要があります。

- 全 stream で `CUNIT1` が一意か
- `n_chan` は最大長の目安として何を採用するか
- 既定 `SpectralAxisUniform` をどの stream から取るか
- どれかの stream が `LSRK` か、または row ごとの `FREQ` ベクトル保存が必要か

特に `SPECSYS='LSRK'` を row 既定値で使う場合、コア writer は `store_freq_column=True` を要求します。

### 17.7 history への書き込み

外部 producer は `writer.add_history()` または `write_sdfits(history=...)` を使って、たとえば次を残すと有用です。

- producer 名
- config 名 / path
- output layout
- db namespace
- telescope / weather / source table 名
- channel slice
- Az/El / boresight / correction 各列の意味
- stream ごとの `fdnum`, `ifnum`, `plnum`, `polariza`, beam model version など

これらは header card よりも `HISTORY` table の方が扱いやすく、後から provenance を追いやすいです。

---

### 18. header keyword の付け方と `HIERARCH`

### 18.1 `set_meta_keyword()` の方針

`sdfits_bintable.py` の `set_meta_keyword()` は FITS header keyword を次の規則で書きます。

- すでに `HIERARCH ` で始まるならそのまま
- 8 文字以内の標準形なら通常 keyword として書く
- それ以外は `HIERARCH <key>` として書く

### 18.2 `COMMENT` / `HISTORY` はここでは扱わない

`set_meta_keyword()` は `COMMENT`, `HISTORY` という reserved card は直接書きません。history は別 Hdu へ分離する方針です。

### 18.3 `shared_meta`

`DatasetInfo.shared_meta` の値は、テーブル header に `HIERARCH` keyword として追加されます。長い key でも round-trip しやすいようにした設計です。

---

### 18.4 `fitsio.write_sdfits()` / `write_scantable()` の追加仕様

解析 I/O 系から直接 SDFITS を書く場合、`fitsio.write_sdfits()` と `write_scantable()` は writer コアと同じ BinTable 推論器を共有します。

`write_scantable()` は、table に列がなく meta にだけある重要キーワードを、書き込み前に定数列として昇格させます。主な対象は次です。

- `RESTFRQ`, `RESTFREQ`
- `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`
- `SPECSYS`, `SSYSOBS`, `VELDEF`
- `VELOSYS`, `VFRAME`
- `TEMPSCAL`, `BEAMEFF`

また `write_sdfits()` は on-disk 温度スケールを明示的に制御できます。`data_scale` をメモリ中の温度スケール、`tempscal` または `out_scale` を保存時スケールとすると、必要に応じて保存時だけ変換を行います。

$$
T_R^* = \frac{T_A^*}{BEAMEFF}
$$

$$
T_A^* = T_R^* BEAMEFF
$$

この変換は非破壊で、入力 `Scantable` 自体を書き換えず、`HISTORY` に `convert_on_write` として記録されます。

---

### 19. 文字列列と単位列の扱い

### 19.1 固定幅文字列

writer は主要文字列列に固定幅を与えます。例:

- `OBJECT`: 32
- `OBSMODE`: 24
- `CALSTAT`: 16
- `POLARIZA`: 8
- `FRONTEND/BACKEND/SAMPLER`: 32
- `CTYPE1`: 8
- `CUNIT1`: 8
- `SPECSYS`: 16
- `VELDEF`: 8

### 19.2 単位の付与

`build_single_dish_table_hdu()` へ `units` マップを渡し、列単位の `TUNITn` を付けます。

重要な例:

- `RA`, `DEC`, `GLON`, `GLAT`: deg
- `EXPOSURE`: s
- `RESTFREQ`, `OBSFREQ`, `LO1FREQ` など: Hz
- `VFRAME`: m/s
- `TAMBIENT`: K
- `PRESSURE`: mmHg
- `FREQ`: Hz

### 19.3 `CRVAL1`, `CDELT1` の unit

`CRVAL1`, `CDELT1` は row ごとの scalar metadata ですが、FITS の `TUNITn` は列ごとに1つしか付けられません。そのため、writer は **dataset 全体で `CUNIT1` が一意であること**を要求します。複数の `CUNIT1` が混在する dataset は 1 つのテーブルとして書けません。

---

### 20. 実装上の重要注意点

#### 20.1 `TIME` はまだ補助列の実装段階

本書の採用仕様では `TIME` は `DATE-OBS` からの秒オフセットですが、現 writer 実装ではまだ 0.0 固定です。絶対時刻は `TIMESTAMP` または `MJD` を優先してください。

#### 20.2 `DATE-OBS`, `DATEOBS`, `TIMESTAMP` は重複情報を持つ

これは互換性と読込側利便性を重視した設計です。

#### 20.3 row ごとの `VELDEF` と header の `VELDEF`

row 側は velocity context で 8 文字形式へ正規化されますが、header 既定値は `SpectralAxisUniform.veldef` をそのまま書くため、表現が一致しない可能性があります。

#### 20.4 weather block は入力単位と保存単位が異なる

- `tamb_c` -> `TAMBIENT[K]`
- `pressure_hpa` -> `PRESSURE[mmHg]`
- `humidity_pct` -> `HUMIDITY[fraction]`

#### 20.5 `source block` と `scan block` は暗黙既定値を許さない

block が active なら `srcframe` / `scanframe` は必須です。

#### 20.6 heterodyne 文脈では `OBSFREQ` と `SIDEBAND` が最低限必要

LO 情報を少しだけ書いて `OBSFREQ` や `SIDEBAND` がない、という中途半端な状態は許しません。

#### 20.7 外部 producer の主系統も現状は `CTYPE1='FREQ'` 主体

速度軸そのものを保存することは可能ですが、現行コア実装で最も安定している主系統は周波数軸保存です。速度窓指定や後段再配置は、`RESTFREQ` と `SPECSYS` を明示した上で reader / baseline / regrid / coadd 側に委ねるのが安全です。

#### 20.8 読込時の温度スケールは非破壊

最新版の `read_scantable()` は、on-disk の `TEMPSCAL` が `TR*` でも自動的に `TA*` へ戻しません。読込時は

- `TEMPSCAL` を尊重して保持する
- `BEAMEFF` 列を補う
- 必要な変換は viewer、coadd、write_sdfits 側で明示的に行う

という方針です。

#### 20.9 legacy 列は読込時に安全移行される

`read_scantable()` は、古い列名を公開意味論へ自動移行します。代表例は次です。

- `V_CORR_KMS [km/s]` -> `VELOSYS [m/s]`, `VFRAME [m/s]`
- `TAMB_K` -> `THOT`
- `TAU` -> `TAU0`

複数の legacy 速度列が不一致なら、読込時点でエラーにします。

#### 20.10 parts mode では `close()` が manifest 出力の完了条件

`chunk_size` を使う parts mode では、part FITS の flush だけでなく、最後に `close()` を呼んで manifest を確定させる必要があります。manifest には少なくとも次が入ります。

- format 識別子
- `n_chan`
- `store_freq_column`
- spectrum column 名
- 出力 part 一覧


---

### 21. 実務上の推奨解釈

#### 21.1 時刻

- 読み取りでは `TIMESTAMP -> MJD -> DATE-OBS + TIME -> 旧 TIME`
- 解析では絶対時刻は `TIMESTAMP` または `MJD` を使う
- `TIME` は補助列とみなす

#### 21.2 偏波

- row の物理偏波は必ず `POLARIZA`
- `PLNUM` は内部番号
- `FDNUM`, `IFNUM`, `PLNUM`, `POLARIZA`, `SAMPLER` の組で row を識別

#### 21.3 周波数 / 速度

- 現 writer の基本は `CTYPE1='FREQ'`, `CUNIT1='Hz'`
- 速度解釈は `RESTFREQ`, `SPECSYS`, `VELDEF` を見て行う
- `SPECSYS='LSRK'` なら `FREQ` row vector を保持する運用が望ましい

#### 21.4 provenance

- converter の設定や beam model、channel slice は `HISTORY` に残す
- file-level 補足は `shared_meta` と `HISTORY` を併用する

---

### 22. 短い使用例

### 22.1 最小の fixed-length 出力

```python
from sd_radio_spectral_fits import (
    SDRadioSpectralSDFITSWriter,
    Site, DatasetInfo, SpectralAxisUniform, Efficiency,
)
import numpy as np

site = Site(lat_deg=34.0, lon_deg=135.0, elev_m=50.0)
axis = SpectralAxisUniform(
    crval1_hz=115.2712018e9,
    cdelt1_hz=1.0e5,
    crpix1=1.0,
    restfreq_hz=115.2712018e9,
    specsys='TOPOCENT',
    ssysobs='TOPOCENT',
    veldef='RADIO',
    ctype1='FREQ',
    cunit1='Hz',
)
info = DatasetInfo(
    telescope='OMU1P85M',
    observer='Observer',
    project='TEST',
    object_name='Orion KL',
    spectral_axis=axis,
    eff=Efficiency(effstat='UNKNOWN'),
)
writer = SDRadioSpectralSDFITSWriter(
    n_chan=4096,
    site=site,
    info=info,
    store_freq_column=False,
)
writer.add_row(
    time_mjd=61000.0,
    scanid=0,
    subscan=0,
    intgrp=0,
    obsmode='ON',
    data=np.zeros(4096, dtype=np.float32),
    exposure_s=0.1,
    polariza='XX',
    fdnum=0,
    ifnum=0,
    plnum=0,
    ra_deg=83.8,
    dec_deg=-5.4,
    glon_deg=209.0,
    glat_deg=-19.4,
)
writer.write('example_fixed.fits')
```

### 22.2 可変長配列での出力

```python
writer = SDRadioSpectralSDFITSWriter(
    n_chan=4096,
    site=site,
    info=info,
    store_freq_column=True,
)

writer.add_row(
    time_mjd=61000.0,
    scanid=0, subscan=0, intgrp=0,
    obsmode='ON',
    data=np.zeros(4096, dtype=np.float32),
    exposure_s=0.1,
    polariza='XX',
    fdnum=0, ifnum=0, plnum=0,
    ra_deg=83.8, dec_deg=-5.4,
    glon_deg=209.0, glat_deg=-19.4,
)

writer.add_row(
    time_mjd=61000.0001,
    scanid=0, subscan=1, intgrp=0,
    obsmode='ON',
    data=np.zeros(2048, dtype=np.float32),
    exposure_s=0.1,
    polariza='YY',
    fdnum=0, ifnum=0, plnum=1,
    ra_deg=83.8, dec_deg=-5.4,
    glon_deg=209.0, glat_deg=-19.4,
)

writer.write('example_vla.fits')
```

この場合、`DATA` は `PE`、`FLAG` は `PL`、`FREQ` は `PD` になります。

### 22.3 chunked parts 出力

```python
writer = SDRadioSpectralSDFITSWriter(
    n_chan=4096,
    site=site,
    info=info,
    chunk_size=100000,
    out_basename='otf_run1',
)

# add_row() を大量に呼ぶ
# chunk_size に達すると自動で part を flush

writer.close()
```

重要な点は次です。

- parts mode では `out_basename` が必須
- `chunk_size` は row 数ベース
- 最後の端数 part と manifest を確定するには `close()` が必須

出力:

- `otf_run1_part0001.fits`
- `otf_run1_part0002.fits`
- ...
- `otf_run1_manifest.json`


---

### 23. `HISTORY` の具体例

converter 的運用では、たとえば次のような history を残すのが有用です。

- `converter`: 使用した converter 名
- `config_name`: spectrometer config 名
- `db_namespace`
- `telescope_tables`
- `channel_slice_cli`
- `azimuth_meaning`
- `bore_meaning`
- `weather_table`
- `stream_0`: stream ごとの `fdnum/ifnum/plnum/polariza/...`

これらは header card よりも `HISTORY` table の方が扱いやすく、後から provenance を追いやすいです。

---

### 24. まとめ

このパッケージが前提とする SDFITS の本質は、次のように要約できます。

1. 中核は `PRIMARY + SINGLE DISH + optional HISTORY` である。  
2. 1 row は 1 本のスペクトルと、その解釈に必要な row metadata 一式である。  
3. `FDNUM`, `IFNUM`, `PLNUM`, `POLARIZA`, `BACKEND`, `SAMPLER` が beam/pol/backend 情報の正規列である。  
4. `POLARIZA` は必須であり、`PLNUM` は内部番号にすぎない。  
5. `DATA` は fixed-length でも VLA でもよく、`FLAG` と `FREQ` もそれに追随する。  
6. `RESTFREQ`, `VELDEF`, `OBSFREQ`, `SIDEBAND` は文脈依存で必須になる。  
7. `source`, `scan`, `pointing`, `weather` は block-optional で、意味があるときだけまとめて出す。  
8. `HISTORY` は dedicated BinTable に保存する。  
9. 時刻は `TIMESTAMP` / `MJD` を主とし、`TIME` は補助列という仕様を採用する。  
10. 現 converter の主系統は `CTYPE1='FREQ'` の周波数軸保存である。  

この理解の上で reader / standardizer / coadd を設計すれば、この writer が出力する SDFITS を無理なく一貫して扱えます。

---

## 第II部: 列一覧 完全表形式
### 1. 本書の位置づけ

本書は `sdfits_writer.py` / `sdfits_bintable.py` / `necst_v4_sdfits_converter.py` の実装を基に、`sd_radio_spectral_fits` が前提とする SDFITS の**列仕様を完全表形式で整理した補遺**です。既存の `sdfits_writer_manual_ja.md` を補うものであり、ここでは特に

- `PRIMARY` HDU のキーワード
- `SINGLE DISH` BinTable の列
- `HISTORY` 拡張の構造
- `DATA` / `FLAG` / `FREQ` の fixed-length / VLA
- `writer.add_row()` に与える引数と、最終的にどの列へ書かれるか

を、列中心に整理します。

---

### 2. FITS 全体構造

この writer は基本的に次の構造を出力します。

1. `PRIMARY` HDU
2. `EXTNAME='SINGLE DISH'` の BinTableHDU
3. 必要に応じて `EXTNAME='HISTORY'` の BinTableHDU

`SINGLE DISH` の 1 row は「1本のスペクトル + その row metadata」です。したがって、同一時刻でも `FDNUM`、`IFNUM`、`PLNUM`、`POLARIZA`、`SAMPLER` などが違えば別 row になります。

---

### 3. 列ポリシーの全体像

| 区分 | 意味 | 現実装での扱い |
|---|---|---|
| always | 常に row ごとに意味を持つ列 | 無条件で出力 |
| context-required | 文脈があるとき必須 | velocity 文脈、heterodyne 文脈で出力 |
| optional | 値が1つでも meaningful なら列を出す | all-null/all-NaN なら列ごと落とす |
| block-optional | 関連する列群をまとめて出す | ブロック内のどれか1つに意味があれば全列出力 |

未設定の原則は次の通りです。

- 数値列: `np.nan`
- optional 文字列列: `None`
- always 文字列列: FITS 書き出し時に `""` へ正規化されることがある

重要なのは、**optional 列に見せかけの既定値を置かない**ことです。特に `srcframe='RADEC'` や `scan_x_deg=0.0` のような擬似既定値は、未設定判定を壊します。

---

### 4. `DATA` / `FLAG` / `FREQ` の配列列仕様

#### 4.1 最重要事項

この writer では `DATA` は**固定長配列でも可変長配列でもよい**という前提です。これは reader / standardizer / coadd の設計とも整合しています。

#### 4.2 fixed-length の場合

全 row のチャンネル数が同じ場合:

- `DATA`: `nE`
- `FLAG`: `nL`
- `FREQ`: `nD`（`store_freq_column=True` のとき）

さらに、`TDIMn=(nchan)` が付きます。

#### 4.3 VLA の場合

row ごとに長さが異なる場合:

- `DATA`: `PE`
- `FLAG`: `PL`
- `FREQ`: `PD`

ここで

- `PE`: 可変長 float32
- `PL`: 可変長 logical
- `PD`: 可変長 float64

です。

#### 4.4 整合条件

| 列 | 条件 |
|---|---|
| `DATA` | 各 row で 1 次元配列 |
| `FLAG` | 各 row で `DATA` と同じ長さ |
| `FREQ` | 各 row で `DATA` と同じ長さ |

`FLAG` が与えられない場合は all-False が自動生成されます。

#### 4.5 `NCHAN`, `NCHANSEL` の意味

| 場合 | 意味 |
|---|---|
| fixed-length | 実チャンネル数 |
| VLA | 各 row の最大長を示すヒント |

VLA のとき、真の長さは `TFORM='P*'` の可変長データから row ごとに判断します。

---

### 5. `writer.add_row()` 引数一覧

以下は `SDRadioSpectralSDFITSWriter.add_row()` が受け取る主要引数です。右端に、最終的に主として反映される列・キーワードを示します。

| 引数 | 型 | 必須 | 既定 | 主な反映先 | 備考 |
|---|---|---:|---|---|---|
| `time_mjd` | float | 必須 | なし | `MJD`, `DATE-OBS`, `DATEOBS`, `TIMESTAMP`, `TIME` | `TIME` は現実装では常に 0.0 |
| `scanid` | int | 必須 | なし | `SCAN` | 0以上 |
| `subscan` | int | 必須 | なし | `SUBSCAN` | 0以上 |
| `intgrp` | int | 必須 | なし | `INTGRP` | 0以上 |
| `obsmode` | str | 必須 | なし | `OBSMODE` | 空文字不可 |
| `data` | ndarray/float | 必須 | なし | `DATA` | 1本のスペクトル |
| `exposure_s` | float | 必須 | なし | `EXPOSURE` | 秒 |
| `polariza` | str | 必須 | なし | `POLARIZA` | 物理偏波ラベル |
| `object_name` | str/None | 任意 | `None` | `OBJECT` | 省略時は `DatasetInfo.object_name` |
| `flagrow` | int | 任意 | `0` | `FLAGROW` | row 全体フラグ |
| `flag` | ndarray/None | 任意 | `None` | `FLAG` | `DATA` と同長 |
| `fdnum` | int | 任意 | `0` | `FDNUM` | beam/feed 番号 |
| `ifnum` | int | 任意 | `0` | `IFNUM` | IF / SPW 番号 |
| `plnum` | int | 任意 | `0` | `PLNUM` | 装置内部偏波番号 |
| `backend` | str/None | 任意 | `None` | `BACKEND` | optional |
| `sampler` | str/None | 任意 | `None` | `SAMPLER` | optional |
| `frontend` | str/None | 任意 | `None` | `FRONTEND` | optional |
| `obsfreq_hz` | float/None | 任意 | `None` | `OBSFREQ` | heterodyne 文脈で必須 |
| `imagfreq_hz` | float/None | 任意 | `None` | `IMAGFREQ` | optional |
| `lo1freq_hz` | float/None | 任意 | `None` | `LO1FREQ` | optional |
| `lo2freq_hz` | float/None | 任意 | `None` | `LO2FREQ` | optional |
| `lo3freq_hz` | float/None | 任意 | `None` | `LO3FREQ` | optional |
| `sideband` | str/None | 任意 | `None` | `SIDEBAND` | heterodyne 文脈で必須 |
| `sb1`,`sb2`,`sb3` | str/None | 任意 | `None` | `SB1`,`SB2`,`SB3` | 変換段 sideband |
| `calstat` | str | 任意 | `RAW` | `CALSTAT` | enum 制約あり |
| `tempscal` | str/None | 任意 | `None` | `TEMPSCAL` | 省略時は `CALSTAT` から推定 |
| `t_cal_k` | float | 任意 | `NaN` | `TCAL` | K |
| `t_hot_k` | float | 任意 | `NaN` | `THOT` | K |
| `tsys_k` | float | 任意 | `NaN` | `TSYS` | K |
| `tau0` | float | 任意 | `NaN` | `TAU0` | 無次元 |
| `ra_deg`,`dec_deg` | float | 条件付き | `NaN` | `RA`,`DEC` | `GLON/GLAT` と相互補完 |
| `glon_deg`,`glat_deg` | float | 条件付き | `NaN` | `GLON`,`GLAT` | `RA/DEC` と相互補完 |
| `srcframe` 以下 | 各種 | 任意 | `None/NaN` | source block | ブロック起動時は frame 明示必須 |
| `scanframe` 以下 | 各種 | 任意 | `None/NaN` | scan block | 同上 |
| `az_center_deg`,`el_center_deg` | float | 任意 | `NaN` | `AZIMUTH`,`ELEVATIO` | beam-center Az/El |
| `az_enc_deg`,`el_enc_deg` | float | 任意 | `NaN` | `AZIMUTH`,`ELEVATIO` | center が無いときの fallback |
| `boresight_az_deg`,`boresight_el_deg` | float | 任意 | `NaN` | `BORE_AZ`,`BORE_EL` | optional |
| `beam_xoff_arcsec`,`beam_yoff_arcsec` | float | 任意 | `NaN` | `BEAMXOFF`,`BEAMYOFF` | 列単位の TUNIT は実装上未設定 |
| `beam_rot_deg` | float | 任意 | `NaN` | `BEAMROT` | deg |
| `az_cmd_deg`,`el_cmd_deg`,`corr_az_deg`,`corr_el_deg`,`calc_refr_deg` | float | 任意 | `NaN` | pointing block | どれか meaningful ならブロック全体 |
| `v_frame_mps` | float | 任意 | `NaN` | `VFRAME` | m/s |
| `f_offset_hz` | float/None | 任意 | `None` | `FOFFSET` | Hz |
| `tamb_c`,`pressure_hpa`,`humidity_pct`,`wind_spd_mps`,`wind_dir_deg` | float | 任意 | `NaN` | weather block | 内部で K/mmHg/fraction へ変換 |
| `freq_hz` | ndarray/None | 任意 | `None` | `FREQ` | `store_freq_column=True` のときのみ |
| `restfreq_hz` | float/None | 条件付き | `None` | `RESTFREQ` | velocity 文脈で必須 |
| `crval1_hz`,`cdelt1_hz`,`crpix1`,`ctype1`,`cunit1`,`specsys`,`veldef` | 各種 | 任意 | `DatasetInfo.spectral_axis` から継承 | WCS core | row 単位で上書き可 |

---

### 6. `SINGLE DISH` 列一覧（always）

#### 6.1 時刻・row 識別・観測基本情報

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 取り得る値 / 規約 | 意味 |
|---|---|---|---|---|---|
| `TIME` | float64 / `D` | s | always | 現実装では全 row `0.0` | `DATE-OBS` 基準の秒オフセット。採用仕様上は補助列 |
| `MJD` | float64 / `D` | d | always | MJD UTC | 絶対時刻 |
| `DATE-OBS` | str / `26A` | – | always | ISO UTC 文字列 | 各 row の絶対 UTC |
| `DATEOBS` | str / `26A` | – | always | `DATE-OBS` と同値 | 互換用別名 |
| `TIMESTAMP` | str / `26A` | – | always | 現実装では `DATE-OBS` と同値 | 読み込み側の内部統一時刻 |
| `SCAN` | int32 / `J` | – | always | 0以上 | scan ID |
| `SUBSCAN` | int32 / `J` | – | always | 0以上 | subscan ID |
| `INTGRP` | int32 / `J` | – | always | 0以上 | integration group |
| `OBJECT` | str | – | always | 任意文字列 | 代表天体名 |
| `OBSMODE` | str | – | always | telescope-specific だが空文字不可 | 観測モード |
| `EXPOSURE` | float32 / `E` | s | always | 正 | 積分時間 |
| `CALSTAT` | str | – | always | `RAW`, `TASTAR`, `TA`, `TMB`, `INTEGRATED` | キャリブレーション状態 |
| `TEMPSCAL` | str | – | always | 正規化された `TA*`, `TR*` など | 温度スケール |
| `FLAGROW` | int32 / `J` | – | always | 0/非0 | row 全体フラグ |

#### 6.2 beam / IF / polarization の中核列

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 取り得る値 / 規約 | 意味 |
|---|---|---|---|---|---|
| `FDNUM` | int32 / `J` | – | always | 0始まり推奨 | beam / feed 番号 |
| `IFNUM` | int32 / `J` | – | always | 0始まり推奨 | IF / spectral window 番号 |
| `PLNUM` | int32 / `J` | – | always | 0始まり推奨 | 装置内部偏波番号 |
| `POLARIZA` | str / `8A` | – | always | `XX`,`YY`,`RR`,`LL`,`XY`,`YX`,`RL`,`LR`,`I`,`Q`,`U`,`V` | 物理偏波ラベル |

#### 6.3 天球位置とスペクトル軸の基本列

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 取り得る値 / 規約 | 意味 |
|---|---|---|---|---|---|
| `RA` | float64 / `D` | deg | always | 有限 | 赤経 |
| `DEC` | float64 / `D` | deg | always | 有限 | 赤緯 |
| `GLON` | float64 / `D` | deg | always | 有限 | 銀経 |
| `GLAT` | float64 / `D` | deg | always | 有限 | 銀緯 |
| `CRVAL1` | float64 / `D` | `CUNIT1` | always | row ごと | 参照値 |
| `CDELT1` | float64 / `D` | `CUNIT1` | always | row ごと | チャンネル間隔 |
| `CRPIX1` | float64 / `D` | – | always | 1-based | 参照ピクセル |
| `CTYPE1` | str / `8A` | – | always | 典型的には `FREQ`, `VRAD` | 軸種別 |
| `CUNIT1` | str / `8A` | – | always | 典型的には `Hz`, `m/s` | 軸単位 |
| `SPECSYS` | str / `16A` | – | always | `TOPOCENT`,`LSRK`,`BARYCENT`,`GEOCENTR`,`HELIOCEN` | スペクトル基準系 |

#### 6.4 配列列

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 取り得る値 / 規約 | 意味 |
|---|---|---|---|---|---|
| `DATA` | `nE` または `PE` | K | always | fixed-length または VLA | スペクトル本体 |
| `FLAG` | `nL` または `PL` | – | always | `DATA` と同長 | チャンネルフラグ |

---

### 7. `SINGLE DISH` 列一覧（context-required）

#### 7.1 velocity 文脈

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 取り得る値 / 規約 | 意味 |
|---|---|---|---|---|---|
| `RESTFREQ` | float64 / `D` | Hz | `CTYPE1 != 'FREQ'` または `SPECSYS='LSRK'` または `VFRAME` あり | 正 | 静止周波数 |
| `VELDEF` | str / `8A` | – | 同上 | `RADI-OBS`,`RADI-LSR`,`OPTI-LSR` など | 速度定義 |

#### 7.2 heterodyne 文脈

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 取り得る値 / 規約 | 意味 |
|---|---|---|---|---|---|
| `OBSFREQ` | float64 / `D` | Hz | LO / sideband 情報のどれかがあるとき必須 | 正 | signal-sideband に対応する空周波数 |
| `SIDEBAND` | str / `8A` | – | 同上で必須 | `USB`,`LSB` | `DATA` が最終的にどちらの sky sideband を表すか |

---

### 8. `SINGLE DISH` 列一覧（optional）

#### 8.1 校正・大気・速度補助

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 意味 |
|---|---|---|---|---|
| `TCAL` | float32 / `E` | K | 値あり | calibration noise / reference temperature |
| `THOT` | float32 / `E` | K | 値あり | hot load 温度 |
| `TSYS` | float32 / `E` | K | 値あり | システム温度 |
| `TAU0` | float32 / `E` | – | 値あり | zenith opacity |
| `VFRAME` | float32 / `E` | m/s | 値あり | observer-to-frame 補正速度 |
| `FOFFSET` | float32 / `E` | Hz | 値あり | 追加の周波数オフセット |

#### 8.2 水平座標・ビーム幾何・装置 provenance

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 意味 |
|---|---|---|---|---|
| `AZIMUTH` | float64 / `D` | deg | 値あり | beam-center Az。無ければ encoder fallback |
| `ELEVATIO` | float64 / `D` | deg | 値あり | beam-center El。無ければ encoder fallback |
| `BORE_AZ` | float64 / `D` | deg | 値あり | boresight Az |
| `BORE_EL` | float64 / `D` | deg | 値あり | boresight El |
| `BEAMXOFF` | float64 / `D` | 実装上 TUNIT なし | 値あり | beam x オフセット。converter では arcsec で供給 |
| `BEAMYOFF` | float64 / `D` | 実装上 TUNIT なし | 値あり | beam y オフセット |
| `BEAMROT` | float64 / `D` | 実装上 TUNIT なし | 値あり | beam 回転角 |
| `FRONTEND` | str / `32A` | – | 値あり | 受信機名 |
| `BACKEND` | str / `32A` | – | 値あり | 分光計装置名 |
| `SAMPLER` | str / `32A` | – | 値あり | 分光計入力系列名 |

#### 8.3 heterodyne 補助列と周波数ベクトル

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 意味 |
|---|---|---|---|---|
| `IMAGFREQ` | float64 / `D` | Hz | 値あり | image-sideband 周波数 |
| `LO1FREQ` | float64 / `D` | Hz | 値あり | 第1 LO |
| `LO2FREQ` | float64 / `D` | Hz | 値あり | 第2 LO |
| `LO3FREQ` | float64 / `D` | Hz | 値あり | 第3 LO |
| `SB1` | str / `8A` | – | 値あり | 第1変換の sideband |
| `SB2` | str / `8A` | – | 値あり | 第2変換の sideband |
| `SB3` | str / `8A` | – | 値あり | 第3変換の sideband |
| `FREQ` | `nD` または `PD` | Hz | `store_freq_column=True` | per-row 周波数ベクトル |

---

### 9. `SINGLE DISH` 列一覧（block-optional）

#### 9.1 source block

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 意味 |
|---|---|---|---|---|
| `SRCFRAME` | str / `16A` | – | ブロック起動時 | source 座標の種別 (`RADEC`,`GALACTIC`,`AZEL`,`AZEL_GEO`) |
| `SRCRDSYS` | str / `16A` | – | `SRCFRAME='RADEC'` のとき有効 | source 座標系 |
| `SRCEQNX` | float64 / `D` | year | ブロック起動時 | source equinox |
| `SRC_LONG` | float64 / `D` | deg | ブロック起動時 | source longitude / RA-like |
| `SRC_LAT` | float64 / `D` | deg | ブロック起動時 | source latitude / Dec-like |

#### 9.2 scan block

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 意味 |
|---|---|---|---|---|
| `SCANFRAM` | str / `16A` | – | ブロック起動時 | scan offset の基準座標種別 |
| `SCANRDSYS` | str / `16A` | – | `SCANFRAM='RADEC'` のとき有効 | scan 座標系 |
| `SCANEQNX` | float64 / `D` | year | ブロック起動時 | scan equinox |
| `SCANX` | float64 / `D` | deg | ブロック起動時 | scan x offset |
| `SCANY` | float64 / `D` | deg | ブロック起動時 | scan y offset |

#### 9.3 pointing diagnostic block

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 意味 |
|---|---|---|---|---|
| `AZ_CMD` | float64 / `D` | deg | ブロック起動時 | commanded Az |
| `EL_CMD` | float64 / `D` | deg | ブロック起動時 | commanded El |
| `CORR_AZ` | float64 / `D` | deg | ブロック起動時 | pointing correction in Az |
| `CORR_EL` | float64 / `D` | deg | ブロック起動時 | pointing correction in El |
| `CALC_REFR` | float64 / `D` | deg | ブロック起動時 | 計算上の refraction 項 |

#### 9.4 weather block

| 列名 | 型 / FITS 形式 | 単位 | 出力条件 | 意味 |
|---|---|---|---|---|
| `TAMBIENT` | float32 / `E` | K | ブロック起動時 | 外気温。入力 `tamb_c` は degC、保存は K |
| `PRESSURE` | float32 / `E` | mmHg | ブロック起動時 | 入力 `pressure_hpa` を mmHg へ変換して保存 |
| `HUMIDITY` | float32 / `E` | – | ブロック起動時 | 入力 `humidity_pct` を 0..1 fraction へ変換 |
| `WINDSPD` | float32 / `E` | m/s | ブロック起動時 | 風速 |
| `WINDDIR` | float32 / `E` | deg | ブロック起動時 | 風向 |

---

### 10. `PRIMARY` HDU キーワード一覧

#### 10.1 基本キーワード

| キーワード | 型 | いつ書くか | 意味 |
|---|---|---|---|
| `SIMPLE`,`BITPIX`,`NAXIS`,`EXTEND` | FITS 基本 | 常時 | FITS 基本構造 |
| `TELESCOP` | str | 常時 | 望遠鏡名 |
| `OBSERVER` | str | 常時 | 観測者 |
| `PROJID` | str | 常時 | project / proposal ID |
| `OBJECT` | str | 常時 | 代表天体名 |
| `DATE` | str | 常時 | ファイル作成 UTC |
| `TIMESYS` | str | 常時 | `UTC` |
| `MJDSTART`,`MJDEND` | float | 常時 | ファイル全体の時刻範囲 |
| `SWNAME`,`SWVER`,`ORIGIN` | str | 常時 | writer provenance |

#### 10.2 観測地点

| キーワード | 型 | 単位 | 意味 |
|---|---|---|---|
| `SITELAT` | float | deg | 緯度 |
| `SITELON` | float | deg | 経度（east+） |
| `SITELONG` | float | deg | `SITELON` alias |
| `SITEELEV` | float | m | 標高 |
| `OBSGEO-X`,`OBSGEO-Y`,`OBSGEO-Z` | float | m | 地心直交座標 |

#### 10.3 既定スペクトル WCS

`DatasetInfo.spectral_axis` が与えられたときに書かれます。

| キーワード | 型 | 単位 | 条件 | 意味 |
|---|---|---|---|---|
| `WCSAXES` | int | – | spectral_axis あり | 1 |
| `CTYPE1` | str | – | 同上 | 既定軸種別 |
| `CUNIT1` | str | – | 同上 | 既定軸単位 |
| `CRVAL1` | float | Hz など | 同上 | 既定参照値 |
| `CDELT1` | float | Hz など | 同上 | 既定チャンネル間隔 |
| `CRPIX1` | float | – | 同上 | 既定参照ピクセル |
| `SPECSYS` | str | – | 同上 | 既定スペクトル基準系 |
| `SSYSOBS` | str | – | 同上 | 観測者基準系 |
| `REFCHAN` | int | – | 同上 | reference channel |
| `RESTFREQ`,`RESTFRQ` | float | Hz | `restfreq_hz > 0` | 既定静止周波数 |
| `VELREF` | int | – | velocity 文脈 | AIPS/casacore 互換 |
| `VELDEF` | str | – | velocity 文脈 | 既定速度定義 |

---

### 11. `SINGLE DISH` 拡張ヘッダの重要キーワード

| キーワード | 型 | いつ書くか | 意味 / 注意 |
|---|---|---|---|
| `EXTNAME='SINGLE DISH'` | str | 常時 | 拡張名 |
| `TIMESYS='UTC'` | str | 常時 | time system |
| `TELESCOP`,`OBSERVER`,`PROJID`,`OBJECT` | str | 常時 | 代表値 |
| `DATE-OBS` | str | 常時 | 先頭 row の代表時刻 |
| `NCHAN`,`NCHANSEL` | int | 常時 | fixed-length では実長、VLA では最大長 |
| `MJDSTART`,`MJDEND` | float | 常時 | テーブル全体範囲 |
| `RADESYS`,`EQUINOX` | str/float | 常時 | row の RA/DEC 解釈用 |
| `HIERARCH SRC_RADESYS`,`HIERARCH SRC_EQUINOX` | str/float | 常時 | source block の既定値 |
| `BMAJ`,`BMIN`,`BPA` | float | 値あり | ビーム形状 |
| `BEAMEFF`,`APEREFF`,`MOONEFF`,`SUNEFF`,`EFFSTAT` | 各種 | 値あり / 常時 | 効率情報 |
| `TEMPSCAL='TA*'` | str | 現実装では常時 | ここは row 列と完全同期していない点に注意 |
| `HIERARCH REFR_INC` | bool | 常時 | `CORR_*` に refraction を含むか |
| `HIERARCH DOPPLER` | bool | 常時 | online Doppler tracking の有無 |
| `TDIMn=(nchan)` | str | fixed-length の `DATA`,`FLAG`,`FREQ` に対して | VLA には付かない |

#### 11.1 `shared_meta` の扱い

`DatasetInfo.shared_meta` に入れた値は、`SINGLE DISH` ヘッダへ `HIERARCH <key>` として追加されます。

---

### 12. `HISTORY` 拡張の仕様

`build_history_hdu()` は `dict` / `list` / `tuple` / `str` から `HISTORY` 拡張を作ります。

#### 12.1 列構造

| 列名 | 型 | 意味 |
|---|---|---|
| `KEY` | 文字列 | 履歴キー |
| `VALUE` | 文字列 | 履歴値 |

#### 12.2 受理される入力形

| 入力形 | 展開例 |
|---|---|
| `dict` | `{a:1,b:2}` -> `KEY=a,b` |
| `list[dict]` | `[{'a':1},{'b':2}]` -> `0:a`, `1:b` |
| `list[str]` | `['x','y']` -> `0`, `1` |
| `str` | `'memo'` -> `0` |

#### 12.3 `add_history()` の挙動

`writer.add_history(key, value)` は file-level history を追記します。converter では、たとえば次のような provenance を記録しています。

- converter 名
- config 名 / path
- DB namespace
- `AZIMUTH/ELEVATIO` など列の意味
- beam model の説明
- stream ごとの `fdnum/ifnum/plnum/polariza`
- `obsfreq_fallback` の由来

したがって、**解析側で重要な意味論は HISTORY に明示的に残す**という運用が推奨されます。

---

### 13. 実装上の注意点

#### 13.1 `POLARIZA` は必須

`PLNUM` は内部番号であり、物理偏波の意味は持ちません。物理偏波は必ず `POLARIZA` で与えます。

#### 13.2 `BEAMNO` / `POLNO` は新規仕様では使わない

この writer の正規列は

- `FDNUM`
- `IFNUM`
- `PLNUM`
- `POLARIZA`
- `BACKEND`
- `SAMPLER`

です。

#### 13.3 `TIME` は採用仕様上は補助列

現 writer は各 row の `DATE-OBS` に絶対時刻を入れ、`TIME` は 0.0 を書きます。reader 側は `TIMESTAMP` / `MJD` / `DATE-OBS + TIME` / legacy `TIME` の順で解決するため、実務上は `TIMESTAMP` または `MJD` を正本と考えるのが安全です。

#### 13.4 `CUNIT1` はデータセット内で一意である必要がある

`CRVAL1` / `CDELT1` に対して 1 つの `TUNITn` しか付けられないため、writer は `CUNIT1` が row ごとに混在しているとエラーにします。

#### 13.5 `BEAMXOFF` / `BEAMYOFF` / `BEAMROT` の単位

converter は `BEAMXOFF/BEAMYOFF` を arcsec、`BEAMROT` を deg として供給しています。しかし現 writer の列 TUNIT はこの 3 列に付けていません。したがって、この意味は**HISTORY と converter 文書で補完することが重要**です。

---

### 14. writer / converter から見た推奨運用

1. `POLARIZA` は必ず確定値で渡す。既定値でごまかさない。  
2. line / velocity 解釈を行うなら `RESTFREQ` を必ず与える。  
3. LO / sideband 情報を出すなら、少なくとも `OBSFREQ` と `SIDEBAND` を必ず出す。  
4. source / scan / weather / pointing は、分からないなら中途半端に埋めず `None/NaN` にする。  
5. `DATA` は fixed-length でも VLA でもよいが、`FLAG` と `FREQ` の長さ整合を必ず守る。  
6. 列の意味論は `HISTORY` に残す。特に beam model、補正式、fallback 規則は残した方がよい。

---

## 第III部: reader / standardizer / coadd 契約
### 1. 本書の目的

本書は `sd_radio_spectral_fits` の reader / standardizer / coadd / baseline / calibration が、入力 SDFITS に対して**何を前提にしているか**を整理した仕様書です。単なる I/O 仕様ではなく、

- 読み込み時にどの列をどう解決するか
- 旧列名をどう公開意味論へ移すか
- 時刻をどの順で解釈するか
- 周波数軸 / 速度軸 / `RESTFREQ` / `SPECSYS` / `VELOSYS` をどう扱うか
- 可変長配列 VLA をどう扱うか
- coadd がどの row を同じグループとして束ねるか
- baseline / coadd / regrid がどの解析列を出力へ追加するか
- 温度スケールと `BEAMEFF` をどう非破壊に扱うか

を、実装に即して明文化します。

対象は主に以下のモジュールです。

- `fitsio.py`
- `rawspec.py`
- `scantable_utils.py`
- `regrid_vlsrk.py`
- `coadd.py`
- `calibrate.py`
- `baseline.py`
- `restfreq.py`
- `tempscale.py`
- `doppler.py`
- `axis.py`
- `atmosphere.py`

---

### 2. reader 側の基本前提

#### 2.1 `SINGLE DISH` テーブルが正本

reader はまず `EXTNAME='SINGLE DISH'` の BinTable を探し、そのテーブルと `PRIMARY` ヘッダを併用して `Scantable` を構成します。

- テーブル列にある値は per-row metadata
- ヘッダは dataset-level 既定値 / 補助情報

という役割分担です。

#### 2.2 VLA を正しく読めることを前提にしている

reader は `DATA` が

- fixed-length 2 次元配列
- row ごとに長さの違う list-of-arrays

の両方を許容します。したがって、writer が `TFORM='PE'/'PL'/'PD'` を使っていてもよい、というのがこのパッケージ全体の前提です。

#### 2.3 `RESTFREQ` と `RESTFRQ` は同義 alias として扱う

多くのモジュールは `RESTFRQ` / `RESTFREQ` の両方を見て値を解決します。したがって、どちらか一方しか無い入力も受理されますが、**出力時には両方を揃える方が安全**です。

#### 2.4 `read_scantable()` は legacy 列を公開意味論へ移す

最新版 `read_scantable()` は、過去の列名を読込時に安全移行します。主な規則は次です。

- `V_CORR_KMS` は km/s とみなし、m/s へ変換して `VFRAME` と `VELOSYS` に移す
- `TAMB_K` は `THOT` へ移す
- `TAU` は `TAU0` へ移す
- `VFRAME` だけが存在し `VELOSYS` が無い場合は `VELOSYS` を mirror する

複数の legacy 速度列が相互に一致しない場合は、読込時点で安全に判断できないのでエラーにします。

#### 2.5 読込時の温度スケールは非破壊

2026 方針では、read 系は `TR*` を勝手に `TA*` へ戻しません。`read_scantable()` は

- `TEMPSCAL` 列を補う
- `BEAMEFF` 列を補う
- ただしデータ値そのものは変換しない

という動作です。したがって、on-disk が `TR*` なら、in-memory でも `TR*` のまま保持されます。

#### 2.6 big-endian を native-endian へ直す

FITS 由来の配列や table が big-endian のままだと、pandas の `iloc`, `query`, `isin`, sort, groupby などで不安定になります。そこで `read_scantable()`、`merge_scantables()`、`filter_scantable()`、`find_scans()` などは numeric 列を native-endian に正規化します。

#### 2.7 `write_scantable()` は重要キーワードを列へ昇格させる

`write_scantable()` は、table に列がなく meta にだけ存在する重要キーワードを、書き込み前に定数列として補います。主な対象は次です。

- `RESTFRQ`, `RESTFREQ`
- `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`
- `SPECSYS`, `SSYSOBS`, `VELDEF`
- `VELOSYS`, `VFRAME`
- `TEMPSCAL`, `BEAMEFF`

このため、解析途中で meta 側だけを書き換えた場合でも、最終書き出し時に row metadata として揃えやすくなっています。

#### 2.8 `write_sdfits()` の保存時温度スケール変換

`write_sdfits()` は `tempscal` または `out_scale` を保存時スケール、`data_scale` をメモリ中のスケールとみなし、必要なときだけ保存時変換を行います。

$$
T_R^* = \frac{T_A^*}{BEAMEFF}
$$

$$
T_A^* = T_R^* BEAMEFF
$$

この変換は非破壊で、history には `convert_on_write` として残ります。

---

### 3. 時刻解決の契約

#### 3.1 優先順位

reader 系ユーティリティは、row 時刻を次の順に解決します。

1. `TIMESTAMP`
2. `MJD`
3. `DATE-OBS` / `DATEOBS` + `TIME`
4. 旧形式 `TIME`

これは `fitsio._resolve_table_timestamps()`、`scantable_utils._resolve_table_timestamps()`、`rawspec._resolve_mapping_timestamps()` に共通した設計です。

#### 3.2 各列の意味

| 列 | reader が期待する意味 |
|---|---|
| `TIMESTAMP` | 各 row の絶対 UTC 時刻。最優先 |
| `MJD` | 同じ時刻を MJD UTC days で表したもの |
| `DATE-OBS` / `DATEOBS` | 絶対 UTC 時刻文字列 |
| `TIME` | 標準意味では `DATE-OBS` からの秒オフセット |
| 旧 `TIME` | `TIMESTAMP/MJD/DATE-OBS` が無い場合のみ、MJD days または Unix 秒として互換解釈 |

#### 3.3 現 writer / write_sdfits との整合

現 writer 系は

- `DATE-OBS` = 各 row の絶対 UTC
- `DATEOBS` = 同値 alias
- `TIMESTAMP` = `DATE-OBS` と同値
- `MJD` = 絶対時刻
- `TIME` = 0.0

という構造で出力します。したがって、reader は `TIMESTAMP` または `MJD` から正しく復元できます。

#### 3.4 ダミー時刻の扱い

`coadd.py` は `_ensure_timestamp_column()` で時刻列を補い、解決できない場合はダミー時刻を入れることがあります。しかし `TOPOCENT` データで速度補正が必要な場合、ダミー時刻では処理できずエラーになります。

---

### 4. 座標解釈の契約

#### 4.1 velocity 補正に必要な座標

`TOPOCENT` 入力から `VELOSYS` / `VFRAME` を計算する場合、coadd と calibration は各 row の sky position を必要とします。通常は

- `RA`,`DEC`
- または frame に応じた対応列

が必要です。

#### 4.2 `RA/DEC` と `GLON/GLAT`

writer は `RA/DEC` だけ、または `GLON/GLAT` だけからでも相互変換して両方を埋めます。reader 側もこの両方があると後段処理が安定します。

#### 4.3 calibration の座標フレーム選択

`run_tastar_calibration()` は座標フレームを次の優先順で解決します。

1. 明示 `coord_frame`
2. `meta['coord_frame']`
3. 使える列の存在から推定

`RA/DEC` と `GLON/GLAT` が両方ある場合でも、metadata が曖昧なら誤解釈の余地があるので、外部 producer 側で frame を明示しておくのが安全です。

#### 4.4 site 情報

Doppler 補正の再計算には観測地点が必要です。`doppler.earth_location_from_meta()` は `SITELAT/SITELONG/SITEELEV` または `OBSGEO-X/Y/Z` を使います。したがって、`TOPOCENT` を扱う入力は site 情報を持っているのが望ましいです。

---

### 5. スペクトル軸解釈の契約

#### 5.1 reader が見る主要列

| 列 / ヘッダ | 役割 |
|---|---|
| `CRVAL1`,`CDELT1`,`CRPIX1` | 1次元 WCS の基本 |
| `CTYPE1` | `FREQ` / `VRAD` など |
| `CUNIT1` | `Hz` / `m/s` など |
| `RESTFREQ` / `RESTFRQ` | 速度軸や line 解釈に必要 |
| `SPECSYS` / `SSYSOBS` | frame 解釈 |
| `VELDEF` | radio / optical / relativistic |
| `VELOSYS` / `VFRAME` | TOPOCENT -> LSRK 変換用補正速度 |
| `FREQ` | row ごとの実周波数ベクトル。特に LSRK では有力 |

#### 5.2 周波数軸の基本式

FITS WCS から channel $i$ の周波数を作る基本式は

$$
f_i = CRVAL1 + ((i+1) - CRPIX1) CDELT1
$$

です。`axis.freq_axis_from_wcs()`、`sdfits_writer` の `freq_hz` 自動生成、`regrid_vlsrk` の row axis 再構成は、この関係を前提にしています。

#### 5.3 radio 速度と frame 補正式

周波数から radio definition の速度を作るときは

$$
v = c \frac{f_{rest} - f}{f_{rest}}
$$

を使います。また relativistic Doppler factor は

$$
\beta = \frac{v}{c}, \quad
k = \sqrt{\frac{1+\beta}{1-\beta}}
$$

で、native 周波数から LSRK 周波数へ直すときは

$$
f_{lsrk} = \frac{f_{obs}}{k}
$$

を使います。`axis.vlsrk_axis_from_freq_meta()` と `doppler.scale_frequency_wcs_by_velocity()` はこの系統の意味論に従います。

#### 5.4 `CTYPE1='FREQ'` の場合

最も単純なケースです。`RESTFREQ` がなくても単なる周波数軸としては読めます。ただし

- 速度窓指定
- Doppler 補正
- coadd / regrid の velocity-grid 化

をするなら `RESTFREQ` が必要です。

#### 5.5 `CTYPE1='VRAD'` など速度軸の場合

この場合 reader / baseline / coadd / regrid は `RESTFREQ` を必須と考えます。`RESTFREQ` が無いと、速度軸の物理的意味が確定できません。

#### 5.6 `CUNIT1` の正規化と混在

reader は `Hz`, `kHz`, `MHz`, `GHz` の表記揺れをある程度正規化しますが、1 つの dataset で `CUNIT1` が混在する入力は安全ではありません。writer は 1 テーブル内の混在を拒否し、coadd / regrid でも最終的に不整合になります。

#### 5.7 `SPECSYS='LSRK'` と `FREQ` 列

writer コアは、row の既定 frame が `LSRK` の場合、`store_freq_column=True` を要求します。理由は、LSRK 周波数軸が時刻と方向に依存し、単一の固定 WCS だけでは厳密に表しにくいからです。

---

### 6. `RESTFREQ` / `VELDEF` / `SPECSYS` / `VELOSYS` の契約

#### 6.1 `RESTFREQ`

次のどれかを行うなら必須です。

- line 観測の解釈
- velocity 軸への変換
- baseline / RMS window を速度指定する処理
- LSRK / BARY / HELIO などの frame 解釈
- coadd / regrid の velocity-grid 化

#### 6.2 `VELDEF`

writer は velocity context で `VELDEF` を要求します。reader / baseline / coadd / regrid もこれを期待します。現行コアの主系統は radio definition なので、`RADIO` ないし `RADI-LSR` 系を前提に考えるのが安全です。

#### 6.3 `SPECSYS`

`SPECSYS` は「現在その軸が属している frame」を表します。

- `TOPOCENT` を含むなら未補正データとして扱い、必要なら `VELOSYS` を計算する
- `LSRK` / `HELIO` / `BARY` などなら、原則として既に補正済みとみなす

#### 6.4 `SSYSOBS`

`SSYSOBS` は「観測時に固定だった frame」を表します。したがって、TOPOCENT 観測データを後段の `regrid` や `coadd` で LSRK へ変換しても、通常は

- `SPECSYS='LSRK'`
- `SSYSOBS='TOPOCENT'`

という関係を保ちます。最新版では `baseline` もこの意味論を壊さず、`SSYSOBS` を `SPECSYS` で上書きしません。

#### 6.5 `VELOSYS` / `VFRAME` の単位と所在

公開意味論では、`VELOSYS` / `VFRAME` は **row ごとの未適用補正速度 [m/s]** です。最新版では、baseline はとくに row-only metadata として扱い、meta 側の古い値を優先しません。

#### 6.6 rest frequency override の契約

`restfreq.apply_restfreq_override()` は、FREQ 軸と VRAD 軸で意味論を分けます。

- `CTYPE1='FREQ'` では WCS は変えず、`RESTFRQ/RESTFREQ` を差し替える
- `CTYPE1='VRAD'` では速度軸の意味が変わるので WCS も更新する

VRAD 軸での更新は、旧 rest を $f_{rest,1}$、新 rest を $f_{rest,2}$ とすると

$$
a = \frac{f_{rest,1}}{f_{rest,2}}, \quad
b = c (1 - \frac{f_{rest,1}}{f_{rest,2}}), \quad
v_2 = a v_1 + b
$$

という affine 変換に従います。

---

### 7. Standardizer と `run_velocity_regrid()`

#### 7.1 役割

`regrid_vlsrk.Standardizer` は、不均質な `Scantable` を共通速度グリッドへ再配置するクラスです。

#### 7.2 AxisSignature

Standardizer は row を次の軸署名でグループ化します。

- `NCHAN`
- `CRVAL1`
- `CDELT1`
- `CRPIX1`
- `CTYPE1`
- `RESTFREQ`
- `CUNIT1`

つまり、**軸が異なる row は別グループとして扱う**のが基本です。

#### 7.3 row 長の扱い

`DATA` が VLA でも、各 row の長さをそのまま見て軸署名を作ります。したがって VLA は Standardizer と矛盾しません。

#### 7.4 速度補正列の参照順

`Standardizer` は `v_corr_col` を持ちますが、内部では row-only に値を探し、候補として

- 指定された `v_corr_col`
- `VFRAME`
- `V_CORR_KMS`
- `VELOSYS`

を順に見ます。`V_CORR_KMS` だけが km/s、それ以外は m/s です。

#### 7.5 同一 dv と coarse dv の処理分岐

最新版では、target dv が native dv とほぼ同じなら補間経路、明確に coarse なら exact overlap-weighted rebin 経路を使います。したがって regrid は単純な 1 次補間一択ではなく、dv の関係で経路が変わります。

#### 7.6 `run_velocity_regrid()` の主引数

`run_velocity_regrid()` は次のような実務用オプションを持ちます。

- `rows`, `exclude_rows`
- `vmin_kms`, `vmax_kms`, `dv_kms`
- `v_corr_col`
- `rest_freq`
- `ch_start`, `ch_stop`
- `max_dumps`
- `drop_allnan_rows`
- `history_tag`

#### 7.7 `run_velocity_regrid()` の出力契約

regrid 出力は共通速度グリッドへ乗った解析生成物です。したがって出力では

- `CTYPE1='VRAD'`
- `CUNIT1='m/s'`
- `SPECSYS='LSRK'`
- `SSYSOBS` は観測時 frame を保持。たとえば TOPOCENT 観測を regrid しても、通常は `SSYSOBS='TOPOCENT'` のまま
- `RESTFRQ`, `RESTFREQ` を再設定

し、未適用補正列は drop します。meta 側でも `VELOSYS`, `VFRAME`, 指定 `v_corr_col` は除去されます。

#### 7.8 heterogeneous `RESTFREQ` はそのままでは許されない

regrid 出力は 1 本の velocity grid を共有するため、row 間で `RESTFRQ/RESTFREQ` が不均一な場合はエラーです。別線ごとに分けるか、`rest_freq=` を明示して揃える必要があります。

---

### 8. coadd が前提とすること

#### 8.1 入力時の正規化

`run_velocity_coadd()` は入力を読んだ後、概ね次の順で正規化します。

1. timestamp 列を解決 / 追加
2. 標準列を補完
3. 必要なら `RESTFREQ` 上書き
4. 周波数 WCS 単位を正規化
5. `SPECSYS` 解決
6. `TOPOCENT` なら `VELOSYS/VFRAME` を解決または計算
7. 全入力の `RESTFREQ` 一致を確認
8. `TEMPSCAL` / `BEAMEFF` を補完

#### 8.2 coadd 入力の必須条件

実務上、velocity coadd に安全に入れるには少なくとも次が必要です。

| 項目 | 必須度 | 理由 |
|---|---|---|
| `DATA` | 必須 | スペクトル本体 |
| `CRVAL1`,`CDELT1`,`CRPIX1`,`CTYPE1`,`CUNIT1` | 必須 | 軸再構成 |
| `RESTFREQ` / `RESTFRQ` | ほぼ必須 | velocity 化に必要 |
| `SPECSYS` または `SSYSOBS` | 必須 | frame 解釈 |
| 有効な時刻列 | `TOPOCENT` では必須級 | `VELOSYS` 計算 |
| 有効な sky coordinate | `TOPOCENT` では必須級 | `VELOSYS` 計算 |

#### 8.3 grouping 規則

`group_mode` は主に次を取ります。

- `position`
- `scan`
- `intgrp`

`group_mode='scan'` のとき、coadd は単に `SCAN` だけで束ねません。実際には

- `SCAN`
- `OBSMODE`
- `FDNUM`
- `IFNUM`
- `PLNUM`
- `_INPUT_ID`

で group を作ります。これにより、異なる beam / IF / polarization / 入力ファイルが同じ `SCAN` 番号を持っていても誤混合しません。

#### 8.4 重みづけ mode

`mode` は主に

- `uniform`
- `rms`
- `rms_weight`

を取ります。`rms` 系では基本的に重みは

$$
w_n \propto \frac{1}{\sigma_n^2}
$$

で、`sigma_n` は `baseline_vwin`、`rms_vwin`、または既存 `BSL_RMS` から評価されます。

#### 8.5 `baseline_vwin` と `rms_vwin` は同時指定不可

最新版では、`baseline_vwin` と `rms_vwin` を同時に指定するとエラーです。意味が異なるためです。

- `baseline_vwin`: ベースライン差し引きと RMS 評価の両方に使う
- `rms_vwin`: ベースライン差し引きはせず、RMS 評価だけに使う

#### 8.6 QC mode

`coadd_qc` は `uniform` mode では使えません。QC は本質的に統計重みづけを伴うため、`mode='rms'` か `mode='rms_weight'` と組み合わせる必要があります。また QC を使うときは `rms_vwin` が必要です。

#### 8.7 `weight_zero_policy`

RMS 由来の重みが 0 または非正になったときの扱いは次の 3 通りです。

- `error`
- `drop`
- `impute_median`

coadd はこの設定に応じて、停止、該当スペクトルの除外、あるいは median 重みの代入を行います。

#### 8.8 `SSYSOBS` の保持と mixed-group 保護

coadd 出力では `SPECSYS` は LSRK へ更新されますが、`SSYSOBS` は観測時 frame を保持します。したがって TOPOCENT 観測由来の行を coadd した出力では、通常

- `SPECSYS='LSRK'`
- `SSYSOBS='TOPOCENT'`

です。

また、1 本の coadd 出力 row に対して `SSYSOBS` 候補が複数混ざる場合は、曖昧な metadata を作らないためエラーで停止します。これは数値平均の仕様変更ではなく、frame metadata の安全化です。

#### 8.9 `normalize_if_mixed`, `out_scale`, `sigma_scale`

同一 coadd group 内で `BEAMEFF` が混在すると、物理的にはそのまま平均するのが危険です。最新版では

- `normalize_if_mixed='auto'`: 必要なら内部で `TR*` 正規化して coadd
- `normalize_if_mixed='never'`: 危険だが正規化せずに進む

です。`out_scale='TR*'` を指定した場合も、最終出力は

$$
T_R^* = \frac{T_A^*}{BEAMEFF}
$$

で作られます。なお `sigma_scale` は現状 `TA*` 以外未実装です。

#### 8.10 post-coadd baseline

coadd には post-coadd baseline 機能があります。主な引数は次です。

- `post_baseline_mode`
- `post_baseline_vwin`
- `post_baseline_poly`
- `post_baseline_iter_max`
- `post_baseline_iter_sigma`

`post_baseline_mode='inherit_all'` を使うと、前段 baseline 設定を継承できます。ただし継承元が必要なので、`baseline_vwin` を与えていないと使えません。

#### 8.10 coadd が追加する代表列

coadd 出力 table には、代表的に次の解析列が追加または再設定されます。

- `TEMPSCAL`
- `BEAMEFF`
- `VELOSYS_OBS`, `VFRAME_OBS`, または指定補正列の `*_OBS`
- `VELOSYS=0`, `VFRAME=0`
- `WGT_MODE`, `WGT_SRC`, `WGT_VWIN`, `WGT_ZERO_POLICY`
- post-baseline を行った場合の `BSL_*`

`VELOSYS` / `VFRAME` を 0 にする理由は、**出力時点で共通 frame へ補正済み** だからです。観測時の補正値は `*_OBS` に退避します。

#### 8.11 coadd 出力の軸

`axis_type='vel'` では出力を `VRAD` / `m/s` / `LSRK` に揃えます。`axis_type='freq'` では周波数軸出力も可能ですが、いずれにせよ共通グリッドへ再投影された解析生成物であることに変わりはありません。

---

### 9. calibration / baseline が前提とすること

#### 9.1 `run_tastar_calibration()` の対象

`run_tastar_calibration()` は raw FITS, `Scantable`, あるいは raw dict を受け取り、HOT / OFF / ON から `Ta*` を生成します。現行の対象は **frequency-axis の TOPOCENT または LSRK 入力** です。

#### 9.2 `run_tastar_calibration()` の主要引数

主な引数は次です。

- `t_hot_k`
- `ch_range`
- `vlsrk_range_kms`
- `coord_frame`
- `vcorr_chunk_sec`
- `dtype`
- `rest_freq`
- `gain_mode`
- `tau_zenith`
- `t_atm_model`
- `t_atm_delta_k`
- `t_atm_eta`

#### 9.3 TOPOCENT 入力の速度補正

TOPOCENT で `VELOSYS` / `VFRAME` が row に無い場合、calibration は時刻・座標・site 情報から補正を再計算します。足りなければ止まります。したがって TOPOCENT raw を扱う producer は、少なくとも次のどちらかを満たす必要があります。

- row ごとの `VELOSYS/VFRAME` を持つ
- それらを再計算できるだけの時刻・座標・site 情報を持つ

#### 9.4 calibration の基本式

basic 1-temperature モデルでは、各 row に対して

$$
G = \frac{T_{cal}}{HOT - OFF}
$$

$$
T_A^* = (ON - OFF) G
$$

という形で `Ta*` を作ります。ここで `HOT`、`OFF`、`ON` は time interpolation 後の各スペクトルです。

#### 9.5 大気モデルを使う場合の `T_cal`

大気モデルを使うときは、`atmosphere.compute_t_cal_array()` により

$$
X = \frac{1}{\sin(El)}
$$

$$
T_{cal} = T_{hot} e^{	au X} - T_{atm} (e^{	au X} - 1) - T_{bg}
$$

を各 row で評価します。`T_atm` は `offset` モデルまたは `ratio` モデルで与えます。

#### 9.6 calibration 出力で保証される列

calibration 後の table では、少なくとも次が整えられます。

- `THOT`
- `TCAL`
- 必要なら `TAU0`
- `TEMPSCAL='TA*'`
- `RESTFRQ`, `RESTFREQ`
- `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`, `SPECSYS`, `SSYSOBS`, `VELDEF`
- TOPOCENT 入力なら `VELOSYS`, `VFRAME`

#### 9.7 `recalibrate_tastar()`

`recalibrate_tastar()` は、既存 `Ta*` に対して新しい `tau` や `T_atm` モデルを適用し直す関数です。

- `new_tau=None` なら basic 1-temperature へ undo 的に戻す
- `new_t_surface_k` が無ければ metadata / table から探す
- 時変外気温があれば配列として引き継ぐ

#### 9.8 `run_baseline_fit()` の主要引数

最新版 `run_baseline_fit()` は次の実務用引数を持ちます。

- `rows`, `exclude_rows`
- `vwin`
- `line_vwin`
- `poly_order`
- `iter_max`, `iter_sigma`
- `max_dumps`
- `ch_start`, `ch_stop`
- `v_corr_col`
- `rest_freq`
- `apply`
- `bsl_overwrite`

#### 9.9 baseline の window 意味論

baseline は `vwin` を line-free window、`line_vwin` を signal window として扱います。`vwin` の解釈 frame は VLSR で、row ごとの補正値は row-only の `VELOSYS` / `VFRAME` から取ります。

#### 9.10 baseline が追加する `BSL_*`

baseline 実行後には、主に次の列が追加されます。

- `BSL_DONE`
- `BSL_APPLIED`
- `BSL_STAGE`
- `BSL_POLY`
- `BSL_WINF`
- `BSL_RMS`
- `BSL_STAT`
- `BSL_NUSED`
- `BSL_COEF`
- `BSL_SCALE`

`BSL_SCALE` は、その baseline metadata がどの温度スケールで評価されたかを示します。

#### 9.11 `apply=False` の意味

`apply=False` にすると、データ値そのものは差し引かず、baseline fit と `BSL_*` だけを記録します。つまり、RMS 評価や品質管理だけ行いたい場合に使えます。

#### 9.12 channel slice 時の WCS 更新

baseline は `ch_start`, `ch_stop` を使うと、実際の data slice に加えて WCS も更新します。基本則は

$$
CRVAL1_{new} = CRVAL1_{old} + ch_{start} CDELT1
$$

で、`CRPIX1` は据え置きです。

---

### 10. helper utilities が前提とすること

#### 10.1 `merge_scantables()`

`merge_scantables()` は複数 `Scantable` を高速結合します。最新版では `shift_scan_id=True` が既定で、入力ごとに `SCAN` が衝突しないよう自動シフトします。

#### 10.2 `update_metadata()`

`update_metadata()` は指定列を row 単位で更新するユーティリティです。`force=False` なら既存値を尊重し、`rows=` で対象行を絞れます。

#### 10.3 `set_beameff()`

`set_beameff()` は `Scantable` の `BEAMEFF` を row 単位または一括で設定します。coadd や write_sdfits の温度スケール変換前に、`BEAMEFF` を明示したいときに使います。

#### 10.4 `calc_mapping_offsets()`

`calc_mapping_offsets()` は座標列から map offset を計算するユーティリティで、`frame`, `projection`, `unit`, `cos_mode` を持ちます。後段で position grouping や map 表示を行うときの補助になります。

---

### 11. reader / standardizer / coadd から見た推奨 SDFITS

#### 11.1 最低限これを満たすと安全

| 項目 | 推奨 |
|---|---|
| 時刻 | `TIMESTAMP` と `MJD` を両方持つ |
| WCS | `CRVAL1`,`CDELT1`,`CRPIX1`,`CTYPE1`,`CUNIT1` を row に持つ |
| line 情報 | `RESTFREQ` と `RESTFRQ` を両方持つ |
| frame | `SPECSYS` と `SSYSOBS` を揃える |
| velocity 補正 | `TOPOCENT` なら `VELOSYS` か有効時刻+座標を持つ |
| beam/pol | `FDNUM`,`IFNUM`,`PLNUM`,`POLARIZA` を持つ |
| 装置 provenance | `BACKEND`,`SAMPLER`,`FRONTEND` を可能なら持つ |
| 温度スケール | `TEMPSCAL` と `BEAMEFF` を持つ |
| 配列 | `DATA` と `FLAG` の長さを必ず一致させる |

#### 11.2 LSRK 入力で特に重要な点

writer 側で `SPECSYS='LSRK'` を使う場合、`store_freq_column=True` を使って `FREQ` を実際の row 周波数として持たせる設計が安全です。現 writer もこの条件を要求します。

#### 11.3 `TEMPSCAL` / `BEAMEFF`

coadd と write_sdfits は `TEMPSCAL` と `BEAMEFF` を整理しながら進みます。`TR*` へ変換する場合、どの row / beam / pol にどの `BEAMEFF` が適用されるかが明確である必要があります。

---

### 12. 失敗しやすいパターン

#### 12.1 `TOPOCENT` なのに時刻が不完全

- `TIMESTAMP` なし
- `MJD` なし
- `DATE-OBS` も信用できない
- `VELOSYS` / `VFRAME` もない

この場合、coadd も calibration も速度補正を計算できません。

#### 12.2 `RESTFREQ` が無いのに velocity 窓を使う

baseline / calibrate / coadd / regrid のいずれかで失敗します。

#### 12.3 multi-beam なのに `FDNUM` が無い

scan grouping、position grouping、`set_beameff()` 適用で曖昧さが残ります。

#### 12.4 `POLARIZA` を書かず `PLNUM` だけで済ませる

後段解析で物理偏波の意味が失われます。

#### 12.5 row ごとに `CUNIT1` を混在させる

writer で弾かれるか、reader / coadd / regrid で扱いにくくなります。

#### 12.6 `TR*` を暗黙に `TA*` とみなす

最新版は読込時に非破壊なので、暗黙変換を仮定すると重み評価や最終出力温度スケールで齟齬が出ます。

#### 12.7 heterogeneous `RESTFREQ` のまま 1 回の velocity regrid を掛ける

regrid は 1 本の velocity grid を作るので、異なる線をそのまま混ぜると失敗します。

---

### 13. 最終結論

このパッケージ全体の観点から見ると、良い SDFITS とは次の条件を満たすものです。

1. `SINGLE DISH` の各 row が自己完結した metadata を持つ。  
2. `DATA` は fixed-length でも VLA でもよいが、WCS と時刻が解決できる。  
3. `RESTFREQ`、`SPECSYS`、`TIMESTAMP/MJD` が明確である。  
4. `FDNUM`、`IFNUM`、`PLNUM`、`POLARIZA` により beam / IF / pol が曖昧でない。  
5. `TOPOCENT` データでは `VELOSYS` を計算できるだけの時刻・座標・site 情報がある。  
6. `TEMPSCAL` と `BEAMEFF` の意味が明確で、必要な変換は非破壊かつ明示的に行われる。  
7. 解析生成物として baseline / coadd / regrid が新しい WCS と解析列を安全に上書きできる。  

この契約を守る限り、reader、standardizer、coadd、baseline、calibration、writer は相互に整合して動作しやすくなります。

---


## 補遺A: `run_velocity_regrid()` 詳細パラメータ補遺

`run_velocity_regrid()` は、`Standardizer` を表に出した CLI / notebook 向け実務ラッパーです。主な引数の意味を、単位・入力・出力の観点で改めて整理します。

### A.1 入力選択

- `input_data`: `Scantable` または FITS path
- `rows`: 対象 row selector
- `exclude_rows`: 除外 row selector
- `max_dumps`: 選択後の先頭からさらに切る上限

`rows` と `exclude_rows` は同時指定できません。`max_dumps` は row 選択後に作用します。

### A.2 目標速度グリッド

- `vmin_kms`, `vmax_kms`: 出力速度範囲 [km/s]
- `dv_kms`: 出力チャネル幅 [km/s]

生成される目標 grid は、内部では `VGrid` として保持され、出力 WCS は

$$
CRVAL1 = 1000 v_{min}
$$

$$
CDELT1 = 1000 dv
$$

$$
CRPIX1 = 1
$$

として `m/s` 単位で保存されます。

### A.3 補正速度列

- `v_corr_col`: TOPO 入力に対して使う row-wise 補正列名

実装上の検索順は、指定列、`VFRAME`、`V_CORR_KMS`、`VELOSYS` です。ただし公開意味論としては `VELOSYS` / `VFRAME` を m/s で持つのが標準です。

### A.4 rest frequency の取り扱い

- `rest_freq`: 出力または解釈用に明示上書きしたい静止周波数 [Hz]

`rest_freq` を与えない場合は row / meta から `RESTFRQ` / `RESTFREQ` を集めます。不均一なら停止します。

### A.5 channel pre-slice

- `ch_start`, `ch_stop`: native channel の先頭 / 末尾を事前に切る

これは「まず元の周波数軸を切ってから velocity regrid する」ことを意味します。WCS は

$$
CRVAL1_{new} = CRVAL1_{old} + ch_{start} CDELT1
$$

で更新されます。

### A.6 `drop_allnan_rows`

regrid 後に全チャネル NaN になった row を drop するかどうかです。map 前処理などで便利ですが、row 数が変わるので provenance 上は注意が必要です。

### A.7 history

- `history_tag`: history 辞書へ書くキー名

history には、入力、grid、row selector、slice、`v_corr_col`、`rest_freq`、`drop_allnan_rows` などがまとまって入ります。

---

## 補遺B: `run_velocity_coadd()` 詳細パラメータ補遺

`run_velocity_coadd()` は、速度再配置、grouping、重みづけ平均、必要なら post-baseline までを一括で行う高機能関数です。ここでは主な引数を機能群ごとに整理します。

### B.1 入力と row 選択

- `inputs`: `Scantable` または path の列
- `rows`, `exclude_rows`
- `output_path`
- `overwrite`

複数入力はまず 1 つの解析対象へ統合されますが、`_INPUT_ID` が内部的に保持されるので、scan grouping で不用意に混ざりません。

### B.2 grouping 関連

- `group_mode`: `position`, `scan`, `intgrp`
- `pos_col`: 位置 ID 列名
- `pos_tol_arcsec`: position grouping の許容距離 [arcsec]

`position` grouping では `pos_col` を優先し、必要に応じて許容距離 grouping を行います。

### B.3 frame / WCS 関連

- `v_corr_col`
- `coord_frame`
- `vcorr_chunk_sec`
- `vmin`, `vmax`, `dv`
- `axis_type`
- `rest_freq`
- `allow_outside_overlap`
- `ch_start`, `ch_stop`
- `block_size`
- `max_dumps`

`vcorr_chunk_sec` は TOPO 補正再計算時の decimation に使われます。時刻列が長い OTF では速度補正計算の高速化に有効です。

### B.4 重みづけ関連

- `mode`
- `rms_vwin`
- `rms_poly`
- `rms_bin`
- `weight_zero_policy`
- `coadd_qc`
- `line_vwin`

`rms_poly` は RMS 評価用の事前 baseline fitting 次数です。coadd 前の本データを必ずしも引くわけではなく、重み算出専用に使われることがあります。

### B.5 baseline 関連

- `baseline_vwin`
- `baseline_poly`
- `baseline_iter_max`
- `baseline_iter_sigma`
- `post_baseline_mode`
- `post_baseline_vwin`
- `post_baseline_poly`
- `post_baseline_iter_max`
- `post_baseline_iter_sigma`

前段 baseline を group 内各スペクトルに掛ける経路と、coadd 後に代表スペクトルへ掛ける経路は別です。後者で作られる `BSL_*` は `BSL_STAGE='post_coadd'` になります。

### B.6 温度スケール関連

- `out_scale`
- `normalize_if_mixed`
- `beameff_tol`
- `sigma_scale`

`out_scale='TR*'` を選んだ場合の最終出力は、row 代表 `BEAMEFF` と共に保存されます。混在 group を `TR*` 正規化した場合、最終出力 `BEAMEFF` は 1 に揃うことがあります。

### B.7 coadd が追加する列の詳細

重みづけ関連では代表的に次が出ます。

- `WGT_MODE`: `uniform`, `rms`, `rms_weight`, `QC:...`
- `WGT_SRC`: `baseline_vwin`, `rms_vwin`, `BSL_RMS`, `uniform`
- `WGT_VWIN`: 実際に使った window 文字列
- `WGT_ZERO_POLICY`: 0 重み対策

補正速度関連では

- `VELOSYS_OBS`
- `VFRAME_OBS`
- 指定補正列の `*_OBS`

が退避列として残ります。

### B.8 coadd 出力の物理的意味

coadd 出力は「観測 dump のそのままの保存」ではなく、「共通 frame と共通 grid に再投影された解析生成物」です。したがって、観測時の raw WCS を完全保存したい用途ではなく、後段解析・可視化・map 作成を主目的とする出力です。

---

## 補遺C: `run_baseline_fit()` 詳細パラメータ補遺

### C.1 row 選択と安全策

- `rows`, `exclude_rows`
- `max_dumps`
- `bsl_overwrite`
- `on_fail`

既に `BSL_*` がある入力に対して `bsl_overwrite='error'` を使うと、重ね掛けを防げます。

### C.2 baseline 窓

- `vwin`: baseline fit に使う line-free window
- `line_vwin`: signal window

実装上は `vwin` が主 window で、`line_vwin` は signal 域の明示や可読性向上に使われます。

### C.3 robust 化

- `iter_max`
- `iter_sigma`

外れ値除去付き baseline は、反復ごとに残差の散らばりを評価し、しきい値超過点を落として再 fit します。

### C.4 `apply` と `BSL_APPLIED`

- `apply=True`: 実データから baseline を差し引く
- `apply=False`: 差し引かず、係数と RMS だけ記録する

したがって `BSL_DONE=True` かつ `BSL_APPLIED=False` という状態は、「fit は済んだがデータ値は変えていない」を意味します。

### C.5 `BSL_COEF` の shape

`BSL_COEF` は vector-in-cell 列として保存されるので、fixed-length でも VLA でも書けます。reader は list-of-arrays として受けてもよい設計です。

### C.6 温度スケールとの関係

baseline metadata は `BSL_SCALE` を持ちます。例えば `TR*` のまま read したデータに baseline を掛ければ、その `BSL_RMS` も `TR*` の散らばりです。`TA*` と仮定して後段で使わないよう注意が必要です。

---

## 補遺D: `run_tastar_calibration()` / `recalibrate_tastar()` 詳細パラメータ補遺

### D.1 `gain_mode`

現行実装では `independent` と `hybrid` の系統があり、`hybrid` では HOT 時刻へ OFF を内挿してゲイン分母をクリーンに作ってから ON 時刻へ再内挿します。定在波や時間変化の影響を抑えたいときに有利です。

### D.2 `vlsrk_range_kms` と `ch_range`

- `vlsrk_range_kms`: LSRK 速度窓で元周波数軸を事前 slice
- `ch_range`: 生 channel index で事前 slice

両方与えた場合、実装は両者を順に反映します。したがって、最終的な `CRVAL1` と `NAXIS1` は切り出し後の値になります。

### D.3 `tau_zenith='auto'`

`tau_zenith='auto'` では、`TAU0`, `TAU`, `OPACITY` の順で table / meta を探します。見つからなければ大気モデルは適用できません。

### D.4 `new_tau=None` の意味

`recalibrate_tastar()` で `new_tau=None` は、旧 `T_cal` を剥がして basic 1-temperature の `T_cal=T_hot` に戻す意味です。つまり「追加大気補正をやめる」方向の再較正です。

---

## 補遺E: 単位・符号・WCS 更新則の要約

### E.1 単位

- `RESTFRQ`, `RESTFREQ`, `OBSFREQ`, `FREQ`, `CRVAL1`, `CDELT1`: Hz 系
- `VELOSYS`, `VFRAME`: m/s
- `vmin_kms`, `vmax_kms`, `dv_kms`: km/s
- `TAMBIENT`: K
- `PRESSURE`: mmHg
- `HUMIDITY`: 0 から 1 の fraction
- `BEAMXOFF`, `BEAMYOFF`: arcsec

### E.2 周波数軸 slice の更新則

channel slice を掛けたときの安全な WCS 更新は

$$
CRVAL1_{new} = CRVAL1_{old} + ch_{start} CDELT1
$$

$$
CRPIX1_{new} = CRPIX1_{old}
$$

です。これは `axis.wcs_slice_channels()`、`baseline`、`regrid` が共有する方針です。

### E.3 速度補正列の符号の意味

公開意味論では、`VELOSYS` / `VFRAME` は「その row にまだ適用されていない observer to rest-frame の補正速度」です。coadd / regrid でこれを消化した後は、出力側で 0 にされるか、削除されます。

### E.4 `RESTFRQ` と `RESTFREQ`

実装上は alias ですが、外部との互換のため両方をそろえて書くのが最も安全です。

---

## 付記

本書は Markdown 版です。今後必要であれば、この統合版をもとに PDF 版または章立てをさらに整理した配布版へ展開できます。
