# sd_radio_spectral_fits SDFITS 最終完全版説明書
## 0. この文書について
本書は、`sd_radio_spectral_fits` が前提とする SDFITS の仕様を、**writer / converter / reader / standardizer / coadd** の全体像が一冊で分かるように統合した最終版です。

対象は主に以下です。

- `sdfits_writer.py`
- `sdfits_bintable.py`
- `necst_v4_sdfits_converter.py`
- `fitsio.py`
- `rawspec.py`
- `scantable_utils.py`
- `restfreq.py`
- `regrid_vlsrk.py`
- `coadd.py`
- `calibrate.py`
- `baseline.py`

本書では、次の3本の文書を統合しています。

1. writer 本体仕様
2. 列一覧の完全表形式版
3. reader / standardizer / coadd 側の前提契約版

重複する内容は完全には削らず、**読み手がどの部から読んでも必要情報に到達できる**ことを優先しています。ただし、構造は次のように整理しました。

- **第I部**: writer / converter が書く SDFITS の本体仕様
- **第II部**: 列・キーワード・履歴の完全表
- **第III部**: reader / standardizer / coadd が期待する契約

特に重要な前提は以下です。

- `DATA` は **fixed-length でも VLA でもよい**
- `FLAG` と `FREQ` は `DATA` と row ごとに整合していなければならない
- `POLARIZA` は物理偏波ラベルとして必須
- `FDNUM` / `IFNUM` / `PLNUM` は row 識別用の中核列
- 時刻解決は `TIMESTAMP -> MJD -> DATE-OBS + TIME -> 旧 TIME` を採用する
- `RESTFREQ` / `VELDEF` は常設ではなく、文脈依存で必須

---

## 第I部: writer が前提とする SDFITS 本体仕様

## sd_radio_spectral_fits が前提とする SDFITS 仕様書

### 1. 目的

本書は、`sd_radio_spectral_fits` パッケージのうち、主に次の3ファイルに基づいて、このパッケージが前提とする SDFITS の構造と、`sdfits_writer.py` が想定している列・ヘッダ・履歴・時刻・偏波・beam・WCS の仕様を、実装ベースで詳細に説明するものです。

- `sdfits_writer.py`
- `sdfits_bintable.py`
- `necst_v4_sdfits_converter.py`

ここでいう SDFITS は、古典的な標準をそのまま再現するというより、**単一鏡スペクトル解析に必要な provenance と row metadata を十分に保持するための SDFITS-like 実装**です。したがって、本書では「一般論としての SDFITS」ではなく、**このパッケージが実際に書く FITS の仕様**を定義します。

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

### 17. `necst_v4_sdfits_converter.py` が writer に与えるもの

converter は、この writer の実務上の代表的利用例です。

### 17.1 stream ごとの設定

converter は TOML から各 spectrometer stream を読みます。各 stream には主に次を持ちます。

- `name`
- `fdnum`
- `ifnum`
- `plnum`
- `polariza`
- `frontend`
- `backend`
- `sampler`
- `beam` ブロック
- `frequency_axis` ブロック
- `local_oscillators` ブロック
- `channel_slice`

### 17.2 `polariza` は必須

`_normalize_stream_block()` は stream 設定に `polariza` がないとエラーにします。つまり converter 側でも `POLARIZA` 必須方針を採っています。

### 17.3 `derive_stream_wcs()`

converter は stream ごとに `StreamWCS` を作ります。現在サポートするのは主として

- `CTYPE1='FREQ'`
- `SPECSYS='TOPOCENT'` または `LSRK`

です。

### 17.4 heterodyne 文脈

LO / sideband 情報があるとき、最終 sideband が導出できなければエラーです。また `obsfreq_hz` が明示されないが `restfreq_hz > 0` なら、互換目的で `RESTFREQ -> OBSFREQ` fallback を許しています。ただしこれは actual sky tuning の厳密表現ではなく、writer 互換用 metadata です。

### 17.5 writer への橋渡し

`write_sdfits()` は最終的に `writer.add_row()` へ、

- `polariza`, `fdnum`, `ifnum`, `plnum`
- `backend`, `sampler`, `frontend`
- `obsfreq_hz`, `imagfreq_hz`, `lo1/2/3freq_hz`, `sideband`, `sb1/2/3`
- `ra/dec/glon/glat`
- `az_center`, `az_enc`, `boresight`, `beam offset`, `command`, `correction`, `calc_refr`
- `weather`
- `restfreq_hz`, `crval1_hz`, `cdelt1_hz`, `crpix1`, `ctype1`, `cunit1`, `specsys`, `veldef`

を渡します。

### 17.6 `make_writer()` の重要仕様

converter の `make_writer()` は

- 全 stream の `CUNIT1` が一致すること
- `n_chan` には stream の最大チャンネル数を使うこと
- 既定 `SpectralAxisUniform` は先頭 stream から取ること
- どれかの stream が `LSRK` か、または `store_freq_column` を要求したら `store_freq_column=True`

という方針をとります。

### 17.7 history への書き込み

converter は `writer.add_history()` を使って、たとえば次を残します。

- converter 名
- config 名 / path
- output layout
- db namespace
- telescope table 名
- channel slice
- Az/El / boresight / correction 各列の意味
- weather table 名
- source / scan block を省略した理由
- stream ごとの `fdnum`, `ifnum`, `plnum`, `polariza`, beam model version など

これは provenance として非常に重要です。

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

#### 20.7 converter は現状 `CTYPE1='FREQ'` 主体

Setsumei にある VRAD / VOPT / REL の議論は後段解析の意味論として重要ですが、現 converter の主系統は周波数軸保存です。

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

## add_row() を大量に呼ぶ
## chunk_size に達すると自動で part を flush

writer.close()
```

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

本書は `sd_radio_spectral_fits` の reader / standardizer / coadd が、入力 SDFITS に対して**何を前提にしているか**を整理した仕様書です。単なる I/O 仕様ではなく、

- 読み込み時にどの列をどう解決するか
- 時刻をどの順で解釈するか
- 速度軸 / 周波数軸 / RESTFREQ / SPECSYS をどう扱うか
- 可変長配列 VLA をどう扱うか
- coadd がどの row を同じグループとして束ねるか
- 出力側へどの列を再設定して返すか

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

---

### 3. 時刻解決の契約

#### 3.1 優先順位

reader 系ユーティリティは、row 時刻を次の順に解決します。

1. `TIMESTAMP`
2. `MJD`
3. `DATE-OBS` / `DATEOBS` + `TIME`
4. 旧形式 `TIME`

これは `fitsio._resolve_table_timestamps()`、`scantable_utils._resolve_table_timestamps()`、`rawspec._resolve_mapping_timestamps()` などに共通した設計です。

#### 3.2 各列の意味

| 列 | reader が期待する意味 |
|---|---|
| `TIMESTAMP` | 各 row の絶対 UTC 時刻。最優先 |
| `MJD` | 同じ時刻を MJD UTC days で表したもの |
| `DATE-OBS` / `DATEOBS` | 絶対 UTC 時刻文字列 |
| `TIME` | 標準意味では `DATE-OBS` からの秒オフセット |
| 旧 `TIME` | `TIMESTAMP/MJD/DATE-OBS` が無い場合のみ、MJD days または Unix 秒として互換解釈 |

#### 3.3 現 writer との整合

現 writer は

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

`TOPOCENT` 入力から `VELOSYS` / `VFRAME` を計算する場合、coadd は各 row の sky position を必要とします。通常は

- `RA`,`DEC`
- または frame に応じた対応列

が必要です。

#### 4.2 `RA/DEC` と `GLON/GLAT`

writer は `RA/DEC` だけ、または `GLON/GLAT` だけからでも相互変換して両方を埋めます。reader 側もこの両方があると後段処理が安定します。

#### 4.3 `RADESYS` / `EQUINOX`

reader / Doppler 補正では、ヘッダまたは table にある座標系情報を使います。RA/DEC の意味論を曖昧にしないため、`RADESYS` と `EQUINOX` は dataset-level に持っておくのが望ましいです。

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

#### 5.2 `CTYPE1='FREQ'` の場合

最も単純なケースです。`RESTFREQ` がなくても単なる周波数軸としては読めます。ただし

- 速度窓指定
- Doppler 補正
- coadd の velocity-grid 化

をするなら `RESTFREQ` が必要です。

#### 5.3 `CTYPE1='VRAD'` など velocity 軸の場合

この場合 reader / baseline / coadd は `RESTFREQ` を必須と考えます。`RESTFREQ` が無いと、速度窓から周波数へ戻せず失敗します。

#### 5.4 `CUNIT1` は混在不可

writer は 1 つのテーブル内で `CUNIT1` が混在すると書き出しを拒否します。reader / coadd 側でも、混在した入力は結局どこかで問題になります。1ファイル1単位を前提とするのが安全です。

---

### 6. `RESTFREQ` / `VELDEF` / `SPECSYS` の契約

#### 6.1 `RESTFREQ`

次のどれかを行うなら必須です。

- line 観測の解釈
- velocity 軸への変換
- baseline / rms window を速度指定する処理
- LSRK / BARY / HELIO などの frame 解釈
- coadd の velocity-grid 化

#### 6.2 `VELDEF`

writer は velocity 文脈で `VELDEF` を要求します。reader / coadd もこれを期待します。`RADIO` 系が前提の実装が多いため、特に `LSRK` を使うなら `RADI-LSR` 相当が安全です。

#### 6.3 `SPECSYS`

coadd は `SPECSYS` を強く見ます。

- `TOPOCENT` を含むなら未補正データとして扱い、必要なら `VELOSYS` を計算する
- `LSRK` / `HELIO` / `BARY` 等なら、原則として既に補正済みとみなす

#### 6.4 `VELOSYS` / `VFRAME`

`TOPOCENT` 入力で velocity coadd を行う場合、

- 既に `VELOSYS` または `VFRAME` が入っていればそれを使う
- 無ければ時刻と座標から計算する
- どちらも無理ならエラー

というのが coadd の契約です。

---

### 7. Standardizer が前提とすること

#### 7.1 役割

`regrid_vlsrk.Standardizer` は、不均質な Scantable を共通速度グリッドへ再配置するクラスです。

#### 7.2 AxisSignature

Standardizer は row を次の軸署名でグループ化します。

- `NCHAN`
- `CRVAL1`
- `CDELT1`
- `CRPIX1`
- `CTYPE1`
- `CUNIT1`
- `RESTFREQ`

つまり、**軸が異なる row は別グループとして扱う**のが基本です。

#### 7.3 row 長の扱い

`DATA` が VLA でも、各 row の長さをそのまま見て軸署名を作ります。したがって VLA は Standardizer と矛盾しません。

#### 7.4 `RESTFREQ` が必要な場面

Standardizer の内部では、観測周波数軸から速度軸を作るとき `RESTFREQ` を必要とします。無い場合、その row グループは正しく扱えません。

#### 7.5 `SPECSYS='TOPOCENT'` とそれ以外

Standardizer は

- `TOPOCENT` なら `VELOSYS/VFRAME` を使って row ごとの補正を掛ける
- それ以外は既に補正済みとみなして 0 扱い

というルールを持っています。

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

#### 8.3 `scan` grouping の契約

`group_mode='scan'` のとき、coadd は単に `SCAN` だけで束ねません。実際には

- `SCAN`
- `OBSMODE`
- `FDNUM`
- `IFNUM`
- `PLNUM`
- さらに `_INPUT_ID`

で group を作ります。

これは、異なる beam / IF / polarization / 入力ファイルが同じ `SCAN` 番号を持っていても、誤って混ざらないようにするためです。

#### 8.4 `FDNUM/IFNUM/PLNUM` が無い場合

後方互換として、これらが無い入力は 0 扱いされます。しかし multi-beam / multi-pol 系では明示しておくべきです。

#### 8.5 出力側で上書きされる列

velocity coadd の出力では、代表的に次が再定義されます。

- `CTYPE1`
- `CUNIT1`
- `CRVAL1`
- `CDELT1`
- `CRPIX1`
- `SPECSYS`
- `SSYSOBS`
- `RESTFRQ`
- `RESTFREQ`
- `VELDEF`
- `VELOSYS` / `VFRAME`（観測時値は `*_OBS` へ退避して 0 にする）

つまり coadd 出力は、**観測生データの記録**というより、**共通グリッドへ再投影された解析生成物**です。

---

### 9. calibration / baseline が前提とすること

#### 9.1 calibration

`calibrate.py` は

- `TIMESTAMP` / `MJD` / `DATE-OBS` / legacy `TIME`
- `RESTFREQ` / `RESTFRQ`
- WCS 基本列

を見てキャリブレーションします。周波数軸ベースの処理では `RESTFREQ` が無いと失敗します。

#### 9.2 baseline

`baseline.py` は velocity window を使うとき、`RESTFRQ/RESTFREQ` と WCS を必要とします。したがって baseline を速度窓指定で使うなら、writer 側で `RESTFREQ` を省略してはいけません。

---

### 10. reader / standardizer / coadd から見た推奨 SDFITS

#### 10.1 最低限これを満たすと安全

| 項目 | 推奨 |
|---|---|
| 時刻 | `TIMESTAMP` と `MJD` を両方持つ |
| WCS | `CRVAL1`,`CDELT1`,`CRPIX1`,`CTYPE1`,`CUNIT1` を row に持つ |
| line 情報 | `RESTFREQ` と `RESTFRQ` を両方持つ |
| frame | `SPECSYS` と `SSYSOBS` を揃える |
| velocity 補正 | `TOPOCENT` なら `VELOSYS` か有効時刻+座標を持つ |
| beam/pol | `FDNUM`,`IFNUM`,`PLNUM`,`POLARIZA` を持つ |
| 装置 provenance | `BACKEND`,`SAMPLER` を可能なら持つ |
| 配列 | `DATA` と `FLAG` の長さを必ず一致させる |

#### 10.2 LSRK 入力で特に重要な点

writer 側で `SPECSYS='LSRK'` を使う場合、`store_freq_column=True` を使って `FREQ` を実際の row 周波数として持たせる設計が安全です。現 writer もこの条件を事実上要求しています。

#### 10.3 `TEMPSCAL` / `BEAMEFF`

coadd は `TEMPSCAL` と `BEAMEFF` を整理しながら進みます。`TR*` へ変換する場合、`BEAMEFF` がどの row / beam / pol に適用されるかが明確である必要があります。

---

### 11. 失敗しやすいパターン

#### 11.1 `TOPOCENT` なのに時刻が不完全

- `TIMESTAMP` なし
- `MJD` なし
- `DATE-OBS` も信用できない
- `VELOSYS` / `VFRAME` もない

この場合、coadd は速度補正を計算できません。

#### 11.2 `RESTFREQ` が無いのに velocity 窓を使う

baseline / calibrate / coadd のいずれかで失敗します。

#### 11.3 multi-beam なのに `FDNUM` が無い

scan grouping や `set_beameff()` 適用で曖昧さが残ります。

#### 11.4 `POLARIZA` を書かず `PLNUM` だけで済ませる

後段解析で物理偏波の意味が失われます。

#### 11.5 row ごとに `CUNIT1` を混在させる

writer で弾かれるか、reader / coadd で扱いにくくなります。

---

### 12. converter を含めた推奨運用

`necst_v4_sdfits_converter.py` は、この契約にかなり忠実に従っています。特に良い点は次の通りです。

- stream ごとに `fdnum/ifnum/plnum/polariza/backend/sampler/frontend` を明示
- `MJD` / `DATE-OBS` / `TIMESTAMP` を揃える
- beam-center Az/El、boresight、command、correction を分けて保存
- `GLON/GLAT` を converter 側で計算して保存
- heterodyne 文脈では `OBSFREQ` / `SIDEBAND` の整合を取る
- 列の意味論を `HISTORY` に詳しく残す

したがって、新しい writer / converter を作る場合も、**reader / standardizer / coadd が困らない列を最初から揃えて出す**ことが最も重要です。

---

### 13. 最終結論

このパッケージ全体の観点から見ると、良い SDFITS とは次の条件を満たすものです。

1. `SINGLE DISH` の各 row が自己完結した metadata を持つ。  
2. `DATA` は fixed-length でも VLA でもよいが、WCS と時刻が解決できる。  
3. `RESTFREQ`、`SPECSYS`、`TIMESTAMP/MJD` が明確である。  
4. `FDNUM`、`IFNUM`、`PLNUM`、`POLARIZA` により beam / IF / pol が曖昧でない。  
5. `TOPOCENT` データでは `VELOSYS` を計算できるだけの時刻・座標情報がある。  
6. 解析生成物として coadd が新しい WCS を安全に上書きできる。  

この契約を守る限り、reader、standardizer、coadd、baseline、calibration は相互に整合して動作しやすくなります。

---

## 付記

本書は Markdown 版です。今後必要であれば、この統合版をもとに PDF 版または章立てをさらに整理した配布版へ展開できます。
