# 単一鏡電波分光解析パッケージ 説明書
## 解析パイプライン・I/O・基準系・温度スケール・coadd・Standardizer 総合ガイド

> 対象  
> 添付された最新版アーカイブに含まれる実装に基づく総合説明書です。主に `fitsio.py`, `rawspec.py`, `atmosphere.py`, `calibrate.py`, `baseline.py`, `regrid_vlsrk.py`, `coadd.py`, `restfreq.py`, `tempscale.py`, `scantable_utils.py`, `axis.py`, `doppler.py`, `ranges.py`, `sdfits_writer.py`, `sdfits_bintable.py` を対象とします。  
> 旧版説明書の情報量は落とさず、最新版で追加・整理された意味論を反映しています。

---

## 1. このパッケージは何をするものか

このパッケージは、単一鏡の分光データを

1. 生データ読み込み  
2. HOT / OFF / ON からの $T_A^*$ キャリブレーション  
3. 必要に応じた静止周波数の上書き  
4. ベースライン推定・減算  
5. 行を保ったままの LSRK 共通速度軸再配置  
6. scan 単位・位置単位・INTGRP 単位の coadd  
7. SDFITS への保存  

まで、一貫したデータ型で扱うための解析基盤です。

設計上の核は `Scantable` で、これは

- `meta`: ファイル全体に共通するメタデータ
- `data`: 各行スペクトル
- `table`: 行ごとの補助情報
- `history`: 処理履歴

をまとめて保持します。

この設計により、スペクトル配列だけでなく、各行がどの基準系・どの静止周波数・どの温度スケールで解釈されるかを、常に同じコンテナの中で追跡できます。

---

## 2. まず押さえるべき全体像

### 2.1 解析の標準的な順番

標準的な流れは次です。

1. `fitsio.read_scantable()` または `rawspec.load_rawspec_auto()` で入力を読む  
2. `calibrate.run_tastar_calibration()` で $T_A^*$ を作る  
3. 必要なら `restfreq.apply_restfreq_override()` で静止周波数を上書きする  
4. `baseline.run_baseline_fit()` でベースラインを評価・減算する  
5. `regrid_vlsrk.run_velocity_regrid()` で各行を共通 LSRK 速度軸へそろえる  
6. `coadd.run_velocity_coadd()` で必要な単位で積分する  
7. `fitsio.write_scantable()` または `fitsio.write_sdfits()` で保存する

### 2.2 典型的な思想

このパッケージは、次の分離を非常に重視します。

- 強度較正は `calibrate.py` の仕事
- ベースライン評価は `baseline.py` の仕事
- 共通速度軸への変換は `regrid_vlsrk.py` と `coadd.py` の仕事
- 温度スケール変換は `tempscale.py` と書き出し時の明示的操作の仕事

重要なのは、読み込み時に勝手に温度スケールを変えないことです。2026 方針では、ファイルに `TEMPSCAL=TR*` と書かれていれば、そのまま `TR*` として保持します。内部配列を自動で $T_A^*$ に戻すことはしません。

同様に、静止周波数の上書きも軸の種類によって扱いを分けます。

- `CTYPE1=FREQ` なら周波数 WCS は保ち、`RESTFRQ` だけを更新する
- `CTYPE1=VRAD` なら速度軸 WCS 自体を厳密に更新する

---

## 3. このパッケージの中核データ構造

## 3.1 Scantable

`fitsio.Scantable` は、解析の主コンテナです。

### 役割

- `meta`: PRIMARY HDU 相当の全体メタデータ
- `data`: スペクトル本体
- `table`: 各行に付随する metadata
- `history`: 処理履歴

### 特徴

- `data` は通常の 2 次元配列でもよい
- 行ごとに長さが異なる場合は `List[np.ndarray]` を使える
- `copy()` は `meta`, `table`, `history` をコピーし、`data` も安全にコピーする
- `write_scantable()` へそのまま渡せる

### Scantable を使うと何が良いか

行ごとの `RESTFRQ`, `CRVAL1`, `CDELT1`, `SPECSYS`, `TEMPSCAL`, `BEAMEFF`, `VELOSYS` などを、スペクトル配列と分離せず保持できます。これにより、

- 同じファイルに複数の IF / 偏波 / ビームが混在しても扱いやすい
- 再グリッドや coadd のときに、各行の意味を失いにくい
- SDFITS の VLA も自然に表現できる

## 3.2 RawSpec

`rawspec.RawSpec` は、HOT / OFF / ON の生データを扱うためのコンテナです。

### 構成

- `meta`
- `hot`
- `off`
- `on`
- `mapping`

### 想定用途

- 観測直後の HOT / OFF / ON をまとめて保持する
- まだ較正されていない生データを、安全に Ta* 化する
- `mapping` に HOT / OFF / ON を含む行 metadata を持たせる

### 重要

最終的に `run_tastar_calibration()` の出力は ON 行主体の `Scantable` になりますが、補間や gain 計算には HOT / OFF の時系列情報が必要です。そのため RawSpec では `mapping` に全観測行を持つ構造が重要です。

---

## 4. 重要仕様 1: スペクトル軸・基準系・速度補正

## 4.1 基本原則

### 軸の物理量

主に次を扱います。

- `CTYPE1 = FREQ`: 周波数軸
- `CTYPE1 = VRAD`: radio definition の速度軸

### 軸の単位

最新版実装では、単位も厳密に見ます。

- 周波数軸は `Hz` を標準とする
- 速度軸は `m/s` または `km/s` を解釈できる
- `Standardizer` は `CUNIT1` を見て、周波数 WCS の単位変換を明示的に行う

### 基準系

- `SPECSYS`: 現在その軸が属している基準系
- `SSYSOBS`: 観測時に固定だった基準系。`regrid` / `coadd` / `baseline` 後も、通常はこの意味で保持される

Raw データが TOPOCENT で保存されているときは、未適用の観測者速度補正を行ごとに持っている可能性があります。逆に LSRK で保存されているときは、その補正は既に消化済みであるべきです。

### 静止周波数

静止周波数は `RESTFRQ` を正本とし、互換用に `RESTFREQ` も同じ値で持ちます。最新版では両者を可能な限り同期させる方針です。

### 観測者速度補正

内部では、行ごとの未適用補正を `VELOSYS` または `VFRAME` で扱います。意味は

- `VELOSYS`, `VFRAME` は原則として $\mathrm{m/s}$
- 旧 `V_CORR_KMS` は $\mathrm{km/s}$ であり、読込時に自動移行される

です。

ドップラー係数は

$$
k = \sqrt{\frac{1+\beta}{1-\beta}}, \qquad \beta = \frac{v}{c}
$$

で与えられ、周波数変換は

$$
\nu_{\rm LSRK} = \frac{\nu_{\rm obs}}{k}
$$

として扱われます。

## 4.2 フェーズごとの意味

### Raw

ドップラー補正をまだ消化していないなら、通常は

- `CTYPE1=FREQ`
- `SPECSYS=TOPOCENT`
- `SSYSOBS=TOPOCENT`

を想定します。

### Calibration

`run_tastar_calibration()` は強度較正が主目的であり、原則として軸の意味そのものは壊しません。したがって

- 周波数軸のまま Ta* を作る
- 必要なら `VFRAME` などの速度補正列を維持する
- `TEMPSCAL` は `TA*` を明示する

という動作になります。

### Coadd

`run_velocity_coadd()` は、入力が TOPOCENT なら行ごとの未適用補正を使って LSRK 共通軸へ揃えます。入力が既に LSRK なら、未適用補正が残っていないことを確認してから処理します。

## 4.3 実務上の注意

### TOPOCENT 入力

TOPOCENT 入力では、各行について LSRK へ移るための補正が必要です。最新版では、次のどちらかが必要です。

1. `VELOSYS` / `VFRAME` が既に行ごとに入っている
2. 有効な時刻列と座標列とサイト情報があり、`compute_vcorr_series()` で再計算できる

これがないと、coadd や regrid は停止します。

### LSRK 入力

LSRK 入力では、未適用補正列が非ゼロで残っていることを不整合として扱います。つまり

- `SPECSYS=LSRK`
- なのに `VELOSYS` や `VFRAME` が非ゼロ

という状態はエラーです。

---

## 5. 重要仕様 2: rest frequency の上書き

## 5.1 FREQ 軸のとき

`CTYPE1=FREQ` の場合、静止周波数を別の分子線へ差し替えても、観測された周波数自体は変わりません。したがって最新版の契約は

- `CRVAL1`, `CDELT1`, `CRPIX1` は変えない
- `RESTFRQ`, `RESTFREQ` を上書きする

です。

これは、同じ raw から $^{13}\mathrm{CO}$ と $\mathrm{C}^{18}\mathrm{O}$ を別々に解釈するような用途で重要です。

## 5.2 VRAD 軸のとき

`CTYPE1=VRAD` では、静止周波数を変えると速度軸そのものの意味が変わります。そこで `apply_restfreq_override()` は厳密なアフィン変換を使います。

静止周波数を $\nu_1$ から $\nu_2$ へ変えるとき、

$$
a = \frac{\nu_1}{\nu_2}
$$

$$
b = c\left(1 - \frac{\nu_1}{\nu_2}\right)
$$

$$
v_2 = a v_1 + b
$$

です。したがって

$$
\mathrm{CRVAL1}_2 = a\,\mathrm{CRVAL1}_1 + b
$$

$$
\mathrm{CDELT1}_2 = a\,\mathrm{CDELT1}_1
$$

で、`CRPIX1` は変えません。

### まとめ

- FREQ 軸: 周波数 WCS は不変、静止周波数だけ更新
- VRAD 軸: 速度 WCS も更新

---

## 6. 重要仕様 3: TEMPSCAL と BEAMEFF

## 6.1 基本方針

2026 方針では、読み込み時に温度スケールを勝手に変換しません。

### つまり

- `TEMPSCAL=TR*` なら、そのまま `TR*` として読む
- `TEMPSCAL=TA*` なら、そのまま `TA*` として読む
- 内部配列を自動で別スケールへ変えない

## 6.2 BEAMEFF の意味

主ビーム効率を $\eta_{\rm mb}$ とすると、

$$
T_R^* = \frac{T_A^*}{\eta_{\rm mb}}
$$

$$
T_A^* = T_R^*\,\eta_{\rm mb}
$$

です。

`tempscale.py` の `ta_to_tr()` と `tr_to_ta()` はこの式をそのまま実装しています。VLA では `convert_rowwise_vla()` が行ごとに同じ変換を適用します。

## 6.3 混在グループへの対処

coadd のグループ内で `BEAMEFF` が混在すると、単純平均は危険です。最新版では

- `normalize_if_mixed="auto"`: 混在グループを内部的に $T_R^*$ へ正規化してからまとめる
- `normalize_if_mixed="never"`: 正規化せず進めるが危険モード

です。

`auto` の場合、グループ内で一度

$$
T_{R,i}^* = \frac{T_{A,i}^*}{\eta_{{\rm mb},i}}
$$

へ変換し、coadd 後に必要な出力スケールへ戻します。

---

## 7. 重要仕様 4: big-endian と pandas の問題

FITS 由来の配列や DataFrame は big-endian で来ることがあります。これをそのまま pandas に渡すと、`iloc`, `query`, `isin`, ソートなどで問題になる場合があります。

最新版では `read_scantable()`, `find_scans()`, `filter_scantable()`, `merge_scantables()` などで native-endian 正規化を行います。

### 実務上の意味

- FITS 読み込み後に DataFrame をそのまま扱いやすい
- ただし外部で別経路から作った DataFrame は、自前で endian をそろえる方が安全

---

# 8. モジュール別の詳細説明

## 8.1 `fitsio.py`

### 主な役割

- `Scantable` の定義
- SDFITS の読込・書込
- 旧形式 DATA / DUMPS と SINGLE DISH の両対応
- 温度スケール変換の書き出し時適用

### 公開 API

- `Scantable`
- `read_scantable(path, tr_input_policy="preserve")`
- `write_scantable(path, scantable, spectrum_column="DATA", overwrite=True, **kwargs)`
- `read_tastar_fits(path)`
- `write_sdfits(out_path, meta, data, table, history=None, ...)`

### `read_scantable()`

#### 何をするか

`read_tastar_fits()` を呼んだあと、最新版の意味論に合わせて列と履歴を正規化し、`Scantable` を返します。

#### 重要

1. legacy 列の自動移行  
   - `V_CORR_KMS` があれば $\mathrm{km/s}$ から $\mathrm{m/s}$ へ直して `VFRAME`, `VELOSYS` へ移す  
   - `TAMB_K` は `THOT` へ rename する  
   - `TAU` は `TAU0` へ rename する

2. `VFRAME` だけあって `VELOSYS` が無い場合は mirror する

3. `TEMPSCAL`, `BEAMEFF` 列をそろえるが、TR* を TA* へ自動変換しない

4. `tr_input_policy` という引数は公開シグネチャに残っていますが、最新版実装の挙動は実質 `preserve` 固定です。つまり、読込時に尺度変換はしません。

### `write_scantable()`

#### 何をするか

`Scantable` を壊さず、書き出しに必要な列をヘッダから table へ昇格させて `write_sdfits()` へ渡します。

#### 主な自動昇格列

- `RESTFRQ`, `RESTFREQ`
- `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`
- `SPECSYS`, `SSYSOBS`, `VELDEF`
- `VELOSYS`, `VFRAME`
- `TEMPSCAL`, `BEAMEFF`

つまり、ヘッダにしか無かった情報も、行ごとに解釈できる形へ寄せてから保存します。

### `write_sdfits()`

#### 何をするか

これが推奨の低レベル保存 API です。`SINGLE DISH` BinTable を構築し、必要なら `HISTORY` 拡張も付けます。

#### 主要引数

- `out_path`: 出力先
- `meta`: PRIMARY HDU 相当のヘッダ辞書
- `data`: 2 次元配列または VLA 用リスト
- `table`: 行 metadata
- `history`: 履歴辞書
- `spectrum_column`: スペクトル列名
- `include_flag`: `FLAG` 列を含めるか
- `string_widths`: 文字列列の幅指定

#### `kwargs` で実質使うもの

実装上、特に重要なのは次です。

- `tempscal` または `out_scale`: 保存時の on-disk スケール
- `data_scale`: 今メモリ上に載っている配列のスケール
- `beameff_on_fail`: 変換に必要な `BEAMEFF` が無いときの扱い

#### 重要

保存時スケール変換は非破壊です。たとえばメモリ中が $T_A^*$ で、保存だけ $T_R^*$ にしたい場合、書き出し時だけ

$$
T_R^* = \frac{T_A^*}{\eta_{\rm mb}}
$$

を適用して保存します。元の `Scantable` は変えません。

また、`TIMESTAMP`, `MJD`, `DATE-OBS`, `DATEOBS`, `TIME` は、手元の列や index から可能な限り補完されます。最新版では各行の `DATE-OBS` をフル UTC 時刻として持つ思想なので、`TIME` は各行 0 秒で埋まる場合があります。

## 8.2 `rawspec.py`

### 公開 API

- `RawSpec`
- `build_rawspec()`
- `save_rawspec()`
- `load_rawspec()`
- `load_rawspec_auto()`
- `load_rawspec_fits()`

### `build_rawspec()`

#### 何をするか

HOT / OFF / ON / meta / mapping を一つの RawSpec にまとめます。

#### 引数の意味

- `hot`, `off`, `on`: それぞれの時系列スペクトル
- `meta`: 全体ヘッダ
- `mapping`: 全観測行の表

#### 注意

`mapping` は ON だけでなく HOT / OFF も含む前提です。時系列補間や gain の作成で必要になります。

### `load_rawspec_auto()`

pickle 形式か FITS 形式かを見て、適切な loader を選びます。

### `load_rawspec_fits()`

SDFITS から HOT / OFF / ON を読み分けて RawSpec を構築します。

## 8.3 `atmosphere.py`

### 公開 API

- `get_airmass()`
- `estimate_t_atm()`
- `compute_t_cal_array()`
- `extract_meta_value()`
- `extract_meta_array()`

### 何をしているか

大気モデル込みの chopper-wheel 補正で必要な、

- エアマス
- 大気有効温度
- 等価較正温度 $T_{\rm cal}$

を計算します。

### 重要仕様

エアマスは

$$
X = \sec z = \frac{1}{\sin({\rm El})}
$$

ですが、実装では低仰角で暴れないよう、$\sin 5^\circ$ で下駄を履かせます。

大気有効温度は二つの簡易モデルを持ちます。

- `offset`: $T_{\rm atm} = T_{\rm surf} - \Delta T$
- `ratio`: $T_{\rm atm} = \eta T_{\rm surf}$

等価較正温度は実装上

$$
T_{\rm cal} = T_{\rm hot} e^{\tau X} - T_{\rm atm}(e^{\tau X}-1) - T_{\rm bg}
$$

です。

`extract_meta_value()` と `extract_meta_array()` は、`TAU0` / `TAU` / `OPACITY` のような alias を許して、table 優先で値を取る補助関数です。

## 8.4 `calibrate.py`

### 公開 API

- `make_tastar_dumps()`
- `tastar_from_rawspec()`
- `run_tastar_calibration()`
- `recalibrate_tastar()`

### `make_tastar_dumps()`

#### 何をするか

RawSpec から ON 行の $T_A^*$ を作ります。行選択、チャンネル切り出し、速度窓からのチャンネル選択、大気モデル、静止周波数の上書き、速度補正列の扱いまで含みます。

#### 入力

- `raw`: RawSpec 互換オブジェクト
- `rows`, `exclude_rows`: 行選択
- `t_hot_k`: ホットロード温度
- `tau_zenith`: 光学的厚み。数値または `"auto"`
- `t_surface_k`: 地表気温
- `t_atm_model`, `t_atm_delta_k`, `t_atm_eta`: ATM モデル
- `gain_mode`: `"hybrid"` または `"independent"`
- `ch_range` または `vlsrk_range_kms`: 使用チャンネル範囲
- `v_corr_col`: 速度補正列の優先名
- `vcorr_chunk_sec`: 速度補正再計算の間引き幅
- `dtype`: 保存 dtype
- `rest_freq`: 静止周波数上書き

#### 主な内部処理

ON 行の時刻を $t_{\rm on}$、補間された OFF を $P_{\rm off}(t_{\rm on})$、gain 計算用分母を $D(t_{\rm on})$ とすると、基本式は

$$
T_A^*(t_{\rm on}) = \left(P_{\rm on}(t_{\rm on}) - P_{\rm off}(t_{\rm on})\right)
\frac{T_{\rm cal}(t_{\rm on})}{D(t_{\rm on})}
$$

です。

`independent` では

$$
D(t_{\rm on}) = P_{\rm hot}(t_{\rm on}) - P_{\rm off}(t_{\rm on})
$$

を ON 時刻で別々に補間します。

`hybrid` では、まず HOT の時刻で OFF を補間して

$$
D(t_{\rm hot}) = P_{\rm hot}(t_{\rm hot}) - P_{\rm off}(t_{\rm hot})
$$

を作り、それを ON 時刻へ再補間します。これにより、HOT 側の分母に乗る定在波や時間変動の影響を少し抑えやすくなります。

#### 全引数の説明

旧説明書で重要だった引数はそのまま有効で、最新版ではさらに次が重要です。

- `rows`, `exclude_rows`: ON/OFF/HOT 全体の入力表に対する選択
- `rest_freq`: FREQ 軸なら静止周波数のみ差し替え、VRAD 軸なら WCS も更新
- `v_corr_col`: 旧 `V_CORR_KMS` を含む候補列の優先順位先頭

#### `gain_mode` の意味

##### `hybrid`（推奨）

HOT 時刻で作ったクリーンな分母を ON 時刻へ補間する方式です。

##### `independent`

従来型で、HOT と OFF をそれぞれ ON 時刻へ補間して差を作ります。

#### ATM 発動条件

`tau_zenith` が `None` でなければ大気モデルを使います。`"auto"` のときは `TAU0` などを table / meta から探します。

#### 出力 table に付く重要列

- `THOT`
- `TAMB_K`（後方互換のため残す）
- `TCAL`
- `CALSTAT`
- `TAU0`（使った場合）
- `CH_START`, `CH_STOP`, `NCHAN_SEL`
- `TEMPSCAL="TA*"`

### `run_tastar_calibration()`

#### 何をするか

入力がファイルでも RawSpec でも受け取り、必要なら保存まで行う高レベル API です。

#### 追加引数

- `output_path`: 保存先
- `spectrum_column`: 保存時の列名
- `store_freq_column`: `FREQ` ベクトル列を持つか
- `overwrite`: 既存ファイル上書き

#### 注意

`rest_freq` の上書きは較正の前段で反映されるため、以後の速度解釈はその静止周波数に従います。

### `recalibrate_tastar()`

#### 何をするか

既に一度作った Ta* に対して、`TAU0` や大気温度モデルだけを変えて再評価します。

#### 原理

保存されている `TCAL` 相当の情報と履歴を用い、強度スケールを再構成します。

#### 引数

- `new_tau`
- `new_t_surface_k`
- `new_t_atm_model`
- `new_t_atm_delta_k`
- `new_t_atm_eta`

#### 挙動

入力 `Scantable` を壊さず、新しい `Scantable` を返します。

#### 注意

元データに必要な履歴や大気情報が無い場合、完全再現はできません。

## 8.5 `baseline.py`

### 公開 API

- `BaselineModel`
- `fit_polynomial_baseline()`
- `run_baseline_fit()`

### `fit_polynomial_baseline()`

#### 何をするか

与えられた軸 `x_axis` とスペクトル `y` に対して、多項式ベースラインをフィットします。

#### 引数

- `base_windows` または `windows`: ベースライン窓
- `line_windows`: 除外したい線窓
- `poly_order`
- `iter_max`, `iter_sigma`: 反復クリップ

#### 返り値

- `coeffs`: 多項式係数
- `baseline`: 評価されたベースライン
- `info`: 追加情報

#### `info` の中身

少なくとも `std` と `mask` が重要です。`std` は最終 line-free 領域の標準偏差で、`BSL_RMS` へ保存されます。

#### 実装上の特徴

line-free 窓は速度窓指定で与えられますが、実際の軸は各行の `SPECSYS`, `RESTFRQ`, `VELOSYS` から構成されます。したがって、baseline 窓は常に「物理軸」で解釈されます。

### `run_baseline_fit()`

#### 何をするか

`Scantable` 全体に対して、行ごとに baseline を評価し、必要なら減算します。VLA にも対応します。

#### 全引数の説明

- `rows`, `exclude_rows`: 行選択
- `vwin`: baseline 窓。必須
- `poly_order`: 多項式次数
- `line_vwin`: 追加の線窓除外
- `iter_max`, `iter_sigma`: 反復クリップ
- `max_dumps`: 先頭何行まで使うか
- `ch_start`, `ch_stop`: 先に切り出すチャネル範囲
- `v_corr_col`: 速度補正列。既定は `VELOSYS`
- `rest_freq`: 静止周波数上書き
- `apply`: 実際に差し引くか。`False` なら評価だけ
- `bsl_overwrite`: 既存 `BSL_*` を置き換えるか
- `on_fail`: `exit` / `warn` / `skip`

#### 重要挙動

- baseline は行ごとの WCS / frame metadata を解釈して動作し、`SSYSOBS` を `SPECSYS` で上書きしません。したがって TOPOCENT 観測を LSRK へ再配置した後でも、通常は `SPECSYS='LSRK'`, `SSYSOBS='TOPOCENT'` の関係を保持します。

1. TOPOCENT 行では、速度窓解釈のために行ごとの速度補正が必要です。  
2. LSRK 行では、未適用補正が非ゼロならエラーです。  
3. `apply=False` でも `BSL_*` は評価結果として保存されます。  
4. channel slice をした場合、WCS もその切り出しに合わせて更新します。

#### 出力で増える主な列

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

`BSL_COEF` は多次元列になり得るため、SDFITS 保存時には vector-in-cell として処理されます。

#### 温度スケール方針

最新版では、baseline fitting 自体は入力スケールのまま行います。`TR*` を見たからといって内部で $T_A^*$ に戻すことはしません。どのスケールで評価したかは `BSL_SCALE` へ記録されます。

## 8.6 `regrid_vlsrk.py`

### 公開 API

- `VGrid`
- `make_vgrid()`
- `get_axis_signature()`
- `Standardizer`
- `vlsrk_axis_for_spectrum()`
- `interp_to_vgrid()`
- `vrange_from_meta_and_vcorr()`
- `run_velocity_regrid()`

### `VGrid`

共通速度グリッドを表す単純なコンテナです。`v0_kms`, `dv_kms`, `nchan`, `crpix1` を持ちます。

### `make_vgrid()`

#### 何をするか

$v_{\min}$, $v_{\max}$, $\Delta v$ から共通グリッドを作ります。

#### 注意

`dv_kms` は正でなければなりません。

### `get_axis_signature()`

#### 何をするか

各行のスペクトル軸を、

- `nchan`
- `CRVAL1`
- `CDELT1`
- `CRPIX1`
- `CTYPE1`
- `CUNIT1`
- `RESTFRQ`
- `SPECSYS`

の組として表し、同じ軸を持つ行をまとめるために使います。

#### なぜ必要か

軸計算を行ごとに毎回やると重いので、同じ signature を持つ行はまとめて処理します。

### `Standardizer`

#### 何をするクラスか

不均質な `Scantable` を、共通 LSRK 速度グリッドを持つ行列へ変換します。

#### 内部の重要思想

1. 行を signature ごとにグループ化する  
2. 元の観測周波数軸をできるだけ正確に復元する  
3. TOPO 行なら relativistic Doppler factor で周波数基準系を直す  
4. その後に radio velocity へ変換する  
5. 目標グリッドへの補間または rebinning を行う

補間だけでなく、粗いグリッドへ落とす場合には overlap-weighted rebin を使う経路が用意されており、same/smaller `dv` と clearly larger `dv` を分けて扱います。

### `Standardizer.auto_determine_grid()`

#### 何をするか

入力 `Scantable` を全部含むマスターグリッドを自動決定します。

#### 引数

- `dv`
- `vmin`
- `vmax`
- `margin`

`__init__` ではさらに

- `auto_grid_mode`
- `auto_grid_anchor_kms`
- `auto_grid_tol_frac`

を指定できます。

#### 決め方

既定の `stable_native` では、

- 全入力で `nchan` がそろい
- `dv` もほぼそろっている

場合、native のチャネル幅と幅をできるだけ保つ安定なグリッドを返します。従来の global min/max 包含だけだと端点の微差が 1--2 channel の差に増幅されやすいので、anchor 量子化を使います。

`legacy_envelope` は旧来型の global envelope 方式です。

### `Standardizer.get_matrix()`

#### 何をするか

`Scantable` 全体を、$(N_{\rm row}, N_{\rm tgt})$ 行列へ変換します。

#### 処理内容

- グリッドを決める
- signature ごとにまとめる
- TOPO 行は行ごとの `v_corr_col` を使って LSRK 化する
- 共通 `v_tgt` へ補間または rebining する

#### 重要

最新版の実装では、周波数軸から LSRK 速度軸を作るときに近似

$$
v_{\rm LSRK} \approx v_{\rm native} + v_{\rm corr}
$$

を使わず、まず周波数を Doppler factor で補正してから速度へ変換します。

### `vlsrk_axis_for_spectrum()`

行ごとの LSRK radio velocity 軸を返します。

### `interp_to_vgrid()`

行の 1 本のスペクトルを目標速度軸へ補間します。

### `vrange_from_meta_and_vcorr()`

WCS と行ごとの補正から、その行がカバーする LSRK 速度範囲を返します。

### `run_velocity_regrid()`

#### 何をするか

coadd はせず、各行を保ったまま共通 LSRK 速度軸へ再配置します。

#### 全引数の意味

- `rows`, `exclude_rows`: 行選択
- `vmin_kms`, `vmax_kms`, `dv_kms`: 目標速度グリッド
- `v_corr_col`: TOPO 行で使う補正列。既定 `VFRAME`
- `rest_freq`: 静止周波数上書き
- `ch_start`, `ch_stop`: 先に元チャネルを切る
- `max_dumps`: 使う最大行数
- `fill_value`: 補間外を何で埋めるか
- `keep_row_order`: 現状は `True` のみ実装
- `drop_allnan_rows`: 完全に NaN の行を落とすか
- `history_tag`: history のキー名

#### 出力仕様

出力は

- `CTYPE1="VRAD"`
- `CUNIT1="m/s"`
- `SPECSYS="LSRK"`
- `SSYSOBS` は観測時 frame を保持。たとえば TOPOCENT 観測を regrid しても、通常は `SSYSOBS="TOPOCENT"` のまま
- `VELDEF="RADIO"`

へ更新されます。さらに

- `REGRID_DONE=True`
- `REGRID_FRAME="LSRK"`

が付き、未適用補正列は落とされます。

## 8.7 `coadd.py`

### 公開 API

- `run_velocity_coadd()`

### QC プリセット

最新版には次の QC プリセットがあります。

- `robust`
- `gentle`
- `hardonly`
- `noclip`

内部パラメータは `hard`, `soft_a`, `soft_p`, `clip_k`, `clip_iter` を持ちます。

### `run_velocity_coadd()`

#### 何をするか

複数入力または単一入力を読み込み、必要なら行ごとに LSRK 共通軸へそろえ、指定したグループ単位で coadd します。

### 全引数の説明

#### 入出力・選択

- `inputs`: path または `Scantable` の列
- `output_path`
- `rows`, `exclude_rows`
- `overwrite`
- `verbose`

#### グループ化

- `group_mode`: `position`, `scan`, `intgrp`
- `pos_col`: 位置 ID 列名
- `pos_tol_arcsec`: 自動位置グルーピングの許容幅

#### 速度補正

- `v_corr_col`: 既定 `VELOSYS`
- `coord_frame`: 速度補正再計算時の座標系
- `vcorr_chunk_sec`: 長い時系列では補正計算を間引く

#### 速度グリッド

- `vmin`, `vmax`, `dv`
- `allow_outside_overlap`
- `axis_type`: `freq` または `vel`
- `rest_freq`

#### 重み付けモード

- `mode`: `uniform`, `rms`, `rms_weight`
- `baseline_vwin`, `baseline_poly`, `baseline_iter_max`, `baseline_iter_sigma`
- `rms_vwin`, `rms_poly`, `rms_bin`
- `weight_zero_policy`: `error`, `drop`, `impute_median`
- `coadd_qc`
- `line_vwin`

#### 入力量制限

- `block_size`
- `max_dumps`
- `ch_start`, `ch_stop`

#### エラーと温度スケール

- `on_fail`
- `out_scale`: 最終出力スケール
- `normalize_if_mixed`
- `beameff_tol`
- `sigma_scale`

### 非常に重要な制約

#### `baseline_vwin` と `rms_vwin` は同時指定不可

最新版ではこれを明示的に禁止しています。前者は baseline subtraction 付き重み付け、後者は subtraction なし RMS 評価です。

#### `coadd_qc` は `mode="uniform"` と併用不可

QC は統計的重み付けを本質的に含むためです。

### グループ化の意味

#### `group_mode="position"`

同じ位置の repeat をまとめます。

#### `group_mode="scan"`

scan 内の dump を積分して 1 本にします。

#### `group_mode="intgrp"`

`INTGRP` を単位にまとめます。

### 重み付けの実際

#### `mode="uniform"`

単純平均です。

#### `mode="rms"` / `"rms_weight"`

実装上は inverse-variance を使います。各入力スペクトルの empirical RMS を $\sigma_i$ とすると、重みは

$$
w_i = \frac{1}{\sigma_i^2}
$$

です。coadd 自体は

$$
T_{\rm out}(k) = \frac{\sum_i w_i T_i(k)}{\sum_i w_i}
$$

で行います。

##### `baseline_vwin` がある

まず各入力行から baseline を引き、その residual の標準偏差を重みに使います。

##### `rms_vwin` がある

baseline subtraction はせず、指定窓で RMS を測って重みに使います。

##### どちらもない

既存の `BSL_RMS` があればそれを入力重みとして利用します。無ければ uniform に近い扱いになります。

#### `rms_bin`

RMS 測定時の spectral binning 幅です。

### QC モード

`coadd_qc` を使うと hard/soft/clip を含む QC weighting を使います。詳細パラメータは `"robust;hard=...;soft=...;clip=..."` のような文字列で上書きできます。

### coadd 後の baseline

最新版では `post_baseline_mode='inherit_all'` が重要です。これは pre-coadd の baseline 設定

- 窓
- 多項式次数
- 反復回数
- 反復 sigma

を post-coadd 側へ継承するモードです。`baseline_vwin` が指定されていないと使えません。

### BEAMEFF 混在時の処理

グループ内で `BEAMEFF` が混在すると、`normalize_if_mixed='auto'` なら内部で TR* 正規化してから coadd します。これは scale history にも記録されます。

### 出力の基準系

coadd 出力は必ず LSRK です。ただし `SSYSOBS` は `SPECSYS` に潰さず、観測時 frame を保持します。したがって TOPOCENT 観測由来の行を coadd した出力では、通常 `SPECSYS='LSRK'`, `SSYSOBS='TOPOCENT'` です。さらに未適用補正の観測時値は退避列として残されます。

- `VELOSYS_OBS`
- `VFRAME_OBS`
- 必要なら `${v_corr_col}_OBS`

また、1 本の coadd 出力 row に対して `SSYSOBS` 候補が複数混ざる場合は、曖昧な metadata を作らないためエラーで停止します。これは数値計算の変更ではなく、frame metadata の安全化です。

#### `axis_type="vel"`

出力軸は

- `CTYPE1="VRAD"`
- `CUNIT1="m/s"`

です。

#### `axis_type="freq"`

LSRK 速度グリッド $v$ から、静止周波数 $\nu_0$ を用いて

$$
\nu = \nu_0\left(1 - \frac{v}{c}\right)
$$

を使い、LSRK の `FREQ` 軸として保存します。

### 出力で増える主な列

重み付け関連:

- `WGT_MODE`
- `WGT_SRC`
- `WGT_VWIN`
- `WGT_POLY`
- `WGT_BIN`
- `WGT_STAT`
- `WGT_INPUT_SUB`
- `WGT_ZERO_POLICY`

グループ情報:

- `N_IN`
- `N_USED`
- `WEIGHT_MODE`

post-baseline を行った場合はさらに `BSL_*` が付きます。

## 8.8 `restfreq.py`

### 公開 API

- `apply_restfreq_override()`

### 何をするか

静止周波数の上書きを、FREQ 軸と VRAD 軸で意味論を分けて安全に実施します。

### 引数

- `meta`
- `table`
- `rest_freq_hz`
- `require_wcs_for_vrad`

### 戻り値

辞書で返し、そこには少なくとも

- `rest1_hz`
- `rest2_hz`
- `axis_vrad`
- `wcs_updated`

が入ります。VRAD 軸のときは `a`, `b`, `vel_unit` も入ります。

## 8.9 `tempscale.py`

### 公開 API

- `normalize_tempscal()`
- `ensure_tempscal_column()`
- `ensure_beameff_column()`
- `beameff_array()`
- `tempscal_array()`
- `require_beameff()`
- `is_beameff_mixed()`
- `representative_beameff()`
- `ta_to_tr()`
- `tr_to_ta()`
- `convert_rowwise_vla()`
- `append_scale_history()`
- `set_beameff()`

### 重要関数

#### `normalize_tempscal()`

`TA`, `TASTAR`, `TR`, `TMB` などの表記揺れを `TA*` / `TR*` の正規形へ変換します。

#### `beameff_array()`

table に `BEAMEFF` があればそれを使い、無ければ meta の `BEAMEFF` を全行へ broadcast します。

#### `tempscal_array()`

table 優先、無ければ meta の `TEMPSCAL` を全行へ展開します。

#### `require_beameff()`

`TR*` 変換に必要な `BEAMEFF` が有限かつ正であることを要求します。

#### `convert_rowwise_vla()`

VLA スペクトルの行ごと変換です。`direction` は `ta_to_tr` または `tr_to_ta` です。

#### `set_beameff()`

`Scantable` に `BEAMEFF` を付与する補助関数です。行選択に対して in-place で書き込みます。

## 8.10 `scantable_utils.py`

### 公開 API

- `describe_columns()`
- `show_scantable()`
- `calc_mapping_offsets()`
- `merge_scantables()`
- `update_metadata()`
- `set_beameff()`
- `find_scans()`
- `filter_scantable()`
- `select_rows()`

### まず押さえるべき設計思想

#### table / data / meta の役割分担

- `meta` はファイル全体の既定値
- `table` は行ごとの真の値
- `data` はスペクトル本体

#### 非破壊か破壊的か

- `filter_scantable()` や `select_rows()` は新しい `Scantable` を返す
- `update_metadata()` や `set_beameff()` は in-place 変更

#### 行指定 `rows` は positional index

行選択は DataFrame index label ではなく位置インデックスです。

#### Big-endian FITS への対策

内部で native-endian 化してから DataFrame 操作します。

#### TIMESTAMP の自動解決

`TIMESTAMP`, `MJD`, `DATE-OBS`, `DATEOBS`, legacy `TIME` を順に見て時刻列を作ります。

### `describe_columns()`

#### 何をするか

列の型や欠損状況の概要を表示します。

#### 使いどころ

入力ファイルの品質確認に便利です。

#### 注意点

VLA や vector-in-cell 列は単純な dtype 表示だけでは本質が見えないことがあります。

### `show_scantable()`

#### 何をするか

`Scantable` の表を見やすく整形して表示します。

#### 入力 `inputs`

単一 `Scantable` でも複数入力でも受けられます。

#### `columns`

`default` のほか、明示的な列リストを指定できます。

#### table に無い列をどう探すか

必要に応じて meta から補って表示します。

#### `rows`

位置インデックスで選びます。

#### `extra_data`

追加計算した列を一緒に表示できます。

#### `ref_coord`, `frame`, `projection`, `unit`

オフセット表示や座標変換のための指定です。

#### 重要な実装上の注意

DataFrame の index と `TIMESTAMP` の関係を保ったまま表示しようとします。

#### 実用例

観測ごとの `SCAN`, `OBJECT`, `TEMPSCAL`, `RESTFRQ` をざっと確認するときに便利です。

#### オフセットをその場で見たい場合

`calc_mapping_offsets()` と組み合わせます。

#### 外部計算列を足して見る

`extra_data` で任意列を差し込めます。

### `calc_mapping_offsets()`

#### 何をするか

天球上の絶対座標から、基準位置に対するオフセット `OFS_LON`, `OFS_LAT` を作ります。

#### 入力座標として認識する列

主に `RA`, `DEC`, `GLON`, `GLAT`, または meta の `OBSRA`, `OBSDEC` です。

#### `frame`

変換先フレームです。`icrs`, `fk5`, `fk4`, `galactic`, `altaz` などを受けます。

#### `ref_coord`

省略時はヘッダまたはデータ平均から決めます。

#### `projection`

- `GLS`, `SFL`: $x$ に $\cos$ 項を掛ける
- `CAR`, `NONE`: そのまま差分

#### 重要: ドキュメント文字列との差異

実装は `SIN` をサポートすると言いながら、分岐では `GLS`, `SFL`, `CAR`, `NONE` のみです。最新版の実際の仕様に従うなら、`SIN` を期待しない方が安全です。

#### `cos_mode`

- `point`: 各点の緯度で $\cos$ を掛ける
- `ref`: 基準緯度で $\cos$ を掛ける

#### `unit`

`arcsec`, `arcmin`, `deg`

#### AltAz データの制限

AltAz から ICRS などへの変換は、site/time が十分無いとできません。

#### 返り値

`OFS_LON`, `OFS_LAT` の DataFrame を返します。

#### 実用例

マッピング座標のずれ確認、位置グループ化前の sanity check に向きます。

#### どんなときに有用か

同じ scan に見えて実は中心がずれている場合などを見抜きやすいです。

### `merge_scantables()`

#### 何をするか

複数の `Scantable` を高速に結合します。VLA と固定長の混在にも対応します。

#### 引数

- `sort_by_time`
- `shift_scan_id`

#### `shift_scan_id=True` が重要な理由

複数ファイルを結合すると `SCAN` が衝突しやすいので、自動でオフセットを足して重複を避けます。

#### 時刻ソートの方法

内部で解決したタイムスタンプを使い、`sort_by_time=True` なら時間順に並べ替えます。

#### meta / history の扱い

meta は最初の入力を基準にし、history は簡潔な merge 情報へ整理されます。

#### 実用例

日をまたいだ観測をまとめる前処理に便利です。

#### 特に有効な場面

scan ID が各ファイルで 0 から振り直されている場合でも安全です。

### `update_metadata()`

#### 何をするか

指定行の table 列を in-place で安全に更新します。

#### 主な引数

- `column`
- `value`
- `rows`
- `force`
- `verbose`

#### 内蔵されている安全装置

##### 1. 危険列のブロック

一部の重要列は `force=False` では簡単に変えられません。

##### 2. enum 値の検証

`SPECSYS`, `TEMPSCAL`, `POLARIZA` などは許容値チェックがあります。

##### 3. `TEMPSCAL='TR*'` と `BEAMEFF` の整合性確認

`BEAMEFF` が無いのに `TR*` を設定しようとすると、通常は止めます。

#### alias 同期

`RESTFREQ` を変えたら `RESTFRQ` も同期します。逆も同様です。

#### 新規列の作成

列が無ければ作ります。

#### 破壊的変更であることに注意

この関数は新しい `Scantable` を返さず、その場で変更します。

#### 実用例: OBJECT の修正

観測対象名の typo 修正に向きます。

#### 実用例: RESTFRQ を強制修正

誤った静止周波数で保存されたファイルに後から補正を入れるときに使えます。

#### どんなときに使うべきか

読み込み後に table のメタ情報だけを安全に直したいときです。

### `set_beameff()`

#### 何をするか

`BEAMEFF` を row-wise に設定します。

#### 実装上の位置づけ

`tempscale.set_beameff()` の util 版です。

#### 重要: 何をしないか

自動でデータ配列を TA* / TR* 変換しません。あくまで metadata 更新です。

#### 典型的な使い方

##### 全行へ同じ効率を付与

単一受信機で全行同じ効率のとき。

##### 一部行だけに付与

ビームや IF によって効率が違うとき。

##### 配列で行ごとに別値を入れる

校正結果が行ごとに違うとき。

#### FDNUM / IFNUM / PLNUM ごとに異なる `BEAMEFF` を入れる

最新版でもこの運用は有効です。`filter_scantable()` や `find_scans()` と組み合わせます。

#### 使ったあとに確認する

`TEMPSCAL`, `BEAMEFF`, `RESTFRQ` を `show_scantable()` で見ておくと安全です。

#### `TEMPSCAL` との関係

`TR*` で扱うなら `BEAMEFF` が必要です。

### `find_scans()`

#### 何をするか

条件に一致する行を検索し、index を返します。

#### 条件指定の方法

##### 1. `query`

pandas query 風の文字列

##### 2. `kwargs`

列名ごとの一致条件

#### `extra_data`

一時的に追加した列を query に使えます。

#### 列名の扱い

大文字小文字の揺れをある程度吸収します。

#### エラー時

不正な query は例外になります。

### `filter_scantable()`

#### 何をするか

条件に合う行だけを残した新しい `Scantable` を返します。

#### 特徴

data, table, history の整合を保ったまま切り出します。

#### `rows` 併用時の挙動

まず query / kwargs で絞り、その後に位置インデックス指定を適用します。

#### 実用例

特定の `SCAN`, `FDNUM`, `IFNUM`, `POLARIZA` のみを抜くとき。

#### どんなときに有用か

coadd 前に不要行を落とす前処理に向きます。

### `select_rows()`

#### 何をするか

位置インデックス指定だけで新しい `Scantable` を作ります。

#### 向いている用途

既に index 列を持っていて、明示的にその行だけ抜きたいときです。

### `scantable_utils.py` をどう使うと効果的か

- まず `describe_columns()` で中身を知る  
- 次に `show_scantable()` で目視確認する  
- 必要なら `filter_scantable()` や `update_metadata()` で整える  
- その後 `baseline` / `regrid` / `coadd` に進む  

という使い方が安全です。

## 8.11 `axis.py`

### 公開 API

- `freq_axis_from_wcs()`
- `radio_velocity_kms()`
- `wcs_slice_channels()`
- `vlsrk_axis_from_freq_meta()`
- `channel_slice_from_vrange_union()`
- `slice_channels()`

### 役割

スペクトル軸の生成と、チャネル切り出し時の WCS 更新を担います。

重要式は

$$
\nu_i = \mathrm{CRVAL1} + \left((i+1)-\mathrm{CRPIX1}\right)\mathrm{CDELT1}
$$

および radio definition

$$
v = c\frac{\nu_{\rm rest}-\nu}{\nu_{\rm rest}}
$$

です。

`vlsrk_axis_from_freq_meta()` は、まず周波数基準系を Doppler factor で補正してから速度へ変換します。

`wcs_slice_channels()` と `slice_channels()` は、チャネルを `[s:e]` で切るとき

$$
\mathrm{CRVAL1}_{\rm new} = \mathrm{CRVAL1}_{\rm old} + s\,\mathrm{CDELT1}
$$

とし、`CRPIX1` はそのまま保つ安全な規約です。

## 8.12 `ranges.py`

### 公開 API

- `parse_windows()`
- `window_to_mask()`
- `windows_to_mask()`

### 使い方

`"-35:-5"` のような文字列を `(-35.0, -5.0)` へ直します。境界は常に小さい方から大きい方へそろえます。

## 8.13 `doppler.py`

### 公開 API

- `earth_location_from_meta()`
- `calc_vlsrk_correction_kms()`
- `compute_vcorr_series()`
- `get_doppler_factor()`
- `scale_frequency_wcs_by_velocity()`

### 役割

観測地点・観測時刻・視線方向から、LSRK 補正速度を計算します。

### `earth_location_from_meta()` の探索順

1. `OBSGEO-X/Y/Z`  
2. `SITELAT`, `SITELONG`, `SITEELEV`  
3. `site_name`  

です。

### `compute_vcorr_series()`

`DatetimeIndex` と各行の座標列から、行ごとの $v_{\rm corr}$ を km/s で返します。

## 8.14 `profile_view`（旧 plotting）

### モジュール

このアーカイブには `profile_view` の実装ファイル自体は含まれていませんが、`__init__.py` では `view_spectra`, `plot_profile_map` を公開 API として import しています。

### 役割分担

 package 全体としては quick-look と可視化を担当する層です。

### 共通思想

軸情報を行ごとに尊重し、解析結果を人間が確認できる形へ落とします。

### 代表コンストラクタ

添付アーカイブに実装本体が無いため、今回の更新ではここを package-level の位置づけ説明に留めます。可視化モジュールの内部仕様を更新したい場合は、そのソースを含む版を別途参照してください。


## 8.15 `sdfits_writer.py`

### 公開 API

- `Site`
- `Efficiency`
- `SpectralAxisUniform`
- `DatasetInfo`
- `SDRadioSpectralSDFITSWriter`

### 位置づけ

`fitsio.write_sdfits()` が解析側の推奨保存 API であるのに対し、`sdfits_writer.py` は観測側・converter 側も含めた低レベル writer の本体です。最新版でも、最終的な BinTable 構築は `sdfits_bintable.py` に委譲されます。

### `Site`

観測地点を

- `lat_deg`
- `lon_deg`
- `elev_m`

で表し、`EarthLocation` と `OBSGEO-X/Y/Z` の両方へ落とせます。

### `Efficiency`

- `beameff`
- `apereff`
- `mooneff`
- `suneff`
- `effstat`

を保持します。`effstat` は `ASSUMED`, `MEASURED`, `UNKNOWN` の enum です。

### `SpectralAxisUniform`

等間隔スペクトル軸の既定値を保持します。周波数軸の基本式は

$$

u_i = \mathrm{CRVAL1} + \left((i+1)-\mathrm{CRPIX1}
ight)\mathrm{CDELT1}
$$

です。`velref()` は古い AIPS / casacore 互換の `VELREF` 整数を返します。

### `DatasetInfo`

望遠鏡名、観測者、プロジェクト名、対象名、座標系、ビーム情報、効率、共通 spectral axis などを一括保持します。

### `SDRadioSpectralSDFITSWriter`

#### 何をするか

行ごとにスペクトルと metadata を蓄積し、必要なら chunk 分割しながら SDFITS を書きます。

#### 重要な初期化引数

- `n_chan`: 固定長の既定チャネル数または最大値の目安
- `site`
- `info`
- `store_freq_column`: `FREQ` ベクトル列を持つか
- `chunk_size`: parts 分割の行数
- `out_basename`: parts 出力の基底名
- `history`

#### 列ポリシー

writer は列を 4 群で扱います。

1. always  
   行の意味やスペクトル軸解釈に必須な列
2. context-required  
   文脈によって必須になる列
3. optional  
   meaningful な値があるときだけ出す列
4. block-optional  
   ブロック内のどれかが meaningful ならブロック全体を出す列

#### always の典型例

- `TIME`, `MJD`, `DATE-OBS`, `DATEOBS`, `TIMESTAMP`
- `SCAN`, `SUBSCAN`, `INTGRP`
- `OBJECT`, `OBSMODE`
- `EXPOSURE`, `CALSTAT`, `TEMPSCAL`, `FLAGROW`
- `RA`, `DEC`, `GLON`, `GLAT`
- `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`, `SPECSYS`
- `FDNUM`, `IFNUM`, `PLNUM`, `POLARIZA`
- `DATA`, `FLAG`

#### context-required の典型例

- `RESTFREQ`, `VELDEF`
- sideband 解釈が必要なときの `OBSFREQ`, `SIDEBAND`

#### optional の典型例

- `TCAL`, `THOT`, `TSYS`, `TAU0`
- `AZIMUTH`, `ELEVATIO`
- `VFRAME`, `FOFFSET`
- `BORE_AZ`, `BORE_EL`, `BEAMXOFF`, `BEAMYOFF`, `BEAMROT`
- `FRONTEND`, `BACKEND`, `SAMPLER`
- `IMAGFREQ`, `LO1FREQ`, `LO2FREQ`, `LO3FREQ`
- `SB1`, `SB2`, `SB3`
- `FREQ`

#### block-optional の典型例

- source block  
  `SRCFRAME`, `SRCRDSYS`, `SRCEQNX`, `SRC_LONG`, `SRC_LAT`
- scan block  
  `SCANFRAM`, `SCANRDSYS`, `SCANEQNX`, `SCANX`, `SCANY`
- pointing block  
  `AZ_CMD`, `EL_CMD`, `CORR_AZ`, `CORR_EL`, `CALC_REFR`
- weather block  
  `TAMBIENT`, `PRESSURE`, `HUMIDITY`, `WINDSPD`, `WINDDIR`

#### parts 出力

`chunk_size` を指定すると、大きな観測を複数の parts へ分割できます。`out_basename` を基底名にして、複数ファイルとして保存する設計です。

## 8.16 `sdfits_bintable.py`

### 公開 API

- `set_meta_keyword()`
- `apply_meta_to_header()`
- `VectorColumnSpec`
- `build_history_hdu()`
- `build_single_dish_table_hdu()`

### 役割

SDFITS の BinTable を、固定長配列と VLA の両方に対応して安全に構築する、Single Source of Truth です。

### `set_meta_keyword()`

非標準キーワードや HIERARCH を含むメタデータを、できるだけ安全に FITS header へ入れます。

### `build_history_hdu()`

`HISTORY` 拡張を key/value テーブルとして作ります。最新版のこのパッケージでは、長い provenance を header card の `HISTORY` に詰め込むより、専用 BinTable に逃がす方針です。

### `build_single_dish_table_hdu()`

#### 何をするか

行 metadata, spectrum, optional vector columns を見て、最終的な `SINGLE DISH` BinTableHDU を作ります。

#### 重要な引数

- `table`
- `data`
- `spectrum_column`
- `include_flag`
- `flag`
- `extra_vector_columns`
- `units`
- `normalize_columns`
- `string_widths`
- `prefer_float64_vectors`
- `allow_vla_vectors`
- `bunit`
- `warn_dropped_vector_columns`
- `extname`

#### fixed-length と VLA の判定

- 全行が同じ長さなら fixed-length vector
- 行ごとに長さが違えば VLA

です。

#### spectrum / flag / freq の整合

- `DATA` と `FLAG` は行ごとに同長でなければならない
- `FREQ` を出すなら `DATA` と行ごとに同長でなければならない

#### vector-in-cell 列

各セルに `list` や `ndarray` が入っている列は、fixed-length か VLA vector 列として保存されます。`BSL_COEF` のような列がこれに該当します。

#### 文字列列

`string_widths` があればその幅を使い、無ければ内容から自動推定します。

# 9. obs_method.txt を一般化した推奨解析フロー

## 9.1 12CO の典型フロー

1. `read_scantable()` または `load_rawspec_auto()` で読む  
2. `run_tastar_calibration()` で Ta* 化する  
3. `run_velocity_coadd(group_mode="scan")` で scan 内積分  
4. `run_baseline_fit()` で baseline 減算  
5. `run_velocity_coadd(group_mode="position")` で位置積分  
6. 必要なら `post_baseline_mode='inherit_all'` 付きで仕上げる  

## 9.2 13CO / C18O のように同じ raw から複数線を取りたい場合

raw の周波数軸は同じでも、静止周波数だけを変えて別々に解釈します。FREQ 軸のときは WCS を変えずに `RESTFRQ` だけを差し替えるのが正しい流儀です。

# 10. よくある設計判断と推奨

## 10.1 baseline は scan coadd の前か後か

scan ごとにまずまとめてから baseline を引く方が S/N は上がりますが、scan 内でベースラインが大きく変わるなら先に個別 baseline が有利です。最新版では両方可能です。

## 10.2 `baseline_vwin` と `rms_vwin` の使い分け

- baseline も引いて、その residual RMS で重み付けしたいなら `baseline_vwin`
- baseline は引かず、line-free 窓の RMS だけ測りたいなら `rms_vwin`

です。同時指定は不可です。

## 10.3 TOPOCENT のまま baseline をしてよいか

はい。ただし velocity 窓で指定する以上、各行の物理軸を正しく作れるだけの `VELOSYS` / `VFRAME` または再計算可能な時刻・座標・サイト情報が必要です。

## 10.4 viewer 用に TR* へ変換したい

読み込み時に自動変換はしません。`write_sdfits(..., tempscal="TR*", data_scale="TA*")` のように保存時だけ変換するか、明示的に `ta_to_tr()` を使います。

# 11. エラーや停止条件の読み方

## 11.1 `TOPOCENT input requires VELOSYS/VFRAME or valid timestamps`

TOPO 入力なのに、未適用補正列も valid timestamp も無い状態です。

## 11.2 `Input has SPECSYS=LSRK but contains non-zero unapplied velocity column`

LSRK と言っているのに未適用補正が残っています。入力 metadata の矛盾です。

## 11.3 `Cannot specify both baseline_vwin and rms_vwin`

重み付け評価法を二重指定しています。

## 11.4 `RESTFRQ missing`

速度窓解釈や FREQ→VRAD 変換に必要な静止周波数がありません。

## 11.5 TR* 変換時の `BEAMEFF` エラー

`BEAMEFF` が無い、非有限、または 0 以下です。

# 12. 実務上の推奨チェックリスト

- `SPECSYS`, `SSYSOBS`, `CTYPE1`, `CUNIT1`, `RESTFRQ` を最初に確認する  
- `TEMPSCAL` と `BEAMEFF` を確認する  
- TOPO 入力なら `VELOSYS` / `VFRAME` または timestamp/coord/site の有無を確認する  
- baseline 前に line-free 窓が本当に line-free かを確かめる  
- coadd 後に `N_IN`, `N_USED`, `WGT_*`, `BSL_*` を見る  
- 保存時に on-disk スケールを変えるなら `data_scale` を必ず明示する  

# 13. 最低限の API 早見表

## 13.1 解析で最もよく使う関数

- `read_scantable()`
- `run_tastar_calibration()`
- `run_baseline_fit()`
- `run_velocity_regrid()`
- `run_velocity_coadd()`
- `write_scantable()`

## 13.2 補助的だが重要

- `apply_restfreq_override()`
- `set_beameff()`
- `calc_mapping_offsets()`
- `merge_scantables()`
- `filter_scantable()`

# 14. 付録 A: 公開 API シグネチャ一覧

## 14.1 I/O

以下は最新版コードから照合した公開シグネチャです。`*` は keyword-only 引数を表します。

```python
class Scantable(meta, data, table, history={})
read_scantable(path, *, tr_input_policy='preserve')
write_scantable(path, scantable, spectrum_column='DATA', overwrite=True, **kwargs)
read_tastar_fits(path)
write_sdfits(out_path, meta, data, table, history=None, *, spectrum_column='DATA', include_flag=True, overwrite=True, string_widths=None, **kwargs)
```

## 14.2 Raw

```python
class RawSpec(meta, hot, off, on, mapping)
build_rawspec(hot, on, off, meta, mapping)
save_rawspec(raw, path)
load_rawspec(path)
load_rawspec_auto(path, prefer=None)
load_rawspec_fits(path)
```

## 14.3 Calibration

```python
make_tastar_dumps(raw, *, rows=None, exclude_rows=None, t_hot_k=None, tau_zenith=None, t_surface_k=None, t_atm_model='offset', t_atm_delta_k=15.0, t_atm_eta=0.95, gain_mode='hybrid', verbose=True, ch_range=None, vlsrk_range_kms=None, v_corr_col='VFRAME', coord_frame=None, vcorr_chunk_sec=None, dtype=None, rest_freq=None)
tastar_from_rawspec(raw, t_hot_k=None, *, ch_range=None, vlsrk_range_kms=None, vcorr_chunk_sec=None, dtype=None, rest_freq=None, v_corr_col='VFRAME', coord_frame=None, rows=None, exclude_rows=None, tau_zenith=None, t_surface_k=None, t_atm_model='offset', t_atm_delta_k=15.0, t_atm_eta=0.95, gain_mode='hybrid', verbose=True)
run_tastar_calibration(input_data, output_path=None, t_hot_k=None, ch_range=None, vlsrk_range_kms=None, coord_frame=None, spectrum_column='DATA', overwrite=False, store_freq_column=False, vcorr_chunk_sec=None, dtype=None, rest_freq=None, v_corr_col='VFRAME', rows=None, exclude_rows=None, tau_zenith=None, t_surface_k=None, t_atm_model='offset', t_atm_delta_k=15.0, t_atm_eta=0.95, gain_mode='hybrid', verbose=True)
recalibrate_tastar(scantable, new_tau=None, new_t_surface_k=None, new_t_atm_model='offset', new_t_atm_delta_k=15.0, new_t_atm_eta=0.95, verbose=True)
```

## 14.4 Baseline

```python
fit_polynomial_baseline(v_kms, y, base_windows=None, *, windows=None, line_windows=None, poly_order=1, order=None, iter_max=0, iter_sigma=3.0)
run_baseline_fit(input_data, output_path=None, *, rows=None, exclude_rows=None, vwin, poly_order=1, line_vwin=None, iter_max=0, iter_sigma=3.0, max_dumps=0, ch_start=None, ch_stop=None, v_corr_col='VELOSYS', rest_freq=None, apply=True, bsl_overwrite='replace', on_fail='exit', overwrite=True)
```

## 14.5 Regrid / Standardizer

```python
class VGrid(v0_kms, dv_kms, nchan, crpix1=1.0)
make_vgrid(vmin_kms, vmax_kms, dv_kms)
get_axis_signature(row, meta, nchan)
class Standardizer(scantable, target_grid=None, v_corr_col='VFRAME', auto_grid_mode='stable_native', auto_grid_anchor_kms=0.0, auto_grid_tol_frac=1e-06)
vlsrk_axis_for_spectrum(meta, *, v_corr_kms, nchan=None)
interp_to_vgrid(v_src, y_src, v_tgt)
vrange_from_meta_and_vcorr(meta, v_corr_kms, *, nchan=None)
run_velocity_regrid(input_data, output_path=None, *, rows=None, exclude_rows=None, vmin_kms, vmax_kms, dv_kms, v_corr_col='VFRAME', rest_freq=None, ch_start=None, ch_stop=None, max_dumps=0, overwrite=True, fill_value=np.nan, keep_row_order=True, drop_allnan_rows=False, history_tag='velocity_regrid')
```

## 14.6 Coadd

```python
run_velocity_coadd(inputs, output_path=None, *, rows=None, exclude_rows=None, mode='uniform', group_mode='position', pos_col='pos_id', pos_tol_arcsec=None, v_corr_col='VELOSYS', coord_frame=None, vcorr_chunk_sec=None, vmin=None, vmax=None, dv=None, allow_outside_overlap=False, axis_type='freq', rest_freq=None, baseline_vwin=None, baseline_poly=0, baseline_iter_max=0, baseline_iter_sigma=3.0, post_baseline_mode=None, post_baseline_vwin='inherit', post_baseline_poly='inherit', post_baseline_iter_max='inherit', post_baseline_iter_sigma='inherit', rms_vwin=None, rms_poly=1, rms_bin=1, weight_zero_policy='error', coadd_qc=None, line_vwin=None, block_size=0, max_dumps=0, ch_start=None, ch_stop=None, on_fail='exit', overwrite=True, out_scale='TA*', normalize_if_mixed='auto', beameff_tol=1e-06, sigma_scale='TA*', verbose=True)
```

## 14.7 Utilities

```python
describe_columns(sc)
show_scantable(inputs, rows=None, columns='default', head=20, show_legend=False, extra_data=None, ref_coord=None, frame='ICRS', projection='GLS', unit='arcsec')
calc_mapping_offsets(sc, ref_coord=None, frame='ICRS', projection='GLS', unit='arcsec', cos_mode='point', verbose=True)
merge_scantables(inputs, sort_by_time=False, shift_scan_id=True)
update_metadata(sc, column, value, rows=None, force=False, verbose=True)
set_beameff(sc, efficiency, rows=None, verbose=True)
find_scans(sc, query=None, extra_data=None, **kwargs)
filter_scantable(sc, query=None, extra_data=None, rows=None, **kwargs)
select_rows(sc, rows)
normalize_tempscal(value, *, default='TA*')
apply_restfreq_override(meta, table, rest_freq_hz, require_wcs_for_vrad=True)
compute_vcorr_series(times, ra_deg, dec_deg, meta, coord_frame='icrs')
```

# 15. 最後に

この最新版で特に重要なのは、

- 読み込み時の温度スケール非破壊化
- legacy 速度列の自動移行
- `run_velocity_regrid()` の追加
- `run_velocity_coadd()` における `weight_zero_policy`, `normalize_if_mixed`, `post_baseline_mode='inherit_all'`
- 保存時だけの明示的 TA* / TR* 変換

です。

旧版からの最大の思想的変化は、「読んだ瞬間に勝手に直す」のではなく、「意味を保ったまま読み、必要な変換は明示的に行う」へ寄ったことです。これは、温度スケール、静止周波数、速度基準系のいずれでも同じです。


# 16. 付録 B: 解析後に現れやすい列の意味

## 16.1 baseline 系

- `BSL_DONE`: baseline 評価が完了したか  
- `BSL_APPLIED`: 実際に減算したか  
- `BSL_STAGE`: `baseline_fit` か `post_coadd` か  
- `BSL_POLY`: 多項式次数  
- `BSL_WINF`: 窓文字列  
- `BSL_RMS`: line-free 領域の標準偏差  
- `BSL_STAT`: 現状は `std`  
- `BSL_NUSED`: フィットに使った点数  
- `BSL_COEF`: 多項式係数列  
- `BSL_SCALE`: その baseline をどの温度スケールで評価したか  

## 16.2 weight 系

- `WGT_MODE`: `inverse_variance` または `uniform`  
- `WGT_SRC`: `baseline_vwin`, `rms_vwin`, `input_bsl_rms`, `uniform`  
- `WGT_VWIN`: 重み評価窓  
- `WGT_POLY`: RMS 推定前に使った baseline 次数  
- `WGT_BIN`: RMS 評価前の spectral binning 幅  
- `WGT_STAT`: 現状は `std`  
- `WGT_INPUT_SUB`: 入力側で baseline subtraction をしたか  
- `WGT_ZERO_POLICY`: 0 または不正重みの扱い  

## 16.3 regrid / coadd 系

- `REGRID_DONE`: 行保持 regrid を通したか  
- `REGRID_FRAME`: 現状は `LSRK`  
- `N_IN`: グループに入った本数  
- `N_USED`: 実際に使えた本数  
- `WEIGHT_MODE`: 人間向けの重みモード文字列  
- `VELOSYS_OBS`, `VFRAME_OBS`: 観測時の未適用補正退避列  

# 17. 付録 C: 式のまとめ

## 17.1 周波数 WCS

$$

u_i = \mathrm{CRVAL1} + \left((i+1)-\mathrm{CRPIX1}
ight)\mathrm{CDELT1}
$$

## 17.2 radio velocity

$$
v = c \frac{\nu_{\mathrm{rest}} - \nu}{\nu_{\mathrm{rest}}}
$$

ここで $c$ は光速、$\nu_{\mathrm{rest}}$ は静止周波数、$\nu$ は観測された周波数です。$\nu < \nu_{\mathrm{rest}}$ なら $v > 0$ で、radio definition では後退を正に取ります。

## 17.3 relativistic Doppler factor

$$
k = \sqrt{\frac{1+\beta}{1-\beta}}, \qquad \beta = \frac{v}{c}
$$

$$
\nu_{\mathrm{LSRK}} = \frac{\nu_{\mathrm{obs}}}{k}
$$

この形は、視線速度 $v$ に対応する相対論的 Doppler factor $k$ を使って、観測周波数 $\nu_{\mathrm{obs}}$ を基準系補正後の周波数 $\nu_{\mathrm{LSRK}}$ へ写す書き方です。

## 17.4 Ta* / Tr* 変換

$$
T_{R}^{*} = \frac{T_{A}^{*}}{\eta_{\mathrm{mb}}}
$$

$$
T_{A}^{*} = T_{R}^{*} \eta_{\mathrm{mb}}
$$

ここで $\eta_{\mathrm{mb}}$ は main-beam efficiency です。

## 17.5 VRAD 軸での静止周波数上書き

$$
a = \frac{\nu_1}{\nu_2}
$$

$$
b = c \left(1 - \frac{\nu_1}{\nu_2}\right)
$$

$$
v_2 = a v_1 + b
$$

これは、同じ周波数軸を異なる静止周波数 $\nu_1$, $\nu_2$ で radio velocity として解釈したとき、速度軸がアフィン変換になることを表しています。

## 17.6 大気モデル付き等価較正温度

$$
X = \frac{1}{\sin(\mathrm{El})}
$$

$$
T_{\mathrm{atm}} = T_{\mathrm{surf}} - \Delta T
$$

または

$$
T_{\mathrm{atm}} = \eta T_{\mathrm{surf}}
$$

で、実装上の等価較正温度は

$$
T_{\mathrm{cal}} = T_{\mathrm{hot}} e^{\tau X} - T_{\mathrm{atm}} \left(e^{\tau X}-1\right) - T_{\mathrm{bg}}
$$

です。

## 17.7 inverse-variance coadd

$$
w_i = \frac{1}{\sigma_i^2}
$$

$$
T_{\mathrm{out}}(k) = \frac{\sum_i w_i T_i(k)}{\sum_i w_i}
$$

# 18. 付録 D: 実用メモとパラメータ解釈の補足

## 18.1 `rows` / `exclude_rows` の考え方

多くの高レベル関数は `rows` と `exclude_rows` を持ちます。これは DataFrame の index label ではなく、現在の `Scantable.table` に対する位置インデックスです。したがって、前段で `filter_scantable()` した後では番号が変わり得ます。

## 18.2 `max_dumps` の意味

`max_dumps` は主にデバッグや quick-look 用で、選択後の先頭から何行まで使うかを制限します。統計量の偏りを避けたい本番解析では、任意の一部分だけ切るより、明示的な `rows` 指定の方が意味が分かりやすいです。

## 18.3 `ch_start` / `ch_stop` の意味

チャネル切り出しは 0-based で、`ch_start` は inclusive, `ch_stop` は exclusive です。切り出し後の WCS は

$$
\mathrm{CRVAL1}_{
m new} = \mathrm{CRVAL1}_{
m old} + ch_{
m start}\,\mathrm{CDELT1}
$$

として更新されます。`CRPIX1` は不変です。

## 18.4 `axis_type='freq'` の意味

`run_velocity_coadd()` の `axis_type='freq'` は、「周波数ベースの入力をそのまま保存する」という意味ではありません。内部では一度 LSRK 共通速度グリッドへそろえてから、その速度グリッドを LSRK の `FREQ` 軸へ戻しています。したがって出力 `SPECSYS` は依然として `LSRK` です。

## 18.5 `sigma_scale` の現状

最新版では `sigma_scale` は公開引数として存在しますが、実装は `TA*` のみを受け付けます。つまり、RMS 評価や重み付けは Ta* スケール基準で行う設計で、最終出力だけ `out_scale='TR*'` へ変換する、というのが現在の正しい読み方です。

## 18.6 `weight_zero_policy` の使い分け

- `error`: 0 または不正重みが出たら止める。最も安全  
- `drop`: 問題のある入力スペクトルをそのグループから除外する  
- `impute_median`: 代表重みで置き換える。探索的解析向け  

特に論文用の最終スペクトルでは、まず `error` で入力の問題を見つける方が安全です。

## 18.7 `normalize_if_mixed` の考え方

異なる `BEAMEFF` を持つ行をそのまま TA* で足すと、実質的に異なるビーム効率でスケールされた量を混ぜることになります。`auto` はこの混合を避けるための安全装置です。出力を TA* にしたい場合でも、内部では一度 TR* へ正規化してから代表効率で戻す、という考え方です。

## 18.8 `post_baseline_mode='inherit_all'` を使うときの注意

これは便利ですが、前段の `baseline_vwin` が妥当であることを前提にします。前段窓が広すぎたり線を食っていると、その設定がそのまま post-coadd にも伝播します。したがって、まず pre-coadd の `BSL_RMS` と `BSL_COEF` を見て、窓設定が妥当か確認するのがよいです。

## 18.9 `rest_freq` 上書きの使いどころ

- 同じ raw から別分子線を解析したい  
- 既存ファイルの `RESTFRQ` が誤っていた  
- velocity window の解釈だけを別線基準に切り替えたい  

というときに使います。FREQ 軸で WCS を書き換えないことが大事です。

## 18.10 `TEMPSCAL` を metadata だけで変えない

`update_metadata(sc, 'TEMPSCAL', 'TR*')` だけでは、スペクトル配列自体は変わりません。配列の物理量と宣言を一致させるには、明示的に Ta* / Tr* 変換を行う必要があります。逆に、保存時だけ変えたいなら `write_sdfits()` の `tempscal` / `data_scale` を使うのが安全です。

# 19. 付録 E: 入力ファイルに対して期待される最低契約

## 19.1 baseline / regrid / coadd の共通契約

最低限、各行に対して次のどれかが必要です。

1. 周波数 WCS  
   `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`
2. 静止周波数  
   `RESTFRQ` または `RESTFREQ`
3. 基準系  
   `SPECSYS` または `SSYSOBS`
4. TOPO 行なら未適用補正列、またはそれを再計算するための時刻・座標・サイト情報

## 19.2 速度補正再計算に必要なもの

- 時刻列  
  `TIMESTAMP`, `MJD`, `DATE-OBS`, `DATEOBS`, legacy `TIME` のどれか
- 座標列  
  `RA`, `DEC` あるいは `GLON`, `GLAT`
- サイト情報  
  `OBSGEO-X/Y/Z` または `SITELAT`, `SITELONG`, `SITEELEV`

## 19.3 温度スケール変換に必要なもの

`TR*` を扱うには `BEAMEFF` が必要です。行列でも header でもよいですが、最終的に行ごとへ展開できる必要があります。

## 19.4 VLA の契約

- `data` は `List[np.ndarray]` でよい  
- 各行の長さが違ってよい  
- ただし、その行を解釈する WCS と静止周波数は行ごとに一意でなければならない  

## 19.5 legacy 列の扱い

- `V_CORR_KMS` は読込時に `VFRAME`, `VELOSYS` へ移行される  
- `TAMB_K` は `THOT` へ rename される  
- `TAU` は `TAU0` へ rename される  

ただし、複数の legacy 列が互いに矛盾しているときは自動移行は停止します。

# 20. 付録 F: 高レベル関数をどう組み合わせるか

## 20.1 最も保守的な流れ

1. `read_scantable()` で読む  
2. `describe_columns()` で列を確認する  
3. `show_scantable()` で `SPECSYS`, `RESTFRQ`, `TEMPSCAL`, `BEAMEFF`, `TIMESTAMP` を見る  
4. 必要なら `update_metadata()` または `set_beameff()` で整える  
5. `run_baseline_fit(apply=False)` で baseline 評価だけ先に行う  
6. `run_velocity_regrid()` で共通軸へそろえた quick-look を見る  
7. 問題なければ `run_velocity_coadd()` を本番設定で実行する  

## 20.2 速度窓だけで channel slice を決めたい場合

まず `vlsrk_range_kms` を使って較正段階で切るか、`run_velocity_regrid()` で明示的な `vmin_kms`, `vmax_kms`, `dv_kms` を与えるのが安全です。WCS の途中で人手でチャネル番号へ変換してしまうより、物理窓をそのまま関数へ渡す方が符号ミスを減らせます。

## 20.3 file ごとに baseline 後、最後に全統合する理由

観測日ごとに gain や baseline の癖が違うとき、全ファイルを先に一つへ混ぜると問題の切り分けが難しくなります。最新版の API は `merge_scantables()` や `run_velocity_coadd(inputs=[...])` に対応しているので、file 内処理と全体統合を分けやすいです。

## 20.4 quick-look と最終保存を分ける

quick-look では `max_dumps`, `drop_allnan_rows`, `fill_value` を使って軽く見て、本番保存ではそれらを保守的設定へ戻す、という使い分けがしやすくなっています。


# 21. 付録 G: ファイル形式と履歴の扱い

## 21.1 `read_tastar_fits()` が読める 2 系統のファイル

最新版の `read_tastar_fits()` は、少なくとも次の 2 系統を扱います。

### A. `SINGLE DISH` テーブルを核にした SDFITS-like 形式

- `PRIMARY`
- `SINGLE DISH`
- 必要なら `HISTORY`

という構成です。スペクトル列は `DATA` または `SPECTRUM` を探します。`FLAG` や `FLAGS` は table から除外し、`table` には scalar 列や vector-in-cell 列を入れます。

### B. 旧内部形式

- `DATA` ImageHDU
- `DUMPS` BinTableHDU
- 必要なら `HISTORY`

です。旧形式も読み続けられるように残していますが、今後の推奨保存先は `SINGLE DISH` 形式です。

## 21.2 `NAXIS1`, `NCHAN`, `NCHANSEL` の意味

固定長配列では、通常

- `NAXIS1`: 行ごとのチャネル数
- `NCHAN`: 選択後チャネル数
- `NCHANSEL`: 実効選択チャネル数

が一致します。

VLA では行ごと長さが違うため、ヘッダ上の `NCHAN` や `NCHANSEL` は代表値または最大長のヒントに近く、真の長さは各行の vector 長を見て判断する必要があります。

## 21.3 `history` の書き方

このパッケージの `history` は辞書、あるいは stage ごとの辞書や list を含む構造です。保存時には `build_history_hdu()` により key/value テーブルへ展開されます。

### baseline の履歴

`run_baseline_fit()` は `baseline_history` という list を持ち、そこへ stage, 作成時刻, 入力名, 窓設定, 反復設定, `v_corr_col`, `rest_freq`, `apply`, `bsl_overwrite` などを書きます。

### regrid の履歴

`run_velocity_regrid()` は既定で `history['velocity_regrid']` に

- `stage='velocity_regrid'`
- `frame='LSRK'`
- `kind='row_preserving_regrid'`
- `vmin_kms`, `vmax_kms`, `dv_kms`
- `source_specsys`
- `source_vcorr_col`
- `drop_allnan_rows`

などを入れます。

### coadd の履歴

`run_velocity_coadd()` は履歴に

- 出力スケール
- `normalize_if_mixed`
- `beameff_tol`
- `sigma_scale`
- `mode`, `group_mode`
- `rms` 設定
- `baseline` 設定
- `post_baseline` 設定
- `weight_zero_policy`

を入れます。さらに混在 `BEAMEFF` を自動正規化した場合は、scale history に `beameff_mixed_detected` が追加されます。

## 21.4 PRIMARY と `SINGLE DISH` の役割分担

- `PRIMARY`: ファイル共通の既定値  
- `SINGLE DISH`: 行ごとの真の値  

です。`write_scantable()` が重要列を meta から table へ昇格させるのは、解析側で「行ごとに解釈可能であること」を優先しているからです。

## 21.5 `HISTORY` 拡張を使う理由

FITS header の `HISTORY` card だけでは、長い JSON 風の provenance や辞書構造を扱いにくく、カード長制限にも当たりやすいです。そのため最新版では専用 `HISTORY` BinTable を使う方針です。

## 21.6 非標準 meta の保存

`write_sdfits()` は header に安全に入らない meta を捨てません。`META:<KEY>` という形で `history` 側へ逃がして diagnostic として残します。これにより、header card 長や FITS キーワード制約により情報が失われることを減らしています。
