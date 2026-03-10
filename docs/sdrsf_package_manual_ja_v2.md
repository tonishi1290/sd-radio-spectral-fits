# 単一鏡電波分光解析パッケージ 説明書
## 解析パイプライン・I/O・基準系・温度スケール・coadd・Standardizer 総合ガイド

> 対象  
> 添付されたソース一式（`fitsio.py`, `rawspec.py`, `calibrate.py`, `baseline.py`, `regrid_vlsrk.py`, `coadd.py`, `restfreq.py`, `tempscale.py`, `scantable_utils.py`, `profile_view/*` ほか）に基づく説明書です。  
> `obs_method.txt` の解析例は、説明書では **一般的で再利用しやすい形** に整理し、**ループを外した単発例** を中心に記述しています。

---

## 1. このパッケージは何をするものか

このパッケージは、単一鏡の分光データを

1. 生データ読み込み  
2. Ta* キャリブレーション  
3. 必要なら rest frequency の上書き  
4. ベースライン処理  
5. LSRK へそろえた速度グリッドへの再配置  
6. scan 単位・position 単位・INTGRP 単位の coadd  
7. 表示・確認・SDFITS 書き出し  

まで、一貫して扱うための解析基盤です。

設計上の重要点は次のとおりです。

- **Scantable** を中心とした一貫したメモリ表現
- **行ごとのメタデータ** を重視し、`RESTFRQ`, `CRVAL1`, `SPECSYS`, `VELOSYS` などを行列と一緒に持てる
- **TOPOCENT → LSRK 補正** の意味論を厳密に分離
- **FREQ 軸と VRAD 軸の扱いを分ける**
- **TEMPSCAL / BEAMEFF** を破壊的に扱わない
- **固定長配列と VLA（可変長スペクトル）** の両対応
- **big-endian FITS** を pandas/Numpy で安全に扱うための変換を内部で実施
- coadd では **Standardizer** により異なる WCS を共通速度グリッドへそろえる
- profile viewer 系では **個別行ごとの軸再計算** を尊重する

---

## 2. まず押さえるべき全体像

### 2.1 解析の標準的な順番

一般的な使用順は次です。

1. `fitsio.read_scantable()` あるいは `rawspec.load_rawspec_auto()` で入力を読む  
2. `calibrate.run_tastar_calibration()` で Ta* を作る  
3. `coadd.run_velocity_coadd(group_mode="scan", ...)` で scan ごとにまとめる  
4. `baseline.run_baseline_fit()` でベースラインを引く  
5. `coadd.run_velocity_coadd(group_mode="position" 相当)` で同位置をまとめる  
6. 必要なら再度 `baseline.run_baseline_fit()`  
7. `profile_view.viewer.view_spectra()` や montage/grid で確認  
8. `fitsio.write_scantable()` / `write_sdfits()` で保存

### 2.2 典型的な思想

- **キャリブレーションは強度の仕事**
- **coadd で基準系補正と再グリッドを消化する**
- **FREQ 軸を持つデータに rest frequency を変えて速度解釈したい場合、WCS をむやみに書き換えない**
- **VRAD 軸を持つデータで rest frequency を変える場合のみ、速度軸 WCS を厳密に更新する**
- **TEMPSCAL が TR* の入力でも、読み込み時に勝手に TA* に変換しない**
- **BEAMEFF が混在しているグループは、必要に応じて TR* 正規化へ逃がす**

---

## 3. このパッケージの中核データ構造

## 3.1 Scantable

`fitsio.Scantable` は、解析済みデータの中核コンテナです。

### 役割
- `meta`: ヘッダ相当の辞書
- `data`: スペクトル本体
- `table`: 行ごとの補助情報
- `history`: 処理履歴

### 特徴
- `data` は通常の 2 次元配列でもよく、VLA 的に **行ごとに長さの違う `List[np.ndarray]`** でもよい
- `copy()` は `meta`, `table`, `history` を浅く、`data` は配列または各行をコピーする
- `write_scantable()` にそのまま渡すと VLA も保存できる

### Scantable を使うと何が良いか
- ベースライン・coadd・viewer が同じデータ型を受け取れる
- `table` に `OBSMODE`, `SCAN`, `RA`, `DEC`, `VELOSYS`, `TEMPSCAL`, `BEAMEFF` などを同居させられる
- 「スペクトル配列」と「その行がどの基準系・どの rest frequency か」が分離されずに済む

---

## 3.2 RawSpec

`rawspec.RawSpec` は、**HOT / OFF / ON の生スペクトル** を扱うコンテナです。

### 構成
- `meta`
- `hot`
- `off`
- `on`
- `mapping`

### 想定用途
- 観測直後の生データを、Ta* キャリブレーション前に扱う
- `OBSMODE=HOT/OFF/ON` が並ぶ SDFITS あるいは pickle 化された内部形式を共通化する
- `mapping` には HOT/OFF/ON を含む全行テーブルを保持する

### 重要
`mapping` は ON だけでなく全行を持つ設計ですが、Ta* キャリブレーションの最終出力は ON 行主体になります。

---

## 4. 重要仕様 1: スペクトル軸・基準系・速度補正

このパッケージで最も重要なのは、`CTYPE1`, `CUNIT1`, `SPECSYS`, `SSYSOBS`, `RESTFRQ`, `VELOSYS`, `VFRAME` の意味づけです。

## 4.1 基本原則

### 軸の物理量
- `CTYPE1 = FREQ` : 周波数軸
- `CTYPE1 = VRAD` : radio definition の速度軸

### 軸の単位
- `FREQ` なら `Hz`
- `VRAD` なら保存時標準は `m/s`

### 基準系
- `SPECSYS` は **現在その軸が属している基準系**
- `SSYSOBS` は **観測時に固定だった基準系**
- Raw / Calibration / Coadd のどこで補正を消化したかはここで判断する

### 静止周波数
- 正本は `RESTFRQ`
- 互換用として `RESTFREQ` も同値で持てる
- 両者は常に一致させる

### 観測者速度補正
- `VELOSYS` は observer → rest-frame の視線速度
- 単位は **m/s**
- `VFRAME` は内部互換用
- `SPECSYS=LSRK` のデータに未消化の `VFRAME` を残してはいけない

---

## 4.2 フェーズごとの意味

### Raw
- ドップラー補正なしなら  
  - `CTYPE1=FREQ`
  - `SPECSYS=TOPOCENT`
  - `SSYSOBS=TOPOCENT`
- 観測時点で LSRK 追尾済みなら  
  - `CTYPE1=FREQ`
  - `SPECSYS=LSRK`
  - `SSYSOBS=LSRK`

### Calibration
- 強度較正のみ
- 軸は原則そのまま
- TOPOCENT 入力なら `VELOSYS` / `VFRAME` を付けてよい
- LSRK 入力なら `VELOSYS`, `VFRAME` を新規付与しない

### Coadd
- 入力が TOPOCENT なら `VELOSYS` を使って LSRK にそろえる
- 入力が LSRK なら追加補正しない
- 出力は必ず **LSRK**

---

## 4.3 実務上の注意

### TOPOCENT 入力
次のどちらかが必要です。

- `VELOSYS` / `VFRAME` がすでにある
- あるいは、時刻・座標・サイト情報があって `compute_vcorr_series()` で再計算できる

これがないと coadd は止まります。

### LSRK 入力
`VELOSYS` や `VFRAME` に **非ゼロの未適用補正** が残っているとエラーになります。  
つまり、「LSRK なのに未適用補正列が残っている」状態を不整合として扱います。

---

## 5. 重要仕様 2: rest frequency の上書き

## 5.1 FREQ 軸のとき

`CTYPE1="FREQ"` のデータに対し、解析上だけ別の rest frequency で速度解釈したい場合は、

- 周波数 WCS (`CRVAL1`, `CDELT1`, `CRPIX1`) はそのまま
- `RESTFRQ` / `RESTFREQ` を新しい線の値に置き換える
- 速度計算は新しい rest frequency で行う

という扱いが基本です。

これは obs_method にあるような「同じ raw から 13CO と C18O を別々に処理する」用途に非常に重要です。

---

## 5.2 VRAD 軸のとき

`CTYPE1="VRAD"` の場合は、rest frequency を変えると **速度軸そのものの意味** が変わるため、`RESTFRQ` だけ変えるのは不整合です。

`restfreq.apply_restfreq_override()` は VRAD 軸に対して、厳密なアフィン変換

- `a = rest1 / rest2`
- `b = c * (1 - rest1 / rest2)`
- `v2 = a * v1 + b`

を使い、

- `CRVAL1`
- `CDELT1`
- table 内の対応列

を更新します。

### まとめ
- **FREQ 軸**: WCS を変えず rest だけ更新
- **VRAD 軸**: WCS も更新

---

## 6. 重要仕様 3: TEMPSCAL と BEAMEFF

## 6.1 基本方針

このパッケージは 2026 方針として、**読み込み時に温度スケールを勝手に変換しません**。

### つまり
- ファイルに `TEMPSCAL=TR*` と書いてあれば、TR* として保持
- 自動的に TA* に戻さない
- viewer や書き出し時の変換は、必要なときだけ明示的に行う

## 6.2 BEAMEFF の意味

- Ta* → Tr* 変換には `BEAMEFF` が必要
- `TR* = TA* / BEAMEFF`
- `TA* = TR* * BEAMEFF`

## 6.3 混在グループへの対処

coadd では同じグループ内で `BEAMEFF` が混在している場合があります。

このとき

- `normalize_if_mixed="auto"`  
  混在を検出したグループを TR* に正規化してまとめる
- `normalize_if_mixed="never"`  
  危険だがそのまま進める

という挙動です。

---

## 7. 重要仕様 4: big-endian と pandas の問題

FITS 由来の配列や table が big-endian のままだと、pandas の `iloc`, `query`, `isin`, 並び替えなどで落ちることがあります。  
このパッケージは `read_scantable()`, `find_scans()`, `filter_scantable()`, `merge_scantables()` などで **native endian へ正規化** する安全策を取っています。

### 実務上の意味
- FITS を読んだ直後の table をそのまま pandas でごりごり処理しても大丈夫な方向へ寄せている
- ただし、自前で別の FITS 読み込みをした場合は、同様の endian 正規化を意識した方がよい

---

# 8. モジュール別の詳細説明

---

## 8.1 `fitsio.py`
Scantable の標準 I/O です。

### 主な役割
- Scantable の読み書き
- SINGLE DISH BinTable 形式の SDFITS 読み込み
- VLA 対応
- `TEMPSCAL`, `BEAMEFF`, `RESTFRQ`, `VELOSYS` などの標準列の取り回し

### 公開 API
```python
def read_scantable(path: str, *, tr_input_policy: str = "preserve") -> Scantable:
def write_scantable(
    path: str,
    scantable: Scantable,
    spectrum_column: str = "DATA",
    overwrite: bool = True,
    **kwargs
) -> None:
def read_tastar_fits(path: str):
def write_sdfits(
    out_path: str,
    meta: dict,
    data: Union[np.ndarray, List[np.ndarray]],
    table: pd.DataFrame,
    history: dict | None = None,
    *,
    spectrum_column: str = "DATA",
    include_flag: bool = True,
    overwrite: bool = True,
    string_widths: dict[str, int] | None = None,
    **kwargs,
) -> None:
```

### `read_scantable()`
#### 何をするか
- FITS を読んで `Scantable` を返す
- `V_CORR_KMS` などの旧列名を `VELOSYS` / `VFRAME` へ移行
- `TAMB_K` → `THOT`, `TAU` → `TAU0` などの旧名も補正
- `TEMPSCAL`, `BEAMEFF` 列を保証
- big-endian table を native endian 化

#### 重要
- `tr_input_policy` は現状ほぼ `preserve` 的に振る舞い、TR* を勝手に TA* に変えない
- `VFRAME` のみあれば `VELOSYS` にミラーする

### `write_scantable()`
#### 何をするか
- `Scantable` をそのまま SDFITS へ保存する高水準ラッパ
- `meta` にだけある重要キーワードを table 列へ昇格させる

#### 主な自動昇格列
- `RESTFRQ`, `RESTFREQ`
- `CRVAL1`, `CDELT1`, `CRPIX1`
- `CTYPE1`, `CUNIT1`
- `SPECSYS`, `SSYSOBS`
- `VELDEF`, `VELOSYS`, `VFRAME`
- `TEMPSCAL`, `BEAMEFF`

### `write_sdfits()`
#### 何をするか
- 実際の SDFITS 書き出し
- 固定長配列と VLA の両対応
- `tempscal` / `data_scale` の指定に応じて、**書き出し時だけ** TA* ↔ TR* 変換できる

#### 主要引数
- `out_path`: 出力先
- `meta`: ヘッダ相当
- `data`: 2D ndarray または `List[np.ndarray]`
- `table`: 行ごとの情報
- `history`: 履歴
- `spectrum_column`: スペクトル列名
- `include_flag`: FLAG 列を出すか
- `overwrite`: 上書き
- `string_widths`: 文字列列幅の指定

#### `kwargs` で実質使うもの
- `tempscal` / `out_scale`: 書き出し時の温度スケール
- `data_scale`: 入力データが今どのスケールか
- `beameff_on_fail`: TR* 変換時に BEAMEFF 欠損をどう扱うか

#### 重要
- データを破壊せず、**保存時だけ変換**
- `TEMPSCAL` は on-disk 宣言として table と meta の両方へ反映

---

## 8.2 `rawspec.py`
HOT/OFF/ON の raw をまとめる層です。

### 公開 API
```python
def build_rawspec(
    *,
    hot: pd.DataFrame,
    on: pd.DataFrame,
    off: pd.DataFrame,
    meta: dict[str, Any] | None,
    mapping: pd.DataFrame | None,
) -> RawSpec:
def save_rawspec(raw: RawSpec, path: str) -> None:
def load_rawspec(path: str) -> RawSpec:
def load_rawspec_auto(path: str, prefer: tuple[str, ...] | None = None) -> RawSpec:
def load_rawspec_fits(path: str) -> RawSpec:
```

### `build_rawspec()`
#### 何をするか
- HOT / OFF / ON の DataFrame から RawSpec を組み立てる
- `mapping` が ON のみなら、HOT / OFF 行を自動生成して全体 table を作る
- 各 DataFrame は `DatetimeIndex` を持つ必要がある

#### 引数の意味
- `hot`, `on`, `off`: 2D DataFrame。行が dump、列が channel
- `meta`: 全体メタ
- `mapping`: ON dump 用の補助 table。`OBSMODE` がなければ ON 用と解釈

#### 注意
- `nchan` 整合性を厳しく確認
- index が naive なら UTC 扱いに直す

### `load_rawspec_auto()`
- FITS か pickle かを自動判定して読み込む
- CLI 的な入口として使いやすい

### `load_rawspec_fits()`
- 主対応は `OBSMODE` と `DATA` / `SPECTRUM` を持つ SDFITS BinTable
- fallback として `HOT`, `ON`, `OFF` の ImageHDU も読む

---

## 8.3 `atmosphere.py`
ATM 補正用の補助関数です。

### 公開 API
```python
def get_airmass(elevation_deg)
def estimate_t_atm(t_surface_k, model="offset", delta_t=15.0, eta=0.95)
def compute_t_cal_array(t_hot_k, t_atm_k, tau_zenith, elevation_deg, t_bg_k=...)
def extract_meta_value(mapping_df, meta_dict, valid_keys)
def extract_meta_array(mapping_df, meta_dict, valid_keys)
```

### 何をしているか
- 仰角からエアマスを計算
- 地表気温から空の有効温度 `T_atm` を推定
- Ulich & Haas 型の式で `T_cal` を計算
- `TAU0`, `WXTEMP`, `TAMBIENT` などの **別名を許容しながら** 値を取得

### 重要仕様
- 仰角 5 度未満は `sec(z)` が暴れないようにクリップ
- `extract_meta_array()` は table の時変列を優先し、一部 NaN は補間する

---

## 8.4 `calibrate.py`
Ta* キャリブレーション本体です。

### 公開 API
```python
def make_tastar_dumps(...)
def tastar_from_rawspec(...)
def run_tastar_calibration(...)
def recalibrate_tastar(...)
```

---

### `make_tastar_dumps()`
#### 何をするか
HOT / OFF / ON の raw から ON dump ごとの Ta* を作り、Scantable を返します。

#### 入力
- `raw["meta"]`
- `raw["on"]`, `raw["off"]`, `raw["hot"]`
- `raw["mapping"]`

#### 主な内部処理
1. `mapping` の時刻を確定
2. 行選択 `rows` / `exclude_rows`
3. `OBSMODE` ごとに時刻列を切る
4. ON 行を主対象に `mapping_on` を作る
5. 必要なら座標列を正規化
6. `SPECSYS=TOPOCENT` なら `VELOSYS` を確保または再計算
7. `ch_range` または `vlsrk_range_kms` でチャンネル範囲を決定
8. Basic 1-Temp または ATM で `T_cal` を作る
9. HOT/OFF を scan ごとに平均
10. OFF を ON 時刻へ内挿
11. gain を計算して Ta* を得る
12. `TCAL`, `CALSTAT`, `THOT`, `TAU0`, `VELOSYS` などを table に詰める

#### 全引数の説明
- `rows`, `exclude_rows`  
  対象 dump の選択。両方同時指定は禁止。
- `t_hot_k`  
  ホットロード温度 [K]。未指定時は metadata / mapping から探す。
- `tau_zenith`  
  `None` なら Basic 1-Temp、数値または `"auto"` なら ATM。
- `t_surface_k`  
  外気温 [K]。未指定時は `WXTEMP`, `TAMBIENT`, `TAMB`, `T_SURF` を探す。
- `t_atm_model`  
  `"offset"` または `"ratio"`。
- `t_atm_delta_k`  
  offset モデルの引き算量。
- `t_atm_eta`  
  ratio モデルの係数。
- `gain_mode`  
  `"hybrid"` または `"independent"`。
- `verbose`  
  ログ出力。
- `ch_range`  
  `(ch_start, ch_stop)`。0-based, stop は exclusive。
- `vlsrk_range_kms`  
  指定速度範囲をカバーする channel slice を自動決定。
- `v_corr_col`  
  未適用速度補正列の優先名。通常は `VELOSYS`。
- `coord_frame`  
  座標列解釈の指定。自動判定に任せることも可能。
- `vcorr_chunk_sec`  
  速度補正を疎サンプリングして補間するときの間隔。長い時系列で有効。
- `dtype`  
  出力データ型。
- `rest_freq`  
  解析上使う静止周波数上書き。

#### `gain_mode` の意味
##### `hybrid`（推奨）
- HOT 取得時刻で OFF を引いて **クリーンな分母** を作る
- その分母を ON 時刻へ内挿して使う
- 定在波や OFF 混入を減らしやすい

##### `independent`
- HOT と OFF をそれぞれ ON 時刻へ内挿して直接分母を作る
- 従来型に近い

#### ATM 発動条件
- `tau_zenith is None`  
  Basic 1-Temp
- `tau_zenith = 数値`  
  指定値で ATM
- `tau_zenith = "auto"` または `"header"`  
  `TAU0` 系列を table / meta から検索

#### 出力 table に付く重要列
- `THOT`
- `TAMB_K`（後方互換）
- `TCAL`
- `CALSTAT`
- `TAU0`（ATM 時）
- `CH_START`, `CH_STOP`, `NCHAN_SEL`
- `TEMPSCAL="TA*"`
- `CRVAL1`, `CDELT1`, `CRPIX1`
- `CTYPE1="FREQ"`
- `CUNIT1="Hz"`
- `RESTFREQ`, `RESTFRQ`
- `SPECSYS`, `SSYSOBS`
- `VELDEF="RADIO"`
- TOPOCENT 入力なら `VELOSYS`, `VFRAME`

---

### `run_tastar_calibration()`
#### 何をするか
`make_tastar_dumps()` の高水準ラッパです。  
入力として

- パス
- Scantable
- raw dict

を受け、必要ならそのまま書き出します。

#### 追加引数
- `input_data`
- `output_path`
- `spectrum_column`
- `overwrite`
- `store_freq_column`

#### 注意
`store_freq_column` は `write_scantable()` へ透過的に渡されます。  
SDFITS writer 側では、LSRK や時間依存軸を厳密に保存したいときに `FREQ` ベクトル列を持たせる文脈があります。

---

### `recalibrate_tastar()`
#### 何をするか
すでに Ta* になっている Scantable に対し、RAW を読み直さず

- ATM 補正のやり直し
- Basic 1-Temp への巻き戻し

を行います。

#### 原理
既存データが

- `TA* = TCAL * CALRATIO`

という比例形で保存されているとみなし、

- `TA*_new = TA*_old * (TCAL_new / TCAL_old)`

で再スケールします。

#### 引数
- `scantable`
- `new_tau`
- `new_t_surface_k`
- `new_t_atm_model`
- `new_t_atm_delta_k`
- `new_t_atm_eta`
- `verbose`

#### 挙動
- `new_tau=None` なら Basic へ戻す
- それ以外なら新しい ATM パラメータで更新
- `TCAL` を差し替え
- `CALSTAT` を更新
- `TAU0` を更新または削除

#### 注意
- `THOT` と `TCAL` が table に必要
- 過去データで `THOT` がなければ `TAMB_K` にフォールバック
- 摂氏っぽい値（100 未満）は K に補正する安全策あり

---

## 8.5 `baseline.py`
ベースラインフィットと差し引きです。

### 公開 API
```python
@dataclass
class BaselineModel:
    poly_order: int
    v_windows_kms: List[Tuple[float, float]]
    convention: str = "radio"
    frame: str = "LSRK"

def fit_polynomial_baseline(...)
def run_baseline_fit(...)
```

### `fit_polynomial_baseline()`
#### 何をするか
1 本のスペクトルに多項式をフィットし、ベースラインを返します。

#### 引数
- `v_kms`  
  速度軸
- `y`  
  スペクトル
- `base_windows`  
  ベースラインに使う窓
- `windows`  
  旧名エイリアス
- `line_windows`  
  除外すべき線窓
- `poly_order` / `order`  
  多項式次数
- `iter_max`
- `iter_sigma`

#### 返り値
- `coeffs`
- `baseline`
- `info`

#### `info` の中身
- `mask`
- `rms`
- `std`
- `poly_order`
- `iter_max`

#### 実装上の特徴
- `np.polynomial.Polynomial.fit()` を使用
- iterative sigma-clipping 対応
- 窓未指定時は「安全側」としてゼロベースライン返却
- 明示的に空窓のときは全域使用に近い挙動

---

### `run_baseline_fit()`
#### 何をするか
Scantable の各行に対してベースラインを引き、残差を新しい data として返します。

#### 全引数の説明
- `input_data`  
  パスまたは Scantable
- `output_path`  
  保存先
- `rows`, `exclude_rows`  
  対象行選択
- `vwin`  
  ベースライン窓。必須。
- `poly_order`
- `line_vwin`
- `iter_max`
- `iter_sigma`
- `max_dumps`
- `ch_start`, `ch_stop`
- `v_corr_col`
- `rest_freq`
- `on_fail`  
  `exit`, `warn`, `skip` の想定
- `overwrite`

#### 重要挙動
- `rest_freq` があれば、まず `apply_restfreq_override()` を適用
- 行ごとに `SPECSYS` と `VELOSYS` の整合をチェック
- TOPOCENT なのに速度補正列がないと停止
- LSR 系なのに非ゼロ `VELOSYS` があると停止
- channel slice をした場合は WCS も更新
- 出力 table に baseline 情報を書き込む

#### 出力で増える主な列
- `BSL_DONE`
- `BSL_POLY`
- `BSL_WINF`
- `BSL_RMS`
- `BSL_COEF`
- `BSL_SCALE`

`BSL_SCALE` は、現在の `TEMPSCAL` と同じ意味で保存されます。  
つまり「この baseline は TA* で引いたのか TR* で引いたのか」を追跡できます。

#### 温度スケール方針
- TR* 入力を読んでも、baseline 前に TA* へ戻さない
- **入力されたスケールのままフィット** する
- viewer 側で表示切替が可能

---

## 8.6 `regrid_vlsrk.py`
異なる WCS を共通速度グリッドへそろえるモジュールです。  
**Standardizer が中核** です。

### 公開 API
```python
@dataclass
class VGrid:
    v0_kms: float
    dv_kms: float
    nchan: int
    crpix1: float = 1.0

def make_vgrid(vmin_kms: float, vmax_kms: float, dv_kms: float) -> VGrid:
def get_axis_signature(row: pd.Series, meta: dict, nchan: int) -> tuple:
class Standardizer:
    def __init__(self, scantable, target_grid: Optional[VGrid] = None, v_corr_col: str = "VELOSYS"):
    def auto_determine_grid(self, dv=None, vmin=None, vmax=None, margin=0.0)
    def get_matrix(self, target_grid=None, dv=None, vmin=None, vmax=None)
def vlsrk_axis_for_spectrum(meta: dict, *, v_corr_kms: float, nchan: Optional[int] = None) -> np.ndarray:
def interp_to_vgrid(v_src: np.ndarray, y_src: np.ndarray, v_tgt: np.ndarray) -> np.ndarray:
def vrange_from_meta_and_vcorr(meta: dict, v_corr_kms: np.ndarray, *, nchan: int | None = None) -> tuple[float, float]:
```

---

### `VGrid`
共通速度グリッドを表す単純なデータクラスです。

- `v0_kms`: 先頭チャンネルの速度
- `dv_kms`: チャンネル幅
- `nchan`: チャンネル数
- `crpix1`: 1-based reference pixel

`axis()` で実際の速度配列を返します。

---

### `make_vgrid()`
#### 何をするか
`vmin`, `vmax`, `dv` から `VGrid` を作ります。

#### 注意
- `dv` は正値でなければならない
- `dv > (vmax-vmin)` なら 1 チャンネルのグリッドを返す

---

### `get_axis_signature()`
#### 何をするか
各行について、

- `nchan`
- `CRVAL1`
- `CDELT1`
- `CRPIX1`
- `CTYPE1`
- `CUNIT1`
- `RESTFREQ`

などから「この行の軸は同型か」を判定する署名を作ります。

#### なぜ必要か
全行ごとに軸を毎回ゼロから作ると遅いので、同じ署名の行はまとめて扱います。

---

### `Standardizer`
#### 何をするクラスか
不均質な Scantable を、**共通の LSRK 速度グリッドを持つ 2 次元行列**へ変換します。

#### 内部の重要思想
- 行ごとに WCS が少し違っていてもよい
- ただし同じ `AxisSignature` を持つ行はグループ化する
- 各グループで基底軸を 1 回だけ作り、キャッシュする
- TOPOCENT 行なら `VELOSYS` を使って行ごとに速度軸をずらす
- LSRK 行なら `v_corr=0` 扱いでそのまま使う

---

### `Standardizer.auto_determine_grid()`
#### 何をするか
Scantable 全体を覆うマスター速度グリッドを自動決定します。

#### 引数
- `dv`  
  明示指定すればそれを採用。`None` なら最小分解能を採用。
- `vmin`, `vmax`  
  明示指定がなければ全行を覆う範囲を自動推定。
- `margin`  
  自動推定範囲に付け足す余裕 [km/s]

#### 決め方
- すべての Signature グループを走査
- 各グループの観測周波数軸から速度分解能を推定
- TOPOCENT なら `v_corr` の最小最大も加味
- 全体を覆う `vmin`, `vmax`, `dv` を決める

---

### `Standardizer.get_matrix()`
#### 何をするか
全行を target grid へ再サンプリングし、

- `(N_row, N_grid)` の行列
- `v_tgt`

を返します。

#### 処理内容
- target grid を確定
- Signature ごとに行を束ねる
- `FREQ` 軸なら基底周波数軸を作る
- `VRAD` 軸なら一度観測周波数へ戻す
- TOPOCENT なら `get_doppler_factor(v_corr)` で各行の LSRK 軸を作る
- `np.interp()` で target grid に載せる

#### 重要
coadd の前段で「入力が全部同じチャンネル軸である」と仮定せずに済むのは、この Standardizer のおかげです。

---

### `vlsrk_axis_for_spectrum()`
1 本のスペクトルについて、指定 `v_corr_kms` を適用した LSRK 速度軸を返します。

### `interp_to_vgrid()`
1 本のスペクトルを target velocity grid に補間します。

### `vrange_from_meta_and_vcorr()`
header と `v_corr` 範囲から、取り得る LSRK 速度域を見積もります。

---

## 8.7 `coadd.py`
速度再グリッドと加算の本体です。  
このパッケージで最も機能が多い部分です。

### 公開 API
```python
def run_velocity_coadd(...)
```

### QC プリセット
```python
QC_PRESETS = {
    "robust":   {"hard": 4.0, "soft_a": 1.0, "soft_p": 4.0, "clip_k": 3.0, "clip_iter": 3},
    "gentle":   {"hard": 0.0, "soft_a": 1.0, "soft_p": 3.0, "clip_k": 3.0, "clip_iter": 2},
    "hardonly": {"hard": 4.0, "soft_a": 0.0, "soft_p": 0.0, "clip_k": 0.0, "clip_iter": 0},
    "noclip":   {"hard": 4.0, "soft_a": 1.0, "soft_p": 4.0, "clip_k": 0.0, "clip_iter": 0},
}
```

---

### `run_velocity_coadd()`
#### 何をするか
複数入力の Scantable を

1. 必要なら TOPOCENT → LSRK 補正
2. Standardizer で共通速度グリッドへ再配置
3. グループ分け
4. 重み付き加算
5. 必要なら post-coadd baseline
6. LSRK 出力

まで行います。

---

### 全引数の説明

#### 入出力・選択
- `inputs`  
  パスまたは Scantable の列。文字列 1 個でも可。
- `output_path`
- `rows`, `exclude_rows`
- `overwrite`

#### グループ化
- `group_mode`  
  `"position"`, `"scan"`, `"intgrp"`
- `pos_col`  
  position grouping の ID 列。通常 `pos_id`
- `pos_tol_arcsec`  
  `pos_col` がないとき、RA/DEC 近傍で grouping する許容距離
- `coord_frame`  
  VELOSYS 再計算時の座標解釈

#### 速度補正
- `v_corr_col`  
  未適用補正列の優先名。通常 `VELOSYS`
- `vcorr_chunk_sec`  
  長い時系列で速度補正を間引き計算する間隔

#### 速度グリッド
- `vmin`, `vmax`, `dv`  
  target velocity grid 指定。未指定なら自動
- `allow_outside_overlap`  
  現状は説明上は「共通重なり外も許す」意図だが、コード理解としては最重要パラメータではない
- `axis_type`  
  出力軸。`"freq"` または `"vel"`
- `rest_freq`  
  解析時の静止周波数上書き

#### 重み付けモード
- `mode`
  - `"uniform"`: 一様重み
  - `"rms"`: RMS から重みを作る
  - `"rms_weight"`: 実質的には RMS ベースの重み付け coadd として使う
- `baseline_vwin`
  baseline subtraction と RMS 評価を同時にやりたい窓
- `baseline_poly`
- `baseline_iter_max`
- `baseline_iter_sigma`
- `rms_vwin`
  データ本体からは引かず、RMS 評価だけに使う窓
- `rms_poly`
- `rms_bin`
- `coadd_qc`
  QC preset / 指定文字列
- `line_vwin`
  baseline / RMS 窓から引き去る線窓

#### 入力量制限
- `block_size`
- `max_dumps`
- `ch_start`, `ch_stop`

#### エラーと温度スケール
- `on_fail`
- `out_scale`
  出力スケール。`"TA*"` または `"TR*"`
- `normalize_if_mixed`
  混在 BEAMEFF 群への対処
- `beameff_tol`
  混在判定の相対閾値
- `sigma_scale`
  現状 `"TA*"` のみ実装
- `verbose`

---

### 非常に重要な制約

#### `baseline_vwin` と `rms_vwin` は同時指定不可
コード上で明示的に禁止されています。  
理由は、「ベースラインを引きながらその窓で RMS を評価するモード」と「引かずに RMS だけ測るモード」を混ぜないためです。

#### `coadd_qc` は `mode="uniform"` と併用不可
QC は統計的重み付けを伴うので uniform と矛盾します。

---

### グループ化の意味

#### `group_mode="position"`
- `pos_col` があればそれでグループ化
- なければ `pos_tol_arcsec` を用いて RA/DEC 近傍でまとめる

#### `group_mode="scan"`
同じ `SCAN` 番号だけでなく、さらに

- `OBSMODE`
- `FDNUM`
- `IFNUM`
- `PLNUM`
- `_INPUT_ID`

まで見て分離します。  
つまり、**同じ scan 番号でも IF や偏波が違うものを勝手に足さない** ようになっています。

#### `group_mode="intgrp"`
`INTGRP` 列でまとめます。

---

### 重み付けの実際

#### `mode="uniform"`
全スペクトル重み 1

#### `mode="rms"` / `"rms_weight"`
ケース分けは次の通りです。

##### `baseline_vwin` がある
- 各スペクトルごとに baseline fit
- baseline を差し引いた残差から RMS を計算
- `1 / RMS^2` を重みにする

##### `rms_vwin` がある
- baseline はデータ本体に引かない
- RMS 評価用にだけフィットして `1 / RMS^2`

##### どちらもない
- 入力 table の `BSL_RMS` を使う
- これも `1 / RMS^2`

#### `rms_bin`
RMS 評価前に箱平均してノイズ推定を安定化する用途です。

---

### QC モード
`coadd_qc="robust"` などを与えると、

- 指定窓で事前に baseline を引いた残差を評価
- hard / soft の基準で重み付け
- 必要なら channel ごとの clipping

を行います。

文字列は例えば

- `robust`
- `gentle`
- `robust;clip=3,2`
- `robust;hard=4;soft=1,4`

のように指定できます。

---

### coadd 後の baseline
`baseline_windows_parsed` が指定されていて QC でない場合、最終 `out_spec` に対して **post-coadd baseline** をもう一度当てます。  
その結果は

- `BSL_DONE`
- `BSL_POLY`
- `BSL_WINF`
- `BSL_RMS`

として出力行へ記録されます。

---

### BEAMEFF 混在時の処理
グループ内 `BEAMEFF` が混在しているとき、

- `out_scale="TR*"` なら TR* で合算
- あるいは `normalize_if_mixed="auto"` なら自動的に TR* 正規化

になります。

このとき出力は
- `TEMPSCAL="TR*"`
- `BEAMEFF=1.0`

相当に寄せる思想です。

---

### 出力の基準系
coadd 出力は常に **LSRK** です。

#### `axis_type="vel"`
- `CTYPE1="VRAD"`
- `CUNIT1="m/s"`
- `SPECSYS="LSRK"`
- `SSYSOBS="LSRK"`
- `VELDEF="RADIO"`

#### `axis_type="freq"`
- `CTYPE1="FREQ"`
- `CUNIT1="Hz"`
- `SPECSYS="LSRK"`
- `SSYSOBS="LSRK"`

また、入力時の `VELOSYS` / `VFRAME` は
- `VELOSYS_OBS`
- `VFRAME_OBS`
のような監査列へ退避し、出力の `VELOSYS` / `VFRAME` 自体は 0 にします。

---

## 8.8 `restfreq.py`
rest frequency 上書き専用です。

### 公開 API
```python
def apply_restfreq_override(
    meta: dict,
    table: Optional[pd.DataFrame],
    rest_freq_hz: float,
    *,
    require_wcs_for_vrad: bool = True,
) -> Dict[str, Any]:
```

### 何をするか
- meta と table に対し、rest1 → rest2 変換を in-place で適用
- FREQ 軸なら WCS は維持
- VRAD 軸なら WCS も厳密更新

### 引数
- `meta`
- `table`
- `rest_freq_hz`
- `require_wcs_for_vrad`  
  VRAD 軸で `CRVAL1`, `CDELT1`, `CRPIX1` がないと困るので、厳密に要求するかどうか

### 戻り値
情報辞書で、例えば
- `rest1_hz`
- `rest2_hz`
- `axis_vrad`
- `wcs_updated`
- `a`, `b`
- `vel_unit`

などを返します。

---

## 8.9 `tempscale.py`
温度スケールと BEAMEFF のユーティリティです。

### 公開 API
```python
def normalize_tempscal(value, *, default="TA*") -> str
def ensure_tempscal_column(df, *, default="TA*") -> pd.DataFrame
def ensure_beameff_column(df, *, default=np.nan) -> pd.DataFrame
def beameff_array(df, meta, n_rows, *, default=np.nan) -> np.ndarray
def tempscal_array(df, meta, n_rows, *, default="TA*") -> np.ndarray
def require_beameff(beameff, *, on_fail="error") -> None
def is_beameff_mixed(beameff, *, tol=1e-6) -> bool
def representative_beameff(beameff) -> float
def ta_to_tr(ta, beameff) -> np.ndarray
def tr_to_ta(tr, beameff) -> np.ndarray
def convert_rowwise_vla(spectra, beameff, *, direction) -> List[np.ndarray]
def append_scale_history(history, entry) -> dict
def set_beameff(sc, efficiency, rows=None, verbose=True) -> None
```

### 重要関数

#### `normalize_tempscal()`
- `"TA*"`, `"TA"`, `"TR*"`, `"TMB"` などを正規化
- `"TMB"` 系もこの実装では `TR*` 側に寄せる

#### `beameff_array()`
- table に `BEAMEFF` があればそれを優先
- なければ meta の `BEAMEFF`
- それもなければ `default`

#### `tempscal_array()`
- table の `TEMPSCAL` → meta → default の順

#### `require_beameff()`
- TR* 変換時に BEAMEFF が有限正値かを厳密確認

#### `convert_rowwise_vla()`
- VLA データにも TA* ↔ TR* 変換を適用

#### `set_beameff()`
- データ自体は変えず、行ごとの `BEAMEFF` を設定するだけ
- viewer のオンザフライ変換や write 時変換の準備として使う

---

## 8.10 `scantable_utils.py`
Scantable の table, data, meta を安全に検査・抽出・補助更新するための実務モジュールです。  
このモジュールは、単なる表示用小物ではなく、実際には以下のような場面で非常に有用です。

- 読み込んだ SDFITS の table をその場で検査する
- `TIME`, `TIMESTAMP`, `RESTFRQ`, `TEMPSCAL`, `BEAMEFF` などの実在状況を確認する
- scan や観測モードごとに行を選別して試験解析する
- 日ごとのファイルを結合するときに `SCAN` 重複を避ける
- マッピング用の相対座標を計算する
- 行ごとに異なる `BEAMEFF` を安全に与える
- 重要列をうっかり破壊しないように制限付きで metadata を更新する

この節では、ソースの実装に基づいて、各関数の役割・引数・戻り値・制約・典型的な使い方を詳しく説明します。

### 公開 API
```python
def describe_columns(sc: Scantable) -> None

def show_scantable(
    inputs: Union[Scantable, str, Sequence[Union[Scantable, str]]],
    rows: Union[str, slice, List[int], None] = None,
    columns: Union[List[str], str] = "default",
    head: Optional[int] = 20,
    show_legend: bool = False,
    extra_data: Optional[pd.DataFrame] = None,
    ref_coord: Optional[Union[str, SkyCoord]] = None,
    frame: str = "ICRS",
    projection: str = "GLS",
    unit: str = "arcsec",
) -> None

def calc_mapping_offsets(
    sc: Scantable,
    ref_coord: Optional[Union[str, SkyCoord]] = None,
    frame: str = "ICRS",
    projection: str = "GLS",
    unit: str = "arcsec",
    cos_mode: str = "point",
    verbose: bool = True,
) -> pd.DataFrame

def merge_scantables(
    inputs: Sequence[Union[Scantable, str]],
    sort_by_time: bool = False,
    shift_scan_id: bool = True,
) -> Scantable

def update_metadata(
    sc: Scantable,
    column: str,
    value: Any,
    rows: Union[str, List[int], None] = None,
    force: bool = False,
    verbose: bool = True,
) -> None

def set_beameff(
    sc: Scantable,
    efficiency: Union[float, np.ndarray, Sequence[float]],
    rows: Union[str, List[int], None] = None,
    verbose: bool = True,
) -> None

def find_scans(
    sc: Scantable,
    query: Optional[str] = None,
    extra_data: Optional[pd.DataFrame] = None,
    **kwargs
) -> np.ndarray

def filter_scantable(
    sc: Scantable,
    query: Optional[str] = None,
    extra_data: Optional[pd.DataFrame] = None,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    **kwargs
) -> Scantable

def select_rows(sc: Scantable, rows: Union[str, slice, Sequence[int], int]) -> Scantable
```

### まず押さえるべき設計思想

#### table / data / meta の役割分担
- `sc.table` は各積分・各ダンプごとの表データです。
- `sc.data` は実スペクトル配列です。
- `sc.meta` はグローバルなヘッダ相当情報です。

`scantable_utils.py` は主に table と meta を扱います。  
ただし `filter_scantable()`, `merge_scantables()`, `select_rows()` は `data` 側も対応する行だけ追随させます。

#### 非破壊か破壊的か
- 非破壊: `describe_columns()`, `show_scantable()`, `calc_mapping_offsets()`, `find_scans()`, `filter_scantable()`, `select_rows()`, `merge_scantables()`
- 破壊的: `update_metadata()`, `set_beameff()`

#### 行指定 `rows` は positional index
多くの関数で使う `rows` は、table の表示 index ではなく **0-based positional index** です。  
負の index は想定されていません。

受理される主な形式:
- `None` : 全行
- `3` : 3 行目のみ
- `slice(0, 10)`
- `[0, 5, 10]`
- `"0:10"`
- `"0:10,50,60:80"`

`"0:10"` は Python slice と同様に stop を含みません。

#### Big-endian FITS への対策
本モジュールには `_df_to_native_endian()` が組み込まれており、FITS 由来の big-endian buffer が Pandas/Numpy の処理で落ちる問題を避けています。  
特に `show_scantable()`, `find_scans()`, `filter_scantable()`, `merge_scantables()` はこの保護の恩恵を受けます。

#### TIMESTAMP の自動解決
`show_scantable()` や `merge_scantables(sort_by_time=True)` では、必要に応じて内部的に時刻列を解決します。優先順位は次の通りです。

1. `TIMESTAMP`
2. `MJD`
3. `DATE-OBS` / `DATEOBS` と `TIME`
4. 旧来の `TIME`
5. DatetimeIndex

したがって、table に `TIMESTAMP` が無くても、表示時には先頭に補われることがあります。

---

### `describe_columns()`
#### 何をするか
Scantable に含まれる table 列の一覧と dtype を表示し、続けて `meta` の全キーを一覧表示します。

#### 使いどころ
- まず何が入っているかを把握したいとき
- `BEAMEFF`, `TEMPSCAL`, `VELOSYS`, `RESTFRQ`, `FDNUM`, `IFNUM`, `PLNUM` の実在確認
- viewer や coadd の前に列の存在を確認したいとき

#### 注意点
- 説明文は `DEFAULT_SHOW_COLS` に登録された列だけに付くため、すべての列に説明が出るわけではありません。
- ただし列名と dtype を俯瞰するには十分有用です。

---

### `show_scantable()`
#### 何をするか
Scantable の table を人間が確認しやすい形で印字します。  
単なる `print(sc.table.head())` よりかなり実務的で、以下の補助機能を持ちます。

- ファイルパスを直接渡せる
- 複数入力を順に表示できる
- `TIMESTAMP` を best-effort で補う
- table に無い列を meta やエイリアスから拾う
- その場で offset を計算して追加できる
- 外部の DataFrame を仮結合できる

#### 入力 `inputs`
- `Scantable`
- ファイルパス文字列
- それらのシーケンス

#### `columns`
- `"default"` : `SCAN`, `TIMESTAMP`, `OBSMODE`, `OBJECT`, `OFS_LON`, `OFS_LAT`, `EXPOSURE`, `TEMPSCAL`
- `"all"` : table 上の全列
- `"SCAN, TIME, OBSMODE"` のようなカンマ区切り文字列
- `['SCAN', 'TIME']` のような list

#### table に無い列をどう探すか
たとえば `columns="RESTFREQ"` を指定したのに table に `RESTFREQ` が無い場合、次の順で探します。

1. table 内の同名列
2. table 内の alias 列 (`RESTFRQ` など)
3. meta の同名キー
4. meta の alias キー

このため、`RESTFRQ` と `RESTFREQ` が環境によって揺れていても確認しやすくなっています。

#### `rows`
`rows` は位置ベース選択です。  
`"0:5,10,20:25"` のような複合指定もできます。

#### `extra_data`
外部 DataFrame を一時的に横結合して表示します。  
行数が一致しない場合は無視されます。

有用な使い方:
- `calc_mapping_offsets()` の結果を足す
- 自前で計算した RMS, OFFS, QC 指標を足す
- basket-weave 補正量や外部フラグを足す

#### `ref_coord`, `frame`, `projection`, `unit`
これらを与えると内部で `calc_mapping_offsets()` を呼び、`OFS_LON`, `OFS_LAT` を追加表示します。

#### 重要な実装上の注意
1. **天体名の名前解決はしていません。**  
   `ref_coord="Orion KL"` のような名前文字列を `SkyCoord.from_name()` で解決する実装ではなく、単に `SkyCoord(ref_coord)` へ渡しています。  
   したがって、確実なのは以下です。
   - 座標文字列を渡す: `"05h35m14.5s -05d22m30s"`
   - `SkyCoord` オブジェクトを渡す

2. `show_scantable()` 自体は data 配列を変更しません。

3. 長大表の全行印字は見づらいので、通常は `head` を使う方が安全です。

#### 実用例
```python
import sd_radio_spectral_fits.scantable_utils as su

su.show_scantable(
    inputs="tastar_dump.fits",
    rows="0:10",
    columns="SCAN, TIMESTAMP, OBSMODE, RESTFRQ, RA, DEC, VELOSYS, TEMPSCAL, BEAMEFF",
    head=20,
)
```

#### オフセットをその場で見たい場合
```python
su.show_scantable(
    inputs=sc,
    rows="0:20",
    columns="TIMESTAMP, RA, DEC, OFS_LON, OFS_LAT",
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
)
```

#### 外部計算列を足して見る
```python
ofs = su.calc_mapping_offsets(
    sc,
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
    cos_mode="ref",
)

su.show_scantable(
    sc,
    rows="0:10",
    columns="SCAN, OFS_LON, OFS_LAT, OBSMODE",
    extra_data=ofs,
)
```

---

### `calc_mapping_offsets()`
#### 何をするか
Scantable の座標列から、参照点に対するマッピングオフセット `OFS_LON`, `OFS_LAT` を計算し、**新しい DataFrame として返します。**  
Scantable 本体は変更しません。

#### 入力座標として認識する列
優先順は次の通りです。

1. `RA`, `DEC`
2. `GLON`, `GLAT`
3. `AZIMUTH`, `ELEVATIO`
4. `AZIMUTH`, `ELEVATION`
5. `AZ`, `EL`

#### `frame`
内部では以下の補助マッピングがあります。
- `j2000` -> `fk5`
- `b1950` -> `fk4`
- `lb` -> `galactic`

それ以外は小文字化して Astropy の frame 名へ渡します。

#### `ref_coord`
- `None` の場合
  - source frame が ICRS 系で、`meta` に `OBSRA`, `OBSDEC` があればそれを優先
  - それが無ければ table の円周平均・平均値から自動設定
- 座標文字列または `SkyCoord` を与えた場合はそれを使用

#### `projection`
現コードで実際にサポートされるのは次です。

- `GLS` または `SFL`
  - `ofs_x = dlon * cos(lat)`
  - `ofs_y = dlat`
- `CAR` または `NONE`
  - `ofs_x = dlon`
  - `ofs_y = dlat`

#### 重要: ドキュメント文字列との差異
関数 docstring には `TAN` や `SIN` に言及がありますが、**現コードでは `TAN` は実装されていません。**  
またエラーメッセージには `SIN` が含まれる一方、分岐実装は `GLS/SFL` と `CAR/NONE` です。  
したがって、現時点で確実に使うべきなのは `GLS` か `CAR` です。

#### `cos_mode`
`GLS/SFL` でだけ意味があります。

- `"point"` : 各点の緯度で `cos(lat)` を評価
- `"ref"` : 参照点の緯度 `lat0` で評価

小領域マップで基準点中心の幾何を強く意識するなら `ref` が分かりやすく、従来互換を重視するなら `point` です。

#### `unit`
- `arcsec`
- `arcmin`
- `deg`

#### AltAz データの制限
source が AltAz 系のとき、現コードは **他 frame への変換を行いません。**  
site/time 情報が必要になるためで、`frame="altaz"` のように同じ frame 内での差分計算だけが許容されます。

#### 返り値
`pd.DataFrame`。列は必ず
- `OFS_LON`
- `OFS_LAT`
です。

#### 実用例
```python
ofs = su.calc_mapping_offsets(
    sc,
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
    cos_mode="ref",
)
```

#### どんなときに有用か
- マップ中心からのオフセットを table レベルで確認したい
- `OFS_LON`, `OFS_LAT` を使って `find_scans()` で位置抽出したい
- viewer や独自グリッド処理の前段確認をしたい

たとえば中心から半径 60 arcsec 以内だけ抜きたいなら、

```python
ofs = su.calc_mapping_offsets(sc, ref_coord="05h35m14.5s -05d22m30s")
idx = su.find_scans(sc, query="OFS_LON**2 + OFS_LAT**2 < 60**2", extra_data=ofs)
sc_center = su.filter_scantable(sc, rows=idx)
```

---

### `merge_scantables()`
#### 何をするか
複数の Scantable またはファイルパスを結合し、単一の Scantable を返します。  
単に table を縦結合するだけではなく、`data` 配列も同順に結合します。

#### 引数
- `inputs` : Scantable またはパスの列
- `sort_by_time=False` : 必要なら時刻順に並べ替える
- `shift_scan_id=True` : 結合時に `SCAN` 重複を避ける

#### `shift_scan_id=True` が重要な理由
日別ファイルや分割ファイルを単純結合すると、各ファイルで `SCAN=0,1,2...` が再利用されていることがあります。  
この状態で後段の coadd を `group_mode="scan"` で動かすと、**別ファイルの scan が同一 scan と誤認される** 危険があります。

この関数は `shift_scan_id=True` のとき、各入力の `SCAN` に累積オフセットを足し、番号衝突を避けます。  
実運用では、この既定値を維持するのが安全です。

#### 時刻ソートの方法
`sort_by_time=True` のとき、内部で `_resolve_table_timestamps()` を使って時刻列を解決し、table と data を同じ順序で並べ替えます。

#### meta / history の扱い
- `meta` は先頭入力のコピーを使います
- `history` には `merged_count` と `shifted_scan_id` が入ります
- 複数 meta をマージしてはいません

つまり、ヘッダ値が入力間で異なる場合、代表値は先頭ファイル依存です。  
重要な meta が混在しているときは、結合後に `show_scantable()` や `update_metadata()` で点検した方が安全です。

#### 実用例
```python
sc_merged = su.merge_scantables(
    inputs=["day1.fits", "day2.fits", "day3.fits"],
    sort_by_time=True,
    shift_scan_id=True,
)
```

#### 特に有効な場面
- 日をまたいだ反復観測の結合
- coadd 前の試験的な multi-file 結合
- profile viewer へ複数ファイルをまとめて渡したいとき

---

### `update_metadata()`
#### 何をするか
指定した列の値を、指定した行に対して **in-place** で更新します。  
便利ですが危険でもあるため、安全装置付きです。

#### 主な引数
- `column` : 更新対象列名
- `value` : スカラーまたは配列
- `rows` : 対象行
- `force` : 安全装置を無効化するか
- `verbose` : ログを出すか

#### 内蔵されている安全装置
##### 1. 危険列のブロック
以下は `force=False` では変更をブロックします。

- `CRVAL1`, `CDELT1`, `CRPIX1`, `CTYPE1`, `CUNIT1`
- `RESTFREQ`, `RESTFRQ`
- `SPECSYS`
- `DATA`, `FLAG`, `SPECTRUM`

これらはスペクトル軸やデータ本体に直結するためです。

##### 2. enum 値の検証
以下の列では、既知の候補以外は `force=False` でブロックされます。

- `OBSMODE`
- `TEMPSCAL`
- `CALSTAT`

##### 3. `TEMPSCAL='TR*'` と `BEAMEFF` の整合性確認
`TEMPSCAL` を `TR*` にしようとした場合、table または meta に有効な `BEAMEFF` が無いと警告して止めます。  
`force=True` なら強行できますが、通常は避けるべきです。

#### alias 同期
`RESTFREQ` と `RESTFRQ` は相互 alias として扱われ、片方を更新すると table 内に存在するもう片方も同期更新されます。

#### 新規列の作成
指定列が無ければ新規作成します。  
したがって、補助列を手で入れる用途にも使えます。

#### 破壊的変更であることに注意
この関数は元の `sc.table` を直接変更します。  
試験的修正なら、事前に `filter_scantable()` や `select_rows()` で小さい Scantable を作ってから試す方が安全です。

#### 実用例: OBJECT の修正
```python
su.update_metadata(
    sc,
    column="OBJECT",
    value="Orion-KL",
    rows="0:100",
)
```

#### 実用例: RESTFRQ を強制修正
```python
su.update_metadata(
    sc,
    column="RESTFRQ",
    value=115.2712018e9,
    force=True,
)
```

#### どんなときに使うべきか
- `OBJECT`, `OBSMODE`, `CALSTAT` など table ラベルの是正
- 手作業で補助列を足す
- どうしても必要なヘッダ修正を、小規模検証の上で限定的に行う

WCS 本体を本格的に作り替える用途には、より専用の関数や再生成処理を優先すべきです。

---

### `set_beameff()`
#### 何をするか
指定した行の `BEAMEFF` 列へビーム効率を設定します。  
**データ配列そのものも、`TEMPSCAL` も変更しません。**  
あくまで、後段の viewer や scale 変換処理が参照するための metadata を埋める関数です。

#### 実装上の位置づけ
同名関数は `tempscale.py` にもあります。  
実装内容はほぼ同趣旨で、現在の `scantable_utils.py` 版も同様に以下を行います。

- `BEAMEFF` 列が無ければ作成
- 指定行へ値を代入
- 値が 0 以下または 1 より大きいと警告
- `TEMPSCAL` は変更しない

#### 重要: 何をしないか
この関数は次をしません。

- `TA* -> TR*` へデータ実体を変換しない
- `TEMPSCAL` を `TR*` に変えない
- `target_scale` のような引数で表示スケールを切り替えない

したがって、ユーザー例として見かける
```python
sd.set_beameff(sc, efficiency=0.45, target_scale="TR*", rows="0:100")
```
の `target_scale` は、**この添付ソースの `set_beameff()` には存在しません。**  
現在のコードに即して書くなら、まず `BEAMEFF` を付与し、その後に
- viewer で TR* 表示へ切り替える
- `coadd.run_velocity_coadd(..., out_scale="TR*")` を使う
- `write_scantable(..., tempscal="TR*")` 側で明示変換する
といった流れになります。

#### 典型的な使い方
##### 全行へ同じ効率を付与
```python
su.set_beameff(sc, efficiency=0.42)
```

##### 一部行だけに付与
```python
su.set_beameff(sc, efficiency=0.42, rows="0:20")
```

##### 配列で行ごとに別値を入れる
```python
rows = [0, 1, 2]
eta = [0.40, 0.41, 0.42]
su.set_beameff(sc, efficiency=eta, rows=rows)
```

#### FDNUM / IFNUM / PLNUM ごとに異なる `BEAMEFF` を入れる
実務上もっとも重要なのはここです。  
多ビーム・多 IF・多偏波のデータでは、`FDNUM`, `IFNUM`, `PLNUM` の組ごとに効率が異なることがあります。  
その場合、`find_scans()` と組み合わせて **組ごとに個別の `BEAMEFF` を入れる** のが自然です。

##### 例1: 特定の 1 組だけ更新
```python
idx = su.find_scans(sc, FDNUM=0, IFNUM=1, PLNUM=0)
su.set_beameff(sc, efficiency=0.45, rows=idx)
```

##### 例2: 複数組に別々の値を入れる
```python
beam_eff_map = {
    (0, 0, 0): 0.46,
    (0, 0, 1): 0.44,
    (1, 0, 0): 0.42,
    (1, 0, 1): 0.41,
}

for (fdnum, ifnum, plnum), eta in beam_eff_map.items():
    idx = su.find_scans(sc, FDNUM=fdnum, IFNUM=ifnum, PLNUM=plnum)
    su.set_beameff(sc, efficiency=eta, rows=idx, verbose=True)
```

##### 例3: query を使ってまとめて選ぶ
```python
idx = su.find_scans(
    sc,
    query="FDNUM == 0 and IFNUM == 1 and PLNUM in [0, 1]"
)
su.set_beameff(sc, efficiency=0.43, rows=idx)
```

#### 使ったあとに確認する
```python
su.show_scantable(
    sc,
    rows="0:20",
    columns="FDNUM, IFNUM, PLNUM, TEMPSCAL, BEAMEFF"
)
```

#### `TEMPSCAL` との関係
- `BEAMEFF` を付けただけでは data は依然として元スケールのままです。
- `TEMPSCAL='TA*'` のデータに `BEAMEFF` を入れておくと、viewer で TR* 表示へ切り替える準備ができます。
- `out_scale='TR*'` で coadd するとき、混在ビーム効率の扱い方が重要になります。詳細は `tempscale.py` と `coadd.py` の節を参照してください。

---

### `find_scans()`
#### 何をするか
条件に合致する行の **0-based positional index** を `np.ndarray` で返します。  
名称は `find_scans` ですが、scan 単位ではなく行レベル選択です。

#### 条件指定の方法
##### 1. `query`
Pandas の `DataFrame.query()` 構文です。

```python
idx = su.find_scans(sc, query="SCAN > 10 and OBSMODE == 'ON'")
```

##### 2. `kwargs`
等価比較です。

```python
idx = su.find_scans(sc, OBSMODE="OFF")
```

list / tuple / ndarray を渡すと `isin` 扱いになり、OR 検索になります。

```python
idx = su.find_scans(sc, OBSMODE=["OFF", "SKY"])
```

#### `extra_data`
行数が一致する外部 DataFrame を追加して、その列も検索に使えます。

```python
ofs = su.calc_mapping_offsets(sc, ref_coord="05h35m14.5s -05d22m30s")
idx = su.find_scans(
    sc,
    query="OFS_LON**2 + OFS_LAT**2 < 120**2",
    extra_data=ofs,
)
```

#### 列名の扱い
`kwargs` では
- まずそのままの列名
- 無ければ大文字化した列名
を探します。

#### エラー時
- 不正 query は `ValueError`
- 存在しない列の kwargs は `KeyError`

---

### `filter_scantable()`
#### 何をするか
`find_scans()` の結果に対応する行だけを取り出し、新しい Scantable を返します。  
`table` だけでなく `data` も対応行だけに切り詰められます。

#### 特徴
- 元の Scantable は変更しない
- `history` に `filter_query` を残す
- `rows` を併用すると、`find_scans()` 結果との **交差** を取る

#### `rows` 併用時の挙動
`rows` が指定されているときは、`rows` の順序を優先した交差になります。  
この仕様は、飛び飛び行を指定しつつ query 条件でも絞りたいときに便利です。

#### 実用例
```python
sc_on = su.filter_scantable(sc, OBSMODE="ON")
```

```python
sc_sub = su.filter_scantable(
    sc,
    query="SCAN >= 10 and SCAN <= 30",
    rows="0:500",
)
```

```python
ofs = su.calc_mapping_offsets(sc, ref_coord="05h35m14.5s -05d22m30s")
sc_center = su.filter_scantable(
    sc,
    query="OFS_LON**2 + OFS_LAT**2 < 60**2",
    extra_data=ofs,
)
```

#### どんなときに有用か
- 特定 scan 群だけで baseline を試す
- `OBSMODE='ON'` だけ抜いて確認する
- 位置中心部だけを抽出して viewer へ渡す
- FDNUM / IFNUM / PLNUM ごとに分けて解析する

---

### `select_rows()`
#### 何をするか
query を使わず、単に位置インデックスで Scantable を切り出すショートカットです。

```python
sc_small = su.select_rows(sc, "0:100")
```

これは
```python
sc_small = su.filter_scantable(sc, rows="0:100")
```
とほぼ同義です。

#### 向いている用途
- まず小さいサブセットで pipeline を試したい
- 不具合調査の最小再現セットを作りたい
- viewer / montage / grid に軽い入力だけ渡したい

---

### `scantable_utils.py` をどう使うと効果的か
以下の順で使うと、単なる補助関数集以上の価値があります。

1. `describe_columns()` / `show_scantable()` で入力確認
2. `find_scans()` / `filter_scantable()` / `select_rows()` で試験サブセット作成
3. 必要なら `calc_mapping_offsets()` で位置情報を追加
4. 複数日ファイルなら `merge_scantables(shift_scan_id=True)`
5. `set_beameff()` でビーム別効率を付与
6. 限定的に `update_metadata()` でラベル是正

この流れにより、後段の `calibrate`, `baseline`, `coadd`, `profile_view` がかなり安全になります。


## 8.11 `axis.py`
低レベルの軸処理です。

### 公開 API
```python
def freq_axis_from_wcs(meta: dict, nchan: int) -> np.ndarray
def radio_velocity_kms(freq_hz: np.ndarray, rest_hz: float) -> np.ndarray
def wcs_slice_channels(meta: dict, ch_start: int, ch_stop: int) -> dict
def channel_slice_from_vrange_union(meta, v_corr_kms, vmin_kms, vmax_kms, rest_hz=None) -> Tuple[int, int]
def slice_channels(meta, data, ch_start, ch_stop)
```

### 役割
- FITS WCS から周波数軸を作る
- radio definition の速度へ変換
- channel slice 後の `CRVAL1` 更新
- dump ごとの `v_corr` を考慮して、指定速度範囲を全 dump で覆う channel range を決める

`vlsrk_range_kms` を使った切り出しの根幹です。

---

## 8.12 `ranges.py`
窓文字列処理です。

### 公開 API
```python
def parse_windows(specs: List[str]) -> List[Tuple[float, float]]
def window_to_mask(x_axis, windows) -> np.ndarray
def windows_to_mask(x, windows) -> np.ndarray
```

### 使い方
- `["-35:-5", "30:55"]` のような指定を数値窓へ変換
- baseline 窓や RMS 窓の指定に使う

---

## 8.13 `doppler.py`
観測者速度補正と Doppler factor です。

### 公開 API
```python
def earth_location_from_meta(meta: dict)
def calc_vlsrk_correction_kms(...)
def compute_vcorr_series(...)
def get_doppler_factor(v_kms)
def scale_frequency_wcs_by_velocity(meta, v_corr_kms) -> dict
```

### 役割
- `meta` から観測地点を復元
- 行ごとの VLSRK 補正を計算
- 相対論的 Doppler factor を返す
- 周波数 WCS を Doppler 補正したコピーを作る

### `earth_location_from_meta()` の探索順
1. `OBSGEO-X/Y/Z`
2. `SITELAT`, `SITELONG`, `SITEELEV`
3. `site_name`

### `compute_vcorr_series()`
- `DatetimeIndex`
- RA/DEC
- meta のサイト情報

から各行の `v_corr` を計算します。

---

## 8.14 `profile_view`（旧 plotting）
viewer 本体は別マニュアル化する予定とのことなので、ここでは概要のみ記します。

### モジュール
- `profile_view/viewer.py`
- `profile_view/montage.py`
- `profile_view/grid.py`
- `profile_view/windowmask.py`
- `profile_view/utils.py`

### 役割分担
- `viewer.py`  
  1 本ずつスペクトルを見る
- `montage.py`  
  複数スペクトルを格子状に並べる
- `grid.py`  
  座標配置に基づいてスペクトルマップ表示を作る
- `windowmask.py`  
  fit windows や RMS windows の表示変換

### 共通思想
- profile_view 系は **個々の行の軸をその行の WCS から再計算** して描く方向を重視している
- 旧来の「全部を先に Standardizer へ通して 1 本化」ではなく、表示の正しさを優先する設計が見える
- `TEMPSCAL` と `BEAMEFF` を参照し、表示スケール切替ができる

### 代表コンストラクタ

#### `SpectralViewer`
```python
def __init__(
    self,
    scantable,
    rows=None,
    exclude_rows=None,
    xaxis="vel",
    axis_type=None,
    rest_freq=None,
    xrange=None,
    yrange=None,
    autoscale_y=True,
    show_fitwin_rms=True,
    show_fit_windows=True,
    smooth_mode="none",
    smooth_width=1,
    box_downsample=False,
    box_policy="trim",
    rms_on="raw",
    show_top_axis=True,
    save_dir=".",
    save_prefix="spectrum",
    save_pdf=None,
    max_pdf_pages=100,
    show=True,
    figsize=None,
    content_aspect=None,
    paper_margins=None,
)
```

#### `ProfileMapMontageViewer`
- `nrows`, `ncols`, `imin`, `imax` で表示ページを制御
- `annotate_rms`, `show_fit_windows`, `rest_freq` などを指定可能

#### `ProfileMapGridViewer`
- `coord`, `projection`, `ref_point`
- `x_grid`, `y_grid`
- `corner_offsets`, `grid_bounds_offsets`, `grid_anchor_offsets`
- `grid_tol`
- `mode`, `combine`
- `show_baseline`, `show_grid`, `offset_unit`

など、空間配置の指定が多いです。

---

# 9. obs_method.txt を一般化した推奨解析フロー

`obs_method.txt` では複数ファイルを loop で回していますが、説明書では **1 ステップずつ意味が分かる形** に直すと次のようになります。

## 9.1 12CO の典型フロー

1. raw FITS を読む  
2. `run_tastar_calibration()` で Ta* を作る  
3. `group_mode="scan"` で scan ごとにまとめる  
4. `run_baseline_fit()`  
5. 同一位置をまとめる  
6. 必要なら再度 baseline  
7. viewer で確認

この順番にする理由は、

- まず短時間の dump を scan 単位で安定化
- その段階で baseline をかける
- その後、同じ位置の複数 repeat を最終的に統合

という二段階 coadd が自然だからです。

## 9.2 13CO / C18O のように同じ raw から複数線を取りたい場合

同じ raw に対し、`rest_freq` を変えて複数回 `run_tastar_calibration()` すればよいです。  
このとき FREQ 軸のまま処理している限り、**元の観測周波数 WCS は不必要に壊さず** に済みます。

---

# 10. よくある設計判断と推奨

## 10.1 baseline は scan coadd の前か後か
推奨は次です。

- dump 生データに直接高次 baseline を乱用しない
- まず scan 単位で coadd
- その後 baseline
- 最終統合後に必要ならもう一度 baseline

## 10.2 `baseline_vwin` と `rms_vwin` の使い分け
- 「引いた後の残差 RMS を重みにしたい」なら `baseline_vwin`
- 「引かずにノイズだけ評価したい」なら `rms_vwin`

## 10.3 TOPOCENT のまま baseline をしてよいか
可能ですが、速度窓 `vwin` を使うなら、その行の `VELOSYS` が必要です。  
このパッケージはその整合を厳しく見ています。

## 10.4 viewer 用に TR* へ変換したい
解析本体は Ta* のまま進め、viewer または write 時に変換する方が安全です。  
特に `BEAMEFF` が行ごとに違う場合、早い段階で一括変換すると追跡が難しくなります。

---

# 11. エラーや停止条件の読み方

## 11.1 `TOPOCENT input requires VELOSYS/VFRAME or valid timestamps`
意味:
- TOPOCENT なのに速度補正列が無く
- しかも時刻が壊れていて再計算もできない

対策:
- `TIMESTAMP` / `MJD` / `DATE-OBS` を整える
- `RA`, `DEC` と site 情報も整える

## 11.2 `Input has SPECSYS=LSRK but contains non-zero unapplied velocity column`
意味:
- すでに LSRK と名乗っているのに、未適用速度補正が残っている

対策:
- `VELOSYS`, `VFRAME` を 0 または削除
- あるいは本当に TOPOCENT なら `SPECSYS` を直す

## 11.3 `Cannot specify both baseline_vwin and rms_vwin`
意味:
- coadd の重み付け方針が曖昧なので禁止

## 11.4 `RESTFRQ missing`
意味:
- 速度解釈や rest override に必要な静止周波数がない

## 11.5 TR* 変換時の `BEAMEFF` エラー
意味:
- `BEAMEFF <= 0` や NaN がある

対策:
- `set_beameff()` で設定
- あるいは `out_scale="TA*"` のまま扱う

---

# 12. 実務上の推奨チェックリスト

- 入力 table に `TIMESTAMP` が正しくあるか
- `SPECSYS` と `SSYSOBS` が妥当か
- `RESTFRQ` / `RESTFREQ` が入っているか
- TOPOCENT なら `RA/DEC` とサイト情報があるか
- `OBSMODE=HOT/OFF/ON` が対応しているか
- baseline 窓と line 窓が衝突していないか
- 13CO / C18O 再解釈時は `rest_freq` を明示したか
- coadd 前に `SCAN`, `INTGRP`, `pos_id` などの grouping 列があるか
- viewer で TR* を見たいなら `BEAMEFF` を設定したか

---

# 13. 最低限の API 早見表

## 13.1 解析で最もよく使う関数
- `fitsio.read_scantable`
- `calibrate.run_tastar_calibration`
- `baseline.run_baseline_fit`
- `coadd.run_velocity_coadd`
- `scantable_utils.show_scantable`
- `profile_view.viewer.view_spectra`

## 13.2 補助的だが重要
- `restfreq.apply_restfreq_override`
- `tempscale.set_beameff`
- `scantable_utils.merge_scantables`
- `scantable_utils.calc_mapping_offsets`
- `regrid_vlsrk.Standardizer`

---

# 14. 付録 A: 公開 API シグネチャ一覧

## 14.1 I/O
```python
def read_scantable(path: str, *, tr_input_policy: str = "preserve") -> Scantable:
def write_scantable(
    path: str,
    scantable: Scantable,
    spectrum_column: str = "DATA",
    overwrite: bool = True,
    **kwargs
) -> None:
def write_sdfits(
    out_path: str,
    meta: dict,
    data: Union[np.ndarray, List[np.ndarray]],
    table: pd.DataFrame,
    history: dict | None = None,
    *,
    spectrum_column: str = "DATA",
    include_flag: bool = True,
    overwrite: bool = True,
    string_widths: dict[str, int] | None = None,
    **kwargs,
) -> None:
```

## 14.2 Raw
```python
def build_rawspec(
    *,
    hot: pd.DataFrame,
    on: pd.DataFrame,
    off: pd.DataFrame,
    meta: dict[str, Any] | None,
    mapping: pd.DataFrame | None,
) -> RawSpec:
def save_rawspec(raw: RawSpec, path: str) -> None:
def load_rawspec(path: str) -> RawSpec:
def load_rawspec_auto(path: str, prefer: tuple[str, ...] | None = None) -> RawSpec:
def load_rawspec_fits(path: str) -> RawSpec:
```

## 14.3 Calibration
```python
def make_tastar_dumps(
    raw,
    *,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    t_hot_k: Optional[float] = None,
    tau_zenith: Union[float, str, None] = None,
    t_surface_k: Optional[float] = None,
    t_atm_model: str = "offset",
    t_atm_delta_k: float = 15.0,
    t_atm_eta: float = 0.95,
    gain_mode: str = "hybrid",
    verbose: bool = True,
    ch_range: Optional[Tuple[int, int]] = None,
    vlsrk_range_kms: Optional[Tuple[float, float]] = None,
    v_corr_col: str = "VELOSYS",
    coord_frame: str | None = None,
    vcorr_chunk_sec: Optional[float] = None,
    dtype: Optional[Union[type, str]] = None,
    rest_freq: Optional[float] = None,
) -> Scantable:

def run_tastar_calibration(
    input_data: Union[str, Scantable, dict],
    output_path: Optional[str] = None,
    t_hot_k: Optional[float] = None,
    ch_range: Optional[Tuple[int, int]] = None,
    vlsrk_range_kms: Optional[Tuple[float, float]] = None,
    coord_frame: str | None = None,
    spectrum_column: str = "DATA",
    overwrite: bool = False,
    store_freq_column: bool = False,
    vcorr_chunk_sec: Optional[float] = None,
    dtype: Optional[Union[type, str]] = None,
    rest_freq: Optional[float] = None,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    tau_zenith: Union[float, str, None] = None,
    t_surface_k: Optional[float] = None,
    t_atm_model: str = "offset",
    t_atm_delta_k: float = 15.0,
    t_atm_eta: float = 0.95,
    gain_mode: str = "hybrid",
    verbose: bool = True,
) -> Scantable:

def recalibrate_tastar(
    scantable: Scantable,
    new_tau: Union[float, None] = None,
    new_t_surface_k: Optional[float] = None,
    new_t_atm_model: str = "offset",
    new_t_atm_delta_k: float = 15.0,
    new_t_atm_eta: float = 0.95,
    verbose: bool = True
) -> Scantable:
```

## 14.4 Baseline
```python
def fit_polynomial_baseline(
    v_kms: np.ndarray,
    y: np.ndarray,
    base_windows: Optional[List[Tuple[float, float]]] = None,
    *,
    windows: Optional[List[Tuple[float, float]]] = None,
    line_windows: Optional[List[Tuple[float, float]]] = None,
    poly_order: int = 1,
    order: Optional[int] = None,
    iter_max: int = 0,
    iter_sigma: float = 3.0
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:

def run_baseline_fit(
    input_data: Union[str, Scantable],
    output_path: Optional[str] = None,
    *,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    vwin: List[str],
    poly_order: int = 1,
    line_vwin: Optional[List[str]] = None,
    iter_max: int = 0,
    iter_sigma: float = 3.0,
    max_dumps: int = 0,
    ch_start: Optional[int] = None,
    ch_stop: Optional[int] = None,
    v_corr_col: str = "VELOSYS",
    rest_freq: Optional[float] = None,
    on_fail: str = "exit",
    overwrite: bool = True
) -> Scantable:
```

## 14.5 Regrid / Standardizer
```python
def make_vgrid(vmin_kms: float, vmax_kms: float, dv_kms: float) -> VGrid:
def vlsrk_axis_for_spectrum(meta: dict, *, v_corr_kms: float, nchan: Optional[int] = None) -> np.ndarray:
def interp_to_vgrid(v_src: np.ndarray, y_src: np.ndarray, v_tgt: np.ndarray) -> np.ndarray:
def vrange_from_meta_and_vcorr(meta: dict, v_corr_kms: np.ndarray, *, nchan: int | None = None) -> tuple[float, float]:
```

## 14.6 Coadd
```python
def run_velocity_coadd(
    inputs: Sequence[Union[str, Scantable]],
    output_path: Optional[str] = None,
    *,
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    mode: str = "uniform",
    group_mode: str = "position",
    pos_col: str = "pos_id",
    pos_tol_arcsec: Optional[float] = None,
    v_corr_col: str = "VELOSYS",
    coord_frame: Optional[str] = None,
    vcorr_chunk_sec: Optional[float] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    dv: Optional[float] = None,
    allow_outside_overlap: bool = False,
    axis_type: str = "freq",
    rest_freq: Optional[float] = None,
    baseline_vwin: Union[str, List[str], None] = None,
    baseline_poly: int = 0,
    baseline_iter_max: int = 0,
    baseline_iter_sigma: float = 3.0,
    rms_vwin: Union[str, List[str], None] = None,
    rms_poly: int = 1,
    rms_bin: int = 1,
    coadd_qc: Optional[str] = None,
    line_vwin: Union[str, List[str], None] = None,
    block_size: int = 0,
    max_dumps: int = 0,
    ch_start: Optional[int] = None,
    ch_stop: Optional[int] = None,
    on_fail: str = "exit",
    overwrite: bool = True,
    out_scale: str = "TA*",
    normalize_if_mixed: str = "auto",
    beameff_tol: float = 1e-6,
    sigma_scale: str = "TA*",
    verbose: bool = True,
) -> Scantable:
```

## 14.7 Utilities
```python
def apply_restfreq_override(
    meta: dict,
    table: Optional[pd.DataFrame],
    rest_freq_hz: float,
    *,
    require_wcs_for_vrad: bool = True,
) -> Dict[str, Any]:

def show_scantable(...)
def calc_mapping_offsets(...)
def merge_scantables(...)
def update_metadata(...)
def find_scans(...)
def filter_scantable(...)
def select_rows(...)
```

---

# 15. 最後に

このパッケージの本質は、単に「スペクトルを足す」ことではなく、

- 基準系の意味を崩さない
- rest frequency の変更を安全に行う
- 温度スケールを破壊的に扱わない
- 行ごとの WCS と metadata を尊重する

という点にあります。

特に **Standardizer**, **coadd**, **restfreq**, **TEMPSCAL/BEAMEFF**, **TOPOCENT/LSRK の意味論** は、解析結果の信頼性に直結します。  
解析例だけを見るより、この説明書の「なぜそう実装されているか」を先に理解しておくと、今後の拡張やデバッグが格段に楽になります。
