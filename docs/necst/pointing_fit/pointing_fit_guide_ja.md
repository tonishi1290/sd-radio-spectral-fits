# pointing-fit 説明書

対象: `pointing_fit`

この文書は、日常的に `pointing-fit` を使うための実務向けガイドです。

---

## 1. このツールでやること

`pointing-fit` は、観測から得た `d_param.csv` を正規化し、
pointing model を fit し、TOML と診断図を書き出します。

基本の流れは次の 4 段階です。

1. `normalize`
   - raw CSV を canonical CSV に変換する
2. `fit`
   - model parameter を推定する
3. `apply`
   - `delta` fit の結果を既存 TOML に加える
   - `fit_result.json` や `absolute` fit の出力は直接は受けない
4. `report`
   - residual 図と summary markdown を作る

---

## 2. 今回の仕様整理で重要な点

### 2.1 raw `dx`, `dy` の意味

通常の pointing fit では、raw `dx`, `dy` は

- **星位置 - CCD中心**

です。

したがって、ユーザー向けの通常運用では

- `measurement_space = "star_offset"`

を使います。

### 2.2 画像の向きと観測プログラムの向きは分けて扱う

画像の向きは

- `positive_az_moves_star`
- `positive_el_moves_star`

で与えます。

器差を観測プログラムへどう入れているかは

- `command_err_mode`

で与えます。

### 2.3 raw 入力列の単位

Az/El は `az_deg`, `el_deg`, `az_rad`, `el_rad` のように列名で単位を明示できます。
raw `dx`, `dy` も `dx_arcsec`, `dy_arcsec`, `dx_deg`, `dy_deg`, `dx_rad`, `dy_rad` のように列名で単位を明示できます。
`dt_cross_id` と `cross_id` は省略可能です。無い場合は normalize が行番号から自動生成します。
Az/El 列名に単位が無い場合だけ `--angle-unit` / manifest の `angle_unit` を使います。raw `dx`, `dy` 列名に単位が無い場合は `--offset-unit` / manifest の `offset_unit` を使います。既定値は legacy `d_param.csv` 互換のため `arcsec` です。

### 2.4 normalize 後の量の意味

normalize 後の

- `delta_err_dx_arcsec`
- `delta_err_dy_arcsec`

は、常に

- **使用中モデルに加えるべき差分**

です。

### 2.5 `sign_convention` は廃止されました

以前は `absolute` fit の解釈に `sign_convention` を使っていましたが、
この版では **削除** です。manifest に残っている場合は normalize が明示的にエラーで停止します。legacy な CSV 列として残っていても fit では無視されます。

---

## 3. 最小実行例

### 3.1 NANTEN2

```bash
pointing-fit normalize \
  --input ./d_param.csv \
  --dataset-id nanten2_run001 \
  --used-param ./pointing_param.toml \
  --telescope nanten2 \
  --model nanten2_actual_v1 \
  --measurement-space star_offset \
  --positive-az-moves-star left \
  --positive-el-moves-star down \
  --command-err-mode subtract \
  --angle-unit deg \
  --output normalized_nanten2.csv

pointing-fit fit \
  --input normalized_nanten2.csv \
  --outdir fit_nanten2_abs \
  --model nanten2_actual_v1 \
  --fit-target absolute
```

### 3.2 OMU 1.85m

```bash
pointing-fit normalize \
  --input ./d_param.csv \
  --dataset-id omu_run001 \
  --used-param ./pointing_param.toml \
  --telescope omu1p85m \
  --model omu1p85m_actual_v1 \
  --measurement-space star_offset \
  --positive-az-moves-star left \
  --positive-el-moves-star down \
  --command-err-mode add \
  --angle-unit deg \
  --output normalized_omu.csv

pointing-fit fit \
  --input normalized_omu.csv \
  --outdir fit_omu_abs \
  --model omu1p85m_actual_v1 \
  --fit-target absolute
```

---

## 4. `command_err_mode` の意味

### 4.1 `subtract`

```text
Az_cmd = Az_true - Err_Az
El_cmd = El_true - Err_El
```

したがって、指令補正 `Δcmd` と器差差分 `ΔErr` の関係は

```text
ΔErr = -Δcmd
```

です。

NANTEN2 は通常こちらです。

### 4.2 `add`

```text
Az_cmd = Az_true + Err_Az
El_cmd = El_true + Err_El
```

したがって

```text
ΔErr = +Δcmd
```

です。

OMU 1.85m は通常こちらです。

---

## 5. normalized CSV の読み方

主に見るべき列は次です。

- `raw_dx_input_arcsec`, `raw_dy_input_arcsec`
  - 入力値を arcsec にそろえたもの
- `delta_err_dx_arcsec`, `delta_err_dy_arcsec`
  - 使用中モデルに加える差分
- `used_model_dx_arcsec`, `used_model_dy_arcsec`
  - 観測時に使っていた TOML に基づく model 値
- `az_deg`, `el_deg`
  - fit に使う Az, El

旧名互換列として

- `measured_dx_arcsec`, `measured_dy_arcsec`
- `applied_model_dx_arcsec`, `applied_model_dy_arcsec`

も残りますが、今後は新名を見る方が安全です。

---

## 6. `absolute` と `delta`

### 6.1 `delta`

`delta` fit は、normalize 後の `delta_err_*` 自体を fit します。

- 何を追加すべきかを知りたいとき
- 既存 TOML を基準に少し更新したいとき

に向いています。

### 6.2 `absolute`

`absolute` fit は、常に

```text
new_model = used_model + delta_err
```

を fit target にします。

したがって、この版では `sign_convention` を考える必要はありません。

---

## 7. manifest の使い方

共通規約を `[defaults]` に置き、各 dataset のファイル名だけを `[[datasets]]` に並べるのが実務上分かりやすいです。

### 7.1 NANTEN2 の例

```toml
[defaults]
angle_unit = "deg"
measurement_space = "star_offset"
positive_az_moves_star = "left"
positive_el_moves_star = "down"
command_err_mode = "subtract"
model = "nanten2_actual_v1"

[[datasets]]
dataset_id = "run001"
input = "d_param.csv"
used_param = "pointing_param.toml"
```

### 7.2 OMU 1.85m の例

```toml
[defaults]
angle_unit = "deg"
measurement_space = "star_offset"
positive_az_moves_star = "left"
positive_el_moves_star = "down"
command_err_mode = "add"
model = "omu1p85m_actual_v1"

[[datasets]]
dataset_id = "run001"
input = "d_param.csv"
used_param = "pointing_param.toml"
```

---

## 8. いつ `command_offset` / `model_delta` を使うか

通常運用では不要です。

- `star_offset`
  - ユーザー向け標準
- `command_offset`
  - upstream がすでに Az/El 正方向の command 補正量を出している特殊ケース
- `model_delta`
  - synthetic 検証や内部確認用

これらは便利ですが、日常運用で主役にしない方が混乱が少なくなります。

---

## 9. 推奨運用

1. raw `dx`, `dy` は常に「星位置 - CCD中心」で保存する
2. manifest では `measurement_space = "star_offset"` を使う
3. 画像の向きは `positive_*` に書く
4. 観測プログラムの器差適用式は `command_err_mode` に書く
5. `sign_convention` は使わない
6. normalize 後は `delta_err_*` を見る

---

## 10. よくある誤解

### 誤解 1
`positive_az_moves_star` は raw `dx` の符号規約そのものを表す。

**違います。**
これは、+Az 指令で星が画像のどちらへ動くか、という装置の向きを表します。

### 誤解 2
`command_err_mode` と旧 `sign_convention` は同じ意味ではありません。

**違います。**
現在の版では、`sign_convention` はこの版で削除しました。実務的に効くのは `command_err_mode` です。

### 誤解 3
`command_offset` の方が常に正しい。

**違います。**
ユーザーが自然に扱うのは `star_offset` です。`command_offset` は内部向け・特殊ケース向けです。


## 10. fixed parameter の実務方針

NANTEN2 では、`g, gg, ggg, gggg` と対応する radio 項を、最初から高次まで一度に free にしないことを推奨します。

- まず低次項だけを free にする
- 高次項は固定したまま残差を見る
- 必要な場合だけ `ggg`, `gggg` と対応する radio 項を段階的に free にする

同梱例:

- `examples/configs/fixed_nanten2_optical.toml`
- `examples/configs/fixed_nanten2_radio.toml`

optical / radio どちらでも、最初から高次まで free にすると他の項を吸収しやすくなるため、この段階的戦略を基本にしてください。
