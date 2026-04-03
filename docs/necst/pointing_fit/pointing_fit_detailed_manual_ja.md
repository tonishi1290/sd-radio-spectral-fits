# pointing-fit 詳細説明書

対象: `pointing_fit`

この文書は、`pointing-fit` CLI の仕様、入力データの意味、各サブコマンドの動作、
出力ファイルの読み方、NANTEN2 / OMU 1.85m での実運用上の注意を詳しくまとめたものです。

---

## 1. 全体像

`pointing-fit` の主な役割は次の通りです。

1. 観測から得られた raw CSV を、fit 用の canonical CSV へ正規化する
2. pointing model を `absolute` または `delta` として fit する
3. 必要なら `delta` を既存 TOML に加える
4. residual 図と summary markdown を出す

本ツールでは、入力 `dx`, `dy` の意味が曖昧だと、fit の符号議論が破綻しやすくなります。
そのため、この版では raw `dx`, `dy` の意味を明示的に固定して扱います。

---

## 2. 用語の定義

### 2.1 raw 入力

- `dx_raw`, `dy_raw`
  - raw CSV に入っている元の `dx`, `dy`
  - 通常運用では **星位置 - CCD中心**

### 2.2 指令補正

- `Δcmd_x`, `Δcmd_y`
  - 星を中心へ戻すために追加で与える command 補正量
  - 正方向は +Az, +El

### 2.3 器差差分

- `ΔErr_x`, `ΔErr_y`
  - 観測時に使っていた器差モデルへ加えるべき差分

### 2.4 normalize 後の列

この版では

- `delta_err_dx_arcsec`
- `delta_err_dy_arcsec`

を、必ず

- **使用中モデルに加えるべき差分**

として解釈します。

### 2.5 使用中モデル

- `used_model_dx_arcsec`
- `used_model_dy_arcsec`

は、観測時に使っていた TOML と Az/El から計算した model 値です。

---

## 3. `measurement_space` の意味

### 3.1 `star_offset`

ユーザー向け標準です。

raw `dx`, `dy` は

- 星位置 - CCD中心

として解釈します。

画像の向きは

- `positive_az_moves_star`
- `positive_el_moves_star`

で指定します。

### 3.2 `command_offset`

upstream がすでに

- +Az に何秒角動かすべきか
- +El に何秒角動かすべきか

という形で補正量を出している場合に使います。

通常の pointing fit 実務では、主役にしない方が混乱が少なくなります。

### 3.3 `model_delta`

raw `dx`, `dy` が最初から model difference として与えられている場合です。
synthetic 検証や内部確認向けです。

---

## 4. 画像の向き: `positive_az_moves_star`, `positive_el_moves_star`

### 4.1 `positive_az_moves_star`

`+Az` 指令を出したとき、星が画像上のどちらへ動くかを表します。

- `left`
- `right`

### 4.2 `positive_el_moves_star`

`+El` 指令を出したとき、星が画像上のどちらへ動くかを表します。

- `up`
- `down`

これは raw `dx`, `dy` の数値そのものではなく、**装置の見え方**を記述するための設定です。

---

## 5. 観測プログラムの器差適用式: `command_err_mode`

### 5.1 `subtract`

```text
Az_cmd = Az_true - Err_Az
El_cmd = El_true - Err_El
```

このとき

```text
ΔErr = -Δcmd
```

です。

NANTEN2 は通常こちらです。

### 5.2 `add`

```text
Az_cmd = Az_true + Err_Az
El_cmd = El_true + Err_El
```

このとき

```text
ΔErr = +Δcmd
```

です。

OMU 1.85m は通常こちらです。

---

## 6. normalize の仕様

### 6.1 何をしているか

`normalize` は、raw CSV を読み、

- Az, El を deg にそろえる
- raw `dx`, `dy` を arcsec にそろえる
- Az/El は列名が `*_rad` なら rad, `*_deg` なら deg と解釈し, それ以外は `angle_unit` を使う
- raw `dx`, `dy` も列名で単位を明示できる。代表例は `dx_arcsec`, `dy_arcsec`, `dx_deg`, `dy_deg`, `dx_rad`, `dy_rad`
- Az/El 列名に単位が無い場合は `angle_unit`、raw `dx`, `dy` 列名に単位が無い場合は `offset_unit` を使う
- `offset_unit` の既定値は legacy `d_param.csv` 互換のため `arcsec`
- `dt_cross_id`, `cross_id` は省略可能。省略時は normalize が行番号から自動生成する
- `measurement_space`, `positive_*`, `command_err_mode` に基づいて内部差分へ変換する
- 使用中 TOML から `used_model_*` を計算する

という処理をします。

### 6.2 出力列

主要列:

- `dataset_id`
- `az_deg`, `el_deg`
- `raw_dx_input_arcsec`, `raw_dy_input_arcsec`
- `delta_err_dx_arcsec`, `delta_err_dy_arcsec`
- `used_model_dx_arcsec`, `used_model_dy_arcsec`

旧互換列:

- `measured_dx_arcsec`, `measured_dy_arcsec`
- `applied_model_dx_arcsec`, `applied_model_dy_arcsec`

旧互換列は残りますが、今後は新名を基準に読む方が安全です。

### 6.3 `sign_convention` は廃止されました

この版では 削除 です。入力や既存 CSV に列が残っていても、
fit の実際の解釈は

```text
new_model = used_model + delta_err
```

に固定されています。

---

## 7. `fit` の仕様

### 7.1 `fit-target = delta`

fit target はそのまま

- `delta_err_dx_arcsec`
- `delta_err_dy_arcsec`

です。

### 7.2 `fit-target = absolute`

fit target は常に

```text
new_model = used_model + delta_err
```

です。

したがって、この版では `absolute` fit に `sign_convention` の分岐はありません。

### 7.3 `solve-mode`

- `joint`
  - dx と dy をまとめて最適化する
- `separate`
  - dx と dy を別々に扱う

実務では、まず `joint` を試し、縮退や不安定さが強い場合だけ `separate` を検討するのがよいです。

### 7.4 robust loss

利用可能な値:

- `linear`
- `soft_l1`
- `huber`
- `cauchy`

通常は `soft_l1` を出発点にするのが無難です。

---

## 8. `apply` の仕様

`apply` は、`delta` fit で得た TOML を既存 TOML に加えるためのコマンドです。`fit_result.json` や `absolute fit` の出力ディレクトリは update source として受けません。

```bash
pointing-fit apply \
  --base-param pointing_param.toml \
  --delta-param fit_delta/delta_params.toml \
  --output pointing_param_updated.toml
```

`absolute` fit の結果は、通常はそのまま新しい候補 TOML として扱います。
`apply` に渡すのは通常は `delta_params.toml`（または separate の `delta_params_dx.toml` / `delta_params_dy.toml` を含む delta fit 出力ディレクトリ）です。`delta_params_arcsec.toml` は旧版互換として読むことはできますが、現行版は新規には出力しません。

---

## 9. `inspect` と `report`

### 9.1 `inspect`

入力 CSV に必要な列がそろっているか、基本的な状態を確認する簡易ツールです。

### 9.2 `report`

`residuals.csv` から summary 図と markdown を作ります。

主な出力:

- `fit_diagnostics.png`
- `fit_diagnostics_azel_map.png`
- `report_summary.md`

---

## 10. manifest の書き方

### 10.1 基本形

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

### 10.2 NANTEN2

```toml
[defaults]
angle_unit = "deg"
measurement_space = "star_offset"
positive_az_moves_star = "left"
positive_el_moves_star = "down"
command_err_mode = "subtract"
model = "nanten2_actual_v1"
```

### 10.3 OMU 1.85m

```toml
[defaults]
angle_unit = "deg"
measurement_space = "star_offset"
positive_az_moves_star = "left"
positive_el_moves_star = "down"
command_err_mode = "add"
model = "omu1p85m_actual_v1"
```

---

## 11. NANTEN2 / 1.85m の違いの見方

### 11.1 NANTEN2

- 器差適用式は通常 `subtract`
- したがって、指令補正と器差差分の符号は反転する
- `delta_err = -Δcmd`

### 11.2 OMU 1.85m

- 器差適用式は通常 `add`
- したがって、指令補正と器差差分の符号は同じ
- `delta_err = +Δcmd`

この違いが `command_err_mode` に集約されます。

---

## 12. よくある詰まりどころ

### 12.1 raw `dx`, `dy` がすでに Az/El 正方向なのに `star_offset` を使ってしまう

この場合は、画像向きの指定と二重になる可能性があります。
upstream が command 補正量を直接出しているなら、`command_offset` を検討してください。

### 12.2 観測時の TOML が残っていない

`absolute` fit の解釈と比較を正しく行うために、`used_param` はできるだけ保存しておくべきです。

### 12.3 旧版 `sign_convention` について

この版では変える必要はありません。削除 です。

---

## 13. 推奨運用

1. raw `dx`, `dy` は「星位置 - CCD中心」で保存する
2. 通常は `measurement_space = "star_offset"`
3. `positive_*` で画像向きを書く
4. `command_err_mode` で観測プログラムの器差適用式を書く
5. normalize 後は `delta_err_*` を見る
6. `sign_convention` は使わない
7. NANTEN2 / 1.85m の違いは `command_err_mode` を主軸に整理する


## 14. fixed parameter の推奨戦略

NANTEN2 では、重力項 `g, gg, ggg, gggg` と、それに対応する radio 項のあいだで縮退が起きやすいため、
最初から高次まで free にすることは推奨しません。

基本方針は次の通りです。

1. まず低次項だけを free にする
2. `ggg`, `gggg` と対応する radio 高次項は固定したまま開始する
3. 残差分布、Az/El 依存、パラメータ相関を確認する
4. 必要な場合だけ次数を一段ずつ上げる

同梱の fixed 例:

- `examples/configs/fixed_nanten2_optical.toml`
- `examples/configs/fixed_nanten2_radio.toml`

optical 用 fixed では、radio 項を基本的に固定し、まず optical の主要項を安定に求めることを狙います。
radio 用 fixed では、既知の optical 項を固定したうえで `g, gg, ggg, gggg = 0` から開始し、必要なときだけ高次項を順に解放します。

この方針の理由は、高次項を最初から free にすると、ゼロ点項や位相項など他の物理的意味を持つ項まで吸収しやすく、
解の安定性と解釈性が落ちるためです。
