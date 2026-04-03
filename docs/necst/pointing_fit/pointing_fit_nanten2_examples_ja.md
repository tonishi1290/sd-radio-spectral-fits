# pointing-fit: NANTEN2 実践例

対象: `pointing_fit`

この文書は、NANTEN2 の optical / radio pointing fit を実務で回すときの、
manifest の書き方、normalize の考え方、fit の使い分けを簡潔にまとめたものです。

---

## 1. NANTEN2 で最初に固定するもの

通常の NANTEN2 optical では、次を前提にします。

- raw `dx`, `dy` は **星位置 - CCD中心**
- `measurement_space = "star_offset"`
- `+Az` で星が画像の左へ動く
- `+El` で星が画像の下へ動く
- 観測プログラムの器差適用式は
  - `Az_cmd = Az_true - Err_Az`
  - `El_cmd = El_true - Err_El`

したがって manifest は次が基本形です。

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

---

## 2. normalize

```bash
pointing-fit normalize \
  --manifest datasets.toml \
  --output normalized.csv
```

あるいは単発指定でも構いません。

```bash
pointing-fit normalize \
  --input d_param.csv \
  --dataset-id run001 \
  --used-param pointing_param.toml \
  --telescope nanten2 \
  --model nanten2_actual_v1 \
  --measurement-space star_offset \
  --positive-az-moves-star left \
  --positive-el-moves-star down \
  --command-err-mode subtract \
  --angle-unit deg \
  --output normalized.csv
```

---

## 3. NANTEN2 での `command_err_mode = subtract` の意味

NANTEN2 では通常

```text
Az_cmd = Az_true - Err_Az
El_cmd = El_true - Err_El
```

です。したがって、画像から得た必要指令補正 `Δcmd` を器差差分に直すと

```text
ΔErr = -Δcmd
```

になります。

normalize 後の `delta_err_*` は、この `ΔErr` です。

---

## 4. absolute fit と delta fit

### 4.1 absolute fit

```bash
pointing-fit fit \
  --input normalized.csv \
  --outdir fit_abs \
  --model nanten2_actual_v1 \
  --fit-target absolute
```

この版では absolute fit は常に

```text
new_model = used_model + delta_err
```

です。

### 4.2 delta fit

```bash
pointing-fit fit \
  --input normalized.csv \
  --outdir fit_delta \
  --model nanten2_actual_v1 \
  --fit-target delta
```

`delta` fit の結果は、既存 TOML に対する差分として読むのが基本です。

---

## 5. fixed parameter 運用

NANTEN2 では縮退しやすい組み合わせがあるため、実務では fixed parameter を使うことが多いです。

```bash
pointing-fit fit \
  --input normalized.csv \
  --outdir fit_abs \
  --model nanten2_actual_v1 \
  --fit-target absolute \
  --fix-file examples/configs/fixed_nanten2_optical.toml
```

添付した fixed 例として、次を同梱しています。

- `examples/configs/fixed_nanten2_optical.toml`
- `examples/configs/fixed_nanten2_radio.toml`

推奨戦略は、`g, gg, ggg, gggg` と対応する radio 項を**できるだけ低次から順に free にする**ことです。
最初から高次項を free にすると、他の項まで吸収してしまう可能性があります。

したがって実務では、まず

1. 高次の重力項を固定した状態で fit する
2. 残差や相関を見る
3. 必要な場合だけ `ggg`, `gggg` や対応する radio 項を順に解放する

という手順を推奨します。

radio pointing では、まず `examples/configs/fixed_nanten2_radio.toml` を出発点にして optical 既知項を固定し、radio 項を段階的に解放してください。

---

## 6. 何を見るべきか

normalize 後は、まず次を確認します。

- `raw_dx_input_arcsec`, `raw_dy_input_arcsec`
- `delta_err_dx_arcsec`, `delta_err_dy_arcsec`
- `used_model_dx_arcsec`, `used_model_dy_arcsec`

fit 後は

- `fit_result.json`

`apply` に渡すのは `fit_result.json` ではなく、delta fit で作られる `delta_params.toml` です。raw CSV の `dt_cross_id` / `cross_id` は省略可能で、無い場合は normalize が自動生成します。
- `parameter_updates.txt`
- `residuals.csv`
- `fit_diagnostics.png`
- `fit_diagnostics_azel_map.png`

を見ます。

---

## 7. 注意点

1. `sign_convention` はこの版で削除しました。NANTEN2 運用では考えなくて構いません。
2. 画像向きが本当に `left/down` かは、装置セットアップごとに必ず一度確認してください。
3. 観測時の `pointing_param.toml` を保存しておくと、absolute fit の比較が明快になります。
