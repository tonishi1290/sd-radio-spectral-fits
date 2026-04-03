# pointing-fit: OMU 1.85m 実践例

対象: `pointing_fit`

この文書は、OMU 1.85m で pointing fit 解析を行うときの、
manifest の書き方、normalize の解釈、fit の進め方を実務向けにまとめたものです。

---

## 1. OMU 1.85m の基本設定

通常の OMU 1.85m optical では、次を基本にします。

- raw `dx`, `dy` は **星位置 - CCD中心**
- `measurement_space = "star_offset"`
- 画像向きは装置に応じて `positive_*` を決める
- 観測プログラムの器差適用式は
  - `Az_cmd = Az_true + Err_Az`
  - `El_cmd = El_true + Err_El`

したがって、典型 manifest は次です。

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

`left/down` は例です。実画像の向きに合わせて `right` や `up` に変えてください。

---

## 2. normalize

```bash
pointing-fit normalize \
  --manifest datasets.toml \
  --output normalized.csv
```

単発指定でも実行できます。

```bash
pointing-fit normalize \
  --input d_param.csv \
  --dataset-id run001 \
  --used-param pointing_param.toml \
  --telescope omu1p85m \
  --model omu1p85m_actual_v1 \
  --measurement-space star_offset \
  --positive-az-moves-star left \
  --positive-el-moves-star down \
  --command-err-mode add \
  --angle-unit deg \
  --output normalized.csv
```

---

## 3. OMU 1.85m での `command_err_mode = add` の意味

OMU 1.85m では通常

```text
Az_cmd = Az_true + Err_Az
El_cmd = El_true + Err_El
```

です。したがって、必要指令補正 `Δcmd` と器差差分 `ΔErr` の関係は

```text
ΔErr = +Δcmd
```

です。

つまり normalize 後の `delta_err_*` は、raw から見て符号をそのまま保ちやすい設定です。

---

## 4. fit

### 4.1 absolute fit

```bash
pointing-fit fit \
  --input normalized.csv \
  --outdir fit_abs \
  --model omu1p85m_actual_v1 \
  --fit-target absolute
```

### 4.2 delta fit

```bash
pointing-fit fit \
  --input normalized.csv \
  --outdir fit_delta \
  --model omu1p85m_actual_v1 \
  --fit-target delta
```

absolute fit はこの版では常に

```text
new_model = used_model + delta_err
```

です。

---

## 5. 旧 TOML をどう扱うか

OMU 1.85m でも、観測時に使っていた `pointing_param.toml` を保存しておくことを強く勧めます。

- absolute fit では、観測時 TOML と新 fit の差が分かる
- delta fit では、その TOML に対する増分として結果を適用できる

---

## 6. 何を見るべきか

normalize 後:

- `raw_dx_input_arcsec`, `raw_dy_input_arcsec`
- `delta_err_dx_arcsec`, `delta_err_dy_arcsec`
- `used_model_dx_arcsec`, `used_model_dy_arcsec`

fit 後:

- `parameter_updates.txt`
- `fit_result.json`

`apply` に渡すのは `fit_result.json` ではなく、delta fit で作られる `delta_params.toml` です。raw CSV の `dt_cross_id` / `cross_id` は省略可能で、無い場合は normalize が自動生成します。
- `residuals.csv`
- `fit_diagnostics.png`
- `fit_diagnostics_azel_map.png`

---

## 7. 注意点

1. `sign_convention` はこの版で削除しました。OMU 1.85m 運用でも考えなくて構いません。
2. `positive_*` は装置の実画像向きに合わせて必ず確認してください。
3. raw `dx`, `dy` が本当に「星位置 - CCD中心」なら、`star_offset` を使うのが自然です。
