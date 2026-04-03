# 測定規約メモ（2026-03-24版）

この版では、raw `dx`, `dy` の意味を次のように整理します。

## 1. 基本方針

ユーザーが与える raw `dx`, `dy` は、原則として

- **星位置 - CCD中心**

です。したがって、通常の pointing fit では

- `measurement_space = "star_offset"`

を使います。

## 2. `positive_az_moves_star`, `positive_el_moves_star`

これらは、画像の向きを装置ごとに指定するための項目です。

- `positive_az_moves_star`
  - `+Az` 指令で、星が画像上のどちらへ動くか
  - `left` または `right`
- `positive_el_moves_star`
  - `+El` 指令で、星が画像上のどちらへ動くか
  - `up` または `down`

## 3. `command_err_mode`

これは、観測プログラムで器差 `Err` をどの向きで指令へ入れているかを表します。

- `subtract`
  - `Az_cmd = Az_true - Err_Az`
  - `El_cmd = El_true - Err_El`
- `add`
  - `Az_cmd = Az_true + Err_Az`
  - `El_cmd = El_true + Err_El`

この違いにより、normalize 後の `delta_err_*` の符号が決まります。

## 4. raw 入力列の単位指定

raw 入力の単位は列名でも明示できます。代表例は

- Az/El: `az_deg`, `el_deg`, `az_rad`, `el_rad`
- offset: `dx_arcsec`, `dy_arcsec`, `dx_deg`, `dy_deg`, `dx_rad`, `dy_rad`

であり、Az/El 列名に単位が無い場合だけ `angle_unit` を使います。offset 列名に単位が無い場合は `offset_unit` を使います。既定値は legacy `d_param.csv` 互換のため `arcsec` です。

## 5. normalize 後の量の意味

この版では、normalize 後の

- `delta_err_dx_arcsec`
- `delta_err_dy_arcsec`

を必ず

- **使用中モデルに加えるべき差分**

と解釈します。

したがって、absolute fit は常に

```text
new_model = used_model + delta_err
```

です。

## 6. `sign_convention`

`sign_convention` は以前の版では absolute fit の target 作成に使っていましたが、
この版では **削除** です。manifest に残っている場合は normalize が明示的にエラーで停止します。legacy な CSV 列として残っていても fit では無視されます。

## 7. 典型例

### NANTEN2

```toml
[defaults]
angle_unit = "deg"
measurement_space = "star_offset"
positive_az_moves_star = "left"
positive_el_moves_star = "down"
command_err_mode = "subtract"
model = "nanten2_actual_v1"
```

### OMU 1.85m

```toml
[defaults]
angle_unit = "deg"
measurement_space = "star_offset"
positive_az_moves_star = "left"
positive_el_moves_star = "down"
command_err_mode = "add"
model = "omu1p85m_actual_v1"
```
