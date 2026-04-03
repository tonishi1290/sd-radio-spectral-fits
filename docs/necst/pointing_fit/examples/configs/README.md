# examples/configs

このディレクトリには、manifest と TOML の最小例を置いています。

## 1. ユーザー向け標準

通常の optical pointing では、raw `dx`, `dy` は

- **星位置 - CCD中心**

として扱います。したがって manifest では通常

- `measurement_space = "star_offset"`

を使います。

## 2. NANTEN2 例

`nanten2_star_offset_manifest_example.toml` は、NANTEN2 で

- `+Az` で星が左へ動く
- `+El` で星が下へ動く
- `Az_cmd = Az_true - Err_Az`
- `El_cmd = El_true - Err_El`

を仮定した最小例です。

## 3. 実行例

```bash
pointing-fit normalize \
  --manifest examples/configs/nanten2_star_offset_manifest_example.toml \
  --output normalized.csv
```

## 4. 注意

- `sign_convention` はこの版で削除しました。manifest に残っていると normalize は明示的にエラーになります。
- 通常のユーザー向け運用では、`command_offset` や `model_delta` を主役にしない方が分かりやすくなります。


## fixed parameter 例

- `fixed_nanten2_optical.toml`
  - NANTEN2 の optical fit 向け fixed 例です。
  - `gggg = 0` と各 radio 項を固定し、まず低次の optical 項から解く出発点として使います。

- `fixed_nanten2_radio.toml`
  - NANTEN2 の radio fit 向け fixed 例です。
  - 既知の optical 項を固定しつつ、`g, gg, ggg, gggg = 0` とし、必要になったときだけ高次項を順に解放する出発点です。

推奨方針は、`g, gg, ggg, gggg` と対応する radio 項を**最初から高次まで free にしない**ことです。
低次から順に free にし、残差と相関を見て必要な場合だけ次数を上げてください。

- `apply` に渡すのは `fit_result.json` ではなく、delta fit で作られる `delta_params.toml` です。

- raw CSV の `dt_cross_id` と `cross_id` は省略できます。省略時は normalize が行番号から自動生成します。
