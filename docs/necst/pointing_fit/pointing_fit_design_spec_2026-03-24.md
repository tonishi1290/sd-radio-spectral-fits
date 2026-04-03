# Optical Pointing CLI 設計整理メモ

作成日: 2026-03-24
対象: `pointing_fit`

## 1. この文書の目的

本メモは、pointing fit 解析における `dx`, `dy` の意味、`command_err_mode` の役割、`sign_convention` の位置付けを整理し、今後の仕様を混乱の少ない形に固めるための設計文書である。

今回の議論で確認できた重要点は次の通りである。

1. ユーザーが自然に与える量は、通常は
   - 星位置 - CCD 中心
   である。
2. その raw 値が画像座標なのか、すでに Az/El 正方向へ解釈済みなのかを明確に分けないと、符号議論が破綻しやすい。
3. 現行実装では `sign_convention` が absolute fit の target 作成に効いているが、仕様を明確にすればこれは原理的には不要である。
4. 実務上は、normalize 後の値を常に
   - 使用中モデルに加えるべき器差差分
   と定義すると、NANTEN2 / 1.85m の違いは `command_err_mode` だけで整理できる。

---

## 2. 現行実装で確認したこと

### 2.1 normalize 側

現行 `normalize.py` では、`measurement_space`, `positive_az_moves_star`, `positive_el_moves_star`, `command_err_mode` から `dx_to_model_sign`, `dy_to_model_sign` を作り、入力 `dx`, `dy` にその符号を掛けて `measured_dx_arcsec`, `measured_dy_arcsec` を作っている。

要するに、normalize の役目は

- raw の `dx`, `dy`
- -> fit が使う内部量

への変換である。

### 2.2 fit 側

現行 `fit.py` では

- `fit_target = "delta"` のときは `target = measured`
- `fit_target = "absolute"` のときは `sign_convention` により
  - `measured + applied_model`
  - または `applied_model - measured`

を選んでいる。

したがって、現行実装における `sign_convention` は normalize ではなく absolute fit の target 作成にだけ効いている。

---

## 3. 今回合意した基本設計

## 3.1 ユーザー向けの raw 入力

ユーザーが与える `dx`, `dy` は、原則として

- 星位置 - CCD 中心

とする。

これは観測者にとって最も自然な定義である。

### ただし重要な分岐

同じ「星位置 - CCD 中心」と言っても、実際の入力値には二通りある。

#### A. 画像座標型

raw `dx`, `dy` が本当に画像上の左右上下で与えられている場合。

このときは

- `measurement_space = "star_offset"`
- `positive_az_moves_star`
- `positive_el_moves_star`

が必要である。

#### B. Az/El 正方向解釈済み型

raw `dx`, `dy` が、すでに

- `+dx` = +Az に動かすべき量
- `+dy` = +El に動かすべき量

として定義されている場合。

このときは本質的には `command_offset` 型であり、`positive_*` は不要である。

---

## 4. 変数の定義を固定する

今後の議論では、以下の記号を固定して使う。

### 4.1 raw 入力

- `dx_raw`, `dy_raw`
  - CSV に入っている元の測定値
  - 単位は最終的に arcsec に統一して扱う

### 4.2 指令補正量

- `Δcmd_x`, `Δcmd_y`
  - 星を CCD 中心へ戻すために、望遠鏡へ追加で与えるべき command 補正
  - 正方向は `+Az`, `+El`

raw がすでに Az/El 正方向で解釈済みなら

- `Δcmd_x = dx_raw`
- `Δcmd_y = dy_raw`

である。

### 4.3 器差差分

- `ΔErr_x`, `ΔErr_y`
  - 現在使っている器差モデルに追加すべき差分
  - 器差モデル `Err` と同じ符号規約で表す

### 4.4 使用中モデル値

- `Err_used_x`, `Err_used_y`
  - 現在の `used_param` から計算されたモデル補正値

### 4.5 新しいモデル値

- `Err_new_x`, `Err_new_y`
  - フィット後に得たい新しいモデル補正値

---

## 5. command_err_mode の意味

`command_err_mode` は、観測プログラムが器差 `Err` を指令へどう適用しているかを表す。

### 5.1 subtract

定義:

- `Az_cmd = Az_true - Err_x`
- `El_cmd = El_true - Err_y`

したがって、必要な指令補正 `Δcmd` と器差差分 `ΔErr` の関係は

- `ΔErr_x = -Δcmd_x`
- `ΔErr_y = -Δcmd_y`

である。

### 5.2 add

定義:

- `Az_cmd = Az_true + Err_x`
- `El_cmd = El_true + Err_y`

したがって

- `ΔErr_x = +Δcmd_x`
- `ΔErr_y = +Δcmd_y`

である。

### 5.3 まとめ

- NANTEN2 は `subtract`
- 1.85m は `add`

で整理する。

---

## 6. normalize 後の内部量の意味を固定する

今後は normalize 後の内部量を必ず

- `delta_err_dx_arcsec = ΔErr_x`
- `delta_err_dy_arcsec = ΔErr_y`

と定義する。

すなわち、normalize 後の量は常に

- 使用中モデルに加えるべき器差差分

である。

この定義にすると、絶対値 fit でも差分 fit でも解釈が一貫する。

---

## 7. absolute fit / delta fit の定義

### 7.1 delta fit

フィット対象はそのまま器差差分である。

- `target_dx = delta_err_dx_arcsec`
- `target_dy = delta_err_dy_arcsec`

### 7.2 absolute fit

フィット対象は新しいモデル値である。

- `target_dx = used_model_dx_arcsec + delta_err_dx_arcsec`
- `target_dy = used_model_dy_arcsec + delta_err_dy_arcsec`

これは式で書けば

- `Err_new = Err_used + ΔErr`

である。

---

## 8. sign_convention についての結論

### 8.1 本質

仕様を明確にした今回の整理では、`sign_convention` は原理的に不要である。

理由は、normalize 前段で

- raw が何を意味するか
- 画像座標 -> command 座標
- command 座標 -> 器差差分

が一意に定まり、その結果 normalize 後の量を

- 使用中モデルに加えるべき差分

と固定できるからである。

### 8.2 現行実装で存在している理由

現行実装では absolute fit で

- `measured + applied_model`
- `applied_model - measured`

の二択が残っているため、`sign_convention` が必要になっている。

しかしこれは、同じ物理量に二重の自由度を与えており、混乱の原因である。

### 8.3 今後の方針

- `sign_convention` は外部仕様から外す
- 少なくとも通常運用ではユーザーに指定させない
- absolute fit は常に
  - `target = used_model + delta_err`
  に固定する

---

## 9. 推奨 manifest

## 9.1 NANTEN2

```toml
[defaults]
angle_unit = "deg"
measurement_space = "star_offset"
positive_az_moves_star = "left"
positive_el_moves_star = "down"
command_err_mode = "subtract"
model = "nanten2_actual_v1"

[[datasets]]
dataset_id = "run1"
input = "d_param.csv"
used_param = "pointing_param.toml"
```

## 9.2 1.85m

```toml
[defaults]
angle_unit = "deg"
measurement_space = "star_offset"
positive_az_moves_star = "left"
positive_el_moves_star = "down"
command_err_mode = "add"
model = "omu1p85m_actual_v1"

[[datasets]]
dataset_id = "run1"
input = "d_param.csv"
used_param = "pointing_param.toml"
```

### 注意

`positive_az_moves_star`, `positive_el_moves_star` は、実際の CCD 画像上で

- `+Az` command で星がどちらへ動くか
- `+El` command で星がどちらへ動くか

を実機ごとに決める必要がある。

ここで例として `left`, `down` を書いているが、装置構成が異なれば変更されうる。

---

## 10. ユーザー向け仕様と内部仕様を分ける

### 10.1 ユーザー向け正式入口

原則として、公開仕様では

- `measurement_space = "star_offset"`

を正式入口とする。

この方がユーザーにとって自然である。

### 10.2 内部仕様

内部では早い段階で

- `star_offset`
- -> `command_offset`
- -> `delta_err`

へ変換してよい。

`command_offset` は内部概念として残してよいが、通常ユーザー向けドキュメントでは主役にしない。

### 10.3 例外

もし upstream がすでに Az/El 正方向の `dx`, `dy` を出している場合は、内部上級者向けオプションとして `command_offset` を残してもよい。

ただしその場合は

- `positive_az_moves_star`
- `positive_el_moves_star`

を禁止または無視する方が安全である。

---

## 11. 変数名・列名の整理案

現行名は、`measured`, `applied_model`, `target` などが抽象的で、符号議論を難しくしている。

以下のように意味の見える名前へ変更すると混乱が減る。

### 11.1 設定名

現行 -> 提案

- `measurement_space` -> `raw_dxdy_space`
- `star_offset` -> `star_center_offset`
- `command_offset` -> `command_delta`
- `model_delta` -> `err_delta`
- `positive_az_moves_star` -> `star_motion_for_positive_az_cmd`
- `positive_el_moves_star` -> `star_motion_for_positive_el_cmd`
- `command_err_mode` -> `command_model_relation`

値も

- `subtract` / `add`

より

- `cmd_eq_true_minus_err`
- `cmd_eq_true_plus_err`

の方が誤解が少ない。

### 11.2 DataFrame 列名

現行 -> 提案

- `measured_dx_arcsec` -> `delta_err_dx_arcsec`
- `measured_dy_arcsec` -> `delta_err_dy_arcsec`
- `applied_model_dx_arcsec` -> `used_model_dx_arcsec`
- `applied_model_dy_arcsec` -> `used_model_dy_arcsec`
- `target_dx_arcsec` -> `fit_target_dx_arcsec`
- `target_dy_arcsec` -> `fit_target_dy_arcsec`
- `dx_to_model_sign` -> `dx_raw_to_delta_err_sign`
- `dy_to_model_sign` -> `dy_raw_to_delta_err_sign`

---

## 12. 実装変更方針

### 12.1 normalize.py

方針:

1. normalize 後の量を常に `delta_err_*` の意味に固定する。
2. `command_offset` のときは `positive_*` を禁止する。
3. `star_offset` のときは `positive_*` を必須にする。
4. 旧 `measured_*` 列は移行期間だけ残すか、alias とする。

### 12.2 fit.py

方針:

1. `sign_convention` への依存をやめる。
2. delta fit は
   - `target = delta_err`
3. absolute fit は
   - `target = used_model + delta_err`
4. `model_minus_measured` 分岐は削除または deprecated とする。

### 12.3 CLI / manifest

方針:

1. `sign_convention` を deprecated にする。
2. 指定されても warning を出して無視する。
3. ドキュメントから `sign_convention` を外す。

---

## 13. 移行戦略

安全に移行するため、以下の段階を推奨する。

### 段階1

- 新しい内部意味を実装する
- 旧列名は互換のため残す
- `sign_convention` 使用時に warning を出す

### 段階2

- 出力 CSV / レポート / ドキュメントを新名称へ切り替える
- `sign_convention` をドキュメントから削除する

### 段階3

- `sign_convention` の分岐そのものを削除する

---

## 14. 最終結論

今回の設計方針は、次の一文にまとめられる。

**ユーザーは常に「星が CCD 中心からどちらへどれだけずれているか」を与える。**
**内部ではそれを「使用中モデルに加えるべき器差差分」へ正規化する。**
**absolute fit は常に `used_model + delta_err` とする。**

この方針にすると

- NANTEN2 / 1.85m の違いは `command_err_mode` で整理できる
- `sign_convention` のような曖昧な自由度を消せる
- ユーザー視点と内部処理の対応が明確になる

以上を、今後の pointing fit CLI の正式設計方針とする。
