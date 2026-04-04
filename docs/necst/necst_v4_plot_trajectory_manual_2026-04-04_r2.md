# `necst_v4_plot_trajectory.py` 詳細説明書

## 1. 概要

`necst_v4_plot_trajectory.py` は、NECST v4 RawData から望遠鏡軌跡を描画するためのスクリプトです。

主な目的は次の 4 点です。

1. 分光配列そのものを使わず、代表 1 stream の時刻列と `OBSMODE` を使って軌跡を描く
2. converter と整合した boresight 補正・beam offset・Az/El からの座標変換を使う
3. CLI と Python 関数の両方から使えるようにする
4. 一度 CSV に書き出した後は、DB を参照せずに再描画できるようにする

本実装は、RawData から直接描画する経路と、`--export-csv` で作った CSV を `--from-csv` で再利用する経路の両方を持ちます。

---

## 2. 設計方針

### 2.1 代表 stream のみを timing に使う

複数 stream が設定されていても、軌跡表示のためには 1 本の時刻列があれば十分です。
そのため、本実装では converter と同様に `use_for_convert` の stream 群から代表 1 本を選び、その spectral table の時刻を `t_spec` として用います。

### 2.2 分光値は使わない

軌跡表示に必要なのは、主に

- 分光計取得時刻 `t_spec`
- `OBSMODE`
- encoder Az/El
- cmd/altaz Az/El
- pointing correction `dlon`, `dlat`
- beam offset
- site 情報

です。

したがって、本実装は spectrum 本体を計算には使いません。
ただし、`necstdb` の実装によっては spectral table の列 subset 読みが使えず、raw record 全体を 1 回読む fallback に入ることがあります。その場合でも spectrum の数値配列は計算には使いません。

### 2.3 config 優先規則

重要なのは、`--telescope` などの CLI オプションを**本当に明示したときだけ** converter 側へ渡すことです。

未指定時は config の `[global]` を優先します。

つまり、最終的な値の決まり方は概念的には

$$
\text{runtime value} =
\begin{cases}
\text{CLI value} & \text{if CLI で明示指定されたとき} \\
\text{config value} & \text{if config に存在するとき} \\
\text{built-in default} & \text{それ以外}
\end{cases}
$$

です。

この規則は特に

- `telescope`
- `db_namespace`
- `encoder_time_col`
- `altaz_time_col`
- `encoder_shift_sec`

で重要です。

---

## 3. 基本用語

### 3.1 時刻

- `t_spec` : 代表 spectral stream の取得時刻を unix 秒に直したもの

### 3.2 生の pointing

- `enc_az_deg`, `enc_el_deg` : encoder の Az, El
- `cmd_az_deg`, `cmd_el_deg` : cmd/altaz の Az, El
- `dlon_deg`, `dlat_deg` : pointing correction

### 3.3 補正後 boresight

encoder 系から補正後 boresight を作るときは、`azel_correction_apply` の規約に従います。

概念的には

$$
(\mathrm{Az}_{\rm bore}, \mathrm{El}_{\rm bore})
=
\mathrm{ApplyCorrection}
\left(
\mathrm{Az}_{\rm enc}, \mathrm{El}_{\rm enc}, d\mathrm{lon}, d\mathrm{lat}
\right)
$$

です。

本実装では CSV 列名として

- `real_boresight_az_deg`
- `real_boresight_el_deg`

を使います。

### 3.4 beam-center

補正後 boresight に beam offset を適用した最終 beam-center は

$$
(\mathrm{Az}_{\rm beam}, \mathrm{El}_{\rm beam})
=
\mathrm{ApplyBeamOffset}
\left(
\mathrm{Az}_{\rm bore}, \mathrm{El}_{\rm bore}
\right)
$$

です。

本実装では CSV 列名として

- `real_beam_az_deg`
- `real_beam_el_deg`

を使います。

### 3.5 corrected cmd

cmd/altaz 側にも同じ規約で `dlon`, `dlat` を適用した量を保存します。

まず corrected cmd boresight を

$$
(\mathrm{Az}_{\rm cmd,bore}, \mathrm{El}_{\rm cmd,bore})
=
\mathrm{ApplyCorrection}
\left(
\mathrm{Az}_{\rm cmd}, \mathrm{El}_{\rm cmd}, d\mathrm{lon}, d\mathrm{lat}
\right)
$$

とし、さらに beam offset を適用して

$$
(\mathrm{Az}_{\rm cmd,beam}, \mathrm{El}_{\rm cmd,beam})
=
\mathrm{ApplyBeamOffset}
\left(
\mathrm{Az}_{\rm cmd,bore}, \mathrm{El}_{\rm cmd,bore}
\right)
$$

を作ります。

CSV 列名は

- `corrected_cmd_boresight_az_deg`
- `corrected_cmd_boresight_el_deg`
- `corrected_cmd_beam_az_deg`
- `corrected_cmd_beam_el_deg`

です。

---

## 4. `coord-source` の意味

`--coord-source` は、描画に使う Az/El の定義を指定します。

正式名は次の 5 つです。

- `beam`
- `boresight`
- `raw_encoder`
- `raw_altaz`
- `corrected_cmd`

### 4.1 正式名の意味

- `beam` : 補正後 encoder boresight に beam offset を載せた最終 beam-center
- `boresight` : 補正後 encoder boresight
- `raw_encoder` : 生の encoder Az/El に beam offset を載せたもの
- `raw_altaz` : 生の cmd/altaz Az/El に beam offset を載せたもの
- `corrected_cmd` : 補正後 cmd boresight に beam offset を載せたもの

### 4.2 alias

CLI と CSV 再読込では、次の alias も使えます。

- `true` -> `beam`
- `encoder` -> `raw_encoder`
- `altaz` -> `raw_altaz`
- `cmd` -> `corrected_cmd`

---

## 5. `coord-frame` の意味

`--coord-frame` は、描画面をどの座標系で表すかを指定します。

- `azel` : Az, El
- `radec` : RA, Dec
- `galactic` : Galactic longitude, latitude

RawData から描画するときは converter 側の Az/El -> RA/Dec 変換を使います。CSV 再描画では、CSV に保存された Az/El と site, `t_spec` を使って local に変換します。

---

## 6. `OBSMODE` と色

軌跡色は `OBSMODE` に基づきます。

- `ON` : 赤
- `OFF` : 青
- `HOT` : 緑
- その他 : 灰

大文字小文字の違いは吸収されます。

---

## 7. OFF block の定義と selection

### 7.1 OFF block

`mode == "OFF"` の連続区間を 1 つの OFF block と定義します。

たとえば、時系列の mode が

```text
ON ON OFF OFF OFF ON HOT OFF OFF ON
```

なら、OFF block は 2 個です。

- block 1 : index `[2, 5)`
- block 2 : index `[7, 9)`

### 7.2 selection

`--selection` には次を指定できます。

- `all`
- `from-first-off`
- `off-block-range`

#### `all`
全サンプルを使います。

#### `from-first-off`
最初の OFF block の先頭から最後までを使います。

#### `off-block-range`
`--off-block-start`, `--off-block-end` を 1-based index として受け取り、

- 第 `start` OFF block の先頭
- 第 `end` OFF block の末尾

で挟まれた連続区間を使います。

内部的には

$$
\mathrm{slice}
=
[s_{\rm start}, e_{\rm end})
$$

なので、**OFF だけではなく、その間の ON/HOT/その他もすべて含む**のが重要です。

---

## 8. 出力規則

`--out` と `--show` の動作は次です。

- `--out` なし、`--show` なし -> 画面表示のみ
- `--out` あり、`--show` なし -> ファイル出力のみ
- `--out` なし、`--show` あり -> 画面表示のみ
- `--out` あり、`--show` あり -> 画面表示 + ファイル出力

つまり、論理的には

$$
\mathrm{do\_show} = (\texttt{--show}) \lor (\texttt{--out が無い})
$$

$$
\mathrm{do\_save} = (\texttt{--out がある})
$$

です。

保存形式は `--out` の拡張子で決まります。たとえば

- `traj.png`
- `traj.pdf`
- `traj.svg`

のように使います。

---

## 9. 方向矢印

### 9.1 有効化

方向矢印は `--show-direction-arrows` で有効になります。

### 9.2 `arrow-mode`

- `segment` : 既定
- `samples`

#### `segment` モード
このモードが既定です。各 mode 連続区間について

- 十分大きく動いているなら、その区間に 1 本
- mode 境界で十分動いているなら、境界付近にも 1 本
- その座標系でほとんど動いていない区間には描かない

という方針です。

「ON の 1 セグメントに 1 本、ON -> OFF 移動に 1 本、でも明確に動いていないときは描かない」という用途を狙っています。

#### `samples` モード
時系列に沿って N サンプルごとに 1 本描きます。N は `--arrow-every` で決まります。

### 9.3 矢印の長さ

矢印は速度ベクトルをそのまま描くのではなく、**方向表示用の固定可視長**で描きます。
図の長辺スパンを $S$ とすると、描画長は

$$
L_{\rm draw} = f S
$$

で、既定では `f = 0.015` です。

したがって `--arrow-length-frac` は「矢印長を図サイズに対してどの程度にするか」を表し、速度の大きさそのものを表しません。

### 9.4 関連オプション

- `--show-direction-arrows`
- `--arrow-mode {segment,samples}`
- `--arrow-every N`
- `--arrow-alpha A`
- `--arrow-width W`
- `--arrow-length-frac F`

---

## 10. CSV 出力と CSV 再描画

### 10.1 `--export-csv`

描画に使った軌跡パラメータを CSV に保存します。

- `--export-csv path.csv` : 指定パスに保存
- `--export-csv` だけ : 既定名に保存

既定名は

- RawData 入力時 : `<rawbasename>_trajectory.csv`
- `--from-csv` 入力時 : `<csvbasename>_plot.csv`

です。

### 10.2 `--from-csv`

以前に `--export-csv` で保存した CSV を読み、DB を参照せずに再描画します。

これにより、

1. 一度だけ RawData から CSV を作る
2. 以後は CSV だけで軽く再描画する

という使い方ができます。

### 10.3 CSV 主要列

CSV には次のような主要列が入ります。

- index / time
  - `sample_index`
  - `unix_time`
  - `iso_time`
- mode
  - `mode`
  - `mode_norm`
  - `off_block_index`
- site
  - `site_lat_deg`
  - `site_lon_deg`
  - `site_elev_m`
- raw pointing
  - `enc_az_deg`, `enc_el_deg`
  - `cmd_az_deg`, `cmd_el_deg`
  - `dlon_deg`, `dlat_deg`
- corrected / beam
  - `real_boresight_az_deg`, `real_boresight_el_deg`
  - `real_beam_az_deg`, `real_beam_el_deg`
- raw alternative sources
  - `raw_encoder_beam_az_deg`, `raw_encoder_beam_el_deg`
  - `raw_altaz_beam_az_deg`, `raw_altaz_beam_el_deg`
- corrected cmd
  - `corrected_cmd_boresight_az_deg`, `corrected_cmd_boresight_el_deg`
  - `corrected_cmd_beam_az_deg`, `corrected_cmd_beam_el_deg`
- beam / pointing error diagnostics
  - `beam_dx_arcsec`, `beam_dy_arcsec`, `beam_rot_deg`
  - `pe_x_arcsec`, `pe_y_arcsec`, `pe_r_arcsec`

---

## 11. Python API

主に次の関数を使います。

- `load_trajectory_samples(...)`
- `load_trajectory_samples_from_csv(...)`
- `export_trajectory_csv(samples, csv_path)`
- `plot_trajectory(samples, ...)`

### 11.1 RawData から読む

```python
from necst_v4_plot_trajectory import load_trajectory_samples, plot_trajectory

samples = load_trajectory_samples(
    rawdata="necst_otf_20260321_201209_orion-kl",
    spectrometer_config="12co-NANTEN2-multi_260331TO.conf",
    stream_names=["1LU"],
    coord_frame="azel",
    coord_source="beam",
    overlay_cmd=True,
    cmd_source="corrected_cmd",
    selection="from-first-off",
)

fig, ax, meta = plot_trajectory(
    samples,
    show_direction_arrows=True,
    arrow_mode="segment",
)
fig.savefig("traj.png", dpi=150, bbox_inches="tight")
```

### 11.2 CSV を書く

```python
from necst_v4_plot_trajectory import export_trajectory_csv

export_trajectory_csv(samples, "traj.csv")
```

### 11.3 CSV から再描画

```python
from necst_v4_plot_trajectory import load_trajectory_samples_from_csv, plot_trajectory

samples = load_trajectory_samples_from_csv(
    "traj.csv",
    coord_frame="radec",
    coord_source="corrected_cmd",
    overlay_cmd=False,
    selection="all",
)

fig, ax, meta = plot_trajectory(samples, show_direction_arrows=True)
fig.savefig("traj_radec.pdf", dpi=150, bbox_inches="tight")
```

---

## 12. CLI 使用例

## 12.1 何も保存せず画面表示

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --stream-name 1LU \
  --coord-frame azel \
  necst_otf_20260321_201209_orion-kl
```

## 12.2 `from-first-off` で corrected cmd を描く

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --coord-source corrected_cmd \
  --stream-name 1LU \
  necst_otf_20260321_201209_orion-kl
```

## 12.3 `off-block-range` を使う

第 2 OFF block の先頭から第 5 OFF block の末尾までを描きます。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection off-block-range \
  --off-block-start 2 \
  --off-block-end 5 \
  --coord-frame azel \
  --stream-name 1LU \
  necst_otf_20260321_201209_orion-kl
```

## 12.4 cmd overlay 付き

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --coord-source beam \
  --overlay-cmd \
  --cmd-source corrected_cmd \
  --cmd-alpha 0.5 \
  --stream-name 1LU \
  necst_otf_20260321_201209_orion-kl
```

## 12.5 RA/Dec で表示

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --coord-frame radec \
  --stream-name 1LU \
  necst_otf_20260321_201209_orion-kl
```

## 12.6 Galactic で表示

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --coord-frame galactic \
  --stream-name 1LU \
  necst_otf_20260321_201209_orion-kl
```

## 12.7 `segment` 矢印

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --coord-source corrected_cmd \
  --stream-name 1LU \
  --show-direction-arrows \
  necst_otf_20260321_201209_orion-kl
```

## 12.8 `samples` 矢印

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  --show-direction-arrows \
  --arrow-mode samples \
  --arrow-every 10 \
  necst_otf_20260321_201209_orion-kl
```

## 12.9 PDF 保存のみ

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  --out traj.pdf \
  necst_otf_20260321_201209_orion-kl
```

## 12.10 画面表示しながら PNG 保存

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  --out traj.png \
  --show \
  necst_otf_20260321_201209_orion-kl
```

## 12.11 CSV を保存

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  --export-csv traj.csv \
  necst_otf_20260321_201209_orion-kl
```

## 12.12 CSV から再描画

```bash
necst_v4_plot_trajectory \
  --from-csv traj.csv \
  --coord-frame radec \
  --coord-source corrected_cmd \
  --show-direction-arrows \
  --show
```

## 12.13 mode summary を表示

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --print-mode-summary \
  necst_otf_20260321_201209_orion-kl
```

---

## 13. 実装上の注意

### 13.1 `coord-source corrected_cmd`

`corrected_cmd` は「cmd に correction を同じ規約で適用し、さらに beam offset を載せたもの」です。したがって、encoder base の `beam` とは一致しません。

### 13.2 `from-csv` の site

CSV 再描画では DB を見ないため、CSV に保存された

- `site_lat_deg`
- `site_lon_deg`
- `site_elev_m`

を使って座標変換します。

### 13.3 `segment` 矢印の見方

`segment` 矢印は「局所的な向きのラベル」です。速度ベクトルや移動量をそのまま表すものではなく、見やすさのために長さを固定しています。

---

## 14. 推奨ワークフロー

### 14.1 最初に RawData から確認

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  necst_otf_20260321_201209_orion-kl
```

### 14.2 CSV を書いて再利用

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  --export-csv traj.csv \
  necst_otf_20260321_201209_orion-kl
```

その後は

```bash
necst_v4_plot_trajectory --from-csv traj.csv --coord-frame radec --show
```

のように軽く再描画します。

### 14.3 overlay と矢印を足す

```bash
necst_v4_plot_trajectory \
  --from-csv traj.csv \
  --coord-frame radec \
  --coord-source beam \
  --overlay-cmd \
  --cmd-source corrected_cmd \
  --show-direction-arrows \
  --show
```

---

## 15. まとめ

この実装の要点は次です。

- 代表 stream の時刻だけを使って軌跡を描く
- config 優先規則を保つ
- `real_boresight`, `real_beam`, `corrected_cmd` を明示的に区別する
- `--out` と `--show` の役割を単純に保つ
- CSV による再利用経路を持つ
- 矢印は既定で `segment` モードにして、見やすさを優先する

この方針により、RawData からの確認、CSV 化、座標系変更、overlay、方向確認までを 1 本のツールで扱えます。
