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
(\mathrm{Az}_{\rm bore}, \mathrm{El}_{\rm bore})=
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
(\mathrm{Az}_{\rm beam}, \mathrm{El}_{\rm beam})=
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
(\mathrm{Az}_{\rm cmd,bore}, \mathrm{El}_{\rm cmd,bore})=
\mathrm{ApplyCorrection}
\left(
\mathrm{Az}_{\rm cmd}, \mathrm{El}_{\rm cmd}, d\mathrm{lon}, d\mathrm{lat}
\right)
$$

とし、さらに beam offset を適用して

$$
(\mathrm{Az}_{\rm cmd,beam}, \mathrm{El}_{\rm cmd,beam})=
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


### 5.1 ON scan color-by 表示の目的

今回の改修版では、通常の `OBSMODE` 色分けとは別に、`--color-by` で **ON scan のみ**を色付き線分として描けます。

この機能の目的は、scan map 上で

- enc と cmd の差がどこで大きいか
- enc 側の frame 上の速度がどこで大きいか

を、`--coord-frame` で指定した座標系そのものの上で可視化することです。

`--color-by` の対象は次です。

- `none`
- `frame_dx`
- `frame_dy`
- `frame_dr`
- `frame_vx`
- `frame_vy`
- `frame_speed`

### 5.2 色付き線分を描く区間

`--color-by` を指定したとき、色付き主軌跡は **ON だけ** を対象にします。さらに、`ON scan A -> ON scan B` の遷移を除くため、線分 $k$ を描く条件は

$$
\\mathrm{same\_on\_run}(k) \\land \\mathrm{same\_scan}(k) \\land \\mathrm{finite\_metric}(k)
$$

です。

ここで

- $\\mathrm{same\_on\_run}(k)$ : 線分 $k$ の両端点が同じ ON 連続区間に属する
- $\\mathrm{same\_scan}(k)$ : 線分 $k$ の両端点が同じ scan group に属する
- $\\mathrm{finite\_metric}(k)$ : 色付けに使う metric 値が有限である

を意味します。

scan group key は次の優先順で決めます。

1. `scan_id`
2. `line_index`
3. `line_label`

つまり、`scan_id` が使えるならそれを最優先し、無い場合だけ `line_index`、さらに無い場合だけ `line_label` を使います。

### 5.3 `frame_dx`, `frame_dy`, `frame_dr` の定義

`--coord-frame` で指定した座標系上で、enc 側の座標を $(x_{\mathrm{enc}}, y_{\mathrm{enc}})$、cmd 側の座標を $(x_{\mathrm{cmd}}, y_{\mathrm{cmd}})$ とします。

ここで

- `azel` なら $(x,y)=(\mathrm{Az},\mathrm{El})$
- `radec` なら $(x,y)=(\mathrm{RA},\mathrm{Dec})$
- `galactic` なら $(x,y)=(l,b)$

です。

$x$ は経度型軸なので、そのまま引くのではなく、局所角距離に直して

$$
\Delta x_{\mathrm{local}}=
\mathrm{wrap}(x_{\mathrm{enc}}-x_{\mathrm{cmd}})
\cos\left(\frac{y_{\mathrm{enc}}+y_{\mathrm{cmd}}}{2}\right)
$$

とします。

また

$$
\Delta y = y_{\mathrm{enc}}-y_{\mathrm{cmd}}
$$

$$
\Delta r = \sqrt{\Delta x_{\mathrm{local}}^2 + \Delta y^2}
$$

です。

表示量は線分値なので、各点値を隣接 2 点で平均して、その線分に割り当てます。たとえば `frame_dx` なら

$$
\Delta x_{\mathrm{seg},k} = \frac{\Delta x_k + \Delta x_{k+1}}{2}
$$

です。同様に `frame_dy`, `frame_dr` も各点値の隣接平均を線分値として使います。

`frame_dx`, `frame_dy`, `frame_dr` の表示単位は arcsec です。

### 5.4 `frame_vx`, `frame_vy`, `frame_speed` の定義

`frame_vx`, `frame_vy`, `frame_speed` は、enc 側の frame 座標に基づく速度です。

まず隣接 2 点から instant な線分速度を

$$
\Delta t_k = t_{k+1}-t_k
$$

$$
\Delta x_{\mathrm{local},k}=
\mathrm{wrap}(x_{k+1}-x_k)
\cos\left(\frac{y_{k+1}+y_k}{2}\right)
$$

$$
\Delta y_k = y_{k+1}-y_k
$$

$$
v_{x,k} = \frac{\Delta x_{\mathrm{local},k}}{\Delta t_k},
\qquad
v_{y,k} = \frac{\Delta y_k}{\Delta t_k}
$$

$$
v_k = \sqrt{v_{x,k}^2 + v_{y,k}^2}
$$

と定義します。

ここで重要なのは、既定値 `--speed-window-points 2` は **点中心速度ではなく、隣接 2 点から作る線分速度**だということです。

したがって、ある「点」の前後点から中心差分で速度を出しているのではなく、各線分に対して速度を割り当てています。線分 $k$ の既定値は、点 $k$ と点 $k+1$ だけを使って作る 2 点差分です。

つまり、既定値では

$$
(v_{x,k}, v_{y,k}, v_k)
$$

は「点 $k$ の速度」ではなく、「点 $k$ と点 $k+1$ を結ぶ線分の代表速度」です。

表示単位は arcsec/s です。

### 5.5 `--speed-window-points` と `--speed-window-mode`

`--speed-window-points N` を 2 より大きくすると、速度量は同一 ON・同一 scan の内部だけで局所平均されます。

- `N = 2` : 隣接 2 点から作る線分速度
- `N > 2` : 同一 ON・同一 scan の内部に限った局所窓平均

平均方法は `--speed-window-mode` で選びます。

#### `secant`

窓の両端点だけを使って平均速度を作ります。

窓の始点を $p_0$、終点を $p_1$ とすると

$$
v_x = \frac{\Delta x_{\mathrm{local}}(p_0,p_1)}{t_{p_1}-t_{p_0}},
\qquad
v_y = \frac{y_{p_1}-y_{p_0}}{t_{p_1}-t_{p_0}}
$$

です。ここで $p_0, p_1$ は**点 index**です。ある線分 $i$ に対して、その近傍窓に含まれる最初の点と最後の点を使います。既定値 $N=2$ では $p_0=i$, $p_1=i+1$ です。

これは「その窓での平均速度」に最も近い定義です。

#### `component_mean`

窓内の instant な線分速度 $v_x, v_y$ をそれぞれ平均します。

$$
\bar v_x = \mathrm{mean}(v_x),
\qquad
\bar v_y = \mathrm{mean}(v_y)
$$

#### `magnitude_mean`

窓内の speed magnitude

$$
v = \sqrt{v_x^2+v_y^2}
$$

を平均します。`frame_speed` ではこの平均値そのものを使います。

一方、`frame_vx` / `frame_vy` では、実装上は

$$
\mathrm{mean}(v)
$$

に対して、それぞれ

$$
\mathrm{sign}(\mathrm{mean}(v_x)),
\qquad
\mathrm{sign}(\mathrm{mean}(v_y))
$$

の符号だけを付けた量を使います。したがって、これは純粋な成分平均ではなく、「速度大きさを成分平均符号で代表させる」挙動です。通常は `frame_speed` に使うのを勧めます。

### 5.6 signed/unsigned と color range

signed 量は

- `frame_dx`
- `frame_dy`
- `frame_vx`
- `frame_vy`

です。unsigned 量は

- `frame_dr`
- `frame_speed`

です。

さらに signed 量については

- `--color-sign-mode signed`
- `--color-sign-mode abs`

を選べます。

`abs` を選ぶと、`frame_dx`, `frame_dy`, `frame_vx`, `frame_vy` も絶対値で表示されます。つまり

$$
|\Delta x|, \quad |\Delta y|, \quad |v_x|, \quad |v_y|
$$

を色付けに使います。この場合、色範囲の扱いも unsigned と同じになります。

#### `--color-vmin`, `--color-vmax`

色範囲は `--color-vmin`, `--color-vmax` で明示できます。

- 両方未指定: 自動
- signed 量で片側だけ指定: その絶対値で 0 対称
- 両側指定: その範囲をそのまま使う

たとえば signed 量で `--color-vmax 20` とした場合は

$$
[v_{\min}, v_{\max}] = [-20, +20]
$$

です。

#### `--color-percentile`

`--color-percentile q` を指定すると、percentile ベースで自動範囲を決めます。

##### signed 量

signed 量では

$$
a = P_q(|v|)
$$

$$
[v_{\min}, v_{\max}] = [-a, +a]
$$

です。

##### unsigned 量の `central`

`--color-percentile-mode central` では、中央 $q$% を残す範囲を使います。

$$
\mathrm{tail} = \frac{100-q}{2}
$$

$$
[v_{\min}, v_{\max}] = [P_{\mathrm{tail}}(v), P_{100-\mathrm{tail}}(v)]
$$

したがって `q=99` なら

$$
[v_{\min}, v_{\max}] = [P_{0.5}(v), P_{99.5}(v)]
$$

です。

##### unsigned 量の `zero_upper`

`--color-percentile-mode zero_upper` では、下端を 0 に固定して上側だけ percentile で切ります。

$$
[v_{\min}, v_{\max}] = [0, P_q(v)]
$$

です。

#### percentile と明示範囲の優先順位

`--color-vmin`, `--color-vmax` を指定した側は、その値を優先します。未指定側だけ percentile 側の自動値を使います。

たとえば unsigned 量で

```bash
--color-percentile 99 --color-percentile-mode zero_upper --color-vmax 30
```

とすると、上端は 30 を使い、下端は `zero_upper` に従って 0 になります。

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
\mathrm{slice}=
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
- scan grouping
  - `scan_id`
  - `line_index`
  - `line_label`
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


これら 3 列は `--color-by` の ON scan 切り分けにも使います。優先順位は

1. `scan_id`
2. `line_index`
3. `line_label`

です。したがって、CSV 再描画でも RawData と同じ scan 切り分け規則を保てます。

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


### 11.4 ON scan を frame 差で色付けする

```python
from necst_v4_plot_trajectory import load_trajectory_samples, plot_trajectory

samples = load_trajectory_samples(
    rawdata="necst_otf_20260321_201209_orion-kl",
    spectrometer_config="12co-NANTEN2-multi_260331TO.conf",
    stream_names=["1LU"],
    coord_frame="radec",
    coord_source="beam",
    overlay_cmd=True,
    cmd_source="corrected_cmd",
    selection="from-first-off",
)

fig, ax, meta = plot_trajectory(
    samples,
    color_by="frame_dr",
    color_percentile=99,
    color_percentile_mode="zero_upper",
    linewidth=1.2,
)
fig.savefig("traj_frame_dr_radec.png", dpi=150, bbox_inches="tight")
```

### 11.5 `|frame_dx|` を表示する

```python
fig, ax, meta = plot_trajectory(
    samples,
    color_by="frame_dx",
    color_sign_mode="abs",
    color_percentile=99,
    color_percentile_mode="zero_upper",
)
fig.savefig("traj_abs_frame_dx.png", dpi=150, bbox_inches="tight")
```

### 11.6 平滑化した frame speed を表示する

```python
fig, ax, meta = plot_trajectory(
    samples,
    color_by="frame_speed",
    speed_window_points=7,
    speed_window_mode="secant",
    color_percentile=99,
    color_percentile_mode="central",
)
fig.savefig("traj_frame_speed_7pt.png", dpi=150, bbox_inches="tight")
```

### 11.7 CSV から color-by を再描画する

```python
samples2 = load_trajectory_samples_from_csv(
    "traj.csv",
    coord_frame="galactic",
    coord_source="beam",
    overlay_cmd=True,
    cmd_source="corrected_cmd",
    selection="all",
)

fig, ax, meta = plot_trajectory(
    samples2,
    color_by="frame_vx",
    color_sign_mode="abs",
    speed_window_points=5,
    speed_window_mode="component_mean",
    color_percentile=99,
    color_percentile_mode="zero_upper",
)
fig.savefig("traj_abs_frame_vx_from_csv.png", dpi=150, bbox_inches="tight")
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


## 12.14 ON scan を frame 差で色付けする

RA/Dec 面で、enc-cmd の位置差大きさ `frame_dr` を色付けします。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame radec \
  --coord-source beam \
  --overlay-cmd \
  --cmd-source corrected_cmd \
  --stream-name 1LU \
  --color-by frame_dr \
  --color-percentile 99 \
  --color-percentile-mode zero_upper \
  --out traj_frame_dr_radec.png \
  necst_otf_20260321_201209_orion-kl
```

## 12.15 `|frame_dx|` を色付けする

Galactic 座標で、signed ではなく絶対値 $|\Delta x|$ を使います。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame galactic \
  --coord-source beam \
  --overlay-cmd \
  --cmd-source corrected_cmd \
  --stream-name 1LU \
  --color-by frame_dx \
  --color-sign-mode abs \
  --color-percentile 99 \
  --color-percentile-mode zero_upper \
  --out traj_abs_frame_dx_gal.png \
  necst_otf_20260321_201209_orion-kl
```

## 12.16 color range を手動指定する

`frame_dx` の表示範囲を $\pm 20$ arcsec に固定します。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame radec \
  --stream-name 1LU \
  --color-by frame_dx \
  --color-vmax 20 \
  --out traj_frame_dx_pm20.png \
  necst_otf_20260321_201209_orion-kl
```

## 12.17 `frame_speed` を 7 点窓で平滑化する

この例の 7 点窓は、各線分の近傍 7 点から代表速度を作るもので、点中心速度ではありません。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --coord-source beam \
  --stream-name 1LU \
  --color-by frame_speed \
  --speed-window-points 7 \
  --speed-window-mode secant \
  --color-percentile 99 \
  --color-percentile-mode central \
  --out traj_frame_speed_7pt.png \
  necst_otf_20260321_201209_orion-kl
```

## 12.18 `frame_vx` を 5 点 `component_mean` で平滑化する

`frame_vx` / `frame_vy` を平滑化したいときは、まず `component_mean` を使うのが自然です。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  --color-by frame_vx \
  --speed-window-points 5 \
  --speed-window-mode component_mean \
  --color-percentile 99 \
  --out traj_frame_vx_5pt.png \
  necst_otf_20260321_201209_orion-kl
```

## 12.19 CSV へ scan 情報も含めて保存する

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  --export-csv traj_with_scan.csv \
  necst_otf_20260321_201209_orion-kl
```

この CSV には `scan_id`, `line_index`, `line_label` も保存されるので、後から `--from-csv` で再描画しても ON scan 境界を保ったまま color-by を使えます。

## 12.20 CSV から color-by を再描画する

```bash
necst_v4_plot_trajectory \
  --from-csv traj_with_scan.csv \
  --coord-frame galactic \
  --coord-source beam \
  --overlay-cmd \
  --cmd-source corrected_cmd \
  --color-by frame_speed \
  --speed-window-points 9 \
  --speed-window-mode secant \
  --color-percentile 99 \
  --color-percentile-mode central \
  --out traj_from_csv_frame_speed_9pt.png
```

## 12.21 `zero_upper` と `central` の使い分け

誤差大きさや速度大きさのように 0 に意味があるときは `zero_upper` が有用です。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --coord-frame radec \
  --stream-name 1LU \
  --color-by frame_dr \
  --color-percentile 99 \
  --color-percentile-mode zero_upper \
  --out traj_frame_dr_zero_upper.png \
  necst_otf_20260321_201209_orion-kl
```

一方、全体コントラストを重視したいときは `central` を使います。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --coord-frame radec \
  --stream-name 1LU \
  --color-by frame_speed \
  --color-percentile 99 \
  --color-percentile-mode central \
  --out traj_frame_speed_central.png \
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

### 13.4 `--color-by` は mode 色分けの置き換えではない

`--color-by` は通常の `OBSMODE` 色分けとは別経路です。通常軌跡全体は mode 色分けの考え方を保ちつつ、色付き主軌跡は ON scan のみを対象にします。`--overlay-cmd` を有効にした場合の cmd overlay は、`--color-by` の ON 制限とは独立です。

### 13.5 `--speed-window-points 2` の意味

`--speed-window-points` を省略したときは既定値 2 です。これは「ある点の前後点から中心差分で速度を出す」という意味ではなく、**隣接 2 点から作る 1 本の線分に速度を割り当てる**という意味です。線分 $k$ に対しては、点 $k$ と点 $k+1$ の 2 点だけを使います。

したがって、`frame_vx`, `frame_vy`, `frame_speed` は点値ではなく線分値です。

### 13.6 `magnitude_mean` の注意

`magnitude_mean` は speed magnitude の平均なので、`frame_speed` を滑らかに見る用途には向きます。

一方、`frame_vx` や `frame_vy` に使うと、実装上は speed magnitude の平均に、それぞれ $\mathrm{sign}(\mathrm{mean}(v_x))$ と $\mathrm{sign}(\mathrm{mean}(v_y))$ の符号だけを付ける挙動になります。つまり、成分そのものの平均ではありません。成分量の可視化では、通常 `secant` か `component_mean` の方が解釈しやすいです。

### 13.7 percentile 範囲と手動範囲の優先関係

`--color-vmin`, `--color-vmax` は `--color-percentile` より優先します。

- 両方未指定: percentile または全値自動
- 片側だけ指定: 指定側は固定、未指定側だけ自動
- 両側指定: 完全に手動

signed 量で片側だけ指定した場合は、その絶対値を使って 0 対称範囲を作ります。

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


### 14.4 まず scan 境界を保った CSV を作る

color-by を繰り返し試すなら、まず RawData から scan 情報付き CSV を作っておくと便利です。

```bash
necst_v4_plot_trajectory \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --selection from-first-off \
  --coord-frame azel \
  --stream-name 1LU \
  --export-csv traj_with_scan.csv \
  necst_otf_20260321_201209_orion-kl
```

### 14.5 誤差 map と速度 map を描き分ける

まず `frame_dr` で enc-cmd の差が大きい位置を見る:

```bash
necst_v4_plot_trajectory \
  --from-csv traj_with_scan.csv \
  --coord-frame radec \
  --coord-source beam \
  --overlay-cmd \
  --cmd-source corrected_cmd \
  --color-by frame_dr \
  --color-percentile 99 \
  --color-percentile-mode zero_upper \
  --out traj_frame_dr.png
```

次に `frame_speed` で scan 上の速度変化を見る:

```bash
necst_v4_plot_trajectory \
  --from-csv traj_with_scan.csv \
  --coord-frame radec \
  --coord-source beam \
  --color-by frame_speed \
  --speed-window-points 7 \
  --speed-window-mode secant \
  --color-percentile 99 \
  --color-percentile-mode central \
  --out traj_frame_speed.png
```

---

## 15. まとめ

この実装の要点は次です。

- 代表 stream の時刻だけを使って軌跡を描く
- config 優先規則を保つ
- `real_boresight`, `real_beam`, `corrected_cmd` を明示的に区別する
- `--out` と `--show` の役割を単純に保つ
- CSV による再利用経路を持つ
- `scan_id` / `line_index` / `line_label` を保持し、CSV 再描画でも ON scan 境界を再現する
- `--color-by` で ON scan のみを frame 差または frame 速度で色付けできる
- `--color-vmin` / `--color-vmax` / `--color-percentile` / `--color-percentile-mode` で color range を制御できる
- `--color-sign-mode abs` で signed metric を absolute 値表示できる
- `--speed-window-points` / `--speed-window-mode` で速度 metric を same-ON same-scan の内部だけで局所平均できる
- 矢印は既定で `segment` モードにして、見やすさを優先する

この方針により、RawData からの確認、CSV 化、座標系変更、overlay、方向確認までを 1 本のツールで扱えます。
