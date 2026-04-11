# `calibrate.py` state-wise QC 詳細説明書

## 0. この文書の位置づけ

この文書は、添付いただいた最新版 `qc_spur_on_prescreen_detailed_manual_ja_2026-04-10_recommended.md` をベースに、**現在の実装済みコード** `calibrate_updated_2026-04-10_statewise.py` に合わせて更新した詳細説明書です。

本書で対象とする内容は次の 4 点です。

1. fixed bad channel の扱い
2. HOT / OFF の state 向け narrow spur QC
3. ON の state 向け高速化付き spur QC
4. 推奨設定と運用例

この版では、以前の文書に含まれていた

- ON pre-screen は提案仕様で未実装
- ON / HOT / OFF にほぼ同じ狭帯域 spur 本判定を当てる

という記述は採用しません。ここでは、**現在の実装に一致する内容だけ**を書きます。

---

## 1. 基本方針

本実装の基本方針は、**共通前処理 + state 別判定器** です。

すなわち、

- 全 state 共通で fixed bad channel を処理する
- HOT / OFF には参照信号向けの row QC を使う
- ON には天体線と spur を分けやすい判定器を使う
- 高速化は主に ON に入れる

という構成です。

記号は次を使います。

- state: $x \in \{\mathrm{ON}, \mathrm{OFF}, \mathrm{HOT}\}$
- raw channel index: $k$
- row index: $r$
- raw spectrum: $S_x(k, r)$
- science window: $W_{\mathrm{sci}}$
- fixed bad channel 集合: $K_{\mathrm{bad}}$
- 有効チャネル集合:

$$
V = W_{\mathrm{sci}} \setminus K_{\mathrm{bad}}
$$

- 有効チャネル数:

$$
N = |V|
$$

- rolling median などの滑らかな成分:

$$
B_x(k, r)
$$

- residual:

$$
R_x(k, r) = S_x(k, r) - B_x(k, r)
$$

- residual の robust 尺度:

$$
\sigma_{R,x}(r) = 1.4826 \, \mathrm{MAD}\left(R_x(k, r),\, k \in V\right)
$$

- 規格化 residual:

$$
z_x(k, r) = \frac{|R_x(k, r)|}{\sigma_{R,x}(r)}
$$

---

## 2. 実装済みの公開パラメータ

現在の `make_tastar_dumps(...)` / `run_tastar_calibration(...)` では、QC 関連として次を受け付けます。

### 2.1 既存 shape QC

- `qc_sigma`
- `qc_apply_flagrow`

### 2.2 狭帯域 spur QC

- `qc_spur_sigma`
- `qc_spur_window_channels`
- `qc_spur_max_consecutive_channels`
- `qc_spur_consider_on`
- `qc_spur_on_apply_flagrow`

### 2.3 fixed bad channel

- `bad_channels`
- `bad_channel_ranges`
- `bad_channel_key_columns`
- `bad_channel_map`
- `bad_channel_policy`
- `bad_channel_interp_max_consecutive_channels`
- `bad_channel_interp_min_valid_neighbors`

---

## 3. 共通前処理: fixed bad channel

### 3.1 基本方針

fixed bad channel は row 異常ではなく、**channel-level の恒常的異常**です。

したがって fixed bad channel は `FLAGROW` では処理せず、channel mask として扱います。

### 3.2 指定方法

#### `bad_channels`

全 group 共通の raw channel index を指定します。

```python
bad_channels=[1536]
```

#### `bad_channel_ranges`

全 group 共通の range を指定します。range は **両端を含む closed interval** です。

```python
bad_channel_ranges=[(1534, 1538)]
```

#### `bad_channel_map`

group ごとに異なる fixed bad channel を指定します。key は `bad_channel_key_columns` の順です。

```python
bad_channel_key_columns=("FDNUM", "IFNUM", "PLNUM")

bad_channel_map = {
    (2, 0, 0): [1536],
    (3, 0, 1): [(1534, 1538)],
}
```

### 3.3 インデックスの解釈

ここで指定する channel index は、**raw の 0-based channel index** です。

`ch_range` や `vlsrk_range_kms` で切り出した後の local index ではありません。

たとえば raw 1536 ch を bad にしたいなら、`1536` を指定します。`ch_range=(1000, 2000)` を使っているからといって `536` を指定してはいけません。

### 3.4 QC 用配列と較正用配列は分ける

実装上、fixed bad channel を処理した配列は 3 つの役割に分けます。

#### spur QC 用配列

狭帯域 spur QC では、bad channel は**常に除外**します。つまり、該当チャネルは NaN にし、`finite` なチャネルだけを見ます。

#### shape QC 用配列

shape QC では、fixed bad channel 自体で row が落ちないように、**局所線形補間済みの専用配列**を使います。

#### 較正計算用配列

較正に実際に使う配列では `bad_channel_policy` に従います。

- `"nan"`: bad channel は NaN のまま
- `"interp_linear"`: 条件を満たす短い区間だけ線形補間

### 3.5 補間の定義

連続 bad 区間を $[k_1, k_2]$ とし、その外側に有効チャネル $k_L < k_1$ と $k_R > k_2$ があるとき、`interp_linear` では

$$
S'_x(k, r)
=
S_x(k_L, r)
+
\frac{k-k_L}{k_R-k_L}
\left(S_x(k_R, r)-S_x(k_L, r)\right)
$$

を使います。

### 3.6 補間条件

`interp_linear` では次を満たすときだけ補間します。

- edge ではない
- 両側に有効チャネルがある
- 連続 bad 区間長が `bad_channel_interp_max_consecutive_channels` 以下
- `bad_channel_interp_min_valid_neighbors` を満たす

満たさない場合は NaN のまま残します。

### 3.7 fixed bad channel 自体では `FLAGROW` を立てない

fixed bad channel は row 異常ではないため、**それ自体では `FLAGROW` を立てません**。

`FLAGROW` は、shape QC や spur QC によって「その row 全体を使わない」と判断されたときだけ立ちます。

---

## 4. bad scan と QC の順序

ここは重要です。

既存 `FLAGROW` によって bad とされている row は、**QC の判定対象から先に外れます**。

ただし実装は「物理的に先に配列を短く切る」方式ではなく、まず

$$
\text{candidate\_mask} = \lnot \text{input\_bad}
$$

を作り、QC 計算では candidate row だけを見る方式です。

その後、

- 既存 `FLAGROW` bad
- shape QC による auto bad
- spur QC による auto bad

をまとめて `FLAGROW` に反映し、最後に残った row だけが平均・内挿・較正に進みます。

したがって、**bad scan を含んだまま他の row を判定する**という順序にはなっていません。

---

## 5. HOT / OFF の QC

### 5.1 基本方針

HOT / OFF には本物の天体線が入らないため、**row 自動 flag の主戦場**にしてよいです。

HOT / OFF の最終判定は

$$
\text{auto\_bad}
=
\text{shape\_bad}
\lor
\text{narrow\_spur\_bad}
$$

です。

### 5.2 broad な異常: shape QC

これは既存の `qc_sigma` 系です。`grouped_affine_shape_qc(...)` によって、同じ group の HOT / OFF の中で broad な shape 外れ値を見ます。

fixed bad channel の影響を減らすため、shape QC では **shape 用補間配列**を入力に使います。

### 5.3 狭帯域異常: narrow spur QC

狭帯域 spur QC では、各 row について rolling median ベースの residual を作り、1 ch から `qc_spur_max_consecutive_channels` ch までの narrow run を探します。

各 narrow run $C$ の local score を

$$
T(C) = \max_{k \in C} z_x(k, r)
$$

とし、その row の代表 score を

$$
T_{\max}(r) = \max_C T(C)
$$

とします。

### 5.4 window サイズ補正

science window をそのまま使うと、チャネル数が多いほど「どこかで 1 回引っかかる」確率が上がります。そこで HOT / OFF の narrow spur QC では、**trial 数補正付き threshold** を使います。

trial 数は単なる $N$ ではなく、run 長を 1 から $L_{\max}$ まで試した総数として

$$
M(N, L_{\max}) = \sum_{\ell=1}^{L_{\max}} (N - \ell + 1)
$$

と定義します。

実装では、`qc_spur_sigma` を **基準有効チャネル数**

$$
N_{\mathrm{ref}} = 512
$$

に対する local threshold とみなし、実際の $N$ に対して有効 threshold を内部で計算します。

このため、512 ch では旧来の固定 threshold とほぼ同じで、数千 ch では少しだけ threshold が上がります。

### 5.5 実装上の意味

- 1000 ch 以下では旧版とかなり近い
- 数千 ch で science window が広いだけで HOT / OFF が落ちやすくなる傾向を抑える

という役割です。

---

## 6. ON の QC

### 6.1 基本方針

ON では本物の天体線が入るため、HOT / OFF と同じ row 自動 flag を主判定にしてはいけません。

現在の実装では ON について、

1. pre-screen
2. narrow spur 本判定
3. 近傍 ON row との孤立性判定

の 3 段階に分けています。

### 6.2 ON pre-screen

各 ON row のスペクトルを $S(i)$ とし、隣接チャネル差分を

$$
D(i) = S(i+1) - S(i)
$$

とします。

差分系列の robust 尺度を

$$
\sigma_D = 1.4826 \, \mathrm{MAD}(D)
$$

とし、差分の最大規格化振幅を

$$
Q_D = \frac{\max_i |D(i)|}{\sigma_D}
$$

と定義します。

実装では内部固定値

$$
T_D = 8
$$

を使い、

$$
Q_D > T_D
$$

の row だけを本判定へ送ります。

ここで重要なのは、**$Q_D$ が大きいから bad ではない**という点です。$Q_D$ はあくまで高速化のための candidate 選別です。

### 6.3 ON narrow spur 本判定

pre-screen を通過した ON row に対してだけ、HOT / OFF と同じ rolling median residual ベースの narrow spur 本判定を行います。

ただし ON では、それだけで最終 auto bad にはしません。

### 6.4 ON の孤立性判定

ON では本物の狭い線もあり得るため、candidate row の周辺 ON row に同じ構造が再現するかを見ます。

candidate row の narrow spur score を $T_{\max}(r)$、近傍 ON row の対応 score を $T_{\mathrm{nbr}}(r)$ とし、孤立性指標を

$$
Q_{\mathrm{iso}}(r) = \frac{T_{\max}(r)}{\max(T_{\mathrm{nbr}}(r), \epsilon)}
$$

で定義します。

実装では内部固定値

$$
Q_{\mathrm{iso,min}} = 3.0
$$

を使い、これを超えたときに「近傍 ON に再現しない孤立した狭帯域異常」とみなします。

### 6.5 ON の `FLAGROW` 反映

ON は天体線を含むため、自動 `FLAGROW` は既定で off です。

- `qc_spur_consider_on=False` なら ON QC 自体を行いません
- `qc_spur_consider_on=True` なら ON QC を実行します
- `qc_spur_on_apply_flagrow=False` が既定で、**report only** として使います
- `qc_spur_on_apply_flagrow=True` のときだけ ON の auto bad を `FLAGROW` に反映します

これにより、ON を見たいが自動で捨てたくはない、という運用ができます。

---

## 7. 1 ch から 3 ch を探す意味

本実装で狙っている狭帯域 spur は、基本的に 1 ch から `qc_spur_max_consecutive_channels` ch です。

もし 1 ch あたりの速度幅を $\Delta v_{\mathrm{ch}}$、本物の最も狭い線幅を $\Delta v_{\mathrm{line,min}}$ とすると、

$$
\Delta v_{\mathrm{line,min}} \gg 3 \, \Delta v_{\mathrm{ch}}
$$

であれば、1 ch から 3 ch の鋭い構造は本物の線では出にくいので、ON にもかなり安全に使えます。

一方で

$$
\Delta v_{\mathrm{line,min}} \lesssim 3 \, \Delta v_{\mathrm{ch}}
$$

なら、本物の狭い線も 1 ch から 3 ch に見え得ます。このため ON では、狭帯域幅だけでなく**孤立性**も併用します。

---

## 8. 高速化

### 8.1 HOT / OFF

HOT / OFF は本数が少ないため、原則として高速化は優先しません。ここは安定性を優先します。

### 8.2 ON

高速化の主対象は ON です。

旧方式では、`qc_spur_consider_on=True` のとき ON 全 row に rolling median 本判定をかけていたため、計算量は概ね

$$
N_{\mathrm{on}} \times N_{\mathrm{ch}} \times W
$$

に比例して増えます。ここで $W$ は `qc_spur_window_channels` です。

新方式では pre-screen 通過率を $f$ とすると、重い本判定は概ね

$$
f \, N_{\mathrm{on}} \times N_{\mathrm{ch}} \times W
$$

になります。通常は $f \ll 1$ を期待するので、かなりの高速化が得られます。

### 8.3 ベンチマーク結果

現時点の簡易シミュレーションでは次でした。

#### HOT / OFF 用 trial 補正 threshold

| 有効チャネル数 $N$ | 新方式の有効 threshold |
|---:|---:|
| 256 | 6.901947 |
| 512 | 6.999994 |
| 768 | 7.056678 |
| 1024 | 7.096611 |
| 2048 | 7.191882 |
| 4096 | 7.285997 |

#### 境界的な 1 ch spike の検出率

| spike 振幅 | $N$ | 旧方式検出率 | 新方式検出率 |
|---:|---:|---:|---:|
| 7.1 | 512 | 0.597 | 0.597 |
| 7.1 | 4096 | 0.633 | 0.527 |
| 7.3 | 512 | 0.663 | 0.663 |
| 7.3 | 4096 | 0.720 | 0.647 |
| 7.5 | 512 | 0.743 | 0.743 |
| 7.5 | 4096 | 0.713 | 0.597 |
| 8.0 | 512 | 0.880 | 0.880 |
| 8.0 | 4096 | 0.860 | 0.793 |

#### ON 高速化

| $N$ | ON row 数 | 注入 spike 数 | pre-screen 通過 row 数 | 旧方式時間 [s] | 新方式時間 [s] | speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 512 | 1200 | 12 | 12 | 1.172 | 0.274 | 4.27 |
| 4096 | 800 | 8 | 8 | 2.298 | 0.374 | 6.14 |

### 8.4 解釈

- 1000 ch 以下では旧方式でも大きな破綻は出にくい
- 512 ch 付近では新旧の threshold はほぼ同じ
- 数千 ch では HOT / OFF の過剰 flag を抑えやすい
- ON は pre-screen により数倍高速化できる

---

## 9. 出力列と記録

### 9.1 fixed bad channel 関連

fixed bad channel を指定したとき、ON table には次が追加されます。

- `BADCHAN_POLICY`
- `BADCHAN_N_GLOBAL`
- `BADCHAN_N_GROUP`
- `BADCHAN_N_TOTAL`
- `BADCHAN_INTERP_APPLIED`
- `BADCHAN_INTERP_NFILLED`
- `BADCHAN_HAS_EDGE_UNFILLED`
- `BADCHAN_SHAPE_INTERP_APPLIED`

### 9.2 ON QC 関連

`qc_spur_sigma is not None` かつ `qc_spur_consider_on=True` のとき、次が追加されます。

- `REF_QC_ON_PRESCREEN_QD`
- `REF_QC_ON_PRESCREEN_CANDIDATE`
- `REF_QC_ON_SPUR_SCORE`
- `REF_QC_ON_SPUR_ABSPEAK`
- `REF_QC_ON_SPUR_RUNLEN`
- `REF_QC_ON_SPUR_THRESHOLD`
- `REF_QC_ON_NEIGHBOR_SCORE`
- `REF_QC_ON_ISO_SCORE`
- `REF_QC_ON_SPUR_AUTO_BAD`

### 9.3 集計列

QC を使ったとき、table には次が入ります。

- `QC_REF_HOT_AUTO_BAD`
- `QC_REF_OFF_AUTO_BAD`
- `QC_REF_HOT_SHAPE_AUTO_BAD`
- `QC_REF_OFF_SHAPE_AUTO_BAD`
- `QC_REF_HOT_SPUR_AUTO_BAD`
- `QC_REF_OFF_SPUR_AUTO_BAD`
- `QC_ON_PRESCREEN_CANDIDATE`
- `QC_ON_SPUR_AUTO_BAD`
- `QC_SIGMA`
- `QC_SPUR_SIGMA`
- `QC_SPUR_WINDOW_CH`
- `QC_SPUR_MAXRUN`
- `QC_SPUR_TRIAL_REF_CH`
- `QC_SPUR_CONSIDER_ON`
- `QC_SPUR_ON_APPLY_FLAGROW`
- `QC_SPUR_ON_PRESCREEN_QD`
- `QC_SPUR_ON_ISO_MIN_RATIO`
- `QC_APPLY_FLAGROW`

### 9.4 `attrs` に入る summary

`table.attrs` には次が記録されます。

- `qc_summary`
- `bad_channel_summary`

これにより、入力 `FLAGROW`、HOT / OFF の shape / spur の件数、ON pre-screen 候補件数、ON auto bad 件数、補間適用行数などを追跡できます。

---

## 10. 推奨設定

### 10.1 まずの推奨

まずは report 的に確認するのが安全です。

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    qc_spur_on_apply_flagrow=False,
    qc_apply_flagrow=False,
)
```

この設定では

- HOT / OFF の shape QC を実行
- HOT / OFF の narrow spur QC を実行
- ON は見ない
- `FLAGROW` にはまだ反映しない

です。

### 10.2 HOT / OFF を実際に落とす標準設定

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    qc_apply_flagrow=True,
)
```

これは最初の標準設定です。

### 10.3 1 ch spur を主に見たい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.5,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=1,
    qc_spur_consider_on=False,
    qc_apply_flagrow=True,
)
```

### 10.4 2 ch から 3 ch spur も見たい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    qc_apply_flagrow=True,
)
```

### 10.5 weak spur を取りこぼしたくない場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=6.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    qc_apply_flagrow=False,
)
```

最初から `qc_apply_flagrow=True` にせず、まず件数を見る方が安全です。

### 10.6 fixed bad channel がある場合

1 ch だけでなく、少し広めの range 指定が有効なことがあります。

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    bad_channel_ranges=[(1534, 1538)],
    bad_channel_policy="interp_linear",
    bad_channel_interp_max_consecutive_channels=5,
    qc_spur_consider_on=False,
    qc_apply_flagrow=True,
)
```

### 10.7 group ごとに異なる fixed bad channel がある場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    bad_channel_key_columns=("FDNUM", "IFNUM", "PLNUM"),
    bad_channel_map={
        (2, 0, 0): [(1534, 1538)],
        (3, 0, 1): [2048],
    },
    bad_channel_policy="interp_linear",
    qc_spur_consider_on=False,
    qc_apply_flagrow=True,
)
```

### 10.8 ON を report only で見る場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=True,
    qc_spur_on_apply_flagrow=False,
    qc_apply_flagrow=True,
)
```

この場合、ON の

- pre-screen
- spur score
- neighbor score
- isolation score

は計算されますが、ON の `FLAGROW` は立ちません。

### 10.9 ON も自動で落としたい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=True,
    qc_spur_on_apply_flagrow=True,
    qc_apply_flagrow=True,
)
```

ただしこれは line-rich source では慎重に使うべきです。

### 10.10 銀河中心や line-rich source の場合

銀河中心のように速度的に広がった本物の emission が強いときは、ON の自動 flag は危険です。

その場合はまず

```python
qc_spur_consider_on=False
```

を基本にするか、ON を見るなら

```python
qc_spur_consider_on=True
qc_spur_on_apply_flagrow=False
```

の report only を推奨します。

### 10.11 速度優先と見逃し優先

#### 速度優先

- `qc_spur_sigma` をやや高めにする
- `qc_spur_max_consecutive_channels` を 1 か 2 にする
- ON は `qc_spur_consider_on=True` かつ `qc_spur_on_apply_flagrow=False`

#### 見逃し優先

- `qc_spur_sigma` を少し下げる
- `qc_spur_max_consecutive_channels` を 3 にする
- まず `qc_apply_flagrow=False` で件数確認

---

## 11. 実務上の推奨手順

1. fixed bad channel が分かっているなら最初に指定する
2. まず `qc_apply_flagrow=False` で件数を見る
3. HOT / OFF の shape / spur 件数を確認する
4. ON を見る必要があるときだけ `qc_spur_consider_on=True` にする
5. ON は最初は `qc_spur_on_apply_flagrow=False` にする
6. 問題なければ `qc_apply_flagrow=True` にする

---

## 12. まとめ

現在の `calibrate_updated_2026-04-10_statewise.py` の QC は、次のように整理されています。

### 共通前処理

- fixed bad channel を channel-level で処理する
- QC 用配列、shape 用配列、較正用配列を分ける

### HOT / OFF

- `qc_sigma` による shape QC
- `qc_spur_sigma` による narrow spur QC
- narrow spur は science window 内で trial 数補正付き threshold を使う

### ON

- `$Q_D$` による pre-screen
- narrow spur 本判定
- 近傍 ON row との孤立性判定
- 既定は report only

この構成により、

- science window をそのまま使いながら
- fixed bad channel と row 異常を分離し
- HOT / OFF と ON にそれぞれ向いた判定を使い
- ON では高速化も得る

という設計になっています。
