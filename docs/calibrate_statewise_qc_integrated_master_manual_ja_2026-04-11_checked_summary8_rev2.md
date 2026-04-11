# `calibrate.py` state-wise QC 統合詳細説明書

## 0. この文書の位置づけ

この文書は、state-wise QC、fixed bad channel、最終 Ta 出力補間を 1 本に統合し、**summary8 実装**へ合わせて更新した統合詳細説明書です。

本書が対象とする実装上の基準コードは、**最新版**

- `calibrate_updated_2026-04-11_statewise_peaksummary8.py`

です。

この版では、以前の文書に残っていた

- `peaksummary5` 前提の記述
- `interp_linear_ta` を既定推奨としていた記述
- `poly2_weighted_ta` 未反映の記述
- summary 列と `history["bad_channel_summary"]` の旧 key 一覧

を見直し、**summary8 の実装と整合する内容**へ更新しています。

また、本書では次を一体として扱います。

1. state-wise QC の基本設計
2. fixed bad channel の扱い
3. HOT / OFF の shape QC と narrow spur QC
4. ON の pre-screen, narrow spur 本判定, 孤立性判定
5. `qc_summary` と `bad_channel_summary` の確認方法
6. `output_bad_channel_fill_policy="none"`, `"interp_linear_ta"`, `"poly2_weighted_ta"` の使い分け
7. 実務向けの推奨設定と多数の使用例
8. 10 ch gap, OTF-like 平均, ripple を踏まえた実務上の注意

---

## 1. 基本方針

本実装の基本方針は、**共通前処理 + state 別判定器** です。

すなわち、

- 全 state 共通で fixed bad channel を処理する
- HOT / OFF には参照信号向けの row QC を使う
- ON には天体線と spur を分けやすい判定器を使う
- 高速化は主に ON に入れる
- 最終的に見せたいスペクトルが Ta なら、補間も原則 Ta 側で行う

という構成です。

### 1.1 記号

本書では次の記号を使います。

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
\sigma_{R,x}(r) = 1.4826 \, \mathrm{MAD}\left(R_x(k, r),\ k \in V\right)
$$

- 規格化 residual:

$$
z_x(k, r) = \frac{|R_x(k, r)|}{\sigma_{R,x}(r)}
$$

---

## 2. 公開 API と実装済み引数

最新版 `make_tastar_dumps(...)`, `tastar_from_rawspec(...)`, `run_tastar_calibration(...)` では、QC / bad channel / 出力補間まわりとして次を受け付けます。

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

### 2.4 最終 Ta 補間

- `output_bad_channel_fill_policy`
- `output_bad_channel_fill_max_consecutive_channels`
- `output_bad_channel_fill_min_valid_neighbors`
- `output_bad_channel_fill_half_window_channels`
- `output_bad_channel_fill_curvature_sigma_threshold`
- `output_bad_channel_fill_edge_consistency_sigma_max`
- `output_bad_channel_fill_fallback_policy`

### 2.5 `output_bad_channel_fill_policy` の選択肢

summary8 では、最終 Ta 補間 policy として次を受け付けます。

- `"none"`
- `"interp_linear_ta"`
- `"poly2_weighted_ta"`

### 2.6 重要な整理

役割は次のように分かれます。

- `bad_channel_policy`
  - raw の ON / OFF / HOT をどう扱うか
- `output_bad_channel_fill_policy`
  - **最終 Ta 出力**の bad channel をどう扱うか

本書の現在の推奨は、**原則として**

- `bad_channel_policy="nan"`
- `output_bad_channel_fill_policy="none"`

です。

理由は単純で、1 ch の isolated spur と違って、**5 ch から 10 ch 程度の連続 bad gap は補間で安全に回復できるとは限らない**からです。特に OTF-like に位置の少し違う row を重み付き平均すると、単一 row では小さく見えた補間誤差が平均後に偏って残り、深い dip として見えることがあります。

したがって、補間は既定ではなく、**どうしても quicklook / viewer 上の連続性が必要なときだけ明示的に使う**方針を本書の標準とします。

---

## 3. 処理全体の流れ

概略フローは次です。

```text
raw ON/OFF/HOT
  -> row selection
  -> 既存 FLAGROW の candidate 除外
  -> science window 切り出し
  -> fixed bad channel 処理
       - spur QC 用配列
       - shape QC 用配列
       - 較正計算用配列
  -> HOT/OFF shape QC
  -> HOT/OFF narrow spur QC
  -> ON pre-screen
  -> ON narrow spur 本判定
  -> ON 孤立性判定
  -> 必要なら FLAGROW 反映
  -> HOT/OFF の平均・時間内挿
  -> Ta 計算
  -> 必要なら最終 Ta 補間
  -> summary / history / table 列へ記録
```

---

## 4. fixed bad channel

### 4.1 基本方針

fixed bad channel は row 異常ではなく、**channel-level の恒常的不良**です。

したがって、fixed bad channel は `FLAGROW` では処理せず、**channel mask として扱います**。

### 4.2 指定方法

#### `bad_channels`

全 group 共通の raw channel index を指定します。

```python
bad_channels=[1536, 3072]
```

#### `bad_channel_ranges`

全 group 共通の range を指定します。range は **両端を含む closed interval** です。

```python
bad_channel_ranges=[(1534, 1538), (3070, 3074)]
```

#### `bad_channel_map`

group ごとに異なる fixed bad channel を指定します。key の順は `bad_channel_key_columns` と一致させます。

```python
bad_channel_key_columns=("FDNUM", "IFNUM", "PLNUM")

bad_channel_map = {
    (2, 0, 0): [1536],
    (3, 0, 1): [(1534, 1538)],
}
```

### 4.3 インデックスの意味

ここで指定する channel index は、**raw の 0-based channel index** です。

`ch_range` や `vlsrk_range_kms` で切り出した後の local index ではありません。

たとえば raw 1536 ch を bad にしたいなら、`1536` を指定します。`ch_range=(1000, 2000)` を使っていても `536` ではありません。

### 4.4 global と group の合成

group $g$ に対する最終 bad channel 集合は、global と group 別指定の和集合です。

$$
K_{\mathrm{bad}}(g) = K_{\mathrm{global}} \cup K_{\mathrm{group}}(g)
$$

### 4.5 bad scan と QC の順序

既存 `FLAGROW` によって bad とされている row は、**QC の判定対象から先に外れます**。

ただし実装は「先に配列を短く切る」のではなく、まず

$$
\mathrm{candidate\_mask} = \lnot \mathrm{input\_bad}
$$

を作って、candidate row だけを QC で見ます。

その後、

- 既存 `FLAGROW` bad
- shape QC による auto bad
- spur QC による auto bad

をまとめて `FLAGROW` に反映し、最後に残った row だけを平均・内挿・較正に進めます。

したがって、**bad scan を含んだまま他の row を判定する**という順序ではありません。

---

## 5. fixed bad channel を処理した 4 種類の配列

実装上、bad channel を処理した配列は 4 つの役割に分かれます。

### 5.1 spur QC 用配列

狭帯域 spur QC では、bad channel は **常に除外**します。該当チャネルは NaN にし、finite なチャネルだけを見ます。

### 5.2 shape QC 用配列

shape QC では、fixed bad channel 自体で row が落ちないように、**局所線形補間済みの専用配列**を使います。

重要なのは、これは science 値の復元ではなく、**QC の安定化のための内部補助配列**だという点です。

### 5.3 較正計算用配列

較正計算に実際に使う raw ON / OFF / HOT の配列です。ここでは `bad_channel_policy` が効きます。

- `"nan"`: bad channel は NaN のまま
- `"interp_linear"`: 短い bad 区間だけ raw で線形補間

### 5.4 最終 Ta 出力配列

最終 Ta を作った後にだけ適用される補間です。ここでは `output_bad_channel_fill_policy` が効きます。

- `"none"`: 何もしない
- `"interp_linear_ta"`: 最終 Ta の bad channel を線形補間
- `"poly2_weighted_ta"`: 最終 Ta の bad channel に対し、gap ごとに局所 2 次 weighted 補間を試し、不適当なら `linear` または `none` に fallback

現在の summary8 実装では、`poly2_weighted_ta` は **無条件 2 次補間ではありません**。gap ごとに

1. 2 次を試す価値があるか
2. 2 次が gap の外側端点と矛盾していないか
3. だめなら `linear` または `none` に戻すか

を判定します。

---

## 6. raw 側の bad channel 補間

### 6.1 `bad_channel_policy="nan"`

これは **科学解析の標準**です。

raw 側では bad channel を使わず、NaN のまま下流へ伝播させます。

長所:

- 余計な値を作らない
- 非線形な較正式を勝手にゆがめない
- どの段階で欠損したかが明確

短所:

- viewer や quicklook では穴が見える
- 一部の下流ツールで NaN を嫌うことがある

### 6.2 `bad_channel_policy="interp_linear"`

raw ON / OFF / HOT を state ごとに線形補間します。

bad な連続区間を $[k_1, k_2]$ とし、その外側に有効チャネル $k_L < k_1$ と $k_R > k_2$ があるとき、各 state について

$$
S'_x(k, r)
=
S_x(k_L, r)
+
\frac{k-k_L}{k_R-k_L}
\left(S_x(k_R, r)-S_x(k_L, r)\right)
$$

で埋めます。

### 6.3 raw 補間の条件

`interp_linear` では次を満たすときだけ補間します。

- edge ではない
- 両側に有効チャネルがある
- 連続 bad 区間長が `bad_channel_interp_max_consecutive_channels` 以下
- `bad_channel_interp_min_valid_neighbors` を満たす

満たさない場合は NaN のまま残します。

### 6.4 なぜ raw 補間を既定推奨にしないか

Ta は概ね

$$
T(k) = \frac{N(k)}{D(k)}
$$

の形です。ここで

$$
N(k) = \mathrm{ON}(k) - \mathrm{OFF}(k)
$$

$$
D(k) = \mathrm{HOT}(k) - \mathrm{OFF}(k)
$$

です。

raw を別々に補間してから比を取ると

$$
\widehat{T}(k)
=
\frac{\widehat{\mathrm{ON}}(k)-\widehat{\mathrm{OFF}}(k)}
     {\widehat{\mathrm{HOT}}(k)-\widehat{\mathrm{OFF}}(k)}
$$

になりますが、一般には

$$
\frac{\widehat N(k)}{\widehat D(k)} \neq \widehat{\left(\frac{N(k)}{D(k)}\right)}
$$

です。

このため、raw で別々に補間すると、bad channel 部分に

- 不自然な dip
- 不自然な bump
- 吸収っぽい構造

が出ることがあります。

したがって、本書では **raw 側は NaN を基本**とし、必要なら **最終 Ta で補間**する方針を推奨します。

---

## 7. 最終 Ta 補間

### 7.1 基本方針

最終的に見たい物理量が Ta なら、補間も Ta に対して行う方が raw 側よりは自然です。

ただし、**Ta 側で補間すること自体が安全という意味ではありません**。

特に 5 ch から 10 ch 程度の連続 gap では、line center 付近や shoulder のように局所曲率が強い場所で、単純線形補間は

- 人工的な dip
- 人工的な bump
- OTF 平均後の負バイアス

を作り得ます。

したがって、現在の基本推奨は

- raw 側: `bad_channel_policy="nan"`
- 最終 Ta 側: `output_bad_channel_fill_policy="none"`

です。

### 7.2 `output_bad_channel_fill_policy`

summary8 では次を受け付けます。

- `"none"`
- `"interp_linear_ta"`
- `"poly2_weighted_ta"`

### 7.3 `"none"`

何もしません。bad channel は最終 Ta にも NaN のまま残ります。

これは見た目には穴が残りますが、**値を勝手に作らない**という意味で最も安全です。

本書では、science 用の基本既定をこれにします。

### 7.4 `"interp_linear_ta"`

bad な連続区間を $[k_1, k_2]$ とし、その外側に有効な Ta の値があるとき、

$$
T'(k, r)
=
T(k_L, r)
+
\frac{k-k_L}{k_R-k_L}
\left(T(k_R, r)-T(k_L, r)\right)
$$

で埋めます。

これは実装が単純で quicklook には便利ですが、局所曲率を表現できません。

### 7.5 `"poly2_weighted_ta"`

summary8 で追加された最終 Ta 側補間です。bad run ごとに gap 中心

$$
k_c = \frac{k_1 + k_2}{2}
$$

を定め、近傍有効点

$$
G = \{ j \notin B \mid |j-k_c| \le h \}
$$

から局所 1 次式と局所 2 次式を作ります。ここで $h$ は半窓幅です。

距離重みは、実装では tri-cube 形で

$$
w(j) = \left(1 - \left|\frac{j-k_c}{h}\right|^3\right)^3
$$

を、$|j-k_c| < h$ の点に使います。

そのうえで、gap 内での 1 次予測と 2 次予測の差

$$
\Delta_{2-1} = \max_{k \in B} |T_2(k) - T_1(k)|
$$

を局所スケール $\sigma_{\mathrm{loc}}$ で割った

$$
S_{\mathrm{curv}} = \frac{\Delta_{2-1}}{\sigma_{\mathrm{loc}}}
$$

を曲率指標とします。

さらに、gap 外側端点で 2 次式が観測値とどれだけ一致するかを

$$
S_{\mathrm{edge}} =
\frac{
\max
\left(
|T_2(k_L)-T(k_L)|,
|T_2(k_R)-T(k_R)|
\right)
}{\sigma_{\mathrm{loc}}}
$$

で見ます。

実装では概念的に

- $S_{\mathrm{curv}}$ が十分大きい
- $S_{\mathrm{edge}}$ が十分小さい

ときだけ 2 次式を採用し、それ以外は fallback へ戻します。

### 7.6 `poly2_weighted_ta` の実装上の重要点

`poly2_weighted_ta` は **常に 2 次で埋めるものではありません**。実装上は次の順で判定します。

1. gap 長が `output_bad_channel_fill_max_consecutive_channels` を超えるなら埋めない
2. 近傍有効点が不足するなら fallback
3. 近傍左右の点が足りないなら fallback
4. 1 次式と 2 次式を両方作る
5. `curvature_score` と `edge_score` を見て 2 次採用を判断
6. だめなら `linear` または `none` へ戻る

したがって、line-free や edge gap では 2 次が一度も使われないことがあります。

### 7.7 `interp_linear_ta` と `poly2_weighted_ta` の使い分け

実務上の整理は次です。

- **science 既定**: `none`
- **quicklook で穴をつなぎたいだけ**: `interp_linear_ta`
- **line center や shoulder で 10 ch 前後の gap による dip を少しでも減らしたい**: `poly2_weighted_ta`

ただし、`poly2_weighted_ta` でも真値再現は保証されません。特に

- gap が長すぎる
- 本物の非常に narrow な線が gap の中に完全に入る
- line-free ではなく ripple 優勢
- 片側近傍しか無い edge gap

では、安全側は依然として `none` です。

### 7.8 viewer での見え方

`output_bad_channel_fill_policy="interp_linear_ta"` あるいは `"poly2_weighted_ta"` を使うと、viewer で見えている Ta 側の該当 channel が補間値になります。

一方、raw 側 `interp_linear` の場合は、viewer に見えるのは補間 raw そのものではなく **その後に較正した Ta** なので、見た目が直感と一致しないことがあります。

### 7.9 本書の推奨結論

- 1 ch isolated spur なら、補間を使っても大崩れしない場合がある
- 5 ch から 10 ch 程度の連続 bad gap は、補間で自然に回復する保証がない
- OTF-like に位置の少し違う row を平均すると、補間誤差が系統的 dip として残り得る

したがって、**基本は NaN, NaN**、すなわち

```python
bad_channel_policy="nan"
output_bad_channel_fill_policy="none"
```

を本書の標準推奨とします。

補間は、どうしても必要な場合にだけ明示的に有効化してください。

---

## 8. HOT / OFF の QC

### 8.1 基本方針

HOT / OFF には本物の天体線が入らないので、**row 自動 flag の主戦場**にしてよいです。

HOT / OFF の最終判定は

$$
\mathrm{auto\_bad} = \mathrm{shape\_bad} \lor \mathrm{narrow\_spur\_bad}
$$

です。

### 8.2 broad な異常: shape QC

これは既存の `qc_sigma` 系です。`grouped_affine_shape_qc(...)` によって、同じ group の HOT / OFF の中で broad な shape 外れ値を見ます。

fixed bad channel の影響を減らすため、shape QC には **shape 用補間配列**を入力に使います。

### 8.3 狭帯域異常: narrow spur QC

各 row について rolling median residual を作り、1 ch から `qc_spur_max_consecutive_channels` ch までの narrow run を探します。

各 narrow run $C$ の local score を

$$
T(C) = \max_{k \in C} z_x(k, r)
$$

とし、その row の代表 score を

$$
T_{\max}(r) = \max_C T(C)
$$

とします。

### 8.4 trial 数補正

science window をそのまま使うと、有効チャネル数が多いほど「どこかで 1 回引っかかる」確率が上がります。

そこで HOT / OFF の narrow spur QC では、**trial 数補正付き threshold** を使います。

試行回数は単純な $N$ ではなく、run 長を 1 から $L_{\max}$ まで試した総数として

$$
M(N, L_{\max}) = \sum_{\ell=1}^{L_{\max}} (N - \ell + 1)
$$

と定義します。

実装では、`qc_spur_sigma` を基準有効チャネル数

$$
N_{\mathrm{ref}} = 512
$$

に対する local threshold とみなし、実際の $N$ に対する有効 threshold を内部で計算します。

### 8.5 最新修正: `qc_spur_sigma` の大値でも単調に効く

以前の途中版では、試行補正の数値計算が極端な tail で underflow し、`qc_spur_sigma` を非常に大きくしても内部 threshold が約 7.94 に張り付く不具合がありました。

最新版コードではこの点を修正し、**大きい `qc_spur_sigma` でも単調に threshold が上がる**ようになっています。

したがって、最新版では

- `qc_spur_sigma = 9`
- `qc_spur_sigma = 12`
- `qc_spur_sigma = 20`
- `qc_spur_sigma = 1e20`

がすべて区別されます。

### 8.6 `qc_spur_sigma` を上げても件数が変わらないとき

それはしばしば

- fixed bad channel が十分に除外されていない
- 1 ch 指定が狭すぎる
- raw index がずれている
- spur が本当に非常に強い

ことを意味します。

この場合は threshold を上げ続けるよりも、まず

- `bad_channel_ranges`
- `bad_channel_map`
- `qc_spur_max_consecutive_channels=1`
- `qc_spur_window_channels=31`

の調整が本筋です。

---

## 9. ON の QC

### 9.1 基本方針

ON では本物の天体線が入るため、HOT / OFF と同じ row 自動 flag を主判定にしてはいけません。

現在の実装では ON について

1. pre-screen
2. narrow spur 本判定
3. 近傍 ON row との孤立性判定

の 3 段に分けています。

### 9.2 ON pre-screen

各 ON row のスペクトルを $S(i)$ とし、隣接チャネル差分を

$$
D(i) = S(i+1) - S(i)
$$

とします。

差分系列の robust 尺度を

$$
\sigma_D = 1.4826 \, \mathrm{MAD}(D)
$$

とし、最大規格化振幅を

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

### 9.3 直感的意味

狭帯域 spur は 1 ch から数 ch の間で急に立ち上がって急に戻ることが多いので、隣接差分に大きく出やすいです。

一方で

- 全体 gain 変化
- 緩やかな bandpass 傾き
- 広い線構造
- ゆっくりした ripple

は、隣接差分には強くは出にくいです。

### 9.4 重要な解釈

$Q_D$ は **bad 判定そのものではありません**。

正確には、

- $Q_D$ が大きい row を本判定へ回す
- $Q_D$ が小さい row には、少なくとも「隣接チャネルで急に跳ぶ型の 1 ch から数 ch の狭帯域 spur」は出にくい

という意味です。

### 9.5 ON narrow spur 本判定

pre-screen を通過した ON row に対してだけ、HOT / OFF と同様の rolling median residual ベース本判定を行います。

### 9.6 ON の孤立性判定

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

### 9.7 ON の `FLAGROW` 反映

- `qc_spur_consider_on=False` なら ON QC をしない
- `qc_spur_consider_on=True` なら ON QC を実行する
- `qc_spur_on_apply_flagrow=False` が既定で、**report only** として使う
- `qc_spur_on_apply_flagrow=True` のときだけ ON の auto bad を `FLAGROW` に反映する

### 9.8 ON に対する注意

ON には本物のスペクトル線が入ります。特に

- 強い狭線
- maser 的に非常に細い線
- 非常に高い S/N の狭い線

では $Q_D$ や本判定 score が大きくなり得ます。

したがって、line-rich source では

- 最初は `qc_spur_consider_on=False`
- あるいは `qc_spur_consider_on=True`, `qc_spur_on_apply_flagrow=False`

から始めるのが安全です。

---

## 10. 1 ch から 3 ch を探す意味

狙う異常幅を $L = 1$ から $3$ ch とするのは、典型的な digital spur を狙うには自然です。

ただし ON では、これが常に安全という意味ではありません。

安全かどうかは、1 ch あたりの速度幅を $\Delta v_{\mathrm{ch}}$、本物の最小線幅を $\Delta v_{\mathrm{line,min}}$ とすると、概ね

$$
\Delta v_{\mathrm{line,min}} \gg 3 \, \Delta v_{\mathrm{ch}}
$$

が言えるかに依存します。

### 10.1 ON に使ってよい場面

- 本物の線は 1-3 ch より十分広い
- 熱的な分子雲線で narrow digital spike だけ落としたい
- 自動 `FLAGROW` ではなく report only で使う

### 10.2 ON に慎重であるべき場面

- maser
- 非常に高分解能
- 鋭い吸収線
- 自己吸収の急峻な肩
- 本物の線が 1-3 ch で見え得る観測

---

## 11. 高速化

### 11.1 HOT / OFF

HOT / OFF は本数が少ないので、基本は高速化不要です。ここでは安定性を優先します。

### 11.2 ON

ON は row 数が多いので、pre-screen による間引きが効きます。

重い本判定部分の計算量は概ね

$$
N_{\mathrm{on}} \times N_{\mathrm{ch}} \times W
$$

に比例します。ここで $W$ は rolling median 窓幅です。

pre-screen を通る割合を $f$ とすると、重い本判定は概ね

$$
f \, N_{\mathrm{on}} \times N_{\mathrm{ch}} \times W
$$

になります。

### 11.3 実装上の重要点

`qc_spur_on_apply_flagrow=False` は **ON を落とさない**だけで、**ON QC の計算自体は実行します**。したがって、速度には直接効きません。

速度を決めるのは主に

- `qc_spur_consider_on`
- pre-screen を何本通過するか

です。

### 11.4 ログの見方

ログに

- `on_spur_prescreen_candidate=...`
- `on_spur_auto_bad=...`

が出ます。

`on_spur_prescreen_candidate` が大きいなら、pre-screen でほとんど間引けていません。line-rich source ではこれが起こり得ます。

---

## 12. 出力列と記録

### 12.1 fixed bad channel 関連

fixed bad channel を指定したとき、ON table には少なくとも次が入ります。

- `BADCHAN_POLICY`
- `BADCHAN_N_GLOBAL`
- `BADCHAN_N_GROUP`
- `BADCHAN_N_TOTAL`
- `BADCHAN_INTERP_APPLIED`
- `BADCHAN_INTERP_NFILLED`
- `BADCHAN_HAS_EDGE_UNFILLED`
- `BADCHAN_SHAPE_INTERP_APPLIED`
- `BADCHAN_OUTFILL_POLICY`
- `BADCHAN_OUTFILL_APPLIED`
- `BADCHAN_OUTFILL_NFILLED`
- `BADCHAN_OUTFILL_HAS_EDGE_UNFILLED`
- `BADCHAN_OUTFILL_POLY2_APPLIED`
- `BADCHAN_OUTFILL_LINEAR_APPLIED`
- `BADCHAN_OUTFILL_POLY2_RUNS`
- `BADCHAN_OUTFILL_LINEAR_RUNS`
- `BADCHAN_OUTFILL_NONE_RUNS`
- `BADCHAN_OUTFILL_CURV_SCORE_MAX`
- `BADCHAN_OUTFILL_EDGE_SCORE_MAX`
- `BADCHAN_OUTFILL_LOCAL_SIGMA_MED`

### 12.2 ON QC 関連

`qc_spur_sigma is not None` かつ `qc_spur_consider_on=True` のとき、次が入ります。

- `REF_QC_ON_PRESCREEN_QD`
- `REF_QC_ON_PRESCREEN_CANDIDATE`
- `REF_QC_ON_SPUR_SCORE`
- `REF_QC_ON_SPUR_ABSPEAK`
- `REF_QC_ON_SPUR_RUNLEN`
- `REF_QC_ON_SPUR_THRESHOLD`
- `REF_QC_ON_NEIGHBOR_SCORE`
- `REF_QC_ON_ISO_SCORE`
- `REF_QC_ON_SPUR_AUTO_BAD`

### 12.3 集計列

QC を使ったとき、table には次が入ります。

- `QC_REF_HOT_AUTO_BAD`
- `QC_REF_OFF_AUTO_BAD`
- `QC_REF_HOT_SHAPE_AUTO_BAD`
- `QC_REF_OFF_SHAPE_AUTO_BAD`
- `QC_REF_HOT_SPUR_AUTO_BAD`
- `QC_REF_OFF_SPUR_AUTO_BAD`
- `QC_REF_HOT_SPUR_TOP_CH_RAW`
- `QC_REF_HOT_SPUR_TOP_COUNT`
- `QC_REF_OFF_SPUR_TOP_CH_RAW`
- `QC_REF_OFF_SPUR_TOP_COUNT`
- `QC_REF_HOT_SPUR_PEAK_COUNTS_RAW_JSON`
- `QC_REF_OFF_SPUR_PEAK_COUNTS_RAW_JSON`
- `QC_REF_HOT_SPUR_PEAK_TOP_RAW_JSON`
- `QC_REF_OFF_SPUR_PEAK_TOP_RAW_JSON`
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

### 12.4 `attrs` に入る summary

`table.attrs` には次が入ります。

- `qc_summary`
- `bad_channel_summary`

### 12.5 `history` に残る summary

最終的な `Scantable` の `history` にも

- `history["qc_summary"]`
- `history["bad_channel_summary"]`

が残ります。

**最終結果の確認は `history` から見るのが最も確実**です。

---

## 13. summary の確認方法

### 13.1 まず `history["qc_summary"]`

```python
qc_summary = sc_cal.history.get("qc_summary", {})

qc_summary.get("hot_input_flagrow_bad", 0)
qc_summary.get("off_input_flagrow_bad", 0)
qc_summary.get("on_input_flagrow_bad", 0)
qc_summary.get("hot_auto_bad_would_flag", 0)
qc_summary.get("off_auto_bad_would_flag", 0)
qc_summary.get("hot_shape_auto_bad_would_flag", 0)
qc_summary.get("off_shape_auto_bad_would_flag", 0)
qc_summary.get("hot_spur_auto_bad_would_flag", 0)
qc_summary.get("off_spur_auto_bad_would_flag", 0)
qc_summary.get("on_spur_prescreen_candidate", 0)
qc_summary.get("on_spur_auto_bad_would_flag", 0)
qc_summary.get("hot_used_rows", 0)
qc_summary.get("off_used_rows", 0)
```

### 13.2 HOT / OFF の spur peak channel

```python
qc_summary.get("hot_spur_peak_channel_counts_raw", {})
qc_summary.get("off_spur_peak_channel_counts_raw", {})
qc_summary.get("hot_spur_peak_channel_top_raw", [])
qc_summary.get("off_spur_peak_channel_top_raw", [])
```

### 13.3 table の top channel 列

```python
sc_cal.table[[
    "QC_REF_HOT_SPUR_TOP_CH_RAW",
    "QC_REF_HOT_SPUR_TOP_COUNT",
    "QC_REF_OFF_SPUR_TOP_CH_RAW",
    "QC_REF_OFF_SPUR_TOP_COUNT",
]].iloc[0]
```

### 13.4 table の JSON 列

```python
sc_cal.table[[
    "QC_REF_HOT_SPUR_PEAK_COUNTS_RAW_JSON",
    "QC_REF_OFF_SPUR_PEAK_COUNTS_RAW_JSON",
    "QC_REF_HOT_SPUR_PEAK_TOP_RAW_JSON",
    "QC_REF_OFF_SPUR_PEAK_TOP_RAW_JSON",
]].iloc[0]
```

### 13.5 `history["bad_channel_summary"]`

```python
bad_summary = sc_cal.history.get("bad_channel_summary", {})

bad_summary.get("policy")
bad_summary.get("interp_max_consecutive_channels")
bad_summary.get("output_fill_policy")
bad_summary.get("output_fill_max_consecutive_channels")
bad_summary.get("output_fill_min_valid_neighbors")
bad_summary.get("output_fill_half_window_channels")
bad_summary.get("output_fill_curvature_sigma_threshold")
bad_summary.get("output_fill_edge_consistency_sigma_max")
bad_summary.get("output_fill_fallback_policy")

bad_summary.get("global_bad_channels_raw_count")
bad_summary.get("shape_interp_applied_rows")

bad_summary.get("on_interp_rows")
bad_summary.get("off_interp_rows")
bad_summary.get("hot_interp_rows")

bad_summary.get("on_interp_nfilled_total")
bad_summary.get("off_interp_nfilled_total")
bad_summary.get("hot_interp_nfilled_total")

bad_summary.get("on_output_fill_rows")
bad_summary.get("on_output_fill_nfilled_total")
bad_summary.get("on_output_fill_has_edge_unfilled_any")
bad_summary.get("on_output_fill_poly2_rows")
bad_summary.get("on_output_fill_linear_rows")
bad_summary.get("on_output_fill_poly2_runs_total")
bad_summary.get("on_output_fill_linear_runs_total")
bad_summary.get("on_output_fill_none_runs_total")

bad_summary.get("bad_channel_map_unused_keys", [])
```

### 13.6 `RAW` の意味

`QC_REF_HOT_SPUR_TOP_CH_RAW` などの `RAW` は、`ch_range` や `vlsrk_range_kms` で切り出す前の **raw 0-based channel index** です。

---

## 14. 推奨設定

### 14.1 まずの標準推奨

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_apply_flagrow=False,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    qc_spur_on_apply_flagrow=False,
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
)
```

意味:

- まず件数だけ見る
- ON は見ない
- fixed bad channel は NaN のまま通す
- 出力側も埋めない

### 14.2 HOT / OFF を実際に落とす標準設定

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    qc_spur_on_apply_flagrow=False,
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
)
```

### 14.3 1 ch spur を主に見たい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.5,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=1,
    qc_spur_consider_on=False,
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
)
```

### 14.4 2 ch から 3 ch spur も見たい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
)
```

### 14.5 弱い spur を取りこぼしたくない場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_apply_flagrow=False,
    qc_spur_sigma=6.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
)
```

### 14.6 fixed bad channel が既知の場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channels=[1536, 3072],
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
)
```

### 14.7 fixed bad channel を range で切りたい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_ranges=[(1534, 1538)],
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=1,
    qc_spur_consider_on=False,
)
```

### 14.8 group ごとに異なる fixed bad channel がある場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_key_columns=("FDNUM", "IFNUM", "PLNUM"),
    bad_channel_map={
        (2, 0, 0): [(1534, 1538)],
        (3, 0, 1): [3072],
    },
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=1,
    qc_spur_consider_on=False,
)
```

### 14.9 本書の基本推奨: NaN, NaN

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_ranges=[(6767, 6776)],
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=2,
    qc_spur_consider_on=False,
)
```

本書ではこれを標準とします。理由は、10 ch 前後の gap は補間で自然に回復できる保証がなく、OTF-like 平均後に dip が目立つことがあるからです。

### 14.10 quicklook 用に最終 Ta を線形で埋めたい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_ranges=[(6767, 6770)],
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="interp_linear_ta",
    output_bad_channel_fill_max_consecutive_channels=4,
    output_bad_channel_fill_min_valid_neighbors=2,
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=2,
    qc_spur_consider_on=False,
)
```

これは **viewer / quicklook 用**です。science の既定ではありません。

### 14.11 10 ch gap に対して `poly2_weighted_ta` を試す標準形

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_ranges=[(6767, 6776)],
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="poly2_weighted_ta",
    output_bad_channel_fill_max_consecutive_channels=10,
    output_bad_channel_fill_min_valid_neighbors=2,
    output_bad_channel_fill_half_window_channels=8,
    output_bad_channel_fill_curvature_sigma_threshold=1.0,
    output_bad_channel_fill_edge_consistency_sigma_max=3.0,
    output_bad_channel_fill_fallback_policy="linear",
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=2,
    qc_spur_consider_on=False,
)
```

これは「10 ch gap をどうしても viewer 上なめらかに見たい」「`interp_linear_ta` の dip を少しでも減らしたい」ときの試験設定です。

### 14.12 より保守的な `poly2_weighted_ta`

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_ranges=[(6767, 6776)],
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="poly2_weighted_ta",
    output_bad_channel_fill_max_consecutive_channels=10,
    output_bad_channel_fill_min_valid_neighbors=2,
    output_bad_channel_fill_half_window_channels=8,
    output_bad_channel_fill_curvature_sigma_threshold=1.5,
    output_bad_channel_fill_edge_consistency_sigma_max=2.0,
    output_bad_channel_fill_fallback_policy="none",
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=2,
    qc_spur_consider_on=False,
)
```

これは「2 次がかなりもっともらしいときだけ使い、だめなら埋めない」設定です。

### 14.13 `interp_linear_ta` と `poly2_weighted_ta` を比較したい場合

```python
sc_cal_lin = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_ranges=[(6767, 6776)],
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="interp_linear_ta",
    output_bad_channel_fill_max_consecutive_channels=10,
    output_bad_channel_fill_min_valid_neighbors=2,
)

sc_cal_poly = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_ranges=[(6767, 6776)],
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="poly2_weighted_ta",
    output_bad_channel_fill_max_consecutive_channels=10,
    output_bad_channel_fill_min_valid_neighbors=2,
    output_bad_channel_fill_half_window_channels=8,
    output_bad_channel_fill_curvature_sigma_threshold=1.0,
    output_bad_channel_fill_edge_consistency_sigma_max=3.0,
    output_bad_channel_fill_fallback_policy="linear",
)
```

### 14.14 raw 側で補間したい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    bad_channel_ranges=[(6767, 6770)],
    bad_channel_policy="interp_linear",
    bad_channel_interp_max_consecutive_channels=4,
    bad_channel_interp_min_valid_neighbors=2,
    output_bad_channel_fill_policy="none",
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=2,
    qc_spur_consider_on=False,
)
```

ただし、吸収っぽい人工構造が見えたら、raw 補間ではなく `bad_channel_policy="nan"` に戻す方が安全です。

### 14.15 ON を report only で見る場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=True,
    qc_spur_on_apply_flagrow=False,
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
)
```

### 14.16 ON も自動で落としたい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=7.5,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=1,
    qc_spur_consider_on=True,
    qc_spur_on_apply_flagrow=True,
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
)
```

### 14.17 銀河中心や line-rich source の場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    vlsrk_range_kms=(-350, 350),
    qc_sigma=5.0,
    qc_apply_flagrow=True,
    qc_spur_sigma=10.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=1,
    qc_spur_consider_on=False,
    bad_channel_policy="nan",
    output_bad_channel_fill_policy="none",
)
```

line-rich source では、まず ON を切るのが安全です。

### 14.18 `qc_spur_sigma` が効きすぎる / 効かなすぎる場合

- 効きすぎるとき
  - `qc_spur_sigma` を上げる
  - `qc_spur_max_consecutive_channels=1` にする
  - `qc_spur_window_channels=31` に下げる
- それでも変わらないとき
  - fixed bad channel 指定が不十分な可能性が高い
  - `bad_channel_ranges` を広げる
  - peak channel summary を確認する

### 14.19 ON pre-screen の内部値

公開引数ではありませんが、summary8 の内部固定値は

- `Q_D` threshold: $T_D = 8$
- ON 孤立性比: $Q_{\mathrm{iso,min}} = 3.0$

です。

### 14.20 `poly2_weighted_ta` の結果確認

```python
cols = [
    "BADCHAN_OUTFILL_POLICY",
    "BADCHAN_OUTFILL_APPLIED",
    "BADCHAN_OUTFILL_NFILLED",
    "BADCHAN_OUTFILL_HAS_EDGE_UNFILLED",
    "BADCHAN_OUTFILL_POLY2_APPLIED",
    "BADCHAN_OUTFILL_LINEAR_APPLIED",
    "BADCHAN_OUTFILL_POLY2_RUNS",
    "BADCHAN_OUTFILL_LINEAR_RUNS",
    "BADCHAN_OUTFILL_NONE_RUNS",
    "BADCHAN_OUTFILL_CURV_SCORE_MAX",
    "BADCHAN_OUTFILL_EDGE_SCORE_MAX",
    "BADCHAN_OUTFILL_LOCAL_SIGMA_MED",
]

sc_cal.table[cols].iloc[0]
```

見方:

- `BADCHAN_OUTFILL_POLY2_APPLIED=True`
  - 少なくとも 1 本の run で 2 次が採用された
- `BADCHAN_OUTFILL_LINEAR_APPLIED=True`
  - fallback linear が使われた
- `BADCHAN_OUTFILL_NONE_RUNS>0`
  - 埋めなかった run がある

---

## 15. 実務上の推奨手順

### Step 1

fixed bad channel が分かっているなら最初に指定します。

### Step 2

まず `qc_apply_flagrow=False`、`bad_channel_policy="nan"`、`output_bad_channel_fill_policy="none"` で件数だけを見ます。

### Step 3

HOT / OFF の shape / spur 件数を確認します。

### Step 4

HOT / OFF が大量に落ちるなら、peak channel summary を見て fixed bad channel を見直します。

### Step 5

ON を見る必要があるときだけ `qc_spur_consider_on=True` にします。

### Step 6

ON は最初は `qc_spur_on_apply_flagrow=False` にします。

### Step 7

science 用では、そのまま `output_bad_channel_fill_policy="none"` を維持します。

### Step 8

どうしても viewer / quicklook の連続性が必要なときだけ、まず `interp_linear_ta` を試します。

### Step 9

10 ch 前後の gap で line center 付近の dip が問題なら、`poly2_weighted_ta` を試します。

### Step 10

`BADCHAN_OUTFILL_*` 列と `history["bad_channel_summary"]` を見て、実際に 2 次が使われたのか、linear に戻ったのか、none のままかを確認します。

### Step 11

設定が妥当なら、最後に `qc_apply_flagrow=True` で本番を回します。

---

## 16. よくある症状と切り分け

### 16.1 `No HOT rows remain after applying ...`

意味:

- HOT が shape QC または spur QC で全滅した

確認:

```python
qc_summary = sc_cal.history.get("qc_summary", {})
qc_summary.get("hot_shape_auto_bad_would_flag", 0)
qc_summary.get("hot_spur_auto_bad_would_flag", 0)
```

多くの場合、犯人は `hot_spur_auto_bad` です。

### 16.2 `qc_spur_sigma` を上げても件数が変わらない

考えられる原因:

- fixed bad channel が十分に除外されていない
- 1 ch 指定が狭すぎる
- raw index がずれている
- peak channel 周辺が非常に強い

対処:

- peak summary を見る
- `bad_channel_ranges` で少し広めに切る
- raw 0-based index を再確認する

### 16.3 viewer で補間が効いていないように見える

確認点:

- `bad_channel_policy` ではなく `output_bad_channel_fill_policy` を見ているか
- `output_bad_channel_fill_max_consecutive_channels` を超えていないか
- edge gap で左右有効点が不足していないか
- `BADCHAN_OUTFILL_APPLIED` が本当に True か

### 16.4 raw 補間で吸収っぽく見える

これは raw ON / OFF / HOT を別々に補間してから非線形比を作ることによる人工構造の可能性があります。

対処:

- `bad_channel_policy="nan"` に戻す
- 出力補間を使うなら Ta 側に限定する

### 16.5 `interp_linear_ta` では 10 ch gap に deep dip が出る

考えられる原因:

- line center や shoulder のように局所曲率が強い場所を直線で橋渡ししている
- OTF-like に位置の少し違う row を平均している
- 各 row の補間誤差が同符号に偏って残っている

対処:

- まず `output_bad_channel_fill_policy="none"` で NaN のまま比較する
- どうしても埋めるなら `poly2_weighted_ta` を試す
- それでも不安なら埋めない

### 16.6 `poly2_weighted_ta` を使っても改善しない

考えられる原因:

- gap が長すぎる
- line-free ではなく ripple 優勢
- gap が edge に近く左右近傍が不足
- `curvature_sigma_threshold` が高すぎて 2 次が採用されていない
- `edge_consistency_sigma_max` が厳しすぎて linear / none に落ちている

確認:

```python
sc_cal.table[[
    "BADCHAN_OUTFILL_POLY2_APPLIED",
    "BADCHAN_OUTFILL_LINEAR_APPLIED",
    "BADCHAN_OUTFILL_POLY2_RUNS",
    "BADCHAN_OUTFILL_LINEAR_RUNS",
    "BADCHAN_OUTFILL_NONE_RUNS",
    "BADCHAN_OUTFILL_CURV_SCORE_MAX",
    "BADCHAN_OUTFILL_EDGE_SCORE_MAX",
]].iloc[0]
```

### 16.7 `poly2_weighted_ta` で line-free がむしろ不自然

line-free + ripple の場所では、2 次が常に良いとは限りません。

対処:

- `output_bad_channel_fill_fallback_policy="none"` にする
- `output_bad_channel_fill_curvature_sigma_threshold` を上げる
- science では `none` に戻す

---

## 17. 実装上の補足

### 17.1 `qc_spur_sigma` の内部基準

trial 補正の基準有効チャネル数は

$$
N_{\mathrm{ref}} = 512
$$

です。したがって、512 ch 付近では旧版とかなり近く、数千 ch でのみ少し厳しくなります。

### 17.2 1000 ch 以下での解釈

1000 ch 以下では旧方式でも大きな破綻は出にくく、最新版との差は大きくなりにくいです。

### 17.3 peak channel summary

最新版では HOT / OFF の narrow spur に対して

- どの raw channel が原因だったか
- その channel が何回出たか

を summary で確認できます。

これは fixed bad channel の切り分けに非常に有効です。

---

## 18. まとめ

最新版 `calibrate_updated_2026-04-11_statewise_peaksummary8.py` の QC は、次のように整理されています。

### 共通前処理

- fixed bad channel を channel-level で処理する
- bad channel の指定は raw 0-based channel index
- QC 用配列、shape 用配列、較正用配列、最終 Ta 補間を分ける

### HOT / OFF

- `qc_sigma` による shape QC
- `qc_spur_sigma` による narrow spur QC
- narrow spur は science window 内で trial 数補正付き threshold を使う
- peak channel summary で犯人 channel を確認できる

### ON

- `$Q_D$` による pre-screen
- narrow spur 本判定
- 近傍 ON row との孤立性判定
- 既定は report only

### bad channel 補間

- raw では `bad_channel_policy` が効く
- 推奨は `bad_channel_policy="nan"`
- science の基本は `output_bad_channel_fill_policy="none"`
- 見た目の連続性だけが必要なら `output_bad_channel_fill_policy="interp_linear_ta"`
- 10 ch 前後の gap で line center の dip 低減を試すなら `output_bad_channel_fill_policy="poly2_weighted_ta"`

この構成により、

- science window をそのまま使いながら
- fixed bad channel と row 異常を分離し
- HOT / OFF と ON にそれぞれ向いた判定を使い
- ON では高速化も得て
- 出力では summary と Ta 補間まで追跡できる

という設計になっています。


---

## 19. `peaksummary8` との完全対応確認で追加した追補

この節は、最新版コード `calibrate_updated_2026-04-11_statewise_peaksummary8.py` と本書を 1 項目ずつ突き合わせた結果、追記が必要だと確認できた事項をまとめたものです。

本節の追加によって、本書は少なくとも次の観点で `peaksummary8` に対応します。

- `run_tastar_calibration(...)` の公開引数
- `history["qc_summary"]` の key
- `history["bad_channel_summary"]` の key
- `table` に保存される summary 列
- 実装上の内部固定値

### 19.1 `run_tastar_calibration(...)` の公開引数の完全列挙

本書の前半では QC 関連引数を中心に説明しているが、最新版 `peaksummary8` の `run_tastar_calibration(...)` は次の公開引数を持つ。

#### 19.1.1 入力・出力・行選択

- `input_data`
- `output_path`
- `rows`
- `exclude_rows`
- `overwrite`
- `verbose`

#### 19.1.2 スペクトル軸・WCS・速度補正

- `ch_range`
- `vlsrk_range_kms`
- `coord_frame`
- `spectrum_column`
- `store_freq_column`
- `vcorr_chunk_sec`
- `dtype`
- `rest_freq`
- `v_corr_col`

#### 19.1.3 温度較正・大気モデル

- `t_hot_k`
- `tau_zenith`
- `t_surface_k`
- `t_atm_model`
- `t_atm_delta_k`
- `t_atm_eta`
- `gain_mode`

#### 19.1.4 QC・bad channel・最終 Ta 補間

- `qc_sigma`
- `qc_spur_sigma`
- `qc_spur_window_channels`
- `qc_spur_max_consecutive_channels`
- `qc_spur_consider_on`
- `qc_spur_on_apply_flagrow`
- `bad_channels`
- `bad_channel_ranges`
- `bad_channel_key_columns`
- `bad_channel_map`
- `bad_channel_policy`
- `bad_channel_interp_max_consecutive_channels`
- `bad_channel_interp_min_valid_neighbors`
- `output_bad_channel_fill_policy`
- `output_bad_channel_fill_max_consecutive_channels`
- `output_bad_channel_fill_min_valid_neighbors`
- `output_bad_channel_fill_half_window_channels`
- `output_bad_channel_fill_curvature_sigma_threshold`
- `output_bad_channel_fill_edge_consistency_sigma_max`
- `output_bad_channel_fill_fallback_policy`
- `qc_apply_flagrow`

### 19.2 `history["qc_summary"]` の全 key

最新版 `peaksummary8` では、`history["qc_summary"]` に次が入る。

- `hot_input_flagrow_bad`
- `off_input_flagrow_bad`
- `on_input_flagrow_bad`
- `hot_auto_bad_would_flag`
- `off_auto_bad_would_flag`
- `hot_shape_auto_bad_would_flag`
- `off_shape_auto_bad_would_flag`
- `hot_spur_auto_bad_would_flag`
- `off_spur_auto_bad_would_flag`
- `hot_spur_peak_channel_counts_raw`
- `off_spur_peak_channel_counts_raw`
- `hot_spur_peak_channel_top_raw`
- `off_spur_peak_channel_top_raw`
- `on_spur_prescreen_candidate`
- `on_spur_auto_bad_would_flag`
- `hot_used_rows`
- `off_used_rows`

#### 19.2.1 代表的な確認コード

```python
qc_summary = sc_cal.history.get("qc_summary", {})

qc_summary.get("hot_input_flagrow_bad", 0)
qc_summary.get("off_input_flagrow_bad", 0)
qc_summary.get("on_input_flagrow_bad", 0)

qc_summary.get("hot_auto_bad_would_flag", 0)
qc_summary.get("off_auto_bad_would_flag", 0)
qc_summary.get("hot_shape_auto_bad_would_flag", 0)
qc_summary.get("off_shape_auto_bad_would_flag", 0)
qc_summary.get("hot_spur_auto_bad_would_flag", 0)
qc_summary.get("off_spur_auto_bad_would_flag", 0)

qc_summary.get("hot_spur_peak_channel_counts_raw", {})
qc_summary.get("off_spur_peak_channel_counts_raw", {})
qc_summary.get("hot_spur_peak_channel_top_raw", [])
qc_summary.get("off_spur_peak_channel_top_raw", [])

qc_summary.get("on_spur_prescreen_candidate", 0)
qc_summary.get("on_spur_auto_bad_would_flag", 0)
qc_summary.get("hot_used_rows", 0)
qc_summary.get("off_used_rows", 0)
```

### 19.3 `history["bad_channel_summary"]` の全 key

最新版 `peaksummary8` では、`history["bad_channel_summary"]` に次が入る。

- `policy`
- `interp_max_consecutive_channels`
- `output_fill_policy`
- `output_fill_max_consecutive_channels`
- `output_fill_min_valid_neighbors`
- `output_fill_half_window_channels`
- `output_fill_curvature_sigma_threshold`
- `output_fill_edge_consistency_sigma_max`
- `output_fill_fallback_policy`
- `global_bad_channels_raw_count`
- `shape_interp_applied_rows`
- `bad_channel_group_cols`
- `bad_channel_key_columns`
- `bad_channel_map_keys`
- `bad_channel_map_unused_keys`
- `on_total_bad_channels_max`
- `off_total_bad_channels_max`
- `hot_total_bad_channels_max`
- `on_interp_rows`
- `off_interp_rows`
- `hot_interp_rows`
- `on_interp_nfilled_total`
- `off_interp_nfilled_total`
- `hot_interp_nfilled_total`
- `on_has_edge_unfilled_any`
- `off_has_edge_unfilled_any`
- `hot_has_edge_unfilled_any`
- `on_output_fill_rows`
- `on_output_fill_nfilled_total`
- `on_output_fill_has_edge_unfilled_any`
- `on_output_fill_poly2_rows`
- `on_output_fill_linear_rows`
- `on_output_fill_poly2_runs_total`
- `on_output_fill_linear_runs_total`
- `on_output_fill_none_runs_total`

#### 19.3.1 代表的な確認コード

```python
bad_summary = sc_cal.history.get("bad_channel_summary", {})

bad_summary.get("policy", "")
bad_summary.get("interp_max_consecutive_channels", 0)
bad_summary.get("output_fill_policy", "")
bad_summary.get("output_fill_max_consecutive_channels", 0)
bad_summary.get("output_fill_min_valid_neighbors", 0)
bad_summary.get("output_fill_half_window_channels")
bad_summary.get("output_fill_curvature_sigma_threshold")
bad_summary.get("output_fill_edge_consistency_sigma_max")
bad_summary.get("output_fill_fallback_policy", "")

bad_summary.get("global_bad_channels_raw_count", 0)
bad_summary.get("shape_interp_applied_rows", 0)

bad_summary.get("bad_channel_group_cols", [])
bad_summary.get("bad_channel_key_columns", [])
bad_summary.get("bad_channel_map_keys", 0)
bad_summary.get("bad_channel_map_unused_keys", [])

bad_summary.get("on_total_bad_channels_max", 0)
bad_summary.get("off_total_bad_channels_max", 0)
bad_summary.get("hot_total_bad_channels_max", 0)

bad_summary.get("on_interp_rows", 0)
bad_summary.get("off_interp_rows", 0)
bad_summary.get("hot_interp_rows", 0)

bad_summary.get("on_interp_nfilled_total", 0)
bad_summary.get("off_interp_nfilled_total", 0)
bad_summary.get("hot_interp_nfilled_total", 0)

bad_summary.get("on_has_edge_unfilled_any", False)
bad_summary.get("off_has_edge_unfilled_any", False)
bad_summary.get("hot_has_edge_unfilled_any", False)

bad_summary.get("on_output_fill_rows", 0)
bad_summary.get("on_output_fill_nfilled_total", 0)
bad_summary.get("on_output_fill_has_edge_unfilled_any", False)
bad_summary.get("on_output_fill_poly2_rows", 0)
bad_summary.get("on_output_fill_linear_rows", 0)
bad_summary.get("on_output_fill_poly2_runs_total", 0)
bad_summary.get("on_output_fill_linear_runs_total", 0)
bad_summary.get("on_output_fill_none_runs_total", 0)
```

### 19.4 `table` に保存される summary 列の完全列挙

本書の前半でも主要列は説明したが、最新版 `peaksummary8` では summary 関連として少なくとも次の列が `table` に保存される。

#### 19.4.1 bad channel 関連

- `BADCHAN_POLICY`
- `BADCHAN_N_GLOBAL`
- `BADCHAN_N_GROUP`
- `BADCHAN_N_TOTAL`
- `BADCHAN_INTERP_APPLIED`
- `BADCHAN_INTERP_NFILLED`
- `BADCHAN_HAS_EDGE_UNFILLED`
- `BADCHAN_SHAPE_INTERP_APPLIED`
- `BADCHAN_OUTFILL_POLICY`
- `BADCHAN_OUTFILL_APPLIED`
- `BADCHAN_OUTFILL_NFILLED`
- `BADCHAN_OUTFILL_HAS_EDGE_UNFILLED`
- `BADCHAN_OUTFILL_POLY2_APPLIED`
- `BADCHAN_OUTFILL_LINEAR_APPLIED`
- `BADCHAN_OUTFILL_POLY2_RUNS`
- `BADCHAN_OUTFILL_LINEAR_RUNS`
- `BADCHAN_OUTFILL_NONE_RUNS`
- `BADCHAN_OUTFILL_CURV_SCORE_MAX`
- `BADCHAN_OUTFILL_EDGE_SCORE_MAX`
- `BADCHAN_OUTFILL_LOCAL_SIGMA_MED`

#### 19.4.2 HOT / OFF summary

- `QC_REF_HOT_AUTO_BAD`
- `QC_REF_OFF_AUTO_BAD`
- `QC_REF_HOT_SHAPE_AUTO_BAD`
- `QC_REF_OFF_SHAPE_AUTO_BAD`
- `QC_REF_HOT_SPUR_AUTO_BAD`
- `QC_REF_OFF_SPUR_AUTO_BAD`
- `QC_REF_HOT_SPUR_TOP_CH_RAW`
- `QC_REF_HOT_SPUR_TOP_COUNT`
- `QC_REF_OFF_SPUR_TOP_CH_RAW`
- `QC_REF_OFF_SPUR_TOP_COUNT`
- `QC_REF_HOT_SPUR_PEAK_COUNTS_RAW_JSON`
- `QC_REF_OFF_SPUR_PEAK_COUNTS_RAW_JSON`
- `QC_REF_HOT_SPUR_PEAK_TOP_RAW_JSON`
- `QC_REF_OFF_SPUR_PEAK_TOP_RAW_JSON`

#### 19.4.3 ON summary

- `QC_ON_PRESCREEN_CANDIDATE`
- `QC_ON_SPUR_AUTO_BAD`
- `REF_QC_ON_PRESCREEN_QD`
- `REF_QC_ON_PRESCREEN_CANDIDATE`
- `REF_QC_ON_SPUR_SCORE`
- `REF_QC_ON_SPUR_ABSPEAK`
- `REF_QC_ON_SPUR_RUNLEN`
- `REF_QC_ON_SPUR_THRESHOLD`
- `REF_QC_ON_NEIGHBOR_SCORE`
- `REF_QC_ON_ISO_SCORE`
- `REF_QC_ON_SPUR_AUTO_BAD`

#### 19.4.4 QC 設定列

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

### 19.5 内部固定値

最新版 `peaksummary8` では、公開引数ではなく内部固定値として次を使う。

- ON pre-screen threshold:
  $$
  T_D = 8
  $$
- ON isolation ratio threshold:
  $$
  Q_{\mathrm{iso,min}} = 3.0
  $$
- HOT / OFF narrow spur の trial 補正基準チャネル数:
  $$
  N_{\mathrm{ref}} = 512
  $$

これらは現在の `peaksummary8` でも公開パラメータではない。

### 19.6 `peaksummary8` 対応確認の結論

最新版 `peaksummary8` と本書を突き合わせた結果、本書の本体部だけでは次が不足しやすかった。

- `run_tastar_calibration(...)` の QC 以外の公開引数の完全列挙
- `history["qc_summary"]` の全 key
- `history["bad_channel_summary"]` の全 key
- `poly2_weighted_ta` 関連の table 列と summary key
- `NaN, NaN` を基本とする最新推奨方針

本節を含めることで、本書は **`peaksummary8` に対して、少なくとも公開引数・summary key・summary 列・内部固定値の観点で完全対応**とする。
