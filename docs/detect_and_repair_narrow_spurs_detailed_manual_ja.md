# `detect_and_repair_narrow_spurs` 詳細説明書

本書は、現在の `spur_verified_v4.py` に実装されている `detect_and_repair_narrow_spurs` の**実装そのもの**を説明する文書です。  
ここでは「理想仕様」ではなく、**今のコードが実際に何を入力とし、何を計算し、どのような条件で補間し、何を出力するか**を明示します。

---

## 1. 目的

`detect_and_repair_narrow_spurs` は、Ta* に較正された dump 群のスペクトルに対して、**1--2 channel 程度の狭いスプリアス**を検出し、必要に応じて**線形内挿で補間**する関数です。

主な前提は次の通りです。

- 入力は **Ta* 較正後**の `Scantable`
- channel 軸は **native channel** のまま
- 対象は **狭い narrow spur**
- dump 全体を bad spectrum として落とす関数ではない
- `Scantable` 自体の構造は変えず、`data`, `table`, `history` に結果を反映する

したがって、図のような**dump 全体にわたる異常スペクトル**や、**HOT/OFF/ON の row 単位 FLAGROW 判定**とは役割が異なります。

---

## 2. 関数の位置づけ

この関数は、論理的には次のような流れの中で使うことを想定しています。

```text
raw spectrum
  -> HOT/OFF を使った calibration
  -> Ta* dump
  -> detect_and_repair_narrow_spurs
  -> velocity regrid / coadd / baseline / gridding
```

重要なのは、**速度補正や regrid の前**に使うことです。  
理由は、backend 固定の narrow spur は native channel では同じ channel に現れやすく、速度軸へ変換した後では位置が見かけ上ずれて見えるためです。

---

## 3. 入力・出力の定義

### 3.1 入力

関数シグネチャは次の通りです。

```python
sc2 = detect_and_repair_narrow_spurs(
    sc,
    spur_mode="local",
    spur_nsigma=7.0,
    spur_max_width=2,
    verbose=True,
)
```

### 3.2 入力の意味

- `sc`: `Scantable`
- `sc.data[i, k] = T_i(k)`
  - $i = 0, 1, \dots, N_{\rm row}-1$: row index
  - $k = 0, 1, \dots, N_{\rm ch}-1$: native channel index
  - 単位は通常 Ta* [K]
- `spur_mode`
  - `"none"`
  - `"local"`
  - `"promoted-all-dumps"`
- `spur_nsigma`
  - 検出しきい値
- `spur_max_width`
  - 自動補間を許す最大連続 channel 数
- `verbose`
  - summary を表示するかどうか

### 3.3 受け付ける mode 文字列

内部では次が同義語として受理されます。

- `none`, `off`, `false` -> `none`
- `local` -> `local`
- `promoted-all-dumps`, `promoted_all_dumps`, `promoted` -> `promoted-all-dumps`

### 3.4 出力

出力は**新しい `Scantable`**です。  
ただし完全に新規構造を作るのではなく、`sc.copy()` を作った上で次を更新します。

- `data`: 補間後スペクトルに置き換える
- `table`: `SPUR_*` 列を追加する
- `history`: spur 処理に関する履歴情報を追加する

元の `Scantable` をその場で破壊的変更するわけではありません。

---

## 4. 現在の実装が前提としているデータ構造

この関数は現在、次を前提にしています。

1. `sc.data` は `numpy.ndarray`
2. `sc.data.ndim == 2`
3. `len(sc.table) == sc.data.shape[0]`

したがって、例えば可変長スペクトルを Python list で持つような `Scantable.data` には対応していません。

異常系では、次のようにエラーになります。

- `sc.data` が ndarray でない -> `TypeError`
- `sc.data` が 2D でない -> `ValueError`
- `table` の行数が一致しない -> `ValueError`

---

## 5. 公開パラメータの意味

### 5.1 `spur_mode`

#### `none`
何もしません。  
内部的には `sc.copy()` を返すだけです。

#### `local`
**その row で検出された channel だけ**を補間します。

#### `promoted-all-dumps`
まず `local` と同じ検出を全 row に対して行い、  
多くの row に共通して現れた channel を「昇格 channel」として選びます。  
その後、その昇格 channel を**全 row に対して補間対象候補**に加えます。

ただし、ここで重要なのは、

- 1回見つかっただけでは昇格しない
- 昇格した channel も、`spur_max_width` を超える run になれば補間されない

という点です。

### 5.2 `spur_nsigma`

row ごとのロバストな residual scale を $\sigma_i$ としたとき、  
判定条件は概念的には

$$
|r_i(k)| > n_\sigma \, \sigma_i
$$

です。  
ここで $n_\sigma$ が `spur_nsigma` です。

### 5.3 `spur_max_width`

連続して検出された channel run の長さを

$$
m = b - a + 1
$$

としたとき、

$$
m \le m_{\max}
$$

を満たす run のみを自動補間します。  
ここで $m_{\max} =$ `spur_max_width` です。

現在の既定は 2 です。

---

## 6. 内部アルゴリズムの全体像

現在の実装の流れは次の通りです。

1. `data_in = np.asarray(sc.data, dtype=float)` を作る
2. 各 row に対して**5点 median residual**を計算する
3. 各 row の residual からロバストスケール $\sigma_i$ を作る
4. row ごとに local spur run を検出する
5. 全 row の local 検出結果を集計する
6. `promoted-all-dumps` のときだけ昇格 channel を選ぶ
7. row ごとに線形内挿で補間する
8. `SPUR_*` 列と `history` を書く
9. summary を表示する

以下で各段階を詳しく説明します。

---

## 7. 5点 median residual の計算

### 7.1 定義

各 row のスペクトルを $T_i(k)$ とします。  
現在の実装では、channel 方向に**半幅 2**の窓を取り、合計 5 点の median を計算します。

概念的には

$$
\tilde{T}_i(k) = \mathrm{median}\bigl(T_i(k-2), T_i(k-1), T_i(k), T_i(k+1), T_i(k+2)\bigr)
$$

として、局所 residual を

$$
r_i(k) = T_i(k) - \tilde{T}_i(k)
$$

で定義します。

実装上は `sliding_window_view` と `np.nanmedian` を用いて 2D 一括で計算しています。

### 7.2 なぜ 5 点なのか

これは現在の実装で**固定**です。  
公開パラメータではありません。

- 1 ch spike を目立たせやすい
- 広い線や緩やかな baseline をある程度打ち消せる
- 計算が単純で高速

という理由で採用されています。

### 7.3 valid 条件

窓内 finite 値が 3 点未満の場所は residual を `NaN` にします。  
また中心 channel 自体が `NaN` の場合も residual は `NaN` になります。

---

## 8. row ごとのロバストスケール $\sigma_i$

各 row の residual 列 $r_i(k)$ から、現在の実装では

$$
\sigma_i = 1.4826 \times \mathrm{MAD}\bigl(r_i(k)\bigr)
$$

を作ります。  
ここで MAD は median absolute deviation です。

より具体的には

$$
\mathrm{MAD}(r_i) = \mathrm{median}\left( \left| r_i(k) - \mathrm{median}(r_i) \right| \right)
$$

です。

もし $\sigma_i$ が有限でない、あるいは $\sigma_i \le 0$ なら、その row の $\sigma_i$ は一旦 `NaN` とします。

---

## 9. local 検出

### 9.1 候補 channel

各 row について、候補 channel 集合を

$$
C_i = \left\{ k \mid |r_i(k)| > n_\sigma \, \sigma_i \right\}
$$

で作ります。

ここで $n_\sigma =$ `spur_nsigma` です。

### 9.2 連続 run へのまとめ

隣接する候補 channel をまとめて run $(a,b)$ にします。  
run 幅は

$$
m = b - a + 1
$$

です。

### 9.3 幅条件

run 幅が

$$
m > \texttt{spur\_max\_width}
$$

なら、その run は捨てます。

### 9.4 bridge residual による再確認

候補 run $(a,b)$ に対して、外側端点 $a-1$ と $b+1$ を結ぶ線形補間を作ります。

$$
\hat{T}_i(k; a,b)
=
T_i(a-1)
+
\frac{k-(a-1)}{(b+1)-(a-1)}
\left[T_i(b+1)-T_i(a-1)\right]
$$

ただし $a \le k \le b$ です。

その residual を

$$
e_i(k; a,b) = T_i(k) - \hat{T}_i(k; a,b)
$$

とします。

run の score は

$$
\mathrm{score}(a,b) = \frac{\max |e_i(k;a,b)|}{\sigma_i}
$$

です。

これが

$$
\mathrm{score}(a,b) > n_\sigma
$$

を満たす run だけを local 検出結果として受理します。

### 9.5 端の channel が補間されない理由

run が band edge にかかっていて $a=0$ または $b=N_{\rm ch}-1$ のような場合、  
外側端点が作れないので `bridge_segment` が失敗し、その run は補間されません。

したがって現在の実装は、**band edge の spur に対して保守的**です。

---

## 10. fallback sigma

row ごとの $\sigma_i$ が `NaN` になる場合があります。  
例えば、非常に平坦で residual の MAD が 0 になった場合です。

そのとき現在の実装では、有限な $\sigma_i$ の中央値

$$
\sigma_{\rm fb} = \mathrm{median}\{\sigma_i \mid \sigma_i \text{ is finite and } \sigma_i > 0\}
$$

を fallback として使います。

bad row については、

- row 自身の $\sigma_i$ の代わりに $\sigma_{\rm fb}$ を使って local 検出をやり直す
- `SPUR_SIGMA` にはその fallback 値を書き込む

という実装です。

---

## 11. `promoted-all-dumps` の意味

これはしばしば誤解されやすいので、実装の意味を明確に書きます。

### 11.1 集計 counts

まず、local 検出で得られた各 row の local channel 集合を $L_i$ とします。  
各 channel $k$ について、出現回数を

$$
C(k) = \sum_i 1[k \in L_i]
$$

で数えます。

### 11.2 昇格しきい値

昇格しきい値は現在の実装では固定式で

$$
N_{\rm promote} = \max\left(3, \left\lceil 0.05 \, N_{\rm row} \right\rceil \right)
$$

です。

### 11.3 昇格 channel

`promoted-all-dumps` のとき、昇格 channel 集合を

$$
G = \{ k \mid C(k) \ge N_{\rm promote} \}
$$

で定義します。

### 11.4 各 row での補間

各 row ではまず local channels を補間します。  
その後、昇格 channel のうち local で既に補間済みでないものを追加で補間します。

ここで重要なのは、**local 補間と promoted 補間を分けている**ことです。  
これにより、たまたま隣接して 3 ch run になり、local の 1 ch まで巻き添えで未補間になる、という以前の問題を避けています。

### 11.5 よくある誤解

`promoted-all-dumps` は

> 1回でも見つかった channel を全 row で補間する

という意味ではありません。  
現在の実装では、**十分多くの row に出た channel だけ**が昇格します。

したがって、`spur_nsigma` を下げても見た目がほとんど変わらない場合、

- local 検出されていない
- 検出はされているが `N_{\rm promote}` に届いていない
- `spur_max_width` を超えていて補間されていない

のどれかである可能性があります。

---

## 12. 補間の方法

各 run $(a,b)$ に対して、現在の実装は**線形内挿**を使います。

$$
T_i^{\rm corr}(k)
=
T_i(a-1)
+
\frac{k-(a-1)}{(b+1)-(a-1)}
\left[T_i(b+1)-T_i(a-1)\right]
$$

ただし $a \le k \le b$ です。

この補間は、

- 左右両端が finite であること
- run 幅が `spur_max_width` 以下であること

を満たすときだけ実行されます。

---

## 13. 出力 `table` の列

現在の実装では `sc.table` に次の列を追加します。

### `SPUR_N`
その row で最終的に補間された channel 数。

### `SPUR_N_LOCAL`
その row で local として補間された channel 数。

### `SPUR_N_PROM`
その row で promoted として補間された channel 数。

### `SPUR_MODE`
その row の最終状態。

- `NONE`
- `LOCAL`
- `PROMOTED`

注意として、`PROMOTED` は

- promoted だけ補間した row
- local と promoted の両方を補間した row

の両方を含みます。

### `SPUR_SIGMA`
その row に対して使われた $\sigma_i$。  
元の row sigma か、fallback sigma のどちらかです。

### `SPUR_CHAN`
最終的に補間された channel 番号配列。  
`np.int32` の object 配列として格納されます。

### `SPUR_DET_LOCAL`
local 検出段階で見つかった channel 番号配列。  
こちらも `np.int32` の object 配列です。

---

## 14. 出力 `history` の内容

現在の実装では `history` に次を書き込みます。

- `spur_enabled`
- `spur_mode`
- `spur_nsigma`
- `spur_max_width`
- `spur_promote_n`
- `spur_local_detected_total`
- `spur_rows_with_local_detections`
- `spur_counts_nonzero`
- `spur_promoted_channels`
- `spur_total_repaired_channels`
- `spur_total_rows_with_repairs`

### `spur_counts_nonzero`
これは JSON 文字列です。  
中身は

```json
{"123": 4, "456": 7}
```

のように、channel ごとの local 検出回数を表します。

### `spur_promoted_channels`
これも JSON 文字列です。  
昇格 channel 番号の配列を表します。

---

## 15. `verbose=True` の summary

`verbose=True` のとき、最後に summary を 1 行 print します。

例:

```text
[SPUR] mode=local, nsigma=7, max_width=2, rows_detected=2/5, detected_channels_total=3, rows_repaired=2/5, repaired_channels_total=3
```

`promoted-all-dumps` のときは追加で

- `promote_n`
- `promoted_channel_count`
- `promoted_channels`

も表示します。

これは、実際に何本見つかったかを素早く把握するための summary です。

---

## 16. `nsigma` を変えても見た目が変わらないときの確認法

これは実運用上かなり重要です。

`spur_nsigma=7` と `2` で見た目がほとんど同じとき、まず確認すべきは

- 本当に 1 点でも補間したか
- local は見つかっているが promoted に昇格していないのか
- 幅条件で落ちているのか

です。

次のように、補間前後差分を直接確認するのが確実です。

```python
sc0 = sc.copy()
sc1 = detect_and_repair_narrow_spurs(
    sc0,
    spur_mode="promoted-all-dumps",
    spur_nsigma=2.0,
    spur_max_width=2,
    verbose=True,
)

print("spur_total_repaired_channels =", sc1.history.get("spur_total_repaired_channels"))
print("spur_total_rows_with_repairs =", sc1.history.get("spur_total_rows_with_repairs"))
print("spur_promote_n =", sc1.history.get("spur_promote_n"))
print("spur_promoted_channels =", sc1.history.get("spur_promoted_channels"))

d = sc1.data - sc0.data
print("max |delta| =", np.nanmax(np.abs(d)))

rows, chans = np.where(np.abs(d) > 0)
print("n_changed_points =", len(rows))
print("first changed rows =", rows[:20])
print("first changed chans =", chans[:20])
```

ここで

- `max |delta| = 0` なら、全く補間していない
- `SPUR_N > 0` の row があるなら、その row では補間している
- `spur_promoted_channels = []` なら、promoted は起きていない

という解釈になります。

---

## 17. 現在の実装の長所

1. **パラメータが少ない**
   - `spur_mode`
   - `spur_nsigma`
   - `spur_max_width`

2. **高速**
   - 5点 median residual を 2D 一括計算している

3. **`Scantable` 無改造**
   - `data/table/history` のみ変更

4. **local と promoted を分離して補間**
   - 以前の巻き添え問題を避けている

5. **何をしたかが `table/history` に残る**

---

## 18. 現在の実装の限界と注意点

### 18.1 dump 全体異常には向かない
この関数は narrow spur 用です。  
図のような**row 全体の変なスペクトル**は、`FLAGROW` ベースの QC で落とすべきです。

### 18.2 線上に重なった狭い異常には限界がある
広い実線の上に narrow spike が重なったケースでは、ある程度は効く可能性がありますが、  
実線の幅や形に依存します。完全自動判定は困難です。

### 18.3 band edge には保守的
端点外側が必要なので、edge にかかった spur は補間されにくいです。

### 18.4 `promoted-all-dumps` は「少数出現」には効かない
昇格条件

$$
C(k) \ge \max\left(3, \left\lceil 0.05 N_{\rm row} \right\rceil \right)
$$

があるため、少数 dump にだけ出るものは promoted されません。

### 18.5 `FLAG` 列とは連動していない
この実装は channel-wise `FLAG` を持ちません。  
補間済みの値を `data` に書き戻し、履歴を `table/history` に記録する方式です。

---

## 19. 使用上の推奨

### 19.1 まずは `local`
最初は

```python
sc2 = detect_and_repair_narrow_spurs(
    sc,
    spur_mode="local",
    spur_nsigma=7.0,
    spur_max_width=2,
    verbose=True,
)
```

から始めるのが安全です。

### 19.2 `promoted-all-dumps` は、同じ native channel に繰り返し出る場合のみ
backend 固定の bad channel が疑われるときに使います。

### 19.3 `SPUR_*` 列を必ず確認する
少なくとも次は見るべきです。

- `SPUR_N`
- `SPUR_MODE`
- `SPUR_SIGMA`
- `SPUR_CHAN`
- `history["spur_promoted_channels"]`

---

## 20. 実装上の注意

### 20.1 dtype
内部計算は `float` で行います。  
最後に元の `sc.data.dtype` が浮動小数ならその dtype に戻します。  
整数 dtype の場合は unsafe cast を避けるため、そのまま float の `data_out` を入れます。

### 20.2 `table` への保存形式
`SPUR_CHAN`, `SPUR_DET_LOCAL` は object 配列で、その中身は `np.int32` 配列です。  
現行 writer はこのような vector-in-cell を扱える前提です。

---

## 21. まとめ

`detect_and_repair_narrow_spurs` は現在、次の関数です。

- 対象: Ta* 後、native channel の narrow spur
- 入力: 2D ndarray の `Scantable.data`
- 検出: 5点 median residual + row ごとの MAD sigma
- 受理: bridge residual でも `nsigma` を超える 1--`max_width` ch run
- 補間: 線形内挿
- mode:
  - `local`
  - `promoted-all-dumps`
- 出力: 補間済み `data` と `SPUR_*` 記録

したがって、これは**局所狭帯域異常の補間関数**であり、  
**HOT/OFF/ON の row 単位 QC や `FLAGROW` 判定とは別物**です。

両者を混同しないことが重要です。

