# FLAGROW ベース QC 詳細説明書

作成日: 2026-04-01

この文書は、現在合意している `FLAGROW` ベースの QC について、目的、責務分離、数式、各関数での挙動、パラメータ、ログの意味、運用上の注意を詳細に説明するための GitHub 向け Markdown 文書である。

本書の対象は、主に以下の 3 系統の QC である。

1. `run_tastar_calibration()` における HOT/OFF QC
2. `run_velocity_coadd()` における ON QC（主に `mode=scan`）
3. `run_baseline_fit()` における OTF baseline QC

本書では、数式は GitHub Markdown で表示できる形にそろえ、複雑な独自マクロは使わない。

---

## 1. この QC の基本思想

### 1.1 何を `FLAGROW` で表すか

この QC では、**1 本の row 全体を使うかどうか**を `FLAGROW` で表す。

- `FLAGROW = 0`: その row は使用可
- `FLAGROW != 0`: その row は使用しない

ここでいう row とは、1 dump に対応する 1 本のスペクトルである。

本 QC は、**channel 単位の細いマスク**ではなく、**row 単位の採否判定**を目的とする。

### 1.2 何を同じ QC で扱わないか

同じ「異常スペクトル」でも、以下は同じ判定量で扱わない。

1. HOT/OFF の参照スペクトル異常
2. coadd 前の ON スペクトル異常
3. OTF baseline 前の gross bad row

理由は、それぞれで「正常な変動」が異なるからである。

- HOT/OFF では、受信機ゲインや空の状態で全体レベルや全体スケールが変わる
- `mode=scan` の ON では、同一 scan 内なら天体信号は大きく変わらないとみなしやすい
- OTF では、scan 中に pointing が変わるので、row 間で天体信号が普通に変わる

したがって、**責務ごとに QC を分ける**必要がある。

### 1.3 自動 QC の実行と `FLAGROW` 反映は分離する

本 QC では、

- 「QC score を計算して診断する」
- 「その結果を `FLAGROW` に反映して実際に除外する」

を分ける。

公開パラメータとしては、以下を想定する。

- `qc_sigma`
- `qc_apply_flagrow`

論理的には、

自動 QC を実行する条件は次である。

$$
qc\_sigma 
eq \mathrm{None}
$$

自動 QC の結果を `FLAGROW` に反映する条件は次である。

$$
qc\_apply\_flagrow = \mathrm{True}
$$

である。

したがって、

- `qc_sigma=5.0, qc_apply_flagrow=False`

なら、QC は計算するが row は落とさない。

一方、

- `qc_sigma=5.0, qc_apply_flagrow=True`

なら、bad と判定された row に `FLAGROW=1` を立て、後段処理から外す。

---

## 2. 記号、変数、単位の定義

### 2.1 基本記号

- $i$: row index
- $k$: spectral channel index
- $v$: 速度軸上の channel 位置または速度

### 2.2 スペクトル

- $H_i(k)$: HOT row の raw spectrum
- $O_i(k)$: OFF row の raw spectrum
- $S_i(k)$: ON row の raw spectrum
- $T_i(k)$: calibration 後の ON spectrum（Ta*）

### 2.3 window

- $V$: `vwin` で指定した window の和集合
- $L$: `line_vwin` で指定した line 除外 window の和集合
- $W$: QC に実際に使う channel 集合

通常、

$$
W =
\begin{cases}
V \setminus L & \text{if } line\_vwin \text{ is given} \\
V & \text{otherwise}
\end{cases}
$$

とする。

### 2.4 robust 幅

本書では、robust な scatter 指標として、例えば median absolute deviation を用いる。

$$
\mathrm{MAD}(x) = \mathrm{median}(|x - \mathrm{median}(x)|)
$$

ガウス分布の標準偏差スケールに合わせたいときは、

$$
\sigma_{\mathrm{robust}} \approx 1.4826\,\mathrm{MAD}(x)
$$

を使う。

---

## 3. 全体の責務分離

### 3.1 HOT/OFF QC

責務:

- calibration の参照として使ってよい HOT/OFF row を選ぶ

実装場所:

- `run_tastar_calibration()`

重要な前提:

- HOT/OFF は全体上下や全体スケール変動が正常に起こりうる
- したがって、**絶対レベル差そのもの**を異常とはみなさない

### 3.2 ON QC for coadd

責務:

- coadd 前に bad ON row を除く
- 特に `mode=scan` で有効

実装場所:

- `run_velocity_coadd()`

重要な前提:

- 同一 scan 内では、基本的に同じ方向を見ているとみなせる
- ただし scan 番号が変われば空間位置も変わる
- したがって、比較は **scan 単位**で行う

### 3.3 OTF baseline QC

責務:

- OTF で gross bad row を baseline / gridding 前に落とす

実装場所:

- `run_baseline_fit()` の内部

重要な前提:

- OTF では row 間で天体信号が変わる
- したがって、template による shape 比較は使いにくい
- その代わり、**line-free / baseline 窓上の軽量 QC**で gross bad row を判定する
- QC のために baseline を 2 回計算してはならない

---

## 4. HOT/OFF QC の詳細

## 4.1 なぜ単純な template 差ではだめか

一見すると、HOT/OFF row $X_i(k)$ と template $\tilde X(k)$ の差

$$
X_i(k) - \tilde X(k)
$$

の大きさを見れば異常を検出できそうに見える。

しかし HOT/OFF では、受信機ゲインや空の状況により、正常でも

- 全体レベルが上がる
- 全体レベルが下がる
- 全体スケールが少し変わる

ということが起こる。

したがって、単純に

$$
\mathrm{RMS}(X_i - \tilde X)
$$

を見ると、正常な全体変動まで異常とみなしてしまう。

## 4.2 HOT/OFF QC の考え方

HOT/OFF では、**offset と scale を吸収した後の残差 shape** を見る。

同一 group の HOT または OFF に対して robust template $\tilde X(k)$ を作る。

各 row $X_i(k)$ に対して、

$$
X_i(k) \approx a_i + b_i\,\tilde X(k)
$$

を近似し、残差を

$$
r_i(k) = X_i(k) - a_i - b_i\,\tilde X(k)
$$

とする。

その row の shape 異常 score を

$$
Q_i^{\mathrm{ref,shape}} = \sigma_{\mathrm{robust}}\bigl(r_i(k)\bigr)
$$

で定義する。

この方法なら、

- 全体上下
- 全体スケール変動

は $a_i, b_i$ に吸収される。

一方、

- 一部 channel だけ壊れている
- 異常な ripple がある
- 狭い異常構造が乗る

といった shape の破綻は $r_i(k)$ に残る。

## 4.3 HOT/OFF QC の group

比較は、少なくとも以下が同じ row 群の中で行う。

- state (`HOT` または `OFF`)
- stream
- beam
- pol
- IF
- spectral setup

## 4.4 HOT/OFF QC の判定式

同一 group 内で $Q_i^{\mathrm{ref,shape}}$ を並べ、

$$
Q_i^{\mathrm{ref,shape}} > \mathrm{median}(Q^{\mathrm{ref,shape}}) + q_c\,\mathrm{MAD}(Q^{\mathrm{ref,shape}})
$$

なら bad HOT/OFF row とする。

ここで $q_c$ は `qc_sigma` である。

## 4.5 HOT/OFF QC の使い方

- `qc_sigma is None`: HOT/OFF auto QC を行わない
- `qc_sigma` が数値: HOT/OFF auto QC を行う
- `qc_apply_flagrow=False`: 診断のみ
- `qc_apply_flagrow=True`: bad HOT/OFF に `FLAGROW=1`

参照平均を作るときは、最終的に `FLAGROW==0` の HOT/OFF のみを使う。

## 4.6 calibrate では絶対閾値が難しい理由

OTF baseline QC では後述の absolute peak 判定を入れているが、同じことを HOT/OFF raw にそのまま入れるのは難しい。

理由は、HOT/OFF raw では正常な level 自体が大きく変わりうるからである。

したがって calibrate では、**絶対 K の閾値**ではなく、**affine 残差の shape 異常**を見るのが自然である。

---

## 5. ON QC for coadd の詳細

## 5.1 この QC の役割

coadd の QC は、**coadd へ入れる前に bad ON row を落とす**ことが目的である。

ここで重要なのは、coadd QC は OTF baseline QC とは別物だという点である。

- coadd QC は、主に `mode=scan` で使う
- OTF baseline QC は、主に OTF の baseline / gridding 前に使う

## 5.2 `mode=position` では既定で必須にしない理由

`mode=position` では、比較できる row 数が少ないことが多い。さらに、位置固定観測では後段での目視や手動 flag の方が安全な場合もある。

したがって、`mode=position` では既定で自動 coadd QC を必須にはしない。

## 5.3 `mode=scan` で shape 比較を使う理由

`mode=scan` では、同一 scan 内では基本的に同じ方向を見ているとみなせる。したがって、同一 scan 内では ON の天体信号 shape は大きく変わらない前提を置ける。

ただし、scan が変われば場所も変わるので、比較は scan ごとに分ける必要がある。

## 5.4 coadd の shape QC

同一 scan group の ON spectrum に対して robust template $\tilde T(k)$ を作る。

各 row $T_i(k)$ に対して、

$$
T_i(k) \approx a_i + b_i\,\tilde T(k)
$$

を当て、残差

$$
r_i(k) = T_i(k) - a_i - b_i\,\tilde T(k)
$$

を作る。

その shape score を

$$
Q_i^{\mathrm{scan,shape}} = \sigma_{\mathrm{robust}}\bigl(r_i(k)\bigr)
$$

で定義する。

同一 group 内で、

$$
Q_i^{\mathrm{scan,shape}} > \mathrm{median}(Q^{\mathrm{scan,shape}}) + q_c\,\mathrm{MAD}(Q^{\mathrm{scan,shape}})
$$

なら shape 異常の候補とする。

## 5.5 coadd の absolute peak QC

shape 比較だけでは、狭いが非常に高い異常ピークを取り逃すことがある。

そこで、window $W$ 内で row ごとの中央値を

$$
m_i = \mathrm{median}\{T_i(v) \mid v \in W\}
$$

とし、最大偏差

$$
A_i = \max_{v \in W} |T_i(v) - m_i|
$$

を定義する。

このとき、

$$
A_i > T_{\mathrm{abs,max}}
$$

なら absolute peak 異常とする。

ここで $T_{\mathrm{abs,max}}$ が `qc_abs_peak_k` である。

## 5.6 coadd での最終 bad 判定

`mode=scan` では、最終的に

shape 判定を $S_i$、abs peak 判定を $P_i$ と書くと、coadd の最終判定は次である。

$$
\mathrm{bad\ row} \iff S_i \lor P_i
$$

とする。

つまり、

- group に対して shape が異常
- あるいは window 内に大きすぎる peak がある

のどちらかなら bad 候補である。

## 5.7 coadd で使う window

coadd で absolute peak を見るときの $W$ は、原則として

$$
W = V \setminus L
$$

である。

ここで、

- $V$: baseline 用または RMS 用の `vwin`
- $L$: `line_vwin`

である。

実際の実装では、典型的には

1. `baseline_vwin - line_vwin`
2. なければ `rms_vwin - line_vwin`

の順に候補を探す。

## 5.8 coadd で `line_vwin` が必要な場合

- `vwin` 自体が line-free なら `line_vwin` は不要
- `vwin` が広く、その中に線があるなら `line_vwin` が必要

なぜなら、本物の線が $W$ に入ると、absolute peak 判定はそれも異常とみなすからである。

## 5.9 coadd での regrid の影響

coadd 側の QC は、実装上、coadd 内部で扱う regridded spectrum 上で評価されることがある。

その場合、1 channel の非常に鋭いスパイクは、regrid により振幅がやや薄まる可能性がある。

したがって、同じ `qc_abs_peak_k` を baseline 側と coadd 側にそのまま共通で使うと、coadd 側では効き方が弱く見えることがある。

このため、運用上は

- baseline 側の `qc_abs_peak_k`
- coadd 側の `qc_abs_peak_k`

を必ずしも同じ値にしない方がよい場合がある。

---

## 6. OTF baseline QC の詳細

## 6.1 この QC の役割

OTF baseline QC は、**OTF の baseline / gridding に入れる前に gross bad row を落とす**ための QC である。

重要なのは、これは **`run_baseline_fit()` の内部で行う**という点である。

つまり、別の前処理関数を増やすのではなく、OTF では `run_baseline_fit()` を必ず通る前提で、その内部で screening する。

## 6.2 なぜ template 比較を使わないか

OTF では scan 中に pointing が連続的に変わる。したがって、row 間で天体信号 shape が変わるのが正常である。

このため、`mode=scan` の coadd のような template shape 比較は、OTF では原則として使わない。

## 6.3 なぜ baseline を 2 回計算しないか

QC のために baseline fit をもう 1 回走らせると、計算量が単純に増えるだけでなく、責務も曖昧になる。

したがって、OTF QC では以下を守る。

- baseline fit は本体で 1 回だけ
- QC のために baseline を再計算しない
- その代わり、軽量な row 指標で gross bad を判定する

## 6.4 OTF の diff-MAD 判定

window $W$ 上で一次差分

$$
\Delta T_i(v) = T_i(v+1) - T_i(v)
$$

を作る。

その robust 幅を

$$
D_i = 1.4826\,\mathrm{MAD}\{\Delta T_i(v) \mid v, v+1 \in W\}
$$

と定義する。

この量は、

- 定数オフセットに不変
- 一次傾きにほぼ不変
- スパイク、ギザギザ、異常ノイズに敏感

という性質を持つ。

同一 group 内で、

$$
D_i > \mathrm{median}(D) + q_c\,\mathrm{MAD}(D)
$$

なら diff-MAD 異常候補とする。

## 6.5 OTF の absolute peak 判定

今回の実装では、diff-MAD だけでは取り逃すケースに対応するため、absolute peak 判定を追加している。

window $W$ 内で row ごとの中央値を

$$
m_i = \mathrm{median}\{T_i(v) \mid v \in W\}
$$

とし、

$$
A_i = \max_{v \in W} |T_i(v) - m_i|
$$

を定義する。

このとき、

$$
A_i > T_{\mathrm{abs,max}}
$$

なら absolute peak 異常とする。

ここで $T_{\mathrm{abs,max}}$ は `qc_abs_peak_k` である。

## 6.6 OTF での最終 bad 判定

OTF baseline QC の最終判定は OR である。ここで、diff-MAD 判定を $D_i$、abs peak 判定を $P_i$ と書く。

$$
\mathrm{bad\ row} \iff D_i \lor P_i
$$

つまり、

- row 全体が荒れている
- または一部に大きすぎる peak がある

のどちらかなら bad 候補とする。

## 6.7 OTF で使う window

OTF では、原則として line-free / baseline 用の窓を使う。

- `vwin` が line-free だけを表すなら、それをそのまま使ってよい
- `vwin` が広く、その中に線があるなら `line_vwin` を与え、

$$
W = V \setminus L
$$

を使う

もし line-free 指定の結果 $W$ が空になるなら、その row / group に対して QC を無理に行わない方が安全である。全帯域へ自動で勝手にフォールバックしてはならない。

## 6.8 OTF と line-free の関係

OTF の absolute peak 判定は、本物の天体線が $W$ に入っていないことを前提にしている。

したがって、OTF では **line-free の指定が重要**である。

line-free が不適切だと、

- 本物の線を異常 peak とみなす
- 本物の線の周辺で diff-MAD が増える

という誤判定につながる。

---

## 7. `vwin` と `line_vwin` の考え方

この QC で非常に重要なのが、`vwin` と `line_vwin` の役割の違いである。

## 7.1 `vwin` が line-free でない場合

例えば、

- `vwin = ["-30:55"]`
- `line_vwin = ["5:15"]`

なら、baseline fit も QC も、実際に使う window は

$$
F = V \setminus L
$$

である。

つまり、fit に使う channel も、QC に使う channel も、`vwin` のうち line を除いた範囲である。

baseline model $B_i(v)$ の係数は

$$
\hat\theta_i = \arg\min_{\theta_i}
\sum_{v \in F,\,\mathrm{finite}}
\left[T_i(v) - B_i(v;\theta_i)\right]^2
$$

で決まり、求まった baseline は全帯域に適用して差し引く。

## 7.2 `vwin` 自体が line-free の場合

例えば、

- `vwin = ["-30:-5", "15:55"]`
- `line_vwin = None`

なら、`vwin` 自体が line-free である。

この場合、baseline fit にも QC にも、そのまま `vwin` を使えばよい。

---

## 8. パラメータ一覧と意味

## 8.1 共通

### `qc_sigma`

- 型: `None` または数値
- 意味: 自動 QC の閾値強度

役割:

- `None`: 自動 QC 無効
- 数値: group 内の robust 外れ判定に使用

## 8.2 共通

### `qc_apply_flagrow`

- 型: `bool`
- 意味: 自動 QC の結果を `FLAGROW` に反映するか

役割:

- `False`: 診断のみ
- `True`: bad row に `FLAGROW=1`

## 8.3 coadd / baseline

### `qc_abs_peak_k`

- 型: `None` または数値 [K]
- 意味: absolute peak 判定の閾値

役割:

- `None`: abs peak 判定を使わない
- 数値: 窓内最大偏差がこの値を超えた row を bad 候補にする

注意:

- baseline 側と coadd 側で同じ値が最適とは限らない
- coadd 側では regrid の影響でピークが少し薄まることがある

---

## 9. ログ、history、出力列の意味

現在の実装では、QC の有無が画面でも追えるように、summary print を出す方針を取っている。

## 9.1 calibrate の出力例

```text
[CALIBRATE] HOT/OFF QC: hot_auto_bad=0, off_auto_bad=16, apply_flagrow=False
```

意味:

- `hot_auto_bad`: HOT で自動 QC により bad 候補になった row 数
- `off_auto_bad`: OFF で自動 QC により bad 候補になった row 数
- `apply_flagrow`: 実際に `FLAGROW` を立てたかどうか

## 9.2 baseline の出力例

```text
[BASELINE] OTF gross QC: auto_bad=12, diffmad_bad=4, abspeak_bad=10, applied=0, dropped_input_flagrow=3, dropped_flagrow=3, window_sources=['baseline_minus_line'], qc_abs_peak_k=100.0, apply_flagrow=False
```

意味:

- `auto_bad`: 自動 QC 全体で bad 候補になった row 数
- `diffmad_bad`: diff-MAD で bad 候補になった row 数
- `abspeak_bad`: abs peak で bad 候補になった row 数
- `applied`: 実際に今回 `FLAGROW` 適用で落とした row 数
- `dropped_input_flagrow`: 入力時点で既に `FLAGROW!=0` のため除外された row 数
- `dropped_flagrow`: 最終的に `FLAGROW` により除外された row 数
- `window_sources`: どの window を QC に使ったか
- `qc_abs_peak_k`: abs peak 閾値
- `apply_flagrow`: 実際に `FLAGROW` に反映したか

## 9.3 coadd の出力例

```text
[COADD] scan QC: would_flag=12, applied=0, dropped_input_flagrow=3, group_mode=scan, scan_groups_used=8, apply_flagrow=False
```

または、実装が進んだ版では abs peak 関連の項目も加わる。

意味:

- `would_flag`: `qc_apply_flagrow=False` でも、もし適用していれば bad だった row 数
- `applied`: 実際に coadd から除外した row 数
- `dropped_input_flagrow`: 入力時点で既に `FLAGROW!=0` の row 数
- `group_mode`: `scan` など
- `scan_groups_used`: shape QC を実際に計算した group 数
- `qc_abs_peak_k`: abs peak 閾値（設定時）
- `scan_qc_window_source`: abs peak 判定に使った窓の種別

## 9.4 table / history に残す意味

`qc_apply_flagrow=False` の report-only 運用でも、

- `would_flag`
- `auto_bad`
- `diffmad_bad`
- `abspeak_bad`

などを history や summary に残すことで、閾値設定の妥当性を後で検証できる。

---

## 10. 実運用の推奨手順

## 10.1 最初は report-only

最初から `FLAGROW` を立てるのではなく、まずは

- `qc_apply_flagrow=False`

で動かし、

- どれくらい `would_flag` が出るか
- 期待した異常だけが拾われているか
- 本物の線を誤って拾っていないか

を確認する。

## 10.2 次に flag 適用

十分に妥当だと判断したら、

- `qc_apply_flagrow=True`

へ進む。

## 10.3 `qc_abs_peak_k` の調整

`qc_abs_peak_k` は、データの単位、regrid の有無、窓の取り方で効き方が変わる。

したがって、

- baseline 側
- coadd 側

で同一値に固定するより、それぞれの出力と実データを見て調整した方がよい。

## 10.4 OTF では line-free の妥当性を先に疑う

OTF で abs peak や diff-MAD が大量に bad を出す場合、まず疑うべきは

- 本当に異常 row が多いのか
- `line_vwin` が不適切で本物の線を見ていないか

である。

---

## 11. この QC が扱わないもの

本 QC は row 単位の採否判定に重点を置いているので、以下は本書の主対象ではない。

1. 1--2 channel の narrow spur 補間
2. channel 単位のフラグ管理
3. calibrate における絶対閾値の強制判定
4. OTF baseline QC のための baseline 再計算

これらは必要に応じて別機能として考えるべきである。

---

## 12. まとめ

この QC の論理構成を一言で言うと、以下である。

- HOT/OFF は **offset/scale を吸収した shape 異常**を見る
- coadd（`mode=scan`）は **scan 単位の shape 異常 + abs peak** を見る
- OTF baseline は **line-free / baseline 窓上の diff-MAD + abs peak** を見る
- 自動 QC 実行と `FLAGROW` 反映は分離する

すなわち、同じ `FLAGROW` を使いながらも、

- calibrate
- coadd
- baseline

で **正常とみなすべき変動**が違うことを明示的に反映した設計である。

