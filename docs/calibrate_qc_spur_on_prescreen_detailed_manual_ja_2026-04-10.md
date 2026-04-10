# `qc_spur_consider_on=True` 用 ON 事前選別機能 詳細説明書

## 0. この文書の位置づけ

この文書は、`calibrate.py` に追加した狭帯域 spur QC に対して、`qc_spur_consider_on=True` のときに ON row の判定を高速化するための **事前選別機能** の詳細仕様書です。

重要な点として、**この文書で説明する ON 事前選別機能は、現時点では提案仕様であり、まだ最新版 `calibrate_updated_2026-04-10.py` には実装されていません。**

現状の実装では、`qc_spur_consider_on=True` のとき、ON 全 row に対して HOT/OFF と同じ狭帯域 spur QC をそのまま実行します。そのため、ON row 数が多い観測では時間がかかります。

本機能の目的は、ON row の大部分を軽い判定で早期に通過させ、**狭帯域 spur がありそうな row にだけ重い rolling median ベースの本判定を行う**ことです。

---

## 1. 背景

現在の狭帯域 spur QC は、各 row のスペクトルを周波数方向に rolling median で平滑化し、その残差の大きさと連続チャネル数を使って狭帯域 spur を検出します。

スペクトルを $S(i)$、チャネル番号を $i$ とすると、平滑成分を $B(i)$、残差を

$$
R(i) = S(i) - B(i)
$$

と定義します。

残差の robust な尺度を

$$
\sigma_R = 1.4826 \, \mathrm{MAD}(R)
$$

とし、

$$
|R(i)| > q \, \sigma_R
$$

を満たすチャネル群を異常候補とします。ここで $q$ は `qc_spur_sigma` です。

さらに、連続して超過したチャネル数を $L$ とし、

$$
1 \le L \le L_{\max}
$$

のときに狭帯域 spur とみなします。ここで $L_{\max}$ は `qc_spur_max_consecutive_channels` です。

この本判定は意味としては妥当ですが、`qc_spur_consider_on=True` の場合、ON row 数が多いと計算量が大きくなります。特に rolling median は row ごとに周波数方向へ実行されるため、ON の全 row に対して毎回行うと遅くなります。

---

## 2. 本機能の目的

本機能の目的は、ON row に対して次の 2 段階判定を導入することです。

1. 軽量な事前選別
2. 重量な本判定

すなわち、

$$
\text{ON row} \rightarrow \text{事前選別} \rightarrow \text{候補 row のみ本判定}
$$

という流れに変更します。

このとき、

- 既存の HOT/OFF 判定ロジックは変えない
- ON の最終判定基準は変えない
- 追加パラメータはできるだけ増やさない

ことを重視します。

---

## 3. 事前選別で見る量

### 3.1 定義

各 ON row のスペクトルを $S(i)$ とし、隣接チャネル差分を

$$
D(i) = S(i+1) - S(i)
$$

と定義します。

さらに、この差分系列の robust な尺度を

$$
\sigma_D = 1.4826 \, \mathrm{MAD}(D)
$$

と定義します。

その上で、差分の最大規格化振幅を

$$
Q_D = \frac{\max_i |D(i)|}{\sigma_D}
$$

と定義します。

### 3.2 直感的意味

狭帯域 spur は、1 ch から数 ch の間で急激に立ち上がって急激に戻ることが多いです。そのため、隣接チャネル差分 $D(i)$ に大きな値が現れやすくなります。

一方で、

- 全体 gain 変化
- 緩やかな bandpass 傾き
- 広い線構造
- ゆっくりした ripple

は隣接チャネル差分には強くは出にくいです。

したがって、$Q_D$ は「この row に鋭い狭帯域異常がありそうか」を軽く見る指標として使えます。

---

## 4. 最小例

### 4.1 1 ch の鋭い spur

例えば

$$
S = [0, 0, 0, 10, 0, 0]
$$

なら、

$$
D = [0, 0, 10, -10, 0]
$$

となります。

このとき $\max |D|$ は大きく、$Q_D$ も大きくなります。

### 4.2 2 ch の spur

$$
S = [0, 0, 8, 8, 0, 0]
$$

なら、

$$
D = [0, 8, 0, -8, 0]
$$

となり、やはり両端で大きな差分が出ます。

### 4.3 滑らかな変化

$$
S = [0, 1, 2, 3, 4, 5]
$$

なら、

$$
D = [1, 1, 1, 1, 1]
$$

となり、急激な跳びはありません。

このように、$Q_D$ は「狭い構造の鋭さ」を簡単に見る量です。

---

## 5. 事前選別の判定規則

事前選別では、各 ON row に対して $Q_D$ を計算し、

$$
Q_D > T_D
$$

を満たす row だけを本判定へ進めます。ここで $T_D$ は事前選別しきい値です。

逆に、

$$
Q_D \le T_D
$$

の row は「狭帯域 spur の可能性が低い」とみなし、本判定を省略します。

### 5.1 重要な解釈

これは

- $Q_D$ が小さいなら spur が絶対に無い

という意味ではありません。

正しくは、

- $Q_D$ が小さい row には、少なくとも「隣接チャネルで急に跳ぶ型の 1 ch から数 ch の狭帯域 spur」は出にくい

という意味です。

したがって、これは本判定そのものではなく、**高速化のための安全側の事前ふるい**です。

---

## 6. 処理フロー

提案する処理フローは次です。

### 6.1 HOT/OFF

HOT/OFF は現行通りです。

$$
\text{HOT/OFF} \rightarrow \text{狭帯域 spur 本判定}
$$

### 6.2 ON

ON だけを次の 2 段階にします。

$$
\text{ON} \rightarrow Q_D \text{ による事前選別} \rightarrow \text{候補 row のみ狭帯域 spur 本判定}
$$

### 6.3 擬似コード

```text
for each ON row:
    compute D(i) = S(i+1) - S(i)
    compute sigma_D = 1.4826 * MAD(D)
    compute Q_D = max(|D|) / sigma_D

    if Q_D <= T_D:
        mark as pre-screen pass
        skip rolling-median spur QC
    else:
        run existing narrow-spur QC
```

---

## 7. 既存本判定との関係

事前選別を通過した row だけに、現行の本判定を実行します。

本判定では、

1. rolling median $B(i)$ を作る
2. 残差 $R(i) = S(i) - B(i)$ を作る
3. $\sigma_R = 1.4826 \, \mathrm{MAD}(R)$ を作る
4. $|R(i)| > q \, \sigma_R$ を見る
5. 連続超過長 $L$ が $1 \le L \le L_{\max}$ なら bad

とします。

したがって、最終的な bad 判定そのものは変えません。

本機能はあくまで、

$$
\text{全 ON row に本判定} \rightarrow \text{候補 ON row のみ本判定}
$$

へ変えるだけです。

---

## 8. パラメータ設計

### 8.1 推奨方針

新しい公開パラメータは増やしすぎない方がよいです。

そのため、本機能ではまず次の方針を推奨します。

- `qc_spur_consider_on=True` のときだけ自動的に ON 事前選別を使う
- 事前選別しきい値 $T_D$ は内部固定値にする
- HOT/OFF には事前選別を入れない

### 8.2 内部固定しきい値

初期値としては

$$
T_D = 8
$$

程度が無難です。

これはかなり保守的な値であり、明らかに鋭い跳びを持つ row だけを本判定へ送る意図です。

### 8.3 将来拡張

将来必要なら `qc_spur_on_prescreen_sigma` のような公開引数を追加することは可能ですが、最初から公開するとパラメータが増えすぎるため、初期実装では非推奨です。

---

## 9. ON に対してのみ使う理由

この事前選別は主に ON 用です。

理由は次の通りです。

- HOT/OFF は本数が少ないので、そのまま本判定しても計算時間が大きな問題になりにくい
- ON は本数が多いため、全 row に rolling median をかけると遅くなりやすい
- ON では本物の天体線があるので、判定の意味を壊さないように本判定は維持したい

したがって、

- HOT/OFF は現行通り
- ON のみ事前選別を追加

が最も自然です。

---

## 10. ON に対する注意点

ON には本物のスペクトル線が入り得ます。特に

- 強い狭線
- maser 的に非常に細い線
- 非常に高い S/N の狭い線

では、$Q_D$ が大きくなることがあります。

そのため、事前選別は「本判定へ送る候補を広めに拾う」役割に限定し、

- $Q_D$ が大きいから bad

とはしません。

bad 判定はあくまで rolling median ベースの本判定で行います。

これは重要です。

本機能の役割は、

$$
\text{bad 判定} \ne \text{事前選別}
$$

です。

事前選別は、

$$
\text{本判定が必要そうな row を逃しにくく選ぶ}
$$

ためのものです。

---

## 11. 計算量の考え方

### 11.1 現状

現状では、`qc_spur_consider_on=True` のとき、ON 全 row に本判定をかけます。

ON row 数を $N_{\mathrm{on}}$、channel 数を $N_{\mathrm{ch}}$、rolling median 窓幅を $W$ とすると、重い部分は概ね

$$
N_{\mathrm{on}} \times N_{\mathrm{ch}} \times W
$$

に比例して増えます。

### 11.2 事前選別導入後

ON row のうち、事前選別を通る割合を $f$ とすると、重い rolling median 本判定は概ね

$$
f \, N_{\mathrm{on}} \times N_{\mathrm{ch}} \times W
$$

になります。

通常は $f \ll 1$ を期待するので、かなりの高速化が見込めます。

---

## 12. 出力と記録

本機能を実装する場合、記録は次のいずれかが考えられます。

### 12.1 最小構成

内部処理だけに使い、table には残さない。

長所:

- 変更が最小
- table 列が増えない

短所:

- どれだけ pre-screen で落ちたか見えない

### 12.2 推奨構成

ON に対して次の列を追加する。

- `REF_QC_ON_PRESCREEN_QD`
- `REF_QC_ON_PRESCREEN_PASS`

意味は

- `REF_QC_ON_PRESCREEN_QD`: 各 ON row の $Q_D$
- `REF_QC_ON_PRESCREEN_PASS`: 本判定を省略したかどうか

です。

ただし、まずは列を増やさず内部だけで使う実装でもよいです。

---

## 13. 推奨設定

この機能の推奨設定は、厳密には観測対象、分光分解能、実際に出る spur の幅、ON に本物の細い線がどれくらいあるかに依存します。ここでは、最初に試すべき設定を目的別に整理します。

### 13.1 基本方針

まず一番大事なのは、ON 事前選別のしきい値 $T_D$ は、

- spur を直接判定する閾値ではない
- 本判定へ送る候補 row を絞るための閾値である

という点です。

したがって、最初は

- 候補を絞りすぎない
- 本判定へ送る row をやや多めに残す

という安全側の運用がよいです。

### 13.2 最初の推奨値

最初は内部固定値として

$$
T_D = 8
$$

を推奨します。

この値は、

- 1 ch から 3 ch の鋭い spur は拾いやすい
- ゆるやかな bandpass 変化では反応しにくい
- ON に本物の広い線がある場合でも過剰には反応しにくい

という意味で、かなり無難です。

### 13.3 安全側の設定

ON に本物の狭い線が入る可能性があり、できるだけ pre-screen による候補増加を抑えたい場合は、

$$
T_D = 10
$$

程度の高め設定が無難です。

この場合、

- 非常に鋭い跳びを持つ row だけ本判定へ回す
- 本判定対象 row 数は少なくなりやすい
- 速度は出やすい
- ただし弱い 2 ch から 3 ch spur を見逃す可能性はやや増える

という傾向になります。

### 13.4 標準設定

一般的には

$$
T_D = 8
$$

を標準設定とするのがよいです。

この場合、

- 1 ch spur に十分反応しやすい
- 2 ch から 3 ch spur にも反応しやすい
- 過剰に候補 row を増やしすぎにくい

というバランスになります。

### 13.5 攻める設定

spur を取りこぼしたくなく、本判定 row 数が増えてもよい場合は、

$$
T_D = 6 \text{ から } 7
$$

を試す余地があります。

正しくは、

$$
T_D = 6 \text{ から } 7
$$

です。

この場合、

- 弱い狭帯域異常でも本判定へ送られやすい
- 速度向上効果は小さくなる
- 本物の細い線を持つ ON row も候補に入りやすい

ため、line-rich な ON には慎重です。

### 13.6 かなり保守的な設定

ON に鋭い astrophysical line が多く、pre-screen はあくまで非常に明瞭な spur だけに限定したいなら、

$$
T_D = 12
$$

程度も候補です。

この場合、

- 明らかな 1 ch spur にはまだ反応する可能性が高い
- 中程度の 2 ch から 3 ch spur は本判定に送られないことが増える
- 高速化は大きくなりやすい

という性質になります。

### 13.7 推奨設定のまとめ

最初に試すべき設定をまとめると、次のようになります。

- 安全側: $T_D = 10$
- 標準: $T_D = 8$
- 取りこぼしを減らす: $T_D = 6$ から $7$
- かなり保守的: $T_D = 12$

最初は

$$
T_D = 8
$$

で始め、

- まだ遅いなら $T_D$ を上げる
- 弱い spur が通り抜けるなら $T_D$ を下げる

の順で調整するのが分かりやすいです。

### 13.8 本判定パラメータとの組み合わせ

ON pre-screen は単独ではなく、既存の狭帯域 spur 本判定と組み合わさって使われます。したがって、`qc_spur_sigma`、`qc_spur_window_channels`、`qc_spur_max_consecutive_channels` との組み合わせも重要です。

#### 13.8.1 まずの推奨

- `qc_spur_sigma = 7.0`
- `qc_spur_window_channels = 31`
- `qc_spur_max_consecutive_channels = 3`
- ON pre-screen は $T_D = 8$

これは一番無難な出発点です。

#### 13.8.2 1 ch spur を主に見たい場合

- `qc_spur_sigma = 7.5`
- `qc_spur_window_channels = 31`
- `qc_spur_max_consecutive_channels = 1`
- ON pre-screen は $T_D = 8$ から $10$

この設定は、針のような 1 ch 異常に敏感で、少し広い構造には反応しにくいです。

#### 13.8.3 2 ch から 3 ch spur も見たい場合

- `qc_spur_sigma = 7.0`
- `qc_spur_window_channels = 31`
- `qc_spur_max_consecutive_channels = 3`
- ON pre-screen は $T_D = 8$

この組み合わせが標準です。

#### 13.8.4 弱い spur を取りこぼしたくない場合

- `qc_spur_sigma = 6.0`
- `qc_spur_window_channels = 31`
- `qc_spur_max_consecutive_channels = 3`
- ON pre-screen は $T_D = 6$ から $7$

この場合、候補 row は増えやすくなり、高速化効果は弱まります。

#### 13.8.5 ON に本物の狭線が多い場合

- `qc_spur_sigma = 7.5` あるいは $8.0$
- `qc_spur_window_channels = 31`
- `qc_spur_max_consecutive_channels = 1` あるいは $2$
- ON pre-screen は $T_D = 10$ あるいは $12$

この設定は、できるだけ astrophysical narrow line に引きずられにくくする方向です。

### 13.9 速度優先と見逃し優先

設定の考え方は次の 2 軸です。

#### 13.9.1 速度優先

- $T_D$ を上げる
- `qc_spur_sigma` もやや高めにする
- `qc_spur_max_consecutive_channels` を 1 あるいは 2 にする

この方向では本判定 row 数が減りやすく、高速です。ただし弱い spur を取りこぼしやすくなります。

#### 13.9.2 見逃し低減優先

- $T_D$ を下げる
- `qc_spur_sigma` もやや低めにする
- `qc_spur_max_consecutive_channels` を 3 にする

この方向では本判定 row 数が増えやすく、速度改善は減りますが、spur の取りこぼしは減りやすいです。

### 13.10 実務上の推奨手順

実データでは、最初から `FLAGROW` に反映するより、まず report 的に動かす方が安全です。したがって、実務上は次の順を推奨します。

1. まず標準設定で試す
2. pre-screen を通る ON row の割合を見る
3. 候補 row が多すぎるなら $T_D$ を上げる
4. 明らかな spur を通しすぎるなら $T_D$ を下げる
5. その後に `qc_apply_flagrow` を使うか検討する

最初の標準設定は次です。

- `qc_spur_sigma = 7.0`
- `qc_spur_window_channels = 31`
- `qc_spur_max_consecutive_channels = 3`
- `qc_spur_consider_on = True`
- ON pre-screen は $T_D = 8$
- 最初は `qc_apply_flagrow = False`

### 13.11 推奨設定の見方

ここでの推奨値は、論理的に一意に決まるものではなく、あくまで観測装置と実データ傾向に依存する出発点です。確定なのは次の関係です。

- $T_D$ を上げるほど速くなりやすい
- $T_D$ を下げるほど候補 row は増えやすい
- `qc_spur_sigma` を上げるほど本判定は厳しくなる
- `qc_spur_max_consecutive_channels` を小さくするほど、細い spur 専用になる

一方で、どの値が最適かは、実データで

- 候補 row 数
- 本判定 row 数
- 実際に見つかった spur の質
- astrophysical line の誤反応

を見て決める必要があります。

---

## 14. 最終 bad 判定との関係

本機能を導入しても、最終 bad 判定は次を維持します。

### 14.1 HOT/OFF

$$
\text{auto\_bad} = \text{shape\_bad} \lor \text{spur\_bad}
$$

### 14.2 ON

`qc_spur_consider_on=True` のとき、ON については

- pre-screen で候補外なら `spur_bad = False`
- pre-screen で候補なら本判定の結果を使う

とします。

したがって、pre-screen は ON の `spur_bad` 計算を省略するための入口条件です。

---

## 15. 実装方針

実装は `calibrate.py` の狭帯域 spur QC 関数とは分けて、例えば次の補助関数を追加するのが自然です。

```text
_spur_prescreen_qd_rows(...)
```

役割は

- 入力: 2-D spectra array
- 出力: `candidate_mask`, `qd_score`

です。

その後、既存の `_narrow_spur_qc_rows(...)` に `candidate_mask` を渡す構成にします。

これにより、既存の本判定関数の責務は保ったまま、ON だけ前段で間引けます。

---

## 16. 推奨初期仕様

最小仕様としては次を推奨します。

1. `qc_spur_consider_on=True` のときだけ ON pre-screen を有効化
2. HOT/OFF は現行通り
3. pre-screen 指標は $Q_D$ のみ
4. pre-screen しきい値 $T_D$ は内部固定値
5. bad 判定は既存の rolling median 本判定のみ
6. 公開パラメータは増やさない

この仕様なら、

- 高速化の効果が見込める
- ユーザー向けパラメータが増えない
- 既存の spur 判定意味を壊さない

という利点があります。

---

## 17. まとめ

本機能は、`qc_spur_consider_on=True` のときだけ ON row に導入する高速化機能です。

核心は、隣接チャネル差分

$$
D(i) = S(i+1) - S(i)
$$

から

$$
Q_D = \frac{\max_i |D(i)|}{1.4826 \, \mathrm{MAD}(D)}
$$

を作り、$Q_D$ が十分大きい row にだけ重い rolling median ベースの本判定を行うことです。

これにより、

- 全 ON row に本判定をかける現在の遅さを改善しつつ
- 最終的な狭帯域 spur 判定の意味は保ちやすい

という設計になります。

本機能は、狭帯域 spur の本判定を置き換えるものではなく、**その前段に入る軽量な候補選別機能**です。


---

## 18. 現行実装と提案仕様の整理

この文書では、内容を次の 2 つに分けて扱う。

### 18.1 現行実装として入っているもの

本チャットでの最新版 `calibrate_updated_2026-04-10_badchan.py` には、少なくとも次が入っている。

- HOT と OFF に対する broad shape QC
- HOT と OFF に対する狭帯域 spur QC
- 必要なら ON に対する狭帯域 spur QC
- `bad_channels` による全体共通の fixed bad channel 指定
- `bad_channel_ranges` による全体共通の fixed bad channel range 指定
- `bad_channel_map` による group ごとの fixed bad channel 指定
- `bad_channel_policy="nan"`
- `bad_channel_policy="interp_linear"`
- QC 計算時に fixed bad channel を除外する処理
- 較正計算時に fixed bad channel を `NaN` または補間値で扱う処理

### 18.2 現時点では提案仕様のもの

次は、現時点では提案仕様として扱う。

- `qc_spur_consider_on=True` のときの ON 事前選別機能
- その事前選別の指標としての $Q_D$
- ON row の一部だけに rolling median ベース本判定を行う高速化

したがって、読み分けは次のようにする。

- Section 1 から Section 17 までは、主に ON 事前選別の提案仕様
- Section 18 以降は、現行実装に入っている fixed bad channel 機能と、その運用指針の追記

---

## 19. fixed bad channel 問題の整理

### 19.1 背景

常に同じ raw channel で大きな spur が出る場合、それは row 全体の異常ではなく、特定 channel の恒常的不良である可能性が高い。

raw channel index を $k$、時刻または dump index を $t$、受信系ごとの group を $g$ とする。

各 state の raw spectrum を

$$
S_{\mathrm{ON}}(k,t), \quad S_{\mathrm{OFF}}(k,t), \quad S_{\mathrm{HOT}}(k,t)
$$

とする。

fixed bad channel 集合を

$$
K_{\mathrm{bad}}(g)
$$

とする。

このとき、問題は

$$
k \in K_{\mathrm{bad}}(g)
$$

である channel が常に不良、ということである。

### 19.2 row QC と channel 処理の違い

row QC は、「その row 全体が使えないかどうか」を決める処理である。

一方、fixed bad channel は「その row の一部 channel だけが不良」である。

したがって、fixed bad channel に対して row 全体を `FLAGROW` で落とし続けるのは損失が大きい。

このため、fixed bad channel は row QC とは分離し、channel-level 処理として扱う。

---

## 20. `bad_channels` / `bad_channel_ranges` / `bad_channel_map`

### 20.1 全体共通指定

全 group 共通で bad とみなす raw channel は、次で指定する。

```python
bad_channels: Optional[Sequence[int]] = None
bad_channel_ranges: Optional[Sequence[Tuple[int, int]]] = None
```

例:

```python
bad_channels = [1536, 3072]
bad_channel_ranges = [(2048, 2050)]
```

ここで `bad_channel_ranges` の `(start, stop)` は、両端を含む inclusive range とする。

### 20.2 group 別指定

group ごとに異なる fixed bad channel を指定したい場合は、次を使う。

```python
bad_channel_key_columns = ("FDNUM", "IFNUM", "PLNUM")
bad_channel_map = {
    (2, 0, 0): [1536],
    (3, 0, 1): [1536, (2048, 2050)],
}
```

ここで key tuple の意味は `bad_channel_key_columns` で決まる。

既定では

$$
g = (\mathrm{FDNUM}, \mathrm{IFNUM}, \mathrm{PLNUM})
$$

である。

### 20.3 global と group の合成

ある group $g$ に対する最終 bad channel 集合は、global 設定と map 設定の和集合とする。

$$
K_{\mathrm{bad}}(g) = K_{\mathrm{global}} \cup K_{\mathrm{map}}(g)
$$

つまり `bad_channel_map` は global 指定を打ち消すものではなく、加算的に追加するものである。

---

## 21. channel index の定義

`bad_channels`、`bad_channel_ranges`、`bad_channel_map` は、すべて **元の raw spectrum の channel index** で指定する。

これは重要である。

- `ch_range` で切り出した後の local index ではない
- `vlsrk_range_kms` で切り出した後の local index ではない
- backend native の raw channel index を使う

内部では、選択後の local channel index へ変換する。

切り出し開始 channel を $k_0$ とすると、local index は

$$
k_{\mathrm{local}} = k_{\mathrm{raw}} - k_0
$$

である。

もし指定した raw channel が現在の切り出し範囲外なら、その channel は単に無視する。

---

## 22. 処理フローにおける位置づけ

現行実装では、まず `ch_range` や `vlsrk_range_kms` に従って `on2_arr`、`off2_arr`、`hot2_arr` を作り、その後に fixed bad channel 処理を適用する。

流れは次の通りである。

1. `on2_arr`、`off2_arr`、`hot2_arr` を作る
2. group ごとに local bad channel 集合を解決する
3. QC 用配列では fixed bad channel を除外する
4. 較正用配列では `bad_channel_policy` に従って `NaN` 化または補間する
5. その後に broad shape QC、狭帯域 spur QC、平均、時間内挿、gain 計算、較正へ進む

したがって、fixed bad channel は

- QC を誤作動させない
- かつ較正本体も汚さない

という 2 つの役割を同時に持つ。

---

## 23. QC 計算での扱い

fixed bad channel は、spur row 判定や shape QC の対象ではない。

したがって、QC の統計量は

$$
k \notin K_{\mathrm{bad}}(g)
$$

の channel だけで計算する。

たとえば residual を

$$
R(k) = S(k) - B(k)
$$

と定義し、robust scale を

$$
\sigma_R = 1.4826 \, \mathrm{MAD}(R)
$$

と定義するときも、`MAD` や最大値の計算には fixed bad channel を含めない。

これにより、同じ channel の恒常 spur が毎回 HOT/OFF row を自動 flag してしまうことを防ぐ。

---

## 24. `bad_channel_policy="nan"`

### 24.1 定義

既定の policy は `"nan"` とする。

このとき、較正用のスペクトルは

$$
S'(k,t) =
\begin{cases}
\mathrm{NaN} & k \in K_{\mathrm{bad}}(g) \\
S(k,t) & \text{otherwise}
\end{cases}
$$

として扱う。

### 24.2 目的

- fixed bad channel を gain 計算から除外する
- 情報を勝手に作らない
- 後段で bad channel が明示的に分かるようにする

### 24.3 特徴

- 最も安全
- 科学解析の標準として自然
- ただし見た目は不連続になりやすい

---

## 25. `bad_channel_policy="interp_linear"`

### 25.1 目的

`interp_linear` は、fixed bad channel を隣接 channel の値から 1 次元線形補間で埋める mode である。

この mode の狙いは、

- 見た目の連続性を改善する
- 1 ch から数 ch の固定 spur による不自然な穴を減らす
- ただし interpolation は最小限に抑える

ことである。

### 25.2 補間の基本式

連続 bad 区間が

$$
[k_1, k_2]
$$

であり、両端の有効 channel が $k_L$ と $k_R$ で、それぞれの値が $S(k_L)$、$S(k_R)$ であるとする。

このとき、$k_1 \le k \le k_2$ に対して

$$
S'(k) = S(k_L) + \frac{k - k_L}{k_R - k_L} \left( S(k_R) - S(k_L) \right)
$$

で埋める。

1 ch の bad channel $k_0$ なら、これは

$$
S'(k_0) = \frac{S(k_0-1) + S(k_0+1)}{2}
$$

に一致する。

### 25.3 補間の条件

補間は次の条件を満たすときだけ行う。

- 連続 bad channel 数が `bad_channel_interp_max_consecutive_channels` 以下
- 左右に少なくとも 1 点ずつ有効 channel がある
- 区間が spectrum の edge に接していない

それ以外は補間せず、その区間は `NaN` のままとする。

### 25.4 なぜ edge を埋めないか

edge では片側の情報しかないため、線形補間ではなく外挿になる。外挿はより不確実で、勝手に構造を作りやすい。

したがって、最小仕様では edge は埋めない。

### 25.5 適用上の注意

`interp_linear` は便利だが、補間値は観測値ではない。したがって、狭い本物の線が bad channel 上にあった場合、その構造は失われる。

このため、科学解析の標準は依然として `"nan"` が安全である。

一方、

- quicklook
- 図示
- baseline の見た目の安定化
- 1 ch 固定 spur が明らかな場合

には `interp_linear` が有用である。

---

## 26. `FLAGROW` との役割分担

fixed bad channel は channel-level の問題であり、それ自体では row 全体を bad にする理由ではない。

したがって、固定 bad channel 自体では `FLAGROW` を立てない。

役割分担は次の通りである。

- fixed bad channel は `bad_channels` / `bad_channel_map` で扱う
- row 全体の異常は broad shape QC や狭帯域 spur QC で扱う
- `FLAGROW` は後者に対して使う

この分離が重要である。

---

## 27. 出力列と記録

最新実装では、fixed bad channel 関連の情報として、少なくとも次の列を出力できる。

- `BADCHAN_POLICY`
- `BADCHAN_N_GLOBAL`
- `BADCHAN_N_GROUP`
- `BADCHAN_N_TOTAL`
- `BADCHAN_INTERP_APPLIED`
- `BADCHAN_INTERP_NFILLED`
- `BADCHAN_HAS_EDGE_UNFILLED`

意味は次の通りである。

- `BADCHAN_POLICY`
  - `nan` か `interp_linear`
- `BADCHAN_N_GLOBAL`
  - global 設定に由来する bad channel 数
- `BADCHAN_N_GROUP`
  - map 設定に由来する bad channel 数
- `BADCHAN_N_TOTAL`
  - 実際にその row へ適用された総 bad channel 数
- `BADCHAN_INTERP_APPLIED`
  - 補間が 1 箇所以上実際に行われたか
- `BADCHAN_INTERP_NFILLED`
  - 補間で埋めた channel 数
- `BADCHAN_HAS_EDGE_UNFILLED`
  - edge に接していて埋められなかった bad channel があるか

また、狭帯域 spur QC 側では引き続き

- `QC_REF_HOT_AUTO_BAD`
- `QC_REF_OFF_AUTO_BAD`
- `QC_REF_HOT_SHAPE_AUTO_BAD`
- `QC_REF_OFF_SHAPE_AUTO_BAD`
- `QC_REF_HOT_SPUR_AUTO_BAD`
- `QC_REF_OFF_SPUR_AUTO_BAD`
- `QC_ON_SPUR_AUTO_BAD`
- `QC_SIGMA`
- `QC_SPUR_SIGMA`
- `QC_SPUR_WINDOW_CH`
- `QC_SPUR_MAXRUN`
- `QC_SPUR_CONSIDER_ON`
- `QC_APPLY_FLAGROW`

などが記録される。

したがって、fixed bad channel と row QC の両方を出力 table 上で追跡できる。

---

## 28. 銀河中心のように広い emission がある場合の考え方

### 28.1 ON に対する注意

銀河中心や line-rich source では、ON に本物の広い emission や複雑な line structure がある。このため、`qc_spur_consider_on=True` のときは、ON の本物の線構造が spur 候補として本判定に送られやすくなる。

したがって、line-rich source では原則として

```python
qc_spur_consider_on = False
```

を基本にする方が安全である。

### 28.2 HOT が大量に flag される場合

HOT が大量に flag される場合、原因は 2 つに分かれる。

1. HOT row 自体に本当に row-level の異常が多い
2. 同じ raw channel の固定 spur が毎回 QC を発火させている

後者なら、問題の本質は row 異常ではなく fixed bad channel である。

この場合は `qc_spur_sigma` の調整で粘るより、まず `bad_channels` あるいは `bad_channel_map` でその channel を明示的に指定すべきである。

---

## 29. 推奨設定

ここでは、現行実装の狭帯域 spur QC と fixed bad channel 機能を、実データでどのように使い始めるかを、なるべく豊富に整理する。

### 29.1 まずの基本方針

最初に大事なのは次の 3 点である。

1. `qc_spur_consider_on` は最初は `False`
2. fixed bad channel が分かっているなら先に `bad_channels` へ入れる
3. `qc_apply_flagrow` は最初は `False`

つまり、最初は「件数と挙動を見るだけ」にする。

### 29.2 最小の標準設定

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    qc_apply_flagrow=False,
)
```

この設定は、

- HOT/OFF の broad shape QC
- HOT/OFF の狭帯域 spur QC
- ON は見ない
- まず report 的に挙動を見る

という意味で、最初の出発点として無難である。

### 29.3 fixed 1 ch spur が既知の場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    bad_channels=[1536],
    bad_channel_policy="nan",
    qc_spur_consider_on=False,
    qc_apply_flagrow=False,
)
```

この設定では、raw channel 1536 を全 group 共通で bad channel として扱う。

科学解析の標準としては、まず `bad_channel_policy="nan"` を推奨する。

### 29.4 fixed 1 ch spur を補間したい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    bad_channels=[1536],
    bad_channel_policy="interp_linear",
    bad_channel_interp_max_consecutive_channels=1,
    qc_spur_consider_on=False,
    qc_apply_flagrow=False,
)
```

1 ch 固定 spur なら、この設定はかなり使いやすい。

### 29.5 2 ch から 3 ch の固定 spur が既知の場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    bad_channel_ranges=[(1536, 1538)],
    bad_channel_policy="interp_linear",
    bad_channel_interp_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    qc_apply_flagrow=False,
)
```

この場合、連続長 3 ch までなら線形補間で埋める。

### 29.6 group ごとに異なる固定 spur がある場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    bad_channel_key_columns=("FDNUM", "IFNUM", "PLNUM"),
    bad_channel_map={
        (2, 0, 0): [1536],
        (3, 0, 1): [1536, (2048, 2050)],
    },
    bad_channel_policy="nan",
    qc_spur_consider_on=False,
    qc_apply_flagrow=False,
)
```

複数ビーム、複数 IF、複数偏波で channel 不良位置が違う場合は、この形が最も自然である。

### 29.7 global と group を併用する場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    bad_channels=[1536],
    bad_channel_map={
        (3, 0, 1): [(2048, 2050)],
    },
    bad_channel_policy="nan",
    qc_spur_consider_on=False,
    qc_apply_flagrow=False,
)
```

この場合、

- 全 group で ch 1536 を除外
- さらに特定 group だけ 2048 から 2050 を除外

となる。

### 29.8 銀河中心のように line-rich な source

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    bad_channel_map={
        (2, 0, 0): [1536],
    },
    bad_channel_policy="nan",
    qc_apply_flagrow=True,
)
```

この場合の考え方は次の通りである。

- ON は本物の line structure が多いので見ない
- HOT/OFF の row QC は使う
- fixed bad channel は map で除外する
- row 全体を落とすのは本当に row-level に異常がある場合だけにする

### 29.9 できるだけ安全側に始める場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.5,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=2,
    qc_spur_consider_on=False,
    bad_channel_policy="nan",
    qc_apply_flagrow=False,
)
```

特徴は、

- 狭帯域 spur 判定をやや保守的にする
- 1 ch から 2 ch の鋭いもの寄りにする
- 最初は何も自動で落とさない

である。

### 29.10 1 ch spur を重点的に見たい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.5,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=1,
    qc_spur_consider_on=False,
    bad_channel_policy="nan",
    qc_apply_flagrow=False,
)
```

### 29.11 2 ch から 3 ch spur も積極的に見たい場合

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=6.5,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=False,
    bad_channel_policy="nan",
    qc_apply_flagrow=False,
)
```

この設定では、弱めの狭帯域 spur も拾いやすくなるが、候補は増えやすい。

### 29.12 ON も見たいが最初は report only にしたい場合

現行実装には `report only` 専用スイッチはないため、最初は `qc_apply_flagrow=False` で見るのが安全である。

```python
sc_cal = cal.run_tastar_calibration(
    input_data=sc_raw,
    qc_sigma=5.0,
    qc_spur_sigma=7.0,
    qc_spur_window_channels=31,
    qc_spur_max_consecutive_channels=3,
    qc_spur_consider_on=True,
    bad_channel_policy="nan",
    qc_apply_flagrow=False,
)
```

ただし line-rich source では、やはり `qc_spur_consider_on=False` を推奨する。

### 29.13 ON 事前選別を将来使う場合の推奨

Section 1 から Section 17 の提案仕様が将来実装された場合、最初の推奨は次である。

- ON pre-screen の内部しきい値は $T_D = 8$
- `qc_spur_sigma = 7.0`
- `qc_spur_window_channels = 31`
- `qc_spur_max_consecutive_channels = 3`
- 最初は `qc_apply_flagrow = False`

ただし、これも line-rich source では注意が必要である。

### 29.14 `nan` と `interp_linear` の使い分け

#### 科学解析の標準

```python
bad_channel_policy = "nan"
```

を推奨する。

#### quicklook や見た目重視

```python
bad_channel_policy = "interp_linear"
bad_channel_interp_max_consecutive_channels = 1
```

から始めるのがよい。

#### 連続 2 ch から 3 ch の固定不良も埋めたい

```python
bad_channel_policy = "interp_linear"
bad_channel_interp_max_consecutive_channels = 3
```

とする。

ただし、長い bad 区間や edge に接する bad 区間は、原理的に不確かさが大きいので、無理に埋めない方が安全である。

---

## 30. 実務上の推奨手順

現実の運用では、次の順で進めるのが分かりやすい。

### Step 1

まず fixed bad channel を known issue として入れる。

- 全 group 共通なら `bad_channels`
- group 別なら `bad_channel_map`

### Step 2

最初は

```python
qc_apply_flagrow = False
```

で走らせる。

### Step 3

出力 table の

- `BADCHAN_*`
- `QC_REF_HOT_*`
- `QC_REF_OFF_*`
- 必要なら `QC_ON_SPUR_*`

を確認する。

### Step 4

fixed bad channel を入れてもなお HOT/OFF が大量に bad なら、その時点で初めて `qc_spur_sigma` や `qc_sigma` を調整する。

### Step 5

挙動に納得してから

```python
qc_apply_flagrow = True
```

へ進む。

この順にすると、

- fixed bad channel と row 異常を混同しにくい
- 何が効いたのかが見えやすい
- いきなり良い row を大量に捨てる事故を減らせる

---

## 31. まとめ

この文書の更新後の要点は次の通りである。

1. ON 事前選別は依然として提案仕様である
2. 現行実装には、fixed bad channel を扱う `bad_channels`、`bad_channel_ranges`、`bad_channel_map` が入っている
3. fixed bad channel は row 異常ではなく channel-level 異常である
4. したがって、fixed bad channel 自体では `FLAGROW` を立てない
5. QC 計算では fixed bad channel を除外する
6. 較正計算では `bad_channel_policy` に従って `NaN` または補間で扱う
7. 科学解析の標準は `bad_channel_policy="nan"`
8. quicklook や見た目重視では `bad_channel_policy="interp_linear"` が有用である
9. line-rich source では原則 `qc_spur_consider_on=False` が安全である
10. fixed bad channel が既知なら、まず `bad_channels` または `bad_channel_map` を入れてから row QC を調整するべきである
