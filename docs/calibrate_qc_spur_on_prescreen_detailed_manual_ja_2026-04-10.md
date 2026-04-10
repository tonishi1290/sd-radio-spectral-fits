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
T_D = 6
text{ から } 7
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
