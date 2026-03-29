# `map_3d` OTF gridding kernel 詳細説明書

Date: 2026-03-29  
Revision: r5  
Target package: `map_3d`  
Scope: OTF gridding における空間 kernel の理論・実装・既定値・使い分け・未実装事項・注意点

---

## 0. この文書の目的

この文書は、`map_3d` パッケージにおける **OTF gridding kernel** の理論と実装を、利用者・開発者の双方が参照できる形で詳細に整理した説明書である。

本書では、次を明確に区別して記述する。

1. **理論**: OTF gridding で kernel が何を決めるのか。  
2. **外部参照に基づく仕様背景**: CASA / casacore / Mangum et al. などが何を採っているか。  
3. **本パッケージの現在実装**: 何が実装済みで、どのように解釈されるか。  
4. **推奨運用**: どの kernel をどの状況で使うべきか。  
5. **未実装・未検証**: 何をまだやっていないか。  
6. **誤解しやすい点**: support, cell, beam, sign, effective beam などの取り違えを防ぐ。

この文書を読めば、少なくとも次が分かるようにする。

- OTF gridding の kernel が、最終マップの **空間分解能**, **ノイズ**, **aliasing**, **欠損耐性** にどう効くか  
- `sf`, `gjinc`, `gauss` が、このパッケージでどう定義されているか  
- `mangum`, `casa`, `legacy`, `mangum2007` がどういう関係にあるか  
- `convsupport`, `truncate`, `gwidth`, `jwidth`, `support_radius_*` が何を意味するか  
- なぜ current default を `sf + convsupport=3` にしたのか  
- どこまで CASA / casacore を参考にし、どこがまだ未完なのか

---

## 1. OTF gridding で kernel は何を決めるのか

単一鏡の OTF (On-The-Fly) 観測では、走査中に得られる離散 sample (dump) を、最終的な 2D map / 3D cube の格子点へ再配置する必要がある。

基本形は

$$
T_{\rm out}(\mathbf{x}) =
\frac{\sum_i w_i(\mathbf{x}) T_i}{\sum_i w_i(\mathbf{x})}
$$

であり、各 sample $T_i$ は、格子点 $\mathbf{x}$ からの距離に応じた kernel 重み $w_i$ を通して寄与する。

したがって kernel の選択は、少なくとも次を決める。

### 1.1 空間分解能

出力マップの実効 PSF は、単純には

$$
\mathrm{PSF}_{\rm out} \approx \mathrm{Beam}_{\rm in} * K
$$

で与えられる。ここで

- $\mathrm{Beam}_{\rm in}$: 望遠鏡自身の入力ビーム  
- $K$: gridding kernel

である。したがって kernel の幅と形は、入力ビームをどれだけ太らせるかを直接決める。

### 1.2 ノイズの空間相関と RMS

- 広い positive kernel は周辺 sample を平均しやすく、ノイズを平滑化しやすい。  
- signed kernel は高域を残しやすく sharp だが、負重みを含むため RMS を増幅しやすい。  
- kernel の tail が長いと、ノイズ相関長も伸びやすい。

### 1.3 aliasing・gap・sampling の粗さへの頑健性

- 支持半径 (support) が狭すぎると、sampling の粗さや欠損に敏感になる。  
- 支持半径が広い positive kernel は安全だが、過度に broad になる可能性がある。  
- kernel の形は、単なる smoothing ではなく、欠損データの埋まり方や局所的なムラにも効く。

### 1.4 stripe との関係

走査線起因の stripe を主に除去するのは basket-weave / plait 側の仕事である。ただし、kernel は

- gap をどう埋めるか  
- scan edge をどれだけ滑らかにするか  
- 局所的なサンプリング密度の違いにどう反応するか

に影響するため、stripe の見え方にも間接的に効く。

---

## 2. 用語・変数・単位の定義

### 2.1 基本量

- $B$: 入力望遠鏡ビーム FWHM [arcsec]  
- $d$: 出力 image の cell size [arcsec/pixel]  
- $r$: dump 位置と格子点の距離 [arcsec]  
- $u = r/d$: cell 基準の無次元半径  
- $q = r/B$: beam 基準の無次元半径

### 2.2 support / cutoff / truncate / convsupport

この文書では、次を明確に区別する。

- **support radius**: kernel を実際に計算に使う最大距離、すなわち kernel を 0 とみなす cutoff 半径。一般概念である。  
- **truncate**: `gauss` / `gjinc` に対して cutoff 半径を指定する public パラメータ名である。  
- **convsupport**: `sf` に対する pixel 単位 cutoff 半径である。  
- **first null**: 主に `gjinc` で使う語で、kernel が最初に 0 になる半径を指す。`GJINC(casa)` では natural な cutoff 候補である。

直感的には、ある dump のまわりに半径 `support radius` の円を描き、その円の中の output pixel にだけ重みを配ると考えるとよい。すなわち

$$
r \le R_{\rm sup}
$$

の範囲だけ kernel を評価し、

$$
r > R_{\rm sup}
$$

ではその dump はその pixel に寄与しない。

本パッケージでは、`sf` の cutoff は通常 `convsupport` で与える。`gjinc` / `gauss` は `truncate` や `support_radius_*` で与える。したがって、**support は一般名、`truncate` と `convsupport` は kernel ごとの public な呼び名**と理解すると混乱しにくい。

### 2.3 入力ビームと実効ビーム

- **入力ビーム**: `beam_fwhm_arcsec` に与える telescope beam FWHM  
- **実効ビーム**: 入力ビームと kernel の畳み込みで決まる出力 PSF

### 2.4 `kernel_sign`

- `signed`: kernel の負値も含めて使う  
- `positive_only`: kernel の正値のみ使う。負値は 0 とみなす  
- `auto`: preset / kernel に応じて自動解決

### 2.5 pixel-based と beam-aware

- **pixel-based**: kernel の幅・support を cell size $d$ 基準で定義する  
- **beam-aware**: kernel の幅・support を beam size $B$ 基準で定義する

**現在の public preset はすべて pixel-based** である。beam-aware 指定子 (`*_beam`, `support_radius_beam`) は後方互換のため残しているが、public default の中心ではない。

---

## 3. このパッケージの基本設計方針

### 3.1 公開 preset はすべて pixel-based

現在の `map_3d` では、公開 preset はすべて pixel-based に統一している。すなわち、kernel の幅や support はまず `cell_arcsec` 基準で決まる。

### 3.2 `beam_fwhm_arcsec` と `cell_arcsec` の役割分担

- `beam_fwhm_arcsec`: 観測系の入力ビームを表す  
- `cell_arcsec`: 画像格子と kernel の幅・support の直接の基準量

両者の役割は異なる。`beam_fwhm_arcsec` は主に

1. 入力ビームの大きさを与える  
2. `cell_arcsec=None` のとき `beam/3` を決める  
3. 実効ビーム推定で入力 PSF として使う

ために必要である。

### 3.3 `cell_arcsec=None` の既定動作

`cell_arcsec=None` のとき、本パッケージは

$$
d = B/3
$$

を採用する。これは CASA `tsdimaging` docs および single-dish imaging notebook にある「primary beam FWHM の 1/3 を cell とする」既定方針と整合している。

現在の clean 実装では、`MapConfig` は **通常の dataclass** であり、初期化は dataclass 自動生成の `__init__` と `__post_init__` で行う。曖昧な位置引数 shorthand は public API に採用していない。

### 3.4 既定 kernel

現在の formal default は

- `kernel='sf'`  
- `convsupport=3`

である。

これは CASA `tsdimaging` docs の formal default (`SF` の `convsupport=3`) に合わせたものである。一方で、ALMA の single-dish imaging notebook では、ALMA data に対して `SF` と `convsupport=6` が推奨例として示されている。したがって、**`3` は formal default、`6` は推奨運用例の一つ**として区別して扱う。

---

## 4. 外部仕様・理論背景

この節は「このパッケージで何を参考にしたか」を、議論経緯ではなく、最終的な整理として書く。

### 4.1 Mangum, Emerson & Greisen (2007)

Mangum et al. (2007) は単一鏡 OTF imaging の理論と実務を整理した主要文献であり、sampling, image formation, gridding kernel の考え方、GJINC / GAUSS 等の式と既定値の基準を与える。

本パッケージでは、`gjinc` の pixel-based 定数

- `jwidth = 1.55 pixel`  
- `gwidth = 2.52 * sqrt(log 2) pixel`

を CASA docs 経由で採用している。CASA docs 自身も、これらの既定値を Mangum et al. 2007 由来と説明している。

### 4.2 CASA `tsdimaging`

CASA `tsdimaging` docs から、本パッケージで直接重要な点は次の通りである。

1. spatial `cell` の default は **primary beam FWHM の 1/3** である。  
2. `SF` の `convsupport` は **pixel 単位の cutoff radius** で、formal default は **3 pixel** である。  
3. `GJINC` の default は  
   - `jwidth = 1.55 pixel`  
   - `gwidth = 2.52 * sqrt(log 2) pixel`  
   - `truncate = first null`  
   である。  
4. `GAUSS` は `gwidth = sqrt(log 2) pixel`, `truncate = 3*HWHM` を既定とする。  
5. CASA task 層では `SF`, `GAUSS`, `GJINC`, `BOX`, `PB` が public kernel 候補として存在する。

本パッケージは CASA 完全複製ではないが、`sf` と `gjinc(casa)` の public parameterization は可能な限りここに合わせている。

### 4.3 CASA single-dish imaging notebook

single-dish imaging notebook には、次の実務情報がある。

1. OTF imaging では、有効分解能は telescope beam と sky sampling function の畳み込みで決まり、理論ビームよりやや悪化する。  
2. sky sampling function は典型的に primary beam の 1/3〜1/5 程度であり、そのため image pixel size として primary beam の 1/3 は安全な oversampling とみなせる。  
3. ALMA data に対しては `SF` kernel が推奨で、`convsupport=6` が推奨例として挙げられている。

本パッケージでは、最後の点を **背景知識としては採るが、そのまま default にはしない**。理由は後述の内部比較で、単一鏡 OTF では `SF(6)` の broadening が大きかったためである。

### 4.4 casacore `Sph_Conv` / `SPHFN`

`sf` 実装の直接の外部参照は `casacore::Sph_Conv` と `MathFunc2.cc` の `SPHFN` である。

#### 4.4.1 `Sph_Conv` の役割

casacore の class reference は `Sph_Conv` を spheroidal function と説明し、「Fred Schwab function を呼ぶ」実装であることを明記している。

#### 4.4.2 family selection

`MathFunc.tcc` にある `Sph_Conv::value(j)` は、cutoff `cut` とパラメータ `sphparm` から

$$
isupp = \mathrm{clip}(\lfloor 2\,cut \rfloor, 4, 8)
$$
$$
ialpha = \mathrm{clip}(\lfloor 2\,sphparm + 1 \rfloor, 1, 5)
$$
$$
j_{\max} = \min(\lfloor cut \rfloor, 7)
$$

を作り、

$$
K(j) = \mathrm{sphfn}(ialpha, isupp, j/j_{\max})
$$

を返す。

ここが重要である。`SF` は単なる「continuous PSWF を `r/convsupport` で評価したもの」ではなく、**cutoff から離散 family を選ぶ実装**になっている。

#### 4.4.3 `SPHFN`

`MathFunc2.cc` のコメントでは、`SPHFN` は Fred Schwab の Fortran `SPHFN.f` 由来であり、support width $m = 4,5,6,7,8$ と weighting exponent $\alpha = 0, 1/2, 1, 3/2, 2$ に対する rational approximation を持つ。

また、コメントには

- `IFLAG <= 0`: gridding に適した branch  
- `IFLAG > 0`: u-v plane convolution に適した Fourier-transform branch

とある。したがって spheroidal approximation は radio astronomy において、**gridding と u-v plane convolution の両方に関係する**が、現在本パッケージが使っているのは **gridding branch** である。

### 4.5 `SF` の背景をどう理解すべきか

以上を踏まえると、`SF` について背景として次のように書ける。

- `SF` は radio astronomy で長く使われてきた Fred Schwab の spheroidal approximation に基づく kernel である。  
- casacore / CASA では、single-dish gridding にも `SF` を用い、`convsupport` で cutoff 半径を pixel 単位に指定する。  
- 同じ source 系は u-v plane convolution に対応する branch も持つ。  

ただし、この文書では **単一鏡 OTF gridding における `SF`** に話を限定する。

---

## 5. 本パッケージで現在実装している kernel 一覧

### 5.1 public kernel 名

現在の public `kernel` は次の 3 種類である。

1. `sf`  
2. `gjinc`  
3. `gauss`

### 5.2 `gjinc` の preset

`kernel='gjinc'` のとき、`kernel_preset` として次を受け付ける。

- `mangum`  
- `casa`

互換 alias:

- `mangum2007 -> mangum`  
- `legacy -> casa`

### 5.3 現在の default

- `kernel='sf'`  
- `convsupport=3`  
- `cell_arcsec=None` のとき `beam_fwhm_arcsec / 3`

### 5.4 `kernel_sign='auto'` の解決

- `sf` -> `positive_only`  
- `gauss` -> `positive_only`  
- `gjinc + mangum` -> `signed`  
- `gjinc + casa` -> `positive_only`

### 5.5 support の既定値

- `sf`: `support = convsupport * cell`  
- `gjinc + mangum`: `support = 3 * cell`  
- `gjinc + casa`: `support = first null of resolved jinc width`  
- `gauss`: `support = 3 * gwidth`

---

## 6. public API と内部解決の対応

### 6.1 `MapConfig` / `GridConfig` の kernel 関連 public パラメータ

現在の実装では `MapConfig` が本体であり、`GridConfig = MapConfig` は互換 alias である。kernel 関連の public パラメータは次である。

- `beam_fwhm_arcsec`  
- `cell_arcsec`  
- `kernel`  
- `kernel_preset`  
- `kernel_sign`  
- `convsupport`  
- `gwidth_pix`, `gwidth_beam`  
- `jwidth_pix`, `jwidth_beam`  
- `truncate`  
- `support_radius_pix`, `support_radius_beam`  
- `estimator`, `n_min_avg`, `n_min_plane`

### 6.2 解決順序

本パッケージでは、幅や support を解決する際、**`gjinc` / `gauss` に限って** おおむね

`explicit pix > explicit beam > preset default`

の順に優先する。`sf` は例外であり、public には `convsupport` だけを使う。

### 6.3 `cell_arcsec=None`

`__post_init__()` と `_validate_and_resolve_config()` で

`__post_init__()` と `_validate_and_resolve_config()` で、`cell_arcsec=None` のとき `cell_arcsec = beam_fwhm_arcsec / 3` へ置換する。現在の clean 実装では、曖昧な手書き `__init__` は使っていない。

### 6.4 `kernel_preset` の正規化

- `None`, `''`, `'none'` -> `gjinc` では `mangum` へ解決  
- `'mangum2007'` -> `mangum`  
- `'legacy'` -> `casa`

### 6.5 `support_radius_*` と `truncate` の関係

- `sf` は public に `convsupport` を使う  
- `sf` に `support_radius_pix`, `support_radius_beam`, `truncate` を与える使い方は**現在実装では禁止**  
- `gjinc` / `gauss` は `truncate` または `support_radius_*` で明示指定できる  
- `sf` / `gauss` に `kernel_sign='signed'` を与えても、現在実装では warning のうえ `positive_only` に正規化される

### 6.6 `n_min_avg` と kernel の見え方

`n_min_avg` は kernel そのものの形を変えるパラメータではないが、**gridding 出力が invalid になる条件**を決めるため、実務上は kernel の見え方に強く影響する。

`estimator='avg'` のとき、各 output pixel に寄与した有効 dump 数を $N_{\rm eff}(p)$ とすると、

$$
N_{\rm eff}(p) < n_{\rm min,avg}
$$

の pixel は invalid とみなされる。現在の formal default は `n_min_avg=2` である。したがって、synthetic な 1 dump 入力では、kernel が正常でも pixel が NaN になることがある。これは kernel バグではなく QC 条件による。

---

## 7. 実装ファイルと役割

### 7.1 `map_3d/config.py`

- `MapConfig` の定義  
- public default の宣言  
- docstring による public 方針の明示

### 7.2 `map_3d/core.py`

kernel 実装の中心。

- config 検証  
- preset 解決  
- `sf`, `gjinc`, `gauss` evaluator  
- support 解決  
- nominal effective beam 推定  
- 実際の gridding 時の重み計算

### 7.3 `map_3d/gridder.py`

- kernel 解決結果のログ表示  
- effective beam summary の人間可読出力

### 7.4 `map_3d/fits_io.py`

- FITS header への kernel 情報出力  
- `KERNEL`, `KPRESET`, `KSIGN`, `CONVSUP` などを書き込む

---

## 8. `SF` kernel の理論と実装

## 8.1 `SF` の数学的イメージ

ここでの `SF` は、一般的な滑らかな positive kernel ではない。radio astronomy で使われてきた Fred Schwab の spheroidal approximation を、casacore `Sph_Conv` と同系統に評価する kernel である。

### 8.2 `SF` を単純化してはいけない点

次の理解は不正確である。

- `convsupport` は単なる切り取り半径である  
- kernel 本体は不変で、`convsupport` で切る位置だけが変わる  
- 連続 PSWF を `r/convsupport` で評価しているだけである

実際には、`convsupport` は cutoff 半径であると同時に、**どの離散 family (`isupp`, `ialpha`, `jmax`) を使うか**にも効く。

### 8.3 本パッケージの `SF` 実装

`core.py` の `_kernel_sf()` は、casacore `Sph_Conv` と同じ family selection を mirror する。

#### 8.3.1 family selection

pixel 単位の support 半径を $R_{\rm sup,pix}$ と書くと、

$$
isupp = \mathrm{clip}(\lfloor 2R_{\rm sup,pix} \rfloor, 4, 8)
$$
$$
ialpha = \mathrm{clip}(\lfloor 2\,sphparm + 1 \rfloor, 1, 5)
$$
$$
j_{\max} = \min(\lfloor R_{\rm sup,pix} \rfloor, 7)
$$

ここで本パッケージは `sphparm=1.0` 相当を採るので、通常 `ialpha=3` である。

#### 8.3.2 kernel 評価

pixel 半径 $r_{pix} = r / d$ を使って

$$
\eta = r_{pix} / j_{\max}
$$

を作り、

$$
K(r) = \mathrm{sphfn}(ialpha, isupp, \eta)
$$

を gridding branch (`IFLAG <= 0`) で評価する。

#### 8.3.3 支持半径外の扱い

$$
r_{\rm pix} > R_{\rm sup,pix}
$$

では 0 にする。

#### 8.3.4 負値の扱い

実装上、`SF` は

```python
out[out < 0] = 0.0
```

としている。すなわち、公開 API 上は **positive_only** として扱う。

### 8.4 `SPHFN` の rational approximation table

本パッケージは `MathFunc2.cc` にある Schwab の coefficient table を Python に移植している。family は

- support width: 4, 5, 6, 7, 8  
- exponent: 0, 1/2, 1, 3/2, 2

であり、piecewise な rational approximation で評価する。

### 8.5 `SF` default の意味

- formal default: `convsupport=3`  
- 実務上の別候補: `convsupport=6`

この二つを混同しないこと。

### 8.6 `convsupport` を大きくすると何が起こるか

`convsupport` を大きくすると、単純な切り取り位置だけでなく family selection と `jmax` も変わる。結果として kernel 自体が広がり、support も伸びる。直感的には「太くなる」でよいが、内部実装としては **family を切り替えたうえで広がる** と理解した方が正確である。

### 8.7 `convsupport < 3` の扱い

コード上は `convsupport >= 1` を通す。しかし実務上は `convsupport < 3` を積極的に推奨しない。

理由:

- casacore 側でも family selection は `isupp >= 4` に clamp される  
- support が短すぎると aliasing / gap / sampling 粗さへの頑健性が落ちやすい  
- default としては CASA formal default の `3` が自然

したがって、**許容される** と **推奨される** は分けて考える。

### 8.8 本パッケージで `SF` に対して実装していないもの

- public な `sphparm` 指定  
- `IFLAG > 0` の Fourier-transform branch  
- gridding correction function の別立て出力  
- CASA runtime と end-to-end で完全一致する保証

後述するように、kernel evaluator 自体は source-level にかなり寄せているが、CASA task 全体を再現したわけではない。

---

## 9. `GJINC` kernel の理論と実装

### 9.1 定義

本パッケージの `GJINC` は

$$
K(r) = \frac{J_1(\pi r / jwidth)}{\pi r / jwidth}
\exp\left[-\ln 2 \left(\frac{r}{gwidth}\right)^2\right]
$$

である。

- `jwidth`: jinc 部の幅  
- `gwidth`: Gaussian taper の HWHM

### 9.2 `r=0` の安全処理

$J_1(x)/x$ は $x=0$ で不定形になるため、コードでは

$$
\lim_{x\to 0} \frac{J_1(x)}{x} = \frac{1}{2}
$$

を使って `eps_u0` 以内を安全に評価している。

### 9.3 `mangum` preset

`gjinc + mangum` は pixel-based の文献系 preset である。

既定値:

- `jwidth_pix = 1.55`  
- `gwidth_pix = 2.52 * sqrt(log 2)`  
- `support = 3 pixel`  
- `kernel_sign = signed` (`auto` のとき)

ここで重要なのは、**現在の `mangum` は beam-aware 版ではなく pixel-based 版**であること。

### 9.4 `casa` preset

`gjinc + casa` は CASA docs に沿った safe side の preset である。

既定値:

- `jwidth_pix = 1.55`  
- `gwidth_pix = 2.52 * sqrt(log 2)`  
- `support = first null of jinc`  
- `kernel_sign = positive_only` (`auto` のとき)

### 9.5 first null

$J_1(\pi r/jwidth)=0$ の第一零点は

$$
x_1 \approx 3.83170597
$$

なので

$$
r_{\rm null} = \frac{x_1}{\pi} jwidth \approx 1.21966989\,jwidth
$$

である。`jwidth=1.55 pixel` では

$$
r_{\rm null} \approx 1.89048833\;\mathrm{pixel}
$$

となる。

### 9.6 `mangum` と `casa` の違い

- `mangum`: 3 pixel まで取り、負ローブを含み得る  
- `casa`: first null で切るため、その範囲では実質 positive-only

### 9.7 `GJINC` の物理的意味

Jinc は 2D sinc に対応する radial kernel であり、Gaussian taper を掛けることで無限 tail を抑える。signed のまま使うと高域成分を比較的保持しやすく、出力ビームの broadening を抑える方向に働くが、その代償として ringing とノイズ増幅が起こりうる。

### 9.8 `casa` preset の位置づけ

`casa` は「signed GJINC の sharpness をそのまま最大化する preset」ではなく、**CASA docs に沿った safe な GJINC** である。

---

## 10. `GAUSS` kernel の理論と実装

### 10.1 定義

$$
K(r) = \exp\left[-\ln 2 \left(\frac{r}{gwidth}\right)^2\right]
$$

### 10.2 default

- `gwidth = sqrt(log 2) pixel`  
- `support = 3 * gwidth`

### 10.3 位置づけ

`gauss` は理論的に特別な singledish optimal kernel として採っているわけではなく、シンプルで扱いやすい参照 kernel である。`sf` と `gjinc` の間の比較基準や、単純 smoothing の基準として有用である。

### 10.4 `gauss` が何をしやすいか

`gauss` は単調な positive kernel なので、直感的には「周囲の sample を滑らかに平均する」方向に働く。したがって

- ringing を出しにくい  
- aliasing を抑えやすい  
- 素朴な positive smoothing として理解しやすい

という利点がある。

### 10.5 `gauss` が第一推奨になりにくい理由

ただし、`gauss` が public default や第一推奨になっていないのは、単に broad だからではない。理由は少なくとも次の三つである。

1. **Gaussian は数学的には無限 support** なので、実装ではどこかで cutoff を決めて打ち切る必要がある。現在実装は `support = 3 * gwidth` を既定にしているが、これは実務的には妥当でも、本質的には実装上の選択である。  
2. **safe positive kernel** としては `sf(convsupport=3)` があり、こちらの方が `convsupport` による support 制御の意味が明確で、CASA/casacore 系の文脈とも対応づけやすい。  
3. **sharp kernel** としては `gjinc(mangum)` があるため、`gauss` は「安全側でも sharp 側でもない中間」に見えやすい。

したがって、`gauss` は悪い kernel なのではなく、**比較用・参照用・素朴な positive smoothing 用途では有用だが、役割の主張が `sf` や `gjinc` ほど強くない** と理解するのがよい。

### 10.6 `gauss` と `sf` の関係

`sf(convsupport=3)` も positive kernel であり、分解能だけを見れば `gauss` と大差ない場合がある。このため、`gauss` を推奨しない理由を「`sf` よりぼけるから」と単純化してはいけない。より正確には、

- `sf`: aliasing 抑制寄りの family として説明しやすく、`convsupport` で support を制御できる  
- `gauss`: 直感的で単純だが、support cutoff の意味づけと理論的な位置づけが相対的に弱い

という違いである。

### 10.7 打ち切り Gaussian を正規化すれば十分か

結論から言うと、**十分ではない**。ただし、何を保存したいかを分けて考える必要がある。

現在の gridding の基本形は

$$
T_{\rm out}(k,p)=
\frac{\sum_i K_{ip}\,q_i\,T_i(k)}{\sum_i K_{ip}\,q_i}
$$

であり、kernel 全体を定数 `c` 倍しても

$$
\frac{\sum_i (cK_{ip})\,q_i\,T_i(k)}{\sum_i (cK_{ip})\,q_i}=
\frac{\sum_i K_{ip}\,q_i\,T_i(k)}{\sum_i K_{ip}\,q_i}
$$

なので、**全体係数は打ち消し合う**。したがって、Gaussian を途中で打ち切ってから再正規化する操作で主に救えるのは

- 定数場の保存  
- 局所平均のスケール  
- DC gain

である。

一方で、再正規化しても変わらないのは

- `gauss` が positive-only であること  
- 単調な low-pass kernel であること  
- 高空間周波数を素直に落としやすいこと  
- point source の peak を下げやすいこと  
- effective beam を太らせやすいこと

である。つまり、Gaussian が第一推奨になりにくい本質は、**打ち切り誤差よりも kernel 形状そのものが smoothing 側であること**にある。

この点を整理すると次のようになる。

| 量 | 何を意味するか | 全体正規化で救えるか | 主に何で決まるか |
|---|---|---:|---|
| 定数場の保存 | 入力が一様なら出力も一様であること | はい | 正規化平均 |
| 局所平均のスケール | 平均値のスケールがずれないこと | はい | 正規化平均 |
| 点源ピーク | 点源中心値がどれだけ下がるか | いいえ | kernel 形状と support |
| effective beam | 出力ビームがどれだけ太るか | いいえ | kernel 形状と support |
| 高空間周波数の保存 | 細かい構造をどれだけ残せるか | いいえ | Fourier 空間での kernel 形 |
| ringing | 点源周囲の波紋 | いいえ | kernel の振動性と符号 |
| aliasing 抑制 | 折り返し偽構造をどれだけ防げるか | いいえ | support と Fourier 空間での減衰 |
| noise correlation | 近傍 pixel 間のノイズ相関 | いいえ | kernel 形状と広がり |

したがって、`gauss` は「打ち切り Gaussian を正規化すれば他と同じになる」わけではない。再正規化しても、それは依然として **positive low-pass kernel** であり、`sf` や `gjinc` と比べたときの役割の違いは残る。

---

## 11. effective beam の考え方とこのパッケージでの扱い

### 11.1 何で決まるか

実効ビームは、入力ビームと kernel の両方で決まる。kernel が pixel-based なので、その物理幅は `cell_arcsec` に依存する。

したがって、実効ビームは

$$
\mathrm{PSF}_{\rm out} \approx \mathrm{Beam}_{\rm in}(B) * K(d)
$$

で決まる。

### 11.2 何が間違いやすいか

誤りやすい理解:

- 実効ビームは beam size だけで決まる  
- 実効ビームは cell size だけで決まる

正しい理解:

- **両方で決まる**  
- ただし cell は kernel の物理幅を通して効く

### 11.3 本パッケージでの nominal effective beam 推定

`core.py` では、入力 Gaussian beam と radial kernel を 2D で畳み込み、中心値で正規化した nominal PSF から radial FWHM を推定する。

この推定は、内部比較と log 出力に有用だが、観測実データの sampling gap, turn-around, flagging, 重み変動をすべて反映した empirical beam とは別物である。

---

## 12. `SF` と `GJINC` の内部比較結果

ここでの数値は、**現在の本パッケージ実装内部での比較**であり、CASA runtime との end-to-end 1:1 比較ではない。

### 12.1 比較条件

- 入力ビーム: $B = 350$ arcsec  
- cell size: $d = B/3 = 116.667$ arcsec/pixel  
- 比較 kernel:  
  - `SF(convsupport=3)`  
  - `SF(convsupport=6)`  
  - `GJINC(casa)`  
  - `GJINC(mangum)`

### 12.2 nominal output PSF FWHM

| kernel | nominal output PSF FWHM [arcsec] | input beam に対する増加率 |
|---|---:|---:|
| `GJINC(mangum)` | 382.839 | +9.383 % |
| `GJINC(casa)`   | 400.666 | +14.476 % |
| `SF(3)`         | 426.225 | +21.779 % |
| `SF(6)`         | 563.252 | +60.929 % |

### 12.3 `SF` が `GJINC` よりどれだけ大きいか

| compare | 差 [arcsec] | 比率 |
|---|---:|---:|
| `SF(3)` vs `GJINC(casa)`   | 25.559  | +6.379 % |
| `SF(6)` vs `GJINC(casa)`   | 162.587 | +40.579 % |
| `SF(3)` vs `GJINC(mangum)` | 43.386  | +11.333 % |
| `SF(6)` vs `GJINC(mangum)` | 180.413 | +47.125 % |

### 12.4 解釈

この条件では、`SF(3)` は `GJINC(casa)` より少し broad だが、default 候補としてはまだ理解可能である。一方 `SF(6)` は support 半径が

$$
R_{cut} = 6d = 2B
$$

となるため、単一鏡 OTF では broadening がかなり大きい。したがって、**`SF(6)` を formal default にしない** という判断は、この internal comparison と整合している。

### 12.5 重要な注意

ここでの比較は package 内部の nominal 比較である。同一 synthetic OTF data を CASA `tsdimaging` に通した end-to-end 比較は、まだ完了していない。

### 12.6 kernel の善し悪しを何で定量化するか

kernel の比較を曖昧な言葉だけで済ませると誤解しやすい。少なくとも次の量を分けて評価するのが望ましい。

1. **effective beam FWHM**  
   出力マップの実効ビームがどれだけ太るか。
2. **peak response**  
   点源や細い構造のピークがどれだけ残るか。
3. **ringing / side-lobe**  
   強い点源や急峻な境界の周囲に波紋状の偽構造が出るか。
4. **aliasing 抑制**  
   Nyquist を超える高空間周波数が低周波へ折り返して偽構造になるのをどれだけ防げるか。
5. **noise correlation / RMS**  
   出力ノイズがどれだけ平均化されるか、近傍画素間でどれだけ相関するか。
6. **欠損・edge・粗い sampling に対する頑健さ**  
   sample 間隔が粗い、coverage が薄い、row が flag で落ちる、edge に近い、といった条件でどれだけ破綻しにくいか。

本パッケージの kernel 比較では、少なくとも `effective beam`, `peak response`, `ringing`, `aliasing`, `noise correlation` を意識すること。

### 12.7 `ringing` とは何か

`ringing` とは、強い点源や急峻な境界の周囲に現れる**波紋状の偽構造**である。典型的には

- 点源の周囲に正負の輪が見える  
- 点源近傍に不自然な負の谷が出る  
- 細い構造の両側に overshoot / undershoot が出る

という見え方をする。

直感的には、**より sharp に復元しようとする kernel ほど、少し引き算を含みやすく、その引き算が波紋に見える** と考えるとよい。したがって

- `sf`, `gauss` は ringing が比較的小さい  
- `gjinc(mangum, signed)` は ringing が出やすい

という理解でよい。

### 12.8 `aliasing` 抑制とは何か

`aliasing` とは、**本来は表現できないほど細かい空間周波数成分**が、低い空間周波数へ折り返して**偽構造**として現れることである。OTF では dump 間隔や pixel 間隔が有限なので、すべての高周波を表現できるわけではない。

したがって kernel には、少なくとも一部、

- 表現不能な高周波を無理に残しすぎない  
- Nyquist 近傍の成分を穏やかに落とす

という役割がある。

この観点では

- `sf`: aliasing 抑制寄り  
- `gauss`: 強い smoothing により aliasing も抑えやすい  
- `gjinc(mangum)`: 分解能を残しやすい代わりに aliasing/ringing 側では注意が必要

と理解するとよい。

### 12.9 `gjinc` は何を保存したい kernel か

`gjinc(mangum, signed)` が主に保存したいのは、

- 入力ビームが本来持っている空間分解能  
- point source の peak response  
- 高空間周波数成分

である。

言い換えると、**望遠鏡が本来持っていた細かい構造を、gridding の段階でなるべく鈍らせずに残したい** kernel である。負ローブを含めるのは、そのために少しデコンボリューション的に働かせるためであり、その代償として

- ringing  
- ノイズ増幅  
- sparse sampling での負アーティファクト

が起こりうる。

### 12.10 粗い sampling のとき何が問題になるか

実際の dump 間隔を $\Delta s$ とする。$\Delta s$ が入力ビームや cell に対して十分細かくないとき、gridding の問題は主に次の二つに集約される。

1. 高空間周波数を十分に表現できず、aliasing が増える。  
2. signed kernel では ringing や負アーティファクトが出やすくなる。

このため、粗い sampling が疑われるときは

- まず `sf + convsupport=3` で safe な map を作る  
- 次に `gjinc + mangum` へ切り替えて sharpness と副作用を比較する

という順序が安全である。

### 12.11 `gauss` の位置づけをどう理解すべきか

`gauss` は

- ringing を出しにくい  
- aliasing を抑えやすい  
- 直感的で単純

という利点があるため、粗い sampling で最も危険な kernel ではない。むしろ `gjinc(mangum)` よりは safe 側である。

一方で、

- 高空間周波数をより素直に落としやすい  
- support cutoff に自然な基準が弱い  
- safe positive kernel としては `sf` が近い役割を持つ

ため、**悪い kernel ではないが第一推奨にもなりにくい** と理解するのが正確である。

---

## 13. quick reference

| kernel 設定 | 符号 | support の概念 | broadening 傾向 | 主な長所 | 主な注意点 | 主な用途 |
|---|---|---|---|---|---|---|
| `sf`, `convsupport=3` | positive-only | `convsupport * cell`。family selection を伴う | 軽〜中程度 | 安全, 破綻しにくい, aliasing に比較的強い | `GJINC` より broad | 1st look, 標準生成 |
| `sf`, `convsupport=6` | positive-only | 同上 | 大 | 強い平滑化, 安定 | 単一鏡では broadening が大きい | 明示指定の特殊用途 |
| `gjinc`, `preset='mangum'` | signed | 3 pixel cutoff | 最も sharp | 分解能寄り | ノイズ増幅, ringing | 最終解析, sharpness 重視 |
| `gjinc`, `preset='casa'` | positive-only | first null cutoff | 中程度 | safe な GJINC | `mangum` より broad | 安全寄り GJINC |
| `gauss` | positive-only | `3*gwidth` | 設定依存 | 単純, ringing が小さい, 参照用 | support cutoff は実装上の選択, safe 側の主役は `sf` | 比較・参照 |

---

## 14. 推奨運用

### 14.1 標準運用（初期解析）

```python
config = MapConfig(
    x0=..., y0=..., nx=..., ny=...,
    beam_fwhm_arcsec=B,
    cell_arcsec=None,      # -> B/3
    kernel='sf',
    convsupport=3,
)
```

意味:

- formal default と整合  
- pixel-based で理解しやすい  
- aliasing / gap / local instability に比較的強い  
- まずは破綻しにくいマップを得る

### 14.2 sharpness 重視の再イメージング

```python
config = MapConfig(
    x0=..., y0=..., nx=..., ny=...,
    beam_fwhm_arcsec=B,
    cell_arcsec=B/3,
    kernel='gjinc',
    kernel_preset='mangum',
    kernel_sign='signed',
)
```

意味:

- 分解能寄り  
- 負ローブを含める  
- ノイズ増幅や ringing を評価しながら使う

### 14.3 safe な GJINC

```python
config = MapConfig(
    x0=..., y0=..., nx=..., ny=...,
    beam_fwhm_arcsec=B,
    cell_arcsec=B/3,
    kernel='gjinc',
    kernel_preset='casa',
)
```

意味:

- `gjinc` family を使いたいが、安全寄りにしたい  
- first-null cutoff を採る  
- `mangum` より broad だが扱いやすい

### 14.4 `convsupport=6` を使う場合

```python
config = MapConfig(
    x0=..., y0=..., nx=..., ny=...,
    beam_fwhm_arcsec=B,
    cell_arcsec=B/3,
    kernel='sf',
    convsupport=6,
)
```

前提:

- broadening を許容する  
- 安定性・平滑性を優先する  
- formal default ではなく明示指定として使う

### 14.5 実務的な 2 段階運用

1. まず `sf + convsupport=3` で全体像とノイズを確認する  
2. sharpness が重要なら `gjinc + mangum` で再イメージングする  
3. `WEIGHT`, `CANCEL` などの診断量を見てトレードオフを判断する

---

## 15. 間違いやすい点

### 15.1 `convsupport` は単なる切り取り位置ではない

`SF` では `convsupport` は cutoff 半径であると同時に family selection にも効く。切る場所だけが変わるわけではない。

### 15.2 `SF(6)` は CASA docs に出てくるが、formal default ではない

CASA docs の formal default は `convsupport=3` である。single-dish imaging notebook の `6` は ALMA data に対する推奨例であり、両者は区別する必要がある。

### 15.3 `mangum` は beam-aware ではない

現在の `mangum` は pixel-based である。以前の beam-aware 解釈は public default には採用していない。

### 15.4 `legacy` は独立 preset ではない

現在 `legacy` は互換 alias であり、実体は `casa` である。主名として使わない。

### 15.5 実効ビームは `beam` と `cell` の両方で決まる

`beam` だけでも `cell` だけでもない。`cell` は kernel の物理幅を通して効く。

### 15.6 `SF` の現在実装は generic continuous PSWF ではない

casacore の離散 Schwab family を意識した実装である。`SPHFN` table に基づく。

### 15.7 `GJINC(casa)` と `GJINC(mangum)` は符号だけが違うわけではない

support の既定値も違う。

- `mangum`: 3 pixel  
- `casa`: first null

### 15.8 `positive_only` は単に「安全」ではあるが、sharpness を落としやすい

負ローブを落とすと ringing は減るが、高域保持も弱くなる。`positive_only` は常に優れているわけではない。

### 15.9 `convsupport < 3` はコード上許容されても、推奨とは限らない

casacore family selection が clamp されるため、support を小さくしすぎると概念上も実務上も不自然になりやすい。

### 15.10 1 dump だけの synthetic 入力で全部 NaN でも、まず `n_min_avg` を疑う

現在の default は `n_min_avg=2` である。したがって 1 dump しか寄与しない pixel は invalid になり、全 kernel で NaN になり得る。これは `sf`, `gjinc`, `gauss` のどれかが壊れていることを直接意味しない。

---


### 15.11 `GridInput.flag` は kernel パラメータではないが、kernel の見え方に影響する

`GridInput.flag` は現在の clean 実装では **1D の row mask** である。

- shape: `(ndump,)`
- 意味: `True = use`, `False = drop`
- 互換入力として `0/1`, `0.0/1.0` は受け付ける
- `-1`, `2`, `0.5` などは invalid として error
- 2D `(ndump, nchan)` flag は現在の公開仕様では扱わない

ここで重要なのは、`flag` は **行全体の採否** を決めるものであり、channel ごとの bad mask ではないという点である。したがって

- dump 行全体を落としたい -> `flag[i] = False`
- 一部 channel だけ落としたい -> `spec[i, bad_channels] = NaN`

という役割分担になる。

この仕様は kernel そのものの数式を変えない。しかし、実際の gridding では

$$
N_{\rm used}(p)
$$

すなわちその pixel に寄与する有効 dump 数を変えるため、特に sparse input や edge では **見かけの kernel broadening / invalid pixel の出方 / cancellation diagnostics** に影響する。

### 15.12 `n_min_avg` と `GridInput.flag` / `spec=NaN` の関係

`n_min_avg` は

$$
N_{\rm used}(p) < n_{\rm min,avg}
$$

の pixel を invalid にする QC 条件である。`GridInput.flag=False` の行や、必要 channel が `NaN` の行は、結果として $N_{\rm used}(p)$ を減らす。

したがって、例えば

- `sf + convsupport=3`
- `gjinc + mangum`
- `gjinc + casa`
- `gauss`

のどの kernel を使っていても、入力が 1 dump しか寄与しない synthetic 条件では、formal default の `n_min_avg=2` により全 pixel が `NaN` になり得る。これは kernel バグではなく QC と入力条件の結果である。

### 15.13 `basketweave` / `ps_gridder` でも同じ `flag` 仕様を使う

現在の clean 実装では、`GridInput.flag` の 1D bool row mask 仕様は `core.grid_otf(...)` だけでなく、`basketweave` の geometry / merge 系、および `ps_gridder` にもそろえて適用している。

したがって、manual 上でも

- `GridInput.flag`: row 単位の採否
- channel 単位の欠損: `spec` の `NaN`

を共通仕様として読んでよい。

### 15.14 `support` という語は一般名であり、`truncate` / `convsupport` はその kernel ごとの呼び方である

`support` は直感的には **kernel を使う最大距離**、すなわち **有効半径**または **打ち切り半径**である。

- `support radius`: 一般概念  
- `truncate`: `gauss` / `gjinc` での public 名  
- `convsupport`: `sf` での public 名  
- `first null`: `gjinc` での自然な cutoff 候補

という関係であり、全部がまったく別物というわけではない。ただし、`sf` の `convsupport` は単なる切り取り位置ではなく、family selection と `jmax` にも効くので、`gauss` のような「無限に広がる kernel をどこで切るか」とは意味が少し異なる。

### 15.15 `gauss` が第一推奨でないのは、途中で打ち切るからだけではない

Gaussian は数学的には無限 support なので、実装ではどこかで打ち切る必要がある。これは事実である。しかし、第一推奨になりにくい本質はそれだけではない。打ち切って再正規化しても、`gauss` は依然として

- positive-only  
- 単調な low-pass  
- 高空間周波数を素直に落としやすい

という性格を持つ。したがって、**定数場や平均値スケールは保てても、point source peak, effective beam, sharpness, aliasing / ringing の性格は kernel 形状で決まり、再正規化では変わらない**。`sf` が default なのは、単に `gauss` より sharp だからではなく、support 制御の意味と aliasing 抑制寄りの位置づけがより明確だからである。



## 16. 何を実装しているか

### 16.1 実装済み

- `sf`, `gjinc`, `gauss` の 3 kernel  
- `gjinc` に対する `mangum` / `casa` preset  
- 互換 alias: `mangum2007 -> mangum`, `legacy -> casa`  
- `sf` に対する `convsupport` public API  
- `cell_arcsec=None -> beam/3` の auto default  
- `MapConfig` を通常の dataclass として維持し、`GridConfig = MapConfig` を互換 alias とする構成  
- nominal effective beam estimation  
- FITS header への `KERNEL`, `KPRESET`, `KSIGN`, `CONVSUP` 書き込み  
- `SF` の `Sph_Conv` family selection と `SPHFN` coefficient table の移植  
- `gjinc` の CASA docs 準拠の pixel-based parameterization

### 16.2 source-level にかなり寄せているもの

- `SF` の `ialpha`, `isupp`, `jmax` 解決  
- `SPHFN` の gridding branch (`IFLAG <= 0`)  
- `GJINC` の public default (`jwidth`, `gwidth`, `truncate`)  
- `cell=beam/3` という default 方針

---

## 17. 何を実装していないか

### 17.1 CASA task 全体の完全複製

このパッケージは CASA task 全体の clone ではない。特に次はそのまま再現していない。

- `BOX` kernel  
- `PB` kernel  
- CASA 独自の文字列 unit parser を含む UI 全体  
- task 層の全オプション  
- CASA image / MS 周辺のエコシステム

### 17.2 `SF` の public `sphparm`

`sphparm` は現在 public に出していない。実質 `1.0` 固定である。

### 17.3 `SPHFN` の Fourier-transform branch

`IFLAG > 0` branch は実装していない。本パッケージは gridding branch のみを使う。

### 17.4 gridding correction function

gridding correction function を独立にユーザーへ公開したり、別出力として扱ったりはしていない。

### 17.5 CASA runtime との end-to-end 1:1 検証

kernel evaluator 自体は source-level にかなり寄せたが、synthetic OTF input を CASA `tsdimaging` に通した end-to-end 比較はまだ完了していない。

### 17.6 `sincgauss`

現時点では未実装。

### 17.7 kernel 最適化の自動探索

`cell/beam` に応じて最適 kernel width や support を自動最適化する機能はない。

---

## 18. このパッケージで何を参考に実装したか

### 18.1 `SF`

- casacore `Sph_Conv` class reference  
- `MathFunc.tcc` の `Sph_Conv::value()`  
- `MathFunc2.cc` の `SPHFN` coefficient table と gridding branch

### 18.2 `GJINC`

- CASA `tsdimaging` docs の `jwidth`, `gwidth`, `truncate` default  
- Mangum et al. 2007 の kernel parameterization

### 18.3 cell default と singledish 背景

- CASA `tsdimaging` docs  
- CASA single-dish imaging notebook の理論説明と ALMA singledish の推奨例

---

## 19. この説明書の利用者への要点

### 19.1 まず安全に作るなら

- `kernel='sf'`  
- `convsupport=3`  
- `cell_arcsec=None` または `beam/3`

### 19.2 分解能を詰めるなら

- `kernel='gjinc'`  
- `kernel_preset='mangum'`  
- `kernel_sign='signed'`

### 19.3 `convsupport=6` を見たら

- それは formal default ではない  
- singledish の実務例としては存在する  
- しかしこの package 内部比較ではかなり broad だった

### 19.4 `legacy` を見たら

- 現在は `casa` alias と読む

### 19.5 `mangum2007` を見たら

- 現在は `mangum` alias と読む

---

## 20. 現時点の最重要結論

1. **公開 preset はすべて pixel-based** である。  
2. **formal default は `sf + convsupport=3`** である。  
3. `convsupport=6` は背景的・実務的には重要だが、formal default ではない。  
4. `SF` は generic smooth kernel ではなく、**casacore の discrete Schwab family** を意識した実装である。  
5. `gjinc + mangum` は sharp だが signed kernel の副作用を理解して使う必要がある。  
6. `gjinc + casa` は safe な GJINC である。  
7. `legacy` と `mangum2007` は互換 alias であり、主名ではない。  
8. まだ未完なのは、**CASA runtime との end-to-end 1:1 比較**である。

---

## 21. 参考資料

### 21.1 文献

- Mangum, Emerson, Greisen (2007), *The On The Fly Imaging Technique*  
- Sawada et al. (2008), *On-The-Fly Observing System of the Nobeyama 45-m and ASTE 10-m Telescopes*

### 21.2 CASA / casacore

- CASAdocs `tsdimaging`  
- CASAdocs `Single-Dish Imaging` notebook  
- casacore `Sph_Conv` class reference  
- casacore `MathFunc.tcc`  
- casacore `MathFunc2.cc`

### 21.3 本パッケージ内資料

- OTF kernel handover 詳細版  
- `SF vs GJINC` 比較メモ  
- kernel 仕様書  
- 現在の patch 一式
