# `estimate_relative_beam_intensities` 詳細説明書

## 1. 概要

`estimate_relative_beam_intensities` は、複数の `OTFBundle` から **ビームごとの相対強度** を推定する関数です。

本関数の主眼は、

- 各 bundle の **共通に比較できる空間領域** だけを使うこと
- OTF gridding 後の `support_mask` / `valid_mask` / `NaN` を尊重すること
- baseline subtraction 後に既に得られている `SIGNAL_MASK_USED` / `SIGNAL_MASK3D_USED` / `MOMENT0` / `BASE_RMS` を優先して使うこと
- 出力が **入力順と 1 対 1 に対応** すること

です。

特に重要なのは、**出力の主識別子は `input_index` と `input_label` であり、`FDNUM/IFNUM/PLNUM` ではない**ことです。

これは、OTF gridding や beam convolution の途中で `FDNUM/IFNUM/PLNUM` を引き継がない設計があり得るためです。したがって、本関数は **入力順を保存したまま** 結果を返します。

---

## 2. 想定している問題設定

各入力 bundle を

$$
B_i \quad (i = 0, 1, \dots, N-1)
$$

とします。

各 bundle から比較用の 2 次元積分強度マップ

$$
M_i(y, x)
$$

を作り、共通比較領域 $Q$ の上で、各ビームの代表強度を求めます。

最終的に返したい量は

$$
R_i = \frac{S_i}{\max_j S_j}
$$

です。

ここで

- $S_i$ : bundle $i$ の生の代表強度
- $R_i$ : 最大値で正規化した相対強度

です。

`method="sum"` では $S_i$ は背景補正後の総和です。

`method="fit"` では $S_i$ は基準 bundle に対する fitted slope です。

---

## 3. 関数シグネチャ

```python
estimate_relative_beam_intensities(
    bundles,
    *,
    input_labels=None,
    method="sum",
    edge_margin_beam=1.5,
    velocity_range_kms=None,
    min_snr=5.0,
    min_pixels=16,
    positive_only=True,
    min_template_fraction=0.2,
    relative_weight_threshold=0.10,
    clip_sigma=4.0,
    max_clip_iter=5,
)
```

返り値は `astropy.table.Table` ではなく、`BeamIntensityEstimateResult` オブジェクトです。

---

## 4. 入力引数の意味

### 4.1 `bundles`

比較したい `OTFBundle` の列です。

要件は次です。

- 空間グリッドが全 bundle で一致していること
- `method` や利用可能な生成物に応じて、必要な `image_ext` や WCS 情報があること

空間グリッド一致の確認では、少なくとも次を見ています。

- `data.shape[1:]`
- `CTYPE1`, `CTYPE2`
- `CRPIX1`, `CRPIX2`
- `CRVAL1`, `CRVAL2`
- `CDELT1`, `CDELT2`
- `CUNIT1`, `CUNIT2`

比較は単位正規化後に行います。たとえば `deg` と `arcsec` の違いは吸収して比較します。

### 4.2 `input_labels`

入力順を人間にわかりやすく表示するラベルです。

例:

```python
input_labels=["bundle2RU_bl", "bundle3LU_bl", "bundle4LU_bl", "bundle5LU_bl"]
```

指定しない場合は自動で

- `bundle_0`
- `bundle_1`
- ...

になります。

### 4.3 `method`

- `"sum"` : 総和ベース
- `"fit"` : robust affine fit ベース

通常運用では、まず `"sum"` を基本とし、確認用に `"fit"` も見るのが推奨です。

### 4.4 `edge_margin_beam`

各 bundle の edge を、ビーム FWHM の何本分だけ内側へ削るかを指定します。

削る半径ピクセル数は

$$
r_{\mathrm{pix}} = \left\lceil \frac{m \times \theta_{\mathrm{beam}}}{\min(\Delta_x, \Delta_y)} \right\rceil
$$

です。

- $m$ : `edge_margin_beam`
- $\theta_{\mathrm{beam}}$ : beam FWHM [arcsec]
- $\Delta_x, \Delta_y$ : 画素サイズ [arcsec/pixel]

beam FWHM は優先的に

- header の `BMAJ`, `BMIN`
- `meta` の `beam_fwhm_arcsec` など

から読みます。

beam 情報が無ければ edge trim は行わず、warning に

- `beam_fwhm_missing:no_edge_trim`

が入ります。

### 4.5 `velocity_range_kms`

比較用マップを cube から作るときの速度範囲です。

これは **`SIGNAL_MASK_USED` / `SIGNAL_MASK3D_USED` が無い場合の明示的 fallback** です。

本関数は、弱い自動推定には落としません。つまり、比較用マップの優先順位は次です。

1. 共通 `SIGNAL_MASK_USED` / `SIGNAL_MASK3D_USED`
2. `velocity_range_kms`
3. 既存 `MOMENT0`

このどれも使えなければエラーです。

### 4.6 `min_snr`

比較領域に採用する pixel の S/N 下限です。

`sigma_map` が全 bundle で構成できる場合にのみ使います。そうでない場合は S/N cut を行いません。

### 4.7 `min_pixels`

最終的な比較領域の最小 pixel 数です。これ未満ならエラーです。

### 4.8 `positive_only`

- `True` : 放射線を仮定し、正の強度を主に扱う
- `False` : 吸収線も想定し、絶対値ベースで比較する

### 4.9 `min_template_fraction`

比較領域を作るとき、template map の 95 パーセンタイル値を基準に、どの程度の強度以上を採用するかを指定します。

閾値は

$$
T(y, x) \ge f \times T_{95}
$$

です。

- $T(y, x)$ : template map
- $T_{95}$ : template の 95 パーセンタイル
- $f$ : `min_template_fraction`

### 4.10 `relative_weight_threshold`

`WEIGHT_SUM` / `WEIGHT` / `WSUM` がある場合、それらの最大値に対して相対的に低すぎる画素を support から除外する閾値です。

### 4.11 `clip_sigma`, `max_clip_iter`

`method="fit"` のときに使う反復外れ値除去の設定です。

---

## 5. 比較に使う情報の優先順位

### 5.1 最優先: 共通 `SIGNAL_MASK`

全 bundle に `SIGNAL_MASK_USED` または `SIGNAL_MASK3D_USED` があれば、それらの **共通 intersection** を使います。

この場合、各 bundle の比較用マップは

$$
M_i(y, x) = \sum_{k \in \mathcal{S}_{\cap}} T_i(k, y, x)\, \Delta v
$$

で作られます。

- $T_i(k, y, x)$ : bundle $i$ の cube
- $\mathcal{S}_{\cap}$ : 全 bundle の共通 signal mask
- $\Delta v$ : チャネル幅 [km/s]

`SIGNAL_MASK` を共通利用する場合、全 bundle の速度軸は一致している必要があります。

確認対象は少なくとも

- `NAXIS3`
- `CTYPE3`
- `CRPIX3`
- `CRVAL3`
- `CDELT3`
- `CUNIT3`

です。

### 5.2 次点: `velocity_range_kms`

`SIGNAL_MASK` が無いときに、ユーザーが明示した速度範囲で比較用マップを作ります。

$$
M_i(y, x) = \sum_{v_{\min} \le v_k \le v_{\max}} T_i(k, y, x)\, \Delta v
$$

### 5.3 最後の fallback: 既存 `MOMENT0`

全 bundle に `MOMENT0` が入っていればそれを使います。

ただしこの場合、bundle ごとに signal integration window が違う可能性があり、gain 差と窓の違いが分離できません。そのため warning に

- `noncommon_signal_mask_risk`

が入ります。

---

## 6. OTF gridding の有効領域の扱い

本関数は配列サイズだけで edge を決めません。

空間的な有効領域は、まず各 bundle ごとに次の順で決めます。

1. `bundle.support_mask`
2. `image_ext["SUPPORT_MASK"]`
3. `image_ext["MASK"]`
4. `np.any(np.isfinite(bundle.data), axis=0)`

さらに `WEIGHT_SUM` / `WEIGHT` / `WSUM` があれば、相対閾値で support を削ります。

最後に `valid_mask` と `NaN` を見て、実際に有効な空間 footprint を作ります。

つまり、各 bundle の空間 footprint は概念的には

$$
F_i = U_i \cap V_i
$$

です。

- $U_i$ : support 側から得た空間領域
- $V_i$ : `valid_mask` と finite data を反映した空間領域

edge trim 後の有効領域を $E_i$ とすると、全 bundle の共通 support は

$$
Q_{\mathrm{support}} = \bigcap_i E_i
$$

です。

---

## 7. 比較領域の作り方

比較領域は `common_support` の中からさらに絞ります。

### 7.1 template map

まず、各 bundle の比較用マップ $M_i(y, x)$ の median を template とします。

$$
T(y, x) = \mathrm{median}_i\, M_i(y, x)
$$

### 7.2 強度閾値

`positive_only=True` のときは、正の template のみを使い、95 パーセンタイル値 $T_{95}$ を基準に

$$
T(y, x) \ge f T_{95}
$$

を満たす pixel だけを残します。

さらに孤立した極端値を避けるため、

$$
T(y, x) \le 5 T_{95}
$$

も課します。

`positive_only=False` のときは $|T(y, x)|$ に対して同じ処理を行います。

### 7.3 S/N 閾値

全 bundle で `sigma_map` がある場合は、各 bundle ごとに

$$
\mathrm{SNR}_i(y, x) = \frac{M_i(y, x)}{\sigma_i(y, x)}
$$

あるいは `positive_only=False` なら

$$
\mathrm{SNR}_i(y, x) = \frac{|M_i(y, x)|}{\sigma_i(y, x)}
$$

を作り、その median を基準に S/N cut を掛けます。

最終比較領域を $Q$ とします。

---

## 8. `sigma_map` の作り方

`BASE_RMS` がある場合、積分強度の不確かさは

$$
\sigma_i(y, x) = \sigma_{\mathrm{base}, i}(y, x)\, \Delta v\, \sqrt{N_i(y, x)}
$$

で近似します。

- $\sigma_{\mathrm{base}, i}(y, x)$ : `BASE_RMS`
- $\Delta v$ : チャネル幅 [km/s]
- $N_i(y, x)$ : 実際に finite だった選択チャネル数

重要なのは、$N_i$ は **mask に入ったチャネル数** ではなく、**mask に入りかつ finite だったチャネル数** であることです。

---

## 9. `method="sum"` の定義

`sum` では、比較領域 $Q$ の外側で、かつ共通 support の内側に背景候補領域を作ります。

$$
B = Q_{\mathrm{support}} \setminus Q
$$

背景候補 pixel 数が `min_pixels` 以上なら、bundle ごとに背景中央値

$$
c_i = \mathrm{median}_{(y,x) \in B} M_i(y, x)
$$

を引きます。足りない場合は

$$
c_i = 0
$$

です。

その上で、生の代表強度は

$$
S_i^{(\mathrm{sum})} = \sum_{(y,x) \in Q} \left[M_i(y, x) - c_i\right]
$$

です。

不確かさは

$$
\delta S_i^{(\mathrm{sum})} = \sqrt{\sum_{(y,x) \in Q} \sigma_i(y, x)^2}
$$

で近似します。

この方法の特徴は次です。

- 実装が比較的単純
- 多 pixel を使うので peak 比より安定
- 小さな加法オフセットに比較的強い

---

## 10. `method="fit"` の定義

`fit` では、まず基準 bundle を 1 本選びます。

基準選択は、比較領域 $Q$ の中で背景補正後の値の 90 パーセンタイルが最大の bundle を選びます。

概念的には

$$
\mathrm{metric}_i = \mathrm{p90}_{(y,x) \in Q} \left(M_i(y, x) - c_i\right)
$$

です。

`positive_only=False` のときは絶対値を使います。

基準 bundle を $r$ とすると、各 bundle $i$ に対して

$$
M_i(y, x) \approx a_i M_r(y, x) + b_i
$$

を robust に fit します。

- $a_i$ : multiplicative scale
- $b_i$ : additive offset

### 10.1 重み

`sigma_map` があれば、近似的に

$$
\sigma_{\mathrm{eff}}^2 \approx \sigma_y^2 + a^2 \sigma_x^2
$$

を使い、

$$
w = \frac{1}{\sigma_{\mathrm{eff}}^2}
$$

とします。

### 10.2 外れ値除去

fit 後の残差に対して median と MAD を使い、

$$
|r - \mathrm{median}(r)| \le \kappa \times 1.4826 \times \mathrm{median}\left(|r - \mathrm{median}(r)|\right)
$$

を満たす点のみ残して反復します。

- $\kappa$ : `clip_sigma`

### 10.3 ratio fallback

基準 map 側の dynamic range がほとんど無く、切片付き fit が退化するときは

$$
M_i(y, x) \approx a_i M_r(y, x)
$$

の ratio fit に落とします。

このとき warning に

- `fit:ratio_fallback`

が入ります。

### 10.4 `fit` の `raw_strength`

`method="fit"` のときの `raw_strength` は、積分総和ではなく fitted slope $a_i$ です。

したがって、`sum` と `fit` では `raw_strength` の意味が異なります。

- `sum` : 背景補正後総和
- `fit` : 基準 bundle に対する slope

---

## 11. 相対強度の定義

最終的に返す相対強度は

$$
R_i = \frac{|S_i|}{\max_j |S_j|}
$$

ではなく、`positive_only` に応じて次のように定義します。

### 11.1 `positive_only=True`

$$
R_i = \frac{S_i}{\max_j S_j}
$$

このとき、最大 raw strength が 0 以下ならエラーです。

### 11.2 `positive_only=False`

$$
R_i = \frac{|S_i|}{\max_j |S_j|}
$$

です。

---

## 12. 返り値

返り値は `BeamIntensityEstimateResult` オブジェクトです。

### 12.1 `print(result)`

summary table だけを表示します。

### 12.2 `result.summary`

主な出力表です。**出力順は入力順のまま** です。

主な列:

- `input_index`
- `input_label`
- `relative_strength`
- `relative_strength_err`
- `raw_strength`
- `raw_strength_err`
- `comparison_pixels`
- `support_pixels`
- `is_strongest`
- `warning_flag`

### 12.3 `result.detail`

詳細情報を含む表です。summary の列に加えて、少なくとも次が入ります。

- `fdnum`
- `ifnum`
- `plnum`
- `stream_name`
- `method`
- `moment_source`
- `weighting_mode`
- `fit_mode`
- `fit_offset`
- `offset_correction`
- `is_reference_for_fit`
- `warning_flag`

ここで `fdnum/ifnum/plnum/stream_name` は **補助情報** です。主識別子ではありません。

### 12.4 `result.config`

実行条件をまとめた辞書です。

主なキー:

- `method`
- `moment_source`
- `input_labels`
- `edge_margin_beam`
- `velocity_range_kms`
- `min_snr`
- `min_pixels`
- `positive_only`
- `min_template_fraction`
- `template_ref`
- `template_upper_cap_factor`
- `relative_weight_threshold`
- `sum_offset_mode`
- `weighting_mode`
- `clip_sigma`
- `max_clip_iter`
- `comparison_pixels`
- `support_pixels`
- `strongest_input_index`
- `strongest_input_label`
- `fit_reference_input_index`
- `fit_reference_input_label`

---

## 13. `warning_flag` の意味

現行実装で summary/detail に残る warning は次のようなものです。

### 13.1 `beam_fwhm_missing:no_edge_trim`

beam FWHM が取れず、edge trim ができなかったことを意味します。

### 13.2 `noncommon_signal_mask_risk`

比較用マップが既存 `MOMENT0` から作られており、beam ごとに積分窓が違う可能性があることを意味します。

### 13.3 `fit:ratio_fallback`

`method="fit"` で切片付き affine fit が退化し、ratio fit に落ちたことを意味します。

---

## 14. 使い方の推奨

### 14.1 通常運用

まずは `method="sum"` を使うのが推奨です。

```python
result = estimate_relative_beam_intensities(
    [bundle2RU_bl, bundle3LU_bl, bundle4LU_bl, bundle5LU_bl],
    input_labels=["bundle2RU_bl", "bundle3LU_bl", "bundle4LU_bl", "bundle5LU_bl"],
    method="sum",
    edge_margin_beam=1.5,
    min_snr=5.0,
    min_pixels=16,
    positive_only=True,
)

print(result)
```

### 14.2 確認用に `fit` も見る

```python
result_fit = estimate_relative_beam_intensities(
    [bundle2RU_bl, bundle3LU_bl, bundle4LU_bl, bundle5LU_bl],
    input_labels=["bundle2RU_bl", "bundle3LU_bl", "bundle4LU_bl", "bundle5LU_bl"],
    method="fit",
    edge_margin_beam=1.5,
    min_snr=5.0,
    min_pixels=16,
    positive_only=True,
)
```

### 14.3 `SIGNAL_MASK` が無く、cube だけある場合

```python
result = estimate_relative_beam_intensities(
    bundles,
    input_labels=labels,
    method="sum",
    velocity_range_kms=(-10.0, 20.0),
    edge_margin_beam=1.5,
)
```

この場合、速度範囲の指定は必須です。

---

## 15. どう解釈すればよいか

### 15.1 `sum` と `fit` が近い

理想的です。比較領域の定義と gain 推定の両方が安定している可能性が高いです。

### 15.2 `sum` と `fit` が大きく違う

次を疑うべきです。

- residual baseline
- `MOMENT0` の積分窓の不一致
- 比較領域が狭すぎる
- 外れ値 pixel の影響
- 基準 bundle の dynamic range 不足

### 15.3 `existing_moment0` しか使えない

比較は可能ですが、`noncommon_signal_mask_risk` が出るなら、beam ごとに異なる signal integration window が混ざっている可能性があります。可能なら共通 `SIGNAL_MASK` か `velocity_range_kms` に戻した方が安全です。

---

## 16. 制約と注意

1. **空間グリッドは全 bundle で一致している必要があります。**
2. **共通 `SIGNAL_MASK` を使うなら速度軸も一致している必要があります。**
3. **弱い自動推定は実装していません。**
4. **`fdnum/ifnum/plnum` は主識別子ではありません。**
5. **出力順は入力順のままです。**
6. `method="fit"` の `raw_strength` は slope であり、`sum` の `raw_strength` とは意味が違います。

---

## 17. まとめ

`estimate_relative_beam_intensities` の設計思想は次の通りです。

- OTF gridding 後の有効領域を尊重する
- ビームごとに edge を削った上で、共通空間領域だけを使う
- 可能なら共通 `SIGNAL_MASK` を使う
- それが無ければ `velocity_range_kms` を明示的に指定する
- 出力は必ず入力順に対応させる
- `sum` を基本、`fit` を確認用または補助診断として併用する

運用上は、まず `method="sum"` を採用し、`method="fit"` が大きく矛盾しないことを確認する、という使い方が最も自然です。
