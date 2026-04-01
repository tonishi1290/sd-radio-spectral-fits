# necst_v4_sdfits_converter 実装詳細メモ

## 0. この文書の目的

この文書は、今回実装・修正した `tools/necst/necst_v4_sdfits_converter.py` のうち、特に以下の追加・変更点を後で全体マニュアルへマージしやすい形で整理した実装説明書である。

主対象は次の 4 点である。

1. `pure_rotation_model` import の堅牢化
2. `restfreq_hz` / `restfreq_ghz` の扱いの整理
3. `--vlsrk-kms-slice` の追加
4. 上記を既存の `channel_slice` / stream WCS / SDFITS writer にどう接続しているか

本書は仕様書というより、**現在の実装が何をどうしているか** を、単位・符号・優先順位・境界条件を明示して記述したものである。

---

## 1. 用語・変数・単位・座標系

### 1.1 stream

本 converter では、`[[spectrometers]]` の各ブロックを 1 本の `stream` として扱う。

各 stream は少なくとも次を持つ。

- `name`
- `fdnum`
- `ifnum`
- `plnum`
- `polariza`
- `beam`
- `frequency_axis`
- `local_oscillators`

内部では `StreamConfig` と `StreamWCS` に正規化される。

### 1.2 周波数軸

各 stream のスペクトル軸は converter 内では基本的に周波数軸として持つ。

- `CTYPE1 = 'FREQ'`
- `CUNIT1 = 'Hz'`
- `SPECSYS = 'TOPOCENT'` または `'LSRK'`
- `RESTFREQ` は Hz 単位で内部保持する

このチャットで追加した `--vlsrk-kms-slice` は、**速度指定で channel 範囲を決めるための入力** であり、converter 出力のスペクトル軸そのものを速度軸へ書き換えるわけではない。

### 1.3 速度定義

`--vlsrk-kms-slice` で最終的に比較に使う速度は、**radio 定義の LSRK 速度** である。

周波数 `nu`、基準周波数 `nu0`、光速 `c` に対して、比較に使う速度は

$$
v_{\mathrm{radio}} = c \left(1 - \frac{\nu}{\nu_0}\right)
$$

である。

ここで

- `nu` は LSRK 系へ変換後の channel 中心周波数
- `nu0` は `RESTFREQ`
- `c = 299792.458 km/s`

である。

### 1.4 相対論補正

観測所速度補正は単なる足し算ではなく、converter 実装では相対論的ドップラー係数

$$
\beta = \frac{v}{c}
$$

$$
k(v) = \sqrt{\frac{1 + \beta}{1 - \beta}}
$$

を使う。

`SPECSYS='TOPOCENT'` の stream では

$$
\nu_{\mathrm{LSRK}} = \frac{\nu_{\mathrm{TOPO}}}{k(v_{\mathrm{corr}})}
$$

として、TOPOCENT の周波数軸を LSRK 側へ 1 回だけ写している。

`SPECSYS='LSRK'` の stream では、この補正は行わない。

---

## 2. 今回の実装変更の概要

### 2.1 `pure_rotation_model` import の堅牢化

以前は `pure_rotation_model` を普通に

```python
from .multibeam_beam_measurement.pure_rotation_model import ...
```

のように import していた。そのため、親 package の `multibeam_beam_measurement/__init__.py` が壊れていると、`pure_rotation_model.py` 自体は正常でも converter 起動時に巻き添えで落ちていた。

今回の実装では `_load_pure_rotation_helpers()` を追加し、以下の順で読む。

1. 通常の relative import
2. 通常の absolute import
3. `pure_rotation_model.py` をファイルパスから直接ロード
4. それも失敗したら converter 内ローカル fallback を使用

したがって、親 package の `__init__.py` が重い・壊れている・余計な import を含む場合でも、converter は起動可能である。

### 2.2 `restfreq_ghz` / `restfreq_hz` の整理

このチャットで、RESTFREQ の指定を次のように整理した。

- CLI では `--restfreq-ghz`
- config では `restfreq_hz` と `restfreq_ghz` の両方を受理
- 内部では最終的に `restfreq_hz` に正規化して `StreamWCS.restfreq_hz` へ渡す

### 2.3 `--vlsrk-kms-slice` の追加

channel index ではなく VLSRK 速度範囲で切り出したい、という要望に対応するため、converter で

```text
--vlsrk-kms-slice "[vmin,vmax]"
```

を追加した。

ただし row ごとに厳密に再計算するのではなく、**各 stream の最初の ON row を基準に 1 回だけ channel 範囲へ変換** し、その結果を既存の `channel_slice` 経路へ流す実装にしている。

この設計にした理由は次である。

- OTF 中の観測所速度変化は通常小さい
- row ごとに完全再計算すると重い
- 既存の `channel_slice_bounds` / WCS 更新 / spectral extract の流れを再利用できる

---

## 3. import fallback 実装

### 3.1 目的

`pure_rotation_model.py` にある次の 2 関数だけが欲しい。

- `rotate_el0_offset`
- `offset_xy_to_beam_azel`

converter は mapping パッケージ全体を使いたいわけではない。beam-center Az/El を計算するために、この 2 関数を使うだけである。

### 3.2 fallback の段階

実装は `_load_pure_rotation_helpers()` にまとまっている。

1. relative import
2. absolute import
3. `importlib.util.spec_from_file_location()` による直接ロード
4. converter 内 fallback 実装

### 3.3 converter 内 fallback の式

fallback の回転は、El=0 で定義されたオフセット `(dx0, dy0)` を elevation に応じて回す簡易式である。

回転角を

$$
\theta = s \cdot \mathrm{El}
$$

とし、`s = +/- 1` を `rotation_sign` とすると、回転後オフセットは

$$
dx' = dx_0 \cos\theta - dy_0 \sin\theta
$$

$$
dy' = dx_0 \sin\theta + dy_0 \cos\theta
$$

である。

その後、beam-center Az/El は

$$
\Delta \mathrm{Az}[\mathrm{deg}] = \frac{dx'[\mathrm{arcsec}]}{3600 \cos(\mathrm{El})}
$$

$$
\Delta \mathrm{El}[\mathrm{deg}] = \frac{dy'[\mathrm{arcsec}]}{3600}
$$

$$
\mathrm{Az}_{\mathrm{beam}} = \mathrm{Az}_{\mathrm{boresight}} + \Delta \mathrm{Az}
$$

$$
\mathrm{El}_{\mathrm{beam}} = \mathrm{El}_{\mathrm{boresight}} + \Delta \mathrm{El}
$$

として作る。

高仰角では `cos(El)` が小さくなるので、実装では `|cos(El)|` が十分小さい領域を避ける安全策を入れている。

---

## 4. RESTFREQ の実装

## 4.1 入力形式

config の `frequency_axis` では次の両方を受ける。

- `restfreq_hz`
- `restfreq_ghz`

例:

```toml
[spectrometers.frequency_axis]
nchan = 32768
definition_mode = "first_center_and_delta"
first_channel_center_hz = 1.500000000e9
channel_spacing_hz = -6.103515625e4
restfreq_ghz = 115.2712018
ctype1 = "FREQ"
cunit1 = "Hz"
specsys = "TOPOCENT"
veldef = "RADIO"
```

### 4.2 正規化規則

内部では `_canonicalize_restfreq_in_frequency_axis()` で正規化する。

- `restfreq_ghz` のみある場合
  - `restfreq_hz = restfreq_ghz * 1e9`
- `restfreq_hz` のみある場合
  - そのまま採用
- 両方ある場合
  - 数値整合を確認
  - 不整合ならエラー

したがって、最終的な内部表現は常に Hz である。

### 4.3 CLI 上書き

CLI の `--restfreq-ghz` が指定された場合、`apply_restfreq_ghz_override()` により、**選択された全 stream の `frequency_axis.restfreq_hz` を同じ値で上書き** する。

このとき config 内の `restfreq_ghz` は削除され、内部辞書は `restfreq_hz` に統一される。

優先順位は

1. CLI `--restfreq-ghz`
2. config `restfreq_hz`
3. config `restfreq_ghz`

である。

### 4.4 legacy single-stream の扱い

`--spectrometer-config` を使わない legacy 単一 stream モードでは、`build_legacy_single_stream_config(args)` が自動設定を作る。

このとき `frequency_axis.restfreq_hz` は

$$
\mathrm{RESTFREQ}[\mathrm{Hz}] = \mathrm{restfreq\_ghz} \times 10^9
$$

で入る。

---

## 5. `--vlsrk-kms-slice` の実装

## 5.1 目的

既存の `--channel-slice` は channel index 指定であり、例えば

```text
--channel-slice "[1000,5000)"
```

のように使う。

今回追加した `--vlsrk-kms-slice` は、これを VLSRK 速度指定で行うためのものである。

例:

```text
--vlsrk-kms-slice "[-20,80]"
```

このオプションは、**converter 実行前半で各 stream ごとに channel 範囲へ変換され、その後は既存の `channel_slice` と同じ経路で処理される**。

## 5.2 入力形式

`_parse_vlsrk_kms_slice_spec()` が受ける形式は次である。

- `"[vmin,vmax]"`
- `"[vmin,vmax)"`
- `"(vmin,vmax]"`
- `"(vmin,vmax)"`
- Python 的な長さ 2 の list / tuple

実装上、文字列パースは指数表記も受ける。

例:

```text
"[-30,50]"
"(1.0e1,8.0e1]"
```

### 5.3 大小反転の扱い

速度範囲は入力で `vmin > vmax` でもよい。

例えば

```text
--vlsrk-kms-slice "[80,-20]"
```

のように与えた場合でも、内部では low/high に正規化する。

ただし、区間端の inclusive / exclusive は**反転前の括弧の意味を保ったまま** 正規化される。

つまり、`[a,b)` を逆順に書いたときでも、実装は low/high 側へ inclusive 情報を正しく移す。

---

## 6. `--vlsrk-kms-slice` から `channel_slice_bounds` を作る流れ

## 6.1 全体の流れ

main の流れは概略的に次である。

1. config を読み、stream を確定する
2. `apply_restfreq_ghz_override()` を適用する
3. `apply_channel_slice_config()` を適用する
4. `--vlsrk-kms-slice` があれば、各 stream で `_compute_stream_vlsrk_channel_slice()` を呼ぶ
5. stream ごとの `channel_slice_bounds` と `stream.wcs` を更新する
6. 以後の spectral 読み出しでは、その bounds を使って切り出す

## 6.2 排他制御

`--vlsrk-kms-slice` は、既存の channel slice 系と同時指定できない。

具体的には次と排他である。

- CLI `--channel-slice`
- config `global.channel_slice`
- config `stream.channel_slice`

理由は、速度範囲から決めた slice と channel 指定の slice が同時に存在すると、どちらを優先するかで曖昧さが生じるためである。

## 6.3 代表 row の選び方

`_compute_stream_vlsrk_channel_slice()` は、各 stream の分光テーブルを最初に 1 回だけ読んで、代表 row を選ぶ。

選び方は次である。

1. spectral table の時刻 `t_spec_unsorted` を得る
2. `position` / `obsmode` / `mode` 等から OBSMODE 相当のラベル列を抽出する
3. 時刻順に並べる
4. 最初の `ON` row を探す
5. `ON` が無ければ、最初の有限時刻 row を fallback とする

この代表 row についてのみ、boresight → beam-center → RA/Dec → site velocity correction を計算し、その結果を stream 全体の代表値として使う。

### 6.4 なぜ row ごとにやらないか

ここでの目的は、厳密な Doppler 補正付き再格子化ではなく、**converter の前段で大まかな帯域切り出しをすること** である。

したがって

- 1 row ごとに `v_corr` を再計算するより軽い
- OTF 中の時間変化は通常小さい
- 切り出し範囲さえ概ね合っていればよい

という理由から、1 stream につき 1 回の代表計算にしている。

---

## 7. 代表 row での VLSRK 補正計算

## 7.1 必要な量

代表 row について必要なのは次である。

- 時刻 `t_ref` [unix s]
- その時刻の beam-center Az/El
- そこから得た RA/Dec
- 観測所 site 情報

### 7.2 観測所速度補正

`_calc_vlsrk_correction_kms_for_site()` は

- `Site.to_earthlocation()`
- `Time(unix, scale='utc')`
- `GCRS -> LSRK`

を使い、観測所速度ベクトルを求めた後、目標方向単位ベクトルへの射影で視線方向速度を得る。

式としては、RA, Dec を使って

$$
\hat{n} = (\cos\delta\cos\alpha,\ \cos\delta\sin\alpha,\ \sin\delta)
$$

とし、LSRK での観測所速度ベクトルを

$$
\vec{v}_{\mathrm{obs}} = (v_x, v_y, v_z)
$$

とすると、実装が使う補正は

$$
v_{\mathrm{corr}} = \vec{v}_{\mathrm{obs}} \cdot \hat{n}
$$

である。

単位は km/s である。

### 7.3 符号の扱い

converter 実装では、TOPOCENT の周波数軸を LSRK 側へ動かすときに

$$
\nu_{\mathrm{LSRK}} = \frac{\nu_{\mathrm{TOPO}}}{k(v_{\mathrm{corr}})}
$$

を使う。

ここで

$$
k(v) = \sqrt{\frac{1 + v/c}{1 - v/c}}
$$

である。

したがって、`v_corr > 0` なら `k > 1` となり、`nu_lsrk < nu_topo` になる。

この符号規約は、converter 内部の「TOPOCENT の channel 軸を LSRK に写してから radio 速度を比較する」という実装に対応している。

---

## 8. LSRK 速度から channel index を決める式

## 8.1 入力

stream WCS `sw` に対して

- `rest_hz = sw.restfreq_hz`
- `freq_axis_hz = _channel_centers_from_wcs(sw)`

を得る。

ここで `sw.ctype1` は `FREQ` でなければならない。

### 8.2 TOPOCENT / LSRK の分岐

- `SPECSYS='TOPOCENT'`
  - 相対論補正で `freq_lsrk_hz` を作る
- `SPECSYS='LSRK'`
  - `freq_lsrk_hz = freq_axis_hz`
- それ以外
  - エラー

### 8.3 radio 速度へ変換

比較用の LSRK radio 速度は

$$
v_{\mathrm{LSRK,radio}} = c \left(1 - \frac{\nu_{\mathrm{LSRK}}}{\nu_0}\right)
$$

である。

ここで

- `nu_LSRK = freq_lsrk_hz`
- `nu0 = rest_hz`

である。

### 8.4 速度範囲から mask を作る

`_mask_from_vlsrk_bounds()` は

- low 側 inclusive / exclusive
- high 側 inclusive / exclusive
- finite check

をまとめて論理積で判定する。

### 8.5 channel index の決め方

mask に入った channel の index を `idx = flatnonzero(mask)` とし、最終的な bounds は

- `start = idx[0]`
- `stop = idx[-1] + 1`

とする。

ここで重要なのは、**軸が昇順か降順かを仮定していない** 点である。

つまり、channel と velocity の増減関係が反転していても、mask で拾った index の最初と最後から作るので安全である。

これは、このチャットで特に注意した点である。

---

## 9. `channel_slice_bounds` 適用後の WCS 更新

## 9.1 方針

一度 `channel_slice_bounds = (start, stop)` が決まったら、既存の `channel_slice` と同じく `_slice_stream_wcs()` で stream の WCS を切り詰める。

### 9.2 更新内容

切り詰め後の `StreamWCS` は

- `nchan = stop - start`
- `crval1_hz = sliced channel centers の先頭`
- `cdelt1_hz = 元の値`
- `crpix1 = 1.0`

になる。

### 9.3 一様軸チェック

切り出し後の中心周波数列が一様でない場合はエラーにする。

これは converter が、切り出した後もなお `SpectralAxisUniform` として扱えることを前提にしているためである。

---

## 10. spectral 配列の切り出し実装

## 10.1 どこで切っているか

実際の 2D スペクトル配列は `_extract_spectral_from_structured()` で読み出すときに `channel_slice_bounds` を渡している。

つまり、full-NCHAN を先に全部作ってから後段で切るのではなく、**可能な限り早い段階で配列を切ってメモリ節約する** 実装である。

### 10.2 object 配列の場合

分光データが object 配列の場合は row ごとに

- full row を 1D 配列へ変換
- `start:stop` で切る
- `float32` へ格納

する。

### 10.3 2D ndarray の場合

すでに `(nrow, nchan)` の 2D 配列であれば、そのまま

```python
arr2[:, start:stop]
```

で切る。

### 10.4 整合チェック

切り出し後の shape は `stream.wcs.nchan` と一致しなければならず、不一致ならエラーになる。

---

## 11. 時刻抽出の実装

## 11.1 目的

spectral table の時刻は、装置や環境によって

- `timestamp` 列の文字列
- `time` 列の数値
- `time_spectrometer`
- `unix_time` / `unixtime`

など、複数形式があり得る。

converter では `_select_spectral_time_from_structured()` でこれを正規化する。

### 11.2 timestamp suffix

文字列 timestamp の末尾 suffix を見て、以下を扱う。

- `UTC`
- `GPS`
- `TAI`
- `PC`

方針は次である。

- `UTC`, `GPS`, `TAI`
  - 変換して unix 秒へ
- `PC`
  - 信頼しないので数値時刻へ fallback
- 空 / NaN / 不明 suffix
  - 数値時刻へ fallback

### 11.3 参照用情報の記録

どの時刻列を使ったかは `time_meta` に記録し、最終的に writer history にも反映される。

---

## 12. site 情報の扱い

converter は site 情報を

- CLI 指定
- RawData 配下の config TOML から自動取得
- 既定値

の優先順位で決める。

`--vlsrk-kms-slice` は site velocity correction を使うため、site 情報が必要である。

site が正しく取れないと、代表 row の `v_corr_kms` が不正になるので注意が必要である。

---

## 13. beam-center Az/El, boresight, RA/Dec の関係

この converter では、Az/El 系の意味を分けて扱う。

- `AZ_ENC/EL_ENC`
  - encoder の実測 Az/El を spectral 時刻へ内挿したもの
- `CORR_AZ/CORR_EL`
  - dlon/dlat
- `BORE_AZ/BORE_EL`
  - corrected boresight Az/El
- `AZIMUTH/ELEVATIO`
  - 最終的な beam-center Az/El

`--vlsrk-kms-slice` の代表 row 計算で最終的に使うのは、beam model 適用後の beam-center Az/El から作った RA/Dec である。

このため、beam offset / beam rotation の定義が変われば、代表 row の `v_corr_kms` もわずかに変わり得る。

---

## 14. FDNUM / IFNUM / PLNUM の扱い

各 stream の `fdnum`, `ifnum`, `plnum` は `StreamConfig` に保存され、各 row を writer に渡すときにも

- `fdnum`
- `ifnum`
- `plnum`

として明示的に渡される。

したがって、multi-stream config を使い `--stream-name` で 1 本だけ選んだ場合も、legacy 単一 stream ではなく、その stream に設定された `fdnum/ifnum/plnum` が使われる。

legacy single-stream モードでのみ既定値 `0,0,0` となる。

---

## 15. CLI と config の優先順位

このチャットで関係する主な優先順位をまとめる。

### 15.1 RESTFREQ

1. CLI `--restfreq-ghz`
2. stream config `restfreq_hz`
3. stream config `restfreq_ghz`

### 15.2 channel slice

1. `--vlsrk-kms-slice` がある場合
   - channel slice 系は同時指定不可
2. CLI `--channel-slice`
3. stream config `channel_slice`
4. global config `channel_slice`
5. slice なし

### 15.3 stream selection

1. `--stream-name` で明示指定されたもの
2. それ以外は `enabled` / `use_for_convert` に従う

---

## 16. 異常系

## 16.1 `--vlsrk-kms-slice` 関連

次の場合はエラーで止める。

- `restfreq_hz <= 0`
- `CTYPE1 != 'FREQ'`
- `SPECSYS` が `TOPOCENT` / `LSRK` 以外
- 指定速度範囲に 1 本も channel が入らない
- 計算された速度軸が finite でない
- `--channel-slice` と同時指定
- config の `channel_slice` と同時指定

### 16.2 beam / pointing 関連

次の場合も明示的に止める。

- beam-center Az/El に NaN が出る
- RA/DEC 用に選んだ Az/El source に NaN が出る
- GLON/GLAT 変換後に NaN が出る

### 16.3 RESTFREQ 関連

- `restfreq_ghz` が有限正値でない
- `restfreq_hz` と `restfreq_ghz` が両方あり不整合

---

## 17. 最小例

### 17.1 config で `restfreq_ghz` を使う例

```toml
schema_version = 1

[global]
telescope = "NANTEN2"
db_namespace = "necst"

[[spectrometers]]
name = "2RU"
fdnum = 2
ifnum = 0
plnum = 0
polariza = "XX"

[spectrometers.frequency_axis]
nchan = 32768
definition_mode = "first_center_and_delta"
first_channel_center_hz = 2.500000000e9
channel_spacing_hz = -7.62939453125e4
restfreq_ghz = 230.5380000
ctype1 = "FREQ"
cunit1 = "Hz"
specsys = "TOPOCENT"
veldef = "RADIO"
```

### 17.2 CLI で全選択 stream の RESTFREQ を上書きする例

```bash
necst_v4_sdfits_converter \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --stream-name 2RU \
  --restfreq-ghz 230.5380000 \
  necst_otf_20260321_201209_orion-kl \
  --out orion-kl-2RU.fits
```

### 17.3 VLSRK 速度範囲で切る例

```bash
necst_v4_sdfits_converter \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --stream-name 2RU \
  --restfreq-ghz 230.5380000 \
  --vlsrk-kms-slice "[-30,60]" \
  necst_otf_20260321_201209_orion-kl \
  --out orion-kl-2RU.fits
```

### 17.4 逆順入力も受理される例

```bash
necst_v4_sdfits_converter \
  --spectrometer-config 12co-NANTEN2-multi_260331TO.conf \
  --stream-name 2RU \
  --vlsrk-kms-slice "[60,-30]" \
  necst_otf_20260321_201209_orion-kl
```

内部では low/high に正規化される。

---

## 18. 実装上の注意

### 18.1 `--vlsrk-kms-slice` は「厳密再格子化」ではない

これは converter の前段での大まかな帯域切り出しである。row ごとに厳密な LSRK 再計算をしているわけではない。

### 18.2 出力軸は速度軸へ変えない

`--vlsrk-kms-slice` は channel 範囲の決定にのみ使う。出力の WCS は依然として周波数軸である。

### 18.3 `SPECSYS='TOPOCENT'` / `'LSRK'` 以外は未対応

今回の実装では、`BARYCENT` や他の spectral frame は扱わない。

### 18.4 `channel_slice` は既存経路を再利用している

今回の実装は、新しいスペクトル切り出し器を別に作ったのではなく、既存の `channel_slice_bounds` の経路へ落としている。このため既存挙動との一貫性が高い。

---

## 19. writer history に残る情報

今回の実装では、次のような情報が history へ残る。

- `restfreq_ghz_cli`
- `channel_slice_cli` または `channel_slice_global`
- `vlsrk_kms_slice_cli`
- stream ごとの `vlsrk_slice_i`
  - stream 名
  - label
  - bounds
  - 代表 row の mode
  - 代表時刻 unix
  - `vcorr_kms`
- beam model の要約
- stream ごとの `fdnum/ifnum/plnum`
- spec time basis / suffix / fallback

したがって、後で「どの stream にどの slice が適用されたか」は FITS 側でもかなり追跡しやすい。

---

## 20. 今後の拡張候補

### 20.1 row ごとの `v_corr` 再計算

現在は 1 stream につき 1 回の代表計算である。必要なら、将来は row ごとの Doppler 補正へ拡張できる。

ただし、その場合は

- 計算コスト
- channel 範囲が row ごとに変わる問題
- 出力 row を同一 NCHAN に保つ方法

を別途設計しなければならない。

### 20.2 他 frame 対応

現状は `TOPOCENT` / `LSRK` のみである。必要なら将来

- `BARYCENT`
- `HELIOCEN`
- 他の frame

も追加できるが、符号規約と周波数→速度変換の責務分離を明確にして進める必要がある。

---

## 21. まとめ

このチャットで入れた実装の本質は、次の 3 点である。

1. **converter 起動を package `__init__.py` に依存させない**
2. **RESTFREQ を Hz/GHz の両方から一貫した内部表現へ正規化する**
3. **VLSRK 速度範囲を各 stream で 1 回だけ channel 範囲へ変換し、既存の `channel_slice` 経路へ落とす**

この設計により、処理系を複雑に増やさずに

- import 問題の回避
- `restfreq_ghz` 対応
- `--vlsrk-kms-slice` 対応
- channel/velocity 反転への安全性

を実現している。

