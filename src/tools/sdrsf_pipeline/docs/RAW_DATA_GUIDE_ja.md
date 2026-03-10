# 生スペクトルデータ作成ガイド (HOT/ON/OFF)

本パイプラインは、まず **dumpごとの生スペクトル** から開始します。

- `HOT` : **常温ロード (Tamb)** を見ているスペクトル
- `ON`  : **天体/空** を見ているスペクトル
- `OFF` : **オフ点 (空)** を見ているスペクトル

いずれも **(Ndump, Nchan)** の2次元データを想定します。
推奨は pandas `DataFrame` で、行インデックスを **UTCの `DatetimeIndex`** にします。

## 1. 最低限必要な情報

### 1.1 HOT/ON/OFF スペクトル表

- 形状: `(Ndump, Nchan)`
- インデックス: dump時刻 (UTC)。`DatetimeIndex` 推奨
- 値: 分光計の生出力 (counts / power 等、線形量を想定)
- チャンネル並び: HOT/ON/OFF で一致していること
  - チャンネルが増えると周波数が下がる場合でも問題ありません。その場合は `CDELT1` を負にします。

### 1.2 meta.json (スペクトルWCS + 基本情報)

`meta.json` はスペクトル軸を定義します。最低限:

- `rest_hz` : 静止周波数 (Hz)
- `crval1_hz` : 参照周波数 (Hz)
- `crpix1` : 参照ピクセル (FITS慣習で1始まり)
- `cdelt1_hz` : チャンネル幅 (Hz/ch)。負になり得ます
- `ctype1` : 通常 `"FREQ"`
- `cunit1` : 通常 `"Hz"`
- `timesys` : `"UTC"`
- `specsys` : `"LSRK"` (ラベル。軸そのものが既にLSRKであることは別問題です)

任意ですが **強く推奨** (速度補正の再現性のため):

- `site_lat_deg`, `site_lon_deg`, `site_height_m` (または `site_name`)

## 2. マッピングテーブル (dumpごと)

座標やスキャン情報がある場合は mapping を与えます。
推奨列 (ON dump と同じ行数・時刻に揃える):

- `timestamp` (UTC) もしくは index が `DatetimeIndex`
- `ra_deg`, `dec_deg` (ICRS/J2000 の度)
- `coord_frame` (文字列): `"icrs"`, `"fk5"`, `"galactic"` 等
- `pos_id` (任意): 位置を表す整数 (同一点の合算などに利用)
- その他の列 (az/el, scan_id 等) があっても構いません

パイプラインは `coord_frame` が想定と一致するか検証し、違えば明示的にエラーにします。

## 3. rawspec コンテナの生成

次段の入力となる rawspec pickle を作ります:

```bash
sdrsf-make-rawspec --hot hot.pkl --on on.pkl --off off.pkl --meta meta.json --mapping mapping.csv --out rawspec.pkl
```

Doppler補正や較正の **前** にデータ量を落とす場合:

```bash
sdrsf-make-rawspec ... --max-dumps 2000 --ch-start 1024 --ch-stop 3072
```

HOT/ON/OFF を同一条件で切り出し、WCSも整合するように (CRVAL1 のシフト、NCHAN更新) 自動で更新します。

## 4. 典型的な落とし穴

- **時刻はUTC**: ローカル時刻が混ざると、OFF内挿や速度補正が破綻します
- **Nchan不一致**: HOT/ON/OFF は切り出し後も必ず一致
- **WCS符号ミス**: 周波数がチャンネルと逆向きなら `cdelt1_hz` を負にする
- **サイト情報不足**: VLSRK 補正には site (lat/lon/height か astropy が認識する site_name) が必須
- **座標系の食い違い**: mapping が Galactic なのに RA/Dec を入れる等は、早めにエラーにして気づける方が安全です

