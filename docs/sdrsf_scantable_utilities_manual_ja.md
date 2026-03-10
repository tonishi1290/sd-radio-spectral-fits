# Scantable Utilities マニュアル

この文書は `scantable_utils.py` と `tempscale.py` の `set_beameff()` を中心に、Scantable の確認・抽出・補助更新の使い方をまとめた補助説明書です。  
対象は添付ソースの現実装です。

## 1. 何のためのモジュールか

Scantable Utilities は、SDFITS / Scantable を処理するときの前処理・点検・部分抽出のための実務モジュールです。  
主に次の用途で役立ちます。

- table にどの列があるか調べる
- `RESTFRQ`, `TEMPSCAL`, `BEAMEFF`, `VELOSYS` などを確認する
- scan や観測モードで一部行だけ切り出す
- マッピング用オフセットを計算する
- 複数ファイルを `SCAN` 重複なく結合する
- `FDNUM`, `IFNUM`, `PLNUM` ごとに `BEAMEFF` を与える

## 2. インポート

添付 ZIP には package root の `__init__.py` が含まれていないため、トップレベル再公開の最終形はこの資料だけでは断定できません。  
したがって、以下のモジュール経由インポートが最も確実です。

```python
import sd_radio_spectral_fits.scantable_utils as su
import sd_radio_spectral_fits.tempscale as ts
```

ユーザー環境で
```python
import sd_radio_spectral_fits as sd
```
として `sd.set_beameff(...)` が使える場合は、その再公開を利用して構いません。ただし、添付ソース上で直接確認できたのは `scantable_utils.py` / `tempscale.py` 側の実装です。

## 3. 関数一覧

```python
describe_columns(sc)
show_scantable(inputs, rows=None, columns="default", head=20, show_legend=False,
               extra_data=None, ref_coord=None, frame="ICRS", projection="GLS", unit="arcsec")
calc_mapping_offsets(sc, ref_coord=None, frame="ICRS", projection="GLS", unit="arcsec",
                     cos_mode="point", verbose=True)
merge_scantables(inputs, sort_by_time=False, shift_scan_id=True)
update_metadata(sc, column, value, rows=None, force=False, verbose=True)
set_beameff(sc, efficiency, rows=None, verbose=True)
find_scans(sc, query=None, extra_data=None, **kwargs)
filter_scantable(sc, query=None, extra_data=None, rows=None, **kwargs)
select_rows(sc, rows)
```

## 4. 実務でのおすすめ順序

1. `describe_columns()` / `show_scantable()` で列確認
2. `find_scans()` / `filter_scantable()` / `select_rows()` で試験サブセット作成
3. 必要に応じて `calc_mapping_offsets()`
4. 複数日なら `merge_scantables()`
5. `set_beameff()` でビーム効率付与
6. どうしても必要なら `update_metadata()`

## 5. `show_scantable()`

### 何が便利か
- table に無い列を meta や alias から拾う
- `TIMESTAMP` を補う
- offset 計算結果や外部 DataFrame をその場で表示に足せる

### 例
```python
su.show_scantable(
    sc,
    rows="0:10",
    columns="SCAN, TIMESTAMP, OBSMODE, RESTFRQ, RA, DEC, TEMPSCAL, BEAMEFF"
)
```

## 6. `calc_mapping_offsets()`

### 現コードで確実に使える投影
- `GLS` / `SFL`
- `CAR` / `NONE`

### 重要
- `TAN` は docstring には出てくるが現コードでは未実装
- `ref_coord="Orion KL"` のような名前解決も実装していない
- 確実なのは座標文字列または `SkyCoord`

### 例
```python
ofs = su.calc_mapping_offsets(
    sc,
    ref_coord="05h35m14.5s -05d22m30s",
    frame="ICRS",
    projection="GLS",
    unit="arcsec",
    cos_mode="ref",
)
```

## 7. `find_scans()` / `filter_scantable()`

### kwargs で簡単に探す
```python
idx = su.find_scans(sc, OBSMODE="ON")
sc_on = su.filter_scantable(sc, OBSMODE="ON")
```

### query で複雑条件
```python
idx = su.find_scans(sc, query="SCAN > 10 and OBSMODE == 'ON'")
```

### offset を使って中心部だけ抜く
```python
ofs = su.calc_mapping_offsets(sc, ref_coord="05h35m14.5s -05d22m30s")
sc_center = su.filter_scantable(
    sc,
    query="OFS_LON**2 + OFS_LAT**2 < 60**2",
    extra_data=ofs,
)
```

## 8. `merge_scantables()`

### 重要ポイント
- `shift_scan_id=True` が既定
- これによりファイル間の `SCAN` 衝突を防ぐ
- `group_mode="scan"` を使う前には特に重要

```python
sc_merged = su.merge_scantables(
    ["day1.fits", "day2.fits"],
    sort_by_time=True,
    shift_scan_id=True,
)
```

## 9. `update_metadata()`

### 使いどころ
- `OBJECT` などのラベル修正
- 補助列追加
- 限定的な `RESTFRQ` 修正

### 注意
- `RESTFRQ`, `RESTFREQ`, `CTYPE1`, `CRVAL1` などは危険列
- `force=True` が必要なことがある
- 破壊的変更

## 10. `set_beameff()`

### 何をするか
- `BEAMEFF` 列を作成または更新する
- data は変えない
- `TEMPSCAL` も変えない

### 何をしないか
- `TA* -> TR*` 変換
- `TEMPSCAL='TR*'` への変更
- `target_scale` 指定

### 全行へ一括設定
```python
su.set_beameff(sc, efficiency=0.42)
```

### 一部行だけ
```python
su.set_beameff(sc, efficiency=0.42, rows="0:20")
```

### FDNUM / IFNUM / PLNUM ごとに設定
```python
beam_eff_map = {
    (0, 0, 0): 0.46,
    (0, 0, 1): 0.44,
    (1, 0, 0): 0.42,
    (1, 0, 1): 0.41,
}

for (fdnum, ifnum, plnum), eta in beam_eff_map.items():
    idx = su.find_scans(sc, FDNUM=fdnum, IFNUM=ifnum, PLNUM=plnum)
    su.set_beameff(sc, efficiency=eta, rows=idx)
```

### 結果確認
```python
su.show_scantable(
    sc,
    rows="0:20",
    columns="FDNUM, IFNUM, PLNUM, TEMPSCAL, BEAMEFF"
)
```

## 11. viewer / coadd との関係

- viewer は `BEAMEFF` を参照して TR* 表示へ切り替えられる
- `coadd.run_velocity_coadd(..., out_scale="TR*")` では `BEAMEFF` の整合が重要
- `BEAMEFF` が行ごとに混在している場合、coadd の `normalize_if_mixed` 設計が効いてくる

## 12. 現場での最小レシピ

```python
import sd_radio_spectral_fits.scantable_utils as su

# 1. まず確認
su.describe_columns(sc)
su.show_scantable(sc, rows="0:5", columns="SCAN, OBSMODE, FDNUM, IFNUM, PLNUM, TEMPSCAL, BEAMEFF")

# 2. ビームごとに効率を付与
idx = su.find_scans(sc, FDNUM=0, IFNUM=1, PLNUM=0)
su.set_beameff(sc, efficiency=0.45, rows=idx)

# 3. 確認
su.show_scantable(sc, rows="0:20", columns="FDNUM, IFNUM, PLNUM, TEMPSCAL, BEAMEFF")
```
