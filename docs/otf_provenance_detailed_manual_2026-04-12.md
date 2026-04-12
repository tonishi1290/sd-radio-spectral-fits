# OTF provenance 実装詳細説明書

## 1. この文書の目的

この文書は、`sd_radio_spectral_fits.map_3d` に実装した provenance 機能について、次の 3 つを将来参照できるようにまとめたものです。

1. **何をどこに保存するか**
2. **なぜその設計にしたか**
3. **bundle / FITS からどう読むか**

対象は、`map_3d` 側の OTF 解析パイプラインです。入力の scan table / SDFITS reader 全体の履歴管理まではまだ扱っていません。

---

## 2. 設計の基本方針

### 2.1 provenance の対象

ここでいう provenance は、**OTF 解析の下流処理履歴**です。

具体的には、たとえば次のような処理を指します。

- `grid_otf_family()` による family cube 作成
- `attach_mosaic_products()` / `attach_mosaic_products_from_mask()` による mosaic 補助量作成
- `coadd_family_cubes()` による coadd
- `mosaic_bundles()` による mosaic combine
- `plait_fft_cubes()` による PLAIT / fallback combine
- `subtract_baseline_from_bundle()` による baseline subtraction
- `make_baseline_viewer_bundle()` による baseline viewer 用 bundle 作成

### 2.2 provenance の対象外

以下は今回の provenance の正本にはしていません。

- 入力 SDFITS の `HISTORY` 自由文
- scan table reader / converter 側の完全履歴
- 巨大な ndarray 全文
- 任意の `bundle.meta` 全体

この切り分けは重要です。

- **SDFITS HISTORY** は上流の人間向け履歴
- **OTF provenance** は下流の機械可読な履歴

役割が違うので、現段階では混ぜていません。

### 2.3 runtime と永続化の分離

provenance の正本は、実行中は Python の構造化データとして持ち、FITS にはその永続化表現を書きます。

- **runtime 正本**: `bundle.meta['provenance']`
- **FITS 要約**: PRIMARY header の `PRV*` keyword
- **FITS 完全版**: `EXTNAME='PROVSTEP'` の Binary Table

この 3 層構成にしている理由は次です。

1. header に長い JSON 本文を押し込むと互換性が悪くなる
2. 実行中は dict/list の方が扱いやすい
3. 保存時は table の方が機械可読で安定

### 2.4 read/write を provenance step にしない

`write_otf_bundle()` と `read_otf_bundle()` は provenance step に入れていません。

理由は、これらは**科学処理ではなく永続化・復元**だからです。

もし step に入れると、保存・再読込だけで provenance が増殖し、処理履歴として不自然になります。

---

## 3. scan table 系との切り分け

### 3.1 現段階で扱うもの

現段階では、scan table 側については **`grid_otf_family()` の入力要約**だけを provenance に入れます。

つまり、入力 `sc_base1LU` そのものの全文を保存するのではなく、

- `nrow`
- `nchan`
- `RESTFREQ` / `RESTFRQ`
- `SPECSYS`
- `CTYPE1`
- `BUNIT`
- `TELESCOP`
- `OBJECT`
- `baseline_subtracted`
- `FDNUM` / `IFNUM` / `PLNUM` / `OBSMODE` / `SCAN` / `IS_TURN`

などを summary として保存します。

### 3.2 現段階で扱わないもの

まだ扱っていないものは次です。

- SDFITS `HISTORY` の全文読取と引き継ぎ
- scan table reader / converter 側の provenance
- 入力ファイル fingerprint / hash の確定管理
- upstream/downstream を跨いだ単一 provenance 連結

したがって、現状の provenance は

- **OTF 側の処理履歴**
- **入力 scan table の要約**

を保存するものであって、scan table 側の完全履歴を置き換えるものではありません。

### 3.3 Python 変数名は保存しない

たとえば

```python
bundle1LU = m3d.grid_otf_family(sc_base1LU, family_label='X', config=cfg)
```

の `sc_base1LU` という**Python 変数名そのもの**は保存しません。

保存するのは、`sc_base1LU` が持っていたデータの構造化 summary です。

これは論理的に自然です。関数に渡された時点で、呼び出し元のローカル変数名は一般には復元できません。

---

## 4. 保存先の仕様

## 4.1 runtime 正本: `bundle.meta`

各 `OTFBundle` には provenance 関連として少なくとも次を持ちます。

```python
bundle.meta['provenance_schema']
bundle.meta['provenance']
bundle.meta['bundle_id']
```

`bundle.meta['provenance']` は record の list です。

### 4.1.1 record の基本構造

各 record は次のキーを持ちます。

```python
{
    'step_id': str,
    'op_id': str,
    'kind': str,
    'aux': str,
    'utc': str,
    'module': str,
    'function': str,
    'input_step_ids': list[str],
    'params_input': dict | list | scalar,
    'params_config': dict | list | scalar,
    'params_resolved': dict | list | scalar,
    'results_summary': dict | list | scalar,
    'duration_sec': float | None,
    'code_version': str,
}
```

### 4.1.2 `kind` と `aux`

- `kind='main'`: 主たる科学処理 step
- `kind='aux'`: 補助的 step

現在の実装では、mosaic 補助量の付与は

- `kind='aux'`
- `aux='mosaic_products'`

として記録します。

### 4.1.3 `step_id`

`step_id` は deterministic hash です。

ハッシュ対象は次です。

- `schema_version`
- `op_id`
- `kind`
- `aux`
- `params_input`
- `params_config`
- `params_resolved`
- `input_step_ids`

`results_summary` は**ハッシュ対象に入れていません**。

理由は、結果 summary は丸め方や集計方法の見直しで揺れやすく、同じ処理の識別子まで変わると困るからです。

### 4.1.4 `bundle_id`

`bundle_id` は最後の `step_id` から導出されます。

現在の実装では概念的に

```python
bundle_id = 'b:' + last_step_id[:32]
```

です。

したがって、bundle の provenance 正本は最後の step によって識別されます。

## 4.2 FITS header の要約

PRIMARY header には provenance の要約だけを書きます。

### 4.2.1 keyword

現在の keyword は次です。

- `PRVSCHEM`: provenance schema
- `PRVNSTEP`: step 数
- `PRVMAIN`: 最後の main step の `op_id`
- `PRVAUX`: 最後の aux step の `op_id`
- `PRVUTC`: 最後の step の UTC
- `BNDLID`: bundle id
- `PRVCODE`: 最後の step の code version

### 4.2.2 目的

header は次の用途だけに使います。

- ファイルを軽く点検する
- provenance の存在をすぐ知る
- 完全版の table を読む前に要約を把握する

header は完全版の正本ではありません。

## 4.3 FITS 完全版: `PROVSTEP`

完全版は `EXTNAME='PROVSTEP'` の `BinTableHDU` に保存します。

### 4.3.1 列名

現在の列は次です。

- `STEPID`
- `OPID`
- `KIND`
- `AUX`
- `UTC`
- `MODULE`
- `FUNC`
- `INSTEP`
- `PINPUT`
- `PCONFIG`
- `PRESOLV`
- `RESULT`
- `CODEVER`
- `DSEC`

### 4.3.2 各列の意味

- `STEPID`: step 識別子
- `OPID`: 安定な処理識別子
- `KIND`: `main` または `aux`
- `AUX`: 補助 step の詳細種別
- `UTC`: step 時刻
- `MODULE`: 実装モジュール名
- `FUNC`: 実装関数名
- `INSTEP`: `input_step_ids` の JSON
- `PINPUT`: `params_input` の JSON
- `PCONFIG`: `params_config` の JSON
- `PRESOLV`: `params_resolved` の JSON
- `RESULT`: `results_summary` の JSON
- `CODEVER`: code version
- `DSEC`: `duration_sec`

### 4.3.3 JSON の方針

`INSTEP`, `PINPUT`, `PCONFIG`, `PRESOLV`, `RESULT` は JSON 文字列です。

方針は次のとおりです。

- ASCII safe
- `ensure_ascii=True`
- `allow_nan=False`
- `sort_keys=True`
- `separators=(',', ':')`

長すぎる payload は切り捨てず、例外にします。

現在の実装上限は

- `MAX_JSON_CHARS = 32760`
- `MAX_TEXT_CHARS = 1024`

です。

### 4.3.4 no-op で step を増やさない

`attach_mosaic_products()` は、既に必要な mosaic product が揃っていて `overwrite=False` かつ更新が不要な場合、何もせず return します。

この no-op では provenance step を増やしません。

これは非常に重要です。no-op 呼び出しのたびに provenance が膨らむと、履歴の意味が壊れます。

---

## 5. code version の仕様

各 step には `code_version` を入れます。

取得優先順位は次です。

1. `importlib.metadata.version('sd-radio-spectral-fits')`
2. `importlib.metadata.version('sd_radio_spectral_fits')`
3. `sd_radio_spectral_fits.__version__`
4. 失敗したら `'unknown'`

格納形式は現在

```python
sd-radio-spectral-fits=<version>
```

です。

たとえば

```python
sd-radio-spectral-fits=0.9.3
```

のように入ります。

---

## 6. `input_scantable_summary` の仕様

## 6.1 どこに入るか

`input_scantable_summary` は現在、`grid_otf_family()` の provenance の

```python
params_resolved['input_scantable_summary']
```

に入ります。

つまり、初回の family gridding step にのみ入ります。

## 6.2 なぜ `params_input` ではなく `params_resolved` か

`params_input` は、原則として**ユーザーが関数に渡した引数**です。

入力 scantable summary は、関数が受け取った Python オブジェクトを解析して作る**内部で解決された情報**なので、`params_resolved` に置いています。

## 6.3 何が入るか

概念的には次のような構造です。

```python
{
  'n_input_scantables': 1,
  'total_nrow': 1234,
  'nchan_values': [4096],
  'per_input_truncated': False,
  'per_input': [
    {
      'input_index': 0,
      'nrow': 1234,
      'nchan': 4096,
      'meta': {
        'RESTFRQ': 230538000000.0,
        'SPECSYS': 'LSRK',
        'CTYPE1': 'VRAD',
        'BUNIT': 'K',
        'TELESCOP': 'NANTEN2',
        'OBJECT': 'Orion',
        'baseline_subtracted': False,
      },
      'columns': {
        'FDNUM': {'present': True, ...},
        'IFNUM': {'present': True, ...},
        'PLNUM': {'present': True, ...},
        'OBSMODE': {'present': True, ...},
        'SCAN': {'present': True, ...},
        'IS_TURN': {'present': True, ...},
      },
    }
  ],
}
```

`sample_values` は全値ではなく、代表値の sample です。

---

## 7. どの関数がどの provenance step を追加するか

## 7.1 operation ID 一覧

現在の実装で使っている `op_id` は次です。

### 主処理 (`kind='main'`)

- `otf.grid.family.v1`
- `otf.coadd.family.v1`
- `otf.mosaic.combine.v1`
- `otf.plait.combine.v1`
- `otf.baseline.subtract.v1`
- `otf.baseline.viewer.v1`

### 補助処理 (`kind='aux'`)

- `otf.mosaic.attach_mask.v1`

## 7.2 各関数と保存内容

### `grid_otf_family()`

主 step:

- `op_id='otf.grid.family.v1'`
- `function='grid_otf_family'`
- `params_input`: family label, coordinate system, projection, velocity range, `otf_input_state` など
- `params_config`: 実効 `MapConfig`
- `params_resolved`: 実際に使われた ref coord, velocity axis, weight mode, `input_scantable_summary`, `otf_scan_summary` など
- `results_summary`: cube shape, support pixel 数, restfreq, specsys など

続いて `attach_mosaic_products()` を呼ぶので、条件に応じて aux step

- `op_id='otf.mosaic.attach_mask.v1'`

が追加されます。

### `coadd_family_cubes()`

- `op_id='otf.coadd.family.v1'`
- 入力 family bundle 群を `input_bundles` として受ける
- その後 `attach_mosaic_products()` を呼ぶ

### `mosaic_bundles()`

- `op_id='otf.mosaic.combine.v1'`
- 入力 bundle 群を `input_bundles` として受ける
- その後 `attach_mosaic_products()` を unity gain/trust で呼ぶ

### `plait_fft_cubes()`

- `op_id='otf.plait.combine.v1'`
- FFT path でも small-map fallback path でも同じ `op_id`
- `params_resolved['selected_method']` に選ばれた方法を入れる

### `subtract_baseline_from_bundle()`

主 step:

- `op_id='otf.baseline.subtract.v1'`
- linefree, ripple, baseline config を保存

`update_mosaic_products=True` で実際に更新が行われた場合のみ aux step:

- `op_id='otf.mosaic.attach_mask.v1'`

を追加します。

### `make_baseline_viewer_bundle()`

- `op_id='otf.baseline.viewer.v1'`
- baseline viewer 用に selected voxel だけ残した cube を作る

---

## 8. ancestry の仕様

## 8.1 direct input と closure

新しい step を作るとき、入力 bundle がある場合は、それらの provenance を closure としてマージします。

ただし、`input_step_ids` に入るのは **direct input の最後の step_id** です。

たとえば 2 つの family bundle から PLAIT を作るとき、

- closure としては両方の履歴を統合
- 新規 step の `input_step_ids` は `[last_x, last_y]`

となります。

## 8.2 direct input の順序は保持する

multi-input のとき、closure は dedupe しますが、`input_step_ids` の direct input 順は保持します。

したがって、

```python
plait_fft_cubes(x_bundle, y_bundle)
```

と

```python
plait_fft_cubes(y_bundle, x_bundle)
```

では、direct input の順は異なります。

これは設計上意図的です。

---

## 9. normalize の仕様

provenance に入れる値は `normalize_provenance_value()` で保守的に正規化します。

主な方針は次です。

- bool / int / float / str / None はそのまま
- `np.nan`, `np.inf`, `-np.inf` は `'nan'`, `'inf'`, `'-inf'`
- dataclass は `asdict()` 後に再帰処理
- dict は key を文字列化して sort
- ndarray は小さいものだけ list 化し、大きいものは `shape` と `dtype` だけ
- set は sort して list 化
- `fits.Header` は key-value dict 化
- `astropy.table.Table` は `colnames` と `nrow` のみ

つまり、巨大データを provenance にそのまま押し込まない設計です。

---

## 10. write / read の仕様

## 10.1 `write_otf_bundle()`

現在の既定値は

```python
write_otf_bundle(bundle, path, checksum=True, output_verify='exception')
```

です。

主な処理は次です。

1. `bundle.header` をコピー
2. `bundle.meta['provenance']` があれば header の `PRV*` を更新
3. `PROVSTEP` table を最後に追加
4. `checksum=True` で書く

### 10.1.1 予約名との関係

`bundle.table_ext` に `PROVSTEP` と同名の table があっても、provenance がある場合は内部で予約名として扱い、そちらは保存しません。

したがって、`PROVSTEP` は provenance 専用拡張名として使います。

## 10.2 `read_otf_bundle()`

`read_otf_bundle()` は、FITS 内に `PROVSTEP` があればそれを読んで

```python
bundle.meta['provenance_schema']
bundle.meta['provenance']
bundle.meta['bundle_id']
```

を復元します。

`PROVSTEP` が無い場合、provenance は空です。

### 10.2.1 重要な注意

今回永続化しているのは provenance 関連だけです。

つまり、`bundle.meta` の arbitrary な他のキーを全部 round-trip しているわけではありません。

これは**以前からそう**で、今回の provenance 実装で悪化したものではありません。

---

## 11. 直接読む方法

## 11.1 `read_otf_bundle()` 経由

もっとも標準的な方法です。

```python
import sd_radio_spectral_fits.map_3d as m3d

bundle = m3d.read_otf_bundle('out.fits')
prov = bundle.meta.get('provenance', [])
summary = bundle.meta.get('provenance_schema')
```

## 11.2 Astropy で `PROVSTEP` を直接読む

```python
from astropy.io import fits
from astropy.table import Table
import json

with fits.open('out.fits') as hdul:
    prov_tab = Table(hdul['PROVSTEP'].data)

for row in prov_tab:
    print(row['OPID'])
    params_input = json.loads(row['PINPUT'])
    params_resolved = json.loads(row['PRESOLV'])
    result = json.loads(row['RESULT'])
    print(params_input)
    print(params_resolved)
    print(result)
```

## 11.3 header だけ読む

```python
from astropy.io import fits

with fits.open('out.fits') as hdul:
    hdr = hdul[0].header
    print(hdr.get('PRVSCHEM'))
    print(hdr.get('PRVNSTEP'))
    print(hdr.get('PRVMAIN'))
    print(hdr.get('PRVAUX'))
    print(hdr.get('PRVUTC'))
    print(hdr.get('BNDLID'))
    print(hdr.get('PRVCODE'))
```

header は要約だけです。詳細は `PROVSTEP` を見てください。

---

## 12. 関数経由で読む方法

現在 `map_3d.__init__` で公開しているユーティリティは次です。

- `get_provenance()`
- `get_provenance_summary()`
- `extract_provenance_table()`
- `format_provenance()`
- `show_provenance()`

## 12.1 `get_provenance()`

source に次を渡せます。

- `OTFBundle`
- FITS path
- `fits.HDUList`

```python
import sd_radio_spectral_fits.map_3d as m3d

prov = m3d.get_provenance('out.fits')
print(len(prov))
print(prov[-1]['op_id'])
```

bundle でも同じです。

```python
prov = m3d.get_provenance(bundle)
```

## 12.2 `get_provenance_summary()`

summary だけ返します。

source に次を渡せます。

- `OTFBundle`
- FITS path
- `fits.HDUList`
- `fits.Header`

```python
summary = m3d.get_provenance_summary('out.fits')
print(summary)
```

返り値は概念的に次です。

```python
{
  'schema': 'otf-prov-v1',
  'nstep': 3,
  'last_main_op': 'otf.plait.combine.v1',
  'last_aux_op': 'otf.mosaic.attach_mask.v1',
  'last_utc': '2026-04-12T12:34:56Z',
  'bundle_id': 'b:....',
  'code_version': 'sd-radio-spectral-fits=0.9.3',
}
```

## 12.3 `extract_provenance_table()`

`astropy.table.Table` として欲しいときに使います。

```python
tab = m3d.extract_provenance_table(bundle)
print(tab.colnames)
print(tab)
```

## 12.4 `format_provenance()`

文字列化して確認したいときに使います。

```python
txt = m3d.format_provenance(bundle, include_params=True, include_results=True)
print(txt)
```

## 12.5 `show_provenance()`

標準出力にそのまま出す版です。

```python
m3d.show_provenance('out.fits', include_params=True, include_results=True)
```

---

## 13. 具体例

## 13.1 family grid の provenance を確認する

```python
import sd_radio_spectral_fits.map_3d as m3d

bundle1LU = m3d.grid_otf_family(
    sc_base1LU,
    config=cfg,
    family_label='X',
    coord_sys='galactic',
    projection='GLS',
    out_scale='TA*',
)

m3d.show_provenance(bundle1LU, include_params=True, include_results=True)
```

ここで見える主な内容は次です。

- `op_id='otf.grid.family.v1'`
- `code_version='sd-radio-spectral-fits=...'`
- `params_input` に family label, coordinate system など
- `params_config` に実効 `MapConfig`
- `params_resolved['input_scantable_summary']`
- `results_summary['cube_shape']`

## 13.2 FITS に保存してから確認する

```python
m3d.write_otf_bundle(bundle1LU, 'bundle1LU.fits', overwrite=True)

# 要約だけ
summary = m3d.get_provenance_summary('bundle1LU.fits')
print(summary)

# 詳細
m3d.show_provenance('bundle1LU.fits', include_params=True, include_results=True)
```

## 13.3 `PROVSTEP` を直接抜き出す

```python
from astropy.io import fits
from astropy.table import Table

with fits.open('bundle1LU.fits') as hdul:
    tab = Table(hdul['PROVSTEP'].data)

print(tab['OPID'])
print(tab['PRESOLV'][0])
```

## 13.4 baseline subtraction 後の provenance を見る

```python
bundle_bl = m3d.subtract_baseline_from_bundle(
    bundle1LU,
    linefree_mode='manual',
    linefree_velocity_windows_kms=['-30:-5', '25:55'],
    update_mosaic_products=True,
)

m3d.show_provenance(bundle_bl, include_params=True, include_results=True)
```

ここでは主に

- `otf.baseline.subtract.v1`
- 条件により `otf.mosaic.attach_mask.v1`

が見えます。

## 13.5 header だけで軽く点検する

```python
from astropy.io import fits

with fits.open('bundle1LU.fits') as hdul:
    hdr = hdul[0].header
    print('schema   =', hdr.get('PRVSCHEM'))
    print('nstep    =', hdr.get('PRVNSTEP'))
    print('main op  =', hdr.get('PRVMAIN'))
    print('aux op   =', hdr.get('PRVAUX'))
    print('utc      =', hdr.get('PRVUTC'))
    print('bundle id=', hdr.get('BNDLID'))
    print('code     =', hdr.get('PRVCODE'))
```

---

## 14. 将来の新規実装時のルール

新しい OTF 処理を追加するときは、以下のルールに従うのが望ましいです。

## 14.1 新しい main 処理を追加する場合

たとえば新しい combine 処理 `foo_combine()` を作るなら、出力 bundle 生成後に

```python
append_bundle_provenance_step(
    out_bundle,
    input_bundles=[...],
    op_id='otf.foo.combine.v1',
    module=__name__,
    function='foo_combine',
    kind='main',
    params_input={...},
    params_config={...},
    params_resolved={...},
    results_summary={...},
)
```

の形で追加します。

### 14.1.1 `op_id` の付け方

推奨は

```text
otf.<category>.<action>.v1
```

です。

例:

- `otf.grid.family.v1`
- `otf.mosaic.combine.v1`
- `otf.baseline.subtract.v1`

将来互換性を壊す実装変更をしたら `v2` に上げます。

## 14.2 aux 処理を追加する場合

mosaic product のような補助 step は

- `kind='aux'`
- `aux='<aux category>'`

を付けます。

main 処理と混同しないことが重要です。

## 14.3 `params_input` / `params_config` / `params_resolved` の使い分け

この 3 つを混ぜないでください。

- `params_input`: ユーザーが渡した引数
- `params_config`: dataclass/config など、処理全体の設定オブジェクト
- `params_resolved`: 実行時に解決された値、内部で決まった値、入力 summary

### 14.3.1 `MapConfig.__post_init__()` の注意

`MapConfig` は `__post_init__()` で内部補完を行うことがあります。

したがって、たとえば `cell_arcsec=None` を自動で `beam/3` に埋めた場合、後から

- ユーザーが省略したのか
- 最初から明示していたのか

は区別できないことがあります。

このような値は、必要なら `params_resolved` 側で「実際に使った値」として保存してください。

## 14.4 巨大配列を入れない

provenance は履歴です。巨大な mask や cube 全体を直接入れないでください。

必要なら

- shape
- dtype
- true 数
- support pixel 数
- representative scalar

などの summary にしてください。

## 14.5 no-op では step を増やさない

入力と出力が実質同じで、何も更新していない処理では provenance を増やさない方がよいです。

これは将来の新規実装でも同じです。

## 14.6 read/write を step にしない

永続化・復元は provenance step にしないでください。

そうしないと、科学処理履歴ではなく I/O 履歴が混ざります。

---

## 15. 現状の制限と今後の拡張候補

## 15.1 現状の制限

- scan table / SDFITS `HISTORY` はまだ取り込まない
- input file fingerprint / hash は未実装
- arbitrary `bundle.meta` 全体は round-trip しない
- `code_version` は実行環境次第で `'unknown'` になりうる
- JSON 長さ制限を超えると例外になる

## 15.2 今後の拡張候補

### A. scan table 側 provenance

将来的には scan table reader / converter 側で

- file path / fingerprint
- upstream processing state
- `HISTORY` summary

を `scantable.meta` に持たせ、OTF 側で参照する設計が考えられます。

### B. `HISTORY` の段階 2 対応

現段階では `HISTORY` を無視ではなく**未対応**としています。

もし扱うなら、OTF provenance に全文を混ぜるのではなく、

- `input_history_count`
- `input_history_hash`
- `input_history_preview`

のような structured summary として追加するのがよいです。

### C. provenance 検査ユーティリティの追加

今後、必要であれば

- `validate_provenance()`
- `diff_provenance(a, b)`
- `find_provenance_steps(source, op_id=...)`

のような補助関数を追加してもよいです。

---

## 16. 実務上の推奨運用

最後に、実際に使うときの推奨をまとめます。

1. **解析途中では `bundle.meta['provenance']` を正本として見る**
2. **保存後の外部確認では `show_provenance()` か `PROVSTEP` を使う**
3. **header の `PRV*` は quick look 専用と考える**
4. **新規処理を追加するときは `op_id` を安定に付ける**
5. **scan table 履歴と OTF 履歴は混ぜず、段階的に統合する**

この方針にしておくと、現在の実装とも整合し、将来 scan table 側の履歴管理を追加するときにも無理なく拡張できます。
