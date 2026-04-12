# Scantable / SDFITS provenance 実装と確認方法

## 1. この文書の対象

この文書は、`sd_radio_spectral_fits` における **Scantable 系列 / SDFITS 系列の provenance 実装**を、現在の公開 API と実装方針に合わせて整理したものです。

ここでいう provenance は、主に次の 3 種類です。

1. **ソフトウェア provenance**
   - どの解析ソフト・どの版がこの `Scantable` / FITS を生成・更新したか
2. **行ごとの provenance**
   - 各 row / dump がどの frontend / backend / beam / IF / polarization / sideband 由来か
3. **処理履歴 provenance**
   - どの処理をどの条件で実行したか

本書では、次をまとめます。

- `Scantable` のどこに何を保存しているか
- SDFITS にどう写像しているか
- 中身の確認方法
- `summarize_provenance()` / `print_provenance()` の使い方
- FITS round-trip 時の制約
- Scantable と FITS の間の非対称性
- 将来整理した方がよい点

---

## 2. 関係する公開 API

主に次の関数・クラスが関係します。

- `fitsio.Scantable`
- `fitsio.read_scantable(path)`
- `fitsio.write_scantable(path, scantable, ...)`
- `fitsio.write_sdfits(out_path, meta, data, table, history, ...)`
- `fitsio.stamp_scantable_code_provenance(scantable, stage=...)`
- `scantable_utils.summarize_provenance(sc_or_path)`
- `scantable_utils.print_provenance(sc_or_path, ...)`

トップレベル import の想定例:

```python
import sd_radio_spectral_fits as sdrsf

sc = sdrsf.read_scantable("example.fits")
prov = sdrsf.summarize_provenance(sc)
sdrsf.print_provenance(sc)
```

---

## 3. まず定義: `Scantable` の provenance 保存場所

`Scantable` は概念的には次の 4 要素を持ちます。

```python
Scantable(
    meta: dict,
    data: np.ndarray | list[np.ndarray],
    table: pandas.DataFrame,
    history: dict,
)
```

provenance の観点では、主役は `meta`, `table`, `history` の 3 つです。

### 3.1 `meta`

**ファイル全体 / Scantable 全体に対する属性**を持ちます。

典型例:

- `SWNAME`, `SWVER`
- `OBJECT`, `TELESCOP`, `INSTRUME`
- `SITELAT`, `SITELONG`, `SITEELEV`
- `OBSGEO-X`, `OBSGEO-Y`, `OBSGEO-Z`
- `RESTFRQ`, `RESTFREQ`
- `SPECSYS`, `SSYSOBS`
- `CTYPE1`, `CRVAL1`, `CDELT1`, `CRPIX1`
- `TEMPSCAL`

### 3.2 `table`

**各 row / dump ごとの provenance** を持ちます。

典型例:

- `FRONTEND`, `BACKEND`, `SAMPLER`
- `FDNUM`, `IFNUM`, `PLNUM`, `POLARIZA`
- `OBSFREQ`, `LO1FREQ`, `LO2FREQ`, `LO3FREQ`
- `SIDEBAND`, `SB1`, `SB2`, `SB3`
- `TEMPSCAL`, `CALSTAT`, `SPECSYS`, `SSYSOBS`
- `SCAN`, `SUBSCAN`, `INTGRP`, `OBSMODE`, `OBJECT`
- `FLAGROW`

### 3.3 `history`

**処理履歴 provenance** を持ちます。

典型例:

- `code_software`, `code_version`, `code_stage`
- `task`, `input`, `created_at_utc`
- `baseline_history`
- `scale_history`
- `velocity_regrid`
- `qc_summary`
- `notes`
- `rows`, `exclude_rows`
- `rest_freq`, `v_corr_col`
- `fitsio_migration`

---

## 4. provenance の役割分担

現在の実装は、provenance を 1 か所に押し込まず、意味に応じて分けています。

### 4.1 ソフトウェア provenance

保存先:

- `sc.meta["SWNAME"]`
- `sc.meta["SWVER"]`
- `sc.history["code_software"]`
- `sc.history["code_version"]`
- `sc.history["code_stage"]` 必要時のみ

意味:

- `SWNAME`, `SWVER`
  - header 的に見えるソフトウェア識別子
- `code_software`, `code_version`
  - history / provenance 空間で見えるソフトウェア識別子
- `code_stage`
  - どの処理段階で現在の `Scantable` が生成・更新されたか

### 4.2 行ごとの provenance

beam / IF / 偏波 / backend / frontend のような **row 単位の由来** は `history` ではなく `table` の列に置きます。

これは次の理由によります。

- row 単位の抽出がしやすい
- `groupby` による集約がしやすい
- FITS でも `SINGLE DISH` table に自然に対応する

### 4.3 処理履歴 provenance

処理履歴は `history` の dict に置きます。

ここで重要なのは、値が flat な文字列だけとは限らず、

- `dict`
- `list[dict]`
- 数値
- 真偽値

になり得る点です。

例:

- `baseline_history` は `list[dict]`
- `velocity_regrid` は `dict`
- `task` は文字列
- `t_hot_k` は数値

---

## 5. `stamp_scantable_code_provenance()` の意味

ソフトウェア provenance を押す中核関数は

```python
stamp_scantable_code_provenance(scantable, stage=...)
```

です。

この関数は、`Scantable` を **生成または更新した側** が呼ぶ想定です。

設定するもの:

- `meta["SWNAME"]`
- `meta["SWVER"]`
- `history["code_software"]`
- `history["code_version"]`
- `history["code_stage"]` 必要時のみ

この関数の意味は、

> この時点の `Scantable` は、どのソフトのどの版が生成・更新したか

を記録することです。

### 5.1 重要: `read_scantable()` では自動 stamp しない

これは今回の議論で最重要だった点です。

**読み込んだだけで producer provenance を書き換えてはいけません。**

例えば、古い版のコードで生成した FITS を新しい版のコードで開いたとき、単に `read_scantable()` しただけで

- `SWVER`
- `code_version`

が新しい版に変わってしまうと、

> どの版で解析された結果なのか

が分からなくなります。

したがって、現在の方針は次です。

- `read_scantable()` は stamp しない
- `filter_scantable()` や `run_velocity_regrid()` など、**新しい `Scantable` を返す処理**で stamp する
- `write_scantable()` は **未設定時だけ** `setdefault()` で `SWVER` / `code_version` を補完する

この「read では上書きしない」が provenance の信頼性にとって重要です。

---

## 6. FITS への対応づけ

現在の実装では、`Scantable` の provenance は SDFITS に次のように写像されます。

### 6.1 `meta` -> `PRIMARY` header 中心

`write_scantable()` -> `write_sdfits()` では、`meta` は基本的に `PRIMARY` header に書かれます。

ただし、全ての key がそのまま FITS keyword として入るわけではありません。

- FITS keyword として書けるものは `PRIMARY` header へ
- FITS keyword として扱いにくいものは、`HISTORY` 側へ `META:<KEY>` として退避されることがあります

### 6.2 `table` -> `SINGLE DISH` BinTable 列

`table` の各列は、基本的に `SINGLE DISH` extension の列になります。

そのため、row provenance は FITS 上でも table として保持されます。

### 6.3 `history` -> `HISTORY` extension

`history` は、専用の `HISTORY` BinTable extension に書かれます。

現在の形式は

- `KEY`
- `VALUE`

の 2 列です。

つまり FITS 上では

```text
HISTORY[KEY, VALUE]
```

という **key-value 表** になります。

---

## 7. `PRIMARY`, `SINGLE DISH`, `HISTORY` の役割分担

### 7.1 `PRIMARY`

主に `meta` が入ります。

典型例:

- `SWNAME`, `SWVER`
- `OBJECT`, `TELESCOP`, `INSTRUME`
- `SITELAT`, `SITELONG`, `SITEELEV`
- `OBSGEO-X`, `OBSGEO-Y`, `OBSGEO-Z`
- `CTYPE1`, `CUNIT1`, `CRVAL1`, `CDELT1`, `CRPIX1`
- `RESTFRQ`, `RESTFREQ`, `SPECSYS`, `SSYSOBS`

### 7.2 `SINGLE DISH`

主に `table` の列が入ります。

典型例:

- `FRONTEND`, `BACKEND`, `SAMPLER`
- `FDNUM`, `IFNUM`, `PLNUM`, `POLARIZA`
- `OBSMODE`, `CALSTAT`, `FLAGROW`
- `OBSFREQ`, `LO1FREQ`, `LO2FREQ`, `LO3FREQ`
- `TEMPSCAL`, `SPECSYS`, `SSYSOBS`
- `DATE-OBS`, `DATEOBS`, `TIMESTAMP`, `MJD`, `TIME`

### 7.3 `SINGLE DISH` header への昇格コピー

`write_sdfits()` では、`PRIMARY` にある一部 key を `SINGLE DISH` header にも複製します。

主な対象:

- `SITELAT`, `SITELONG`, `SITELON`, `SITEELEV`
- `OBSGEO-X`, `OBSGEO-Y`, `OBSGEO-Z`
- `CTYPE1`, `CUNIT1`, `CRVAL1`, `CDELT1`, `CRPIX1`
- `RESTFRQ`, `RESTFREQ`, `SPECSYS`, `SSYSOBS`, `VELDEF`, `VELOSYS`
- `TIMESYS`, `RADESYS`, `EQUINOX`
- `NCHAN`, `NCHANSEL`

これは downstream reader が `SINGLE DISH` 側だけを読む場合への互換性配慮です。

### 7.4 `HISTORY`

主に `history` の内容が入ります。

典型例:

- `code_software`, `code_version`, `code_stage`
- `task`, `input`, `created_at_utc`
- `baseline_history`
- `scale_history`
- `velocity_regrid`
- `qc_summary`
- `fitsio_migration`
- `META:<KEY>` 形式で退避された非 FITS 型 `meta`

---

## 8. 現在の write 時の補完動作

`write_scantable()` は、元の `Scantable` を壊さないように

- `table.copy()`
- `meta = dict(scantable.meta or {})`
- `history = dict(scantable.history or {})`

で作業した上で、未設定時のみ

```python
meta.setdefault("SWNAME", ...)
meta.setdefault("SWVER", ...)
history.setdefault("code_software", ...)
history.setdefault("code_version", ...)
```

を入れます。

意味は次の通りです。

- 解析途中で既に `stamp_scantable_code_provenance()` が入れていた provenance があれば、それを尊重する
- provenance が全く無い古い `Scantable` でも、少なくとも write 時点のソフト名・版は残す

したがって、**write-time fallback は補完であって上書きではありません。**

---

## 9. `read_scantable()` / `read_tastar_fits()` が返す provenance

### 9.1 `meta`

`read_tastar_fits()` は、SDFITS なら

- `PRIMARY` header
- `SINGLE DISH` header

を読み、必要な正規化をした上で `meta` にまとめます。

このため `meta` には、もともと `PRIMARY` にあった key に加えて、`SINGLE DISH` header 側の WCS / site 情報が入ることがあります。

### 9.2 `table`

`SINGLE DISH` table から `DATA` / `SPECTRUM` を除いた列を `pandas.DataFrame` として返します。

row provenance はここを見ます。

### 9.3 `history`

`HISTORY` extension があれば、`KEY`, `VALUE` を読んで `dict` にします。

ただしここで重要なのは、現在の実装では `VALUE` は **まず文字列として読まれる** という点です。

つまり、書き込み時に `dict` や `list` だったものも、読み込み後は多くの場合

- `"{'stage': 'baseline_fit', ...}"`
- `"[{'stage': 'write_sdfits', ...}]"`

のような **文字列** として `history` に入ります。

この制約は後で詳しく説明します。

---

## 10. 実際の確認方法

### 10.1 まずは `Scantable` として読む

```python
import sd_radio_spectral_fits as sdrsf

sc = sdrsf.read_scantable("your_file.fits")
```

### 10.2 ソフトウェア provenance を見る

```python
print(sc.meta.get("SWNAME"))
print(sc.meta.get("SWVER"))
print(sc.history.get("code_software"))
print(sc.history.get("code_version"))
print(sc.history.get("code_stage"))
```

### 10.3 row provenance を見る

```python
cols = [c for c in [
    "FRONTEND", "BACKEND", "SAMPLER",
    "FDNUM", "IFNUM", "PLNUM", "POLARIZA",
    "OBSFREQ", "LO1FREQ", "LO2FREQ", "LO3FREQ",
    "SIDEBAND", "SB1", "SB2", "SB3",
    "TEMPSCAL", "CALSTAT", "SPECSYS", "SSYSOBS",
] if c in sc.table.columns]

print(cols)
print(sc.table[cols].head())
```

### 10.4 `history` の key 一覧を見る

```python
print(sorted(sc.history.keys()))
```

### 10.5 `history` の個別項目を見る

```python
print(sc.history.get("task"))
print(sc.history.get("baseline_history"))
print(sc.history.get("scale_history"))
print(sc.history.get("velocity_regrid"))
print(sc.history.get("qc_summary"))
print(sc.history.get("fitsio_migration"))
```

---

## 11. `summarize_provenance()` の使い方

`summary = summarize_provenance(sc_or_path)` は provenance を人間が追いやすい形に整理した dict を返します。

例:

```python
import sd_radio_spectral_fits as sdrsf

prov = sdrsf.summarize_provenance("your_file.fits")
print(prov.keys())
```

主な返り値:

- `source`
- `n_rows`
- `meta_summary`
- `meta_keys_all`
- `row_provenance_columns`
- `row_provenance_unique`
- `row_provenance_groups`
- `history_keys_all`
- `history_focus`
- `history_parsed`

### 11.1 `meta_summary`

`meta` のうち、よく見る key を優先して抜き出したものです。

```python
print(prov["meta_summary"])
```

### 11.2 `row_provenance_columns`

provenance として重要な列のうち、実際に存在する列の一覧です。

```python
print(prov["row_provenance_columns"])
```

### 11.3 `row_provenance_unique`

各 provenance 列について、非 null 件数・unique 数・代表値を要約します。

```python
print(prov["row_provenance_unique"]["FDNUM"])
print(prov["row_provenance_unique"]["FRONTEND"])
```

### 11.4 `row_provenance_groups`

複数列を組み合わせた provenance key ごとに、何行あるかを集約します。

```python
print(prov["row_provenance_groups"]["group_columns"])
for row in prov["row_provenance_groups"]["groups"][:5]:
    print(row)
```

これは、例えば

- `(FRONTEND, BACKEND, SAMPLER, FDNUM, IFNUM, PLNUM)`

の組ごとに `n_rows` を見る用途に向いています。

### 11.5 `history_focus`

`history` のうち、よく見る項目だけを抜き出したものです。

```python
print(prov["history_focus"])
```

現在優先している主な key:

- `task`, `input`, `created_at_utc`
- `baseline_history`
- `scale_history`
- `velocity_regrid`
- `qc_summary`
- `notes`
- `stage`, `mode`, `group_mode`
- `rows`, `exclude_rows`
- `t_hot_k`, `tau_zenith`
- `rest_freq`, `v_corr_col`
- `fitsio_migration`

### 11.6 `history_parsed`

`history` の値に対して、可能なら `dict` / `list` / 数値 / 真偽値へ再解釈したものです。

```python
print(type(prov["history_parsed"].get("baseline_history")))
print(prov["history_parsed"].get("baseline_history"))
```

---

## 12. `print_provenance()` の使い方

```python
import sd_radio_spectral_fits as sdrsf

sdrsf.print_provenance("your_file.fits")
```

これで

- `meta_summary`
- `row_provenance_columns`
- `row_provenance_groups`
- `history_keys_all`
- `history_focus` または `history_parsed`

をまとめて表示できます。

### 12.1 `show_all_history`

```python
sdrsf.print_provenance("your_file.fits", show_all_history=True)
```

とすると、`history_focus` ではなく `history_parsed` 全体を表示します。

### 12.2 `... (truncated) ...` の意味

`print_provenance()` では、表示文字数が `max_history_chars` を超えると

```text
... (truncated) ...
```

が付きます。

これは **画面表示だけを省略した** という意味であり、provenance 自体が失われたわけではありません。

全部見たい場合は、例えば

```python
sdrsf.print_provenance(
    "your_file.fits",
    show_all_history=True,
    max_history_chars=100000,
)
```

のように `max_history_chars` を大きくします。

あるいは `summarize_provenance()` の返り値をそのまま見れば、省略なしで扱えます。

---

## 13. `astropy.io.fits` で直接確認する方法

### 13.1 HDU 構成を見る

```python
from astropy.io import fits

with fits.open("your_file.fits") as hdul:
    hdul.info()
```

典型的には

- `PRIMARY`
- `SINGLE DISH`
- `HISTORY`

が見えます。

### 13.2 `PRIMARY` header の software provenance を見る

```python
from astropy.io import fits

with fits.open("your_file.fits") as hdul:
    hdr = hdul[0].header
    print(hdr.get("SWNAME"))
    print(hdr.get("SWVER"))
```

### 13.3 `SINGLE DISH` table の row provenance を見る

```python
from astropy.io import fits

with fits.open("your_file.fits") as hdul:
    tab = hdul["SINGLE DISH"].data
    for name in ["FRONTEND", "BACKEND", "SAMPLER", "FDNUM", "IFNUM", "PLNUM"]:
        if name in tab.names:
            print(name, tab[name][:5])
```

### 13.4 `HISTORY` extension を見る

```python
from astropy.io import fits

with fits.open("your_file.fits") as hdul:
    htab = hdul["HISTORY"].data
    history = {
        str(k).strip(): str(v).strip()
        for k, v in zip(htab["KEY"], htab["VALUE"])
    }

print(history.get("code_version"))
print(history.get("baseline_history"))
print(history.get("scale_history"))
```

---

## 14. よく使う確認例

### 14.1 どのソフト版で解析されたかを見る

```python
sc = sdrsf.read_scantable("your_file.fits")

print("SWNAME =", sc.meta.get("SWNAME"))
print("SWVER  =", sc.meta.get("SWVER"))
print("code_software =", sc.history.get("code_software"))
print("code_version  =", sc.history.get("code_version"))
print("code_stage    =", sc.history.get("code_stage"))
```

### 14.2 どの beam / IF / 偏波が含まれているかを見る

```python
keycols = [c for c in ["FRONTEND", "BACKEND", "SAMPLER", "FDNUM", "IFNUM", "PLNUM"] if c in sc.table.columns]
print(sc.table[keycols].drop_duplicates())
```

### 14.3 `baseline_history` を構造として見たい

`read_scantable()` 後は文字列になっていることがあるので、

```python
prov = sdrsf.summarize_provenance(sc)
print(type(prov["history_parsed"].get("baseline_history")))
print(prov["history_parsed"].get("baseline_history"))
```

とするのが楽です。

### 14.4 どの provenance 組が何行あるか見る

```python
prov = sdrsf.summarize_provenance(sc)
for row in prov["row_provenance_groups"]["groups"]:
    print(row)
```

---

## 15. 重要な制約: `HISTORY` は FITS round-trip で型が崩れる

これは現在の実装上、最も重要な制約です。

### 15.1 何が起こるか

`history` には本来

- `dict`
- `list[dict]`
- 数値
- 真偽値

を入れられますが、FITS に書くとき `HISTORY` extension は `KEY`, `VALUE` の 2 列なので、最終的に多くのものは **文字列化** されます。

そのため、例えば

```python
history["baseline_history"] = [
    {"stage": "baseline_fit", "poly_order": 3}
]
```

のようなものも、FITS を保存して再読込した後は

```python
"[{'stage': 'baseline_fit', 'poly_order': 3}]"
```

のような文字列になることがあります。

### 15.2 現在どう対処しているか

`scantable_utils.summarize_provenance()` は、文字列が

- JSON
- Python literal
- bool / None / number

らしければ再解釈を試みます。

そのため、**読むときは `history` 生の dict より `history_parsed` を見る方が便利**です。

### 15.3 ただし完全ではない

文字列の再解釈は best-effort です。

- 文字列の表現が曖昧
- 非 JSON 的な表現
- 複雑な型

では、元の構造を完全に復元できないことがあります。

---

## 16. Scantable と FITS の非対称性

今回の議論で重要だったのは、`Scantable` と FITS の間にはいくつか **完全には対称でない部分** があることです。

### 16.1 `history` の型対称性が崩れる

これは前節の通りです。

- `Scantable.history` では構造化データを持てる
- FITS `HISTORY` では `KEY`, `VALUE` の文字列表になる

### 16.2 `meta` の一部が `SINGLE DISH` header に複製される

write 時に downstream compatibility のため、`PRIMARY` の一部 key を `SINGLE DISH` header にも複製しています。

そのため、論理的には 1 つの値でも、FITS 上では複数箇所に見えることがあります。

### 16.3 `meta` の一部は `META:<KEY>` として `HISTORY` へ逃がされる

FITS keyword にしにくい `meta` は `HISTORY` へ退避されることがあります。

つまり `meta` は write 後に

- `PRIMARY` keyword
- `SINGLE DISH` header
- `HISTORY` の `META:<KEY>`

へ分散し得ます。

### 16.4 read では producer provenance を更新しない

これは意図的な非対称性です。

- write / processing 側は provenance を更新する
- read 側は provenance を保持する

こうしないと「誰が作ったか」と「今誰が読んでいるか」が混ざってしまいます。

---

## 17. 将来整理した方がよい点

### 17.1 `HISTORY` を typed / JSON 化する

現在の `KEY`, `VALUE` 2 列方式は簡単ですが、構造化 provenance の round-trip に弱いです。

将来的には例えば

- `KEY`, `JSON_VALUE`, `TYPE`

のような table にするか、あるいは

- 1 つの `HISTORY_JSON` / `PROVENANCE_JSON` extension に丸ごと JSON で保存する

方が型保存の観点では堅牢です。

### 17.2 software provenance の表現を一本化する

現在は

- `meta["SWNAME"]`, `meta["SWVER"]`
- `history["code_software"]`, `history["code_version"]`

の両方を持っています。

これは実用上便利ですが、役割分担をさらに明確にする余地があります。

例:

- `SWNAME`, `SWVER` は header 用の短い識別子
- `history` 側はより詳細な provenance 専用

### 17.3 `summarize_provenance()` の先頭表示を software provenance 中心にする

現在の要約は `meta_summary` と `history_focus` が中心です。

今後、`SWVER` / `code_version` / `code_stage` をより前面に出す専用要約ブロックを追加すると、

- どの版で解析したか
- どの stage で生成したか

が一目で分かりやすくなります。

### 17.4 in-place 更新系と copy-return 系の stamp 方針を統一する

現在は `Scantable` を新しく返す関数を中心に stamp しています。

将来的には、

- 新しい `Scantable` を返す関数
- in-place で内容を更新する関数

の provenance policy を文書化して統一すると、`code_stage` の意味がさらに明確になります。

---

## 18. 実務上の推奨ワークフロー

### 18.1 まずファイルを読む

```python
sc = sdrsf.read_scantable("your_file.fits")
```

### 18.2 まず provenance 全体をざっと見る

```python
sdrsf.print_provenance(sc)
```

### 18.3 詳しく見る

```python
prov = sdrsf.summarize_provenance(sc)

print(prov["meta_summary"])
print(prov["row_provenance_groups"])
print(prov["history_focus"])
```

### 18.4 software version を明示確認する

```python
print(sc.meta.get("SWVER"))
print(sc.history.get("code_version"))
print(sc.history.get("code_stage"))
```

### 18.5 `history` の中の構造化項目は `history_parsed` を見る

```python
print(prov["history_parsed"].get("baseline_history"))
print(prov["history_parsed"].get("scale_history"))
print(prov["history_parsed"].get("velocity_regrid"))
```

---

## 19. まとめ

現在の provenance 実装は、次の方針で整理されています。

- **ファイル全体の属性**は `meta`
- **row ごとの由来**は `table`
- **処理履歴**は `history`
- **software provenance** は `meta` と `history` の両方に持つ
- **read では producer provenance を上書きしない**
- **write は未設定時だけ補完**する
- FITS では `history` が `KEY`, `VALUE` table になるため、**構造化データは文字列化されやすい**

したがって、現状で provenance を確認するときの実用的な順序は

1. `read_scantable()` で読む
2. `print_provenance()` で全体を見る
3. `summarize_provenance()` の `history_parsed` を使って中身を追う
4. 必要なら `astropy.io.fits` で `PRIMARY`, `SINGLE DISH`, `HISTORY` を直接確認する

です。

この方針は、今回ここで議論した重要点、特に

- どの版で解析したのかを追跡できること
- 読んだだけで provenance を壊さないこと
- FITS round-trip の制約を明示すること

を満たすためのものです。
