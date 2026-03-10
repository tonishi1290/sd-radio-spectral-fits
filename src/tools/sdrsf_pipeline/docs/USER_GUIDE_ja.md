# USER GUIDE（日本語版）: sdrsf_pipeline

この文書は「電波スペクトル解析に不慣れな人でも勘違いしにくい」ことを最優先に、あえて冗長に書いています。

---

## 0. このツールがやること／やらないこと

### やること
- HOT/ON/OFF（Ndump×Nchan）を **rawspec.pkl** に成形して保存
- 1温度 chopper wheel で **Ta*（dump毎）** を作る
- OFF/HOT は **時刻で内挿**（前後があれば線形、片側しか無ければ最近傍）
- チャンネル切り出し（ch範囲／VLSRK範囲）でデータ量削減  
  ※切り出し後も周波数軸WCSが壊れないように更新
- baseline fitting（複数速度窓、多項式）を再利用可能な形で実施  
  → 履歴（窓の個数・範囲・次数・係数）を **JSONとしてFITSに保存**（後で読める）
- 同じ座標点のスペクトルを、uniform平均またはRMS重み平均で加算
- Doppler補正が異なる複数FITSを、**共通VLSRK速度軸に内挿してから加算**（ヘッダーも整合）
- WCSが完全一致するFITSは、ドップラー補正なしで「そのまま」連結できる

### やらないこと（現段階）
- 望遠鏡ログを全部解析して自動でHOT/ON/OFFを分類する  
  → すでに「HOT/ON/OFFとして指定されている」前提を採用
- 観測所（設置場所）をCLIで毎回指定させる運用  
  → **生データ側（meta/FITS）に含まれている前提**で扱う

---

## 1. 入力データのモデル（最重要）

### 1.1 HOT/ON/OFF スペクトル配列
- それぞれ shape = **(Ndump, Nchan)**
- NumPy配列でも pandas.DataFrame でも良い
- 強く推奨: **DataFrame + UTCのDatetimeIndex**  
  → 「いつのdumpか」を確実に保持でき、内挿が安全になります

### 1.2 周波数軸（WCS）を meta に入れる理由
チャンネルを切り出したとき、周波数軸の計算を誤る事故が多いです。  
そこで、FITS互換のWCS（周波数軸）で統一します。

metaに必須:
- `nchan` : チャンネル数
- `crval1_hz` : 基準ピクセルの周波数（Hz）
- `crpix1` : 基準ピクセル番号（**1始まり**）
- `cdelt1_hz` : 1chあたり周波数刻み（Hz）。**LSB等なら負でOK**
- `rest_hz` : 静止周波数（Hz, RESTFRQ）
- `timesys="UTC"`

周波数は配列index i(0始まり)に対して
f(i) = crval1_hz + ((i+1) - crpix1) * cdelt1_hz

#### チャンネル切り出し時のWCS更新（事故防止の肝）
データを `A[:, ch_start:ch_stop]` に切り出したら、
- `nchan_new = ch_stop - ch_start`
- `crval1_hz_new = crval1_hz + ch_start * cdelt1_hz`
- `crpix1` はそのまま
とします（本ツールはこの規則で更新します）。

---

## 2. 観測所（設置場所）情報は「生データ側に含まれる前提」
VLSRK補正を“情報から計算”するには、観測所が必要です。

meta（raw段階）にどちらかを入れてください：
- `site_name`（Astropyのサイトレジストリ名）
- もしくは `site_lat_deg`, `site_lon_deg`, `site_height_m`

本ツールが書き出すTa* FITSには、ヘッダーへ以下も保存します：
- `SITELAT`, `SITELON`, `SITEELEV`, `SITENAME`

---

## 3. マッピング情報（同一点の定義）
ON dumpごとに、少なくとも以下のどちらかが必要です。

推奨（確実）:
- `pos_id`（整数ID）

代替（IDが無い場合）:
- `ra_deg`, `dec_deg`（度）
  + その上で `--pos-tol-arcsec`（許容差）でグルーピング

---

## 4. ステージA: rawspec.pkl を作る（rawは成形して保存するだけ）
- `examples/make_rawspec_example.py` を参照
- `sdrsf-make-rawspec` も利用可能

---

## 5. ステージB: chopper wheel（1温度）で Ta*（dump毎）を作る
### 5.1 Ta*の式（一般的な1温度）
- gain = Tamb / (HOT - OFF)
- Ta*  = (ON - OFF) * gain
Tambは `--tamb`（デフォルト300K）

### 5.2 OFF/HOTの内挿ルール（要件通り）
ONの時刻tに対して
- OFF(t): 前後のOFFがあれば線形内挿
- 片側しか無ければ最近傍
HOTも同様です

### 5.3 データ量削減：ch範囲／VLSRK範囲で切り出し
- ch範囲: `--ch-start/--ch-stop`
- VLSRK範囲: `--vmin/--vmax`
  - VLSRK補正は、観測時刻(UTC)＋観測所＋(RA,Dec) から自動計算します
  - そのVLSRK軸で指定した速度範囲を含むように **ch範囲を自動決定**します
  - 決定したch範囲に基づいてWCSを安全に更新します

---

## 6. ステージC: baseline fitting（複数窓・再利用）
`--vwin "vmin:vmax"` を **複数回**指定できます。
例:
- `--vwin "-100:-50" --vwin "50:100"`

次数:
- `--poly 0`（定数）
- `--poly 1`（直線）
- `--poly 2`（2次）…

履歴はFITSの `HISTORY` 拡張に JSON として保存され、
後からプログラムで読めます（ログ文字列ではありません）。

保存される主な情報:
- baseline窓の個数と範囲（速度）
- 多項式次数
- 係数
- 作成日時

---

## 7. ステージD: 同一点で加算（uniform / RMS重み）
- uniform: 単純平均
- rms_weight: 1/rms^2 で重み付け

RMSの定義（要件通り）:
1) 指定範囲（速度窓）を抽出
2) 多項式でfit（次数指定）
3) 残差residについて  
   RMS = sqrt( mean( (resid - mean(resid))^2 ) )

さらに、平均前に baseline を引くオプションもあります。

---

## 8. ステージE: Doppler補正が異なる複数FITSを加算
`--vmin/--vmax/--dv` で共通VLSRKグリッドを定義し、
各スペクトルをその速度軸に内挿してから加算します。  
チャネルがずれるので内挿が必須、という要件を満たします。

---

## 9. FITS連結（そのまま連結したい場合）
WCSが完全一致する場合だけ、`sdrsf-concat-fits` で安全に連結できます。  
WCSが異なるなら、必ず regrid/coadd 側を使ってください。

---

## 10. 生データスペクトルを作る際の指針（チェックリスト）
1. 各dumpに正しいUTC時刻を付ける（DataFrame index推奨）
2. HOT/ON/OFFの分類は厳密に（あなたの前提）
3. Nchanは全ブロックで一致させる
4. WCS（crval/crpix/cdelt/rest）を必ず埋める  
   ※切り出し事故防止のため
5. 観測所情報（site_name か lat/lon/height）を入れる  
   ※VLSRK補正を自動計算するため
6. ON dumpにRA/Decを入れる（pos_idがあるとさらに安全）
7. まずは examples と pytest で動作確認する



---

## 9b. 明示的な「内挿して連結」: sdrsf-regrid-and-concat（周波数軸）
`WCSが一致しないけれど、どうしても“連結して一つにまとめたい”` という場合に使います。

- **最初のFITS（先頭）**を基準（参照WCS）とし、
- 2本目以降は **周波数軸を先頭に合わせて線形内挿**してから、
- 行方向に連結します。

注意（重要）:
- これは **データを変更**します（内挿＝再サンプリング）。
- “そのまま連結”ではありません。
- RESTFRQが違う場合はデフォルトでエラー（線の同一性が崩れるため）。
  理解した上で続行する場合だけ `--allow-rest-mismatch` を使ってください。

例:
```bash
sdrsf-regrid-and-concat --out concat_on_ref.fits a.fits b.fits c.fits
```


---

## 付記: sdrsf-coadd-fits-vlsrk の既定のVLSRKグリッド挙動
- `--vmin/--vmax/--dv` を指定した場合はそれを使用します。
- 省略がある場合は自動決定します：
  - `dv` は **先頭FITSの周波数WCSとRESTFRQから導かれる速度刻み**を既定値にします。
  - `vmin/vmax` は **全入力がカバーする共通部分（overlap）**を使います
    （どの入力も寄与できる速度範囲だけを採用）。
- 決まったグリッドは、出力FITSのHISTORY(JSON)に保存されます。


---

## 11. コマンド別オプション（リファレンス）

### 11.1 共通オプション（多くのコマンドで共通）
- `--max-dumps N`  
  先頭から N 行（dump/point）だけ処理します。大容量データの素早い試行に便利です。  
  **重要**: これはドップラー補正や内挿の前に適用されます（要件）。

- `--ch-start i --ch-stop j`  
  チャンネル範囲 `[i, j)` を切り出します（0始まり）。  
  **重要**: 切り出し後も周波数軸WCSが破綻しないよう `crval1_hz` 等を更新します（要件）。

- `--on-fail {exit,skip,mask}`  
  計算不能が生じたときの動作です。  
  - `exit`: その時点でエラー終了（再現性・事故防止の観点で推奨）  
  - `skip`: そのスペクトル/点を捨てて続行  
  - `mask`: NaN等で埋めて続行（「逃げ道」）

### 11.2 baseline窓を複数個設定する方法（最重要）
baseline窓（`--vwin` / `--baseline-vwin` など）は **複数回指定**できます。

例（2つの速度窓を使う）:
```bash
--vwin "-100:-50" --vwin "50:100"
```

line窓（`--line-vwin`）も同様に複数回指定できます:
```bash
--line-vwin "-10:10" --line-vwin "120:130"
```

内部では
- baseline窓 = (指定した baseline窓) から (line窓) を差し引いたもの
として扱います。  
**line窓がbaseline窓と重ならない場合は、そのままです。**

### 11.3 反復フィット（補助オプション）
`--iter-max > 0`（または `--baseline-iter-max > 0`）で、
baselineフィット時に「外れ値（線の残り等）をσクリップして再フィット」を実行できます。

- 反復は **あくまで補助**です。  
  まずは `--line-vwin` などで線領域を明示的に除外するのが安全です。

- スペクトルが線だらけ等で、baseline窓が実質的に残らない場合は計算不能になります。  
  このときは既定では `exit` で止まります（事故防止）。必要なら `--on-fail mask` などを使ってください。

### 11.4 --block-size の目安（大容量I/O）
`--block-size` は「1つの位置グループ内で、何行ずつ処理するか」です。  
大きいデータでメモリを節約したいときに使います。

目安（概算）:
- 1行のスペクトルが `Nchan`、float64 なら 1行 ≈ `8*Nchan` byte  
- `block_size` 行を同時に処理するなら ≈ `8*Nchan*block_size` byte

例:
- Nchan=32768 のとき  
  `block_size=512` → 約 8*32768*512 ≈ 134 MB 程度

安全側に倒すなら:
- まず `block_size=256` から開始し、余裕があれば増やす、が無難です。

---

## 12. 連結（concat）と内挿（regrid）の方針
- **そのまま連結（concat）できるのは、WCSが完全一致する場合のみ**が安全です。  
- ただし `sdrsf-concat-fits --regrid-mismatch-to-first` を有効にすると、  
  WCSが一致しない入力も「先頭のWCSへ線形内挿してから」連結できます。  
  これは `sdrsf-regrid-and-concat` と同等の挙動で、**データが変更されます**。

推奨:
- 解析用に「そのまま並べたい」だけなら WCS一致のものだけ concat
- どうしても一つにまとめたい（行方向に結合したい）なら regrid-and-concat
- 「加算して1本にしたい」なら coadd（必要ならVLSRKグリッドへ）

---

## 13. sdrsf-coadd-fits-vlsrk の仕様（質問への明確な回答）
- 入力ごとにWCSが異なっても、行ごとに `v_lsrk` 軸を構築して共通グリッドへ内挿します。
- 速度グリッドは
  - 何も指定しなければ **入力群の共通重なり範囲（overlap）**を自動採用
  - `--vmin/--vmax/--dv` を指定したらそれを採用（ただし overlap 外は既定でエラー）
- `--vmin/--vmax` を overlap の外に出したい場合は `--allow-outside-overlap` が必要です
  （その場合、一部入力は NaN を含む可能性があります）



## Ndump/Nchan の切り出しに関する注意

多くのスクリプトで **重い処理の前にデータ量を削減** できます。

- `--max-dumps N` : 先頭 N 本のスペクトル(行)だけを使用
- `--ch-start S --ch-stop E` : チャンネルを切り出してから処理
  (WCS も整合するように CRVAL1 のシフトと NCHAN 更新を行います)

例:

```bash
sdrsf-make-rawspec ... --max-dumps 2000 --ch-start 1024 --ch-stop 3072
```

同様の考え方は `sdrsf-regrid-and-concat` や `sdrsf-baseline-fit-fits` でも利用できます。


## baseline fitting と VLSRK 速度窓

`sdrsf-baseline-fit-fits` は baseline 窓を **km/s** (radio 定義) で指定します。
FITS テーブルに dump ごとの速度補正列 (既定 `v_corr_kms`) がある場合は、

- `v_lsrk(ch,row) = v_radio(ch) + v_corr_kms[row]`

として各dumpで VLSRK 軸を構成し、同一の速度窓が時間変化しても破綻しにくいようにします。

補正列が無い場合は `v_corr_kms=0` とみなして処理を継続します (履歴に明示します)。
