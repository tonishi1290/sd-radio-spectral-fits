# sdrsf_pipeline（独立ツール群）

このフォルダは、`sd_radio_spectral_fits` 本体を変更せずに使える「観測スペクトル解析パイプライン」ツール群です。

対応フロー（最小）:
1. HOT/ON/OFF（Ndump×Nchan）を **rawspec.pkl** に成形して保存
2. 1温度 chopper wheel により **Ta*（dump毎）** を作成（OFF/HOTは時間内挿）
3. baseline fitting（複数速度窓・多項式）
4. 同一座標点で加算（uniform または RMS重み）
5. Doppler補正が異なる複数FITSを **共通VLSRK速度軸に内挿**して加算
6. WCSが完全に一致するFITSの単純連結

重要:
- `sd_radio_spectral_fits` のソースには一切手を入れません（`tools/` 以下に閉じます）
- 設置場所（観測所情報）は **生データ側（raw meta / FITSヘッダー）に含まれる**ことを前提にします
  - 本ツールが書き出すFITSには `SITELAT/SITELON/SITEELEV/SITENAME` を保存します

インストール:
```bash
pip install -e "tools/sdrsf_pipeline[dev]"
```

ドキュメント:
- 詳細日本語ガイド: `docs/USER_GUIDE_ja.md`
- Detailed English guide: `docs/USER_GUIDE.md`


### 連結まわり
- `sdrsf-concat-fits` : WCS完全一致のみ（内挿なし）
- `sdrsf-regrid-and-concat` : 先頭FITSの周波数軸へ線形内挿してから連結（明示的にデータが変わる）


## 0.1.2 の更新点

- site 情報が欠損/None の場合でも不意に `float(None)` で落ちないように堅牢化しました (必要な場合は明示的にエラー)。
- `sdrsf-make-rawspec` で `--max-dumps` と `--ch-start/--ch-stop` による切り出しをサポートし、WCS を整合するよう更新します。
- `sdrsf-regrid-and-concat` でも `--max-dumps` とチャンネル切り出しを追加しました。
- `sdrsf-baseline-fit-fits` で dump ごとの速度補正列 `v_corr_kms` を用いた速度窓指定をサポートしました (`--v-corr-col`)。
