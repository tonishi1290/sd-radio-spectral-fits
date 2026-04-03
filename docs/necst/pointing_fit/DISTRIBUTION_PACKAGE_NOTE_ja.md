# 配布用パッケージについて

この配布用パッケージでは、容量削減のため `validation/` 配下の生成済み出力ディレクトリ（`*_out`）を削除しています。

残しているもの:
- 実行コード本体 `pointing_fit/`
- 実例生成用スクリプト `examples/`
- 説明書 `docs/`
- validation 再実行用のスクリプト類と補助資料

削除したもの:
- validation の大量生成物（図、CSV、json、report などの出力結果）
- `__pycache__`, `*.pyc`

必要な場合は、`validation/` 配下のスクリプトを実行して検証出力を再生成してください。
