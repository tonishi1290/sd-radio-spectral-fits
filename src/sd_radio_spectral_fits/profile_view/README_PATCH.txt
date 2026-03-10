このパッチは plotting/windowmask.py を新設し、viewer/grid/montage を共通実装へ寄せました。
目的: rest_freq 上書きで Native→Display の速度シフトを行った場合でも、
  - Fit window(BSL_WINF) の描画
  - RMS(win) の計算マスク（raw/display 両モード）
が必ず同じ座標変換を通り、整合するようにする。

適用方法:
  1) リポジトリ root で、このパッチの 'src/' をあなたの 'src/' に上書きコピー（同名ファイルは置換）
  2) 追加ファイル windowmask.py は plotting パッケージ配下に入るだけです（他依存の追加なし）

含まれるファイル:
  src/sd_radio_spectral_fits/plotting/windowmask.py  (NEW)
  src/sd_radio_spectral_fits/plotting/viewer.py      (UPDATED)
  src/sd_radio_spectral_fits/plotting/grid.py        (UPDATED)
  src/sd_radio_spectral_fits/plotting/montage.py     (UPDATED)

windowmask.py の新API（他の plot 系へ横展開する際はこれを使用）:
  - fit_windows_xaxis(...)
      BSL_WINF (Native vel) → Display vel(rest_freq補正) → 現在のxaxis(vel/freq/chan) の窓へ変換
  - compute_rms_win(...)
      residual(y) の RMS(win) を、xr_vel(速度) と fit windows を両方考慮して一発で計算
      ※ raw/display モードでRMS用グリッドが変わっても、必ず同じ座標系でマスク生成

他の plot モジュールを同じ方針で寄せる手順（推奨）:
  1) リポジトリ全体で BSL_WINF と window 処理を検索
       rg -n "BSL_WINF|win_str|fit_windows_disp_vel|vel_windows_to_xaxis|recalculate_velocity_axis" src/
  2) Fit window 描画は
       wins_x = fit_windows_xaxis(win_str, rest_hz_meta, rest_hz_user, xaxis, vel_disp, freq_ghz, chan)
       for x0, x1 in wins_x: ax.axvspan(...)
     に統一
  3) RMS(win) 相当は
       rms, n, use = compute_rms_win(...)
     に置換（range+xaxis変換+windowシフト+raw/displayを丸ごと共通化）
