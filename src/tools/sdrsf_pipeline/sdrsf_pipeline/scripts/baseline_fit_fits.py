# sdrsf_pipeline/scripts/baseline_fit_fits.py
from __future__ import annotations

import argparse
# 新しいライブラリ関数をインポート
from sd_radio_spectral_fits.baseline import run_baseline_fit


def main():
    ap = argparse.ArgumentParser(description="Baseline-fit spectra in a Ta* FITS (multi-window polynomial).")
    ap.add_argument("--in", dest="inp", required=True, help="Input FITS (DATA: Nspec x Nchan)")
    ap.add_argument("--out", required=True, help="Output FITS (baseline-subtracted)")
    ap.add_argument("--max-dumps", type=int, default=0, help="Use only first N spectra (0=all)")
    ap.add_argument("--ch-start", type=int, default=None, help="Pre-slice channels (0-based inclusive)")
    ap.add_argument("--ch-stop", type=int, default=None, help="Pre-slice channels (0-based exclusive)")
    ap.add_argument("--v-corr-col", default="v_corr_kms", help="Velocity correction column in DUMPS table (km/s).")
    ap.add_argument("--vwin", action="append", required=True, help="Baseline velocity window vmin:vmax (km/s).")
    ap.add_argument("--poly", type=int, default=1, help="Baseline polynomial order (default 1).")
    ap.add_argument("--line-vwin", action="append", default=[], help="Line mask window vmin:vmax (km/s).")
    ap.add_argument("--iter-max", type=int, default=0, help="Sigma-clip iterations for baseline fit.")
    ap.add_argument("--iter-sigma", type=float, default=3.0, help="Sigma threshold for iterative rejection.")
    ap.add_argument("--on-fail", choices=["exit","skip","mask"], default="exit", help="If fit fails: exit, skip row, or mask with NaNs.")
    args = ap.parse_args()

    # run_baseline_fit を呼び出し (内部で Scantable の読み書き・処理が完結)
    run_baseline_fit(
        input_data=args.inp,
        output_path=args.out,
        vwin=args.vwin,
        poly_order=args.poly,
        line_vwin=args.line_vwin,
        iter_max=args.iter_max,
        iter_sigma=args.iter_sigma,
        max_dumps=args.max_dumps,
        ch_start=args.ch_start,
        ch_stop=args.ch_stop,
        v_corr_col=args.v_corr_col,
        on_fail=args.on_fail,
        overwrite=True
    )

if __name__ == "__main__":
    main()
