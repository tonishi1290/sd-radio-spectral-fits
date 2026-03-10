# sdrsf_pipeline/scripts/coadd_fits_vlsrk.py
from __future__ import annotations

import argparse
# 新しいライブラリ関数をインポート
from sd_radio_spectral_fits.coadd import run_velocity_coadd


def main():
    ap = argparse.ArgumentParser(description="Coadd Ta* FITS (Correct handling of FREQ+LSRK).")
    ap.add_argument("--out", required=True, help="Output FITS")
    ap.add_argument("--mode", choices=["uniform", "rms_weight"], default="uniform")
    ap.add_argument("--pos-col", default="pos_id")
    ap.add_argument("--pos-tol-arcsec", type=float, default=None)
    ap.add_argument("--v-corr-col", default="v_corr_kms")
    ap.add_argument("--target-frame", default="icrs")
    ap.add_argument("--vmin", type=float, required=False)
    ap.add_argument("--vmax", type=float, required=False)
    ap.add_argument("--dv", type=float, required=False)
    ap.add_argument("--allow-outside-overlap", action="store_true")
    ap.add_argument("--axis", choices=["freq", "vel"], default="freq", help="Output spectral axis type.")
    ap.add_argument("--rest-freq", type=float, default=None, help="Force rest frequency (Hz).")
    ap.add_argument("--baseline-vwin", action="append", default=[])
    ap.add_argument("--baseline-poly", type=int, default=0)
    ap.add_argument("--baseline-iter-max", type=int, default=0)
    ap.add_argument("--baseline-iter-sigma", type=float, default=3.0)
    ap.add_argument("--rms-vwin", action="append", default=[])
    ap.add_argument("--rms-poly", type=int, default=0)
    ap.add_argument("--line-vwin", action="append", default=[])
    ap.add_argument("--block-size", type=int, default=0)
    ap.add_argument("--max-dumps", type=int, default=0)
    ap.add_argument("--ch-start", type=int, default=None)
    ap.add_argument("--ch-stop", type=int, default=None)
    ap.add_argument("--on-fail", choices=["exit", "skip", "mask"], default="exit")
    ap.add_argument("inputs", nargs="+")
    args = ap.parse_args()

    # run_velocity_coadd を呼び出し (内部で Scantable の読み書き・処理が完結)
    run_velocity_coadd(
        inputs=args.inputs,
        output_path=args.out,
        mode=args.mode,
        pos_col=args.pos_col,
        pos_tol_arcsec=args.pos_tol_arcsec,
        v_corr_col=args.v_corr_col,
        target_frame=args.target_frame,
        vmin=args.vmin,
        vmax=args.vmax,
        dv=args.dv,
        allow_outside_overlap=args.allow_outside_overlap,
        axis_type=args.axis,
        rest_freq=args.rest_freq,
        baseline_vwin=args.baseline_vwin,
        baseline_poly=args.baseline_poly,
        baseline_iter_max=args.baseline_iter_max,
        baseline_iter_sigma=args.baseline_iter_sigma,
        rms_vwin=args.rms_vwin,
        rms_poly=args.rms_poly,
        line_vwin=args.line_vwin,
        block_size=args.block_size,
        max_dumps=args.max_dumps,
        ch_start=args.ch_start,
        ch_stop=args.ch_stop,
        on_fail=args.on_fail,
        overwrite=True
    )

if __name__ == "__main__":
    main()
