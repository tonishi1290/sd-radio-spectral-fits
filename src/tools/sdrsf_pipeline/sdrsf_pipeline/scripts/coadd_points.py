from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from ..fitsio import read_tastar_fits, write_tastar_sdfits
from ..coadd import coadd_by_position
from ..baseline import BaselineModel
from ..ranges import parse_windows
from ..utils import subtract_windows
from ..axis import slice_channels

def main():
    ap = argparse.ArgumentParser(description="Coadd Ta* dumps by position (uniform or rms_weight).")
    ap.add_argument("--in", dest="inp", required=True, help="Input Ta* dumps FITS")
    ap.add_argument("--out", required=True, help="Output coadded FITS (Npos x Nchan)")
    ap.add_argument("--max-dumps", type=int, default=0, help="Use only first N dumps (0=all)")
    ap.add_argument("--ch-start", type=int, default=None, help="Pre-slice channels (0-based inclusive)")
    ap.add_argument("--ch-stop", type=int, default=None, help="Pre-slice channels (0-based exclusive)")
    ap.add_argument("--line-vwin", action="append", default=[], help="Line windows (VLSRK km/s) excluded from baseline/RMS; repeatable")
    ap.add_argument("--mode", choices=["uniform", "rms_weight"], default="uniform")
    ap.add_argument("--pos-col", default="pos_id", help="Column name for position ID in DUMPS table")
    ap.add_argument("--pos-tol-arcsec", type=float, default=None, help="If pos_id absent, group by RA/Dec within this tolerance")
    ap.add_argument("--rms-vwin", action="append", default=[], help="RMS velocity window 'vmin:vmax' (repeatable)")
    ap.add_argument("--rms-poly", type=int, default=0, help="Polynomial order for RMS estimation (default 0)")
    ap.add_argument("--baseline-vwin", action="append", default=[], help="Baseline velocity window 'vmin:vmax' (repeatable)")
    ap.add_argument("--baseline-poly", type=int, default=1, help="Baseline polynomial order (default 1)")
    ap.add_argument("--v-corr-col", default="v_corr_kms", help="Velocity correction column in DUMPS table")
    args = ap.parse_args()

    meta, data, table, hist_in = read_tastar_fits(args.inp)
    if args.max_dumps and args.max_dumps > 0:
        data = data[: int(args.max_dumps), :]
        table = table.iloc[: int(args.max_dumps)].copy()
    meta, data = slice_channels(meta, np.asarray(data, float), args.ch_start, args.ch_stop)
    if isinstance(table.index, pd.DatetimeIndex):
        if table.index.name is None:
            table.index.name = "timestamp"
        table = table.reset_index(drop=False)
    else:
        table = table.reset_index(drop=True)

    rms_vwins = parse_windows(args.rms_vwin) if args.rms_vwin else None
    base_vwins = parse_windows(args.baseline_vwin) if args.baseline_vwin else None
    line_vwins = parse_windows(args.line_vwin) if args.line_vwin else []
    if line_vwins and rms_vwins:
        rms_vwins = subtract_windows(rms_vwins, line_vwins)
    if line_vwins and base_vwins:
        base_vwins = subtract_windows(base_vwins, line_vwins)
    base_model = None
    if base_vwins:
        base_model = BaselineModel(poly_order=int(args.baseline_poly), v_windows_kms=base_vwins)

    res = coadd_by_position(
        meta=meta,
        data=data,
        table=table,
        mode=args.mode,
        pos_col=args.pos_col,
        tol_arcsec=args.pos_tol_arcsec,
        rms_vwins=rms_vwins,
        line_vwins=line_vwins,
        ch_slice=dict(ch_start=args.ch_start, ch_stop=args.ch_stop),
        max_dumps=int(args.max_dumps),
        rms_poly_order=int(args.rms_poly),
        baseline_model=base_model,
        v_corr_col=args.v_corr_col,
    )

    history = dict(
        stage="coadd_points",
        input=args.inp,
        mode=args.mode,
        pos_col=args.pos_col,
        pos_tol_arcsec=args.pos_tol_arcsec,
        rms_vwins=rms_vwins,
        line_vwins=line_vwins,
        ch_slice=dict(ch_start=args.ch_start, ch_stop=args.ch_stop),
        max_dumps=int(args.max_dumps),
        rms_poly_order=int(args.rms_poly),
        baseline_model=None if base_model is None else dict(poly_order=base_model.poly_order, v_windows_kms=base_model.v_windows_kms),
        from_history=hist_in,
        per_group=res.baseline_histories,
    )

    # write coadded as FITS: reuse DUMPS table schema as "POINTS"
    write_tastar_sdfits(
        args.out,
        meta=res.meta,
        data=res.data,
        table=res.table,
        history=history,
    )
    print(f"saved: {args.out}")

if __name__ == "__main__":
    main()
