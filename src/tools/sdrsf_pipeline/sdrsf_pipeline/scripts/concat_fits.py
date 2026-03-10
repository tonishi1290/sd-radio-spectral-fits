from __future__ import annotations

import argparse
import numpy as np
import pandas as pd

from ..fitsio import read_tastar_fits, write_tastar_sdfits
from ..axis import freq_axis_from_wcs, slice_channels

def _wcs_key(meta: dict):
    return (
        meta.get("ctype1",""),
        meta.get("cunit1",""),
        float(meta["crval1_hz"]),
        float(meta["crpix1"]),
        float(meta["cdelt1_hz"]),
        float(meta["rest_hz"]),
        meta.get("specsys",""),
        meta.get("timesys",""),
        int(meta.get("nchan",0)),
    )


def _regrid_to_ref(meta_ref: dict, y: np.ndarray, meta_src: dict) -> np.ndarray:
    """Regrid 1D spectrum y from meta_src frequency axis to meta_ref frequency axis (linear interp).

    This is used only when --regrid-mismatch-to-first is enabled.
    The interpolation is done in frequency (Hz) on the axis defined by FITS-like WCS.
    """
    f_ref = freq_axis_from_wcs(meta_ref, nchan=int(meta_ref.get("nchan", y.size)))
    f_src = freq_axis_from_wcs(meta_src, nchan=int(meta_src.get("nchan", y.size)))
    # np.interp requires ascending x; handle descending axes
    if f_src[0] > f_src[-1]:
        f_src2 = f_src[::-1]
        y2 = y[::-1]
    else:
        f_src2 = f_src
        y2 = y
    if f_ref[0] > f_ref[-1]:
        f_ref2 = f_ref[::-1]
        rev = True
    else:
        f_ref2 = f_ref
        rev = False
    yy = np.interp(f_ref2, f_src2, y2, left=np.nan, right=np.nan)
    if rev:
        yy = yy[::-1]
    return yy

def main():
    ap = argparse.ArgumentParser(description="Concatenate Ta* FITS along spectrum axis (rows). Requires identical spectral WCS by default.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-dumps", type=int, default=0, help="Use only first N rows from each input (0=all)")
    ap.add_argument("--ch-start", type=int, default=None, help="Pre-slice channels (0-based inclusive)")
    ap.add_argument("--ch-stop", type=int, default=None, help="Pre-slice channels (0-based exclusive)")
    ap.add_argument("--regrid-mismatch-to-first", action="store_true", help="If WCS differs, regrid each spectrum to the first file's frequency axis before concatenation.")
    ap.add_argument("--allow-rest-mismatch", action="store_true", help="Allow RESTFRQ mismatch when regridding (NOT recommended).")
    ap.add_argument("inputs", nargs="+")
    args = ap.parse_args()

    metas = []
    datas = []
    tables = []
    histories = []

    for p in args.inputs:
        meta, data, table, hist = read_tastar_fits(p)
        data = np.asarray(data, float)
        if args.max_dumps and args.max_dumps > 0:
            data = data[: int(args.max_dumps), :]
            table = table.iloc[: int(args.max_dumps)].copy()
        meta, data = slice_channels(meta, data, args.ch_start, args.ch_stop)

        metas.append(meta)
        datas.append(data)
        tables.append(table)
        histories.append(hist)

    k0 = _wcs_key(metas[0])
    for i, m in enumerate(metas[1:], start=1):
        if _wcs_key(m) != k0:
            if not args.regrid_mismatch_to_first:
                raise SystemExit(f"WCS mismatch for input #{i}: {args.inputs[i]} (use sdrsf-regrid-and-concat or enable --regrid-mismatch-to-first)")
            # else: allow; will regrid below

    # If enabled, regrid spectra to the first file's frequency axis when WCS differs
    if args.regrid_mismatch_to_first:
        ref = metas[0]
        if (not args.allow_rest_mismatch):
            for i, m in enumerate(metas[1:], start=1):
                if float(m["rest_hz"]) != float(ref["rest_hz"]):
                    raise SystemExit(f"RESTFRQ mismatch for input #{i}: {args.inputs[i]} (use --allow-rest-mismatch to override)")
        datas2 = []
        k0 = _wcs_key(ref)
        for i, (m, d) in enumerate(zip(metas, datas)):
            if _wcs_key(m) == k0:
                datas2.append(np.asarray(d, float))
            else:
                # regrid each row
                A = np.asarray(d, float)
                B = np.vstack([_regrid_to_ref(ref, A[j], m) for j in range(A.shape[0])])
                datas2.append(B)
        data2 = np.vstack(datas2)
    else:
        data2 = np.vstack(datas)
    table2 = pd.concat(tables, axis=0, ignore_index=True)

    history = dict(stage="concat_fits", inputs=args.inputs, from_histories=histories,
                   max_dumps=int(args.max_dumps), ch_slice=dict(ch_start=args.ch_start, ch_stop=args.ch_stop),
                   regrid_mismatch_to_first=bool(args.regrid_mismatch_to_first), allow_rest_mismatch=bool(args.allow_rest_mismatch))
    write_tastar_sdfits(args.out, meta=metas[0], data=data2, table=table2, history=history)
    print(f"saved: {args.out}")

if __name__ == "__main__":
    main()
