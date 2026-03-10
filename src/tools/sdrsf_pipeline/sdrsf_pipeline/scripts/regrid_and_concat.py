from __future__ import annotations

import argparse
import datetime
import numpy as np
import pandas as pd

from ..fitsio import read_tastar_fits, write_tastar_sdfits
from ..axis import freq_axis_from_wcs, slice_channels


def _interp_to_freq_grid(freq_src: np.ndarray, y_src: np.ndarray, freq_tgt: np.ndarray) -> np.ndarray:
    """Interpolate one spectrum onto the target frequency grid (Hz).

    - Handles ascending/descending axes robustly.
    - Uses linear interpolation (np.interp).
    - Out-of-range is filled with NaN.
    """
    x = np.asarray(freq_src, float)
    y = np.asarray(y_src, float)

    # np.interp requires increasing x
    if x.size <= 1:
        return np.full(freq_tgt.size, np.nan, dtype=float)
    if x[0] > x[-1]:
        x2 = x[::-1]
        y2 = y[::-1]
    else:
        x2 = x
        y2 = y

    xt = np.asarray(freq_tgt, float)
    rev = False
    if xt[0] > xt[-1]:
        xt2 = xt[::-1]
        rev = True
    else:
        xt2 = xt

    yy = np.interp(xt2, x2, y2, left=np.nan, right=np.nan)
    if rev:
        yy = yy[::-1]
    return yy


def _wcs_key(meta: dict):
    """Tuple key for WCS equality check.

    Notes:
      - float-cast numeric keys to avoid string-formatting mismatches.
      - Includes RESTFRQ/specsys/timesys; these MUST match for a strict 'same WCS'.
    """
    return (
        str(meta.get("ctype1", "")),
        str(meta.get("cunit1", "")),
        float(meta.get("crval1_hz")),
        float(meta.get("crpix1")),
        float(meta.get("cdelt1_hz")),
        float(meta.get("rest_hz")),
        str(meta.get("specsys", "")),
        str(meta.get("timesys", "")),
        int(meta.get("nchan", 0)),
    )


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Regrid-and-concatenate Ta* FITS.\n"
            "\n"
            "Policy:\n"
            "  - The first input defines the *reference* frequency axis.\n"
            "  - If subsequent inputs have identical WCS (including RESTFRQ/specsys/timesys), they are appended as-is.\n"
            "  - Otherwise, they are regridded (linear interpolation) to the reference frequency axis and then appended.\n"
            "\n"
            "This tool is primarily for *frequency-axis* Ta* dump FITS (CTYPE1='FREQ').\n"
            "For Doppler-aligned coadding on a velocity grid, use sdrsf-coadd-fits-vlsrk.\n"
        )
    )
    ap.add_argument("--out", required=True, help="Output FITS (same WCS as the first input)")
    ap.add_argument("--allow-rest-mismatch", action="store_true", help="Allow REST frequency mismatch (NOT recommended)")
    ap.add_argument("--max-dumps", type=int, default=0, help="Use only first N rows from each input (0=all)")
    ap.add_argument("--ch-start", type=int, default=None, help="Pre-slice channels (0-based inclusive)")
    ap.add_argument("--ch-stop", type=int, default=None, help="Pre-slice channels (0-based exclusive)")
    ap.add_argument("inputs", nargs="+", help="Input FITS files (Ta* dumps)")

    args = ap.parse_args()
    if not args.inputs:
        raise SystemExit("No input FITS given.")

    # ------------------------------------------------------------
    # Read reference (first) file
    # ------------------------------------------------------------
    meta0, data0, tab0, hist0 = read_tastar_fits(args.inputs[0])
    data0 = np.asarray(data0, float)
    if args.max_dumps and args.max_dumps > 0:
        data0 = data0[: int(args.max_dumps), :]
        tab0 = tab0.iloc[: int(args.max_dumps)].copy()
    meta0, data0 = slice_channels(meta0, data0, args.ch_start, args.ch_stop)

    freq_tgt = freq_axis_from_wcs(meta0, nchan=data0.shape[1])

    out_datas = [np.asarray(data0, float)]
    out_tabs = [tab0.reset_index(drop=True)]
    histories = [hist0]

    regrid_log = [{"file": args.inputs[0], "action": "keep", "reason": "reference"}]

    # ------------------------------------------------------------
    # Process remaining inputs
    # ------------------------------------------------------------
    for fn in args.inputs[1:]:
        meta, data, tab, hist = read_tastar_fits(fn)
        data = np.asarray(data, float)

        if args.max_dumps and args.max_dumps > 0:
            data = data[: int(args.max_dumps), :]
            tab = tab.iloc[: int(args.max_dumps)].copy()

        meta, data = slice_channels(meta, data, args.ch_start, args.ch_stop)

        # REST frequency check (important for line identity)
        if (float(meta["rest_hz"]) != float(meta0["rest_hz"])) and (not args.allow_rest_mismatch):
            raise SystemExit(
                f"REST frequency mismatch: {fn} rest_hz={meta['rest_hz']} vs ref={meta0['rest_hz']}. "
                "If you really want to proceed, use --allow-rest-mismatch (NOT recommended)."
            )

        # Check whether WCS is identical (after any pre-slicing)
        same = (_wcs_key(meta) == _wcs_key(meta0))

        A = np.asarray(data, float)
        if same and (A.shape[1] == out_datas[0].shape[1]):
            out_datas.append(A)
            out_tabs.append(tab.reset_index(drop=True))
            regrid_log.append({"file": fn, "action": "keep", "reason": "WCS identical"})
        else:
            # Regrid each row in frequency domain to target freq grid
            freq_src = freq_axis_from_wcs(meta, nchan=A.shape[1])
            B = np.empty((A.shape[0], freq_tgt.size), dtype=float)
            for i in range(A.shape[0]):
                B[i] = _interp_to_freq_grid(freq_src, A[i], freq_tgt)

            out_datas.append(B)
            out_tabs.append(tab.reset_index(drop=True))
            regrid_log.append({"file": fn, "action": "regrid", "reason": "WCS mismatch"})

        histories.append(hist)

    data_out = np.vstack(out_datas)
    tab_out = pd.concat(out_tabs, axis=0, ignore_index=True)

    history = dict(
        stage="regrid_and_concat",
        created_at_utc=datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        inputs=list(args.inputs),
        allow_rest_mismatch=bool(args.allow_rest_mismatch),
        max_dumps=int(args.max_dumps),
        ch_slice=dict(ch_start=args.ch_start, ch_stop=args.ch_stop),
        regrid_log=regrid_log,
        from_histories=histories,
    )
    write_tastar_sdfits(
        args.out,
        meta=meta0,
        data=data_out,
        table=tab_out,
        history=history,
        overwrite=False,
    )
    print(f"saved: {args.out}")


if __name__ == "__main__":
    main()
