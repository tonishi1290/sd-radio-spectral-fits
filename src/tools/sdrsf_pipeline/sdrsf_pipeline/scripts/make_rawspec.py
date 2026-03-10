from __future__ import annotations

import argparse
import json
import pathlib

import pandas as pd

from ..rawspec import build_rawspec, save_rawspec
from ..axis import wcs_slice_channels


def _load_df(path: str) -> pd.DataFrame:
    p = pathlib.Path(path)
    if p.suffix.lower() in (".pkl", ".pickle"):
        return pd.read_pickle(p)
    if p.suffix.lower() in (".csv",):
        df = pd.read_csv(p)
        if "timestamp" in df.columns:
            df["timestamp"] = pd.to_datetime(df["timestamp"])
            df = df.set_index("timestamp")
        return df
    raise ValueError(f"Unsupported input format: {p}")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Build rawspec pickle from HOT/ON/OFF and meta/mapping.\n"
            "\n"
            "HOT/ON/OFF are expected to be (Ndump, Nchan) tables (pandas DataFrame recommended).\n"
            "The output is a pickle container used by subsequent pipeline stages.\n"
        )
    )
    ap.add_argument("--hot", required=True, help="HOT spectra DataFrame (pkl/csv). DatetimeIndex recommended.")
    ap.add_argument("--on", required=True, help="ON  spectra DataFrame (pkl/csv). DatetimeIndex recommended.")
    ap.add_argument("--off", required=True, help="OFF spectra DataFrame (pkl/csv). DatetimeIndex recommended.")
    ap.add_argument("--meta", required=True, help="Meta JSON file (must include spectral WCS keywords)")
    ap.add_argument("--mapping", default=None, help="Mapping table (csv or pkl), indexed by timestamp (recommended)")

    ap.add_argument("--out", required=True, help="Output rawspec pickle")

    # ------------------------------------------------------------
    # Ndump/Nchan reduction (before any Doppler correction)
    # ------------------------------------------------------------
    ap.add_argument(
        "--max-dumps",
        type=int,
        default=0,
        help=(
            "Use only the first N rows from HOT/ON/OFF (0=all). "
            "Mapping is aligned (reindexed) to the truncated ON timestamps when possible."
        ),
    )
    ap.add_argument("--max-hot", type=int, default=0, help="Override --max-dumps only for HOT (0=use --max-dumps).")
    ap.add_argument("--max-on", type=int, default=0, help="Override --max-dumps only for ON  (0=use --max-dumps).")
    ap.add_argument("--max-off", type=int, default=0, help="Override --max-dumps only for OFF (0=use --max-dumps).")
    ap.add_argument("--ch-start", type=int, default=None, help="Slice channels before saving (0-based inclusive).")
    ap.add_argument("--ch-stop", type=int, default=None, help="Slice channels before saving (0-based exclusive).")

    args = ap.parse_args()

    hot = _load_df(args.hot)
    on = _load_df(args.on)
    off = _load_df(args.off)

    meta = json.loads(pathlib.Path(args.meta).read_text(encoding="utf-8"))
    mapping = _load_df(args.mapping) if args.mapping else None

    # ------------------------------------------------------------
    # Optional: truncate dumps (rows) before building rawspec
    # ------------------------------------------------------------
    def _nmax(v: int) -> int:
        return int(v) if (v is not None and int(v) > 0) else 0

    n_default = _nmax(args.max_dumps)
    n_hot = _nmax(args.max_hot) if _nmax(args.max_hot) else n_default
    n_on = _nmax(args.max_on) if _nmax(args.max_on) else n_default
    n_off = _nmax(args.max_off) if _nmax(args.max_off) else n_default

    if n_hot:
        hot = hot.iloc[:n_hot, :].copy()
    if n_on:
        on = on.iloc[:n_on, :].copy()
    if n_off:
        off = off.iloc[:n_off, :].copy()

    # Align mapping to ON timestamps (recommended). If mapping is absent, build_rawspec will create empty mapping.
    if mapping is not None:
        try:
            if isinstance(on.index, pd.DatetimeIndex):
                mp = mapping.copy()
                if not isinstance(mp.index, pd.DatetimeIndex):
                    if "timestamp" in mp.columns:
                        mp["timestamp"] = pd.to_datetime(mp["timestamp"])
                        mp = mp.set_index("timestamp")
                mapping = mp.reindex(on.index)
        except Exception:
            # Keep original mapping if alignment fails; build_rawspec will still validate frame if possible.
            pass

    # ------------------------------------------------------------
    # Optional: channel slice BEFORE saving rawspec
    #   - Important: update meta WCS consistently (CRVAL1 shift, NCHAN)
    # ------------------------------------------------------------
    if (args.ch_start is not None) or (args.ch_stop is not None):
        nchan0 = int(on.shape[1])
        s = 0 if args.ch_start is None else int(args.ch_start)
        e = nchan0 if args.ch_stop is None else int(args.ch_stop)
        if not (0 <= s < e <= nchan0):
            raise SystemExit(f"Invalid channel slice {s}:{e} for Nchan={nchan0}")

        hot = hot.iloc[:, s:e].copy()
        on = on.iloc[:, s:e].copy()
        off = off.iloc[:, s:e].copy()

        meta = wcs_slice_channels(meta, s, e)

    raw = build_rawspec(hot=hot, on=on, off=off, meta=meta, mapping=mapping)
    save_rawspec(raw, args.out)
    print(f"saved: {args.out}")


if __name__ == "__main__":
    main()
