#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Export pipeline internal FITS to CASA-friendly SDFITS-like FITS.

The pipeline uses an internal working format:
  PRIMARY + ImageHDU("DATA") + BinTableHDU("DUMPS") (+ optional "HISTORY")

CASA/ASAP typically expects SDFITS:
  PRIMARY + BinTableHDU(EXTNAME="SINGLE DISH") containing a vector spectrum column

This script converts the internal format to an SDFITS-like output. It does *not*
change the spectral WCS; it simply repackages the spectra and per-dump metadata.
"""

from __future__ import annotations

import argparse

from ..fitsio import read_tastar_fits, write_tastar_sdfits


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sdrsf-export-sdfits",
        description="Export internal Ta* FITS (DATA+DUMPS) to SDFITS-like FITS (SINGLE DISH BinTable).",
    )
    p.add_argument("input_fits", help="Input FITS (pipeline internal format; or already SDFITS).")
    p.add_argument("-o", "--out", required=True, help="Output FITS path (SDFITS-like).")
    p.add_argument(
        "--spectrum-column",
        default="DATA",
        help="Name of the spectrum column in the output SINGLE DISH table (default: DATA).",
    )
    p.add_argument(
        "--no-flag",
        action="store_true",
        help="Do not write FLAG column (default: write all-False FLAG).",
    )
    p.add_argument("--overwrite", action="store_true", help="Overwrite output file if exists.")
    return p


def main(argv: list[str] | None = None) -> int:
    p = build_parser()
    args = p.parse_args(argv)

    meta, data, table, history = read_tastar_fits(args.input_fits)
    write_tastar_sdfits(
        args.out,
        meta,
        data,
        table,
        history=history,
        spectrum_column=args.spectrum_column,
        include_flag=not args.no_flag,
        overwrite=args.overwrite,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
