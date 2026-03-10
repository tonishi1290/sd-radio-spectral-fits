#!/usr/bin/env python3
from __future__ import annotations

import argparse
from typing import Tuple

# ライブラリから高レベル関数をインポート
from sd_radio_spectral_fits.calibrate import run_tastar_calibration


def _parse_range_int(s: str) -> Tuple[int, int]:
    """Parse 'start:stop' into (start, stop)."""
    parts = s.split(":")
    if len(parts) != 2:
        raise ValueError("range must be 'start:stop'")
    a, b = sorted(parts)
    return int(a), int(b)


def _parse_range_float(s: str) -> Tuple[float, float]:
    """Parse 'vmin:vmax' into (vmin, vmax)."""
    parts = s.split(":")
    if len(parts) != 2:
        raise ValueError("range must be 'vmin:vmax'")
    a, b = sorted(parts)
    return float(a), float(b)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Calibrate raw spectra to Ta* (1-load chopper wheel)."
    )
    ap.add_argument("input_fits", help="Raw FITS input")
    ap.add_argument("-o", "--output", required=True, help="Output SDFITS path")
    ap.add_argument("--tamb-k", type=float, default=300.0, help="Ambient load temp [K]")
    ap.add_argument("--ch-range", help="Channel trim 'start:stop'")
    ap.add_argument("--vlsrk-range-kms", help="Velocity window 'vmin:vmax'")
    ap.add_argument("--target-frame", default="icrs", help="Target frame (default: icrs)")
    ap.add_argument("--store-freq-column", action="store_true", help="Store FREQ column (optional)")
    ap.add_argument("--spectrum-column", default="DATA", help="Output column name")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite output")

    args = ap.parse_args()

    # 引数のパース処理
    ch_range = _parse_range_int(args.ch_range) if args.ch_range else None
    vlsrk_range = _parse_range_float(args.vlsrk_range_kms) if args.vlsrk_range_kms else None

    # ライブラリ関数の呼び出し
    run_tastar_calibration(
        input_data=args.input_fits,
        output_path=args.output,
        tamb_k=args.tamb_k,
        ch_range=ch_range,
        vlsrk_range_kms=vlsrk_range,
        target_frame=args.target_frame,
        spectrum_column=args.spectrum_column,
        overwrite=args.overwrite,
        store_freq_column=args.store_freq_column
    )

if __name__ == "__main__":
    main()
