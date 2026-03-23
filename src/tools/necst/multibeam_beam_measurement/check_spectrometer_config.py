from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from .config_io import validate_spectrometer_config



def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Validate a converter-compatible spectrometer config for multi-beam analysis.")
    ap.add_argument("spectrometer_config", help="Path to spectrometer TOML")
    ap.add_argument("--stream-name", dest="stream_names", action="append", default=None, help="Explicit fit stream selection (repeatable)")
    ap.add_argument("--out-csv", default=None, help="Optional CSV path for the stream table")
    args = ap.parse_args(argv)
    result = validate_spectrometer_config(Path(args.spectrometer_config).expanduser().resolve(), explicit_stream_names=args.stream_names)
    print(result.stream_table.to_string(index=False))
    if result.primary_streams:
        print(f"\nfit_streams = {', '.join(result.primary_streams)}")
    if result.duplicate_beam_ids:
        print(f"duplicate_beam_ids (info) = {', '.join(result.duplicate_beam_ids)}")
    if result.warnings:
        print("\nwarnings:")
        for w in result.warnings:
            print(f"  - {w}")
    if args.out_csv:
        out = Path(args.out_csv).expanduser().resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        result.stream_table.to_csv(out, index=False)
        print(f"\n[done] wrote {out}")


if __name__ == "__main__":
    main()
