from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from . import sunscan_extract_multibeam as extract_mod
from . import sunscan_fit_multibeam as fit_mod
from . import synthetic_multibeam as pseudo_mod
from .public_api import (
    check_spectrometer_config,
    run_multibeam_extract,
    run_multibeam_fit,
    run_pseudo_multibeam,
)



def main(argv: Optional[Sequence[str]] = None) -> None:
    if argv is None:
        import sys
        argv = list(sys.argv[1:])
    ap = argparse.ArgumentParser(description="Multi-beam sun_scan-compatible extraction, dry-run generation, and beam fitting.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_extract = sub.add_parser("extract", help="run sun_scan-compatible extraction for every selected stream")
    extract_mod.add_extract_arguments(ap_extract)

    ap_fit = sub.add_parser("fit", help="fit beam geometry from aggregated summary CSV files")
    fit_mod.add_fit_arguments(ap_fit)

    ap_pseudo = sub.add_parser("pseudo", help="generate pseudo multi-beam summaries from single-beam summary CSV files")
    pseudo_mod.add_arguments(ap_pseudo)

    ap_check = sub.add_parser("check-config", help="validate a converter-compatible spectrometer config")
    ap_check.add_argument("spectrometer_config")
    ap_check.add_argument("--stream-name", dest="stream_names", action="append", default=None)
    ap_check.add_argument("--out-csv", default=None)

    args = ap.parse_args(argv)
    if args.cmd == "extract":
        cfg = extract_mod.config_from_args(args, argv=argv)
        outputs = run_multibeam_extract(
            rawdata_path=cfg.input.rawdata_path,
            spectrometer_config=Path(args.spectrometer_config).expanduser().resolve(),
            base_config=cfg,
            outdir=cfg.report.outdir,
            run_id=args.run_id,
            stream_names=args.stream_names,
        )
        for key in ["all_summary_csv", "manifest_csv", "stream_table_csv", "config_snapshot_json"]:
            print(f"[done] {key}: {outputs[key]}")
        return
    if args.cmd == "fit":
        outputs = run_multibeam_fit(
            summary_paths=[Path(p).expanduser().resolve() for p in args.summary_csv],
            spectrometer_config=Path(args.spectrometer_config).expanduser().resolve(),
            outdir=Path(args.outdir).expanduser().resolve(),
            center_beam_id=args.center_beam_id,
            stream_names=args.stream_names,
            model=args.model,
            reference_angle_deg=args.reference_angle_deg,
            sigma_clip=(None if args.sigma_clip is None or float(args.sigma_clip) <= 0 else float(args.sigma_clip)),
            clip_iters=int(args.clip_iters),
            min_points_per_beam=int(args.min_points_per_beam),
            min_scans_per_beam=int(args.min_scans_per_beam),
        )
        for key, path in outputs.items():
            print(f"[done] {key}: {path}")
        return
    if args.cmd == "pseudo":
        outputs = run_pseudo_multibeam(
            singlebeam_summary_paths=[Path(p).expanduser().resolve() for p in args.singlebeam_summary_csv],
            spectrometer_config=Path(args.spectrometer_config).expanduser().resolve(),
            outdir=Path(args.outdir).expanduser().resolve(),
            stream_names=args.stream_names,
            noise_arcsec=float(args.noise_arcsec),
            seed=int(args.seed),
            tag=args.tag,
        )
        for key, path in outputs.items():
            print(f"[done] {key}: {path}")
        return
    if args.cmd == "check-config":
        result = check_spectrometer_config(Path(args.spectrometer_config).expanduser().resolve(), stream_names=args.stream_names)
        print(result.stream_table.to_string(index=False))
        if result.primary_streams:
            print(f"\nprimary_streams = {', '.join(result.primary_streams)}")
        if result.duplicate_beam_ids:
            print(f"duplicate_beam_ids = {', '.join(result.duplicate_beam_ids)}")
        if result.warnings:
            print("\nwarnings:")
            for w in result.warnings:
                print(f"  - {w}")
        if args.out_csv:
            out = Path(args.out_csv).expanduser().resolve()
            out.parent.mkdir(parents=True, exist_ok=True)
            result.stream_table.to_csv(out, index=False)
            print(f"\n[done] wrote {out}")
        return
    raise SystemExit(f"unsupported command: {args.cmd!r}")


if __name__ == "__main__":
    main()
