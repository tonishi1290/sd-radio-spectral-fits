from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Sequence

from . import sunscan_extract_multibeam as extract_mod
from . import sunscan_fit_multibeam as fit_mod
from . import synthetic_multibeam as pseudo_mod


def _load_public_api():
    """Load runtime-heavy public API only for commands that need real analysis."""

    from .public_api import (  # type: ignore
        check_spectrometer_config,
        run_multibeam_extract,
        run_multibeam_extract_many,
        run_multibeam_fit,
        run_pseudo_multibeam,
    )
    return {
        "check_spectrometer_config": check_spectrometer_config,
        "run_multibeam_extract": run_multibeam_extract,
        "run_multibeam_extract_many": run_multibeam_extract_many,
        "run_multibeam_fit": run_multibeam_fit,
        "run_pseudo_multibeam": run_pseudo_multibeam,
    }



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
    ap_check.add_argument("--config-loader", default="legacy", choices=["legacy", "adapter"], help="Loader used for spectrometer_config validation")
    ap_check.add_argument("--sunscan-analysis-config", "--sunscan-analysis", default=None, help="Standalone sunscan_analysis TOML to apply before validation")
    ap_check.add_argument("--analysis-stream-selection", default=None, action="append", help="Standalone analysis_stream_selection TOML to apply before validation. May be repeated.")
    ap_check.add_argument("--out-csv", default=None)

    ap_inspect = sub.add_parser("inspect-sidecars", help="inspect DB-embedded spectral-recording sidecars and exit")
    ap_inspect.add_argument("rawdata", nargs="+", help="One or more RawData directories to inspect")

    args = ap.parse_args(argv)
    if args.cmd == "inspect-sidecars":
        rawdata_paths = [Path(p).expanduser().resolve() for p in args.rawdata]
        print(json.dumps(extract_mod.inspect_rawdata_sidecars(rawdata_paths), ensure_ascii=False, indent=2))
        return
    if args.cmd == "extract":
        rawdata_paths = [Path(p).expanduser().resolve() for p in args.rawdata]
        if bool(getattr(args, "inspect_sidecars", False)):
            print(json.dumps(extract_mod.inspect_rawdata_sidecars(rawdata_paths), ensure_ascii=False, indent=2))
            return
        cfg = extract_mod.config_from_args(args, argv=argv)
        api = _load_public_api()
        if len(rawdata_paths) == 1:
            outputs = api["run_multibeam_extract"](
                rawdata_path=cfg.input.rawdata_path,
                spectrometer_config=Path(args.spectrometer_config).expanduser().resolve() if getattr(args, "spectrometer_config", None) else None,
                base_config=cfg,
                outdir=cfg.report.outdir,
                run_id=args.run_id,
                stream_names=args.stream_names,
                config_loader=getattr(args, "config_loader", "legacy"),
                sunscan_analysis_config=Path(args.sunscan_analysis_config).expanduser().resolve() if getattr(args, "sunscan_analysis_config", None) else None,
                analysis_stream_selection=[Path(p).expanduser().resolve() for p in list(getattr(args, "analysis_stream_selection", None) or [])],
                spectral_recording_snapshot=Path(args.spectral_recording_snapshot).expanduser().resolve() if getattr(args, "spectral_recording_snapshot", None) else None,
                beam_model_path=Path(args.beam_model).expanduser().resolve() if getattr(args, "beam_model", None) else None,
                allow_beam_model_override=bool(getattr(args, "allow_beam_model_override", False)),
            )
            for key in ["all_summary_csv", "manifest_csv", "stream_table_csv", "config_snapshot_json"]:
                print(f"[done] {key}: {outputs[key]}")
            return
        outputs = api["run_multibeam_extract_many"](
            rawdata_paths=rawdata_paths,
            spectrometer_config=Path(args.spectrometer_config).expanduser().resolve() if getattr(args, "spectrometer_config", None) else None,
            base_config=cfg,
            outdir=cfg.report.outdir,
            run_ids=None,
            stream_names=args.stream_names,
            merged_tag=args.run_id,
            config_loader=getattr(args, "config_loader", "legacy"),
            sunscan_analysis_config=Path(args.sunscan_analysis_config).expanduser().resolve() if getattr(args, "sunscan_analysis_config", None) else None,
            analysis_stream_selection=[Path(p).expanduser().resolve() for p in list(getattr(args, "analysis_stream_selection", None) or [])],
            spectral_recording_snapshot=Path(args.spectral_recording_snapshot).expanduser().resolve() if getattr(args, "spectral_recording_snapshot", None) else None,
            beam_model_path=Path(args.beam_model).expanduser().resolve() if getattr(args, "beam_model", None) else None,
            allow_beam_model_override=bool(getattr(args, "allow_beam_model_override", False)),
        )
        for key in ["all_summary_csv", "manifest_csv", "run_table_csv"]:
            print(f"[done] {key}: {outputs[key]}")
        return
    if args.cmd == "fit":
        outputs = _load_public_api()["run_multibeam_fit"](
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
            config_loader=getattr(args, "config_loader", "legacy"),
            analysis_stream_selection=[Path(p).expanduser().resolve() for p in list(getattr(args, "analysis_stream_selection", None) or [])],
        )
        for key, path in outputs.items():
            print(f"[done] {key}: {path}")
        return
    if args.cmd == "pseudo":
        outputs = _load_public_api()["run_pseudo_multibeam"](
            singlebeam_summary_paths=[Path(p).expanduser().resolve() for p in args.singlebeam_summary_csv],
            spectrometer_config=Path(args.spectrometer_config).expanduser().resolve(),
            outdir=Path(args.outdir).expanduser().resolve(),
            stream_names=args.stream_names,
            noise_arcsec=float(args.noise_arcsec),
            seed=int(args.seed),
            tag=args.tag,
            rep_el_degs=args.rep_el_degs,
            config_loader=getattr(args, "config_loader", "legacy"),
            analysis_stream_selection=[Path(p).expanduser().resolve() for p in list(getattr(args, "analysis_stream_selection", None) or [])],
        )
        for key, path in outputs.items():
            print(f"[done] {key}: {path}")
        return
    if args.cmd == "check-config":
        result = _load_public_api()["check_spectrometer_config"](
            Path(args.spectrometer_config).expanduser().resolve(),
            stream_names=args.stream_names,
            config_loader=getattr(args, "config_loader", "legacy"),
            sunscan_analysis_config=Path(args.sunscan_analysis_config).expanduser().resolve() if getattr(args, "sunscan_analysis_config", None) else None,
            analysis_stream_selection=[Path(p).expanduser().resolve() for p in list(getattr(args, "analysis_stream_selection", None) or [])],
        )
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
