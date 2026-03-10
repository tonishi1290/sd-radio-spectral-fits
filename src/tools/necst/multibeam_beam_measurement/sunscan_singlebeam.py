from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from .public_api import run_singlebeam
from .sunscan_config import SunScanAnalysisConfig
from .sunscan_report import build_summary_lines



def add_singlebeam_arguments(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("rawdata", help="RawData directory containing the NECST database")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--spectral-name", default="xffts-board1", help="Spectral stream name")

    ap.add_argument("--azel-source", choices=["encoder", "altaz"], default="encoder")
    ap.add_argument("--altaz-apply", choices=["none", "minus", "plus"], default="none")
    ap.add_argument("--encoder-shift-sec", type=float, default=0.0)
    ap.add_argument("--encoder-vavg-sec", type=float, default=0.0)
    ap.add_argument("--no-chopper-wheel", dest="chopper_wheel", action="store_false")
    ap.add_argument("--tamb-k", type=float, default=None)
    ap.add_argument("--chopper-win-sec", type=float, default=5.0)
    ap.add_argument("--chopper-stat", choices=["median", "mean"], default="median")
    ap.set_defaults(chopper_wheel=True)

    ap.add_argument("--profile-xlim-deg", type=float, default=1.0)
    ap.add_argument("--ripple-no-remove", dest="ripple_remove", action="store_false")
    ap.add_argument("--ripple-preset", choices=["auto", "safe", "normal", "strong"], default="auto")
    ap.add_argument("--ripple-model", choices=["auto", "add", "mul"], default="auto")
    ap.add_argument("--ripple-target-hz", type=float, default=1.2)
    ap.add_argument("--ripple-search-hz", type=float, default=0.3)
    ap.add_argument("--ripple-bw-hz", type=float, default=None)
    ap.add_argument("--ripple-max-harm", type=int, default=None)
    ap.add_argument("--ripple-order", type=int, default=None)
    ap.add_argument("--ripple-notch-pass", type=int, default=None)
    ap.add_argument("--ripple-trend-win-sec", type=float, default=None)
    ap.add_argument("--ripple-resample-dt-sec", type=float, default=None)
    ap.add_argument("--ripple-eval-band-hz", type=float, default=None)
    ap.set_defaults(ripple_remove=True)

    ap.add_argument("--no-edge-fit", dest="edge_fit", action="store_false")
    ap.add_argument("--edge-fit-win-deg", type=float, default=0.15)
    ap.add_argument("--edge-fit-threshold", type=float, default=0.20)
    ap.add_argument("--hpbw-init-arcsec", type=float, default=324.0)
    ap.add_argument("--edge-fit-plot-max-scans", type=int, default=3)
    ap.set_defaults(edge_fit=True)

    ap.add_argument("--no-trim-scan", dest="trim_scan", action="store_false")
    ap.add_argument("--trim-vfrac", type=float, default=0.20)
    ap.add_argument("--trim-vmin", type=float, default=1e-4)
    ap.add_argument("--trim-gap", type=int, default=10)
    ap.add_argument("--trim-min-samples", type=int, default=100)
    ap.add_argument("--trim-dominant-axis", dest="trim_dominant_axis", action="store_true")
    ap.add_argument("--trim-no-dominant-axis", dest="trim_dominant_axis", action="store_false")
    ap.add_argument("--trim-axis-ratio-min", type=float, default=3.0)
    ap.add_argument("--trim-vpercentile", type=float, default=95.0)
    ap.add_argument("--trim-no-steady-scan", dest="trim_steady_scan", action="store_false")
    ap.add_argument("--trim-include-hotoff", dest="trim_use_on_only", action="store_false")
    ap.add_argument("--trim-scan-speed-min-arcsec", type=float, default=20.0)
    ap.add_argument("--trim-xwin-factor", type=float, default=1.2)
    ap.add_argument("--trim-cross-offset-max-deg", type=float, default=0.5)
    ap.add_argument("--trim-steady-cv-max", type=float, default=0.8)
    ap.set_defaults(trim_scan=True)
    ap.set_defaults(trim_dominant_axis=True)
    ap.set_defaults(trim_steady_scan=True)
    ap.set_defaults(trim_use_on_only=True)

    ap.add_argument("--strict-deriv", dest="strict_deriv", action="store_true")
    ap.add_argument("--no-strict-deriv", dest="strict_deriv", action="store_false")
    ap.set_defaults(strict_deriv=True)

    ap.add_argument("--continue-on-error", action="store_true", default=False)
    ap.add_argument("--debug-plot", action="store_true", default=False)
    ap.add_argument("--dish-diameter-m", type=float, default=1.85)
    ap.add_argument("--hpbw-factor", type=float, default=1.2)



def config_from_args(args: argparse.Namespace) -> SunScanAnalysisConfig:
    cfg = SunScanAnalysisConfig.default(
        rawdata_path=Path(args.rawdata).expanduser().resolve(),
        spectral_name=str(args.spectral_name),
        outdir=Path(args.outdir).expanduser().resolve(),
    )
    cfg.input.azel_source = str(args.azel_source)
    cfg.input.altaz_apply = str(args.altaz_apply)
    cfg.input.encoder_shift_sec = float(args.encoder_shift_sec)
    cfg.input.encoder_vavg_sec = float(args.encoder_vavg_sec)
    cfg.calibration.chopper_wheel = bool(args.chopper_wheel)
    cfg.calibration.tamb_k = args.tamb_k
    cfg.calibration.chopper_win_sec = float(args.chopper_win_sec)
    cfg.calibration.chopper_stat = str(args.chopper_stat)
    cfg.profile.profile_xlim_deg = float(args.profile_xlim_deg)
    cfg.ripple.enabled = bool(args.ripple_remove)
    cfg.ripple.preset = str(args.ripple_preset)
    cfg.ripple.model = str(args.ripple_model)
    cfg.ripple.target_hz = float(args.ripple_target_hz)
    cfg.ripple.search_hz = float(args.ripple_search_hz)
    cfg.ripple.bw_hz = args.ripple_bw_hz
    cfg.ripple.max_harm = args.ripple_max_harm
    cfg.ripple.order = args.ripple_order
    cfg.ripple.notch_passes = args.ripple_notch_pass
    cfg.ripple.trend_win_sec = args.ripple_trend_win_sec
    cfg.ripple.resample_dt_sec = args.ripple_resample_dt_sec
    cfg.ripple.eval_band_hz = args.ripple_eval_band_hz
    cfg.edge_fit.enabled = bool(args.edge_fit)
    cfg.edge_fit.fit_win_deg = float(args.edge_fit_win_deg)
    cfg.edge_fit.fit_threshold = float(args.edge_fit_threshold)
    cfg.edge_fit.hpbw_init_arcsec = float(args.hpbw_init_arcsec)
    cfg.edge_fit.strict_deriv = bool(args.strict_deriv)
    cfg.report.edge_fit_plot_max_scans = int(args.edge_fit_plot_max_scans)
    cfg.trim.enabled = bool(args.trim_scan)
    cfg.trim.vfrac = float(args.trim_vfrac)
    cfg.trim.vmin = float(args.trim_vmin)
    cfg.trim.gap_fill = int(args.trim_gap)
    cfg.trim.min_samples = int(args.trim_min_samples)
    cfg.trim.dominant_axis = bool(args.trim_dominant_axis)
    cfg.trim.ratio_min = float(args.trim_axis_ratio_min)
    cfg.trim.vpercentile = float(args.trim_vpercentile)
    cfg.trim.steady_scan = bool(args.trim_steady_scan)
    cfg.trim.use_on_only = bool(args.trim_use_on_only)
    cfg.trim.xwin_factor = float(args.trim_xwin_factor)
    cfg.trim.cross_offset_max_deg = float(args.trim_cross_offset_max_deg)
    cfg.trim.speed_min_deg_s = float(args.trim_scan_speed_min_arcsec) / 3600.0
    cfg.trim.steady_cv_max = float(args.trim_steady_cv_max)
    cfg.runtime.continue_on_error = bool(args.continue_on_error)
    cfg.report.debug_plot = bool(args.debug_plot)
    cfg.dish_diameter_m = float(args.dish_diameter_m)
    cfg.hpbw_factor = float(args.hpbw_factor)
    return cfg



def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Run single-stream sun_scan-compatible analysis via the shared package API.")
    add_singlebeam_arguments(ap)
    args = ap.parse_args(argv)
    cfg = config_from_args(args)
    result = run_singlebeam(cfg)

    if bool(result.config.edge_fit.enabled) and result.results:
        lines = build_summary_lines(
            result.config.resolved_tag(),
            result.ytitle,
            result.config,
            result.scan_ids,
            result.results,
        )
        print("\n".join(lines))

    if result.report_paths.derivative_pngs:
        page_pat = f"sun_scan_derivative_fits_{result.config.resolved_tag()}_pXXX.png"
        print(f"[fit] wrote {len(result.report_paths.derivative_pngs)} page(s) of {page_pat}")
    if result.report_paths.summary_csv:
        print(f"[fit] wrote {result.report_paths.summary_csv.name}")

    print("[done]")
    print(f"  outputs in: {result.config.report.outdir}")
    print(f"  df rows: {len(result.df)}  az_scans: {len(result.az_scans)}  el_scans: {len(result.el_scans)}")


if __name__ == "__main__":
    main()
