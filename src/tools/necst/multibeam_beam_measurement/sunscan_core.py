from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List

from . import sunscan_legacy_compat as compat
from .sunscan_config import SunScanAnalysisConfig
from .sunscan_io import build_dataframe_from_config
from .sunscan_report import SingleBeamAnalysisResult


def make_trim_params(config: SunScanAnalysisConfig) -> Dict[str, object]:
    return {
        "vfrac": float(config.trim.vfrac),
        "vmin": float(config.trim.vmin),
        "gap_fill": int(config.trim.gap_fill),
        "min_samples": int(config.trim.min_samples),
        "ratio_min": float(config.trim.ratio_min),
        "vpercentile": float(config.trim.vpercentile),
        "dominant_axis": bool(config.trim.dominant_axis),
        "steady_scan": bool(config.trim.steady_scan),
        "use_on_only": bool(config.trim.use_on_only),
        "profile_xlim_deg": float(config.profile.profile_xlim_deg),
        "xwin_factor": float(config.trim.xwin_factor),
        "cross_offset_max_deg": float(config.trim.cross_offset_max_deg),
        "speed_min_deg_s": float(config.trim.speed_min_deg_s),
        "steady_cv_max": float(config.trim.steady_cv_max),
    }


def make_ripple_policy(config: SunScanAnalysisConfig) -> Dict[str, object]:
    args = SimpleNamespace(
        ripple_preset=config.ripple.preset,
        ripple_model=config.ripple.model,
        ripple_target_hz=config.ripple.target_hz,
        ripple_search_hz=config.ripple.search_hz,
        ripple_bw_hz=config.ripple.bw_hz,
        ripple_max_harm=config.ripple.max_harm,
        ripple_order=config.ripple.order,
        ripple_notch_pass=config.ripple.notch_passes,
        ripple_trend_win_sec=config.ripple.trend_win_sec,
        ripple_resample_dt_sec=config.ripple.resample_dt_sec,
        ripple_eval_band_hz=config.ripple.eval_band_hz,
    )
    return compat.make_ripple_policy_from_args(args)


def _annotate_results_metadata(results: Dict[int, dict], config: SunScanAnalysisConfig, ycol: str, df_attrs: Dict[str, object] | None = None) -> None:
    tag = config.resolved_tag()
    df_attrs = dict(df_attrs or {})
    for sid, row in list(results.items()):
        if not isinstance(sid, int) or not isinstance(row, dict):
            continue
        row.setdefault("data_tag", tag)
        row.setdefault("y_axis", ycol)
        row.setdefault("spec_time_basis", df_attrs.get("spec_time_basis"))
        row.setdefault("spec_time_suffix", df_attrs.get("spec_time_suffix"))
        row.setdefault("spec_time_fallback_field", df_attrs.get("spec_time_fallback_field"))
        row.setdefault("spec_time_example", df_attrs.get("spec_time_example"))
        row.setdefault("azel_source", config.input.azel_source)
        resolved_azel_apply = df_attrs.get("azel_correction_apply")
        if resolved_azel_apply is None:
            resolved_azel_apply = config.input.azel_correction_apply
        row.setdefault("azel_correction_apply", resolved_azel_apply)
        resolved_altaz_apply = df_attrs.get("azel_correction_apply")
        if resolved_altaz_apply is None:
            resolved_altaz_apply = config.input.altaz_apply
        row.setdefault("altaz_apply", resolved_altaz_apply)
        row.setdefault("spectrometer_time_offset_sec", float(config.input.spectrometer_time_offset_sec))
        row.setdefault("encoder_shift_sec", float(config.input.encoder_shift_sec))
        row.setdefault("encoder_az_time_offset_sec", float(config.input.encoder_az_time_offset_sec))
        row.setdefault("encoder_el_time_offset_sec", float(config.input.encoder_el_time_offset_sec))
        row.setdefault("encoder_vavg_sec", float(config.input.encoder_vavg_sec))
        row.setdefault("chopper_wheel", bool(config.calibration.chopper_wheel))
        row.setdefault("ripple_remove", bool(config.ripple.enabled))
        row.setdefault("ripple_preset", str(config.ripple.preset))
        row.setdefault("edge_fit_win_deg", float(config.edge_fit.fit_win_deg))
        row.setdefault("edge_fit_threshold", float(config.edge_fit.fit_threshold))
        row.setdefault("hpbw_init_arcsec", float(config.edge_fit.hpbw_init_arcsec))
        row.setdefault("trim_scan", bool(config.trim.enabled))
        row.setdefault("profile_xlim_deg", float(config.profile.profile_xlim_deg))
        row.setdefault("db_namespace", str(config.input.db_namespace))
        row.setdefault("telescope", str(config.input.telescope))
        row.setdefault("tel_loaddata", str(config.input.tel_loaddata))
        row.setdefault("planet", str(config.input.planet))
        row.setdefault("encoder_table_resolved", df_attrs.get("encoder_table_resolved"))
        row.setdefault("altaz_table_resolved", df_attrs.get("altaz_table_resolved"))
        row.setdefault("weather_inside_table_resolved", df_attrs.get("weather_inside_table_resolved"))
        row.setdefault("weather_outside_table_resolved", df_attrs.get("weather_outside_table_resolved"))
        row.setdefault("encoder_time_col", df_attrs.get("encoder_time_col"))
        row.setdefault("altaz_time_col", df_attrs.get("altaz_time_col"))
        row.setdefault("weather_inside_time_col", df_attrs.get("weather_inside_time_col"))
        row.setdefault("weather_outside_time_col", df_attrs.get("weather_outside_time_col"))
        row.setdefault("tamb_source", df_attrs.get("tamb_source"))
        row.setdefault("outside_meteo_source", df_attrs.get("outside_meteo_source"))
        row.setdefault("outside_default_used_temperature", df_attrs.get("outside_default_used_temperature"))
        row.setdefault("outside_default_used_pressure", df_attrs.get("outside_default_used_pressure"))
        row.setdefault("outside_default_used_humidity", df_attrs.get("outside_default_used_humidity"))


def run_debug_plots(config: SunScanAnalysisConfig, az_scans: Dict[int, object], el_scans: Dict[int, object], trim_params: Dict[str, object]) -> List[Path]:
    outdir = Path(config.report.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tag = config.resolved_tag()
    out_paths: List[Path] = []

    for i, seg in az_scans.items():
        try:
            if bool(config.report.debug_plot):
                p1 = outdir / f"AZ_dxdy_{i}_{tag}.png"
                p2 = outdir / f"AZ_mount_azel_{i}_{tag}.png"
                compat.plot_dxdy_simple(seg, title=f"AZ scan {i}: daz/d_el and Δ per sample", out_png=p1, trim_enabled=bool(config.trim.enabled), trim_params=trim_params, which_xcol="daz")
                compat.plot_mount_azel_simple(seg, title=f"AZ scan {i}: Az_true/El_true and Δ per sample", out_png=p2)
                out_paths.extend([p1, p2])
        except Exception as exc:
            print(f"[error] AZ debug plot scan {i}: {exc}")
            if not config.runtime.continue_on_error:
                raise

    for i, seg in el_scans.items():
        try:
            if bool(config.report.debug_plot):
                p1 = outdir / f"EL_dxdy_{i}_{tag}.png"
                p2 = outdir / f"EL_mount_azel_{i}_{tag}.png"
                compat.plot_dxdy_simple(seg, title=f"EL scan {i}: daz/d_el and Δ per sample", out_png=p1, trim_enabled=bool(config.trim.enabled), trim_params=trim_params, which_xcol="d_el")
                compat.plot_mount_azel_simple(seg, title=f"EL scan {i}: Az_true/El_true and Δ per sample", out_png=p2)
                out_paths.extend([p1, p2])
        except Exception as exc:
            print(f"[error] EL debug plot scan {i}: {exc}")
            if not config.runtime.continue_on_error:
                raise
    return out_paths


def analyze_single_stream(config: SunScanAnalysisConfig) -> SingleBeamAnalysisResult:
    df, az_scans, el_scans = build_dataframe_from_config(config)
    trim_params = make_trim_params(config)
    ripple_policy = make_ripple_policy(config)
    ycol = "ta_star" if bool(config.calibration.chopper_wheel) else "tp1"
    ytitle = "Ta*" if bool(config.calibration.chopper_wheel) else "tp1"

    debug_pngs = run_debug_plots(config, az_scans, el_scans, trim_params)

    result = SingleBeamAnalysisResult(
        config=config,
        df=df,
        az_scans=az_scans,
        el_scans=el_scans,
        ycol=ycol,
        ytitle=ytitle,
        trim_params=trim_params,
        ripple_policy=ripple_policy,
    )
    result.report_paths.debug_pngs = debug_pngs

    if bool(config.edge_fit.enabled):
        try:
            az_deriv, el_deriv, az_fit, el_fit, results = compat.run_edge_fits_and_summarize(
                az_scans=az_scans,
                el_scans=el_scans,
                ycol=ycol,
                trim_enabled=bool(config.trim.enabled),
                trim_params=trim_params,
                strict_deriv=bool(config.edge_fit.strict_deriv),
                ripple_enabled=bool(config.ripple.enabled),
                ripple_policy=ripple_policy,
                profile_xlim_deg=float(config.profile.profile_xlim_deg),
                fit_win_deg=float(config.edge_fit.fit_win_deg),
                fit_threshold=float(config.edge_fit.fit_threshold),
                hpbw_init_arcsec=float(config.edge_fit.hpbw_init_arcsec),
            )
            scan_ids = results.pop("__scan_ids__", sorted(set(az_scans.keys()) | set(el_scans.keys()))) if isinstance(results, dict) else []
            az_deriv_raw = results.pop("__az_deriv_raw__", {}) if isinstance(results, dict) else {}
            el_deriv_raw = results.pop("__el_deriv_raw__", {}) if isinstance(results, dict) else {}
            az_profile_raw = results.pop("__az_profile_raw__", {}) if isinstance(results, dict) else {}
            el_profile_raw = results.pop("__el_profile_raw__", {}) if isinstance(results, dict) else {}
            az_profile = results.pop("__az_profile__", {}) if isinstance(results, dict) else {}
            el_profile = results.pop("__el_profile__", {}) if isinstance(results, dict) else {}
            _annotate_results_metadata(results, config, ycol, df_attrs=getattr(df, "attrs", {}))
            result.az_deriv = az_deriv
            result.el_deriv = el_deriv
            result.az_fit = az_fit
            result.el_fit = el_fit
            result.results = results
            result.scan_ids = list(scan_ids)
            result.az_deriv_raw = az_deriv_raw
            result.el_deriv_raw = el_deriv_raw
            result.az_profile_raw = az_profile_raw
            result.el_profile_raw = el_profile_raw
            result.az_profile = az_profile
            result.el_profile = el_profile
        except Exception as exc:
            print(f"[warn] edge fitting failed: {exc}")

    return result
