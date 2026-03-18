from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, Optional, Sequence
import math
import sys

from .config_io import load_spectrometer_config, resolve_stream_usage_policy, restfreq_hz_for_stream, stream_extras_by_name
from .public_api import run_singlebeam
from .sunscan_config import SunScanAnalysisConfig
from .sunscan_report import build_summary_lines


def _argv_has_option(argv: Optional[Sequence[str]], *opts: str) -> bool:
    sargv = [str(x) for x in (argv or [])]
    for a in sargv:
        for opt in opts:
            if a == opt or a.startswith(opt + "="):
                return True
    return False


def _choose_float_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: float = 0.0) -> float:
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt) or (global_cfg.get(global_key) is None):
        return float(getattr(args, attr, default))
    return float(global_cfg.get(global_key))


def _choose_str_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: str) -> str:
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt) or (global_cfg.get(global_key) is None):
        return str(getattr(args, attr, default))
    return str(global_cfg.get(global_key))

def _choose_optional_str_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: Optional[str] = None) -> Optional[str]:
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt) or (global_cfg.get(global_key) is None):
        value = getattr(args, attr, default)
    else:
        value = global_cfg.get(global_key)
    if value is None:
        return default
    s = str(value).strip()
    return s if s else default


def _normalize_legacy_azel_correction_apply(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    s = str(value).strip().lower()
    if not s:
        return None
    if s == "minus":
        return "subtract"
    if s == "plus":
        return "add"
    if s in {"subtract", "add", "none"}:
        return s
    raise ValueError(f"unsupported azel correction mode: {value!r}")


def _resolve_azel_correction_apply(args: argparse.Namespace, argv: Optional[Sequence[str]], global_cfg: Dict[str, Any], azel_source: str) -> str:
    raw = None
    if _argv_has_option(argv, "--azel-correction-apply"):
        raw = getattr(args, "azel_correction_apply", None)
    elif _argv_has_option(argv, "--altaz-apply"):
        raw = getattr(args, "altaz_apply", None)
    elif global_cfg.get("azel_correction_apply") is not None:
        raw = global_cfg.get("azel_correction_apply")
    elif global_cfg.get("altaz_apply") is not None:
        raw = global_cfg.get("altaz_apply")
    norm = _normalize_legacy_azel_correction_apply(raw)
    if norm is not None:
        return norm
    return "subtract" if str(azel_source).strip().lower() == "encoder" else "none"


def _estimate_hpbw_init_arcsec(restfreq_hz: Optional[float], dish_diameter_m: float, hpbw_factor: float, fallback_arcsec: float) -> float:
    if restfreq_hz is None or not math.isfinite(restfreq_hz) or restfreq_hz <= 0:
        return float(fallback_arcsec)
    c_m_s = 299792458.0
    return float(hpbw_factor) * (c_m_s / float(restfreq_hz)) / float(dish_diameter_m) * 206265.0


def _select_stream_from_config(config_dict: Dict[str, Any], *, stream_name: Optional[str], spectral_name: Optional[str]) -> Any:
    streams = list(config_dict.get("streams", []) or [])
    if not streams:
        raise ValueError("spectrometer config contains no streams")
    if stream_name:
        wanted = str(stream_name)
        matches = [s for s in streams if str(getattr(s, "name", "")) == wanted]
        if not matches:
            raise ValueError(f"stream_name={wanted!r} was not found in spectrometer config")
        if len(matches) != 1:
            raise ValueError(f"stream_name={wanted!r} matched multiple streams")
        return matches[0]
    if spectral_name:
        wanted = str(spectral_name)
        matches = [s for s in streams if str(getattr(s, "name", "")) == wanted or str(getattr(s, "db_stream_name", "")) == wanted]
        if len(matches) == 1:
            return matches[0]
    if len(streams) == 1:
        return streams[0]
    names = [str(getattr(s, "name", "")) for s in streams]
    raise ValueError(f"multiple streams are defined in spectrometer config; specify --stream-name. available={names}")



def add_singlebeam_arguments(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("rawdata", help="RawData directory containing the NECST database")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--spectrometer-config", default=None, help="converter-compatible spectrometer config TOML")
    ap.add_argument("--stream-name", default=None, help="Select one stream from --spectrometer-config")
    ap.add_argument("--db-namespace", default="necst", help="Database namespace/prefix used in NECST DB table names")
    ap.add_argument("--telescope", default="OMU1P85M", help="Telescope name used in NECST DB table names")
    ap.add_argument("--tel-loaddata", default="OMU1p85m", help="Telescope name passed to nercst loaddb()")
    ap.add_argument("--planet", default="sun", help="Target body name passed to astropy.get_body()")
    ap.add_argument("--spectral-name", default="xffts-board1", help="Spectral stream name")

    ap.add_argument("--azel-source", choices=["encoder", "altaz"], default="encoder")
    ap.add_argument("--azel-correction-apply", choices=["none", "subtract", "add"], default=None)
    ap.add_argument("--altaz-apply", choices=["none", "minus", "plus"], default=None, help="Legacy alias for --azel-correction-apply")
    ap.add_argument("--spectrometer-time-offset-sec", type=float, default=0.0)
    ap.add_argument("--encoder-shift-sec", type=float, default=0.0)
    ap.add_argument("--encoder-az-time-offset-sec", type=float, default=0.0)
    ap.add_argument("--encoder-el-time-offset-sec", type=float, default=0.0)
    ap.add_argument("--encoder-vavg-sec", type=float, default=0.0)
    ap.add_argument("--encoder-table", default=None)
    ap.add_argument("--encoder-table-suffix", default=None)
    ap.add_argument("--altaz-table", default=None)
    ap.add_argument("--altaz-table-suffix", default=None)
    ap.add_argument("--encoder-time-col", default="time")
    ap.add_argument("--altaz-time-col", default="time")
    ap.add_argument("--weather-inside-table", default=None)
    ap.add_argument("--weather-inside-table-suffix", default=None)
    ap.add_argument("--weather-inside-time-col", default=None)
    ap.add_argument("--weather-outside-table", default=None)
    ap.add_argument("--weather-outside-table-suffix", default=None)
    ap.add_argument("--weather-outside-time-col", default=None)
    ap.add_argument("--weather-table", default=None, help="Legacy alias applied to both inside/outside weather tables")
    ap.add_argument("--weather-time-col", default=None, help="Legacy alias applied to both inside/outside weather time columns")
    ap.add_argument("--no-chopper-wheel", dest="chopper_wheel", action="store_false")
    ap.add_argument("--tamb-k", type=float, default=None)
    ap.add_argument("--chopper-win-sec", type=float, default=5.0)
    ap.add_argument("--chopper-stat", choices=["median", "mean"], default="median")
    ap.set_defaults(chopper_wheel=True)
    ap.add_argument("--tamb-default-k", type=float, default=300.0)
    ap.add_argument("--tamb-min-k", type=float, default=250.0)
    ap.add_argument("--tamb-max-k", type=float, default=330.0)
    ap.add_argument("--outside-default-temperature-c", type=float, default=0.0)
    ap.add_argument("--outside-default-pressure-hpa", type=float, default=760.0)
    ap.add_argument("--outside-default-humidity-pct", type=float, default=30.0)
    ap.add_argument("--outside-temperature-min-c", type=float, default=-50.0)
    ap.add_argument("--outside-temperature-max-c", type=float, default=50.0)
    ap.add_argument("--outside-pressure-min-hpa", type=float, default=400.0)
    ap.add_argument("--outside-pressure-max-hpa", type=float, default=1100.0)
    ap.add_argument("--outside-humidity-min-pct", type=float, default=0.0)
    ap.add_argument("--outside-humidity-max-pct", type=float, default=100.0)

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



def config_from_args(args: argparse.Namespace, argv: Optional[Sequence[str]] = None) -> SunScanAnalysisConfig:
    cfg = SunScanAnalysisConfig.default(
        rawdata_path=Path(args.rawdata).expanduser().resolve(),
        spectral_name=str(args.spectral_name),
        outdir=Path(args.outdir).expanduser().resolve(),
    )
    global_cfg: Dict[str, Any] = {}
    selected_stream = None
    if args.spectrometer_config:
        config_path = Path(args.spectrometer_config).expanduser().resolve()
        config_dict = load_spectrometer_config(config_path)
        global_cfg = dict(config_dict.get("global", {}) or {})
        selected_stream = _select_stream_from_config(config_dict, stream_name=getattr(args, "stream_name", None), spectral_name=getattr(args, "spectral_name", None))
        explicit_stream_override = bool(getattr(args, "stream_name", None)) or _argv_has_option(argv, "--spectral-name")
        if (selected_stream is not None) and (not explicit_stream_override):
            usage_policy = resolve_stream_usage_policy(stream_extras_by_name(config_path).get(str(getattr(selected_stream, "name", "")), {}))
            if not bool(usage_policy.get("use_for_sunscan", True)):
                raise ValueError(
                    f"selected stream {getattr(selected_stream, 'name', None)!r} is disabled for sunscan by enabled/use_for_sunscan in the spectrometer config. "
                    "Pass --stream-name explicitly if you intentionally want to override this for a one-off run."
                )

    cfg.input.db_namespace = _choose_str_setting(args, argv, "--db-namespace", "db_namespace", global_cfg, "db_namespace", "necst")
    cfg.input.telescope = _choose_str_setting(args, argv, "--telescope", "telescope", global_cfg, "telescope", "OMU1P85M")
    cfg.input.tel_loaddata = _choose_str_setting(args, argv, "--tel-loaddata", "tel_loaddata", global_cfg, "tel_loaddata", "OMU1p85m")
    cfg.input.planet = _choose_str_setting(args, argv, "--planet", "planet", global_cfg, "planet", "sun")
    cfg.input.azel_source = _choose_str_setting(args, argv, "--azel-source", "azel_source", global_cfg, "azel_source", "encoder")
    cfg.input.azel_correction_apply = _resolve_azel_correction_apply(args, argv, global_cfg, cfg.input.azel_source)
    cfg.input.altaz_apply = cfg.input.azel_correction_apply
    cfg.input.encoder_table = _choose_optional_str_setting(args, argv, "--encoder-table", "encoder_table", global_cfg, "encoder_table", None)
    cfg.input.encoder_table_suffix = _choose_optional_str_setting(args, argv, "--encoder-table-suffix", "encoder_table_suffix", global_cfg, "encoder_table_suffix", "ctrl-antenna-encoder") or "ctrl-antenna-encoder"
    cfg.input.altaz_table = _choose_optional_str_setting(args, argv, "--altaz-table", "altaz_table", global_cfg, "altaz_table", None)
    cfg.input.altaz_table_suffix = _choose_optional_str_setting(args, argv, "--altaz-table-suffix", "altaz_table_suffix", global_cfg, "altaz_table_suffix", "ctrl-antenna-altaz") or "ctrl-antenna-altaz"
    cfg.input.encoder_time_col = _choose_str_setting(args, argv, "--encoder-time-col", "encoder_time_col", global_cfg, "encoder_time_col", "time")
    cfg.input.altaz_time_col = _choose_str_setting(args, argv, "--altaz-time-col", "altaz_time_col", global_cfg, "altaz_time_col", "time")
    cfg.input.spectrometer_time_offset_sec = _choose_float_setting(args, argv, "--spectrometer-time-offset-sec", "spectrometer_time_offset_sec", global_cfg, "spectrometer_time_offset_sec", 0.0)
    cfg.input.encoder_shift_sec = _choose_float_setting(args, argv, "--encoder-shift-sec", "encoder_shift_sec", global_cfg, "encoder_shift_sec", 0.0)
    cfg.input.encoder_az_time_offset_sec = _choose_float_setting(args, argv, "--encoder-az-time-offset-sec", "encoder_az_time_offset_sec", global_cfg, "encoder_az_time_offset_sec", 0.0)
    cfg.input.encoder_el_time_offset_sec = _choose_float_setting(args, argv, "--encoder-el-time-offset-sec", "encoder_el_time_offset_sec", global_cfg, "encoder_el_time_offset_sec", 0.0)
    cfg.input.encoder_vavg_sec = _choose_float_setting(args, argv, "--encoder-vavg-sec", "encoder_vavg_sec", global_cfg, "encoder_vavg_sec", 0.0)
    if selected_stream is not None:
        cfg.input.spectral_name = str(getattr(selected_stream, "db_stream_name", None) or getattr(selected_stream, "name", cfg.input.spectral_name))
        cfg.beam_override.stream_name = str(getattr(selected_stream, "name", cfg.beam_override.stream_name or "")) or cfg.beam_override.stream_name
        cfg.beam_override.beam_id = str(getattr(getattr(selected_stream, "beam", None), "beam_id", cfg.beam_override.beam_id or "")) or cfg.beam_override.beam_id
        cfg.beam_override.restfreq_hz = restfreq_hz_for_stream(selected_stream)
        cfg.beam_override.polariza = str(getattr(selected_stream, "polariza", cfg.beam_override.polariza or "")) or cfg.beam_override.polariza
        cfg.beam_override.fdnum = int(getattr(selected_stream, "fdnum", cfg.beam_override.fdnum or 0))
        cfg.beam_override.ifnum = int(getattr(selected_stream, "ifnum", cfg.beam_override.ifnum or 0))
        cfg.beam_override.plnum = int(getattr(selected_stream, "plnum", cfg.beam_override.plnum or 0))
        cfg.beam_override.sampler = (str(getattr(selected_stream, "sampler", "")).strip() or cfg.beam_override.sampler)

    cfg.calibration.chopper_wheel = bool(args.chopper_wheel)
    cfg.calibration.tamb_k = args.tamb_k
    cfg.calibration.tamb_default_k = _choose_float_setting(args, argv, "--tamb-default-k", "tamb_default_k", global_cfg, "tamb_default_k", 300.0)
    cfg.calibration.tamb_min_k = _choose_float_setting(args, argv, "--tamb-min-k", "tamb_min_k", global_cfg, "tamb_min_k", 250.0)
    cfg.calibration.tamb_max_k = _choose_float_setting(args, argv, "--tamb-max-k", "tamb_max_k", global_cfg, "tamb_max_k", 330.0)
    legacy_weather_table = _choose_optional_str_setting(args, argv, "--weather-table", "weather_table", global_cfg, "weather_table", None)
    legacy_weather_time_col = _choose_optional_str_setting(args, argv, "--weather-time-col", "weather_time_col", global_cfg, "weather_time_col", None)
    cfg.calibration.weather_inside_table = _choose_optional_str_setting(args, argv, "--weather-inside-table", "weather_inside_table", global_cfg, "weather_inside_table", legacy_weather_table)
    cfg.calibration.weather_inside_table_suffix = _choose_optional_str_setting(args, argv, "--weather-inside-table-suffix", "weather_inside_table_suffix", global_cfg, "weather_inside_table_suffix", "weather-ambient") or "weather-ambient"
    cfg.calibration.weather_inside_time_col = _choose_optional_str_setting(args, argv, "--weather-inside-time-col", "weather_inside_time_col", global_cfg, "weather_inside_time_col", legacy_weather_time_col or "time") or (legacy_weather_time_col or "time")
    cfg.refraction.weather_outside_table = _choose_optional_str_setting(args, argv, "--weather-outside-table", "weather_outside_table", global_cfg, "weather_outside_table", legacy_weather_table)
    cfg.refraction.weather_outside_table_suffix = _choose_optional_str_setting(args, argv, "--weather-outside-table-suffix", "weather_outside_table_suffix", global_cfg, "weather_outside_table_suffix", "weather-ambient") or "weather-ambient"
    cfg.refraction.weather_outside_time_col = _choose_optional_str_setting(args, argv, "--weather-outside-time-col", "weather_outside_time_col", global_cfg, "weather_outside_time_col", legacy_weather_time_col or "time") or (legacy_weather_time_col or "time")
    legacy_weather_table_cli = _argv_has_option(argv, "--weather-table")
    legacy_weather_time_cli = _argv_has_option(argv, "--weather-time-col")
    if legacy_weather_table_cli:
        if not _argv_has_option(argv, "--weather-inside-table"):
            cfg.calibration.weather_inside_table = legacy_weather_table
        if not _argv_has_option(argv, "--weather-outside-table"):
            cfg.refraction.weather_outside_table = legacy_weather_table
    if legacy_weather_time_cli:
        if not _argv_has_option(argv, "--weather-inside-time-col"):
            cfg.calibration.weather_inside_time_col = legacy_weather_time_col or "time"
        if not _argv_has_option(argv, "--weather-outside-time-col"):
            cfg.refraction.weather_outside_time_col = legacy_weather_time_col or "time"
    cfg.refraction.outside_default_temperature_c = _choose_float_setting(args, argv, "--outside-default-temperature-c", "outside_default_temperature_c", global_cfg, "outside_default_temperature_c", 0.0)
    cfg.refraction.outside_default_pressure_hpa = _choose_float_setting(args, argv, "--outside-default-pressure-hpa", "outside_default_pressure_hpa", global_cfg, "outside_default_pressure_hpa", 760.0)
    cfg.refraction.outside_default_humidity_pct = _choose_float_setting(args, argv, "--outside-default-humidity-pct", "outside_default_humidity_pct", global_cfg, "outside_default_humidity_pct", 30.0)
    cfg.refraction.outside_temperature_min_c = _choose_float_setting(args, argv, "--outside-temperature-min-c", "outside_temperature_min_c", global_cfg, "outside_temperature_min_c", -50.0)
    cfg.refraction.outside_temperature_max_c = _choose_float_setting(args, argv, "--outside-temperature-max-c", "outside_temperature_max_c", global_cfg, "outside_temperature_max_c", 50.0)
    cfg.refraction.outside_pressure_min_hpa = _choose_float_setting(args, argv, "--outside-pressure-min-hpa", "outside_pressure_min_hpa", global_cfg, "outside_pressure_min_hpa", 400.0)
    cfg.refraction.outside_pressure_max_hpa = _choose_float_setting(args, argv, "--outside-pressure-max-hpa", "outside_pressure_max_hpa", global_cfg, "outside_pressure_max_hpa", 1100.0)
    cfg.refraction.outside_humidity_min_pct = _choose_float_setting(args, argv, "--outside-humidity-min-pct", "outside_humidity_min_pct", global_cfg, "outside_humidity_min_pct", 0.0)
    cfg.refraction.outside_humidity_max_pct = _choose_float_setting(args, argv, "--outside-humidity-max-pct", "outside_humidity_max_pct", global_cfg, "outside_humidity_max_pct", 100.0)
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
    if (selected_stream is not None) and (not _argv_has_option(argv, "--hpbw-init-arcsec")):
        cfg.edge_fit.hpbw_init_arcsec = _estimate_hpbw_init_arcsec(
            restfreq_hz=cfg.beam_override.restfreq_hz,
            dish_diameter_m=float(args.dish_diameter_m),
            hpbw_factor=float(args.hpbw_factor),
            fallback_arcsec=float(cfg.edge_fit.hpbw_init_arcsec),
        )
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
    if argv is None:
        argv = list(sys.argv[1:])
    ap = argparse.ArgumentParser(description="Run single-stream sun_scan-compatible analysis via the shared package API.")
    add_singlebeam_arguments(ap)
    args = ap.parse_args(argv)
    cfg = config_from_args(args, argv=argv)
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
