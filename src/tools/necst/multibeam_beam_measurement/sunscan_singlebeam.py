from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, Optional, Sequence
import math
import sys

from .config_io import load_spectrometer_config, resolve_stream_usage_policy, restfreq_hz_for_stream, stream_extras_by_name
from .public_api import run_singlebeam, run_singlebeam_many
from .sunscan_config import SunScanAnalysisConfig
from .sunscan_report import build_summary_lines


def _argv_has_option(argv: Optional[Sequence[str]], *opts: str) -> bool:
    sargv = [str(x) for x in (argv or [])]
    for a in sargv:
        for opt in opts:
            if a == opt or a.startswith(opt + "="):
                return True
    return False


def _choose_float_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: float = 0.0, *, stream_cfg: Optional[Dict[str, Any]] = None, stream_key: Optional[str] = None) -> float:
    key = stream_key or global_key
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt):
        value = getattr(args, attr, default)
    elif stream_cfg is not None and stream_cfg.get(key) is not None:
        value = stream_cfg.get(key)
    elif global_cfg.get(global_key) is not None:
        value = global_cfg.get(global_key)
    else:
        value = getattr(args, attr, default)
    return float(value)


def _choose_str_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: str, *, stream_cfg: Optional[Dict[str, Any]] = None, stream_key: Optional[str] = None) -> str:
    key = stream_key or global_key
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt):
        value = getattr(args, attr, default)
    elif stream_cfg is not None and stream_cfg.get(key) is not None:
        value = stream_cfg.get(key)
    elif global_cfg.get(global_key) is not None:
        value = global_cfg.get(global_key)
    else:
        value = getattr(args, attr, default)
    return str(value)


def _choose_optional_str_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: Optional[str] = None, *, stream_cfg: Optional[Dict[str, Any]] = None, stream_key: Optional[str] = None) -> Optional[str]:
    key = stream_key or global_key
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt):
        value = getattr(args, attr, default)
    elif stream_cfg is not None and stream_cfg.get(key) is not None:
        value = stream_cfg.get(key)
    elif global_cfg.get(global_key) is not None:
        value = global_cfg.get(global_key)
    else:
        value = getattr(args, attr, default)
    if value is None:
        return default
    s = str(value).strip()
    return s if s else default


def _coerce_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise ValueError(f"cannot coerce to bool: {value!r}")


def _canonicalize_boresight_settings(d: Dict[str, Any]) -> Dict[str, Any]:
    if d.get("azel_source") is None and d.get("boresight_source") is not None:
        d["azel_source"] = d.get("boresight_source")
    if d.get("azel_correction_apply") is None and d.get("boresight_correction_apply") is not None:
        d["azel_correction_apply"] = d.get("boresight_correction_apply")
    return d


def _stream_override_dict(stream: Any) -> Dict[str, Any]:
    raw = getattr(stream, "override", None)
    override: Dict[str, Any] = {}
    if isinstance(raw, dict):
        override.update(raw)
    elif raw is not None:
        try:
            override.update(dict(raw))
        except Exception:
            pass
    return _canonicalize_boresight_settings(override)


def _choose_bool_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: Any, attr: str, global_cfg: Dict[str, Any], global_key: str, default: bool = False, *, stream_cfg: Optional[Dict[str, Any]] = None, stream_key: Optional[str] = None) -> bool:
    key = stream_key or global_key
    cli_opts = tuple(cli_opt) if isinstance(cli_opt, (list, tuple)) else (str(cli_opt),)
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, *cli_opts):
        value = getattr(args, attr, default)
    elif stream_cfg is not None and stream_cfg.get(key) is not None:
        value = stream_cfg.get(key)
    elif global_cfg.get(global_key) is not None:
        value = global_cfg.get(global_key)
    else:
        value = getattr(args, attr, default)
    return bool(_coerce_bool(value))


def _choose_int_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: int = 0, *, stream_cfg: Optional[Dict[str, Any]] = None, stream_key: Optional[str] = None) -> int:
    key = stream_key or global_key
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt):
        value = getattr(args, attr, default)
    elif stream_cfg is not None and stream_cfg.get(key) is not None:
        value = stream_cfg.get(key)
    elif global_cfg.get(global_key) is not None:
        value = global_cfg.get(global_key)
    else:
        value = getattr(args, attr, default)
    return int(value)


def _choose_optional_float_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: Optional[float] = None, *, stream_cfg: Optional[Dict[str, Any]] = None, stream_key: Optional[str] = None) -> Optional[float]:
    key = stream_key or global_key
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt):
        value = getattr(args, attr, default)
    elif stream_cfg is not None and stream_cfg.get(key) is not None:
        value = stream_cfg.get(key)
    elif global_cfg.get(global_key) is not None:
        value = global_cfg.get(global_key)
    else:
        value = getattr(args, attr, default)
    if value is None:
        return default
    return float(value)


def _choose_optional_int_setting(args: argparse.Namespace, argv: Optional[Sequence[str]], cli_opt: str, attr: str, global_cfg: Dict[str, Any], global_key: str, default: Optional[int] = None, *, stream_cfg: Optional[Dict[str, Any]] = None, stream_key: Optional[str] = None) -> Optional[int]:
    key = stream_key or global_key
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, cli_opt):
        value = getattr(args, attr, default)
    elif stream_cfg is not None and stream_cfg.get(key) is not None:
        value = stream_cfg.get(key)
    elif global_cfg.get(global_key) is not None:
        value = global_cfg.get(global_key)
    else:
        value = getattr(args, attr, default)
    if value is None:
        return default
    return int(value)


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


def _resolve_azel_correction_apply(args: argparse.Namespace, argv: Optional[Sequence[str]], global_cfg: Dict[str, Any], azel_source: str, *, stream_cfg: Optional[Dict[str, Any]] = None) -> str:
    raw = None
    if _argv_has_option(argv, "--azel-correction-apply", "--boresight-correction-apply"):
        raw = getattr(args, "azel_correction_apply", None)
    elif _argv_has_option(argv, "--altaz-apply"):
        raw = getattr(args, "altaz_apply", None)
    elif stream_cfg is not None and stream_cfg.get("azel_correction_apply") is not None:
        raw = stream_cfg.get("azel_correction_apply")
    elif stream_cfg is not None and stream_cfg.get("boresight_correction_apply") is not None:
        raw = stream_cfg.get("boresight_correction_apply")
    elif stream_cfg is not None and stream_cfg.get("altaz_apply") is not None:
        raw = stream_cfg.get("altaz_apply")
    elif global_cfg.get("azel_correction_apply") is not None:
        raw = global_cfg.get("azel_correction_apply")
    elif global_cfg.get("boresight_correction_apply") is not None:
        raw = global_cfg.get("boresight_correction_apply")
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
    ap.add_argument("rawdata", nargs="+", help="One or more RawData directories containing the NECST database")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--spectrometer-config", default=None, help="converter-compatible spectrometer config TOML")
    ap.add_argument("--run-id", default=None, help="Merged output tag for multi-run mode; also used as output tag in single-run mode when explicitly set")
    ap.add_argument("--stream-name", default=None, help="Select one stream from --spectrometer-config")
    ap.add_argument("--db-namespace", default="necst", help="Database namespace/prefix used in NECST DB table names")
    ap.add_argument("--telescope", default="OMU1P85M", help="Telescope name used in NECST DB table names")
    ap.add_argument("--tel-loaddata", default="OMU1p85m", help="Telescope name passed to nercst loaddb()")
    ap.add_argument("--planet", default="sun", help="Target body name passed to astropy.get_body()")
    ap.add_argument("--spectral-name", default="xffts-board1", help="Spectral stream name")

    ap.add_argument("--azel-source", "--boresight-source", dest="azel_source", choices=["encoder", "altaz"], default="encoder", help="Raw Az/El source used to construct boresight before sunscan fitting")
    ap.add_argument("--azel-correction-apply", "--boresight-correction-apply", dest="azel_correction_apply", choices=["none", "subtract", "add"], default=None, help="How to combine raw boresight Az/El with dlon/dlat when constructing corrected boresight")
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
    ap.add_argument("--chopper-wheel", dest="chopper_wheel", action="store_true")
    ap.add_argument("--no-chopper-wheel", dest="chopper_wheel", action="store_false")
    ap.add_argument("--tamb-k", type=float, default=None)
    ap.add_argument("--chopper-win-sec", type=float, default=5.0)
    ap.add_argument("--chopper-stat", choices=["median", "mean"], default="mean")
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
    ap.add_argument("--ripple-remove", dest="ripple_remove", action="store_true")
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

    ap.add_argument("--edge-fit", dest="edge_fit", action="store_true")
    ap.add_argument("--no-edge-fit", dest="edge_fit", action="store_false")
    ap.add_argument("--edge-fit-win-deg", type=float, default=0.15)
    ap.add_argument("--edge-fit-threshold", type=float, default=0.20)
    ap.add_argument("--hpbw-init-arcsec", type=float, default=324.0)
    ap.add_argument("--edge-fit-plot-max-scans", type=int, default=3)
    ap.set_defaults(edge_fit=True)

    ap.add_argument("--trim-scan", dest="trim_scan", action="store_true")
    ap.add_argument("--no-trim-scan", dest="trim_scan", action="store_false")
    ap.add_argument("--trim-vfrac", type=float, default=0.20)
    ap.add_argument("--trim-vmin", type=float, default=1e-4)
    ap.add_argument("--trim-gap", type=int, default=10)
    ap.add_argument("--trim-min-samples", type=int, default=100)
    ap.add_argument("--trim-dominant-axis", dest="trim_dominant_axis", action="store_true")
    ap.add_argument("--trim-no-dominant-axis", dest="trim_dominant_axis", action="store_false")
    ap.add_argument("--trim-axis-ratio-min", type=float, default=3.0)
    ap.add_argument("--trim-vpercentile", type=float, default=95.0)
    ap.add_argument("--trim-steady-scan", dest="trim_steady_scan", action="store_true")
    ap.add_argument("--trim-no-steady-scan", dest="trim_steady_scan", action="store_false")
    ap.add_argument("--trim-use-on-only", dest="trim_use_on_only", action="store_true")
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

    ap.add_argument("--continue-on-error", dest="continue_on_error", action="store_true", default=False)
    ap.add_argument("--no-continue-on-error", dest="continue_on_error", action="store_false")
    ap.add_argument("--debug-plot", dest="debug_plot", action="store_true", default=False)
    ap.add_argument("--no-debug-plot", dest="debug_plot", action="store_false")
    ap.add_argument("--dish-diameter-m", type=float, default=1.85)
    ap.add_argument("--hpbw-factor", type=float, default=1.2)



def config_from_args(args: argparse.Namespace, argv: Optional[Sequence[str]] = None) -> SunScanAnalysisConfig:
    resolved_spectral_name = getattr(args, "spectral_name", None)

    global_cfg: Dict[str, Any] = {}
    selected_stream = None
    stream_cfg: Dict[str, Any] = {}
    if args.spectrometer_config:
        config_path = Path(args.spectrometer_config).expanduser().resolve()
        config_dict = load_spectrometer_config(config_path)
        global_cfg = _canonicalize_boresight_settings(dict(config_dict.get("global", {}) or {}))
        if (not _argv_has_option(argv, "--spectral-name")) and (global_cfg.get("spectral_name") is not None):
            gs = str(global_cfg.get("spectral_name")).strip()
            if gs:
                resolved_spectral_name = gs
        selected_stream = _select_stream_from_config(config_dict, stream_name=getattr(args, "stream_name", None), spectral_name=resolved_spectral_name)
        stream_cfg = _stream_override_dict(selected_stream) if selected_stream is not None else {}
        explicit_stream_override = bool(getattr(args, "stream_name", None)) or _argv_has_option(argv, "--spectral-name") or (global_cfg.get("spectral_name") is not None)
        if (selected_stream is not None) and (not explicit_stream_override):
            usage_policy = resolve_stream_usage_policy(stream_extras_by_name(config_path).get(str(getattr(selected_stream, "name", "")), {}))
            if not bool(usage_policy.get("use_for_sunscan", True)):
                raise ValueError(
                    f"selected stream {getattr(selected_stream, 'name', None)!r} is disabled for sunscan by enabled/use_for_sunscan in the spectrometer config. "
                    "Pass --stream-name explicitly if you intentionally want to override this for a one-off run."
                )

    rawdata_values = list(args.rawdata) if isinstance(args.rawdata, (list, tuple)) else [args.rawdata]
    first_rawdata = Path(rawdata_values[0]).expanduser().resolve()
    cfg = SunScanAnalysisConfig.default(rawdata_path=first_rawdata, spectral_name=(resolved_spectral_name or args.spectral_name), outdir=Path(args.outdir).expanduser().resolve())

    cfg.input.db_namespace = _choose_str_setting(args, argv, "--db-namespace", "db_namespace", global_cfg, "db_namespace", "necst")
    cfg.input.telescope = _choose_str_setting(args, argv, "--telescope", "telescope", global_cfg, "telescope", "OMU1P85M")
    cfg.input.tel_loaddata = _choose_str_setting(args, argv, "--tel-loaddata", "tel_loaddata", global_cfg, "tel_loaddata", "OMU1p85m")
    cfg.input.planet = _choose_str_setting(args, argv, "--planet", "planet", global_cfg, "planet", "sun")
    cfg.runtime.azel_source_explicit = bool(
        _argv_has_option(argv, "--azel-source", "--boresight-source")
        or (stream_cfg.get("azel_source") is not None)
        or (stream_cfg.get("boresight_source") is not None)
        or (global_cfg.get("azel_source") is not None)
        or (global_cfg.get("boresight_source") is not None)
    )
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, "--azel-source", "--boresight-source"):
        cfg.input.azel_source = str(getattr(args, "azel_source", "encoder"))
    elif stream_cfg.get("azel_source") is not None:
        cfg.input.azel_source = str(stream_cfg.get("azel_source"))
    elif global_cfg.get("azel_source") is not None:
        cfg.input.azel_source = str(global_cfg.get("azel_source"))
    else:
        cfg.input.azel_source = str(getattr(args, "azel_source", "encoder"))
    cfg.runtime.azel_correction_apply_explicit = bool(
        _argv_has_option(argv, "--azel-correction-apply", "--boresight-correction-apply", "--altaz-apply")
        or (stream_cfg.get("azel_correction_apply") is not None)
        or (stream_cfg.get("boresight_correction_apply") is not None)
        or (stream_cfg.get("altaz_apply") is not None)
        or (global_cfg.get("azel_correction_apply") is not None)
        or (global_cfg.get("boresight_correction_apply") is not None)
        or (global_cfg.get("altaz_apply") is not None)
    )
    cfg.input.azel_correction_apply = _resolve_azel_correction_apply(args, argv, global_cfg, cfg.input.azel_source, stream_cfg=stream_cfg)
    cfg.input.altaz_apply = cfg.input.azel_correction_apply
    cfg.input.encoder_table = _choose_optional_str_setting(args, argv, "--encoder-table", "encoder_table", global_cfg, "encoder_table", None, stream_cfg=stream_cfg)
    cfg.input.encoder_table_suffix = _choose_optional_str_setting(args, argv, "--encoder-table-suffix", "encoder_table_suffix", global_cfg, "encoder_table_suffix", "ctrl-antenna-encoder", stream_cfg=stream_cfg) or "ctrl-antenna-encoder"
    cfg.input.altaz_table = _choose_optional_str_setting(args, argv, "--altaz-table", "altaz_table", global_cfg, "altaz_table", None, stream_cfg=stream_cfg)
    cfg.input.altaz_table_suffix = _choose_optional_str_setting(args, argv, "--altaz-table-suffix", "altaz_table_suffix", global_cfg, "altaz_table_suffix", "ctrl-antenna-altaz", stream_cfg=stream_cfg) or "ctrl-antenna-altaz"
    cfg.input.encoder_time_col = _choose_str_setting(args, argv, "--encoder-time-col", "encoder_time_col", global_cfg, "encoder_time_col", "time", stream_cfg=stream_cfg)
    cfg.input.altaz_time_col = _choose_str_setting(args, argv, "--altaz-time-col", "altaz_time_col", global_cfg, "altaz_time_col", "time", stream_cfg=stream_cfg)
    cfg.input.spectrometer_time_offset_sec = _choose_float_setting(args, argv, "--spectrometer-time-offset-sec", "spectrometer_time_offset_sec", global_cfg, "spectrometer_time_offset_sec", 0.0, stream_cfg=stream_cfg)
    cfg.input.encoder_shift_sec = _choose_float_setting(args, argv, "--encoder-shift-sec", "encoder_shift_sec", global_cfg, "encoder_shift_sec", 0.0, stream_cfg=stream_cfg)
    cfg.input.encoder_az_time_offset_sec = _choose_float_setting(args, argv, "--encoder-az-time-offset-sec", "encoder_az_time_offset_sec", global_cfg, "encoder_az_time_offset_sec", 0.0, stream_cfg=stream_cfg)
    cfg.input.encoder_el_time_offset_sec = _choose_float_setting(args, argv, "--encoder-el-time-offset-sec", "encoder_el_time_offset_sec", global_cfg, "encoder_el_time_offset_sec", 0.0, stream_cfg=stream_cfg)
    cfg.input.encoder_vavg_sec = _choose_float_setting(args, argv, "--encoder-vavg-sec", "encoder_vavg_sec", global_cfg, "encoder_vavg_sec", 0.0, stream_cfg=stream_cfg)
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

    cfg.calibration.chopper_wheel = _choose_bool_setting(args, argv, ("--chopper-wheel", "--no-chopper-wheel"), "chopper_wheel", global_cfg, "chopper_wheel", True, stream_cfg=stream_cfg)
    cfg.calibration.tamb_k = _choose_optional_float_setting(args, argv, "--tamb-k", "tamb_k", global_cfg, "tamb_k", None, stream_cfg=stream_cfg)
    cfg.calibration.tamb_default_k = _choose_float_setting(args, argv, "--tamb-default-k", "tamb_default_k", global_cfg, "tamb_default_k", 300.0, stream_cfg=stream_cfg)
    cfg.calibration.tamb_min_k = _choose_float_setting(args, argv, "--tamb-min-k", "tamb_min_k", global_cfg, "tamb_min_k", 250.0, stream_cfg=stream_cfg)
    cfg.calibration.tamb_max_k = _choose_float_setting(args, argv, "--tamb-max-k", "tamb_max_k", global_cfg, "tamb_max_k", 330.0, stream_cfg=stream_cfg)
    legacy_weather_table = _choose_optional_str_setting(args, argv, "--weather-table", "weather_table", global_cfg, "weather_table", None, stream_cfg=stream_cfg)
    legacy_weather_time_col = _choose_optional_str_setting(args, argv, "--weather-time-col", "weather_time_col", global_cfg, "weather_time_col", None, stream_cfg=stream_cfg)
    cfg.calibration.weather_inside_table = _choose_optional_str_setting(args, argv, "--weather-inside-table", "weather_inside_table", global_cfg, "weather_inside_table", legacy_weather_table, stream_cfg=stream_cfg)
    cfg.calibration.weather_inside_table_suffix = _choose_optional_str_setting(args, argv, "--weather-inside-table-suffix", "weather_inside_table_suffix", global_cfg, "weather_inside_table_suffix", "weather-ambient", stream_cfg=stream_cfg) or "weather-ambient"
    cfg.calibration.weather_inside_time_col = _choose_optional_str_setting(args, argv, "--weather-inside-time-col", "weather_inside_time_col", global_cfg, "weather_inside_time_col", legacy_weather_time_col or "time", stream_cfg=stream_cfg) or (legacy_weather_time_col or "time")
    cfg.refraction.weather_outside_table = _choose_optional_str_setting(args, argv, "--weather-outside-table", "weather_outside_table", global_cfg, "weather_outside_table", legacy_weather_table, stream_cfg=stream_cfg)
    cfg.refraction.weather_outside_table_suffix = _choose_optional_str_setting(args, argv, "--weather-outside-table-suffix", "weather_outside_table_suffix", global_cfg, "weather_outside_table_suffix", "weather-ambient", stream_cfg=stream_cfg) or "weather-ambient"
    cfg.refraction.weather_outside_time_col = _choose_optional_str_setting(args, argv, "--weather-outside-time-col", "weather_outside_time_col", global_cfg, "weather_outside_time_col", legacy_weather_time_col or "time", stream_cfg=stream_cfg) or (legacy_weather_time_col or "time")
    # Legacy weather aliases should work symmetrically for the selected stream
    # override as well as for CLI/global settings.  This matters when a stream
    # wants to specify a single weather table/time column and let it fan out to
    # both inside/outside uses without repeating both specific keys.
    legacy_stream_weather_table = stream_cfg.get("weather_table") if stream_cfg else None
    legacy_stream_weather_time_col = stream_cfg.get("weather_time_col") if stream_cfg else None
    if legacy_stream_weather_table is not None:
        if stream_cfg.get("weather_inside_table") is None:
            cfg.calibration.weather_inside_table = str(legacy_stream_weather_table)
        if stream_cfg.get("weather_outside_table") is None:
            cfg.refraction.weather_outside_table = str(legacy_stream_weather_table)
    if legacy_stream_weather_time_col is not None:
        if stream_cfg.get("weather_inside_time_col") is None:
            cfg.calibration.weather_inside_time_col = str(legacy_stream_weather_time_col)
        if stream_cfg.get("weather_outside_time_col") is None:
            cfg.refraction.weather_outside_time_col = str(legacy_stream_weather_time_col)

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
    cfg.refraction.outside_default_temperature_c = _choose_float_setting(args, argv, "--outside-default-temperature-c", "outside_default_temperature_c", global_cfg, "outside_default_temperature_c", 0.0, stream_cfg=stream_cfg)
    cfg.refraction.outside_default_pressure_hpa = _choose_float_setting(args, argv, "--outside-default-pressure-hpa", "outside_default_pressure_hpa", global_cfg, "outside_default_pressure_hpa", 760.0, stream_cfg=stream_cfg)
    cfg.refraction.outside_default_humidity_pct = _choose_float_setting(args, argv, "--outside-default-humidity-pct", "outside_default_humidity_pct", global_cfg, "outside_default_humidity_pct", 30.0, stream_cfg=stream_cfg)
    cfg.refraction.outside_temperature_min_c = _choose_float_setting(args, argv, "--outside-temperature-min-c", "outside_temperature_min_c", global_cfg, "outside_temperature_min_c", -50.0, stream_cfg=stream_cfg)
    cfg.refraction.outside_temperature_max_c = _choose_float_setting(args, argv, "--outside-temperature-max-c", "outside_temperature_max_c", global_cfg, "outside_temperature_max_c", 50.0, stream_cfg=stream_cfg)
    cfg.refraction.outside_pressure_min_hpa = _choose_float_setting(args, argv, "--outside-pressure-min-hpa", "outside_pressure_min_hpa", global_cfg, "outside_pressure_min_hpa", 400.0, stream_cfg=stream_cfg)
    cfg.refraction.outside_pressure_max_hpa = _choose_float_setting(args, argv, "--outside-pressure-max-hpa", "outside_pressure_max_hpa", global_cfg, "outside_pressure_max_hpa", 1100.0, stream_cfg=stream_cfg)
    cfg.refraction.outside_humidity_min_pct = _choose_float_setting(args, argv, "--outside-humidity-min-pct", "outside_humidity_min_pct", global_cfg, "outside_humidity_min_pct", 0.0, stream_cfg=stream_cfg)
    cfg.refraction.outside_humidity_max_pct = _choose_float_setting(args, argv, "--outside-humidity-max-pct", "outside_humidity_max_pct", global_cfg, "outside_humidity_max_pct", 100.0, stream_cfg=stream_cfg)
    cfg.calibration.chopper_win_sec = _choose_float_setting(args, argv, "--chopper-win-sec", "chopper_win_sec", global_cfg, "chopper_win_sec", 5.0, stream_cfg=stream_cfg)
    cfg.calibration.chopper_stat = _choose_str_setting(args, argv, "--chopper-stat", "chopper_stat", global_cfg, "chopper_stat", "mean", stream_cfg=stream_cfg)
    cfg.profile.profile_xlim_deg = _choose_float_setting(args, argv, "--profile-xlim-deg", "profile_xlim_deg", global_cfg, "profile_xlim_deg", 1.0, stream_cfg=stream_cfg)
    cfg.ripple.enabled = _choose_bool_setting(args, argv, ("--ripple-remove", "--ripple-no-remove"), "ripple_remove", global_cfg, "ripple_remove", True, stream_cfg=stream_cfg)
    cfg.ripple.preset = _choose_str_setting(args, argv, "--ripple-preset", "ripple_preset", global_cfg, "ripple_preset", "auto", stream_cfg=stream_cfg)
    cfg.ripple.model = _choose_str_setting(args, argv, "--ripple-model", "ripple_model", global_cfg, "ripple_model", "auto", stream_cfg=stream_cfg)
    cfg.ripple.target_hz = _choose_float_setting(args, argv, "--ripple-target-hz", "ripple_target_hz", global_cfg, "ripple_target_hz", 1.2, stream_cfg=stream_cfg)
    cfg.ripple.search_hz = _choose_float_setting(args, argv, "--ripple-search-hz", "ripple_search_hz", global_cfg, "ripple_search_hz", 0.3, stream_cfg=stream_cfg)
    cfg.ripple.bw_hz = _choose_optional_float_setting(args, argv, "--ripple-bw-hz", "ripple_bw_hz", global_cfg, "ripple_bw_hz", None, stream_cfg=stream_cfg)
    cfg.ripple.max_harm = _choose_optional_int_setting(args, argv, "--ripple-max-harm", "ripple_max_harm", global_cfg, "ripple_max_harm", None, stream_cfg=stream_cfg)
    cfg.ripple.order = _choose_optional_int_setting(args, argv, "--ripple-order", "ripple_order", global_cfg, "ripple_order", None, stream_cfg=stream_cfg)
    cfg.ripple.notch_passes = _choose_optional_int_setting(args, argv, "--ripple-notch-pass", "ripple_notch_pass", global_cfg, "ripple_notch_pass", None, stream_cfg=stream_cfg)
    cfg.ripple.trend_win_sec = _choose_optional_float_setting(args, argv, "--ripple-trend-win-sec", "ripple_trend_win_sec", global_cfg, "ripple_trend_win_sec", None, stream_cfg=stream_cfg)
    cfg.ripple.resample_dt_sec = _choose_optional_float_setting(args, argv, "--ripple-resample-dt-sec", "ripple_resample_dt_sec", global_cfg, "ripple_resample_dt_sec", None, stream_cfg=stream_cfg)
    cfg.ripple.eval_band_hz = _choose_optional_float_setting(args, argv, "--ripple-eval-band-hz", "ripple_eval_band_hz", global_cfg, "ripple_eval_band_hz", None, stream_cfg=stream_cfg)
    cfg.edge_fit.enabled = _choose_bool_setting(args, argv, ("--edge-fit", "--no-edge-fit"), "edge_fit", global_cfg, "edge_fit", True, stream_cfg=stream_cfg)
    cfg.edge_fit.fit_win_deg = _choose_float_setting(args, argv, "--edge-fit-win-deg", "edge_fit_win_deg", global_cfg, "edge_fit_win_deg", 0.15, stream_cfg=stream_cfg)
    cfg.edge_fit.fit_threshold = _choose_float_setting(args, argv, "--edge-fit-threshold", "edge_fit_threshold", global_cfg, "edge_fit_threshold", 0.20, stream_cfg=stream_cfg)
    cfg.edge_fit.hpbw_init_arcsec = _choose_float_setting(args, argv, "--hpbw-init-arcsec", "hpbw_init_arcsec", global_cfg, "hpbw_init_arcsec", 324.0, stream_cfg=stream_cfg)
    cfg.runtime.hpbw_init_explicit = bool(_argv_has_option(argv, "--hpbw-init-arcsec") or (global_cfg.get("hpbw_init_arcsec") is not None) or (stream_cfg.get("hpbw_init_arcsec") is not None))
    cfg.edge_fit.strict_deriv = _choose_bool_setting(args, argv, ("--strict-deriv", "--no-strict-deriv"), "strict_deriv", global_cfg, "strict_deriv", True, stream_cfg=stream_cfg)
    cfg.report.edge_fit_plot_max_scans = _choose_int_setting(args, argv, "--edge-fit-plot-max-scans", "edge_fit_plot_max_scans", global_cfg, "edge_fit_plot_max_scans", 3, stream_cfg=stream_cfg)
    cfg.trim.enabled = _choose_bool_setting(args, argv, ("--trim-scan", "--no-trim-scan"), "trim_scan", global_cfg, "trim_scan", True, stream_cfg=stream_cfg)
    cfg.trim.vfrac = _choose_float_setting(args, argv, "--trim-vfrac", "trim_vfrac", global_cfg, "trim_vfrac", 0.20, stream_cfg=stream_cfg)
    cfg.trim.vmin = _choose_float_setting(args, argv, "--trim-vmin", "trim_vmin", global_cfg, "trim_vmin", 1e-4, stream_cfg=stream_cfg)
    cfg.trim.gap_fill = _choose_int_setting(args, argv, "--trim-gap", "trim_gap", global_cfg, "trim_gap", 10, stream_cfg=stream_cfg)
    cfg.trim.min_samples = _choose_int_setting(args, argv, "--trim-min-samples", "trim_min_samples", global_cfg, "trim_min_samples", 100, stream_cfg=stream_cfg)
    cfg.trim.dominant_axis = _choose_bool_setting(args, argv, ("--trim-dominant-axis", "--trim-no-dominant-axis"), "trim_dominant_axis", global_cfg, "trim_dominant_axis", True, stream_cfg=stream_cfg)
    cfg.trim.ratio_min = _choose_float_setting(args, argv, "--trim-axis-ratio-min", "trim_axis_ratio_min", global_cfg, "trim_axis_ratio_min", 3.0, stream_cfg=stream_cfg)
    cfg.trim.vpercentile = _choose_float_setting(args, argv, "--trim-vpercentile", "trim_vpercentile", global_cfg, "trim_vpercentile", 95.0, stream_cfg=stream_cfg)
    cfg.trim.steady_scan = _choose_bool_setting(args, argv, ("--trim-steady-scan", "--trim-no-steady-scan"), "trim_steady_scan", global_cfg, "trim_steady_scan", True, stream_cfg=stream_cfg)
    cfg.trim.use_on_only = _choose_bool_setting(args, argv, ("--trim-use-on-only", "--trim-include-hotoff"), "trim_use_on_only", global_cfg, "trim_use_on_only", True, stream_cfg=stream_cfg)
    cfg.trim.xwin_factor = _choose_float_setting(args, argv, "--trim-xwin-factor", "trim_xwin_factor", global_cfg, "trim_xwin_factor", 1.2, stream_cfg=stream_cfg)
    cfg.trim.cross_offset_max_deg = _choose_float_setting(args, argv, "--trim-cross-offset-max-deg", "trim_cross_offset_max_deg", global_cfg, "trim_cross_offset_max_deg", 0.5, stream_cfg=stream_cfg)
    cfg.trim.speed_min_deg_s = _choose_float_setting(args, argv, "--trim-scan-speed-min-arcsec", "trim_scan_speed_min_arcsec", global_cfg, "trim_scan_speed_min_arcsec", 20.0, stream_cfg=stream_cfg) / 3600.0
    cfg.trim.steady_cv_max = _choose_float_setting(args, argv, "--trim-steady-cv-max", "trim_steady_cv_max", global_cfg, "trim_steady_cv_max", 0.8, stream_cfg=stream_cfg)
    cfg.runtime.continue_on_error = _choose_bool_setting(args, argv, ("--continue-on-error", "--no-continue-on-error"), "continue_on_error", global_cfg, "continue_on_error", False, stream_cfg=stream_cfg)
    cfg.report.debug_plot = _choose_bool_setting(args, argv, ("--debug-plot", "--no-debug-plot"), "debug_plot", global_cfg, "debug_plot", False, stream_cfg=stream_cfg)
    cfg.dish_diameter_m = _choose_float_setting(args, argv, "--dish-diameter-m", "dish_diameter_m", global_cfg, "dish_diameter_m", 1.85, stream_cfg=stream_cfg)
    cfg.hpbw_factor = _choose_float_setting(args, argv, "--hpbw-factor", "hpbw_factor", global_cfg, "hpbw_factor", 1.2, stream_cfg=stream_cfg)

    if getattr(args, "run_id", None) and len(rawdata_values) == 1:
        cfg.report.tag = str(args.run_id).strip() or None

    if selected_stream is not None and not cfg.runtime.hpbw_init_explicit:
        cfg.edge_fit.hpbw_init_arcsec = _estimate_hpbw_init_arcsec(
            restfreq_hz=cfg.beam_override.restfreq_hz,
            dish_diameter_m=float(cfg.dish_diameter_m),
            hpbw_factor=float(cfg.hpbw_factor),
            fallback_arcsec=float(cfg.edge_fit.hpbw_init_arcsec),
        )
    return cfg


def main(argv: Optional[Sequence[str]] = None) -> None:
    if argv is None:
        argv = list(sys.argv[1:])
    ap = argparse.ArgumentParser(description="Run single-stream sun_scan-compatible analysis via the shared package API.")
    add_singlebeam_arguments(ap)
    args = ap.parse_args(argv)
    cfg = config_from_args(args, argv=argv)
    rawdata_paths = [Path(p).expanduser().resolve() for p in args.rawdata]

    if len(rawdata_paths) == 1:
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
        return

    outputs = run_singlebeam_many(
        rawdata_paths=rawdata_paths,
        base_config=cfg,
        outdir=Path(args.outdir).expanduser().resolve(),
        run_id=getattr(args, "run_id", None),
    )
    run_table_df = outputs["run_table_df"]
    success_count = int((run_table_df["status"] == "success").sum()) if (not run_table_df.empty and "status" in run_table_df.columns) else 0
    error_count = int((run_table_df["status"] == "error").sum()) if (not run_table_df.empty and "status" in run_table_df.columns) else 0
    if outputs.get("merged_summary_csv") is not None:
        print(f"[fit] wrote {Path(outputs['merged_summary_csv']).name}")
    print(f"[done] wrote {Path(outputs['run_table_csv']).name}")
    print(f"  outputs in: {Path(args.outdir).expanduser().resolve()}")
    print(f"  runs: total={len(rawdata_paths)} success={success_count} error={error_count}")


if __name__ == "__main__":
    main()
