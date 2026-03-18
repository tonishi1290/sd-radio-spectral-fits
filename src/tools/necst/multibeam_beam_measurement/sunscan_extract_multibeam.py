from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence
import json
import math
import traceback
import sys

import pandas as pd

from .config_io import (
    filter_streams_for_purpose,
    load_spectrometer_config,
    resolve_stream_usage_policy,
    restfreq_hz_for_stream,
    stream_extras_by_name,
    stream_table_from_config,
    validate_spectrometer_config,
)
from .sunscan_config import SunScanAnalysisConfig
from .sunscan_core import analyze_single_stream
from .sunscan_report import results_to_summary_dataframe, write_singlebeam_outputs


def _argv_has_option(argv: Optional[Sequence[str]], *opts: str) -> bool:
    sargv = [str(x) for x in (argv or [])]
    for a in sargv:
        for opt in opts:
            if a == opt or a.startswith(opt + "="):
                return True
    return False


def _load_global_cfg(spectrometer_config: Optional[str]) -> Dict[str, Any]:
    if not spectrometer_config:
        return {}
    cfg = load_spectrometer_config(Path(spectrometer_config).expanduser().resolve())
    return dict(cfg.get("global", {}) or {})


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



def estimate_hpbw_init_arcsec(restfreq_hz: Optional[float], dish_diameter_m: float, hpbw_factor: float, fallback_arcsec: float) -> float:
    if restfreq_hz is None or not math.isfinite(restfreq_hz) or restfreq_hz <= 0:
        return float(fallback_arcsec)
    c_m_s = 299792458.0
    return float(hpbw_factor) * (c_m_s / float(restfreq_hz)) / float(dish_diameter_m) * 206265.0



def stream_output_tag(base_tag: str, stream_name: str) -> str:
    safe = str(stream_name).replace("/", "_").replace(" ", "_")
    return f"{base_tag}__{safe}"



def _collect_singlebeam_rows(result, run_id: str, source_db_path: Path, stream: Any) -> pd.DataFrame:
    df = results_to_summary_dataframe(result.scan_ids, result.results)
    if df.empty:
        return df
    if "scan_id" in df.columns:
        df.insert(0, "scan_id", df.pop("scan_id"))
    df.insert(0, "beam_id", str(stream.beam.beam_id))
    df.insert(0, "stream_name", str(stream.name))
    df.insert(0, "tag", result.config.resolved_tag())
    df.insert(0, "run_id", str(run_id))
    if "center_az_deg" in df.columns:
        df["x_arcsec"] = pd.to_numeric(df["center_az_deg"], errors="coerce").astype(float) * 3600.0
    else:
        df["x_arcsec"] = float("nan")
    if "center_el_deg" in df.columns:
        df["y_arcsec"] = pd.to_numeric(df["center_el_deg"], errors="coerce").astype(float) * 3600.0
    else:
        df["y_arcsec"] = float("nan")
    df["source_db_path"] = str(source_db_path)
    df["spectral_name"] = str(result.config.input.spectral_name)
    df["restfreq_hz"] = restfreq_hz_for_stream(stream)
    df["analysis_config_digest"] = result.config.config_digest()
    df["polariza"] = str(stream.polariza)
    df["fdnum"] = int(stream.fdnum)
    df["ifnum"] = int(stream.ifnum)
    df["plnum"] = int(stream.plnum)
    df["sampler"] = stream.sampler
    if getattr(stream, "db_stream_name", None) is not None:
        df["db_stream_name"] = str(stream.db_stream_name)
    return df



def _empty_manifest_row(stream: Any, spectral_name: str, error: str, *, status: str = "error", skip_reason: Optional[str] = None, usage_policy: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    usage_policy = dict(usage_policy or {})
    return {
        "stream_name": str(stream.name),
        "beam_id": str(stream.beam.beam_id),
        "spectral_name": spectral_name,
        "summary_csv": None,
        "derivative_png_count": 0,
        "debug_png_count": 0,
        "spec_time_basis": None,
        "spec_time_suffix": None,
        "spec_time_fallback_field": None,
        "spec_time_example": None,
        "status": status,
        "skip_reason": skip_reason,
        "enabled": usage_policy.get("enabled"),
        "use_for_convert": usage_policy.get("use_for_convert"),
        "use_for_sunscan": usage_policy.get("use_for_sunscan"),
        "use_for_fit": usage_policy.get("use_for_fit"),
        "beam_fit_use": usage_policy.get("beam_fit_use"),
        "error": error,
    }



def _write_config_snapshot(path: Path, *, base_config: SunScanAnalysisConfig, validation, stream_usage: Dict[str, Dict[str, Any]]) -> Path:
    payload = {
        "rawdata_path": str(base_config.input.rawdata_path),
        "resolved_outdir": str(base_config.report.outdir),
        "config_digest": base_config.config_digest(),
        "analysis": {
            "db_namespace": base_config.input.db_namespace,
            "telescope": base_config.input.telescope,
            "tel_loaddata": base_config.input.tel_loaddata,
            "planet": base_config.input.planet,
            "azel_source": base_config.input.azel_source,
            "azel_correction_apply": base_config.input.azel_correction_apply,
            "altaz_apply": base_config.input.altaz_apply,
            "encoder_table": base_config.input.encoder_table,
            "encoder_table_suffix": base_config.input.encoder_table_suffix,
            "altaz_table": base_config.input.altaz_table,
            "altaz_table_suffix": base_config.input.altaz_table_suffix,
            "weather_inside_table": base_config.calibration.weather_inside_table,
            "weather_inside_table_suffix": base_config.calibration.weather_inside_table_suffix,
            "weather_inside_time_col": base_config.calibration.weather_inside_time_col,
            "weather_outside_table": base_config.refraction.weather_outside_table,
            "weather_outside_table_suffix": base_config.refraction.weather_outside_table_suffix,
            "weather_outside_time_col": base_config.refraction.weather_outside_time_col,
            "spectrometer_time_offset_sec": base_config.input.spectrometer_time_offset_sec,
            "encoder_shift_sec": base_config.input.encoder_shift_sec,
            "encoder_az_time_offset_sec": base_config.input.encoder_az_time_offset_sec,
            "encoder_el_time_offset_sec": base_config.input.encoder_el_time_offset_sec,
            "encoder_vavg_sec": base_config.input.encoder_vavg_sec,
            "chopper_wheel": base_config.calibration.chopper_wheel,
            "tamb_default_k": base_config.calibration.tamb_default_k,
            "tamb_min_k": base_config.calibration.tamb_min_k,
            "tamb_max_k": base_config.calibration.tamb_max_k,
            "outside_default_temperature_c": base_config.refraction.outside_default_temperature_c,
            "outside_default_pressure_hpa": base_config.refraction.outside_default_pressure_hpa,
            "outside_default_humidity_pct": base_config.refraction.outside_default_humidity_pct,
            "outside_temperature_min_c": base_config.refraction.outside_temperature_min_c,
            "outside_temperature_max_c": base_config.refraction.outside_temperature_max_c,
            "outside_pressure_min_hpa": base_config.refraction.outside_pressure_min_hpa,
            "outside_pressure_max_hpa": base_config.refraction.outside_pressure_max_hpa,
            "outside_humidity_min_pct": base_config.refraction.outside_humidity_min_pct,
            "outside_humidity_max_pct": base_config.refraction.outside_humidity_max_pct,
            "ripple_enabled": base_config.ripple.enabled,
            "trim_enabled": base_config.trim.enabled,
            "edge_fit_enabled": base_config.edge_fit.enabled,
            "hpbw_init_arcsec": base_config.edge_fit.hpbw_init_arcsec,
            "dish_diameter_m": base_config.dish_diameter_m,
            "hpbw_factor": base_config.hpbw_factor,
        },
        "validation": {
            "duplicate_beam_ids": validation.duplicate_beam_ids,
            "primary_streams": validation.primary_streams,
            "warnings": validation.warnings,
        },
        "stream_usage": stream_usage,
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path



def run_extract(
    rawdata_path: Path,
    spectrometer_config: Path,
    base_config: SunScanAnalysisConfig,
    *,
    outdir: Path,
    run_id: Optional[str] = None,
    stream_names: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    config = load_spectrometer_config(spectrometer_config)
    validation = validate_spectrometer_config(spectrometer_config, explicit_stream_names=stream_names)
    all_streams = list(config.get("streams", []) or [])
    extras_by_name = stream_extras_by_name(spectrometer_config)
    streams = filter_streams_for_purpose(spectrometer_config, config, "sunscan", explicit_stream_names=stream_names)
    if not streams:
        raise ValueError("no streams selected for extraction")

    base_tag = rawdata_path.name
    run_id = str(run_id or rawdata_path.name)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_rows: List[pd.DataFrame] = []
    manifests: List[Dict[str, Any]] = []

    stream_table = stream_table_from_config(spectrometer_config, config)
    if stream_names:
        stream_table = stream_table.loc[stream_table["stream_name"].astype(str).isin([str(s) for s in stream_names])].copy()
    stream_table_csv = outdir / f"spectrometer_stream_table_{base_tag}.csv"
    stream_table.to_csv(stream_table_csv, index=False)
    config_snapshot_json = outdir / f"analysis_config_snapshot_{base_tag}.json"
    _write_config_snapshot(config_snapshot_json, base_config=base_config, validation=validation, stream_usage=extras_by_name)

    selected_names = {str(stream.name) for stream in streams}
    if not stream_names:
        for stream in all_streams:
            stream_name = str(stream.name)
            if stream_name in selected_names:
                continue
            usage_policy = resolve_stream_usage_policy(extras_by_name.get(stream_name, {}))
            manifests.append(_empty_manifest_row(
                stream,
                str(getattr(stream, "db_stream_name", None) or stream.name),
                "",
                status="skipped",
                skip_reason="disabled by enabled/use_for_sunscan flags",
                usage_policy=usage_policy,
            ))

    for stream in streams:
        per_stream_outdir = outdir / "per_stream" / str(stream.name)
        per_stream_outdir.mkdir(parents=True, exist_ok=True)
        restfreq_hz = restfreq_hz_for_stream(stream)
        hpbw_init_arcsec = estimate_hpbw_init_arcsec(
            restfreq_hz,
            dish_diameter_m=float(base_config.dish_diameter_m),
            hpbw_factor=float(base_config.hpbw_factor),
            fallback_arcsec=float(base_config.edge_fit.hpbw_init_arcsec),
        )
        spectral_name = str(getattr(stream, "db_stream_name", None) or stream.name)
        cfg_stream = base_config.with_stream_override(
            spectral_name=spectral_name,
            stream_name=str(stream.name),
            beam_id=str(stream.beam.beam_id),
            restfreq_hz=restfreq_hz,
            hpbw_init_arcsec=hpbw_init_arcsec,
            polariza=str(stream.polariza),
            fdnum=int(stream.fdnum),
            ifnum=int(stream.ifnum),
            plnum=int(stream.plnum),
            sampler=str(stream.sampler) if stream.sampler is not None else None,
            outdir=per_stream_outdir,
            tag=stream_output_tag(base_tag, str(stream.name)),
        )
        try:
            result = analyze_single_stream(cfg_stream)
            result = write_singlebeam_outputs(result)
            rows = _collect_singlebeam_rows(result, run_id=run_id, source_db_path=rawdata_path, stream=stream)
            all_rows.append(rows)
            usage_policy = resolve_stream_usage_policy(extras_by_name.get(str(stream.name), {}))
            manifests.append({
                "stream_name": str(stream.name),
                "beam_id": str(stream.beam.beam_id),
                "spectral_name": str(cfg_stream.input.spectral_name),
                "restfreq_hz": restfreq_hz,
                "hpbw_init_arcsec": hpbw_init_arcsec,
                "summary_csv": str(result.report_paths.summary_csv) if result.report_paths.summary_csv else None,
                "derivative_png_count": int(len(result.report_paths.derivative_pngs)),
                "debug_png_count": int(len(result.report_paths.debug_pngs)),
                "spec_time_basis": result.df.attrs.get("spec_time_basis"),
                "spec_time_suffix": result.df.attrs.get("spec_time_suffix"),
                "spec_time_fallback_field": result.df.attrs.get("spec_time_fallback_field"),
                "spec_time_example": result.df.attrs.get("spec_time_example"),
                "status": "ok",
                "skip_reason": None,
                "enabled": usage_policy.get("enabled"),
                "use_for_convert": usage_policy.get("use_for_convert"),
                "use_for_sunscan": usage_policy.get("use_for_sunscan"),
                "use_for_fit": usage_policy.get("use_for_fit"),
                "beam_fit_use": usage_policy.get("beam_fit_use"),
                "error": None,
            })
        except Exception as exc:
            manifests.append(_empty_manifest_row(
                stream,
                spectral_name,
                f"{type(exc).__name__}: {exc}",
                usage_policy=resolve_stream_usage_policy(extras_by_name.get(str(stream.name), {})),
            ))
            (per_stream_outdir / "analysis_error.txt").write_text(traceback.format_exc(), encoding="utf-8")
            if not base_config.runtime.continue_on_error:
                raise

    all_df = pd.concat(all_rows, axis=0, ignore_index=True) if all_rows else pd.DataFrame()
    preferred_prefix = [
        "run_id", "tag", "stream_name", "beam_id", "scan_id",
        "x_arcsec", "y_arcsec",
        "source_db_path", "spectral_name", "db_stream_name", "restfreq_hz", "analysis_config_digest",
        "polariza", "fdnum", "ifnum", "plnum", "sampler",
    ]
    ordered_cols = [c for c in preferred_prefix if c in all_df.columns] + [c for c in all_df.columns if c not in preferred_prefix]
    if ordered_cols:
        all_df = all_df.loc[:, ordered_cols]
    all_csv = outdir / f"sunscan_multibeam_scan_summary_{base_tag}.csv"
    all_df.to_csv(all_csv, index=False)
    manifest_df = pd.DataFrame(manifests)
    manifest_csv = outdir / f"sunscan_multibeam_manifest_{base_tag}.csv"
    manifest_df.to_csv(manifest_csv, index=False)
    return {
        "all_summary_csv": all_csv,
        "manifest_csv": manifest_csv,
        "stream_table_csv": stream_table_csv,
        "config_snapshot_json": config_snapshot_json,
        "all_df": all_df,
        "manifest_df": manifest_df,
    }



def add_extract_arguments(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("rawdata", help="RawData directory containing the NECST database")
    ap.add_argument("--spectrometer-config", required=True, help="converter-compatible spectrometer config TOML")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--run-id", default=None, help="Override run_id (default: rawdata directory name)")
    ap.add_argument("--db-namespace", default="necst", help="Database namespace/prefix used in NECST DB table names")
    ap.add_argument("--telescope", default="OMU1P85M", help="Telescope name used in NECST DB table names")
    ap.add_argument("--tel-loaddata", default="OMU1p85m", help="Telescope name passed to nercst loaddb()")
    ap.add_argument("--planet", default="sun", help="Target body name passed to astropy.get_body()")
    ap.add_argument("--stream-name", dest="stream_names", action="append", default=None, help="Select a specific stream name (repeatable)")

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
    cfg = SunScanAnalysisConfig.default(rawdata_path=Path(args.rawdata).expanduser().resolve(), spectral_name="__unused__", outdir=Path(args.outdir).expanduser().resolve())
    global_cfg = _load_global_cfg(getattr(args, "spectrometer_config", None))

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
    ap = argparse.ArgumentParser(description="Run sun_scan-compatible analysis for every configured stream in a RawData directory.")
    add_extract_arguments(ap)
    args = ap.parse_args(argv)
    cfg = config_from_args(args, argv=argv)
    outputs = run_extract(
        rawdata_path=Path(args.rawdata).expanduser().resolve(),
        spectrometer_config=Path(args.spectrometer_config).expanduser().resolve(),
        base_config=cfg,
        outdir=Path(args.outdir).expanduser().resolve(),
        run_id=args.run_id,
        stream_names=args.stream_names,
    )
    print(f"[done] wrote {outputs['all_summary_csv']}")
    print(f"[done] wrote {outputs['manifest_csv']}")
    print(f"[done] wrote {outputs['stream_table_csv']}")
    print(f"[done] wrote {outputs['config_snapshot_json']}")


if __name__ == "__main__":
    main()
