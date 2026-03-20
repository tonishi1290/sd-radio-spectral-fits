from __future__ import annotations

import argparse
import copy
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence
import json
import math
import traceback
import sys
import re

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
    return _canonicalize_boresight_settings(dict(cfg.get("global", {}) or {}))


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



def estimate_hpbw_init_arcsec(restfreq_hz: Optional[float], dish_diameter_m: float, hpbw_factor: float, fallback_arcsec: float) -> float:
    if restfreq_hz is None or not math.isfinite(restfreq_hz) or restfreq_hz <= 0:
        return float(fallback_arcsec)
    c_m_s = 299792458.0
    return float(hpbw_factor) * (c_m_s / float(restfreq_hz)) / float(dish_diameter_m) * 206265.0



def stream_output_tag(base_tag: str, stream_name: str) -> str:
    safe = str(stream_name).replace("/", "_").replace(" ", "_")
    return f"{base_tag}__{safe}"



def _sanitize_path_tag_component(value: Any) -> str:
    text = str(value).strip()
    if not text:
        return "unnamed"
    safe = re.sub(r"[^0-9A-Za-z._-]+", "_", text)
    safe = safe.strip("._-")
    return safe or "unnamed"



def _default_multi_base_tag(rawdata_path: Path) -> str:
    name = _sanitize_path_tag_component(rawdata_path.name)
    parent = _sanitize_path_tag_component(rawdata_path.parent.name)
    if parent and parent != "unnamed" and parent != name:
        return f"{parent}__{name}"
    return name



def _deduplicate_tags(tags: Sequence[str]) -> List[str]:
    counts: Dict[str, int] = {}
    unique: List[str] = []
    for tag in tags:
        base = _sanitize_path_tag_component(tag)
        n = counts.get(base, 0) + 1
        counts[base] = n
        unique.append(base if n == 1 else f"{base}__{n:02d}")
    return unique



def _prepare_multi_run_entries(rawdata_paths: Sequence[Path], run_ids: Optional[Sequence[Optional[str]]] = None) -> List[Dict[str, Any]]:
    normalized_paths = [Path(p).expanduser().resolve() for p in rawdata_paths]
    if run_ids is not None and len(run_ids) != len(normalized_paths):
        raise ValueError("run_ids must have the same length as rawdata_paths")
    name_counts: Dict[str, int] = {}
    for path in normalized_paths:
        key = str(path.name)
        name_counts[key] = name_counts.get(key, 0) + 1
    candidate_tags: List[str] = []
    for idx, path in enumerate(normalized_paths):
        explicit_run_id = None if run_ids is None else run_ids[idx]
        if explicit_run_id is not None and str(explicit_run_id).strip():
            candidate_tags.append(str(explicit_run_id))
        elif name_counts.get(str(path.name), 0) > 1:
            candidate_tags.append(_default_multi_base_tag(path))
        else:
            candidate_tags.append(_sanitize_path_tag_component(path.name))
    unique_tags = _deduplicate_tags(candidate_tags)
    entries: List[Dict[str, Any]] = []
    for path, tag in zip(normalized_paths, unique_tags):
        entries.append({
            "rawdata_path": path,
            "base_tag": tag,
            "run_id": tag,
        })
    return entries



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
    base_tag: Optional[str] = None,
) -> Dict[str, Any]:
    config = load_spectrometer_config(spectrometer_config)
    validation = validate_spectrometer_config(spectrometer_config, explicit_stream_names=stream_names)
    all_streams = list(config.get("streams", []) or [])
    extras_by_name = stream_extras_by_name(spectrometer_config)
    streams = filter_streams_for_purpose(spectrometer_config, config, "sunscan", explicit_stream_names=stream_names)
    if not streams:
        raise ValueError("no streams selected for extraction")

    base_tag = str(base_tag or rawdata_path.name)
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
        stream_cfg = _stream_override_dict(stream)
        restfreq_hz = restfreq_hz_for_stream(stream)
        dish_diameter_m = float(stream_cfg.get("dish_diameter_m", base_config.dish_diameter_m))
        hpbw_factor = float(stream_cfg.get("hpbw_factor", base_config.hpbw_factor))
        if stream_cfg.get("hpbw_init_arcsec") is not None:
            hpbw_init_arcsec = float(stream_cfg.get("hpbw_init_arcsec"))
        elif base_config.runtime.hpbw_init_explicit:
            hpbw_init_arcsec = float(base_config.edge_fit.hpbw_init_arcsec)
        else:
            hpbw_init_arcsec = estimate_hpbw_init_arcsec(
                restfreq_hz,
                dish_diameter_m=dish_diameter_m,
                hpbw_factor=hpbw_factor,
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
        cfg_stream.dish_diameter_m = dish_diameter_m
        cfg_stream.hpbw_factor = hpbw_factor
        stream_azel_source = str(stream_cfg.get("azel_source", stream_cfg.get("boresight_source", base_config.input.azel_source)))
        stream_corr_raw = None
        if stream_cfg.get("azel_correction_apply") is not None:
            stream_corr_raw = stream_cfg.get("azel_correction_apply")
        elif stream_cfg.get("altaz_apply") is not None:
            stream_corr_raw = stream_cfg.get("altaz_apply")
        if stream_corr_raw is not None:
            stream_azel_correction_apply = _normalize_legacy_azel_correction_apply(stream_corr_raw)
        elif base_config.runtime.azel_correction_apply_explicit:
            stream_azel_correction_apply = str(base_config.input.azel_correction_apply)
        else:
            stream_azel_correction_apply = "subtract" if stream_azel_source.strip().lower() == "encoder" else "none"
        cfg_stream.input.azel_source = stream_azel_source
        cfg_stream.input.azel_correction_apply = stream_azel_correction_apply
        cfg_stream.input.altaz_apply = stream_azel_correction_apply
        # Input/timing-related per-stream overrides.  These matter for real runs
        # because different spectrometers can have different timing offsets or
        # even different supporting NECST tables.
        input_override_specs = [
            ("encoder_table", "input", "encoder_table", str),
            ("encoder_table_suffix", "input", "encoder_table_suffix", str),
            ("altaz_table", "input", "altaz_table", str),
            ("altaz_table_suffix", "input", "altaz_table_suffix", str),
            ("encoder_time_col", "input", "encoder_time_col", str),
            ("altaz_time_col", "input", "altaz_time_col", str),
            ("spectrometer_time_offset_sec", "input", "spectrometer_time_offset_sec", float),
            ("encoder_shift_sec", "input", "encoder_shift_sec", float),
            ("encoder_az_time_offset_sec", "input", "encoder_az_time_offset_sec", float),
            ("encoder_el_time_offset_sec", "input", "encoder_el_time_offset_sec", float),
            ("encoder_vavg_sec", "input", "encoder_vavg_sec", float),
        ]
        override_specs = [
            ("chopper_wheel", "calibration", "chopper_wheel", _coerce_bool),
            ("tamb_k", "calibration", "tamb_k", float),
            ("tamb_default_k", "calibration", "tamb_default_k", float),
            ("tamb_min_k", "calibration", "tamb_min_k", float),
            ("tamb_max_k", "calibration", "tamb_max_k", float),
            ("weather_inside_table", "calibration", "weather_inside_table", str),
            ("weather_inside_table_suffix", "calibration", "weather_inside_table_suffix", str),
            ("weather_inside_time_col", "calibration", "weather_inside_time_col", str),
            ("weather_outside_table", "refraction", "weather_outside_table", str),
            ("weather_outside_table_suffix", "refraction", "weather_outside_table_suffix", str),
            ("weather_outside_time_col", "refraction", "weather_outside_time_col", str),
            ("outside_default_temperature_c", "refraction", "outside_default_temperature_c", float),
            ("outside_default_pressure_hpa", "refraction", "outside_default_pressure_hpa", float),
            ("outside_default_humidity_pct", "refraction", "outside_default_humidity_pct", float),
            ("outside_temperature_min_c", "refraction", "outside_temperature_min_c", float),
            ("outside_temperature_max_c", "refraction", "outside_temperature_max_c", float),
            ("outside_pressure_min_hpa", "refraction", "outside_pressure_min_hpa", float),
            ("outside_pressure_max_hpa", "refraction", "outside_pressure_max_hpa", float),
            ("outside_humidity_min_pct", "refraction", "outside_humidity_min_pct", float),
            ("outside_humidity_max_pct", "refraction", "outside_humidity_max_pct", float),
            ("chopper_win_sec", "calibration", "chopper_win_sec", float),
            ("chopper_stat", "calibration", "chopper_stat", str),
            ("profile_xlim_deg", "profile", "profile_xlim_deg", float),
            ("ripple_remove", "ripple", "enabled", _coerce_bool),
            ("ripple_preset", "ripple", "preset", str),
            ("ripple_model", "ripple", "model", str),
            ("ripple_target_hz", "ripple", "target_hz", float),
            ("ripple_search_hz", "ripple", "search_hz", float),
            ("ripple_bw_hz", "ripple", "bw_hz", float),
            ("ripple_max_harm", "ripple", "max_harm", int),
            ("ripple_order", "ripple", "order", int),
            ("ripple_notch_pass", "ripple", "notch_passes", int),
            ("ripple_trend_win_sec", "ripple", "trend_win_sec", float),
            ("ripple_resample_dt_sec", "ripple", "resample_dt_sec", float),
            ("ripple_eval_band_hz", "ripple", "eval_band_hz", float),
            ("edge_fit", "edge_fit", "enabled", _coerce_bool),
            ("edge_fit_win_deg", "edge_fit", "fit_win_deg", float),
            ("edge_fit_threshold", "edge_fit", "fit_threshold", float),
            ("strict_deriv", "edge_fit", "strict_deriv", _coerce_bool),
            ("edge_fit_plot_max_scans", "report", "edge_fit_plot_max_scans", int),
            ("trim_scan", "trim", "enabled", _coerce_bool),
            ("trim_vfrac", "trim", "vfrac", float),
            ("trim_vmin", "trim", "vmin", float),
            ("trim_gap", "trim", "gap_fill", int),
            ("trim_min_samples", "trim", "min_samples", int),
            ("trim_dominant_axis", "trim", "dominant_axis", _coerce_bool),
            ("trim_axis_ratio_min", "trim", "ratio_min", float),
            ("trim_vpercentile", "trim", "vpercentile", float),
            ("trim_steady_scan", "trim", "steady_scan", _coerce_bool),
            ("trim_use_on_only", "trim", "use_on_only", _coerce_bool),
            ("trim_xwin_factor", "trim", "xwin_factor", float),
            ("trim_cross_offset_max_deg", "trim", "cross_offset_max_deg", float),
            ("trim_scan_speed_min_arcsec", "trim", "speed_min_deg_s", lambda v: float(v) / 3600.0),
            ("trim_steady_cv_max", "trim", "steady_cv_max", float),
            ("continue_on_error", "runtime", "continue_on_error", _coerce_bool),
            ("debug_plot", "report", "debug_plot", _coerce_bool),
            ("encoder_vavg_sec", "input", "encoder_vavg_sec", float),
        ]
        # Legacy/symmetric aliases for per-stream weather settings.
        legacy_stream_weather_table = stream_cfg.get("weather_table")
        legacy_stream_weather_time_col = stream_cfg.get("weather_time_col")
        if legacy_stream_weather_table is not None:
            if stream_cfg.get("weather_inside_table") is None:
                cfg_stream.calibration.weather_inside_table = str(legacy_stream_weather_table)
            if stream_cfg.get("weather_outside_table") is None:
                cfg_stream.refraction.weather_outside_table = str(legacy_stream_weather_table)
        if legacy_stream_weather_time_col is not None:
            if stream_cfg.get("weather_inside_time_col") is None:
                cfg_stream.calibration.weather_inside_time_col = str(legacy_stream_weather_time_col)
            if stream_cfg.get("weather_outside_time_col") is None:
                cfg_stream.refraction.weather_outside_time_col = str(legacy_stream_weather_time_col)

        for key, section_name, attr_name, caster in input_override_specs + override_specs:
            if stream_cfg.get(key) is None:
                continue
            section_obj = getattr(cfg_stream, section_name)
            setattr(section_obj, attr_name, caster(stream_cfg.get(key)))
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



def run_extract_many(
    rawdata_paths: Sequence[Path],
    spectrometer_config: Path,
    base_config: SunScanAnalysisConfig,
    *,
    outdir: Path,
    run_ids: Optional[Sequence[Optional[str]]] = None,
    stream_names: Optional[Sequence[str]] = None,
    merged_tag: Optional[str] = None,
) -> Dict[str, Any]:
    entries = _prepare_multi_run_entries(rawdata_paths, run_ids=run_ids)
    if not entries:
        raise ValueError("rawdata_paths must not be empty")

    outdir = Path(outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    per_run_root = outdir / "per_run"
    merged_root = outdir / "merged"
    per_run_root.mkdir(parents=True, exist_ok=True)
    merged_root.mkdir(parents=True, exist_ok=True)

    merged_label = _sanitize_path_tag_component(merged_tag or "all")
    all_rows: List[pd.DataFrame] = []
    manifest_rows: List[pd.DataFrame] = []
    run_rows: List[Dict[str, Any]] = []
    per_run_outputs: List[Dict[str, Any]] = []

    for seq, entry in enumerate(entries, start=1):
        rawdata_path = Path(entry["rawdata_path"]).expanduser().resolve()
        base_tag = str(entry["base_tag"])
        run_id = str(entry["run_id"])
        per_run_outdir = per_run_root / base_tag
        per_run_outdir.mkdir(parents=True, exist_ok=True)
        cfg_run = copy.deepcopy(base_config)
        cfg_run.input.rawdata_path = rawdata_path
        cfg_run.report.outdir = per_run_outdir
        try:
            outputs = run_extract(
                rawdata_path=rawdata_path,
                spectrometer_config=Path(spectrometer_config).expanduser().resolve(),
                base_config=cfg_run,
                outdir=per_run_outdir,
                run_id=run_id,
                stream_names=stream_names,
                base_tag=base_tag,
            )
            all_df = outputs.get("all_df", pd.DataFrame()).copy()
            if not all_df.empty:
                if "rawdata_path" not in all_df.columns:
                    all_df.insert(0, "rawdata_path", str(rawdata_path))
                if "base_tag" not in all_df.columns:
                    all_df.insert(0, "base_tag", base_tag)
                if "run_id" not in all_df.columns:
                    all_df.insert(0, "run_id", run_id)
            manifest_df = outputs.get("manifest_df", pd.DataFrame()).copy()
            if manifest_df.empty:
                manifest_df = pd.DataFrame()
            if "rawdata_path" not in manifest_df.columns:
                manifest_df.insert(0, "rawdata_path", str(rawdata_path))
            if "base_tag" not in manifest_df.columns:
                manifest_df.insert(0, "base_tag", base_tag)
            if "run_id" not in manifest_df.columns:
                manifest_df.insert(0, "run_id", run_id)
            all_rows.append(all_df)
            manifest_rows.append(manifest_df)
            outputs["base_tag"] = base_tag
            outputs["run_id"] = run_id
            outputs["rawdata_path"] = rawdata_path
            outputs["per_run_outdir"] = per_run_outdir
            per_run_outputs.append(outputs)
            run_rows.append({
                "sequence_index": seq,
                "run_id": run_id,
                "base_tag": base_tag,
                "rawdata_path": str(rawdata_path),
                "per_run_outdir": str(per_run_outdir),
                "status": "ok",
                "error": None,
                "all_summary_csv": str(outputs.get("all_summary_csv")) if outputs.get("all_summary_csv") is not None else None,
                "manifest_csv": str(outputs.get("manifest_csv")) if outputs.get("manifest_csv") is not None else None,
                "stream_table_csv": str(outputs.get("stream_table_csv")) if outputs.get("stream_table_csv") is not None else None,
                "config_snapshot_json": str(outputs.get("config_snapshot_json")) if outputs.get("config_snapshot_json") is not None else None,
                "n_summary_rows": int(len(all_df)),
                "n_manifest_rows": int(len(manifest_df)),
            })
        except Exception as exc:
            run_rows.append({
                "sequence_index": seq,
                "run_id": run_id,
                "base_tag": base_tag,
                "rawdata_path": str(rawdata_path),
                "per_run_outdir": str(per_run_outdir),
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
                "all_summary_csv": None,
                "manifest_csv": None,
                "stream_table_csv": None,
                "config_snapshot_json": None,
                "n_summary_rows": 0,
                "n_manifest_rows": 0,
            })
            (per_run_outdir / "run_extract_many_error.txt").write_text(traceback.format_exc(), encoding="utf-8")
            if not base_config.runtime.continue_on_error:
                raise

    merged_all_df = pd.concat(all_rows, axis=0, ignore_index=True) if all_rows else pd.DataFrame()
    preferred_prefix = [
        "run_id", "base_tag", "rawdata_path", "tag", "stream_name", "beam_id", "scan_id",
        "x_arcsec", "y_arcsec",
        "source_db_path", "spectral_name", "db_stream_name", "restfreq_hz", "analysis_config_digest",
        "polariza", "fdnum", "ifnum", "plnum", "sampler",
    ]
    ordered_cols = [c for c in preferred_prefix if c in merged_all_df.columns] + [c for c in merged_all_df.columns if c not in preferred_prefix]
    if ordered_cols:
        merged_all_df = merged_all_df.loc[:, ordered_cols]
    merged_manifest_df = pd.concat(manifest_rows, axis=0, ignore_index=True) if manifest_rows else pd.DataFrame()
    run_table_df = pd.DataFrame(run_rows)

    all_csv = merged_root / f"sunscan_multibeam_scan_summary_{merged_label}.csv"
    manifest_csv = merged_root / f"sunscan_multibeam_manifest_{merged_label}.csv"
    run_table_csv = merged_root / f"sunscan_multibeam_run_table_{merged_label}.csv"
    merged_all_df.to_csv(all_csv, index=False)
    merged_manifest_df.to_csv(manifest_csv, index=False)
    run_table_df.to_csv(run_table_csv, index=False)

    return {
        "all_summary_csv": all_csv,
        "manifest_csv": manifest_csv,
        "run_table_csv": run_table_csv,
        "merged_outdir": merged_root,
        "per_run_root": per_run_root,
        "all_df": merged_all_df,
        "manifest_df": merged_manifest_df,
        "run_table_df": run_table_df,
        "per_run_outputs": per_run_outputs,
    }



def add_extract_arguments(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("rawdata", nargs="+", help="One or more RawData directories containing the NECST database")
    ap.add_argument("--spectrometer-config", required=True, help="converter-compatible spectrometer config TOML")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--run-id", default=None, help="Override run_id for a single RawData input. When multiple RawData directories are given, this value is used as the merged output tag.")
    ap.add_argument("--db-namespace", default="necst", help="Database namespace/prefix used in NECST DB table names")
    ap.add_argument("--telescope", default="OMU1P85M", help="Telescope name used in NECST DB table names")
    ap.add_argument("--tel-loaddata", default="OMU1p85m", help="Telescope name passed to nercst loaddb()")
    ap.add_argument("--planet", default="sun", help="Target body name passed to astropy.get_body()")
    ap.add_argument("--stream-name", dest="stream_names", action="append", default=None, help="Select a specific stream name (repeatable)")

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
    rawdata_value = args.rawdata[0] if isinstance(args.rawdata, (list, tuple)) else args.rawdata
    cfg = SunScanAnalysisConfig.default(rawdata_path=Path(rawdata_value).expanduser().resolve(), spectral_name="__unused__", outdir=Path(args.outdir).expanduser().resolve())
    global_cfg = _load_global_cfg(getattr(args, "spectrometer_config", None))

    cfg.input.db_namespace = _choose_str_setting(args, argv, "--db-namespace", "db_namespace", global_cfg, "db_namespace", "necst")
    cfg.input.telescope = _choose_str_setting(args, argv, "--telescope", "telescope", global_cfg, "telescope", "OMU1P85M")
    cfg.input.tel_loaddata = _choose_str_setting(args, argv, "--tel-loaddata", "tel_loaddata", global_cfg, "tel_loaddata", "OMU1p85m")
    cfg.input.planet = _choose_str_setting(args, argv, "--planet", "planet", global_cfg, "planet", "sun")
    cfg.runtime.azel_source_explicit = bool(_argv_has_option(argv, "--azel-source", "--boresight-source") or (global_cfg.get("azel_source") is not None) or (global_cfg.get("boresight_source") is not None))
    if (not getattr(args, "spectrometer_config", None)) or _argv_has_option(argv, "--azel-source", "--boresight-source"):
        cfg.input.azel_source = str(getattr(args, "azel_source", "encoder"))
    elif global_cfg.get("azel_source") is not None:
        cfg.input.azel_source = str(global_cfg.get("azel_source"))
    else:
        cfg.input.azel_source = str(getattr(args, "azel_source", "encoder"))
    cfg.runtime.azel_correction_apply_explicit = bool(
        _argv_has_option(argv, "--azel-correction-apply", "--boresight-correction-apply", "--altaz-apply")
        or (global_cfg.get("azel_correction_apply") is not None)
        or (global_cfg.get("boresight_correction_apply") is not None)
        or (global_cfg.get("altaz_apply") is not None)
    )
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
    cfg.calibration.chopper_wheel = _choose_bool_setting(args, argv, ("--chopper-wheel", "--no-chopper-wheel"), "chopper_wheel", global_cfg, "chopper_wheel", True)
    cfg.calibration.tamb_k = _choose_optional_float_setting(args, argv, "--tamb-k", "tamb_k", global_cfg, "tamb_k", None)
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
    cfg.calibration.chopper_win_sec = _choose_float_setting(args, argv, "--chopper-win-sec", "chopper_win_sec", global_cfg, "chopper_win_sec", 5.0)
    cfg.calibration.chopper_stat = _choose_str_setting(args, argv, "--chopper-stat", "chopper_stat", global_cfg, "chopper_stat", "mean")
    cfg.profile.profile_xlim_deg = _choose_float_setting(args, argv, "--profile-xlim-deg", "profile_xlim_deg", global_cfg, "profile_xlim_deg", 1.0)
    cfg.ripple.enabled = _choose_bool_setting(args, argv, ("--ripple-remove", "--ripple-no-remove"), "ripple_remove", global_cfg, "ripple_remove", True)
    cfg.ripple.preset = _choose_str_setting(args, argv, "--ripple-preset", "ripple_preset", global_cfg, "ripple_preset", "auto")
    cfg.ripple.model = _choose_str_setting(args, argv, "--ripple-model", "ripple_model", global_cfg, "ripple_model", "auto")
    cfg.ripple.target_hz = _choose_float_setting(args, argv, "--ripple-target-hz", "ripple_target_hz", global_cfg, "ripple_target_hz", 1.2)
    cfg.ripple.search_hz = _choose_float_setting(args, argv, "--ripple-search-hz", "ripple_search_hz", global_cfg, "ripple_search_hz", 0.3)
    cfg.ripple.bw_hz = _choose_optional_float_setting(args, argv, "--ripple-bw-hz", "ripple_bw_hz", global_cfg, "ripple_bw_hz", None)
    cfg.ripple.max_harm = _choose_optional_int_setting(args, argv, "--ripple-max-harm", "ripple_max_harm", global_cfg, "ripple_max_harm", None)
    cfg.ripple.order = _choose_optional_int_setting(args, argv, "--ripple-order", "ripple_order", global_cfg, "ripple_order", None)
    cfg.ripple.notch_passes = _choose_optional_int_setting(args, argv, "--ripple-notch-pass", "ripple_notch_pass", global_cfg, "ripple_notch_pass", None)
    cfg.ripple.trend_win_sec = _choose_optional_float_setting(args, argv, "--ripple-trend-win-sec", "ripple_trend_win_sec", global_cfg, "ripple_trend_win_sec", None)
    cfg.ripple.resample_dt_sec = _choose_optional_float_setting(args, argv, "--ripple-resample-dt-sec", "ripple_resample_dt_sec", global_cfg, "ripple_resample_dt_sec", None)
    cfg.ripple.eval_band_hz = _choose_optional_float_setting(args, argv, "--ripple-eval-band-hz", "ripple_eval_band_hz", global_cfg, "ripple_eval_band_hz", None)
    cfg.edge_fit.enabled = _choose_bool_setting(args, argv, ("--edge-fit", "--no-edge-fit"), "edge_fit", global_cfg, "edge_fit", True)
    cfg.edge_fit.fit_win_deg = _choose_float_setting(args, argv, "--edge-fit-win-deg", "edge_fit_win_deg", global_cfg, "edge_fit_win_deg", 0.15)
    cfg.edge_fit.fit_threshold = _choose_float_setting(args, argv, "--edge-fit-threshold", "edge_fit_threshold", global_cfg, "edge_fit_threshold", 0.20)
    cfg.edge_fit.hpbw_init_arcsec = _choose_float_setting(args, argv, "--hpbw-init-arcsec", "hpbw_init_arcsec", global_cfg, "hpbw_init_arcsec", 324.0)
    cfg.runtime.hpbw_init_explicit = bool(_argv_has_option(argv, "--hpbw-init-arcsec") or (global_cfg.get("hpbw_init_arcsec") is not None))
    cfg.edge_fit.strict_deriv = _choose_bool_setting(args, argv, ("--strict-deriv", "--no-strict-deriv"), "strict_deriv", global_cfg, "strict_deriv", True)
    cfg.report.edge_fit_plot_max_scans = _choose_int_setting(args, argv, "--edge-fit-plot-max-scans", "edge_fit_plot_max_scans", global_cfg, "edge_fit_plot_max_scans", 3)
    cfg.trim.enabled = _choose_bool_setting(args, argv, ("--trim-scan", "--no-trim-scan"), "trim_scan", global_cfg, "trim_scan", True)
    cfg.trim.vfrac = _choose_float_setting(args, argv, "--trim-vfrac", "trim_vfrac", global_cfg, "trim_vfrac", 0.20)
    cfg.trim.vmin = _choose_float_setting(args, argv, "--trim-vmin", "trim_vmin", global_cfg, "trim_vmin", 1e-4)
    cfg.trim.gap_fill = _choose_int_setting(args, argv, "--trim-gap", "trim_gap", global_cfg, "trim_gap", 10)
    cfg.trim.min_samples = _choose_int_setting(args, argv, "--trim-min-samples", "trim_min_samples", global_cfg, "trim_min_samples", 100)
    cfg.trim.dominant_axis = _choose_bool_setting(args, argv, ("--trim-dominant-axis", "--trim-no-dominant-axis"), "trim_dominant_axis", global_cfg, "trim_dominant_axis", True)
    cfg.trim.ratio_min = _choose_float_setting(args, argv, "--trim-axis-ratio-min", "trim_axis_ratio_min", global_cfg, "trim_axis_ratio_min", 3.0)
    cfg.trim.vpercentile = _choose_float_setting(args, argv, "--trim-vpercentile", "trim_vpercentile", global_cfg, "trim_vpercentile", 95.0)
    cfg.trim.steady_scan = _choose_bool_setting(args, argv, ("--trim-steady-scan", "--trim-no-steady-scan"), "trim_steady_scan", global_cfg, "trim_steady_scan", True)
    cfg.trim.use_on_only = _choose_bool_setting(args, argv, ("--trim-use-on-only", "--trim-include-hotoff"), "trim_use_on_only", global_cfg, "trim_use_on_only", True)
    cfg.trim.xwin_factor = _choose_float_setting(args, argv, "--trim-xwin-factor", "trim_xwin_factor", global_cfg, "trim_xwin_factor", 1.2)
    cfg.trim.cross_offset_max_deg = _choose_float_setting(args, argv, "--trim-cross-offset-max-deg", "trim_cross_offset_max_deg", global_cfg, "trim_cross_offset_max_deg", 0.5)
    cfg.trim.speed_min_deg_s = _choose_float_setting(args, argv, "--trim-scan-speed-min-arcsec", "trim_scan_speed_min_arcsec", global_cfg, "trim_scan_speed_min_arcsec", 20.0) / 3600.0
    cfg.trim.steady_cv_max = _choose_float_setting(args, argv, "--trim-steady-cv-max", "trim_steady_cv_max", global_cfg, "trim_steady_cv_max", 0.8)
    cfg.runtime.continue_on_error = _choose_bool_setting(args, argv, ("--continue-on-error", "--no-continue-on-error"), "continue_on_error", global_cfg, "continue_on_error", False)
    cfg.report.debug_plot = _choose_bool_setting(args, argv, ("--debug-plot", "--no-debug-plot"), "debug_plot", global_cfg, "debug_plot", False)
    cfg.dish_diameter_m = _choose_float_setting(args, argv, "--dish-diameter-m", "dish_diameter_m", global_cfg, "dish_diameter_m", 1.85)
    cfg.hpbw_factor = _choose_float_setting(args, argv, "--hpbw-factor", "hpbw_factor", global_cfg, "hpbw_factor", 1.2)
    return cfg


def main(argv: Optional[Sequence[str]] = None) -> None:
    if argv is None:
        argv = list(sys.argv[1:])
    ap = argparse.ArgumentParser(description="Run sun_scan-compatible analysis for every configured stream in one or more RawData directories.")
    add_extract_arguments(ap)
    args = ap.parse_args(argv)
    cfg = config_from_args(args, argv=argv)
    rawdata_paths = [Path(p).expanduser().resolve() for p in args.rawdata]
    if len(rawdata_paths) == 1:
        outputs = run_extract(
            rawdata_path=rawdata_paths[0],
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
        return
    outputs = run_extract_many(
        rawdata_paths=rawdata_paths,
        spectrometer_config=Path(args.spectrometer_config).expanduser().resolve(),
        base_config=cfg,
        outdir=Path(args.outdir).expanduser().resolve(),
        run_ids=None,
        stream_names=args.stream_names,
        merged_tag=args.run_id,
    )
    print(f"[done] wrote {outputs['all_summary_csv']}")
    print(f"[done] wrote {outputs['manifest_csv']}")
    print(f"[done] wrote {outputs['run_table_csv']}")


if __name__ == "__main__":
    main()
