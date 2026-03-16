from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence
import math
import tomllib

import pandas as pd

from .legacy_loaders import load_converter_module


ALLOWED_POLARIZA = {"RR", "LL", "RL", "LR", "XX", "YY", "XY", "YX"}


@dataclass
class LightBeam:
    beam_id: str
    az_offset_arcsec: float = 0.0
    el_offset_arcsec: float = 0.0
    rotation_mode: str = "none"
    reference_angle_deg: float = 0.0
    rotation_sign: float = 1.0
    dewar_angle_deg: float = 0.0


@dataclass
class LightWCS:
    restfreq_hz: float = float("nan")


@dataclass
class LightStream:
    name: str
    fdnum: int
    ifnum: int
    plnum: int
    polariza: str
    beam: LightBeam
    db_stream_name: str
    db_table_name: Optional[str] = None
    sampler: Optional[str] = None
    frontend: Optional[str] = None
    backend: Optional[str] = None
    frequency_axis: Dict[str, Any] = field(default_factory=dict)
    local_oscillators: Dict[str, Any] = field(default_factory=dict)
    override: Dict[str, Any] = field(default_factory=dict)
    enabled: bool = True
    use_for_convert: bool = True
    use_for_sunscan: bool = True
    use_for_fit: bool = True
    beam_fit_use: Optional[bool] = None
    wcs: Optional[LightWCS] = None
    wcs_full: Optional[LightWCS] = None


@dataclass
class SpectrometerValidationResult:
    stream_table: pd.DataFrame
    duplicate_beam_ids: List[str]
    primary_streams: List[str]
    warnings: List[str]



def load_raw_toml(config_path: Path) -> Dict[str, Any]:
    with open(config_path, "rb") as fh:
        return tomllib.load(fh)


def _coerce_optional_bool(value: Any) -> Optional[bool]:
    if value is None:
        return None
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


def resolve_stream_usage_policy(extra: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    extra = dict(extra or {})
    enabled_raw = _coerce_optional_bool(extra.get("enabled", None))
    enabled = True if enabled_raw is None else bool(enabled_raw)

    def _resolve_simple(key: str) -> bool:
        raw = _coerce_optional_bool(extra.get(key, None))
        if raw is None:
            return bool(enabled)
        return bool(enabled) and bool(raw)

    beam_fit_use = _coerce_optional_bool(extra.get("beam_fit_use", None))
    use_for_fit_raw = _coerce_optional_bool(extra.get("use_for_fit", None))
    if use_for_fit_raw is None:
        if beam_fit_use is None:
            use_for_fit = bool(enabled)
        else:
            use_for_fit = bool(enabled) and bool(beam_fit_use)
    else:
        use_for_fit = bool(enabled) and bool(use_for_fit_raw)

    return {
        "enabled": bool(enabled),
        "use_for_convert": _resolve_simple("use_for_convert"),
        "use_for_sunscan": _resolve_simple("use_for_sunscan"),
        "use_for_fit": bool(use_for_fit),
        "beam_fit_use": beam_fit_use,
    }


def stream_enabled_for_purpose(extra: Optional[Dict[str, Any]], purpose: str) -> bool:
    resolved = resolve_stream_usage_policy(extra)
    purpose_key = {
        "convert": "use_for_convert",
        "sunscan": "use_for_sunscan",
        "fit": "use_for_fit",
    }.get(str(purpose))
    if purpose_key is None:
        raise ValueError(f"unsupported purpose={purpose!r}")
    return bool(resolved[purpose_key])


def filter_streams_for_purpose(raw_config_path: Path, stream_cfg: Dict[str, Any], purpose: str, *, explicit_stream_names: Optional[Sequence[str]] = None) -> List[Any]:
    streams = list(stream_cfg.get("streams", []) or [])
    if explicit_stream_names:
        wanted = {str(s) for s in explicit_stream_names}
        selected = [stream for stream in streams if str(getattr(stream, "name", "")) in wanted]
        missing = sorted(wanted - {str(getattr(stream, "name", "")) for stream in selected})
        if missing:
            raise ValueError(f"stream(s) not found in spectrometer config: {missing}")
        return selected
    extras = stream_extras_by_name(raw_config_path)
    return [stream for stream in streams if stream_enabled_for_purpose(extras.get(str(getattr(stream, "name", "")), {}), purpose)]



def _fallback_load_spectrometer_config(config_path: Path) -> Dict[str, Any]:
    raw = load_raw_toml(config_path)
    blocks = list(raw.get("spectrometers", []) or [])
    streams: List[LightStream] = []
    for block in blocks:
        name = str(block.get("name", "")).strip()
        if not name:
            raise ValueError("Each [[spectrometers]] block must have non-empty name")
        polariza = str(block.get("polariza", "")).strip().upper()
        if polariza not in ALLOWED_POLARIZA:
            raise ValueError(f"stream {name!r} has unsupported polariza={polariza!r}")
        fdnum = int(block.get("fdnum", 0))
        ifnum = int(block.get("ifnum", 0))
        plnum = int(block.get("plnum", 0))
        beam_id = str(block.get("beam_id", f"B{fdnum:02d}")).strip() or f"B{fdnum:02d}"
        beam_block = dict(block.get("beam", {}) or {})
        beam = LightBeam(
            beam_id=beam_id,
            az_offset_arcsec=float(beam_block.get("az_offset_arcsec", 0.0)),
            el_offset_arcsec=float(beam_block.get("el_offset_arcsec", 0.0)),
            rotation_mode=str(beam_block.get("rotation_mode", "none")),
            reference_angle_deg=float(beam_block.get("reference_angle_deg", beam_block.get("reference_el_deg", 0.0))),
            rotation_sign=float(beam_block.get("rotation_sign", 1.0)),
            dewar_angle_deg=float(beam_block.get("dewar_angle_deg", 0.0)),
        )
        freq_axis = dict(block.get("frequency_axis", {}) or {})
        restfreq_hz = freq_axis.get("restfreq_hz", float("nan"))
        try:
            restfreq_hz = float(restfreq_hz)
        except Exception:
            restfreq_hz = float("nan")
        wcs = LightWCS(restfreq_hz=restfreq_hz)
        usage_policy = resolve_stream_usage_policy(block)
        streams.append(
            LightStream(
                name=name,
                fdnum=fdnum,
                ifnum=ifnum,
                plnum=plnum,
                polariza=polariza,
                beam=beam,
                db_stream_name=str(block.get("db_stream_name", name)).strip() or name,
                db_table_name=(str(block.get("db_table_name")).strip() if block.get("db_table_name") is not None else None),
                sampler=(str(block.get("sampler")).strip() if block.get("sampler") is not None else None),
                frontend=(str(block.get("frontend")).strip() if block.get("frontend") is not None else None),
                backend=(str(block.get("backend")).strip() if block.get("backend") is not None else None),
                frequency_axis=freq_axis,
                local_oscillators=dict(block.get("local_oscillators", {}) or {}),
                override=dict(block.get("override", {}) or {}),
                wcs=wcs,
                wcs_full=wcs,
                enabled=usage_policy.get("enabled", True),
                use_for_convert=usage_policy.get("use_for_convert", True),
                use_for_sunscan=usage_policy.get("use_for_sunscan", True),
                use_for_fit=usage_policy.get("use_for_fit", True),
                beam_fit_use=usage_policy.get("beam_fit_use", None),
            )
        )
    return {
        "schema_version": int(raw.get("schema_version", raw.get("config_version", 1))),
        "global": dict(raw.get("global", {}) or {}),
        "provenance": dict(raw.get("provenance", {}) or {}),
        "streams": streams,
    }



def load_spectrometer_config(config_path: Path):
    try:
        converter = load_converter_module()
        return converter.load_spectrometer_config(config_path)
    except Exception:
        return _fallback_load_spectrometer_config(Path(config_path))



def stream_extras_by_name(config_path: Path) -> Dict[str, Dict[str, Any]]:
    raw = load_raw_toml(config_path)
    extras: Dict[str, Dict[str, Any]] = {}
    for block in list(raw.get("spectrometers", []) or []):
        name = str(block.get("name", "")).strip()
        if not name:
            continue
        extras[name] = {
            "enabled": block.get("enabled", None),
            "use_for_convert": block.get("use_for_convert", None),
            "use_for_sunscan": block.get("use_for_sunscan", None),
            "use_for_fit": block.get("use_for_fit", None),
            "beam_fit_use": block.get("beam_fit_use", None),
        }
    return extras



def restfreq_hz_for_stream(stream: Any) -> Optional[float]:
    for obj in [getattr(stream, "wcs", None), getattr(stream, "wcs_full", None)]:
        if obj is not None:
            try:
                val = float(getattr(obj, "restfreq_hz"))
                if math.isfinite(val) and val > 0:
                    return val
            except Exception:
                pass
    try:
        val = float((getattr(stream, "frequency_axis", {}) or {}).get("restfreq_hz"))
        if math.isfinite(val) and val > 0:
            return val
    except Exception:
        pass
    return None



def stream_table_from_config(config_path: Path, stream_cfg: Optional[Dict[str, Any]] = None) -> pd.DataFrame:
    cfg = stream_cfg if stream_cfg is not None else load_spectrometer_config(config_path)
    extras = stream_extras_by_name(config_path)
    rows: List[Dict[str, Any]] = []
    for stream in list(cfg.get("streams", []) or []):
        usage_policy = resolve_stream_usage_policy(extras.get(str(stream.name), {}))
        rows.append(
            {
                "stream_name": str(stream.name),
                "beam_id": str(stream.beam.beam_id),
                "db_stream_name": str(getattr(stream, "db_stream_name", None) or stream.name),
                "restfreq_hz": restfreq_hz_for_stream(stream),
                "fdnum": int(stream.fdnum),
                "ifnum": int(stream.ifnum),
                "plnum": int(stream.plnum),
                "polariza": str(stream.polariza),
                "sampler": stream.sampler,
                "frontend": getattr(stream, "frontend", None),
                "backend": getattr(stream, "backend", None),
                "enabled": usage_policy.get("enabled", True),
                "use_for_convert": usage_policy.get("use_for_convert", True),
                "use_for_sunscan": usage_policy.get("use_for_sunscan", True),
                "use_for_fit": usage_policy.get("use_for_fit", True),
                "beam_fit_use": usage_policy.get("beam_fit_use", None),
                "rotation_mode": str(getattr(stream.beam, "rotation_mode", "none")),
                "reference_angle_deg": float(getattr(stream.beam, "reference_angle_deg", 0.0)),
                "rotation_sign": float(getattr(stream.beam, "rotation_sign", 1.0)),
                "dewar_angle_deg": float(getattr(stream.beam, "dewar_angle_deg", 0.0)),
                "az_offset_arcsec": float(getattr(stream.beam, "az_offset_arcsec", 0.0)),
                "el_offset_arcsec": float(getattr(stream.beam, "el_offset_arcsec", 0.0)),
            }
        )
    return pd.DataFrame(rows)



def resolve_primary_stream_names(raw_config_path: Path, stream_cfg: Dict[str, Any], explicit_stream_names: Optional[Sequence[str]] = None) -> List[str]:
    streams = filter_streams_for_purpose(raw_config_path, stream_cfg, "fit", explicit_stream_names=explicit_stream_names)
    if explicit_stream_names:
        return [str(s) for s in explicit_stream_names]
    extras = stream_extras_by_name(raw_config_path)
    primary: List[str] = []
    by_beam: Dict[str, List[Any]] = {}
    for stream in streams:
        by_beam.setdefault(str(stream.beam.beam_id), []).append(stream)
    for beam_id, group in by_beam.items():
        chosen = []
        for stream in group:
            policy = resolve_stream_usage_policy(extras.get(str(stream.name), {}))
            if policy.get("beam_fit_use") is True:
                chosen.append(stream.name)
        if len(chosen) > 1:
            raise ValueError(f"multiple primary streams declared for beam_id={beam_id!r}: {chosen}")
        if len(group) > 1 and not chosen:
            names = [s.name for s in group]
            raise ValueError(
                f"beam_id={beam_id!r} has multiple fit-enabled streams {names} but no primary stream was specified. "
                "Add use_for_fit=true to the desired streams and beam_fit_use=true to exactly one stream, or pass --fit-stream-name."
            )
        primary.extend(chosen if chosen else [group[0].name])
    return primary



def validate_spectrometer_config(config_path: Path, *, explicit_stream_names: Optional[Sequence[str]] = None) -> SpectrometerValidationResult:
    cfg = load_spectrometer_config(config_path)
    table = stream_table_from_config(config_path, cfg)
    warnings: List[str] = []
    if table.empty:
        raise ValueError("spectrometer config contains no streams")
    duplicate_beam_ids = sorted([str(b) for b, n in table["beam_id"].value_counts().items() if int(n) > 1])
    try:
        primary_streams = resolve_primary_stream_names(config_path, cfg, explicit_stream_names=explicit_stream_names)
    except Exception as exc:
        warnings.append(str(exc))
        primary_streams = list(explicit_stream_names or [])
    for _, row in table.iterrows():
        if str(row["polariza"]).upper() not in ALLOWED_POLARIZA:
            warnings.append(f"unsupported polariza for stream {row['stream_name']!r}: {row['polariza']!r}")
        restfreq = row.get("restfreq_hz")
        if restfreq is None or (isinstance(restfreq, float) and not math.isfinite(restfreq)):
            warnings.append(f"stream {row['stream_name']!r} has no finite restfreq_hz")
        rot_sign = row.get("rotation_sign")
        try:
            rot_sign_f = float(rot_sign)
            if rot_sign_f not in (-1.0, 0.0, 1.0):
                warnings.append(
                    f"stream {row['stream_name']!r} has rotation_sign={rot_sign_f!r}; "
                    "converter-compatible configurations should use -1, 0, or +1"
                )
        except Exception:
            warnings.append(f"stream {row['stream_name']!r} has non-numeric rotation_sign={rot_sign!r}")

    if primary_streams:
        primary_table = table.loc[table["stream_name"].astype(str).isin([str(s) for s in primary_streams])].copy()
        if not primary_table.empty:
            geom_mag = (
                pd.to_numeric(primary_table["az_offset_arcsec"], errors="coerce").fillna(0.0).abs()
                + pd.to_numeric(primary_table["el_offset_arcsec"], errors="coerce").fillna(0.0).abs()
            )
            if bool((geom_mag <= 0.0).all()):
                warnings.append(
                    "all selected primary streams have zero nominal beam offsets; "
                    "pseudo multi-beam dry-runs will not constrain rotation unless non-zero az/el offsets are set in [spectrometers.beam]"
                )
            rot_modes = primary_table["rotation_mode"].astype(str).str.lower().str.strip()
            if bool((rot_modes == "none").all()):
                warnings.append(
                    "all selected primary streams have rotation_mode='none'; pseudo multi-beam dry-runs will not exercise EL-dependent beam rotation"
                )
    return SpectrometerValidationResult(
        stream_table=table,
        duplicate_beam_ids=duplicate_beam_ids,
        primary_streams=primary_streams,
        warnings=warnings,
    )



def format_toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        raise ValueError("None is not representable in TOML")
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            raise ValueError("NaN is not representable in TOML")
        if math.isinf(value):
            raise ValueError("inf is not representable in TOML")
        return repr(float(value))
    if isinstance(value, str):
        esc = value.replace("\\", "\\\\").replace('"', '\\"')
        return f'"{esc}"'
    if isinstance(value, list):
        return "[" + ", ".join(format_toml_value(v) for v in value) + "]"
    raise TypeError(f"unsupported TOML value type: {type(value)!r}")



def _write_table_lines(lines: List[str], prefix: str, table: Dict[str, Any]) -> None:
    scalar_items = []
    nested_items = []
    for key, value in table.items():
        if isinstance(value, dict):
            nested_items.append((key, value))
        else:
            scalar_items.append((key, value))
    if prefix:
        lines.append(f"[{prefix}]")
    for key, value in scalar_items:
        if value is None:
            continue
        lines.append(f"{key} = {format_toml_value(value)}")
    if scalar_items or prefix:
        lines.append("")
    for key, value in nested_items:
        _write_table_lines(lines, f"{prefix}.{key}" if prefix else key, value)



def write_beam_model_toml(output_path: Path, raw_config_path: Path, beam_rows_df, model_name: str) -> Path:
    raw = load_raw_toml(raw_config_path)
    beam_map = {str(row["beam_id"]): row for _, row in beam_rows_df.iterrows()}
    lines: List[str] = []

    top_level = {k: v for k, v in raw.items() if k != "spectrometers"}
    _write_table_lines(lines, "", top_level)

    spectrometers = list(raw.get("spectrometers", []) or [])
    for block in spectrometers:
        lines.append("[[spectrometers]]")
        beam_id = str(block.get("beam_id", f"B{int(block.get('fdnum', 0)):02d}"))
        beam_row = beam_map.get(beam_id)
        block_copy = dict(block)
        beam_block = dict(block_copy.get("beam", {}) or {})
        if beam_row is not None:
            beam_block.update({
                "az_offset_arcsec": float(beam_row["az_offset_arcsec"]),
                "el_offset_arcsec": float(beam_row["el_offset_arcsec"]),
                "rotation_mode": str(beam_row["rotation_mode"]),
                "reference_angle_deg": float(beam_row["reference_angle_deg"]),
                "rotation_sign": float(beam_row["rotation_sign"]),
                "dewar_angle_deg": float(beam_row["dewar_angle_deg"]),
                "beam_model_version": f"sunscan_multibeam_{model_name}",
            })
        block_copy["beam"] = beam_block

        scalar_keys = [
            "name", "fdnum", "ifnum", "plnum", "polariza", "beam_id",
            "frontend", "backend", "sampler", "db_stream_name", "db_table_name",
            "enabled", "use_for_convert", "use_for_sunscan", "use_for_fit", "beam_fit_use",
        ]
        for key in scalar_keys:
            if key in block_copy and block_copy[key] is not None:
                if key == "polariza" and str(block_copy[key]).upper() not in ALLOWED_POLARIZA:
                    raise ValueError(f"unsupported polariza={block_copy[key]!r} for stream {block_copy.get('name')!r}")
                lines.append(f"{key} = {format_toml_value(block_copy[key])}")
        lines.append("")
        for sub_key in ["frequency_axis", "local_oscillators", "override", "beam"]:
            sub = block_copy.get(sub_key)
            if isinstance(sub, dict) and sub:
                _write_table_lines(lines, f"spectrometers.{sub_key}", sub)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return output_path
