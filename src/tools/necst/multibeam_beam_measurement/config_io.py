from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence
import math
import tomllib

import pandas as pd

from .legacy_loaders import load_converter_module


ALLOWED_POLARIZA = {"RR", "LL", "RL", "LR", "XX", "YY", "XY", "YX"}

PURE_ROTATION_MODEL = "pure_rotation_v1"
PURE_ROTATION_LEGACY_KEYS = {
    "az_offset_arcsec",
    "el_offset_arcsec",
    "rotation_mode",
    "reference_angle_deg",
    "reference_el_deg",
    "rotation_sign",
    "rotation_slope_deg_per_deg",
    "dewar_angle_deg",
}


def _normalize_beam_model(value: Any) -> str:
    if value is None:
        return "legacy"
    s = str(value).strip()
    if not s:
        return "legacy"
    if s == PURE_ROTATION_MODEL:
        return PURE_ROTATION_MODEL
    return s


def _normalize_rotation_mode(value: Any) -> str:
    s = str(value if value is not None else "none").strip().lower()
    if s not in {"none", "elevation"}:
        raise ValueError(f"rotation_mode must be one of ['none', 'elevation'], got {value!r}")
    return s


def _build_light_beam(beam_id: str, beam_block: Dict[str, Any]) -> LightBeam:
    beam_block = dict(beam_block or {})
    beam_model = _normalize_beam_model(beam_block.get("model", None))
    beam_model_version = (str(beam_block.get("beam_model_version")).strip() if beam_block.get("beam_model_version") is not None else None)
    if beam_model == PURE_ROTATION_MODEL:
        legacy_keys = sorted(k for k in PURE_ROTATION_LEGACY_KEYS if k in beam_block and beam_block.get(k) is not None)
        if legacy_keys:
            raise ValueError(
                f"beam.model={PURE_ROTATION_MODEL!r} forbids legacy beam keys: {legacy_keys}"
            )
        pure_block = dict(beam_block.get("pure_rotation", {}) or {})
        for req in ("offset_x_el0_arcsec", "offset_y_el0_arcsec", "rotation_sign"):
            if pure_block.get(req) is None:
                raise ValueError(f"beam.model={PURE_ROTATION_MODEL!r} requires beam.pure_rotation.{req}")
        pure_sign = float(pure_block.get("rotation_sign"))
        if pure_sign not in (-1.0, 1.0):
            raise ValueError(f"beam.pure_rotation.rotation_sign must be +/-1, got {pure_sign!r}")
        return LightBeam(
            beam_id=beam_id,
            beam_model=PURE_ROTATION_MODEL,
            beam_model_version=beam_model_version,
            rotation_sign=pure_sign,
            pure_rotation_offset_x_el0_arcsec=float(pure_block.get("offset_x_el0_arcsec")),
            pure_rotation_offset_y_el0_arcsec=float(pure_block.get("offset_y_el0_arcsec")),
            pure_rotation_sign=pure_sign,
        )
    return LightBeam(
        beam_id=beam_id,
        beam_model=beam_model,
        beam_model_version=beam_model_version,
        az_offset_arcsec=float(beam_block.get("az_offset_arcsec", 0.0)),
        el_offset_arcsec=float(beam_block.get("el_offset_arcsec", 0.0)),
        rotation_mode=_normalize_rotation_mode(beam_block.get("rotation_mode", "none")),
        reference_angle_deg=float(beam_block.get("reference_angle_deg", beam_block.get("reference_el_deg", 0.0))),
        rotation_sign=float(beam_block.get("rotation_sign", 1.0)),
        rotation_slope_deg_per_deg=(
            float(beam_block.get("rotation_slope_deg_per_deg"))
            if beam_block.get("rotation_slope_deg_per_deg") is not None else None
        ),
        dewar_angle_deg=float(beam_block.get("dewar_angle_deg", 0.0)),
    )



@dataclass
class LightBeam:
    beam_id: str
    beam_model: str = "legacy"
    beam_model_version: Optional[str] = None
    az_offset_arcsec: float = 0.0
    el_offset_arcsec: float = 0.0
    rotation_mode: str = "none"
    reference_angle_deg: float = 0.0
    rotation_sign: float = 1.0
    rotation_slope_deg_per_deg: Optional[float] = None
    dewar_angle_deg: float = 0.0
    pure_rotation_offset_x_el0_arcsec: Optional[float] = None
    pure_rotation_offset_y_el0_arcsec: Optional[float] = None
    pure_rotation_sign: Optional[float] = None


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
        beam = _build_light_beam(beam_id, beam_block)
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
                "beam_model": str(getattr(stream.beam, "beam_model", "legacy")),
                "beam_model_version": getattr(stream.beam, "beam_model_version", None),
                "rotation_mode": str(getattr(stream.beam, "rotation_mode", "none")),
                "reference_angle_deg": float(getattr(stream.beam, "reference_angle_deg", 0.0)),
                "rotation_sign": float(getattr(stream.beam, "rotation_sign", 1.0)),
                "rotation_slope_deg_per_deg": (
                    float(getattr(stream.beam, "rotation_slope_deg_per_deg"))
                    if getattr(stream.beam, "rotation_slope_deg_per_deg", None) is not None else None
                ),
                "dewar_angle_deg": float(getattr(stream.beam, "dewar_angle_deg", 0.0)),
                "az_offset_arcsec": float(getattr(stream.beam, "az_offset_arcsec", 0.0)),
                "el_offset_arcsec": float(getattr(stream.beam, "el_offset_arcsec", 0.0)),
                "offset_x_el0_arcsec": getattr(stream.beam, "pure_rotation_offset_x_el0_arcsec", None),
                "offset_y_el0_arcsec": getattr(stream.beam, "pure_rotation_offset_y_el0_arcsec", None),
                "pure_rotation_sign": getattr(stream.beam, "pure_rotation_sign", None),
            }
        )
    return pd.DataFrame(rows)



def resolve_primary_stream_names(raw_config_path: Path, stream_cfg: Dict[str, Any], explicit_stream_names: Optional[Sequence[str]] = None) -> List[str]:
    """Return fit-enabled stream names.

    The historical name is kept for CLI compatibility, but pure_rotation_v1 does
    not require a single primary stream per beam. Duplicate beam_id entries are
    therefore allowed and all selected fit streams are returned.
    """
    streams = filter_streams_for_purpose(raw_config_path, stream_cfg, "fit", explicit_stream_names=explicit_stream_names)
    if explicit_stream_names:
        return [str(s) for s in explicit_stream_names]
    return [str(getattr(stream, "name", "")) for stream in streams]



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
    if duplicate_beam_ids:
        for beam_id in duplicate_beam_ids:
            dup = table.loc[table["beam_id"].astype(str) == str(beam_id)].copy()
            pure_mask = dup.get("beam_model", pd.Series(["legacy"] * len(dup), index=dup.index)).astype(str) == PURE_ROTATION_MODEL
            if bool(pure_mask.all()):
                pure = dup.loc[pure_mask].copy()
                sig = pd.to_numeric(pure["pure_rotation_sign"], errors="coerce")
                x0 = pd.to_numeric(pure["offset_x_el0_arcsec"], errors="coerce")
                y0 = pd.to_numeric(pure["offset_y_el0_arcsec"], errors="coerce")
                if not (sig.nunique(dropna=True) <= 1 and x0.nunique(dropna=True) <= 1 and y0.nunique(dropna=True) <= 1):
                    warnings.append(
                        f"beam_id={beam_id!r} appears in multiple streams but pure_rotation_v1 parameters are not identical across those streams"
                    )
            elif bool(pure_mask.any()):
                warnings.append(
                    f"beam_id={beam_id!r} mixes pure_rotation_v1 and legacy beam models across streams"
                )

    for _, row in table.iterrows():
        if str(row["polariza"]).upper() not in ALLOWED_POLARIZA:
            warnings.append(f"unsupported polariza for stream {row['stream_name']!r}: {row['polariza']!r}")
        restfreq = row.get("restfreq_hz")
        if restfreq is None or (isinstance(restfreq, float) and not math.isfinite(restfreq)):
            warnings.append(f"stream {row['stream_name']!r} has no finite restfreq_hz")
        rot_sign = row.get("pure_rotation_sign") if str(row.get("beam_model", "legacy")) == PURE_ROTATION_MODEL else row.get("rotation_sign")
        try:
            rot_sign_f = float(rot_sign)
            allowed = (-1.0, 1.0) if str(row.get("beam_model", "legacy")) == PURE_ROTATION_MODEL else (-1.0, 0.0, 1.0)
            if rot_sign_f not in allowed:
                warnings.append(
                    f"stream {row['stream_name']!r} has rotation sign={rot_sign_f!r}; "
                    f"beam model {row.get('beam_model', 'legacy')!r} expects one of {allowed}"
                )
        except Exception:
            warnings.append(f"stream {row['stream_name']!r} has non-numeric rotation sign={rot_sign!r}")

    if primary_streams:
        primary_table = table.loc[table["stream_name"].astype(str).isin([str(s) for s in primary_streams])].copy()
        if not primary_table.empty:
            geom_mag = (
                pd.to_numeric(primary_table["az_offset_arcsec"], errors="coerce").fillna(0.0).abs()
                + pd.to_numeric(primary_table["el_offset_arcsec"], errors="coerce").fillna(0.0).abs()
                + pd.to_numeric(primary_table.get("offset_x_el0_arcsec"), errors="coerce").fillna(0.0).abs()
                + pd.to_numeric(primary_table.get("offset_y_el0_arcsec"), errors="coerce").fillna(0.0).abs()
            )
            if bool((geom_mag <= 0.0).all()):
                warnings.append(
                    "all selected fit streams have zero nominal beam offsets; pseudo multi-beam dry-runs will not constrain rotation unless non-zero beam offsets are set in [spectrometers.beam]"
                )
            beam_models = primary_table.get("beam_model", pd.Series(["legacy"] * len(primary_table), index=primary_table.index)).astype(str).str.strip()
            pure_mask = beam_models == PURE_ROTATION_MODEL
            if bool(pure_mask.any()):
                pure_signs = pd.to_numeric(primary_table.loc[pure_mask, "pure_rotation_sign"], errors="coerce")
                if bool((pure_signs.abs() != 1.0).any()):
                    warnings.append("some selected fit streams use pure_rotation_v1 but do not have rotation_sign=+/-1")
            elif "rotation_mode" in primary_table.columns:
                rot_modes = primary_table["rotation_mode"].astype(str).str.lower().str.strip()
                if bool((rot_modes == "none").all()):
                    warnings.append(
                        "all selected fit streams have rotation_mode='none'; pseudo multi-beam dry-runs will not exercise EL-dependent beam rotation"
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
    missing_required: List[str] = []

    top_level = {k: v for k, v in raw.items() if k != "spectrometers"}
    _write_table_lines(lines, "", top_level)

    spectrometers = list(raw.get("spectrometers", []) or [])
    for block in spectrometers:
        lines.append("[[spectrometers]]")
        beam_id = str(block.get("beam_id", f"B{int(block.get('fdnum', 0)):02d}"))
        beam_row = beam_map.get(beam_id)
        block_copy = dict(block)
        beam_block = dict(block_copy.get("beam", {}) or {})
        usage_policy = resolve_stream_usage_policy(block)
        requires_geometry = bool(usage_policy.get("use_for_convert", True) or usage_policy.get("use_for_fit", True))
        if beam_row is not None:
            has_pure = all(col in beam_row.index for col in ["offset_x_el0_arcsec", "offset_y_el0_arcsec", "rotation_sign"])
            if has_pure and pd.notna(beam_row["offset_x_el0_arcsec"]) and pd.notna(beam_row["offset_y_el0_arcsec"]):
                for legacy_key in [
                    "az_offset_arcsec", "el_offset_arcsec", "rotation_mode", "reference_angle_deg", "reference_el_deg",
                    "rotation_sign", "rotation_slope_deg_per_deg", "dewar_angle_deg",
                ]:
                    beam_block.pop(legacy_key, None)
                beam_block["model"] = PURE_ROTATION_MODEL
                beam_block["beam_model_version"] = f"sunscan_multibeam_{model_name}"
                beam_block["pure_rotation"] = {
                    "offset_x_el0_arcsec": float(beam_row["offset_x_el0_arcsec"]),
                    "offset_y_el0_arcsec": float(beam_row["offset_y_el0_arcsec"]),
                    "rotation_sign": float(beam_row["rotation_sign"]),
                }
            else:
                beam_block.update({
                    "az_offset_arcsec": float(beam_row["az_offset_arcsec"]),
                    "el_offset_arcsec": float(beam_row["el_offset_arcsec"]),
                    "rotation_mode": str(beam_row["rotation_mode"]),
                    "reference_angle_deg": float(beam_row["reference_angle_deg"]),
                    "rotation_sign": float(beam_row["rotation_sign"]),
                    "dewar_angle_deg": float(beam_row["dewar_angle_deg"]),
                    "beam_model_version": f"sunscan_multibeam_{model_name}",
                })
                if "rotation_slope_deg_per_deg" in beam_row and pd.notna(beam_row["rotation_slope_deg_per_deg"]):
                    beam_block["rotation_slope_deg_per_deg"] = float(beam_row["rotation_slope_deg_per_deg"])
        elif requires_geometry:
            missing_required.append(f"stream={block.get('name', '')!r} beam_id={beam_id!r}")
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

    if missing_required:
        raise ValueError(
            "cannot write beam model TOML because fit results are missing required beam geometry for: "
            + ", ".join(missing_required)
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return output_path
