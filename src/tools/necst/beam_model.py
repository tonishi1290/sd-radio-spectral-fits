"""Shared NECST beam-model loader and override helpers.

PR2d foundation for the spectral-recording redesign.

Definitions
-----------
* ``beam_model.toml`` is the standalone geometry truth for converter,
  sunscan, and multibeam beam measurement override workflows.
* Streams reference beams by ``beam_id`` only.  This module validates that every
  stream beam_id is present in the selected beam model.
* External beam models are explicit overrides.  They must not rewrite the DB
  snapshot in place; callers should record the returned provenance in FITS
  HISTORY, JSON sidecars, or CSV manifests.

The module is dependency-light and intentionally avoids numpy/pandas/astropy so
it can be imported by CLI tests under ``python -S``.
"""

from __future__ import annotations

from dataclasses import dataclass
import copy
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:  # Python >= 3.11
    import tomllib as _toml_reader
except ModuleNotFoundError:  # pragma: no cover
    try:
        import tomli as _toml_reader  # type: ignore[no-redef]
    except ModuleNotFoundError:  # pragma: no cover
        _toml_reader = None  # type: ignore[assignment]


CANONICAL_BEAM_MODEL_NAME = "beam_model.toml"
PURE_ROTATION_MODEL = "pure_rotation_v1"


class BeamModelError(ValueError):
    """Raised for invalid beam-model or beam-id resolution state."""


@dataclass(frozen=True)
class BeamModelDocument:
    """Normalized beam-model document.

    ``beams`` is keyed by beam_id and contains only JSON/TOML scalar values that
    downstream code can copy into its own BeamConfig/LightBeam dataclass.
    """

    beams: Dict[str, Dict[str, Any]]
    source_path: Optional[Path] = None
    source_kind: str = "beam_model_toml"
    input_file_sha256: Optional[str] = None
    beam_model_version: Optional[str] = None
    warnings: Tuple[str, ...] = ()

    def as_dict(self) -> Dict[str, Any]:
        return {
            "beams": copy.deepcopy(self.beams),
            "source_path": None if self.source_path is None else str(self.source_path),
            "source_kind": self.source_kind,
            "input_file_sha256": self.input_file_sha256,
            "beam_model_version": self.beam_model_version,
            "warnings": list(self.warnings),
        }


def sha256_file(path: os.PathLike[str] | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _read_toml(path: Path) -> Dict[str, Any]:
    if _toml_reader is None:  # pragma: no cover
        raise BeamModelError("No TOML reader is available. Use Python >= 3.11, or install tomli.")
    text = path.read_text(encoding="utf-8")
    if hasattr(_toml_reader, "loads"):
        return _toml_reader.loads(text)  # type: ignore[no-any-return]
    raise BeamModelError("Configured TOML reader has no loads() function")


def _as_float(value: Any, *, key: str, default: Optional[float] = None) -> float:
    if value is None:
        if default is None:
            raise BeamModelError(f"beam field {key!r} is required")
        return float(default)
    try:
        out = float(value)
    except Exception as exc:
        raise BeamModelError(f"beam field {key!r} must be numeric, got {value!r}") from exc
    if not math.isfinite(out):
        raise BeamModelError(f"beam field {key!r} must be finite, got {value!r}")
    return out


def _normalize_rotation_mode(value: Any) -> str:
    text = str(value if value is not None else "none").strip().lower()
    if text in {"", "legacy", "fixed"}:
        return "none"
    if text not in {"none", "elevation", PURE_ROTATION_MODEL}:
        raise BeamModelError(f"rotation_mode must be none, elevation, or {PURE_ROTATION_MODEL!r}; got {value!r}")
    return text


def _normalize_model(value: Any, rotation_mode: str) -> str:
    text = str(value if value is not None else "legacy").strip().lower()
    if text in {"", "none"}:
        text = "legacy"
    if rotation_mode == PURE_ROTATION_MODEL:
        text = PURE_ROTATION_MODEL
    return text


def _normalize_pure_sign(value: Any) -> float:
    sign = _as_float(value, key="pure_rotation_sign/rotation_sign")
    if sign not in (-1.0, 1.0):
        raise BeamModelError(f"pure_rotation_v1 requires sign +/-1, got {value!r}")
    return sign


def normalize_beam_block(beam_id: str, raw_block: Mapping[str, Any]) -> Dict[str, Any]:
    """Normalize one ``[beams.<beam_id>]`` block.

    Supports both the standalone flat schema used by the v26 spec and the older
    nested ``pure_rotation`` block produced by legacy sunscan config exports.
    """

    if not str(beam_id).strip():
        raise BeamModelError("beam_id must be non-empty")
    block = dict(raw_block or {})
    rotation_mode = _normalize_rotation_mode(block.get("rotation_mode", "none"))
    model = _normalize_model(block.get("model", None), rotation_mode)
    version = block.get("beam_model_version", block.get("version", None))
    version_s = None if version is None else str(version).strip()

    if model == PURE_ROTATION_MODEL:
        nested = dict(block.get("pure_rotation", {}) or {})

        def _nonzero_legacy_value(key: str) -> bool:
            if key not in block or block.get(key) is None:
                return False
            return abs(_as_float(block.get(key), key=key)) > 0.0

        legacy_conflicts = [
            key for key in ("az_offset_arcsec", "el_offset_arcsec", "reference_angle_deg", "reference_el_deg")
            if _nonzero_legacy_value(key)
        ]
        if "rotation_slope_deg_per_deg" in block and block.get("rotation_slope_deg_per_deg") is not None:
            if abs(_as_float(block.get("rotation_slope_deg_per_deg"), key="rotation_slope_deg_per_deg")) > 0.0:
                legacy_conflicts.append("rotation_slope_deg_per_deg")
        if legacy_conflicts:
            raise BeamModelError(
                f"beam {beam_id!r}: {PURE_ROTATION_MODEL} cannot be mixed with non-zero legacy beam fields: "
                f"{sorted(set(legacy_conflicts))}"
            )

        def _first_consistent_float(key: str, *values: Any) -> Any:
            present = [v for v in values if v is not None]
            if not present:
                return None
            first = _as_float(present[0], key=key)
            for value in present[1:]:
                other = _as_float(value, key=key)
                if abs(other - first) > 1e-9:
                    raise BeamModelError(
                        f"beam {beam_id!r}: inconsistent {key} values in flat/nested pure_rotation schema: "
                        f"{present!r}"
                    )
            return first

        x0 = _first_consistent_float(
            "pure_rotation_offset_x_el0_arcsec",
            block.get("pure_rotation_offset_x_el0_arcsec"),
            block.get("offset_x_el0_arcsec"),
            nested.get("offset_x_el0_arcsec"),
        )
        y0 = _first_consistent_float(
            "pure_rotation_offset_y_el0_arcsec",
            block.get("pure_rotation_offset_y_el0_arcsec"),
            block.get("offset_y_el0_arcsec"),
            nested.get("offset_y_el0_arcsec"),
        )
        sign = _first_consistent_float(
            "pure_rotation_sign/rotation_sign",
            block.get("pure_rotation_sign"),
            block.get("rotation_sign"),
            nested.get("rotation_sign"),
        )
        return {
            "beam_id": str(beam_id),
            "model": PURE_ROTATION_MODEL,
            "rotation_mode": PURE_ROTATION_MODEL,
            "beam_model_version": version_s,
            "az_offset_arcsec": 0.0,
            "el_offset_arcsec": 0.0,
            "reference_angle_deg": 0.0,
            "rotation_sign": _normalize_pure_sign(sign),
            "rotation_slope_deg_per_deg": None,
            "dewar_angle_deg": _as_float(block.get("dewar_angle_deg", 0.0), key="dewar_angle_deg", default=0.0),
            "pure_rotation_offset_x_el0_arcsec": _as_float(x0, key="pure_rotation_offset_x_el0_arcsec"),
            "pure_rotation_offset_y_el0_arcsec": _as_float(y0, key="pure_rotation_offset_y_el0_arcsec"),
            "pure_rotation_sign": _normalize_pure_sign(sign),
        }

    return {
        "beam_id": str(beam_id),
        "model": model,
        "rotation_mode": rotation_mode,
        "beam_model_version": version_s,
        "az_offset_arcsec": _as_float(block.get("az_offset_arcsec", 0.0), key="az_offset_arcsec", default=0.0),
        "el_offset_arcsec": _as_float(block.get("el_offset_arcsec", 0.0), key="el_offset_arcsec", default=0.0),
        "reference_angle_deg": _as_float(block.get("reference_angle_deg", block.get("reference_el_deg", 0.0)), key="reference_angle_deg", default=0.0),
        "rotation_sign": _as_float(block.get("rotation_sign", 1.0), key="rotation_sign", default=1.0),
        "rotation_slope_deg_per_deg": (
            _as_float(block.get("rotation_slope_deg_per_deg"), key="rotation_slope_deg_per_deg")
            if block.get("rotation_slope_deg_per_deg") is not None else None
        ),
        "dewar_angle_deg": _as_float(block.get("dewar_angle_deg", 0.0), key="dewar_angle_deg", default=0.0),
        "pure_rotation_offset_x_el0_arcsec": None,
        "pure_rotation_offset_y_el0_arcsec": None,
        "pure_rotation_sign": None,
    }


def _beams_from_legacy_spectrometers(raw: Mapping[str, Any]) -> Tuple[Dict[str, Dict[str, Any]], Tuple[str, ...]]:
    beams: Dict[str, Dict[str, Any]] = {}
    warnings: List[str] = []
    for idx, block in enumerate(list(raw.get("spectrometers", []) or [])):
        if not isinstance(block, Mapping):
            continue
        beam_id = str(block.get("beam_id", f"B{int(block.get('fdnum', idx)):02d}")).strip()
        if not beam_id:
            raise BeamModelError(f"legacy spectrometers[{idx}] has empty beam_id")
        beam_block = dict(block.get("beam", {}) or {})
        normalized = normalize_beam_block(beam_id, beam_block)
        previous = beams.get(beam_id)
        if previous is not None and previous != normalized:
            raise BeamModelError(f"legacy config contains conflicting beam geometry for beam_id={beam_id!r}")
        beams[beam_id] = normalized
    if beams:
        warnings.append("beam model was derived from legacy [[spectrometers]].beam blocks")
    return beams, tuple(warnings)


def normalize_beam_model(raw: Mapping[str, Any], *, source_path: Optional[Path] = None, source_kind: str = "beam_model_toml", input_file_sha256: Optional[str] = None) -> BeamModelDocument:
    raw = dict(raw or {})
    warnings: List[str] = []
    if isinstance(raw.get("beams"), Mapping):
        beams_raw = dict(raw.get("beams", {}) or {})
        beams = {str(beam_id): normalize_beam_block(str(beam_id), block) for beam_id, block in beams_raw.items()}
    else:
        beams, legacy_warnings = _beams_from_legacy_spectrometers(raw)
        warnings.extend(legacy_warnings)
    if not beams:
        raise BeamModelError("beam model contains no [beams.<beam_id>] entries")
    versions = sorted({str(b.get("beam_model_version")) for b in beams.values() if b.get("beam_model_version")})
    version = str(raw.get("beam_model_version", raw.get("version", ""))).strip() or (versions[0] if len(versions) == 1 else None)
    return BeamModelDocument(
        beams=beams,
        source_path=source_path,
        source_kind=source_kind,
        input_file_sha256=input_file_sha256,
        beam_model_version=version,
        warnings=tuple(warnings),
    )


def load_beam_model(path_or_db_dir: os.PathLike[str] | str) -> BeamModelDocument:
    """Load a standalone beam model or discover ``beam_model.toml`` in a DB dir.

    If a directory is supplied, only the DB-sidecar ``beam_model.toml`` is read.
    The DB snapshot is not modified and is not used as an implicit replacement
    for a missing standalone beam model.
    """

    path = Path(path_or_db_dir).expanduser().resolve()
    if path.is_dir():
        candidate = path / CANONICAL_BEAM_MODEL_NAME
        if not candidate.exists():
            raise BeamModelError(f"No {CANONICAL_BEAM_MODEL_NAME} found in {path}")
        path = candidate
    if not path.exists():
        raise BeamModelError(f"beam model file does not exist: {path}")
    raw = _read_toml(path)
    return normalize_beam_model(raw, source_path=path, source_kind="beam_model_toml", input_file_sha256=sha256_file(path))


def beam_model_from_snapshot(snapshot: Mapping[str, Any]) -> BeamModelDocument:
    return normalize_beam_model({"beams": dict(snapshot.get("beams", {}) or {})}, source_path=None, source_kind="snapshot_embedded_beams")


def _stream_name(entry_key: str, stream_entry: Mapping[str, Any]) -> str:
    return str(stream_entry.get("stream_id", entry_key)).strip() or str(entry_key)


def validate_beam_ids(stream_table: Mapping[str, Any] | Sequence[Any], beam_model: BeamModelDocument | Mapping[str, Any]) -> List[str]:
    """Validate that every stream beam_id exists in a beam model.

    Returns sorted used beam_ids.  Raises ``BeamModelError`` for missing/empty
    beam_id entries.
    """

    bm = beam_model if isinstance(beam_model, BeamModelDocument) else normalize_beam_model(beam_model)
    beams = set(str(k) for k in bm.beams.keys())
    missing: List[str] = []
    used: List[str] = []

    if isinstance(stream_table, Mapping):
        iterable = list(stream_table.items())
        for key, entry_any in iterable:
            entry = dict(entry_any or {}) if isinstance(entry_any, Mapping) else {}
            beam_id = str(entry.get("beam_id", "")).strip()
            name = _stream_name(str(key), entry)
            if not beam_id:
                missing.append(f"{name}:<empty>")
            elif beam_id not in beams:
                missing.append(f"{name}:{beam_id}")
            else:
                used.append(beam_id)
    else:
        for idx, stream in enumerate(stream_table):
            beam = getattr(stream, "beam", None)
            beam_id = str(getattr(beam, "beam_id", getattr(stream, "beam_id", ""))).strip()
            name = str(getattr(stream, "name", getattr(stream, "stream_id", idx)))
            if not beam_id:
                missing.append(f"{name}:<empty>")
            elif beam_id not in beams:
                missing.append(f"{name}:{beam_id}")
            else:
                used.append(beam_id)
    if missing:
        raise BeamModelError("stream beam_id(s) not found in beam model: " + ", ".join(sorted(missing)))
    return sorted(set(used))


def resolve_stream_beam_geometry(stream_table: Mapping[str, Any] | Sequence[Any], beam_model: BeamModelDocument | Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    bm = beam_model if isinstance(beam_model, BeamModelDocument) else normalize_beam_model(beam_model)
    validate_beam_ids(stream_table, bm)
    out: Dict[str, Dict[str, Any]] = {}
    if isinstance(stream_table, Mapping):
        for key, entry_any in stream_table.items():
            entry = dict(entry_any or {}) if isinstance(entry_any, Mapping) else {}
            name = _stream_name(str(key), entry)
            beam_id = str(entry.get("beam_id", "")).strip()
            out[name] = copy.deepcopy(bm.beams[beam_id])
    else:
        for idx, stream in enumerate(stream_table):
            name = str(getattr(stream, "name", getattr(stream, "stream_id", idx)))
            beam = getattr(stream, "beam", None)
            beam_id = str(getattr(beam, "beam_id", getattr(stream, "beam_id", ""))).strip()
            out[name] = copy.deepcopy(bm.beams[beam_id])
    return out


def apply_beam_model_to_snapshot(snapshot: Mapping[str, Any], beam_model: BeamModelDocument | Mapping[str, Any]) -> Dict[str, Any]:
    bm = beam_model if isinstance(beam_model, BeamModelDocument) else normalize_beam_model(beam_model)
    out = copy.deepcopy(dict(snapshot))
    validate_beam_ids(dict(out.get("streams", {}) or {}), bm)
    out["beams"] = copy.deepcopy(bm.beams)
    prov = dict(out.get("provenance", {}) or {})
    prov["beam_model_source"] = bm.source_kind
    if bm.source_path is not None:
        prov["beam_model_path"] = str(bm.source_path)
    if bm.input_file_sha256:
        prov["beam_model_input_file_sha256"] = bm.input_file_sha256
    prov["beam_model_override_used"] = bool(bm.source_kind != "snapshot_embedded_beams")
    out["provenance"] = prov
    return out


def beam_model_provenance_record(beam_model: BeamModelDocument | Mapping[str, Any], *, override_used: bool) -> Dict[str, Any]:
    bm = beam_model if isinstance(beam_model, BeamModelDocument) else normalize_beam_model(beam_model)
    return {
        "beam_model_source_kind": bm.source_kind,
        "beam_model_source_path": None if bm.source_path is None else str(bm.source_path),
        "beam_model_input_file_sha256": bm.input_file_sha256,
        "beam_model_version": bm.beam_model_version,
        "beam_model_override_used": bool(override_used),
        "beam_ids": sorted(bm.beams.keys()),
        "warnings": list(bm.warnings),
    }


def write_beam_model_provenance(output: Any, beam_model: BeamModelDocument | Mapping[str, Any], *, source: Optional[str] = None, override_used: bool) -> Dict[str, Any]:
    """Write provenance to common output objects and return the record.

    Supported outputs:
    * object with ``add_history(key, value)`` method, e.g. SDFITS writer
    * list, receiving ``(key, value)`` tuples
    * dict, receiving string values by key
    """

    record = beam_model_provenance_record(beam_model, override_used=override_used)
    if source is not None:
        record["beam_model_provenance_source"] = str(source)
    pairs = [
        ("beam_model_source", record["beam_model_source_kind"]),
        ("beam_model_path", record.get("beam_model_source_path") or ""),
        ("beam_model_sha256", record.get("beam_model_input_file_sha256") or ""),
        ("beam_model_version", record.get("beam_model_version") or ""),
        ("beam_model_override", str(bool(record.get("beam_model_override_used"))).lower()),
        ("beam_model_ids", ",".join(record.get("beam_ids", []))),
    ]
    if hasattr(output, "add_history"):
        for key, value in pairs:
            output.add_history(key, value)
    elif isinstance(output, list):
        output.extend(pairs)
    elif isinstance(output, dict):
        for key, value in pairs:
            output[key] = value
    else:
        raise TypeError("output must provide add_history(), or be a list/dict")
    return record


def dumps_beam_model(beam_model: BeamModelDocument | Mapping[str, Any]) -> str:
    bm = beam_model if isinstance(beam_model, BeamModelDocument) else normalize_beam_model(beam_model)
    lines: List[str] = []
    if bm.beam_model_version:
        lines.append(f"beam_model_version = {json.dumps(str(bm.beam_model_version))}")
        lines.append("")
    for beam_id in sorted(bm.beams):
        block = bm.beams[beam_id]
        lines.append(f"[beams.{beam_id}]")
        for key in [
            "beam_model_version", "model", "az_offset_arcsec", "el_offset_arcsec", "rotation_mode",
            "reference_angle_deg", "rotation_sign", "rotation_slope_deg_per_deg", "dewar_angle_deg",
            "pure_rotation_offset_x_el0_arcsec", "pure_rotation_offset_y_el0_arcsec", "pure_rotation_sign",
        ]:
            value = block.get(key)
            if value is None:
                continue
            if isinstance(value, str):
                text = json.dumps(value)
            elif isinstance(value, bool):
                text = "true" if value else "false"
            else:
                text = repr(float(value)) if isinstance(value, float) else str(value)
            lines.append(f"{key} = {text}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
