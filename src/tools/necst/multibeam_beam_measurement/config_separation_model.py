from __future__ import annotations

"""Config-separation internal model and legacy adapter for NECST/XFFTS tools.

This module is intentionally side-effect free and does not import heavy runtime
dependencies such as necstdb, astropy, numpy, or pandas.  It provides the first
implementation step of the converter/sunscan config separation:

* read the current all-in-one legacy ``--spectrometer-config`` TOML;
* normalize it into a purpose-separated ``ResolvedConfigBundle``;
* materialize the bundle back into the existing converter ``StreamConfig`` or
  sunscan ``LightStream`` objects.

The existing converter/sunscan entry points are not changed by this module.
Keeping this layer independent makes it possible to add regression tests before
switching the tools to the new config stack.
"""

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import copy
import importlib.util
import math
import re
import sys

try:
    import tomllib  # Python 3.11+
except Exception:  # pragma: no cover
    tomllib = None  # type: ignore[assignment]


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


@dataclass
class ResolvedBeamModel:
    """Beam geometry separated from stream selection and fit policy."""

    beam_id: str
    model: str = "legacy"
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
    raw: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ResolvedStream:
    """Stream truth needed to bind DB rows to physical/SDFITS metadata."""

    stream_id: str
    db_stream_name: str
    fdnum: int
    ifnum: int
    plnum: int
    polariza: str
    beam_id: str
    stream_index: int = 0
    aliases: List[str] = field(default_factory=list)
    spectrometer_key: Optional[str] = None
    board_id: Optional[int] = None
    db_table_name: Optional[str] = None
    sampler: Optional[str] = None
    frontend: Optional[str] = None
    backend: Optional[str] = None
    frequency_axis: Dict[str, Any] = field(default_factory=dict)
    local_oscillators: Dict[str, Any] = field(default_factory=dict)
    override: Dict[str, Any] = field(default_factory=dict)
    channel_slice_spec: Optional[Any] = None
    raw: Dict[str, Any] = field(default_factory=dict)


@dataclass
class StreamSelectionPolicy:
    """Legacy use flags normalized as analysis-time selection policy."""

    enabled: bool = True
    use_for_convert: bool = True
    use_for_sunscan: bool = True
    use_for_fit: bool = True
    beam_fit_use: Optional[bool] = None


@dataclass
class AnalysisStreamSelection:
    """Analysis stream selection derived from old use_for_* flags."""

    policies: Dict[str, StreamSelectionPolicy] = field(default_factory=dict)
    convert_stream_ids: List[str] = field(default_factory=list)
    sunscan_extract_stream_ids: List[str] = field(default_factory=list)
    sunscan_fit_stream_ids: List[str] = field(default_factory=list)

    def streams_for(self, purpose: str) -> List[str]:
        key = str(purpose).strip().lower()
        if key in {"convert", "converter"}:
            return list(self.convert_stream_ids)
        if key in {"sunscan", "sunscan_extract", "extract"}:
            return list(self.sunscan_extract_stream_ids)
        if key in {"fit", "sunscan_fit"}:
            return list(self.sunscan_fit_stream_ids)
        raise ValueError(f"unsupported stream-selection purpose: {purpose!r}")


@dataclass
class DbTimeEnvironmentConfig:
    """Global DB/time/weather/temperature bindings preserved from legacy config."""

    values: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ConverterAnalysisConfig:
    """Converter-specific analysis settings separated from stream truth.

    In the legacy all-in-one config, per-stream ``channel_slice`` means an
    export-time converter slice.  It must not be confused with an observation
    recording window.
    """

    values: Dict[str, Any] = field(default_factory=dict)
    export_slices_by_stream_id: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SunScanAnalysisConfigModel:
    """Sunscan-specific analysis settings and per-stream legacy overrides."""

    values: Dict[str, Any] = field(default_factory=dict)
    per_stream_overrides: Dict[str, Dict[str, Any]] = field(default_factory=dict)


@dataclass
class ResolvedConfigBundle:
    """Purpose-separated internal model produced by config adapters."""

    schema_version: int
    config_name: Optional[str] = None
    config_description: Optional[str] = None
    source_format: str = "legacy_spectrometer_config"
    source_path: Optional[str] = None
    global_config: Dict[str, Any] = field(default_factory=dict)
    provenance: Dict[str, Any] = field(default_factory=dict)
    streams: List[ResolvedStream] = field(default_factory=list)
    beams: Dict[str, ResolvedBeamModel] = field(default_factory=dict)
    selection: AnalysisStreamSelection = field(default_factory=AnalysisStreamSelection)
    db_time_environment: DbTimeEnvironmentConfig = field(default_factory=DbTimeEnvironmentConfig)
    converter_analysis: ConverterAnalysisConfig = field(default_factory=ConverterAnalysisConfig)
    sunscan_analysis: SunScanAnalysisConfigModel = field(default_factory=SunScanAnalysisConfigModel)

    def stream_by_id(self) -> Dict[str, ResolvedStream]:
        return {s.stream_id: s for s in self.streams}

    def alias_map(self) -> Dict[str, str]:
        out: Dict[str, str] = {}
        for stream in self.streams:
            for alias in stream.aliases:
                if alias:
                    out[str(alias)] = stream.stream_id
        return out

    def resolve_stream_id(self, stream_name_or_alias: str) -> str:
        key = str(stream_name_or_alias)
        amap = self.alias_map()
        if key not in amap:
            raise KeyError(f"unknown stream id/alias: {stream_name_or_alias!r}")
        return amap[key]


def load_raw_toml(path: Path | str) -> Dict[str, Any]:
    if tomllib is None:  # pragma: no cover
        raise RuntimeError("tomllib is required; use Python 3.11 or newer")
    with open(Path(path), "rb") as fh:
        data = tomllib.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"TOML root must be a table: {path}")
    return data


def _nonempty_str(value: Any, default: Optional[str] = None) -> Optional[str]:
    if value is None:
        return default
    s = str(value).strip()
    return s if s else default


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


def resolve_legacy_selection_policy(block: Mapping[str, Any]) -> StreamSelectionPolicy:
    """Resolve old enabled/use_for_*/beam_fit_use flags without changing behavior."""

    enabled_raw = _coerce_optional_bool(block.get("enabled"))
    enabled = True if enabled_raw is None else bool(enabled_raw)

    def _purpose_flag(key: str) -> bool:
        raw = _coerce_optional_bool(block.get(key))
        if raw is None:
            return bool(enabled)
        return bool(enabled) and bool(raw)

    beam_fit_use = _coerce_optional_bool(block.get("beam_fit_use"))
    use_for_fit_raw = _coerce_optional_bool(block.get("use_for_fit"))
    if use_for_fit_raw is None:
        use_for_fit = bool(enabled) if beam_fit_use is None else bool(enabled) and bool(beam_fit_use)
    else:
        use_for_fit = bool(enabled) and bool(use_for_fit_raw)

    return StreamSelectionPolicy(
        enabled=bool(enabled),
        use_for_convert=_purpose_flag("use_for_convert"),
        use_for_sunscan=_purpose_flag("use_for_sunscan"),
        use_for_fit=bool(use_for_fit),
        beam_fit_use=beam_fit_use,
    )


def _normalize_rotation_mode(value: Any) -> str:
    s = str(value if value is not None else "none").strip().lower()
    if s not in {"none", "elevation"}:
        raise ValueError(f"rotation_mode must be one of ['none', 'elevation'], got {value!r}")
    return s


def _normalize_beam_model(value: Any) -> str:
    s = str(value if value is not None else "legacy").strip()
    return s or "legacy"


def _normalize_pure_rotation_sign(value: Any) -> float:
    sign = float(value)
    if sign not in (-1.0, 1.0):
        raise ValueError(f"pure_rotation rotation_sign must be +/-1, got {value!r}")
    return sign


def _build_resolved_beam(beam_id: str, beam_block: Mapping[str, Any]) -> ResolvedBeamModel:
    block = dict(beam_block or {})
    model = _normalize_beam_model(block.get("model"))
    rotation_mode_raw = block.get("rotation_mode")
    if str(rotation_mode_raw if rotation_mode_raw is not None else "").strip().lower() == PURE_ROTATION_MODEL:
        model = PURE_ROTATION_MODEL
    version = _nonempty_str(block.get("beam_model_version"))
    if model == PURE_ROTATION_MODEL:
        # Accept both pure_rotation schemas:
        #   1. legacy nested form: [beam.pure_rotation] offset_x_el0_arcsec = ...
        #   2. normalized flat form produced by tools.necst.beam_model:
        #      pure_rotation_offset_x_el0_arcsec = ...
        #
        # The normalized flat form deliberately also contains harmless legacy
        # compatibility fields such as rotation_mode, rotation_sign, and
        # dewar_angle_deg.  These must not be treated as conflicts.  Only
        # non-zero legacy az/el/reference offsets or an explicit legacy slope
        # would make the interpretation ambiguous.
        pure = dict(block.get("pure_rotation", {}) or {})

        def _first_consistent_float(key: str, *values: Any) -> Any:
            present = [v for v in values if v is not None]
            if not present:
                return None
            try:
                first = float(present[0])
            except Exception as exc:
                raise ValueError(f"beam_id={beam_id!r}: {key} must be numeric, got {present[0]!r}") from exc
            for value in present[1:]:
                try:
                    other = float(value)
                except Exception as exc:
                    raise ValueError(f"beam_id={beam_id!r}: {key} must be numeric, got {value!r}") from exc
                if abs(other - first) > 1e-9:
                    raise ValueError(
                        f"beam_id={beam_id!r}: inconsistent {key} values in flat/nested pure_rotation schema: "
                        f"{present!r}"
                    )
            return first

        x0 = _first_consistent_float(
            "pure_rotation_offset_x_el0_arcsec",
            block.get("pure_rotation_offset_x_el0_arcsec"),
            block.get("offset_x_el0_arcsec"),
            pure.get("offset_x_el0_arcsec"),
        )
        y0 = _first_consistent_float(
            "pure_rotation_offset_y_el0_arcsec",
            block.get("pure_rotation_offset_y_el0_arcsec"),
            block.get("offset_y_el0_arcsec"),
            pure.get("offset_y_el0_arcsec"),
        )
        sign_raw = _first_consistent_float(
            "pure_rotation_sign/rotation_sign",
            block.get("pure_rotation_sign"),
            block.get("rotation_sign"),
            pure.get("rotation_sign"),
        )
        for req, value in (
            ("pure_rotation_offset_x_el0_arcsec", x0),
            ("pure_rotation_offset_y_el0_arcsec", y0),
            ("pure_rotation_sign/rotation_sign", sign_raw),
        ):
            if value is None:
                raise ValueError(f"beam_id={beam_id!r}: {req} is required for {PURE_ROTATION_MODEL}")
        legacy_conflicts: List[str] = []
        for key in ("az_offset_arcsec", "el_offset_arcsec", "reference_angle_deg", "reference_el_deg"):
            value = block.get(key)
            if value is None:
                continue
            try:
                if abs(float(value)) > 0.0:
                    legacy_conflicts.append(key)
            except Exception:
                legacy_conflicts.append(key)
        slope = block.get("rotation_slope_deg_per_deg")
        if slope is not None:
            try:
                if abs(float(slope)) > 0.0:
                    legacy_conflicts.append("rotation_slope_deg_per_deg")
            except Exception:
                legacy_conflicts.append("rotation_slope_deg_per_deg")
        if legacy_conflicts:
            raise ValueError(
                f"beam_id={beam_id!r}: {PURE_ROTATION_MODEL!r} cannot be mixed with non-zero legacy beam keys: "
                f"{sorted(set(legacy_conflicts))}"
            )
        sign = _normalize_pure_rotation_sign(sign_raw)
        return ResolvedBeamModel(
            beam_id=str(beam_id),
            model=PURE_ROTATION_MODEL,
            beam_model_version=version,
            rotation_mode=PURE_ROTATION_MODEL,
            rotation_sign=sign,
            dewar_angle_deg=float(block.get("dewar_angle_deg", 0.0)),
            pure_rotation_offset_x_el0_arcsec=float(x0),
            pure_rotation_offset_y_el0_arcsec=float(y0),
            pure_rotation_sign=sign,
            raw=copy.deepcopy(block),
        )
    return ResolvedBeamModel(
        beam_id=str(beam_id),
        model=model,
        beam_model_version=version,
        az_offset_arcsec=float(block.get("az_offset_arcsec", 0.0)),
        el_offset_arcsec=float(block.get("el_offset_arcsec", 0.0)),
        rotation_mode=_normalize_rotation_mode(block.get("rotation_mode", "none")),
        reference_angle_deg=float(block.get("reference_angle_deg", block.get("reference_el_deg", 0.0))),
        rotation_sign=float(block.get("rotation_sign", 1.0)),
        rotation_slope_deg_per_deg=(
            float(block.get("rotation_slope_deg_per_deg"))
            if block.get("rotation_slope_deg_per_deg") is not None
            else None
        ),
        dewar_angle_deg=float(block.get("dewar_angle_deg", 0.0)),
        raw=copy.deepcopy(block),
    )


def _beam_signature(beam: ResolvedBeamModel) -> Dict[str, Any]:
    data = asdict(beam)
    data.pop("raw", None)
    return data


def _infer_board_id(block: Mapping[str, Any], db_stream_name: str) -> Optional[int]:
    raw = block.get("board_id", None)
    if raw is not None:
        return int(raw)
    match = re.search(r"(?:board|xffts-board|xffts_board)[^\d]*(\d+)", str(db_stream_name), flags=re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None


def _unique_aliases(*items: Any) -> List[str]:
    out: List[str] = []
    seen = set()
    for item in items:
        if item is None:
            continue
        if isinstance(item, (list, tuple)):
            values = item
        else:
            values = [item]
        for value in values:
            s = str(value).strip()
            if s and s not in seen:
                out.append(s)
                seen.add(s)
    return out


def _extract_channel_slice(block: Mapping[str, Any], frequency_axis: Dict[str, Any]) -> Optional[Any]:
    channel_slice = block.get("channel_slice", None)
    if channel_slice is None and "channel_slice" in frequency_axis:
        channel_slice = frequency_axis.pop("channel_slice")
    return channel_slice


def _finite_positive_float_or_nan(raw: Any, *, context: str) -> float:
    if raw is None:
        return float("nan")
    try:
        value = float(raw)
    except Exception as exc:
        raise ValueError(f"{context} must be a positive finite scalar, got {raw!r}") from exc
    if (not math.isfinite(value)) or value <= 0.0:
        raise ValueError(f"{context} must be positive and finite, got {raw!r}")
    return value


def _extract_restfreq_hz_from_frequency_axis(freq_cfg: Mapping[str, Any], *, context: str) -> float:
    hz_raw = freq_cfg.get("restfreq_hz", None)
    ghz_raw = freq_cfg.get("restfreq_ghz", None)
    hz = _finite_positive_float_or_nan(hz_raw, context=f"{context}.restfreq_hz")
    ghz = _finite_positive_float_or_nan(ghz_raw, context=f"{context}.restfreq_ghz")
    hz_from_ghz = ghz * 1.0e9 if math.isfinite(ghz) else float("nan")

    if (hz_raw is not None) and (ghz_raw is not None):
        atol = max(1.0e-6, abs(hz) * 1.0e-12)
        if not math.isclose(hz, hz_from_ghz, rel_tol=1.0e-12, abs_tol=atol):
            raise ValueError(
                f"{context} has inconsistent restfreq_hz={hz_raw!r} and restfreq_ghz={ghz_raw!r}"
            )
        return float(hz)
    if hz_raw is not None:
        return float(hz)
    if ghz_raw is not None:
        return float(hz_from_ghz)
    return float("nan")


def _canonicalize_frequency_axis_restfreq(freq_cfg: Mapping[str, Any], *, context: str) -> Dict[str, Any]:
    out = dict(freq_cfg or {})
    restfreq_hz = _extract_restfreq_hz_from_frequency_axis(out, context=context)
    out.pop("restfreq_ghz", None)
    if math.isfinite(restfreq_hz):
        out["restfreq_hz"] = float(restfreq_hz)
    else:
        out.pop("restfreq_hz", None)
    return out


def load_legacy_spectrometer_config_bundle(
    config_path: Path | str,
    *,
    strict_converter_identity: bool = False,
) -> ResolvedConfigBundle:
    """Load old all-in-one ``--spectrometer-config`` into the new internal model.

    Parameters
    ----------
    strict_converter_identity:
        When True, enforce the converter's historical duplicate
        ``(fdnum, ifnum, plnum, polariza, sampler)`` check.  The default is
        False because the sunscan fallback loader historically did not enforce
        this check; converter integration can enable it explicitly.
    """

    path = Path(config_path)
    raw = load_raw_toml(path)
    version_raw = raw.get("schema_version", raw.get("config_version", 1))
    schema_version = int(version_raw)
    if schema_version != 1:
        raise ValueError(f"unsupported schema_version/config_version={schema_version!r}; expected 1")

    blocks = list(raw.get("spectrometers", []) or [])
    if not blocks:
        raise ValueError("legacy spectrometer config requires at least one [[spectrometers]] block")

    streams: List[ResolvedStream] = []
    beams: Dict[str, ResolvedBeamModel] = {}
    policies: Dict[str, StreamSelectionPolicy] = {}
    convert_ids: List[str] = []
    sunscan_ids: List[str] = []
    fit_ids: List[str] = []
    export_slices: Dict[str, Any] = {}
    per_stream_overrides: Dict[str, Dict[str, Any]] = {}

    seen_stream_ids = set()
    seen_sdfits_keys = set()

    for idx, raw_block in enumerate(blocks):
        block = dict(raw_block or {})
        stream_id = _nonempty_str(block.get("name"))
        if stream_id is None:
            raise ValueError(f"[[spectrometers]] block index={idx} has no non-empty name")
        if stream_id in seen_stream_ids:
            raise ValueError(f"duplicate stream id/name: {stream_id!r}")
        seen_stream_ids.add(stream_id)

        polariza = str(block.get("polariza", "")).strip().upper()
        if polariza not in ALLOWED_POLARIZA:
            raise ValueError(f"stream {stream_id!r}: unsupported polariza={polariza!r}")
        fdnum = int(block.get("fdnum", 0))
        ifnum = int(block.get("ifnum", 0))
        plnum = int(block.get("plnum", 0))
        sampler = _nonempty_str(block.get("sampler"))
        sdfits_key = (fdnum, ifnum, plnum, polariza, sampler)
        if strict_converter_identity and sdfits_key in seen_sdfits_keys:
            raise ValueError(f"duplicate (fdnum, ifnum, plnum, polariza, sampler) combination: {sdfits_key!r}")
        seen_sdfits_keys.add(sdfits_key)

        db_stream_name = _nonempty_str(block.get("db_stream_name"), default=stream_id)
        assert db_stream_name is not None
        beam_id = _nonempty_str(block.get("beam_id"), default=f"B{fdnum:02d}")
        assert beam_id is not None
        beam = _build_resolved_beam(beam_id, dict(block.get("beam", {}) or {}))
        if beam_id in beams:
            if _beam_signature(beams[beam_id]) != _beam_signature(beam):
                raise ValueError(f"beam_id={beam_id!r} has inconsistent definitions across streams")
        else:
            beams[beam_id] = beam

        frequency_axis = copy.deepcopy(dict(block.get("frequency_axis", {}) or {}))
        channel_slice = _extract_channel_slice(block, frequency_axis)
        frequency_axis = _canonicalize_frequency_axis_restfreq(
            frequency_axis,
            context=f"spectrometers[{idx}].frequency_axis",
        )
        local_oscillators = copy.deepcopy(dict(block.get("local_oscillators", {}) or {}))
        override = copy.deepcopy(dict(block.get("override", {}) or {}))

        policy = resolve_legacy_selection_policy(block)
        policies[stream_id] = policy
        if policy.use_for_convert:
            convert_ids.append(stream_id)
        if policy.use_for_sunscan:
            sunscan_ids.append(stream_id)
        if policy.use_for_fit:
            fit_ids.append(stream_id)
        if channel_slice is not None:
            export_slices[stream_id] = copy.deepcopy(channel_slice)
        if override:
            per_stream_overrides[stream_id] = copy.deepcopy(override)

        aliases = _unique_aliases(stream_id, db_stream_name, block.get("legacy_aliases"))
        streams.append(
            ResolvedStream(
                stream_id=stream_id,
                db_stream_name=db_stream_name,
                db_table_name=_nonempty_str(block.get("db_table_name")),
                fdnum=fdnum,
                ifnum=ifnum,
                plnum=plnum,
                polariza=polariza,
                beam_id=beam_id,
                stream_index=int(idx),
                aliases=aliases,
                spectrometer_key=_nonempty_str(block.get("spectrometer_key")),
                board_id=_infer_board_id(block, db_stream_name),
                sampler=sampler,
                frontend=_nonempty_str(block.get("frontend")),
                backend=_nonempty_str(block.get("backend")),
                frequency_axis=frequency_axis,
                local_oscillators=local_oscillators,
                override=override,
                channel_slice_spec=copy.deepcopy(channel_slice),
                raw=copy.deepcopy(block),
            )
        )

    global_cfg = copy.deepcopy(dict(raw.get("global", {}) or {}))
    provenance = copy.deepcopy(dict(raw.get("provenance", {}) or {}))
    selection = AnalysisStreamSelection(
        policies=policies,
        convert_stream_ids=convert_ids,
        sunscan_extract_stream_ids=sunscan_ids,
        sunscan_fit_stream_ids=fit_ids,
    )

    return ResolvedConfigBundle(
        schema_version=schema_version,
        config_name=_nonempty_str(raw.get("config_name")),
        config_description=_nonempty_str(raw.get("config_description")),
        source_format="legacy_spectrometer_config",
        source_path=str(path),
        global_config=global_cfg,
        provenance=provenance,
        streams=streams,
        beams=beams,
        selection=selection,
        db_time_environment=DbTimeEnvironmentConfig(values=copy.deepcopy(global_cfg)),
        converter_analysis=ConverterAnalysisConfig(
            values=copy.deepcopy(global_cfg),
            export_slices_by_stream_id=export_slices,
        ),
        sunscan_analysis=SunScanAnalysisConfigModel(
            values=copy.deepcopy(global_cfg),
            per_stream_overrides=per_stream_overrides,
        ),
    )


def _restfreq_hz_from_frequency_axis(frequency_axis: Mapping[str, Any]) -> float:
    for key, factor in (("restfreq_hz", 1.0), ("rest_frequency_hz", 1.0), ("restfreq_ghz", 1.0e9), ("rest_frequency_ghz", 1.0e9)):
        if key in frequency_axis and frequency_axis.get(key) is not None:
            try:
                val = float(frequency_axis.get(key)) * factor
                if math.isfinite(val) and val > 0:
                    return val
            except Exception:
                pass
    return float("nan")


def materialize_converter_config(bundle: ResolvedConfigBundle, converter_module: Any) -> Dict[str, Any]:
    """Build the current converter config dict from ``ResolvedConfigBundle``.

    ``converter_module`` is passed in explicitly so this module does not import
    the heavy converter at import time.  The object must provide ``BeamConfig``,
    ``StreamConfig`` and ``derive_stream_wcs``.
    """

    streams = []
    for resolved in bundle.streams:
        beam_model = bundle.beams[resolved.beam_id]
        beam = converter_module.BeamConfig(
            beam_id=beam_model.beam_id,
            model=beam_model.model,
            az_offset_arcsec=beam_model.az_offset_arcsec,
            el_offset_arcsec=beam_model.el_offset_arcsec,
            rotation_mode=beam_model.rotation_mode,
            reference_angle_deg=beam_model.reference_angle_deg,
            rotation_sign=beam_model.rotation_sign,
            rotation_slope_deg_per_deg=beam_model.rotation_slope_deg_per_deg,
            dewar_angle_deg=beam_model.dewar_angle_deg,
            beam_model_version=beam_model.beam_model_version,
            pure_rotation_offset_x_el0_arcsec=beam_model.pure_rotation_offset_x_el0_arcsec,
            pure_rotation_offset_y_el0_arcsec=beam_model.pure_rotation_offset_y_el0_arcsec,
            pure_rotation_sign=beam_model.pure_rotation_sign,
        )
        policy = bundle.selection.policies.get(resolved.stream_id, StreamSelectionPolicy())
        stream = converter_module.StreamConfig(
            name=resolved.stream_id,
            fdnum=int(resolved.fdnum),
            ifnum=int(resolved.ifnum),
            plnum=int(resolved.plnum),
            polariza=str(resolved.polariza),
            beam=beam,
            enabled=bool(policy.enabled),
            use_for_convert=bool(policy.use_for_convert),
            frontend=resolved.frontend,
            backend=resolved.backend,
            sampler=resolved.sampler,
            db_stream_name=resolved.db_stream_name,
            db_table_name=resolved.db_table_name,
            frequency_axis=copy.deepcopy(resolved.frequency_axis),
            local_oscillators=copy.deepcopy(resolved.local_oscillators),
            override=copy.deepcopy(resolved.override),
            channel_slice_spec=copy.deepcopy(resolved.channel_slice_spec),
            stream_index=int(resolved.stream_index),
        )
        stream.wcs_full = converter_module.derive_stream_wcs(stream)
        stream.wcs = stream.wcs_full
        streams.append(stream)

    return {
        "schema_version": int(bundle.schema_version),
        "config_name": bundle.config_name,
        "config_description": bundle.config_description,
        "global": copy.deepcopy(bundle.global_config),
        "provenance": copy.deepcopy(bundle.provenance),
        "streams": streams,
    }


def materialize_light_config(bundle: ResolvedConfigBundle, config_io_module: Any) -> Dict[str, Any]:
    """Build the current sunscan/config_io style dict from ``ResolvedConfigBundle``."""

    streams = []
    for resolved in bundle.streams:
        beam_model = bundle.beams[resolved.beam_id]
        beam = config_io_module.LightBeam(
            beam_id=beam_model.beam_id,
            beam_model=beam_model.model,
            beam_model_version=beam_model.beam_model_version,
            az_offset_arcsec=beam_model.az_offset_arcsec,
            el_offset_arcsec=beam_model.el_offset_arcsec,
            rotation_mode=beam_model.rotation_mode,
            reference_angle_deg=beam_model.reference_angle_deg,
            rotation_sign=beam_model.rotation_sign,
            rotation_slope_deg_per_deg=beam_model.rotation_slope_deg_per_deg,
            dewar_angle_deg=beam_model.dewar_angle_deg,
            pure_rotation_offset_x_el0_arcsec=beam_model.pure_rotation_offset_x_el0_arcsec,
            pure_rotation_offset_y_el0_arcsec=beam_model.pure_rotation_offset_y_el0_arcsec,
            pure_rotation_sign=beam_model.pure_rotation_sign,
        )
        policy = bundle.selection.policies.get(resolved.stream_id, StreamSelectionPolicy())
        wcs = config_io_module.LightWCS(restfreq_hz=_restfreq_hz_from_frequency_axis(resolved.frequency_axis))
        streams.append(
            config_io_module.LightStream(
                name=resolved.stream_id,
                fdnum=int(resolved.fdnum),
                ifnum=int(resolved.ifnum),
                plnum=int(resolved.plnum),
                polariza=str(resolved.polariza),
                beam=beam,
                db_stream_name=resolved.db_stream_name,
                db_table_name=resolved.db_table_name,
                sampler=resolved.sampler,
                frontend=resolved.frontend,
                backend=resolved.backend,
                frequency_axis=copy.deepcopy(resolved.frequency_axis),
                local_oscillators=copy.deepcopy(resolved.local_oscillators),
                override=copy.deepcopy(resolved.override),
                enabled=bool(policy.enabled),
                use_for_convert=bool(policy.use_for_convert),
                use_for_sunscan=bool(policy.use_for_sunscan),
                use_for_fit=bool(policy.use_for_fit),
                beam_fit_use=policy.beam_fit_use,
                wcs=wcs,
                wcs_full=wcs,
            )
        )

    return {
        "schema_version": int(bundle.schema_version),
        "config_name": bundle.config_name,
        "config_description": bundle.config_description,
        "global": copy.deepcopy(bundle.global_config),
        "provenance": copy.deepcopy(bundle.provenance),
        "streams": streams,
    }


def bundle_to_summary_dict(bundle: ResolvedConfigBundle) -> Dict[str, Any]:
    """Return a JSON/TOML-friendly summary useful for tests and diagnostics."""

    return {
        "schema_version": bundle.schema_version,
        "config_name": bundle.config_name,
        "source_format": bundle.source_format,
        "source_path": bundle.source_path,
        "streams": [
            {
                "stream_id": s.stream_id,
                "aliases": list(s.aliases),
                "db_stream_name": s.db_stream_name,
                "db_table_name": s.db_table_name,
                "board_id": s.board_id,
                "fdnum": s.fdnum,
                "ifnum": s.ifnum,
                "plnum": s.plnum,
                "polariza": s.polariza,
                "beam_id": s.beam_id,
                "channel_slice_spec": s.channel_slice_spec,
            }
            for s in bundle.streams
        ],
        "beams": {k: _beam_signature(v) for k, v in bundle.beams.items()},
        "selection": {
            "convert": list(bundle.selection.convert_stream_ids),
            "sunscan_extract": list(bundle.selection.sunscan_extract_stream_ids),
            "sunscan_fit": list(bundle.selection.sunscan_fit_stream_ids),
        },
        "export_slices_by_stream_id": copy.deepcopy(bundle.converter_analysis.export_slices_by_stream_id),
        "per_stream_overrides": copy.deepcopy(bundle.sunscan_analysis.per_stream_overrides),
    }




# -----------------------------------------------------------------------------
# Standalone analysis-config overlays
# -----------------------------------------------------------------------------

def _load_analysis_toml(path: Path | str, *, expected_kind: Optional[str] = None) -> Dict[str, Any]:
    """Load a converter/sunscan analysis TOML.

    The accepted files are intentionally shallow and explicit.  They do not
    redefine stream truth (LO, IF, beam identity, frequency axis); they only
    provide analysis defaults, export slices, stream-selection flags, and
    per-stream sunscan overrides.
    """

    raw = load_raw_toml(Path(path))
    kind = _nonempty_str(raw.get("kind"), default=None)
    if expected_kind is not None and kind is not None and kind != expected_kind:
        raise ValueError(f"{path}: kind={kind!r} is not {expected_kind!r}")
    return raw


def _merge_into_global(config_dict: Dict[str, Any], values: Mapping[str, Any]) -> None:
    if not values:
        return
    global_cfg = dict(config_dict.get("global", {}) or {})
    for key, value in dict(values).items():
        if value is not None:
            global_cfg[str(key)] = copy.deepcopy(value)
    config_dict["global"] = global_cfg


def _stream_map(config_dict: Dict[str, Any]) -> Dict[str, Any]:
    return {str(getattr(stream, "name", "")): stream for stream in list(config_dict.get("streams", []) or [])}


def _stream_slice_value(raw_value: Any) -> Any:
    if isinstance(raw_value, Mapping):
        for key in ("channel_slice", "export_slice", "slice", "value"):
            if raw_value.get(key) is not None:
                return copy.deepcopy(raw_value.get(key))
        raise ValueError(f"stream export slice table must contain channel_slice/export_slice/slice/value: {raw_value!r}")
    return copy.deepcopy(raw_value)


def _extract_export_slices(raw: Mapping[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for top_key in ("export_slices_by_stream_id", "export_slices", "channel_slices_by_stream_id", "channel_slices"):
        table = raw.get(top_key)
        if isinstance(table, Mapping):
            for stream_id, value in table.items():
                out[str(stream_id)] = _stream_slice_value(value)

    converter = raw.get("converter_analysis")
    if isinstance(converter, Mapping):
        for nested_key in ("export_slices_by_stream_id", "export_slices", "channel_slices_by_stream_id", "channel_slices"):
            table = converter.get(nested_key)
            if isinstance(table, Mapping):
                for stream_id, value in table.items():
                    out[str(stream_id)] = _stream_slice_value(value)
    return out


def _extract_analysis_values(raw: Mapping[str, Any], section_name: str) -> Dict[str, Any]:
    """Return scalar analysis values from [global] and [section_name].

    Nested tables that have special meaning are removed here and handled by the
    caller.  This keeps TOML files readable while preserving legacy code that
    still consumes a single ``config_dict["global"]`` namespace.
    """

    values: Dict[str, Any] = {}
    global_section = raw.get("global")
    if isinstance(global_section, Mapping):
        values.update(copy.deepcopy(dict(global_section)))

    section = raw.get(section_name)
    if isinstance(section, Mapping):
        for key, value in dict(section).items():
            if key in {
                "export_slices",
                "export_slices_by_stream_id",
                "channel_slices",
                "channel_slices_by_stream_id",
                "per_stream_overrides",
                "streams",
                "analysis_selection",
            }:
                continue
            values[str(key)] = copy.deepcopy(value)
    return values


_SELECTION_KEYS = {"enabled", "use_for_convert", "use_for_sunscan", "use_for_fit", "beam_fit_use"}


def _coerce_optional_bool_local(value: Any) -> Optional[bool]:
    if value is None:
        return None
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value)
    s = str(value).strip().lower()
    if s in {"1", "true", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "no", "n", "off"}:
        return False
    raise ValueError(f"cannot coerce to bool: {value!r}")


def _apply_stream_selection_table(config_dict: Dict[str, Any], selection_raw: Any) -> None:
    """Apply stream-selection flags from a TOML table.

    Accepted forms:

    [analysis_selection]
    convert_stream_ids = ["A"]
    sunscan_extract_stream_ids = ["A", "B"]
    sunscan_fit_stream_ids = ["A"]

    [analysis_selection.streams.A]
    enabled = true
    use_for_convert = true
    use_for_sunscan = false

    [streams.A]
    use_for_convert = true
    """

    if not isinstance(selection_raw, Mapping):
        return

    by_name = _stream_map(config_dict)

    def require_known(name: str) -> Any:
        if name not in by_name:
            raise ValueError(f"analysis selection references unknown stream_id={name!r}")
        return by_name[name]

    def set_bool_attr(stream: Any, attr: str, value: Any) -> None:
        flag = _coerce_optional_bool_local(value)
        if flag is not None:
            setattr(stream, attr, bool(flag))

    # Explicit lists are convenient for cookbook-style selection.
    list_specs = {
        "convert_stream_ids": "use_for_convert",
        "converter_stream_ids": "use_for_convert",
        "sunscan_extract_stream_ids": "use_for_sunscan",
        "sunscan_stream_ids": "use_for_sunscan",
        "sunscan_fit_stream_ids": "use_for_fit",
        "fit_stream_ids": "use_for_fit",
    }
    for key, attr in list_specs.items():
        if selection_raw.get(key) is None:
            continue
        wanted = {str(x) for x in list(selection_raw.get(key) or [])}
        for name, stream in by_name.items():
            setattr(stream, attr, name in wanted)

    streams_table = selection_raw.get("streams")
    if isinstance(streams_table, Mapping):
        for name, values in streams_table.items():
            stream = require_known(str(name))
            if not isinstance(values, Mapping):
                raise ValueError(f"analysis_selection.streams.{name} must be a table")
            for key, value in values.items():
                if key in _SELECTION_KEYS:
                    set_bool_attr(stream, str(key), value)

    # Also allow [analysis_selection.<stream_id>] for compact files, but ignore
    # known list keys and scalar metadata.
    for name, values in selection_raw.items():
        if name in set(list_specs) | {"streams"}:
            continue
        if not isinstance(values, Mapping):
            continue
        stream = require_known(str(name))
        for key, value in values.items():
            if key in _SELECTION_KEYS:
                set_bool_attr(stream, str(key), value)


def _apply_top_level_stream_selection(config_dict: Dict[str, Any], raw: Mapping[str, Any]) -> None:
    if isinstance(raw.get("analysis_selection"), Mapping):
        _apply_stream_selection_table(config_dict, raw.get("analysis_selection"))

    streams_table = raw.get("streams")
    if isinstance(streams_table, Mapping):
        by_name = _stream_map(config_dict)
        for name, values in streams_table.items():
            if not isinstance(values, Mapping):
                continue
            if not any(key in values for key in _SELECTION_KEYS):
                continue
            if str(name) not in by_name:
                raise ValueError(f"streams.{name} selection references unknown stream_id={name!r}")
            _apply_stream_selection_table(config_dict, {"streams": {str(name): dict(values)}})


def apply_converter_analysis_config(
    config_dict: Dict[str, Any],
    analysis_path: Path | str,
) -> Dict[str, Any]:
    """Overlay a standalone converter-analysis TOML onto materialized config.

    This implements the public separation promised by the handover: stream truth
    remains in spectrometer/snapshot-derived config, while converter-only
    analysis settings live in a separate file.
    """

    raw = _load_analysis_toml(analysis_path, expected_kind="converter_analysis")
    _merge_into_global(config_dict, _extract_analysis_values(raw, "converter_analysis"))

    converter = raw.get("converter_analysis") if isinstance(raw.get("converter_analysis"), Mapping) else {}
    if isinstance(converter, Mapping):
        for key in ("channel_slice", "export_slice"):
            if converter.get(key) is not None:
                config_dict.setdefault("global", {})["channel_slice"] = copy.deepcopy(converter.get(key))

    stream_slices = _extract_export_slices(raw)
    if stream_slices:
        by_name = _stream_map(config_dict)
        for stream_id, slice_spec in stream_slices.items():
            if stream_id not in by_name:
                raise ValueError(f"{analysis_path}: export slice references unknown stream_id={stream_id!r}")
            setattr(by_name[stream_id], "channel_slice_spec", copy.deepcopy(slice_spec))

    _apply_top_level_stream_selection(config_dict, raw)

    provenance = dict(config_dict.get("provenance", {}) or {})
    provenance["converter_analysis_config_path"] = str(Path(analysis_path))
    provenance["converter_analysis_config_kind"] = "converter_analysis"
    config_dict["provenance"] = provenance
    return config_dict


def apply_sunscan_analysis_config(
    config_dict: Dict[str, Any],
    analysis_path: Path | str,
) -> Dict[str, Any]:
    """Overlay a standalone sunscan-analysis TOML onto materialized config."""

    raw = _load_analysis_toml(analysis_path, expected_kind="sunscan_analysis")
    _merge_into_global(config_dict, _extract_analysis_values(raw, "sunscan_analysis"))

    def merge_override(stream_id: str, values: Mapping[str, Any]) -> None:
        by_name = _stream_map(config_dict)
        if stream_id not in by_name:
            raise ValueError(f"{analysis_path}: per-stream override references unknown stream_id={stream_id!r}")
        stream = by_name[stream_id]
        override = dict(getattr(stream, "override", {}) or {})
        for key, value in dict(values).items():
            if key in _SELECTION_KEYS:
                continue
            override[str(key)] = copy.deepcopy(value)
        setattr(stream, "override", override)

    for top_key in ("per_stream_overrides", "stream_overrides"):
        table = raw.get(top_key)
        if isinstance(table, Mapping):
            for stream_id, values in table.items():
                if not isinstance(values, Mapping):
                    raise ValueError(f"{analysis_path}: {top_key}.{stream_id} must be a table")
                merge_override(str(stream_id), values)

    section = raw.get("sunscan_analysis")
    if isinstance(section, Mapping):
        for nested_key in ("per_stream_overrides", "stream_overrides"):
            table = section.get(nested_key)
            if isinstance(table, Mapping):
                for stream_id, values in table.items():
                    if not isinstance(values, Mapping):
                        raise ValueError(f"{analysis_path}: sunscan_analysis.{nested_key}.{stream_id} must be a table")
                    merge_override(str(stream_id), values)

    streams_table = raw.get("streams")
    if isinstance(streams_table, Mapping):
        for stream_id, values in streams_table.items():
            if not isinstance(values, Mapping):
                continue
            # Selection flags are handled separately; remaining keys are
            # per-stream sunscan overrides.
            non_selection = {str(k): copy.deepcopy(v) for k, v in values.items() if k not in _SELECTION_KEYS}
            if non_selection:
                merge_override(str(stream_id), non_selection)

    _apply_top_level_stream_selection(config_dict, raw)

    provenance = dict(config_dict.get("provenance", {}) or {})
    provenance["sunscan_analysis_config_path"] = str(Path(analysis_path))
    provenance["sunscan_analysis_config_kind"] = "sunscan_analysis"
    config_dict["provenance"] = provenance
    return config_dict


def apply_analysis_stream_selection_config(
    config_dict: Dict[str, Any],
    selection_path: Path | str,
) -> Dict[str, Any]:
    """Apply standalone stream-selection TOML to converter/sunscan config."""

    raw = _load_analysis_toml(selection_path, expected_kind="analysis_stream_selection")
    if isinstance(raw.get("analysis_selection"), Mapping):
        _apply_stream_selection_table(config_dict, raw.get("analysis_selection"))
    else:
        _apply_stream_selection_table(config_dict, raw)

    if isinstance(raw.get("streams"), Mapping):
        _apply_top_level_stream_selection(config_dict, raw)

    provenance = dict(config_dict.get("provenance", {}) or {})
    paths = list(provenance.get("analysis_stream_selection_config_paths", []) or [])
    paths.append(str(Path(selection_path)))
    provenance["analysis_stream_selection_config_paths"] = paths
    config_dict["provenance"] = provenance
    return config_dict

# -----------------------------------------------------------------------------
# New-DB spectral-recording snapshot adapter
# -----------------------------------------------------------------------------

def _load_necst_helper_module(module_basename: str) -> Any:
    """Load a sibling ``tools.necst`` helper module robustly.

    ``config_separation_model.py`` is imported in three modes in practice:
    package-relative imports, absolute ``tools.necst`` imports, and direct
    file-based imports from the converter.  Keeping helper import fallback here
    avoids adding heavy dependencies or requiring the caller to know the import
    mode.
    """

    errors: List[str] = []
    try:
        from .. import spectral_recording_snapshot as _srsnap  # type: ignore
        from .. import beam_model as _beam_model  # type: ignore
        from .. import analysis_stream_selection as _selection  # type: ignore
        mapping = {
            "spectral_recording_snapshot": _srsnap,
            "beam_model": _beam_model,
            "analysis_stream_selection": _selection,
        }
        if module_basename in mapping:
            return mapping[module_basename]
    except Exception as exc:  # pragma: no cover - depends on import mode
        errors.append(f"relative: {exc}")

    try:
        if module_basename == "spectral_recording_snapshot":
            from tools.necst import spectral_recording_snapshot as mod  # type: ignore
        elif module_basename == "beam_model":
            from tools.necst import beam_model as mod  # type: ignore
        elif module_basename == "analysis_stream_selection":
            from tools.necst import analysis_stream_selection as mod  # type: ignore
        else:
            raise ValueError(module_basename)
        return mod
    except Exception as exc:  # pragma: no cover - depends on import mode
        errors.append(f"absolute: {exc}")

    helper_path = Path(__file__).resolve().parents[1] / f"{module_basename}.py"
    try:
        spec = importlib.util.spec_from_file_location(f"_necst_{module_basename}", str(helper_path))
        if spec is None or spec.loader is None:
            raise RuntimeError(f"cannot create import spec for {helper_path}")
        mod = importlib.util.module_from_spec(spec)
        sys.modules.setdefault(spec.name, mod)
        spec.loader.exec_module(mod)
        return mod
    except Exception as exc:  # pragma: no cover - fallback path
        errors.append(f"direct: {exc}")

    raise RuntimeError(f"failed to load helper module {module_basename!r}: " + " | ".join(errors))


def _snapshot_numeric(raw: Any, *, default: float = float("nan")) -> float:
    if raw is None:
        return default
    try:
        out = float(raw)
    except Exception as exc:
        raise ValueError(f"snapshot numeric field must be float-compatible, got {raw!r}") from exc
    if not math.isfinite(out):
        return default
    return out


def _snapshot_frequency_alias_hz(entry: Mapping[str, Any], keys_and_factors: Sequence[Tuple[str, float]]) -> float:
    values: List[Tuple[str, float]] = []
    for key, factor in keys_and_factors:
        if key in entry and entry.get(key) is not None:
            val = _snapshot_numeric(entry.get(key), default=float("nan")) * factor
            if math.isfinite(val) and val > 0.0:
                values.append((key, val))
    if not values:
        return float("nan")
    ref_key, ref_val = values[0]
    for key, val in values[1:]:
        tol = max(1.0e-6, abs(ref_val) * 1.0e-12)
        if abs(val - ref_val) > tol:
            raise ValueError(
                f"snapshot has inconsistent frequency aliases: "
                f"{ref_key}={ref_val} Hz but {key}={val} Hz"
            )
    return float(ref_val)


def _snapshot_axis_to_legacy_frequency_axis(
    *,
    stream_id: str,
    stream_entry: Mapping[str, Any],
    layout: Any,
    snapshot: Mapping[str, Any],
) -> Dict[str, Any]:
    """Convert snapshot frequency-axis truth into legacy explicit-WCS form.

    The new snapshot records the full input axis and the saved channel window.
    Converter/sunscan receive already-saved spectra, so the materialized axis is
    local to the saved data: local channel 0 corresponds to full channel
    ``saved_ch_start``.
    """

    axis_id = _nonempty_str(stream_entry.get("frequency_axis_id"))
    axes = snapshot.get("frequency_axes", {})
    if not axis_id or not isinstance(axes, Mapping) or axis_id not in axes:
        raise ValueError(f"snapshot stream {stream_id!r} has no valid frequency_axis_id")
    axis = dict(axes[axis_id] or {})

    full_start = int(getattr(layout, "saved_ch_start"))
    saved_nchan = int(getattr(layout, "saved_nchan"))
    if saved_nchan <= 0:
        raise ValueError(f"snapshot stream {stream_id!r} has invalid saved_nchan={saved_nchan!r}")

    if "sky_freq_at_full_ch0_hz" in axis and "sky_freq_step_hz" in axis:
        full_ch0_hz = _snapshot_numeric(axis.get("sky_freq_at_full_ch0_hz"))
        step_hz = _snapshot_numeric(axis.get("sky_freq_step_hz"))
        axis_source = "snapshot_sky_frequency"
    else:
        if0 = _snapshot_numeric(axis.get("if_freq_at_full_ch0_hz"))
        if_step = _snapshot_numeric(axis.get("if_freq_step_hz"))
        lo_chain_id = _nonempty_str(stream_entry.get("lo_chain"))
        lo_chains = snapshot.get("lo_chains", {})
        signed_sum_hz = None
        if lo_chain_id and isinstance(lo_chains, Mapping) and lo_chain_id in lo_chains:
            signed_sum_hz = _snapshot_numeric(dict(lo_chains[lo_chain_id]).get("signed_lo_sum_hz"))
        if_sign = 1.0
        if signed_sum_hz is None or not math.isfinite(signed_sum_hz):
            signed_sum_hz = 0.0
            axis_source = "snapshot_if_frequency_no_lo_sum"
        else:
            axis_source = "snapshot_signed_lo_sum_plus_if"
            try:
                if_sign = float(dict(lo_chains[lo_chain_id]).get("if_frequency_sign", 1.0))
            except Exception:
                if_sign = 1.0
        full_ch0_hz = signed_sum_hz + if_sign * if0
        step_hz = if_sign * if_step

    if not math.isfinite(full_ch0_hz) or not math.isfinite(step_hz) or step_hz == 0.0:
        raise ValueError(f"snapshot stream {stream_id!r} has invalid frequency axis")
    crval1_hz = float(full_ch0_hz + full_start * step_hz)

    restfreq_hz = _snapshot_frequency_alias_hz(
        stream_entry,
        (
            ("rest_frequency_hz", 1.0),
            ("restfreq_hz", 1.0),
            ("default_rest_frequency_hz", 1.0),
            ("rest_frequency_ghz", 1.0e9),
            ("restfreq_ghz", 1.0e9),
            ("default_rest_frequency_ghz", 1.0e9),
            ("rest_frequency_mhz", 1.0e6),
            ("restfreq_mhz", 1.0e6),
            ("default_rest_frequency_mhz", 1.0e6),
        ),
    )
    if not math.isfinite(restfreq_hz):
        restfreq_hz = _snapshot_frequency_alias_hz(
            axis,
            (
                ("rest_frequency_hz", 1.0),
                ("restfreq_hz", 1.0),
                ("rest_frequency_ghz", 1.0e9),
                ("restfreq_ghz", 1.0e9),
                ("rest_frequency_mhz", 1.0e6),
                ("restfreq_mhz", 1.0e6),
            ),
        )

    out = {
        "definition_mode": "explicit_wcs",
        "nchan": saved_nchan,
        "crval1_hz": crval1_hz,
        "cdelt1_hz": float(step_hz),
        "crpix1": 1.0,
        "ctype1": str(axis.get("ctype1", "FREQ")),
        "cunit1": str(axis.get("cunit1", "Hz")),
        "specsys": str(axis.get("specsys", "TOPOCENT")),
        "veldef": str(axis.get("veldef", "RADIO")),
        "restfreq_hz": restfreq_hz,
        "store_freq_column": axis.get("store_freq_column", "auto"),
        "snapshot_frequency_axis_id": axis_id,
        "snapshot_full_nchan": int(getattr(layout, "full_nchan")),
        "snapshot_saved_ch_start": int(getattr(layout, "saved_ch_start")),
        "snapshot_saved_ch_stop": int(getattr(layout, "saved_ch_stop")),
        "snapshot_axis_source": axis_source,
    }
    return _canonicalize_frequency_axis_restfreq(out, context=f"snapshot stream {stream_id!r} frequency_axis")


def _snapshot_lo_chain_to_local_oscillators(
    *,
    stream_entry: Mapping[str, Any],
    snapshot: Mapping[str, Any],
) -> Dict[str, Any]:
    """Return lightweight LO provenance for materialized legacy configs.

    The actual spectral axis is explicit-WCS, so these fields are not used to
    re-resolve channel frequencies.  They are retained only as analysis
    provenance where they can be represented safely.
    """

    lo_chain_id = _nonempty_str(stream_entry.get("lo_chain"))
    if not lo_chain_id:
        return {}
    lo_chains = snapshot.get("lo_chains", {})
    if not isinstance(lo_chains, Mapping) or lo_chain_id not in lo_chains:
        return {"lo_chain": lo_chain_id}
    chain = dict(lo_chains[lo_chain_id] or {})
    out: Dict[str, Any] = {
        "lo_chain": lo_chain_id,
        "formula_version": chain.get("formula_version"),
        "signed_lo_sum_hz": chain.get("signed_lo_sum_hz"),
        "if_frequency_sign": chain.get("if_frequency_sign"),
        "lo_roles": list(chain.get("lo_roles", []) or []),
        "signs": list(chain.get("signs", []) or []),
        "physical_lo_frequencies_hz": list(chain.get("physical_lo_frequencies_hz", []) or []),
    }
    legacy = chain.get("legacy_local_oscillators")
    if isinstance(legacy, Mapping):
        out.update(dict(legacy))
    return {k: v for k, v in out.items() if v is not None}



def _beam_position_signature_for_pointing(raw: Mapping[str, Any]) -> Tuple[Any, ...]:
    """Return geometry fields that affect pointing-reference correction."""

    d = dict(raw or {})
    model = str(d.get("model", d.get("beam_model", "legacy")) or "legacy")
    rotation_mode = str(d.get("rotation_mode", "none") or "none")
    if model == PURE_ROTATION_MODEL or rotation_mode == PURE_ROTATION_MODEL:
        fields = [
            ("model", PURE_ROTATION_MODEL),
            ("rotation_mode", PURE_ROTATION_MODEL),
            ("pure_rotation_offset_x_el0_arcsec", None),
            ("pure_rotation_offset_y_el0_arcsec", None),
            ("pure_rotation_sign", None),
            ("dewar_angle_deg", 0.0),
        ]
    else:
        fields = [
            ("model", model),
            ("rotation_mode", rotation_mode),
            ("az_offset_arcsec", 0.0),
            ("el_offset_arcsec", 0.0),
            ("reference_angle_deg", 0.0),
            ("reference_el_deg", None),
            ("rotation_sign", 1.0),
            ("rotation_slope_deg_per_deg", None),
            ("dewar_angle_deg", 0.0),
        ]
    out: List[Any] = []
    for key, default in fields:
        value = d.get(key, default)
        if isinstance(value, (int, float)):
            out.append(round(float(value), 9))
        elif value is None:
            out.append(None)
        else:
            out.append(str(value))
    return tuple(out)

def _validate_pointing_reference_beam_override(
    *,
    snapshot: Mapping[str, Any],
    external_beams: Mapping[str, ResolvedBeamModel],
) -> None:
    """Validate external beam_model override against the observation pointing beam.

    Only the beam used for observation-time pointing correction must match.  The
    other beam positions may legitimately change for later converter/sunscan
    analysis after a new multibeam sunscan solution.
    """

    pr = snapshot.get("pointing_reference_beam", {})
    if not isinstance(pr, Mapping):
        return
    ref_id = str(pr.get("pointing_reference_beam_id", "B00") or "B00")
    if not ref_id or ref_id == "B00":
        return
    snapshot_beams = snapshot.get("beams", {})
    if not isinstance(snapshot_beams, Mapping) or ref_id not in snapshot_beams:
        raise ValueError(f"snapshot pointing_reference_beam_id={ref_id!r} is not present in snapshot beams")
    if ref_id not in external_beams:
        raise ValueError(f"external beam_model is missing pointing_reference_beam_id={ref_id!r}")
    snap_sig = _beam_position_signature_for_pointing(dict(snapshot_beams[ref_id] or {}))
    ext_sig = _beam_position_signature_for_pointing(asdict(external_beams[ref_id]))
    if snap_sig != ext_sig:
        raise ValueError(
            f"external beam_model changes the observation-time pointing_reference_beam_id={ref_id!r}. "
            "This would break boresight/beam-center consistency; use the snapshot beam for that ID or "
            "rerun with an explicitly validated override workflow."
        )


def _resolved_beams_from_snapshot(snapshot: Mapping[str, Any]) -> Dict[str, ResolvedBeamModel]:
    helper = _load_necst_helper_module("beam_model")
    doc = helper.beam_model_from_snapshot(snapshot)
    beams: Dict[str, ResolvedBeamModel] = {}
    for beam_id, raw_beam in doc.beams.items():
        beams[str(beam_id)] = _build_resolved_beam(str(beam_id), raw_beam)
    return beams


def _analysis_selection_from_snapshot(
    *,
    stream_ids: Sequence[str],
    analysis_selection: Optional[Mapping[str, Any]],
) -> AnalysisStreamSelection:
    selection_helper = _load_necst_helper_module("analysis_stream_selection")
    available = [str(s) for s in stream_ids]
    convert = list(selection_helper.resolve_stream_selection(
        available,
        purpose="convert",
        analysis_selection=analysis_selection,
        default_streams=available,
    ).selected)
    extract = list(selection_helper.resolve_stream_selection(
        available,
        purpose="sunscan.extract",
        analysis_selection=analysis_selection,
        default_streams=available,
    ).selected)
    fit = list(selection_helper.resolve_stream_selection(
        available,
        purpose="sunscan.fit",
        analysis_selection=analysis_selection,
        default_streams=available,
    ).selected)

    policies: Dict[str, StreamSelectionPolicy] = {}
    convert_set = set(convert)
    extract_set = set(extract)
    fit_set = set(fit)
    for stream_id in available:
        policies[stream_id] = StreamSelectionPolicy(
            enabled=True,
            use_for_convert=stream_id in convert_set,
            use_for_sunscan=stream_id in extract_set,
            use_for_fit=stream_id in fit_set,
            beam_fit_use=None,
        )
    return AnalysisStreamSelection(
        policies=policies,
        convert_stream_ids=convert,
        sunscan_extract_stream_ids=extract,
        sunscan_fit_stream_ids=fit,
    )


def load_spectral_recording_snapshot_config_bundle(
    snapshot_path: Path | str,
    *,
    analysis_stream_selection_path: Optional[Path | str] = None,
    beam_model_path: Optional[Path | str] = None,
    allow_beam_model_override: bool = False,
    stream_ids: Optional[Sequence[str]] = None,
) -> ResolvedConfigBundle:
    """Load a new-DB spectral-recording snapshot into ``ResolvedConfigBundle``.

    This is the PR8g-compatible bridge from observation-time snapshot truth into
    the new separated converter/sunscan internal model.  TP streams are not
    materialized as normal spectra.  If an external beam model is supplied, it
    is treated as an explicit analysis-time override and requires
    ``allow_beam_model_override=True``.
    """

    srsnap = _load_necst_helper_module("spectral_recording_snapshot")
    selection_helper = _load_necst_helper_module("analysis_stream_selection")
    beam_helper = _load_necst_helper_module("beam_model")

    path = Path(snapshot_path)
    snapshot = srsnap.load_spectral_recording_snapshot(path)

    if beam_model_path is not None:
        if not allow_beam_model_override:
            raise ValueError(
                "external beam_model_path is an analysis-time override; "
                "pass allow_beam_model_override=True to use it"
            )
        beam_doc = beam_helper.load_beam_model(beam_model_path)
        external_beams = {str(beam_id): _build_resolved_beam(str(beam_id), raw_beam) for beam_id, raw_beam in beam_doc.beams.items()}
        _validate_pointing_reference_beam_override(snapshot=snapshot, external_beams=external_beams)
        snapshot = beam_helper.apply_beam_model_to_snapshot(snapshot, beam_doc)

    layouts_all = list(srsnap.extract_stream_layouts(snapshot))
    spectral_layouts = [layout for layout in layouts_all if not bool(getattr(layout, "is_tp"))]
    if stream_ids is not None:
        wanted = {str(x) for x in stream_ids}
        by_name = {str(layout.stream_id): layout for layout in spectral_layouts}
        missing = sorted(wanted - set(by_name))
        if missing:
            raise ValueError(f"requested snapshot stream_id(s) are not available as spectral streams: {missing}")
        spectral_layouts = [layout for layout in spectral_layouts if str(layout.stream_id) in wanted]
    if not spectral_layouts:
        raise ValueError("snapshot has no spectral streams available for converter/sunscan materialization")

    raw_streams = dict(snapshot.get("streams", {}) or {})
    beams = _resolved_beams_from_snapshot(snapshot)
    streams: List[ResolvedStream] = []
    for idx, layout in enumerate(spectral_layouts):
        stream_id = str(layout.stream_id)
        raw_entry = dict(raw_streams.get(stream_id, {}) or {})
        frequency_axis = _snapshot_axis_to_legacy_frequency_axis(
            stream_id=stream_id,
            stream_entry=raw_entry,
            layout=layout,
            snapshot=snapshot,
        )
        local_oscillators = _snapshot_lo_chain_to_local_oscillators(
            stream_entry=raw_entry,
            snapshot=snapshot,
        )
        # Preserve the snapshot DB binding.  ``db_stream_name`` is an alias used
        # by old configs, while ``db_table_path`` is the new-DB truth.  Do not
        # collapse ``data/spectral/xffts/board1`` to just ``board1``: doing so
        # loses the namespace-relative path needed by snapshot-aware readers.
        db_table_name = _nonempty_str(raw_entry.get("db_table_name"))
        if db_table_name is None:
            db_table_name = _nonempty_str(getattr(layout, "db_table_path", None))
        aliases = _unique_aliases(
            stream_id,
            getattr(layout, "db_stream_name", None),
            raw_entry.get("legacy_aliases"),
        )
        streams.append(
            ResolvedStream(
                stream_id=stream_id,
                db_stream_name=str(layout.db_stream_name),
                db_table_name=db_table_name,
                fdnum=int(layout.fdnum),
                ifnum=int(layout.ifnum),
                plnum=int(layout.plnum),
                polariza=str(layout.polariza).upper(),
                beam_id=str(layout.beam_id),
                stream_index=idx,
                aliases=aliases,
                spectrometer_key=str(layout.spectrometer_key),
                board_id=int(layout.board_id),
                sampler=_nonempty_str(raw_entry.get("sampler")),
                frontend=_nonempty_str(raw_entry.get("frontend")),
                backend=_nonempty_str(raw_entry.get("backend")),
                frequency_axis=frequency_axis,
                local_oscillators=local_oscillators,
                override={
                    key: copy.deepcopy(raw_entry[key])
                    for key in (
                        "source_stream_id",
                        "source_db_stream_name",
                        "recorded_stream_id",
                        "recorded_db_stream_name",
                        "window_id",
                        "line_name",
                    )
                    if key in raw_entry
                },
                channel_slice_spec=None,
                raw=copy.deepcopy(raw_entry),
            )
        )

    missing_beams = sorted({s.beam_id for s in streams} - set(beams))
    if missing_beams:
        raise ValueError(f"snapshot stream(s) reference missing beam_id(s): {missing_beams}")

    analysis_selection = selection_helper.load_analysis_stream_selection(analysis_stream_selection_path)
    selection = _analysis_selection_from_snapshot(
        stream_ids=[s.stream_id for s in streams],
        analysis_selection=analysis_selection,
    )

    provenance = copy.deepcopy(dict(snapshot.get("provenance", {}) or {}))
    provenance.update(
        {
            "config_loader": "snapshot_adapter",
            "config_source_format": "spectral_recording_snapshot",
            "config_source_path": str(path),
            "setup_id": snapshot.get("setup_id"),
            "canonical_snapshot_sha256": snapshot.get("canonical_snapshot_sha256"),
            "analysis_stream_selection_path": None if analysis_stream_selection_path is None else str(Path(analysis_stream_selection_path)),
            "beam_model_override_used": bool(beam_model_path is not None),
            "beam_model_override_path": None if beam_model_path is None else str(Path(beam_model_path)),
            "excluded_tp_stream_ids": [
                str(layout.stream_id) for layout in layouts_all if bool(getattr(layout, "is_tp"))
            ],
        }
    )

    return ResolvedConfigBundle(
        schema_version=1,
        config_name=_nonempty_str(snapshot.get("setup_id")),
        config_description=_nonempty_str(snapshot.get("snapshot_version")),
        source_format="spectral_recording_snapshot",
        source_path=str(path),
        global_config={},
        provenance=provenance,
        streams=streams,
        beams=beams,
        selection=selection,
        db_time_environment=DbTimeEnvironmentConfig(values={}),
        converter_analysis=ConverterAnalysisConfig(values={}, export_slices_by_stream_id={}),
        sunscan_analysis=SunScanAnalysisConfigModel(values={}, per_stream_overrides={}),
    )

