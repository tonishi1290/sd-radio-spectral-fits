"""NECST spectral-recording snapshot reader and DB sidecar discovery helpers.

This module is the PR2a foundation for the NECST/XFFTS spectral recording
redesign.  It deliberately does *not* change the converter runtime path yet;
PR2b and later PRs will call these helpers from ``necst_v4_sdfits_converter.py``
and the sunscan/multibeam adapters.

Definitions used here are the same as the v26-aligned specification:

* ``spectral_recording_snapshot.toml`` is the new-DB truth for frequency axes,
  saved channel ranges, TP/spectrum mode, and stream identity.
* ``db_table_path`` is the canonical NECSTDB logical path from the snapshot.
* ``db_stream_name`` is a legacy-compatible stream/table name used by the older
  converter/sunscan view.
* DB sidecar TOMLs other than the snapshot are validation/display material; they
  must not silently re-resolve saved channel windows when a snapshot exists.

The code is dependency-light and safe to import in tests that do not have
Astropy, necstdb, or the sd_radio_spectral_fits package installed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import copy
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:  # Python >= 3.11
    import tomllib as _toml_reader  # type: ignore[attr-defined]
except ModuleNotFoundError:  # pragma: no cover - depends on runtime Python
    try:
        import tomli as _toml_reader  # type: ignore[no-redef]
    except ModuleNotFoundError:  # pragma: no cover - optional fallback
        try:
            import toml as _toml_legacy  # type: ignore[no-redef]
        except ModuleNotFoundError as _toml_import_error:  # pragma: no cover
            _toml_reader = None  # type: ignore[assignment]
            _TOML_IMPORT_ERROR = _toml_import_error
        else:  # pragma: no cover
            _toml_reader = _toml_legacy  # type: ignore[assignment]
            _TOML_IMPORT_ERROR = None
    else:
        _TOML_IMPORT_ERROR = None
else:
    _TOML_IMPORT_ERROR = None


SNAPSHOT_SCHEMA_VERSION = "spectral_recording_snapshot_v1"
CANONICAL_SNAPSHOT_NAME = "spectral_recording_snapshot.toml"
CANONICAL_SIDECAR_NAMES = {
    "snapshot": CANONICAL_SNAPSHOT_NAME,
    "lo_profile": "lo_profile.toml",
    "recording_window_setup": "recording_window_setup.toml",
    "beam_model": "beam_model.toml",
    "pointing_param": "pointing_param.toml",
}
CONFIG_SUFFIXES = ("config.toml", "_config.toml")
OBS_SUFFIX = ".obs"

LAYOUT_LEGACY_DB = "LEGACY_DB"
LAYOUT_INCOMPLETE_NEW_DB = "INCOMPLETE_NEW_DB"
LAYOUT_NEW_FULL_SPECTRUM = "NEW_FULL_SPECTRUM"
LAYOUT_NEW_SLICED_SPECTRUM = "NEW_SLICED_SPECTRUM"
LAYOUT_NEW_MIXED_SPECTRUM = "NEW_MIXED_SPECTRUM"
LAYOUT_NEW_TP = "NEW_TP"
LAYOUT_NEW_MIXED_WITH_TP = "NEW_MIXED_WITH_TP"


class SpectralRecordingSnapshotError(ValueError):
    """Raised for invalid snapshot or sidecar discovery state."""


class SnapshotHashMismatchError(SpectralRecordingSnapshotError):
    """Raised when ``canonical_snapshot_sha256`` does not match the contents."""


@dataclass(frozen=True)
class SidecarDiscovery:
    """Result of scanning a NECST RawData/DB directory for setup sidecars."""

    db_dir: Path
    snapshot_path: Optional[Path] = None
    lo_profile_path: Optional[Path] = None
    recording_window_setup_path: Optional[Path] = None
    beam_model_path: Optional[Path] = None
    pointing_param_path: Optional[Path] = None
    obs_paths: Tuple[Path, ...] = ()
    config_paths: Tuple[Path, ...] = ()
    extra_toml_paths: Tuple[Path, ...] = ()
    duplicate_canonical_paths: Mapping[str, Tuple[Path, ...]] = field(default_factory=dict)
    warnings: Tuple[str, ...] = ()

    def has_new_sidecars(self) -> bool:
        """Return True only for spectral-recording sidecars that imply new setup.

        ``.obs`` files and generic ``config.toml`` / ``*_config.toml`` files are
        intentionally *not* counted here.  Those files may also be present in
        historical RawData directories; treating them as incomplete new spectral
        recording DBs would incorrectly block the legacy converter fallback.
        """

        return any(
            p is not None
            for p in (
                self.snapshot_path,
                self.lo_profile_path,
                self.recording_window_setup_path,
                self.beam_model_path,
                self.pointing_param_path,
            )
        )

    def as_dict(self) -> Dict[str, Any]:
        return {
            "db_dir": str(self.db_dir),
            "snapshot_path": _path_to_str(self.snapshot_path),
            "lo_profile_path": _path_to_str(self.lo_profile_path),
            "recording_window_setup_path": _path_to_str(self.recording_window_setup_path),
            "beam_model_path": _path_to_str(self.beam_model_path),
            "pointing_param_path": _path_to_str(self.pointing_param_path),
            "obs_paths": [str(p) for p in self.obs_paths],
            "config_paths": [str(p) for p in self.config_paths],
            "extra_toml_paths": [str(p) for p in self.extra_toml_paths],
            "duplicate_canonical_paths": {
                str(k): [str(p) for p in tuple(v)]
                for k, v in dict(self.duplicate_canonical_paths or {}).items()
            },
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class StreamLayout:
    """Per-stream layout extracted from a spectral-recording snapshot."""

    stream_id: str
    recording_mode: str
    recording_table_kind: str
    spectrometer_key: str
    board_id: int
    db_table_path: str
    db_stream_name: str
    fdnum: int
    ifnum: int
    plnum: int
    polariza: str
    beam_id: str
    full_nchan: int
    input_data_nchan: int
    saved_ch_start: int
    saved_ch_stop: int
    saved_nchan: int
    # Backward-compatible legacy fields may exist in older snapshots, but PR8g
    # analysis selection is stream_id based and does not treat use_for_* as
    # observation-time truth.
    use_for_convert: bool = True
    use_for_sunscan: bool = True
    use_for_fit: bool = True
    enabled: bool = True
    analysis_defaults: Mapping[str, Any] = field(default_factory=dict)

    @property
    def is_tp(self) -> bool:
        return self.recording_mode == "tp" or self.recording_table_kind == "tp"

    @property
    def is_full_spectrum(self) -> bool:
        return (
            not self.is_tp
            and self.saved_ch_start == 0
            and self.saved_ch_stop == self.full_nchan
            and self.saved_nchan == self.full_nchan
        )

    @property
    def is_sliced_spectrum(self) -> bool:
        return (not self.is_tp) and (not self.is_full_spectrum)

    def local_channel_to_full_channel(self, c_local: int) -> int:
        if c_local < 0 or c_local >= self.saved_nchan:
            raise SpectralRecordingSnapshotError(
                f"local channel {c_local} is outside [0,{self.saved_nchan}) for stream {self.stream_id!r}"
            )
        return int(self.saved_ch_start) + int(c_local)


@dataclass(frozen=True)
class LayoutSummary:
    """Overall new/legacy DB layout classification."""

    layout_kind: str
    stream_layouts: Tuple[StreamLayout, ...] = ()
    snapshot_path: Optional[Path] = None
    warnings: Tuple[str, ...] = ()

    @property
    def has_snapshot(self) -> bool:
        return self.snapshot_path is not None

    @property
    def has_tp(self) -> bool:
        return any(s.is_tp for s in self.stream_layouts)

    @property
    def has_sliced_spectrum(self) -> bool:
        return any(s.is_sliced_spectrum for s in self.stream_layouts)

    @property
    def has_full_spectrum(self) -> bool:
        return any(s.is_full_spectrum for s in self.stream_layouts)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "layout_kind": self.layout_kind,
            "snapshot_path": _path_to_str(self.snapshot_path),
            "warnings": list(self.warnings),
            "stream_layouts": [s.__dict__.copy() for s in self.stream_layouts],
        }


def _path_to_str(path: Optional[Path]) -> Optional[str]:
    return None if path is None else str(path)


def _toml_loads(text: str) -> Dict[str, Any]:
    if _toml_reader is None:  # pragma: no cover
        raise SpectralRecordingSnapshotError(
            "No TOML reader is available. Use Python >= 3.11, or install tomli/toml."
        ) from _TOML_IMPORT_ERROR
    if hasattr(_toml_reader, "loads"):
        return _toml_reader.loads(text)  # type: ignore[no-any-return]
    raise SpectralRecordingSnapshotError("Configured TOML reader has no loads() function")


def read_toml_file(path: os.PathLike[str] | str) -> Dict[str, Any]:
    p = Path(path)
    try:
        return _toml_loads(p.read_text(encoding="utf-8"))
    except OSError as exc:
        raise SpectralRecordingSnapshotError(f"Cannot read TOML file: {p}") from exc
    except Exception as exc:
        if isinstance(exc, SpectralRecordingSnapshotError):
            raise
        raise SpectralRecordingSnapshotError(f"Cannot parse TOML file: {p}: {exc}") from exc


def sha256_file(path: os.PathLike[str] | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _format_float(value: float) -> str:
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return format(float(value), ".17g")


def _toml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return _format_float(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if value is None:
        return json.dumps("")
    if isinstance(value, list):
        return "[" + ", ".join(_toml_scalar(v) for v in value) + "]"
    raise TypeError(f"Unsupported TOML scalar value: {value!r}")


def dumps_toml(data: Mapping[str, Any]) -> str:
    """Return the deterministic TOML subset used by the NECST PR1 resolver."""

    lines: List[str] = []

    def emit_table(path: List[str], table: Mapping[str, Any]) -> None:
        scalar_keys: List[str] = []
        table_keys: List[str] = []
        array_table_keys: List[str] = []
        for key in sorted(table):
            value = table[key]
            if isinstance(value, Mapping):
                table_keys.append(key)
            elif isinstance(value, list) and all(isinstance(x, Mapping) for x in value):
                array_table_keys.append(key)
            else:
                scalar_keys.append(key)
        if path:
            if lines and lines[-1] != "":
                lines.append("")
            lines.append("[" + ".".join(path) + "]")
        for key in scalar_keys:
            lines.append(f"{key} = {_toml_scalar(table[key])}")
        for key in table_keys:
            emit_table(path + [key], table[key])  # type: ignore[arg-type]
        for key in array_table_keys:
            for item in table[key]:
                if lines and lines[-1] != "":
                    lines.append("")
                lines.append("[[" + ".".join(path + [key]) + "]]")
                emit_array_item(path + [key], item)

    def emit_array_item(path: List[str], table: Mapping[str, Any]) -> None:
        scalar_keys: List[str] = []
        table_keys: List[str] = []
        for key in sorted(table):
            value = table[key]
            if isinstance(value, Mapping):
                table_keys.append(key)
            else:
                scalar_keys.append(key)
        for key in scalar_keys:
            lines.append(f"{key} = {_toml_scalar(table[key])}")
        for key in table_keys:
            emit_table(path + [key], table[key])  # type: ignore[arg-type]

    emit_table([], data)
    return "\n".join(lines).rstrip() + "\n"


def canonical_snapshot_sha256(snapshot: Mapping[str, Any]) -> str:
    """Recompute ``canonical_snapshot_sha256`` from a loaded snapshot.

    The hash is defined over the deterministic TOML representation with the
    ``canonical_snapshot_sha256`` field itself blanked, matching PR1.
    """

    tmp = copy.deepcopy(dict(snapshot))
    tmp["canonical_snapshot_sha256"] = ""
    return sha256_text(dumps_toml(tmp))


def _utf8_len(value: Any) -> int:
    return len(str(value).encode("utf-8"))


def _require_fixed_string_limit(value: Any, limit: int, *, field: str) -> None:
    length = _utf8_len(value)
    if length > limit:
        raise SpectralRecordingSnapshotError(
            f"{field} is too long for fixed schema: {length} bytes > {limit} bytes"
        )


def _validate_fixed_schema_metadata(snapshot: Mapping[str, Any]) -> None:
    _require_fixed_string_limit(str(snapshot.get("setup_id", "")), 64, field="snapshot.setup_id")
    setup_hash = str(snapshot.get("canonical_snapshot_sha256", "") or snapshot.get("setup_hash", ""))
    if setup_hash:
        _require_fixed_string_limit(setup_hash, 64, field="snapshot.canonical_snapshot_sha256")
    streams = snapshot.get("streams")
    if not isinstance(streams, Mapping):
        return
    for table_key, stream in streams.items():
        if not isinstance(stream, Mapping):
            continue
        _require_fixed_string_limit(str(stream.get("stream_id", table_key)), 64, field=f"streams.{table_key}.stream_id")
        _require_fixed_string_limit(str(stream.get("spectrometer_key", "")), 32, field=f"streams.{table_key}.spectrometer_key")
        _require_fixed_string_limit(str(stream.get("polariza", "")), 8, field=f"streams.{table_key}.polariza")
        _require_fixed_string_limit(str(stream.get("beam_id", "")), 32, field=f"streams.{table_key}.beam_id")
        _require_fixed_string_limit(str(stream.get("saved_window_policy", "")), 32, field=f"streams.{table_key}.saved_window_policy")
        if str(stream.get("recording_mode", "spectrum")) == "tp":
            _require_fixed_string_limit(str(stream.get("tp_stat", "sum_mean")), 16, field=f"streams.{table_key}.tp_stat")
            windows_json = str(stream.get("computed_windows_json", "") or "")
            if not windows_json:
                start = int(stream.get("saved_ch_start", 0))
                stop = int(stream.get("saved_ch_stop", 0))
                windows_json = json.dumps(
                    [{"kind": "channel", "computed_ch_start": start, "computed_ch_stop": stop}],
                    sort_keys=True,
                    separators=(",", ":"),
                )
            _require_fixed_string_limit(windows_json, 2048, field=f"streams.{table_key}.tp_windows_json")


def _validate_channel_order_metadata(snapshot: Mapping[str, Any]) -> None:
    axes = snapshot.get("frequency_axes", {})
    if not isinstance(axes, Mapping):
        return
    increasing = {"increasing_if", "increasing_if_frequency", "increasing"}
    decreasing = {"decreasing_if", "decreasing_if_frequency", "decreasing"}
    for axis_id, axis in axes.items():
        if not isinstance(axis, Mapping):
            continue
        order = str(axis.get("channel_order", "increasing_if"))
        if "if_freq_step_hz" not in axis:
            continue
        step = float(axis.get("if_freq_step_hz"))
        if order in increasing and step < 0:
            raise SpectralRecordingSnapshotError(
                f"frequency_axes.{axis_id}.channel_order={order!r} requires if_freq_step_hz > 0, got {step}"
            )
        if order in decreasing and step > 0:
            raise SpectralRecordingSnapshotError(
                f"frequency_axes.{axis_id}.channel_order={order!r} requires if_freq_step_hz < 0, got {step}"
            )


def validate_snapshot(snapshot: Mapping[str, Any], *, verify_hash: bool = True) -> None:
    """Validate the snapshot fields required by PR2a and downstream adapters."""

    if snapshot.get("schema_version") != SNAPSHOT_SCHEMA_VERSION:
        raise SpectralRecordingSnapshotError(
            f"Unsupported schema_version: {snapshot.get('schema_version')!r}"
        )
    _validate_channel_order_metadata(snapshot)
    _validate_fixed_schema_metadata(snapshot)

    streams = snapshot.get("streams")
    if not isinstance(streams, Mapping) or not streams:
        raise SpectralRecordingSnapshotError("snapshot.streams must be a non-empty table")
    seen_stream: set[str] = set()
    seen_raw: Dict[Tuple[str, int], str] = {}
    seen_db: Dict[str, str] = {}
    seen_sdfits: Dict[Tuple[int, int, int, str], str] = {}
    for table_key, raw_stream in streams.items():
        if not isinstance(raw_stream, Mapping):
            raise SpectralRecordingSnapshotError(f"snapshot.streams.{table_key} must be a table")
        layout = stream_layout_from_snapshot_entry(str(table_key), raw_stream)
        if layout.stream_id != str(table_key):
            raise SpectralRecordingSnapshotError(
                f"snapshot stream table key {table_key!r} does not match stream_id={layout.stream_id!r}"
            )
        if layout.stream_id in seen_stream:
            raise SpectralRecordingSnapshotError(f"Duplicate stream_id: {layout.stream_id}")
        seen_stream.add(layout.stream_id)
        raw_key = (layout.spectrometer_key, layout.board_id)
        is_recorded_product = bool(raw_stream.get("source_stream_id")) or bool(raw_stream.get("recorded_stream_id"))
        if raw_key in seen_raw and not is_recorded_product:
            raise SpectralRecordingSnapshotError(
                f"Duplicate raw input key {raw_key}: {seen_raw[raw_key]} and {layout.stream_id}"
            )
        seen_raw.setdefault(raw_key, layout.stream_id)
        if layout.db_table_path in seen_db:
            raise SpectralRecordingSnapshotError(
                f"Duplicate db_table_path {layout.db_table_path!r}: {seen_db[layout.db_table_path]} and {layout.stream_id}"
            )
        seen_db[layout.db_table_path] = layout.stream_id
        sdfits_key = (layout.fdnum, layout.ifnum, layout.plnum, layout.polariza)
        if sdfits_key in seen_sdfits and not is_recorded_product:
            raise SpectralRecordingSnapshotError(
                f"Duplicate SDFITS key {sdfits_key}: {seen_sdfits[sdfits_key]} and {layout.stream_id}"
            )
        seen_sdfits.setdefault(sdfits_key, layout.stream_id)

    if verify_hash:
        expected = str(snapshot.get("canonical_snapshot_sha256", "")).strip()
        if not expected:
            raise SnapshotHashMismatchError("snapshot has no canonical_snapshot_sha256")
        actual = canonical_snapshot_sha256(snapshot)
        if actual != expected:
            raise SnapshotHashMismatchError(
                f"canonical_snapshot_sha256 mismatch: expected {expected}, recomputed {actual}"
            )


def load_spectral_recording_snapshot(
    path: os.PathLike[str] | str,
    *,
    verify_hash: bool = True,
) -> Dict[str, Any]:
    """Load and validate ``spectral_recording_snapshot.toml``."""

    snapshot = read_toml_file(path)
    validate_snapshot(snapshot, verify_hash=verify_hash)
    return snapshot


def _as_bool(value: Any, *, default: bool = True) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return bool(value)
    s = str(value).strip().lower()
    if s in ("1", "true", "yes", "on"):
        return True
    if s in ("0", "false", "no", "off"):
        return False
    raise SpectralRecordingSnapshotError(f"Expected bool-compatible value, got {value!r}")


def _as_int(value: Any, *, field: str) -> int:
    if isinstance(value, bool):
        raise SpectralRecordingSnapshotError(f"{field} must be int, got bool")
    try:
        return int(value)
    except Exception as exc:
        raise SpectralRecordingSnapshotError(f"{field} must be int, got {value!r}") from exc


def _as_str(value: Any, *, field: str) -> str:
    if value is None:
        raise SpectralRecordingSnapshotError(f"{field} must be string, got None")
    s = str(value).strip()
    if not s:
        raise SpectralRecordingSnapshotError(f"{field} must be non-empty string")
    return s


def stream_layout_from_snapshot_entry(table_key: str, stream: Mapping[str, Any]) -> StreamLayout:
    """Convert one ``snapshot.streams.<stream_id>`` table into ``StreamLayout``."""

    stream_id = _as_str(stream.get("stream_id", table_key), field=f"streams.{table_key}.stream_id")
    recording_mode = _as_str(stream.get("recording_mode"), field=f"streams.{table_key}.recording_mode")
    recording_table_kind = _as_str(
        stream.get("recording_table_kind", "tp" if recording_mode == "tp" else "spectral"),
        field=f"streams.{table_key}.recording_table_kind",
    )
    if recording_mode not in ("spectrum", "tp"):
        raise SpectralRecordingSnapshotError(
            f"streams.{table_key}.recording_mode must be 'spectrum' or 'tp', got {recording_mode!r}"
        )
    if recording_table_kind not in ("spectral", "tp"):
        raise SpectralRecordingSnapshotError(
            f"streams.{table_key}.recording_table_kind must be 'spectral' or 'tp', got {recording_table_kind!r}"
        )
    full_nchan = _as_int(stream.get("full_nchan"), field=f"streams.{table_key}.full_nchan")
    saved_ch_start = _as_int(stream.get("saved_ch_start"), field=f"streams.{table_key}.saved_ch_start")
    saved_ch_stop = _as_int(stream.get("saved_ch_stop"), field=f"streams.{table_key}.saved_ch_stop")
    saved_nchan = _as_int(stream.get("saved_nchan"), field=f"streams.{table_key}.saved_nchan")
    if full_nchan <= 0:
        raise SpectralRecordingSnapshotError(f"streams.{table_key}.full_nchan must be positive")
    if not (0 <= saved_ch_start < saved_ch_stop <= full_nchan):
        raise SpectralRecordingSnapshotError(
            f"streams.{table_key}: invalid saved range [{saved_ch_start}:{saved_ch_stop}] for full_nchan={full_nchan}"
        )
    if saved_nchan != saved_ch_stop - saved_ch_start:
        raise SpectralRecordingSnapshotError(
            f"streams.{table_key}: saved_nchan={saved_nchan} does not equal saved_ch_stop-saved_ch_start"
        )
    if recording_mode == "tp" and recording_table_kind != "tp":
        raise SpectralRecordingSnapshotError(
            f"streams.{table_key}: TP stream must have recording_table_kind='tp'"
        )
    if recording_table_kind == "tp" and recording_mode != "tp":
        raise SpectralRecordingSnapshotError(
            f"streams.{table_key}: recording_table_kind='tp' requires recording_mode='tp'"
        )
    db_table_path = _as_str(stream.get("db_table_path"), field=f"streams.{table_key}.db_table_path")
    if recording_table_kind == "tp" and not db_table_path.startswith("data/tp/"):
        raise SpectralRecordingSnapshotError(
            f"streams.{table_key}: recording_table_kind='tp' requires db_table_path under 'data/tp/', got {db_table_path!r}"
        )
    if recording_table_kind == "spectral" and not db_table_path.startswith("data/spectral/"):
        raise SpectralRecordingSnapshotError(
            f"streams.{table_key}: recording_table_kind='spectral' requires db_table_path under 'data/spectral/', got {db_table_path!r}"
        )

    return StreamLayout(
        stream_id=stream_id,
        recording_mode=recording_mode,
        recording_table_kind=recording_table_kind,
        spectrometer_key=_as_str(stream.get("spectrometer_key"), field=f"streams.{table_key}.spectrometer_key"),
        board_id=_as_int(stream.get("board_id"), field=f"streams.{table_key}.board_id"),
        db_table_path=db_table_path,
        db_stream_name=_as_str(stream.get("db_stream_name"), field=f"streams.{table_key}.db_stream_name"),
        fdnum=_as_int(stream.get("fdnum"), field=f"streams.{table_key}.fdnum"),
        ifnum=_as_int(stream.get("ifnum"), field=f"streams.{table_key}.ifnum"),
        plnum=_as_int(stream.get("plnum"), field=f"streams.{table_key}.plnum"),
        polariza=_as_str(stream.get("polariza"), field=f"streams.{table_key}.polariza"),
        beam_id=_as_str(stream.get("beam_id"), field=f"streams.{table_key}.beam_id"),
        full_nchan=full_nchan,
        input_data_nchan=_as_int(stream.get("input_data_nchan", full_nchan), field=f"streams.{table_key}.input_data_nchan"),
        saved_ch_start=saved_ch_start,
        saved_ch_stop=saved_ch_stop,
        saved_nchan=saved_nchan,
        use_for_convert=_as_bool(stream.get("use_for_convert"), default=True),
        use_for_sunscan=_as_bool(stream.get("use_for_sunscan"), default=True),
        use_for_fit=_as_bool(stream.get("use_for_fit"), default=True),
        enabled=_as_bool(stream.get("enabled"), default=True),
        analysis_defaults=dict(stream.get("analysis_defaults", {}) or {}),
    )


def extract_stream_layouts(snapshot: Mapping[str, Any]) -> Tuple[StreamLayout, ...]:
    streams = snapshot.get("streams")
    if not isinstance(streams, Mapping):
        raise SpectralRecordingSnapshotError("snapshot.streams must be a table")
    return tuple(
        stream_layout_from_snapshot_entry(str(stream_id), streams[stream_id])
        for stream_id in sorted(streams)
    )


def classify_snapshot_layout(snapshot: Mapping[str, Any]) -> LayoutSummary:
    """Classify a loaded snapshot without looking at DB files."""

    validate_snapshot(snapshot, verify_hash=False)
    layouts = extract_stream_layouts(snapshot)
    if not layouts:
        raise SpectralRecordingSnapshotError("snapshot has no stream layouts")
    if all(s.is_tp for s in layouts):
        kind = LAYOUT_NEW_TP
    elif any(s.is_tp for s in layouts):
        kind = LAYOUT_NEW_MIXED_WITH_TP
    elif all(s.is_full_spectrum for s in layouts):
        kind = LAYOUT_NEW_FULL_SPECTRUM
    elif all(s.is_sliced_spectrum for s in layouts):
        kind = LAYOUT_NEW_SLICED_SPECTRUM
    else:
        kind = LAYOUT_NEW_MIXED_SPECTRUM
    return LayoutSummary(layout_kind=kind, stream_layouts=layouts)


def select_stream_layouts(
    snapshot: Mapping[str, Any],
    *,
    purpose: str = "convert",
    stream_names: Optional[Sequence[str]] = None,
) -> Tuple[StreamLayout, ...]:
    """Return stream layouts for ``convert``, ``sunscan`` or ``fit``.

    PR8g separates observation-time stream registry from analysis-time stream
    selection.  Therefore this helper no longer filters by ``use_for_*`` flags.
    Explicit ``stream_names`` selects stream_id values; otherwise all non-TP
    spectral streams are returned for convert/sunscan/fit.  ``purpose='all'``
    returns every stream including TP.
    """

    purpose = str(purpose).strip().lower()
    if purpose not in ("convert", "sunscan", "fit", "all"):
        raise SpectralRecordingSnapshotError(
            f"purpose must be convert, sunscan, fit, or all; got {purpose!r}"
        )
    requested = None if stream_names is None else {str(x) for x in stream_names}
    layouts = extract_stream_layouts(snapshot)
    by_name = {s.stream_id: s for s in layouts}
    if requested is not None:
        missing = sorted(requested - set(by_name))
        if missing:
            raise SpectralRecordingSnapshotError(f"requested unknown stream_name(s): {missing}")
    out: List[StreamLayout] = []
    for layout in layouts:
        if requested is not None and layout.stream_id not in requested:
            continue
        if purpose != "all" and layout.is_tp:
            continue
        out.append(layout)
    if requested is not None and purpose != "all":
        unsupported = sorted(name for name in requested if name in by_name and by_name[name].is_tp)
        if unsupported:
            raise SpectralRecordingSnapshotError(
                f"requested TP stream(s) are not valid for purpose {purpose!r}: {unsupported}"
            )
    return tuple(out)


def _choose_unique_named_file(
    paths: Sequence[Path],
    *,
    kind: str,
    warnings: List[str],
    strict_duplicates: bool = True,
) -> Optional[Path]:
    if not paths:
        return None
    ordered = sorted(paths, key=lambda p: (len(p.parts), str(p)))
    if len(ordered) > 1:
        candidates = ", ".join(str(p) for p in ordered)
        message = (
            f"multiple {kind} sidecars found; refusing to choose automatically. "
            f"candidate_paths=[{candidates}]. "
            "recommended_action=keep exactly one canonical sidecar in the RawData tree "
            "or pass the intended file explicitly via the corresponding CLI option."
        )
        if strict_duplicates:
            raise SpectralRecordingSnapshotError(message)
        warnings.append(message + f" using_first_for_inspection_only={ordered[0]}")
    return ordered[0]


def discover_db_sidecar_files(
    db_dir: os.PathLike[str] | str,
    *,
    recursive: bool = True,
    strict_duplicates: bool = True,
) -> SidecarDiscovery:
    """Discover snapshot and setup sidecar files under a RawData/DB directory.

    The scanner intentionally accepts both flat DB directories and DBs where
    sidecars were copied under a metadata/config subdirectory.  Exact canonical
    names are preferred. By default, duplicate canonical sidecars are a
    strict error because one DB directory must have one authoritative snapshot
    and one unambiguous set of setup sidecars.
    """

    root = Path(db_dir).expanduser().resolve()
    if not root.exists() or not root.is_dir():
        raise SpectralRecordingSnapshotError(f"DB directory does not exist or is not a directory: {root}")
    iterator = root.rglob("*") if recursive else root.glob("*")
    files = [p for p in iterator if p.is_file()]
    by_name: Dict[str, List[Path]] = {}
    for path in files:
        by_name.setdefault(path.name, []).append(path)

    duplicate_canonical_paths: Dict[str, Tuple[Path, ...]] = {}
    for label, canonical_name in CANONICAL_SIDECAR_NAMES.items():
        ordered = tuple(sorted(by_name.get(canonical_name, []), key=lambda p: (len(p.parts), str(p))))
        if len(ordered) > 1:
            duplicate_canonical_paths[label] = ordered

    warnings: List[str] = []
    snapshot_path = _choose_unique_named_file(
        by_name.get(CANONICAL_SNAPSHOT_NAME, []),
        kind="snapshot",
        warnings=warnings,
        strict_duplicates=strict_duplicates,
    )
    lo_profile_path = _choose_unique_named_file(
        by_name.get("lo_profile.toml", []),
        kind="lo_profile",
        warnings=warnings,
        strict_duplicates=strict_duplicates,
    )
    recording_window_setup_path = _choose_unique_named_file(
        by_name.get("recording_window_setup.toml", []),
        kind="recording_window_setup",
        warnings=warnings,
        strict_duplicates=strict_duplicates,
    )
    beam_model_path = _choose_unique_named_file(
        by_name.get("beam_model.toml", []),
        kind="beam_model",
        warnings=warnings,
        strict_duplicates=strict_duplicates,
    )
    pointing_param_path = _choose_unique_named_file(
        by_name.get("pointing_param.toml", []),
        kind="pointing_param",
        warnings=warnings,
        strict_duplicates=strict_duplicates,
    )

    obs_paths = tuple(sorted((p for p in files if p.name.endswith(OBS_SUFFIX)), key=str))
    config_paths = tuple(sorted((p for p in files if p.name.endswith(CONFIG_SUFFIXES)), key=str))
    recognized = {
        p
        for p in (
            snapshot_path,
            lo_profile_path,
            recording_window_setup_path,
            beam_model_path,
            pointing_param_path,
        )
        if p is not None
    }
    recognized.update(obs_paths)
    recognized.update(config_paths)
    extra_toml_paths = tuple(sorted((p for p in files if p.suffix == ".toml" and p not in recognized), key=str))

    if snapshot_path is None and any((lo_profile_path, recording_window_setup_path, beam_model_path)):
        warnings.append(
            "setup sidecar files exist but spectral_recording_snapshot.toml is absent; treat as incomplete new DB unless legacy mode is explicit"
        )

    return SidecarDiscovery(
        db_dir=root,
        snapshot_path=snapshot_path,
        lo_profile_path=lo_profile_path,
        recording_window_setup_path=recording_window_setup_path,
        beam_model_path=beam_model_path,
        pointing_param_path=pointing_param_path,
        obs_paths=obs_paths,
        config_paths=config_paths,
        extra_toml_paths=extra_toml_paths,
        duplicate_canonical_paths=duplicate_canonical_paths,
        warnings=tuple(warnings),
    )


def inspect_db_sidecars(
    db_dir: os.PathLike[str] | str,
    *,
    recursive: bool = True,
) -> Dict[str, Any]:
    """Return a JSON-serializable dry-run report for RawData sidecar discovery.

    This is intentionally non-mutating and uses ``strict_duplicates=False`` so
    that a duplicate snapshot state can be reported with all candidate paths
    instead of aborting before the operator sees what must be fixed.
    """

    sidecars = discover_db_sidecar_files(db_dir, recursive=recursive, strict_duplicates=False)
    payload = sidecars.as_dict()
    duplicates = dict(payload.get("duplicate_canonical_paths", {}) or {})
    has_snapshot = bool(payload.get("snapshot_path"))
    has_new_sidecars = bool(sidecars.has_new_sidecars())
    if duplicates:
        status = "error"
        layout_kind = "DUPLICATE_CANONICAL_SIDECARS"
        recommended_action = (
            "Do not let converter/sunscan choose automatically. Keep exactly one "
            "spectral_recording_snapshot.toml under the RawData tree, or pass the "
            "intended snapshot explicitly with --spectral-recording-snapshot."
        )
    elif has_snapshot:
        status = "ok"
        layout_kind = "NEW_DB_WITH_SNAPSHOT"
        recommended_action = (
            "Default converter/sunscan invocation may use this snapshot automatically "
            "when no explicit --spectrometer-config or --spectral-recording-snapshot is supplied."
        )
    elif has_new_sidecars:
        status = "error"
        layout_kind = LAYOUT_INCOMPLETE_NEW_DB
        recommended_action = (
            "RawData contains new spectral-recording setup sidecars but no authoritative "
            "spectral_recording_snapshot.toml. Regenerate or copy the snapshot, or pass a "
            "legacy --spectrometer-config explicitly if this is intentionally legacy data."
        )
    else:
        status = "ok"
        layout_kind = LAYOUT_LEGACY_DB
        recommended_action = (
            "No spectral_recording_snapshot.toml was found. Default converter legacy "
            "fallback is allowed; sunscan still needs an explicit config unless using "
            "a command path that supports legacy defaults."
        )
    payload.update(
        {
            "status": status,
            "layout_kind": layout_kind,
            "has_snapshot": has_snapshot,
            "has_new_sidecars": has_new_sidecars,
            "recommended_action": recommended_action,
        }
    )
    return payload


def classify_db_layout(
    db_dir: os.PathLike[str] | str,
    *,
    verify_hash: bool = True,
    recursive: bool = True,
    strict_duplicates: bool = True,
) -> LayoutSummary:
    """Classify a DB directory as legacy, incomplete-new, or new snapshot layout."""

    sidecars = discover_db_sidecar_files(db_dir, recursive=recursive, strict_duplicates=strict_duplicates)
    warnings = list(sidecars.warnings)
    if sidecars.snapshot_path is None:
        kind = LAYOUT_INCOMPLETE_NEW_DB if sidecars.has_new_sidecars() else LAYOUT_LEGACY_DB
        return LayoutSummary(layout_kind=kind, snapshot_path=None, warnings=tuple(warnings))
    snapshot = load_spectral_recording_snapshot(sidecars.snapshot_path, verify_hash=verify_hash)
    summary = classify_snapshot_layout(snapshot)
    return LayoutSummary(
        layout_kind=summary.layout_kind,
        stream_layouts=summary.stream_layouts,
        snapshot_path=sidecars.snapshot_path,
        warnings=tuple(warnings),
    )


def infer_expected_data_nchan_for_stream(layout: StreamLayout) -> int:
    """Expected stored array length for the normal spectrum converter.

    For TP streams, the normal spectrum converter must not try to read a
    spectrum; PR2b will detect ``NEW_TP`` / ``NEW_MIXED_WITH_TP`` and stop or
    route to a TP reader.  This helper therefore raises for TP.
    """

    if layout.is_tp:
        raise SpectralRecordingSnapshotError(
            f"stream {layout.stream_id!r} is TP; normal spectrum converter must not infer spectral nchan"
        )
    return int(layout.saved_nchan)


def saved_local_channel_to_full_channel(layout: StreamLayout, c_local: int) -> int:
    """Map local saved channel index to original full-channel index."""

    return layout.local_channel_to_full_channel(c_local)


def if_frequency_at_full_channel(axis: Mapping[str, Any], c_full: int) -> float:
    """IF frequency at zero-based full channel index."""

    return float(axis["if_freq_at_full_ch0_hz"]) + int(c_full) * float(axis["if_freq_step_hz"])


def if_frequency_at_saved_local_channel(
    axis: Mapping[str, Any],
    layout: StreamLayout,
    c_local: int,
) -> float:
    """IF frequency at zero-based local saved channel index."""

    return if_frequency_at_full_channel(axis, layout.local_channel_to_full_channel(c_local))


__all__ = [
    "CANONICAL_SNAPSHOT_NAME",
    "CANONICAL_SIDECAR_NAMES",
    "LAYOUT_INCOMPLETE_NEW_DB",
    "LAYOUT_LEGACY_DB",
    "LAYOUT_NEW_FULL_SPECTRUM",
    "LAYOUT_NEW_MIXED_SPECTRUM",
    "LAYOUT_NEW_MIXED_WITH_TP",
    "LAYOUT_NEW_SLICED_SPECTRUM",
    "LAYOUT_NEW_TP",
    "LayoutSummary",
    "SidecarDiscovery",
    "SnapshotHashMismatchError",
    "SpectralRecordingSnapshotError",
    "StreamLayout",
    "canonical_snapshot_sha256",
    "classify_db_layout",
    "classify_snapshot_layout",
    "discover_db_sidecar_files",
    "dumps_toml",
    "extract_stream_layouts",
    "if_frequency_at_full_channel",
    "if_frequency_at_saved_local_channel",
    "inspect_db_sidecars",
    "infer_expected_data_nchan_for_stream",
    "load_spectral_recording_snapshot",
    "read_toml_file",
    "saved_local_channel_to_full_channel",
    "select_stream_layouts",
    "sha256_file",
    "sha256_text",
    "stream_layout_from_snapshot_entry",
    "validate_snapshot",
]
