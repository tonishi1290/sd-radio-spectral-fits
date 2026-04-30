"""Analysis-time stream selection for NECST spectral-recording snapshots.

This module deliberately belongs to the post-observation analysis layer.  The
spectral-recording snapshot records what was observed and how it was stored;
this file reads a separate TOML that says which stream_id values should be used
for conversion, sunscan extraction, or sunscan fitting.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:  # Python >= 3.11
    import tomllib  # type: ignore[attr-defined]
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore[no-redef]


@dataclass(frozen=True)
class StreamSelectionResult:
    """Resolved stream_id selection for one analysis purpose."""

    selected: Tuple[str, ...]
    reason: str
    explicit: bool = False
    requested: Tuple[str, ...] = ()
    excluded: Tuple[str, ...] = ()


def _load_toml(path: Path) -> Dict[str, Any]:
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def load_analysis_stream_selection(path: Optional[Path | str]) -> Dict[str, Any]:
    """Load an analysis stream-selection TOML file.

    Supported canonical form::

        [convert]
        streams = ["xffts_board3", "xffts_board4"]

        [sunscan.extract]
        streams = ["xffts_board3", "xffts_board4"]

        [sunscan.fit]
        streams = ["xffts_board3"]
        exclude_streams = ["xffts_board5"]

    Backward-friendly aliases are also accepted under ``[sunscan]``:
    ``extract_streams`` and ``fit_streams``.
    """

    if path is None:
        return {}
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"analysis stream-selection file does not exist: {p}")
    data = _load_toml(p)
    if not isinstance(data, dict):
        raise ValueError(f"analysis stream-selection file must contain TOML tables: {p}")
    data = dict(data)
    data["_selection_path"] = str(p)
    return data


def _as_string_list(value: Any, *, field: str) -> Optional[List[str]]:
    if value is None:
        return None
    if isinstance(value, str):
        items = [value]
    elif isinstance(value, (list, tuple)):
        items = list(value)
    else:
        raise ValueError(f"{field} must be a string or list of strings")
    out: List[str] = []
    for item in items:
        text = str(item).strip()
        if not text:
            raise ValueError(f"{field} contains an empty stream name")
        out.append(text)
    return out


def _purpose_blocks(selection: Mapping[str, Any], purpose: str) -> List[Tuple[str, Mapping[str, Any]]]:
    purpose = str(purpose).strip().lower()
    blocks: List[Tuple[str, Mapping[str, Any]]] = []
    if purpose == "convert":
        block = selection.get("convert", {})
        if isinstance(block, Mapping):
            blocks.append(("convert", block))
        return blocks
    if purpose in ("sunscan", "sunscan.extract", "extract"):
        sun = selection.get("sunscan", {})
        if isinstance(sun, Mapping):
            block = sun.get("extract", {})
            if isinstance(block, Mapping):
                blocks.append(("sunscan.extract", block))
            # alias: [sunscan] extract_streams = [...]
            alias: Dict[str, Any] = {}
            if "extract_streams" in sun:
                alias["streams"] = sun.get("extract_streams")
            if "exclude_extract_streams" in sun:
                alias["exclude_streams"] = sun.get("exclude_extract_streams")
            if alias:
                blocks.append(("sunscan.extract alias", alias))
        return blocks
    if purpose in ("fit", "sunscan.fit"):
        sun = selection.get("sunscan", {})
        if isinstance(sun, Mapping):
            block = sun.get("fit", {})
            if isinstance(block, Mapping):
                blocks.append(("sunscan.fit", block))
            alias: Dict[str, Any] = {}
            if "fit_streams" in sun:
                alias["streams"] = sun.get("fit_streams")
            if "exclude_fit_streams" in sun:
                alias["exclude_streams"] = sun.get("exclude_fit_streams")
            if alias:
                blocks.append(("sunscan.fit alias", alias))
        return blocks
    raise ValueError(f"unsupported analysis selection purpose: {purpose!r}")


def _selection_from_file(selection: Mapping[str, Any], purpose: str) -> Tuple[Optional[List[str]], List[str], List[str]]:
    requested: Optional[List[str]] = None
    excluded: List[str] = []
    sources: List[str] = []
    for label, block in _purpose_blocks(selection, purpose):
        streams = _as_string_list(block.get("streams"), field=f"{label}.streams")
        excludes = _as_string_list(block.get("exclude_streams"), field=f"{label}.exclude_streams") or []
        if streams is not None:
            if requested is not None and streams != requested:
                raise ValueError(
                    f"analysis stream-selection has conflicting stream lists for {purpose!r}: {requested!r} vs {streams!r}"
                )
            requested = streams
            sources.append(label)
        excluded.extend(excludes)
    return requested, excluded, sources


def _dedup_keep_order(names: Iterable[str]) -> List[str]:
    out: List[str] = []
    seen = set()
    for name in names:
        text = str(name)
        if text not in seen:
            seen.add(text)
            out.append(text)
    return out


def resolve_stream_selection(
    available_streams: Sequence[str],
    *,
    purpose: str,
    explicit_stream_names: Optional[Sequence[str]] = None,
    analysis_selection: Optional[Mapping[str, Any]] = None,
    default_streams: Optional[Sequence[str]] = None,
) -> StreamSelectionResult:
    """Resolve stream_id selection.

    Precedence:

    1. explicit CLI/API stream list;
    2. analysis stream-selection TOML;
    3. caller-provided default stream list;
    4. all available streams.
    """

    available = _dedup_keep_order([str(s) for s in available_streams])
    available_set = set(available)
    default_list = _dedup_keep_order([str(s) for s in (default_streams if default_streams is not None else available)])

    if explicit_stream_names:
        requested = _dedup_keep_order([str(s) for s in explicit_stream_names])
        missing = [s for s in requested if s not in available_set]
        if missing:
            raise ValueError(f"requested stream_id(s) are not available for {purpose}: {missing}")
        return StreamSelectionResult(tuple(requested), reason="explicit-stream-name", explicit=True, requested=tuple(requested))

    requested_from_file: Optional[List[str]] = None
    excluded: List[str] = []
    sources: List[str] = []
    if analysis_selection:
        requested_from_file, excluded, sources = _selection_from_file(analysis_selection, purpose)

    if requested_from_file is not None:
        requested = _dedup_keep_order(requested_from_file)
        missing = [s for s in requested if s not in available_set]
        if missing:
            raise ValueError(f"analysis stream-selection references unavailable stream_id(s) for {purpose}: {missing}")
        selected = requested
        reason = "analysis-stream-selection:" + ",".join(sources)
    else:
        selected = list(default_list)
        reason = "default"

    if excluded:
        exclude_set = set(_dedup_keep_order(excluded))
        missing_excluded = sorted(exclude_set - available_set)
        if missing_excluded:
            raise ValueError(f"analysis stream-selection exclude_streams references unavailable stream_id(s) for {purpose}: {missing_excluded}")
        selected = [s for s in selected if s not in exclude_set]
        reason += "+exclude"
    return StreamSelectionResult(tuple(selected), reason=reason, requested=tuple(requested_from_file or ()), excluded=tuple(_dedup_keep_order(excluded)))


__all__ = [
    "StreamSelectionResult",
    "load_analysis_stream_selection",
    "resolve_stream_selection",
]
