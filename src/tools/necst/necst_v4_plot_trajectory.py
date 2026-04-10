#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
necst_v4_plot_trajectory.py

Plot telescope trajectory from NECST v4 RawData without materializing spectral arrays.

Design goals
------------
- Reuse converter semantics for time normalization, table naming, site resolution,
  boresight correction, beam offset, and Az/El -> RA/Dec conversion.
- Read only the timing / OBSMODE side of the spectral table when possible.
- Use one representative stream for timing if multiple streams exist.
- Provide both a Python API and a CLI.

The script defaults to encoder-based beam-center coordinates. Optionally it can
overlay the command trajectory with partial transparency.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import importlib.util
import inspect
import pathlib
import sys
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# Headless-safe default. When run as a script, use Agg only for save-only mode
# inferred from the raw argv. This keeps the default no-option behavior interactive.
import matplotlib
if __name__ == "__main__":
    _argv0 = [str(a) for a in sys.argv[1:]]
    _has_show0 = any((a == "--show") or a.startswith("--show=") for a in _argv0)
    _has_out0 = any((a == "--out") or a.startswith("--out=") for a in _argv0)
    if (not _has_show0) and _has_out0:
        matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D



MODE_COLORS = {
    "ON": "red",
    "OFF": "blue",
    "HOT": "green",
    "OTHER": "gray",
}


@dataclass
class SimpleSite:
    lat_deg: float
    lon_deg: float
    elev_m: float

    def to_earthlocation(self):
        from astropy.coordinates import EarthLocation
        from astropy import units as u
        return EarthLocation(
            lat=float(self.lat_deg) * u.deg,
            lon=float(self.lon_deg) * u.deg,
            height=float(self.elev_m) * u.m,
        )


@dataclass
class TrajectorySamples:
    t_spec: np.ndarray
    mode_str: np.ndarray
    stream_name: str
    spec_table_name: str
    spec_time_field: str
    spec_time_basis: str
    spec_time_suffix: Optional[str]
    spec_time_fallback_field: Optional[str]
    spec_time_example: Optional[str]
    pointing: Dict[str, np.ndarray]
    coord_source: str
    coord_meaning: str
    x: np.ndarray
    y: np.ndarray
    x_label: str
    y_label: str
    frame_name: str
    coord_meta: Dict[str, Any]
    site: Optional[SimpleSite] = None
    cmd_x: Optional[np.ndarray] = None
    cmd_y: Optional[np.ndarray] = None
    cmd_source: Optional[str] = None
    cmd_meaning: Optional[str] = None
    off_blocks: Optional[List[Tuple[int, int]]] = None
    scan_id: Optional[np.ndarray] = None
    line_index: Optional[np.ndarray] = None
    line_label: Optional[np.ndarray] = None
    scan_field_name: Optional[str] = None
    line_index_field_name: Optional[str] = None
    line_label_field_name: Optional[str] = None


# -----------------------------------------------------------------------------
# Converter module loading
# -----------------------------------------------------------------------------
def _load_module_from_path(path: pathlib.Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module spec from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_converter_module(converter_path: Optional[str] = None):
    if converter_path is None:
        converter_path = pathlib.Path(__file__).with_name("necst_v4_sdfits_converter.py")
    else:
        converter_path = pathlib.Path(converter_path).expanduser().resolve()
    if not converter_path.exists():
        raise FileNotFoundError(f"converter module not found: {converter_path}")
    return _load_module_from_path(converter_path, "_necst_v4_sdfits_converter_for_plot")


# -----------------------------------------------------------------------------
# Lightweight spectral-timing extraction
# -----------------------------------------------------------------------------
def _pick_supported_kwargs(func, **kwargs):
    try:
        sig = inspect.signature(func)
    except Exception:
        return kwargs
    params = sig.parameters
    out = {}
    for key, value in kwargs.items():
        if key in params:
            out[key] = value
    return out


def _try_read_structured_subset(table, columns: Sequence[str]):
    """
    Best-effort subset read. Different necstdb versions may support different
    keyword names; we try a small set and fall back to None.
    """
    read = getattr(table, "read", None)
    if read is None:
        return None

    candidates = [
        {"astype": "array", "columns": list(columns)},
        {"astype": "array", "fields": list(columns)},
        {"astype": "array", "include_columns": list(columns)},
        {"columns": list(columns), "astype": "array"},
        {"fields": list(columns), "astype": "array"},
    ]
    for cand in candidates:
        try:
            used = _pick_supported_kwargs(read, **cand)
            if not used:
                continue
            arr = read(**used)
            if isinstance(arr, np.ndarray) and arr.dtype.names is not None:
                names = set(arr.dtype.names)
                if set(columns).issubset(names):
                    return arr
        except TypeError:
            continue
        except Exception:
            continue

    read_column = getattr(table, "read_column", None)
    if callable(read_column):
        data = {}
        try:
            for col in columns:
                data[col] = np.asarray(read_column(col))
            return np.rec.fromarrays([data[c] for c in columns], names=list(columns))
        except Exception:
            return None
    return None


def _select_spectral_time_without_spectra(conv, arr):
    names = list(arr.dtype.names or [])
    timestamp_field = conv._pick_field_name(names, None, ["timestamp", "time_spectrometer"])
    fallback_field = conv._pick_field_name(names, None, ["time", "t", "unix_time", "unixtime"])

    time_meta = {
        "applied": None,
        "timestamp_field": timestamp_field,
        "fallback_field": fallback_field,
        "suffix": None,
        "first_timestamp_text": None,
        "reason": None,
    }

    if timestamp_field is not None:
        first_text = conv._find_first_nonempty_timestamp_text(arr[timestamp_field])
        time_meta["first_timestamp_text"] = first_text
        suffix = conv._extract_timestamp_suffix(first_text)
        time_meta["suffix"] = suffix
        if suffix in ("UTC", "GPS", "TAI"):
            try:
                t_spec = conv._timestamp_texts_to_unix(arr[timestamp_field], suffix)
                if np.all(np.isfinite(t_spec)):
                    time_meta["applied"] = f"{timestamp_field}:{suffix}->unix"
                    time_meta["reason"] = "timestamp suffix accepted"
                    return np.asarray(t_spec, dtype=float), time_meta
                time_meta["reason"] = f"non-finite values after {suffix} conversion"
            except Exception as exc:
                time_meta["reason"] = f"{suffix} conversion failed: {exc}"
        elif suffix == "PC":
            time_meta["reason"] = "timestamp suffix is PC; falling back to numeric time"
        elif first_text is None:
            time_meta["reason"] = "timestamp field is empty/NaN; falling back to numeric time"
        else:
            time_meta["reason"] = f"unknown timestamp suffix in {first_text!r}; falling back to numeric time"

    if fallback_field is not None:
        t_spec = np.asarray(arr[fallback_field], dtype=float)
        time_meta["applied"] = fallback_field
        if time_meta["reason"] is None:
            time_meta["reason"] = "timestamp field unavailable; using numeric time"
        return t_spec, time_meta

    if timestamp_field is not None:
        raise RuntimeError(
            "timestamp-like field exists but no usable numeric fallback field is available; "
            f"timestamp_field={timestamp_field} available={names} reason={time_meta['reason']}"
        )
    raise RuntimeError(f"no usable spectral time field. available={names}")


def _extract_scan_related_fields(conv, arr, names: Sequence[str], nrow: int) -> Dict[str, Any]:
    scan_field = conv._pick_field_name(
        list(names),
        None,
        ["scan", "scan_id", "scanid", "scan_no", "scan_num", "scan_number", "scan_index", "scanindex", "subscan", "subscan_id"],
    )
    line_index_field = conv._pick_field_name(
        list(names),
        None,
        ["line_index", "lineindex", "line_no", "line_num", "line_number", "line_id", "lineid", "subscan_index"],
    )
    line_label_field = conv._pick_field_name(
        list(names),
        None,
        ["line_label", "linelabel", "line_name", "linename", "scan_label", "scanlabel", "subscan_label"],
    )

    def _take(name: Optional[str]):
        if name is None:
            return None
        try:
            out = np.asarray(arr[name])
        except Exception:
            return None
        if out.ndim != 1:
            out = np.ravel(out)
        if out.size != int(nrow):
            return None
        return out

    line_label = _take(line_label_field)
    if line_label is not None:
        line_label = np.asarray([
            "" if (v is None or str(v).strip().lower() in ("", "nan", "none")) else str(v).strip()
            for v in np.asarray(line_label, dtype=object)
        ], dtype=object)

    return {
        "scan_id": _take(scan_field),
        "line_index": _take(line_index_field),
        "line_label": line_label,
        "scan_field_name": scan_field,
        "line_index_field_name": line_index_field,
        "line_label_field_name": line_label_field,
    }



def read_spectral_timing_only(conv, db, spec_table_name: str):
    """
    Return only spectral timing, OBSMODE-like labels, and optional scan/line identifiers.

    Fast path: ask necstdb for a subset of columns if supported.
    Fallback: read raw records once and access only needed fields without
    materializing the spectral data column.
    """
    table = db.open_table(spec_table_name)
    dtype = conv._get_table_dtype(table)
    if dtype is None or dtype.names is None:
        arr_full = conv._read_structured_array_tolerant(db, spec_table_name)
        t_spec, time_meta = _select_spectral_time_without_spectra(conv, arr_full)
        mode_str, pos_field = conv._extract_pos_labels_from_structured(arr_full, len(t_spec))
        scan_meta = _extract_scan_related_fields(conv, arr_full, list(arr_full.dtype.names or []), len(t_spec))
        return t_spec, mode_str, time_meta, pos_field, scan_meta

    names = list(dtype.names)
    pos_field = conv._pick_field_name(names, None, ["position", "obsmode", "obs_mode", "obsMode", "mode", "state", "status", "label"])
    timestamp_field = conv._pick_field_name(names, None, ["timestamp", "time_spectrometer"])
    fallback_field = conv._pick_field_name(names, None, ["time", "t", "unix_time", "unixtime"])
    scan_meta_fields = _extract_scan_related_fields(conv, np.empty(0, dtype=dtype), names, 0)

    needed_cols = []
    for name in (timestamp_field, fallback_field, pos_field,
                 scan_meta_fields.get("scan_field_name"),
                 scan_meta_fields.get("line_index_field_name"),
                 scan_meta_fields.get("line_label_field_name")):
        if name is not None and name not in needed_cols:
            needed_cols.append(name)

    arr_subset = None
    if needed_cols:
        arr_subset = _try_read_structured_subset(table, needed_cols)

    if arr_subset is None:
        raw = conv._read_table_raw_bytes(table)
        itemsize = int(dtype.itemsize)
        nrec = len(raw) // itemsize
        if nrec <= 0:
            raise RuntimeError(
                f"raw buffer too small for spectral table {spec_table_name}: len={len(raw)} itemsize={itemsize}"
            )
        arr_subset = np.frombuffer(raw[: nrec * itemsize], dtype=dtype)

    t_spec, time_meta = _select_spectral_time_without_spectra(conv, arr_subset)
    mode_str, pos_field = conv._extract_pos_labels_from_structured(arr_subset, len(t_spec))
    scan_meta = _extract_scan_related_fields(conv, arr_subset, list(arr_subset.dtype.names or []), len(t_spec))
    return t_spec, mode_str, time_meta, pos_field, scan_meta


# -----------------------------------------------------------------------------
# Selection / plotting helpers
# -----------------------------------------------------------------------------
def find_off_blocks(mode_str: Sequence[Any]) -> List[Tuple[int, int]]:
    arr = np.asarray(mode_str, dtype=object)
    n = int(arr.size)
    blocks: List[Tuple[int, int]] = []
    i = 0
    while i < n:
        if str(arr[i]).strip().upper() != "OFF":
            i += 1
            continue
        j = i + 1
        while j < n and str(arr[j]).strip().upper() == "OFF":
            j += 1
        blocks.append((i, j))
        i = j
    return blocks


def _normalize_segment_mode(mode: Any) -> str:
    s = str(mode).strip().upper()
    if s == "ON":
        return "ON"
    if s == "OFF":
        return "OFF"
    if s == "HOT":
        return "HOT"
    return "OTHER"


def _slice_by_selection(n: int, off_blocks: Sequence[Tuple[int, int]], selection: str,
                        off_block_start: Optional[int] = None,
                        off_block_end: Optional[int] = None) -> slice:
    selection = str(selection).strip().lower()
    if selection == "all":
        return slice(0, n)
    if selection == "from-first-off":
        if not off_blocks:
            raise ValueError("selection=from-first-off was requested but no OFF block exists")
        return slice(int(off_blocks[0][0]), n)
    if selection == "off-block-range":
        if not off_blocks:
            raise ValueError("selection=off-block-range was requested but no OFF block exists")
        if off_block_start is None:
            raise ValueError("selection=off-block-range requires --off-block-start")
        i0 = int(off_block_start) - 1
        i1 = i0 if off_block_end is None else int(off_block_end) - 1
        if i0 < 0 or i1 < 0 or i0 >= len(off_blocks) or i1 >= len(off_blocks):
            raise ValueError(
                f"OFF block range out of bounds: requested {off_block_start}:{off_block_end}, available 1..{len(off_blocks)}"
            )
        if i1 < i0:
            raise ValueError(f"OFF block range must satisfy start <= end, got {off_block_start}:{off_block_end}")
        return slice(int(off_blocks[i0][0]), int(off_blocks[i1][1]))
    raise ValueError(f"unsupported selection={selection!r}")


def _apply_slice_to_samples(samples: TrajectorySamples, sl: slice) -> TrajectorySamples:
    def maybe_take(x):
        if x is None:
            return None
        return np.asarray(x)[sl]

    pointing = {k: maybe_take(v) if isinstance(v, np.ndarray) else v for k, v in samples.pointing.items()}
    out = TrajectorySamples(
        t_spec=np.asarray(samples.t_spec)[sl],
        mode_str=np.asarray(samples.mode_str, dtype=object)[sl],
        stream_name=samples.stream_name,
        spec_table_name=samples.spec_table_name,
        spec_time_field=samples.spec_time_field,
        spec_time_basis=samples.spec_time_basis,
        spec_time_suffix=samples.spec_time_suffix,
        spec_time_fallback_field=samples.spec_time_fallback_field,
        spec_time_example=samples.spec_time_example,
        pointing=pointing,
        coord_source=samples.coord_source,
        coord_meaning=samples.coord_meaning,
        x=np.asarray(samples.x)[sl],
        y=np.asarray(samples.y)[sl],
        x_label=samples.x_label,
        y_label=samples.y_label,
        frame_name=samples.frame_name,
        coord_meta=dict(samples.coord_meta),
        site=samples.site,
        cmd_x=maybe_take(samples.cmd_x),
        cmd_y=maybe_take(samples.cmd_y),
        cmd_source=samples.cmd_source,
        cmd_meaning=samples.cmd_meaning,
        off_blocks=find_off_blocks(np.asarray(samples.mode_str, dtype=object)[sl]),
        scan_id=maybe_take(samples.scan_id),
        line_index=maybe_take(samples.line_index),
        line_label=maybe_take(samples.line_label),
        scan_field_name=samples.scan_field_name,
        line_index_field_name=samples.line_index_field_name,
        line_label_field_name=samples.line_label_field_name,
    )
    return out


def select_samples_by_off_blocks(samples: TrajectorySamples, selection: str = "all",
                                 off_block_start: Optional[int] = None,
                                 off_block_end: Optional[int] = None) -> TrajectorySamples:
    off_blocks = find_off_blocks(samples.mode_str)
    sl = _slice_by_selection(len(samples.t_spec), off_blocks, selection, off_block_start, off_block_end)
    out = _apply_slice_to_samples(samples, sl)
    out.off_blocks = find_off_blocks(out.mode_str)
    return out


def _make_segments(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    pts = np.column_stack([x, y])
    if pts.shape[0] < 2:
        return np.empty((0, 2, 2), dtype=float)
    return np.stack([pts[:-1], pts[1:]], axis=1)


def _make_mode_masks(mode_str: Sequence[Any]) -> Dict[str, np.ndarray]:
    mode_arr = np.asarray([_normalize_segment_mode(m) for m in mode_str], dtype=object)
    if mode_arr.size <= 1:
        base = np.zeros(0, dtype=bool)
        return {"ON": base, "OFF": base, "HOT": base, "OTHER": base}
    seg_mode = mode_arr[:-1]
    return {key: (seg_mode == key) for key in ("ON", "OFF", "HOT", "OTHER")}


def _unwrap_if_needed(x: np.ndarray, frame_name: str, axis_name: str, unwrap_az: bool) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if not unwrap_az:
        return arr
    if str(frame_name).lower() != "azel":
        return arr
    if str(axis_name).lower() != "az":
        return arr
    return np.rad2deg(np.unwrap(np.deg2rad(arr)))


def _wrap_delta_deg(delta_deg: np.ndarray) -> np.ndarray:
    d = np.asarray(delta_deg, dtype=float)
    return ((d + 180.0) % 360.0) - 180.0


def _frame_default_reverse_x(frame_name: str) -> bool:
    return str(frame_name).strip().lower() in ("radec", "galactic")


def _circular_center_deg(x_deg: np.ndarray) -> float:
    x = np.asarray(x_deg, dtype=float)
    finite = np.isfinite(x)
    if not np.any(finite):
        return np.nan
    theta = np.deg2rad(np.mod(x[finite], 360.0))
    z = np.nanmean(np.exp(1j * theta))
    if (not np.isfinite(z.real)) or (not np.isfinite(z.imag)) or (abs(z) < 1e-12):
        return float(x[finite][0])
    return float(np.mod(np.rad2deg(np.angle(z)), 360.0))


def _prepare_display_xy(
    x: np.ndarray,
    y: np.ndarray,
    frame_name: str,
    *,
    equal_aspect: bool,
    display_lon0_deg: Optional[float] = None,
    longitude_continuous: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    frame = str(frame_name).strip().lower()
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    meta: Dict[str, Any] = {
        "reverse_x": _frame_default_reverse_x(frame),
        "equal_aspect": bool(equal_aspect),
        "lat0_deg": None,
        "lon0_deg": None,
        "cos_lat0": None,
        "aspect_ratio": 1.0,
        "longitude_continuous": bool(longitude_continuous),
    }

    x_disp = x.copy()
    y_disp = y.copy()

    finite = np.isfinite(x) & np.isfinite(y)
    if not np.any(finite):
        return x_disp, y_disp, meta

    if _frame_longitude_x(frame) and (not bool(longitude_continuous)):
        if display_lon0_deg is None or (not np.isfinite(display_lon0_deg)):
            lon0 = _circular_center_deg(x[finite])
        else:
            lon0 = float(display_lon0_deg)
        x_disp = lon0 + _wrap_delta_deg(x - lon0)
    else:
        lon0 = float(np.nanmedian(x[finite]))

    lat0 = float(np.nanmedian(y[finite]))
    cos_lat0 = float(np.cos(np.deg2rad(lat0)))
    if (not np.isfinite(cos_lat0)) or (abs(cos_lat0) < 1e-12):
        cos_lat0 = 1.0

    meta["lat0_deg"] = lat0
    meta["lon0_deg"] = lon0
    meta["cos_lat0"] = cos_lat0
    meta["aspect_ratio"] = 1.0 / cos_lat0 if bool(equal_aspect) else 1.0
    return x_disp, y_disp, meta


def _format_utc_range_from_unix(t_spec: np.ndarray) -> str:
    t = np.asarray(t_spec, dtype=float)
    finite = np.isfinite(t)
    if not np.any(finite):
        return "UTC unknown"
    t0 = float(np.nanmin(t[finite]))
    t1 = float(np.nanmax(t[finite]))
    dt0 = datetime.fromtimestamp(t0, tz=timezone.utc)
    dt1 = datetime.fromtimestamp(t1, tz=timezone.utc)
    return f"UTC {dt0.strftime('%Y-%m-%d %H:%M:%S')} .. {dt1.strftime('%Y-%m-%d %H:%M:%S')}"


def _compute_on_run_ids(mode_str: Sequence[Any]) -> np.ndarray:
    mode_arr = np.asarray([_normalize_segment_mode(m) for m in mode_str], dtype=object)
    out = np.full(mode_arr.size, -1, dtype=int)
    run_id = -1
    i = 0
    while i < mode_arr.size:
        if str(mode_arr[i]) != "ON":
            i += 1
            continue
        j = i + 1
        while j < mode_arr.size and str(mode_arr[j]) == "ON":
            j += 1
        run_id += 1
        out[i:j] = run_id
        i = j
    return out



def _group_key_valid_mask(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=object)
    out = np.zeros(arr.shape, dtype=bool)
    for i, v in enumerate(arr.tolist()):
        if v is None:
            out[i] = False
        elif isinstance(v, (np.floating, float)):
            out[i] = bool(np.isfinite(v))
        elif isinstance(v, (np.integer, int)):
            out[i] = True
        else:
            s = str(v).strip()
            out[i] = (s != "") and (s.lower() not in ("nan", "none"))
    return out



def _select_scan_group_key(samples: TrajectorySamples) -> Tuple[Optional[np.ndarray], Optional[str]]:
    candidates = [
        (samples.scan_id, "scan_id"),
        (samples.line_index, "line_index"),
        (samples.line_label, "line_label"),
    ]
    for values, name in candidates:
        if values is None:
            continue
        arr = np.asarray(values, dtype=object)
        if arr.size <= 0:
            continue
        if np.any(_group_key_valid_mask(arr)):
            return arr, name
    return None, None


def _frame_longitude_x(frame_name: str) -> bool:
    return str(frame_name).strip().lower() in ("azel", "radec", "galactic")


def _segment_mid_lat_deg(y0: np.ndarray, y1: np.ndarray) -> np.ndarray:
    return 0.5 * (np.asarray(y0, dtype=float) + np.asarray(y1, dtype=float))


def _dx_local_deg(x0: np.ndarray, x1: np.ndarray, y_ref_deg: np.ndarray, frame_name: str) -> np.ndarray:
    dx = np.asarray(x0, dtype=float) - np.asarray(x1, dtype=float)
    if _frame_longitude_x(frame_name):
        dx = _wrap_delta_deg(dx)
    return dx * np.cos(np.deg2rad(np.asarray(y_ref_deg, dtype=float)))


def _segment_value_from_point_value(point_value: np.ndarray) -> np.ndarray:
    p = np.asarray(point_value, dtype=float)
    if p.size < 2:
        return np.empty(0, dtype=float)
    return 0.5 * (p[:-1] + p[1:])


def _compute_colored_on_segments(
    x: np.ndarray,
    y: np.ndarray,
    t: np.ndarray,
    on_run_ids: np.ndarray,
    frame_name: str,
    color_by: str,
    cmd_x: Optional[np.ndarray] = None,
    cmd_y: Optional[np.ndarray] = None,
    scan_group_key: Optional[np.ndarray] = None,
    color_sign_mode: str = "signed",
    speed_window_points: int = 2,
    speed_window_mode: str = "secant",
) -> Tuple[np.ndarray, np.ndarray, str, bool]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    on_run_ids = np.asarray(on_run_ids, dtype=int)
    metric = str(color_by).strip().lower()
    sign_mode = str(color_sign_mode).strip().lower()
    if sign_mode not in ("signed", "abs"):
        raise ValueError(f"unsupported color_sign_mode={color_sign_mode!r}")
    speed_mode = str(speed_window_mode).strip().lower()
    if speed_mode not in ("secant", "component_mean", "magnitude_mean"):
        raise ValueError(f"unsupported speed_window_mode={speed_window_mode!r}")
    try:
        speed_win_pts = int(speed_window_points)
    except Exception as exc:
        raise ValueError(f"speed_window_points must be an integer, got {speed_window_points!r}") from exc
    if speed_win_pts < 2:
        raise ValueError(f"speed_window_points must be >= 2, got {speed_window_points!r}")

    segments = _make_segments(x, y)
    if segments.shape[0] <= 0:
        return segments, np.empty(0, dtype=float), "", False

    same_on_run = (on_run_ids[:-1] >= 0) & (on_run_ids[1:] == on_run_ids[:-1])
    same_scan = np.ones(same_on_run.shape, dtype=bool)
    if scan_group_key is not None:
        scan_group_key = np.asarray(scan_group_key, dtype=object)
        if scan_group_key.size != x.size:
            raise ValueError("scan_group_key must have the same length as x/y/t")
        valid_group = _group_key_valid_mask(scan_group_key)
        same_scan = valid_group[:-1] & valid_group[1:] & (scan_group_key[:-1] == scan_group_key[1:])
    seg_mask_base = same_on_run & same_scan

    def _contiguous_true_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
        mask = np.asarray(mask, dtype=bool)
        runs: List[Tuple[int, int]] = []
        i = 0
        nmask = int(mask.size)
        while i < nmask:
            if not bool(mask[i]):
                i += 1
                continue
            j = i + 1
            while j < nmask and bool(mask[j]):
                j += 1
            runs.append((i, j))
            i = j
        return runs

    if metric in ("frame_dx", "frame_dy", "frame_dr"):
        if cmd_x is None or cmd_y is None:
            raise ValueError(f"color_by={metric} requires command-frame coordinates")
        cmd_x = np.asarray(cmd_x, dtype=float)
        cmd_y = np.asarray(cmd_y, dtype=float)
        y_ref = 0.5 * (y + cmd_y)
        point_dx_deg = _dx_local_deg(x, cmd_x, y_ref, frame_name)
        point_dy_deg = y - cmd_y
        point_dr_deg = np.hypot(point_dx_deg, point_dy_deg)
        if metric == "frame_dx":
            seg_value = _segment_value_from_point_value(point_dx_deg) * 3600.0
            label = "frame dx (enc-cmd) [arcsec]"
            signed = True
        elif metric == "frame_dy":
            seg_value = _segment_value_from_point_value(point_dy_deg) * 3600.0
            label = "frame dy (enc-cmd) [arcsec]"
            signed = True
        else:
            seg_value = _segment_value_from_point_value(point_dr_deg) * 3600.0
            label = "frame dr (enc-cmd) [arcsec]"
            signed = False
    elif metric in ("frame_vx", "frame_vy", "frame_speed"):
        if x.size < 2:
            return np.empty((0, 2, 2), dtype=float), np.empty(0, dtype=float), "", False
        dt = np.diff(t)
        mid_lat = _segment_mid_lat_deg(y[:-1], y[1:])
        dx_deg = _dx_local_deg(x[1:], x[:-1], mid_lat, frame_name)
        dy_deg = np.diff(y)
        good_dt = np.isfinite(dt) & (dt > 0.0)
        vx_inst = np.full(dt.shape, np.nan, dtype=float)
        vy_inst = np.full(dt.shape, np.nan, dtype=float)
        vx_inst[good_dt] = dx_deg[good_dt] / dt[good_dt]
        vy_inst[good_dt] = dy_deg[good_dt] / dt[good_dt]

        vx = vx_inst.copy()
        vy = vy_inst.copy()
        if speed_win_pts > 2:
            runs = _contiguous_true_runs(seg_mask_base & good_dt)
            vx[:] = np.nan
            vy[:] = np.nan
            win_seg = speed_win_pts - 1
            left = max(0, (win_seg - 1) // 2)
            right = max(0, win_seg - 1 - left)
            for a, b in runs:
                for i in range(a, b):
                    s0 = max(a, i - left)
                    s1 = min(b - 1, i + right)
                    if s1 < s0:
                        continue
                    if speed_mode == "secant":
                        p0 = s0
                        p1 = s1 + 1
                        dtw = t[p1] - t[p0]
                        if np.isfinite(dtw) and (dtw > 0.0):
                            yref = 0.5 * (y[p0] + y[p1])
                            vx[i] = _dx_local_deg(x[p1], x[p0], yref, frame_name) / dtw
                            vy[i] = (y[p1] - y[p0]) / dtw
                    elif speed_mode == "component_mean":
                        sl = slice(s0, s1 + 1)
                        if np.any(np.isfinite(vx_inst[sl])):
                            vx[i] = float(np.nanmean(vx_inst[sl]))
                        if np.any(np.isfinite(vy_inst[sl])):
                            vy[i] = float(np.nanmean(vy_inst[sl]))
                    else:
                        sl = slice(s0, s1 + 1)
                        sp = np.hypot(vx_inst[sl], vy_inst[sl])
                        if np.any(np.isfinite(sp)):
                            sp_mean = float(np.nanmean(sp))
                            if metric == "frame_vx":
                                if np.any(np.isfinite(vx_inst[sl])):
                                    sgn = float(np.sign(np.nanmean(vx_inst[sl])))
                                    vx[i] = sgn * sp_mean
                            elif metric == "frame_vy":
                                if np.any(np.isfinite(vy_inst[sl])):
                                    sgn = float(np.sign(np.nanmean(vy_inst[sl])))
                                    vy[i] = sgn * sp_mean
                            else:
                                vx[i] = sp_mean
                                vy[i] = 0.0
        if metric == "frame_vx":
            seg_value = vx * 3600.0
            label = "frame vx (enc) [arcsec/s]"
            signed = True
        elif metric == "frame_vy":
            seg_value = vy * 3600.0
            label = "frame vy (enc) [arcsec/s]"
            signed = True
        else:
            seg_value = np.hypot(vx, vy) * 3600.0
            label = "frame speed (enc) [arcsec/s]"
            signed = False
            if speed_win_pts > 2:
                label = f"frame speed (enc, {speed_win_pts}-point {speed_mode}) [arcsec/s]"
    else:
        raise ValueError(f"unsupported color_by={color_by!r}")

    if metric in ("frame_dx", "frame_dy", "frame_vx", "frame_vy") and sign_mode == "abs":
        seg_value = np.abs(seg_value)
        signed = False
        label = f"|{label.split(' [', 1)[0]}| [{label.split(' [', 1)[1]}" if " [" in label else f"|{label}|"

    finite_metric = np.isfinite(seg_value)
    seg_mask = seg_mask_base & finite_metric
    return segments[seg_mask], np.asarray(seg_value[seg_mask], dtype=float), label, signed


def _make_metric_norm(metric_values: np.ndarray, *, signed: bool,
                      color_vmin: Optional[float] = None,
                      color_vmax: Optional[float] = None,
                      color_percentile: Optional[float] = None,
                      color_percentile_mode: str = "central"):
    vals = np.asarray(metric_values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size <= 0:
        return None

    auto_vmin = float(np.nanmin(vals))
    auto_vmax = float(np.nanmax(vals))
    if not np.isfinite(auto_vmin) or not np.isfinite(auto_vmax):
        return None

    user_vmin = None if color_vmin is None else float(color_vmin)
    user_vmax = None if color_vmax is None else float(color_vmax)
    if user_vmin is not None and not np.isfinite(user_vmin):
        raise ValueError(f"color_vmin must be finite, got {color_vmin!r}")
    if user_vmax is not None and not np.isfinite(user_vmax):
        raise ValueError(f"color_vmax must be finite, got {color_vmax!r}")

    pct = None if color_percentile is None else float(color_percentile)
    if pct is not None:
        if not np.isfinite(pct):
            raise ValueError(f"color_percentile must be finite, got {color_percentile!r}")
        if not (0.0 < pct <= 100.0):
            raise ValueError(f"color_percentile must satisfy 0 < q <= 100, got {pct}")
    pct_mode = str(color_percentile_mode).strip().lower()
    if pct_mode not in ("central", "zero_upper"):
        raise ValueError(f"unsupported color_percentile_mode={color_percentile_mode!r}")

    pct_auto_vmin = auto_vmin
    pct_auto_vmax = auto_vmax
    if pct is not None:
        if signed:
            pct_auto_vmax = float(np.nanpercentile(np.abs(vals), pct))
            pct_auto_vmin = -pct_auto_vmax
        else:
            if pct_mode == "central":
                tail = max(0.0, 0.5 * (100.0 - pct))
                pct_auto_vmin = float(np.nanpercentile(vals, tail))
                pct_auto_vmax = float(np.nanpercentile(vals, 100.0 - tail))
            else:
                pct_auto_vmin = 0.0
                pct_auto_vmax = float(np.nanpercentile(vals, pct))

    if signed:
        if user_vmin is None and user_vmax is None:
            vmax_abs = float(max(abs(pct_auto_vmin), abs(pct_auto_vmax)))
            if vmax_abs <= 0.0:
                return None
            return mcolors.TwoSlopeNorm(vmin=-vmax_abs, vcenter=0.0, vmax=vmax_abs)

        if user_vmin is None and user_vmax is not None:
            vmax_abs = abs(float(user_vmax))
            if vmax_abs <= 0.0:
                return None
            return mcolors.TwoSlopeNorm(vmin=-vmax_abs, vcenter=0.0, vmax=vmax_abs)

        if user_vmax is None and user_vmin is not None:
            vmax_abs = abs(float(user_vmin))
            if vmax_abs <= 0.0:
                return None
            return mcolors.TwoSlopeNorm(vmin=-vmax_abs, vcenter=0.0, vmax=vmax_abs)

        vmin = float(user_vmin)
        vmax = float(user_vmax)
        if vmax <= vmin:
            raise ValueError(f"color range must satisfy vmin < vmax, got {vmin} >= {vmax}")
        if vmin < 0.0 < vmax:
            return mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
        return mcolors.Normalize(vmin=vmin, vmax=vmax)

    vmin = pct_auto_vmin if user_vmin is None else float(user_vmin)
    vmax = pct_auto_vmax if user_vmax is None else float(user_vmax)
    if vmax <= vmin:
        raise ValueError(f"color range must satisfy vmin < vmax, got {vmin} >= {vmax}")
    return mcolors.Normalize(vmin=vmin, vmax=vmax)


# -----------------------------------------------------------------------------
# Coordinate conversion
# -----------------------------------------------------------------------------
def _convert_frame(conv, site: SimpleSite, t_spec: np.ndarray,
                   az_deg: np.ndarray, el_deg: np.ndarray,
                   frame_name: str) -> Tuple[np.ndarray, np.ndarray, str, str, Dict[str, Any]]:
    frame_name = str(frame_name).strip().lower()
    meta: Dict[str, Any] = {}
    if frame_name == "azel":
        return np.asarray(az_deg, dtype=float), np.asarray(el_deg, dtype=float), "Az [deg]", "El [deg]", meta
    if frame_name == "radec":
        radec = conv.normalize_radec(
            wcs_df=None,
            t_spec=np.asarray(t_spec, dtype=float),
            beam_az_deg=np.asarray(az_deg, dtype=float),
            beam_el_deg=np.asarray(el_deg, dtype=float),
            site=site,
            radec_method="azel",
            apply_refraction=False,
            press_hpa=None,
            temp_c=None,
            humid_pct=None,
            obswl_um=None,
        )
        meta["used_wcs"] = False
        return np.asarray(radec["ra_deg"], dtype=float), np.asarray(radec["dec_deg"], dtype=float), "RA [deg]", "Dec [deg]", meta
    if frame_name == "galactic":
        radec = conv.normalize_radec(
            wcs_df=None,
            t_spec=np.asarray(t_spec, dtype=float),
            beam_az_deg=np.asarray(az_deg, dtype=float),
            beam_el_deg=np.asarray(el_deg, dtype=float),
            site=site,
            radec_method="azel",
            apply_refraction=False,
            press_hpa=None,
            temp_c=None,
            humid_pct=None,
            obswl_um=None,
        )
        gal = conv._gal_from_radec(radec["ra_deg"], radec["dec_deg"])
        meta["used_wcs"] = False
        return np.asarray(gal["glon_deg"], dtype=float), np.asarray(gal["glat_deg"], dtype=float), "l [deg]", "b [deg]", meta
    raise ValueError(f"unsupported coord frame={frame_name!r}")



def _normalize_coord_source_local(coord_source: str) -> str:
    src = str(coord_source).strip().lower()
    mapping = {
        "true": "beam",
        "encoder": "raw_encoder",
        "altaz": "raw_altaz",
        "cmd": "corrected_cmd",
    }
    return mapping.get(src, src)


def _coord_source_to_columns(coord_source: str) -> Tuple[str, str, str]:
    src = _normalize_coord_source_local(coord_source)
    if src == "beam":
        return "real_beam_az_deg", "real_beam_el_deg", "beam-center Az/El from corrected encoder boresight + beam offset"
    if src == "boresight":
        return "real_boresight_az_deg", "real_boresight_el_deg", "corrected encoder boresight Az/El without beam offset"
    if src == "raw_encoder":
        return "raw_encoder_beam_az_deg", "raw_encoder_beam_el_deg", "beam-offset Az/El from raw encoder lon/lat"
    if src == "raw_altaz":
        return "raw_altaz_beam_az_deg", "raw_altaz_beam_el_deg", "beam-offset Az/El from raw altaz/cmd lon/lat"
    if src == "corrected_cmd":
        return "corrected_cmd_beam_az_deg", "corrected_cmd_beam_el_deg", "beam-offset Az/El from corrected altaz/cmd boresight"
    raise ValueError(f"unsupported coord source={coord_source!r}")


def _csv_float_or_nan(value: Any) -> float:
    if value is None:
        return np.nan
    s = str(value).strip()
    if s == "" or s.lower() == "nan":
        return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan


def _compute_off_block_indices(mode_str: Sequence[Any]) -> np.ndarray:
    mode_arr = np.asarray(mode_str, dtype=object)
    out = np.zeros(mode_arr.size, dtype=int)
    block_id = 0
    i = 0
    while i < mode_arr.size:
        if str(mode_arr[i]).strip().upper() != "OFF":
            i += 1
            continue
        j = i + 1
        while j < mode_arr.size and str(mode_arr[j]).strip().upper() == "OFF":
            j += 1
        block_id += 1
        out[i:j] = block_id
        i = j
    return out


def _csv_any_or_empty(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (np.floating, float)) and (not np.isfinite(value)):
        return ""
    s = str(value).strip()
    if s.lower() in ("nan", "none"):
        return ""
    return s


def export_trajectory_csv(samples: TrajectorySamples, csv_path: str) -> pathlib.Path:
    if samples.site is None:
        raise ValueError("samples.site is required for CSV export")

    path = pathlib.Path(csv_path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    site = samples.site
    t = np.asarray(samples.t_spec, dtype=float)
    mode = np.asarray(samples.mode_str, dtype=object)
    p = samples.pointing

    def arr(name: str) -> np.ndarray:
        x = p.get(name)
        if x is None:
            return np.full(t.shape, np.nan, dtype=float)
        return np.asarray(x, dtype=float)

    off_idx = _compute_off_block_indices(mode)

    columns = [
        "sample_index", "unix_time", "iso_time", "mode", "mode_norm", "off_block_index",
        "scan_id", "line_index", "line_label",
        "stream_name", "spec_table_name",
        "site_lat_deg", "site_lon_deg", "site_elev_m",
        "enc_az_deg", "enc_el_deg",
        "cmd_az_deg", "cmd_el_deg",
        "dlon_deg", "dlat_deg",
        "real_boresight_az_deg", "real_boresight_el_deg",
        "real_beam_az_deg", "real_beam_el_deg",
        "raw_encoder_beam_az_deg", "raw_encoder_beam_el_deg",
        "raw_altaz_beam_az_deg", "raw_altaz_beam_el_deg",
        "corrected_cmd_boresight_az_deg", "corrected_cmd_boresight_el_deg",
        "corrected_cmd_beam_az_deg", "corrected_cmd_beam_el_deg",
        "beam_dx_arcsec", "beam_dy_arcsec", "beam_rot_deg",
        "pe_x_arcsec", "pe_y_arcsec", "pe_r_arcsec",
    ]

    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        for i in range(t.size):
            writer.writerow([
                i,
                f"{float(t[i]):.9f}",
                str(np.datetime64(int(round(float(t[i]) * 1e9)), "ns")),
                str(mode[i]),
                _normalize_segment_mode(mode[i]),
                int(off_idx[i]),
                _csv_any_or_empty(None if samples.scan_id is None else np.asarray(samples.scan_id, dtype=object)[i]),
                _csv_any_or_empty(None if samples.line_index is None else np.asarray(samples.line_index, dtype=object)[i]),
                _csv_any_or_empty(None if samples.line_label is None else np.asarray(samples.line_label, dtype=object)[i]),
                samples.stream_name,
                samples.spec_table_name,
                f"{float(site.lat_deg):.12f}",
                f"{float(site.lon_deg):.12f}",
                f"{float(site.elev_m):.6f}",
                f"{float(arr('az_enc_t')[i]):.12f}",
                f"{float(arr('el_enc_t')[i]):.12f}",
                f"{float(arr('az_cmd_t')[i]):.12f}",
                f"{float(arr('el_cmd_t')[i]):.12f}",
                f"{float(arr('dlon_t')[i]):.12f}",
                f"{float(arr('dlat_t')[i]):.12f}",
                f"{float(arr('real_boresight_az_deg')[i]):.12f}",
                f"{float(arr('real_boresight_el_deg')[i]):.12f}",
                f"{float(arr('real_beam_az_deg')[i]):.12f}",
                f"{float(arr('real_beam_el_deg')[i]):.12f}",
                f"{float(arr('raw_encoder_beam_az_deg')[i]):.12f}",
                f"{float(arr('raw_encoder_beam_el_deg')[i]):.12f}",
                f"{float(arr('raw_altaz_beam_az_deg')[i]):.12f}",
                f"{float(arr('raw_altaz_beam_el_deg')[i]):.12f}",
                f"{float(arr('corrected_cmd_boresight_az_deg')[i]):.12f}",
                f"{float(arr('corrected_cmd_boresight_el_deg')[i]):.12f}",
                f"{float(arr('corrected_cmd_beam_az_deg')[i]):.12f}",
                f"{float(arr('corrected_cmd_beam_el_deg')[i]):.12f}",
                f"{float(arr('beam_dx_arcsec')[i]):.12f}",
                f"{float(arr('beam_dy_arcsec')[i]):.12f}",
                f"{float(arr('beam_rot_deg')[i]):.12f}",
                f"{float(arr('pe_x_arcsec')[i]):.12f}",
                f"{float(arr('pe_y_arcsec')[i]):.12f}",
                f"{float(arr('pe_r_arcsec')[i]):.12f}",
            ])
    return path


def _convert_frame_local(site: SimpleSite, t_spec: np.ndarray,
                         az_deg: np.ndarray, el_deg: np.ndarray,
                         frame_name: str) -> Tuple[np.ndarray, np.ndarray, str, str, Dict[str, Any]]:
    frame_name = str(frame_name).strip().lower()
    meta: Dict[str, Any] = {"used_wcs": False, "converter": False}
    az = np.asarray(az_deg, dtype=float)
    el = np.asarray(el_deg, dtype=float)
    t = np.asarray(t_spec, dtype=float)
    if frame_name == "azel":
        return az, el, "Az [deg]", "El [deg]", meta

    from astropy.time import Time
    from astropy.coordinates import AltAz, SkyCoord
    from astropy import units as u

    altaz_frame = AltAz(obstime=Time(t, format="unix"), location=site.to_earthlocation())
    sc_altaz = SkyCoord(az=az * u.deg, alt=el * u.deg, frame=altaz_frame)
    icrs = sc_altaz.icrs
    if frame_name == "radec":
        return np.asarray(icrs.ra.deg, dtype=float), np.asarray(icrs.dec.deg, dtype=float), "RA [deg]", "Dec [deg]", meta
    if frame_name == "galactic":
        gal = icrs.galactic
        return np.asarray(gal.l.deg, dtype=float), np.asarray(gal.b.deg, dtype=float), "l [deg]", "b [deg]", meta
    raise ValueError(f"unsupported coord frame={frame_name!r}")


def load_trajectory_samples_from_csv(csv_path: str,
                                     coord_frame: str = "azel",
                                     coord_source: str = "beam",
                                     overlay_cmd: bool = False,
                                     cmd_source: str = "corrected_cmd",
                                     selection: str = "all",
                                     off_block_start: Optional[int] = None,
                                     off_block_end: Optional[int] = None) -> TrajectorySamples:
    path = pathlib.Path(csv_path).expanduser().resolve()
    rows = []
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    if not rows:
        raise ValueError(f"CSV is empty: {path}")

    def col(name: str, numeric: bool = True):
        if numeric:
            return np.asarray([_csv_float_or_nan(r.get(name)) for r in rows], dtype=float)
        return np.asarray([str(r.get(name, "")) for r in rows], dtype=object)

    def obj_col(name: str) -> Optional[np.ndarray]:
        if rows and (name not in rows[0]):
            return None
        arr = np.asarray([_csv_any_or_empty(r.get(name)) for r in rows], dtype=object)
        if not np.any(_group_key_valid_mask(arr)):
            return None
        return arr

    t_spec = col("unix_time")
    mode_str = col("mode", numeric=False)
    stream_name = str(rows[0].get("stream_name", "CSV"))
    spec_table_name = str(rows[0].get("spec_table_name", path.name))
    site = SimpleSite(
        lat_deg=float(_csv_float_or_nan(rows[0].get("site_lat_deg"))),
        lon_deg=float(_csv_float_or_nan(rows[0].get("site_lon_deg"))),
        elev_m=float(_csv_float_or_nan(rows[0].get("site_elev_m"))),
    )

    pointing = {
        "az_enc_t": col("enc_az_deg"),
        "el_enc_t": col("enc_el_deg"),
        "az_cmd_t": col("cmd_az_deg"),
        "el_cmd_t": col("cmd_el_deg"),
        "dlon_t": col("dlon_deg"),
        "dlat_t": col("dlat_deg"),
        "real_boresight_az_deg": col("real_boresight_az_deg"),
        "real_boresight_el_deg": col("real_boresight_el_deg"),
        "boresight_az": col("real_boresight_az_deg"),
        "boresight_el": col("real_boresight_el_deg"),
        "real_beam_az_deg": col("real_beam_az_deg"),
        "real_beam_el_deg": col("real_beam_el_deg"),
        "beam_az_deg": col("real_beam_az_deg"),
        "beam_el_deg": col("real_beam_el_deg"),
        "raw_encoder_beam_az_deg": col("raw_encoder_beam_az_deg"),
        "raw_encoder_beam_el_deg": col("raw_encoder_beam_el_deg"),
        "raw_altaz_beam_az_deg": col("raw_altaz_beam_az_deg"),
        "raw_altaz_beam_el_deg": col("raw_altaz_beam_el_deg"),
        "corrected_cmd_boresight_az_deg": col("corrected_cmd_boresight_az_deg"),
        "corrected_cmd_boresight_el_deg": col("corrected_cmd_boresight_el_deg"),
        "corrected_cmd_beam_az_deg": col("corrected_cmd_beam_az_deg"),
        "corrected_cmd_beam_el_deg": col("corrected_cmd_beam_el_deg"),
        "beam_dx_arcsec": col("beam_dx_arcsec"),
        "beam_dy_arcsec": col("beam_dy_arcsec"),
        "beam_rot_deg": col("beam_rot_deg"),
        "pe_x_arcsec": col("pe_x_arcsec"),
        "pe_y_arcsec": col("pe_y_arcsec"),
        "pe_r_arcsec": col("pe_r_arcsec"),
    }

    coord_source_norm = _normalize_coord_source_local(coord_source)
    az_col, el_col, meaning = _coord_source_to_columns(coord_source_norm)
    x, y, x_label, y_label, coord_meta = _convert_frame_local(site, t_spec, pointing[az_col], pointing[el_col], coord_frame)

    csv_scan_id = obj_col("scan_id")
    csv_line_index = obj_col("line_index")
    csv_line_label = obj_col("line_label")

    cmd_x = cmd_y = None
    cmd_meaning = None
    cmd_source_norm = None
    if overlay_cmd:
        cmd_source_norm = _normalize_coord_source_local(cmd_source)
        cmd_az_col, cmd_el_col, cmd_meaning = _coord_source_to_columns(cmd_source_norm)
        cmd_x, cmd_y, _, _, _ = _convert_frame_local(site, t_spec, pointing[cmd_az_col], pointing[cmd_el_col], coord_frame)

    samples = TrajectorySamples(
        t_spec=t_spec,
        mode_str=mode_str,
        stream_name=stream_name,
        spec_table_name=spec_table_name,
        spec_time_field="unix_time",
        spec_time_basis="csv:unix_time",
        spec_time_suffix=None,
        spec_time_fallback_field=None,
        spec_time_example=None,
        pointing=pointing,
        coord_source=str(coord_source_norm),
        coord_meaning=meaning,
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        x_label=str(x_label),
        y_label=str(y_label),
        frame_name=str(coord_frame).strip().lower(),
        coord_meta=dict(coord_meta),
        site=site,
        cmd_x=None if cmd_x is None else np.asarray(cmd_x, dtype=float),
        cmd_y=None if cmd_y is None else np.asarray(cmd_y, dtype=float),
        cmd_source=cmd_source_norm,
        cmd_meaning=cmd_meaning,
        off_blocks=find_off_blocks(mode_str),
        scan_id=csv_scan_id,
        line_index=csv_line_index,
        line_label=csv_line_label,
        scan_field_name="scan_id" if csv_scan_id is not None else None,
        line_index_field_name="line_index" if csv_line_index is not None else None,
        line_label_field_name="line_label" if csv_line_label is not None else None,
    )
    return select_samples_by_off_blocks(samples, selection=selection, off_block_start=off_block_start, off_block_end=off_block_end)


# -----------------------------------------------------------------------------
# High-level API
# -----------------------------------------------------------------------------
def load_trajectory_samples(rawdata: str,
                            spectrometer_config: Optional[str] = None,
                            stream_names: Optional[Sequence[str]] = None,
                            telescope: str = "OMU1P85M",
                            db_namespace: Optional[str] = None,
                            coord_frame: str = "azel",
                            coord_source: str = "beam",
                            overlay_cmd: bool = False,
                            cmd_source: str = "corrected_cmd",
                            azel_correction_apply: Optional[str] = None,
                            encoder_time_col: str = "time",
                            altaz_time_col: str = "time",
                            encoder_shift_sec: float = 0.0,
                            interp_extrap: str = "hold",
                            site_lat: Optional[float] = None,
                            site_lon: Optional[float] = None,
                            site_elev: Optional[float] = None,
                            converter_path: Optional[str] = None,
                            selection: str = "all",
                            off_block_start: Optional[int] = None,
                            off_block_end: Optional[int] = None,
                            telescope_explicit: bool = False,
                            db_namespace_explicit: bool = False,
                            encoder_time_col_explicit: bool = False,
                            altaz_time_col_explicit: bool = False,
                            encoder_shift_sec_explicit: bool = False) -> TrajectorySamples:
    conv = load_converter_module(converter_path)

    # Build an argparse-like namespace because the converter helpers expect args.*.
    args = argparse.Namespace(
        rawdata=str(rawdata),
        spectrometer_config=spectrometer_config,
        stream_names=list(stream_names) if stream_names is not None else None,
        telescope=str(telescope),
        db_namespace=db_namespace,
        encoder_table=None,
        encoder_table_suffix=None,
        altaz_table=None,
        altaz_table_suffix=None,
        weather_inside_table=None,
        weather_inside_table_suffix=None,
        weather_inside_time_col=None,
        weather_outside_table=None,
        weather_outside_table_suffix=None,
        weather_outside_time_col=None,
        weather_table="auto",
        weather_time_col="time",
        encoder_time_col=str(encoder_time_col),
        altaz_time_col=str(altaz_time_col),
        site_lat=site_lat,
        site_lon=site_lon,
        site_elev=site_elev,
        spectral="xffts-board1",
        nchan=2**15,
        if0_ghz=0.0,
        if1_ghz=2.5,
        lo1_ghz=109.8,
        lo2_ghz=4.0,
        restfreq_ghz=None,
        channel_slice=None,
        vlsrk_kms_slice=None,
        radec_method=None,
        radec_azel_source=coord_source,
        azel_correction_apply=azel_correction_apply,
        encoder_shift_sec=float(encoder_shift_sec),
        met_pressure_hpa=None,
        met_temperature_c=10.0,
        met_humidity_pct=50.0,
        object=None,
        out=None,
        strict_config=False,
        thot_k=None,
        thot_default_k=None,
        thot_min_k=None,
        thot_max_k=None,
        tamb_k=None,
        tamb_default_k=None,
        tamb_min_k=None,
        tamb_max_k=None,
        inside_default_temperature_c=None,
        inside_temperature_min_c=None,
        inside_temperature_max_c=None,
        outside_default_temperature_c=None,
        outside_default_pressure_hpa=None,
        outside_default_humidity_pct=None,
        outside_temperature_min_c=None,
        outside_temperature_max_c=None,
        outside_pressure_min_hpa=None,
        outside_pressure_max_hpa=None,
        outside_humidity_min_pct=None,
        outside_humidity_max_pct=None,
    )
    argv: List[str] = []

    rawdata_path = pathlib.Path(rawdata).expanduser().resolve()
    site_info = conv._resolve_site_from_cli_or_rawdata(args, rawdata_path)

    if spectrometer_config:
        config_dict = conv.load_spectrometer_config(spectrometer_config)
    else:
        config_dict = conv.build_legacy_single_stream_config(args)

    global_cfg = conv._canonicalize_global_coordinate_settings(dict(config_dict.get("global", {}) or {}))
    config_dict["global"] = global_cfg

    coord_source_norm = conv._normalize_output_azel_source(coord_source)
    if azel_correction_apply is None:
        azel_correction_apply_norm = conv._normalize_legacy_correction_mode(
            global_cfg.get("azel_correction_apply", global_cfg.get("boresight_correction_apply", "subtract"))
        )
    else:
        azel_correction_apply_norm = conv._normalize_legacy_correction_mode(azel_correction_apply)

    argv = []
    if bool(telescope_explicit):
        argv.extend(["--telescope", str(telescope)])
    if bool(encoder_time_col_explicit):
        argv.extend(["--encoder-time-col", str(encoder_time_col)])
    if bool(altaz_time_col_explicit):
        argv.extend(["--altaz-time-col", str(altaz_time_col)])
    if bool(encoder_shift_sec_explicit):
        argv.extend(["--encoder-shift-sec", str(float(encoder_shift_sec))])
    if bool(db_namespace_explicit) and (db_namespace is not None):
        argv.extend(["--db-namespace", str(db_namespace)])

    runtime = conv._resolve_runtime_naming(args, config_dict, argv=argv)
    db_namespace = runtime["db_namespace"]
    telescope_name = runtime["telescope"]
    encoder_time_col = runtime["encoder_time_col"]
    altaz_time_col = runtime["altaz_time_col"]
    encoder_shift_sec = float(runtime["encoder_shift_sec"])

    selected_streams = conv._select_streams_for_convert(config_dict.get("streams", []), explicit_stream_names=args.stream_names)
    if not selected_streams:
        raise RuntimeError("no enabled/use_for_convert stream is available")
    if len(selected_streams) > 1:
        print(
            "[info] multiple streams are available for timing; using the first selected stream: {}".format(
                getattr(selected_streams[0], "name", "UNKNOWN")
            )
        )
    stream = selected_streams[0]
    spec_table_name = conv._resolve_spec_table_name(db_namespace, telescope_name, stream)

    db = conv.necstdb.opendb(str(rawdata_path))
    enc_table = conv._resolve_prefixed_table_name(
        db_namespace, telescope_name, runtime.get("encoder_table"), runtime.get("encoder_table_suffix"), "ctrl-antenna-encoder"
    )
    alt_table = conv._resolve_prefixed_table_name(
        db_namespace, telescope_name, runtime.get("altaz_table"), runtime.get("altaz_table_suffix"), "ctrl-antenna-altaz"
    )
    arr_enc = conv._read_structured_array_tolerant(db, enc_table)
    arr_alt = conv._read_structured_array_tolerant(db, alt_table)

    t_spec_unsorted, mode_unsorted, time_meta, pos_field, scan_meta = read_spectral_timing_only(conv, db, spec_table_name)
    order = np.argsort(np.asarray(t_spec_unsorted, dtype=float))
    t_spec = np.asarray(t_spec_unsorted, dtype=float)[order]
    mode_str = np.asarray(mode_unsorted, dtype=object)[order]
    scan_id = None if scan_meta.get("scan_id") is None else np.asarray(scan_meta["scan_id"])[order]
    line_index = None if scan_meta.get("line_index") is None else np.asarray(scan_meta["line_index"])[order]
    line_label = None if scan_meta.get("line_label") is None else np.asarray(scan_meta["line_label"], dtype=object)[order]

    pointing = conv.normalize_pointing(
        t_spec=t_spec,
        arr_enc=arr_enc,
        arr_alt=arr_alt,
        encoder_time_col=encoder_time_col,
        altaz_time_col=altaz_time_col,
        extrap=str(interp_extrap),
        encoder_shift_sec=float(encoder_shift_sec),
        azel_correction_apply=str(azel_correction_apply_norm),
    )

    beam_applied = conv.apply_beam_offset(pointing["boresight_az"], pointing["boresight_el"], stream.beam)
    pointing["beam_az_deg"] = np.asarray(beam_applied["beam_az_deg"], dtype=float)
    pointing["beam_el_deg"] = np.asarray(beam_applied["beam_el_deg"], dtype=float)
    pointing["beam_dx_arcsec"] = np.asarray(beam_applied["beam_dx_arcsec"], dtype=float)
    pointing["beam_dy_arcsec"] = np.asarray(beam_applied["beam_dy_arcsec"], dtype=float)
    pointing["beam_rot_deg"] = np.asarray(beam_applied["beam_rot_deg"], dtype=float)
    pointing["real_boresight_az_deg"] = np.asarray(pointing["boresight_az"], dtype=float)
    pointing["real_boresight_el_deg"] = np.asarray(pointing["boresight_el"], dtype=float)
    pointing["real_beam_az_deg"] = np.asarray(pointing["beam_az_deg"], dtype=float)
    pointing["real_beam_el_deg"] = np.asarray(pointing["beam_el_deg"], dtype=float)

    raw_encoder_applied = conv.apply_beam_offset(pointing["az_enc_t"], pointing["el_enc_t"], stream.beam)
    pointing["raw_encoder_beam_az_deg"] = np.asarray(raw_encoder_applied["beam_az_deg"], dtype=float)
    pointing["raw_encoder_beam_el_deg"] = np.asarray(raw_encoder_applied["beam_el_deg"], dtype=float)

    raw_altaz_applied = conv.apply_beam_offset(pointing["az_cmd_t"], pointing["el_cmd_t"], stream.beam)
    pointing["raw_altaz_beam_az_deg"] = np.asarray(raw_altaz_applied["beam_az_deg"], dtype=float)
    pointing["raw_altaz_beam_el_deg"] = np.asarray(raw_altaz_applied["beam_el_deg"], dtype=float)

    corrected_cmd_boresight_az, corrected_cmd_boresight_el = conv.apply_azel_correction(
        np.asarray(pointing["az_cmd_t"], dtype=float),
        np.asarray(pointing["el_cmd_t"], dtype=float),
        np.asarray(pointing["dlon_t"], dtype=float),
        np.asarray(pointing["dlat_t"], dtype=float),
        str(azel_correction_apply_norm),
    )
    pointing["corrected_cmd_boresight_az_deg"] = np.asarray(corrected_cmd_boresight_az, dtype=float)
    pointing["corrected_cmd_boresight_el_deg"] = np.asarray(corrected_cmd_boresight_el, dtype=float)
    corrected_cmd_applied = conv.apply_beam_offset(corrected_cmd_boresight_az, corrected_cmd_boresight_el, stream.beam)
    pointing["corrected_cmd_beam_az_deg"] = np.asarray(corrected_cmd_applied["beam_az_deg"], dtype=float)
    pointing["corrected_cmd_beam_el_deg"] = np.asarray(corrected_cmd_applied["beam_el_deg"], dtype=float)

    coord_src = conv.select_radec_azel_source(
        pointing=pointing,
        beam=stream.beam,
        radec_azel_source=coord_source_norm,
        azel_correction_apply=str(azel_correction_apply_norm),
    )

    site = SimpleSite(
        lat_deg=float(site_info["lat_deg"]),
        lon_deg=float(site_info["lon_deg"]),
        elev_m=float(site_info["elev_m"]),
    )
    x, y, x_label, y_label, coord_meta = _convert_frame(
        conv=conv,
        site=site,
        t_spec=t_spec,
        az_deg=np.asarray(coord_src["az_deg"], dtype=float),
        el_deg=np.asarray(coord_src["el_deg"], dtype=float),
        frame_name=coord_frame,
    )

    cmd_x = cmd_y = None
    cmd_meaning = None
    cmd_source_norm = None
    if overlay_cmd:
        cmd_source_norm = conv._normalize_output_azel_source(cmd_source)
        cmd_src = conv.select_radec_azel_source(
            pointing=pointing,
            beam=stream.beam,
            radec_azel_source=cmd_source_norm,
            azel_correction_apply=str(azel_correction_apply_norm),
        )
        cmd_x, cmd_y, _, _, _ = _convert_frame(
            conv=conv,
            site=site,
            t_spec=t_spec,
            az_deg=np.asarray(cmd_src["az_deg"], dtype=float),
            el_deg=np.asarray(cmd_src["el_deg"], dtype=float),
            frame_name=coord_frame,
        )
        cmd_meaning = str(cmd_src.get("meaning", ""))

    samples = TrajectorySamples(
        t_spec=t_spec,
        mode_str=mode_str,
        stream_name=str(stream.name),
        spec_table_name=str(spec_table_name),
        spec_time_field=str(time_meta.get("timestamp_field") or time_meta.get("fallback_field") or "UNKNOWN"),
        spec_time_basis=str(time_meta.get("applied") or "UNKNOWN"),
        spec_time_suffix=time_meta.get("suffix"),
        spec_time_fallback_field=time_meta.get("fallback_field"),
        spec_time_example=time_meta.get("first_timestamp_text"),
        pointing=pointing,
        coord_source=str(coord_src.get("source", coord_source_norm)),
        coord_meaning=str(coord_src.get("meaning", "")),
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        x_label=str(x_label),
        y_label=str(y_label),
        frame_name=str(coord_frame).strip().lower(),
        coord_meta=dict(coord_meta),
        site=site,
        cmd_x=None if cmd_x is None else np.asarray(cmd_x, dtype=float),
        cmd_y=None if cmd_y is None else np.asarray(cmd_y, dtype=float),
        cmd_source=cmd_source_norm,
        cmd_meaning=cmd_meaning,
        off_blocks=find_off_blocks(mode_str),
        scan_id=scan_id,
        line_index=line_index,
        line_label=line_label,
        scan_field_name=scan_meta.get("scan_field_name"),
        line_index_field_name=scan_meta.get("line_index_field_name"),
        line_label_field_name=scan_meta.get("line_label_field_name"),
    )
    return select_samples_by_off_blocks(samples, selection=selection, off_block_start=off_block_start, off_block_end=off_block_end)


# -----------------------------------------------------------------------------
# Plotting API
# -----------------------------------------------------------------------------
def _find_mode_runs(mode: Sequence[Any]) -> List[Tuple[int, int, str]]:
    arr = np.asarray(mode, dtype=object)
    runs: List[Tuple[int, int, str]] = []
    n = int(arr.size)
    if n <= 0:
        return runs
    i = 0
    while i < n:
        key = _normalize_segment_mode(arr[i])
        j = i + 1
        while j < n and _normalize_segment_mode(arr[j]) == key:
            j += 1
        runs.append((i, j, key))
        i = j
    return runs


def _draw_quiver_arrows(ax: plt.Axes, x0: np.ndarray, y0: np.ndarray, dx: np.ndarray, dy: np.ndarray,
                        *, color: str, alpha: float, width: float) -> bool:
    x0 = np.asarray(x0, dtype=float)
    y0 = np.asarray(y0, dtype=float)
    dx = np.asarray(dx, dtype=float)
    dy = np.asarray(dy, dtype=float)
    good = np.isfinite(x0) & np.isfinite(y0) & np.isfinite(dx) & np.isfinite(dy) & ((dx * dx + dy * dy) > 0.0)
    if not np.any(good):
        return False
    ax.quiver(
        x0[good],
        y0[good],
        dx[good],
        dy[good],
        angles='xy',
        scale_units='xy',
        scale=1.0,
        color=str(color),
        alpha=float(alpha),
        width=float(width),
        headwidth=4.0,
        headlength=5.0,
        headaxislength=4.5,
        pivot='mid',
        zorder=4,
    )
    return True


def _draw_direction_arrows_samples(ax: plt.Axes, x: np.ndarray, y: np.ndarray, *, arrow_every: int,
                                   arrow_alpha: float, arrow_color: str, arrow_width: float,
                                   arrow_length_frac: float) -> bool:
    if x.size < 2:
        return False
    step_arrow = max(int(arrow_every), 1)
    arrow_idx = np.arange(0, max(x.size - 1, 0), step_arrow, dtype=int)
    if arrow_idx.size <= 0:
        return False
    dx = x[arrow_idx + 1] - x[arrow_idx]
    dy = y[arrow_idx + 1] - y[arrow_idx]
    norm = np.hypot(dx, dy)
    good = np.isfinite(dx) & np.isfinite(dy) & np.isfinite(norm) & (norm > 0.0)
    if not np.any(good):
        return False
    span_x = float(np.nanmax(x) - np.nanmin(x)) if x.size else 0.0
    span_y = float(np.nanmax(y) - np.nanmin(y)) if y.size else 0.0
    span = max(span_x, span_y)
    if (not np.isfinite(span)) or (span <= 0.0):
        span = 1.0
    arrow_len = max(float(arrow_length_frac), 0.0) * span
    if not np.isfinite(arrow_len) or arrow_len <= 0.0:
        arrow_len = 0.015 * span
    ux = dx[good] / norm[good]
    uy = dy[good] / norm[good]
    return _draw_quiver_arrows(
        ax,
        x[arrow_idx][good],
        y[arrow_idx][good],
        ux * arrow_len,
        uy * arrow_len,
        color=arrow_color,
        alpha=arrow_alpha,
        width=arrow_width,
    )


def _draw_direction_arrows_segment(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    mode: Sequence[Any],
    *,
    arrow_alpha: float,
    arrow_color: str,
    arrow_width: float,
    arrow_length_frac: float,
) -> bool:
    """
    Draw one arrow per contiguous run of *drawn line segments* with the same mode.

    Notes
    -----
    - The colored trajectory uses segment-mode = mode[:-1].
    - Therefore segment arrows must also be built from segment runs, not point runs.
    - No extra arrows are drawn at mode boundaries.
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mode_arr = np.asarray(mode, dtype=object)

    npts = int(x.size)
    if npts < 2 or mode_arr.size < 2:
        return False

    # Segment i is the visible line from point i -> i+1,
    # and its mode/color is determined by mode[i].
    seg_mode = np.asarray([_normalize_segment_mode(m) for m in mode_arr[:-1]], dtype=object)
    nseg = int(seg_mode.size)
    if nseg <= 0:
        return False

    span_x = float(np.nanmax(x) - np.nanmin(x)) if npts else 0.0
    span_y = float(np.nanmax(y) - np.nanmin(y)) if npts else 0.0
    span = max(span_x, span_y)
    if (not np.isfinite(span)) or (span <= 0.0):
        span = 1.0

    draw_len = max(float(arrow_length_frac), 0.0) * span
    if (not np.isfinite(draw_len)) or (draw_len <= 0.0):
        draw_len = 0.015 * span

    # Reject only truly degenerate local motion.
    move_min = max(1.0e-12 * span, 1.0e-15)

    # Build contiguous runs in segment-index space: [s, e)
    runs = []
    i = 0
    while i < nseg:
        key = str(seg_mode[i])
        j = i + 1
        while j < nseg and str(seg_mode[j]) == key:
            j += 1
        runs.append((i, j, key))
        i = j

    starts_x = []
    starts_y = []
    vec_x = []
    vec_y = []

    # One arrow per run, aligned with a local drawn segment near the middle.
    for s, e, _key in runs:
        mid_seg = int((s + e - 1) // 2)
        if mid_seg < 0 or mid_seg >= nseg:
            continue

        dx = float(x[mid_seg + 1] - x[mid_seg])
        dy = float(y[mid_seg + 1] - y[mid_seg])
        dist = float(np.hypot(dx, dy))
        if (not np.isfinite(dist)) or (dist <= move_min):
            continue

        starts_x.append(0.5 * float(x[mid_seg] + x[mid_seg + 1]))
        starts_y.append(0.5 * float(y[mid_seg] + y[mid_seg + 1]))
        vec_x.append(draw_len * dx / dist)
        vec_y.append(draw_len * dy / dist)

    if len(starts_x) <= 0:
        return False

    return _draw_quiver_arrows(
        ax,
        np.asarray(starts_x, dtype=float),
        np.asarray(starts_y, dtype=float),
        np.asarray(vec_x, dtype=float),
        np.asarray(vec_y, dtype=float),
        color=arrow_color,
        alpha=arrow_alpha,
        width=arrow_width,
    )

def plot_trajectory(samples: TrajectorySamples,
                    title: Optional[str] = None,
                    ax: Optional[plt.Axes] = None,
                    overlay_cmd: Optional[bool] = None,
                    cmd_alpha: float = 0.5,
                    linewidth: float = 1.0,
                    cmd_linewidth: Optional[float] = None,
                    thin_every: int = 1,
                    unwrap_az: bool = False,
                    equal_aspect: bool = False,
                    legend: bool = True,
                    show_points: bool = False,
                    point_size: float = 4.0,
                    show_direction_arrows: bool = False,
                    arrow_mode: str = "segment",
                    arrow_every: int = 50,
                    arrow_alpha: float = 0.25,
                    arrow_color: str = "black",
                    arrow_width: float = 0.0025,
                    arrow_length_frac: float = 0.015,
                    color_by: str = "none",
                    color_vmin: Optional[float] = None,
                    color_vmax: Optional[float] = None,
                    color_percentile: Optional[float] = None,
                    color_percentile_mode: str = "central",
                    color_sign_mode: str = "signed",
                    speed_window_points: int = 2,
                    speed_window_mode: str = "secant") -> Tuple[plt.Figure, plt.Axes, Dict[str, Any]]:
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    else:
        fig = ax.figure

    if cmd_linewidth is None:
        cmd_linewidth = max(0.7, 0.9 * float(linewidth))

    step = max(int(thin_every), 1)
    idx = slice(None, None, step)
    color_metric = str(color_by).strip().lower()

    full_on_run_ids = _compute_on_run_ids(samples.mode_str)
    plot_index = np.arange(int(samples.t_spec.size), dtype=int)[idx]

    x_raw = _unwrap_if_needed(np.asarray(samples.x)[idx], samples.frame_name, samples.x_label.split()[0], unwrap_az)
    y_raw = np.asarray(samples.y)[idx]
    t_raw = np.asarray(samples.t_spec, dtype=float)[idx]
    mode = np.asarray(samples.mode_str, dtype=object)[idx]
    on_run_ids = np.asarray(full_on_run_ids, dtype=int)[idx]
    scan_group_key_raw, scan_group_source = _select_scan_group_key(samples)
    if scan_group_key_raw is not None:
        scan_group_key_raw = np.asarray(scan_group_key_raw, dtype=object)[idx]

    finite = np.isfinite(x_raw) & np.isfinite(y_raw) & np.isfinite(t_raw)

    cmd_x_raw = None
    cmd_y_raw = None
    do_overlay_cmd = bool(samples.cmd_x is not None and samples.cmd_y is not None) if overlay_cmd is None else bool(overlay_cmd)
    need_cmd_metric = color_metric in ("frame_dx", "frame_dy", "frame_dr")
    if samples.cmd_x is not None and samples.cmd_y is not None and (do_overlay_cmd or need_cmd_metric):
        cmd_x_raw = _unwrap_if_needed(np.asarray(samples.cmd_x)[idx], samples.frame_name, samples.x_label.split()[0], unwrap_az)
        cmd_y_raw = np.asarray(samples.cmd_y)[idx]
        finite = finite & np.isfinite(cmd_x_raw) & np.isfinite(cmd_y_raw)

    x_raw = x_raw[finite]
    y_raw = y_raw[finite]
    t_raw = t_raw[finite]
    mode = mode[finite]
    on_run_ids = on_run_ids[finite]
    plot_index = plot_index[finite]
    if cmd_x_raw is not None and cmd_y_raw is not None:
        cmd_x_raw = cmd_x_raw[finite]
        cmd_y_raw = cmd_y_raw[finite]
    if scan_group_key_raw is not None:
        scan_group_key_raw = scan_group_key_raw[finite]

    axis_name0 = samples.x_label.split()[0].strip().lower()
    longitude_continuous = bool(unwrap_az) and (str(samples.frame_name).strip().lower() == "azel") and (axis_name0 == "az")

    x, y, display_meta = _prepare_display_xy(
        x_raw, y_raw, samples.frame_name,
        equal_aspect=bool(equal_aspect),
        longitude_continuous=longitude_continuous,
    )

    drew_direction_arrows = False
    metric_label = None
    metric_signed = False
    n_metric_segments = 0

    if color_metric == "none":
        segments = _make_segments(x, y)
        masks = _make_mode_masks(mode)

        for key in ("ON", "OFF", "HOT", "OTHER"):
            mask = masks[key]
            if mask.size == 0 or not np.any(mask):
                continue
            coll = LineCollection(
                segments[mask],
                colors=MODE_COLORS[key],
                linewidths=float(linewidth),
                alpha=1.0,
                label=key,
            )
            ax.add_collection(coll)
    else:
        metric_segments, metric_values, metric_label, metric_signed = _compute_colored_on_segments(
            x=x_raw,
            y=y_raw,
            t=t_raw,
            on_run_ids=on_run_ids,
            frame_name=samples.frame_name,
            color_by=color_metric,
            cmd_x=cmd_x_raw,
            cmd_y=cmd_y_raw,
            scan_group_key=scan_group_key_raw,
            color_sign_mode=color_sign_mode,
            speed_window_points=speed_window_points,
            speed_window_mode=speed_window_mode,
        )
        n_metric_segments = int(metric_segments.shape[0])
        if n_metric_segments > 0:
            metric_segments_disp = metric_segments.copy()
            metric_x_flat, metric_y_flat, _ = _prepare_display_xy(
                metric_segments[:, :, 0].reshape(-1),
                metric_segments[:, :, 1].reshape(-1),
                samples.frame_name,
                equal_aspect=bool(equal_aspect),
                display_lon0_deg=display_meta.get("lon0_deg"),
                longitude_continuous=longitude_continuous,
            )
            metric_segments_disp[:, :, 0] = np.asarray(metric_x_flat, dtype=float).reshape(metric_segments.shape[0], metric_segments.shape[1])
            metric_segments_disp[:, :, 1] = np.asarray(metric_y_flat, dtype=float).reshape(metric_segments.shape[0], metric_segments.shape[1])
            metric_norm = _make_metric_norm(metric_values, signed=bool(metric_signed),
                                            color_vmin=color_vmin, color_vmax=color_vmax,
                                            color_percentile=color_percentile,
                                            color_percentile_mode=color_percentile_mode)
            metric_coll = LineCollection(
                metric_segments_disp,
                linewidths=float(linewidth),
                alpha=1.0,
                norm=metric_norm,
            )
            metric_coll.set_array(np.asarray(metric_values, dtype=float))
            ax.add_collection(metric_coll)
            cbar = fig.colorbar(metric_coll, ax=ax)
            cbar.set_label(metric_label or color_metric)
        else:
            ax.text(0.02, 0.98, f"No valid ON segments for color_by={color_metric}", transform=ax.transAxes,
                    ha="left", va="top")

    if do_overlay_cmd and cmd_x_raw is not None and cmd_y_raw is not None:
        cmd_x, cmd_y, _ = _prepare_display_xy(
            cmd_x_raw, cmd_y_raw, samples.frame_name,
            equal_aspect=bool(equal_aspect),
            display_lon0_deg=display_meta.get("lon0_deg"),
            longitude_continuous=longitude_continuous,
        )
        cmd_segments = _make_segments(cmd_x, cmd_y)
        cmd_masks = _make_mode_masks(mode)
        for key in ("ON", "OFF", "HOT", "OTHER"):
            mask = cmd_masks[key]
            if mask.size == 0 or not np.any(mask):
                continue
            coll = LineCollection(
                cmd_segments[mask],
                colors=MODE_COLORS[key],
                linewidths=float(cmd_linewidth),
                alpha=float(cmd_alpha),
            )
            ax.add_collection(coll)

    if show_points and x.size > 0:
        if color_metric == "none":
            color_arr = [MODE_COLORS[_normalize_segment_mode(m)] for m in mode]
            ax.scatter(x, y, s=float(point_size), c=color_arr, alpha=0.8, linewidths=0)
        else:
            on_mask = np.asarray([_normalize_segment_mode(m) == "ON" for m in mode], dtype=bool)
            if np.any(on_mask):
                ax.scatter(x[on_mask], y[on_mask], s=float(point_size), c="black", alpha=0.4, linewidths=0)

    arrow_mode_norm = str(arrow_mode).strip().lower()
    if bool(show_direction_arrows) and x.size >= 2 and color_metric == "none":
        if arrow_mode_norm == "samples":
            drew_direction_arrows = _draw_direction_arrows_samples(
                ax,
                x,
                y,
                arrow_every=max(int(arrow_every), 1),
                arrow_alpha=float(arrow_alpha),
                arrow_color=str(arrow_color),
                arrow_width=float(arrow_width),
                arrow_length_frac=float(arrow_length_frac),
            )
        else:
            drew_direction_arrows = _draw_direction_arrows_segment(
                ax,
                x,
                y,
                mode,
                arrow_alpha=float(arrow_alpha),
                arrow_color=str(arrow_color),
                arrow_width=float(arrow_width),
                arrow_length_frac=float(arrow_length_frac),
            )

    ax.autoscale()
    if bool(display_meta.get("reverse_x")):
        ax.invert_xaxis()

    if bool(equal_aspect):
        aspect_ratio = float(display_meta.get("aspect_ratio", 1.0) or 1.0)
        if (not np.isfinite(aspect_ratio)) or (aspect_ratio <= 0.0):
            aspect_ratio = 1.0
        ax.set_aspect(aspect_ratio, adjustable="box")
    ax.set_xlabel(samples.x_label)
    ax.set_ylabel(samples.y_label)
    ax.grid(True, alpha=0.3)

    if title is None:
        metric_title = "" if color_metric == "none" else f", color_by={color_metric}"
        if samples.t_spec.size > 0:
            title = (
                f"Trajectory ({samples.frame_name}, source={samples.coord_source}, stream={samples.stream_name}{metric_title})\n"
                f"{_format_utc_range_from_unix(samples.t_spec)}"
            )
        else:
            title = f"Trajectory ({samples.frame_name}, source={samples.coord_source}, stream={samples.stream_name}{metric_title})"
    ax.set_title(title)

    if legend:
        handles = []
        if color_metric == "none":
            handles.extend([
                Line2D([0], [0], color=MODE_COLORS["ON"], lw=2, label="ON"),
                Line2D([0], [0], color=MODE_COLORS["OFF"], lw=2, label="OFF"),
                Line2D([0], [0], color=MODE_COLORS["HOT"], lw=2, label="HOT"),
                Line2D([0], [0], color=MODE_COLORS["OTHER"], lw=2, label="OTHER"),
            ])
        else:
            handles.append(Line2D([0], [0], color="black", lw=2, label=f"ON colored ({color_metric})"))
        if do_overlay_cmd and cmd_x_raw is not None and cmd_y_raw is not None:
            handles.append(Line2D([0], [0], color="black", lw=2, alpha=float(cmd_alpha), label=f"cmd overlay ({samples.cmd_source})"))
        if drew_direction_arrows:
            handles.append(Line2D([0], [0], color=str(arrow_color), lw=1.2, alpha=float(arrow_alpha), label="direction"))
        if handles:
            ax.legend(handles=handles, loc="best")

    meta = {
        "stream_name": samples.stream_name,
        "spec_table_name": samples.spec_table_name,
        "spec_time_field": samples.spec_time_field,
        "spec_time_basis": samples.spec_time_basis,
        "spec_time_suffix": samples.spec_time_suffix,
        "coord_source": samples.coord_source,
        "coord_meaning": samples.coord_meaning,
        "frame_name": samples.frame_name,
        "n_samples": int(samples.t_spec.size),
        "off_blocks": list(samples.off_blocks or []),
        "cmd_source": samples.cmd_source,
        "cmd_meaning": samples.cmd_meaning,
        "direction_arrows": bool(drew_direction_arrows),
        "arrow_mode": arrow_mode_norm,
        "arrow_every": max(int(arrow_every), 1),
        "arrow_length_frac": float(arrow_length_frac),
        "reverse_x": bool(display_meta.get("reverse_x")),
        "equal_aspect": bool(equal_aspect),
        "display_lat0_deg": display_meta.get("lat0_deg"),
        "display_lon0_deg": display_meta.get("lon0_deg"),
        "display_cos_lat0": display_meta.get("cos_lat0"),
        "color_by": color_metric,
        "metric_label": metric_label,
        "metric_signed": bool(metric_signed),
        "n_metric_segments": int(n_metric_segments),
        "scan_group_source": scan_group_source,
        "scan_field_name": samples.scan_field_name,
        "line_index_field_name": samples.line_index_field_name,
        "line_label_field_name": samples.line_label_field_name,
    }
    return fig, ax, meta



# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def parse_args(argv: Optional[Sequence[str]] = None):
    p = argparse.ArgumentParser(
        prog=pathlib.Path(__file__).name,
        description="Plot NECST v4 telescope trajectory from RawData using a representative spectral stream for timing.",
    )
    p.add_argument("rawdata", nargs="?", default=None, help="RawData folder path (necstdb + nercst)")
    p.add_argument("--spectrometer-config", default=None, help="TOML file describing spectrometers / beams / LO / WCS")
    p.add_argument("--from-csv", default=None, help="Read previously exported trajectory CSV instead of querying RawData DB")
    p.add_argument("--stream-name", dest="stream_names", action="append", default=None,
                   help="Representative stream name used for timing. When omitted, the first enabled/use_for_convert stream is used.")
    p.add_argument("--telescope", default="OMU1P85M", help="Telescope name used in table names")
    p.add_argument("--db-namespace", default=None, help="Database namespace/prefix used in table names")
    p.add_argument("--site-lat", type=float, default=None, help="Site latitude [deg]. Default: auto from RawData config.")
    p.add_argument("--site-lon", type=float, default=None, help="Site longitude [deg]. Default: auto from RawData config.")
    p.add_argument("--site-elev", type=float, default=None, help="Site elevation [m]. Default: auto from RawData config.")
    p.add_argument("--coord-frame", default="azel", choices=["azel", "radec", "galactic"], help="Coordinate frame for plotting")
    p.add_argument("--coord-source", default="beam",
                   choices=["beam", "boresight", "raw_encoder", "raw_altaz", "corrected_cmd", "true", "encoder", "altaz", "cmd"],
                   help="Az/El source before frame conversion. Default: encoder-based beam center.")
    p.add_argument("--overlay-cmd", action="store_true", help="Overlay command trajectory")
    p.add_argument("--cmd-source", default="corrected_cmd",
                   choices=["beam", "boresight", "raw_encoder", "raw_altaz", "corrected_cmd", "true", "encoder", "altaz", "cmd"],
                   help="Command-like source used for overlay")
    p.add_argument("--cmd-alpha", type=float, default=0.5, help="Alpha for command overlay")
    p.add_argument("--azel-correction-apply", "--boresight-correction-apply", dest="azel_correction_apply",
                   default=None, choices=["none", "subtract", "add"],
                   help="How dlon/dlat are applied when forming corrected boresight")
    p.add_argument("--encoder-time-col", default="time", help="Encoder time column (recommended: time)")
    p.add_argument("--altaz-time-col", default="time", help="Altaz time column (recommended: time)")
    p.add_argument("--encoder-shift-sec", type=float, default=0.0, help="Shift applied to encoder timestamps before interpolation")
    p.add_argument("--interp-extrap", default="hold", choices=["nan", "hold"], help="Extrapolation policy for interpolation")
    p.add_argument("--selection", default="all", choices=["all", "from-first-off", "off-block-range"],
                   help="Which part of the trajectory to plot")
    p.add_argument("--off-block-start", type=int, default=None, help="1-based OFF block index used with selection=off-block-range")
    p.add_argument("--off-block-end", type=int, default=None, help="1-based OFF block end index used with selection=off-block-range")
    p.add_argument("--unwrap-az", action="store_true", help="Unwrap Az when coord-frame=azel")
    p.add_argument("--equal-aspect", action="store_true", help="Keep the original coordinate values on the axes, but match the local angular scales by applying a cos(lat0) correction to the x/y display scale. RA/Dec and Galactic are also shown with the conventional left-right reversal.")
    p.add_argument("--thin-every", type=int, default=1, help="Plot every Nth sample")
    p.add_argument("--linewidth", type=float, default=1.0, help="Main trajectory linewidth")
    p.add_argument("--show-points", action="store_true", help="Overlay sample points")
    p.add_argument("--point-size", type=float, default=4.0, help="Point size for --show-points")
    p.add_argument("--show-direction-arrows", action="store_true", help="Overlay thin arrows showing the main-trajectory direction")
    p.add_argument("--arrow-mode", default="segment", choices=["segment", "samples"], help="Direction-arrow placement: segment=one per moving mode segment and major mode transition; samples=every N samples")
    p.add_argument("--arrow-every", type=int, default=50, help="Draw one direction arrow every N plotted samples when --arrow-mode=samples")
    p.add_argument("--arrow-alpha", type=float, default=0.25, help="Alpha for direction arrows")
    p.add_argument("--arrow-width", type=float, default=0.0025, help="Quiver width for direction arrows (axes fraction)")
    p.add_argument("--arrow-length-frac", type=float, default=0.015, help="Arrow length as a fraction of the larger plot span; arrows show direction, not speed")
    p.add_argument("--color-by", default="none", choices=["none", "frame_dx", "frame_dy", "frame_dr", "frame_vx", "frame_vy", "frame_speed"],
                   help="Color ON scan segments by frame-space enc-cmd difference or by frame-space encoder speed. Transitions between ON scans are excluded.")
    p.add_argument("--color-vmin", type=float, default=None,
                   help="Optional lower bound for the color range. For signed metrics, specifying only one of --color-vmin/--color-vmax makes the range symmetric about zero using the absolute value of the specified bound.")
    p.add_argument("--color-vmax", type=float, default=None,
                   help="Optional upper bound for the color range. For signed metrics, specifying only one of --color-vmin/--color-vmax makes the range symmetric about zero using the absolute value of the specified bound.")
    p.add_argument("--color-percentile", type=float, default=None,
                   help="Optional percentile-based auto-range for --color-by. For signed metrics, q means ±Pq(|value|). For unsigned metrics, the mode is controlled by --color-percentile-mode. Explicit --color-vmin/--color-vmax override the percentile side(s) they specify.")
    p.add_argument("--color-percentile-mode", default="central", choices=["central", "zero_upper"],
                   help="How --color-percentile is interpreted for unsigned metrics: central keeps the central q percent of values using [P((100-q)/2), P(100-(100-q)/2)], while zero_upper uses [0, Pq]. Ignored for signed metrics.")
    p.add_argument("--color-sign-mode", default="signed", choices=["signed", "abs"],
                   help="For signed metrics (frame_dx, frame_dy, frame_vx, frame_vy), choose signed values or their absolute values. Unsigned metrics are unchanged.")
    p.add_argument("--speed-window-points", type=int, default=2,
                   help="Number of points used for velocity metrics (frame_vx, frame_vy, frame_speed). 2 means per-segment instantaneous velocity; values >2 use a same-ON same-scan local window.")
    p.add_argument("--speed-window-mode", default="secant", choices=["secant", "component_mean", "magnitude_mean"],
                   help="How to average velocity metrics when --speed-window-points > 2: secant uses the window endpoints, component_mean averages vx and vy separately, magnitude_mean averages the speed magnitude.")
    p.add_argument("--hide-legend", action="store_true", help="Do not show legend")
    p.add_argument("--title", default=None, help="Optional plot title")
    p.add_argument("--out", default=None, help="Output image path (png/pdf/svg). Default: rawbasename_trajectory_<frame>.png")
    p.add_argument("--export-csv", nargs="?", const="", default=None,
                   help="Write trajectory parameters to CSV. When no path is given, use rawbasename_trajectory.csv or csvbasename_plot.csv")
    p.add_argument("--dpi", type=float, default=150.0, help="Savefig DPI")
    p.add_argument("--show", action="store_true", help="Show plot interactively")
    p.add_argument("--converter-path", default=None, help="Path to necst_v4_sdfits_converter.py when it is not adjacent to this script")
    p.add_argument("--print-mode-summary", action="store_true", help="Print all encountered OBSMODE labels and counts")
    return p.parse_args(argv)


def _default_output_path(rawdata_path: pathlib.Path, coord_frame: str) -> pathlib.Path:
    stem = rawdata_path.name if rawdata_path.is_dir() else rawdata_path.stem
    return rawdata_path.with_name(f"{stem}_trajectory_{coord_frame}.png")


def _default_csv_path(raw_path: pathlib.Path, from_csv: bool = False) -> pathlib.Path:
    stem = raw_path.stem if from_csv else raw_path.name
    suffix = "_plot.csv" if from_csv else "_trajectory.csv"
    return raw_path.with_name(f"{stem}{suffix}")


def _print_mode_summary(mode_str: Sequence[Any]):
    uniq, cnt = np.unique(np.asarray(mode_str, dtype=object), return_counts=True)
    print("[info] mode summary:")
    for u, c in zip(uniq.tolist(), cnt.tolist()):
        print(f"  {u}: {c}")


def _argv_has_option_local(argv_list: Sequence[str], opt: str) -> bool:
    opt = str(opt)
    prefix = opt + "="
    for a in argv_list:
        s = str(a)
        if s == opt or s.startswith(prefix):
            return True
    return False


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv_list = list(sys.argv[1:] if argv is None else argv)
    args = parse_args(argv_list)

    if (args.rawdata is None) and (args.from_csv is None):
        raise SystemExit("either RAWDATA or --from-csv is required")
    if (args.rawdata is not None) and (args.from_csv is not None):
        raise SystemExit("use either RAWDATA or --from-csv, not both")

    need_cmd_frame = bool(args.overlay_cmd) or str(args.color_by).strip().lower() in ("frame_dx", "frame_dy", "frame_dr")

    if args.from_csv is not None:
        samples = load_trajectory_samples_from_csv(
            csv_path=args.from_csv,
            coord_frame=args.coord_frame,
            coord_source=args.coord_source,
            overlay_cmd=bool(need_cmd_frame),
            cmd_source=args.cmd_source,
            selection=args.selection,
            off_block_start=args.off_block_start,
            off_block_end=args.off_block_end,
        )
        base_path = pathlib.Path(args.from_csv).expanduser().resolve()
    else:
        samples = load_trajectory_samples(
            rawdata=args.rawdata,
            spectrometer_config=args.spectrometer_config,
            stream_names=args.stream_names,
            telescope=args.telescope,
            db_namespace=args.db_namespace,
            coord_frame=args.coord_frame,
            coord_source=args.coord_source,
            overlay_cmd=bool(need_cmd_frame),
            cmd_source=args.cmd_source,
            azel_correction_apply=args.azel_correction_apply,
            encoder_time_col=args.encoder_time_col,
            altaz_time_col=args.altaz_time_col,
            encoder_shift_sec=float(args.encoder_shift_sec),
            interp_extrap=args.interp_extrap,
            site_lat=args.site_lat,
            site_lon=args.site_lon,
            site_elev=args.site_elev,
            converter_path=args.converter_path,
            selection=args.selection,
            off_block_start=args.off_block_start,
            off_block_end=args.off_block_end,
            telescope_explicit=_argv_has_option_local(argv_list, "--telescope"),
            db_namespace_explicit=_argv_has_option_local(argv_list, "--db-namespace"),
            encoder_time_col_explicit=_argv_has_option_local(argv_list, "--encoder-time-col"),
            altaz_time_col_explicit=_argv_has_option_local(argv_list, "--altaz-time-col"),
            encoder_shift_sec_explicit=_argv_has_option_local(argv_list, "--encoder-shift-sec"),
        )
        base_path = pathlib.Path(args.rawdata).expanduser().resolve()

    if args.print_mode_summary:
        _print_mode_summary(samples.mode_str)
        if samples.off_blocks is not None:
            print(f"[info] OFF blocks (1-based): {[(i + 1, a, b) for i, (a, b) in enumerate(samples.off_blocks)]}")

    if args.export_csv is not None:
        if args.export_csv == "":
            csv_path = _default_csv_path(base_path, from_csv=(args.from_csv is not None))
        else:
            csv_path = pathlib.Path(args.export_csv).expanduser().resolve()
        wrote_csv = export_trajectory_csv(samples, str(csv_path))
        print(f"[info] wrote CSV {wrote_csv}")

    do_show = bool(args.show) or (args.out is None)
    do_save = (args.out is not None)
    if not do_show:
        try:
            plt.switch_backend("Agg")
        except Exception:
            pass

    fig, ax, meta = plot_trajectory(
        samples,
        title=args.title,
        overlay_cmd=bool(args.overlay_cmd),
        cmd_alpha=float(args.cmd_alpha),
        linewidth=float(args.linewidth),
        thin_every=max(int(args.thin_every), 1),
        unwrap_az=bool(args.unwrap_az),
        equal_aspect=bool(args.equal_aspect),
        legend=not bool(args.hide_legend),
        show_points=bool(args.show_points),
        point_size=float(args.point_size),
        show_direction_arrows=bool(args.show_direction_arrows),
        arrow_mode=str(args.arrow_mode),
        arrow_every=max(int(args.arrow_every), 1),
        arrow_alpha=float(args.arrow_alpha),
        arrow_width=float(args.arrow_width),
        arrow_length_frac=float(args.arrow_length_frac),
        color_by=str(args.color_by),
        color_vmin=args.color_vmin,
        color_vmax=args.color_vmax,
        color_percentile=args.color_percentile,
        color_percentile_mode=str(args.color_percentile_mode),
        color_sign_mode=str(args.color_sign_mode),
        speed_window_points=int(args.speed_window_points),
        speed_window_mode=str(args.speed_window_mode),
    )

    if do_save:
        out_path = pathlib.Path(args.out).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(out_path), dpi=float(args.dpi), bbox_inches="tight")
        print(f"[info] wrote {out_path}")

    print(f"[info] stream={meta['stream_name']} spec_table={meta['spec_table_name']} time_basis={meta['spec_time_basis']}")
    print(f"[info] coord_source={meta['coord_source']} frame={meta['frame_name']} n_samples={meta['n_samples']}")

    if do_show:
        plt.show()
    else:
        plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
