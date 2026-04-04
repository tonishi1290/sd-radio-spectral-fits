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


def read_spectral_timing_only(conv, db, spec_table_name: str):
    """
    Return only spectral timing and OBSMODE-like labels.

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
        return t_spec, mode_str, time_meta, pos_field

    names = list(dtype.names)
    pos_field = conv._pick_field_name(names, None, ["position", "obsmode", "obs_mode", "obsMode", "mode", "state", "status", "label"])
    timestamp_field = conv._pick_field_name(names, None, ["timestamp", "time_spectrometer"])
    fallback_field = conv._pick_field_name(names, None, ["time", "t", "unix_time", "unixtime"])

    needed_cols = []
    for name in (timestamp_field, fallback_field, pos_field):
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
    return t_spec, mode_str, time_meta, pos_field


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

    t_spec_unsorted, mode_unsorted, time_meta, pos_field = read_spectral_timing_only(conv, db, spec_table_name)
    order = np.argsort(np.asarray(t_spec_unsorted, dtype=float))
    t_spec = np.asarray(t_spec_unsorted, dtype=float)[order]
    mode_str = np.asarray(mode_unsorted, dtype=object)[order]

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
        pivot='tail',
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


def _draw_direction_arrows_segment(ax: plt.Axes, x: np.ndarray, y: np.ndarray, mode: Sequence[Any], *,
                                   arrow_alpha: float, arrow_color: str, arrow_width: float,
                                   arrow_length_frac: float) -> bool:
    if x.size < 2:
        return False
    span_x = float(np.nanmax(x) - np.nanmin(x)) if x.size else 0.0
    span_y = float(np.nanmax(y) - np.nanmin(y)) if y.size else 0.0
    span = max(span_x, span_y)
    if (not np.isfinite(span)) or (span <= 0.0):
        return False
    draw_len = max(float(arrow_length_frac), 0.0) * span
    if not np.isfinite(draw_len) or draw_len <= 0.0:
        draw_len = 0.015 * span
    move_min = 0.01 * span

    runs = _find_mode_runs(mode)
    starts_x = []
    starts_y = []
    vec_x = []
    vec_y = []

    for s, e, _key in runs:
        if e - s < 2:
            continue
        dx = float(x[e - 1] - x[s])
        dy = float(y[e - 1] - y[s])
        dist = float(np.hypot(dx, dy))
        if (not np.isfinite(dist)) or (dist < move_min):
            continue
        mid = int((s + e - 1) // 2)
        starts_x.append(float(x[mid]))
        starts_y.append(float(y[mid]))
        vec_x.append(draw_len * dx / dist)
        vec_y.append(draw_len * dy / dist)

    for (s0, e0, _k0), (s1, e1, _k1) in zip(runs[:-1], runs[1:]):
        if s1 <= 0 or e0 <= 0:
            continue
        dx = float(x[s1] - x[e0 - 1])
        dy = float(y[s1] - y[e0 - 1])
        dist = float(np.hypot(dx, dy))
        if (not np.isfinite(dist)) or (dist < move_min):
            continue
        x0 = 0.5 * float(x[e0 - 1] + x[s1])
        y0 = 0.5 * float(y[e0 - 1] + y[s1])
        starts_x.append(x0)
        starts_y.append(y0)
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
                    legend: bool = True,
                    show_points: bool = False,
                    point_size: float = 4.0,
                    show_direction_arrows: bool = False,
                    arrow_mode: str = "segment",
                    arrow_every: int = 50,
                    arrow_alpha: float = 0.25,
                    arrow_color: str = "black",
                    arrow_width: float = 0.0025,
                    arrow_length_frac: float = 0.015) -> Tuple[plt.Figure, plt.Axes, Dict[str, Any]]:
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    else:
        fig = ax.figure

    if cmd_linewidth is None:
        cmd_linewidth = max(0.7, 0.9 * float(linewidth))

    step = max(int(thin_every), 1)
    idx = slice(None, None, step)

    x = _unwrap_if_needed(np.asarray(samples.x)[idx], samples.frame_name, samples.x_label.split()[0], unwrap_az)
    y = np.asarray(samples.y)[idx]
    mode = np.asarray(samples.mode_str, dtype=object)[idx]

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    mode = mode[finite]

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

    do_overlay_cmd = bool(samples.cmd_x is not None and samples.cmd_y is not None) if overlay_cmd is None else bool(overlay_cmd)
    if do_overlay_cmd:
        cmd_x = _unwrap_if_needed(np.asarray(samples.cmd_x)[idx], samples.frame_name, samples.x_label.split()[0], unwrap_az)
        cmd_y = np.asarray(samples.cmd_y)[idx]
        cmd_finite = np.isfinite(cmd_x) & np.isfinite(cmd_y)
        cmd_x = cmd_x[cmd_finite]
        cmd_y = cmd_y[cmd_finite]
        cmd_mode = np.asarray(samples.mode_str, dtype=object)[idx][cmd_finite]
        cmd_segments = _make_segments(cmd_x, cmd_y)
        cmd_masks = _make_mode_masks(cmd_mode)
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
        color_arr = [MODE_COLORS[_normalize_segment_mode(m)] for m in mode]
        ax.scatter(x, y, s=float(point_size), c=color_arr, alpha=0.8, linewidths=0)

    drew_direction_arrows = False
    arrow_mode_norm = str(arrow_mode).strip().lower()
    if bool(show_direction_arrows) and x.size >= 2:
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
    ax.set_xlabel(samples.x_label)
    ax.set_ylabel(samples.y_label)
    ax.grid(True, alpha=0.3)

    if title is None:
        if samples.t_spec.size > 0:
            t0 = float(np.nanmin(samples.t_spec))
            t1 = float(np.nanmax(samples.t_spec))
            title = (
                f"Trajectory ({samples.frame_name}, source={samples.coord_source}, stream={samples.stream_name})\n"
                f"unix={t0:.3f} .. {t1:.3f}"
            )
        else:
            title = f"Trajectory ({samples.frame_name}, source={samples.coord_source}, stream={samples.stream_name})"
    ax.set_title(title)

    if legend:
        handles = [
            Line2D([0], [0], color=MODE_COLORS["ON"], lw=2, label="ON"),
            Line2D([0], [0], color=MODE_COLORS["OFF"], lw=2, label="OFF"),
            Line2D([0], [0], color=MODE_COLORS["HOT"], lw=2, label="HOT"),
            Line2D([0], [0], color=MODE_COLORS["OTHER"], lw=2, label="OTHER"),
        ]
        if do_overlay_cmd:
            handles.append(Line2D([0], [0], color="black", lw=2, alpha=float(cmd_alpha), label=f"cmd overlay ({samples.cmd_source})"))
        if drew_direction_arrows:
            handles.append(Line2D([0], [0], color=str(arrow_color), lw=1.2, alpha=float(arrow_alpha), label="direction"))
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

    if args.from_csv is not None:
        samples = load_trajectory_samples_from_csv(
            csv_path=args.from_csv,
            coord_frame=args.coord_frame,
            coord_source=args.coord_source,
            overlay_cmd=bool(args.overlay_cmd),
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
            overlay_cmd=bool(args.overlay_cmd),
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
        legend=not bool(args.hide_legend),
        show_points=bool(args.show_points),
        point_size=float(args.point_size),
        show_direction_arrows=bool(args.show_direction_arrows),
        arrow_mode=str(args.arrow_mode),
        arrow_every=max(int(args.arrow_every), 1),
        arrow_alpha=float(args.arrow_alpha),
        arrow_width=float(args.arrow_width),
        arrow_length_frac=float(args.arrow_length_frac),
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
