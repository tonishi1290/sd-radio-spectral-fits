# src/sd_radio_spectral_fits/sdfits_bintable.py
from __future__ import annotations

"""
sdfits_bintable.py
==================

Single Source of Truth for building the SDFITS (SINGLE DISH) BinTableHDU.

Design goals
------------
- **Permissive**: accept extra analysis columns (including vector columns stored as
  "array-in-a-cell" in a pandas DataFrame) and write them back to FITS.
- **Robust**: automatically select fixed-length vectors vs variable-length arrays (VLA: 'P*'/'Q*')
  based on per-row lengths.
- **Minimal policy**: this module does *not* decide how to "promote" header keywords into
  per-row columns. That policy belongs to wrappers (e.g., fitsio.write_scantable()).

This module is used by:
- fitsio.py (analysis I/O; DataFrame-driven)
- sdfits_writer.py (observing/convert writer; row-buffer-driven)
"""

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import re
import warnings

from astropy.io import fits

# -----------------------------------------------------------------------------
# FITS keyword helpers (HIERARCH unification)
# -----------------------------------------------------------------------------

_FITS_STD_KEY_RE = re.compile(r"^[A-Z0-9_-]{1,8}$")


def set_meta_keyword(
    hdr: "fits.Header",
    key: str,
    value: Any,
    comment: str | None = None,
    *,
    prefer_hierarch_for_nonstandard: bool = True,
) -> None:
    """Set a FITS header keyword robustly, unifying HIERARCH handling.

    Rules (minimal, predictable):
    - If `key` already starts with 'HIERARCH ' (case-insensitive), it is used as-is.
    - If `key` is a standard FITS keyword (<=8 chars, A-Z0-9_-), set directly.
    - Otherwise, if `prefer_hierarch_for_nonstandard=True`, write as 'HIERARCH <key>'.
      This keeps long or nonstandard keys round-trippable.

    `value` may be a numpy scalar; it will be converted to a plain Python scalar when possible.
    """
    if key is None:
        return
    k = str(key).strip()
    if not k:
        return
    ku = k.upper()

    # normalize numpy scalar
    try:
        if isinstance(value, np.generic):
            value = value.item()
    except Exception:
        pass

    # skip reserved commentary cards
    if ku in ("COMMENT", "HISTORY"):
        return

    out_key = k
    if ku.startswith("HIERARCH "):
        out_key = k  # already explicit
    elif _FITS_STD_KEY_RE.match(ku) is not None:
        out_key = ku
    else:
        out_key = f"HIERARCH {k}" if prefer_hierarch_for_nonstandard else k

    if comment is None:
        hdr[out_key] = value
    else:
        hdr[out_key] = (value, str(comment))


def apply_meta_to_header(
    hdr: "fits.Header",
    meta: Mapping[str, Any],
    *,
    prefer_hierarch_for_nonstandard: bool = True,
) -> Dict[str, str]:
    """Apply `meta` dict to header and return a dict of failures (key -> reason).

    This is intentionally permissive: failures are returned (not raised) so callers may
    decide how to record them (e.g., into HISTORY).
    """
    failures: Dict[str, str] = {}
    for k, v in (meta or {}).items():
        try:
            if isinstance(v, tuple) and len(v) >= 2:
                val, com = v[0], v[1]
            else:
                val, com = v, None
            set_meta_keyword(
                hdr,
                str(k),
                val,
                comment=(str(com) if com is not None else None),
                prefer_hierarch_for_nonstandard=prefer_hierarch_for_nonstandard,
            )
        except Exception as e:
            failures[str(k)] = str(e)
    return failures



# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _is_null(x: Any) -> bool:
    if x is None:
        return True
    try:
        # numpy/pandas NaN
        return bool(pd.isna(x))
    except Exception:
        return False


def _first_nonnull(values: Iterable[Any]) -> Any:
    for v in values:
        if not _is_null(v):
            return v
    return None


def _dedup_case_insensitive_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate columns that would collide in FITS after uppercasing.

    - FITS column names are case-insensitive.
    - If duplicates exist, keep the earliest occurrence (stable), unless a preferred
      name is present (e.g., keep 'TIMESTAMP' over 'timestamp').
    """
    cols = list(df.columns)
    prefer_keep = {"TIMESTAMP": "TIMESTAMP"}

    drop = set()
    up_to_originals: Dict[str, List[str]] = {}
    for c in cols:
        up = str(c).upper()
        up_to_originals.setdefault(up, []).append(c)

    for up, keep_name in prefer_keep.items():
        if up in up_to_originals and keep_name in up_to_originals[up]:
            for c in up_to_originals[up]:
                if c != keep_name:
                    drop.add(c)

    ordered_keep: List[str] = []
    kept_upper = set()
    for c in cols:
        if c in drop:
            continue
        up = str(c).upper()
        if up in kept_upper:
            continue
        kept_upper.add(up)
        ordered_keep.append(c)

    return df.loc[:, ordered_keep].copy()


def _normalize_columns_upper(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(c).strip().upper() for c in out.columns]
    out = _dedup_case_insensitive_columns(out)
    # If exact duplicates remain, keep last (pandas behavior previously used in fitsio.py)
    if out.columns.duplicated().any():
        out = out.loc[:, ~out.columns.duplicated(keep="last")]
    return out


def _dtype_to_scalar_format(arr: np.ndarray) -> Tuple[str, np.ndarray]:
    """
    Map numpy dtype to FITS scalar column format and a clean array.
    Uses conservative, widely-supported formats.
    """
    a = np.asarray(arr)
    if a.dtype.kind in ("S", "U", "O"):
        # strings: handled elsewhere
        raise TypeError("string dtype should be handled by _make_string_column")
    if a.dtype.kind in ("b",):
        # preserve historical behavior: store as int8 'B' (not logical 'L')
        return "B", a.astype(np.int8, copy=False)
    if a.dtype.kind in ("i", "u"):
        # choose by itemsize
        if a.dtype.itemsize <= 2:
            return "I", a.astype(np.int16, copy=False)
        if a.dtype.itemsize <= 4:
            return "J", a.astype(np.int32, copy=False)
        return "K", a.astype(np.int64, copy=False)
    if a.dtype.kind == "f":
        if a.dtype.itemsize <= 4:
            return "E", a.astype(np.float32, copy=False)
        return "D", a.astype(np.float64, copy=False)
    # fallback: stringify
    raise TypeError(f"unsupported dtype for scalar: {a.dtype}")


def _make_string_column(name: str, values: Sequence[Any], width: Optional[int] = None) -> fits.Column:
    s = ["" if _is_null(x) else str(x) for x in values]
    if width is None:
        width = int(max(1, max((len(v) for v in s), default=1)))
    # Use bytes storage for FITS ASCII columns
    arr = np.asarray(s, dtype=f"S{width}")
    return fits.Column(name=name, array=arr, format=f"{width}A")


def _as_1d_array(x: Any, dtype: Optional[np.dtype] = None) -> Optional[np.ndarray]:
    if _is_null(x):
        return None
    a = np.asarray(x)
    if a.ndim == 0:
        # scalar; treat as length-1 vector (rare)
        a = a.reshape(1)
    else:
        a = a.reshape(-1)
    if dtype is not None:
        a = a.astype(dtype, copy=False)
    return a


@dataclass(frozen=True)
class VectorColumnSpec:
    name: str
    format: str
    array: np.ndarray  # 2D fixed or 1D object array for VLA
    unit: Optional[str] = None


def _infer_vector_column(
    name: str,
    series: pd.Series,
    *,
    prefer_float64: bool = True,
    allow_vla: bool = True,
) -> Optional[VectorColumnSpec]:
    """
    Infer a FITS vector column from a Series containing array-like entries.

    - If all non-null rows have the same length and no nulls -> fixed-length 'nD'/'nE'/'nJ' etc
    - Otherwise -> VLA 'P*' (object array) if allow_vla, else return None.
    """
    vals = list(series.values)
    sample = _first_nonnull(vals)
    if sample is None:
        return None

    # Determine base dtype from sample
    samp_arr = _as_1d_array(sample)
    if samp_arr is None:
        return None

    # Choose dtype/format letter
    if samp_arr.dtype.kind == "f":
        if prefer_float64:
            base_dtype = np.float64
            letter = "D"
        else:
            base_dtype = np.float32
            letter = "E"
    elif samp_arr.dtype.kind in ("i", "u"):
        base_dtype = np.int64
        letter = "K"
    elif samp_arr.dtype.kind == "b":
        base_dtype = np.bool_
        letter = "L"
    else:
        # vector of strings or unknown -> not supported here
        return None

    arrays: List[Optional[np.ndarray]] = []
    lengths: List[int] = []
    has_null = False
    for v in vals:
        a = _as_1d_array(v, dtype=base_dtype)
        if a is None:
            has_null = True
            arrays.append(None)
            lengths.append(0)
        else:
            arrays.append(a)
            lengths.append(int(a.size))

    # If any nulls exist, force VLA (to avoid inventing fill values)
    uniform = (len(set(lengths)) <= 1) and (not has_null)
    n = int(lengths[0]) if lengths else 0

    if uniform and n > 0:
        mat = np.vstack([a for a in arrays if a is not None])
        fmt = f"{n}{letter}"
        return VectorColumnSpec(name=name, format=fmt, array=mat)
    if not allow_vla:
        return None

    # VLA: represent null as empty vector
    obj = np.empty(len(vals), dtype=object)
    for i, a in enumerate(arrays):
        if a is None:
            obj[i] = np.asarray([], dtype=base_dtype)
        else:
            obj[i] = a
    fmt = f"P{letter}"
    return VectorColumnSpec(name=name, format=fmt, array=obj)


def _spectrum_and_flag_columns(
    data: Union[np.ndarray, List[np.ndarray]],
    *,
    spectrum_column: str = "DATA",
    include_flag: bool = True,
    flag: Optional[Union[np.ndarray, List[np.ndarray]]] = None,
) -> Tuple[List[fits.Column], int, bool]:
    """
    Build spectrum (and optional flag) columns and return:
      columns, nchan_hint, is_vla
    """
    # Normalize data into either 2D float32 (fixed) or object array (VLA)
    is_vla = False
    nchan_hint = 0

    if isinstance(data, np.ndarray) and (data.dtype != object):
        d2 = np.asarray(data)
        if d2.ndim != 2:
            raise ValueError(f"data must be 2D (ndumps, nchan), got shape={d2.shape}")
        nd, nchan = d2.shape
        nchan_hint = int(nchan)
        cols: List[fits.Column] = [
            fits.Column(name=spectrum_column, array=d2.astype(np.float32, copy=False), format=f"{nchan}E")
        ]
        if include_flag:
            if flag is None:
                f2 = np.zeros((nd, nchan), dtype=bool)
            else:
                f2 = np.asarray(flag, dtype=bool)
                if f2.shape != (nd, nchan):
                    raise ValueError(f"FLAG shape {f2.shape} != DATA shape {(nd, nchan)}")
            cols.append(fits.Column(name="FLAG", array=f2, format=f"{nchan}L"))
        return cols, nchan_hint, False

    # list or object array
    if isinstance(data, np.ndarray) and data.dtype == object:
        rows = list(data)
    else:
        rows = list(data)

    lengths = [len(np.asarray(x).reshape(-1)) for x in rows]
    if len(lengths) == 0:
        raise ValueError("empty data (no rows)")

    nchan_max = int(max(lengths)) if lengths else 0
    nchan_hint = nchan_max

    if len(set(lengths)) == 1:
        # fixed, stack
        nchan = int(lengths[0])
        d2 = np.vstack([np.asarray(x, dtype=np.float32).reshape(-1) for x in rows])
        nd = d2.shape[0]
        cols = [fits.Column(name=spectrum_column, array=d2, format=f"{nchan}E")]
        if include_flag:
            if flag is None:
                f2 = np.zeros((nd, nchan), dtype=bool)
            else:
                if isinstance(flag, np.ndarray) and flag.dtype != object:
                    f2 = np.asarray(flag, dtype=bool)
                else:
                    f2 = np.vstack([np.asarray(x, dtype=bool).reshape(-1) for x in list(flag)])  # type: ignore[arg-type]
                if f2.shape != (nd, nchan):
                    raise ValueError(f"FLAG shape {f2.shape} != DATA shape {(nd, nchan)}")
            cols.append(fits.Column(name="FLAG", array=f2, format=f"{nchan}L"))
        return cols, nchan, False

    # VLA
    is_vla = True
    nd = len(rows)

    data_obj = np.empty(nd, dtype=object)
    data_obj[:] = [np.asarray(x, dtype=np.float32).reshape(-1) for x in rows]
    cols = [fits.Column(name=spectrum_column, array=data_obj, format="PE")]

    if include_flag:
        if flag is None:
            flag_obj = np.empty(nd, dtype=object)
            flag_obj[:] = [np.zeros(len(np.asarray(x).reshape(-1)), dtype=bool) for x in rows]
        else:
            frows = list(flag) if not (isinstance(flag, np.ndarray) and flag.dtype != object) else list(flag)  # type: ignore[list-item]
            if len(frows) != nd:
                raise ValueError("FLAG rows != DATA rows")
            flag_obj = np.empty(nd, dtype=object)
            for i, fr in enumerate(frows):
                fa = np.asarray(fr, dtype=bool).reshape(-1)
                if fa.size != len(np.asarray(rows[i]).reshape(-1)):
                    raise ValueError(f"FLAG length {fa.size} != DATA length {len(np.asarray(rows[i]).reshape(-1))} at row {i}")
                flag_obj[i] = fa
        cols.append(fits.Column(name="FLAG", array=flag_obj, format="PL"))

    return cols, nchan_hint, True


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# HISTORY extension helper (shared by fitsio.py and sdfits_writer.py)
# -----------------------------------------------------------------------------

def build_history_hdu(history: Any) -> fits.BinTableHDU | None:
    """Build a HISTORY BinTableHDU from dict/list/str-like history.

    Accepted forms:
    - dict: {key: value}
    - list/tuple: sequence of dicts (merged with prefix 'i:key') or scalars
    - other: stringified into key '0'

    Returns None if `history` is falsy or results in an empty table.
    """
    if not history:
        return None

    hist_flat: Dict[str, str] = {}
    if isinstance(history, dict):
        for k, v in history.items():
            hist_flat[str(k)] = str(v)
    elif isinstance(history, (list, tuple)):
        for i, h in enumerate(history):
            if h is None:
                continue
            if not isinstance(h, dict):
                hist_flat[f"{i}"] = str(h)
                continue
            for k, v in h.items():
                hist_flat[f"{i}:{k}"] = str(v)
    else:
        hist_flat["0"] = str(history)

    if not hist_flat:
        return None

    keys = list(hist_flat.keys())
    vals = [str(hist_flat[k]) for k in keys]

    w_k = int(max(1, max((len(str(k)) for k in keys), default=1)))
    w_v = int(max(1, max((len(str(v)) for v in vals), default=1)))

    col_k = fits.Column(name="KEY", format=f"{w_k}A", array=keys)
    col_v = fits.Column(name="VALUE", format=f"{w_v}A", array=vals)
    return fits.BinTableHDU.from_columns([col_k, col_v], name="HISTORY")

def build_single_dish_table_hdu(
    *,
    table: pd.DataFrame,
    data: Union[np.ndarray, List[np.ndarray]],
    spectrum_column: str = "DATA",
    include_flag: bool = True,
    flag: Optional[Union[np.ndarray, List[np.ndarray]]] = None,
    extra_vector_columns: Optional[Mapping[str, Sequence[np.ndarray]]] = None,
    units: Optional[Mapping[str, str]] = None,
    normalize_columns: bool = True,
    string_widths: Optional[Mapping[str, int]] = None,
    prefer_float64_vectors: bool = True,
    allow_vla_vectors: bool = True,
    bunit: Optional[str] = "K",
    warn_dropped_vector_columns: bool = False,
    extname: str = "SINGLE DISH",
) -> Tuple[fits.BinTableHDU, int]:
    """
    Build the SINGLE DISH BinTableHDU from a DataFrame and spectrum data.

    Parameters
    ----------
    table:
        Per-row metadata table (pandas DataFrame). Columns may include scalar values or
        vector-in-cell values (list/tuple/np.ndarray).
    data:
        Spectrum array per row. Either:
          - 2D ndarray: (nrow, nchan) fixed-length
          - list of 1D arrays: ragged OK (VLA written if needed)
    spectrum_column:
        FITS column name for the spectrum (default 'DATA').
    include_flag / flag:
        If include_flag=True, write a FLAG column. If flag is None, FLAG is created as zeros.
        If provided, flag must match DATA shape/lengths.
    extra_vector_columns:
        Additional vector columns not present in the DataFrame, e.g. {'FREQ': list_of_freq_arrays}.
        This is useful for row-buffer writers.
    normalize_columns:
        If True, columns are uppercased and deduplicated in FITS-safe manner.
    string_widths:
        Optional fixed widths for specific string columns (keys are case-insensitive).
    prefer_float64_vectors:
        If True, vector columns inferred as floats are written as float64 ('D'). Otherwise float32 ('E').
    allow_vla_vectors:
        If True, vector columns with ragged rows are written as VLA ('P*').
    bunit:
        Set table header BUNIT (and/or spectrum unit). Default 'K'.
    extname:
        Extension name for the table (default 'SINGLE DISH').
    warn_dropped_vector_columns:
        If True, emit a RuntimeWarning when a column is dropped because its values are all null/NaN
        (i.e., the format cannot be inferred).


    Returns
    -------
    (hdu, nchan_hint)
        hdu is the BinTableHDU.
        nchan_hint is the fixed channel count if fixed, else the maximum length across rows.
    """
    if normalize_columns:
        df = _normalize_columns_upper(table)
    else:
        df = table.copy()

    widths = {str(k).upper(): int(v) for k, v in (string_widths or {}).items()}

    reserved = {spectrum_column.upper(), "SPECTRUM", "FLAG", "FLAGS"}

    columns: List[fits.Column] = []
    dropped_cols: List[str] = []  # columns with unknown type (all-null), dropped

    # 1) metadata columns from DataFrame
    for name in df.columns:
        uname = str(name).strip().upper()
        if uname in reserved:
            continue

        ser = df[name]

        sample = _first_nonnull(ser.values)
        if sample is None:
            # no type info -> drop (all-null / all-NaN optional column)
            dropped_cols.append(uname)
            continue

        # vector-in-cell?
        if isinstance(sample, (list, tuple, np.ndarray)) and not (np.asarray(sample).ndim == 0):
            spec = _infer_vector_column(
                uname,
                ser,
                prefer_float64=prefer_float64_vectors,
                allow_vla=allow_vla_vectors,
            )
            if spec is None:
                # fallback: stringify scalar-ish representation per row
                columns.append(_make_string_column(uname, ser.values, width=widths.get(uname)))
            else:
                col = fits.Column(name=spec.name, format=spec.format, array=spec.array)
                if spec.unit:
                    col.unit = spec.unit
                columns.append(col)
            continue

        # scalar
        try:
            arr = ser.to_numpy()
            if hasattr(arr, "ndim") and arr.ndim > 1:
                # unexpected; skip
                continue

            if np.asarray(arr).dtype.kind in ("S", "U", "O"):
                columns.append(_make_string_column(uname, arr, width=widths.get(uname)))
            else:
                fmt, a2 = _dtype_to_scalar_format(np.asarray(arr))
                columns.append(fits.Column(name=uname, format=fmt, array=a2))
        except Exception:
            # last resort: stringify
            columns.append(_make_string_column(uname, ser.values, width=widths.get(uname)))

    # 2) extra vector columns (e.g. FREQ)
    if extra_vector_columns:
        for k, vlist in extra_vector_columns.items():
            uname = str(k).strip().upper()
            if uname in reserved:
                continue
            ser = pd.Series(list(vlist))
            spec = _infer_vector_column(
                uname,
                ser,
                prefer_float64=prefer_float64_vectors,
                allow_vla=allow_vla_vectors,
            )
            if spec is None:
                # Can't infer; drop silently
                dropped_cols.append(uname)
                continue
            col = fits.Column(name=spec.name, format=spec.format, array=spec.array)
            columns.append(col)

    # 3) spectrum + flag
    spec_cols, nchan_hint, is_vla = _spectrum_and_flag_columns(
        data,
        spectrum_column=spectrum_column.upper(),
        include_flag=include_flag,
        flag=flag,
    )
    # Assign units to spectrum column (and possibly others)
    if bunit:
        for c in spec_cols:
            if c.name.upper() == spectrum_column.upper():
                c.unit = bunit
    columns.extend(spec_cols)

    units_map = {str(k).upper(): str(v) for k, v in (units or {}).items()}
    if units_map:
        for c in columns:
            ukey = str(c.name).upper()
            if ukey in units_map and getattr(c, 'unit', None) in (None, ''):
                c.unit = units_map[ukey]

    hdu = fits.BinTableHDU.from_columns(columns, name=extname)
    hdu.name = extname

    # Header cosmetics
    eh = hdu.header
    eh["EXTNAME"] = (extname, "Single dish table (SDFITS-like)")
    if bunit:
        eh["BUNIT"] = (str(bunit), "Brightness unit for spectrum")

    # If fixed-length, write TDIM for commonly-used vector columns
    if not is_vla:
        # Find nchan for DATA/FLAG/FREQ etc
        for i in range(1, int(eh.get("TFIELDS", 0)) + 1):
            ttype = (eh.get(f"TTYPE{i}", "") or "").strip().upper()
            tform = (eh.get(f"TFORM{i}", "") or "").strip().upper()
            if ttype in ("DATA", "FLAG", "FREQ") and (not tform.startswith(("P", "Q"))):
                eh[f"TDIM{i}"] = (f"({int(nchan_hint)})", "Vector length per row")


    if warn_dropped_vector_columns and dropped_cols:
        warnings.warn(
            "Dropped columns with unknown type (all values null/NaN): " + ", ".join(sorted(set(dropped_cols))),
            RuntimeWarning,
        )

    return hdu, int(nchan_hint)
