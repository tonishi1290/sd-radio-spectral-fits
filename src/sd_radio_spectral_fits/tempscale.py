# src/sd_radio_spectral_fits/tempscale.py
from __future__ import annotations

import datetime
import warnings
from typing import Any, Iterable, Optional, Tuple, Union, List, TYPE_CHECKING, Sequence
from collections.abc import Mapping

import numpy as np
import pandas as pd

# [MODIFIED] Avoid circular import by using TYPE_CHECKING
if TYPE_CHECKING:
    from .fitsio import Scantable

# -----------------------------------------------------------------------------
# Temperature scale utilities
# -----------------------------------------------------------------------------
# Policy:
# - Internal spectral arrays are treated as Ta*.
# - Tr* is derived as Ta* / BEAMEFF (main-beam efficiency).
# - The on-disk scale is declared by TEMPSCAL ("TA*" or "TR*").
# - Never switch scales silently: if conversion is performed, record it in HISTORY/log.
# -----------------------------------------------------------------------------


def normalize_tempscal(value: Any, *, default: str = "TA*") -> str:
    """Normalize various inputs into canonical 'TA*' or 'TR*'."""
    if value is None:
        return str(default).upper()
    s = str(value).strip().upper()
    if not s:
        return str(default).upper()
    # Common variants
    if s in ("TA*", "TA", "TASTAR", "TASTAR*", "T_A*", "T_A", "TA-STAR", "TASTAR-STAR"):
        return "TA*"
    if s in ("TR*", "TR", "TMB", "T_MB", "TMAINBEAM", "MAINBEAM", "TMB*"):
        return "TR*"
    # Unknown -> default (safer than propagating nonsense)
    return str(default).upper()


def ensure_tempscal_column(df: pd.DataFrame, *, default: str = "TA*") -> pd.DataFrame:
    """Ensure df has a TEMPSCAL column (string)."""
    if df is None:
        return df
    if "TEMPSCAL" not in df.columns:
        df = df.copy()
        df["TEMPSCAL"] = normalize_tempscal(default, default=default)
        return df
    # normalize in-place on a copy (avoid modifying caller unexpectedly)
    out = df.copy()
    try:
        out["TEMPSCAL"] = [normalize_tempscal(v, default=default) for v in out["TEMPSCAL"].to_numpy()]
    except Exception:
        out["TEMPSCAL"] = normalize_tempscal(default, default=default)
    return out


def ensure_beameff_column(df: pd.DataFrame, *, default: float = np.nan) -> pd.DataFrame:
    """Ensure df has a BEAMEFF column (float)."""
    if df is None:
        return df
    if "BEAMEFF" not in df.columns:
        df = df.copy()
        df["BEAMEFF"] = float(default)
        return df
    out = df.copy()
    try:
        out["BEAMEFF"] = pd.to_numeric(out["BEAMEFF"], errors="coerce").astype(float)
    except Exception:
        out["BEAMEFF"] = float(default)
    return out


def beameff_array(
    df: pd.DataFrame,
    meta: dict | None,
    n_rows: int,
    *,
    default: float = np.nan,
) -> np.ndarray:
    """Return BEAMEFF per row (length n_rows)."""
    if df is not None and "BEAMEFF" in df.columns:
        try:
            arr = pd.to_numeric(df["BEAMEFF"], errors="coerce").to_numpy(dtype=float)
            if arr.size >= n_rows:
                return arr[:n_rows]
            if arr.size > 0:
                out = np.full(int(n_rows), float(default), dtype=float)
                out[:arr.size] = arr
                return out
        except Exception:
            pass

    # fallback: global keyword
    if meta is not None:
        for k in ("BEAMEFF",):
            if k in meta and meta[k] not in (None, ""):
                try:
                    v = float(meta[k])
                    return np.full(int(n_rows), v, dtype=float)
                except Exception:
                    pass

    return np.full(int(n_rows), float(default), dtype=float)


def tempscal_array(
    df: pd.DataFrame,
    meta: dict | None,
    n_rows: int,
    *,
    default: str = "TA*",
) -> np.ndarray:
    """Return TEMPSCAL per row (length n_rows) as normalized strings."""
    if df is not None and "TEMPSCAL" in df.columns:
        try:
            col = df["TEMPSCAL"].to_numpy()
            arr = np.array([normalize_tempscal(v, default=default) for v in col], dtype="U4")
            if arr.size >= n_rows:
                return arr[:n_rows]
            out = np.full(int(n_rows), normalize_tempscal(default, default=default), dtype="U4")
            out[:arr.size] = arr
            return out
        except Exception:
            pass

    # fallback: header keyword
    if meta is not None:
        for k in ("TEMPSCAL", "TEMP_SCAL", "TSCALE"):
            if k in meta and meta[k] not in (None, ""):
                return np.full(int(n_rows), normalize_tempscal(meta[k], default=default), dtype="U4")
    return np.full(int(n_rows), normalize_tempscal(default, default=default), dtype="U4")


def require_beameff(beameff: np.ndarray, *, on_fail: str = "error") -> None:
    """
    Validate BEAMEFF (main-beam efficiency).

    - Must be finite and > 0 for Tr* conversions.
    - Internal spectral arrays are treated as Ta*; Tr* is derived as Ta* / BEAMEFF.
    """
    b = np.asarray(beameff, dtype=float)
    bad = (~np.isfinite(b)) | (b <= 0.0)
    if np.any(bad):
        msg = "BEAMEFF is required (finite and > 0) for TR* conversions, but missing/invalid values were found."
        if on_fail == "warn":
            warnings.warn(msg, RuntimeWarning)
            return
        raise ValueError(msg)


def is_beameff_mixed(beameff: np.ndarray, *, tol: float = 1e-6) -> bool:
    """
    Check whether BEAMEFF is mixed beyond tolerance.

    `tol` is a *relative* tolerance:
        mixed if (max - min) > tol * median(|beameff|).

    Any non-positive or non-finite BEAMEFF is treated as mixed (because TR* conversion would be unsafe).
    """
    b = np.asarray(beameff, dtype=float).ravel()
    b = b[np.isfinite(b)]
    if b.size <= 1:
        return False
    bmin = float(np.min(b))
    bmax = float(np.max(b))
    if bmin <= 0.0:
        return True
    bref = float(np.median(b))
    scale = max(abs(bref), 1e-12)
    return (bmax - bmin) > (float(tol) * scale)


def representative_beameff(beameff: np.ndarray) -> float:
    """Return a representative BEAMEFF (median of finite positive values)."""
    b = np.asarray(beameff, dtype=float).ravel()
    b = b[np.isfinite(b) & (b > 0.0)]
    if b.size == 0:
        return float("nan")
    return float(np.median(b))


def ta_to_tr(ta: np.ndarray, beameff: np.ndarray) -> np.ndarray:
    """Convert Ta* -> Tr* (vectorized): Tr* = Ta* / BEAMEFF."""
    require_beameff(beameff, on_fail="error")
    ta = np.asarray(ta, dtype=float)
    b = np.asarray(beameff, dtype=float)

    if ta.ndim == 2 and b.ndim == 1 and ta.shape[0] == b.shape[0]:
        out = np.empty_like(ta, dtype=float)
        np.divide(ta, b[:, None], out=out, where=(b[:, None] > 0.0))
        return out.astype(np.float32, copy=False)

    out = ta / b
    return np.asarray(out, dtype=np.float32)


def tr_to_ta(tr: np.ndarray, beameff: np.ndarray) -> np.ndarray:
    """Convert Tr* -> Ta* (vectorized): Ta* = Tr* * BEAMEFF."""
    require_beameff(beameff, on_fail="error")
    tr = np.asarray(tr, dtype=float)
    b = np.asarray(beameff, dtype=float)

    if tr.ndim == 2 and b.ndim == 1 and tr.shape[0] == b.shape[0]:
        out = (tr * b[:, None]).astype(np.float32, copy=False)
        return out

    out = tr * b
    return np.asarray(out, dtype=np.float32)


def convert_rowwise_vla(
    spectra: List[np.ndarray],
    beameff: np.ndarray,
    *,
    direction: str,
) -> List[np.ndarray]:
    """Row-wise Ta*<->Tr* conversion for ragged (variable-length) spectra.

    Parameters
    ----------
    spectra
        List of 1D numpy arrays (one per row). Each element may have different length.
    beameff
        Per-row BEAMEFF array. If length==1, it is broadcast to all rows.
    direction
        "ta_to_tr" or "tr_to_ta".

    Returns
    -------
    list of np.ndarray
        Converted spectra as float32 arrays, preserving original per-row shapes.
    """
    if not isinstance(spectra, list):
        spectra = list(spectra)

    b = np.asarray(beameff, dtype=float).ravel()
    if b.size == 1 and len(spectra) > 1:
        b = np.full(len(spectra), float(b[0]), dtype=float)
    if len(spectra) != int(b.size):
        raise ValueError(f"convert_rowwise_vla: len(spectra)={len(spectra)} != len(BEAMEFF)={int(b.size)}")

    require_beameff(b, on_fail="error")

    d = str(direction).strip().lower()
    out: List[np.ndarray] = []
    if d in ("ta_to_tr", "ta2tr", "ta->tr"):
        for spec, bi in zip(spectra, b):
            arr = np.asarray(spec, dtype=float)
            out.append((arr / float(bi)).astype(np.float32, copy=False))
        return out
    if d in ("tr_to_ta", "tr2ta", "tr->ta"):
        for spec, bi in zip(spectra, b):
            arr = np.asarray(spec, dtype=float)
            out.append((arr * float(bi)).astype(np.float32, copy=False))
        return out

    raise ValueError("convert_rowwise_vla: direction must be 'ta_to_tr' or 'tr_to_ta'")



def append_scale_history(history: dict | None, entry: dict) -> dict:
    """Append a scale-related history entry (non-breaking)."""
    if history is None or not isinstance(history, dict):
        history = {} if history is None else {"_history": str(history)}
    hist = dict(history)
    key = "scale_history"
    if key not in hist or not isinstance(hist.get(key), list):
        hist[key] = []
    entry2 = dict(entry)
    entry2.setdefault("created_at_utc", datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"))
    hist[key].append(entry2)
    return hist



# -----------------------------------------------------------------------------
# Row-wise metadata / intensity scaling helpers
# -----------------------------------------------------------------------------

_SCALE_DERIVED_CLEAR_COLS: Tuple[str, ...] = (
    "BSL_RMS",
    "BSL_COEF",
    "BSL_SCALE",
    "MOMENT0",
    "MOM0",
    "INTEGRATED",
)


def _ensure_intscale_column(df: pd.DataFrame) -> None:
    """Ensure an in-place cumulative intensity scale column exists."""
    if "INTSCALE" not in df.columns:
        df["INTSCALE"] = np.ones(len(df), dtype=np.float64)
        return
    vals = pd.to_numeric(df["INTSCALE"], errors="coerce").to_numpy(dtype=float)
    bad = ~np.isfinite(vals)
    if np.any(bad):
        vals[bad] = 1.0
    df["INTSCALE"] = vals.astype(np.float64, copy=False)


def _validate_positive_factors(values: np.ndarray, *, context: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0:
        raise ValueError(f"{context}: no factor values were provided.")
    bad = (~np.isfinite(arr)) | (arr <= 0.0)
    if np.any(bad):
        examples = arr[bad][:5]
        raise ValueError(
            f"{context}: factors must be finite and > 0. Invalid examples: {examples.tolist()}"
        )
    return arr


def _selected_index_labels(df: pd.DataFrame, idxs: np.ndarray) -> pd.Index:
    return df.index[np.asarray(idxs, dtype=int)]


def _normalize_key_columns(key_columns: Union[str, Sequence[str]]) -> List[str]:
    cols = [str(key_columns)] if isinstance(key_columns, str) else [str(c) for c in key_columns]
    if len(cols) == 0:
        raise ValueError("key_columns must not be empty.")
    return cols


def _resolve_row_factors(
    df: pd.DataFrame,
    idxs: np.ndarray,
    scale: Union[float, Sequence[float], np.ndarray, Mapping[Any, float]],
    *,
    key_columns: Union[str, Sequence[str]] = ("FDNUM", "IFNUM", "PLNUM"),
    strict: bool = False,
    context: str,
) -> Tuple[np.ndarray, str, dict]:
    idxs = np.asarray(idxs, dtype=int)
    if idxs.size == 0:
        return np.array([], dtype=float), "empty", {}

    if isinstance(scale, Mapping):
        cols = _normalize_key_columns(key_columns)
        missing_cols = [c for c in cols if c not in df.columns]
        if missing_cols:
            raise KeyError(f"{context}: key columns not found: {missing_cols}")

        normalized_map: dict[Any, float] = {}
        for raw_key, raw_fac in scale.items():
            fac = float(_validate_positive_factors(np.array([raw_fac], dtype=float), context=f"{context}(mapping)")[0])
            if len(cols) == 1:
                if isinstance(raw_key, tuple):
                    if len(raw_key) != 1:
                        raise ValueError(
                            f"{context}: key {raw_key!r} does not match single key column {cols[0]!r}."
                        )
                    key = raw_key[0]
                else:
                    key = raw_key
            else:
                if not isinstance(raw_key, tuple) or len(raw_key) != len(cols):
                    raise ValueError(
                        f"{context}: key {raw_key!r} must be a tuple of length {len(cols)} for key_columns={cols}."
                    )
                key = tuple(raw_key)
            normalized_map[key] = fac

        subset = df.iloc[idxs]
        used_keys: set[Any] = set()

        if len(cols) == 1:
            key_series = subset[cols[0]]
            factors_series = key_series.map(normalized_map)
            factors = pd.to_numeric(factors_series, errors="coerce").to_numpy(dtype=float)
            finite = np.isfinite(factors)
            if finite.any():
                used_keys = set(pd.unique(key_series[finite]))
            if not finite.all():
                miss_pos = np.flatnonzero(~finite)
                missing_rows = [(int(idxs[m]), key_series.iloc[m]) for m in miss_pos[:5]]
                raise ValueError(
                    f"{context}: {int((~finite).sum())} selected rows have no matching factor. Examples: {missing_rows}"
                )
        else:
            map_rows = [tuple(k) + (v,) for k, v in normalized_map.items()]
            map_df = pd.DataFrame(map_rows, columns=cols + ["__factor__"])
            merged = subset[cols].reset_index(drop=True).merge(map_df, on=cols, how="left", sort=False)
            factors = pd.to_numeric(merged["__factor__"], errors="coerce").to_numpy(dtype=float)
            finite = np.isfinite(factors)
            if finite.any():
                used_keys = set(merged.loc[finite, cols].itertuples(index=False, name=None))
            if not finite.all():
                miss_pos = np.flatnonzero(~finite)
                missing_rows = []
                for m in miss_pos[:5]:
                    key = tuple(subset.iloc[int(m)][c] for c in cols)
                    missing_rows.append((int(idxs[m]), key))
                raise ValueError(
                    f"{context}: {int((~finite).sum())} selected rows have no matching factor. Examples: {missing_rows}"
                )

        used_key_counts: dict[Any, int] = {}
        if len(cols) == 1:
            if finite.any():
                used_vals, counts = np.unique(key_series[finite].to_numpy(), return_counts=True)
                used_key_counts = {used_vals[i]: int(counts[i]) for i in range(len(used_vals))}
        else:
            if finite.any():
                finite_keys = list(merged.loc[finite, cols].itertuples(index=False, name=None))
                for key in finite_keys:
                    used_key_counts[key] = used_key_counts.get(key, 0) + 1

        unused_keys = [k for k in normalized_map.keys() if k not in used_keys]
        if unused_keys:
            msg = f"{context}: {len(unused_keys)} mapping keys were unused. Examples: {unused_keys[:5]}"
            if strict:
                raise ValueError(msg)
            warnings.warn(msg)

        applied_summary = []
        for key in normalized_map.keys():
            n_rows = int(used_key_counts.get(key, 0))
            if n_rows <= 0:
                continue
            applied_summary.append({
                "key": key,
                "factor": float(normalized_map[key]),
                "n_rows": n_rows,
            })

        return factors, "mapping", {
            "key_columns": cols,
            "n_map": len(normalized_map),
            "n_unused": len(unused_keys),
            "applied_summary": applied_summary,
        }

    arr = np.asarray(scale, dtype=float)
    if arr.ndim == 0:
        fac = float(_validate_positive_factors(np.array([arr.item()], dtype=float), context=f"{context}(scalar)")[0])
        return np.full(idxs.size, fac, dtype=float), "scalar", {"factor": fac}

    vals = _validate_positive_factors(arr, context=f"{context}(sequence)")
    if vals.size != idxs.size:
        raise ValueError(
            f"{context}: factor length ({int(vals.size)}) must equal number of selected rows ({int(idxs.size)})."
        )
    return vals.astype(float, copy=False), "sequence", {"n_factors": int(vals.size)}


def _scale_row_array_inplace_or_replace(arr: np.ndarray, factor: float) -> np.ndarray:
    a = np.asarray(arr)
    if np.issubdtype(a.dtype, np.floating):
        try:
            a *= np.asarray(factor, dtype=a.dtype)
            return a
        except Exception:
            pass
    out = np.asarray(a, dtype=np.float32) * np.float32(factor)
    return np.asarray(out, dtype=np.float32)


def _scale_scantable_data_inplace(sc: "Scantable", idxs: np.ndarray, row_factors: np.ndarray) -> None:
    idxs = np.asarray(idxs, dtype=int)
    fac = np.asarray(row_factors, dtype=float).ravel()
    if fac.size != idxs.size:
        raise ValueError("_scale_scantable_data_inplace: len(row_factors) must equal len(idxs).")

    data = sc.data
    scalar_like = bool(fac.size > 0 and np.allclose(fac, fac[0], rtol=0.0, atol=0.0))
    if isinstance(data, np.ndarray):
        if data.ndim != 2:
            raise ValueError(f"_scale_scantable_data_inplace: expected 2D ndarray, got shape={data.shape}")
        if not np.issubdtype(data.dtype, np.floating):
            data = np.asarray(data, dtype=np.float32)
            sc.data = data
        if idxs.size == data.shape[0]:
            if scalar_like:
                data *= np.asarray(float(fac[0]), dtype=data.dtype)
            else:
                data *= fac.astype(data.dtype, copy=False)[:, None]
            return
        if scalar_like:
            data[idxs, :] *= np.asarray(float(fac[0]), dtype=data.dtype)
        else:
            data[idxs, :] = data[idxs, :] * fac.astype(data.dtype, copy=False)[:, None]
        return

    if not isinstance(data, list):
        raise TypeError(f"_scale_scantable_data_inplace: unsupported data type {type(data)}")

    if scalar_like:
        f0 = float(fac[0])
        for i in idxs:
            data[int(i)] = _scale_row_array_inplace_or_replace(data[int(i)], f0)
        return

    for i, f in zip(idxs, fac):
        data[int(i)] = _scale_row_array_inplace_or_replace(data[int(i)], float(f))


def _invalidate_scale_dependent_columns(df: pd.DataFrame, idxs: np.ndarray) -> List[str]:
    cleared: List[str] = []
    labels = _selected_index_labels(df, idxs)
    for col in _SCALE_DERIVED_CLEAR_COLS:
        if col not in df.columns:
            continue
        series = df[col]
        if getattr(series.dtype, "kind", "") in ("f", "i", "u"):
            df.loc[labels, col] = np.nan
        else:
            df.loc[labels, col] = None
        cleared.append(col)
    return cleared


def _apply_scale_common(
    sc: "Scantable",
    row_factors: np.ndarray,
    idxs: np.ndarray,
    *,
    action: str,
    factor_mode: str,
    history_extra: Optional[dict] = None,
    invalidate_derived: bool = True,
    verbose: bool = True,
) -> None:
    df = sc.table
    idxs = np.asarray(idxs, dtype=int)
    row_factors = np.asarray(row_factors, dtype=float).ravel()
    if idxs.size == 0:
        if verbose:
            print("No rows selected.")
        return

    _scale_scantable_data_inplace(sc, idxs, row_factors)
    _ensure_intscale_column(df)
    # pandas may return a read-only NumPy view here; take an explicit writable copy.
    current = pd.to_numeric(df["INTSCALE"], errors="coerce").to_numpy(dtype=float, copy=True)
    bad = ~np.isfinite(current)
    if np.any(bad):
        current[bad] = 1.0
    current[idxs] *= row_factors
    df["INTSCALE"] = current.astype(np.float64, copy=False)

    cleared = _invalidate_scale_dependent_columns(df, idxs) if invalidate_derived else []

    entry = {
        "stage": "tempscale",
        "action": action,
        "factor_mode": factor_mode,
        "rows_selected": int(idxs.size),
        "invalidate_derived": bool(invalidate_derived),
        "cleared_columns": cleared,
        "tempscal_preserved": True,
        "record_column": "INTSCALE",
    }
    if row_factors.size > 0:
        entry.update({
            "factor_min": float(np.min(row_factors)),
            "factor_max": float(np.max(row_factors)),
        })
    if history_extra:
        entry.update(history_extra)
    sc.history = append_scale_history(sc.history, entry)

    if verbose:
        current_scale = df["TEMPSCAL"].iloc[int(idxs[0])] if "TEMPSCAL" in df.columns else "Unknown"
        print(f"Applied {action} to {len(idxs)} rows (mode={factor_mode}).")
        print(f"Note: Data remain labeled as '{current_scale}'. Cumulative factor is recorded in 'INTSCALE'.")
        if factor_mode == "mapping" and history_extra:
            key_columns = history_extra.get("key_columns")
            applied_summary = history_extra.get("applied_summary") or []
            if key_columns and applied_summary:
                key_label = ", ".join(str(c) for c in key_columns)
                print(f"Mapping summary by ({key_label}):")
                for item in applied_summary:
                    key = item.get("key")
                    factor = float(item.get("factor"))
                    n_rows = int(item.get("n_rows", 0))
                    print(f"  key={key!r} -> factor={factor:.15g}, rows={n_rows}")
        if cleared:
            print(f"Invalidated scale-dependent columns: {', '.join(cleared)}")


def set_beameff(
    sc: "Scantable",
    efficiency: Union[float, Sequence[float], np.ndarray, Mapping[Any, float]],
    rows: Union[str, List[int], int, slice, None] = None,
    *,
    key_columns: Union[str, Sequence[str]] = ("FDNUM", "IFNUM", "PLNUM"),
    strict: bool = False,
    verbose: bool = True,
) -> None:
    """
    Set per-row BEAMEFF values in-place.

    This function only updates the ``BEAMEFF`` column. It does **not** modify
    the spectral data array and does **not** change ``TEMPSCAL``.

    Accepted efficiency forms
    -------------------------
    1. scalar
       Apply one value to all selected rows.
    2. 1-D sequence / ndarray
       Length must exactly match the number of selected rows.
    3. mapping
       Apply values by key columns such as ``(FDNUM, IFNUM, PLNUM)``.
    """
    from .scantable_utils import _parse_row_selector

    df = sc.table
    idxs = np.asarray(_parse_row_selector(rows, len(df)), dtype=int)
    if idxs.size == 0:
        if verbose:
            print("No rows selected.")
        return

    if "BEAMEFF" not in df.columns:
        df["BEAMEFF"] = np.nan

    values, mode, info = _resolve_row_factors(
        df,
        idxs,
        efficiency,
        key_columns=key_columns,
        strict=strict,
        context="set_beameff",
    )
    if np.any(values > 1.0):
        bad = values[values > 1.0][:5]
        raise ValueError(
            "set_beameff: BEAMEFF values must satisfy 0 < eta <= 1. "
            f"Invalid examples: {bad.tolist()}"
        )
    df.loc[_selected_index_labels(df, idxs), "BEAMEFF"] = values.astype(float, copy=False)

    if verbose:
        current_scale = df["TEMPSCAL"].iloc[int(idxs[0])] if "TEMPSCAL" in df.columns else "Unknown"
        print(f"Updated BEAMEFF for {len(idxs)} rows (mode={mode}).")
        print(f"Note: Data remain in '{current_scale}'. Use Viewer / coadd / write_sdfits to interpret TR*.")

    hist = {"stage": "tempscale", "action": "set_beameff", "mode": mode, "rows_selected": int(idxs.size)}
    hist.update(info)
    sc.history = append_scale_history(sc.history, hist)


def apply_relative_scale(
    sc: "Scantable",
    scale: Union[float, Sequence[float], np.ndarray, Mapping[Any, float]],
    rows: Union[str, List[int], int, slice, None] = None,
    *,
    key_columns: Union[str, Sequence[str]] = ("FDNUM", "IFNUM", "PLNUM"),
    strict: bool = False,
    invalidate_derived: bool = True,
    verbose: bool = True,
) -> None:
    """
    Apply row-wise relative intensity scale factors in-place.

    This function multiplies the spectral data themselves and preserves ``TEMPSCAL``.
    It is intended for beam-to-beam relative scaling, typically immediately after
    calibration and before baseline / RMS / moment / coadd products are computed.

    The cumulative applied factor is recorded in ``INTSCALE``.
    By default, known scale-dependent derived columns are invalidated.
    """
    from .scantable_utils import _parse_row_selector

    df = sc.table
    idxs = np.asarray(_parse_row_selector(rows, len(df)), dtype=int)
    if idxs.size == 0:
        if verbose:
            print("No rows selected.")
        return

    factors, mode, info = _resolve_row_factors(
        df,
        idxs,
        scale,
        key_columns=key_columns,
        strict=strict,
        context="apply_relative_scale",
    )
    _apply_scale_common(
        sc,
        factors,
        idxs,
        action="apply_relative_scale",
        factor_mode=mode,
        history_extra=info,
        invalidate_derived=invalidate_derived,
        verbose=verbose,
    )


def apply_global_scale(
    sc: "Scantable",
    factor: float,
    rows: Union[str, List[int], int, slice, None] = None,
    *,
    invalidate_derived: bool = True,
    verbose: bool = True,
) -> None:
    """
    Apply one multiplicative global intensity factor in-place.

    Typical use: after beam-to-beam relative alignment has already been done,
    multiply all selected rows by one final absolute scale factor.
    ``TEMPSCAL`` is preserved and the cumulative factor is recorded in ``INTSCALE``.
    """
    from .scantable_utils import _parse_row_selector

    df = sc.table
    idxs = np.asarray(_parse_row_selector(rows, len(df)), dtype=int)
    if idxs.size == 0:
        if verbose:
            print("No rows selected.")
        return

    factors, mode, info = _resolve_row_factors(
        df,
        idxs,
        float(factor),
        key_columns=("FDNUM", "IFNUM", "PLNUM"),
        strict=False,
        context="apply_global_scale",
    )
    _apply_scale_common(
        sc,
        factors,
        idxs,
        action="apply_global_scale",
        factor_mode=mode,
        history_extra=info,
        invalidate_derived=invalidate_derived,
        verbose=verbose,
    )
