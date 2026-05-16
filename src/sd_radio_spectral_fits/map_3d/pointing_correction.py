"""Runtime pointing-correction helpers for moving-body maps.

The correction is intentionally lightweight and has no Astropy/SPICE imports.
It is applied in the Moon/Sun sky-offset tangent plane before optional Moon
surface intersection, so callers may pass it only when needed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd


_DX_KEYS = ("dx_arcsec", "dx", "x_arcsec", "dra_arcsec", "daz_arcsec")
_DY_KEYS = ("dy_arcsec", "dy", "y_arcsec", "ddec_arcsec", "del_arcsec")


def _is_scalar_number(value) -> bool:
    if isinstance(value, np.generic):
        value = value.item()
    return isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)


def _is_pair(value) -> bool:
    if isinstance(value, (str, bytes, bytearray, Mapping)):
        return False
    if not isinstance(value, Sequence):
        return False
    if len(value) != 2:
        return False
    return _is_scalar_number(value[0]) and _is_scalar_number(value[1])


def _mapping_get_any(mapping: Mapping, keys: tuple[str, ...], *, default=None):
    for key in keys:
        for cand in (key, key.upper(), key.lower()):
            if cand in mapping:
                return mapping[cand]
    return default


def _coerce_pair(value, *, default=(0.0, 0.0)) -> tuple[float, float]:
    """Coerce a correction object to ``(dx_arcsec, dy_arcsec)``."""
    if value is None:
        return (float(default[0]), float(default[1]))
    if isinstance(value, Mapping):
        dx = _mapping_get_any(value, _DX_KEYS, default=default[0])
        dy = _mapping_get_any(value, _DY_KEYS, default=default[1])
        return (float(dx), float(dy))
    if _is_pair(value):
        return (float(value[0]), float(value[1]))
    raise ValueError(
        "pointing_correction must be None, a (dx_arcsec, dy_arcsec) pair, "
        "a dict with dx_arcsec/dy_arcsec, a list of such entries, or a dict keyed by scantable index."
    )


def select_pointing_correction_for_index(pointing_correction, index: int | None):
    """Select one correction entry for a scantable index.

    Accepted forms:

    - ``None``
    - ``(dx_arcsec, dy_arcsec)``
    - ``{"dx_arcsec": dx, "dy_arcsec": dy}``
    - ``[(dx0, dy0), (dx1, dy1), ...]``
    - ``{0: {"dx_arcsec": ...}, 1: (...), "default": (...)}``
    """
    if pointing_correction is None or index is None:
        return pointing_correction
    if isinstance(pointing_correction, Mapping):
        if any(k in pointing_correction for k in _DX_KEYS + _DY_KEYS):
            return pointing_correction
        for key in (index, str(index), f"scantable_{index}", f"table_{index}"):
            if key in pointing_correction:
                return pointing_correction[key]
        for key in ("default", "DEFAULT", "*"):
            if key in pointing_correction:
                return pointing_correction[key]
        return None
    if isinstance(pointing_correction, Sequence) and not isinstance(pointing_correction, (str, bytes, bytearray)):
        if _is_pair(pointing_correction):
            return pointing_correction
        if 0 <= int(index) < len(pointing_correction):
            return pointing_correction[int(index)]
        return None
    return pointing_correction




def select_pointing_component_for_index(value, index: int | None, *, default: float = 0.0):
    """Select one scalar dx/dy component for a scantable index.

    This helper is used when a public wrapper has filtered or re-ordered a
    list of scantables before passing them to a lower-level routine that will
    enumerate the effective list from zero.  Scalars are returned unchanged;
    sequences and mappings are resolved against the original input index.
    """
    if value is None:
        return None
    if _is_scalar_number(value):
        return float(value)
    if index is None:
        return _component_for_index(value, None, default=default)
    return _component_for_index(value, int(index), default=default)


def _component_for_index(value, index: int | None, *, default: float = 0.0) -> float:
    if value is None:
        return float(default)
    if _is_scalar_number(value):
        return float(value)
    if index is None:
        # A non-scalar component without an index is ambiguous.
        raise ValueError(
            "pointing_dx_arcsec/pointing_dy_arcsec sequence or mapping requires a scantable index "
            "or an OTF_TABLE_INDEX column."
        )
    if isinstance(value, Mapping):
        for key in (index, str(index), f"scantable_{index}", f"table_{index}"):
            if key in value:
                return float(value[key])
        for key in ("default", "DEFAULT", "*"):
            if key in value:
                return float(value[key])
        return float(default)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if 0 <= int(index) < len(value):
            return float(value[int(index)])
        return float(default)
    return float(value)


def _row_indices_from_table(table: pd.DataFrame, *, explicit_index: int | None = None) -> np.ndarray:
    n = int(len(table))
    if explicit_index is not None:
        return np.full(n, int(explicit_index), dtype=np.int64)
    for col in ("OTF_TABLE_INDEX", "TABLE_INDEX", "SCANTABLE_INDEX", "INPUT_INDEX"):
        if col in table.columns:
            arr = pd.to_numeric(table[col], errors="coerce").to_numpy(dtype=float)
            out = np.full(n, -1, dtype=np.int64)
            finite = np.isfinite(arr)
            out[finite] = arr[finite].astype(np.int64)
            return out
    return np.full(n, -1, dtype=np.int64)


def resolve_pointing_offsets_arcsec(
    table: pd.DataFrame,
    *,
    pointing_correction=None,
    pointing_dx_arcsec=None,
    pointing_dy_arcsec=None,
    table_index: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return per-row pointing correction arrays in arcsec.

    Sign convention:

    ``corrected sky offset = original sky offset + (dx_arcsec, dy_arcsec)``.

    Thus, if a Moon map made in ``moon_sky_offset`` shows the fitted lunar
    center at ``(+12, -8) arcsec``, the correction that recenters it is
    ``dx_arcsec=-12`` and ``dy_arcsec=+8``.
    """
    n = int(len(table))
    if n == 0:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)

    row_idx = _row_indices_from_table(table, explicit_index=table_index)
    dx = np.zeros(n, dtype=float)
    dy = np.zeros(n, dtype=float)

    if pointing_dx_arcsec is not None or pointing_dy_arcsec is not None:
        for idx in np.unique(row_idx):
            mask = row_idx == idx
            idx_use = None if int(idx) < 0 else int(idx)
            dx[mask] = _component_for_index(pointing_dx_arcsec, idx_use, default=0.0)
            dy[mask] = _component_for_index(pointing_dy_arcsec, idx_use, default=0.0)
    elif pointing_correction is not None:
        for idx in np.unique(row_idx):
            mask = row_idx == idx
            idx_use = None if int(idx) < 0 else int(idx)
            corr = select_pointing_correction_for_index(pointing_correction, idx_use)
            dx_val, dy_val = _coerce_pair(corr, default=(0.0, 0.0))
            dx[mask] = dx_val
            dy[mask] = dy_val

    if not (np.all(np.isfinite(dx)) and np.all(np.isfinite(dy))):
        raise ValueError("pointing correction contains NaN or non-finite values.")
    return dx, dy


def has_nonzero_pointing_correction(
    table: pd.DataFrame,
    *,
    pointing_correction=None,
    pointing_dx_arcsec=None,
    pointing_dy_arcsec=None,
    table_index: int | None = None,
) -> bool:
    dx, dy = resolve_pointing_offsets_arcsec(
        table,
        pointing_correction=pointing_correction,
        pointing_dx_arcsec=pointing_dx_arcsec,
        pointing_dy_arcsec=pointing_dy_arcsec,
        table_index=table_index,
    )
    return bool(np.any(dx != 0.0) or np.any(dy != 0.0))
