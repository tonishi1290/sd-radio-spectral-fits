from __future__ import annotations

import builtins
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import warnings


QC_MAD_SCALE = 1.4826


def ensure_flagrow_column(table: pd.DataFrame) -> pd.DataFrame:
    out = table.copy()
    if "FLAGROW" not in out.columns:
        out["FLAGROW"] = np.zeros(len(out), dtype=np.int32)
        return out
    vals = pd.to_numeric(out["FLAGROW"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    out["FLAGROW"] = vals.astype(np.int32)
    return out


def flagrow_mask(table: pd.DataFrame) -> np.ndarray:
    if "FLAGROW" not in table.columns:
        return np.zeros(len(table), dtype=bool)
    vals = pd.to_numeric(table["FLAGROW"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    return vals != 0


def apply_flagrow(table: pd.DataFrame, bad_mask: np.ndarray, *, code: int = 1) -> pd.DataFrame:
    out = ensure_flagrow_column(table)
    bad = np.asarray(bad_mask, dtype=bool)
    if bad.size != len(out):
        raise ValueError("bad_mask length does not match table length")
    vals = pd.to_numeric(out["FLAGROW"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    vals[bad] = builtins.max(int(code), 1)
    out["FLAGROW"] = vals.astype(np.int32)
    return out


def robust_mad_scale(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan
    med = float(np.median(arr))
    mad = float(np.median(np.abs(arr - med)))
    return QC_MAD_SCALE * mad


def robust_center(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def affine_shape_score(row: np.ndarray, template: np.ndarray) -> float:
    y = np.asarray(row, dtype=float)
    t = np.asarray(template, dtype=float)
    mask = np.isfinite(y) & np.isfinite(t)
    if np.count_nonzero(mask) < 3:
        return np.nan
    A = np.column_stack([np.ones(np.count_nonzero(mask), dtype=float), t[mask]])
    coef, *_ = np.linalg.lstsq(A, y[mask], rcond=None)
    resid = y[mask] - (A @ coef)
    return robust_mad_scale(resid)


def grouped_affine_shape_qc(
    table: pd.DataFrame,
    data: np.ndarray,
    *,
    qc_sigma: float | None,
    candidate_mask: np.ndarray | None = None,
    group_cols: Sequence[str] | None = None,
    min_group_size: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    n = len(table)
    scores = np.full(n, np.nan, dtype=float)
    bad = np.zeros(n, dtype=bool)
    if qc_sigma is None or n == 0:
        return bad, scores
    if data.shape[0] != n:
        raise ValueError("data row count does not match table length")

    cand = np.ones(n, dtype=bool) if candidate_mask is None else np.asarray(candidate_mask, dtype=bool)
    if cand.size != n:
        raise ValueError("candidate_mask length does not match table length")

    if group_cols:
        work = table.loc[:, list(group_cols)].copy()
        for col in work.columns:
            if pd.api.types.is_numeric_dtype(work[col]):
                work[col] = pd.to_numeric(work[col], errors="coerce").fillna(-999999)
            else:
                work[col] = work[col].astype(str).fillna("")
        grouped = work.groupby(list(work.columns), sort=False, dropna=False).indices
        group_lists = [np.asarray(idxs, dtype=int) for idxs in grouped.values()]
    else:
        group_lists = [np.arange(n, dtype=int)]

    for idxs in group_lists:
        idxs = idxs[cand[idxs]]
        if idxs.size < builtins.max(2, int(min_group_size)):
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            template = np.nanmedian(np.asarray(data[idxs], dtype=float), axis=0)
        local_scores = np.full(idxs.size, np.nan, dtype=float)
        for j, ridx in enumerate(idxs):
            local_scores[j] = affine_shape_score(data[ridx], template)
        scores[idxs] = local_scores
        finite = np.isfinite(local_scores)
        if np.count_nonzero(finite) < builtins.max(2, int(min_group_size)):
            continue
        med = float(np.median(local_scores[finite]))
        scl = robust_mad_scale(local_scores[finite])
        if not np.isfinite(scl) or scl <= 0:
            continue
        bad[idxs] = np.isfinite(local_scores) & (local_scores > med + float(qc_sigma) * scl)
    return bad, scores


def diffmad_score(row: np.ndarray, window_mask: np.ndarray | None = None) -> float:
    y = np.asarray(row, dtype=float)
    if window_mask is None:
        valid = np.isfinite(y)
    else:
        valid = np.isfinite(y) & np.asarray(window_mask, dtype=bool)
    idx = np.flatnonzero(valid)
    if idx.size < 4:
        return np.nan
    dy = np.diff(y[idx])
    consec = np.diff(idx) == 1
    if dy.size == 0 or not np.any(consec):
        return np.nan
    return robust_mad_scale(dy[consec])


def abspeak_score(row: np.ndarray, window_mask: np.ndarray | None = None) -> float:
    y = np.asarray(row, dtype=float)
    if window_mask is None:
        vals = y[np.isfinite(y)]
    else:
        valid = np.isfinite(y) & np.asarray(window_mask, dtype=bool)
        vals = y[valid]
    if vals.size == 0:
        return np.nan
    med = float(np.median(vals))
    return float(np.max(np.abs(vals - med)))


def grouped_diffmad_qc(
    table: pd.DataFrame,
    data: np.ndarray,
    *,
    qc_sigma: float | None,
    candidate_mask: np.ndarray | None = None,
    group_cols: Sequence[str] | None = None,
    window_mask: np.ndarray | None = None,
    min_group_size: int = 4,
) -> tuple[np.ndarray, np.ndarray]:
    n = len(table)
    scores = np.full(n, np.nan, dtype=float)
    bad = np.zeros(n, dtype=bool)
    if qc_sigma is None or n == 0:
        return bad, scores
    if data.shape[0] != n:
        raise ValueError("data row count does not match table length")

    cand = np.ones(n, dtype=bool) if candidate_mask is None else np.asarray(candidate_mask, dtype=bool)
    if cand.size != n:
        raise ValueError("candidate_mask length does not match table length")

    if group_cols:
        work = table.loc[:, list(group_cols)].copy()
        for col in work.columns:
            if pd.api.types.is_numeric_dtype(work[col]):
                work[col] = pd.to_numeric(work[col], errors="coerce").fillna(-999999)
            else:
                work[col] = work[col].astype(str).fillna("")
        grouped = work.groupby(list(work.columns), sort=False, dropna=False).indices
        group_lists = [np.asarray(idxs, dtype=int) for idxs in grouped.values()]
    else:
        group_lists = [np.arange(n, dtype=int)]

    for idxs in group_lists:
        idxs = idxs[cand[idxs]]
        if idxs.size < builtins.max(2, int(min_group_size)):
            continue
        local_scores = np.full(idxs.size, np.nan, dtype=float)
        for j, ridx in enumerate(idxs):
            local_scores[j] = diffmad_score(data[ridx], window_mask=window_mask)
        scores[idxs] = local_scores
        finite = np.isfinite(local_scores)
        if np.count_nonzero(finite) < builtins.max(2, int(min_group_size)):
            continue
        med = float(np.median(local_scores[finite]))
        scl = robust_mad_scale(local_scores[finite])
        if not np.isfinite(scl) or scl <= 0:
            continue
        bad[idxs] = np.isfinite(local_scores) & (local_scores > med + float(qc_sigma) * scl)
    return bad, scores


def choose_existing_group_cols(table: pd.DataFrame, preferred: Sequence[str]) -> list[str]:
    return [c for c in preferred if c in table.columns]
