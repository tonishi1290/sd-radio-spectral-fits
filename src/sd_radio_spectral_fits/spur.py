from __future__ import annotations

import builtins
import json
import warnings
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .fitsio import Scantable, stamp_scantable_code_provenance


_VALID_MODES = {
    "none": "none",
    "off": "none",
    "false": "none",
    "local": "local",
    "promoted-all-dumps": "promoted-all-dumps",
    "promoted_all_dumps": "promoted-all-dumps",
    "promoted": "promoted-all-dumps",
}


@dataclass(frozen=True)
class SpurConfig:
    mode: str = "none"
    nsigma: float = 7.0
    max_width: int = 2

    @classmethod
    def normalize(
        cls,
        *,
        mode: str | None = "none",
        nsigma: float = 7.0,
        max_width: int = 2,
    ) -> "SpurConfig":
        mode_key = str(mode or "none").strip().lower()
        if mode_key not in _VALID_MODES:
            raise ValueError(
                f"Invalid spur_mode={mode!r}. Must be one of: none, local, promoted-all-dumps."
            )
        mode_norm = _VALID_MODES[mode_key]

        ns = float(nsigma)
        if not np.isfinite(ns) or ns <= 0.0:
            raise ValueError(f"spur_nsigma must be > 0 and finite; got {nsigma!r}.")

        mw = int(max_width)
        if mw < 1:
            raise ValueError(f"spur_max_width must be >= 1; got {max_width!r}.")

        return cls(mode=mode_norm, nsigma=ns, max_width=mw)


def _group_contiguous(indices: np.ndarray) -> list[tuple[int, int]]:
    if indices.size == 0:
        return []
    idx = np.array(indices, dtype=int, copy=True)
    idx.sort()
    runs: list[tuple[int, int]] = []
    start = int(idx[0])
    prev = int(idx[0])
    for val in idx[1:]:
        v = int(val)
        if v == prev + 1:
            prev = v
            continue
        runs.append((start, prev))
        start = v
        prev = v
    runs.append((start, prev))
    return runs


def _mad_sigma_rows(residual_2d: np.ndarray) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered", category=RuntimeWarning)
        med = np.nanmedian(residual_2d, axis=1)
        mad = np.nanmedian(np.abs(residual_2d - med[:, None]), axis=1)
    sigma = 1.4826 * mad
    bad = (~np.isfinite(sigma)) | (sigma <= 0.0)
    sigma[bad] = np.nan
    return sigma.astype(float, copy=False)


def _compute_local_median_residual_2d(
    data: np.ndarray,
    *,
    half_window: int = 2,
) -> np.ndarray:
    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError(f"Expected 2D array; got shape={data.shape!r}.")

    pad = int(half_window)
    if pad < 0:
        raise ValueError(f"half_window must be >= 0; got {half_window!r}.")

    padded = np.pad(data, ((0, 0), (pad, pad)), mode="constant", constant_values=np.nan)
    windows = np.lib.stride_tricks.sliding_window_view(
        padded,
        window_shape=(2 * pad + 1),
        axis=1,
    )
    finite_counts = np.isfinite(windows).sum(axis=-1)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered", category=RuntimeWarning)
        window_median = np.nanmedian(windows, axis=-1)

    resid = data - window_median
    valid = (finite_counts >= 3) & np.isfinite(data)
    resid[~valid] = np.nan
    return resid


def _bridge_segment(spec: np.ndarray, a: int, b: int) -> tuple[np.ndarray | None, np.ndarray | None]:
    nchan = int(spec.size)
    if a <= 0 or b >= (nchan - 1) or b < a:
        return None, None

    left = float(spec[a - 1])
    right = float(spec[b + 1])
    if not (np.isfinite(left) and np.isfinite(right)):
        return None, None

    x = np.arange(a, b + 1, dtype=float)
    xa = float(a - 1)
    xb = float(b + 1)
    interp = left + ((x - xa) / (xb - xa)) * (right - left)
    orig = spec[a:b + 1].astype(float, copy=False)
    resid = orig - interp
    return interp, resid


def _detect_local_runs_from_resid(
    spec: np.ndarray,
    median_resid: np.ndarray,
    sigma: float,
    cfg: SpurConfig,
) -> list[tuple[int, int]]:
    if not np.isfinite(sigma) or sigma <= 0.0:
        return []

    candidate_idx = np.flatnonzero(np.abs(median_resid) > (cfg.nsigma * sigma))
    if candidate_idx.size == 0:
        return []

    accepted: list[tuple[int, int]] = []
    for a, b in _group_contiguous(candidate_idx):
        width = b - a + 1
        if width > cfg.max_width:
            continue
        _, seg_resid = _bridge_segment(spec, a, b)
        if seg_resid is None:
            continue
        finite = seg_resid[np.isfinite(seg_resid)]
        if finite.size == 0:
            continue
        score = float(np.nanmax(np.abs(finite)) / sigma)
        if score > cfg.nsigma:
            accepted.append((a, b))
    return accepted


def _runs_to_channels(runs: Sequence[tuple[int, int]]) -> np.ndarray:
    if not runs:
        return np.empty(0, dtype=int)
    parts = [np.arange(a, b + 1, dtype=int) for a, b in runs]
    return np.unique(np.concatenate(parts)).astype(int, copy=False)


def _repair_channels(spec: np.ndarray, channels: np.ndarray, cfg: SpurConfig) -> tuple[np.ndarray, np.ndarray]:
    out = np.asarray(spec, dtype=float).copy()
    ch = np.unique(np.asarray(channels, dtype=int))
    if ch.size == 0:
        return out, np.empty(0, dtype=int)

    repaired_runs = _group_contiguous(ch)
    repaired_channels: list[np.ndarray] = []
    for a, b in repaired_runs:
        width = b - a + 1
        if width > cfg.max_width:
            continue
        interp, _ = _bridge_segment(out, a, b)
        if interp is None:
            continue
        out[a:b + 1] = interp
        repaired_channels.append(np.arange(a, b + 1, dtype=int))

    if not repaired_channels:
        return out, np.empty(0, dtype=int)
    merged = np.unique(np.concatenate(repaired_channels)).astype(int, copy=False)
    return out, merged


def _as_object_array(rows: list[np.ndarray]) -> np.ndarray:
    out = np.empty(len(rows), dtype=object)
    for i, arr in enumerate(rows):
        out[i] = np.asarray(arr, dtype=np.int32)
    return out


def _promote_threshold(n_rows: int) -> int:
    return int(builtins.max(3, int(np.ceil(0.05 * float(n_rows)))))


def detect_and_repair_narrow_spurs(
    scantable: Scantable,
    *,
    spur_mode: str | None = "none",
    spur_nsigma: float = 7.0,
    spur_max_width: int = 2,
    verbose: bool = True,
) -> Scantable:
    """Detect and repair narrow 1-2 channel spurs on Ta* dumps.

    This is the Scantable-unmodified implementation. The output data are
    in-place repaired in a copied Scantable, while per-row bookkeeping is
    stored in table/history.
    """
    cfg = SpurConfig.normalize(mode=spur_mode, nsigma=spur_nsigma, max_width=spur_max_width)
    if cfg.mode == "none":
        return scantable.copy()

    if not isinstance(scantable.data, np.ndarray):
        raise TypeError("Narrow-spur repair currently supports only 2D numpy-array Scantable.data.")
    if scantable.data.ndim != 2:
        raise ValueError(f"Expected 2D Scantable.data; got shape={scantable.data.shape!r}.")

    sc = scantable.copy()
    data_in = np.asarray(sc.data, dtype=float)
    nrow, nchan = data_in.shape

    if len(sc.table) != nrow:
        raise ValueError(
            f"Row mismatch: len(table)={len(sc.table)} but data.shape[0]={nrow}."
        )

    median_resid_2d = _compute_local_median_residual_2d(data_in, half_window=2)
    sigma_per_row = _mad_sigma_rows(median_resid_2d)

    local_channels_per_row: list[np.ndarray] = []
    counts = np.zeros(nchan, dtype=np.int32)

    for i in range(nrow):
        runs = _detect_local_runs_from_resid(data_in[i], median_resid_2d[i], sigma_per_row[i], cfg)
        local_channels = _runs_to_channels(runs)
        local_channels_per_row.append(local_channels)
        if local_channels.size > 0:
            counts[local_channels] += 1

    finite_sigma = sigma_per_row[np.isfinite(sigma_per_row) & (sigma_per_row > 0.0)]
    fallback_sigma = float(np.median(finite_sigma)) if finite_sigma.size > 0 else np.nan
    if np.isfinite(fallback_sigma):
        bad_rows = np.flatnonzero((~np.isfinite(sigma_per_row)) | (sigma_per_row <= 0.0))
        for i in bad_rows:
            runs = _detect_local_runs_from_resid(data_in[i], median_resid_2d[i], fallback_sigma, cfg)
            sigma_per_row[i] = fallback_sigma
            local_channels_per_row[i] = _runs_to_channels(runs)

    counts[:] = 0
    for local_channels in local_channels_per_row:
        if local_channels.size > 0:
            counts[local_channels] += 1

    n_rows_local_detected = int(sum(1 for x in local_channels_per_row if len(x) > 0))
    n_local_detected = int(sum(len(x) for x in local_channels_per_row))

    promote_n = _promote_threshold(nrow)
    if cfg.mode == "promoted-all-dumps":
        promoted_channels = np.flatnonzero(counts >= promote_n).astype(int, copy=False)
    else:
        promoted_channels = np.empty(0, dtype=int)

    data_out = data_in.copy()
    repaired_channels_rows: list[np.ndarray] = []
    mode_rows: list[str] = []
    local_only_rows = np.zeros(nrow, dtype=np.int32)
    promoted_only_rows = np.zeros(nrow, dtype=np.int32)

    for i in range(nrow):
        local_channels = local_channels_per_row[i]
        repaired_spec, repaired_local = _repair_channels(data_in[i], local_channels, cfg)
        repaired_promoted = np.empty(0, dtype=int)

        if cfg.mode == "promoted-all-dumps" and promoted_channels.size > 0:
            promoted_only = np.setdiff1d(promoted_channels, repaired_local, assume_unique=False).astype(int, copy=False)
            repaired_spec, repaired_promoted = _repair_channels(repaired_spec, promoted_only, cfg)

        repaired_channels = np.union1d(repaired_local, repaired_promoted).astype(int, copy=False)
        data_out[i] = repaired_spec
        repaired_channels_rows.append(repaired_channels)

        local_only_rows[i] = int(repaired_local.size)
        promoted_only_rows[i] = int(repaired_promoted.size)

        if repaired_channels.size == 0:
            mode_rows.append("NONE")
        elif repaired_promoted.size > 0:
            mode_rows.append("PROMOTED")
        else:
            mode_rows.append("LOCAL")

    table = sc.table.copy()
    table["SPUR_N"] = np.asarray([len(x) for x in repaired_channels_rows], dtype=np.int32)
    table["SPUR_N_LOCAL"] = local_only_rows
    table["SPUR_N_PROM"] = promoted_only_rows
    table["SPUR_MODE"] = mode_rows
    table["SPUR_SIGMA"] = sigma_per_row
    table["SPUR_CHAN"] = _as_object_array(repaired_channels_rows)
    table["SPUR_DET_LOCAL"] = _as_object_array(local_channels_per_row)

    hist = dict(sc.history or {})
    hist["spur_enabled"] = True
    hist["spur_mode"] = cfg.mode
    hist["spur_nsigma"] = float(cfg.nsigma)
    hist["spur_max_width"] = int(cfg.max_width)
    hist["spur_promote_n"] = int(promote_n)
    hist["spur_local_detected_total"] = int(n_local_detected)
    hist["spur_rows_with_local_detections"] = int(n_rows_local_detected)
    hist["spur_counts_nonzero"] = json.dumps(
        {int(k): int(v) for k, v in enumerate(counts.tolist()) if v > 0},
        sort_keys=True,
    )
    hist["spur_promoted_channels"] = json.dumps([int(x) for x in promoted_channels.tolist()])
    hist["spur_total_repaired_channels"] = int(sum(len(x) for x in repaired_channels_rows))
    hist["spur_total_rows_with_repairs"] = int(sum(1 for x in repaired_channels_rows if len(x) > 0))

    orig_dtype = np.asarray(sc.data).dtype
    if np.issubdtype(orig_dtype, np.floating):
        sc.data = data_out.astype(orig_dtype, copy=False)
    else:
        sc.data = data_out
    sc.table = table
    sc.history = hist

    if verbose:
        n_rows_hit = int(sum(1 for x in repaired_channels_rows if len(x) > 0))
        n_rep = int(sum(len(x) for x in repaired_channels_rows))
        msg = (
            f"[SPUR] mode={cfg.mode}, nsigma={cfg.nsigma:g}, max_width={cfg.max_width}, "
            f"rows_detected={n_rows_local_detected}/{nrow}, detected_channels_total={n_local_detected}, "
            f"rows_repaired={n_rows_hit}/{nrow}, repaired_channels_total={n_rep}"
        )
        if cfg.mode == "promoted-all-dumps":
            msg += (
                f", promote_n={promote_n}, promoted_channel_count={int(promoted_channels.size)}, "
                f"promoted_channels={promoted_channels.tolist()}"
            )
        print(msg)

    return stamp_scantable_code_provenance(sc, stage="detect_and_repair_narrow_spurs")
