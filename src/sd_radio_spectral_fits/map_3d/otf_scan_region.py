from __future__ import annotations

from dataclasses import dataclass
import builtins
from pathlib import Path
import warnings

import numpy as np
import pandas as pd


@dataclass
class OTFScanRegionResult:
    is_turn: np.ndarray
    effscan: np.ndarray
    scan_dir_deg: np.ndarray
    summary: dict


_MODE_SETTINGS = {
    "loose": {
        "motion_fraction": 0.10,
        "straightness_min": 0.18,
        "gap_fill": 2,
        "jump_factor": 12.0,
        "backtrack_fraction": 0.90,
        "turn_angle_deg": 70.0,
    },
    "auto": {
        "motion_fraction": 0.14,
        "straightness_min": 0.24,
        "gap_fill": 2,
        "jump_factor": 9.0,
        "backtrack_fraction": 0.70,
        "turn_angle_deg": 58.0,
    },
    "normal": {
        "motion_fraction": 0.18,
        "straightness_min": 0.30,
        "gap_fill": 2,
        "jump_factor": 7.0,
        "backtrack_fraction": 0.55,
        "turn_angle_deg": 48.0,
    },
    "strict": {
        "motion_fraction": 0.25,
        "straightness_min": 0.40,
        "gap_fill": 1,
        "jump_factor": 5.0,
        "backtrack_fraction": 0.35,
        "turn_angle_deg": 38.0,
    },
}

_FORMAL_STREAM_COLS = (
    "FDNUM",
    "IFNUM",
    "PLNUM",
)

_FALLBACK_STREAM_GROUP_COLS = (
    "BEAM",
    "POL",
    "IF",
    "ARRAYID",
    "FEED",
    "STREAM",
    "SIDEBAND",
    "BAND",
    "RX",
    "SAMPLER",
    "BOARD",
    "XFFTSBOARD",
)

_VALID_OTF_INPUT_STATES = (
    "with_turnarounds",
    "scan_only",
    "use_existing_labels",
)


def normalize_otf_input_state(input_state) -> str:
    text = "with_turnarounds" if input_state is None else str(input_state).strip().lower()
    if text not in _VALID_OTF_INPUT_STATES:
        raise ValueError(
            "otf_input_state must be one of 'with_turnarounds', 'scan_only', or 'use_existing_labels'."
        )
    return text


def _normalize_mode(mode) -> str:
    if mode is True:
        return "auto"
    if mode in (None, False):
        return "auto"
    text = str(mode).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return "auto"
    if text in _MODE_SETTINGS:
        return text
    raise ValueError(f"Unsupported otf_scan_region mode: {mode!r}")


def _normalize_existing_turn_labels(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text not in {"ignore", "respect"}:
        raise ValueError("existing_turn_labels must be 'ignore' or 'respect'.")
    return text


def _normalize_legacy_existing_is_turn_mode(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text == "prefer":
        return "respect"
    if text == "ignore":
        return "ignore"
    raise ValueError("otf_scan_existing_is_turn must be 'prefer' or 'ignore'.")


def _resolve_existing_turn_labels(*, existing_turn_labels=None, existing_is_turn_mode=None, input_state: str | None = None) -> str | None:
    state = normalize_otf_input_state(input_state)
    new_val = _normalize_existing_turn_labels(existing_turn_labels)
    legacy_given = existing_is_turn_mode is not None
    legacy_val = _normalize_legacy_existing_is_turn_mode(existing_is_turn_mode)
    if state == "use_existing_labels":
        if new_val is not None or legacy_given:
            raise ValueError(
                "otf_input_state='use_existing_labels' does not accept existing_turn_labels or otf_scan_existing_is_turn. "
                "IS_TURN is required and used directly in this mode."
            )
        return None
    if new_val is None and not legacy_given:
        return "ignore"
    if new_val is None:
        return legacy_val
    if legacy_given and legacy_val != new_val:
        raise ValueError(
            "existing_turn_labels conflicts with deprecated otf_scan_existing_is_turn."
        )
    return new_val


def _sorted_group_indices(group_indices: np.ndarray, time_arr: np.ndarray) -> np.ndarray:
    if group_indices.size <= 1:
        return np.asarray(group_indices, dtype=np.int64)
    t = np.asarray(time_arr[group_indices], dtype=float)
    order_local = np.arange(group_indices.size, dtype=np.int64)
    finite = np.isfinite(t)
    if np.count_nonzero(finite) >= 2 and np.unique(t[finite]).size >= 2:
        t_key = np.where(finite, t, np.inf)
        order = np.lexsort((order_local, t_key))
        return np.asarray(group_indices[order], dtype=np.int64)
    return np.asarray(group_indices, dtype=np.int64)


def _rolling_sum_centered(arr: np.ndarray, half_window: int) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    if a.size == 0:
        return np.array([], dtype=float)
    if half_window <= 0:
        return a.copy()
    half_window = int(min(int(half_window), max(0, (a.size - 1) // 2)))
    if half_window <= 0:
        return a.copy()
    full = np.convolve(a, np.ones(2 * half_window + 1, dtype=float), mode="same")
    return np.asarray(full, dtype=float)


def _fill_small_false_gaps(mask: np.ndarray, allowed: np.ndarray, max_gap: int) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    allowed = np.asarray(allowed, dtype=bool)
    if max_gap <= 0 or out.size == 0:
        return out
    i = 0
    n = out.size
    while i < n:
        if out[i]:
            i += 1
            continue
        j = i
        while j < n and (not out[j]):
            j += 1
        if i > 0 and j < n and out[i - 1] and out[j] and (j - i) <= max_gap and np.all(allowed[i:j]):
            out[i:j] = True
        i = j
    return out


def _iter_true_runs(mask: np.ndarray):
    m = np.asarray(mask, dtype=bool)
    n = m.size
    i = 0
    while i < n:
        while i < n and not m[i]:
            i += 1
        if i >= n:
            break
        j = i + 1
        while j < n and m[j]:
            j += 1
        yield i, j
        i = j


def _iter_segment_slices_from_breaks(run_size: int, breaks: np.ndarray):
    """Yield contiguous [start, stop) slices split at break positions, without dropping points.

    breaks[j] == True means sample j starts a new segment. The boundary sample stays in the
    later segment instead of being discarded.
    """
    n = int(run_size)
    if n <= 0:
        return
    br = np.asarray(breaks, dtype=bool)
    if br.shape[0] != n:
        raise ValueError('breaks must have the same length as the run.')
    start = 0
    for j in range(1, n):
        if br[j]:
            if start < j:
                yield start, j
            start = j
    if start < n:
        yield start, n


def _local_tangent_vectors(xs: np.ndarray, ys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    n = xs.size
    tx = np.full(n, np.nan, dtype=float)
    ty = np.full(n, np.nan, dtype=float)
    if n < 2:
        return tx, ty
    tx[0] = xs[1] - xs[0]
    ty[0] = ys[1] - ys[0]
    tx[-1] = xs[-1] - xs[-2]
    ty[-1] = ys[-1] - ys[-2]
    if n > 2:
        tx[1:-1] = xs[2:] - xs[:-2]
        ty[1:-1] = ys[2:] - ys[:-2]
    return tx, ty

def _dominant_direction(xs: np.ndarray, ys: np.ndarray) -> np.ndarray | None:
    if xs.size < 2:
        return None
    dx = float(xs[-1] - xs[0])
    dy = float(ys[-1] - ys[0])
    norm = float(np.hypot(dx, dy))
    if np.isfinite(norm) and norm > 0:
        return np.array([dx / norm, dy / norm], dtype=float)
    arr = np.column_stack([xs - np.nanmean(xs), ys - np.nanmean(ys)])
    if arr.shape[0] < 2 or not np.isfinite(arr).all():
        return None
    try:
        _, _, vh = np.linalg.svd(arr, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    vec = np.asarray(vh[0], dtype=float)
    norm = float(np.hypot(vec[0], vec[1]))
    if not np.isfinite(norm) or norm <= 0:
        return None
    return vec / norm


def _direction_to_angle_deg(vec: np.ndarray | None) -> float:
    if vec is None:
        return np.nan
    ang = float(np.degrees(np.arctan2(vec[1], vec[0])))
    ang = np.mod(ang, 180.0)
    if ang < 0:
        ang += 180.0
    return ang


def _coerce_group_value(value):
    if isinstance(value, np.generic):
        value = value.item()
    if pd.isna(value):
        return "<NA>"
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8", errors="replace")
        except Exception:
            return repr(value)
    if isinstance(value, (int, float, str, bool)):
        return value
    return str(value)


def _candidate_group_columns(table: pd.DataFrame) -> tuple[list[str], list[str], list[str], list[str]]:
    cols = list(getattr(table, "columns", []))
    formal_present = [c for c in _FORMAL_STREAM_COLS if c in cols]
    formal_missing = [c for c in _FORMAL_STREAM_COLS if c not in cols]
    fallback_cols: list[str] = []
    if formal_missing:
        fallback_cols = [c for c in _FALLBACK_STREAM_GROUP_COLS if c in cols and c not in formal_present]
    group_cols = formal_present + fallback_cols
    return group_cols, formal_present, formal_missing, fallback_cols


def _build_groups(table: pd.DataFrame, base_scan_id: np.ndarray | None) -> tuple[dict[tuple, list[int]], dict]:
    n = len(table)
    group_cols, formal_present, formal_missing, fallback_cols = _candidate_group_columns(table)
    if formal_missing:
        if fallback_cols:
            warnings.warn(
                'Formal OTF stream columns FDNUM/IFNUM/PLNUM are incomplete; falling back to '                 + ', '.join(str(c) for c in fallback_cols) + '.',
                RuntimeWarning,
                stacklevel=2,
            )
        else:
            warnings.warn(
                'Formal OTF stream columns FDNUM/IFNUM/PLNUM are incomplete and no fallback stream key exists; '
                'grouping will use SCAN only.',
                RuntimeWarning,
                stacklevel=2,
            )
    groups: dict[tuple, list[int]] = {}
    if base_scan_id is None:
        base = np.zeros(n, dtype=float)
        finite_scan = np.ones(n, dtype=bool)
    else:
        base = np.asarray(base_scan_id, dtype=float)
        finite_scan = np.isfinite(base)
    for i in range(n):
        if not finite_scan[i]:
            continue
        key = [int(round(float(base[i])))]
        for col in group_cols:
            key.append(_coerce_group_value(table.iloc[i][col]))
        groups.setdefault(tuple(key), []).append(i)
    info = {
        "stream_group_columns": list(group_cols),
        "group_key_columns": ["SCAN"] + list(group_cols),
        "formal_stream_columns_present": list(formal_present),
        "missing_formal_stream_columns": list(formal_missing),
        "fallback_stream_columns": list(fallback_cols),
        "used_fallback_stream_key": bool(fallback_cols),
        "stream_group_strategy": (
            "formal_only" if not formal_missing and formal_present else
            "formal_plus_fallback" if formal_present or fallback_cols else
            "scan_only"
        ),
    }
    return groups, info

def _build_existing_label_stream_groups(table: pd.DataFrame) -> tuple[dict[tuple, list[int]], dict]:
    cols = list(getattr(table, "columns", []))
    missing_formal = [c for c in _FORMAL_STREAM_COLS if c not in cols]
    if missing_formal:
        raise ValueError(
            "otf_input_state='use_existing_labels' requires formal stream columns FDNUM/IFNUM/PLNUM for stream separation. "
            f"Missing columns: {', '.join(str(c) for c in missing_formal)}."
        )
    groups: dict[tuple, list[int]] = {}
    for i in range(len(table)):
        key = tuple(_coerce_group_value(table.iloc[i][col]) for col in _FORMAL_STREAM_COLS)
        groups.setdefault(key, []).append(i)
    info = {
        "stream_group_columns": list(_FORMAL_STREAM_COLS),
        "group_key_columns": list(_FORMAL_STREAM_COLS) + ["SCAN"],
        "formal_stream_columns_present": list(_FORMAL_STREAM_COLS),
        "missing_formal_stream_columns": [],
        "fallback_stream_columns": [],
        "used_fallback_stream_key": False,
        "stream_group_strategy": "formal_only",
    }
    return groups, info


def _unique_ints_in_order(values: np.ndarray) -> list[int]:
    seen = set()
    out: list[int] = []
    for v in np.asarray(values, dtype=float):
        if not np.isfinite(v):
            continue
        iv = int(round(float(v)))
        if iv not in seen:
            seen.add(iv)
            out.append(iv)
    return out


def _write_region_png(
    png_path: str,
    *,
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    local_motion: np.ndarray,
    accepted_mask: np.ndarray,
    final_is_turn: np.ndarray,
    effscan: np.ndarray,
    threshold: float,
):
    try:
        import matplotlib.lines as mlines
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.figure import Figure
        from matplotlib.cm import get_cmap
    except Exception as e:
        warnings.warn(f'Failed to import matplotlib for OTF scan-region PNG output: {e}', RuntimeWarning, stacklevel=2)
        return

    path = Path(png_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    x = np.asarray(x_arcsec, dtype=float)
    y = np.asarray(y_arcsec, dtype=float)
    accepted = np.asarray(accepted_mask, dtype=bool)
    rejected = np.asarray(final_is_turn, dtype=bool)
    eff = np.asarray(effscan, dtype=int)
    idx_all = np.arange(x.size, dtype=int)

    fig = Figure(figsize=(14, 7), constrained_layout=True)
    FigureCanvas(fig)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 1.0])
    ax_path = fig.add_subplot(gs[:, 0])
    ax_x = fig.add_subplot(gs[0, 1])
    ax_y = fig.add_subplot(gs[1, 1], sharex=ax_x)

    finite_xy = np.isfinite(x) & np.isfinite(y)
    finite_acc = finite_xy & accepted
    finite_rej = finite_xy & rejected

    # Full path first, so the accepted scan region stands out clearly.
    if np.any(finite_xy):
        ax_path.plot(x[finite_xy], y[finite_xy], linewidth=0.7, alpha=0.35, color='0.75', zorder=1)

    accepted_ids = np.unique(eff[eff >= 0])
    n_runs = int(accepted_ids.size)
    if n_runs > 0:
        cmap = get_cmap('turbo')
        label_every = max(1, int(np.ceil(n_runs / 80.0)))
        for kk, sid in enumerate(accepted_ids):
            idx = np.flatnonzero(eff == sid)
            if idx.size == 0:
                continue
            color = cmap(kk / max(1, n_runs - 1))
            ax_path.plot(x[idx], y[idx], linewidth=2.4, alpha=0.95, color=color, zorder=3)
            ax_path.scatter(x[idx], y[idx], s=7, color=color, alpha=0.95, zorder=4)
            if (kk % label_every) == 0 or kk == 0 or kk == n_runs - 1:
                mid = idx[idx.size // 2]
                if np.isfinite(x[mid]) and np.isfinite(y[mid]):
                    ax_path.text(
                        x[mid],
                        y[mid],
                        str(int(sid)),
                        fontsize=6,
                        ha='center',
                        va='center',
                        bbox=dict(boxstyle='round,pad=0.12', facecolor='white', edgecolor='none', alpha=0.75),
                        zorder=5,
                    )
    else:
        label_every = 1

    if np.any(finite_rej):
        ax_path.scatter(
            x[finite_rej],
            y[finite_rej],
            s=12,
            marker='x',
            linewidths=0.8,
            color='tab:blue',
            alpha=0.85,
            zorder=2,
        )

    ax_path.set_title(f"Accepted OTF scan region (runs={n_runs})")
    ax_path.set_xlabel("X [arcsec]")
    ax_path.set_ylabel("Y [arcsec]")
    ax_path.set_aspect("equal", adjustable="box")
    ax_path.grid(alpha=0.18)

    text_lines = [
        "Gray thin line: full path",
        "Colored thick lines: accepted scan runs",
        "Blue x: rejected / turnaround / line-change points",
    ]
    if n_runs > 80:
        text_lines.append(f"Run labels shown every {label_every} runs")
    ax_path.text(
        0.015,
        0.985,
        "\n".join(text_lines),
        transform=ax_path.transAxes,
        ha='left',
        va='top',
        fontsize=8,
        bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='0.8', alpha=0.92),
        zorder=6,
    )

    legend_handles = [
        mlines.Line2D([], [], color='0.75', linewidth=1.0, label='Full path'),
        mlines.Line2D([], [], color='tab:red', linewidth=2.4, label='Accepted scan run'),
        mlines.Line2D([], [], color='tab:blue', marker='x', linestyle='None', markersize=6, label='Rejected / turnaround'),
    ]
    ax_path.legend(handles=legend_handles, loc='lower left', fontsize=8, framealpha=0.92)

    def _shade_accepted_spans(ax):
        if not np.any(accepted):
            return
        start = None
        for i, flag in enumerate(accepted):
            if flag and start is None:
                start = i
            elif (not flag) and start is not None:
                ax.axvspan(start, i - 1, alpha=0.16, color='tab:green')
                start = None
        if start is not None:
            ax.axvspan(start, len(accepted) - 1, alpha=0.16, color='tab:green')

    finite_x = np.isfinite(x)
    finite_y = np.isfinite(y)
    ax_x.plot(idx_all[finite_x], x[finite_x], linewidth=1.0, color='black')
    _shade_accepted_spans(ax_x)
    ax_x.set_title("Row-order diagnostic: X during accepted spans")
    ax_x.set_ylabel("X [arcsec]")
    ax_x.grid(alpha=0.18)

    ax_y.plot(idx_all[finite_y], y[finite_y], linewidth=1.0, color='black')
    _shade_accepted_spans(ax_y)
    ax_y.set_title("Row-order diagnostic: Y during accepted spans")
    ax_y.set_xlabel("Row index")
    ax_y.set_ylabel("Y [arcsec]")
    ax_y.grid(alpha=0.18)

    if np.isfinite(threshold):
        finite_lm = np.isfinite(local_motion)
        lm_pos = np.asarray(local_motion[finite_lm], dtype=float)
        if lm_pos.size:
            txt = f"Local-motion threshold used internally: {threshold:.1f} arcsec"
            ax_y.text(
                0.99,
                0.02,
                txt,
                transform=ax_y.transAxes,
                ha='right',
                va='bottom',
                fontsize=8,
                bbox=dict(boxstyle='round,pad=0.18', facecolor='white', edgecolor='0.85', alpha=0.9),
            )

    try:
        fig.savefig(path, dpi=160)
    except Exception as e:
        warnings.warn(f'Failed to write OTF scan-region PNG to {path}: {e}', RuntimeWarning, stacklevel=2)


def identify_otf_scan_regions(
    table,
    *,
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    flag_on: np.ndarray,
    time_arr: np.ndarray,
    base_scan_id: np.ndarray | None,
    existing_is_turn: np.ndarray | None = None,
    mode: str = "auto",
    existing_turn_labels: str | None = None,
    existing_is_turn_mode: str | None = None,
    input_state: str = "with_turnarounds",
    png_path: str | None = None,
) -> OTFScanRegionResult:
    table = pd.DataFrame(table)
    x = np.asarray(x_arcsec, dtype=float)
    y = np.asarray(y_arcsec, dtype=float)
    flag_on = np.asarray(flag_on, dtype=bool)
    time_arr = np.asarray(time_arr, dtype=float)
    n = len(table)
    if x.shape[0] != n or y.shape[0] != n or flag_on.shape[0] != n or time_arr.shape[0] != n:
        raise ValueError("x_arcsec, y_arcsec, flag_on, time_arr, and table must have the same length.")
    if base_scan_id is not None and np.asarray(base_scan_id).shape[0] != n:
        raise ValueError("base_scan_id must have the same length as table.")
    if existing_is_turn is not None and np.asarray(existing_is_turn).shape[0] != n:
        raise ValueError("existing_is_turn must have the same length as table.")

    mode_norm = _normalize_mode(mode)
    input_state_norm = normalize_otf_input_state(input_state)
    existing_turn_policy = _resolve_existing_turn_labels(
        existing_turn_labels=existing_turn_labels,
        existing_is_turn_mode=existing_is_turn_mode,
        input_state=input_state_norm,
    )
    settings = _MODE_SETTINGS[mode_norm]

    finite_xy = np.isfinite(x) & np.isfinite(y)
    existing_turn = np.zeros(n, dtype=bool)
    if existing_is_turn is not None and existing_turn_policy == "respect":
        existing_turn = np.asarray(existing_is_turn, dtype=bool).copy()

    if input_state_norm == "use_existing_labels":
        if base_scan_id is None:
            raise ValueError("otf_input_state='use_existing_labels' requires base_scan_id / SCAN labels.")
        if existing_is_turn is None:
            raise ValueError("otf_input_state='use_existing_labels' requires existing_is_turn / IS_TURN labels.")
        scan_arr = np.asarray(base_scan_id, dtype=float)
        turn_arr = np.asarray(existing_is_turn, dtype=bool)
        groups, group_info = _build_existing_label_stream_groups(table)
        accepted = finite_xy & flag_on & np.isfinite(scan_arr) & (~turn_arr)
        final_is_turn = np.asarray(turn_arr, dtype=bool).copy()
        effscan = np.full(n, -1, dtype=np.int64)
        next_global = 0
        for _, idx_list in groups.items():
            idx = np.asarray(idx_list, dtype=np.int64)
            usable = idx[accepted[idx]]
            if usable.size == 0:
                continue
            local_ids = _unique_ints_in_order(scan_arr[usable])
            mapping = {sid: int(next_global + j) for j, sid in enumerate(local_ids)}
            for row in usable:
                effscan[int(row)] = mapping[int(round(float(scan_arr[int(row)])))]
            next_global += len(local_ids)
        scan_dir = np.full(n, np.nan, dtype=float)
        accepted_scan_ids = np.unique(effscan[effscan >= 0]).astype(np.int64) if np.any(effscan >= 0) else np.array([], dtype=np.int64)
        summary = {
            "summary_version": 2,
            "kind": "single_table",
            "reference_scan_id_semantics": "local for a single GridInput; global after merged basketweave inputs",
            "mode": None,
            "input_state": str(input_state_norm),
            "existing_turn_labels": None,
            "existing_is_turn_mode": None,
            **group_info,
            "num_groups": int(len(groups)),
            "num_runs": int(next_global),
            "num_points_total": int(n),
            "num_points_flag_on": int(np.count_nonzero(flag_on)),
            "num_points_existing_turn": int(np.count_nonzero(turn_arr)),
            "num_points_scan": int(np.count_nonzero(accepted)),
            "num_points_turn": int(np.count_nonzero(final_is_turn)),
            "num_points_unresolved": int(np.count_nonzero(np.isfinite(scan_arr) & (~turn_arr) & (~accepted))),
            "num_points_excluded": int(np.count_nonzero(~accepted)),
            "has_existing_is_turn": True,
            "accepted_scan_id_min": None if accepted_scan_ids.size == 0 else int(accepted_scan_ids.min()),
            "accepted_scan_id_max": None if accepted_scan_ids.size == 0 else int(accepted_scan_ids.max()),
            "accepted_scan_id_count": int(accepted_scan_ids.size),
            "geometry_mask_applied": False,
        }
        if png_path:
            _write_region_png(
                str(png_path),
                x_arcsec=x,
                y_arcsec=y,
                local_motion=np.full(n, np.nan, dtype=float),
                accepted_mask=accepted,
                final_is_turn=final_is_turn,
                effscan=effscan,
                threshold=np.nan,
            )
        return OTFScanRegionResult(
            is_turn=np.asarray(final_is_turn, dtype=bool),
            effscan=np.asarray(effscan, dtype=np.int64),
            scan_dir_deg=np.asarray(scan_dir, dtype=float),
            summary=summary,
        )

    groups, group_info = _build_groups(table, base_scan_id)

    final_is_turn = np.ones(n, dtype=bool)
    final_is_turn[existing_turn] = True
    effscan = np.full(n, -1, dtype=np.int64)
    scan_dir = np.full(n, np.nan, dtype=float)
    local_motion_global = np.full(n, np.nan, dtype=float)

    next_effscan = 0
    accepted_runs = 0
    accepted_points = 0

    for _, group_list in groups.items():
        group_idx = np.asarray(group_list, dtype=np.int64)
        if group_idx.size == 0:
            continue
        ordered = _sorted_group_indices(group_idx, time_arr)
        xs = x[ordered]
        ys = y[ordered]
        finite_group = finite_xy[ordered]
        candidate = finite_group & flag_on[ordered]
        if existing_turn_policy == "respect":
            candidate &= ~existing_turn[ordered]
        if np.count_nonzero(candidate) < 2:
            continue

        time_group = np.asarray(time_arr[ordered], dtype=float)
        dt_days = np.diff(time_group)
        dt_sec = dt_days * 86400.0
        dt_pos = dt_sec[np.isfinite(dt_sec) & (dt_sec > 0)]
        if dt_pos.size:
            half_window = int(np.clip(np.rint(0.75 / np.nanmedian(dt_pos)), 1, 8))
        else:
            half_window = 2

        step = np.zeros(ordered.size, dtype=float)
        step[1:] = np.hypot(np.diff(xs), np.diff(ys))
        step_valid = np.zeros_like(step)
        step_valid[1:] = step[1:] * (candidate[1:] & candidate[:-1])
        local_motion = _rolling_sum_centered(step_valid, half_window)
        local_motion_global[ordered] = local_motion

        positive_motion = local_motion[candidate & np.isfinite(local_motion) & (local_motion > 0)]
        move_threshold = np.nan
        if input_state_norm == "scan_only":
            moving = candidate.copy()
        else:
            if positive_motion.size == 0:
                continue
            q90 = float(np.nanpercentile(positive_motion, 90.0))
            if not np.isfinite(q90) or q90 <= 0:
                continue
            move_threshold = float(settings["motion_fraction"] * q90)
            moving = candidate & np.isfinite(local_motion) & (local_motion >= move_threshold)
            moving = _fill_small_false_gaps(moving, candidate, int(settings["gap_fill"]))

        positive_steps = step_valid[(step_valid > 0) & np.isfinite(step_valid)]
        median_step = float(np.nanmedian(positive_steps)) if positive_steps.size else np.nan
        if not np.isfinite(median_step) or median_step <= 0:
            median_step = 1.0
        min_points_per_run = int(max(5, 2 * half_window + 1))
        min_path_length = float(max(10.0, 3.0 * median_step))

        segmented_runs: list[np.ndarray] = []
        for start, stop in list(_iter_true_runs(moving)):
            run_idx = np.arange(start, stop, dtype=np.int64)
            if run_idx.size < 2:
                segmented_runs.append(run_idx)
                continue
            u = _dominant_direction(xs[run_idx], ys[run_idx])
            if u is None:
                segmented_runs.append(run_idx)
                continue
            dx = np.diff(xs[run_idx])
            dy = np.diff(ys[run_idx])
            step_run = np.hypot(dx, dy)
            positive_step_run = step_run[(step_run > 0) & np.isfinite(step_run)]
            run_step_med = float(np.nanmedian(positive_step_run)) if positive_step_run.size else median_step
            if not np.isfinite(run_step_med) or run_step_med <= 0:
                run_step_med = median_step
            tol = float(settings["backtrack_fraction"] * run_step_med)
            jump_limit = float(settings["jump_factor"] * run_step_med)
            breaks = np.zeros(run_idx.size, dtype=bool)
            if run_idx.size > 1:
                breaks[1:] = (step_run > jump_limit)

            tx, ty = _local_tangent_vectors(xs[run_idx], ys[run_idx])
            tnorm = np.hypot(tx, ty)
            tnorm_safe = np.where(tnorm > 0, tnorm, np.nan)
            tux = tx / tnorm_safe
            tuy = ty / tnorm_safe
            turn_limit = float(np.cos(np.deg2rad(float(settings["turn_angle_deg"]))))
            anchor_sum = None
            anchor_count = 0
            for jj in range(run_idx.size):
                if not np.isfinite(tux[jj]) or not np.isfinite(tuy[jj]):
                    continue
                cur_u = np.array([tux[jj], tuy[jj]], dtype=float)
                if jj > 0 and breaks[jj]:
                    anchor_sum = cur_u.copy()
                    anchor_count = 1
                    continue
                if anchor_sum is None:
                    anchor_sum = cur_u.copy()
                    anchor_count = 1
                    continue
                anchor_norm = float(np.hypot(anchor_sum[0], anchor_sum[1]))
                if not np.isfinite(anchor_norm) or anchor_norm <= 0:
                    anchor_sum = cur_u.copy()
                    anchor_count = 1
                    continue
                anchor_u = anchor_sum / anchor_norm
                dot = float(np.clip(cur_u[0] * anchor_u[0] + cur_u[1] * anchor_u[1], -1.0, 1.0))
                if jj > 0 and anchor_count >= 3 and dot < turn_limit:
                    breaks[jj] = True
                    anchor_sum = cur_u.copy()
                    anchor_count = 1
                else:
                    if anchor_count < 4:
                        anchor_sum = anchor_sum + cur_u
                        anchor_count += 1

            for seg_start, seg_stop in _iter_segment_slices_from_breaks(run_idx.size, breaks):
                segmented_runs.append(run_idx[seg_start:seg_stop])

        for run_idx in segmented_runs:
            if run_idx.size == 0:
                continue
            if input_state_norm == "scan_only":
                if run_idx.size >= 2:
                    u_run = _dominant_direction(xs[run_idx], ys[run_idx])
                else:
                    u_run = None
            else:
                if run_idx.size < min_points_per_run:
                    continue
                step_run = np.zeros(run_idx.size, dtype=float)
                if run_idx.size > 1:
                    step_run[1:] = np.hypot(np.diff(xs[run_idx]), np.diff(ys[run_idx]))
                path_length = float(np.nansum(step_run))
                if not np.isfinite(path_length) or path_length < min_path_length:
                    continue
                net_disp = float(np.hypot(xs[run_idx[-1]] - xs[run_idx[0]], ys[run_idx[-1]] - ys[run_idx[0]]))
                straightness = net_disp / path_length if path_length > 0 else 0.0
                if not np.isfinite(straightness) or straightness < float(settings["straightness_min"]):
                    continue
                u_run = _dominant_direction(xs[run_idx], ys[run_idx])
                if u_run is None:
                    continue
                if run_idx.size > 1:
                    dxr = np.diff(xs[run_idx])
                    dyr = np.diff(ys[run_idx])
                    along_run = dxr * u_run[0] + dyr * u_run[1]
                    finite_along = np.isfinite(along_run)
                    if np.any(finite_along):
                        neg_fraction = float(np.count_nonzero(along_run[finite_along] < -tol) / np.count_nonzero(finite_along))
                        if neg_fraction > 0.25:
                            continue
            out_idx = ordered[run_idx]
            final_is_turn[out_idx] = False
            effscan[out_idx] = next_effscan
            scan_dir[out_idx] = _direction_to_angle_deg(u_run)
            next_effscan += 1
            accepted_runs += 1
            accepted_points += int(out_idx.size)

    accepted_mask = effscan >= 0
    final_is_turn[accepted_mask] = False
    if existing_turn_policy == "respect":
        final_is_turn |= existing_turn
        effscan[existing_turn] = -1
        scan_dir[existing_turn] = np.nan

    finite_scan_mask = np.isfinite(np.asarray(base_scan_id, dtype=float)) if base_scan_id is not None else np.ones(n, dtype=bool)
    unresolved = finite_xy & flag_on & finite_scan_mask & (~accepted_mask) & (~final_is_turn)
    final_is_turn[unresolved] = True

    accepted_scan_ids = np.unique(effscan[effscan >= 0]).astype(np.int64)
    summary = {
        "summary_version": 2,
        "kind": "single_table",
        "reference_scan_id_semantics": "local for a single GridInput; global after merged basketweave inputs",
        "mode": mode_norm,
        "input_state": input_state_norm,
        "existing_turn_labels": existing_turn_policy,
        "existing_is_turn_mode": None if existing_turn_policy is None else ("prefer" if existing_turn_policy == "respect" else "ignore"),
        **group_info,
        "num_groups": int(len(groups)),
        "num_runs": int(accepted_runs),
        "num_points_total": int(n),
        "num_points_flag_on": int(np.count_nonzero(flag_on)),
        "num_points_existing_turn": int(np.count_nonzero(existing_turn)),
        "num_points_scan": int(accepted_points),
        "num_points_turn": int(np.count_nonzero(final_is_turn)),
        "num_points_unresolved": int(np.count_nonzero(unresolved)),
        "num_points_excluded": int(np.count_nonzero(~accepted_mask)),
        "has_existing_is_turn": bool(existing_is_turn is not None),
        "accepted_scan_id_min": None if accepted_scan_ids.size == 0 else int(accepted_scan_ids.min()),
        "accepted_scan_id_max": None if accepted_scan_ids.size == 0 else int(accepted_scan_ids.max()),
        "accepted_scan_id_count": int(accepted_scan_ids.size),
        "geometry_mask_applied": bool(input_state_norm == "with_turnarounds"),
    }

    if png_path:
        finite_motion = local_motion_global[np.isfinite(local_motion_global) & (local_motion_global > 0)]
        if input_state_norm == "scan_only":
            threshold_for_plot = np.nan
        else:
            threshold_for_plot = float(np.nanpercentile(finite_motion, 90.0) * settings["motion_fraction"]) if finite_motion.size else np.nan
        _write_region_png(
            str(png_path),
            x_arcsec=x,
            y_arcsec=y,
            local_motion=local_motion_global,
            accepted_mask=accepted_mask,
            final_is_turn=final_is_turn,
            effscan=effscan,
            threshold=threshold_for_plot,
        )

    return OTFScanRegionResult(
        is_turn=np.asarray(final_is_turn, dtype=bool),
        effscan=np.asarray(effscan, dtype=np.int64),
        scan_dir_deg=np.asarray(scan_dir, dtype=float),
        summary=summary,
    )


def format_otf_scan_summary(summary: dict | None) -> str:
    if not summary:
        return 'otf_scan_region: summary unavailable'
    kind = str(summary.get('kind', 'single_table'))
    if kind == 'multi_input':
        n_inputs = summary.get('num_inputs', '?')
        runs = summary.get('num_runs_total', '?')
        kept = summary.get('num_points_scan_total', '?')
        total = summary.get('num_points_total', '?')
        ranges = summary.get('global_scan_id_ranges') or []
        if ranges:
            mins = [r.get('global_scan_id_min') for r in ranges if r.get('global_scan_id_min') is not None]
            maxs = [r.get('global_scan_id_max') for r in ranges if r.get('global_scan_id_max') is not None]
            if mins and maxs:
                return f'otf_scan_region: inputs={n_inputs} runs={runs} kept={kept}/{total} global_scan_ids={builtins.min(mins)}-{builtins.max(maxs)}'
        return f'otf_scan_region: inputs={n_inputs} runs={runs} kept={kept}/{total}'
    groups = summary.get('num_groups', '?')
    runs = summary.get('num_runs', '?')
    kept = summary.get('num_points_scan', '?')
    total = summary.get('num_points_total', '?')
    strategy = summary.get('stream_group_strategy', 'unknown')
    fallback_cols = summary.get('fallback_stream_columns') or []
    scan_min = summary.get('accepted_scan_id_min')
    scan_max = summary.get('accepted_scan_id_max')
    text = f'otf_scan_region: groups={groups} runs={runs} kept={kept}/{total} strategy={strategy}'
    if scan_min is not None and scan_max is not None:
        text += f' scan_ids={scan_min}-{scan_max}'
    if fallback_cols:
        text += ' fallback=' + ','.join(str(v) for v in fallback_cols)
    return text
