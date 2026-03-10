# src/sd_radio_spectral_fits/plotting/windowmask.py
from __future__ import annotations

import builtins
from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..ranges import parse_windows
from .utils import _interp_extrap, _process_spectrum, _rms_from_mean, recalculate_velocity_axis


# -----------------------------------------------------------------------------
# Fit window conversion helpers
# -----------------------------------------------------------------------------

def fit_windows_disp_vel(
    win_str: Optional[str],
    rest_hz_meta: Optional[float],
    rest_hz_user: Optional[float],
    *,
    shift_threshold_hz: float = 1.0,
) -> List[Tuple[float, float]]:
    """Parse BSL_WINF and return windows in *Display* velocity [km/s].

    Notes
    -----
    - BSL_WINF in the DB is treated as 'Native' velocity windows.
    - If user rest frequency differs from meta rest frequency, we shift window
      boundaries to the Display velocity frame using `recalculate_velocity_axis`.

    Returns
    -------
    list of (v0, v1)
        Window edges in Display velocity (km/s), each edge pair ordered (min,max).
    """
    if win_str is None:
        return []
    s = str(win_str).strip()
    if s.lower() in ("", "nan", "none"):
        return []

    raw_wins = [w.strip() for w in s.split(";") if w.strip()]
    if not raw_wins:
        return []

    do_shift = (
        rest_hz_meta is not None
        and rest_hz_user is not None
        and np.isfinite(rest_hz_meta)
        and np.isfinite(rest_hz_user)
        and abs(float(rest_hz_meta) - float(rest_hz_user)) > float(shift_threshold_hz)
    )

    out: List[Tuple[float, float]] = []
    for w0, w1 in parse_windows(raw_wins):
        v_edges = np.array([w0, w1], dtype=float)
        if do_shift:
            v_edges = recalculate_velocity_axis(v_edges, float(rest_hz_meta), float(rest_hz_user))
        lo = float(builtins.min(v_edges[0], v_edges[1]))
        hi = float(builtins.max(v_edges[0], v_edges[1]))
        out.append((lo, hi))
    return out


def vel_windows_to_xaxis(
    wins_vel: Sequence[Tuple[float, float]],
    xaxis: str,
    vel_disp: np.ndarray,
    freq_ghz: Optional[np.ndarray],
    chan: Optional[np.ndarray],
) -> List[Tuple[float, float]]:
    """Convert Display-velocity windows to current xaxis windows."""
    if not wins_vel:
        return []

    xa = str(xaxis).lower().strip() if xaxis is not None else "chan"
    if xa == "vel":
        return [(float(builtins.min(a, b)), float(builtins.max(a, b))) for (a, b) in wins_vel]

    vref = np.asarray(vel_disp, dtype=float)

    if xa == "freq" and freq_ghz is not None:
        xref = np.asarray(freq_ghz, dtype=float)
    else:
        if chan is None:
            xref = np.arange(vref.size, dtype=float)
        else:
            xref = np.asarray(chan, dtype=float)

    out: List[Tuple[float, float]] = []
    for v0, v1 in wins_vel:
        x_edges = _interp_extrap(vref, xref, np.array([v0, v1], dtype=float))
        lo = float(builtins.min(x_edges[0], x_edges[1]))
        hi = float(builtins.max(x_edges[0], x_edges[1]))
        out.append((lo, hi))
    return out


def fit_windows_xaxis(
    win_str: Optional[str],
    rest_hz_meta: Optional[float],
    rest_hz_user: Optional[float],
    xaxis: str,
    vel_disp: np.ndarray,
    freq_ghz: Optional[np.ndarray],
    chan: Optional[np.ndarray],
    *,
    shift_threshold_hz: float = 1.0,
) -> List[Tuple[float, float]]:
    """Convenience: BSL_WINF -> (Native vel) -> (Display vel) -> (current xaxis)."""
    wins_v = fit_windows_disp_vel(
        win_str,
        rest_hz_meta,
        rest_hz_user,
        shift_threshold_hz=shift_threshold_hz,
    )
    return vel_windows_to_xaxis(wins_v, xaxis, vel_disp, freq_ghz, chan)


def vel_range_to_xbounds(
    xaxis: str,
    vrange: Optional[Tuple[float, float]],
    vel_disp: np.ndarray,
    freq_ghz: Optional[np.ndarray],
    chan: Optional[np.ndarray],
) -> Optional[Tuple[float, float]]:
    """Convert a Display-velocity range (vmin, vmax) into (xmin, xmax) on current xaxis."""
    if vrange is None:
        return None
    v0, v1 = float(vrange[0]), float(vrange[1])
    wins = vel_windows_to_xaxis([(v0, v1)], xaxis, vel_disp, freq_ghz, chan)
    if not wins:
        return None
    return wins[0]


# -----------------------------------------------------------------------------
# Mask helpers
# -----------------------------------------------------------------------------

def mask_from_xbounds(x: np.ndarray, xbounds: Tuple[float, float]) -> np.ndarray:
    """Mask selecting x within [min(xbounds), max(xbounds)]."""
    lo = float(builtins.min(xbounds[0], xbounds[1]))
    hi = float(builtins.max(xbounds[0], xbounds[1]))
    xx = np.asarray(x, dtype=float)
    return (xx >= lo) & (xx <= hi)


def mask_from_windows(x: np.ndarray, wins_x: Sequence[Tuple[float, float]]) -> np.ndarray:
    """OR-mask of multiple x windows."""
    xx = np.asarray(x, dtype=float)
    m = np.zeros_like(xx, dtype=bool)
    for a, b in wins_x:
        lo = float(builtins.min(a, b))
        hi = float(builtins.max(a, b))
        m |= (xx >= lo) & (xx <= hi)
    return m


# -----------------------------------------------------------------------------
# RMS helpers
# -----------------------------------------------------------------------------

def prepare_rms_xy(
    x_raw: np.ndarray,
    y_raw: np.ndarray,
    rms_on: str,
    smooth_mode: str,
    smooth_width: int,
    box_downsample: bool,
    box_policy: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (x_for_rms, y_for_rms) on the grid used for RMS computation."""
    mode = str(rms_on).lower().strip() if rms_on is not None else "raw"
    if mode == "display":
        x_for_rms, y_for_rms = _process_spectrum(
            np.asarray(x_raw),
            np.asarray(y_raw),
            smooth_mode,
            smooth_width,
            box_downsample,
            box_policy,
        )
        return x_for_rms, y_for_rms
    return np.asarray(x_raw), np.asarray(y_raw)


def compute_rms_win(
    *,
    x_raw: np.ndarray,
    y_resid: np.ndarray,
    rms_on: str,
    smooth_mode: str,
    smooth_width: int,
    box_downsample: bool,
    box_policy: str,
    # axis / conversion context
    xaxis: str,
    vel_disp: np.ndarray,
    freq_ghz: Optional[np.ndarray],
    chan: Optional[np.ndarray],
    # selection inputs
    xr_vel: Optional[Tuple[float, float]],
    win_str: Optional[str],
    rest_hz_meta: Optional[float],
    rest_hz_user: Optional[float],
    # logic control
    has_bsl: bool,
    use_fit_windows: bool = True,
) -> Tuple[float, int, np.ndarray]:
    """Compute RMS on residual spectrum with consistent restfreq/window handling.

    Semantics (matches viewer/grid):
    - Always respect `xr_vel` (Display-velocity range).
    - If `has_bsl` and `win_str` is usable and `use_fit_windows` is True, compute
      RMS over (fit_windows & xr_vel). Otherwise compute RMS over xr_vel only.

    Returns
    -------
    rms : float
    n   : int
    use_mask : np.ndarray[bool]
        The final mask on the RMS grid.
    """
    x_for_rms, y_for_rms = prepare_rms_xy(
        x_raw,
        y_resid,
        rms_on,
        smooth_mode,
        smooth_width,
        box_downsample,
        box_policy,
    )

    mx_range = np.ones_like(y_for_rms, dtype=bool)
    if xr_vel is not None:
        xbounds = vel_range_to_xbounds(xaxis, xr_vel, vel_disp, freq_ghz, chan)
        if xbounds is not None:
            mx_range = mask_from_xbounds(x_for_rms, xbounds)

    if not use_fit_windows:
        use = mx_range
        rms = _rms_from_mean(y_for_rms[use])
        return rms, int(np.count_nonzero(use)), use

    wins_x = fit_windows_xaxis(
        win_str,
        rest_hz_meta,
        rest_hz_user,
        xaxis,
        vel_disp,
        freq_ghz,
        chan,
    )
    mask_win = mask_from_windows(x_for_rms, wins_x) if wins_x else np.zeros_like(y_for_rms, dtype=bool)

    use = (mask_win & mx_range) if (has_bsl and np.any(mask_win)) else mx_range
    rms = _rms_from_mean(y_for_rms[use])
    return rms, int(np.count_nonzero(use)), use
