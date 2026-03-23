from __future__ import annotations

import numpy as np
from typing import Tuple


def _as_float_array(value):
    return np.asarray(value, dtype=float)


def _rotation_theta_deg(boresight_el_deg, rotation_sign: float):
    return float(rotation_sign) * _as_float_array(boresight_el_deg)


def rotate_el0_offset(offset_x_el0_arcsec, offset_y_el0_arcsec, boresight_el_deg, rotation_sign) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rotate projected El=0 beam offsets by sigma * El.

    Returns
    -------
    beam_dx_arcsec, beam_dy_arcsec, theta_deg
        Projected beam offset on the sky where
        x = dAz * cos(El), y = dEl.
    """
    theta_deg = _rotation_theta_deg(boresight_el_deg, rotation_sign)
    theta = np.deg2rad(theta_deg)
    dx0 = _as_float_array(offset_x_el0_arcsec)
    dy0 = _as_float_array(offset_y_el0_arcsec)
    dx = dx0 * np.cos(theta) - dy0 * np.sin(theta)
    dy = dx0 * np.sin(theta) + dy0 * np.cos(theta)
    return _as_float_array(dx), _as_float_array(dy), _as_float_array(theta_deg)


def offset_xy_to_beam_azel(boresight_az_deg, boresight_el_deg, beam_dx_arcsec, beam_dy_arcsec) -> Tuple[np.ndarray, np.ndarray]:
    """Convert projected offsets (x=dAz*cosEl, y=dEl) to beam-center Az/El."""
    az0 = _as_float_array(boresight_az_deg)
    el0 = _as_float_array(boresight_el_deg)
    dx = _as_float_array(beam_dx_arcsec)
    dy = _as_float_array(beam_dy_arcsec)

    cos_el = np.cos(np.deg2rad(el0))
    tiny = np.cos(np.deg2rad(89.9))
    safe = np.abs(cos_el) >= tiny

    d_az_deg = np.full_like(el0, np.nan, dtype=float)
    d_az_deg[safe] = dx[safe] / (3600.0 * cos_el[safe])
    d_el_deg = dy / 3600.0
    return (az0 + d_az_deg) % 360.0, el0 + d_el_deg


def beam_azel_to_offset_xy(boresight_az_deg, boresight_el_deg, beam_az_deg, beam_el_deg) -> Tuple[np.ndarray, np.ndarray]:
    """Recover projected offsets from boresight and beam-center Az/El."""
    az0 = _as_float_array(boresight_az_deg)
    el0 = _as_float_array(boresight_el_deg)
    azb = _as_float_array(beam_az_deg)
    elb = _as_float_array(beam_el_deg)
    d_az_deg = ((azb - az0 + 180.0) % 360.0) - 180.0
    dx = 3600.0 * d_az_deg * np.cos(np.deg2rad(el0))
    dy = 3600.0 * (elb - el0)
    return dx, dy


def recover_el0_offset_from_xy(beam_dx_arcsec, beam_dy_arcsec, boresight_el_deg, rotation_sign) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Undo sigma * El rotation and recover the El=0 projected offset."""
    return rotate_el0_offset(beam_dx_arcsec, beam_dy_arcsec, boresight_el_deg, -float(rotation_sign))
