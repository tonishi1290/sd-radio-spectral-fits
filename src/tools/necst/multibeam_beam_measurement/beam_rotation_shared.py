from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Tuple
import numpy as np

from .pure_rotation_model import rotate_el0_offset


@dataclass
class BeamRotationConfig:
    beam_model: str = "legacy"
    az_offset_arcsec: float = 0.0
    el_offset_arcsec: float = 0.0
    rotation_mode: str = "none"
    reference_angle_deg: float = 0.0
    rotation_sign: float = 1.0
    rotation_slope_deg_per_deg: float | None = None
    dewar_angle_deg: float = 0.0
    offset_x_el0_arcsec: float | None = None
    offset_y_el0_arcsec: float | None = None
    pure_rotation_sign: float | None = None


def beam_rotation_angle_deg(boresight_el_deg: np.ndarray | float, beam: Any) -> np.ndarray:
    """Beam rotation angle in degrees.

    For ``pure_rotation_v1`` this is simply ``rotation_sign * El``.
    Legacy mode mirrors necst_v4_sdfits_converter.py.
    """
    el = np.asarray(boresight_el_deg, dtype=float)
    beam_model = str(getattr(beam, "beam_model", getattr(beam, "model", "legacy")) or "legacy").strip().lower()
    if beam_model == "pure_rotation_v1":
        pure_sign = getattr(beam, "pure_rotation_sign", None)
        if pure_sign is None:
            pure_sign = getattr(beam, "rotation_sign", 1.0)
        return np.asarray(float(pure_sign) * el, dtype=float)
    rotation_mode = str(getattr(beam, "rotation_mode", "none") or "none").lower().strip()
    reference_angle_deg = float(getattr(beam, "reference_angle_deg", 0.0))
    rotation_sign = float(getattr(beam, "rotation_sign", 1.0))
    rotation_slope_deg_per_deg = getattr(beam, "rotation_slope_deg_per_deg", None)
    slope = float(rotation_slope_deg_per_deg) if rotation_slope_deg_per_deg is not None else rotation_sign
    dewar_angle_deg = float(getattr(beam, "dewar_angle_deg", 0.0))
    if rotation_mode == "none":
        return np.full_like(el, dewar_angle_deg, dtype=float)
    if rotation_mode == "elevation":
        return slope * (el - reference_angle_deg) + dewar_angle_deg
    raise ValueError(f"unsupported rotation_mode={rotation_mode!r}")


def rotate_offset_arcsec(dx0_arcsec: np.ndarray | float, dy0_arcsec: np.ndarray | float, theta_deg: np.ndarray | float) -> Tuple[np.ndarray, np.ndarray]:
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))
    dx0 = np.asarray(dx0_arcsec, dtype=float)
    dy0 = np.asarray(dy0_arcsec, dtype=float)
    dx = dx0 * np.cos(theta) - dy0 * np.sin(theta)
    dy = dx0 * np.sin(theta) + dy0 * np.cos(theta)
    return np.asarray(dx, dtype=float), np.asarray(dy, dtype=float)


def apply_beam_offset_arcsec(boresight_el_deg: np.ndarray | float, beam: Any) -> Tuple[np.ndarray, np.ndarray]:
    beam_model = str(getattr(beam, "beam_model", getattr(beam, "model", "legacy")) or "legacy").strip().lower()
    if beam_model == "pure_rotation_v1":
        dx, dy, _ = rotate_el0_offset(
            getattr(beam, "offset_x_el0_arcsec", getattr(beam, "pure_rotation_offset_x_el0_arcsec", 0.0)),
            getattr(beam, "offset_y_el0_arcsec", getattr(beam, "pure_rotation_offset_y_el0_arcsec", 0.0)),
            boresight_el_deg,
            getattr(beam, "pure_rotation_sign", getattr(beam, "rotation_sign", 1.0)),
        )
        return dx, dy
    theta_deg = beam_rotation_angle_deg(boresight_el_deg, beam)
    return rotate_offset_arcsec(
        getattr(beam, "az_offset_arcsec", 0.0),
        getattr(beam, "el_offset_arcsec", 0.0),
        theta_deg,
    )
