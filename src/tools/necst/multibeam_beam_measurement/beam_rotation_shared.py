from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Tuple
import numpy as np


@dataclass
class BeamRotationConfig:
    az_offset_arcsec: float = 0.0
    el_offset_arcsec: float = 0.0
    rotation_mode: str = "none"
    reference_angle_deg: float = 0.0
    rotation_sign: float = 1.0
    dewar_angle_deg: float = 0.0


def beam_rotation_angle_deg(boresight_el_deg: np.ndarray | float, beam: Any) -> np.ndarray:
    """Converter-compatible beam rotation angle.

    This intentionally mirrors the current formula in necst_v4_sdfits_converter.py.
    """
    el = np.asarray(boresight_el_deg, dtype=float)
    rotation_mode = str(getattr(beam, "rotation_mode", "none") or "none").lower().strip()
    reference_angle_deg = float(getattr(beam, "reference_angle_deg", 0.0))
    rotation_sign = float(getattr(beam, "rotation_sign", 1.0))
    dewar_angle_deg = float(getattr(beam, "dewar_angle_deg", 0.0))
    if rotation_mode == "none":
        return np.full_like(el, dewar_angle_deg, dtype=float)
    if rotation_mode == "elevation":
        return rotation_sign * (el - reference_angle_deg) + dewar_angle_deg
    raise ValueError(f"unsupported rotation_mode={rotation_mode!r}")


def rotate_offset_arcsec(dx0_arcsec: np.ndarray | float, dy0_arcsec: np.ndarray | float, theta_deg: np.ndarray | float) -> Tuple[np.ndarray, np.ndarray]:
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))
    dx0 = np.asarray(dx0_arcsec, dtype=float)
    dy0 = np.asarray(dy0_arcsec, dtype=float)
    dx = dx0 * np.cos(theta) - dy0 * np.sin(theta)
    dy = dx0 * np.sin(theta) + dy0 * np.cos(theta)
    return np.asarray(dx, dtype=float), np.asarray(dy, dtype=float)


def apply_beam_offset_arcsec(boresight_el_deg: np.ndarray | float, beam: Any) -> Tuple[np.ndarray, np.ndarray]:
    theta_deg = beam_rotation_angle_deg(boresight_el_deg, beam)
    return rotate_offset_arcsec(
        getattr(beam, "az_offset_arcsec", 0.0),
        getattr(beam, "el_offset_arcsec", 0.0),
        theta_deg,
    )
