from __future__ import annotations

import math
import numpy as np

from .base import ParameterDef, PointingModelBase

ARCSEC_PER_DEG = 3600.0


def _phase_component_from_c1_c2(
    params: dict[str, float],
    covariance_deg: dict[tuple[str, str], float] | None = None,
    fixed_names: set[str] | None = None,
    *,
    axis: str | None,
    theta_label: str = "az-el",
) -> dict[str, dict[str, object]]:
    if "c1" not in params or "c2" not in params:
        return {}
    c1 = float(params["c1"])
    c2 = float(params["c2"])
    amp_deg = float(math.hypot(c1, c2))
    amp_arcsec = amp_deg * ARCSEC_PER_DEG
    fixed_names = set(fixed_names or set())
    basis_fixed = {"c1": bool("c1" in fixed_names), "c2": bool("c2" in fixed_names)}
    entry: dict[str, object] = {
        "basis": f"R*sin({theta_label} - phi) = c1*sin({theta_label}) + c2*cos({theta_label})",
        "theta": theta_label,
        "basis_parameters": ["c1", "c2"],
        "axis": axis or "joint",
        "amplitude_deg": amp_deg,
        "amplitude_arcsec": amp_arcsec,
        "phase_defined": bool(amp_deg > 0.0),
        "fixed": bool(basis_fixed["c1"] and basis_fixed["c2"]),
        "basis_fixed": basis_fixed,
    }
    if amp_deg > 0.0:
        phi_deg = math.degrees(math.atan2(-c2, c1))
        entry["phase_deg"] = float(phi_deg)
        entry["phase_rad"] = float(math.radians(phi_deg))
    else:
        entry["phase_deg"] = None
        entry["phase_rad"] = None

    if covariance_deg is not None:
        v11 = float(covariance_deg.get(("c1", "c1"), 0.0))
        v22 = float(covariance_deg.get(("c2", "c2"), 0.0))
        v12 = float(covariance_deg.get(("c1", "c2"), covariance_deg.get(("c2", "c1"), 0.0)))
        cov = np.asarray([[v11, v12], [v12, v22]], dtype=float)
        if np.all(np.isfinite(cov)):
            if amp_deg > 0.0:
                grad_amp = np.asarray([c1 / amp_deg, c2 / amp_deg], dtype=float)
                var_amp = float(grad_amp @ cov @ grad_amp)
                entry["amplitude_stderr_deg"] = float(math.sqrt(max(var_amp, 0.0)))
                entry["amplitude_stderr_arcsec"] = float(math.sqrt(max(var_amp, 0.0)) * ARCSEC_PER_DEG)
                r2 = c1 * c1 + c2 * c2
                grad_phi_rad = np.asarray([c2 / r2, -c1 / r2], dtype=float)
                var_phi_rad = float(grad_phi_rad @ cov @ grad_phi_rad)
                entry["phase_stderr_deg"] = float(math.degrees(math.sqrt(max(var_phi_rad, 0.0))))
                entry["phase_stderr_rad"] = float(math.sqrt(max(var_phi_rad, 0.0)))
            else:
                entry["amplitude_stderr_deg"] = None
                entry["amplitude_stderr_arcsec"] = None
                entry["phase_stderr_deg"] = None
                entry["phase_stderr_rad"] = None
        else:
            entry["amplitude_stderr_deg"] = None
            entry["amplitude_stderr_arcsec"] = None
            entry["phase_stderr_deg"] = None
            entry["phase_stderr_rad"] = None
    return {"phase_component_az_minus_el": entry}


class Omu1p85mOpticalBasic(PointingModelBase):
    name = "omu1p85m_optical_basic"

    def axis_parameter_names(self, axis: str) -> list[str]:
        if axis == "dx":
            return ["a1", "a2", "a3", "b1", "b2"]
        if axis == "dy":
            return ["b1", "b2", "b3", "g1"]
        raise ValueError("axis must be 'dx' or 'dy'")

    def parameter_defs(self) -> list[ParameterDef]:
        return [
            ParameterDef("a1", "deg", 0.0),
            ParameterDef("a2", "deg", 0.0),
            ParameterDef("a3", "deg", 0.0),
            ParameterDef("b1", "deg", 0.0),
            ParameterDef("b2", "deg", 0.0),
            ParameterDef("b3", "deg", 0.0),
            ParameterDef("g1", "1", 0.0),
        ]

    def predict_dx_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el = np.deg2rad(np.asarray(el_deg, dtype=float))
        p = params
        return (
            p["a1"] * np.sin(el)
            + p["a2"]
            + p["a3"] * np.cos(el)
            + p["b1"] * np.sin(az) * np.sin(el)
            - p["b2"] * np.cos(az) * np.sin(el)
        )

    def predict_dy_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el = np.asarray(el_deg, dtype=float)
        p = params
        return p["b1"] * np.cos(az) + p["b2"] * np.sin(az) + p["b3"] + p["g1"] * el


class Omu1p85mRadioBasic(PointingModelBase):
    name = "omu1p85m_radio_basic"

    def axis_parameter_names(self, axis: str) -> list[str]:
        if axis == "dx":
            return ["c1", "c2", "d1"]
        if axis == "dy":
            return ["c1", "c2", "d2"]
        raise ValueError("axis must be 'dx' or 'dy'")

    def parameter_defs(self) -> list[ParameterDef]:
        return [
            ParameterDef("c1", "deg", 0.0),
            ParameterDef("c2", "deg", 0.0),
            ParameterDef("d1", "deg", 0.0),
            ParameterDef("d2", "deg", 0.0),
        ]

    def initial_guess(self, data):
        guess = super().initial_guess(data)
        if data is not None:
            if "measured_dx_arcsec" in data:
                guess["d1"] = float(np.nanmean(np.asarray(data["measured_dx_arcsec"], dtype=float))) / 3600.0
            if "measured_dy_arcsec" in data:
                guess["d2"] = float(np.nanmean(np.asarray(data["measured_dy_arcsec"], dtype=float))) / 3600.0
        return guess

    def predict_dx_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el = np.deg2rad(np.asarray(el_deg, dtype=float))
        p = params
        return p["c1"] * np.sin(az - el) + p["c2"] * np.cos(az - el) + p["d1"]

    def predict_dy_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el = np.deg2rad(np.asarray(el_deg, dtype=float))
        p = params
        return p["c1"] * np.cos(az - el) - p["c2"] * np.sin(az - el) + p["d2"]

    def derived_parameters(self, params_deg, covariance_deg=None, fixed_names=None, axis=None):
        return _phase_component_from_c1_c2(params_deg, covariance_deg, fixed_names, axis=axis, theta_label="az-el")


class Omu1p85mCombinedV1(PointingModelBase):
    name = "omu1p85m_combined_v1"

    def axis_parameter_names(self, axis: str) -> list[str]:
        if axis == "dx":
            return ["a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1"]
        if axis == "dy":
            return ["b1", "b2", "b3", "g1", "c1", "c2", "d2"]
        raise ValueError("axis must be 'dx' or 'dy'")

    def parameter_defs(self) -> list[ParameterDef]:
        return [
            ParameterDef("a1", "deg", 0.0),
            ParameterDef("a2", "deg", 0.0),
            ParameterDef("a3", "deg", 0.0),
            ParameterDef("b1", "deg", 0.0),
            ParameterDef("b2", "deg", 0.0),
            ParameterDef("b3", "deg", 0.0),
            ParameterDef("g1", "1", 0.0),
            ParameterDef("c1", "deg", 0.0),
            ParameterDef("c2", "deg", 0.0),
            ParameterDef("d1", "deg", 0.0),
            ParameterDef("d2", "deg", 0.0),
        ]

    def initial_guess(self, data):
        guess = super().initial_guess(data)
        if data is not None:
            if "measured_dx_arcsec" in data:
                guess["d1"] = float(np.nanmean(np.asarray(data["measured_dx_arcsec"], dtype=float))) / 3600.0
            if "measured_dy_arcsec" in data:
                guess["d2"] = float(np.nanmean(np.asarray(data["measured_dy_arcsec"], dtype=float))) / 3600.0
        return guess

    def predict_dx_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el = np.deg2rad(np.asarray(el_deg, dtype=float))
        p = params
        return (
            p["a1"] * np.sin(el)
            + p["a2"]
            + p["a3"] * np.cos(el)
            + p["b1"] * np.sin(az) * np.sin(el)
            - p["b2"] * np.cos(az) * np.sin(el)
            + p["c1"] * np.sin(az - el)
            + p["c2"] * np.cos(az - el)
            + p["d1"]
        )

    def predict_dy_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el_deg_arr = np.asarray(el_deg, dtype=float)
        el = np.deg2rad(el_deg_arr)
        p = params
        return (
            p["b1"] * np.cos(az)
            + p["b2"] * np.sin(az)
            + p["b3"]
            + p["g1"] * el_deg_arr
            + p["c1"] * np.cos(az - el)
            - p["c2"] * np.sin(az - el)
            + p["d2"]
        )

    def derived_parameters(self, params_deg, covariance_deg=None, fixed_names=None, axis=None):
        return _phase_component_from_c1_c2(params_deg, covariance_deg, fixed_names, axis=axis, theta_label="az-el")


class Omu1p85mCombinedReducedV1(PointingModelBase):
    """Identifiable reduced variant of the combined model.

    The original combined_v1 contains exact duplicate constant terms:
    - dx constant appears in both a2 and d1
    - dy constant appears in both b3 and d2

    This reduced form keeps a2 and b3 as the per-axis constants and removes d1/d2.
    It is numerically safer for simulation studies and future production use.
    """

    name = "omu1p85m_combined_reduced_v1"

    def axis_parameter_names(self, axis: str) -> list[str]:
        if axis == "dx":
            return ["a1", "a2", "a3", "b1", "b2", "c1", "c2"]
        if axis == "dy":
            return ["b1", "b2", "b3", "g1", "c1", "c2"]
        raise ValueError("axis must be 'dx' or 'dy'")

    def parameter_defs(self) -> list[ParameterDef]:
        return [
            ParameterDef("a1", "deg", 0.0),
            ParameterDef("a2", "deg", 0.0),
            ParameterDef("a3", "deg", 0.0),
            ParameterDef("b1", "deg", 0.0),
            ParameterDef("b2", "deg", 0.0),
            ParameterDef("b3", "deg", 0.0),
            ParameterDef("g1", "1", 0.0),
            ParameterDef("c1", "deg", 0.0),
            ParameterDef("c2", "deg", 0.0),
        ]

    def initial_guess(self, data):
        guess = super().initial_guess(data)
        if data is not None:
            if "measured_dx_arcsec" in data:
                guess["a2"] = float(np.nanmean(np.asarray(data["measured_dx_arcsec"], dtype=float))) / 3600.0
            if "measured_dy_arcsec" in data:
                guess["b3"] = float(np.nanmean(np.asarray(data["measured_dy_arcsec"], dtype=float))) / 3600.0
        return guess

    def predict_dx_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el = np.deg2rad(np.asarray(el_deg, dtype=float))
        p = params
        return (
            p["a1"] * np.sin(el)
            + p["a2"]
            + p["a3"] * np.cos(el)
            + p["b1"] * np.sin(az) * np.sin(el)
            - p["b2"] * np.cos(az) * np.sin(el)
            + p["c1"] * np.sin(az - el)
            + p["c2"] * np.cos(az - el)
        )

    def predict_dy_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el_deg_arr = np.asarray(el_deg, dtype=float)
        el = np.deg2rad(el_deg_arr)
        p = params
        return (
            p["b1"] * np.cos(az)
            + p["b2"] * np.sin(az)
            + p["b3"]
            + p["g1"] * el_deg_arr
            + p["c1"] * np.cos(az - el)
            - p["c2"] * np.sin(az - el)
        )

    def derived_parameters(self, params_deg, covariance_deg=None, fixed_names=None, axis=None):
        return _phase_component_from_c1_c2(params_deg, covariance_deg, fixed_names, axis=axis, theta_label="az-el")


class Omu1p85mActualV1(PointingModelBase):
    """OMU1P85M model matching the telescope source implementation.

    This includes e1/e2 radio-axis offset terms. For pointing-fit fits,
    keep e1=e2=0 using --fix unless they are intentionally being solved.
    """

    name = "omu1p85m_actual_v1"

    def axis_parameter_names(self, axis: str) -> list[str]:
        if axis == "dx":
            return ["a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1", "e1", "e2"]
        if axis == "dy":
            return ["b1", "b2", "b3", "g1", "c1", "c2", "d2", "e1", "e2"]
        raise ValueError("axis must be 'dx' or 'dy'")

    def parameter_defs(self) -> list[ParameterDef]:
        return [
            ParameterDef("a1", "deg", 0.0),
            ParameterDef("a2", "deg", 0.0),
            ParameterDef("a3", "deg", 0.0),
            ParameterDef("b1", "deg", 0.0),
            ParameterDef("b2", "deg", 0.0),
            ParameterDef("b3", "deg", 0.0),
            ParameterDef("g1", "1", 0.0),
            ParameterDef("c1", "deg", 0.0),
            ParameterDef("c2", "deg", 0.0),
            ParameterDef("d1", "deg", 0.0),
            ParameterDef("d2", "deg", 0.0),
            ParameterDef("e1", "deg", 0.0),
            ParameterDef("e2", "deg", 0.0),
        ]

    def initial_guess(self, data):
        guess = super().initial_guess(data)
        if data is not None:
            if "measured_dx_arcsec" in data:
                guess["d1"] = float(np.nanmean(np.asarray(data["measured_dx_arcsec"], dtype=float))) / 3600.0
            if "measured_dy_arcsec" in data:
                guess["d2"] = float(np.nanmean(np.asarray(data["measured_dy_arcsec"], dtype=float))) / 3600.0
        return guess

    def predict_dx_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el = np.deg2rad(np.asarray(el_deg, dtype=float))
        p = params
        return (
            p["a1"] * np.sin(el)
            + p["a2"]
            + p["a3"] * np.cos(el)
            + p["b1"] * np.sin(az) * np.sin(el)
            - p["b2"] * np.cos(az) * np.sin(el)
            + p["c1"] * np.sin(az - el)
            + p["c2"] * np.cos(az - el)
            + p["d1"]
            + p["e1"] * np.cos(el)
            - p["e2"] * np.sin(el)
        )

    def predict_dy_deg(self, params: dict[str, float], az_deg, el_deg):
        az = np.deg2rad(np.asarray(az_deg, dtype=float))
        el_deg_arr = np.asarray(el_deg, dtype=float)
        el = np.deg2rad(el_deg_arr)
        p = params
        return (
            p["b1"] * np.cos(az)
            + p["b2"] * np.sin(az)
            + p["b3"]
            + p["g1"] * el_deg_arr
            + p["c1"] * np.cos(az - el)
            - p["c2"] * np.sin(az - el)
            + p["d2"]
            + p["e1"] * np.sin(el)
            + p["e2"] * np.cos(el)
        )

    def derived_parameters(self, params_deg, covariance_deg=None, fixed_names=None, axis=None):
        return _phase_component_from_c1_c2(params_deg, covariance_deg, fixed_names, axis=axis, theta_label="az-el")


# legacy aliases
Omu1p85mOpticalBasic.predict_dx_deg = Omu1p85mOpticalBasic.predict_dx_deg
Omu1p85mOpticalBasic.predict_dy_deg = Omu1p85mOpticalBasic.predict_dy_deg
Omu1p85mRadioBasic.predict_dx_deg = Omu1p85mRadioBasic.predict_dx_deg
Omu1p85mRadioBasic.predict_dy_deg = Omu1p85mRadioBasic.predict_dy_deg
Omu1p85mCombinedV1.predict_dx_deg = Omu1p85mCombinedV1.predict_dx_deg
Omu1p85mCombinedV1.predict_dy_deg = Omu1p85mCombinedV1.predict_dy_deg
Omu1p85mCombinedReducedV1.predict_dx_deg = Omu1p85mCombinedReducedV1.predict_dx_deg
Omu1p85mCombinedReducedV1.predict_dy_deg = Omu1p85mCombinedReducedV1.predict_dy_deg

Omu1p85mActualV1.predict_dx_deg = Omu1p85mActualV1.predict_dx_deg
Omu1p85mActualV1.predict_dy_deg = Omu1p85mActualV1.predict_dy_deg
