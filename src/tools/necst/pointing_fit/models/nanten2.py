from __future__ import annotations

import numpy as np

from .base import ParameterDef, PointingModelBase

ARCSEC_PER_DEG = 3600.0


class Nanten2ActualV1(PointingModelBase):
    """NANTEN2 model matching the telescope source implementation.

    Parameters g/gg/ggg/gggg and their radio counterparts are kept in the same
    source-side units as the telescope code: arcsec/deg^n. The predictor
    converts them to degree-valued offsets internally by dividing the evaluated
    polynomial by 3600.
    """

    name = "nanten2_actual_v1"

    def axis_parameter_names(self, axis: str) -> list[str]:
        if axis == "dx":
            return [
                "chi_Az", "omega_Az", "eps", "chi2_Az", "omega2_Az",
                "dAz", "de", "cor_v", "cor_p", "de_radio",
            ]
        if axis == "dy":
            return [
                "chi_El", "omega_El", "chi2_El", "omega2_El",
                "g", "gg", "ggg", "gggg",
                "dEl", "g_radio", "gg_radio", "ggg_radio", "gggg_radio",
                "cor_v", "cor_p", "del_radio",
            ]
        raise ValueError("axis must be 'dx' or 'dy'")

    def parameter_defs(self) -> list[ParameterDef]:
        return [
            ParameterDef("dAz", "deg", 0.0),
            ParameterDef("de", "deg", 0.0),
            ParameterDef("chi_Az", "deg", 0.0),
            ParameterDef("omega_Az", "deg", 0.0),
            ParameterDef("eps", "deg", 0.0),
            ParameterDef("chi2_Az", "deg", 0.0),
            ParameterDef("omega2_Az", "deg", 0.0),
            ParameterDef("chi_El", "deg", 0.0),
            ParameterDef("omega_El", "deg", 0.0),
            ParameterDef("chi2_El", "deg", 0.0),
            ParameterDef("omega2_El", "deg", 0.0),
            ParameterDef("g", "arcsec/deg", 0.0),
            ParameterDef("gg", "arcsec/deg^2", 0.0),
            ParameterDef("ggg", "arcsec/deg^3", 0.0),
            ParameterDef("gggg", "arcsec/deg^4", 0.0),
            ParameterDef("dEl", "deg", 0.0),
            ParameterDef("de_radio", "deg", 0.0),
            ParameterDef("del_radio", "deg", 0.0),
            ParameterDef("cor_v", "deg", 0.0),
            ParameterDef("cor_p", "deg", 0.0),
            ParameterDef("g_radio", "arcsec/deg", 0.0),
            ParameterDef("gg_radio", "arcsec/deg^2", 0.0),
            ParameterDef("ggg_radio", "arcsec/deg^3", 0.0),
            ParameterDef("gggg_radio", "arcsec/deg^4", 0.0),
        ]

    def initial_guess(self, data):
        guess = super().initial_guess(data)
        if data is not None:
            if "measured_dx_arcsec" in data:
                guess["de"] = float(np.nanmean(np.asarray(data["measured_dx_arcsec"], dtype=float))) / ARCSEC_PER_DEG
            if "measured_dy_arcsec" in data:
                guess["dEl"] = float(np.nanmean(np.asarray(data["measured_dy_arcsec"], dtype=float))) / ARCSEC_PER_DEG
        return guess

    @staticmethod
    def _grav_deg(el_deg_arr: np.ndarray, p: dict[str, float], radio: bool = False) -> np.ndarray:
        if radio:
            g = float(p["g_radio"])
            gg = float(p["gg_radio"])
            ggg = float(p["ggg_radio"])
            gggg = float(p["gggg_radio"])
        else:
            g = float(p["g"])
            gg = float(p["gg"])
            ggg = float(p["ggg"])
            gggg = float(p["gggg"])
        arcsec = g * el_deg_arr + gg * el_deg_arr**2 + ggg * el_deg_arr**3 + gggg * el_deg_arr**4
        return arcsec / ARCSEC_PER_DEG



    def _analytic_jacobian_dx_deg(self, params: dict[str, float], az_deg, el_deg, names: list[str] | None = None) -> np.ndarray:
        az_deg_arr = np.asarray(az_deg, dtype=float)
        el_deg_arr = np.asarray(el_deg, dtype=float)
        el = np.deg2rad(el_deg_arr)
        sin_el = np.sin(el)
        cos_el = np.cos(el)
        omega_az = np.deg2rad(float(params["omega_Az"]) - az_deg_arr)
        omega2_az = 2.0 * np.deg2rad(float(params["omega2_Az"]) - az_deg_arr)
        psi = el + np.deg2rad(float(params["cor_p"]))
        names = list(names or self.axis_parameter_names("dx"))
        cols: list[np.ndarray] = []
        for name in names:
            if name == "chi_Az":
                col = np.sin(omega_az) * sin_el
            elif name == "omega_Az":
                col = float(params["chi_Az"]) * np.cos(omega_az) * sin_el * np.deg2rad(1.0)
            elif name == "eps":
                col = sin_el
            elif name == "chi2_Az":
                col = np.sin(omega2_az) * sin_el
            elif name == "omega2_Az":
                col = float(params["chi2_Az"]) * np.cos(omega2_az) * sin_el * (2.0 * np.deg2rad(1.0))
            elif name == "dAz":
                col = cos_el
            elif name == "de":
                col = np.ones_like(az_deg_arr, dtype=float)
            elif name == "cor_v":
                col = np.cos(psi)
            elif name == "cor_p":
                col = -float(params["cor_v"]) * np.sin(psi) * np.deg2rad(1.0)
            elif name == "de_radio":
                col = np.ones_like(az_deg_arr, dtype=float)
            else:
                col = np.zeros_like(az_deg_arr, dtype=float)
            cols.append(np.asarray(col, dtype=float))
        return np.column_stack(cols) if cols else np.zeros((az_deg_arr.size, 0), dtype=float)

    def _analytic_jacobian_dy_deg(self, params: dict[str, float], az_deg, el_deg, names: list[str] | None = None) -> np.ndarray:
        az_deg_arr = np.asarray(az_deg, dtype=float)
        el_deg_arr = np.asarray(el_deg, dtype=float)
        el = np.deg2rad(el_deg_arr)
        omega_el = np.deg2rad(float(params["omega_El"]) - az_deg_arr)
        omega2_el = 2.0 * np.deg2rad(float(params["omega2_El"]) - az_deg_arr)
        psi = el + np.deg2rad(float(params["cor_p"]))
        names = list(names or self.axis_parameter_names("dy"))
        cols: list[np.ndarray] = []
        for name in names:
            if name == "chi_El":
                col = -np.cos(omega_el)
            elif name == "omega_El":
                col = float(params["chi_El"]) * np.sin(omega_el) * np.deg2rad(1.0)
            elif name == "chi2_El":
                col = -np.cos(omega2_el)
            elif name == "omega2_El":
                col = float(params["chi2_El"]) * np.sin(omega2_el) * (2.0 * np.deg2rad(1.0))
            elif name == "g":
                col = el_deg_arr / ARCSEC_PER_DEG
            elif name == "gg":
                col = el_deg_arr**2 / ARCSEC_PER_DEG
            elif name == "ggg":
                col = el_deg_arr**3 / ARCSEC_PER_DEG
            elif name == "gggg":
                col = el_deg_arr**4 / ARCSEC_PER_DEG
            elif name == "dEl":
                col = np.ones_like(az_deg_arr, dtype=float)
            elif name == "g_radio":
                col = el_deg_arr / ARCSEC_PER_DEG
            elif name == "gg_radio":
                col = el_deg_arr**2 / ARCSEC_PER_DEG
            elif name == "ggg_radio":
                col = el_deg_arr**3 / ARCSEC_PER_DEG
            elif name == "gggg_radio":
                col = el_deg_arr**4 / ARCSEC_PER_DEG
            elif name == "cor_v":
                col = -np.sin(psi)
            elif name == "cor_p":
                col = -float(params["cor_v"]) * np.cos(psi) * np.deg2rad(1.0)
            elif name == "del_radio":
                col = np.ones_like(az_deg_arr, dtype=float)
            else:
                col = np.zeros_like(az_deg_arr, dtype=float)
            cols.append(np.asarray(col, dtype=float))
        return np.column_stack(cols) if cols else np.zeros((az_deg_arr.size, 0), dtype=float)

    def supports_linearized_delta(self) -> bool:
        return True

    def predict_dx_deg(self, params: dict[str, float], az_deg, el_deg):
        az_deg_arr = np.asarray(az_deg, dtype=float)
        el_deg_arr = np.asarray(el_deg, dtype=float)
        az = np.deg2rad(az_deg_arr)
        el = np.deg2rad(el_deg_arr)
        p = params
        return (
            p["chi_Az"] * np.sin(np.deg2rad(float(p["omega_Az"]) - az_deg_arr)) * np.sin(el)
            + p["eps"] * np.sin(el)
            + p["chi2_Az"] * np.sin(2.0 * np.deg2rad(float(p["omega2_Az"]) - az_deg_arr)) * np.sin(el)
            + p["dAz"] * np.cos(el)
            + p["de"]
            + p["cor_v"] * np.cos(el + np.deg2rad(float(p["cor_p"])))
            + p["de_radio"]
        )

    def predict_dy_deg(self, params: dict[str, float], az_deg, el_deg):
        az_deg_arr = np.asarray(az_deg, dtype=float)
        el_deg_arr = np.asarray(el_deg, dtype=float)
        az = np.deg2rad(az_deg_arr)
        el = np.deg2rad(el_deg_arr)
        p = params
        return (
            -1.0 * p["chi_El"] * np.cos(np.deg2rad(float(p["omega_El"]) - az_deg_arr))
            - p["chi2_El"] * np.cos(2.0 * np.deg2rad(float(p["omega2_El"]) - az_deg_arr))
            + self._grav_deg(el_deg_arr, p)
            + p["dEl"]
            + self._grav_deg(el_deg_arr, p, radio=True)
            - p["cor_v"] * np.sin(el + np.deg2rad(float(p["cor_p"])))
            + p["del_radio"]
        )
