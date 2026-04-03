from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class ParameterDef:
    name: str
    unit: str
    default: float = 0.0
    bounds: tuple[float, float] | None = None
    description: str = ""
    external_name: str | None = None


class PointingModelBase:
    name: str = "base"

    def parameter_defs(self) -> list[ParameterDef]:
        raise NotImplementedError

    def parameter_names(self) -> list[str]:
        return [p.name for p in self.parameter_defs()]

    def axis_parameter_names(self, axis: str) -> list[str]:
        return self.parameter_names()

    def initial_guess(self, data: Any) -> dict[str, float]:
        return {p.name: p.default for p in self.parameter_defs()}

    def predict_dx_deg(self, params: dict[str, float], az_deg, el_deg):
        if hasattr(self, "predict_dx_deg"):
            return self.predict_dx_deg(params, az_deg, el_deg)  # type: ignore[attr-defined]
        raise NotImplementedError

    def predict_dy_deg(self, params: dict[str, float], az_deg, el_deg):
        if hasattr(self, "predict_dy_deg"):
            return self.predict_dy_deg(params, az_deg, el_deg)  # type: ignore[attr-defined]
        raise NotImplementedError

    # legacy aliases
    def predict_dx_deg(self, params: dict[str, float], az_deg, el_deg):
        return self.predict_dx_deg(params, az_deg, el_deg)

    def predict_dy_deg(self, params: dict[str, float], az_deg, el_deg):
        return self.predict_dy_deg(params, az_deg, el_deg)

    def to_external_params(self, internal_params: dict[str, float]) -> dict[str, float]:
        return dict(internal_params)

    def from_external_params(self, external_params: dict[str, float]) -> dict[str, float]:
        out = self.initial_guess(None)
        for k, v in external_params.items():
            if k in out:
                out[k] = float(v)
        return out

    def derived_parameters(
        self,
        params_deg: dict[str, float],
        covariance_deg: dict[tuple[str, str], float] | None = None,
        fixed_names: set[str] | None = None,
        axis: str | None = None,
    ) -> dict[str, dict[str, Any]]:
        return {}
