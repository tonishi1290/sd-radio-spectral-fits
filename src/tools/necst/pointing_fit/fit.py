from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .__init__ import __version__
from .io_utils import (
    ARCSEC_PER_DEG,
    dump_json,
    parse_param_file_sections,
    parse_pointing_param_toml,
    read_csv,
    write_csv,
    write_pointing_param_toml,
    write_pointing_param_toml_from_template,
    parse_param_key_styles,
)
from .models import canonical_model_name, resolve_model


@dataclass
class FitOptions:
    input_csv: str
    outdir: str
    model: str
    fit_target: str = "absolute"
    solve_mode: str = "joint"
    robust_loss: str = "soft_l1"
    f_scale: str = "auto"
    clip_sigma: float = 0.0
    clip_max_iter: int = 0
    bootstrap: int = 0
    seed: int = 0
    fixed_param_specs: list[str] | None = None
    fixed_param_file: str | None = None


@dataclass
class FitCore:
    model: Any
    all_names: list[str]
    free_names: list[str]
    fixed_params: dict[str, float]
    residual_fun: Any
    x_to_params: Any
    target_deg: np.ndarray
    init_params: dict[str, float]
    axis: str | None = None
    predict_fun: Any | None = None
    base_params: dict[str, float] | None = None
    linearized_delta: bool = False


@dataclass
class SimpleFitResult:
    x: np.ndarray
    fun: np.ndarray
    jac: np.ndarray
    cost: float
    success: bool
    message: str
    nfev: int = 1
    njev: int = 1


EXPECTED_COLS = [
    "dataset_id",
    "dt_cross_id",
    "cross_id",
    "az_deg",
    "el_deg",
]


_DELTA_ERR_DX_CANDIDATES = ["delta_err_dx_arcsec", "measured_dx_arcsec"]
_DELTA_ERR_DY_CANDIDATES = ["delta_err_dy_arcsec", "measured_dy_arcsec"]
_USED_MODEL_DX_CANDIDATES = ["used_model_dx_arcsec", "applied_model_dx_arcsec"]
_USED_MODEL_DY_CANDIDATES = ["used_model_dy_arcsec", "applied_model_dy_arcsec"]


def _first_present_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for name in candidates:
        if name in df.columns:
            return name
    return None


def _require_first_present_column(df: pd.DataFrame, candidates: list[str], purpose: str) -> str:
    name = _first_present_column(df, candidates)
    if name is None:
        raise KeyError(f"Missing required column for {purpose}: tried {candidates}")
    return name


def _first_existing_used_param_file(df: pd.DataFrame) -> str | None:
    if "used_param_file" not in df.columns:
        return None
    for value in pd.Series(df["used_param_file"]):
        if pd.isna(value):
            continue
        path = str(value).strip()
        if not path:
            continue
        if Path(path).exists():
            return path
    return None



def infer_output_model_label(df: pd.DataFrame, requested_model: str) -> str:
    requested = str(requested_model or "").strip()
    requested_canonical = canonical_model_name(requested)
    used_param_file = _first_existing_used_param_file(df)
    if used_param_file:
        try:
            raw_model, _ = parse_pointing_param_toml(used_param_file)
        except Exception:
            raw_model = None
        if raw_model and canonical_model_name(str(raw_model)) == requested_canonical:
            return str(raw_model)
    return requested or requested_canonical



def _template_for_full_param_output(df: pd.DataFrame, fixed_param_file: str | None) -> str | None:
    used_param_file = _first_existing_used_param_file(df)
    if used_param_file:
        return used_param_file
    if fixed_param_file and Path(fixed_param_file).exists():
        _, _, sec = parse_param_file_sections(fixed_param_file, ("pointing_params", "fixed_params"))
        if sec in {"pointing_params", "top_level"}:
            return fixed_param_file
    return None


def _template_for_fixed_param_output(fixed_param_file: str | None) -> str | None:
    if fixed_param_file and Path(fixed_param_file).exists():
        return fixed_param_file
    return None


def _write_param_output_with_optional_template(
    output_path: str | Path,
    *,
    template_path: str | None,
    params_deg: dict[str, float],
    model_name: str,
    section_name: str = "pointing_params",
    key_style_map: dict[str, dict[str, Any]] | None = None,
) -> None:
    if template_path and Path(template_path).exists():
        write_pointing_param_toml_from_template(
            template_path,
            output_path,
            params_deg=params_deg,
            model_name=model_name,
            section_names=(section_name, "pointing_params", "fixed_params"),
        )
        return
    write_pointing_param_toml(
        output_path,
        model_name=model_name,
        params_deg=params_deg,
        unit="preserve",
        metadata_extra={"model": model_name},
        section_name=section_name,
        key_style_map=key_style_map,
    )

def radial_rms_arcsec(r_daz_deg: np.ndarray, r_del_deg: np.ndarray) -> dict[str, float]:
    r_daz = r_daz_deg * ARCSEC_PER_DEG
    r_del = r_del_deg * ARCSEC_PER_DEG
    rr = np.sqrt(r_daz**2 + r_del**2)
    return {
        "dx": float(np.sqrt(np.nanmean(r_daz**2))),
        "dy": float(np.sqrt(np.nanmean(r_del**2))),
        "radial": float(np.sqrt(np.nanmean(rr**2))),
    }


def mad_scale(x: np.ndarray) -> float:
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    sigma = 1.4826 * mad
    return float(sigma if sigma > 0 and np.isfinite(sigma) else max(float(np.nanstd(x)), 1e-9))


def prepare_targets(df: pd.DataFrame, fit_target: str) -> pd.DataFrame:
    out = df.copy()
    delta_dx_col = _require_first_present_column(out, _DELTA_ERR_DX_CANDIDATES, "delta_err dx")
    delta_dy_col = _require_first_present_column(out, _DELTA_ERR_DY_CANDIDATES, "delta_err dy")
    used_dx_col = _first_present_column(out, _USED_MODEL_DX_CANDIDATES)
    used_dy_col = _first_present_column(out, _USED_MODEL_DY_CANDIDATES)

    out["measured_dx_arcsec"] = pd.to_numeric(out[delta_dx_col], errors="coerce").astype(float)
    out["measured_dy_arcsec"] = pd.to_numeric(out[delta_dy_col], errors="coerce").astype(float)
    out["delta_err_dx_arcsec"] = out["measured_dx_arcsec"]
    out["delta_err_dy_arcsec"] = out["measured_dy_arcsec"]
    if used_dx_col is not None:
        out["used_model_dx_arcsec"] = pd.to_numeric(out[used_dx_col], errors="coerce").astype(float)
        out["applied_model_dx_arcsec"] = out["used_model_dx_arcsec"]
    if used_dy_col is not None:
        out["used_model_dy_arcsec"] = pd.to_numeric(out[used_dy_col], errors="coerce").astype(float)
        out["applied_model_dy_arcsec"] = out["used_model_dy_arcsec"]

    if "sign_convention" in out.columns:
        warnings.warn(
            "Input CSV contains legacy column 'sign_convention'. It is ignored and will be dropped; absolute fit always uses used_model + delta_err.",
            DeprecationWarning,
            stacklevel=2,
        )
        out = out.drop(columns=["sign_convention"])
    out["absolute_target_rule"] = "used_model_plus_delta_err"

    if fit_target == "delta":
        out["target_dx_arcsec"] = out["delta_err_dx_arcsec"].astype(float)
        out["target_dy_arcsec"] = out["delta_err_dy_arcsec"].astype(float)
        return out
    if fit_target != "absolute":
        raise ValueError("fit_target must be 'absolute' or 'delta'.")
    if used_dx_col is None or used_dy_col is None:
        raise KeyError("absolute fit requires used_model/applied_model columns in the normalized CSV.")
    out["target_dx_arcsec"] = out["used_model_dx_arcsec"].astype(float) + out["delta_err_dx_arcsec"].astype(float)
    out["target_dy_arcsec"] = out["used_model_dy_arcsec"].astype(float) + out["delta_err_dy_arcsec"].astype(float)
    return out


def _parse_assignment(spec: str) -> tuple[str, float]:
    if "=" not in spec:
        raise ValueError(f"Invalid fixed parameter specification: {spec}")
    lhs, rhs = spec.split("=", 1)
    lhs = lhs.strip()
    rhs = rhs.strip()
    if not lhs:
        raise ValueError(f"Invalid fixed parameter specification: {spec}")
    unit = None
    name = lhs
    if "[" in lhs and lhs.endswith("]"):
        name, unit = lhs[:-1].split("[", 1)
        name = name.strip()
        unit = unit.strip()
    value = float(rhs)
    if unit == "arcsec":
        value = value / ARCSEC_PER_DEG
    return name, value


def resolve_fixed_params(model_name: str, fixed_param_specs: list[str] | None, fixed_param_file: str | None) -> tuple[dict[str, float], dict[str, str]]:
    model = resolve_model(model_name)
    valid = set(model.parameter_names())
    fixed: dict[str, float] = {}
    sources: dict[str, str] = {}
    if fixed_param_file:
        _, params, sec = parse_param_file_sections(fixed_param_file, ("fixed_params", "pointing_params"))
        for key, val in params.items():
            if key not in valid:
                raise KeyError(f"Unknown fixed parameter in {fixed_param_file}: {key}")
            fixed[key] = float(val)
            sources[key] = f"file:{Path(fixed_param_file).name}:{sec or 'unknown'}"
    for spec in fixed_param_specs or []:
        key, val = _parse_assignment(spec)
        if key not in valid:
            raise KeyError(f"Unknown fixed parameter in --fix: {key}")
        fixed[key] = float(val)
        sources[key] = "cli"
    return fixed, sources


def _default_params_for_model(model_name: str) -> dict[str, float]:
    model = resolve_model(model_name)
    return {name: 0.0 for name in model.parameter_names()}


def _param_file_model_is_compatible(file_model: str | None, file_params: dict[str, float], model_name: str) -> bool:
    if not file_model:
        return True
    requested = canonical_model_name(model_name).split("__", 1)[0]
    parsed = canonical_model_name(str(file_model)).split("__", 1)[0]
    if parsed == requested:
        return True
    try:
        valid = set(resolve_model(model_name).parameter_names())
    except Exception:
        return False
    return set(file_params).issubset(valid)


def resolve_delta_base_params(df: pd.DataFrame, model_name: str, fixed_param_file: str | None) -> dict[str, float] | None:
    model_name = canonical_model_name(model_name)
    model = resolve_model(model_name)
    base = _default_params_for_model(model_name)

    if "used_param_hash" in df.columns:
        hashes = sorted({str(x) for x in df["used_param_hash"].astype(str) if str(x) and str(x) != "nan"})
        if len(hashes) > 1:
            raise ValueError("delta fit with Jacobian linearization requires a single used_param_hash in the normalized CSV.")
    if "used_param_file" in df.columns:
        paths = sorted({str(x) for x in df["used_param_file"].astype(str) if str(x) and str(x) != "nan"})
        if len(paths) == 1:
            used_path = Path(paths[0])
            if used_path.exists():
                used_model, params, _ = parse_param_file_sections(used_path, ("pointing_params", "fixed_params"))
                if not _param_file_model_is_compatible(used_model, params, model_name):
                    raise ValueError(f"used_param model mismatch for delta fit: {used_model} vs {model_name}")
                for key, val in params.items():
                    if key in base:
                        base[key] = float(val)
                return base
        if len(paths) > 1:
            raise ValueError("delta fit with Jacobian linearization requires a single used_param_file in the normalized CSV.")

    if fixed_param_file:
        file_model, params, sec = parse_param_file_sections(fixed_param_file, ("pointing_params", "fixed_params"))
        if sec == "pointing_params":
            if not _param_file_model_is_compatible(file_model, params, model_name):
                raise ValueError(f"fixed_param_file model mismatch for delta fit: {file_model} vs {model_name}")
            for key, val in params.items():
                if key in base:
                    base[key] = float(val)
            return base

    if getattr(model, "supports_linearized_delta", lambda: False)():
        return base
    return None


def _free_index_lookup(free_names: list[str]) -> dict[str, int]:
    return {name: idx for idx, name in enumerate(free_names)}


def build_axis_core(
    df: pd.DataFrame,
    model_name: str,
    axis: str,
    fixed_params: dict[str, float],
    fit_target: str,
    delta_base_params: dict[str, float] | None = None,
) -> FitCore:
    if axis not in {"dx", "dy"}:
        raise ValueError("axis must be 'dx' or 'dy'")
    model = resolve_model(model_name)
    all_names = model.parameter_names()
    init_params = dict(model.initial_guess(df))
    active_names = list(model.axis_parameter_names(axis)) if hasattr(model, "axis_parameter_names") else list(all_names)
    free_names = [n for n in active_names if n not in fixed_params]
    az = df["az_deg"].to_numpy(dtype=float)
    el = df["el_deg"].to_numpy(dtype=float)
    target_deg = df[f"target_{axis}_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG

    use_linearized_delta = (
        fit_target == "delta"
        and delta_base_params is not None
        and getattr(model, "supports_linearized_delta", lambda: False)()
        and hasattr(model, f"_analytic_jacobian_{axis}_deg")
    )

    if use_linearized_delta:
        base_params = dict(_default_params_for_model(model_name))
        base_params.update(delta_base_params)
        for k, v in fixed_params.items():
            base_params[k] = float(v)
        jac_fun = getattr(model, f"_analytic_jacobian_{axis}_deg")
        jac_full = np.asarray(jac_fun(base_params, az, el, active_names), dtype=float)
        free_lookup = _free_index_lookup(free_names)
        col_indices = [active_names.index(name) for name in free_names]
        jac_free = jac_full[:, col_indices] if col_indices else np.zeros((az.size, 0), dtype=float)

        def x_to_params(x_free: np.ndarray) -> dict[str, float]:
            p = {k: 0.0 for k in all_names}
            for k, v in zip(free_names, x_free):
                p[k] = float(v)
            return p

        def predict_axis(x_free: np.ndarray) -> np.ndarray:
            x_arr = np.asarray(x_free, dtype=float)
            if jac_free.shape[1] == 0:
                return np.zeros_like(target_deg, dtype=float)
            return np.asarray(jac_free @ x_arr, dtype=float)

        def residuals_axis(x_free: np.ndarray) -> np.ndarray:
            return target_deg - predict_axis(x_free)

        return FitCore(
            model=model,
            all_names=all_names,
            free_names=free_names,
            fixed_params={k: v for k, v in fixed_params.items() if k in active_names},
            residual_fun=residuals_axis,
            x_to_params=x_to_params,
            target_deg=target_deg,
            init_params={k: 0.0 for k in all_names},
            axis=axis,
            predict_fun=predict_axis,
            base_params=base_params,
            linearized_delta=True,
        )

    def x_to_params(x_free: np.ndarray) -> dict[str, float]:
        p = dict(init_params)
        for key in all_names:
            p.setdefault(key, 0.0)
        for k, v in fixed_params.items():
            p[k] = float(v)
        for k, v in zip(free_names, x_free):
            p[k] = float(v)
        return p

    def predict_axis(x_free: np.ndarray) -> np.ndarray:
        p = x_to_params(x_free)
        if axis == "dx":
            return np.asarray(model.predict_dx_deg(p, az, el), dtype=float)
        return np.asarray(model.predict_dy_deg(p, az, el), dtype=float)

    def residuals_axis(x_free: np.ndarray) -> np.ndarray:
        return target_deg - predict_axis(x_free)

    return FitCore(
        model=model,
        all_names=all_names,
        free_names=free_names,
        fixed_params={k: v for k, v in fixed_params.items() if k in active_names},
        residual_fun=residuals_axis,
        x_to_params=x_to_params,
        target_deg=target_deg,
        init_params=init_params,
        axis=axis,
        predict_fun=predict_axis,
        base_params=None,
        linearized_delta=False,
    )


def build_joint_core(
    df: pd.DataFrame,
    model_name: str,
    fixed_params: dict[str, float],
    fit_target: str,
    delta_base_params: dict[str, float] | None = None,
) -> FitCore:
    model = resolve_model(model_name)
    all_names = model.parameter_names()
    init_params = dict(model.initial_guess(df))
    free_names = [n for n in all_names if n not in fixed_params]
    az = df["az_deg"].to_numpy(dtype=float)
    el = df["el_deg"].to_numpy(dtype=float)
    y_daz = df["target_dx_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG
    y_del = df["target_dy_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG

    use_linearized_delta = (
        fit_target == "delta"
        and delta_base_params is not None
        and getattr(model, "supports_linearized_delta", lambda: False)()
        and hasattr(model, "_analytic_jacobian_dx_deg")
        and hasattr(model, "_analytic_jacobian_dy_deg")
    )

    if use_linearized_delta:
        base_params = dict(_default_params_for_model(model_name))
        base_params.update(delta_base_params)
        for k, v in fixed_params.items():
            base_params[k] = float(v)
        jdx_full = np.asarray(model._analytic_jacobian_dx_deg(base_params, az, el, all_names), dtype=float)
        jdy_full = np.asarray(model._analytic_jacobian_dy_deg(base_params, az, el, all_names), dtype=float)
        col_indices = [all_names.index(name) for name in free_names]
        jdx = jdx_full[:, col_indices] if col_indices else np.zeros((az.size, 0), dtype=float)
        jdy = jdy_full[:, col_indices] if col_indices else np.zeros((az.size, 0), dtype=float)

        def x_to_params(x_free: np.ndarray) -> dict[str, float]:
            p = {k: 0.0 for k in all_names}
            for k, v in zip(free_names, x_free):
                p[k] = float(v)
            return p

        def predict_joint(x_free: np.ndarray) -> np.ndarray:
            x_arr = np.asarray(x_free, dtype=float)
            if len(free_names) == 0:
                pred_dx = np.zeros_like(y_daz, dtype=float)
                pred_dy = np.zeros_like(y_del, dtype=float)
            else:
                pred_dx = np.asarray(jdx @ x_arr, dtype=float)
                pred_dy = np.asarray(jdy @ x_arr, dtype=float)
            return np.concatenate([pred_dx, pred_dy])

        def residuals_joint(x_free: np.ndarray) -> np.ndarray:
            pred = predict_joint(x_free)
            return np.concatenate([y_daz, y_del]) - pred

        return FitCore(
            model=model,
            all_names=all_names,
            free_names=free_names,
            fixed_params=dict(fixed_params),
            residual_fun=residuals_joint,
            x_to_params=x_to_params,
            target_deg=np.concatenate([y_daz, y_del]),
            init_params={k: 0.0 for k in all_names},
            axis=None,
            predict_fun=predict_joint,
            base_params=base_params,
            linearized_delta=True,
        )

    def x_to_params(x_free: np.ndarray) -> dict[str, float]:
        p = dict(init_params)
        for k in all_names:
            p.setdefault(k, 0.0)
        for k, v in fixed_params.items():
            p[k] = float(v)
        for k, v in zip(free_names, x_free):
            p[k] = float(v)
        return p

    def predict_joint(x_free: np.ndarray) -> np.ndarray:
        p = x_to_params(x_free)
        pred_daz = np.asarray(model.predict_dx_deg(p, az, el), dtype=float)
        pred_del = np.asarray(model.predict_dy_deg(p, az, el), dtype=float)
        return np.concatenate([pred_daz, pred_del])

    def residuals_joint(x_free: np.ndarray) -> np.ndarray:
        return np.concatenate([y_daz, y_del]) - predict_joint(x_free)

    return FitCore(
        model=model,
        all_names=all_names,
        free_names=free_names,
        fixed_params=dict(fixed_params),
        residual_fun=residuals_joint,
        x_to_params=x_to_params,
        target_deg=np.concatenate([y_daz, y_del]),
        init_params=init_params,
        axis=None,
        predict_fun=predict_joint,
        base_params=None,
        linearized_delta=False,
    )


def fit_with_core(core: FitCore, robust_loss: str, f_scale: str | float):
    x0_dict = dict(core.init_params)
    x0 = np.array([x0_dict.get(n, 0.0) for n in core.free_names], dtype=float)
    resid0 = core.residual_fun(x0)
    scale = mad_scale(resid0) if f_scale == "auto" else float(f_scale)
    if len(core.free_names) == 0:
        fun = core.residual_fun(x0)
        jac = np.zeros((fun.size, 0), dtype=float)
        result = SimpleFitResult(
            x=x0,
            fun=fun,
            jac=jac,
            cost=0.5 * float(np.sum(fun**2)),
            success=True,
            message="No free parameters: evaluated fixed model only.",
        )
        return result, scale
    result = least_squares(core.residual_fun, x0=x0, loss=robust_loss, f_scale=scale)
    return result, scale


def params_to_dict(names: list[str], x: np.ndarray) -> dict[str, float]:
    return {k: float(v) for k, v in zip(names, x)}


def fit_once(
    df: pd.DataFrame,
    model_name: str,
    solve_mode: str,
    robust_loss: str,
    f_scale: str | float,
    fixed_params: dict[str, float],
    fit_target: str,
    delta_base_params: dict[str, float] | None = None,
):
    if solve_mode == "joint":
        core = build_joint_core(df, model_name, fixed_params, fit_target, delta_base_params)
        result, scale = fit_with_core(core, robust_loss, f_scale)
        return core, result, scale
    if solve_mode == "separate":
        cores: dict[str, FitCore] = {}
        results: dict[str, Any] = {}
        scales: dict[str, float] = {}
        for axis in ("dx", "dy"):
            core_axis = build_axis_core(df, model_name, axis, fixed_params, fit_target, delta_base_params)
            result_axis, scale_axis = fit_with_core(core_axis, robust_loss, f_scale)
            cores[axis] = core_axis
            results[axis] = result_axis
            scales[axis] = scale_axis
        return cores, results, scales
    raise ValueError("solve_mode must be 'joint' or 'separate'.")


def compute_covariance(result, n_obs: int, n_params: int) -> tuple[np.ndarray | None, np.ndarray | None]:
    if n_params == 0:
        return np.zeros((0, 0), dtype=float), np.zeros((0,), dtype=float)
    if result.jac is None:
        return None, None
    jac = np.asarray(result.jac, dtype=float)
    dof = max(n_obs - n_params, 1)
    try:
        jtj_inv = np.linalg.pinv(jac.T @ jac)
    except np.linalg.LinAlgError:
        return None, None
    rss = float(np.sum(np.asarray(result.fun, dtype=float) ** 2))
    cov = jtj_inv * (rss / dof)
    stderr = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))
    return cov, stderr


def covariance_lookup(free_names: list[str], cov: np.ndarray | None) -> dict[tuple[str, str], float] | None:
    if cov is None:
        return None
    arr = np.asarray(cov, dtype=float)
    out: dict[tuple[str, str], float] = {}
    for i, ni in enumerate(free_names):
        for j, nj in enumerate(free_names):
            out[(ni, nj)] = float(arr[i, j])
    return out


def derive_model_parameters(model, params: dict[str, float], free_names: list[str], fixed_params: dict[str, float], cov: np.ndarray | None, axis: str | None = None) -> dict[str, Any]:
    cov_lookup = covariance_lookup(free_names, cov)
    try:
        return model.derived_parameters(params, covariance_deg=cov_lookup, fixed_names=set(fixed_params), axis=axis)
    except Exception as exc:
        return {"__error__": {"message": f"derived parameter computation failed: {exc}"}}


def analyze_jacobian(result, n_free_params: int, n_total_params: int) -> dict[str, Any]:
    if result.jac is None:
        return {"available": False, "n_total_params": int(n_total_params), "n_free_params": int(n_free_params)}
    jac = np.asarray(result.jac, dtype=float)
    if jac.shape[1] == 0:
        return {
            "available": True,
            "shape": [int(jac.shape[0]), int(jac.shape[1])],
            "rank": 0,
            "n_free_params": int(n_free_params),
            "n_total_params": int(n_total_params),
            "tolerance": 0.0,
            "singular_values": [],
            "condition_number": float("nan"),
            "warnings": [],
        }
    _, s, _ = np.linalg.svd(jac, full_matrices=False)
    tol = max(jac.shape) * np.finfo(float).eps * float(s[0])
    rank = int(np.sum(s > tol))
    cond = float(s[0] / s[rank - 1]) if rank > 0 and s[rank - 1] > 0 else float("inf")
    warnings: list[str] = []
    if rank < n_free_params:
        warnings.append(
            f"Jacobian rank deficient: rank={rank} < n_free_params={n_free_params}. Some free parameters are not independently identifiable."
        )
    elif cond > 1.0e8:
        warnings.append(f"Jacobian is ill-conditioned: condition_number={cond:.3e}.")
    return {
        "available": True,
        "shape": [int(jac.shape[0]), int(jac.shape[1])],
        "rank": rank,
        "n_free_params": int(n_free_params),
        "n_total_params": int(n_total_params),
        "tolerance": float(tol),
        "singular_values": [float(v) for v in s.tolist()],
        "condition_number": cond,
        "warnings": warnings,
    }


def apply_sigma_clip_blocks(df: pd.DataFrame, residual_daz_deg: np.ndarray, residual_del_deg: np.ndarray, sigma: float) -> pd.Series:
    r_arcsec = np.sqrt((residual_daz_deg * ARCSEC_PER_DEG) ** 2 + (residual_del_deg * ARCSEC_PER_DEG) ** 2)
    tmp = pd.DataFrame({"dt_cross_id": df["dt_cross_id"].values, "r_arcsec": r_arcsec})
    block = tmp.groupby("dt_cross_id", dropna=False)["r_arcsec"].median()
    thr = np.nanmedian(block.values) + sigma * mad_scale(block.values)
    bad_ids = set(block.index[block > thr].tolist())
    return df["dt_cross_id"].isin(bad_ids)


def bootstrap_parameters_joint(df: pd.DataFrame, opts: FitOptions, all_names: list[str], fixed_params: dict[str, float], delta_base_params: dict[str, float] | None = None) -> dict[str, dict[str, float]]:
    if opts.bootstrap <= 0:
        return {}
    rng = np.random.default_rng(opts.seed)
    block_ids = df["dt_cross_id"].astype(str).unique()
    samples: list[dict[str, float]] = []
    for _ in range(opts.bootstrap):
        chosen = rng.choice(block_ids, size=len(block_ids), replace=True)
        sub = pd.concat([df[df["dt_cross_id"].astype(str) == bid] for bid in chosen], ignore_index=True)
        core, result, _ = fit_once(sub, opts.model, "joint", opts.robust_loss, opts.f_scale, fixed_params, opts.fit_target, delta_base_params)
        samples.append(core.x_to_params(np.asarray(result.x, dtype=float)))
    out: dict[str, dict[str, float]] = {}
    for name in all_names:
        arr = np.asarray([s[name] for s in samples], dtype=float)
        out[name] = {
            "p16_deg": float(np.percentile(arr, 16)),
            "p50_deg": float(np.percentile(arr, 50)),
            "p84_deg": float(np.percentile(arr, 84)),
            "fixed": bool(name in fixed_params),
        }
    return out


def bootstrap_parameters_separate(df: pd.DataFrame, opts: FitOptions, cores: dict[str, FitCore], fixed_params: dict[str, float], delta_base_params: dict[str, float] | None = None) -> dict[str, Any]:
    if opts.bootstrap <= 0:
        return {}
    rng = np.random.default_rng(opts.seed)
    block_ids = df["dt_cross_id"].astype(str).unique()
    samples_by_axis: dict[str, list[dict[str, float]]] = {"dx": [], "dy": []}
    for _ in range(opts.bootstrap):
        chosen = rng.choice(block_ids, size=len(block_ids), replace=True)
        sub = pd.concat([df[df["dt_cross_id"].astype(str) == bid] for bid in chosen], ignore_index=True)
        core_map, result_map, _ = fit_once(sub, opts.model, "separate", opts.robust_loss, opts.f_scale, fixed_params, opts.fit_target, delta_base_params)
        for axis in ("dx", "dy"):
            samples_by_axis[axis].append(core_map[axis].x_to_params(np.asarray(result_map[axis].x, dtype=float)))
    out: dict[str, Any] = {}
    for axis in ("dx", "dy"):
        out[axis] = {}
        for name in cores[axis].all_names:
            arr = np.asarray([s.get(name, np.nan) for s in samples_by_axis[axis]], dtype=float)
            valid = np.isfinite(arr)
            if not np.any(valid):
                continue
            arrv = arr[valid]
            out[axis][name] = {
                "p16_deg": float(np.percentile(arrv, 16)),
                "p50_deg": float(np.percentile(arrv, 50)),
                "p84_deg": float(np.percentile(arrv, 84)),
                "fixed": bool(name in cores[axis].fixed_params),
                "active": bool(name in getattr(cores[axis].model, "axis_parameter_names")(axis)),
            }
    return out






def _preferred_display_unit(name: str, native: str) -> str:
    native = (native or "").strip()
    if native != "deg":
        return native
    if str(name).startswith("omega") or str(name) == "cor_p":
        return "deg"
    return "arcsec"

def _model_unit_map(model_name: str) -> dict[str, str]:
    model = resolve_model(model_name.split("__", 1)[0])
    return {p.name: p.unit for p in model.parameter_defs()}


def _format_param_display(model_name: str, name: str, value: float, stderr: float | None = None) -> tuple[str, str]:
    native = _model_unit_map(model_name).get(name, "deg")
    if native in {"", "1"}:
        val = f"{value:.6g}"
        if stderr is not None:
            val += f" ± {stderr:.6g}"
        return val, ""
    preferred = _preferred_display_unit(name, native)
    if native == "deg" and preferred == "arcsec":
        val = f"{value * ARCSEC_PER_DEG:.6f}"
        if stderr is not None:
            val += f" ± {stderr * ARCSEC_PER_DEG:.6g}"
        return val, "arcsec"
    if native == "deg" and preferred == "deg":
        val = f"{value:.9g}"
        if stderr is not None:
            val += f" ± {stderr:.6g}"
        return val, "deg"
    val = f"{value:.6g}"
    if stderr is not None:
        val += f" ± {stderr:.6g}"
    return val, native


def _external_param_entry(model_name: str, name: str, value: float, stderr: float | None, fixed: bool) -> dict[str, Any]:
    native = _model_unit_map(model_name).get(name, "deg")
    if native in {"", "1"}:
        return {"value": float(value), "stderr": None if stderr is None else float(stderr), "fixed": bool(fixed), "unit": native or "1"}
    return {"value": float(value), "stderr": None if stderr is None else float(stderr), "fixed": bool(fixed), "unit": native}


def _external_param_key(model_name: str, name: str) -> str:
    native = str(_model_unit_map(model_name).get(name, "deg") or "").strip()
    if native in {"", "1"}:
        return name
    return f"{name}[{native}]"
def _fixed_meta_dict(model_name: str, fixed_params: dict[str, float], fixed_sources: dict[str, str]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in fixed_params.items():
        entry = _external_param_entry(model_name, k, v, None, True)
        entry["source"] = fixed_sources.get(k, "")
        out[k] = entry
    return out
def build_parameter_text_joint(
    path: Path,
    fit_target: str,
    model_name: str,
    params: dict[str, float],
    stderr: dict[str, float] | None,
    fixed_params: dict[str, float],
    fixed_sources: dict[str, str],
    free_names: list[str],
    jacobian_diagnostics: dict[str, Any] | None = None,
    derived_parameters: dict[str, Any] | None = None,
) -> None:
    lines = [
        f"Fit target : {fit_target}",
        f"Solve mode : joint",
        f"Model      : {model_name}",
        "",
    ]
    if fixed_params:
        lines.append("Fixed parameters")
        for k in sorted(fixed_params):
            vtxt, unit_txt = _format_param_display(model_name, k, fixed_params[k], None)
            src = fixed_sources.get(k, "")
            lines.append(f"  {k} = {vtxt}{(' ' + unit_txt) if unit_txt else ''}  ({src})")
        lines.append("")
    lines.append("Fitted parameters" if free_names else "No fitted parameters")
    for k in free_names:
        stderr_k = None if not stderr or k not in stderr else stderr[k]
        vtxt, unit_txt = _format_param_display(model_name, k, params[k], stderr_k)
        if fit_target == "delta":
            lines.append(f"  {k} += {vtxt}{(' ' + unit_txt) if unit_txt else ''}")
        else:
            lines.append(f"  {k} = {vtxt}{(' ' + unit_txt) if unit_txt else ''}")
    derived_parameters = derived_parameters or {}
    printable_derived = {k: v for k, v in derived_parameters.items() if not str(k).startswith("__")}
    if printable_derived:
        lines.extend(["", "Derived parameters"])
        for name, info in printable_derived.items():
            amp = info.get("amplitude_arcsec")
            amp_err = info.get("amplitude_stderr_arcsec")
            phi = info.get("phase_deg")
            phi_err = info.get("phase_stderr_deg")
            basis = info.get("basis", "")
            if amp is not None:
                msg = f"  {name}: amplitude = {amp:.6f} arcsec"
                if amp_err is not None:
                    msg += f" ± {amp_err:.6g}"
                if phi is None:
                    msg += "; phase = undefined"
                else:
                    msg += f"; phase = {phi:.6f} deg"
                    if phi_err is not None:
                        msg += f" ± {phi_err:.6g}"
                lines.append(msg)
            else:
                lines.append(f"  {name}: {info}")
            if basis:
                lines.append(f"    basis: {basis}")
    if jacobian_diagnostics and jacobian_diagnostics.get("warnings"):
        lines.extend(["", "Warnings"])
        for msg in jacobian_diagnostics["warnings"]:
            lines.append(f"  - {msg}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parameter_text_separate(
    path: Path,
    fit_target: str,
    model_name: str,
    axis_params: dict[str, dict[str, float]],
    axis_stderr: dict[str, dict[str, float] | None],
    fixed_params: dict[str, float],
    fixed_sources: dict[str, str],
    axis_free_names: dict[str, list[str]],
    axis_jacobian_diagnostics: dict[str, dict[str, Any]],
    derived_parameters_by_axis: dict[str, dict[str, Any]] | None = None,
) -> None:
    lines = [
        f"Fit target : {fit_target}",
        f"Solve mode : separate",
        f"Model      : {model_name}",
        "",
        "Separate fit produces independent dx-side and dy-side parameter estimates.",
        "Shared parameter names may therefore take different fitted values in the two axis sections below.",
        "",
    ]
    if fixed_params:
        lines.append("User-requested fixed parameters")
        for k in sorted(fixed_params):
            src = fixed_sources.get(k, "")
            vtxt, unit_txt = _format_param_display(model_name, k, fixed_params[k], None)
            lines.append(f"  {k} = {vtxt}{(' ' + unit_txt) if unit_txt else ''}  ({src})")
        lines.append("")
    for axis in ("dx", "dy"):
        label = "dx-axis fit" if axis == "dx" else "dy-axis fit"
        lines.append(label)
        lines.append("-" * len(label))
        active_names = sorted(axis_params[axis].keys())
        fixed_axis = sorted([k for k in active_names if k in fixed_params])
        if fixed_axis:
            lines.append("  Fixed in this axis:")
            for k in fixed_axis:
                vtxt, unit_txt = _format_param_display(model_name, k, fixed_params[k], None)
                lines.append(f"    {k} = {vtxt}{(' ' + unit_txt) if unit_txt else ''}")
        free_names = axis_free_names[axis]
        if free_names:
            lines.append("  Fitted in this axis:")
            for k in free_names:
                stderr_map = axis_stderr.get(axis) or {}
                stderr_k = None if k not in stderr_map else stderr_map[k]
                vtxt, unit_txt = _format_param_display(model_name, k, axis_params[axis][k], stderr_k)
                op = "+=" if fit_target == "delta" else "="
                lines.append(f"    {k} {op} {vtxt}{(' ' + unit_txt) if unit_txt else ''}")
        else:
            lines.append("  No free parameters in this axis.")
        derived_axis = (derived_parameters_by_axis or {}).get(axis, {})
        printable_derived = {k: v for k, v in derived_axis.items() if not str(k).startswith("__")}
        if printable_derived:
            lines.append("  Derived parameters:")
            for name, info in printable_derived.items():
                amp = info.get("amplitude_arcsec")
                amp_err = info.get("amplitude_stderr_arcsec")
                phi = info.get("phase_deg")
                phi_err = info.get("phase_stderr_deg")
                basis = info.get("basis", "")
                msg = f"    {name}: amplitude = {amp:.6f} arcsec" if amp is not None else f"    {name}: {info}"
                if amp is not None:
                    if amp_err is not None:
                        msg += f" ± {amp_err:.6g}"
                    if phi is None:
                        msg += "; phase = undefined"
                    else:
                        msg += f"; phase = {phi:.6f} deg"
                        if phi_err is not None:
                            msg += f" ± {phi_err:.6g}"
                lines.append(msg)
                if basis:
                    lines.append(f"      basis: {basis}")
        warns = axis_jacobian_diagnostics.get(axis, {}).get("warnings", [])
        if warns:
            lines.append("  Warnings:")
            for msg in warns:
                lines.append(f"    - {msg}")
        lines.append("")
    shared = sorted(set(axis_params["dx"]).intersection(axis_params["dy"]))
    shared_lines = []
    for k in shared:
        vd = axis_params["dx"].get(k)
        ve = axis_params["dy"].get(k)
        if vd is None or ve is None or not np.isfinite(vd) or not np.isfinite(ve):
            continue
        if abs(vd - ve) == 0:
            continue
        vd_txt, vd_unit = _format_param_display(model_name, k, vd, None)
        ve_txt, ve_unit = _format_param_display(model_name, k, ve, None)
        diff_txt, diff_unit = _format_param_display(model_name, k, abs(vd - ve), None)
        shared_lines.append(
            f"  {k}: dx-side={vd_txt}{(' ' + vd_unit) if vd_unit else ''}, dy-side={ve_txt}{(' ' + ve_unit) if ve_unit else ''}, |difference|={diff_txt}{(' ' + diff_unit) if diff_unit else ''}"
        )
    if shared_lines:
        lines.append("Shared-name comparison")
        lines.extend(shared_lines)
        lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_residual_table_joint(
    df: pd.DataFrame,
    model_name: str,
    best_params: dict[str, float],
    fit_target: str,
    core: FitCore | None = None,
    x_free: np.ndarray | None = None,
) -> pd.DataFrame:
    model = resolve_model(model_name)
    out = df.copy()
    use_core_predict = core is not None and core.predict_fun is not None and x_free is not None
    if use_core_predict:
        pred_joint = np.asarray(core.predict_fun(np.asarray(x_free, dtype=float)), dtype=float)
        n = len(out)
        if pred_joint.size == 2 * n:
            pred_daz_deg = pred_joint[:n]
            pred_del_deg = pred_joint[n:]
        else:
            pred_daz_deg = np.asarray(model.predict_dx_deg(best_params, out["az_deg"], out["el_deg"]), dtype=float)
            pred_del_deg = np.asarray(model.predict_dy_deg(best_params, out["az_deg"], out["el_deg"]), dtype=float)
    else:
        pred_daz_deg = np.asarray(model.predict_dx_deg(best_params, out["az_deg"], out["el_deg"]), dtype=float)
        pred_del_deg = np.asarray(model.predict_dy_deg(best_params, out["az_deg"], out["el_deg"]), dtype=float)
    out["predicted_dx_arcsec"] = pred_daz_deg * ARCSEC_PER_DEG
    out["predicted_dy_arcsec"] = pred_del_deg * ARCSEC_PER_DEG
    out["residual_dx_arcsec"] = out["target_dx_arcsec"].astype(float) - out["predicted_dx_arcsec"]
    out["residual_dy_arcsec"] = out["target_dy_arcsec"].astype(float) - out["predicted_dy_arcsec"]
    out["radial_residual_arcsec"] = np.sqrt(out["residual_dx_arcsec"] ** 2 + out["residual_dy_arcsec"] ** 2)
    if fit_target == "absolute" and "used_model_dx_arcsec" in out.columns:
        out["delta_to_used_dx_arcsec"] = out["predicted_dx_arcsec"] - out["used_model_dx_arcsec"].astype(float)
        out["delta_to_used_dy_arcsec"] = out["predicted_dy_arcsec"] - out["used_model_dy_arcsec"].astype(float)
    elif fit_target == "delta" and "used_model_dx_arcsec" in out.columns:
        out["delta_to_used_dx_arcsec"] = out["predicted_dx_arcsec"]
        out["delta_to_used_dy_arcsec"] = out["predicted_dy_arcsec"]
        out["predicted_total_dx_arcsec"] = out["used_model_dx_arcsec"].astype(float) + out["predicted_dx_arcsec"]
        out["predicted_total_dy_arcsec"] = out["used_model_dy_arcsec"].astype(float) + out["predicted_dy_arcsec"]
    return out


def make_residual_table_separate(
    df: pd.DataFrame,
    model_name: str,
    best_params_by_axis: dict[str, dict[str, float]],
    fit_target: str,
    core_map: dict[str, FitCore] | None = None,
    x_map: dict[str, np.ndarray] | None = None,
) -> pd.DataFrame:
    model = resolve_model(model_name)
    out = df.copy()
    use_core_predict = core_map is not None and x_map is not None and all(core_map[a].predict_fun is not None for a in ("dx", "dy"))
    if use_core_predict:
        pred_daz_deg = np.asarray(core_map["dx"].predict_fun(np.asarray(x_map["dx"], dtype=float)), dtype=float)
        pred_del_deg = np.asarray(core_map["dy"].predict_fun(np.asarray(x_map["dy"], dtype=float)), dtype=float)
        if pred_daz_deg.size != len(out) or pred_del_deg.size != len(out):
            pred_daz_deg = np.asarray(model.predict_dx_deg(best_params_by_axis["dx"], out["az_deg"], out["el_deg"]), dtype=float)
            pred_del_deg = np.asarray(model.predict_dy_deg(best_params_by_axis["dy"], out["az_deg"], out["el_deg"]), dtype=float)
    else:
        pred_daz_deg = np.asarray(model.predict_dx_deg(best_params_by_axis["dx"], out["az_deg"], out["el_deg"]), dtype=float)
        pred_del_deg = np.asarray(model.predict_dy_deg(best_params_by_axis["dy"], out["az_deg"], out["el_deg"]), dtype=float)
    out["predicted_dx_arcsec"] = pred_daz_deg * ARCSEC_PER_DEG
    out["predicted_dy_arcsec"] = pred_del_deg * ARCSEC_PER_DEG
    out["residual_dx_arcsec"] = out["target_dx_arcsec"].astype(float) - out["predicted_dx_arcsec"]
    out["residual_dy_arcsec"] = out["target_dy_arcsec"].astype(float) - out["predicted_dy_arcsec"]
    out["radial_residual_arcsec"] = np.sqrt(out["residual_dx_arcsec"] ** 2 + out["residual_dy_arcsec"] ** 2)
    if fit_target == "absolute" and "used_model_dx_arcsec" in out.columns:
        out["delta_to_used_dx_arcsec"] = out["predicted_dx_arcsec"] - out["used_model_dx_arcsec"].astype(float)
        out["delta_to_used_dy_arcsec"] = out["predicted_dy_arcsec"] - out["used_model_dy_arcsec"].astype(float)
    elif fit_target == "delta" and "used_model_dx_arcsec" in out.columns:
        out["delta_to_used_dx_arcsec"] = out["predicted_dx_arcsec"]
        out["delta_to_used_dy_arcsec"] = out["predicted_dy_arcsec"]
        out["predicted_total_dx_arcsec"] = out["used_model_dx_arcsec"].astype(float) + out["predicted_dx_arcsec"]
        out["predicted_total_dy_arcsec"] = out["used_model_dy_arcsec"].astype(float) + out["predicted_dy_arcsec"]
    return out


def run_fit(opts: FitOptions) -> dict[str, Any]:
    outdir = Path(opts.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = read_csv(opts.input_csv)
    for col in EXPECTED_COLS:
        if col not in df.columns:
            raise KeyError(f"Missing required column: {col}")
    df = prepare_targets(df, opts.fit_target)
    fixed_params, fixed_sources = resolve_fixed_params(opts.model, opts.fixed_param_specs, opts.fixed_param_file)
    delta_base_params = resolve_delta_base_params(df, opts.model, opts.fixed_param_file) if opts.fit_target == "delta" else None

    if opts.solve_mode == "joint":
        core, result, scale = fit_once(df, opts.model, "joint", opts.robust_loss, opts.f_scale, fixed_params, opts.fit_target, delta_base_params)
        residual_table = make_residual_table_joint(df, opts.model, core.x_to_params(np.asarray(result.x, dtype=float)), opts.fit_target, core=core, x_free=np.asarray(result.x, dtype=float))
        if opts.clip_sigma > 0 and opts.clip_max_iter > 0:
            current = df.copy()
            for _ in range(opts.clip_max_iter):
                residual_table_iter = make_residual_table_joint(current, opts.model, core.x_to_params(np.asarray(result.x, dtype=float)), opts.fit_target, core=core, x_free=np.asarray(result.x, dtype=float))
                bad_iter = apply_sigma_clip_blocks(
                    current,
                    residual_table_iter["residual_dx_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
                    residual_table_iter["residual_dy_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
                    opts.clip_sigma,
                )
                if not bad_iter.any():
                    break
                current = current.loc[~bad_iter].reset_index(drop=True)
                core, result, scale = fit_once(current, opts.model, "joint", opts.robust_loss, opts.f_scale, fixed_params, opts.fit_target, delta_base_params)
            full = make_residual_table_joint(df, opts.model, core.x_to_params(np.asarray(result.x, dtype=float)), opts.fit_target, core=core, x_free=np.asarray(result.x, dtype=float))
            bad_full = apply_sigma_clip_blocks(
                df,
                full["residual_dx_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
                full["residual_dy_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
                opts.clip_sigma,
            )
            full["is_outlier_block"] = bad_full.to_numpy(dtype=bool)
            residual_table = full
        else:
            residual_table["is_outlier_block"] = False

        best_params = core.x_to_params(np.asarray(result.x, dtype=float))
        cov, stderr_arr = compute_covariance(result, n_obs=result.fun.size, n_params=len(core.free_names))
        stderr_dict = {k: float(v) for k, v in zip(core.free_names, stderr_arr)} if stderr_arr is not None else None
        jacobian_diagnostics = analyze_jacobian(result, len(core.free_names), len(core.all_names))
        derived_parameters = derive_model_parameters(core.model, best_params, core.free_names, fixed_params, cov, axis=None)
        bootstrap_summary = bootstrap_parameters_joint(df.loc[~residual_table["is_outlier_block"]].reset_index(drop=True), opts, core.all_names, fixed_params, delta_base_params)

        rms = radial_rms_arcsec(
            residual_table.loc[~residual_table["is_outlier_block"], "residual_dx_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
            residual_table.loc[~residual_table["is_outlier_block"], "residual_dy_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
        )

        residual_table["robust_weight"] = 1.0
        residual_csv = outdir / "residuals.csv"
        outlier_csv = outdir / "outliers.csv"
        write_csv(residual_table, residual_csv)
        write_csv(residual_table.loc[residual_table["is_outlier_block"]], outlier_csv)

        output_model_label = infer_output_model_label(df, opts.model)
        param_template_file = _template_for_full_param_output(df, opts.fixed_param_file)
        fixed_template_file = _template_for_fixed_param_output(opts.fixed_param_file)
        key_style_map: dict[str, dict[str, Any]] = {}
        if param_template_file and Path(param_template_file).exists():
            try:
                key_style_map = parse_param_key_styles(param_template_file, section_names=("pointing_params", "fixed_params"))
            except Exception:
                key_style_map = {}
        fixed_key_style_map: dict[str, dict[str, Any]] = {}
        if fixed_template_file and Path(fixed_template_file).exists():
            try:
                fixed_key_style_map = parse_param_key_styles(fixed_template_file, section_names=("fixed_params", "pointing_params"))
            except Exception:
                fixed_key_style_map = {}
        param_out = outdir / ("absolute_params.toml" if opts.fit_target == "absolute" else "delta_params.toml")
        _write_param_output_with_optional_template(
            param_out,
            template_path=param_template_file,
            params_deg=best_params,
            model_name=output_model_label,
            section_name="pointing_params",
            key_style_map=key_style_map,
        )
        if fixed_params:
            _write_param_output_with_optional_template(
                outdir / "resolved_fixed_params.toml",
                template_path=fixed_template_file,
                params_deg=fixed_params,
                model_name=output_model_label,
                section_name="fixed_params",
                key_style_map=fixed_key_style_map or key_style_map,
            )

        build_parameter_text_joint(
            outdir / "parameter_updates.txt",
            opts.fit_target,
            opts.model,
            best_params,
            stderr_dict,
            fixed_params,
            fixed_sources,
            core.free_names,
            jacobian_diagnostics,
            derived_parameters,
        )

        result_json: dict[str, Any] = {
            "tool": {"name": "pointing-fit", "version": __version__},
            "run": {
                "fit_target": opts.fit_target,
                "solve_mode": opts.solve_mode,
                "model_name": opts.model,
                "robust_loss": opts.robust_loss,
                "f_scale": scale,
                "bootstrap": opts.bootstrap,
                "absolute_target_rule": "used_model_plus_delta_err",
                "fixed_param_file": opts.fixed_param_file or "",
                "fixed_param_specs": list(opts.fixed_param_specs or []),
                "delta_linearized": bool(core.linearized_delta),
            },
            "input": {
                "input_csv": opts.input_csv,
                "n_rows": int(len(df)),
                "dataset_ids": sorted(df["dataset_id"].astype(str).unique().tolist()),
            },
            "quality": {
                "n_used": int((~residual_table["is_outlier_block"]).sum()),
                "n_outlier_blocks": int(residual_table["is_outlier_block"].sum()),
                "residual_rms_arcsec": rms,
                "cost": float(result.cost),
                "success": bool(result.success),
                "message": str(result.message),
            },
            "jacobian_diagnostics": jacobian_diagnostics,
            "free_parameter_names": list(core.free_names),
            "fixed_parameters": _fixed_meta_dict(opts.model, fixed_params, fixed_sources),
            "parameters_internal": {
                k: {
                    "value_deg": float(best_params[k]),
                    "stderr_deg": None if stderr_dict is None or k not in stderr_dict else float(stderr_dict[k]),
                    "fixed": bool(k in fixed_params),
                }
                for k in core.all_names
            },
            "parameters_external": {
                _external_param_key(opts.model, k): _external_param_entry(opts.model, k, best_params[k], None if stderr_dict is None or k not in stderr_dict else stderr_dict[k], k in fixed_params)
                for k in core.all_names
            },
            "derived_parameters": derived_parameters,
            "bootstrap_summary": bootstrap_summary,
            "files": {
                "residuals_csv": str(residual_csv),
                "outliers_csv": str(outlier_csv),
                "param_toml": str(param_out),
                "resolved_fixed_params_toml": str(outdir / "resolved_fixed_params.toml") if fixed_params else "",
            },
        }
        if cov is not None:
            result_json["covariance"] = {"parameter_order": list(core.free_names), "matrix": np.asarray(cov, dtype=float).tolist()}
        if core.base_params is not None:
            result_json["delta_base_parameters_internal"] = {k: float(v) for k, v in core.base_params.items()}
        dump_json(outdir / "fit_result.json", result_json)
        return result_json

    # separate mode
    core_map, result_map, scale_map = fit_once(df, opts.model, "separate", opts.robust_loss, opts.f_scale, fixed_params, opts.fit_target, delta_base_params)
    assert isinstance(core_map, dict)
    residual_table = make_residual_table_separate(
        df,
        opts.model,
        {
            "dx": core_map["dx"].x_to_params(np.asarray(result_map["dx"].x, dtype=float)),
            "dy": core_map["dy"].x_to_params(np.asarray(result_map["dy"].x, dtype=float)),
        },
        opts.fit_target,
        core_map=core_map,
        x_map={axis: np.asarray(result_map[axis].x, dtype=float) for axis in ("dx", "dy")},
    )
    if opts.clip_sigma > 0 and opts.clip_max_iter > 0:
        current = df.copy()
        for _ in range(opts.clip_max_iter):
            residual_table_iter = make_residual_table_separate(
                current,
                opts.model,
                {
                    "dx": core_map["dx"].x_to_params(np.asarray(result_map["dx"].x, dtype=float)),
                    "dy": core_map["dy"].x_to_params(np.asarray(result_map["dy"].x, dtype=float)),
                },
                opts.fit_target,
                core_map=core_map,
                x_map={axis: np.asarray(result_map[axis].x, dtype=float) for axis in ("dx", "dy")},
            )
            bad_iter = apply_sigma_clip_blocks(
                current,
                residual_table_iter["residual_dx_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
                residual_table_iter["residual_dy_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
                opts.clip_sigma,
            )
            if not bad_iter.any():
                break
            current = current.loc[~bad_iter].reset_index(drop=True)
            core_map, result_map, scale_map = fit_once(current, opts.model, "separate", opts.robust_loss, opts.f_scale, fixed_params, opts.fit_target, delta_base_params)
        full = make_residual_table_separate(
            df,
            opts.model,
            {
                "dx": core_map["dx"].x_to_params(np.asarray(result_map["dx"].x, dtype=float)),
                "dy": core_map["dy"].x_to_params(np.asarray(result_map["dy"].x, dtype=float)),
            },
            opts.fit_target,
            core_map=core_map,
            x_map={axis: np.asarray(result_map[axis].x, dtype=float) for axis in ("dx", "dy")},
        )
        bad_full = apply_sigma_clip_blocks(
            df,
            full["residual_dx_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
            full["residual_dy_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
            opts.clip_sigma,
        )
        full["is_outlier_block"] = bad_full.to_numpy(dtype=bool)
        residual_table = full
    else:
        residual_table["is_outlier_block"] = False

    best_params_by_axis = {
        axis: core_map[axis].x_to_params(np.asarray(result_map[axis].x, dtype=float)) for axis in ("dx", "dy")
    }
    cov_by_axis: dict[str, np.ndarray | None] = {}
    stderr_by_axis: dict[str, dict[str, float] | None] = {}
    jacobian_by_axis: dict[str, dict[str, Any]] = {}
    for axis in ("dx", "dy"):
        cov_axis, stderr_arr_axis = compute_covariance(result_map[axis], n_obs=result_map[axis].fun.size, n_params=len(core_map[axis].free_names))
        cov_by_axis[axis] = cov_axis
        stderr_by_axis[axis] = {k: float(v) for k, v in zip(core_map[axis].free_names, stderr_arr_axis)} if stderr_arr_axis is not None else None
        jacobian_by_axis[axis] = analyze_jacobian(result_map[axis], len(core_map[axis].free_names), len(core_map[axis].all_names))

    derived_by_axis = {
        axis: derive_model_parameters(core_map[axis].model, best_params_by_axis[axis], core_map[axis].free_names, core_map[axis].fixed_params, cov_by_axis[axis], axis=axis)
        for axis in ("dx", "dy")
    }

    bootstrap_summary = bootstrap_parameters_separate(
        df.loc[~residual_table["is_outlier_block"]].reset_index(drop=True),
        opts,
        core_map,
        fixed_params,
        delta_base_params,
    )

    rms = radial_rms_arcsec(
        residual_table.loc[~residual_table["is_outlier_block"], "residual_dx_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
        residual_table.loc[~residual_table["is_outlier_block"], "residual_dy_arcsec"].to_numpy(dtype=float) / ARCSEC_PER_DEG,
    )

    residual_table["robust_weight"] = 1.0
    residual_csv = outdir / "residuals.csv"
    outlier_csv = outdir / "outliers.csv"
    write_csv(residual_table, residual_csv)
    write_csv(residual_table.loc[residual_table["is_outlier_block"]], outlier_csv)

    output_model_label = infer_output_model_label(df, opts.model)
    param_template_file = _template_for_full_param_output(df, opts.fixed_param_file)
    fixed_template_file = _template_for_fixed_param_output(opts.fixed_param_file)
    key_style_map: dict[str, dict[str, Any]] = {}
    if param_template_file and Path(param_template_file).exists():
        try:
            key_style_map = parse_param_key_styles(param_template_file, section_names=("pointing_params", "fixed_params"))
        except Exception:
            key_style_map = {}
    fixed_key_style_map: dict[str, dict[str, Any]] = {}
    if fixed_template_file and Path(fixed_template_file).exists():
        try:
            fixed_key_style_map = parse_param_key_styles(fixed_template_file, section_names=("fixed_params", "pointing_params"))
        except Exception:
            fixed_key_style_map = {}
    prefix = "absolute" if opts.fit_target == "absolute" else "delta"
    axis_param_files = {}
    for axis in ("dx", "dy"):
        param_out = outdir / f"{prefix}_params_{axis}.toml"
        _write_param_output_with_optional_template(
            param_out,
            template_path=param_template_file,
            params_deg=best_params_by_axis[axis],
            model_name=output_model_label,
            section_name="pointing_params",
            key_style_map=key_style_map,
        )
        axis_param_files[axis] = str(param_out)

    if fixed_params:
        _write_param_output_with_optional_template(
            outdir / "resolved_fixed_params.toml",
            template_path=fixed_template_file,
            params_deg=fixed_params,
            model_name=output_model_label,
            section_name="fixed_params",
            key_style_map=fixed_key_style_map or key_style_map,
        )

    build_parameter_text_separate(
        outdir / "parameter_updates.txt",
        opts.fit_target,
        opts.model,
        axis_params=best_params_by_axis,
        axis_stderr=stderr_by_axis,
        fixed_params=fixed_params,
        fixed_sources=fixed_sources,
        axis_free_names={axis: list(core_map[axis].free_names) for axis in ("dx", "dy")},
        axis_jacobian_diagnostics=jacobian_by_axis,
        derived_parameters_by_axis=derived_by_axis,
    )

    shared_parameter_difference: dict[str, Any] = {}
    for name in sorted(set(best_params_by_axis["dx"]).intersection(best_params_by_axis["dy"])):
        vd = best_params_by_axis["dx"][name]
        ve = best_params_by_axis["dy"][name]
        entry_dx = _external_param_entry(opts.model, name, float(vd), None, False)
        entry_dy = _external_param_entry(opts.model, name, float(ve), None, False)
        diff_entry = _external_param_entry(opts.model, name, float(vd - ve), None, False)
        shared_parameter_difference[name] = {
            "dx_value_deg": float(vd),
            "dy_value_deg": float(ve),
            "difference_deg": float(vd - ve),
            "display_unit": diff_entry.get("unit", ""),
            "dx_value_display": entry_dx.get("value"),
            "dy_value_display": entry_dy.get("value"),
            "difference_display": diff_entry.get("value"),
        }

    derived_difference: dict[str, Any] = {}
    common_derived = sorted(set(derived_by_axis["dx"]).intersection(derived_by_axis["dy"]))
    for name in common_derived:
        info_d = derived_by_axis["dx"].get(name, {})
        info_e = derived_by_axis["dy"].get(name, {})
        if not isinstance(info_d, dict) or not isinstance(info_e, dict):
            continue
        if "amplitude_deg" not in info_d or "amplitude_deg" not in info_e:
            continue
        amp_d = info_d.get("amplitude_deg")
        amp_e = info_e.get("amplitude_deg")
        phi_d = info_d.get("phase_deg")
        phi_e = info_e.get("phase_deg")
        phase_diff = None
        if phi_d is not None and phi_e is not None:
            phase_diff = float(((float(phi_d) - float(phi_e) + 180.0) % 360.0) - 180.0)
        derived_difference[name] = {
            "dx_amplitude_arcsec": None if amp_d is None else float(float(amp_d) * ARCSEC_PER_DEG),
            "dy_amplitude_arcsec": None if amp_e is None else float(float(amp_e) * ARCSEC_PER_DEG),
            "difference_amplitude_arcsec": None if amp_d is None or amp_e is None else float((float(amp_d) - float(amp_e)) * ARCSEC_PER_DEG),
            "dx_phase_deg": None if phi_d is None else float(phi_d),
            "dy_phase_deg": None if phi_e is None else float(phi_e),
            "difference_phase_deg_wrapped": phase_diff,
        }

    quality = {
        "n_used": int((~residual_table["is_outlier_block"]).sum()),
        "n_outlier_blocks": int(residual_table["is_outlier_block"].sum()),
        "residual_rms_arcsec": rms,
        "success": bool(result_map["dx"].success and result_map["dy"].success),
        "message": f"dx: {result_map['dx'].message} | dy: {result_map['dy'].message}",
        "cost_by_axis": {axis: float(result_map[axis].cost) for axis in ("dx", "dy")},
    }

    result_json = {
        "tool": {"name": "pointing-fit", "version": __version__},
        "run": {
            "fit_target": opts.fit_target,
            "solve_mode": opts.solve_mode,
            "model_name": opts.model,
            "robust_loss": opts.robust_loss,
            "f_scale": {axis: float(scale_map[axis]) for axis in ("dx", "dy")},
            "bootstrap": opts.bootstrap,
            "fixed_param_file": opts.fixed_param_file or "",
            "fixed_param_specs": list(opts.fixed_param_specs or []),
            "delta_linearized_by_axis": {axis: bool(core_map[axis].linearized_delta) for axis in ("dx", "dy")},
        },
        "input": {
            "input_csv": opts.input_csv,
            "n_rows": int(len(df)),
            "dataset_ids": sorted(df["dataset_id"].astype(str).unique().tolist()),
        },
        "quality": quality,
        "jacobian_diagnostics_by_axis": jacobian_by_axis,
        "free_parameter_names_by_axis": {axis: list(core_map[axis].free_names) for axis in ("dx", "dy")},
        "fixed_parameters": _fixed_meta_dict(opts.model, fixed_params, fixed_sources),
        "parameters_by_axis_internal": {
            axis: {
                k: {
                    "value_deg": float(best_params_by_axis[axis][k]),
                    "stderr_deg": None if stderr_by_axis[axis] is None or k not in (stderr_by_axis[axis] or {}) else float(stderr_by_axis[axis][k]),
                    "fixed": bool(k in core_map[axis].fixed_params),
                    "active": bool(k in core_map[axis].model.axis_parameter_names(axis)),
                }
                for k in core_map[axis].all_names
            }
            for axis in ("dx", "dy")
        },
        "parameters_by_axis_external": {
            axis: {
                _external_param_key(opts.model, k): _external_param_entry(
                    opts.model,
                    k,
                    best_params_by_axis[axis][k],
                    None if stderr_by_axis[axis] is None or k not in (stderr_by_axis[axis] or {}) else stderr_by_axis[axis][k],
                    k in core_map[axis].fixed_params,
                )
                for k in core_map[axis].all_names
            }
            for axis in ("dx", "dy")
        },
        "derived_parameters_by_axis": derived_by_axis,
        "shared_parameter_difference": shared_parameter_difference,
        "derived_parameter_difference": derived_difference,
        "bootstrap_summary_by_axis": bootstrap_summary,
        "files": {
            "residuals_csv": str(residual_csv),
            "outliers_csv": str(outlier_csv),
            "param_toml_by_axis": axis_param_files,
            "resolved_fixed_params_toml": str(outdir / "resolved_fixed_params.toml") if fixed_params else "",
        },
    }
    cov_json = {}
    for axis in ("dx", "dy"):
        if cov_by_axis[axis] is not None:
            cov_json[axis] = {"parameter_order": list(core_map[axis].free_names), "matrix": np.asarray(cov_by_axis[axis], dtype=float).tolist()}
    if cov_json:
        result_json["covariance_by_axis"] = cov_json
    base_by_axis = {axis: core_map[axis].base_params for axis in ("dx", "dy") if core_map[axis].base_params is not None}
    if base_by_axis:
        result_json["delta_base_parameters_internal_by_axis"] = {
            axis: {k: float(v) for k, v in base.items()} for axis, base in base_by_axis.items()
        }
    dump_json(outdir / "fit_result.json", result_json)
    return result_json
