from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from .io_utils import read_csv


@dataclass
class ReportOptions:
    input_csv: str
    outdir: str
    title: str = "pointing-fit report"


def _load_fit_result_near(csv_path: Path):
    candidate = csv_path.parent / "fit_result.json"
    if not candidate.exists():
        return None
    try:
        return json.loads(candidate.read_text(encoding="utf-8"))
    except Exception:
        return None


def _derived_lines(fit_result: dict | None) -> list[str]:
    if not fit_result:
        return []
    lines: list[str] = []
    if fit_result.get("run", {}).get("solve_mode") == "joint":
        derived = fit_result.get("derived_parameters", {}) or {}
        printable = {k: v for k, v in derived.items() if not str(k).startswith("__")}
        if printable:
            lines.extend(["", "## Derived parameters"])
            for name, info in printable.items():
                amp = info.get("amplitude_arcsec")
                amp_err = info.get("amplitude_stderr_arcsec")
                phi = info.get("phase_deg")
                phi_err = info.get("phase_stderr_deg")
                basis = info.get("basis", "")
                line = f"- `{name}`: amplitude = {amp:.6f} arcsec" if amp is not None else f"- `{name}`: {info}"
                if amp is not None:
                    if amp_err is not None:
                        line += f" ± {amp_err:.6g}"
                    if phi is None:
                        line += "; phase = undefined"
                    else:
                        line += f"; phase = {phi:.6f} deg"
                        if phi_err is not None:
                            line += f" ± {phi_err:.6g}"
                lines.append(line)
                if basis:
                    lines.append(f"  - basis: `{basis}`")
    else:
        derived_by_axis = fit_result.get("derived_parameters_by_axis", {}) or {}
        if derived_by_axis:
            lines.extend(["", "## Derived parameters by axis"])
            for axis in ("dx", "dy"):
                derived = {k: v for k, v in (derived_by_axis.get(axis, {}) or {}).items() if not str(k).startswith("__")}
                if not derived:
                    continue
                lines.append(f"### {axis}")
                for name, info in derived.items():
                    amp = info.get("amplitude_arcsec")
                    amp_err = info.get("amplitude_stderr_arcsec")
                    phi = info.get("phase_deg")
                    phi_err = info.get("phase_stderr_deg")
                    basis = info.get("basis", "")
                    line = f"- `{name}`: amplitude = {amp:.6f} arcsec" if amp is not None else f"- `{name}`: {info}"
                    if amp is not None:
                        if amp_err is not None:
                            line += f" ± {amp_err:.6g}"
                        if phi is None:
                            line += "; phase = undefined"
                        else:
                            line += f"; phase = {phi:.6f} deg"
                            if phi_err is not None:
                                line += f" ± {phi_err:.6g}"
                    lines.append(line)
                    if basis:
                        lines.append(f"  - basis: `{basis}`")
            diff = fit_result.get("derived_parameter_difference", {}) or {}
            if diff:
                lines.extend(["", "## Derived parameter axis difference"])
                for name, info in diff.items():
                    lines.append(
                        f"- `{name}`: Δamplitude = {info.get('difference_amplitude_arcsec')} arcsec, Δphase = {info.get('difference_phase_deg_wrapped')} deg"
                    )
    return lines


def _nan_limits(*arrays: np.ndarray) -> tuple[float, float]:
    vals = []
    for arr in arrays:
        a = np.asarray(arr, dtype=float).ravel()
        a = a[np.isfinite(a)]
        if a.size:
            vals.append(a)
    if not vals:
        return -1.0, 1.0
    merged = np.concatenate(vals)
    vmin = float(np.nanmin(merged))
    vmax = float(np.nanmax(merged))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        return -1.0, 1.0
    if abs(vmax - vmin) < 1e-12:
        pad = max(abs(vmax), 1.0) * 1.0e-6 + 1.0e-6
        return vmin - pad, vmax + pad
    pad = 0.08 * (vmax - vmin)
    return vmin - pad, vmax + pad


def _symmetric_limit(arrays: list[np.ndarray | None]) -> tuple[float, float]:
    vals = []
    for arr in arrays:
        if arr is None:
            continue
        a = np.asarray(arr, dtype=float).ravel()
        a = a[np.isfinite(a)]
        if a.size:
            vals.append(a)
    if not vals:
        return -1.0, 1.0
    vmax = float(np.nanmax(np.abs(np.concatenate(vals))))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    vmax *= 1.08
    return -vmax, vmax


def _infer_fit_target(df, fit_result: dict | None) -> str:
    target = str(fit_result.get("run", {}).get("fit_target", "")).strip().lower() if fit_result else ""
    if target in {"absolute", "delta"}:
        return target
    if "predicted_total_dx_arcsec" in df.columns or "predicted_total_dy_arcsec" in df.columns:
        return "delta"
    return "absolute"


def _infer_solve_mode(fit_result: dict | None) -> str:
    if not fit_result:
        return ""
    return str(fit_result.get("run", {}).get("solve_mode", "")).strip()


def _finite_array(df, col: str | None) -> np.ndarray | None:
    if col is None or col not in df.columns:
        return None
    return np.asarray(df[col], dtype=float)


def _pick_first_available(df, *cols: str) -> np.ndarray | None:
    for col in cols:
        arr = _finite_array(df, col)
        if arr is not None:
            return arr
    return None


def _used_mask(df, az: np.ndarray, el: np.ndarray) -> np.ndarray:
    mask = np.isfinite(az) & np.isfinite(el)
    if "is_outlier_block" in df.columns:
        mask &= ~df["is_outlier_block"].astype(bool).to_numpy()
    if "robust_weight" in df.columns:
        w = np.asarray(df["robust_weight"], dtype=float)
        mask &= np.isfinite(w) & (w > 0.0)
    return mask


def _outlier_mask(df, az: np.ndarray, el: np.ndarray) -> np.ndarray:
    mask = np.isfinite(az) & np.isfinite(el)
    if "is_outlier_block" in df.columns:
        return mask & df["is_outlier_block"].astype(bool).to_numpy()
    return np.zeros_like(mask, dtype=bool)


def _report_arrays(df):
    return {
        "used_dx": _finite_array(df, "used_model_dx_arcsec"),
        "used_dy": _finite_array(df, "used_model_dy_arcsec"),
        "delta_dx": _pick_first_available(df, "delta_err_dx_arcsec", "measured_dx_arcsec"),
        "delta_dy": _pick_first_available(df, "delta_err_dy_arcsec", "measured_dy_arcsec"),
        "target_dx": _finite_array(df, "target_dx_arcsec"),
        "target_dy": _finite_array(df, "target_dy_arcsec"),
        "pred_dx": _finite_array(df, "predicted_dx_arcsec"),
        "pred_dy": _finite_array(df, "predicted_dy_arcsec"),
        "res_dx": _finite_array(df, "residual_dx_arcsec"),
        "res_dy": _finite_array(df, "residual_dy_arcsec"),
        "upd_dx": _finite_array(df, "delta_to_used_dx_arcsec"),
        "upd_dy": _finite_array(df, "delta_to_used_dy_arcsec"),
        "pred_total_dx": _finite_array(df, "predicted_total_dx_arcsec"),
        "pred_total_dy": _finite_array(df, "predicted_total_dy_arcsec"),
        "radial_res": _finite_array(df, "radial_residual_arcsec"),
    }


def _profile_bin_count(n_valid: int) -> int:
    if n_valid <= 0:
        return 8
    return int(np.clip(np.round(np.sqrt(max(n_valid, 1))), 8, 18))


def _binned_profile(x: np.ndarray, y: np.ndarray | None, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if y is None:
        return np.array([], dtype=float), np.array([], dtype=float)
    valid = mask & np.isfinite(x) & np.isfinite(y)
    if not np.any(valid):
        return np.array([], dtype=float), np.array([], dtype=float)
    xv = np.asarray(x[valid], dtype=float)
    yv = np.asarray(y[valid], dtype=float)
    if xv.size <= 2:
        order = np.argsort(xv)
        return xv[order], yv[order]
    xmin = float(np.nanmin(xv))
    xmax = float(np.nanmax(xv))
    if not np.isfinite(xmin) or not np.isfinite(xmax):
        return np.array([], dtype=float), np.array([], dtype=float)
    if np.isclose(xmin, xmax):
        return np.array([float(np.nanmedian(xv))]), np.array([float(np.nanmedian(yv))])
    nbins = _profile_bin_count(xv.size)
    edges = np.linspace(xmin, xmax, nbins + 1)
    idx = np.clip(np.digitize(xv, edges) - 1, 0, nbins - 1)
    xout: list[float] = []
    yout: list[float] = []
    min_count = 1 if xv.size < 40 else 2
    for i in range(nbins):
        sel = idx == i
        if int(np.count_nonzero(sel)) < min_count:
            continue
        xout.append(float(np.nanmedian(xv[sel])))
        yout.append(float(np.nanmedian(yv[sel])))
    if not xout:
        order = np.argsort(xv)
        return xv[order], yv[order]
    return np.asarray(xout, dtype=float), np.asarray(yout, dtype=float)


def _style_series_axis(ax, xlabel: str, ylabel: str, title: str, ylim: tuple[float, float] | None = None):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True)
    ax.axhline(0.0, color="0.65", lw=0.9, zorder=0)
    if ylim is not None:
        ax.set_ylim(*ylim)


def _plot_scatter_points(ax, x: np.ndarray, y: np.ndarray | None, mask: np.ndarray, label: str, color: str, alpha: float = 0.25, s: float = 12.0, zorder: int = 1, marker: str = "o"):
    if y is None:
        return
    valid = mask & np.isfinite(x) & np.isfinite(y)
    if not np.any(valid):
        return
    ax.scatter(x[valid], y[valid], s=s, color=color, alpha=alpha, edgecolors="none", label=label, zorder=zorder, marker=marker)


def _plot_outliers(ax, x: np.ndarray, y: np.ndarray | None, outlier_mask: np.ndarray):
    if y is None:
        return
    valid = outlier_mask & np.isfinite(x) & np.isfinite(y)
    if not np.any(valid):
        return
    ax.scatter(x[valid], y[valid], s=28, color="tab:red", marker="x", linewidths=1.0, label="outlier", zorder=4)


def _plot_binned_curve(ax, x: np.ndarray, y: np.ndarray | None, mask: np.ndarray, label: str, color: str, linestyle: str = "-", marker: str = "o", lw: float = 2.0, alpha: float = 1.0, zorder: int = 3):
    xb, yb = _binned_profile(x, y, mask)
    if xb.size == 0:
        return
    ax.plot(xb, yb, linestyle=linestyle, marker=marker, color=color, lw=lw, ms=4.5, alpha=alpha, label=label, zorder=zorder)


def _legend_first_panel(ax, loc: str = "best"):
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        uniq_h = []
        uniq_l = []
        seen = set()
        for h, l in zip(handles, labels):
            if l in seen or not l:
                continue
            seen.add(l)
            uniq_h.append(h)
            uniq_l.append(l)
        ax.legend(uniq_h, uniq_l, loc=loc, fontsize=8)


def _panel_text(ax, lines: list[str]):
    ax.text(
        0.02,
        0.98,
        "\n".join(lines),
        ha="left",
        va="top",
        transform=ax.transAxes,
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.8"),
    )


def _plot_1d_main_panel(ax, x: np.ndarray, x_label: str, y_label: str, title: str, y_target: np.ndarray | None, y_fit: np.ndarray | None, used_mask: np.ndarray, outlier_mask: np.ndarray, ylim: tuple[float, float], y_used: np.ndarray | None = None, target_label: str = "target", fit_label: str = "fit", used_label: str = "used model"):
    _style_series_axis(ax, x_label, y_label, title, ylim)
    _plot_scatter_points(ax, x, y_target, used_mask, f"{target_label} points", "tab:blue", alpha=0.22, s=12.0, zorder=1)
    _plot_outliers(ax, x, y_target, outlier_mask)
    if y_used is not None:
        _plot_binned_curve(ax, x, y_used, used_mask, used_label, "0.45", linestyle="--", marker=None, lw=1.7, alpha=1.0, zorder=2)
    _plot_binned_curve(ax, x, y_target, used_mask, f"{target_label} median", "tab:blue", linestyle="-", marker="o", lw=1.8, alpha=0.95, zorder=3)
    _plot_binned_curve(ax, x, y_fit, used_mask, f"{fit_label} median", "tab:red", linestyle="-", marker="o", lw=2.2, alpha=0.95, zorder=4)


def _plot_1d_compare_panel(ax, x: np.ndarray, x_label: str, y_label: str, title: str, y_obs: np.ndarray | None, y_model: np.ndarray | None, used_mask: np.ndarray, outlier_mask: np.ndarray, ylim: tuple[float, float], obs_label: str, model_label: str, y_ref: np.ndarray | None = None, ref_label: str | None = None):
    _style_series_axis(ax, x_label, y_label, title, ylim)
    _plot_scatter_points(ax, x, y_obs, used_mask, obs_label + " points", "tab:blue", alpha=0.22, s=12.0, zorder=1)
    _plot_outliers(ax, x, y_obs, outlier_mask)
    _plot_binned_curve(ax, x, y_obs, used_mask, obs_label + " median", "tab:blue", linestyle="-", marker="o", lw=1.9, alpha=0.95, zorder=3)
    _plot_binned_curve(ax, x, y_model, used_mask, model_label + " median", "tab:red", linestyle="-", marker="o", lw=2.2, alpha=0.95, zorder=4)
    if y_ref is not None and ref_label:
        _plot_binned_curve(ax, x, y_ref, used_mask, ref_label, "0.45", linestyle="--", marker=None, lw=1.6, alpha=1.0, zorder=2)


def _plot_1d_residual_panel(ax, x: np.ndarray, x_label: str, y_label: str, title: str, y_res: np.ndarray | None, used_mask: np.ndarray, outlier_mask: np.ndarray, ylim: tuple[float, float]):
    _style_series_axis(ax, x_label, y_label, title, ylim)
    _plot_scatter_points(ax, x, y_res, used_mask, "residual points", "tab:blue", alpha=0.25, s=12.0, zorder=1)
    _plot_outliers(ax, x, y_res, outlier_mask)
    _plot_binned_curve(ax, x, y_res, used_mask, "residual median", "tab:red", linestyle="-", marker="o", lw=2.1, alpha=0.95, zorder=3)


def _plot_coverage_panel(ax, df, az: np.ndarray, el: np.ndarray, fit_target: str, solve_mode: str):
    used = _used_mask(df, az, el)
    outlier = _outlier_mask(df, az, el)
    base = np.isfinite(az) & np.isfinite(el)
    dataset_col = "dataset_id" if "dataset_id" in df.columns else None
    ax.set_xlabel("Az [deg]")
    ax.set_ylabel("El [deg]")
    ax.set_title("Coverage / sampling")
    ax.grid(True)
    if dataset_col is not None:
        dataset = df[dataset_col].astype(str).to_numpy()
        uniq = [u for u in np.unique(dataset[base]) if str(u) != "nan"]
        cmap = plt.get_cmap("tab10")
        if 0 < len(uniq) <= 8:
            for idx, name in enumerate(uniq):
                trace = base & (dataset == name)
                if np.any(trace):
                    ax.plot(az[trace], el[trace], color=cmap(idx), lw=0.8, alpha=0.25, zorder=1)
                sel = used & (dataset == name)
                if np.any(sel):
                    ax.scatter(az[sel], el[sel], s=18, color=cmap(idx), alpha=0.9, label=name, zorder=3)
        else:
            if np.any(base):
                ax.plot(az[base], el[base], color="0.55", lw=0.8, alpha=0.30, zorder=1)
            ax.scatter(az[used], el[used], s=18, color="tab:blue", alpha=0.90, label="used", zorder=3)
    else:
        if np.any(base):
            ax.plot(az[base], el[base], color="0.55", lw=0.8, alpha=0.30, zorder=1)
        ax.scatter(az[used], el[used], s=18, color="tab:blue", alpha=0.90, label="used", zorder=3)
    if np.any(outlier):
        ax.scatter(az[outlier], el[outlier], s=32, color="tab:red", marker="x", linewidths=1.1, label="outlier", zorder=4)
    unused = base & ~used & ~outlier
    if np.any(unused):
        ax.scatter(az[unused], el[unused], s=20, facecolors="none", edgecolors="0.55", marker="o", linewidths=0.9, label="unused", zorder=2)
    n_total = int(np.count_nonzero(base))
    n_used = int(np.count_nonzero(used))
    n_out = int(np.count_nonzero(outlier))
    lines = [
        f"fit_target = {fit_target}",
        f"solve_mode = {solve_mode or 'unknown'}",
        f"rows = {n_total}",
        f"used = {n_used}",
    ]
    if n_out > 0:
        lines.append(f"outlier = {n_out}")
    _panel_text(ax, lines)
    _legend_first_panel(ax)


def _main_limits_absolute(arrs) -> tuple[tuple[float, float], tuple[float, float]]:
    dx_lim = _nan_limits(*(a for a in [arrs["used_dx"], arrs["target_dx"], arrs["pred_dx"]] if a is not None))
    dy_lim = _nan_limits(*(a for a in [arrs["used_dy"], arrs["target_dy"], arrs["pred_dy"]] if a is not None))
    return dx_lim, dy_lim


def _main_limits_delta(arrs) -> tuple[tuple[float, float], tuple[float, float]]:
    dx_lim = _nan_limits(*(a for a in [arrs["delta_dx"], arrs["pred_dx"]] if a is not None))
    dy_lim = _nan_limits(*(a for a in [arrs["delta_dy"], arrs["pred_dy"]] if a is not None))
    return dx_lim, dy_lim


def _save_overview_figure(outdir: Path, title: str, df, fit_target: str, solve_mode: str, az, el, arrs) -> Path:
    used_mask = _used_mask(df, az, el)
    outlier_mask = _outlier_mask(df, az, el)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    if fit_target == "absolute":
        dx_lim, dy_lim = _main_limits_absolute(arrs)
        _plot_1d_main_panel(axes[0, 0], az, "Az [deg]", "dx [arcsec]", "dx vs Az", arrs["target_dx"], arrs["pred_dx"], used_mask, outlier_mask, dx_lim, y_used=arrs["used_dx"], target_label="pointing absolute target", fit_label="fitted pointing model", used_label="observation-time pointing model")
        _plot_1d_main_panel(axes[0, 1], az, "Az [deg]", "dy [arcsec]", "dy vs Az", arrs["target_dy"], arrs["pred_dy"], used_mask, outlier_mask, dy_lim, y_used=arrs["used_dy"], target_label="pointing absolute target", fit_label="fitted pointing model", used_label="observation-time pointing model")
        _plot_1d_main_panel(axes[1, 0], el, "El [deg]", "dx [arcsec]", "dx vs El", arrs["target_dx"], arrs["pred_dx"], used_mask, outlier_mask, dx_lim, y_used=arrs["used_dx"], target_label="pointing absolute target", fit_label="fitted pointing model", used_label="observation-time pointing model")
        _plot_1d_main_panel(axes[1, 1], el, "El [deg]", "dy [arcsec]", "dy vs El", arrs["target_dy"], arrs["pred_dy"], used_mask, outlier_mask, dy_lim, y_used=arrs["used_dy"], target_label="pointing absolute target", fit_label="fitted pointing model", used_label="observation-time pointing model")
        fig.suptitle(f"{title} : 1D main diagnostics (pointing absolute fit)", fontsize=14)
        _panel_text(axes[0, 0], [f"fit_target = {fit_target}", f"solve_mode = {solve_mode or 'unknown'}", "gray dashed = observation-time pointing model", "blue = pointing absolute target (= observation-time model + observed correction)", "red = fitted pointing model"]) 
    else:
        dx_lim, dy_lim = _main_limits_delta(arrs)
        _plot_1d_main_panel(axes[0, 0], az, "Az [deg]", "dx [arcsec]", "dx vs Az", arrs["delta_dx"], arrs["pred_dx"], used_mask, outlier_mask, dx_lim, y_used=None, target_label="observed correction", fit_label="fitted correction model")
        _plot_1d_main_panel(axes[0, 1], az, "Az [deg]", "dy [arcsec]", "dy vs Az", arrs["delta_dy"], arrs["pred_dy"], used_mask, outlier_mask, dy_lim, y_used=None, target_label="observed correction", fit_label="fitted correction model")
        _plot_1d_main_panel(axes[1, 0], el, "El [deg]", "dx [arcsec]", "dx vs El", arrs["delta_dx"], arrs["pred_dx"], used_mask, outlier_mask, dx_lim, y_used=None, target_label="observed correction", fit_label="fitted correction model")
        _plot_1d_main_panel(axes[1, 1], el, "El [deg]", "dy [arcsec]", "dy vs El", arrs["delta_dy"], arrs["pred_dy"], used_mask, outlier_mask, dy_lim, y_used=None, target_label="observed correction", fit_label="fitted correction model")
        fig.suptitle(f"{title} : 1D main diagnostics (delta fit)", fontsize=14)
        _panel_text(axes[0, 0], [f"fit_target = {fit_target}", f"solve_mode = {solve_mode or 'unknown'}", "blue = observed correction relative to observation-time model", "red = fitted correction model"]) 
    _legend_first_panel(axes[0, 0], loc="best")
    out = outdir / "fit_diagnostics_azel_map.png"
    fig.savefig(out, dpi=150)
    compat = locals().get("compat_out")
    if compat is not None and compat != out:
        fig.savefig(compat, dpi=150)
    plt.close(fig)
    return out


def _save_component_figure(outdir: Path, title: str, fit_target: str, df, az, el, arrs) -> Path:
    used_mask = _used_mask(df, az, el)
    outlier_mask = _outlier_mask(df, az, el)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    res_dx_lim, res_dy_lim = _symmetric_limit([arrs["res_dx"]]), _symmetric_limit([arrs["res_dy"]])
    _plot_1d_residual_panel(axes[0, 0], az, "Az [deg]", "residual dx [arcsec]", "Residual dx vs Az", arrs["res_dx"], used_mask, outlier_mask, res_dx_lim)
    _plot_1d_residual_panel(axes[0, 1], az, "Az [deg]", "residual dy [arcsec]", "Residual dy vs Az", arrs["res_dy"], used_mask, outlier_mask, res_dy_lim)
    _plot_1d_residual_panel(axes[1, 0], el, "El [deg]", "residual dx [arcsec]", "Residual dx vs El", arrs["res_dx"], used_mask, outlier_mask, res_dx_lim)
    _plot_1d_residual_panel(axes[1, 1], el, "El [deg]", "residual dy [arcsec]", "Residual dy vs El", arrs["res_dy"], used_mask, outlier_mask, res_dy_lim)
    fig.suptitle(f"{title} : 1D residual diagnostics ({fit_target})", fontsize=14)
    _legend_first_panel(axes[0, 0], loc="best")
    out = outdir / "fit_diagnostics_azel_components.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _save_update_or_total_figure(outdir: Path, title: str, fit_target: str, df, az, el, arrs) -> Path | None:
    used_mask = _used_mask(df, az, el)
    outlier_mask = _outlier_mask(df, az, el)
    if fit_target == "absolute":
        if arrs["upd_dx"] is None or arrs["upd_dy"] is None:
            return None
        fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
        dx_lim = _nan_limits(*(a for a in [arrs["delta_dx"], arrs["upd_dx"]] if a is not None))
        dy_lim = _nan_limits(*(a for a in [arrs["delta_dy"], arrs["upd_dy"]] if a is not None))
        _plot_1d_compare_panel(axes[0, 0], az, "Az [deg]", "dx [arcsec]", "Observed correction vs fitted correction: dx vs Az", arrs["delta_dx"], arrs["upd_dx"], used_mask, outlier_mask, dx_lim, "observed correction", "fitted correction")
        _plot_1d_compare_panel(axes[0, 1], az, "Az [deg]", "dy [arcsec]", "Observed correction vs fitted correction: dy vs Az", arrs["delta_dy"], arrs["upd_dy"], used_mask, outlier_mask, dy_lim, "observed correction", "fitted correction")
        _plot_1d_compare_panel(axes[1, 0], el, "El [deg]", "dx [arcsec]", "Observed correction vs fitted correction: dx vs El", arrs["delta_dx"], arrs["upd_dx"], used_mask, outlier_mask, dx_lim, "observed correction", "fitted correction")
        _plot_1d_compare_panel(axes[1, 1], el, "El [deg]", "dy [arcsec]", "Observed correction vs fitted correction: dy vs El", arrs["delta_dy"], arrs["upd_dy"], used_mask, outlier_mask, dy_lim, "observed correction", "fitted correction")
        fig.suptitle(f"{title} : observed correction vs fitted correction relative to observation-time pointing model", fontsize=14)
        _legend_first_panel(axes[0, 0], loc="best")
        out = outdir / "fit_diagnostics_azel_correction_relative_to_used.png"
        compat_out = outdir / "fit_diagnostics_azel_update.png"
    else:
        if arrs["pred_total_dx"] is None or arrs["pred_total_dy"] is None:
            return None
        fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
        dx_lim = _nan_limits(*(a for a in [arrs["used_dx"], arrs["pred_total_dx"]] if a is not None))
        dy_lim = _nan_limits(*(a for a in [arrs["used_dy"], arrs["pred_total_dy"]] if a is not None))
        _plot_1d_compare_panel(axes[0, 0], az, "Az [deg]", "dx [arcsec]", "Observation-time pointing model vs total model after delta fit: dx vs Az", arrs["pred_total_dx"], arrs["pred_total_dx"], used_mask, outlier_mask, dx_lim, "predicted total", "predicted total", y_ref=arrs["used_dx"], ref_label="observation-time pointing model")
        _plot_1d_compare_panel(axes[0, 1], az, "Az [deg]", "dy [arcsec]", "Observation-time pointing model vs total model after delta fit: dy vs Az", arrs["pred_total_dy"], arrs["pred_total_dy"], used_mask, outlier_mask, dy_lim, "predicted total", "predicted total", y_ref=arrs["used_dy"], ref_label="observation-time pointing model")
        _plot_1d_compare_panel(axes[1, 0], el, "El [deg]", "dx [arcsec]", "Observation-time pointing model vs total model after delta fit: dx vs El", arrs["pred_total_dx"], arrs["pred_total_dx"], used_mask, outlier_mask, dx_lim, "predicted total", "predicted total", y_ref=arrs["used_dx"], ref_label="observation-time pointing model")
        _plot_1d_compare_panel(axes[1, 1], el, "El [deg]", "dy [arcsec]", "Observation-time pointing model vs total model after delta fit: dy vs El", arrs["pred_total_dy"], arrs["pred_total_dy"], used_mask, outlier_mask, dy_lim, "predicted total", "predicted total", y_ref=arrs["used_dy"], ref_label="observation-time pointing model")
        fig.suptitle(f"{title} : total pointing model after delta fit", fontsize=14)
        _legend_first_panel(axes[0, 0], loc="best")
        out = outdir / "fit_diagnostics_azel_total.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _prepare_3d_arrays(az: np.ndarray, el: np.ndarray, z_obs: np.ndarray | None, z_fit: np.ndarray | None, used_mask: np.ndarray):
    if z_obs is None or z_fit is None:
        return None
    valid = used_mask & np.isfinite(az) & np.isfinite(el) & np.isfinite(z_obs) & np.isfinite(z_fit)
    if not np.any(valid):
        return None
    azv = np.asarray(az[valid], dtype=float)
    elv = np.asarray(el[valid], dtype=float)
    zobsv = np.asarray(z_obs[valid], dtype=float)
    zfitv = np.asarray(z_fit[valid], dtype=float)
    return azv, elv, zobsv, zfitv


def _set_3d_common(ax, title: str, zlabel: str, elev: float, azim: float):
    ax.set_title(title)
    ax.set_xlabel("Az [deg]")
    ax.set_ylabel("El [deg]")
    ax.set_zlabel(zlabel)
    ax.view_init(elev=elev, azim=azim)


def _plot_3d_surface_overview_panel(ax, az: np.ndarray, el: np.ndarray, z_obs: np.ndarray | None, z_fit: np.ndarray | None, used_mask: np.ndarray, title: str, zlabel: str, elev: float, azim: float):
    _set_3d_common(ax, title, zlabel, elev, azim)
    prepared = _prepare_3d_arrays(az, el, z_obs, z_fit, used_mask)
    if prepared is None:
        ax.text2D(0.5, 0.5, "not available", transform=ax.transAxes, ha="center", va="center")
        return
    azv, elv, zobsv, zfitv = prepared
    plotted_surface = False
    if azv.size >= 3:
        try:
            surf = ax.plot_trisurf(
                azv,
                elv,
                zfitv,
                cmap="viridis",
                alpha=0.70,
                linewidth=0.55,
                edgecolor="0.20",
                antialiased=True,
                shade=True,
            )
            surf.set_clim(float(np.nanmin(zfitv)), float(np.nanmax(zfitv)))
            plotted_surface = True
        except Exception:
            plotted_surface = False
    if not plotted_surface:
        ax.scatter(azv, elv, zfitv, s=12, alpha=0.50, color="tab:green", depthshade=True)
    ax.scatter(azv, elv, zobsv, s=16, color="k", alpha=0.42, edgecolors="none", depthshade=False)
    ax.text2D(
        0.02,
        0.98,
        "surface = fitted model\nblack points = observed values\n(no signed-distance coloring in this figure)",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.82, edgecolor="0.8"),
    )

def _plot_3d_signed_distance_panel(ax, az: np.ndarray, el: np.ndarray, z_obs: np.ndarray | None, z_fit: np.ndarray | None, used_mask: np.ndarray, title: str, zlabel: str, elev: float, azim: float):
    _set_3d_common(ax, title, zlabel, elev, azim)
    prepared = _prepare_3d_arrays(az, el, z_obs, z_fit, used_mask)
    if prepared is None:
        ax.text2D(0.5, 0.5, "not available", transform=ax.transAxes, ha="center", va="center")
        return
    azv, elv, zobsv, zfitv = prepared
    diff = zobsv - zfitv
    max_abs = float(np.nanmax(np.abs(diff))) if diff.size else 1.0
    if not np.isfinite(max_abs) or max_abs <= 0.0:
        max_abs = 1.0
    norm = Normalize(vmin=-max_abs, vmax=max_abs)
    cmap = plt.get_cmap("coolwarm")
    plotted_surface = False
    if azv.size >= 3:
        try:
            ax.plot_trisurf(
                azv,
                elv,
                zfitv,
                color="0.70",
                alpha=0.18,
                linewidth=0.12,
                edgecolor="0.70",
                antialiased=True,
                shade=False,
            )
            plotted_surface = True
        except Exception:
            plotted_surface = False
    if not plotted_surface:
        ax.scatter(azv, elv, zfitv, s=10, alpha=0.18, color="0.65", depthshade=False)
    for xi, yi, zoi, zfi, di in zip(azv, elv, zobsv, zfitv, diff):
        ax.plot([xi, xi], [yi, yi], [zfi, zoi], color=cmap(norm(di)), alpha=0.90, lw=1.8)
    ax.scatter(azv, elv, zobsv, s=34, c=diff, cmap=cmap, norm=norm, alpha=0.98, edgecolors="k", linewidths=0.35, depthshade=False)
    ax.text2D(
        0.02,
        0.98,
        "gray surface = fitted model\nred point/segment: observed value above fitted surface\nblue point/segment: observed value below fitted surface",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.82, edgecolor="0.8"),
    )

def _save_3d_figure(outdir: Path, title: str, fit_target: str, df, az, el, arrs) -> Path | None:
    used_mask = _used_mask(df, az, el)
    if fit_target == "absolute":
        obs_dx, obs_dy = arrs["target_dx"], arrs["target_dy"]
        fit_dx, fit_dy = arrs["pred_dx"], arrs["pred_dy"]
        obs_label = "pointing absolute target"
        fit_label = "fitted pointing model"
    else:
        obs_dx, obs_dy = arrs["delta_dx"], arrs["delta_dy"]
        fit_dx, fit_dy = arrs["pred_dx"], arrs["pred_dy"]
        obs_label = "observed correction relative to observation-time model"
        fit_label = "fitted correction model"
    if obs_dx is None or obs_dy is None or fit_dx is None or fit_dy is None:
        return None
    target_name = "pointing absolute fit" if fit_target == "absolute" else "delta fit"

    fig = plt.figure(figsize=(18, 12), constrained_layout=True)
    axes = [fig.add_subplot(2, 2, i, projection="3d") for i in range(1, 5)]
    _plot_3d_surface_overview_panel(axes[0], az, el, obs_dx, fit_dx, used_mask, f"dx 3D surface overview: {obs_label} vs {fit_label} (view A)", "dx [arcsec]", elev=26.0, azim=-58.0)
    _plot_3d_surface_overview_panel(axes[1], az, el, obs_dy, fit_dy, used_mask, f"dy 3D surface overview: {obs_label} vs {fit_label} (view A)", "dy [arcsec]", elev=26.0, azim=-58.0)
    _plot_3d_surface_overview_panel(axes[2], az, el, obs_dx, fit_dx, used_mask, f"dx 3D surface overview: {obs_label} vs {fit_label} (view B)", "dx [arcsec]", elev=26.0, azim=32.0)
    _plot_3d_surface_overview_panel(axes[3], az, el, obs_dy, fit_dy, used_mask, f"dy 3D surface overview: {obs_label} vs {fit_label} (view B)", "dy [arcsec]", elev=26.0, azim=32.0)
    fig.suptitle(f"{title} : supplementary 3D surface overview ({target_name})", fontsize=14)
    out_surface = outdir / "fit_diagnostics_3d.png"
    fig.savefig(out_surface, dpi=150)
    plt.close(fig)

    fig = plt.figure(figsize=(18, 12), constrained_layout=True)
    axes = [fig.add_subplot(2, 2, i, projection="3d") for i in range(1, 5)]
    _plot_3d_signed_distance_panel(axes[0], az, el, obs_dx, fit_dx, used_mask, f"dx 3D signed distance: {obs_label} vs {fit_label} (view A)", "dx [arcsec]", elev=26.0, azim=-58.0)
    _plot_3d_signed_distance_panel(axes[1], az, el, obs_dy, fit_dy, used_mask, f"dy 3D signed distance: {obs_label} vs {fit_label} (view A)", "dy [arcsec]", elev=26.0, azim=-58.0)
    _plot_3d_signed_distance_panel(axes[2], az, el, obs_dx, fit_dx, used_mask, f"dx 3D signed distance: {obs_label} vs {fit_label} (view B)", "dx [arcsec]", elev=26.0, azim=32.0)
    _plot_3d_signed_distance_panel(axes[3], az, el, obs_dy, fit_dy, used_mask, f"dy 3D signed distance: {obs_label} vs {fit_label} (view B)", "dy [arcsec]", elev=26.0, azim=32.0)
    fig.suptitle(f"{title} : supplementary 3D signed-distance diagnostics ({target_name})", fontsize=14)
    out_signed = outdir / "fit_diagnostics_3d_signed_distance.png"
    fig.savefig(out_signed, dpi=150)
    plt.close(fig)
    return out_signed

def run_report(opts: ReportOptions) -> dict[str, str]:
    outdir = Path(opts.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = Path(opts.input_csv)
    df = read_csv(csv_path)
    fit_result = _load_fit_result_near(csv_path)

    fig = plt.figure(figsize=(14, 8))
    ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]

    ax[0].plot(df["az_deg"], df["el_deg"], ".")
    ax[0].set_xlabel("Az [deg]")
    ax[0].set_ylabel("El [deg]")

    ax[1].plot(df["az_deg"], df["residual_dx_arcsec"], ".")
    ax[1].set_xlabel("Az [deg]")
    ax[1].set_ylabel("Residual dx [arcsec]")

    ax[2].plot(df["el_deg"], df["residual_dx_arcsec"], ".")
    ax[2].set_xlabel("El [deg]")
    ax[2].set_ylabel("Residual dx [arcsec]")

    ax[3].plot(df["residual_dx_arcsec"], df["residual_dy_arcsec"], ".")
    th = np.linspace(0.0, 2.0 * np.pi, 360)
    r = 10.0
    ax[3].plot(r * np.cos(th), r * np.sin(th), color="tab:orange", lw=1.4)
    ax[3].set_xlabel("Residual dx [arcsec]")
    ax[3].set_ylabel("Residual dy [arcsec]")
    lim = _symmetric_limit([np.asarray(df["residual_dx_arcsec"], dtype=float), np.asarray(df["residual_dy_arcsec"], dtype=float), np.array([r, -r], dtype=float)])
    ax[3].set_xlim(*lim)
    ax[3].set_ylim(*lim)
    ax[3].set_aspect("equal", adjustable="box")

    ax[4].plot(df["az_deg"], df["residual_dy_arcsec"], ".")
    ax[4].set_xlabel("Az [deg]")
    ax[4].set_ylabel("Residual dy [arcsec]")

    ax[5].plot(df["el_deg"], df["residual_dy_arcsec"], ".")
    ax[5].set_xlabel("El [deg]")
    ax[5].set_ylabel("Residual dy [arcsec]")

    for a in ax:
        a.grid(True)
    fig.suptitle(opts.title)
    png = outdir / "fit_diagnostics.png"
    fig.tight_layout()
    fig.savefig(png, dpi=150)
    plt.close(fig)

    az = np.asarray(df["az_deg"], dtype=float)
    el = np.asarray(df["el_deg"], dtype=float)
    fit_target = _infer_fit_target(df, fit_result)
    solve_mode = _infer_solve_mode(fit_result)
    arrs = _report_arrays(df)
    fit_map_png = _save_overview_figure(outdir, opts.title, df, fit_target, solve_mode, az, el, arrs)
    _save_update_or_total_figure(outdir, opts.title, fit_target, df, az, el, arrs)
    _save_3d_figure(outdir, opts.title, fit_target, df, az, el, arrs)

    md = outdir / "report_summary.md"
    radial_rms = float(np.sqrt(np.nanmean(df["radial_residual_arcsec"].to_numpy(dtype=float) ** 2)))
    text = [
        f"# {opts.title}",
        "",
        f"- rows: {len(df)}",
        f"- radial RMS [arcsec]: {radial_rms:.3f}",
        f"- diagnostics: `{png.name}`",
        f"- Az/El fit maps: `{fit_map_png.name}`",
    ]
    if fit_result:
        text.extend([
            f"- solve mode: `{fit_result.get('run', {}).get('solve_mode', '')}`",
            f"- model: `{fit_result.get('run', {}).get('model_name', '')}`",
        ])
    text.extend(_derived_lines(fit_result))
    md.write_text("\n".join(text) + "\n", encoding="utf-8")
    return {"png": str(png), "fit_map_png": str(fit_map_png), "markdown": str(md)}
