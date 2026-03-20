from __future__ import annotations
import argparse
import builtins
import math
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REQUIRED_COLUMNS = ["x_rel_arcsec", "y_rel_arcsec", "pred_x_arcsec", "pred_y_arcsec"]
OPTIONAL_NUMERIC_COLUMNS = [
    "rep_el_deg",
    "rotation_sign",
    "reference_angle_deg",
    "dewar_angle_deg",
]


def _rotate_xy(x_arcsec: float, y_arcsec: float, theta_deg: float) -> Tuple[float, float]:
    th = math.radians(float(theta_deg))
    c = math.cos(th)
    s = math.sin(th)
    return float(c * x_arcsec - s * y_arcsec), float(s * x_arcsec + c * y_arcsec)


def _inverse_rotate_xy(x_arcsec: float, y_arcsec: float, theta_deg: float) -> Tuple[float, float]:
    return _rotate_xy(x_arcsec, y_arcsec, -float(theta_deg))


def _coerce_numeric_columns(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def _theta_from_row(row: pd.Series, rep_el_deg: float) -> Optional[float]:
    if not np.isfinite(rep_el_deg):
        return None
    sign = float(row["rotation_sign"]) if "rotation_sign" in row and np.isfinite(row["rotation_sign"]) else 0.0
    ref = float(row["reference_angle_deg"]) if "reference_angle_deg" in row and np.isfinite(row["reference_angle_deg"]) else 0.0
    dew = float(row["dewar_angle_deg"]) if "dewar_angle_deg" in row and np.isfinite(row["dewar_angle_deg"]) else 0.0
    return float(sign * (rep_el_deg - ref) + dew)


def _group_column(df: pd.DataFrame) -> str:
    if "beam_id" in df.columns:
        return "beam_id"
    if "stream_name" in df.columns:
        return "stream_name"
    return "__group__"


def _estimate_template_by_group(df: pd.DataFrame, group_col: str) -> Dict[str, Tuple[float, float]]:
    templates: Dict[str, Tuple[float, float]] = {}
    needed = {"rep_el_deg", "rotation_sign", "pred_x_arcsec", "pred_y_arcsec"}
    if not needed.issubset(df.columns):
        return templates
    for group_name, sub in df.groupby(group_col, sort=True):
        vals = []
        for _, row in sub.iterrows():
            rep_el = float(row["rep_el_deg"])
            theta = _theta_from_row(row, rep_el)
            if theta is None:
                continue
            x0, y0 = _inverse_rotate_xy(float(row["pred_x_arcsec"]), float(row["pred_y_arcsec"]), theta)
            vals.append((x0, y0))
        if vals:
            arr = np.asarray(vals, dtype=float)
            templates[str(group_name)] = (float(np.nanmean(arr[:, 0])), float(np.nanmean(arr[:, 1])))
    return templates


def _label_text(group_name: str, sub: pd.DataFrame) -> str:
    lines = [str(group_name)]
    if "stream_name" in sub.columns:
        streams = sorted({str(v) for v in sub["stream_name"].dropna().astype(str).tolist()})
        if len(streams) == 1 and streams[0] != str(group_name):
            lines.append(streams[0])
    if "rep_el_deg" in sub.columns:
        vals = pd.to_numeric(sub["rep_el_deg"], errors="coerce")
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            lines.append(f"El {float(vals.min()):.0f}°→{float(vals.max()):.0f}°")
    return "\n".join(lines)


def _choose_label_anchor(x_ref: float, y_ref: float, x_rel: np.ndarray, y_rel: np.ndarray) -> Tuple[float, float]:
    if np.isfinite(x_ref) and np.isfinite(y_ref):
        return float(x_ref), float(y_ref)
    r = np.hypot(x_rel, y_rel)
    if len(r) == 0:
        return 0.0, 0.0
    idx = int(np.nanargmax(r))
    return float(x_rel[idx]), float(y_rel[idx])


def _choose_label_position(anchor_x: float, anchor_y: float, span: float, *, radial_scale: float = 0.080, tangential_scale: float = 0.040, tangential_sign: Optional[float] = None) -> Tuple[float, float, str, str]:
    r = math.hypot(anchor_x, anchor_y)
    base_offset = radial_scale * span
    tangential_offset = tangential_scale * span
    if r < 1e-9:
        ux, uy = 1.0, 0.0
    else:
        ux, uy = anchor_x / r, anchor_y / r
    tx, ty = -uy, ux
    if tangential_sign is None:
        tangential_sign = 1.0 if anchor_y >= 0 else -1.0
    lx = anchor_x + base_offset * ux + tangential_sign * tangential_offset * tx
    ly = anchor_y + base_offset * uy + tangential_sign * tangential_offset * ty
    ha = "left" if ux >= 0 else "right"
    va = "bottom" if uy >= 0 else "top"
    return lx, ly, ha, va


def _friendly_title(title: Optional[str], output_path: Path) -> str:
    if title is None:
        stem = output_path.stem
        if stem.startswith("beam_fit_residuals_") and stem.endswith("_xy"):
            model_name = stem[len("beam_fit_residuals_"):-len("_xy")]
            return f"{model_name}: measured beam tracks and fitted model"
        return stem
    fixed = title.replace("rel vs pred", "measured beam tracks and fitted model")
    fixed = fixed.replace("Rel vs Pred", "measured beam tracks and fitted model")
    return fixed


def write_beam_fit_xy_png(
    residual_df: pd.DataFrame,
    output_path: Path,
    *,
    title: Optional[str] = None,
    connect_mode: str = "line",
    dpi: int = 150,
) -> Path:
    missing = [c for c in REQUIRED_COLUMNS if c not in residual_df.columns]
    if missing:
        raise ValueError(f"missing required columns for residual plot: {missing}")

    df = residual_df.copy()
    if "__group__" not in df.columns:
        df["__group__"] = "all"

    df = _coerce_numeric_columns(df, list(REQUIRED_COLUMNS) + OPTIONAL_NUMERIC_COLUMNS)
    group_col = _group_column(df)

    mask = np.ones(len(df), dtype=bool)
    for col in REQUIRED_COLUMNS:
        mask &= np.isfinite(df[col].to_numpy(dtype=float))
    df = df.loc[mask].copy()
    if df.empty:
        raise ValueError("no finite rows are available for residual plotting")

    if connect_mode not in {"none", "line", "arrow"}:
        raise ValueError(f"unsupported connect_mode={connect_mode!r}")

    templates = _estimate_template_by_group(df, group_col)

    fig, ax = plt.subplots(figsize=(9.2, 9.2))

    color_cycle = plt.rcParams.get("axes.prop_cycle", None)
    if color_cycle is not None:
        colors = color_cycle.by_key().get("color", ["C0", "C1", "C2", "C3", "C4"])
    else:
        colors = ["C0", "C1", "C2", "C3", "C4"]

    if connect_mode in {"line", "arrow"}:
        for _, row in df.iterrows():
            x0 = float(row["pred_x_arcsec"])
            y0 = float(row["pred_y_arcsec"])
            x1 = float(row["x_rel_arcsec"])
            y1 = float(row["y_rel_arcsec"])
            if connect_mode == "line":
                ax.plot([x0, x1], [y0, y1], "-", color="0.84", alpha=0.26, linewidth=0.70, zorder=1)
            else:
                ax.annotate(
                    "",
                    xy=(x1, y1),
                    xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", lw=0.65, color="0.84", alpha=0.26),
                    zorder=1,
                )

    x_all = np.concatenate([df["x_rel_arcsec"].to_numpy(float), df["pred_x_arcsec"].to_numpy(float)])
    y_all = np.concatenate([df["y_rel_arcsec"].to_numpy(float), df["pred_y_arcsec"].to_numpy(float)])

    if templates and {"rep_el_deg", "rotation_sign"}.issubset(df.columns):
        for group_name, sub in df.groupby(group_col, sort=True):
            key = str(group_name)
            if key not in templates:
                continue
            rep = sub.iloc[0]
            xs90 = []
            ys90 = []
            for el in np.linspace(0.0, 90.0, 181):
                theta = _theta_from_row(rep, float(el))
                if theta is None:
                    continue
                tx, ty = _rotate_xy(templates[key][0], templates[key][1], theta)
                xs90.append(tx)
                ys90.append(ty)
            if xs90:
                x_all = np.append(x_all, np.asarray(xs90, dtype=float))
                y_all = np.append(y_all, np.asarray(ys90, dtype=float))

    xmin = builtins.min(x_all)
    xmax = builtins.max(x_all)
    ymin = builtins.min(y_all)
    ymax = builtins.max(y_all)
    dx = xmax - xmin
    dy = ymax - ymin
    span = builtins.max(dx, dy)
    if not np.isfinite(span) or span <= 0:
        span = 1.0

    margin = 0.18 * span
    xmid = 0.5 * (xmin + xmax)
    ymid = 0.5 * (ymin + ymax)
    ax.set_xlim(xmid - 0.5 * span - margin, xmid + 0.5 * span + margin)
    ax.set_ylim(ymid - 0.5 * span - margin, ymid + 0.5 * span + margin)

    dummy_zero90 = None

    for ig, (group_name, sub) in enumerate(df.groupby(group_col, sort=True)):
        color = colors[ig % len(colors)]
        sub = sub.copy()
        if "rep_el_deg" in sub.columns:
            sub = sub.sort_values("rep_el_deg", kind="mergesort")
        else:
            sub = sub.sort_index()

        x_rel = sub["x_rel_arcsec"].to_numpy(dtype=float)
        y_rel = sub["y_rel_arcsec"].to_numpy(dtype=float)
        x_pred = sub["pred_x_arcsec"].to_numpy(dtype=float)
        y_pred = sub["pred_y_arcsec"].to_numpy(dtype=float)

        key = str(group_name)
        x_ref = np.nan
        y_ref = np.nan

        if key in templates and "rep_el_deg" in sub.columns and "rotation_sign" in sub.columns:
            rep = sub.iloc[0]
            xs90 = []
            ys90 = []
            for el in np.linspace(0.0, 90.0, 181):
                theta = _theta_from_row(rep, float(el))
                if theta is None:
                    continue
                tx, ty = _rotate_xy(templates[key][0], templates[key][1], theta)
                xs90.append(tx)
                ys90.append(ty)
            if xs90:
                line90, = ax.plot(xs90, ys90, linestyle=":", linewidth=1.1, color=color, alpha=0.18, zorder=1)
                if dummy_zero90 is None:
                    dummy_zero90 = line90

            theta0 = _theta_from_row(rep, 0.0)
            if theta0 is not None:
                x_ref, y_ref = _rotate_xy(templates[key][0], templates[key][1], theta0)
                ax.scatter(
                    [x_ref],
                    [y_ref],
                    s=130.0,
                    marker="*",
                    facecolors=color,
                    edgecolors="black",
                    linewidths=0.9,
                    alpha=0.98,
                    zorder=6,
                )

        ax.plot(x_pred, y_pred, linestyle="--", linewidth=1.0, color=color, alpha=0.58, zorder=2)
        ax.plot(x_rel, y_rel, linestyle="-", linewidth=1.7, color=color, alpha=0.95, zorder=3)

        ax.scatter(
            x_pred,
            y_pred,
            s=28.0,
            marker="x",
            color=color,
            alpha=0.76,
            linewidths=0.95,
            zorder=4,
        )
        ax.scatter(
            x_rel,
            y_rel,
            s=48.0,
            marker="o",
            facecolors=color,
            edgecolors="white",
            alpha=0.98,
            linewidths=0.8,
            zorder=5,
        )

        if len(x_rel) >= 1:
            ax.scatter(
                [x_rel[0]],
                [y_rel[0]],
                s=72.0,
                marker="v",
                color=color,
                edgecolors="black",
                linewidths=0.7,
                zorder=6,
            )

        if len(x_rel) >= 2:
            ax.annotate(
                "",
                xy=(x_rel[1], y_rel[1]),
                xytext=(x_rel[0], y_rel[0]),
                arrowprops=dict(arrowstyle="->", color=color, lw=1.0, alpha=0.95),
                zorder=6,
            )

        anchor_x, anchor_y = _choose_label_anchor(x_ref, y_ref, x_rel, y_rel)
        tangential_sign = 1.0 if (ig % 2 == 0) else -1.0
        lx, ly, ha, va = _choose_label_position(anchor_x, anchor_y, span, radial_scale=0.095, tangential_scale=0.055, tangential_sign=tangential_sign)
        ax.annotate(
            _label_text(key, sub),
            xy=(anchor_x, anchor_y),
            xytext=(lx, ly),
            textcoords="data",
            color=color,
            fontsize=9.0,
            ha=ha,
            va=va,
            bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor=color, alpha=0.84, linewidth=0.75),
            arrowprops=dict(arrowstyle="-", color=color, lw=0.8, alpha=0.80, shrinkA=3, shrinkB=3),
            zorder=7,
        )
        if np.isfinite(x_ref) and np.isfinite(y_ref):
            rx, ry, rha, rva = _choose_label_position(x_ref, y_ref, span, radial_scale=0.040, tangential_scale=0.030, tangential_sign=-tangential_sign)
            ax.annotate(
                "ref (El=0)",
                xy=(x_ref, y_ref),
                xytext=(rx, ry),
                textcoords="data",
                color=color,
                fontsize=7.6,
                ha=rha,
                va=rva,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor=color, alpha=0.78, linewidth=0.55),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.7, alpha=0.70, shrinkA=2, shrinkB=2),
                xycoords="data",
                zorder=7,
            )

    ax.axhline(0.0, color="0.55", linewidth=0.8, alpha=0.6, zorder=0)
    ax.axvline(0.0, color="0.55", linewidth=0.8, alpha=0.6, zorder=0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [arcsec]")
    ax.set_ylabel("y [arcsec]")
    ax.set_title(_friendly_title(title, Path(output_path)))
    ax.grid(True, alpha=0.22)

    legend_handles = [
        mlines.Line2D([], [], color="0.45", marker="o", linestyle="-", markersize=6, label="measured beam centers"),
        mlines.Line2D([], [], color="0.45", marker="x", linestyle="--", markersize=6, label="fitted model at sampled El"),
        mlines.Line2D([], [], color="0.25", marker="v", linestyle="None", markersize=7, label="lower-El end of measured track"),
        mlines.Line2D([], [], color="0.25", marker=None, linestyle="-", linewidth=1.0, label="arrow shows increasing El"),
        mlines.Line2D([], [], color="0.25", marker="*", linestyle="None", markersize=10, label="reference position / model point at El = 0°"),
    ]
    if dummy_zero90 is not None:
        legend_handles.append(
            mlines.Line2D([], [], color="0.45", linestyle=":", linewidth=1.1, alpha=0.35, label="thin model track for El = 0° to 90°")
        )
    ax.legend(handles=legend_handles, loc="best", fontsize=8)


    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=int(dpi))
    plt.close(fig)
    return output_path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Plot beam-fit residuals in the x-y plane with beam labels and El-direction cues.")
    p.add_argument("csv", help="Input residual CSV, e.g. beam_fit_residuals_center_beam.csv")
    p.add_argument("--output", default=None, help="Output PNG path. Default: <input_stem>_xy.png")
    p.add_argument("--title", default=None, help="Figure title")
    p.add_argument("--dpi", type=int, default=150, help="PNG DPI")
    p.add_argument("--connect-mode", choices=["none", "line", "arrow"], default="line", help="How to connect pred and rel for each row")
    return p


def main() -> int:
    args = _build_parser().parse_args()
    csv_path = Path(args.csv).expanduser().resolve()
    if not csv_path.exists():
        raise FileNotFoundError(f"input CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    output = Path(args.output).expanduser().resolve() if args.output is not None else csv_path.with_name(f"{csv_path.stem}_xy.png")
    write_beam_fit_xy_png(df, output, title=args.title, connect_mode=args.connect_mode, dpi=int(args.dpi))
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
