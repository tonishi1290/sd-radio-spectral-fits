from __future__ import annotations

from dataclasses import dataclass

from .io_utils import read_csv


@dataclass
class InspectOptions:
    input_csv: str


def _unique_nonempty(df, col: str) -> list[str]:
    if col not in df.columns:
        return []
    vals = []
    for x in df[col].astype(str).unique().tolist():
        sx = str(x)
        if sx and sx != "nan":
            vals.append(sx)
    return sorted(set(vals))


def run_inspect(opts: InspectOptions) -> str:
    df = read_csv(opts.input_csv)
    lines = []
    lines.append(f"Rows                : {len(df)}")
    lines.append(f"Columns             : {len(df.columns)}")
    lines.append(f"Column names        : {', '.join(df.columns)}")
    if "dataset_id" in df.columns:
        counts = df["dataset_id"].astype(str).value_counts().sort_index()
        lines.append("Dataset counts:")
        for k, v in counts.items():
            lines.append(f"  {k}: {v}")
    if "dt_cross_id" in df.columns:
        lines.append(f"Unique dt_cross_id  : {df['dt_cross_id'].nunique(dropna=False)}")
    if "used_param_hash" in df.columns:
        vals = _unique_nonempty(df, "used_param_hash")
        lines.append(f"used_param_hash set : {len(vals)}")
    if "input_angle_unit" in df.columns:
        vals = _unique_nonempty(df, "input_angle_unit")
        lines.append(f"input_angle_unit    : {', '.join(vals)}")
        lines.append("canonical Az/El     : stored internally as az_deg/el_deg")
    for col, label in [
        ("measurement_space", "measurement_space  "),
        ("positive_az_moves_star", "+Az moves star     "),
        ("positive_el_moves_star", "+El moves star     "),
        ("command_err_mode", "command_err_mode   "),
    ]:
        vals = _unique_nonempty(df, col)
        if vals:
            lines.append(f"{label}: {', '.join(vals)}")
    if "dx_to_model_sign" in df.columns or "dy_to_model_sign" in df.columns:
        dxv = _unique_nonempty(df, "dx_to_model_sign")
        dyv = _unique_nonempty(df, "dy_to_model_sign")
        if dxv:
            lines.append(f"dx_to_model_sign   : {', '.join(dxv)}")
        if dyv:
            lines.append(f"dy_to_model_sign   : {', '.join(dyv)}")
    absolute_ok = all(c in df.columns for c in ["applied_model_dx_arcsec", "applied_model_dy_arcsec"])
    lines.append(f"Absolute fit ready  : {absolute_ok}")
    missing = df.isna().mean().sort_values(ascending=False)
    lines.append("Missing fraction:")
    for k, v in missing.items():
        lines.append(f"  {k}: {v:.3f}")
    return "\n".join(lines)
