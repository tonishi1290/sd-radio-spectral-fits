from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence
import builtins
import math

import pandas as pd

from .beam_geometry import beam_rows_to_dataframe, fit_beam_model
from .config_io import load_spectrometer_config, resolve_primary_stream_names, write_beam_model_toml
from .beam_rotation_shared import rotate_offset_arcsec
from .plot_beam_fit_residuals_xy import write_beam_fit_xy_png


MODEL_CHOICES = ["both", "center_beam", "virtual_center"]


def _load_summary_tables(paths: Sequence[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in paths]
    if not frames:
        raise ValueError("no input summary tables were provided")
    return pd.concat(frames, axis=0, ignore_index=True)



def _coerce_bool_series(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip().str.lower()
    return text.isin(["1", "true", "t", "yes", "y"])



def _usable_rows(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    for col in ["fit_ok_az", "fit_ok_el"]:
        if col in work.columns:
            work = work.loc[_coerce_bool_series(work[col])]
    for col in ["x_arcsec", "y_arcsec", "rep_el_deg", "beam_id"]:
        if col not in work.columns:
            raise ValueError(f"required column missing from summary table: {col!r}")
    work["beam_id"] = work["beam_id"].astype(str)
    work["stream_name"] = work.get("stream_name", pd.Series([None] * len(work))).astype(str)
    return work



def _singlebeam_compat_outputs(work: pd.DataFrame, *, outdir: Path, spectrometer_config: Path) -> Dict[str, Path]:
    outputs: Dict[str, Path] = {}
    selected_csv = outdir / "selected_rows_singlebeam_compat.csv"
    work.to_csv(selected_csv, index=False)
    outputs["selected_rows_csv"] = selected_csv
    txt = outdir / "fit_summary.txt"
    beam_ids = ", ".join(sorted(work["beam_id"].astype(str).unique()))
    lines = [
        "Multi-beam beam-geometry fit summary",
        "",
        "single-beam compatibility mode was triggered because fewer than two beam IDs were present.",
        f"input rows = {len(work)}",
        f"beam_ids = {beam_ids}",
        f"spectrometer_config = {spectrometer_config}",
        "",
        "No geometry fit was executed. This is expected when validating the package with single-beam data.",
    ]
    txt.write_text("\n".join(lines) + "\n", encoding="utf-8")
    outputs["fit_summary_txt"] = txt
    return outputs



def _theta_map_from_rep_el(rep_el_by_scan: Dict[str, float], rotation_sign: float) -> Dict[str, float]:
    sign = float(rotation_sign)
    return {str(scan_key): sign * float(rep_el) for scan_key, rep_el in rep_el_by_scan.items()}



def _template_from_points(points_df: pd.DataFrame, theta_map_deg: Dict[str, float], *, center_beam_id: Optional[str], model_name: str) -> Dict[str, complex]:
    template: Dict[str, complex] = {}
    for beam_id, subset in points_df.groupby(points_df["beam_id"].astype(str), sort=True):
        vals = []
        for row in subset.itertuples(index=False):
            z = complex(float(row.x_rel_arcsec), float(row.y_rel_arcsec))
            theta_deg = float(theta_map_deg.get(str(row.scan_key), 0.0))
            theta_rad = math.radians(theta_deg)
            vals.append(z * complex(math.cos(-theta_rad), math.sin(-theta_rad)))
        if vals:
            template[str(beam_id)] = sum(vals) / float(len(vals))
    if model_name == "center_beam" and center_beam_id is not None:
        template[str(center_beam_id)] = 0.0 + 0.0j
    return template



def _recompute_converter_compatible_outputs(
    model_result,
    beam_df: pd.DataFrame,
    shifts_df: pd.DataFrame,
    residual_df: pd.DataFrame,
    *,
    center_beam_id: Optional[str],
    model_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    raw_reference = float(model_result.reference_angle_deg)
    raw_sign = float(model_result.rotation_sign)
    raw_slope = float(model_result.rotation_slope_deg_per_deg)
    raw_dewar = float(model_result.dewar_angle_deg)

    points = residual_df.copy()
    rep_el_by_scan: Dict[str, float] = {}
    for scan_key, subset in points.groupby(points["scan_key"].astype(str), sort=False):
        rep_el_by_scan[str(scan_key)] = float(pd.to_numeric(subset["rep_el_deg"], errors="coerce").iloc[0])

    theta_map_deg = _theta_map_from_rep_el(rep_el_by_scan, raw_sign)
    template = _template_from_points(points, theta_map_deg, center_beam_id=center_beam_id, model_name=model_name)

    norm_beam_df = beam_df.copy()
    norm_beam_df["raw_reference_angle_deg"] = raw_reference
    norm_beam_df["raw_rotation_slope_deg_per_deg"] = raw_slope
    norm_beam_df["raw_dewar_angle_deg"] = raw_dewar

    beam_counts = points.groupby(points["beam_id"].astype(str)).size().to_dict()
    beam_scan_counts = points.groupby(points["beam_id"].astype(str))["scan_key"].nunique().to_dict()
    beam_run_counts = points.groupby(points["beam_id"].astype(str))["run_id"].nunique().to_dict()

    for idx, row in norm_beam_df.iterrows():
        beam_id = str(row["beam_id"])
        z0 = template.get(beam_id, 0.0 + 0.0j)
        az = float(z0.real)
        el = float(z0.imag)
        subset = points.loc[points["beam_id"].astype(str) == beam_id].copy()
        resids = []
        for prow in subset.itertuples(index=False):
            theta_deg = float(theta_map_deg.get(str(prow.scan_key), 0.0))
            pred_x, pred_y = rotate_offset_arcsec(az, el, theta_deg)
            pred_x = float(pred_x)
            pred_y = float(pred_y)
            resids.append(math.hypot(float(prow.x_rel_arcsec) - pred_x, float(prow.y_rel_arcsec) - pred_y))
        norm_beam_df.at[idx, "fit_n"] = int(beam_counts.get(beam_id, len(subset)))
        norm_beam_df.at[idx, "n_scan"] = int(beam_scan_counts.get(beam_id, subset["scan_key"].nunique())) if not subset.empty else 0
        norm_beam_df.at[idx, "n_run"] = int(beam_run_counts.get(beam_id, subset["run_id"].nunique())) if not subset.empty else 0
        norm_beam_df.at[idx, "az_offset_arcsec"] = az
        norm_beam_df.at[idx, "el_offset_arcsec"] = el
        norm_beam_df.at[idx, "radius_arcsec"] = float(abs(z0))
        norm_beam_df.at[idx, "phase_deg"] = float(math.degrees(math.atan2(el, az))) if abs(z0) > 0 else 0.0
        norm_beam_df.at[idx, "reference_angle_deg"] = 0.0
        norm_beam_df.at[idx, "rotation_sign"] = raw_sign
        norm_beam_df.at[idx, "rotation_slope_deg_per_deg"] = raw_sign
        norm_beam_df.at[idx, "dewar_angle_deg"] = 0.0
        norm_beam_df.at[idx, "rms_resid_arcsec"] = float((sum(r * r for r in resids) / len(resids)) ** 0.5) if resids else float("nan")
        norm_beam_df.at[idx, "max_resid_arcsec"] = float(max(resids)) if resids else float("nan")

    norm_residual_df = residual_df.copy()
    norm_residual_df["raw_pred_x_arcsec"] = norm_residual_df.get("pred_x_arcsec")
    norm_residual_df["raw_pred_y_arcsec"] = norm_residual_df.get("pred_y_arcsec")
    norm_residual_df["raw_resid_arcsec"] = norm_residual_df.get("resid_arcsec")
    norm_residual_df["raw_reference_angle_deg"] = raw_reference
    norm_residual_df["raw_rotation_slope_deg_per_deg"] = raw_slope
    norm_residual_df["raw_dewar_angle_deg"] = raw_dewar

    beam_lookup = {
        str(row["beam_id"]): (float(row["az_offset_arcsec"]), float(row["el_offset_arcsec"]))
        for _, row in norm_beam_df.iterrows()
    }
    pred_x_vals = []
    pred_y_vals = []
    resid_vals = []
    for row in norm_residual_df.itertuples(index=False):
        az, el = beam_lookup[str(row.beam_id)]
        theta_deg = float(theta_map_deg.get(str(row.scan_key), 0.0))
        pred_x, pred_y = rotate_offset_arcsec(az, el, theta_deg)
        pred_x = float(pred_x)
        pred_y = float(pred_y)
        pred_x_vals.append(pred_x)
        pred_y_vals.append(pred_y)
        resid_vals.append(math.hypot(float(row.x_rel_arcsec) - pred_x, float(row.y_rel_arcsec) - pred_y))
    norm_residual_df["pred_x_arcsec"] = pred_x_vals
    norm_residual_df["pred_y_arcsec"] = pred_y_vals
    norm_residual_df["resid_arcsec"] = resid_vals
    norm_residual_df["reference_angle_deg"] = 0.0
    norm_residual_df["rotation_sign"] = raw_sign
    norm_residual_df["rotation_slope_deg_per_deg"] = raw_sign
    norm_residual_df["dewar_angle_deg"] = 0.0

    norm_shifts_df = shifts_df.copy()
    if "theta_deg" in norm_shifts_df.columns:
        norm_shifts_df["raw_theta_deg"] = norm_shifts_df["theta_deg"]
    norm_shifts_df["raw_reference_angle_deg"] = raw_reference
    norm_shifts_df["raw_rotation_slope_deg_per_deg"] = raw_slope
    norm_shifts_df["raw_dewar_angle_deg"] = raw_dewar
    norm_shifts_df["theta_deg"] = [float(theta_map_deg.get(str(scan_key), 0.0)) for scan_key in norm_shifts_df["scan_key"].astype(str)]
    norm_shifts_df["reference_angle_deg"] = 0.0
    norm_shifts_df["rotation_sign"] = raw_sign
    norm_shifts_df["rotation_slope_deg_per_deg"] = raw_sign
    norm_shifts_df["dewar_angle_deg"] = 0.0

    info = {
        "raw_reference_angle_deg": raw_reference,
        "raw_rotation_slope_deg_per_deg": raw_slope,
        "raw_dewar_angle_deg": raw_dewar,
        "normalized_reference_angle_deg": 0.0,
        "normalized_dewar_angle_deg": 0.0,
        "normalized_rotation_slope_deg_per_deg": raw_sign,
        "rotation_sign": raw_sign,
    }
    return norm_beam_df, norm_shifts_df, norm_residual_df, info



def run_fit(
    summary_paths: Sequence[Path],
    spectrometer_config: Path,
    *,
    outdir: Path,
    center_beam_id: Optional[str] = None,
    stream_names: Optional[Sequence[str]] = None,
    model: str = "both",
    reference_angle_deg: Optional[float] = None,
    sigma_clip: Optional[float] = 4.5,
    clip_iters: int = 2,
    min_points_per_beam: int = 2,
    min_scans_per_beam: int = 2,
) -> Dict[str, Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary_df = _load_summary_tables(summary_paths)
    summary_df = _usable_rows(summary_df)
    cfg = load_spectrometer_config(spectrometer_config)
    primary_streams = resolve_primary_stream_names(Path(spectrometer_config), cfg, explicit_stream_names=stream_names)
    work = summary_df.loc[summary_df["stream_name"].astype(str).isin([str(s) for s in primary_streams])].copy()
    if work.empty:
        raise ValueError("no usable rows remain after applying stream selection")

    selected_rows_csv = outdir / "selected_rows_for_fit.csv"
    work.to_csv(selected_rows_csv, index=False)
    outputs: Dict[str, Path] = {"selected_rows_csv": selected_rows_csv}

    per_beam_point_counts = work.groupby("beam_id").size() if "beam_id" in work.columns else pd.Series(dtype=int)
    per_beam_scan_counts = work.groupby("beam_id")["scan_id"].nunique() if "scan_id" in work.columns else pd.Series(dtype=int)
    effective_min_points_per_beam = int(min_points_per_beam)
    effective_min_scans_per_beam = int(min_scans_per_beam)
    auto_relaxed_min_points = False
    auto_relaxed_min_scans = False
    if not per_beam_point_counts.empty:
        max_points = int(per_beam_point_counts.max())
        if effective_min_points_per_beam > max_points:
            effective_min_points_per_beam = builtins.max(1, max_points)
            auto_relaxed_min_points = True
    if not per_beam_scan_counts.empty:
        max_scans = int(per_beam_scan_counts.max())
        if effective_min_scans_per_beam > max_scans:
            effective_min_scans_per_beam = builtins.max(1, max_scans)
            auto_relaxed_min_scans = True

    if int(work["beam_id"].astype(str).nunique()) < 2:
        compat_outputs = _singlebeam_compat_outputs(work, outdir=outdir, spectrometer_config=spectrometer_config)
        outputs.update(compat_outputs)
        return outputs

    if model not in MODEL_CHOICES:
        raise ValueError(f"unsupported model={model!r}; expected one of {MODEL_CHOICES}")
    if model == "both":
        model_names = ["virtual_center"] if center_beam_id is None else ["center_beam", "virtual_center"]
    elif model == "center_beam":
        if center_beam_id is None:
            raise ValueError("model='center_beam' requires center_beam_id")
        model_names = ["center_beam"]
    else:
        model_names = ["virtual_center"]

    best_model_name: Optional[str] = None
    best_metric = float("inf")
    fit_summary_lines: List[str] = []
    fit_summary_lines.append("Multi-beam beam-geometry fit summary")
    fit_summary_lines.append(f"input rows = {len(work)}")
    fit_summary_lines.append(f"selected primary streams = {', '.join(primary_streams)}")
    fit_summary_lines.append(f"center_beam_id = {center_beam_id}")
    fit_summary_lines.append(f"sigma_clip = {sigma_clip}")
    fit_summary_lines.append(f"clip_iters = {clip_iters}")
    fit_summary_lines.append(f"min_points_per_beam = {min_points_per_beam}")
    fit_summary_lines.append(f"effective_min_points_per_beam = {effective_min_points_per_beam}")
    fit_summary_lines.append(f"min_scans_per_beam = {min_scans_per_beam}")
    fit_summary_lines.append(f"effective_min_scans_per_beam = {effective_min_scans_per_beam}")
    if auto_relaxed_min_points:
        fit_summary_lines.append("auto_relax_note_points = minimum point count threshold was relaxed because the input contains fewer points per beam than requested")
    if auto_relaxed_min_scans:
        fit_summary_lines.append("auto_relax_note_scans = minimum scan count threshold was relaxed because the input contains fewer scans per beam than requested")

    for model_name in model_names:
        model_result = fit_beam_model(
            work,
            model_name=model_name,
            center_beam_id=center_beam_id,
            reference_angle_deg=reference_angle_deg,
            sigma_clip=sigma_clip,
            clip_iters=clip_iters,
            min_points_per_beam=effective_min_points_per_beam,
            min_scans_per_beam=effective_min_scans_per_beam,
        )
        raw_beam_df = beam_rows_to_dataframe(model_result.beam_rows)
        raw_shifts_df = model_result.shifts_df.copy()
        raw_residual_df = model_result.usable_df.copy()
        rejected_df = model_result.rejected_df.copy()

        beam_df = raw_beam_df.copy()
        shifts_df = raw_shifts_df.copy()
        residual_df = raw_residual_df.copy()
        norm_info = {
            "raw_reference_angle_deg": float(model_result.reference_angle_deg),
            "raw_rotation_slope_deg_per_deg": float(model_result.rotation_slope_deg_per_deg),
            "raw_dewar_angle_deg": float(model_result.dewar_angle_deg),
            "normalized_reference_angle_deg": float(model_result.reference_angle_deg),
            "normalized_rotation_slope_deg_per_deg": float(model_result.rotation_slope_deg_per_deg),
            "normalized_dewar_angle_deg": float(model_result.dewar_angle_deg),
            "rotation_sign": float(model_result.rotation_sign),
        }
        if not rejected_df.empty:
            rejected_df = rejected_df.copy()
            rejected_df["reference_angle_deg"] = float(model_result.reference_angle_deg)
            rejected_df["rotation_sign"] = float(model_result.rotation_sign)
            rejected_df["rotation_slope_deg_per_deg"] = float(model_result.rotation_slope_deg_per_deg)
            rejected_df["dewar_angle_deg"] = float(model_result.dewar_angle_deg)

        beam_csv = outdir / f"beam_fit_results_{model_name}.csv"
        shift_csv = outdir / f"run_shift_results_{model_name}.csv"
        resid_csv = outdir / f"beam_fit_residuals_{model_name}.csv"
        reject_csv = outdir / f"beam_fit_rejected_{model_name}.csv"
        png_xy = outdir / f"beam_fit_residuals_{model_name}_xy.png"
        beam_df.to_csv(beam_csv, index=False)
        shifts_df.to_csv(shift_csv, index=False)
        residual_df.to_csv(resid_csv, index=False)
        rejected_df.to_csv(reject_csv, index=False)
        write_beam_fit_xy_png(residual_df, png_xy, title=f"{model_name}: measured beam tracks and fitted model")
        outputs[f"beam_csv_{model_name}"] = beam_csv
        outputs[f"shift_csv_{model_name}"] = shift_csv
        outputs[f"resid_csv_{model_name}"] = resid_csv
        outputs[f"rejected_csv_{model_name}"] = reject_csv
        outputs[f"resid_xy_png_{model_name}"] = png_xy

        toml_path = outdir / f"beam_model_{model_name}.toml"
        write_beam_model_toml(toml_path, Path(spectrometer_config), beam_df, model_name=model_name)
        outputs[f"beam_toml_{model_name}"] = toml_path

        fit_summary_lines.append("")
        fit_summary_lines.append(f"[{model_name}]")
        fit_summary_lines.append(f"raw_reference_angle_deg = {model_result.reference_angle_deg:.6f}")
        fit_summary_lines.append(f"raw_rotation_sign = {model_result.rotation_sign:.6f}")
        fit_summary_lines.append(f"raw_rotation_slope_deg_per_deg = {model_result.rotation_slope_deg_per_deg:.6f}")
        fit_summary_lines.append(f"raw_dewar_angle_deg = {model_result.dewar_angle_deg:.6f}")
        fit_summary_lines.append("output_normalization = raw_fit_parameters_are_exported_without_converter_renormalization")
        fit_summary_lines.append(f"reference_angle_deg = {norm_info['normalized_reference_angle_deg']:.6f}")
        fit_summary_lines.append(f"rotation_sign = {norm_info['rotation_sign']:.6f}")
        fit_summary_lines.append(f"rotation_slope_deg_per_deg = {norm_info['normalized_rotation_slope_deg_per_deg']:.6f}")
        fit_summary_lines.append(f"dewar_angle_deg = {norm_info['normalized_dewar_angle_deg']:.6f}")
        fit_summary_lines.append(f"n_beams = {len(beam_df)}")
        fit_summary_lines.append(f"n_rows_used = {len(residual_df)}")
        fit_summary_lines.append(f"n_rows_rejected = {len(rejected_df)}")
        fit_summary_lines.append(f"mean_resid_arcsec = {pd.to_numeric(residual_df['resid_arcsec'], errors='coerce').mean():.6f}")
        fit_summary_lines.append(f"max_resid_arcsec = {pd.to_numeric(residual_df['resid_arcsec'], errors='coerce').max():.6f}")
        metric = float(pd.to_numeric(residual_df["resid_arcsec"], errors="coerce").mean())
        if metric < best_metric:
            best_metric = metric
            best_model_name = model_name

    if best_model_name is not None:
        fit_summary_lines.append("")
        fit_summary_lines.append(f"best_model = {best_model_name}")
        fit_summary_lines.append(f"best_metric_mean_resid_arcsec = {best_metric:.6f}")

    summary_txt = outdir / "fit_summary.txt"
    summary_txt.write_text("\n".join(fit_summary_lines) + "\n", encoding="utf-8")
    outputs["fit_summary_txt"] = summary_txt
    return outputs



def add_fit_arguments(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("summary_csv", nargs="+", help="one or more all-stream summary CSV files")
    ap.add_argument("--spectrometer-config", required=True, help="converter-compatible spectrometer config TOML")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--center-beam-id", default=None, help="Center beam ID for center-beam model")
    ap.add_argument("--fit-stream-name", dest="stream_names", action="append", default=None, help="Select a stream name for fitting (repeatable)")
    ap.add_argument("--model", choices=MODEL_CHOICES, default="both", help="Fit model selection")
    ap.add_argument("--reference-angle-deg", type=float, default=None, help="Override raw fit reference elevation angle in degrees before output normalization")
    ap.add_argument("--sigma-clip", type=float, default=4.5, help="Sigma threshold for residual clipping; <=0 disables clipping")
    ap.add_argument("--clip-iters", type=int, default=2, help="Number of sigma-clipping iterations")
    ap.add_argument("--min-points-per-beam", type=int, default=2, help="Minimum retained points per beam")
    ap.add_argument("--min-scans-per-beam", type=int, default=2, help="Minimum retained scans per beam")



def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Fit a multi-beam geometry model from aggregated sun_scan summaries.")
    add_fit_arguments(ap)
    args = ap.parse_args(argv)
    outputs = run_fit(
        summary_paths=[Path(p).expanduser().resolve() for p in args.summary_csv],
        spectrometer_config=Path(args.spectrometer_config).expanduser().resolve(),
        outdir=Path(args.outdir).expanduser().resolve(),
        center_beam_id=args.center_beam_id,
        stream_names=args.stream_names,
        model=str(args.model),
        reference_angle_deg=args.reference_angle_deg,
        sigma_clip=(None if args.sigma_clip is None or float(args.sigma_clip) <= 0 else float(args.sigma_clip)),
        clip_iters=int(args.clip_iters),
        min_points_per_beam=int(args.min_points_per_beam),
        min_scans_per_beam=int(args.min_scans_per_beam),
    )
    for key, path in outputs.items():
        print(f"[done] {key}: {path}")


if __name__ == "__main__":
    main()
