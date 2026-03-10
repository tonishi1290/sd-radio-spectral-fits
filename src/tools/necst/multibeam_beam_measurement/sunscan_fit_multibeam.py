from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence
import builtins

import pandas as pd

from .beam_geometry import beam_rows_to_dataframe, fit_beam_model
from .config_io import load_spectrometer_config, resolve_primary_stream_names, write_beam_model_toml


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
        beam_df = beam_rows_to_dataframe(model_result.beam_rows)
        shifts_df = model_result.shifts_df.copy()
        residual_df = model_result.usable_df.copy()
        rejected_df = model_result.rejected_df.copy()

        beam_csv = outdir / f"beam_fit_results_{model_name}.csv"
        shift_csv = outdir / f"run_shift_results_{model_name}.csv"
        resid_csv = outdir / f"beam_fit_residuals_{model_name}.csv"
        reject_csv = outdir / f"beam_fit_rejected_{model_name}.csv"
        beam_df.to_csv(beam_csv, index=False)
        shifts_df.to_csv(shift_csv, index=False)
        residual_df.to_csv(resid_csv, index=False)
        rejected_df.to_csv(reject_csv, index=False)
        outputs[f"beam_csv_{model_name}"] = beam_csv
        outputs[f"shift_csv_{model_name}"] = shift_csv
        outputs[f"resid_csv_{model_name}"] = resid_csv
        outputs[f"rejected_csv_{model_name}"] = reject_csv

        toml_path = outdir / f"beam_model_{model_name}.toml"
        write_beam_model_toml(toml_path, Path(spectrometer_config), beam_df, model_name=model_name)
        outputs[f"beam_toml_{model_name}"] = toml_path

        metric = float(beam_df["rms_resid_arcsec"].mean()) if not beam_df.empty else float("inf")
        if metric < best_metric:
            best_metric = metric
            best_model_name = model_name

        fit_summary_lines.append("")
        fit_summary_lines.append(f"[{model_name}]")
        fit_summary_lines.append(f"reference_angle_deg = {model_result.reference_angle_deg:.6f}")
        fit_summary_lines.append(f"rotation_sign = {model_result.rotation_sign:.6f}")
        fit_summary_lines.append(f"rotation_slope_deg_per_deg = {model_result.rotation_slope_deg_per_deg:.6f}")
        fit_summary_lines.append(f"dewar_angle_deg = {model_result.dewar_angle_deg:.6f}")
        fit_summary_lines.append(f"mean_rms_resid_arcsec = {metric:.6f}")
        fit_summary_lines.append(f"beam_rows = {len(beam_df)}")
        fit_summary_lines.append(f"selected_rows = {len(residual_df)}")
        fit_summary_lines.append(f"rejected_rows = {len(rejected_df)}")

    fit_summary_lines.append("")
    fit_summary_lines.append(f"recommended_model = {best_model_name}")
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
    ap.add_argument("--reference-angle-deg", type=float, default=None, help="Override reference elevation angle in degrees")
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
