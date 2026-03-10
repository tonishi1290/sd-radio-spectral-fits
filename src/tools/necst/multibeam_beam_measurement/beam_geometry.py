from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import builtins
import math

import numpy as np
import pandas as pd

from .beam_rotation_shared import rotate_offset_arcsec


@dataclass
class BeamFitRow:
    model_name: str
    beam_id: str
    fit_n: int
    n_scan: int
    n_run: int
    az_offset_arcsec: float
    el_offset_arcsec: float
    radius_arcsec: float
    phase_deg: float
    rotation_mode: str
    reference_angle_deg: float
    rotation_sign: float
    rotation_slope_deg_per_deg: float
    dewar_angle_deg: float
    rms_resid_arcsec: float
    max_resid_arcsec: float


@dataclass
class BeamFitModelResult:
    model_name: str
    reference_angle_deg: float
    rotation_sign: float
    rotation_slope_deg_per_deg: float
    dewar_angle_deg: float
    beam_rows: List[BeamFitRow]
    shifts_df: pd.DataFrame
    usable_df: pd.DataFrame
    selected_df: pd.DataFrame
    rejected_df: pd.DataFrame



def _safe_float_array(series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").to_numpy(float)



def _estimate_scan_rotation(template: Dict[str, complex], scan_points: Dict[str, complex]) -> float:
    keys = [k for k in scan_points.keys() if k in template and np.isfinite(template[k].real) and np.isfinite(scan_points[k].real)]
    if not keys:
        return 0.0
    acc = 0.0 + 0.0j
    for key in keys:
        acc += scan_points[key] * np.conjugate(template[key])
    if abs(acc) <= 0:
        return 0.0
    return float(np.rad2deg(np.angle(acc)))



def _linear_fit_theta(rep_el_deg: np.ndarray, theta_deg: np.ndarray, reference_angle_deg: float) -> Tuple[float, float]:
    x = np.asarray(rep_el_deg, dtype=float) - float(reference_angle_deg)
    y = np.asarray(theta_deg, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if int(np.count_nonzero(mask)) >= 2:
        slope, intercept = np.polyfit(x[mask], y[mask], 1)
        return float(slope), float(intercept)
    if int(np.count_nonzero(mask)) == 1:
        return 0.0, float(y[mask][0])
    return 0.0, 0.0


def _quantize_rotation_sign_and_intercept(rep_el_deg: np.ndarray, theta_deg: np.ndarray, reference_angle_deg: float, slope: float, intercept: float) -> Tuple[float, float]:
    x = np.asarray(rep_el_deg, dtype=float) - float(reference_angle_deg)
    y = np.asarray(theta_deg, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if int(np.count_nonzero(mask)) < 2:
        return 0.0, float(intercept)

    if float(np.nanmax(x[mask]) - np.nanmin(x[mask])) <= 0.0:
        return 0.0, float(intercept)

    sign = 1.0 if float(slope) >= 0.0 else -1.0
    dewar = float(np.nanmean(y[mask] - sign * x[mask]))
    return sign, dewar



def _iterative_rigid_rotation_fit(points_df: pd.DataFrame, reference_angle_deg: float) -> Tuple[Dict[str, complex], pd.DataFrame, float, float]:
    scans = sorted(points_df["scan_key"].unique())
    beam_ids = sorted(points_df["beam_id"].unique())
    theta_map = {scan_key: 0.0 for scan_key in scans}
    slope = 0.0
    intercept = 0.0

    for _ in range(12):
        template: Dict[str, complex] = {}
        for beam_id in beam_ids:
            subset = points_df.loc[points_df["beam_id"] == beam_id]
            vals: List[complex] = []
            for _, row in subset.iterrows():
                z = complex(float(row["x_rel_arcsec"]), float(row["y_rel_arcsec"]))
                theta = math.radians(theta_map[str(row["scan_key"])])
                vals.append(z * complex(math.cos(-theta), math.sin(-theta)))
            if vals:
                template[str(beam_id)] = sum(vals) / float(len(vals))

        est_rows: List[dict] = []
        for scan_key in scans:
            subset = points_df.loc[points_df["scan_key"] == scan_key]
            scan_points = {
                str(row["beam_id"]): complex(float(row["x_rel_arcsec"]), float(row["y_rel_arcsec"]))
                for _, row in subset.iterrows()
            }
            theta_est = _estimate_scan_rotation(template, scan_points)
            est_rows.append({
                "scan_key": str(scan_key),
                "theta_est_deg": theta_est,
                "rep_el_deg": float(subset["rep_el_deg"].iloc[0]),
            })
        theta_df = pd.DataFrame(est_rows)
        slope, intercept = _linear_fit_theta(_safe_float_array(theta_df["rep_el_deg"]), _safe_float_array(theta_df["theta_est_deg"]), reference_angle_deg)
        for _, row in theta_df.iterrows():
            theta_map[str(row["scan_key"])] = float(slope) * (float(row["rep_el_deg"]) - float(reference_angle_deg)) + float(intercept)

    final_template: Dict[str, complex] = {}
    for beam_id in beam_ids:
        subset = points_df.loc[points_df["beam_id"] == beam_id]
        vals: List[complex] = []
        for _, row in subset.iterrows():
            theta = math.radians(theta_map[str(row["scan_key"])])
            z = complex(float(row["x_rel_arcsec"]), float(row["y_rel_arcsec"]))
            vals.append(z * complex(math.cos(-theta), math.sin(-theta)))
        if vals:
            final_template[str(beam_id)] = sum(vals) / float(len(vals))

    shift_rows = []
    for scan_key in scans:
        subset = points_df.loc[points_df["scan_key"] == scan_key]
        shift_rows.append({
            "scan_key": str(scan_key),
            "rep_el_deg": float(subset["rep_el_deg"].iloc[0]),
            "theta_deg": float(theta_map[str(scan_key)]),
        })
    return final_template, pd.DataFrame(shift_rows), float(slope), float(intercept)



def _build_relative_points(df: pd.DataFrame, model_name: str, center_beam_id: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    records: List[dict] = []
    shift_records: List[dict] = []
    grouped = df.groupby(["run_id", "scan_id"], sort=True)
    for (run_id, scan_id), grp in grouped:
        grp = grp.copy()
        grp = grp.loc[np.isfinite(_safe_float_array(grp["x_arcsec"])) & np.isfinite(_safe_float_array(grp["y_arcsec"])) & np.isfinite(_safe_float_array(grp["rep_el_deg"]))]
        if grp.empty:
            continue
        scan_key = f"{run_id}::{int(scan_id)}"
        rep_el_deg = float(np.nanmedian(_safe_float_array(grp["rep_el_deg"])))
        if model_name == "center_beam":
            if center_beam_id is None:
                raise ValueError("center_beam model requires center_beam_id")
            center_rows = grp.loc[grp["beam_id"].astype(str) == str(center_beam_id)]
            if center_rows.empty:
                continue
            x0 = float(_safe_float_array(center_rows["x_arcsec"])[0])
            y0 = float(_safe_float_array(center_rows["y_arcsec"])[0])
        elif model_name == "virtual_center":
            x0 = float(np.nanmean(_safe_float_array(grp["x_arcsec"])))
            y0 = float(np.nanmean(_safe_float_array(grp["y_arcsec"])))
        else:
            raise ValueError(f"unsupported model_name={model_name!r}")
        shift_records.append({
            "model_name": model_name,
            "run_id": str(run_id),
            "scan_id": int(scan_id),
            "scan_key": scan_key,
            "rep_el_deg": rep_el_deg,
            "shift_x_arcsec": x0,
            "shift_y_arcsec": y0,
        })
        for _, row in grp.iterrows():
            x_val = float(row["x_arcsec"])
            y_val = float(row["y_arcsec"])
            records.append({
                "model_name": model_name,
                "run_id": str(run_id),
                "scan_id": int(scan_id),
                "scan_key": scan_key,
                "beam_id": str(row["beam_id"]),
                "stream_name": row.get("stream_name", None),
                "rep_el_deg": rep_el_deg,
                "x_rel_arcsec": x_val - x0,
                "y_rel_arcsec": y_val - y0,
                "x_arcsec": x_val,
                "y_arcsec": y_val,
            })
    return pd.DataFrame(records), pd.DataFrame(shift_records)



def _residual_table(points_df: pd.DataFrame, template: Dict[str, complex], theta_df: pd.DataFrame, rotation_sign: float, rotation_slope_deg_per_deg: float, intercept: float, reference_angle_deg: float, model_name: str) -> pd.DataFrame:
    theta_map = {str(row["scan_key"]): float(row["theta_deg"]) for _, row in theta_df.iterrows()}
    rows: List[dict] = []
    for _, row in points_df.iterrows():
        beam_id = str(row["beam_id"])
        z0 = template.get(beam_id, 0.0 + 0.0j)
        theta_deg = theta_map[str(row["scan_key"])]
        pred_x, pred_y = rotate_offset_arcsec(z0.real, z0.imag, theta_deg)
        dx = float(row["x_rel_arcsec"]) - float(np.asarray(pred_x))
        dy = float(row["y_rel_arcsec"]) - float(np.asarray(pred_y))
        resid = float(math.hypot(dx, dy))
        rows.append({
            "model_name": model_name,
            "beam_id": beam_id,
            "stream_name": row.get("stream_name", None),
            "scan_key": str(row["scan_key"]),
            "run_id": row["run_id"],
            "scan_id": row["scan_id"],
            "rep_el_deg": row["rep_el_deg"],
            "x_arcsec": row["x_arcsec"],
            "y_arcsec": row["y_arcsec"],
            "x_rel_arcsec": row["x_rel_arcsec"],
            "y_rel_arcsec": row["y_rel_arcsec"],
            "pred_x_arcsec": float(np.asarray(pred_x)),
            "pred_y_arcsec": float(np.asarray(pred_y)),
            "resid_arcsec": resid,
            "reference_angle_deg": float(reference_angle_deg),
            "rotation_sign": float(rotation_sign),
            "rotation_slope_deg_per_deg": float(rotation_slope_deg_per_deg),
            "dewar_angle_deg": float(intercept),
        })
    return pd.DataFrame(rows)



def _sigma_clip_mask(residual_df: pd.DataFrame, sigma_clip: Optional[float]) -> np.ndarray:
    if sigma_clip is None or not math.isfinite(float(sigma_clip)) or float(sigma_clip) <= 0:
        return np.ones(len(residual_df), dtype=bool)
    resid = _safe_float_array(residual_df["resid_arcsec"])
    med = float(np.nanmedian(resid)) if np.any(np.isfinite(resid)) else 0.0
    mad = float(np.nanmedian(np.abs(resid - med))) if np.any(np.isfinite(resid)) else 0.0
    scale = 1.4826 * mad
    if not math.isfinite(scale) or scale <= 0:
        return np.ones(len(residual_df), dtype=bool)
    thresh = med + float(sigma_clip) * scale
    return np.asarray(resid <= thresh, dtype=bool)



def fit_beam_model(
    df: pd.DataFrame,
    model_name: str,
    center_beam_id: Optional[str] = None,
    *,
    reference_angle_deg: Optional[float] = None,
    sigma_clip: Optional[float] = 4.5,
    clip_iters: int = 2,
    min_points_per_beam: int = 2,
    min_scans_per_beam: int = 2,
) -> BeamFitModelResult:
    points_df, shifts_df = _build_relative_points(df, model_name=model_name, center_beam_id=center_beam_id)
    if points_df.empty:
        raise ValueError(f"no usable rows for model '{model_name}'")

    if reference_angle_deg is None:
        reference_angle_deg = float(np.nanmedian(_safe_float_array(points_df["rep_el_deg"])))
    reference_angle_deg = float(reference_angle_deg)

    current_points = points_df.copy()
    rejected_frames: List[pd.DataFrame] = []
    template: Dict[str, complex] = {}
    theta_df = pd.DataFrame()
    slope = 0.0
    intercept = 0.0
    rotation_sign = 0.0
    residual_df = pd.DataFrame()

    for _ in range(builtins.max(int(clip_iters), 1)):
        beam_counts = current_points.groupby("beam_id").size().rename("fit_n")
        beam_scan_counts = current_points.groupby("beam_id")["scan_key"].nunique().rename("n_scan")
        valid_beams = set(
            beam_counts.loc[beam_counts >= int(min_points_per_beam)].index.astype(str)
        ) & set(
            beam_scan_counts.loc[beam_scan_counts >= int(min_scans_per_beam)].index.astype(str)
        )
        current_points = current_points.loc[current_points["beam_id"].astype(str).isin(valid_beams)].copy()
        if current_points.empty:
            raise ValueError(f"no beams remain for model '{model_name}' after minimum-count filtering")

        template, theta_df, slope, intercept = _iterative_rigid_rotation_fit(current_points, reference_angle_deg=reference_angle_deg)
        rotation_sign, intercept = _quantize_rotation_sign_and_intercept(
            _safe_float_array(theta_df["rep_el_deg"]),
            _safe_float_array(theta_df["theta_deg"]),
            reference_angle_deg=float(reference_angle_deg),
            slope=float(slope),
            intercept=float(intercept),
        )
        residual_df = _residual_table(current_points, template, theta_df, rotation_sign, slope, intercept, reference_angle_deg, model_name)
        keep_mask = _sigma_clip_mask(residual_df, sigma_clip)
        if bool(np.all(keep_mask)):
            break
        rejected = residual_df.loc[~keep_mask].copy()
        if not rejected.empty:
            rejected_frames.append(rejected)
        keep_pairs = set(zip(residual_df.loc[keep_mask, "scan_key"].astype(str), residual_df.loc[keep_mask, "beam_id"].astype(str)))
        current_points = current_points.loc[
            [
                (str(r.scan_key), str(r.beam_id)) in keep_pairs
                for r in current_points[["scan_key", "beam_id"]].itertuples(index=False)
            ]
        ].copy()
        if current_points.empty:
            raise ValueError(f"all rows were rejected for model '{model_name}' during sigma clipping")

    if model_name == "center_beam" and center_beam_id is not None and str(center_beam_id) not in template:
        template[str(center_beam_id)] = 0.0 + 0.0j

    theta_map = {str(row["scan_key"]): float(row["theta_deg"]) for _, row in theta_df.iterrows()}
    shifts_df = shifts_df.merge(theta_df, on=["scan_key", "rep_el_deg"], how="left")
    beam_rows: List[BeamFitRow] = []
    for beam_id in sorted(current_points["beam_id"].astype(str).unique()):
        subset = residual_df.loc[residual_df["beam_id"].astype(str) == beam_id].copy()
        z0 = template.get(str(beam_id), 0.0 + 0.0j)
        if model_name == "center_beam" and center_beam_id is not None and str(beam_id) == str(center_beam_id):
            z0 = 0.0 + 0.0j
        resids = _safe_float_array(subset["resid_arcsec"])
        beam_rows.append(
            BeamFitRow(
                model_name=model_name,
                beam_id=str(beam_id),
                fit_n=int(len(subset)),
                n_scan=int(subset["scan_key"].nunique()),
                n_run=int(subset["run_id"].nunique()),
                az_offset_arcsec=float(z0.real),
                el_offset_arcsec=float(z0.imag),
                radius_arcsec=float(abs(z0)),
                phase_deg=float(np.rad2deg(np.angle(z0))) if abs(z0) > 0 else 0.0,
                rotation_mode="elevation",
                reference_angle_deg=float(reference_angle_deg),
                rotation_sign=float(rotation_sign),
                rotation_slope_deg_per_deg=float(slope),
                dewar_angle_deg=float(intercept),
                rms_resid_arcsec=float(np.sqrt(np.mean(np.square(resids)))) if resids.size else float("nan"),
                max_resid_arcsec=float(np.nanmax(resids)) if resids.size else float("nan"),
            )
        )

    selected_df = residual_df.copy()
    rejected_df = pd.concat(rejected_frames, axis=0, ignore_index=True) if rejected_frames else pd.DataFrame(columns=residual_df.columns)
    return BeamFitModelResult(
        model_name=model_name,
        reference_angle_deg=float(reference_angle_deg),
        rotation_sign=float(rotation_sign),
        rotation_slope_deg_per_deg=float(slope),
        dewar_angle_deg=float(intercept),
        beam_rows=beam_rows,
        shifts_df=shifts_df,
        usable_df=residual_df,
        selected_df=selected_df,
        rejected_df=rejected_df,
    )



def beam_rows_to_dataframe(rows: Sequence[BeamFitRow]) -> pd.DataFrame:
    return pd.DataFrame([
        {
            "model_name": r.model_name,
            "beam_id": r.beam_id,
            "fit_n": r.fit_n,
            "n_scan": r.n_scan,
            "n_run": r.n_run,
            "az_offset_arcsec": r.az_offset_arcsec,
            "el_offset_arcsec": r.el_offset_arcsec,
            "radius_arcsec": r.radius_arcsec,
            "phase_deg": r.phase_deg,
            "rotation_mode": r.rotation_mode,
            "reference_angle_deg": r.reference_angle_deg,
            "rotation_sign": r.rotation_sign,
            "rotation_slope_deg_per_deg": r.rotation_slope_deg_per_deg,
            "dewar_angle_deg": r.dewar_angle_deg,
            "rms_resid_arcsec": r.rms_resid_arcsec,
            "max_resid_arcsec": r.max_resid_arcsec,
        }
        for r in rows
    ])
