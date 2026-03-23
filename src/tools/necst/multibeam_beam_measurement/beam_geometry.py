from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import builtins
import math

import numpy as np
import pandas as pd

from .pure_rotation_model import recover_el0_offset_from_xy, rotate_el0_offset


@dataclass
class BeamFitRow:
    model_name: str
    beam_id: str
    fit_n: int
    n_scan: int
    n_run: int
    az_offset_arcsec: float
    el_offset_arcsec: float
    observed_template_x_arcsec: float
    observed_template_y_arcsec: float
    radius_arcsec: float
    phase_deg: float
    rotation_mode: str
    reference_angle_deg: float
    rotation_sign: float
    rotation_slope_deg_per_deg: float
    dewar_angle_deg: float
    rms_resid_arcsec: float
    max_resid_arcsec: float
    offset_x_el0_arcsec: float
    offset_y_el0_arcsec: float


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


def observed_rel_to_physical_beam_offset(x_arcsec: float, y_arcsec: float) -> Tuple[float, float]:
    """Convert sunscan observed relative offsets to boresight->beam offsets."""
    return -float(x_arcsec), -float(y_arcsec)


def physical_beam_offset_to_observed_rel(x_arcsec: float, y_arcsec: float) -> Tuple[float, float]:
    """Convert boresight->beam offsets into sunscan observed relative offsets."""
    return -float(x_arcsec), -float(y_arcsec)


def _aggregate_beam_measurements(df: pd.DataFrame) -> pd.DataFrame:
    records: List[dict] = []
    work = df.copy()
    work["beam_id"] = work["beam_id"].astype(str)
    work["stream_name"] = work.get("stream_name", pd.Series([None] * len(work), index=work.index)).astype(str)
    for (run_id, scan_id, beam_id), grp in work.groupby(["run_id", "scan_id", "beam_id"], sort=True):
        x_vals = _safe_float_array(grp["x_arcsec"])
        y_vals = _safe_float_array(grp["y_arcsec"])
        el_vals = _safe_float_array(grp["rep_el_deg"])
        finite = np.isfinite(x_vals) & np.isfinite(y_vals) & np.isfinite(el_vals)
        if not bool(np.any(finite)):
            continue
        streams = sorted({str(v) for v in grp.loc[finite, "stream_name"].dropna().astype(str).tolist() if str(v) and str(v) != "None"})
        records.append(
            {
                "run_id": str(run_id),
                "scan_id": int(scan_id),
                "beam_id": str(beam_id),
                "stream_names": ",".join(streams),
                "n_stream_rows": int(np.count_nonzero(finite)),
                "rep_el_deg": float(np.nanmedian(el_vals[finite])),
                "x_arcsec": float(np.nanmean(x_vals[finite])),
                "y_arcsec": float(np.nanmean(y_vals[finite])),
            }
        )
    return pd.DataFrame(records)


def _build_relative_points(df: pd.DataFrame, model_name: str, center_beam_id: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    agg = _aggregate_beam_measurements(df)
    if agg.empty:
        return pd.DataFrame(), pd.DataFrame()

    records: List[dict] = []
    shift_records: List[dict] = []
    grouped = agg.groupby(["run_id", "scan_id"], sort=True)
    for (run_id, scan_id), grp in grouped:
        grp = grp.copy()
        grp = grp.loc[
            np.isfinite(_safe_float_array(grp["x_arcsec"]))
            & np.isfinite(_safe_float_array(grp["y_arcsec"]))
            & np.isfinite(_safe_float_array(grp["rep_el_deg"]))
        ]
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
        shift_records.append(
            {
                "model_name": model_name,
                "run_id": str(run_id),
                "scan_id": int(scan_id),
                "scan_key": scan_key,
                "rep_el_deg": rep_el_deg,
                "shift_x_arcsec": x0,
                "shift_y_arcsec": y0,
            }
        )
        for _, row in grp.iterrows():
            x_val = float(row["x_arcsec"])
            y_val = float(row["y_arcsec"])
            records.append(
                {
                    "model_name": model_name,
                    "run_id": str(run_id),
                    "scan_id": int(scan_id),
                    "scan_key": scan_key,
                    "beam_id": str(row["beam_id"]),
                    "stream_names": str(row.get("stream_names", "")),
                    "n_stream_rows": int(row.get("n_stream_rows", 1)),
                    "rep_el_deg": rep_el_deg,
                    "x_rel_arcsec": x_val - x0,
                    "y_rel_arcsec": y_val - y0,
                    "x_arcsec": x_val,
                    "y_arcsec": y_val,
                }
            )
    return pd.DataFrame(records), pd.DataFrame(shift_records)


def _compute_template_and_residuals(
    points_df: pd.DataFrame,
    *,
    rotation_sign: float,
    model_name: str,
    center_beam_id: Optional[str],
) -> Tuple[Dict[str, complex], pd.DataFrame, float]:
    sign = float(rotation_sign)
    obs_x = _safe_float_array(points_df["x_rel_arcsec"])
    obs_y = _safe_float_array(points_df["y_rel_arcsec"])
    phys_x = -obs_x
    phys_y = -obs_y
    inv_x, inv_y, _ = recover_el0_offset_from_xy(
        phys_x,
        phys_y,
        _safe_float_array(points_df["rep_el_deg"]),
        sign,
    )
    work = points_df.copy()
    work["offset_x_el0_est_arcsec"] = inv_x
    work["offset_y_el0_est_arcsec"] = inv_y

    template: Dict[str, complex] = {}
    for beam_id, subset in work.groupby(work["beam_id"].astype(str), sort=True):
        if model_name == "center_beam" and center_beam_id is not None and str(beam_id) == str(center_beam_id):
            template[str(beam_id)] = 0.0 + 0.0j
            continue
        template[str(beam_id)] = complex(
            float(pd.to_numeric(subset["offset_x_el0_est_arcsec"], errors="coerce").mean()),
            float(pd.to_numeric(subset["offset_y_el0_est_arcsec"], errors="coerce").mean()),
        )

    rows: List[dict] = []
    sq = []
    for row in work.itertuples(index=False):
        z0 = template.get(str(row.beam_id), 0.0 + 0.0j)
        pred_phys_x, pred_phys_y, theta_deg = rotate_el0_offset(z0.real, z0.imag, float(row.rep_el_deg), sign)
        pred_phys_x = float(np.asarray(pred_phys_x))
        pred_phys_y = float(np.asarray(pred_phys_y))
        pred_x = -pred_phys_x
        pred_y = -pred_phys_y
        resid_x = float(row.x_rel_arcsec) - pred_x
        resid_y = float(row.y_rel_arcsec) - pred_y
        resid = float(math.hypot(resid_x, resid_y))
        sq.append(resid_x * resid_x + resid_y * resid_y)
        rows.append(
            {
                "model_name": model_name,
                "beam_id": str(row.beam_id),
                "stream_names": getattr(row, "stream_names", ""),
                "scan_key": str(row.scan_key),
                "run_id": row.run_id,
                "scan_id": row.scan_id,
                "rep_el_deg": row.rep_el_deg,
                "x_arcsec": row.x_arcsec,
                "y_arcsec": row.y_arcsec,
                "x_rel_arcsec": row.x_rel_arcsec,
                "y_rel_arcsec": row.y_rel_arcsec,
                "pred_x_arcsec": pred_x,
                "pred_y_arcsec": pred_y,
                "pred_phys_x_arcsec": pred_phys_x,
                "pred_phys_y_arcsec": pred_phys_y,
                "resid_x_arcsec": resid_x,
                "resid_y_arcsec": resid_y,
                "resid_arcsec": resid,
                "offset_x_el0_arcsec": float(z0.real),
                "offset_y_el0_arcsec": float(z0.imag),
                "observed_template_x_arcsec": float(-z0.real),
                "observed_template_y_arcsec": float(-z0.imag),
                "rotation_sign": sign,
                "theta_deg": float(np.asarray(theta_deg)),
                "n_stream_rows": getattr(row, "n_stream_rows", 1),
            }
        )
    residual_df = pd.DataFrame(rows)
    metric = float(np.sqrt(np.mean(np.asarray(sq, dtype=float)))) if sq else float("inf")
    return template, residual_df, metric


def _choose_best_sign(points_df: pd.DataFrame, model_name: str, center_beam_id: Optional[str]) -> Tuple[float, Dict[str, complex], pd.DataFrame]:
    best_sign = 1.0
    best_template: Dict[str, complex] = {}
    best_residual_df = pd.DataFrame()
    best_metric = float("inf")
    for sign in (+1.0, -1.0):
        template, residual_df, metric = _compute_template_and_residuals(
            points_df,
            rotation_sign=sign,
            model_name=model_name,
            center_beam_id=center_beam_id,
        )
        if metric < best_metric:
            best_metric = metric
            best_sign = sign
            best_template = template
            best_residual_df = residual_df
    return best_sign, best_template, best_residual_df


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

    current_points = points_df.copy()
    rejected_frames: List[pd.DataFrame] = []
    residual_df = pd.DataFrame()
    template: Dict[str, complex] = {}
    rotation_sign = 1.0

    for _ in range(builtins.max(int(clip_iters), 1)):
        beam_counts = current_points.groupby("beam_id").size().rename("fit_n")
        beam_scan_counts = current_points.groupby("beam_id")["scan_key"].nunique().rename("n_scan")
        valid_beams = set(beam_counts.loc[beam_counts >= int(min_points_per_beam)].index.astype(str)) & set(
            beam_scan_counts.loc[beam_scan_counts >= int(min_scans_per_beam)].index.astype(str)
        )
        current_points = current_points.loc[current_points["beam_id"].astype(str).isin(valid_beams)].copy()
        if current_points.empty:
            raise ValueError(f"no beams remain for model '{model_name}' after minimum-count filtering")

        rotation_sign, template, residual_df = _choose_best_sign(current_points, model_name=model_name, center_beam_id=center_beam_id)
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
                observed_template_x_arcsec=float(-z0.real),
                observed_template_y_arcsec=float(-z0.imag),
                radius_arcsec=float(abs(z0)),
                phase_deg=float(np.rad2deg(np.angle(z0))) if abs(z0) > 0 else 0.0,
                rotation_mode="pure_rotation_v1",
                reference_angle_deg=0.0,
                rotation_sign=float(rotation_sign),
                rotation_slope_deg_per_deg=float("nan"),
                dewar_angle_deg=float("nan"),
                rms_resid_arcsec=float(np.sqrt(np.mean(np.square(resids)))) if resids.size else float("nan"),
                max_resid_arcsec=float(np.nanmax(resids)) if resids.size else float("nan"),
                offset_x_el0_arcsec=float(z0.real),
                offset_y_el0_arcsec=float(z0.imag),
            )
        )

    selected_df = residual_df.copy()
    rejected_df = pd.concat(rejected_frames, axis=0, ignore_index=True) if rejected_frames else pd.DataFrame(columns=residual_df.columns)
    if not shifts_df.empty:
        shifts_df = shifts_df.copy()
        shifts_df["rotation_sign"] = float(rotation_sign)
        shifts_df["theta_deg"] = float("nan")
    return BeamFitModelResult(
        model_name=model_name,
        reference_angle_deg=0.0,
        rotation_sign=float(rotation_sign),
        rotation_slope_deg_per_deg=float("nan"),
        dewar_angle_deg=float("nan"),
        beam_rows=beam_rows,
        shifts_df=shifts_df,
        usable_df=residual_df,
        selected_df=selected_df,
        rejected_df=rejected_df,
    )


def beam_rows_to_dataframe(rows: Sequence[BeamFitRow]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model_name": r.model_name,
                "beam_id": r.beam_id,
                "fit_n": r.fit_n,
                "n_scan": r.n_scan,
                "n_run": r.n_run,
                "az_offset_arcsec": r.az_offset_arcsec,
                "el_offset_arcsec": r.el_offset_arcsec,
                "observed_template_x_arcsec": r.observed_template_x_arcsec,
                "observed_template_y_arcsec": r.observed_template_y_arcsec,
                "radius_arcsec": r.radius_arcsec,
                "phase_deg": r.phase_deg,
                "rotation_mode": r.rotation_mode,
                "reference_angle_deg": r.reference_angle_deg,
                "rotation_sign": r.rotation_sign,
                "rotation_slope_deg_per_deg": r.rotation_slope_deg_per_deg,
                "dewar_angle_deg": r.dewar_angle_deg,
                "rms_resid_arcsec": r.rms_resid_arcsec,
                "max_resid_arcsec": r.max_resid_arcsec,
                "offset_x_el0_arcsec": r.offset_x_el0_arcsec,
                "offset_y_el0_arcsec": r.offset_y_el0_arcsec,
            }
            for r in rows
        ]
    )
