from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd

from . import sunscan_legacy_compat as compat
from .sunscan_config import SunScanAnalysisConfig


CANONICAL_SUMMARY_COLUMNS = [
    "scan_id",
    "rep_az_deg", "rep_el_deg",
    "center_az_deg", "center_el_deg",
    "hpbw_az_arcsec", "hpbw_el_arcsec",
    "sun_az_deg", "sun_el_deg",
    "speed_az_arcsec_s", "speed_el_arcsec_s",
    "fit_ok_az", "fit_ok_el",
    "az_error", "el_error",
    "n_az", "n_az_used", "n_el", "n_el_used",
    "az_left_amp", "az_left_center_deg", "az_left_sigma_deg",
    "az_right_amp", "az_right_center_deg", "az_right_sigma_deg",
    "el_low_amp", "el_low_center_deg", "el_low_sigma_deg",
    "el_high_amp", "el_high_center_deg", "el_high_sigma_deg",
    "ripple_applied_az", "ripple_applied_el",
    "az_track_t_unix", "el_track_t_unix",
    "az_track_az_deg", "az_track_el_deg", "az_track_main_offset_deg", "az_track_cross_offset_deg",
    "el_track_az_deg", "el_track_el_deg", "el_track_main_offset_deg", "el_track_cross_offset_deg",
    "data_tag", "y_axis",
    "azel_source", "altaz_apply",
    "encoder_shift_sec", "encoder_vavg_sec",
    "chopper_wheel",
    "ripple_remove", "ripple_preset",
    "edge_fit_win_deg", "edge_fit_threshold", "hpbw_init_arcsec",
    "trim_scan", "profile_xlim_deg",
]


@dataclass
class SingleBeamReportPaths:
    summary_csv: Optional[Path] = None
    derivative_pngs: List[Path] = field(default_factory=list)
    summary_text_png: Optional[Path] = None
    debug_pngs: List[Path] = field(default_factory=list)


@dataclass
class SingleBeamAnalysisResult:
    config: SunScanAnalysisConfig
    df: pd.DataFrame
    az_scans: Dict[int, pd.DataFrame]
    el_scans: Dict[int, pd.DataFrame]
    ycol: str
    ytitle: str
    trim_params: Dict[str, object]
    ripple_policy: Dict[str, object]
    az_deriv: Dict[int, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    el_deriv: Dict[int, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    az_fit: Dict[int, Tuple[Optional[np.ndarray], Optional[np.ndarray]]] = field(default_factory=dict)
    el_fit: Dict[int, Tuple[Optional[np.ndarray], Optional[np.ndarray]]] = field(default_factory=dict)
    results: Dict[int, dict] = field(default_factory=dict)
    scan_ids: List[int] = field(default_factory=list)
    az_deriv_raw: Dict[int, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    el_deriv_raw: Dict[int, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    az_profile_raw: Dict[int, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    el_profile_raw: Dict[int, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    az_profile: Dict[int, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    el_profile: Dict[int, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    report_paths: SingleBeamReportPaths = field(default_factory=SingleBeamReportPaths)


def results_to_summary_dataframe(scan_ids: List[int], results: Dict[int, dict]) -> pd.DataFrame:
    rows: List[dict] = []
    for sid in sorted([int(s) for s in scan_ids]):
        row = dict(results.get(sid, {}))
        row = {k: v for k, v in row.items() if not str(k).startswith("__")}
        row["scan_id"] = int(sid)
        rows.append(row)
    df = pd.DataFrame(rows)
    cols = [c for c in CANONICAL_SUMMARY_COLUMNS if c in df.columns] + [c for c in df.columns if c not in CANONICAL_SUMMARY_COLUMNS]
    if cols:
        df = df.loc[:, cols]
    return df


def build_summary_lines(tag: str, ytitle: str, config: SunScanAnalysisConfig, scan_ids: List[int], results: Dict[int, dict]) -> List[str]:
    lines: List[str] = []
    lines.append("=" * 112)
    lines.append(f"Sun scan limb-fit summary: {tag}")
    lines.append(
        f"  y-axis: {ytitle}  ripple_remove={bool(config.ripple.enabled)}  "
        f"ripple_preset={str(config.ripple.preset)}  chopper_wheel={bool(config.calibration.chopper_wheel)}"
    )
    lines.append(
        f"  edge_fit_win={float(config.edge_fit.fit_win_deg):.3f} deg  "
        f"threshold={float(config.edge_fit.fit_threshold):.2f}  "
        f"HPBW_init={float(config.edge_fit.hpbw_init_arcsec):.1f} arcsec"
    )
    lines.append("-" * 112)
    lines.append("per-scan results (rep Az/El = mean of the AZ/EL tracking points nearest zero main-axis offset):")
    lines.append("  scan   rep_Az   rep_El | center_az  HPBW_az  v_az | center_el  HPBW_el  v_el | fitAZ fitEL")
    lines.append("  ----  -------  ------- | ---------  -------  ---- | ---------  -------  ---- | ----- -----")
    for sid in sorted([int(s) for s in scan_ids]):
        rr = results.get(sid, {})
        caz = rr.get("center_az_deg", float("nan")) * 3600.0
        haz = rr.get("hpbw_az_arcsec", float("nan"))
        vaz = rr.get("speed_az_arcsec_s", float("nan"))
        cel = rr.get("center_el_deg", float("nan")) * 3600.0
        hel = rr.get("hpbw_el_arcsec", float("nan"))
        vel = rr.get("speed_el_arcsec_s", float("nan"))
        raz = rr.get("rep_az_deg", float("nan"))
        rel = rr.get("rep_el_deg", float("nan"))
        faz = "Y" if rr.get("fit_ok_az", False) else "N"
        fel = "Y" if rr.get("fit_ok_el", False) else "N"
        lines.append(
            f"  {sid:4d}  {raz:7.3f}  {rel:7.3f} | {caz:9.1f}  {haz:7.1f}  {vaz:4.0f} | "
            f"{cel:9.1f}  {hel:7.1f}  {vel:4.0f} |   {faz}     {fel}"
        )

    def _nanmedian(key: str) -> float:
        vals = [results.get(s, {}).get(key, float("nan")) for s in scan_ids]
        arr = np.asarray(vals, dtype=float)
        return float(np.nanmedian(arr)) if np.any(np.isfinite(arr)) else float("nan")

    center_az_med = _nanmedian("center_az_deg") * 3600.0
    center_el_med = _nanmedian("center_el_deg") * 3600.0
    hpbw_az_med = _nanmedian("hpbw_az_arcsec")
    hpbw_el_med = _nanmedian("hpbw_el_arcsec")
    rep_az_med = _nanmedian("rep_az_deg")
    rep_el_med = _nanmedian("rep_el_deg")

    lines.append("-" * 112)
    lines.append("median over scans:")
    lines.append(f"  rep_Az = {rep_az_med:7.3f} deg   rep_El = {rep_el_med:7.3f} deg")
    lines.append(f"  center_az = {center_az_med:7.1f} arcsec   HPBW_az = {hpbw_az_med:6.1f} arcsec")
    lines.append(f"  center_el = {center_el_med:7.1f} arcsec   HPBW_el = {hpbw_el_med:6.1f} arcsec")
    lines.append("=" * 112)
    return lines


def write_singlebeam_outputs(result: SingleBeamAnalysisResult) -> SingleBeamAnalysisResult:
    outdir = Path(result.config.report.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tag = result.config.resolved_tag()
    paths = result.report_paths

    if result.config.edge_fit.enabled and result.results:
        lines = build_summary_lines(tag, result.ytitle, result.config, result.scan_ids, result.results)
        csv_path = outdir / f"sun_scan_summary_{tag}.csv"
        compat.write_scan_summary_csv(csv_path, result.scan_ids, result.results)
        paths.summary_csv = csv_path

        if result.config.report.debug_plot:
            summary_png = outdir / f"summary_text_{tag}.png"
            compat.save_text_summary("\n".join(lines), summary_png)
            paths.summary_text_png = summary_png

        pngs = compat.plot_derivative_fits_paginated(
            outdir,
            tag=tag,
            ycol=result.ycol,
            scan_ids=result.scan_ids,
            az_profile_raw=result.az_profile_raw,
            el_profile_raw=result.el_profile_raw,
            az_profile=result.az_profile,
            el_profile=result.el_profile,
            az_deriv_raw=result.az_deriv_raw,
            el_deriv_raw=result.el_deriv_raw,
            az_deriv=result.az_deriv,
            el_deriv=result.el_deriv,
            az_fit=result.az_fit,
            el_fit=result.el_fit,
            results=result.results,
            xlim_deg=float(result.config.profile.profile_xlim_deg),
            scans_per_page=int(result.config.report.edge_fit_plot_max_scans),
        )
        paths.derivative_pngs = list(pngs)

    result.report_paths = paths
    return result
