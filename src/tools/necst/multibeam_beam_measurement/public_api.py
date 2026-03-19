from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import copy
import re

import pandas as pd

from .config_io import validate_spectrometer_config
from .sunscan_config import SunScanAnalysisConfig
from .sunscan_core import analyze_single_stream
from .sunscan_report import write_singlebeam_outputs
from .sunscan_extract_multibeam import run_extract, run_extract_many
from .sunscan_fit_multibeam import run_fit
from .synthetic_multibeam import make_pseudo_multibeam



def run_singlebeam(config: SunScanAnalysisConfig):
    result = analyze_single_stream(config)
    return write_singlebeam_outputs(result)



def _sanitize_path_tag_component(value: Any) -> str:
    text = str(value).strip()
    if not text:
        return "unnamed"
    safe = re.sub(r"[^0-9A-Za-z._-]+", "_", text)
    safe = safe.strip("._-")
    return safe or "unnamed"



def _default_multi_base_tag(rawdata_path: Path) -> str:
    name = _sanitize_path_tag_component(rawdata_path.name)
    parent = _sanitize_path_tag_component(rawdata_path.parent.name)
    if parent and parent != "unnamed" and parent != name:
        return f"{parent}__{name}"
    return name



def _deduplicate_tags(tags: Sequence[str]) -> List[str]:
    counts: Dict[str, int] = {}
    unique: List[str] = []
    for tag in tags:
        base = _sanitize_path_tag_component(tag)
        n = counts.get(base, 0) + 1
        counts[base] = n
        unique.append(base if n == 1 else f"{base}__{n:02d}")
    return unique



def _prepare_multi_run_entries(rawdata_paths: Sequence[Path]) -> List[Dict[str, Any]]:
    normalized_paths = [Path(p).expanduser().resolve() for p in rawdata_paths]
    name_counts: Dict[str, int] = {}
    for path in normalized_paths:
        key = str(path.name)
        name_counts[key] = name_counts.get(key, 0) + 1
    candidate_tags: List[str] = []
    for path in normalized_paths:
        if name_counts.get(str(path.name), 0) > 1:
            candidate_tags.append(_default_multi_base_tag(path))
        else:
            candidate_tags.append(_sanitize_path_tag_component(path.name))
    unique_tags = _deduplicate_tags(candidate_tags)
    return [
        {
            "rawdata_path": path,
            "base_tag": tag,
            "run_id": tag,
        }
        for path, tag in zip(normalized_paths, unique_tags)
    ]



def _insert_or_replace_front(df: pd.DataFrame, column: str, value: Any) -> None:
    if column in df.columns:
        df.pop(column)
    df.insert(0, column, value)



def run_singlebeam_many(rawdata_paths: Sequence[Path], base_config: SunScanAnalysisConfig, *, outdir: Path, run_id: Optional[str] = None) -> Dict[str, Any]:
    entries = _prepare_multi_run_entries(rawdata_paths)
    resolved_outdir = Path(outdir).expanduser().resolve()
    per_run_root = resolved_outdir / "per_run"
    merged_root = resolved_outdir / "merged"
    per_run_root.mkdir(parents=True, exist_ok=True)
    merged_root.mkdir(parents=True, exist_ok=True)

    summary_frames: List[pd.DataFrame] = []
    run_rows: List[Dict[str, Any]] = []
    results: List[Any] = []

    for entry in entries:
        rawdata_path = Path(entry["rawdata_path"]).expanduser().resolve()
        base_tag = str(entry["base_tag"])
        this_run_id = str(entry["run_id"])
        per_run_outdir = per_run_root / base_tag
        per_run_outdir.mkdir(parents=True, exist_ok=True)

        cfg = copy.deepcopy(base_config)
        cfg.input.rawdata_path = rawdata_path
        cfg.report.outdir = per_run_outdir
        cfg.report.tag = base_tag

        try:
            result = run_singlebeam(cfg)
            results.append(result)
            summary_csv = result.report_paths.summary_csv
            derivative_pngs = list(result.report_paths.derivative_pngs or [])
            debug_png = result.report_paths.summary_text_png
            run_rows.append({
                "run_id": this_run_id,
                "base_tag": base_tag,
                "rawdata_path": str(rawdata_path),
                "per_run_outdir": str(per_run_outdir),
                "status": "success",
                "summary_csv": str(summary_csv) if summary_csv else None,
                "derivative_png_count": len(derivative_pngs),
                "summary_text_png": str(debug_png) if debug_png else None,
                "error": None,
            })
            if summary_csv is not None and Path(summary_csv).exists():
                df = pd.read_csv(summary_csv)
                _insert_or_replace_front(df, "rawdata_path", str(rawdata_path))
                _insert_or_replace_front(df, "base_tag", base_tag)
                _insert_or_replace_front(df, "run_id", this_run_id)
                summary_frames.append(df)
        except Exception as exc:
            run_rows.append({
                "run_id": this_run_id,
                "base_tag": base_tag,
                "rawdata_path": str(rawdata_path),
                "per_run_outdir": str(per_run_outdir),
                "status": "error",
                "summary_csv": None,
                "derivative_png_count": 0,
                "summary_text_png": None,
                "error": f"{type(exc).__name__}: {exc}",
            })
            if not base_config.runtime.continue_on_error:
                raise

    merged_label = _sanitize_path_tag_component(run_id) if run_id else "all"
    merged_summary_csv = None
    if summary_frames:
        merged_df = pd.concat(summary_frames, ignore_index=True, sort=False)
        if "scan_id" in merged_df.columns:
            scan_col = merged_df.pop("scan_id")
            merged_df.insert(3 if len(merged_df.columns) >= 3 else len(merged_df.columns), "scan_id", scan_col)
        merged_summary_csv = merged_root / f"sun_scan_summary_{merged_label}.csv"
        merged_df.to_csv(merged_summary_csv, index=False)
    run_table_df = pd.DataFrame(run_rows)
    run_table_csv = merged_root / f"sun_scan_run_table_{merged_label}.csv"
    run_table_df.to_csv(run_table_csv, index=False)

    return {
        "merged_summary_csv": merged_summary_csv,
        "run_table_csv": run_table_csv,
        "run_table_df": run_table_df,
        "results": results,
        "summary_csvs": [Path(p) for p in run_table_df["summary_csv"].dropna().astype(str).tolist()] if (not run_table_df.empty and "summary_csv" in run_table_df.columns) else [],
    }



def run_multibeam_extract(rawdata_path: Path, spectrometer_config: Path, base_config: SunScanAnalysisConfig, *, outdir: Path, run_id: Optional[str] = None, stream_names: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    return run_extract(rawdata_path=rawdata_path, spectrometer_config=spectrometer_config, base_config=base_config, outdir=outdir, run_id=run_id, stream_names=stream_names)



def run_multibeam_extract_many(rawdata_paths: Sequence[Path], spectrometer_config: Path, base_config: SunScanAnalysisConfig, *, outdir: Path, run_ids: Optional[Sequence[Optional[str]]] = None, stream_names: Optional[Sequence[str]] = None, merged_tag: Optional[str] = None) -> Dict[str, Any]:
    return run_extract_many(
        rawdata_paths=rawdata_paths,
        spectrometer_config=spectrometer_config,
        base_config=base_config,
        outdir=outdir,
        run_ids=run_ids,
        stream_names=stream_names,
        merged_tag=merged_tag,
    )



def run_multibeam_fit(summary_paths: Sequence[Path], spectrometer_config: Path, *, outdir: Path, center_beam_id: Optional[str] = None, stream_names: Optional[Sequence[str]] = None, model: str = "both", reference_angle_deg: Optional[float] = None, sigma_clip: Optional[float] = 4.5, clip_iters: int = 2, min_points_per_beam: int = 2, min_scans_per_beam: int = 2) -> Dict[str, Path]:
    return run_fit(
        summary_paths=summary_paths,
        spectrometer_config=spectrometer_config,
        outdir=outdir,
        center_beam_id=center_beam_id,
        stream_names=stream_names,
        model=model,
        reference_angle_deg=reference_angle_deg,
        sigma_clip=sigma_clip,
        clip_iters=clip_iters,
        min_points_per_beam=min_points_per_beam,
        min_scans_per_beam=min_scans_per_beam,
    )



def run_pseudo_multibeam(singlebeam_summary_paths: Sequence[Path], spectrometer_config: Path, *, outdir: Path, stream_names: Optional[Sequence[str]] = None, noise_arcsec: float = 0.0, seed: int = 0, tag: Optional[str] = None, rep_el_degs: Optional[Sequence[float]] = None) -> Dict[str, Path]:
    return make_pseudo_multibeam(
        singlebeam_summary_paths=singlebeam_summary_paths,
        spectrometer_config=spectrometer_config,
        outdir=outdir,
        stream_names=stream_names,
        noise_arcsec=noise_arcsec,
        seed=seed,
        tag=tag,
        rep_el_degs=rep_el_degs,
    )



def check_spectrometer_config(spectrometer_config: Path, *, stream_names: Optional[Sequence[str]] = None):
    return validate_spectrometer_config(spectrometer_config, explicit_stream_names=stream_names)
