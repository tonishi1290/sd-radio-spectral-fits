from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .beam_rotation_shared import apply_beam_offset_arcsec
from .config_io import load_spectrometer_config, restfreq_hz_for_stream, stream_table_from_config



def _source_tag(path: Path, df: pd.DataFrame) -> str:
    if "tag" in df.columns and len(df) > 0:
        v = str(df["tag"].iloc[0]).strip()
        if v:
            return v
    return path.stem



def _source_run_id(path: Path, df: pd.DataFrame) -> str:
    if "run_id" in df.columns and len(df) > 0:
        v = str(df["run_id"].iloc[0]).strip()
        if v:
            return v
    return path.stem



def _ensure_xy(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    if "x_arcsec" not in work.columns:
        work["x_arcsec"] = pd.to_numeric(work.get("center_az_deg"), errors="coerce") * 3600.0
    if "y_arcsec" not in work.columns:
        work["y_arcsec"] = pd.to_numeric(work.get("center_el_deg"), errors="coerce") * 3600.0
    return work



def _expanded_source_rows(src: pd.DataFrame, rep_el_degs: Optional[Sequence[float]]) -> pd.DataFrame:
    work = src.copy()
    if not rep_el_degs:
        if "source_scan_id" not in work.columns:
            work["source_scan_id"] = work.get("scan_id", pd.Series(np.arange(len(work), dtype=int)))
        return work

    rep_el_vals = [float(v) for v in rep_el_degs]
    if len(rep_el_vals) == 0:
        if "source_scan_id" not in work.columns:
            work["source_scan_id"] = work.get("scan_id", pd.Series(np.arange(len(work), dtype=int)))
        return work

    max_scan = 0
    if "scan_id" in work.columns:
        try:
            max_scan = int(np.nanmax(pd.to_numeric(work["scan_id"], errors="coerce").to_numpy(float)))
            if not np.isfinite(max_scan):
                max_scan = 0
        except Exception:
            max_scan = 0
    scan_stride = max(int(max_scan) + 1, len(work), 1)
    expanded: List[pd.DataFrame] = []
    for j, rep_el in enumerate(rep_el_vals):
        clone = work.copy()
        clone["source_scan_id"] = clone.get("scan_id", pd.Series(np.arange(len(clone), dtype=int)))
        clone["rep_el_deg"] = float(rep_el)
        if "scan_id" in clone.columns:
            base_scan = pd.to_numeric(clone["scan_id"], errors="coerce")
            fill_series = pd.Series(np.arange(len(clone), dtype=int), index=clone.index)
            base_scan = base_scan.where(np.isfinite(base_scan), fill_series).astype(int)
        else:
            base_scan = pd.Series(np.arange(len(clone), dtype=int), index=clone.index)
        clone["scan_id"] = base_scan + j * scan_stride
        expanded.append(clone)
    return pd.concat(expanded, axis=0, ignore_index=True)



def make_pseudo_multibeam(
    singlebeam_summary_paths: Sequence[Path],
    spectrometer_config: Path,
    *,
    outdir: Path,
    stream_names: Optional[Sequence[str]] = None,
    noise_arcsec: float = 0.0,
    seed: int = 0,
    tag: Optional[str] = None,
    rep_el_degs: Optional[Sequence[float]] = None,
) -> Dict[str, Path]:
    cfg = load_spectrometer_config(spectrometer_config)
    streams = list(cfg.get("streams", []) or [])
    if stream_names:
        wanted = {str(s) for s in stream_names}
        streams = [s for s in streams if str(s.name) in wanted]
    if not streams:
        raise ValueError("no streams selected for pseudo multibeam generation")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    rows: List[pd.DataFrame] = []
    manifest_rows: List[Dict[str, Any]] = []

    for path in singlebeam_summary_paths:
        src = _ensure_xy(pd.read_csv(path))
        src = _expanded_source_rows(src, rep_el_degs)
        src_tag = _source_tag(path, src)
        src_run_id = _source_run_id(path, src)
        base_tag = str(tag or src_tag)
        for stream in streams:
            beam_dx, beam_dy = apply_beam_offset_arcsec(pd.to_numeric(src["rep_el_deg"], errors="coerce").to_numpy(float), stream.beam)
            stream_df = src.copy()
            jitter_x = rng.normal(0.0, float(noise_arcsec), size=len(stream_df)) if float(noise_arcsec) > 0 else 0.0
            jitter_y = rng.normal(0.0, float(noise_arcsec), size=len(stream_df)) if float(noise_arcsec) > 0 else 0.0
            stream_df["run_id"] = str(src_run_id)
            stream_df["tag"] = f"{base_tag}__pseudo__{stream.name}"
            stream_df["stream_name"] = str(stream.name)
            stream_df["beam_id"] = str(stream.beam.beam_id)
            stream_df["spectral_name"] = str(getattr(stream, "db_stream_name", None) or stream.name)
            stream_df["db_stream_name"] = str(getattr(stream, "db_stream_name", None) or stream.name)
            stream_df["restfreq_hz"] = restfreq_hz_for_stream(stream)
            stream_df["polariza"] = str(stream.polariza)
            stream_df["fdnum"] = int(stream.fdnum)
            stream_df["ifnum"] = int(stream.ifnum)
            stream_df["plnum"] = int(stream.plnum)
            stream_df["sampler"] = stream.sampler
            # Synthetic summaries should mimic the fitted source-relative centers
            # seen by sunscan. A beam physically offset from the boresight by
            # (+dx, +dy) peaks when the boresight is moved by (-dx, -dy), so the
            # synthetic fitted centers must subtract the converter-style offsets.
            stream_df["x_arcsec"] = pd.to_numeric(stream_df["x_arcsec"], errors="coerce").to_numpy(float) - np.asarray(beam_dx, dtype=float) + jitter_x
            stream_df["y_arcsec"] = pd.to_numeric(stream_df["y_arcsec"], errors="coerce").to_numpy(float) - np.asarray(beam_dy, dtype=float) + jitter_y
            stream_df["center_az_deg"] = stream_df["x_arcsec"] / 3600.0
            stream_df["center_el_deg"] = stream_df["y_arcsec"] / 3600.0
            rows.append(stream_df)
            manifest_rows.append({
                "source_summary_csv": str(path),
                "run_id": str(src_run_id),
                "tag": str(base_tag),
                "stream_name": str(stream.name),
                "beam_id": str(stream.beam.beam_id),
                "restfreq_hz": restfreq_hz_for_stream(stream),
                "rep_el_mode": "override" if rep_el_degs else "source",
                "rep_el_values": ",".join(f"{float(v):.6f}" for v in rep_el_degs) if rep_el_degs else None,
            })

    all_df = pd.concat(rows, axis=0, ignore_index=True) if rows else pd.DataFrame()
    preferred_prefix = [
        "run_id", "tag", "stream_name", "beam_id", "scan_id", "source_scan_id",
        "x_arcsec", "y_arcsec", "restfreq_hz",
        "polariza", "fdnum", "ifnum", "plnum", "sampler",
    ]
    ordered_cols = [c for c in preferred_prefix if c in all_df.columns] + [c for c in all_df.columns if c not in preferred_prefix]
    if ordered_cols:
        all_df = all_df.loc[:, ordered_cols]

    base_name = str(tag or "pseudo_multibeam")
    summary_csv = outdir / f"pseudo_multibeam_summary_{base_name}.csv"
    manifest_csv = outdir / f"pseudo_multibeam_manifest_{base_name}.csv"
    stream_table_csv = outdir / f"pseudo_multibeam_stream_table_{base_name}.csv"
    all_df.to_csv(summary_csv, index=False)
    pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)
    stream_table_from_config(spectrometer_config, cfg).to_csv(stream_table_csv, index=False)
    return {
        "summary_csv": summary_csv,
        "manifest_csv": manifest_csv,
        "stream_table_csv": stream_table_csv,
    }



def add_arguments(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("singlebeam_summary_csv", nargs="+", help="One or more single-beam sun_scan summary CSV files")
    ap.add_argument("--spectrometer-config", required=True, help="converter-compatible spectrometer config TOML")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--stream-name", dest="stream_names", action="append", default=None, help="Select a specific stream name (repeatable)")
    ap.add_argument("--noise-arcsec", type=float, default=0.0, help="Optional Gaussian jitter for dry-run testing")
    ap.add_argument("--seed", type=int, default=0, help="Random seed for optional jitter")
    ap.add_argument("--tag", default=None, help="Override base tag for output file names")
    ap.add_argument("--rep-el-deg", dest="rep_el_degs", action="append", type=float, default=None, help="Override representative EL values for synthetic scans (repeatable)")



def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Generate pseudo multi-beam summary CSVs from single-beam sun_scan summaries.")
    add_arguments(ap)
    args = ap.parse_args(argv)
    outputs = make_pseudo_multibeam(
        singlebeam_summary_paths=[Path(p).expanduser().resolve() for p in args.singlebeam_summary_csv],
        spectrometer_config=Path(args.spectrometer_config).expanduser().resolve(),
        outdir=Path(args.outdir).expanduser().resolve(),
        stream_names=args.stream_names,
        noise_arcsec=float(args.noise_arcsec),
        seed=int(args.seed),
        tag=args.tag,
        rep_el_degs=args.rep_el_degs,
    )
    for key, path in outputs.items():
        print(f"[done] {key}: {path}")


if __name__ == "__main__":
    main()
