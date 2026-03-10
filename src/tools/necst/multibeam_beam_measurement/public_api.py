from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from .config_io import validate_spectrometer_config
from .sunscan_config import SunScanAnalysisConfig
from .sunscan_core import analyze_single_stream
from .sunscan_report import write_singlebeam_outputs
from .sunscan_extract_multibeam import run_extract
from .sunscan_fit_multibeam import run_fit
from .synthetic_multibeam import make_pseudo_multibeam



def run_singlebeam(config: SunScanAnalysisConfig):
    result = analyze_single_stream(config)
    return write_singlebeam_outputs(result)



def run_multibeam_extract(rawdata_path: Path, spectrometer_config: Path, base_config: SunScanAnalysisConfig, *, outdir: Path, run_id: Optional[str] = None, stream_names: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    return run_extract(rawdata_path=rawdata_path, spectrometer_config=spectrometer_config, base_config=base_config, outdir=outdir, run_id=run_id, stream_names=stream_names)



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



def run_pseudo_multibeam(singlebeam_summary_paths: Sequence[Path], spectrometer_config: Path, *, outdir: Path, stream_names: Optional[Sequence[str]] = None, noise_arcsec: float = 0.0, seed: int = 0, tag: Optional[str] = None) -> Dict[str, Path]:
    return make_pseudo_multibeam(
        singlebeam_summary_paths=singlebeam_summary_paths,
        spectrometer_config=spectrometer_config,
        outdir=outdir,
        stream_names=stream_names,
        noise_arcsec=noise_arcsec,
        seed=seed,
        tag=tag,
    )



def check_spectrometer_config(spectrometer_config: Path, *, stream_names: Optional[Sequence[str]] = None):
    return validate_spectrometer_config(spectrometer_config, explicit_stream_names=stream_names)
