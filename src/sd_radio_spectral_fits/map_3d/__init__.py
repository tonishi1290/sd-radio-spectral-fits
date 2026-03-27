# src/sd_radio_spectral_fits/map_3d/__init__.py
"""
Gridding, mapping, and 3D cube analysis tools for single-dish radio astronomy data.
"""

# ==========================================
# 1. Configurations & Data Structures
# ==========================================
from .config import MapConfig, GridInput, GridResult
from .ps_gridder import PSMapConfig

# ==========================================
# 2. Mapping Pipelines (OTF & PS)
# ==========================================
# OTF (On-The-Fly)
from .gridder import run_mapping_pipeline, create_grid_input
from .core import grid_otf
from .basketweave import solve_basket_weave_offsets, apply_basket_weave_correction, basket_weave_inplace, basketweave_cubes, basketweave_fits, run_otf_plait_pipeline
from .otf_bundle import OTFBundle
from .otf_bundle_io import read_otf_bundle, write_otf_bundle, gridresult_to_otf_bundle, validate_otf_bundle
from .otf_family_grid import grid_otf_family
from .cube_coadd import coadd_family_cubes
from .mosaic import attach_mosaic_products, attach_mosaic_products_from_mask, mosaic_bundles, mosaic_fits
from .plait_fft import plait_fft_cubes
from .otf_pipeline import run_otf_full_pipeline, run_otf_full_pipeline_multi
from .otf_scan_region import identify_otf_scan_regions, OTFScanRegionResult, format_otf_scan_summary

try:
    from .baseline_subtraction import subtract_baseline_from_bundle, subtract_baseline_from_fits, make_baseline_viewer_bundle
    _baseline_import_error = None
except Exception as _baseline_exc:
    _baseline_import_error = _baseline_exc

    def _raise_baseline_import_error(*args, **kwargs):
        raise ModuleNotFoundError(
            "baseline_subtraction features require optional dependencies such as astropy/scipy. "
            f"Original import error: {_baseline_import_error}"
        )

    subtract_baseline_from_bundle = _raise_baseline_import_error
    subtract_baseline_from_fits = _raise_baseline_import_error
    make_baseline_viewer_bundle = _raise_baseline_import_error

# PS (Position Switch)
from .ps_gridder import run_ps_mapping_pipeline, grid_ps

# ==========================================
# 3. 3D Cube Analysis & Moment Generation
# ==========================================
try:
    from .cube_analysis import (
        make_3d_mask_for_existing_fits,
        estimate_robust_rms,
        generate_cube_mask,
        append_analysis_hdus_to_fits,
    )
    _cube_analysis_import_error = None
except Exception as _cube_analysis_exc:
    _cube_analysis_import_error = _cube_analysis_exc

    def _raise_cube_analysis_import_error(*args, **kwargs):
        raise ModuleNotFoundError(
            "cube_analysis features require optional dependencies such as spectral_cube. "
            f"Original import error: {_cube_analysis_import_error}"
        )

    make_3d_mask_for_existing_fits = _raise_cube_analysis_import_error
    estimate_robust_rms = _raise_cube_analysis_import_error
    generate_cube_mask = _raise_cube_analysis_import_error
    append_analysis_hdus_to_fits = _raise_cube_analysis_import_error

# 明示的に公開するAPIのリスト (IDEの補完候補をクリーンに保つため)
__all__ = [
    # Configs
    "MapConfig",
    "PSMapConfig",
    "GridInput",
    "GridResult",
    
    # Pipelines
    "run_mapping_pipeline",
    "create_grid_input",
    "grid_otf",
    "run_ps_mapping_pipeline",
    "grid_ps",
    "solve_basket_weave_offsets",
    "apply_basket_weave_correction",
    "basket_weave_inplace",
    "basketweave_cubes",
    "basketweave_fits",
    "run_otf_plait_pipeline",
    "OTFBundle",
    "read_otf_bundle",
    "write_otf_bundle",
    "gridresult_to_otf_bundle",
    "validate_otf_bundle",
    "grid_otf_family",
    "coadd_family_cubes",
    "attach_mosaic_products",
    "attach_mosaic_products_from_mask",
    "mosaic_bundles",
    "mosaic_fits",
    "subtract_baseline_from_bundle",
    "subtract_baseline_from_fits",
    "make_baseline_viewer_bundle",
    "plait_fft_cubes",
    "run_otf_full_pipeline",
    "run_otf_full_pipeline_multi",
    "identify_otf_scan_regions",
    "OTFScanRegionResult",
    "format_otf_scan_summary",

    # Cube Analysis
    "estimate_robust_rms",
    "generate_cube_mask",
    "append_analysis_hdus_to_fits",
    "make_3d_mask_for_existing_fits",
    #
]
