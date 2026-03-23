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
from .basketweave import solve_basket_weave_offsets, apply_basket_weave_correction, basket_weave_inplace
from .otf_scan_region import identify_otf_scan_regions, OTFScanRegionResult
# PS (Position Switch)
from .ps_gridder import run_ps_mapping_pipeline, grid_ps

# ==========================================
# 3. 3D Cube Analysis & Moment Generation
# ==========================================
from .cube_analysis import (
    make_3d_mask_for_existing_fits,
    estimate_robust_rms,
    generate_cube_mask,
    append_analysis_hdus_to_fits
)

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
    "identify_otf_scan_regions",
    "OTFScanRegionResult",

    # Cube Analysis
    "estimate_robust_rms",
    "generate_cube_mask",
    "append_analysis_hdus_to_fits",
    "make_3d_mask_for_existing_fits",
    #
]
