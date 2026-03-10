# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map.cube_baseline

Iterative baseline subtraction framework for spectral cubes.

Recommended module layout
-------------------------
- sd_radio_spectral_fits.map.baseline_subtraction:
    Core baseline subtraction utilities (line-free detection, multi-ripple frequency estimation,
    poly + multi-sin batch least squares, FITS wrapper).
- sd_radio_spectral_fits.map.cube_baseline:
    Session + orchestrator layer for iterative operation, QC collection, and a thin CLI pipeline.

Axis convention
---------------
This subpackage standardizes cube arrays to (nchan, ny, nx) internally.
For backward-compatibility, BaselineSession can accept (ny, nx, nchan) too and will convert
based on v_axis length.

Public API
----------
- BaselineSession
- LineFreeConfig, RippleConfig, BaselineConfig
- run_one_iteration
- run_cli_pipeline
"""

from .session import BaselineSession
from .orchestrator import run_one_iteration
from .cli import run_cli_pipeline

from ..baseline_subtraction import LineFreeConfig, RippleConfig, BaselineConfig

__all__ = [
    "BaselineSession",
    "LineFreeConfig",
    "RippleConfig",
    "BaselineConfig",
    "run_one_iteration",
    "run_cli_pipeline",
]
