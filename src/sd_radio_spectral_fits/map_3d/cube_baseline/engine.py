# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map.cube_baseline.engine

Vectorized, chunked baseline fitting engine (poly + multi-sin with fixed frequencies).

This engine delegates the heavy lifting to sd_radio_spectral_fits.map.baseline_subtraction.subtract_baseline_cube
and returns the fitted baseline cube + QC maps.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from sd_radio_spectral_fits.map_3d.baseline_subtraction import subtract_baseline_cube, BaselineConfig


def fit_cube_baseline(
    cube_original: np.ndarray,
    *,
    linefree_mask_1d: np.ndarray,
    baseline_cfg: BaselineConfig,
    ripple_freqs: Optional[Sequence[float]] = None,
    target_mask_2d: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Fit baseline model and return baseline cube + QC stats.

    Returns
    -------
    baseline_cube : (nchan, ny, nx) float32
    stats : dict
        rms_map (ny,nx), flag_map (ny,nx), fit_ok_map (ny,nx), n_linefree (scalar array)
    """
    data = np.asarray(cube_original, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(f"cube_original must be 3D (nchan,ny,nx), got shape={data.shape}")
    nchan, ny, nx = data.shape

    lf = np.asarray(linefree_mask_1d, dtype=bool)
    if lf.shape != (nchan,):
        raise ValueError(f"linefree_mask_1d shape mismatch: {lf.shape} vs ({nchan},)")

    out_cube, rms_map, flag_map = subtract_baseline_cube(
        data,
        linefree_mask=lf,
        bcfg=baseline_cfg,
        ripple_freqs=ripple_freqs,
        return_qc=True,
    )
    baseline_cube = np.zeros_like(data, dtype=np.float32)
    finite_both = np.isfinite(data) & np.isfinite(out_cube)
    baseline_cube[finite_both] = (data[finite_both] - out_cube[finite_both]).astype(np.float32, copy=False)

    rms_map = np.asarray(rms_map, dtype=np.float32)
    flag_map = np.asarray(flag_map, dtype=np.uint8)
    fit_ok = (flag_map == 0) & np.isfinite(rms_map)

    if target_mask_2d is not None:
        tm = np.asarray(target_mask_2d, dtype=bool)
        if tm.shape != (ny, nx):
            raise ValueError(f"target_mask_2d shape mismatch: {tm.shape} vs ({ny},{nx})")
        keep = tm

        baseline_cube2 = np.zeros_like(baseline_cube)
        baseline_cube2[:, keep] = baseline_cube[:, keep]
        baseline_cube = baseline_cube2

        rms2 = np.full((ny, nx), np.nan, dtype=np.float32)
        rms2[keep] = rms_map[keep]
        rms_map = rms2

        flg2 = np.zeros((ny, nx), dtype=np.uint8)
        flg2[keep] = flag_map[keep]
        flag_map = flg2

        ok2 = np.zeros((ny, nx), dtype=bool)
        ok2[keep] = fit_ok[keep]
        fit_ok = ok2

    stats: Dict[str, np.ndarray] = {
        "rms_map": rms_map,
        "flag_map": flag_map,
        "fit_ok_map": fit_ok,
        "n_linefree": np.asarray(int(lf.sum()), dtype=np.int32),
    }
    return baseline_cube.astype(np.float32, copy=False), stats
