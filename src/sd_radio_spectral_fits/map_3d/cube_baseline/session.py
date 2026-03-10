# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map.cube_baseline.session

BaselineSession holds the state of iterative baseline subtraction.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np


def _standardize_cube_axis(
    cube: np.ndarray,
    v_axis: Optional[np.ndarray] = None,
    axis_order_hint: Optional[str] = None,
) -> Tuple[np.ndarray, Optional[np.ndarray], str]:
    """
    Standardize cube axis order to (nchan, ny, nx).

    Accepts:
    - (nchan, ny, nx)  (preferred)
    - (ny, nchan, nx)  (legacy spectral axis in the middle)
    - (ny, nx, nchan)  (legacy spectral axis last)

    Uses v_axis length (if provided) to disambiguate.

    Returns
    -------
    cube_std, v_axis_std, axis_order_in
    """
    arr = np.asarray(cube)
    if arr.ndim != 3:
        raise ValueError(f"cube must be 3D, got shape={arr.shape}")

    if axis_order_hint is not None:
        hint = str(axis_order_hint)
        if hint == "v_y_x":
            return arr, v_axis, hint
        if hint == "y_v_x":
            return np.transpose(arr, (1, 0, 2)), v_axis, hint
        if hint == "y_x_v":
            return np.transpose(arr, (2, 0, 1)), v_axis, hint
        raise ValueError(f"Unknown axis_order_hint: {axis_order_hint}")

    if v_axis is not None:
        nchan = int(np.asarray(v_axis).size)
        if arr.shape[0] == nchan:
            return arr, v_axis, "v_y_x"
        if arr.shape[1] == nchan:
            return np.transpose(arr, (1, 0, 2)), v_axis, "y_v_x"
        if arr.shape[2] == nchan:
            return np.transpose(arr, (2, 0, 1)), v_axis, "y_x_v"
        raise ValueError(
            f"Cannot match v_axis length={nchan} to cube shape={arr.shape}. "
            "Expected one cube axis to equal nchan."
        )

    return arr, v_axis, "unknown"


@dataclass
class IterationRecord:
    action: str
    params: Dict[str, Any]


class BaselineSession:
    """
    State container for iterative baseline subtraction.

    Parameters
    ----------
    cube_original : np.ndarray
        3D cube. Either (nchan, ny, nx) or (ny, nx, nchan).
    v_axis : np.ndarray
        1D axis in km/s (recommended) or any consistent unit. Length must equal nchan.
    """

    def __init__(self, cube_original: np.ndarray, v_axis: np.ndarray, *, axis_order_hint: Optional[str] = None):
        cube_std, v_std, axis_order = _standardize_cube_axis(cube_original, v_axis=v_axis, axis_order_hint=axis_order_hint)
        if v_std is None:
            raise ValueError("v_axis must be provided")

        self.axis_order_in = axis_order

        self.cube_original = np.asarray(cube_std, dtype=np.float32)
        try:
            self.cube_original.flags.writeable = False
        except Exception:
            pass

        self.v_axis = np.asarray(v_std, dtype=np.float64)

        nchan, ny, nx = self.cube_original.shape
        self._shape = (int(nchan), int(ny), int(nx))

        # baseline model and cache
        self.baseline_cube = np.zeros(self._shape, dtype=np.float32)
        self._cube_work_cache: Optional[np.ndarray] = None

        # prior state (loaded from existing FITS or previous completed run)
        self.linefree_mask_1d_prior: Optional[np.ndarray] = None
        self.ripple_freqs_prior: Optional[np.ndarray] = None
        self.prior_source: Optional[str] = None

        # current state
        self.linefree_mask_1d_current: Optional[np.ndarray] = None  # bool (nchan,)
        self.ripple_freqs_current: Optional[np.ndarray] = None      # float (nfreq,)
        self.target_mask_2d: Optional[np.ndarray] = None            # bool (ny, nx)

        # QC maps / stats from last iteration
        self.fit_stats: Dict[str, np.ndarray] = {}

        # history of actions
        self.edit_history: List[IterationRecord] = []

    @property
    def shape(self) -> Tuple[int, int, int]:
        return self._shape

    @property
    def nchan(self) -> int:
        return int(self._shape[0])

    @property
    def ny(self) -> int:
        return int(self._shape[1])

    @property
    def nx(self) -> int:
        return int(self._shape[2])

    def get_full_cube_work(self) -> np.ndarray:
        """
        Return corrected cube (cube_original - baseline_cube), cached.
        """
        if self._cube_work_cache is None:
            self._cube_work_cache = self.cube_original - self.baseline_cube
        return self._cube_work_cache

    def set_prior(
        self,
        *,
        linefree_mask_1d: Optional[np.ndarray] = None,
        ripple_freqs: Optional[Sequence[float]] = None,
        source: Optional[str] = None,
    ) -> None:
        """Store prior baseline products loaded from an external FITS or previous run."""
        if linefree_mask_1d is not None:
            lf = np.asarray(linefree_mask_1d, dtype=bool)
            if lf.shape != (self.nchan,):
                raise ValueError(f"linefree_mask_1d prior shape mismatch: {lf.shape} vs {(self.nchan,)}")
            self.linefree_mask_1d_prior = lf

        if ripple_freqs is not None:
            self.ripple_freqs_prior = np.asarray(list(ripple_freqs), dtype=float)

        if source is not None:
            self.prior_source = str(source)

    def promote_current_to_prior(self, *, source: Optional[str] = None) -> None:
        """Promote the current baseline solution to prior products for future reuse."""
        if self.linefree_mask_1d_current is not None:
            self.linefree_mask_1d_prior = np.asarray(self.linefree_mask_1d_current, dtype=bool).copy()
        if self.ripple_freqs_current is not None:
            self.ripple_freqs_prior = np.asarray(self.ripple_freqs_current, dtype=float).copy()
        if source is not None:
            self.prior_source = str(source)

    def update_baseline(
        self,
        baseline_cube_new: np.ndarray,
        *,
        target_mask_2d: Optional[np.ndarray],
        stats: Dict[str, np.ndarray],
        params: Dict[str, Any],
        linefree_mask_1d: Optional[np.ndarray] = None,
        ripple_freqs: Optional[Sequence[float]] = None,
    ) -> None:
        """
        Update baseline model and QC maps in-place.

        baseline_cube_new must be (nchan, ny, nx).
        Only pixels in target_mask_2d are updated if provided; else full cube.
        """
        base = np.asarray(baseline_cube_new, dtype=np.float32)
        if base.shape != self._shape:
            raise ValueError(f"baseline_cube_new shape mismatch: {base.shape} vs {self._shape}")

        if target_mask_2d is None:
            target_mask_2d = np.ones((self.ny, self.nx), dtype=bool)
        else:
            target_mask_2d = np.asarray(target_mask_2d, dtype=bool)
            if target_mask_2d.shape != (self.ny, self.nx):
                raise ValueError(f"target_mask_2d shape mismatch: {target_mask_2d.shape} vs {(self.ny, self.nx)}")

        self.baseline_cube[:, target_mask_2d] = base[:, target_mask_2d]

        for key, new_map in stats.items():
            new_arr = np.asarray(new_map)
            if new_arr.ndim == 2 and new_arr.shape == (self.ny, self.nx):
                if key not in self.fit_stats:
                    fill = np.nan if np.issubdtype(new_arr.dtype, np.floating) else 0
                    self.fit_stats[key] = np.full((self.ny, self.nx), fill, dtype=new_arr.dtype)
                self.fit_stats[key][target_mask_2d] = new_arr[target_mask_2d]
            elif new_arr.ndim == 1 and new_arr.shape == (self.nchan,):
                self.fit_stats[key] = new_arr.copy()
            else:
                self.fit_stats[key] = new_arr.copy()

        if linefree_mask_1d is not None:
            lf = np.asarray(linefree_mask_1d, dtype=bool)
            if lf.shape != (self.nchan,):
                raise ValueError(f"linefree_mask_1d shape mismatch: {lf.shape} vs {(self.nchan,)}")
            self.linefree_mask_1d_current = lf

        if ripple_freqs is not None:
            self.ripple_freqs_current = np.asarray(list(ripple_freqs), dtype=float)

        self.target_mask_2d = target_mask_2d
        self._cube_work_cache = None

        self.edit_history.append(IterationRecord(action="fit_baseline", params=params))
