from __future__ import annotations

import builtins
from dataclasses import dataclass
from typing import Dict, Tuple, Optional

import numpy as np

C_KMS = 299792.458

def freq_axis_from_wcs(meta: dict, nchan: int) -> np.ndarray:
    # Compatibility: accept internal keys or FITS WCS keys
    crval = float(meta.get("CRVAL1"))
    cdelt = float(meta.get("CDELT1"))
    crpix = float(meta.get("CRPIX1", 1.0))

    if not np.isfinite(crval) or not np.isfinite(cdelt) or not np.isfinite(crpix):
        raise KeyError("freq_axis_from_wcs: missing/invalid CRVAL1/CDELT1/CRPIX1 (or internal equivalents).")

    i = np.arange(int(nchan), dtype=float)  # 0..nchan-1
    # FITS WCS: world = CRVAL + ((i+1) - CRPIX)*CDELT
    return crval + ((i + 1.0) - crpix) * cdelt
    
def radio_velocity_kms(freq_hz: np.ndarray, rest_hz: float) -> np.ndarray:
    """Radio definition: v = c * (rest - freq) / rest."""
    rest = float(rest_hz)
    return C_KMS * (rest - freq_hz) / rest

def wcs_slice_channels(meta: dict, ch_start: int, ch_stop: int) -> dict:
    """Return a new meta dict updated for channel slicing [ch_start:ch_stop].

    Policy (safe):
      - keep CRPIX1 unchanged
      - shift CRVAL1 by ch_start * CDELT1
    """
    if ch_start < 0 or ch_stop <= ch_start:
        raise ValueError(f"invalid slice: {ch_start}:{ch_stop}")
    meta2 = dict(meta)
    meta2["NAXIS1"] = int(ch_stop - ch_start)
    meta2["CRVAL1"] = float(meta["CRVAL1"]) + float(ch_start) * float(meta["CDELT1"])
    # CRPIX1 unchanged
    return meta2

def channel_slice_from_vrange_union(
    meta: dict,
    v_corr_kms: np.ndarray,
    vmin_kms: float,
    vmax_kms: float,
    rest_hz: Optional[float] = None,
) -> Tuple[int, int]:
    """Determine a single channel slice that covers [vmin, vmax] for ALL dumps.

    This uses:
      v_lsrk(ch, dump) = v_radio(ch) + v_corr_kms[dump]

    Returns (ch_start, ch_stop) with stop exclusive.

    Notes:
      - Requires v_corr_kms per dump (can be constant array).
      - Uses union over dumps, to ensure requested velocity range is included for all times.
    """
    import builtins
    rest = float(meta.get("RESTFREQ") if rest_hz is None else rest_hz)
    nchan = int(meta.get("NAXIS1", 0))
    if nchan <= 0:
        raise KeyError("channel_slice_from_vrange_union: meta lacks nchan (or NCHAN).")
    freq = freq_axis_from_wcs(meta, nchan=nchan)
    v_radio = radio_velocity_kms(freq, rest)  # shape (nchan,)
    v1, v2 = float(vmin_kms), float(vmax_kms)
    if v2 < v1:
        v1, v2 = v2, v1

    starts = []
    stops = []
    for vc in np.asarray(v_corr_kms, float):
        v = v_radio + vc
        m = (v >= v1) & (v <= v2)
        if not np.any(m):
            raise ValueError("Requested velocity window does not overlap spectral axis.")
        idx = np.where(m)[0]
        starts.append(int(idx[0]))
        stops.append(int(idx[-1]) + 1)  # exclusive

    return builtins.min(starts), builtins.max(stops)


def slice_channels(meta: dict, data: np.ndarray, ch_start: int | None, ch_stop: int | None):
    """Slice channels and update frequency WCS safely.

    - ch_start: inclusive, 0-based
    - ch_stop : exclusive, 0-based

    Updates:
      nchan and crval1_hz (CRVAL1_new = CRVAL1 + ch_start * CDELT1).
    """
    import copy
    if ch_start is None and ch_stop is None:
        return meta, data
    nchan = data.shape[1]
    s = 0 if ch_start is None else int(ch_start)
    e = nchan if ch_stop is None else int(ch_stop)
    if not (0 <= s < e <= nchan):
        raise ValueError(f"Invalid channel slice {s}:{e} for nchan={nchan}")
    out = np.asarray(data[:, s:e], float)
    m = copy.deepcopy(meta)
    m["NAXIS1"] = int(e - s)
    m["CRVAL1"] = float(m["CRVAL1"]) + float(s) * float(m["CDELT1"])
    return m, out
