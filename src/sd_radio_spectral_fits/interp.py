from __future__ import annotations

import numpy as np
import pandas as pd

def _as_int_ns(tidx: pd.DatetimeIndex) -> np.ndarray:
    # pandas Timestamp to int64 ns
    return tidx.view("i8")

def interp_dataframe_in_time(df: pd.DataFrame, target_times: pd.DatetimeIndex) -> np.ndarray:
    """Interpolate spectra in time.

    - Linear interpolation if bracketing samples exist.
    - If only one side exists, use the nearest sample.
    - Returns ndarray (Ntarget, Nchan).
    """
    if not isinstance(df.index, pd.DatetimeIndex):
        raise ValueError("df must have DatetimeIndex for time interpolation")

    src_t = _as_int_ns(df.index)
    tgt_t = _as_int_ns(target_times)

    A = df.to_numpy(dtype=float)
    nsrc, nchan = A.shape

    out = np.empty((tgt_t.size, nchan), dtype=float)

    # indices where to insert to keep order
    idx = np.searchsorted(src_t, tgt_t, side="left")

    for k, j in enumerate(idx):
        if j <= 0:
            out[k] = A[0]
            continue
        if j >= nsrc:
            out[k] = A[-1]
            continue
        t = tgt_t[k]
        t0 = src_t[j - 1]
        t1 = src_t[j]
        if t1 == t0:
            out[k] = A[j]
            continue
        w = (t - t0) / (t1 - t0)
        # if exact match at j, w may be 1.0; that's fine
        out[k] = (1.0 - w) * A[j - 1] + w * A[j]
    return out
