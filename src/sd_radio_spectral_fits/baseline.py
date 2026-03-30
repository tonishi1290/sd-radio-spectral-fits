# src/sd_radio_spectral_fits/baseline.py
from __future__ import annotations

import datetime
from dataclasses import dataclass
from typing import List, Tuple, Optional, Union, Sequence, Dict, Any

import numpy as np
import pandas as pd

from .fitsio import Scantable, read_scantable, write_scantable
# [MODIFIED] Use robust axis helpers from regrid_vlsrk
from .regrid_vlsrk import vlsrk_axis_for_spectrum
from .ranges import parse_windows, window_to_mask
from .utils import FailPolicy, subtract_windows, in_any_windows
from .tempscale import (
    normalize_tempscal, tempscal_array, beameff_array, require_beameff, append_scale_history
)
# [ADDED] Row selector and Endian utilities
from .scantable_utils import _parse_row_selector, _df_to_native_endian

from .restfreq import apply_restfreq_override
from .axis import wcs_slice_channels

# =========================================================
# 0. Data Structures
# =========================================================

@dataclass
class BaselineModel:
    """
    Model parameters for baseline fitting.
    Used by coadd.py and other modules.
    """
    poly_order: int
    v_windows_kms: List[Tuple[float, float]]
    convention: str = "radio"
    frame: str = "LSRK"


BSL_METADATA_COLUMNS = [
    "BSL_DONE",
    "BSL_APPLIED",
    "BSL_STAGE",
    "BSL_POLY",
    "BSL_WINF",
    "BSL_RMS",
    "BSL_STAT",
    "BSL_NUSED",
    "BSL_COEF",
    "BSL_SCALE",
]


def _drop_bsl_metadata_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=[c for c in BSL_METADATA_COLUMNS if c in df.columns], errors="ignore")

# =========================================================
# 1. Low-level Fitting Logic (Core)
# =========================================================

def fit_polynomial_baseline(
    v_kms: np.ndarray,
    y: np.ndarray,
    base_windows: Optional[List[Tuple[float, float]]] = None,
    *,
    # Aliases / Legacy Support
    windows: Optional[List[Tuple[float, float]]] = None, 
    line_windows: Optional[List[Tuple[float, float]]] = None,
    poly_order: int = 1, 
    order: Optional[int] = None,
    # Iterative parameters
    iter_max: int = 0,
    iter_sigma: float = 3.0
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Fit a polynomial baseline to the spectrum, optionally with iterative sigma-clipping.
    """
    # 1. Parameter Normalization
    v = np.asarray(v_kms, float)
    y = np.asarray(y, float)
    deg = int(order if order is not None else poly_order)

    # Resolve Windows
    # base_windows takes precedence, then windows
    wins = base_windows if base_windows is not None else windows
    
    # If absolutely no window is provided, usually implies "use everything" 
    # OR "cannot fit". Here we default to "cannot fit" to be safe, 
    # unless user explicitly passed empty list? 
    # Let's assume None = "No Fit", [] = "Use All" (Legacy behavior varies)
    # -> Safe bet: If None, return zero baseline.
    if wins is None:
        return np.zeros(deg + 1), np.zeros_like(y), {"mask": np.zeros_like(y, dtype=bool), "rms": np.nan, "std": np.nan}

    # Subtract line windows if any
    line_wins = list(line_windows or [])
    use_windows = subtract_windows(wins, line_wins) if line_wins else list(wins)
    
    # 2. Build Initial Mask
    if not use_windows:
        # If explicit empty list or subtraction resulted in empty, 
        # usually means "use full range" in some contexts, but here we treat strictly.
        # Fallback: if no windows defined, use valid data points? 
        # Better: use everything (continuum mode)
        mask_window = np.ones_like(y, dtype=bool)
    else:
        mask_window = in_any_windows(v, use_windows)
        
    mask_valid = np.isfinite(y) & mask_window
    
    # Check points count
    if np.count_nonzero(mask_valid) <= deg:
        return np.zeros(deg + 1), np.zeros_like(y), {"mask": mask_valid, "rms": np.nan, "std": np.nan}

    # 3. Fitting Loop (Iterative)
    current_mask = mask_valid.copy()
    coeffs = np.zeros(deg + 1)
    
    # [追加] フィッティング完了時のベースライン配列を初期化
    baseline = np.zeros_like(y)
    
    # Numerical stability: normalize x to [-1, 1] within the valid range
    # But polyfit/polyval works on raw x. For high order, normalization is recommended.
    # Here we stick to raw x for simplicity and consistency with stored coefficients.
    # (If stability issues arise, we should switch to Chebyshev or normalized domain)

    for i in range(iter_max + 1):
        if np.count_nonzero(current_mask) <= deg:
            break
            
        try:
            # [変更] np.polyfit から np.polynomial.Polynomial.fit へ移行
            p = np.polynomial.Polynomial.fit(v[current_mask], y[current_mask], deg)
            
            # [変更] np.polyval から p() へ移行
            baseline = p(v)
            
            # FITS保存用に、従来の np.polyfit と同じ形式（高次->低次）に変換
            coeffs = p.convert().coef[::-1]
            
        except np.linalg.LinAlgError:
            break
            
        if i == iter_max:
            break
            
        # Clipping
        resid = y - baseline
        # Calc sigma only on current mask
        resid_valid = resid[current_mask]
        if resid_valid.size < 2:
            break
            
        sigma = np.std(resid_valid, ddof=1)
        if sigma <= 0:
            break
            
        # Outlier detection
        # Note: We only reject points that are currently in the mask
        # (Once rejected, they stay rejected? Or do we allow re-inclusion?)
        # Standard sigma clipping usually allows re-inclusion if fit changes, 
        # but here we strictly shrink mask for stability.
        
        abs_dev = np.abs(resid)
        # Check against threshold
        is_outlier = (abs_dev > iter_sigma * sigma)
        
        # Update mask: Keep points that are valid AND not outliers
        new_mask = mask_valid & (~is_outlier)
        
        if np.array_equal(new_mask, current_mask):
            break # Converged
            
        current_mask = new_mask

    # 4. Final Calc
    # [削除] 元の baseline = np.polyval(coeffs, v) は削除
    
    # Calc RMS on final used points
    resid_final = (y - baseline)[current_mask]
    if resid_final.size > 1:
        std_val = float(np.std(resid_final, ddof=1))
        rms_val = float(np.sqrt(np.mean(resid_final**2)))
    else:
        std_val = np.nan
        rms_val = np.nan

    info = {
        "mask": current_mask,
        "rms": rms_val,
        "std": std_val,
        "poly_order": deg,
        "iter_max": iter_max
    }
    
    return coeffs, baseline, info
    

# =========================================================
# 2. High-level Interface (VLA Support)
# =========================================================

def _apply_baseline_to_table(
    df: pd.DataFrame,
    coeffs_list: List[Optional[np.ndarray]],
    rms_list: List[float],
    poly_orders: Union[int, List[int]],
    win_f_list: Union[str, List[str]] = "",
    *,
    applied: bool = True,
    stage: str = "baseline_fit",
    done_list: Optional[List[bool]] = None,
    nused_list: Optional[List[int]] = None,
) -> pd.DataFrame:
    """Add baseline results to the dataframe."""
    n = len(df)
    df = _drop_bsl_metadata_columns(df.copy())

    if done_list is None:
        done_arr = np.ones(n, dtype=bool)
    else:
        done_arr = np.asarray(done_list, dtype=bool)
        if done_arr.size != n:
            raise ValueError(f"done_list size mismatch: {done_arr.size} != {n}")

    if nused_list is None:
        nused_arr = np.zeros(n, dtype=np.int32)
    else:
        nused_arr = np.asarray(nused_list, dtype=np.int32)
        if nused_arr.size != n:
            raise ValueError(f"nused_list size mismatch: {nused_arr.size} != {n}")

    df["BSL_DONE"] = done_arr
    df["BSL_APPLIED"] = (done_arr & bool(applied))
    df["BSL_STAGE"] = [str(stage)] * n
    df["BSL_RMS"] = np.asarray(rms_list, dtype=np.float32)
    df["BSL_STAT"] = ["std"] * n
    df["BSL_NUSED"] = nused_arr

    if applied:
        if isinstance(poly_orders, (list, np.ndarray)):
            poly_arr = np.asarray(poly_orders, dtype=np.int16)
            if poly_arr.size != n:
                raise ValueError(f"poly_orders size mismatch: {poly_arr.size} != {n}")
        else:
            poly_arr = np.full(n, int(poly_orders), dtype=np.int16)
        df["BSL_POLY"] = poly_arr

        if isinstance(win_f_list, list):
            if len(win_f_list) != n:
                raise ValueError(f"win_f_list size mismatch: {len(win_f_list)} != {n}")
            win_list = [str(v) for v in win_f_list]
        else:
            win_list = [str(win_f_list)] * n
        df["BSL_WINF"] = win_list

        coef_out: list[Optional[np.ndarray]] = []
        for c, done in zip(coeffs_list, done_arr):
            if (not done) or c is None:
                coef_out.append(None)
            else:
                coef_out.append(np.asarray(c, dtype=float))
        df["BSL_COEF"] = coef_out

        if "TEMPSCAL" in df.columns:
            df["BSL_SCALE"] = df["TEMPSCAL"]
        else:
            df["BSL_SCALE"] = "TA*"

    return df

def _row_or_meta(row: pd.Series, meta: dict, key: str, default: Any = None) -> Any:
    if key in row.index:
        val = row[key]
        if val is not None and not pd.isna(val):
            return val
    if key in meta:
        val = meta.get(key)
        if val is not None and not pd.isna(val):
            return val
    return default


def _row_only_value(row: pd.Series, key: str, default: Any = None) -> Any:
    """Return a row value only; ignore meta/header fallbacks.

    VELOSYS / VFRAME are row-wise metadata in this package. We intentionally do not
    fall back to global meta/header values because those can only represent a single
    value and are therefore ambiguous for multi-row data.
    """
    if key in row.index:
        val = row[key]
        if val is not None and not pd.isna(val):
            return val
    return default


def _resolve_specsys(row: pd.Series, meta: dict) -> str:
    for key in ("SPECSYS", "SSYSOBS"):
        val = _row_or_meta(row, meta, key, None)
        if val not in (None, ""):
            return str(val).strip().upper()
    return ""


def _resolve_ssysobs(row: pd.Series, meta: dict) -> str:
    for key in ("SSYSOBS", "SPECSYS"):
        val = _row_or_meta(row, meta, key, None)
        if val not in (None, ""):
            return str(val).strip().upper()
    return ""


def _vcorr_scale_to_kms(key: str) -> float:
    ku = str(key or "").strip().upper()
    if ku.endswith("_KMS") or ku == "V_CORR_KMS":
        return 1.0
    # Standard policy: VELOSYS / VFRAME are stored in m/s.
    return 1.0e-3


def _vcorr_candidate_names(preferred_key: str = "VELOSYS") -> list[str]:
    candidates: list[str] = []
    pk = str(preferred_key or "").strip().upper()
    for key in (pk, "VELOSYS", "VFRAME", "V_CORR_KMS"):
        if key and key not in candidates:
            candidates.append(key)
    return candidates


def _extract_vcorr_kms(row: pd.Series, preferred_key: str = "VELOSYS") -> tuple[float, str | None]:
    """Return row-wise unapplied velocity correction in km/s.

    Policy:
      - VELOSYS / VFRAME / V_CORR_KMS are row-only metadata.
      - Global meta/header values must not be used here because a single header value
        cannot represent row-varying corrections safely.
      - 0.0 is a valid stored correction and must not be treated as missing.
    """
    candidates = _vcorr_candidate_names(preferred_key)

    found: list[tuple[str, float]] = []
    for key in candidates:
        val = _row_only_value(row, key, None)
        if val is None or pd.isna(val):
            continue
        found.append((key, float(val) * _vcorr_scale_to_kms(key)))

    if not found:
        return 0.0, None

    ref_key, ref_val = found[0]
    for key, val in found[1:]:
        if np.isfinite(ref_val) and np.isfinite(val) and abs(val - ref_val) > 1.0e-9:
            raise ValueError(f"Velocity columns disagree for a row: {ref_key}={ref_val} km/s, {key}={val} km/s")
    return ref_val, ref_key


def _row_wcs_meta(row: pd.Series, meta: dict, nchan: int) -> dict:
    row_meta = dict(meta)
    # VELOSYS / VFRAME / V_CORR_KMS are intentionally row-only metadata.
    # Remove any stale global values first so they can never leak into the per-row WCS.
    for key in ("VELOSYS", "VFRAME", "V_CORR_KMS"):
        row_meta.pop(key, None)

    for key in (
        "CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1",
        "RESTFRQ", "RESTFREQ", "SPECSYS", "SSYSOBS", "VELDEF",
        "NAXIS1"
    ):
        val = _row_or_meta(row, meta, key, None)
        if val is not None and not pd.isna(val):
            row_meta[key] = val

    for key in ("VELOSYS", "VFRAME", "V_CORR_KMS"):
        val = _row_only_value(row, key, None)
        if val is not None and not pd.isna(val):
            row_meta[key] = val

    row_meta["NAXIS1"] = int(nchan)

    # Keep RESTFRQ/RESTFREQ synchronized
    rest = None
    for key in ("RESTFRQ", "RESTFREQ"):
        val = row_meta.get(key, None)
        if val not in (None, ""):
            try:
                rest = float(val)
                break
            except Exception:
                pass
    if rest is not None:
        row_meta["RESTFRQ"] = rest
        row_meta["RESTFREQ"] = rest

    specsys = _resolve_specsys(row, row_meta)
    if specsys:
        row_meta["SPECSYS"] = specsys

    ssysobs = _resolve_ssysobs(row, row_meta)
    if ssysobs:
        row_meta["SSYSOBS"] = ssysobs
    elif specsys:
        row_meta["SSYSOBS"] = specsys

    return row_meta


def _physical_axis_for_row(row: pd.Series, meta: dict, nchan_full: int, *, v_corr_kms: float, ch_slice: slice | None = None) -> np.ndarray:
    """Return the per-row auxiliary VLSR axis used for BSL_WINF interpretation.

    This helper does *not* regrid or rewrite the spectrum WCS. It only computes the
    row-specific VLSR velocity axis needed to interpret baseline windows when the
    stored spectrum is still in native FREQ/TOPOCENT.
    """
    row_meta = _row_wcs_meta(row, meta, nchan_full)
    axis = vlsrk_axis_for_spectrum(row_meta, v_corr_kms=float(v_corr_kms), nchan=int(nchan_full))
    if ch_slice is not None:
        axis = axis[ch_slice]
    return np.asarray(axis, dtype=float)







def _update_wcs_after_channel_slice(table: pd.DataFrame, meta: dict, data_list: list[np.ndarray], ch_start: Optional[int], ch_stop: Optional[int], sliced_flags: Optional[list[bool]] = None) -> tuple[pd.DataFrame, dict]:
    """
    If run_baseline_fit outputs channel-sliced spectra, update row/header WCS so that
    data length and WCS remain consistent.

    Policy:
      - keep the current "output is cropped" semantics
      - update per-row CRVAL1 / NAXIS1 using the same rule as axis.wcs_slice_channels()
      - keep CRPIX1 unchanged (consistent with axis.py)
    """
    if ch_start is None and ch_stop is None:
        return table, meta

    out_table = table.copy()
    out_meta = dict(meta)

    s = 0 if ch_start is None else int(ch_start)

    # Update row-wise NAXIS1 and CRVAL1 when possible.
    if len(out_table) != len(data_list):
        raise ValueError(f"Internal error: table/data length mismatch after baseline fit ({len(out_table)} vs {len(data_list)}).")

    if sliced_flags is None:
        sliced_flags = [True] * len(data_list)
    if len(sliced_flags) != len(data_list):
        raise ValueError(f"Internal error: sliced_flags/data length mismatch ({len(sliced_flags)} vs {len(data_list)}).")

    if "NAXIS1" not in out_table.columns:
        out_table["NAXIS1"] = np.nan

    for i, spec_i in enumerate(data_list):
        nout = int(len(spec_i))
        out_table.at[i, "NAXIS1"] = nout

        if not bool(sliced_flags[i]):
            continue

        row = out_table.iloc[i]
        row_has_wcs = all(k in row.index and pd.notna(row[k]) for k in ("CRVAL1", "CDELT1"))
        if row_has_wcs:
            row_meta = {
                "CRVAL1": float(row["CRVAL1"]),
                "CDELT1": float(row["CDELT1"]),
                "CRPIX1": float(row["CRPIX1"]) if ("CRPIX1" in row.index and pd.notna(row["CRPIX1"])) else 1.0,
                "NAXIS1": int(row["NAXIS1"]) if pd.notna(row["NAXIS1"]) else int(len(spec_i)),
            }
            row_meta2 = wcs_slice_channels(row_meta, s, s + nout)
            out_table.at[i, "CRVAL1"] = float(row_meta2["CRVAL1"])
            out_table.at[i, "CRPIX1"] = float(row_meta2.get("CRPIX1", row_meta["CRPIX1"]))
            out_table.at[i, "NAXIS1"] = int(row_meta2["NAXIS1"])

    # Update global meta only when every output row was channel-sliced, the output is
    # uniform-length, and the global meta has a usable frequency WCS.
    try:
        lengths = [int(len(x)) for x in data_list]
        if (
            lengths
            and all(bool(f) for f in sliced_flags)
            and all(n == lengths[0] for n in lengths)
            and ("CRVAL1" in out_meta)
            and ("CDELT1" in out_meta)
        ):
            out_meta = wcs_slice_channels(
                {
                    "CRVAL1": float(out_meta["CRVAL1"]),
                    "CDELT1": float(out_meta["CDELT1"]),
                    "CRPIX1": float(out_meta.get("CRPIX1", 1.0)),
                    "NAXIS1": int(out_meta.get("NAXIS1", lengths[0] + s)),
                },
                s,
                s + lengths[0],
            ) | {k: v for k, v in out_meta.items() if k not in ("CRVAL1", "CDELT1", "CRPIX1", "NAXIS1")}
    except Exception:
        # Keep row-wise WCS as the source of truth even if global meta update fails.
        pass

    return out_table, out_meta


def run_baseline_fit(
    input_data: Union[str, Scantable],
    output_path: Optional[str] = None,
    *,
    # [ADDED] rows and exclude_rows (mutually exclusive)
    rows: Union[str, slice, Sequence[int], int, None] = None,
    exclude_rows: Union[str, slice, Sequence[int], int, None] = None,
    vwin: List[str],
    poly_order: int = 1,
    line_vwin: Optional[List[str]] = None,
    iter_max: int = 0,
    iter_sigma: float = 3.0,
    max_dumps: int = 0,
    ch_start: Optional[int] = None, # ch_slice for VLA is tricky, usually ignored or applied if uniform
    ch_stop: Optional[int] = None,
    # Standard policy: VELOSYS is the canonical unapplied velocity column [m/s].
    v_corr_col: str = "VELOSYS",
    rest_freq: Optional[float] = None,
    apply: bool = True,
    bsl_overwrite: str = "replace",
    on_fail: str = "exit",
    overwrite: bool = True
) -> Scantable:
    """
    Apply or evaluate a polynomial baseline on a Scantable (VLA supported).
    """
    policy = FailPolicy(on_fail)

    apply = bool(apply)
    bsl_overwrite = str(bsl_overwrite).strip().lower()
    if bsl_overwrite not in {"replace", "error"}:
        raise ValueError("bsl_overwrite must be 'replace' or 'error'.")

    # 1. Load Data
    if isinstance(input_data, str):
        sc = read_scantable(input_data)
        input_name = input_data
    else:
        sc = input_data
        input_name = "memory_object"

    meta = dict(sc.meta)
    # data can be np.ndarray OR List[np.ndarray]
    data_source = sc.data 
    table = sc.table.copy()      # [MODIFIED] コピーして汚染防止
    hist = dict(sc.history)      # [MODIFIED] コピーして汚染防止

    if bsl_overwrite == "error" and any(str(c).startswith("BSL_") for c in table.columns):
        raise ValueError("BSL_* columns already exist; bsl_overwrite='error' forbids replacing them.")
    
    # [ADDED] Ensure table is native endian immediately to prevent iloc issues
    table = _df_to_native_endian(table)

    # [ADDED] REST frequency override (rest1->rest2) for velocity interpretation
    if rest_freq is not None:
        apply_restfreq_override(meta, table, float(rest_freq), require_wcs_for_vrad=True)


    total_len = len(table)

    # [ADDED] Row selection / exclusion logic
    if rows is not None and exclude_rows is not None:
        raise ValueError("Cannot specify both 'rows' and 'exclude_rows'.")

    final_idxs = None
    if rows is not None:
        final_idxs = _parse_row_selector(rows, total_len)
    elif exclude_rows is not None:
        exclude_idxs = _parse_row_selector(exclude_rows, total_len)
        # Difference set (np.setdiff1d returns unique and sorted)
        final_idxs = np.setdiff1d(np.arange(total_len), exclude_idxs)
    
    # [MODIFIED] Ensure sorted and unique (Time order preservation) for ALL cases
    if final_idxs is not None:
        final_idxs = np.unique(final_idxs)
        
        if len(final_idxs) == 0:
            raise ValueError("Selector resulted in empty table.")
        
        # Safe slicing (table is already native endian)
        table = table.iloc[final_idxs].reset_index(drop=True)
        
        if isinstance(data_source, list):
            data_source = [data_source[i] for i in final_idxs]
        else:
            data_source = data_source[final_idxs]
        
        # Update n_rows for subsequent max_dumps logic
        total_len = len(table)

    # Handle max_dumps
    n_rows = total_len
    if max_dumps > 0 and n_rows > max_dumps:
        n_rows = int(max_dumps)
        table = table.iloc[:n_rows].reset_index(drop=True)
        # Slice data
        if isinstance(data_source, list):
            data_source = data_source[:n_rows]
        else:
            data_source = data_source[:n_rows]

    # ---------------------------------------------------------------------
    # Temperature scale policy (New 2026: Non-destructive)
    # - Just ensure TEMPSCAL/BEAMEFF columns exist for downstream consistency.
    # - Do NOT convert data even if it says "TR*". We fit on whatever it is.
    # ---------------------------------------------------------------------
    if "TEMPSCAL" not in table.columns:
        table = table.copy()
        table["TEMPSCAL"] = normalize_tempscal(meta.get("TEMPSCAL", "TA*"), default="TA*")
    else:
        # Normalize strings but don't change meaning
        table = table.copy()
        try:
            table["TEMPSCAL"] = [normalize_tempscal(v, default="TA*") for v in table["TEMPSCAL"].to_numpy()]
        except Exception:
            table["TEMPSCAL"] = normalize_tempscal("TA*", default="TA*")

    # Ensure BEAMEFF exists (as NaN is fine)
    if "BEAMEFF" not in table.columns and "BEAMEFF" in meta:
        try:
            table["BEAMEFF"] = float(meta["BEAMEFF"])
        except Exception:
            table["BEAMEFF"] = np.nan

    # [REMOVED] The logic that forced conversion from TR* to TA* has been deleted.
    # We now trust the input data as-is.

    # VLA Data Access Helper
    def get_spec(idx):
        if isinstance(data_source, list):
            return data_source[idx]
        else:
            return data_source[idx]

    # 2. Parse Windows
    base_windows_req = parse_windows(vwin)
    line_windows_req = parse_windows(line_vwin) if line_vwin else []
    
    # 3. Fit Loop
    new_data_list = []
    coeffs_list = []
    rms_list = []
    nused_list = []
    done_list = []
    sliced_flags = []
    keep_mask = [] # For skipping rows if policy=skip

    for i in range(n_rows):
        spec = get_spec(i)
        
        did_apply_channel_slice = False
        try:
            row_i = table.iloc[i]
            full_spec = np.asarray(spec, dtype=float)
            n_chan_full = len(full_spec)

            row_meta = _row_wcs_meta(row_i, meta, n_chan_full)
            specsys = _resolve_specsys(row_i, row_meta)
            if not specsys:
                raise ValueError("Missing SPECSYS/SSYSOBS; baseline fitting requires an explicit spectral reference frame.")
            rest_hz = row_meta.get("RESTFRQ", row_meta.get("RESTFREQ", None))
            if rest_hz in (None, ""):
                raise ValueError("Missing RESTFRQ/RESTFREQ; baseline fitting in velocity windows requires a rest frequency.")
            v_corr_kms, v_corr_source = _extract_vcorr_kms(row_i, preferred_key=v_corr_col)

            if "TOPO" in specsys and v_corr_source is None:
                cand_txt = ", ".join(_vcorr_candidate_names(v_corr_col))
                raise ValueError(
                    f"Row {i} has SPECSYS={specsys} but no usable row-wise velocity correction column "
                    f"({cand_txt})."
                )
            if "LSR" in specsys and abs(v_corr_kms) > 1.0e-9:
                raise ValueError(
                    f"Row {i} has SPECSYS={specsys} but contains non-zero unapplied velocity correction "
                    f"({v_corr_source}={v_corr_kms} km/s)."
                )

            if n_chan_full <= 0:
                raise ValueError("Encountered an empty spectrum.")

            ch_slice = slice(ch_start, ch_stop) if (ch_start is not None or ch_stop is not None) else None
            x_axis = _physical_axis_for_row(
                row_i,
                meta,
                n_chan_full,
                v_corr_kms=float(v_corr_kms),
                ch_slice=ch_slice,
            )

            if ch_slice is not None:
                spec = full_spec[ch_slice]
                did_apply_channel_slice = True
            else:
                spec = full_spec

            if len(spec) != len(x_axis):
                raise ValueError(
                    f"Axis/data length mismatch after channel selection: len(spec)={len(spec)}, len(axis)={len(x_axis)}"
                )
            if not np.all(np.isfinite(x_axis)):
                raise ValueError("Generated spectral axis contains non-finite values.")

            # Fit
            coeffs, baseline, info = fit_polynomial_baseline(
                x_axis, spec, 
                base_windows=base_windows_req,
                line_windows=line_windows_req,
                poly_order=int(poly_order),
                iter_max=int(iter_max),
                iter_sigma=float(iter_sigma),
            )
            
            # Result
            resid = spec - baseline
            rms_val = info.get("std", np.nan) # Use STD as RMS

            out_spec = resid if apply else spec
            new_data_list.append(np.asarray(out_spec, dtype=np.float32))
            coeffs_list.append(np.asarray(coeffs, dtype=float))
            rms_list.append(rms_val)
            nused_list.append(int(np.count_nonzero(info.get("mask", np.zeros_like(spec, dtype=bool)))))
            done_list.append(True)
            sliced_flags.append(did_apply_channel_slice)
            keep_mask.append(True)

        except Exception as e:
            if policy.mode == "exit":
                raise ValueError(f"Baseline fit failed at row {i}: {e}") from e
            else:
                # Fill Dummy
                # Keep row but mark as failed? Or skip?
                # If "skip", we remove it later.
                # If "warn", we keep original data.
                keep_mask.append(False if policy.mode == "skip" else True)
                new_data_list.append(np.asarray(spec, dtype=np.float32)) # No subtract
                coeffs_list.append(None)
                rms_list.append(np.nan)
                nused_list.append(0)
                done_list.append(False)
                sliced_flags.append(did_apply_channel_slice)
                if policy.mode == "warn":
                    print(f"Warning: Baseline fit failed row {i}, kept original.")

    # 4. Handle Skipped Rows
    if policy.mode == "skip":
        keep_bool = np.array(keep_mask, dtype=bool)
        if not np.all(keep_bool):
            table = table.iloc[keep_bool].reset_index(drop=True)
            new_data_list = [d for d, k in zip(new_data_list, keep_bool) if k]
            coeffs_list = [c for c, k in zip(coeffs_list, keep_bool) if k]
            rms_list = [r for r, k in zip(rms_list, keep_bool) if k]
            nused_list = [u for u, k in zip(nused_list, keep_bool) if k]
            done_list = [d for d, k in zip(done_list, keep_bool) if k]
            sliced_flags = [f for f, k in zip(sliced_flags, keep_bool) if k]
            n_rows = len(table)

    # 5. Reconstruct Scantable Data
    # If original was ndarray and lengths are uniform, try to keep as ndarray
    # Otherwise leave as list (VLA)
    is_vla = False
    if new_data_list:
        l0 = len(new_data_list[0])
        if any(len(d) != l0 for d in new_data_list):
            is_vla = True
    
    if is_vla:
        final_data = new_data_list
    else:
        if new_data_list:
            final_data = np.vstack(new_data_list)
        else:
            final_data = np.zeros((0, 0)) # Empty

    # 6. Update WCS for channel-sliced output (if any)
    table, meta = _update_wcs_after_channel_slice(table, meta, new_data_list, ch_start, ch_stop, sliced_flags)

    # 7. Update Table & History
    table = _drop_bsl_metadata_columns(table)
    vwin_str = ";".join(vwin) if vwin else ""
    table = _apply_baseline_to_table(
        table, 
        coeffs_list, 
        rms_list, 
        poly_orders=poly_order,
        win_f_list=vwin_str,
        applied=apply,
        stage="baseline_fit",
        done_list=done_list,
        nused_list=nused_list,
    )

    history_entry = {
        "stage": "baseline_fit",
        "created_at_utc": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "input": input_name,
        "parameters": {
            "poly_order": int(poly_order),
            "rows": str(rows),
            "exclude_rows": str(exclude_rows),
            "vwin": vwin,
            "line_vwin": line_vwin,
            "iter_max": int(iter_max),
            "iter_sigma": float(iter_sigma),
            "on_fail": str(on_fail),
            "max_dumps": int(max_dumps),
            "ch_start": None if ch_start is None else int(ch_start),
            "ch_stop": None if ch_stop is None else int(ch_stop),
            "v_corr_col": str(v_corr_col),
            "rest_freq": None if rest_freq is None else float(rest_freq),
            "apply": bool(apply),
            "bsl_overwrite": str(bsl_overwrite),
            "window_frame": "VLSR",
            "vcorr_source_policy": "row_only",
        }
    }
    
    new_history = hist.copy()
    if "baseline_history" not in new_history:
        new_history["baseline_history"] = []
    
    if isinstance(new_history["baseline_history"], list):
        new_history["baseline_history"].append(history_entry)
    
    res = Scantable(meta=meta, data=final_data, table=table, history=new_history)

    if output_path:
        write_scantable(output_path, res, overwrite=overwrite)

    return res
