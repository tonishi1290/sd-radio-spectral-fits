# src/sd_radio_spectral_fits/tempscale.py
from __future__ import annotations

import datetime
import warnings
from typing import Any, Iterable, Optional, Tuple, Union, List, TYPE_CHECKING

import numpy as np
import pandas as pd

# [MODIFIED] Avoid circular import by using TYPE_CHECKING
if TYPE_CHECKING:
    from .fitsio import Scantable

# -----------------------------------------------------------------------------
# Temperature scale utilities
# -----------------------------------------------------------------------------
# Policy:
# - Internal spectral arrays are treated as Ta*.
# - Tr* is derived as Ta* / BEAMEFF (main-beam efficiency).
# - The on-disk scale is declared by TEMPSCAL ("TA*" or "TR*").
# - Never switch scales silently: if conversion is performed, record it in HISTORY/log.
# -----------------------------------------------------------------------------


def normalize_tempscal(value: Any, *, default: str = "TA*") -> str:
    """Normalize various inputs into canonical 'TA*' or 'TR*'."""
    if value is None:
        return str(default).upper()
    s = str(value).strip().upper()
    if not s:
        return str(default).upper()
    # Common variants
    if s in ("TA*", "TA", "TASTAR", "TASTAR*", "T_A*", "T_A", "TA-STAR", "TASTAR-STAR"):
        return "TA*"
    if s in ("TR*", "TR", "TMB", "T_MB", "TMAINBEAM", "MAINBEAM", "TMB*"):
        return "TR*"
    # Unknown -> default (safer than propagating nonsense)
    return str(default).upper()


def ensure_tempscal_column(df: pd.DataFrame, *, default: str = "TA*") -> pd.DataFrame:
    """Ensure df has a TEMPSCAL column (string)."""
    if df is None:
        return df
    if "TEMPSCAL" not in df.columns:
        df = df.copy()
        df["TEMPSCAL"] = normalize_tempscal(default, default=default)
        return df
    # normalize in-place on a copy (avoid modifying caller unexpectedly)
    out = df.copy()
    try:
        out["TEMPSCAL"] = [normalize_tempscal(v, default=default) for v in out["TEMPSCAL"].to_numpy()]
    except Exception:
        out["TEMPSCAL"] = normalize_tempscal(default, default=default)
    return out


def ensure_beameff_column(df: pd.DataFrame, *, default: float = np.nan) -> pd.DataFrame:
    """Ensure df has a BEAMEFF column (float)."""
    if df is None:
        return df
    if "BEAMEFF" not in df.columns:
        df = df.copy()
        df["BEAMEFF"] = float(default)
        return df
    out = df.copy()
    try:
        out["BEAMEFF"] = pd.to_numeric(out["BEAMEFF"], errors="coerce").astype(float)
    except Exception:
        out["BEAMEFF"] = float(default)
    return out


def beameff_array(
    df: pd.DataFrame,
    meta: dict | None,
    n_rows: int,
    *,
    default: float = np.nan,
) -> np.ndarray:
    """Return BEAMEFF per row (length n_rows)."""
    if df is not None and "BEAMEFF" in df.columns:
        try:
            arr = pd.to_numeric(df["BEAMEFF"], errors="coerce").to_numpy(dtype=float)
            if arr.size >= n_rows:
                return arr[:n_rows]
            if arr.size > 0:
                out = np.full(int(n_rows), float(default), dtype=float)
                out[:arr.size] = arr
                return out
        except Exception:
            pass

    # fallback: global keyword
    if meta is not None:
        for k in ("BEAMEFF",):
            if k in meta and meta[k] not in (None, ""):
                try:
                    v = float(meta[k])
                    return np.full(int(n_rows), v, dtype=float)
                except Exception:
                    pass

    return np.full(int(n_rows), float(default), dtype=float)


def tempscal_array(
    df: pd.DataFrame,
    meta: dict | None,
    n_rows: int,
    *,
    default: str = "TA*",
) -> np.ndarray:
    """Return TEMPSCAL per row (length n_rows) as normalized strings."""
    if df is not None and "TEMPSCAL" in df.columns:
        try:
            col = df["TEMPSCAL"].to_numpy()
            arr = np.array([normalize_tempscal(v, default=default) for v in col], dtype="U4")
            if arr.size >= n_rows:
                return arr[:n_rows]
            out = np.full(int(n_rows), normalize_tempscal(default, default=default), dtype="U4")
            out[:arr.size] = arr
            return out
        except Exception:
            pass

    # fallback: header keyword
    if meta is not None:
        for k in ("TEMPSCAL", "TEMP_SCAL", "TSCALE"):
            if k in meta and meta[k] not in (None, ""):
                return np.full(int(n_rows), normalize_tempscal(meta[k], default=default), dtype="U4")
    return np.full(int(n_rows), normalize_tempscal(default, default=default), dtype="U4")


def require_beameff(beameff: np.ndarray, *, on_fail: str = "error") -> None:
    """
    Validate BEAMEFF (main-beam efficiency).

    - Must be finite and > 0 for Tr* conversions.
    - Internal spectral arrays are treated as Ta*; Tr* is derived as Ta* / BEAMEFF.
    """
    b = np.asarray(beameff, dtype=float)
    bad = (~np.isfinite(b)) | (b <= 0.0)
    if np.any(bad):
        msg = "BEAMEFF is required (finite and > 0) for TR* conversions, but missing/invalid values were found."
        if on_fail == "warn":
            warnings.warn(msg, RuntimeWarning)
            return
        raise ValueError(msg)


def is_beameff_mixed(beameff: np.ndarray, *, tol: float = 1e-6) -> bool:
    """
    Check whether BEAMEFF is mixed beyond tolerance.

    `tol` is a *relative* tolerance:
        mixed if (max - min) > tol * median(|beameff|).

    Any non-positive or non-finite BEAMEFF is treated as mixed (because TR* conversion would be unsafe).
    """
    b = np.asarray(beameff, dtype=float).ravel()
    b = b[np.isfinite(b)]
    if b.size <= 1:
        return False
    bmin = float(np.min(b))
    bmax = float(np.max(b))
    if bmin <= 0.0:
        return True
    bref = float(np.median(b))
    scale = max(abs(bref), 1e-12)
    return (bmax - bmin) > (float(tol) * scale)


def representative_beameff(beameff: np.ndarray) -> float:
    """Return a representative BEAMEFF (median of finite positive values)."""
    b = np.asarray(beameff, dtype=float).ravel()
    b = b[np.isfinite(b) & (b > 0.0)]
    if b.size == 0:
        return float("nan")
    return float(np.median(b))


def ta_to_tr(ta: np.ndarray, beameff: np.ndarray) -> np.ndarray:
    """Convert Ta* -> Tr* (vectorized): Tr* = Ta* / BEAMEFF."""
    require_beameff(beameff, on_fail="error")
    ta = np.asarray(ta, dtype=float)
    b = np.asarray(beameff, dtype=float)

    if ta.ndim == 2 and b.ndim == 1 and ta.shape[0] == b.shape[0]:
        out = np.empty_like(ta, dtype=float)
        np.divide(ta, b[:, None], out=out, where=(b[:, None] > 0.0))
        return out.astype(np.float32, copy=False)

    out = ta / b
    return np.asarray(out, dtype=np.float32)


def tr_to_ta(tr: np.ndarray, beameff: np.ndarray) -> np.ndarray:
    """Convert Tr* -> Ta* (vectorized): Ta* = Tr* * BEAMEFF."""
    require_beameff(beameff, on_fail="error")
    tr = np.asarray(tr, dtype=float)
    b = np.asarray(beameff, dtype=float)

    if tr.ndim == 2 and b.ndim == 1 and tr.shape[0] == b.shape[0]:
        out = (tr * b[:, None]).astype(np.float32, copy=False)
        return out

    out = tr * b
    return np.asarray(out, dtype=np.float32)


def convert_rowwise_vla(
    spectra: List[np.ndarray],
    beameff: np.ndarray,
    *,
    direction: str,
) -> List[np.ndarray]:
    """Row-wise Ta*<->Tr* conversion for ragged (variable-length) spectra.

    Parameters
    ----------
    spectra
        List of 1D numpy arrays (one per row). Each element may have different length.
    beameff
        Per-row BEAMEFF array. If length==1, it is broadcast to all rows.
    direction
        "ta_to_tr" or "tr_to_ta".

    Returns
    -------
    list of np.ndarray
        Converted spectra as float32 arrays, preserving original per-row shapes.
    """
    if not isinstance(spectra, list):
        spectra = list(spectra)

    b = np.asarray(beameff, dtype=float).ravel()
    if b.size == 1 and len(spectra) > 1:
        b = np.full(len(spectra), float(b[0]), dtype=float)
    if len(spectra) != int(b.size):
        raise ValueError(f"convert_rowwise_vla: len(spectra)={len(spectra)} != len(BEAMEFF)={int(b.size)}")

    require_beameff(b, on_fail="error")

    d = str(direction).strip().lower()
    out: List[np.ndarray] = []
    if d in ("ta_to_tr", "ta2tr", "ta->tr"):
        for spec, bi in zip(spectra, b):
            arr = np.asarray(spec, dtype=float)
            out.append((arr / float(bi)).astype(np.float32, copy=False))
        return out
    if d in ("tr_to_ta", "tr2ta", "tr->ta"):
        for spec, bi in zip(spectra, b):
            arr = np.asarray(spec, dtype=float)
            out.append((arr * float(bi)).astype(np.float32, copy=False))
        return out

    raise ValueError("convert_rowwise_vla: direction must be 'ta_to_tr' or 'tr_to_ta'")



def append_scale_history(history: dict | None, entry: dict) -> dict:
    """Append a scale-related history entry (non-breaking)."""
    if history is None or not isinstance(history, dict):
        history = {} if history is None else {"_history": str(history)}
    hist = dict(history)
    key = "scale_history"
    if key not in hist or not isinstance(hist.get(key), list):
        hist[key] = []
    entry2 = dict(entry)
    entry2.setdefault("created_at_utc", datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z")
    hist[key].append(entry2)
    return hist


# [MODIFIED] Added safe metadata setter based on New Spec
def set_beameff(
    sc: "Scantable", # Use string forward reference or TYPE_CHECKING alias
    efficiency: Union[float, np.ndarray, Sequence[float]],
    rows: Union[str, List[int], None] = None,
    verbose: bool = True
) -> None:
    """
    BEAMEFF (ビーム効率) を設定するユーティリティ。
    
    注意:
    この関数はデータの実体や TEMPSCAL (温度スケールラベル) を変更しません。
    単に BEAMEFF カラムに値をセットし、Viewer等でのオンザフライ変換 (Ta* <-> TR*) を可能にするためのメタデータを供給します。

    Parameters
    ----------
    sc : Scantable
        対象のScantable
    efficiency : float or array
        設定するビーム効率 (0.0 < eta <= 1.0)
    rows : selector
        適用する行。Noneの場合は全行。
    verbose : bool
        詳細表示。
    """
    df = sc.table
    
    # 1. 効率値のバリデーション
    eff_arr = np.array(efficiency, dtype=float)
    if np.any((eff_arr <= 0) | (eff_arr > 1.0)):
        # 物理的にありえない値ですが、計算エラーを防ぐため警告にとどめる
        warnings.warn("Efficiency values out of physical range (0.0 < eta <= 1.0). Check your inputs.")
    
    # 2. 行の特定 (簡易実装フォールバック付き)
    try:
        from .scantable_utils import _parse_row_selector
        idxs = _parse_row_selector(rows, len(df))
    except ImportError:
        if rows is None:
            idxs = df.index
        elif isinstance(rows, (list, np.ndarray)):
            idxs = rows
        elif isinstance(rows, slice):
            idxs = np.arange(len(df))[rows]
        else:
            idxs = df.index 

    if len(idxs) == 0:
        return

    # 3. BEAMEFF の更新
    
    # BEAMEFFカラムがなければ作成
    if "BEAMEFF" not in df.columns:
        df["BEAMEFF"] = np.nan

    # 値の代入
    try:
        # Pandasのloc代入。eff_arrがスカラーならブロードキャスト、配列ならサイズ一致が必要
        df.loc[df.index[idxs], "BEAMEFF"] = eff_arr
        
        if verbose:
            print(f"Updated BEAMEFF for {len(idxs)} rows.")
            # ユーザーへの案内（TEMPSCALは変わっていないことを明示）
            current_scale = df["TEMPSCAL"].iloc[idxs[0]] if "TEMPSCAL" in df.columns else "Unknown"
            print(f"Note: Data remains in '{current_scale}'. Use Viewer to toggle TR* display.")
            
    except Exception as e:
        print(f"Error setting BEAMEFF: {e}")
        return
