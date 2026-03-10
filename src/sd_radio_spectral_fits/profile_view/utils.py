# src/sd_radio_spectral_fits/plotting/utils.py
from __future__ import annotations

import ast
import builtins
from dataclasses import dataclass
from typing import Optional, Tuple, List, Union, Any, Callable
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd

# 親ディレクトリのモジュールを相対インポート
from ..axis import freq_axis_from_wcs, radio_velocity_kms
from ..ranges import parse_windows
from ..doppler import C_KMS # [MODIFIED] Use shared constant

# [ADDED] Paper Sizes (inch)
# NOTE:
#   - Default keys (e.g., "A4", "A3", "LETTER") are PORTRAIT.
#   - Use *P / *L suffix to force Portrait/Landscape explicitly (e.g., A4P, A4L).
#   - Backward-compatible aliases with underscores (e.g., A4_P) are also accepted.
PAPER_SIZES = {
    # ISO A-series
    "A4":  (8.27, 11.69),  # Portrait (default)
    "A4P": (8.27, 11.69),  # Portrait
    "A4L": (11.69, 8.27),  # Landscape
    "A4_P": (8.27, 11.69), # alias
    "A4_L": (11.69, 8.27), # alias

    "A3":  (11.69, 16.54), # Portrait (default)
    "A3P": (11.69, 16.54), # Portrait
    "A3L": (16.54, 11.69), # Landscape
    "A3_P": (11.69, 16.54),# alias
    "A3_L": (16.54, 11.69),# alias

    # US Letter
    "LETTER":  (8.5, 11.0),   # Portrait (default)
    "LETTERP": (8.5, 11.0),   # Portrait
    "LETTERL": (11.0, 8.5),   # Landscape
    "LETTER_P": (8.5, 11.0),  # alias
    "LETTER_L": (11.0, 8.5),  # alias
}



# -------------------------
# Axis bundle (Common Data Structure for Legacy Viewers)
# -------------------------
@dataclass(frozen=True)
class AxisBundle:
    chan: np.ndarray
    freq_hz: Optional[np.ndarray]
    specsys: str
    ctype1: str

    @staticmethod
    def build(meta: dict, nchan: int) -> "AxisBundle":
        ctype1 = str(meta.get("CTYPE1", "")).upper()
        specsys = str(meta.get("SPECSYS", "")).upper()
        
        # Determine if native axis is velocity
        is_velo = ctype1.startswith(("VRAD", "VELO"))
        
        if is_velo:
            # Native velocity: Freq conversion might be ambiguous without RESTFREQ
            # For visualization, we keep freq_hz=None unless explicitly needed
            return AxisBundle(chan=np.arange(nchan, dtype=float), freq_hz=None, specsys=specsys, ctype1=ctype1)
        else:
            # Native frequency
            f = freq_axis_from_wcs(meta, nchan)
            return AxisBundle(chan=np.arange(nchan, dtype=float), freq_hz=np.asarray(f, dtype=float), specsys=specsys, ctype1=ctype1)

    def vel_uncorr(self, rest_hz: Optional[float]) -> Optional[np.ndarray]:
        """Calculate velocity from freq without additional v_corr (Topocentric)."""
        if self.ctype1.startswith(("VRAD", "VELO")):
            return None # Already velocity (values depend on CRVAL/CDELT)
        if rest_hz is None or rest_hz <= 0 or self.freq_hz is None:
            return None
        return radio_velocity_kms(self.freq_hz, rest_hz)

# -------------------------
# QC / Policy Helpers
# -------------------------
@dataclass
class FailPolicy:
    mode: str  # "exit", "warn", "skip"

    def __init__(self, mode: str):
        m = str(mode).lower().strip()
        if m not in ("exit", "warn", "skip"):
            m = "exit"
        self.mode = m

# -------------------------
# Helper Functions
# -------------------------
def _norm_range(r) -> Optional[Tuple[float, float]]:
    if r is None:
        return None
    if not isinstance(r, (tuple, list)) or len(r) != 2:
        # raise ValueError("xrange/yrange must be (min,max)")
        return None
    a, b = r
    if a is None or b is None:
        return None
    a = float(a); b = float(b)
    return (a, b) if a <= b else (b, a)

def _finite(a: np.ndarray) -> np.ndarray:
    return np.isfinite(a)

def _rms_from_mean(y: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return float("nan")
    mu = float(np.mean(y))
    return float(np.sqrt(np.mean((y - mu) ** 2)))

def _parse_list_like(obj: Any) -> Optional[np.ndarray]:
    """
    Parse a DataFrame cell that might contain a list, numpy array, or string representation.
    Used for BSL_COEF.
    """
    if obj is None:
        return None
    
    # Already an array?
    if isinstance(obj, np.ndarray):
        return obj.astype(float, copy=False) if obj.size > 0 else None
        
    # List or Tuple?
    if isinstance(obj, (list, tuple)):
        if len(obj) == 0: return None
        return np.asarray(obj, dtype=float)
        
    # String representation? (Legacy FITS string column)
    if isinstance(obj, str):
        s = obj.strip()
        if (not s) or s.lower() in ("nan", "none", "null"):
            return None
        try:
            # Safer than eval, handles "[1, 2, 3]"
            v = ast.literal_eval(s)
            if isinstance(v, (list, tuple)):
                return np.asarray(v, dtype=float)
        except Exception:
            pass
            
    # Float (single value)?
    if isinstance(obj, (float, int, np.number)):
        if np.isnan(obj): return None
        return np.array([float(obj)])
        
    return None

def _interp_extrap(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    """Linear interpolation/extrapolation. NaNs in input propagate."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xq = np.asarray(xq, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]; y = y[ok]
    if x.size < 2:
        return np.full_like(xq, np.nan, dtype=float)

    idx = np.argsort(x)
    x = x[idx]; y = y[idx]

    yq = np.interp(xq, x, y, left=np.nan, right=np.nan)

    # Extrapolate left
    if abs(x[1] - x[0]) > 1e-12:
        m0 = (y[1] - y[0]) / (x[1] - x[0])
        left = (xq < x[0])
        if np.any(left):
            yq[left] = y[0] + m0 * (xq[left] - x[0])

    # Extrapolate right
    if abs(x[-1] - x[-2]) > 1e-12:
        m1 = (y[-1] - y[-2]) / (x[-1] - x[-2])
        right = (xq > x[-1])
        if np.any(right):
            yq[right] = y[-1] + m1 * (xq[right] - x[-1])

    return yq

def _vel_from_freq_ghz(freq_ghz: np.ndarray, rest_hz: float) -> np.ndarray:
    f = np.asarray(freq_ghz, dtype=float) * 1e9
    return C_KMS * (1.0 - f / float(rest_hz))

def _freq_ghz_from_vel(vel_kms: np.ndarray, rest_hz: float) -> np.ndarray:
    v = np.asarray(vel_kms, dtype=float)
    f = float(rest_hz) * (1.0 - v / C_KMS)
    return f * 1e-9

def recalculate_velocity_axis(vel_kms_old: np.ndarray, rest_hz_old: float, rest_hz_new: float) -> np.ndarray:
    """
    Recalculate velocity axis based on a new rest frequency.
    
    Logic:
        f_obs = rest_old * (1 - v_old / c)
        v_new = c * (1 - f_obs / rest_new)
              = c * (1 - (rest_old/rest_new) * (1 - v_old/c))
    """
    if rest_hz_old <= 0 or rest_hz_new <= 0:
        return vel_kms_old
    if abs(rest_hz_old - rest_hz_new) < 1.0: # Identical
        return vel_kms_old

    r_ratio = rest_hz_old / rest_hz_new
    # Formula simplified:
    # v_new = c * (1 - r_ratio) + v_old * r_ratio
    return C_KMS * (1.0 - r_ratio) + vel_kms_old * r_ratio


def _running_mean(y: np.ndarray, w: int) -> np.ndarray:
    w = int(w)
    if w <= 1:
        return y
    y = np.asarray(y, dtype=float)
    kern = np.ones(w, dtype=float)
    y0 = np.nan_to_num(y, nan=0.0)
    ww = np.isfinite(y).astype(float)
    
    # Avoid boundary effects if possible, but 'same' is standard
    num = np.convolve(y0, kern, mode="same")
    den = np.convolve(ww, kern, mode="same")
    
    with np.errstate(divide='ignore', invalid='ignore'):
        out = num / den
    out[den == 0] = np.nan
    return out

def _box_average_downsample(x: np.ndarray, y: np.ndarray, w: int, policy: str = "trim"):
    w = int(w)
    if w <= 1:
        return x, y
    n = int(y.size)
    
    if policy == "pad":
        m = int(np.ceil(n / w))
        pad = m * w - n
        if pad > 0:
            x = np.concatenate([x, np.full(pad, np.nan)])
            y = np.concatenate([y, np.full(pad, np.nan)])
        n = int(y.size) # Updated n
    else:
        # Trim
        m = n // w
        n2 = m * w
        x = x[:n2]; y = y[:n2]
    
    if m <= 0:
        return np.array([]), np.array([])
        
    xx = x.reshape(m, w)
    yy = y.reshape(m, w)
    return np.nanmean(xx, axis=1), np.nanmean(yy, axis=1)

def _process_spectrum(x, y, smooth_mode, smooth_width, box_downsample, box_policy):
    mode = str(smooth_mode).lower()
    w = int(smooth_width)
    if mode in ("none", "") or w <= 1:
        return x, y
    if mode in ("running", "run", "mean"):
        return x, _running_mean(y, w)
    if mode in ("boxcar", "box", "avg", "average"):
        if box_downsample:
            return _box_average_downsample(x, y, w, policy=box_policy)
        return x, _running_mean(y, w)
    return x, y

def in_any_windows(v_axis: np.ndarray, windows: List[Tuple[float, float]]) -> np.ndarray:
    """Return boolean mask for v_axis in any of the windows."""
    mask = np.zeros_like(v_axis, dtype=bool)
    for (w0, w1) in windows:
        lo, hi = (w0, w1) if w0 <= w1 else (w1, w0)
        mask |= (v_axis >= lo) & (v_axis <= hi)
    return mask

def subtract_windows(
    base_windows: List[Tuple[float, float]], 
    exclude_windows: List[Tuple[float, float]]
) -> List[Tuple[float, float]]:
    """
    Subtract exclude_windows ranges from base_windows.
    Returns a new list of disjoint windows.
    Simple implementation: works on continuum logic (float ranges).
    """
    # This can be complex to do geometrically.
    # An easier way for spectral fitting context is to rely on masks,
    # but if we need the window list back, we do 1D boolean geometry.
    # Here is a simplified approach:
    # 1. We won't strictly compute geometry here to keep util simple.
    #    Instead, baseline.py uses `in_any_windows` to build a positive mask,
    #    and then `in_any_windows(exclude)` to build a negative mask.
    #    So this function might be redundant if we handle it via masks.
    #    However, if we strictly need windows:
    # NOTE: This helper is currently a no-op for backward compatibility.
    # Prefer mask-based logic (in_any_windows) when fitting baselines.
    return base_windows

def _fit_baseline_poly(row: pd.Series, vel: np.ndarray, base_arr: np.ndarray, finite_mask: np.ndarray):
    """
    Helper to re-fit a polynomial if only the baseline array (not coeffs) is available.
    Used by montage viewer for legacy data.
    """
    deg = None
    for k in ("BSL_DEG", "BSL_ORDER", "BSL_NPOLY", "BSL_POLYORDER"):
        if k in row.index:
            try:
                deg = int(row.get(k))
                break
            except Exception:
                pass
    if deg is None:
        deg = 3 # Default guess
        
    n = int(np.count_nonzero(finite_mask))
    if n <= deg + 1:
        return np.array([0.0], dtype=float)
        
    # Limit degree to avoid overfitting on small segments
    deg = builtins.max(0, builtins.min(int(deg), n - 1, 9))
    
    try:
        co = np.polyfit(vel[finite_mask], base_arr[finite_mask], deg=deg)
        return co.astype(float, copy=False)
    except Exception:
        return np.array([0.0], dtype=float)

def _baseline_extrapolated_line(
    row: pd.Series,
    vel: np.ndarray,
    chan: np.ndarray,
    freq_ghz: Optional[np.ndarray],
    rest_hz: Optional[float],
    vel_range: Tuple[float, float],
    xaxis: str
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reconstruct the baseline curve for display.
    Handles BSL_COEF being either polynomial coefficients (new) or full array (legacy).
    """
    if "BSL_COEF" not in row.index:
        return np.array([]), np.array([])
        
    obj = row.get("BSL_COEF")
    arr = _parse_list_like(obj)
    if arr is None:
        return np.array([]), np.array([])
        
    v0, v1 = vel_range
    nline = 400
    vline = np.linspace(v0, v1, nline)
    yline = np.zeros_like(vline)
    
    # Case A: arr is Coefficients (Standard for VLA/New pipeline)
    # Heuristic: Coefficients usually have small size (e.g. < 20) vs Spectrum (> 100)
    if arr.size < vel.size: 
        try:
            yline = np.polyval(arr, vline).astype(float, copy=False)
        except Exception:
            return np.array([]), np.array([])
            
    # Case B: arr is Full Baseline Array (Legacy)
    else:
        # We need to fit a poly to this array to extrapolate it to the view range
        base = arr.astype(float, copy=True)
        # Ensure sizes match
        if base.size != vel.size:
             # Size mismatch -> cannot align -> give up
             return np.array([]), np.array([])
             
        finite_mask = np.isfinite(base)
        co = _fit_baseline_poly(row, vel, base, finite_mask)
        yline = np.polyval(co, vline).astype(float, copy=False)

    # Convert x-axis for display
    xa = str(xaxis).lower()
    x_out = vline
    
    if xa == "freq":
        if rest_hz is not None and rest_hz > 0:
            x_out = _freq_ghz_from_vel(vline, rest_hz)
        elif freq_ghz is not None and np.any(np.isfinite(freq_ghz)):
            x_out = _interp_extrap(vel, freq_ghz, vline)
    elif xa == "chan":
        x_out = _interp_extrap(vel, chan, vline)
        
    return x_out, yline


def _build_fitwin_mask_on_vel(vel: np.ndarray, win_str: str) -> np.ndarray:
    """
    BSL_WINF文字列から速度軸上のマスクを生成する。
    """
    mask = np.zeros_like(vel, dtype=bool)
    if not win_str:
        return mask
        
    s = str(win_str).strip()
    if (not s) or s.lower() in ("nan", "none"):
        return mask
        
    try:
        # 1. セミコロンで区切られた文字列をリスト化
        raw_wins = [w.strip() for w in s.split(";") if w.strip()]
        # 2. ranges.py の parse_windows で数値タプルのリストに変換
        parsed_wins = parse_windows(raw_wins)
        # 3. ranges.py の window_to_mask で一括してマスクを作成（推奨）
        from ..ranges import window_to_mask
        return window_to_mask(vel, parsed_wins)
    except Exception:
        # パースに失敗した場合は全て False のマスクを返す
        return np.zeros_like(vel, dtype=bool)
        
        
def _convert_vel_range_to_xrange(xaxis: str, vel_range: Tuple[float, float],
                                vel: np.ndarray, chan: np.ndarray,
                                freq_ghz: Optional[np.ndarray], rest_hz: Optional[float]) -> Tuple[float, float]:
    v0, v1 = vel_range
    vq = np.array([v0, v1], dtype=float)
    xa = str(xaxis).lower()
    if xa == "vel":
        return (v0, v1)
    if xa == "freq":
        if rest_hz is not None and rest_hz > 0:
            fq = _freq_ghz_from_vel(vq, rest_hz)
            return (float(fq[0]), float(fq[1]))
        if freq_ghz is not None and np.any(np.isfinite(freq_ghz)):
            xq = _interp_extrap(vel, freq_ghz, vq)
            return (float(xq[0]), float(xq[1]))
        return (v0, v1)
    # chan
    chq = _interp_extrap(vel, chan, vq)
    lo = float(np.floor(builtins.min(chq)))
    hi = float(np.ceil(builtins.max(chq)))
    return (lo, hi) if (chq[0] <= chq[1]) else (hi, lo)

def _convert_user_xrange_to_vel_range(
    xaxis: str,
    user_xrange: Tuple[float, float],
    vel: np.ndarray,
    chan: np.ndarray,
    freq_ghz: Optional[np.ndarray],
    rest_hz: Optional[float],
    v_corr_kms: float,
) -> Tuple[float, float]:
    a, b = user_xrange
    xa = str(xaxis).lower()
    if xa == "vel":
        return (float(a), float(b))
    if xa == "freq":
        if rest_hz is not None and rest_hz > 0:
            vq = _vel_from_freq_ghz(np.array([a, b], dtype=float), rest_hz)
            return (float(vq[0]), float(vq[1]))
        if freq_ghz is not None and np.any(np.isfinite(freq_ghz)):
            vq = _interp_extrap(freq_ghz, vel, np.array([a, b], dtype=float))
            return (float(vq[0]), float(vq[1]))
        return (float(a), float(b))
    # chan
    vq = _interp_extrap(chan, vel, np.array([a, b], dtype=float))
    return (float(vq[0]), float(vq[1]))


def parse_figsize(
    size_spec: Union[str, Tuple[float, float], None],
    default: Optional[Tuple[float, float]] = None
) -> Tuple[Optional[Tuple[float, float]], Optional[str]]:
    """
    Parse figsize argument.

    Returns
    -------
    (figsize, bbox_inches_mode)
        - figsize: (width_in, height_in) in inches, or None
        - bbox_inches_mode: None for fixed paper size (do not trim), "tight" otherwise

    Policy
    ------
    - If size_spec is a known paper key (e.g., A4, A4P, A4L) or an explicit (w,h),
      return bbox_inches_mode=None to *respect the fixed size*.
    - If size_spec is None, use default and bbox_inches_mode="tight".
    - If size_spec is an unknown string, try to parse "W,H" / "W x H" style sizes.
      *Units*:
        - If "mm" or "cm" appears anywhere in the string, that unit is used.
        - Else if "in" appears, inches.
        - Else numbers are interpreted as inches (Matplotlib convention).

    Notes
    -----
    - Paper keys are case-insensitive.
    - Common separators are ignored, so "A4_P" and "A4P" are treated the same.
    """
    if size_spec is None:
        return default, "tight"

    # --- string: paper key or size pair ---
    if isinstance(size_spec, str):
        raw = size_spec.strip()
        # Normalize: remove non-alphanumerics so "A4_P" -> "A4P", "letter-l" -> "LETTERL"
        key = re.sub(r"[^A-Z0-9]+", "", raw.upper())
        if key in PAPER_SIZES:
            return PAPER_SIZES[key], None

        # Try to parse explicit size like "8,6", "8x6", "210mm,297mm"
        s = raw.lower()
        if re.search(r"\d", s):
            # Replace common separators with spaces and extract numbers
            s2 = s.replace(",", " ").replace(";", " ")
            s2 = s2.replace("×", " ")
            s2 = re.sub(r"[xX]", " ", s2)
            nums = re.findall(r"[+-]?(?:\d+\.?\d*|\d*\.?\d+)", s2)
            if len(nums) == 2:
                w = float(nums[0]); h = float(nums[1])
                unit = "in"
                if "mm" in s:
                    unit = "mm"
                elif "cm" in s:
                    unit = "cm"
                elif "in" in s or "inch" in s:
                    unit = "in"
                return (_to_inches(w, unit), _to_inches(h, unit)), None

        # Keep behavior: unknown key -> default + tight, but use warnings (library-friendly)
        warnings.warn(
            f"Unknown paper size '{size_spec}'. Using default figsize.",
            UserWarning,
            stacklevel=2,
        )
        return default, "tight"

    # --- tuple/list: explicit (w,h) inches ---
    if isinstance(size_spec, (tuple, list)) and len(size_spec) == 2:
        return (float(size_spec[0]), float(size_spec[1])), None

    return default, "tight"

def drive_pdf_generation(
    filename: Union[str, Path],
    fig: Any,
    updater_func: Callable[[int], Any],
    count: int,
    max_pages: Optional[int] = 100,
    bbox_inches: Optional[str] = "tight",
    pad_inches: float = 0.1, # [ADDED]
    verbose: bool = True
) -> None:
    """
    Helper to generate a multipage PDF by driving an updater function.
    
    Parameters
    ----------
    filename : str or Path
        Output PDF filename.
    fig : matplotlib.figure.Figure
        The figure object to save.
    updater_func : Callable[[int], Any]
        Function that accepts an index (or page number) and updates the figure.
    count : int
        Total number of pages/items to loop through.
    max_pages : int, optional
        Safety limit for maximum pages to output. Default 100. Set to None to disable.
    bbox_inches : str or None
        Passed to savefig. Use 'tight' for screen-optimized, None for fixed paper size.
    pad_inches : float
        Padding around the figure when bbox_inches is 'tight'.
    verbose : bool
        If True, print progress.
    """
    from matplotlib.backends.backend_pdf import PdfPages
    
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    limit = count
    if max_pages is not None and count > max_pages:
        if verbose:
            print(f"⚠️ PDF output limit reached: Outputting first {max_pages} of {count} pages.")
            print(f"   (Use max_pdf_pages=None to disable limit)")
        limit = int(max_pages)
        
    if verbose:
        print(f"Starting PDF generation: {path}")
    
    try:
        with PdfPages(path) as pdf:
            for i in range(limit):
                updater_func(i)
                # [MODIFIED] Use explicit bbox_inches and pad_inches
                pdf.savefig(fig, bbox_inches=bbox_inches, pad_inches=pad_inches)
                if verbose:
                    print(f"Saved PDF page {i+1}/{limit}", end="\r")
        if verbose:
            print(f"\nCompleted. PDF saved to: {path}")
    except Exception as e:
        print(f"\n❌ PDF generation failed: {e}")


# -----------------------------------------------------------------------------
# Paper layout utilities (A4, margins in mm, content_aspect as Axes box aspect)
# -----------------------------------------------------------------------------
_MM_PER_INCH = 25.4

def _to_inches(value: float, unit: str) -> float:
    unit = unit.lower().strip()
    if unit in ("in", "inch", "inches"):
        return float(value)
    if unit in ("mm",):
        return float(value) / _MM_PER_INCH
    if unit in ("cm",):
        return float(value) * 10.0 / _MM_PER_INCH
    raise ValueError(f"Unsupported unit: {unit!r}")

def parse_paper_margins(
    spec: Union[str, Tuple[float, float, float, float], List[float], None],
    default: Union[str, Tuple[float, float, float, float]] = "20mm,15mm,20mm,30mm",
) -> Tuple[float, float, float, float]:
    """Parse paper margins and return (left,right,bottom,top) in inches.

    Accepts:
      - None -> default
      - "20mm,15mm,20mm,25mm" (L,R,B,T)
      - "0.8in,0.6in,0.8in,1.0in"
      - "20,15,20,25" -> interpreted as mm (by design)
      - (20,15,20,25) -> interpreted as mm (by design)
    """
    if spec is None:
        spec = default

    def _parse_one(tok: str) -> float:
        t = tok.strip().lower()
        m = re.fullmatch(r"([+-]?(?:\d+\.?\d*|\d*\.?\d+))(mm|cm|in)?", t)
        if not m:
            raise ValueError(f"Invalid margin token: {tok!r}")
        val = float(m.group(1))
        unit = m.group(2) or "mm"  # numbers without unit -> mm
        return _to_inches(val, unit)

    if isinstance(spec, str):
        parts = [p.strip() for p in spec.replace(" ", "").split(",") if p.strip() != ""]
        if len(parts) != 4:
            raise ValueError("paper_margins must have 4 values: left,right,bottom,top")
        l, r, b, t = (_parse_one(p) for p in parts)
        return (float(l), float(r), float(b), float(t))

    if isinstance(spec, (tuple, list)) and len(spec) == 4:
        # numeric sequence -> mm by design
        l, r, b, t = (float(x) / _MM_PER_INCH for x in spec)
        return (float(l), float(r), float(b), float(t))

    raise ValueError(f"Unsupported paper_margins spec: {spec!r}")

def parse_content_aspect(spec: Union[str, float, int, Tuple[float, float], None]) -> Optional[float]:
    """Parse content_aspect and return width/height as float.

    Accepts:
      - None
      - "4:3", "1:1"
      - float (already width/height)
      - (w,h)
    """
    if spec is None:
        return None
    if isinstance(spec, (int, float)):
        v = float(spec)
        if not np.isfinite(v) or v <= 0:
            raise ValueError(f"Invalid content_aspect: {spec!r}")
        return v
    if isinstance(spec, (tuple, list)) and len(spec) == 2:
        w = float(spec[0]); h = float(spec[1])
        if not np.isfinite(w) or not np.isfinite(h) or w <= 0 or h <= 0:
            raise ValueError(f"Invalid content_aspect: {spec!r}")
        return w / h
    if isinstance(spec, str):
        s = spec.strip()
        # allow "w:h"
        if ":" in s:
            a,b = s.split(":", 1)
            w = float(a); h = float(b)
            if w <= 0 or h <= 0:
                raise ValueError(f"Invalid content_aspect: {spec!r}")
            return w / h
        # allow "w/h"
        if "/" in s:
            a,b = s.split("/",1)
            w = float(a); h = float(b)
            if w <= 0 or h <= 0:
                raise ValueError(f"Invalid content_aspect: {spec!r}")
            return w / h
        # allow float-like
        v = float(s)
        if v <= 0:
            raise ValueError(f"Invalid content_aspect: {spec!r}")
        return v
    raise ValueError(f"Unsupported content_aspect spec: {spec!r}")

def paper_inner_rect_frac(
    paper_size_in: Tuple[float, float],
    margins_in: Tuple[float, float, float, float],
    *,
    strict: bool = False,
) -> Tuple[float, float, float, float]:
    """Return inner rect (left,bottom,width,height) in figure fraction.

    Parameters
    ----------
    paper_size_in:
        (width_in, height_in) in inches.
    margins_in:
        (left,right,bottom,top) in inches.
    strict:
        If True, raise ValueError when margins make the inner area non-positive.
        If False, keep legacy behavior and clamp to a tiny positive size.
    """
    pw, ph = float(paper_size_in[0]), float(paper_size_in[1])
    l, r, b, t = (float(margins_in[0]), float(margins_in[1]), float(margins_in[2]), float(margins_in[3]))

    pw = builtins.max(1e-9, pw)
    ph = builtins.max(1e-9, ph)

    left = l / pw
    right = 1.0 - (r / pw)
    bottom = b / ph
    top = 1.0 - (t / ph)

    w_raw = right - left
    h_raw = top - bottom

    if strict and (w_raw <= 0.0 or h_raw <= 0.0):
        raise ValueError(
            "Invalid paper_margins: inner area is non-positive. "
            f"paper_size_in=({pw:.6g},{ph:.6g}) "
            f"margins_in(L,R,B,T)=({l:.6g},{r:.6g},{b:.6g},{t:.6g})"
        )

    w = builtins.max(1e-9, w_raw)
    h = builtins.max(1e-9, h_raw)
    return (float(left), float(bottom), float(w), float(h))

def place_single_axes_in_rect(
    ax: Any,
    rect: Tuple[float, float, float, float],
    content_aspect: Optional[float],
) -> Tuple[float, float, float, float]:
    """Place a single Axes inside rect, maximizing size while keeping box aspect (width/height)."""
    x0, y0, w, h = (float(rect[0]), float(rect[1]), float(rect[2]), float(rect[3]))
    if content_aspect is None:
        ax.set_position([x0, y0, w, h])
        return (x0, y0, w, h)

    aspect = float(content_aspect)
    aspect = builtins.max(1e-9, aspect)

    # Max rectangle with given aspect inside (w,h)
    rw = w
    rh = rw / aspect
    if rh > h:
        rh = h
        rw = rh * aspect

    xx = x0 + (w - rw) / 2.0
    yy = y0 + (h - rh) / 2.0
    ax.set_position([xx, yy, rw, rh])
    try:
        ax.set_box_aspect(1.0 / aspect)  # height/width
    except Exception:
        pass
    return (float(xx), float(yy), float(rw), float(rh))

def place_grid_axes_in_rect(
    axes_2d: Any,
    rect: Tuple[float, float, float, float],
    nrows: int,
    ncols: int,
    pad: float = 0.08,
    content_aspect: Optional[float] = None,
) -> None:
    """Manually place a grid of axes inside rect.
    Each cell is centered; if content_aspect is provided, each axes box uses width/height=content_aspect.
    pad is a fraction of each axes size (like box_padding).
    """
    x0, y0, w, h = (float(rect[0]), float(rect[1]), float(rect[2]), float(rect[3]))
    nrows = int(nrows); ncols = int(ncols)
    if nrows <= 0 or ncols <= 0:
        return

    pad = float(pad)
    pad = builtins.max(0.0, pad)

    denom_w = float(ncols) + float(builtins.max(ncols - 1, 0)) * pad
    denom_h = float(nrows) + float(builtins.max(nrows - 1, 0)) * pad
    denom_w = builtins.max(1e-9, denom_w)
    denom_h = builtins.max(1e-9, denom_h)

    if content_aspect is None:
        cell_w = w / denom_w
        cell_h = h / denom_h
    else:
        aspect = builtins.max(1e-9, float(content_aspect))
        # choose cell_h that satisfies both width and height constraints
        cell_h = builtins.min(w / (aspect * denom_w), h / denom_h)
        cell_h = builtins.max(1e-9, float(cell_h))
        cell_w = cell_h * aspect

    gap_w = cell_w * pad
    gap_h = cell_h * pad

    grid_w = (float(ncols) * cell_w) + float(builtins.max(ncols - 1, 0)) * gap_w
    grid_h = (float(nrows) * cell_h) + float(builtins.max(nrows - 1, 0)) * gap_h

    gx0 = x0 + (w - grid_w) / 2.0
    gy0 = y0 + (h - grid_h) / 2.0

    # axes_2d indexing: [row, col] with row 0 at top (matplotlib subplots convention)
    for r in range(nrows):
        for c in range(ncols):
            left = gx0 + float(c) * (cell_w + gap_w)
            bottom = gy0 + float((nrows - 1) - r) * (cell_h + gap_h)
            ax = axes_2d[r, c]
            ax.set_position([left, bottom, cell_w, cell_h])
            if content_aspect is not None:
                try:
                    ax.set_box_aspect(1.0 / float(content_aspect))
                except Exception:
                    pass
            try:
                ax.set_anchor("C")
            except Exception:
                pass

# -----------------------------------------------------------------------------
# PAPER / Layout helpers for figure-level labels (commonized)
# -----------------------------------------------------------------------------
def axes_grid_bbox_frac(
    axes_2d: Any,
    *,
    fallback: Tuple[float, float, float, float] = (0.95, 0.05, 0.10, 0.90),
) -> Tuple[float, float, float, float]:
    """Return (top_y1, bot_y0, left_x0, right_x1) in figure fraction for a grid of Axes.

    axes_2d:
        2D-like container (e.g., np.ndarray) of matplotlib Axes.
    fallback:
        Used when positions cannot be computed for any reason.
    """
    try:
        arr = np.asarray(axes_2d)
        if arr.ndim != 2:
            # Try to reshape if possible
            # (This is best-effort; callers should pass proper 2D.)
            arr = arr.reshape((-1, 1)) if arr.size else arr
        if arr.size == 0:
            return fallback

        top_y1 = -999.0
        bot_y0 = 999.0
        left_x0 = 999.0
        right_x1 = -999.0

        for ax in arr.ravel():
            if ax is None:
                continue
            try:
                pos = ax.get_position()
            except Exception:
                continue
            top_y1 = pos.y1 if pos.y1 > top_y1 else top_y1
            bot_y0 = pos.y0 if pos.y0 < bot_y0 else bot_y0
            left_x0 = pos.x0 if pos.x0 < left_x0 else left_x0
            right_x1 = pos.x1 if pos.x1 > right_x1 else right_x1

        if top_y1 == -999.0 or bot_y0 == 999.0 or left_x0 == 999.0 or right_x1 == -999.0:
            return fallback

        # Clamp into [0,1] just in case
        top_y1 = float(builtins.min(1.0, builtins.max(0.0, top_y1)))
        bot_y0 = float(builtins.min(1.0, builtins.max(0.0, bot_y0)))
        left_x0 = float(builtins.min(1.0, builtins.max(0.0, left_x0)))
        right_x1 = float(builtins.min(1.0, builtins.max(0.0, right_x1)))

        return (top_y1, bot_y0, left_x0, right_x1)
    except Exception:
        return fallback


def adjacent_label_positions_from_bbox(
    *,
    top_row_y1: float,
    bot_row_y0: float,
    left_col_x0: float,
    title_pad: float = 0.07,
    xlabel_pad: float = 0.07,
    ylabel_pad: float = 0.08,
    shrink: float = 0.90,
    min_x: float = 0.01,
) -> Tuple[float, float, float]:
    """Compute figure-level label positions adjacent to the content bbox.

    IMPORTANT:
      - This places labels relative to the *content* (e.g., the grid of axes),
        not relative to the page edges.
      - If the remaining margin is too small, pads are automatically reduced
        to avoid clipping in fixed-paper outputs.

    Returns:
      (title_y, xlabel_y, ylabel_x) in figure fraction.
    """
    # Remaining spaces to figure edges
    top_space = builtins.max(0.0, 1.0 - float(top_row_y1))
    bot_space = builtins.max(0.0, float(bot_row_y0))
    left_space = builtins.max(0.0, float(left_col_x0))

    # Reduce pads if needed (keep a small safety margin using `shrink`)
    title_off = builtins.min(float(title_pad), float(shrink) * top_space) if top_space > 0 else 0.0
    xlabel_off = builtins.min(float(xlabel_pad), float(shrink) * bot_space) if bot_space > 0 else 0.0
    ylabel_off = builtins.min(float(ylabel_pad), float(shrink) * left_space) if left_space > 0 else 0.0

    title_y = float(top_row_y1) + float(title_off)
    xlabel_y = float(bot_row_y0) - float(xlabel_off)
    ylabel_x = builtins.max(float(min_x), float(left_col_x0) - float(ylabel_off))

    # Final clamp
    title_y = float(builtins.min(1.0, builtins.max(0.0, title_y)))
    xlabel_y = float(builtins.min(1.0, builtins.max(0.0, xlabel_y)))
    ylabel_x = float(builtins.min(1.0, builtins.max(0.0, ylabel_x)))

    return title_y, xlabel_y, ylabel_x
