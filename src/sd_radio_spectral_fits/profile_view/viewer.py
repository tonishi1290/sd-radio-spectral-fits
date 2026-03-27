# src/sd_radio_spectral_fits/plotting/viewer.py
from __future__ import annotations

import builtins
from pathlib import Path
from typing import Optional, Union, Tuple, Sequence, List

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from ..fitsio import Scantable, read_scantable
from ..ranges import parse_windows
# [MODIFIED] Import helpers from scantable_utils
from ..scantable_utils import _parse_row_selector, _df_to_native_endian
from .utils import (
AxisBundle, _norm_range, _parse_list_like, _interp_extrap, _freq_ghz_from_vel,
    _convert_user_xrange_to_vel_range, _convert_vel_range_to_xrange, _build_fitwin_mask_on_vel,
    _process_spectrum, _baseline_extrapolated_line, _rms_from_mean,
    _running_mean, _box_average_downsample, recalculate_velocity_axis,
    drive_pdf_generation, # [ADDED]
    parse_figsize, # [ADDED]
    parse_paper_margins, # [ADDED PAPER]
    parse_content_aspect, # [ADDED PAPER]
    paper_inner_rect_frac, # [ADDED PAPER]
    place_single_axes_in_rect, # [ADDED PAPER]
)
from .windowmask import fit_windows_xaxis, compute_rms_win
# [MODIFIED] Use robust axis generation
from ..regrid_vlsrk import vlsrk_axis_for_spectrum, _generate_obs_freq_axis
# [MODIFIED] Import tempscale helpers
from ..tempscale import (
    normalize_tempscal,
    beameff_array,
    ta_to_tr,
    tr_to_ta
)



def _row_has_display_baseline(row: pd.Series) -> bool:
    applied = row.get("BSL_APPLIED", False)
    try:
        if pd.isna(applied):
            applied = False
    except Exception:
        pass
    if not bool(applied):
        return False
    arr = _parse_list_like(row.get("BSL_COEF")) if "BSL_COEF" in row.index else None
    return arr is not None

# -----------------------------------------------------------------------------
# Row Selection Helper
# -----------------------------------------------------------------------------
def _filter_scantable_by_rows(st: Scantable, rows=None, exclude_rows=None) -> Scantable:
    """
    Apply row selection/exclusion logic with endian safety and synchronization.
    Enforces unique, sorted indices.
    """
    if rows is None and exclude_rows is None:
        return st
        
    if rows is not None and exclude_rows is not None:
        raise ValueError("Cannot specify both 'rows' and 'exclude_rows'.")
        
    n_total = len(st.table)
    
    if rows is not None:
        # Parse and normalize (unique + sorted)
        raw_idx = _parse_row_selector(rows, n_total)
        final_idx = np.unique(raw_idx)
    else:
        # Exclude logic
        exclude_idx = _parse_row_selector(exclude_rows, n_total)
        all_idx = np.arange(n_total, dtype=int)
        final_idx = np.setdiff1d(all_idx, exclude_idx) # setdiff1d returns unique & sorted
        
    if len(final_idx) == 0:
        raise ValueError(f"Row selection resulted in empty Scantable. (rows={rows}, exclude_rows={exclude_rows})")
    
    # Check bounds
    if len(final_idx) > 0 and (final_idx.min() < 0 or final_idx.max() >= n_total):
         raise ValueError(f"Row selector contains indices out of bounds [0, {n_total-1}].")
         
    # Apply to Table (Safe Endian fix for FITS compatibility)
    safe_table = _df_to_native_endian(st.table)
    
    # [FIX] Preserve original global row index before reset
    new_table = safe_table.iloc[final_idx].copy()
    if "ORIG_ROW_ID" not in new_table.columns:
        new_table["ORIG_ROW_ID"] = final_idx
    new_table.reset_index(drop=True, inplace=True)
    
    # Apply to Data (Sync with table)
    if isinstance(st.data, list):
        new_data = [st.data[i] for i in final_idx]
    else:
        new_data = st.data[final_idx]
        
    # Metadata
    new_meta = st.meta.copy() if hasattr(st.meta, 'copy') else dict(st.meta)
    hist_key = 'HISTORY_ROWS' if rows is not None else 'HISTORY_EXCLUDE_ROWS'
    hist_val = str(rows) if rows is not None else str(exclude_rows)
    new_meta[hist_key] = hist_val

    # [FIX] History (Append to list to avoid overwrite)
    new_hist = st.history.copy() if hasattr(st, "history") and st.history else {}
    if "selection_log" not in new_hist:
        new_hist["selection_log"] = []
    
    if not isinstance(new_hist["selection_log"], list):
         new_hist["selection_log"] = [str(new_hist["selection_log"])]
         
    new_hist["selection_log"].append(f"Select: {hist_key}={hist_val}")

    return Scantable(meta=new_meta, table=new_table, data=new_data, history=new_hist)


def _merge_scantables(inputs: Union[Scantable, str, Sequence[Union[Scantable, str]]]) -> Scantable:
    """Helper to merge multiple inputs into a single VLA-style Scantable."""
    if isinstance(inputs, (Scantable, str)):
        inputs = [inputs]
    
    sc_list = []
    for inp in inputs:
        if isinstance(inp, str):
            sc_list.append(read_scantable(inp))
        else:
            sc_list.append(inp)

    if not sc_list:
        raise ValueError("No input data provided.")
    
    if len(sc_list) == 1:
        return sc_list[0]

    # Merge
    # Use meta from first
    meta = sc_list[0].meta
    
    tabs = []
    data_all = []
    
    for sc in sc_list:
        # Table
        tabs.append(_df_to_native_endian(sc.table))
        
        # Data (Convert to list of arrays to support heterogeneous lengths)
        d = sc.data
        if isinstance(d, list):
            data_all.extend(d)
        else:
            # 2D array -> list of arrays
            data_all.extend(list(d))
            
    table_all = pd.concat(tabs, axis=0, ignore_index=True)
    
    # Init history for merged
    merged_hist = {"merged_count": len(sc_list), "sources": [str(x) for x in inputs]}
    
    # Return new Scantable with list-based data (VLA ready)
    return Scantable(meta=meta, data=data_all, table=table_all, history=merged_hist)

def _row_specsys(row: Optional[pd.Series], meta: dict) -> str:
    """Return effective spectral reference frame from row/header."""
    for key in ("SPECSYS", "SSYSOBS"):
        val = None
        if row is not None:
            try:
                if key in row.index:
                    val = row.get(key)
            except Exception:
                val = None
        if val is None and isinstance(meta, dict):
            val = meta.get(key)
        if val is None:
            continue
        s = str(val).strip().upper()
        if s and s not in ("NAN", "NONE"):
            return s
    return ""


def _row_meta_for_axis(row: Optional[pd.Series], meta: dict) -> dict:
    """Build effective per-row spectral metadata for axis generation."""
    merged = meta.copy() if isinstance(meta, dict) else dict(meta)
    if row is None:
        return merged
    for key in (
        "CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1",
        "RESTFRQ", "RESTFREQ", "SPECSYS", "SSYSOBS"
    ):
        try:
            if key in row.index:
                val = row.get(key)
                if val is not None and not pd.isna(val):
                    merged[key] = val
        except Exception:
            pass
    return merged


def _row_rest_hz_meta(row: Optional[pd.Series], meta: dict) -> Optional[float]:
    eff = _row_meta_for_axis(row, meta)
    for key in ("RESTFRQ", "RESTFREQ"):
        val = eff.get(key)
        try:
            f = float(val)
        except Exception:
            continue
        if np.isfinite(f) and f > 0:
            return float(f)
    return None


def _row_vcorr_kms(row: Optional[pd.Series], meta: dict) -> float:
    """
    Resolve per-row LSRK correction in km/s.

    New spec:
      - VELOSYS / VFRAME : m/s
      - V_CORR_KMS / v_corr_kms : km/s (legacy)
    If SPECSYS/SSYSOBS explicitly indicates a non-TOPO frame, no additional
    correction is applied. If the frame is unknown, legacy velocity columns are
    still honored for backward compatibility in plotting.
    """
    specsys = _row_specsys(row, meta)
    if specsys and ("TOPO" not in specsys):
        return 0.0

    def _get_num(key: str):
        val = None
        if row is not None:
            try:
                if key in row.index:
                    val = row.get(key)
            except Exception:
                val = None
        if val is None and isinstance(meta, dict):
            val = meta.get(key)
        try:
            f = float(val)
        except Exception:
            return np.nan
        return f if np.isfinite(f) else np.nan

    for key in ("VELOSYS", "VFRAME"):
        val = _get_num(key)
        if np.isfinite(val):
            return float(val) * 1.0e-3

    for key in ("V_CORR_KMS", "v_corr_kms"):
        val = _get_num(key)
        if np.isfinite(val):
            return float(val)

    return 0.0





def _apply_freq_axis_format(ax, labelsize: Optional[str] = None) -> None:
    """Render GHz ticks plainly without scientific offset text."""
    try:
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, pos: f"{x:.3f}"))
        ax.xaxis.offsetText.set_visible(False)
    except Exception:
        pass
    if labelsize is not None:
        try:
            ax.tick_params(axis="x", labelsize=labelsize)
        except Exception:
            pass

class SpectralViewer:
    """
    Interactive Viewer for Single Dish Spectra.
    Supports both Fixed-Length (2D Array) and VLA (List of Arrays) Scantables.
    
    Parameters
    ----------
    scantable : Union[Scantable, str, Sequence[Union[Scantable, str]]]
        Input data.
    rows : Union[int, slice, List[int], str], optional
        Select specific rows (global indices) to display.
        Cannot be used with exclude_rows.
        Indices are sorted and duplicates removed.
    exclude_rows : Union[int, slice, List[int], str], optional
        Exclude specific rows (global indices) from display.
        Cannot be used with rows.
    """
    def __init__(
        self,
        scantable: Union[Scantable, str, Sequence[Union[Scantable, str]]],
        rows: Optional[Union[int, slice, List[int], str]] = None,
        exclude_rows: Optional[Union[int, slice, List[int], str]] = None,
        xaxis: str = "vel",
        axis_type: Optional[str] = None,
        rest_freq: Optional[float] = None,
        xrange: Optional[Tuple[float, float]] = None,
        yrange: Optional[Tuple[float, float]] = None,
        autoscale_y: bool = True,
        show_fitwin_rms: bool = True,
        show_fit_windows: bool = True,
        smooth_mode: str = "none",
        smooth_width: int = 1,
        box_downsample: bool = False,
        box_policy: str = "trim",
        rms_on: str = "raw",
        show_top_axis: bool = True,
        save_dir: Union[str, Path] = ".",
        save_prefix: str = "spectrum",
        save_pdf: Optional[Union[str, Path]] = None, # [ADDED]
        max_pdf_pages: Optional[int] = 100, # [ADDED]
        show: bool = True,
        figsize: Optional[Union[str, Tuple[float, float]]] = None, # [ADDED]
        content_aspect: Optional[Union[str, float, Tuple[float, float]]] = None, # [ADDED PAPER]
        paper_margins: Optional[Union[str, Tuple[float, float, float, float], List[float]]] = None, # [ADDED PAPER]
    ):
        if axis_type is not None and (xaxis is None or str(xaxis).lower() == "vel"):
            xaxis = axis_type

        # [MODIFIED] Handle multiple inputs / files and Row Selection
        merged_st = _merge_scantables(scantable)
        self.st = _filter_scantable_by_rows(merged_st, rows=rows, exclude_rows=exclude_rows)
        
        self.idx = 0
        self.xaxis = str(xaxis).lower()
        if self.xaxis not in ("vel", "freq", "chan"):
            self.xaxis = "vel"

        # VLA Handling: Check if data is list or array
        self.is_vla = isinstance(self.st.data, list)
        self.n_spec = len(self.st.table)
        
        # [MODIFIED] Check first row to set default display (but not strictly bound to it)
        tempscal_col = self.st.table["TEMPSCAL"].iloc[0] if "TEMPSCAL" in self.st.table.columns else "TA*"
        self.default_scale = normalize_tempscal(tempscal_col)
        self.display_scale = self.default_scale

        # Meta defaults
        self.rest_hz_user = float(rest_freq) if rest_freq is not None else None
        
        # Get Header Rest Freq (Fallback if user doesn't specify, or needed for conversion)
        val = self.st.meta.get("RESTFREQ") or self.st.meta.get("RESTFRQ")
        self.rest_hz_meta = float(val) if val is not None else None

        # Effective Rest Freq for Display
        if self.rest_hz_user is not None and self.rest_hz_user > 0:
            self.rest_hz = self.rest_hz_user
        else:
            self.rest_hz = self.rest_hz_meta

        self.xrange_user = _norm_range(xrange)
        self.yrange = _norm_range(yrange)
        self.autoscale_y = bool(autoscale_y)
        self.show_fitwin_rms = bool(show_fitwin_rms)
        self.show_fit_windows = bool(show_fit_windows)
        self.smooth_mode = str(smooth_mode).lower()
        self.smooth_width = int(smooth_width)
        self.box_downsample = bool(box_downsample)
        self.box_policy = str(box_policy).lower()
        self.rms_on = str(rms_on).lower()
        self.show_top_axis = bool(show_top_axis)
        self.save_dir = Path(save_dir)
        self.save_prefix = str(save_prefix)

        # Note: AxisBundle from utils is mostly for fixed length.
        # For VLA, we will regenerate axis per row.
        # However, we keep a "default" nchan for initial fallback
        if not self.is_vla:
            self.default_nchan = int(self.st.data.shape[1])
        else:
            self.default_nchan = len(self.st.data[0]) if len(self.st.data) > 0 else 0
            
        self.view_mode = "original"

        self.fig = None
        self.ax = None
        self._ax_top = None
        self._warned_topo = False
        self.xr_vel = None # Will be set on first plot

        # Pre-calc initial velocity range if possible (using first row)
        if self.xrange_user is not None and self.n_spec > 0:
            self.xr_vel = self._user_xrange_to_vel(self.xaxis, self.xrange_user)

        # [MODIFIED] Determine Figure Size
        default_figsize = (10, 6)
        fig_size, bbox_mode = parse_figsize(figsize, default=default_figsize)
        # [ADDED PAPER] Paper mode settings
        self._paper_mode = (bbox_mode is None and fig_size is not None)
        self._paper_size_in = fig_size if self._paper_mode else None
        self._paper_margins_in = parse_paper_margins(paper_margins) if self._paper_mode else None
        self._content_aspect = parse_content_aspect(content_aspect) if self._paper_mode else None
        self._paper_title_y = None

        # [MODIFIED] Create Figure
        if show or save_pdf:
            self.fig, self.ax = plt.subplots(figsize=fig_size)
            # [ADDED PAPER] Layout policy
            if self._paper_mode and self._paper_size_in is not None and self._paper_margins_in is not None:
                # Fixed paper size: place axes inside margins; (A) box aspect via content_aspect
                inner = paper_inner_rect_frac(self._paper_size_in, self._paper_margins_in, strict=True)
                place_single_axes_in_rect(self.ax, inner, self._content_aspect)
                # Figure-level title position (center of top margin)
                pw, ph = float(self._paper_size_in[0]), float(self._paper_size_in[1])
                top_in = float(self._paper_margins_in[3])
                self._paper_title_y = 1.0 - (top_in / builtins.max(1e-9, ph)) * 0.5
            else:
                # Screen optimized
                self.fig.subplots_adjust(top=0.85, bottom=0.15, left=0.10, right=0.90)

            self.cid = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
            self._print_help()
            
            # Initial plot
            self.update_plot()

            # [MODIFIED] PDF Generation
            if save_pdf:
                def _updater(i):
                    self.idx = i
                    self.update_plot()
                
                # [MODIFIED] Enforce tight layout with padding for PDF
                drive_pdf_generation(
                    save_pdf,
                    self.fig,
                    _updater,
                    count=self.n_spec,
                    max_pages=max_pdf_pages,
                    bbox_inches=(None if self._paper_mode else "tight"),
                    pad_inches=(0.0 if self._paper_mode else 0.75),
                )
                
                # [FIX] State Reset
                self.idx = 0
                self.update_plot()
                if not show:
                    plt.close(self.fig)
                    self.fig = None

            if show and self.fig is not None:
                plt.show()

    def _print_help(self):
        print("\n--- Controls ---")
        print(" n / p : Next / Previous Spectrum")
        print(" b     : Toggle Baseline View (Original <-> Baseline Added)")
        print(" x     : Cycle x-axis")
        print(" t     : Toggle Temp Scale (TA* <-> TR*)")
        print(" m     : Cycle smoothing")
        print(" [ / ] : Change smooth width")
        print(" d     : Toggle downsample")
        print(" r / R : Toggle rms_on (raw <-> display)")
        print(" T     : Toggle top axis")
        print(" s     : Save PDF")
        print(" q     : Quit")
        print("----------------")

    def _cycle_xaxis(self):
        order = ["vel", "freq", "chan"]
        i = order.index(self.xaxis) if self.xaxis in order else 0
        self.xaxis = order[(i + 1) % len(order)]
        
    def _toggle_scale(self):
        if self.display_scale == "TA*":
            # [MODIFIED] Check if BEAMEFF is available before toggling to TR*
            idx = int(self.idx) % self.n_spec
            row = self.st.table.iloc[idx]
            beameff = float(row.get("BEAMEFF", np.nan))
            
            # Fallback to metadata check if row has no BEAMEFF
            if not (np.isfinite(beameff) and beameff > 0):
                 val = self.st.meta.get("BEAMEFF")
                 if val is not None:
                     try: beameff = float(val)
                     except: pass
            
            if not (np.isfinite(beameff) and beameff > 0):
                print("⚠️ Cannot toggle to TR*: BEAMEFF is missing or invalid for this spectrum.")
                return

            self.display_scale = "TR*"
        else:
            self.display_scale = "TA*"
        print(f"Scale toggled to: {self.display_scale}")

    def _get_data(self, idx: int) -> np.ndarray:
        """Retrieve spectrum data safely (VLA aware)."""
        if self.is_vla:
            return np.asarray(self.st.data[idx], dtype=float)
        else:
            return np.asarray(self.st.data[idx], dtype=float)

    def _make_top_axis(self, vel: np.ndarray, chan: np.ndarray, freq_ghz: Optional[np.ndarray], xaxis: str):
        if not self.show_top_axis:
            return None
        xa = str(xaxis).lower()
        ax_top = self.ax.twiny()
        # [追加] メイン軸と位置・サイズを完全同期させる
        ax_top.set_position(self.ax.get_position())
        
        ax_top.set_xlim(self.ax.get_xlim())

        # Simplified top axis logic:
        x0, x1 = self.ax.get_xlim()
        
        # Define mapping function based on current X axis
        if xa == "chan":
            # Top axis: Velocity (if available) else Freq
            if vel is not None:
                def x_from_bottom(v_bottom): 
                    return _interp_extrap(chan, vel, np.asarray(v_bottom, dtype=float))
                label = "Velocity [km/s]"
            else:
                return None
        elif xa == "vel":
            # Top axis: Channel
            def x_from_bottom(v_bottom):
                return _interp_extrap(vel, chan, np.asarray(v_bottom, dtype=float))
            label = "Channel"
        else: # freq
            # Top axis: Velocity
            def x_from_bottom(v_bottom):
                if freq_ghz is None: return v_bottom
                return _interp_extrap(freq_ghz, vel, np.asarray(v_bottom, dtype=float))
            label = "Velocity [km/s]"

        # Create ticks
        # limits map
        l_bot, r_bot = self.ax.get_xlim()
        l_top = x_from_bottom(np.array([l_bot]))[0]
        r_top = x_from_bottom(np.array([r_bot]))[0]
        
        ax_top.set_xlim(l_top, r_top)
        ax_top.set_xlabel(label, fontsize="small")
        return ax_top

    def _row_vel_lsrk_and_oriented(self, row: pd.Series, y: np.ndarray):
        """
        Generate physical axes for the current row.
        Handles VLA (variable nchan).
        Includes Logic to RECALCULATE velocity if user rest_freq != header rest_freq.
        """
        nchan = len(y)
        
        # 1. Get Velocity Axis (LSRK) using robust helper.
        #    Row-level WCS / rest-frequency override global metadata when present.
        eff_meta = _row_meta_for_axis(row, self.st.meta)
        row_rest_hz_meta = _row_rest_hz_meta(row, self.st.meta)
        v_corr = _row_vcorr_kms(row, self.st.meta)
        vel_native = vlsrk_axis_for_spectrum(eff_meta, v_corr_kms=float(v_corr), nchan=nchan)
        
        # 2. Check for Rest Frequency Override (Shift)
        # Only do this if we have a valid header restfreq AND a valid user restfreq
        vel_disp = vel_native
        if (row_rest_hz_meta is not None and row_rest_hz_meta > 0 and 
            self.rest_hz_user is not None and self.rest_hz_user > 0 and
            abs(row_rest_hz_meta - self.rest_hz_user) > 1.0):
            
            # Recalculate: V_disp = f(V_native, Rest_old, Rest_new)
            vel_disp = recalculate_velocity_axis(vel_native, row_rest_hz_meta, self.rest_hz_user)

        # 3. Get Frequency Axis (LSRK) for display
        # Always use the EFFECTIVE rest frequency (self.rest_hz)
        if self.rest_hz is not None and self.rest_hz > 0:
            c_kms = 299792.458
            freq_ghz = (self.rest_hz * (1.0 - vel_disp / c_kms)) * 1e-9
        else:
            # Try to get obs freq (approx) from WCS
            try:
                f_obs = _generate_obs_freq_axis(eff_meta, nchan)
                freq_ghz = f_obs * 1e-9
            except:
                freq_ghz = None

        chan = np.arange(nchan, dtype=float)
        
        # Orient (ensure increasing velocity if preferred, or standard logic)
        if vel_disp.size >= 2 and vel_disp[0] > vel_disp[-1]:
            # Flip to increasing velocity
            vel_disp = vel_disp[::-1]
            chan = chan[::-1]
            y = y[::-1]
            if freq_ghz is not None: freq_ghz = freq_ghz[::-1]
            
        return vel_disp, chan, freq_ghz, v_corr, y

    def _user_xrange_to_vel(self, initial_xaxis: str, user_xrange: Tuple[float, float]) -> Tuple[float, float]:
        # Use current index to resolve
        y0 = self._get_data(int(self.idx))
        row0 = self.st.table.iloc[int(self.idx)]
        vel, chan, freq_ghz, v_corr, _ = self._row_vel_lsrk_and_oriented(row0, y0)
        
        v0, v1 = _convert_user_xrange_to_vel_range(
            initial_xaxis, user_xrange, 
            vel, chan, freq_ghz, self.rest_hz, v_corr
        )
        return (v0, v1)

    def _xlim_for_current_axis(self, vel: np.ndarray, chan: np.ndarray, freq_ghz: Optional[np.ndarray]) -> Optional[Tuple[float, float]]:
        if self.xr_vel is None: return None
        return _convert_vel_range_to_xrange(self.xaxis, self.xr_vel, vel, chan, freq_ghz, self.rest_hz)

    def update_plot(self):
        if self.fig is None or self.ax is None: return
        self.ax.clear()
        
        # [MODIFIED] Explicitly remove twin axis to prevent ghosting
        if self._ax_top is not None:
            try:
                self._ax_top.remove()
            except Exception:
                pass
            self._ax_top = None

        if self.n_spec <= 0: return
        idx = int(self.idx) % self.n_spec
        row = self.st.table.iloc[idx]
        
        # Get Data (VLA safe)
        current_data = self._get_data(idx)
        
        # [MODIFIED] Get Per-Row Scale (instead of assuming homogeneous)
        row_scale = normalize_tempscal(row.get("TEMPSCAL", "TA*"))
        
        # [MODIFIED] Scale Conversion (Non-destructive / On-the-fly)
        if self.display_scale != row_scale:
            beameff = float(row.get("BEAMEFF", np.nan))
            # Fallback to header BEAMEFF if row BEAMEFF is invalid
            if not (np.isfinite(beameff) and beameff > 0):
                 val = self.st.meta.get("BEAMEFF")
                 if val is not None:
                     try: beameff = float(val)
                     except: pass
            
            # Prepare arrays for conversion (handle broadcasting properly)
            # Reshape 1D spectrum to (1, N) to ensure ta_to_tr returns correct shape
            # This fixes "TypeError: object of type 'numpy.float32' has no len()"
            data_2d = current_data.reshape(1, -1)
            eff_arr = np.array([beameff])
            
            # Apply conversion if BEAMEFF is valid
            if np.isfinite(beameff) and beameff > 0:
                if self.display_scale == "TR*" and row_scale == "TA*":
                    # ta_to_tr returns (1, N), [0] takes the first row -> (N,)
                    current_data = ta_to_tr(data_2d, eff_arr)[0]
                elif self.display_scale == "TA*" and row_scale == "TR*":
                    current_data = tr_to_ta(data_2d, eff_arr)[0]

        # Generate Axes (Includes Velocity Shift Logic)
        vel, chan, freq_ghz, v_corr, spectrum = self._row_vel_lsrk_and_oriented(row, current_data)
        
        # Setup X axis for plot
        if self.xaxis == "vel": 
            x_raw = vel; xlabel = "Velocity (LSRK) [km/s]"
        elif self.xaxis == "freq" and freq_ghz is not None:
            x_raw = freq_ghz; xlabel = "Frequency [GHz]"
        else: 
            x_raw = chan; xlabel = "Channel"

        if self.xr_vel is None: self.xr_vel = (float(np.min(vel)), float(np.max(vel)))
        xlim = self._xlim_for_current_axis(vel, chan, freq_ghz)

        # Baseline Logic
        has_bsl = _row_has_display_baseline(row)
        win_str = str(row.get("BSL_WINF", "")) if has_bsl else ""
        arr = _parse_list_like(row.get("BSL_COEF")) if has_bsl and "BSL_COEF" in row.index else None
        
        if self.view_mode == "baseline_added" and not has_bsl:
            self.view_mode = "original"

        baseline_val = None
        reconstructed_raw = None
        
        if has_bsl:
            if arr.size == spectrum.size:
                # Vector baseline (Legacy)
                base_raw = arr.astype(float, copy=True)
                # Check if we flipped in _row_vel...
                if chan[0] > chan[-1]: 
                     base_raw = base_raw[::-1]
                baseline_val = base_raw
            else:
                # Poly coefficients (New)
                co = arr.astype(float, copy=False)
                try:
                    baseline_val = np.polyval(co, vel).astype(float, copy=False)
                except Exception:
                    has_bsl = False

            if baseline_val is not None:
                # [MODIFIED] Apply scale conversion to baseline too
                # Must respect the BSL_SCALE if available, else assume same as data input scale
                bsl_scale = normalize_tempscal(row.get("BSL_SCALE", row_scale))
                
                if self.display_scale != bsl_scale:
                    beameff = float(row.get("BEAMEFF", np.nan))
                    # Use header fallback if needed
                    if not (np.isfinite(beameff) and beameff > 0):
                         val = self.st.meta.get("BEAMEFF")
                         if val is not None:
                             try: beameff = float(val)
                             except: pass

                    if np.isfinite(beameff) and beameff > 0:
                        # Use same reshaping logic for robustness
                        base_2d = baseline_val.reshape(1, -1)
                        eff_arr = np.array([beameff])
                        
                        if self.display_scale == "TR*" and bsl_scale == "TA*":
                             baseline_val = ta_to_tr(base_2d, eff_arr)[0]
                        elif self.display_scale == "TA*" and bsl_scale == "TR*":
                             baseline_val = tr_to_ta(base_2d, eff_arr)[0]
                         
                reconstructed_raw = spectrum + baseline_val

        # Choose Data
        if self.view_mode == "baseline_added" and reconstructed_raw is not None:
            y_plot = reconstructed_raw
            label_main = "Raw Data (Restored)"
        else:
            y_plot = spectrum
            label_main = "Original (Clean)"

        # Smooth & Plot
        x_disp, y_disp = _process_spectrum(x_raw, y_plot, self.smooth_mode, self.smooth_width, self.box_downsample, self.box_policy)
        self.ax.step(x_disp, y_disp, where="mid", color="black", lw=0.8, label=label_main)
        
        # Plot Baseline Curve (if Added mode) or Zero Line
        if self.view_mode == "baseline_added" and baseline_val is not None:
             if has_bsl:
                 x_bl, y_bl = _process_spectrum(x_raw, baseline_val, self.smooth_mode, self.smooth_width, self.box_downsample, self.box_policy)
                 self.ax.plot(x_bl, y_bl, color="red", lw=1.2, alpha=0.85, label="Baseline Model")
        else:
            self.ax.axhline(0, color="gray", lw=0.8, ls="--")

        # Fit Windows
        if self.show_fit_windows and win_str and str(win_str).lower() not in ("nan", "none", ""):
            try:
                wins_x = fit_windows_xaxis(
                    win_str,
                    self.rest_hz_meta,
                    self.rest_hz_user,
                    self.xaxis,
                    vel,
                    freq_ghz,
                    chan,
                )
                for x0, x1 in wins_x:
                    self.ax.axvspan(float(x0), float(x1), color="green", alpha=0.10)
            except Exception:
                pass



        if xlim is not None: self.ax.set_xlim(*xlim)

                # RMS & Labels
        rms_bsl = row.get("BSL_RMS", np.nan)

        rms_win, n_win, _ = compute_rms_win(
            x_raw=x_raw,
            y_resid=spectrum,
            rms_on=self.rms_on,
            smooth_mode=self.smooth_mode,
            smooth_width=self.smooth_width,
            box_downsample=self.box_downsample,
            box_policy=self.box_policy,
            xaxis=self.xaxis,
            vel_disp=vel,
            freq_ghz=freq_ghz,
            chan=chan,
            xr_vel=self.xr_vel,
            win_str=win_str,
            rest_hz_meta=self.rest_hz_meta,
            rest_hz_user=self.rest_hz_user,
            has_bsl=bool(has_bsl),
            use_fit_windows=True,
        )

        smooth_tag = f"{self.smooth_mode}(w={int(self.smooth_width)})"
        if self.smooth_mode in ("boxcar", "box", "avg", "average") and self.box_downsample: smooth_tag += ",downsample"
        
        mode_label = "Baseline Added" if self.view_mode == "baseline_added" else "Original"
        line1 = f"Index: {idx} | Mode: {mode_label} | Scale: {self.display_scale} (Native: {row_scale})"
        line2 = f"Smooth: {smooth_tag} | rms_on: {self.rms_on}"
        if np.isfinite(rms_bsl): line2 += f" | RMS(bsl): {float(rms_bsl):.4f} K"
        if self.show_fitwin_rms: line2 += f" | RMS(win): {rms_win:.4f} K"
        
        # [PAPER] Title placement
        if getattr(self, "_paper_mode", False) and getattr(self, "_paper_title_y", None) is not None:
            # Use figure-level title inside the top margin
            #self.fig.suptitle(line1 + "\n" + line2, fontsize=10, y=float(self._paper_title_y))
            #self.ax.set_title("")
            # [変更後] 常にAxes基準でタイトルを表示し、padで距離を調整する
            pad_title = 35 if self.show_top_axis else 10
            self.ax.set_title(line1 + "\n" + line2, fontsize=10, pad=pad_title)
        else:
            # Screen default
            self.ax.set_title(line1 + "\n" + line2, fontsize=10, y=1.1)
            

        
        self.ax.set_xlabel(xlabel)
        # [MODIFIED] Dynamic Y-label
        self.ax.set_ylabel(f"Intensity ({self.display_scale}) [{self.st.meta.get('BUNIT', 'K')}]")

        # Auto Scale
        if self.yrange is not None:
            self.ax.set_ylim(*self.yrange)
        elif self.autoscale_y:
            l, r = self.ax.get_xlim()
            mx = (x_disp >= min(l,r)) & (x_disp <= max(l,r))
            yy = y_disp[mx]
            yy = yy[np.isfinite(yy)]
            if yy.size:
                ymin, ymax = float(np.min(yy)), float(np.max(yy))
                yr = (ymax - ymin) or 1.0
                self.ax.set_ylim(ymin - 0.1 * yr, ymax + 0.1 * yr)

        self.ax.legend(loc="upper right", fontsize="small")
        self.ax.grid(True, linestyle=":", alpha=0.5)
        if self.xaxis == "freq":
            _apply_freq_axis_format(self.ax)
        
        # [MODIFIED] Correctly assign the new top axis to self._ax_top
        if self.show_top_axis:
            self._ax_top = self._make_top_axis(vel=vel, chan=chan, freq_ghz=freq_ghz, xaxis=self.xaxis)
            
        if hasattr(self.fig.canvas, "draw_idle"):
            self.fig.canvas.draw_idle()
        else:
            self.fig.canvas.draw()

    def on_key(self, event):
        if event is None or event.key is None: return
        key = event.key
        
        if key == "n": self.idx = (self.idx + 1) % self.n_spec; self.update_plot()
        elif key == "p": self.idx = (self.idx - 1 + self.n_spec) % self.n_spec; self.update_plot()
        elif key == "b":
            row = self.st.table.iloc[int(self.idx)]
            has_bsl = _row_has_display_baseline(row)
            if not has_bsl:
                print("⚠️ BSL_APPLIED=True の baseline 情報が無いため Baseline View は切替できません。")
                self.view_mode = "original"
            else:
                self.view_mode = "baseline_added" if self.view_mode == "original" else "original"
            self.update_plot()
        elif key == "x": self._cycle_xaxis(); self.update_plot()
        elif key == "t": self._toggle_scale(); self.update_plot()
        elif key == "m":
            # [MODIFIED] Removed "running" from cycle
            order = ["none", "boxcar"]
            i = order.index(self.smooth_mode) if self.smooth_mode in order else 0
            self.smooth_mode = order[(i + 1) % len(order)]
            self.update_plot()
        elif key == "[": self.smooth_width = builtins.max(1, int(self.smooth_width) - 1); self.update_plot()
        elif key == "]": self.smooth_width = int(self.smooth_width) + 1; self.update_plot()
        elif key == "d": self.box_downsample = not self.box_downsample; self.update_plot()
        elif key in ("r", "R"): self.rms_on = "display" if self.rms_on == "raw" else "raw"; self.update_plot()
        elif key == "T": self.show_top_axis = not self.show_top_axis; self.update_plot()
        elif key == "s":
            save_dir = Path(getattr(self, "save_dir", "."))
            save_dir.mkdir(parents=True, exist_ok=True)
            prefix = getattr(self, "save_prefix", "spectrum")
            # [FIXED] Corrected f-string syntax (added missing opening brace)
            fn = f"{prefix}_{int(self.idx):04d}_{self.view_mode}_{self.xaxis}_{self.display_scale[:2]}.pdf"
            out = save_dir / fn
            self.fig.savefig(out, bbox_inches="tight", pad_inches=0.5)
            print(f"Saved: {out}")
        elif key == "q": plt.close(self.fig)
        
def view_spectra(input_data: Union[Scantable, str, Sequence[Union[Scantable, str]]], **kwargs):
    return SpectralViewer(input_data, **kwargs)

