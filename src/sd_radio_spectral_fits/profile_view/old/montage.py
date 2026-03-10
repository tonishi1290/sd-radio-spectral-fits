# src/sd_radio_spectral_fits/plotting/montage.py
from __future__ import annotations

import builtins
from pathlib import Path
from typing import Optional, Tuple, Union, Sequence, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..fitsio import Scantable, read_scantable
# [MODIFIED] Added vlsrk_axis_for_spectrum for per-row robustness
from ..regrid_vlsrk import vlsrk_axis_for_spectrum
from ..axis import freq_axis_from_wcs, radio_velocity_kms
from ..ranges import parse_windows
# [MODIFIED] Import helpers from scantable_utils
from ..scantable_utils import _parse_row_selector, _df_to_native_endian
from .utils import (
_norm_range, _process_spectrum, _rms_from_mean,
    _parse_list_like, _running_mean, recalculate_velocity_axis,
    _interp_extrap, drive_pdf_generation, # [ADDED]
    parse_figsize, # [ADDED]
    parse_paper_margins, # [ADDED PAPER]
    parse_content_aspect, # [ADDED PAPER]
    paper_inner_rect_frac, # [ADDED PAPER]
    place_grid_axes_in_rect, # [ADDED PAPER]
    axes_grid_bbox_frac, # [PAPER COMMON]
    adjacent_label_positions_from_bbox, # [PAPER COMMON]
)
from .windowmask import fit_windows_xaxis
from ..tempscale import (
    normalize_tempscal,
    beameff_array,
    tempscal_array,
    ta_to_tr,
    tr_to_ta
)

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
    meta = sc_list[0].meta
    
    tabs = []
    data_to_merge = []
    
    for sc in sc_list:
        tabs.append(_df_to_native_endian(sc.table))
        data_to_merge.append(sc.data)
            
    table_all = pd.concat(tabs, axis=0, ignore_index=True)
    
    # --- オリジナルのシンプルなリスト化処理 ---
    data_all = []
    for d in data_to_merge:
        if isinstance(d, list):
            data_all.extend(d)
        else:
            data_all.extend(list(d))
    
    # Init history for merged
    merged_hist = {"merged_count": len(sc_list), "sources": [str(x) for x in inputs]}
    
    return Scantable(meta=meta, data=data_all, table=table_all, history=merged_hist)    

class ProfileMapMontageViewer:
    """
    Robust Montage Viewer.
    Uses 'Hybrid' axis generation:
      - If data is regular 2D array, calculates axis PER ROW (like Viewer) to handle TOPO/v_corr.
      - If data is VLA (ragged), uses Standardizer to regrid to common axis.
    
    Parameters
    ----------
    scantable : Union[Scantable, str, Sequence[Union[Scantable, str]]]
        Input data.
    rows : Union[int, slice, List[int], str], optional
        Select specific rows (global indices) to display.
        Examples: 0, "0,2,5", "10:20". Cannot be used with exclude_rows.
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
        imin: int = 0,
        imax: Optional[int] = None,
        nrows: int = 5,
        ncols: int = 4,
        xrange: Optional[Tuple[float, float]] = None,
        yrange: Optional[Tuple[float, float]] = None,
        autoscale_y: bool = True,
        annotate_rms: bool = True,
        show_fit_windows: bool = True,
        smooth_mode: str = "none",
        smooth_width: int = 1,
        box_downsample: bool = False,
        box_policy: str = "trim",
        rms_on: str = "raw",
        save_dir: Union[str, Path] = ".",
        save_prefix: str = "montage",
        save_pdf: Optional[Union[str, Path]] = None, # [ADDED]
        max_pdf_pages: Optional[int] = 100, # [ADDED]
        show: bool = True,
        figsize: Optional[Union[str, Tuple[float, float]]] = None, # [ADDED]
        content_aspect: Optional[Union[str, float, Tuple[float, float]]] = None, # [ADDED PAPER]
        paper_margins: Optional[Union[str, Tuple[float, float, float, float], List[float]]] = None, # [ADDED PAPER] （left, right, bottom, top）
    ):
        if axis_type is not None and (xaxis is None or str(xaxis).lower() == "vel"):
            xaxis = axis_type
        
        # Handle multiple inputs and Row Selection
        merged_st = _merge_scantables(scantable)
        self.st = _filter_scantable_by_rows(merged_st, rows=rows, exclude_rows=exclude_rows)

        self.xaxis = str(xaxis).lower()
        if self.xaxis not in ("vel", "freq", "chan"):
            self.xaxis = "vel"

        # 1. Prepare Data & Axes
        print("Preparing data for Montage View...")
        data = self.st.data
        meta = self.st.meta
        
        # --- [MODIFIED] 常に独立処理ルートを使用（Standardizerを完全排除） ---
        self.use_per_row_axis = True
        
        # データをリスト（長さが異なる場合）または2次元配列のまま保持
        if isinstance(data, list):
            self.matrix = data
            self.n_spec = len(data)
            nchan_nominal = len(data[0]) if self.n_spec > 0 else 0
        else:
            self.matrix = np.asarray(data, dtype=float)
            self.n_spec = self.matrix.shape[0]
            nchan_nominal = self.matrix.shape[1] if self.matrix.ndim > 1 else 0

        # 代表値（初期表示用）としてのチャンネル数と軸をセット
        self.nchan = nchan_nominal
        self.ch_axis = np.arange(self.nchan, dtype=float)

        # 全体基準用のNOMINAL軸（個別の描画時には update_plot 内で再計算されます）
        try:
            self.v_axis_native = vlsrk_axis_for_spectrum(meta, v_corr_kms=0.0, nchan=nchan_nominal)
        except Exception:
            self.v_axis_native = np.arange(nchan_nominal, dtype=float)
            
            
        # Cache BEAMEFF/TEMPSCAL
        n_spec = self.n_spec
        self.beameffs = beameff_array(self.st.table, self.st.meta, n_spec)
        self.tempscals = tempscal_array(self.st.table, self.st.meta, n_spec)
        
        first_ts = self.tempscals[0] if len(self.tempscals) > 0 else "TA*"
        self.display_scale = first_ts

        # 2. Setup Frequency Axis
        self.rest_hz_user = float(rest_freq) if rest_freq is not None else None
        val = self.st.meta.get("RESTFREQ") or self.st.meta.get("RESTFRQ")
        self.rest_hz_meta = float(val) if val is not None else None
        
        if self.rest_hz_user is not None and self.rest_hz_user > 0:
            self.rest_hz = self.rest_hz_user
        else:
            self.rest_hz = self.rest_hz_meta
        
        # Velocity Recalculation (Global reference)
        self.v_axis_disp = self.v_axis_native
        if (self.rest_hz_meta is not None and self.rest_hz_meta > 0 and 
            self.rest_hz_user is not None and self.rest_hz_user > 0 and
            abs(self.rest_hz_meta - self.rest_hz_user) > 1.0):
            self.v_axis_disp = recalculate_velocity_axis(self.v_axis_native, self.rest_hz_meta, self.rest_hz_user)
        
        if self.rest_hz is not None and self.rest_hz > 0:
            c_kms = 299792.458
            self.f_axis = self.rest_hz * (1.0 - self.v_axis_disp / c_kms)
            self.f_axis_ghz = self.f_axis * 1e-9
        else:
            self.f_axis_ghz = None

        # 3. Visualization Params
        self.xrange = _norm_range(xrange)
        self.yrange = _norm_range(yrange)
        self.autoscale_y = bool(autoscale_y)
        self.annotate_rms = bool(annotate_rms)
        self.show_fit_windows = bool(show_fit_windows)
        self.smooth_mode = str(smooth_mode).lower()
        self.smooth_width = int(smooth_width)
        self.box_downsample = bool(box_downsample)
        self.box_policy = str(box_policy).lower()
        self.rms_on = str(rms_on).lower()
        self.save_dir = Path(save_dir)
        self.save_prefix = str(save_prefix)
        self.view_mode = "original" 

        imin = builtins.max(0, int(imin))
        if imax is None: imax = self.n_spec
        imax = builtins.min(self.n_spec, int(imax))
        self.indices = list(range(imin, imax))
        self.page = 0
        self.nrows = int(nrows)
        self.ncols = int(ncols)
        self.page_size = builtins.max(1, self.nrows * self.ncols)

        # [MODIFIED] Determine Figure Size
        default_w = self.ncols * 3.0
        default_h = self.nrows * 2.2
        fig_size, bbox_mode = parse_figsize(figsize, default=(default_w, default_h))
        # [ADDED PAPER] Paper mode settings
        self._paper_mode = (bbox_mode is None and fig_size is not None)
        self._paper_size_in = fig_size if self._paper_mode else None
        self._paper_margins_in = parse_paper_margins(paper_margins) if self._paper_mode else None
        self._content_aspect = parse_content_aspect(content_aspect) if self._paper_mode else None
        self._paper_title_y = None
        self._paper_xlabel_y = None
        self._paper_ylabel_x = None
        self._paper_grid_pad = 0.10

        # [MODIFIED] Create Figure
        if show or save_pdf:
            self.fig, self.axes = plt.subplots(
                self.nrows, self.ncols,
                figsize=fig_size,
                sharex=True, sharey=True,
                constrained_layout=(bbox_mode is not None) # Use CL only if not fixed paper
            )
            
                        # [PAPER] Apply paper margins + (A) panel box aspect
            if self._paper_mode and self._paper_size_in is not None and self._paper_margins_in is not None:
                inner = paper_inner_rect_frac(self._paper_size_in, self._paper_margins_in, strict=True)
                axes_arr_paper = np.asarray(self.axes).reshape(self.nrows, self.ncols)
                place_grid_axes_in_rect(axes_arr_paper, inner, self.nrows, self.ncols, pad=self._paper_grid_pad, content_aspect=self._content_aspect)
                pw, ph = float(self._paper_size_in[0]), float(self._paper_size_in[1])
                l_in, r_in, b_in, t_in = self._paper_margins_in
                self._paper_title_y = 1.0 - (float(t_in) / builtins.max(1e-9, ph)) * 0.5
                self._paper_xlabel_y = (float(b_in) / builtins.max(1e-9, ph)) * 0.5
                self._paper_ylabel_x = (float(l_in) / builtins.max(1e-9, pw)) * 0.5
            else:
                # Screen / tight mode
                self.fig.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.15, wspace=0.2, hspace=0.2)

            axes_arr = np.asarray(self.axes).reshape(self.nrows, self.ncols)
            self.axes_flat = axes_arr[::-1, :].ravel()
            self.cid = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
            
            self._print_help()
            self.update_plot()
            
            # [ADDED] PDF Generation
            if save_pdf:
                total_items = len(self.indices)
                n_pages = int(np.ceil(total_items / self.page_size))
                
                def _updater(p):
                    self.page = p
                    self.update_plot()
                
                # [MODIFIED] Enforce tight layout with padding for PDF
                drive_pdf_generation(
                    save_pdf,
                    self.fig,
                    _updater,
                    count=n_pages,
                    max_pages=max_pdf_pages,
                    bbox_inches=(None if self._paper_mode else "tight"),
                    pad_inches=(0.0 if self._paper_mode else 0.75),
                )

                # [FIX] State Reset
                self.page = 0
                self.update_plot()
                
        # [FIX] show=False の場合は自動表示を防ぐためにFigureを閉じる
                if show:
                    # [FIX] Reset to first page before showing
                    self.page = 0
                    self.update_plot()
                else:
                    plt.close(self.fig)
                    self.fig = None
                    
            if show and self.fig is not None:
                plt.show()

    def _print_help(self):
        print("\n--- Montage Controls ---")
        print(" n / p : Next / Previous Page")
        print(" b     : Toggle Baseline View (Residual <-> Raw/Restored)")
        print(" w     : Toggle Fit Windows")
        print(" x     : Cycle x-axis (Vel/Freq/Chan)")
        print(" t     : Toggle Temp Scale (TA* <-> TR*)")
        print(" m     : Cycle smoothing")
        print(" [ / ] : Change smooth width")
        print(" d     : Toggle downsample")
        print(" r / R : Toggle rms_on (raw <-> display)")
        print(" s     : Save Page PDF")
        print(" q     : Quit")
        print("------------------------")

    def _cycle_xaxis(self):
        order = ["vel", "freq", "chan"]
        i = order.index(self.xaxis) if self.xaxis in order else 0
        self.xaxis = order[(i + 1) % len(order)]
        
    def _toggle_scale(self):
        if self.display_scale == "TA*":
            self.display_scale = "TR*"
        else:
            self.display_scale = "TA*"
        print(f"Scale toggled to: {self.display_scale}")

    def _get_current_xaxis_data(self):
        # NOTE: This returns the 'nominal' global axis. 
        if self.xaxis == "vel":
            return self.v_axis_disp, "Velocity (LSRK) [km/s]"
        elif self.xaxis == "freq" and self.f_axis_ghz is not None:
            return self.f_axis_ghz, "Frequency [GHz]"
        else:
            return self.ch_axis, "Channel"

    def update_plot(self):
        total = len(self.indices)
        if total == 0: return
        n_pages = int(np.ceil(total / self.page_size))
        self.page = self.page % builtins.max(n_pages, 1)

        start = self.page * self.page_size
        page_idx = self.indices[start:start + self.page_size]

        for ax in self.axes_flat:
            ax.clear(); ax.axis("off")

        ymins, ymaxs = [], []
        # Get global reference label
        _, xlabel = self._get_current_xaxis_data()

        for k, ax in enumerate(self.axes_flat):
            if k >= len(page_idx): continue
            i = page_idx[k]
            row = self.st.table.iloc[i]
            
            # 1. Get Spectrum Data
            spec_val = self.matrix[i]
            n_points = len(spec_val)

            # 2. Determine Axis for THIS spectrum
            if self.use_per_row_axis:
                # [MODIFIED] Calculate exact axis for this row (like viewer.py)
                v_corr = row.get("v_corr_kms", 0.0)
                if not np.isfinite(v_corr): v_corr = 0.0
                
                # Native LSRK axis (handles TOPO+v_corr)
                vel_native_row = vlsrk_axis_for_spectrum(self.st.meta, v_corr_kms=float(v_corr), nchan=n_points)
                
                # Recalculate if user restfreq is different
                vel_disp_row = vel_native_row
                if (self.rest_hz_meta is not None and self.rest_hz_meta > 0 and 
                    self.rest_hz_user is not None and self.rest_hz_user > 0 and
                    abs(self.rest_hz_meta - self.rest_hz_user) > 1.0):
                    vel_disp_row = recalculate_velocity_axis(vel_native_row, self.rest_hz_meta, self.rest_hz_user)
                
                # Freq axis
                if self.rest_hz is not None and self.rest_hz > 0:
                    c_kms = 299792.458
                    f_axis_ghz_row = (self.rest_hz * (1.0 - vel_disp_row / c_kms)) * 1e-9
                else:
                    f_axis_ghz_row = None
                
                # Select correct axis based on self.xaxis
                if self.xaxis == "vel":
                    x_plot_base = vel_disp_row
                elif self.xaxis == "freq":
                    x_plot_base = f_axis_ghz_row if f_axis_ghz_row is not None else np.arange(n_points)
                else:
                    x_plot_base = np.arange(n_points)
                    
                # Store these for helpers that need them (like baseline/rms)
                current_vel = vel_disp_row
                current_native_vel = vel_native_row
                current_freq = f_axis_ghz_row
                
            else:
                # Use global axis (Standardizer case)
                if self.xaxis == "vel":
                    x_plot_base = self.v_axis_disp
                elif self.xaxis == "freq" and self.f_axis_ghz is not None:
                    x_plot_base = self.f_axis_ghz
                else:
                    x_plot_base = self.ch_axis
                
                current_vel = self.v_axis_disp
                current_native_vel = self.v_axis_native
                current_freq = self.f_axis_ghz
            
            # [MODIFIED] Scale Conversion (Row-wise check)
            row_scale = self.tempscals[i]
            
            if self.display_scale != row_scale:
                eff = self.beameffs[i]
                if self.display_scale == "TR*" and row_scale == "TA*":
                    spec_val = ta_to_tr(spec_val, np.array([eff]))[0]
                elif self.display_scale == "TA*" and row_scale == "TR*":
                    spec_val = tr_to_ta(spec_val, np.array([eff]))[0]

            # Baseline Reconstruction
            arr = _parse_list_like(row.get("BSL_COEF")) if "BSL_COEF" in row.index else None
            has_bsl = arr is not None
            
            baseline_curve = None
            
            if has_bsl:
                try:
                    # Convert to float array
                    coeffs = np.array(arr, dtype=float)
                    # Polyval needs NATIVE velocity (before restfreq shift)
                    baseline_curve = np.polyval(coeffs, current_native_vel)
                    
                    # [MODIFIED] Scale baseline curve too if converting
                    bsl_scale = normalize_tempscal(row.get("BSL_SCALE", row_scale))
                    
                    if self.display_scale != bsl_scale:
                        eff = self.beameffs[i]
                        if self.display_scale == "TR*" and bsl_scale == "TA*":
                            baseline_curve = ta_to_tr(baseline_curve, np.array([eff]))[0]
                        elif self.display_scale == "TA*" and bsl_scale == "TR*":
                            baseline_curve = tr_to_ta(baseline_curve, np.array([eff]))[0]
                            
                except Exception:
                    baseline_curve = None

            is_residual = bool(row.get("BSL_DONE", False))
            
            if self.view_mode == "baseline_added":
                # User wants to see "Raw" data and the baseline curve
                if is_residual and baseline_curve is not None:
                    y_plot = spec_val + baseline_curve
                else:
                    y_plot = spec_val
            else:
                # User wants to see "Current" data (Residual)
                y_plot = spec_val

            # Process / Smooth
            x_disp, y_disp = _process_spectrum(
                x_plot_base, y_plot, 
                self.smooth_mode, self.smooth_width, 
                self.box_downsample, self.box_policy
            )
            
            ax.axis("on")
            ax.step(x_disp, y_disp, where="mid", lw=0.6, color="black")
            
            # Plot Baseline Curve if requested
            if self.view_mode == "baseline_added" and baseline_curve is not None:
                # Need to process baseline curve same as data to match coordinates
                x_bl, y_bl = _process_spectrum(
                    x_plot_base, baseline_curve,
                    self.smooth_mode, self.smooth_width,
                    self.box_downsample, self.box_policy
                )
                ax.plot(x_bl, y_bl, color="red", lw=0.8, alpha=0.8)
            else:
                # Zero line
                ax.axhline(0, lw=0.4, color="gray", ls="--")
                
            # [MODIFIED] Draw Fit Windows
            win_str = str(row.get("BSL_WINF", ""))
            if self.show_fit_windows and win_str and str(win_str).lower() not in ("nan", "none", ""):
                try:
                    wins_x = fit_windows_xaxis(
                        win_str,
                        self.rest_hz_meta,
                        self.rest_hz_user,
                        self.xaxis,
                        current_vel, # Pass correct axis
                        current_freq,
                        self.ch_axis,
                    )
                    for x0, x1 in wins_x:
                        ax.axvspan(float(x0), float(x1), color="green", alpha=0.10, zorder=0)
                except Exception:
                    pass

            ax.set_title(f"#{i}", fontsize=9)

            if self.annotate_rms:
                # [FIXED] RMS calculation with windows logic for Montage
                
                # 1. Determine base data (smoothed vs raw)
                data_rms = y_disp if self.rms_on == "display" else spec_val
                
                # Determine X axis for RMS masking
                x_for_win = x_disp if self.rms_on == "display" else x_plot_base

                # 2. Build Mask
                # (A) X-range mask
                if self.xrange is not None:
                    lo, hi = self.xrange
                    mask_base = (x_for_win >= lo) & (x_for_win <= hi)
                else:
                    mask_base = np.ones_like(data_rms, dtype=bool)

                # (B) Window mask
                mask_win = np.zeros_like(data_rms, dtype=bool)
                has_win = False
                if win_str and str(win_str).lower() not in ("nan", "none", ""):
                    try:
                        wins_x_rms = fit_windows_xaxis(
                            win_str,
                            self.rest_hz_meta,
                            self.rest_hz_user,
                            self.xaxis,
                            current_vel, # Pass correct axis
                            current_freq,
                            self.ch_axis,
                        )
                        if wins_x_rms:
                            has_win = True
                            for x0_w, x1_w in wins_x_rms:
                                mask_win |= (x_for_win >= x0_w) & (x_for_win <= x1_w)
                    except Exception:
                        pass
                
                if has_win:
                    final_mask = mask_base & mask_win
                else:
                    final_mask = mask_base

                rms = _rms_from_mean(data_rms[final_mask])

                if np.isfinite(rms): 
                    ax.text(0.03, 0.88, f"rms={rms:.3g}", transform=ax.transAxes, fontsize=7, va="top")

            if self.xrange is not None: ax.set_xlim(*self.xrange)
            
            mx = np.ones_like(y_disp, dtype=bool)
            if self.xrange is not None:
                lo, hi = self.xrange
                mx = (x_disp >= lo) & (x_disp <= hi)
            
            yy = y_disp[mx]
            yy = yy[np.isfinite(yy)]
            if yy.size: 
                ymins.append(float(np.min(yy)))
                ymaxs.append(float(np.max(yy)))

        if self.yrange is not None:
            for ax in self.axes_flat:
                if ax.has_data(): ax.set_ylim(*self.yrange)
        elif self.autoscale_y and ymins and ymaxs:
            ymin = builtins.min(ymins)
            ymax = builtins.max(ymaxs)
            yr = (ymax - ymin) or 1.0
            lo, hi = ymin - 0.1 * yr, ymax + 0.1 * yr
            for ax in self.axes_flat:
                if ax.has_data(): ax.set_ylim(lo, hi)

        mode_str = "Residual" if self.view_mode == "original" else "Raw + Baseline"
        smooth_tag = f"{self.smooth_mode}(w={int(self.smooth_width)})"
        if self.box_downsample: smooth_tag += "*"
        
                # [PAPER] Figure-level labels within margins
        if getattr(self, "_paper_mode", False) and getattr(self, "_paper_title_y", None) is not None:
            # [変更] グリッドの実体位置（バウンディングボックス）を取得して配置基準にする
            axes_arr = np.asarray(self.axes).reshape(self.nrows, self.ncols)

            # NOTE: タイトル/ラベルは「ページ端」ではなく「図(グリッド)の直上/直下/直左」に置く。
            #       ただし余白が小さい場合に保存時に切れないよう、オフセットを利用可能余白に合わせて縮める。
            top_row_y1, bot_row_y0, left_col_x0, _right_col_x1 = axes_grid_bbox_frac(
                axes_arr,
                fallback=(0.95, 0.05, 0.10, 0.90),
            )
            title_y, label_y, label_x = adjacent_label_positions_from_bbox(
                top_row_y1=top_row_y1,
                bot_row_y0=bot_row_y0,
                left_col_x0=left_col_x0,
                title_pad=0.07,
                xlabel_pad=0.07,
                ylabel_pad=0.08,
                shrink=0.90,
                min_x=0.01,
            )
# fig.suptitle / supxlabel を動的座標で設定
            self.fig.suptitle(f"Montage p.{self.page+1}/{builtins.max(n_pages,1)} | {mode_str} | Scale={self.display_scale} | {smooth_tag}", fontsize=10, y=title_y)
            self.fig.supxlabel(xlabel if xlabel else "X", y=label_y)
            # [変更] x=0.02 固定をやめ、計算した label_x を使用
            self.fig.supylabel(f"Intensity ({self.display_scale}) [{self.st.meta.get('BUNIT','K')}]", x=label_x)
        else:
            self.fig.suptitle(f"Montage p.{self.page+1}/{builtins.max(n_pages,1)} | {mode_str} | Scale={self.display_scale} | {smooth_tag}", fontsize=10)
            self.fig.supxlabel(xlabel if xlabel else "X")
            self.fig.supylabel(f"Intensity ({self.display_scale}) [{self.st.meta.get('BUNIT','K')}] ")
            
        if hasattr(self.fig.canvas, "draw_idle"):
            self.fig.canvas.draw_idle()
        else:
            self.fig.canvas.draw()

    def on_key(self, event):
        if event is None or event.key is None: return
        key = event.key
        total = len(self.indices)
        n_pages = int(np.ceil(total / self.page_size)) if total else 1

        if key == "n": 
            self.page = (self.page + 1) % builtins.max(n_pages, 1)
            self.update_plot()
        elif key == "p": 
            self.page = (self.page - 1 + builtins.max(n_pages, 1)) % builtins.max(n_pages, 1)
            self.update_plot()
        elif key == "b":
            self.view_mode = "baseline_added" if self.view_mode == "original" else "original"
            self.update_plot()
        elif key == "w":
            self.show_fit_windows = not self.show_fit_windows
            self.update_plot()
        elif key == "x": 
            self._cycle_xaxis()
            self.update_plot()
        elif key == "t":
            self._toggle_scale()
            self.update_plot()
        elif key == "m":
            # [MODIFIED] Removed "running" from cycle
            order = ["none", "boxcar"]
            i = order.index(self.smooth_mode) if self.smooth_mode in order else 0
            self.smooth_mode = order[(i + 1) % len(order)]
            self.update_plot()
        elif key == "[": 
            self.smooth_width = builtins.max(1, int(self.smooth_width) - 1)
            self.update_plot()
        elif key == "]": 
            self.smooth_width = int(self.smooth_width) + 1
            self.update_plot()
        elif key == "d": 
            self.box_downsample = not self.box_downsample
            self.update_plot()
        elif key in ("r", "R"): 
            self.rms_on = "display" if self.rms_on == "raw" else "raw"
            self.update_plot()
        elif key == "s":
            save_dir = Path(getattr(self, "save_dir", "."))
            save_dir.mkdir(parents=True, exist_ok=True)
            prefix = getattr(self, "save_prefix", "montage")
            fn = f"{prefix}_page{self.page+1:03d}_{self.view_mode}_{self.display_scale[:2]}.pdf"
            out = save_dir / fn
            self.fig.savefig(out, bbox_inches="tight", pad_inches=0.5)
            print(f"Saved: {out}")
        elif key == "q": 
            plt.close(self.fig)
