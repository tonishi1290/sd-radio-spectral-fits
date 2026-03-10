# src/sd_radio_spectral_fits/plotting/grid.py
from __future__ import annotations

import builtins
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple, List, Union, Sequence

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from ..fitsio import Scantable, read_scantable
from ..regrid_vlsrk import vlsrk_axis_for_spectrum
from ..ranges import parse_windows
# [MODIFIED] Import helpers from scantable_utils
from ..scantable_utils import _parse_row_selector, _df_to_native_endian
from .utils import (
    AxisBundle, _norm_range, _process_spectrum, _rms_from_mean,
    _parse_list_like, _interp_extrap,
    _convert_user_xrange_to_vel_range, _convert_vel_range_to_xrange,
    recalculate_velocity_axis,
    drive_pdf_generation,
    parse_figsize, # [ADDED]
    parse_paper_margins, # [ADDED PAPER]
    parse_content_aspect, # [ADDED PAPER]
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
    
    # Ensure it's a list (safety)
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
    data_all = []
    
    for sc in sc_list:
        tab = _df_to_native_endian(sc.table)
        # Preserve original row IDs per input table when available; otherwise create them
        if "ORIG_ROW_ID" not in tab.columns:
            tab = tab.copy()
            tab["ORIG_ROW_ID"] = np.arange(len(tab), dtype=int)
        tabs.append(tab)
        d = sc.data
        if isinstance(d, list):
            data_all.extend(d)
        else:
            data_all.extend(list(d))
            
    table_all = pd.concat(tabs, axis=0, ignore_index=True)
    
    # Init history for merged
    merged_hist = {"merged_count": len(sc_list), "sources": [str(x) for x in inputs]}
    
    return Scantable(meta=meta, data=data_all, table=table_all, history=merged_hist)


# -------------------------
# Gridding Logic (Helpers)
# -------------------------
def _find_col_ci(df: pd.DataFrame, *candidates: str) -> Optional[str]:
    cols = {str(c).strip().lower(): c for c in df.columns}
    for k in candidates:
        kk = str(k).strip().lower()
        if kk in cols: return cols[kk]
    return None

def _unwrap_deg(a_deg: np.ndarray) -> np.ndarray:
    a = np.asarray(a_deg, dtype=float)
    m = np.isfinite(a)
    if m.sum() < 2: return a
    rad = np.deg2rad(a[m])
    ref = np.nanmedian(rad)
    rad2 = np.unwrap(rad - ref) + ref
    out = a.copy()
    out[m] = np.rad2deg(rad2)
    return out

def _maybe_hours_to_deg(a: np.ndarray, assume_hours_if_le24: bool = True) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    if not assume_hours_if_le24: return a
    m = np.isfinite(a)
    if m.sum() == 0: return a
    mx = np.nanmax(a[m]); mn = np.nanmin(a[m])
    if (mn >= -1.0) and (mx <= 24.5): return a * 15.0
    return a

def _estimate_step(vals: np.ndarray) -> Optional[float]:
    v = np.asarray(vals, dtype=float); v = v[np.isfinite(v)]
    if v.size < 2: return None
    s = np.sort(v); d = np.abs(np.diff(s)); d = d[d > 0]
    if d.size == 0: return None
    lo = np.percentile(d, 5); hi = np.percentile(d, 50)
    d2 = d[(d >= lo) & (d <= hi)]
    step = float(np.median(d2)) if d2.size >= 3 else float(np.median(d))
    return step if (np.isfinite(step) and step > 0) else None

@dataclass(frozen=True)
class GridSolution:
    coord: str
    projection: str
    x_label: str
    y_label: str
    ix: np.ndarray
    iy: np.ndarray
    nx: int
    ny: int
    dx: float
    dy: float
    x0: float
    y0: float
    x_deg: np.ndarray
    y_deg: np.ndarray
    score: float

def _wrap_delta_deg(delta_deg: np.ndarray) -> np.ndarray:
    d = np.asarray(delta_deg, dtype=float)
    return ((d + 180.0) % 360.0) - 180.0

def _fmt_ra(deg: float) -> str:
    val = float(deg)
    if val < 0: val += 360.0
    val = val % 360.0
    hours = val / 15.0
    h = int(hours)
    m = int((hours - h) * 60.0)
    s = (hours - h - m/60.0) * 3600.0
    return f"{h:02d}:{m:02d}:{s:04.1f}"

def _fmt_dec(deg: float) -> str:
    val = float(deg)
    sign = '-' if val < 0 else '+'
    val = abs(val)
    d = int(val)
    m = int((val - d) * 60.0)
    s = (val - d - m/60.0) * 3600.0
    return f"{sign}{d:02d}:{m:02d}:{s:04.1f}"

def _solve_grid_manual_from_lonlat(
    lon_deg: np.ndarray, 
    lat_deg: np.ndarray, 
    coord: str, 
    ref_point: Tuple[float, float], 
    dx_arcsec: float, 
    dy_arcsec: float, 
    corner_offsets_arcsec: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None, 
    grid_anchor_offsets_arcsec: Optional[Tuple[float, float]] = None,
    tol_factor: float = 0.5, 
    invert_x: bool = True, 
    projection: str = "SFL"
) -> Optional[GridSolution]:
    """
    Solve grid mapping manually.
    """
    lon = np.asarray(lon_deg, dtype=float); lat = np.asarray(lat_deg, dtype=float)
    n = int(lon.size)
    if n == 0 or lat.size != n: return None
    m = np.isfinite(lon) & np.isfinite(lat)
    if int(m.sum()) < 1: return None
    
    lon0, lat0 = float(ref_point[0]), float(ref_point[1])
    dx = float(dx_arcsec); dy = float(dy_arcsec)
    if not np.isfinite(dx) or dx <= 0 or not np.isfinite(dy) or dy <= 0: return None
    
    lon_mod = np.mod(lon, 360.0); lon0_mod = float(np.mod(lon0, 360.0))
    dlon = _wrap_delta_deg(lon_mod - lon0_mod)
    dlat = (lat - lat0)
    sgn = -1.0 if invert_x else 1.0
    
    # --- Projection Logic ---
    proj_mode = str(projection).upper()
    if proj_mode in ("TAN", "SIN", "GNOMONIC"):
        # (2) Reference Point Dec: Orthogonal grid (Keeps square pixels)
        cosf = float(np.cos(np.deg2rad(lat0)))
    elif proj_mode in ("GLS", "SFL", "SINE"):
        # (1) Local Dec: Equal Area (GLS), standard for large single-dish maps
        cosf = np.cos(np.deg2rad(lat))
        print(lat)
    else:
        # NONE or raw
        cosf = 1.0
        
    x_arcsec = sgn * dlon * cosf * 3600.0
    y_arcsec = dlat * 3600.0

    use_anchor = (grid_anchor_offsets_arcsec is not None)
    if use_anchor:
        xa, ya = grid_anchor_offsets_arcsec
        xa = float(xa)
        ya = float(ya)
    else:
        xa = None
        ya = None

    if corner_offsets_arcsec is not None:
        (x1, y1), (x2, y2) = corner_offsets_arcsec
        xmin_req, xmax_req = sorted([float(x1), float(x2)])
        ymin_req, ymax_req = sorted([float(y1), float(y2)])
    else:
        xmin_req = float(np.nanmin(x_arcsec[m])); xmax_req = float(np.nanmax(x_arcsec[m]))
        ymin_req = float(np.nanmin(y_arcsec[m])); ymax_req = float(np.nanmax(y_arcsec[m]))

    if use_anchor:
        if corner_offsets_arcsec is not None:
            ix_min = int(np.ceil((xmin_req - xa) / dx))
            ix_max = int(np.floor((xmax_req - xa) / dx))
            iy_min = int(np.ceil((ymin_req - ya) / dy))
            iy_max = int(np.floor((ymax_req - ya) / dy))
        else:
            # Auto range: expand to anchor-aligned grid that covers all data
            ix_min = int(np.floor((xmin_req - xa) / dx))
            ix_max = int(np.ceil((xmax_req - xa) / dx))
            iy_min = int(np.floor((ymin_req - ya) / dy))
            iy_max = int(np.ceil((ymax_req - ya) / dy))

        nx = int(ix_max - ix_min + 1)
        ny = int(iy_max - iy_min + 1)
        if nx <= 0 or ny <= 0: return None

        xmin = float(xa + ix_min * dx)
        ymin = float(ya + iy_min * dy)

        ix_abs = np.rint((x_arcsec - xa) / dx).astype(int)
        iy_abs = np.rint((y_arcsec - ya) / dy).astype(int)
        ix0 = ix_abs - ix_min
        iy0 = iy_abs - iy_min

        xg = xa + ix_abs * dx
        yg = ya + iy_abs * dy
    else:
        xmin = float(xmin_req)
        ymin = float(ymin_req)
        if corner_offsets_arcsec is not None:
            xmax = float(xmax_req)
            ymax = float(ymax_req)
            nx = int(np.ceil((xmax - xmin) / dx)) + 1
            ny = int(np.ceil((ymax - ymin) / dy)) + 1
        else:
            xmax = float(xmax_req)
            ymax = float(ymax_req)
            nx = int(np.rint((xmax - xmin) / dx)) + 1
            ny = int(np.rint((ymax - ymin) / dy)) + 1

        if nx <= 0 or ny <= 0: return None

        ix0 = np.rint((x_arcsec - xmin) / dx).astype(int)
        iy0 = np.rint((y_arcsec - ymin) / dy).astype(int)

        xg = xmin + ix0 * dx
        yg = ymin + iy0 * dy
    tolx = abs(float(tol_factor)) * dx; toly = abs(float(tol_factor)) * dy
    
    good = m & (ix0 >= 0) & (ix0 < nx) & (iy0 >= 0) & (iy0 < ny) & (np.abs(x_arcsec - xg) <= tolx) & (np.abs(y_arcsec - yg) <= toly)
    
    ix = np.full(n, -1, dtype=int); iy_disp = np.full(n, -1, dtype=int)
    ix[good] = ix0[good]; iy_disp[good] = (ny - 1) - iy0[good]
    
    if int(good.sum()) > 0:
        cells = (iy_disp[good].astype(np.int64) * np.int64(nx)) + ix[good].astype(np.int64)
        unique = int(np.unique(cells).size); total = int(good.sum())
        dup = int(total - unique); occ = float(unique) / float(nx * ny)
        dup_rate = float(dup) / float(total) if total > 0 else 0.0
        score = float(occ - 0.35 * dup_rate)
    else: score = -1.0
    
    x_label, y_label = ("RA", "DEC") if str(coord).lower() in ("radec", "ra-dec", "ra_dec") else ("GLON", "GLAT")
    return GridSolution(
        coord=str(coord).lower(), projection=proj_mode,
        x_label=x_label, y_label=y_label, 
        ix=ix, iy=iy_disp, nx=nx, ny=ny, dx=float(dx), dy=float(dy), 
        x0=float(lon0), y0=float(lat0), x_deg=_unwrap_deg(lon), y_deg=lat, score=float(score)
    )


class ProfileMapGridViewer:
    """
    Profile Map Viewer using manual gridding and Standardizer.
    
    Parameters
    ----------
    st : Union[Scantable, str, Sequence[Union[Scantable, str]]]
        Input data.
    rows : Union[int, slice, List[int], str], optional
        Select specific rows (global indices) to display.
        Examples: 0, "0,2,5", "10:20". Cannot be used with exclude_rows.
    exclude_rows : Union[int, slice, List[int], str], optional
        Exclude specific rows (global indices) from display.
        Cannot be used with rows.
    mode : str
        "coadd" to average spectra in a cell.
        "all" (default) to show all spectra in a cell individually.
    """
    def __init__(
        self,
        st: Union[Scantable, str, Sequence[Union[Scantable, str]]],
        rows: Optional[Union[int, slice, List[int], str]] = None,
        exclude_rows: Optional[Union[int, slice, List[int], str]] = None,
        mode: str = "all",
        combine: str = "mean",
        coord: str = "radec",
        projection: str = "TAN", 
        ref_point: Optional[Tuple[float, float]] = None,
        x_grid: float = 30.0,
        y_grid: float = 30.0,
        corner_offsets: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
        grid_bounds_offsets: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
        grid_anchor_offsets: Optional[Tuple[float, float]] = None,
        grid_tol: float = 0.5,
        invert_x: bool = True,
        show_rms: bool = True,
        show_coord_text: bool = True,
        nrows: Optional[int] = None,
        ncols: Optional[int] = None,
        xaxis: str = "vel",
        axis_type: Optional[str] = None,
        rest_freq: Optional[float] = None,
        imin: int = 0,
        imax: Optional[int] = None,
        xrange: Optional[Tuple[float, float]] = None,
        yrange: Optional[Tuple[float, float]] = None,
        annotate_rms: bool = True,
        smooth_mode: str = "none",
        smooth_width: int = 1,
        box_downsample: bool = False,
        box_policy: str = "trim",
        rms_on: str = "raw",
        show: bool = True,
        square_aspect: bool = True,
        box_padding: float = 0.0,
        show_baseline: bool = True,
        show_id: bool = True, # [ADDED] ID表示の切り替え
        show_fit_windows: bool = True,
        show_grid: Union[bool, str] = "auto",
        offset_unit: str = "arcsec",
        save_dir: Union[str, Path] = ".",
        save_prefix: str = "grid",
        save_pdf: Optional[Union[str, Path]] = None, # [ADDED]
        max_pdf_pages: Optional[int] = 100, # [ADDED]
        figsize: Optional[Union[str, Tuple[float, float]]] = None, # [ADDED]
        content_aspect: Optional[Union[str, float, Tuple[float, float]]] = None, # [ADDED PAPER]
        paper_margins: Optional[Union[str, Tuple[float, float, float, float], List[float]]] = None, # [ADDED PAPER]
    ):
        if axis_type is not None and (xaxis is None or str(xaxis).lower() == "vel"):
            xaxis = axis_type
        
        # Handle multiple inputs and Row Selection
        merged_st = _merge_scantables(st)
        self.st = _filter_scantable_by_rows(merged_st, rows=rows, exclude_rows=exclude_rows)
        
        self.mode = str(mode).lower() 

        self.xaxis = str(xaxis).lower()
        if self.xaxis not in ("vel", "freq", "chan"): self.xaxis = "vel"

        # [MODIFIED] Data Preparation (Standardizer完全排除)
        print("Preparing data for Grid View...")
        data = self.st.data
        
        # リストで渡された場合も、チャンネル数が同じなら綺麗な2次元配列として扱う
        if isinstance(data, list):
            try:
                self.matrix = np.array(data, dtype=float)
            except ValueError:
                self.matrix = data # 万が一長さが違った場合はリストのまま保持
        else:
            self.matrix = np.asarray(data, dtype=float)

        n_spec = len(self.matrix)
        if n_spec > 0:
            nchan_nominal = len(self.matrix[0]) if isinstance(self.matrix, list) else self.matrix.shape[1]
        else:
            nchan_nominal = 0
            
        self.nchan = nchan_nominal

        # 代表値としての軸をセット
        try:
            self.v_axis_native = vlsrk_axis_for_spectrum(self.st.meta, v_corr_kms=0.0, nchan=nchan_nominal)
        except Exception:
            self.v_axis_native = np.arange(nchan_nominal, dtype=float)

        # Cache BEAMEFF and TEMPSCAL
        n_rows = len(self.st.table)
        self.beameffs = beameff_array(self.st.table, self.st.meta, n_rows)
        self.tempscals = tempscal_array(self.st.table, self.st.meta, n_rows)
        
        # Default display scale
        first_ts = self.tempscals[0] if len(self.tempscals) > 0 else "TA*"
        self.display_scale = first_ts 

        # Setup Freq Axis
        self.rest_hz_user = float(rest_freq) if rest_freq is not None else None
        
        val = self.st.meta.get("RESTFREQ") or self.st.meta.get("RESTFRQ")
        self.rest_hz_meta = float(val) if val is not None else None
        
        if self.rest_hz_user is not None and self.rest_hz_user > 0:
            self.rest_hz = self.rest_hz_user
        else:
            self.rest_hz = self.rest_hz_meta

        # Velocity Recalculation
        self.v_axis_disp = self.v_axis_native
        if (self.rest_hz_meta is not None and self.rest_hz_meta > 0 and 
            self.rest_hz_user is not None and self.rest_hz_user > 0 and
            abs(self.rest_hz_meta - self.rest_hz_user) > 1.0):
            self.v_axis_disp = recalculate_velocity_axis(self.v_axis_native, self.rest_hz_meta, self.rest_hz_user)
        
        if self.rest_hz is not None and self.rest_hz > 0:
            c_kms = 299792.458
            self.f_axis_ghz = (self.rest_hz * (1.0 - self.v_axis_disp / c_kms)) * 1e-9
        else:
            self.f_axis_ghz = None
            
        self.ch_axis = np.arange(self.nchan)

        n_spec = len(self.st.data)
        imin = builtins.max(0, int(imin))
        if imax is None: imax = n_spec
        imax = builtins.min(n_spec, int(imax))
        self.indices = np.arange(imin, imax, dtype=int)
        self._nrows_req = None if nrows is None else int(nrows)
        self._ncols_req = None if ncols is None else int(ncols)
        if self._nrows_req is not None and self._nrows_req <= 0:
            raise ValueError("nrows must be >= 1 when specified.")
        if self._ncols_req is not None and self._ncols_req <= 0:
            raise ValueError("ncols must be >= 1 when specified.")

        # 1. Coordinate Gridding
        coord_l = str(coord).lower()
        if coord_l == "auto": raise ValueError("coord='auto' is disabled. Please specify coord='radec' or coord='gal'.")

        if coord_l in ("radec", "ra-dec", "ra_dec"):
            lon_col = _find_col_ci(self.st.table, "RA", "RA_DEG", "RA2000", "RAJ2000")
            lat_col = _find_col_ci(self.st.table, "DEC", "DEC_DEG", "DEC2000", "DECJ2000")
            if lon_col is None or lat_col is None: raise ValueError("RA/DEC columns not found for coord='radec'.")
            lon_raw = self.st.table[lon_col].to_numpy(float)
            lat_deg = self.st.table[lat_col].to_numpy(float)
            lon_deg = _maybe_hours_to_deg(lon_raw)
            if ref_point is None:
                lon0 = float(np.nanmedian(lon_deg[np.isfinite(lon_deg)]))
                lat0 = float(np.nanmedian(lat_deg[np.isfinite(lat_deg)]))
            else:
                lon0 = float(ref_point[0]); lat0 = float(ref_point[1])
                if (np.nanmax(np.abs(lon_raw[np.isfinite(lon_raw)])) <= 24.5 and abs(lon0) <= 24.5): lon0 *= 15.0
            coord_use = "radec"
        elif coord_l in ("gal", "galactic", "glon-glat", "glon_glat"):
            lon_col = _find_col_ci(self.st.table, "GLON", "L", "GAL_L", "LON")
            lat_col = _find_col_ci(self.st.table, "GLAT", "B", "GAL_B", "LAT")
            if lon_col is None or lat_col is None: raise ValueError("GLON/GLAT columns not found for coord='gal'.")
            lon_deg = np.asarray(self.st.table[lon_col].to_numpy(float), dtype=float)
            lat_deg = np.asarray(self.st.table[lat_col].to_numpy(float), dtype=float)
            if ref_point is None:
                lon0 = float(np.nanmedian(lon_deg[np.isfinite(lon_deg)]))
                lat0 = float(np.nanmedian(lat_deg[np.isfinite(lat_deg)]))
            else:
                lon0 = float(ref_point[0]); lat0 = float(ref_point[1])
            coord_use = "gal"
        else: raise ValueError("coord must be 'radec' or 'gal'.")

        # Range parameter alias handling (corner_offsets kept for backward compatibility)
        if corner_offsets is not None and grid_bounds_offsets is not None:
            if tuple(corner_offsets) != tuple(grid_bounds_offsets):
                raise ValueError("Specify only one of 'corner_offsets' or 'grid_bounds_offsets' (or pass identical values).")
            grid_bounds_use = corner_offsets
        elif grid_bounds_offsets is not None:
            grid_bounds_use = grid_bounds_offsets
        else:
            grid_bounds_use = corner_offsets

        self.sol = _solve_grid_manual_from_lonlat(
            lon_deg=lon_deg, lat_deg=lat_deg, coord=coord_use, ref_point=(lon0, lat0),
            dx_arcsec=float(x_grid), dy_arcsec=float(y_grid), corner_offsets_arcsec=grid_bounds_use,
            grid_anchor_offsets_arcsec=grid_anchor_offsets,
            tol_factor=float(grid_tol), invert_x=bool(invert_x), projection=str(projection),
        )
        if self.sol is None: raise ValueError("Failed to build grid solution. Check coord/ref_point/grid settings.")

        # Display tiling: if nrows/ncols are omitted, show the full grid on one page
        self.nrows = int(self.sol.ny if self._nrows_req is None else self._nrows_req)
        self.ncols = int(self.sol.nx if self._ncols_req is None else self._ncols_req)

        self.grid = [[[] for _ in range(self.sol.nx)] for _ in range(self.sol.ny)]
        for idx in self.indices:
            if idx < 0 or idx >= n_spec: continue
            ix = int(self.sol.ix[idx]); iy = int(self.sol.iy[idx])
            if 0 <= ix < self.sol.nx and 0 <= iy < self.sol.ny:
                self.grid[iy][ix].append(int(idx))

        # 2. Visualization Setup
        self._xaxis_init = str(self.xaxis)
        self.xrange_user = _norm_range(xrange)
        self.xrange = self.xrange_user
        self.xr_vel = None
        if self.xrange_user is not None:
            try:
                self.xr_vel = _convert_user_xrange_to_vel_range(
                    self._xaxis_init,
                    self.xrange_user,
                    self.v_axis_disp,
                    self.ch_axis,
                    self.f_axis_ghz,
                    self.rest_hz,
                    0.0,
                )
            except Exception:
                self.xr_vel = None
        self.yrange = _norm_range(yrange)
        self.annotate_rms = bool(annotate_rms)
        self.smooth_mode = str(smooth_mode).lower()
        self.smooth_width = int(smooth_width)
        self.box_downsample = bool(box_downsample)
        self.box_policy = str(box_policy).lower()
        self.rms_on = str(rms_on).lower()
        self.combine = str(combine).lower()
        self.show_coord_text = bool(show_coord_text)
        
        self.square_aspect = bool(square_aspect)
        self.box_padding = float(box_padding)
        self.show_baseline = bool(show_baseline)
        self.show_fit_windows = bool(show_fit_windows)
        self.show_grid = show_grid
        self.show_id = bool(show_id) # [ADDED]
        self.offset_unit = str(offset_unit).strip().lower()
        self.save_dir = Path(save_dir)
        self.save_prefix = str(save_prefix)
        self.grid_anchor_offsets = None if grid_anchor_offsets is None else (float(grid_anchor_offsets[0]), float(grid_anchor_offsets[1]))
        self.grid_bounds_offsets = None if grid_bounds_use is None else (
            (float(grid_bounds_use[0][0]), float(grid_bounds_use[0][1])),
            (float(grid_bounds_use[1][0]), float(grid_bounds_use[1][1])),
        )
        self.n_win_x = int(np.ceil(self.sol.nx / self.ncols))
        self.n_win_y = int(np.ceil(self.sol.ny / self.nrows))

        self._page_map: List[Tuple[int, int]] = []
        for _py in range(self.n_win_y):
            for _px in range(self.n_win_x):
                _x0 = int(_px * self.ncols)
                _y0 = int(_py * self.nrows)
                _has = False
                for _gy in range(_y0, builtins.min(_y0 + self.nrows, self.sol.ny)):
                    for _gx in range(_x0, builtins.min(_x0 + self.ncols, self.sol.nx)):
                        if self.grid[_gy][_gx]:
                            _has = True
                            break
                    if _has:
                        break
                if _has:
                    self._page_map.append((int(_px), int(_py)))

        if not self._page_map:
            self._page_map = [(0, 0)]

        self.n_pages = int(len(self._page_map))
        self.page = 0

        # --- Layout & Figure Setup ---
        
        # [MODIFIED] Determine Figure Size and BBox Policy
        # Determine Margins (Fixed safety margins)
        self._marg_left_in = 0.8
        self._marg_right_in = 0.2
        self._marg_bottom_in = 0.8
        self._marg_top_in = 0.75
        pad = builtins.max(0.0, float(self.box_padding))
        
        # Parse user preference
        user_size, bbox_mode = parse_figsize(figsize, default=None)
        # [ADDED PAPER] Paper mode settings
        self._paper_mode = (bbox_mode is None and user_size is not None)
        self._paper_size_in = user_size if self._paper_mode else None
        self._paper_margins_in = parse_paper_margins(paper_margins) if self._paper_mode else None
        self._content_aspect = parse_content_aspect(content_aspect) if self._paper_mode else None
        if self._paper_mode and self._paper_margins_in is not None:
            ml, mr, mb, mt = self._paper_margins_in
            self._marg_left_in = float(ml)
            self._marg_right_in = float(mr)
            self._marg_bottom_in = float(mb)
            self._marg_top_in = float(mt)

        if user_size is not None:
            # --- [Fit-to-Paper Mode] ---
            # Fixed Figure Size -> Calculate Panel Size
            fig_w_in, fig_h_in = user_size
            
            # Effective area for panels
            eff_w = fig_w_in - (self._marg_left_in + self._marg_right_in)
            eff_h = fig_h_in - (self._marg_bottom_in + self._marg_top_in)

            # [ADDED] Fail-fast if margins consume the whole page
            if eff_w <= 0 or eff_h <= 0:
                raise ValueError(
                    "Invalid paper_margins: margins exceed paper size. "
                    f"paper_size_in=({fig_w_in:.6g},{fig_h_in:.6g}) "
                    f"margins_in(L,R,B,T)=({self._marg_left_in:.6g},{self._marg_right_in:.6g},{self._marg_bottom_in:.6g},{self._marg_top_in:.6g})"
                )
            
            # Calculate required space for gaps (ncols-1 gaps)
            # Gap is defined relative to panel size: gap = panel * pad
            # Width = ncols * panel + (ncols-1) * (panel * pad)
            #       = panel * (ncols + (ncols-1)*pad)
            denom_w = float(self.ncols) + float(builtins.max(self.ncols - 1, 0)) * pad
            denom_h = float(self.nrows) + float(builtins.max(self.nrows - 1, 0)) * pad
            
            # Calculate max possible panel size to fit width and height
            panel_w_raw = eff_w / denom_w
            panel_h_raw = eff_h / denom_h
            
            # Determine actual panel size (Square or Rect)
            if self.square_aspect:
                panel_size = builtins.min(panel_w_raw, panel_h_raw)
                panel_w_in = panel_size
                panel_h_in = panel_size
            else:
                panel_w_in = panel_w_raw
                panel_h_in = panel_h_raw
                
            # Note: We don't change fig_w_in/fig_h_in here, we use the paper size.
            # Layout logic in _apply_square_layout handles centering.
            
        else:
            # --- [Stacking Mode (Legacy)] ---
            # Fixed Panel Size -> Calculate Figure Size
            panel_w_in = 3.0
            panel_h_in = 3.0 if self.square_aspect else 2.0
            
            bbox_mode = "tight" # Default behavior if no size specified

            gap_w_in = panel_w_in * pad
            gap_h_in = panel_h_in * pad

            content_w_in = (self.ncols * panel_w_in) + (builtins.max(self.ncols - 1, 0) * gap_w_in)
            content_h_in = (self.nrows * panel_h_in) + (builtins.max(self.nrows - 1, 0) * gap_h_in)

            fig_w_in = content_w_in + self._marg_left_in + self._marg_right_in
            fig_h_in = content_h_in + self._marg_bottom_in + self._marg_top_in

        # [MODIFIED] Create Figure
        if show or save_pdf:
            self.fig = plt.figure(figsize=(fig_w_in, fig_h_in))

            axs = []
            n_axes = int(self.nrows) * int(self.ncols)
            for i in range(n_axes):
                if i == 0:
                    ax = self.fig.add_axes([0.1, 0.1, 0.1, 0.1])
                else:
                    ax = self.fig.add_axes([0.1, 0.1, 0.1, 0.1], sharex=axs[0], sharey=axs[0])
                axs.append(ax)

            self.axes = np.array(axs, dtype=object).reshape(self.nrows, self.ncols)
            self.axes_flat = np.atleast_1d(self.axes).flatten()

            self._apply_square_layout() # This uses panel_w_in calculated above
            self._rid = self.fig.canvas.mpl_connect("resize_event", self._on_resize)
            self.cid = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
            self._print_help()
            
            self.update_plot()
            
            # [ADDED] PDF Generation
            if save_pdf:
                def _updater(p):
                    self.page = p
                    self.update_plot()
                
                # [MODIFIED] Enforce tight layout with padding for PDF
                drive_pdf_generation(
                    save_pdf,
                    self.fig,
                    _updater,
                    count=self.n_pages,
                    max_pages=max_pdf_pages,
                    bbox_inches=(None if self._paper_mode else "tight"),
                    pad_inches=(0.0 if self._paper_mode else 0.75)
                )
                
                # [FIX] State reset and cleanup
                self.page = 0
                self.update_plot()

                if not show:
                    plt.close(self.fig)
                    self.fig = None
            if show and self.fig is not None:
                plt.show()


    # -------------------------------------------------------------------------
    # Layout helpers
    # -------------------------------------------------------------------------
    def _apply_square_layout(self) -> None:
        if not hasattr(self, "fig") or not hasattr(self, "axes"): return

        fig_w_in, fig_h_in = self.fig.get_size_inches()
        fig_w_in = float(fig_w_in); fig_h_in = float(fig_h_in)
        if not np.isfinite(fig_w_in) or not np.isfinite(fig_h_in) or fig_w_in <= 0 or fig_h_in <= 0:
            return

        ml = float(getattr(self, "_marg_left_in", 1.2))
        mr = float(getattr(self, "_marg_right_in", 0.2))
        mb = float(getattr(self, "_marg_bottom_in", 1.0))
        mt = float(getattr(self, "_marg_top_in", 0.8))

        aw = fig_w_in - (ml + mr)
        ah = fig_h_in - (mb + mt)
        aw = builtins.max(1e-6, aw)
        ah = builtins.max(1e-6, ah)

        pad = builtins.max(0.0, float(getattr(self, "box_padding", 0.0)))

        denom_w = float(self.ncols) + float(builtins.max(self.ncols - 1, 0)) * pad
        denom_h = float(self.nrows) + float(builtins.max(self.nrows - 1, 0)) * pad
        denom_w = builtins.max(1e-9, denom_w)
        denom_h = builtins.max(1e-9, denom_h)

        if bool(getattr(self, "square_aspect", True)):
            cell = builtins.min(aw / denom_w, ah / denom_h)
            cell = builtins.max(1e-6, float(cell))
            gap = cell * pad

            grid_w = (float(self.ncols) * cell) + float(builtins.max(self.ncols - 1, 0)) * gap
            grid_h = (float(self.nrows) * cell) + float(builtins.max(self.nrows - 1, 0)) * gap

            x0_in = ml + (aw - grid_w) / 2.0
            y0_in = mb + (ah - grid_h) / 2.0

            for r in range(self.nrows):
                for c in range(self.ncols):
                    left_in = x0_in + float(c) * (cell + gap)
                    bottom_in = y0_in + float((self.nrows - 1) - r) * (cell + gap)
                    self.axes[r, c].set_position([
                        left_in / fig_w_in,
                        bottom_in / fig_h_in,
                        cell / fig_w_in,
                        cell / fig_h_in,
                    ])
        else:
            # Non-square panels. If content_aspect (width/height) is specified, enforce it.
            aspect = getattr(self, "_content_aspect", None)
            if aspect is not None and np.isfinite(float(aspect)) and float(aspect) > 0:
                a = builtins.max(1e-9, float(aspect))
                cell_h = builtins.min(aw / (a * denom_w), ah / denom_h)
                cell_h = builtins.max(1e-6, float(cell_h))
                cell_w = builtins.max(1e-6, float(cell_h * a))
            else:
                cell_w = builtins.max(1e-6, float(aw / denom_w))
                cell_h = builtins.max(1e-6, float(ah / denom_h))
            gap_w = cell_w * pad
            gap_h = cell_h * pad

            grid_w = (float(self.ncols) * cell_w) + float(builtins.max(self.ncols - 1, 0)) * gap_w
            grid_h = (float(self.nrows) * cell_h) + float(builtins.max(self.nrows - 1, 0)) * gap_h

            x0_in = ml + (aw - grid_w) / 2.0
            y0_in = mb + (ah - grid_h) / 2.0

            for r in range(self.nrows):
                for c in range(self.ncols):
                    left_in = x0_in + float(c) * (cell_w + gap_w)
                    bottom_in = y0_in + float((self.nrows - 1) - r) * (cell_h + gap_h)
                    self.axes[r, c].set_position([
                        left_in / fig_w_in,
                        bottom_in / fig_h_in,
                        cell_w / fig_w_in,
                        cell_h / fig_h_in,
                    ])

    def _on_resize(self, event) -> None:
        try:
            self._apply_square_layout()
            if hasattr(self, "fig") and hasattr(self.fig, "canvas"):
                self.fig.canvas.draw_idle()
        except Exception:
            return

    def _print_help(self):
        print("\n--- Controls (Grid) ---")
        print(" n / p : Next / Previous Page")
        print(" x     : Cycle x-axis (vel -> freq -> chan)")
        print(" m     : Cycle smoothing")
        print(" [ / ] : Change smooth width")
        print(" d     : Toggle downsample")
        print(" t     : Toggle Temp Scale (TA* <-> TR*)")
        print(" r / R : Toggle rms_on (raw <-> display)")
        print(" b     : Toggle fit windows (BSL_WINF)")
        print(" w     : Toggle fit windows (alias of b)")
        print(" g     : Toggle grid lines (auto -> off -> on)")
        print(" a     : Toggle square aspect (and re-layout)")
        print(" s     : Save PDF")
        print(" h     : Print this help")
        print(" q     : Quit")
        print("-----------------------")

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
        if self.xaxis == "vel":
            return self.v_axis_disp, "Velocity (LSRK) [km/s]"
        elif self.xaxis == "freq" and self.f_axis_ghz is not None:
            return self.f_axis_ghz, "Frequency [GHz]"
        else:
            return self.ch_axis, "Channel"

    def _xlim_for_current_axis(self, x: Optional[np.ndarray] = None) -> Optional[Tuple[float, float]]:
        if getattr(self, "xr_vel", None) is not None:
            try:
                return _convert_vel_range_to_xrange(
                    self.xaxis,
                    self.xr_vel,
                    self.v_axis_disp,
                    self.ch_axis,
                    self.f_axis_ghz,
                    self.rest_hz,
                )
            except Exception:
                pass

        if getattr(self, "xrange_user", None) is not None and str(getattr(self, "_xaxis_init", self.xaxis)) == str(self.xaxis):
            return self.xrange_user

        return None

    def _cell_spectrum(self, idx_list: List[int]) -> Optional[Tuple[pd.Series, np.ndarray]]:
        if not idx_list: return None
        
        target = self.display_scale
        
        # Get data and scales for all indices in this cell
        # Optimization: if len=1, simpler path.
        if len(idx_list) == 1:
            i0 = idx_list[0]
            y = self.matrix[i0]
            row = self.st.table.iloc[i0]
            scale = self.tempscals[i0]
            
            if scale != target:
                eff = self.beameffs[i0]
                if target == "TR*" and scale == "TA*":
                    y = ta_to_tr(y, np.array([eff]))[0]
                elif target == "TA*" and scale == "TR*":
                    y = tr_to_ta(y, np.array([eff]))[0]
            return row, y

        # Average case
        if isinstance(self.matrix, list):
            data_sub = np.array([self.matrix[i] for i in idx_list], dtype=float)
        else:
            data_sub = self.matrix[idx_list].copy() # (N, Ch)
            
        scales_sub = self.tempscals[idx_list] # (N,)
        effs_sub = self.beameffs[idx_list] # (N,)
        
        # We want to convert everything to 'target' scale before averaging
        
        # Mask of rows that need TA->TR
        if target == "TR*":
            # Convert TA* rows to TR*
            mask_ta = (scales_sub == "TA*")
            if np.any(mask_ta):
                # [FIXED] Add [:, None] to broadcast correctly (N,) -> (N,1)
                data_sub[mask_ta] = ta_to_tr(data_sub[mask_ta], effs_sub[mask_ta, None])
                
        elif target == "TA*":
            # Convert TR* rows to TA*
            mask_tr = (scales_sub == "TR*")
            if np.any(mask_tr):
                # [FIXED] Add [:, None] to broadcast correctly
                data_sub[mask_tr] = tr_to_ta(data_sub[mask_tr], effs_sub[mask_tr, None])
        
        y = np.nanmean(data_sub, axis=0)
        row = self.st.table.iloc[idx_list[0]]
        return row, y

    def update_plot(self):
        # [MODIFIED] Construct info string for the figure title (like viewer.py)
        smooth_tag = f"{self.smooth_mode}(w={int(self.smooth_width)})"
        if self.smooth_mode in ("boxcar", "box", "avg", "average") and self.box_downsample: 
            smooth_tag += ",downsample"
        
        title_str = f"Scale: {self.display_scale} | Smooth: {smooth_tag} | RMS Mode: {self.rms_on}"
        
        if getattr(self, "_paper_mode", False) and getattr(self, "_paper_size_in", None) is not None:
            # [変更] グリッドの実体位置（バウンディングボックス）を取得して配置基準にする
            axes_arr = np.asarray(self.axes).reshape(self.nrows, self.ncols)

            # NOTE: タイトルは「ページ端」ではなく「図(グリッド)の直上」に置く。
            #       ただし余白が小さい場合に保存時に切れないよう、オフセットを利用可能余白に合わせて縮める。
            top_row_y1, bot_row_y0, left_col_x0, _right_col_x1 = axes_grid_bbox_frac(
                axes_arr,
                fallback=(0.95, 0.05, 0.10, 0.90),
            )
            title_y, _xlabel_y, _ylabel_x = adjacent_label_positions_from_bbox(
                top_row_y1=top_row_y1,
                bot_row_y0=bot_row_y0,
                left_col_x0=left_col_x0,
                title_pad=0.07,
                xlabel_pad=0.07,
                ylabel_pad=0.08,
                shrink=0.90,
                min_x=0.01,
            )
            self.fig.suptitle(title_str, fontsize=10, y=title_y)

        else:
            self.fig.suptitle(title_str, fontsize=10, y=0.98)


        if hasattr(self, "_page_map") and getattr(self, "_page_map", None):
            px, py = self._page_map[int(self.page) % builtins.max(1, int(len(self._page_map)))]
        else:
            px = int(self.page % self.n_win_x)
            py = int(self.page // self.n_win_x)

        x0 = int(px) * self.ncols
        y0 = int(py) * self.nrows

        for ax in self.axes_flat:
            ax.clear(); ax.axis("off")

        ymins, ymaxs = [], []
        x, xlabel_text = self._get_current_xaxis_data()
        ylabel_text = f"Intensity ({self.display_scale}) [{self.st.meta.get('BUNIT','K')}]"
        xlim = self._xlim_for_current_axis(x=x)

        ref_x, ref_y = self.sol.x0, self.sol.y0
        cos_dec_ref = 1.0
        if self.sol.coord == "radec":
            cos_dec_ref = np.cos(np.deg2rad(ref_y))
        
        do_grid = False
        if self.show_grid == "auto":
            do_grid = True
        elif isinstance(self.show_grid, bool):
            do_grid = self.show_grid
        elif str(self.show_grid).lower() not in ("none", "false", "0"):
            do_grid = True

        unit_scale = 1.0
        unit_label = '"'
        if self.offset_unit.startswith("arcmin"):
            unit_scale = 1.0 / 60.0
            unit_label = "'"
        elif self.offset_unit.startswith("deg"):
            unit_scale = 1.0 / 3600.0
            unit_label = "deg"

        for r in range(self.nrows):
            for c in range(self.ncols):
                ax = self.axes[r, c] if hasattr(self.axes, "__getitem__") else self.axes_flat[r * self.ncols + c]
                gx = x0 + c
                gy = y0 + r
                
                if self.square_aspect:
                    ax.set_box_aspect(1)
                else:
                    ax.set_box_aspect(None)

                # Determine if this is the "Main" bottom-left anchor axis
                is_anchor_ax = (r == self.nrows - 1) and (c == 0)

                idx_list = []
                if gx < self.sol.nx and gy < self.sol.ny:
                    idx_list = self.grid[gy][gx]

                # [MODIFIED] Force axis ON if it is the anchor, even if empty
                if not idx_list and not is_anchor_ax:
                    continue
                
                ax.axis("on")
                
                # [MODIFIED] Specific Labels/Ticks for bottom-left axis ONLY
                if is_anchor_ax:
                    ax.tick_params(labelbottom=True, labelleft=True)
                    ax.set_xlabel(xlabel_text if xlabel_text else "X", fontsize="small")
                    ax.set_ylabel(ylabel_text, fontsize="small")
                else:
                    ax.tick_params(labelbottom=False, labelleft=False)

                if do_grid:
                    ax.grid(True, linestyle=':', linewidth=0.5, alpha=0.7)

                if not idx_list:
                    # Empty anchor case: just draw frame/ticks, no data
                    continue

                # [FIXED] Logic for multiple spectra
                draw_items = []
                                
                if self.mode == "coadd":
                    # Average mode
                    got = self._cell_spectrum(idx_list)
                    if got:
                        draw_items.append((got[0], got[1], "black", 1.0, 0.8, True))
                else:
                    # All/Spaghetti mode: draw each
                    is_multi_cell = (len(idx_list) > 1)
                    line_color = "gray" if is_multi_cell else "black"
                    line_alpha = 0.6 if is_multi_cell else 1.0
                    line_lw = 0.5 if is_multi_cell else 0.8
                    for idx in idx_list:
                        got = self._cell_spectrum([idx])
                        if got:
                            draw_items.append((got[0], got[1], line_color, line_alpha, line_lw, False))

                if not draw_items: continue

                for row, y, color, alpha, lw, is_main in draw_items:
                    x_plot, y_plot = _process_spectrum(x, y, self.smooth_mode, self.smooth_width, self.box_downsample, self.box_policy)
                    
                    mx = np.ones_like(y_plot, dtype=bool)
                    if xlim is not None:
                        lo, hi = xlim
                        mx = (x_plot >= builtins.min(lo, hi)) & (x_plot <= builtins.max(lo, hi))

                    ax.step(x_plot, y_plot, where="mid", lw=lw, color=color, alpha=alpha)

                    yy = y_plot[mx]
                    yy = yy[np.isfinite(yy)]
                    if yy.size: 
                        ymins.append(float(np.min(yy)))
                        ymaxs.append(float(np.max(yy)))

                    # Only draw annotations (RMS, windows) for the "main" spectrum (coadd)
                    # or the last one if in "all" mode (to avoid clutter)
                    if is_main or (self.mode != "coadd" and row is draw_items[-1][0]):
                        ax.axhline(0, lw=0.3, color="gray", ls="--")
                        win_str = str(row.get("BSL_WINF", ""))

                        # Windows Display
                        if self.show_fit_windows and win_str and str(win_str).lower() not in ("nan", "none", ""):
                            try:
                                wins_x = fit_windows_xaxis(
                                    win_str,
                                    self.rest_hz_meta,
                                    self.rest_hz_user,
                                    self.xaxis,
                                    self.v_axis_disp,
                                    self.f_axis_ghz,
                                    self.ch_axis,
                                )
                                for x0_w, x1_w in wins_x:
                                    ax.axvspan(float(x0_w), float(x1_w), color="green", alpha=0.10, zorder=0)
                            except Exception:
                                pass

                        if self.annotate_rms:
                            # [FIXED] RMS calculation respects Fit Windows now
                            # 1. Determine base data (smoothed vs raw)
                            data_rms = y_plot if self.rms_on == "display" else y
                            
                            # 2. Build mask from X-limits and Windows
                            # (A) Start with X-limit mask
                            if xlim is not None:
                                lo, hi = xlim
                                # If display (already processed x_plot), mask on x_plot
                                # If raw (raw y), we need raw x (self.v_axis_disp etc).
                                # But `y` is raw, so we must use `x` (global raw axis)
                                if self.rms_on == "display":
                                    mask_base = (x_plot >= builtins.min(lo, hi)) & (x_plot <= builtins.max(lo, hi))
                                    x_for_win = x_plot
                                else:
                                    mask_base = (x >= builtins.min(lo, hi)) & (x <= builtins.max(lo, hi))
                                    x_for_win = x
                            else:
                                mask_base = np.ones_like(data_rms, dtype=bool)
                                x_for_win = x_plot if self.rms_on == "display" else x

                            # (B) Apply Windows Mask if available
                            mask_win = np.zeros_like(data_rms, dtype=bool)
                            has_win = False
                            if win_str and str(win_str).lower() not in ("nan", "none", ""):
                                try:
                                    wins_x_rms = fit_windows_xaxis(
                                        win_str,
                                        self.rest_hz_meta,
                                        self.rest_hz_user,
                                        self.xaxis,
                                        self.v_axis_disp,
                                        self.f_axis_ghz,
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
                                # Fallback: if no windows defined, use entire range (or xlim)
                                final_mask = mask_base

                            rms = _rms_from_mean(data_rms[final_mask])
                            
                            if np.isfinite(rms): 
                                # [MODIFIED] RMS text moved to upper-left, below coordinates
                                ax.text(0.03, 0.82, f"rms={rms:.3g}", transform=ax.transAxes, fontsize=6, va="top", ha="left")

                if self.show_coord_text and len(idx_list) > 0:
                    i0 = idx_list[0]
                    curr_x = float(self.sol.x_deg[i0])
                    curr_y = float(self.sol.y_deg[i0])
                    
                    dx_deg = curr_x - ref_x
                    if dx_deg > 180: dx_deg -= 360
                    elif dx_deg < -180: dx_deg += 360
                    
                    if self.sol.coord == "radec":
                        dx_arcsec = dx_deg * cos_dec_ref * 3600.0
                        dy_arcsec = (curr_y - ref_y) * 3600.0
                    else:
                        dx_arcsec = dx_deg * 3600.0
                        dy_arcsec = (curr_y - ref_y) * 3600.0

                    val_x = dx_arcsec * unit_scale
                    val_y = dy_arcsec * unit_scale
                    
                    # [MODIFIED] Format (dx, dy) and move to upper-left
                    label_text = f"({val_x:+.1f}{unit_label}, {val_y:+.1f}{unit_label})"
                    ax.text(0.03, 0.95, label_text, transform=ax.transAxes, fontsize=6, va="top", ha="left")

                # [ADDED] IDを右上に青字で全て表示 (はみ出し防止で4つごとに改行)
                if getattr(self, "show_id", False) and len(idx_list) > 0:
                    if "ORIG_ROW_ID" in self.st.table.columns:
                        orig_vals = self.st.table.iloc[idx_list]["ORIG_ROW_ID"].tolist()
                        chunks = [str(v) for v in orig_vals]
                    else:
                        chunks = [str(x) for x in idx_list]
                    lines = []
                    for i in range(0, len(chunks), 4):
                        lines.append(",".join(chunks[i:i+4]))
                    id_text = "id:\n" + "\n".join(lines)
                    ax.text(0.97, 0.95, id_text, transform=ax.transAxes, fontsize=5, va="top", ha="right", color="blue", alpha=0.7)
                    

                if xlim is not None: ax.set_xlim(*xlim)

        # Determine limits
        final_ylim = None
        if self.yrange is not None:
            final_ylim = self.yrange
        else:
            if ymins and ymaxs:
                ymin = builtins.min(ymins); ymax = builtins.max(ymaxs)
                yr = (ymax - ymin) or 1.0
                final_ylim = (ymin - 0.1 * yr, ymax + 0.1 * yr)
        
        # Apply limits to ALL active axes (including empty bottom-left)
        if final_ylim:
            for ax in self.axes_flat:
                # Only apply if axis is "on" (meaning it has data OR is the anchor)
                if ax.lines or (ax == self.axes[self.nrows-1, 0]):
                    ax.set_ylim(*final_ylim)
        
        # Helper for tick format of freq
        if self.xaxis == "freq":
             for ax in self.axes_flat:
                 if ax.lines or (ax == self.axes[self.nrows-1, 0]):
                     ax.ticklabel_format(axis="x", style="plain", useOffset=False)
                     ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))

        if hasattr(self.fig.canvas, "draw_idle"):
            self.fig.canvas.draw_idle()
        else:
            self.fig.canvas.draw()

    def on_key(self, event):
        if event is None or event.key is None: return
        key = event.key
        if key == "n":
            self.page = (self.page + 1) % builtins.max(self.n_pages, 1); self.update_plot()
        elif key == "p":
            self.page = (self.page - 1 + builtins.max(self.n_pages, 1)) % builtins.max(self.n_pages, 1); self.update_plot()
        elif key == "x":
            self._cycle_xaxis(); self.update_plot()
        elif key == "t":
            self._toggle_scale(); self.update_plot()
        elif key == "m":
            # [MODIFIED] Removed "running" from cycle
            order = ["none", "boxcar"]
            i = order.index(self.smooth_mode) if self.smooth_mode in order else 0
            self.smooth_mode = order[(i + 1) % len(order)]
            self.update_plot()
        elif key in ("[", "]"):
            step = -1 if key == "[" else 1
            w = int(self.smooth_width) + step
            w = builtins.max(1, w)
            self.smooth_width = w
            self.update_plot()
        elif key == "d":
            self.box_downsample = not self.box_downsample; self.update_plot()
        elif key in ("r", "R"):
            self.rms_on = "display" if self.rms_on == "raw" else "raw"; self.update_plot()
        elif key in ("b", "B"):
            self.show_fit_windows = not self.show_fit_windows; self.update_plot()
        elif key in ("w", "W"):
            self.show_fit_windows = not self.show_fit_windows; self.update_plot()
        elif key in ("g", "G"):
            if self.show_grid == "auto": self.show_grid = False
            elif self.show_grid is False: self.show_grid = True
            else: self.show_grid = "auto"
            self.update_plot()
        elif key in ("a", "A"):
            self.square_aspect = not self.square_aspect
            self._apply_square_layout(); self.update_plot()
        elif key in ("h", "H"):
            self._print_help()
        elif key == "s":
            try:
                save_dir = Path(getattr(self, "save_dir", "."))
                save_dir.mkdir(parents=True, exist_ok=True)
                prefix = str(getattr(self, "save_prefix", "grid"))
                fn = f"{prefix}_pg{int(self.page)+1:03d}_{self.xaxis}_{self.display_scale[:2]}.pdf"
                out = save_dir / fn
                self.fig.savefig(out, bbox_inches="tight", pad_inches=0.5)
                print(f"Saved: {out}")
            except Exception as e:
                print(f"Save failed: {e}")
        elif key == "q":
            plt.close(self.fig)
