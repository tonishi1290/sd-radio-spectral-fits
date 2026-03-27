from __future__ import annotations
import builtins
import logging
from pathlib import Path
from dataclasses import dataclass
from typing import Sequence
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree
from .config import GridInput
from .gridder import create_grid_input, _resolve_otf_processing_policy, _resolve_projection_reference_for_scantables
from ..ranges import parse_windows, window_to_mask


C_KMS = 299792.458


_FORMAL_STREAM_COLS_BW = ('FDNUM', 'IFNUM', 'PLNUM')
_FALLBACK_STREAM_COLS_BW = ('BEAM', 'POL', 'IF', 'ARRAYID', 'FEED', 'STREAM', 'SIDEBAND', 'BAND', 'RX', 'SAMPLER', 'BOARD', 'XFFTSBOARD')


def _normalize_otf_scan_png_path(path_like, *, default_path: str | None = None):
    if path_like in (None, False):
        return None
    if path_like is True:
        return str(default_path if default_path is not None else 'otf_scan_region.png')
    return str(path_like)


def _indexed_otf_scan_png_path(path_like, index: int, total: int, *, default_path: str | None = None):
    resolved = _normalize_otf_scan_png_path(path_like, default_path=default_path)
    if resolved is None:
        return None
    path = Path(str(resolved))
    if total <= 1:
        return str(path)
    return str(path.with_name(f'{path.stem}.{index:03d}{path.suffix}'))


def _resolve_otf_scan_png_sequence(otf_scan_png, count: int, *, default_path: str | None = None):
    if count <= 0:
        return []
    if isinstance(otf_scan_png, (list, tuple)):
        if len(otf_scan_png) != count:
            raise ValueError('otf_scan_png sequence length must match the number of inputs.')
        out = []
        for i, item in enumerate(otf_scan_png):
            if item is True:
                out.append(_indexed_otf_scan_png_path(item, i, count, default_path=default_path))
            else:
                out.append(_normalize_otf_scan_png_path(item, default_path=default_path))
        return out
    return [
        _indexed_otf_scan_png_path(otf_scan_png, i, count, default_path=default_path)
        for i in range(count)
    ]


def _ensure_not_all_empty_inputs(dataset_or_input_data):
    items = dataset_or_input_data if isinstance(dataset_or_input_data, (list, tuple)) else [dataset_or_input_data]
    lengths = []
    for idx, item in enumerate(items):
        if isinstance(item, GridInput):
            lengths.append(int(len(np.asarray(item.x))))
            continue
        table = getattr(item, 'table', None)
        if table is None:
            continue
        lengths.append(int(len(table)))
    if lengths and not any(n > 0 for n in lengths):
        raise ValueError('All input scantables/GridInputs are empty; at least one non-empty input is required for basketweave.')

def _input_row_count_for_basketweave(item) -> int | None:
    if isinstance(item, GridInput):
        return int(len(np.asarray(item.x)))
    table = getattr(item, 'table', None)
    if table is None:
        return None
    return int(len(table))

def _filter_nonempty_inputs_for_basketweave(dataset_or_input_data):
    items = dataset_or_input_data if isinstance(dataset_or_input_data, (list, tuple)) else [dataset_or_input_data]
    kept = []
    skipped = []
    for idx, item in enumerate(items):
        nrow = _input_row_count_for_basketweave(item)
        if nrow is None or int(nrow) > 0:
            kept.append((idx, item))
        else:
            skipped.append(int(idx))
    if skipped:
        import warnings
        warnings.warn(
            f'basketweave is skipping empty inputs at indices {skipped}.',
            RuntimeWarning,
            stacklevel=2,
        )
    return kept, skipped



def _ensure_safe_basketweave_input(dataset_or_input_data, *, otf_input_state=None, otf_scan_region=None):
    state, _ = _resolve_otf_processing_policy(otf_input_state=otf_input_state, otf_scan_region=otf_scan_region)
    if state != "use_existing_labels":
        return
    items = dataset_or_input_data if isinstance(dataset_or_input_data, (list, tuple)) else [dataset_or_input_data]
    required_cols = {"SCAN", "IS_TURN", "FDNUM", "IFNUM", "PLNUM"}
    for idx, item in enumerate(items):
        if isinstance(item, GridInput):
            continue
        table = getattr(item, 'table', None)
        cols = set(getattr(table, 'columns', [])) if table is not None else set()
        missing = sorted(required_cols - cols)
        if missing:
            raise ValueError(
                "otf_input_state='use_existing_labels' requires SCAN, IS_TURN, FDNUM, IFNUM, and PLNUM in every raw scantable "
                f"before basketweave. Input {idx} is missing: {', '.join(str(c) for c in missing)}."
            )

def _normalize_unit_local(unit_str: str) -> str:
    u = str(unit_str).strip().lower()
    if 'hz' in u:
        return 'hz'
    if u in ('m/s', 'm s-1', 'ms-1', 'meter/sec', 'm/sec'):
        return 'm/s'
    if u in ('km/s', 'km s-1', 'kms-1', 'kilometer/sec', 'km/sec'):
        return 'km/s'
    return u

def _freq_scale_to_hz_local(unit_str: str) -> float:
    u = str(unit_str).strip().lower()
    if u in ('', 'hz'):
        return 1.0
    if u == 'khz':
        return 1.0e3
    if u == 'mhz':
        return 1.0e6
    if u == 'ghz':
        return 1.0e9
    if u == 'thz':
        return 1.0e12
    return 1.0

def _get_restfreq_local(meta: dict) -> float:
    for k in ('RESTFRQ', 'RESTFREQ', 'rest_hz', 'restfrq_hz', 'restfreq_hz'):
        if k in meta and meta[k] not in (None, ''):
            val = float(meta[k])
            if np.isfinite(val) and val > 0:
                return val
    raise ValueError('RESTFRQ/RESTFREQ is required to convert spectral axis to velocity.')

def _radio_velocity_kms_from_freq(freq_hz: np.ndarray, rest_hz: float) -> np.ndarray:
    return C_KMS * (1.0 - np.asarray(freq_hz, dtype=float) / float(rest_hz))

def _infer_nchan_from_dataset(dataset) -> int:
    data = np.asarray(dataset.data)
    if data.ndim == 1:
        return int(data.shape[0])
    if data.ndim >= 2:
        return int(data.shape[1])
    raise ValueError(f'Cannot infer NCHAN from dataset.data with shape={data.shape}')

def _dataset_meta_like(dataset) -> dict:
    meta = getattr(dataset, 'meta', None)
    if meta is None:
        raise ValueError('Input dataset does not provide .meta; cannot infer spectral axis automatically.')
    if isinstance(meta, dict):
        return dict(meta)
    try:
        return dict(meta)
    except Exception as exc:
        raise ValueError('Input dataset .meta is not dict-like; cannot infer spectral axis automatically.') from exc

def _validate_non_topocentric_specsys(meta: dict) -> None:
    specsys = str(meta.get('SPECSYS', meta.get('SSYSOBS', ''))).strip().upper()
    if specsys == 'TOPOCENT':
        raise ValueError('basketweave automatic velocity-axis inference does not support SPECSYS=TOPOCENT. Regrid/coadd to a non-topocentric spectral frame first, or pass v_axis/channel_mask explicitly.')

def _velocity_axis_kms_from_meta_simple(meta: dict, *, nchan: int) -> np.ndarray:
    _validate_non_topocentric_specsys(meta)
    raw_ctype = str(meta.get('CTYPE1', 'FREQ')).strip()
    ctype = raw_ctype.upper()
    raw_cunit = str(meta.get('CUNIT1', ''))
    unit_norm = _normalize_unit_local(raw_cunit)
    if 'CRVAL1' not in meta or 'CDELT1' not in meta:
        raise ValueError('CRVAL1/CDELT1 are required for automatic velocity-axis inference.')
    crval = float(meta['CRVAL1'])
    cdelt = float(meta['CDELT1'])
    crpix = float(meta.get('CRPIX1', 1.0))
    i = np.arange(int(nchan), dtype=float)
    axis_native = crval + (i + 1.0 - crpix) * cdelt
    if ctype.startswith('FREQ') or unit_norm == 'hz':
        freq_hz = axis_native * _freq_scale_to_hz_local(raw_cunit)
        rest_hz = _get_restfreq_local(meta)
        return _radio_velocity_kms_from_freq(freq_hz, rest_hz)
    scale_to_kms = 0.001 if unit_norm == 'm/s' else 1.0
    v_axis_kms = axis_native * scale_to_kms
    if ctype.startswith('VRAD'):
        return v_axis_kms
    if ctype.startswith('VOPT'):
        rest_hz = _get_restfreq_local(meta)
        freq_hz = rest_hz / (v_axis_kms / C_KMS + 1.0)
        return _radio_velocity_kms_from_freq(freq_hz, rest_hz)
    if ctype.startswith('VELO'):
        rest_hz = _get_restfreq_local(meta)
        beta = np.clip(v_axis_kms / C_KMS, -0.999, 0.999)
        freq_hz = rest_hz * np.sqrt((1.0 - beta) / (1.0 + beta))
        return _radio_velocity_kms_from_freq(freq_hz, rest_hz)
    if unit_norm in ('m/s', 'km/s'):
        return v_axis_kms
    raise ValueError(f'Unsupported spectral axis for automatic velocity inference: CTYPE1={raw_ctype!r}, CUNIT1={raw_cunit!r}')

def _infer_velocity_axis_kms_from_dataset(dataset) -> np.ndarray:
    if isinstance(dataset, GridInput):
        raise ValueError('Automatic velocity-axis inference is not available for GridInput input. Pass v_axis or channel_mask explicitly.')
    meta = _dataset_meta_like(dataset)
    nchan = _infer_nchan_from_dataset(dataset)
    return np.asarray(_velocity_axis_kms_from_meta_simple(meta, nchan=nchan), dtype=float)

def _infer_velocity_axis_kms_for_input(dataset_or_input_data) -> np.ndarray | None:
    if isinstance(dataset_or_input_data, (list, tuple)):
        _ensure_not_all_empty_inputs(dataset_or_input_data)
        _validate_grid_input_contract(dataset_or_input_data)
        axes = []
        for item in dataset_or_input_data:
            if isinstance(item, GridInput):
                raise ValueError('Automatic velocity-axis inference for a list/tuple input is only supported when all elements are scantable-like objects. For GridInput elements, pass v_axis or channel_mask explicitly.')
            axes.append(_infer_velocity_axis_kms_from_dataset(item))
        if len(axes) == 0:
            return None
        ref = np.asarray(axes[0], dtype=float)
        for k, ax in enumerate(axes[1:], start=1):
            ax = np.asarray(ax, dtype=float)
            if ax.shape != ref.shape or not np.allclose(ax, ref, rtol=0.0, atol=1.0e-12, equal_nan=True):
                raise ValueError(f'Automatic velocity-axis inference requires all inputs to share the same spectral axis. Input 0 and input {k} differ; pass v_axis or channel_mask explicitly.')
        return ref
    if isinstance(dataset_or_input_data, GridInput):
        return None
    return _infer_velocity_axis_kms_from_dataset(dataset_or_input_data)

@dataclass
class _PairSelectionResult:
    i_idx: np.ndarray
    j_idx: np.ndarray
    scan_i: np.ndarray
    scan_j: np.ndarray
    n_pairs_total: int
    n_pairs_cross_scan: int
    n_pairs_used: int
    n_pairs_direction_rejected: int
    direction_filter_applied: bool
    direction_filter_fallback_used: bool

@dataclass
class _DirectionInferenceResult:
    uniq_scan_ids: np.ndarray
    angles_deg: np.ndarray
    straightness: np.ndarray
    from_scan_dir: np.ndarray

@dataclass
class _SamplingScaleResult:
    along_step_arcsec: float
    cross_step_arcsec: float

def _safe_bool_array(values, default: bool=False) -> np.ndarray:
    """Convert mixed boolean-like values into a strict bool ndarray."""
    a = np.asarray(values)
    if a.dtype.kind == 'b':
        return a.astype(bool, copy=False)
    out = np.full(a.shape, bool(default), dtype=bool)
    if a.dtype.kind in 'iuf':
        finite = np.isfinite(a)
        out[finite] = a[finite] != 0
        return out
    flat = a.astype(object, copy=False).ravel()
    flat_out = out.ravel()
    true_vals = {'true', 't', '1', 'yes', 'y', 'on'}
    false_vals = {'false', 'f', '0', 'no', 'n', 'off', '', 'nan', 'none', 'null', '<na>'}
    for i, v in enumerate(flat):
        if v is None:
            flat_out[i] = bool(default)
            continue
        s = str(v).strip().lower()
        if s in true_vals:
            flat_out[i] = True
        elif s in false_vals:
            flat_out[i] = False
        else:
            flat_out[i] = bool(default)
    return out

def _safe_scan_ids_int64(values) -> np.ndarray:
    """Convert scan ids to int64, mapping invalid entries to -1."""
    arr = np.asarray(values)
    if arr.size == 0:
        return np.array([], dtype=np.int64)
    num = np.empty(arr.shape, dtype=np.float64)
    num.fill(np.nan)
    try:
        num = arr.astype(np.float64)
    except (TypeError, ValueError):
        flat_src = arr.astype(object, copy=False).ravel()
        flat_dst = num.ravel()
        for i, v in enumerate(flat_src):
            try:
                flat_dst[i] = float(v)
            except (TypeError, ValueError):
                flat_dst[i] = np.nan
    out = np.full(num.shape, -1, dtype=np.int64)
    finite = np.isfinite(num)
    out[finite] = num[finite].astype(np.int64)
    return out

def _zero_offsets_from_scan_ids(scan_ids_all: np.ndarray) -> np.ndarray:
    max_scan = int(np.max(scan_ids_all)) if np.any(scan_ids_all >= 0) else -1
    return np.zeros(max_scan + 1, dtype=float)

def _sorted_pairs_array(pairs_arr) -> np.ndarray:
    arr = np.asarray(pairs_arr, dtype=np.int64)
    if arr.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    arr = arr.reshape(-1, 2)
    order = np.lexsort((arr[:, 1], arr[:, 0]))
    return arr[order]

def _base_geometry_mask(input_data: GridInput, *, x_all: np.ndarray, y_all: np.ndarray, scan_ids_all: np.ndarray) -> np.ndarray:
    good = np.ones_like(x_all, dtype=bool)
    if input_data.flag is not None:
        good &= np.asarray(input_data.flag) > 0
    if input_data.is_turnaround is not None:
        good &= ~_safe_bool_array(input_data.is_turnaround, default=False)
    good &= np.isfinite(x_all) & np.isfinite(y_all)
    good &= scan_ids_all >= 0
    return good

def _resolve_linefree_velocity_windows(*, linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None, v_windows_kms: list[str] | list[tuple[float, float]] | None) -> list[str] | list[tuple[float, float]] | None:
    if linefree_velocity_windows_kms is not None and v_windows_kms is not None:
        raise ValueError("Specify only one of 'linefree_velocity_windows_kms' or the deprecated alias 'v_windows_kms'.")
    if linefree_velocity_windows_kms is not None:
        return linefree_velocity_windows_kms
    if v_windows_kms is not None:
        logging.warning("'v_windows_kms' is deprecated. Use 'linefree_velocity_windows_kms' instead.")
        return v_windows_kms
    return None

def _build_channel_mask(*, spec: np.ndarray, v_axis: np.ndarray | None, linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None, channel_mask: np.ndarray | None) -> np.ndarray | None:
    resolved_mask = channel_mask
    if linefree_velocity_windows_kms is not None:
        if v_axis is None:
            raise ValueError('v_axis must be provided when linefree_velocity_windows_kms is specified.')
        if len(linefree_velocity_windows_kms) > 0 and isinstance(linefree_velocity_windows_kms[0], str):
            windows = parse_windows(linefree_velocity_windows_kms)
        else:
            windows = linefree_velocity_windows_kms
        resolved_mask = window_to_mask(np.asarray(v_axis, dtype=float), windows)
        if not np.any(resolved_mask):
            logging.error('Basket-weave aborted: No channels found in windows %s.', linefree_velocity_windows_kms)
            logging.error('(Available velocity range is approx %.1f to %.1f km/s)', np.nanmin(v_axis), np.nanmax(v_axis))
            raise ValueError('Invalid linefree_velocity_windows_kms: No matching channels. Please check your velocity range.')
    if resolved_mask is not None:
        resolved_mask = np.asarray(resolved_mask, dtype=bool)
        if spec.ndim != 2:
            raise ValueError('channel_mask can only be used when input_data.spec is 2D.')
        if resolved_mask.ndim != 1 or resolved_mask.shape[0] != spec.shape[1]:
            raise ValueError(f'channel_mask must have shape ({spec.shape[1]},), got {resolved_mask.shape}')
        if not np.any(resolved_mask):
            raise ValueError('channel_mask selected zero channels.')
    return resolved_mask

def _aggregate_spectrum_for_offsets(spec: np.ndarray, channel_mask: np.ndarray | None) -> np.ndarray:
    if spec.ndim == 2:
        if channel_mask is not None:
            return np.nanmean(spec[:, channel_mask], axis=1)
        logging.warning("Basket-weaving is averaging across ALL channels. Strong emission lines might bias the baseline offset estimation. Consider specifying 'linefree_velocity_windows_kms' for better accuracy.")
        return np.nanmean(spec, axis=1)
    return np.asarray(spec, dtype=float)

def _orientation_difference_deg(a_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    """Return orientation difference in [0, 90] deg for 180-deg periodic angles."""
    diff = np.abs(a_deg - b_deg)
    diff = np.mod(diff, 180.0)
    diff = np.where(diff > 90.0, 180.0 - diff, diff)
    return diff

def _normalize_scan_dir_value_to_angle_deg(value) -> float:
    """Normalize one scan_dir-like value into an orientation angle in [0, 180)."""
    if value is None:
        return np.nan
    if isinstance(value, bytes):
        try:
            value = value.decode('utf-8', errors='ignore')
        except Exception:
            return np.nan
    if isinstance(value, (bool, np.bool_)):
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        val = float(value)
        if not np.isfinite(val):
            return np.nan
        return float(val % 180.0)
    s = str(value).strip().lower()
    if s == '':
        return np.nan
    aliases_x = {'x', 'scan_x', 'xscan', 'along_x', '+x', '-x', 'lon', 'longitude', 'ra', 'az', 'horizontal', 'h', 'row', 'rows'}
    aliases_y = {'y', 'scan_y', 'yscan', 'along_y', '+y', '-y', 'lat', 'latitude', 'dec', 'el', 'vertical', 'v', 'col', 'cols', 'column', 'columns'}
    if s in aliases_x:
        return 0.0
    if s in aliases_y:
        return 90.0
    for prefix in ('scan_', 'dir_', 'axis_'):
        if s.startswith(prefix):
            base = s[len(prefix):]
            if base in aliases_x:
                return 0.0
            if base in aliases_y:
                return 90.0
    try:
        val = float(s)
    except ValueError:
        return np.nan
    if not np.isfinite(val):
        return np.nan
    return float(val % 180.0)

def _infer_compact_scan_directions_from_scan_dir(scan_dir, inv_scan: np.ndarray, *, uniq_scan_ids: np.ndarray, min_fraction: float=0.6, angle_consistency_tol_deg: float=20.0) -> tuple[np.ndarray, np.ndarray]:
    """Infer per-scan orientation from input_data.scan_dir when available."""
    num_scans = int(len(uniq_scan_ids))
    angles_deg = np.full(num_scans, np.nan, dtype=float)
    confidence = np.full(num_scans, np.nan, dtype=float)
    if scan_dir is None:
        return (angles_deg, confidence)
    scan_dir_arr = np.asarray(scan_dir, dtype=object)
    if scan_dir_arr.ndim != 1:
        scan_dir_arr = scan_dir_arr.ravel()
    normalized = np.array([_normalize_scan_dir_value_to_angle_deg(v) for v in scan_dir_arr], dtype=float)
    for compact_scan_index in range(num_scans):
        idx = np.flatnonzero(inv_scan == compact_scan_index)
        if idx.size == 0:
            continue
        vals = normalized[idx]
        finite = np.isfinite(vals)
        if not np.any(finite):
            continue
        vals = vals[finite] % 180.0
        if vals.size == 1:
            angles_deg[compact_scan_index] = float(vals[0])
            confidence[compact_scan_index] = 1.0
            continue
        bins = np.array([0.0, 90.0], dtype=float)
        diff_to_bins = np.stack([_orientation_difference_deg(vals, b) for b in bins], axis=1)
        nearest_bin_idx = np.argmin(diff_to_bins, axis=1)
        counts = np.bincount(nearest_bin_idx, minlength=2)
        dominant_idx = int(np.argmax(counts))
        dominant_fraction = float(counts[dominant_idx] / vals.size)
        dominant_bin = bins[dominant_idx]
        dominant_vals = vals[nearest_bin_idx == dominant_idx]
        dominant_diff = _orientation_difference_deg(dominant_vals, dominant_bin)
        if dominant_fraction >= float(min_fraction) and np.nanmax(dominant_diff) <= float(angle_consistency_tol_deg):
            angles_deg[compact_scan_index] = float(dominant_bin)
            confidence[compact_scan_index] = dominant_fraction
            continue
        circ = np.exp(1j * np.deg2rad(2.0 * vals))
        mean_vec = np.nanmean(circ)
        if np.isfinite(mean_vec.real) and np.isfinite(mean_vec.imag):
            angle = 0.5 * np.degrees(np.angle(mean_vec))
            angle = float(angle % 180.0)
            diffs = _orientation_difference_deg(vals, angle)
            if np.nanmax(diffs) <= float(angle_consistency_tol_deg):
                angles_deg[compact_scan_index] = angle
                confidence[compact_scan_index] = float(np.abs(mean_vec))
    return (angles_deg, confidence)

def _infer_compact_scan_directions(x: np.ndarray, y: np.ndarray, inv_scan: np.ndarray, *, uniq_scan_ids: np.ndarray, scan_dir=None, min_points_per_scan: int=3, min_straightness: float=0.85) -> _DirectionInferenceResult:
    """
    Infer scan orientation angle [deg] for each compact scan index.

    Preference order:
      1. per-sample scan_dir, when it can be interpreted consistently within a scan
      2. geometric inference from x/y sampling

    The angle is defined modulo 180 deg. 0 deg means approximately along +x/-x,
    and 90 deg means approximately along +y/-y.
    """
    num_scans = int(len(uniq_scan_ids))
    angles_deg = np.full(num_scans, np.nan, dtype=float)
    straightness = np.full(num_scans, np.nan, dtype=float)
    from_scan_dir = np.zeros(num_scans, dtype=bool)
    dir_angles, dir_confidence = _infer_compact_scan_directions_from_scan_dir(scan_dir, inv_scan, uniq_scan_ids=uniq_scan_ids)
    finite_dir = np.isfinite(dir_angles)
    angles_deg[finite_dir] = dir_angles[finite_dir]
    straightness[finite_dir] = dir_confidence[finite_dir]
    from_scan_dir[finite_dir] = True
    for compact_scan_index in range(num_scans):
        if np.isfinite(angles_deg[compact_scan_index]):
            continue
        idx = np.flatnonzero(inv_scan == compact_scan_index)
        if idx.size < 2:
            continue
        xs = np.asarray(x[idx], dtype=float)
        ys = np.asarray(y[idx], dtype=float)
        finite = np.isfinite(xs) & np.isfinite(ys)
        xs = xs[finite]
        ys = ys[finite]
        if xs.size < builtins.max(2, min_points_per_scan):
            continue
        pts = np.column_stack([xs, ys])
        pts -= np.nanmean(pts, axis=0, keepdims=True)
        try:
            _, svals, vh = np.linalg.svd(pts, full_matrices=False)
            principal = vh[0]
            if svals.size >= 2:
                denom = float(svals[0] ** 2 + svals[1] ** 2)
                strength = float(svals[0] ** 2 / denom) if denom > 0 else np.nan
            else:
                strength = 1.0
        except np.linalg.LinAlgError:
            principal = None
            strength = np.nan
        if principal is None or not np.all(np.isfinite(principal)):
            dx = float(xs[-1] - xs[0])
            dy = float(ys[-1] - ys[0])
            if not np.isfinite(dx) or not np.isfinite(dy) or (dx == 0.0 and dy == 0.0):
                continue
            principal = np.array([dx, dy], dtype=float)
        angle = float(np.degrees(np.arctan2(principal[1], principal[0])) % 180.0)
        if np.isfinite(strength) and strength < min_straightness:
            pass
        angles_deg[compact_scan_index] = angle
        straightness[compact_scan_index] = strength
    return _DirectionInferenceResult(uniq_scan_ids=np.asarray(uniq_scan_ids, dtype=np.int64), angles_deg=angles_deg, straightness=straightness, from_scan_dir=from_scan_dir)

def _select_pair_mask_by_direction(scan_i: np.ndarray, scan_j: np.ndarray, scan_angles_deg: np.ndarray | None, *, orthogonality_tolerance_deg: float, fallback_to_all_cross_scan_pairs: bool) -> tuple[np.ndarray, int, bool, bool]:
    if scan_angles_deg is None:
        used_mask = np.ones_like(scan_i, dtype=bool)
        return (used_mask, 0, False, False)
    ang_i = np.asarray(scan_angles_deg[scan_i], dtype=float)
    ang_j = np.asarray(scan_angles_deg[scan_j], dtype=float)
    finite = np.isfinite(ang_i) & np.isfinite(ang_j)
    if not np.any(finite):
        if fallback_to_all_cross_scan_pairs:
            return (np.ones_like(scan_i, dtype=bool), 0, True, True)
        return (np.zeros_like(scan_i, dtype=bool), int(len(scan_i)), True, False)
    diff = _orientation_difference_deg(ang_i, ang_j)
    orth_mask = finite & (np.abs(diff - 90.0) <= float(orthogonality_tolerance_deg))
    if np.any(orth_mask):
        rejected = int(len(scan_i) - np.count_nonzero(orth_mask))
        return (orth_mask, rejected, True, False)
    if fallback_to_all_cross_scan_pairs:
        return (np.ones_like(scan_i, dtype=bool), 0, True, True)
    return (np.zeros_like(scan_i, dtype=bool), int(len(scan_i)), True, False)

def _find_cross_scan_pairs(x: np.ndarray, y: np.ndarray, inv_scan: np.ndarray, *, search_radius_arcsec: float, scan_angles_deg: np.ndarray | None=None, cross_direction_only: bool=True, orthogonality_tolerance_deg: float=30.0, fallback_to_all_cross_scan_pairs: bool=True) -> _PairSelectionResult:
    tree = cKDTree(np.c_[x, y])
    try:
        pairs_arr = tree.query_pairs(r=float(search_radius_arcsec), output_type='ndarray')
    except TypeError:
        pairs_set = tree.query_pairs(r=float(search_radius_arcsec))
        if not pairs_set:
            empty = np.empty(0, dtype=np.int64)
            return _PairSelectionResult(i_idx=empty, j_idx=empty, scan_i=empty, scan_j=empty, n_pairs_total=0, n_pairs_cross_scan=0, n_pairs_used=0, n_pairs_direction_rejected=0, direction_filter_applied=False, direction_filter_fallback_used=False)
        pairs_arr = np.array(list(pairs_set), dtype=np.int64)
    pairs_arr = _sorted_pairs_array(pairs_arr)
    n_total = int(len(pairs_arr))
    if n_total == 0:
        empty = np.empty(0, dtype=np.int64)
        return _PairSelectionResult(i_idx=empty, j_idx=empty, scan_i=empty, scan_j=empty, n_pairs_total=0, n_pairs_cross_scan=0, n_pairs_used=0, n_pairs_direction_rejected=0, direction_filter_applied=False, direction_filter_fallback_used=False)
    i_idx = pairs_arr[:, 0]
    j_idx = pairs_arr[:, 1]
    scan_i = inv_scan[i_idx]
    scan_j = inv_scan[j_idx]
    cross_scan_mask = scan_i != scan_j
    i_idx = i_idx[cross_scan_mask]
    j_idx = j_idx[cross_scan_mask]
    scan_i = scan_i[cross_scan_mask]
    scan_j = scan_j[cross_scan_mask]
    n_cross = int(len(i_idx))
    if n_cross == 0:
        empty = np.empty(0, dtype=np.int64)
        return _PairSelectionResult(i_idx=empty, j_idx=empty, scan_i=empty, scan_j=empty, n_pairs_total=n_total, n_pairs_cross_scan=0, n_pairs_used=0, n_pairs_direction_rejected=0, direction_filter_applied=False, direction_filter_fallback_used=False)
    if cross_direction_only:
        used_mask, n_rejected, direction_filter_applied, fallback_used = _select_pair_mask_by_direction(scan_i, scan_j, scan_angles_deg, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs)
    else:
        used_mask = np.ones_like(scan_i, dtype=bool)
        n_rejected = 0
        direction_filter_applied = False
        fallback_used = False
    i_idx = i_idx[used_mask]
    j_idx = j_idx[used_mask]
    scan_i = scan_i[used_mask]
    scan_j = scan_j[used_mask]
    return _PairSelectionResult(i_idx=i_idx, j_idx=j_idx, scan_i=scan_i, scan_j=scan_j, n_pairs_total=n_total, n_pairs_cross_scan=n_cross, n_pairs_used=int(len(i_idx)), n_pairs_direction_rejected=int(n_rejected), direction_filter_applied=bool(direction_filter_applied), direction_filter_fallback_used=bool(fallback_used))

def _estimate_sampling_scales_from_arrays(x: np.ndarray, y: np.ndarray, inv_scan: np.ndarray) -> _SamplingScaleResult:
    """Estimate representative along-scan and cross-scan spacings [arcsec]."""
    along_dists_list: list[np.ndarray] = []
    uniq_compact = np.unique(inv_scan)
    for compact_scan_index in uniq_compact:
        idx = np.flatnonzero(inv_scan == compact_scan_index)
        if idx.size < 2:
            continue
        dx = np.diff(x[idx])
        dy = np.diff(y[idx])
        dist = np.hypot(dx, dy)
        finite = np.isfinite(dist) & (dist > 0)
        if np.any(finite):
            along_dists_list.append(dist[finite])
    along_step = np.nan
    if along_dists_list:
        along_step = float(np.nanmedian(np.concatenate(along_dists_list)))
    cross_step = np.nan
    coords = np.c_[x, y]
    n = len(coords)
    if n >= 2 and uniq_compact.size >= 2:
        tree = cKDTree(coords)
        max_k = builtins.min(64, n)
        k = builtins.min(8, n)
        cross_nn = np.full(n, np.nan, dtype=float)
        while True:
            dists, indices = tree.query(coords, k=k)
            dists = np.asarray(dists, dtype=float)
            indices = np.asarray(indices, dtype=np.int64)
            if dists.ndim == 1:
                dists = dists[:, np.newaxis]
                indices = indices[:, np.newaxis]
            unresolved = ~np.isfinite(cross_nn)
            if np.any(unresolved):
                unresolved_idx = np.flatnonzero(unresolved)
                for row in unresolved_idx:
                    for dist, nbr in zip(dists[row], indices[row]):
                        if nbr < 0 or nbr == row:
                            continue
                        if not np.isfinite(dist) or dist <= 0:
                            continue
                        if inv_scan[nbr] == inv_scan[row]:
                            continue
                        cross_nn[row] = float(dist)
                        break
            if np.all(np.isfinite(cross_nn)) or k >= max_k:
                break
            k = builtins.min(max_k, k * 2)
        finite_cross = np.isfinite(cross_nn) & (cross_nn > 0)
        if np.any(finite_cross):
            cross_step = float(np.nanmedian(cross_nn[finite_cross]))
    return _SamplingScaleResult(along_step_arcsec=along_step, cross_step_arcsec=cross_step)

def estimate_basket_weave_sampling_scales_arcsec(input_data: GridInput) -> _SamplingScaleResult:
    """Estimate representative along-scan and cross-scan spacings [arcsec] from GridInput."""
    if input_data.scan_id is None:
        raise ValueError('scan_id must be provided in GridInput to estimate sampling scales.')
    x_all = np.asarray(input_data.x, dtype=float)
    y_all = np.asarray(input_data.y, dtype=float)
    scan_ids_all = _safe_scan_ids_int64(input_data.scan_id)
    good = _base_geometry_mask(input_data, x_all=x_all, y_all=y_all, scan_ids_all=scan_ids_all)
    if np.count_nonzero(good) < 2:
        return _SamplingScaleResult(along_step_arcsec=np.nan, cross_step_arcsec=np.nan)
    x = x_all[good]
    y = y_all[good]
    scan_ids = scan_ids_all[good]
    _, inv_scan = np.unique(scan_ids, return_inverse=True)
    return _estimate_sampling_scales_from_arrays(x, y, inv_scan)

def _warn_if_search_radius_exceeds_cross_scan_spacing(*, search_radius_arcsec: float, cross_step_arcsec: float, emit_warning: bool=True) -> bool:
    exceeds = bool(np.isfinite(search_radius_arcsec) and np.isfinite(cross_step_arcsec) and (cross_step_arcsec > 0) and (search_radius_arcsec >= cross_step_arcsec))
    if exceeds and emit_warning:
        logging.warning('search_radius_arcsec=%.3f is >= estimated cross-scan spacing %.3f arcsec. This may admit near-parallel scan pairs and bias basket-weave offsets.', search_radius_arcsec, cross_step_arcsec)
    return exceeds

def _compute_constraint_residual_rms(rhs_diff: np.ndarray, scan_i: np.ndarray, scan_j: np.ndarray, offsets_compact: np.ndarray, *, weights: np.ndarray | None=None) -> tuple[float | None, float | None]:
    if len(scan_i) == 0:
        return (None, None)
    before = np.asarray(rhs_diff, dtype=float)
    after = before - np.asarray(offsets_compact[scan_i] - offsets_compact[scan_j], dtype=float)

    def _wrms(values: np.ndarray, w: np.ndarray | None) -> float | None:
        finite = np.isfinite(values)
        if w is not None:
            ww = np.asarray(w, dtype=float)
            finite &= np.isfinite(ww) & (ww > 0)
        if not np.any(finite):
            return None
        if w is None:
            return float(np.sqrt(np.nanmean(values[finite] ** 2)))
        ww = np.asarray(w, dtype=float)[finite]
        vv = values[finite]
        denom = float(np.sum(ww))
        if denom <= 0:
            return None
        return float(np.sqrt(np.sum(ww * vv * vv) / denom))
    before_rms = _wrms(before, weights)
    after_rms = _wrms(after, weights)
    return (before_rms, after_rms)

def _build_scan_connectivity_diagnostics(scan_i: np.ndarray, scan_j: np.ndarray, uniq_scan_ids: np.ndarray) -> tuple[int, list[int], list[int], list[tuple[int, int]]]:
    num_scans = int(len(uniq_scan_ids))
    if num_scans == 0:
        return (0, [], [], [])
    if len(scan_i) == 0:
        isolated = [int(v) for v in np.asarray(uniq_scan_ids, dtype=np.int64).tolist()]
        return (num_scans, [1] * num_scans, isolated, [])
    edges_arr = np.stack([scan_i, scan_j], axis=1).astype(np.int64, copy=False)
    edges_arr = np.sort(edges_arr, axis=1)
    edges_arr = np.unique(edges_arr, axis=0)
    neighbors = [set() for _ in range(num_scans)]
    for a, b in edges_arr:
        ai = int(a)
        bi = int(b)
        neighbors[ai].add(bi)
        neighbors[bi].add(ai)
    visited = np.zeros(num_scans, dtype=bool)
    component_sizes: list[int] = []
    for start in range(num_scans):
        if visited[start]:
            continue
        stack = [start]
        visited[start] = True
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            for nbr in neighbors[node]:
                if not visited[nbr]:
                    visited[nbr] = True
                    stack.append(nbr)
        component_sizes.append(int(size))
    isolated_scan_ids = [int(uniq_scan_ids[i]) for i in range(num_scans) if len(neighbors[i]) == 0]
    edges_scan_ids = [(int(uniq_scan_ids[a]), int(uniq_scan_ids[b])) for a, b in edges_arr.tolist()]
    return (int(len(component_sizes)), component_sizes, isolated_scan_ids, edges_scan_ids)


def _compute_connectivity_components(scan_i: np.ndarray, scan_j: np.ndarray, num_scans: int) -> list[np.ndarray]:
    if num_scans <= 0:
        return []
    neighbors = [set() for _ in range(int(num_scans))]
    if len(scan_i) > 0:
        edges_arr = np.stack([scan_i, scan_j], axis=1).astype(np.int64, copy=False)
        edges_arr = np.sort(edges_arr, axis=1)
        edges_arr = np.unique(edges_arr, axis=0)
        for a, b in edges_arr:
            ai = int(a)
            bi = int(b)
            if 0 <= ai < num_scans and 0 <= bi < num_scans and ai != bi:
                neighbors[ai].add(bi)
                neighbors[bi].add(ai)
    visited = np.zeros(int(num_scans), dtype=bool)
    components: list[np.ndarray] = []
    for start in range(int(num_scans)):
        if visited[start]:
            continue
        stack = [int(start)]
        visited[start] = True
        comp = []
        while stack:
            node = stack.pop()
            comp.append(int(node))
            for nbr in neighbors[node]:
                if not visited[nbr]:
                    visited[nbr] = True
                    stack.append(int(nbr))
        components.append(np.asarray(comp, dtype=np.int64))
    return components

def _reference_compact_scan_indices(*, uniq_scan_ids: np.ndarray, scan_angles_deg: np.ndarray | None, reference_mode_applied: str | None, reference_direction_deg_applied: float | None, reference_scan_id_applied: int | None, orthogonality_tolerance_deg: float) -> np.ndarray:
    mode = None if reference_mode_applied is None else str(reference_mode_applied).strip().lower()
    if mode == 'direction_zero' and reference_direction_deg_applied is not None and scan_angles_deg is not None:
        diffs = _orientation_difference_deg(np.asarray(scan_angles_deg, dtype=float), float(reference_direction_deg_applied))
        return np.flatnonzero(np.isfinite(diffs) & (diffs <= float(orthogonality_tolerance_deg))).astype(np.int64)
    if mode == 'scan_zero' and reference_scan_id_applied is not None:
        return np.flatnonzero(np.asarray(uniq_scan_ids, dtype=np.int64) == int(reference_scan_id_applied)).astype(np.int64)
    return np.empty(0, dtype=np.int64)

def _apply_componentwise_gauge_fix(coeff_compact: np.ndarray, components: list[np.ndarray], *, reference_mode_applied: str | None, referenced_compact_indices: np.ndarray | None=None, reference_constrain_terms: str='offset_only') -> tuple[np.ndarray, int]:
    coeff = np.asarray(coeff_compact, dtype=float).copy()
    if coeff.size == 0:
        return coeff, 0
    mode = None if reference_mode_applied is None else str(reference_mode_applied).strip().lower()
    constrain_mode = str(reference_constrain_terms).strip().lower()
    if constrain_mode in ('offset_only', 'offset', 'constant_only'):
        term_indices = np.array([0], dtype=np.int64)
    elif constrain_mode in ('all_terms', 'all', 'full'):
        term_indices = np.arange(coeff.shape[1], dtype=np.int64)
    else:
        raise ValueError("reference_constrain_terms must be 'offset_only' or 'all_terms'.")
    referenced = set(int(v) for v in np.asarray(referenced_compact_indices if referenced_compact_indices is not None else np.empty(0, dtype=np.int64), dtype=np.int64).tolist())
    reanchored = 0
    if not components:
        components = [np.arange(coeff.shape[0], dtype=np.int64)]
    for comp in components:
        comp = np.asarray(comp, dtype=np.int64)
        if comp.size == 0:
            continue
        if mode not in (None, '', 'mean_zero') and any(int(v) in referenced for v in comp.tolist()):
            continue
        mean_vals = np.nanmean(coeff[np.ix_(comp, term_indices)], axis=0, keepdims=True)
        coeff[np.ix_(comp, term_indices)] -= mean_vals
        if mode not in (None, '', 'mean_zero'):
            reanchored += 1
    return coeff, reanchored

def estimate_basket_weave_search_radius_arcsec(input_data: GridInput, *, cross_scan_factor: float=0.45, along_scan_factor: float=1.5, fallback_arcsec: float=3.0) -> float:
    """
    Estimate a practical search radius [arcsec] from sampling geometry.

    The estimate is based on two characteristic scales:
      1. median along-scan dump spacing
      2. median nearest-neighbor distance to a point on a different scan

    The returned radius is conservatively chosen as
        min(along_scan_factor * along_step, cross_scan_factor * cross_step)
    whenever both scales are available.
    """
    if cross_scan_factor <= 0:
        raise ValueError('cross_scan_factor must be positive')
    if along_scan_factor <= 0:
        raise ValueError('along_scan_factor must be positive')
    if fallback_arcsec <= 0:
        raise ValueError('fallback_arcsec must be positive')
    scales = estimate_basket_weave_sampling_scales_arcsec(input_data)
    along_step = scales.along_step_arcsec
    cross_step = scales.cross_step_arcsec
    candidates: list[float] = []
    if np.isfinite(along_step) and along_step > 0:
        candidates.append(float(along_scan_factor * along_step))
    if np.isfinite(cross_step) and cross_step > 0:
        candidates.append(float(cross_scan_factor * cross_step))
    if candidates:
        radius = float(np.min(np.asarray(candidates, dtype=float)))
        if np.isfinite(radius) and radius > 0:
            return radius
    return float(fallback_arcsec)

def _as_grid_input(dataset_or_input_data, *, projection: str, ref_coord, frame: str, otf_input_state=None, otf_scan_region=None, otf_scan_png=None, existing_turn_labels: str | None=None, otf_scan_existing_is_turn: str | None=None) -> GridInput:
    if isinstance(dataset_or_input_data, GridInput):
        return dataset_or_input_data
    return create_grid_input(dataset_or_input_data, ref_coord=ref_coord, frame=frame, projection=projection, otf_input_state=otf_input_state, otf_scan_region=otf_scan_region, otf_scan_png=otf_scan_png, existing_turn_labels=existing_turn_labels, otf_scan_existing_is_turn=otf_scan_existing_is_turn)

def _concatenate_optional_numeric(arrays: list[np.ndarray | None], lengths: list[int]) -> np.ndarray | None:
    if all((a is None for a in arrays)):
        return None
    out_list: list[np.ndarray] = []
    for arr, n in zip(arrays, lengths):
        if arr is None:
            out_list.append(np.full(n, np.nan, dtype=float))
        else:
            out_list.append(np.asarray(arr, dtype=float))
    return np.concatenate(out_list, axis=0)

def _concatenate_optional_bool(arrays: list[np.ndarray | None], lengths: list[int], *, default: bool) -> np.ndarray | None:
    if all((a is None for a in arrays)):
        return None
    out_list: list[np.ndarray] = []
    for arr, n in zip(arrays, lengths):
        if arr is None:
            out_list.append(np.full(n, bool(default), dtype=bool))
        else:
            out_list.append(_safe_bool_array(arr, default=default))
    return np.concatenate(out_list, axis=0)


def _summarize_global_scan_ranges(scan_id_chunks: list[np.ndarray], inputs: Sequence[GridInput]) -> list[dict]:
    out: list[dict] = []
    for idx, (sid, inp) in enumerate(zip(scan_id_chunks, inputs)):
        sid = np.asarray(sid, dtype=np.int64)
        valid = sid[sid >= 0]
        uniq = np.unique(valid) if valid.size else np.empty(0, dtype=np.int64)
        source_input_index = getattr(inp, '_source_input_index', idx)
        entry = {
            'input_index': int(source_input_index),
            'global_scan_id_min': None if uniq.size == 0 else int(uniq.min()),
            'global_scan_id_max': None if uniq.size == 0 else int(uniq.max()),
            'num_global_scan_ids': int(uniq.size),
            'num_points_total': int(sid.size),
            'num_points_scan': int(valid.size),
            'source_summary': getattr(inp, '_otf_scan_summary', None),
        }
        out.append(entry)
    return out


def _validate_grid_input_contract(dataset_or_input_data):
    items = dataset_or_input_data if isinstance(dataset_or_input_data, (list, tuple)) else [dataset_or_input_data]
    for idx, item in enumerate(items):
        if not isinstance(item, GridInput):
            continue
        if item.scan_id is None:
            raise ValueError(f'GridInput input {idx} is missing scan_id. When passing GridInput directly to basketweave, the caller must materialize scan_id / is_turnaround / scan_dir consistently in advance.')
        if item.is_turnaround is None:
            raise ValueError(f'GridInput input {idx} is missing is_turnaround. When passing GridInput directly to basketweave, the caller must materialize scan_id / is_turnaround / scan_dir consistently in advance.')


def _validate_reference_scan_id_requested(reference_mode, reference_scan_id: int | None, uniq_scan_ids: np.ndarray):
    mode = None if reference_mode is None else str(reference_mode).strip().lower()
    if mode != 'scan_zero' or reference_scan_id is None:
        return
    uniq = np.asarray(uniq_scan_ids, dtype=np.int64)
    if uniq.size == 0:
        raise ValueError('reference_scan_id was given, but no valid scan_id exists in the basketweave input.')
    ref_id = int(reference_scan_id)
    if not np.any(uniq == ref_id):
        lo = int(uniq.min())
        hi = int(uniq.max())
        raise ValueError(f'reference_scan_id={ref_id} is not present in the basketweave input. Valid global scan_id range is {lo}..{hi}.')


def _merge_grid_inputs(inputs: Sequence[GridInput]) -> GridInput:
    if len(inputs) == 0:
        raise ValueError('At least one GridInput is required.')
    if len(inputs) == 1:
        return inputs[0]
    lengths = [int(len(np.asarray(inp.x))) for inp in inputs]
    spec_arrays = [np.asarray(inp.spec) for inp in inputs]
    first_spec = spec_arrays[0]
    first_ndim = first_spec.ndim
    if first_ndim not in (1, 2):
        raise ValueError(f'spec must be 1D or 2D, got ndim={first_ndim}')
    if first_ndim == 2:
        nchan = first_spec.shape[1]
        for arr in spec_arrays[1:]:
            if arr.ndim != 2 or arr.shape[1] != nchan:
                raise ValueError('All input spectra must have the same channel count when merging GridInput objects.')
    else:
        for arr in spec_arrays[1:]:
            if arr.ndim != 1:
                raise ValueError('Cannot merge 1D and 2D spectra in basket-weave input.')
    x = np.concatenate([np.asarray(inp.x, dtype=float) for inp in inputs], axis=0)
    y = np.concatenate([np.asarray(inp.y, dtype=float) for inp in inputs], axis=0)
    spec = np.concatenate(spec_arrays, axis=0)
    flag_arrays = []
    for inp, n in zip(inputs, lengths):
        if inp.flag is None:
            flag_arrays.append(np.ones(n, dtype=bool))
        else:
            flag_arrays.append(np.asarray(inp.flag) > 0)
    flag = np.concatenate(flag_arrays, axis=0)
    time_arrays = []
    for inp, n in zip(inputs, lengths):
        if inp.time is None:
            time_arrays.append(np.full(n, np.nan, dtype=float))
        else:
            time_arrays.append(np.asarray(inp.time, dtype=float))
    time = np.concatenate(time_arrays, axis=0)
    rms = _concatenate_optional_numeric([inp.rms for inp in inputs], lengths)
    tint = _concatenate_optional_numeric([inp.tint for inp in inputs], lengths)
    tsys = _concatenate_optional_numeric([inp.tsys for inp in inputs], lengths)
    is_turnaround = _concatenate_optional_bool([inp.is_turnaround for inp in inputs], lengths, default=False)
    if any((inp.scan_id is None for inp in inputs)):
        raise ValueError('scan_id must be present for every input when merging basket-weave inputs.')
    scan_id_chunks: list[np.ndarray] = []
    next_scan_id = 0
    for inp in inputs:
        sid = _safe_scan_ids_int64(inp.scan_id)
        out_sid = np.full_like(sid, -1, dtype=np.int64)
        valid = sid >= 0
        if np.any(valid):
            uniq = np.unique(sid[valid])
            mapping = {int(old): int(next_scan_id + i) for i, old in enumerate(uniq)}
            mapped = np.array([mapping[int(v)] for v in sid[valid]], dtype=np.int64)
            out_sid[valid] = mapped
            next_scan_id += len(uniq)
        scan_id_chunks.append(out_sid)
    scan_id = np.concatenate(scan_id_chunks, axis=0)
    global_scan_ranges = _summarize_global_scan_ranges(scan_id_chunks, inputs)
    subscan_id = None
    if any((inp.subscan_id is not None for inp in inputs)):
        subscan_list = []
        for inp, n in zip(inputs, lengths):
            if inp.subscan_id is None:
                subscan_list.append(np.full(n, np.nan, dtype=float))
            else:
                subscan_list.append(np.asarray(inp.subscan_id, dtype=float))
        subscan_id = np.concatenate(subscan_list, axis=0)
    scan_dir = None
    if any((inp.scan_dir is not None for inp in inputs)):
        scan_dir_list = []
        for inp, n in zip(inputs, lengths):
            if inp.scan_dir is None:
                scan_dir_list.append(np.full(n, None, dtype=object))
            else:
                scan_dir_list.append(np.asarray(inp.scan_dir, dtype=object))
        scan_dir = np.concatenate(scan_dir_list, axis=0)
    merged = GridInput(x=x, y=y, spec=spec, flag=flag, time=time, rms=rms, tint=tint, tsys=tsys, scan_id=scan_id, subscan_id=subscan_id, scan_dir=scan_dir, is_turnaround=is_turnaround)
    child_summaries = [getattr(inp, '_otf_scan_summary', None) for inp in inputs]
    merged_summary = {
        'summary_version': 2,
        'kind': 'multi_input',
        'reference_scan_id_semantics': 'global after merge',
        'num_inputs': int(len(inputs)),
        'num_points_total': int(len(scan_id)),
        'num_points_scan_total': int(np.count_nonzero(scan_id >= 0)),
        'num_runs_total': int(len(np.unique(scan_id[scan_id >= 0]))) if np.any(scan_id >= 0) else 0,
        'global_scan_id_ranges': global_scan_ranges,
        'inputs': child_summaries,
    }
    input_states = sorted({str(s.get('input_state')) for s in child_summaries if isinstance(s, dict) and s.get('input_state') is not None})
    if len(input_states) == 1:
        merged_summary['input_state'] = input_states[0]
    setattr(merged, '_otf_scan_summary', merged_summary)
    return merged

def _prepare_input_data(dataset_or_input_data, *, projection: str, ref_coord, frame: str, otf_input_state=None, otf_scan_region=None, otf_scan_png=None, existing_turn_labels: str | None=None, otf_scan_existing_is_turn: str | None=None) -> tuple[GridInput, object | None, bool, list[int] | None]:
    """
    Normalize input into a single GridInput.

    Returns
    -------
    input_data : GridInput
    original_target : object | None
        Original single dataset/GridInput for write-back, or the original list/tuple for merged inputs.
    writeback_to_dataset : bool
        True when a single non-GridInput dataset was supplied.
    segment_lengths : list[int] | None
        Lengths of the merged segments when the input was a list/tuple.
    """
    if isinstance(dataset_or_input_data, (list, tuple)):
        _ensure_not_all_empty_inputs(dataset_or_input_data)
        kept_items, skipped_empty_indices = _filter_nonempty_inputs_for_basketweave(dataset_or_input_data)
        used_items = [item for _, item in kept_items]
        _validate_grid_input_contract(used_items)
        raw_items = [item for item in used_items if not isinstance(item, GridInput)]
        grid_items = [item for item in used_items if isinstance(item, GridInput)]
        ref_coord_use = ref_coord
        if ref_coord_use is None and raw_items:
            if grid_items:
                raise ValueError(
                    'When multi-input basketweave mixes GridInput and scantable objects, an explicit ref_coord is required '
                    'so that all inputs share the same projection reference.'
                )
            lon0, lat0 = _resolve_projection_reference_for_scantables(raw_items, frame=frame)
            ref_coord_use = (float(lon0), float(lat0))
        png_paths_all = _resolve_otf_scan_png_sequence(otf_scan_png, len(dataset_or_input_data))
        grid_inputs = []
        for orig_idx, item in kept_items:
            gi = _as_grid_input(
                item,
                projection=projection,
                ref_coord=ref_coord_use,
                frame=frame,
                otf_input_state=otf_input_state,
                otf_scan_region=otf_scan_region,
                otf_scan_png=png_paths_all[orig_idx],
                existing_turn_labels=existing_turn_labels,
                otf_scan_existing_is_turn=otf_scan_existing_is_turn,
            )
            setattr(gi, '_source_input_index', int(orig_idx))
            grid_inputs.append(gi)
        lengths = [int(len(np.asarray(gi.x))) for gi in grid_inputs]
        merged = _merge_grid_inputs(grid_inputs)
        if not hasattr(merged, '_otf_scan_summary'):
            merged_summary = {
                'summary_version': 2,
                'kind': 'multi_input',
                'reference_scan_id_semantics': 'global after merge',
                'num_inputs': int(len(grid_inputs)),
                'inputs': [getattr(gi, '_otf_scan_summary', None) for gi in grid_inputs],
                'num_inputs_original': int(len(dataset_or_input_data)),
                'num_inputs_used': int(len(grid_inputs)),
                'used_input_indices': [int(v) for v, _ in kept_items],
                'skipped_empty_input_indices': [int(v) for v in skipped_empty_indices],
            }
            setattr(merged, '_otf_scan_summary', merged_summary)
        return (merged, dataset_or_input_data, False, lengths)
    if isinstance(dataset_or_input_data, GridInput):
        _validate_grid_input_contract(dataset_or_input_data)
        return (dataset_or_input_data, None, False, None)
    input_data = _as_grid_input(dataset_or_input_data, projection=projection, ref_coord=ref_coord, frame=frame, otf_input_state=otf_input_state, otf_scan_region=otf_scan_region, otf_scan_png=otf_scan_png, existing_turn_labels=existing_turn_labels, otf_scan_existing_is_turn=otf_scan_existing_is_turn)
    return (input_data, dataset_or_input_data, True, None)

def _sort_scan_sample_indices(indices: np.ndarray, time_values: np.ndarray | None) -> np.ndarray:
    indices = np.asarray(indices, dtype=np.int64)
    if indices.size <= 1:
        return indices
    if time_values is None:
        return indices
    t = np.asarray(time_values[indices], dtype=float)
    if np.count_nonzero(np.isfinite(t)) < 2:
        return indices
    sort_key = np.where(np.isfinite(t), t, np.inf)
    order = np.argsort(sort_key, kind='stable')
    return indices[order]

def _cross2d(ax: np.ndarray, ay: np.ndarray, bx: np.ndarray, by: np.ndarray) -> np.ndarray:
    return np.asarray(ax, dtype=float) * np.asarray(by, dtype=float) - np.asarray(ay, dtype=float) * np.asarray(bx, dtype=float)

@dataclass
class BasketWeaveSolution:
    """Stage3 basket-weave model solution."""
    offset_model: str
    term_names: tuple[str, ...]
    search_radius_arcsec: float
    pair_mode_requested: str
    pair_mode_used: str
    coefficients_full: np.ndarray
    compact_scan_ids: np.ndarray
    compact_coefficients: np.ndarray
    reference_mode_requested: str | None
    reference_mode_applied: str | None
    reference_direction_deg_applied: float | None
    reference_scan_id_applied: int | None
    residuals_before: np.ndarray | None
    residuals_after: np.ndarray | None
    constraint_x_arcsec: np.ndarray | None
    constraint_y_arcsec: np.ndarray | None

@dataclass
class BasketWeaveResult:
    """High-level basket-weave execution result."""
    offsets: np.ndarray
    model_coefficients: np.ndarray
    offset_model_requested: str
    offset_model_used: str
    model_term_names: tuple[str, ...]
    reference_mode_requested: str | None
    reference_mode_applied: str | None
    reference_direction_deg_applied: float | None
    reference_scan_id_applied: int | None
    search_radius_arcsec: float
    pair_mode_requested: str
    pair_mode_used: str
    weighted_solver_used: bool
    num_scans: int
    num_scans_with_direction: int
    num_good_samples: int
    num_constraints: int
    num_pairs_total: int
    num_pairs_cross_scan: int
    num_pairs_used: int
    num_pairs_direction_rejected: int
    num_segment_candidates: int
    num_segment_intersections: int
    direction_filter_applied: bool
    direction_filter_fallback_used: bool
    used_channel_count: int | None
    used_all_channels: bool
    estimated_along_scan_step_arcsec: float | None
    estimated_cross_scan_step_arcsec: float | None
    search_radius_exceeds_cross_scan_spacing: bool
    residual_rms_before: float | None
    residual_rms_after: float | None
    connectivity_num_components: int
    connectivity_component_sizes: list[int]
    connectivity_isolated_scan_ids: list[int]
    connectivity_edges_scan_ids: list[tuple[int, int]]
    constraint_residuals_before: np.ndarray | None = None
    constraint_residuals_after: np.ndarray | None = None
    constraint_x_arcsec: np.ndarray | None = None
    constraint_y_arcsec: np.ndarray | None = None
    search_radius_relevant_for_pair_mode: bool = True
    connectivity_component_reanchored_count: int = 0

@dataclass
class _ConstraintSystemResult:
    scan_i: np.ndarray
    scan_j: np.ndarray
    alpha_i: np.ndarray
    alpha_j: np.ndarray
    rhs_diff: np.ndarray
    weights: np.ndarray
    constraint_x: np.ndarray
    constraint_y: np.ndarray
    n_pairs_total: int
    n_pairs_cross_scan: int
    n_pairs_used: int
    n_pairs_direction_rejected: int
    n_segment_candidates: int
    n_segment_intersections: int
    direction_filter_applied: bool
    direction_filter_fallback_used: bool
    pair_mode_used: str

@dataclass
class _SegmentCollection:
    scan_idx: np.ndarray
    p0x: np.ndarray
    p0y: np.ndarray
    p1x: np.ndarray
    p1y: np.ndarray
    d0: np.ndarray
    d1: np.ndarray
    a0: np.ndarray
    a1: np.ndarray
    seg_angle_deg: np.ndarray
    seg_length_arcsec: np.ndarray

def _resolve_offset_model(offset_model: str) -> tuple[str, tuple[str, ...]]:
    mode = str(offset_model).strip().lower()
    if mode == 'constant':
        return (mode, ('offset',))
    if mode == 'linear':
        return (mode, ('offset', 'slope'))
    raise ValueError("offset_model must be 'constant' or 'linear'.")

def _normalize_reference_direction_deg(reference_direction) -> float | None:
    if reference_direction is None:
        return None
    if isinstance(reference_direction, bytes):
        try:
            reference_direction = reference_direction.decode('utf-8', errors='ignore')
        except Exception:
            return None
    if isinstance(reference_direction, str):
        s = reference_direction.strip().lower()
        aliases_x = {'x', 'scan_x', 'xscan', 'along_x', 'ra', 'az', 'horizontal', 'h'}
        aliases_y = {'y', 'scan_y', 'yscan', 'along_y', 'dec', 'el', 'vertical', 'v'}
        if s in aliases_x:
            return 0.0
        if s in aliases_y:
            return 90.0
        if s == 'auto' or s == '':
            return None
        try:
            val = float(s)
        except ValueError:
            return None
        return float(val % 180.0) if np.isfinite(val) else None
    try:
        val = float(reference_direction)
    except (TypeError, ValueError):
        return None
    return float(val % 180.0) if np.isfinite(val) else None

def _choose_reference_direction_deg(scan_angles_deg: np.ndarray | None) -> float | None:
    if scan_angles_deg is None:
        return None
    ang = np.asarray(scan_angles_deg, dtype=float)
    finite = np.isfinite(ang)
    if not np.any(finite):
        return None
    bins = np.array([0.0, 90.0], dtype=float)
    diffs = np.stack([_orientation_difference_deg(ang[finite], b) for b in bins], axis=1)
    nearest = np.argmin(diffs, axis=1)
    counts = np.bincount(nearest, minlength=2)
    return float(bins[int(np.argmax(counts))])

def _compute_sample_model_coordinate(x: np.ndarray, y: np.ndarray, inv_scan: np.ndarray, *, time_values: np.ndarray | None=None) -> np.ndarray:
    """Per-sample normalized along-scan coordinate in roughly [-1, 1]."""
    coord = np.zeros(len(x), dtype=float)
    for compact_scan_index in np.unique(inv_scan):
        idx = np.flatnonzero(inv_scan == compact_scan_index)
        if idx.size == 0:
            continue
        sidx = _sort_scan_sample_indices(idx, time_values)
        xs = np.asarray(x[sidx], dtype=float)
        ys = np.asarray(y[sidx], dtype=float)
        if sidx.size == 1:
            coord[sidx] = 0.0
            continue
        dx = np.diff(xs)
        dy = np.diff(ys)
        ds = np.hypot(dx, dy)
        ds[~np.isfinite(ds)] = 0.0
        s = np.concatenate([[0.0], np.cumsum(ds)])
        total = float(s[-1]) if s.size else 0.0
        if not np.isfinite(total) or total <= 0:
            local = np.linspace(-1.0, 1.0, sidx.size, dtype=float)
        else:
            center = 0.5 * total
            scale = 0.5 * total
            if not np.isfinite(scale) or scale <= 0:
                local = np.linspace(-1.0, 1.0, sidx.size, dtype=float)
            else:
                local = (s - center) / scale
        coord[sidx] = np.asarray(local, dtype=float)
    return coord

def _compute_full_sample_model_coordinate(input_data: GridInput) -> tuple[np.ndarray, np.ndarray]:
    """Per-sample normalized along-scan coordinate for the full dataset.

    This is the *application* coordinate used by linear basket-weave models.
    It is intentionally computed from the full scan geometry (finite x/y, valid
    scan_id, excluding turnaround) rather than from the solver's filtered sample
    subset. This keeps the along-scan coordinate system consistent between
    solving and applying even when rows contain random flagged/missing spectra
    inside a scan.
    """
    if input_data.scan_id is None:
        raise ValueError('scan_id is required to compute model coordinates.')
    x = np.asarray(input_data.x, dtype=float)
    y = np.asarray(input_data.y, dtype=float)
    scan_ids_all = _safe_scan_ids_int64(input_data.scan_id)
    coord = np.zeros(len(x), dtype=float)
    usable = np.isfinite(x) & np.isfinite(y) & (scan_ids_all >= 0)
    if input_data.is_turnaround is not None:
        usable &= ~_safe_bool_array(input_data.is_turnaround, default=False)
    if np.count_nonzero(usable) == 0:
        return coord, usable
    _, inv = np.unique(scan_ids_all[usable], return_inverse=True)
    time_values = None
    if input_data.time is not None:
        time_values = np.asarray(input_data.time, dtype=float)[usable]
    coord_usable = _compute_sample_model_coordinate(x[usable], y[usable], inv, time_values=time_values)
    coord[usable] = coord_usable
    return coord, usable


def _build_scan_segments(x: np.ndarray, y: np.ndarray, d: np.ndarray, inv_scan: np.ndarray, *, sample_alpha: np.ndarray, time_values: np.ndarray | None=None, min_segment_length_arcsec: float=0.0) -> _SegmentCollection:
    scan_idx_list = []
    p0x_list = []
    p0y_list = []
    p1x_list = []
    p1y_list = []
    d0_list = []
    d1_list = []
    a0_list = []
    a1_list = []
    angle_list = []
    length_list = []
    for compact_scan_index in np.unique(inv_scan):
        idx = np.flatnonzero(inv_scan == compact_scan_index)
        if idx.size < 2:
            continue
        idx = _sort_scan_sample_indices(idx, time_values)
        xs = np.asarray(x[idx], dtype=float)
        ys = np.asarray(y[idx], dtype=float)
        ds = np.asarray(d[idx], dtype=float)
        al = np.asarray(sample_alpha[idx], dtype=float)
        for k in range(idx.size - 1):
            x0 = float(xs[k])
            y0 = float(ys[k])
            x1 = float(xs[k + 1])
            y1 = float(ys[k + 1])
            d0 = float(ds[k])
            d1 = float(ds[k + 1])
            a0 = float(al[k])
            a1 = float(al[k + 1])
            if not (np.isfinite(x0) and np.isfinite(y0) and np.isfinite(x1) and np.isfinite(y1)):
                continue
            if not (np.isfinite(d0) and np.isfinite(d1) and np.isfinite(a0) and np.isfinite(a1)):
                continue
            dx = x1 - x0
            dy = y1 - y0
            seg_len = float(np.hypot(dx, dy))
            if not np.isfinite(seg_len) or seg_len <= float(min_segment_length_arcsec):
                continue
            angle = float(np.degrees(np.arctan2(dy, dx)) % 180.0)
            scan_idx_list.append(int(compact_scan_index))
            p0x_list.append(x0)
            p0y_list.append(y0)
            p1x_list.append(x1)
            p1y_list.append(y1)
            d0_list.append(d0)
            d1_list.append(d1)
            a0_list.append(a0)
            a1_list.append(a1)
            angle_list.append(angle)
            length_list.append(seg_len)
    if not scan_idx_list:
        empty_f = np.empty(0, dtype=float)
        empty_i = np.empty(0, dtype=np.int64)
        return _SegmentCollection(scan_idx=empty_i, p0x=empty_f, p0y=empty_f, p1x=empty_f, p1y=empty_f, d0=empty_f, d1=empty_f, a0=empty_f, a1=empty_f, seg_angle_deg=empty_f, seg_length_arcsec=empty_f)
    return _SegmentCollection(scan_idx=np.asarray(scan_idx_list, dtype=np.int64), p0x=np.asarray(p0x_list, dtype=float), p0y=np.asarray(p0y_list, dtype=float), p1x=np.asarray(p1x_list, dtype=float), p1y=np.asarray(p1y_list, dtype=float), d0=np.asarray(d0_list, dtype=float), d1=np.asarray(d1_list, dtype=float), a0=np.asarray(a0_list, dtype=float), a1=np.asarray(a1_list, dtype=float), seg_angle_deg=np.asarray(angle_list, dtype=float), seg_length_arcsec=np.asarray(length_list, dtype=float))

def _build_constraint_system_from_point_pairs(x: np.ndarray, y: np.ndarray, d: np.ndarray, inv_scan: np.ndarray, *, sample_alpha: np.ndarray, search_radius_arcsec: float, scan_angles_deg: np.ndarray | None, cross_direction_only: bool, orthogonality_tolerance_deg: float, fallback_to_all_cross_scan_pairs: bool) -> _ConstraintSystemResult:
    pair_result = _find_cross_scan_pairs(x, y, inv_scan, search_radius_arcsec=search_radius_arcsec, scan_angles_deg=scan_angles_deg, cross_direction_only=cross_direction_only, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs)
    empty_f = np.empty(0, dtype=float)
    empty_i = np.empty(0, dtype=np.int64)
    if pair_result.n_pairs_used == 0:
        return _ConstraintSystemResult(scan_i=empty_i, scan_j=empty_i, alpha_i=empty_f, alpha_j=empty_f, rhs_diff=empty_f, weights=empty_f, constraint_x=empty_f, constraint_y=empty_f, n_pairs_total=pair_result.n_pairs_total, n_pairs_cross_scan=pair_result.n_pairs_cross_scan, n_pairs_used=0, n_pairs_direction_rejected=pair_result.n_pairs_direction_rejected, n_segment_candidates=0, n_segment_intersections=0, direction_filter_applied=pair_result.direction_filter_applied, direction_filter_fallback_used=pair_result.direction_filter_fallback_used, pair_mode_used='points')
    rhs_diff = np.asarray(d[pair_result.i_idx] - d[pair_result.j_idx], dtype=float)
    weights = np.ones_like(rhs_diff, dtype=float)
    alpha_i = np.asarray(sample_alpha[pair_result.i_idx], dtype=float)
    alpha_j = np.asarray(sample_alpha[pair_result.j_idx], dtype=float)
    cx = 0.5 * (np.asarray(x[pair_result.i_idx], dtype=float) + np.asarray(x[pair_result.j_idx], dtype=float))
    cy = 0.5 * (np.asarray(y[pair_result.i_idx], dtype=float) + np.asarray(y[pair_result.j_idx], dtype=float))
    return _ConstraintSystemResult(scan_i=np.asarray(pair_result.scan_i, dtype=np.int64), scan_j=np.asarray(pair_result.scan_j, dtype=np.int64), alpha_i=alpha_i, alpha_j=alpha_j, rhs_diff=rhs_diff, weights=weights, constraint_x=cx, constraint_y=cy, n_pairs_total=pair_result.n_pairs_total, n_pairs_cross_scan=pair_result.n_pairs_cross_scan, n_pairs_used=pair_result.n_pairs_used, n_pairs_direction_rejected=pair_result.n_pairs_direction_rejected, n_segment_candidates=0, n_segment_intersections=0, direction_filter_applied=pair_result.direction_filter_applied, direction_filter_fallback_used=pair_result.direction_filter_fallback_used, pair_mode_used='points')

def _build_constraint_system_from_segments(x: np.ndarray, y: np.ndarray, d: np.ndarray, inv_scan: np.ndarray, *, sample_alpha: np.ndarray, time_values: np.ndarray | None, search_radius_arcsec: float, scan_angles_deg: np.ndarray | None, cross_direction_only: bool, orthogonality_tolerance_deg: float, fallback_to_all_cross_scan_pairs: bool, min_segment_length_arcsec: float=0.0, segment_parallel_epsilon: float=1e-12, endpoint_weight_floor: float=0.25) -> _ConstraintSystemResult:
    seg = _build_scan_segments(x, y, d, inv_scan, sample_alpha=sample_alpha, time_values=time_values, min_segment_length_arcsec=min_segment_length_arcsec)
    empty_f = np.empty(0, dtype=float)
    empty_i = np.empty(0, dtype=np.int64)
    if len(seg.scan_idx) == 0:
        return _ConstraintSystemResult(scan_i=empty_i, scan_j=empty_i, alpha_i=empty_f, alpha_j=empty_f, rhs_diff=empty_f, weights=empty_f, constraint_x=empty_f, constraint_y=empty_f, n_pairs_total=0, n_pairs_cross_scan=0, n_pairs_used=0, n_pairs_direction_rejected=0, n_segment_candidates=0, n_segment_intersections=0, direction_filter_applied=False, direction_filter_fallback_used=False, pair_mode_used='segments')
    unique_scans = np.unique(seg.scan_idx)
    seg_indices_by_scan = {int(s): np.flatnonzero(seg.scan_idx == s) for s in unique_scans}
    scan_i_list = []
    scan_j_list = []
    alpha_i_list = []
    alpha_j_list = []
    rhs_list = []
    weight_list = []
    cx_list = []
    cy_list = []
    n_pairs_total = 0
    n_pairs_cross_scan = 0
    n_pairs_used = 0
    n_pairs_direction_rejected = 0
    n_segment_candidates = 0
    n_segment_intersections = 0
    direction_filter_applied = False
    direction_filter_fallback_used = False
    pair_scan_i = []
    pair_scan_j = []
    pair_allow = []
    for ia in range(len(unique_scans)):
        for ib in range(ia + 1, len(unique_scans)):
            a = int(unique_scans[ia])
            b = int(unique_scans[ib])
            if cross_direction_only:
                used_mask, _, applied, fallback_used = _select_pair_mask_by_direction(np.array([a], dtype=np.int64), np.array([b], dtype=np.int64), scan_angles_deg, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs)
                allow = bool(used_mask[0])
                direction_filter_applied = direction_filter_applied or bool(applied)
                direction_filter_fallback_used = direction_filter_fallback_used or bool(fallback_used)
            else:
                allow = True
            pair_scan_i.append(a)
            pair_scan_j.append(b)
            pair_allow.append(allow)
    for a, b, allow in zip(pair_scan_i, pair_scan_j, pair_allow):
        idx_a = seg_indices_by_scan[a]
        idx_b = seg_indices_by_scan[b]
        if idx_a.size == 0 or idx_b.size == 0:
            continue
        ax0 = seg.p0x[idx_a][:, None]
        ay0 = seg.p0y[idx_a][:, None]
        ax1 = seg.p1x[idx_a][:, None]
        ay1 = seg.p1y[idx_a][:, None]
        bx0 = seg.p0x[idx_b][None, :]
        by0 = seg.p0y[idx_b][None, :]
        bx1 = seg.p1x[idx_b][None, :]
        by1 = seg.p1y[idx_b][None, :]
        axmin = np.minimum(ax0, ax1)
        axmax = np.maximum(ax0, ax1)
        aymin = np.minimum(ay0, ay1)
        aymax = np.maximum(ay0, ay1)
        bxmin = np.minimum(bx0, bx1)
        bxmax = np.maximum(bx0, bx1)
        bymin = np.minimum(by0, by1)
        bymax = np.maximum(by0, by1)
        overlap = (axmax >= bxmin) & (bxmax >= axmin) & (aymax >= bymin) & (bymax >= aymin)
        cand_ai, cand_bj = np.where(overlap)
        n_cand = int(cand_ai.size)
        n_pairs_total += n_cand
        n_pairs_cross_scan += n_cand
        n_segment_candidates += n_cand
        if n_cand == 0:
            continue
        if not allow:
            n_pairs_direction_rejected += n_cand
            continue
        sidx_a = idx_a[cand_ai]
        sidx_b = idx_b[cand_bj]
        p0x = seg.p0x[sidx_a]
        p0y = seg.p0y[sidx_a]
        p1x = seg.p1x[sidx_a]
        p1y = seg.p1y[sidx_a]
        q0x = seg.p0x[sidx_b]
        q0y = seg.p0y[sidx_b]
        q1x = seg.p1x[sidx_b]
        q1y = seg.p1y[sidx_b]
        rx = p1x - p0x
        ry = p1y - p0y
        sx = q1x - q0x
        sy = q1y - q0y
        qmpx = q0x - p0x
        qmpy = q0y - p0y
        denom = _cross2d(rx, ry, sx, sy)
        nz = np.abs(denom) > float(segment_parallel_epsilon)
        if not np.any(nz):
            continue
        t = np.full_like(denom, np.nan, dtype=float)
        u = np.full_like(denom, np.nan, dtype=float)
        t[nz] = _cross2d(qmpx[nz], qmpy[nz], sx[nz], sy[nz]) / denom[nz]
        u[nz] = _cross2d(qmpx[nz], qmpy[nz], rx[nz], ry[nz]) / denom[nz]
        inside = nz & (t >= -1e-10) & (t <= 1.0 + 1e-10) & (u >= -1e-10) & (u <= 1.0 + 1e-10)
        if not np.any(inside):
            continue
        t = np.clip(t[inside], 0.0, 1.0)
        u = np.clip(u[inside], 0.0, 1.0)
        sidx_a = sidx_a[inside]
        sidx_b = sidx_b[inside]
        da = (1.0 - t) * seg.d0[sidx_a] + t * seg.d1[sidx_a]
        db = (1.0 - u) * seg.d0[sidx_b] + u * seg.d1[sidx_b]
        aa = (1.0 - t) * seg.a0[sidx_a] + t * seg.a1[sidx_a]
        ab = (1.0 - u) * seg.a0[sidx_b] + u * seg.a1[sidx_b]
        rhs = np.asarray(da - db, dtype=float)
        ix = np.asarray(p0x[inside] + t * rx[inside], dtype=float)
        iy = np.asarray(p0y[inside] + t * ry[inside], dtype=float)
        angle_diff = _orientation_difference_deg(seg.seg_angle_deg[sidx_a], seg.seg_angle_deg[sidx_b])
        angle_w = np.sin(np.deg2rad(angle_diff)) ** 2
        end_w = np.sqrt(np.clip(4.0 * t * (1.0 - t), 0.0, 1.0) * np.clip(4.0 * u * (1.0 - u), 0.0, 1.0))
        weights = np.asarray(np.clip(angle_w * (endpoint_weight_floor + (1.0 - endpoint_weight_floor) * end_w), 1e-06, None), dtype=float)
        n_use = int(rhs.size)
        n_pairs_used += n_use
        n_segment_intersections += n_use
        scan_i_list.append(np.full(n_use, a, dtype=np.int64))
        scan_j_list.append(np.full(n_use, b, dtype=np.int64))
        alpha_i_list.append(np.asarray(aa, dtype=float))
        alpha_j_list.append(np.asarray(ab, dtype=float))
        rhs_list.append(rhs)
        weight_list.append(weights)
        cx_list.append(ix)
        cy_list.append(iy)
    if rhs_list:
        return _ConstraintSystemResult(scan_i=np.concatenate(scan_i_list, axis=0), scan_j=np.concatenate(scan_j_list, axis=0), alpha_i=np.concatenate(alpha_i_list, axis=0), alpha_j=np.concatenate(alpha_j_list, axis=0), rhs_diff=np.concatenate(rhs_list, axis=0), weights=np.concatenate(weight_list, axis=0), constraint_x=np.concatenate(cx_list, axis=0), constraint_y=np.concatenate(cy_list, axis=0), n_pairs_total=n_pairs_total, n_pairs_cross_scan=n_pairs_cross_scan, n_pairs_used=n_pairs_used, n_pairs_direction_rejected=n_pairs_direction_rejected, n_segment_candidates=n_segment_candidates, n_segment_intersections=n_segment_intersections, direction_filter_applied=direction_filter_applied, direction_filter_fallback_used=direction_filter_fallback_used, pair_mode_used='segments')
    point_result = _build_constraint_system_from_point_pairs(x, y, d, inv_scan, sample_alpha=sample_alpha, search_radius_arcsec=search_radius_arcsec, scan_angles_deg=scan_angles_deg, cross_direction_only=cross_direction_only, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs)
    point_result.n_pairs_total = builtins.max(point_result.n_pairs_total, n_pairs_total)
    point_result.n_pairs_cross_scan = builtins.max(point_result.n_pairs_cross_scan, n_pairs_cross_scan)
    point_result.n_pairs_direction_rejected = builtins.max(point_result.n_pairs_direction_rejected, n_pairs_direction_rejected)
    point_result.n_segment_candidates = n_segment_candidates
    point_result.n_segment_intersections = 0
    point_result.direction_filter_applied = point_result.direction_filter_applied or direction_filter_applied
    point_result.direction_filter_fallback_used = True
    point_result.pair_mode_used = 'points_fallback'
    return point_result

def _build_constraint_system(x: np.ndarray, y: np.ndarray, d: np.ndarray, inv_scan: np.ndarray, *, sample_alpha: np.ndarray, time_values: np.ndarray | None, search_radius_arcsec: float, scan_angles_deg: np.ndarray | None, pair_mode: str, cross_direction_only: bool, orthogonality_tolerance_deg: float, fallback_to_all_cross_scan_pairs: bool) -> _ConstraintSystemResult:
    mode = str(pair_mode).strip().lower()
    if mode == 'points':
        return _build_constraint_system_from_point_pairs(x, y, d, inv_scan, sample_alpha=sample_alpha, search_radius_arcsec=search_radius_arcsec, scan_angles_deg=scan_angles_deg, cross_direction_only=cross_direction_only, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs)
    if mode == 'segments':
        return _build_constraint_system_from_segments(x, y, d, inv_scan, sample_alpha=sample_alpha, time_values=time_values, search_radius_arcsec=search_radius_arcsec, scan_angles_deg=scan_angles_deg, cross_direction_only=cross_direction_only, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs)
    raise ValueError("pair_mode must be 'points' or 'segments'.")

def _compute_model_prediction(coeff_compact: np.ndarray, scan_i: np.ndarray, scan_j: np.ndarray, alpha_i: np.ndarray, alpha_j: np.ndarray) -> np.ndarray:
    coeff_compact = np.asarray(coeff_compact, dtype=float)
    pred = coeff_compact[scan_i, 0] - coeff_compact[scan_j, 0]
    if coeff_compact.shape[1] >= 2:
        pred = pred + coeff_compact[scan_i, 1] * np.asarray(alpha_i, dtype=float) - coeff_compact[scan_j, 1] * np.asarray(alpha_j, dtype=float)
    return np.asarray(pred, dtype=float)

def _compute_constraint_residual_rms_model(rhs_diff: np.ndarray, scan_i: np.ndarray, scan_j: np.ndarray, alpha_i: np.ndarray, alpha_j: np.ndarray, coeff_compact: np.ndarray, *, weights: np.ndarray | None=None) -> tuple[float | None, float | None, np.ndarray | None, np.ndarray | None]:
    if len(scan_i) == 0:
        return (None, None, None, None)
    before = np.asarray(rhs_diff, dtype=float)
    pred = _compute_model_prediction(coeff_compact, scan_i, scan_j, alpha_i, alpha_j)
    after = before - pred

    def _wrms(values: np.ndarray, w: np.ndarray | None) -> float | None:
        finite = np.isfinite(values)
        if w is not None:
            ww = np.asarray(w, dtype=float)
            finite &= np.isfinite(ww) & (ww > 0)
        if not np.any(finite):
            return None
        if w is None:
            return float(np.sqrt(np.nanmean(values[finite] ** 2)))
        ww = np.asarray(w, dtype=float)[finite]
        vv = values[finite]
        denom = float(np.sum(ww))
        if denom <= 0:
            return None
        return float(np.sqrt(np.sum(ww * vv * vv) / denom))
    return (_wrms(before, weights), _wrms(after, weights), before, after)

def _append_reference_constraints(row_indices: np.ndarray, col_indices: np.ndarray, data_vals: np.ndarray, rhs_w: np.ndarray, *, num_rows_existing: int, coeff_terms: int, uniq_scan_ids: np.ndarray, scan_angles_deg: np.ndarray | None, reference_mode: str | None, reference_direction, reference_scan_id: int | None, reference_constraint_weight: float, orthogonality_tolerance_deg: float, reference_constrain_terms: str='offset_only') -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, str | None, float | None, int | None]:
    mode = None if reference_mode is None else str(reference_mode).strip().lower()
    if mode in (None, '', 'mean_zero'):
        return (row_indices, col_indices, data_vals, rhs_w, mode or 'mean_zero', None, None)
    if reference_constraint_weight <= 0:
        raise ValueError('reference_constraint_weight must be positive when reference constraints are used.')
    selected_compact = None
    applied_direction = None
    applied_scan_id = None
    if mode == 'direction_zero':
        ref_deg = _normalize_reference_direction_deg(reference_direction)
        if ref_deg is None:
            ref_deg = _choose_reference_direction_deg(scan_angles_deg)
        if ref_deg is not None and scan_angles_deg is not None:
            diffs = _orientation_difference_deg(np.asarray(scan_angles_deg, dtype=float), ref_deg)
            selected_compact = np.flatnonzero(np.isfinite(diffs) & (diffs <= float(orthogonality_tolerance_deg)))
            applied_direction = float(ref_deg)
        if selected_compact is None or selected_compact.size == 0:
            return (row_indices, col_indices, data_vals, rhs_w, 'mean_zero', None, None)
    elif mode == 'scan_zero':
        if reference_scan_id is None:
            for candidate in uniq_scan_ids:
                if int(candidate) >= 0:
                    reference_scan_id = int(candidate)
                    break
        if reference_scan_id is None:
            return (row_indices, col_indices, data_vals, rhs_w, 'mean_zero', None, None)
        selected_compact = np.flatnonzero(np.asarray(uniq_scan_ids, dtype=np.int64) == int(reference_scan_id))
        if selected_compact.size == 0:
            return (row_indices, col_indices, data_vals, rhs_w, 'mean_zero', None, None)
        applied_scan_id = int(reference_scan_id)
    else:
        raise ValueError("reference_mode must be one of None, 'mean_zero', 'direction_zero', or 'scan_zero'.")
    constrain_mode = str(reference_constrain_terms).strip().lower()
    if constrain_mode in ('offset_only', 'offset', 'constant_only'):
        term_indices = [0]
    elif constrain_mode in ('all_terms', 'all', 'full'):
        term_indices = list(range(int(coeff_terms)))
    else:
        raise ValueError("reference_constrain_terms must be 'offset_only' or 'all_terms'.")
    start_row = int(num_rows_existing)
    ref_rows = []
    ref_cols = []
    ref_data = []
    rhs_extra = []
    row_counter = start_row
    for s in np.asarray(selected_compact, dtype=np.int64):
        for term in term_indices:
            ref_rows.append(row_counter)
            ref_cols.append(int(s) * int(coeff_terms) + int(term))
            ref_data.append(float(reference_constraint_weight))
            rhs_extra.append(0.0)
            row_counter += 1
    if not ref_rows:
        return (row_indices, col_indices, data_vals, rhs_w, 'mean_zero', None, None)
    row_indices = np.concatenate([row_indices, np.asarray(ref_rows, dtype=np.int64)])
    col_indices = np.concatenate([col_indices, np.asarray(ref_cols, dtype=np.int64)])
    data_vals = np.concatenate([data_vals, np.asarray(ref_data, dtype=float)])
    rhs_w = np.concatenate([rhs_w, np.asarray(rhs_extra, dtype=float)])
    return (row_indices, col_indices, data_vals, rhs_w, mode, applied_direction, applied_scan_id)

def _solve_basket_weave_core(input_data: GridInput, search_radius_arcsec: float | str=3.0, damp: float=0.01, v_axis: np.ndarray | None=None, linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None=None, channel_mask: np.ndarray | None=None, v_windows_kms: list[str] | list[tuple[float, float]] | None=None, cross_direction_only: bool=True, orthogonality_tolerance_deg: float=30.0, fallback_to_all_cross_scan_pairs: bool=True, pair_mode: str='segments', offset_model: str='constant', reference_mode: str | None='mean_zero', reference_direction=None, reference_scan_id: int | None=None, reference_constraint_weight: float=1000.0, reference_constrain_terms: str='offset_only') -> tuple[BasketWeaveSolution, dict]:
    resolved_windows = _resolve_linefree_velocity_windows(linefree_velocity_windows_kms=linefree_velocity_windows_kms, v_windows_kms=v_windows_kms)
    offset_model, term_names = _resolve_offset_model(offset_model)
    if isinstance(search_radius_arcsec, str):
        if search_radius_arcsec.strip().lower() != 'auto':
            raise ValueError("search_radius_arcsec must be a positive float or 'auto'.")
        resolved_search_radius_arcsec = estimate_basket_weave_search_radius_arcsec(input_data)
    else:
        resolved_search_radius_arcsec = float(search_radius_arcsec)
    if resolved_search_radius_arcsec <= 0:
        raise ValueError('search_radius_arcsec must be positive')
    if damp < 0:
        raise ValueError('damp must be non-negative')
    if orthogonality_tolerance_deg < 0 or orthogonality_tolerance_deg > 90:
        raise ValueError('orthogonality_tolerance_deg must be in [0, 90].')
    spec = np.asarray(input_data.spec)
    resolved_channel_mask = _build_channel_mask(spec=spec, v_axis=v_axis, linefree_velocity_windows_kms=resolved_windows, channel_mask=channel_mask)
    x_all = np.asarray(input_data.x, dtype=float)
    y_all = np.asarray(input_data.y, dtype=float)
    if input_data.scan_id is None:
        raise ValueError('scan_id must be provided in GridInput for basket-weaving.')
    scan_ids_all = _safe_scan_ids_int64(input_data.scan_id)
    d_all = np.asarray(_aggregate_spectrum_for_offsets(spec, resolved_channel_mask), dtype=float)
    good = _base_geometry_mask(input_data, x_all=x_all, y_all=y_all, scan_ids_all=scan_ids_all)
    good &= np.isfinite(d_all)
    num_good_samples = int(np.count_nonzero(good))
    sampling_scales = estimate_basket_weave_sampling_scales_arcsec(input_data)
    search_radius_exceeds = _warn_if_search_radius_exceeds_cross_scan_spacing(search_radius_arcsec=resolved_search_radius_arcsec, cross_step_arcsec=sampling_scales.cross_step_arcsec, emit_warning=False)
    diagnostics = dict(offset_model_requested=str(offset_model), offset_model_used=str(offset_model), search_radius_arcsec=float(resolved_search_radius_arcsec), num_good_samples=num_good_samples, num_scans=0, num_scans_with_direction=0, num_constraints=0, num_pairs_total=0, num_pairs_cross_scan=0, num_pairs_used=0, num_pairs_direction_rejected=0, num_segment_candidates=0, num_segment_intersections=0, direction_filter_applied=False, direction_filter_fallback_used=False, pair_mode_requested=str(pair_mode), pair_mode_used=str(pair_mode), estimated_along_scan_step_arcsec=float(sampling_scales.along_step_arcsec) if np.isfinite(sampling_scales.along_step_arcsec) else None, estimated_cross_scan_step_arcsec=float(sampling_scales.cross_step_arcsec) if np.isfinite(sampling_scales.cross_step_arcsec) else None, search_radius_exceeds_cross_scan_spacing=bool(search_radius_exceeds), residual_rms_before=None, residual_rms_after=None, connectivity_num_components=0, connectivity_component_sizes=[], connectivity_isolated_scan_ids=[], connectivity_edges_scan_ids=[], constraint_result=None, resolved_channel_mask=resolved_channel_mask, search_radius_relevant_for_pair_mode=True, connectivity_component_reanchored_count=0)
    if np.count_nonzero(scan_ids_all >= 0) > 0:
        num_scans_full = int(np.max(scan_ids_all[scan_ids_all >= 0])) + 1
    else:
        num_scans_full = 0
    coeff_full = np.zeros((num_scans_full, len(term_names)), dtype=float)
    if num_good_samples < 2:
        sol = BasketWeaveSolution(offset_model=offset_model, term_names=term_names, search_radius_arcsec=float(resolved_search_radius_arcsec), pair_mode_requested=str(pair_mode), pair_mode_used=str(pair_mode), coefficients_full=coeff_full, compact_scan_ids=np.empty(0, dtype=np.int64), compact_coefficients=np.empty((0, len(term_names)), dtype=float), reference_mode_requested=reference_mode, reference_mode_applied='mean_zero', reference_direction_deg_applied=None, reference_scan_id_applied=None, residuals_before=None, residuals_after=None, constraint_x_arcsec=None, constraint_y_arcsec=None)
        return (sol, diagnostics)
    x = x_all[good]
    y = y_all[good]
    d = d_all[good]
    scan_ids = scan_ids_all[good]
    uniq_scan_ids, inv_scan = np.unique(scan_ids, return_inverse=True)
    _validate_reference_scan_id_requested(reference_mode, reference_scan_id, uniq_scan_ids)
    diagnostics['num_scans'] = int(len(uniq_scan_ids))
    scan_dir_good = None
    if input_data.scan_dir is not None:
        scan_dir_good = np.asarray(input_data.scan_dir, dtype=object)[good]
    time_good = None
    if input_data.time is not None:
        time_good = np.asarray(input_data.time, dtype=float)[good]
    direction_info = _infer_compact_scan_directions(x, y, inv_scan, uniq_scan_ids=uniq_scan_ids, scan_dir=scan_dir_good)
    diagnostics['num_scans_with_direction'] = int(np.count_nonzero(np.isfinite(direction_info.angles_deg)))
    full_sample_alpha, full_alpha_usable = _compute_full_sample_model_coordinate(input_data)
    sample_alpha = np.asarray(full_sample_alpha[good], dtype=float)
    constraint_result = _build_constraint_system(x, y, d, inv_scan, sample_alpha=sample_alpha, time_values=time_good, search_radius_arcsec=resolved_search_radius_arcsec, scan_angles_deg=direction_info.angles_deg, pair_mode=pair_mode, cross_direction_only=cross_direction_only, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs)
    diagnostics['constraint_result'] = constraint_result
    diagnostics['pair_mode_used'] = str(constraint_result.pair_mode_used)
    diagnostics['num_constraints'] = int(len(constraint_result.rhs_diff))
    diagnostics['num_pairs_total'] = int(constraint_result.n_pairs_total)
    diagnostics['num_pairs_cross_scan'] = int(constraint_result.n_pairs_cross_scan)
    diagnostics['num_pairs_used'] = int(constraint_result.n_pairs_used)
    diagnostics['num_pairs_direction_rejected'] = int(constraint_result.n_pairs_direction_rejected)
    diagnostics['num_segment_candidates'] = int(constraint_result.n_segment_candidates)
    diagnostics['num_segment_intersections'] = int(constraint_result.n_segment_intersections)
    diagnostics['direction_filter_applied'] = bool(constraint_result.direction_filter_applied)
    diagnostics['direction_filter_fallback_used'] = bool(constraint_result.direction_filter_fallback_used)
    diagnostics['connectivity_num_components'], diagnostics['connectivity_component_sizes'], diagnostics['connectivity_isolated_scan_ids'], diagnostics['connectivity_edges_scan_ids'] = _build_scan_connectivity_diagnostics(constraint_result.scan_i, constraint_result.scan_j, uniq_scan_ids)
    components = _compute_connectivity_components(constraint_result.scan_i, constraint_result.scan_j, len(uniq_scan_ids))
    if int(diagnostics['connectivity_num_components']) > 1:
        logging.warning('Basket-weave constraint graph has %d connected components. Relative zero-points between components are not constrained; applying per-component gauge fixing where needed.', int(diagnostics['connectivity_num_components']))
    if str(constraint_result.pair_mode_used).startswith('segments'):
        diagnostics['search_radius_relevant_for_pair_mode'] = False
        diagnostics['search_radius_exceeds_cross_scan_spacing'] = False
    if constraint_result.n_pairs_used == 0:
        logging.warning('No usable basket-weave constraints found. Returning zero model coefficients.')
        sol = BasketWeaveSolution(offset_model=offset_model, term_names=term_names, search_radius_arcsec=float(resolved_search_radius_arcsec), pair_mode_requested=str(pair_mode), pair_mode_used=str(constraint_result.pair_mode_used), coefficients_full=coeff_full, compact_scan_ids=np.asarray(uniq_scan_ids, dtype=np.int64), compact_coefficients=np.zeros((len(uniq_scan_ids), len(term_names)), dtype=float), reference_mode_requested=reference_mode, reference_mode_applied='mean_zero', reference_direction_deg_applied=None, reference_scan_id_applied=None, residuals_before=None, residuals_after=None, constraint_x_arcsec=None, constraint_y_arcsec=None)
        return (sol, diagnostics)
    rhs_diff = np.asarray(constraint_result.rhs_diff, dtype=float)
    weights = np.asarray(constraint_result.weights, dtype=float)
    w_sqrt = np.sqrt(np.clip(weights, 0.0, None))
    n_eq = len(rhs_diff)
    n_terms = len(term_names)
    num_scans_compact = int(len(uniq_scan_ids))
    if offset_model == 'constant':
        row_indices = np.repeat(np.arange(n_eq, dtype=np.int64), 2)
        col_indices = np.stack([constraint_result.scan_i, constraint_result.scan_j], axis=1).ravel()
        data_vals = np.empty(n_eq * 2, dtype=float)
        data_vals[0::2] = w_sqrt
        data_vals[1::2] = -w_sqrt
    else:
        row_indices = np.repeat(np.arange(n_eq, dtype=np.int64), 4)
        col_indices = np.empty(n_eq * 4, dtype=np.int64)
        col_indices[0::4] = constraint_result.scan_i * 2
        col_indices[1::4] = constraint_result.scan_i * 2 + 1
        col_indices[2::4] = constraint_result.scan_j * 2
        col_indices[3::4] = constraint_result.scan_j * 2 + 1
        data_vals = np.empty(n_eq * 4, dtype=float)
        data_vals[0::4] = w_sqrt
        data_vals[1::4] = w_sqrt * np.asarray(constraint_result.alpha_i, dtype=float)
        data_vals[2::4] = -w_sqrt
        data_vals[3::4] = -w_sqrt * np.asarray(constraint_result.alpha_j, dtype=float)
    rhs_w = rhs_diff * w_sqrt
    effective_reference_mode = reference_mode
    effective_reference_direction = reference_direction
    effective_reference_scan_id = reference_scan_id
    effective_reference_constrain_terms = reference_constrain_terms
    ref_mode_lower = None if reference_mode is None else str(reference_mode).strip().lower()
    ref_terms_lower = str(reference_constrain_terms).strip().lower()
    if offset_model == 'linear' and ref_terms_lower in ('offset_only', 'offset', 'constant_only'):
        logging.warning("offset_model='linear' with reference_constrain_terms=%r leaves slope gauges weakly constrained. Automatically upgrading to 'all_terms'.", reference_constrain_terms)
        effective_reference_constrain_terms = 'all_terms'
    if offset_model == 'linear' and ref_mode_lower in (None, '', 'mean_zero'):
        if np.any(np.isfinite(direction_info.angles_deg)):
            logging.warning("offset_model='linear' with reference_mode='mean_zero' is only weakly anchored. Automatically upgrading to reference_mode='direction_zero' with reference_constrain_terms=%r.", reference_constrain_terms)
            effective_reference_mode = 'direction_zero'
            if effective_reference_direction is None:
                effective_reference_direction = 'auto'
        else:
            logging.warning("offset_model='linear' with reference_mode='mean_zero' is only weakly anchored and no scan directions were inferred. Automatically upgrading to reference_mode='scan_zero' (first valid scan) with reference_constrain_terms=%r.", reference_constrain_terms)
            effective_reference_mode = 'scan_zero'
            if effective_reference_scan_id is None and len(uniq_scan_ids) > 0:
                effective_reference_scan_id = int(uniq_scan_ids[0])
    row_indices, col_indices, data_vals, rhs_w, reference_mode_applied, ref_dir_applied, ref_scan_applied = _append_reference_constraints(row_indices, col_indices, data_vals, rhs_w, num_rows_existing=n_eq, coeff_terms=n_terms, uniq_scan_ids=uniq_scan_ids, scan_angles_deg=direction_info.angles_deg, reference_mode=effective_reference_mode, reference_direction=effective_reference_direction, reference_scan_id=effective_reference_scan_id, reference_constraint_weight=float(reference_constraint_weight), orthogonality_tolerance_deg=float(orthogonality_tolerance_deg), reference_constrain_terms=effective_reference_constrain_terms)
    if offset_model == 'linear' and reference_mode_applied == 'mean_zero' and len(uniq_scan_ids) > 0:
        # Final safety net: if a requested anchor could not be realized, anchor to the first valid scan.
        row_indices, col_indices, data_vals, rhs_w, reference_mode_applied, ref_dir_applied, ref_scan_applied = _append_reference_constraints(row_indices[:n_eq * (2 if n_terms == 1 else 4)], col_indices[:n_eq * (2 if n_terms == 1 else 4)], data_vals[:n_eq * (2 if n_terms == 1 else 4)], rhs_w[:n_eq], num_rows_existing=n_eq, coeff_terms=n_terms, uniq_scan_ids=uniq_scan_ids, scan_angles_deg=direction_info.angles_deg, reference_mode='scan_zero', reference_direction=None, reference_scan_id=int(uniq_scan_ids[0]), reference_constraint_weight=float(reference_constraint_weight), orthogonality_tolerance_deg=float(orthogonality_tolerance_deg), reference_constrain_terms=effective_reference_constrain_terms)
    A = sp.coo_matrix((data_vals, (row_indices, col_indices)), shape=(int(row_indices.max()) + 1 if row_indices.size else 0, num_scans_compact * n_terms))
    res = lsqr(A, rhs_w, damp=damp)
    coeff_compact = np.asarray(res[0], dtype=float).reshape(num_scans_compact, n_terms)
    referenced_compact = _reference_compact_scan_indices(
        uniq_scan_ids=np.asarray(uniq_scan_ids, dtype=np.int64),
        scan_angles_deg=direction_info.angles_deg,
        reference_mode_applied=reference_mode_applied,
        reference_direction_deg_applied=ref_dir_applied,
        reference_scan_id_applied=ref_scan_applied,
        orthogonality_tolerance_deg=float(orthogonality_tolerance_deg),
    )
    coeff_compact, reanchored_count = _apply_componentwise_gauge_fix(
        coeff_compact,
        components,
        reference_mode_applied=reference_mode_applied,
        referenced_compact_indices=referenced_compact,
        reference_constrain_terms=effective_reference_constrain_terms,
    )
    diagnostics['connectivity_component_reanchored_count'] = int(reanchored_count)
    if reference_mode_applied in (None, ''):
        reference_mode_applied = 'mean_zero'
    coeff_full = np.zeros((num_scans_full, n_terms), dtype=float)
    if num_scans_full > 0:
        coeff_full[np.asarray(uniq_scan_ids, dtype=np.int64)] = coeff_compact
    before_rms, after_rms, resid_before, resid_after = _compute_constraint_residual_rms_model(rhs_diff, constraint_result.scan_i, constraint_result.scan_j, constraint_result.alpha_i, constraint_result.alpha_j, coeff_compact, weights=weights)
    diagnostics['residual_rms_before'] = before_rms
    diagnostics['residual_rms_after'] = after_rms
    sol = BasketWeaveSolution(offset_model=offset_model, term_names=term_names, search_radius_arcsec=float(resolved_search_radius_arcsec), pair_mode_requested=str(pair_mode), pair_mode_used=str(constraint_result.pair_mode_used), coefficients_full=coeff_full, compact_scan_ids=np.asarray(uniq_scan_ids, dtype=np.int64), compact_coefficients=np.asarray(coeff_compact, dtype=float), reference_mode_requested=reference_mode, reference_mode_applied=reference_mode_applied, reference_direction_deg_applied=ref_dir_applied, reference_scan_id_applied=ref_scan_applied, residuals_before=resid_before, residuals_after=resid_after, constraint_x_arcsec=np.asarray(constraint_result.constraint_x, dtype=float), constraint_y_arcsec=np.asarray(constraint_result.constraint_y, dtype=float))
    return (sol, diagnostics)

def solve_basket_weave_solution(input_data: GridInput, search_radius_arcsec: float | str=3.0, damp: float=0.01, v_axis: np.ndarray | None=None, linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None=None, channel_mask: np.ndarray | None=None, v_windows_kms: list[str] | list[tuple[float, float]] | None=None, cross_direction_only: bool=True, orthogonality_tolerance_deg: float=30.0, fallback_to_all_cross_scan_pairs: bool=True, pair_mode: str='segments', offset_model: str='constant', reference_mode: str | None='mean_zero', reference_direction=None, reference_scan_id: int | None=None, reference_constraint_weight: float=1000.0, reference_constrain_terms: str='offset_only') -> BasketWeaveSolution:
    sol, _ = _solve_basket_weave_core(input_data, search_radius_arcsec=search_radius_arcsec, damp=damp, v_axis=v_axis, linefree_velocity_windows_kms=linefree_velocity_windows_kms, channel_mask=channel_mask, v_windows_kms=v_windows_kms, cross_direction_only=cross_direction_only, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs, pair_mode=pair_mode, offset_model=offset_model, reference_mode=reference_mode, reference_direction=reference_direction, reference_scan_id=reference_scan_id, reference_constraint_weight=reference_constraint_weight, reference_constrain_terms=reference_constrain_terms)
    return sol

def solve_basket_weave_offsets(input_data: GridInput, search_radius_arcsec: float | str=3.0, damp: float=0.01, v_axis: np.ndarray | None=None, linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None=None, channel_mask: np.ndarray | None=None, v_windows_kms: list[str] | list[tuple[float, float]] | None=None, cross_direction_only: bool=True, orthogonality_tolerance_deg: float=30.0, fallback_to_all_cross_scan_pairs: bool=True, pair_mode: str='segments', offset_model: str='constant') -> np.ndarray:
    """Compatibility wrapper returning only constant offsets.

    For ``offset_model='linear'``, use :func:`solve_basket_weave_solution` or
    :func:`basket_weave_inplace`, because a single offset per scan is no longer
    sufficient to represent the correction.
    """
    mode, _ = _resolve_offset_model(offset_model)
    if mode != 'constant':
        raise ValueError("solve_basket_weave_offsets() only returns per-scan constants. Use solve_basket_weave_solution() or basket_weave_inplace() for offset_model='linear'.")
    sol = solve_basket_weave_solution(input_data, search_radius_arcsec=search_radius_arcsec, damp=damp, v_axis=v_axis, linefree_velocity_windows_kms=linefree_velocity_windows_kms, channel_mask=channel_mask, v_windows_kms=v_windows_kms, cross_direction_only=cross_direction_only, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs, pair_mode=pair_mode, offset_model='constant', reference_mode='mean_zero')
    return np.asarray(sol.coefficients_full[:, 0], dtype=float)

def _compute_correction_vector_from_solution(input_data: GridInput, solution: BasketWeaveSolution) -> np.ndarray:
    scan_ids = _safe_scan_ids_int64(input_data.scan_id)
    n = len(scan_ids)
    corr = np.zeros(n, dtype=float)
    valid = (scan_ids >= 0) & (scan_ids < solution.coefficients_full.shape[0])
    if not np.any(valid):
        return corr
    coeff = np.asarray(solution.coefficients_full, dtype=float)
    corr[valid] = coeff[scan_ids[valid], 0]
    if coeff.shape[1] >= 2:
        alpha_full, alpha_usable = _compute_full_sample_model_coordinate(input_data)
        valid &= np.asarray(alpha_usable, dtype=bool)
        corr[valid] = corr[valid] + coeff[scan_ids[valid], 1] * np.asarray(alpha_full[valid], dtype=float)
    return corr

def apply_basket_weave_correction(input_data: GridInput, offsets_or_solution) -> None:
    """Apply constant or stage3 linear basket-weave correction in-place.

    Both 1D aggregated inputs and 2D spectral inputs are supported.
    """
    if input_data.scan_id is None:
        raise ValueError('scan_id is required in GridInput to apply corrections.')
    spec = np.asarray(input_data.spec)
    if spec.ndim not in (1, 2):
        raise ValueError(f'input_data.spec must be 1D or 2D, got shape={spec.shape}')
    if not np.issubdtype(spec.dtype, np.floating):
        spec = spec.astype(np.float64)
        input_data.spec = spec
    if isinstance(offsets_or_solution, BasketWeaveSolution):
        correction_vector = _compute_correction_vector_from_solution(input_data, offsets_or_solution).astype(spec.dtype, copy=False)
    else:
        offsets = np.asarray(offsets_or_solution, dtype=float)
        scan_ids = _safe_scan_ids_int64(input_data.scan_id)
        valid_offset_mask = (scan_ids >= 0) & (scan_ids < len(offsets))
        correction_vector = np.zeros(len(scan_ids), dtype=spec.dtype)
        correction_vector[valid_offset_mask] = offsets[scan_ids[valid_offset_mask]].astype(spec.dtype, copy=False)
    if spec.ndim == 1:
        input_data.spec -= correction_vector
    else:
        input_data.spec -= correction_vector[:, np.newaxis]

def basket_weave_inplace(dataset_or_input_data, *, projection: str='SFL', ref_coord=None, frame: str='ICRS', search_radius_arcsec: float | str='auto', damp: float=0.01, v_axis: np.ndarray | None=None, linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None=None, channel_mask: np.ndarray | None=None, v_windows_kms: list[str] | list[tuple[float, float]] | None=None, cross_direction_only: bool=True, orthogonality_tolerance_deg: float=30.0, fallback_to_all_cross_scan_pairs: bool=True, pair_mode: str='segments', offset_model: str='constant', reference_mode: str | None='mean_zero', reference_direction=None, reference_scan_id: int | None=None, reference_constraint_weight: float=1000.0, reference_constrain_terms: str='offset_only', otf_input_state=None, otf_scan_region=None, otf_scan_png=None, existing_turn_labels: str | None=None, otf_scan_existing_is_turn: str | None=None) -> BasketWeaveResult:
    _ensure_safe_basketweave_input(dataset_or_input_data, otf_input_state=otf_input_state, otf_scan_region=otf_scan_region)
    effective_v_axis = None if v_axis is None else np.asarray(v_axis, dtype=float)
    resolved_windows = _resolve_linefree_velocity_windows(linefree_velocity_windows_kms=linefree_velocity_windows_kms, v_windows_kms=v_windows_kms)
    if effective_v_axis is None and channel_mask is None and resolved_windows is not None:
        inferred_v_axis = _infer_velocity_axis_kms_for_input(dataset_or_input_data)
        if inferred_v_axis is not None:
            effective_v_axis = np.asarray(inferred_v_axis, dtype=float)
    input_data, original_target, writeback_to_dataset, segment_lengths = _prepare_input_data(dataset_or_input_data, projection=projection, ref_coord=ref_coord, frame=frame, otf_input_state=otf_input_state, otf_scan_region=otf_scan_region, otf_scan_png=otf_scan_png, existing_turn_labels=existing_turn_labels, otf_scan_existing_is_turn=otf_scan_existing_is_turn)
    spec = np.asarray(input_data.spec)
    solution, diagnostics = _solve_basket_weave_core(input_data, search_radius_arcsec=search_radius_arcsec, damp=damp, v_axis=effective_v_axis, linefree_velocity_windows_kms=linefree_velocity_windows_kms, channel_mask=channel_mask, v_windows_kms=v_windows_kms, cross_direction_only=cross_direction_only, orthogonality_tolerance_deg=orthogonality_tolerance_deg, fallback_to_all_cross_scan_pairs=fallback_to_all_cross_scan_pairs, pair_mode=pair_mode, offset_model=offset_model, reference_mode=reference_mode, reference_direction=reference_direction, reference_scan_id=reference_scan_id, reference_constraint_weight=reference_constraint_weight, reference_constrain_terms=reference_constrain_terms)
    apply_basket_weave_correction(input_data, solution)
    if writeback_to_dataset and original_target is not None:
        original_target.data = input_data.spec
    elif segment_lengths is not None and original_target is not None:
        start = 0
        for item, seg_len in zip(original_target, segment_lengths):
            stop = start + int(seg_len)
            chunk = input_data.spec[start:stop]
            if isinstance(item, GridInput):
                item.spec = np.asarray(chunk)
            else:
                item.data = np.asarray(chunk)
            start = stop
    resolved_channel_mask = diagnostics['resolved_channel_mask']
    used_channel_count = None
    used_all_channels = False
    if spec.ndim == 2:
        if resolved_channel_mask is None:
            used_all_channels = True
            used_channel_count = int(spec.shape[1])
        else:
            used_channel_count = int(np.count_nonzero(resolved_channel_mask))
    offsets_full = np.asarray(solution.coefficients_full[:, 0], dtype=float)
    return BasketWeaveResult(offsets=offsets_full, model_coefficients=np.asarray(solution.coefficients_full, dtype=float), offset_model_requested=str(offset_model), offset_model_used=str(solution.offset_model), model_term_names=tuple(solution.term_names), reference_mode_requested=reference_mode, reference_mode_applied=solution.reference_mode_applied, reference_direction_deg_applied=solution.reference_direction_deg_applied, reference_scan_id_applied=solution.reference_scan_id_applied, search_radius_arcsec=float(solution.search_radius_arcsec), pair_mode_requested=str(pair_mode), pair_mode_used=str(solution.pair_mode_used), weighted_solver_used=True, num_scans=int(diagnostics['num_scans']), num_scans_with_direction=int(diagnostics['num_scans_with_direction']), num_good_samples=int(diagnostics['num_good_samples']), num_constraints=int(diagnostics['num_constraints']), num_pairs_total=int(diagnostics['num_pairs_total']), num_pairs_cross_scan=int(diagnostics['num_pairs_cross_scan']), num_pairs_used=int(diagnostics['num_pairs_used']), num_pairs_direction_rejected=int(diagnostics['num_pairs_direction_rejected']), num_segment_candidates=int(diagnostics['num_segment_candidates']), num_segment_intersections=int(diagnostics['num_segment_intersections']), direction_filter_applied=bool(diagnostics['direction_filter_applied']), direction_filter_fallback_used=bool(diagnostics['direction_filter_fallback_used']), used_channel_count=used_channel_count, used_all_channels=bool(used_all_channels), estimated_along_scan_step_arcsec=diagnostics['estimated_along_scan_step_arcsec'], estimated_cross_scan_step_arcsec=diagnostics['estimated_cross_scan_step_arcsec'], search_radius_exceeds_cross_scan_spacing=bool(diagnostics['search_radius_exceeds_cross_scan_spacing']), residual_rms_before=diagnostics['residual_rms_before'], residual_rms_after=diagnostics['residual_rms_after'], connectivity_num_components=int(diagnostics['connectivity_num_components']), connectivity_component_sizes=[int(v) for v in diagnostics['connectivity_component_sizes']], connectivity_isolated_scan_ids=[int(v) for v in diagnostics['connectivity_isolated_scan_ids']], connectivity_edges_scan_ids=[(int(a), int(b)) for a, b in diagnostics['connectivity_edges_scan_ids']], constraint_residuals_before=None if solution.residuals_before is None else np.asarray(solution.residuals_before, dtype=float), constraint_residuals_after=None if solution.residuals_after is None else np.asarray(solution.residuals_after, dtype=float), constraint_x_arcsec=None if solution.constraint_x_arcsec is None else np.asarray(solution.constraint_x_arcsec, dtype=float), constraint_y_arcsec=None if solution.constraint_y_arcsec is None else np.asarray(solution.constraint_y_arcsec, dtype=float), search_radius_relevant_for_pair_mode=bool(diagnostics['search_radius_relevant_for_pair_mode']), connectivity_component_reanchored_count=int(diagnostics['connectivity_component_reanchored_count']))

def plot_basket_weave_constraint_points(result: BasketWeaveResult, *, ax=None):
    if result.constraint_x_arcsec is None or result.constraint_y_arcsec is None:
        raise ValueError('No constraint-point diagnostics are available in result.')
    import matplotlib.pyplot as plt
    created = ax is None
    if created:
        fig, ax = plt.subplots()
    ax.scatter(result.constraint_x_arcsec, result.constraint_y_arcsec, s=6, alpha=0.6)
    ax.set_xlabel('x [arcsec]')
    ax.set_ylabel('y [arcsec]')
    ax.set_title('Basket-weave constraint points')
    return fig if created else ax

def plot_basket_weave_constraint_residual_histogram(result: BasketWeaveResult, *, after: bool=True, ax=None):
    values = result.constraint_residuals_after if after else result.constraint_residuals_before
    if values is None:
        raise ValueError('No constraint residual diagnostics are available in result.')
    import matplotlib.pyplot as plt
    created = ax is None
    if created:
        fig, ax = plt.subplots()
    ax.hist(np.asarray(values, dtype=float), bins=50)
    ax.set_xlabel('Residual')
    ax.set_ylabel('Count')
    ax.set_title('Constraint residuals ({})'.format('after' if after else 'before'))
    return fig if created else ax

def plot_basket_weave_scan_coefficients(result: BasketWeaveResult, *, term: int=0, ax=None):
    coeff = np.asarray(result.model_coefficients, dtype=float)
    if coeff.ndim != 2 or coeff.shape[1] == 0:
        raise ValueError('No model coefficients are available in result.')
    if term < 0 or term >= coeff.shape[1]:
        raise ValueError(f'term must be in [0, {coeff.shape[1] - 1}]')
    import matplotlib.pyplot as plt
    created = ax is None
    if created:
        fig, ax = plt.subplots()
    x = np.arange(coeff.shape[0], dtype=int)
    ax.plot(x, coeff[:, term], marker='o', linestyle='none')
    label = result.model_term_names[term] if term < len(result.model_term_names) else f'term {term}'
    ax.set_xlabel('scan_id')
    ax.set_ylabel(label)
    ax.set_title('Basket-weave scan coefficients')
    return fig if created else ax

# -----------------------------------------------------------------------------
# FFT/PLAIT public wrappers (new standard path)
# -----------------------------------------------------------------------------
from .otf_bundle import OTFBundle  # noqa: E402
from .otf_bundle_io import read_otf_bundle, write_otf_bundle  # noqa: E402
from .otf_family_grid import grid_otf_family  # noqa: E402
from .cube_coadd import coadd_family_cubes  # noqa: E402
from .plait_fft import plait_fft_cubes  # noqa: E402


def basketweave_cubes(
    x_bundle: OTFBundle,
    y_bundle: OTFBundle,
    *,
    linefree_velocity_windows_kms,
    output_fits: str | None = None,
    overwrite: bool = False,
    **kwargs,
) -> OTFBundle:
    """User-facing FFT/PLAIT basket-weave wrapper."""
    out = plait_fft_cubes(
        x_bundle,
        y_bundle,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        **kwargs,
    )
    if output_fits is not None:
        write_otf_bundle(out, output_fits, overwrite=overwrite)
    return out


def basketweave_fits(
    x_fits: str,
    y_fits: str,
    output_fits: str,
    *,
    linefree_velocity_windows_kms,
    overwrite: bool = False,
    **kwargs,
) -> OTFBundle:
    """FITS wrapper for FFT/PLAIT basket-weaving."""
    x_bundle = read_otf_bundle(x_fits)
    y_bundle = read_otf_bundle(y_fits)
    return basketweave_cubes(
        x_bundle,
        y_bundle,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        output_fits=output_fits,
        overwrite=overwrite,
        **kwargs,
    )


def run_otf_plait_pipeline(
    x_inputs,
    y_inputs,
    config,
    *,
    linefree_velocity_windows_kms,
    output_fits: str | None = None,
    x_family_label: str = "X",
    y_family_label: str = "Y",
    x_output_fits: str | None = None,
    y_output_fits: str | None = None,
    overwrite: bool = False,
    grid_kwargs: dict | None = None,
    plait_kwargs: dict | None = None,
    **legacy_kwargs,
) -> OTFBundle:
    """High-level X/Y family gridding + FFT/PLAIT basket-weave pipeline."""
    plait_only_keys = {
        'noise_mode', 'plait_noise_mode', 'pad_frac', 'apodize', 'apodize_alpha', 'support_taper',
        'support_taper_width_pix', 'science_mask_mode', 'fft_workers', 'dtype',
        'diagnostics', 'min_plait_size_pix', 'small_map_policy',
        'quality_gate_mode', 'min_improvement_frac',
    }
    grid_kwargs = {} if grid_kwargs is None else dict(grid_kwargs)
    plait_kwargs = {} if plait_kwargs is None else dict(plait_kwargs)
    for key, value in legacy_kwargs.items():
        if key in plait_only_keys:
            plait_kwargs.setdefault(key, value)
        else:
            grid_kwargs.setdefault(key, value)

    x_bundle = grid_otf_family(
        x_inputs,
        config,
        family_label=x_family_label,
        output_fits=x_output_fits,
        overwrite=overwrite,
        **grid_kwargs,
    )
    y_bundle = grid_otf_family(
        y_inputs,
        config,
        family_label=y_family_label,
        output_fits=y_output_fits,
        overwrite=overwrite,
        **grid_kwargs,
    )
    return basketweave_cubes(
        x_bundle,
        y_bundle,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        output_fits=output_fits,
        overwrite=overwrite,
        **plait_kwargs,
    )
