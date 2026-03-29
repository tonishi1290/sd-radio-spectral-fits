# src/sd_radio_spectral_fits/map/gridder.py
import builtins
import contextlib
import hashlib
import io
import json
import os
import platform
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

from ..regrid_vlsrk import Standardizer
from .core import grid_otf
from .wcs_proj import project_to_plane
from .config import MapConfig, GridInput, normalize_row_flag_mask
from .fits_io import save_map_fits
from ..tempscale import beameff_array, tempscal_array, ta_to_tr
from ..scantable_utils import _df_to_native_endian


class _ConfigView:
    """Read-only config proxy with local runtime overrides."""

    def __init__(self, base, **overrides):
        self._base = base
        self._overrides = {k: v for k, v in overrides.items() if v is not None}

    def __getattr__(self, name):
        if name in self._overrides:
            return self._overrides[name]
        return getattr(self._base, name)


def _coerce_optional_bool(value):
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y", "on"}


def _resolve_weight_mode(config) -> str:
    mode = str(getattr(config, "weight_mode", "uniform")).strip().lower()
    if mode not in {"uniform", "rms"}:
        raise ValueError(f"Unknown weight_mode={mode!r}. Use 'uniform' or 'rms'.")
    return mode


def _effective_config(config, reproducible_mode=None, workers=None, sort_neighbors=None, verbose=None):
    overrides = {}
    rep = _coerce_optional_bool(reproducible_mode)
    if rep is not None:
        overrides["_reproducible_mode"] = rep
    srt = _coerce_optional_bool(sort_neighbors)
    if srt is not None:
        overrides["_sort_neighbors"] = srt
    if workers is not None:
        overrides["_workers"] = int(workers)
    vb = _coerce_optional_bool(verbose)
    if vb is not None:
        overrides["verbose"] = vb
    if not overrides:
        return config
    return _ConfigView(config, **overrides)


def _sha256_file(path: str, chunk_size: int = 1024 * 1024) -> str | None:
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    break
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


def _hash_ndarray(arr: np.ndarray, *, chunk_bytes: int = 1024 * 1024) -> str:
    a = np.ascontiguousarray(np.asarray(arr))
    h = hashlib.sha256()
    h.update(str(a.shape).encode("utf-8"))
    h.update(str(a.dtype).encode("utf-8"))
    view = memoryview(a).cast("B")
    for off in range(0, len(view), chunk_bytes):
        h.update(view[off:off + chunk_bytes])
    return h.hexdigest()


def _nan_stats(arr: np.ndarray) -> dict:
    a = np.asarray(arr, dtype=float)
    finite = np.isfinite(a)
    out = {
        "shape": list(a.shape),
        "dtype": str(a.dtype),
        "finite_count": int(np.count_nonzero(finite)),
        "nan_count": int(np.count_nonzero(~finite)),
    }
    if np.any(finite):
        af = a[finite]
        out.update({
            "min": float(np.nanmin(af)),
            "max": float(np.nanmax(af)),
            "mean": float(np.nanmean(af)),
            "median": float(np.nanmedian(af)),
            "std": float(np.nanstd(af)),
        })
    return out


def _first_last_values(arr: np.ndarray, n: int = 5) -> dict:
    a = np.asarray(arr)
    if a.ndim != 1:
        a = a.ravel()
    return {
        "first": [float(x) if np.isfinite(x) else None for x in a[:n]],
        "last": [float(x) if np.isfinite(x) else None for x in a[-n:]],
    }


def _capture_numpy_runtime() -> str:
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            np.show_runtime()
    except Exception:
        try:
            with contextlib.redirect_stdout(buf):
                np.show_config()
        except Exception:
            pass
    return buf.getvalue()


def _safe_meta_value(meta: dict, key: str):
    try:
        v = meta.get(key)
        if isinstance(v, (np.generic,)):
            return v.item()
        return v
    except Exception:
        return None


def _nanmean_no_warn(arr: np.ndarray, axis=None) -> np.ndarray:
    """nanmean without RuntimeWarning on all-NaN slices."""
    a = np.asarray(arr, dtype=np.float64)
    finite = np.isfinite(a)
    cnt = np.sum(finite, axis=axis)
    summed = np.sum(np.where(finite, a, 0.0), axis=axis, dtype=np.float64)
    out = np.full(np.shape(cnt), np.nan, dtype=np.float64)
    np.divide(summed, cnt, out=out, where=(cnt > 0))
    return out


def _collect_otf_diagnostics(
    *,
    scantable,
    table,
    config,
    coord_sys,
    projection,
    out_scale_norm,
    ref_lon,
    ref_lat,
    full_matrix,
    v_tgt,
    grid_input,
    grid_res,
    output_fits,
) -> dict:
    try:
        import scipy
        scipy_version = scipy.__version__
    except Exception:
        scipy_version = None
    try:
        import astropy
        astropy_version = astropy.__version__
    except Exception:
        astropy_version = None

    meta = getattr(scantable, "meta", {}) or {}
    spec_in = np.asarray(scantable.data)
    agg_in = _nanmean_no_warn(np.asarray(full_matrix, dtype=np.float64), axis=0)
    agg_out = _nanmean_no_warn(np.asarray(grid_res.cube, dtype=np.float64), axis=(0, 1))

    diag = {
        "schema": "otf_diag_v1",
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy_version,
            "astropy": astropy_version,
            "numpy_runtime": _capture_numpy_runtime(),
        },
        "pipeline": {
            "coord_sys": str(coord_sys),
            "projection": str(projection),
            "out_scale": str(out_scale_norm),
            "output_fits": str(output_fits),
        },
        "config": {
            "nx": int(config.nx),
            "ny": int(config.ny),
            "cell_arcsec": float(config.cell_arcsec),
            "beam_fwhm_arcsec": float(config.beam_fwhm_arcsec),
            "kernel": str(config.kernel),
            "kernel_preset": getattr(config, "kernel_preset", None),
            "convsupport": None if getattr(config, "kernel", None) != "sf" else getattr(config, "convsupport", None),
            "kernel_sign": str(getattr(config, "kernel_sign", "auto")),
            "dtype": str(config.dtype),
            "backend": str(config.backend),
            "chunk_ch": int(config.chunk_ch),
            "exclude_turnaround": bool(config.exclude_turnaround),
            "weight_mode": str(_resolve_weight_mode(config)),
            "alpha_rms": float(config.alpha_rms),
            "beta_tint": float(config.beta_tint),
            "dv_kms": None if getattr(config, "dv_kms", None) is None else float(config.dv_kms),
            "vmin_kms": None if getattr(config, "vmin_kms", None) is None else float(config.vmin_kms),
            "vmax_kms": None if getattr(config, "vmax_kms", None) is None else float(config.vmax_kms),
        },
        "scantable_input": {
            "nrow": int(len(table)),
            "table_ncol": int(len(getattr(table, "columns", []))),
            "table_columns": [str(c) for c in getattr(table, "columns", [])],
            "data_shape": list(spec_in.shape),
            "data_dtype": str(spec_in.dtype),
            "meta": {
                "CTYPE1": _safe_meta_value(meta, "CTYPE1"),
                "CUNIT1": _safe_meta_value(meta, "CUNIT1"),
                "CRVAL1": _safe_meta_value(meta, "CRVAL1"),
                "CDELT1": _safe_meta_value(meta, "CDELT1"),
                "CRPIX1": _safe_meta_value(meta, "CRPIX1"),
                "RESTFREQ": _safe_meta_value(meta, "RESTFREQ"),
                "RESTFRQ": _safe_meta_value(meta, "RESTFRQ"),
                "SPECSYS": _safe_meta_value(meta, "SPECSYS"),
            },
        },
        "standardizer": {
            "full_matrix_shape": list(full_matrix.shape),
            "full_matrix_dtype": str(np.asarray(full_matrix).dtype),
            "full_matrix_hash": _hash_ndarray(np.asarray(full_matrix, dtype=np.float32)),
            "agg_spec_hash": _hash_ndarray(np.asarray(agg_in, dtype=np.float64)),
            "agg_spec_stats": _nan_stats(agg_in),
            "v_tgt_nchan": int(len(v_tgt)),
            "v_tgt_hash": _hash_ndarray(np.asarray(v_tgt, dtype=np.float64)),
            "v_tgt_stats": _nan_stats(v_tgt),
            "v_tgt_edges": _first_last_values(v_tgt, n=5),
            "v_tgt_diff_stats": _nan_stats(np.diff(v_tgt) if len(v_tgt) > 1 else np.array([], dtype=float)),
        },
        "projection": {
            "ref_lon_deg": float(ref_lon),
            "ref_lat_deg": float(ref_lat),
            "x_stats": _nan_stats(grid_input.x),
            "y_stats": _nan_stats(grid_input.y),
        },
        "grid_input": {
            "flag_on_count": int(np.count_nonzero(normalize_row_flag_mask(grid_input.flag, ndump=len(grid_input.x), allow_none=True, name='flag'))),
            "time_stats": _nan_stats(grid_input.time),
            "rms_stats": None if grid_input.rms is None else _nan_stats(grid_input.rms),
            "tint_stats": None if grid_input.tint is None else _nan_stats(grid_input.tint),
            "tsys_stats": None if grid_input.tsys is None else _nan_stats(grid_input.tsys),
            "scan_id_stats": None if grid_input.scan_id is None else _nan_stats(pd.to_numeric(np.asarray(grid_input.scan_id), errors="coerce")),
            "is_turnaround_count": None if grid_input.is_turnaround is None else int(np.count_nonzero(np.asarray(grid_input.is_turnaround, dtype=bool))),
            "otf_scan_summary": getattr(grid_input, "_otf_scan_summary", None),
        },
        "grid_output": {
            "cube_shape": list(np.asarray(grid_res.cube).shape),
            "cube_dtype": str(np.asarray(grid_res.cube).dtype),
            "cube_hash": _hash_ndarray(np.asarray(grid_res.cube, dtype=np.float32)),
            "agg_spec_hash": _hash_ndarray(np.asarray(agg_out, dtype=np.float64)),
            "agg_spec_stats": _nan_stats(agg_out),
            "weight_stats": _nan_stats(grid_res.weight_map),
            "hit_stats": _nan_stats(grid_res.hit_map),
            "nsamp_stats": None if getattr(grid_res, "nsamp_map", None) is None else _nan_stats(grid_res.nsamp_map),
            "wsum_stats": None if getattr(grid_res, "wsum_map", None) is None else _nan_stats(grid_res.wsum_map),
            "wabs_stats": None if getattr(grid_res, "wabs_map", None) is None else _nan_stats(grid_res.wabs_map),
            "mask_true": int(np.count_nonzero(np.asarray(grid_res.mask_map, dtype=bool))),
            "time_map_stats": None if getattr(grid_res, "time_map", None) is None else _nan_stats(grid_res.time_map),
            "rms_map_stats": None if getattr(grid_res, "rms_map", None) is None else _nan_stats(grid_res.rms_map),
            "tint_map_stats": None if getattr(grid_res, "tint_map", None) is None else _nan_stats(grid_res.tint_map),
            "tsys_map_stats": None if getattr(grid_res, "tsys_map", None) is None else _nan_stats(grid_res.tsys_map),
            "meta": getattr(grid_res, "meta", None),
        },
        "output_file": {
            "path": str(output_fits),
            "exists": os.path.exists(output_fits),
            "size_bytes": int(os.path.getsize(output_fits)) if os.path.exists(output_fits) else None,
            "sha256": _sha256_file(output_fits) if os.path.exists(output_fits) else None,
        },
    }
    if os.path.exists(output_fits):
        try:
            from astropy.io import fits
            with fits.open(output_fits, memmap=True) as hdul:
                hdr = hdul[0].header
                data = hdul[0].data
                diag["output_fits_header"] = {
                    "NAXIS": hdr.get("NAXIS"),
                    "NAXIS1": hdr.get("NAXIS1"),
                    "NAXIS2": hdr.get("NAXIS2"),
                    "NAXIS3": hdr.get("NAXIS3"),
                    "CTYPE1": hdr.get("CTYPE1"),
                    "CTYPE2": hdr.get("CTYPE2"),
                    "CTYPE3": hdr.get("CTYPE3"),
                    "CUNIT1": hdr.get("CUNIT1"),
                    "CUNIT2": hdr.get("CUNIT2"),
                    "CUNIT3": hdr.get("CUNIT3"),
                    "CRPIX1": hdr.get("CRPIX1"),
                    "CRPIX2": hdr.get("CRPIX2"),
                    "CRPIX3": hdr.get("CRPIX3"),
                    "CRVAL1": hdr.get("CRVAL1"),
                    "CRVAL2": hdr.get("CRVAL2"),
                    "CRVAL3": hdr.get("CRVAL3"),
                    "CDELT1": hdr.get("CDELT1"),
                    "CDELT2": hdr.get("CDELT2"),
                    "CDELT3": hdr.get("CDELT3"),
                    "RESTFRQ": hdr.get("RESTFRQ"),
                    "SPECSYS": hdr.get("SPECSYS"),
                    "data_shape": list(data.shape) if data is not None else None,
                    "data_dtype": str(data.dtype) if data is not None else None,
                }
        except Exception as e:
            diag["output_fits_header_error"] = repr(e)
    return diag


def _format_otf_diag_text(diag: dict) -> str:
    lines = []
    lines.append(f"schema: {diag.get('schema')}")
    rt = diag.get("runtime", {})
    lines.append(f"python: {rt.get('python')}")
    lines.append(f"platform: {rt.get('platform')}")
    lines.append(f"machine: {rt.get('machine')}")
    lines.append(f"numpy: {rt.get('numpy')}  pandas: {rt.get('pandas')}  scipy: {rt.get('scipy')}  astropy: {rt.get('astropy')}")
    pipe = diag.get("pipeline", {})
    lines.append(f"coord_sys: {pipe.get('coord_sys')}  projection: {pipe.get('projection')}  out_scale: {pipe.get('out_scale')}")
    cfg = diag.get("config", {})
    lines.append(f"config: nx={cfg.get('nx')} ny={cfg.get('ny')} cell={cfg.get('cell_arcsec')} beam={cfg.get('beam_fwhm_arcsec')} kernel={cfg.get('kernel')} dtype={cfg.get('dtype')} backend={cfg.get('backend')} chunk_ch={cfg.get('chunk_ch')}")
    std = diag.get("standardizer", {})
    lines.append(f"full_matrix_shape: {std.get('full_matrix_shape')} dtype={std.get('full_matrix_dtype')} hash={std.get('full_matrix_hash')}")
    lines.append(f"v_tgt_nchan: {std.get('v_tgt_nchan')}  v_tgt_hash: {std.get('v_tgt_hash')}")
    edges = std.get("v_tgt_edges", {})
    lines.append(f"v_tgt_first5: {edges.get('first')}")
    lines.append(f"v_tgt_last5 : {edges.get('last')}")
    gout = diag.get("grid_output", {})
    lines.append(f"cube_shape: {gout.get('cube_shape')} dtype={gout.get('cube_dtype')} cube_hash={gout.get('cube_hash')}")
    lines.append(f"output agg_spec_hash: {gout.get('agg_spec_hash')}")
    if gout.get('weight_stats') is not None:
        lines.append(f"weight_stats: {gout.get('weight_stats')}")
    if gout.get('nsamp_stats') is not None:
        lines.append(f"nsamp_stats: {gout.get('nsamp_stats')}")
    if gout.get('wsum_stats') is not None:
        lines.append(f"wsum_stats: {gout.get('wsum_stats')}")
    if gout.get('wabs_stats') is not None:
        lines.append(f"wabs_stats: {gout.get('wabs_stats')}")
    lines.append(f"mask_true: {gout.get('mask_true')}")
    meta = gout.get("meta", {}) or {}
    if meta:
        if meta.get("kernel") is not None:
            lines.append(
                f"kernel_resolved: kernel={meta.get('kernel')} preset={meta.get('kernel_preset')} sign={meta.get('kernel_sign')} "
                f"convsupport={meta.get('convsupport')} "
                f"gwidth={meta.get('gwidth_arcsec', np.nan):.3f} arcsec "
                f"jwidth={meta.get('jwidth_arcsec', np.nan):.3f} arcsec "
                f"support={meta.get('support_radius_arcsec', np.nan):.3f} arcsec "
                f"beam_pix={meta.get('beam_pix', np.nan):.3f} "
                f"cell_over_beam={meta.get('cell_over_beam', np.nan):.4f}"
            )
        if np.isfinite(meta.get("nominal_radial_fwhm_arcsec", np.nan)):
            lines.append(
                f"beam_nominal: radial={meta.get('nominal_radial_fwhm_arcsec'):.3f} arcsec "
                f"bmaj={meta.get('bmaj_nominal_arcsec', np.nan):.3f} arcsec "
                f"bmin={meta.get('bmin_nominal_arcsec', np.nan):.3f} arcsec "
                f"bpa={meta.get('bpa_nominal_deg', np.nan):.2f} deg"
            )
        if np.isfinite(meta.get("empirical_radial_fwhm_arcsec", np.nan)):
            lines.append(
                f"beam_empirical_center: radial={meta.get('empirical_radial_fwhm_arcsec'):.3f} arcsec "
                f"bmaj={meta.get('bmaj_empirical_arcsec', np.nan):.3f} arcsec "
                f"bmin={meta.get('bmin_empirical_arcsec', np.nan):.3f} arcsec "
                f"bpa={meta.get('bpa_empirical_deg', np.nan):.2f} deg"
            )
    of = diag.get("output_file", {})
    lines.append(f"output_file: {of.get('path')}")
    lines.append(f"output_size_bytes: {of.get('size_bytes')}  sha256: {of.get('sha256')}")
    oh = diag.get("output_fits_header", {})
    if oh:
        lines.append(
            f"output_header: NAXIS=({oh.get('NAXIS1')}, {oh.get('NAXIS2')}, {oh.get('NAXIS3')}) "
            f"CTYPE3={oh.get('CTYPE3')} CUNIT3={oh.get('CUNIT3')} "
            f"CRVAL3={oh.get('CRVAL3')} CDELT3={oh.get('CDELT3')} CRPIX3={oh.get('CRPIX3')}"
        )
    return "\n".join(lines) + "\n"


def _print_effective_beam_summary(meta: dict | None, *, indent: str = "   ", print_kernel: bool = True) -> None:
    meta = meta or {}
    if print_kernel and meta:
        print(
            f"{indent}kernel resolved: "
            f"kernel={meta.get('kernel')} "
            f"preset={meta.get('kernel_preset')} "
            f"sign={meta.get('kernel_sign')} "
            f"gwidth={meta.get('gwidth_arcsec', np.nan):.2f}\" "
            f"jwidth={meta.get('jwidth_arcsec', np.nan):.2f}\" "
            f"support={meta.get('support_radius_arcsec', np.nan):.2f}\""
        )
    if np.isfinite(meta.get("nominal_radial_fwhm_arcsec", np.nan)):
        print(
            f"{indent}nominal effective beam: "
            f"radial={meta.get('nominal_radial_fwhm_arcsec'):.2f}\" "
            f"bmaj={meta.get('bmaj_nominal_arcsec', np.nan):.2f}\" "
            f"bmin={meta.get('bmin_nominal_arcsec', np.nan):.2f}\" "
            f"bpa={meta.get('bpa_nominal_deg', np.nan):.2f} deg"
        )
    if np.isfinite(meta.get("empirical_radial_fwhm_arcsec", np.nan)):
        print(
            f"{indent}empirical center beam: "
            f"radial={meta.get('empirical_radial_fwhm_arcsec'):.2f}\" "
            f"bmaj={meta.get('bmaj_empirical_arcsec', np.nan):.2f}\" "
            f"bmin={meta.get('bmin_empirical_arcsec', np.nan):.2f}\" "
            f"bpa={meta.get('bpa_empirical_deg', np.nan):.2f} deg"
        )


def _write_otf_diagnostics(diag: dict, output_fits: str, diagnostics_prefix: str | None = None) -> tuple[str, str]:
    out_path = Path(output_fits)
    if diagnostics_prefix:
        prefix = Path(diagnostics_prefix)
    else:
        prefix = out_path.with_suffix(out_path.suffix + ".otf_diag")
    json_path = str(prefix) + ".json"
    txt_path = str(prefix) + ".txt"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(diag, f, ensure_ascii=False, indent=2)
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write(_format_otf_diag_text(diag))
    return json_path, txt_path


def _time_series_to_mjd_utc(series: pd.Series) -> np.ndarray:
    """Convert a timestamp-like Series to MJD UTC days (float), unit-agnostic."""
    ts = pd.to_datetime(series, errors="coerce", utc=True)
    if len(ts) == 0:
        return np.array([], dtype=float)

    # Avoid relying on the internal datetime64 unit (ns/us/ms/...) because
    # pandas 3.0 may infer lower resolutions such as datetime64[us].
    epoch = pd.Timestamp("1970-01-01", tz="UTC")
    dt_days = (pd.DatetimeIndex(ts) - epoch) / pd.Timedelta(days=1)
    return np.asarray(dt_days, dtype=np.float64) + 40587.0


def _safe_bool_array(values, default: bool = False) -> np.ndarray:
    """Convert mixed boolean-like values into a strict bool ndarray."""
    if values is None:
        return None

    ser = pd.Series(values, copy=False)
    if ser.empty:
        return np.array([], dtype=bool)

    if pd.api.types.is_bool_dtype(ser.dtype):
        return ser.fillna(default).to_numpy(dtype=bool)

    lowered = ser.astype("string").str.strip().str.lower()
    true_vals = {"true", "t", "1", "yes", "y", "on"}
    false_vals = {"false", "f", "0", "no", "n", "off", "", "nan", "none", "null", "<na>"}

    out = np.full(len(ser), bool(default), dtype=bool)
    valid = lowered.notna().to_numpy()
    if np.any(valid):
        arr = lowered.to_numpy(dtype=object)
        m_true = np.array([(x in true_vals) for x in arr], dtype=bool) & valid
        m_false = np.array([(x in false_vals) for x in arr], dtype=bool) & valid
        out[m_true] = True
        out[m_false] = False
    return out


def _safe_scan_id_array(values) -> np.ndarray | None:
    """Convert scan ids to float array while preserving NaN for invalid entries."""
    if values is None:
        return None
    return pd.to_numeric(pd.Series(values, copy=False), errors="coerce").to_numpy(dtype=float)


def _extract_restfreq(scantable, table: pd.DataFrame) -> float:
    """Resolve RESTFREQ/RESTFRQ from meta first, then table columns."""
    meta = getattr(scantable, "meta", {}) or {}
    for key in ("RESTFREQ", "RESTFRQ"):
        val = meta.get(key, None)
        try:
            valf = float(val)
        except (TypeError, ValueError):
            valf = np.nan
        if np.isfinite(valf) and valf > 0:
            return valf

    for key in ("RESTFREQ", "RESTFRQ"):
        if key in table.columns:
            arr = pd.to_numeric(table[key], errors="coerce").to_numpy(dtype=float)
            arr = arr[np.isfinite(arr) & (arr > 0)]
            if arr.size:
                return float(arr[0])

    return 0.0


def _resolve_time_array(table: pd.DataFrame) -> np.ndarray:
    """Resolve per-row time to MJD UTC days from preferred columns."""
    if "MJD" in table.columns:
        mjd = pd.to_numeric(table["MJD"], errors="coerce").to_numpy(float)
        if np.isfinite(mjd).any():
            return mjd

    if "TIMESTAMP" in table.columns:
        mjd = _time_series_to_mjd_utc(table["TIMESTAMP"])
        if np.isfinite(mjd).any():
            return mjd

    date_col = None
    if "DATE-OBS" in table.columns:
        date_col = "DATE-OBS"
    elif "DATEOBS" in table.columns:
        date_col = "DATEOBS"

    if date_col is not None:
        base = _time_series_to_mjd_utc(table[date_col])
        if np.isfinite(base).any():
            if "TIME" in table.columns:
                sec = pd.to_numeric(table["TIME"], errors="coerce").to_numpy(float)
                if len(sec) == len(base):
                    finite_base = base[np.isfinite(base)]
                    finite_sec = sec[np.isfinite(sec)]
                    base_span_sec = (np.nanmax(finite_base) - np.nanmin(finite_base)) * 86400.0 if finite_base.size else np.nan
                    sec_span = np.nanmax(finite_sec) - np.nanmin(finite_sec) if finite_sec.size else np.nan
                    # DATE-OBS が全行ほぼ同一のときだけ TIME をオフセット秒として解釈する。
                    if np.isfinite(base_span_sec) and base_span_sec < 1.0 and np.isfinite(sec_span) and sec_span > 0.0:
                        return base + np.where(np.isfinite(sec), sec / 86400.0, np.nan)
            return base

    # Legacy fallback: TIME stored directly as MJD day.
    if "TIME" in table.columns:
        t = pd.to_numeric(table["TIME"], errors="coerce").to_numpy(float)
        if np.isfinite(t).any():
            finite = t[np.isfinite(t)]
            if finite.size and np.nanmedian(finite) > 1e4:
                return t

    return np.arange(len(table), dtype=float)




def _frame_is_galactic(frame: str) -> bool:
    return str(frame).strip().lower() in {"gal", "galactic"}


def _extract_lon_lat_arrays(table: pd.DataFrame, frame: str) -> tuple[np.ndarray, np.ndarray]:
    is_galactic = _frame_is_galactic(frame)
    lon_col = "GLON" if is_galactic else "RA"
    lat_col = "GLAT" if is_galactic else "DEC"
    if lon_col not in table.columns or lat_col not in table.columns:
        raise ValueError(f"Required coordinate columns not found for frame={frame!r}: {lon_col}, {lat_col}")
    lon_deg = pd.to_numeric(table[lon_col], errors="coerce").to_numpy(float)
    lat_deg = pd.to_numeric(table[lat_col], errors="coerce").to_numpy(float)
    return lon_deg, lat_deg


def _angle_like_to_deg(value) -> float:
    if isinstance(value, np.generic):
        value = value.item()
    if value is None:
        return np.nan
    for attr in ("deg", "degree"):
        if hasattr(value, attr):
            try:
                return float(getattr(value, attr))
            except Exception:
                pass
    if hasattr(value, "to_value"):
        try:
            return float(value.to_value("deg"))
        except Exception:
            pass
    try:
        return float(value)
    except Exception:
        return np.nan


def _resolve_reference_lonlat(ref_coord, frame: str, lon_deg: np.ndarray, lat_deg: np.ndarray, *, ref_lon=None, ref_lat=None) -> tuple[float, float]:
    valid_mask = np.isfinite(lon_deg) & np.isfinite(lat_deg)
    if not np.any(valid_mask):
        raise ValueError("No finite coordinates found.")

    if ref_lon is not None and ref_lat is not None:
        return float(ref_lon), float(ref_lat)

    if ref_coord is not None:
        if isinstance(ref_coord, (tuple, list, np.ndarray)) and len(ref_coord) >= 2:
            lon0 = _angle_like_to_deg(ref_coord[0])
            lat0 = _angle_like_to_deg(ref_coord[1])
            if np.isfinite(lon0) and np.isfinite(lat0):
                return float(lon0), float(lat0)
        try:
            from astropy.coordinates import SkyCoord
        except Exception:
            SkyCoord = None
        if SkyCoord is not None and isinstance(ref_coord, SkyCoord):
            if _frame_is_galactic(frame):
                return float(ref_coord.galactic.l.deg), float(ref_coord.galactic.b.deg)
            return float(ref_coord.icrs.ra.deg), float(ref_coord.icrs.dec.deg)
        if _frame_is_galactic(frame):
            for lon_name, lat_name in (("l", "b"), ("lon", "lat")):
                if hasattr(ref_coord, lon_name) and hasattr(ref_coord, lat_name):
                    lon0 = _angle_like_to_deg(getattr(ref_coord, lon_name))
                    lat0 = _angle_like_to_deg(getattr(ref_coord, lat_name))
                    if np.isfinite(lon0) and np.isfinite(lat0):
                        return float(lon0), float(lat0)
        else:
            for lon_name, lat_name in (("ra", "dec"), ("lon", "lat")):
                if hasattr(ref_coord, lon_name) and hasattr(ref_coord, lat_name):
                    lon0 = _angle_like_to_deg(getattr(ref_coord, lon_name))
                    lat0 = _angle_like_to_deg(getattr(ref_coord, lat_name))
                    if np.isfinite(lon0) and np.isfinite(lat0):
                        return float(lon0), float(lat0)
        raise ValueError("Could not interpret ref_coord as a reference sky position.")

    return float(np.nanmedian(lon_deg[valid_mask])), float(np.nanmedian(lat_deg[valid_mask]))


def _prepare_projection_table_and_coords(
    scantable,
    *,
    frame: str = "ICRS",
    projection: str = "SFL",
    ref_coord=None,
    ref_lon: float | None = None,
    ref_lat: float | None = None,
):
    table = _df_to_native_endian(scantable.table).copy()
    lon_deg, lat_deg = _extract_lon_lat_arrays(table, frame)
    lon0, lat0 = _resolve_reference_lonlat(ref_coord, frame, lon_deg, lat_deg, ref_lon=ref_lon, ref_lat=ref_lat)
    x_arcsec, y_arcsec = project_to_plane(lon_deg, lat_deg, lon0, lat0, projection=projection)
    return table, x_arcsec, y_arcsec, lon0, lat0


def _resolve_projection_reference_for_scantables(
    scantables,
    *,
    frame: str = "ICRS",
    ref_coord=None,
    ref_lon: float | None = None,
    ref_lat: float | None = None,
) -> tuple[float, float]:
    if scantables is None:
        raise ValueError('scantables must not be None.')
    lon_parts = []
    lat_parts = []
    for idx, scantable in enumerate(scantables):
        table = getattr(scantable, 'table', None)
        if table is None:
            raise ValueError(f'scantable at index {idx} does not have a .table attribute.')
        table_df = _df_to_native_endian(table)
        lon_deg, lat_deg = _extract_lon_lat_arrays(table_df, frame)
        valid = np.isfinite(lon_deg) & np.isfinite(lat_deg)
        if np.any(valid):
            lon_parts.append(np.asarray(lon_deg[valid], dtype=float))
            lat_parts.append(np.asarray(lat_deg[valid], dtype=float))
    if not lon_parts:
        raise ValueError('No finite coordinates found across the supplied scantables.')
    lon_all = np.concatenate(lon_parts)
    lat_all = np.concatenate(lat_parts)
    return _resolve_reference_lonlat(ref_coord, frame, lon_all, lat_all, ref_lon=ref_lon, ref_lat=ref_lat)


_VALID_OTF_INPUT_STATES = (
    "with_turnarounds",
    "scan_only",
    "use_existing_labels",
)


def _normalize_otf_input_state(value):
    text = None if value is None else str(value).strip().lower()
    if text is None:
        return None
    if text not in _VALID_OTF_INPUT_STATES:
        raise ValueError(
            "otf_input_state must be one of 'with_turnarounds', 'scan_only', or 'use_existing_labels'."
        )
    return text




def _normalize_existing_turn_labels(value):
    text = None if value is None else str(value).strip().lower()
    if text is None:
        return None
    if text not in {'ignore', 'respect'}:
        raise ValueError("existing_turn_labels must be 'ignore' or 'respect'.")
    return text


def _normalize_legacy_existing_is_turn(value):
    text = None if value is None else str(value).strip().lower()
    if text is None:
        return None
    if text == 'prefer':
        return 'respect'
    if text == 'ignore':
        return 'ignore'
    raise ValueError("otf_scan_existing_is_turn must be 'prefer' or 'ignore'.")


def _resolve_existing_turn_labels_policy(*, otf_state: str, existing_turn_labels=None, otf_scan_existing_is_turn=None):
    new_val = _normalize_existing_turn_labels(existing_turn_labels)
    legacy_given = otf_scan_existing_is_turn is not None
    legacy_val = _normalize_legacy_existing_is_turn(otf_scan_existing_is_turn)
    if otf_state == 'use_existing_labels':
        if new_val is not None or legacy_given:
            raise ValueError(
                "otf_input_state='use_existing_labels' does not accept existing_turn_labels or deprecated otf_scan_existing_is_turn."
            )
        return None
    if new_val is None and not legacy_given:
        return 'ignore'
    if new_val is None:
        warnings.warn(
            "otf_scan_existing_is_turn is deprecated; use existing_turn_labels='respect'|'ignore' instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        return legacy_val
    if legacy_given:
        warnings.warn(
            "otf_scan_existing_is_turn is deprecated; use existing_turn_labels instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        if legacy_val != new_val:
            raise ValueError("existing_turn_labels conflicts with deprecated otf_scan_existing_is_turn.")
    return new_val
def _resolve_otf_processing_policy(*, otf_input_state=None, otf_scan_region=None):
    state = _normalize_otf_input_state(otf_input_state)
    legacy_mode = None
    legacy_explicit = otf_scan_region is not None
    if otf_scan_region not in (None, False):
        legacy_mode = _normalize_otf_scan_region_mode(otf_scan_region)

    if state is None:
        if not legacy_explicit:
            return ("with_turnarounds", "auto")
        if otf_scan_region is False:
            warnings.warn(
                "otf_scan_region is deprecated; use otf_input_state='use_existing_labels' instead.",
                DeprecationWarning,
                stacklevel=3,
            )
            return ("use_existing_labels", None)
        warnings.warn(
            "otf_scan_region is deprecated; use otf_input_state instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        return ("with_turnarounds", legacy_mode or "auto")

    if state == "use_existing_labels":
        if legacy_explicit and otf_scan_region not in (None, False):
            raise ValueError(
                "otf_input_state='use_existing_labels' conflicts with geometry-enabled otf_scan_region."
            )
        if otf_scan_region is False:
            warnings.warn(
                "otf_scan_region is deprecated; use otf_input_state='use_existing_labels' only.",
                DeprecationWarning,
                stacklevel=3,
            )
        return (state, None)

    if otf_scan_region is False:
        raise ValueError(
            f"otf_input_state={state!r} conflicts with otf_scan_region=False. Use otf_input_state='use_existing_labels' instead."
        )
    if legacy_explicit:
        warnings.warn(
            "otf_scan_region is deprecated; use otf_input_state instead.",
            DeprecationWarning,
            stacklevel=3,
        )
    return (state, legacy_mode or "auto")


def _build_existing_labels_summary(*, base_scan_id, existing_is_turn, flag_on, n_total, input_state: str):
    scan_arr = np.asarray(base_scan_id, dtype=float)
    turn_arr = np.asarray(existing_is_turn, dtype=bool)
    flag_arr = np.asarray(flag_on, dtype=bool)
    accepted = np.isfinite(scan_arr) & (~turn_arr) & flag_arr
    accepted_ids = np.unique(scan_arr[accepted]).astype(np.int64) if np.any(accepted) else np.array([], dtype=np.int64)
    return {
        "summary_version": 2,
        "kind": "single_table",
        "reference_scan_id_semantics": "local for a single GridInput; global after merged basketweave inputs",
        "mode": None,
        "input_state": str(input_state),
        "existing_turn_labels": None,
        "existing_is_turn_mode": None,
        "stream_group_columns": [],
        "group_key_columns": ["SCAN"],
        "formal_stream_columns_present": [],
        "missing_formal_stream_columns": [],
        "fallback_stream_columns": [],
        "used_fallback_stream_key": False,
        "stream_group_strategy": "existing_labels",
        "num_groups": int(len(accepted_ids)),
        "num_runs": int(len(accepted_ids)),
        "num_points_total": int(n_total),
        "num_points_flag_on": int(np.count_nonzero(flag_arr)),
        "num_points_existing_turn": int(np.count_nonzero(turn_arr)),
        "num_points_scan": int(np.count_nonzero(accepted)),
        "num_points_turn": int(np.count_nonzero(turn_arr)),
        "num_points_unresolved": 0,
        "num_points_excluded": int(np.count_nonzero(~accepted)),
        "has_existing_is_turn": True,
        "accepted_scan_id_min": None if accepted_ids.size == 0 else int(accepted_ids.min()),
        "accepted_scan_id_max": None if accepted_ids.size == 0 else int(accepted_ids.max()),
        "accepted_scan_id_count": int(accepted_ids.size),
        "geometry_mask_applied": False,
    }


def _normalize_otf_scan_region_mode(value):
    if value in (None, False):
        return None
    if value is True:
        return "auto"
    text = str(value).strip().lower()
    if text in {"0", "false", "f", "no", "n", "off"}:
        return None
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return "auto"
    if text in {"auto", "loose", "normal", "strict"}:
        return text
    raise ValueError(f"Unsupported otf_scan_region mode: {value!r}")


def _resolve_otf_scan_png_path(otf_scan_png, *, default_path: str | None = None) -> str | None:
    if otf_scan_png in (None, False):
        return None
    if otf_scan_png is True:
        return default_path if default_path is not None else "otf_scan_region.png"
    return str(otf_scan_png)


def _assemble_grid_input(
    *,
    scantable,
    table: pd.DataFrame,
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    spec: np.ndarray,
    out_scale_norm: str | None = None,
    beameff_vec: np.ndarray | None = None,
    weight_mode: str = "uniform",
    otf_input_state=None,
    otf_scan_region=None,
    otf_scan_png=None,
    existing_turn_labels: str | None = None,
    otf_scan_existing_is_turn: str | None = None,
):
    obsmode = table["OBSMODE"].astype(str).str.upper().to_numpy() if "OBSMODE" in table.columns else np.full(len(table), "ON")
    flag_on = (obsmode == "ON")
    time_arr = _resolve_time_array(table)
    base_scan_id = _safe_scan_id_array(table["SCAN"]) if "SCAN" in table.columns else None
    existing_is_turn = _safe_bool_array(table["IS_TURN"], default=False) if "IS_TURN" in table.columns else None

    scan_id = base_scan_id
    is_turn = existing_is_turn
    scan_dir = None
    otf_summary = None

    otf_state, otf_mode = _resolve_otf_processing_policy(
        otf_input_state=otf_input_state,
        otf_scan_region=otf_scan_region,
    )
    existing_turn_policy = _resolve_existing_turn_labels_policy(
        otf_state=otf_state,
        existing_turn_labels=existing_turn_labels,
        otf_scan_existing_is_turn=otf_scan_existing_is_turn,
    )
    from .otf_scan_region import identify_otf_scan_regions

    otf_result = identify_otf_scan_regions(
        table,
        x_arcsec=x_arcsec,
        y_arcsec=y_arcsec,
        flag_on=flag_on,
        time_arr=time_arr,
        base_scan_id=base_scan_id,
        existing_is_turn=existing_is_turn,
        mode=otf_mode,
        existing_turn_labels=existing_turn_policy,
        existing_is_turn_mode=otf_scan_existing_is_turn,
        input_state=otf_state,
        png_path=_resolve_otf_scan_png_path(otf_scan_png),
    )
    scan_id = otf_result.effscan
    is_turn = otf_result.is_turn
    scan_dir = otf_result.scan_dir_deg
    otf_summary = otf_result.summary

    return GridInput(
        x=np.asarray(x_arcsec, dtype=float),
        y=np.asarray(y_arcsec, dtype=float),
        spec=np.asarray(spec, dtype=float),
        flag=np.asarray(flag_on, dtype=bool),
        time=np.asarray(time_arr, dtype=float),
        rms=_extract_rms_for_weighting(table, out_scale_norm, beameff_vec, weight_mode=weight_mode),
        tint=_extract_meta_col(table, ("EXPOSURE", "INTTIME", "DUR")),
        tsys=_extract_tsys_scalar(table),
        scan_id=scan_id,
        scan_dir=scan_dir,
        is_turnaround=is_turn,
    ), otf_summary

def run_mapping_pipeline(
    scantable,
    config: MapConfig,
    output_fits: str,
    coord_sys: str = "icrs",
    projection: str = "SFL",
    out_scale: str = "TA*",
    dv_kms: float = None,
    vmin_kms: float = None,
    vmax_kms: float = None,
    ref_coord=None,
    ref_lon: float = None,
    ref_lat: float = None,
    reproducible_mode: bool | None = None,
    workers: int | None = None,
    sort_neighbors: bool | None = None,
    verbose: bool | None = None,
    write_diagnostics: bool = False,
    diagnostics_prefix: str | None = None,
    otf_input_state=None,
    otf_scan_region=None,
    otf_scan_png=None,
    existing_turn_labels: str | None = None,
    otf_scan_existing_is_turn: str | None = None,
):
    """
    Scantable から 3D FITS キューブを生成する汎用統合パイプライン。
    旧 run_otf_pipeline の全機能を包含し、PS観測にも対応可能。
    """
    runtime_config = _effective_config(
        config,
        reproducible_mode=reproducible_mode,
        workers=workers,
        sort_neighbors=sort_neighbors,
        verbose=verbose,
    )

    # 0. SIG_* を除去（coadd後にVRADへ確定している場合は残骸になる）
    t = scantable.table
    if hasattr(t, "columns"):
        meta_ctype1 = str(getattr(scantable, "meta", {}).get("CTYPE1", "")).upper()
        if ("VEL" in meta_ctype1) or ("VRAD" in meta_ctype1):
            sig_cols = [c for c in t.columns if str(c).startswith("SIG_")]
            if sig_cols:
                if getattr(runtime_config, "verbose", False):
                    print(f"[info] drop SIG_* columns: {sig_cols}")
                scantable.table = t.drop(columns=sig_cols)

    # 1. 速度軸標準化 (Standardizer)
    print("1. Regridding velocity axis (Standardizer)...")
    std = Standardizer(scantable, v_corr_col="VFRAME")
    dv_use = dv_kms if dv_kms is not None else getattr(config, "dv_kms", None)
    vmin_use = vmin_kms if vmin_kms is not None else getattr(config, "vmin_kms", None)
    vmax_use = vmax_kms if vmax_kms is not None else getattr(config, "vmax_kms", None)
    full_matrix, v_tgt = std.get_matrix(dv=dv_use, vmin=vmin_use, vmax=vmax_use)

    if getattr(runtime_config, "verbose", False) and (vmin_use is not None or vmax_use is not None):
        print(f"   requested velocity range: {vmin_use} to {vmax_use} km/s")

    if len(v_tgt) > 1:
        dv = float(np.nanmedian(np.diff(v_tgt)))
        dv_nonlin = float(np.nanmax(np.abs(np.diff(v_tgt) - dv)))
        if getattr(runtime_config, "verbose", False):
            print("v_tgt[0:3] =", v_tgt[:3])
            print("v_tgt[-3:] =", v_tgt[-3:])
            print("dv(median) =", dv, " range =", (np.nanmin(v_tgt), np.nanmax(v_tgt)))
            print("nonlinear check (max|d-dmed|) =", dv_nonlin)
        if np.isfinite(dv_nonlin) and np.isfinite(dv) and abs(dv) > 0 and dv_nonlin > builtins.max(1e-6, abs(dv) * 1e-6):
            warnings.warn(
                f"Velocity axis is not perfectly linear (max|Δv-Δv_med|={dv_nonlin:g} km/s). "
                "The FITS WCS will use a linear approximation.",
                RuntimeWarning,
                stacklevel=2,
            )

    # 2. 座標投影 (deg -> arcsec)
    print("2. Projecting coordinates to local plane...")
    table, x_arcsec, y_arcsec, lon0, lat0 = _prepare_projection_table_and_coords(
        scantable,
        frame=coord_sys,
        projection=projection,
        ref_coord=ref_coord,
        ref_lon=ref_lon,
        ref_lat=ref_lat,
    )

    # 3. 温度スケールと BEAMEFF の処理
    print("3. Handling BEAMEFF and Temperature Scale...")
    out_scale_norm = str(out_scale).strip().upper()
    n_rows = len(table)

    beameff_vec = beameff_array(table, scantable.meta, n_rows)
    tempscal_vec = tempscal_array(table, scantable.meta, n_rows)

    rep_beameff = np.nan
    if out_scale_norm == "TR*":
        mask_ta = (tempscal_vec == "TA*")
        if np.any(mask_ta):
            full_matrix[mask_ta] = ta_to_tr(full_matrix[mask_ta], beameff_vec[mask_ta, None])
            print(f" -> Converted {np.count_nonzero(mask_ta)} spectra to TR*.")
        rep_beameff = 1.0  # TR*時は効率1.0として扱う
    else:
        valid_beameff = beameff_vec[np.isfinite(beameff_vec)]
        rep_beameff = float(np.median(valid_beameff)) if len(valid_beameff) > 0 else np.nan

    # 4. GridInput の詳細な組み立て
    weight_mode = _resolve_weight_mode(runtime_config)
    grid_input, otf_scan_summary = _assemble_grid_input(
        scantable=scantable,
        table=table,
        x_arcsec=x_arcsec,
        y_arcsec=y_arcsec,
        spec=full_matrix,
        out_scale_norm=out_scale_norm,
        beameff_vec=beameff_vec,
        weight_mode=weight_mode,
        otf_input_state=otf_input_state,
        otf_scan_region=otf_scan_region,
        otf_scan_png=_resolve_otf_scan_png_path(
            otf_scan_png,
            default_path=str(Path(output_fits).with_suffix(Path(output_fits).suffix + ".scan_region.png")),
        ),
        existing_turn_labels=existing_turn_labels,
        otf_scan_existing_is_turn=otf_scan_existing_is_turn,
    )
    setattr(grid_input, "_otf_scan_summary", otf_scan_summary)
    if getattr(runtime_config, "verbose", False) and otf_scan_summary is not None:
        state_txt = otf_scan_summary.get("input_state")
        if state_txt:
            extra = ''
            if state_txt != 'use_existing_labels':
                etl = otf_scan_summary.get('existing_turn_labels')
                if etl is not None:
                    extra = f", existing_turn_labels={etl}"
            print(f"   otf_input_state={state_txt}{extra}")
        try:
            from .otf_scan_region import format_otf_scan_summary
            print(f"   {format_otf_scan_summary(otf_scan_summary)}")
        except Exception:
            print(f"   otf_scan_region: runs={otf_scan_summary.get('num_runs', '?')} kept={otf_scan_summary.get('num_points_scan', '?')}/{otf_scan_summary.get('num_points_total', '?')}")

    # 5. グリッディング実行
    print("4. Executing Gridding Engine...")
    res = grid_otf(grid_input, runtime_config)

    if not res.meta:
        res.meta = {}
    res.meta["RESTFREQ"] = _extract_restfreq(scantable, table)
    res.meta["SPECSYS"] = "LSRK"

    beam_meta = res.meta or {}
    _print_effective_beam_summary(beam_meta)


    # 6. FITS 書き出し
    print("5. Writing Multi-Extension FITS...")
    save_map_fits(
        res,
        v_tgt,
        coord_sys,
        projection,
        lon0,
        lat0,
        runtime_config,
        output_fits,
        out_scale=out_scale_norm,
        rep_beameff=rep_beameff,
    )

    if write_diagnostics or getattr(runtime_config, "write_diagnostics", False):
        diag = _collect_otf_diagnostics(
            scantable=scantable,
            table=table,
            config=runtime_config,
            coord_sys=coord_sys,
            projection=projection,
            out_scale_norm=out_scale_norm,
            ref_lon=lon0,
            ref_lat=lat0,
            full_matrix=full_matrix,
            v_tgt=v_tgt,
            grid_input=grid_input,
            grid_res=res,
            output_fits=output_fits,
        )
        json_path, txt_path = _write_otf_diagnostics(
            diag,
            output_fits=output_fits,
            diagnostics_prefix=diagnostics_prefix or getattr(runtime_config, "diagnostics_prefix", None),
        )
        print(f"OTF diagnostics written to {json_path} and {txt_path}")

    print(f"Done! Saved to {output_fits}")
    return res


# --- 内部ヘルパー ---

def _extract_rms(table, out_scale, beameff_vec):
    if "BSL_RMS" not in table.columns:
        return np.ones(len(table))
    rms = pd.to_numeric(table["BSL_RMS"], errors="coerce").to_numpy(float, copy=True)
    rms[~np.isfinite(rms) | (rms <= 0)] = 1.0
    if out_scale == "TR*":
        rms /= np.where(beameff_vec > 0, beameff_vec, 1.0)
    return rms


def _extract_rms_for_weighting(table, out_scale, beameff_vec, *, weight_mode: str):
    mode = str(weight_mode).strip().lower()
    if mode == "uniform":
        return None
    if mode != "rms":
        raise ValueError(f"Unknown weight_mode={mode!r}. Use 'uniform' or 'rms'.")
    if "BSL_RMS" not in table.columns:
        raise ValueError(
            "weight_mode='rms' requires a BSL_RMS column for every dump, but the input table has no BSL_RMS column. "
            "Run baseline statistics first, or set weight_mode='uniform'."
        )
    rms = pd.to_numeric(table["BSL_RMS"], errors="coerce").to_numpy(float, copy=True)
    bad = (~np.isfinite(rms)) | (rms <= 0)
    if np.any(bad):
        n_bad = int(np.count_nonzero(bad))
        raise ValueError(
            "weight_mode='rms' requires valid positive BSL_RMS for every dump. "
            f"Found {n_bad}/{len(rms)} invalid rows in BSL_RMS. "
            "This often means some scantables lack BSL_RMS or contain NaN/non-positive values. "
            "Provide BSL_RMS for all inputs, or set weight_mode='uniform'."
        )
    if str(out_scale).strip().upper() == "TR*":
        if beameff_vec is None:
            raise ValueError(
                "weight_mode='rms' with out_scale='TR*' requires BEAMEFF information to scale BSL_RMS consistently."
            )
        rms /= np.where(beameff_vec > 0, beameff_vec, 1.0)
    return rms


def _extract_tsys_scalar(table):
    if "TSYS" not in table.columns:
        return None

    out = np.full(len(table), np.nan, dtype=float)
    for i, x in enumerate(table["TSYS"]):
        if isinstance(x, (np.ndarray, list, tuple)):
            arr = np.asarray(x, dtype=float)
            if arr.size and np.isfinite(arr).any():
                out[i] = np.nanmean(arr)
            else:
                out[i] = np.nan
            continue
        try:
            out[i] = float(x)
        except (TypeError, ValueError):
            out[i] = np.nan
    return out


def _extract_meta_col(table, names):
    for n in names:
        if n in table.columns:
            return pd.to_numeric(table[n], errors="coerce").to_numpy(float)
    return None


def create_grid_input(
    scantable,
    ref_coord=None,
    ref_lon: float | None = None,
    ref_lat: float | None = None,
    frame="ICRS",
    projection="SFL",
    weight_mode: str = "uniform",
    otf_input_state=None,
    otf_scan_region=None,
    otf_scan_png=None,
    existing_turn_labels: str | None = None,
    otf_scan_existing_is_turn: str | None = None,
):
    """
    Scantable から GridInput を安全に生成する。
    Basket-weave と最終 gridding で同じ座標投影・同じ OTF 領域判定を使う。
    """
    table, x_arcsec, y_arcsec, _, _ = _prepare_projection_table_and_coords(
        scantable,
        frame=frame,
        projection=projection,
        ref_coord=ref_coord,
        ref_lon=ref_lon,
        ref_lat=ref_lat,
    )
    grid_input, otf_scan_summary = _assemble_grid_input(
        scantable=scantable,
        table=table,
        x_arcsec=x_arcsec,
        y_arcsec=y_arcsec,
        spec=np.asarray(scantable.data, dtype=float),
        out_scale_norm=None,
        beameff_vec=None,
        weight_mode=weight_mode,
        otf_input_state=otf_input_state,
        otf_scan_region=otf_scan_region,
        otf_scan_png=_resolve_otf_scan_png_path(otf_scan_png),
        existing_turn_labels=existing_turn_labels,
        otf_scan_existing_is_turn=otf_scan_existing_is_turn,
    )
    setattr(grid_input, "_otf_scan_summary", otf_scan_summary)
    return grid_input
