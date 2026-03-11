# src/sd_radio_spectral_fits/map/gridder.py
import builtins
import warnings
import numpy as np
import pandas as pd

from ..regrid_vlsrk import Standardizer
from .core import grid_otf
from .wcs_proj import project_to_plane
from .config import MapConfig, GridInput
from .fits_io import save_map_fits
from ..tempscale import beameff_array, tempscal_array, ta_to_tr
from ..scantable_utils import calc_mapping_offsets, _df_to_native_endian


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


def _effective_config(config, reproducible_mode=None, workers=None, sort_neighbors=None):
    overrides = {}
    rep = _coerce_optional_bool(reproducible_mode)
    if rep is not None:
        overrides["_reproducible_mode"] = rep
    srt = _coerce_optional_bool(sort_neighbors)
    if srt is not None:
        overrides["_sort_neighbors"] = srt
    if workers is not None:
        overrides["_workers"] = int(workers)
    if not overrides:
        return config
    return _ConfigView(config, **overrides)


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


def run_mapping_pipeline(
    scantable,
    config: MapConfig,
    output_fits: str,
    coord_sys: str = "icrs",
    projection: str = "SFL",
    out_scale: str = "TA*",
    dv_kms: float = None,
    ref_lon: float = None,
    ref_lat: float = None,
    reproducible_mode: bool | None = None,
    workers: int | None = None,
    sort_neighbors: bool | None = None,
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
    full_matrix, v_tgt = std.get_matrix(dv=dv_use)

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
    table = _df_to_native_endian(scantable.table).copy()
    is_galactic = coord_sys.lower() in ("galactic", "gal")
    lon_deg = pd.to_numeric(table["GLON" if is_galactic else "RA"], errors="coerce").to_numpy(float)
    lat_deg = pd.to_numeric(table["GLAT" if is_galactic else "DEC"], errors="coerce").to_numpy(float)

    valid_mask = np.isfinite(lon_deg) & np.isfinite(lat_deg)
    if not np.any(valid_mask):
        raise ValueError("No finite coordinates found.")

    lon0 = float(ref_lon) if ref_lon is not None else float(np.nanmedian(lon_deg[valid_mask]))
    lat0 = float(ref_lat) if ref_lat is not None else float(np.nanmedian(lat_deg[valid_mask]))

    x_arcsec, y_arcsec = project_to_plane(lon_deg, lat_deg, lon0, lat0, projection=projection)

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
    grid_input = GridInput(
        x=x_arcsec,
        y=y_arcsec,
        spec=full_matrix,
        flag=(table["OBSMODE"].astype(str).str.upper() == "ON").to_numpy(dtype=bool) if "OBSMODE" in table.columns else np.ones(len(table), dtype=bool),
        time=_resolve_time_array(table),
        rms=_extract_rms(table, out_scale_norm, beameff_vec),
        tint=_extract_meta_col(table, ("EXPOSURE", "INTTIME", "DUR")),
        tsys=_extract_tsys_scalar(table),
        scan_id=_safe_scan_id_array(table["SCAN"]) if "SCAN" in table.columns else None,
        is_turnaround=_safe_bool_array(table["IS_TURN"], default=False) if "IS_TURN" in table.columns else None,
    )

    # 5. グリッディング実行
    print("4. Executing Gridding Engine...")
    res = grid_otf(grid_input, runtime_config)

    if not res.meta:
        res.meta = {}
    res.meta["RESTFREQ"] = _extract_restfreq(scantable, table)
    res.meta["SPECSYS"] = "LSRK"

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


def create_grid_input(scantable, ref_coord=None, frame="ICRS", projection="SFL"):
    """
    ScantableからGridInputを安全に生成する。
    ※ 将来的に src/map/gridder.py に移設することを推奨します。
    """
    table = _df_to_native_endian(scantable.table).copy()

    offsets_df = calc_mapping_offsets(
        scantable,
        ref_coord=ref_coord,
        frame=frame,
        projection=projection,
        unit="arcsec",
        verbose=False,
    )

    x_arcsec = offsets_df["OFS_LON"].to_numpy()
    y_arcsec = offsets_df["OFS_LAT"].to_numpy()

    obsmode = table["OBSMODE"].astype(str).str.upper().to_numpy() if "OBSMODE" in table.columns else np.full(len(table), "ON")
    time_arr = _resolve_time_array(table)
    scan_id = _safe_scan_id_array(table["SCAN"]) if "SCAN" in table.columns else None
    is_turn = _safe_bool_array(table["IS_TURN"], default=False) if "IS_TURN" in table.columns else None

    return GridInput(
        x=x_arcsec,
        y=y_arcsec,
        spec=np.asarray(scantable.data, dtype=float),
        flag=(obsmode == "ON"),
        time=time_arr,
        scan_id=scan_id,
        is_turnaround=is_turn,
    )
