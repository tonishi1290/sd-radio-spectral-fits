# -*- coding: utf-8 -*-
"""
sd_radio_spectral_fits.map.cube_baseline.cli

CLI-like pipeline function for iterative baseline subtraction.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from astropy.io import fits
import astropy.units as u

try:
    import psutil  # type: ignore
except Exception:  # pragma: no cover
    psutil = None  # type: ignore

from sd_radio_spectral_fits.map_3d.baseline_subtraction import (
    LineFreeConfig,
    RippleConfig,
    build_baseline_diagnostic_payload,
    write_baseline_diagnostic_files,
)
from .session import BaselineSession
from .orchestrator import run_one_iteration

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def _profile_rss_gb() -> Optional[float]:
    if psutil is None:
        return None
    try:
        return float(psutil.Process(os.getpid()).memory_info().rss) / 1.0e9
    except Exception:
        return None


def _profile_record(store: List[Dict[str, Any]], stage: str, t_prev: float, t0: float, **extra: Any) -> float:
    now = time.perf_counter()
    rec: Dict[str, Any] = {"stage": str(stage), "dt_stage_sec": float(now - t_prev), "dt_total_sec": float(now - t0)}
    rss = _profile_rss_gb()
    if rss is not None:
        rec["rss_gb"] = float(rss)
    rec.update(extra)
    store.append(rec)
    return now


def _write_profile_sidecar(profile_path_prefix: str, records: List[Dict[str, Any]]) -> Tuple[str, str]:
    json_path = profile_path_prefix + ".baseline_profile.json"
    txt_path = profile_path_prefix + ".baseline_profile.txt"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump({"records": records}, f, ensure_ascii=False, indent=2)
    lines = ["# cli baseline profile"]
    for rec in records:
        rss = rec.get("rss_gb")
        rss_txt = f"  rss_gb={rss:.6f}" if isinstance(rss, (int, float)) else ""
        lines.append(f"{rec['stage']}: dt_stage_sec={rec['dt_stage_sec']:.6f}  dt_total_sec={rec['dt_total_sec']:.6f}{rss_txt}")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return json_path, txt_path


def _get_cube_from_hdul(
    hdul: fits.HDUList,
    cube_ext: Optional[Union[int, str]] = None,
) -> Tuple[np.ndarray, fits.Header, int]:
    image_types = (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)
    excluded_names = {
        "MASK3D", "BASESUP3D", "LINECAND3D",
        "RMS", "BASE_RMS", "WEIGHT", "HIT", "MASK", "TSYS", "TINT", "TIME",
        "MOMENT0", "MOM0_BASESUP", "MOM0_LINECAND",
    }
    excluded_btypes = {
        "SIGNALMASK", "BASELINESUPPORTMASK", "LINECANDIDATEMASK",
        "MOMENT0", "MOMENT0BASELINESUPPORT", "MOMENT0LINECANDIDATE",
        "WEIGHT", "HITCOUNT", "VALIDMASK", "SYSTEMTEMP", "INTEGRATIONTIME", "OBSERVATIONTIME",
        "BASELINERESIDRMS", "BASELINERMS",
    }

    def _is_3d_image(hdu: object) -> bool:
        return isinstance(hdu, image_types) and getattr(hdu, "data", None) is not None and np.ndim(hdu.data) == 3

    def _looks_like_analysis_hdu(hdu: object) -> bool:
        name = str(getattr(hdu, "name", "") or "").upper()
        btype = str(getattr(getattr(hdu, "header", None), "get", lambda *_: "")("BTYPE", "") or "").upper()
        return (name in excluded_names) or (btype in excluded_btypes)

    def _is_cli_usable(hdu: object) -> bool:
        hdr = getattr(hdu, "header", None)
        if hdr is None:
            return False
        return _is_cli_usable_spectral_hdu(hdr)

    if cube_ext is None:
        if _is_3d_image(hdul[0]) and (not _looks_like_analysis_hdu(hdul[0])) and _is_cli_usable(hdul[0]):
            return np.asarray(hdul[0].data), hdul[0].header.copy(), 0
        for i, hdu in enumerate(hdul):
            if _is_3d_image(hdu) and (not _looks_like_analysis_hdu(hdu)) and _is_cli_usable(hdu):
                return np.asarray(hdu.data), hdu.header.copy(), i
        raise ValueError(
            "No suitable CLI-usable 3D spectral cube found in FITS. Need a velocity-like axis or a frequency-like axis with RESTFRQ/RESTFREQ."
        )

    hdu = hdul[cube_ext]
    if not _is_3d_image(hdu):
        shape = getattr(getattr(hdu, "data", None), "shape", None)
        raise ValueError(f"Selected cube_ext={cube_ext!r} is not a 3D image HDU (shape={shape}).")
    if _looks_like_analysis_hdu(hdu):
        raise ValueError(f"Selected cube_ext={cube_ext!r} points to an analysis/product HDU, not an input data cube.")
    if not _is_cli_usable(hdu):
        hdr = getattr(hdu, "header", None)
        kind = _cli_axis_kind(hdr) if hdr is not None else None
        if kind == "wavelength":
            raise ValueError(
                f"Selected cube_ext={cube_ext!r} uses a wavelength-like spectral axis (WAVE/AWAV); "
                "this CLI only accepts velocity-like axes or frequency-like axes with RESTFRQ/RESTFREQ."
            )
        raise ValueError(
            f"Selected cube_ext={cube_ext!r} does not provide a CLI-usable spectral axis "
            "(need velocity-like axis or frequency-like axis with RESTFRQ/RESTFREQ)."
        )
    idx = int(cube_ext) if isinstance(cube_ext, int) else hdul.index_of(cube_ext)
    return np.asarray(hdu.data), hdu.header.copy(), idx


def _read_linefree_prior(hdul: fits.HDUList) -> Optional[np.ndarray]:
    for name in ("LINEFREE", "LINEFREE_USED", "LINEFREE_PRIOR"):
        if name in hdul:
            arr = np.asarray(hdul[name].data)
            arr = np.squeeze(arr)
            if arr.ndim == 1:
                return arr.astype(bool)
    return None


def _read_ripple_prior(hdul: fits.HDUList) -> Optional[np.ndarray]:
    for name in ("RIPFREQ", "RIPFREQ_USED", "RIPFREQ_PRIOR"):
        if name not in hdul:
            continue
        tbl = hdul[name].data
        if tbl is None:
            continue
        if "FREQ_CYC_PER_CH" in getattr(tbl, "names", []):
            arr = np.asarray(tbl["FREQ_CYC_PER_CH"], dtype=float)
            if arr.size > 0:
                return arr
    return None


def _find_spectral_fits_axis(header: fits.Header) -> Optional[int]:
    """Return the FITS axis number (1-based) that appears to be spectral."""
    for ax in (1, 2, 3):
        ctype = str(header.get(f"CTYPE{ax}", "")).upper()
        if any(tok in ctype for tok in ("FREQ", "VRAD", "VELO", "VOPT", "WAVE", "AWAV")):
            return ax
    return None


def _cli_axis_kind(header: fits.Header, spectral_fits_axis: Optional[int] = None) -> Optional[str]:
    """Return spectral axis kind for CLI use: 'velocity', 'frequency', 'wavelength', or None."""
    ax = _find_spectral_fits_axis(header) if spectral_fits_axis is None else int(spectral_fits_axis)
    if ax is None:
        return None
    ctype = str(header.get(f"CTYPE{ax}", "")).upper()
    if ("VRAD" in ctype) or ("VELO" in ctype) or ("VOPT" in ctype):
        return "velocity"
    if "FREQ" in ctype:
        return "frequency"
    if ("WAVE" in ctype) or ("AWAV" in ctype):
        return "wavelength"
    return None


def _is_cli_usable_spectral_hdu(header: fits.Header) -> bool:
    """
    Return True only when the CLI can *actually* build a km/s axis.

    This is intentionally stricter than a simple CTYPE/RESTFRQ check. In particular,
    auto-selection must reject cubes whose spectral metadata looks plausible but whose
    CUNIT is invalid or incompatible, because otherwise the CLI can pick that cube first
    and fail before reaching a later valid candidate.
    """
    fits_ax = _find_spectral_fits_axis(header)
    if fits_ax is None:
        return False

    kind = _cli_axis_kind(header, spectral_fits_axis=fits_ax)
    if kind is None or kind == "wavelength":
        return False

    if kind == "frequency":
        restfreq = header.get("RESTFREQ", header.get("RESTFRQ", 0.0))
        try:
            if float(restfreq) <= 0.0:
                return False
        except Exception:
            return False

    # Validate the actual unit/conversion path rather than trusting CTYPE alone.
    try:
        build_v_axis_kms_from_header(header, 1, spectral_fits_axis=fits_ax)
        return True
    except Exception:
        return False


def _try_make_session(cube_raw: np.ndarray, header: fits.Header) -> Tuple[BaselineSession, np.ndarray]:
    """
    Build a BaselineSession while being robust to non-standard FITS axis order.

    Preference order for nchan/spectral-axis inference:
      1) header spectral axis (CTYPE*)

    Important:
      We do not rely on shape matching alone because a spatial dimension can equal nchan.
      If the spectral axis cannot be identified from the header, we fail loudly rather than guess.
    """
    shape = tuple(int(v) for v in cube_raw.shape)
    ndim = len(shape)
    spectral_axis = _find_spectral_fits_axis(header)
    if spectral_axis is None:
        raise ValueError(
            "Could not identify the spectral FITS axis from CTYPE1..3. "
            "Refusing to guess axis order from shape alone."
        )

    fits_axis_candidates: List[int] = [int(spectral_axis)]

    errors: List[str] = []
    for fits_ax in fits_axis_candidates:
        np_ax = ndim - int(fits_ax)
        if not (0 <= np_ax < ndim):
            continue
        nchan = int(shape[np_ax])
        try:
            v_axis = build_v_axis_kms_from_header(header, nchan, spectral_fits_axis=fits_ax)
            if np_ax == 0:
                axis_order_hint = "v_y_x"
            elif np_ax == 1:
                axis_order_hint = "y_v_x"
            elif np_ax == 2:
                axis_order_hint = "y_x_v"
            else:
                raise ValueError(f"Could not map numpy axis {np_ax} for shape={shape}")
            session = BaselineSession(cube_raw, v_axis, axis_order_hint=axis_order_hint)
            logging.info(
                "Initialized session with spectral FITS axis=%d, nchan=%d (axis_order_in=%s).",
                fits_ax,
                nchan,
                session.axis_order_in,
            )
            return session, v_axis
        except Exception as exc:
            errors.append(f"fits_axis={fits_ax}, nchan={nchan}: {exc}")

    raise ValueError("Could not initialize BaselineSession from header/cube shape. " + " | ".join(errors))


def _convert_axis_values(values: np.ndarray, cunit_raw: str, target_unit: u.Unit) -> np.ndarray:
    """Convert axis values from FITS CUNIT into target_unit; refuse ambiguous nonblank units."""
    txt = (cunit_raw or '').strip()
    if txt == '':
        return np.asarray(values, dtype=float)

    alias = {
        'm/s': 'm / s', 'm s-1': 'm / s', 'ms-1': 'm / s', 'm/sec': 'm / s',
        'km/s': 'km / s', 'km s-1': 'km / s', 'kms-1': 'km / s', 'kmsec-1': 'km / s',
        'hz': 'Hz', 'khz': 'kHz', 'mhz': 'MHz', 'ghz': 'GHz',
    }
    unit_txt = alias.get(txt.lower(), txt)
    try:
        unit_in = u.Unit(unit_txt)
    except Exception as e:
        raise ValueError(f"Unrecognized spectral CUNIT={cunit_raw!r}; refusing to guess conversion to {target_unit!r}.") from e
    if not unit_in.is_equivalent(target_unit):
        raise ValueError(f"Spectral CUNIT={cunit_raw!r} is not compatible with target unit {target_unit!r}.")
    return np.asarray((np.asarray(values, dtype=float) * unit_in).to_value(target_unit), dtype=float)


def build_v_axis_kms_from_header(header: fits.Header, nchan: int, *, spectral_fits_axis: int = 3) -> np.ndarray:
    """
    Build velocity axis in km/s from the spectral FITS axis keywords.

    - If the spectral axis already carries velocity units, converts into km/s.
    - If the spectral axis is frequency-like, first converts into Hz (respecting CUNIT when possible)
      and then converts to radio velocity using RESTFRQ/RESTFREQ.
    """
    axnum = int(spectral_fits_axis)
    if axnum not in (1, 2, 3):
        raise ValueError(f"spectral_fits_axis must be 1, 2, or 3; got {spectral_fits_axis!r}")

    crval = header.get(f"CRVAL{axnum}", 0.0)
    cdelt = header.get(f"CDELT{axnum}", 1.0)
    crpix = header.get(f"CRPIX{axnum}", 1.0)
    ax = crval + (np.arange(nchan) + 1 - crpix) * cdelt

    ctype = str(header.get(f"CTYPE{axnum}", "")).upper()
    cunit_raw = str(header.get(f"CUNIT{axnum}", "")).strip()
    cunit = cunit_raw.lower().replace(" ", "")

    velocity_like = ("VRAD" in ctype) or ("VELO" in ctype) or ("VOPT" in ctype)
    if velocity_like or cunit in ("km/s", "kms-1", "kmsec-1", "m/s", "ms-1", "m/sec"):
        # FITS CUNIT default is SI for dimensional coordinates; for velocity-like axes
        # a blank CUNIT is therefore interpreted as m/s, not km/s.
        ax_kms = _convert_axis_values(ax, cunit_raw or 'm/s', u.km / u.s)
        return np.asarray(ax_kms, dtype=float)

    if ("FREQ" in ctype) or cunit in ("hz", "khz", "mhz", "ghz"):
        ax_hz = _convert_axis_values(ax, cunit_raw or 'Hz', u.Hz)
        restfreq = header.get("RESTFREQ", header.get("RESTFRQ", 0.0))
        if restfreq <= 0:
            raise ValueError("Spectral axis is frequency-like but RESTFRQ/RESTFREQ is missing; cannot build km/s axis.")
        v = (float(restfreq) - ax_hz) / float(restfreq) * 299792.458
        return np.asarray(v, dtype=float)

    if ("WAVE" in ctype) or ("AWAV" in ctype):
        raise ValueError(
            "Spectral axis is wavelength-like (WAVE/AWAV). This CLI expects velocity-like or frequency-like axes."
        )

    raise ValueError(
        f"Unsupported spectral axis for km/s construction: CTYPE{axnum}={ctype!r}, CUNIT{axnum}={cunit_raw!r}."
    )


def _read_target_mask_2d(hdul: fits.HDUList, ext: Optional[Union[int, str]]) -> Optional[np.ndarray]:
    if ext is None:
        return None
    hdu = hdul[ext]
    arr = np.asarray(hdu.data)
    arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"target mask HDU must be 2D after squeeze, got shape={arr.shape}")
    return np.asarray(arr != 0, dtype=bool)


def _append_linefree_hdu(hdul_out: fits.HDUList, name: str, data_1d: Sequence[bool], *, comment: str) -> None:
    hdr_lf = fits.Header()
    hdr_lf["BTYPE"] = "LineFreeMask"
    hdr_lf["COMMENT"] = comment
    hdul_out.append(fits.ImageHDU(data=np.asarray(data_1d, dtype=np.uint8), header=hdr_lf, name=name))


def _append_ripple_hdu(hdul_out: fits.HDUList, name: str, freqs: Sequence[float], *, comment: str) -> None:
    arr = np.asarray(freqs, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        period = np.where(arr != 0.0, 1.0 / arr, np.inf)
    cols = [
        fits.Column(name="FREQ_CYC_PER_CH", format="D", array=arr),
        fits.Column(name="PERIOD_CH", format="D", array=period),
    ]
    tbhdu = fits.BinTableHDU.from_columns(cols, name=name)
    tbhdu.header["COMMENT"] = comment
    hdul_out.append(tbhdu)


def _replace_or_append_hdu(hdul: fits.HDUList, hdu_new: fits.hdu.base.ExtensionHDU) -> None:
    """Replace an existing HDU with the same EXTNAME, or append if absent."""
    name = str(getattr(hdu_new, "name", "") or "").upper()
    if not name:
        hdul.append(hdu_new)
        return
    for i, hdu in enumerate(hdul):
        if str(getattr(hdu, "name", "") or "").upper() == name:
            hdul[i] = hdu_new
            return
    hdul.append(hdu_new)


def _remove_named_hdus(hdul: fits.HDUList, names: Sequence[str], *, protect_hdu: object | None = None) -> None:
    """Remove all HDUs whose EXTNAME matches one of `names` (case-insensitive)."""
    names_up = {str(n).upper() for n in names}
    remove_idx = []
    for i, hdu in enumerate(hdul):
        if i == 0:
            continue
        if protect_hdu is not None and hdu is protect_hdu:
            continue
        hname = str(getattr(hdu, 'name', '') or '').upper()
        if hname in names_up:
            remove_idx.append(i)
    for i in reversed(remove_idx):
        del hdul[i]


def _strip_checksum_all_hdus(hdul: fits.HDUList) -> None:
    """Remove CHECKSUM/DATASUM cards from all HDUs in-place."""
    for hdu in hdul:
        hdr = getattr(hdu, "header", None)
        if hdr is None:
            continue
        for key in ("CHECKSUM", "DATASUM", "ZHECKSUM", "ZDATASUM"):
            while key in hdr:
                del hdr[key]


def run_cli_pipeline(
    input_fits: str,
    output_fits: str,
    *,
    cube_ext: Optional[Union[int, str]] = None,
    iterations: int = 2,
    poly_order: int = 1,
    auto_linefree: bool = True,
    linefree_cfg: LineFreeConfig = LineFreeConfig(),
    manual_v_windows: Optional[List[str]] = None,
    linefree_mode: str = "auto",
    linefree_scope: str = "global",
    target_mask_ext: Optional[Union[int, str]] = None,
    enable_ripple: bool = True,
    ripple_cfg: RippleConfig = RippleConfig(),
    ripple_mode: str = "auto",
    ripple_scope: str = "global",
    robust: bool = False,
    chunk_pix: int = 65536,
    reproducible_mode: bool = False,
    load_prior_from_input: bool = True,
    run_cube_analysis: bool = False,
    cube_analysis_method: str = "smooth_mask_lite",
    cube_analysis_write_provisional: bool = True,
    write_diagnostics: bool = False,
    diagnostics_prefix: Optional[str] = None,
    write_profile: bool = False,
    profile_prefix: Optional[str] = None,
) -> None:
    """
    Pipeline:
      read FITS -> iterate baseline fit -> write baselined FITS with QC HDUs
      optionally run cube_analysis to append provisional masks, MASK3D, and moments.
    """
    _profile_records: List[Dict[str, Any]] = []
    _t0 = time.perf_counter()
    _tprev = _t0
    if write_profile:
        _tprev = _profile_record(_profile_records, "start", _tprev, _t0, input_fits=str(input_fits), output_fits=str(output_fits), iterations=int(iterations))
    logging.info("Loading FITS: %s", input_fits)
    with fits.open(input_fits, memmap=True) as hdul:
        cube_raw, header, cube_idx = _get_cube_from_hdul(hdul, cube_ext=cube_ext)
        lf_prior = _read_linefree_prior(hdul) if load_prior_from_input else None
        rf_prior = _read_ripple_prior(hdul) if load_prior_from_input else None
        target_mask_2d = _read_target_mask_2d(hdul, target_mask_ext)

    if cube_raw.ndim != 3:
        raise ValueError(f"Cube must be 3D, got shape={cube_raw.shape}")

    session, _v_axis = _try_make_session(cube_raw, header)
    if write_profile:
        _tprev = _profile_record(_profile_records, "read_and_init_session", _tprev, _t0, cube_shape=list(np.asarray(cube_raw).shape), session_shape=[int(session.nchan), int(session.ny), int(session.nx)])

    if target_mask_2d is not None and target_mask_2d.shape != (session.ny, session.nx):
        raise ValueError(
            f"target_mask_2d shape mismatch: {target_mask_2d.shape} vs {(session.ny, session.nx)}"
        )

    if lf_prior is not None and np.asarray(lf_prior).shape != (session.nchan,):
        logging.warning(
            "Ignoring prior LINEFREE because shape=%s does not match nchan=%d.",
            np.asarray(lf_prior).shape,
            session.nchan,
        )
        lf_prior = None

    if load_prior_from_input and (lf_prior is not None or rf_prior is not None):
        session.set_prior(linefree_mask_1d=lf_prior, ripple_freqs=rf_prior, source=input_fits)
        logging.info(
            "Loaded prior baseline info from input FITS: LINEFREE=%s RIPFREQ=%s",
            lf_prior is not None,
            rf_prior is not None,
        )

    for it in range(int(iterations)):
        logging.info("--- Iteration %d / %d ---", it + 1, int(iterations))
        run_one_iteration(
            session=session,
            target_mask_2d=target_mask_2d,
            auto_linefree=auto_linefree,
            linefree_cfg=linefree_cfg,
            manual_v_windows=manual_v_windows,
            linefree_mode=linefree_mode,
            linefree_scope=linefree_scope,
            enable_ripple=enable_ripple,
            ripple_cfg=ripple_cfg,
            ripple_mode=ripple_mode,
            ripple_scope=ripple_scope,
            poly_order=int(poly_order),
            robust=robust,
            chunk_pix=int(chunk_pix),
            reproducible_mode=bool(reproducible_mode),
        )
        rms_map = session.fit_stats.get("rms_map", None)
        if rms_map is not None:
            logging.info("  -> Median BASE_RMS: %.4f", float(np.nanmedian(rms_map)))
        if write_profile:
            _tprev = _profile_record(_profile_records, f"iteration_{it + 1}", _tprev, _t0, median_base_rms=(None if rms_map is None else float(np.nanmedian(rms_map))))

    logging.info("Writing baselined cube: %s", output_fits)
    cube_out = session.get_full_cube_work(cache=False).astype(np.float32, copy=False)
    if write_profile:
        _tprev = _profile_record(_profile_records, "assemble_output_cube", _tprev, _t0)
    axis_order_in = session.axis_order_in

    linefree_mask_1d_prior = None if session.linefree_mask_1d_prior is None else np.asarray(session.linefree_mask_1d_prior, dtype=bool).copy()
    linefree_mask_1d_current = None if session.linefree_mask_1d_current is None else np.asarray(session.linefree_mask_1d_current, dtype=bool).copy()
    ripple_freqs_prior = None if session.ripple_freqs_prior is None else np.asarray(session.ripple_freqs_prior, dtype=float).copy()
    ripple_freqs_current = None if session.ripple_freqs_current is None else np.asarray(session.ripple_freqs_current, dtype=float).copy()
    rms_map = session.fit_stats.get("rms_map", None)
    if rms_map is not None:
        rms_map = np.asarray(rms_map, dtype=np.float32)
    flag_map = session.fit_stats.get("flag_map", None)
    if flag_map is not None:
        flag_map = np.asarray(flag_map, dtype=np.uint8)

    if axis_order_in == "y_x_v":
        cube_out_write = np.transpose(cube_out, (1, 2, 0))
    elif axis_order_in == "y_v_x":
        cube_out_write = np.transpose(cube_out, (1, 0, 2))
    else:
        cube_out_write = cube_out

    session.clear_work_cache()
    del session

    with fits.open(input_fits, memmap=True) as hdul_in:
        hdul_out = fits.HDUList([hdu.copy() for hdu in hdul_in])
    if write_profile:
        _tprev = _profile_record(_profile_records, "copy_output_hdus", _tprev, _t0, nhdu=int(len(hdul_out)))

    target_cube_hdu = hdul_out[cube_idx]
    target_cube_name = str(getattr(target_cube_hdu, "name", "") or "").strip()

    hdul_out[cube_idx].data = cube_out_write

    # Remove stale baseline/provenance and downstream analysis products from copied input FITS.
    _remove_named_hdus(
        hdul_out,
        [
            "LINEFREE",
            "LINEFREE_USED",
            "LINEFREE_PRIOR",
            "RIPFREQ",
            "RIPFREQ_USED",
            "RIPFREQ_PRIOR",
            "BASE_RMS",
            "BASE_FLG",
            "BASEFLAG",
            "BASE_DIAG",
            "RMS",
            "MASK3D",
            "MOMENT0",
            "BASESUP3D",
            "LINECAND3D",
            "MOM0_BASESUP",
            "MOM0_LINECAND",
        ],
        protect_hdu=target_cube_hdu,
    )

    # Resolve the output cube locator after any HDU removals. Integer ext indices can
    # shift when stale HDUs are deleted, so we track the actual copied cube HDU object
    # and then re-discover its post-cleanup location.
    cube_idx_out = None
    for i, hdu in enumerate(hdul_out):
        if hdu is target_cube_hdu:
            cube_idx_out = i
            break

    # Fallback: identity should normally survive HDU deletions, but if not, try the
    # original EXTNAME before giving up.
    if cube_idx_out is None and target_cube_name:
        try:
            cube_idx_out = hdul_out.index_of(target_cube_name)
        except Exception:
            cube_idx_out = None

    if cube_idx_out is None:
        raise RuntimeError("Could not resolve output cube HDU after cleanup.")

    cube_ext_out: Optional[Union[int, str]]
    if cube_idx_out == 0:
        cube_ext_out = None
    else:
        # Use the resolved integer index for an unambiguous hand-off to cube_analysis.
        cube_ext_out = cube_idx_out

    _strip_checksum_all_hdus(hdul_out)

    if linefree_mask_1d_prior is not None:
        hdr = fits.Header()
        hdr["BTYPE"] = "LineFreeMask"
        hdr["COMMENT"] = "1=line-free prior loaded from input FITS; 0=non-prior/line candidate."
        _replace_or_append_hdu(
            hdul_out,
            fits.ImageHDU(data=np.asarray(linefree_mask_1d_prior, dtype=np.uint8), header=hdr, name="LINEFREE_PRIOR"),
        )
    if linefree_mask_1d_current is not None:
        hdr = fits.Header()
        hdr["BTYPE"] = "LineFreeMask"
        hdr["COMMENT"] = "1=line-free channel used for baseline fitting; 0=line/ignored."
        _replace_or_append_hdu(
            hdul_out,
            fits.ImageHDU(data=np.asarray(linefree_mask_1d_current, dtype=np.uint8), header=hdr, name="LINEFREE"),
        )
        _replace_or_append_hdu(
            hdul_out,
            fits.ImageHDU(data=np.asarray(linefree_mask_1d_current, dtype=np.uint8), header=hdr.copy(), name="LINEFREE_USED"),
        )

    if enable_ripple and ripple_freqs_prior is not None and len(ripple_freqs_prior) > 0:
        arr = np.asarray(ripple_freqs_prior, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            period = np.where(arr != 0.0, 1.0 / arr, np.inf)
        cols = [
            fits.Column(name="FREQ_CYC_PER_CH", format="D", array=arr),
            fits.Column(name="PERIOD_CH", format="D", array=period),
        ]
        tbhdu = fits.BinTableHDU.from_columns(cols, name="RIPFREQ_PRIOR")
        tbhdu.header["COMMENT"] = "Ripple frequencies loaded from input FITS as priors."
        _replace_or_append_hdu(hdul_out, tbhdu)

    freqs = ripple_freqs_current
    if enable_ripple and (freqs is not None) and (len(freqs) > 0):
        arr = np.asarray(freqs, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            period = np.where(arr != 0.0, 1.0 / arr, np.inf)
        cols = [
            fits.Column(name="FREQ_CYC_PER_CH", format="D", array=arr),
            fits.Column(name="PERIOD_CH", format="D", array=period),
        ]
        tbhdu = fits.BinTableHDU.from_columns(cols, name="RIPFREQ")
        tbhdu.header["COMMENT"] = "Ripple frequencies used in baseline model."
        _replace_or_append_hdu(hdul_out, tbhdu)
        tbhdu_used = fits.BinTableHDU.from_columns(cols, name="RIPFREQ_USED")
        tbhdu_used.header["COMMENT"] = "Ripple frequencies used in baseline model."
        _replace_or_append_hdu(hdul_out, tbhdu_used)

    if rms_map is not None:
        hdr = fits.Header()
        hdr["BTYPE"] = "BaselineResidRMS"
        hdr["BUNIT"] = "K"
        _replace_or_append_hdu(
            hdul_out,
            fits.ImageHDU(data=np.asarray(rms_map, dtype=np.float32), header=hdr, name="BASE_RMS"),
        )

    if flag_map is not None:
        hdr = fits.Header()
        hdr["BTYPE"] = "BaselineFitFlag"
        hdr["COMMENT"] = "0=OK, 1=insufficient finite line-free points, 2=lstsq failed (spectrum left unchanged), 3=all-NaN spectrum."
        _replace_or_append_hdu(
            hdul_out,
            fits.ImageHDU(data=np.asarray(flag_map, dtype=np.uint8), header=hdr, name="BASE_FLG"),
        )

    if write_profile:
        _tprev = _profile_record(_profile_records, "prepare_output_hdus", _tprev, _t0)
    hdul_out.writeto(output_fits, overwrite=True)
    if write_profile:
        _tprev = _profile_record(_profile_records, "write_fits", _tprev, _t0)
    try:
        hdul_out.close()
    except Exception:
        pass
    del cube_out
    del cube_out_write

    if write_diagnostics:
        payload = build_baseline_diagnostic_payload(
            linefree_mask=(np.ones((0,), dtype=bool) if linefree_mask_1d_current is None else linefree_mask_1d_current),
            ripple_freqs=ripple_freqs_current,
            rms_map=rms_map,
            flag_map=flag_map,
            extra={
                "input_fits": input_fits,
                "output_fits": output_fits,
                "cube_ext": cube_ext,
                "iterations": int(iterations),
                "poly_order": int(poly_order),
                "linefree_mode": linefree_mode,
                "linefree_scope": linefree_scope,
                "manual_v_windows": manual_v_windows,
                "ripple_mode": ripple_mode,
                "ripple_scope": ripple_scope,
                "enable_ripple": bool(enable_ripple),
                "robust": bool(robust),
                "chunk_pix": int(chunk_pix),
                "reproducible_mode": bool(reproducible_mode),
                "linefree_cfg": linefree_cfg.__dict__,
                "ripple_cfg": ripple_cfg.__dict__,
            },
        )
        prefix = diagnostics_prefix if diagnostics_prefix else output_fits
        json_path, txt_path = write_baseline_diagnostic_files(prefix, payload)
        logging.info("Wrote baseline diagnostics: %s | %s", json_path, txt_path)

    if write_profile:
        prefix = profile_prefix if profile_prefix else output_fits
        json_path, txt_path = _write_profile_sidecar(prefix, _profile_records)
        logging.info("Wrote baseline profile: %s | %s", json_path, txt_path)

    if run_cube_analysis:
        from sd_radio_spectral_fits.map.cube_analysis import make_3d_mask_for_existing_fits

        logging.info("Running cube_analysis on baselined cube (method=%s).", cube_analysis_method)
        make_3d_mask_for_existing_fits(
            output_fits,
            output_fits=output_fits,
            cube_ext=cube_ext_out,
            method=cube_analysis_method,
            file_format="fits",
            require_kms=False,
            rms_ext="auto",
            write_provisional_masks=bool(cube_analysis_write_provisional),
        )
        logging.info("cube_analysis done.")

    logging.info("Done.")


def _main() -> None:
    p = argparse.ArgumentParser(description="Iterative cube baseline subtraction (poly + multi-ripple).")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--cube-ext", default=None)
    p.add_argument("--iter", type=int, default=2)
    p.add_argument("--poly", type=int, default=1)
    p.add_argument("--no-auto-linefree", action="store_true")
    p.add_argument("--manual-vwin", nargs="*", default=None)
    p.add_argument("--linefree-mode", default="auto", choices=["auto", "prior", "current", "or"])
    p.add_argument("--linefree-scope", default="global", choices=["global", "target", "both"])
    p.add_argument("--target-mask-ext", default=None, help="2D FITS HDU (name or index) used as target_mask_2d for target/both scopes.")
    p.add_argument("--no-ripple", action="store_true")
    p.add_argument("--ripple-mode", default="auto", choices=["auto", "prior", "current"])
    p.add_argument("--ripple-scope", default="global", choices=["global", "target", "both"])
    p.add_argument("--no-load-prior", action="store_true")
    p.add_argument("--nfreq", type=int, default=2)
    p.add_argument("--period", nargs=2, type=float, default=(20.0, 400.0), metavar=("LO", "HI"))
    p.add_argument("--linefree-conservative-q", type=float, default=None, help="Enable conservative line-candidate detection using an additional quantile spectrum, e.g. 0.9.")
    p.add_argument("--linefree-one-tail", action="store_true", help="When using --linefree-conservative-q, inspect only the upper quantile (skip the symmetric lower quantile).")
    p.add_argument("--robust", action="store_true")
    p.add_argument("--chunk-pix", type=int, default=65536)
    p.add_argument("--reproducible", action="store_true", help="Use deterministic baseline settings (float64, normalized polynomial basis, stride sampling, stable peak ranking).")
    p.add_argument("--run-analysis", action="store_true")
    p.add_argument("--no-provisional", action="store_true")
    p.add_argument("--write-diagnostics", action="store_true", help="Write <prefix>.baseline_diag.json/txt for cross-environment comparison.")
    p.add_argument("--diagnostics-prefix", default=None, help="Prefix (or FITS path) for diagnostic sidecar files.")
    p.add_argument("--write-profile", action="store_true", help="Write <prefix>.baseline_profile.json/txt with stage timings and RSS.")
    p.add_argument("--profile-prefix", default=None, help="Prefix (or FITS path) for profiling sidecar files.")
    p.add_argument(
        "--analysis-method",
        default="smooth_mask_lite",
        choices=["smooth_mask", "smooth_mask_lite", "simple", "derivative"],
    )
    args = p.parse_args()

    cube_ext = None
    if args.cube_ext is not None:
        try:
            cube_ext = int(args.cube_ext)
        except Exception:
            cube_ext = str(args.cube_ext)

    lcfg = LineFreeConfig(
        conservative_quantile=(None if args.linefree_conservative_q is None else float(args.linefree_conservative_q)),
        conservative_both_tails=(not bool(args.linefree_one_tail)),
    )
    rcfg = RippleConfig(nfreq=int(args.nfreq), period_range_chan=(float(args.period[0]), float(args.period[1])))

    target_mask_ext = None
    if args.target_mask_ext is not None:
        try:
            target_mask_ext = int(args.target_mask_ext)
        except Exception:
            target_mask_ext = str(args.target_mask_ext)

    run_cli_pipeline(
        args.input,
        args.output,
        cube_ext=cube_ext,
        iterations=int(args.iter),
        poly_order=int(args.poly),
        auto_linefree=not args.no_auto_linefree,
        linefree_cfg=lcfg,
        manual_v_windows=args.manual_vwin if args.manual_vwin else None,
        linefree_mode=str(args.linefree_mode),
        linefree_scope=str(args.linefree_scope),
        target_mask_ext=target_mask_ext,
        enable_ripple=not args.no_ripple,
        ripple_cfg=rcfg,
        ripple_mode=str(args.ripple_mode),
        ripple_scope=str(args.ripple_scope),
        robust=bool(args.robust),
        chunk_pix=int(args.chunk_pix),
        reproducible_mode=bool(args.reproducible),
        load_prior_from_input=not args.no_load_prior,
        run_cube_analysis=bool(args.run_analysis),
        cube_analysis_method=str(args.analysis_method),
        cube_analysis_write_provisional=not args.no_provisional,
        write_diagnostics=bool(args.write_diagnostics),
        diagnostics_prefix=args.diagnostics_prefix,
        write_profile=bool(args.write_profile),
        profile_prefix=args.profile_prefix,
    )


if __name__ == "__main__":
    _main()
