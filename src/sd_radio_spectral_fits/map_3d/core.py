# src/sd_radio_spectral_fits/otf/core.py
from __future__ import annotations

import builtins
import warnings
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import j1
from scipy.signal import fftconvolve

from .config import MapConfig, GridInput, GridResult


SQRT_LN2 = float(np.sqrt(np.log(2.0)))
FIRST_JINC_NULL_OVER_PI = 3.8317059702075125 / np.pi

# Mangum et al. (2007) / Sawada et al. (2008) と整合する beam-aware default
# cell = beam / 3 のとき、CASA / NOSTAR でよく使われる pixel 既定値に戻る。
#
# GJINC (Mangum/Sawada):
#   jwidth = 1.55 * (beam / 3)
#   gwidth = 2.52 * sqrt(log(2)) * (beam / 3)
#   support radius = beam FWHM
#
# GAUSS (Sawada):
#   c(r) = exp(-(r/a)^2),  a = beam / 3
#   実装上の gwidth は HWHM なので gwidth = sqrt(log(2)) * a
#   support radius = 3a = 3 * gwidth / sqrt(log(2))
DEFAULT_GJINC_JWIDTH_BEAM = 1.55 / 3.0
DEFAULT_GJINC_GWIDTH_BEAM = 2.52 * SQRT_LN2 / 3.0
DEFAULT_GAUSS_GWIDTH_BEAM = SQRT_LN2 / 3.0
RECOMMENDED_MAX_CELL_OVER_BEAM = 1.0 / 3.0


# ==========================================
# 1. Validation / Normalization Helpers
# ==========================================

def _ensure_1d(name: str, arr: np.ndarray, ndump: int) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim != 1 or a.shape[0] != ndump:
        raise ValueError(f"{name} must have shape ({ndump},), got {a.shape}")
    return a


def _safe_bool_array(values, default: bool = False) -> np.ndarray:
    """Convert mixed boolean-like values into a strict bool ndarray."""
    a = np.asarray(values)
    if a.dtype.kind == "b":
        return a.astype(bool, copy=False)

    out = np.full(a.shape, bool(default), dtype=bool)
    if a.dtype.kind in "iuf":
        finite = np.isfinite(a)
        out[finite] = a[finite] != 0
        return out

    flat = a.astype(object, copy=False).ravel()
    flat_out = out.ravel()
    true_vals = {"true", "t", "1", "yes", "y", "on"}
    false_vals = {"false", "f", "0", "no", "n", "off", "", "nan", "none", "null", "<na>"}
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


def _validate_input(input_data: GridInput) -> None:
    x = np.asarray(input_data.x)
    y = np.asarray(input_data.y)
    spec = np.asarray(input_data.spec)
    flag = np.asarray(input_data.flag)
    time = np.asarray(input_data.time)

    if spec.ndim != 2:
        raise ValueError(f"spec must be 2D (ndump, nchan), got shape={spec.shape}")
    ndump = spec.shape[0]

    _ensure_1d("x", x, ndump)
    _ensure_1d("y", y, ndump)
    _ensure_1d("flag", flag, ndump)
    _ensure_1d("time", time, ndump)

    if input_data.rms is not None:
        _ensure_1d("rms", np.asarray(input_data.rms), ndump)
    if input_data.tint is not None:
        _ensure_1d("tint", np.asarray(input_data.tint), ndump)
    if input_data.tsys is not None:
        _ensure_1d("tsys", np.asarray(input_data.tsys), ndump)
    if input_data.is_turnaround is not None:
        _ensure_1d("is_turnaround", np.asarray(input_data.is_turnaround), ndump)


def _coerce_optional_bool(value, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y", "on"}


def _runtime_options(config: MapConfig) -> dict[str, object]:
    reproducible = _coerce_optional_bool(
        getattr(config, "_reproducible_mode", getattr(config, "reproducible_mode", getattr(config, "reproducible", False))),
        default=False,
    )

    if hasattr(config, "_sort_neighbors"):
        sort_neighbors = _coerce_optional_bool(getattr(config, "_sort_neighbors"), default=reproducible)
    else:
        sort_neighbors = True if reproducible else _coerce_optional_bool(getattr(config, "sort_neighbors", False), default=False)

    if reproducible:
        dtype_name = "float64"
        workers = 1
    else:
        dtype_name = str(getattr(config, "dtype", "float32"))
        workers = int(getattr(config, "_workers", getattr(config, "workers", -1)))

    return {
        "reproducible": reproducible,
        "sort_neighbors": sort_neighbors,
        "dtype_name": dtype_name,
        "workers": workers,
    }


def _normalize_backend_and_kernel(config: MapConfig) -> tuple[str, str]:
    backend = str(config.backend)
    kernel = str(config.kernel).lower()

    # バックエンドとカーネルのエイリアス解決
    if backend == "numpy_gjinc":
        backend = "numpy"
        if kernel != "gjinc":
            raise ValueError("backend='numpy_gjinc' requires kernel='gjinc'")
    elif backend == "numpy_gauss":
        backend = "numpy"
        if kernel != "gauss":
            raise ValueError("backend='numpy_gauss' requires kernel='gauss'")
    elif backend == "cygrid_gauss":
        backend = "cygrid"
        if kernel != "gauss":
            raise ValueError("backend='cygrid_gauss' requires kernel='gauss'")

    if kernel not in ("gjinc", "gauss"):
        raise ValueError(f"Unknown kernel: {config.kernel!r}. Use 'gjinc' or 'gauss'.")

    return backend, kernel



def _validate_and_resolve_config(config: MapConfig) -> dict[str, float | str | bool]:
    backend_impl, kernel_name = _normalize_backend_and_kernel(config)
    runtime = _runtime_options(config)

    if config.estimator == "plane":
        raise NotImplementedError("estimator='plane' is not implemented yet. Please use 'avg'.")

    if backend_impl != "numpy":
        raise NotImplementedError(
            f"backend={backend_impl!r} is not implemented in this build. "
            "Please use backend='numpy' or one of its aliases."
        )

    # 基本的な正の数チェック
    if config.nx <= 0 or config.ny <= 0:
        raise ValueError("nx and ny must be positive")
    if config.cell_arcsec <= 0 or config.beam_fwhm_arcsec <= 0:
        raise ValueError("cell_arcsec and beam_fwhm_arcsec must be positive")
    if config.chunk_ch <= 0:
        raise ValueError("chunk_ch must be positive")
    if runtime["dtype_name"] not in ("float32", "float64"):
        raise ValueError("dtype must be 'float32' or 'float64'")
    if config.weight_clip_quantile is not None and not (0.0 <= float(config.weight_clip_quantile) <= 1.0):
        raise ValueError("weight_clip_quantile must be within [0, 1]")
    if config.weight_clip_max is not None and float(config.weight_clip_max) <= 0:
        raise ValueError("weight_clip_max must be positive when specified")
    if float(config.eps_weight_sum) <= 0:
        raise ValueError("eps_weight_sum must be positive")
    if float(getattr(config, 'min_abs_weight_ratio', 0.0)) < 0:
        raise ValueError("min_abs_weight_ratio must be >= 0")
    if float(getattr(config, 'min_cancel_ratio', 0.0)) < 0 or float(getattr(config, 'min_cancel_ratio', 0.0)) > 1:
        raise ValueError("min_cancel_ratio must be within [0, 1]")

    beam_pix = float(config.beam_fwhm_arcsec / config.cell_arcsec)
    cell_over_beam = float(config.cell_arcsec / config.beam_fwhm_arcsec)
    cell_is_coarse = bool(cell_over_beam > RECOMMENDED_MAX_CELL_OVER_BEAM)

    if getattr(config, "warn_if_cell_coarse", True) and cell_is_coarse:
        warnings.warn(
            "cell_arcsec is coarser than beam_fwhm_arcsec/3. "
            "This may increase aliasing and degrade the gridded beam. "
            f"Got cell={float(config.cell_arcsec):.6g} arcsec, "
            f"beam={float(config.beam_fwhm_arcsec):.6g} arcsec "
            f"(cell/beam={cell_over_beam:.4f}).",
            RuntimeWarning,
            stacklevel=2,
        )

    kernel_preset = str(getattr(config, "kernel_preset", "mangum2007")).lower()
    if kernel_preset not in ("mangum2007", "legacy"):
        raise ValueError("kernel_preset must be 'mangum2007' or 'legacy'")

    kernel_sign = str(getattr(config, "kernel_sign", "auto")).lower()
    if kernel_sign not in ("auto", "signed", "positive_only"):
        raise ValueError("kernel_sign must be 'auto', 'signed', or 'positive_only'")
    if kernel_sign == "auto":
        kernel_sign = "signed" if kernel_preset == "mangum2007" else "positive_only"
    if kernel_preset == "mangum2007" and kernel_sign != "signed":
        warnings.warn(
            "kernel_preset='mangum2007' is literature-oriented and is normally used with kernel_sign='signed'. "
            "You explicitly requested a non-standard combination.",
            RuntimeWarning,
            stacklevel=2,
        )
    if kernel_preset == "legacy" and kernel_sign != "positive_only":
        warnings.warn(
            "kernel_preset='legacy' is normally used with kernel_sign='positive_only'. "
            "You explicitly requested a non-standard combination.",
            RuntimeWarning,
            stacklevel=2,
        )

    def _resolve_support_radius(
        *,
        default_pix: float,
        default_source: str,
        allow_first_null: bool,
        jwidth_pix: float | None,
    ) -> tuple[float, str]:
        if config.support_radius_pix is not None:
            support = float(config.support_radius_pix)
            source = "support_radius_pix"
        elif getattr(config, "support_radius_beam", None) is not None:
            support = float(config.support_radius_beam) * beam_pix
            source = "support_radius_beam"
        elif config.truncate is None:
            support = float(default_pix)
            source = default_source
        elif isinstance(config.truncate, (int, float)):
            support = float(config.truncate)
            source = "truncate_numeric"
        elif config.truncate == "first_null":
            if not allow_first_null or jwidth_pix is None or not np.isfinite(jwidth_pix):
                raise ValueError("truncate='first_null' is only valid for kernel='gjinc'")
            support = float(FIRST_JINC_NULL_OVER_PI * jwidth_pix)
            source = "truncate_first_null"
        else:
            raise ValueError("Invalid truncate mode.")
        if support <= 0:
            raise ValueError("support_radius must be positive")
        return support, source

    # 空間カーネルのパラメータ解決
    if kernel_name == "gjinc":
        default_gwidth_beam = DEFAULT_GJINC_GWIDTH_BEAM if kernel_preset == "mangum2007" else None
        default_jwidth_beam = DEFAULT_GJINC_JWIDTH_BEAM if kernel_preset == "mangum2007" else None
        default_gwidth_pix = None if kernel_preset == "mangum2007" else (2.52 * SQRT_LN2)
        default_jwidth_pix = None if kernel_preset == "mangum2007" else 1.55
    else:
        default_gwidth_beam = DEFAULT_GAUSS_GWIDTH_BEAM if kernel_preset == "mangum2007" else None
        default_jwidth_beam = None
        default_gwidth_pix = None if kernel_preset == "mangum2007" else SQRT_LN2
        default_jwidth_pix = np.nan

    if config.gwidth_pix is not None:
        gwidth_pix = float(config.gwidth_pix)
        gwidth_source = "gwidth_pix"
    elif config.gwidth_beam is not None:
        gwidth_pix = float(config.gwidth_beam) * beam_pix
        gwidth_source = "gwidth_beam"
    elif default_gwidth_beam is not None:
        gwidth_pix = float(default_gwidth_beam * beam_pix)
        gwidth_source = "preset_beam"
    else:
        gwidth_pix = float(default_gwidth_pix)
        gwidth_source = "preset_pix"

    if gwidth_pix <= 0:
        raise ValueError("gwidth must be positive")

    if kernel_name == "gjinc":
        if config.jwidth_pix is not None:
            jwidth_pix = float(config.jwidth_pix)
            jwidth_source = "jwidth_pix"
        elif config.jwidth_beam is not None:
            jwidth_pix = float(config.jwidth_beam) * beam_pix
            jwidth_source = "jwidth_beam"
        elif default_jwidth_beam is not None:
            jwidth_pix = float(default_jwidth_beam * beam_pix)
            jwidth_source = "preset_beam"
        else:
            jwidth_pix = float(default_jwidth_pix)
            jwidth_source = "preset_pix"

        if jwidth_pix <= 0:
            raise ValueError("jwidth must be positive for kernel='gjinc'")
    else:
        jwidth_pix = np.nan
        jwidth_source = "not_used"

    widths_explicit = any(
        v is not None for v in (config.gwidth_pix, config.gwidth_beam, config.jwidth_pix, config.jwidth_beam)
    )

    if kernel_name == "gjinc":
        if kernel_preset == "mangum2007":
            default_support_pix = float(beam_pix)
            default_support_source = "preset_beam_fwhm"
            if widths_explicit and config.support_radius_pix is None and getattr(config, "support_radius_beam", None) is None and config.truncate is None:
                warnings.warn(
                    "kernel='gjinc' with kernel_preset='mangum2007' uses support radius = beam FWHM by default. "
                    "You changed gwidth/jwidth explicitly while leaving support implicit; if you intend a different cutoff radius, set support_radius_* or truncate explicitly.",
                    RuntimeWarning,
                    stacklevel=2,
                )
        else:
            default_support_pix = float(FIRST_JINC_NULL_OVER_PI * jwidth_pix)
            default_support_source = "preset_first_null"
    else:
        if kernel_preset == "mangum2007":
            default_support_pix = float(3.0 * gwidth_pix / SQRT_LN2)
            default_support_source = "preset_3a"
        else:
            default_support_pix = float(3.0 * gwidth_pix)
            default_support_source = "preset_3hwhm"

    support_radius_pix, support_source = _resolve_support_radius(
        default_pix=float(default_support_pix),
        default_source=default_support_source,
        allow_first_null=(kernel_name == "gjinc"),
        jwidth_pix=None if not np.isfinite(jwidth_pix) else float(jwidth_pix),
    )

    return {
        "backend_impl": backend_impl,
        "kernel_name": kernel_name,
        "kernel_preset": kernel_preset,
        "kernel_sign": kernel_sign,
        "gwidth_pix": float(gwidth_pix),
        "jwidth_pix": float(jwidth_pix),
        "support_radius_pix": float(support_radius_pix),
        "beam_pix": float(beam_pix),
        "cell_over_beam": float(cell_over_beam),
        "cell_is_coarse": bool(cell_is_coarse),
        "gwidth_source": gwidth_source,
        "jwidth_source": jwidth_source,
        "support_source": support_source,
    }




def _gaussian_beam_response(r_arcsec: np.ndarray, beam_fwhm_arcsec: float) -> np.ndarray:
    return np.exp(-4.0 * np.log(2.0) * (r_arcsec / beam_fwhm_arcsec) ** 2)


def _estimate_fwhm_from_profile(r_arcsec: np.ndarray, profile: np.ndarray) -> float | None:
    r = np.asarray(r_arcsec, dtype=float)
    p = np.asarray(profile, dtype=float)
    m = np.isfinite(r) & np.isfinite(p) & (r >= 0)
    if np.count_nonzero(m) < 2:
        return None
    r = r[m]
    p = p[m]
    if p.size == 0 or np.nanmax(p) <= 0:
        return None
    order = np.argsort(r)
    r = r[order]
    p = p[order]
    half = 0.5 * np.nanmax(p)
    above = p >= half
    if not np.any(above):
        return None
    i_last = int(np.max(np.where(above)[0]))
    if i_last >= len(r) - 1:
        return None
    r0, r1 = float(r[i_last]), float(r[i_last + 1])
    p0, p1 = float(p[i_last]), float(p[i_last + 1])
    if not np.isfinite(p0) or not np.isfinite(p1) or p0 == p1:
        return None
    frac = (half - p0) / (p1 - p0)
    if not np.isfinite(frac):
        return None
    rh = r0 + frac * (r1 - r0)
    return 2.0 * rh


def _fit_gaussian_core_2d(image: np.ndarray, x2d: np.ndarray, y2d: np.ndarray, min_frac: float = 0.2) -> dict[str, float] | None:
    img = np.asarray(image, dtype=float)
    if img.ndim != 2 or img.size == 0:
        return None
    finite = np.isfinite(img)
    if not np.any(finite):
        return None

    peak_idx = np.nanargmax(np.where(finite, img, -np.inf))
    iy0, ix0 = np.unravel_index(int(peak_idx), img.shape)
    x_ref = float(x2d[iy0, ix0])
    y_ref = float(y2d[iy0, ix0])

    img_peak = float(img[iy0, ix0])
    if not np.isfinite(img_peak) or img_peak <= 0:
        return None

    work = np.where(finite, img / img_peak, np.nan)
    min_frac = float(min_frac)
    min_frac = min(max(min_frac, 1e-3), 0.95)
    thresholds = [min_frac, 0.1, 0.05, 0.02]
    fit = None

    for thr in thresholds:
        sel = np.isfinite(work) & (work >= thr) & (work < 1.0)
        if np.count_nonzero(sel) < 6:
            continue
        x = np.asarray(x2d[sel] - x_ref, dtype=float)
        y = np.asarray(y2d[sel] - y_ref, dtype=float)
        z = np.asarray(work[sel], dtype=float)
        z = np.clip(z, 1e-12, None)
        A = np.column_stack([x * x, x * y, y * y, x, y, np.ones_like(x)])
        b = np.log(z)
        try:
            coeff, *_ = np.linalg.lstsq(A, b, rcond=None)
        except np.linalg.LinAlgError:
            continue
        c_xx, c_xy, c_yy, c_x, c_y, c_0 = coeff
        Q = -np.array([[2.0 * c_xx, c_xy], [c_xy, 2.0 * c_yy]], dtype=float)
        if not np.all(np.isfinite(Q)):
            continue
        try:
            eigvals, eigvecs = np.linalg.eigh(Q)
        except np.linalg.LinAlgError:
            continue
        if np.any(~np.isfinite(eigvals)) or np.any(eigvals <= 0):
            continue
        H = np.array([[2.0 * c_xx, c_xy], [c_xy, 2.0 * c_yy]], dtype=float)
        try:
            center = -np.linalg.solve(H, np.array([c_x, c_y], dtype=float))
        except np.linalg.LinAlgError:
            center = np.array([0.0, 0.0], dtype=float)
        fwhm_axes = 2.0 * np.sqrt(2.0 * np.log(2.0) / eigvals)
        order = np.argsort(fwhm_axes)[::-1]
        fwhm_axes = fwhm_axes[order]
        eigvecs = eigvecs[:, order]
        pa_deg = float(np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0])))
        while pa_deg >= 180.0:
            pa_deg -= 180.0
        while pa_deg < 0.0:
            pa_deg += 180.0
        fit = {
            'bmaj_arcsec': float(fwhm_axes[0]),
            'bmin_arcsec': float(fwhm_axes[1]),
            'bpa_deg': pa_deg,
            'x_center_arcsec': float(x_ref + center[0]),
            'y_center_arcsec': float(y_ref + center[1]),
            'fit_threshold_frac': float(thr),
            'peak_value': img_peak,
        }
        break

    if fit is None:
        return None
    return fit


def _estimate_nominal_effective_beam(
    *,
    beam_fwhm_arcsec: float,
    kernel_name: str,
    kernel_sign: str,
    gwidth_arcsec: float,
    jwidth_arcsec: float,
    r_trunc_arcsec: float,
    eps_u0: float,
    cell_arcsec: float,
    eps_weight_sum: float,
) -> dict[str, float] | None:
    beam = float(beam_fwhm_arcsec)
    support = float(r_trunc_arcsec)
    if not np.isfinite(beam) or beam <= 0 or not np.isfinite(support) or support <= 0:
        return None

    step = min(cell_arcsec / 8.0, beam / 40.0, support / 25.0)
    step = max(step, beam / 200.0, 1e-3)
    radius = max(3.0 * beam, 3.0 * support)
    n = int(np.ceil((2.0 * radius) / step)) + 1
    n = max(129, min(n, 513))
    if n % 2 == 0:
        n += 1
    coords = (np.arange(n, dtype=float) - n // 2) * step
    xx, yy = np.meshgrid(coords, coords)
    rr = np.hypot(xx, yy)

    beam_img = _gaussian_beam_response(rr, beam)
    if kernel_name == 'gjinc':
        kern = _kernel_gjinc(rr, gwidth_arcsec, jwidth_arcsec, eps_u0)
    else:
        kern = _kernel_gauss(rr, gwidth_arcsec)
    kern = np.where(rr <= support, kern, 0.0)
    if kernel_sign == 'positive_only':
        kern = np.where(np.isfinite(kern) & (kern > 0), kern, 0.0)
    else:
        kern = np.where(np.isfinite(kern), kern, 0.0)

    beam_sum = float(np.sum(beam_img, dtype=float))
    kern_sum = float(np.sum(kern, dtype=float))
    if beam_sum <= 0 or (not np.isfinite(kern_sum)) or (abs(kern_sum) <= float(eps_weight_sum)):
        return None
    beam_img = beam_img / beam_sum
    kern = kern / kern_sum
    conv = fftconvolve(beam_img, kern, mode='same')
    if not np.all(np.isfinite(conv)) or np.nanmax(conv) <= 0:
        return None
    conv /= np.nanmax(conv)

    fit = _fit_gaussian_core_2d(conv, xx, yy, min_frac=0.2)
    radial = None
    try:
        radial = _estimate_fwhm_from_profile(np.abs(coords[n // 2:]), conv[n // 2, n // 2:])
    except Exception:
        radial = None
    if fit is None and radial is None:
        return None

    out = {
        'nominal_radial_fwhm_arcsec': float(radial) if radial is not None else np.nan,
    }
    if fit is not None:
        out.update({
            'bmaj_nominal_arcsec': float(fit['bmaj_arcsec']),
            'bmin_nominal_arcsec': float(fit['bmin_arcsec']),
            'bpa_nominal_deg': float(fit['bpa_deg']),
            'beam_nominal_fit_threshold_frac': float(fit['fit_threshold_frac']),
        })
    if radial is not None and ('bmaj_nominal_arcsec' not in out or not np.isfinite(out['bmaj_nominal_arcsec'])):
        out['bmaj_nominal_arcsec'] = float(radial)
        out['bmin_nominal_arcsec'] = float(radial)
        out['bpa_nominal_deg'] = 0.0
    return out


def _estimate_empirical_center_beam(
    *,
    x: np.ndarray,
    y: np.ndarray,
    q: np.ndarray,
    tree: cKDTree,
    config: MapConfig,
    kernel_name: str,
    kernel_sign: str,
    gwidth_arcsec: float,
    jwidth_arcsec: float,
    r_trunc_arcsec: float,
    xg_1d: np.ndarray,
    yg_1d: np.ndarray,
) -> dict[str, float] | None:
    beam = float(config.beam_fwhm_arcsec)
    fit_radius_arcsec = float(getattr(config, 'beam_fit_radius_beam', 1.5)) * beam
    fit_radius_arcsec = max(fit_radius_arcsec, 1.5 * beam)
    min_frac = float(getattr(config, 'beam_fit_min_frac', 0.2))

    x_src = float(xg_1d[len(xg_1d) // 2])
    y_src = float(yg_1d[len(yg_1d) // 2])
    src_r = np.hypot(x - x_src, y - y_src)
    src_amp = _gaussian_beam_response(src_r, beam)

    mx = np.abs(xg_1d - x_src) <= fit_radius_arcsec
    my = np.abs(yg_1d - y_src) <= fit_radius_arcsec
    if not np.any(mx) or not np.any(my):
        return None
    xloc = np.asarray(xg_1d[mx], dtype=float)
    yloc = np.asarray(yg_1d[my], dtype=float)
    xx, yy = np.meshgrid(xloc, yloc)
    pts = np.column_stack((xx.ravel(), yy.ravel()))
    nbrs_list = _query_ball_neighbors(tree, pts, float(r_trunc_arcsec), workers=1, sort_neighbors=True)

    psf = np.full(pts.shape[0], np.nan, dtype=float)
    for i, nbrs_list_i in enumerate(nbrs_list):
        if len(nbrs_list_i) < config.n_min_avg:
            continue
        nbrs = np.asarray(nbrs_list_i, dtype=np.int64)
        dx = x[nbrs] - pts[i, 0]
        dy = y[nbrs] - pts[i, 1]
        rr = np.hypot(dx, dy)
        keep = rr <= float(r_trunc_arcsec)
        if not np.any(keep):
            continue
        nbrs = nbrs[keep]
        rr = rr[keep]
        if nbrs.size < config.n_min_avg:
            continue
        if kernel_name == 'gjinc':
            k = _kernel_gjinc(rr, float(gwidth_arcsec), float(jwidth_arcsec), float(config.eps_u0)).astype(float, copy=False)
        else:
            k = _kernel_gauss(rr, float(gwidth_arcsec)).astype(float, copy=False)
        w = k * q[nbrs]
        m_finite = np.isfinite(w)
        if not np.any(m_finite):
            continue
        if kernel_sign == 'positive_only':
            m_use = m_finite & (w > 0)
        else:
            m_use = m_finite
        if not np.any(m_use):
            continue
        nbrs = nbrs[m_use]
        w = w[m_use]
        if nbrs.size < config.n_min_avg:
            continue
        wsum = np.sum(w, dtype=float)
        if not np.isfinite(wsum) or abs(wsum) <= float(config.eps_weight_sum):
            continue
        psf[i] = np.sum(w * src_amp[nbrs], dtype=float) / wsum

    psf_map = psf.reshape(yy.shape)
    if not np.isfinite(psf_map).any() or np.nanmax(psf_map) <= 0:
        return None
    psf_map = psf_map / np.nanmax(psf_map)
    fit = _fit_gaussian_core_2d(psf_map, xx, yy, min_frac=min_frac)
    radial = None
    try:
        iy0 = np.argmin(np.abs(yloc - y_src))
        ix0 = np.argmin(np.abs(xloc - x_src))
        radial = _estimate_fwhm_from_profile(np.abs(xloc[ix0:] - x_src), psf_map[iy0, ix0:])
    except Exception:
        radial = None
    if fit is None and radial is None:
        return None
    out = {
        'beam_center_x_arcsec': x_src,
        'beam_center_y_arcsec': y_src,
        'beam_empirical_fit_threshold_frac': float(min_frac),
        'empirical_radial_fwhm_arcsec': float(radial) if radial is not None else np.nan,
    }
    if fit is not None:
        out.update({
            'bmaj_empirical_arcsec': float(fit['bmaj_arcsec']),
            'bmin_empirical_arcsec': float(fit['bmin_arcsec']),
            'bpa_empirical_deg': float(fit['bpa_deg']),
            'beam_empirical_peak_x_arcsec': float(fit['x_center_arcsec']),
            'beam_empirical_peak_y_arcsec': float(fit['y_center_arcsec']),
        })
    if radial is not None and ('bmaj_empirical_arcsec' not in out or not np.isfinite(out['bmaj_empirical_arcsec'])):
        out['bmaj_empirical_arcsec'] = float(radial)
        out['bmin_empirical_arcsec'] = float(radial)
        out['bpa_empirical_deg'] = 0.0
    return out

# ==========================================
# 2. Kernels
# ==========================================

def _kernel_gauss(r_arcsec: np.ndarray, gwidth_arcsec: float) -> np.ndarray:
    """Gaussian kernel with HWHM parameter gwidth_arcsec."""
    return np.exp(-np.log(2.0) * (r_arcsec / gwidth_arcsec) ** 2)


def _kernel_gjinc(r_arcsec: np.ndarray, gwidth_arcsec: float, jwidth_arcsec: float, eps_u0: float) -> np.ndarray:
    """GJINC kernel (CASA-like form) with safe handling at r=0."""
    u = np.pi * r_arcsec / jwidth_arcsec
    jinc = np.empty_like(u)

    # r=0 付近の発散を防ぐ安全処理
    m0 = np.abs(u) < eps_u0
    jinc[m0] = 0.5

    mu = ~m0
    if np.any(mu):
        jinc[mu] = j1(u[mu]) / u[mu]

    gauss = np.exp(-np.log(2.0) * (r_arcsec / gwidth_arcsec) ** 2)
    return jinc * gauss


# ==========================================
# 3. Main API & Implementation
# ==========================================

def grid_otf(input_data: GridInput, config: MapConfig) -> GridResult:
    """
    OTF Gridding コアエンジン。指定された空間カーネルと設定に基づいてグリッディングを実行する。
    """
    _validate_input(input_data)
    kres = _validate_and_resolve_config(config)

    return _backend_numpy_avg(input_data, config, kres)


def _query_ball_neighbors(tree, grid_pts, r_trunc_arcsec: float, workers: int, sort_neighbors: bool):
    query_kwargs = {"r": float(r_trunc_arcsec)}
    if workers is not None:
        query_kwargs["workers"] = int(workers)
    if sort_neighbors:
        query_kwargs["return_sorted"] = True

    try:
        return tree.query_ball_point(grid_pts, **query_kwargs)
    except TypeError:
        query_kwargs.pop("return_sorted", None)
        try:
            return tree.query_ball_point(grid_pts, **query_kwargs)
        except TypeError:
            query_kwargs.pop("workers", None)
            return tree.query_ball_point(grid_pts, **query_kwargs)


def _backend_numpy_avg(input_data: GridInput, config: MapConfig, kernel_resolved: dict[str, float | str]) -> GridResult:
    """NumPyベースの加重平均グリッディング実装"""
    runtime = _runtime_options(config)
    dtype_np = np.float64 if runtime["dtype_name"] == "float64" else np.float32

    x_all = np.asarray(input_data.x)
    y_all = np.asarray(input_data.y)
    spec_all = np.asarray(input_data.spec)
    flag_all = np.asarray(input_data.flag)
    time_all = np.asarray(input_data.time)
    ndump, nchan = spec_all.shape

    # --- 1. データの前処理・フィルタリング ---
    valid = np.asarray(flag_all > 0, dtype=bool)
    if config.exclude_turnaround and input_data.is_turnaround is not None:
        valid &= ~_safe_bool_array(input_data.is_turnaround, default=False)

    # 空間・スペクトルの有効性チェック
    valid &= np.isfinite(x_all) & np.isfinite(y_all)
    valid &= np.any(np.isfinite(spec_all), axis=1)

    rms_all, tint_all, tsys_all = None, None, None
    if input_data.rms is not None:
        rms_all = np.asarray(input_data.rms)
        valid &= np.isfinite(rms_all) & (rms_all > 0)
    if input_data.tint is not None:
        tint_all = np.asarray(input_data.tint)
        valid &= np.isfinite(tint_all) & (tint_all >= 0)
    if input_data.tsys is not None:
        tsys_all = np.asarray(input_data.tsys)

    if not np.any(valid):
        raise ValueError("No valid dumps remain after filtering.")

    # フィルタリングの適用
    x = x_all[valid].astype(dtype_np, copy=False)
    y = y_all[valid].astype(dtype_np, copy=False)
    spec = spec_all[valid].astype(dtype_np, copy=False)
    time_vals = time_all[valid].astype(dtype_np, copy=False)
    rms = rms_all[valid].astype(dtype_np, copy=False) if rms_all is not None else None
    tint = tint_all[valid].astype(dtype_np, copy=False) if tint_all is not None else None
    tsys = tsys_all[valid].astype(dtype_np, copy=False) if tsys_all is not None else None

    # 重み (q) の計算
    q = np.ones(len(x), dtype=dtype_np)
    if rms is not None:
        q *= (1.0 / (rms * rms + dtype_np(1e-12))) ** dtype_np(config.alpha_rms)
    if tint is not None and config.beta_tint > 0:
        q *= tint ** dtype_np(config.beta_tint)

    # 重みのクリップ（異常に重い1点が結果を支配するのを防ぐ）
    if config.weight_clip_quantile is not None and q.size > 0:
        qclip = np.quantile(q, config.weight_clip_quantile)
        q = np.minimum(q, dtype_np(qclip))
    if config.weight_clip_max is not None:
        q = np.minimum(q, dtype_np(config.weight_clip_max))

    # --- 2. 幾何学パラメータとグリッドの準備 ---
    kernel_name = str(kernel_resolved["kernel_name"])
    kernel_sign = str(kernel_resolved.get("kernel_sign", "signed"))
    backend_impl = str(kernel_resolved["backend_impl"])
    gwidth_arcsec = dtype_np(float(kernel_resolved["gwidth_pix"]) * config.cell_arcsec)
    jwidth_arcsec = dtype_np(float(kernel_resolved["jwidth_pix"]) * config.cell_arcsec) if kernel_name == "gjinc" else dtype_np(np.nan)
    r_trunc_arcsec = dtype_np(float(kernel_resolved["support_radius_pix"]) * config.cell_arcsec)

    # 出力グリッド座標 [arcsec] の生成
    xg_1d = (config.x0 + np.arange(config.nx, dtype=dtype_np) * dtype_np(config.cell_arcsec)).astype(dtype_np)
    yg_1d = (config.y0 + np.arange(config.ny, dtype=dtype_np) * dtype_np(config.cell_arcsec)).astype(dtype_np)
    xg2d, yg2d = np.meshgrid(xg_1d, yg_1d)
    grid_pts = np.column_stack((xg2d.ravel(), yg2d.ravel()))
    ngrid = grid_pts.shape[0]

    # KDTree で近傍探索
    tree = cKDTree(np.column_stack((x, y)))
    neighbors_list = _query_ball_neighbors(
        tree,
        grid_pts,
        r_trunc_arcsec=float(r_trunc_arcsec),
        workers=int(runtime["workers"]),
        sort_neighbors=bool(runtime["sort_neighbors"]),
    )

    # --- 3. 出力配列の初期化 ---
    if config.fill_nan_for_invalid:
        cube_flat = np.full((ngrid, nchan), np.nan, dtype=dtype_np)
    else:
        cube_flat = np.zeros((ngrid, nchan), dtype=dtype_np)

    weight_flat = np.zeros(ngrid, dtype=dtype_np)
    hit_flat = np.zeros(ngrid, dtype=np.int32)
    mask_flat = np.zeros(ngrid, dtype=bool)
    nsamp_flat = np.zeros(ngrid, dtype=np.int32)
    wsum_raw_flat = np.zeros(ngrid, dtype=dtype_np)
    wabs_raw_flat = np.zeros(ngrid, dtype=dtype_np)
    cancel_flat = np.full(ngrid, np.nan, dtype=dtype_np)
    weight_rel_flat = np.full(ngrid, np.nan, dtype=dtype_np)

    xeff_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_diag_maps else None
    yeff_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_diag_maps else None
    dr_eff_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_diag_maps else None
    neff_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_neff_map else None
    rms_flat = np.full(ngrid, np.nan, dtype=dtype_np) if (config.emit_rms_map and rms is not None) else None
    time_flat = np.full(ngrid, np.nan, dtype=dtype_np) if config.emit_time_map else None
    tint_flat = np.full(ngrid, np.nan, dtype=dtype_np) if (config.emit_tint_map and tint is not None) else None
    tsys_flat = np.full(ngrid, np.nan, dtype=dtype_np) if (config.emit_tsys_map and tsys is not None) else None

    # --- 4. メインループ (グリッド点ごと) ---
    for g_idx, nbrs_list in enumerate(neighbors_list):
        if len(nbrs_list) < config.n_min_avg:
            continue

        nbrs = np.asarray(nbrs_list, dtype=np.int64)
        if runtime["sort_neighbors"]:
            nbrs.sort()

        dx = x[nbrs] - grid_pts[g_idx, 0]
        dy = y[nbrs] - grid_pts[g_idx, 1]
        r = np.sqrt(dx * dx + dy * dy)

        m_keep = r <= r_trunc_arcsec
        if not np.any(m_keep):
            continue
        if not np.all(m_keep):
            nbrs = nbrs[m_keep]
            r = r[m_keep]

        if nbrs.size < config.n_min_avg:
            continue

        if kernel_name == "gjinc":
            k = _kernel_gjinc(r, float(gwidth_arcsec), float(jwidth_arcsec), float(config.eps_u0)).astype(dtype_np, copy=False)
        else:
            k = _kernel_gauss(r, float(gwidth_arcsec)).astype(dtype_np, copy=False)

        w_raw = k * q[nbrs]
        m_finite = np.isfinite(w_raw)
        if not np.any(m_finite):
            continue

        nsamp_total = int(np.count_nonzero(m_finite))
        wsum_raw = np.sum(w_raw[m_finite], dtype=dtype_np)
        wabs_raw = np.sum(np.abs(w_raw[m_finite]), dtype=dtype_np)
        nsamp_flat[g_idx] = nsamp_total
        wsum_raw_flat[g_idx] = wsum_raw
        wabs_raw_flat[g_idx] = wabs_raw
        if np.isfinite(wabs_raw) and (wabs_raw > 0):
            cancel_flat[g_idx] = np.abs(wsum_raw) / wabs_raw

        if kernel_sign == 'positive_only':
            m_use = m_finite & (w_raw > 0)
        else:
            m_use = m_finite

        used_count = int(np.count_nonzero(m_use))
        hit_flat[g_idx] = used_count

        if used_count == 0:
            continue

        if not np.all(m_use):
            nbrs = nbrs[m_use]
            w = w_raw[m_use]
        else:
            w = w_raw[m_use]

        if nbrs.size < config.n_min_avg:
            continue

        W_sum = np.sum(w, dtype=dtype_np)
        if (not np.isfinite(W_sum)) or (np.abs(W_sum) <= dtype_np(config.eps_weight_sum)):
            continue

        weight_flat[g_idx] = W_sum
        mask_flat[g_idx] = bool((used_count >= config.n_min_avg) and (np.abs(W_sum) > dtype_np(config.eps_weight_sum)))
        if not mask_flat[g_idx]:
            continue

        if config.emit_diag_maps:
            xeff = np.sum(w * x[nbrs], dtype=dtype_np) / W_sum
            yeff = np.sum(w * y[nbrs], dtype=dtype_np) / W_sum
            xeff_flat[g_idx] = xeff
            yeff_flat[g_idx] = yeff
            dr_eff_flat[g_idx] = np.sqrt((xeff - grid_pts[g_idx, 0]) ** 2 + (yeff - grid_pts[g_idx, 1]) ** 2) / dtype_np(config.cell_arcsec)

        if config.emit_neff_map:
            w2_sum = np.sum(w * w, dtype=dtype_np)
            if w2_sum > 0:
                neff_flat[g_idx] = (W_sum * W_sum) / w2_sum

        if rms_flat is not None:
            rms_flat[g_idx] = np.sqrt(np.sum((w * w) * (rms[nbrs] * rms[nbrs]), dtype=dtype_np) / (W_sum * W_sum))

        if time_flat is not None:
            time_flat[g_idx] = np.sum(w * time_vals[nbrs], dtype=dtype_np) / W_sum

        if tint_flat is not None:
            if neff_flat is not None:
                neff_val = neff_flat[g_idx]
            else:
                w2_sum = np.sum(w * w, dtype=dtype_np)
                neff_val = (W_sum * W_sum) / w2_sum if w2_sum > 0 else 1.0
            tint_flat[g_idx] = (np.sum(w * tint[nbrs], dtype=dtype_np) / W_sum) * neff_val

        if tsys_flat is not None:
            w_tsys = w * tsys[nbrs]
            m_tsys = np.isfinite(w_tsys)
            if np.any(m_tsys):
                denom = np.sum(w[m_tsys], dtype=dtype_np)
                if np.isfinite(denom) and (np.abs(denom) > dtype_np(config.eps_weight_sum)):
                    tsys_flat[g_idx] = np.sum(w_tsys[m_tsys], dtype=dtype_np) / denom

        spec_nbrs = spec[nbrs, :]
        for ch0 in range(0, nchan, config.chunk_ch):
            ch1 = builtins.min(ch0 + config.chunk_ch, nchan)
            spec_chunk = spec_nbrs[:, ch0:ch1]
            valid_data = np.isfinite(spec_chunk)
            w_2d = w[:, None] * valid_data
            W_sum_chan = np.sum(w_2d, axis=0, dtype=dtype_np)
            spec_sum_chan = np.sum(w_2d * np.where(valid_data, spec_chunk, 0.0), axis=0, dtype=dtype_np)
            m_chan = np.abs(W_sum_chan) > dtype_np(config.eps_weight_sum)
            cube_view = cube_flat[g_idx, ch0:ch1]
            cube_view[m_chan] = spec_sum_chan[m_chan] / W_sum_chan[m_chan]
            cube_flat[g_idx, ch0:ch1] = cube_view

    # --- 4. post-gridding 妥当性判定（相対 weight / cancellation） ---
    provisional_valid = mask_flat.copy()
    median_abs_weight = np.nan
    valid_abs_weight = np.abs(weight_flat[provisional_valid])
    valid_abs_weight = valid_abs_weight[np.isfinite(valid_abs_weight)]
    if valid_abs_weight.size > 0:
        median_abs_weight = float(np.nanmedian(valid_abs_weight))
        if np.isfinite(median_abs_weight) and median_abs_weight > 0:
            weight_rel_flat[:] = np.abs(weight_flat) / dtype_np(median_abs_weight)

    low_weight_invalid_count = 0
    low_cancel_invalid_count = 0

    post_valid = provisional_valid.copy()
    min_abs_weight_ratio = float(getattr(config, 'min_abs_weight_ratio', 0.0))
    if min_abs_weight_ratio > 0 and np.isfinite(median_abs_weight) and median_abs_weight > 0:
        low_weight = provisional_valid & np.isfinite(weight_rel_flat) & (weight_rel_flat < dtype_np(min_abs_weight_ratio))
        low_weight_invalid_count = int(np.count_nonzero(low_weight))
        post_valid[low_weight] = False

    min_cancel_ratio = float(getattr(config, 'min_cancel_ratio', 0.0))
    if min_cancel_ratio > 0:
        low_cancel = provisional_valid & np.isfinite(cancel_flat) & (cancel_flat < dtype_np(min_cancel_ratio))
        low_cancel_invalid_count = int(np.count_nonzero(low_cancel & post_valid))
        post_valid[low_cancel] = False

    if config.fill_nan_for_invalid and np.any(~post_valid):
        cube_flat[~post_valid, :] = np.nan
    mask_flat = post_valid

    if (min_abs_weight_ratio > 0) and (low_weight_invalid_count > 0):
        warnings.warn(
            f"{low_weight_invalid_count} pixels were masked by min_abs_weight_ratio={min_abs_weight_ratio:g}.",
            RuntimeWarning,
            stacklevel=2,
        )
    if (min_cancel_ratio > 0) and (low_cancel_invalid_count > 0):
        warnings.warn(
            f"{low_cancel_invalid_count} pixels were masked by min_cancel_ratio={min_cancel_ratio:g}.",
            RuntimeWarning,
            stacklevel=2,
        )

    ny, nx = config.ny, config.nx

    n_warn = None
    if dr_eff_flat is not None:
        n_warn = int(np.sum(np.isfinite(dr_eff_flat) & (dr_eff_flat > dtype_np(config.dr_eff_warn_pix))))

    beam_meta = {}
    if getattr(config, 'estimate_effective_beam', True):
        nominal_beam = _estimate_nominal_effective_beam(
            beam_fwhm_arcsec=float(config.beam_fwhm_arcsec),
            kernel_name=kernel_name,
            kernel_sign=kernel_sign,
            gwidth_arcsec=float(gwidth_arcsec),
            jwidth_arcsec=float(jwidth_arcsec),
            r_trunc_arcsec=float(r_trunc_arcsec),
            eps_u0=float(config.eps_u0),
            cell_arcsec=float(config.cell_arcsec),
            eps_weight_sum=float(config.eps_weight_sum),
        )
        if nominal_beam is not None:
            beam_meta.update(nominal_beam)

        empirical_beam = _estimate_empirical_center_beam(
            x=x,
            y=y,
            q=q,
            tree=tree,
            config=config,
            kernel_name=kernel_name,
            kernel_sign=kernel_sign,
            gwidth_arcsec=float(gwidth_arcsec),
            jwidth_arcsec=float(jwidth_arcsec),
            r_trunc_arcsec=float(r_trunc_arcsec),
            xg_1d=xg_1d,
            yg_1d=yg_1d,
        )
        if empirical_beam is not None:
            beam_meta.update(empirical_beam)

    bmaj_eff_arcsec = float(
        beam_meta.get('bmaj_empirical_arcsec', beam_meta.get('bmaj_nominal_arcsec', np.nan))
    )

    return GridResult(
        cube=cube_flat.reshape(ny, nx, nchan),
        weight_map=weight_flat.reshape(ny, nx),
        hit_map=hit_flat.reshape(ny, nx),
        mask_map=mask_flat.reshape(ny, nx),
        nsamp_map=nsamp_flat.reshape(ny, nx),
        wsum_map=wsum_raw_flat.reshape(ny, nx),
        wabs_map=wabs_raw_flat.reshape(ny, nx),
        cancel_map=cancel_flat.reshape(ny, nx),
        weight_rel_map=weight_rel_flat.reshape(ny, nx),
        xeff_map=xeff_flat.reshape(ny, nx) if xeff_flat is not None else None,
        yeff_map=yeff_flat.reshape(ny, nx) if yeff_flat is not None else None,
        dr_eff_map_pix=dr_eff_flat.reshape(ny, nx) if dr_eff_flat is not None else None,
        neff_map=neff_flat.reshape(ny, nx) if neff_flat is not None else None,
        rms_map=rms_flat.reshape(ny, nx) if rms_flat is not None else None,
        time_map=time_flat.reshape(ny, nx) if time_flat is not None else None,
        tint_map=tint_flat.reshape(ny, nx) if tint_flat is not None else None,
        tsys_map=tsys_flat.reshape(ny, nx) if tsys_flat is not None else None,
        meta={
            "kernel": kernel_name,
            "kernel_preset": str(kernel_resolved.get("kernel_preset", getattr(config, "kernel_preset", "mangum2007"))),
            "kernel_sign": str(kernel_sign),
            "backend": backend_impl,
            "cell_arcsec": float(config.cell_arcsec),
            "beam_fwhm_arcsec": float(config.beam_fwhm_arcsec),
            "beam_pix": float(kernel_resolved.get("beam_pix", config.beam_fwhm_arcsec / config.cell_arcsec)),
            "cell_over_beam": float(kernel_resolved.get("cell_over_beam", config.cell_arcsec / config.beam_fwhm_arcsec)),
            "cell_is_coarse": bool(kernel_resolved.get("cell_is_coarse", False)),
            "bmaj_eff_arcsec": bmaj_eff_arcsec,
            "gwidth_arcsec": float(gwidth_arcsec),
            "jwidth_arcsec": float(jwidth_arcsec) if kernel_name == "gjinc" else np.nan,
            "support_radius_arcsec": float(r_trunc_arcsec),
            "weight_map_semantics": "gridding_denominator_sum_after_sign_filter",
            "hit_map_semantics": "samples_used_in_gridding_denominator_after_sign_filter",
            "nsamp_map_semantics": "finite_support_sample_count_before_sign_filter",
            "wsum_map_semantics": "signed_sum_all_finite_support_weights_before_sign_filter",
            "wabs_map_semantics": "absolute_sum_all_finite_support_weights_before_sign_filter",
            "cancel_map_semantics": "absolute(WSUM) / WABS before sign filter",
            "weight_rel_map_semantics": "absolute(WEIGHT) / median_absolute_WEIGHT_over_provisional_valid_pixels",
            "gwidth_source": str(kernel_resolved.get("gwidth_source", "")),
            "jwidth_source": str(kernel_resolved.get("jwidth_source", "")),
            "support_source": str(kernel_resolved.get("support_source", "")),
            "mask_map_semantics": "spatial_denominator_valid_after_sign_filter_before_channel_nan_mask",
            "alpha_rms": float(config.alpha_rms),
            "beta_tint": float(config.beta_tint),
            "support_radius_pix": float(kernel_resolved["support_radius_pix"]),
            "median_abs_weight": median_abs_weight,
            "min_abs_weight_ratio": float(getattr(config, 'min_abs_weight_ratio', 0.0)),
            "min_cancel_ratio": float(getattr(config, 'min_cancel_ratio', 0.0)),
            "low_weight_invalid_count": int(low_weight_invalid_count),
            "low_cancel_invalid_count": int(low_cancel_invalid_count),
            "dr_eff_warn_count": n_warn,
            "reproducible_mode": bool(runtime["reproducible"]),
            "gridding_workers": int(runtime["workers"]),
            "sorted_neighbors": bool(runtime["sort_neighbors"]),
            **beam_meta,
        },
    )
