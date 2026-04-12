# src/sd_radio_spectral_fits/otf/core.py
from __future__ import annotations

import builtins
import warnings
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import j1, pro_ang1
from scipy.signal import fftconvolve

from .config import MapConfig, GridInput, GridResult, normalize_row_flag_mask


SQRT_LN2 = float(np.sqrt(np.log(2.0)))
FIRST_JINC_NULL_OVER_PI = 3.8317059702075125 / np.pi

# Public presets are pixel-based.
#
# GJINC (Mangum/Sawada-like):
#   jwidth = 1.55 pixel
#   gwidth = 2.52 * sqrt(log(2)) pixel
#   support radius = 3 pixel
#
# GJINC (CASA-like):
#   jwidth = 1.55 pixel
#   gwidth = 2.52 * sqrt(log(2)) pixel
#   support radius = first null of the resolved jinc width
#
# SF (casacore Sph_Conv-like):
#   convsupport = cutoff radius in pixel
#   support radius = convsupport * cell
#
# GAUSS:
#   gwidth = sqrt(log(2)) pixel (HWHM)
#   support radius = 3 * gwidth
DEFAULT_GJINC_JWIDTH_PIX = 1.55
DEFAULT_GJINC_GWIDTH_PIX = 2.52 * SQRT_LN2
DEFAULT_GAUSS_GWIDTH_PIX = SQRT_LN2
DEFAULT_SF_CONVSUPPORT = 3
DEFAULT_SF_BEAM_MATCH_CONVSUPPORT = 3
SF_BEAM_MATCH_SHAPE_C = {
    3: 10.685961961746216,
    4: 13.737787961959839,
    5: 13.719022989273071,
    6: 13.70369029045105,
}
DEFAULT_GJINC_BEAM_JWIDTH_BEAM = DEFAULT_GJINC_JWIDTH_PIX / 3.0
DEFAULT_GJINC_BEAM_GWIDTH_BEAM = DEFAULT_GJINC_GWIDTH_PIX / 3.0
DEFAULT_GJINC_BEAM_SUPPORT_BEAM_MANGUM = 1.0
DEFAULT_GJINC_BEAM_SUPPORT_BEAM_CASA = FIRST_JINC_NULL_OVER_PI * DEFAULT_GJINC_BEAM_JWIDTH_BEAM
RECOMMENDED_MAX_CELL_OVER_BEAM = 1.0 / 3.0
CELL_OVER_BEAM_TOL = 1.0e-12


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
    time = np.asarray(input_data.time)

    if spec.ndim != 2:
        raise ValueError(f"spec must be 2D (ndump, nchan), got shape={spec.shape}")
    ndump = spec.shape[0]
    flag = normalize_row_flag_mask(input_data.flag, ndump=ndump, allow_none=True, name='flag')

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
    elif backend == "numpy_sf":
        backend = "numpy"
        if kernel not in ("sf", "sf_beam"):
            raise ValueError("backend='numpy_sf' requires kernel='sf' or 'sf_beam'")
    elif backend == "numpy_gjinc_beam":
        backend = "numpy"
        if kernel != "gjinc_beam":
            raise ValueError("backend='numpy_gjinc_beam' requires kernel='gjinc_beam'")
    elif backend == "numpy_gauss":
        backend = "numpy"
        if kernel != "gauss":
            raise ValueError("backend='numpy_gauss' requires kernel='gauss'")
    elif backend == "cygrid_gauss":
        backend = "cygrid"
        if kernel != "gauss":
            raise ValueError("backend='cygrid_gauss' requires kernel='gauss'")

    if kernel not in ("sf", "sf_beam", "gjinc", "gjinc_beam", "gauss"):
        raise ValueError(f"Unknown kernel: {config.kernel!r}. Use 'sf', 'sf_beam', 'gjinc', 'gjinc_beam', or 'gauss'.")

    return backend, kernel



def _parse_exact_positive_int(value, *, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer >= 1, not bool")
    if isinstance(value, (int, np.integer)):
        iv = int(value)
    elif isinstance(value, (float, np.floating)):
        fv = float(value)
        if not np.isfinite(fv) or not fv.is_integer():
            raise ValueError(f"{name} must be an integer >= 1")
        iv = int(fv)
    elif isinstance(value, str):
        s = value.strip()
        if not s or s.startswith('+'):
            s = s.lstrip('+')
        if not s.isdigit():
            raise ValueError(f"{name} must be an integer >= 1")
        iv = int(s)
    else:
        raise ValueError(f"{name} must be an integer >= 1")
    if iv <= 0:
        raise ValueError(f"{name} must be >= 1")
    return iv


def _validate_and_resolve_config(config: MapConfig) -> dict[str, object]:
    backend_impl, kernel_name = _normalize_backend_and_kernel(config)
    runtime = _runtime_options(config)

    if config.estimator == "plane":
        raise NotImplementedError("estimator='plane' is not implemented yet. Please use 'avg'.")

    if backend_impl != "numpy":
        raise NotImplementedError(
            f"backend={backend_impl!r} is not implemented in this build. "
            "Please use backend='numpy' or one of its aliases."
        )

    if config.nx <= 0 or config.ny <= 0:
        raise ValueError("nx and ny must be positive")
    if getattr(config, 'beam_fwhm_arcsec', None) is None or float(config.beam_fwhm_arcsec) <= 0:
        raise ValueError("beam_fwhm_arcsec must be positive")
    if getattr(config, 'cell_arcsec', None) is None:
        config.cell_arcsec = float(config.beam_fwhm_arcsec) / 3.0
    if float(config.cell_arcsec) <= 0:
        raise ValueError("cell_arcsec must be positive")
    if config.chunk_ch <= 0:
        raise ValueError("chunk_ch must be positive")
    if runtime["dtype_name"] not in ("float32", "float64"):
        raise ValueError("dtype must be 'float32' or 'float64'")
    weight_mode = str(getattr(config, 'weight_mode', 'uniform')).strip().lower()
    if weight_mode not in ('uniform', 'rms'):
        raise ValueError("weight_mode must be 'uniform' or 'rms'")
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
    cell_is_coarse = bool(cell_over_beam > (RECOMMENDED_MAX_CELL_OVER_BEAM * (1.0 + CELL_OVER_BEAM_TOL)))

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

    def _normalize_kernel_preset(value, *, kernel_name: str) -> str | None:
        raw = None if value is None else str(value).strip().lower()
        if kernel_name in ('gjinc', 'gjinc_beam'):
            if raw in (None, '', 'none'):
                return 'mangum'
            if raw in ('mangum', 'mangum2007'):
                return 'mangum'
            if raw in ('casa', 'legacy'):
                return 'casa'
            raise ValueError(
                "kernel_preset for kernel='gjinc'/'gjinc_beam' must be one of: None, 'mangum', 'mangum2007', 'casa', 'legacy'"
            )
        if raw not in (None, '', 'none', 'default'):
            warnings.warn(
                f"kernel_preset={value!r} is ignored for kernel={kernel_name!r}.",
                RuntimeWarning,
                stacklevel=2,
            )
        return None

    kernel_preset = _normalize_kernel_preset(getattr(config, 'kernel_preset', None), kernel_name=kernel_name)

    kernel_sign = str(getattr(config, 'kernel_sign', 'auto')).lower()
    if kernel_sign not in ('auto', 'signed', 'positive_only'):
        raise ValueError("kernel_sign must be 'auto', 'signed', or 'positive_only'")
    if kernel_sign == 'auto':
        if kernel_name in ('gjinc', 'gjinc_beam') and kernel_preset == 'mangum':
            kernel_sign = 'signed'
        else:
            kernel_sign = 'positive_only'
    if kernel_name in ('sf', 'sf_beam', 'gauss') and kernel_sign != 'positive_only':
        warnings.warn(
            f"kernel={kernel_name!r} only supports positive weights; overriding kernel_sign={kernel_sign!r} -> 'positive_only'.",
            RuntimeWarning,
            stacklevel=2,
        )
        kernel_sign = 'positive_only'
    if kernel_name in ('gjinc', 'gjinc_beam') and kernel_preset == 'mangum' and kernel_sign != 'signed':
        warnings.warn(
            "kernel_preset='mangum' is literature-oriented and is normally used with kernel_sign='signed'. "
            "You explicitly requested a non-standard combination.",
            RuntimeWarning,
            stacklevel=2,
        )
    if kernel_name in ('gjinc', 'gjinc_beam') and kernel_preset == 'casa' and kernel_sign != 'positive_only':
        warnings.warn(
            "kernel_preset='casa' is normally used with kernel_sign='positive_only'. "
            "You explicitly requested a non-standard combination.",
            RuntimeWarning,
            stacklevel=2,
        )

    convsupport = None
    sf_beam_match_convsupport = None
    sf_beam_support_beam = np.nan
    sf_beam_shape_c = np.nan
    gjinc_beam_support_beam = np.nan
    gjinc_beam_gwidth_beam = np.nan
    gjinc_beam_jwidth_beam = np.nan
    if kernel_name == 'sf':
        convsupport = _parse_exact_positive_int(getattr(config, 'convsupport', DEFAULT_SF_CONVSUPPORT), name='convsupport')
        if convsupport < 3:
            warnings.warn(
                f"convsupport={convsupport} is smaller than the formal CASA/casacore default 3 for kernel='sf'. "
                "This may increase aliasing and is not recommended for routine OTF mapping.",
                RuntimeWarning,
                stacklevel=2,
            )
        if any(getattr(config, name) is not None for name in ('support_radius_pix', 'support_radius_beam', 'truncate')):
            raise ValueError(
                "kernel='sf' uses convsupport only. Do not set support_radius_pix, support_radius_beam, or truncate."
            )
        if any(getattr(config, name) is not None for name in ('gwidth_pix', 'gwidth_beam', 'jwidth_pix', 'jwidth_beam')):
            warnings.warn(
                "kernel='sf' ignores gwidth_*/jwidth_* parameters; they will be ignored.",
                RuntimeWarning,
                stacklevel=2,
            )

    def _resolve_support_radius(*, default_pix: float, default_source: str, allow_first_null: bool, jwidth_pix: float | None) -> tuple[float, str]:
        if config.support_radius_pix is not None:
            support = float(config.support_radius_pix)
            source = 'support_radius_pix'
        elif getattr(config, 'support_radius_beam', None) is not None:
            support = float(config.support_radius_beam) * beam_pix
            source = 'support_radius_beam'
        elif config.truncate is None:
            support = float(default_pix)
            source = default_source
        elif isinstance(config.truncate, (int, float)) and not isinstance(config.truncate, (bool, np.bool_)):
            support = float(config.truncate)
            source = 'truncate_numeric'
        elif config.truncate == 'first_null':
            if not allow_first_null or jwidth_pix is None or not np.isfinite(jwidth_pix):
                raise ValueError("truncate='first_null' is only valid for kernel='gjinc'")
            support = float(FIRST_JINC_NULL_OVER_PI * jwidth_pix)
            source = 'truncate_first_null'
        else:
            raise ValueError('Invalid truncate mode.')
        if support <= 0:
            raise ValueError('support_radius must be positive')
        return support, source

    if kernel_name == 'sf':
        gwidth_pix = np.nan
        gwidth_source = 'not_used'
        jwidth_pix = np.nan
        jwidth_source = 'not_used'
        support_radius_pix = float(convsupport)
        support_source = 'convsupport'
    elif kernel_name == 'sf_beam':
        if any(getattr(config, name) is not None for name in ('support_radius_pix', 'support_radius_beam', 'truncate')):
            raise ValueError(
                "kernel='sf_beam' uses sf_beam_support_beam / sf_beam_shape_c (or sf_beam_match_convsupport). "
                "Do not set support_radius_pix, support_radius_beam, or truncate."
            )
        if any(getattr(config, name) is not None for name in ('gwidth_pix', 'gwidth_beam', 'jwidth_pix', 'jwidth_beam')):
            warnings.warn(
                "kernel='sf_beam' ignores gwidth_*/jwidth_* parameters; they will be ignored.",
                RuntimeWarning,
                stacklevel=2,
            )
        if getattr(config, 'convsupport', DEFAULT_SF_CONVSUPPORT) != DEFAULT_SF_CONVSUPPORT:
            warnings.warn(
                "kernel='sf_beam' ignores convsupport directly. Use sf_beam_match_convsupport or explicit sf_beam_* parameters.",
                RuntimeWarning,
                stacklevel=2,
            )
        raw_match = getattr(config, 'sf_beam_match_convsupport', DEFAULT_SF_BEAM_MATCH_CONVSUPPORT)
        if raw_match is not None:
            sf_beam_match_convsupport = _parse_exact_positive_int(raw_match, name='sf_beam_match_convsupport')
            if sf_beam_match_convsupport not in SF_BEAM_MATCH_SHAPE_C:
                raise ValueError("sf_beam_match_convsupport must be one of 3, 4, 5, 6")
        raw_support_beam = getattr(config, 'sf_beam_support_beam', None)
        if raw_support_beam is None:
            if sf_beam_match_convsupport is None:
                raise ValueError("kernel='sf_beam' requires sf_beam_support_beam or sf_beam_match_convsupport")
            sf_beam_support_beam = float(sf_beam_match_convsupport) / 3.0
            support_source = 'sf_beam_match_convsupport'
        else:
            sf_beam_support_beam = float(raw_support_beam)
            support_source = 'sf_beam_support_beam'
        if (not np.isfinite(sf_beam_support_beam)) or sf_beam_support_beam <= 0:
            raise ValueError('sf_beam_support_beam must be positive')
        raw_shape_c = getattr(config, 'sf_beam_shape_c', None)
        if raw_shape_c is None:
            if sf_beam_match_convsupport is None:
                raise ValueError("kernel='sf_beam' requires sf_beam_shape_c or sf_beam_match_convsupport")
            sf_beam_shape_c = float(SF_BEAM_MATCH_SHAPE_C[int(sf_beam_match_convsupport)])
        else:
            sf_beam_shape_c = float(raw_shape_c)
        if (not np.isfinite(sf_beam_shape_c)) or sf_beam_shape_c <= 0:
            raise ValueError('sf_beam_shape_c must be positive')
        gwidth_pix = np.nan
        gwidth_source = 'not_used'
        jwidth_pix = np.nan
        jwidth_source = 'not_used'
        support_radius_pix = float(sf_beam_support_beam * beam_pix)
    elif kernel_name == 'gjinc_beam':
        if any(getattr(config, name) is not None for name in ('support_radius_pix', 'support_radius_beam', 'truncate')):
            raise ValueError(
                "kernel='gjinc_beam' uses gjinc_beam_* parameters and kernel_preset only. "
                "Do not set support_radius_pix, support_radius_beam, or truncate."
            )
        if any(getattr(config, name) is not None for name in ('gwidth_pix', 'gwidth_beam', 'jwidth_pix', 'jwidth_beam')):
            warnings.warn(
                "kernel='gjinc_beam' ignores generic gwidth_*/jwidth_* parameters; use gjinc_beam_* parameters instead.",
                RuntimeWarning,
                stacklevel=2,
            )
        if getattr(config, 'convsupport', DEFAULT_SF_CONVSUPPORT) != DEFAULT_SF_CONVSUPPORT:
            warnings.warn(
                "kernel='gjinc_beam' ignores convsupport directly. Use kernel_preset and/or gjinc_beam_* parameters.",
                RuntimeWarning,
                stacklevel=2,
            )
        raw_jwidth_beam = getattr(config, 'gjinc_beam_jwidth_beam', None)
        if raw_jwidth_beam is None:
            gjinc_beam_jwidth_beam = float(DEFAULT_GJINC_BEAM_JWIDTH_BEAM)
            jwidth_source = 'gjinc_beam_preset_default'
        else:
            gjinc_beam_jwidth_beam = float(raw_jwidth_beam)
            jwidth_source = 'gjinc_beam_jwidth_beam'
        if (not np.isfinite(gjinc_beam_jwidth_beam)) or gjinc_beam_jwidth_beam <= 0:
            raise ValueError('gjinc_beam_jwidth_beam must be positive')
        raw_gwidth_beam = getattr(config, 'gjinc_beam_gwidth_beam', None)
        if raw_gwidth_beam is None:
            gjinc_beam_gwidth_beam = float(DEFAULT_GJINC_BEAM_GWIDTH_BEAM)
            gwidth_source = 'gjinc_beam_preset_default'
        else:
            gjinc_beam_gwidth_beam = float(raw_gwidth_beam)
            gwidth_source = 'gjinc_beam_gwidth_beam'
        if (not np.isfinite(gjinc_beam_gwidth_beam)) or gjinc_beam_gwidth_beam <= 0:
            raise ValueError('gjinc_beam_gwidth_beam must be positive')
        raw_support_beam = getattr(config, 'gjinc_beam_support_beam', None)
        if raw_support_beam is None:
            if kernel_preset == 'mangum':
                gjinc_beam_support_beam = float(DEFAULT_GJINC_BEAM_SUPPORT_BEAM_MANGUM)
                support_source = 'gjinc_beam_preset_mangum'
            elif kernel_preset == 'casa':
                gjinc_beam_support_beam = float(FIRST_JINC_NULL_OVER_PI * gjinc_beam_jwidth_beam)
                support_source = 'gjinc_beam_preset_casa_first_null'
            else:
                raise ValueError(f"Unsupported kernel_preset={kernel_preset!r} for kernel='gjinc_beam'")
        else:
            gjinc_beam_support_beam = float(raw_support_beam)
            support_source = 'gjinc_beam_support_beam'
        if (not np.isfinite(gjinc_beam_support_beam)) or gjinc_beam_support_beam <= 0:
            raise ValueError('gjinc_beam_support_beam must be positive')
        jwidth_pix = float(gjinc_beam_jwidth_beam * beam_pix)
        gwidth_pix = float(gjinc_beam_gwidth_beam * beam_pix)
        support_radius_pix = float(gjinc_beam_support_beam * beam_pix)
    else:
        if kernel_name == 'gjinc':
            default_gwidth_pix = DEFAULT_GJINC_GWIDTH_PIX
            default_jwidth_pix = DEFAULT_GJINC_JWIDTH_PIX
        else:
            default_gwidth_pix = DEFAULT_GAUSS_GWIDTH_PIX
            default_jwidth_pix = np.nan

        if config.gwidth_pix is not None:
            gwidth_pix = float(config.gwidth_pix)
            gwidth_source = 'gwidth_pix'
        elif config.gwidth_beam is not None:
            gwidth_pix = float(config.gwidth_beam) * beam_pix
            gwidth_source = 'gwidth_beam'
        else:
            gwidth_pix = float(default_gwidth_pix)
            gwidth_source = 'preset_pix'
        if gwidth_pix <= 0:
            raise ValueError('gwidth must be positive')

        if kernel_name == 'gjinc':
            if config.jwidth_pix is not None:
                jwidth_pix = float(config.jwidth_pix)
                jwidth_source = 'jwidth_pix'
            elif config.jwidth_beam is not None:
                jwidth_pix = float(config.jwidth_beam) * beam_pix
                jwidth_source = 'jwidth_beam'
            else:
                jwidth_pix = float(default_jwidth_pix)
                jwidth_source = 'preset_pix'
            if jwidth_pix <= 0:
                raise ValueError("jwidth must be positive for kernel='gjinc'")
        else:
            jwidth_pix = np.nan
            jwidth_source = 'not_used'
            if any(getattr(config, name) is not None for name in ('jwidth_pix', 'jwidth_beam')):
                warnings.warn(
                    "kernel='gauss' ignores jwidth_*/truncate='first_null' semantics for jinc kernels.",
                    RuntimeWarning,
                    stacklevel=2,
                )

        if kernel_name == 'gjinc':
            if kernel_preset == 'mangum':
                default_support_pix = 3.0
                default_support_source = 'preset_3pix'
            else:
                default_support_pix = float(FIRST_JINC_NULL_OVER_PI * jwidth_pix)
                default_support_source = 'preset_first_null'
        else:
            default_support_pix = float(3.0 * gwidth_pix)
            default_support_source = 'preset_3hwhm'

        support_radius_pix, support_source = _resolve_support_radius(
            default_pix=float(default_support_pix),
            default_source=default_support_source,
            allow_first_null=(kernel_name == 'gjinc'),
            jwidth_pix=None if not np.isfinite(jwidth_pix) else float(jwidth_pix),
        )

    return {
        'backend_impl': backend_impl,
        'kernel_name': kernel_name,
        'kernel_preset': kernel_preset,
        'kernel_sign': kernel_sign,
        'convsupport': convsupport,
        'gwidth_pix': float(gwidth_pix),
        'jwidth_pix': float(jwidth_pix),
        'support_radius_pix': float(support_radius_pix),
        'sf_beam_match_convsupport': sf_beam_match_convsupport,
        'sf_beam_support_beam': float(sf_beam_support_beam),
        'sf_beam_shape_c': float(sf_beam_shape_c),
        'gjinc_beam_support_beam': float(gjinc_beam_support_beam),
        'gjinc_beam_gwidth_beam': float(gjinc_beam_gwidth_beam),
        'gjinc_beam_jwidth_beam': float(gjinc_beam_jwidth_beam),
        'beam_pix': float(beam_pix),
        'cell_over_beam': float(cell_over_beam),
        'cell_is_coarse': bool(cell_is_coarse),
        'gwidth_source': gwidth_source,
        'jwidth_source': jwidth_source,
        'support_source': support_source,
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
    sf_beam_support_beam: float,
    sf_beam_shape_c: float,
    sf_convsupport_exact: float | None,
    gjinc_beam_support_beam: float,
    gjinc_beam_gwidth_beam: float,
    gjinc_beam_jwidth_beam: float,
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
    kern = _evaluate_kernel(
        rr,
        kernel_name=kernel_name,
        gwidth_arcsec=gwidth_arcsec,
        jwidth_arcsec=jwidth_arcsec,
        support_radius_pix=(float(sf_convsupport_exact) if kernel_name == 'sf' and sf_convsupport_exact is not None and np.isfinite(float(sf_convsupport_exact)) else float(support / cell_arcsec)),
        cell_arcsec=float(cell_arcsec),
        beam_fwhm_arcsec=float(beam),
        sf_beam_support_beam=float(sf_beam_support_beam),
        sf_beam_shape_c=float(sf_beam_shape_c),
        gjinc_beam_support_beam=float(gjinc_beam_support_beam),
        gjinc_beam_gwidth_beam=float(gjinc_beam_gwidth_beam),
        gjinc_beam_jwidth_beam=float(gjinc_beam_jwidth_beam),
        eps_u0=eps_u0,
    )
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
            'bmaj_nominal_gaussian_fit_arcsec': float(fit['bmaj_arcsec']),
            'bmin_nominal_gaussian_fit_arcsec': float(fit['bmin_arcsec']),
            'bpa_nominal_gaussian_fit_deg': float(fit['bpa_deg']),
            'beam_nominal_fit_threshold_frac': float(fit['fit_threshold_frac']),
        })
    if radial is not None and kernel_name in ('sf', 'sf_beam'):
        out['bmaj_nominal_arcsec'] = float(radial)
        out['bmin_nominal_arcsec'] = float(radial)
        out['bpa_nominal_deg'] = 0.0
    elif fit is not None:
        out['bmaj_nominal_arcsec'] = float(fit['bmaj_arcsec'])
        out['bmin_nominal_arcsec'] = float(fit['bmin_arcsec'])
        out['bpa_nominal_deg'] = float(fit['bpa_deg'])
    elif radial is not None:
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
    sf_beam_support_beam: float,
    sf_beam_shape_c: float,
    sf_convsupport_exact: float | None,
    gjinc_beam_support_beam: float,
    gjinc_beam_gwidth_beam: float,
    gjinc_beam_jwidth_beam: float,
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
        k = _evaluate_kernel(
            rr,
            kernel_name=kernel_name,
            gwidth_arcsec=float(gwidth_arcsec),
            jwidth_arcsec=float(jwidth_arcsec),
            support_radius_pix=(float(sf_convsupport_exact) if kernel_name == 'sf' and sf_convsupport_exact is not None and np.isfinite(float(sf_convsupport_exact)) else float(r_trunc_arcsec / config.cell_arcsec)),
            cell_arcsec=float(config.cell_arcsec),
            beam_fwhm_arcsec=float(config.beam_fwhm_arcsec),
            sf_beam_support_beam=float(sf_beam_support_beam),
            sf_beam_shape_c=float(sf_beam_shape_c),
            gjinc_beam_support_beam=float(gjinc_beam_support_beam),
            gjinc_beam_gwidth_beam=float(gjinc_beam_gwidth_beam),
            gjinc_beam_jwidth_beam=float(gjinc_beam_jwidth_beam),
            eps_u0=float(config.eps_u0),
        ).astype(float, copy=False)
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
            'bmaj_empirical_gaussian_fit_arcsec': float(fit['bmaj_arcsec']),
            'bmin_empirical_gaussian_fit_arcsec': float(fit['bmin_arcsec']),
            'bpa_empirical_gaussian_fit_deg': float(fit['bpa_deg']),
            'beam_empirical_peak_x_arcsec': float(fit['x_center_arcsec']),
            'beam_empirical_peak_y_arcsec': float(fit['y_center_arcsec']),
        })
    if radial is not None and kernel_name in ('sf', 'sf_beam'):
        out['bmaj_empirical_arcsec'] = float(radial)
        out['bmin_empirical_arcsec'] = float(radial)
        out['bpa_empirical_deg'] = 0.0
    elif fit is not None:
        out['bmaj_empirical_arcsec'] = float(fit['bmaj_arcsec'])
        out['bmin_empirical_arcsec'] = float(fit['bmin_arcsec'])
        out['bpa_empirical_deg'] = float(fit['bpa_deg'])
    elif radial is not None:
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


# Fred Schwab / casacore spheroidal rational-approximation tables.
# Shapes follow the original Fortran ordering.
_SF_ALPHA = np.asarray([0.0, 0.5, 1.0, 1.5, 2.0], dtype=float)
_SF_P4 = np.asarray([
    [0.01584774, -0.1269612, 0.2333851, -0.1636744, 0.05014648],
    [0.03101855, -0.1641253, 0.23855, -0.1417069, 0.03773226],
    [0.050079, -0.1971357, 0.2363775, -0.1215569, 0.02853104],
    [0.0720126, -0.225158, 0.2293715, -0.1038359, 0.02174211],
    [0.09585932, -0.2481381, 0.2194469, -0.08862132, 0.01672243],
], dtype=float)
_SF_Q4 = np.asarray([
    [0.4845581, 0.07457381],
    [0.4514531, 0.0645864],
    [0.4228767, 0.05655715],
    [0.3978515, 0.04997164],
    [0.3756999, 0.044488],
], dtype=float)
_SF_P5 = np.asarray([
    [0.003722238, -0.04991683, 0.1658905, -0.238724, 0.1877469, -0.08159855, 0.03051959],
    [0.008182649, -0.07325459, 0.1945697, -0.2396387, 0.1667832, -0.06620786, 0.02224041],
    [0.01466325, -0.09858686, 0.2180684, -0.2347118, 0.1464354, -0.05350728, 0.01624782],
    [0.02314317, -0.1246383, 0.2362036, -0.2257366, 0.1275895, -0.04317874, 0.01193168],
    [0.03346886, -0.1503778, 0.2492826, -0.2142055, 0.1106482, -0.03486024, 0.008821107],
], dtype=float)
_SF_Q5 = np.asarray([0.241882, 0.2291233, 0.2177793, 0.2075784, 0.1983358], dtype=float)
_SF_P6L = np.asarray([
    [0.05613913, -0.3019847, 0.6256387, -0.6324887, 0.3303194],
    [0.06843713, -0.3342119, 0.6302307, -0.5829747, 0.27657],
    [0.08203343, -0.3644705, 0.627866, -0.5335581, 0.2312756],
    [0.09675562, -0.3922489, 0.6197133, -0.485747, 0.1934013],
    [0.1124069, -0.4172349, 0.6069622, -0.4405326, 0.1618978],
], dtype=float)
_SF_Q6L = np.asarray([
    [0.9077644, 0.2535284],
    [0.8626056, 0.22914],
    [0.8212018, 0.2078043],
    [0.7831755, 0.1890848],
    [0.7481828, 0.1726085],
], dtype=float)
_SF_P6U = np.asarray([
    [8.531865e-4, -0.01616105, 0.06888533, -0.1109391, 0.07747182],
    [0.00206076, -0.02558954, 0.08595213, -0.1170228, 0.07094106],
    [0.004028559, -0.03697768, 0.1021332, -0.1201436, 0.06412774],
    [0.006887946, -0.04994202, 0.1168451, -0.1207733, 0.0574421],
    [0.01071895, -0.06404749, 0.1297386, -0.1194208, 0.05112822],
], dtype=float)
_SF_Q6U = np.asarray([
    [1.10127, 0.3858544],
    [1.025431, 0.3337648],
    [0.9599102, 0.2918724],
    [0.9025276, 0.2575336],
    [0.851747, 0.2289667],
], dtype=float)
_SF_P7L = np.asarray([
    [0.02460495, -0.1640964, 0.434011, -0.5705516, 0.4418614],
    [0.03070261, -0.1879546, 0.4565902, -0.5544891, 0.389279],
    [0.03770526, -0.2121608, 0.4746423, -0.5338058, 0.3417026],
    [0.04559398, -0.236267, 0.4881998, -0.5098448, 0.2991635],
    [0.054325, -0.2598752, 0.4974791, -0.4837861, 0.2614838],
], dtype=float)
_SF_Q7L = np.asarray([
    [1.124957, 0.3784976],
    [1.07542, 0.3466086],
    [1.029374, 0.3181219],
    [0.9865496, 0.2926441],
    [0.9466891, 0.2698218],
], dtype=float)
_SF_P7U = np.asarray([
    [1.924318e-4, -0.005044864, 0.02979803, -0.06660688, 0.06792268],
    [5.030909e-4, -0.008639332, 0.04018472, -0.07595456, 0.06696215],
    [0.001059406, -0.01343605, 0.0513536, -0.08386588, 0.06484517],
    [0.001941904, -0.01943727, 0.06288221, -0.09021607, 0.06193],
    [0.003224785, -0.02657664, 0.07438627, -0.09500554, 0.05850884],
], dtype=float)
_SF_Q7U = np.asarray([
    [1.45073, 0.6578685],
    [1.353872, 0.5724332],
    [1.269924, 0.5032139],
    [1.196177, 0.4460948],
    [1.130719, 0.3982785],
], dtype=float)
_SF_P8L = np.asarray([
    [0.0137803, -0.1097846, 0.3625283, -0.6522477, 0.6684458, -0.4703556],
    [0.01721632, -0.1274981, 0.3917226, -0.6562264, 0.6305859, -0.4067119],
    [0.02121871, -0.1461891, 0.4185427, -0.6543539, 0.590466, -0.3507098],
    [0.02580565, -0.1656048, 0.4426283, -0.6473472, 0.5494752, -0.3018936],
    [0.03098251, -0.1854823, 0.4637398, -0.6359482, 0.5086794, -0.2595588],
], dtype=float)
_SF_Q8L = np.asarray([
    [1.076975, 0.3394154],
    [1.036132, 0.3145673],
    [0.9978025, 0.2920529],
    [0.9617584, 0.2715949],
    [0.9278774, 0.2530051],
], dtype=float)
_SF_P8U = np.asarray([
    [4.29046e-5, -0.001508077, 0.01233763, -0.0409127, 0.06547454, -0.05664203],
    [1.201008e-4, -0.002778372, 0.01797999, -0.05055048, 0.07125083, -0.05469912],
    [2.698511e-4, -0.004628815, 0.0247089, -0.06017759, 0.07566434, -0.05202678],
    [5.259595e-4, -0.007144198, 0.03238633, -0.06946769, 0.07873067, -0.0488949],
    [9.255826e-4, -0.01038126, 0.04083176, -0.07815954, 0.08054087, -0.04552077],
], dtype=float)
_SF_Q8U = np.asarray([
    [1.379457, 0.5786953],
    [1.300303, 0.5135748],
    [1.230436, 0.4593779],
    [1.168075, 0.4135871],
    [1.111893, 0.3744076],
], dtype=float)


def _sphconv_family_params(support_radius_pix: float, sphparm: float = 1.0) -> tuple[int, int, int]:
    """Mirror casacore::Sph_Conv family selection."""
    cut = float(support_radius_pix)
    isupp = int(2.0 * cut)
    if isupp < 4:
        isupp = 4
    if isupp > 8:
        isupp = 8

    ialpha = int(2.0 * float(sphparm) + 1.0)
    if ialpha < 1:
        ialpha = 1
    if ialpha > 5:
        ialpha = 5

    jmax = int(cut)
    if jmax > 7:
        jmax = 7
    if jmax < 1:
        raise ValueError('Spheroidal kernel requires int(convsupport) >= 1')
    return ialpha, isupp, jmax


def _eval_schwab_sphfn_gridding(ialpha: int, isupp: int, eta: np.ndarray) -> np.ndarray:
    """Evaluate the Fred Schwab spheroidal approximation used by casacore.

    This matches the gridding branch (IFLAG <= 0) of the original SPHFN routine.
    """
    ax = np.abs(np.asarray(eta, dtype=float))
    out = np.zeros_like(ax, dtype=float)

    valid = np.isfinite(ax) & (ax <= 1.0)
    if not np.any(valid):
        return out

    row = int(ialpha) - 1
    eta2 = ax[valid] * ax[valid]

    if isupp == 4:
        x = eta2 - 1.0
        p = _SF_P4[row]
        q = _SF_Q4[row]
        psi = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * p[4])))) / (x * (q[0] + x * q[1]) + 1.0)
    elif isupp == 5:
        x = eta2 - 1.0
        p = _SF_P5[row]
        q = _SF_Q5[row]
        psi = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * p[6])))))) / (x * q + 1.0)
    elif isupp == 6:
        psi = np.empty_like(eta2)
        m = ax[valid] <= 0.75
        if np.any(m):
            x = eta2[m] - 0.5625
            p = _SF_P6L[row]
            q = _SF_Q6L[row]
            psi[m] = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * p[4])))) / (x * (q[0] + x * q[1]) + 1.0)
        if np.any(~m):
            x = eta2[~m] - 1.0
            p = _SF_P6U[row]
            q = _SF_Q6U[row]
            psi[~m] = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * p[4])))) / (x * (q[0] + x * q[1]) + 1.0)
    elif isupp == 7:
        psi = np.empty_like(eta2)
        m = ax[valid] <= 0.775
        if np.any(m):
            x = eta2[m] - 0.600625
            p = _SF_P7L[row]
            q = _SF_Q7L[row]
            psi[m] = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * p[4])))) / (x * (q[0] + x * q[1]) + 1.0)
        if np.any(~m):
            x = eta2[~m] - 1.0
            p = _SF_P7U[row]
            q = _SF_Q7U[row]
            psi[~m] = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * p[4])))) / (x * (q[0] + x * q[1]) + 1.0)
    elif isupp == 8:
        psi = np.empty_like(eta2)
        m = ax[valid] <= 0.775
        if np.any(m):
            x = eta2[m] - 0.600625
            p = _SF_P8L[row]
            q = _SF_Q8L[row]
            psi[m] = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * p[5]))))) / (x * (q[0] + x * q[1]) + 1.0)
        if np.any(~m):
            x = eta2[~m] - 1.0
            p = _SF_P8U[row]
            q = _SF_Q8U[row]
            psi[~m] = (p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * p[5]))))) / (x * (q[0] + x * q[1]) + 1.0)
    else:
        raise ValueError(f'Unsupported spheroidal family isupp={isupp}')

    # casacore MathFunc2.cc gridding branch: apply (1-eta^2)^alpha except for
    # alpha=0 (ialpha==1) and eta==0 or |eta|==1 special cases.
    if ialpha != 1:
        edge = np.clip(1.0 - eta2, 0.0, None)
        alpha = _SF_ALPHA[row]
        psi = np.where(edge > 0.0, np.power(edge, alpha) * psi, 0.0)

    out[valid] = psi
    out[~np.isfinite(out)] = 0.0
    return out


def _kernel_sf(r_arcsec: np.ndarray, *, cell_arcsec: float, support_radius_pix: float) -> np.ndarray:
    """casacore-like spheroidal kernel."""
    if not np.isfinite(cell_arcsec) or cell_arcsec <= 0:
        raise ValueError("cell_arcsec must be positive for kernel='sf'")
    ialpha, isupp, jmax = _sphconv_family_params(support_radius_pix)
    r_pix = np.asarray(r_arcsec, dtype=float) / float(cell_arcsec)
    eta = r_pix / float(jmax)
    out = _eval_schwab_sphfn_gridding(ialpha, isupp, eta)
    out = np.where(r_pix <= float(support_radius_pix), out, 0.0)
    out[out < 0] = 0.0
    return out


def _evaluate_kernel(
    r_arcsec: np.ndarray,
    *,
    kernel_name: str,
    gwidth_arcsec: float,
    jwidth_arcsec: float,
    support_radius_pix: float,
    cell_arcsec: float,
    beam_fwhm_arcsec: float,
    sf_beam_support_beam: float,
    sf_beam_shape_c: float,
    gjinc_beam_support_beam: float,
    gjinc_beam_gwidth_beam: float,
    gjinc_beam_jwidth_beam: float,
    eps_u0: float,
) -> np.ndarray:
    if kernel_name == 'gjinc':
        return _kernel_gjinc(r_arcsec, gwidth_arcsec, jwidth_arcsec, eps_u0)
    if kernel_name == 'gjinc_beam':
        return _kernel_gjinc_beam(
            r_arcsec,
            beam_fwhm_arcsec=beam_fwhm_arcsec,
            gwidth_beam=gjinc_beam_gwidth_beam,
            jwidth_beam=gjinc_beam_jwidth_beam,
            eps_u0=eps_u0,
        )
    if kernel_name == 'gauss':
        return _kernel_gauss(r_arcsec, gwidth_arcsec)
    if kernel_name == 'sf':
        return _kernel_sf(r_arcsec, cell_arcsec=cell_arcsec, support_radius_pix=support_radius_pix)
    if kernel_name == 'sf_beam':
        return _kernel_sf_beam(
            r_arcsec,
            beam_fwhm_arcsec=beam_fwhm_arcsec,
            support_beam=sf_beam_support_beam,
            shape_c=sf_beam_shape_c,
        )
    raise ValueError(f"Unknown kernel_name={kernel_name!r}")


def _kernel_sf_beam(
    r_arcsec: np.ndarray,
    *,
    beam_fwhm_arcsec: float,
    support_beam: float,
    shape_c: float,
) -> np.ndarray:
    """Continuous beam-aware spheroidal proof-of-concept kernel."""
    if not np.isfinite(beam_fwhm_arcsec) or beam_fwhm_arcsec <= 0:
        raise ValueError("beam_fwhm_arcsec must be positive for kernel='sf_beam'")
    if not np.isfinite(support_beam) or support_beam <= 0:
        raise ValueError("support_beam must be positive for kernel='sf_beam'")
    if not np.isfinite(shape_c) or shape_c <= 0:
        raise ValueError("shape_c must be positive for kernel='sf_beam'")
    support_arcsec = float(support_beam) * float(beam_fwhm_arcsec)
    r = np.abs(np.asarray(r_arcsec, dtype=float))
    out = np.zeros_like(r, dtype=float)
    inside = r < support_arcsec
    if not np.any(inside):
        return out
    x = np.clip(r[inside] / support_arcsec, 0.0, 1.0 - 1.0e-12)
    v0 = float(pro_ang1(0, 0, float(shape_c), 0.0)[0])
    if (not np.isfinite(v0)) or abs(v0) <= 0:
        raise ValueError("pro_ang1 returned an invalid normalization for kernel='sf_beam'")
    vx = np.asarray(pro_ang1(0, 0, float(shape_c), x)[0], dtype=float)
    vals = vx / v0
    vals = np.where(np.isfinite(vals) & (vals > 0.0), vals, 0.0)
    out[inside] = vals
    return out


def _kernel_gjinc_beam(
    r_arcsec: np.ndarray,
    *,
    beam_fwhm_arcsec: float,
    gwidth_beam: float,
    jwidth_beam: float,
    eps_u0: float,
) -> np.ndarray:
    """Beam-aware GJINC kernel with widths expressed in beam-FWHM units."""
    if not np.isfinite(beam_fwhm_arcsec) or beam_fwhm_arcsec <= 0:
        raise ValueError("beam_fwhm_arcsec must be positive for kernel='gjinc_beam'")
    if not np.isfinite(gwidth_beam) or gwidth_beam <= 0:
        raise ValueError("gwidth_beam must be positive for kernel='gjinc_beam'")
    if not np.isfinite(jwidth_beam) or jwidth_beam <= 0:
        raise ValueError("jwidth_beam must be positive for kernel='gjinc_beam'")
    return _kernel_gjinc(
        np.asarray(r_arcsec, dtype=float),
        float(gwidth_beam) * float(beam_fwhm_arcsec),
        float(jwidth_beam) * float(beam_fwhm_arcsec),
        eps_u0=float(eps_u0),
    )


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
    time_all = np.asarray(input_data.time)
    ndump, nchan = spec_all.shape
    flag_all = normalize_row_flag_mask(input_data.flag, ndump=ndump, allow_none=True, name='flag')

    # --- 1. データの前処理・フィルタリング ---
    valid = flag_all.astype(bool, copy=False)
    if config.exclude_turnaround and input_data.is_turnaround is not None:
        valid &= ~_safe_bool_array(input_data.is_turnaround, default=False)

    # 空間・スペクトルの有効性チェック
    valid &= np.isfinite(x_all) & np.isfinite(y_all)
    valid &= np.any(np.isfinite(spec_all), axis=1)

    weight_mode = str(getattr(config, 'weight_mode', 'uniform')).strip().lower()
    rms_all, tint_all, tsys_all = None, None, None
    if input_data.rms is not None:
        rms_all = np.asarray(input_data.rms)
        if weight_mode == 'rms':
            bad_rms = (~np.isfinite(rms_all)) | (rms_all <= 0)
            if np.any(bad_rms):
                n_bad = int(np.count_nonzero(bad_rms))
                raise ValueError(
                    "weight_mode='rms' requires finite positive GridInput.rms for every dump. "
                    f"Found {n_bad}/{len(rms_all)} invalid dumps. "
                    "Provide valid BSL_RMS for all inputs, or set weight_mode='uniform'."
                )
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
    if weight_mode == 'rms':
        if rms is None:
            raise ValueError(
                "weight_mode='rms' requires GridInput.rms for every dump. "
                "Provide valid BSL_RMS through the OTF pipeline, or set weight_mode='uniform'."
            )
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
    sf_beam_support_beam = dtype_np(float(kernel_resolved.get("sf_beam_support_beam", np.nan)))
    sf_beam_shape_c = dtype_np(float(kernel_resolved.get("sf_beam_shape_c", np.nan)))
    gjinc_beam_support_beam = dtype_np(float(kernel_resolved.get("gjinc_beam_support_beam", np.nan)))
    gjinc_beam_gwidth_beam = dtype_np(float(kernel_resolved.get("gjinc_beam_gwidth_beam", np.nan)))
    gjinc_beam_jwidth_beam = dtype_np(float(kernel_resolved.get("gjinc_beam_jwidth_beam", np.nan)))
    gwidth_arcsec = dtype_np(float(kernel_resolved["gwidth_pix"]) * config.cell_arcsec) if np.isfinite(float(kernel_resolved["gwidth_pix"])) else dtype_np(np.nan)
    jwidth_arcsec = dtype_np(float(kernel_resolved["jwidth_pix"]) * config.cell_arcsec) if kernel_name in ("gjinc", "gjinc_beam") else dtype_np(np.nan)
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

        k = _evaluate_kernel(
            r,
            kernel_name=kernel_name,
            gwidth_arcsec=float(gwidth_arcsec),
            jwidth_arcsec=float(jwidth_arcsec),
            support_radius_pix=(float(kernel_resolved['convsupport']) if kernel_name == 'sf' and kernel_resolved.get('convsupport') is not None else float(r_trunc_arcsec / config.cell_arcsec)),
            cell_arcsec=float(config.cell_arcsec),
            beam_fwhm_arcsec=float(config.beam_fwhm_arcsec),
            sf_beam_support_beam=float(sf_beam_support_beam),
            sf_beam_shape_c=float(sf_beam_shape_c),
            gjinc_beam_support_beam=float(gjinc_beam_support_beam),
            gjinc_beam_gwidth_beam=float(gjinc_beam_gwidth_beam),
            gjinc_beam_jwidth_beam=float(gjinc_beam_jwidth_beam),
            eps_u0=float(config.eps_u0),
        ).astype(dtype_np, copy=False)

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
            sf_beam_support_beam=float(sf_beam_support_beam),
            sf_beam_shape_c=float(sf_beam_shape_c),
            sf_convsupport_exact=(float(kernel_resolved['convsupport']) if kernel_resolved.get('convsupport') is not None else None),
            gjinc_beam_support_beam=float(gjinc_beam_support_beam),
            gjinc_beam_gwidth_beam=float(gjinc_beam_gwidth_beam),
            gjinc_beam_jwidth_beam=float(gjinc_beam_jwidth_beam),
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
            sf_beam_support_beam=float(sf_beam_support_beam),
            sf_beam_shape_c=float(sf_beam_shape_c),
            sf_convsupport_exact=(float(kernel_resolved['convsupport']) if kernel_resolved.get('convsupport') is not None else None),
            gjinc_beam_support_beam=float(gjinc_beam_support_beam),
            gjinc_beam_gwidth_beam=float(gjinc_beam_gwidth_beam),
            gjinc_beam_jwidth_beam=float(gjinc_beam_jwidth_beam),
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
            "kernel_preset": kernel_resolved.get("kernel_preset"),
            "kernel_sign": str(kernel_sign),
            "backend": backend_impl,
            "cell_arcsec": float(config.cell_arcsec),
            "beam_fwhm_arcsec": float(config.beam_fwhm_arcsec),
            "beam_pix": float(kernel_resolved.get("beam_pix", config.beam_fwhm_arcsec / config.cell_arcsec)),
            "cell_over_beam": float(kernel_resolved.get("cell_over_beam", config.cell_arcsec / config.beam_fwhm_arcsec)),
            "cell_is_coarse": bool(kernel_resolved.get("cell_is_coarse", False)),
            "bmaj_eff_arcsec": bmaj_eff_arcsec,
            "gwidth_arcsec": float(gwidth_arcsec) if np.isfinite(gwidth_arcsec) else np.nan,
            "jwidth_arcsec": float(jwidth_arcsec) if kernel_name == "gjinc" else np.nan,
            "support_radius_arcsec": float(r_trunc_arcsec),
            "convsupport": kernel_resolved.get("convsupport"),
            "weight_map_semantics": "denom_after_sign",
            "hit_map_semantics": "used_samples_after_sign",
            "nsamp_map_semantics": "support_samples_before_sign",
            "wsum_map_semantics": "signed_wsum_before_sign",
            "wabs_map_semantics": "abs_wsum_before_sign",
            "cancel_map_semantics": "abs_wsum_over_wabs",
            "weight_rel_map_semantics": "abs_weight_over_medabs",
            "gwidth_source": str(kernel_resolved.get("gwidth_source", "")),
            "jwidth_source": str(kernel_resolved.get("jwidth_source", "")),
            "support_source": str(kernel_resolved.get("support_source", "")),
            "mask_map_semantics": "spatial_valid_pre_chanmask",
            "weight_mode": str(getattr(config, 'weight_mode', 'uniform')).lower(),
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
