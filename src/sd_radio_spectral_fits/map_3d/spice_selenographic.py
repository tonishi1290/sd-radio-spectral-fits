"""SPICE-backed Moon selenographic coordinate helpers.

This module is intentionally imported only from the ``moon_selenographic``
coordinate path.  The heavy optional dependency ``spiceypy`` is imported lazily
inside helper functions, so normal RA/Dec, Galactic, and dAz/dEl mapping does
not require SPICE.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from ..scantable_utils import _resolve_table_timestamps
from ..doppler import earth_location_from_meta
from .pointing_correction import resolve_pointing_offsets_arcsec


_DEFAULT_MOON_FRAME = "IAU_MOON"
_DEFAULT_ABCORR = "NONE"


def _require_spiceypy():
    try:
        import spiceypy as spice  # type: ignore
        return spice
    except Exception as exc:
        raise ImportError(
            "coord_sys='moon_selenographic' requires SpiceyPy. "
            "Install it with: pip install spiceypy. "
            "This dependency is imported only for SPICE-backed Moon surface mapping."
        ) from exc


def _require_astropy():
    try:
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astropy.time import Time
        return u, SkyCoord, Time
    except Exception as exc:
        raise ImportError(
            "coord_sys='moon_selenographic' requires astropy as well as SpiceyPy, "
            "because the input RA/DEC or GLON/GLAT beam coordinates and observatory "
            "position are converted to inertial vectors before SPICE surface intersection."
        ) from exc


def _meta_get_any(meta: dict, *keys):
    for key in keys:
        if key in meta and meta[key] not in (None, ""):
            return meta[key]
    return None


def _combined_meta_for_location(scantable, table: pd.DataFrame) -> dict:
    meta = dict(getattr(scantable, "meta", {}) or {})
    for key in (
        "OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z", "OBSGEO_X", "OBSGEO_Y", "OBSGEO_Z",
        "SITELAT", "SITELONG", "SITELON", "SITEELEV", "SITE_LAT", "SITE_LON", "SITE_ELEV",
        "site_lat_deg", "site_lon_deg", "site_height_m", "site_name", "SITENAME", "OBS_SITE", "SITE",
    ):
        if key in meta and meta[key] not in (None, ""):
            continue
        if key in table.columns:
            arr = pd.Series(table[key]).dropna()
            if len(arr):
                meta[key] = arr.iloc[0]
    return meta


def _resolve_time_for_astropy(table: pd.DataFrame):
    ts = _resolve_table_timestamps(table)
    if ts is None or len(ts) == 0 or ts.isna().all():
        raise ValueError(
            "SPICE Moon selenographic mapping requires real UTC timestamps. "
            "Provide TIMESTAMP, MJD, or DATE-OBS/DATEOBS(+TIME) columns."
        )
    return ts


def _as_kernel_list(value) -> list[str]:
    if value in (None, ""):
        return []
    if isinstance(value, (str, os.PathLike)):
        text = os.fspath(value)
        # Accept either os.pathsep-separated strings or newline/comma separated metadata.
        parts: list[str] = []
        for chunk in text.replace("\n", os.pathsep).replace(",", os.pathsep).split(os.pathsep):
            chunk = chunk.strip()
            if chunk:
                parts.append(chunk)
        return parts
    if isinstance(value, Iterable):
        out = []
        for item in value:
            out.extend(_as_kernel_list(item))
        return out
    return [str(value)]


def _load_spice_kernels_if_requested(spice, scantable, table: pd.DataFrame) -> list[str]:
    """Load kernels listed in metadata or environment, if supplied.

    The normal SPICE workflow, where the caller has already run ``spice.furnsh``
    before the mapping pipeline, is also supported.  In that case this function
    simply loads nothing.
    """
    meta = dict(getattr(scantable, "meta", {}) or {})
    env_value = os.environ.get("SDRSF_SPICE_KERNELS")
    candidates = []
    candidates.extend(_as_kernel_list(env_value))
    for key in (
        "SPICE_KERNELS", "SPICEKERNELS", "SPICE_KERNEL", "SPICEKRN", "SPICE_KRN",
        "KERNELS", "KERNEL_LIST",
    ):
        candidates.extend(_as_kernel_list(_meta_get_any(meta, key)))
    if "SPICE_KERNELS" in table.columns:
        vals = pd.Series(table["SPICE_KERNELS"]).dropna().astype(str)
        if len(vals):
            candidates.extend(_as_kernel_list(vals.iloc[0]))

    loaded = []
    seen = set()
    for kernel in candidates:
        path = str(Path(kernel).expanduser())
        if path in seen:
            continue
        seen.add(path)
        spice.furnsh(path)
        loaded.append(path)
    return loaded


def _moon_frame_from_meta(scantable) -> str:
    meta = getattr(scantable, "meta", {}) or {}
    value = _meta_get_any(meta, "SPICE_MOON_FRAME", "MOON_FRAME", "SELFRAME")
    if value in (None, ""):
        return _DEFAULT_MOON_FRAME
    return str(value).strip()


def _abcorr_from_meta(scantable) -> str:
    meta = getattr(scantable, "meta", {}) or {}
    value = _meta_get_any(meta, "SPICE_ABCORR", "ABCORR")
    if value in (None, ""):
        return _DEFAULT_ABCORR
    text = str(value).strip().upper()
    # The current implementation is a geometric ray/ellipsoid intersection.  It
    # accepts the keyword for provenance, but the actual manual intersection is
    # defined at the observation epoch and uses geometric Earth/Moon states.
    if text not in {"NONE", "LT", "LT+S", "CN", "CN+S"}:
        raise ValueError(f"Unsupported SPICE_ABCORR={text!r}.")
    return text


def _table_or_meta_get_any(scantable, table: pd.DataFrame, *keys, default=None):
    meta = getattr(scantable, "meta", {}) or {}
    for key in keys:
        if key in meta and meta[key] not in (None, ""):
            return meta[key]
        for cand in (key, key.upper(), key.lower()):
            if cand in table.columns:
                vals = pd.Series(table[cand]).dropna()
                if len(vals):
                    return vals.iloc[0]
    return default


def _moving_body_radec_input_mode(scantable, table: pd.DataFrame) -> str:
    """Return how RA/DEC columns should be interpreted for Moon SPICE rays."""
    value = _table_or_meta_get_any(
        scantable,
        table,
        "MOVING_BODY_RADEC_FRAME",
        "BODY_RADEC_FRAME",
        "BODY_OFFSET_RADEC_FRAME",
        "MOVBODY_RADEC_FRAME",
        default="icrs",
    )
    text = str(value).strip().lower().replace("-", "_")
    aliases = {
        "": "icrs",
        "icrs": "icrs",
        "j2000": "icrs",
        "fk5": "fk5",
        "fk4": "fk4",
        "native": "apparent",
        "apparent": "apparent",
        "topocentric": "apparent",
        "topocentric_apparent": "apparent",
        "gcrs": "apparent",
    }
    if text not in aliases:
        raise ValueError(
            "Unsupported MOVING_BODY_RADEC_FRAME for moon_selenographic "
            f"{value!r}. Use 'ICRS' (default), 'apparent', 'FK5', or 'FK4'."
        )
    return aliases[text]


def _beam_coord_in_moon_apparent_frame(scantable, table: pd.DataFrame, *, obstime, location):
    u, SkyCoord, _Time = _require_astropy()
    from astropy.coordinates import FK4, FK5, get_body

    body_coord = get_body("moon", obstime, location=location)
    body_frame = body_coord.frame.replicate_without_data()
    if {"RA", "DEC"}.issubset(table.columns):
        ra = pd.to_numeric(table["RA"], errors="coerce").to_numpy(dtype=float)
        dec = pd.to_numeric(table["DEC"], errors="coerce").to_numpy(dtype=float)
        mode = _moving_body_radec_input_mode(scantable, table)
        if mode == "apparent":
            return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=body_frame), body_coord, mode
        if mode == "icrs":
            return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs").transform_to(body_frame), body_coord, mode
        equinox = _table_or_meta_get_any(scantable, table, "EQUINOX", "RADECEQNX", default=2000.0)
        if mode == "fk5":
            return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=FK5(equinox=f"J{float(equinox):.6f}")).transform_to(body_frame), body_coord, mode
        if mode == "fk4":
            return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=FK4(equinox=f"B{float(equinox):.6f}")).transform_to(body_frame), body_coord, mode
    elif {"GLON", "GLAT"}.issubset(table.columns):
        glon = pd.to_numeric(table["GLON"], errors="coerce").to_numpy(dtype=float)
        glat = pd.to_numeric(table["GLAT"], errors="coerce").to_numpy(dtype=float)
        return SkyCoord(l=glon * u.deg, b=glat * u.deg, frame="galactic").transform_to(body_frame), body_coord, "galactic"
    else:
        raise ValueError(
            "coord_sys='moon_selenographic' requires per-beam RA/DEC or GLON/GLAT columns. "
            "BORE_AZ/BORE_EL are deliberately not used."
        )
    raise RuntimeError("unreachable Moon SPICE beam-coordinate branch")


def _unit_vectors_from_table(scantable, table: pd.DataFrame, *, obstime, location):
    beam_apparent, _body_coord, input_mode = _beam_coord_in_moon_apparent_frame(
        scantable, table, obstime=obstime, location=location
    )
    coord = beam_apparent.icrs
    table["BODY_RADEC_INPUT"] = input_mode
    cart = coord.cartesian.xyz.to_value()
    return np.asarray(cart, dtype=float).T


def _unit_vectors_from_table_with_optional_pointing_correction(
    scantable,
    table: pd.DataFrame,
    *,
    obstime,
    location,
    pointing_correction=None,
    pointing_dx_arcsec=None,
    pointing_dy_arcsec=None,
    pointing_table_index: int | None = None,
):
    """Return J2000/ICRS unit vectors after optional Moon-sky pointing correction."""
    corr_dx, corr_dy = resolve_pointing_offsets_arcsec(
        table,
        pointing_correction=pointing_correction,
        pointing_dx_arcsec=pointing_dx_arcsec,
        pointing_dy_arcsec=pointing_dy_arcsec,
        table_index=pointing_table_index,
    )
    if not (np.any(corr_dx != 0.0) or np.any(corr_dy != 0.0)):
        unit = _unit_vectors_from_table(scantable, table, obstime=obstime, location=location)
        table["POINT_DX_ARCSEC"] = corr_dx
        table["POINT_DY_ARCSEC"] = corr_dy
        return unit

    u, _SkyCoord, _Time = _require_astropy()

    # The correction is defined in the apparent Moon-centered sky tangent plane.
    # Build beam directions in the same apparent frame as astropy.get_body('moon'),
    # add (dx, dy), then convert the corrected direction back to ICRS/J2000 unit vectors.
    beam_apparent, body_coord, input_mode = _beam_coord_in_moon_apparent_frame(
        scantable, table, obstime=obstime, location=location
    )
    table["BODY_RADEC_INPUT"] = input_mode

    dx0, dy0 = body_coord.spherical_offsets_to(beam_apparent)
    corrected = body_coord.spherical_offsets_by(
        dx0 + corr_dx * u.arcsec,
        dy0 + corr_dy * u.arcsec,
    )
    corrected_icrs = corrected.icrs
    table["POINT_DX_ARCSEC"] = corr_dx
    table["POINT_DY_ARCSEC"] = corr_dy
    table["RA_POINT_CORR"] = corrected_icrs.ra.to_value(u.deg)
    table["DEC_POINT_CORR"] = corrected_icrs.dec.to_value(u.deg)
    cart = corrected_icrs.cartesian.xyz.to_value()
    return np.asarray(cart, dtype=float).T

def _observer_geocentric_vectors_km(location, obstime) -> np.ndarray:
    u, _SkyCoord, _Time = _require_astropy()
    pos, _vel = location.get_gcrs_posvel(obstime)
    return np.asarray(pos.xyz.to_value(u.km), dtype=float).T


def ray_ellipsoid_intersections(
    observer_body_fixed_km: np.ndarray,
    ray_body_fixed_unit: np.ndarray,
    radii_km: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Intersect observer rays with a tri-axial ellipsoid.

    Parameters
    ----------
    observer_body_fixed_km
        Observer position relative to the body center in body-fixed coordinates,
        shape ``(N, 3)``.
    ray_body_fixed_unit
        Unit line-of-sight vectors in the same body-fixed frame, shape ``(N, 3)``.
    radii_km
        Ellipsoid radii ``(a, b, c)`` in km.

    Returns
    -------
    points, valid
        ``points`` has shape ``(N, 3)`` and contains NaN rows for rays that do
        not hit the ellipsoid in the forward direction.  ``valid`` is a boolean
        mask.
    """
    p = np.asarray(observer_body_fixed_km, dtype=float)
    d = np.asarray(ray_body_fixed_unit, dtype=float)
    r = np.asarray(radii_km, dtype=float)
    if p.shape != d.shape or p.ndim != 2 or p.shape[1] != 3:
        raise ValueError("observer_body_fixed_km and ray_body_fixed_unit must both have shape (N, 3).")
    if r.shape != (3,) or np.any(~np.isfinite(r)) or np.any(r <= 0):
        raise ValueError("radii_km must contain three positive finite radii.")

    inv_r2 = 1.0 / (r * r)
    a = np.sum(d * d * inv_r2, axis=1)
    b = 2.0 * np.sum(p * d * inv_r2, axis=1)
    c = np.sum(p * p * inv_r2, axis=1) - 1.0
    disc = b * b - 4.0 * a * c
    valid = np.isfinite(disc) & (disc >= 0.0) & np.isfinite(a) & (a > 0.0)

    t = np.full(p.shape[0], np.nan, dtype=float)
    if np.any(valid):
        sqrt_disc = np.sqrt(np.maximum(disc[valid], 0.0))
        a_v = a[valid]
        b_v = b[valid]
        t1 = (-b_v - sqrt_disc) / (2.0 * a_v)
        t2 = (-b_v + sqrt_disc) / (2.0 * a_v)
        t_near = np.where((t1 > 0.0) & np.isfinite(t1), t1, t2)
        idx = np.where(valid)[0]
        ok = (t_near > 0.0) & np.isfinite(t_near)
        t[idx[ok]] = t_near[ok]

    hit = np.isfinite(t)
    points = np.full_like(p, np.nan, dtype=float)
    points[hit] = p[hit] + t[hit, None] * d[hit]
    return points, hit


def lonlat_from_body_fixed_points_deg(points_body_fixed_km: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return east-positive longitude and planetocentric latitude in degrees."""
    p = np.asarray(points_body_fixed_km, dtype=float)
    lon = np.full(p.shape[0], np.nan, dtype=float)
    lat = np.full(p.shape[0], np.nan, dtype=float)
    finite = np.all(np.isfinite(p), axis=1)
    if np.any(finite):
        x = p[finite, 0]
        y = p[finite, 1]
        z = p[finite, 2]
        rho = np.sqrt(x * x + y * y + z * z)
        lon[finite] = np.degrees(np.arctan2(y, x)) % 360.0
        lat[finite] = np.degrees(np.arcsin(np.clip(z / rho, -1.0, 1.0)))
    return lon, lat


def compute_moon_selenographic_lonlat(
    scantable,
    table: pd.DataFrame,
    *,
    pointing_correction=None,
    pointing_dx_arcsec=None,
    pointing_dy_arcsec=None,
    pointing_table_index: int | None = None,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, dict]:
    """Compute Moon surface lon/lat for each beam row using SPICE.

    The caller must provide either already loaded SPICE kernels or kernel paths
    via ``scantable.meta['SPICE_KERNELS']`` or the ``SDRSF_SPICE_KERNELS``
    environment variable.  Required kernels are at least a leapseconds kernel,
    an ephemeris SPK, and a PCK or binary PCK providing Moon radii/frame data.
    """
    spice = _require_spiceypy()
    u, _SkyCoord, Time = _require_astropy()

    table = pd.DataFrame(table).copy()
    loaded_kernels = _load_spice_kernels_if_requested(spice, scantable, table)
    moon_frame = _moon_frame_from_meta(scantable)
    abcorr = _abcorr_from_meta(scantable)

    timestamps = _resolve_time_for_astropy(table)
    obstime = Time(timestamps.to_pydatetime(), scale="utc")
    et = np.asarray([spice.utc2et(t) for t in obstime.utc.isot], dtype=float)

    try:
        _dim, radii = spice.bodvrd("MOON", "RADII", 3)
    except Exception as exc:
        raise RuntimeError(
            "SPICE Moon radii are unavailable. Load a planetary constants kernel, "
            "for example pck00011.tpc, before using coord_sys='moon_selenographic'."
        ) from exc
    radii = np.asarray(radii, dtype=float)

    location = earth_location_from_meta(_combined_meta_for_location(scantable, table))
    beam_unit_j2000 = _unit_vectors_from_table_with_optional_pointing_correction(
        scantable,
        table,
        obstime=obstime,
        location=location,
        pointing_correction=pointing_correction,
        pointing_dx_arcsec=pointing_dx_arcsec,
        pointing_dy_arcsec=pointing_dy_arcsec,
        pointing_table_index=pointing_table_index,
    )
    obs_geo_j2000 = _observer_geocentric_vectors_km(location, obstime)

    points = np.full((len(table), 3), np.nan, dtype=float)
    valid = np.zeros(len(table), dtype=bool)
    for i, eti in enumerate(et):
        if not np.all(np.isfinite(beam_unit_j2000[i])) or not np.all(np.isfinite(obs_geo_j2000[i])):
            continue
        try:
            # Geometric Earth center relative to Moon center.  The current manual
            # ray/ellipsoid intersection is explicitly geometric at the observation
            # epoch.  ``abcorr`` is still recorded as provenance so users know the
            # chosen convention, and future sincpt/custom-observer support can use it.
            earth_from_moon, _lt = spice.spkpos("EARTH", float(eti), "J2000", "NONE", "MOON")
            j2moon = spice.pxform("J2000", moon_frame, float(eti))
        except Exception as exc:
            raise RuntimeError(
                "SPICE state/frame transform failed. Check that the loaded kernels cover "
                f"the observation time and define the Moon frame {moon_frame!r}."
            ) from exc
        obs_from_moon_j2000 = np.asarray(earth_from_moon, dtype=float) + obs_geo_j2000[i]
        obs_body = np.asarray(j2moon, dtype=float) @ obs_from_moon_j2000
        ray_body = np.asarray(j2moon, dtype=float) @ beam_unit_j2000[i]
        norm = np.linalg.norm(ray_body)
        if not np.isfinite(norm) or norm <= 0:
            continue
        ray_body = ray_body / norm
        hit_point, hit_mask = ray_ellipsoid_intersections(obs_body[None, :], ray_body[None, :], radii)
        if bool(hit_mask[0]):
            points[i] = hit_point[0]
            valid[i] = True

    lon_deg, lat_deg = lonlat_from_body_fixed_points_deg(points)
    table["SELON"] = lon_deg
    table["SELAT"] = lat_deg
    table["SPICE_HIT"] = valid
    table["SPICE_X_KM"] = points[:, 0]
    table["SPICE_Y_KM"] = points[:, 1]
    table["SPICE_Z_KM"] = points[:, 2]

    meta = {
        "spice_loaded_kernels": loaded_kernels,
        "spice_moon_frame": moon_frame,
        "spice_abcorr": abcorr,
        "spice_body": "MOON",
        "spice_surface": "ELLIPSOID",
        "spice_lon_positive": "EAST",
        "spice_latitude": "PLANETOCENTRIC",
        "spice_valid_count": int(np.count_nonzero(valid)),
        "spice_pointing_correction_applied": bool(np.any(table.get("POINT_DX_ARCSEC", 0) != 0) or np.any(table.get("POINT_DY_ARCSEC", 0) != 0)),
        "spice_row_count": int(len(valid)),
    }
    return table, lon_deg, lat_deg, meta
