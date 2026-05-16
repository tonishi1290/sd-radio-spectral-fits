import re
import numpy as np


_AZEL_OFFSET_GENERIC_ALIASES = {
    "azel_offset",
    "altaz_offset",
    "dazel",
    "dazdel",
    "az_el_offset",
    "az-el-offset",
    "az_el",
}

_AZEL_OFFSET_BODY_ALIASES = {
    "moon": "moon",
    "luna": "moon",
    "tsuki": "moon",
    "sun": "sun",
    "solar": "sun",
}


def _norm_token(value) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def normalize_moving_body_name(value):
    """Return 'moon' or 'sun' from a loose body name, or None."""
    if value is None:
        return None
    token = _norm_token(value)
    if token in _AZEL_OFFSET_BODY_ALIASES:
        return _AZEL_OFFSET_BODY_ALIASES[token]
    # SDFITS OBJECT names are often free-form, so allow containment checks.
    if "moon" in token or "luna" in token:
        return "moon"
    if token == "sun" or "solar" in token:
        return "sun"
    return None


def normalize_azel_offset_coord_sys(coord_sys, *, ref_coord=None, object_name=None):
    """Normalize moving-object Az/El-offset coord_sys names.

    Returns
    -------
    (is_azel_offset, body, canonical_coord_sys)
        body is 'moon', 'sun', or None for generic offset frames.

    Accepted canonical names are 'moon_azel_offset' and 'sun_azel_offset'.
    The generic aliases 'azel_offset'/'dazel' are accepted when the body can be
    inferred from ref_coord or OBJECT metadata.
    """
    raw = "" if coord_sys is None else str(coord_sys).strip()
    token = _norm_token(raw)
    token_compact = token.replace("_", "")

    for body in ("moon", "sun"):
        patterns = {
            f"{body}_azel_offset",
            f"{body}_altaz_offset",
            f"{body}_dazel",
            f"{body}_dazdel",
            f"{body}_az_el_offset",
            f"{body}_az_el",
            f"{body}_offset",
        }
        compact_patterns = {p.replace("_", "") for p in patterns}
        if token in patterns or token_compact in compact_patterns:
            return True, body, f"{body}_azel_offset"

    if token in _AZEL_OFFSET_GENERIC_ALIASES or token_compact in {a.replace("_", "") for a in _AZEL_OFFSET_GENERIC_ALIASES}:
        body = normalize_moving_body_name(ref_coord) or normalize_moving_body_name(object_name)
        canonical = f"{body}_azel_offset" if body is not None else "azel_offset"
        return True, body, canonical

    return False, None, raw


_SKY_OFFSET_GENERIC_ALIASES = {
    "sky_offset",
    "radec_offset",
    "ra_dec_offset",
    "equatorial_offset",
    "apparent_sky_offset",
}


def normalize_sky_offset_coord_sys(coord_sys, *, ref_coord=None, object_name=None):
    """Normalize moving-object sky-plane offset coord_sys names.

    Returns
    -------
    (is_sky_offset, body, canonical_coord_sys)
        body is 'moon', 'sun', or None for generic sky-offset aliases.

    Accepted canonical names are 'moon_sky_offset' and 'sun_sky_offset'.
    These are Moon/Sun-centered apparent sky-plane offsets, not AltAz offsets.
    For the Moon this is the practical frame that removes the parallactic/AltAz
    field rotation from dAz/dEl scans while keeping the apparent disk about
    30 arcmin across.
    """
    raw = "" if coord_sys is None else str(coord_sys).strip()
    token = _norm_token(raw)
    token_compact = token.replace("_", "")

    for body in ("moon", "sun"):
        patterns = {
            f"{body}_sky_offset",
            f"{body}_radec_offset",
            f"{body}_ra_dec_offset",
            f"{body}_equatorial_offset",
            f"{body}_apparent_offset",
            f"{body}_apparent_sky_offset",
        }
        compact_patterns = {p.replace("_", "") for p in patterns}
        if token in patterns or token_compact in compact_patterns:
            return True, body, f"{body}_sky_offset"

    aliases_compact = {a.replace("_", "") for a in _SKY_OFFSET_GENERIC_ALIASES}
    if token in _SKY_OFFSET_GENERIC_ALIASES or token_compact in aliases_compact:
        body = normalize_moving_body_name(ref_coord) or normalize_moving_body_name(object_name)
        canonical = f"{body}_sky_offset" if body is not None else "sky_offset"
        return True, body, canonical

    return False, None, raw



_MOON_DISK_OFFSET_ALIASES = {
    "moon_disk_offset",
    "lunar_disk_offset",
    "moon_body_offset",
    "lunar_body_offset",
    "moon_oriented_offset",
    "lunar_oriented_offset",
    "moon_apparent_disk_offset",
    "lunar_apparent_disk_offset",
    "moon_disk",
    "lunar_disk",
}

_MOON_DISK_GENERIC_ALIASES = {
    "disk_offset",
    "body_offset",
    "oriented_offset",
    "apparent_disk_offset",
    "disk",
}


def normalize_moon_disk_offset_coord_sys(coord_sys, *, ref_coord=None, object_name=None):
    """Normalize apparent Moon-disk offset coordinate-system names.

    Returns
    -------
    (is_moon_disk_offset, body, canonical_coord_sys)
        body is currently always ``"moon"`` when true.

    Notes
    -----
    ``moon_disk_offset`` is not a selenographic longitude/latitude grid.  It is
    an apparent 30-arcmin lunar disk map: first compute the Moon-centered
    apparent sky-plane offset, then rotate the local axes so +Y is along the
    apparent lunar north pole and +X is 90 deg eastward from that direction on
    the sky.  This removes AltAz field rotation while keeping the Moon as a
    small apparent disk.
    """
    raw = "" if coord_sys is None else str(coord_sys).strip()
    token = _norm_token(raw)
    token_compact = token.replace("_", "")
    aliases_compact = {a.replace("_", "") for a in _MOON_DISK_OFFSET_ALIASES}
    if token in _MOON_DISK_OFFSET_ALIASES or token_compact in aliases_compact:
        return True, "moon", "moon_disk_offset"
    generic_compact = {a.replace("_", "") for a in _MOON_DISK_GENERIC_ALIASES}
    if token in _MOON_DISK_GENERIC_ALIASES or token_compact in generic_compact:
        body = normalize_moving_body_name(ref_coord) or normalize_moving_body_name(object_name)
        if body == "moon":
            return True, "moon", "moon_disk_offset"
        canonical = f"{body}_disk_offset" if body is not None else "disk_offset"
        return True, body, canonical
    return False, None, raw


_SELENOGRAPHIC_ALIASES = {
    "selenographic",
    "moon_selenographic",
    "lunar_selenographic",
    "moon_spice_selenographic",
    "spice_selenographic",
    "moon_fixed",
    "lunar_fixed",
    "moon_body_fixed",
    "lunar_body_fixed",
}


def normalize_selenographic_coord_sys(coord_sys, *, ref_coord=None, object_name=None):
    """Normalize Moon body-fixed / selenographic coordinate-system names.

    Returns
    -------
    (is_selenographic, body, canonical_coord_sys)
        body is currently always ``"moon"`` when true.

    Notes
    -----
    This is deliberately separate from ``normalize_azel_offset_coord_sys``.
    ``moon_azel_offset`` is an observer-topocentric offset plane, while
    ``moon_selenographic`` is a Moon body-fixed surface coordinate grid.
    """
    raw = "" if coord_sys is None else str(coord_sys).strip()
    token = _norm_token(raw)
    token_compact = token.replace("_", "")
    aliases_compact = {a.replace("_", "") for a in _SELENOGRAPHIC_ALIASES}
    if token in _SELENOGRAPHIC_ALIASES or token_compact in aliases_compact:
        return True, "moon", "moon_selenographic"
    if token in {"body_fixed", "fixed", "surface"}:
        body = normalize_moving_body_name(ref_coord) or normalize_moving_body_name(object_name)
        if body == "moon":
            return True, "moon", "moon_selenographic"
    return False, None, raw


def is_selenographic_coord_sys(coord_sys) -> bool:
    return normalize_selenographic_coord_sys(coord_sys)[0]

def is_azel_offset_coord_sys(coord_sys) -> bool:
    return normalize_azel_offset_coord_sys(coord_sys)[0]


def is_sky_offset_coord_sys(coord_sys) -> bool:
    return normalize_sky_offset_coord_sys(coord_sys)[0]


def is_moon_disk_offset_coord_sys(coord_sys) -> bool:
    return normalize_moon_disk_offset_coord_sys(coord_sys)[0]


def project_to_plane(
    lon_deg: np.ndarray, 
    lat_deg: np.ndarray, 
    lon0: float, 
    lat0: float, 
    projection: str = "SFL", 
    invert_x: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    """球面座標を局所平面 (x, y) [arcsec] に投影する。サポートは CAR と SFL のみ。"""
    lon_mod = np.mod(lon_deg, 360.0)
    lon0_mod = np.mod(lon0, 360.0)
    
    dlon = ((lon_mod - lon0_mod + 180.0) % 360.0) - 180.0
    dlat = lat_deg - lat0
    sgn = -1.0 if invert_x else 1.0
    proj = str(projection).upper()
    
    if proj in ("CAR", "PLATE_CARREE", "NONE"):
        cosf = 1.0
    elif proj in ("GLS", "SFL", "SINE"):
        cosf = np.cos(np.deg2rad(lat_deg))
    else:
        # TANの記述を削除し、エラーメッセージから除外
        raise ValueError(f"Unknown projection: {proj}. Supported: CAR, SFL.")
        
    x_arcsec = sgn * dlon * cosf * 3600.0
    y_arcsec = dlat * 3600.0
    
    return x_arcsec, y_arcsec


def project_azel_offset_to_plane(
    az_deg: np.ndarray,
    el_deg: np.ndarray,
    body_az_deg: np.ndarray,
    body_el_deg: np.ndarray,
    projection: str = "SFL",
) -> tuple[np.ndarray, np.ndarray]:
    """Project per-row telescope AltAz minus moving-body AltAz to arcsec.

    x is positive toward increasing azimuth.  For SFL/GLS/SINE the azimuth
    difference is multiplied by cos(El_body), giving the usual cross-elevation
    angular offset.  For CAR/NONE the raw dAz is used.
    """
    az = np.asarray(az_deg, dtype=float)
    el = np.asarray(el_deg, dtype=float)
    baz = np.asarray(body_az_deg, dtype=float)
    bel = np.asarray(body_el_deg, dtype=float)
    daz = ((az - baz + 180.0) % 360.0) - 180.0
    del_ = el - bel
    proj = str(projection).upper()
    if proj in ("CAR", "PLATE_CARREE", "NONE"):
        cosf = 1.0
    elif proj in ("GLS", "SFL", "SINE"):
        cosf = np.cos(np.deg2rad(bel))
    else:
        raise ValueError(f"Unknown projection: {proj}. Supported: CAR, SFL.")
    return daz * cosf * 3600.0, del_ * 3600.0




def project_sky_offset_to_plane(
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    body_lon_deg: np.ndarray,
    body_lat_deg: np.ndarray,
    projection: str = "SFL",
) -> tuple[np.ndarray, np.ndarray]:
    """Project per-row beam RA/Dec minus moving-body RA/Dec to arcsec.

    x is positive toward increasing right ascension / sky longitude.  For
    SFL/GLS/SINE, the longitude difference is multiplied by cos(Dec_body),
    yielding the small-angle sky-plane offset that removes AltAz field rotation.
    For CAR/NONE, the raw longitude difference is used.
    """
    lon = np.asarray(lon_deg, dtype=float)
    lat = np.asarray(lat_deg, dtype=float)
    blon = np.asarray(body_lon_deg, dtype=float)
    blat = np.asarray(body_lat_deg, dtype=float)
    dlon = ((lon - blon + 180.0) % 360.0) - 180.0
    dlat = lat - blat
    proj = str(projection).upper()
    if proj in ("CAR", "PLATE_CARREE", "NONE"):
        cosf = 1.0
    elif proj in ("GLS", "SFL", "SINE"):
        cosf = np.cos(np.deg2rad(blat))
    else:
        raise ValueError(f"Unknown projection: {proj}. Supported: CAR, SFL.")
    return dlon * cosf * 3600.0, dlat * 3600.0


def build_spatial_wcs_dict(
    coord_sys: str, projection: str, 
    lon0: float, lat0: float, 
    config_x0: float, config_y0: float, 
    cell_arcsec: float, 
    nx: int, ny: int
) -> dict:
    """FITS出力用の空間WCSヘッダ辞書を構築する (CAR/SFL 厳密対応版)"""
    is_offset, body, canonical = normalize_azel_offset_coord_sys(coord_sys)
    is_sky_offset, sky_body, sky_canonical = normalize_sky_offset_coord_sys(coord_sys)
    is_disk_offset, disk_body, disk_canonical = normalize_moon_disk_offset_coord_sys(coord_sys)
    sys_upper = str(disk_canonical if is_disk_offset else (sky_canonical if is_sky_offset else canonical)).upper()
    proj_upper = projection.upper()
    
    # 【修正】未知の投影法は弾く（CARへの無言フォールバックをやめる）
    if proj_upper in ("GLS", "SFL", "SINE"):
        proj_fits = "SFL"
    elif proj_upper in ("CAR", "PLATE_CARREE", "NONE"):
        proj_fits = "CAR"
    else:
        raise ValueError(f"Unknown projection for WCS header: {proj_upper}. Supported: CAR, SFL.")

    is_seleno, seleno_body, seleno_canonical = normalize_selenographic_coord_sys(coord_sys)

    # Moon selenographic maps are fixed body-surface longitude/latitude grids.
    # Use generic LON/LAT projected WCS rather than faking RA/Dec; semantic
    # keywords document that the axes are Moon-fixed selenographic coordinates.
    if is_seleno:
        cdelt = float(cell_arcsec) / 3600.0
        pix_target_x = 1.0 - (float(config_x0) / float(cell_arcsec))
        pix_target_y = 1.0 - (float(config_y0) / float(cell_arcsec))
        common = {
            "WCSNAME": "MOON-SELENO",
            "COORDSYS": "SELENOGRAPHIC",
            "BODY": "MOON",
            "SELONDIR": "EAST",
            "SELLAT": "PLANETOCENTRIC",
            "PROJTYPE": proj_fits,
            "SELON0": float(lon0),
            "SELAT0": float(lat0),
        }
        if proj_fits == "CAR":
            # Use ordinary linear non-celestial WCS axes.  Astropy/WCSLIB and
            # generic FITS viewers accept arbitrary linear CTYPE strings such as
            # SELON/SELAT, whereas pseudo-celestial LON--CAR/LAT--CAR is not a
            # recognized FITS-WCS celestial pair.
            out = {
                "CTYPE1": "SELON", "CRVAL1": float(lon0), "CDELT1": cdelt, "CRPIX1": float(pix_target_x), "CUNIT1": "deg",
                "CTYPE2": "SELAT", "CRVAL2": float(lat0), "CDELT2": cdelt, "CRPIX2": float(pix_target_y), "CUNIT2": "deg",
            }
        else:
            # SFL gridding coordinates are projected surface offsets, not a
            # simple linear longitude/latitude WCS.  Keep the FITS interoperable
            # by labeling them as linear offsets and document the definition.
            out = {
                "CTYPE1": "SELOFFX", "CRVAL1": 0.0, "CDELT1": cdelt, "CRPIX1": float(pix_target_x), "CUNIT1": "deg",
                "CTYPE2": "SELOFFY", "CRVAL2": 0.0, "CDELT2": cdelt, "CRPIX2": float(pix_target_y), "CUNIT2": "deg",
                "OFFXDEF": "(SELON-SELON0)*cos(SELAT)",
                "OFFYDEF": "SELAT-SELAT0",
            }
        out.update(common)
        return out

    # Apparent Moon-disk offsets keep the Moon as a small apparent disk but rotate
    # the axes so +Y is along the apparent lunar north pole.  This is a linear
    # non-celestial WCS, not RA/Dec and not selenographic longitude/latitude.
    if is_disk_offset:
        cdelt = float(cell_arcsec) / 3600.0
        pix_target_x = 1.0 - (float(config_x0) / float(cell_arcsec))
        pix_target_y = 1.0 - (float(config_y0) / float(cell_arcsec))
        return {
            "CTYPE1": "OFFSETX", "CRVAL1": 0.0, "CDELT1": cdelt, "CRPIX1": float(pix_target_x), "CUNIT1": "deg",
            "CTYPE2": "OFFSETY", "CRVAL2": 0.0, "CDELT2": cdelt, "CRPIX2": float(pix_target_y), "CUNIT2": "deg",
            "WCSNAME": "MOON-DISK",
            "COORDSYS": "MOON_DISK_OFFSET",
            "OFFSYS": "MOON_DISK",
            "OFFBODY": "MOON",
            "OFFREF": "BODYCTR",
            "OFFXDEF": "PA_NORTH_PLUS_90",
            "OFFYDEF": "LUNAR_NORTH",
            "DISKREF": "APPARENT",
            "DISKYPOS": "LUNAR_NORTH",
            "DISKXPOS": "PA_NORTH_PLUS_90",
            "PAPOLE": "ROW_DEP",
            "PROJTYPE": proj_fits,
            "REFOFFX": (0.0, "Reference offset X (deg)"),
            "REFOFFY": (0.0, "Reference offset Y (deg)"),
        }

    # Moving-object sky-plane offsets are centered on the apparent body center,
    # but the axes are equatorial tangent-plane offsets rather than a fixed
    # RA/Dec image.  Do not fake RA---/DEC--; each row used a different body
    # center.  A linear OFFSET WCS keeps generic FITS viewers interoperable.
    if is_sky_offset:
        cdelt = float(cell_arcsec) / 3600.0
        pix_target_x = 1.0 - (float(config_x0) / float(cell_arcsec))
        pix_target_y = 1.0 - (float(config_y0) / float(cell_arcsec))
        body_name = "" if sky_body is None else sky_body.upper()
        return {
            "CTYPE1": "OFFSETX", "CRVAL1": 0.0, "CDELT1": cdelt, "CRPIX1": float(pix_target_x), "CUNIT1": "deg",
            "CTYPE2": "OFFSETY", "CRVAL2": 0.0, "CDELT2": cdelt, "CRPIX2": float(pix_target_y), "CUNIT2": "deg",
            "WCSNAME": "BODY-SKY",
            "COORDSYS": "SKY_OFFSET",
            "OFFSYS": "RADEC",
            "OFFBODY": body_name,
            "OFFREF": "BODYCTR",
            "OFFXDEF": "dRA*cos(Dec_body)" if proj_fits == "SFL" else "dRA",
            "OFFYDEF": "dDec",
            "PROJTYPE": proj_fits,
            "REFOFFX": (0.0, "Reference offset X (deg)"),
            "REFOFFY": (0.0, "Reference offset Y (deg)"),
        }

    # Moving-object dAz/dEl maps are not fixed celestial images.  Do not fake
    # RA---TAN/DEC--TAN.  Use a valid FITS linear WCS with units of degrees so
    # generic FITS viewers can display it, while custom keywords document the
    # physical meaning.
    if is_offset:
        cdelt = float(cell_arcsec) / 3600.0
        pix_target_x = 1.0 - (float(config_x0) / float(cell_arcsec))
        pix_target_y = 1.0 - (float(config_y0) / float(cell_arcsec))
        body_name = "" if body is None else body.upper()
        return {
            "CTYPE1": "OFFSETX", "CRVAL1": 0.0, "CDELT1": cdelt, "CRPIX1": float(pix_target_x), "CUNIT1": "deg",
            "CTYPE2": "OFFSETY", "CRVAL2": 0.0, "CDELT2": cdelt, "CRPIX2": float(pix_target_y), "CUNIT2": "deg",
            "WCSNAME": "AZEL-OFF",
            "COORDSYS": "AZEL_OFFSET",
            "OFFSYS": "AZEL",
            "OFFBODY": body_name,
            "OFFREF": "BODYCTR",
            "OFFXDEF": "dAz*cos(El_body)" if proj_fits == "SFL" else "dAz",
            "OFFYDEF": "dEl",
            "PROJTYPE": proj_fits,
            "REFOFFX": (0.0, "Reference offset X (deg)"),
            "REFOFFY": (0.0, "Reference offset Y (deg)"),
        }
        
    if sys_upper in ("RADEC", "ICRS"):
        ctype1, ctype2 = f"RA---{proj_fits}", f"DEC--{proj_fits}"
    elif sys_upper in ("GALACTIC", "GAL"):
        ctype1, ctype2 = f"GLON-{proj_fits}", f"GLAT-{proj_fits}"
    else:
        ctype1, ctype2 = f"LON--{proj_fits}", f"LAT--{proj_fits}"

    cdelt1 = -cell_arcsec / 3600.0
    cdelt2 = cell_arcsec / 3600.0

    # ターゲット座標(lon0, lat0)が、局所平面で x=0, y=0 になるFITSのピクセル位置(1-based)
    pix_target_x = 1.0 - (config_x0 / cell_arcsec)
    pix_target_y = 1.0 - (config_y0 / cell_arcsec)

    if proj_fits == "SFL":
        crval1 = float(lon0)
        crval2 = 0.0  # SFLの鉄則：基準は必ず赤道に固定する
        crpix1 = float(pix_target_x)
        # CRVAL2が0なので、ターゲット緯度(lat0)が正しい位置に来るようにCRPIX2をシフトする
        crpix2 = float(pix_target_y - (lat0 / cdelt2))
    else:
        crval1 = float(lon0)
        crval2 = float(lat0)
        crpix1 = float(pix_target_x)
        crpix2 = float(pix_target_y)

    return {
        "CTYPE1": ctype1, "CRVAL1": crval1, "CDELT1": cdelt1, "CRPIX1": crpix1, "CUNIT1": "deg",
        "CTYPE2": ctype2, "CRVAL2": crval2, "CDELT2": cdelt2, "CRPIX2": crpix2, "CUNIT2": "deg",
    }
