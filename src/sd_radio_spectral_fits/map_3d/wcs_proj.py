import numpy as np

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

def build_spatial_wcs_dict(
    coord_sys: str, projection: str, 
    lon0: float, lat0: float, 
    config_x0: float, config_y0: float, 
    cell_arcsec: float, 
    nx: int, ny: int
) -> dict:
    """FITS出力用の空間WCSヘッダ辞書を構築する (CAR/SFL 厳密対応版)"""
    sys_upper = coord_sys.upper()
    proj_upper = projection.upper()
    
    # 【修正】未知の投影法は弾く（CARへの無言フォールバックをやめる）
    if proj_upper in ("GLS", "SFL", "SINE"):
        proj_fits = "SFL"
    elif proj_upper in ("CAR", "PLATE_CARREE", "NONE"):
        proj_fits = "CAR"
    else:
        raise ValueError(f"Unknown projection for WCS header: {proj_upper}. Supported: CAR, SFL.")
        
    
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
