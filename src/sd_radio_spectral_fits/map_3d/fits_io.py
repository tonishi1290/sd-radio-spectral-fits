# src/map/fits_io.py
import numpy as np
from astropy.io import fits
from .wcs_proj import build_spatial_wcs_dict


def save_map_fits(
    grid_res,
    v_tgt,
    coord_sys,
    projection,
    lon0,
    lat0,
    config,
    out_path,
    out_scale="TA*",
    rep_beameff=1.0,
):
    """GridResult を 3D Multi-Extension FITS として書き出す。"""
    ny, nx, nchan = grid_res.cube.shape
    wcs_hdr = build_spatial_wcs_dict(
        coord_sys, projection, lon0, lat0, config.x0, config.y0, config.cell_arcsec, nx, ny
    )

    header = fits.Header()
    header.update(wcs_hdr)

    _fill_header_metadata(header, coord_sys, lon0, lat0, out_scale, rep_beameff, grid_res)

    header["CTYPE3"] = "VRAD"
    header["CUNIT3"] = "km/s"
    header["CRPIX3"] = 1.0
    header["CRVAL3"] = float(v_tgt[0])
    header["CDELT3"] = float(v_tgt[1] - v_tgt[0]) if len(v_tgt) > 1 else 1.0

    if getattr(grid_res, "meta", None):
        if grid_res.meta.get("RESTFREQ", 0) > 0:
            header["RESTFRQ"] = (grid_res.meta["RESTFREQ"], "Rest Frequency (Hz)")
        if "SPECSYS" in grid_res.meta:
            header["SPECSYS"] = (grid_res.meta["SPECSYS"], "Spectral reference frame")
        header["VELDEF"] = ("RADI-LSR", "Radio velocity definition")

    cube_fits = np.transpose(grid_res.cube, (2, 0, 1))
    hdu_primary = fits.PrimaryHDU(data=cube_fits.astype(np.float32), header=header)
    hdul = fits.HDUList([hdu_primary])

    _add_diagnostic_hdus(hdul, grid_res, header)

    hdul.writeto(out_path, overwrite=True)


def _fill_header_metadata(header, coord_sys, lon0, lat0, out_scale, rep_beameff, grid_res):
    """ヘッダーへの観測情報の書き込み"""
    is_galactic = coord_sys.lower() in ("galactic", "gal")
    if is_galactic:
        header["OBSGLON"] = (float(lon0), "Ref Galactic Longitude (deg)")
        header["OBSGLAT"] = (float(lat0), "Ref Galactic Latitude (deg)")
    else:
        header["OBSRA"] = (float(lon0), "Ref Right Ascension (deg)")
        header["OBSDEC"] = (float(lat0), "Ref Declination (deg)")

    # FITS BUNIT is kept machine-readable; the temperature scale is split out.
    header["BUNIT"] = ("K", "Unit of data")
    header["TEMPSCAL"] = (str(out_scale), "Temperature scale of primary cube")
    header["BTYPE"] = ("BrightnessTemperature", "Physical type of primary cube")
    if np.isfinite(rep_beameff):
        header["BEAMEFF"] = (float(rep_beameff), "Main beam efficiency")

    if getattr(grid_res, "meta", None):
        meta = grid_res.meta
        header["BMAJ"] = (float(meta.get("bmaj_eff_arcsec", 0)) / 3600.0, "deg")


def _add_diagnostic_hdus(hdul, grid_res, base_header):
    """2D診断マップ（HDU）の動的追加（漏れなく全マップを追加）"""
    header_2d = base_header.copy()
    for k in ["CTYPE3", "CUNIT3", "CRPIX3", "CRVAL3", "CDELT3"]:
        if k in header_2d:
            del header_2d[k]

    maps = [
        ("WEIGHT", grid_res.weight_map, "Weight", ""),
        ("HIT", grid_res.hit_map, "HitCount", ""),
        ("MASK", grid_res.mask_map.astype(np.float32), "ValidMask", ""),
        ("TSYS", grid_res.tsys_map, "SystemTemp", "K"),
        ("TINT", grid_res.tint_map, "IntegrationTime", "s"),
        ("TIME", grid_res.time_map, "ObservationTime", "MJD"),
        ("RMS", grid_res.rms_map, "BaselineRMS", "K"),
        ("NEFF", getattr(grid_res, "neff_map", None), "EffectiveSamples", ""),
        ("XEFF", getattr(grid_res, "xeff_map", None), "EffectiveX", "arcsec"),
        ("YEFF", getattr(grid_res, "yeff_map", None), "EffectiveY", "arcsec"),
        ("BIAS_PIX", getattr(grid_res, "dr_eff_map_pix", None), "PositionBias", "pixel"),
    ]

    for name, data, btype, bunit in maps:
        if data is not None:
            h = header_2d.copy()
            h["BTYPE"], h["BUNIT"] = btype, bunit
            hdul.append(fits.ImageHDU(data=data.astype(np.float32), header=h, name=name))
