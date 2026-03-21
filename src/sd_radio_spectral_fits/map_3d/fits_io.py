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
        if "reproducible_mode" in grid_res.meta:
            header["REPRODUC"] = (bool(grid_res.meta["reproducible_mode"]), "Reproducible mode")
        if "gridding_workers" in grid_res.meta:
            header["GRIDWORK"] = (int(grid_res.meta["gridding_workers"]), "KDTree query workers")
        if "sorted_neighbors" in grid_res.meta:
            header["SORTNBR"] = (bool(grid_res.meta["sorted_neighbors"]), "Neighbor sort fixed")

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

    header["BUNIT"] = ("K", "Unit of data")
    header["TEMPSCAL"] = (str(out_scale), "Temperature scale of primary cube")
    header["BTYPE"] = ("BrightnessTemperature", "Physical type of primary cube")
    if np.isfinite(rep_beameff):
        header["BEAMEFF"] = (float(rep_beameff), "Main beam efficiency")

    if getattr(grid_res, "meta", None):
        meta = grid_res.meta
        if "kernel" in meta:
            header["KERNEL"] = (str(meta["kernel"]), "Gridding kernel")
        if "kernel_preset" in meta:
            header["KPRESET"] = (str(meta["kernel_preset"]), "Kernel default preset")
        if "kernel_sign" in meta:
            header["KSIGN"] = (str(meta["kernel_sign"]), "Kernel sign handling mode")
        if np.isfinite(meta.get("gwidth_arcsec", np.nan)):
            header["GWIDTH"] = (float(meta["gwidth_arcsec"]) / 3600.0, "Gaussian HWHM [deg]")
        if np.isfinite(meta.get("jwidth_arcsec", np.nan)):
            header["JWIDTH"] = (float(meta["jwidth_arcsec"]) / 3600.0, "Jinc width [deg]")
        if np.isfinite(meta.get("support_radius_arcsec", np.nan)):
            header["SUPPORTR"] = (float(meta["support_radius_arcsec"]) / 3600.0, "Kernel support radius [deg]")
        if "weight_map_semantics" in meta:
            header["WGTSEM"] = (str(meta["weight_map_semantics"]), "WEIGHT semantics")
        if "hit_map_semantics" in meta:
            header["HITSEM"] = (str(meta["hit_map_semantics"]), "HIT semantics")
        if "nsamp_map_semantics" in meta:
            header["NSAMPSEM"] = (str(meta["nsamp_map_semantics"]), "NSAMP semantics")
        if "mask_map_semantics" in meta:
            header["MASKSEM"] = (str(meta["mask_map_semantics"]), "MASK semantics")
        if "wsum_map_semantics" in meta:
            header["WSUMSEM"] = (str(meta["wsum_map_semantics"]), "WSUM semantics")
        if "wabs_map_semantics" in meta:
            header["WABSSEM"] = (str(meta["wabs_map_semantics"]), "WABS semantics")
        if "cancel_map_semantics" in meta:
            header["CANCSEM"] = (str(meta["cancel_map_semantics"]), "CANCEL semantics")
        if "weight_rel_map_semantics" in meta:
            header["WRELSEM"] = (str(meta["weight_rel_map_semantics"]), "WREL semantics")
        if np.isfinite(meta.get("cell_over_beam", np.nan)):
            header["CELLBEAM"] = (float(meta["cell_over_beam"]), "Cell / beam FWHM")
        if "cell_is_coarse" in meta:
            header["CELLCOAR"] = (bool(meta["cell_is_coarse"]), "cell > beam/3 warning")
        if np.isfinite(meta.get("median_abs_weight", np.nan)):
            header["MEDAWGT"] = (float(meta["median_abs_weight"]), "Median |WEIGHT|")
        if "min_abs_weight_ratio" in meta:
            header["MINWREL"] = (float(meta["min_abs_weight_ratio"]), "Min |WEIGHT|/median")
        if "min_cancel_ratio" in meta:
            header["MINCANC"] = (float(meta["min_cancel_ratio"]), "Min cancel ratio")

        bmaj = meta.get("bmaj_empirical_arcsec", meta.get("bmaj_nominal_arcsec", meta.get("bmaj_eff_arcsec", np.nan)))
        bmin = meta.get("bmin_empirical_arcsec", meta.get("bmin_nominal_arcsec", bmaj))
        bpa = meta.get("bpa_empirical_deg", meta.get("bpa_nominal_deg", 0.0))
        if np.isfinite(bmaj):
            header["BMAJ"] = (float(bmaj) / 3600.0, "Eff beam major [deg]")
        if np.isfinite(bmin):
            header["BMIN"] = (float(bmin) / 3600.0, "Eff beam minor [deg]")
        if np.isfinite(bpa):
            header["BPA"] = (float(bpa), "Eff beam PA [deg]")

        if np.isfinite(meta.get("bmaj_nominal_arcsec", np.nan)):
            header["BMAJNOM"] = (float(meta["bmaj_nominal_arcsec"]) / 3600.0, "Nom beam major [deg]")
        if np.isfinite(meta.get("bmin_nominal_arcsec", np.nan)):
            header["BMINNOM"] = (float(meta["bmin_nominal_arcsec"]) / 3600.0, "Nom beam minor [deg]")
        if np.isfinite(meta.get("bpa_nominal_deg", np.nan)):
            header["BPANOM"] = (float(meta["bpa_nominal_deg"]), "Nom beam PA [deg]")
        if np.isfinite(meta.get("nominal_radial_fwhm_arcsec", np.nan)):
            header["BEAMNOM"] = (float(meta["nominal_radial_fwhm_arcsec"]) / 3600.0, "Nom radial FWHM [deg]")

        if np.isfinite(meta.get("bmaj_empirical_arcsec", np.nan)):
            header["BMAJCEN"] = (float(meta["bmaj_empirical_arcsec"]) / 3600.0, "Ctr emp beam maj [deg]")
        if np.isfinite(meta.get("bmin_empirical_arcsec", np.nan)):
            header["BMINCEN"] = (float(meta["bmin_empirical_arcsec"]) / 3600.0, "Ctr emp beam min [deg]")
        if np.isfinite(meta.get("bpa_empirical_deg", np.nan)):
            header["BPACEN"] = (float(meta["bpa_empirical_deg"]), "Ctr emp beam PA [deg]")
        if np.isfinite(meta.get("empirical_radial_fwhm_arcsec", np.nan)):
            header["BEAMCEN"] = (float(meta["empirical_radial_fwhm_arcsec"]) / 3600.0, "Ctr emp radial FWHM [deg]")


def _add_diagnostic_hdus(hdul, grid_res, base_header):
    """2D診断マップ（HDU）の動的追加（漏れなく全マップを追加）"""
    header_2d = base_header.copy()
    for k in ["CTYPE3", "CUNIT3", "CRPIX3", "CRVAL3", "CDELT3"]:
        if k in header_2d:
            del header_2d[k]

    maps = [
        ("WEIGHT", grid_res.weight_map, "GriddingDenominator", ""),
        ("HIT", grid_res.hit_map, "UsedSampleCount", ""),
        ("NSAMP", getattr(grid_res, "nsamp_map", None), "SupportSampleCount", ""),
        ("WSUM", getattr(grid_res, "wsum_map", None), "SignedWeightSum", ""),
        ("WABS", getattr(grid_res, "wabs_map", None), "AbsoluteWeightSum", ""),
        ("CANCEL", getattr(grid_res, "cancel_map", None), "CancellationRatio", ""),
        ("WREL", getattr(grid_res, "weight_rel_map", None), "RelativeAbsWeight", ""),
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
