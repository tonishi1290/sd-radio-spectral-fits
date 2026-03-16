from __future__ import annotations

from typing import Dict, Tuple
import pandas as pd

from . import sunscan_legacy_compat as compat
from .sunscan_config import SunScanAnalysisConfig


def build_dataframe_from_config(config: SunScanAnalysisConfig):
    return compat.build_dataframe(
        config.input.rawdata_path,
        config.input.spectral_name,
        db_namespace=str(config.input.db_namespace),
        telescope=str(config.input.telescope),
        tel_loaddata=str(config.input.tel_loaddata),
        planet=str(config.input.planet),
        azel_source=str(config.input.azel_source),
        altaz_apply=str(config.input.altaz_apply),
        spectrometer_time_offset_sec=float(config.input.spectrometer_time_offset_sec),
        encoder_shift_sec=float(config.input.encoder_shift_sec),
        encoder_az_time_offset_sec=float(config.input.encoder_az_time_offset_sec),
        encoder_el_time_offset_sec=float(config.input.encoder_el_time_offset_sec),
        encoder_vavg_sec=float(config.input.encoder_vavg_sec),
        chopper_wheel=bool(config.calibration.chopper_wheel),
        tamb_k=config.calibration.tamb_k,
        chopper_win_sec=float(config.calibration.chopper_win_sec),
        chopper_stat=str(config.calibration.chopper_stat),
    )


def dataframe_bundle_from_config(config: SunScanAnalysisConfig) -> Dict[str, object]:
    df, az_scans, el_scans = build_dataframe_from_config(config)
    return {
        "df": df,
        "az_scans": az_scans,
        "el_scans": el_scans,
    }
