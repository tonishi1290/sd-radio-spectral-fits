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
        azel_correction_apply=(None if config.input.azel_correction_apply is None else str(config.input.azel_correction_apply)),
        altaz_apply=(None if config.input.altaz_apply is None else str(config.input.altaz_apply)),
        encoder_table=config.input.encoder_table,
        encoder_table_suffix=str(config.input.encoder_table_suffix),
        altaz_table=config.input.altaz_table,
        altaz_table_suffix=str(config.input.altaz_table_suffix),
        encoder_time_col=str(config.input.encoder_time_col),
        altaz_time_col=str(config.input.altaz_time_col),
        spectrometer_time_offset_sec=float(config.input.spectrometer_time_offset_sec),
        encoder_shift_sec=float(config.input.encoder_shift_sec),
        encoder_az_time_offset_sec=float(config.input.encoder_az_time_offset_sec),
        encoder_el_time_offset_sec=float(config.input.encoder_el_time_offset_sec),
        encoder_vavg_sec=float(config.input.encoder_vavg_sec),
        chopper_wheel=bool(config.calibration.chopper_wheel),
        tamb_k=config.calibration.tamb_k,
        tamb_default_k=float(config.calibration.tamb_default_k),
        tamb_min_k=float(config.calibration.tamb_min_k),
        tamb_max_k=float(config.calibration.tamb_max_k),
        weather_inside_table=config.calibration.weather_inside_table,
        weather_inside_table_suffix=str(config.calibration.weather_inside_table_suffix),
        weather_inside_time_col=str(config.calibration.weather_inside_time_col),
        weather_outside_table=config.refraction.weather_outside_table,
        weather_outside_table_suffix=str(config.refraction.weather_outside_table_suffix),
        weather_outside_time_col=str(config.refraction.weather_outside_time_col),
        outside_default_temperature_c=float(config.refraction.outside_default_temperature_c),
        outside_default_pressure_hpa=float(config.refraction.outside_default_pressure_hpa),
        outside_default_humidity_pct=float(config.refraction.outside_default_humidity_pct),
        outside_temperature_min_c=float(config.refraction.outside_temperature_min_c),
        outside_temperature_max_c=float(config.refraction.outside_temperature_max_c),
        outside_pressure_min_hpa=float(config.refraction.outside_pressure_min_hpa),
        outside_pressure_max_hpa=float(config.refraction.outside_pressure_max_hpa),
        outside_humidity_min_pct=float(config.refraction.outside_humidity_min_pct),
        outside_humidity_max_pct=float(config.refraction.outside_humidity_max_pct),
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
