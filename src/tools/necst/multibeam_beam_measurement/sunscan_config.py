from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import copy
from pathlib import Path
from typing import Any, Dict, Optional
import hashlib
import json


DEFAULT_DB_NAMESPACE = "necst"
DEFAULT_TELESCOPE = "OMU1P85M"
DEFAULT_TEL_LOADDATA = "OMU1p85m"
DEFAULT_PLANET = "sun"
DEFAULT_SPECTRAL_NAME = "xffts-board1"
DEFAULT_OUTDIR = "."
DEFAULT_AZEL_SOURCE = "encoder"
DEFAULT_AZEL_CORRECTION_APPLY = None
DEFAULT_ALTAZ_APPLY = DEFAULT_AZEL_CORRECTION_APPLY
DEFAULT_ENCODER_TABLE_SUFFIX = "ctrl-antenna-encoder"
DEFAULT_ALTAZ_TABLE_SUFFIX = "ctrl-antenna-altaz"
DEFAULT_WEATHER_INSIDE_TABLE_SUFFIX = "weather-ambient"
DEFAULT_WEATHER_OUTSIDE_TABLE_SUFFIX = "weather-ambient"
DEFAULT_TIME_COL = "time"
DEFAULT_SPECTROMETER_TIME_OFFSET_SEC = 0.0
DEFAULT_ENCODER_SHIFT_SEC = 0.0
DEFAULT_ENCODER_AZ_TIME_OFFSET_SEC = 0.0
DEFAULT_ENCODER_EL_TIME_OFFSET_SEC = 0.0
DEFAULT_ENCODER_VAVG_SEC = 0.0
DEFAULT_CHOPPER_WHEEL = True
DEFAULT_TAMB_FALLBACK_K = 300.0
DEFAULT_CHOPPER_WIN_SEC = 5.0
DEFAULT_CHOPPER_STAT = "mean"
DEFAULT_RIPPLE_REMOVE = True
DEFAULT_RIPPLE_MODEL = "auto"
DEFAULT_RIPPLE_TARGET_HZ = 1.2
DEFAULT_RIPPLE_SEARCH_HZ = 0.3
DEFAULT_RIPPLE_BW_HZ = 0.7
DEFAULT_RIPPLE_MAX_HARM = 8
DEFAULT_RIPPLE_ORDER = 2
DEFAULT_RIPPLE_NOTCH_PASS = 4
DEFAULT_RIPPLE_TREND_WIN_SEC = 5.0
DEFAULT_RIPPLE_RESAMPLE_DT_SEC = 0.0
DEFAULT_RIPPLE_PRESET = "auto"
DEFAULT_PROFILE_XLIM_DEG = 1.0
DEFAULT_DEBUG_PLOT = False
DEFAULT_RIPPLE_EVAL_BAND_HZ = 0.0
DEFAULT_TRIM_SCAN = True
DEFAULT_TRIM_VFRAC = 0.20
DEFAULT_TRIM_VMIN = 1e-4
DEFAULT_TRIM_GAP = 10
DEFAULT_TRIM_MIN_SAMPLES = 100
DEFAULT_TRIM_DOMINANT_AXIS = True
DEFAULT_TRIM_AXIS_RATIO_MIN = 3.0
DEFAULT_TRIM_VPERCENTILE = 95.0
DEFAULT_TRIM_STEADY_SCAN = True
DEFAULT_TRIM_USE_ON_ONLY = True
DEFAULT_TRIM_SCAN_SPEED_MIN_ARCSEC_S = 20.0
DEFAULT_TRIM_XWIN_FACTOR = 1.2
DEFAULT_TRIM_CROSS_OFFSET_MAX_DEG = 0.5
DEFAULT_TRIM_STEADY_CV_MAX = 0.8
DEFAULT_STRICT_DERIV = True
DEFAULT_CONTINUE_ON_ERROR = False
DEFAULT_EDGE_FIT = True
DEFAULT_EDGE_FIT_WIN_DEG = 0.15
DEFAULT_EDGE_FIT_THRESHOLD = 0.20
DEFAULT_HPBW_INIT_ARCSEC = 324.0
DEFAULT_EDGE_FIT_PLOT_MAX_SCANS = 3
DEFAULT_HPBW_FACTOR = 1.2
DEFAULT_DISH_DIAMETER_M = 1.85


@dataclass
class InputConfig:
    rawdata_path: Path
    db_namespace: str = DEFAULT_DB_NAMESPACE
    telescope: str = DEFAULT_TELESCOPE
    tel_loaddata: str = DEFAULT_TEL_LOADDATA
    planet: str = DEFAULT_PLANET
    spectral_name: str = DEFAULT_SPECTRAL_NAME
    azel_source: str = DEFAULT_AZEL_SOURCE
    azel_correction_apply: Optional[str] = DEFAULT_AZEL_CORRECTION_APPLY
    altaz_apply: Optional[str] = DEFAULT_ALTAZ_APPLY
    encoder_table: Optional[str] = None
    encoder_table_suffix: str = DEFAULT_ENCODER_TABLE_SUFFIX
    altaz_table: Optional[str] = None
    altaz_table_suffix: str = DEFAULT_ALTAZ_TABLE_SUFFIX
    encoder_time_col: str = DEFAULT_TIME_COL
    altaz_time_col: str = DEFAULT_TIME_COL
    spectrometer_time_offset_sec: float = DEFAULT_SPECTROMETER_TIME_OFFSET_SEC
    encoder_shift_sec: float = DEFAULT_ENCODER_SHIFT_SEC
    encoder_az_time_offset_sec: float = DEFAULT_ENCODER_AZ_TIME_OFFSET_SEC
    encoder_el_time_offset_sec: float = DEFAULT_ENCODER_EL_TIME_OFFSET_SEC
    encoder_vavg_sec: float = DEFAULT_ENCODER_VAVG_SEC


@dataclass
class CalibrationConfig:
    chopper_wheel: bool = DEFAULT_CHOPPER_WHEEL
    tamb_k: Optional[float] = None
    tamb_default_k: float = DEFAULT_TAMB_FALLBACK_K
    tamb_min_k: float = 250.0
    tamb_max_k: float = 330.0
    weather_inside_table: Optional[str] = None
    weather_inside_table_suffix: str = DEFAULT_WEATHER_INSIDE_TABLE_SUFFIX
    weather_inside_time_col: str = DEFAULT_TIME_COL
    chopper_win_sec: float = DEFAULT_CHOPPER_WIN_SEC
    chopper_stat: str = DEFAULT_CHOPPER_STAT


@dataclass
class RefractionConfig:
    weather_outside_table: Optional[str] = None
    weather_outside_table_suffix: str = DEFAULT_WEATHER_OUTSIDE_TABLE_SUFFIX
    weather_outside_time_col: str = DEFAULT_TIME_COL
    outside_default_temperature_c: float = 0.0
    outside_default_pressure_hpa: float = 760.0
    outside_default_humidity_pct: float = 30.0
    outside_temperature_min_c: float = -50.0
    outside_temperature_max_c: float = 50.0
    outside_pressure_min_hpa: float = 400.0
    outside_pressure_max_hpa: float = 1100.0
    outside_humidity_min_pct: float = 0.0
    outside_humidity_max_pct: float = 100.0


@dataclass
class RippleConfig:
    enabled: bool = DEFAULT_RIPPLE_REMOVE
    preset: str = DEFAULT_RIPPLE_PRESET
    model: str = DEFAULT_RIPPLE_MODEL
    target_hz: float = DEFAULT_RIPPLE_TARGET_HZ
    search_hz: float = DEFAULT_RIPPLE_SEARCH_HZ
    bw_hz: Optional[float] = None
    max_harm: Optional[int] = None
    order: Optional[int] = None
    notch_passes: Optional[int] = None
    trend_win_sec: Optional[float] = None
    resample_dt_sec: Optional[float] = None
    eval_band_hz: Optional[float] = None


@dataclass
class ProfileConfig:
    profile_xlim_deg: float = DEFAULT_PROFILE_XLIM_DEG


@dataclass
class TrimConfig:
    enabled: bool = DEFAULT_TRIM_SCAN
    vfrac: float = DEFAULT_TRIM_VFRAC
    vmin: float = DEFAULT_TRIM_VMIN
    gap_fill: int = DEFAULT_TRIM_GAP
    min_samples: int = DEFAULT_TRIM_MIN_SAMPLES
    dominant_axis: bool = DEFAULT_TRIM_DOMINANT_AXIS
    ratio_min: float = DEFAULT_TRIM_AXIS_RATIO_MIN
    vpercentile: float = DEFAULT_TRIM_VPERCENTILE
    steady_scan: bool = DEFAULT_TRIM_STEADY_SCAN
    use_on_only: bool = DEFAULT_TRIM_USE_ON_ONLY
    xwin_factor: float = DEFAULT_TRIM_XWIN_FACTOR
    cross_offset_max_deg: float = DEFAULT_TRIM_CROSS_OFFSET_MAX_DEG
    speed_min_deg_s: float = DEFAULT_TRIM_SCAN_SPEED_MIN_ARCSEC_S / 3600.0
    steady_cv_max: float = DEFAULT_TRIM_STEADY_CV_MAX


@dataclass
class EdgeFitConfig:
    enabled: bool = DEFAULT_EDGE_FIT
    strict_deriv: bool = DEFAULT_STRICT_DERIV
    fit_win_deg: float = DEFAULT_EDGE_FIT_WIN_DEG
    fit_threshold: float = DEFAULT_EDGE_FIT_THRESHOLD
    hpbw_init_arcsec: float = DEFAULT_HPBW_INIT_ARCSEC


@dataclass
class ReportConfig:
    outdir: Path = Path(DEFAULT_OUTDIR)
    debug_plot: bool = DEFAULT_DEBUG_PLOT
    edge_fit_plot_max_scans: int = DEFAULT_EDGE_FIT_PLOT_MAX_SCANS
    tag: Optional[str] = None


@dataclass
class RuntimeConfig:
    continue_on_error: bool = DEFAULT_CONTINUE_ON_ERROR
    hpbw_init_explicit: bool = False
    azel_source_explicit: bool = False
    azel_correction_apply_explicit: bool = False


@dataclass
class BeamOverride:
    stream_name: Optional[str] = None
    beam_id: Optional[str] = None
    restfreq_hz: Optional[float] = None
    hpbw_init_arcsec: Optional[float] = None
    polariza: Optional[str] = None
    fdnum: Optional[int] = None
    ifnum: Optional[int] = None
    plnum: Optional[int] = None
    sampler: Optional[str] = None


@dataclass
class SunScanAnalysisConfig:
    input: InputConfig
    calibration: CalibrationConfig = field(default_factory=CalibrationConfig)
    refraction: RefractionConfig = field(default_factory=RefractionConfig)
    ripple: RippleConfig = field(default_factory=RippleConfig)
    profile: ProfileConfig = field(default_factory=ProfileConfig)
    trim: TrimConfig = field(default_factory=TrimConfig)
    edge_fit: EdgeFitConfig = field(default_factory=EdgeFitConfig)
    report: ReportConfig = field(default_factory=ReportConfig)
    runtime: RuntimeConfig = field(default_factory=RuntimeConfig)
    beam_override: BeamOverride = field(default_factory=BeamOverride)
    dish_diameter_m: float = DEFAULT_DISH_DIAMETER_M
    hpbw_factor: float = DEFAULT_HPBW_FACTOR

    @classmethod
    def default(cls, rawdata_path: Path, spectral_name: str = DEFAULT_SPECTRAL_NAME, outdir: Path | str = DEFAULT_OUTDIR) -> "SunScanAnalysisConfig":
        rawdata_path = Path(rawdata_path)
        outdir_path = Path(outdir)
        return cls(
            input=InputConfig(rawdata_path=rawdata_path, spectral_name=spectral_name),
            report=ReportConfig(outdir=outdir_path),
        )

    def resolved_tag(self) -> str:
        return str(self.report.tag) if self.report.tag else self.input.rawdata_path.name

    def with_updates(self, **kwargs: Any) -> "SunScanAnalysisConfig":
        return replace(self, **kwargs)

    def with_stream_override(
        self,
        *,
        spectral_name: Optional[str] = None,
        stream_name: Optional[str] = None,
        beam_id: Optional[str] = None,
        restfreq_hz: Optional[float] = None,
        hpbw_init_arcsec: Optional[float] = None,
        polariza: Optional[str] = None,
        fdnum: Optional[int] = None,
        ifnum: Optional[int] = None,
        plnum: Optional[int] = None,
        sampler: Optional[str] = None,
        outdir: Optional[Path] = None,
        tag: Optional[str] = None,
    ) -> "SunScanAnalysisConfig":
        new_cfg = copy.deepcopy(self)
        if spectral_name is not None:
            new_cfg.input.spectral_name = spectral_name
        if outdir is not None:
            new_cfg.report.outdir = Path(outdir)
        if tag is not None:
            new_cfg.report.tag = tag
        if hpbw_init_arcsec is not None:
            new_cfg.edge_fit.hpbw_init_arcsec = float(hpbw_init_arcsec)
            new_cfg.beam_override.hpbw_init_arcsec = hpbw_init_arcsec
        if stream_name is not None:
            new_cfg.beam_override.stream_name = stream_name
        if beam_id is not None:
            new_cfg.beam_override.beam_id = beam_id
        if restfreq_hz is not None:
            new_cfg.beam_override.restfreq_hz = restfreq_hz
        if polariza is not None:
            new_cfg.beam_override.polariza = polariza
        if fdnum is not None:
            new_cfg.beam_override.fdnum = fdnum
        if ifnum is not None:
            new_cfg.beam_override.ifnum = ifnum
        if plnum is not None:
            new_cfg.beam_override.plnum = plnum
        if sampler is not None:
            new_cfg.beam_override.sampler = sampler
        return new_cfg

    def config_digest(self) -> str:
        def _normalize(value: Any) -> Any:
            if isinstance(value, Path):
                return str(value)
            if isinstance(value, dict):
                return {str(k): _normalize(v) for k, v in value.items()}
            if isinstance(value, (list, tuple)):
                return [_normalize(v) for v in value]
            return value

        payload: Dict[str, Any] = _normalize(asdict(self))
        blob = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(blob).hexdigest()[:16]
