"""Multi-beam beam measurement package.

The package keeps imports intentionally light at module-import time so that
fit-only / config-only workflows do not require the full NECST runtime stack.
Heavy analysis modules are imported lazily by the CLI wrappers and the public
API helpers.
"""

from .sunscan_config import (
    InputConfig,
    CalibrationConfig,
    RippleConfig,
    ProfileConfig,
    TrimConfig,
    EdgeFitConfig,
    ReportConfig,
    RuntimeConfig,
    SunScanAnalysisConfig,
)

__all__ = [
    "InputConfig",
    "CalibrationConfig",
    "RippleConfig",
    "ProfileConfig",
    "TrimConfig",
    "EdgeFitConfig",
    "ReportConfig",
    "RuntimeConfig",
    "SunScanAnalysisConfig",
    "run_singlebeam",
    "run_multibeam_extract",
    "run_multibeam_fit",
    "run_pseudo_multibeam",
    "check_spectrometer_config",
]

__version__ = "v3.8"


def __getattr__(name):
    if name in {
        "run_singlebeam",
        "run_multibeam_extract",
        "run_multibeam_fit",
        "run_pseudo_multibeam",
        "check_spectrometer_config",
    }:
        from .public_api import (
            run_singlebeam,
            run_multibeam_extract,
            run_multibeam_fit,
            run_pseudo_multibeam,
            check_spectrometer_config,
        )
        return {
            "run_singlebeam": run_singlebeam,
            "run_multibeam_extract": run_multibeam_extract,
            "run_multibeam_fit": run_multibeam_fit,
            "run_pseudo_multibeam": run_pseudo_multibeam,
            "check_spectrometer_config": check_spectrometer_config,
        }[name]
    raise AttributeError(name)
