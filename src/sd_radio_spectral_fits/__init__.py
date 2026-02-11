"""sd_radio_spectral_fits

Single-dish radio spectral SDFITS writer.

Public API is intentionally small. For advanced/low-level symbols,
import from :mod:`sd_radio_spectral_fits.sdfits_writer`.
"""

from .version import __version__, SWNAME
from .sdfits_writer import (
    Site,
    Efficiency,
    SpectralAxisUniform,
    DatasetInfo,
    SDRadioSpectralSDFITSWriter,
)

__all__ = [
    "__version__",
    "SWNAME",
    "Site",
    "Efficiency",
    "SpectralAxisUniform",
    "DatasetInfo",
    "SDRadioSpectralSDFITSWriter",
]
