"""Compatibility import path.

This module is kept for older code that used:
    from sd_radio_spectral_fits.writer import ...

New code should prefer:
    from sd_radio_spectral_fits import SDRadioSpectralSDFITSWriter
"""

from .sdfits_writer import *  # noqa: F401,F403
