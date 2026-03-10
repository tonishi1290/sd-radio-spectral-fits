# src/sd_radio_spectral_fits/__init__.py
__version__ = "1.0.0"

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


from .fitsio import (
    Scantable,
    read_scantable,
    write_scantable,
)
from .calibrate import (
    run_tastar_calibration,
)
from .coadd import (
    run_velocity_coadd,
)
from .rawspec import (
    load_rawspec_auto,
)

from .scantable_utils import (
    show_scantable,
    describe_columns,
    update_metadata,
    merge_scantables,
    filter_scantable,
    find_scans,
    calc_mapping_offsets,
)

from .atmosphere import (
    get_airmass,
    compute_t_cal_array,
    extract_meta_value,
    estimate_t_atm
)


from .profile_view import view_spectra, plot_profile_map

from .baseline import run_baseline_fit

from .regrid_vlsrk import VGrid, make_vgrid

__all__ = [
    "__version__",
    "SWNAME",
    "Site",
    "Efficiency",
    "SpectralAxisUniform",
    "DatasetInfo",
    "SDRadioSpectralSDFITSWriter",
    "show_scantable",
    "describe_columns",
    "update_metadata",
    "merge_scantables",
    "filter_scantable",
    "find_scans",
    "calc_mapping_offsets",
]
