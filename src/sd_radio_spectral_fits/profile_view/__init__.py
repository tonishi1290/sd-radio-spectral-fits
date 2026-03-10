# src/sd_radio_spectral_fits/profile_view/__init__.py
from __future__ import annotations
from typing import Union

from ..fitsio import Scantable, read_scantable

# viewer.py から view_spectra 関数と SpectralViewer クラスを読み込む
from .viewer import view_spectra, SpectralViewer

# grid.py からはクラスだけ読み込む（関数はここで定義するため）
from .grid import ProfileMapGridViewer
from .montage import ProfileMapMontageViewer

def plot_profile_map(input_data: Union[str, Scantable], **kwargs):
    """
    Launch profile map viewer.
    
    Parameters
    ----------
    input_data : str or Scantable
    mode : "montage" (default) or "grid"
    projection : "TAN" (default) or "GLS" (for grid mode)
    """
    st = read_scantable(input_data) if isinstance(input_data, str) else input_data
    
    # alias checks
    if "axis_type" in kwargs and "xaxis" not in kwargs: kwargs["xaxis"] = kwargs.pop("axis_type")
    if "xranges" in kwargs and "xrange" not in kwargs: kwargs["xrange"] = kwargs.pop("xranges")
    if "yranges" in kwargs and "yrange" not in kwargs: kwargs["yrange"] = kwargs.pop("yranges")

    mode = str(kwargs.pop("mode", "montage")).lower()
    if mode in ("montage", "mosaic"):
        return ProfileMapMontageViewer(st, **kwargs)
    if mode in ("grid", "gridding"):
        return ProfileMapGridViewer(st, **kwargs)

    raise ValueError("mode must be 'montage' or 'grid'")

__all__ = [
    "view_spectra",
    "plot_profile_map",
    "SpectralViewer",
    "ProfileMapMontageViewer",
    "ProfileMapGridViewer",
]
