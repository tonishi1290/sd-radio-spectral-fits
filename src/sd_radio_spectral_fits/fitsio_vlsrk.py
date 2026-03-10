from __future__ import annotations

import builtins
from typing import Optional

import numpy as np
import pandas as pd

def _require_astropy():
    try:
        from astropy.io import fits  # noqa
        return fits
    except Exception as e:
        raise ImportError("This function requires astropy. Install: pip install 'sdrsf-pipeline[fits]'") from e

def write_vlsrk_fits(path: str, *, meta: dict, data: np.ndarray, table: pd.DataFrame, history: Optional[dict] = None) -> None:
    """Write spectra on a common VLSRK velocity grid.

    Header WCS (axis 1):
      CTYPE1='VRAD', CUNIT1='km/s', CRVAL1/CRPIX1/CDELT1 defined on VLSRK axis.
      RESTFRQ is kept (Hz). SPECSYS='LSRK'.
    """
    fits = _require_astropy()
    hdr = fits.Header()
    hdr["CTYPE1"] = str(meta.get("CTYPE1", "VRAD"))
    hdr["CUNIT1"] = str(meta.get("CUNIT1", "km/s"))
    hdr["CRVAL1"] = float(meta["CRVAL1"])
    hdr["CRPIX1"] = float(meta.get("CRPIX1", 1.0))
    hdr["CDELT1"] = float(meta["CDELT1"])
    hdr["RESTFRQ"] = float(meta["RESTFREQ"])
    hdr["SPECSYS"] = str(meta.get("SPECSYS", "LSRK"))
    hdr["TIMESYS"] = str(meta.get("TIMESYS", "UTC"))
    hdr["BUNIT"] = str(meta.get("BUNIT", "K"))
    hdr["NAXIS1"] = int(meta.get("NAXIS1", data.shape[1]))

    if "TELESCOPE" in meta and meta["TELESCOPE"]:
        hdr["TELESCOP"] = str(meta["TELESCOPE"])
    if "BACKEND" in meta and meta["BACKEND"]:
        hdr["BACKEND"] = str(meta["BACKEND"])

    hdu_data = fits.ImageHDU(data.astype(np.float32), header=hdr, name="DATA")

    tab = table.copy()
    cols = []
    for col in tab.columns:
        arr = tab[col].to_numpy()
        if arr.dtype.kind in "if":
            fmt = "D"
        elif arr.dtype.kind in "i":
            fmt = "K"
        else:
            arr = arr.astype(str)
            mx = int(np.max([len(s) for s in arr]) if arr.size else 1)
            fmt = f"A{builtins.max(1, mx)}"
        cols.append(fits.Column(name=str(col).upper(), format=fmt, array=arr))
    hdu_tab = fits.BinTableHDU.from_columns(cols, name="POINTS")

    hdus = [fits.PrimaryHDU(), hdu_data, hdu_tab]

    if history is not None:
        import json
        js = json.dumps(history, ensure_ascii=False)
        col = fits.Column(name="HISTORY_JSON", format=f"A{builtins.max(1,len(js))}", array=np.asarray([js], dtype=object))
        hdu_hist = fits.BinTableHDU.from_columns([col], name="HISTORY")
        hdus.append(hdu_hist)

    fits.HDUList(hdus).writeto(path, overwrite=True)
