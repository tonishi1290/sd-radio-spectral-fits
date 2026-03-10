from __future__ import annotations

from typing import Optional, Dict, Any

import pandas as pd

# Speed of light
C_KMS = 299792.458
C_MS = 299792458.0


def _u(s: str) -> str:
    return str(s or "").strip().upper()


def _normalize_unit(u: str) -> str:
    s = str(u or "").strip().lower().replace(" ", "")
    if s in ("hz",):
        return "hz"
    if s in ("km/s", "kms-1", "km.s-1", "km*s-1", "kmsec-1"):
        return "km/s"
    if s in ("m/s", "ms-1", "m.s-1", "m*s-1", "msec-1"):
        return "m/s"
    return s


def _get_restfreq_hz(meta: dict, table: Optional[pd.DataFrame] = None) -> Optional[float]:
    for k in ("RESTFREQ", "RESTFRQ", "rest_hz", "restfreq_hz", "restfrq_hz"):
        v = meta.get(k, None)
        if v not in (None, ""):
            try:
                return float(v)
            except Exception:
                pass
    if table is not None:
        for k in ("RESTFREQ", "RESTFRQ"):
            if k in table.columns:
                s = pd.to_numeric(table[k], errors="coerce")
                if s.notna().any():
                    return float(s.dropna().iloc[0])
    return None


def _is_vrad_axis(meta: dict) -> bool:
    specsys = _u(meta.get("SPECSYS", ""))
    ctype1 = _u(meta.get("CTYPE1", ""))
    return (specsys == "VRAD") or ("VRAD" in ctype1)


def apply_restfreq_override(
    meta: dict,
    table: Optional[pd.DataFrame],
    rest_freq_hz: float,
    *,
    require_wcs_for_vrad: bool = True,
) -> Dict[str, Any]:
    """Apply rest1->rest2 policy in-place to meta (and table columns if present).

    - FREQ axis: keep frequency WCS unchanged; stamp RESTFREQ/RESTFRQ to rest2.
    - VRAD axis (SYSFRAME=VRAD or CTYPE1 contains VRAD): update linear WCS by exact affine transform
        v2 = a*v1 + b
        a = rest1/rest2
        b = c*(1 - rest1/rest2)
      and stamp RESTFREQ/RESTFRQ to rest2.

    Table handling:
      - If table has RESTFREQ/RESTFRQ columns, they are overwritten to rest2.
      - If VRAD axis and table has CRVAL1/CDELT1 columns, they are updated consistently.

    This function does NOT add new columns; it only overwrites existing ones.
    """
    rest2 = float(rest_freq_hz)
    if not (rest2 > 0.0):
        raise ValueError("rest_freq_hz must be > 0")

    info: Dict[str, Any] = {
        "rest1_hz": None,
        "rest2_hz": rest2,
        "axis_vrad": False,
        "wcs_updated": False,
    }

    axis_vrad = _is_vrad_axis(meta)
    info["axis_vrad"] = axis_vrad

    rest1 = _get_restfreq_hz(meta, table)
    if rest1 is not None:
        info["rest1_hz"] = float(rest1)
        # preserve originals once
        meta.setdefault("RESTFREQ_ORIG", float(rest1))
        meta.setdefault("RESTFRQ_ORIG", float(rest1))
    else:
        if axis_vrad:
            raise ValueError("VRAD axis rest conversion requires existing RESTFREQ/RESTFRQ (rest1).")

    # VRAD axis: strict affine WCS update if rest actually changes
    if axis_vrad and (rest1 is not None) and (float(rest1) != rest2):
        # Choose c consistent with units of CRVAL1/CDELT1
        unit = _normalize_unit(meta.get("CUNIT1", "km/s"))
        c = C_MS if unit == "m/s" else C_KMS
        a = float(rest1) / rest2
        b = c * (1.0 - float(rest1) / rest2)

        if require_wcs_for_vrad:
            for k in ("CRVAL1", "CDELT1", "CRPIX1"):
                if meta.get(k, None) in (None, ""):
                    raise ValueError(f"VRAD axis rest conversion requires meta['{k}'].")

        # Update meta WCS if present
        if meta.get("CRVAL1", None) not in (None, ""):
            meta["CRVAL1"] = a * float(meta["CRVAL1"]) + b
        if meta.get("CDELT1", None) not in (None, ""):
            meta["CDELT1"] = a * float(meta["CDELT1"])
        # CRPIX1 unchanged

        # Update table WCS columns if present
        if table is not None:
            if "CRVAL1" in table.columns:
                s = pd.to_numeric(table["CRVAL1"], errors="coerce")
                if s.isna().all() and table["CRVAL1"].notna().any():
                    raise ValueError("Cannot parse table['CRVAL1'] as numeric for VRAD WCS update.")
                table["CRVAL1"] = a * s + b
            if "CDELT1" in table.columns:
                s = pd.to_numeric(table["CDELT1"], errors="coerce")
                if s.isna().all() and table["CDELT1"].notna().any():
                    raise ValueError("Cannot parse table['CDELT1'] as numeric for VRAD WCS update.")
                table["CDELT1"] = a * s

        info.update({"wcs_updated": True, "a": a, "b": b, "vel_unit": unit})

    # Stamp rest2 into meta
    meta["RESTFREQ"] = rest2
    meta["RESTFRQ"] = rest2
    meta["USER_RESTFREQ"] = rest2
    # Audit/intent (non-invasive)
    meta.setdefault("RESTCONV_VELDEF", "VRAD")
    meta.setdefault("RESTCONV_SPECSYS", str(meta.get("SPECSYS", "LSRK")))

    # Overwrite restfreq columns if present in table (do not add new columns)
    if table is not None:
        if "RESTFREQ" in table.columns:
            table["RESTFREQ"] = rest2
        if "RESTFRQ" in table.columns:
            table["RESTFRQ"] = rest2

    return info
