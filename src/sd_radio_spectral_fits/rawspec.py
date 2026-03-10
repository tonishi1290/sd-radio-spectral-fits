# tools/sdrsf_pipeline/sdrsf_pipeline/rawspec.py
from __future__ import annotations

from dataclasses import dataclass
import builtins
from typing import Any, Optional

import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# Endianness / timestamp helpers
# -----------------------------------------------------------------------------
def _to_native_endian_array(a: np.ndarray) -> np.ndarray:
    """Return a native-endian view/copy of a NumPy array when needed."""
    a = np.asarray(a)
    dt = a.dtype
    if not isinstance(dt, np.dtype):
        return a
    if dt.byteorder == "|" or dt.isnative:
        return a
    return a.byteswap().view(dt.newbyteorder("="))


def _df_to_native_endian(df: pd.DataFrame) -> pd.DataFrame:
    """Convert numeric big-endian columns to native-endian for safe pandas ops."""
    if df is None or len(df.columns) == 0:
        return df
    need = False
    for col in df.columns:
        dt = getattr(df[col], "dtype", None)
        if isinstance(dt, np.dtype) and dt.kind in ("i", "u", "f", "c") and (not dt.isnative):
            need = True
            break
    if not need:
        return df
    out = df.copy(deep=False)
    for col in out.columns:
        s = out[col]
        dt = getattr(s, "dtype", None)
        if isinstance(dt, np.dtype) and dt.kind in ("i", "u", "f", "c") and (not dt.isnative):
            out[col] = _to_native_endian_array(s.to_numpy(copy=False))
    return out


def _resolve_mapping_timestamps(mapping: pd.DataFrame) -> pd.DatetimeIndex | None:
    """Resolve UTC timestamps from standard/legacy time columns.

    Priority:
      1. TIMESTAMP
      2. MJD
      3. DATE-OBS / DATEOBS (+ TIME seconds offset if present)
      4. legacy TIME (MJD days or Unix seconds)
    """
    if mapping is None or len(mapping) == 0:
        return pd.DatetimeIndex([])

    if "TIMESTAMP" in mapping.columns:
        ts = pd.to_datetime(mapping["TIMESTAMP"], utc=True, errors="coerce")
        if ts.notna().any():
            return pd.DatetimeIndex(ts)

    if "MJD" in mapping.columns:
        arr = pd.to_numeric(mapping["MJD"], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(arr).any():
            ts = pd.to_datetime(arr, unit="D", origin=pd.Timestamp("1858-11-17"), utc=True, errors="coerce")
            return pd.DatetimeIndex(ts)

    for dcol in ("DATE-OBS", "DATEOBS"):
        if dcol in mapping.columns:
            base = pd.to_datetime(mapping[dcol].astype(str), utc=True, errors="coerce")
            if base.notna().any():
                ts = pd.DatetimeIndex(base)
                if "TIME" in mapping.columns:
                    offs = pd.to_numeric(mapping["TIME"], errors="coerce").to_numpy(dtype=float)
                    if np.isfinite(offs).any():
                        ts = pd.DatetimeIndex(ts + pd.to_timedelta(np.nan_to_num(offs, nan=0.0), unit="s"))
                return ts

    if "TIME" in mapping.columns:
        arr = pd.to_numeric(mapping["TIME"], errors="coerce").to_numpy(dtype=float)
        finite = arr[np.isfinite(arr)]
        if finite.size:
            if float(np.nanmedian(finite)) < 1.0e6:
                ts = pd.to_datetime(arr, unit="D", origin=pd.Timestamp("1858-11-17"), utc=True, errors="coerce")
            else:
                ts = pd.to_datetime(arr, unit="s", utc=True, errors="coerce")
            return pd.DatetimeIndex(ts)

    return None


@dataclass
class RawSpec:
    """
    Raw single-dish spectra container.

    Supports:
      - Pickle-based internal format (load_rawspec/save_rawspec)
      - FITS SDFITS-style (BinTable rows; OBSMODE=HOT/OFF/ON; spectrum in vector columns)
      - Fallback ImageHDU-style (HOT/ON/OFF ImageHDUs)

    Notes:
      - For SDFITS BinTable style, we keep scalar metadata columns in `mapping`
        (vector columns like DATA/SPECTRUM/FLAG are dropped).
      - We normalize essential axis keys to pipeline internal names:
          crval1_hz, cdelt1_hz, crpix1, rest_hz/restfrq_hz, nchan
      - Site keys in SDFITS table header are captured into meta:
          SITELAT, SITELONG, SITEELEV, OBSGEO-X/Y/Z
    """
    meta: dict[str, Any]
    hot: np.ndarray
    on: np.ndarray
    off: np.ndarray
    mapping: pd.DataFrame

    # dict-like compatibility
    def __getitem__(self, key: str):
        k = str(key)
        if k == "meta":
            return self.meta
        if k.lower() == "hot":
            return self.hot
        if k.lower() == "on":
            return self.on
        if k.lower() == "off":
            return self.off
        if k.lower() in ("mapping", "table"):
            return self.mapping
        raise KeyError(k)

    def get(self, key: str, default=None):
        try:
            return self[key]
        except KeyError:
            return default


def build_rawspec(
    *,
    hot: pd.DataFrame,
    on: pd.DataFrame,
    off: pd.DataFrame,
    meta: dict[str, Any] | None,
    mapping: pd.DataFrame | None,
) -> RawSpec:
    """Build a RawSpec object from HOT/ON/OFF spectra DataFrames.

    Notes
    -----
    - `hot`, `on`, `off` must be pandas DataFrame with a DatetimeIndex (UTC recommended).
    - `mapping` is optional. If provided and does NOT contain an `OBSMODE` column,
      it is interpreted as an ON-only mapping table (e.g., position IDs, scan IDs, etc.).
      HOT/OFF rows will be created automatically.
    - Internally we always store a single `mapping_all` table that includes all rows
      (HOT/OFF/ON) with mandatory columns: `TIMESTAMP` (datetime64[ns, UTC]) and `OBSMODE`.
    """

    if meta is None:
        meta = {}

    # ---- basic validation (nchan)
    def _nchan(df: pd.DataFrame) -> int:
        if df.ndim != 2:
            raise ValueError('spectra must be a 2-D DataFrame')
        return int(df.shape[1])

    n_hot = _nchan(hot)
    n_on  = _nchan(on)
    n_off = _nchan(off)

    n_meta = meta.get('nchan', None)
    if n_meta is not None:
        n_meta = int(n_meta)
        if (n_hot != n_meta) or (n_on != n_meta) or (n_off != n_meta):
            raise ValueError(
                f"nchan mismatch: meta nchan={n_meta}, hot={n_hot}, on={n_on}, off={n_off}"
            )
    else:
        # infer
        meta['nchan'] = int(n_on)
        n_meta = int(n_on)

    # ---- normalize / build mapping_all
    def _ensure_dtindex(df: pd.DataFrame, name: str) -> None:
        if not isinstance(df.index, pd.DatetimeIndex):
            raise ValueError(f"{name} index must be a pandas.DatetimeIndex")

    _ensure_dtindex(hot, 'hot')
    _ensure_dtindex(off, 'off')
    _ensure_dtindex(on,  'on')

    # make sure indices are tz-aware (assume UTC if naive)
    def _tz_utc(idx: pd.DatetimeIndex) -> pd.DatetimeIndex:
        if idx.tz is None:
            return idx.tz_localize('UTC')
        return idx.tz_convert('UTC')

    hot_idx = _tz_utc(hot.index)
    off_idx = _tz_utc(off.index)
    on_idx  = _tz_utc(on.index)

    # ON mapping: user-supplied mapping is interpreted as ON-only unless it already
    # contains OBSMODE.
    if mapping is None:
        on_map = pd.DataFrame(index=on_idx)
    else:
        on_map = mapping.copy()
        # if mapping has a timestamp column but not datetime index, normalize
        if not isinstance(on_map.index, pd.DatetimeIndex):
            if 'timestamp' in on_map.columns:
                on_map = on_map.set_index('timestamp')
            elif 'TIMESTAMP' in on_map.columns:
                on_map = on_map.set_index('TIMESTAMP')
            else:
                raise ValueError('mapping must have a DatetimeIndex or a timestamp/TIMESTAMP column')
        # align length
        if len(on_map) != len(on_idx):
            raise ValueError(f"mapping length ({len(on_map)}) does not match ON dumps ({len(on_idx)})")
        # force index to ON times (assume row-order match)
        on_map = on_map.copy()
        on_map.index = on_idx

    def _mk_map(idx: pd.DatetimeIndex, mode: str) -> pd.DataFrame:
        m = pd.DataFrame(index=idx)
        m['TIMESTAMP'] = idx
        m['OBSMODE'] = mode
        return m

    hot_map = _mk_map(hot_idx, 'HOT')
    off_map = _mk_map(off_idx, 'OFF')

    # ON map: keep user columns, but ensure mandatory columns exist
    on_map = on_map.copy()
    on_map['TIMESTAMP'] = on_idx
    on_map['OBSMODE'] = 'ON'

    mapping_all = pd.concat([hot_map, off_map, on_map], axis=0, ignore_index=True, sort=False)

    # ---- finalize RawSpec
    return RawSpec(
        meta=dict(meta),
        hot=hot.to_numpy(float),
        on=on.to_numpy(float),
        off=off.to_numpy(float),
        mapping=_df_to_native_endian(mapping_all),
    )

def save_rawspec(raw: RawSpec, path: str) -> None:
    import pickle

    with open(path, "wb") as f:
        pickle.dump(
            {
                "meta": raw.meta,
                "hot": raw.hot,
                "on": raw.on,
                "off": raw.off,
                "mapping": raw.mapping,
            },
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )


def load_rawspec(path: str) -> RawSpec:
    import pickle

    with open(path, "rb") as f:
        d = pickle.load(f)

    mapping = d.get("mapping", pd.DataFrame())
    if isinstance(mapping, dict):
        mapping = pd.DataFrame(mapping)
    mapping = _df_to_native_endian(mapping)
    ts = _resolve_mapping_timestamps(mapping)
    if ts is not None and len(ts) == len(mapping):
        mapping = mapping.copy()
        mapping["TIMESTAMP"] = ts
    return RawSpec(
        meta=dict(d.get("meta", {})),
        hot=np.asarray(d["hot"], float),
        on=np.asarray(d["on"], float),
        off=np.asarray(d["off"], float),
        mapping=mapping,
    )


def load_rawspec_auto(path: str, prefer: tuple[str, ...] | None = None) -> RawSpec:
    """Auto-load raw spectra.

    Notes
    -----
    - This function is used by CLI scripts.
    - `prefer` is accepted for forward/backward compatibility.
      (Current loader supports FITS SDFITS-like input and pickled RawSpec.)
    """
    # FITS判定: 拡張子 or 先頭マジック
    head = b""
    try:
        with open(path, "rb") as f:
            head = f.read(16)
    except Exception:
        pass
    is_fits = path.lower().endswith((".fits", ".fit", ".fts")) or head.startswith(b"SIMPLE")
    if is_fits:
        return load_rawspec_fits(path)
    return load_rawspec(path)


def load_rawspec_fits(path: str) -> RawSpec:
    """
    Load raw spectra from FITS.

    Primary supported layout (SDFITS-like):
      - A BinTableHDU (often EXTNAME='SINGLE DISH') with columns:
          OBSMODE (HOT/OFF/ON), and spectrum vector column DATA (preferred) or SPECTRUM.
      - WCS / site keywords are typically stored in the BinTable header:
          CTYPE1/CRVAL1/CDELT1/CRPIX1, RESTFRQ, SITELAT/SITELONG/SITEELEV, OBSGEO-X/Y/Z, ...

    Fallback:
      - ImageHDUs named HOT/ON/OFF (or LOAD) with 2D arrays (Ndump,Nchan).
    """
    from astropy.io import fits

    hdul = fits.open(path)

    # --- SDFITS BinTable detection ---
    table_hdu = None
    for hh in hdul:
        if type(hh).__name__ != "BinTableHDU":
            continue
        cols = set(hh.columns.names or [])
        if "OBSMODE" not in cols:
            continue
        if ("DATA" in cols) or ("SPECTRUM" in cols):
            table_hdu = hh
            break

    if table_hdu is not None:
        h = table_hdu
        tab = h.data
        cols = list(h.columns.names)

        # Build meta from PRIMARY + BinTable header (BinTable overwrites)
        meta: dict[str, Any] = dict(hdul[0].header)
        meta.update(dict(h.header))

        # Ensure site keywords are captured from BinTable header explicitly
        for k in ("SITELAT", "SITELONG", "SITEELEV", "OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"):
            if k in h.header and h.header[k] not in (None, ""):
                meta[k] = h.header[k]

        # Choose spectrum column
        spec_col = "DATA" if "DATA" in cols else "SPECTRUM"

        obsmode = np.char.upper(np.char.strip(np.asarray(tab["OBSMODE"]).astype(str)))

        def pick(mode: str) -> np.ndarray:
            m = (obsmode == mode)
            if not np.any(m):
                raise KeyError(f"Missing rows with OBSMODE={mode} in BinTable '{h.name}'.")
            arr = np.asarray(tab[spec_col][m], float)
            if arr.ndim != 2:
                raise ValueError(f"Spectrum column '{spec_col}' must be 2-D (Ndump,Nchan), got shape={arr.shape}")
            return arr

        hot = pick("HOT")
        on = pick("ON")
        off = pick("OFF")

        # mapping: keep only scalar columns (drop vector columns like DATA/SPECTRUM/FLAG)
        drop_vec = {"DATA", "SPECTRUM", "FLAG"}
        mapping_dict: dict[str, Any] = {}
        for c in cols:
            if c in drop_vec:
                continue
            v = tab[c]
            if getattr(v, "ndim", 1) != 1:
                continue
            mapping_dict[c] = np.asarray(v)
        mapping = _df_to_native_endian(pd.DataFrame(mapping_dict))

        # Mark provenance
        mapping["_SPEC_COL"] = str(spec_col)
        mapping["_TABLE_HDU"] = str(h.name)

        # Resolve robust UTC timestamps from standard / legacy time columns.
        ts = _resolve_mapping_timestamps(mapping)
        if ts is not None and len(ts) == len(mapping):
            mapping["TIMESTAMP"] = ts

        meta["CTYPE1"] = str(meta["CTYPE1"]).strip()

        # REST frequency: keep multiple aliases for compatibility
        rest = None
        if "RESTFRQ" in meta and meta["RESTFRQ"] not in (None, ""):
            rest = float(meta["RESTFRQ"])
        elif "RESTFREQ" in meta and meta["RESTFREQ"] not in (None, ""):
            rest = float(meta["RESTFREQ"])
        if rest is not None:
            meta["RESTFRQ"] = rest
            meta["RESTFREQ"] = rest

        meta["NAXIS1"] = int(on.shape[1])

        # Site internal aliases (optional but useful)
        if "SITELAT" in meta and "SITELONG" in meta:
            meta["site_lat_deg"] = float(meta["SITELAT"])
            meta["site_lon_deg"] = float(meta["SITELONG"])
            if "SITEELEV" in meta and meta["SITEELEV"] not in (None, ""):
                meta["site_height_m"] = float(meta["SITEELEV"])
        if "OBSGEO-X" in meta and "OBSGEO-Y" in meta and "OBSGEO-Z" in meta:
            meta["obsgeo_x_m"] = float(meta["OBSGEO-X"])
            meta["obsgeo_y_m"] = float(meta["OBSGEO-Y"])
            meta["obsgeo_z_m"] = float(meta["OBSGEO-Z"])

        return RawSpec(meta=meta, hot=hot, on=on, off=off, mapping=mapping)

    # ----------------------------
    # Fallback: ImageHDU-style
    # ----------------------------
    meta = dict(hdul[0].header)

    def get_img(names):
        for nm in names:
            if nm in hdul and getattr(hdul[nm], "data", None) is not None:
                return np.asarray(hdul[nm].data, float)
        return None

    hot = get_img(["HOT", "hot", "LOAD", "load", "AMBIENT", "ambient", "TAMB", "Tamb"])
    on = get_img(["ON", "on"])
    off = get_img(["OFF", "off"])

    if hot is None or on is None or off is None:
        raise KeyError(
            "Raw FITS does not contain a BinTable with OBSMODE+DATA/SPECTRUM, "
            "and also missing ImageHDUs for HOT/ON/OFF."
        )

    mapping = pd.DataFrame()
    if "MAPPING" in hdul and hdul["MAPPING"].data is not None:
        tab = hdul["MAPPING"].data
        mapping = pd.DataFrame({n: np.asarray(tab[n]) for n in tab.names})

    # Normalize axis keys if present in PRIMARY
    rest = None
    if "RESTFRQ" in meta and meta["RESTFRQ"] not in (None, ""):
        rest = float(meta["RESTFRQ"])
    elif "RESTFREQ" in meta and meta["RESTFREQ"] not in (None, ""):
        rest = float(meta["RESTFREQ"])
    if rest is not None:
        meta["RESTFRQ"] = rest
        meta["RESTFREQ"] = rest

    meta["NAXIS1"] = int(on.shape[1])

    return RawSpec(meta=meta, hot=hot, on=on, off=off, mapping=mapping)
