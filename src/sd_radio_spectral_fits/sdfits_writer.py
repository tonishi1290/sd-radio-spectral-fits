# src/sd_radio_spectral_fits/sdfits_writer.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
互換性最優先 SDFITS Writer（Single-Dish Spectral FITS）
==========================================================

更新点 (v1.5.0)
------------------------
- 可変長配列 (Variable Length Array: VLA) の書き込みに対応。
  チャンネル数が行ごとに異なる場合、自動的に TFORM='PE'/'PD'/'PL' を使用します。
- 周波数/WCS関連パラメータをデータ列(Column)として出力する機能を維持。
"""

from __future__ import annotations

import json
import math
import builtins
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from .tempscale import normalize_tempscal

# NOTE: Single Source of Truth for FITS BinTable writing
try:
    from .sdfits_bintable import build_single_dish_table_hdu, build_history_hdu, set_meta_keyword
except Exception:  # pragma: no cover
    from sdfits_bintable import build_single_dish_table_hdu, build_history_hdu, set_meta_keyword

from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import EarthLocation, SkyCoord, FK5, FK4, Galactic

__version__ = "1.5.1"
SWNAME = "SDFITS_WRITER"


Number = Union[int, float, np.number]

# -----------------------------------------------------------------------------
# Column policy (final)
# -----------------------------------------------------------------------------
# ALWAYS: row の意味やスペクトル軸解釈に必須な列
# CONTEXT-REQUIRED: 文脈に応じて必須（例: 速度解釈時の RESTFREQ / VELDEF）
# OPTIONAL: 値があるときだけ列として出す
# BLOCK-OPTIONAL: 関連列群のどれか1つでも意味があれば、ブロック全体を出す
ALWAYS_COLUMNS = {
    "TIME", "MJD", "DATE-OBS", "DATEOBS", "TIMESTAMP", "SCAN", "SUBSCAN", "INTGRP", "OBJECT", "OBSMODE",
    "EXPOSURE", "CALSTAT", "TEMPSCAL", "FLAGROW",
    "RA", "DEC", "GLON", "GLAT",
    "CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1", "SPECSYS",
    "FDNUM", "IFNUM", "PLNUM", "POLARIZA",
    "DATA", "FLAG",
}

CONTEXT_REQUIRED_COLUMNS = {
    "RESTFREQ", "VELDEF",
    # LO / sideband interpretation becomes mandatory only when the row carries
    # heterodyne sideband semantics.
    "OBSFREQ", "SIDEBAND",
}

OPTIONAL_COLUMNS = {
    "TCAL", "THOT", "TSYS", "TAU0",
    "AZIMUTH", "ELEVATIO",
    "VFRAME", "FOFFSET",
    "BORE_AZ", "BORE_EL", "BEAMXOFF", "BEAMYOFF", "BEAMROT",
    "FRONTEND", "BACKEND", "SAMPLER",
    "IMAGFREQ", "LO1FREQ", "LO2FREQ", "LO3FREQ",
    "SB1", "SB2", "SB3",
    "FREQ",
}

BLOCK_OPTIONAL = {
    "SOURCE": {"SRCFRAME", "SRCRDSYS", "SRCEQNX", "SRC_LONG", "SRC_LAT"},
    "SCAN": {"SCANFRAM", "SCANRDSYS", "SCANEQNX", "SCANX", "SCANY"},
    "POINTING": {"AZ_CMD", "EL_CMD", "CORR_AZ", "CORR_EL", "CALC_REFR"},
    "WEATHER": {"TAMBIENT", "PRESSURE", "HUMIDITY", "WINDSPD", "WINDDIR"},
}


def _is_meaningful_str(x: Any) -> bool:
    """True if x is a non-empty string-like value with semantic content."""
    return x is not None and str(x).strip() != ""


def _is_meaningful_scalar(x: Any) -> bool:
    """True if x is not None / NaN / empty-string.

    This is the core helper used to decide whether an optional column should be
    emitted into the FITS table.
    """
    if x is None:
        return False
    if isinstance(x, str):
        return x.strip() != ""
    try:
        return not bool(pd.isna(x))
    except Exception:
        return True


def _series_has_meaningful(xs: Sequence[Any]) -> bool:
    """True if any element in the sequence carries meaningful information."""
    return any(_is_meaningful_scalar(x) for x in xs)


def _block_has_meaningful(buffer: Any, names: Sequence[str]) -> bool:
    """True if any column inside a block-optional group has data worth writing."""
    return any(_series_has_meaningful(getattr(buffer, name)) for name in names)

# =============================================================================
# 1) 列挙値
# =============================================================================
RADESYS_VALUES = ("ICRS", "FK5", "FK4", "FK4-NO-E", "GALACTIC")
SPECSYS_VALUES = ("TOPOCENT", "LSRK", "BARYCENT", "GEOCENTR", "HELIOCEN")
VELDEF_VALUES  = ("RADIO", "OPTICAL", "RELATIVISTIC", "RADI", "OPTI", "RELA")
VELDEF_FRAME_CODES = {
    "TOPOCENT": "-OBS",
    "LSRK": "-LSR",
    "HELIOCEN": "-HEL",
    "BARYCENT": "-BAR",
    "GEOCENTR": "-GEO",
}
COORDTYPE_VALUES = ("RADEC", "GALACTIC", "AZEL", "AZEL_GEO")

OBSMODE_VALUES = ("ON", "OFF", "HOT", "SKY", "FS_SIG", "FS_REF", "CAL_ON", "CAL_OFF", "UNKNOWN")
CALSTAT_VALUES = ("RAW", "TASTAR", "TA", "TMB", "INTEGRATED")
EFFSTAT_VALUES = ("ASSUMED", "MEASURED", "UNKNOWN")
POLARIZA_VALUES = ("XX", "YY", "RR", "LL", "XY", "YX", "RL", "LR", "I", "Q", "U", "V")
SIDEBAND_VALUES = ("USB", "LSB")

# =============================================================================
# 2) 仕様一覧（SCHEMA）
# =============================================================================
SCHEMA: Dict[str, Dict[str, Any]] = {
    "writer.n_chan": {"type": "int", "required": True, "meaning": "チャンネル数 (デフォルト/最大値)"},
    "writer.store_freq_column": {"type": "bool", "required": False, "meaning": "Trueなら FREQ(Hzベクトル) を各行に保持"},
    "writer.chunk_size": {"type": "Optional[int]", "required": False, "meaning": "巨大OTF用 parts 分割行数"},
    "writer.out_basename": {"type": "Optional[str]", "required": "chunk_size is not None", "meaning": "parts基底名"},
    "writer.spectrum_column": {"type": "str", "required": True, "meaning": "スペクトル配列カラム名"},
}

# =============================================================================
# 3) dataclasses
# =============================================================================

@dataclass(frozen=True)
class Site:
    """観測地点（EarthLocationへ変換）。"""
    lat_deg: float
    lon_deg: float
    elev_m: float

    def to_earthlocation(self) -> EarthLocation:
        return EarthLocation.from_geodetic(lon=self.lon_deg, lat=self.lat_deg, height=self.elev_m)

    def obsgeo_xyz_m(self) -> Tuple[float, float, float]:
        loc = self.to_earthlocation()
        x, y, z = loc.to_geocentric()
        return float(x.to_value("m")), float(y.to_value("m")), float(z.to_value("m"))


@dataclass
class Efficiency:
    beameff: float = np.nan
    apereff: float = np.nan
    mooneff: float = np.nan
    suneff: float = np.nan
    effstat: str = "UNKNOWN"


@dataclass
class SpectralAxisUniform:
    """
    等間隔周波数WCS（デフォルト値として使用）。
    FREQ[i] = CRVAL1 + ((i+1)-CRPIX1)*CDELT1
    """
    crval1_hz: float
    cdelt1_hz: float
    crpix1: float = 1.0
    restfreq_hz: float = 0.0
    specsys: str = "TOPOCENT"
    ssysobs: str = "TOPOCENT"
    veldef: str = "RADIO"
    ctype1: str = "FREQ"
    cunit1: str = "Hz"
    refchan: int = 1

    def validate(self) -> None:
        if self.specsys not in SPECSYS_VALUES:
            raise ValueError(f"SPECSYS must be one of {SPECSYS_VALUES}, got {self.specsys}")
        if self.ssysobs not in SPECSYS_VALUES:
            raise ValueError(f"SSYSOBS must be one of {SPECSYS_VALUES}, got {self.ssysobs}")
        if self.veldef not in VELDEF_VALUES:
            raise ValueError(f"VELDEF must be one of {VELDEF_VALUES}, got {self.veldef}")
        if self.ctype1 != "FREQ":
            raise ValueError("CTYPE1 should be 'FREQ'")
        if self.cunit1 != "Hz":
            raise ValueError("CUNIT1 should be 'Hz'")
        if self.refchan < 1:
            raise ValueError("REFCHAN must be >=1")
        if self.restfreq_hz < 0:
            raise ValueError("RESTFREQ must be >=0")

    def velref(self) -> int:
        """Return AIPS/casacore-style VELREF integer. Add +256 for RADIO."""
        base_map = {
            "LSRK": 1,
            "BARYCENT": 2,
            "HELIOCEN": 2,
            "TOPOCENT": 3,
            "GEOCENTR": 5,
        }
        base = base_map.get(self.specsys, 0)
        if base == 0:
            return 0
        if self.veldef == "RADIO":
            base += 256
        return int(base)


@dataclass
class DatasetInfo:
    telescope: str
    observer: str = "Unknown"
    project: str = "ProjectID"
    object_name: str = "Target"

    radesys: str = "ICRS"
    equinox: float = 2000.0

    src_radesys: str = "ICRS"
    src_equinox: float = 2000.0

    bmaj_deg: float = np.nan
    bmin_deg: float = np.nan
    bpa_deg: float = np.nan

    eff: Efficiency = field(default_factory=Efficiency)

    refr_included_in_corr: bool = True
    doppler_tracking_applied: bool = False

    spectral_axis: Optional[SpectralAxisUniform] = None
    shared_meta: Dict[str, Any] = field(default_factory=dict)

    def validate(self, store_freq_column: bool) -> None:
        if self.radesys not in RADESYS_VALUES:
            raise ValueError(f"RADESYS must be one of {RADESYS_VALUES}, got {self.radesys}")
        if self.src_radesys not in RADESYS_VALUES:
            raise ValueError(f"SRC_RADESYS must be one of {RADESYS_VALUES}, got {self.src_radesys}")
        if self.eff.effstat not in EFFSTAT_VALUES:
            raise ValueError(f"EFFSTAT must be one of {EFFSTAT_VALUES}, got {self.eff.effstat}")
        if (not store_freq_column) and (self.spectral_axis is None):
            raise ValueError("store_freq_column=False requires info.spectral_axis (uniform WCS).")
        if self.spectral_axis is not None:
            self.spectral_axis.validate()
            if (self.spectral_axis.specsys == "LSRK") and (not store_freq_column):
                raise ValueError("SPECSYS='LSRK' requires store_freq_column=True (time-dependent axis).")


# =============================================================================
# 4) ColumnBuffer（列ごとに蓄積）
# =============================================================================

@dataclass
class _ColumnBuffer:
    n_chan: int
    store_freq_column: bool

    TIME: List[float] = field(default_factory=list)      # Seconds since DATE-OBS (standard SDFITS meaning)
    MJD: List[float] = field(default_factory=list)       # Absolute UTC time in MJD days (custom convenience column)
    DATEOBS: List[str] = field(default_factory=list)     # ISO UTC (DATE-OBS string)

    SCAN: List[int] = field(default_factory=list)
    SUBSCAN: List[int] = field(default_factory=list)
    INTGRP: List[int] = field(default_factory=list)

    OBJECT: List[str] = field(default_factory=list)
    OBSMODE: List[str] = field(default_factory=list)

    # beam / pol / IF の正規列（BEAMNO/POLNO は廃止）
    FDNUM: List[int] = field(default_factory=list)      # Feed / beam number
    IFNUM: List[int] = field(default_factory=list)      # IF / spectral-window number
    PLNUM: List[int] = field(default_factory=list)      # Internal polarization index
    POLARIZA: List[str] = field(default_factory=list)   # Physical polarization label (XX/YY/RR/LL/...)

    # 分光計 / 受信機 provenance（値が分かるときだけ列を出す）
    FRONTEND: List[Optional[str]] = field(default_factory=list)
    BACKEND: List[Optional[str]] = field(default_factory=list)
    SAMPLER: List[Optional[str]] = field(default_factory=list)

    # Heterodyne LO / sideband metadata.
    # OBSFREQ/SIDEBAND are context-required when sideband semantics are present.
    OBSFREQ: List[float] = field(default_factory=list)
    IMAGFREQ: List[float] = field(default_factory=list)
    LO1FREQ: List[float] = field(default_factory=list)
    LO2FREQ: List[float] = field(default_factory=list)
    LO3FREQ: List[float] = field(default_factory=list)
    SIDEBAND: List[Optional[str]] = field(default_factory=list)
    SB1: List[Optional[str]] = field(default_factory=list)
    SB2: List[Optional[str]] = field(default_factory=list)
    SB3: List[Optional[str]] = field(default_factory=list)

    EXPOSURE: List[float] = field(default_factory=list)
    CALSTAT: List[str] = field(default_factory=list)
    TEMPSCAL: List[str] = field(default_factory=list)
    TCAL: List[float] = field(default_factory=list)
    THOT: List[float] = field(default_factory=list)
    TSYS: List[float] = field(default_factory=list)
    TAU0: List[float] = field(default_factory=list)  # ★ TAU -> TAU0 に変更

    FLAGROW: List[int] = field(default_factory=list)

    RA: List[float] = field(default_factory=list)
    DEC: List[float] = field(default_factory=list)
    GLON: List[float] = field(default_factory=list)
    GLAT: List[float] = field(default_factory=list)

    AZIMUTH: List[float] = field(default_factory=list)
    ELEVATIO: List[float] = field(default_factory=list)
    BORE_AZ: List[float] = field(default_factory=list)
    BORE_EL: List[float] = field(default_factory=list)
    BEAMXOFF: List[float] = field(default_factory=list)
    BEAMYOFF: List[float] = field(default_factory=list)
    BEAMROT: List[float] = field(default_factory=list)

    # Source-coordinate block (block-optional)
    SRCFRAME: List[Optional[str]] = field(default_factory=list)
    SRCRDSYS: List[Optional[str]] = field(default_factory=list)
    SRCEQNX: List[float] = field(default_factory=list)
    SRC_LONG: List[float] = field(default_factory=list)
    SRC_LAT: List[float] = field(default_factory=list)

    # Scan-coordinate block (block-optional)
    SCANFRAM: List[Optional[str]] = field(default_factory=list)
    SCANRDSYS: List[Optional[str]] = field(default_factory=list)
    SCANEQNX: List[float] = field(default_factory=list)
    SCANX: List[float] = field(default_factory=list)
    SCANY: List[float] = field(default_factory=list)

    AZ_CMD: List[float] = field(default_factory=list)
    EL_CMD: List[float] = field(default_factory=list)
    CORR_AZ: List[float] = field(default_factory=list)
    CORR_EL: List[float] = field(default_factory=list)
    CALC_REFR: List[float] = field(default_factory=list)

    VFRAME: List[float] = field(default_factory=list)
    FOFFSET: List[float] = field(default_factory=list)

    TAMBIENT: List[float] = field(default_factory=list)
    PRESSURE: List[float] = field(default_factory=list)
    HUMIDITY: List[float] = field(default_factory=list)
    WINDSPD: List[float] = field(default_factory=list)
    WINDDIR: List[float] = field(default_factory=list)

    # --- WCS / Frequency Parameters per row ---
    RESTFREQ: List[float] = field(default_factory=list)
    CRVAL1: List[float] = field(default_factory=list)
    CDELT1: List[float] = field(default_factory=list)
    CRPIX1: List[float] = field(default_factory=list)
    CTYPE1: List[str] = field(default_factory=list)
    CUNIT1: List[str] = field(default_factory=list)
    SPECSYS: List[str] = field(default_factory=list)
    VELDEF: List[Optional[str]] = field(default_factory=list)
    # ------------------------------------------------

    _DATA: List[np.ndarray] = field(default_factory=list)
    _FLAG: List[np.ndarray] = field(default_factory=list)
    _FREQ: List[np.ndarray] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.TIME)

    def clear(self) -> None:
        for v in self.__dict__.values():
            if isinstance(v, list):
                v.clear()

    def append_spectrum(self, spec: Union[np.ndarray, float], flag: Optional[np.ndarray]) -> None:
        # [MODIFIED] Relax n_chan check for VLA support
        x = np.asarray(spec, dtype=np.float32).reshape(-1)
        self._DATA.append(x)
        
        current_len = x.size
        
        if flag is None:
            self._FLAG.append(np.zeros(current_len, dtype=bool))
        else:
            f = np.asarray(flag, dtype=bool).reshape(-1)
            if f.size != current_len:
                raise ValueError(f"FLAG length {f.size} != DATA length {current_len}")
            self._FLAG.append(f)

    def append_freq(self, freq_hz: np.ndarray) -> None:
        if not self.store_freq_column:
            return
        f = np.asarray(freq_hz, dtype=np.float64).reshape(-1)
        # [MODIFIED] Relax check, match with DATA length
        expected = self._DATA[-1].size if self._DATA else self.n_chan
        if f.size != expected:
            raise ValueError(f"FREQ length {f.size} != DATA length {expected}")
        self._FREQ.append(f)


# =============================================================================
# 5) Writer本体
# =============================================================================

class SDRadioSpectralSDFITSWriter:
    """
    CASA互換を優先した SDFITS-like FITS Writer (VLA対応)

    - スペクトル配列カラムは 'DATA' のみ（SPECTRUM は出さない）
    - チャンネル数が変動する場合、TFORM='PE'/'PL'/'PD' (Variable Length Array) を使用
    """

    STR_OBJECT  = 32
    STR_OBSMODE = 24
    STR_CALSTAT = 16
    STR_FRAME   = 16
    STR_CTYPE   = 8
    STR_CUNIT   = 8
    STR_SPECSYS = 16
    STR_VELDEF  = 8

    SPECTRUM_COLUMN = "DATA"

    def __init__(
        self,
        *,
        n_chan: int,
        site: Site,
        info: DatasetInfo,
        store_freq_column: bool = False,
        chunk_size: Optional[int] = None,
        out_basename: Optional[str] = None,
        history: Any | None = None,
    ):
        self.n_chan = int(n_chan) # デフォルト/最大の目安として保持
        if self.n_chan < 1:
            raise ValueError("n_chan must be >= 1")

        self.site = site
        self.store_freq_column = bool(store_freq_column)

        info.validate(self.store_freq_column)
        self.info = info

        self.chunk_size = int(chunk_size) if chunk_size is not None else None
        if self.chunk_size is not None and self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive int or None")
        self.out_basename = out_basename
        # Optional file-level HISTORY (dict/list/str).
        if isinstance(history, dict):
            self._history: Any | None = dict(history)
        elif isinstance(history, (list, tuple)):
            self._history = list(history)
        else:
            self._history = history

        if self.chunk_size is not None and not self.out_basename:
            raise ValueError("parts mode requires out_basename, e.g. 'otf_run1'.")

        self._buf = _ColumnBuffer(n_chan=self.n_chan, store_freq_column=self.store_freq_column)
        self._part_index = 0
        self._manifest: Dict[str, Any] = {
            "format": "sd-radio-spectral-fits parts manifest",
            "n_chan": self.n_chan,
            "store_freq_column": self.store_freq_column,
            "spectrum_column": self.SPECTRUM_COLUMN,
            "parts": [],
        }

    # ... (schema helper, coordinate helpers, freq helpers: 省略なしで保持) ...
    # [この部分は既存コードと同じため、変更点のみ記述します]
    @staticmethod
    def print_schema() -> None:
        keys = sorted(SCHEMA.keys())
        for k in keys:
            d = SCHEMA[k]
            print(f"\n[{k}] {d.get('meaning','')}")
            for kk in ("type", "required"):
                if kk in d:
                    print(f"  {kk}: {d[kk]}")

    def _radec_frame_for_columns(self):
        rs = self.info.radesys
        if rs == "ICRS":
            return "icrs"
        if rs == "FK5":
            return FK5(equinox=Time(self.info.equinox, format="jyear"))
        if rs in ("FK4", "FK4-NO-E"):
            return FK4(equinox=Time(self.info.equinox, format="jyear"))
        if rs == "GALACTIC":
            return Galactic()
        return "icrs"

    def _fill_radec_glon_glat(
        self,
        ra_deg: float,
        dec_deg: float,
        glon_deg: float,
        glat_deg: float,
    ) -> Tuple[float, float, float, float]:
        have_radec = np.isfinite(ra_deg) and np.isfinite(dec_deg)
        have_gal = np.isfinite(glon_deg) and np.isfinite(glat_deg)

        if (not have_radec) and (not have_gal):
            raise ValueError("Provide either (ra_deg, dec_deg) or (glon_deg, glat_deg).")

        if have_radec and (not have_gal):
            sc = SkyCoord(ra=float(ra_deg)*u.deg, dec=float(dec_deg)*u.deg, frame=self._radec_frame_for_columns())
            gal = sc.galactic
            glon_deg = float(gal.l.to_value(u.deg))
            glat_deg = float(gal.b.to_value(u.deg))
            return float(ra_deg), float(dec_deg), glon_deg, glat_deg

        if have_gal and (not have_radec):
            gal = SkyCoord(l=float(glon_deg)*u.deg, b=float(glat_deg)*u.deg, frame=Galactic())
            radec = gal.transform_to(self._radec_frame_for_columns())
            if hasattr(radec, "ra") and hasattr(radec, "dec"):
                ra_deg = float(radec.ra.to_value(u.deg))
                dec_deg = float(radec.dec.to_value(u.deg))
            else:
                ra_deg = float(gal.l.to_value(u.deg))
                dec_deg = float(gal.b.to_value(u.deg))
            return ra_deg, dec_deg, float(glon_deg), float(glat_deg)

        return float(ra_deg), float(dec_deg), float(glon_deg), float(glat_deg)

    def _freq_from_uniform_axis(self) -> np.ndarray:
        ax = self.info.spectral_axis
        if ax is None:
            raise ValueError("No uniform spectral_axis to generate freq.")
        pix = (np.arange(self.n_chan, dtype=np.float64) + 1.0)  # 1-based
        return ax.crval1_hz + (pix - float(ax.crpix1)) * ax.cdelt1_hz

    def _calc_vobs_kms_user(self, ra_deg: float, dec_deg: float, time_mjd: float) -> u.Quantity:
        loc = self.site.to_earthlocation()
        unixtime = Time(float(time_mjd), format="mjd", scale="utc").unix
        tobs = Time(unixtime, format="unix", scale="utc")

        target_frame = self._radec_frame_for_columns()
        target = SkyCoord(
            float(ra_deg) * u.deg,
            float(dec_deg) * u.deg,
            frame=target_frame,
            obstime=tobs,
            location=loc,
        )

        vobs = SkyCoord(loc.get_gcrs(tobs)).transform_to(coord.LSRK()).velocity

        vcorrection = (
            vobs.d_x * np.cos(target.icrs.dec) * np.cos(target.icrs.ra)
            + vobs.d_y * np.cos(target.icrs.dec) * np.sin(target.icrs.ra)
            + vobs.d_z * np.sin(target.icrs.dec)
        )
        return vcorrection.to(u.km/u.s)

    def _topo_freq_to_lsrk_freq(self, freq_topo_hz: np.ndarray, time_mjd: float, ra_deg: float, dec_deg: float) -> np.ndarray:
        ax = self.info.spectral_axis
        if ax is None:
            raise ValueError("LSRK conversion requires info.spectral_axis (RESTFREQ needed).")
        if ax.veldef != "RADIO":
            raise ValueError("This writer assumes VELDEF='RADIO' for LSRK conversion.")
        if ax.restfreq_hz <= 0:
            raise ValueError("LSRK conversion requires RESTFREQ > 0.")

        rest = ax.restfreq_hz * u.Hz
        f_topo = np.asarray(freq_topo_hz, dtype=np.float64).reshape(-1) * u.Hz
        v_topo = f_topo.to(u.m/u.s, equivalencies=u.doppler_radio(rest))
        vcor_kms = self._calc_vobs_kms_user(ra_deg=float(ra_deg), dec_deg=float(dec_deg), time_mjd=float(time_mjd))
        vcor = vcor_kms.to(u.m/u.s)
        v_lsrk = v_topo + vcor
        f_lsrk = v_lsrk.to(u.Hz, equivalencies=u.doppler_radio(rest))
        return np.asarray(f_lsrk.to_value(u.Hz), dtype=np.float64)


    @staticmethod
    def _normalize_veldef(value: Optional[str], specsys: str) -> Optional[str]:
        """Normalize user-friendly velocity definitions to standard 8-char VELDEF.

        Examples
        --------
        RADIO + TOPOCENT -> RADI-OBS
        OPTICAL + LSRK   -> OPTI-LSR
        RELATIVISTIC + HELIOCEN -> RELA-HEL
        """
        if value is None:
            return None
        v = str(value).strip().upper().replace('_', '-')
        if v == '':
            return None

        prefix_map = {
            'RADIO': 'RADI', 'RADI': 'RADI',
            'OPTICAL': 'OPTI', 'OPTI': 'OPTI',
            'RELATIVISTIC': 'RELA', 'RELA': 'RELA',
        }
        frame_codes = {'-OBS', '-LSR', '-HEL', '-BAR', '-GEO'}

        # Already in standard 8-char form? Keep if recognized.
        if len(v) == 8 and v[:4] in {'RADI', 'OPTI', 'RELA'} and v[4:] in frame_codes:
            return v

        if v not in prefix_map:
            raise ValueError(
                "veldef must be one of RADIO/OPTICAL/RELATIVISTIC or a standard 8-char form like RADI-LSR"
            )
        frame = VELDEF_FRAME_CODES.get(str(specsys).strip().upper())
        if frame is None:
            raise ValueError(f"Cannot derive VELDEF frame suffix from SPECSYS={specsys!r}")
        return prefix_map[v] + frame

    @staticmethod
    def _normalize_obsmode(val: str) -> str:
        """Normalize OBSMODE without over-constraining telescope-specific values."""
        s = str(val).strip().upper()
        if not s:
            raise ValueError('obsmode must be a non-empty string')
        return s

    @staticmethod
    def _validate_enum(name: str, val: str, options: Sequence[str]) -> str:
        if val not in options:
            raise ValueError(f"{name} must be one of {options}, got {val}")
        return val

    @staticmethod
    def _validate_str(name: str, val: str, maxlen: int) -> str:
        if not isinstance(val, str):
            raise TypeError(f"{name} must be str")
        if len(val) > maxlen:
            raise ValueError(f"{name} too long (len={len(val)} > {maxlen}). val={val}")
        return val

    # -------------------- public: add_row --------------------

    def add_history(self, key: str, value: Any = "") -> None:
        """Append a simple (key, value) entry to the file-level HISTORY.

        Notes
        -----
        - HISTORY is stored as a dedicated BinTable extension (KEY, VALUE) for round-trip safety.
        - This method is optional; you may also pass `history=` at initialization.
        """
        if self._history is None:
            self._history = {}
        if isinstance(self._history, dict):
            self._history[str(key)] = str(value)
            return
        if isinstance(self._history, list):
            self._history.append({str(key): str(value)})
            return
        # fallback: stringify existing and start a dict
        self._history = {"0": str(self._history), str(key): str(value)}

    def add_row(
        self,
        *,
        time_mjd: float,
        scanid: int,
        subscan: int,
        intgrp: int,
        obsmode: str,
        data: Union[np.ndarray, float],
        exposure_s: float,
        polariza: str,  # Physical polarization label: XX/YY/RR/LL/... (required)
        object_name: Optional[str] = None,
        flagrow: int = 0,
        flag: Optional[np.ndarray] = None,
        fdnum: int = 0,  # Feed / beam number
        ifnum: int = 0,  # IF / spectral-window number
        plnum: int = 0,  # Internal polarization index
        backend: Optional[str] = None,  # Spectrometer backend name (optional)
        sampler: Optional[str] = None,  # Spectrometer input / sampler name (optional)
        frontend: Optional[str] = None,  # Frontend / receiver name (optional)
        obsfreq_hz: Optional[float] = None,  # Signal-sideband sky/reference frequency [Hz]
        imagfreq_hz: Optional[float] = None,  # Image-sideband sky/reference frequency [Hz]
        lo1freq_hz: Optional[float] = None,  # 1st local oscillator frequency [Hz]
        lo2freq_hz: Optional[float] = None,  # 2nd local oscillator frequency [Hz]
        lo3freq_hz: Optional[float] = None,  # 3rd local oscillator frequency [Hz]
        sideband: Optional[str] = None,  # Final sky sideband represented by DATA: USB / LSB
        sb1: Optional[str] = None,  # 1st conversion sideband: USB / LSB
        sb2: Optional[str] = None,  # 2nd conversion sideband: USB / LSB
        sb3: Optional[str] = None,  # 3rd conversion sideband: USB / LSB
        calstat: str = "RAW",
        tempscal: Optional[str] = None,
        t_cal_k: float = np.nan,
        t_hot_k: float = np.nan,
        tsys_k: float = np.nan,
        tau0: float = np.nan,
        ra_deg: float = np.nan,
        dec_deg: float = np.nan,
        glon_deg: float = np.nan,
        glat_deg: float = np.nan,
        srcframe: Optional[str] = None,
        src_radesys: Optional[str] = None,
        src_equinox: Optional[float] = None,
        src_long_deg: Optional[float] = None,
        src_lat_deg: Optional[float] = None,
        scanframe: Optional[str] = None,
        scan_radesys: Optional[str] = None,
        scan_equinox: Optional[float] = None,
        scan_x_deg: Optional[float] = None,
        scan_y_deg: Optional[float] = None,
        az_center_deg: float = np.nan,
        el_center_deg: float = np.nan,
        az_enc_deg: float = np.nan,
        el_enc_deg: float = np.nan,
        boresight_az_deg: float = np.nan,
        boresight_el_deg: float = np.nan,
        beam_xoff_arcsec: float = np.nan,
        beam_yoff_arcsec: float = np.nan,
        beam_rot_deg: float = np.nan,
        az_cmd_deg: float = np.nan,
        el_cmd_deg: float = np.nan,
        corr_az_deg: float = np.nan,
        corr_el_deg: float = np.nan,
        calc_refr_deg: float = np.nan,
        v_frame_mps: float = np.nan,
        f_offset_hz: Optional[float] = None,
        tamb_c: float = np.nan,
        pressure_hpa: float = np.nan,
        humidity_pct: float = np.nan,
        wind_spd_mps: float = np.nan,
        wind_dir_deg: float = np.nan,
        freq_hz: Optional[np.ndarray] = None,
        restfreq_hz: Optional[float] = None,
        crval1_hz: Optional[float] = None,
        cdelt1_hz: Optional[float] = None,
        crpix1: Optional[float] = None,
        ctype1: Optional[str] = None,
        cunit1: Optional[str] = None,
        specsys: Optional[str] = None,
        veldef: Optional[str] = None,
    ) -> None:
        # Check finite MJD, positive exposure, non-negative IDs, valid enums... (Same as before)
        if not np.isfinite(time_mjd):
            raise ValueError("time_mjd must be finite float (MJD UTC days)")
        if exposure_s <= 0:
            raise ValueError("exposure_s must be > 0")
        if scanid < 0 or subscan < 0 or intgrp < 0:
            raise ValueError("scanid/subscan/intgrp must be >= 0")
        
        if fdnum < 0 or ifnum < 0 or plnum < 0:
            raise ValueError("fdnum/ifnum/plnum must be >= 0")

        # Validate enums and normalize key labels.
        obsmode = self._normalize_obsmode(obsmode)
        calstat = self._validate_enum("calstat", str(calstat), CALSTAT_VALUES)

        # Physical polarization is a required, first-class row attribute.
        if not _is_meaningful_str(polariza):
            raise ValueError("polariza is required and must be a non-empty string")
        polariza = str(polariza).strip().upper()
        if polariza not in POLARIZA_VALUES:
            raise ValueError(f"polariza must be one of {POLARIZA_VALUES}, got {polariza}")

        # Optional provenance strings: keep None if unknown so optional-column logic works.
        backend = (str(backend).strip() if _is_meaningful_str(backend) else None)
        sampler = (str(sampler).strip() if _is_meaningful_str(sampler) else None)
        frontend = (str(frontend).strip() if _is_meaningful_str(frontend) else None)

        # Optional sideband labels. Normalize to upper case when present.
        sideband = (str(sideband).strip().upper() if _is_meaningful_str(sideband) else None)
        sb1 = (str(sb1).strip().upper() if _is_meaningful_str(sb1) else None)
        sb2 = (str(sb2).strip().upper() if _is_meaningful_str(sb2) else None)
        sb3 = (str(sb3).strip().upper() if _is_meaningful_str(sb3) else None)
        for name, sb in (("sideband", sideband), ("sb1", sb1), ("sb2", sb2), ("sb3", sb3)):
            if sb is not None and sb not in SIDEBAND_VALUES:
                raise ValueError(f"{name} must be one of {SIDEBAND_VALUES}, got {sb}")

        # Heterodyne sideband context: if any LO/sideband/image metadata are present,
        # the row should also define the signal-sideband reference frequency and final
        # sideband interpretation of DATA. FRONTEND is strongly recommended but not
        # enforced because some existing pipelines may not know it at row-build time.
        heterodyne_context_active = any([
            _is_meaningful_scalar(obsfreq_hz),
            _is_meaningful_scalar(imagfreq_hz),
            _is_meaningful_scalar(lo1freq_hz),
            _is_meaningful_scalar(lo2freq_hz),
            _is_meaningful_scalar(lo3freq_hz),
            _is_meaningful_str(sideband),
            _is_meaningful_str(sb1),
            _is_meaningful_str(sb2),
            _is_meaningful_str(sb3),
        ])
        if heterodyne_context_active:
            if not _is_meaningful_scalar(obsfreq_hz):
                raise ValueError("heterodyne LO/sideband metadata require obsfreq_hz")
            if not _is_meaningful_str(sideband):
                raise ValueError("heterodyne LO/sideband metadata require sideband ('USB' or 'LSB')")

        ra_deg, dec_deg, glon_deg, glat_deg = self._fill_radec_glon_glat(
            ra_deg=ra_deg, dec_deg=dec_deg, glon_deg=glon_deg, glat_deg=glat_deg
        )

        # Effective per-row spectral WCS values.
        ax = self.info.spectral_axis
        eff_restfreq_hz = restfreq_hz if restfreq_hz is not None else (ax.restfreq_hz if ax else np.nan)
        if not np.isfinite(eff_restfreq_hz) or float(eff_restfreq_hz) <= 0:
            eff_restfreq_hz = np.nan

        eff_crval1_hz = crval1_hz if crval1_hz is not None else (ax.crval1_hz if ax else 0.0)
        eff_cdelt1_hz = cdelt1_hz if cdelt1_hz is not None else (ax.cdelt1_hz if ax else 1.0)
        eff_crpix1 = crpix1 if crpix1 is not None else (ax.crpix1 if ax else 1.0)
        eff_ctype1 = str(ctype1 if ctype1 is not None else (ax.ctype1 if ax else "FREQ")).strip().upper()
        eff_cunit1 = str(cunit1 if cunit1 is not None else (ax.cunit1 if ax else "Hz")).strip()
        eff_specsys = str(specsys if specsys is not None else (ax.specsys if ax else "TOPOCENT")).strip().upper()
        eff_veldef = veldef if veldef is not None else (ax.veldef if ax else None)
        eff_veldef = None if eff_veldef is None else str(eff_veldef).strip().upper()

        # Context-required columns: only mandatory when the row carries velocity semantics.
        needs_velocity_context = (
            eff_ctype1 != "FREQ"
            or eff_specsys == "LSRK"
            or np.isfinite(v_frame_mps)
        )
        if needs_velocity_context:
            if not np.isfinite(eff_restfreq_hz) or float(eff_restfreq_hz) <= 0:
                raise ValueError("velocity-context rows require RESTFREQ > 0")
            if eff_veldef is None:
                raise ValueError("velocity-context rows require VELDEF")
            eff_veldef = self._normalize_veldef(eff_veldef, eff_specsys)
        elif eff_veldef is not None and eff_veldef != "":
            eff_veldef = self._normalize_veldef(eff_veldef, eff_specsys)

        eff_specsys = self._validate_enum("specsys", eff_specsys, SPECSYS_VALUES)

        # Block-optional groups: activate only when at least one field is meaningful.
        source_block_active = any([
            _is_meaningful_str(srcframe),
            _is_meaningful_str(src_radesys),
            _is_meaningful_scalar(src_equinox),
            _is_meaningful_scalar(src_long_deg),
            _is_meaningful_scalar(src_lat_deg),
        ])
        scan_block_active = any([
            _is_meaningful_str(scanframe),
            _is_meaningful_str(scan_radesys),
            _is_meaningful_scalar(scan_equinox),
            _is_meaningful_scalar(scan_x_deg),
            _is_meaningful_scalar(scan_y_deg),
        ])
        pointing_block_active = any([
            _is_meaningful_scalar(az_cmd_deg),
            _is_meaningful_scalar(el_cmd_deg),
            _is_meaningful_scalar(corr_az_deg),
            _is_meaningful_scalar(corr_el_deg),
            _is_meaningful_scalar(calc_refr_deg),
        ])
        weather_block_active = any([
            _is_meaningful_scalar(tamb_c),
            _is_meaningful_scalar(pressure_hpa),
            _is_meaningful_scalar(humidity_pct),
            _is_meaningful_scalar(wind_spd_mps),
            _is_meaningful_scalar(wind_dir_deg),
        ])

        if source_block_active:
            if not _is_meaningful_str(srcframe):
                raise ValueError(
                    "source block requires explicit srcframe; do not rely on an implicit default like RADEC"
                )
            eff_srcframe = str(srcframe).strip().upper()
            eff_srcframe = self._validate_enum("srcframe", eff_srcframe, COORDTYPE_VALUES)
            if eff_srcframe == "RADEC":
                eff_src_radesys = (
                    str(src_radesys).strip().upper()
                    if _is_meaningful_str(src_radesys)
                    else str(self.info.src_radesys).upper()
                )
                eff_src_radesys = self._validate_enum("src_radesys", eff_src_radesys, RADESYS_VALUES)
                eff_src_equinox = (
                    float(src_equinox)
                    if _is_meaningful_scalar(src_equinox)
                    else float(self.info.src_equinox)
                )
            else:
                if _is_meaningful_str(src_radesys) or _is_meaningful_scalar(src_equinox):
                    raise ValueError(
                        "src_radesys/src_equinox are only valid when srcframe='RADEC'"
                    )
                eff_src_radesys = None
                eff_src_equinox = np.nan
        else:
            eff_srcframe = None
            eff_src_radesys = None
            eff_src_equinox = np.nan

        if scan_block_active:
            if not _is_meaningful_str(scanframe):
                raise ValueError(
                    "scan block requires explicit scanframe; do not rely on an implicit default like RADEC"
                )
            eff_scanframe = str(scanframe).strip().upper()
            eff_scanframe = self._validate_enum("scanframe", eff_scanframe, COORDTYPE_VALUES)
            if eff_scanframe == "RADEC":
                eff_scan_radesys = (
                    str(scan_radesys).strip().upper()
                    if _is_meaningful_str(scan_radesys)
                    else str(self.info.radesys).upper()
                )
                eff_scan_radesys = self._validate_enum("scan_radesys", eff_scan_radesys, RADESYS_VALUES)
                eff_scan_equinox = (
                    float(scan_equinox)
                    if _is_meaningful_scalar(scan_equinox)
                    else float(self.info.equinox)
                )
            else:
                if _is_meaningful_str(scan_radesys) or _is_meaningful_scalar(scan_equinox):
                    raise ValueError(
                        "scan_radesys/scan_equinox are only valid when scanframe='RADEC'"
                    )
                eff_scan_radesys = None
                eff_scan_equinox = np.nan
        else:
            eff_scanframe = None
            eff_scan_radesys = None
            eff_scan_equinox = np.nan

        # Buffer addition logic
        b = self._buf
        dateobs_iso = Time(time_mjd, format="mjd", scale="utc").isot
        b.TIME.append(0.0)
        b.MJD.append(float(time_mjd))
        b.DATEOBS.append(dateobs_iso)
        b.SCAN.append(int(scanid))
        b.SUBSCAN.append(int(subscan))
        b.INTGRP.append(int(intgrp))
        b.OBJECT.append(str(object_name if object_name else self.info.object_name))
        b.OBSMODE.append(str(obsmode))
        b.FDNUM.append(int(fdnum))
        b.IFNUM.append(int(ifnum))
        b.PLNUM.append(int(plnum))
        b.POLARIZA.append(str(polariza))

        # Receiver / spectrometer provenance.
        b.FRONTEND.append(frontend)
        b.BACKEND.append(backend)
        b.SAMPLER.append(sampler)

        # Heterodyne LO / sideband metadata.  These are optional row attributes,
        # but when any of them are present the signal-sideband reference
        # frequency (OBSFREQ) and the final sideband interpretation (SIDEBAND)
        # become mandatory for that row.
        b.OBSFREQ.append(float(obsfreq_hz) if _is_meaningful_scalar(obsfreq_hz) else np.nan)
        b.IMAGFREQ.append(float(imagfreq_hz) if _is_meaningful_scalar(imagfreq_hz) else np.nan)
        b.LO1FREQ.append(float(lo1freq_hz) if _is_meaningful_scalar(lo1freq_hz) else np.nan)
        b.LO2FREQ.append(float(lo2freq_hz) if _is_meaningful_scalar(lo2freq_hz) else np.nan)
        b.LO3FREQ.append(float(lo3freq_hz) if _is_meaningful_scalar(lo3freq_hz) else np.nan)
        b.SIDEBAND.append(sideband)
        b.SB1.append(sb1)
        b.SB2.append(sb2)
        b.SB3.append(sb3)

        b.EXPOSURE.append(float(exposure_s))
        b.CALSTAT.append(str(calstat))
        # Temperature scale (SDFITS/CASA). Chopper-wheel calibrated data are Ta*.
        try:
            if tempscal is not None:
                b.TEMPSCAL.append(normalize_tempscal(tempscal))
            else:
                cs = str(calstat).upper()
                if "TMB" in cs or "TR" in cs:
                    b.TEMPSCAL.append("TR*")
                else:
                    b.TEMPSCAL.append("TA*")
        except Exception:
            b.TEMPSCAL.append("TA*")
        
        b.TCAL.append(float(t_cal_k))
        b.THOT.append(float(t_hot_k))
        b.TSYS.append(float(tsys_k))
        b.TAU0.append(float(tau0))
        b.FLAGROW.append(int(flagrow))
        b.RA.append(float(ra_deg))
        b.DEC.append(float(dec_deg))
        b.GLON.append(float(glon_deg))
        b.GLAT.append(float(glat_deg))
        eff_az = float(az_center_deg) if _is_meaningful_scalar(az_center_deg) else float(az_enc_deg)
        eff_el = float(el_center_deg) if _is_meaningful_scalar(el_center_deg) else float(el_enc_deg)
        b.AZIMUTH.append(eff_az)
        b.ELEVATIO.append(eff_el)
        b.BORE_AZ.append(float(boresight_az_deg) if _is_meaningful_scalar(boresight_az_deg) else np.nan)
        b.BORE_EL.append(float(boresight_el_deg) if _is_meaningful_scalar(boresight_el_deg) else np.nan)
        b.BEAMXOFF.append(float(beam_xoff_arcsec) if _is_meaningful_scalar(beam_xoff_arcsec) else np.nan)
        b.BEAMYOFF.append(float(beam_yoff_arcsec) if _is_meaningful_scalar(beam_yoff_arcsec) else np.nan)
        b.BEAMROT.append(float(beam_rot_deg) if _is_meaningful_scalar(beam_rot_deg) else np.nan)
        
        # Source-coordinate block: emit only if at least one source coordinate field is meaningful.
        b.SRCFRAME.append(eff_srcframe)
        b.SRCRDSYS.append(eff_src_radesys)
        b.SRCEQNX.append(float(eff_src_equinox) if np.isfinite(eff_src_equinox) else np.nan)
        b.SRC_LONG.append(float(src_long_deg) if _is_meaningful_scalar(src_long_deg) else np.nan)
        b.SRC_LAT.append(float(src_lat_deg) if _is_meaningful_scalar(src_lat_deg) else np.nan)

        # Scan-coordinate block: emit only if scan offsets / scan-frame metadata exist.
        b.SCANFRAM.append(eff_scanframe)
        b.SCANRDSYS.append(eff_scan_radesys)
        b.SCANEQNX.append(float(eff_scan_equinox) if np.isfinite(eff_scan_equinox) else np.nan)
        b.SCANX.append(float(scan_x_deg) if _is_meaningful_scalar(scan_x_deg) else np.nan)
        b.SCANY.append(float(scan_y_deg) if _is_meaningful_scalar(scan_y_deg) else np.nan)

        # Pointing-diagnostic block (optional as a coherent group).
        b.AZ_CMD.append(float(az_cmd_deg) if pointing_block_active and _is_meaningful_scalar(az_cmd_deg) else np.nan)
        b.EL_CMD.append(float(el_cmd_deg) if pointing_block_active and _is_meaningful_scalar(el_cmd_deg) else np.nan)
        b.CORR_AZ.append(float(corr_az_deg) if pointing_block_active and _is_meaningful_scalar(corr_az_deg) else np.nan)
        b.CORR_EL.append(float(corr_el_deg) if pointing_block_active and _is_meaningful_scalar(corr_el_deg) else np.nan)
        b.CALC_REFR.append(float(calc_refr_deg) if pointing_block_active and _is_meaningful_scalar(calc_refr_deg) else np.nan)

        # Optional scalar metadata.
        b.VFRAME.append(float(v_frame_mps) if _is_meaningful_scalar(v_frame_mps) else np.nan)
        b.FOFFSET.append(float(f_offset_hz) if _is_meaningful_scalar(f_offset_hz) else np.nan)

        # Weather block (optional as a coherent group).
        # Store SDFITS-like standard units where practical:
        #   TAMBIENT -> K (input here is degC)
        #   PRESSURE -> mmHg (input here is hPa)
        #   HUMIDITY -> fraction 0..1 (input here is percent)
        tamb_k = (float(tamb_c) + 273.15) if weather_block_active and _is_meaningful_scalar(tamb_c) else np.nan
        pressure_mmhg = (float(pressure_hpa) * 0.750061683) if weather_block_active and _is_meaningful_scalar(pressure_hpa) else np.nan
        humidity_frac = (float(humidity_pct) / 100.0) if weather_block_active and _is_meaningful_scalar(humidity_pct) else np.nan
        b.TAMBIENT.append(tamb_k)
        b.PRESSURE.append(pressure_mmhg)
        b.HUMIDITY.append(humidity_frac)
        b.WINDSPD.append(float(wind_spd_mps) if weather_block_active and _is_meaningful_scalar(wind_spd_mps) else np.nan)
        b.WINDDIR.append(float(wind_dir_deg) if weather_block_active and _is_meaningful_scalar(wind_dir_deg) else np.nan)

        # WCS Parameters (RESTFREQ / VELDEF are context-required, not unconditional).
        b.RESTFREQ.append(float(eff_restfreq_hz) if np.isfinite(eff_restfreq_hz) else np.nan)
        b.CRVAL1.append(float(eff_crval1_hz))
        b.CDELT1.append(float(eff_cdelt1_hz))
        b.CRPIX1.append(float(eff_crpix1))
        b.CTYPE1.append(str(eff_ctype1))
        b.CUNIT1.append(str(eff_cunit1))
        b.SPECSYS.append(str(eff_specsys))
        b.VELDEF.append(None if eff_veldef is None else str(eff_veldef))

        # Arrays
        b.append_spectrum(data, flag)

        if self.store_freq_column:
            if freq_hz is None:
                # Need explicit WCS for freq generation in VLA mode
                # Assume provided/default WCS is valid for generating the freq array
                c_v = b.CRVAL1[-1]
                c_d = b.CDELT1[-1]
                c_p = b.CRPIX1[-1]
                
                # Determine length from data just added
                row_len = b._DATA[-1].size
                
                pix = (np.arange(row_len, dtype=np.float64) + 1.0)
                freq_hz = c_v + (pix - c_p) * c_d

            # If the stored per-row spectral frame is LSRK, convert the generated
            # topocentric frequency grid into the LSRK frequency grid.
            if eff_specsys == "LSRK":
                freq_hz = self._topo_freq_to_lsrk_freq(freq_hz, time_mjd=float(time_mjd), ra_deg=float(ra_deg), dec_deg=float(dec_deg))

            b.append_freq(freq_hz)

        if self.chunk_size is not None and len(b) >= self.chunk_size:
            self._flush_part()

    # -------------------- FITS construction --------------------
    def _build_primary_hdu(self, mjdstart: float, mjdend: float) -> fits.PrimaryHDU:
        # Same as before
        hdr = fits.Header()
        hdr["SIMPLE"] = True
        hdr["BITPIX"] = 8
        hdr["NAXIS"] = 0
        hdr["EXTEND"] = True
        hdr["TELESCOP"] = (self.info.telescope, "Telescope name")
        hdr["OBSERVER"] = (self.info.observer, "Observer")
        hdr["PROJID"]   = (self.info.project,  "Project/Proposal ID")
        hdr["OBJECT"]   = (self.info.object_name, "Representative target name")
        hdr["DATE"]     = (Time.now().utc.isot, "File creation time (UTC)")
        hdr["TIMESYS"]  = ("UTC", "Time system used in TIME column")
        hdr["MJDSTART"] = (float(mjdstart), "Start MJD(UTC)")
        hdr["MJDEND"]   = (float(mjdend), "End MJD(UTC)")
        hdr["SWNAME"]   = (SWNAME, "Software name")
        hdr["SWVER"]    = (__version__, "Software version")
        hdr["ORIGIN"]   = (SWNAME, "Writer ID")
        self._write_site_keywords(hdr)
        self._write_spectral_keywords(hdr)
        return fits.PrimaryHDU(header=hdr)

    def _write_site_keywords(self, hdr: fits.Header) -> None:
        x, y, z = self.site.obsgeo_xyz_m()
        hdr["SITELAT"]  = (float(self.site.lat_deg), "Site latitude (deg)")
        hdr["SITELON"]  = (float(self.site.lon_deg), "Site longitude (deg, east+)")
        hdr["SITELONG"] = (float(self.site.lon_deg), "Site longitude (deg, east+) (alias)")
        hdr["SITEELEV"] = (float(self.site.elev_m),  "Site elevation (m)")
        hdr["OBSGEO-X"] = (float(x), "Observatory geocenter X (m)")
        hdr["OBSGEO-Y"] = (float(y), "Observatory geocenter Y (m)")
        hdr["OBSGEO-Z"] = (float(z), "Observatory geocenter Z (m)")

    def _write_spectral_keywords(self, hdr: fits.Header) -> None:
        if self.info.spectral_axis is None:
            return
        ax = self.info.spectral_axis
        hdr["WCSAXES"] = (1, "Number of WCS axes (spectral only)")
        hdr["CTYPE1"]  = (ax.ctype1, "Spectral axis type (default)")
        hdr["CUNIT1"]  = (ax.cunit1, "Spectral axis unit")
        hdr["CRVAL1"]  = (float(ax.crval1_hz), "Ref freq at CRPIX1 (Hz) (default)")
        hdr["CDELT1"]  = (float(ax.cdelt1_hz), "Channel spacing (Hz) (default)")
        hdr["CRPIX1"]  = (float(ax.crpix1), "Ref pixel (1-based) (default)")
        hdr["SPECSYS"]  = (str(ax.specsys), "Spectral ref frame")
        hdr["SSYSOBS"]  = (str(ax.ssysobs), "Observer frame")
        hdr["REFCHAN"]  = (int(ax.refchan), "Reference channel index (1-based)")

        # Context-required defaults: only write them when they carry actual
        # velocity / line semantics.  This avoids meaningless header defaults.
        velocity_context = (
            str(ax.ctype1).upper() != "FREQ"
            or str(ax.specsys).upper() == "LSRK"
            or float(ax.restfreq_hz) > 0
        )
        if float(ax.restfreq_hz) > 0:
            hdr["RESTFREQ"] = (float(ax.restfreq_hz), "Rest frequency (Hz) (default)")
            hdr["RESTFRQ"]  = (float(ax.restfreq_hz), "Rest frequency (Hz) (alias)")
        if velocity_context:
            hdr["VELREF"]   = (int(ax.velref()), "AIPS/casacore VELREF")
            hdr["VELDEF"]   = (str(ax.veldef), "Velocity definition")

    def _build_table_hdu(self) -> fits.BinTableHDU:
        b = self._buf
        if len(b) == 0:
            raise RuntimeError("No data to write.")

        def _to_str_list(xs: Sequence[Any]) -> list[str]:
            """Convert always-present string columns to plain strings for FITS writing."""
            out: list[str] = []
            for x in xs:
                if x is None:
                    out.append("")
                elif isinstance(x, (bytes, bytearray)):
                    out.append(x.decode("utf-8", errors="ignore"))
                else:
                    out.append(str(x))
            return out

        # ---------------------------------------------------------------------
        # Build a DataFrame in a stable column order.  ALWAYS columns are added
        # unconditionally; OPTIONAL columns only when meaningful; BLOCK-OPTIONAL
        # groups are emitted as a whole if any member carries information.
        # ---------------------------------------------------------------------
        df = pd.DataFrame({
            # Basic row identity / observation info (ALWAYS)
            "TIME":     np.array(b.TIME, dtype=np.float64),
            "MJD":      np.array(b.MJD, dtype=np.float64),
            "DATE-OBS": _to_str_list(b.DATEOBS),
            "DATEOBS":  _to_str_list(b.DATEOBS),
            "TIMESTAMP": _to_str_list(b.DATEOBS),
            "SCAN":     np.array(b.SCAN, dtype=np.int32),
            "SUBSCAN":  np.array(b.SUBSCAN, dtype=np.int32),
            "INTGRP":   np.array(b.INTGRP, dtype=np.int32),
            "OBJECT":   _to_str_list(b.OBJECT),
            "OBSMODE":  _to_str_list(b.OBSMODE),
            "FDNUM":    np.array(b.FDNUM, dtype=np.int32),
            "IFNUM":    np.array(b.IFNUM, dtype=np.int32),
            "PLNUM":    np.array(b.PLNUM, dtype=np.int32),
            "POLARIZA": _to_str_list(b.POLARIZA),
            "EXPOSURE": np.array(b.EXPOSURE, dtype=np.float32),
            "CALSTAT":  _to_str_list(b.CALSTAT),
            "TEMPSCAL": _to_str_list(getattr(b, "TEMPSCAL", ["TA*"] * len(b.CALSTAT))),
            "FLAGROW":  np.array(b.FLAGROW, dtype=np.int32),

            # Coordinates (ALWAYS)
            "RA":       np.array(b.RA, dtype=np.float64),
            "DEC":      np.array(b.DEC, dtype=np.float64),
            "GLON":     np.array(b.GLON, dtype=np.float64),
            "GLAT":     np.array(b.GLAT, dtype=np.float64),

            # Spectral WCS core (ALWAYS)
            "CRVAL1":   np.array(b.CRVAL1, dtype=np.float64),
            "CDELT1":   np.array(b.CDELT1, dtype=np.float64),
            "CRPIX1":   np.array(b.CRPIX1, dtype=np.float64),
            "CTYPE1":   _to_str_list(b.CTYPE1),
            "CUNIT1":   _to_str_list(b.CUNIT1),
            "SPECSYS":  _to_str_list(b.SPECSYS),
        })

        # Optional scalar columns: include only if at least one row has a meaningful value.
        optional_float32 = {
            "TCAL": b.TCAL,
            "THOT": b.THOT,
            "TSYS": b.TSYS,
            "TAU0": b.TAU0,
            "VFRAME": b.VFRAME,
            "FOFFSET": b.FOFFSET,
        }
        optional_float64 = {
            "AZIMUTH": b.AZIMUTH,
            "ELEVATIO": b.ELEVATIO,
            "BORE_AZ": b.BORE_AZ,
            "BORE_EL": b.BORE_EL,
            "BEAMXOFF": b.BEAMXOFF,
            "BEAMYOFF": b.BEAMYOFF,
            "BEAMROT": b.BEAMROT,
            "IMAGFREQ": b.IMAGFREQ,
            "LO1FREQ": b.LO1FREQ,
            "LO2FREQ": b.LO2FREQ,
            "LO3FREQ": b.LO3FREQ,
        }
        optional_object = {
            "FRONTEND": b.FRONTEND,
            "BACKEND": b.BACKEND,
            "SAMPLER": b.SAMPLER,
            "SB1": b.SB1,
            "SB2": b.SB2,
            "SB3": b.SB3,
        }

        for name, values in optional_float32.items():
            if _series_has_meaningful(values):
                df[name] = np.array(values, dtype=np.float32)
        for name, values in optional_float64.items():
            if _series_has_meaningful(values):
                df[name] = np.array(values, dtype=np.float64)
        for name, values in optional_object.items():
            if _series_has_meaningful(values):
                df[name] = pd.Series(list(values), dtype="object")

        # Context-required columns: emit only when the overall dataset carries
        # velocity / line semantics.  If emitted, add the whole column.
        velocity_context_active = (
            any(str(x).strip().upper() != "FREQ" for x in b.CTYPE1)
            or any(str(x).strip().upper() == "LSRK" for x in b.SPECSYS)
            or _series_has_meaningful(b.VFRAME)
        )
        if velocity_context_active:
            df["RESTFREQ"] = np.array(b.RESTFREQ, dtype=np.float64)
            df["VELDEF"] = pd.Series(list(b.VELDEF), dtype="object")

        # Heterodyne LO / sideband context.  OBSFREQ and SIDEBAND are the
        # minimum pair required to interpret DATA as a specific sky-sideband
        # spectrum.  Additional LO / image-sideband information remains optional.
        heterodyne_context_active = (
            _series_has_meaningful(b.OBSFREQ)
            or _series_has_meaningful(b.IMAGFREQ)
            or _series_has_meaningful(b.LO1FREQ)
            or _series_has_meaningful(b.LO2FREQ)
            or _series_has_meaningful(b.LO3FREQ)
            or _series_has_meaningful(b.SIDEBAND)
            or _series_has_meaningful(b.SB1)
            or _series_has_meaningful(b.SB2)
            or _series_has_meaningful(b.SB3)
        )
        if heterodyne_context_active:
            if not _series_has_meaningful(b.OBSFREQ):
                raise ValueError("heterodyne-context dataset requires OBSFREQ")
            if not _series_has_meaningful(b.SIDEBAND):
                raise ValueError("heterodyne-context dataset requires SIDEBAND")
            df["OBSFREQ"] = np.array(b.OBSFREQ, dtype=np.float64)
            df["SIDEBAND"] = pd.Series(list(b.SIDEBAND), dtype="object")

        # Block-optional groups: if any member has data, emit the whole block.
        if _block_has_meaningful(b, sorted(BLOCK_OPTIONAL["SOURCE"])):
            df["SRCFRAME"] = pd.Series(list(b.SRCFRAME), dtype="object")
            df["SRCRDSYS"] = pd.Series(list(b.SRCRDSYS), dtype="object")
            df["SRCEQNX"] = np.array(b.SRCEQNX, dtype=np.float64)
            df["SRC_LONG"] = np.array(b.SRC_LONG, dtype=np.float64)
            df["SRC_LAT"] = np.array(b.SRC_LAT, dtype=np.float64)

        if _block_has_meaningful(b, sorted(BLOCK_OPTIONAL["SCAN"])):
            df["SCANFRAM"] = pd.Series(list(b.SCANFRAM), dtype="object")
            df["SCANRDSYS"] = pd.Series(list(b.SCANRDSYS), dtype="object")
            df["SCANEQNX"] = np.array(b.SCANEQNX, dtype=np.float64)
            df["SCANX"] = np.array(b.SCANX, dtype=np.float64)
            df["SCANY"] = np.array(b.SCANY, dtype=np.float64)

        if _block_has_meaningful(b, sorted(BLOCK_OPTIONAL["POINTING"])):
            df["AZ_CMD"] = np.array(b.AZ_CMD, dtype=np.float64)
            df["EL_CMD"] = np.array(b.EL_CMD, dtype=np.float64)
            df["CORR_AZ"] = np.array(b.CORR_AZ, dtype=np.float64)
            df["CORR_EL"] = np.array(b.CORR_EL, dtype=np.float64)
            df["CALC_REFR"] = np.array(b.CALC_REFR, dtype=np.float64)

        if _block_has_meaningful(b, sorted(BLOCK_OPTIONAL["WEATHER"])):
            df["TAMBIENT"] = np.array(b.TAMBIENT, dtype=np.float32)
            df["PRESSURE"] = np.array(b.PRESSURE, dtype=np.float32)
            df["HUMIDITY"] = np.array(b.HUMIDITY, dtype=np.float32)
            df["WINDSPD"] = np.array(b.WINDSPD, dtype=np.float32)
            df["WINDDIR"] = np.array(b.WINDDIR, dtype=np.float32)

        # Fixed widths for SDFITS-style string columns (stability + readability)
        string_widths = {
            "DATE-OBS": 26,
            "DATEOBS": 26,
            "TIMESTAMP": 26,
            "OBJECT": int(self.STR_OBJECT),
            "OBSMODE": int(self.STR_OBSMODE),
            "CALSTAT": int(self.STR_CALSTAT),
            "POLARIZA": 8,
            "FRONTEND": 32,
            "BACKEND": 32,
            "SAMPLER": 32,
            "SIDEBAND": 8,
            "SB1": 8,
            "SB2": 8,
            "SB3": 8,
            "SRCFRAME": int(self.STR_FRAME),
            "SRCRDSYS": int(self.STR_FRAME),
            "SCANFRAM": int(self.STR_FRAME),
            "SCANRDSYS": int(self.STR_FRAME),
            "CTYPE1": int(self.STR_CTYPE),
            "CUNIT1": int(self.STR_CUNIT),
            "SPECSYS": int(self.STR_SPECSYS),
            "VELDEF": int(self.STR_VELDEF),
        }

        units = {
            "TIME": "s",
            "MJD": "d",
            "EXPOSURE": "s",
            "TCAL": "K",
            "THOT": "K",
            "TSYS": "K",
            "RA": "deg",
            "DEC": "deg",
            "GLON": "deg",
            "GLAT": "deg",
            "AZIMUTH": "deg",
            "ELEVATIO": "deg",
            "RESTFREQ": "Hz",
            "OBSFREQ": "Hz",
            "IMAGFREQ": "Hz",
            "LO1FREQ": "Hz",
            "LO2FREQ": "Hz",
            "LO3FREQ": "Hz",
            "SRCEQNX": "year",
            "SCANEQNX": "year",
            "SRC_LONG": "deg",
            "SRC_LAT": "deg",
            "SCANX": "deg",
            "SCANY": "deg",
            "AZ_CMD": "deg",
            "EL_CMD": "deg",
            "CORR_AZ": "deg",
            "CORR_EL": "deg",
            "CALC_REFR": "deg",
            "VFRAME": "m/s",
            "FOFFSET": "Hz",
            "TAMBIENT": "K",
            "PRESSURE": "mmHg",
            "WINDSPD": "m/s",
            "WINDDIR": "deg",
            "FREQ": "Hz",
        }

        # CRVAL1/CDELT1 are scalar row metadata, but FITS TUNITn is attached per
        # column (not per row).  Therefore the dataset must use one consistent
        # CUNIT1 if we want to assign a physically correct unit to those columns.
        cunit1_values = [str(x).strip() for x in b.CUNIT1 if _is_meaningful_str(x)]
        if cunit1_values:
            cunit1_unique = {x.upper(): x for x in cunit1_values}
            if len(cunit1_unique) != 1:
                raise ValueError(
                    "CUNIT1 varies across rows; FITS cannot assign a single TUNIT to CRVAL1/CDELT1. "
                    "Please split the dataset or normalize CUNIT1 before writing."
                )
            dataset_cunit1 = builtins.next(iter(cunit1_unique.values()))
            units["CRVAL1"] = dataset_cunit1
            units["CDELT1"] = dataset_cunit1

        extra_vec = {"FREQ": b._FREQ} if self.store_freq_column else None

        hdu, nchan_hint = build_single_dish_table_hdu(
            table=df,
            data=b._DATA,
            spectrum_column=self.SPECTRUM_COLUMN,
            include_flag=True,
            flag=b._FLAG,
            extra_vector_columns=extra_vec,
            units=units,
            string_widths=string_widths,
            normalize_columns=True,
            prefer_float64_vectors=True,
            allow_vla_vectors=True,
            bunit="K",
            extname="SINGLE DISH",
        )
        hdu.name = "SINGLE DISH"

        eh = hdu.header
        eh["EXTNAME"]  = ("SINGLE DISH", "Single dish table (SDFITS-like)")
        eh["TIMESYS"]  = ("UTC", "Time system for DATE-OBS/TIME")
        eh["TELESCOP"] = (self.info.telescope, "Telescope name")
        eh["OBSERVER"] = (self.info.observer,  "Observer")
        eh["PROJID"]   = (self.info.project,   "Project/Proposal ID")
        eh["OBJECT"]   = (self.info.object_name, "Representative target name")
        if len(b.DATEOBS) > 0:
            eh["DATE-OBS"] = (str(b.DATEOBS[0]), "Representative row start time (UTC)")
        eh["SWNAME"]   = (SWNAME, "Software name")
        eh["SWVER"]    = (__version__, "Software version")

        # Determine whether DATA is VLA (P/Q) or fixed-length
        is_var = False
        for i in range(1, int(eh.get("TFIELDS", 0)) + 1):
            ttype = (eh.get(f"TTYPE{i}", "") or "").strip().upper()
            if ttype == self.SPECTRUM_COLUMN.upper():
                tform = (eh.get(f"TFORM{i}", "") or "").strip().upper()
                is_var = tform.startswith(("P", "Q"))
                break

        eh["NCHAN"]    = (int(nchan_hint), "Channel count (fixed) or max length (VLA)")
        eh["NCHANSEL"] = (int(nchan_hint), "Channel count (fixed) or max length (VLA)")

        mjd0 = float(builtins.min(b.MJD))
        mjd1 = float(builtins.max(b.MJD))
        eh["TIMESYS"]  = ("UTC", "Time system for TIME column")
        eh["MJDSTART"] = (mjd0, "Start MJD(UTC)")
        eh["MJDEND"]   = (mjd1, "End MJD(UTC)")

        eh["RADESYS"]  = (self.info.radesys, "RA/Dec frame")
        eh["EQUINOX"]  = (float(self.info.equinox), "Equinox (FK4/FK5)")
        eh["HIERARCH SRC_RADESYS"] = (self.info.src_radesys, "SRC RA/Dec frame")
        eh["HIERARCH SRC_EQUINOX"] = (float(self.info.src_equinox), "SRC equinox (FK4/FK5)")

        self._write_site_keywords(eh)

        if np.isfinite(self.info.bmaj_deg):
            eh["BMAJ"] = (float(self.info.bmaj_deg), "Beam major axis (deg)")
        if np.isfinite(self.info.bmin_deg):
            eh["BMIN"] = (float(self.info.bmin_deg), "Beam minor axis (deg)")
        if np.isfinite(self.info.bpa_deg):
            eh["BPA"]  = (float(self.info.bpa_deg), "Beam position angle (deg)")

        if np.isfinite(self.info.eff.beameff):
            eh["BEAMEFF"] = (float(self.info.eff.beameff), "Main-beam efficiency eta_mb")
        eh["TEMPSCAL"] = ("TA*", "Temperature scale of DATA (TA* or TR*)")
        if np.isfinite(self.info.eff.apereff):
            eh["APEREFF"] = (float(self.info.eff.apereff), "Aperture efficiency eta_a")
        if np.isfinite(self.info.eff.mooneff):
            eh["MOONEFF"] = (float(self.info.eff.mooneff), "Moon efficiency")
        if np.isfinite(self.info.eff.suneff):
            eh["SUNEFF"]  = (float(self.info.eff.suneff), "Sun efficiency")
        eh["EFFSTAT"] = (str(self.info.eff.effstat), "Efficiency status")

        eh["HIERARCH REFR_INC"] = (bool(self.info.refr_included_in_corr), "CORR_* includes refraction term")
        eh["HIERARCH DOPPLER"]  = (bool(self.info.doppler_tracking_applied), "Online Doppler tracking applied")

        self._write_spectral_keywords(eh)

        # TDIM: keep the convention for fixed-length vectors (skip for VLA)
        if not is_var:
            for i in range(1, int(eh.get("TFIELDS", 0)) + 1):
                ttype = (eh.get(f"TTYPE{i}", "") or "").strip().upper()
                tform = (eh.get(f"TFORM{i}", "") or "").strip().upper()
                if ttype in ("DATA", "FLAG", "FREQ") and (not tform.startswith(("P", "Q"))):
                    eh[f"TDIM{i}"] = (f"({int(nchan_hint)})", "Vector length per row")

        for k, v in self.info.shared_meta.items():
                set_meta_keyword(eh, f"HIERARCH {k}", v)

        return hdu

    def write(self, filename: str, overwrite: bool = True) -> None:
        b = self._buf
        if len(b) == 0:
            raise RuntimeError("No data to write.")
        mjd0 = float(builtins.min(b.MJD))
        mjd1 = float(builtins.max(b.MJD))
        pri = self._build_primary_hdu(mjd0, mjd1)
        tab = self._build_table_hdu()
        hdus = [pri, tab]
        hist_hdu = build_history_hdu(self._history)
        if hist_hdu is not None:
            hdus.append(hist_hdu)
        fits.HDUList(hdus).writeto(filename, overwrite=overwrite)
        print(f"[{SWNAME} {__version__}] Saved FITS: {filename} (rows={len(b)})")

    def _flush_part(self) -> None:
        if self.chunk_size is None:
            return
        b = self._buf
        if len(b) == 0:
            return
        self._part_index += 1
        partname = f"{self.out_basename}_part{self._part_index:04d}.fits"

        mjd0 = float(builtins.min(b.MJD))
        mjd1 = float(builtins.max(b.MJD))
        pri = self._build_primary_hdu(mjd0, mjd1)
        tab = self._build_table_hdu()
        hdus = [pri, tab]
        hist_hdu = build_history_hdu(self._history)
        if hist_hdu is not None:
            hdus.append(hist_hdu)
        fits.HDUList(hdus).writeto(partname, overwrite=True)

        self._manifest["parts"].append({
            "file": partname,
            "rows": len(b),
            "mjdstart": mjd0,
            "mjdend": mjd1,
        })
        print(f"[{SWNAME} {__version__}] Flushed part: {partname} (rows={len(b)})")
        b.clear()

    def close(self) -> None:
        if self.chunk_size is None:
            return
        self._flush_part()
        manifest_name = f"{self.out_basename}_manifest.json"
        with open(manifest_name, "w", encoding="utf-8") as f:
            json.dump(self._manifest, f, ensure_ascii=False, indent=2)
        print(f"[{SWNAME} {__version__}] Wrote manifest: {manifest_name}")

    @staticmethod
    def airmass_from_el_deg(el_deg: float) -> float:
        el = math.radians(el_deg)
        s = math.sin(el)
        if s <= 0:
            return float("nan")
        return 1.0 / s
