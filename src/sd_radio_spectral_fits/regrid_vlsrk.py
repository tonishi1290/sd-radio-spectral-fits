# src/sd_radio_spectral_fits/regrid_vlsrk.py
from __future__ import annotations

import builtins
from dataclasses import dataclass
from typing import Optional, Tuple, Union, List, Dict, Any

import numpy as np
import pandas as pd

from .axis import freq_axis_from_wcs, radio_velocity_kms
from .doppler import get_doppler_factor, C_KMS

# =========================================================
# 1. データ構造 & ヘルパー
# =========================================================

@dataclass
class VGrid:
    """共通速度グリッド定義"""
    v0_kms: float
    dv_kms: float
    nchan: int
    crpix1: float = 1.0  # 1-based

    def axis(self) -> np.ndarray:
        i = np.arange(self.nchan, dtype=float)
        return self.v0_kms + ((i + 1.0) - self.crpix1) * self.dv_kms

    def to_dict(self) -> dict:
        return {
            "v0_kms": self.v0_kms,
            "dv_kms": self.dv_kms,
            "nchan": self.nchan,
            "crpix1": self.crpix1
        }


def make_vgrid(vmin_kms: float, vmax_kms: float, dv_kms: float) -> VGrid:
    """範囲と分解能からVGridを作成"""
    v1, v2 = float(vmin_kms), float(vmax_kms)
    if v2 < v1:
        v1, v2 = v2, v1
    dv = abs(float(dv_kms))
    if dv <= 0:
        raise ValueError("dv_kms must be positive")
    
    if dv > (v2 - v1):
        # 1 channel covering the center
        return VGrid(v0_kms=(v1+v2)/2.0, dv_kms=dv, nchan=1, crpix1=1.0)

    n = int(np.floor((v2 - v1) / dv)) + 1
    return VGrid(v0_kms=v1, dv_kms=dv, nchan=n, crpix1=1.0)


def _get_val(row: pd.Series, meta: dict, key: str, default: Any = None) -> Any:
    """行データを優先し、無ければメタデータを参照する。"""
    # 1) Row override: plain <KEY>
    if key in row.index:
        val = row[key]
        if val is not None and not pd.isna(val):
            return val

    # 2) Meta fallback
    if key in meta:
        val = meta.get(key)
        if val is not None and not pd.isna(val):
            return val
    return default


def _get_specsys(row: pd.Series, meta: dict) -> str:
    for key in ("SPECSYS", "SSYSOBS"):
        val = _get_val(row, meta, key, None)
        if val not in (None, ""):
            return str(val).strip().upper()
    return ""


def _vcorr_scale_to_kms(key: str) -> float:
    ku = str(key or "").strip().upper()
    if ku.endswith("_KMS") or ku == "V_CORR_KMS":
        return 1.0
    # Standard policy: VELOSYS / VFRAME are stored in m/s.
    return 1.0e-3


def _get_vcorr_kms(row: pd.Series, meta: dict, preferred_key: str = "VELOSYS") -> float:
    candidates = ["VELOSYS"]
    if preferred_key and preferred_key not in candidates:
        candidates.append(preferred_key)
    if "VFRAME" not in candidates:
        candidates.append("VFRAME")
    if "V_CORR_KMS" not in candidates:
        candidates.append("V_CORR_KMS")

    for key in candidates:
        val = _get_val(row, meta, key, None)
        if val is None or pd.isna(val):
            continue
        return float(val) * _vcorr_scale_to_kms(key)
    return 0.0


def get_axis_signature(row: pd.Series, meta: dict, nchan: int) -> tuple:
    """
    スペクトルの軸特性を一意に識別する「指紋」を生成する。

    (NCHAN, CRVAL1, CDELT1, CRPIX1, CTYPE1, RESTFREQ)

    注意:
    - 正本は通常の row/header WCS とする。
    - REST 周波数は RESTFREQ / RESTFRQ の両方を許容する。
    """
    restfreq = _get_val(row, meta, "RESTFREQ", None)
    if restfreq is None or (isinstance(restfreq, float) and not np.isfinite(restfreq)) or float(restfreq) <= 0:
        restfreq = _get_val(row, meta, "RESTFRQ", 0.0)
    try:
        restfreq_f = float(restfreq)
    except Exception:
        restfreq_f = 0.0

    return (
        int(nchan),
        float(_get_val(row, meta, "CRVAL1", 0.0)),
        float(_get_val(row, meta, "CDELT1", 1.0)),
        float(_get_val(row, meta, "CRPIX1", 1.0)),
        str(_get_val(row, meta, "CTYPE1", "FREQ")).strip().upper(),
        str(_get_val(row, meta, "CUNIT1", "Hz")).strip().lower(), # 追加
        restfreq_f,
    )


# =========================================================
# 2. Standardizer (標準化エンジン)
# =========================================================

class Standardizer:
    """
    不均質なScantableを、共通の速度グリッドを持つ単一の行列に変換するクラス。
    AxisSignatureによるグルーピングを行い、WCS計算コストを最小化する。
    """
    # [MODIFIED] Added v_corr_col parameter
    def __init__(self, scantable, target_grid: Optional[VGrid] = None, v_corr_col: str = "VELOSYS"):
        self.st = scantable
        self.target_grid = target_grid
        self.v_corr_col = v_corr_col
        
        # キャッシュ: signature -> base_axis_array (Hz or km/s)
        self._base_axis_cache = {}
        
    def _get_base_axis(self, sig: tuple) -> np.ndarray:
        """
        Signatureに対応するベース軸を作成（キャッシュ対応）。
        """
        if sig in self._base_axis_cache:
            return self._base_axis_cache[sig]

        (nchan, crval, cdelt, crpix, ctype, cunit, restfreq) = sig
        
        pix = np.arange(nchan, dtype=float) + 1.0
        val_axis = crval + (pix - crpix) * cdelt
        
        obs_freq = val_axis # Default assumption: CTYPE=FREQ, CUNIT=Hz
        
        # CTYPEによる補正
        if "VEL" in ctype:
            # c_kms = 299792.458 # unused here, used in _generate_obs_freq_axis
            if restfreq > 0:
                dummy_meta = {
                    "CRVAL1": crval, "CDELT1": cdelt, "CRPIX1": crpix,
                    "CTYPE1": ctype, "RESTFREQ": restfreq, "CUNIT1": cunit
                }
                obs_freq = _generate_obs_freq_axis(dummy_meta, nchan)
                
        self._base_axis_cache[sig] = obs_freq
        return obs_freq

    def _group_by_signature(self) -> Dict[tuple, List[int]]:
        """Scantableの全行をAxisSignatureでグルーピングする"""
        df = self.st.table
        meta = self.st.meta
        is_list_data = isinstance(self.st.data, list)
        
        groups = {}
        for i in range(len(df)):
            if is_list_data:
                nc = len(self.st.data[i])
            else:
                nc = self.st.data.shape[1]
                
            sig = get_axis_signature(df.iloc[i], meta, nc)
            if sig not in groups:
                groups[sig] = []
            groups[sig].append(i)
        return groups

    def auto_determine_grid(
        self, 
        dv: Optional[float] = None, 
        vmin: Optional[float] = None, 
        vmax: Optional[float] = None, 
        margin: float = 0.0
    ) -> VGrid:
        """
        Scantable内の全データを包含するマスターグリッドを自動決定する。
        
        Parameters
        ----------
        dv : float, optional
            速度分解能 (km/s)。Noneの場合、全データ中の最小分解能を採用。
        vmin : float, optional
            速度範囲の下限 (km/s)。Noneの場合、全データの最小値 - margin。
        vmax : float, optional
            速度範囲の上限 (km/s)。Noneの場合、全データの最大値 + margin。
        margin : float
            自動範囲決定時に追加するマージン (km/s)。
        """
        # 全て指定されている場合は即作成
        if dv is not None and vmin is not None and vmax is not None:
            return make_vgrid(vmin, vmax, dv)

        groups = self._group_by_signature()
        
        global_min_dv = float("inf")
        global_vmin = float("inf")
        global_vmax = float("-inf")
        
        found_any = False
        
        # 全Signatureを走査して min_dv と min/max range を探す
        for sig, idxs in groups.items():
            f_obs = self._get_base_axis(sig)
            if len(f_obs) < 2: continue
            
            (_, _, _, _, _, _, restfreq) = sig
            if restfreq <= 0: continue
            
            # dv 推定 (km/s)
            df_obs = abs(f_obs[1] - f_obs[0])
            dv_est = C_KMS * df_obs / restfreq
            if dv_est < global_min_dv:
                global_min_dv = dv_est
            
            # Range 推定 (v_corr込み)
            # [MODIFIED] Check SPECSYS using strict TOPO rule
            row_ref = self.st.table.iloc[idxs[0]]
            specsys = _get_specsys(row_ref, self.st.meta)
            
            # 従来は "LSR" in specsys で判定していたが、
            # "TOPO" または "TOPOCENT" が含まれる場合のみ補正を行うように厳格化
            if "TOPO" in specsys:
                v_corrs = np.asarray([_get_vcorr_kms(self.st.table.iloc[j], self.st.meta, self.v_corr_col) for j in idxs], dtype=float)
                if len(v_corrs) == 0: v_c_min, v_c_max = 0.0, 0.0
                else: v_c_min, v_c_max = float(np.min(v_corrs)), float(np.max(v_corrs))
            else:
                # LSRK, HELIO, BARY等は補正済み(=0)とみなす
                v_c_min, v_c_max = 0.0, 0.0
            
            # v_topo (v_corr=0) range
            v_topo_axis = C_KMS * (1.0 - f_obs / restfreq)
            vt_min, vt_max = np.min(v_topo_axis), np.max(v_topo_axis)
            
            # v_lsrk range
            g_min = vt_min + v_c_min
            g_max = vt_max + v_c_max
            
            if g_min < global_vmin: global_vmin = g_min
            if g_max > global_vmax: global_vmax = g_max
            
            found_any = True
            
        if not found_any:
            # Fallback
            return make_vgrid(-100, 100, 0.1 if dv is None else dv)
            
        # 決定されたパラメータ (Noneなら自動値を使用)
        target_dv = dv if dv is not None else global_min_dv
        target_vmin = vmin if vmin is not None else (global_vmin - margin)
        target_vmax = vmax if vmax is not None else (global_vmax + margin)
        
        return make_vgrid(target_vmin, target_vmax, target_dv)

    def get_matrix(
        self, 
        target_grid: Optional[VGrid] = None,
        dv: Optional[float] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        全データをリグリッドし、(N_row, N_grid) の行列と、その速度軸を返す。
        
        Parameters
        ----------
        target_grid : VGrid, optional
            定義済みのグリッドオブジェクト。
        dv, vmin, vmax : float, optional
            target_grid が指定されていない場合に、自動決定ロジックへ渡すパラメータ。
            指定された値は優先的に使用される。
        """
        # 1. ターゲットグリッドの確定
        if target_grid is None:
            # ユーザー指定パラメータがある場合は、それを使って新規作成
            if dv is not None or vmin is not None or vmax is not None:
                target_grid = self.auto_determine_grid(dv=dv, vmin=vmin, vmax=vmax)
            # 指定がなく、キャッシュがあれば使う
            elif self.target_grid is not None:
                target_grid = self.target_grid
            # キャッシュもなければフルオートで作成
            else:
                self.target_grid = self.auto_determine_grid()
                target_grid = self.target_grid
            
        v_tgt = target_grid.axis()
        n_tgt = len(v_tgt)
        n_rows = len(self.st.table)
        
        # 2. 結果格納用行列 (NaNで初期化)
        out_matrix = np.full((n_rows, n_tgt), np.nan, dtype=np.float32)
        
        # 3. グルーピング
        groups = self._group_by_signature()
        
        # 4. グループごとに一括処理
        rest_freq_global = float(self.st.meta.get("RESTFRQ", self.st.meta.get("RESTFREQ", 0)))
        is_list_data = isinstance(self.st.data, list)
        
        for sig, idxs in groups.items():
            (_, _, _, _, _, _, sig_rest) = sig
            rest_freq = sig_rest if sig_rest > 0 else rest_freq_global
            if rest_freq <= 0: continue

            f_obs_axis = self._get_base_axis(sig)
            
            if is_list_data:
                raw_spectra = [self.st.data[i] for i in idxs]
                raw_block = np.vstack(raw_spectra)
            else:
                raw_block = self.st.data[idxs]
            
            # [MODIFIED] Check SPECSYS using strict TOPO rule
            row_ref = self.st.table.iloc[idxs[0]]
            specsys = _get_specsys(row_ref, self.st.meta)
            
            if "TOPO" in specsys:
                v_corrs = np.asarray([_get_vcorr_kms(self.st.table.iloc[j], self.st.meta, self.v_corr_col) for j in idxs], dtype=float)
                v_corrs = np.nan_to_num(v_corrs, nan=0.0)
            else:
                v_corrs = np.zeros(len(idxs), dtype=float)
                
            k_factors = get_doppler_factor(v_corrs)
            k_inv_arr = 1.0 / k_factors
            term_b = (C_KMS * f_obs_axis) / rest_freq


            if len(term_b) > 1:
                # 速度軸の向き（昇順か降順か）はグループ内で不変のため、最初に1度だけ判定する
                # np.interp は x が昇順である必要がある
                is_reversed = (C_KMS - term_b[0] * k_inv_arr[0]) > (C_KMS - term_b[-1] * k_inv_arr[0])
                
                if is_reversed:
                    term_b_rev = term_b[::-1]
                    for local_i, global_i in enumerate(idxs):
                        v_axis_row = C_KMS - term_b_rev * k_inv_arr[local_i]
                        out_matrix[global_i] = np.interp(v_tgt, v_axis_row, raw_block[local_i, ::-1], left=np.nan, right=np.nan)
                else:
                    for local_i, global_i in enumerate(idxs):
                        v_axis_row = C_KMS - term_b * k_inv_arr[local_i]
                        out_matrix[global_i] = np.interp(v_tgt, v_axis_row, raw_block[local_i], left=np.nan, right=np.nan)
            else:
                # チャンネル数が1以下の場合の安全なフォールバック
                for local_i, global_i in enumerate(idxs):
                    v_axis_row = C_KMS - term_b * k_inv_arr[local_i]
                    out_matrix[global_i] = interp_to_vgrid(v_axis_row, raw_block[local_i], v_tgt)
                
        return out_matrix, v_tgt


# =========================================================
# 3. 既存ヘルパー (Legacy Support)
# =========================================================

def _get_restfreq(meta: dict) -> float:
    for k in ("RESTFRQ", "RESTFREQ", "rest_hz"):
        if k in meta and meta[k] not in (None, ""):
            return float(meta[k])
    raise ValueError("RESTFREQ missing.")

def _normalize_unit_local(unit_str: str) -> str:
    u = str(unit_str).strip().lower()
    if "hz" in u: return "hz"
    if u in ("m/s", "m s-1", "ms-1", "meter/sec", "m/sec"): return "m/s"
    if u in ("km/s", "km s-1", "kms-1", "kilometer/sec", "km/sec"): return "km/s"
    return u


def _freq_scale_to_hz_local(unit_str: str) -> float:
    u = str(unit_str).strip().lower()
    if u in ("", "hz"):
        return 1.0
    if u == "khz":
        return 1.0e3
    if u == "mhz":
        return 1.0e6
    if u == "ghz":
        return 1.0e9
    if u == "thz":
        return 1.0e12
    return 1.0

def _generate_obs_freq_axis(meta: dict, nchan: int) -> np.ndarray:
    raw_ctype = str(meta.get("CTYPE1", "FREQ")).strip()
    ctype = raw_ctype.upper()
    raw_cunit = str(meta.get("CUNIT1", "")) 
    unit_norm = _normalize_unit_local(raw_cunit)

    if ctype.startswith("FREQ") or (unit_norm == "hz"):
        freq = freq_axis_from_wcs(meta, nchan=nchan)
        return freq * _freq_scale_to_hz_local(raw_cunit)

    if "CRVAL1" not in meta or "CDELT1" not in meta:
        return np.arange(nchan, dtype=float)

    crval = float(meta["CRVAL1"])
    cdelt = float(meta["CDELT1"])
    crpix = float(meta.get("CRPIX1", 1.0))
    scale_to_kms = 0.001 if unit_norm == "m/s" else 1.0
    
    i = np.arange(nchan, dtype=float)
    v_axis_kms = (crval + (i + 1.0 - crpix) * cdelt) * scale_to_kms
    
    try:
        rf = _get_restfreq(meta)
    except ValueError:
        return v_axis_kms 

    if ctype.startswith("VRAD"):
        return rf * (1.0 - v_axis_kms / C_KMS)
    elif ctype.startswith("VOPT"):
        return rf / (v_axis_kms / C_KMS + 1.0)
    elif ctype.startswith("VELO"):
        beta = np.clip(v_axis_kms / C_KMS, -0.999, 0.999)
        return rf * np.sqrt((1.0 - beta) / (1.0 + beta))
    else:
        return rf * (1.0 - v_axis_kms / C_KMS)

def vlsrk_axis_for_spectrum(meta: dict, *, v_corr_kms: float, nchan: Optional[int] = None) -> np.ndarray:
    if nchan is None: nchan = int(meta.get("NAXIS1", 0))
    f_obs = _generate_obs_freq_axis(meta, nchan)
    k = get_doppler_factor(v_corr_kms)
    f_lsrk = f_obs / k
    try:
        rf = _get_restfreq(meta)
        return radio_velocity_kms(f_lsrk, rf)
    except ValueError:
        return f_lsrk

def interp_to_vgrid(v_src: np.ndarray, y_src: np.ndarray, v_tgt: np.ndarray) -> np.ndarray:
    v_src = np.asarray(v_src, float)
    y_src = np.asarray(y_src, float)
    v_tgt = np.asarray(v_tgt, float)
    if v_src.size == 0: return np.full_like(v_tgt, np.nan)
    if v_src[0] > v_src[-1]: v, y = v_src[::-1], y_src[::-1]
    else: v, y = v_src, y_src
    m = np.isfinite(y)
    if np.count_nonzero(m) < 2: return np.full_like(v_tgt, np.nan)
    return np.interp(v_tgt, v[m], y[m], left=np.nan, right=np.nan)

def vrange_from_meta_and_vcorr(meta: dict, v_corr_kms: np.ndarray, *, nchan: int | None = None) -> tuple[float, float]:
    n = int(nchan) if nchan is not None else int(meta.get("NAXIS1", 1))
    f_obs = _generate_obs_freq_axis(meta, n)
    f_min, f_max = np.min(f_obs), np.max(f_obs)
    vc_min = float(np.nanmin(v_corr_kms)) if len(v_corr_kms) > 0 else 0.0
    vc_max = float(np.nanmax(v_corr_kms)) if len(v_corr_kms) > 0 else 0.0
    k_min = get_doppler_factor(vc_min)
    k_max = get_doppler_factor(vc_max)
    f_l_min = f_min / k_max
    f_l_max = f_max / k_min
    try:
        rf = _get_restfreq(meta)
        v1 = C_KMS * (1.0 - f_l_min / rf)
        v2 = C_KMS * (1.0 - f_l_max / rf)
        return (float(min(v1, v2)), float(max(v1, v2)))
    except ValueError:
        return (0.0, 1.0)
