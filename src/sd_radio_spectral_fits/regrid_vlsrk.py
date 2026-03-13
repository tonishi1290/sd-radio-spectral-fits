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


def _get_vcorr_kms(row: pd.Series, meta: dict, preferred_key: str = "VFRAME") -> float:
    preferred = str(preferred_key or "").strip().upper()
    candidates: List[str] = []
    for key in (preferred, "VFRAME", "V_CORR_KMS", "VELOSYS"):
        if key and key not in candidates:
            candidates.append(key)

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
        str(_get_val(row, meta, "CUNIT1", "Hz")).strip().lower(),
        restfreq_f,
        _get_specsys(row, meta),
    )


# =========================================================
# 2. Standardizer (標準化エンジン)
# =========================================================

class Standardizer:
    """
    不均質なScantableを、共通の速度グリッドを持つ単一の行列に変換するクラス。
    AxisSignatureによるグルーピングを行い、WCS計算コストを最小化する。
    """
    def __init__(
        self,
        scantable,
        target_grid: Optional[VGrid] = None,
        v_corr_col: str = "VFRAME",
        auto_grid_mode: str = "stable_native",
        auto_grid_anchor_kms: float = 0.0,
        auto_grid_tol_frac: float = 1.0e-6,
    ):
        self.st = scantable
        self.target_grid = target_grid
        self.v_corr_col = str(v_corr_col).strip().upper() if v_corr_col is not None else "VFRAME"
        self.auto_grid_mode = str(auto_grid_mode).strip().lower()
        self.auto_grid_anchor_kms = float(auto_grid_anchor_kms)
        self.auto_grid_tol_frac = float(auto_grid_tol_frac)
        # キャッシュ: signature -> base_axis_array (Hz or km/s)
        self._base_axis_cache = {}
        
    def _get_base_axis(self, sig: tuple) -> np.ndarray:
        """
        Signatureに対応するベース軸を作成（キャッシュ対応）。
        """
        if sig in self._base_axis_cache:
            return self._base_axis_cache[sig]

        (nchan, crval, cdelt, crpix, ctype, cunit, restfreq, _specsys) = sig
        
        pix = np.arange(nchan, dtype=float) + 1.0
        val_axis = crval + (pix - crpix) * cdelt

        # CTYPE/CUNIT を尊重して観測周波数軸(Hz)へ正規化する。
        # これまでは CTYPE=FREQ なら val_axis をそのまま使っており、
        # CUNIT1 が kHz/MHz/GHz などでも暗黙に Hz 扱いされていた。
        # 基礎部分なので、FREQ 系も _generate_obs_freq_axis() に統一する。
        obs_freq = val_axis
        unit_norm = _normalize_unit_local(cunit)
        if ctype.startswith("FREQ") or unit_norm in ("hz",):
            dummy_meta = {
                "CRVAL1": crval, "CDELT1": cdelt, "CRPIX1": crpix,
                "CTYPE1": ctype, "RESTFREQ": restfreq, "CUNIT1": cunit
            }
            obs_freq = _generate_obs_freq_axis(dummy_meta, nchan)
        elif "VEL" in ctype and restfreq > 0:
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

    def _quantize_to_anchor(self, value_kms: float, dv_kms: float) -> float:
        dv = abs(float(dv_kms))
        if dv <= 0:
            return float(value_kms)
        anchor = float(self.auto_grid_anchor_kms)
        return anchor + np.round((float(value_kms) - anchor) / dv) * dv

    def _compute_axis_stats(self, sig: tuple, idxs: List[int]) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        f_obs = self._get_base_axis(sig)
        if len(f_obs) < 2:
            raise ValueError("nchan must be >= 2 for auto_determine_grid")
        (_, _, _, _, _, _, restfreq, specsys) = sig
        if restfreq <= 0:
            raise ValueError("REST frequency must be positive")

        df_obs = abs(f_obs[1] - f_obs[0])
        dv_est = C_KMS * df_obs / restfreq

        if "TOPO" in str(specsys).upper():
            v_corrs = np.asarray([
                _get_vcorr_kms(self.st.table.iloc[j], self.st.meta, self.v_corr_col) for j in idxs
            ], dtype=float)
            v_corrs = np.nan_to_num(v_corrs, nan=0.0)
        else:
            v_corrs = np.zeros(len(idxs), dtype=float)

        term_b = (C_KMS * f_obs) / restfreq
        if len(term_b) == 0:
            raise ValueError("empty spectral axis")

        v_first = C_KMS - term_b[0] / get_doppler_factor(v_corrs)
        v_last = C_KMS - term_b[-1] / get_doppler_factor(v_corrs)
        v_lo = np.minimum(v_first, v_last)
        v_hi = np.maximum(v_first, v_last)
        centers = 0.5 * (v_lo + v_hi)
        half_widths = 0.5 * (v_hi - v_lo)
        return float(dv_est), centers.astype(float), half_widths.astype(float), np.array([len(f_obs)], dtype=int)

    def _stable_native_grid(self, dv_kms: float, centers: np.ndarray, half_widths: np.ndarray, native_nchan: int, margin_kms: float = 0.0) -> VGrid:
        dv = abs(float(dv_kms))
        if dv <= 0:
            raise ValueError("dv_kms must be positive")
        centers = np.asarray(centers, dtype=float)
        half_widths = np.asarray(half_widths, dtype=float)
        if centers.size == 0 or half_widths.size == 0:
            return make_vgrid(-100.0, 100.0, dv)

        c_ref = self._quantize_to_anchor(float(np.median(centers)), dv)
        tol = builtins.max(1.0e-9, abs(dv) * builtins.max(self.auto_grid_tol_frac, 0.0))

        if native_nchan < 1:
            native_nchan = 1
        # Preserve native width when the input family is homogeneous.
        n_half = 0.5 * (native_nchan - 1) * dv
        v0 = c_ref - n_half

        # Expand only if needed to guarantee coverage of all rows.
        need_left = 0.0
        need_right = 0.0
        current_left = c_ref - n_half
        current_right = c_ref + n_half
        if centers.size > 0:
            row_left = float(np.min(centers - half_widths) - float(margin_kms))
            row_right = float(np.max(centers + half_widths) + float(margin_kms))
            need_left = builtins.max(0.0, current_left - row_left)
            need_right = builtins.max(0.0, row_right - current_right)

        extra_left = int(np.ceil((need_left - tol) / dv)) if need_left > tol else 0
        extra_right = int(np.ceil((need_right - tol) / dv)) if need_right > tol else 0
        nchan = int(native_nchan + extra_left + extra_right)
        v0 = float(v0 - extra_left * dv)
        return VGrid(v0_kms=v0, dv_kms=dv, nchan=nchan, crpix1=1.0)

    def auto_determine_grid(
        self,
        dv: Optional[float] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        margin: float = 0.0,
    ) -> VGrid:
        """
        Scantable内の全データを包含するマスターグリッドを自動決定する。

        既定では、同質な入力（同一 nchan / ほぼ同一 dv）の場合に、
        native のチャネル幅・幅をできるだけ保ちながら、代表中心を
        anchor 基準で量子化した安定なグリッドを返す。
        比較実験で端点の微差が 1--2 ch へ増幅されるのを抑えるため、
        legacy の global min/max envelope は明示モード時のみ使う。
        """
        if dv is not None and vmin is not None and vmax is not None:
            return make_vgrid(vmin, vmax, dv)

        groups = self._group_by_signature()
        if len(groups) == 0:
            return make_vgrid(-100.0, 100.0, 0.1 if dv is None else dv)

        dv_candidates: List[float] = []
        all_centers: List[np.ndarray] = []
        all_half_widths: List[np.ndarray] = []
        all_native_nchan: List[int] = []
        global_vmin = float("inf")
        global_vmax = float("-inf")

        for sig, idxs in groups.items():
            try:
                dv_est, centers, half_widths, nchan_arr = self._compute_axis_stats(sig, idxs)
            except Exception:
                continue
            dv_candidates.append(float(dv_est))
            all_centers.append(centers)
            all_half_widths.append(half_widths)
            all_native_nchan.extend([int(x) for x in nchan_arr])
            global_vmin = builtins.min(global_vmin, float(np.min(centers - half_widths)))
            global_vmax = builtins.max(global_vmax, float(np.max(centers + half_widths)))

        if len(dv_candidates) == 0 or len(all_centers) == 0:
            return make_vgrid(-100.0, 100.0, 0.1 if dv is None else dv)

        target_dv = abs(float(dv)) if dv is not None else float(np.min(dv_candidates))
        centers_all = np.concatenate(all_centers).astype(float)
        half_widths_all = np.concatenate(all_half_widths).astype(float)

        if vmin is not None or vmax is not None:
            target_vmin = float(vmin) if vmin is not None else float(global_vmin - margin)
            target_vmax = float(vmax) if vmax is not None else float(global_vmax + margin)
            return make_vgrid(target_vmin, target_vmax, target_dv)

        mode = self.auto_grid_mode
        if mode == "legacy_envelope":
            return make_vgrid(float(global_vmin - margin), float(global_vmax + margin), target_dv)

        homogeneous_nchan = (len(set(all_native_nchan)) == 1)
        homogeneous_dv = bool(np.allclose(np.asarray(dv_candidates, dtype=float), target_dv, rtol=1.0e-6, atol=builtins.max(1.0e-9, abs(target_dv) * builtins.max(self.auto_grid_tol_frac, 0.0))))
        if mode == "stable_native" and homogeneous_nchan and homogeneous_dv:
            native_nchan = int(all_native_nchan[0])
            return self._stable_native_grid(target_dv, centers_all, half_widths_all, native_nchan, margin_kms=margin)

        # Fallback: anchor-quantized envelope. More stable than raw global min/max.
        left = float(np.min(centers_all - half_widths_all) - margin)
        right = float(np.max(centers_all + half_widths_all) + margin)
        left_q = self._quantize_to_anchor(left, target_dv)
        right_q = self._quantize_to_anchor(right, target_dv)
        if right_q < right:
            right_q += target_dv
        if left_q > left:
            left_q -= target_dv
        return make_vgrid(left_q, right_q, target_dv)

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
            (_, _, _, _, _, _, sig_rest, sig_specsys) = sig
            rest_freq = sig_rest if sig_rest > 0 else rest_freq_global
            if rest_freq <= 0: continue

            f_obs_axis = self._get_base_axis(sig)
            
            if is_list_data:
                raw_spectra = [self.st.data[i] for i in idxs]
                raw_block = np.vstack(raw_spectra)
            else:
                raw_block = self.st.data[idxs]
            
            specsys = str(sig_specsys).upper()

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
