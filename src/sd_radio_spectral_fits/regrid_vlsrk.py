from __future__ import annotations

import builtins
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any

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
            "crpix1": self.crpix1,
        }


def make_vgrid(vmin_kms: float, vmax_kms: float, dv_kms: float) -> VGrid:
    """範囲と分解能からVGridを作成する legacy helper。"""
    v1, v2 = float(vmin_kms), float(vmax_kms)
    if v2 < v1:
        v1, v2 = v2, v1
    dv = abs(float(dv_kms))
    if dv <= 0:
        raise ValueError("dv_kms must be positive")

    if dv > (v2 - v1):
        return VGrid(v0_kms=(v1 + v2) / 2.0, dv_kms=dv, nchan=1, crpix1=1.0)

    n = int(np.floor((v2 - v1) / dv)) + 1
    return VGrid(v0_kms=v1, dv_kms=dv, nchan=n, crpix1=1.0)


def make_vgrid_from_center(center_kms: float, dv_kms: float, nchan: int) -> VGrid:
    """中心速度・分解能・チャンネル数から VGrid を作成する。"""
    dv = abs(float(dv_kms))
    n = int(nchan)
    if dv <= 0:
        raise ValueError("dv_kms must be positive")
    if n < 1:
        raise ValueError("nchan must be >= 1")
    v0 = float(center_kms) - 0.5 * (n - 1) * dv
    return VGrid(v0_kms=v0, dv_kms=dv, nchan=n, crpix1=1.0)


def make_vgrid_anchored(
    vmin_kms: float,
    vmax_kms: float,
    dv_kms: float,
    *,
    anchor_kms: float = 0.0,
    tol_frac: float = 1.0e-9,
) -> VGrid:
    """
    anchor 付きの決定論的グリッド作成。

    floor/ceil の境界で極小の丸め差が出ても結果がぶれにくいよう、
    anchor 基準のインデックスへ量子化する。
    """
    v1, v2 = float(vmin_kms), float(vmax_kms)
    if v2 < v1:
        v1, v2 = v2, v1
    dv = abs(float(dv_kms))
    if dv <= 0:
        raise ValueError("dv_kms must be positive")

    if dv > (v2 - v1):
        return VGrid(v0_kms=(v1 + v2) / 2.0, dv_kms=dv, nchan=1, crpix1=1.0)

    tol = abs(float(tol_frac))
    start_idx = int(np.floor(((v1 - anchor_kms) / dv) + tol))
    end_idx = int(np.ceil(((v2 - anchor_kms) / dv) - tol))
    n = builtins.max(1, end_idx - start_idx + 1)
    v0 = float(anchor_kms) + start_idx * dv
    return VGrid(v0_kms=v0, dv_kms=dv, nchan=n, crpix1=1.0)


def _quantize_to_anchor(value_kms: float, dv_kms: float, anchor_kms: float) -> float:
    dv = abs(float(dv_kms))
    if dv <= 0:
        raise ValueError("dv_kms must be positive")
    idx = int(np.round((float(value_kms) - float(anchor_kms)) / dv))
    return float(anchor_kms) + idx * dv


def _cover_nchan_from_center(
    center_kms: float,
    vmin_kms: float,
    vmax_kms: float,
    dv_kms: float,
    *,
    tol_frac: float = 1.0e-9,
) -> int:
    """center を固定したまま [vmin, vmax] を包含する最小 nchan を決定論的に求める。"""
    dv = abs(float(dv_kms))
    if dv <= 0:
        raise ValueError("dv_kms must be positive")
    v1, v2 = float(vmin_kms), float(vmax_kms)
    if v2 < v1:
        v1, v2 = v2, v1
    tol = abs(float(tol_frac))
    half_left_idx = int(np.ceil(((float(center_kms) - v1) / dv) - tol))
    half_right_idx = int(np.ceil(((v2 - float(center_kms)) / dv) - tol))
    half_idx = builtins.max(half_left_idx, half_right_idx, 0)
    return 2 * half_idx + 1


def _get_val(row: pd.Series, meta: dict, key: str, default: Any = None) -> Any:
    """行データを優先し、無ければメタデータを参照する。"""
    if key in row.index:
        val = row[key]
        if val is not None and not pd.isna(val):
            return val
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
    return 1.0e-3


def _get_vcorr_kms(row: pd.Series, meta: dict, preferred_key: str = "VFRAME") -> float:
    preferred = str(preferred_key or "VFRAME").strip().upper()
    candidates: List[str] = []
    for key in (preferred, "VFRAME", "V_CORR_KMS", "VELOSYS"):
        ku = str(key).strip().upper()
        if ku and ku not in candidates:
            candidates.append(ku)
    for key in candidates:
        val = _get_val(row, meta, key, None)
        if val is None or pd.isna(val):
            continue
        return float(val) * _vcorr_scale_to_kms(key)
    return 0.0


def get_axis_signature(row: pd.Series, meta: dict, nchan: int) -> tuple:
    """スペクトル軸の識別 signature。"""
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
    不均質な Scantable を共通速度グリッドへ変換するクラス。

    Parameters
    ----------
    auto_grid_mode : {"centered_native", "legacy_envelope"}
        centered_native
            各 row の native スペクトル幅を尊重しつつ、中心速度を anchor 基準へ量子化した
            決定論的なグリッドを自動生成する。環境差に強い推奨モード。
        legacy_envelope
            従来どおり全 row の速度範囲 envelope を用いて自動生成する。
    auto_grid_anchor_kms : float
        centered_native および anchored envelope で使う基準点。
    auto_grid_tol_frac : float
        legacy_envelope + anchored grid の境界丸め用微小許容値。
    """

    def __init__(
        self,
        scantable,
        target_grid: Optional[VGrid] = None,
        v_corr_col: str = "VFRAME",
        *,
        auto_grid_mode: str = "centered_native",
        auto_grid_anchor_kms: float = 0.0,
        auto_grid_tol_frac: float = 1.0e-9,
    ):
        self.st = scantable
        self.target_grid = target_grid
        self.v_corr_col = str(v_corr_col or "VFRAME").strip().upper() or "VFRAME"
        self.auto_grid_mode = str(auto_grid_mode or "centered_native").strip().lower()
        self.auto_grid_anchor_kms = float(auto_grid_anchor_kms)
        self.auto_grid_tol_frac = float(auto_grid_tol_frac)
        self._base_axis_cache: Dict[tuple, np.ndarray] = {}

    def _get_base_axis(self, sig: tuple) -> np.ndarray:
        if sig in self._base_axis_cache:
            return self._base_axis_cache[sig]

        (nchan, crval, cdelt, crpix, ctype, cunit, restfreq, _specsys) = sig
        pix = np.arange(nchan, dtype=float) + 1.0
        val_axis = crval + (pix - crpix) * cdelt
        obs_freq = val_axis

        if "VEL" in ctype and restfreq > 0:
            dummy_meta = {
                "CRVAL1": crval,
                "CDELT1": cdelt,
                "CRPIX1": crpix,
                "CTYPE1": ctype,
                "RESTFREQ": restfreq,
                "CUNIT1": cunit,
            }
            obs_freq = _generate_obs_freq_axis(dummy_meta, nchan)

        self._base_axis_cache[sig] = obs_freq
        return obs_freq

    def _group_by_signature(self) -> Dict[tuple, List[int]]:
        df = self.st.table
        meta = self.st.meta
        is_list_data = isinstance(self.st.data, list)
        groups: Dict[tuple, List[int]] = {}
        for i in range(len(df)):
            nc = len(self.st.data[i]) if is_list_data else self.st.data.shape[1]
            sig = get_axis_signature(df.iloc[i], meta, nc)
            groups.setdefault(sig, []).append(i)
        return groups

    def _row_axis_bounds_and_center(
        self,
        f_obs_axis: np.ndarray,
        restfreq: float,
        v_corr_kms: float,
    ) -> Tuple[float, float, float]:
        """単一 row の LSRK 速度軸 min/max/center を返す。"""
        if f_obs_axis.size == 0:
            return 0.0, 0.0, 0.0
        k_inv = 1.0 / get_doppler_factor(float(v_corr_kms))
        term_b0 = (C_KMS * float(f_obs_axis[0])) / restfreq
        term_b1 = (C_KMS * float(f_obs_axis[-1])) / restfreq
        v0 = C_KMS - term_b0 * k_inv
        v1 = C_KMS - term_b1 * k_inv
        vmin = builtins.min(v0, v1)
        vmax = builtins.max(v0, v1)
        center = 0.5 * (vmin + vmax)
        return float(vmin), float(vmax), float(center)

    def auto_determine_grid(
        self,
        dv: Optional[float] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        margin: float = 0.0,
    ) -> VGrid:
        if dv is not None and vmin is not None and vmax is not None:
            # 明示指定はそのまま尊重する。
            return make_vgrid(vmin, vmax, dv)

        groups = self._group_by_signature()
        global_min_dv = float("inf")
        envelope_vmin = float("inf")
        envelope_vmax = float("-inf")
        native_required_nchan = 0
        centers: List[float] = []
        found_any = False

        for sig, idxs in groups.items():
            (nchan_sig, _crval, _cdelt, _crpix, _ctype, _cunit, restfreq, specsys) = sig
            f_obs = self._get_base_axis(sig)
            if len(f_obs) < 2 or restfreq <= 0:
                continue

            df_obs = abs(float(f_obs[1]) - float(f_obs[0]))
            dv_est = C_KMS * df_obs / restfreq
            if np.isfinite(dv_est) and dv_est > 0 and dv_est < global_min_dv:
                global_min_dv = float(dv_est)

            specsys_u = str(specsys).upper()
            for j in idxs:
                row = self.st.table.iloc[j]
                v_corr = _get_vcorr_kms(row, self.st.meta, self.v_corr_col) if "TOPO" in specsys_u else 0.0
                row_vmin, row_vmax, row_center = self._row_axis_bounds_and_center(f_obs, restfreq, v_corr)
                envelope_vmin = builtins.min(envelope_vmin, row_vmin)
                envelope_vmax = builtins.max(envelope_vmax, row_vmax)
                centers.append(row_center)
                found_any = True

            target_dv_tmp = abs(float(dv)) if dv is not None else abs(float(dv_est))
            if target_dv_tmp > 0:
                native_span = abs(float(nchan_sig - 1)) * abs(float(dv_est))
                n_req = int(np.ceil(native_span / target_dv_tmp)) + 1
                native_required_nchan = builtins.max(native_required_nchan, n_req)

        if not found_any:
            target_dv = 0.1 if dv is None else float(dv)
            return make_vgrid_anchored(
                -100.0,
                100.0,
                target_dv,
                anchor_kms=self.auto_grid_anchor_kms,
                tol_frac=self.auto_grid_tol_frac,
            )

        target_dv = abs(float(dv)) if dv is not None else float(global_min_dv)
        if not np.isfinite(target_dv) or target_dv <= 0:
            target_dv = 0.1

        # 明示的に範囲指定された場合は、その拘束を最優先する。
        if vmin is not None or vmax is not None:
            vmin_eff = float(vmin) if vmin is not None else (envelope_vmin - float(margin))
            vmax_eff = float(vmax) if vmax is not None else (envelope_vmax + float(margin))
            return make_vgrid(vmin_eff, vmax_eff, target_dv)

        if self.auto_grid_mode == "legacy_envelope":
            # 従来挙動をなるべくそのまま保つ。
            return make_vgrid(
                envelope_vmin - float(margin),
                envelope_vmax + float(margin),
                target_dv,
            )

        # 推奨: centered_native
        if len(centers) == 0:
            center_ref = 0.5 * (envelope_vmin + envelope_vmax)
        else:
            center_ref = float(np.median(np.asarray(centers, dtype=float)))
        center_q = _quantize_to_anchor(center_ref, target_dv, self.auto_grid_anchor_kms)

        # native スペクトル幅だけでなく、row ごとの中心ずれも必ず包含する。
        native_cover_nchan = _cover_nchan_from_center(
            center_q,
            envelope_vmin,
            envelope_vmax,
            target_dv,
            tol_frac=self.auto_grid_tol_frac,
        )

        if native_required_nchan < 1:
            native_required_nchan = native_cover_nchan
        native_required_nchan = builtins.max(int(native_required_nchan), int(native_cover_nchan), 1)

        # 念のため margin をチャンネル数へ反映する（既定 0 なら互換的）。
        if margin > 0:
            extra = int(np.ceil(float(margin) / target_dv))
            native_required_nchan += 2 * extra

        return make_vgrid_from_center(center_q, target_dv, native_required_nchan)

    def get_matrix(
        self,
        target_grid: Optional[VGrid] = None,
        dv: Optional[float] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        if target_grid is None:
            if dv is not None or vmin is not None or vmax is not None:
                target_grid = self.auto_determine_grid(dv=dv, vmin=vmin, vmax=vmax)
            elif self.target_grid is not None:
                target_grid = self.target_grid
            else:
                self.target_grid = self.auto_determine_grid()
                target_grid = self.target_grid

        v_tgt = target_grid.axis()
        n_tgt = len(v_tgt)
        n_rows = len(self.st.table)
        out_matrix = np.full((n_rows, n_tgt), np.nan, dtype=np.float32)

        groups = self._group_by_signature()
        rest_freq_global = float(self.st.meta.get("RESTFRQ", self.st.meta.get("RESTFREQ", 0)))
        is_list_data = isinstance(self.st.data, list)

        for sig, idxs in groups.items():
            (_nchan, _crval, _cdelt, _crpix, _ctype, _cunit, sig_rest, specsys) = sig
            rest_freq = sig_rest if sig_rest > 0 else rest_freq_global
            if rest_freq <= 0:
                continue

            f_obs_axis = self._get_base_axis(sig)
            raw_block = np.vstack([self.st.data[i] for i in idxs]) if is_list_data else self.st.data[idxs]

            if "TOPO" in str(specsys).upper():
                v_corrs = np.asarray(
                    [_get_vcorr_kms(self.st.table.iloc[j], self.st.meta, self.v_corr_col) for j in idxs],
                    dtype=float,
                )
                v_corrs = np.nan_to_num(v_corrs, nan=0.0)
            else:
                v_corrs = np.zeros(len(idxs), dtype=float)

            k_inv_arr = 1.0 / get_doppler_factor(v_corrs)
            term_b = (C_KMS * f_obs_axis) / rest_freq

            if len(term_b) > 1:
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
    if "hz" in u:
        return "hz"
    if u in ("m/s", "m s-1", "ms-1", "meter/sec", "m/sec"):
        return "m/s"
    if u in ("km/s", "km s-1", "kms-1", "kilometer/sec", "km/sec"):
        return "km/s"
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
    if nchan is None:
        nchan = int(meta.get("NAXIS1", 0))
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
    if v_src.size == 0:
        return np.full_like(v_tgt, np.nan)
    if v_src[0] > v_src[-1]:
        v, y = v_src[::-1], y_src[::-1]
    else:
        v, y = v_src, y_src
    m = np.isfinite(y)
    if np.count_nonzero(m) < 2:
        return np.full_like(v_tgt, np.nan)
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
