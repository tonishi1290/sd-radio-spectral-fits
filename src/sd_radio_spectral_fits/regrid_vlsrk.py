# src/sd_radio_spectral_fits/regrid_vlsrk.py
from __future__ import annotations

import builtins
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional, Tuple, Union, List, Dict, Any

import numpy as np
import pandas as pd

from .axis import freq_axis_from_wcs, radio_velocity_kms
from .doppler import get_doppler_factor, C_KMS
from .fitsio import Scantable, stamp_scantable_code_provenance

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


def _get_ssysobs(row: pd.Series, meta: dict) -> str:
    val = _get_val(row, meta, "SSYSOBS", None)
    if val not in (None, ""):
        return str(val).strip().upper()

    val = _get_val(row, meta, "SPECSYS", None)
    if val not in (None, ""):
        return str(val).strip().upper()

    return ""


def _vcorr_scale_to_kms(key: str) -> float:
    ku = str(key or "").strip().upper()
    if ku.endswith("_KMS") or ku == "V_CORR_KMS":
        return 1.0
    # Standard policy: VELOSYS / VFRAME are stored in m/s.
    return 1.0e-3


def _get_vcorr_kms(row: pd.Series, meta: dict, preferred_key: str = "VFRAME", *, row_only: bool = False) -> Optional[float]:
    """Return velocity correction in km/s.

    When ``row_only=True``, only table-row columns are consulted. This is the
    strict policy for TOPO row-wise corrections and intentionally ignores
    meta/header fallback values.
    """
    preferred = str(preferred_key or "").strip().upper()
    candidates: List[str] = []
    for key in (preferred, "VFRAME", "V_CORR_KMS", "VELOSYS"):
        if key and key not in candidates:
            candidates.append(key)

    for key in candidates:
        if row_only:
            if key not in row.index:
                continue
            val = row[key]
        else:
            val = _get_val(row, meta, key, None)
        if val is None or pd.isna(val):
            continue
        return float(val) * _vcorr_scale_to_kms(key)
    return None


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
        # 行単位の前計算キャッシュ（1回の public call の中だけで共有）
        self._row_nchan_cache = None
        self._row_signature_cache = None
        self._row_vcorr_cache = None
        self._group_cache = None
        self._cache_session_depth = 0

    def clear_caches(self, clear_target_grid: bool = False) -> None:
        """
        内部キャッシュを破棄する。通常利用では呼び出し不要。
        """
        self._base_axis_cache = {}
        self._invalidate_runtime_caches()
        if clear_target_grid:
            self.target_grid = None

    def _invalidate_runtime_caches(self) -> None:
        self._row_nchan_cache = None
        self._row_signature_cache = None
        self._row_vcorr_cache = None
        self._group_cache = None

    @contextmanager
    def _cache_session(self):
        outermost = (self._cache_session_depth == 0)
        if outermost:
            self._invalidate_runtime_caches()
        self._cache_session_depth += 1
        try:
            yield
        finally:
            self._cache_session_depth -= 1
            if outermost:
                self._invalidate_runtime_caches()

    def _ensure_row_cache(self) -> None:
        if self._group_cache is not None:
            return

        df = self.st.table
        meta = self.st.meta
        n_rows = len(df)
        is_list_data = isinstance(self.st.data, list)

        if is_list_data:
            nchan_cache = np.fromiter((len(self.st.data[i]) for i in range(n_rows)), dtype=int, count=n_rows)
        else:
            nchan_fixed = self.st.data.shape[1]
            nchan_cache = np.full(n_rows, nchan_fixed, dtype=int)

        signature_cache: List[tuple] = [()] * n_rows
        vcorr_cache = np.zeros(n_rows, dtype=float)
        groups: Dict[tuple, List[int]] = {}

        for i in range(n_rows):
            row = df.iloc[i]
            sig = get_axis_signature(row, meta, int(nchan_cache[i]))
            signature_cache[i] = sig
            # 元実装では v_corr は TOPO 系グループでのみ参照される。
            # 非 TOPO 行まで前計算すると、従来は無視されていた不正な
            # VFRAME/VELOSYS 値で失敗し得るため、ここでも同じ条件に揃える。
            if "TOPO" in str(sig[-1]).upper():
                vc = _get_vcorr_kms(row, meta, self.v_corr_col, row_only=True)
                if vc is None:
                    cand_txt = ", ".join([k for k in (self.v_corr_col, "VFRAME", "V_CORR_KMS", "VELOSYS") if k])
                    raise ValueError(
                        f"Row {i} has SPECSYS={str(sig[-1]).upper()} but no usable row-wise velocity correction column ({cand_txt})."
                    )
                vcorr_cache[i] = float(vc)
            else:
                vcorr_cache[i] = 0.0
            if sig not in groups:
                groups[sig] = []
            groups[sig].append(i)

        self._row_nchan_cache = nchan_cache
        self._row_signature_cache = tuple(signature_cache)
        self._row_vcorr_cache = np.nan_to_num(vcorr_cache, nan=0.0)
        self._group_cache = groups

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
        is_velocity_axis = (
            ctype.startswith("VRAD")
            or ctype.startswith("VOPT")
            or ctype.startswith("VELO")
            or (unit_norm in ("m/s", "km/s") and not ctype.startswith("FREQ"))
        )
        if ctype.startswith("FREQ") or unit_norm in ("hz",):
            dummy_meta = {
                "CRVAL1": crval, "CDELT1": cdelt, "CRPIX1": crpix,
                "CTYPE1": ctype, "RESTFREQ": restfreq, "CUNIT1": cunit
            }
            obs_freq = _generate_obs_freq_axis(dummy_meta, nchan)
        elif is_velocity_axis and restfreq > 0:
            dummy_meta = {
                "CRVAL1": crval, "CDELT1": cdelt, "CRPIX1": crpix,
                "CTYPE1": ctype, "RESTFREQ": restfreq, "CUNIT1": cunit
            }
            obs_freq = _generate_obs_freq_axis(dummy_meta, nchan)
                
        self._base_axis_cache[sig] = obs_freq
        return obs_freq

    def _group_by_signature(self) -> Dict[tuple, List[int]]:
        """Scantableの全行をAxisSignatureでグルーピングする"""
        self._ensure_row_cache()
        return self._group_cache if self._group_cache is not None else {}

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
            self._ensure_row_cache()
            v_corrs = self._row_vcorr_cache[np.asarray(idxs, dtype=int)]
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
        with self._cache_session():
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
        with self._cache_session():
            # 1. ターゲットグリッドの確定
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
            dv_tgt_abs = abs(float(target_grid.dv_kms))

            # 2. 結果格納用行列 (NaNで初期化)
            out_matrix = np.full((n_rows, n_tgt), np.nan, dtype=np.float32)

            # 3. グルーピング
            self._ensure_row_cache()
            groups = self._group_by_signature()

            # 4. グループごとに一括処理
            rest_freq_global = float(self.st.meta.get("RESTFRQ", self.st.meta.get("RESTFREQ", 0)))
            is_list_data = isinstance(self.st.data, list)

            # Internal policy:
            # - same/smaller dv -> preserve original np.interp path
            # - clearly larger dv -> exact overlap-weighted rebin
            coarse_rtol = 5.0e-3
            min_valid_coverage_frac = 0.999999

            for sig, idxs in groups.items():
                (_, _, _, _, _, _, sig_rest, sig_specsys) = sig
                rest_freq = sig_rest if sig_rest > 0 else rest_freq_global
                if rest_freq <= 0:
                    continue

                f_obs_axis = self._get_base_axis(sig)

                if is_list_data:
                    raw_spectra = [self.st.data[i] for i in idxs]
                    raw_block = np.vstack(raw_spectra)
                else:
                    raw_block = self.st.data[idxs]

                specsys = str(sig_specsys).upper()

                if "TOPO" in specsys:
                    v_corrs = self._row_vcorr_cache[np.asarray(idxs, dtype=int)]
                else:
                    v_corrs = np.zeros(len(idxs), dtype=float)

                k_factors = get_doppler_factor(v_corrs)
                k_inv_arr = 1.0 / k_factors
                term_b = (C_KMS * f_obs_axis) / rest_freq

                if len(term_b) > 1:
                    is_reversed = (C_KMS - term_b[0] * k_inv_arr[0]) > (C_KMS - term_b[-1] * k_inv_arr[0])

                    if is_reversed:
                        term_axis = term_b[::-1]
                    else:
                        term_axis = term_b

                    # row-wise source dv differs only through k_inv. Use this only to decide
                    # whether the target dv is clearly coarser than the native row dv.
                    dv_src_abs_arr = abs(float(term_axis[1] - term_axis[0])) * k_inv_arr

                    group_interp = bool(np.all(dv_tgt_abs <= dv_src_abs_arr * (1.0 + coarse_rtol)))
                    group_rebin = bool(np.all(dv_tgt_abs > dv_src_abs_arr * (1.0 + coarse_rtol)))

                    if group_interp:
                        if is_reversed:
                            for local_i, global_i in enumerate(idxs):
                                v_axis_row = C_KMS - term_axis * k_inv_arr[local_i]
                                out_matrix[global_i] = np.interp(v_tgt, v_axis_row, raw_block[local_i, ::-1], left=np.nan, right=np.nan)
                        else:
                            for local_i, global_i in enumerate(idxs):
                                v_axis_row = C_KMS - term_axis * k_inv_arr[local_i]
                                out_matrix[global_i] = np.interp(v_tgt, v_axis_row, raw_block[local_i], left=np.nan, right=np.nan)
                    elif group_rebin:
                        raw_block_use = raw_block[:, ::-1] if is_reversed else raw_block
                        v0_src_arr = C_KMS - term_axis[0] * k_inv_arr
                        dv_src_arr = (C_KMS - term_axis[1] * k_inv_arr) - v0_src_arr
                        coarse_factor_med = float(np.median(dv_tgt_abs / dv_src_abs_arr)) if len(dv_src_abs_arr) else 0.0
                        use_block_rebin = (len(idxs) >= 32 and raw_block_use.shape[1] <= 4096 and coarse_factor_med >= 3.0)
                        if use_block_rebin:
                            out_matrix[np.asarray(idxs, dtype=int)] = _rebin_uniform_block_to_regular_target(
                                raw_block_use,
                                v0_src_arr,
                                np.abs(dv_src_arr),
                                v0_tgt=float(target_grid.v0_kms),
                                dv_tgt=float(target_grid.dv_kms),
                                n_tgt=n_tgt,
                                min_valid_coverage_frac=min_valid_coverage_frac,
                            )
                        else:
                            for local_i, global_i in enumerate(idxs):
                                k_inv = float(k_inv_arr[local_i])
                                v0_src = float(C_KMS - term_axis[0] * k_inv)
                                dv_src = float((C_KMS - term_axis[1] * k_inv) - v0_src)
                                y_row = raw_block_use[local_i]
                                out_matrix[global_i] = _rebin_uniform_to_regular_target(
                                    y_row,
                                    v0_src,
                                    dv_src,
                                    v0_tgt=float(target_grid.v0_kms),
                                    dv_tgt=float(target_grid.dv_kms),
                                    n_tgt=n_tgt,
                                    min_valid_coverage_frac=min_valid_coverage_frac,
                                ).astype(np.float32, copy=False)
                    else:
                        for local_i, global_i in enumerate(idxs):
                            k_inv = float(k_inv_arr[local_i])
                            dv_src_abs = float(dv_src_abs_arr[local_i])
                            if dv_tgt_abs <= dv_src_abs * (1.0 + coarse_rtol):
                                v_axis_row = C_KMS - term_axis * k_inv
                                y_row = raw_block[local_i, ::-1] if is_reversed else raw_block[local_i]
                                out_matrix[global_i] = np.interp(v_tgt, v_axis_row, y_row, left=np.nan, right=np.nan)
                            else:
                                v0_src = float(C_KMS - term_axis[0] * k_inv)
                                dv_src = float((C_KMS - term_axis[1] * k_inv) - v0_src)
                                y_row = raw_block[local_i, ::-1] if is_reversed else raw_block[local_i]
                                out_matrix[global_i] = _rebin_uniform_to_regular_target(
                                    y_row,
                                    v0_src,
                                    dv_src,
                                    v0_tgt=float(target_grid.v0_kms),
                                    dv_tgt=float(target_grid.dv_kms),
                                    n_tgt=n_tgt,
                                    min_valid_coverage_frac=min_valid_coverage_frac,
                                ).astype(np.float32, copy=False)
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


# -----------------------------------------------------------------
# Internal exact rebin helpers (coarse dv only)
# -----------------------------------------------------------------

def _normalize_uniform_axis(y_src: np.ndarray, v0_src: float, dv_src: float) -> tuple[np.ndarray, float, float]:
    y = np.asarray(y_src, dtype=float)
    v0 = float(v0_src)
    dv = float(dv_src)
    if y.ndim != 1:
        raise ValueError("y_src must be 1D")
    if y.size <= 1 or dv == 0.0:
        return y, v0, abs(dv)
    if dv > 0.0:
        return y, v0, dv
    return y[::-1], v0 + (y.size - 1) * dv, -dv


def _build_piecewise_constant_prefix(weights: np.ndarray, dx: float) -> np.ndarray:
    w = np.asarray(weights, dtype=float)
    prefix = np.empty(w.size + 1, dtype=float)
    prefix[0] = 0.0
    if w.size > 0:
        np.cumsum(w * float(dx), dtype=float, out=prefix[1:])
    return prefix


def _eval_piecewise_constant_prefix(prefix: np.ndarray, weights: np.ndarray, x: np.ndarray, x0: float, dx: float) -> np.ndarray:
    q = np.asarray(x, dtype=float)
    n = int(np.asarray(weights).size)
    if n <= 0:
        return np.zeros_like(q, dtype=float)
    u = (q - float(x0)) / float(dx)
    np.clip(u, 0.0, float(n), out=u)
    k = np.floor(u).astype(np.int64)
    frac = u - k
    out = prefix[k].astype(float, copy=False)
    m = (k >= 0) & (k < n)
    if np.any(m):
        out[m] += frac[m] * float(dx) * np.asarray(weights, dtype=float)[k[m]]
    return out


def _rebin_uniform_to_regular_target(
    y_src: np.ndarray,
    v0_src: float,
    dv_src: float,
    *,
    v0_tgt: float,
    dv_tgt: float,
    n_tgt: int,
    min_valid_coverage_frac: float = 0.999999,
) -> np.ndarray:
    y, v0_asc, dv_asc = _normalize_uniform_axis(y_src, v0_src, dv_src)
    out = np.full(int(n_tgt), np.nan, dtype=float)
    if out.size == 0 or y.size == 0:
        return out
    dv_tgt_abs = abs(float(dv_tgt))
    if y.size <= 1 or dv_asc == 0.0 or dv_tgt_abs == 0.0:
        v_src = v0_asc + np.arange(y.size, dtype=float) * (dv_asc if dv_asc != 0.0 else 1.0)
        v_tgt = float(v0_tgt) + np.arange(int(n_tgt), dtype=float) * float(dv_tgt)
        return interp_to_vgrid(v_src, y, v_tgt)

    tgt_left = (float(v0_tgt) - 0.5 * dv_tgt_abs) + np.arange(int(n_tgt), dtype=float) * dv_tgt_abs
    tgt_right = tgt_left + dv_tgt_abs

    src_edge0 = float(v0_asc - 0.5 * dv_asc)
    src_edge1 = float(src_edge0 + y.size * dv_asc)

    finite = np.isfinite(y)
    all_finite = bool(np.all(finite))

    if all_finite:
        weights_sig = y
        prefix_sig = _build_piecewise_constant_prefix(weights_sig, dv_asc)
        left_sig = _eval_piecewise_constant_prefix(prefix_sig, weights_sig, tgt_left, src_edge0, dv_asc)
        right_sig = _eval_piecewise_constant_prefix(prefix_sig, weights_sig, tgt_right, src_edge0, dv_asc)
        valid_width = np.minimum(tgt_right, src_edge1) - np.maximum(tgt_left, src_edge0)
        np.maximum(valid_width, 0.0, out=valid_width)
    else:
        y_fill = np.where(finite, y, 0.0)
        valid_fill = finite.astype(float, copy=False)
        prefix_sig = _build_piecewise_constant_prefix(y_fill, dv_asc)
        prefix_val = _build_piecewise_constant_prefix(valid_fill, dv_asc)
        left_sig = _eval_piecewise_constant_prefix(prefix_sig, y_fill, tgt_left, src_edge0, dv_asc)
        right_sig = _eval_piecewise_constant_prefix(prefix_sig, y_fill, tgt_right, src_edge0, dv_asc)
        left_val = _eval_piecewise_constant_prefix(prefix_val, valid_fill, tgt_left, src_edge0, dv_asc)
        right_val = _eval_piecewise_constant_prefix(prefix_val, valid_fill, tgt_right, src_edge0, dv_asc)
        valid_width = right_val - left_val

    good = valid_width >= float(min_valid_coverage_frac) * dv_tgt_abs
    if np.any(good):
        out[good] = (right_sig[good] - left_sig[good]) / dv_tgt_abs
    return out



def _rebin_uniform_block_to_regular_target(
    y_block: np.ndarray,
    v0_src_arr: np.ndarray,
    dv_src_arr: np.ndarray,
    *,
    v0_tgt: float,
    dv_tgt: float,
    n_tgt: int,
    min_valid_coverage_frac: float = 0.999999,
) -> np.ndarray:
    """Exact coarse rebin for a block of rows with row-wise affine source axes.

    This private helper is used only on the coarse path. Each row must already be
    ordered so that the source velocity axis is ascending.
    """
    y2d = np.asarray(y_block, dtype=float)
    if y2d.ndim != 2:
        raise ValueError('y_block must be 2D')
    n_rows, n_src = y2d.shape
    out = np.full((n_rows, int(n_tgt)), np.nan, dtype=np.float32)
    if n_rows == 0 or n_src == 0 or int(n_tgt) <= 0:
        return out

    dv_src_abs = np.abs(np.asarray(dv_src_arr, dtype=float))
    v0_src = np.asarray(v0_src_arr, dtype=float)
    dv_tgt_abs = abs(float(dv_tgt))
    if n_src <= 1 or dv_tgt_abs == 0.0 or np.any(dv_src_abs == 0.0):
        for i in range(n_rows):
            out[i] = _rebin_uniform_to_regular_target(
                y2d[i],
                float(v0_src[i]),
                float(dv_src_abs[i]),
                v0_tgt=float(v0_tgt),
                dv_tgt=float(dv_tgt),
                n_tgt=int(n_tgt),
                min_valid_coverage_frac=float(min_valid_coverage_frac),
            ).astype(np.float32, copy=False)
        return out

    finite = np.isfinite(y2d)
    y_fill = np.where(finite, y2d, 0.0)
    valid_fill = finite.astype(float)

    dv2 = dv_src_abs[:, None]
    src_edge0 = v0_src[:, None] - 0.5 * dv2

    prefix_sig = np.zeros((n_rows, n_src + 1), dtype=float)
    np.cumsum(y_fill * dv2, axis=1, out=prefix_sig[:, 1:])

    prefix_val = np.zeros((n_rows, n_src + 1), dtype=float)
    np.cumsum(valid_fill * dv2, axis=1, out=prefix_val[:, 1:])

    j = np.arange(int(n_tgt), dtype=float)
    left_1d = (float(v0_tgt) - 0.5 * dv_tgt_abs) + j * dv_tgt_abs
    right_1d = left_1d + dv_tgt_abs

    def eval_prefix(prefix_2d: np.ndarray, weights_2d: np.ndarray, q_1d: np.ndarray) -> np.ndarray:
        q2 = q_1d[None, :]
        u = (q2 - src_edge0) / dv2
        np.clip(u, 0.0, float(n_src), out=u)
        k = np.floor(u).astype(np.int64)
        frac = u - k
        val = np.take_along_axis(prefix_2d, k, axis=1)
        k_safe = np.clip(k, 0, n_src - 1)
        w = np.take_along_axis(weights_2d, k_safe, axis=1)
        val += frac * dv2 * w
        return val

    left_sig = eval_prefix(prefix_sig, y_fill, left_1d)
    right_sig = eval_prefix(prefix_sig, y_fill, right_1d)
    left_val = eval_prefix(prefix_val, valid_fill, left_1d)
    right_val = eval_prefix(prefix_val, valid_fill, right_1d)

    valid_width = right_val - left_val
    good = valid_width >= float(min_valid_coverage_frac) * dv_tgt_abs
    if np.any(good):
        tmp = (right_sig - left_sig) / dv_tgt_abs
        out[good] = tmp[good].astype(np.float32, copy=False)
    return out

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



def _slice_row_wcs_columns(table: pd.DataFrame, ch_start: Optional[int], ch_stop: Optional[int], new_naxis1: int) -> pd.DataFrame:
    out = table.copy()
    if ch_start is None and ch_stop is None:
        return out
    s = 0 if ch_start is None else int(ch_start)
    if 'CRVAL1' in out.columns and 'CDELT1' in out.columns:
        out['CRVAL1'] = pd.to_numeric(out['CRVAL1'], errors='coerce') + float(s) * pd.to_numeric(out['CDELT1'], errors='coerce')
    if 'NAXIS1' in out.columns:
        out['NAXIS1'] = int(new_naxis1)
    if 'CRPIX1' not in out.columns:
        out['CRPIX1'] = 1.0
    return out


def _slice_meta_wcs_if_present(meta: dict, ch_start: Optional[int], ch_stop: Optional[int], new_naxis1: int) -> dict:
    m = dict(meta)
    if ch_start is None and ch_stop is None:
        return m
    s = 0 if ch_start is None else int(ch_start)
    # Respect row-first policy: only update meta WCS if those keys actually exist.
    if 'CRVAL1' in m and 'CDELT1' in m:
        m['CRVAL1'] = float(m['CRVAL1']) + float(s) * float(m['CDELT1'])
    if 'NAXIS1' in m or new_naxis1 is not None:
        m['NAXIS1'] = int(new_naxis1)
    return m


def _collect_restfreq_values(table: pd.DataFrame, meta: dict) -> list[float]:
    vals: list[float] = []
    for key in ('RESTFRQ', 'RESTFREQ'):
        if key in table.columns:
            s = pd.to_numeric(table[key], errors='coerce')
            vals.extend([float(x) for x in s[np.isfinite(s)] if float(x) > 0])
    for key in ('RESTFRQ', 'RESTFREQ'):
        if key in meta:
            try:
                v = float(meta[key])
            except Exception:
                continue
            if np.isfinite(v) and v > 0:
                vals.append(v)
    return vals


def _resolve_uniform_restfreq(table: pd.DataFrame, meta: dict, *, rtol: float = 1.0e-12, atol: float = 0.0) -> float:
    vals = _collect_restfreq_values(table, meta)
    if len(vals) == 0:
        raise ValueError('RESTFRQ/RESTFREQ is required for velocity regrid output.')
    ref = float(vals[0])
    for v in vals[1:]:
        if not np.isclose(float(v), ref, rtol=rtol, atol=atol):
            raise ValueError('velocity_regrid does not support heterogeneous RESTFRQ/RESTFREQ across rows. Please split by line or pass rest_freq=... explicitly.')
    return ref


# =========================================================
# 4. Public API: row-preserving velocity regrid
# =========================================================

def _slice_data_and_rows(sc, idxs: np.ndarray, ch_start: Optional[int], ch_stop: Optional[int]):
    """Return a row-selected / channel-sliced shallow copy of Scantable-like contents.

    - preserves row order given by idxs
    - supports both ndarray and list-of-arrays storage
    - updates row-wise CRVAL1/CRPIX1/NAXIS1 when channel slicing is requested
    - does not require meta/header WCS when row WCS exists
    """
    from .fitsio import Scantable

    idxs = np.asarray(idxs, dtype=int)
    table = sc.table.iloc[idxs].reset_index(drop=True).copy()
    meta = dict(sc.meta)

    if isinstance(sc.data, list):
        data_sel = [np.asarray(sc.data[i], dtype=float) for i in idxs]
        if ch_start is not None or ch_stop is not None:
            sliced = []
            for j, arr in enumerate(data_sel):
                nchan = int(arr.shape[0])
                s = 0 if ch_start is None else int(ch_start)
                e = nchan if ch_stop is None else int(ch_stop)
                if not (0 <= s < e <= nchan):
                    raise ValueError(f"Invalid channel slice {s}:{e} for row {j} with nchan={nchan}")
                sliced.append(np.asarray(arr[s:e], dtype=float))
            data_sel = sliced
            nchan0 = int(len(data_sel[0])) if len(data_sel) > 0 else 0
            table = _slice_row_wcs_columns(table, ch_start, ch_stop, nchan0)
            meta = _slice_meta_wcs_if_present(meta, ch_start, ch_stop, nchan0)
        return stamp_scantable_code_provenance(
            Scantable(meta=meta, data=data_sel, table=table, history=dict(sc.history)),
            stage="select_rows_and_channels",
        )

    data_arr = np.asarray(sc.data)
    data_sel = np.asarray(data_arr[idxs], dtype=float)
    if ch_start is not None or ch_stop is not None:
        nchan = int(data_sel.shape[1])
        s = 0 if ch_start is None else int(ch_start)
        e = nchan if ch_stop is None else int(ch_stop)
        if not (0 <= s < e <= nchan):
            raise ValueError(f"Invalid channel slice {s}:{e} for nchan={nchan}")
        data_sel = np.asarray(data_sel[:, s:e], dtype=float)
        nchan_new = int(data_sel.shape[1])
        table = _slice_row_wcs_columns(table, ch_start, ch_stop, nchan_new)
        meta = _slice_meta_wcs_if_present(meta, ch_start, ch_stop, nchan_new)
    return stamp_scantable_code_provenance(
        Scantable(meta=meta, data=data_sel, table=table, history=dict(sc.history)),
        stage="select_rows_and_channels",
    )


def _prepare_row_selection(n_rows: int, rows=None, exclude_rows=None, max_dumps: int = 0) -> np.ndarray:
    from .scantable_utils import _parse_row_selector
    idxs = _parse_row_selector(rows, n_rows)
    if exclude_rows is not None:
        ex = set(_parse_row_selector(exclude_rows, n_rows).tolist())
        idxs = np.asarray([i for i in idxs.tolist() if i not in ex], dtype=int)
    if int(max_dumps or 0) > 0:
        idxs = idxs[: int(max_dumps)]
    return idxs


def _drop_velocity_correction_columns(table: pd.DataFrame, meta: dict, preferred_key: str) -> tuple[pd.DataFrame, dict]:
    out = table.copy()
    m = dict(meta)
    related_cols = []
    preferred = str(preferred_key or '').strip().upper()
    for key in (preferred, 'VELOSYS', 'VFRAME', 'V_CORR_KMS'):
        if key and key not in related_cols:
            related_cols.append(key)
    for key in related_cols:
        if key in out.columns:
            out = out.drop(columns=[key])
        m.pop(key, None)
    return out, m


def run_velocity_regrid(
    input_data,
    output_path: Optional[str] = None,
    *,
    rows=None,
    exclude_rows=None,
    vmin_kms: float,
    vmax_kms: float,
    dv_kms: float,
    v_corr_col: str = 'VFRAME',
    rest_freq: Optional[float] = None,
    ch_start: Optional[int] = None,
    ch_stop: Optional[int] = None,
    max_dumps: int = 0,
    overwrite: bool = True,
    fill_value: float = np.nan,
    keep_row_order: bool = True,
    drop_allnan_rows: bool = False,
    history_tag: str = 'velocity_regrid',
):
    """Regrid each row onto a common LSRK velocity axis without coadding.

    Parameters
    ----------
    input_data : str | Scantable
        Input scantable or path.
    output_path : str, optional
        Output FITS path.
    rows, exclude_rows : selector, optional
        Row selection compatible with ``_parse_row_selector``.
    vmin_kms, vmax_kms, dv_kms : float
        Target common LSRK velocity grid in km/s.
    v_corr_col : str, default 'VFRAME'
        Row-wise correction column for TOPO input. ``VELOSYS`` / ``VFRAME`` are
        interpreted as m/s; ``*_KMS`` columns are interpreted as km/s.
    rest_freq : float, optional
        Override rest frequency [Hz].
    ch_start, ch_stop : int, optional
        Optional pre-slice of native spectral channels before regridding.
    max_dumps : int, default 0
        If >0, keep only first ``max_dumps`` selected rows.
    overwrite : bool, default True
        Passed to ``write_scantable`` when ``output_path`` is given.
    fill_value : float, default np.nan
        Fill value outside interpolation support.
    keep_row_order : bool, default True
        Currently row order is always preserved; kept for explicit API.
    drop_allnan_rows : bool, default False
        Drop rows that are all NaN after regridding.
    history_tag : str, default 'velocity_regrid'
        Stored in history.
    """
    from .fitsio import Scantable, read_scantable, write_scantable
    from .restfreq import apply_restfreq_override

    if dv_kms is None or float(dv_kms) <= 0:
        raise ValueError('dv_kms must be positive')

    if not bool(keep_row_order):
        raise NotImplementedError('keep_row_order=False is not implemented; row order is preserved.')

    if isinstance(input_data, Scantable):
        sc_in = input_data.copy()
    elif isinstance(input_data, str):
        sc_in = read_scantable(input_data)
    else:
        raise TypeError('input_data must be a path or Scantable')

    idxs = _prepare_row_selection(len(sc_in.table), rows=rows, exclude_rows=exclude_rows, max_dumps=max_dumps)
    if idxs.size == 0:
        raise ValueError('No rows selected for velocity regrid.')

    sc_sel = _slice_data_and_rows(sc_in, idxs, ch_start=ch_start, ch_stop=ch_stop)

    # Optional REST frequency override before Standardizer so source->target conversion is consistent.
    if rest_freq is not None:
        apply_restfreq_override(sc_sel.meta, sc_sel.table, float(rest_freq), require_wcs_for_vrad=False)

    target_grid = make_vgrid(float(vmin_kms), float(vmax_kms), float(dv_kms))
    std = Standardizer(sc_sel, target_grid=target_grid, v_corr_col=v_corr_col)
    out_matrix, v_tgt = std.get_matrix(target_grid=target_grid)

    keep = None
    if drop_allnan_rows:
        keep = np.any(np.isfinite(out_matrix), axis=1)
        out_matrix = out_matrix[keep]

    if not np.isnan(fill_value):
        out_matrix = np.where(np.isfinite(out_matrix), out_matrix, float(fill_value)).astype(np.float32, copy=False)

    out_table = sc_sel.table.copy().reset_index(drop=True)
    if keep is not None:
        out_table = out_table.loc[keep].reset_index(drop=True)
    meta_out = dict(sc_sel.meta)

    # Preserve unrelated metadata/columns, update only spectral-axis/frame related fields.
    out_table, meta_out = _drop_velocity_correction_columns(out_table, meta_out, preferred_key=v_corr_col)

    rest0 = _resolve_uniform_restfreq(sc_sel.table, sc_sel.meta)

    resolved_ssysobs = [
        _get_ssysobs(out_table.iloc[i], sc_sel.meta)
        for i in range(len(out_table))
    ]
    ssysobs_vals = pd.Series(resolved_ssysobs, dtype=object)
    ssysobs_vals = ssysobs_vals.dropna().astype(str).str.strip()
    ssysobs_vals = ssysobs_vals[ssysobs_vals != '']
    ssysobs_unique = pd.unique(ssysobs_vals)

    n_tgt = int(len(v_tgt))
    out_table['CTYPE1'] = 'VRAD'
    out_table['CUNIT1'] = 'm/s'
    out_table['CRVAL1'] = 1000.0 * (float(v_tgt[0]) if n_tgt > 0 else float(vmin_kms))
    out_table['CDELT1'] = 1000.0 * float(target_grid.dv_kms)
    out_table['CRPIX1'] = float(target_grid.crpix1)
    out_table['NAXIS1'] = n_tgt
    out_table['SPECSYS'] = 'LSRK'
    out_table['SSYSOBS'] = resolved_ssysobs
    out_table['RESTFRQ'] = float(rest0)
    out_table['RESTFREQ'] = float(rest0)
    out_table['VELDEF'] = 'RADIO'
    out_table['REGRID_DONE'] = True
    out_table['REGRID_FRAME'] = 'LSRK'

    meta_out.update(dict(
        CTYPE1='VRAD',
        CUNIT1='m/s',
        CRVAL1=1000.0 * (float(v_tgt[0]) if n_tgt > 0 else float(vmin_kms)),
        CDELT1=1000.0 * float(target_grid.dv_kms),
        CRPIX1=float(target_grid.crpix1),
        NAXIS1=n_tgt,
        SPECSYS='LSRK',
        RESTFRQ=float(rest0),
        RESTFREQ=float(rest0),
        VELDEF='RADIO',
        TIMESYS=str(meta_out.get('TIMESYS', 'UTC')),
    ))
    if len(ssysobs_unique) == 1:
        meta_out['SSYSOBS'] = str(ssysobs_unique[0]).upper()
    else:
        meta_out.pop('SSYSOBS', None)
    meta_out.pop('VELOSYS', None)
    meta_out.pop('VFRAME', None)
    meta_out.pop(str(v_corr_col or '').strip().upper(), None)

    hist = dict(sc_sel.history)
    hist[history_tag] = dict(
        stage='velocity_regrid',
        frame='LSRK',
        kind='row_preserving_regrid',
        vmin_kms=float(vmin_kms),
        vmax_kms=float(vmax_kms),
        dv_kms=float(dv_kms),
        fill_value=None if np.isnan(fill_value) else float(fill_value),
        source_specsys=sorted({str(x).upper() for x in pd.Series(sc_sel.table.get('SPECSYS', sc_sel.meta.get('SPECSYS', ''))).dropna().astype(str).tolist()}) if len(sc_sel.table) else [str(sc_sel.meta.get('SPECSYS', '')).upper()],
        source_vcorr_col=str(v_corr_col),
        rows=str(rows),
        exclude_rows=str(exclude_rows),
        ch_start=None if ch_start is None else int(ch_start),
        ch_stop=None if ch_stop is None else int(ch_stop),
        max_dumps=int(max_dumps or 0),
        keep_row_order=bool(keep_row_order),
        drop_allnan_rows=bool(drop_allnan_rows),
    )

    res = stamp_scantable_code_provenance(
        Scantable(meta=meta_out, data=np.asarray(out_matrix, dtype=np.float32), table=out_table, history=hist),
        stage="run_velocity_regrid",
    )
    if output_path:
        write_scantable(output_path, res, overwrite=overwrite)
    return res
