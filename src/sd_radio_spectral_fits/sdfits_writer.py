#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
互換性最優先 SDFITS Writer（Single-Dish Spectral FITS）
==========================================================

狙い
----
- 汎用（asap系/単一鏡スペクトル処理系）で「読めて」「解析に使える」ことを最優先にした SDFITS 風FITSを書き出します。
- 先の版で良かった点（全パラメータの意味・単位・具体例がコード内で自明、厳密バリデーション、巨大OTFのparts分割）を維持し、
  互換性の観点で「SDFITSっぽいカラム名・ヘッダ名」を採用します。

重要（互換性のための実装方針）
------------------------------
(1) 1行=1スペクトル（dump/未積分/積分済も同一形式）
(2) スペクトル配列カラムは SDFITS/GBT系でよくある 'SPECTRUM' を採用し、
    互換保険として同一内容の別名 'DATA' も同時に持てるオプションを用意（既定は両方出力）。
    ※冗長だが「読めない」を避けるための安全策。
(3) 座標・観測モード・TSYS・積分時間・周波数WCS（またはFREQベクトル）を必ず保持。
(4) 追加で、観測意図座標/オフセット/ハードウェア生値/補正/τ 等は
    既存機能を省略せず、SDFITS標準外は HIERARCH ではなく「通常のカラム」として残す。

巨大OTF対策
-----------
- FITS Binary table は安全な追記が難しいため、「parts分割（chunk_size指定）」を推奨。
- partsモードでは base_part0001.fits ... を出し、base_manifest.json を作る。

依存
----
numpy, astropy

使い方
------
1) site/info を作る
2) writer = SDRadioSpectralSDFITSWriter(...)
3) writer.add_row(...)
4) smallなら writer.write("out.fits")
   large(parts)なら writer.close()（残りflush+manifest）

このコードは「見れば仕様が分かる」ことを重視し、全項目を SCHEMA に列挙します。
writer.print_schema() で一覧表示できます。
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from astropy.io import fits
from .version import __version__, SWNAME
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, FK5, FK4, Galactic
import astropy.units as u


Number = Union[int, float, np.number]


# =============================================================================
# 1) 列挙値（取り得る値をすべて明示）
# =============================================================================
RADESYS_VALUES = ("ICRS", "FK5", "FK4", "FK4-NO-E", "GALACTIC")
SPECSYS_VALUES = ("TOPOCENT", "LSRK", "BARYCENT", "GEOCENTR", "HELIOCEN")
VELDEF_VALUES  = ("RADIO", "OPTICAL")

COORDTYPE_VALUES = ("RADEC", "GALACTIC", "AZEL", "AZEL_GEO")

# 標準的な単一鏡解析/単一鏡で扱いやすい短い固定コード推奨
OBSMODE_VALUES = ("ON", "OFF", "HOT", "SKY", "FS_SIG", "FS_REF", "CAL_ON", "CAL_OFF", "UNKNOWN")

CALSTAT_VALUES = ("RAW", "TASTAR", "TA", "TMB", "INTEGRATED")

EFFSTAT_VALUES = ("ASSUMED", "MEASURED", "UNKNOWN")


# =============================================================================
# 2) 仕様一覧（SCHEMA）— 省略せず全項目
# =============================================================================
SCHEMA: Dict[str, Dict[str, Any]] = {
    # --- writer ---
    "writer.n_chan": {
        "where": "SDRadioSpectralSDFITSWriter.__init__",
        "type": "int",
        "unit": "",
        "required": True,
        "meaning": "チャンネル数。n_chan==1 なら連続波/総電力の1点扱い。",
        "examples": [1024, 1],
        "notes": "n_chanが混在するデータは同一FITSに入れない。"
    },
    "writer.store_freq_column": {
        "where": "SDRadioSpectralSDFITSWriter.__init__",
        "type": "bool",
        "unit": "",
        "required": False,
        "meaning": "Trueなら各行に FREQ(Hzベクトル) を持たせる。FalseならWCS(ヘッダ)のみ。",
        "examples": [False, True],
        "notes": "周波数軸が行ごとに変わる可能性がある場合は True を推奨。"
    },
    "writer.duplicate_data_columns": {
        "where": "SDRadioSpectralSDFITSWriter.__init__",
        "type": "bool",
        "unit": "",
        "required": False,
        "meaning": "Trueならスペクトル配列を SPECTRUM と DATA の両方に同内容で出力（互換保険）。",
        "examples": [True, False],
        "notes": "互換性最優先では True 推奨（既定True）。サイズは増える。"
    },
    "writer.chunk_size": {
        "where": "SDRadioSpectralSDFITSWriter.__init__",
        "type": "Optional[int]",
        "unit": "rows",
        "required": False,
        "meaning": "巨大OTF用。指定すると chunkごとに partXXXX.fits に分割書き出し。",
        "examples": [None, 50000],
        "notes": "FITS追記の不確実性を避けるため parts方式を推奨。"
    },
    "writer.out_basename": {
        "where": "SDRadioSpectralSDFITSWriter.__init__",
        "type": "Optional[str]",
        "unit": "",
        "required": "chunk_size is not None",
        "meaning": "parts方式の基底名。例 'otf' -> otf_part0001.fits を生成。",
        "examples": ["otf_run1"],
        "notes": ""
    },

    # --- site ---
    "site.lat_deg": {"where": "Site", "type": "float", "unit": "deg", "required": True,
                     "meaning": "観測地点緯度（測地）。", "examples": [35.9400], "notes": ""},
    "site.lon_deg": {"where": "Site", "type": "float", "unit": "deg (east+)", "required": True,
                     "meaning": "観測地点経度（東経+）。", "examples": [138.4720], "notes": ""},
    "site.elev_m": {"where": "Site", "type": "float", "unit": "m", "required": True,
                    "meaning": "観測地点標高。", "examples": [1350.0], "notes": ""},

    # --- dataset info ---
    "info.telescope": {"where": "DatasetInfo", "type": "str", "unit": "", "required": True,
                       "meaning": "望遠鏡名（TELESCOP）。", "examples": ["MyDish", "NRO45m"], "notes": ""},
    "info.observer": {"where": "DatasetInfo", "type": "str", "unit": "", "required": False,
                      "meaning": "観測者（OBSERVER）。", "examples": ["A. Observer"], "notes": ""},
    "info.project": {"where": "DatasetInfo", "type": "str", "unit": "", "required": False,
                     "meaning": "プロジェクトID（PROJID/PROJECT）。", "examples": ["P001"], "notes": ""},
    "info.object_name": {"where": "DatasetInfo", "type": "str", "unit": "", "required": False,
                         "meaning": "代表天体名（OBJECT）。行ごとに別OBJECTも持てる。", "examples": ["OrionKL"], "notes": ""},

    # --- RA/DEC frame ---
    "info.radesys": {"where": "DatasetInfo", "type": "str (enum)", "unit": "", "required": False,
                     "meaning": "RA/DEC列の基準系（RADESYS）。", "options": list(RADESYS_VALUES),
                     "examples": ["ICRS", "FK5"], "notes": "推奨ICRS。FK5/FK4ならequinox明示。"},
    "info.equinox": {"where": "DatasetInfo", "type": "float", "unit": "year", "required": False,
                     "meaning": "FK5/FK4の分点（EQUINOX）。", "examples": [2000.0, 1950.0],
                     "notes": "ICRSでも2000.0としておく運用は問題ない。"},
    "info.src_radesys": {"where": "DatasetInfo", "type": "str (enum)", "unit": "", "required": False,
                         "meaning": "SRCFRAME='RADEC' の SRC_LONG/LAT の基準系（HIERARCH SRC_RADESYS）。",
                         "options": list(RADESYS_VALUES), "examples": ["ICRS", "FK4"], "notes": ""},
    "info.src_equinox": {"where": "DatasetInfo", "type": "float", "unit": "year", "required": False,
                         "meaning": "SRCFRAME='RADEC' の分点（HIERARCH SRC_EQUINOX）。",
                         "examples": [2000.0, 1950.0], "notes": "B1950ならFK4/1950.0。"},
    "info.bmaj_deg": {"where": "DatasetInfo", "type": "float", "unit": "deg", "required": False,
                      "meaning": "ビーム長径（BMAJ）。未知ならNaN。", "examples": [0.002, np.nan],
                      "notes": "arcsecなら/3600。"},
    "info.bmin_deg": {"where": "DatasetInfo", "type": "float", "unit": "deg", "required": False,
                      "meaning": "ビーム短径（BMIN）。未知ならNaN。", "examples": [0.002, np.nan], "notes": ""},
    "info.bpa_deg": {"where": "DatasetInfo", "type": "float", "unit": "deg", "required": False,
                     "meaning": "ビームPA（BPA）。未知ならNaN。", "examples": [0.0, np.nan], "notes": ""},

    # --- efficiency ---
    "info.eff.beameff": {"where": "Efficiency", "type": "float", "unit": "", "required": False,
                         "meaning": "主ビーム能率 η_mb（BEAMEFF）。", "examples": [0.45, np.nan],
                         "notes": "不明ならNaN+effstat。"},
    "info.eff.apereff": {"where": "Efficiency", "type": "float", "unit": "", "required": False,
                         "meaning": "開口能率 η_a（APEREFF）。", "examples": [0.30, np.nan], "notes": ""},
    "info.eff.mooneff": {"where": "Efficiency", "type": "float", "unit": "", "required": False,
                         "meaning": "月能率（MOONEFF）。", "examples": [0.85, np.nan], "notes": ""},
    "info.eff.suneff": {"where": "Efficiency", "type": "float", "unit": "", "required": False,
                        "meaning": "太陽能率（SUNEFF）。", "examples": [0.90, np.nan], "notes": ""},
    "info.eff.effstat": {"where": "Efficiency", "type": "str (enum)", "unit": "", "required": False,
                         "meaning": "能率の信頼度（EFFSTAT）。", "options": list(EFFSTAT_VALUES),
                         "examples": ["UNKNOWN", "ASSUMED", "MEASURED"], "notes": ""},

    # --- interpretation flags ---
    "info.refr_included_in_corr": {"where": "DatasetInfo", "type": "bool", "unit": "", "required": False,
                                   "meaning": "CORR_* が屈折を含むか（HIERARCH REFR_INC）。", "examples": [True, False], "notes": ""},
    "info.doppler_tracking_applied": {"where": "DatasetInfo", "type": "bool", "unit": "", "required": False,
                                      "meaning": "オンラインドップラー追尾を行ったか（HIERARCH DOPPLER）。", "examples": [True, False], "notes": ""},

    # --- spectral axis (uniform WCS) ---
    "info.spectral_axis.crval1_hz": {"where": "SpectralAxisUniform", "type": "float", "unit": "Hz",
                                     "required": "store_freq_column==False",
                                     "meaning": "基準ピクセル周波数（CRVAL1）。", "examples": [110.201354e9],
                                     "notes": "推奨：TOPOCENT観測周波数で保存。"},
    "info.spectral_axis.cdelt1_hz": {"where": "SpectralAxisUniform", "type": "float", "unit": "Hz",
                                     "required": "store_freq_column==False",
                                     "meaning": "チャンネル間隔（CDELT1）。", "examples": [61e3, -61e3], "notes": ""},
    "info.spectral_axis.crpix1": {"where": "SpectralAxisUniform", "type": "float", "unit": "pixel(1-based)",
                                  "required": False, "meaning": "基準ピクセル（CRPIX1）。", "examples": [1.0], "notes": ""},
    "info.spectral_axis.restfreq_hz": {"where": "SpectralAxisUniform", "type": "float", "unit": "Hz",
                                       "required": True,
                                       "meaning": "静止周波数（RESTFREQ/RESTFRQ）。速度軸変換に必須。", "examples": [110.201354e9],
                                       "notes": "未確定でも可だが互換性/再現性のため確定推奨。"},
    "info.spectral_axis.specsys": {"where": "SpectralAxisUniform", "type": "str(enum)", "unit": "",
                                   "required": False, "meaning": "スペクトル基準系（SPECSYS）。", "options": list(SPECSYS_VALUES),
                                   "examples": ["TOPOCENT", "LSRK"], "notes": "推奨TOPOCENT。"},
    "info.spectral_axis.ssysobs": {"where": "SpectralAxisUniform", "type": "str(enum)", "unit": "",
                                   "required": False, "meaning": "観測者の基準系（SSYSOBS）。通常TOPOCENT。", "options": list(SPECSYS_VALUES),
                                   "examples": ["TOPOCENT"], "notes": "SPECSYSがLSRKでも、SSYSOBSはTOPOCENTとするのが一般的。"},
    "info.spectral_axis.veldef": {"where": "SpectralAxisUniform", "type": "str(enum)", "unit": "",
                                  "required": False, "meaning": "速度定義（VELDEF）。", "options": list(VELDEF_VALUES),
                                  "examples": ["RADIO"], "notes": "電波ではRADIO推奨。"},
    "info.spectral_axis.ctype1": {"where": "SpectralAxisUniform", "type": "str", "unit": "",
                                  "required": False, "meaning": "WCS軸種別（CTYPE1）。推奨FREQ。", "examples": ["FREQ"], "notes": ""},
    "info.spectral_axis.cunit1": {"where": "SpectralAxisUniform", "type": "str", "unit": "",
                                  "required": False, "meaning": "WCS単位（CUNIT1）。推奨Hz。", "examples": ["Hz"], "notes": ""},
    "info.spectral_axis.refchan": {"where": "SpectralAxisUniform", "type": "int", "unit": "channel(1-based)",
                                   "required": False, "meaning": "参照チャンネル（REFCHAN）。", "examples": [1, 512], "notes": ""},

    # --- shared meta ---
    "info.shared_meta": {"where": "DatasetInfo", "type": "dict[str, scalar]", "unit": "", "required": False,
                         "meaning": "共通メタを HIERARCH KEY=VALUE として保存。", "examples": [{"TAU_DEF": "zenith tau"}],
                         "notes": "値は FITS scalar推奨（str/int/float/bool）。"},

    # --- row core ---
    "row.time_mjd": {"where": "add_row", "type": "float", "unit": "day (MJD UTC)", "required": True,
                     "meaning": "MJD(UTC)。TIME列。", "examples": [60200.0], "notes": ""},
    "row.scanid": {"where": "add_row", "type": "int", "unit": "", "required": True,
                   "meaning": "較正参照グループ。OFF更新で増やす。SCAN列。", "examples": [1, 2], "notes": ""},
    "row.subscan": {"where": "add_row", "type": "int", "unit": "", "required": True,
                    "meaning": "幾何学スキャンID（OTFライン等）。SUBSCAN列。", "examples": [10, 11], "notes": ""},
    "row.intgrp": {"where": "add_row", "type": "int", "unit": "", "required": True,
                   "meaning": "dump束など“同一位置のまとまり”。INTGRP列。", "examples": [100, 101], "notes": ""},
    "row.obsmode": {"where": "add_row", "type": "str(enum)", "unit": "", "required": True,
                    "meaning": "観測モード（OBSMODE列）。", "options": list(OBSMODE_VALUES),
                    "examples": ["ON", "OFF", "HOT"], "notes": "長名は避ける（固定長事故回避）。"},
    "row.object": {"where": "add_row", "type": "str", "unit": "", "required": False,
                   "meaning": "行ごとの天体名（OBJECT列）。省略時は info.object_name。", "examples": ["OrionKL"], "notes": ""},

    "row.data": {"where": "add_row", "type": "np.ndarray(float) or float", "unit": "K or counts (calstat依存)",
                 "required": True, "meaning": "スペクトル配列（SPECTRUM列、互換保険でDATA列も可）。",
                 "examples": ["np.ndarray(n_chan)", 12.3], "notes": ""},
    "row.flag": {"where": "add_row", "type": "np.ndarray(bool)", "unit": "", "required": False,
                 "meaning": "チャンネルマスク（FLAG列）。未指定は全False。", "examples": ["mask[500:510]=True"], "notes": ""},
    "row.flagrow": {"where": "add_row", "type": "int(0/1)", "unit": "", "required": False,
                    "meaning": "行全体フラグ（FLAGROW列）。", "examples": [0, 1], "notes": ""},

    "row.exposure_s": {"where": "add_row", "type": "float", "unit": "s", "required": True,
                       "meaning": "積分時間（EXPOSURE列）。", "examples": [0.1, 20.0], "notes": ""},
    "row.calstat": {"where": "add_row", "type": "str(enum)", "unit": "", "required": False,
                    "meaning": "較正状態（CALSTAT列）。", "options": list(CALSTAT_VALUES),
                    "examples": ["RAW", "TASTAR"], "notes": ""},
    "row.tsys_k": {"where": "add_row", "type": "float", "unit": "K", "required": False,
                   "meaning": "Tsys（TSYS列）。不明ならNaN。", "examples": [200.0, np.nan], "notes": ""},
    "row.tau": {"where": "add_row", "type": "float", "unit": "", "required": False,
                "meaning": "天頂光学的厚み τ（TAU列）。", "examples": [0.10, np.nan],
                "notes": "定義を info.shared_meta['TAU_DEF'] に明記推奨。"},
    "row.el_enc_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": False,
                       "meaning": "仰角（ELEVATIO列）。airmass計算に有用。", "examples": [45.0, np.nan], "notes": ""},

    # --- layer1 ---
    "row.ra_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": "ra/dec or glon/glat required",
                   "meaning": "最終絶対RA（RA列）。", "examples": [83.82208], "notes": "基準系は info.radesys/equinox。"},
    "row.dec_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": "ra/dec or glon/glat required",
                    "meaning": "最終絶対Dec（DEC列）。", "examples": [-5.39111], "notes": ""},
    "row.glon_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": "ra/dec or glon/glat required",
                     "meaning": "最終絶対GLON（GLON列）。", "examples": [209.0], "notes": "銀河経度 l。"},
    "row.glat_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": "ra/dec or glon/glat required",
                     "meaning": "最終絶対GLAT（GLAT列）。", "examples": [-19.4], "notes": "銀河緯度 b。"},

    # --- layer2 ---
    "row.srcframe": {"where": "add_row", "type": "str(enum)", "unit": "", "required": True,
                     "meaning": "SRC_LONG/LATの意味（SRCFRAME列）。", "options": list(COORDTYPE_VALUES),
                     "examples": ["RADEC", "GALACTIC"], "notes": "RADECの分点は SRCRDSYS/SRCEQNX（既定はinfo.src_*）。"},
    "row.src_radesys": {"where": "add_row", "type": "str(enum)", "unit": "", "required": "srcframe=='RADEC'",
                        "meaning": "SRCFRAME='RADEC' の基準系（SRCRDSYS列）。", "options": list(RADESYS_VALUES),
                        "examples": ["ICRS", "FK5", "FK4"], "notes": "省略時は info.src_radesys を用いる。"},
    "row.src_equinox": {"where": "add_row", "type": "float", "unit": "year", "required": "srcframe=='RADEC'",
                        "meaning": "SRCFRAME='RADEC' の分点（SRCEQNX列）。", "examples": [2000.0, 1950.0],
                        "notes": "省略時は info.src_equinox を用いる。"},
    "row.src_long_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": True,
                         "meaning": "意図座標Longitude（SRC_LONG列）。", "examples": [83.82208, 110.0], "notes": ""},
    "row.src_lat_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": True,
                        "meaning": "意図座標Latitude（SRC_LAT列）。", "examples": [-5.39111, 20.0], "notes": ""},
    "row.scanframe": {"where": "add_row", "type": "str(enum)", "unit": "", "required": True,
                      "meaning": "SCAN_X/Yの意味（SCANFRAM列）。", "options": list(COORDTYPE_VALUES),
                      "examples": ["AZEL_GEO", "RADEC"], "notes": ""},
    "row.scan_radesys": {"where": "add_row", "type": "str(enum)", "unit": "", "required": "scanframe=='RADEC'",
                         "meaning": "SCANFRAM='RADEC' の基準系（SCANRDSYS列）。", "options": list(RADESYS_VALUES),
                         "examples": ["ICRS", "FK5", "FK4"], "notes": ""},
    "row.scan_equinox": {"where": "add_row", "type": "float", "unit": "year", "required": "scanframe=='RADEC'",
                         "meaning": "SCANFRAM='RADEC' の分点（SCANEQNX列）。", "examples": [2000.0, 1950.0], "notes": ""},
    "row.scan_x_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": True,
                       "meaning": "オフセットX（SCANX列）。", "examples": [0.0, 0.001], "notes": ""},
    "row.scan_y_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": True,
                       "meaning": "オフセットY（SCANY列）。", "examples": [0.0, -0.001], "notes": ""},

    # --- layer3 hardware/correction ---
    "row.az_enc_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": False,
                       "meaning": "方位角（AZIMUTH列）。", "examples": [120.0, np.nan], "notes": ""},
    "row.az_cmd_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": False,
                       "meaning": "方位指令（AZ_CMD列）。", "examples": [120.1, np.nan], "notes": ""},
    "row.el_cmd_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": False,
                       "meaning": "仰角指令（EL_CMD列）。", "examples": [45.1, np.nan], "notes": ""},
    "row.corr_az_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": False,
                        "meaning": "適用補正Az（CORR_AZ列）。", "examples": [0.002, np.nan], "notes": "屈折含むかはREFFR_INC。"},
    "row.corr_el_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": False,
                        "meaning": "適用補正El（CORR_EL列）。", "examples": [-0.001, np.nan], "notes": ""},
    "row.calc_refr_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": False,
                          "meaning": "屈折モデル値（CALC_REFR列）。", "examples": [0.03, np.nan], "notes": ""},

    # --- doppler / LO ---
    "row.v_frame_mps": {"where": "add_row", "type": "float", "unit": "m/s", "required": False,
                        "meaning": "適用速度シフト等（VFRAME列）。", "examples": [-40000.0, np.nan], "notes": ""},
    "row.f_offset_hz": {"where": "add_row", "type": "float", "unit": "Hz", "required": False,
                        "meaning": "周波数オフセット（FOFFSET列）。", "examples": [0.0, 10e6], "notes": ""},

    # --- weather ---
    "row.tamb_c": {"where": "add_row", "type": "float", "unit": "degC", "required": False,
                   "meaning": "外気温（TAMBIENT列）。", "examples": [10.0, np.nan], "notes": ""},
    "row.pressure_hpa": {"where": "add_row", "type": "float", "unit": "hPa", "required": False,
                         "meaning": "気圧（PRESSURE列）。", "examples": [800.0, np.nan], "notes": ""},
    "row.humidity_pct": {"where": "add_row", "type": "float", "unit": "%", "required": False,
                         "meaning": "湿度（HUMIDITY列）。", "examples": [30.0, np.nan], "notes": ""},
    "row.wind_spd_mps": {"where": "add_row", "type": "float", "unit": "m/s", "required": False,
                         "meaning": "風速（WINDSPD列）。", "examples": [5.0, np.nan], "notes": ""},
    "row.wind_dir_deg": {"where": "add_row", "type": "float", "unit": "deg", "required": False,
                         "meaning": "風向（WINDDIR列）。", "examples": [180.0, np.nan], "notes": ""},

    # --- per-row frequency ---
    "row.freq_hz": {"where": "add_row", "type": "Optional[np.ndarray(float)]", "unit": "Hz", "required": "store_freq_column==True",
                    "meaning": "行ごとの周波数ベクトル（FREQ列）。未指定ならuniform axisから自動生成可。",
                    "examples": ["np.ndarray(n_chan)"], "notes": ""},
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
    等間隔周波数WCS（全行共通の場合に推奨）。
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
            raise ValueError("CTYPE1 should be 'FREQ' for 標準的な単一鏡解析/asap interoperability")
        if self.cunit1 != "Hz":
            raise ValueError("CUNIT1 should be 'Hz'")
        if self.refchan < 1:
            raise ValueError("REFCHAN must be >=1")
        if self.restfreq_hz < 0:
            raise ValueError("RESTFREQ must be >=0")
    def velref(self) -> int:
        """Return AIPS/casacore-style VELREF integer.

        - 1: LSRK
        - 2: Barycentric/heliocentric (historically "HEL")
        - 3: Topocentric (historically "OBS")
        - 5: Geocentric

        Add +256 when VELDEF is RADIO.
        """
        base_map = {
            "LSRK": 1,
            "BARYCENT": 2,
            "HELIOCEN": 2,
            "TOPOCENT": 3,
            "GEOCENTR": 5,
        }
        base = base_map.get(self.specsys)
        if base is None:
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
                raise ValueError("SPECSYS='LSRK' requires store_freq_column=True (LSRK axis is generally time/pointing dependent).")


# =============================================================================
# 4) ColumnBuffer（Rowを溜めず列ごとに蓄積）
# =============================================================================

@dataclass
class _ColumnBuffer:
    n_chan: int
    store_freq_column: bool

    # --- SDFITS-ish core columns ---
    TIME: List[float] = field(default_factory=list)      # MJD UTC [day]
    DATEOBS: List[str] = field(default_factory=list)     # ISO UTC (string) - interoperability aid

    SCAN: List[int] = field(default_factory=list)        # scanid
    SUBSCAN: List[int] = field(default_factory=list)
    INTGRP: List[int] = field(default_factory=list)

    OBJECT: List[str] = field(default_factory=list)      # per-row object name
    OBSMODE: List[str] = field(default_factory=list)

    BEAMNO: List[int] = field(default_factory=list)
    POLNO: List[int] = field(default_factory=list)

    EXPOSURE: List[float] = field(default_factory=list)  # [s]
    CALSTAT: List[str] = field(default_factory=list)
    TSYS: List[float] = field(default_factory=list)      # [K]
    TAU: List[float] = field(default_factory=list)       # dimensionless (zenith)

    FLAGROW: List[int] = field(default_factory=list)

    # --- sky coordinates (Layer1) ---
    RA: List[float] = field(default_factory=list)        # deg
    DEC: List[float] = field(default_factory=list)       # deg
    GLON: List[float] = field(default_factory=list)      # deg
    GLAT: List[float] = field(default_factory=list)      # deg

    # --- pointing/hardware (common single-dish columns) ---
    AZIMUTH: List[float] = field(default_factory=list)   # deg
    ELEVATIO: List[float] = field(default_factory=list)  # deg (common name in many SDFITS variants)

    # --- extras (Layer2/Layer3) keep all features ---
    SRCFRAME: List[str] = field(default_factory=list)
    SRCRDSYS: List[str] = field(default_factory=list)    # SRC frame RADESYS for SRCFRAME='RADEC' (per-row)
    SRCEQNX: List[float] = field(default_factory=list)   # SRC frame EQUINOX for SRCFRAME='RADEC' (per-row)
    SRC_LONG: List[float] = field(default_factory=list)
    SRC_LAT: List[float] = field(default_factory=list)

    SCANFRAM: List[str] = field(default_factory=list)    # 互換性のため 8文字にしない（列名はfitsでOK）
    SCANRDSYS: List[str] = field(default_factory=list)   # SCAN frame RADESYS for SCANFRAM='RADEC'
    SCANEQNX: List[float] = field(default_factory=list)  # SCAN frame EQUINOX for SCANFRAM='RADEC'
    SCANX: List[float] = field(default_factory=list)
    SCANY: List[float] = field(default_factory=list)

    AZ_CMD: List[float] = field(default_factory=list)
    EL_CMD: List[float] = field(default_factory=list)
    CORR_AZ: List[float] = field(default_factory=list)
    CORR_EL: List[float] = field(default_factory=list)
    CALC_REFR: List[float] = field(default_factory=list)

    VFRAME: List[float] = field(default_factory=list)    # m/s
    FOFFSET: List[float] = field(default_factory=list)   # Hz

    TAMBIENT: List[float] = field(default_factory=list)
    PRESSURE: List[float] = field(default_factory=list)
    HUMIDITY: List[float] = field(default_factory=list)
    WINDSPD: List[float] = field(default_factory=list)
    WINDDIR: List[float] = field(default_factory=list)

    # vectors
    _SPEC: List[np.ndarray] = field(default_factory=list)  # (n_chan,) float32 or (1,)
    _FLAG: List[np.ndarray] = field(default_factory=list)  # (n_chan,) bool
    _FREQ: List[np.ndarray] = field(default_factory=list)  # (n_chan,) float64

    def __len__(self) -> int:
        return len(self.TIME)

    def clear(self) -> None:
        for v in self.__dict__.values():
            if isinstance(v, list):
                v.clear()

    def append_spectrum(self, spec: Union[np.ndarray, float], flag: Optional[np.ndarray]) -> None:
        if self.n_chan == 1:
            x = float(np.asarray(spec).reshape(-1)[0]) if isinstance(spec, np.ndarray) else float(spec)
            self._SPEC.append(np.array([x], dtype=np.float32))
            if flag is None:
                self._FLAG.append(np.array([False], dtype=bool))
            else:
                f = np.asarray(flag, dtype=bool).reshape(-1)
                if f.size != 1:
                    raise ValueError("FLAG must have length 1 for n_chan=1")
                self._FLAG.append(np.array([bool(f[0])], dtype=bool))
        else:
            x = np.asarray(spec, dtype=np.float32).reshape(-1)
            if x.size != self.n_chan:
                raise ValueError(f"SPECTRUM must have length {self.n_chan}, got {x.size}")
            self._SPEC.append(x)
            if flag is None:
                self._FLAG.append(np.zeros(self.n_chan, dtype=bool))
            else:
                f = np.asarray(flag, dtype=bool).reshape(-1)
                if f.size != self.n_chan:
                    raise ValueError(f"FLAG must have length {self.n_chan}, got {f.size}")
                self._FLAG.append(f)

    def append_freq(self, freq_hz: np.ndarray) -> None:
        if not self.store_freq_column:
            return
        f = np.asarray(freq_hz, dtype=np.float64).reshape(-1)
        need = 1 if self.n_chan == 1 else self.n_chan
        if f.size != need:
            raise ValueError(f"FREQ must have length {need}, got {f.size}")
        self._FREQ.append(f)


# =============================================================================
# 5) Writer本体（互換性最優先）
# =============================================================================

class SDRadioSpectralSDFITSWriter:
    """
    互換性最優先の SDFITS-like FITS Writer

    - 主要カラム名を SDFITS/単一鏡で一般的なものに寄せる：
        TIME, DATE-OBS(列としてDATEOBS), SCAN, OBJECT, OBSMODE,
        RA, DEC, AZIMUTH, ELEVATIO, EXPOSURE, TSYS, SPECTRUM, FLAG, FLAGROW など
    - 追加情報（意図座標/オフセット/補正等）は機能省略せず保持。

    注：FITS標準には“単一鏡SDFITSの唯一の正解”が無いので、互換保険として
        * スペクトルを SPECTRUM と DATA の両方に同内容で出す（既定）
      を採用します。
    """

    # fixed-length strings (generous to avoid truncation)
    STR_OBJECT  = 32
    STR_OBSMODE = 24
    STR_CALSTAT = 16
    STR_FRAME   = 16

    def __init__(
        self,
        *,
        n_chan: int,
        site: Site,
        info: DatasetInfo,
        store_freq_column: bool = False,
        duplicate_data_columns: bool = True,
        chunk_size: Optional[int] = None,
        out_basename: Optional[str] = None,
    ):
        self.n_chan = int(n_chan)
        if self.n_chan < 1:
            raise ValueError("n_chan must be >= 1")

        self.site = site
        self.store_freq_column = bool(store_freq_column)
        self.duplicate_data_columns = bool(duplicate_data_columns)

        info.validate(self.store_freq_column)
        self.info = info

        self.chunk_size = int(chunk_size) if chunk_size is not None else None
        if self.chunk_size is not None and self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive int or None")
        self.out_basename = out_basename
        if self.chunk_size is not None and not self.out_basename:
            raise ValueError("parts mode requires out_basename, e.g. 'otf_run1'.")

        self._buf = _ColumnBuffer(n_chan=self.n_chan, store_freq_column=self.store_freq_column)
        self._part_index = 0
        self._manifest: Dict[str, Any] = {
            "format": "sd-radio-spectral-fits parts manifest",
            "n_chan": self.n_chan,
            "store_freq_column": self.store_freq_column,
            "duplicate_data_columns": self.duplicate_data_columns,
            "parts": [],
        }

    # -------------------- schema helper --------------------
    @staticmethod
    def print_schema() -> None:
        keys = sorted(SCHEMA.keys())
        for k in keys:
            d = SCHEMA[k]
            print(f"\n[{k}]")
            for kk in ("where", "type", "unit", "required", "meaning"):
                if kk in d:
                    print(f"  {kk}: {d[kk]}")
            if "options" in d:
                print(f"  options: {d['options']}")
            if "examples" in d:
                print(f"  examples: {d['examples']}")
            if "notes" in d and d["notes"]:
                print(f"  notes: {d['notes']}")

    # -------------------- coordinate helpers --------------------
    def _radec_frame_for_columns(self):
        rs = self.info.radesys
        if rs == "ICRS":
            return "icrs"
        if rs == "FK5":
            return FK5(equinox=Time(self.info.equinox, format="jyear"))
        if rs in ("FK4", "FK4-NO-E"):
            return FK4(equinox=Time(self.info.equinox, format="jyear"))
        if rs == "GALACTIC":
            # RA/DEC列にGALACTICを採る運用は基本想定しないが、変換先として許容
            return Galactic()
        return "icrs"

    def _skycoord_for_rvcor(self, ra_deg: float, dec_deg: float) -> SkyCoord:
        rs = self.info.radesys
        if rs == "ICRS":
            return SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs")
        if rs == "FK5":
            eq = Time(self.info.equinox, format="jyear")
            return SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame=FK5(equinox=eq))
        if rs in ("FK4", "FK4-NO-E"):
            eq = Time(self.info.equinox, format="jyear")
            return SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame=FK4(equinox=eq))
        return SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs")

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

        # If we have RA/DEC -> compute Galactic
        if have_radec and (not have_gal):
            sc = SkyCoord(ra=float(ra_deg)*u.deg, dec=float(dec_deg)*u.deg, frame=self._radec_frame_for_columns())
            gal = sc.galactic
            glon_deg = float(gal.l.to_value(u.deg))
            glat_deg = float(gal.b.to_value(u.deg))
            return float(ra_deg), float(dec_deg), glon_deg, glat_deg

        # If we have Galactic -> compute RA/DEC in info.radesys/equinox
        if have_gal and (not have_radec):
            gal = SkyCoord(l=float(glon_deg)*u.deg, b=float(glat_deg)*u.deg, frame=Galactic())
            radec = gal.transform_to(self._radec_frame_for_columns())
            # transform_to may return Galactic() when info.radesys == 'GALACTIC'
            if hasattr(radec, "ra") and hasattr(radec, "dec"):
                ra_deg = float(radec.ra.to_value(u.deg))
                dec_deg = float(radec.dec.to_value(u.deg))
            else:
                # fallback: if columns are Galactic (non-typical), map l,b to RA,DEC fields anyway
                ra_deg = float(gal.l.to_value(u.deg))
                dec_deg = float(gal.b.to_value(u.deg))
            return ra_deg, dec_deg, float(glon_deg), float(glat_deg)

        # both provided -> keep as-is
        return float(ra_deg), float(dec_deg), float(glon_deg), float(glat_deg)

    # -------------------- frequency helpers --------------------
    def _freq_from_uniform_axis(self) -> np.ndarray:
        ax = self.info.spectral_axis
        if ax is None:
            raise ValueError("No uniform spectral_axis to generate freq.")
        pix = (np.arange(self.n_chan, dtype=np.float64) + 1.0)  # 1-based
        return ax.crval1_hz + (pix - float(ax.crpix1)) * ax.cdelt1_hz

    def _topo_freq_to_lsrk_freq(self, freq_topo_hz: np.ndarray, time_mjd: float, ra_deg: float, dec_deg: float) -> np.ndarray:
        ax = self.info.spectral_axis
        if ax is None:
            raise ValueError("LSRK conversion requires info.spectral_axis (RESTFREQ needed).")
        if ax.veldef != "RADIO":
            raise ValueError("LSRK conversion in this writer assumes VELDEF='RADIO'.")
        if ax.restfreq_hz <= 0:
            raise ValueError("LSRK conversion requires RESTFREQ > 0.")

        rest = ax.restfreq_hz * u.Hz
        f_topo = np.asarray(freq_topo_hz, dtype=np.float64).reshape(-1) * u.Hz

        t = Time(float(time_mjd), format="mjd", scale="utc")
        loc = self.site.to_earthlocation()
        sc = self._skycoord_for_rvcor(ra_deg, dec_deg)

        # NOTE: astropy's radial_velocity_correction supports barycentric/heliocentric kinds.
        # Here we keep frequency vector as stored, and VLSRK should be carried separately if needed by upstream.
        v_topo = f_topo.to(u.m/u.s, equivalencies=u.doppler_radio(rest))
        v_corr = sc.radial_velocity_correction(kind="barycentric", obstime=t, location=loc)
        v_lsrk_like = v_topo + v_corr
        f_lsrk_like = v_lsrk_like.to(u.Hz, equivalencies=u.doppler_radio(rest))
        return np.asarray(f_lsrk_like.to_value(u.Hz), dtype=np.float64)

    # -------------------- validation helpers --------------------
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
            raise ValueError(f"{name} too long (len={len(val)} > {maxlen}). Use short code. val={val}")
        return val

    # -------------------- public: add_row --------------------
    def add_row(
        self,
        *,
        # core
        time_mjd: float,     # MJD UTC (days)
        scanid: int,         # SCAN
        subscan: int,        # SUBSCAN
        intgrp: int,         # INTGRP
        obsmode: str,        # OBSMODE enum
        data: Union[np.ndarray, float],  # spectrum
        exposure_s: float,   # EXPOSURE

        # per-row object
        object_name: Optional[str] = None,

        # QA
        flagrow: int = 0,
        flag: Optional[np.ndarray] = None,

        # identifiers
        beamno: int = 0,
        polno: int = 0,

        # calibration
        calstat: str = "RAW",
        tsys_k: float = np.nan,
        tau: float = np.nan,

        # layer1 (either RA/DEC or GLON/GLAT required; both will be stored)
        ra_deg: float = np.nan,
        dec_deg: float = np.nan,
        glon_deg: float = np.nan,
        glat_deg: float = np.nan,

        # layer2 (required)
        srcframe: str = "RADEC",
        src_radesys: Optional[str] = None,
        src_equinox: Optional[float] = None,
        src_long_deg: float = np.nan,
        src_lat_deg: float = np.nan,
        scanframe: str = "RADEC",
        scan_radesys: str = "ICRS",
        scan_equinox: float = 2000.0,
        scan_x_deg: float = 0.0,
        scan_y_deg: float = 0.0,

        # layer3
        az_enc_deg: float = np.nan,
        el_enc_deg: float = np.nan,
        az_cmd_deg: float = np.nan,
        el_cmd_deg: float = np.nan,
        corr_az_deg: float = np.nan,
        corr_el_deg: float = np.nan,
        calc_refr_deg: float = np.nan,

        # doppler/LO
        v_frame_mps: float = np.nan,
        f_offset_hz: float = 0.0,

        # weather
        tamb_c: float = np.nan,
        pressure_hpa: float = np.nan,
        humidity_pct: float = np.nan,
        wind_spd_mps: float = np.nan,
        wind_dir_deg: float = np.nan,

        # per-row freq
        freq_hz: Optional[np.ndarray] = None,
    ) -> None:
        # basic checks
        if not np.isfinite(time_mjd):
            raise ValueError("time_mjd must be finite float (MJD UTC days)")
        if exposure_s <= 0:
            raise ValueError("exposure_s must be > 0")
        if scanid < 0 or subscan < 0 or intgrp < 0:
            raise ValueError("scanid/subscan/intgrp must be >= 0")
        if flagrow not in (0, 1):
            raise ValueError("flagrow must be 0 or 1")

        # enums
        obsmode = self._validate_enum("obsmode", obsmode, OBSMODE_VALUES)
        calstat = self._validate_enum("calstat", calstat, CALSTAT_VALUES)
        srcframe = self._validate_enum("srcframe", srcframe, COORDTYPE_VALUES)
        scanframe = self._validate_enum("scanframe", scanframe, COORDTYPE_VALUES)

        # strings length
        obsmode = self._validate_str("obsmode", obsmode, self.STR_OBSMODE)
        calstat = self._validate_str("calstat", calstat, self.STR_CALSTAT)
        srcframe = self._validate_str("srcframe", srcframe, self.STR_FRAME)
        scanframe = self._validate_str("scanframe", scanframe, self.STR_FRAME)

        obj = self.info.object_name if object_name is None else object_name
        obj = self._validate_str("object_name", obj, self.STR_OBJECT)

        # require finite key coordinates for 標準的な単一鏡解析 usability (fill missing side)
        ra_deg, dec_deg, glon_deg, glat_deg = self._fill_radec_glon_glat(
            ra_deg=ra_deg, dec_deg=dec_deg, glon_deg=glon_deg, glat_deg=glat_deg
        )
        for nm, vv in (("ra_deg", ra_deg), ("dec_deg", dec_deg), ("glon_deg", glon_deg), ("glat_deg", glat_deg),
                       ("src_long_deg", src_long_deg), ("src_lat_deg", src_lat_deg)):
            if not np.isfinite(vv):
                raise ValueError(f"{nm} must be finite float (deg). Provide actual values, not NaN.")

        # SRCFRAME='RADEC' needs per-row RADESYS/EQUINOX (stored as columns)
        # 省略時は info.src_radesys/info.src_equinox を用いる
        if srcframe == "RADEC":
            if src_radesys is None:
                src_radesys = self.info.src_radesys
            if src_equinox is None:
                src_equinox = float(self.info.src_equinox)
            src_radesys = self._validate_enum("src_radesys", str(src_radesys), RADESYS_VALUES)
            src_radesys = self._validate_str("src_radesys", str(src_radesys), self.STR_FRAME)
            if (src_equinox is None) or (not np.isfinite(float(src_equinox))):
                raise ValueError("src_equinox must be finite float (year) when srcframe=='RADEC'")
            src_equinox = float(src_equinox)
        else:
            # non-RADEC SRCFRAME: keep placeholders (stored, but not required/used)
            src_radesys = self._validate_str("src_radesys", "" if src_radesys is None else str(src_radesys), self.STR_FRAME)
            src_equinox = float("nan") if src_equinox is None else float(src_equinox)

        # SCANFRAM='RADEC' needs per-row RADESYS/EQUINOX (stored as columns)
        if scanframe == "RADEC":
            scan_radesys = self._validate_enum("scan_radesys", scan_radesys, RADESYS_VALUES)
            scan_radesys = self._validate_str("scan_radesys", scan_radesys, self.STR_FRAME)
            if not np.isfinite(scan_equinox):
                raise ValueError("scan_equinox must be finite float (year) when scanframe=='RADEC'")
        else:
            scan_radesys = self._validate_str("scan_radesys", str(scan_radesys), self.STR_FRAME)

        # buffer append
        b = self._buf
        b.TIME.append(float(time_mjd))
        b.DATEOBS.append(Time(time_mjd, format="mjd", scale="utc").isot)

        b.SCAN.append(int(scanid))
        b.SUBSCAN.append(int(subscan))
        b.INTGRP.append(int(intgrp))

        b.OBJECT.append(obj)
        b.OBSMODE.append(obsmode)

        b.BEAMNO.append(int(beamno))
        b.POLNO.append(int(polno))

        b.EXPOSURE.append(float(exposure_s))
        b.CALSTAT.append(calstat)
        b.TSYS.append(float(tsys_k))
        b.TAU.append(float(tau))

        b.FLAGROW.append(int(flagrow))

        b.RA.append(float(ra_deg))
        b.DEC.append(float(dec_deg))
        b.GLON.append(float(glon_deg))
        b.GLAT.append(float(glat_deg))

        b.AZIMUTH.append(float(az_enc_deg))
        b.ELEVATIO.append(float(el_enc_deg))

        b.SRCFRAME.append(srcframe)
        b.SRCRDSYS.append(str(src_radesys))
        b.SRCEQNX.append(float(src_equinox))
        b.SRC_LONG.append(float(src_long_deg))
        b.SRC_LAT.append(float(src_lat_deg))

        b.SCANFRAM.append(scanframe)
        b.SCANRDSYS.append(str(scan_radesys))
        b.SCANEQNX.append(float(scan_equinox))
        b.SCANX.append(float(scan_x_deg))
        b.SCANY.append(float(scan_y_deg))

        b.AZ_CMD.append(float(az_cmd_deg))
        b.EL_CMD.append(float(el_cmd_deg))
        b.CORR_AZ.append(float(corr_az_deg))
        b.CORR_EL.append(float(corr_el_deg))
        b.CALC_REFR.append(float(calc_refr_deg))

        b.VFRAME.append(float(v_frame_mps))
        b.FOFFSET.append(float(f_offset_hz))

        b.TAMBIENT.append(float(tamb_c))
        b.PRESSURE.append(float(pressure_hpa))
        b.HUMIDITY.append(float(humidity_pct))
        b.WINDSPD.append(float(wind_spd_mps))
        b.WINDDIR.append(float(wind_dir_deg))

        b.append_spectrum(data, flag)

        if self.store_freq_column:
            if freq_hz is None:
                if self.info.spectral_axis is None:
                    raise ValueError("store_freq_column=True and freq_hz=None requires info.spectral_axis for auto generation.")
                freq_hz = np.array([self.info.spectral_axis.crval1_hz], dtype=np.float64) if self.n_chan == 1 else self._freq_from_uniform_axis()

            # --- Case C: store per-row "LSRK-like" frequency vector while keeping CTYPE1='FREQ' ---
            if (self.info.spectral_axis is not None) and (self.info.spectral_axis.specsys == "LSRK"):
                freq_hz = self._topo_freq_to_lsrk_freq(freq_hz, time_mjd=float(time_mjd), ra_deg=float(ra_deg), dec_deg=float(dec_deg))

            b.append_freq(freq_hz)

        # parts flush
        if self.chunk_size is not None and len(b) >= self.chunk_size:
            self._flush_part()

    # -------------------- FITS construction --------------------
    def _build_primary_hdu(self, mjdstart: float, mjdend: float) -> fits.PrimaryHDU:
        hdr = fits.Header()
        hdr["SIMPLE"] = True
        hdr["BITPIX"] = 8
        hdr["NAXIS"] = 0
        hdr["EXTEND"] = True

        # Primary-level identifiers (common in many SDFITS)
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
        return fits.PrimaryHDU(header=hdr)

    def _build_table_hdu(self) -> fits.BinTableHDU:
        b = self._buf
        if len(b) == 0:
            raise RuntimeError("No data to write.")

        nrow = len(b)
        nchan = 1 if self.n_chan == 1 else self.n_chan

        spec2 = np.vstack(b._SPEC).astype(np.float32)  # (nrow,nchan)
        flag2 = np.vstack(b._FLAG).astype(np.bool_)    # (nrow,nchan)

        cols: List[fits.Column] = []

        # --- core time/ids ---
        cols.append(fits.Column(name="TIME",    format="D",  unit="d",   array=np.array(b.TIME, dtype=np.float64)))  # MJD
        cols.append(fits.Column(name="DATEOBS", format="26A",          array=np.array(b.DATEOBS, dtype="S26")))    # ISO UTC string (aid)
        cols.append(fits.Column(name="SCAN",    format="J",            array=np.array(b.SCAN, dtype=np.int32)))
        cols.append(fits.Column(name="SUBSCAN", format="J",            array=np.array(b.SUBSCAN, dtype=np.int32)))
        cols.append(fits.Column(name="INTGRP",  format="J",            array=np.array(b.INTGRP, dtype=np.int32)))

        # --- object/mode ---
        cols.append(fits.Column(name="OBJECT",  format=f"{self.STR_OBJECT}A",
                                array=np.array(b.OBJECT, dtype=f"S{self.STR_OBJECT}")))
        cols.append(fits.Column(name="OBSMODE", format=f"{self.STR_OBSMODE}A",
                                array=np.array(b.OBSMODE, dtype=f"S{self.STR_OBSMODE}")))

        cols.append(fits.Column(name="BEAMNO",  format="J",            array=np.array(b.BEAMNO, dtype=np.int32)))
        cols.append(fits.Column(name="POLNO",   format="J",            array=np.array(b.POLNO, dtype=np.int32)))

        # --- integration/cal ---
        cols.append(fits.Column(name="EXPOSURE", format="E", unit="s", array=np.array(b.EXPOSURE, dtype=np.float32)))
        cols.append(fits.Column(name="CALSTAT",  format=f"{self.STR_CALSTAT}A",
                                array=np.array(b.CALSTAT, dtype=f"S{self.STR_CALSTAT}")))
        cols.append(fits.Column(name="TSYS",     format="E", unit="K", array=np.array(b.TSYS, dtype=np.float32)))
        cols.append(fits.Column(name="TAU",      format="E", unit="",  array=np.array(b.TAU, dtype=np.float32)))

        cols.append(fits.Column(name="FLAGROW",  format="J",           array=np.array(b.FLAGROW, dtype=np.int32)))

        # --- coordinates / pointing ---
        cols.append(fits.Column(name="RA",       format="D", unit="deg", array=np.array(b.RA, dtype=np.float64)))
        cols.append(fits.Column(name="DEC",      format="D", unit="deg", array=np.array(b.DEC, dtype=np.float64)))
        cols.append(fits.Column(name="GLON",     format="D", unit="deg", array=np.array(b.GLON, dtype=np.float64)))
        cols.append(fits.Column(name="GLAT",     format="D", unit="deg", array=np.array(b.GLAT, dtype=np.float64)))

        cols.append(fits.Column(name="AZIMUTH",  format="D", unit="deg", array=np.array(b.AZIMUTH, dtype=np.float64)))
        cols.append(fits.Column(name="ELEVATIO", format="D", unit="deg", array=np.array(b.ELEVATIO, dtype=np.float64)))

        # --- spectrum + flags ---
        cols.append(fits.Column(name="SPECTRUM", format=f"{nchan}E", unit="K", array=spec2))
        cols.append(fits.Column(name="FLAG",     format=f"{nchan}L",          array=flag2))

        # compatibility insurance: duplicate spectrum to 'DATA' if requested
        if self.duplicate_data_columns:
            cols.append(fits.Column(name="DATA", format=f"{nchan}E", unit="K", array=spec2))

        # --- keep all extra features (Layer2/Layer3) as columns ---
        cols.append(fits.Column(name="SRCFRAME", format=f"{self.STR_FRAME}A",
                                array=np.array(b.SRCFRAME, dtype=f"S{self.STR_FRAME}")))
        cols.append(fits.Column(name="SRCRDSYS", format=f"{self.STR_FRAME}A",
                                array=np.array(b.SRCRDSYS, dtype=f"S{self.STR_FRAME}")))
        cols.append(fits.Column(name="SRCEQNX",  format="D", unit="year",
                                array=np.array(b.SRCEQNX, dtype=np.float64)))
        cols.append(fits.Column(name="SRC_LONG", format="D", unit="deg", array=np.array(b.SRC_LONG, dtype=np.float64)))
        cols.append(fits.Column(name="SRC_LAT",  format="D", unit="deg", array=np.array(b.SRC_LAT, dtype=np.float64)))

        cols.append(fits.Column(name="SCANFRAM", format=f"{self.STR_FRAME}A",
                                array=np.array(b.SCANFRAM, dtype=f"S{self.STR_FRAME}")))
        cols.append(fits.Column(name="SCANRDSYS", format=f"{self.STR_FRAME}A",
                                array=np.array(b.SCANRDSYS, dtype=f"S{self.STR_FRAME}")))
        cols.append(fits.Column(name="SCANEQNX", format="D", unit="year", array=np.array(b.SCANEQNX, dtype=np.float64)))

        cols.append(fits.Column(name="SCANX",    format="D", unit="deg", array=np.array(b.SCANX, dtype=np.float64)))
        cols.append(fits.Column(name="SCANY",    format="D", unit="deg", array=np.array(b.SCANY, dtype=np.float64)))

        cols.append(fits.Column(name="AZ_CMD",   format="D", unit="deg", array=np.array(b.AZ_CMD, dtype=np.float64)))
        cols.append(fits.Column(name="EL_CMD",   format="D", unit="deg", array=np.array(b.EL_CMD, dtype=np.float64)))
        cols.append(fits.Column(name="CORR_AZ",  format="D", unit="deg", array=np.array(b.CORR_AZ, dtype=np.float64)))
        cols.append(fits.Column(name="CORR_EL",  format="D", unit="deg", array=np.array(b.CORR_EL, dtype=np.float64)))
        cols.append(fits.Column(name="CALC_REFR",format="D", unit="deg", array=np.array(b.CALC_REFR, dtype=np.float64)))

        cols.append(fits.Column(name="VFRAME",   format="E", unit="m/s", array=np.array(b.VFRAME, dtype=np.float32)))
        cols.append(fits.Column(name="FOFFSET",  format="E", unit="Hz",  array=np.array(b.FOFFSET, dtype=np.float32)))

        cols.append(fits.Column(name="TAMBIENT", format="E", unit="degC", array=np.array(b.TAMBIENT, dtype=np.float32)))
        cols.append(fits.Column(name="PRESSURE", format="E", unit="hPa",  array=np.array(b.PRESSURE, dtype=np.float32)))
        cols.append(fits.Column(name="HUMIDITY", format="E", unit="percent", array=np.array(b.HUMIDITY, dtype=np.float32)))
        cols.append(fits.Column(name="WINDSPD",  format="E", unit="m/s",  array=np.array(b.WINDSPD, dtype=np.float32)))
        cols.append(fits.Column(name="WINDDIR",  format="E", unit="deg",  array=np.array(b.WINDDIR, dtype=np.float32)))

        if self.store_freq_column:
            freq2 = np.vstack(b._FREQ).astype(np.float64)  # (nrow,nchan)
            cols.append(fits.Column(name="FREQ", format=f"{nchan}D", unit="Hz", array=freq2))

        hdu = fits.BinTableHDU.from_columns(cols)
        hdu.name = "SINGLE DISH"

        # -------- extension header (SDFITS-like metadata) --------
        eh = hdu.header
        eh["EXTNAME"]  = ("SINGLE DISH", "Single dish table (SDFITS-like)")
        eh["TELESCOP"] = (self.info.telescope, "Telescope name")
        eh["OBSERVER"] = (self.info.observer,  "Observer")
        eh["PROJID"]   = (self.info.project,   "Project/Proposal ID")
        eh["OBJECT"]   = (self.info.object_name, "Representative target name")
        eh["SWNAME"]   = (SWNAME, "Software name")
        eh["SWVER"]    = (__version__, "Software version")

        mjd0 = float(min(b.TIME))
        mjd1 = float(max(b.TIME))
        eh["TIMESYS"]  = ("UTC", "Time system for TIME column")
        eh["MJDSTART"] = (mjd0, "Start MJD(UTC)")
        eh["MJDEND"]   = (mjd1, "End MJD(UTC)")

        # coordinate frames
        eh["RADESYS"]  = (self.info.radesys, "RA/Dec frame")
        eh["EQUINOX"]  = (float(self.info.equinox), "Equinox (FK4/FK5)")
        eh["HIERARCH SRC_RADESYS"] = (self.info.src_radesys, "SRC RA/Dec frame")
        eh["HIERARCH SRC_EQUINOX"] = (float(self.info.src_equinox), "SRC equinox (FK4/FK5)")
        
        # site (OBSGEO is widely used)
        x, y, z = self.site.obsgeo_xyz_m()
        eh["SITELAT"]  = (float(self.site.lat_deg), "Site latitude (deg)")
        eh["SITELONG"] = (float(self.site.lon_deg), "Site longitude (deg, east+)")
        eh["SITEELEV"] = (float(self.site.elev_m),  "Site elevation (m)")
        eh["OBSGEO-X"] = (float(x), "Observatory geocenter X (m)")
        eh["OBSGEO-Y"] = (float(y), "Observatory geocenter Y (m)")
        eh["OBSGEO-Z"] = (float(z), "Observatory geocenter Z (m)")

        # beam (optional)
        if np.isfinite(self.info.bmaj_deg):
            eh["BMAJ"] = (float(self.info.bmaj_deg), "Beam major axis (deg)")
        if np.isfinite(self.info.bmin_deg):
            eh["BMIN"] = (float(self.info.bmin_deg), "Beam minor axis (deg)")
        if np.isfinite(self.info.bpa_deg):
            eh["BPA"]  = (float(self.info.bpa_deg), "Beam position angle (deg)")

        # efficiencies (optional)
        if np.isfinite(self.info.eff.beameff):
            eh["BEAMEFF"] = (float(self.info.eff.beameff), "Main-beam efficiency eta_mb")
        if np.isfinite(self.info.eff.apereff):
            eh["APEREFF"] = (float(self.info.eff.apereff), "Aperture efficiency eta_a")
        if np.isfinite(self.info.eff.mooneff):
            eh["MOONEFF"] = (float(self.info.eff.mooneff), "Moon efficiency")
        if np.isfinite(self.info.eff.suneff):
            eh["SUNEFF"]  = (float(self.info.eff.suneff), "Sun efficiency")
        eh["EFFSTAT"] = (str(self.info.eff.effstat), "Efficiency status")
        
        # interpretation flags
        eh["HIERARCH REFR_INC"] = (bool(self.info.refr_included_in_corr), "CORR_* includes refraction term")
        eh["HIERARCH DOPPLER"]  = (bool(self.info.doppler_tracking_applied), "Online Doppler tracking applied")

        # spectral WCS (uniform axis)
        if self.info.spectral_axis is not None:
            ax = self.info.spectral_axis

            # --- Astropy WCS friendliness (minimal, no 標準的な単一鏡解析 downside) ---
            eh["WCSAXES"] = (1, "Number of WCS axes (spectral only)")
            eh["PC1_1"]   = (1.0, "Coordinate transformation matrix element")
            eh["WCSNAME"] = ("SPECTRAL", "WCS name")

            # widely recognized keywords
            eh["CTYPE1"]   = (ax.ctype1, "Spectral axis type (recommended: FREQ)")
            eh["CUNIT1"]   = (ax.cunit1, "Spectral axis unit (recommended: Hz)")
            eh["CRVAL1"]   = (float(ax.crval1_hz), "Reference frequency at CRPIX1 (Hz)")
            eh["CDELT1"]   = (float(ax.cdelt1_hz), "Channel spacing (Hz)")
            eh["CRPIX1"]   = (float(ax.crpix1), "Reference pixel (1-based)")
            # both spellings are seen in the wild; include both for compatibility
            eh["RESTFREQ"] = (float(ax.restfreq_hz), "Rest frequency (Hz)")
            eh["RESTFRQ"]  = (float(ax.restfreq_hz), "Rest frequency (Hz) (alias)")
            eh["SPECSYS"]  = (str(ax.specsys), "Spectral ref frame")
            eh["SSYSOBS"] = (str(ax.ssysobs), "Observer frame for SSYSOBS (often TOPOCENT)")
            eh["VELREF"]  = (int(ax.velref()), "AIPS/casacore VELREF (+256 for RADIO)")
            eh["VELDEF"]   = (str(ax.veldef),  "Velocity definition")
            eh["REFCHAN"]  = (int(ax.refchan), "Reference channel index (1-based)")

        # vector dimensions (TDIMi)
        for i in range(1, eh["TFIELDS"] + 1):
            ttype = (eh.get(f"TTYPE{i}", "") or "").strip()
            if ttype in ("SPECTRUM", "FLAG", "DATA", "FREQ"):
                eh[f"TDIM{i}"] = (f"({nchan})", "Vector length per row")

        # shared meta
        for k, v in self.info.shared_meta.items():
            eh[f"HIERARCH {k}"] = v

        return hdu

    # -------------------- write modes --------------------
    def write(self, filename: str, overwrite: bool = True) -> None:
        b = self._buf
        if len(b) == 0:
            raise RuntimeError("No data to write.")
        mjd0 = float(min(b.TIME))
        mjd1 = float(max(b.TIME))
        pri = self._build_primary_hdu(mjd0, mjd1)
        tab = self._build_table_hdu()
        fits.HDUList([pri, tab]).writeto(filename, overwrite=overwrite)
        print(f"[{SWNAME} {__version__}] Saved FITS: {filename} (rows={len(b)}, n_chan={self.n_chan})")

    def _flush_part(self) -> None:
        if self.chunk_size is None:
            return
        b = self._buf
        if len(b) == 0:
            return
        self._part_index += 1
        partname = f"{self.out_basename}_part{self._part_index:04d}.fits"

        mjd0 = float(min(b.TIME))
        mjd1 = float(max(b.TIME))
        pri = self._build_primary_hdu(mjd0, mjd1)
        tab = self._build_table_hdu()
        fits.HDUList([pri, tab]).writeto(partname, overwrite=True)

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

    # convenience
    @staticmethod
    def airmass_from_el_deg(el_deg: float) -> float:
        """簡易 airmass：X ≈ 1/sin(El)"""
        el = math.radians(el_deg)
        s = math.sin(el)
        if s <= 0:
            return float("nan")
        return 1.0 / s


# =============================================================================
# 6) 例（必要なら削除）
# =============================================================================
def _example_small():
    site = Site(lat_deg=35.9400, lon_deg=138.4720, elev_m=1350.0)

    axis = SpectralAxisUniform(
        crval1_hz=110.201354e9,
        cdelt1_hz=61e3,
        crpix1=1.0,
        restfreq_hz=110.201354e9,
        specsys="LSRK",
        veldef="RADIO",
        refchan=1,
    )

    info = DatasetInfo(
        telescope="MyDish",
        observer="Observer",
        project="TEST",
        object_name="OrionKL",
        radesys="ICRS",
        equinox=2000.0,
        src_radesys="ICRS",
        src_equinox=2000.0,
        eff=Efficiency(beameff=np.nan, apereff=np.nan, effstat="UNKNOWN"),
        refr_included_in_corr=True,
        doppler_tracking_applied=False,
        spectral_axis=axis,
        shared_meta={
            "CALMETHOD": "CHOPPER_WHEEL",
            "TAU_DEF": "zenith tau at observing frequency (dimensionless)",
            "RECEIVER": "RX_A",
            "BACKEND": "XFFTS",
        },
    )

    w = SDRadioSpectralSDFITSWriter(
        n_chan=1024,
        site=site,
        info=info,
        store_freq_column=True,
        duplicate_data_columns=True,  # 互換保険
    )

    mjd0 = 60200.0
    for i in range(10):
        spec = np.random.normal(0, 1, 1024).astype(np.float32)
        w.add_row(
            time_mjd=mjd0 + i * (0.1 / 86400.0),
            scanid=1, subscan=10, intgrp=100,
            obsmode="ON",
            data=spec,
            exposure_s=0.1,

            calstat="TASTAR",
            tsys_k=200.0,
            tau=0.10,

            ra_deg=83.82208, dec_deg=-5.39111,
            srcframe="RADEC", src_long_deg=83.82208, src_lat_deg=-5.39111,
            scanframe="AZEL_GEO", scan_x_deg=0.0, scan_y_deg=0.0,

            az_enc_deg=120.0,
            el_enc_deg=45.0,
        )

    w.write("example_casa_sdfits.fits", overwrite=True)


def _example_parts():
    site = Site(35.9400, 138.4720, 1350.0)
    axis = SpectralAxisUniform(
        crval1_hz=110.201354e9,
        cdelt1_hz=61e3,
        crpix1=1.0,
        restfreq_hz=110.201354e9,
        specsys="LSRK",
        veldef="RADIO",
        refchan=1,
    )
    info = DatasetInfo(telescope="MyDish", project="OTF", object_name="Map", spectral_axis=axis, shared_meta={"TAU_DEF": "zenith tau"})
    w = SDRadioSpectralSDFITSWriter(
        n_chan=1024,
        site=site,
        info=info,
        store_freq_column=True,
        duplicate_data_columns=True,
        chunk_size=5000,
        out_basename="otf_run1",
    )

    mjd0 = 60200.0
    for i in range(12000):
        spec = np.random.normal(0, 1, 1024).astype(np.float32)
        w.add_row(
            time_mjd=mjd0 + i * (0.1 / 86400.0),
            scanid=1,
            subscan=int(i // 1000),
            intgrp=int(i // 200),        # 例：200dumpを1束扱い
            obsmode="ON",
            data=spec,
            exposure_s=0.1,
            ra_deg=83.0 + 0.001*(i % 1000),
            dec_deg=-5.0 + 0.001*(i // 1000),
            srcframe="RADEC",
            src_long_deg=83.0,
            src_lat_deg=-5.0,
            scanframe="RADEC",
            scan_radesys="ICRS",
            scan_equinox=2000.0,
            scan_x_deg=0.001*(i % 1000),
            scan_y_deg=0.001*(i // 1000),
            tsys_k=200.0,
            tau=0.10,
            el_enc_deg=45.0,
            az_enc_deg=120.0,
        )

    w.close()


if __name__ == "__main__":
    # SDRadioSpectralSDFITSWriter.print_schema()
    _example_small()
    # _example_parts()
