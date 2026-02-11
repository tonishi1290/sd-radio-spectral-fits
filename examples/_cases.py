# =============================================================================
# 7) main: 観測モード別の書き出し例
# =============================================================================
import numpy as np
from astropy.time import Time

from sd_radio_spectral_fits import (
    SDRadioSpectralSDFITSWriter,
    Site,
    DatasetInfo,
    SpectralAxisUniform,
    Efficiency,
)


def _mk_site_info_axis_line_115GHz(*, specsys: str, store_freq_column: bool):
    """
    例として 115 GHz 帯（CO(1-0)等）を想定した Site/DatasetInfo/SpectralAxis の雛形を返す。
    specsys: "TOPOCENT" or "LSRK" etc.
    store_freq_column: True の場合、行ごとFREQ列を持てる（OTF/ドップラー補正後などに強い）。
    """
    site = Site(lat_deg=35.9400, lon_deg=138.4720, elev_m=1350.0)

    # 例：1024ch, 61kHz spacing（適宜置換）
    axis = SpectralAxisUniform(
        crval1_hz=115.2712018e9,   # 例: CO(1-0) 115.2712018 GHz 付近（観測中心周波数として仮）
        cdelt1_hz=61e3,
        crpix1=1.0,
        restfreq_hz=115.2712018e9, # 速度変換の基準（必須級）
        specsys=specsys,
        veldef="RADIO",
        refchan=1,
    )

    info = DatasetInfo(
        telescope="MyDish",
        observer="Observer",
        project="DEMO",
        object_name="Target",
        radesys="ICRS",
        equinox=2000.0,

        # SRCFRAME='RADEC' の「SRC_LONG/LAT の基準系」
        src_radesys="ICRS",
        src_equinox=2000.0,

        eff=Efficiency(effstat="UNKNOWN"),
        refr_included_in_corr=True,
        doppler_tracking_applied=False,  # オンライン追尾をしたなら True

        spectral_axis=axis,
        shared_meta={
            "CALMETHOD": "CHOPPER_WHEEL",
            "TAU_DEF": "zenith tau at observing frequency (dimensionless)",
            "BACKEND": "XFFTS",
        },
    )

    # store_freq_column=False なら WCS はヘッダにのみ載る（axisは必須）
    # store_freq_column=True なら FREQ 列（ベクトル）を各行に持てる
    return site, info


def _rand_spec(n_chan: int, seed: int = 0, scale: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (rng.normal(0.0, 1.0, n_chan).astype(np.float32)) * float(scale)


def _integrate_spectra(spec_list: list[np.ndarray]) -> np.ndarray:
    """
    単純平均で積分（例）。実際は重み付き・フラグ・基線等を考慮してください。
    """
    a = np.stack(spec_list, axis=0).astype(np.float64)
    return np.mean(a, axis=0).astype(np.float32)


def _ta_star_from_on_off_hot(on: np.ndarray, off: np.ndarray, hot: np.ndarray,
                            Thot_K: float = 300.0) -> np.ndarray:
    """
    かなり単純化した例：
      Ta* ≈ Thot * (ON - OFF) / (HOT - OFF)
    実運用のあなたの較正式に置換してください（gain・spillover・atm等）。
    """
    den = (hot - off).astype(np.float64)
    den[np.abs(den) < 1e-12] = np.nan
    ta = Thot_K * (on - off).astype(np.float64) / den
    return ta.astype(np.float32)


def _write_case1_raw_on_off_hot():
    """
    1. ON/OFF/HOT の生データをそのまま保存
       - CALSTAT="RAW"
       - OBSMODE を "ON"/"OFF"/"HOT" で区別
       - スペクトル軸は通常 TOPOCENT を推奨（生に近い意味）
    """
    site, info = _mk_site_info_axis_line_115GHz(specsys="TOPOCENT", store_freq_column=False)

    w = SDRadioSpectralSDFITSWriter(
        n_chan=1024,
        site=site,
        info=info,
        store_freq_column=False,
        duplicate_data_columns=True,
    )

    mjd0 = 60200.0
    ra, dec = 83.82208, -5.39111  # 例
    glon, glat = np.nan, np.nan   # RA/DECだけ与えても writer がGLON/GLATを埋める

    # 例：1サイクル分（ON/OFF/HOT）を同一scanidに入れる
    seq = [
        ("OFF", 1.0),
        ("HOT", 30.0),
        ("ON",  3.0),
    ]
    for i, (mode, amp) in enumerate(seq):
        spec = _rand_spec(1024, seed=100+i, scale=amp)

        w.add_row(
            time_mjd=mjd0 + i*(0.1/86400.0),
            scanid=1, subscan=0, intgrp=0,
            obsmode=mode,
            data=spec,
            exposure_s=0.1,

            calstat="RAW",
            tsys_k=np.nan,
            tau=np.nan,

            ra_deg=ra, dec_deg=dec,
            glon_deg=glon, glat_deg=glat,

            # SRC: その行が狙っている中心（たとえば同一ターゲット中心）
            srcframe="RADEC",
            src_long_deg=ra, src_lat_deg=dec,

            # SCAN: 今回は「オフセット無し」扱い
            scanframe="RADEC",
            scan_radesys="ICRS",
            scan_equinox=2000.0,
            scan_x_deg=0.0, scan_y_deg=0.0,

            az_enc_deg=120.0, el_enc_deg=45.0,
        )

    w.write("case1_raw_on_off_hot.fits", overwrite=True)


def _write_case2_tastar_dump_multi_points():
    """
    2. dumpデータから求めた Ta*（積分前）。
       - CALSTAT="TASTAR"
       - OBSMODE は基本 "ON"（較正後スペクトルとして扱う）
       - 複数観測点（点A, 点B）を同一FITSに混在：intgrp で束ねる運用例
    """
    site, info = _mk_site_info_axis_line_115GHz(specsys="TOPOCENT", store_freq_column=False)
    info.object_name = "MultiPoint_TaStar"

    w = SDRadioSpectralSDFITSWriter(
        n_chan=1024, site=site, info=info,
        store_freq_column=False, duplicate_data_columns=True
    )

    mjd0 = 60200.1

    # 例：2点観測（Point A / Point B）、各点で dump が数本ある
    points = [
        ("A", 83.82208, -5.39111),
        ("B", 83.82500, -5.39000),
    ]

    row_i = 0
    for pidx, (pname, ra, dec) in enumerate(points):
        # 例：各点で OFF/HOT を1回ずつ取り、ONを複数dump取ったと仮定
        off = _rand_spec(1024, seed=2000+pidx, scale=1.0)
        hot = _rand_spec(1024, seed=2100+pidx, scale=30.0)

        for d in range(5):  # 5 dumps
            on = _rand_spec(1024, seed=2200+pidx*10+d, scale=3.0)
            ta = _ta_star_from_on_off_hot(on, off, hot, Thot_K=300.0)

            # 観測点を束ねる例：intgrp=点ID、subscan=点ID、scanidは観測ブロック
            w.add_row(
                time_mjd=mjd0 + row_i*(0.1/86400.0),
                scanid=10,
                subscan=pidx,
                intgrp=100+pidx,
                obsmode="ON",
                data=ta,
                exposure_s=0.1,

                object_name=f"Point{pname}",
                calstat="TASTAR",
                tsys_k=200.0,
                tau=0.10,

                ra_deg=ra, dec_deg=dec,
                srcframe="RADEC", src_long_deg=ra, src_lat_deg=dec,

                # 点観測なので scanframe=RADEC、オフセット0
                scanframe="RADEC",
                scan_radesys="ICRS",
                scan_equinox=2000.0,
                scan_x_deg=0.0, scan_y_deg=0.0,

                az_enc_deg=120.0, el_enc_deg=45.0,
            )
            row_i += 1

    w.write("case2_tastar_dump_multi_points.fits", overwrite=True)


def _write_case3_integrated_tastar_multi_points():
    """
    3. 積分した Ta*（複数観測点）。
       - CALSTAT="TASTAR"（もしくは "TA" にしても良いが、ここではTa*として統一）
       - exposure_s は「積分後の有効積分時間」を入れる
    """
    site, info = _mk_site_info_axis_line_115GHz(specsys="TOPOCENT", store_freq_column=False)
    info.object_name = "MultiPoint_TaStar_INT"

    w = SDRadioSpectralSDFITSWriter(
        n_chan=1024, site=site, info=info,
        store_freq_column=False, duplicate_data_columns=True
    )

    mjd0 = 60200.2

    points = [
        ("A", 83.82208, -5.39111),
        ("B", 83.82500, -5.39000),
        ("C", 83.83000, -5.39500),
    ]

    row_i = 0
    for pidx, (pname, ra, dec) in enumerate(points):
        off = _rand_spec(1024, seed=3000+pidx, scale=1.0)
        hot = _rand_spec(1024, seed=3100+pidx, scale=30.0)

        dumps = []
        for d in range(8):
            on = _rand_spec(1024, seed=3200+pidx*100+d, scale=3.0)
            dumps.append(_ta_star_from_on_off_hot(on, off, hot, Thot_K=300.0))

        ta_int = _integrate_spectra(dumps)

        w.add_row(
            time_mjd=mjd0 + row_i*(1.0/86400.0),  # 例：点ごとに時刻が進む
            scanid=20,
            subscan=pidx,
            intgrp=200+pidx,
            obsmode="ON",
            data=ta_int,
            exposure_s=0.1 * len(dumps),  # dump積分時間×本数

            object_name=f"Point{pname}",
            calstat="TASTAR",
            tsys_k=200.0,
            tau=0.10,

            ra_deg=ra, dec_deg=dec,
            srcframe="RADEC", src_long_deg=ra, src_lat_deg=dec,

            scanframe="RADEC",
            scan_radesys="ICRS",
            scan_equinox=2000.0,
            scan_x_deg=0.0, scan_y_deg=0.0,

            az_enc_deg=120.0, el_enc_deg=45.0,
        )
        row_i += 1

    w.write("case3_tastar_integrated_multi_points.fits", overwrite=True)


def _write_case4_integrated_doppler_corrected():
    """
    4. 解析してドップラー補正した積分データ（複数観測点）。
       - 典型運用：
         (A) 生の周波数は TOPOCENT のまま保持し、VFRAME 等に補正量を持つ
         (B) 解析後のスペクトル軸（FREQ列）を LSRK 等へ変換して持つ
       ここでは (B) を例示：
         - store_freq_column=True にして各行FREQを持つ
         - info.spectral_axis.specsys="LSRK" を宣言（解析後の基準系）
    """
    site, info = _mk_site_info_axis_line_115GHz(specsys="LSRK", store_freq_column=True)
    info.object_name = "MultiPoint_TaStar_INT_DOPCOR"
    # オフライン補正であることを明記（自由に鍵追加可能）
    info.shared_meta["DOPPLER_CORR"] = "applied offline (analysis stage)"

    w = SDRadioSpectralSDFITSWriter(
        n_chan=1024, site=site, info=info,
        store_freq_column=True, duplicate_data_columns=True
    )

    mjd0 = 60200.3
    points = [
        ("A", 83.82208, -5.39111),
        ("B", 83.82500, -5.39000),
    ]

    for pidx, (pname, ra, dec) in enumerate(points):
        # 解析済みの積分スペクトル（例：既にTa*で積分済み）
        spec = _rand_spec(1024, seed=4000+pidx, scale=2.0)

        # 周波数ベクトル：ここでは「uniformから生成したTOPO」を渡し、
        # writer側で（現状実装の）LSRK-like変換を適用してFREQ列に入れる想定。
        # ※実運用では、あなたの確定したVLSRK変換で作ったFREQをそのまま渡してOKです。
        freq_topo = w._freq_from_uniform_axis()

        w.add_row(
            time_mjd=mjd0 + pidx*(2.0/86400.0),
            scanid=30, subscan=pidx, intgrp=300+pidx,
            obsmode="ON",
            data=spec,
            exposure_s=10.0,

            object_name=f"Point{pname}",
            calstat="TASTAR",
            tsys_k=200.0,
            tau=0.10,

            ra_deg=ra, dec_deg=dec,
            srcframe="RADEC", src_long_deg=ra, src_lat_deg=dec,

            scanframe="RADEC",
            scan_radesys="ICRS", scan_equinox=2000.0,
            scan_x_deg=0.0, scan_y_deg=0.0,

            # 解析で求めた補正速度を別途格納したい場合（任意）
            v_frame_mps=np.nan,

            # 解析後FREQを入れる（store_freq_column=True）
            freq_hz=freq_topo,
        )

    w.write("case4_integrated_doppler_corrected.fits", overwrite=True)


def _write_case5_otf_raw_hot_off_continuous_on():
    """
    5. OTFの生データ（HOT, OFF, 連続するON）
       - 巨大化するので chunk_size parts を推奨
       - scanid: 観測ブロック
       - subscan: 走査ライン番号
       - intgrp: 較正参照グループ（OFF更新ごと等）
       - scanframe: 走査オフセットの定義（例：RADECオフセット or AZEL_GEOオフセット）
    """
    site, info = _mk_site_info_axis_line_115GHz(specsys="TOPOCENT", store_freq_column=False)
    info.object_name = "OTF_RAW"

    w = SDRadioSpectralSDFITSWriter(
        n_chan=1024, site=site, info=info,
        store_freq_column=False, duplicate_data_columns=True,
        chunk_size=5000, out_basename="case5_otf_raw"
    )

    mjd0 = 60200.4
    ra0, dec0 = 83.8, -5.39  # マップ中心（例）

    row_i = 0
    n_lines = 3
    dumps_per_line = 2000  # 例（実際はもっと多い）

    # 走査前に HOT/OFF を入れる例
    hot = _rand_spec(1024, seed=5000, scale=30.0)
    off = _rand_spec(1024, seed=5001, scale=1.0)

    # HOT
    w.add_row(
        time_mjd=mjd0 + row_i*(0.1/86400.0),
        scanid=100, subscan=0, intgrp=0,
        obsmode="HOT", data=hot, exposure_s=0.1,
        calstat="RAW",
        ra_deg=ra0, dec_deg=dec0,
        srcframe="RADEC", src_long_deg=ra0, src_lat_deg=dec0,
        scanframe="RADEC", scan_radesys="ICRS", scan_equinox=2000.0,
        scan_x_deg=0.0, scan_y_deg=0.0,
        az_enc_deg=120.0, el_enc_deg=45.0,
    )
    row_i += 1

    # OFF
    w.add_row(
        time_mjd=mjd0 + row_i*(0.1/86400.0),
        scanid=100, subscan=0, intgrp=0,
        obsmode="OFF", data=off, exposure_s=0.1,
        calstat="RAW",
        ra_deg=ra0, dec_deg=dec0,
        srcframe="RADEC", src_long_deg=ra0, src_lat_deg=dec0,
        scanframe="RADEC", scan_radesys="ICRS", scan_equinox=2000.0,
        scan_x_deg=0.0, scan_y_deg=0.0,
        az_enc_deg=120.0, el_enc_deg=45.0,
    )
    row_i += 1

    # 連続ON（OTF）
    for line in range(n_lines):
        for j in range(dumps_per_line):
            # 例：RA方向にスキャン、DECはラインごとに変える
            dra = (j - dumps_per_line/2) * (1.0/3600.0)  # 1 arcsec step相当（例）
            ddec = (line - (n_lines-1)/2) * (30.0/3600.0) # ライン間隔30 arcsec（例）
            ra = ra0 + dra
            dec = dec0 + ddec

            spec = _rand_spec(1024, seed=5100 + line*100000 + j, scale=3.0)

            w.add_row(
                time_mjd=mjd0 + row_i*(0.1/86400.0),
                scanid=100,
                subscan=line,          # ライン番号
                intgrp=line,           # OFF更新単位などに合わせて設計
                obsmode="ON",
                data=spec,
                exposure_s=0.1,
                calstat="RAW",

                ra_deg=ra, dec_deg=dec,
                srcframe="RADEC", src_long_deg=ra0, src_lat_deg=dec0,  # マップ中心を意図座標に入れる例

                # SCAN: オフセットはRADECで表現（中心からの差）
                scanframe="RADEC",
                scan_radesys="ICRS",
                scan_equinox=2000.0,
                scan_x_deg=(ra - ra0),   # 単純差分（厳密にはcos(dec)等を考慮する設計も可）
                scan_y_deg=(dec - dec0),

                az_enc_deg=120.0, el_enc_deg=45.0,
            )
            row_i += 1

    w.close()  # parts flush + manifest


def _write_case6_sun_scan_raw_radec_target_azel_offset():
    """
    6. 太陽スキャン生データ
       - 太陽の中心（狙い）は RA/DEC (J2000) で与える
       - 走査オフセットは Az/El（= AZEL_GEO）で与える
       - したがって：
         srcframe="RADEC", src_long/lat = (Sun RA, Sun Dec)
         scanframe="AZEL_GEO", scan_x/y = (dAz, dEl)  [deg]
       - 行の絶対座標（RA/DEC, GLON/GLAT）は、各dumpの実際の指向（推定）を入れるのが理想。
         ここでは例として、中心RA/DECに近い値を入れる（オフセットはSCAN側に保持）。
    """
    site, info = _mk_site_info_axis_line_115GHz(specsys="TOPOCENT", store_freq_column=False)
    info.object_name = "SUN_SCAN_RAW"

    w = SDRadioSpectralSDFITSWriter(
        n_chan=1024, site=site, info=info,
        store_freq_column=False, duplicate_data_columns=True
    )

    mjd0 = 60200.5

    # 例：太陽中心（観測計算で得たJ2000中心；ここではダミー）
    sun_ra_j2000 = 300.0
    sun_dec_j2000 = -20.0

    # 例：Az方向に往復スキャン（±1 deg）、Elは0固定
    n = 200
    az_span_deg = 2.0
    for i in range(n):
        # -1..+1 deg
        daz = (i/(n-1) - 0.5) * az_span_deg
        delv = 0.0

        spec = _rand_spec(1024, seed=6000+i, scale=50.0)  # 太陽なので強い例

        w.add_row(
            time_mjd=mjd0 + i*(0.1/86400.0),
            scanid=500, subscan=0, intgrp=0,
            obsmode="ON",
            data=spec,
            exposure_s=0.1,

            object_name="SUN",
            calstat="RAW",
            tsys_k=np.nan,
            tau=np.nan,

            # 絶対座標（例：中心付近を入れる）
            ra_deg=sun_ra_j2000,
            dec_deg=sun_dec_j2000,

            # SRC: 太陽中心（RADEC/J2000）
            srcframe="RADEC",
            src_long_deg=sun_ra_j2000,
            src_lat_deg=sun_dec_j2000,

            # SCAN: オフセットは AZEL_GEO
            scanframe="AZEL_GEO",
            scan_x_deg=daz,   # dAz
            scan_y_deg=delv,  # dEl

            # エンコーダ（例）
            az_enc_deg=180.0 + daz,
            el_enc_deg=45.0 + delv,
        )

    w.write("case6_sun_scan_raw.fits", overwrite=True)


def _write_case7_sun_scan_total_power_only():
    """
    7. 太陽スキャンの TP（total power）データのみ
       - n_chan=1 として 1点値を SPECTRUM[0] に入れる
       - 周波数軸が意味を持たない（あるいは固定）場合でも、
         writer仕様上 store_freq_column=False なら spectral_axis は必要なので、
         info.spectral_axis はダミーでも良い（restfreq=0でも可）。
       - 走査の座標設計は case6 と同様：SRCはRADEC、SCANはAZEL_GEO等
    """
    # ここだけ continuum 扱いのため、axisは「形式上」持たせる（値は用途に合わせて）
    site = Site(lat_deg=35.9400, lon_deg=138.4720, elev_m=1350.0)
    axis = SpectralAxisUniform(
        crval1_hz=0.0, cdelt1_hz=1.0, crpix1=1.0,
        restfreq_hz=0.0,
        specsys="TOPOCENT",
        veldef="RADIO",
        refchan=1,
    )
    info = DatasetInfo(
        telescope="MyDish",
        observer="Observer",
        project="DEMO",
        object_name="SUN_SCAN_TP",
        radesys="ICRS",
        equinox=2000.0,
        src_radesys="ICRS",
        src_equinox=2000.0,
        spectral_axis=axis,
        shared_meta={"DATATYPE": "TOTAL_POWER"},
    )

    w = SDRadioSpectralSDFITSWriter(
        n_chan=1,
        site=site,
        info=info,
        store_freq_column=False,
        duplicate_data_columns=True,
    )

    mjd0 = 60200.6
    sun_ra_j2000 = 300.0
    sun_dec_j2000 = -20.0

    n = 200
    az_span_deg = 2.0
    for i in range(n):
        daz = (i/(n-1) - 0.5) * az_span_deg
        delv = 0.0

        # TP値（例）
        tp = float(1e3 + 200*np.cos(2*np.pi*i/(n-1)))  # 適当な変動例

        w.add_row(
            time_mjd=mjd0 + i*(0.1/86400.0),
            scanid=700, subscan=0, intgrp=0,
            obsmode="ON",
            data=tp,                # n_chan=1 なら float OK
            exposure_s=0.1,

            object_name="SUN",
            calstat="RAW",

            ra_deg=sun_ra_j2000,
            dec_deg=sun_dec_j2000,

            srcframe="RADEC",
            src_long_deg=sun_ra_j2000,
            src_lat_deg=sun_dec_j2000,

            scanframe="AZEL_GEO",
            scan_x_deg=daz,
            scan_y_deg=delv,

            az_enc_deg=180.0 + daz,
            el_enc_deg=45.0 + delv,
        )

    w.write("case7_sun_scan_tp_only.fits", overwrite=True)


def main():
    """
    観測モード別：FITS書き出しの具体例
    """
    _write_case1_raw_on_off_hot()
    _write_case2_tastar_dump_multi_points()
    _write_case3_integrated_tastar_multi_points()
    _write_case4_integrated_doppler_corrected()
    _write_case5_otf_raw_hot_off_continuous_on()
    _write_case6_sun_scan_raw_radec_target_azel_offset()
    _write_case7_sun_scan_total_power_only()
