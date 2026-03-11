#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plotting.py の cookbook 的な使用例集

このファイルでは、以下をまとめて示します。
1. 2D FITS をそのまま描く
2. 2D extension を指定して描く
3. 2D extension を読む場合と、3D cube から再計算する場合の違い
4. (header, data) から描く
5. 3D cube から moment0 / channel_sum / channel_slice を作る
6. provisional moment と final moment を作る
7. 複数 contour を重ねる
8. 銀河座標 (GLON/GLAT) の map を描く
9. RGB 合成を 2D / 3D から作る
10. normalize 手法を切り替える
11. color scale の下限・上限 (cmin/cmax) を指定する
12. HPBW が header に無い場合に orig_hpbw_arcsec を明示する

注意:
- ここではファイル名は例です。手元の実データに合わせて書き換えてください。
- すべての例を一気に実行する必要はありません。
- plotting_checked_v7.py が同じディレクトリにある想定ですが、
  パッケージに組み込んだ後は import 部分を変えればそのまま使えます。
"""

from pathlib import Path

import numpy as np
from astropy.io import fits

try:
    # ローカル単独ファイルとして使う場合
    from plotting import (
        build_normalize,
        make_2d_map,
        make_final_moment,
        make_provisional_moment,
        make_rgb_map,
        plot_map,
        plot_rgb,
        resolve_map_input,
    )
except ImportError:
    # パッケージへ入れた後の想定
    from sd_radio_spectral_fits.map_3d.plotting import (
        build_normalize,
        make_2d_map,
        make_final_moment,
        make_provisional_moment,
        make_rgb_map,
        plot_map,
        plot_rgb,
        resolve_map_input,
    )


# ---------------------------------------------------------------------
# 手元のデータに合わせて編集する場所
# ---------------------------------------------------------------------

DATA = {
    # 2D FITS の例
    "map2d_eq": "moment0_equatorial.fits",
    "map2d_gal": "moment0_galactic.fits",

    # 3D cube の例
    "cube_eq": "ori-kl-12co_bsub.fits",
    "cube_gal": "galactic_cube.fits",

    # mask extension を含む cube の例
    "cube_with_masks": "cube_with_masks.fits",

    # 別データの contour 用
    "contour_cube": "ori-kl-13co_bsub.fits",
}

OUTDIR = Path("plotting_examples")
OUTDIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------
# まず重要な考え方
# ---------------------------------------------------------------------

EXPLANATION_AUTO_PROVISIONAL = r"""
auto で provisional moment を作るときの意味
-------------------------------------------

provisional moment は、「本番の MASK3D が無い、または作る前の段階で、
まず信号らしい場所だけを使って仮の moment map を作る」ためのものです。

auto の探索順は次の通りです。

1. LINECAND3D があればそれを使う
   - True = この voxel は線信号候補
   - つまり、True の場所をそのまま積分します。
   - これが最も直接的に「ここが信号」と書いてある情報です。

2. LINECAND3D が無く、BASESUP3D があればその補集合を使う
   - BASESUP3D の True = baseline を決めるのに使ってよい領域
   - つまり、baseline 側の領域です。
   - provisional moment では、その反対側（補集合）を「信号候補」とみなして積分します。

3. さらに BASESUP3D も無く、LINEFREE があればその補集合を使う
   - LINEFREE の True = line-free channel（線が無いチャネル）
   - つまり、baseline 用チャネルです。
   - provisional moment では、その反対側（line-free ではないチャネル）を積分します。

ひとことで言うと、
「信号がある場所を直接示すマスクがあればそれを使う。
それが無ければ、『ここは baseline 側です』『ここは line-free です』
と書かれた情報の反対側を、仮の signal 領域として使う」
という意味です。

実務的には次の理解で十分です。
- LINECAND3D がある: 「ここを積分すればよい」
- BASESUP3D しかない: 「baseline に使った場所以外を積分する」
- LINEFREE しかない: 「line-free ではないチャネルを積分する」
"""


def print_auto_provisional_explanation():
    print(EXPLANATION_AUTO_PROVISIONAL)


def save_result(result, filename: str, dpi: int = 250):
    out = OUTDIR / filename
    result["fig"].savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"[saved] {out}")


# ---------------------------------------------------------------------
# 1) 2D FITS をそのまま描く
# ---------------------------------------------------------------------

def example_2d_identity():
    """
    最も単純な 2D 表示。
    2D FITS を resolve -> plot する例。
    """
    map2d = make_2d_map(
        DATA["map2d_eq"],
        mode="identity",
    )

    result = plot_map(
        map2d,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        stretch_a=0.08,
        beam="header",
        title="2D FITS identity example",
        show=False,
    )
    save_result(result, "ex01_identity_2d.png")


# ---------------------------------------------------------------------
# 2) 2D extension を指定して描く
# ---------------------------------------------------------------------

def example_2d_extension():
    """
    1つの FITS に複数 extension がある場合。
    例: PRIMARY は空で、MOMENT0 や RMS が extension に入っているケース。
    """
    result = plot_map(
        source=DATA["cube_with_masks"],
        ext="MOMENT0",
        cmap="viridis",
        norm_mode="linear",
        norm_percentile=(1.0, 99.0),
        title="2D extension example: MOMENT0",
        show=False,
    )
    save_result(result, "ex02_extension_moment0.png")

    # RMS は linear/log を切り替えたいことが多い
    result = plot_map(
        source=DATA["cube_with_masks"],
        ext="BASE_RMS",
        cmap="magma",
        norm_mode="log",
        norm_percentile=(1.0, 99.5),
        title="2D extension example: BASE_RMS (log)",
        show=False,
    )
    save_result(result, "ex02_extension_base_rms_log.png")


# ---------------------------------------------------------------------
# 2.5) 2D extension を読む場合と、3D cube から再計算する場合の違い
# ---------------------------------------------------------------------

def example_precomputed_moment0_vs_recomputed_moment0():
    """
    ext="MOMENT0" は既に 2D です。
    したがって vel_range を与えても、その MOMENT0 を再計算するわけではありません。

    - ext="MOMENT0" を読む: 既存の 2D map をそのまま読む
    - 3D cube を指定して mode="moment0": 指定した速度範囲で新しく積分する
    """
    # A) 既に用意されている MOMENT0 extension をそのまま読む
    precomputed = make_2d_map(
        DATA["cube_with_masks"],
        ext="MOMENT0",
        mode="map",
        target_hpbw_arcsec=400.0,
        orig_hpbw_arcsec=350.0,
    )
    result = plot_map(
        precomputed,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        stretch_a=0.10,
        beam="header",
        title="Precomputed MOMENT0 extension (2D passthrough)",
        show=False,
    )
    save_result(result, "ex02b_precomputed_moment0_passthrough.png")

    # B) 元の 3D cube から速度範囲を指定して新しく moment0 を作る
    remade = make_2d_map(
        DATA["cube_eq"],
        ext=0,
        mode="moment0",
        vel_range=(0.0, 20.0),
        spectral_unit="km/s",
        target_hpbw_arcsec=400.0,
        orig_hpbw_arcsec=350.0,
    )
    result = plot_map(
        remade,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        stretch_a=0.10,
        beam="header",
        title="Moment0 remade from 3D cube over 0-20 km/s",
        show=False,
    )
    save_result(result, "ex02b_remade_moment0_from_3dcube.png")


# ---------------------------------------------------------------------
# 3) (header, data) から描く
# ---------------------------------------------------------------------

def example_header_data_input():
    """
    既に Python 上で header, data を持っている場合。
    例えば別処理で moment0 を作ってから描くときの形。
    """
    data = fits.getdata(DATA["map2d_eq"])
    header = fits.getheader(DATA["map2d_eq"])

    result = plot_map(
        data=data,
        header=header,
        cmap="turbo",
        norm_mode="sqrt",
        norm_percentile=(1.0, 99.5),
        title="(header, data) input example",
        show=False,
    )
    save_result(result, "ex03_header_data.png")


# ---------------------------------------------------------------------
# 4) 3D cube -> moment0
# ---------------------------------------------------------------------

def example_moment0_with_target_hpbw():
    """
    ユーザーが最初に示した典型ユースケースに近い例。
    target_hpbw_arcsec を使う場合、元 HPBW が header に無いなら
    orig_hpbw_arcsec を明示する。
    """
    map2d = make_2d_map(
        DATA["cube_eq"],
        mode="moment0",
        vel_range=(0.0, 20.0),
        spectral_unit="km/s",
        target_hpbw_arcsec=200.0,
        orig_hpbw_arcsec=100.0,   # header に BMAJ/BMIN が無い場合は明示
    )

    result = plot_map(
        map2d,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        stretch_a=0.10,
        contours=[
            {
                "levels": "auto",
                "colors": "white",
                "linewidths": 0.8,
                "alpha": 0.6,
            }
        ],
        beam="header",
        title="Moment 0 with target HPBW",
        show=False,
    )
    save_result(result, "ex04_moment0_target_hpbw.png")


# ---------------------------------------------------------------------
# 5) channel_sum / channel_mean / channel_slice
# ---------------------------------------------------------------------

def example_channel_based_maps():
    """
    velocity ではなく channel 指定で 2D map を作る例。
    """
    channel_sum = make_2d_map(
        DATA["cube_eq"],
        mode="channel_sum",
        chan_range=(20, 50),
        smooth_fwhm_arcsec=60.0,
    )
    result = plot_map(
        channel_sum,
        cmap="viridis",
        norm_mode="linear",
        norm_percentile=(1.0, 99.5),
        title="Channel sum (20:50)",
        show=False,
    )
    save_result(result, "ex05_channel_sum.png")

    channel_mean = make_2d_map(
        DATA["cube_eq"],
        mode="channel_mean",
        chan_range=(20, 50),
    )
    result = plot_map(
        channel_mean,
        cmap="viridis",
        norm_mode="sqrt",
        norm_percentile=(1.0, 99.5),
        title="Channel mean (20:50)",
        show=False,
    )
    save_result(result, "ex05_channel_mean.png")

    channel_slice = make_2d_map(
        DATA["cube_eq"],
        mode="channel_slice",
        chan_range=(35, 35),
    )
    result = plot_map(
        channel_slice,
        cmap="plasma",
        norm_mode="linear",
        norm_percentile=(1.0, 99.5),
        title="Single channel slice (35)",
        show=False,
    )
    save_result(result, "ex05_channel_slice.png")


# ---------------------------------------------------------------------
# 6) provisional moment
# ---------------------------------------------------------------------

def example_provisional_auto_and_manual():
    """
    provisional moment の自動選択と手動選択。
    自動では LINECAND3D -> BASESUP3D -> LINEFREE補集合 の順に探索する。
    """
    pmap_auto = make_provisional_moment(
        DATA["cube_with_masks"],
        prefer="auto",
        vel_range=(0.0, 20.0),
        spectral_unit="km/s",
        nan_fill=True,
    )
    result = plot_map(
        pmap_auto,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        title="Provisional moment (auto)",
        show=False,
    )
    save_result(result, "ex06_provisional_auto.png")

    pmap_linecand = make_provisional_moment(
        DATA["cube_with_masks"],
        prefer="linecand3d",
        vel_range=(0.0, 20.0),
        spectral_unit="km/s",
        zero_fill=True,
        nan_fill=False,
    )
    result = plot_map(
        pmap_linecand,
        cmap="viridis",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        title="Provisional moment (LINECAND3D, zero fill)",
        show=False,
    )
    save_result(result, "ex06_provisional_linecand_zero_fill.png")


# ---------------------------------------------------------------------
# 7) final moment
# ---------------------------------------------------------------------

def example_final_moment():
    fmap = make_final_moment(
        DATA["cube_with_masks"],
        vel_range=(0.0, 20.0),
        spectral_unit="km/s",
        smooth_fwhm_arcsec=45.0,
    )
    result = plot_map(
        fmap,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        beam="header",
        title="Final moment from MASK3D",
        show=False,
    )
    save_result(result, "ex07_final_moment.png")


# ---------------------------------------------------------------------
# 8) 複数 contour を重ねる
# ---------------------------------------------------------------------

def example_multiple_contours():
    """
    pseudo color と contour を別データで重ねる例。
    同じ cube から別速度範囲を contour にしてもよいし、
    別の FITS を contour source にしてもよい。
    """
    base = make_2d_map(
        DATA["cube_eq"],
        mode="moment0",
        vel_range=(0.0, 20.0),
        spectral_unit="km/s",
    )

    contours = [
        {
            # 同じ cube の別速度範囲
            "source": DATA["cube_eq"],
            "mode": "moment0",
            "vel_range": (5.0, 10.0),
            "spectral_unit": "km/s",
            "levels": "auto",
            "colors": "white",
            "linewidths": 0.8,
            "alpha": 0.8,
            "label": "12CO 5-10 km/s",
        },
        {
            # 別データセット
            "source": DATA["contour_cube"],
            "mode": "moment0",
            "vel_range": (0.0, 20.0),
            "spectral_unit": "km/s",
            "levels": "auto",
            "colors": "cyan",
            "linewidths": 0.8,
            "linestyles": "dashed",
            "alpha": 0.8,
            "label": "13CO 0-20 km/s",
        },
    ]

    result = plot_map(
        base,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        contours=contours,
        title="Multiple contour overlays",
        show=False,
    )
    save_result(result, "ex08_multiple_contours.png")


# ---------------------------------------------------------------------
# 9) 銀河座標の map
# ---------------------------------------------------------------------

def example_galactic_coordinates():
    """
    header の CTYPE が GLON/GLAT であれば、軸ラベルは自動で
    Galactic Longitude / Galactic Latitude になる。
    """
    gmap = make_2d_map(
        DATA["cube_gal"],
        mode="moment0",
        vel_range=(-20.0, 20.0),
        spectral_unit="km/s",
    )
    result = plot_map(
        gmap,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        title="Galactic coordinate moment map",
        show=False,
    )
    save_result(result, "ex09_galactic_coords.png")


# ---------------------------------------------------------------------
# 10) normalize を切り替える
# ---------------------------------------------------------------------

def example_normalize_variants():
    """
    同じ map を linear / sqrt / log / asinh / power で描き比べる例。
    """
    map2d = make_2d_map(
        DATA["cube_eq"],
        mode="moment0",
        vel_range=(0.0, 20.0),
        spectral_unit="km/s",
    )

    # external norm を build_normalize で作って渡す例
    custom_norm = build_normalize(
        map2d.data,
        mode="power",
        percentile=(1.0, 99.5),
        power_gamma=0.5,
    )
    result = plot_map(
        map2d,
        cmap="turbo",
        norm=custom_norm,
        title="Power-law normalize (external norm)",
        show=False,
    )
    save_result(result, "ex10_power_norm_external.png")

    for mode in ["linear", "sqrt", "log", "asinh", "power"]:
        kwargs = dict(
            cmap="turbo",
            norm_mode=mode,
            norm_percentile=(1.0, 99.5),
            title=f"Normalize mode: {mode}",
            show=False,
        )
        if mode == "power":
            kwargs["power_gamma"] = 0.6
        if mode == "asinh":
            kwargs["stretch_a"] = 0.08
        try:
            result = plot_map(map2d, **kwargs)
            save_result(result, f"ex10_norm_{mode}.png")
        except ValueError as exc:
            print(f"[skip] normalize={mode}: {exc}")


# ---------------------------------------------------------------------
# 11) RGB 合成: 2D map 3枚から
# ---------------------------------------------------------------------

def example_rgb_from_2d():
    """
    2D map を 3枚用意して RGB 合成する例。
    """
    red = resolve_map_input(DATA["map2d_eq"])
    green = resolve_map_input(DATA["map2d_eq"])
    blue = resolve_map_input(DATA["map2d_eq"])

    rgb = make_rgb_map(
        red,
        green,
        blue,
        red_norm={"norm_mode": "asinh", "percentile": (1.0, 99.5), "stretch_a": 0.08},
        green_norm={"norm_mode": "sqrt", "percentile": (1.0, 99.5)},
        blue_norm={"norm_mode": "linear", "percentile": (1.0, 99.5)},
    )
    result = plot_rgb(
        rgb_source=rgb,
        title="RGB from 2D maps",
        show=False,
    )
    save_result(result, "ex11_rgb_from_2d.png")


# ---------------------------------------------------------------------
# 12) RGB 合成: 3D cube から速度範囲を色分け
# ---------------------------------------------------------------------

def example_rgb_from_3d():
    """
    1つの cube から 3つの速度範囲を作って RGB に割り当てる例。
    """
    rgb = make_rgb_map(
        red={
            "source": DATA["cube_eq"],
            "mode": "moment0",
            "vel_range": (-5.0, 2.0),
            "spectral_unit": "km/s",
        },
        green={
            "source": DATA["cube_eq"],
            "mode": "moment0",
            "vel_range": (2.0, 8.0),
            "spectral_unit": "km/s",
        },
        blue={
            "source": DATA["cube_eq"],
            "mode": "moment0",
            "vel_range": (8.0, 15.0),
            "spectral_unit": "km/s",
        },
        red_norm={"norm_mode": "asinh", "percentile": (1.0, 99.5), "stretch_a": 0.08},
        green_norm={"norm_mode": "asinh", "percentile": (1.0, 99.5), "stretch_a": 0.08},
        blue_norm={"norm_mode": "asinh", "percentile": (1.0, 99.5), "stretch_a": 0.08},
    )
    result = plot_rgb(
        rgb_source=rgb,
        title="RGB from 3D cube velocity ranges",
        show=False,
    )
    save_result(result, "ex12_rgb_from_3d.png")


# ---------------------------------------------------------------------
# 13) 2D extension と contour の混在
# ---------------------------------------------------------------------

def example_2d_base_with_3d_contours():
    """
    pseudo color は 2D extension、contour は 3D cube から moment0 を作る例。
    """
    base = resolve_map_input(
        source=DATA["cube_with_masks"],
        ext="MOMENT0",
    )
    result = plot_map(
        base,
        cmap="turbo",
        norm_mode="asinh",
        norm_percentile=(1.0, 99.5),
        contours=[
            {
                "source": DATA["cube_eq"],
                "mode": "moment0",
                "vel_range": (0.0, 20.0),
                "spectral_unit": "km/s",
                "levels": "auto",
                "colors": "white",
                "linewidths": 0.8,
            }
        ],
        title="2D base + 3D contour overlay",
        show=False,
    )
    save_result(result, "ex13_2d_base_with_3d_contour.png")


# ---------------------------------------------------------------------
# 実行例
# ---------------------------------------------------------------------

def main():
    """
    必要なものだけ呼んでください。
    最初は ex02b, ex04, ex06, ex08, ex12 あたりが分かりやすいです。
    """
    print_auto_provisional_explanation()

    # まずは基本
    # example_2d_identity()
    # example_2d_extension()
    # example_precomputed_moment0_vs_recomputed_moment0()
    # example_header_data_input()
    # example_moment0_with_target_hpbw()

    # channel / mask 系
    # example_channel_based_maps()
    # example_provisional_auto_and_manual()
    # example_final_moment()

    # overlay / galactic / normalize / RGB
    # example_multiple_contours()
    # example_galactic_coordinates()
    # example_normalize_variants()
    # example_rgb_from_2d()
    # example_rgb_from_3d()
    # example_2d_base_with_3d_contours()

    # デモとして 1 つだけ有効にしておく
    example_precomputed_moment0_vs_recomputed_moment0()


if __name__ == "__main__":
    main()
