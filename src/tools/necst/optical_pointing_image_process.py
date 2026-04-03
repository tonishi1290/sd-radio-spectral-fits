#!/usr/bin/env python3
"""CLI tool for optical pointing CSV generation and image-based offset analysis.

This script refactors the workflow in `optical_pointing_make_csv.ipynb` into a
reproducible command-line tool.

Main workflow
-------------
1. Read capture metadata either from a NECST DB table or from an existing CSV.
2. For each image, threshold the grayscale image and choose the largest contour.
3. Convert centroid pixel offsets to arcsec using a fixed pixel scale.
4. Save per-frame results, legacy-compatible CSV/DAT outputs, annotated images,
   and summary plots.

Notes
-----
- The default behavior is intentionally close to the notebook:
  fixed reference center = (4800/2, 3200/2), threshold=50, pixel scale=0.8808.
- If image sizes differ from the expected size, use `--center-mode image` to
  measure offsets with respect to each image center, or `--strict-image-shape`
  to stop immediately.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import cv2
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_DB_TABLE = "necst-NANTEN2-rx-spectra_meta"
DEFAULT_IMAGE_EXT = ".JPG"
DEFAULT_NPIX_X = 4800
DEFAULT_NPIX_Y = 3200
DEFAULT_PIXEL_SCALE_ARCSEC = 0.8808
DEFAULT_THRESHOLD = 50


@dataclass
class Config:
    datetime_tag: str
    datadir: Path
    datadate: str | None
    dbname: str | None
    meta_csv: Path | None
    pic_dir: Path
    output_dir: Path
    figs_dir: Path
    marks_dir: Path
    npix_x: int = DEFAULT_NPIX_X
    npix_y: int = DEFAULT_NPIX_Y
    pixel_scale_arcsec: float = DEFAULT_PIXEL_SCALE_ARCSEC
    threshold: int = DEFAULT_THRESHOLD
    image_ext: str = DEFAULT_IMAGE_EXT
    db_table: str = DEFAULT_DB_TABLE
    save_meta_csv: bool = True
    save_mark_images: bool = True
    dpi: int = 200
    verbose: bool = False
    strict_image_shape: bool = False
    center_mode: str = "fixed"


@dataclass
class DetectionResult:
    file: str
    image_path: str
    cap_az: float
    cap_el: float
    ra: float | None
    dec: float | None
    cap_time: Any
    detected: bool
    reason: str
    image_width: float
    image_height: float
    x_center_ref: float
    y_center_ref: float
    pix_x: float
    pix_y: float
    contour_area: float
    dpix_x: float
    dpix_y: float
    dx_arcsec: float
    dy_arcsec: float


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Optical pointing analysis CLI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--datetime", required=True, help="Image directory timestamp tag, e.g. 20260318_083530")
    parser.add_argument("--datadate", default=None, help="DB timestamp tag, e.g. 20260318_083519")
    parser.add_argument("--datadir", default="./data", help="Base data directory")
    parser.add_argument("--dbname", default=None, help="DB name. Overrides auto-generated necst_opticalpointing_<datadate>")
    parser.add_argument("--meta-csv", default=None, help="Use an existing metadata CSV instead of reading necstdb")
    parser.add_argument("--pic-dir", default=None, help="Directory containing image files. Default: <datadir>/<datetime>")
    parser.add_argument("--output-dir", default=None, help="Output directory. Default: <datadir>/<datetime>/analysis_cli")
    parser.add_argument("--figs-dir", default=None, help="Directory for plots. Default: <output-dir>/figs")
    parser.add_argument("--marks-dir", default=None, help="Directory for annotated images. Default: <output-dir>/marks")
    parser.add_argument("--db-table", default=DEFAULT_DB_TABLE, help="NECST DB table name")
    parser.add_argument("--npix-x", type=int, default=DEFAULT_NPIX_X, help="Expected image width in pixels")
    parser.add_argument("--npix-y", type=int, default=DEFAULT_NPIX_Y, help="Expected image height in pixels")
    parser.add_argument("--pixel-scale-arcsec", type=float, default=DEFAULT_PIXEL_SCALE_ARCSEC, help="Pixel scale in arcsec/pixel")
    parser.add_argument("--threshold", type=int, default=DEFAULT_THRESHOLD, help="Binary threshold for contour extraction")
    parser.add_argument("--image-ext", default=DEFAULT_IMAGE_EXT, help="Preferred image extension when metadata file names have no suffix")
    parser.add_argument("--center-mode", choices=["fixed", "image"], default="fixed", help="Reference center for offset calculation")
    parser.add_argument("--strict-image-shape", action="store_true", help="Stop if an image shape differs from --npix-x/--npix-y")
    parser.add_argument("--no-save-meta-csv", action="store_true", help="Do not write metadata CSV")
    parser.add_argument("--no-save-mark-images", action="store_true", help="Do not save centroid-marked images")
    parser.add_argument("--dpi", type=int, default=200, help="DPI for saved figures")
    parser.add_argument("--verbose", action="store_true", help="Print progress")
    return parser.parse_args(argv)


def build_config(args: argparse.Namespace) -> Config:
    datadir = Path(args.datadir)
    pic_dir = Path(args.pic_dir) if args.pic_dir else datadir / args.datetime

    if args.dbname:
        dbname = args.dbname
    elif args.datadate:
        dbname = f"necst_opticalpointing_{args.datadate}"
    else:
        dbname = None

    output_dir = Path(args.output_dir) if args.output_dir else (datadir / args.datetime / "analysis_cli")
    figs_dir = Path(args.figs_dir) if args.figs_dir else (output_dir / "figs")
    marks_dir = Path(args.marks_dir) if args.marks_dir else (output_dir / "marks")

    return Config(
        datetime_tag=args.datetime,
        datadir=datadir,
        datadate=args.datadate,
        dbname=dbname,
        meta_csv=Path(args.meta_csv) if args.meta_csv else None,
        pic_dir=pic_dir,
        output_dir=output_dir,
        figs_dir=figs_dir,
        marks_dir=marks_dir,
        npix_x=args.npix_x,
        npix_y=args.npix_y,
        pixel_scale_arcsec=args.pixel_scale_arcsec,
        threshold=args.threshold,
        image_ext=normalize_image_ext(args.image_ext),
        db_table=args.db_table,
        save_meta_csv=not args.no_save_meta_csv,
        save_mark_images=not args.no_save_mark_images,
        dpi=args.dpi,
        verbose=args.verbose,
        strict_image_shape=args.strict_image_shape,
        center_mode=args.center_mode,
    )


def ensure_dirs(cfg: Config) -> None:
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.figs_dir.mkdir(parents=True, exist_ok=True)
    if cfg.save_mark_images:
        cfg.marks_dir.mkdir(parents=True, exist_ok=True)


def log(msg: str, cfg: Config) -> None:
    if cfg.verbose:
        print(msg)


def decode_maybe_bytes(value: Any) -> str:
    if isinstance(value, (bytes, bytearray)):
        return value.decode(errors="replace")
    return str(value)


def normalize_image_ext(ext: str) -> str:
    ext = str(ext).strip()
    if ext == "":
        return DEFAULT_IMAGE_EXT
    if not ext.startswith("."):
        ext = "." + ext
    return ext


def normalize_file_cell(value: Any) -> str:
    """Normalize metadata `file` cell to a clean stem-or-name string."""
    text = decode_maybe_bytes(value).strip()
    if text == "" or text.lower() == "nan":
        return ""
    return text


def normalize_metadata(meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()

    required = ["cap_az", "cap_el", "file"]
    optional = ["ra", "dec", "cap_time"]

    for col in required:
        if col not in meta.columns:
            raise ValueError(f"Metadata is missing required column: {col}")

    for col in optional:
        if col not in meta.columns:
            meta[col] = np.nan

    meta["cap_az"] = pd.to_numeric(meta["cap_az"], errors="coerce")
    meta["cap_el"] = pd.to_numeric(meta["cap_el"], errors="coerce")
    meta["file"] = meta["file"].map(normalize_file_cell)

    if "ra" in meta.columns:
        meta["ra"] = pd.to_numeric(meta["ra"], errors="coerce")
    if "dec" in meta.columns:
        meta["dec"] = pd.to_numeric(meta["dec"], errors="coerce")

    if meta["file"].eq("").any():
        bad_rows = meta.index[meta["file"].eq("")].tolist()
        raise ValueError(f"Metadata contains empty file names at rows: {bad_rows[:10]}")

    if meta["cap_az"].isna().any() or meta["cap_el"].isna().any():
        bad_rows = meta.index[meta["cap_az"].isna() | meta["cap_el"].isna()].tolist()
        raise ValueError(f"Metadata contains invalid cap_az/cap_el at rows: {bad_rows[:10]}")

    return meta.reset_index(drop=True)


def load_metadata_from_db(cfg: Config) -> pd.DataFrame:
    if cfg.dbname is None:
        raise ValueError("DB mode requires either --datadate or --dbname.")

    try:
        import necstdb  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "Failed to import necstdb. Use --meta-csv or install necstdb in this environment."
        ) from exc

    db_path = cfg.datadir / cfg.dbname
    log(f"Opening DB: {db_path}", cfg)
    db = necstdb.opendb(str(db_path))
    df_meta = db.open_table(cfg.db_table).read(astype="pandas")

    cap_az: list[float] = []
    cap_el: list[float] = []
    ra: list[float] = []
    dec: list[float] = []
    cap_time: list[Any] = []
    file_name: list[str] = []

    for i in range(len(df_meta)):
        row_data = df_meta["data"].iloc[i]
        cap_time.append(df_meta["time"].iloc[i])
        cap_az.append(float(row_data[0]))
        cap_el.append(float(row_data[1]))
        ra.append(float(row_data[2]))
        dec.append(float(row_data[3]))
        name = decode_maybe_bytes(df_meta["id"].iloc[i]) + decode_maybe_bytes(df_meta["position"].iloc[i])
        file_name.append(name.replace(" ", ""))

    return pd.DataFrame(
        {
            "cap_az": cap_az,
            "cap_el": cap_el,
            "ra": ra,
            "dec": dec,
            "cap_time": cap_time,
            "file": file_name,
        }
    )


def load_metadata(cfg: Config) -> pd.DataFrame:
    if cfg.meta_csv is not None:
        log(f"Reading metadata CSV: {cfg.meta_csv}", cfg)
        meta = pd.read_csv(cfg.meta_csv)
    else:
        meta = load_metadata_from_db(cfg)

    return normalize_metadata(meta)


def save_metadata_csv(meta: pd.DataFrame, cfg: Config) -> Path:
    out = cfg.output_dir / f"necst_opticalpointing_{cfg.datetime_tag}_meta.csv"
    meta.to_csv(out, index=False)
    return out


def find_contours_compat(binary_image: np.ndarray) -> list[np.ndarray]:
    """OpenCV 3/4 compatible contour extraction."""
    res = cv2.findContours(binary_image, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    if len(res) == 2:
        contours, _hier = res
    elif len(res) == 3:
        _img, contours, _hier = res
    else:
        raise RuntimeError(f"Unexpected cv2.findContours return length: {len(res)}")
    return list(contours)


def choose_largest_contour(contours: Iterable[np.ndarray]) -> tuple[float, float, float]:
    stars: list[tuple[float, float]] = []
    areas: list[float] = []

    for cnt in contours:
        if cnt is None or len(cnt) == 0:
            continue
        M = cv2.moments(cnt)
        if M["m00"] != 0:
            cx = float(M["m10"] / M["m00"])
            cy = float(M["m01"] / M["m00"])
        else:
            cx, cy = map(float, cnt[0][0])
        area = float(cv2.contourArea(cnt))
        stars.append((cx, cy))
        areas.append(area)

    if not areas:
        raise ValueError("no contours")

    idx = int(np.argmax(areas))
    return stars[idx][0], stars[idx][1], areas[idx]


def candidate_image_names(file_value: str, preferred_ext: str) -> list[str]:
    p = Path(file_value)
    names: list[str] = []
    ext0 = normalize_image_ext(preferred_ext)
    exts = [ext0, ext0.lower(), ext0.upper(), ".JPG", ".jpg", ".JPEG", ".jpeg", ".PNG", ".png"]
    seen: set[str] = set()

    def add(name: str) -> None:
        if name not in seen:
            names.append(name)
            seen.add(name)

    if p.suffix:
        add(p.name)
        add(p.stem + p.suffix.lower())
        add(p.stem + p.suffix.upper())
    else:
        for ext in exts:
            add(file_value + ext)

    return names


def find_image_by_casefold_stem(pic_dir: Path, file_value: str) -> Path | None:
    target = Path(file_value).stem.casefold()
    target_suffix = Path(file_value).suffix.casefold()
    supported_exts = {".jpg", ".jpeg", ".png"}
    for path in pic_dir.iterdir():
        if not path.is_file():
            continue
        if path.stem.casefold() != target:
            continue
        if target_suffix and path.suffix.casefold() != target_suffix:
            continue
        if not target_suffix and path.suffix.casefold() not in supported_exts:
            continue
        return path
    return None


def safe_output_name(text: str) -> str:
    value = decode_maybe_bytes(text).strip()
    if value == "":
        return "unnamed"
    for ch in ("/", "\\", ":", "\n", "\r", "\t"):
        value = value.replace(ch, "_")
    return value


def resolve_image_path(pic_dir: Path, file_value: str, preferred_ext: str) -> Path:
    raw_path = Path(file_value)
    if raw_path.is_absolute() and raw_path.exists():
        return raw_path
    if raw_path.exists():
        return raw_path

    candidates = candidate_image_names(file_value, preferred_ext)
    for name in candidates:
        path = pic_dir / name
        if path.exists():
            return path

    fallback = find_image_by_casefold_stem(pic_dir, file_value)
    if fallback is not None:
        return fallback

    return pic_dir / candidates[0]


def get_reference_center(image_width: int, image_height: int, cfg: Config) -> tuple[float, float]:
    if cfg.center_mode == "image":
        return image_width / 2.0, image_height / 2.0
    return cfg.npix_x / 2.0, cfg.npix_y / 2.0


def annotate_detection(
    image_color_flipped: np.ndarray,
    x_center: float,
    y_center: float,
    title: str,
    output_path: Path,
    cfg: Config,
    detected_x: float | None = None,
    detected_y: float | None = None,
    message: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.imshow(image_color_flipped, vmin=0, vmax=256)
    ax.set_xlim(0, image_color_flipped.shape[1])
    ax.set_ylim(0, image_color_flipped.shape[0])
    ax.axvline(x_center, color="white", linestyle="--")
    ax.axhline(y_center, color="white", linestyle="--")
    if detected_x is not None and detected_y is not None:
        ax.plot(detected_x, detected_y, "+", color="red", markersize=14)
    if message:
        ax.text(100, 100, message, color="yellow" if "NOT FOUND" in message else "red", fontsize=14)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=cfg.dpi)
    plt.close(fig)


def analyze_images(meta: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    rows: list[DetectionResult] = []

    for i, row in meta.iterrows():
        file_value = str(row["file"])
        image_path = resolve_image_path(cfg.pic_dir, file_value, cfg.image_ext)
        cap_az = float(row["cap_az"])
        cap_el = float(row["cap_el"])
        ra = None if pd.isna(row.get("ra", np.nan)) else float(row["ra"])
        dec = None if pd.isna(row.get("dec", np.nan)) else float(row["dec"])
        cap_time = row.get("cap_time", np.nan)

        log(f"[{i + 1}/{len(meta)}] {image_path.name}", cfg)

        detected = False
        reason = ""
        pix_x = np.nan
        pix_y = np.nan
        contour_area = np.nan
        image_width = float(cfg.npix_x)
        image_height = float(cfg.npix_y)
        image_for_plot = np.zeros((cfg.npix_y, cfg.npix_x, 3), dtype=np.uint8)

        img_color = None
        img_gray = None
        if image_path.exists():
            img_color = cv2.imread(str(image_path))
            img_gray = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)

        if img_color is None or img_gray is None:
            reason = "image not found"
            if image_path.exists():
                reason = "image load failed"
            x_center_ref = cfg.npix_x / 2.0
            y_center_ref = cfg.npix_y / 2.0
        else:
            image_height = float(img_gray.shape[0])
            image_width = float(img_gray.shape[1])

            if int(image_width) != cfg.npix_x or int(image_height) != cfg.npix_y:
                shape_msg = (
                    f"image shape mismatch: got {int(image_width)}x{int(image_height)}, "
                    f"expected {cfg.npix_x}x{cfg.npix_y}"
                )
                if cfg.strict_image_shape:
                    raise ValueError(f"{image_path}: {shape_msg}")
                reason = shape_msg
            else:
                reason = "ok"

            x_center_ref, y_center_ref = get_reference_center(int(image_width), int(image_height), cfg)
            image_for_plot = np.flipud(img_color)
            img_gray_flipped = np.flipud(img_gray)

            try:
                _thr, binary = cv2.threshold(img_gray_flipped, cfg.threshold, 255, cv2.THRESH_BINARY)
                contours = find_contours_compat(binary)
                pix_x, pix_y, contour_area = choose_largest_contour(contours)
                detected = True
                if reason == "ok":
                    reason = "ok"
                else:
                    reason = reason + "; detection ok"
            except Exception as exc:
                if reason == "ok":
                    reason = f"detection failed: {exc}"
                else:
                    reason = reason + f"; detection failed: {exc}"

        dpix_x = pix_x - x_center_ref if np.isfinite(pix_x) else np.nan
        dpix_y = pix_y - y_center_ref if np.isfinite(pix_y) else np.nan
        dx_arcsec = dpix_x * cfg.pixel_scale_arcsec if np.isfinite(dpix_x) else np.nan
        dy_arcsec = dpix_y * cfg.pixel_scale_arcsec if np.isfinite(dpix_y) else np.nan

        rows.append(
            DetectionResult(
                file=file_value,
                image_path=str(image_path),
                cap_az=cap_az,
                cap_el=cap_el,
                ra=ra,
                dec=dec,
                cap_time=cap_time,
                detected=detected,
                reason=reason,
                image_width=image_width,
                image_height=image_height,
                x_center_ref=float(x_center_ref),
                y_center_ref=float(y_center_ref),
                pix_x=float(pix_x) if np.isfinite(pix_x) else np.nan,
                pix_y=float(pix_y) if np.isfinite(pix_y) else np.nan,
                contour_area=float(contour_area) if np.isfinite(contour_area) else np.nan,
                dpix_x=float(dpix_x) if np.isfinite(dpix_x) else np.nan,
                dpix_y=float(dpix_y) if np.isfinite(dpix_y) else np.nan,
                dx_arcsec=float(dx_arcsec) if np.isfinite(dx_arcsec) else np.nan,
                dy_arcsec=float(dy_arcsec) if np.isfinite(dy_arcsec) else np.nan,
            )
        )

        if cfg.save_mark_images:
            msg = None
            if "image not found" in reason:
                msg = "IMAGE NOT FOUND"
            elif "detection failed" in reason:
                msg = "DETECTION ERROR"
            annotate_detection(
                image_color_flipped=image_for_plot,
                x_center=x_center_ref,
                y_center=y_center_ref,
                title=f"Az : {cap_az} EL : {cap_el}",
                output_path=cfg.marks_dir / f"{safe_output_name(file_value)}.mark.png",
                cfg=cfg,
                detected_x=pix_x if detected else None,
                detected_y=pix_y if detected else None,
                message=msg,
            )

    return pd.DataFrame([asdict(r) for r in rows])


def finite_rms(values: np.ndarray) -> float:
    finite = np.isfinite(values)
    if not np.any(finite):
        return float("nan")
    return float(np.sqrt(np.mean(values[finite] ** 2)))


def finite_sigma(values: np.ndarray) -> float:
    finite = np.isfinite(values)
    if not np.any(finite):
        return float("nan")
    return float(np.nanstd(values[finite]))


def scatter_plot(
    x: np.ndarray,
    y: np.ndarray,
    xlabel: tuple[str, str],
    ylabel: tuple[str, str],
    d_rms: float,
    output_path: Path,
    circle_radius: float = 5.0,
) -> None:
    fig, ax = plt.subplots()
    ax.scatter(x, y, s=10)

    if (xlabel[0] == "dAz" and ylabel[0] == "dEl") or (xlabel[0] == "zansa_dAz" and ylabel[0] == "zansa_dEl"):
        ax.set_title(f"{xlabel[0]}_vs_{ylabel[0]}\nrms = {d_rms:0.2f}[arcsec]")
        ax.set_aspect("equal", "datalim")
        deg = np.linspace(-180, 180, 360)
        X = [circle_radius * math.sin(math.radians(num)) for num in deg]
        Y = [circle_radius * math.cos(math.radians(num)) for num in deg]
        ax.plot(X, Y, "r")
    else:
        ax.set_title(f"{xlabel[0]}_vs_{ylabel[0]}")

    ax.set_xlabel(f"{xlabel[0]} [{xlabel[1]}]")
    ax.set_ylabel(f"{ylabel[0]} [{ylabel[1]}]")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def create_summary_plots(results: pd.DataFrame, cfg: Config) -> dict[str, Path]:
    az = results["cap_az"].to_numpy(dtype=float)
    el = results["cap_el"].to_numpy(dtype=float)
    dx = results["dx_arcsec"].to_numpy(dtype=float)
    dy = results["dy_arcsec"].to_numpy(dtype=float)

    d_x_rms = finite_rms(dx)
    d_y_rms = finite_rms(dy)
    d_rms = float(np.sqrt(d_x_rms ** 2 + d_y_rms ** 2)) if np.isfinite(d_x_rms) and np.isfinite(d_y_rms) else float("nan")

    outputs = {
        "az_vs_el": cfg.figs_dir / "Az_vs_El.png",
        "az_vs_dx": cfg.figs_dir / "Az_vs_dAz.png",
        "az_vs_dy": cfg.figs_dir / "Az_vs_dEl.png",
        "el_vs_dx": cfg.figs_dir / "El_vs_dAz.png",
        "el_vs_dy": cfg.figs_dir / "El_vs_dEl.png",
        "dx_vs_dy": cfg.figs_dir / "dAz_vs_dEl.png",
        "summary": cfg.figs_dir / "optical_pointing_results.png",
    }

    scatter_plot(az, el, ("Az", "degree"), ("El", "degree"), d_rms, outputs["az_vs_el"])
    scatter_plot(az, dx, ("Az", "degree"), ("dAz", "arcsec"), d_rms, outputs["az_vs_dx"])
    scatter_plot(az, dy, ("Az", "degree"), ("dEl", "arcsec"), d_rms, outputs["az_vs_dy"])
    scatter_plot(el, dx, ("El", "degree"), ("dAz", "arcsec"), d_rms, outputs["el_vs_dx"])
    scatter_plot(el, dy, ("El", "degree"), ("dEl", "arcsec"), d_rms, outputs["el_vs_dy"])
    scatter_plot(dx, dy, ("dAz", "arcsec"), ("dEl", "arcsec"), d_rms, outputs["dx_vs_dy"], circle_radius=10.0)

    fig = plt.figure(figsize=(15, 8))
    ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]

    ax[0].plot(az, el, "o")
    ax[0].set_xlabel("Az [deg]")
    ax[0].set_ylabel("El [deg]")

    ax[1].plot(az, dx, "o")
    ax[1].set_xlabel("Az [deg]")
    ax[1].set_ylabel("dx [arcsec]")

    ax[2].plot(el, dx, "o")
    ax[2].set_xlabel("El [deg]")
    ax[2].set_ylabel("dx [arcsec]")

    ax[3].plot(dx, dy, "o")
    ax[3].set_xlabel("dx [arcsec]")
    ax[3].set_ylabel("dy [arcsec]")
    deg = np.linspace(-180, 180, 360)
    X = [10.0 * math.sin(math.radians(num)) for num in deg]
    Y = [10.0 * math.cos(math.radians(num)) for num in deg]
    ax[3].plot(X, Y, "r")

    ax[4].plot(az, dy, "o")
    ax[4].set_xlabel("Az [deg]")
    ax[4].set_ylabel("dy [arcsec]")

    ax[5].plot(el, dy, "o")
    ax[5].set_xlabel("El [deg]")
    ax[5].set_ylabel("dy [arcsec]")

    for a in ax:
        a.grid(True)

    fig.suptitle("optical_pointing")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outputs["summary"])
    plt.close(fig)

    return outputs


def save_results(results: pd.DataFrame, cfg: Config) -> dict[str, Path]:
    legacy_dat = cfg.output_dir / "Az_El_dAz_dEl.dat"
    finite_results = results[["cap_az", "cap_el", "dx_arcsec", "dy_arcsec"]].to_numpy(dtype=float)
    np.savetxt(legacy_dat, finite_results, delimiter=", ")

    detailed_csv = cfg.output_dir / f"{cfg.datetime_tag}_optical_pointing_results_detailed.csv"
    results.to_csv(detailed_csv, index=False)

    formatted = pd.DataFrame(
        {
            "dt_cross_id": [f"c{i:03d}" for i in range(len(results))],
            "cross_id": range(len(results)),
            "Az": results["cap_az"],
            "El": results["cap_el"],
            "dx": results["dx_arcsec"],
            "dy": results["dy_arcsec"],
        }
    )
    formatted_csv = cfg.output_dir / f"{cfg.datetime_tag}_optical_pointing_results_formatted.csv"
    formatted.to_csv(formatted_csv, index=False)

    return {
        "legacy_dat": legacy_dat,
        "detailed_csv": detailed_csv,
        "formatted_csv": formatted_csv,
    }


def build_summary(meta: pd.DataFrame, results: pd.DataFrame, cfg: Config) -> dict[str, Any]:
    dx = results["dx_arcsec"].to_numpy(dtype=float)
    dy = results["dy_arcsec"].to_numpy(dtype=float)
    ok = results["detected"].to_numpy(dtype=bool)

    reasons = results["reason"].fillna("").astype(str)
    n_missing = int(reasons.str.contains("image not found", regex=False).sum())
    n_load_failed = int(reasons.str.contains("image load failed", regex=False).sum())
    n_detection_failed = int(reasons.str.contains("detection failed", regex=False).sum())
    n_shape_mismatch = int(reasons.str.contains("image shape mismatch", regex=False).sum())

    dx_rms = finite_rms(dx)
    dy_rms = finite_rms(dy)
    dx_sigma = finite_sigma(dx)
    dy_sigma = finite_sigma(dy)

    summary = {
        "n_total": int(len(meta)),
        "n_detected": int(np.count_nonzero(ok)),
        "n_failed": int(len(meta) - np.count_nonzero(ok)),
        "n_missing_image": n_missing,
        "n_image_load_failed": n_load_failed,
        "n_detection_failed": n_detection_failed,
        "n_shape_mismatch": n_shape_mismatch,
        "dx_rms_arcsec": dx_rms,
        "dy_rms_arcsec": dy_rms,
        "d_rms_arcsec": float(np.sqrt(dx_rms ** 2 + dy_rms ** 2)) if np.isfinite(dx_rms) and np.isfinite(dy_rms) else float("nan"),
        "dx_sigma_arcsec": dx_sigma,
        "dy_sigma_arcsec": dy_sigma,
        "d_sigma_arcsec": float(np.sqrt(dx_sigma ** 2 + dy_sigma ** 2)) if np.isfinite(dx_sigma) and np.isfinite(dy_sigma) else float("nan"),
        "pixel_scale_arcsec": cfg.pixel_scale_arcsec,
        "threshold": cfg.threshold,
        "npix_x": cfg.npix_x,
        "npix_y": cfg.npix_y,
        "center_mode": cfg.center_mode,
        "strict_image_shape": cfg.strict_image_shape,
        "pic_dir": str(cfg.pic_dir),
        "output_dir": str(cfg.output_dir),
    }
    return summary


def write_summary_json(summary: dict[str, Any], cfg: Config) -> Path:
    path = cfg.output_dir / f"{cfg.datetime_tag}_summary.json"
    with path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    return path


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cfg = build_config(args)

    try:
        ensure_dirs(cfg)

        meta = load_metadata(cfg)
        meta_csv_path = None
        if cfg.save_meta_csv:
            meta_csv_path = save_metadata_csv(meta, cfg)

        results = analyze_images(meta, cfg)
        result_paths = save_results(results, cfg)
        plot_paths = create_summary_plots(results, cfg)
        summary = build_summary(meta, results, cfg)
        summary_path = write_summary_json(summary, cfg)
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1

    print(f"Processed {summary['n_total']} frames")
    print(f"Detected   {summary['n_detected']} frames")
    print(f"Failed     {summary['n_failed']} frames")
    print(f"Missing    {summary['n_missing_image']} images")
    print(f"LoadFail   {summary['n_image_load_failed']} frames")
    print(f"DetFail    {summary['n_detection_failed']} frames")
    print(f"ShapeMis   {summary['n_shape_mismatch']} frames")
    print(f"rms        {summary['d_rms_arcsec']:.2f} [arcsec]")
    print(f"sigma      {summary['d_sigma_arcsec']:.2f} [arcsec]")
    if meta_csv_path is not None:
        print(f"meta_csv   {meta_csv_path}")
    for name, path in {**result_paths, **plot_paths, "summary_json": summary_path}.items():
        print(f"{name:24s} {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
