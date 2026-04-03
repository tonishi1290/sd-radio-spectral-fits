from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from .io_utils import (
    ARCSEC_PER_DEG,
    DEG_PER_RAD,
    load_manifest,
    parse_pointing_param_toml,
    read_csv,
    sha256_file,
    write_csv,
)
from .models import resolve_model


@dataclass
class NormalizeOptions:
    input_csv: str | None = None
    manifest: str | None = None
    output_csv: str = "normalized.csv"
    pair_mode: str = "union"
    dataset_id: str | None = None
    used_param: str | None = None
    telescope: str | None = None
    model: str = "omu1p85m_combined_v1"
    angle_unit: str = "deg"
    offset_unit: str = "arcsec"
    measurement_space: str = "star_offset"
    positive_az_moves_star: str | None = None
    positive_el_moves_star: str | None = None
    command_err_mode: str | None = None


CANONICAL_COLUMNS = [
    "dataset_id",
    "dt_cross_id",
    "cross_id",
    "az_deg",
    "el_deg",
    "input_angle_unit",
    "raw_dx_input_arcsec",
    "raw_dy_input_arcsec",
    "measured_dx_arcsec",
    "measured_dy_arcsec",
    "delta_err_dx_arcsec",
    "delta_err_dy_arcsec",
    "source_dir",
    "used_param_file",
    "used_param_hash",
    "telescope",
    "model_name",
    "measurement_space",
    "positive_az_moves_star",
    "positive_el_moves_star",
    "command_err_mode",
    "dx_to_model_sign",
    "dy_to_model_sign",
    "applied_model_dx_arcsec",
    "applied_model_dy_arcsec",
    "used_model_dx_arcsec",
    "used_model_dy_arcsec",
    "absolute_target_rule",
]

RAW_COLUMN_CANDIDATES = {
    "az_deg": ["az_deg", "Az", "az", "AZ", "az_rad", "Az_rad", "AZ_RAD"],
    "el_deg": ["el_deg", "El", "el", "EL", "el_rad", "El_rad", "EL_RAD"],
    "measured_dx_arcsec": [
        "measured_dx_arcsec",
        "dx_arcsec",
        "dx_deg",
        "dx_rad",
        "dx",
        "dX",
        "DX",
        "dAz_arcsec",
        "dAz_deg",
        "dAz_rad",
        "dAz",
        "daz_arcsec",
        "daz_deg",
        "daz_rad",
        "daz",
    ],
    "measured_dy_arcsec": [
        "measured_dy_arcsec",
        "dy_arcsec",
        "dy_deg",
        "dy_rad",
        "dy",
        "dY",
        "DY",
        "dEl_arcsec",
        "dEl_deg",
        "dEl_rad",
        "dEl",
        "del_arcsec",
        "del_deg",
        "del_rad",
        "del",
    ],
    "dt_cross_id": ["dt_cross_id"],
    "cross_id": ["cross_id"],
    "source_dir": ["source_dir"],
}

_ALLOWED_MEASUREMENT_SPACE = {"model_delta", "command_offset", "star_offset"}
_ALLOWED_STAR_X_MOTION = {"left", "right"}
_ALLOWED_STAR_Y_MOTION = {"up", "down"}
_ALLOWED_COMMAND_ERR_MODE = {"subtract", "add"}


def first_existing(df: pd.DataFrame, names: list[str]) -> str | None:
    for n in names:
        if n in df.columns:
            return n
    return None


def _convert_angle_to_deg(series: pd.Series, src_name: str, angle_unit: str) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce").astype(float)
    src = str(src_name).strip().lower()
    if src.endswith("rad") or src.endswith("_rad") or src == "az_rad" or src == "el_rad":
        return vals * DEG_PER_RAD
    if src.endswith("deg") or src.endswith("_deg"):
        return vals
    if angle_unit == "deg":
        return vals
    if angle_unit == "rad":
        return vals * DEG_PER_RAD
    raise ValueError("angle_unit must be 'deg' or 'rad'.")


def _detect_explicit_angle_unit(src_name: str) -> str | None:
    src = str(src_name).strip().lower()
    if src.endswith("rad") or src.endswith("_rad") or src in {"az_rad", "el_rad"}:
        return "rad"
    if src.endswith("deg") or src.endswith("_deg"):
        return "deg"
    return None


def _resolve_input_angle_unit(az_src_name: str, el_src_name: str, fallback_angle_unit: str) -> str:
    az_unit = _detect_explicit_angle_unit(az_src_name) or fallback_angle_unit
    el_unit = _detect_explicit_angle_unit(el_src_name) or fallback_angle_unit
    if az_unit != el_unit:
        raise ValueError(
            f"Az/El input units are inconsistent: az source '{az_src_name}' resolves to {az_unit}, "
            f"but el source '{el_src_name}' resolves to {el_unit}."
        )
    return az_unit


def _convert_offset_to_arcsec(series: pd.Series, src_name: str, offset_unit: str) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce").astype(float)
    src = str(src_name).strip().lower()
    if "arcsec" in src:
        return vals
    if src.endswith("deg") or src.endswith("_deg"):
        return vals * ARCSEC_PER_DEG
    if src.endswith("rad") or src.endswith("_rad"):
        return vals * DEG_PER_RAD * ARCSEC_PER_DEG

    if offset_unit == "arcsec":
        warnings.warn(
            f"Offset column '{src_name}' has no unit in its name; assuming arcsec via offset_unit='arcsec'. "
            "Use an explicit name such as dx_arcsec/dy_arcsec, dx_deg/dy_deg, or dx_rad/dy_rad to avoid ambiguity.",
            UserWarning,
            stacklevel=2,
        )
        return vals
    if offset_unit == "deg":
        warnings.warn(
            f"Offset column '{src_name}' has no unit in its name; assuming degrees via offset_unit='deg'. "
            "Use an explicit name such as dx_arcsec/dy_arcsec, dx_deg/dy_deg, or dx_rad/dy_rad to avoid ambiguity.",
            UserWarning,
            stacklevel=2,
        )
        return vals * ARCSEC_PER_DEG
    if offset_unit == "rad":
        warnings.warn(
            f"Offset column '{src_name}' has no unit in its name; assuming radians via offset_unit='rad'. "
            "Use an explicit name such as dx_arcsec/dy_arcsec, dx_deg/dy_deg, or dx_rad/dy_rad to avoid ambiguity.",
            UserWarning,
            stacklevel=2,
        )
        return vals * DEG_PER_RAD * ARCSEC_PER_DEG
    raise ValueError("offset_unit must be 'arcsec', 'deg', or 'rad'.")


def standardize_raw_table(df: pd.DataFrame, angle_unit: str, offset_unit: str) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    source_names: dict[str, str] = {}
    n = len(df)
    for canonical, names in RAW_COLUMN_CANDIDATES.items():
        src = first_existing(df, names)
        if src is None:
            if canonical == "source_dir":
                continue
            if canonical == "dt_cross_id":
                out[canonical] = np.arange(n, dtype=int)
                source_names[canonical] = "<auto:index>"
                continue
            if canonical == "cross_id":
                if "dt_cross_id" in out:
                    out[canonical] = out["dt_cross_id"]
                    source_names[canonical] = "<auto:dt_cross_id>"
                else:
                    out[canonical] = np.arange(n, dtype=int)
                    source_names[canonical] = "<auto:index>"
                continue
            raise KeyError(f"Could not find a source column for {canonical}. Available columns: {list(df.columns)}")
        source_names[canonical] = src
        out[canonical] = df[src]
    out["az_deg"] = _convert_angle_to_deg(out["az_deg"], source_names["az_deg"], angle_unit)
    out["el_deg"] = _convert_angle_to_deg(out["el_deg"], source_names["el_deg"], angle_unit)
    out["measured_dx_arcsec"] = _convert_offset_to_arcsec(out["measured_dx_arcsec"], source_names["measured_dx_arcsec"], offset_unit)
    out["measured_dy_arcsec"] = _convert_offset_to_arcsec(out["measured_dy_arcsec"], source_names["measured_dy_arcsec"], offset_unit)
    out["raw_dx_input_arcsec"] = out["measured_dx_arcsec"]
    out["raw_dy_input_arcsec"] = out["measured_dy_arcsec"]
    if "source_dir" not in out:
        out["source_dir"] = ""
    out["input_angle_unit"] = _resolve_input_angle_unit(source_names["az_deg"], source_names["el_deg"], angle_unit)
    return out


def build_union(df: pd.DataFrame) -> pd.DataFrame:
    xpart = df.dropna(subset=["measured_dx_arcsec"]).copy()
    ypart = df.dropna(subset=["measured_dy_arcsec"]).copy()
    merged = pd.merge(
        xpart[["dt_cross_id", "cross_id", "az_deg", "el_deg", "input_angle_unit", "measured_dx_arcsec", "raw_dx_input_arcsec", "source_dir"]],
        ypart[["dt_cross_id", "measured_dy_arcsec", "raw_dy_input_arcsec"]],
        on="dt_cross_id",
        how="inner",
    )
    return merged


def _normalize_convention_settings(
    *,
    measurement_space: str,
    positive_az_moves_star: str | None,
    positive_el_moves_star: str | None,
    command_err_mode: str | None,
) -> tuple[str, str, str, int, int]:
    ms = str(measurement_space or "star_offset").strip()
    if ms not in _ALLOWED_MEASUREMENT_SPACE:
        raise ValueError(f"measurement_space must be one of {sorted(_ALLOWED_MEASUREMENT_SPACE)}")

    px = ""
    py = ""
    cem = ""
    if ms == "model_delta":
        if positive_az_moves_star or positive_el_moves_star or command_err_mode:
            warnings.warn(
                "measurement_space='model_delta' ignores positive_az_moves_star, positive_el_moves_star, and command_err_mode.",
                UserWarning,
                stacklevel=2,
            )
        return ms, px, py, 1, 1

    cem = str(command_err_mode or "").strip()
    if cem not in _ALLOWED_COMMAND_ERR_MODE:
        raise ValueError("command_err_mode must be provided as 'subtract' or 'add' when measurement_space is not 'model_delta'.")
    err_from_command = -1 if cem == "subtract" else 1

    if ms == "command_offset":
        if positive_az_moves_star or positive_el_moves_star:
            raise ValueError("positive_az_moves_star and positive_el_moves_star must not be provided when measurement_space='command_offset'.")
        return ms, px, py, err_from_command, err_from_command

    px = str(positive_az_moves_star or "").strip()
    py = str(positive_el_moves_star or "").strip()
    if px not in _ALLOWED_STAR_X_MOTION:
        raise ValueError("positive_az_moves_star must be 'left' or 'right' when measurement_space='star_offset'.")
    if py not in _ALLOWED_STAR_Y_MOTION:
        raise ValueError("positive_el_moves_star must be 'up' or 'down' when measurement_space='star_offset'.")

    command_from_star_x = 1 if px == "left" else -1
    command_from_star_y = 1 if py == "down" else -1
    dx_to_model_sign = err_from_command * command_from_star_x
    dy_to_model_sign = err_from_command * command_from_star_y
    return ms, px, py, dx_to_model_sign, dy_to_model_sign


def apply_measurement_convention(
    df: pd.DataFrame,
    *,
    measurement_space: str,
    positive_az_moves_star: str | None,
    positive_el_moves_star: str | None,
    command_err_mode: str | None,
) -> pd.DataFrame:
    ms, px, py, dx_to_model_sign, dy_to_model_sign = _normalize_convention_settings(
        measurement_space=measurement_space,
        positive_az_moves_star=positive_az_moves_star,
        positive_el_moves_star=positive_el_moves_star,
        command_err_mode=command_err_mode,
    )
    out = df.copy()
    out["measured_dx_arcsec"] = pd.to_numeric(out["measured_dx_arcsec"], errors="coerce") * float(dx_to_model_sign)
    out["measured_dy_arcsec"] = pd.to_numeric(out["measured_dy_arcsec"], errors="coerce") * float(dy_to_model_sign)
    out["delta_err_dx_arcsec"] = out["measured_dx_arcsec"]
    out["delta_err_dy_arcsec"] = out["measured_dy_arcsec"]
    out["measurement_space"] = ms
    out["positive_az_moves_star"] = px
    out["positive_el_moves_star"] = py
    out["command_err_mode"] = command_err_mode or ""
    out["dx_to_model_sign"] = int(dx_to_model_sign)
    out["dy_to_model_sign"] = int(dy_to_model_sign)
    out["absolute_target_rule"] = "used_model_plus_delta_err"
    return out


def annotate_used_model(df: pd.DataFrame, used_param: str | None, model_name: str) -> pd.DataFrame:
    out = df.copy()
    out["used_param_file"] = used_param or ""
    out["used_param_hash"] = ""
    out["applied_model_dx_arcsec"] = np.nan
    out["applied_model_dy_arcsec"] = np.nan
    out["used_model_dx_arcsec"] = np.nan
    out["used_model_dy_arcsec"] = np.nan
    if not used_param:
        return out
    _, params = parse_pointing_param_toml(used_param)
    model = resolve_model(model_name)
    out["used_param_hash"] = sha256_file(used_param)
    dx_deg = model.predict_dx_deg(params, out["az_deg"], out["el_deg"])
    dy_deg = model.predict_dy_deg(params, out["az_deg"], out["el_deg"])
    out["applied_model_dx_arcsec"] = np.asarray(dx_deg, dtype=float) * ARCSEC_PER_DEG
    out["applied_model_dy_arcsec"] = np.asarray(dy_deg, dtype=float) * ARCSEC_PER_DEG
    out["used_model_dx_arcsec"] = out["applied_model_dx_arcsec"]
    out["used_model_dy_arcsec"] = out["applied_model_dy_arcsec"]
    return out


def normalize_one(
    *,
    input_csv: str,
    output_rows: list[pd.DataFrame],
    dataset_id: str | None,
    used_param: str | None,
    telescope: str | None,
    model_name: str,
    pair_mode: str,
    angle_unit: str,
    offset_unit: str,
    measurement_space: str,
    positive_az_moves_star: str | None,
    positive_el_moves_star: str | None,
    command_err_mode: str | None,
) -> None:
    raw = read_csv(input_csv)
    base = standardize_raw_table(raw, angle_unit=angle_unit, offset_unit=offset_unit)
    if pair_mode == "union":
        base = build_union(base)
    elif pair_mode != "raw":
        raise ValueError("pair_mode must be 'union' or 'raw'.")
    base = apply_measurement_convention(
        base,
        measurement_space=measurement_space,
        positive_az_moves_star=positive_az_moves_star,
        positive_el_moves_star=positive_el_moves_star,
        command_err_mode=command_err_mode,
    )
    base["dataset_id"] = dataset_id or Path(input_csv).stem
    base["telescope"] = telescope or ""
    base["model_name"] = model_name
    base = annotate_used_model(base, used_param, model_name)
    string_like = {
        "dataset_id",
        "telescope",
        "model_name",
        "source_dir",
        "input_angle_unit",
        "used_param_file",
        "used_param_hash",
        "measurement_space",
        "positive_az_moves_star",
        "positive_el_moves_star",
        "command_err_mode",
        "absolute_target_rule",
    }
    integer_like = {"dx_to_model_sign", "dy_to_model_sign"}
    for col in CANONICAL_COLUMNS:
        if col not in base.columns:
            if col in string_like:
                base[col] = ""
            elif col in integer_like:
                base[col] = np.nan
            else:
                base[col] = np.nan
    output_rows.append(base[CANONICAL_COLUMNS])


def _reject_removed_manifest_keys(ds: dict) -> None:
    if "sign_convention" in ds:
        raise ValueError(
            "Manifest key 'sign_convention' has been removed. Normalized values are always interpreted as delta_err to be added to the used model; remove this key from the manifest."
        )


def run_normalize(opts: NormalizeOptions) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    if opts.manifest:
        for ds in load_manifest(opts.manifest):
            _reject_removed_manifest_keys(ds)
            normalize_one(
                input_csv=str(ds["input_csv"]),
                output_rows=rows,
                dataset_id=ds.get("dataset_id"),
                used_param=ds.get("used_param"),
                telescope=ds.get("telescope"),
                model_name=ds.get("model", opts.model),
                pair_mode=ds.get("pair_mode", opts.pair_mode),
                angle_unit=ds.get("angle_unit", opts.angle_unit),
                offset_unit=ds.get("offset_unit", opts.offset_unit),
                measurement_space=ds.get("measurement_space", opts.measurement_space),
                positive_az_moves_star=ds.get("positive_az_moves_star", opts.positive_az_moves_star),
                positive_el_moves_star=ds.get("positive_el_moves_star", opts.positive_el_moves_star),
                command_err_mode=ds.get("command_err_mode", opts.command_err_mode),
            )
    elif opts.input_csv:
        normalize_one(
            input_csv=opts.input_csv,
            output_rows=rows,
            dataset_id=opts.dataset_id,
            used_param=opts.used_param,
            telescope=opts.telescope,
            model_name=opts.model,
            pair_mode=opts.pair_mode,
            angle_unit=opts.angle_unit,
            offset_unit=opts.offset_unit,
            measurement_space=opts.measurement_space,
            positive_az_moves_star=opts.positive_az_moves_star,
            positive_el_moves_star=opts.positive_el_moves_star,
            command_err_mode=opts.command_err_mode,
        )
    else:
        raise ValueError("Either input_csv or manifest must be provided.")
    out = pd.concat(rows, ignore_index=True)
    write_csv(out, opts.output_csv)
    return out
