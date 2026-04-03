from __future__ import annotations

import csv
import hashlib
import json
import math
import tomllib

import tomlkit
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ARCSEC_PER_DEG = 3600.0
RAD_PER_DEG = math.pi / 180.0
DEG_PER_RAD = 180.0 / math.pi


def _preferred_angle_output_unit(param_name: str, native_unit: str, requested_unit: str) -> str:
    native_unit = (native_unit or "").strip()
    if native_unit != "deg":
        return native_unit
    if str(param_name).startswith("omega") or str(param_name) == "cor_p":
        return "deg"
    return requested_unit


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_toml(path: str | Path) -> dict[str, Any]:
    with open(path, "rb") as f:
        return tomllib.load(f)


def _parse_param_key(key: str) -> tuple[str, str | None]:
    key = str(key).strip()
    if "[" in key and key.endswith("]"):
        name, unit = key[:-1].split("[", 1)
        return name.strip(), unit.strip()
    return key, None


def parse_param_mapping(raw: dict[str, Any]) -> dict[str, float]:
    params: dict[str, float] = {}
    for key, val in raw.items():
        name, unit = _parse_param_key(key)
        fv = float(val)
        if unit == "arcsec":
            params[name] = fv / ARCSEC_PER_DEG
        elif unit == "rad":
            params[name] = fv * DEG_PER_RAD
        else:
            params[name] = fv
    return params


def parse_pointing_param_toml(path: str | Path) -> tuple[str | None, dict[str, float]]:
    doc = read_toml(path)
    model = doc.get("metadata", {}).get("model")
    raw = doc.get("pointing_params", {})
    return model, parse_param_mapping(raw)


def parse_param_key_styles(path: str | Path, section_names: tuple[str, ...] = ("pointing_params", "fixed_params")) -> dict[str, dict[str, Any]]:
    doc = read_toml(path)
    styles: dict[str, dict[str, Any]] = {}
    for sec in section_names:
        raw = doc.get(sec)
        if not isinstance(raw, dict):
            continue
        for raw_key in raw.keys():
            name, explicit_unit = _parse_param_key(str(raw_key))
            styles[name] = {
                "original_key": str(raw_key),
                "explicit_unit": explicit_unit,
            }
        if styles:
            return styles
    return styles


def parse_param_file_sections(path: str | Path, section_names: tuple[str, ...] = ("fixed_params", "pointing_params")) -> tuple[str | None, dict[str, float], str | None]:
    doc = read_toml(path)
    model = doc.get("metadata", {}).get("model")
    for sec in section_names:
        raw = doc.get(sec)
        if isinstance(raw, dict):
            return model, parse_param_mapping(raw), sec
    top_level = {k: v for k, v in doc.items() if k != "metadata" and not isinstance(v, dict)}
    if top_level:
        return model, parse_param_mapping(top_level), "top_level"
    return model, {}, None



def _resolve_param_defs_for_model(model_name: str) -> dict[str, str]:
    try:
        from .models import resolve_model  # local import to avoid circular import at module load time

        base_name = model_name.split("__", 1)[0]
        model = resolve_model(base_name)
        return {p.name: p.unit for p in model.parameter_defs()}
    except Exception:
        return {}


def _format_param_for_output(value: float, param_name: str, native_unit: str, requested_unit: str, explicit_unit: str | None = None) -> tuple[float, str | None]:
    native_unit = (native_unit or "").strip()
    if native_unit in {"", "1"}:
        return float(value), None
    if requested_unit == "preserve" and explicit_unit:
        actual_unit = explicit_unit
    else:
        actual_unit = _preferred_angle_output_unit(param_name, native_unit, requested_unit)
    if native_unit == "deg":
        if actual_unit == "arcsec":
            return float(value) * ARCSEC_PER_DEG, "arcsec"
        if actual_unit == "rad":
            return float(value) * RAD_PER_DEG, "rad"
        return float(value), "deg"
    # Mixed units such as arcsec/deg^n are kept in their native representation.
    return float(value), native_unit


def write_pointing_param_toml(
    path: str | Path,
    *,
    model_name: str,
    params_deg: dict[str, float],
    unit: str = "deg",
    metadata_extra: dict[str, Any] | None = None,
    section_name: str = "pointing_params",
    key_style_map: dict[str, dict[str, Any]] | None = None,
) -> None:
    p = Path(path)
    ensure_parent(p)
    meta = dict(metadata_extra or {})
    meta.setdefault("model", model_name)
    param_units = _resolve_param_defs_for_model(model_name)
    preferred_order = list(param_units) or ["a1", "a2", "a3", "b1", "b2", "b3", "g1", "c1", "c2", "d1", "d2", "e1", "e2"]
    with p.open("w", encoding="utf-8") as f:
        f.write("[metadata]\n")
        for k, v in meta.items():
            f.write(f"{k} = {toml_scalar(v)}\n")
        f.write(f"\n[{section_name}]\n")
        for key in preferred_order + sorted(set(params_deg) - set(preferred_order)):
            if key not in params_deg:
                continue
            val = float(params_deg[key])
            native_unit = param_units.get(key, "deg")
            style = (key_style_map or {}).get(key, {})
            explicit_unit = style.get("explicit_unit")
            out, out_unit = _format_param_for_output(val, key, native_unit, unit, explicit_unit=explicit_unit)
            if out_unit is None:
                out_key = key
            elif explicit_unit is None and style:
                out_key = key
            elif explicit_unit == out_unit and style.get("original_key"):
                out_key = str(style.get("original_key"))
            else:
                out_key = f"{key}[{out_unit}]"
            if "[" in out_key and out_key.endswith("]"):
                f.write(f'"{out_key}" = {out}\n')
            else:
                f.write(f'{out_key} = {out}\n')




def _render_value_for_template(param_name: str, value: float, explicit_unit: str | None, model_name: str) -> float:
    param_units = _resolve_param_defs_for_model(model_name)
    native_unit = param_units.get(param_name, "deg")
    out, _ = _format_param_for_output(value, param_name, native_unit, "preserve", explicit_unit=explicit_unit)
    return float(out)


def _find_param_container_in_template(doc: Any, section_names: tuple[str, ...]) -> tuple[str | None, Any]:
    for sec in section_names:
        if sec in doc and isinstance(doc[sec], dict):
            return sec, doc[sec]
    return None, doc


def write_pointing_param_toml_from_template(
    template_path: str | Path,
    output_path: str | Path,
    *,
    params_deg: dict[str, float],
    model_name: str | None = None,
    section_names: tuple[str, ...] = ("pointing_params", "fixed_params"),
) -> None:
    template_path = Path(template_path)
    output_path = Path(output_path)
    ensure_parent(output_path)
    text = template_path.read_text(encoding="utf-8")
    doc = tomlkit.parse(text)
    _, container = _find_param_container_in_template(doc, section_names)

    for raw_key in list(container.keys()):
        name, explicit_unit = _parse_param_key(str(raw_key))
        if name not in params_deg:
            continue
        container[raw_key] = _render_value_for_template(
            name,
            float(params_deg[name]),
            explicit_unit,
            str(model_name or doc.get("metadata", {}).get("model") or ""),
        )

    with output_path.open("w", encoding="utf-8") as f:
        f.write(tomlkit.dumps(doc))

def toml_scalar(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            raise ValueError("NaN/Inf cannot be written to TOML")
        return repr(v)
    return json.dumps(str(v), ensure_ascii=False)


def ensure_parent(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def dump_json(path: str | Path, obj: Any) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def _resolve_manifest_relative_path(value: Any, manifest_dir: Path) -> Any:
    if not isinstance(value, str):
        return value
    s = value.strip()
    if not s:
        return value
    p = Path(s).expanduser()
    if p.is_absolute():
        return str(p)
    return str((manifest_dir / p).resolve())


def load_manifest(manifest_path: str | Path) -> list[dict[str, Any]]:
    manifest_path = Path(manifest_path)
    manifest_dir = manifest_path.resolve().parent
    doc = read_toml(manifest_path)
    defaults = doc.get("defaults", {})
    if defaults is None:
        defaults = {}
    if not isinstance(defaults, dict):
        raise ValueError("Manifest [defaults] must be a table when provided.")
    datasets = doc.get("dataset")
    if datasets is None:
        datasets = doc.get("datasets")
    if not isinstance(datasets, list) or not datasets:
        raise ValueError("Manifest must contain one or more [[dataset]] or [[datasets]] entries.")
    path_like_keys = {"input", "input_csv", "used_param", "fix_file", "fixed_param_file", "source_model_csv"}
    merged: list[dict[str, Any]] = []
    for ds in datasets:
        if not isinstance(ds, dict):
            raise ValueError("Each manifest dataset entry must be a table.")
        row = dict(defaults)
        row.update(ds)
        if "input_csv" not in row and "input" in row:
            row["input_csv"] = row["input"]
        for key in path_like_keys:
            if key in row:
                row[key] = _resolve_manifest_relative_path(row[key], manifest_dir)
        merged.append(row)
    return merged


def read_csv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path)


def write_csv(df: pd.DataFrame, path: str | Path) -> None:
    ensure_parent(path)
    df.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def numeric_series(df: pd.DataFrame, name: str) -> np.ndarray:
    return pd.to_numeric(df[name], errors="coerce").to_numpy(dtype=float)
