from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .io_utils import parse_pointing_param_toml, write_pointing_param_toml, write_pointing_param_toml_from_template, parse_param_key_styles
from .models import canonical_model_name


@dataclass
class ApplyOptions:
    base_param: str
    delta_param: str
    output: str


def _load_update_param_source(path: str | Path) -> tuple[str | None, dict[str, float]]:
    p = Path(path)
    if p.is_dir():
        for name in ("delta_params_arcsec.toml", "delta_params.toml"):
            cand = p / name
            if cand.exists():
                return parse_pointing_param_toml(cand)
        pair_candidates = [
            (p / "delta_params_dx_arcsec.toml", p / "delta_params_dy_arcsec.toml"),
            (p / "delta_params_dx.toml", p / "delta_params_dy.toml"),
        ]
        for dx_path, dy_path in pair_candidates:
            if dx_path.exists() and dy_path.exists():
                dx_model, dx_params = parse_pointing_param_toml(dx_path)
                dy_model, dy_params = parse_pointing_param_toml(dy_path)
                dx_canon = canonical_model_name(str(dx_model or "").split("__", 1)[0])
                dy_canon = canonical_model_name(str(dy_model or "").split("__", 1)[0])
                if dx_canon != dy_canon:
                    raise ValueError(f"Model mismatch between separate delta files: {dx_model} vs {dy_model}")
                merged: dict[str, float] = {}
                for key in sorted(set(dx_params) | set(dy_params)):
                    merged[key] = float(dx_params.get(key, 0.0)) + float(dy_params.get(key, 0.0))
                return (str(dx_model or dy_model).split("__", 1)[0] if (dx_model or dy_model) else None, merged)
        if (p / "fit_result.json").exists() or (p / "absolute_params.toml").exists() or (p / "absolute_params_arcsec.toml").exists():
            raise ValueError(
                "apply accepts delta-fit outputs only. Pass a delta_params.toml file or a delta fit output directory; absolute fit outputs and fit_result.json are not valid update sources."
            )
        raise FileNotFoundError(
            f"No delta parameter TOML found in directory: {p}. apply accepts delta_params.toml (or delta_params_arcsec.toml) or a delta-fit output directory."
        )
    if p.suffix.lower() == ".json" or p.name == "fit_result.json":
        raise ValueError(
            "apply does not accept fit_result.json directly. Pass delta_params.toml, delta_params_arcsec.toml, or a delta fit output directory."
        )
    if p.name.startswith("absolute_params"):
        raise ValueError(
            "apply accepts delta parameter TOML only. absolute_params.toml is already a full model candidate and should not be added onto a base file."
        )
    return parse_pointing_param_toml(p)


def run_apply(opts: ApplyOptions) -> dict[str, float]:
    base_model, base = parse_pointing_param_toml(opts.base_param)
    delta_model, delta = _load_update_param_source(opts.delta_param)
    if canonical_model_name(str(base_model or "")) != canonical_model_name(str(delta_model or "")):
        raise ValueError(f"Model mismatch: {base_model} vs {delta_model}")
    out: dict[str, float] = {}
    for key in sorted(set(base) | set(delta)):
        out[key] = float(base.get(key, 0.0)) + float(delta.get(key, 0.0))
    key_style_map = {}
    try:
        key_style_map = parse_param_key_styles(opts.base_param)
    except Exception:
        key_style_map = {}
    output_model_label = str(base_model or delta_model or canonical_model_name(str(base_model or delta_model or "")))
    if Path(opts.base_param).exists():
        write_pointing_param_toml_from_template(
            opts.base_param,
            opts.output,
            params_deg=out,
            model_name=output_model_label,
            section_names=("pointing_params", "fixed_params"),
        )
    else:
        write_pointing_param_toml(
            opts.output,
            model_name=output_model_label,
            params_deg=out,
            unit="preserve",
            metadata_extra={"model": output_model_label},
            key_style_map=key_style_map,
        )
    return out
