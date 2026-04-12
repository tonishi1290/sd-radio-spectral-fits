from __future__ import annotations

import builtins
import copy
import time
from pathlib import Path
from types import SimpleNamespace
import warnings

import numpy as np
import pandas as pd

from ..regrid_vlsrk import Standardizer
from ..tempscale import beameff_array, tempscal_array, ta_to_tr
from .config import MapConfig
from .core import grid_otf
from .gridder import (
    _effective_config,
    _prepare_projection_table_and_coords,
    _resolve_weight_mode,
    _assemble_grid_input,
    _resolve_otf_scan_png_path,
    _extract_restfreq,
    _print_effective_beam_summary,
)
from .otf_bundle_io import gridresult_to_otf_bundle, write_otf_bundle
from .provenance import append_bundle_provenance_step, summarize_scantable_input


def _deepcopy_meta(meta):
    try:
        return copy.deepcopy(meta)
    except Exception:
        return dict(meta or {})


def _merge_scantables(scantables) -> object:
    items = list(scantables) if isinstance(scantables, (list, tuple)) else [scantables]
    if len(items) == 0:
        raise ValueError("No scantables were supplied.")
    if len(items) == 1:
        return items[0]

    tables = []
    specs = []
    base_meta = _deepcopy_meta(getattr(items[0], "meta", {}) or {})
    nchan0 = None
    for idx, sc in enumerate(items):
        table = getattr(sc, "table", None)
        data = getattr(sc, "data", None)
        if table is None or data is None:
            raise ValueError(f"Input scantable at index {idx} must have .table and .data attributes.")
        arr = np.asarray(data)
        if arr.ndim != 2:
            raise ValueError(f"scantable.data must be 2D; got shape={arr.shape} at index {idx}")
        if nchan0 is None:
            nchan0 = int(arr.shape[1])
        elif int(arr.shape[1]) != nchan0:
            raise ValueError(
                "All scantables must currently have the same input channel count before Standardizer. "
                f"Got first nchan={nchan0}, but index {idx} has nchan={arr.shape[1]}."
            )
        tables.append(pd.DataFrame(table).copy())
        specs.append(np.asarray(arr, dtype=float))
    merged_table = pd.concat(tables, ignore_index=True)
    merged_data = np.vstack(specs)
    merged = SimpleNamespace(table=merged_table, data=merged_data, meta=base_meta)
    merged.meta.setdefault("MERGED_NINPUT", len(items))
    return merged


def grid_otf_family(
    scantables,
    config: MapConfig,
    *,
    family_label: str,
    coord_sys: str = "icrs",
    projection: str = "SFL",
    out_scale: str = "TA*",
    dv_kms: float | None = None,
    vmin_kms: float | None = None,
    vmax_kms: float | None = None,
    linefree_velocity_windows_kms: list[str] | list[tuple[float, float]] | None = None,
    ref_coord=None,
    ref_lon: float | None = None,
    ref_lat: float | None = None,
    reproducible_mode: bool | None = None,
    workers: int | None = None,
    sort_neighbors: bool | None = None,
    verbose: bool | None = None,
    otf_input_state=None,
    otf_scan_region=None,
    otf_scan_png=None,
    existing_turn_labels: str | None = None,
    otf_scan_existing_is_turn: str | None = None,
    output_fits: str | None = None,
    overwrite: bool = False,
):
    """Grid one family (e.g. X or Y) of OTF scantables into an OTFBundle."""
    _t0_prov = time.perf_counter()
    input_scantable_summary = summarize_scantable_input(scantables)
    runtime_config = _effective_config(
        config,
        reproducible_mode=reproducible_mode,
        workers=workers,
        sort_neighbors=sort_neighbors,
        verbose=verbose,
    )
    scantable = _merge_scantables(scantables)

    t = scantable.table
    if hasattr(t, "columns"):
        meta_ctype1 = str(getattr(scantable, "meta", {}).get("CTYPE1", "")).upper()
        if ("VEL" in meta_ctype1) or ("VRAD" in meta_ctype1):
            sig_cols = [c for c in t.columns if str(c).startswith("SIG_")]
            if sig_cols:
                if getattr(runtime_config, "verbose", False):
                    print(f"[info] drop SIG_* columns: {sig_cols}")
                scantable.table = t.drop(columns=sig_cols)

    std = Standardizer(scantable, v_corr_col="VFRAME")
    dv_use = dv_kms if dv_kms is not None else getattr(config, "dv_kms", None)
    vmin_use = vmin_kms if vmin_kms is not None else getattr(config, "vmin_kms", None)
    vmax_use = vmax_kms if vmax_kms is not None else getattr(config, "vmax_kms", None)
    full_matrix, v_tgt = std.get_matrix(dv=dv_use, vmin=vmin_use, vmax=vmax_use)

    if getattr(runtime_config, "verbose", False) and (vmin_use is not None or vmax_use is not None):
        print(f"   requested velocity range: {vmin_use} to {vmax_use} km/s")
    if len(v_tgt) > 1:
        dv = float(np.nanmedian(np.diff(v_tgt)))
        dv_nonlin = float(np.nanmax(np.abs(np.diff(v_tgt) - dv)))
        if getattr(runtime_config, "verbose", False):
            print("v_tgt[0:3] =", v_tgt[:3])
            print("v_tgt[-3:] =", v_tgt[-3:])
            print("dv(median) =", dv, " range =", (np.nanmin(v_tgt), np.nanmax(v_tgt)))
            print("nonlinear check (max|d-dmed|) =", dv_nonlin)
        if np.isfinite(dv_nonlin) and np.isfinite(dv) and abs(dv) > 0 and dv_nonlin > builtins.max(1e-6, abs(dv) * 1e-6):
            warnings.warn(
                f"Velocity axis is not perfectly linear (max|Δv-Δv_med|={dv_nonlin:g} km/s). "
                "The bundle header WCS will use a linear approximation.",
                RuntimeWarning,
                stacklevel=2,
            )

    table, x_arcsec, y_arcsec, lon0, lat0 = _prepare_projection_table_and_coords(
        scantable,
        frame=coord_sys,
        projection=projection,
        ref_coord=ref_coord,
        ref_lon=ref_lon,
        ref_lat=ref_lat,
    )

    out_scale_norm = str(out_scale).strip().upper()
    n_rows = len(table)
    beameff_vec = beameff_array(table, scantable.meta, n_rows)
    tempscal_vec = tempscal_array(table, scantable.meta, n_rows)

    rep_beameff = np.nan
    if out_scale_norm == "TR*":
        mask_ta = (tempscal_vec == "TA*")
        if np.any(mask_ta):
            full_matrix[mask_ta] = ta_to_tr(full_matrix[mask_ta], beameff_vec[mask_ta, None])
            if getattr(runtime_config, "verbose", False):
                print(f" -> Converted {np.count_nonzero(mask_ta)} spectra to TR*.")
        rep_beameff = 1.0
    else:
        valid_beameff = beameff_vec[np.isfinite(beameff_vec)]
        rep_beameff = float(np.median(valid_beameff)) if len(valid_beameff) > 0 else np.nan

    weight_mode = _resolve_weight_mode(runtime_config)
    grid_input, otf_scan_summary = _assemble_grid_input(
        scantable=scantable,
        table=table,
        x_arcsec=x_arcsec,
        y_arcsec=y_arcsec,
        spec=full_matrix,
        out_scale_norm=out_scale_norm,
        beameff_vec=beameff_vec,
        weight_mode=weight_mode,
        otf_input_state=otf_input_state,
        otf_scan_region=otf_scan_region,
        otf_scan_png=_resolve_otf_scan_png_path(
            otf_scan_png,
            default_path=None if output_fits is None else str(Path(output_fits).with_suffix(Path(output_fits).suffix + ".scan_region.png")),
        ),
        existing_turn_labels=existing_turn_labels,
        otf_scan_existing_is_turn=otf_scan_existing_is_turn,
    )
    setattr(grid_input, "_otf_scan_summary", otf_scan_summary)

    if getattr(runtime_config, "verbose", False) and otf_scan_summary is not None:
        state_txt = otf_scan_summary.get("input_state")
        if state_txt:
            print(f"   otf_input_state={state_txt}")

    res = grid_otf(grid_input, runtime_config)
    if getattr(runtime_config, "verbose", False):
        _print_effective_beam_summary(res.meta or {})
    if not res.meta:
        res.meta = {}
    res.meta["RESTFREQ"] = _extract_restfreq(scantable, table)
    res.meta["SPECSYS"] = "LSRK"
    res.meta["family_label"] = str(family_label)
    res.meta["otf_scan_summary"] = otf_scan_summary
    if isinstance(scantables, (list, tuple)):
        res.meta["merged_inputs"] = int(len(scantables))

    bundle = gridresult_to_otf_bundle(
        res,
        v_tgt=v_tgt,
        coord_sys=coord_sys,
        projection=projection,
        lon0=lon0,
        lat0=lat0,
        config=runtime_config,
        out_scale=out_scale_norm,
        rep_beameff=rep_beameff,
        family_label=family_label,
        extra_meta={
            "coord_sys": str(coord_sys),
            "projection": str(projection),
            "ref_lon": float(lon0),
            "ref_lat": float(lat0),
            "out_scale": str(out_scale_norm),
            "vmin_kms": None if vmin_use is None else float(vmin_use),
            "vmax_kms": None if vmax_use is None else float(vmax_use),
            "dv_kms": None if dv_use is None else float(dv_use),
            "linefree_velocity_windows_kms": list(linefree_velocity_windows_kms or []) if linefree_velocity_windows_kms is not None else None,
            "baseline_subtracted": False,
        },
    )
    append_bundle_provenance_step(
        bundle,
        input_bundles=[],
        op_id="otf.grid.family.v1",
        module=__name__,
        function="grid_otf_family",
        kind="main",
        params_input={
            "family_label": str(family_label),
            "coord_sys": str(coord_sys),
            "projection": str(projection),
            "out_scale": str(out_scale),
            "dv_kms": dv_kms,
            "vmin_kms": vmin_kms,
            "vmax_kms": vmax_kms,
            "linefree_velocity_windows_kms": list(linefree_velocity_windows_kms or []) if linefree_velocity_windows_kms is not None else None,
            "ref_coord": ref_coord,
            "ref_lon": ref_lon,
            "ref_lat": ref_lat,
            "reproducible_mode": reproducible_mode,
            "workers": workers,
            "sort_neighbors": sort_neighbors,
            "verbose": verbose,
            "otf_input_state": otf_input_state,
            "otf_scan_region": otf_scan_region,
            "existing_turn_labels": existing_turn_labels,
            "otf_scan_existing_is_turn": otf_scan_existing_is_turn,
        },
        params_config=runtime_config,
        params_resolved={
            "family_label": str(family_label),
            "coord_sys": str(coord_sys),
            "projection": str(projection),
            "ref_lon_used_deg": float(lon0),
            "ref_lat_used_deg": float(lat0),
            "out_scale_used": str(out_scale_norm),
            "dv_kms_used": None if dv_use is None else float(dv_use),
            "vmin_kms_used": None if vmin_use is None else float(vmin_use),
            "vmax_kms_used": None if vmax_use is None else float(vmax_use),
            "weight_mode_used": str(weight_mode),
            "rep_beameff_used": float(rep_beameff) if np.isfinite(rep_beameff) else None,
            "otf_scan_summary": otf_scan_summary,
            "input_scantable_summary": input_scantable_summary,
        },
        results_summary={
            "cube_shape": [int(v) for v in bundle.data.shape],
            "support_npix": int(np.count_nonzero(bundle.support_mask)) if bundle.support_mask is not None else None,
            "variance_source": bundle.meta.get("variance_source"),
            "restfreq_hz": bundle.meta.get("RESTFREQ"),
            "specsys": bundle.meta.get("SPECSYS"),
            "meta_keys": sorted(str(k) for k in (res.meta or {}).keys()),
        },
        duration_sec=float(time.perf_counter() - _t0_prov),
    )
    from .mosaic import attach_mosaic_products

    attach_mosaic_products(
        bundle,
        linefree_velocity_windows_kms=linefree_velocity_windows_kms,
        overwrite=True,
        in_place=True,
        _record_provenance=True,
    )
    if output_fits is not None:
        write_otf_bundle(bundle, output_fits, overwrite=overwrite)
    return bundle
