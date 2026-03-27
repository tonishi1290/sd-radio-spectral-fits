from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from .gridder import run_mapping_pipeline, create_grid_input, _resolve_projection_reference_for_scantables
from .basketweave import (
    basket_weave_inplace,
    _merge_grid_inputs,
)
from .otf_scan_region import format_otf_scan_summary


_BASKET_KW_NAMES = {
    "search_radius_arcsec",
    "damp",
    "v_axis",
    "v_windows_kms",
    "channel_mask",
    "cross_direction_only",
    "orthogonality_tolerance_deg",
    "fallback_to_all_cross_scan_pairs",
    "pair_mode",
    "offset_model",
    "reference_mode",
    "reference_direction",
    "reference_scan_id",
    "reference_constraint_weight",
    "reference_constrain_terms",
    "linefree_velocity_windows_kms",
}

_RUNTIME_KW_NAMES = {
    "reproducible_mode",
    "workers",
    "sort_neighbors",
    "verbose",
}

_OTF_REGION_KW_NAMES = {
    "otf_input_state",
    "otf_scan_region",
    "otf_scan_png",
    "existing_turn_labels",
    "otf_scan_existing_is_turn",
}


def _validate_multi_scantable_sequence(scantables, *, caller: str):
    if not isinstance(scantables, (list, tuple)) or len(scantables) == 0:
        raise ValueError(f'{caller} requires a non-empty list/tuple of scantables.')
    if all(isinstance(item, (str, Path)) for item in scantables):
        raise ValueError(f'{caller} does not accept a file-path list. Pass loaded scantable objects instead.')


def _count_scantable_rows(scantable) -> int:
    table = getattr(scantable, 'table', None)
    if table is None:
        raise ValueError('Each scantable must have a .table attribute.')
    return int(len(table))


def _ensure_not_all_empty_scantables(scantables, *, caller: str):
    nrows = [_count_scantable_rows(sc) for sc in scantables]
    if len(nrows) == 0 or not any(n > 0 for n in nrows):
        raise ValueError(f'{caller} received only empty scantables; at least one non-empty scantable is required.')
    return nrows

def _filter_nonempty_scantables(scantables, *, caller: str):
    nrows = _ensure_not_all_empty_scantables(scantables, caller=caller)
    kept = []
    skipped = []
    for idx, (sc, nrow) in enumerate(zip(scantables, nrows)):
        if int(nrow) > 0:
            kept.append((idx, sc))
        else:
            skipped.append(int(idx))
    if skipped:
        warnings.warn(
            f"{caller} is skipping empty scantable inputs at indices {skipped}.",
            RuntimeWarning,
            stacklevel=2,
        )
    return kept, skipped



def _extract_named_kwargs(kwargs: dict, names: set[str]) -> dict:
    return {k: kwargs.pop(k) for k in list(kwargs.keys()) if k in names}


def _normalize_single_otf_scan_png_path(path_like, *, default_path: str | None = None):
    if path_like in (None, False):
        return None
    if path_like is True:
        return str(default_path if default_path is not None else 'otf_scan_region.png')
    return str(path_like)


def _indexed_otf_scan_png_path(path_like, index: int, total: int, *, default_path: str | None = None):
    resolved = _normalize_single_otf_scan_png_path(path_like, default_path=default_path)
    if resolved is None:
        return None
    path = Path(str(resolved))
    if total <= 1:
        return str(path)
    return str(path.with_name(f'{path.stem}.{index:03d}{path.suffix}'))


def _resolve_otf_scan_png_sequence(otf_scan_png, count: int, *, default_path: str | None = None):
    if count <= 0:
        return []
    if isinstance(otf_scan_png, (list, tuple)):
        if len(otf_scan_png) != count:
            raise ValueError('otf_scan_png sequence length must match the number of scantables.')
        out = []
        for i, item in enumerate(otf_scan_png):
            if item is True:
                out.append(_indexed_otf_scan_png_path(item, i, count, default_path=default_path))
            else:
                out.append(_normalize_single_otf_scan_png_path(item, default_path=default_path))
        return out
    return [
        _indexed_otf_scan_png_path(otf_scan_png, i, count, default_path=default_path)
        for i in range(count)
    ]


def _concat_scantable_data(scantables) -> np.ndarray:
    arrays = [np.asarray(sc.data) for sc in scantables]
    if not arrays:
        raise ValueError('At least one scantable is required.')
    first = arrays[0]
    if first.ndim not in (1, 2):
        raise ValueError(f'scantable.data must be 1D or 2D, got ndim={first.ndim}')
    if first.ndim == 2:
        nchan = int(first.shape[1])
        for idx, arr in enumerate(arrays[1:], start=1):
            if arr.ndim != 2 or int(arr.shape[1]) != nchan:
                raise ValueError(
                    f'All input scantables must have the same spectral channel count for multi-OTF gridding. '
                    f'Input 0 has {nchan} channels, but input {idx} has shape {arr.shape}.'
                )
    else:
        for idx, arr in enumerate(arrays[1:], start=1):
            if arr.ndim != 1:
                raise ValueError('Cannot mix 1D and 2D scantable.data in run_otf_full_pipeline_multi().')
    return np.concatenate(arrays, axis=0)


def _materialize_single_otf_table(
    scantable,
    *,
    frame: str,
    projection: str,
    ref_coord=None,
    ref_lon: float | None = None,
    ref_lat: float | None = None,
    otf_region_kwargs: dict,
    table_index: int,
):
    grid_input = create_grid_input(
        scantable,
        ref_coord=ref_coord,
        ref_lon=ref_lon,
        ref_lat=ref_lat,
        frame=frame,
        projection=projection,
        **otf_region_kwargs,
    )
    table = pd.DataFrame(scantable.table).copy()
    scan_id = None if grid_input.scan_id is None else np.asarray(grid_input.scan_id, dtype=np.int64)
    is_turn = None if grid_input.is_turnaround is None else np.asarray(grid_input.is_turnaround, dtype=bool)
    scan_dir = None if grid_input.scan_dir is None else np.asarray(grid_input.scan_dir, dtype=float)
    if scan_id is not None:
        table['OTF_EFFSCAN'] = scan_id
        table['SCAN'] = scan_id
    if is_turn is not None:
        table['OTF_IS_TURN'] = is_turn
        table['IS_TURN'] = is_turn
    if scan_dir is not None:
        table['OTF_SCAN_DIR'] = scan_dir
    table['OTF_TABLE_INDEX'] = int(table_index)
    summary = getattr(grid_input, '_otf_scan_summary', None)
    if isinstance(summary, dict) and int(summary.get('num_runs', 0)) == 0:
        warnings.warn(
            f'OTF scan-region found 0 accepted runs for scantable index {int(table_index)}.',
            RuntimeWarning,
            stacklevel=2,
        )
    return grid_input, table, summary


def _build_multi_scantable_for_mapping(
    scantables,
    *,
    frame: str,
    projection: str,
    ref_coord=None,
    ref_lon: float | None = None,
    ref_lat: float | None = None,
    otf_region_kwargs: dict,
    kept_inputs: list[tuple[int, object]] | None = None,
    skipped_empty_indices: list[int] | None = None,
):
    _validate_multi_scantable_sequence(scantables, caller='run_otf_full_pipeline_multi()')
    if kept_inputs is None or skipped_empty_indices is None:
        kept_inputs, skipped_empty_indices = _filter_nonempty_scantables(
            scantables,
            caller='run_otf_full_pipeline_multi()',
        )
    else:
        kept_inputs = [(int(idx), sc) for idx, sc in kept_inputs]
        skipped_empty_indices = [int(v) for v in skipped_empty_indices]
    effective_scantables = [sc for _, sc in kept_inputs]
    if ref_coord is None and (ref_lon is None or ref_lat is None):
        ref_lon, ref_lat = _resolve_projection_reference_for_scantables(
            effective_scantables,
            frame=frame,
            ref_coord=ref_coord,
            ref_lon=ref_lon,
            ref_lat=ref_lat,
        )
        ref_coord = (float(ref_lon), float(ref_lat))
    png_seq_all = _resolve_otf_scan_png_sequence(
        otf_region_kwargs.get('otf_scan_png', None),
        len(scantables),
    )
    grid_inputs = []
    tables = []
    summaries = []
    used_input_indices = []
    for orig_idx, sc in kept_inputs:
        item_otf_kwargs = dict(otf_region_kwargs)
        item_otf_kwargs['otf_scan_png'] = png_seq_all[orig_idx]
        gi, tbl, summary = _materialize_single_otf_table(
            sc,
            frame=frame,
            projection=projection,
            ref_coord=ref_coord,
            ref_lon=ref_lon,
            ref_lat=ref_lat,
            otf_region_kwargs=item_otf_kwargs,
            table_index=orig_idx,
        )
        setattr(gi, '_source_input_index', int(orig_idx))
        grid_inputs.append(gi)
        tables.append(tbl)
        summaries.append(summary)
        used_input_indices.append(int(orig_idx))

    merged_gi = _merge_grid_inputs(grid_inputs)
    start = 0
    for idx, (tbl, gi) in enumerate(zip(tables, grid_inputs)):
        n = int(len(np.asarray(gi.x)))
        stop = start + n
        sid = np.asarray(merged_gi.scan_id[start:stop], dtype=np.int64)
        tbl['OTF_EFFSCAN_GLOBAL'] = sid
        tbl['SCAN'] = sid
        if merged_gi.is_turnaround is not None:
            turn = np.asarray(merged_gi.is_turnaround[start:stop], dtype=bool)
            tbl['OTF_IS_TURN_GLOBAL'] = turn
            tbl['IS_TURN'] = turn
        if merged_gi.scan_dir is not None:
            tbl['OTF_SCAN_DIR_GLOBAL'] = np.asarray(merged_gi.scan_dir[start:stop], dtype=float)
        start = stop

    merged_table = pd.concat(tables, ignore_index=True, sort=False)
    merged_data = _concat_scantable_data(effective_scantables)
    merged_meta = dict(getattr(effective_scantables[0], 'meta', {}) or {})
    merged_scantable = SimpleNamespace(table=merged_table, data=merged_data, meta=merged_meta)

    source_summary = getattr(merged_gi, '_otf_scan_summary', None)
    if len(scantables) > 1 or skipped_empty_indices:
        if isinstance(source_summary, dict) and source_summary.get('kind') == 'multi_input':
            merged_summary = dict(source_summary)
        else:
            merged_summary = {
                'summary_version': 2,
                'kind': 'multi_input',
                'reference_scan_id_semantics': 'global after merge',
                'num_inputs': int(len(grid_inputs)),
                'inputs': summaries,
            }
            if isinstance(source_summary, dict):
                merged_summary['source_summary'] = source_summary
        input_states = sorted({str(s.get('input_state')) for s in summaries if isinstance(s, dict) and s.get('input_state') is not None})
        if len(input_states) == 1:
            merged_summary['input_state'] = input_states[0]
    else:
        merged_summary = dict(source_summary) if isinstance(source_summary, dict) else {
            'summary_version': 2,
            'kind': 'single_table',
            'inputs': summaries,
        }
    merged_summary['num_inputs_original'] = int(len(scantables))
    merged_summary['num_inputs_used'] = int(len(grid_inputs))
    merged_summary['used_input_indices'] = [int(v) for v in used_input_indices]
    merged_summary['skipped_empty_input_indices'] = [int(v) for v in skipped_empty_indices]
    setattr(merged_scantable, '_otf_scan_summary', merged_summary)
    return merged_scantable, merged_summary


def run_otf_full_pipeline(
    scantable,
    config,
    output_fits: str,
    do_basket_weave: bool = True,
    write_diagnostics: bool = False,
    diagnostics_prefix: str | None = None,
    **kwargs
):
    basket_kwargs = _extract_named_kwargs(kwargs, _BASKET_KW_NAMES)
    runtime_kwargs = _extract_named_kwargs(kwargs, _RUNTIME_KW_NAMES)
    otf_region_kwargs = _extract_named_kwargs(kwargs, _OTF_REGION_KW_NAMES)
    ref_coord = kwargs.pop('ref_coord', None)
    ref_lon = kwargs.pop('ref_lon', None)
    ref_lat = kwargs.pop('ref_lat', None)

    frame = str(kwargs.get("coord_sys", "ICRS")).upper()
    projection = kwargs.get("projection", "SFL")
    verbose = bool(getattr(config, 'verbose', False) or runtime_kwargs.get('verbose', False))
    common_ref_coord = ref_coord
    common_ref_lon = ref_lon
    common_ref_lat = ref_lat
    if common_ref_coord is None:
        if common_ref_lon is not None and common_ref_lat is not None:
            common_ref_coord = (float(common_ref_lon), float(common_ref_lat))
        else:
            common_ref_lon, common_ref_lat = _resolve_projection_reference_for_scantables(
                [scantable],
                frame=frame,
                ref_coord=ref_coord,
                ref_lon=ref_lon,
                ref_lat=ref_lat,
            )
            common_ref_coord = (float(common_ref_lon), float(common_ref_lat))

    if do_basket_weave:
        print("Executing Basket-weave correction...")
        bw_result = basket_weave_inplace(
            scantable,
            projection=projection,
            ref_coord=common_ref_coord,
            frame=frame,
            **basket_kwargs,
            **otf_region_kwargs,
        )
        if verbose:
            print(
                '   basketweave: '
                f'scans={bw_result.num_scans} '
                f'pairs={bw_result.num_pairs_used}/{bw_result.num_pairs_total} '
                f'used_channels={bw_result.used_channel_count if bw_result.used_channel_count is not None else "?"}'
            )

    print("Delegating to General Mapping Engine...")
    otf_region_kwargs_for_mapping = dict(otf_region_kwargs)
    if do_basket_weave:
        otf_region_kwargs_for_mapping['otf_scan_png'] = None
    return run_mapping_pipeline(
        scantable=scantable,
        config=config,
        output_fits=output_fits,
        write_diagnostics=write_diagnostics,
        diagnostics_prefix=diagnostics_prefix,
        ref_coord=common_ref_coord,
        ref_lon=common_ref_lon,
        ref_lat=common_ref_lat,
        **runtime_kwargs,
        **otf_region_kwargs_for_mapping,
        **kwargs,
    )


def run_otf_full_pipeline_multi(
    scantables,
    config,
    output_fits: str,
    do_basket_weave: bool = True,
    write_diagnostics: bool = False,
    diagnostics_prefix: str | None = None,
    **kwargs,
):
    _validate_multi_scantable_sequence(scantables, caller='run_otf_full_pipeline_multi()')
    kept_inputs, skipped_empty_indices = _filter_nonempty_scantables(
        scantables,
        caller='run_otf_full_pipeline_multi()',
    )
    effective_scantables = [sc for _, sc in kept_inputs]

    basket_kwargs = _extract_named_kwargs(kwargs, _BASKET_KW_NAMES)
    runtime_kwargs = _extract_named_kwargs(kwargs, _RUNTIME_KW_NAMES)
    otf_region_kwargs = _extract_named_kwargs(kwargs, _OTF_REGION_KW_NAMES)
    ref_coord = kwargs.pop('ref_coord', None)
    ref_lon = kwargs.pop('ref_lon', None)
    ref_lat = kwargs.pop('ref_lat', None)

    frame = str(kwargs.get('coord_sys', 'ICRS')).upper()
    projection = kwargs.get('projection', 'SFL')
    verbose = bool(getattr(config, 'verbose', False) or runtime_kwargs.get('verbose', False))
    default_png_base = str(Path(output_fits).with_suffix(Path(output_fits).suffix + '.scan_region.png'))
    common_ref_lon, common_ref_lat = _resolve_projection_reference_for_scantables(
        effective_scantables, frame=frame, ref_coord=ref_coord, ref_lon=ref_lon, ref_lat=ref_lat
    )
    common_ref_coord = (float(common_ref_lon), float(common_ref_lat))
    otf_region_kwargs_for_multi = dict(otf_region_kwargs)
    png_seq = _resolve_otf_scan_png_sequence(
        otf_region_kwargs.get('otf_scan_png', None),
        len(scantables),
        default_path=default_png_base,
    )

    if skipped_empty_indices and verbose:
        print(f'   skipping empty scantables at indices {skipped_empty_indices}')

    if do_basket_weave:
        print('1. Basket-weave correction across multiple scantables...')
        bw_region_kwargs = dict(otf_region_kwargs)
        bw_region_kwargs['otf_scan_png'] = png_seq
        bw_result = basket_weave_inplace(
            effective_scantables,
            projection=projection,
            ref_coord=common_ref_coord,
            frame=frame,
            **basket_kwargs,
            **bw_region_kwargs,
        )
        if verbose:
            print(
                '   basketweave: '
                f'scans={bw_result.num_scans} '
                f'pairs={bw_result.num_pairs_used}/{bw_result.num_pairs_total} '
                f'used_channels={bw_result.used_channel_count if bw_result.used_channel_count is not None else "?"}'
            )

    print('2. Materializing safe OTF columns for merged gridding...')
    if do_basket_weave:
        otf_region_kwargs_for_multi['otf_scan_png'] = None
    else:
        otf_region_kwargs_for_multi['otf_scan_png'] = png_seq
    merged_scantable, merged_summary = _build_multi_scantable_for_mapping(
        scantables,
        frame=frame,
        projection=projection,
        ref_coord=common_ref_coord,
        ref_lon=common_ref_lon,
        ref_lat=common_ref_lat,
        otf_region_kwargs=otf_region_kwargs_for_multi,
        kept_inputs=kept_inputs,
        skipped_empty_indices=skipped_empty_indices,
    )
    if verbose and merged_summary is not None:
        print(f'   {format_otf_scan_summary(merged_summary)}')

    print('3. Delegating merged data to General Mapping Engine...')
    return run_mapping_pipeline(
        scantable=merged_scantable,
        config=config,
        output_fits=output_fits,
        write_diagnostics=write_diagnostics,
        diagnostics_prefix=diagnostics_prefix,
        ref_coord=common_ref_coord,
        ref_lon=common_ref_lon,
        ref_lat=common_ref_lat,
        otf_input_state='use_existing_labels',
        otf_scan_region=None,
        otf_scan_png=None,
        existing_turn_labels=None,
        otf_scan_existing_is_turn=None,
        **runtime_kwargs,
        **kwargs,
    )
