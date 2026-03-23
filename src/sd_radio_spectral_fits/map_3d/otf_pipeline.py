from .gridder import run_mapping_pipeline, create_grid_input
from .basketweave import solve_basket_weave_offsets, apply_basket_weave_correction


_BASKET_KW_NAMES = {
    "search_radius_arcsec",
    "damp",
    "v_axis",
    "v_windows_kms",
    "channel_mask",
}

_RUNTIME_KW_NAMES = {
    "reproducible_mode",
    "workers",
    "sort_neighbors",
}

_OTF_REGION_KW_NAMES = {
    "otf_scan_region",
    "otf_scan_png",
    "otf_scan_existing_is_turn",
}


def run_otf_full_pipeline(
    scantable,
    config,
    output_fits: str,
    do_basket_weave: bool = True,
    write_diagnostics: bool = False,
    diagnostics_prefix: str | None = None,
    **kwargs
):
    basket_kwargs = {k: kwargs.pop(k) for k in list(kwargs.keys()) if k in _BASKET_KW_NAMES}
    runtime_kwargs = {k: kwargs.pop(k) for k in list(kwargs.keys()) if k in _RUNTIME_KW_NAMES}
    otf_region_kwargs = {k: kwargs.pop(k) for k in list(kwargs.keys()) if k in _OTF_REGION_KW_NAMES}

    if do_basket_weave:
        print("Executing Basket-weave correction...")

        frame = str(kwargs.get("coord_sys", "ICRS")).upper()
        projection = kwargs.get("projection", "SFL")

        grid_input = create_grid_input(
            scantable,
            frame=frame,
            projection=projection,
            **otf_region_kwargs,
        )

        offsets = solve_basket_weave_offsets(grid_input, **basket_kwargs)

        apply_basket_weave_correction(grid_input, offsets)
        scantable.data = grid_input.spec

    print("Delegating to General Mapping Engine...")
    return run_mapping_pipeline(
        scantable=scantable,
        config=config,
        output_fits=output_fits,
        write_diagnostics=write_diagnostics,
        diagnostics_prefix=diagnostics_prefix,
        **runtime_kwargs,
        **otf_region_kwargs,
        **kwargs,
    )
