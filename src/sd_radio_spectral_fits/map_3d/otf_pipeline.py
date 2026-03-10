from .gridder import run_mapping_pipeline, create_grid_input
from .basketweave import solve_basket_weave_offsets, apply_basket_weave_correction
import numpy as np
import pandas as pd

def run_otf_full_pipeline(
    scantable,
    config,
    output_fits: str,
    do_basket_weave: bool = True,
    **kwargs
):
    if do_basket_weave:
        print("Executing Basket-weave correction...")
        # 1. ScantableからGridInputを生成
        grid_input = create_grid_input(scantable)
        
        # 2. オフセット計算
        offsets = solve_basket_weave_offsets(grid_input)
        
        # 3. Scantableのデータ実体へ直接オフセットを適用
        if "SCAN" in scantable.table.columns:
            scan_ids = pd.to_numeric(scantable.table["SCAN"], errors="coerce").to_numpy(dtype=np.int64)
            valid = (scan_ids >= 0) & (scan_ids < len(offsets))
            corr_vec = np.zeros(len(scan_ids))
            corr_vec[valid] = offsets[scan_ids[valid]]
            
            if isinstance(scantable.data, list):
                for i in range(len(scantable.data)):
                    scantable.data[i] -= corr_vec[i]
            else:
                scantable.data -= corr_vec[:, np.newaxis]

    print("Delegating to General Mapping Engine...")
    return run_mapping_pipeline(
        scantable=scantable,
        config=config,
        output_fits=output_fits,
        **kwargs
    )
