import numpy as np
from sdrsf_pipeline.baseline_extra import fit_baseline_poly

def test_baseline_poly_with_line_mask_and_iter():
    v = np.linspace(-100, 100, 401)
    # baseline: 1 + 0.01*v ; line: gaussian near 0
    y = 1.0 + 0.01*v + 3.0*np.exp(-(v/5.0)**2)
    base_windows = [(-100, -20), (20, 100)]
    line_windows = [(-10, 10)]
    coeff, model, info = fit_baseline_poly(v, y, base_windows, line_windows=line_windows, poly_order=1, iter_max=2, iter_sigma=3.0)
    # baseline near ends should match line-free baseline
    assert abs((y[0] - model[0]) ) < 0.1
    assert info["poly_order"] == 1
