from __future__ import annotations

import numpy as np
from sdrsf_pipeline.baseline import fit_polynomial_baseline, rms_after_polyfit

def test_baseline_fit_linear():
    x = np.linspace(-50, 50, 101)
    y = 0.1 * x + 2.0
    # add a line feature outside baseline windows
    y2 = y.copy()
    y2[(x>-5)&(x<5)] += 10.0

    windows = [(-50,-10), (10,50)]
    base, coeffs, mask = fit_polynomial_baseline(x, y2, windows=windows, order=1)
    # baseline should match original linear trend
    assert np.allclose(base[mask], y[mask], atol=1e-6)

def test_rms_after_polyfit_zero_noise():
    x = np.linspace(-50, 50, 101)
    y = 0.2 * x - 1.0
    windows = [(-50,50)]
    r = rms_after_polyfit(x, y, windows=windows, order=1)
    assert abs(r) < 1e-10
