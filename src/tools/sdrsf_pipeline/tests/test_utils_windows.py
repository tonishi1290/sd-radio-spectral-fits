from __future__ import annotations
from sdrsf_pipeline.utils import parse_windows, subtract_windows

def test_parse_windows():
    w = parse_windows(["1:2", "5:3"])
    assert w[0] == (1.0, 2.0)
    assert w[1] == (3.0, 5.0)

def test_subtract_windows():
    base = [(0.0, 10.0)]
    sub = [(2.0, 3.0), (7.0, 20.0)]
    out = subtract_windows(base, sub)
    assert out == [(0.0, 2.0), (3.0, 7.0)]
