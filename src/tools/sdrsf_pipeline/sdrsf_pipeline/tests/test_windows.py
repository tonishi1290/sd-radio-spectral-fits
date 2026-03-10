from sdrsf_pipeline.utils import subtract_windows

def test_subtract_windows_multiple():
    base = [(-100, -50), (50, 100)]
    line = [(-70, -60), (60, 70)]
    out = subtract_windows(base, line)
    assert out == [(-100, -70), (-60, -50), (50, 60), (70, 100)]
