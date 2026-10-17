import importlib.util, pathlib
spec = importlib.util.spec_from_file_location("nnk", pathlib.Path(__file__).parents[1] / "src" / "05_nnk.py")
nnk = importlib.util.module_from_spec(spec); spec.loader.exec_module(nnk)

def test_sizes():
    assert nnk.library_size(2) == 1024
    assert nnk.library_size(3) == 32768

def test_coverage():
    n = nnk.clones_for_coverage(1024)
    assert 3000 < n < 3200  # ~3x library size
