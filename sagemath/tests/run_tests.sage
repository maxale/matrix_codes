"""Run the complete SageMath regression suite."""

from pathlib import Path
import runpy
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

suite = unittest.defaultTestLoader.discover(str(TESTS), pattern="test_python*.py")
result = unittest.TextTestRunner(verbosity=2).run(suite)
if not result.wasSuccessful():
    raise SystemExit(1)

for name in (
    "test_sage_row_column.sage",
    "test_sage_increasing_rows.sage",
    "test_sage_krylov.sage",
    "test_sage_fiduccia_integration.sage",
    "test_sage_symmetry_reductions.sage",
    "test_sage_circulant_bridge.sage",
):
    print(f"\n=== {name} ===")
    runpy.run_path(str(TESTS / name), run_name="__main__")

print("\n=== paper_examples.sage ===")
runpy.run_path(str(ROOT / "examples" / "paper_examples.sage"), run_name="__main__")

print("\n=== symmetry_examples.sage ===")
runpy.run_path(str(ROOT / "examples" / "symmetry_examples.sage"), run_name="__main__")

print("\n=== circulant_bridge_examples.sage ===")
runpy.run_path(str(ROOT / "examples" / "circulant_bridge_examples.sage"), run_name="__main__")

print("\nALL TOEPLITZ SAGEMATH TESTS: PASS")
