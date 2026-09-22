from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_CODE = [
    ROOT / "src" / "toeplitz_recurrences.py",
    ROOT / "src" / "toeplitz_symmetry_reductions.py",
    ROOT / "src" / "toeplitz_recurrences.sage",
]
SAGE_SURFACES = list((ROOT / "tests").glob("*.sage")) + list((ROOT / "examples").glob("*.sage"))
DOCS = [ROOT / "README.md"]


class NoSymbolicRingTests(unittest.TestCase):
    def test_no_sage_symbolic_ring_dependency(self):
        offenders = []
        for path in PUBLIC_CODE + SAGE_SURFACES + DOCS:
            text = path.read_text(encoding="utf-8")
            for needle in (
                "sage.SR",
                "from sage.all import SR",
                "SR(",
                " var(",
                "= var(",
                'var("',
                "var('",
            ):
                if needle in text:
                    offenders.append((path.relative_to(ROOT).as_posix(), needle))
        self.assertFalse(
            offenders,
            f"Sage symbolic-ring dependencies remain: {offenders}",
        )


if __name__ == "__main__":
    unittest.main()
