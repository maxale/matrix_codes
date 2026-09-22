import ast
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "toeplitz_recurrences.py"
README = ROOT / "README.md"
PAPER_TITLE = "Constructive recurrences for determinants and permanents of banded Toeplitz matrices"

EXPECTED_DEFAULTS = {
    "fiduccia_pol_squarings": {"modulus": "None", "term_count": "1"},
    "laplace_d": {"diagonal_values": "None", "position_dependent": "False", "band_entry_function": "None"},
    "laplace_p": {"diagonal_values": "None", "position_dependent": "False", "band_entry_function": "None"},
    "lin_rec_for_dtm": {
        "diagonal_values": "None", "scalar_recurrence": "'auto'",
        "characteristic_polynomial_variable": "'x'", "max_characteristic_order": "20",
        "verbose": "True", "position_dependent": "False", "band_entry_function": "None",
        "scalar_recurrence_method": "'characteristic_polynomial'",
    },
    "lin_rec_for_ptm": {
        "diagonal_values": "None", "scalar_recurrence": "'auto'",
        "characteristic_polynomial_variable": "'x'", "max_characteristic_order": "20",
        "verbose": "True", "position_dependent": "False", "band_entry_function": "None",
        "scalar_recurrence_method": "'characteristic_polynomial'",
    },
    "lin_rec_for_dtm_increasing_rows": {
        "diagonal_values": "None", "ell_range": "'forward'",
        "characteristic_polynomial_variable": "'x'", "verbose": "True",
    },
    "lin_rec_for_ptm_increasing_rows": {
        "diagonal_values": "None", "ell_range": "'forward'",
        "characteristic_polynomial_variable": "'x'", "verbose": "True",
    },
    "term_for_dtm": {
        "diagonal_values": "None", "scalar_recurrence_method": "'characteristic_polynomial'",
        "modulus": "None", "term_count": "1", "verbose": "False", "position_dependent": "False",
    },
    "term_for_ptm": {
        "diagonal_values": "None", "scalar_recurrence_method": "'characteristic_polynomial'",
        "modulus": "None", "term_count": "1", "verbose": "False", "position_dependent": "False",
    },
    "sage_cfinite_term": {"term_count": "1"},
}


def source_defaults():
    mod = ast.parse(SOURCE.read_text(encoding="utf-8"))
    out = {}
    for node in mod.body:
        if isinstance(node, ast.FunctionDef) and node.name in EXPECTED_DEFAULTS:
            out[node.name] = {
                arg.arg: ast.unparse(default)
                for arg, default in zip(node.args.kwonlyargs, node.args.kw_defaults)
                if default is not None
            }
    return out


class ReleaseContractTests(unittest.TestCase):
    def test_required_layout(self):
        required = [
            "README.md", "CITATION.cff", "CHANGELOG.md", "PACKAGE_INFO.txt", "VERIFICATION.md",
            ".gitignore", "src/toeplitz_recurrences.py", "src/toeplitz_symmetry_reductions.py", "src/toeplitz_recurrences.sage",
            "examples/paper_examples.sage", "examples/symmetry_examples.sage", "examples/circulant_bridge_examples.sage", "examples/compare_fiduccia.sage",
            "tests/run_tests.sage", "tests/test_sage_symmetry_reductions.sage", "tests/test_sage_circulant_bridge.sage",
            "tests/test_python_symmetry_logic.py", "tools/verify_release.py",
        ]
        self.assertEqual([p for p in required if not (ROOT / p).is_file()], [])

    def test_source_defaults_are_expected(self):
        got = source_defaults()
        for function, defaults in EXPECTED_DEFAULTS.items():
            self.assertIn(function, got)
            for key, value in defaults.items():
                self.assertEqual(got[function].get(key), value, (function, key))

    def test_readme_documents_public_defaults(self):
        text = README.read_text(encoding="utf-8")
        tokens = [
            "## Public API and defaults",
            "`diagonal_values` | `None`", "`scalar_recurrence` | `\"auto\"`",
            "`characteristic_polynomial_variable` | `\"x\"`", "`max_characteristic_order` | `20`",
            "`verbose` | `True`", "`position_dependent` | `False`", "`band_entry_function` | `None`",
            "`scalar_recurrence_method` | `\"characteristic_polynomial\"`",
            "`ell_range` | `\"forward\"`", "`modulus` | `None`", "`term_count` | `1`",
            "`verbose` | `False`", "`sage_cfinite_term",
            "symmetric_determinant_recurrence", "skew_symmetric_determinant_recurrence",
            "two_step_observable_order", "full_annihilator_order", "Toeplitz-circulant bridge",
        ]
        for token in tokens:
            self.assertIn(token, text)

    def test_no_manuscript_revision_leakage(self):
        pattern = re.compile(
            r"toeplitz_recurrences_rev\d+\.tex|paper\s+rev|manuscript\s+rev|" + "rev" + r"(?:36|009)",
            re.I,
        )
        for path in ROOT.rglob("*"):
            if path.is_file() and path.suffix.lower() in {".py", ".sage", ".md", ".txt", ".cff"}:
                self.assertIsNone(pattern.search(path.read_text(encoding="utf-8")), str(path.relative_to(ROOT)))

    def test_ascii_public_text(self):
        for path in ROOT.rglob("*"):
            if path.is_file() and path.suffix.lower() in {".py", ".sage", ".md", ".txt", ".cff"}:
                self.assertTrue(all(b < 128 for b in path.read_bytes()), str(path.relative_to(ROOT)))

    def test_paper_and_fiduccia_provenance(self):
        source = SOURCE.read_text(encoding="utf-8")
        readme = README.read_text(encoding="utf-8")
        self.assertIn(PAPER_TITLE, source)
        self.assertIn(PAPER_TITLE, readme)
        self.assertIn("Efficient Computation of Terms of Linear Recurrence Sequences of Any Order", readme)
        self.assertIn("INTEGERS 18", readme)
        self.assertIn("10.1137/0214007", readme)
        self.assertIn("hardcoded modular-polynomial-squaring formulas", readme)


if __name__ == "__main__":
    unittest.main(verbosity=2)
