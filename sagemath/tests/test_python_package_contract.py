import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
import toeplitz_recurrences as tr


class PackageContractTests(unittest.TestCase):
    def test_release_metadata(self):
        self.assertEqual(tr.SOFTWARE_VERSION, "5")
        self.assertEqual(
            tr.PAPER_TITLE,
            "Constructive recurrences for determinants and permanents of banded Toeplitz matrices",
        )

    def test_internal_remote_polynomial_variable_names_are_fresh(self):
        first = tr._fresh_polynomial_variable_name()
        second = tr._fresh_polynomial_variable_name()
        self.assertNotEqual(first, second)
        self.assertTrue(first.startswith("toeplitz_lambda_"))
        self.assertRegex(first, r"^[A-Za-z][A-Za-z0-9_]*$")

    def test_public_api_is_complete(self):
        expected = {
            "m_reduce", "laplace_d", "laplace_p",
            "lin_rec_for_dtm", "lin_rec_for_ptm",
            "lin_rec_for_dtm_increasing_rows", "lin_rec_for_ptm_increasing_rows",
            "fiduccia_pol_squarings", "term_for_dtm", "term_for_ptm",
            "sage_cfinite_term",
        }
        self.assertTrue(expected.issubset(set(tr.__all__)))

    def test_fraction_simplifier_does_not_replace_value_by_factorization(self):
        class FactorableOnly:
            def factor(self):
                return ("factorization-object",)
        value = FactorableOnly()
        self.assertIs(tr._simplify_fraction(value), value)

    def test_distribution_scaffolding_exists(self):
        for name in [
            "README.md", "CITATION.cff", "CHANGELOG.md", "PACKAGE_INFO.txt",
            "VERIFICATION.md", "src/toeplitz_recurrences.py", "src/toeplitz_symmetry_reductions.py",
            "src/toeplitz_recurrences.sage", "tests/run_tests.sage",
            "examples/compare_fiduccia.sage", "examples/paper_examples.sage",
            "examples/symmetry_examples.sage", "tests/test_sage_symmetry_reductions.sage",
            "tools/verify_release.py",
        ]:
            with self.subTest(name=name):
                self.assertTrue((ROOT / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
