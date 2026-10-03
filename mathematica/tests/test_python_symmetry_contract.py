import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class MathematicaSymmetryModuleContractTests(unittest.TestCase):
    def test_optional_symmetry_module_is_loaded_and_exposes_both_apis(self):
        module = ROOT / "src" / "ToeplitzSymmetryReductions.wl"
        self.assertTrue(module.is_file(), "optional symmetry module is missing")
        text = module.read_text(encoding="ascii")
        self.assertIn("SymmetricDeterminantRecurrence::usage", text)
        self.assertIn("SkewSymmetricDeterminantRecurrence::usage", text)
        self.assertIn("trSymmetricStraighteningData", text)
        self.assertIn("trKrylovScalarization", text)

        loader = (ROOT / "src" / "ToeplitzRecurrencesAll.wl").read_text(encoding="ascii")
        self.assertIn('"ToeplitzSymmetryReductions.wl"', loader)


    def test_skew_implementation_uses_two_step_transfer_krylov(self):
        text = (ROOT / "src" / "ToeplitzSymmetryReductions.wl").read_text(encoding="ascii")
        self.assertIn("trSkewTwoStepKrylov", text)
        body = text.split("SkewSymmetricDeterminantRecurrence[m_Integer?Positive", 1)[1]
        self.assertIn("TwoStepTransferMatrixExpression", body)
        self.assertNotIn("trEvenPolynomial[poly", body)


    def test_circulant_bridge_verification_is_registered(self):
        bridge = ROOT / "tests" / "ToeplitzCirculantBridgeExamples.wl"
        example = ROOT / "examples" / "CirculantBridgeExamples.wl"
        self.assertTrue(bridge.is_file())
        self.assertTrue(example.is_file())
        runner = (ROOT / "tests" / "RunToeplitzAllTests.wl").read_text(encoding="ascii")
        self.assertIn("RunToeplitzCirculantBridgeTests", runner)

    def test_symmetry_runtime_suite_is_registered(self):
        test_file = ROOT / "tests" / "ToeplitzSymmetryReductionsExamples.wl"
        self.assertTrue(test_file.is_file(), "symmetry runtime test suite is missing")
        runner = (ROOT / "tests" / "RunToeplitzAllTests.wl").read_text(encoding="ascii")
        self.assertIn("RunToeplitzSymmetryReductionTests", runner)

    def test_doset_coefficient_keys_are_raw_exponent_vectors(self):
        text = (ROOT / "src" / "ToeplitzSymmetryReductions.wl").read_text(encoding="ascii")
        helper = text.split("trHeldCoefficientAssociation[poly_, vars_List] :=", 1)[1]
        helper = helper.split("(* Universal straightening data", 1)[0]
        self.assertNotIn("HoldComplete[parts[[1]]]", helper)
        self.assertIn("parts[[1]] -> parts[[2]]", helper)


if __name__ == "__main__":
    unittest.main()
