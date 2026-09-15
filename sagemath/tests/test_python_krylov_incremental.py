import sys
import unittest
import random
from fractions import Fraction

import sympy as sp
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from toeplitz_recurrences import _IncrementalRowBasis


class IncrementalRowBasisTests(unittest.TestCase):
    def make_basis(self, width):
        return _IncrementalRowBasis(width, Fraction(0), Fraction(1))

    def test_reports_first_dependence_in_original_row_coordinates(self):
        basis = self.make_basis(3)
        self.assertIsNone(basis.add([Fraction(1), Fraction(0), Fraction(0)]))
        self.assertIsNone(basis.add([Fraction(0), Fraction(1), Fraction(0)]))
        relation = basis.add([Fraction(2), Fraction(3), Fraction(0)])
        self.assertEqual(relation, [Fraction(-2), Fraction(-3), Fraction(1)])

    def test_handles_new_pivot_inserted_before_an_existing_pivot(self):
        basis = self.make_basis(3)
        self.assertIsNone(basis.add([Fraction(0), Fraction(1), Fraction(0)]))
        self.assertIsNone(basis.add([Fraction(1), Fraction(1), Fraction(0)]))
        relation = basis.add([Fraction(2), Fraction(3), Fraction(0)])
        self.assertEqual(relation, [Fraction(-1), Fraction(-2), Fraction(1)])

    def test_tracks_full_rank_then_detects_relation(self):
        basis = self.make_basis(3)
        rows = [
            [Fraction(1), Fraction(2), Fraction(0)],
            [Fraction(0), Fraction(1), Fraction(1)],
            [Fraction(1), Fraction(0), Fraction(1)],
        ]
        for row in rows:
            self.assertIsNone(basis.add(row))
        self.assertEqual(basis.rank, 3)
        relation = basis.add([
            rows[0][0] + 2 * rows[1][0] - rows[2][0],
            rows[0][1] + 2 * rows[1][1] - rows[2][1],
            rows[0][2] + 2 * rows[1][2] - rows[2][2],
        ])
        self.assertEqual(
            relation,
            [Fraction(-1), Fraction(-2), Fraction(1), Fraction(1)],
        )

    def test_random_exact_dependence_certificates(self):
        rng = random.Random(20260912)
        for width in range(2, 8):
            basis = self.make_basis(width)
            rows = []
            for i in range(width):
                row = [Fraction(0)] * i + [Fraction(1)]
                row += [Fraction(rng.randint(-4, 4)) for _ in range(width-i-1)]
                rows.append(row)
                self.assertIsNone(basis.add(row))
            coeffs = [Fraction(rng.randint(-3, 3)) for _ in range(width)]
            dependent = [
                sum(coeffs[i] * rows[i][j] for i in range(width))
                for j in range(width)
            ]
            relation = basis.add(dependent)
            self.assertEqual(relation, [-c for c in coeffs] + [Fraction(1)])

    def test_krylov_dependence_matches_sympy_nullspace(self):
        rng = random.Random(20260913)
        for width in range(2, 7):
            for _ in range(8):
                matrix = [
                    [Fraction(rng.randint(-3, 3)) for _ in range(width)]
                    for _ in range(width)
                ]
                row = [Fraction(1)] + [Fraction(0)] * (width - 1)
                rows = [row]
                basis = self.make_basis(width)
                self.assertIsNone(basis.add(row))
                for _step in range(width):
                    nxt = [
                        sum(rows[-1][i] * matrix[i][j] for i in range(width))
                        for j in range(width)
                    ]
                    trial_rows = rows + [nxt]
                    relation = basis.add(nxt)
                    trial = sp.Matrix(
                        [[sp.Rational(v.numerator, v.denominator) for v in r]
                         for r in trial_rows]
                    ).T
                    has_newest_relation = any(v[-1] != 0 for v in trial.nullspace())
                    self.assertEqual(relation is not None, has_newest_relation)
                    if relation is not None:
                        rel = sp.Matrix([sp.Rational(v.numerator, v.denominator) for v in relation])
                        self.assertEqual(trial * rel, sp.zeros(width, 1))
                        self.assertEqual(relation[-1], Fraction(1))
                        break
                    rows.append(nxt)
                else:
                    self.fail("Krylov sequence did not become dependent within dimension")


if __name__ == "__main__":
    unittest.main()
