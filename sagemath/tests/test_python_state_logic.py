import math
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
import toeplitz_recurrences as tr


def raw_closure(m1, m2, quantity):
    k = sp.symbols('k')
    a = {
        s: sp.Symbol(('a_m' + str(-s)) if s < 0 else ('a_p' + str(s) if s > 0 else 'a_0'))
        for s in range(-m1, m2 + 1)
    }
    entry = lambda s, t: a.get(s, 0)
    seen = [((1,), (1,))]
    answer = []
    i = 0
    while i < len(seen):
        data = tr._laplace_expansion(m1, m2, seen[i], quantity, entry, k)
        answer.append(data)
        for _, target in data[1]:
            tp = tr.m_reduce(*tr._final_signature(target))
            if tp not in seen:
                seen.append(tp)
        i += 1
    finals = [d[0] for d in answer]
    shifted = [
        (tuple(v - 1 for v in rows), tuple(v - 1 for v in cols))
        for rows, cols in finals
    ]
    M = sp.zeros(len(answer))
    for r, (_, expansion) in enumerate(answer):
        for coeff, target in expansion:
            M[r, shifted.index(target)] = sp.expand(coeff)
    return finals, M, a

class ScalarPolicyTests(unittest.TestCase):
    def test_invalid_scalar_policy_is_rejected_before_backend_dispatch(self):
        with self.assertRaises(ValueError):
            tr._validate_scalar_policy('banana')

    def test_valid_scalar_policies_are_accepted(self):
        for value in (True, False, 'auto', 'automatic'):
            with self.subTest(value=value):
                self.assertIsNone(tr._validate_scalar_policy(value))


class RowColumnLogicTests(unittest.TestCase):
    def test_position_dependent_template_uses_exact_integer_shifts(self):
        closure = tr._row_column_closure_data(
            2, 2, "Determinant", tr._formal_band_entry, 0
        )
        self.assertEqual(closure["entries"][(0, 0)], tr._BandTerm(0, 0, 1))
        self.assertEqual(closure["entries"][(0, 1)], tr._BandTerm(-1, -1, -1))
        self.assertEqual(closure["entries"][(0, 2)], tr._BandTerm(-2, -2, 1))
        self.assertEqual(closure["entries"][(4, 0)], tr._BandTerm(0, 1, 1))

    def test_pentadiagonal_state_order_and_determinant_matrix(self):
        states, M, a = raw_closure(2, 2, 'Determinant')
        expected_states = [
            ((1,), (1,)),
            ((2, 1), (2, 0)),
            ((2, 1), (2, -1)),
            ((2, 0), (2, 1)),
            ((2, 0), (2, 0)),
            ((3, 2, 1), (3, 1, 0)),
        ]
        self.assertEqual(states, expected_states)
        expected = sp.Matrix([
            [a[0], -a[-1], a[-2], 0, 0, 0],
            [a[1], 0, 0, -a[2], 0, 0],
            [0, a[1], 0, 0, -a[2], 0],
            [a[-1], -a[-2], 0, 0, 0, 0],
            [a[0], 0, 0, 0, 0, -a[-2]],
            [a[2], 0, 0, 0, 0, 0],
        ])
        self.assertEqual(M, expected)

    def test_pentadiagonal_characteristic_polynomial_matches_paper_formula(self):
        _, M, a = raw_closure(2, 2, 'Determinant')
        x = sp.symbols('x')
        got = sp.expand(M.charpoly(x).as_expr())
        expected = (
            x**6 - a[0]*x**5 + (a[-1]*a[1] - a[-2]*a[2])*x**4
            + (-a[-2]*a[1]**2 - a[-1]**2*a[2] + 2*a[-2]*a[0]*a[2])*x**3
            + (a[-2]*a[-1]*a[1]*a[2] - a[-2]**2*a[2]**2)*x**2
            - a[-2]**2*a[0]*a[2]**2*x + a[-2]**3*a[2]**3
        )
        self.assertEqual(sp.expand(got - expected), 0)

    def test_permanent_uses_same_state_graph(self):
        d_states, _, _ = raw_closure(2, 2, 'Determinant')
        p_states, P, a = raw_closure(2, 2, 'Permanent')
        self.assertEqual(p_states, d_states)
        self.assertEqual(P[0, 1], a[-1])
        self.assertEqual(P[1, 3], a[2])

    def test_asymmetric_state_counts(self):
        for m1, m2 in [(3, 2), (2, 3), (4, 1), (1, 4)]:
            with self.subTest(m1=m1, m2=m2):
                states, _, _ = raw_closure(m1, m2, 'Determinant')
                self.assertEqual(len(states), math.comb(m1 + m2, m1))


if __name__ == '__main__':
    unittest.main()
