import sys
import unittest
import random
import sympy as sp
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from toeplitz_recurrences import m_reduce, fiduccia_pol_squarings


class MReduceTests(unittest.TestCase):
    def test_strips_terminal_pairs_in_offset_form(self):
        self.assertEqual(m_reduce([0, 1, 2], [0, 1, 2]), ((0,), (0,)))

    def test_leaves_minimal_signature_unchanged(self):
        self.assertEqual(m_reduce([1, -1], [1, 0]), ((1, -1), (1, 0)))


class FiducciaTests(unittest.TestCase):

    def test_supplied_symbolic_next_term(self):
        p1, p2, p3 = sp.symbols('p1 p2 p3')
        got = fiduccia_pol_squarings([1, 23, 3], [p1, p2, p3], 3)
        self.assertEqual(sp.expand(got - (3*p1 + 23*p2 + p3)), 0)

    def test_supplied_symbolic_three_term_window(self):
        p1, p2, p3 = sp.symbols('p1 p2 p3')
        got = fiduccia_pol_squarings(
            [1, 23, 3], [p1, p2, p3], 3, term_count=3
        )
        expected = [
            3*p1 + 23*p2 + p3,
            3*p1**2 + 3*p2 + 23*p1*p2 + 23*p3 + p1*p3,
            3*p1**3 + 6*p1*p2 + 23*p1**2*p2 + 23*p2**2
            + 3*p3 + 23*p1*p3 + p1**2*p3 + p2*p3,
        ]
        self.assertTrue(all(sp.expand(a-b) == 0 for a, b in zip(got, expected)))

    def test_random_modular_cases_against_companion_power(self):
        rng = random.Random(20260906)

        def matmul(A, B, m):
            nr, nk, nc = len(A), len(B), len(B[0])
            return [[sum(A[i][k] * B[k][j] for k in range(nk)) % m
                     for j in range(nc)] for i in range(nr)]

        def matpow(A, e, m):
            n = len(A)
            R = [[int(i == j) for j in range(n)] for i in range(n)]
            while e:
                if e & 1:
                    R = matmul(R, A, m)
                e >>= 1
                if e:
                    A = matmul(A, A, m)
            return R

        def reference(initial, recurrence, n, m):
            d = len(initial)
            if n < d:
                return initial[n] % m
            C = [[0] * d for _ in range(d)]
            C[0] = [c % m for c in recurrence]
            for i in range(1, d):
                C[i][i-1] = 1
            P = matpow(C, n - (d - 1), m)
            state = [[initial[d-1-i] % m] for i in range(d)]
            return matmul(P, state, m)[0][0]

        for _ in range(250):
            d = rng.randint(2, 8)
            mod = rng.choice([97, 101, 123, 257])
            initial = [rng.randint(-20, 20) for _ in range(d)]
            recurrence = [rng.randint(-10, 10) for _ in range(d)]
            n = rng.randint(d, 10**6)
            got = fiduccia_pol_squarings(initial, recurrence, n, modulus=mod)
            self.assertEqual(got, reference(initial, recurrence, n, mod))

    def test_random_exact_cases_against_sequential_reference(self):
        rng = random.Random(20260907)
        for _ in range(100):
            d = rng.randint(2, 8)
            initial = [rng.randint(-5, 5) for _ in range(d)]
            recurrence = [rng.randint(-3, 3) for _ in range(d)]
            n = rng.randint(d, 40)
            seq = initial[:]
            while len(seq) <= n + 2:
                seq.append(sum(recurrence[i] * seq[-1-i] for i in range(d)))
            got = fiduccia_pol_squarings(initial, recurrence, n, term_count=3)
            self.assertEqual(got, seq[n:n+3])

    def test_supplied_numeric_sequence(self):
        got = [fiduccia_pol_squarings([1, 23, 3], [4, 5, 6], i) for i in range(7)]
        self.assertEqual(got, [1, 23, 3, 133, 685, 3423, 17915])

    def test_moderate_terms_match_sequential_reference(self):
        initial = [2, -1, 5, 4]
        recurrence = [3, -2, 7, 1]
        seq = initial[:]
        while len(seq) < 41:
            seq.append(sum(recurrence[i] * seq[-1-i] for i in range(4)))
        got = [fiduccia_pol_squarings(initial, recurrence, i) for i in range(41)]
        self.assertEqual(got, seq)

    def test_supplied_modular_n100_window(self):
        self.assertEqual(
            fiduccia_pol_squarings([1, 23, 3], [111, 3332, 12], 100,
                                   modulus=123, term_count=3),
            [36, 11, 30],
        )

    def test_supplied_modular_far_windows(self):
        cases = [
            (10**7, [84, 62, 93]),
            (10**100, [117, 20, 117]),
            (10**1000, [24, 5, 33]),
            (10**10000, [105, 56, 42]),
        ]
        for n, expected in cases:
            with self.subTest(n_bit_length=n.bit_length()):
                self.assertEqual(
                    fiduccia_pol_squarings([1, 23, 3], [111, 3332, 12], n,
                                           modulus=123, term_count=3),
                    expected,
                )

    def test_supplied_modular_order_five_window(self):
        self.assertEqual(
            fiduccia_pol_squarings(
                [1, 23, 3, -11, 133], [4, 5, 6, -3, 6], 10**7,
                modulus=123, term_count=5,
            ),
            [100, 57, 68, 107, 66],
        )

    def test_window_starting_inside_initial_block_extends_correctly(self):
        self.assertEqual(
            fiduccia_pol_squarings([1, 23, 3], [4, 5, 6], 1, term_count=5),
            [23, 3, 133, 685, 3423],
        )

    def test_rejects_mismatched_lengths(self):
        with self.assertRaises(ValueError):
            fiduccia_pol_squarings([1, 2], [3], 10)

    def test_rejects_negative_index(self):
        with self.assertRaises(ValueError):
            fiduccia_pol_squarings([1, 2], [3, 4], -1)

    def test_rejects_invalid_modulus(self):
        with self.assertRaises(ValueError):
            fiduccia_pol_squarings([1, 2], [3, 4], 10, modulus=1)

    def test_rejects_zero_term_count(self):
        with self.assertRaises(ValueError):
            fiduccia_pol_squarings([1, 2], [3, 4], 10, term_count=0)


if __name__ == "__main__":
    unittest.main()
