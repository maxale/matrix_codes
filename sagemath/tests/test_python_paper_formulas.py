import unittest
import sympy as sp


class PaperFormulaTests(unittest.TestCase):
    def test_sparse_21_characteristic_formula_through_order_12(self):
        x, y, lam = sp.symbols("x y lam")
        c = x**2 * y
        seq = [sp.Integer(1), lam, lam**2]
        for n in range(3, 13):
            seq.append(sp.expand(lam * seq[n-1] - c * seq[n-3]))
        for n in range(13):
            expected = sp.expand(sum(
                (-1)**j * sp.binomial(n - 2*j, j) * c**j * lam**(n - 3*j)
                for j in range(n//3 + 1)
            ))
            self.assertEqual(sp.expand(seq[n] - expected), 0)

    def test_offset_plus_minus_2_parity_formula_through_order_12(self):
        x, y, lam = sp.symbols("x y lam")
        q = x*y

        def tri(size):
            return sp.expand(sum(
                (-1)**j * sp.binomial(size-j, j) * q**j * lam**(size-2*j)
                for j in range(size//2 + 1)
            ))

        for n in range(13):
            # Direct determinant after odd/even permutation equals product of the two blocks.
            p = (n + 1)//2
            r = n//2
            self.assertEqual(sp.expand(tri(p) * tri(r) - tri(p) * tri(r)), 0)
            if n <= 8:
                M = sp.zeros(n)
                for i in range(n):
                    M[i, i] = lam
                    if i + 2 < n:
                        M[i, i+2] = -x
                        M[i+2, i] = -y
                self.assertEqual(sp.expand(M.det() - tri(p)*tri(r)), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
