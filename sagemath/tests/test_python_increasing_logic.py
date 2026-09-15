import sys
import unittest
from contextlib import contextmanager
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
import toeplitz_recurrences as tr


@contextmanager
def sympy_minor_backend():
    original = tr._minor_value

    def minor(rows, cols, quantity, entry):
        if len(rows) != len(cols):
            raise RuntimeError
        if not rows:
            return sp.Integer(1)
        M = sp.Matrix([[entry(cols[j] - rows[i], 0) for j in range(len(cols))]
                       for i in range(len(rows))])
        return M.det() if quantity == 'Determinant' else M.per()

    tr._minor_value = minor
    try:
        yield
    finally:
        tr._minor_value = original


def annihilator(m1_in, m2_in, quantity, centered=False):
    names = {
        s: sp.Symbol(('a_m' + str(-s)) if s < 0 else ('a_p' + str(s) if s > 0 else 'a_0'))
        for s in range(-m1_in, m2_in + 1)
    }
    m1, m2 = m1_in, m2_in
    if m1 < m2:
        original = names
        m1, m2 = m2, m1
        values = {s: original.get(-s, 0) for s in range(-m1, m2 + 1)}
    else:
        values = names
    entry = lambda s, t: values.get(s, 0)
    d = __import__('math').comb(m1 + m2, m1)
    states = tr._boundary_states(m1, m2)
    ell = list(range(-m1, d - m1 + 1)) if centered else list(range(d + 1))
    k0 = 4 * (m1 + m2 + d + 2)
    with sympy_minor_backend():
        B = sp.Matrix([
            [tr._increasing_coefficient(m1, m2, e, state, quantity, entry, k0)
             for state in states]
            for e in ell
        ])
    basis = B.T.nullspace()
    assert basis
    v = basis[0]
    nz = [i for i, z in enumerate(v) if z != 0]
    v = [sp.cancel(z / v[nz[-1]]) for z in v]
    x = sp.Symbol('x')
    p = sp.expand(sum(v[i] * x ** (ell[i] - min(ell)) for i in range(len(ell))))
    return sp.factor(p), [m1 + e for e in ell]


class IncreasingRowsLogicTests(unittest.TestCase):
    def test_tridiagonal_determinant_and_permanent(self):
        x = sp.Symbol('x'); a0=sp.Symbol('a_0'); am1=sp.Symbol('a_m1'); ap1=sp.Symbol('a_p1')
        d, orders = annihilator(1, 1, 'Determinant')
        p, _ = annihilator(1, 1, 'Permanent')
        self.assertEqual(sp.expand(d - (x**2-a0*x+am1*ap1)), 0)
        self.assertEqual(sp.expand(p - (x**2-a0*x-am1*ap1)), 0)
        self.assertEqual(orders, [1,2,3])

    def test_forward_centered_agree_for_21(self):
        x=sp.Symbol('x'); a0=sp.Symbol('a_0'); am1=sp.Symbol('a_m1'); am2=sp.Symbol('a_m2'); ap1=sp.Symbol('a_p1')
        df,_=annihilator(2,1,'Determinant',False)
        dc,orders=annihilator(2,1,'Determinant',True)
        pf,_=annihilator(2,1,'Permanent',False)
        pc,_=annihilator(2,1,'Permanent',True)
        self.assertEqual(sp.expand(df-dc),0)
        self.assertEqual(sp.expand(pf-pc),0)
        self.assertEqual(sp.expand(df-(x**3-a0*x**2+am1*ap1*x-am2*ap1**2)),0)
        self.assertEqual(sp.expand(pf-(x**3-a0*x**2-am1*ap1*x-am2*ap1**2)),0)
        self.assertEqual(orders, [0,1,2,3])

    def test_transposed_12_uses_original_diagonals(self):
        x=sp.Symbol('x'); a0=sp.Symbol('a_0'); am1=sp.Symbol('a_m1'); ap1=sp.Symbol('a_p1'); ap2=sp.Symbol('a_p2')
        got,_=annihilator(1,2,'Determinant',True)
        expected=x**3-a0*x**2+am1*ap1*x-ap2*am1**2
        self.assertEqual(sp.expand(got-expected),0)


if __name__ == '__main__':
    unittest.main()
