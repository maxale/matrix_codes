import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from sage.all import PolynomialRing, QQ, matrix
from sage.matrix.matrix_misc import permanental_minor_polynomial
import toeplitz_recurrences as tr
from toeplitz_recurrences import (
    fiduccia_pol_squarings,
    term_for_dtm,
    term_for_ptm,
    sage_cfinite_term,
)


def toeplitz_matrix(m1, m2, diagonals, n):
    return matrix(n, n, lambda i, j: diagonals[j-i+m1] if -m1 <= j-i <= m2 else 0)


# SageMath requires polynomial variable names to begin with a letter.
_tmp_name = tr._fresh_polynomial_variable_name()
assert _tmp_name[0].isalpha()
PolynomialRing(QQ, _tmp_name)


# Symbolic preservation of the supplied order-three example.
R = PolynomialRing(QQ, names=('p1', 'p2', 'p3'))
p1, p2, p3 = R.gens()
assert fiduccia_pol_squarings([1, 23, 3], [p1, p2, p3], 3) == 3*p1 + 23*p2 + p3

# Native Sage CFiniteSequences is an independent reference over QQ/ZZ.
initial = [2, -1, 5, 4]
recurrence = [3, -2, 7, 1]
for n in [0, 1, 4, 10, 40, 137]:
    assert fiduccia_pol_squarings(initial, recurrence, n) == sage_cfinite_term(initial, recurrence, n)

# Toeplitz remote determinant/permanent integration.
d50 = term_for_dtm(1, 1, 50, diagonal_values=[2, 1, 3])
p50 = term_for_ptm(1, 1, 50, diagonal_values=[2, 1, 3])
assert d50 == sage_cfinite_term([1, 1], [1, -6], 50)
assert p50 == sage_cfinite_term([1, 1], [1, 6], 50)

assert term_for_dtm(1, 1, 10**7, diagonal_values=[2, 1, 3], modulus=101) == fiduccia_pol_squarings([1, 1], [1, -6], 10**7, modulus=101)
assert term_for_ptm(1, 1, 10**7, diagonal_values=[2, 1, 3], modulus=101) == fiduccia_pol_squarings([1, 1], [1, 6], 10**7, modulus=101)

assert term_for_dtm(1, 1, 30, diagonal_values=[2, 1, 3], term_count=3) == [
    sage_cfinite_term([1, 1], [1, -6], n) for n in range(30, 33)
]

d22_diags = [2, -1, 3, 4, 1]
d22_direct = toeplitz_matrix(2, 2, d22_diags, 8).det()
assert term_for_dtm(2, 2, 8, diagonal_values=d22_diags) == d22_direct

p21_diags = [1, 2, 3, 4]
p21_matrix = toeplitz_matrix(2, 1, p21_diags, 6)
p21_direct = permanental_minor_polynomial(p21_matrix, permanent_only=True)
assert term_for_ptm(2, 1, 6, diagonal_values=p21_diags) == p21_direct

# Characteristic and Krylov scalarizations feed the same Fiduccia evaluator.
assert term_for_dtm(1, 1, 50, diagonal_values=[2, 1, 3], scalar_recurrence_method='krylov') == d50

try:
    term_for_dtm(1, 1, 10, diagonal_values=[2, 1, 3], position_dependent=True)
    raise AssertionError('position-dependent remote evaluation was not rejected')
except ValueError:
    pass

print('Fiduccia/Toeplitz Sage integration tests: PASS')
