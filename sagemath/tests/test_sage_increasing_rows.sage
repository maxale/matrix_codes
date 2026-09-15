import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from sage.all import PolynomialRing, QQ, matrix
from sage.matrix.matrix_misc import permanental_minor_polynomial
from toeplitz_recurrences import (
    lin_rec_for_dtm_increasing_rows,
    lin_rec_for_ptm_increasing_rows,
)

R = PolynomialRing(QQ, names=("a_m2", "a_m1", "a_0", "a_p1", "a_p2"))
a_m2, a_m1, a_0, a_p1, a_p2 = R.gens()
PX = PolynomialRing(R, "x")
x = PX.gen()


def same_poly(p, q):
    return p == p.parent()(q)


def toeplitz_matrix(m1, m2, values, n):
    return matrix(n, n, lambda i, j: values[j-i+m1] if -m1 <= j-i <= m2 else 0)


def sequence(quantity, m1, m2, values, nmax):
    out = [R.one()]
    for n in range(1, nmax + 1):
        M = toeplitz_matrix(m1, m2, values, n)
        value = M.det() if quantity == "Determinant" else permanental_minor_polynomial(M, permanent_only=True)
        out.append(value)
    return out


def polynomial_annihilates(poly, seq):
    degree = int(poly.degree())
    if len(seq) < degree + 1:
        return False
    for n in range(len(seq) - degree):
        value = sum(poly[j] * seq[n+j] for j in range(degree + 1))
        if value != 0:
            return False
    return True


d11 = [a_m1, a_0, a_p1]
df = lin_rec_for_dtm_increasing_rows(1, 1, diagonal_values=d11, ell_range="forward", verbose=False)
dc = lin_rec_for_dtm_increasing_rows(1, 1, diagonal_values=d11, ell_range="centered", verbose=False)
pf = lin_rec_for_ptm_increasing_rows(1, 1, diagonal_values=d11, ell_range="forward", verbose=False)
pc = lin_rec_for_ptm_increasing_rows(1, 1, diagonal_values=d11, ell_range="centered", verbose=False)

assert same_poly(df["annihilating_polynomial"], x**2 - a_0*x + a_m1*a_p1)
assert same_poly(pf["annihilating_polynomial"], x**2 - a_0*x - a_m1*a_p1)
assert same_poly(df["annihilating_polynomial"], dc["annihilating_polynomial"])
assert same_poly(pf["annihilating_polynomial"], pc["annihilating_polynomial"])
assert df["coefficient_minor_orders"] == [1, 2, 3]
assert dc["coefficient_minor_orders"] == [0, 1, 2]

d21 = [a_m2, a_m1, a_0, a_p1]
d21f = lin_rec_for_dtm_increasing_rows(2, 1, diagonal_values=d21, ell_range="forward", verbose=False)
d21c = lin_rec_for_dtm_increasing_rows(2, 1, diagonal_values=d21, ell_range="centered", verbose=False)
p21f = lin_rec_for_ptm_increasing_rows(2, 1, diagonal_values=d21, ell_range="forward", verbose=False)
p21c = lin_rec_for_ptm_increasing_rows(2, 1, diagonal_values=d21, ell_range="centered", verbose=False)
assert same_poly(d21f["annihilating_polynomial"], x**3 - a_0*x**2 + a_m1*a_p1*x - a_m2*a_p1**2)
assert same_poly(p21f["annihilating_polynomial"], x**3 - a_0*x**2 - a_m1*a_p1*x - a_m2*a_p1**2)
assert same_poly(d21f["annihilating_polynomial"], d21c["annihilating_polynomial"])
assert same_poly(p21f["annihilating_polynomial"], p21c["annihilating_polynomial"])

direct_d21 = sequence("Determinant", 2, 1, d21, 6)
direct_p21 = sequence("Permanent", 2, 1, d21, 5)
assert polynomial_annihilates(d21c["annihilating_polynomial"], direct_d21)
assert polynomial_annihilates(p21c["annihilating_polynomial"], direct_p21)

tr = lin_rec_for_dtm_increasing_rows(
    1, 2,
    diagonal_values=[a_m1, a_0, a_p1, a_p2],
    ell_range="centered", verbose=False,
)
assert tr["transposed"] is True
assert same_poly(tr["annihilating_polynomial"], x**3 - a_0*x**2 + a_m1*a_p1*x - a_p2*a_m1**2)

print("increasing-rows Sage tests: PASS")
