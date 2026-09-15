import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from sage.all import PolynomialRing, QQ, zero_matrix
from toeplitz_recurrences import lin_rec_for_dtm, lin_rec_for_ptm

R5 = PolynomialRing(QQ, names=("a_m2", "a_m1", "a_0", "a_p1", "a_p2"))
a_m2, a_m1, a_0, a_p1, a_p2 = R5.gens()
P5 = PolynomialRing(R5, "x")
x = P5.gen()

R3 = PolynomialRing(QQ, names=("A", "B", "C"))
A, B, C = R3.gens()
P3 = PolynomialRing(R3, "x")
y = P3.gen()


def same_poly(p, q):
    return p == p.parent()(q)


def same_element(p, q):
    return p == p.parent()(q)


base22 = [a_m2, a_m1, a_0, a_p1, a_p2]
k11 = lin_rec_for_dtm(1, 1, diagonal_values=[a_m1, a_0, a_p1], scalar_recurrence=True,
                      scalar_recurrence_method="krylov", verbose=False)
k22 = lin_rec_for_dtm(2, 2, diagonal_values=base22, scalar_recurrence=True,
                      scalar_recurrence_method="krylov", verbose=False)
kp22 = lin_rec_for_ptm(2, 2, diagonal_values=base22, scalar_recurrence=True,
                       scalar_recurrence_method="krylov", verbose=False)
ks22 = lin_rec_for_dtm(2, 2, diagonal_values=[C, B, A, B, C],
                       scalar_recurrence=True,
                       scalar_recurrence_method="krylov", verbose=False)
skipped = lin_rec_for_dtm(2, 2, diagonal_values=base22, scalar_recurrence=False,
                          scalar_recurrence_method="krylov", verbose=False)
pos22 = lin_rec_for_dtm(2, 2, position_dependent=True,
                        scalar_recurrence=False,
                        scalar_recurrence_method="krylov", verbose=False)

expected22 = (
    x**6 - a_0*x**5 + (a_m1*a_p1 - a_m2*a_p2)*x**4
    + (-a_m2*a_p1**2 - a_m1**2*a_p2 + 2*a_m2*a_0*a_p2)*x**3
    + (a_m2*a_m1*a_p1*a_p2 - a_m2**2*a_p2**2)*x**2
    - a_m2**2*a_0*a_p2**2*x + a_m2**3*a_p2**3
)
expected_sym5 = (
    y**5 - (A-C)*y**4 - (A*C-B**2)*y**3
    + C*(A*C-B**2)*y**2 + C**3*(A-C)*y - C**5
)

assert k11["scalar_recurrence_method"] == "krylov"
assert k11["observable_order"] == 2
assert same_poly(k11["scalar_recurrence_polynomial"], x**2-a_0*x+a_m1*a_p1)
assert k22["state_ordering"] == "Paper first-discovery order"
assert k22["observable_coordinate"] == 0
assert k22["observable_order"] == 6
assert same_poly(k22["scalar_recurrence_polynomial"], expected22)
assert same_element(
    k22["observable_matrix"].change_ring(R5).det(),
    -a_p2**5 * a_m2**7 * (a_p1**2*a_m2 - a_p2*a_m1**2),
)
assert kp22["observable_order"] == 6
assert ks22["observable_order"] == 5
assert same_poly(ks22["scalar_recurrence_polynomial"], expected_sym5)
assert ks22["similarity_matrix"] is None
assert skipped["scalar_recurrence_polynomial"] is None
assert pos22["scalar_recurrence_polynomial"] is None

O = k22["observable_matrix"]
M = k22["transfer_matrix_expression"]
K = k22["krylov_companion_matrix"]
assert O*M - K*O == zero_matrix(O.base_ring(), 6, 6)
assert k22["similarity_matrix"] == O

print("Krylov Sage tests: PASS")
