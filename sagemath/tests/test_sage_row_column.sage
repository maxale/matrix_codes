import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from sage.all import PolynomialRing, QQ, binomial, matrix
from toeplitz_recurrences import laplace_d, laplace_p, lin_rec_for_dtm, lin_rec_for_ptm

R = PolynomialRing(QQ, names=("a_m2", "a_m1", "a_0", "a_p1", "a_p2", "u", "v", "w", "q"))
a_m2, a_m1, a_0, a_p1, a_p2, u, v, w, q = R.gens()
PX = PolynomialRing(R, "x")
x = PX.gen()


def same_poly(p, qexpr):
    return p == p.parent()(qexpr)


def same_matrix(p, qmat):
    return p.change_ring(R) == qmat.change_ring(R)


ld11 = laplace_d(1, 1, ((1,), (1,)), diagonal_values=[a_m1, a_0, a_p1])
lp11 = laplace_p(1, 1, ((1,), (1,)), diagonal_values=[a_m1, a_0, a_p1])
assert ld11[0] == ((1,), (1,))
assert ld11[1][0] == (a_0, ((0,), (0,)))
assert ld11[1][1] == (-a_m1, ((1, 0), (1, -1)))
assert lp11[1][1] == (a_m1, ((1, 0), (1, -1)))

base22 = [a_m2, a_m1, a_0, a_p1, a_p2]
det11 = lin_rec_for_dtm(1, 1, diagonal_values=[a_m1, a_0, a_p1], scalar_recurrence=True, verbose=False)
det22 = lin_rec_for_dtm(2, 2, diagonal_values=base22, scalar_recurrence=True, verbose=False)
per22 = lin_rec_for_ptm(2, 2, diagonal_values=base22, scalar_recurrence=True, verbose=False)
det32 = lin_rec_for_dtm(3, 2, scalar_recurrence=False, verbose=False)
per23 = lin_rec_for_ptm(2, 3, scalar_recurrence=False, verbose=False)

assert same_matrix(det11["transfer_matrix"], matrix(R, [[a_0, -a_m1], [a_p1, 0]], sparse=True))
assert same_poly(det11["characteristic_polynomial"], x**2 - a_0*x + a_m1*a_p1)

expected_states22 = [
    ((1,), (1,)),
    ((2, 1), (2, 0)),
    ((2, 1), (2, -1)),
    ((2, 0), (2, 1)),
    ((2, 0), (2, 0)),
    ((3, 2, 1), (3, 1, 0)),
]
assert det22["states"] == expected_states22
assert per22["states"] == expected_states22
assert det22["state_ordering"] == "Paper first-discovery order"

expected_det22 = matrix(R, [
    [a_0, -a_m1, a_m2, 0, 0, 0],
    [a_p1, 0, 0, -a_p2, 0, 0],
    [0, a_p1, 0, 0, -a_p2, 0],
    [a_m1, -a_m2, 0, 0, 0, 0],
    [a_0, 0, 0, 0, 0, -a_m2],
    [a_p2, 0, 0, 0, 0, 0],
], sparse=True)
expected_per22 = matrix(R, [
    [a_0, a_m1, a_m2, 0, 0, 0],
    [a_p1, 0, 0, a_p2, 0, 0],
    [0, a_p1, 0, 0, a_p2, 0],
    [a_m1, a_m2, 0, 0, 0, 0],
    [a_0, 0, 0, 0, 0, a_m2],
    [a_p2, 0, 0, 0, 0, 0],
], sparse=True)
assert same_matrix(det22["transfer_matrix"], expected_det22)
assert same_matrix(per22["transfer_matrix"], expected_per22)
assert det32["state_order"] == binomial(5, 3)
assert per23["state_order"] == binomial(5, 2)

sweet22 = (
    x**6 - a_0*x**5 + (a_m1*a_p1 - a_m2*a_p2)*x**4
    + (-a_m2*a_p1**2 - a_m1**2*a_p2 + 2*a_m2*a_0*a_p2)*x**3
    + (a_m2*a_m1*a_p1*a_p2 - a_m2**2*a_p2**2)*x**2
    - a_m2**2*a_0*a_p2**2*x + a_m2**3*a_p2**3
)
assert same_poly(det22["characteristic_polynomial"], sweet22)

custom11 = lin_rec_for_dtm(1, 1, diagonal_values=[u, v, w], scalar_recurrence=False, verbose=False)
assert same_matrix(custom11["transfer_matrix"], matrix(R, [[v, -u], [w, 0]], sparse=True))
assert custom11["characteristic_polynomial"] is None

sparse22 = lin_rec_for_dtm(2, 2, diagonal_values=[0, q, 1, q, 0], scalar_recurrence=False, verbose=False)
expected_sparse22 = matrix(R, [
    [1, -q, 0, 0, 0, 0],
    [q, 0, 0, 0, 0, 0],
    [0, q, 0, 0, 0, 0],
    [q, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
], sparse=True)
assert same_matrix(sparse22["transfer_matrix"], expected_sparse22)

one_param22 = lin_rec_for_dtm(2, 2, diagonal_values=[q**2, q, 1, q, q**2], scalar_recurrence=True, verbose=False)
assert one_param22["characteristic_polynomial"] is not None

assert lin_rec_for_dtm(1, 1, scalar_recurrence="auto", max_characteristic_order=0, verbose=False)["characteristic_polynomial"] is None
assert lin_rec_for_dtm(1, 1, scalar_recurrence=True, max_characteristic_order=0, verbose=False)["characteristic_polynomial"] is not None

# Default nonautonomous output now uses exact polynomial indeterminates for
# the finitely many shifted entries appearing in Q[0].
pos22 = lin_rec_for_dtm(2, 2, position_dependent=True, scalar_recurrence=False, verbose=False)
posp22 = lin_rec_for_ptm(2, 2, position_dependent=True, scalar_recurrence=False, verbose=False)
assert pos22["state_order"] == 6
assert pos22["characteristic_polynomial"] is None
assert posp22["characteristic_polynomial"] is None
assert pos22["transfer_matrix_expression"][0, 1] == -posp22["transfer_matrix_expression"][0, 1]
assert str(pos22["transfer_matrix_expression"][0, 0]) == "a_0_at_k_0"
assert str(pos22["transfer_matrix_expression"][0, 1]) == "-a_m1_at_k_m1"
assert str(pos22["transfer_matrix_expression"][0, 2]) == "a_m2_at_k_m2"
assert str(pos22["transfer_matrix_expression"][4, 0]) == "a_0_at_k_p1"
q5 = pos22["transfer_matrix_function"](5)
assert str(q5[0, 1]) == "-a_m1_at_k_p4"
assert str(q5[4, 0]) == "a_0_at_k_p6"

# A custom position-dependent entry callback receives exact integer positions.
def custom_entry(offset, position):
    return 1000 * offset + position

custom_pos = lin_rec_for_dtm(
    1, 1, position_dependent=True,
    band_entry_function=custom_entry,
    scalar_recurrence=False, verbose=False,
)
assert custom_pos["transfer_matrix_expression"][0, 0] == custom_entry(0, 0)
assert custom_pos["transfer_matrix_expression"][0, 1] == -custom_entry(-1, -1)
assert custom_pos["transfer_matrix_function"](5)[0, 1] == -custom_entry(-1, 4)

print("row-column Sage tests: PASS")
