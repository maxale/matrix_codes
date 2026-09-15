"""Paper-facing examples for the Toeplitz recurrence algorithms SageMath code."""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sage.all import PolynomialRing, QQ, binomial, ceil, floor
from toeplitz_recurrences import (
    fiduccia_pol_squarings,
    lin_rec_for_dtm,
    lin_rec_for_ptm,
    term_for_dtm,
)


def _same(a, b):
    if hasattr(a, "parent"):
        try:
            return a == a.parent()(b)
        except (TypeError, ValueError):
            pass
    return a == b


def _tridiagonal_char(size, lam, q):
    return sum(
        (-1) ** j * binomial(size - j, j) * q ** j * lam ** (size - 2 * j)
        for j in range(size // 2 + 1)
    )


def run_toeplitz_paper_examples():
    R = PolynomialRing(QQ, names=("x", "y", "lam"))
    x, y, lam = R.gens()

    det22 = lin_rec_for_dtm(2, 2, scalar_recurrence=True, verbose=False)
    per22 = lin_rec_for_ptm(2, 2, scalar_recurrence=True, verbose=False)

    det33 = lin_rec_for_dtm(3, 3, scalar_recurrence=False, verbose=False)
    nonzero33 = len(det33["transfer_matrix_expression"].dict())

    sparse_diags = [-y, 0, lam, -x]
    sparse21 = lin_rec_for_dtm(
        2, 1,
        diagonal_values=sparse_diags,
        scalar_recurrence=True,
        scalar_recurrence_method="krylov",
        characteristic_polynomial_variable="z",
        verbose=False,
    )
    sparse_terms = [
        term_for_dtm(
            2, 1, n,
            diagonal_values=sparse_diags,
            scalar_recurrence_method="krylov",
        )
        for n in range(13)
    ]
    sparse_formula = [
        sum(
            (-1) ** j * binomial(n - 2 * j, j) * (x ** 2 * y) ** j * lam ** (n - 3 * j)
            for j in range(n // 3 + 1)
        )
        for n in range(13)
    ]

    parity_diags = [-y, 0, lam, 0, -x]
    parity_terms = [
        term_for_dtm(
            2, 2, n,
            diagonal_values=parity_diags,
            scalar_recurrence_method="krylov",
        )
        for n in range(13)
    ]
    parity_formula = []
    for n in range(13):
        p = int(ceil(n / 2))
        q = int(floor(n / 2))
        parity_formula.append(
            _tridiagonal_char(p, lam, x * y) * _tridiagonal_char(q, lam, x * y)
        )

    krylov22 = lin_rec_for_dtm(
        2, 2,
        scalar_recurrence=True,
        scalar_recurrence_method="krylov",
        verbose=False,
    )
    position22 = lin_rec_for_dtm(
        2, 2,
        position_dependent=True,
        scalar_recurrence=False,
        verbose=False,
    )

    recurrence_parent = sparse21["scalar_recurrence_polynomial"].parent()
    z = recurrence_parent.gen()
    expected_sparse_recurrence = z ** 3 - recurrence_parent(lam) * z ** 2 + recurrence_parent(x ** 2 * y)

    result = {
        "pentadiagonal_determinant": {
            "transfer_matrix": det22["transfer_matrix_expression"],
            "characteristic_polynomial": det22["characteristic_polynomial"],
        },
        "pentadiagonal_permanent": {
            "transfer_matrix": per22["transfer_matrix_expression"],
            "characteristic_polynomial": per22["characteristic_polynomial"],
        },
        "balanced_33": {
            "state_count": det33["state_order"],
            "nonzero_transition_count": nonzero33,
        },
        "sparse_21_spectrum": {
            "scalar_recurrence_polynomial": sparse21["scalar_recurrence_polynomial"],
            "expected_polynomial": expected_sparse_recurrence,
            "terms_0_through_12": sparse_terms,
            "explicit_formula_0_through_12": sparse_formula,
            "explicit_formula_check": all(_same(a, b) for a, b in zip(sparse_terms, sparse_formula)),
        },
        "offset_plus_minus_2_parity_spectrum": {
            "terms_0_through_12": parity_terms,
            "parity_formula_0_through_12": parity_formula,
            "parity_split_check": all(_same(a, b) for a, b in zip(parity_terms, parity_formula)),
        },
        "krylov_22": {
            "observable_order": krylov22["observable_order"],
            "scalar_recurrence_polynomial": krylov22["scalar_recurrence_polynomial"],
        },
        "position_dependent_22": {
            "states": position22["states"],
            "transfer_matrix_expression": position22["transfer_matrix_expression"],
        },
        "fiduccia_example": fiduccia_pol_squarings(
            [1, 23, 3], [111, 3332, 12], 100, modulus=123, term_count=3
        ),
    }

    assert result["balanced_33"] == {"state_count": 20, "nonzero_transition_count": 50}
    assert result["sparse_21_spectrum"]["scalar_recurrence_polynomial"] == expected_sparse_recurrence
    assert result["sparse_21_spectrum"]["explicit_formula_check"]
    assert result["offset_plus_minus_2_parity_spectrum"]["parity_split_check"]
    assert result["fiduccia_example"] == [36, 11, 30]
    return result


if __name__ == "__main__":
    data = run_toeplitz_paper_examples()
    print("Paper examples: PASS")
    print(data)
