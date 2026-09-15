# Toeplitz recurrence algorithms - SageMath

This directory contains the SageMath implementation accompanying
**Constructive recurrences for determinants and permanents of banded Toeplitz matrices**
by Max A. Alekseyev and Dmitry I. Khomovsky.

The code constructs finite recurrences for determinants and permanents of fixed-band
Toeplitz matrices, supports the increasing-rows and row-column constructions developed in
the paper, adds optional Krylov scalarization, handles the position-dependent row-column
cocycle, and provides fast distant-term evaluation once a homogeneous constant-coefficient
recurrence is known. The Sage implementation uses Sage's native matrix, polynomial, kernel,
and C-finite sequence facilities wherever they fit the algorithms directly.

Public repository: `https://github.com/maxale/matrix_codes`

## Requirements

- SageMath 10.x or later is recommended.
- No external packages are required for the Sage implementation itself.
- Python 3 can run the backend-independent Fiduccia/state tests; those tests use SymPy.
- Python 3 is also used by `tools/verify_release.py` for static release checks.

## Repository layout

```text
toeplitz_recurrences_sagemath/
|-- README.md
|-- CITATION.cff
|-- CHANGELOG.md
|-- PACKAGE_INFO.txt
|-- VERIFICATION.md
|-- MANIFEST.sha256
|-- src/
|   |-- toeplitz_recurrences.py
|   `-- toeplitz_recurrences.sage
|-- examples/
|   |-- paper_examples.sage
|   `-- compare_fiduccia.sage
|-- tests/
|   |-- run_tests.sage
|   |-- test_python_core.py
|   |-- test_python_increasing_logic.py
|   |-- test_python_package_contract.py
|   |-- test_python_paper_formulas.py
|   |-- test_python_release_contract.py
|   |-- test_python_state_logic.py
|   |-- test_python_no_symbolic_ring.py
|   |-- test_sage_row_column.sage
|   |-- test_sage_increasing_rows.sage
|   |-- test_sage_krylov.sage
|   `-- test_sage_fiduccia_integration.sage
`-- tools/
    `-- verify_release.py
```

## Quick start

From the repository root in Sage:

```sage
load("src/toeplitz_recurrences.sage")
```

or in Sage/Python code:

```python
import sys
sys.path.insert(0, "src")
from toeplitz_recurrences import *
```

The public convention is always

```text
m1 = number of subdiagonals   (offsets -m1,...,-1)
m2 = number of superdiagonals (offsets 1,...,m2)
```

and `diagonal_values` are supplied in offset order `-m1,...,m2`.

## Paper-to-code map

| Paper construction or application | Main Sage interface | Notes |
|---|---|---|
| Increasing-rows construction | `lin_rec_for_dtm_increasing_rows`, `lin_rec_for_ptm_increasing_rows` | Supports `ell_range="forward"`, `ell_range="centered"`, or an explicit level list; transparently transposes when `m1 < m2`. |
| Row-column transfer construction | `lin_rec_for_dtm`, `lin_rec_for_ptm` | Uses the paper's first-discovery state ordering. |
| Characteristic-polynomial scalarization | `scalar_recurrence=True` | Uses Sage matrix `charpoly`. |
| Observable/Krylov scalarization | `scalar_recurrence_method="krylov"` | Can produce a smaller observable annihilator. |
| Position-dependent row-column cocycle | `position_dependent=True` | Returns a position-dependent transfer; constant scalarization is intentionally disabled. |
| Distant-term evaluation | `fiduccia_pol_squarings`, `term_for_dtm`, `term_for_ptm` | Applies only after a homogeneous constant-coefficient recurrence has been obtained. |
| Sage-native recurrence comparison | `sage_cfinite_term` | Uses Sage `CFiniteSequences` as an independent reference path for rational/integer recurrences. |
| Sparse `(2,1)` Toeplitz-Hessenberg spectral example | `lin_rec_for_dtm(2,1,...)`, `term_for_dtm(2,1,...)` | `examples/paper_examples.sage` checks the characteristic recurrence and explicit finite-section formula behind the threefold spectrum. |
| Offset `+/-2` parity-split comparison | `term_for_dtm(2,2,...)` | `examples/paper_examples.sage` checks the decomposition into two tridiagonal Toeplitz blocks. |

The code does **not** currently provide dedicated implementations of the cyclic-closure
results, the proposed finite-corner-defect state extension, or the nonautonomous increasing-rows theorem. The compound/Widom similarity theorem is likewise a mathematical
identification in the paper rather than a separate runtime constructor.

## Public API and defaults

Sage uses Python keyword arguments rather than Wolfram options. `None` plays the role of
an automatic/default object for several inputs, while `scalar_recurrence="auto"` is the
explicit automatic scalarization policy.

### Row-column recurrence constructors

`lin_rec_for_dtm(m1, m2, ...)` and `lin_rec_for_ptm(m1, m2, ...)` accept:

| Keyword | Default | Meaning |
|---|---|---|
| `diagonal_values` | `None` | Constant diagonal values in offset order `-m1,...,m2`; `None` creates exact polynomial indeterminates over `QQ`. |
| `scalar_recurrence` | `"auto"` | `True` requests scalarization, `False` skips it, and `"auto"` scalarizes only when the transfer order does not exceed `max_characteristic_order`. |
| `characteristic_polynomial_variable` | `"x"` | Variable name used for the returned scalar recurrence polynomial. |
| `max_characteristic_order` | `20` | Transfer-order cutoff used only with `scalar_recurrence="auto"`. |
| `verbose` | `True` | Controls interactive printing. |
| `position_dependent` | `False` | Selects the nonautonomous row-column cocycle. |
| `band_entry_function` | `None` | Optional callable `(offset, position) -> entry` in position-dependent mode; `None` creates exact polynomial indeterminates for the finitely many shifted entries in each transfer matrix. |
| `scalar_recurrence_method` | `"characteristic_polynomial"` | Choose `"characteristic_polynomial"` or `"krylov"`. |

`position_dependent=True` cannot be combined with explicit constant `diagonal_values`.
Constant characteristic/Krylov scalarization is intentionally disabled in position-dependent
mode.

### Increasing-rows constructors

`lin_rec_for_dtm_increasing_rows(m1, m2, ...)` and
`lin_rec_for_ptm_increasing_rows(m1, m2, ...)` accept:

| Keyword | Default | Meaning |
|---|---|---|
| `diagonal_values` | `None` | Constant diagonal values in offset order `-m1,...,m2`. |
| `ell_range` | `"forward"` | Use `"forward"`, `"centered"`, or an explicit strictly increasing list of the required `d+1` levels. |
| `characteristic_polynomial_variable` | `"x"` | Variable name used for the returned annihilating polynomial. |
| `verbose` | `True` | Controls interactive printing. |

### One-step Laplace helpers

`laplace_d(m1, m2, state, ...)` and `laplace_p(m1, m2, state, ...)` accept:

| Keyword | Default | Meaning |
|---|---|---|
| `diagonal_values` | `None` | Constant diagonal values in offset order `-m1,...,m2`. |
| `position_dependent` | `False` | Select constant or position-dependent band entries. |
| `band_entry_function` | `None` | Optional callable `(offset, position) -> entry` in position-dependent mode. |

### Distant-term evaluation

`fiduccia_pol_squarings(initial, recurrence, n, ...)` accepts:

| Keyword | Default | Meaning |
|---|---|---|
| `modulus` | `None` | Exact arithmetic in the supplied coefficient parent by default; an integer greater than 1 selects arithmetic modulo that integer. |
| `term_count` | `1` | Return one term by default; a larger value returns consecutive terms starting at zero-based index `n`. |

`term_for_dtm(m1, m2, n, ...)` and `term_for_ptm(m1, m2, n, ...)` accept:

| Keyword | Default | Meaning |
|---|---|---|
| `diagonal_values` | `None` | Constant Toeplitz diagonal values in offset order `-m1,...,m2`. |
| `scalar_recurrence_method` | `"characteristic_polynomial"` | Choose the full transfer characteristic polynomial or the observable Krylov recurrence. |
| `modulus` | `None` | Exact arithmetic or integer arithmetic modulo the supplied integer. |
| `term_count` | `1` | One distant term or a consecutive window of terms. |
| `verbose` | `False` | Controls printing by the internal recurrence constructor. |
| `position_dependent` | `False` | Must remain `False`; distant-term evaluation requires a homogeneous constant-coefficient recurrence. |

`sage_cfinite_term(initial, recurrence, n, term_count=1)` is a Sage-native reference path
for integer/rational recurrences. It reverses the public recurrence list because Sage's
`CFiniteSequences` interface uses the opposite coefficient order.

## Sage functionality used directly

The implementation deliberately relies on existing Sage functionality rather than
reimplementing it:

- dense and sparse Sage matrices for transfer systems;
- matrix `.det()` and Sage's permanent helper for minors and direct checks;
- matrix `.charpoly(...)` for characteristic-polynomial scalarization;
- exact polynomial rings and fraction fields for symbolic transfers;
- matrix kernels for increasing-rows relations;
- incremental fraction-free row elimination for observable/Krylov scalarization;
- `CFiniteSequences(QQ).from_recurrence(...)` as an independent recurrence evaluator.

The row-column state closure and increasing-rows boundary combinatorics remain explicit
because they are the paper's constructive algorithms. The hardcoded Fiduccia/Khomovsky
squaring formulas also remain explicit so they can be compared directly with Sage's native
C-finite machinery.

## Exact coefficient-ring policy

The constant-coefficient algorithms deliberately avoid Sage's general symbolic ring. Generic
Toeplitz diagonals are created in multivariate polynomial rings over `QQ`. Krylov rows and their
incremental fraction-free elimination stay in that exact base ring while they remain independent;
the corresponding fraction field is introduced only when the first dependence relation must be
normalized. Characteristic and annihilating polynomials are formed over those exact coefficient
parents. User-supplied
`diagonal_values` retain their exact Sage parent whenever coercion permits, including `ZZ`,
`QQ`, finite fields, polynomial rings, fraction fields, and algebraic extensions.

The nonautonomous row-column mode is exact as well. Instead of symbolic functions of an
unspecified variable, `Q[r]` is represented in a finite polynomial ring whose independent
generators stand for the shifted entries actually occurring at integer step `r`. This keeps
coercion, equality, factorization, and linear algebra in explicit algebraic parents.

## Row-column transfer examples

```python
d = lin_rec_for_dtm(2, 2, scalar_recurrence=True, verbose=False)
print(d["states"])
print(d["transfer_matrix"])
print(d["characteristic_polynomial"])

p = lin_rec_for_ptm(2, 2, scalar_recurrence=True, verbose=False)
```

For the balanced `(3,3)` case the paper predicts 20 states and 50 nonzero transitions:

```python
d33 = lin_rec_for_dtm(3, 3, scalar_recurrence=False, verbose=False)
print(d33["state_order"])
print(len(d33["transfer_matrix_expression"].dict()))
```

## Increasing rows

```python
a = lin_rec_for_dtm_increasing_rows(2, 1, ell_range="forward", verbose=False)
b = lin_rec_for_dtm_increasing_rows(2, 1, ell_range="centered", verbose=False)
print(a["annihilating_polynomial"])
```

The centered range keeps the coefficient minors smaller while producing the same
annihilating recurrence.

## Krylov scalarization

```python
dk = lin_rec_for_dtm(
    2, 2,
    scalar_recurrence=True,
    scalar_recurrence_method="krylov",
    verbose=False,
)
print(dk["observable_order"])
print(dk["scalar_recurrence_polynomial"])
```

The Krylov implementation maintains the observable row space incrementally instead of
recomputing a full matrix kernel after every new row. Pivot elimination is fraction-free until the
first dependence is found, which is particularly helpful for multivariate polynomial Toeplitz
parameters.

## Position-dependent band weights

```python
q = lin_rec_for_dtm(
    2, 2,
    position_dependent=True,
    scalar_recurrence=False,
    verbose=False,
)
print(q["transfer_matrix_expression"])
```

The displayed matrix is `Q[0]`. Its coefficients live in an exact multivariate polynomial ring with generator names such as `a_m1_at_k_m1` and `a_0_at_k_p1`. `q["transfer_matrix_function"](r)` constructs the exact matrix `Q[r]` with integer-shifted generators. A custom entry function can be supplied through `band_entry_function`; it receives an offset and an integer position as `(offset, position)` and should return elements of an exact Sage parent (or exact integers/rationals).

## Fiduccia distant-term evaluation

The recurrence convention is

```text
a[k] = c[0] a[k-1] + c[1] a[k-2] + ... + c[d-1] a[k-d].
```

Examples:

```python
fiduccia_pol_squarings([1, 23, 3], [4, 5, 6], 10**7)

fiduccia_pol_squarings(
    [1, 23, 3], [111, 3332, 12], 100,
    modulus=123, term_count=3,
)
# [36, 11, 30]
```

The function preserves the hardcoded modular-polynomial-squaring formulas described by
D. I. Khomovsky in *Efficient Computation of Terms of Linear Recurrence Sequences of Any Order*
(INTEGERS 18 (2018), A39), within the broader fast recurrence-evaluation framework
associated with Fiduccia. The mathematical transition logic is intentionally kept explicit.

For comparison with Sage's native C-finite implementation over rational/integer data:

```python
sage_cfinite_term([2, -1, 5, 4], [3, -2, 7, 1], 137)
```

Run a reproducible comparison with:

```bash
sage examples/compare_fiduccia.sage
```

## Remote Toeplitz determinants and permanents

```python
term_for_dtm(2, 2, 10**6, diagonal_values=[2, -1, 3, 4, 1])

term_for_ptm(
    1, 1, 10**8,
    diagonal_values=[2, 1, 3],
    modulus=1000003,
)
```

To request the observable recurrence when it is smaller:

```python
from sage.all import PolynomialRing, QQ
R = PolynomialRing(QQ, names=("A", "B", "C"))
A, B, C = R.gens()
term_for_dtm(
    2, 2, 10**6,
    diagonal_values=[C, B, A, B, C],
    scalar_recurrence_method="krylov",
)
```

## Reproducing paper-facing examples

Run:

```bash
sage examples/paper_examples.sage
```

The example file collects the pentadiagonal determinant/permanent transfers, the `(3,3)`
state/sparsity count, the sparse `(2,1)` characteristic recurrence and explicit formula,
the offset `+/-2` parity-split comparison, a Krylov example, a position-dependent transfer,
and a Fiduccia distant-term calculation.

## Tests

With Sage installed, run from the repository root:

```bash
sage tests/run_tests.sage
```

Backend-independent tests can also be run under CPython (with SymPy installed):

```bash
python3 -m unittest discover -s tests -p 'test_python*.py' -v
```

Static packaging and API/default checks do not require Sage:

```bash
python3 tools/verify_release.py .
```

See `VERIFICATION.md` for the validation status of this exported bundle.

## References for distant-term evaluation

- C. M. Fiduccia, "An Efficient Formula for Linear Recurrences", SIAM Journal on
  Computing 14 (1985), 106-112. DOI: 10.1137/0214007. `https://doi.org/10.1137/0214007`
- D. I. Khomovsky, "Efficient Computation of Terms of Linear Recurrence Sequences of Any Order",
  INTEGERS 18 (2018), A39.
  `https://math.colgate.edu/~integers/s39/s39.pdf`

## Citation

Please cite the accompanying paper **Constructive recurrences for determinants and
permanents of banded Toeplitz matrices**. Machine-readable citation metadata are provided
in `CITATION.cff`.

## License

This directory is intended to live inside the `matrix_codes` repository and does not choose
a separate license. The repository-level license, if present, governs this subdirectory.
