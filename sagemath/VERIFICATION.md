# Verification

This release keeps the tested generic recurrence core unchanged and synchronizes the optional
symmetry-aware determinant module with **Symmetry reductions and recurrence degrees for banded
Toeplitz determinants and permanents**.

## Symmetry-module checks

The symmetric constructor continues to follow the paper-facing state-generation pseudocode:
newly generated balanced row-column minors are decoded as subset pairs, immediately
straightened into the doset basis, and only then entered into the reduced transfer.  The
backend-independent tests check Catalan dimensions `2, 5, 14, 42` for semibandwidths
`1, 2, 3, 4`, a nontrivial `m=4` straightening relation, and round-trip conversion between
subset-pair and normalized signature coordinates.

The skew constructor now follows the corrected two-step route.  It keeps the normalized
row-column transfer `Q`, forms `Q^2`, applies principal Krylov scalarization directly to that
two-step transfer, and lifts the resulting polynomial `R(y)` to the full one-step annihilator
`R(x^2)`.  The value `binomial(2*m,m)/2` is reported as the theoretical Hodge-half dimension;
`hodge_half_constructed` is false because no local quotient of the normalized row-column
states is claimed.  A backend-independent contract test verifies that the production helper
squares the transfer before invoking Krylov.

The Sage runtime symmetry suite checks skew semibandwidths `m=1,2,3` for two-step/full orders
`1/2`, `3/6`, and `9/18`.  It also independently recomputes the direct full-sequence Krylov
polynomial and checks agreement with the lifted `R(x^2)`.

## Toeplitz-circulant bridge check

`tests/test_sage_circulant_bridge.sage` verifies exact rational examples of the fixed-size
Jacobi correction

```text
det(T_n) = det(C_(n+m)) * det(C_(n+m)^(-1)[S,S]).
```

The companion `examples/circulant_bridge_examples.sage` exposes one exact paper-facing
instance.  These files are verification artifacts, not a general cyclic-closure API.

## Verified in the export environment

- CPython backend-independent regression suite passes when run from the release tree.
- Generic `src/toeplitz_recurrences.py` is byte-for-byte unchanged from the tested v5 baseline.
- The exact-ring policy remains in force: generic constant parameters use exact polynomial
  rings, Krylov division is delayed to the coefficient fraction field, and the public
  code/test/example surfaces contain no Sage Symbolic Ring (`SR`) dependency.
- Static release verification checks layout, API/default documentation, symmetry-module
  documentation, paper-facing examples, citation metadata, ASCII public text, and manifest
  integrity.

## SageMath runtime check still required

The export environment does not provide a `sage` executable, so the Sage-specific runtime
suite cannot be executed here.  On a SageMath installation, run from the repository root:

```bash
sage tests/run_tests.sage
```

This runs the generic row-column, increasing-rows, Krylov, Fiduccia, symmetry, and
Toeplitz-circulant bridge suites, then evaluates all paper-facing examples.
