# Verification status

This release contains the Wolfram Language implementation accompanying
**Constructive recurrences for determinants and permanents of banded Toeplitz matrices**,
with an optional determinant-symmetry module accompanying
**Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents**.

The generic row-column, optimized incremental Krylov, and Fiduccia production sources remain
unchanged in this release.

The v14 symmetry fix uses the exponent vectors produced by `CoefficientRules` directly as
association keys in the doset-straightening coefficient map.  The runtime suite retains the
first nontrivial `m=4` relation and additionally checks that the returned coefficients reconstruct
the original generic symmetric minor exactly.

## Symmetry-module checks

For symmetric balanced bands, `SymmetricDeterminantRecurrence` retains the direct doset/Catalan
construction: every generated boundary minor is normalized, decoded, immediately straightened,
and entered only in standard doset coordinates.

For skew-symmetric bands, `SkewSymmetricDeterminantRecurrence` now keeps the normalized
row-column transfer `Q`, forms `Q^2`, and applies principal Krylov scalarization directly to the
two-step transfer.  The returned two-step polynomial `R(y)` is lifted to the full one-step
annihilator `R(x^2)`.  `StructuralHodgeHalfOrder` records the theoretical
`Binomial[2m,m]/2` Hodge-half dimension, while `HodgeHalfConstructed -> False` makes explicit
that no local quotient of the normalized row-column states is being claimed.

The runtime symmetry tests cover semibandwidths `m=1,2,3`, checking two-step/full orders
`1/2`, `3/6`, and `9/18`, checking that the stored two-step matrix is exactly `Q.Q`, and
independently comparing the lifted annihilator with a direct full-sequence Krylov computation.

## Toeplitz-circulant bridge check

`tests/ToeplitzCirculantBridgeExamples.wl` verifies exact integer/rational instances of the
fixed-size Jacobi correction

```text
det(T_n) = det(C_(n+m)) det(C_(n+m)^(-1)[S,S]).
```

`examples/CirculantBridgeExamples.wl` provides a paper-facing exact example.  These are
verification files, not a new cyclic-closure constructor.

## Static verification

Run from the repository directory:

```bash
python3 tools/verify_release.py .
```

The verifier checks required files, ASCII public documentation, paper-title alignment,
public API/default documentation, symmetry and bridge examples/tests, Wolfram delimiter and
comment balance, and SHA-256 manifest integrity.

## Mathematica runtime verification

A Wolfram kernel is required for the actual `TestReport` suites.  From the repository root:

```wl
Get["tests/RunToeplitzAllTests.wl"];
RunToeplitzAllTests[]
```

Then evaluate the publication-facing examples:

```wl
Get["examples/PaperExamples.wl"];
RunToeplitzPaperExamples[]
Get["examples/SymmetryExamples.wl"];
RunToeplitzSymmetryExamples[]
Get["examples/CirculantBridgeExamples.wl"];
RunToeplitzCirculantBridgeExample[]
```

## Environment note

The export environment used to build this archive does not contain Mathematica,
`WolframKernel`, or `wolframscript`.  Static verification is therefore not presented as a
substitute for a fresh-kernel `TestReport` run.
