# Changelog

## v7 - 2026-09-19

- Synchronized the optional symmetry module with the corrected skew two-step formulation in
  *Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents*.
- Replaced the production skew route based on extracting `R` from a direct full-sequence
  polynomial `P(x)=R(x^2)` by the constructive route `Q -> Q^2 -> principal Krylov R(y)`,
  followed by the full one-step annihilator `R(x^2)`.
- Reported `binomial(2*m,m)/2` explicitly as the theoretical Hodge-half dimension rather than
  as a locally constructed quotient of normalized row-column states.
- Retained the direct full-sequence Krylov computation only in small-width runtime tests as an
  independent cross-check of the lifted two-step annihilator.
- Added exact Toeplitz-circulant bridge verification scripts and a README reproducibility map
  separating the companion recurrence-paper code from the symmetry-paper code.
- The tested v5 generic recurrence core and the v6 symmetric doset/Catalan construction are
  unchanged.

## v6 - 2026-09-16

- Added `src/toeplitz_symmetry_reductions.py` as an optional layer accompanying
  *Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents*.
- Implemented the symmetric balanced determinant pseudocode by straightening every newly
  generated row-column minor immediately into the doset basis, assembling the Catalan
  transfer, and applying the existing optimized principal-coordinate Krylov scalarization.
- Added a skew-symmetric determinant constructor that computes the exact full observable
  polynomial and its even form `P(x)=R(x^2)`, with an autonomous companion transfer for the
  even-size subsequence and the structural half-state dimension reported separately.
- Added backend-independent doset/straightening tests plus Sage runtime comparisons against
  the unrestricted principal Krylov polynomial, symmetry examples, and release-verifier coverage.
- The v5 generic recurrence core remains algorithmically unchanged and continues to avoid
  Sage's Symbolic Ring.

## v5 - 2026-09-13

- Fixed the internal fresh polynomial-variable name used by `term_for_dtm` and `term_for_ptm`: temporary names now begin with a letter (`toeplitz_lambda_N`), as required by SageMath `PolynomialRing`.
- Added a backend-independent contract regression for Sage-compatible temporary variable names and a Sage runtime construction check in the Fiduccia integration suite.
- Retained the v4 incremental fraction-free Krylov optimization unchanged.
- Mathematical recurrence formulas and public API semantics are unchanged.

## v4 - 2026-09-12

- Replaced repeated full right-kernel recomputation in Krylov scalarization with an incremental exact row-echelon reducer.
- The reducer uses fraction-free pivot elimination while Krylov rows are independent and introduces the coefficient fraction field only when the first dependence relation is extracted.
- Added backend-independent dependence-certificate tests, including randomized comparison against SymPy nullspaces.
- Reduced the Sage Krylov regression workload without weakening coverage: the generic five-parameter and symmetric three-parameter examples now use separate polynomial rings, the duplicated characteristic-polynomial computation was removed, and the observable determinant is checked back in its polynomial ring.
- Public Krylov outputs and recurrence formulas are unchanged.

## v3 - 2026-09-12

- Eliminated the Sage Symbolic Ring from production code, tests, and paper-facing examples.
- Generic constant Toeplitz coefficients now live in explicit multivariate polynomial rings
  over `QQ`; Krylov elimination uses the corresponding fraction field only when needed.
- Reworked position-dependent row-column transfers to use exact finite polynomial rings for
  the shifted entries occurring in each integer-step transfer `Q[r]`. Custom
  `band_entry_function` callbacks now receive exact integer positions.
- Removed the symbolic fallback from the increasing-rows annihilator construction: an
  unsupported coefficient parent now raises a targeted error instead of silently switching
  algebra systems.
- Rewrote Sage regression tests and paper-facing examples to construct parameters with
  `PolynomialRing` and to compare objects in compatible exact parents.
- Added a release regression that rejects Sage Symbolic Ring dependencies.
- Mathematical recurrence, state-closure, Krylov, and Fiduccia formulas are otherwise
  unchanged.

## v2 - 2026-09-12

- Reorganized the SageMath supplement into a GitHub-ready stable layout with `src/`,
  `examples/`, `tests/`, and `tools/` directories.
- Aligned public documentation with *Constructive recurrences for determinants and
  permanents of banded Toeplitz matrices*.
- Added complete public API/default documentation for row-column, increasing-rows,
  Laplace, Krylov, Fiduccia, remote-term, and Sage C-finite interfaces.
- Added paper-facing examples for the `(2,2)` transfers, balanced `(3,3)` state/sparsity
  count, sparse `(2,1)` threefold spectral recurrence, offset `+/-2` parity split,
  Krylov scalarization, position-dependent transfer, and distant-term evaluation.
- Added `CITATION.cff`, release metadata, an ASCII-safe README, and static release
  verification including API/default checks and manifest integrity.
- Clarified that the hardcoded modular-polynomial squaring formulas in
  `fiduccia_pol_squarings` follow Khomovsky (INTEGERS 18 (2018), A39) within the
  broader Fiduccia fast recurrence-evaluation framework.
- The recurrence, state-closure, Krylov, increasing-rows, and Fiduccia mathematical
  algorithms are unchanged from the preceding Sage implementation.
