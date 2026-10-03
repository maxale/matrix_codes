# Changelog

## v14 - 2026-09-26 - Symmetric doset straightening fix

- Fixed the Mathematica coefficient-association helper used by symmetric doset straightening.
  Exponent vectors returned by `CoefficientRules` are now used directly as association keys rather
  than being wrapped around an unevaluated `Part` expression.
- Corrected the first nontrivial `m=4` doset relation to
  `Delta_{23,14} = -Delta_{12,34} + Delta_{13,24}`.
- Added a reconstruction regression that substitutes the computed straightening coefficients back
  into the corresponding generic symmetric minor.
- No generic row-column, Krylov, Fiduccia, skew two-step, or circulant-bridge production algorithm
  was changed.

## v13 - 2026-09-19 - Symmetry-paper synchronization

- Synchronized `ToeplitzSymmetryReductions.wl` with the corrected skew two-step formulation in
  *Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents*.
- Replaced the production skew route based on extracting `R` from a direct full-sequence
  polynomial `P(x)=R(x^2)` by `Q -> Q^2 -> principal Krylov R(y)`, followed by the full
  one-step annihilator `R(x^2)`.
- Reported `Binomial[2m,m]/2` explicitly as the theoretical Hodge-half dimension rather than as
  a locally constructed quotient of normalized row-column states.
- Kept the direct full-sequence Krylov calculation only in small-width regression tests as an
  independent cross-check.
- Added exact Toeplitz-circulant bridge verification files and a README reproducibility map for
  the recurrence and symmetry papers.
- The generic row-column, incremental Krylov, and Fiduccia production sources are unchanged.

## v12 - 2026-09-16 - Optional determinant symmetry reductions

- Added `ToeplitzSymmetryReductions.wl` as an optional layer accompanying
  *Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents*.
- Implemented the symmetric balanced determinant pseudocode by straightening every newly
  generated row-column minor immediately into the doset basis, assembling the Catalan
  transfer, and applying the existing optimized principal-coordinate Krylov scalarization.
- Added a skew-symmetric determinant constructor that computes the exact full observable
  polynomial and its even form `P(x)=R(x^2)`, with an autonomous companion transfer for the
  even-size subsequence and the structural half-state dimension reported separately.
- Integrated the preceding incremental fraction-free Krylov update: the generic Krylov engine
  no longer recomputes a full symbolic nullspace after every new Krylov row.
- Added symmetry examples, smoke/regression tests, API documentation, and release-verifier
  coverage. The generic row-column and Fiduccia algorithms remain unchanged.

## v10 - Public option reference and Fiduccia provenance

- Added a `Public API and options` section to `README.md` listing the supported options and
  default values for the row-column constructors, increasing-rows constructors, one-step
  Laplace helpers, `FiducciaPolSquarings`, `TermForDTM`, and `TermForPTM`.
- Documented the semantics and compatibility restrictions of `DiagonalValues`,
  `ScalarRecurrence`, `ScalarRecurrenceMethod`, `PositionDependent`, `BandEntryFunction`,
  `Modulus`, `TermCount`, and the remaining public options.
- Clarified the provenance of the hardcoded modular-polynomial squaring implementation: it
  follows D. I. Khomovsky, "Efficient Computation of Terms of Linear Recurrence Sequences of
  Any Order", INTEGERS 18 (2018), A39, in the Fiduccia-style distant-term evaluation setting.
- Added references to Khomovsky (2018) and C. M. Fiduccia, "An Efficient Formula for Linear
  Recurrences", SIAM Journal on Computing 14 (1985), 106-112, DOI 10.1137/0214007.
- Extended `tools/verify_release.py` so a release fails if the public API/default-option
  documentation or the distant-term provenance references are missing.
- Updated release metadata to v10. No recurrence construction, transfer rule, Krylov
  scalarization, Fiduccia identity, or distant-term arithmetic was changed.

## v9 - ASCII-portable documentation and current paper alignment

- Replaced the Unicode em dash and box-drawing characters in `README.md` with ASCII-only
  equivalents so the documentation renders cleanly in editors that do not detect UTF-8.
- Extended the static release verifier to require ASCII-only public documentation and metadata.
- Synchronized `CITATION.cff` with code release v9.
- Updated the paper-to-code description of the sparse `(2,1)` specialization to match the
  current paper: the code reproduces the recurrence and explicit characteristic-polynomial
  formula underlying the threefold root-of-unity geometry, radial interlacing, and
  asymptotically sharp outer radius.
- Extended `examples/PaperExamples.wl` with a finite-section check of the explicit sparse
  characteristic-polynomial formula.
- Added the paper's offset `+/-2` parity-split spectral comparison as a second finite-section
  characteristic-polynomial check.
- No transfer construction, Laplace formula, Krylov algorithm, Fiduccia squaring identity, or
  recurrence-evaluation logic was changed in this release.

## v8 - GitHub packaging and paper-title alignment

- Reorganized the Wolfram implementation into `src/`, `examples/`, `tests/`, and `tools/`.
- Replaced manuscript-version language in public documentation and source comments by the title
  **Constructive recurrences for determinants and permanents of banded Toeplitz matrices**.
- Added a paper-to-code map and explicit scope boundaries to `README.md`.
- Added `examples/PaperExamples.wl` for the pentadiagonal, balanced `(3,3)`, sparse `(2,1)`,
  position-dependent, Krylov, and distant-term examples.
- Added `CITATION.cff`, `PACKAGE_INFO.txt`, `.gitignore`, a release manifest, and a standalone
  static release verifier.
- Updated test loaders so all tests run correctly from the new directory layout.
- No transfer construction, Laplace formula, Krylov algorithm, Fiduccia squaring identity, or
  recurrence-evaluation logic was changed in this release.

## Earlier development fixes retained

- Default symbolic output was isolated from user assignments to ordinary `a`, `x`, and `k` by
  using Wolfram formal symbols.
- Public Fiduccia options and messages were initialized in a safe order so that `ClearAll` cannot
  erase their definitions after installation.
- Exact and modular distant-term routines were consolidated into `FiducciaPolSquarings`, with
  `TermForDTM` and `TermForPTM` as Toeplitz-facing wrappers.
- Optional observable/Krylov scalarization was added without changing the default
  characteristic-polynomial route.
