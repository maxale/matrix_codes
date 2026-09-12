# Changelog

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
