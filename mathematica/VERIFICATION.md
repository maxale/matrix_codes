# Verification status

This export is the GitHub-facing Wolfram Language implementation accompanying
**Constructive recurrences for determinants and permanents of banded Toeplitz matrices**.

## What was changed for this export

The production algorithms are carried forward unchanged. This release adds a complete README
reference for the public options and their default values, documents the Khomovsky (2018)
provenance of the hardcoded modular-polynomial squaring evaluator, and strengthens the release
verifier so those API/default and provenance statements cannot silently disappear. The
repository layout remains stable.

## Static verification

Run from the repository directory:

```bash
python3 tools/verify_release.py .
```

The verifier checks:

- required GitHub-facing files and directory layout;
- ASCII-only public documentation and metadata;
- use of the paper title rather than manuscript filenames or manuscript revision labels;
- absence of legacy `Global` symbolic placeholders in production source;
- presence of the public interfaces corresponding to increasing rows, row-column transfer,
  Krylov scalarization, position-dependent transfer, and Fiduccia evaluation;
- README documentation of the public options and their default values for all option-bearing
  public procedures;
- README provenance for the hardcoded modular-polynomial squaring implementation, including
  Khomovsky (2018) and the original Fiduccia (1985) reference;
- paper-facing examples for `(2,2)`, `(3,3)`, sparse `(2,1)` (including its explicit
  characteristic-polynomial coefficient formula), the offset `+/-2` parity-split product formula,
  position-dependent transfer, and distant-term evaluation;
- Wolfram delimiter/comment/string balance;
- SHA-256 manifest integrity and coverage.

## Mathematica runtime verification

A Wolfram kernel is required for the actual `TestReport` suites. From the repository root,
run in a fresh kernel:

```wl
Get["tests/smoke/KrylovDispatchSmokeTest.wl"]
Get["tests/smoke/DispatchSmokeTest.wl"]
Get["tests/smoke/FiducciaSmokeTest.wl"]
Get["tests/smoke/ScopeSmokeTest.wl"];
RunToeplitzScopeSmokeTest[]
```

Then run the complete regression suite:

```wl
Get["tests/RunToeplitzAllTests.wl"];
RunToeplitzAllTests[]
```

Finally, the publication-facing examples can be evaluated with

```wl
Get["examples/PaperExamples.wl"];
RunToeplitzPaperExamples[]
```

The expected `"Balanced33"` summary is a state count of `20` and a nonzero transition count
of `50`. The expected sparse `(2,1)` scalar recurrence polynomial for
`det(lambda I-A_n)` is

```text
X^3 - lambda X^2 + x^2 y.
```

## Environment note

The export environment used to build this archive does not contain Mathematica,
`WolframKernel`, or `wolframscript`. Consequently the Wolfram runtime suite cannot be
executed here; successful static verification is not presented as a substitute for the
fresh-kernel `TestReport` run above.

## Export-time results

The final static verifier reports `PASS`: 21 required paths, 13 Wolfram files, public
API/default-option documentation, distant-term provenance references, ASCII portability,
source structure, and full manifest integrity were checked.

The packaged regression sources contain 43 base tests, 20 Krylov tests, 30 Fiduccia tests,
and 8 namespace-isolation tests. These Wolfram tests are included but were not executed in
the export environment because no Wolfram runtime is installed.

Relative to the preceding GitHub export, `ToeplitzRecurrences.wl`,
`ToeplitzRecurrencesKrylov.wl`, and `ToeplitzRecurrencesAll.wl` are byte-for-byte unchanged.
`ToeplitzRecurrencesFiduccia.wl` differs only in its leading documentation comment; after
Wolfram comments are stripped, the source is identical. Thus no computational definition
changed in this release.
