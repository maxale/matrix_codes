# Toeplitz recurrence algorithms - Wolfram Language

This directory contains the Mathematica / Wolfram Language implementation accompanying
**Constructive recurrences for determinants and permanents of banded Toeplitz matrices**
by Max A. Alekseyev and Dmitry I. Khomovsky.

The code constructs finite recurrences for determinants and permanents of fixed-band
Toeplitz matrices, supports the two Laplace constructions developed in the paper, adds
optional Krylov scalarization, handles the position-dependent row-column cocycle, and
provides Fiduccia-style evaluation of distant terms once a constant recurrence is known.
No external Wolfram packages are required.

Public repository: `https://github.com/maxale/matrix_codes`

## Requirements

- Mathematica / Wolfram Language 2021.1 or later.
- No external packages.
- Python 3 is optional and is used only by `tools/verify_release.py` for static release checks.

## Repository layout

```text
toeplitz_recurrences_mathematica/
|-- README.md
|-- CITATION.cff
|-- CHANGELOG.md
|-- PACKAGE_INFO.txt
|-- VERIFICATION.md
|-- MANIFEST.sha256
|-- src/
|   |-- ToeplitzRecurrences.wl
|   |-- ToeplitzRecurrencesKrylov.wl
|   |-- ToeplitzRecurrencesFiduccia.wl
|   `-- ToeplitzRecurrencesAll.wl
|-- examples/
|   `-- PaperExamples.wl
|-- tests/
|   |-- RunToeplitzAllTests.wl
|   |-- ToeplitzRecurrencesExamples.wl
|   |-- ToeplitzRecurrencesKrylovExamples.wl
|   |-- ToeplitzRecurrencesFiducciaExamples.wl
|   `-- smoke/
|       |-- DispatchSmokeTest.wl
|       |-- KrylovDispatchSmokeTest.wl
|       |-- FiducciaSmokeTest.wl
|       `-- ScopeSmokeTest.wl
`-- tools/
    `-- verify_release.py
```

## Quick start

From the root of this directory:

```wl
Get["src/ToeplitzRecurrencesAll.wl"];
```

The base layer alone can be loaded with

```wl
Get["src/ToeplitzRecurrences.wl"];
```

The public convention is always

```text
m1 = number of subdiagonals    (offsets -m1,...,-1)
m2 = number of superdiagonals (offsets 1,...,m2)
```

and `DiagonalValues` are supplied in offset order `-m1,...,m2`.

## Paper-to-code map

| Paper construction or application | Main Wolfram interface | Notes |
|---|---|---|
| Increasing-rows construction | `LinRecForDTMIncreasingRows`, `LinRecForPTMIncreasingRows` | Supports `EllRange -> "Forward"` and `EllRange -> "Centered"`; transparently transposes when `m1 < m2`. |
| Row-column transfer construction | `LinRecForDTM`, `LinRecForPTM` | Uses the same first-discovery state ordering as the paper. |
| Characteristic-polynomial scalarization | `ScalarRecurrence -> True` | Default scalarization when the transfer is constant. |
| Observable/Krylov scalarization | `ScalarRecurrenceMethod -> "Krylov"` | Can produce an annihilator smaller than the full transfer characteristic polynomial. |
| Position-dependent row-column cocycle | `PositionDependent -> True` | Returns a position-dependent transfer; constant scalarization is intentionally disabled. |
| Distant-term evaluation | `FiducciaPolSquarings`, `TermForDTM`, `TermForPTM` | Applies only after a homogeneous constant-coefficient recurrence has been obtained. |
| Sparse `(2,1)` Toeplitz-Hessenberg spectral example | `LinRecForDTM[2,1,...]`, `TermForDTM[2,1,...]` | `examples/PaperExamples.wl` constructs `det(lambda I-A_n)`, checks the explicit coefficient formula, and exposes the recurrence behind the threefold spectrum. |
| Offset `+/-2` parity-split comparison | `TermForDTM[2,2,...]` | Reorders odd/even indices into two tridiagonal blocks and checks the product of the corresponding Lucas-U characteristic polynomials. |

The code does **not** currently provide dedicated implementations of the cyclic-closure
results, the proposed finite-corner-defect state extension, or the nonautonomous
increasing-rows theorem. The compound/Widom similarity theorem is likewise a mathematical
identification in the paper rather than a separate runtime constructor in this package.

## Public API and options

The tables below assume the recommended full loader
`Get["src/ToeplitzRecurrencesAll.wl"]`. In particular, the full loader installs the optional
`ScalarRecurrenceMethod` selector on `LinRecForDTM` and `LinRecForPTM`.

### Row-column recurrence constructors

`LinRecForDTM[m1,m2,...]` and `LinRecForPTM[m1,m2,...]` accept the following options.

| Option | Default | Meaning |
|---|---|---|
| `DiagonalValues` | `Automatic` | Constant diagonal values in offset order `-m1,...,m2`. `Automatic` uses `\[FormalA][-m1],...,\[FormalA][m2]`. |
| `ScalarRecurrence` | `Automatic` | `True` requests scalarization, `False` skips it, and `Automatic` scalarizes only when the transfer order does not exceed `MaxCharacteristicOrder`. |
| `CharacteristicPolynomialVariable` | `\[FormalX]` | Polynomial variable used for the returned scalar recurrence polynomial. |
| `MaxCharacteristicOrder` | `20` | Transfer-order cutoff used only when `ScalarRecurrence -> Automatic`. |
| `Verbose` | `True` | Controls interactive printing of the transfer order, matrix, and scalar polynomial. |
| `PositionDependent` | `False` | `True` constructs the nonautonomous row-column cocycle instead of a constant transfer. |
| `BandEntryFunction` | `Automatic` | Used only with `PositionDependent -> True`. A custom value has the form `Function[{s,t}, ...]`; `Automatic` uses formal entries `\[FormalA][s][t]`. |
| `ScalarRecurrenceMethod` | `"CharacteristicPolynomial"` | With the Krylov extension loaded, choose `"CharacteristicPolynomial"` or `"Krylov"`. |

`PositionDependent -> True` cannot be combined with explicit `DiagonalValues`. Constant
scalarization is disabled in position-dependent mode.

### Increasing-rows constructors

`LinRecForDTMIncreasingRows[m1,m2,...]` and
`LinRecForPTMIncreasingRows[m1,m2,...]` accept:

| Option | Default | Meaning |
|---|---|---|
| `DiagonalValues` | `Automatic` | Constant diagonal values in offset order `-m1,...,m2`. |
| `EllRange` | `"Forward"` | Use `"Forward"`, `"Centered"`, or an explicit strictly increasing list of the required `d+1` levels. |
| `CharacteristicPolynomialVariable` | `\[FormalX]` | Variable used for the returned recurrence polynomial. |
| `Verbose` | `True` | Controls interactive printing. |

### One-step Laplace helpers

`LaplaceD[m1,m2,S,...]` and `LaplaceP[m1,m2,S,...]` accept:

| Option | Default | Meaning |
|---|---|---|
| `DiagonalValues` | `Automatic` | Constant diagonal values in offset order `-m1,...,m2`. |
| `PositionDependent` | `False` | Select constant or position-dependent band entries. |
| `BandEntryFunction` | `Automatic` | Custom `Function[{s,t}, ...]` in position-dependent mode. |

### Distant-term evaluation

`FiducciaPolSquarings[initial,recurrence,N,...]` accepts:

| Option | Default | Meaning |
|---|---|---|
| `Modulus` | `None` | Exact/symbolic arithmetic by default; an integer greater than 1 selects arithmetic modulo that integer. |
| `TermCount` | `1` | Return one term by default; `TermCount -> r` returns `r` consecutive terms starting at the zero-based index `N`. |

`TermForDTM[m1,m2,N,...]` and `TermForPTM[m1,m2,N,...]` accept:

| Option | Default | Meaning |
|---|---|---|
| `DiagonalValues` | `Automatic` | Constant Toeplitz diagonal values in offset order `-m1,...,m2`. |
| `ScalarRecurrenceMethod` | `"CharacteristicPolynomial"` | Choose the full transfer characteristic polynomial or, when loaded, `"Krylov"`. |
| `Modulus` | `None` | Exact/symbolic arithmetic or integer arithmetic modulo the supplied integer. |
| `TermCount` | `1` | One distant term or a consecutive window of terms. |
| `Verbose` | `False` | Controls printing by the internal recurrence constructor. |
| `PositionDependent` | `False` | Must remain `False`; Fiduccia evaluation requires a homogeneous constant-coefficient recurrence. |

The current option defaults can also be queried directly in Mathematica, for example:

```wl
Options[LinRecForDTM]
Options[LinRecForDTMIncreasingRows]
Options[FiducciaPolSquarings]
Options[TermForDTM]
```

For the sparse `(2,1)` double-band specialization, the paper uses the recurrence to recover
the classical threefold root-of-unity geometry and then refines it with explicit radial roots,
strict radial interlacing, and an asymptotically sharp outer radius. The package does not
reimplement the general zero-location theorems; instead, `PaperExamples.wl` reproduces the
characteristic-polynomial recurrence and checks its explicit finite-section formula.
The same example file also records the paper's offset `+/-2` parity-split comparison: after
permuting odd and even coordinates, the matrix becomes a direct sum of two ordinary
tridiagonal Toeplitz blocks, and the code checks the resulting product formula.

## Row-column transfer examples

A symbolic pentadiagonal determinant transfer and its scalar recurrence:

```wl
r = LinRecForDTM[2, 2,
  ScalarRecurrence -> True,
  Verbose -> False
];

Normal[r["TransferMatrix"]]
r["CharacteristicPolynomial"]
```

The permanent uses the same state graph with signless Laplace weights:

```wl
rp = LinRecForPTM[2, 2,
  ScalarRecurrence -> True,
  Verbose -> False
];
```

For the balanced `(3,3)` case, the paper predicts 20 states and 50 nonzero
transitions. The following constructs that transfer without requesting a degree-20
characteristic polynomial:

```wl
r33 = LinRecForDTM[3, 3,
  ScalarRecurrence -> False,
  Verbose -> False
];

Length[r33["States"]]
Count[Flatten[Normal[r33["TransferMatrix"]]], z_ /; ! TrueQ[z === 0]]
```

## Increasing rows

```wl
LinRecForDTMIncreasingRows[2, 1,
  EllRange -> "Forward",
  Verbose -> False
]

LinRecForDTMIncreasingRows[2, 1,
  EllRange -> "Centered",
  Verbose -> False
]
```

The centered range keeps the coefficient minors smaller while producing the same
annihilating recurrence.

## Krylov scalarization

The default row-column scalarization is the full characteristic polynomial:

```wl
LinRecForDTM[2, 2,
  ScalarRecurrence -> True,
  ScalarRecurrenceMethod -> "CharacteristicPolynomial",
  Verbose -> False
]
```

The observable route is selected by

```wl
LinRecForDTM[2, 2,
  ScalarRecurrence -> True,
  ScalarRecurrenceMethod -> "Krylov",
  Verbose -> False
]
```

The returned association includes `"ObservableOrder"` and
`"ScalarRecurrencePolynomial"` when Krylov scalarization is computed.

## Position-dependent band weights

```wl
q = LinRecForDTM[2, 2,
  PositionDependent -> True,
  ScalarRecurrence -> False,
  Verbose -> False
];

q["TransferMatrixExpression"]
q["TransferMatrixFunction"]
```

A custom band entry can be supplied with

```wl
BandEntryFunction -> Function[{s, t}, b[s, t]]
```

In this mode the transfer is a cocycle. `TermForDTM` and `TermForPTM` therefore reject
`PositionDependent -> True`, because Fiduccia evaluation requires a constant-coefficient
recurrence.

## Fiduccia distant-term evaluation

The public name `FiducciaPolSquarings` reflects the Fiduccia-style fast evaluation setting.
The hardcoded polynomial-squaring transitions implemented here follow Dmitry I. Khomovsky,
*Efficient Computation of Terms of Linear Recurrence Sequences of Any Order*, INTEGERS 18
(2018), A39. That paper explicitly develops modular-polynomial-squaring formulas as an
approach distinct from a generic polynomial multiplication/reduction implementation of
Fiduccia's 1985 algorithm; the present code preserves those hardcoded identities.

For a recurrence

```text
a_n = p1 a_(n-1) + ... + pd a_(n-d),
```

use

```wl
FiducciaPolSquarings[
  {1, 23, 3}, {111, 3332, 12}, 100,
  Modulus -> 123,
  TermCount -> 3
]
```

which returns `{36, 11, 30}`.

For Toeplitz determinants and permanents the wrappers first construct a scalar recurrence
and then invoke the same evaluator:

```wl
TermForDTM[2, 2, 10^6,
  DiagonalValues -> {2, -1, 3, 4, 1}]

TermForPTM[1, 1, 10^8,
  DiagonalValues -> {2, 1, 3},
  Modulus -> 1000003]
```

The Krylov recurrence may be requested when it is smaller:

```wl
TermForDTM[2, 2, 10^3,
  DiagonalValues -> {c, B, A, B, c},
  ScalarRecurrenceMethod -> "Krylov"]
```

## Reproducing paper-facing examples

```wl
Get["examples/PaperExamples.wl"];
RunToeplitzPaperExamples[]
```

This collects the pentadiagonal determinant/permanent transfers, the `(3,3)` state and
sparsity counts, the sparse `(2,1)` characteristic recurrence and explicit formula, the offset `+/-2`
parity-split comparison, a position-dependent transfer example, and a distant-term calculation
in one association.

## Tests

Run the full Mathematica regression suite from the repository root:

```wl
Get["tests/RunToeplitzAllTests.wl"];
RunToeplitzAllTests[]
```

Focused smoke tests are in `tests/smoke/`. In particular, the namespace-isolation
regression deliberately assigns values to ordinary `a`, `x`, and `k` and verifies that
the package's protected formal symbols remain unaffected:

```wl
Get["tests/smoke/ScopeSmokeTest.wl"];
RunToeplitzScopeSmokeTest[]
```

Static packaging and provenance checks, which do not require Mathematica, can be run with

```bash
python3 tools/verify_release.py .
```

See `VERIFICATION.md` for the validation status of this exported bundle.

## Symbol hygiene

Default symbolic parameters use Wolfram formal symbols `\[FormalA]`, `\[FormalX]`, and
`\[FormalK]`; internal computational variables are localized. Thus prior user assignments
such as `a = 9` do not alter the symbolic transfer matrices. Explicit symbols or values can
still be supplied through `DiagonalValues`, `BandEntryFunction`, and
`CharacteristicPolynomialVariable`.

## References for distant-term evaluation

- C. M. Fiduccia, "An Efficient Formula for Linear Recurrences", SIAM Journal on
  Computing 14 (1985), 106-112. DOI: 10.1137/0214007. `https://doi.org/10.1137/0214007`
- D. I. Khomovsky, "Efficient Computation of Terms of Linear Recurrence Sequences of Any
  Order", INTEGERS 18 (2018), A39.
  `https://math.colgate.edu/~integers/s39/s39.pdf`

## Citation

Please cite the accompanying paper **Constructive recurrences for determinants and
permanents of banded Toeplitz matrices**. Machine-readable citation metadata are provided
in `CITATION.cff`.

## License

This directory is intended to live inside the `matrix_codes` repository and does not choose
a separate license. The repository-level license, if present, governs this subdirectory.
