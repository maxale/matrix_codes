(* ::Package:: *)

(* Reproducible examples accompanying
   "Constructive recurrences for determinants and permanents of banded Toeplitz matrices".

   Load this file from any working directory.  RunToeplitzPaperExamples[] returns
   a single Association collecting the examples; no Global symbols are required.
*)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "src", "ToeplitzRecurrencesAll.wl"}]];

ClearAll[RunToeplitzPaperExamples];

RunToeplitzPaperExamples[] := Module[
  {fx = \[FormalX], fy = \[FormalY], fl = \[FormalL],
   det22, per22, det33, sparse21, sparseTerms, sparseFormula,
   parityTerms, parityFormula, krylov22, position22, nonzero33},

  det22 = LinRecForDTM[2, 2,
    ScalarRecurrence -> True,
    Verbose -> False
  ];

  per22 = LinRecForPTM[2, 2,
    ScalarRecurrence -> True,
    Verbose -> False
  ];

  det33 = LinRecForDTM[3, 3,
    ScalarRecurrence -> False,
    Verbose -> False
  ];
  nonzero33 = Count[
    Flatten[Normal[det33["TransferMatrix"]]],
    z_ /; !TrueQ[z === 0]
  ];

  (* For A_n with only a_-2 = y and a_1 = x nonzero, construct the
     determinant recurrence of lambda I - A_n.  The expected scalar
     recurrence polynomial is X^3 - lambda X^2 + x^2 y. *)
  sparse21 = LinRecForDTM[2, 1,
    DiagonalValues -> {-fy, 0, fl, -fx},
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    CharacteristicPolynomialVariable -> \[FormalX],
    Verbose -> False
  ];

  (* Finite-section characteristic polynomials for the same sparse family.
     This checks chi_n(lambda) = Sum_j (-1)^j Binomial[n-2j,j]
     (x^2 y)^j lambda^(n-3j), the formula used in the spectral corollary. *)
  sparseTerms = Table[
    Expand[TermForDTM[2, 1, n,
      DiagonalValues -> {-fy, 0, fl, -fx},
      ScalarRecurrenceMethod -> "Krylov"
    ]],
    {n, 0, 12}
  ];
  sparseFormula = Table[
    Expand[Sum[
      (-1)^j Binomial[n - 2 j, j] (fx^2 fy)^j fl^(n - 3 j),
      {j, 0, Floor[n/3]}
    ]],
    {n, 0, 12}
  ];

  (* Classical offset +/-2 comparison from the paper.  Splitting odd and
     even indices gives tridiagonal blocks of sizes Ceiling[n/2] and
     Floor[n/2], so the characteristic polynomial is the product of two
     Lucas-U polynomials with parameters (lambda, x y). *)
  parityTerms = Table[
    Expand[TermForDTM[2, 2, n,
      DiagonalValues -> {-fy, 0, fl, 0, -fx},
      ScalarRecurrenceMethod -> "Krylov"
    ]],
    {n, 0, 12}
  ];
  parityFormula = Table[
    With[{p = Ceiling[n / 2], q = Floor[n / 2]},
      Expand[
        Sum[
          (-1)^j Binomial[p - j, j] (fx fy)^j fl^(p - 2 j),
          {j, 0, Floor[p / 2]}
        ]
        Sum[
          (-1)^j Binomial[q - j, j] (fx fy)^j fl^(q - 2 j),
          {j, 0, Floor[q / 2]}
        ]
      ]
    ],
    {n, 0, 12}
  ];

  krylov22 = LinRecForDTM[2, 2,
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    Verbose -> False
  ];

  position22 = LinRecForDTM[2, 2,
    PositionDependent -> True,
    ScalarRecurrence -> False,
    Verbose -> False
  ];

  <|
    "PentadiagonalDeterminant" -> <|
      "TransferMatrix" -> Normal[det22["TransferMatrix"]],
      "CharacteristicPolynomial" -> det22["CharacteristicPolynomial"]
    |>,
    "PentadiagonalPermanent" -> <|
      "TransferMatrix" -> Normal[per22["TransferMatrix"]],
      "CharacteristicPolynomial" -> per22["CharacteristicPolynomial"]
    |>,
    "Balanced33" -> <|
      "StateCount" -> Length[det33["States"]],
      "NonzeroTransitionCount" -> nonzero33
    |>,
    "Sparse21Spectrum" -> <|
      "ScalarRecurrencePolynomial" -> sparse21["ScalarRecurrencePolynomial"],
      "ExpectedPolynomial" -> \[FormalX]^3 - fl \[FormalX]^2 + fx^2 fy,
      "Terms0Through12" -> sparseTerms,
      "ExplicitFormula0Through12" -> sparseFormula,
      "ExplicitFormulaCheck" -> And @@ MapThread[TrueQ[#1 === #2] &, {sparseTerms, sparseFormula}]
    |>,
    "OffsetPlusMinus2ParitySpectrum" -> <|
      "Terms0Through12" -> parityTerms,
      "ParityFormula0Through12" -> parityFormula,
      "ParitySplitCheck" -> And @@ MapThread[TrueQ[#1 === #2] &, {parityTerms, parityFormula}]
    |>,
    "Krylov22" -> <|
      "ObservableOrder" -> krylov22["ObservableOrder"],
      "ScalarRecurrencePolynomial" -> krylov22["ScalarRecurrencePolynomial"]
    |>,
    "PositionDependent22" -> <|
      "States" -> position22["States"],
      "TransferMatrixExpression" -> position22["TransferMatrixExpression"]
    |>,
    "FiducciaExample" -> FiducciaPolSquarings[
      {1, 23, 3}, {111, 3332, 12}, 100,
      Modulus -> 123, TermCount -> 3
    ]
  |>
];
