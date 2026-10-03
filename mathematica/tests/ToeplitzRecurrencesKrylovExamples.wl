(* ::Package:: *)

(* Regression checks for the optional observable/Krylov scalarization.
   The suite includes an explicit dispatch regression for a previously
   observed Mathematica 14.2.1 dispatch failure. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "src", "ToeplitzRecurrences.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "src", "ToeplitzRecurrencesKrylov.wl"}]];

ClearAll[TRKrylovSamePolynomialQ, RunToeplitzKrylovTests];
TRKrylovSamePolynomialQ[p_, q_] := TrueQ[Expand[p - q] === 0];

RunToeplitzKrylovTests[] := Module[
  {x = \[FormalX], aa = \[FormalA],
   A, B, C, k11, k22, kp22, ks22, c22, skipped, pos22, manualKrylov, scopeKrylov, scopeGlobals, tests,
   expected22, expectedSym5},

  k11 = LinRecForDTM[1, 1,
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    Verbose -> False];

  k22 = LinRecForDTM[2, 2,
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    Verbose -> False];

  c22 = LinRecForDTM[2, 2,
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "CharacteristicPolynomial",
    Verbose -> False];

  kp22 = LinRecForPTM[2, 2,
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    Verbose -> False];

  ks22 = LinRecForDTM[2, 2,
    DiagonalValues -> {C, B, A, B, C},
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    Verbose -> False];

  skipped = LinRecForDTM[2, 2,
    ScalarRecurrence -> False,
    ScalarRecurrenceMethod -> "Krylov",
    Verbose -> False];

  pos22 = LinRecForDTM[2, 2,
    PositionDependent -> True,
    ScalarRecurrence -> False,
    ScalarRecurrenceMethod -> "Krylov",
    Verbose -> False];

  (* Direct private-level regression with a non-first observable coordinate.
     This exercises the pivot ordering used by the incremental elimination. *)
  manualKrylov = ToeplitzRecurrences`Private`trKrylovScalarization[
    {{1, 2, 0}, {0, 3, 4}, {5, 0, 6}}, x, 2
  ];

  {scopeKrylov, scopeGlobals} = Block[
    {Global`a = 9, Global`x = 17, Global`k = 23},
    {
      LinRecForDTM[1, 1,
        ScalarRecurrence -> True,
        ScalarRecurrenceMethod -> "Krylov",
        Verbose -> False],
      {Global`a, Global`x, Global`k}
    }
  ];

  expected22 =
    x^6 - aa[0] x^5
      + (aa[-1] aa[1] - aa[-2] aa[2]) x^4
      + (-aa[-2] aa[1]^2 - aa[-1]^2 aa[2]
          + 2 aa[-2] aa[0] aa[2]) x^3
      + (aa[-2] aa[-1] aa[1] aa[2] - aa[-2]^2 aa[2]^2) x^2
      - aa[-2]^2 aa[0] aa[2]^2 x
      + aa[-2]^3 aa[2]^3;

  expectedSym5 =
    x^5 - (A - C) x^4 - (A C - B^2) x^3
      + C (A C - B^2) x^2 + C^3 (A - C) x - C^5;

  tests = {
    VerificationTest[
      AssociationQ[k11],
      True,
      TestID -> "explicit Krylov dispatch returns an Association"
    ],

    VerificationTest[
      k11["ScalarRecurrenceMethod"],
      "Krylov",
      TestID -> "Krylov method is reported"
    ],

    VerificationTest[
      k11["ObservableOrder"],
      2,
      TestID -> "tridiagonal observable order is two"
    ],

    VerificationTest[
      TRKrylovSamePolynomialQ[
        k11["ScalarRecurrencePolynomial"],
        x^2 - aa[0] x + aa[-1] aa[1]],
      True,
      TestID -> "tridiagonal Krylov recurrence"
    ],

    VerificationTest[
      k22["StateOrdering"],
      "Paper first-discovery order",
      TestID -> "Krylov extension preserves paper state order"
    ],

    VerificationTest[
      k22["ObservableCoordinate"],
      1,
      TestID -> "principal coordinate is first in paper state order"
    ],

    VerificationTest[
      k22["ObservableOrder"],
      6,
      TestID -> "generic pentadiagonal observable order is six"
    ],

    VerificationTest[
      TRKrylovSamePolynomialQ[k22["ScalarRecurrencePolynomial"], expected22],
      True,
      TestID -> "pentadiagonal Krylov recurrence equals Sweet polynomial"
    ],

    VerificationTest[
      Factor[Det[k22["ObservableMatrix"]]],
      -aa[2]^5 aa[-2]^7 (aa[1]^2 aa[-2] - aa[2] aa[-1]^2),
      TestID -> "pentadiagonal observability determinant"
    ],

    VerificationTest[
      Simplify[
        k22["ObservableMatrix"].k22["TransferMatrixExpression"] -
        Normal[k22["KrylovCompanionMatrix"]].k22["ObservableMatrix"]
      ],
      ConstantArray[0, {6, 6}],
      TestID -> "full-order Krylov intertwining identity"
    ],

    VerificationTest[
      k22["SimilarityMatrix"],
      k22["ObservableMatrix"],
      TestID -> "full observable matrix is the Krylov similarity matrix"
    ],

    VerificationTest[
      TRKrylovSamePolynomialQ[
        k22["ScalarRecurrencePolynomial"],
        c22["CharacteristicPolynomial"]],
      True,
      TestID -> "generic pentadiagonal Krylov and characteristic methods agree"
    ],

    VerificationTest[
      kp22["ObservableOrder"],
      6,
      TestID -> "permanent pentadiagonal Krylov scalarization is available"
    ],

    VerificationTest[
      ks22["ObservableOrder"],
      5,
      TestID -> "symmetric pentadiagonal observable order drops to five"
    ],

    VerificationTest[
      TRKrylovSamePolynomialQ[ks22["ScalarRecurrencePolynomial"], expectedSym5],
      True,
      TestID -> "symmetric pentadiagonal Krylov recurrence has degree five"
    ],

    VerificationTest[
      ks22["SimilarityMatrix"],
      Missing["ObservableQuotient"],
      TestID -> "order drop returns an observable quotient"
    ],

    VerificationTest[
      skipped["ScalarRecurrencePolynomial"],
      Missing["NotComputed"],
      TestID -> "Krylov scalarization obeys ScalarRecurrence False"
    ],

    VerificationTest[
      pos22["ScalarRecurrencePolynomial"],
      Missing["PositionDependent"],
      TestID -> "Krylov leaves position-dependent transfer unsimplified"
    ],
    VerificationTest[
      {manualKrylov["ObservableCoordinate"], manualKrylov["ObservableOrder"]},
      {2, 3},
      TestID -> "incremental Krylov supports a non-first observable coordinate"
    ],
    VerificationTest[
      TRKrylovSamePolynomialQ[
        manualKrylov["ScalarRecurrencePolynomial"],
        x^3 - 10 x^2 + 27 x - 58],
      True,
      TestID -> "incremental Krylov non-first-coordinate recurrence"
    ],
    VerificationTest[
      TRKrylovSamePolynomialQ[
        scopeKrylov["ScalarRecurrencePolynomial"],
        x^2 - aa[0] x + aa[-1] aa[1]],
      True,
      TestID -> "Krylov defaults are isolated from Global a/x/k assignments"
    ],
    VerificationTest[
      scopeGlobals, {9, 17, 23},
      TestID -> "Krylov constructor does not modify user Global values"
    ]
  };

  TestReport[tests]
];
