(* ::Package:: *)

(* Regression checks for the optional symmetry-aware determinant reductions. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "src", "ToeplitzRecurrencesAll.wl"}]];

ClearAll[TRSymmetrySamePolynomialQ, RunToeplitzSymmetryReductionTests];
TRSymmetrySamePolynomialQ[p_, q_] := TrueQ[Expand[p - q] === 0];

RunToeplitzSymmetryReductionTests[] := Module[
  {x = \[FormalX], A, B, C, s2, full2, s3, full3, k1, k2, k3,
   k1direct, k2direct, k3direct, straight4, u4, x4, straight4Reconstructs, tests},

  s2 = SymmetricDeterminantRecurrence[2,
    DiagonalValues -> {A, B, C},
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];
  full2 = LinRecForDTM[2, 2,
    DiagonalValues -> {C, B, A, B, C},
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];

  s3 = SymmetricDeterminantRecurrence[3,
    DiagonalValues -> {2, 3, 5, 7},
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];
  full3 = LinRecForDTM[3, 3,
    DiagonalValues -> {7, 5, 3, 2, 3, 5, 7},
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];

  straight4 = ToeplitzRecurrences`Private`trSymmetricStraightenPair[
    4, {{2, 3}, {1, 4}}
  ];

  u4 = Table[Unique["trTestSym$"], {10}];
  x4 = {
    {u4[[1]], u4[[2]], u4[[3]], u4[[4]]},
    {u4[[2]], u4[[5]], u4[[6]], u4[[7]]},
    {u4[[3]], u4[[6]], u4[[8]], u4[[9]]},
    {u4[[4]], u4[[7]], u4[[9]], u4[[10]]}
  };
  straight4Reconstructs = TrueQ[
    Expand[
      Det[x4[[{2, 3}, {1, 4}]]] -
      Total[
        (#[[2]] Det[x4[[#[[1, 1]], #[[1, 2]]]]]) & /@
          (List @@@ straight4)
      ]
    ] === 0
  ];

  k1 = SkewSymmetricDeterminantRecurrence[1,
    DiagonalValues -> {3},
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];
  k2 = SkewSymmetricDeterminantRecurrence[2,
    DiagonalValues -> {3, 4},
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];
  k3 = SkewSymmetricDeterminantRecurrence[3,
    DiagonalValues -> {3, 4, 5},
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];

  k1direct = LinRecForDTM[1, 1,
    DiagonalValues -> {-3, 0, 3},
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];
  k2direct = LinRecForDTM[2, 2,
    DiagonalValues -> {-4, -3, 0, 3, 4},
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];
  k3direct = LinRecForDTM[3, 3,
    DiagonalValues -> {-5, -4, -3, 0, 3, 4, 5},
    ScalarRecurrence -> True,
    ScalarRecurrenceMethod -> "Krylov",
    CharacteristicPolynomialVariable -> x,
    Verbose -> False
  ];

  tests = {
    VerificationTest[
      {s2["UnrestrictedStateOrder"], s2["ReducedStateOrder"],
       s2["ObservableOrder"], s2["PredictedGenericDegree"]},
      {6, 5, 5, 5},
      TestID -> "symmetric pentadiagonal reduction is 6 to 5 to 5"
    ],
    VerificationTest[
      s2["ReducedLevelCounts"],
      {1, 3, 1},
      TestID -> "symmetric m2 doset levels are Narayana counts"
    ],
    VerificationTest[
      TRSymmetrySamePolynomialQ[
        s2["ScalarRecurrencePolynomial"],
        full2["ScalarRecurrencePolynomial"]
      ],
      True,
      TestID -> "symmetric m2 reduced and unrestricted principal Krylov polynomials agree"
    ],
    VerificationTest[
      {s3["UnrestrictedStateOrder"], s3["ReducedStateOrder"],
       s3["ObservableOrder"], s3["PredictedGenericDegree"]},
      {20, 14, 14, 14},
      TestID -> "symmetric m3 reduction is 20 to 14 to 14"
    ],
    VerificationTest[
      s3["ReducedLevelCounts"],
      {1, 6, 6, 1},
      TestID -> "symmetric m3 doset levels are Narayana counts"
    ],
    VerificationTest[
      TRSymmetrySamePolynomialQ[
        s3["ScalarRecurrencePolynomial"],
        full3["ScalarRecurrencePolynomial"]
      ],
      True,
      TestID -> "symmetric m3 reduced and unrestricted principal Krylov polynomials agree"
    ],
    VerificationTest[
      straight4,
      {
        {{1, 2}, {3, 4}} -> -1,
        {{1, 3}, {2, 4}} -> 1
      },
      TestID -> "nontrivial m4 doset straightening relation"
    ],
    VerificationTest[
      straight4Reconstructs,
      True,
      TestID -> "m4 doset straightening reconstructs the original minor"
    ],
    VerificationTest[
      {k1["Method"], k1["HodgeHalfConstructed"],
       k1["FullAnnihilatorOrder"], k1["TwoStepObservableOrder"]},
      {"SkewTwoStepKrylov", False, 2, 1},
      TestID -> "skew m1 uses Q2 Krylov and has orders 2 and 1"
    ],
    VerificationTest[
      {k2["FullAnnihilatorOrder"], k2["TwoStepObservableOrder"]},
      {6, 3},
      TestID -> "skew m2 full and two-step observable orders"
    ],
    VerificationTest[
      {k3["FullAnnihilatorOrder"], k3["TwoStepObservableOrder"]},
      {18, 9},
      TestID -> "skew m3 full and two-step observable orders"
    ],
    VerificationTest[
      And[
        Normal[k1["TwoStepTransferMatrixExpression"]] ===
          Normal[k1["TransferMatrixExpression"].k1["TransferMatrixExpression"]],
        Normal[k2["TwoStepTransferMatrixExpression"]] ===
          Normal[k2["TransferMatrixExpression"].k2["TransferMatrixExpression"]],
        Normal[k3["TwoStepTransferMatrixExpression"]] ===
          Normal[k3["TransferMatrixExpression"].k3["TransferMatrixExpression"]]
      ],
      True,
      TestID -> "skew two-step transfer is Q squared"
    ],
    VerificationTest[
      And[
        TRSymmetrySamePolynomialQ[k1["FullAnnihilatorPolynomial"], k1direct["ScalarRecurrencePolynomial"]],
        TRSymmetrySamePolynomialQ[k2["FullAnnihilatorPolynomial"], k2direct["ScalarRecurrencePolynomial"]],
        TRSymmetrySamePolynomialQ[k3["FullAnnihilatorPolynomial"], k3direct["ScalarRecurrencePolynomial"]]
      ],
      True,
      TestID -> "skew lifted Q2 annihilator agrees with direct full Krylov"
    ],
    VerificationTest[
      And[
        k1["FullAttainsPredictedDegree"], k1["TwoStepAttainsPredictedDegree"],
        k2["FullAttainsPredictedDegree"], k2["TwoStepAttainsPredictedDegree"],
        k3["FullAttainsPredictedDegree"], k3["TwoStepAttainsPredictedDegree"]
      ],
      True,
      TestID -> "skew numerical witnesses attain predicted generic degrees"
    ]
  };

  TestReport[tests]
];
