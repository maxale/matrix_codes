(* ::Package:: *)

(* Paper-facing examples for the optional determinant symmetry reductions. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "src", "ToeplitzRecurrencesAll.wl"}]];

ClearAll[RunToeplitzSymmetryExamples];
RunToeplitzSymmetryExamples[] := Module[{s2, s4, k3},
  s2 = SymmetricDeterminantRecurrence[2,
    DiagonalValues -> {2, 3, 5},
    Verbose -> False
  ];

  s4 = SymmetricDeterminantRecurrence[4,
    DiagonalValues -> {2, 3, 5, 7, 11},
    Verbose -> False
  ];

  k3 = SkewSymmetricDeterminantRecurrence[3,
    DiagonalValues -> {3, 4, 5},
    Verbose -> False
  ];

  <|
    "SymmetricM2" -> <|
      "Unrestricted" -> s2["UnrestrictedStateOrder"],
      "Catalan" -> s2["ReducedStateOrder"],
      "Observable" -> s2["ObservableOrder"],
      "Predicted" -> s2["PredictedGenericDegree"]
    |>,
    "SymmetricM4" -> <|
      "Unrestricted" -> s4["UnrestrictedStateOrder"],
      "Catalan" -> s4["ReducedStateOrder"],
      "Observable" -> s4["ObservableOrder"],
      "Predicted" -> s4["PredictedGenericDegree"]
    |>,
    "SkewM3" -> <|
      "Unrestricted" -> k3["UnrestrictedStateOrder"],
      "TheoreticalHodgeHalf" -> k3["StructuralHodgeHalfOrder"],
      "TwoStepObservable" -> k3["TwoStepObservableOrder"],
      "FullAnnihilator" -> k3["FullAnnihilatorOrder"],
      "PredictedTwoStep" -> k3["PredictedTwoStepGenericDegree"],
      "PredictedFull" -> k3["PredictedGenericDegree"]
    |>
  |>
];
