Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "src", "ToeplitzRecurrencesAll.wl"}]];

s = SymmetricDeterminantRecurrence[2,
  DiagonalValues -> {2, 3, 5},
  Verbose -> False
];
If[!AssociationQ[s] || s["ReducedStateOrder"] =!= 5 || s["ObservableOrder"] =!= 5,
  Print["SymmetryReductionSmokeTest: FAIL"];
  Exit[1]
];

k = SkewSymmetricDeterminantRecurrence[2,
  DiagonalValues -> {3, 4},
  Verbose -> False
];
If[!AssociationQ[k] || k["ObservableOrder"] =!= 6 || k["EvenObservableOrder"] =!= 3,
  Print["SymmetryReductionSmokeTest: FAIL"];
  Exit[1]
];

Print["SymmetryReductionSmokeTest: PASS"];
