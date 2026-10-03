(* Minimal regression reproducer for the Krylov scalarization path. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "src", "ToeplitzRecurrences.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "src", "ToeplitzRecurrencesKrylov.wl"}]];

rDefault = LinRecForDTM[1, 1, ScalarRecurrence -> False, Verbose -> False];
rKrylov = LinRecForDTM[1, 1,
  ScalarRecurrence -> True,
  ScalarRecurrenceMethod -> "Krylov",
  Verbose -> False];

Print["$Version = ", $Version];
Print["default dispatch AssociationQ = ", AssociationQ[rDefault]];
Print["Krylov dispatch AssociationQ = ", AssociationQ[rKrylov]];
If[AssociationQ[rKrylov],
  Print["Krylov method = ", rKrylov["ScalarRecurrenceMethod"]];
  Print["observable order = ", rKrylov["ObservableOrder"]];
];

If[
  AssociationQ[rDefault] && AssociationQ[rKrylov] &&
  rKrylov["ScalarRecurrenceMethod"] === "Krylov" &&
  rKrylov["ObservableOrder"] === 2,
  Print["SMOKE TEST: PASS"],
  Print["SMOKE TEST: FAIL"];
  $Failed
]
