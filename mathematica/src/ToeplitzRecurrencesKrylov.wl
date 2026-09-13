(* ::Package:: *)

(* Observable/Krylov scalarization extension for ToeplitzRecurrences.wl,
   accompanying "Constructive recurrences for determinants and permanents of banded Toeplitz matrices".
   Load the base package first.  This version deliberately avoids OptionValue
   inside Condition (/;) dispatch rules, which is not robust across kernels.
   The public API remains
       ScalarRecurrenceMethod -> "CharacteristicPolynomial" | "Krylov".
*)

BeginPackage["ToeplitzRecurrences`"];

ScalarRecurrenceMethod::usage =
  "ScalarRecurrenceMethod selects the constant-transfer scalarization method: \"CharacteristicPolynomial\" (default) or \"Krylov\".";

Begin["`Private`"];

LinRecForDTM::scalarmethod =
  "ScalarRecurrenceMethod must be \"CharacteristicPolynomial\" or \"Krylov\".";
LinRecForPTM::scalarmethod = LinRecForDTM::scalarmethod;
LinRecForDTM::scalar =
  "ScalarRecurrence must be Automatic, True, or False.";
LinRecForPTM::scalar = LinRecForDTM::scalar;
LinRecForDTM::cutoff =
  "MaxCharacteristicOrder must be a nonnegative integer when ScalarRecurrence->Automatic.";
LinRecForPTM::cutoff = LinRecForDTM::cutoff;
LinRecForDTM::krylovpos =
  "Krylov scalarization is available only for constant transfer matrices; returning the position-dependent transfer system without scalarization.";
LinRecForPTM::krylovpos = LinRecForDTM::krylovpos;

(* Add the selector to the public option sets. *)
Options[LinRecForDTM] = Join[
  DeleteCases[Options[LinRecForDTM], ScalarRecurrenceMethod -> _],
  {ScalarRecurrenceMethod -> "CharacteristicPolynomial"}
];
Options[LinRecForPTM] = Join[
  DeleteCases[Options[LinRecForPTM], ScalarRecurrenceMethod -> _],
  {ScalarRecurrenceMethod -> "CharacteristicPolynomial"}
];

ClearAll[trKrylovShouldScalarize];
trKrylovShouldScalarize[policy_, ord_, cutoff_] := Switch[policy,
  True, True,
  False, False,
  Automatic, IntegerQ[cutoff] && cutoff >= 0 && ord <= cutoff,
  _, $Failed
];

ClearAll[trPrincipalCoordinate];
trPrincipalCoordinate[states_List] := Module[{p},
  p = FirstPosition[
    states,
    s_ /; MatchQ[s, {_List, _List}] && Length[s[[1]]] === 1 && Length[s[[2]]] === 1,
    Missing["NotFound"]
  ];
  If[MissingQ[p], 1, First[p]]
];

ClearAll[trKrylovZeroQ];
trKrylovZeroQ[expr_] := TrueQ[expr === 0] || TrueQ[PossibleZeroQ[expr]];

(* Incremental fraction-free row elimination with relation certificates.
   The pivot order starts with the observable coordinate, so the initial
   coordinate row is the first echelon pivot even if that coordinate is not
   column 1.  Accepted rows remain denominator-free whenever the transfer
   entries are polynomial. *)
ClearAll[trKrylovIncrementalInsert];
trKrylovIncrementalInsert[
    row_List, basisRows_List, basisReps_List, pivots_List,
    rowCount_Integer, columnOrder_List] := Module[
  {work = row, reps, rep, factor, pivotValue, pivot, i},

  reps = (PadRight[#, rowCount + 1, 0] &) /@ basisReps;
  rep = UnitVector[rowCount + 1, rowCount + 1];

  Do[
    factor = work[[pivots[[i]]]];
    If[!trKrylovZeroQ[factor],
      pivotValue = basisRows[[i, pivots[[i]]]];
      work = MapThread[
        Expand[pivotValue #1 - factor #2] &,
        {work, basisRows[[i]]}
      ];
      rep = MapThread[
        Expand[pivotValue #1 - factor #2] &,
        {rep, reps[[i]]}
      ]
    ],
    {i, Length[pivots]}
  ];

  pivot = SelectFirst[
    columnOrder,
    !trKrylovZeroQ[work[[#]]] &,
    Missing["NotFound"]
  ];

  If[MissingQ[pivot],
    <|
      "Dependent" -> True,
      "Relation" -> rep,
      "BasisRows" -> basisRows,
      "BasisRepresentations" -> reps,
      "Pivots" -> pivots
    |>,
    <|
      "Dependent" -> False,
      "BasisRows" -> Append[basisRows, work],
      "BasisRepresentations" -> Append[reps, rep],
      "Pivots" -> Append[pivots, pivot]
    |>
  ]
];

ClearAll[trKrylovScalarization];
trKrylovScalarization[matrix_, var_, coordinate_Integer] := Module[
  {n, row, rows, next, basisRows, basisReps, pivots, columnOrder,
   insert, rel, coeffs, r, poly, companionRules, companion,
   observable, similarity, result = $Failed},

  n = First[Dimensions[matrix]];
  row = UnitVector[n, coordinate];
  rows = {row};

  (* The first Krylov row is already an echelon row.  Subsequent rows are
     inserted incrementally instead of recomputing NullSpace on the whole
     growing Krylov matrix at every step. *)
  basisRows = {row};
  basisReps = {{1}};
  pivots = {coordinate};
  columnOrder = Join[{coordinate}, DeleteCases[Range[n], coordinate]];

  Do[
    next = Expand /@ Normal[Last[rows].matrix];
    insert = trKrylovIncrementalInsert[
      next, basisRows, basisReps, pivots, Length[rows], columnOrder
    ];

    If[TrueQ[insert["Dependent"]],
      rel = insert["Relation"];
      If[Length[rel] === Length[rows] + 1 &&
          !trKrylovZeroQ[Last[rel]],
        coeffs = (Cancel[Together[#]] &) /@ (-Most[rel]/Last[rel]);
        r = Length[coeffs];
        poly = Expand[var^r - Sum[coeffs[[j + 1]] var^j, {j, 0, r - 1}]];
        observable = rows;
        companionRules = Join[
          Table[{j, j + 1} -> 1, {j, 1, r - 1}],
          Table[{r, j} -> coeffs[[j]], {j, 1, r}]
        ];
        companion = SparseArray[companionRules, {r, r}];
        similarity = If[r === n, observable, Missing["ObservableQuotient"]];

        result = <|
          "ObservableOrder" -> r,
          "ObservableCoordinate" -> coordinate,
          "ObservableMatrix" -> observable,
          "KrylovRelationCoefficients" -> coeffs,
          "ScalarRecurrencePolynomial" -> poly,
          "KrylovCompanionMatrix" -> companion,
          "SimilarityMatrix" -> similarity
        |>
      ];
      Break[]
    ];

    basisRows = insert["BasisRows"];
    basisReps = insert["BasisRepresentations"];
    pivots = insert["Pivots"];
    AppendTo[rows, next],
    {n}
  ];

  result
];

(* Call the already-loaded base engine directly.  In particular, do not call
   LinRecForDTM/PTM recursively with a selector option. *)
ClearAll[trKrylovWrapper];
trKrylovWrapper[head_, quantity_, m1_, m2_, diag_, policy_, var_, cutoff_, verbose_,
    posdep_, bandfun_] := Module[
  {base, ord, should, coord, kd, result},

  base = trRowColumnEngine[head, quantity, m1, m2,
    DiagonalValues -> diag,
    ScalarRecurrence -> False,
    CharacteristicPolynomialVariable -> var,
    MaxCharacteristicOrder -> cutoff,
    Verbose -> False,
    PositionDependent -> posdep,
    BandEntryFunction -> bandfun
  ];
  If[base === $Failed || FailureQ[base], Return[base]];

  ord = base["StateOrder"];
  If[!MemberQ[{Automatic, True, False}, policy],
    Message[head::scalar];
    Return[$Failed]
  ];
  If[policy === Automatic && (!IntegerQ[cutoff] || cutoff < 0),
    Message[head::cutoff];
    Return[$Failed]
  ];
  should = trKrylovShouldScalarize[policy, ord, cutoff];

  If[TrueQ[posdep],
    If[policy === True, Message[head::krylovpos]];
    result = Join[base, <|
      "ScalarRecurrenceMethod" -> "Krylov",
      "CharacteristicPolynomial" -> Missing["PositionDependent"],
      "ScalarRecurrencePolynomial" -> Missing["PositionDependent"],
      "ObservableOrder" -> Missing["PositionDependent"],
      "ObservableCoordinate" -> Missing["PositionDependent"],
      "ObservableMatrix" -> Missing["PositionDependent"],
      "KrylovRelationCoefficients" -> Missing["PositionDependent"],
      "KrylovCompanionMatrix" -> Missing["PositionDependent"],
      "SimilarityMatrix" -> Missing["PositionDependent"]
    |>];
    If[TrueQ[verbose],
      Print["The transfer-matrix order is ", ord, "."];
      Print["The position-dependent transfer matrix Q[k] is:"];
      Print[MatrixForm[Normal[base["TransferMatrixExpression"]]]];
      Print["Krylov scalarization was not performed because the transfer matrix depends on k."]
    ];
    Return[result]
  ];

  If[!TrueQ[should],
    result = Join[base, <|
      "ScalarRecurrenceMethod" -> "Krylov",
      "CharacteristicPolynomial" -> Missing["NotComputed"],
      "ScalarRecurrencePolynomial" -> Missing["NotComputed"],
      "ObservableOrder" -> Missing["NotComputed"],
      "ObservableCoordinate" -> trPrincipalCoordinate[base["States"]],
      "ObservableMatrix" -> Missing["NotComputed"],
      "KrylovRelationCoefficients" -> Missing["NotComputed"],
      "KrylovCompanionMatrix" -> Missing["NotComputed"],
      "SimilarityMatrix" -> Missing["NotComputed"]
    |>];
    If[TrueQ[verbose],
      Print["The transfer-matrix order is ", ord, "."];
      Print["The transfer matrix is:"];
      Print[MatrixForm[Normal[base["TransferMatrixExpression"]]]];
      Print["The Krylov scalar recurrence was not computed."]
    ];
    Return[result]
  ];

  coord = trPrincipalCoordinate[base["States"]];
  kd = trKrylovScalarization[base["TransferMatrixExpression"], var, coord];
  If[kd === $Failed, Return[$Failed]];

  result = Join[base, kd, <|
    "ScalarRecurrenceMethod" -> "Krylov",
    "CharacteristicPolynomial" -> Missing["NotComputedByKrylov"]
  |>];

  If[TrueQ[verbose],
    Print["The transfer-matrix order is ", ord, "."];
    Print["The transfer matrix is:"];
    Print[MatrixForm[Normal[base["TransferMatrixExpression"]]]];
    Print["The observable Krylov order is ", result["ObservableOrder"], "."];
    Print["The scalar recurrence polynomial from the determinant/permanent coordinate is:"];
    Print[result["ScalarRecurrencePolynomial"]]
  ];

  result
];

(* Replace only the two public wrappers.  Options and the base private engine
   remain intact.  Keeping all option dispatch in the RHS avoids the 14.2.1
   failure caused by OptionValue inside /; conditions. *)
DownValues[LinRecForDTM] = {};
DownValues[LinRecForPTM] = {};

LinRecForDTM[m1_Integer?NonNegative, m2_Integer?NonNegative,
    OptionsPattern[]] := Module[{method = OptionValue[ScalarRecurrenceMethod]},
  Switch[method,
    "CharacteristicPolynomial",
      trRowColumnEngine[LinRecForDTM, "Determinant", m1, m2,
        DiagonalValues -> OptionValue[DiagonalValues],
        ScalarRecurrence -> OptionValue[ScalarRecurrence],
        CharacteristicPolynomialVariable -> OptionValue[CharacteristicPolynomialVariable],
        MaxCharacteristicOrder -> OptionValue[MaxCharacteristicOrder],
        Verbose -> OptionValue[Verbose],
        PositionDependent -> OptionValue[PositionDependent],
        BandEntryFunction -> OptionValue[BandEntryFunction]
      ],
    "Krylov",
      trKrylovWrapper[LinRecForDTM, "Determinant", m1, m2,
        OptionValue[DiagonalValues],
        OptionValue[ScalarRecurrence],
        OptionValue[CharacteristicPolynomialVariable],
        OptionValue[MaxCharacteristicOrder],
        OptionValue[Verbose],
        OptionValue[PositionDependent],
        OptionValue[BandEntryFunction]
      ],
    _,
      Message[LinRecForDTM::scalarmethod];
      $Failed
  ]
];

LinRecForPTM[m1_Integer?NonNegative, m2_Integer?NonNegative,
    OptionsPattern[]] := Module[{method = OptionValue[ScalarRecurrenceMethod]},
  Switch[method,
    "CharacteristicPolynomial",
      trRowColumnEngine[LinRecForPTM, "Permanent", m1, m2,
        DiagonalValues -> OptionValue[DiagonalValues],
        ScalarRecurrence -> OptionValue[ScalarRecurrence],
        CharacteristicPolynomialVariable -> OptionValue[CharacteristicPolynomialVariable],
        MaxCharacteristicOrder -> OptionValue[MaxCharacteristicOrder],
        Verbose -> OptionValue[Verbose],
        PositionDependent -> OptionValue[PositionDependent],
        BandEntryFunction -> OptionValue[BandEntryFunction]
      ],
    "Krylov",
      trKrylovWrapper[LinRecForPTM, "Permanent", m1, m2,
        OptionValue[DiagonalValues],
        OptionValue[ScalarRecurrence],
        OptionValue[CharacteristicPolynomialVariable],
        OptionValue[MaxCharacteristicOrder],
        OptionValue[Verbose],
        OptionValue[PositionDependent],
        OptionValue[BandEntryFunction]
      ],
    _,
      Message[LinRecForPTM::scalarmethod];
      $Failed
  ]
];

End[];
EndPackage[];
