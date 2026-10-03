(* ::Package:: *)

(* Optional symmetry-aware determinant reductions accompanying
   "Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents"
   by Max A. Alekseyev and Dmitry I. Khomovsky.

   This module deliberately leaves the generic row-column and Krylov engines
   unchanged.  In the symmetric balanced case, each generated boundary minor
   is straightened immediately into the doset basis before the reduced transfer
   is assembled.  In the skew-symmetric balanced case, the normalized
   row-column transfer Q is retained and the optimized principal Krylov
   calculation is applied directly to Q^2; its polynomial R(y) is then lifted
   to the full one-step annihilator R(x^2).
*)

BeginPackage["ToeplitzRecurrences`"];

SymmetricDeterminantRecurrence::usage =
  "SymmetricDeterminantRecurrence[m,opts] constructs the doset/Catalan row-column transfer for a symmetric balanced (m,m)-banded Toeplitz determinant and computes its principal Krylov recurrence.";
SkewSymmetricDeterminantRecurrence::usage =
  "SkewSymmetricDeterminantRecurrence[m,opts] keeps the full normalized skew row-column transfer Q, applies principal Krylov scalarization to Q^2, and lifts the resulting two-step polynomial R(y) to the full annihilator R(x^2).";

Begin["`Private`"];

SymmetricDeterminantRecurrence::diag =
  "DiagonalValues must be Automatic, a compact list {a0,a1,...,am}, or a full symmetric list in offset order -m,...,m.";
SymmetricDeterminantRecurrence::notsym =
  "The supplied full DiagonalValues list is not symmetric.";
SkewSymmetricDeterminantRecurrence::diag =
  "DiagonalValues must be Automatic, a compact list {a1,...,am}, or a full skew-symmetric list in offset order -m,...,m.";
SkewSymmetricDeterminantRecurrence::notskew =
  "The supplied full DiagonalValues list is not skew-symmetric with zero main diagonal.";
SkewSymmetricDeterminantRecurrence::evenodd =
  "The principal observable polynomial was not even, so the skew even-subsequence reduction could not be formed.";

Options[SymmetricDeterminantRecurrence] = {
  DiagonalValues -> Automatic,
  CharacteristicPolynomialVariable -> \[FormalX],
  Verbose -> True
};
Options[SkewSymmetricDeterminantRecurrence] = Options[SymmetricDeterminantRecurrence];

ClearAll[trSymmetryZeroQ];
trSymmetryZeroQ[expr_] := TrueQ[expr === 0] || TrueQ[PossibleZeroQ[expr]];

ClearAll[trPairKey];
trPairKey[pair_] := ToString[pair, InputForm];

ClearAll[trDosetPairQ];
trDosetPairQ[{A_List, B_List}] := And @@ MapThread[LessEqual, {A, B}];

ClearAll[trAllStatePairs, trDosetPairs];
trAllStatePairs[m_Integer?Positive] := Flatten[
  Table[
    With[{subs = Subsets[Range[m], {j}]}, Tuples[{subs, subs}]],
    {j, 0, m}
  ],
  1
];
trDosetPairs[m_Integer?Positive] := Select[trAllStatePairs[m], trDosetPairQ];

(* Pair labels <-> the normalized balanced row-column signatures used by the
   base row-column Laplace engine.  These are the Sigma(A,B) coordinates of
   the symmetry paper. *)
ClearAll[trStateSignatureFromPair, trStatePairFromSignature];
trStateSignatureFromPair[{A_List, B_List}, m_Integer?Positive] := Module[
  {j, ell, ahat, numericRows, numericCols},
  j = Length[A];
  ell = j + 1;
  ahat = If[
    MemberQ[B, 1],
    Sort@Join[DeleteCases[A, 1], If[MemberQ[A, 1], {m + 1}, {}]],
    A
  ];
  numericRows = Prepend[ell - B, ell];
  numericCols = Prepend[ell - ahat, ell];
  {
    (\[FormalK] + #) & /@ numericRows,
    (\[FormalK] + #) & /@ numericCols
  }
];

trStatePairFromSignature[state : {_List, _List}, m_Integer?Positive] := Module[
  {rows, cols, ell, numericRows, numericCols, xs, ys, A, B},
  rows = state[[1]];
  cols = state[[2]];
  ell = Length[rows];
  If[Length[cols] =!= ell, Return[$Failed]];
  numericRows = Simplify[# - \[FormalK]] & /@ rows;
  numericCols = Simplify[# - \[FormalK]] & /@ cols;
  xs = Sort[Rest[numericRows]];
  ys = Sort[Rest[numericCols]];
  B = Sort[ell - # & /@ xs];
  A = Sort[Mod[ell - # - 1, m] + 1 & /@ ys];
  {A, B}
];

ClearAll[trHeldCoefficientAssociation];
trHeldCoefficientAssociation[poly_, vars_List] := Association[
  Map[
    Function[rule,
      With[{parts = List @@ rule}, parts[[1]] -> parts[[2]]]
    ],
    CoefficientRules[Expand[poly], vars]
  ]
];

(* Universal straightening data for minors of a generic symmetric m x m
   matrix.  At each level j, standard pairs A<=B form the doset basis. *)
ClearAll[trSymmetricLevelStraightening, trSymmetricStraighteningData];
trSymmetricLevelStraightening[m_Integer?Positive, j_Integer?NonNegative] := Module[
  {subs, pairs, basis, vars, varPairs, varPos, xmat, polys, coeffAssocs,
   monomials, vectors, basisPositions, basisMatrix, leftInverse,
   coeffVectors, coeffAssoc, i, coeffs, rules},

  subs = Subsets[Range[m], {j}];
  pairs = Tuples[{subs, subs}];
  basis = Select[pairs, trDosetPairQ];

  varPairs = Flatten[Table[{r, c}, {r, 1, m}, {c, r, m}], 1];
  vars = Table[Unique["trSym$"], {Length[varPairs]}];
  varPos = AssociationThread[(ToString[#, InputForm] & /@ varPairs) -> Range[Length[varPairs]]];
  xmat = Table[
    vars[[varPos[ToString[{Min[r, c], Max[r, c]}, InputForm]]]],
    {r, 1, m}, {c, 1, m}
  ];

  polys = If[
    j === 0,
    ConstantArray[1, Length[pairs]],
    Det[xmat[[#[[1]], #[[2]]]]] & /@ pairs
  ];
  coeffAssocs = trHeldCoefficientAssociation[#, vars] & /@ polys;
  monomials = Union @@ (Keys /@ coeffAssocs);
  vectors = (Lookup[#, monomials, 0] &) /@ coeffAssocs;

  basisPositions = Flatten[Position[pairs, p_ /; trDosetPairQ[p]]];
  basisMatrix = Transpose[vectors[[basisPositions]]];
  leftInverse = LinearSolve[
    Transpose[basisMatrix].basisMatrix,
    Transpose[basisMatrix]
  ];
  coeffVectors = (Together /@ (leftInverse.#)) & /@ vectors;

  coeffAssoc = <||>;
  Do[
    coeffs = coeffVectors[[i]];
    rules = DeleteCases[
      Table[
        If[trSymmetryZeroQ[coeffs[[k]]], Nothing, basis[[k]] -> coeffs[[k]]],
        {k, Length[basis]}
      ],
      Nothing
    ];
    AssociateTo[coeffAssoc, trPairKey[pairs[[i]]] -> rules],
    {i, Length[pairs]}
  ];

  <|
    "Basis" -> basis,
    "Coefficients" -> coeffAssoc,
    "LevelCount" -> Length[basis]
  |>
];

trSymmetricStraighteningData[m_Integer?Positive] :=
  trSymmetricStraighteningData[m] = Module[
    {levels, basis, coeffs, levelCounts, expected},
    levels = Table[trSymmetricLevelStraightening[m, j], {j, 0, m}];
    basis = Flatten[levels[[All, "Basis"]], 1];
    coeffs = Join @@ (levels[[All, "Coefficients"]]);
    levelCounts = levels[[All, "LevelCount"]];
    expected = Binomial[2 m + 2, m + 1]/(m + 2);
    If[Length[basis] =!= expected,
      Return[Failure["CatalanDimensionMismatch", <|
        "Actual" -> Length[basis], "Expected" -> expected
      |>]]
    ];
    <|
      "Basis" -> basis,
      "Coefficients" -> coeffs,
      "LevelCounts" -> levelCounts,
      "CatalanDimension" -> expected
    |>
  ];

ClearAll[trSymmetricStraightenPair];
trSymmetricStraightenPair[m_Integer?Positive, pair : {_List, _List}] := Module[
  {data = trSymmetricStraighteningData[m]},
  If[FailureQ[data], Return[data]];
  Lookup[data["Coefficients"], trPairKey[pair], Missing["UnknownPair"]]
];

ClearAll[trSymmetricDiagonalSetup];
trSymmetricDiagonalSetup[m_Integer?Positive, diag_] := Module[
  {values, compact, full, assoc, entry, symmetricQ},
  values = If[
    diag === Automatic,
    Join[{\[FormalA][0]}, Table[\[FormalA][s], {s, 1, m}]],
    diag
  ];
  If[!ListQ[values],
    Message[SymmetricDeterminantRecurrence::diag];
    Return[$Failed]
  ];
  Which[
    Length[values] === m + 1,
      compact = values,
    Length[values] === 2 m + 1,
      symmetricQ = And @@ Table[
        trSymmetryZeroQ[values[[m + 1 - s]] - values[[m + 1 + s]]],
        {s, 1, m}
      ];
      If[!TrueQ[symmetricQ],
        Message[SymmetricDeterminantRecurrence::notsym];
        Return[$Failed]
      ];
      compact = Join[{values[[m + 1]]}, Table[values[[m + 1 + s]], {s, 1, m}]],
    True,
      Message[SymmetricDeterminantRecurrence::diag];
      Return[$Failed]
  ];
  full = Join[Reverse[Rest[compact]], compact];
  assoc = AssociationThread[Range[-m, m] -> full];
  entry = With[{aassoc = assoc}, Function[{s, t}, Lookup[aassoc, s, 0]]];
  <|"Compact" -> compact, "Full" -> full, "EntryFunction" -> entry|>
];

ClearAll[trSkewDiagonalSetup];
trSkewDiagonalSetup[m_Integer?Positive, diag_] := Module[
  {values, compact, full, skewQ},
  values = If[
    diag === Automatic,
    Table[\[FormalA][s], {s, 1, m}],
    diag
  ];
  If[!ListQ[values],
    Message[SkewSymmetricDeterminantRecurrence::diag];
    Return[$Failed]
  ];
  Which[
    Length[values] === m,
      compact = values,
    Length[values] === 2 m + 1,
      skewQ = trSymmetryZeroQ[values[[m + 1]]] && And @@ Table[
        trSymmetryZeroQ[values[[m + 1 - s]] + values[[m + 1 + s]]],
        {s, 1, m}
      ];
      If[!TrueQ[skewQ],
        Message[SkewSymmetricDeterminantRecurrence::notskew];
        Return[$Failed]
      ];
      compact = Table[values[[m + 1 + s]], {s, 1, m}],
    True,
      Message[SkewSymmetricDeterminantRecurrence::diag];
      Return[$Failed]
  ];
  full = Join[-Reverse[compact], {0}, compact];
  <|"Compact" -> compact, "Full" -> full|>
];

(* Direct symmetry-aware row-column construction.  Only doset states are
   stored.  Each Laplace target is normalized, decoded to (A,B), and
   immediately straightened before its transition is entered in the matrix. *)
ClearAll[trSymmetricReducedTransfer];
trSymmetricReducedTransfer[m_Integer?Positive, entry_] := Module[
  {straightening, basis, basisIndex, rawRules = {}, sourceState, data,
   term, targetState, targetPair, expansion, i, r, standardPair, coeff,
   pos, groupedPositions, finalRules, matrix},

  straightening = trSymmetricStraighteningData[m];
  If[FailureQ[straightening], Return[straightening]];
  basis = straightening["Basis"];
  basisIndex = AssociationThread[(trPairKey /@ basis) -> Range[Length[basis]]];

  Do[
    sourceState = trStateSignatureFromPair[basis[[i]], m];
    data = trLaplaceExpansion[
      SymmetricDeterminantRecurrence, m, m, sourceState,
      "Determinant", entry
    ];
    If[data === $Failed, Return[$Failed]];

    Do[
      term = data[[2, r]];
      targetState = trFinalSignature[term[[2]]];
      targetState = MReduce[targetState[[1]], targetState[[2]]];
      targetPair = trStatePairFromSignature[targetState, m];
      If[targetPair === $Failed, Return[$Failed]];
      expansion = trSymmetricStraightenPair[m, targetPair];
      If[MissingQ[expansion] || FailureQ[expansion], Return[$Failed]];
      Do[
        standardPair = First[expansion[[pos]]];
        coeff = Last[expansion[[pos]]];
        AppendTo[
          rawRules,
          {i, basisIndex[trPairKey[standardPair]]} -> Expand[term[[1]] coeff]
        ],
        {pos, Length[expansion]}
      ],
      {r, Length[data[[2]]]}
    ],
    {i, Length[basis]}
  ];

  groupedPositions = DeleteDuplicates[First /@ rawRules];
  finalRules = Table[
    pos -> Expand@Total[Last /@ Select[rawRules, First[#] === pos &]],
    {pos, groupedPositions}
  ];
  matrix = SparseArray[finalRules, {Length[basis], Length[basis]}];

  <|
    "Basis" -> basis,
    "LevelCounts" -> straightening["LevelCounts"],
    "CatalanDimension" -> straightening["CatalanDimension"],
    "TransferMatrixExpression" -> matrix
  |>
];

SymmetricDeterminantRecurrence[m_Integer?Positive, OptionsPattern[]] := Module[
  {setup, reduced, basis, principal, var, kd, predicted, result, verbose},
  setup = trSymmetricDiagonalSetup[m, OptionValue[DiagonalValues]];
  If[setup === $Failed, Return[$Failed]];
  var = OptionValue[CharacteristicPolynomialVariable];
  verbose = TrueQ[OptionValue[Verbose]];

  reduced = trSymmetricReducedTransfer[m, setup["EntryFunction"]];
  If[reduced === $Failed || FailureQ[reduced], Return[reduced]];
  basis = reduced["Basis"];
  principal = First@FirstPosition[basis, {{}, {}}];
  kd = trKrylovScalarization[
    reduced["TransferMatrixExpression"], var, principal
  ];
  If[kd === $Failed, Return[$Failed]];
  predicted = (3^m + 1)/2;

  result = Join[
    <|
      "Method" -> "SymmetricDosetRowColumn",
      "SymmetryClass" -> "Symmetric",
      "PaperTitle" -> "Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents",
      "Semibandwidth" -> m,
      "DiagonalValues" -> setup["Compact"],
      "FullDiagonalValues" -> setup["Full"],
      "UnrestrictedStateOrder" -> Binomial[2 m, m],
      "ReducedBasisType" -> "Doset",
      "ReducedStates" -> basis,
      "ReducedLevelCounts" -> reduced["LevelCounts"],
      "ReducedStateOrder" -> reduced["CatalanDimension"],
      "CatalanDimension" -> reduced["CatalanDimension"],
      "TransferMatrixExpression" -> reduced["TransferMatrixExpression"],
      "PrincipalCoordinate" -> principal,
      "PredictedGenericDegree" -> predicted
    |>,
    kd
  ];
  result = Join[result, <|
    "AttainsPredictedDegree" -> TrueQ[result["ObservableOrder"] === predicted]
  |>];

  If[verbose,
    Print["The unrestricted row-column order is ", result["UnrestrictedStateOrder"], "."];
    Print["The symmetric doset/Catalan transfer order is ", result["ReducedStateOrder"], "."];
    Print["The observable Krylov order is ", result["ObservableOrder"], "."];
    Print["The generic symmetric degree predicted by the theorem is ", predicted, "."]
  ];
  result
];

ClearAll[trSkewTwoStepKrylov, trLiftTwoStepPolynomial];
trSkewTwoStepKrylov[matrix_, var_, principal_Integer?Positive] := Module[
  {twoStep, kd},
  twoStep = matrix.matrix;
  kd = trKrylovScalarization[twoStep, var, principal];
  If[kd === $Failed || FailureQ[kd], Return[kd]];
  <|
    "TwoStepTransferMatrixExpression" -> twoStep,
    "KrylovData" -> kd
  |>
];

trLiftTwoStepPolynomial[poly_, var_] := Module[{degree},
  degree = Exponent[poly, var];
  Expand@Sum[Coefficient[poly, var, j] var^(2 j), {j, 0, degree}]
];

ClearAll[trEvenPolynomial, trCompanionFromMonicPolynomial];
trEvenPolynomial[poly_, var_] := Module[{degree, oddExponents, ypoly},
  degree = Exponent[poly, var];
  oddExponents = Select[Range[1, degree, 2], !trSymmetryZeroQ[Coefficient[poly, var, #]] &];
  If[oddExponents =!= {}, Return[$Failed]];
  ypoly = Expand@Sum[
    Coefficient[poly, var, 2 j] var^j,
    {j, 0, Floor[degree/2]}
  ];
  ypoly
];

trCompanionFromMonicPolynomial[poly_, var_] := Module[
  {r, leading, coeffs, rules, companion},
  r = Exponent[poly, var];
  leading = Coefficient[poly, var, r];
  If[r < 1 || !trSymmetryZeroQ[leading - 1], Return[$Failed]];
  coeffs = -Table[Coefficient[poly, var, j], {j, 0, r - 1}];
  rules = Join[
    Table[{j, j + 1} -> 1, {j, 1, r - 1}],
    Table[{r, j} -> coeffs[[j]], {j, 1, r}]
  ];
  companion = SparseArray[rules, {r, r}];
  <|"Coefficients" -> coeffs, "Companion" -> companion|>
];

SkewSymmetricDeterminantRecurrence[m_Integer?Positive, OptionsPattern[]] := Module[
  {setup, var, verbose, base, transfer, principal, twoStepData, twoStepKrylov,
   twoStepPoly, fullPoly, fullData, predictedFull, predictedTwoStep,
   twoStepOrder, fullOrder, hodgeHalfOrder, result},

  setup = trSkewDiagonalSetup[m, OptionValue[DiagonalValues]];
  If[setup === $Failed, Return[$Failed]];
  var = OptionValue[CharacteristicPolynomialVariable];
  verbose = TrueQ[OptionValue[Verbose]];

  base = LinRecForDTM[m, m,
    DiagonalValues -> setup["Full"],
    ScalarRecurrence -> False,
    ScalarRecurrenceMethod -> "CharacteristicPolynomial",
    CharacteristicPolynomialVariable -> var,
    MaxCharacteristicOrder -> Binomial[2 m, m],
    Verbose -> False
  ];
  If[base === $Failed || FailureQ[base], Return[base]];

  transfer = base["TransferMatrixExpression"];
  principal = trPrincipalCoordinate[base["States"]];
  twoStepData = trSkewTwoStepKrylov[transfer, var, principal];
  If[twoStepData === $Failed || FailureQ[twoStepData], Return[twoStepData]];
  twoStepKrylov = twoStepData["KrylovData"];
  twoStepPoly = twoStepKrylov["ScalarRecurrencePolynomial"];
  fullPoly = trLiftTwoStepPolynomial[twoStepPoly, var];
  fullData = trCompanionFromMonicPolynomial[fullPoly, var];
  If[fullData === $Failed, Return[$Failed]];

  predictedFull = 2 3^(m - 1);
  predictedTwoStep = 3^(m - 1);
  twoStepOrder = twoStepKrylov["ObservableOrder"];
  fullOrder = Exponent[fullPoly, var];
  hodgeHalfOrder = Binomial[2 m, m]/2;

  result = <|
    "Method" -> "SkewTwoStepKrylov",
    "SymmetryClass" -> "SkewSymmetric",
    "PaperTitle" -> "Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents",
    "Semibandwidth" -> m,
    "DiagonalValues" -> setup["Compact"],
    "FullDiagonalValues" -> setup["Full"],
    "UnrestrictedStateOrder" -> Binomial[2 m, m],
    "StructuralHodgeHalfOrder" -> hodgeHalfOrder,
    "StructuralParityStateOrder" -> hodgeHalfOrder,
    "HodgeHalfConstructed" -> False,
    "TransferMatrixExpression" -> transfer,
    "TwoStepTransferMatrixExpression" -> twoStepData["TwoStepTransferMatrixExpression"],
    "PrincipalCoordinate" -> principal,
    "TwoStepObservableOrder" -> twoStepOrder,
    "TwoStepScalarRecurrencePolynomial" -> twoStepPoly,
    "TwoStepObservableMatrix" -> twoStepKrylov["ObservableMatrix"],
    "TwoStepKrylovRelationCoefficients" -> twoStepKrylov["KrylovRelationCoefficients"],
    "TwoStepKrylovCompanionMatrix" -> twoStepKrylov["KrylovCompanionMatrix"],
    "FullAnnihilatorOrder" -> fullOrder,
    "FullAnnihilatorPolynomial" -> fullPoly,
    "FullAnnihilatorRelationCoefficients" -> fullData["Coefficients"],
    "FullAnnihilatorCompanionMatrix" -> fullData["Companion"],
    "ObservableOrder" -> fullOrder,
    "ScalarRecurrencePolynomial" -> fullPoly,
    "KrylovCompanionMatrix" -> fullData["Companion"],
    "EvenObservableOrder" -> twoStepOrder,
    "EvenScalarRecurrencePolynomial" -> twoStepPoly,
    "EvenKrylovRelationCoefficients" -> twoStepKrylov["KrylovRelationCoefficients"],
    "EvenObservableTransfer" -> twoStepKrylov["KrylovCompanionMatrix"],
    "PredictedGenericDegree" -> predictedFull,
    "PredictedEvenGenericDegree" -> predictedTwoStep,
    "PredictedTwoStepGenericDegree" -> predictedTwoStep,
    "FullAttainsPredictedDegree" -> TrueQ[fullOrder === predictedFull],
    "EvenAttainsPredictedDegree" -> TrueQ[twoStepOrder === predictedTwoStep],
    "TwoStepAttainsPredictedDegree" -> TrueQ[twoStepOrder === predictedTwoStep]
  |>;

  If[verbose,
    Print["The unrestricted row-column order is ", result["UnrestrictedStateOrder"], "."];
    Print["The theoretical Hodge-half dimension is ", hodgeHalfOrder,
      " (not constructed as a local row-column quotient)."];
    Print["The two-step principal Krylov order is ", twoStepOrder, "."];
    Print["The lifted full one-step annihilator order is ", fullOrder, "."];
    Print["The generic skew degrees predicted by the theorem are ",
      predictedFull, " (full) and ", predictedTwoStep, " (two-step/even subsequence)."]
  ];
  result
];
End[];
EndPackage[];
