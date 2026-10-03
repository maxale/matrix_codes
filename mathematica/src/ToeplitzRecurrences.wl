(* ::Package:: *)

(*
  ToeplitzRecurrences.wl
  Wolfram Language implementation accompanying
  "Constructive recurrences for determinants and permanents of banded Toeplitz matrices".

  Public convention: m1 is the number of subdiagonals (negative offsets)
  and m2 is the number of superdiagonals (positive offsets).
  Row-column states are numbered in first-discovery order.
  Intended for Mathematica 2021.1 or later; no external packages are used.
*)

BeginPackage["ToeplitzRecurrences`"];

MReduce::usage =
  "MReduce[rows, cols] normalizes a row-column boundary-minor signature.";
LaplaceD::usage =
  "LaplaceD[m1,m2,S] performs one determinant Laplace expansion for a normalized signature S. m1 counts subdiagonals and m2 superdiagonals.";
LaplaceP::usage =
  "LaplaceP[m1,m2,S] performs the signless permanent analogue of LaplaceD.";
LinRecForDTM::usage =
  "LinRecForDTM[m1,m2,opts] constructs the row-column transfer matrix for determinants of an (m1,m2)-banded Toeplitz family.";
LinRecForPTM::usage =
  "LinRecForPTM[m1,m2,opts] constructs the row-column transfer matrix for permanents of an (m1,m2)-banded Toeplitz family.";
LinRecForDTMIncreasingRows::usage =
  "LinRecForDTMIncreasingRows[m1,m2,opts] implements Algorithm 1 (increasing rows) for determinants.";
LinRecForPTMIncreasingRows::usage =
  "LinRecForPTMIncreasingRows[m1,m2,opts] implements Algorithm 1 (increasing rows) for permanents.";

DiagonalValues::usage =
  "DiagonalValues is an option specifying constant diagonal entries in offset order -m1,...,m2. Automatic uses \[FormalA][-m1],...,\[FormalA][m2].";
ScalarRecurrence::usage =
  "ScalarRecurrence controls scalar recurrence extraction: Automatic, True, or False.";
CharacteristicPolynomialVariable::usage =
  "CharacteristicPolynomialVariable specifies the polynomial variable used for scalar recurrences.";
MaxCharacteristicOrder::usage =
  "MaxCharacteristicOrder is the transfer dimension cutoff used when ScalarRecurrence->Automatic.";
PositionDependent::usage =
  "PositionDependent->True makes row-column transition entries depend on the position along each diagonal.";
BandEntryFunction::usage =
  "BandEntryFunction specifies f[s,t] for the entry on offset s at position t=min(i,j) in position-dependent mode.";
EllRange::usage =
  "EllRange specifies the levels used by the increasing-rows method: \"Forward\", \"Centered\", or an increasing list of d+1 integers.";
Verbose::usage =
  "Verbose controls interactive printing by the recurrence constructors.";

Begin["`Private`"];


trMessages::diaglen = "`1`: DiagonalValues must contain exactly `2` entries, in offset order `3`.";
trMessages::posdiag = "`1`: DiagonalValues cannot be supplied together with PositionDependent->True; use BandEntryFunction instead.";
trMessages::scalar = "`1`: ScalarRecurrence must be Automatic, True, or False.";
trMessages::cutoff = "`1`: MaxCharacteristicOrder must be a nonnegative integer.";
trMessages::posscalar = "`1`: a constant characteristic polynomial is not defined for position-dependent weights; returning the transfer matrix only.";
trMessages::laplace = "`1`: the boundary signature `2` did not match an expansion case.";
trMessages::target = "`1`: a Laplace target signature could not be matched to a closed state: `2`.";
trMessages::collision = "`1`: two Laplace terms produced the same normalized transition from state `2` to state `3`.";
trMessages::ell = "`1`: EllRange must be \"Forward\", \"Centered\", or a strictly increasing list of `2` integers with minimum at least `3`.";
trMessages::kernel = "`1`: no nonzero left-kernel vector was found for the increasing-rows coefficient matrix.";
trMessages::minor = "`1`: internal coefficient-minor row and column sets have different sizes (`2` and `3`).";

trMessage[head_, "diaglen", n_, offsets_] := Message[trMessages::diaglen, head, n, offsets];
trMessage[head_, "posdiag"] := Message[trMessages::posdiag, head];
trMessage[head_, "scalar"] := Message[trMessages::scalar, head];
trMessage[head_, "cutoff"] := Message[trMessages::cutoff, head];
trMessage[head_, "posscalar"] := Message[trMessages::posscalar, head];
trMessage[head_, "laplace", sig_] := Message[trMessages::laplace, head, sig];
trMessage[head_, "target", sig_] := Message[trMessages::target, head, sig];
trMessage[head_, "collision", i_, j_] := Message[trMessages::collision, head, i, j];
trMessage[head_, "ell", n_, min_] := Message[trMessages::ell, head, n, min];
trMessage[head_, "kernel"] := Message[trMessages::kernel, head];
trMessage[head_, "minor", nr_, nc_] := Message[trMessages::minor, head, nr, nc];

$TRRowColumnOptions = {
  DiagonalValues -> Automatic,
  ScalarRecurrence -> Automatic,
  CharacteristicPolynomialVariable -> \[FormalX],
  MaxCharacteristicOrder -> 20,
  Verbose -> True,
  PositionDependent -> False,
  BandEntryFunction -> Automatic
};

Options[LinRecForDTM] = $TRRowColumnOptions;
Options[LinRecForPTM] = $TRRowColumnOptions;
Options[trRowColumnEngine] = $TRRowColumnOptions;
Options[LaplaceD] = {
  DiagonalValues -> Automatic,
  PositionDependent -> False,
  BandEntryFunction -> Automatic
};
Options[LaplaceP] = Options[LaplaceD];

$TRIncreasingRowsOptions = {
  DiagonalValues -> Automatic,
  EllRange -> "Forward",
  CharacteristicPolynomialVariable -> \[FormalX],
  Verbose -> True
};

Options[LinRecForDTMIncreasingRows] = $TRIncreasingRowsOptions;
Options[LinRecForPTMIncreasingRows] = $TRIncreasingRowsOptions;
Options[trIncreasingRowsEngine] = $TRIncreasingRowsOptions;

MReduce[l1_List, l2_List] := Module[{L1, L2},
  L1 = Reverse[Sort[l1]];
  L2 = Reverse[Sort[l2]];
  While[
    Length[L1] >= 2 && Length[L2] >= 2 &&
    TrueQ[Simplify[L1[[2]] == L2[[2]]]] &&
    TrueQ[Simplify[L1[[1]] == L2[[1]]]] &&
    TrueQ[Simplify[L1[[1]] == L2[[2]] + 1]],
    L1 = Rest[L1];
    L2 = Rest[L2];
  ];
  {L1, L2}
];

ClearAll[trFinalSignature];
trFinalSignature[S : {_List, _List}] := Module[{len, delta},
  len = Length[S[[1]]];
  delta = Simplify[len - (S[[1, 1]] - \[FormalK])];
  S /. \[FormalK] -> \[FormalK] + delta
];

ClearAll[trDefaultDiagonalValues, trBandEntrySetup];
trDefaultDiagonalValues[m1_, m2_] := (\[FormalA][#] & /@ Range[-m1, m2]);

trBandEntrySetup[head_, m1_, m2_, diag_, pos_, bandfun_] := Module[
  {values, assoc, entry},
  If[TrueQ[pos],
    If[diag =!= Automatic,
      trMessage[head, "posdiag"];
      Return[$Failed]
    ];
    entry = If[
      bandfun === Automatic,
      With[{aa = \[FormalA]}, Function[{s, t}, aa[s][t]]],
      bandfun
    ];
    Return[<|
      "EntryFunction" -> entry,
      "DiagonalValues" -> Missing["PositionDependent"]
    |>]
  ];

  values = If[diag === Automatic, trDefaultDiagonalValues[m1, m2], diag];
  If[!ListQ[values] || Length[values] =!= m1 + m2 + 1,
    trMessage[head, "diaglen", m1 + m2 + 1, Row[Range[-m1, m2], ","]];
    Return[$Failed]
  ];
  assoc = AssociationThread[Range[-m1, m2] -> values];
  entry = With[{aassoc = assoc}, Function[{s, t}, Lookup[aassoc, s, 0]]];
  <|"EntryFunction" -> entry, "DiagonalValues" -> values|>
];

ClearAll[trLaplaceExpansion];
trLaplaceExpansion[head_, m1_, m2_, S : {_List, _List}, quantity_, entry_] := Module[
  {expan = {}, len, Sf, j, s, row, col, coeff, target, gap},

  len = Length[S[[1]]];
  Sf = trFinalSignature[S];

  If[len === 1,
    Do[
      coeff = If[quantity === "Determinant", (-1)^j, 1] entry[-j, \[FormalK] - j];
      target = MReduce[
        Append[Sf[[1]], \[FormalK]],
        Append[Sf[[2]], \[FormalK] - j]
      ];
      expan = Append[expan, {coeff, target}],
      {j, 0, m1}
    ];
    Return[{Sf, expan}]
  ];

  gap = Simplify[Sf[[1, 1]] - Sf[[1, 2]]];

  If[TrueQ[gap === 1],
    s = 0;
    Do[
      row = Sf[[1, 1]] - 1 - j;
      col = Sf[[2, 1]] - 1;
      If[!MemberQ[Sf[[1]], row],
        coeff = If[quantity === "Determinant", (-1)^s, 1] entry[j, row];
        target = MReduce[
          Append[Sf[[1]], row],
          Append[Sf[[2]], col]
        ];
        expan = Append[expan, {coeff, target}];
        s++;
      ],
      {j, 0, m2}
    ];
    Return[{Sf, expan}]
  ];

  If[IntegerQ[gap] && gap > 1,
    s = 0;
    Do[
      row = Sf[[1, 1]] - 1;
      col = Sf[[2, 1]] - 1 - j;
      If[!MemberQ[Sf[[2]], col],
        coeff = If[quantity === "Determinant", (-1)^s, 1] entry[-j, col];
        target = MReduce[
          Append[Sf[[1]], row],
          Append[Sf[[2]], col]
        ];
        expan = Append[expan, {coeff, target}];
        s++;
      ],
      {j, 0, m1}
    ];
    Return[{Sf, expan}]
  ];

  trMessage[head, "laplace", S];
  $Failed
];

LaplaceD[m1_Integer?NonNegative, m2_Integer?NonNegative,
    S : {_List, _List}, OptionsPattern[]] := Module[{setup},
  setup = trBandEntrySetup[LaplaceD, m1, m2,
    OptionValue[DiagonalValues], OptionValue[PositionDependent],
    OptionValue[BandEntryFunction]];
  If[setup === $Failed, Return[$Failed]];
  trLaplaceExpansion[LaplaceD, m1, m2, S, "Determinant", setup["EntryFunction"]]
];

LaplaceP[m1_Integer?NonNegative, m2_Integer?NonNegative,
    S : {_List, _List}, OptionsPattern[]] := Module[{setup},
  setup = trBandEntrySetup[LaplaceP, m1, m2,
    OptionValue[DiagonalValues], OptionValue[PositionDependent],
    OptionValue[BandEntryFunction]];
  If[setup === $Failed, Return[$Failed]];
  trLaplaceExpansion[LaplaceP, m1, m2, S, "Permanent", setup["EntryFunction"]]
];

ClearAll[trRowColumnClosure];
trRowColumnClosure[head_, m1_, m2_, quantity_, entry_] := Module[
  {seen, ans = {}, data, tp, stateIndex = 1,
   maxStates, finals, shiftedFinals, ord, rules = {}, pos, i, j, term, key,
   occupied = <||>, matrix},

  seen = {{{\[FormalK] + 1}, {\[FormalK] + 1}}};
  maxStates = 2 Binomial[m1 + m2, m1] + m1 + m2 + 10;

  While[stateIndex <= Length[seen],
    If[stateIndex > maxStates,
      Return[Failure["StateClosureLimit", <|
        "Message" -> "Row-column state closure exceeded the safety state bound.",
        "SeenStates" -> seen
      |>]]
    ];

    data = trLaplaceExpansion[head, m1, m2, seen[[stateIndex]], quantity, entry];
    If[data === $Failed, Return[$Failed]];
    ans = Append[ans, data];

    Do[
      tp = trFinalSignature[data[[2, j, 2]]];
      tp = MReduce[tp[[1]], tp[[2]]];
      If[!MemberQ[seen, tp], seen = Append[seen, tp]],
      {j, Length[data[[2]]]}
    ];

    stateIndex++;
  ];

  finals = ans[[All, 1]];
  shiftedFinals = (# /. \[FormalK] -> \[FormalK] - 1) & /@ finals;
  ord = Length[ans];

  Do[
    Do[
      term = ans[[i, 2, j]];
      pos = FirstPosition[shiftedFinals, term[[2]], Missing["NotFound"]];
      If[MissingQ[pos],
        trMessage[head, "target", term[[2]]];
        Return[$Failed]
      ];
      key = ToString[{i, First[pos]}, InputForm];
      If[KeyExistsQ[occupied, key],
        trMessage[head, "collision", i, First[pos]];
        Return[$Failed]
      ];
      AssociateTo[occupied, key -> True];
      rules = Append[rules, {i, First[pos]} -> term[[1]]],
      {j, Length[ans[[i, 2]]]}
    ],
    {i, ord}
  ];

  matrix = SparseArray[rules, {ord, ord}];
  <|
    "ExpansionData" -> ans,
    "States" -> finals,
    "ShiftedStates" -> shiftedFinals,
    "StateOrder" -> ord,
    "StateOrdering" -> "Paper first-discovery order",
    "TransferMatrixExpression" -> matrix
  |>
];

ClearAll[trShouldScalarize];
trShouldScalarize[head_, policy_, ord_, cutoff_] := Module[{},
  If[!MemberQ[{Automatic, True, False}, policy],
    trMessage[head, "scalar"];
    Return[$Failed]
  ];
  If[policy === Automatic && (!IntegerQ[cutoff] || cutoff < 0),
    trMessage[head, "cutoff"];
    Return[$Failed]
  ];
  Switch[policy,
    True, True,
    False, False,
    Automatic, ord <= cutoff
  ]
];

trRowColumnEngine[head_, quantity_, m1_, m2_, opts : OptionsPattern[]] := Module[
  {diag, policy, var, cutoff, verbose, posdep, bandfun, setup, closure,
   ord, matrixExpr, transfer, qfun, should, charpoly, values, result},

  diag = OptionValue[DiagonalValues];
  policy = OptionValue[ScalarRecurrence];
  var = OptionValue[CharacteristicPolynomialVariable];
  cutoff = OptionValue[MaxCharacteristicOrder];
  verbose = TrueQ[OptionValue[Verbose]];
  posdep = TrueQ[OptionValue[PositionDependent]];
  bandfun = OptionValue[BandEntryFunction];

  setup = trBandEntrySetup[head, m1, m2, diag, posdep, bandfun];
  If[setup === $Failed, Return[$Failed]];
  values = setup["DiagonalValues"];

  closure = trRowColumnClosure[head, m1, m2, quantity, setup["EntryFunction"]];
  If[closure === $Failed || FailureQ[closure], Return[closure]];

  ord = closure["StateOrder"];
  matrixExpr = closure["TransferMatrixExpression"];

  If[!MemberQ[{Automatic, True, False}, policy],
    trMessage[head, "scalar"];
    Return[$Failed]
  ];

  If[posdep,
    If[policy === True, trMessage[head, "posscalar"]];
    charpoly = Missing["PositionDependent"];
    qfun = Function[{z}, Evaluate[matrixExpr /. \[FormalK] -> z]];
    transfer = qfun,

    should = trShouldScalarize[head, policy, ord, cutoff];
    If[should === $Failed, Return[$Failed]];
    charpoly = If[
      TrueQ[should],
      Collect[CharacteristicPolynomial[matrixExpr, var], var],
      Missing["NotComputed"]
    ];
    qfun = Missing["ConstantMatrix"];
    transfer = matrixExpr
  ];

  result = <|
    "Method" -> "RowColumn",
    "Quantity" -> quantity,
    "m1" -> m1,
    "m2" -> m2,
    "OffsetConvention" -> "Paper: offsets -m1,...,m2",
    "States" -> closure["States"],
    "ShiftedStates" -> closure["ShiftedStates"],
    "StateOrder" -> ord,
    "StateOrdering" -> closure["StateOrdering"],
    "ScalarRecurrenceMethod" -> "CharacteristicPolynomial",
    "TransferMatrix" -> transfer,
    "TransferMatrixExpression" -> matrixExpr,
    "TransferMatrixFunction" -> qfun,
    "DiagonalValues" -> values,
    "PositionDependent" -> posdep,
    "CharacteristicPolynomialVariable" -> var,
    "CharacteristicPolynomial" -> charpoly
  |>;

  If[verbose,
    Print["The transfer-matrix order is ", ord, "."];
    If[posdep,
      Print["The position-dependent transfer matrix Q[k] is:"];
      Print[MatrixForm[Normal[matrixExpr]]],
      Print["The transfer matrix is:"];
      Print[MatrixForm[Normal[matrixExpr]]];
      If[MissingQ[charpoly],
        Print["The characteristic polynomial was not computed."],
        Print["The characteristic polynomial is:"];
        Print[charpoly]
      ]
    ]
  ];

  result
];

LinRecForDTM[m1_Integer?NonNegative, m2_Integer?NonNegative,
    OptionsPattern[]] :=
  trRowColumnEngine[LinRecForDTM, "Determinant", m1, m2,
    DiagonalValues -> OptionValue[DiagonalValues],
    ScalarRecurrence -> OptionValue[ScalarRecurrence],
    CharacteristicPolynomialVariable -> OptionValue[CharacteristicPolynomialVariable],
    MaxCharacteristicOrder -> OptionValue[MaxCharacteristicOrder],
    Verbose -> OptionValue[Verbose],
    PositionDependent -> OptionValue[PositionDependent],
    BandEntryFunction -> OptionValue[BandEntryFunction]
  ];

LinRecForPTM[m1_Integer?NonNegative, m2_Integer?NonNegative,
    OptionsPattern[]] :=
  trRowColumnEngine[LinRecForPTM, "Permanent", m1, m2,
    DiagonalValues -> OptionValue[DiagonalValues],
    ScalarRecurrence -> OptionValue[ScalarRecurrence],
    CharacteristicPolynomialVariable -> OptionValue[CharacteristicPolynomialVariable],
    MaxCharacteristicOrder -> OptionValue[MaxCharacteristicOrder],
    Verbose -> OptionValue[Verbose],
    PositionDependent -> OptionValue[PositionDependent],
    BandEntryFunction -> OptionValue[BandEntryFunction]
  ];

ClearAll[trPermanentRyser];
trPermanentRyser[m_List] := Module[{n, subsets},
  n = Length[m];
  If[n === 0, Return[1]];
  subsets = Rest[Subsets[Range[n]]];
  (-1)^n Total[
    Table[
      (-1)^Length[S] Times @@ Table[Total[m[[i, S]]], {i, n}],
      {S, subsets}
    ]
  ]
];

ClearAll[trMinorValue];
trMinorValue[head_, rows_List, cols_List, quantity_, entry_] := Module[{mat},
  If[Length[rows] =!= Length[cols],
    trMessage[head, "minor", Length[rows], Length[cols]];
    Return[$Failed]
  ];
  If[rows === {}, Return[1]];
  mat = Table[entry[cols[[j]] - rows[[i]], 0],
    {i, Length[rows]}, {j, Length[cols]}];
  If[quantity === "Determinant", Det[mat], trPermanentRyser[mat]]
];

ClearAll[trBoundaryStates];
trBoundaryStates[m1_, m2_] := Reverse /@ Subsets[
  Range[m1 - m2, 2 m1 - 1], {m2}
];

ClearAll[trIncreasingCoefficient];
trIncreasingCoefficient[head_, m1_, m2_, ell_, S_, quantity_, entry_, k0_] := Module[
  {Rell, CellS, H, J, initial, tail},

  If[ell >= 0,
    Rell = Range[k0 + 1 - m1, k0 + ell];
    initial = (k0 - #) & /@ S;
    tail = Range[k0 + 1 + m2 - m1, k0 + ell];
    CellS = Union[initial, tail];
    Return[trMinorValue[head, Rell, CellS, quantity, entry]]
  ];

  H = (k0 - #) & /@ S;
  J = Complement[Range[1, k0 + m2 - m1], H];
  If[Complement[J, Range[1, k0 + ell]] =!= {}, Return[0]];
  Rell = Range[k0 + 1 - m1, k0 + ell];
  CellS = Complement[Range[1, k0 + ell], J];
  trMinorValue[head, Rell, CellS, quantity, entry]
];

ClearAll[trNormalizeKernelVector];
trNormalizeKernelVector[v_List] := Module[{nz, lead},
  nz = Select[Range[Length[v]], v[[#]] =!= 0 &];
  If[nz === {}, Return[v]];
  lead = v[[Last[nz]]];
  (Cancel[Together[#/lead]] &) /@ v
];

trIncreasingRowsEngine[head_, quantity_, m1in_, m2in_, opts : OptionsPattern[]] := Module[
  {diag, ellOpt, var, verbose, m1 = m1in, m2 = m2in, transposed = False,
   setup, values, entry, originalAssoc, internalValues, internalAssoc,
   states, d, ellValues, k0, B, basis, c, poly, orders, result},

  diag = OptionValue[DiagonalValues];
  ellOpt = OptionValue[EllRange];
  var = OptionValue[CharacteristicPolynomialVariable];
  verbose = TrueQ[OptionValue[Verbose]];

  setup = trBandEntrySetup[head, m1in, m2in, diag, False, Automatic];
  If[setup === $Failed, Return[$Failed]];
  values = setup["DiagonalValues"];

  If[m1 < m2,
    transposed = True;
    originalAssoc = AssociationThread[Range[-m1, m2] -> values];
    {m1, m2} = {m2, m1};
    internalValues = Table[Lookup[originalAssoc, -s, 0], {s, -m1, m2}],
    internalValues = values
  ];

  internalAssoc = AssociationThread[Range[-m1, m2] -> internalValues];
  entry = With[{aassoc = internalAssoc}, Function[{s, t}, Lookup[aassoc, s, 0]]];

  d = Binomial[m1 + m2, m1];
  states = trBoundaryStates[m1, m2];

  ellValues = Which[
    ellOpt === "Forward", Range[0, d],
    ellOpt === "Centered", Range[-m1, d - m1],
    ListQ[ellOpt] && Length[ellOpt] === d + 1 &&
      And @@ (IntegerQ /@ ellOpt) && DuplicateFreeQ[ellOpt] && OrderedQ[ellOpt] &&
      Min[ellOpt] >= -m1,
      ellOpt,
    True,
      trMessage[head, "ell", d + 1, -m1];
      Return[$Failed]
  ];

  k0 = 4 (m1 + m2 + d + 2);
  B = Table[
    trIncreasingCoefficient[head, m1, m2, ellValues[[i]], states[[j]],
      quantity, entry, k0],
    {i, Length[ellValues]}, {j, Length[states]}
  ];
  If[!FreeQ[B, $Failed], Return[$Failed]];

  basis = NullSpace[Transpose[B]];
  If[basis === {},
    trMessage[head, "kernel"];
    Return[$Failed]
  ];
  c = trNormalizeKernelVector[First[basis]];
  poly = Expand[Sum[
    c[[i]] var^(ellValues[[i]] - Min[ellValues]),
    {i, Length[ellValues]}
  ]];
  orders = m1 + ellValues;

  result = <|
    "Method" -> "IncreasingRows",
    "Quantity" -> quantity,
    "m1" -> m1in,
    "m2" -> m2in,
    "Internalm1" -> m1,
    "Internalm2" -> m2,
    "Transposed" -> transposed,
    "OffsetConvention" -> "Paper: offsets -m1,...,m2",
    "DiagonalValues" -> values,
    "InternalDiagonalValues" -> internalValues,
    "BoundaryStates" -> states,
    "EllValues" -> ellValues,
    "CoefficientMinorOrders" -> orders,
    "CoefficientMatrix" -> B,
    "KernelVector" -> c,
    "CharacteristicPolynomialVariable" -> var,
    "AnnihilatingPolynomial" -> poly
  |>;

  If[verbose,
    Print["The increasing-rows boundary-state count is ", d, "."];
    Print["The ell values are ", ellValues, "."];
    Print["The coefficient-minor orders are ", orders, "."];
    Print["An annihilating polynomial is:"];
    Print[poly]
  ];

  result
];

LinRecForDTMIncreasingRows[m1_Integer?NonNegative, m2_Integer?NonNegative,
    OptionsPattern[]] :=
  trIncreasingRowsEngine[LinRecForDTMIncreasingRows, "Determinant", m1, m2,
    DiagonalValues -> OptionValue[DiagonalValues],
    EllRange -> OptionValue[EllRange],
    CharacteristicPolynomialVariable -> OptionValue[CharacteristicPolynomialVariable],
    Verbose -> OptionValue[Verbose]
  ];

LinRecForPTMIncreasingRows[m1_Integer?NonNegative, m2_Integer?NonNegative,
    OptionsPattern[]] :=
  trIncreasingRowsEngine[LinRecForPTMIncreasingRows, "Permanent", m1, m2,
    DiagonalValues -> OptionValue[DiagonalValues],
    EllRange -> OptionValue[EllRange],
    CharacteristicPolynomialVariable -> OptionValue[CharacteristicPolynomialVariable],
    Verbose -> OptionValue[Verbose]
  ];

End[];
EndPackage[];
