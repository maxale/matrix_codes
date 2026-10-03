(* ::Package:: *)

(* Regression checks for the base implementation accompanying
   "Constructive recurrences for determinants and permanents of banded Toeplitz matrices".
   Load this file in Mathematica and evaluate RunToeplitzRecurrenceTests[]. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "src", "ToeplitzRecurrences.wl"}]];

ClearAll[
  TRSamePolynomialQ, TRToeplitzMatrix, TRPermanentForTest,
  TRSequenceForTest, TRPolynomialAnnihilatesQ, RunToeplitzRecurrenceTests
];

TRSamePolynomialQ[p_, q_, z_] := TrueQ[Expand[p - q] === 0];

TRToeplitzMatrix[n_Integer?NonNegative, m1_, m2_, values_List] := Module[{assoc},
  assoc = AssociationThread[Range[-m1, m2] -> values];
  Table[Lookup[assoc, j - i, 0], {i, n}, {j, n}]
];

TRPermanentForTest[m_List] := Module[{n = Length[m], perms},
  If[n === 0, Return[1]];
  perms = Permutations[Range[n]];
  Total[Table[Product[m[[i, p[[i]]]], {i, n}], {p, perms}]]
];

TRSequenceForTest[quantity_, m1_, m2_, values_List, nmax_Integer?NonNegative] :=
  Table[
    If[n === 0, 1,
      With[{mat = TRToeplitzMatrix[n, m1, m2, values]},
        If[quantity === "Determinant", Det[mat], TRPermanentForTest[mat]]
      ]
    ],
    {n, 0, nmax}
  ];

TRPolynomialAnnihilatesQ[poly_, z_, seq_List] := Module[{c, d},
  c = CoefficientList[Expand[poly], z];
  d = Length[c] - 1;
  If[Length[seq] < d + 1, Return[False]];
  And @@ Table[
    TrueQ[Expand[Sum[c[[j + 1]] seq[[n + j + 1]], {j, 0, d}]] === 0],
    {n, 0, Length[seq] - d - 1}
  ]
];

RunToeplitzRecurrenceTests[] := Module[
  {x = \[FormalX], k = \[FormalK], aa = \[FormalA],
   u, v, w, q, b,
   det11, det22, det32, per22, per23, generic11, custom11, sparse22,
   oneParam22, pos22, posP22, customPos22,
   incDForward, incDCentered, incPForward, incPCentered,
   incD21Forward, incD21Centered, incP21Forward, incP21Centered,
   incTranspose, scopeDet, scopePer, scopePos, scopeInc, scopeGlobals,
   expectedStates22, expectedDet11, expectedDet22, expectedPer22,
   expectedPos22, expectedPosP22, sweet22, tests},

  det11 = LinRecForDTM[1, 1, ScalarRecurrence -> True, Verbose -> False];
  det22 = LinRecForDTM[2, 2, ScalarRecurrence -> True, Verbose -> False];
  det32 = LinRecForDTM[3, 2, ScalarRecurrence -> False, Verbose -> False];
  per22 = LinRecForPTM[2, 2, ScalarRecurrence -> True, Verbose -> False];
  per23 = LinRecForPTM[2, 3, ScalarRecurrence -> False, Verbose -> False];

  generic11 = LinRecForDTM[1, 1, ScalarRecurrence -> False, Verbose -> False];
  custom11 = LinRecForDTM[1, 1,
    DiagonalValues -> {u, v, w}, ScalarRecurrence -> False, Verbose -> False];
  sparse22 = LinRecForDTM[2, 2,
    DiagonalValues -> {0, q, 1, q, 0}, ScalarRecurrence -> False, Verbose -> False];
  oneParam22 = LinRecForDTM[2, 2,
    DiagonalValues -> {q^2, q, 1, q, q^2},
    ScalarRecurrence -> True, Verbose -> False];

  pos22 = LinRecForDTM[2, 2,
    PositionDependent -> True, ScalarRecurrence -> False, Verbose -> False];
  posP22 = LinRecForPTM[2, 2,
    PositionDependent -> True, ScalarRecurrence -> False, Verbose -> False];
  customPos22 = LinRecForDTM[2, 2,
    PositionDependent -> True,
    BandEntryFunction -> Function[{s, t}, b[s, t]],
    ScalarRecurrence -> False, Verbose -> False];

  incDForward = LinRecForDTMIncreasingRows[1, 1,
    EllRange -> "Forward", Verbose -> False];
  incDCentered = LinRecForDTMIncreasingRows[1, 1,
    EllRange -> "Centered", Verbose -> False];
  incPForward = LinRecForPTMIncreasingRows[1, 1,
    EllRange -> "Forward", Verbose -> False];
  incPCentered = LinRecForPTMIncreasingRows[1, 1,
    EllRange -> "Centered", Verbose -> False];

  incD21Forward = LinRecForDTMIncreasingRows[2, 1,
    EllRange -> "Forward", Verbose -> False];
  incD21Centered = LinRecForDTMIncreasingRows[2, 1,
    EllRange -> "Centered", Verbose -> False];
  incP21Forward = LinRecForPTMIncreasingRows[2, 1,
    EllRange -> "Forward", Verbose -> False];
  incP21Centered = LinRecForPTMIncreasingRows[2, 1,
    EllRange -> "Centered", Verbose -> False];
  incTranspose = LinRecForDTMIncreasingRows[1, 2,
    EllRange -> "Centered", Verbose -> False];

  {scopeDet, scopePer, scopePos, scopeInc, scopeGlobals} = Block[
    {Global`a = 9, Global`x = 17, Global`k = 23},
    {
      LinRecForDTM[1, 1, ScalarRecurrence -> True, Verbose -> False],
      LinRecForPTM[1, 1, ScalarRecurrence -> False, Verbose -> False],
      LinRecForDTM[1, 1, PositionDependent -> True,
        ScalarRecurrence -> False, Verbose -> False],
      LinRecForDTMIncreasingRows[1, 1, Verbose -> False],
      {Global`a, Global`x, Global`k}
    }
  ];

  expectedStates22 = {
    {{k + 1}, {k + 1}},
    {{k + 2, k + 1}, {k + 2, k}},
    {{k + 2, k + 1}, {k + 2, k - 1}},
    {{k + 2, k}, {k + 2, k + 1}},
    {{k + 2, k}, {k + 2, k}},
    {{k + 3, k + 2, k + 1}, {k + 3, k + 1, k}}
  };

  expectedDet11 = {
    {aa[0], -aa[-1]},
    {aa[1], 0}
  };

  expectedDet22 = {
    {aa[0], -aa[-1], aa[-2], 0, 0, 0},
    {aa[1], 0, 0, -aa[2], 0, 0},
    {0, aa[1], 0, 0, -aa[2], 0},
    {aa[-1], -aa[-2], 0, 0, 0, 0},
    {aa[0], 0, 0, 0, 0, -aa[-2]},
    {aa[2], 0, 0, 0, 0, 0}
  };

  expectedPer22 = {
    {aa[0], aa[-1], aa[-2], 0, 0, 0},
    {aa[1], 0, 0, aa[2], 0, 0},
    {0, aa[1], 0, 0, aa[2], 0},
    {aa[-1], aa[-2], 0, 0, 0, 0},
    {aa[0], 0, 0, 0, 0, aa[-2]},
    {aa[2], 0, 0, 0, 0, 0}
  };

  expectedPos22 = {
    {aa[0][k], -aa[-1][k - 1], aa[-2][k - 2], 0, 0, 0},
    {aa[1][k], 0, 0, -aa[2][k - 1], 0, 0},
    {0, aa[1][k], 0, 0, -aa[2][k - 1], 0},
    {aa[-1][k], -aa[-2][k - 1], 0, 0, 0, 0},
    {aa[0][k + 1], 0, 0, 0, 0, -aa[-2][k - 1]},
    {aa[2][k], 0, 0, 0, 0, 0}
  };

  expectedPosP22 = {
    {aa[0][k], aa[-1][k - 1], aa[-2][k - 2], 0, 0, 0},
    {aa[1][k], 0, 0, aa[2][k - 1], 0, 0},
    {0, aa[1][k], 0, 0, aa[2][k - 1], 0},
    {aa[-1][k], aa[-2][k - 1], 0, 0, 0, 0},
    {aa[0][k + 1], 0, 0, 0, 0, aa[-2][k - 1]},
    {aa[2][k], 0, 0, 0, 0, 0}
  };

  sweet22 =
    x^6 - aa[0] x^5
      + (aa[-1] aa[1] - aa[-2] aa[2]) x^4
      + (-aa[-2] aa[1]^2 - aa[-1]^2 aa[2]
          + 2 aa[-2] aa[0] aa[2]) x^3
      + (aa[-2] aa[-1] aa[1] aa[2] - aa[-2]^2 aa[2]^2) x^2
      - aa[-2]^2 aa[0] aa[2]^2 x
      + aa[-2]^3 aa[2]^3;

  tests = {
    VerificationTest[
      MReduce[{k, k + 1, k + 2}, {k, k + 1, k + 2}],
      {{k}, {k}},
      TestID -> "MReduce strips terminal pairs"
    ],
    VerificationTest[
      MReduce[{k + 1, k - 1}, {k + 1, k}],
      {{k + 1, k - 1}, {k + 1, k}},
      TestID -> "MReduce leaves minimal signature unchanged"
    ],
    VerificationTest[
      Normal[det11["TransferMatrix"]], expectedDet11,
      TestID -> "determinant (1,1) transfer matrix"
    ],
    VerificationTest[
      TRSamePolynomialQ[det11["CharacteristicPolynomial"],
        x^2 - aa[0] x + aa[-1] aa[1], x],
      True,
      TestID -> "determinant (1,1) characteristic polynomial"
    ],
    VerificationTest[
      det22["StateOrdering"], "Paper first-discovery order",
      TestID -> "paper state-order convention is reported"
    ],
    VerificationTest[
      det22["States"], expectedStates22,
      TestID -> "determinant (2,2) states in paper order"
    ],
    VerificationTest[
      Normal[det22["TransferMatrix"]], expectedDet22,
      TestID -> "paper pentadiagonal determinant matrix"
    ],
    VerificationTest[
      TRSamePolynomialQ[det22["CharacteristicPolynomial"], sweet22, x],
      True,
      TestID -> "pentadiagonal determinant characteristic polynomial"
    ],
    VerificationTest[
      per22["States"], expectedStates22,
      TestID -> "permanent uses the same first-discovery state graph"
    ],
    VerificationTest[
      Normal[per22["TransferMatrix"]], expectedPer22,
      TestID -> "pentadiagonal permanent matrix in paper state order"
    ],
    VerificationTest[
      det32["StateOrder"], Binomial[5, 3],
      TestID -> "asymmetric determinant exact state count (3,2)"
    ],
    VerificationTest[
      per23["StateOrder"], Binomial[5, 2],
      TestID -> "asymmetric permanent exact state count (2,3)"
    ],
    VerificationTest[
      Normal[custom11["TransferMatrix"]],
      Normal[generic11["TransferMatrix"]] /. {aa[-1] -> u, aa[0] -> v, aa[1] -> w},
      TestID -> "custom diagonal values are applied before scalarization"
    ],
    VerificationTest[
      Normal[sparse22["TransferMatrix"]],
      Normal[det22["TransferMatrix"]] /. {aa[-2] -> 0, aa[-1] -> q, aa[0] -> 1,
        aa[1] -> q, aa[2] -> 0},
      TestID -> "sparse band equals generic transfer specialization"
    ],
    VerificationTest[
      MissingQ[oneParam22["CharacteristicPolynomial"]], False,
      TestID -> "one-parameter band supports forced characteristic polynomial"
    ],
    VerificationTest[
      generic11["CharacteristicPolynomial"], Missing["NotComputed"],
      TestID -> "ScalarRecurrence False skips scalarization"
    ],
    VerificationTest[
      LinRecForDTM[1, 1, ScalarRecurrence -> Automatic,
        MaxCharacteristicOrder -> 0, Verbose -> False]["CharacteristicPolynomial"],
      Missing["NotComputed"],
      TestID -> "automatic scalarization respects order cutoff"
    ],
    VerificationTest[
      MissingQ[LinRecForDTM[1, 1, ScalarRecurrence -> True,
        MaxCharacteristicOrder -> 0, Verbose -> False]["CharacteristicPolynomial"]],
      False,
      TestID -> "forced scalarization overrides order cutoff"
    ],
    VerificationTest[
      Normal[pos22["TransferMatrixExpression"]], expectedPos22,
      TestID -> "position-dependent determinant matrix in paper state order"
    ],
    VerificationTest[
      Normal[posP22["TransferMatrixExpression"]], expectedPosP22,
      TestID -> "position-dependent permanent matrix removes Laplace signs"
    ],
    VerificationTest[
      Normal[customPos22["TransferMatrixFunction"][k]],
      Normal[pos22["TransferMatrixExpression"]] /. aa[s_][t_] :> b[s, t],
      TestID -> "custom position-dependent band-entry function"
    ],
    VerificationTest[
      pos22["CharacteristicPolynomial"], Missing["PositionDependent"],
      TestID -> "position-dependent mode does not claim a constant characteristic polynomial"
    ],
    VerificationTest[
      TRSamePolynomialQ[incDForward["AnnihilatingPolynomial"],
        x^2 - aa[0] x + aa[-1] aa[1], x],
      True,
      TestID -> "increasing-rows determinant tridiagonal recurrence"
    ],
    VerificationTest[
      TRSamePolynomialQ[incDForward["AnnihilatingPolynomial"],
        incDCentered["AnnihilatingPolynomial"], x],
      True,
      TestID -> "determinant forward and centered increasing-rows agree"
    ],
    VerificationTest[
      TRSamePolynomialQ[incPForward["AnnihilatingPolynomial"],
        x^2 - aa[0] x - aa[-1] aa[1], x],
      True,
      TestID -> "increasing-rows permanent tridiagonal recurrence"
    ],
    VerificationTest[
      TRSamePolynomialQ[incPForward["AnnihilatingPolynomial"],
        incPCentered["AnnihilatingPolynomial"], x],
      True,
      TestID -> "permanent forward and centered increasing-rows agree"
    ],
    VerificationTest[
      TRSamePolynomialQ[incD21Forward["AnnihilatingPolynomial"],
        x^3 - aa[0] x^2 + aa[-1] aa[1] x - aa[-2] aa[1]^2, x],
      True,
      TestID -> "increasing-rows determinant (2,1) recurrence"
    ],
    VerificationTest[
      TRSamePolynomialQ[incD21Forward["AnnihilatingPolynomial"],
        incD21Centered["AnnihilatingPolynomial"], x],
      True,
      TestID -> "determinant (2,1) forward and centered agree"
    ],
    VerificationTest[
      TRSamePolynomialQ[incP21Forward["AnnihilatingPolynomial"],
        x^3 - aa[0] x^2 - aa[-1] aa[1] x - aa[-2] aa[1]^2, x],
      True,
      TestID -> "increasing-rows permanent (2,1) recurrence"
    ],
    VerificationTest[
      TRSamePolynomialQ[incP21Forward["AnnihilatingPolynomial"],
        incP21Centered["AnnihilatingPolynomial"], x],
      True,
      TestID -> "permanent (2,1) forward and centered agree"
    ],
    VerificationTest[
      TRPolynomialAnnihilatesQ[
        incD21Centered["AnnihilatingPolynomial"], x,
        TRSequenceForTest["Determinant", 2, 1,
          {aa[-2], aa[-1], aa[0], aa[1]}, 6]],
      True,
      TestID -> "centered determinant recurrence annihilates direct Toeplitz determinants"
    ],
    VerificationTest[
      TRPolynomialAnnihilatesQ[
        incP21Centered["AnnihilatingPolynomial"], x,
        TRSequenceForTest["Permanent", 2, 1,
          {aa[-2], aa[-1], aa[0], aa[1]}, 5]],
      True,
      TestID -> "centered permanent recurrence annihilates direct Toeplitz permanents"
    ],
    VerificationTest[
      incDForward["CoefficientMinorOrders"], {1, 2, 3},
      TestID -> "forward increasing-rows minor orders"
    ],
    VerificationTest[
      incDCentered["CoefficientMinorOrders"], {0, 1, 2},
      TestID -> "centered increasing-rows minor orders"
    ],
    VerificationTest[
      incTranspose["Transposed"], True,
      TestID -> "increasing-rows transparently transposes when m1<m2"
    ],
    VerificationTest[
      TRSamePolynomialQ[incTranspose["AnnihilatingPolynomial"],
        x^3 - aa[0] x^2 + aa[-1] aa[1] x - aa[2] aa[-1]^2, x],
      True,
      TestID -> "transposed increasing-rows recurrence uses original diagonals"
    ],
    VerificationTest[
      Normal[scopeDet["TransferMatrix"]], expectedDet11,
      TestID -> "Global a assignment does not affect determinant defaults"
    ],
    VerificationTest[
      Normal[scopePer["TransferMatrix"]],
      {{aa[0], aa[-1]}, {aa[1], 0}},
      TestID -> "Global a assignment does not affect permanent defaults"
    ],
    VerificationTest[
      scopeDet["CharacteristicPolynomialVariable"], x,
      TestID -> "Global x assignment does not affect characteristic variable"
    ],
    VerificationTest[
      scopeDet["States"][[1]], {{k + 1}, {k + 1}},
      TestID -> "Global k assignment does not affect state parameter"
    ],
    VerificationTest[
      Normal[scopePos["TransferMatrixExpression"]][[1, 1]], aa[0][k],
      TestID -> "Global a and k assignments do not affect position-dependent defaults"
    ],
    VerificationTest[
      scopeInc["CharacteristicPolynomialVariable"], x,
      TestID -> "Global x assignment does not affect increasing-rows variable"
    ],
    VerificationTest[
      scopeGlobals, {9, 17, 23},
      TestID -> "constructors do not modify user Global values"
    ]
  };

  TestReport[tests]
];

(* Small interactive examples. *)
(* LinRecForDTM[2, 2] *)
(* LinRecForPTM[2, 2] *)
(* LinRecForDTM[3, 3, ScalarRecurrence -> False] *)
(* LinRecForDTM[3, 3, ScalarRecurrence -> True] *)
(* LinRecForDTM[2, 2, DiagonalValues -> {0, q, 1, q, 0}] *)
(* LinRecForDTM[2, 2, PositionDependent -> True] *)
(* LinRecForDTMIncreasingRows[2, 1, EllRange -> "Forward"] *)
(* LinRecForDTMIncreasingRows[2, 1, EllRange -> "Centered"] *)
