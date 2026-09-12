(* ::Package:: *)

(* Regression tests for the Fiduccia remote-term extension.
   Intended load order: ToeplitzRecurrences.wl, optional Krylov extension,
   then ToeplitzRecurrencesFiduccia.wl. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "src", "ToeplitzRecurrencesAll.wl"}]];

ClearAll[TRFiducciaSamePolynomialQ, TRFiducciaModMatMul,
  TRFiducciaModMatPow, TRFiducciaSecondOrderMod, TRFiducciaToeplitzMatrix,
  RunFiducciaOptionTests, RunFiducciaRecurrenceTests,
  RunToeplitzFiducciaIntegrationTests, RunToeplitzFiducciaTests];

TRFiducciaSamePolynomialQ[p_, q_] := TrueQ[Expand[p - q] === 0];

TRFiducciaModMatMul[a_, b_, modulus_Integer] :=
  Mod[a.b, modulus];

TRFiducciaModMatPow[matrix_, exponent_Integer?NonNegative, modulus_Integer] :=
  Module[{base = Mod[matrix, modulus], e = exponent, result = IdentityMatrix[Length[matrix]]},
    While[e > 0,
      If[OddQ[e], result = TRFiducciaModMatMul[result, base, modulus]];
      e = Quotient[e, 2];
      If[e > 0, base = TRFiducciaModMatMul[base, base, modulus]]
    ];
    Mod[result, modulus]
  ];


TRFiducciaToeplitzMatrix[m1_Integer?NonNegative, m2_Integer?NonNegative, diagonals_List, n_Integer?NonNegative] :=
  Table[
    If[-m1 <= j - i <= m2, diagonals[[j - i + m1 + 1]], 0],
    {i, 1, n}, {j, 1, n}
  ];

RunFiducciaOptionTests[] := TestReport[{
  VerificationTest[
    Options[FiducciaPolSquarings],
    {Modulus -> None, TermCount -> 1},
    TestID -> "Fiduccia option defaults survive package initialization"
  ],

  VerificationTest[
    Options[TermForDTM],
    {
      DiagonalValues -> Automatic,
      ScalarRecurrenceMethod -> "CharacteristicPolynomial",
      Modulus -> None,
      TermCount -> 1,
      Verbose -> False,
      PositionDependent -> False
    },
    TestID -> "TermForDTM option defaults survive package initialization"
  ],

  VerificationTest[
    Options[TermForPTM],
    Options[TermForDTM],
    TestID -> "TermForPTM option defaults survive package initialization"
  ],

  VerificationTest[
    And[
      StringQ[FiducciaPolSquarings::mod],
      StringQ[FiducciaPolSquarings::count],
      StringQ[TermForDTM::method],
      StringQ[TermForPTM::method]
    ],
    True,
    TestID -> "Fiduccia and Toeplitz wrapper messages survive package initialization"
  ],

  VerificationTest[
    FiducciaPolSquarings[
      {1, 23, 3}, {111, 3332, 12}, 100,
      Modulus -> 123, TermCount -> 3
    ],
    {36, 11, 30},
    TestID -> "Fiduccia option smoke test"
  ],

  VerificationTest[
    TermForDTM[
      1, 1, 20,
      DiagonalValues -> {2, 1, 3},
      ScalarRecurrenceMethod -> "CharacteristicPolynomial",
      Modulus -> 101
    ],
    95,
    TestID -> "TermForDTM option smoke test"
  ]
}];

TRFiducciaSecondOrderMod[{c1_, c2_}, {a0_, a1_}, n_Integer?NonNegative, modulus_Integer] :=
  Which[
    n === 0, Mod[a0, modulus],
    n === 1, Mod[a1, modulus],
    True, First[
      Mod[
        TRFiducciaModMatPow[{{c1, c2}, {1, 0}}, n - 1, modulus].{a1, a0},
        modulus
      ]
    ]
  ];

RunFiducciaRecurrenceTests[] := Module[
  {p1, p2, p3, symbolic3, symbolicWindow, numericTerms, linearCheck,
   modularWindow, modularFar, modularHuge, modularHuge1000, modularHuge10000,
   modularOrder5, tests},

  symbolic3 = FiducciaPolSquarings[{1, 23, 3}, {p1, p2, p3}, 3];
  symbolicWindow = FiducciaPolSquarings[
    {1, 23, 3}, {p1, p2, p3}, 3, TermCount -> 3];

  numericTerms = Table[
    FiducciaPolSquarings[{1, 23, 3}, {4, 5, 6}, i],
    {i, 0, 6}
  ];

  linearCheck = Table[
    FiducciaPolSquarings[{2, -1, 5, 4}, {3, -2, 7, 1}, i],
    {i, 0, 40}
  ];

  modularWindow = FiducciaPolSquarings[
    {1, 23, 3}, {111, 3332, 12}, 100,
    Modulus -> 123, TermCount -> 3];

  modularFar = FiducciaPolSquarings[
    {1, 23, 3}, {111, 3332, 12}, 10^7,
    Modulus -> 123, TermCount -> 3];

  modularHuge = FiducciaPolSquarings[
    {1, 23, 3}, {111, 3332, 12}, 10^100,
    Modulus -> 123, TermCount -> 3];

  modularHuge1000 = FiducciaPolSquarings[
    {1, 23, 3}, {111, 3332, 12}, 10^1000,
    Modulus -> 123, TermCount -> 3];

  modularHuge10000 = FiducciaPolSquarings[
    {1, 23, 3}, {111, 3332, 12}, 10^10000,
    Modulus -> 123, TermCount -> 3];

  modularOrder5 = FiducciaPolSquarings[
    {1, 23, 3, -11, 133}, {4, 5, 6, -3, 6}, 10^7,
    Modulus -> 123, TermCount -> 5];

  tests = {
    VerificationTest[
      TRFiducciaSamePolynomialQ[symbolic3, 3 p1 + 23 p2 + p3],
      True,
      TestID -> "Fiduccia symbolic next term"
    ],

    VerificationTest[
      And @@ MapThread[TRFiducciaSamePolynomialQ, {
        symbolicWindow,
        {
          3 p1 + 23 p2 + p3,
          3 p1^2 + 3 p2 + 23 p1 p2 + 23 p3 + p1 p3,
          3 p1^3 + 6 p1 p2 + 23 p1^2 p2 + 23 p2^2 + 3 p3 +
            23 p1 p3 + p1^2 p3 + p2 p3
        }
      }],
      True,
      TestID -> "Fiduccia symbolic three-term window"
    ],

    VerificationTest[
      numericTerms,
      {1, 23, 3, 133, 685, 3423, 17915},
      TestID -> "Fiduccia supplied numeric sequence"
    ],

    VerificationTest[
      linearCheck,
      LinearRecurrence[{3, -2, 7, 1}, {2, -1, 5, 4}, 41],
      TestID -> "Fiduccia agrees with LinearRecurrence at moderate indices"
    ],

    VerificationTest[
      modularWindow,
      {36, 11, 30},
      TestID -> "Fiduccia supplied modular N=100 window"
    ],

    VerificationTest[
      FiducciaPolSquarings[
        {1, 23, 3}, {111, 3332, 12}, 100, Modulus -> 123],
      36,
      TestID -> "Fiduccia supplied modular N=100 scalar"
    ],

    VerificationTest[
      modularFar,
      {84, 62, 93},
      TestID -> "Fiduccia supplied modular N=10^7 window"
    ],

    VerificationTest[
      modularHuge,
      {117, 20, 117},
      TestID -> "Fiduccia supplied modular N=10^100 window"
    ],

    VerificationTest[
      modularHuge1000,
      {24, 5, 33},
      TestID -> "Fiduccia supplied modular N=10^1000 window"
    ],

    VerificationTest[
      modularHuge10000,
      {105, 56, 42},
      TestID -> "Fiduccia supplied modular N=10^10000 window"
    ],

    VerificationTest[
      modularOrder5,
      {100, 57, 68, 107, 66},
      TestID -> "Fiduccia supplied modular order-five window"
    ],

    VerificationTest[
      FiducciaPolSquarings[{1, 2}, {3}, 10],
      $Failed,
      {FiducciaPolSquarings::data},
      TestID -> "Fiduccia rejects mismatched data lengths"
    ],

    VerificationTest[
      FiducciaPolSquarings[{1, 2}, {3, 4}, -1],
      $Failed,
      {FiducciaPolSquarings::index},
      TestID -> "Fiduccia rejects negative index"
    ],

    VerificationTest[
      FiducciaPolSquarings[{1, 2}, {3, 4}, 10, Modulus -> 1],
      $Failed,
      {FiducciaPolSquarings::mod},
      TestID -> "Fiduccia rejects invalid modulus"
    ],

    VerificationTest[
      FiducciaPolSquarings[{1, 2}, {3, 4}, 10, TermCount -> 0],
      $Failed,
      {FiducciaPolSquarings::count},
      TestID -> "Fiduccia rejects zero TermCount"
    ]
  };

  TestReport[tests]
];

RunToeplitzFiducciaIntegrationTests[] := Module[
  {d50, p50, dmod, pmod, dwin, dk, d22, p21, directD22, directP21, tests},

  d50 = TermForDTM[1, 1, 50,
    DiagonalValues -> {2, 1, 3},
    ScalarRecurrenceMethod -> "CharacteristicPolynomial"];
  p50 = TermForPTM[1, 1, 50,
    DiagonalValues -> {2, 1, 3},
    ScalarRecurrenceMethod -> "CharacteristicPolynomial"];

  dmod = TermForDTM[1, 1, 10^7,
    DiagonalValues -> {2, 1, 3},
    ScalarRecurrenceMethod -> "CharacteristicPolynomial",
    Modulus -> 101];
  pmod = TermForPTM[1, 1, 10^7,
    DiagonalValues -> {2, 1, 3},
    ScalarRecurrenceMethod -> "CharacteristicPolynomial",
    Modulus -> 101];

  dwin = TermForDTM[1, 1, 30,
    DiagonalValues -> {2, 1, 3},
    TermCount -> 3];

  d22 = TermForDTM[2, 2, 8,
    DiagonalValues -> {2, -1, 3, 4, 1},
    ScalarRecurrenceMethod -> "CharacteristicPolynomial"];
  directD22 = Det[TRFiducciaToeplitzMatrix[2, 2, {2, -1, 3, 4, 1}, 8]];

  p21 = TermForPTM[2, 1, 6,
    DiagonalValues -> {1, 2, 3, 4},
    ScalarRecurrenceMethod -> "CharacteristicPolynomial"];
  directP21 = Total[
    Times @@ MapIndexed[
      TRFiducciaToeplitzMatrix[2, 1, {1, 2, 3, 4}, 6][[First[#2], #1]] &,
      #
    ] & /@ Permutations[Range[6]]
  ];

  dk = If[MemberQ[First /@ Options[LinRecForDTM], ScalarRecurrenceMethod],
    TermForDTM[1, 1, 50,
      DiagonalValues -> {2, 1, 3},
      ScalarRecurrenceMethod -> "Krylov"],
    Missing["KrylovNotLoaded"]
  ];

  tests = {
    VerificationTest[
      d50,
      Last[LinearRecurrence[{1, -6}, {1, 1}, 51]],
      TestID -> "TermForDTM tridiagonal remote determinant"
    ],

    VerificationTest[
      p50,
      Last[LinearRecurrence[{1, 6}, {1, 1}, 51]],
      TestID -> "TermForPTM tridiagonal remote permanent"
    ],

    VerificationTest[
      dmod,
      TRFiducciaSecondOrderMod[{1, -6}, {1, 1}, 10^7, 101],
      TestID -> "TermForDTM modular remote determinant"
    ],

    VerificationTest[
      pmod,
      TRFiducciaSecondOrderMod[{1, 6}, {1, 1}, 10^7, 101],
      TestID -> "TermForPTM modular remote permanent"
    ],

    VerificationTest[
      dwin,
      Take[LinearRecurrence[{1, -6}, {1, 1}, 33], -3],
      TestID -> "TermForDTM three-term window"
    ],

    VerificationTest[
      d22,
      directD22,
      TestID -> "TermForDTM pentadiagonal agrees with direct determinant"
    ],

    VerificationTest[
      p21,
      directP21,
      TestID -> "TermForPTM asymmetric band agrees with direct permanent"
    ],

    VerificationTest[
      If[MissingQ[dk], True, dk === d50],
      True,
      TestID -> "Krylov and characteristic remote determinant agree"
    ],

    VerificationTest[
      TermForDTM[1, 1, 10,
        DiagonalValues -> {2, 1, 3}, PositionDependent -> True],
      $Failed,
      {TermForDTM::pos},
      TestID -> "TermForDTM rejects position-dependent input"
    ]
  };

  TestReport[tests]
];

RunToeplitzFiducciaTests[] := <|
  "Options" -> RunFiducciaOptionTests[],
  "Recurrence" -> RunFiducciaRecurrenceTests[],
  "ToeplitzIntegration" -> RunToeplitzFiducciaIntegrationTests[]
|>;
