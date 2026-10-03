(* ::Package:: *)

(* Exact regression checks for the fixed-size Toeplitz-circulant bridge. *)

ClearAll[TRBridgeToeplitzMatrix, TRBridgeCirculant, RunToeplitzCirculantBridgeTests];

TRBridgeToeplitzMatrix[m_Integer?NonNegative, values_List, n_Integer?Positive] :=
  Table[
    With[{d = j - i}, If[-m <= d <= m, values[[d + m + 1]], 0]],
    {i, n}, {j, n}
  ];

TRBridgeCirculant[m_Integer?NonNegative, values_List, n_Integer?Positive] := Module[
  {N = n + m},
  Table[
    With[{d = Mod[j - i, N]},
      Which[
        d <= m, values[[d + m + 1]],
        d >= N - m, values[[d - N + m + 1]],
        True, 0
      ]
    ],
    {i, N}, {j, N}
  ]
];

RunToeplitzCirculantBridgeTests[] := Module[
  {cases, tests},
  cases = {
    {1, 4, {-2, 7, 3}},
    {2, 5, {2, 3, 11, 5, 7}},
    {3, 6, {2, -1, 3, 13, 5, 7, -2}}
  };
  tests = Map[
    Function[case,
      With[
        {m = case[[1]], n = case[[2]], values = case[[3]]},
        With[
          {T = TRBridgeToeplitzMatrix[case[[1]], case[[3]], case[[2]]],
           C = TRBridgeCirculant[case[[1]], case[[3]], case[[2]]],
           N = case[[2]] + case[[1]]},
          VerificationTest[
            And[
              C[[1 ;; n, 1 ;; n]] === T,
              Det[C] =!= 0,
              Together[Det[T] - Det[C] Det[Inverse[C][[n + 1 ;; N, n + 1 ;; N]]]] === 0
            ],
            True,
            TestID -> StringJoin["circulant bridge m=", ToString[m], ", n=", ToString[n]]
          ]
        ]
      ]
    ],
    cases
  ];
  TestReport[tests]
];
