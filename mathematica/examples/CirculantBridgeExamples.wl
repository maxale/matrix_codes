(* ::Package:: *)

(* Exact paper-facing example for the fixed-size Toeplitz-circulant correction. *)

ClearAll[RunToeplitzCirculantBridgeExample];
RunToeplitzCirculantBridgeExample[] := Module[
  {m = 2, n = 5, values = {2, 3, 11, 5, 7}, N, T, C, G},
  N = n + m;
  T = Table[
    With[{d = j - i}, If[-m <= d <= m, values[[d + m + 1]], 0]],
    {i, n}, {j, n}
  ];
  C = Table[
    With[{d = Mod[j - i, N]},
      Which[
        d <= m, values[[d + m + 1]],
        d >= N - m, values[[d - N + m + 1]],
        True, 0
      ]
    ],
    {i, N}, {j, N}
  ];
  G = Inverse[C][[n + 1 ;; N, n + 1 ;; N]];
  <|
    "ToeplitzDeterminant" -> Det[T],
    "CirculantDeterminant" -> Det[C],
    "CorrectionDeterminant" -> Det[G],
    "BridgeVerified" -> TrueQ[Together[Det[T] - Det[C] Det[G]] === 0]
  |>
];
