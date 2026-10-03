(* ::Package:: *)

(* Remote-term evaluation for constant-coefficient recurrences, accompanying
   "Constructive recurrences for determinants and permanents of banded Toeplitz matrices".

   The core Fiduccia routine is a programming cleanup of the co-author's
   implementation with hardcoded modular-polynomial squarings from
   D. I. Khomovsky, "Efficient Computation of Terms of Linear Recurrence
   Sequences of Any Order", INTEGERS 18 (2018), A39.  The mathematical
   identities and binary transition logic are preserved; exact and modular
   arithmetic are handled by one code path.

   Load after ToeplitzRecurrences.wl.  For
   ScalarRecurrenceMethod -> "Krylov", also load
   ToeplitzRecurrencesKrylov.wl first. *)

BeginPackage["ToeplitzRecurrences`"];

(* Clear public function symbols before assigning usage strings, messages,
   options, and definitions.  A previous initialization-order bug showed that
   clearing these symbols later erases their Options and message texts. *)
ClearAll[FiducciaPolSquarings];
ClearAll[TermForDTM];
ClearAll[TermForPTM];

FiducciaPolSquarings::usage =
  "FiducciaPolSquarings[initial, recurrence, N] computes the N-th term (zero-based) of a constant-coefficient linear recurrence using Fiduccia's binary method with hardcoded polynomial-squaring reductions. Use Modulus -> m for integer arithmetic modulo m and TermCount -> r to return r consecutive terms starting at N.";

TermForDTM::usage =
  "TermForDTM[m1, m2, N] constructs a constant scalar recurrence for the determinant sequence of the (m1,m2)-banded Toeplitz family and evaluates its N-th term with FiducciaPolSquarings. Here m1 counts subdiagonals and m2 counts superdiagonals.";

TermForPTM::usage =
  "TermForPTM[m1, m2, N] constructs a constant scalar recurrence for the permanent sequence of the (m1,m2)-banded Toeplitz family and evaluates its N-th term with FiducciaPolSquarings. Here m1 counts subdiagonals and m2 counts superdiagonals.";

TermCount::usage =
  "TermCount is an option for remote-term routines. TermCount -> r returns r consecutive terms starting at the requested zero-based index; the default is 1.";

(* Keep the same option name as the optional Krylov extension.  Repeating the
   usage assignment is harmless when that extension has already been loaded. *)
ScalarRecurrenceMethod::usage =
  "ScalarRecurrenceMethod selects the constant-transfer scalarization method: \"CharacteristicPolynomial\" (default) or \"Krylov\".";

Begin["`Private`"];

FiducciaPolSquarings::data =
  "The initial-value list and recurrence-coefficient list must have the same length, at least 2.";
FiducciaPolSquarings::index =
  "The term index `1` must be a nonnegative integer.";
FiducciaPolSquarings::mod =
  "Modulus must be None or an integer greater than 1.";
FiducciaPolSquarings::moddata =
  "With a finite Modulus, all initial values and recurrence coefficients must be integers.";
FiducciaPolSquarings::count =
  "TermCount must be a positive integer.";

Options[FiducciaPolSquarings] = {
  Modulus -> None,
  TermCount -> 1
};

FiducciaPolSquarings[
    initial_List, recurrence_List, index_, OptionsPattern[]] := Module[
  {
    n = Length[initial], modulus = OptionValue[Modulus],
    termCount = OptionValue[TermCount], normalize, normalizeInput,
    e, bits, X, x, f, u, p = recurrence, a = initial,
    halfFloor, halfCeiling, root, weights, termAt,
    bitPosition, i, j, jt, l, r, sequence, requested
  },

  If[n =!= Length[recurrence] || n < 2,
    Message[FiducciaPolSquarings::data];
    Return[$Failed]
  ];

  If[!IntegerQ[index] || index < 0,
    Message[FiducciaPolSquarings::index, index];
    Return[$Failed]
  ];

  If[!(modulus === None || (IntegerQ[modulus] && modulus > 1)),
    Message[FiducciaPolSquarings::mod];
    Return[$Failed]
  ];

  If[!IntegerQ[termCount] || termCount < 1,
    Message[FiducciaPolSquarings::count];
    Return[$Failed]
  ];

  If[modulus =!= None && !AllTrue[Join[initial, recurrence], IntegerQ],
    Message[FiducciaPolSquarings::moddata];
    Return[$Failed]
  ];

  normalize = If[
    modulus === None,
    Function[value, Expand[value]],
    Function[value, Mod[value, modulus]]
  ];
  normalizeInput = If[
    modulus === None,
    Identity,
    Function[value, Mod[value, modulus]]
  ];

  (* If the requested window begins inside the supplied initial block, a
     short sequential extension is simpler and avoids entering the binary
     jump before it is needed. *)
  If[index < n,
    sequence = normalizeInput /@ initial;
    While[Length[sequence] < index + termCount,
      sequence = Append[
        sequence,
        normalize[
          Sum[
            p[[i]] sequence[[Length[sequence] + 1 - i]],
            {i, 1, n}
          ]
        ]
      ]
    ];
    requested = Take[sequence, {index + 1, index + termCount}];
    Return[If[termCount === 1, First[requested], requested]]
  ];

  e = Mod[n, 2];
  bits = IntegerDigits[index, 2];
  halfFloor = Floor[n/2];
  halfCeiling = Ceiling[n/2];

  (* Initial polynomial residues from the co-author's implementation. *)
  Do[X[i] = 0, {i, 1, n - 2}];
  X[n - 1] = 1;
  X[n] = p[[1]];

  If[n === 2,
    (* Hardcoded order-two squaring identities. *)
    f[1] = u[1] (-p[[1]] u[1] + 2 u[2]);
    f[2] = p[[2]] u[1]^2 + u[2]^2;

    (* Preserve the original difference-of-squares optimization when
       -p2 is a perfect square integer. *)
    root = Sqrt[-p[[2]]];
    If[IntegerQ[root],
      f[2] = (u[2] - root u[1]) (u[2] + root u[1])
    ];

    f[3] = u[2] (p[[1]] u[2] + 2 p[[2]] u[1]),

    (* Extend the symbolic u-sequence just far enough for the hardcoded
       even/odd polynomial-squaring formulas below. *)
    Do[
      u[n + i] = If[
        modulus === None,
        Sum[p[[j + 1]] u[n + i - 1 - j], {j, 0, n - 1}],
        Mod[
          Sum[p[[j + 1]] u[n + i - 1 - j], {j, 0, n - 1}],
          modulus
        ]
      ],
      {i, 1, halfFloor}
    ];

    (* X_(2 k). *)
    f[1] = If[
      modulus === None,
      e u[(n + 1)/2]^2 +
        Sum[
          u[halfFloor - i] * Expand[
            2 u[halfCeiling + 1 + i] -
              p[[2 i + e + 1]] u[halfFloor - i] -
              2 Sum[
                p[[j + 1]] u[halfCeiling + i - j],
                {j, 0, 2 i - 1 + e}
              ]
          ],
          {i, 0, halfFloor - 1}
        ],
      Mod[
        e u[(n + 1)/2]^2 +
          Sum[
            u[halfFloor - i] * Expand[
              2 u[halfCeiling + 1 + i] -
                p[[2 i + e + 1]] u[halfFloor - i] -
                2 Sum[
                  p[[j + 1]] u[halfCeiling + i - j],
                  {j, 0, 2 i - 1 + e}
                ]
            ],
            {i, 0, halfFloor - 1}
          ],
        modulus
      ]
    ];

    (* X_(2(k+jt)), 1 <= jt <= floor(n/2). *)
    Do[
      f[2 jt + 1] =
        e u[(n + 1)/2 + jt]^2 +
          Sum[
            u[halfFloor - i + jt] * Expand[
              2 u[halfCeiling + 1 + i + jt] -
                p[[2 i + e + 1]] u[halfFloor - i + jt] -
                2 Sum[
                  p[[j + 1]] u[halfCeiling + i - j + jt],
                  {j, 0, 2 i - 1 + e}
                ]
            ],
            {i, 0, halfFloor - 1}
          ],
      {jt, 1, halfFloor}
    ];

    (* X_(2 k + 1). *)
    f[2] =
      p[[n]] u[1]^2 + (1 - e) u[(n + 2)/2]^2 +
        Sum[
          u[halfCeiling - i] * Expand[
            2 u[halfFloor + 2 + i] -
              p[[2 i + 2 - e]] u[halfCeiling - i] -
              2 Sum[
                p[[j + 1]] u[halfFloor + 1 + i - j],
                {j, 0, 2 i - e}
              ]
          ],
          {i, 0, halfCeiling - 2}
        ];

    (* X_(2 k + 2 jt + 1), 1 <= jt <= floor(n/2). *)
    Do[
      f[2 jt + 2] =
        p[[n]] u[1 + jt]^2 +
          (1 - e) u[(n + 2)/2 + jt]^2 +
          Sum[
            u[halfCeiling - i + jt] * Expand[
              2 u[halfFloor + 2 + i + jt] -
                p[[2 i + 2 - e]] u[halfCeiling - i + jt] -
                2 Sum[
                  p[[j + 1]] u[halfFloor + 1 + i - j + jt],
                  {j, 0, 2 i - e}
                ]
            ],
            {i, 0, halfCeiling - 2}
          ],
      {jt, 1, halfFloor}
    ]
  ];

  (* Read the binary expansion after its leading 1.  The two branches are
     exactly the original X_(2k+1) and X_(2k) transitions. *)
  Do[
    Do[x[i] = normalize[X[i]], {i, 1, n}];

    If[bits[[bitPosition]] === 1,
      Do[
        X[i] = normalize[f[i + 1] /. u -> x],
        {i, 1, n - 1}
      ];
      X[n] = normalize[f[n + 1] /. u -> x],

      Do[
        X[i] = normalize[f[i] /. u -> x],
        {i, 1, n}
      ]
    ],
    {bitPosition, 2, Length[bits]}
  ];

  (* The final linear form from the supplied initial values. *)
  weights = Table[
    a[[n - j]] -
      Sum[
        a[[n - j - 1 - i]] p[[i + 1]],
        {i, 0, n - j - 2}
      ],
    {j, 0, n - 1}
  ];

  termAt[r_Integer] := normalize[
    Sum[weights[[j + 1]] X[j + r], {j, 0, n - 1}]
  ];

  If[termCount === 1,
    Return[termAt[1]]
  ];

  (* Extend the X-sequence only as far as the requested output window. *)
  Do[
    X[n + 1 + l] = normalize[
      Sum[X[n + l - i] p[[i + 1]], {i, 0, n - 1}]
    ],
    {l, 0, termCount - 2}
  ];

  Table[termAt[r], {r, 1, termCount}]
];

(* ---------------------------------------------------------------------- *)
(* Toeplitz determinant/permanent wrappers                                *)
(* ---------------------------------------------------------------------- *)

TermForDTM::pos =
  "Fiduccia remote-term evaluation requires a constant-coefficient recurrence; PositionDependent -> True is not supported.";
TermForPTM::pos = TermForDTM::pos;
TermForDTM::method =
  "ScalarRecurrenceMethod must be \"CharacteristicPolynomial\" or \"Krylov\".";
TermForPTM::method = TermForDTM::method;
TermForDTM::krylov =
  "ScalarRecurrenceMethod -> \"Krylov\" requires ToeplitzRecurrencesKrylov.wl to be loaded before this extension.";
TermForPTM::krylov = TermForDTM::krylov;
TermForDTM::constructor =
  "The scalar recurrence constructor failed for (`1`,`2`).";
TermForPTM::constructor = TermForDTM::constructor;
TermForDTM::scalar =
  "No usable constant scalar recurrence polynomial was returned for (`1`,`2`).";
TermForPTM::scalar = TermForDTM::scalar;
TermForDTM::order =
  "The returned scalar recurrence has order `1`; FiducciaPolSquarings requires order at least 2.";
TermForPTM::order = TermForDTM::order;
TermForDTM::diag =
  "DiagonalValues must be Automatic or a list of length m1+m2+1.";
TermForPTM::diag = TermForDTM::diag;
TermForDTM::moddiag =
  "With a finite Modulus, the Toeplitz diagonal values must all be integers.";
TermForPTM::moddiag = TermForDTM::moddiag;

Options[TermForDTM] = {
  DiagonalValues -> Automatic,
  ScalarRecurrenceMethod -> "CharacteristicPolynomial",
  Modulus -> None,
  TermCount -> 1,
  Verbose -> False,
  PositionDependent -> False
};
Options[TermForPTM] = Options[TermForDTM];

ClearAll[trFiducciaOptionSupportedQ];
trFiducciaOptionSupportedQ[constructor_, option_] :=
  !FreeQ[Options[constructor], HoldPattern[option -> _]];

ClearAll[trFiducciaToeplitzMatrix];
trFiducciaToeplitzMatrix[m1_Integer, m2_Integer, diagonals_List, size_Integer] :=
  Table[
    If[
      -m1 <= j - i <= m2,
      diagonals[[j - i + m1 + 1]],
      0
    ],
    {i, 1, size}, {j, 1, size}
  ];

ClearAll[trFiducciaPermanent];
trFiducciaPermanent[matrix_List, modulus_] := Module[
  {n = Length[matrix], normalize, subsets, terms},

  If[n === 0, Return[1]];

  normalize = If[
    modulus === None,
    Function[value, Expand[value]],
    Function[value, Mod[value, modulus]]
  ];

  subsets = Subsets[Range[n]];
  terms = Map[
    Function[subset,
      (-1)^(n - Length[subset]) *
        Product[
          normalize[Total[matrix[[i, subset]]]],
          {i, 1, n}
        ]
    ],
    subsets
  ];

  normalize[Total[terms]]
];

ClearAll[trFiducciaInitialToeplitzTerms];
trFiducciaInitialToeplitzTerms[
    quantity_, m1_Integer, m2_Integer, diagonals_List,
    order_Integer, modulus_] := Module[
  {normalize, valueAt},

  normalize = If[
    modulus === None,
    Function[value, Expand[value]],
    Function[value, Mod[value, modulus]]
  ];

  valueAt[0] = 1;
  valueAt[size_Integer?Positive] := valueAt[size] = Module[{matrix},
    matrix = trFiducciaToeplitzMatrix[m1, m2, diagonals, size];
    normalize[
      Switch[quantity,
        "Determinant", Det[matrix],
        "Permanent", trFiducciaPermanent[matrix, modulus]
      ]
    ]
  ];

  Table[valueAt[size], {size, 0, order - 1}]
];

ClearAll[trFiducciaRecurrenceCoefficients];
trFiducciaRecurrenceCoefficients[polynomial_, variable_] := Module[
  {degree, leading, monic},

  degree = Exponent[polynomial, variable];
  If[!IntegerQ[degree] || degree < 1, Return[$Failed]];

  leading = Coefficient[polynomial, variable, degree];
  If[TrueQ[leading === 0], Return[$Failed]];

  monic = Expand[Cancel[Together[polynomial/leading]]];
  Table[
    -Coefficient[monic, variable, degree - i],
    {i, 1, degree}
  ]
];

ClearAll[trFiducciaTermForToeplitz];
trFiducciaTermForToeplitz[
    constructor_, quantity_, messageHead_,
    m1_Integer?NonNegative, m2_Integer?NonNegative, index_,
    diagonalValues_, method_, modulus_, termCount_, verbose_,
    positionDependent_] := Module[
  {
    variable = Unique["x$"], constructorRules, data,
    polynomial, coefficients, order, diagonals, initials
  },

  If[TrueQ[positionDependent],
    Message[messageHead::pos];
    Return[$Failed]
  ];

  If[!MemberQ[{"CharacteristicPolynomial", "Krylov"}, method],
    Message[messageHead::method];
    Return[$Failed]
  ];

  If[!(diagonalValues === Automatic ||
       (ListQ[diagonalValues] && Length[diagonalValues] === m1 + m2 + 1)),
    Message[messageHead::diag];
    Return[$Failed]
  ];

  If[method === "Krylov" &&
     !trFiducciaOptionSupportedQ[constructor, ScalarRecurrenceMethod],
    Message[messageHead::krylov];
    Return[$Failed]
  ];

  constructorRules = {
    DiagonalValues -> diagonalValues,
    ScalarRecurrence -> True,
    CharacteristicPolynomialVariable -> variable,
    Verbose -> verbose,
    PositionDependent -> False
  };

  If[trFiducciaOptionSupportedQ[constructor, ScalarRecurrenceMethod],
    constructorRules = Append[
      constructorRules,
      ScalarRecurrenceMethod -> method
    ]
  ];

  constructorRules = FilterRules[constructorRules, Options[constructor]];
  data = constructor[m1, m2, Sequence @@ constructorRules];

  If[data === $Failed || !AssociationQ[data],
    Message[messageHead::constructor, m1, m2];
    Return[$Failed]
  ];

  polynomial = If[
    method === "Krylov",
    Lookup[data, "ScalarRecurrencePolynomial", Missing["NotComputed"]],
    Lookup[data, "CharacteristicPolynomial", Missing["NotComputed"]]
  ];

  If[MissingQ[polynomial] || polynomial === $Failed,
    Message[messageHead::scalar, m1, m2];
    Return[$Failed]
  ];

  coefficients = trFiducciaRecurrenceCoefficients[polynomial, variable];
  If[coefficients === $Failed,
    Message[messageHead::scalar, m1, m2];
    Return[$Failed]
  ];

  order = Length[coefficients];
  If[order < 2,
    Message[messageHead::order, order];
    Return[$Failed]
  ];

  diagonals = Lookup[data, "DiagonalValues", diagonalValues];
  If[!ListQ[diagonals] || Length[diagonals] =!= m1 + m2 + 1,
    Message[messageHead::diag];
    Return[$Failed]
  ];

  If[modulus =!= None && !AllTrue[diagonals, IntegerQ],
    Message[messageHead::moddiag];
    Return[$Failed]
  ];

  initials = trFiducciaInitialToeplitzTerms[
    quantity, m1, m2, diagonals, order, modulus
  ];

  FiducciaPolSquarings[
    initials, coefficients, index,
    Modulus -> modulus,
    TermCount -> termCount
  ]
];

TermForDTM[
    m1_Integer?NonNegative, m2_Integer?NonNegative, index_,
    OptionsPattern[]] :=
  trFiducciaTermForToeplitz[
    LinRecForDTM, "Determinant", TermForDTM,
    m1, m2, index,
    OptionValue[DiagonalValues],
    OptionValue[ScalarRecurrenceMethod],
    OptionValue[Modulus],
    OptionValue[TermCount],
    OptionValue[Verbose],
    OptionValue[PositionDependent]
  ];

TermForPTM[
    m1_Integer?NonNegative, m2_Integer?NonNegative, index_,
    OptionsPattern[]] :=
  trFiducciaTermForToeplitz[
    LinRecForPTM, "Permanent", TermForPTM,
    m1, m2, index,
    OptionValue[DiagonalValues],
    OptionValue[ScalarRecurrenceMethod],
    OptionValue[Modulus],
    OptionValue[TermCount],
    OptionValue[Verbose],
    OptionValue[PositionDependent]
  ];

End[];
EndPackage[];
