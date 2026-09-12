(* ::Package:: *)

(* Regression test for namespace isolation of default symbolic parameters.
   Run in a fresh Mathematica kernel from the release directory. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "src", "ToeplitzRecurrencesAll.wl"}]];

ClearAll[RunToeplitzScopeSmokeTest];

RunToeplitzScopeSmokeTest[] := Block[
  {Global`a = 9, Global`x = 17, Global`k = 23},
  Module[{det, per, pos, inc, kry, expectedDet11, tests},
    det = LinRecForDTM[1, 1,
      ScalarRecurrence -> True,
      Verbose -> False];
    per = LinRecForPTM[1, 1,
      ScalarRecurrence -> False,
      Verbose -> False];
    pos = LinRecForDTM[1, 1,
      PositionDependent -> True,
      ScalarRecurrence -> False,
      Verbose -> False];
    inc = LinRecForDTMIncreasingRows[1, 1,
      Verbose -> False];
    kry = LinRecForDTM[1, 1,
      ScalarRecurrence -> True,
      ScalarRecurrenceMethod -> "Krylov",
      Verbose -> False];

    expectedDet11 = {
      {\[FormalA][0], -\[FormalA][-1]},
      {\[FormalA][1], 0}
    };

    tests = {
      VerificationTest[
        Normal[det["TransferMatrix"]],
        expectedDet11,
        TestID -> "global a value does not affect default determinant diagonals"
      ],
      VerificationTest[
        Normal[per["TransferMatrix"]],
        {{\[FormalA][0], \[FormalA][-1]}, {\[FormalA][1], 0}},
        TestID -> "global a value does not affect default permanent diagonals"
      ],
      VerificationTest[
        det["CharacteristicPolynomialVariable"],
        \[FormalX],
        TestID -> "global x value does not affect characteristic variable"
      ],
      VerificationTest[
        det["States"][[1]],
        {{\[FormalK] + 1}, {\[FormalK] + 1}},
        TestID -> "global k value does not affect boundary-state parameter"
      ],
      VerificationTest[
        Normal[pos["TransferMatrixExpression"]][[1, 1]],
        \[FormalA][0][\[FormalK]],
        TestID -> "global a and k values do not affect position-dependent defaults"
      ],
      VerificationTest[
        inc["CharacteristicPolynomialVariable"],
        \[FormalX],
        TestID -> "increasing-rows defaults use formal characteristic variable"
      ],
      VerificationTest[
        kry["ScalarRecurrencePolynomial"],
        \[FormalX]^2 - \[FormalA][0] \[FormalX] + \[FormalA][-1] \[FormalA][1],
        TestID -> "Krylov defaults are isolated from Global symbols"
      ],
      VerificationTest[
        {Global`a, Global`x, Global`k},
        {9, 17, 23},
        TestID -> "recurrence constructors do not modify user Global values"
      ]
    };

    TestReport[tests]
  ]
];

RunToeplitzScopeSmokeTest[]
