(* Minimal smoke test for the Fiduccia integration. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "src", "ToeplitzRecurrencesAll.wl"}]];

optionsOK =
  Options[FiducciaPolSquarings] === {Modulus -> None, TermCount -> 1} &&
  MemberQ[Options[TermForDTM], ScalarRecurrenceMethod -> "CharacteristicPolynomial"] &&
  MemberQ[Options[TermForDTM], Modulus -> None] &&
  MemberQ[Options[TermForDTM], TermCount -> 1];

r1 = FiducciaPolSquarings[
  {1, 23, 3}, {111, 3332, 12}, 100,
  Modulus -> 123, TermCount -> 3
];
r2 = TermForDTM[1, 1, 20,
  DiagonalValues -> {2, 1, 3}, Modulus -> 101
];

If[optionsOK && r1 === {36, 11, 30} && r2 === 95,
  Print["FIDUCCIA SMOKE TEST: PASS"],
  Print["FIDUCCIA SMOKE TEST: FAIL"];
  Print[<|"Options" -> optionsOK, "Recurrence" -> r1, "ToeplitzDeterminant" -> r2|>]
];
