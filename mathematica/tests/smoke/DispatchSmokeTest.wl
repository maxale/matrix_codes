(* Focused regression test for ScalarRecurrenceMethod dispatch. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "src", "ToeplitzRecurrencesAll.wl"}]];

TestReport[{
  VerificationTest[
    AssociationQ[LinRecForDTM[1, 1, Verbose -> False]],
    True,
    TestID -> "DTM default method returns Association"
  ],
  VerificationTest[
    AssociationQ[LinRecForDTM[1, 1,
      ScalarRecurrenceMethod -> "CharacteristicPolynomial",
      Verbose -> False]],
    True,
    TestID -> "DTM explicit characteristic method returns Association"
  ],
  VerificationTest[
    AssociationQ[LinRecForDTM[1, 1,
      ScalarRecurrenceMethod -> "Krylov",
      ScalarRecurrence -> True,
      Verbose -> False]],
    True,
    TestID -> "DTM explicit Krylov method returns Association"
  ],
  VerificationTest[
    AssociationQ[LinRecForPTM[1, 1, Verbose -> False]],
    True,
    TestID -> "PTM default method returns Association"
  ],
  VerificationTest[
    AssociationQ[LinRecForPTM[1, 1,
      ScalarRecurrenceMethod -> "Krylov",
      ScalarRecurrence -> True,
      Verbose -> False]],
    True,
    TestID -> "PTM explicit Krylov method returns Association"
  ],
  VerificationTest[
    Quiet[LinRecForDTM[1, 1,
      ScalarRecurrenceMethod -> "NotAMethod",
      Verbose -> False], LinRecForDTM::scalarmethod],
    $Failed,
    TestID -> "invalid method is rejected"
  ]
}]
