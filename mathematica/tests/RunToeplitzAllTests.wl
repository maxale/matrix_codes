(* Load all regression suites and run them. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "ToeplitzRecurrencesExamples.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "ToeplitzRecurrencesKrylovExamples.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "ToeplitzRecurrencesFiducciaExamples.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "ToeplitzSymmetryReductionsExamples.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "ToeplitzCirculantBridgeExamples.wl"}]];

ClearAll[RunToeplitzAllTests];
RunToeplitzAllTests[] := <|
  "Base" -> RunToeplitzRecurrenceTests[],
  "Krylov" -> RunToeplitzKrylovTests[],
  "Fiduccia" -> RunToeplitzFiducciaTests[],
  "Symmetry" -> RunToeplitzSymmetryReductionTests[],
  "CirculantBridge" -> RunToeplitzCirculantBridgeTests[]
|>;
