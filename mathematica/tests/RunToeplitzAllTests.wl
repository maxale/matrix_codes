(* Load all regression suites and run them. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "ToeplitzRecurrencesExamples.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "ToeplitzRecurrencesKrylovExamples.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "ToeplitzRecurrencesFiducciaExamples.wl"}]];

ClearAll[RunToeplitzAllTests];
RunToeplitzAllTests[] := <|
  "Base" -> RunToeplitzRecurrenceTests[],
  "Krylov" -> RunToeplitzKrylovTests[],
  "Fiduccia" -> RunToeplitzFiducciaTests[]
|>;
