#!/usr/bin/env python3
"""Static integrity checks for the Toeplitz recurrence Wolfram release."""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

PAPER_TITLE = "Constructive recurrences for determinants and permanents of banded Toeplitz matrices"

REQUIRED = [
    "README.md",
    "CITATION.cff",
    "CHANGELOG.md",
    "PACKAGE_INFO.txt",
    "VERIFICATION.md",
    ".gitignore",
    "MANIFEST.sha256",
    "src/ToeplitzRecurrences.wl",
    "src/ToeplitzRecurrencesKrylov.wl",
    "src/ToeplitzRecurrencesFiduccia.wl",
    "src/ToeplitzRecurrencesAll.wl",
    "examples/PaperExamples.wl",
    "tests/RunToeplitzAllTests.wl",
    "tests/ToeplitzRecurrencesExamples.wl",
    "tests/ToeplitzRecurrencesKrylovExamples.wl",
    "tests/ToeplitzRecurrencesFiducciaExamples.wl",
    "tests/smoke/DispatchSmokeTest.wl",
    "tests/smoke/KrylovDispatchSmokeTest.wl",
    "tests/smoke/FiducciaSmokeTest.wl",
    "tests/smoke/ScopeSmokeTest.wl",
    "tools/verify_release.py",
]

FORBIDDEN_MANUSCRIPT_MARKERS = [
    re.compile(r"toeplitz_recurrences_rev\d+", re.I),
    re.compile(r"paper\s+rev(?:ision)?\s*\d+", re.I),
    re.compile(r"paper\s+revision\s*:", re.I),
    re.compile(r"rev29|rev36", re.I),
]

PUBLIC_APIS = [
    "LinRecForDTM",
    "LinRecForPTM",
    "LinRecForDTMIncreasingRows",
    "LinRecForPTMIncreasingRows",
    "ScalarRecurrenceMethod",
    "PositionDependent",
    "FiducciaPolSquarings",
    "TermForDTM",
    "TermForPTM",
]

ASCII_PUBLIC_SUFFIXES = {".md", ".txt", ".cff", ".py"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check_wolfram_balance(path: Path) -> list[str]:
    """Check (), [], {}, strings, and nested Wolfram comments without evaluating code."""
    text = path.read_text(encoding="utf-8")
    stack: list[tuple[str, int]] = []
    pairs = {")": "(", "]": "[", "}": "{"}
    openers = set(pairs.values())
    errors: list[str] = []
    i = 0
    in_string = False
    escaped = False
    comment_depth = 0
    while i < len(text):
        ch = text[i]
        nxt = text[i + 1] if i + 1 < len(text) else ""
        if comment_depth:
            if ch == "(" and nxt == "*":
                comment_depth += 1
                i += 2
                continue
            if ch == "*" and nxt == ")":
                comment_depth -= 1
                i += 2
                continue
            i += 1
            continue
        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == '"':
                in_string = False
            i += 1
            continue
        if ch == "(" and nxt == "*":
            comment_depth = 1
            i += 2
            continue
        if ch == '"':
            in_string = True
            i += 1
            continue
        if ch in openers:
            stack.append((ch, i))
        elif ch in pairs:
            if not stack or stack[-1][0] != pairs[ch]:
                errors.append(f"{path}: unmatched {ch!r} near byte {i}")
                return errors
            stack.pop()
        i += 1
    if comment_depth:
        errors.append(f"{path}: unterminated comment")
    if in_string:
        errors.append(f"{path}: unterminated string")
    if stack:
        errors.append(f"{path}: unclosed delimiter {stack[-1][0]!r} near byte {stack[-1][1]}")
    return errors


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    errors: list[str] = []

    for rel in REQUIRED:
        if not (root / rel).is_file():
            errors.append(f"missing required file: {rel}")

    text_files = [
        p for p in root.rglob("*")
        if p.is_file() and p.name not in {"MANIFEST.sha256", "verify_release.py"}
        and p.suffix.lower() in {".wl", ".md", ".txt", ".cff", ".py"}
    ]
    for path in text_files:
        text = path.read_text(encoding="utf-8")
        for pattern in FORBIDDEN_MANUSCRIPT_MARKERS:
            if pattern.search(text):
                errors.append(f"manuscript filename/revision marker in {path.relative_to(root)}")
                break

    for path in root.rglob("*"):
        if path.is_file() and path.suffix.lower() in ASCII_PUBLIC_SUFFIXES:
            data = path.read_bytes()
            if any(byte >= 128 for byte in data):
                errors.append(f"non-ASCII byte in public text file: {path.relative_to(root)}")

    for rel in ["README.md", "PACKAGE_INFO.txt", "CITATION.cff", "src/ToeplitzRecurrences.wl"]:
        path = root / rel
        if path.is_file() and PAPER_TITLE not in path.read_text(encoding="utf-8"):
            errors.append(f"paper title missing from {rel}")

    cff = (root / "CITATION.cff").read_text(encoding="utf-8") if (root / "CITATION.cff").is_file() else ""
    if 'version: "10"' not in cff:
        errors.append("CITATION.cff code version is not 10")
    if 'title: "Toeplitz recurrence algorithms - Wolfram Language"' not in cff:
        errors.append("CITATION.cff software title is stale")

    source = "\n".join(
        p.read_text(encoding="utf-8") for p in sorted((root / "src").glob("*.wl"))
    ) if (root / "src").is_dir() else ""
    for api in PUBLIC_APIS:
        if api not in source:
            errors.append(f"public interface not found in src/: {api}")
    for bad in ["Global`a", "Global`x", "Global`k"]:
        if bad in source:
            errors.append(f"production source contains legacy global placeholder {bad}")

    source_option_tokens = {
        "src/ToeplitzRecurrences.wl": [
            "DiagonalValues -> Automatic",
            "ScalarRecurrence -> Automatic",
            "CharacteristicPolynomialVariable -> \\[FormalX]",
            "MaxCharacteristicOrder -> 20",
            "Verbose -> True",
            "PositionDependent -> False",
            "BandEntryFunction -> Automatic",
            "EllRange -> \"Forward\"",
        ],
        "src/ToeplitzRecurrencesKrylov.wl": [
            "ScalarRecurrenceMethod -> \"CharacteristicPolynomial\"",
        ],
        "src/ToeplitzRecurrencesFiduccia.wl": [
            "Modulus -> None",
            "TermCount -> 1",
            "ScalarRecurrenceMethod -> \"CharacteristicPolynomial\"",
            "Verbose -> False",
            "PositionDependent -> False",
        ],
    }
    for rel, tokens in source_option_tokens.items():
        path = root / rel
        text = path.read_text(encoding="utf-8") if path.is_file() else ""
        for token in tokens:
            if token not in text:
                errors.append(f"source option/default missing from {rel}: {token}")

    readme = (root / "README.md").read_text(encoding="utf-8") if (root / "README.md").is_file() else ""
    for token in [
        "increasing-rows", "row-column", "Krylov", "PositionDependent",
        "FiducciaPolSquarings", "TermForDTM", "TermForPTM",
        "threefold root-of-unity geometry", "radial interlacing",
        "asymptotically sharp outer radius", "offset `+/-2` parity-split comparison",
    ]:
        if token.lower() not in readme.lower():
            errors.append(f"README missing paper-to-code topic: {token}")
    for scope_phrase in ["cyclic-closure", "finite-corner-defect", "nonautonomous"]:
        if scope_phrase.lower() not in readme.lower():
            errors.append(f"README missing explicit scope boundary: {scope_phrase}")

    option_doc_tokens = [
        "## Public API and options",
        "`LinRecForDTM[m1,m2,...]`",
        "`LinRecForPTM[m1,m2,...]`",
        "`LinRecForDTMIncreasingRows[m1,m2,...]`",
        "`LinRecForPTMIncreasingRows[m1,m2,...]`",
        "`LaplaceD[m1,m2,S,...]`",
        "`LaplaceP[m1,m2,S,...]`",
        "`FiducciaPolSquarings[initial,recurrence,N,...]`",
        "`TermForDTM[m1,m2,N,...]`",
        "`TermForPTM[m1,m2,N,...]`",
        "`DiagonalValues` | `Automatic`",
        "`ScalarRecurrence` | `Automatic`",
        "`CharacteristicPolynomialVariable` | `\\[FormalX]`",
        "`MaxCharacteristicOrder` | `20`",
        "`EllRange` | `\"Forward\"`",
        "`ScalarRecurrenceMethod` | `\"CharacteristicPolynomial\"`",
        "`Modulus` | `None`",
        "`TermCount` | `1`",
        "`Verbose` | `False`",
        "`PositionDependent` | `False`",
    ]
    for token in option_doc_tokens:
        if token not in readme:
            errors.append(f"README missing public option/default documentation: {token}")

    for token in [
        "Efficient Computation of Terms of Linear Recurrence Sequences of Any Order",
        "INTEGERS 18",
        "10.1137/0214007",
        "modular-polynomial-squaring formulas",
    ]:
        if token not in readme:
            errors.append(f"README missing distant-term reference/provenance: {token}")

    examples = (root / "examples/PaperExamples.wl").read_text(encoding="utf-8") if (root / "examples/PaperExamples.wl").is_file() else ""
    for token in [
        "LinRecForDTM[2, 2",
        "LinRecForPTM[2, 2",
        "LinRecForDTM[3, 3",
        "LinRecForDTM[2, 1",
        "ScalarRecurrenceMethod -> \"Krylov\"",
        "PositionDependent -> True",
        "FiducciaPolSquarings",
        "ExplicitFormulaCheck",
        "Binomial[n - 2 j, j]",
        "ParitySplitCheck",
        "Ceiling[n / 2]",
        "Floor[n / 2]",
    ]:
        if token not in examples:
            errors.append(f"PaperExamples missing: {token}")

    for path in root.rglob("*.wl"):
        errors.extend(check_wolfram_balance(path))

    manifest_path = root / "MANIFEST.sha256"
    if manifest_path.is_file():
        listed: set[str] = set()
        for lineno, line in enumerate(manifest_path.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip():
                continue
            try:
                digest, rel = line.split("  ", 1)
            except ValueError:
                errors.append(f"MANIFEST.sha256:{lineno}: malformed line")
                continue
            listed.add(rel)
            path = root / rel
            if not path.is_file():
                errors.append(f"manifest entry missing: {rel}")
            elif sha256(path) != digest:
                errors.append(f"manifest hash mismatch: {rel}")
        expected = {
            str(p.relative_to(root)) for p in root.rglob("*")
            if p.is_file() and p.name != "MANIFEST.sha256"
        }
        if listed != expected:
            for rel in sorted(expected - listed):
                errors.append(f"file absent from manifest: {rel}")
            for rel in sorted(listed - expected):
                errors.append(f"manifest lists unexpected file: {rel}")

    if errors:
        print(f"FAIL ({len(errors)} issue(s))")
        for error in errors:
            print(f"- {error}")
        return 1

    print("PASS")
    print(f"Checked {len(REQUIRED)} required paths, {len(list(root.rglob('*.wl')))} Wolfram files, and manifest integrity.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
