#!/usr/bin/env python3
"""Static release verifier for the SageMath Toeplitz recurrence package."""

from __future__ import annotations

import ast
import hashlib
import re
import sys
from pathlib import Path

PAPER_TITLE = "Constructive recurrences for determinants and permanents of banded Toeplitz matrices"
REQUIRED = [
    "README.md", "CITATION.cff", "CHANGELOG.md", "PACKAGE_INFO.txt", "VERIFICATION.md", ".gitignore",
    "src/toeplitz_recurrences.py", "src/toeplitz_recurrences.sage",
    "examples/paper_examples.sage", "examples/compare_fiduccia.sage",
    "tests/run_tests.sage", "tests/test_python_core.py", "tests/test_python_increasing_logic.py",
    "tests/test_python_package_contract.py", "tests/test_python_krylov_incremental.py", "tests/test_python_paper_formulas.py",
    "tests/test_python_release_contract.py", "tests/test_python_state_logic.py", "tests/test_sage_row_column.sage",
    "tests/test_sage_increasing_rows.sage", "tests/test_sage_krylov.sage",
    "tests/test_sage_fiduccia_integration.sage", "tests/test_python_no_symbolic_ring.py", "tools/verify_release.py",
]
TEXT_SUFFIXES = {".py", ".sage", ".md", ".txt", ".cff"}
FORBIDDEN = [
    re.compile(r"toeplitz_recurrences_rev\d+\.tex", re.I),
    re.compile(r"paper\s+rev\d+", re.I),
    re.compile(r"manuscript\s+rev\d+", re.I),
    re.compile("rev" + r"(?:36|009)", re.I),
]
PUBLIC_APIS = [
    "m_reduce", "laplace_d", "laplace_p", "lin_rec_for_dtm", "lin_rec_for_ptm",
    "lin_rec_for_dtm_increasing_rows", "lin_rec_for_ptm_increasing_rows",
    "fiduccia_pol_squarings", "term_for_dtm", "term_for_ptm", "sage_cfinite_term",
]
EXPECTED_DEFAULTS = {
    "fiduccia_pol_squarings": {"modulus": "None", "term_count": "1"},
    "laplace_d": {"diagonal_values": "None", "position_dependent": "False", "band_entry_function": "None"},
    "laplace_p": {"diagonal_values": "None", "position_dependent": "False", "band_entry_function": "None"},
    "lin_rec_for_dtm": {
        "diagonal_values": "None", "scalar_recurrence": "'auto'", "characteristic_polynomial_variable": "'x'",
        "max_characteristic_order": "20", "verbose": "True", "position_dependent": "False",
        "band_entry_function": "None", "scalar_recurrence_method": "'characteristic_polynomial'",
    },
    "lin_rec_for_ptm": {
        "diagonal_values": "None", "scalar_recurrence": "'auto'", "characteristic_polynomial_variable": "'x'",
        "max_characteristic_order": "20", "verbose": "True", "position_dependent": "False",
        "band_entry_function": "None", "scalar_recurrence_method": "'characteristic_polynomial'",
    },
    "lin_rec_for_dtm_increasing_rows": {
        "diagonal_values": "None", "ell_range": "'forward'", "characteristic_polynomial_variable": "'x'", "verbose": "True",
    },
    "lin_rec_for_ptm_increasing_rows": {
        "diagonal_values": "None", "ell_range": "'forward'", "characteristic_polynomial_variable": "'x'", "verbose": "True",
    },
    "term_for_dtm": {
        "diagonal_values": "None", "scalar_recurrence_method": "'characteristic_polynomial'",
        "modulus": "None", "term_count": "1", "verbose": "False", "position_dependent": "False",
    },
    "term_for_ptm": {
        "diagonal_values": "None", "scalar_recurrence_method": "'characteristic_polynomial'",
        "modulus": "None", "term_count": "1", "verbose": "False", "position_dependent": "False",
    },
    "sage_cfinite_term": {"term_count": "1"},
}


def is_generated_artifact(path: Path) -> bool:
    """Return True for local test/build artifacts excluded from the release manifest."""
    return (
        "__pycache__" in path.parts
        or ".pytest_cache" in path.parts
        or path.suffix.lower() in {".pyc", ".pyo"}
        or path.name.endswith(".sage.py")
    )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_defaults(source_path: Path):
    mod = ast.parse(source_path.read_text(encoding="utf-8"))
    functions = {}
    names = set()
    for node in mod.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            names.add(node.name)
            functions[node.name] = {
                arg.arg: ast.unparse(default)
                for arg, default in zip(node.args.kwonlyargs, node.args.kw_defaults)
                if default is not None
            }
    return names, functions


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    errors = []

    for rel in REQUIRED:
        if not (root / rel).is_file():
            errors.append(f"missing required file: {rel}")

    for path in root.rglob("*"):
        if (
            not path.is_file()
            or path.name == "MANIFEST.sha256"
            or is_generated_artifact(path)
        ):
            continue
        if path.suffix.lower() in TEXT_SUFFIXES:
            data = path.read_bytes()
            if any(b >= 128 for b in data):
                errors.append(f"non-ASCII byte in public text file: {path.relative_to(root)}")
            text = data.decode("ascii", errors="replace")
            for pattern in FORBIDDEN:
                if pattern.search(text):
                    errors.append(f"manuscript filename/revision marker in {path.relative_to(root)}")
                    break

    source_path = root / "src/toeplitz_recurrences.py"
    if source_path.is_file():
        source = source_path.read_text(encoding="utf-8")
        if PAPER_TITLE not in source:
            errors.append("paper title missing from Sage source header")
        try:
            names, defaults = parse_defaults(source_path)
        except SyntaxError as exc:
            errors.append(f"Python source syntax error: {exc}")
            names, defaults = set(), {}
        for name in PUBLIC_APIS:
            if name not in names:
                errors.append(f"public interface missing from source: {name}")
        for func, expected in EXPECTED_DEFAULTS.items():
            actual = defaults.get(func, {})
            for key, value in expected.items():
                if actual.get(key) != value:
                    errors.append(f"source default mismatch: {func}.{key} expected {value}, got {actual.get(key)}")

    exact_surfaces = [source_path, root / "src/toeplitz_recurrences.sage", root / "README.md"]
    exact_surfaces.extend((root / "tests").glob("*.sage"))
    exact_surfaces.extend((root / "examples").glob("*.sage"))
    symbolic_markers = ("sage.SR", "from sage.all import SR", "SR(", "var(")
    for path in exact_surfaces:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for marker in symbolic_markers:
            if marker in text:
                errors.append(f"Sage Symbolic Ring dependency in {path.relative_to(root)}: {marker}")

    readme_path = root / "README.md"
    readme = readme_path.read_text(encoding="ascii") if readme_path.is_file() else ""
    required_readme = [
        PAPER_TITLE, "## Public API and defaults", "## Paper-to-code map",
        "lin_rec_for_dtm", "lin_rec_for_ptm", "lin_rec_for_dtm_increasing_rows",
        "lin_rec_for_ptm_increasing_rows", "laplace_d", "laplace_p",
        "fiduccia_pol_squarings", "term_for_dtm", "term_for_ptm", "sage_cfinite_term",
        "threefold spectrum", "offset `+/-2` parity-split comparison",
        "cyclic-closure", "finite-corner-defect", "nonautonomous increasing-rows",
        "Efficient Computation of Terms of Linear Recurrence Sequences of Any Order",
        "INTEGERS 18", "10.1137/0214007", "hardcoded modular-polynomial-squaring formulas",
        "## Exact coefficient-ring policy", "multivariate polynomial rings over `QQ`",
        "fraction field", "integer position",
    ]
    for token in required_readme:
        if token.lower() not in readme.lower():
            errors.append(f"README missing required topic: {token}")

    cff_path = root / "CITATION.cff"
    cff = cff_path.read_text(encoding="ascii") if cff_path.is_file() else ""
    for token in [
        'title: "Toeplitz recurrence algorithms - SageMath"', 'version: "5"', PAPER_TITLE,
        'repository-code: "https://github.com/maxale/matrix_codes"',
    ]:
        if token not in cff:
            errors.append(f"CITATION.cff missing/stale: {token}")

    examples_path = root / "examples/paper_examples.sage"
    examples = examples_path.read_text(encoding="ascii") if examples_path.is_file() else ""
    for token in [
        "lin_rec_for_dtm(2, 2", "lin_rec_for_ptm(2, 2", "lin_rec_for_dtm(3, 3",
        "lin_rec_for_dtm(\n        2, 1", "scalar_recurrence_method=\"krylov\"",
        "position_dependent=True", "fiduccia_pol_squarings", "explicit_formula_check",
        "parity_split_check", "binomial(n - 2 * j, j)",
    ]:
        if token not in examples:
            errors.append(f"paper_examples.sage missing: {token}")

    manifest = root / "MANIFEST.sha256"
    if manifest.is_file():
        listed = set()
        for lineno, line in enumerate(manifest.read_text(encoding="ascii").splitlines(), 1):
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
            str(path.relative_to(root)) for path in root.rglob("*")
            if (
                path.is_file()
                and path.name != "MANIFEST.sha256"
                and not is_generated_artifact(path)
            )
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
    print(f"Checked {len(REQUIRED)} required paths, source signatures/defaults, public docs, examples, and manifest integrity.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
