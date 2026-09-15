# Verification

This release keeps the exact-ring policy of v3 and the incremental Krylov optimization of v4, and fixes SageMath compatibility of internal temporary polynomial-variable names.

## Verified in the export environment

- CPython backend-independent regression suite: 45 tests passed, including dedicated tests for the incremental Krylov row reducer.
- The fresh remote-term polynomial variable is regression-tested to start with a letter (`toeplitz_lambda_N`), avoiding SageMath `PolynomialRing` rejection of the former leading-underscore name.
- Random exact dependence certificates are cross-checked against SymPy nullspaces.
- The Krylov implementation no longer recomputes `right_kernel()` at every step; it uses incremental fraction-free row elimination and normalizes only the first dependence relation.
- Generic constant Toeplitz parameters remain in exact polynomial rings over `QQ`; the coefficient fraction field is introduced only when a Krylov dependence relation requires division.
- `tests/test_sage_krylov.sage` uses separate five-parameter and three-parameter polynomial rings, avoids a duplicated characteristic-polynomial computation already covered by `test_sage_row_column.sage`, and checks the observable determinant in its polynomial ring.
- The no-symbolic-ring regression continues to cover production source, Sage test scripts, paper examples, and README.
- Python-compatible syntax compilation passed for 18 `.py`/`.sage` source, test, example, and verifier files.
- Release verification checks layout, API/default documentation, exact-ring documentation, absence of symbolic-ring dependencies, citation metadata, ASCII public text, and manifest integrity.

## SageMath runtime check still required

The export environment does not provide a `sage` executable, so the Sage-specific runtime suite could not be executed here. On a SageMath installation, run from the repository root:

```bash
sage tests/run_tests.sage
```

This is the decisive runtime check for Sage matrix coercions, polynomial/fraction-field Krylov calculations, position-dependent exact transfer construction, direct determinant/permanent comparisons, and the paper-facing examples.
