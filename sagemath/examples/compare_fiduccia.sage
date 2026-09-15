"""Compare the preserved Fiduccia implementation with Sage CFiniteSequences.

Run with

    sage compare_fiduccia.sage

The documented Sage C-finite API is used as an independent reference path.
The script intentionally reports construction+evaluation time for both public
wrappers, so the comparison is reproducible and easy to modify.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


from time import perf_counter

from toeplitz_recurrences import fiduccia_pol_squarings, sage_cfinite_term


def timed_call(function, *args, repetitions=3, **kwargs):
    best = None
    value = None
    for _ in range(repetitions):
        start = perf_counter()
        candidate = function(*args, **kwargs)
        elapsed = perf_counter() - start
        if best is None or elapsed < best:
            best = elapsed
            value = candidate
    return best, value


cases = [
    ("order 3", [1, 23, 3], [4, 5, 6], 1000),
    ("order 4", [2, -1, 5, 4], [3, -2, 7, 1], 2000),
]

print("Fiduccia hardcoded squarings vs Sage CFiniteSequences")
print("Times are best of three and include public-wrapper setup.\n")

for label, initial, recurrence, n in cases:
    tf, vf = timed_call(fiduccia_pol_squarings, initial, recurrence, n)
    ts, vs = timed_call(sage_cfinite_term, initial, recurrence, n)
    if vf != vs:
        raise AssertionError(f"{label}: implementations disagree at n={n}")
    print(f"{label:8s}  n={n:6d}  hardcoded={tf:.6f}s  CFinite={ts:.6f}s")

print("\nCorrectness checks passed for all benchmark cases.")
print("For modular huge-index timings, benchmark fiduccia_pol_squarings(..., modulus=m) separately;")
print("the Sage CFiniteSequences reference used here is over QQ rather than residue rings.")
