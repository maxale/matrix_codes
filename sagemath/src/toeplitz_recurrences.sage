# Convenience loader for Sage's load(...) workflow.
# The implementation itself is ordinary Sage-compatible Python.
import sys as _sys
from pathlib import Path as _Path
_SRC = _Path(__file__).resolve().parent
if str(_SRC) not in _sys.path:
    _sys.path.insert(0, str(_SRC))
from toeplitz_recurrences import *
del _SRC, _Path, _sys
