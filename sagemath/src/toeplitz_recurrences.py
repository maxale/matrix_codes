"""SageMath implementation accompanying
"Constructive recurrences for determinants and permanents of banded Toeplitz matrices"
by Max A. Alekseyev and Dmitry I. Khomovsky.

This module is Sage-first but keeps the boundary-state combinatorics and the
hardcoded-squaring distant-term evaluator backend-independent. Sage imports are
performed lazily so those pieces can also be regression-tested with ordinary
CPython.

Paper convention: ``m1`` counts subdiagonals (negative offsets) and ``m2``
counts superdiagonals (positive offsets). Row-column states use the paper's
first-discovery order.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, count
from math import comb, isqrt
from operator import index as integer_index
from typing import Iterable, Sequence


SOFTWARE_VERSION = "5"
PAPER_TITLE = "Constructive recurrences for determinants and permanents of banded Toeplitz matrices"
_REMOTE_VARIABLE_COUNTER = count()


def _fresh_polynomial_variable_name():
    """Return a collision-resistant temporary variable name for scalarization."""
    return f"toeplitz_lambda_{next(_REMOTE_VARIABLE_COUNTER)}"


__all__ = [
    "m_reduce",
    "laplace_d",
    "laplace_p",
    "lin_rec_for_dtm",
    "lin_rec_for_ptm",
    "lin_rec_for_dtm_increasing_rows",
    "lin_rec_for_ptm_increasing_rows",
    "fiduccia_pol_squarings",
    "term_for_dtm",
    "term_for_ptm",
    "sage_cfinite_term",
]


def _as_int(value, name: str) -> int:
    """Return an integer-like value as a Python int or raise TypeError."""
    try:
        return int(integer_index(value))
    except Exception as exc:  # pragma: no cover - exact exception is type-dependent
        raise TypeError(f"{name} must be an integer") from exc


def _exact_integer_value(value):
    """Return an exact integer value as Python int, rejecting nonintegers."""
    try:
        return int(integer_index(value))
    except Exception:
        pass
    denominator = getattr(value, "denominator", None)
    numerator = getattr(value, "numerator", None)
    try:
        den = denominator() if callable(denominator) else denominator
        num = numerator() if callable(numerator) else numerator
        if den == 1 and num is not None:
            return int(num)
    except Exception:
        pass
    raise TypeError(f"{value!r} is not an exact integer")


def _expand_if_available(value):
    """Mimic Wolfram Expand without imposing a symbolic backend."""
    expand = getattr(value, "expand", None)
    if callable(expand):
        try:
            return expand()
        except Exception:
            pass
    return value


def _normalize(value, modulus):
    if modulus is None:
        return _expand_if_available(value)
    return value % modulus


def m_reduce(rows: Iterable, cols: Iterable):
    """Normalize a row-column boundary-minor signature.

    The Sage port stores a state by integer offsets relative to the symbolic
    index ``k``.  The operation is exactly the Wolfram ``MReduce`` operation:
    sort both sides in descending order and repeatedly strip a common leading
    consecutive pair.
    """
    left = sorted(tuple(rows), reverse=True)
    right = sorted(tuple(cols), reverse=True)
    while (
        len(left) >= 2
        and len(right) >= 2
        and left[1] == right[1]
        and left[0] == right[0]
        and left[0] == left[1] + 1
    ):
        left = left[1:]
        right = right[1:]
    return tuple(left), tuple(right)


def _fiduccia_step(x, recurrence, bit, modulus):
    """One preserved hardcoded polynomial-squaring transition.

    ``x`` represents the current residue coefficients X[1],...,X[n].  This is
    the direct arithmetic form of the co-author's symbolic f[i]/u[i]
    construction: extend u by the recurrence, evaluate the hardcoded even and
    odd squaring identities, and choose X_(2k) or X_(2k+1) from the next binary
    digit.
    """
    p = recurrence
    n = len(p)

    # 1-based storage mirrors the source formulas and substantially reduces
    # transcription risk.
    u = [None] + [_normalize(v, modulus) for v in x]

    if n == 2:
        u1, u2 = u[1], u[2]
        f1 = u1 * (-p[0] * u1 + 2 * u2)
        f2 = p[1] * u1 * u1 + u2 * u2

        # Preserve the supplied difference-of-squares micro-optimization when
        # it applies over ordinary integers.  It is algebraically identical.
        try:
            q = -integer_index(p[1])
            if q >= 0:
                root = isqrt(q)
                if root * root == q:
                    f2 = (u2 - root * u1) * (u2 + root * u1)
        except Exception:
            pass

        f3 = u2 * (p[0] * u2 + 2 * p[1] * u1)
        f = [None, _normalize(f1, modulus), _normalize(f2, modulus), _normalize(f3, modulus)]
        if bit:
            return [f[2], f[3]]
        return [f[1], f[2]]

    e = n % 2
    half_floor = n // 2
    half_ceil = (n + 1) // 2

    # Extend u[n+1],...,u[n+floor(n/2)] as in the original code.
    for shift in range(1, half_floor + 1):
        idx = n + shift
        value = sum(p[j] * u[idx - 1 - j] for j in range(n))
        u.append(_normalize(value, modulus))

    def even_formula(jt: int):
        value = 0
        if e:
            value += u[(n + 1) // 2 + jt] ** 2
        for i in range(half_floor):
            left = u[half_floor - i + jt]
            inner = 2 * u[half_ceil + 1 + i + jt]
            inner -= p[2 * i + e] * left
            upper = 2 * i - 1 + e
            if upper >= 0:
                inner -= 2 * sum(
                    p[j] * u[half_ceil + i - j + jt]
                    for j in range(upper + 1)
                )
            value += left * inner
        return _normalize(value, modulus)

    def odd_formula(jt: int):
        value = p[n - 1] * u[1 + jt] ** 2
        if not e:
            value += u[(n + 2) // 2 + jt] ** 2
        for i in range(half_ceil - 1):
            left = u[half_ceil - i + jt]
            inner = 2 * u[half_floor + 2 + i + jt]
            inner -= p[2 * i + 1 - e] * left
            upper = 2 * i - e
            if upper >= 0:
                inner -= 2 * sum(
                    p[j] * u[half_floor + 1 + i - j + jt]
                    for j in range(upper + 1)
                )
            value += left * inner
        return _normalize(value, modulus)

    # f is 1-based.  The source computes one harmless extra even/odd entry in
    # some parities; here we compute exactly what the transition can consume.
    f = [None] * (n + 2)
    for jt in range(half_floor + 1):
        idx = 2 * jt + 1
        if idx <= n + 1:
            f[idx] = even_formula(jt)
        idx = 2 * jt + 2
        if idx <= n + 1:
            f[idx] = odd_formula(jt)

    if bit:
        return [_normalize(f[i + 1], modulus) for i in range(1, n + 1)]
    return [_normalize(f[i], modulus) for i in range(1, n + 1)]


def fiduccia_pol_squarings(
    initial: Sequence,
    recurrence: Sequence,
    n,
    *,
    modulus=None,
    term_count=1,
):
    """Jump to a remote term of a constant-coefficient recurrence.

    The recurrence convention is

    ``a[k] = recurrence[0]*a[k-1] + ... + recurrence[d-1]*a[k-d]``.

    The implementation preserves the hardcoded modular-polynomial squaring
    identities and binary transition logic used by D. I. Khomovsky in
    "Efficient Computation of Terms of Linear Recurrence Sequences of Any Order"
    (INTEGERS 18 (2018), A39), within the fast recurrence-evaluation framework
    associated with C. M. Fiduccia. ``n`` is zero-based. With
    ``term_count > 1`` a list of consecutive terms starting at ``n`` is returned.
    """
    initial = list(initial)
    recurrence = list(recurrence)
    order = len(initial)
    if order != len(recurrence) or order < 2:
        raise ValueError(
            "initial and recurrence must have the same length, at least 2"
        )

    n = _as_int(n, "n")
    if n < 0:
        raise ValueError("n must be nonnegative")

    term_count = _as_int(term_count, "term_count")
    if term_count < 1:
        raise ValueError("term_count must be positive")

    if modulus is not None:
        modulus = _as_int(modulus, "modulus")
        if modulus <= 1:
            raise ValueError("modulus must be None or an integer greater than 1")
        try:
            initial = [_exact_integer_value(value) % modulus for value in initial]
            recurrence = [_exact_integer_value(value) % modulus for value in recurrence]
        except TypeError as exc:
            raise TypeError(
                "with a finite modulus, initial values and recurrence coefficients must be integers"
            ) from exc

    # The source enters the binary jump only after the supplied initial block.
    # For a requested window that starts inside it, sequential extension is the
    # simpler equivalent path.
    if n < order:
        seq = list(initial)
        target_len = n + term_count
        while len(seq) < target_len:
            value = sum(recurrence[i] * seq[-1 - i] for i in range(order))
            seq.append(_normalize(value, modulus))
        out = seq[n:target_len]
        return out[0] if term_count == 1 else out

    bits = bin(n)[2:]
    x = [0] * (order - 2) + [1, recurrence[0]]
    for digit in bits[1:]:
        x = _fiduccia_step(x, recurrence, digit == "1", modulus)

    weights = []
    for j in range(order):
        value = initial[order - 1 - j]
        value -= sum(
            initial[order - 2 - j - i] * recurrence[i]
            for i in range(order - j - 1)
        )
        weights.append(_normalize(value, modulus))

    # 1-based X dictionary for the final linear form and optional short tail.
    X = {i + 1: _normalize(x[i], modulus) for i in range(order)}

    if term_count > 1:
        for ell in range(term_count - 1):
            idx = order + 1 + ell
            value = sum(
                X[order + ell - i] * recurrence[i]
                for i in range(order)
            )
            X[idx] = _normalize(value, modulus)

    def term_at(r):
        return _normalize(
            sum(weights[j] * X[j + r] for j in range(order)), modulus
        )

    if term_count == 1:
        return term_at(1)
    return [term_at(r) for r in range(1, term_count + 1)]


# ---------------------------------------------------------------------------
# Sage-dependent Toeplitz transfer construction
# ---------------------------------------------------------------------------


def _require_sage():
    """Import Sage lazily and give a targeted diagnostic under plain Python."""
    try:
        import sage.all as sage
    except ImportError as exc:  # pragma: no cover - depends on host environment
        raise RuntimeError(
            "This operation requires SageMath. Run it with `sage` or `sage -python`."
        ) from exc
    return sage


def _validate_bandwidth(m1, m2):
    m1 = _as_int(m1, "m1")
    m2 = _as_int(m2, "m2")
    if m1 < 0 or m2 < 0:
        raise ValueError("m1 and m2 must be nonnegative")
    return m1, m2


def _diag_name(offset: int) -> str:
    if offset < 0:
        return f"a_m{-offset}"
    if offset > 0:
        return f"a_p{offset}"
    return "a_0"


def _default_diagonal_values(m1: int, m2: int):
    sage = _require_sage()
    names = [_diag_name(s) for s in range(-m1, m2 + 1)]
    ring = sage.PolynomialRing(sage.QQ, names=names)
    return list(ring.gens())


def _constant_band_entry_setup(m1: int, m2: int, *, diagonal_values=None):
    """Prepare exact constant diagonal data without using Sage's symbolic ring."""
    values = (
        _default_diagonal_values(m1, m2)
        if diagonal_values is None
        else list(diagonal_values)
    )
    if len(values) != m1 + m2 + 1:
        raise ValueError(
            f"diagonal_values must have length {m1 + m2 + 1}, "
            "in offset order -m1,...,m2"
        )
    assoc = {offset: values[offset + m1] for offset in range(-m1, m2 + 1)}

    def entry(offset, _position):
        return assoc.get(offset, 0)

    return {
        "entry_function": entry,
        "diagonal_values": values,
        "k": 0,
    }


@dataclass(frozen=True)
class _BandTerm:
    """A signed formal band entry a_offset(k+shift) used before ring creation."""

    offset: int
    shift: int
    sign: int = 1

    def __neg__(self):
        return _BandTerm(self.offset, self.shift, -self.sign)

    def __mul__(self, scalar):
        if scalar == 1:
            return self
        if scalar == -1:
            return -self
        return NotImplemented

    __rmul__ = __mul__


def _formal_band_entry(offset, position):
    return _BandTerm(int(offset), _as_int(position, "position"), 1)


def _position_suffix(position: int) -> str:
    if position < 0:
        return f"m{-position}"
    if position > 0:
        return f"p{position}"
    return "0"


def _formal_position_ring(terms, base_position):
    """Create a finite exact polynomial ring for the shifted entries in one Q[k]."""
    sage = _require_sage()
    pairs = sorted({
        (term.offset, base_position + term.shift)
        for term in terms
    })
    if not pairs:
        return sage.QQ, {}
    names = [
        f"{_diag_name(offset)}_at_k_{_position_suffix(position)}"
        for offset, position in pairs
    ]
    ring = sage.PolynomialRing(sage.QQ, names=names)
    return ring, dict(zip(pairs, ring.gens()))


def _instantiate_band_terms(terms, base_position=0, band_entry_function=None):
    """Evaluate signed shifted band terms in an exact Sage parent."""
    base_position = _as_int(base_position, "position")
    if band_entry_function is None:
        _, lookup = _formal_position_ring(terms, base_position)

        def value(term):
            return term.sign * lookup[(term.offset, base_position + term.shift)]
    else:
        def value(term):
            return term.sign * band_entry_function(
                term.offset, base_position + term.shift
            )
    return [value(term) for term in terms]


def _final_signature(state):
    rows, cols = state
    length = len(rows)
    delta = length - rows[0]
    return (
        tuple(v + delta for v in rows),
        tuple(v + delta for v in cols),
    )


def _laplace_expansion(m1, m2, state, quantity, entry, k):
    rows, cols = state
    sf = _final_signature((tuple(rows), tuple(cols)))
    srows, scols = sf
    length = len(srows)
    expansion = []

    if length == 1:
        for j in range(m1 + 1):
            coeff = entry(-j, k - j)
            if quantity == "Determinant":
                coeff *= -1 if j % 2 else 1
            target = m_reduce(srows + (0,), scols + (-j,))
            expansion.append((coeff, target))
        return sf, expansion

    gap = srows[0] - srows[1]
    if gap == 1:
        sign_index = 0
        for j in range(m2 + 1):
            row = srows[0] - 1 - j
            col = scols[0] - 1
            if row not in srows:
                coeff = entry(j, k + row)
                if quantity == "Determinant" and sign_index % 2:
                    coeff = -coeff
                target = m_reduce(srows + (row,), scols + (col,))
                expansion.append((coeff, target))
                sign_index += 1
        return sf, expansion

    if isinstance(gap, int) and gap > 1:
        sign_index = 0
        for j in range(m1 + 1):
            row = srows[0] - 1
            col = scols[0] - 1 - j
            if col not in scols:
                coeff = entry(-j, k + col)
                if quantity == "Determinant" and sign_index % 2:
                    coeff = -coeff
                target = m_reduce(srows + (row,), scols + (col,))
                expansion.append((coeff, target))
                sign_index += 1
        return sf, expansion

    raise RuntimeError(f"boundary signature {state!r} did not match an expansion case")


def _instantiate_position_expansion(data, *, base_position=0, band_entry_function=None):
    """Instantiate one structural Laplace expansion at an integer position."""
    signature, expansion = data
    terms = [term for term, _ in expansion]
    values = _instantiate_band_terms(
        terms, base_position=base_position,
        band_entry_function=band_entry_function,
    )
    return signature, [
        (value, target)
        for value, (_, target) in zip(values, expansion)
    ]


def laplace_d(
    m1,
    m2,
    state,
    *,
    diagonal_values=None,
    position_dependent=False,
    band_entry_function=None,
):
    """One determinant Laplace expansion for a normalized boundary state.

    In position-dependent mode the default formal expansion is anchored at
    integer position 0 and uses an exact multivariate polynomial ring whose
    generators encode the required shifted band entries.  A custom
    ``band_entry_function(offset, position)`` receives integer positions.
    """
    m1, m2 = _validate_bandwidth(m1, m2)
    if position_dependent:
        if diagonal_values is not None:
            raise ValueError(
                "diagonal_values cannot be supplied with position_dependent=True; "
                "use band_entry_function instead"
            )
        data = _laplace_expansion(
            m1, m2, state, "Determinant", _formal_band_entry, 0
        )
        return _instantiate_position_expansion(
            data, base_position=0, band_entry_function=band_entry_function
        )

    setup = _constant_band_entry_setup(
        m1, m2, diagonal_values=diagonal_values
    )
    return _laplace_expansion(
        m1, m2, state, "Determinant", setup["entry_function"], 0
    )


def laplace_p(
    m1,
    m2,
    state,
    *,
    diagonal_values=None,
    position_dependent=False,
    band_entry_function=None,
):
    """Signless permanent analogue of :func:`laplace_d`."""
    m1, m2 = _validate_bandwidth(m1, m2)
    if position_dependent:
        if diagonal_values is not None:
            raise ValueError(
                "diagonal_values cannot be supplied with position_dependent=True; "
                "use band_entry_function instead"
            )
        data = _laplace_expansion(
            m1, m2, state, "Permanent", _formal_band_entry, 0
        )
        return _instantiate_position_expansion(
            data, base_position=0, band_entry_function=band_entry_function
        )

    setup = _constant_band_entry_setup(
        m1, m2, diagonal_values=diagonal_values
    )
    return _laplace_expansion(
        m1, m2, state, "Permanent", setup["entry_function"], 0
    )


def _row_column_closure_data(m1, m2, quantity, entry, k):
    """Return closed row-column state data before choosing a coefficient parent."""
    seen = [((1,), (1,))]
    answer = []
    max_states = 2 * comb(m1 + m2, m1) + m1 + m2 + 10
    state_index = 0

    while state_index < len(seen):
        if state_index >= max_states:
            raise RuntimeError("row-column state closure exceeded its safety bound")
        data = _laplace_expansion(
            m1, m2, seen[state_index], quantity, entry, k
        )
        answer.append(data)
        for _, target in data[1]:
            tp = m_reduce(*_final_signature(target))
            if tp not in seen:
                seen.append(tp)
        state_index += 1

    finals = [item[0] for item in answer]
    shifted_finals = [
        (tuple(v - 1 for v in rows), tuple(v - 1 for v in cols))
        for rows, cols in finals
    ]
    order = len(answer)
    entries = {}
    for i, (_, expansion) in enumerate(answer):
        for coeff, target in expansion:
            try:
                j = shifted_finals.index(target)
            except ValueError as exc:
                raise RuntimeError(
                    f"Laplace target {target!r} did not match a closed state"
                ) from exc
            key = (i, j)
            if key in entries:
                raise RuntimeError(
                    f"two Laplace terms collided in transition {i + 1}->{j + 1}"
                )
            entries[key] = coeff

    return {
        "expansion_data": answer,
        "states": finals,
        "shifted_states": shifted_finals,
        "state_order": order,
        "state_ordering": "Paper first-discovery order",
        "entries": entries,
    }


def _matrix_from_entry_dict(order, entries):
    sage = _require_sage()
    return sage.matrix(
        nrows=order, ncols=order, entries=entries, sparse=True
    )


def _position_matrix_from_template(
    template, base_position=0, band_entry_function=None
):
    """Instantiate an exact position-dependent transfer matrix Q[base_position]."""
    base_position = _as_int(base_position, "position")
    entries = template["entries"]
    terms = list(entries.values())
    if band_entry_function is None:
        _, lookup = _formal_position_ring(terms, base_position)

        def evaluate(term):
            return term.sign * lookup[
                (term.offset, base_position + term.shift)
            ]
    else:
        def evaluate(term):
            return term.sign * band_entry_function(
                term.offset, base_position + term.shift
            )

    instantiated = {key: evaluate(term) for key, term in entries.items()}
    return _matrix_from_entry_dict(template["state_order"], instantiated)


def _validate_scalar_policy(policy):
    if policy in (True, False):
        return None
    if isinstance(policy, str) and policy.lower() in {"auto", "automatic"}:
        return None
    raise ValueError("scalar_recurrence must be True, False, or 'auto'")


def _should_scalarize(policy, order, cutoff):
    _validate_scalar_policy(policy)
    if policy in (True, False):
        return policy
    cutoff = _as_int(cutoff, "max_characteristic_order")
    if cutoff < 0:
        raise ValueError("max_characteristic_order must be nonnegative")
    return order <= cutoff


def _normalize_scalar_method(method):
    key = str(method).replace("-", "_").replace(" ", "_").lower()
    if key in {"characteristicpolynomial", "characteristic_polynomial"}:
        return "characteristic_polynomial"
    if key == "krylov":
        return "krylov"
    raise ValueError(
        "scalar_recurrence_method must be 'characteristic_polynomial' or 'krylov'"
    )


def _row_column_engine(
    quantity,
    m1,
    m2,
    *,
    diagonal_values=None,
    scalar_recurrence="auto",
    characteristic_polynomial_variable="x",
    max_characteristic_order=20,
    verbose=True,
    position_dependent=False,
    band_entry_function=None,
):
    m1, m2 = _validate_bandwidth(m1, m2)
    _validate_scalar_policy(scalar_recurrence)

    if position_dependent:
        if diagonal_values is not None:
            raise ValueError(
                "diagonal_values cannot be supplied with position_dependent=True; "
                "use band_entry_function instead"
            )
        closure = _row_column_closure_data(
            m1, m2, quantity, _formal_band_entry, 0
        )

        def transfer_function(position):
            return _position_matrix_from_template(
                closure,
                base_position=position,
                band_entry_function=band_entry_function,
            )

        matrix_expr = transfer_function(0)
        transfer = transfer_function
        characteristic = None
        diagonal_values_out = None
    else:
        setup = _constant_band_entry_setup(
            m1, m2, diagonal_values=diagonal_values
        )
        closure = _row_column_closure_data(
            m1, m2, quantity, setup["entry_function"], 0
        )
        matrix_expr = _matrix_from_entry_dict(
            closure["state_order"], closure["entries"]
        )
        should = _should_scalarize(
            scalar_recurrence, closure["state_order"], max_characteristic_order
        )
        characteristic = (
            matrix_expr.charpoly(str(characteristic_polynomial_variable))
            if should
            else None
        )
        transfer_function = None
        transfer = matrix_expr
        diagonal_values_out = setup["diagonal_values"]

    order = closure["state_order"]
    result = {
        "method": "RowColumn",
        "quantity": quantity,
        "m1": m1,
        "m2": m2,
        "offset_convention": "Paper: offsets -m1,...,m2",
        "states": closure["states"],
        "shifted_states": closure["shifted_states"],
        "state_order": order,
        "state_ordering": closure["state_ordering"],
        "scalar_recurrence_method": "characteristic_polynomial",
        "transfer_matrix": transfer,
        "transfer_matrix_expression": matrix_expr,
        "transfer_matrix_function": transfer_function,
        "diagonal_values": diagonal_values_out,
        "position_dependent": bool(position_dependent),
        "position_origin": 0 if position_dependent else None,
        "position_variable_convention": (
            "Exact polynomial generators a_<offset>_at_k_<integer>; "
            "transfer_matrix_function(r) returns Q[r]"
            if position_dependent and band_entry_function is None
            else None
        ),
        "characteristic_polynomial_variable": str(characteristic_polynomial_variable),
        "characteristic_polynomial": characteristic,
    }

    if verbose:
        print(f"The transfer-matrix order is {order}.")
        print(
            "The exact position-dependent transfer matrix Q[0] is:"
            if position_dependent
            else "The transfer matrix is:"
        )
        print(matrix_expr)
        if not position_dependent:
            if characteristic is None:
                print("The characteristic polynomial was not computed.")
            else:
                print("The characteristic polynomial is:")
                print(characteristic)
    return result


def lin_rec_for_dtm(
    m1,
    m2,
    *,
    diagonal_values=None,
    scalar_recurrence="auto",
    characteristic_polynomial_variable="x",
    max_characteristic_order=20,
    verbose=True,
    position_dependent=False,
    band_entry_function=None,
    scalar_recurrence_method="characteristic_polynomial",
):
    """Construct the determinant row-column transfer system."""
    method = _normalize_scalar_method(scalar_recurrence_method)
    if method == "krylov":
        return _krylov_wrapper(
            "Determinant", m1, m2,
            diagonal_values=diagonal_values,
            scalar_recurrence=scalar_recurrence,
            characteristic_polynomial_variable=characteristic_polynomial_variable,
            max_characteristic_order=max_characteristic_order,
            verbose=verbose,
            position_dependent=position_dependent,
            band_entry_function=band_entry_function,
        )
    return _row_column_engine(
        "Determinant",
        m1,
        m2,
        diagonal_values=diagonal_values,
        scalar_recurrence=scalar_recurrence,
        characteristic_polynomial_variable=characteristic_polynomial_variable,
        max_characteristic_order=max_characteristic_order,
        verbose=verbose,
        position_dependent=position_dependent,
        band_entry_function=band_entry_function,
    )


def lin_rec_for_ptm(
    m1,
    m2,
    *,
    diagonal_values=None,
    scalar_recurrence="auto",
    characteristic_polynomial_variable="x",
    max_characteristic_order=20,
    verbose=True,
    position_dependent=False,
    band_entry_function=None,
    scalar_recurrence_method="characteristic_polynomial",
):
    """Construct the permanent row-column transfer system."""
    method = _normalize_scalar_method(scalar_recurrence_method)
    if method == "krylov":
        return _krylov_wrapper(
            "Permanent", m1, m2,
            diagonal_values=diagonal_values,
            scalar_recurrence=scalar_recurrence,
            characteristic_polynomial_variable=characteristic_polynomial_variable,
            max_characteristic_order=max_characteristic_order,
            verbose=verbose,
            position_dependent=position_dependent,
            band_entry_function=band_entry_function,
        )
    return _row_column_engine(
        "Permanent",
        m1,
        m2,
        diagonal_values=diagonal_values,
        scalar_recurrence=scalar_recurrence,
        characteristic_polynomial_variable=characteristic_polynomial_variable,
        max_characteristic_order=max_characteristic_order,
        verbose=verbose,
        position_dependent=position_dependent,
        band_entry_function=band_entry_function,
    )


# ---------------------------------------------------------------------------
# Observable/Krylov scalarization
# ---------------------------------------------------------------------------


class _IncrementalRowBasis:
    """Incremental fraction-free row elimination with certificates.

    Accepted rows are kept in row-echelon form.  To eliminate a pivot ``p``
    from a new row ``w`` using a basis row ``b``, the update is
    ``b[p] * w - w[p] * b``.  This avoids introducing rational-function
    denominators at every Krylov step when the transfer entries are
    polynomial.

    ``add(row)`` returns ``None`` when the row increases the rank.  When the
    row is dependent, it returns normalized coefficients ``r`` in the
    coordinates of all rows supplied so far such that
    ``sum(r[i] * row_i) == 0`` and the newest coefficient is one.
    """

    def __init__(self, width, zero, one):
        self.width = int(width)
        self.zero = zero
        self.one = one
        self.pivots = []
        self.rows = []
        self.representations = []
        self.row_count = 0

    @property
    def rank(self):
        return len(self.pivots)

    def add(self, row):
        work = list(row)
        if len(work) != self.width:
            raise ValueError("row has the wrong width")

        for rep in self.representations:
            rep.append(self.zero)

        rep = [self.zero] * self.row_count + [self.one]
        self.row_count += 1

        # Fraction-free elimination against the current echelon basis.
        for pivot, basis_row, basis_rep in zip(
            self.pivots, self.rows, self.representations
        ):
            factor = work[pivot]
            if factor != self.zero:
                pivot_value = basis_row[pivot]
                work = [
                    pivot_value * a - factor * b
                    for a, b in zip(work, basis_row)
                ]
                rep = [
                    pivot_value * a - factor * b
                    for a, b in zip(rep, basis_rep)
                ]

        pivot = next(
            (j for j, value in enumerate(work) if value != self.zero),
            None,
        )
        if pivot is None:
            newest = rep[-1]
            return [value / newest for value in rep]

        insert_at = 0
        while insert_at < len(self.pivots) and self.pivots[insert_at] < pivot:
            insert_at += 1
        self.pivots.insert(insert_at, pivot)
        self.rows.insert(insert_at, work)
        self.representations.insert(insert_at, rep)
        return None


def _principal_coordinate(states):
    for i, state in enumerate(states):
        rows, cols = state
        if len(rows) == 1 and len(cols) == 1:
            return i
    return 0


def _fraction_field_for_matrix(mat):
    ring = mat.base_ring()
    try:
        if ring.is_field():
            return ring
    except Exception:
        pass
    try:
        return ring.fraction_field()
    except Exception:
        return ring


def _krylov_scalarization(matrix_obj, variable, coordinate):
    sage = _require_sage()
    n = matrix_obj.nrows()
    ring = matrix_obj.base_ring()
    work_matrix = matrix_obj

    row = sage.vector(ring, [1 if i == coordinate else 0 for i in range(n)])
    rows = [row]
    reducer = _IncrementalRowBasis(n, ring.zero(), ring.one())
    reducer.add(row)

    for _ in range(n):
        nxt = rows[-1] * work_matrix
        relation = reducer.add(nxt)
        if relation is not None:
            field = _fraction_field_for_matrix(matrix_obj)
            relation = [field(value) for value in relation]
            coeffs = [
                _simplify_fraction(-relation[j] / relation[-1])
                for j in range(len(relation) - 1)
            ]
            r = len(coeffs)
            poly_ring = sage.PolynomialRing(field, str(variable))
            x = poly_ring.gen()
            poly = x**r - sum(poly_ring(coeffs[j]) * x**j for j in range(r))
            observable = sage.matrix(field, rows)
            companion_entries = {}
            for j in range(r - 1):
                companion_entries[(j, j + 1)] = field.one()
            for j, coeff in enumerate(coeffs):
                if coeff != 0:
                    companion_entries[(r - 1, j)] = coeff
            companion = sage.matrix(
                field, r, r, companion_entries, sparse=True
            )
            similarity = observable if r == n else None
            return {
                "observable_order": r,
                "observable_coordinate": coordinate,
                "observable_matrix": observable,
                "krylov_relation_coefficients": coeffs,
                "scalar_recurrence_polynomial": poly,
                "krylov_companion_matrix": companion,
                "similarity_matrix": similarity,
            }
        rows.append(nxt)

    raise RuntimeError("no usable Krylov relation found")


def _krylov_wrapper(
    quantity,
    m1,
    m2,
    *,
    diagonal_values=None,
    scalar_recurrence="auto",
    characteristic_polynomial_variable="x",
    max_characteristic_order=20,
    verbose=True,
    position_dependent=False,
    band_entry_function=None,
):
    base = _row_column_engine(
        quantity,
        m1,
        m2,
        diagonal_values=diagonal_values,
        scalar_recurrence=False,
        characteristic_polynomial_variable=characteristic_polynomial_variable,
        max_characteristic_order=max_characteristic_order,
        verbose=False,
        position_dependent=position_dependent,
        band_entry_function=band_entry_function,
    )
    order = base["state_order"]
    should = _should_scalarize(
        scalar_recurrence, order, max_characteristic_order
    )
    coordinate = _principal_coordinate(base["states"])

    if position_dependent:
        result = dict(base)
        result.update({
            "scalar_recurrence_method": "krylov",
            "characteristic_polynomial": None,
            "scalar_recurrence_polynomial": None,
            "observable_order": None,
            "observable_coordinate": None,
            "observable_matrix": None,
            "krylov_relation_coefficients": None,
            "krylov_companion_matrix": None,
            "similarity_matrix": None,
        })
        if verbose:
            print(f"The transfer-matrix order is {order}.")
            print("Krylov scalarization was not performed because the transfer matrix depends on k.")
        return result

    if not should:
        result = dict(base)
        result.update({
            "scalar_recurrence_method": "krylov",
            "characteristic_polynomial": None,
            "scalar_recurrence_polynomial": None,
            "observable_order": None,
            "observable_coordinate": coordinate,
            "observable_matrix": None,
            "krylov_relation_coefficients": None,
            "krylov_companion_matrix": None,
            "similarity_matrix": None,
        })
        if verbose:
            print(f"The transfer-matrix order is {order}.")
            print("The Krylov scalar recurrence was not computed.")
        return result

    data = _krylov_scalarization(
        base["transfer_matrix_expression"],
        characteristic_polynomial_variable,
        coordinate,
    )
    result = dict(base)
    result.update(data)
    result["scalar_recurrence_method"] = "krylov"
    result["characteristic_polynomial"] = None

    if verbose:
        print(f"The transfer-matrix order is {order}.")
        print("The transfer matrix is:")
        print(base["transfer_matrix_expression"])
        print(f"The observable Krylov order is {result['observable_order']}.")
        print("The scalar recurrence polynomial from the determinant/permanent coordinate is:")
        print(result["scalar_recurrence_polynomial"])
    return result


def _permanent_sage(mat):
    sage = _require_sage()
    from sage.matrix.matrix_misc import permanental_minor_polynomial

    if mat.nrows() == 0:
        return sage.ZZ.one()
    return permanental_minor_polynomial(mat, permanent_only=True)


def _minor_value(rows, cols, quantity, entry):
    sage = _require_sage()
    if len(rows) != len(cols):
        raise RuntimeError(
            f"coefficient-minor row/column sizes differ: {len(rows)} vs {len(cols)}"
        )
    if not rows:
        return sage.ZZ.one()
    data = [[entry(cols[j] - rows[i], 0) for j in range(len(cols))]
            for i in range(len(rows))]
    mat = sage.matrix(data)
    return mat.det() if quantity == "Determinant" else _permanent_sage(mat)


def _boundary_states(m1, m2):
    universe = range(m1 - m2, 2 * m1)
    return [tuple(reversed(c)) for c in combinations(universe, m2)]


def _integer_range(a, b):
    """Inclusive integer range matching Mathematica Range[a,b]."""
    if b < a:
        return []
    return list(range(a, b + 1))


def _increasing_coefficient(m1, m2, ell, state, quantity, entry, k0):
    if ell >= 0:
        r_ell = _integer_range(k0 + 1 - m1, k0 + ell)
        initial = [k0 - s for s in state]
        tail = _integer_range(k0 + 1 + m2 - m1, k0 + ell)
        cell_s = sorted(set(initial).union(tail))
        return _minor_value(r_ell, cell_s, quantity, entry)

    h = [k0 - s for s in state]
    universe = set(_integer_range(1, k0 + m2 - m1))
    jset = sorted(universe.difference(h))
    if set(jset).difference(_integer_range(1, k0 + ell)):
        return 0
    r_ell = _integer_range(k0 + 1 - m1, k0 + ell)
    cell_s = sorted(set(_integer_range(1, k0 + ell)).difference(jset))
    return _minor_value(r_ell, cell_s, quantity, entry)


def _simplify_fraction(value):
    for method_name in ("cancel", "simplify_rational", "simplify_full"):
        method = getattr(value, method_name, None)
        if callable(method):
            try:
                return method()
            except Exception:
                pass
    return value


def _kernel_over_fraction_field(mat):
    ring = mat.base_ring()
    try:
        is_field = bool(ring.is_field())
    except Exception:
        is_field = False
    work = mat
    if not is_field:
        try:
            work = mat.change_ring(ring.fraction_field())
        except Exception:
            pass
    return work.right_kernel().basis()


def _normalize_kernel_vector(v):
    nz = [i for i, value in enumerate(v) if value != 0]
    if not nz:
        return list(v)
    lead = v[nz[-1]]
    return [_simplify_fraction(value / lead) for value in v]


def _increasing_rows_engine(
    quantity,
    m1_in,
    m2_in,
    *,
    diagonal_values=None,
    ell_range="forward",
    characteristic_polynomial_variable="x",
    verbose=True,
):
    sage = _require_sage()
    m1_in, m2_in = _validate_bandwidth(m1_in, m2_in)
    setup = _constant_band_entry_setup(
        m1_in, m2_in, diagonal_values=diagonal_values
    )
    values = setup["diagonal_values"]

    m1, m2 = m1_in, m2_in
    transposed = False
    if m1 < m2:
        transposed = True
        original = {s: values[s + m1] for s in range(-m1, m2 + 1)}
        m1, m2 = m2, m1
        internal_values = [original.get(-s, 0) for s in range(-m1, m2 + 1)]
    else:
        internal_values = list(values)

    assoc = {s: internal_values[s + m1] for s in range(-m1, m2 + 1)}

    def entry(s, t):
        return assoc.get(s, 0)

    d = comb(m1 + m2, m1)
    states = _boundary_states(m1, m2)

    if isinstance(ell_range, str):
        key = ell_range.lower()
        if key == "forward":
            ell_values = list(range(0, d + 1))
        elif key == "centered":
            ell_values = list(range(-m1, d - m1 + 1))
        else:
            raise ValueError("ell_range must be 'forward', 'centered', or an explicit list")
    else:
        ell_values = [_as_int(v, "ell_range entry") for v in ell_range]
        if (
            len(ell_values) != d + 1
            or len(set(ell_values)) != len(ell_values)
            or ell_values != sorted(ell_values)
            or min(ell_values) < -m1
        ):
            raise ValueError(
                f"explicit ell_range must be a strictly increasing list of {d + 1} "
                f"integers with minimum at least {-m1}"
            )

    k0 = 4 * (m1 + m2 + d + 2)
    b_rows = [
        [
            _increasing_coefficient(m1, m2, ell, state, quantity, entry, k0)
            for state in states
        ]
        for ell in ell_values
    ]
    B = sage.matrix(b_rows)
    basis = _kernel_over_fraction_field(B.transpose())
    if not basis:
        raise RuntimeError("no nonzero left-kernel vector found")
    c = _normalize_kernel_vector(basis[0])

    coefficient_ring = c[0].parent() if hasattr(c[0], "parent") else sage.QQ
    try:
        poly_ring = sage.PolynomialRing(
            coefficient_ring, str(characteristic_polynomial_variable)
        )
        x = poly_ring.gen()
        poly = sum(
            poly_ring(c[i]) * x ** (ell_values[i] - min(ell_values))
            for i in range(len(ell_values))
        )
    except Exception as exc:
        raise TypeError(
            "the coefficient parent does not support an exact polynomial ring; "
            "supply diagonal_values in an exact Sage ring"
        ) from exc

    orders = [m1 + ell for ell in ell_values]
    result = {
        "method": "IncreasingRows",
        "quantity": quantity,
        "m1": m1_in,
        "m2": m2_in,
        "internal_m1": m1,
        "internal_m2": m2,
        "transposed": transposed,
        "offset_convention": "Paper: offsets -m1,...,m2",
        "diagonal_values": values,
        "internal_diagonal_values": internal_values,
        "boundary_states": states,
        "ell_values": ell_values,
        "coefficient_minor_orders": orders,
        "coefficient_matrix": B,
        "kernel_vector": c,
        "characteristic_polynomial_variable": str(characteristic_polynomial_variable),
        "annihilating_polynomial": poly,
    }
    if verbose:
        print(f"The increasing-rows boundary-state count is {d}.")
        print(f"The ell values are {ell_values}.")
        print(f"The coefficient-minor orders are {orders}.")
        print("An annihilating polynomial is:")
        print(poly)
    return result


def lin_rec_for_dtm_increasing_rows(
    m1,
    m2,
    *,
    diagonal_values=None,
    ell_range="forward",
    characteristic_polynomial_variable="x",
    verbose=True,
):
    """Algorithm 1 (increasing rows) for determinants."""
    return _increasing_rows_engine(
        "Determinant", m1, m2,
        diagonal_values=diagonal_values,
        ell_range=ell_range,
        characteristic_polynomial_variable=characteristic_polynomial_variable,
        verbose=verbose,
    )


def lin_rec_for_ptm_increasing_rows(
    m1,
    m2,
    *,
    diagonal_values=None,
    ell_range="forward",
    characteristic_polynomial_variable="x",
    verbose=True,
):
    """Algorithm 1 (increasing rows) for permanents."""
    return _increasing_rows_engine(
        "Permanent", m1, m2,
        diagonal_values=diagonal_values,
        ell_range=ell_range,
        characteristic_polynomial_variable=characteristic_polynomial_variable,
        verbose=verbose,
    )


# ---------------------------------------------------------------------------
# Remote Toeplitz terms and Sage-native C-finite comparison
# ---------------------------------------------------------------------------


def _toeplitz_matrix(m1, m2, diagonals, size):
    sage = _require_sage()
    return sage.matrix(
        size,
        size,
        lambda i, j: diagonals[j - i + m1] if -m1 <= j - i <= m2 else 0,
    )


def _initial_toeplitz_terms(quantity, m1, m2, diagonals, order, modulus):
    values = [1]
    for size in range(1, order):
        mat = _toeplitz_matrix(m1, m2, diagonals, size)
        value = mat.det() if quantity == "Determinant" else _permanent_sage(mat)
        values.append(_normalize(value, modulus))
    return values


def _recurrence_coefficients(polynomial):
    try:
        degree = int(polynomial.degree())
    except Exception as exc:
        raise ValueError("scalar recurrence polynomial has no usable degree") from exc
    if degree < 1:
        raise ValueError("scalar recurrence polynomial must have positive degree")
    try:
        leading = polynomial.leading_coefficient()
    except Exception:
        leading = polynomial[degree]
    if leading == 0:
        raise ValueError("scalar recurrence polynomial has zero leading coefficient")
    return [
        _simplify_fraction(-polynomial[degree - i] / leading)
        for i in range(1, degree + 1)
    ]


def _term_for_toeplitz(
    quantity,
    m1,
    m2,
    n,
    *,
    diagonal_values=None,
    scalar_recurrence_method="characteristic_polynomial",
    modulus=None,
    term_count=1,
    verbose=False,
    position_dependent=False,
):
    m1, m2 = _validate_bandwidth(m1, m2)
    if position_dependent:
        raise ValueError(
            "Fiduccia remote-term evaluation requires a constant-coefficient recurrence; "
            "position_dependent=True is not supported"
        )
    method = _normalize_scalar_method(scalar_recurrence_method)
    constructor = lin_rec_for_dtm if quantity == "Determinant" else lin_rec_for_ptm
    data = constructor(
        m1,
        m2,
        diagonal_values=diagonal_values,
        scalar_recurrence=True,
        characteristic_polynomial_variable=_fresh_polynomial_variable_name(),
        verbose=verbose,
        position_dependent=False,
        scalar_recurrence_method=method,
    )
    polynomial = (
        data.get("scalar_recurrence_polynomial")
        if method == "krylov"
        else data.get("characteristic_polynomial")
    )
    if polynomial is None:
        raise RuntimeError("no usable constant scalar recurrence polynomial was returned")
    recurrence = _recurrence_coefficients(polynomial)
    order = len(recurrence)
    if order < 2:
        raise ValueError(
            f"returned scalar recurrence has order {order}; Fiduccia requires order at least 2"
        )
    diagonals = data.get("diagonal_values")
    if diagonals is None or len(diagonals) != m1 + m2 + 1:
        raise RuntimeError("constructor did not return usable diagonal values")

    if modulus is not None:
        modulus = _as_int(modulus, "modulus")
        if modulus <= 1:
            raise ValueError("modulus must be greater than 1")
        try:
            diagonals = [_exact_integer_value(v) for v in diagonals]
            recurrence = [_exact_integer_value(v) for v in recurrence]
        except TypeError as exc:
            raise TypeError(
                "with a finite modulus, Toeplitz diagonals and recurrence coefficients must be integers"
            ) from exc

    initials = _initial_toeplitz_terms(
        quantity, m1, m2, diagonals, order, modulus
    )
    return fiduccia_pol_squarings(
        initials,
        recurrence,
        n,
        modulus=modulus,
        term_count=term_count,
    )


def term_for_dtm(
    m1,
    m2,
    n,
    *,
    diagonal_values=None,
    scalar_recurrence_method="characteristic_polynomial",
    modulus=None,
    term_count=1,
    verbose=False,
    position_dependent=False,
):
    """Remote determinant term for a constant banded Toeplitz family."""
    return _term_for_toeplitz(
        "Determinant", m1, m2, n,
        diagonal_values=diagonal_values,
        scalar_recurrence_method=scalar_recurrence_method,
        modulus=modulus,
        term_count=term_count,
        verbose=verbose,
        position_dependent=position_dependent,
    )


def term_for_ptm(
    m1,
    m2,
    n,
    *,
    diagonal_values=None,
    scalar_recurrence_method="characteristic_polynomial",
    modulus=None,
    term_count=1,
    verbose=False,
    position_dependent=False,
):
    """Remote permanent term for a constant banded Toeplitz family."""
    return _term_for_toeplitz(
        "Permanent", m1, m2, n,
        diagonal_values=diagonal_values,
        scalar_recurrence_method=scalar_recurrence_method,
        modulus=modulus,
        term_count=term_count,
        verbose=verbose,
        position_dependent=position_dependent,
    )


def sage_cfinite_term(initial, recurrence, n, *, term_count=1):
    """Evaluate a recurrence with Sage's native :class:`CFiniteSequence`.

    This is intentionally a comparison/reference path, not a replacement for
    :func:`fiduccia_pol_squarings`.  Sage's C-finite implementation accepts
    rational/integer coefficients.  Our recurrence list is reversed before it
    is passed to Sage because Sage uses the convention
    ``a[n+d] = c[0] a[n] + ... + c[d-1] a[n+d-1]``.
    """
    sage = _require_sage()
    initial = list(initial)
    recurrence = list(recurrence)
    if len(initial) != len(recurrence) or len(initial) < 1:
        raise ValueError("initial and recurrence must have the same positive length")
    n = _as_int(n, "n")
    term_count = _as_int(term_count, "term_count")
    if n < 0 or term_count < 1:
        raise ValueError("n must be nonnegative and term_count positive")
    try:
        q_initial = [sage.QQ(v) for v in initial]
        q_recurrence = [sage.QQ(v) for v in recurrence]
    except Exception as exc:
        raise TypeError("sage_cfinite_term requires integer/rational data") from exc
    C = sage.CFiniteSequences(sage.QQ, names=("z",))
    seq = C.from_recurrence(list(reversed(q_recurrence)), q_initial)
    if term_count == 1:
        return seq[n]
    return [seq[n + i] for i in range(term_count)]
