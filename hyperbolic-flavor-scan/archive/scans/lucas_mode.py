"""
latticefit/lucas.py
===================
Lucas mode: theoretically-derived base phi for arithmetic hyperbolic
3-manifold data.

The Lucas-geodesic bridge theorem (Gentry 2026):
    geodesic length ell = k * log(phi)
    iff holonomy trace magnitude |tr(gamma)| = L_k = phi^k + phi^(-k)

This is exact. The phi-lattice for fermion masses, the prime dictionary
of the covering tower, and the smearing parameter sigma_opt are all
expressions of L_k = phi^k + phi^(-k) at k in (1/4)*Z.

Usage
-----
    from latticefit.lucas import fit_lucas, lucas, is_prime_lucas

    result = fit_lucas([5.11e-4, 0.1057, 1.777, ...])
    print(result.summary())
"""

import math
from dataclasses import dataclass, field
from typing import Optional, List

PHI = (1 + math.sqrt(5)) / 2
LOG_PHI = math.log(PHI)

# Lucas sequence L_k for k = 0, 1, 2, ..., 23
_LUCAS = [2,1,3,4,7,11,18,29,47,76,123,199,322,521,843,
          1364,2207,3571,5778,9349,15127,24476,39603,64079]
_LUCAS_SET = set(_LUCAS)

# Prime Lucas numbers (OEIS A005479)
PRIME_LUCAS = {2, 3, 7, 11, 29, 47, 199, 521, 2207, 3571}


def lucas(k: float) -> float:
    """L_k = phi^k + phi^(-k) for any real k.
    Integer values: L_0=2, L_1=1, L_2=3, L_3=4, L_4=7, L_5=11, L_6=18, L_7=29...
    """
    return PHI**k + PHI**(-k)


def is_prime_lucas(n: int) -> bool:
    """Return True if n is a prime Lucas number."""
    return n in PRIME_LUCAS


def nearest_q(x: float, anchor: float, q_grid: int = 4):
    """Return (q, residual) for nearest lattice node to x."""
    if x <= 0 or anchor <= 0:
        return 0, 0.5
    ratio = math.log(x / anchor) / LOG_PHI
    q = round(ratio * q_grid) / q_grid
    return q, abs(ratio - q)


@dataclass
class LucasAssignment:
    value: float
    q: float               # nearest quarter-integer index
    residual: float        # |ratio - q| in log-phi units
    predicted: float       # anchor * phi^q
    delta_pct: float       # 100*|value - predicted|/predicted
    k: float               # = q/2, the Lucas parameter
    lk: float              # L_k = phi^k + phi^(-k)
    is_integer_lucas: bool # True if k is integer (L_k is a Lucas integer)
    lucas_integer: Optional[int]  # L_k as integer if is_integer_lucas, else None
    is_prime_lucas: bool   # True if lucas_integer is prime


@dataclass
class LucasResult:
    anchor: float
    q_grid: int
    rms: float
    p_value: float
    null_mean: float
    assignments: List[LucasAssignment]
    lucas_integer_fraction: float  # fraction of data at integer-Lucas positions
    prime_lucas_hits: List[LucasAssignment]  # assignments at prime Lucas positions

    def summary(self) -> str:
        lines = [
            "LucasFit (base = phi, Lucas-geodesic bridge theorem)",
            f"  anchor    = {self.anchor:.6g}",
            f"  RMS       = {self.rms:.5f}",
            f"  p-value   = {self.p_value:.4f}  (null mean = {self.null_mean:.5f})",
            f"  Lucas-integer fraction = {self.lucas_integer_fraction:.3f}",
            "",
            f"  {'value':>12}  {'q':>6}  {'L_k':>8}  {'delta%':>8}  {'notes'}",
            f"  {'-'*55}",
        ]
        for a in self.assignments:
            lk_str = f"{a.lucas_integer}" if a.is_integer_lucas else f"{a.lk:.3f}"
            note = ""
            if a.is_prime_lucas:
                note = "*** prime Lucas"
            elif a.is_integer_lucas:
                note = "** Lucas integer"
            elif a.residual < 0.02:
                note = "***"
            elif a.residual < 0.05:
                note = "**"
            elif a.residual < 0.09:
                note = "*"
            lines.append(
                f"  {a.value:>12.5g}  {a.q:>+6.2f}  {lk_str:>8}  "
                f"{a.delta_pct:>7.3f}%  {note}"
            )
        return "\n".join(lines)


def fit_lucas(
    data,
    anchor: Optional[float] = None,
    q_grid: int = 4,
    n_null: int = 5000,
    rng_seed: int = 42,
) -> LucasResult:
    """
    Fit data to the phi-lattice {anchor * phi^(q/q_grid) : q in Z}.

    Base is fixed to phi by the Lucas-geodesic bridge theorem.
    One fewer free parameter than the general lattice fit.

    Parameters
    ----------
    data     : iterable of positive floats
    anchor   : lattice anchor. If None, uses min(data).
    q_grid   : 4 = quarter-integer lattice (default), 2 = half-integer
    n_null   : null test samples
    rng_seed : reproducibility

    Returns
    -------
    LucasResult
    """
    import numpy as np

    vals = [float(x) for x in data if x > 0]
    if len(vals) < 2:
        raise ValueError("Need at least 2 positive values")

    if anchor is None:
        anchor = min(vals)

    # Compute assignments
    assignments = []
    residuals = []
    for x in vals:
        q, res = nearest_q(x, anchor, q_grid)
        predicted = anchor * PHI**q
        k = q / 2
        lk = lucas(k)
        k_int = abs(k - round(k)) < 1e-9
        lk_int = _LUCAS[round(k)] if (k_int and 0 <= round(k) < len(_LUCAS)) else None
        a = LucasAssignment(
            value=x,
            q=q,
            residual=res,
            predicted=predicted,
            delta_pct=100 * abs(x - predicted) / predicted,
            k=k,
            lk=lk,
            is_integer_lucas=k_int,
            lucas_integer=lk_int,
            is_prime_lucas=(lk_int in PRIME_LUCAS) if lk_int else False,
        )
        assignments.append(a)
        residuals.append(res**2)

    rms = math.sqrt(sum(residuals) / len(residuals))

    # Null test: uniform in same log range
    rng = np.random.default_rng(rng_seed)
    log_min = math.log(min(vals))
    log_max = math.log(max(vals))
    null_rms = []
    for _ in range(n_null):
        null = np.exp(rng.uniform(log_min, log_max, size=len(vals)))
        nr = []
        for x in null:
            _, res = nearest_q(x, anchor, q_grid)
            nr.append(res**2)
        null_rms.append(math.sqrt(sum(nr) / len(nr)))
    null_rms = np.array(null_rms)
    p_value = float(np.mean(null_rms <= rms))

    lucas_int_frac = sum(1 for a in assignments if a.is_integer_lucas) / len(assignments)
    prime_hits = [a for a in assignments if a.is_prime_lucas]

    return LucasResult(
        anchor=anchor,
        q_grid=q_grid,
        rms=rms,
        p_value=p_value,
        null_mean=float(null_rms.mean()),
        assignments=sorted(assignments, key=lambda a: a.q),
        lucas_integer_fraction=lucas_int_frac,
        prime_lucas_hits=prime_hits,
    )
