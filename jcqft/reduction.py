"""C*-equivariant normal form and the dimensionally-reduced Keller condition.

In the scaling invariants w = x*y, v = x^2*z, every equivariant map of the
Alpöge–Mathew weight system has the form F = (P(w,v)/x^2, Q(w,v)/x, x*R(w,v)),
and the 3D Keller condition det DF = kappa collapses to the two-variable
identity  J2(P*R^2, Q*R) = kappa * R^2.
"""

import sympy as sp

from jcqft.core import F, PHI

x, y, z = PHI
w, v = sp.symbols("w v")

KAPPA = -2


def extract(F3, extra_syms=()):
    """(P, Q, R) of an equivariant map F = (P/x^2, Q/x, x*R)."""
    allowed = {w, v} | set(extra_syms)
    out = []
    for expr, xpow in zip(F3, (2, 1, -1)):
        e = sp.cancel(sp.expand(expr.subs({y: w / x, z: v / x**2})) * x**xpow)
        assert e.free_symbols <= allowed, f"not equivariant: {e}"
        out.append(sp.expand(e))
    return tuple(out)


def j2(A, B):
    """Jacobian of (A, B) with respect to (w, v)."""
    return sp.expand(sp.diff(A, w) * sp.diff(B, v) - sp.diff(A, v) * sp.diff(B, w))


def keller_residual(P, Q, R, kappa):
    """Zero iff F = (P/x^2, Q/x, x*R) has det DF = kappa identically."""
    return sp.expand(j2(P * R**2, Q * R) - kappa * R**2)


P0, Q0, R0 = extract(F)
assert keller_residual(P0, Q0, R0, KAPPA) == 0
