"""C*-equivariant reduction for the whole family of weight systems (1,-1,-m).

Generalizes jcqft.reduction (the Alpöge–Mathew case m = 2) to arbitrary
integer m >= 1.  For the source action

    lam . (x, y, z) = (lam*x, y/lam, z/lam^m),

the invariant ring is the FREE polynomial ring C[w, v] with

    w = x*y,   v = x^m*z.

A component of C*-weight d is exactly x^d * S(w, v).  Invertibility of the
linear part DF(0) forces the three component weights to be a permutation of
the source weights {1, -1, -m} (a linear monomial x/y/z has weight 1/-1/-m,
and each variable must appear linearly in some component); we fix the order

    F = ( P(w,v)/x^m,  Q(w,v)/x,  x*R(w,v) ),

so d = (-m, -1, 1), sum(d) = -m.

REDUCED KELLER IDENTITY (exact, proved symbolically in
scripts/reduction_113.py).  In coordinates (x, w, v),
det d(x,w,v)/d(x,y,z) = x^{m+1}, and

    det DF = x^{sum(d) + m + 1} * det M = det M,
    M = [[ d_i S_i,  dS_i/dw,  dS_i/dv ]],   (S_1,S_2,S_3) = (P,Q,R),

a function of (w, v) alone.  Expanding det M and comparing with the
two-variable Jacobian J2(A,B) = A_w B_v - A_v B_w gives the compact form

    J2(P * R^m, Q * R)  =  kappa * R^m        <=>      det DF = kappa.

(m = 2 is jcqft.reduction's identity J2(P R^2, Q R) = kappa R^2.)

POLYNOMIALITY.  x^d * w^j * v^k = x^{d + j + m k} y^j z^k, so

    P:  j + m*k >= m,     Q:  j + m*k >= 1,     R:  no constraint,

and DF(0) invertible additionally requires v in P (gives z), w in Q
(gives y), 1 in R (gives x).
"""

import sympy as sp

x, y, z = sp.symbols("x y z")
w, v = sp.symbols("w v")


def invariants(m):
    """(w, v) as polynomials in (x, y, z) for the weight system (1,-1,-m)."""
    return x * y, x**m * z


def assemble(P, Q, R, m):
    """The equivariant map F = (P/x^m, Q/x, x*R) in the variables (x,y,z)."""
    sub = {w: x * y, v: x**m * z}
    return (sp.cancel(P.subs(sub) / x**m),
            sp.cancel(Q.subs(sub) / x),
            sp.expand(x * R.subs(sub)))


def extract(F3, m, extra_syms=()):
    """(P, Q, R) of an equivariant map F = (P/x^m, Q/x, x*R)."""
    allowed = {w, v} | set(extra_syms)
    out = []
    for expr, xpow in zip(F3, (m, 1, -1)):
        e = sp.cancel(sp.expand(expr.subs({y: w / x, z: v / x**m})) * x**xpow)
        assert e.free_symbols <= allowed, f"not equivariant: {e}"
        out.append(sp.expand(e))
    return tuple(out)


def j2(A, B):
    """Jacobian of (A, B) with respect to (w, v)."""
    return sp.expand(sp.diff(A, w) * sp.diff(B, v) - sp.diff(A, v) * sp.diff(B, w))


def keller_residual(P, Q, R, kappa, m):
    """Zero iff F = (P/x^m, Q/x, x*R) has det DF = kappa identically."""
    return sp.expand(j2(P * R**m, Q * R) - kappa * R**m)


def det_m(P, Q, R, m):
    """det M: the reduced det DF as a polynomial in (w, v).

    Equal to det DF of the assembled 3D map (proved in
    scripts/reduction_113.py); equal to J2(P R^m, Q R)/R^m off {R = 0}."""
    d = (-m, -1, 1)
    M = sp.Matrix([[di * Si, sp.diff(Si, w), sp.diff(Si, v)]
                   for di, Si in zip(d, (P, Q, R))])
    return sp.expand(M.det())


def monomial_box(component, m, jmax, kmax):
    """Admissible monomials w^j v^k (j <= jmax, k <= kmax) for one component
    of the ansatz, respecting polynomiality of the assembled 3D map."""
    lo = {"P": m, "Q": 1, "R": -10**9}[component]
    return [w**j * v**k for j in range(jmax + 1) for k in range(kmax + 1)
            if j + m * k >= lo]
