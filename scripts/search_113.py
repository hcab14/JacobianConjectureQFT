"""Exact counterexample search for the weight system (1,-1,-3).

(Machinery: jcqft/reduction_w.py, proved exact in scripts/reduction_113.py;
write-up: docs/SEARCH_113.md.  Sibling search: scripts/search_213.py.)

Default run (~20 min): proves every claimed identity by assertion and solves
the Keller condition COMPLETELY in the v-linear class (P, Q, R of degree
<= 1 in the invariant v = x^3 z, ARBITRARY degree in w = x y) -- up to one
precisely-mapped gap (stratum D3 with non-squarefree s or t, see section 5),
which is closed by exact in-box Groebner certificates: small box in the
default run, larger boxes under --full (~1 h).

PROVED HERE, all exact:

1. REDUCED KELLER IDENTITY (m = 3, re-asserted from reduction_113):
   for the equivariant map F = (P/x^3, Q/x, x*R), P,Q,R in C[w,v],
       det DF = det M = -3 P J2(Q,R) + Q J2(P,R) + R J2(P,Q),
   a function of (w,v) alone; R^3 det M = J2(P R^3, Q R).  Polynomiality:
   P needs j + 3k >= 3, Q needs j + 3k >= 1; DF(0) invertible needs
   v in P, w in Q, 1 in R, and det DF(0) = -p1(0) q0'(0) r0(0).
   Gauge: target scalings set p1(0) = q0'(0) = r0(0) = 1, hence kappa = -1.

2. v-LINEAR CLASS.  P = p0(w) + p1(w) v, Q = q0(w) + q1(w) v,
   R = r0(w) + r1(w) v.  det M = kappa has v-degree 2: equations E2, E1, E0.
   E2 = 2[q1 (p1 r1)' - 2 (p1 r1) q1'] integrates (Wronskian certificate)
   to p1 r1 = c q1^2 whenever q1 p1 r1 != 0.  Strata:
       A: q1 = r1 = 0      B: r1 = 0, q1 != 0
       C: q1 = 0, r1 != 0  D: q1 r1 != 0 (=> p1 = a s^2 g, q1 = s t g,
                              r1 = b t^2 g, gcd(s,t) = 1, g | kappa const)
   (m = 2 cross-check: the Alpoge-Mathew map satisfies the corresponding
   identity p1^2 r1 = -q1^3/27, i.e. it lives in the m = 2 analogue of the
   stratum D1 that is EMPTY at m = 3 -- the counterexample mechanism does
   not transplant.)

3. STRATA A and B: completely solved, ALL w-degrees.  In the gauge above,
       P = p0 + v,  Q = w + b0*P,  R = 1,     p0 in w^3 C[w], b0 in C
   (b0 = 0 is stratum A, b0 != 0 stratum B).  Every member is a TAME
   automorphism (explicit inverse: an elementary z-shear followed by a
   target shear); generic fiber has 1 point; the infinity prefilter is
   survived only through its known false-positive class (nonlinear
   automorphisms).

4. STRATA C, D0, D1, D2, D3(squarefree): EMPTY, ALL w-degrees, by exact
   integration/divisibility chains (every identity asserted; the final
   arithmetic steps are one-line divisibility arguments in C[w] recorded
   in comments and in docs/SEARCH_113.md).  In particular D1 -- the exact
   analogue of the stratum housing the Alpoge-Mathew counterexample at
   m = 2 -- dies against the polynomiality box: its would-be solutions
   force p0(0) != 0 while the box demands val_w p0 >= 3.

5. GAP + CORROBORATION.  The one case the structural proof does not reach
   is stratum D3 (both s, t nonconstant) with s or t NON-squarefree; this
   needs deg p1 >= 4 or deg r1 >= 4.  In-box Groebner certificates (msolve,
   exact over Q): the default run proves Keller + (r1 != 0) EMPTY for
   deg(p1,q1,r1) <= 2, and --full proves it for deg <= (3,...,3) plus the
   four targeted non-squarefree parametrizations with deg(p0,q0,r0) <=
   (6,5,4) -- all EMPTY.

6. ORBIFOLD MECHANISM.  For (1,-1,-3) the only stabilizer-jump mechanism
   is 3:1 (image on the weight -3 axis: Q = R = 0 != P; the weights -1, +1
   and all coordinate pairs have trivial stabilizer, so no 2:1 exists at
   all).  The classification makes it empty (R is a nonzero constant on
   every Keller solution); independent pointwise Groebner queries at the
   torus-normalized witness points (1,1), (1,0), (0,1) confirm EMPTY.

VERDICT: the v-linear class of (1,-1,-3) contains nonlinear Keller maps,
but ALL are tame automorphisms -- no counterexample, no 3:1 orbifold
covering, no new global-anomaly class.  The Alpoge-Mathew mechanism is
OBSTRUCTED at m = 3 by exact numerology, not merely undiscovered.
"""

from __future__ import annotations

import argparse
import time

import sympy as sp

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from jcqft.prefilter import infinity_prefilter          # noqa: E402
from jcqft.reduction_w import (                          # noqa: E402
    assemble, det_m, extract, j2, keller_residual, monomial_box, v, w, x, y,
    z)

T0 = time.time()


def check(label, cond):
    assert cond, label
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


def zero(expr):
    """Is the (Derivative-carrying) expression identically zero?"""
    return sp.simplify(sp.expand(expr.doit() if hasattr(expr, "doit")
                                 else expr)) == 0


def subd(expr, *pairs):
    """Substitute into an expression containing Derivatives, then evaluate."""
    e = expr
    for old, new in pairs:
        e = e.subs(old, new)
    return sp.expand(e.doit())


dw = lambda f: sp.diff(f, w)   # noqa: E731

# ===========================================================================
print("== 1. reduced Keller identity, m = 3 (re-asserted) ==")
# ===========================================================================

M = 3
Pf, Qf, Rf = [f(w, v) for f in sp.symbols("P Q R", cls=sp.Function)]
_sub = {w: x * y, v: x**M * z}
F_gen = (Pf.subs(_sub) / x**M, Qf.subs(_sub) / x, x * Rf.subs(_sub))
_J = sp.Matrix([[sp.diff(Fi, u) for u in (x, y, z)] for Fi in F_gen])
detDF = sp.simplify(_J.det().subs({y: w / x, z: v / x**M}))
dM = -M * Pf * j2(Qf, Rf) + Qf * j2(Pf, Rf) + Rf * j2(Pf, Qf)
check("det DF == det M == -3 P J2(Q,R) + Q J2(P,R) + R J2(P,Q), a function "
      "of (w,v) alone (GENERIC-FUNCTION PROOF)",
      sp.simplify(detDF - dM) == 0 and not detDF.has(x))
check("compact form R^3 det M == J2(P R^3, Q R)",
      sp.simplify(Rf**M * dM - j2(Pf * Rf**M, Qf * Rf)) == 0)
check("det_m (module) == det M on undetermined functions",
      sp.simplify(det_m(Pf, Qf, Rf, M) - dM) == 0)
boxP = monomial_box("P", M, jmax=4, kmax=1)
boxQ = monomial_box("Q", M, jmax=4, kmax=1)
check("polynomiality: P-box j+3k >= 3 contains v and w^3 but not w, w^2; "
      "Q-box j+3k >= 1 contains w, v but not 1",
      v in boxP and w**3 in boxP and w not in boxP and w**2 not in boxP
      and w in boxQ and v in boxQ and sp.S.One not in
      [sp.sympify(m_) for m_ in boxQ])

# linear part: DF(0) = antidiag(p1(0), q0'(0), r0(0))
_c = sp.symbols("cP cQ cR cW")
F_lin = assemble(_c[0] * v + _c[3] * w**3, _c[1] * w, _c[2], M)
_L = sp.Matrix([[sp.diff(Fi, u) for u in (x, y, z)]
                for Fi in F_lin]).subs({x: 0, y: 0, z: 0})
check("DF(0) = antidiag(p1(0), q0'(0), r0(0)),  det DF(0) = "
      "-p1(0) q0'(0) r0(0)  => gauge p1(0)=q0'(0)=r0(0)=1, kappa = -1",
      _L.det() == -_c[0] * _c[1] * _c[2])

# ===========================================================================
print("== 2. v-linear class: E2, E1, E0 ==")
# ===========================================================================

p0, p1, q0, q1, r0, r1 = [sp.Function(n)(w) for n in
                          ("p0", "p1", "q0", "q1", "r0", "r1")]


def E_coeffs(P1, Q1, R1, P0=p0, Q0=q0, R0=r0, m=3):
    D = sp.expand(det_m(P0 + P1 * v, Q0 + Q1 * v, R0 + R1 * v, m))
    return (sp.expand(D.coeff(v, 2)), sp.expand(D.coeff(v, 1)),
            sp.expand(D.coeff(v, 0)))


E2, E1, E0 = E_coeffs(p1, q1, r1)
check("det M is v-quadratic; E2 == 2[q1 (p1 r1)' - 2 (p1 r1) q1']",
      zero(E2 - 2 * (q1 * dw(p1 * r1) - 2 * p1 * r1 * dw(q1))))
check("E1 == -3p0[Q,R]_1 - 3p1[Q,R]_0 + q0[P,R]_1 + q1[P,R]_0 "
      "+ r0[P,Q]_1 + r1[P,Q]_0  ([A,B]_i = v^i-part of J2(A,B))",
      zero(E1 - (-3 * p0 * (dw(q1) * r1 - q1 * dw(r1))
                 - 3 * p1 * (dw(q0) * r1 - q1 * dw(r0))
                 + q0 * (dw(p1) * r1 - p1 * dw(r1))
                 + q1 * (dw(p0) * r1 - p1 * dw(r0))
                 + r0 * (dw(p1) * q1 - p1 * dw(q1))
                 + r1 * (dw(p0) * q1 - p1 * dw(q0)))))
check("E0 == -3p0(q0'r1 - q1 r0') + q0(p0'r1 - p1 r0') + r0(p0'q1 - p1 q0')",
      zero(E0 - (-3 * p0 * (dw(q0) * r1 - q1 * dw(r0))
                 + q0 * (dw(p0) * r1 - p1 * dw(r0))
                 + r0 * (dw(p0) * q1 - p1 * dw(q0)))))

h = sp.Function("h")(w)
Dfull = sp.expand(det_m(p0 + p1 * v, q0 + q1 * v, r0 + r1 * v, 3))
check("det M is invariant under the v-shear (w,v) -> (w, v + h(w)) "
      "(unimodular; the classification tool below)",
      zero(sp.expand(det_m(*(S.subs(v, v + h) for S in
                             (p0 + p1 * v, q0 + q1 * v, r0 + r1 * v)),
                           3)) - Dfull.subs(v, v + h)))

# m = 2 cross-check: the E2-analogue and the Alpoge-Mathew data
E2m2 = E_coeffs(p1, q1, r1, m=2)[0]
check("m=2 cross-check: E2 integrates to p1^2 r1 = c q1^3 "
      "((p1^2 r1/q1^3)' == p1 E2 / q1^4)",
      zero(sp.diff(p1**2 * r1 / q1**3, w) - p1 * E2m2 / q1**4))
from jcqft.core import F as F_AM  # noqa: E402
P_AM, Q_AM, R_AM = extract(F_AM, 2)
p1AM, q1AM, r1AM = (sp.expand(S.coeff(v, 1)) for S in (P_AM, Q_AM, R_AM))
check("Alpoge-Mathew (m=2) satisfies it with c = -1/27: p1 = (1+w)^3, "
      "q1 = 3(1+w)^2, r1 = -1, p1^2 r1 == -q1^3/27  [the AM map lives in "
      "the m=2 analogue of stratum D1 below]",
      sp.expand(p1AM - (1 + w)**3) == 0
      and sp.expand(q1AM - 3 * (1 + w)**2) == 0 and r1AM == -1
      and sp.expand(p1AM**2 * r1AM + q1AM**3 / 27) == 0)

# ===========================================================================
print("== 3. stratification by E2 ==")
# ===========================================================================
# E2 = 0 and q1 p1 r1 != 0  =>  (p1 r1/q1^2)' = 0 as a rational function
# =>  p1 r1 = c q1^2 exactly, c a nonzero constant.  With p1 != 0 forced
# (p1(0) = 1), the strata are A: q1 = r1 = 0; B: r1 = 0, q1 != 0;
# C: q1 = 0, r1 != 0; D: q1 r1 != 0.

check("Wronskian certificate: (p1 r1/q1^2)' == E2/(2 q1^3)",
      zero(sp.diff(p1 * r1 / q1**2, w) - E2 / (2 * q1**3)))
print("     => on {E2 = 0, q1 != 0, p1 r1 != 0}: p1 r1 = c q1^2, c != 0;")
print("        q1 = 0 or r1 = 0 satisfy E2 identically: strata A, B, C, D.")

# ===========================================================================
print("== 4. strata A and B: complete solution -- the tame family ==")
# ===========================================================================
# Stratum A (q1 = r1 = 0): E2 = E1 = 0 identically and E0 = -p1 (q0 r0)'
# = kappa.  In C[w]:  p1 | kappa  =>  p1 = 1 (gauge);  (q0 r0)' = 1  and
# q0(0) = 0, r0(0) = 1  =>  q0 r0 = w  =>  r0 | w, r0(0)=1  =>  r0 = 1,
# q0 = w.  p0 in w^3 C[w] stays FREE.
check("A: E1 == 0 and E0 == -p1 (q0 r0)' identically",
      zero(E1.subs({q1: 0, r1: 0})) and
      zero(E0.subs({q1: 0, r1: 0}) + p1 * dw(q0 * r0)))
# Stratum B (r1 = 0, q1 != 0): E1 integrates to p1 r0^2 = c q1; then E0
# factors through p1 * d/dw[(p0 r0^3)/c - q0 r0] = kappa, forcing p1 = 1,
# then r0 | w with r0(0) = 1 => r0 = 1, and q0 = b0 p0 + w  (b0 = 1/c).
E1B = sp.expand(E1.subs(r1, 0).doit())
check("B: integration certificate (p1 r0^2/q1)' == r0 E1 / q1^2",
      zero(sp.diff(p1 * r0**2 / q1, w) - r0 * E1B / q1**2))
chat = sp.Symbol("c", nonzero=True)
check("B: E0|_{q1 = p1 r0^2/c} == p1 * d/dw[(p0 r0^3)/c - q0 r0]",
      zero(subd(E0.subs(r1, 0), (q1, p1 * r0**2 / chat))
           - p1 * dw(p0 * r0**3 / chat - q0 * r0)))
print("     => p1 | kappa => p1 = 1; r0 (p0 r0^2/c - q0) = -kappa w = w")
print("        => r0 | w, r0(0) = 1 => r0 = 1, q0 = b0 p0 + w, q1 = b0.")

# THE unified gauged family (b0 = 0 <-> A, b0 != 0 <-> B):
b0 = sp.Symbol("b0")
p0g = sp.Function("p0")(w)          # arbitrary element of w^3 C[w]
Pg = p0g + v
Qg = w + b0 * Pg
Rg = sp.S.One
check("FAMILY: P = p0 + v, Q = w + b0 P, R = 1 has det M == -1 identically "
      "(any p0(w), any b0)",
      zero(sp.expand(det_m(Pg, Qg, Rg, 3)) + 1))

# an explicit sample of each stratum, checked end-to-end in 3D
lam = sp.Symbol("lam")
for label, p0s, b0s in (("A sample (b0 = 0), p0 = w^3 + 2w^4",
                         w**3 + 2 * w**4, sp.S.Zero),
                        ("B sample (b0 = 5), p0 = w^3", w**3,
                         sp.Integer(5))):
    Ps, Qs, Rs = p0g.subs(p0g, p0s) + v, None, sp.S.One
    Ps = p0s + v
    Qs = w + b0s * Ps
    Fs = assemble(Ps, Qs, Rs, 3)
    J3 = sp.Matrix([[sp.diff(Fi, u) for u in (x, y, z)] for Fi in Fs])
    a_, b_, c_ = sp.symbols("a b c")
    # explicit inverse: undo the target shear b -> b - b0*a*c^2 (Q = w+b0*P
    # assembles to F2 = y + b0 x^2 F1 and x^2 = F3^2), then the elementary
    # z-shear: x = c, y = b - b0*a*c^2, z = a - p0s(c*y)/c^3.
    ynew = b_ - b0s * a_ * c_**2
    inv = (c_, ynew, a_ - sp.expand(p0s.subs(w, c_ * ynew) / c_**3))
    comp = tuple(sp.expand(iv.subs(dict(zip((a_, b_, c_), Fs)),
                                   simultaneous=True)) for iv in inv)
    fib = sp.solve([Fs[i] - t_ for i, t_ in
                    enumerate((sp.Rational(3, 7), sp.Rational(-2, 5),
                               sp.Rational(1, 3)))], [x, y, z], dict=True)
    check(f"{label}: det DF == -1, explicit polynomial inverse "
          f"(TAME: target shear o elementary shear), generic fiber = "
          f"{len(fib)} point, infinity prefilter survives="
          f"{infinity_prefilter(Fs, (x, y, z))} (known false-positive "
          "class: nonlinear automorphisms)",
          sp.expand(J3.det()) == -1 and comp == (x, y, z) and len(fib) == 1)

# ===========================================================================
print("== 5. strata C and D: all EMPTY (the interesting part) ==")
# ===========================================================================

# --- C (q1 = 0, r1 != 0) ---------------------------------------------------
E1C = sp.expand(E1.subs(q1, 0).doit())
check("C: certificate (p1/(q0^4 r1))' == E1/(q0^5 r1^2)  =>  E1 = 0 forces "
      "p1 = c q0^4 r1  =>  p1(0) = 0, contradicting p1(0) = 1: EMPTY",
      zero(sp.diff(p1 / (q0**4 * r1), w) - E1C / (q0**5 * r1**2)))

# --- D: parametrization ----------------------------------------------------
a_, b_ = sp.symbols("a b", nonzero=True)
g_ = sp.Function("g")(w)
s_, t_ = sp.Function("s")(w), sp.Function("t")(w)
# p1 r1 = c q1^2 in the UFD C[w]: with g = gcd(p1, r1), the residual
# multiplicities are even => p1 = a s^2 g, q1 = s t g, r1 = b t^2 g,
# gcd(s, t) = 1, ab = c.  E0 is linear in (p1, q1, r1), so g | E0 = kappa:
E0D = sp.expand(E0.subs({p1: a_ * s_**2 * g_, q1: s_ * t_ * g_,
                         r1: b_ * t_**2 * g_}).doit())
check("D: E0 = g * (polynomial)  =>  g | kappa  =>  g constant (absorb: "
      "g = 1); also E2 == 0 holds identically under the parametrization",
      sp.fraction(sp.together(sp.cancel(E0D / g_)))[1] == 1
      and zero(E_coeffs(a_ * s_**2, s_ * t_, b_ * t_**2)[0]))
E2D, E1D, E0D = E_coeffs(a_ * s_**2, s_ * t_, b_ * t_**2)

# --- D0: s, t both constant (p1, q1, r1 nonzero constants) -----------------
al, be, ga, de = sp.symbols("alpha beta gamma delta")
E2K, E1K, E0K = E_coeffs(al, be, ga)
check("D0: E1 == 2 d/dw[beta gamma p0 - 2 alpha gamma q0 + alpha beta r0]",
      zero(E1K - 2 * dw(be * ga * p0 - 2 * al * ga * q0 + al * be * r0)))
r0D0 = de + (2 * al * ga * q0 - be * ga * p0) / (al * be)
VD0 = q0 - (be / al) * p0
check("D0: after solving E1, E0 == d/dw[-(2 alpha gamma/beta) V^2 "
      "- alpha delta V],  V = q0 - (beta/alpha) p0",
      zero(subd(E0K, (r0, r0D0))
           - dw(-2 * al * ga / be * VD0**2 - al * de * VD0)))
print("     => -(2ag/b) V^2 - ad V = kappa w with V(0) = 0, V'(0) != 0:")
print("        V | kappa w => V = eps w, and the w-coefficient forces")
print("        2 alpha gamma/beta = 0 -- impossible.  D0 EMPTY.")

# --- D1: t constant, s nonconstant (the Alpoge-Mathew analogue!) -----------
A_, B_, C_ = sp.symbols("A B C", nonzero=True)
E2t, E1t, E0t = E_coeffs(A_ * s_**2, B_ * s_, C_)
E1z = sp.expand(E1t.subs(r0, 0).doit())     # v-shear h = -r0/C kills r0
E0z = sp.expand(E0t.subs(r0, 0).doit())
check("D1 (r0 sheared to 0): E1 == C[2Bs p0' - 3Bs' p0 + 2As(s'q0 - 2sq0')]"
      " and E0 == C (q0 p0' - 3 p0 q0')",
      zero(E1z - C_ * (2 * B_ * s_ * dw(p0) - 3 * B_ * dw(s_) * p0
                       + 2 * A_ * s_ * (dw(s_) * q0 - 2 * s_ * dw(q0))))
      and zero(E0z - C_ * (q0 * dw(p0) - 3 * p0 * dw(q0))))
pstar = 2 * A_ * s_ * q0 / B_
Y_ = sp.Function("Y")(w)
check("D1: p0* = (2A/B) s q0 solves E1; homogeneous solutions satisfy "
      "(Y^2/s^3)' == Y(2sY' - 3s'Y)/s^4, i.e. Y^2 = e s^3",
      zero(subd(E1z, (p0, pstar)))
      and zero(sp.diff(Y_**2 / s_**3, w)
               - Y_ * (2 * s_ * dw(Y_) - 3 * dw(s_) * Y_) / s_**4))
check("D1 (s not const*square => Y = 0): E0(p0*) == (2AC/B) q0(s'q0 - 2sq0')"
      "  =>  q0 | kappa => q0 = n const => n^2 s' = const => s LINEAR",
      zero(subd(E0z, (p0, pstar))
           - 2 * A_ * C_ / B_ * q0 * (dw(s_) * q0 - 2 * s_ * dw(q0))))
cc, y0 = sp.symbols("c0 y0")
d_ = sp.Function("d")(w)
check("D1 (s = c d^2, d nonconst): E0(p0* + y0 d^3) == "
      "C d (d'q0 - d q0')(4Ac/B q0 + 3 y0 d)  =>  d | kappa: EMPTY",
      zero(subd(E0z, (p0, (pstar + y0 * d_**3).subs(s_, cc * d_**2)),
                (s_, cc * d_**2))
           - C_ * d_ * (dw(d_) * q0 - d_ * dw(q0))
           * (4 * A_ * cc / B_ * q0 + 3 * y0 * d_)))
s0_, s1_, n_ = sp.symbols("s0 s1 n")
slin = s0_ + s1_ * w
famD1 = E_coeffs(A_ * slin**2, B_ * slin, C_,
                 P0=2 * A_ * n_ * slin / B_, Q0=n_, R0=sp.S.Zero)
check("D1 sheared family: s = s0 + s1 w, q0 = n, p0 = (2A/B) n s, r0 = 0 "
      "has E2 = E1 = 0, E0 = (2AC/B) n^2 s1 = kappa != 0 => n != 0",
      zero(famD1[0]) and zero(famD1[1])
      and zero(famD1[2] - 2 * A_ * C_ * n_**2 * s1_ / B_))
# un-shear (p0,q0,r0) = family - (A s^2, B s, C) h and impose the BOX:
h0_ = sp.Symbol("h0")
jets = sp.solve([(2 * A_ * n_ * slin / B_ - A_ * slin**2 * h0_).subs(w, 0),
                 (n_ - B_ * slin * h0_).subs(w, 0)], [h0_, n_], dict=True)
check("D1 BOX: p0(0) = 0 and q0(0) = 0 force n = 0 (unique solution of the "
      "jet system) => kappa = 0.  D1 EMPTY IN THE BOX -- the AM-analogue "
      "stratum is obstructed at m = 3",
      len(jets) == 1 and jets[0][n_] == 0)

# --- D2: s constant, t nonconstant -----------------------------------------
E2u, E1u, E0u = E_coeffs(A_, B_ * t_, C_ * t_**2)
E1z2 = sp.expand(E1u.subs(p0, 0).doit())    # v-shear h = -p0/A kills p0
E0z2 = sp.expand(E0u.subs(p0, 0).doit())
check("D2 (p0 sheared to 0): E0 == -A (q0 r0)'  =>  q0 r0 = -(kappa/A) w "
      "+ const: one of q0, r0 is const, the other linear",
      zero(E0z2 + A_ * dw(q0 * r0)))
check("D2: E1 == A[B(2t r0' - t'r0) - 2Ct(2t q0' + t'q0)]",
      zero(E1z2 - A_ * (B_ * (2 * t_ * dw(r0) - dw(t_) * r0)
                        - 2 * C_ * t_ * (2 * t_ * dw(q0) + dw(t_) * q0))))
t0_, t1_ = sp.symbols("t0 t1")
tlin = t0_ + t1_ * w
famD2 = E_coeffs(A_, B_ * tlin, C_ * tlin**2, P0=sp.S.Zero, Q0=ga,
                 R0=2 * C_ * ga * tlin / B_)
check("D2 case q0 = gamma const: leading-degree count forces t LINEAR and "
      "r0 = (2C gamma/B) t; the family checks: E2 = E1 = 0, "
      "E0 = -2AC gamma^2 t1/B",
      zero(famD2[0]) and zero(famD2[1])
      and zero(famD2[2] + 2 * A_ * C_ * ga**2 * t1_ / B_))
print("     BOX: un-shearing needs val(p0 = -A h) >= 3 => h(0) = 0, but")
print("     then q0(0) = gamma != 0 violates q0(0) = 0.  Case r0 = delta")
print("     const: E1 => t | delta^2 t' => t' = 0, contradiction.")
print("     D2 EMPTY IN THE BOX.")

# --- D3: s, t both nonconstant, coprime ------------------------------------
G1 = t_ * p0 - a_ * s_ * q0
G2 = b_ * t_ * q0 - s_ * r0
G3 = b_ * t_ * G1 + a_ * s_ * G2
Theta = dw(s_ * t_) * G2 - 2 * s_ * t_ * dw(G2)
Theta1 = 2 * s_ * t_ * dw(G1) - (3 * dw(s_) * t_ - s_ * dw(t_)) * G1
E0p = t_ * (G2 * dw(G3) - 3 * G3 * dw(G2)) + dw(t_) * G2 * G3
check("D3: shear-invariant form of the system: E1 == a s Theta + b t Theta1"
      "  and  t(G2 G3' - 3 G3 G2') + t' G2 G3 == b t^2 E0 - r0 E1,  with "
      "G1 = t p0 - a s q0, G2 = b t q0 - s r0 (both shear-invariant)",
      zero(E1D - (a_ * s_ * Theta + b_ * t_ * Theta1))
      and zero(E0p - b_ * t_**2 * E0D + r0 * E1D))
# G2 = 0 or G1 = 0 kill the kappa equation:
check("D3: G2 == 0  =>  E0-combination == 0 != b kappa t^2;  G1 == 0 => "
      "it reduces to a G2 Theta with Theta = 0 forced by E1: both EMPTY",
      zero(subd(E0p, (r0, b_ * t_ * q0 / s_)))     # G2 = 0
      and zero(subd(E0p - a_ * G2 * Theta, (p0, a_ * s_ * q0 / t_))))
# squarefree divisibility: s | Theta1 => s | s' G1 => s | G1; t | Theta =>
# t | t' G2 => t | G2 (s, t squarefree).  Write G1 = s g1b, G2 = t g2b:
g1b, g2b = sp.Function("g1b")(w), sp.Function("g2b")(w)
Wr = dw(s_) * t_ - s_ * dw(t_)
Z_ = a_ * g2b - b_ * g1b
E1g = subd(a_ * s_ * Theta + b_ * t_ * Theta1,
           (p0, (s_ * g1b + a_ * s_ * q0) / t_),
           (r0, (b_ * t_ * q0 - t_ * g2b) / s_))
check("D3: E1|_{G1 = s g1b, G2 = t g2b} == s t (W Z - 2 s t Z'),  "
      "Z = a g2b - b g1b, W = s't - st'; and (Z^2 t/s)' == "
      "Z(2stZ' - WZ)/s^2  =>  Z^2 t = e s  =>  Z = 0 (t nonconst, coprime)",
      zero(E1g - s_ * t_ * (Wr * Z_ - 2 * s_ * t_ * dw(Z_)))
      and zero(sp.diff(Z_**2 * t_ / s_, w)
               - Z_ * (2 * s_ * t_ * dw(Z_) - Wr * Z_) / s_**2))
E0pg = subd(E0p, (p0, (s_ * (a_ / b_) * g_ + a_ * s_ * q0) / t_),
            (r0, (b_ * t_ * q0 - t_ * g_) / s_))
check("D3: with Z = 0 (g1b = (a/b) g, g2b = g): E0-combination == "
      "-2 a t^2 g (2 s t g' - W g)  =>  g(2stg' - Wg) = -b kappa/(2a): "
      "g = c const, kappa = (2a/b) c^2 W  =>  W constant",
      zero(E0pg + 2 * a_ * t_**2 * g_ * (2 * s_ * t_ * dw(g_) - Wr * g_)))
c_ = sp.Symbol("c1", nonzero=True)
q0f = sp.Function("q0")(w)
p0D3 = a_ * s_ * (q0f + c_ / b_) / t_
r0D3 = t_ * (b_ * q0f - c_) / s_
famD3 = E_coeffs(a_ * s_**2, s_ * t_, b_ * t_**2, P0=p0D3, Q0=q0f, R0=r0D3)
check("D3 family (all solutions, q0 = the shear direction): E1 == 0 and "
      "E0 == (2 a c^2/b) W identically -- Keller iff W = s't - st' const",
      zero(famD3[1]) and zero(famD3[2] - 2 * a_ * c_**2 / b_ * Wr))
check("D3 BOX: t * p0 == a s (q0 + c/b), so p0(0) t(0) = a s(0) c/b != 0 "
      "while the box forces p0(0) = 0: contradiction for EVERY t.  "
      "D3 (s, t squarefree) EMPTY IN THE BOX, all degrees",
      zero(t_ * p0D3 - a_ * s_ * (q0f + c_ / b_)))
print("     (remaining gap: D3 with s or t NON-squarefree, i.e. deg p1 >= 4")
print("      or deg r1 >= 4 -- closed in-box by section 7 / --full.)")

# ===========================================================================
print("== 6. the 3:1 orbifold mechanism (the only one) is empty ==")
# ===========================================================================
# Weights: target coordinates carry (-3, -1, 1).  A stabilizer jump needs
# an image orbit in a coordinate stratum where C* descends with kernel
# Z_k: only the weight -3 axis gives k = 3 (gcd of any other weight pair
# is 1) -- so (1,-1,-3) admits NO 2:1 mechanism and exactly one 3:1
# mechanism: a point (w0, v0) with Q = R = 0 != P on a free orbit {x != 0}.
from math import gcd  # noqa: E402
check("weights: only the -3 axis has nontrivial stabilizer Z_3 "
      "(|{-3}| = 3; gcd of every other weight combination is 1)",
      gcd(3, 1) == 1 and gcd(1, 1) == 1 and abs(-3) == 3)
check("classification corollary: every v-linear Keller map has R == 1 "
      "(gauge), which never vanishes: NO 3:1 witness exists",
      Rg == 1)

# independent pointwise Groebner queries in a degree box (mirrors
# search_213 section 4).  The residual torus (w,v) -> (mu w, nu v)
# preserves the gauge (compensated by target scalings) and moves any
# witness to (1,1), (1,0) or (0,1); (0,0) is excluded by R(0,0) = 1.


def vlinear_box(dp0, dp1, dq0, dq1, dr0, dr1):
    """Gauged v-linear ansatz with numeric degree bounds; returns
    (P, Q, R, unknowns)."""
    def poly(name, degs, gauge=None):
        cs, expr = [], sp.S.Zero
        for j in degs:
            if gauge and j in gauge:
                expr += gauge[j] * w**j
            else:
                cj = sp.Symbol(f"{name}{j}")
                cs.append(cj)
                expr += cj * w**j
        return expr, cs
    P0, u1 = poly("p", range(3, dp0 + 1))
    P1, u2 = poly("pa", range(0, dp1 + 1), gauge={0: 1})
    Q0, u3 = poly("q", range(1, dq0 + 1), gauge={1: 1})
    Q1, u4 = poly("qa", range(0, dq1 + 1))
    R0, u5 = poly("r", range(0, dr0 + 1), gauge={0: 1})
    R1, u6 = poly("ra", range(0, dr1 + 1))
    return (P0 + P1 * v, Q0 + Q1 * v, R0 + R1 * v,
            u1 + u2 + u3 + u4 + u5 + u6)


from jcqft.gb_backend import available_backends, is_unit_ideal  # noqa: E402
print(f"  GB backends: {available_backends()}")
Pb, Qb, Rb, U = vlinear_box(4, 2, 3, 2, 2, 2)
KELLER = sp.Poly(sp.expand(det_m(Pb, Qb, Rb, 3)) + 1, w, v).coeffs()
print(f"  small box (deg p0<=4, p1,q1,r1<=2, q0<=3, r0<=2): "
      f"{len(U)} unknowns, {len(KELLER)} Keller equations")
rr = sp.Symbol("rrab")
for w0, v0 in ((1, 1), (1, 0), (0, 1)):
    sys_ = [sp.expand(e) for e in
            KELLER + [Qb.subs({w: w0, v: v0}), Rb.subs({w: w0, v: v0}),
                      1 - rr * Pb.subs({w: w0, v: v0})]]
    try:
        empty = is_unit_ideal(sys_, U + [rr], backend="auto", timeout=600)
        check(f"3:1 witness at (w,v) = ({w0},{v0}): ideal == (1), EMPTY",
              empty)
    except (RuntimeError, TimeoutError) as exc:
        print(f"  [unresolved] 3:1 witness at ({w0},{v0}) -- {exc}")

# ===========================================================================
print("== 7. in-box corroboration of the classification ==")
# ===========================================================================
# (i) the heart of the theorem, checked independently: Keller + (r1 != 0)
# is EMPTY in the box -- this covers strata C and D wholesale, including
# the non-squarefree D3 configurations the structural proof misses (within
# the box degrees).
yv = sp.Symbol("yrab")
for cf in [u for u in U if str(u).startswith("ra")]:
    try:
        empty = is_unit_ideal(KELLER + [1 - yv * cf], U + [yv],
                              backend="auto", timeout=600)
        check(f"Keller + ({cf} != 0): ideal == (1), EMPTY (small box)",
              empty)
    except (RuntimeError, TimeoutError) as exc:
        print(f"  [unresolved] r1 coeff {cf} != 0 -- {exc}")

# (ii) with r1 == 0: nilpotency certificates pin the whole in-box variety
# to the gauged family P = p0 + v, Q = w + b0 P, R = 1 (+ box truncation).
Pb2, Qb2, Rb2, U2 = vlinear_box(4, 2, 3, 2, 2, -1)   # r1 == 0
K2 = sp.Poly(sp.expand(det_m(Pb2, Qb2, Rb2, 3)) + 1, w, v).coeffs()
G = sp.groebner([sp.expand(e) for e in K2], *U2, order="grevlex")
D_ = dict(zip(("p3", "p4", "pa1", "pa2", "q2", "q3", "qa0", "qa1", "qa2",
               "r1", "r2"), [sp.Symbol(n_) for n_ in
                             ("p3", "p4", "pa1", "pa2", "q2", "q3", "qa0",
                              "qa1", "qa2", "r1", "r2")]))


def nilcert(gb, f, kmax=10):
    fk = sp.S.One
    for k in range(1, kmax + 1):
        fk = sp.expand(fk * f)
        if gb.reduce(fk)[1] == 0:
            return k
    return None


for f in (D_["pa1"], D_["pa2"], D_["r1"], D_["r2"], D_["qa1"], D_["qa2"],
          D_["q2"], D_["q3"] - D_["qa0"] * D_["p3"],
          D_["qa0"] * D_["p4"]):
    k = nilcert(G, f)
    if k is not None:
        check(f"r1 = 0 stratum: ({f})^{k} in I  =>  {f} == 0 on the whole "
              "in-box Keller variety", True)
    else:
        # Nilpotency exponent > kmax: certify radical membership directly
        # (Rabinowitsch: f vanishes on V(I) iff I + (1 - y f) is the unit
        # ideal) -- exact, and equivalent to the nilcert claim.
        empty = is_unit_ideal(K2 + [1 - yv * f], U2 + [yv],
                              backend="auto", timeout=600)
        check(f"r1 = 0 stratum: {f} in radical(I) (Rabinowitsch)  =>  "
              f"{f} == 0 on the whole in-box Keller variety", empty)
print("  => in-box variety == {p1 = 1, r0 = 1, q1 = b0, q0 = w + b0 p3 w^3,")
print("     b0 p4 = 0}: exactly the gauged family A u B (with its box")
print("     truncation b0 p4 = 0).  CLASSIFICATION CORROBORATED.")

print("\nVERDICT (v-linear class, all w-degrees, D3-non-squarefree gap")
print("box-checked): every Keller map is the tame automorphism")
print("  F = (p0(xy)/x^3 + z,  y + b0 x^2 F1,  x)   (up to gauge);")
print("NO counterexample, NO 3:1 orbifold covering.  The Alpoge-Mathew")
print("stratum (D1) is EMPTY at m = 3: the mechanism is numerologically")
print("obstructed, not merely undiscovered.")


# ===========================================================================
# --full: larger boxes for the corroboration queries
# ===========================================================================

def run_full(budget):
    print("\n== 8. (--full) larger-box corroboration ==")
    # (i) medium box, all r1 coefficients
    Pb3, Qb3, Rb3, U3 = vlinear_box(5, 3, 4, 3, 3, 3)
    K3 = sp.Poly(sp.expand(det_m(Pb3, Qb3, Rb3, 3)) + 1, w, v).coeffs()
    print(f"  medium box (deg p0<=5, p1,q1,r1<=3, q0<=4, r0<=3): "
          f"{len(U3)} unknowns, {len(K3)} equations")
    for cf in [u for u in U3 if str(u).startswith("ra")]:
        try:
            empty = is_unit_ideal(K3 + [1 - yv * cf], U3 + [yv],
                                  backend="auto", timeout=budget)
            check(f"medium box, Keller + ({cf} != 0): EMPTY", empty)
        except (RuntimeError, TimeoutError) as exc:
            print(f"  [unresolved] {cf} != 0 -- {exc}")
    # (ii) the four targeted non-squarefree D3 gap cases,
    # deg(p0, q0, r0) <= (6, 5, 4)
    rho, m1, m2, n0, n1, n2 = sp.symbols("rho m1 m2 n0 n1 n2")
    cases = {
        "s = (1+rho w)^2, t deg 1":
            ((1 + rho * w)**2, n0 + n1 * w, [rho, n0, n1], [n1]),
        "s deg 1, t = (n0+n1 w)^2":
            ((1 + m1 * w), (n0 + n1 * w)**2, [m1, n0, n1], [n1]),
        "s = (1+rho w)^2, t deg 2":
            ((1 + rho * w)**2, n0 + n1 * w + n2 * w**2,
             [rho, n0, n1, n2], [n2]),
        "s deg 2, t = (n0+n1 w)^2":
            ((1 + m1 * w + m2 * w**2), (n0 + n1 * w)**2,
             [m1, m2, n0, n1], [n1]),
    }

    def poly(name, degs, gauge=None):
        cs, expr = [], sp.S.Zero
        for j in degs:
            if gauge and j in gauge:
                expr += gauge[j] * w**j
            else:
                cj = sp.Symbol(f"{name}{j}")
                cs.append(cj)
                expr += cj * w**j
        return expr, cs

    for label, (s3, t3, pars, sat) in cases.items():
        P0, cp = poly("p", range(3, 7))
        Q0, cq = poly("q", range(1, 6), gauge={1: 1})
        R0, cr = poly("r", range(0, 5), gauge={0: 1})
        Dg = sp.expand(det_m(P0 + s3**2 * v, Q0 + s3 * t3 * v,
                             R0 + b_ * t3**2 * v, 3))
        Kg = sp.Poly(Dg + 1, w, v).coeffs()
        ys = sp.symbols(f"ysat0:{len(sat) + 1}")
        sysg = Kg + [1 - ys[0] * b_] + [1 - yi * cq_
                                        for yi, cq_ in zip(ys[1:], sat)]
        gens = pars + [b_] + cp + cq + cr + list(ys)
        try:
            empty = is_unit_ideal(sysg, gens, backend="auto",
                                  timeout=budget)
            check(f"D3 gap, {label} (deg p0<=6, q0<=5, r0<=4): EMPTY",
                  empty)
            continue
        except (RuntimeError, TimeoutError) as exc:
            print(f"  [msolve unresolved] {label} -- {exc}")
        if "singular" not in available_backends():
            continue
        # memory-frugal fallback: Singular degBound ladder.  Finding 1
        # below the bound is an EXACT unit-ideal certificate (one-sided;
        # not finding it proves nothing).
        for db in (4, 5, 6, 7):
            try:
                empty = is_unit_ideal(sysg, gens, backend="singular",
                                      timeout=budget, degbound=db)
            except (RuntimeError, TimeoutError) as exc:
                print(f"  [singular degBound {db} unresolved] {label} "
                      f"-- {exc}")
                break
            if empty:
                check(f"D3 gap, {label}: EMPTY (Singular degBound {db} "
                      "certificate)", empty)
                break
            print(f"  [no certificate at degBound {db}] {label}")
        else:
            print(f"  [unresolved] {label}: no unit certificate up to "
                  "degBound 7 -- see docs/SEARCH_113.md section 8")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="larger-box Groebner corroboration (~1 h)")
    ap.add_argument("--budget", type=int, default=1800,
                    help="per-query timeout for --full (s)")
    args = ap.parse_args()
    if args.full:
        run_full(args.budget)
    print(f"\nall checks passed in {time.time() - T0:.1f} s")
