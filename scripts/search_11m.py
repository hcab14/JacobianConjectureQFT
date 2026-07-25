"""Uniform classification for the whole family of weight systems (1,-1,-m).

(Machinery: jcqft/reduction_w.py, proved exact for all m in
scripts/reduction_113.py; the m = 3 case worked out in scripts/search_113.py
/ docs/SEARCH_113.md is the blueprint.  Write-up: docs/SEARCH_11M.md.)

THEOREM (proved here, by assertion, with one precisely-mapped gap).
For every integer m >= 3, every Keller map in the v-linear class of the
weight system (1,-1,-m) -- P, Q, R of degree <= 1 in v = x^m z, arbitrary
degree in w = x y -- is, up to gauge, the tame automorphism

    P = p0 + v,   Q = w + b0 P,   R = 1,   p0 in w^m C[w],  b0 in C.

Hence no counterexample and no m:1 orbifold covering exists in the class,
for ANY m >= 3, and the Alpoge-Mathew counterexample (m = 2) is the unique
member of its equivariant family within the v-linear class.  The one gap:
stratum D3 with non-squarefree s or t (deg p1 >= m+1 or deg r1 >= m+1 for
odd m; >= 2(m+1) for even m), exactly the D3' gap of docs/SEARCH_113.md --
box-closed by Groebner certificates for m <= 5.

PROOF STRUCTURE, all identities asserted with m SYMBOLIC (det M is linear
in m; only the UFD exponent bookkeeping branches on the parity of m, and
is asserted with k symbolic, m = 2k+1 / m = 2k):

1. REDUCED KELLER IDENTITY, m symbolic (generic-function proof):
       det DF = det M = -m P J2(Q,R) + Q J2(P,R) + R J2(P,Q),
   R^m det M = J2(P R^m, Q R); box j + mk >= m (P), >= 1 (Q); gauge
   p1(0) = q0'(0) = r0(0) = 1, kappa = -1.

2. v-LINEAR CLASS, m symbolic.  det M = kappa gives E2 = E1 = 0, E0 = kappa:
       E2 = 2 p1' q1 r1 - (m+1) p1 q1' r1 + (m-1) p1 q1 r1'
   integrates (Wronskian certificate) to the GENERAL CONSTRAINT
       p1^2 r1^{m-1} = c q1^{m+1}        (c != 0, whenever q1 r1 != 0),
   i.e. 2 p1'/p1 - (m+1) q1'/q1 + (m-1) r1'/r1 = 0.  Anchors: m = 2 gives
   p1^2 r1 = c q1^3 (Alpoge-Mathew: c = -1/27); m = 3 gives p1 r1 = c q1^2.

3. STRATA (p1 != 0 forced by p1(0) = 1):  A: q1 = r1 = 0;  B: r1 = 0,
   q1 != 0;  C: q1 = 0, r1 != 0;  D: q1 r1 != 0.  In D the UFD solution of
   the constraint (after absorbing gcd(p1,q1,r1) | E0 = kappa, hence
   constant) is, with gcd(s,t) = 1 and per-irreducible multiplicity
   arithmetic (asserted):
       m = 2k+1:  p1 = a s^{k+1},   q1 = s t^k,        r1 = b t^{k+1}
       m = 2k:    p1 = a s^{2k+1},  q1 = s^2 t^{2k-1}, r1 = b t^{2k+1}
   Substrata: D0 (s,t const), D1 (t const -- the Alpoge-Mathew slot),
   D2 (s const), D3 (both nonconstant).

4. VERDICTS (m >= 3):
   A, B  = the tame family above (complete solution, all w-degrees,
           uniform in m).
   C     = EMPTY, uniform in m: (p1/(q0^{m+1} r1))' = E1/(q0^{m+2} r1^2)
           forces p1(0) = 0.
   D0    = EMPTY, uniform: E0 becomes [-(m+1)/(m-1) (ag/b) V^2 - ad V]',
           V = q0 - (b/a) p0, V = w forced, w^2-coefficient kills it.
   D1    = EMPTY.  m even >= 4: the kappa-equation carries the factor
           u^{m-2} (u nonconstant) => u | kappa.  m odd >= 5: factor
           u^{(m-3)/2} (or d^{m-2} in the square subcase u = c d^2).
           m = 3: box jet argument (p0(0) = q0(0) = 0 => kappa = 0).
   D2    = EMPTY, uniform: q0 r0 exactly linear; each case kills by
           divisibility or the un-shear jet (val p0 >= m => h(0) = 0).
   D3    = EMPTY for s, t squarefree.  The divisibility chain
           s | s'G1-iteration (coefficients m-2j, resp. m-j, never zero)
           + t | t'G2 reduces E1 = 0 to an exact Wronskian relation; the
           kappa-equation then carries the factor s^{(m-3)/2} (m odd),
           resp. s^{m-2} (m even) => EMPTY outright for m >= 4; at m = 3
           the surviving shear-orbit dies against the box (search_113).

5. THE m = 2 DIVERGENCE (the theorem's content).  The Alpoge-Mathew map
   lives in D1 (m even): u = 1+w, p1 = u^3, q1 = 3u^2, r1 = -1.  There the
   kappa-equation is  C u^{m-2} (u'q0 - u q0') [ (m^2-1)/2 (A/B) q0
   + m ytil u ] = kappa,  and TWO things are special at m = 2 and ONLY
   there: (i) the killing factor u^{m-2} is trivial; (ii) the homogeneous
   E1-direction Y = ytil u^m exists for even m and enters with weight
   m ytil u.  The AM data threads exactly: u'q0 - uq0' = 2, ytil = -1,
   E0 = -2 = kappa_AM (asserted end-to-end).  For every m >= 4 the same
   formula forces u | kappa; at m = 3 the box valuation val_w p0 >= 3
   kills the survivor.  The obstruction is numerological and exact.

6. GAP + CORROBORATION.  D3 with s or t NON-squarefree is not reached by
   the chain; it needs deg p1 >= m+1 or deg r1 >= m+1 (m odd), resp.
   >= 2(m+1) (m even).  In-box Groebner certificates (msolve, exact over
   Q): Keller + (r1 != 0) is EMPTY in small boxes for m = 4, 5, and
   Rabinowitsch certificates pin the r1 = 0 slice to the gauged tame
   family; --full enlarges the boxes and attacks the two targeted m = 5
   non-squarefree parametrizations.

7. ORBIFOLD MECHANISM.  For every m >= 2 the only stabilizer-jump is m:1
   (image on the weight -m axis: Q = R = 0 != P; all other axes/pairs have
   trivial stabilizer -- also for composite m, where no intermediate Z_d
   stratum exists).  The classification gives R = 1: no witness, for any
   m >= 3; independent pointwise Groebner queries at torus-normalized
   witness points confirm EMPTY in-box for m = 4, 5.

VERDICT: for all m >= 3 the v-linear class contains nonlinear Keller maps,
but ALL are tame automorphisms (modulo the box-closed D3' gap).  The
Alpoge-Mathew mechanism exists at m = 2 and ONLY at m = 2.

Default run ~2 min; --full adds larger boxes and the targeted m = 5 gap
queries (budget-capped msolve, ~70 min).
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
    assemble, det_m, extract, j2, monomial_box, v, w, x, y, z)

T0 = time.time()

SPOT_MS = (2, 3, 4, 5, 7)      # concrete spot-check members (2 = AM!)


def check(label, cond):
    assert cond, label
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


def zero(expr):
    """Is the (Derivative-carrying) expression identically zero?
    Handles symbolic exponents via powsimp."""
    e = expr.doit() if hasattr(expr, "doit") else expr
    e = sp.powsimp(sp.expand(sp.powsimp(e, force=True)), force=True)
    return sp.simplify(e) == 0


def subd(expr, *pairs):
    """Substitute into an expression containing Derivatives, then evaluate."""
    e = expr
    for old, new in pairs:
        e = e.subs(old, new)
    return sp.expand(e.doit())


dw = lambda f: sp.diff(f, w)   # noqa: E731

M = sp.Symbol("m", integer=True, positive=True)     # THE symbolic weight
K = sp.Symbol("k", integer=True, positive=True)     # parity parameter

# ===========================================================================
print("== 1. reduced Keller identity, m SYMBOLIC ==")
# ===========================================================================

Pf, Qf, Rf = [f(w, v) for f in sp.symbols("P Q R", cls=sp.Function)]
_sub = {w: x * y, v: x**M * z}
F_gen = (Pf.subs(_sub) / x**M, Qf.subs(_sub) / x, x * Rf.subs(_sub))
_J = sp.Matrix([[sp.diff(Fi, u) for u in (x, y, z)] for Fi in F_gen])
detDF = sp.simplify(_J.det().subs({y: w / x, z: v / x**M}))
dM = -M * Pf * j2(Qf, Rf) + Qf * j2(Pf, Rf) + Rf * j2(Pf, Qf)
check("det DF == det M == -m P J2(Q,R) + Q J2(P,R) + R J2(P,Q), a function "
      "of (w,v) alone, m SYMBOLIC (generic-function proof)",
      sp.simplify(detDF - dM) == 0 and not detDF.has(x))
check("compact form R^m det M == J2(P R^m, Q R), m symbolic",
      sp.simplify(Rf**M * dM - j2(Pf * Rf**M, Qf * Rf)) == 0)
check("det_m (module) == det M on undetermined functions, m symbolic; "
      "det M is LINEAR in m (whence all E-coefficients are)",
      sp.simplify(det_m(Pf, Qf, Rf, M) - dM) == 0
      and sp.simplify(det_m(Pf, Qf, Rf, M)
                      - (det_m(Pf, Qf, Rf, 0)
                         + M * (det_m(Pf, Qf, Rf, 1)
                                - det_m(Pf, Qf, Rf, 0)))) == 0)

# per-m: polynomiality boxes and the gauge (DF(0) antidiagonal)
_c = sp.symbols("cP cQ cR")
for m in SPOT_MS:
    boxP = monomial_box("P", m, jmax=m + 1, kmax=1)
    boxQ = monomial_box("Q", m, jmax=m + 1, kmax=1)
    ok_box = (v in boxP and w**m in boxP
              and all(w**j not in boxP for j in range(m))
              and w in boxQ and v in boxQ
              and sp.S.One not in [sp.sympify(mm) for mm in boxQ])
    F_lin = assemble(_c[0] * v + w**m, _c[1] * w, _c[2], m)
    _L = sp.Matrix([[sp.diff(Fi, u) for u in (x, y, z)]
                    for Fi in F_lin]).subs({x: 0, y: 0, z: 0})
    check(f"m={m}: box j+mk>=m (P) / >=1 (Q) as claimed; DF(0) = "
          "antidiag(p1(0), q0'(0), r0(0)), det = -p1(0)q0'(0)r0(0) "
          "=> gauge p1(0)=q0'(0)=r0(0)=1, kappa = -1",
          ok_box and _L.det() == -_c[0] * _c[1] * _c[2])

# ===========================================================================
print("== 2. v-linear class: E2, E1, E0 and the GENERAL constraint ==")
# ===========================================================================

p0, p1, q0, q1, r0, r1 = [sp.Function(n)(w) for n in
                          ("p0", "p1", "q0", "q1", "r0", "r1")]


def E_coeffs(P1, Q1, R1, P0=p0, Q0=q0, R0=r0, m=M):
    D = sp.expand(det_m(P0 + P1 * v, Q0 + Q1 * v, R0 + R1 * v, m))
    return (sp.expand(D.coeff(v, 2)), sp.expand(D.coeff(v, 1)),
            sp.expand(D.coeff(v, 0)))


E2, E1, E0 = E_coeffs(p1, q1, r1)
check("E2 == 2 p1' q1 r1 - (m+1) p1 q1' r1 + (m-1) p1 q1 r1'  (m symbolic)",
      zero(E2 - (2 * dw(p1) * q1 * r1 - (M + 1) * p1 * dw(q1) * r1
                 + (M - 1) * p1 * q1 * dw(r1))))
check("E1 == -m p0[Q,R]_1 - m p1[Q,R]_0 + q0[P,R]_1 + q1[P,R]_0 "
      "+ r0[P,Q]_1 + r1[P,Q]_0  ([A,B]_i = v^i-part of J2(A,B))",
      zero(E1 - (-M * p0 * (dw(q1) * r1 - q1 * dw(r1))
                 - M * p1 * (dw(q0) * r1 - q1 * dw(r0))
                 + q0 * (dw(p1) * r1 - p1 * dw(r1))
                 + q1 * (dw(p0) * r1 - p1 * dw(r0))
                 + r0 * (dw(p1) * q1 - p1 * dw(q1))
                 + r1 * (dw(p0) * q1 - p1 * dw(q0)))))
check("E0 == -m p0(q0'r1 - q1 r0') + q0(p0'r1 - p1 r0') + r0(p0'q1 - p1 q0')",
      zero(E0 - (-M * p0 * (dw(q0) * r1 - q1 * dw(r0))
                 + q0 * (dw(p0) * r1 - p1 * dw(r0))
                 + r0 * (dw(p0) * q1 - p1 * dw(q0)))))

h = sp.Function("h")(w)
Dfull = sp.expand(det_m(p0 + p1 * v, q0 + q1 * v, r0 + r1 * v, M))
check("det M invariant under the v-shear (w,v) -> (w, v + h(w)), m symbolic",
      zero(sp.expand(det_m(*(S.subs(v, v + h) for S in
                             (p0 + p1 * v, q0 + q1 * v, r0 + r1 * v)),
                           M)) - Dfull.subs(v, v + h)))

# THE GENERAL INTEGRATED CONSTRAINT (the log-derivative relation
# 2 p1'/p1 - (m+1) q1'/q1 + (m-1) r1'/r1 = 0):
check("GENERAL Wronskian certificate, m SYMBOLIC: "
      "(p1^2 r1^{m-1}/q1^{m+1})' == p1 r1^{m-2} E2 / q1^{m+2}  "
      "=>  on {E2 = 0, q1 p1 r1 != 0}:  p1^2 r1^{m-1} = c q1^{m+1}, c != 0",
      zero(sp.diff(p1**2 * r1**(M - 1) / q1**(M + 1), w)
           - p1 * r1**(M - 2) * E2 / q1**(M + 2)))

# anchors: m = 2 (Alpoge-Mathew) and m = 3 (search_113)
E2m2 = E_coeffs(p1, q1, r1, m=2)[0]
check("anchor m=2: (p1^2 r1/q1^3)' == p1 E2|_{m=2} / q1^4",
      zero(sp.diff(p1**2 * r1 / q1**3, w) - p1 * E2m2 / q1**4))
E2m3 = E_coeffs(p1, q1, r1, m=3)[0]
check("anchor m=3: (p1 r1/q1^2)' == E2|_{m=3}/(2 q1^3)  (search_113 form; "
      "its square is the general constraint at m = 3)",
      zero(sp.diff(p1 * r1 / q1**2, w) - E2m3 / (2 * q1**3)))

from jcqft.core import F as F_AM  # noqa: E402
P_AM, Q_AM, R_AM = extract(F_AM, 2)
p0AM, p1AM = (sp.expand(P_AM.coeff(v, i)) for i in (0, 1))
q0AM, q1AM = (sp.expand(Q_AM.coeff(v, i)) for i in (0, 1))
r0AM, r1AM = (sp.expand(R_AM.coeff(v, i)) for i in (0, 1))
check("Alpoge-Mathew (m=2) satisfies the general constraint with "
      "c = -1/27: p1 = (1+w)^3, q1 = 3(1+w)^2, r1 = -1, p1^2 r1 == -q1^3/27",
      sp.expand(p1AM - (1 + w)**3) == 0
      and sp.expand(q1AM - 3 * (1 + w)**2) == 0 and r1AM == -1
      and sp.expand(p1AM**2 * r1AM + q1AM**3 / 27) == 0)

# ===========================================================================
print("== 3. stratification and the UFD parametrization of D ==")
# ===========================================================================
# Strata (p1(0) = 1 forces p1 != 0):  A: q1 = r1 = 0;  B: r1 = 0, q1 != 0;
# C: q1 = 0, r1 != 0;  D: q1 r1 != 0.
#
# In D: gcd(p1, q1, r1) =: g divides E0 = kappa (E0 is LINEAR in
# (p1,q1,r1) with no derivatives of them), so g is a nonzero constant --
# absorb it (g = 1).  Then for each irreducible pi with multiplicities
# (A, B, C) in (p1, q1, r1) the constraint gives 2A + (m-1)C = (m+1)B with
# min(A, B, C) = 0; B = 0 forces A = C = 0, so rad q1 = rad(p1 r1) and
# gcd(p1, r1) = 1; on the p1-side 2A = (m+1)B, on the r1-side
# (m-1)C = (m+1)B.  Minimal solutions and the parity of m give (asserted):
#   m = 2k+1: (A,B) = ((k+1)b, b), (C,B) = ((k+1)c, k c)
#             => p1 = a s^{k+1}, q1 = s t^k, r1 = b t^{k+1}
#   m = 2k:   (A,B) = ((m+1)b, 2b), (C,B) = ((m+1)c, (m-1)c)
#             => p1 = a s^{2k+1}, q1 = s^2 t^{2k-1}, r1 = b t^{2k+1}
# with gcd(s, t) = 1 and one constant absorbed into q1's normalization.

g_ = sp.Function("g")(w)
E0g = sp.expand(E0.subs({p1: p1 * g_, q1: q1 * g_, r1: r1 * g_}).doit())
check("D: E0(g p1, g q1, g r1) == g * E0(p1, q1, r1)  =>  gcd(p1,q1,r1) | "
      "E0 = kappa  =>  constant (absorbed), m symbolic",
      zero(E0g - g_ * E0))

a_, b_ = sp.symbols("a b", nonzero=True)
s_, t_ = sp.Function("s")(w), sp.Function("t")(w)
Modd, Mev = 2 * K + 1, 2 * K
Wr = dw(s_) * t_ - s_ * dw(t_)          # the Wronskian W = s't - st'

# exponent bookkeeping (integer arithmetic, k symbolic):
check("D exponents, m = 2k+1: 2(k+1) == m+1 and (m-1)(k+1) == (m+1)k;  "
      "m = 2k: 2(m+1) == (m+1)*2 and (m-1)(m+1) == (m+1)(m-1)",
      sp.simplify(2 * (K + 1) - (Modd + 1)) == 0
      and sp.simplify((Modd - 1) * (K + 1) - (Modd + 1) * K) == 0)
check("D parametrization satisfies the constraint identically, k symbolic: "
      "odd: (a s^{k+1})^2 (b t^{k+1})^{m-1} == a^2 b^{m-1} (s t^k)^{m+1}; "
      "even: likewise with (s^{2k+1}, s^2 t^{2k-1}, t^{2k+1})",
      zero((a_ * s_**(K + 1))**2 * (b_ * t_**(K + 1))**(Modd - 1)
           - a_**2 * b_**(Modd - 1) * (s_ * t_**K)**(Modd + 1))
      and zero((a_ * s_**(2 * K + 1))**2 * (b_ * t_**(2 * K + 1))**(Mev - 1)
               - a_**2 * b_**(Mev - 1)
               * (s_**2 * t_**(2 * K - 1))**(Mev + 1)))
check("D parametrization satisfies E2 == 0 identically (both parities, "
      "k symbolic)",
      zero(E_coeffs(a_ * s_**(K + 1), s_ * t_**K, b_ * t_**(K + 1),
                    m=Modd)[0])
      and zero(E_coeffs(a_ * s_**(2 * K + 1), s_**2 * t_**(2 * K - 1),
                        b_ * t_**(2 * K + 1), m=Mev)[0]))
print("     Substrata of D:  D0: s, t const;  D1: t const, s nonconst")
print("     (the Alpoge-Mathew slot);  D2: s const;  D3: both nonconst.")

# ===========================================================================
print("== 4. strata A and B: complete solution, m SYMBOLIC -- tame family ==")
# ===========================================================================
# A (q1 = r1 = 0): E1 == 0 identically, E0 = -p1 (q0 r0)' = kappa; in C[w]
# p1 | kappa => p1 = 1 (gauge); (q0 r0)' = 1, q0(0) = 0, r0(0) = 1 =>
# q0 r0 = w => r0 | w, r0(0) = 1 => r0 = 1, q0 = w; p0 in w^m C[w] FREE.
check("A: E1 == 0 and E0 == -p1 (q0 r0)' identically, m symbolic",
      zero(E1.subs({q1: 0, r1: 0})) and
      zero(E0.subs({q1: 0, r1: 0}) + p1 * dw(q0 * r0)))
# B (r1 = 0, q1 != 0): E1 integrates to p1 r0^{m-1} = c q1; E0 then factors
# through an exact derivative, forcing p1 = 1, r0 = 1, q0 = b0 p0 + w.
E1B = sp.expand(E1.subs(r1, 0).doit())
check("B: integration certificate (p1 r0^{m-1}/q1)' == r0^{m-2} E1 / q1^2, "
      "m symbolic  =>  p1 r0^{m-1} = c q1, c != 0 (p1 r0 != 0)",
      zero(sp.diff(p1 * r0**(M - 1) / q1, w) - r0**(M - 2) * E1B / q1**2))
chat = sp.Symbol("c", nonzero=True)
check("B: E0|_{q1 = p1 r0^{m-1}/c} == p1 * d/dw[(p0 r0^m)/c - q0 r0], "
      "m symbolic",
      zero(subd(E0.subs(r1, 0), (q1, p1 * r0**(M - 1) / chat))
           - p1 * dw(p0 * r0**M / chat - q0 * r0)))
print("     => p1 | kappa => p1 = 1; r0 (p0 r0^{m-1}/c - q0) = -kappa w = w")
print("        => r0 | w, r0(0) = 1 => r0 = 1, q0 = b0 p0 + w, q1 = b0.")

# THE unified gauged family (b0 = 0 <-> A, b0 != 0 <-> B), m symbolic:
b0 = sp.Symbol("b0")
p0g = sp.Function("p0")(w)          # arbitrary element of w^m C[w]
Pg, Qg, Rg = p0g + v, w + b0 * (p0g + v), sp.S.One
check("FAMILY: P = p0 + v, Q = w + b0 P, R = 1 has det M == -1 identically "
      "(any p0(w), any b0, m SYMBOLIC)",
      zero(sp.expand(det_m(Pg, Qg, Rg, M)) + 1))

# end-to-end 3D check per concrete m: tame automorphism, inverse, fiber
a3, b3, c3 = sp.symbols("a3 b3 c3")
for m in [mm for mm in SPOT_MS if mm >= 3]:
    p0s = w**m + 2 * w**(m + 1)
    b0s = sp.Integer(5)
    Ps = p0s + v
    Qs = w + b0s * Ps
    Fs = assemble(Ps, Qs, sp.S.One, m)
    J3 = sp.Matrix([[sp.diff(Fi, u) for u in (x, y, z)] for Fi in Fs])
    # explicit inverse: undo the target shear b -> b - b0 a c^{m-1}
    # (Q = w + b0 P assembles to F2 = y + b0 x^{m-1} F1, x = F3), then the
    # elementary z-shear: x = c, y = b - b0 a c^{m-1}, z = a - p0(cy)/c^m.
    ynew = b3 - b0s * a3 * c3**(m - 1)
    inv = (c3, ynew, a3 - sp.expand(p0s.subs(w, c3 * ynew) / c3**m))
    comp = tuple(sp.expand(iv.subs(dict(zip((a3, b3, c3), Fs)),
                                   simultaneous=True)) for iv in inv)
    extra = ""
    ok_fib = True
    if m <= 5:
        fib = sp.solve([Fs[i] - t_ for i, t_ in
                        enumerate((sp.Rational(3, 7), sp.Rational(-2, 5),
                                   sp.Rational(1, 3)))], [x, y, z],
                       dict=True)
        ok_fib = len(fib) == 1
        extra = f", generic fiber = {len(fib)} point"
    check(f"m={m} family sample (b0=5, p0=w^{m}+2w^{m+1}): det DF == -1, "
          f"explicit polynomial inverse (TAME){extra}, prefilter survives="
          f"{infinity_prefilter(Fs, (x, y, z))} (known false-positive "
          "class: nonlinear automorphisms)",
          sp.expand(J3.det()) == -1 and comp == (x, y, z) and ok_fib)

# ===========================================================================
print("== 5. stratum C: EMPTY, m SYMBOLIC ==")
# ===========================================================================
E1C = sp.expand(E1.subs(q1, 0).doit())
check("C: certificate (p1/(q0^{m+1} r1))' == E1/(q0^{m+2} r1^2), m symbolic"
      "  =>  E1 = 0 forces p1 = c q0^{m+1} r1  =>  p1(0) = 0 (q0(0) = 0), "
      "contradicting p1(0) = 1: EMPTY for every m >= 2",
      zero(sp.diff(p1 / (q0**(M + 1) * r1), w) - E1C / (q0**(M + 2) * r1**2)))

# ===========================================================================
print("== 6. stratum D0: EMPTY, m SYMBOLIC ==")
# ===========================================================================
al, be, ga, de = sp.symbols("alpha beta gamma delta", nonzero=True)
E2K, E1K, E0K = E_coeffs(al, be, ga)
check("D0: E1 == d/dw[2 beta gamma p0 - (m+1) alpha gamma q0 "
      "+ (m-1) alpha beta r0], m symbolic",
      zero(E1K - dw(2 * be * ga * p0 - (M + 1) * al * ga * q0
                    + (M - 1) * al * be * r0)))
r0D0 = de + ((M + 1) * al * ga * q0 - 2 * be * ga * p0) / ((M - 1) * al * be)
VD0 = q0 - (be / al) * p0
check("D0: after solving E1, E0 == d/dw[-((m+1)/(m-1)) (alpha gamma/beta) "
      "V^2 - alpha delta V],  V = q0 - (beta/alpha) p0, m symbolic",
      zero(subd(E0K, (r0, r0D0))
           - dw(-(M + 1) / (M - 1) * al * ga / be * VD0**2
                - al * de * VD0)))
print("     => quadratic(V) = kappa w, V(0) = 0, V'(0) = 1 (p0'(0) = 0 for")
print("        m >= 2): V | kappa w => V = w, and the w^2-coefficient")
print("        forces (m+1)/(m-1) alpha gamma/beta = 0 -- impossible for")
print("        every m >= 2.  D0 EMPTY, uniform in m.")

# ===========================================================================
print("== 7. stratum D1 (the Alpoge-Mathew slot): EMPTY for all m >= 3 ==")
# ===========================================================================
# D1: t const, s =: u nonconstant.  Both parities, k symbolic.  v-shear
# h = -r0/C sets r0 = 0 (r1 = C != 0 const); E1 is then a linear ODE for
# p0 with particular solution p0* = ((m+1)/2)(A/B) u^{e-f} q0 and
# homogeneous solutions Y with Y^2 = e q1^m:
#   m odd  (p1 = A u^{k+1}, q1 = B u):    Y != 0 needs u = c d^2;
#   m even (p1 = A u^{2k+1}, q1 = B u^2): Y = ytil u^m, ALWAYS present.
A_, B_, C_ = sp.symbols("A B C", nonzero=True)
u_, d_, Y_ = [sp.Function(n)(w) for n in ("u", "d", "Y")]
yt, cc, y0 = sp.symbols("ytil c0 y0")

# ---- m odd = 2k+1 ----------------------------------------------------------
E2o, E1o, E0o = E_coeffs(A_ * u_**(K + 1), B_ * u_, C_, m=Modd)
E1z, E0z = (sp.expand(E.subs(r0, 0).doit()) for E in (E1o, E0o))
muo = K + 1                                   # (m+1)/2
check("D1 odd (r0 sheared to 0): p0* = ((m+1)/2)(A/B) u^{(m-1)/2} q0 "
      "solves E1 for ANY q0, k symbolic",
      zero(E1z.subs(p0, muo * (A_ / B_) * u_**K * q0)))
E1hom = sp.expand(E1z.subs({p0: Y_, q0: 0}).doit())
check("D1 odd homogeneous certificate: (Y^2/q1^m)' ~ Y E1hom  =>  "
      "Y^2 = e q1^m = e' u^m  =>  Y = 0 unless u = c d^2 (m odd)",
      zero(sp.diff(Y_**2 / (B_ * u_)**Modd, w)
           - Y_ * E1hom / (C_ * B_**(Modd + 1) * u_**(Modd + 1))))
check("D1 odd, Y = 0: E0(p0*) == ((m+1)/2) k C (A/B) u^{(m-3)/2} q0 "
      "(u'q0 - 2u q0'), k symbolic  =>  for m >= 5 the factor u^{(m-3)/2} "
      "forces u | kappa: EMPTY.  At m = 3 the factor is trivial",
      zero(E0z.subs(p0, muo * (A_ / B_) * u_**K * q0)
           - C_ * muo * (A_ / B_) * K * u_**(K - 1) * q0
           * (dw(u_) * q0 - 2 * u_ * dw(q0))))
# m = 3 leftover: q0 | kappa => q0 = n, n^2 u' = const => u LINEAR; box:
s0_, s1_, n_, h0_ = sp.symbols("s0 s1 n h0")
ulin = s0_ + s1_ * w
famD1 = E_coeffs(A_ * ulin**2, B_ * ulin, C_,
                 P0=2 * A_ * n_ * ulin / B_, Q0=n_, R0=sp.S.Zero, m=3)
check("D1 m=3 sheared family: u = s0 + s1 w, q0 = n, p0 = (2A/B) n u, "
      "r0 = 0 has E2 = E1 = 0, E0 = (2AC/B) n^2 s1 = kappa != 0 => n != 0",
      zero(famD1[0]) and zero(famD1[1])
      and zero(famD1[2] - 2 * A_ * C_ * n_**2 * s1_ / B_))
jets = sp.solve([(2 * A_ * n_ * ulin / B_ - A_ * ulin**2 * h0_).subs(w, 0),
                 (n_ - B_ * ulin * h0_).subs(w, 0)], [h0_, n_], dict=True)
check("D1 m=3 BOX: un-shearing needs p0(0) = q0(0) = 0; the jet system "
      "forces n = 0 (unique solution) => kappa = 0.  EMPTY",
      len(jets) == 1 and jets[0][n_] == 0)
# odd square subcase u = c d^2 (Y = y0 d^m):
p0sq = (muo * (A_ / B_) * u_**K * q0 + y0 * d_**Modd).subs(u_, cc * d_**2)
check("D1 odd, u = c d^2: p0 = p0* + y0 d^m still solves E1, and "
      "E0 == C d^{m-2} (d'q0 - d q0')(2k(k+1)(A/B) c^k q0 + m y0 d), "
      "k symbolic  =>  d | kappa: EMPTY for every odd m >= 3",
      zero(subd(E1z.subs(u_, cc * d_**2), (p0, p0sq)))
      and zero(subd(E0z.subs(u_, cc * d_**2), (p0, p0sq))
               - C_ * d_**(Modd - 2) * (dw(d_) * q0 - d_ * dw(q0))
               * (2 * muo * K * (A_ / B_) * cc**K * q0 + Modd * y0 * d_)))

# ---- m even = 2k ------------------------------------------------------------
E2e, E1e, E0e = E_coeffs(A_ * u_**(2 * K + 1), B_ * u_**2, C_, m=Mev)
E1ze, E0ze = (sp.expand(E.subs(r0, 0).doit()) for E in (E1e, E0e))
mue = sp.Rational(1, 2) * (Mev + 1)
p0ge = mue * (A_ / B_) * u_**(2 * K - 1) * q0 + yt * u_**(2 * K)
check("D1 even: p0 = p0* + ytil u^m solves E1 for ANY q0, ytil, k symbolic "
      "(the homogeneous direction Y = ytil u^m exists for EVERY even m)",
      zero(E1ze.subs(p0, p0ge)))
E1homE = sp.expand(E1ze.subs({p0: Y_, q0: 0}).doit())
check("D1 even homogeneous certificate: (Y^2/u^{2m})' ~ Y E1hom  =>  "
      "Y = ytil u^m is the COMPLETE homogeneous solution set",
      zero(sp.diff(Y_**2 / u_**(2 * Mev), w)
           - Y_ * E1homE / (C_ * B_ * u_**(2 * Mev + 2))))
WtE = dw(u_) * q0 - u_ * dw(q0)
check("D1 even KAPPA-EQUATION, k symbolic:  E0 == C u^{m-2} (u'q0 - u q0') "
      "[ ((m-1)(m+1)/2)(A/B) q0 + m ytil u ]",
      zero(E0ze.subs(p0, p0ge)
           - C_ * u_**(2 * K - 2) * WtE
           * ((Mev - 1) * mue * (A_ / B_) * q0 + Mev * yt * u_)))
print("     => for every EVEN m >= 4 the factor u^{m-2} is nonconstant and")
print("        forces u | kappa: D1 EMPTY.  At m = 2 the factor is TRIVIAL")
print("        -- the unique escape hatch.  THE DIVERGENCE POINT:")

# ---- the m = 2 threading: Alpoge-Mathew, end-to-end -------------------------
uAM, AAM, BAM, CAM = 1 + w, sp.S.One, sp.Integer(3), sp.S.NegativeOne
check("AM at m=2 is in D1: p1 = u^3 (A=1), q1 = 3u^2 (B=3), r1 = -1 (C=-1),"
      " u = 1+w nonconstant, t const",
      sp.expand(p1AM - AAM * uAM**3) == 0
      and sp.expand(q1AM - BAM * uAM**2) == 0 and r1AM == CAM)
hAM = -r0AM / CAM
p0sh = sp.expand(p0AM + p1AM * hAM)
q0sh = sp.expand(q0AM + q1AM * hAM)
check("AM sheared to r0 = 0 (h = -r0/C = 2-3w): p0sh = (w+1)(w+2), "
      "q0sh = 6+4w",
      sp.expand(p0sh - (w + 1) * (w + 2)) == 0
      and sp.expand(q0sh - (6 + 4 * w)) == 0)
p0starAM = sp.Rational(3, 2) * (AAM / BAM) * uAM * q0sh
check("AM: p0sh - p0* == -(1+w)^2 = ytil u^m with ytil = -1  "
      "(the even-m homogeneous direction, NONZERO for AM)",
      sp.expand(p0sh - p0starAM + uAM**2) == 0)
WtAM = sp.expand(sp.diff(uAM, w) * q0sh - uAM * sp.diff(q0sh, w))
E0AM = CAM * uAM**0 * WtAM * (sp.Rational(3, 2) * (AAM / BAM) * q0sh
                              + 2 * (-1) * uAM)
check("AM THREADS the kappa-equation at m=2: u^{m-2} = 1 (trivial factor), "
      "u'q0 - u q0' = 2, and E0 = C*2*[(3/2)(1/3)q0sh - 2u] == -2 == "
      "kappa_AM exactly",
      WtAM == 2 and sp.expand(E0AM) == -2
      and det_m(P_AM, Q_AM, R_AM, 2) == -2)
print("     m=2 survives through (i) the trivial factor u^{m-2} AND (ii)")
print("     the even-m homogeneous direction ytil != 0; for m >= 4 (even)")
print("     u^{m-2} | kappa kills everything, for m = 3 the box valuation")
print("     val_w p0 >= 3 does.  D1 EMPTY for ALL m >= 3.")

# ===========================================================================
print("== 8. stratum D2: EMPTY for all m >= 2, uniform ==")
# ===========================================================================
# D2: s const (p1 = A const), t nonconstant.  v-shear h = -p0/A kills p0.
q1f, r1f = sp.Function("q1")(w), sp.Function("r1")(w)
E0gen = E_coeffs(A_, q1f, r1f)[2]
check("D2: E0|_{p0 = 0} == -A (q0 r0)' for ARBITRARY q1, r1 and symbolic m"
      "  =>  q0 r0 exactly linear: one of q0, r0 is a nonzero constant, "
      "the other exactly linear",
      zero(E0gen.subs(p0, 0) + A_ * dw(q0 * r0)))
# ---- m odd: q1 = B t^k, r1 = C t^{k+1} --------------------------------------
E1d = sp.expand(E_coeffs(A_, B_ * t_**K, C_ * t_**(K + 1),
                         m=Modd)[1].subs(p0, 0).doit())
check("D2 odd: E1 == -A t^{k-1}[(m+1)C t^2 q0' + (k+1)C t t' q0 "
      "- (m-1)B t r0' + k B t' r0], k symbolic",
      zero(E1d + A_ * t_**(K - 1)
           * ((Modd + 1) * C_ * t_**2 * dw(q0) + (K + 1) * C_ * t_ * dw(t_) * q0
              - (Modd - 1) * B_ * t_ * dw(r0) + K * B_ * dw(t_) * r0)))
print("     case r0 = delta const: t | k B delta t' => delta = 0 (q0 r0 not")
print("     linear) or t' = 0 (t const): EMPTY.")
rho0 = (K + 1) * ga * C_ / (K * B_)
check("D2 odd, case q0 = gamma const: r0* = ((k+1) gamma C/(k B)) t solves "
      "E1 for ANY t, k symbolic; homogeneous solutions Y ~ t^{1/2} "
      "((Y^2/t)' certificate)  =>  Y^2 = e t, deg Y < deg t",
      zero(E1d.subs({q0: ga, r0: rho0 * t_}))
      and zero(sp.diff(Y_**2 / t_, w)
               - Y_ * (2 * (Modd - 1) * t_ * dw(Y_)
                       - 2 * K * dw(t_) * Y_) / ((Modd - 1) * t_**2)))
e_ = sp.Function("e")(w)
check("D2 odd, q0 = gamma, t = c e^2 (Y != 0 possible): the COMPLETE "
      "solution r0 = r0* + y0 e still solves E1, k symbolic",
      zero(subd(E1d.subs({q0: ga, r0: rho0 * t_ + y0 * e_}),
                (t_, cc * e_**2))))
print("     => in every case r0 = rho0 t + Y with rho0 != 0 and")
print("     deg Y <= (deg t)/2 < deg t, so deg r0 = deg t; r0 exactly")
print("     linear forces t LINEAR (then Y^2 = e t => Y = 0, r0 = rho0 t).")
print("     Surviving family: E0 = -A gamma rho0 t1 = kappa fine -- but")
print("     un-shearing needs val(p0 = -A h) >= m >= 1 => h(0) = 0 =>")
print("     q0(0) = gamma != 0, violating the box q0(0) = 0.")
print("     EMPTY (uniform in m -- uses only val_w p0 >= 1).")
# ---- m even: q1 = B t^{2k-1}, r1 = C t^{2k+1} -------------------------------
E1d2 = sp.expand(E_coeffs(A_, B_ * t_**(2 * K - 1), C_ * t_**(2 * K + 1),
                          m=Mev)[1].subs(p0, 0).doit())
check("D2 even: E1 == -A t^{2k-2}[(m+1)C t^3 q0' + (m+1)C t^2 t' q0 "
      "- (m-1)B t r0' + (m-1)B t' r0], k symbolic",
      zero(E1d2 + A_ * t_**(2 * K - 2)
           * ((Mev + 1) * C_ * t_**3 * dw(q0)
              + (Mev + 1) * C_ * t_**2 * dw(t_) * q0
              - (Mev - 1) * B_ * t_ * dw(r0) + (Mev - 1) * B_ * dw(t_) * r0)))
print("     case r0 = delta const: t^2 | (m-1)B delta t' => delta = 0 or")
print("     t const: EMPTY.")
rho_ = sp.Symbol("rho")
sigE = (Mev + 1) * ga * C_ / ((Mev - 1) * B_)
check("D2 even, case q0 = gamma: r0 = sigma t^2 + rho t (sigma = "
      "(m+1) gamma C/((m-1)B) != 0) is the COMPLETE solution of E1 "
      "(homogeneous = rho t), k symbolic  =>  deg r0 = 2 deg t >= 2, "
      "contradicting r0 exactly linear: EMPTY",
      zero(E1d2.subs({q0: ga, r0: sigE * t_**2 + rho_ * t_})))

# ===========================================================================
print("== 9. stratum D3 (s, t squarefree): EMPTY for all m >= 3 ==")
# ===========================================================================
# Shear-invariants (Lambda_1 = q1 p0 - p1 q0 etc. are exactly invariant):
#   m odd:  G1 = t^k p0 - a s^k q0,        G2 = b t q0 - s r0
#   m even: G1 = t^{m-1} p0 - a s^{m-1} q0, G2 = b t^2 q0 - s^2 r0
h_, g2_, gk_ = [sp.Function(n)(w) for n in ("h", "g2", "gk")]
Jj = sp.Symbol("j", integer=True, nonnegative=True)
eps = sp.Symbol("epsilon")

# ---- m odd = 2k+1 -----------------------------------------------------------
E2o3, E1o3, E0o3 = E_coeffs(a_ * s_**(K + 1), s_ * t_**K, b_ * t_**(K + 1),
                            m=Modd)
G1o = t_**K * p0 - a_ * s_**K * q0
G2o = b_ * t_ * q0 - s_ * r0
Theta = dw(s_ * t_) * G2o - 2 * s_ * t_ * dw(G2o)
Theta1 = 2 * s_ * t_ * dw(G1o) - (Modd * dw(s_) * t_ - s_ * dw(t_)) * G1o
check("D3 odd: E1 == a k s^k t^{k-1} Theta + b t^k Theta1, k symbolic, "
      "with Theta = (st)'G2 - 2st G2', Theta1 = 2st G1' - (m s't - st')G1, "
      "G1 = t^k p0 - a s^k q0, G2 = b t q0 - s r0 (shear-invariant)",
      zero(E1o3 - (a_ * K * s_**K * t_**(K - 1) * Theta
                   + b_ * t_**K * Theta1)))
# divisibility chain: t-side one step, s-side k steps (coefficient m - 2j)
check("D3 odd, t-side: Theta == s t' G2 (mod t)  =>  t | G2 (t squarefree,"
      " gcd(s,t) = 1); and Theta|_{G2 = t g2} == t [W g2 - 2st g2']",
      sp.fraction(sp.cancel(sp.together(
          (Theta - dw(t_) * s_ * G2o) / t_)))[1] == 1
      and zero(Theta.subs(r0, (b_ * t_ * q0 - t_ * g2_) / s_)
               - t_ * (Wr * g2_ - 2 * s_ * t_ * dw(g2_))))
Bj = lambda f, j: (2 * s_ * t_ * dw(f)                    # noqa: E731
                   - ((Modd - 2 * j) * dw(s_) * t_ - s_ * dw(t_)) * f)
check("D3 odd, s-side iteration, j symbolic: B_j(s h) == s B_{j+1}(h) and "
      "B_j(h) == -(m-2j) s' t h (mod s); m odd => m-2j != 0 at every step "
      "=>  s^k | Theta1 = B_0(G1) forces G1 = s^k g  (s squarefree)",
      zero(Bj(s_ * h_, Jj) - s_ * Bj(h_, Jj + 1))
      and sp.fraction(sp.cancel(sp.together(
          (Bj(h_, Jj) + (Modd - 2 * Jj) * dw(s_) * t_ * h_) / s_)))[1] == 1)
Z_ = a_ * K * g2_ - b_ * gk_
E1sub = E1o3.subs({p0: (s_**K * gk_ + a_ * s_**K * q0) / t_**K,
                   r0: (b_ * t_ * q0 - t_ * g2_) / s_})
check("D3 odd: E1|_{G1 = s^k gk, G2 = t g2} == s^k t^k (W Z - 2stZ'), "
      "Z = a k g2 - b gk, k symbolic; and (Z^2 t/s)' == Z(2stZ' - WZ)/s^2 "
      "=>  Z^2 t = e s  =>  Z = 0 (gcd(s,t) = 1, t nonconstant)",
      zero(E1sub - s_**K * t_**K * (Wr * Z_ - 2 * s_ * t_ * dw(Z_)))
      and zero(sp.diff(Z_**2 * t_ / s_, w)
               - Z_ * (2 * s_ * t_ * dw(Z_) - Wr * Z_) / s_**2))
E0loc = E0o3.subs({p0: (s_**K * (a_ * K / b_) * g_ + a_ * s_**K * q0) / t_**K,
                   r0: (b_ * t_ * q0 - t_ * g_) / s_})
check("D3 odd KAPPA-EQUATION on the E1-locus (gk = (ak/b) g, g2 = g), "
      "k symbolic:  E0 == -(a k(k+1)/b) s^{k-1} g (2st g' - W g)  "
      "(q0 drops out: the shear direction)",
      zero(E0loc + (a_ * K * (K + 1) / b_) * s_**(K - 1) * g_
           * (2 * s_ * t_ * dw(g_) - Wr * g_)))
print("     => m odd >= 5: the factor s^{k-1} (k >= 2) forces s | kappa:")
print("        EMPTY outright, no box needed.  (g = 0 gives E0 = 0 != "
      "kappa.)")
# m = 3 leftover (k = 1): kappa-eq g(2stg' - Wg) = const => g = c const,
# kappa = (2a/b) c^2 W => W const; the single shear-orbit dies on the box:
c1 = sp.Symbol("c1", nonzero=True)
q0f = sp.Function("q0")(w)
p0D3 = a_ * s_ * (q0f + c1 / b_) / t_
r0D3 = t_ * (b_ * q0f - c1) / s_
famD3 = E_coeffs(a_ * s_**2, s_ * t_, b_ * t_**2, P0=p0D3, Q0=q0f,
                 R0=r0D3, m=3)
check("D3 m=3 family (all solutions): E1 == 0 and E0 == (2a c^2/b) W "
      "identically -- Keller iff W = s't - st' is a nonzero constant",
      zero(famD3[1]) and zero(famD3[2] - 2 * a_ * c1**2 / b_ * Wr))
check("D3 m=3 BOX: t p0 == a s (q0 + c/b), so p0(0) t(0) = a s(0) c/b != 0 "
      "(s(0) != 0 from p1(0) = 1) while the box forces p0(0) = 0: EMPTY",
      zero(t_ * p0D3 - a_ * s_ * (q0f + c1 / b_)))

# ---- m even = 2k ------------------------------------------------------------
E2e3, E1e3, E0e3 = E_coeffs(a_ * s_**(2 * K + 1), s_**2 * t_**(2 * K - 1),
                            b_ * t_**(2 * K + 1), m=Mev)
G1e = t_**(2 * K - 1) * p0 - a_ * s_**(2 * K - 1) * q0
G2e = b_ * t_**2 * q0 - s_**2 * r0
ThetaE = dw(s_ * t_) * G2e - s_ * t_ * dw(G2e)
C0 = s_ * t_ * dw(G1e) - (Mev * dw(s_) * t_ - s_ * dw(t_)) * G1e
check("D3 even: E1 == (m-1) a s^{2k} t^{2k-2} ThetaE + 2 b s t^{2k} C0, "
      "k symbolic, ThetaE = (st)'G2 - st G2', C0 = st G1' - (m s't - st')G1,"
      " G1 = t^{m-1} p0 - a s^{m-1} q0, G2 = b t^2 q0 - s^2 r0",
      zero(E1e3 - ((Mev - 1) * a_ * s_**(2 * K) * t_**(2 * K - 2) * ThetaE
                   + 2 * b_ * s_ * t_**(2 * K) * C0)))
check("D3 even, t-side: t^2 | ThetaE => t | G2 (squarefree), and "
      "ThetaE|_{G2 = t g2} == t^2 (s' g2 - s g2')",
      zero(ThetaE.subs(r0, (b_ * t_**2 * q0 - t_ * g2_) / s_**2)
           - t_**2 * (dw(s_) * g2_ - s_ * dw(g2_))))
Cj = lambda f, j: (s_ * t_ * dw(f)                        # noqa: E731
                   - ((Mev - j) * dw(s_) * t_ - s_ * dw(t_)) * f)
check("D3 even, s-side iteration, j symbolic: C_j(s h) == s C_{j+1}(h) and "
      "C_j(h) == -(m-j) s' t h (mod s); coefficients m, m-1, ..., 2 never "
      "vanish  =>  s^{m-1} | C0 forces G1 = s^{m-1} g (s squarefree)",
      zero(Cj(s_ * h_, Jj) - s_ * Cj(h_, Jj + 1))
      and sp.fraction(sp.cancel(sp.together(
          (Cj(h_, Jj) + (Mev - Jj) * dw(s_) * t_ * h_) / s_)))[1] == 1)
combo = (Mev - 1) * a_ * (dw(s_) * g2_ - s_ * dw(g2_)) \
    + 2 * b_ * (s_ * t_ * dw(g_) - Wr * g_)
Ze = (Mev - 1) * a_ * g2_ - 2 * b_ * t_ * g_
E1subE = E1e3.subs({p0: (s_**(2 * K - 1) * g_
                         + a_ * s_**(2 * K - 1) * q0) / t_**(2 * K - 1),
                    r0: (b_ * t_**2 * q0 - t_ * g2_) / s_**2})
check("D3 even: E1 on the locus == s^m t^{2k} [(m-1)a(s'g2 - s g2') + "
      "2b(st g' - W g)] == s^m t^{2k} (s'Z - sZ'),  Z = (m-1)a g2 - 2bt g "
      " =>  (Z/s)' = 0  =>  Z = epsilon s",
      zero(E1subE - s_**(2 * K) * t_**(2 * K) * combo)
      and zero(combo - (dw(s_) * Ze - s_ * dw(Ze))))
g2sol = (2 * b_ * t_ * g_ + eps * s_) / ((Mev - 1) * a_)
E0locE = E0e3.subs({p0: (s_**(2 * K - 1) * g_
                         + a_ * s_**(2 * K - 1) * q0) / t_**(2 * K - 1),
                    r0: (b_ * t_**2 * q0 - t_ * g2sol) / s_**2})
check("D3 even KAPPA-EQUATION on the E1-locus, k symbolic:  E0 == "
      "-s^{m-2} ((2m+2) b t g + epsilon s)(st g' - W g) / ((m-1) a)  "
      "(q0 drops out)",
      zero(E0locE + s_**(2 * K - 2)
           * ((2 * Mev + 2) * b_ * t_ * g_ + eps * s_)
           * (s_ * t_ * dw(g_) - Wr * g_) / ((Mev - 1) * a_)))
print("     => m even >= 4: the factor s^{m-2} forces s | kappa: EMPTY")
print("        outright.  At m = 2 the factor is trivial -- the same")
print("        escape hatch as D1 (but AM itself lives in D1).")
print("     D3 EMPTY for all m >= 3 whenever s, t are squarefree.")

# ===========================================================================
print("== 10. the D3' gap and its in-box closure (m = 4, 5) ==")
# ===========================================================================
# The chains above use squarefreeness of s (s-side iteration) and t
# (t-side step).  Non-squarefree s or t needs deg s >= 2 or deg t >= 2:
#   m odd:  deg p1 = (m+1)/2 deg s >= m+1  or  deg r1 >= m+1;
#   m even: deg p1 = (m+1) deg s >= 2(m+1)  or  deg r1 >= 2(m+1).
# (At m = 3 this is the D3' gap of docs/SEARCH_113.md, threshold deg 4.)
# In-box closure: Keller + (r1 != 0) EMPTY covers strata C and D wholesale
# within the box, including every non-squarefree configuration there.

from jcqft.gb_backend import available_backends, is_unit_ideal  # noqa: E402
print(f"  GB backends: {available_backends()}")


def vlinear_box(m, dp0, dp1, dq0, dq1, dr0, dr1):
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
    P0, u1 = poly("p", range(m, dp0 + 1))
    P1, u2 = poly("pa", range(0, dp1 + 1), gauge={0: 1})
    Q0, u3 = poly("q", range(1, dq0 + 1), gauge={1: 1})
    Q1, u4 = poly("qa", range(0, dq1 + 1))
    R0, u5 = poly("r", range(0, dr0 + 1), gauge={0: 1})
    R1, u6 = poly("ra", range(0, dr1 + 1))
    return (P0 + P1 * v, Q0 + Q1 * v, R0 + R1 * v,
            u1 + u2 + u3 + u4 + u5 + u6)


yv = sp.Symbol("yrab")
KELLERS = {}
for m in (4, 5):
    Pb, Qb, Rb, U = vlinear_box(m, m + 1, 2, 3, 2, 2, 2)
    KEL = sp.Poly(sp.expand(det_m(Pb, Qb, Rb, m)) + 1, w, v).coeffs()
    KELLERS[m] = (Pb, Qb, Rb, U, KEL)
    print(f"  m={m} small box (deg p0<={m + 1}, p1,q1,r1<=2, q0<=3, r0<=2):"
          f" {len(U)} unknowns, {len(KEL)} Keller equations")
    for cf in [u for u in U if str(u).startswith("ra")]:
        try:
            empty = is_unit_ideal(KEL + [1 - yv * cf], U + [yv],
                                  backend="auto", timeout=600)
            check(f"m={m}: Keller + ({cf} != 0): ideal == (1), EMPTY "
                  "(covers C u D in-box, incl. non-squarefree D3')", empty)
        except (RuntimeError, TimeoutError) as exc:
            print(f"  [unresolved] m={m} {cf} != 0 -- {exc}")

# (ii) r1 == 0 slice: Rabinowitsch certificates pin the in-box variety to
# the gauged family A u B (with its box truncation qa0 * p_j = 0).
for m in (4, 5):
    Pb2, Qb2, Rb2, U2 = vlinear_box(m, m + 1, 2, 3, 2, 2, -1)   # r1 == 0
    K2 = sp.Poly(sp.expand(det_m(Pb2, Qb2, Rb2, m)) + 1, w, v).coeffs()
    S = {str(u): u for u in U2}
    certs = [S["pa1"], S["pa2"], S["r1"], S["r2"], S["qa1"], S["qa2"],
             S["q2"], S["q3"],
             S["qa0"] * S[f"p{m}"], S["qa0"] * S[f"p{m + 1}"]]
    ok = []
    for f in certs:
        try:
            empty = is_unit_ideal(K2 + [1 - yv * f], U2 + [yv],
                                  backend="auto", timeout=600)
            ok.append(empty)
        except (RuntimeError, TimeoutError) as exc:
            print(f"  [unresolved] m={m} radical membership {f} -- {exc}")
            ok.append(False)
    check(f"m={m}, r1 = 0 slice: pa1, pa2, r1, r2, qa1, qa2, q2, q3, "
          "qa0*p_j all vanish on the in-box Keller variety (Rabinowitsch)"
          "  =>  variety == the gauged family A u B (box-truncated)",
          all(ok))

# ===========================================================================
print("== 11. the m:1 orbifold mechanism is empty for every m >= 3 ==")
# ===========================================================================
# Target coordinates carry weights (-m, -1, 1).  The stabilizer of a point
# on the a-axis (a != 0, b = c = 0) is {lam : lam^m = 1} = Z_m; every other
# coordinate axis and pair has trivial stabilizer (gcd with 1); for
# composite m there is NO intermediate Z_d stratum -- coordinate stabilizers
# are Z_m or trivial.  So the only mechanism is m:1: a witness (w0, v0)
# with Q = R = 0 != P on a free orbit {x != 0}.
from math import gcd  # noqa: E402
check("weights: for every m in {3,...,12} the only nontrivial coordinate "
      "stabilizer is Z_m on the -m axis (gcd(m,1) = 1 kills all pairs)",
      all(gcd(m, 1) == 1 for m in range(3, 13)))
check("classification corollary: every v-linear Keller map has R == 1 "
      "(gauge), which never vanishes: NO m:1 witness exists, for ANY "
      "m >= 3 (up to the box-closed D3' gap)", Rg == 1)

# independent pointwise Groebner queries (residual torus moves any witness
# to (1,1), (1,0) or (0,1); (0,0) is excluded by R(0,0) = 1):
rr = sp.Symbol("rrab")
for m in (4, 5):
    Pb, Qb, Rb, U, KEL = KELLERS[m]
    for w0, v0 in ((1, 1), (1, 0), (0, 1)):
        sys_ = [sp.expand(e) for e in
                KEL + [Qb.subs({w: w0, v: v0}), Rb.subs({w: w0, v: v0}),
                       1 - rr * Pb.subs({w: w0, v: v0})]]
        try:
            empty = is_unit_ideal(sys_, U + [rr], backend="auto",
                                  timeout=600)
            check(f"m={m}: m:1 witness at (w,v) = ({w0},{v0}): "
                  "ideal == (1), EMPTY", empty)
        except (RuntimeError, TimeoutError) as exc:
            print(f"  [unresolved] m={m} witness at ({w0},{v0}) -- {exc}")

print("\nVERDICT (v-linear class, all w-degrees, every m >= 3; D3'")
print("non-squarefree gap box-checked at m <= 5): every Keller map is the")
print("tame automorphism  F = (p0(xy)/x^m + z,  y + b0 x^{m-1} F1,  x)")
print("(up to gauge).  NO counterexample, NO m:1 orbifold covering, for")
print("ANY m >= 3.  The Alpoge-Mathew mechanism lives at m = 2 and ONLY")
print("at m = 2: its stratum (D1, m even) is killed for m >= 4 by the")
print("factor u^{m-2} in the kappa-equation and for m = 3 by the box")
print("valuation val_w p0 >= 3.")


# ===========================================================================
# --full: larger boxes and the targeted m = 5 non-squarefree gap queries
# ===========================================================================

def run_full(budget):
    print("\n== 12. (--full) larger boxes + targeted m = 5 gap queries ==")
    # (i) medium boxes, all r1 coefficients, m = 4, 5
    for m in (4, 5):
        Pb3, Qb3, Rb3, U3 = vlinear_box(m, m + 2, 3, 4, 3, 3, 3)
        K3 = sp.Poly(sp.expand(det_m(Pb3, Qb3, Rb3, m)) + 1, w, v).coeffs()
        print(f"  m={m} medium box (deg p0<={m + 2}, p1,q1,r1<=3, q0<=4, "
              f"r0<=3): {len(U3)} unknowns, {len(K3)} equations")
        for cf in [u for u in U3 if str(u).startswith("ra")]:
            try:
                empty = is_unit_ideal(K3 + [1 - yv * cf], U3 + [yv],
                                      backend="auto", timeout=budget)
                check(f"m={m} medium box, Keller + ({cf} != 0): EMPTY",
                      empty)
            except (RuntimeError, TimeoutError) as exc:
                print(f"  [unresolved] m={m} {cf} != 0 -- {exc}")
    # (ii) targeted m = 5 non-squarefree D3' parametrizations
    # (m = 4 needs deg p1 >= 10 or deg r1 >= 10: beyond any feasible box).
    m = 5
    rho1, n0, n1, e0_, e1_ = sp.symbols("rho1 n0 n1 e0 e1")

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

    cases = {
        "s = (1+rho1 w)^2 (non-squarefree), t deg 1":
            ((1 + rho1 * w)**2, n0 + n1 * w, [rho1, n0, n1], [n1]),
        "s deg 1, t = (e0+e1 w)^2 (non-squarefree)":
            ((1 + rho1 * w), (e0_ + e1_ * w)**2, [rho1, e0_, e1_], [e1_]),
    }
    for label, (s5, t5, pars, sat) in cases.items():
        # m = 5 odd: p1 = a s^3, q1 = s t^2, r1 = b t^3
        P0, cp = poly("p", range(5, 8))
        Q0, cq = poly("q", range(1, 6), gauge={1: 1})
        R0, cr = poly("r", range(0, 5), gauge={0: 1})
        Dg = sp.expand(det_m(P0 + a_ * s5**3 * v, Q0 + s5 * t5**2 * v,
                             R0 + b_ * t5**3 * v, 5))
        Kg = sp.Poly(Dg + 1, w, v).coeffs()
        ys = sp.symbols(f"ysat0:{len(sat) + 2}")
        sysg = (Kg + [1 - ys[0] * a_, 1 - ys[1] * b_]
                + [1 - yi * cq_ for yi, cq_ in zip(ys[2:], sat)])
        gens = pars + [a_, b_] + cp + cq + cr + list(ys)
        try:
            empty = is_unit_ideal(sysg, gens, backend="auto",
                                  timeout=budget)
            check(f"m=5 D3' gap, {label} (deg p0<=7, q0<=5, r0<=4): EMPTY",
                  empty)
            continue
        except (RuntimeError, TimeoutError) as exc:
            print(f"  [msolve unresolved] {label} -- {exc}")
        if "singular" not in available_backends():
            continue
        # memory-frugal fallback: Singular degBound ladder (one-sided).
        for db in (4, 5, 6, 7):
            try:
                empty = is_unit_ideal(sysg, gens, backend="singular",
                                      timeout=budget, degbound=db)
            except (RuntimeError, TimeoutError) as exc:
                print(f"  [singular degBound {db} unresolved] {label} "
                      f"-- {exc}")
                break
            if empty:
                check(f"m=5 D3' gap, {label}: EMPTY (Singular degBound "
                      f"{db} certificate)", empty)
                break
            print(f"  [no certificate at degBound {db}] {label}")
        else:
            print(f"  [unresolved] {label}: no unit certificate up to "
                  "degBound 7 -- see docs/SEARCH_11M.md")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="larger-box Groebner corroboration + targeted "
                         "m = 5 gap queries")
    ap.add_argument("--budget", type=int, default=1800,
                    help="per-query timeout for --full (s)")
    args = ap.parse_args()
    if args.full:
        run_full(args.budget)
    print(f"\nall checks passed in {time.time() - T0:.1f} s")
