"""Exact reduced Keller identity for the weight systems (1,-1,-m), m >= 1.

(OPEN_QUESTIONS B3, step 1: before any new-counterexample search in the
weight system (1,-1,-3) — or any (1,-1,-m) — the dimensional reduction of
the Keller condition must be exact.  Module: jcqft/reduction_w.py;
write-up folded into docs/NEW_COUNTEREXAMPLES.md.)

Statements proved here, all exact:

1.  GENERIC-FUNCTION PROOF (m = 1, 2, 3, 4): for arbitrary smooth P, Q, R
    of the invariants w = x*y, v = x^m*z, the equivariant map
        F = (P/x^m, Q/x, x*R)
    has
        det DF  =  det M  =  -m*P*J2(Q,R) + Q*J2(P,R) + R*J2(P,Q),
    a function of (w, v) ALONE, where M = [[d_i S_i, S_i_w, S_i_v]],
    d = (-m,-1,1), S = (P,Q,R).  (Sympy applied to undetermined functions
    P(w,v), Q(w,v), R(w,v): a bona fide identity of differential
    polynomials, not a spot check.)
2.  COMPACT FORM:  R^m * det M = J2(P*R^m, Q*R)  identically in
    undetermined functions.  Hence the REDUCED KELLER IDENTITY
        det DF = kappa   <=>   J2(P*R^m, Q*R) = kappa*R^m
    (m = 2 is exactly jcqft.reduction's identity).
3.  CONSISTENCY, m = 2: reduction_w reproduces jcqft.reduction on the
    Alpöge–Mathew map (same P, Q, R; residual 0 at kappa = -2).
4.  POLYNOMIALITY BOXES: assemble() of every admissible ansatz monomial is
    polynomial in (x,y,z), and monomials violating j + m*k >= m (P) or
    >= 1 (Q) assemble to non-polynomial maps (checked exhaustively in a
    box for m = 3).
5.  m = 3 SANITY: a random polynomial triple (P, Q, R) in the admissible
    box has det DF (3D, symbolic) == det M (2D) == J2(P R^3, Q R)/R^3.
"""

import time

import sympy as sp

from jcqft.reduction_w import (assemble, det_m, extract, j2, keller_residual,
                               monomial_box, v, w, x, y, z)

t0 = time.time()


def check(label, cond):
    assert cond, label
    print(f"  [ok] {label}   ({time.time() - t0:.1f} s)")


print("== 1+2. generic-function proof of the reduction, m = 1..4 ==")
Pf, Qf, Rf = [f(w, v) for f in sp.symbols("P Q R", cls=sp.Function)]
for m in (1, 2, 3, 4):
    sub = {w: x * y, v: x**m * z}
    F3 = (Pf.subs(sub) / x**m, Qf.subs(sub) / x, x * Rf.subs(sub))
    J = sp.Matrix([[sp.diff(Fi, u) for u in (x, y, z)] for Fi in F3])
    detDF = sp.simplify(J.det())
    # back to invariant coordinates: eliminate (y, z) via (w, v)
    detDF_wv = sp.simplify(detDF.subs({y: w / x, z: v / x**m}))
    dM = -m * Pf * j2(Qf, Rf) + Qf * j2(Pf, Rf) + Rf * j2(Pf, Qf)
    check(f"m={m}: det DF == det M  (function of (w,v) alone)",
          sp.simplify(detDF_wv - dM) == 0 and not detDF_wv.has(x))
    check(f"m={m}: R^m * det M == J2(P R^m, Q R)  (compact form)",
          sp.simplify(Rf**m * dM - j2(Pf * Rf**m, Qf * Rf)) == 0)

print("== 3. m = 2 reproduces jcqft.reduction on Alpöge–Mathew ==")
from jcqft.core import F as F_AM  # noqa: E402
from jcqft import reduction as red2  # noqa: E402

P0, Q0, R0 = extract(F_AM, 2)
check("extract (m=2) matches jcqft.reduction.extract",
      (P0, Q0, R0) == (red2.P0, red2.Q0, red2.R0))
check("Keller residual == 0 at kappa = -2",
      keller_residual(P0, Q0, R0, -2, 2) == 0)
check("det_m == -2 exactly", det_m(P0, Q0, R0, 2) == -2)

print("== 4. polynomiality boxes, m = 3 ==")
m = 3


def is_poly_map(F3):
    return all(sp.cancel(Fi).as_numer_denom()[1].is_number for Fi in F3)


boxP = monomial_box("P", m, jmax=4, kmax=2)
boxQ = monomial_box("Q", m, jmax=4, kmax=2)
boxR = monomial_box("R", m, jmax=4, kmax=2)
check(f"P-box ({len(boxP)}) all assemble to polynomials",
      all(is_poly_map(assemble(mon, sp.S.Zero + w * 0 + 1, sp.S.One, m)[:1])
          for mon in boxP))
check(f"Q-box ({len(boxQ)}) all assemble to polynomials",
      all(is_poly_map(assemble(sp.S.Zero, mon, sp.S.One, m)[1:2])
          for mon in boxQ))
rejP = [w**j * v**k for j in range(5) for k in range(3) if j + m * k < m]
rejQ = [w**j * v**k for j in range(5) for k in range(3) if j + m * k < 1]
check("every P-monomial with j + 3k < 3 assembles NON-polynomially",
      all(not is_poly_map(assemble(mon, sp.S.Zero, sp.S.One, m)[:1])
          for mon in rejP))
check("every Q-monomial with j + 3k < 1 assembles NON-polynomially",
      all(not is_poly_map(assemble(sp.S.Zero, mon, sp.S.One, m)[1:2])
          for mon in rejQ))
check("linear part needs: v in P-box, w in Q-box, 1 in R-box",
      v in boxP and w in boxQ and sp.S.One in [sp.sympify(mm) for mm in boxR])

print("== 5. m = 3 random-ansatz consistency (3D det vs reduced forms) ==")
import random  # noqa: E402

random.seed(7)


def rand_poly(box):
    return sp.Add(*[random.randint(-3, 3) * mon for mon in box])


P1, Q1, R1 = rand_poly(boxP) + v, rand_poly(boxQ) + w, rand_poly(boxR) + 1
F1 = assemble(P1, Q1, R1, m)
check("assembled map is polynomial", is_poly_map(F1))
J1 = sp.Matrix([[sp.diff(Fi, u) for u in (x, y, z)] for Fi in F1])
det3 = sp.cancel(J1.det())
dM1 = det_m(P1, Q1, R1, m)
check("det DF (3D) == det M (2D) after substituting invariants",
      sp.simplify(det3 - dM1.subs({w: x * y, v: x**m * z})) == 0)
check("R^3 * det M == J2(P R^3, Q R)",
      sp.expand(R1**3 * dM1 - j2(P1 * R1**3, Q1 * R1)) == 0)

print(f"\nall checks passed in {time.time() - t0:.1f} s")
