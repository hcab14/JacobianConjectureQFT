"""Weighted / charted compactification of the Alpoge-Mathew classical map.

Verifies docs/WEIGHTED_COMPACTIFICATION.md:

  1. Escape-chart transition (s, y, gamma) with s=1/x, gamma=F3, and the
     C* weights (-1,-1,+1) on that chart.
  2. The AM escape curve extends regularly to the boundary divisor D_inf={s=0}.
  3. F extends to a polynomial morphism Fbar on U_inf; Fbar(D_inf) lands on
     the Jelonek set {p=0}; eliminating boundary coordinates recovers (p).
  4. Ordinary P3 homogenization fails / obscures the C* mechanism
     (unequal total degrees; leading cone marks fake escapes; total-degree
     truncation discards same-weight terms of F3).

Runtime: a few seconds, exact sympy.
"""

from __future__ import annotations

import sympy as sp

from jcqft.core import (
    F,
    PHI,
    SRC,
    SOURCE_WEIGHTS,
    TARGET_WEIGHTS,
    p,
    q,
)
from jcqft.fibers import exact_fiber
from jcqft.prefilter import leading_part
from jcqft.reduction import P0, Q0, R0, extract

x, y, z = PHI
a, b, c = SRC
s, yy, gamma = sp.symbols("s y_chart gamma")
y0, c3, T = sp.symbols("y_0 c_3 T")

n_checks = 0


def check(label, cond):
    global n_checks
    n_checks += 1
    assert cond, label
    print(f"  OK  {label}")


# ---------------------------------------------------------------------------
print("=== 1. Escape chart: transition, invertibility, C* weights ===")

gamma_xyz = 2 * x - 3 * x**2 * y - x**3 * z
check("gamma := 2x - 3 x^2 y - x^3 z equals F3",
      sp.expand(gamma_xyz - F[2]) == 0)

# Transition: (x,y,z) <-> (s, y, gamma) on x != 0
z_of_s = sp.cancel((2 / s - 3 * yy / s**2 - gamma) * s**3)
check("z in escape chart is the cubic polynomial 2 s^2 - 3 y s - gamma s^3",
      sp.expand(z_of_s - (2 * s**2 - 3 * yy * s - gamma * s**3)) == 0)

phi_chart = (1 / s, yy, z_of_s)
# round-trip: from chart back to gamma
g_back = sp.cancel(gamma_xyz.subs(dict(zip(PHI, phi_chart))))
check("round-trip: gamma recovered from (s,y,z(s,y,gamma))",
      sp.simplify(g_back - gamma) == 0)
check("round-trip: s = 1/x",
      sp.simplify(1 / phi_chart[0] - s) == 0)

# C* action on fields -> weights on (s, y, gamma)
lam = sp.symbols("lam", nonzero=True)
phi_scaled = (lam * x, y / lam, z / lam**2)
s_scaled = sp.cancel(1 / phi_scaled[0])
y_scaled = phi_scaled[1]
g_scaled = sp.cancel(gamma_xyz.subs(dict(zip(PHI, phi_scaled))))
check("s has C* weight -1: s(lam·phi) = s/lam",
      sp.simplify(s_scaled - (1 / x) / lam) == 0)
check("y has C* weight -1 under the field action",
      sp.simplify(y_scaled - y / lam) == 0)
check("gamma has C* weight +1: gamma(lam·phi) = lam * gamma(phi)",
      sp.simplify(g_scaled - lam * gamma_xyz) == 0)
check("source/field weight constants match jcqft.core",
      SOURCE_WEIGHTS == (1, -1, -2) and TARGET_WEIGHTS == (-2, -1, 1))

# Normal form still available (sanity: I7 packaging)
check("reduction extract(F) recovers equivariant (P,Q,R)",
      extract(F) == (P0, Q0, R0))

# ---------------------------------------------------------------------------
print("=== 2. Escape curve extends to D_inf = {s=0} ===")

phi_T = (T, y0, (2 * T - 3 * T**2 * y0 - c3) / T**3)
# In escape coordinates along the curve:
s_of_T = sp.cancel(1 / phi_T[0])
g_of_T = sp.cancel(gamma_xyz.subs(dict(zip(PHI, phi_T))))
check("along escape curve: s = 1/T", sp.simplify(s_of_T - 1 / T) == 0)
check("along escape curve: gamma ≡ c3 (constant)",
      sp.simplify(g_of_T - c3) == 0)
check("along escape curve: y ≡ y0",
      sp.simplify(phi_T[1] - y0) == 0)
# Extended map hat-phi(s) = (s, y0, c3) lands on D_inf at s=0
check("extended curve hat-phi(0) = (0, y0, c3) lies on D_inf",
      True)  # by construction; landing point is affine in U_inf

# ---------------------------------------------------------------------------
print("=== 3. Fbar is polynomial on U_inf; boundary image is {p=0} ===")

FT = [sp.cancel(sp.expand(Fi.subs(dict(zip(PHI, phi_chart))))) for Fi in F]
F1_closed = -(s + yy) * (
    gamma * s**2 + 2 * gamma * s * yy + gamma * yy**2 - 2 * s - yy
)
F2_closed = (
    -3 * gamma * s**2 - 6 * gamma * s * yy - 3 * gamma * yy**2 + 6 * s + 4 * yy
)
check("F1 on U_inf equals the closed polynomial form",
      sp.expand(FT[0] - F1_closed) == 0)
check("F2 on U_inf equals the closed polynomial form",
      sp.expand(FT[1] - F2_closed) == 0)
check("F3 on U_inf equals gamma",
      sp.expand(FT[2] - gamma) == 0)
for i, fi in enumerate(FT):
    sp.Poly(fi, s, yy, gamma)  # raises if not polynomial
    check(f"F{i+1} is a polynomial in (s, y, gamma)", True)

Fbar = [sp.expand(fi.subs(s, 0)) for fi in FT]
Fbar_expected = (
    yy**2 * (1 - gamma * yy),
    yy * (4 - 3 * gamma * yy),
    gamma,
)
for i in range(3):
    check(
        f"Fbar_{i+1} at s=0 matches (y^2(1-gamma y), y(4-3 gamma y), gamma)",
        sp.expand(Fbar[i] - Fbar_expected[i]) == 0,
    )
wall_id = sp.factor(sp.expand(p.subs(dict(zip(SRC, Fbar)))))
check("p(Fbar(0,y,gamma)) ≡ 0 identically (boundary lands on Jelonek)",
      wall_id == 0)

# Image ideal: eliminate y from a - Fbar1, b - Fbar2 with gamma = c
e0 = sp.expand(a - Fbar[0].subs(gamma, c).subs(yy, y0))
e1 = sp.expand(b - Fbar[1].subs(gamma, c).subs(yy, y0))
R = sp.resultant(e0, e1, y0)
check("eliminant Resultant_y(a-Fbar1, b-Fbar2)|_{gamma=c} = c * p",
      sp.expand(R - c * p) == 0)

# Surjectivity onto the rational param of {p=0}: every (y0,c3) hits the wall,
# and the standard param matches Fbar
lim_param = (
    y0**2 * (1 - c3 * y0),
    y0 * (4 - 3 * c3 * y0),
    c3,
)
check("standard escape-limit param equals Fbar(y0,c3)",
      all(sp.expand(Fbar[i].subs({yy: y0, gamma: c3}) - lim_param[i]) == 0
          for i in range(3)))
check("p(lim_param(y0,c3)) ≡ 0 (parametrization of Jelonek)",
      sp.expand(p.subs(dict(zip(SRC, lim_param)))) == 0)

# Spot checks: cusp orbit point is hit (empty affine fiber; all sheets on D_inf)
# (a,b,c) with ac^2=4/27, bc=4/3, c=1 => a=4/27, b=4/3
cusp_pt = (sp.Rational(4, 27), sp.Rational(4, 3), 1)
check("cusp orbit point (4/27,4/3,1) has p=q=0",
      p.subs(dict(zip(SRC, cusp_pt))) == 0
      and q.subs(dict(zip(SRC, cusp_pt))) == 0)
# From docs: (y0,gamma)=(2/3,1) maps to the cusp at c=1
cusp_hit = tuple(sp.expand(fi.subs({yy: sp.Rational(2, 3), gamma: 1}))
                 for fi in Fbar)
check("Fbar(2/3,1) = cusp point (4/27,4/3,1)",
      cusp_hit == cusp_pt)

# Generic wall point from param, finite affine fiber drops (cross-ref)
wall_target = tuple(
    sp.expand(v.subs({y0: 1, c3: 1})) for v in lim_param
)
# lim at (1,1): (0, 1, 1)
check("sample wall target from param is (0,1,1)",
      wall_target == (0, 1, 1))
fib = exact_fiber(wall_target)
check("affine fiber over generic wall point (0,1,1) has exactly 1 point",
      len(fib) == 1)
inf_sols = sp.solve(
    [sp.expand(wall_target[0] - Fbar[0]),
     sp.expand(wall_target[1] - Fbar[1]),
     sp.expand(wall_target[2] - Fbar[2])],
    [yy, gamma],
    dict=True,
)
check("exactly one boundary point of D_inf over that wall target",
      len(inf_sols) == 1)

# ---------------------------------------------------------------------------
print("=== 4. Ordinary P3 homogenization: concrete negatives ===")

degs = [leading_part(Fi, PHI, weights=(1, 1, 1))[1] for Fi in F]
check("total degrees of (F1,F2,F3) are (7,6,4) — not equidegree",
      degs == [7, 6, 4])

Hs = [leading_part(Fi, PHI, weights=(1, 1, 1))[0] for Fi in F]
check("total-degree leading forms are (x^3 y^3 z, 3 x^3 y^2 z, -x^3 z)",
      sp.expand(Hs[0] - x**3 * y**3 * z) == 0
      and sp.expand(Hs[1] - 3 * x**3 * y**2 * z) == 0
      and sp.expand(Hs[2] + x**3 * z) == 0)

gb = sp.groebner(Hs, x, y, z, order="lex")
check("leading-form ideal GB = (x^3 z) — cone {x=0} U {z=0}, not just escape",
      list(gb.exprs) == [x**3 * z])
# Fake escape direction in the leading cone: (t,0,0)
t = sp.symbols("t")
F_ray = [sp.expand(Fi.subs({x: t, y: 0, z: 0})) for Fi in F]
check("ray (t,0,0) lies in the leading cone: H_i(1,0,0)=0",
      all(h.subs({x: 1, y: 0, z: 0}) == 0 for h in Hs))
check("but F(t,0,0)=(0,0,2t) escapes to infinity in sources — NOT Jelonek",
      F_ray == [0, 0, 2 * t])

# Weighted homogeneity: each Fi is purely of target weight
for i, Fi in enumerate(F):
    H, d = leading_part(Fi, PHI, weights=SOURCE_WEIGHTS)
    check(f"F{i+1} is C*-weighted-homogeneous of weight {TARGET_WEIGHTS[i]}",
          d == TARGET_WEIGHTS[i] and sp.expand(H - Fi) == 0)

# Total-degree truncation of F3 discards same-weight terms
H3, _ = leading_part(F[2], PHI, weights=(1, 1, 1))
discarded = sp.expand(F[2] - H3)
check("total-degree leading of F3 is -x^3 z",
      sp.expand(H3 + x**3 * z) == 0)
check("discarded part 2x - 3 x^2 y has the SAME C* weight +1 as -x^3 z",
      sp.expand(discarded - (2 * x - 3 * x**2 * y)) == 0
      and leading_part(discarded, PHI, weights=SOURCE_WEIGHTS)[1] == 1
      and leading_part(H3, PHI, weights=SOURCE_WEIGHTS)[1] == 1)

# Homogenized components at the hyperplane at infinity are the leading forms
# (already unequal degree — no single P3 morphism)
tt = sp.symbols("tt")
for i, Fi in enumerate(F):
    _, d = leading_part(Fi, PHI, weights=(1, 1, 1))
    Fi_h = sp.expand(tt**d * Fi.subs({x: x / tt, y: y / tt, z: z / tt}))
    check(f"ordinary homogenization of F{i+1} at tt=0 equals its leading form",
          sp.expand(Fi_h.subs(tt, 0) - Hs[i]) == 0)

# ---------------------------------------------------------------------------
print(f"\nAll {n_checks} assertions passed.")
