"""Exact and numeric fiber computation for the counterexample map.

The lex Groebner basis of the ideal (F(phi) - J) yields:
  * the x-eliminant  p X^3 + q X + r  (re-derived and checked against core),
  * a rational parametrization of the rest of the fiber,
        y = -B/A,  z = -D/C,  A = 2*D0,  C = 8*D0,
    valid off the x-collision locus {D0 = 0}.
"""

import sympy as sp

from jcqft.core import D0, F, PHI, SRC, X, cubic, p, q, r

x, y, z = PHI
a, b, c = SRC

SYSTEM = [F[0] - a, F[1] - b, F[2] - c]

_gb = sp.groebner(SYSTEM, y, z, x, order="lex")
g_y = next(g for g in _gb.exprs if sp.degree(g, y) == 1)
g_z = next(g for g in _gb.exprs if sp.degree(g, z) == 1)
g_x = next(g for g in _gb.exprs if not g.has(y) and not g.has(z))
A, B = sp.Poly(g_y, y).all_coeffs()  # A*y + B = 0,  A independent of x
C, D = sp.Poly(g_z, z).all_coeffs()  # C*z + D = 0,  C independent of x

assert sp.expand(g_x.subs(x, X) - cubic) == 0, "eliminant mismatch with core"
assert sp.expand(A - 2 * D0) == 0 and sp.expand(C - 8 * D0) == 0

y_of_x = -B / A
z_of_x = -D / C

# numeric helpers
F_num = sp.lambdify(PHI, F, "numpy")
p_num = sp.lambdify(SRC, p, "numpy")
coef_num = sp.lambdify(SRC, (p, q, r), "numpy")  # cubic is p*X^3 + q*X + r
yz_num = sp.lambdify((x,) + SRC, (y_of_x, z_of_x), "numpy")


def exact_fiber(target):
    """All finite preimages of a rational/symbolic target, exactly.

    Uses the cubic + rational parametrization; falls back to sp.solve on the
    x-collision locus where the parametrization denominator vanishes."""
    sub = dict(zip(SRC, target))
    if D0.subs(sub) != 0:
        pts = []
        for xr in sp.Poly(g_x.subs(sub), x).all_roots():
            pts.append(tuple(sp.simplify(v) for v in
                             (xr, y_of_x.subs(sub).subs(x, xr),
                              z_of_x.subs(sub).subs(x, xr))))
        return pts
    sols = sp.solve([eq.subs(sub) for eq in SYSTEM], [x, y, z], dict=True)
    return [tuple(s[v] for v in (x, y, z)) for s in sols]
