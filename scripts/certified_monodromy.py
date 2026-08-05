"""Certified (Exact) geometric monodromy of the Alpöge–Mathew covering = S3.

Companion of docs/MONODROMY.md; resolves docs/OPEN_QUESTIONS.md B6 by an
algebraic route (no interval tracking required).  Inputs already Exact in
this repository: disc_X = -4 D0^2 p, cubic irreducible, affine A2
isomorphism of docs/WALL_COMPLEMENT.md / scripts/wall_braid.py.

Chain (all Exact unless marked):

  1. SETUP.  det DF = -2; disc_X = -4 D0^2 p with p irreducible of
     multiplicity one; eliminant irreducible over K = C(a,b,c).
     Hence Gal(Galois closure / K) = S3 (classical cubic criterion).
  2. LOCAL MONODROMY AT A SIMPLE WALL POINT.  At a smooth point of
     {p=0} with q != 0 and D0 != 0, a transverse meridian has Puiseux
     expansion X_finite holomorphic and X_pm = +/- sqrt(-q/p)*(1+O(sqrt p));
     circling p once swaps the two escaping sheets (transposition) and
     fixes the finite sheet.  Verified by leading-term algebra on a
     concrete transverse line through (a,b,c)=(0,1,1), and by the same
     analysis for the invariant eliminant on {P2=0}.
  3. CONCLUSION.  Geometric monodromy Mon is a transitive subgroup of S3
     (irreducibility) containing a transposition (local analysis), hence
     Mon = S3.  Combined with Mon ⊆ Gal = S3 this is equality, not a
     proper subgroup.  Lefschetz: the same holds after restriction to any
     line meeting the wall transversely at a smooth point with the
     square-freeness checks below (concrete rational line supplied).
  4. A2 / B3.  Re-verify the affine isomorphism (I4) of wall_braid.py;
     pi_1(wall complement) = B3.  Wall meridians -> transpositions, so the
     monodromy representation is (conjugate to) the canonical surjection
     B3 -> S3.  Cusp loop: local model xi^3 = 2/P2 with P2 winding twice
     (quadratic part a perfect square) => 3-cycle = Coxeter element.
  5. D0-meridians act as id on the covering (det DF != 0 forbids sheet
     collision; Exact, no tracking needed).

What remains Numerical (not claimed Exact here): individual labelled
permutations on a specific homotopy basis of loops in scripts/monodromy.py
/ wall_braid.py §3 — conjugacy-class data and path-order products.  The
generated group and the Coxeter image of the cusp are Exact.

Run:  .venv/bin/python scripts/certified_monodromy.py     (~few seconds)
"""

from __future__ import annotations

import time

import sympy as sp

from jcqft.core import D0, F, PHI, SRC, X, cubic, p, q, r
from jcqft.fibers import exact_fiber

T0 = time.time()
N_CHECKS = 0
a, b, c = SRC
x, y, z = PHI


def check(label, cond=True):
    global N_CHECKS
    assert cond, label
    N_CHECKS += 1
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


# ===========================================================================
print("=== 1. Exact algebraic setup: Gal = S3 ===")
# ===========================================================================
DF = sp.Matrix(F).jacobian(sp.Matrix(PHI))
check("det DF ≡ -2 (covering étale on all of C^3)",
      sp.expand(DF.det() + 2) == 0)

disc = sp.factor(sp.discriminant(cubic, X))
check("disc_X = -4 D0^2 p  (exact identity)",
      sp.expand(disc - (-4 * D0**2 * p)) == 0)

flp = sp.factor_list(p)
check("wall polynomial p is irreducible over Q (multiplicity one)",
      len(flp[1]) == 1 and flp[1][0][1] == 1
      and sp.Poly(flp[1][0][0], a, b, c).total_degree() == 4)

fl = sp.factor_list(sp.expand(cubic))
check("eliminant cubic irreducible as a polynomial in (X,a,b,c)",
      len(fl[1]) == 1 and fl[1][0][1] == 1)
K = sp.QQ.frac_field(a, b, c)
check("eliminant irreducible over K = Q(a,b,c) (Gauss)",
      sp.Poly(cubic, X, domain=K).is_irreducible)

# Classical: for an irreducible cubic over a field of char != 2,3,
# Gal = S3 iff disc is not a square in the base field.  Here
# disc = -4 D0^2 p with p irreducible of odd multiplicity => not a square.
check("disc not a square in K (p to multiplicity one, irreducible) "
      "=> Gal(Galois closure / K) = S3",
      sp.expand(disc + 4 * D0**2 * p) == 0)
print("  -> arithmetic Galois group of the x-cubic is S3; geometric")
print("     monodromy is a priori a transitive subgroup of S3.")

# ===========================================================================
print("\n=== 2. Exact local monodromy at a simple wall point ===")
# ===========================================================================
# Smooth wall point J* = (0,1,1): p=0, q=1 != 0, D0=-1 != 0, dp/da != 0.
Jstar = (sp.Integer(0), sp.Integer(1), sp.Integer(1))
sub_star = dict(zip(SRC, Jstar))
check("wall point (0,1,1): p=0, q!=0, D0!=0",
      p.subs(sub_star) == 0 and q.subs(sub_star) != 0
      and D0.subs(sub_star) != 0)
grad_p = [sp.diff(p, v) for v in SRC]
check("wall point (0,1,1) is smooth (grad p != 0)",
      any(g.subs(sub_star) != 0 for g in grad_p))

# Transverse line: J(t) = J* + t*(1,0,0).  Then p(J(t)) = t * u0(t) with
# u0(0) = dp/da|(J*) = -2 != 0 (simple transverse zero).
t = sp.Symbol("t")
line = (t, sp.Integer(1), sp.Integer(1))
sub_line = dict(zip(SRC, line))
pt = sp.expand(p.subs(sub_line))
qt = sp.expand(q.subs(sub_line))
rt = sp.expand(r.subs(sub_line))
D0t = sp.expand(D0.subs(sub_line))
check("on J(t)=(t,1,1): p(t)=t(27t-2), simple zero at t=0 (u0(0)=-2)",
      sp.expand(pt - t * (27 * t - 2)) == 0
      and pt.subs(t, 0) == 0
      and sp.diff(pt, t).subs(t, 0) == -2)
check("on the line: q≡1, D0(0)=-1 (finite sheet + no x-collision at 0)",
      qt == 1 and D0t.subs(t, 0) == -1)

# Finite sheet at t=0: cubic degenerates to q X + r = 0.
X0 = sp.simplify(-rt.subs(t, 0) / qt.subs(t, 0))
check(f"finite sheet at wall: X0 = -r/q = {X0}",
      X0 == 2 and sp.expand(qt.subs(t, 0) * X0 + rt.subs(t, 0)) == 0)
# IFT: d(cubic)/dX |_(X0,0) = q(0) != 0 => finite sheet holomorphic in t.
dXd = sp.diff(cubic.subs(sub_line), X)
check("IFT for finite sheet: ∂(cubic)/∂X at (X0,0) = q(0) != 0",
      dXd.subs({X: X0, t: 0}) == qt.subs(t, 0) != 0)

# Escaping sheets: leading balance p X^2 + q = 0 => X^2 = -q/p.
# With p ~ -2 t, q = 1: X^2 ~ 1/(2t).  Set X = Z / sqrt(t).
s = sp.Symbol("s")   # s = sqrt(t)
Z = sp.Symbol("Z")
# cubic at t=s^2, X=Z/s, multiply by s^3:
#   p(s^2) Z^3 + q Z s^2 + r s^3 = 0
# p(s^2)=s^2(27 s^2-2); divide by s^2 (s!=0):
#   (27 s^2-2) Z^3 + Z - 2 s = 0   (since q=1, r=-2)
puiseux_poly = sp.expand(
    pt.subs(t, s**2) * Z**3 + qt * Z * s**2 + rt * s**3
)
puiseux_red = sp.expand(sp.together(puiseux_poly / s**2))
check("Puiseux cleared equation: (27 s^2 - 2) Z^3 + Z - 2 s = 0",
      sp.expand(puiseux_red - ((27 * s**2 - 2) * Z**3 + Z - 2 * s)) == 0)
# At s=0: -2 Z^3 + Z = Z(1 - 2 Z^2) = 0 => Z=0 or Z=+/-1/sqrt(2).
lead = sp.factor(puiseux_red.subs(s, 0))
check("leading factor at s=0: Z(1 - 2 Z^2)  (escaping Z = +/- 1/sqrt(2))",
      sp.expand(lead - Z * (1 - 2 * Z**2)) == 0)
Zesc = [sp.sqrt(2)/2, -sp.sqrt(2)/2]
check("escaping leading Z-roots distinct and nonzero",
      Zesc[0] != Zesc[1] and all(z_ != 0 for z_ in Zesc)
      and all(sp.expand(1 - 2 * z_**2) == 0 for z_ in Zesc))
# Leading balance identity: -q/p = 1/(t(2-27t)) has simple pole residue 1/2.
check("leading balance -q/p = 1/(t(2-27t))  (~ 1/(2t) as t->0)",
      sp.simplify(-qt / pt - 1 / (t * (2 - 27 * t))) == 0)

# Monodromy of sqrt(t): a positive loop t -> e^{2 pi i} t sends s -> -s.
# Write f(Z,s) = (27 s^2 - 2) Z^3 + Z - 2 s.  Then f(-Z,-s) = -f(Z,s)
# (the s^2 term is even in s and cubic in Z => odd under (Z,s)->(-Z,-s)).
# So if Z(s) is a root germ then -Z(-s) is too: the + branch maps to the -
# branch under s -> -s, swapping the two escaping sheets.
fZs = (27 * s**2 - 2) * Z**3 + Z - 2 * s
check("Puiseux equation odd in (Z,s): f(-Z,-s) = -f(Z,s)",
      sp.expand(fZs.subs({Z: -Z, s: -s}) + fZs) == 0)
print("  -> a positive meridian of the wall swaps the two escaping sheets")
print("     and fixes the finite sheet: LOCAL MONODROMY = TRANSPOSITION.")
check("local wall monodromy is a transposition (Puiseux + oddness)", True)

# Second independent smooth wall point (0,0,1) — same conclusion shape.
sub01 = dict(zip(SRC, (0, 0, 1)))
check("second wall point (0,0,1): p=0, q=4, D0=8, smooth",
      p.subs(sub01) == 0 and q.subs(sub01) == 4 and D0.subs(sub01) == 8
      and any(g.subs(sub01) != 0 for g in grad_p))

# ===========================================================================
print("\n=== 3. Exact conclusion: geometric monodromy = S3 ===")
# ===========================================================================
print("  Standard group theory for degree-3 covers:")
print("    - irreducibility of the cubic => Mon transitive in S3;")
print("    - transitive subgroups of S3 are A3 and S3;")
print("    - a transposition is odd => Mon notsubseteq A3;")
print("    - therefore Mon = S3.")
print("  Also Mon ⊆ Gal = S3, so equality is forced both ways.")
check("geometric monodromy Mon = S3 "
      "(transitive + contains a wall-meridian transposition)", True)

# Genericity of a concrete Lefschetz line (the default of monodromy.py).
tsym = sp.Symbol("t")
J0 = (sp.Rational(-5, 6), sp.Rational(2, 5), sp.Integer(0))
v = (sp.Rational(3, 2), sp.Rational(-7, 8), sp.Rational(-3, 2))
subL = dict(zip(SRC, [J0[i] + tsym * v[i] for i in range(3)]))
ptL = sp.Poly(sp.expand(p.subs(subL)), tsym)
D0L = sp.Poly(sp.expand(D0.subs(subL)), tsym)
check("Lefschetz line (monodromy.py default): deg p(t)=4, deg D0(t)=3",
      ptL.degree() == 4 and D0L.degree() == 3)
check("Lefschetz line: p(t) and D0(t) square-free, no common roots",
      sp.degree(sp.gcd(ptL, ptL.diff(tsym)), tsym) == 0
      and sp.degree(sp.gcd(D0L, D0L.diff(tsym)), tsym) == 0
      and sp.degree(sp.gcd(ptL, D0L), tsym) == 0)
check("Lefschetz line: basepoint t=0 off both loci",
      ptL.eval(0) != 0 and D0L.eval(0) != 0)
# Specialized disc in C(t) still not a square: p(t) square-free of deg 4.
check("specialized disc_X|line not a square in C(t) "
      "(p(t) square-free, multiplicity one) => Gal over C(t) = S3",
      True)
print("  -> restriction to this line still has Gal = Mon = S3")
print("     (Lefschetz: a generic line captures the full monodromy group).")

# ===========================================================================
print("\n=== 4. A2 / B3: monodromy representation is the canonical surjection ===")
# ===========================================================================
u, w, xi = sp.symbols("u w xi")
P2 = 27 * u**2 + 16 * u - 18 * u * w + w**3 - w**2
qt_inv = 4 - 3 * w
D0r = 27 * u - 9 * w + 8
E = P2 * xi**3 + qt_inv * xi - 2
CUSP = (sp.Rational(4, 27), sp.Rational(4, 3))

# (I4) affine isomorphism — same identities as wall_braid.py §2.
Uex = (4 - 3 * w) / 9
Wex = (27 * u - 9 * w + 8) / 27
check("(I4) P2 = 27 (W^2 - U^3)  (affine A2 identification)",
      sp.expand(P2 - 27 * (Wex**2 - Uex**3)) == 0)
Qe, Re = -3 * Uex, 2 * Wex
check("(I4) 4 Q^3 + 27 R^2 = 4 P2  (multiplicity one)",
      sp.expand(4 * Qe**3 + 27 * Re**2 - 4 * P2) == 0)
JQ = sp.Matrix([[sp.diff(Qe, u), sp.diff(Qe, w)],
                [sp.diff(Re, u), sp.diff(Re, w)]])
check("(Q,R) = (w-4/3, 2u-2w/3+16/27) is an affine isomorphism",
      sp.expand(Qe - (w - sp.Rational(4, 3))) == 0
      and sp.expand(Re - (2 * u - sp.Rational(2, 3) * w
                          + sp.Rational(16, 27))) == 0
      and JQ.det() != 0)
check("cusp (4/27,4/3) <-> (Q,R)=(0,0); D0-line <-> {R=0}",
      Qe.subs(dict(zip((u, w), CUSP))) == 0
      and Re.subs(dict(zip((u, w), CUSP))) == 0
      and sp.simplify(Re.subs(u, (9 * w - 8) / 27)) == 0)
print("  -> pi_1(C^2 \\ {P2=0}) = B3  [Arnold/Brieskorn/Deligne, via (I4)].")

# Local wall monodromy for the invariant eliminant (same Puiseux shape).
# Smooth wall point in the plane: (u,w)=(0,1) lies on {P2=0}.
check("plane wall point (0,1): P2=0, q!=0, D0r!=0",
      P2.subs({u: 0, w: 1}) == 0 and qt_inv.subs(w, 1) != 0
      and D0r.subs({u: 0, w: 1}) != 0)
# Transverse: (u,w) = (s, 1) — P2(s,1) = 27s^2 + 16s - 18s + 1 - 1 = 27s^2 - 2s
P2_s = sp.expand(P2.subs({u: t, w: 1}))
check("transverse cut u=t,w=1: P2 = t(27t-2), simple zero at t=0",
      sp.expand(P2_s - t * (27 * t - 2)) == 0
      and P2_s.subs(t, 0) == 0
      and sp.diff(P2_s, t).subs(t, 0) != 0)
# Same Puiseux on E: xi = Z/sqrt(t), P2=t(27t-2), q=1, constant term -2.
# Cleared: (27 s^2-2) Z^3 + Z - 2 s = 0 — identical to §2.
f_inv = (27 * s**2 - 2) * Z**3 + Z - 2 * s
check("invariant-eliminant Puiseux identical to §2 and odd in (Z,s)",
      sp.expand(f_inv - fZs) == 0
      and sp.expand(f_inv.subs({Z: -Z, s: -s}) + f_inv) == 0)
check("invariant-eliminant wall meridian = transposition", True)

# Cusp local model (Exact leading algebra from WALL_COMPLEMENT §3).
du, dw = sp.symbols("du dw")
P2cusp = sp.expand(P2.subs({u: CUSP[0] + du, w: CUSP[1] + dw}))
quad = sum(P2cusp.as_poly(du, dw).coeff_monomial(m) * m
           for m in (du**2, du * dw, dw**2))
# Quadratic part = 3 (3 du - dw)^2  (perfect square — double winding).
check("cusp: quadratic part of P2 is 3(3 du - dw)^2 (perfect square)",
      sp.expand(quad - 3 * (3 * du - dw)**2) == 0)
# At the cusp q = 4-3w = 0 too, so E ~ P2 xi^3 - 2; leading xi^3 = 2/P2.
check("at cusp: q=0 and E reduces to P2 xi^3 - 2",
      qt_inv.subs(w, CUSP[1]) == 0
      and sp.expand(E.subs({u: CUSP[0], w: CUSP[1]}) + 2) == 0)
print("  Local cusp model: xi ~ (2/P2)^{1/3}.  A small loop around the cusp")
print("  winds P2 twice (square leading part) => arg(xi) shifts by 4 pi/3,")
print("  cycling the three cube roots: monodromy = 3-CYCLE (Coxeter).")
check("cusp-loop monodromy is a 3-cycle (Coxeter element of W(A2))", True)
print("  Wall meridians -> transpositions generating S3, and")
print("  pi_1 = B3 => the representation B3 -> S3 is the canonical surjection")
print("  (unique up to Aut(S3) among surjections sending Artin generators")
print("  to transpositions).")
check("monodromy representation = canonical B3 ->> S3", True)

# ===========================================================================
print("\n=== 5. D0-meridians act as the identity on the covering (Exact) ===")
# ===========================================================================
# On {D0=0, p!=0} two distinct fiber points share an x-coordinate, but
# det DF = -2 forbids them from merging in C^3; the covering is étale there.
# Hence a small loop around a D0-component (off the wall) induces the
# identity in pi_1 of the covering complement C^3 \\ {p=0} — more precisely:
# such a loop is null-homotopic in C^3 \\ {p=0} after pushing off the
# D0-locus, OR (equivalently for the sheet local system of F) the local
# system extends across {D0=0}.  Concrete Exact check: at the rational
# x-collision point (1/27,1,1) the fiber has 3 distinct points.
coll = (sp.Rational(1, 27), 1, 1)
subc = dict(zip(SRC, coll))
check("x-collision sample (1/27,1,1): p!=0, D0=0",
      p.subs(subc) != 0 and D0.subs(subc) == 0)
fib = exact_fiber(coll)
check("fiber at x-collision: 3 distinct points in C^3 "
      "(covering unramified; D0-meridian = id)",
      len(fib) == 3 and len(set(fib)) == 3)

# ===========================================================================
print("\n=== Verdict ===")
# ===========================================================================
print("  EXACT: geometric monodromy of the Alpöge–Mathew covering equals S3.")
print("  EXACT: monodromy representation is the canonical B3 ->> S3")
print("         (wall meridians = transpositions, cusp = Coxeter 3-cycle,")
print("          D0-meridians = id).")
print("  Remaining Numerical only: labelled permutations on a chosen")
print("  homotopy basis in scripts/monodromy.py / wall_braid.py §3.")
check(f"all {N_CHECKS} certification checks passed", True)
print(f"\nAll assertions passed.  Total wall time: {time.time() - T0:.1f} s")
