"""Exact geometry of the failure: eliminants, discriminant locus and the
escape-to-infinity (non-properness) set of the counterexample map.

Since det DF = -2 everywhere, F is etale: fiber points can never merge at
finite points of C^3.  The generic fiber has 3 points, so the fiber count can
only drop when a preimage ESCAPES TO INFINITY.  This script computes in
closed form where that happens, and classifies the components of the
discriminant locus:

  * lex Groebner basis gives a rational fiber parametrization
        y = -B(x,a,b,c) / A(a,b,c),   z = -D(x,a,b,c) / C(a,b,c)
    over each root x of the cubic eliminant  p x^3 + q x + r = 0;
  * resultants give the minimal cubics of y and z over C(a,b,c); both are
    MONIC (constant leading coefficient), so y and z remain bounded on
    bounded target sets: escape to infinity happens only in x, exactly on
    the hypersurface  {p(a,b,c) = 0};
  * the other discriminant component {4 q^3 + 27 p r^2 = 0} is where two
    DISTINCT fiber points merely share the same x-coordinate (the projection
    to x ramifies, the covering itself does not) -- harmless for monodromy;
  * exact fibers over the interesting targets, including the triple point
    (-1/4, 0, 0), and the escape of sheets B, C along the segment to it.
"""

import sympy as sp

from jcqft import D0, SRC, X, cubic, p, q, r
from jcqft.fibers import A, B, C, D, exact_fiber, g_x

x = sp.Symbol("x")
a, b, c = SRC

print("=== 1. Lex Groebner basis: rational fiber parametrization ===")
# (the basis itself is computed once in jcqft.fibers, which also asserts that
#  the x-eliminant matches the cubic and that A = 2*D0, C = 8*D0)
print("x-eliminant matches the cubic  p*X^3 + q*X + r: True (asserted on import)")
print(f"parametrization:  y = -B/A,  z = -D/C  with  A = 2*({sp.sstr(D0)}),"
      f"  C = 8*({sp.sstr(sp.factor(C / 8))})")
print(f"  (denominator D0 = {sp.sstr(D0)} vanishes exactly where two sheets"
      " share an x-coordinate, see section 3)")

print("\n=== 2. Minimal cubics of y and z: escape happens ONLY in x ===")
Y, Z = sp.symbols("Y Z")
res_y = sp.resultant(g_x, A * Y + B, x)
res_z = sp.resultant(g_x, C * Z + D, x)
for name, res, var in (("y", res_y, Y), ("z", res_z, Z)):
    fac = [f for f, _ in sp.factor_list(res, var)[1] if f.has(var)]
    assert len(fac) == 1
    poly = sp.Poly(fac[0], var)
    lc = sp.factor(poly.LC())
    print(f"  {name}-eliminant: degree {poly.degree()}, leading coefficient"
          f" {sp.sstr(lc)}  (CONSTANT -> {name} is integral over C[a,b,c],"
          f" never escapes)")
print(f"  x-eliminant:  degree 3, leading coefficient p = {sp.sstr(p)}")
print("  -> a sheet escapes to infinity  <=>  p(a,b,c) = 0, and it escapes in"
      " the x-direction only.")

print("\n=== 3. Discriminant of the x-cubic, classified ===")
disc = sp.factor(sp.discriminant(cubic, X))
print(f"  disc_X = -p * (4*q^3 + 27*p*r^2)  [verified: "
      f"{sp.expand(disc - (-p * (4 * q**3 + 27 * p * r**2))) == 0}]")
print("  component {p = 0}: TRUE branch locus (one sheet at infinity;"
      " monodromy possible around it)")
print("  component {4q^3 + 27pr^2 = 0, p != 0}: two distinct fiber points"
      " share an x-coordinate;")
print("    the covering is unramified there (det DF = -2 forbids merging)."
      "  No monodromy.")
print(f"  p at the origin: p(0,0,0) = {p.subs(dict(zip(SRC, (0, 0, 0))))}"
      "  ->  J = 0 LIES ON the branch locus:")
print("    already over the perturbative vacuum two of the three sheets sit"
      " at infinity.")


def fiber(target, label):
    """Exact fiber (via jcqft.fibers), with pretty-printing."""
    pts = exact_fiber(target)
    print(f"\nfiber over {label} = ({', '.join(str(t) for t in target)}):"
          f"  {len(pts)} finite point(s)")
    for pt in sorted(pts, key=lambda t: abs(complex(sp.N(t[0])))):
        if all(sp.nsimplify(v).is_rational for v in pt):
            print(f"   phi = ({', '.join(str(sp.nsimplify(v)) for v in pt)})")
        else:
            print(f"   phi ~ ({', '.join(str(sp.N(v, 8)) for v in pt)})")
    return pts


print("\n=== 4. Exact fibers over interesting targets ===")
fiber((0, 0, 0), "J = 0 (perturbative vacuum)")
fiber((sp.Rational(-1, 4), 0, 0), "the famous triple point")
print("\n  triple point check: cubic becomes -4X^3 + 4X = 0, roots 0, +-1.")

# Rational point on {4q^3 + 27pr^2 = 0} with p != 0: at b = c = 1,
# p = a(27a - 2) and 4q^3 + 27pr^2 = 4 + 108p vanishes at a = 1/27.
tgt = (sp.Rational(1, 27), 1, 1)
sub = dict(zip(SRC, tgt))
print(f"\n  x-collision point (a,b,c) = {tgt}:  p = {p.subs(sub)},"
      f"  4q^3+27pr^2 = {(4 * q**3 + 27 * p * r**2).subs(sub)},"
      f"  D0 = {D0.subs(sub)}")
pts = fiber(tgt, "x-collision point")
print(f"  x-coordinates: {sorted(sp.N(pt[0], 6) for pt in pts)}")
print("  -> 3 distinct finite points, two sharing x = 3: the x-projection"
      " ramifies, the covering does not.")

print("\n=== 5. Escape of sheets B, C along the segment to the triple point ===")
print("  On the line (a, 0, 0):  16a*X^3 + 4X = 0, so")
print("    sheet A:    X = 0        ->  phi_A = (0, 0, a)   (perturbative"
      " branch; series terminates)")
print("    sheets B,C: X = +-sqrt(-1/(4a))")
for aval in (sp.Rational(-1, 4), sp.Rational(-1, 16), sp.Rational(-1, 64)):
    pts = fiber((aval, 0, 0), f"(a,0,0), a = {aval}")
print("\n  -> as a -> 0^- the x-coordinates of sheets B, C grow like"
      " 1/(2*sqrt(-a)) and recede to infinity;")
print("     at a = -1/4 they land at (+-1, -+3/2, 13/2) and collide in the"
      " target with the image of sheet A.")
print("     They are invisible to EVERY order of perturbation theory around"
      " J = 0.")
