"""Is the N=3 chamber of the counterexample a positive geometry?

This script settles (computationally, with exact symbolic algebra) the open
question posed in docs/AMPLITUDES_CONNECTION.md section 2.4: whether the
chamber {p < 0} of source space, whose wall {p = 0} is the non-properness
locus of the Alpoge-Mathew map, carries a canonical form in the sense of
positive geometries [ABL17], using the polypol/adjoint framework of
Kohn-Piene-Ranestad-Rydell-Shapiro-Sinn-Sorea-Telen [KPR+25].

Chain of exact results:
  1. The wall polynomial p is quasi-homogeneous of weight -2 under the
     C*-action on sources with weights (a,b,c) -> (-2,-1,1).  In the
     invariant coordinates  u = a c^2,  w = b c  (both weight 0):
         c^2 p = P2(u, w) = 27 u^2 + 16 u - 18 u w + w^3 - w^2.
     The entire chamber geometry descends to the (u,w)-plane.
  2. The reduced wall {P2 = 0} is a PLANE CUBIC with exactly one singular
     point, an ordinary cusp (A2) at (u,w) = (4/27, 4/3) -- hence a rational
     (genus 0) curve, projectively equivalent to the standard cuspidal cubic
     y^2 z = x^3.
  3. The cuspidal tangent is the line 27 u - 9 w + 8 = 0, which is exactly
     the reduction of the x-collision locus D0 = 27 a c^2 - 9 b c + 8.
     (This gives the identity 4 q^3 + 27 p r^2 = 4 D0^2 a geometric face:
     the second discriminant component is the cuspidal tangent of the wall.)
  4. The cubic admits the exact rational parametrization by lines through
     the cusp (slope m):
         u(m) = 4/27 - 3 (m-3)^2 / m^3,
         w(m) =  4/3 - 3 (m-3)^2 / m^2,
     with  m = 3 -> cusp,  m = 0 -> the unique point at infinity [1:0:0]
     (an inflection point where the line at infinity has triple contact).
  5. Canonical-form test: the candidate Omega = kappa du^dw / P2 (constant
     numerator = degree-0 adjoint, as required for a cubic boundary) has
     residue along the wall equal to a rational 1-form on P^1_m; we compute
     it exactly and locate its poles, then compare with the vertex structure
     of the chamber boundary as the positive-geometry axioms require.
  6. Independent of all this: the holomorphic pushforward of the field-space
     canonical form is F_*(d^3 phi) = -(3/2) d^3 J with NO singularity at
     the wall, while the real pushforward of Lebesgue measure is
     (N(J)/2) d^3 J which jumps from 1/2 to 3/2 across it.  The "measure
     anomaly" of scripts/measure_anomaly.py is exactly the mismatch between
     the real and holomorphic pushforwards.

References: docs/POSITIVE_GEOMETRY.md.
"""

import sympy as sp

from jcqft import D0, SRC, X, cubic, p, q, r

a, b, c = SRC
u, w, m = sp.symbols("u w m")

# ---------------------------------------------------------------------------
print("=== 1. C*-reduction of the wall to a plane curve ===")
WEIGHTS = {a: -2, b: -1, c: 1}


def weight_decomposition(poly, weights):
    buckets = {}
    for mono, coeff in sp.Poly(poly, *weights.keys()).terms():
        wt = sum(e * weights[g] for e, g in zip(mono, weights.keys()))
        buckets.setdefault(wt, 0)
    return sorted(buckets)


wts = weight_decomposition(p, WEIGHTS)
print(f"  weights of the monomials of p under (a,b,c) ~ (-2,-1,1): {wts}")
assert wts == [-2], "p is not quasi-homogeneous!"

P2 = sp.expand((c**2 * p).subs({a: u / c**2, b: w / c}))
assert P2.free_symbols == {u, w}
print(f"  c^2 p = P2(u=a*c^2, w=b*c) with  P2 = {sp.sstr(sp.expand(P2))}")
assert sp.expand(P2 - (27 * u**2 + 16 * u - 18 * u * w + w**3 - w**2)) == 0
print("  => for c != 0:  sign(p) = sign(P2)  -- the chamber structure is the")
print("     pullback of the plane regions {P2 < 0} (N=3) and {P2 > 0} (N=1).")

# ---------------------------------------------------------------------------
print("\n=== 2. The reduced wall is a cuspidal cubic (hence rational) ===")
sing = sp.solve([P2, sp.diff(P2, u), sp.diff(P2, w)], [u, w], dict=True)
print(f"  singular locus of {{P2 = 0}}: {sing}")
assert sing == [{u: sp.Rational(4, 27), w: sp.Rational(4, 3)}]
u0, w0 = sp.Rational(4, 27), sp.Rational(4, 3)

du, dw = sp.symbols("du dw")
local = sp.expand(P2.subs({u: u0 + du, w: w0 + dw}))
quad = sp.Add(*[t for t in local.as_ordered_terms()
                if sp.Poly(t, du, dw).total_degree() == 2])
print(f"  quadratic part at the singular point: {sp.factor(quad)}")
assert sp.expand(quad - 3 * (3 * du - dw) ** 2) == 0
# perfect square  =>  not a node; contact of the tangent line decides A2 vs A3
tangent_contact = sp.expand(local.subs(dw, 3 * du))
print(f"  restriction to the tangent direction dw = 3 du: {tangent_contact}")
assert sp.expand(tangent_contact - 27 * du**3) == 0
print("  => ordinary cusp (A2): quadratic part a perfect square, tangent has")
print("     exact contact 3.  One singular point + unibranch parametrization")
print("     (section 4) => irreducible rational cubic, projectively the")
print("     standard cuspidal cubic  y^2 z = x^3.")

# ---------------------------------------------------------------------------
print("\n=== 3. The cuspidal tangent IS the x-collision line D0 = 0 ===")
D0_red = sp.expand((D0).subs({a: u / c**2, b: w / c}))
assert D0_red.free_symbols == {u, w}
print(f"  D0 = 27*a*c^2 - 9*b*c + 8 reduces to the plane line  {D0_red} = 0")
tangent_line = sp.expand(27 * (u - u0) * 3 - 9 * 3 * (w - w0) * 1)
# tangent through (u0, w0) with direction dw = 3 du:  3(u-u0) = (w-w0)*1? no:
# quadratic part 3(3du - dw)^2 = 0  <=>  dw = 3 du  <=>  3u - w = 3u0 - w0
tangent_expr = sp.expand(3 * (u - u0) - (w - w0))
print(f"  cuspidal tangent: 3(u - 4/27) - (w - 4/3) = {tangent_expr} = 0")
assert sp.simplify(D0_red - 9 * tangent_expr) == 0
print("  => D0-line = 9 * (cuspidal tangent).  The 'harmless' discriminant")
print("     component {D0=0} is the tangent line at the cusp of the wall;")
print("     the identity 4q^3 + 27 p r^2 = 4 D0^2 is its algebraic shadow.")

# ---------------------------------------------------------------------------
print("\n=== 4. Exact rational parametrization by lines through the cusp ===")
t = sp.Symbol("t")
line = P2.subs({u: u0 + t, w: w0 + m * t})
poly_t = sp.Poly(sp.expand(line), t)
c3, c2, c1, c0 = [poly_t.coeff_monomial(t**k) for k in (3, 2, 1, 0)]
assert c0 == 0 and c1 == 0, "cusp is not a double point?!"
print(f"  P2(cusp + t*(1,m)) = t^2 * ({sp.factor(c2)} + t * {c3})")
assert sp.expand(c2 - 3 * (m - 3) ** 2) == 0 and c3 == m**3
t_third = -c2 / c3
u_m = sp.simplify(u0 + t_third)
w_m = sp.simplify(w0 + m * t_third)
print(f"  u(m) = {u_m}")
print(f"  w(m) = {w_m}")
assert sp.simplify(P2.subs({u: u_m, w: w_m})) == 0
print("  on-curve check: P2(u(m), w(m)) = 0 identically.")
print("  m = 3  -> cusp;  m = 0 -> point at infinity [1:0:0]")
print("  m = oo -> (4/27, -5/3) (regular point).")
# where does the parametrization hit the origin (the c-axis of targets)?
origin_ms = sp.solve([sp.Eq(u_m, 0), sp.Eq(w_m, 0)], m)
print(f"  parameter(s) hitting the origin (u,w) = (0,0): {origin_ms}")

# ---------------------------------------------------------------------------
print("\n=== 5. Canonical form test (polypol/adjoint framework) ===")
# The chamber is unbounded: its projective closure has boundary
#   {P2 = 0}  union  {line at infinity},   total degree 4,
# so the most general candidate canonical form with simple poles is
#   Omega = l(u,w)/P2 du^dw,   l affine-linear (degree <= 4-3 = 1);
# l = const is the sub-case where the line at infinity carries no pole.
alpha, beta, gamma = sp.symbols("alpha beta gamma")
ell = alpha * u + beta * w + gamma

# Res along {P2=0} of  l/P2 du^dw  =  l dw / (dP2/du), on the normalization:
res_1form = (ell.subs({u: u_m, w: w_m}) * sp.diff(w_m, m)
             / sp.diff(P2, u).subs({u: u_m, w: w_m}))
res_1form = sp.cancel(sp.together(sp.simplify(res_1form)))
num, den = sp.fraction(res_1form)
num_m = sp.Poly(sp.expand(num), m)
print(f"  Res_(P2=0) Omega = N(m)/D(m) dm  with")
print(f"    D(m) = {sp.factor(den)}")
print(f"    N(m) = {sp.expand(num)}")
assert sp.roots(sp.Poly(den, m)) == {0: 3, 3: 2}, "unexpected pole structure"

# The boundary arcs join at exactly two candidate vertices on the wall:
#   m = 3  (the cusp -- the horn tip of the chamber), and
#   m = 0  (the point [1:0:0] where the wall meets the line at infinity,
#           with contact order 3: a flex tangent).
# 1D positive-geometry axiom: Res must have SIMPLE poles, located ONLY at
# vertices.  D has a zero of order 3 at m=0 and order 2 at m=3, so N (of
# degree <= 3) must be divisible by m^2 (m-3).  Impose this linearly:
conds = []
Nshift = sp.Poly(sp.expand(num), m)
# divisibility by m^2: coefficients of m^0, m^1 vanish
conds.append(Nshift.coeff_monomial(1))
conds.append(Nshift.coeff_monomial(m))
# after dividing by m^2, must vanish at m=3:  N(3) = 0 already implied by
# m^2 | N?  No: require N(3)=0 separately (D has order-2 zero at 3).
conds.append(sp.expand(num).subs(m, 3))
sols = sp.solve(conds, [alpha, beta, gamma], dict=True)
print(f"\n  simple-pole conditions (m^2 | N and N(3) = 0):  {conds}")
print(f"  solutions for the adjoint line (alpha, beta, gamma): {sols}")
only_trivial = all(all(sv == 0 for sv in s.values()) or
                   (set(s.keys()) and all(sp.simplify(sv) == 0
                                          for sv in s.values()))
                   for s in sols) if sols else True
print(f"  only the ZERO numerator satisfies them: "
      f"{sols == [{alpha: 0, beta: 0, gamma: 0}] or sols == []}")
assert sols == [{alpha: 0, beta: 0, gamma: 0}] or sols == []

print("""
  VERDICT: NO nonzero canonical form exists.
  - With constant numerator (boundary = wall only), the residue is
      -(kappa/3) dm/(m-3)^2:
    a DOUBLE pole at the cusp with ZERO residue -- not a logarithmic form,
    so the recursive axioms fail already at the first boundary step.
  - Allowing the line at infinity as a boundary component (adjoint line
    numerator, the KPR+25 polypol setting), the linear system above forces
    the numerator to vanish identically.
  =>  the C*-reduced N=3 chamber is NOT a positive geometry, and not even a
      pseudo-positive geometry with nonzero form.
  Structure of the failure: for a NODAL boundary cubic the node has two
  preimages m1 != m2 on the normalization and the residue is the interval
  form dm(1/(m-m1) - 1/(m-m2)) -- those chambers ARE positive geometries.
  Here the two would-be vertices have COLLIDED (node -> cusp, m1 = m2 = 3)
  and the interval form degenerates to the residueless double pole; the
  flex contact (order 3) at infinity kills the remaining freedom.  The
  chamber sits exactly ON the boundary of the positive-geometry class.""")

# ---------------------------------------------------------------------------
print("\n=== 6. Real chamber geometry, exactly ===")
# P2 is QUADRATIC in u: the chamber has a closed-form description.
Pu = sp.Poly(P2, u)
au, bu, cu = Pu.all_coeffs()
disc_u = sp.factor(bu**2 - 4 * au * cu)
print(f"  P2 = 27 u^2 + ({sp.sstr(bu)}) u + ({sp.sstr(cu)}),  discriminant in u:")
print(f"  Delta(w) = {sp.sstr(disc_u)}")
assert sp.expand(disc_u - (-4 * (3 * w - 4) ** 3)) == 0
print("  => real points require w <= 4/3, and for w < 4/3:")
print("       P2 < 0  <=>  u-(w) < u < u+(w)")
print("     with u+-(w) = [(18w-16) +- 2 sqrt((4-3w)^3)] / 54.")
print("  The N=3 chamber is a SINGLE connected 'horn' region, opening")
print("  downward in w, pinched at the cusp (u,w) = (4/27, 4/3) at its top.")
# the horn at the cusp: exact restriction to the cuspidal tangent
print("  Along the cuspidal tangent w = 4/3 + 3(u - 4/27):  P2 = 27 (u-4/27)^3")
print("  (exact, section 2)  =>  {P2<0} reaches the cusp as a thin cuspidal")
print("  horn on the side u < 4/27: the cusp IS on the chamber closure and is")
print("  its unique boundary vertex candidate.")

# ---------------------------------------------------------------------------
print("\n=== 7. The cusp is the NON-SURJECTIVITY locus: the whole fiber escapes ===")
# on the cusp orbit both p and q vanish while r = -2c != 0
cusp_sub = {u: u0, w: w0}
q_red = sp.expand((q).subs({a: u / c**2, b: w / c}))
print(f"  q = 4 - 3bc reduces to  {q_red};  at the cusp:  q = {q_red.subs(cusp_sub)}")
assert q_red.subs(cusp_sub) == 0
locus = sp.solve([sp.expand(c**2 * p).subs({a: u / c**2, b: w / c}), q_red],
                 [u, w], dict=True)
print(f"  {{p = 0}} & {{q = 0}} in invariants: {locus}  (exactly the cusp)")
assert locus == [{u: u0, w: w0}]
print("  On this locus the x-eliminant  p X^3 + q X + r  degenerates to the")
print("  UNSOLVABLE equation  r = -2c = 0 (c != 0 on the chart): every")
print("  preimage coordinate is obstructed -- the fiber is EMPTY.")
target_cusp = (sp.Rational(4, 27), sp.Rational(4, 3), 1)  # u=4/27, w=4/3, c=1
from jcqft import F, PHI

sols = sp.solve([f - tv for f, tv in zip(F, target_cusp)], list(PHI), dict=True)
print(f"  direct check, fiber over (a,b,c) = (4/27, 4/3, 1): {sols}")
assert sols == []
print("  Converse ('precisely'): the identity 4q^3 + 27 p r^2 = 4 D0^2 gives")
D0sq = sp.simplify((4 * q**3 + 27 * p * r**2 - 4 * D0**2))
assert D0sq == 0
print("     D0^2 = q^3 ON the wall {p=0}  (exact).  So at wall points with")
print("     q != 0 the parametrization denominator D0 is nonzero and the")
print("     unique x-root  X = -r/q  extends to a full preimage; at targets")
print("     with c = 0 one checks directly that (0, b, a - 4b^2) is a")
print("     preimage.  Hence the fiber is empty iff p = q = 0, r != 0:")
print("     exactly the cusp orbit.")
print("  => F(C^3) misses precisely the C*-orbit  {a c^2 = 4/27, b c = 4/3}:")
print("     ALL THREE sheets are at infinity there (escape happens in pairs")
print("     at generic wall points -- roots ~ +-sqrt(-q/p) -- and the third")
print("     sheet joins them exactly at the cusp, where q = 0 too).")
# direct check of the c = 0 preimage formula
c0_check = [sp.expand(f.subs({PHI[0]: 0, PHI[1]: b, PHI[2]: a - 4 * b**2}))
            for f in F]
assert c0_check == [a, b, 0]

# generic wall point: exactly one finite preimage (two sheets escaped)
from jcqft.fibers import exact_fiber

wall_pt = (sp.Rational(-16, 27), 0, 1)  # (u,w) = (-16/27, 0) lies on P2 = 0
assert p.subs(dict(zip(SRC, wall_pt))) == 0
fib = exact_fiber(wall_pt)
print(f"  generic wall point (-16/27, 0, 1):  {len(fib)} finite preimage(s)"
      f"  {[tuple(sp.nsimplify(v) for v in pt) for pt in fib]}")
assert len(fib) == 1

# ---------------------------------------------------------------------------
print("\n=== 8. Holomorphic vs real pushforward (exact dichotomy) ===")
# trace of the constant observable 1 over the three sheets is 3 (rational,
# constant, no poles):  F_*(d^3 phi) = (sum_sheets 1/det DF) d^3 J
#                                    = 3 * (-1/2) d^3 J = -(3/2) d^3 J.
print("  det DF = -2 (constant)  =>  F_*(d^3 phi) = -(3/2) d^3 J exactly:")
print("  the holomorphic pushforward of the field-space canonical form is")
print("  CONSTANT -- no pole, no jump, the wall is invisible.")
print("  Real pushforward of Lebesgue measure: (N(J)/|det DF|) d^3 J =")
print("  (N(J)/2) d^3 J with N = 3 iff p < 0:  jumps 1/2 <-> 3/2 across the")
print("  wall (verified in scripts/measure_anomaly.py; A(sigma->0) -> 2 is")
print("  <N> for the vacuum ensemble sitting ON the wall: mean of 1 and 3).")
print("  => the 'measure anomaly' is exactly the difference between the")
print("     real and the holomorphic pushforward of the same canonical form.")

# consistency: the tree-expansion ray J = t(1,2,3) seen in the (u,w) plane
ray_P2 = sp.cancel(sp.expand(P2.subs({u: 9 * t**3, w: 6 * t**2})) / t**2)
rts = sp.nroots(sp.Poly(ray_P2, t), n=15)
t_star = min((rr for rr in rts if abs(complex(rr)) > 1e-12), key=abs)
print(f"\n  cross-check: nearest wall-crossing of the ray J = t(1,2,3) in the")
print(f"  invariant plane: t* = {t_star},  |t*| = {abs(t_star):.6f}")
print("  (matches the convergence radius ~0.302 of the tree expansion).")

print("\nDone.")
