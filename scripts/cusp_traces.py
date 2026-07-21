"""A1: exact trace asymptotics at the wall and at the cusp.

Settles open question A1 of docs/OPEN_QUESTIONS.md (= POSITIVE_GEOMETRY.md
section 5 Q3): the behavior of the trace (sum-over-sheets) observables as the
source approaches the non-properness wall {p=0} and its cusp stratum, where
the fiber is empty.

Results (all exact):
  1. Newton power sums of the x-sheets: S_1 = 0 identically; S_2 = -2q/p,
     S_3 = -3r/p, S_4 = 2q^2/p^2, S_5 = 5qr/p^2, S_6 = -2q^3/p^3 + 3r^2/p^2.
     Every S_k is a polynomial in q/p and r/p; poles ONLY on the wall.
  2. GENERIC WALL POINT (p -> 0, q != 0): two sheets escape as a
     square-root pair  x ~ +-sqrt(-q/p);  S_k diverges with pole order
     floor(k/2) in p.  The odd-k leading terms CANCEL between the pair.
  3. CUSP, GENERIC APPROACH (p ~ eps^2, q ~ eps, r = O(1)): all three
     sheets escape as CUBE ROOTS  x_j ~ omega^j (-r/p)^{1/3}.  Because
     1 + omega + omega^2 = 0, the naive divergence eps^{-2k/3} of S_k is
     reduced whenever 3 does not divide k.  Exact law: the eps-order of
     S_k is  max{ a + 2b : 2a + 3b = k, a,b >= 0 }  (from S_k as a
     quasi-homogeneous polynomial in e_2 = q/p ~ eps^{-1} and
     e_3 = -r/p ~ eps^{-2}), always <= 2k/3, with equality iff 3 | k.
  4. CUSP, TANGENT APPROACH (along the cuspidal tangent = the D0-line):
     the theory is EXACTLY SOLVABLE.  With c = 1 and
     (u,w) = cusp + eps*(1,3):  P2 = 27 eps^3 exactly, q = -9 eps,
     r = -2, and the eliminant factors for every eps:
         27 eps^3 X^3 - 9 eps X - 2 = (3 eps X - 2)(3 eps X + 1)^2 .
     The fiber x-coordinates are exactly {2/(3 eps), -1/(3 eps) (double)}:
     the x-collision of the D0-locus persists all the way into the cusp,
     and all sheets escape like 1/eps.  No cube-root regime on the
     tangent: there the pair structure survives.
  5. CORRECTION to POSITIVE_GEOMETRY.md section 5 Q3 as originally posed:
     e_2 = q/p and e_3 = -r/p are NOT finite at the cusp -- they diverge
     (all traces have their poles on the whole wall, cusp included).  What
     is true and sharp is: S_1 = 0 exactly everywhere; the divergence
     RATES at the cusp are anomalously small for 3 ∤ k (cube-root
     cancellation); and the tangent approach is exactly solvable.
"""

import mpmath as mp
import numpy as np
import sympy as sp

from jcqft import SRC, p, q, r

a, b, c = SRC
eps = sp.Symbol("epsilon", positive=True)

# ---------------------------------------------------------------------------
print("=== 1. Exact power sums of the x-sheets (Newton, monic form) ===")
# monic cubic: x^3 + e2p*x + e3p with e2p = q/p, e3p = -(-r/p)... careful:
# x^3 + (q/p) x + (r/p) = 0  =>  e1 = 0, e2 = q/p, e3 = -r/p.
e2, e3 = q / p, -r / p
S = {0: sp.Integer(3), 1: sp.Integer(0)}
S[2] = -2 * e2
for k in range(3, 13):
    S[k] = sp.expand(-e2 * S[k - 2] + e3 * S[k - 3])  # e1 = 0

expected = {2: -2 * q / p, 3: -3 * r / p, 4: 2 * q**2 / p**2,
            5: 5 * q * r / p**2, 6: -2 * q**3 / p**3 + 3 * r**2 / p**2}
for k, ex in expected.items():
    assert sp.simplify(S[k] - ex) == 0
    print(f"  S_{k} = {sp.sstr(sp.factor(ex))}")
print("  S_1 = 0 identically: the sheet-sum of the field's x-component")
print("  vanishes EXACTLY, at every source, to all orders.")

# numeric cross-check of the Newton recursion at a random rational point
mp.mp.dps = 30
to_mp = lambda v: mp.mpf(int(v.p)) / mp.mpf(int(v.q))
pt = {a: sp.Rational(3, 7), b: sp.Rational(-2, 5), c: sp.Rational(1, 3)}
coeffs = [sp.Rational(v.subs(pt)) for v in (p, sp.Integer(0), q, r)]
roots = mp.polyroots([to_mp(v) for v in coeffs], maxsteps=200)
for k in (2, 3, 4, 5, 6):
    num = sum(rt**k for rt in roots)
    sym = to_mp(sp.Rational(S[k].subs(pt)))
    assert abs(num - sym) < mp.mpf("1e-18") * (1 + abs(sym))
print("  (Newton recursion cross-checked numerically at a rational point.)")

# generic power sums as functions of abstract (pp, qq, rr), for line restrictions
pp, qq, rr = sp.symbols("pp qq rr")
Sgen = {1: sp.Integer(0), 2: -2 * qq / pp, 3: 3 * (-rr / pp)}
for k in range(4, 13):
    Sgen[k] = sp.expand(-(qq / pp) * Sgen[k - 2] + (-rr / pp) * Sgen[k - 3])


def pole_order_in_eps(expr):
    """leading pole order of a rational function of eps as eps -> 0."""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    lo = lambda e: 0 if not e.has(eps) else min(m[0] for m in sp.Poly(e, eps).monoms())
    return lo(den) - lo(num)


# ---------------------------------------------------------------------------
print("\n=== 2. Generic wall point: square-root pair escape, floor(k/2) law ===")
# transverse line through the wall point (-16/27, 0, 1) (p = 0, q = 4 there)
line = {a: sp.Rational(-16, 27) + eps, b: sp.Integer(0), c: sp.Integer(1)}
p_line, q_line, r_line = (sp.expand(v.subs(line)) for v in (p, q, r))
print(f"  along J = (-16/27,0,1) + eps*(1,0,0):  p = {sp.sstr(p_line)},  "
      f"q = {q_line},  r = {r_line}")
for k in range(2, 9):
    Sk = Sgen[k].subs({pp: p_line, qq: q_line, rr: r_line})
    lead = pole_order_in_eps(Sk)
    print(f"  S_{k}: pole order in eps = {lead}   (floor(k/2) = {k // 2})")
    assert lead == k // 2

# ---------------------------------------------------------------------------
print("\n=== 3. Cusp, generic approach: cube-root escape and omega-cancellation ===")
u, w = sp.symbols("u w")
P2 = 27 * u**2 + 16 * u - 18 * u * w + w**3 - w**2
u0, w0 = sp.Rational(4, 27), sp.Rational(4, 3)
# generic direction (du, dw) = (1, 1); chart c = 1 so p = P2, q = 4 - 3w, r = -2
p_c = sp.expand(P2.subs({u: u0 + eps, w: w0 + eps}))
q_c = sp.expand((4 - 3 * w).subs(w, w0 + eps))
r_c = sp.Integer(-2)
print(f"  p(eps) = {sp.sstr(p_c)}   (leading 3*(3*du-dw)^2 eps^2 = 12 eps^2)")
print(f"  q(eps) = {q_c},  r = {r_c}")
assert sp.Poly(p_c, eps).all_coeffs()[-1] == 0  # p(0) = 0: on the wall

print("  eps-order of S_k (actual vs naive 2k/3 vs law max{a+2b : 2a+3b=k}):")
for k in range(2, 10):
    Sk = Sgen[k].subs({pp: p_c, qq: q_c, rr: r_c})
    actual = pole_order_in_eps(Sk)
    law = max(aa + 2 * bb for bb in range(k // 3 + 1)
              for aa in ((k - 3 * bb) // 2,) if 2 * aa + 3 * bb == k)
    naive = sp.Rational(2 * k, 3)
    tag = "(cancellation)" if actual < naive else "(no cancellation, 3|k)"
    print(f"    k={k}:  actual {actual}   naive {naive}   law {law}   {tag}")
    assert actual == law and (actual < naive) == (k % 3 != 0)

# numeric confirmation that the roots really are ~ omega^j (-r/p)^(1/3):
# solve the RESCALED cubic in Y = X/s, s = (2/p)^(1/3), which is
# well-conditioned:  2 Y^3 + (q s) Y - 2 = 0.
ee = 1e-8
pv = mp.mpf(ee) ** 3 + 12 * mp.mpf(ee) ** 2
qv = -3 * mp.mpf(ee)
s = (2 / pv) ** (mp.mpf(1) / 3)
rts_Y = mp.polyroots([2, 0, qv * s, -2], maxsteps=200)
rts = [s * y for y in rts_Y]
scale = abs(rts[0])
phases = sorted(round(float(mp.arg(rt / rts[0]) / mp.pi) * 3) / 3 for rt in rts)
print(f"  numeric at eps=1e-8: |x_j|/|x_0| = "
      f"{[round(float(abs(rt) / scale), 4) for rt in rts]}, relative phases "
      f"{phases} pi  -> cube roots of unity")
assert max(abs(abs(rt) / scale - 1) for rt in rts) < 1e-2
assert abs(sum(rts)) < 1e-4 * scale  # S_1 = 0 numerically too

# ---------------------------------------------------------------------------
print("\n=== 4. Cusp, tangent approach: EXACT solvability on the D0-line ===")
X = sp.Symbol("X")
elim_tan = 27 * eps**3 * X**3 - 9 * eps * X - 2
fact = sp.factor(elim_tan)
print(f"  eliminant along (u,w) = cusp + eps*(1,3), c=1:")
print(f"    27 eps^3 X^3 - 9 eps X - 2 = {sp.sstr(fact)}")
assert sp.expand(elim_tan - (3 * eps * X - 2) * (3 * eps * X + 1) ** 2) == 0
print("  fiber x-coordinates EXACTLY: { 2/(3 eps),  -1/(3 eps) (double) }")
print("  -> the x-collision structure of the D0-locus persists into the cusp;")
print("     on the tangent there is NO cube-root regime: q*X balances p*X^3.")
# full 3D check: the two x-degenerate points are distinct in (y,z)
from jcqft.fibers import exact_fiber

eps_val = sp.Rational(1, 10)
target = (sp.Rational(4, 27) + eps_val, sp.Rational(4, 3) + 3 * eps_val, 1)
fib = exact_fiber(target)
xs = sorted(sp.nsimplify(pt[0]) for pt in fib)
print(f"  exact fiber over cusp+eps*(1,3), eps=1/10: x-coords = {xs}")
assert xs == [sp.Rational(-10, 3), sp.Rational(-10, 3), sp.Rational(20, 3)]
distinct = len({tuple(pt) for pt in fib})
print(f"  number of DISTINCT points in the fiber: {distinct} "
      "(collided x, distinct (y,z): etale covering unbroken)")
assert distinct == 3

# ---------------------------------------------------------------------------
print("\n=== 5. Conclusion (corrects POSITIVE_GEOMETRY.md §5 Q3 as posed) ===")
print("""  e_2 = q/p and e_3 = -r/p DIVERGE at the cusp (as everywhere on the
  wall); the original parenthetical claiming removable behavior was wrong.
  The true exact statements:
    (i)   S_1 = 0 identically (the only trace that stays finite).
    (ii)  Pole orders: floor(k/2) in p at generic wall points (pair
          escape, odd-power cancellation).
    (iii) At the cusp the divergence rate of S_k drops below the naive
          2k/3 exactly when 3 does not divide k: the three escaping
          sheets are asymptotic cube roots of unity and their
          contributions cancel at leading order (the 'omega-cancellation'
          -- the trace remembers that the escaping configuration is a
          Z_3-symmetric triple even though no fiber point exists at the
          cusp itself).
    (iv)  Along the cuspidal tangent (the D0-line) the model is exactly
          solvable for all eps; the escaping fiber keeps the collided
          pair structure of the D0-locus.""")

print("Done.")
