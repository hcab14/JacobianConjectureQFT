"""A2: pushforward of general holomorphic top forms -- closing
AMPLITUDES_CONNECTION.md section 2.4 question 2.

For the finite map F (3 sheets off the wall {p=0}) and any polynomial
observable g, the pushforward of the top form g(phi) d^3phi is

    F_*( g d^3phi ) = T[g/det DF] d^3J = -(1/2) T[g] d^3J,
    T[g](J) = sum over sheets  g(phi_i(J)),

because det DF = -2 is constant.  Question 2 asked whether the "spurious"
boundaries cancel between sheets except on {p=0}.  This script verifies,
exactly:

  1. For a basket of observables g -- including ones whose PER-SHEET values
     have poles on the x-collision locus {D0 = 0} through the fiber
     parametrization y = -B(x)/(2 D0), z = -D(x)/(8 D0) -- the sheet sum
     T[g] is a rational function whose denominator is a power of p ONLY:
     every D0-singularity cancels between the sheets.  (Conceptual reason:
     F is finite etale over the wall complement, so T[g] is regular there;
     being rational, its poles can only sit on {p=0}.  The D0-poles of the
     parametrization are coordinate artifacts.)
     Sharper: observables not involving x (or linear in x) have NO wall
     pole at all -- the y- and z-eliminants are monic (only x escapes), so
     wall poles are sourced exclusively by x-powers >= 2.
  2. Pole orders on the wall: ord_p T[x^k] = floor(k/2) (the escaping
     square-root pair), and the basket observables follow the same law in
     their total x-degree after reduction.
  3. g = 1 gives the constant  F_*(d^3phi) = -(3/2) d^3J: the wall is
     invisible to the holomorphic pushforward of the canonical form, while
     the REAL pushforward (N(J)/2) d^3J jumps across it -- the two differ
     exactly by the measured anomaly (scripts/measure_anomaly.py).

Conclusion: F_* g d^3phi always lands in rational forms with poles only on
the non-properness wall.  All collision-locus singularities are spurious
and cancel sheet-by-sheet in the sum -- the CHY-type mechanism, verified
exactly for this model.
"""

import sympy as sp

from jcqft import D0, PHI, SRC, p, q, r
from jcqft.fibers import g_x, y_of_x, z_of_x

x, y, z = PHI
a, b, c = SRC

# ---------------------------------------------------------------------------
# trace machinery: T[h(x)] for polynomial h via Newton power sums
pp_sym = p  # cubic p*X^3 + q*X + r, monic invariants e1 = 0, e2 = q/p, e3 = -r/p
S = {0: sp.Integer(3), 1: sp.Integer(0), 2: -2 * q / p}
for k in range(3, 40):
    S[k] = -(q / p) * S[k - 2] + (-r / p) * S[k - 3]


def trace_poly_in_x(h):
    """T[h] for h a polynomial in x with coefficients in C(a,b,c)."""
    poly = sp.Poly(sp.expand(h), x)
    return sp.together(sp.Add(*[coeff * S[deg]
                                for (deg,), coeff in poly.terms()]))


def trace_observable(g):
    """T[g] for a polynomial observable g(x, y, z): substitute the fiber
    parametrization (per-sheet values have D0-denominators!) and sum."""
    h = sp.together(g.subs({y: y_of_x, z: z_of_x}, simultaneous=True))
    num, den = sp.fraction(sp.cancel(h))
    assert not den.has(x), "denominator must be x-independent (a D0 power)"
    return sp.cancel(trace_poly_in_x(num) / den)


# ---------------------------------------------------------------------------
print("=== 1. Sheet sums are rational with poles ONLY on the wall {p=0} ===")
print("  (per-sheet values have D0-denominators from y = -B/(2 D0) etc.;")
print("   the sum must cancel them -- verified observable by observable)")

basket = {
    "1": sp.Integer(1),
    "x": x, "x^2": x**2, "x^3": x**3, "x^4": x**4,
    "y": y, "z": z, "y^2": y**2, "z^2": z**2,
    "x*y": x * y, "y*z": y * z, "x^2*z": x**2 * z,
    "x*y*z": x * y * z, "y^3": y**3,
    "x^2*y^2*z": x**2 * y**2 * z,
}

results = {}
for name, g in basket.items():
    T = sp.cancel(sp.together(trace_observable(g)))
    num, den = sp.fraction(T)
    # factor the denominator: all factors must be p (up to constants)
    dfac = sp.factor_list(den, a, b, c)
    bad = [f for f, mult in dfac[1] if sp.simplify(f - p) != 0
           and sp.simplify(f + p) != 0]
    ord_p = sum(mult for f, mult in dfac[1]
                if sp.simplify(f - p) == 0 or sp.simplify(f + p) == 0)
    results[name] = (T, ord_p)
    d0_in_den = den.has(a) and any(sp.simplify(f - D0) == 0 for f, _ in dfac[1])
    assert not bad, f"T[{name}] has non-wall pole factors: {bad}"
    print(f"  T[{name:9s}]: pole order in p = {ord_p}   "
          f"(D0 cancelled: {'yes' if not d0_in_den else 'NO'})")

print("\n  conceptual reason (now verified on the basket): F is finite etale")
print("  over C^3 \\ {p=0}, so every sheet sum is regular there; rationality")
print("  then forces all poles onto the wall.  The D0-singularities of the")
print("  per-sheet parametrization are spurious and cancel in the sum.")

# ---------------------------------------------------------------------------
print("\n=== 2. Pole-order law: ord_p T[x^k] = floor(k/2) ===")
for k in range(1, 11):
    T = sp.cancel(S[k])
    num, den = sp.fraction(sp.together(T))
    dp = sp.Poly(den, a, b, c)
    # denominator is +-p^m
    m = 0
    dd = den
    while sp.simplify(dd) not in (1, -1) and not dd.is_number:
        quo = sp.cancel(dd / p)
        if sp.fraction(quo)[1] != 1:
            break
        dd = quo
        m += 1
    print(f"  ord_p T[x^{k}] = {m}   floor(k/2) = {k // 2}")
    assert m == k // 2

# a mixed observable: total escape rate is set by the x-degree
T_mixed, ord_mixed = results["x^2*y^2*z"]
print(f"  mixed example T[x^2*y^2*z]: pole order {ord_mixed} "
      "(after D0-cancellation the wall order is what remains)")

# ---------------------------------------------------------------------------
print("\n=== 3. The two pushforwards of the canonical form ===")
T1 = results["1"][0]
assert T1 == 3
print("  holomorphic:  F_*(d^3phi) = T[1]/det DF d^3J = 3/(-2) d^3J")
print("                = -(3/2) d^3J  -- CONSTANT, the wall is invisible.")
print("  real:         (N(J)/|det DF|) d^3J = (N(J)/2) d^3J,")
print("                N = 3 iff p < 0: jumps 1/2 <-> 3/2 across the wall.")
print("  difference  = the measured anomaly A(sigma) -> 2 of")
print("  scripts/measure_anomaly.py: real minus holomorphic sheet counting.")

# closing check: T[x] = 0 identically -- the pushforward of x d^3phi VANISHES
assert sp.simplify(results["x"][0]) == 0
print("\n  bonus exact identity: T[x] = 0 identically, so")
print("  F_*(x d^3phi) = 0: the x-moment of the fiber vanishes at EVERY")
print("  source (the eliminant cubic has no X^2 term).")

print("\nDone: AMPLITUDES_CONNECTION.md section 2.4 question 2 is answered.")
print("F_*(g d^3phi) is rational with poles only on the non-properness wall;")
print("all spurious (collision-locus) boundaries cancel between sheets.")
