"""The field-redefinition "measure anomaly" of the counterexample, made
quantitative over the reals.

F restricts to a real map R^3 -> R^3 with constant Jacobian -2: a local
diffeomorphism everywhere.  Standard practice would treat phi -> F(phi) as a
legitimate change of variables ("equivalence theorem": the Jacobian is a
harmless constant).  Globally this is wrong: the change-of-variables formula
for a non-injective local diffeomorphism reads

    integral f(F(phi)) |det DF| d^3phi  =  integral f(J) N(J) d^3J ,

where N(J) is the number of real preimages.  This script establishes:

  1. CHAMBER RULE (exact):  N(J) = 3 if p(a,b,c) < 0, and N(J) = 1 if
     p > 0.  Proof: the real preimage count equals the number of real roots
     of the x-eliminant p X^3 + q X + r (y, z are real rational functions of
     the root), and the monic discriminant is -(4q^3 + 27 p r^2)/p^3 =
     -4 D0^2 / p^3 with D0 = 27ac^2 - 9bc + 8, so its sign is -sign(p).
  2. Numerical spot-checks of the rule against direct preimage computation.
  3. Monte Carlo measurement of the anomaly factor
        A(sigma) = <N(J)> = 1 + 2 P[p < 0]
     for Gaussian source ensembles of width sigma.  Because the perturbative
     vacuum J = 0 lies ON the wall {p = 0} (p has no constant term and its
     linear part is 16a), the anomaly does NOT switch off near the vacuum:
     A(sigma) -> 2 as sigma -> 0.
  4. Monodromy-invariant ("sum over all sheets") observables are RATIONAL:
     e1 = x1+x2+x3 = 0, e2 = q/p, e3 = -r/p, with poles exactly on {p = 0}.
"""

import numpy as np
import sympy as sp

from counterexample import F, PHI, SRC, X, cubic, p, q, r

a, b, c = SRC
D0 = 27 * a * c**2 - 9 * b * c + 8

# ---------------------------------------------------------------------------
print("=== 1. Exact chamber rule for the real preimage count N(J) ===")
identity_ok = sp.expand(4 * q**3 + 27 * p * r**2 - 4 * D0**2) == 0
print(f"  4q^3 + 27 p r^2 = 4 D0^2 (symbolic): {identity_ok}")
print("  monic discriminant of the x-cubic: -4 D0^2 / p^3  ->  sign = -sign(p)")
print("  =>  N(J) = 3  iff  p(a,b,c) < 0;   N(J) = 1  iff  p > 0")
print(f"  p has no constant term (p(0)= {p.subs(dict(zip(SRC,(0,0,0))))}), "
      "linear part 16a: the vacuum J = 0 sits ON the wall {p=0}.")

# ---------------------------------------------------------------------------
print("\n=== 2. Spot-checks: chamber rule vs direct real-preimage count ===")
F_num = sp.lambdify(PHI, F, "numpy")
p_num = sp.lambdify(SRC, p, "numpy")
coef_num = sp.lambdify(SRC, (p, q, r), "numpy")

# rational y,z parametrization from the lex Groebner basis
gb = sp.groebner([F[0] - a, F[1] - b, F[2] - c],
                 PHI[1], PHI[2], PHI[0], order="lex")
g_y = next(g for g in gb.exprs if sp.degree(g, PHI[1]) == 1)
g_z = next(g for g in gb.exprs if sp.degree(g, PHI[2]) == 1)
A_, B_ = sp.Poly(g_y, PHI[1]).all_coeffs()
C_, D_ = sp.Poly(g_z, PHI[2]).all_coeffs()
yz_num = sp.lambdify((PHI[0],) + SRC, (-B_ / A_, -D_ / C_), "numpy")

rng = np.random.default_rng(0)
mismatches = 0
for _ in range(300):
    J = rng.normal(0, 1.0, 3)
    pc, qc, rc = coef_num(*J)
    roots = np.roots([pc, 0.0, qc, rc])
    real_roots = [rt.real for rt in roots if abs(rt.imag) < 1e-9]
    n_direct = 0
    for xr in real_roots:
        yr, zr = yz_num(xr, *J)
        resid = np.linalg.norm(np.array(F_num(xr, yr, zr)) - J)
        assert resid < 1e-6, f"parametrization residual {resid}"
        n_direct += 1
    n_rule = 3 if pc < 0 else 1
    if n_direct != n_rule:
        mismatches += 1
print(f"  300 random targets (sigma = 1): direct count vs rule mismatches: "
      f"{mismatches}")
assert mismatches == 0

# ---------------------------------------------------------------------------
print("\n=== 3. Monte Carlo: anomaly factor A(sigma) = <N(J)> = 1 + 2 P[p<0] ===")
print("  (naive change of variables predicts A = 1 identically)")
NSAMP = 10_000_000
print(f"  {NSAMP:.0e} Gaussian samples per width:")
for sigma in (10.0, 1.0, 0.1, 0.01, 0.001):
    Jv = rng.normal(0, sigma, (3, NSAMP))
    frac = float(np.mean(p_num(*Jv) < 0))
    A_sigma = 1 + 2 * frac
    err = 2 * np.sqrt(frac * (1 - frac) / NSAMP)
    print(f"    sigma = {sigma:7.3f}:  A = {A_sigma:.4f} +- {err:.4f}"
          f"   (P[p<0] = {frac:.4f})")
print("  -> A -> 2 as sigma -> 0: the defect is O(1) for sources concentrated")
print("     ARBITRARILY CLOSE to the perturbative vacuum, because J = 0 lies")
print("     on the escape wall {p = 0}.  It fades (A -> 1) only for wide")
print("     ensembles dominated by the quartic term 27 a^2 c^2 > 0 in p.")

# ---------------------------------------------------------------------------
print("\n=== 4. Sum over all sheets restores single-valuedness ===")
print("  Elementary symmetric functions of the three x-roots (all sheets):")
print("    e1 = x1 + x2 + x3 = 0   (no X^2 term in the cubic)")
print(f"    e2 = q/p = {sp.sstr(sp.simplify(q / p))}")
print(f"    e3 = -r/p = {sp.sstr(sp.simplify(-r / p))}")
print("  Monodromy-invariant 'sum over all vacua' observables are RATIONAL in")
print("  the sources, with poles exactly on the non-properness locus {p=0} --")
print("  the single perturbative branch is multivalued, their sum is not.")
