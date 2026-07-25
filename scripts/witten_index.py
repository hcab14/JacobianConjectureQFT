"""The 0D Witten index of the counterexample: signed solution count,
Mathai-Quillen localization, and index jumping at the non-properness wall
(write-up: docs/WITTEN_INDEX.md).

Main results:

1. EXACT: signed count = Brouwer degree per chamber.  det DF = -2
   identically, so sign det DF(phi) = -1 at EVERY real solution of
   F(phi) = J, and the signed solution count ("0D Witten index") is
       deg(F, J) = sum over F^{-1}(J) of sign det DF = -N(J)
   with N the real preimage count: deg = -1 on {p > 0}, deg = -3 on
   {p < 0}.  Asserted via exact fibers at rational points per chamber
   (each preimage verified to map back EXACTLY, reduction mod the minimal
   polynomial of the root).

2. EXACT CERTIFICATE OF NON-PROPERNESS: every J is a regular value
   (det DF = -2 never vanishes) and every fiber is finite, so deg(F, .)
   is defined everywhere; for a PROPER map R^3 -> R^3 it would be a
   single constant (extension to S^3 -> S^3).  The jump -1 -> -3 across
   {p = 0} therefore certifies non-properness, independently of the
   escape-curve certificate of scripts/missing_observables.py (re-asserted
   here in one line).  Converse link: {p = 0} IS the Jelonek set
   (scripts/branch_locus.py).

3. EXACT: the Mathai-Quillen / bosonic partition function
       Z_sigma(J) = (2 pi sigma^2)^{-3/2} int det DF e^{-|F-J|^2/(2 sigma^2)}
   has the closed form  Z_sigma(J) = -E[N(J + sigma xi)]  (xi ~ N(0,1_3)),
   hence is FINITE for every J and sigma > 0, with -3 < Z_sigma(J) < -1,
   and off the wall Z_sigma(J) -> deg(F, J) at the exact Gaussian rate
   e^{-dist(J, wall)^2/(2 sigma^2)} (two-sided bound).  On the wall the
   limit is -2, at the empty-fiber cusp orbit -1: boundary contributions
   from solutions at infinity, exactly quantified.

4. EXACT SUSY STRUCTURE:  DF is NOT symmetric (one-line assert), so no
   superpotential W with F = grad W exists (no Parisi-Sourlas model).
   The Mathai-Quillen completion exists for arbitrary F: fields
   (phi, psi, chibar, B), BRST charge delta(phi) = psi, delta(psi) = 0,
   delta(chibar) = B, delta(B) = 0 (nilpotent), action
       S = delta[ chibar . ( i(F(phi) - J) + (sigma^2/2) B ) ]
         = i B.(F - J) + (sigma^2/2) |B|^2 - i chibar . DF psi .
   Verified with an explicit 6-generator Grassmann algebra: delta^2 = 0,
   S = delta(Psi), delta(S) = 0, Berezin integral of e^{i chibar M psi}
   = -i det M (generic M), Gaussian B-integral, and the normalization
   N = i (2 pi)^{-3} that reduces the MQ integral to the bosonic Z_sigma.

5. NUMERICAL EVIDENCE: direct phi-space quadrature of the MQ integral
   (no change of variables) matches the closed form; convergence tables
   Z_sigma -> -1 / -3 per chamber; measured decay rate matches
   dist(J, wall)^2; wall / vacuum / cusp limits -2 / -2 / -1; the
   near-wall crossover matches the flat-wall profile -1 - 2 Phi(eps/sigma).

6. EXACT CONTRAST WITH C: the complex fiber count is 3 for every J off
   {p = 0} (no jump); the discriminant sign rule (monic disc =
   -4 D0^2 / p^3) decides 1 vs 3 REAL roots -- the wall is a real
   phenomenon, invisible to the complex count.

Exact claims carry asserts; numerics are labelled.  Runtime ~10 s.
Usage:  .venv/bin/python scripts/witten_index.py
"""

import math
import time

import numpy as np
import sympy as sp

from jcqft.core import D0, F, PHI, SRC, X, cubic, p, q, r
from jcqft.fibers import exact_fiber

T0 = time.time()
x, y, z = PHI
a, b, c = SRC

# ===========================================================================
print("=== 1. EXACT: the signed count is -N(J) (chamber-wise Brouwer degree) ===")

DF = sp.Matrix([[sp.diff(Fi, v) for v in PHI] for Fi in F])
assert sp.simplify(DF.det()) == -2
print("  det DF = -2 identically  =>  sign det DF(phi) = -1 at EVERY point;")
print("  every J is a regular value; each real solution contributes -1.")


def _zero_mod_minpoly(expr):
    """Exact zero test for a polynomial expression in (at most) one CRootOf:
    reduce the numerator modulo the root's minimal polynomial."""
    expr = sp.together(sp.expand(expr))
    roots = expr.atoms(sp.CRootOf)
    if not roots:
        return sp.simplify(expr) == 0
    assert len(roots) == 1
    (al,) = roots
    t = sp.Dummy("t")
    num, den = sp.fraction(expr.subs(al, t))
    assert not den.has(t)
    m = al.poly.as_expr().subs(al.poly.gens[0], t)
    return sp.expand(sp.rem(sp.expand(num), m, t)) == 0


# Rational sample targets, two per chamber, all with D0 != 0 (fiber
# parametrization valid) and p != 0 (off the wall).
SAMPLES = [
    ((sp.Rational(-1, 4), 0, 0), -1),   # p = -4  < 0
    ((0, 2, 0), -1),                    # p = -4  < 0
    ((1, 0, 0), +1),                    # p = 16  > 0
    ((2, 1, 1), +1),                    # p = 104 > 0
]
degrees = {}
for J, sgn_p in SAMPLES:
    sub = dict(zip(SRC, J))
    pJ, d0J = p.subs(sub), D0.subs(sub)
    assert sp.sign(pJ) == sgn_p and d0J != 0
    fib = exact_fiber(J)
    assert len(fib) == 3, "complex fiber count"
    # every preimage maps back EXACTLY (reduction mod minimal polynomial)
    for pt in fib:
        psub = dict(zip(PHI, pt))
        assert all(_zero_mod_minpoly(Fi.subs(psub) - Ji)
                   for Fi, Ji in zip(F, J))
    n_real = sum(1 for pt in fib if all(v.is_real for v in pt))
    deg = -n_real                        # each solution contributes sign(-2)
    n_rule = 3 if pJ < 0 else 1
    assert n_real == n_rule, (J, n_real, n_rule)
    degrees[J] = deg
    print(f"  J = {sp.sstr(J):>14s}:  p = {sp.sstr(pJ):>4s},  D0 = "
          f"{sp.sstr(d0J):>2s},  complex fiber = 3,  real fiber = {n_real},"
          f"  deg(F, J) = {deg}")
assert sorted(set(degrees.values())) == [-3, -1]
print("  =>  deg(F, J) = -1 on {p > 0},  deg(F, J) = -3 on {p < 0}.")
print("  (Fiber exhaustiveness is exact: the eliminant cubic lies in the")
print("   ideal (F - J), and off {D0 = 0} the parametrization y = -B/2D0,")
print("   z = -D/8D0 is forced -- jcqft/fibers.py, MISSING_OBSERVABLES.md §4.)")

# ===========================================================================
print("\n=== 2. EXACT: the index jump certifies non-properness ===")
print("  DEGREE THEORY (standard; Outerelo-Ruiz GSM 108, Milnor): a PROPER")
print("  C^1 map f: R^3 -> R^3 extends to S^3 -> S^3 (infinity -> infinity),")
print("  and at every regular value the signed count equals the mapping")
print("  degree of the extension -- a single integer, independent of J.")
print("  Here every J is regular and every fiber is finite, yet the signed")
print("  count takes BOTH values -1 and -3  =>  F is NOT proper.  QED.")
print("  This certificate is independent of the escape-curve certificate")
print("  of scripts/missing_observables.py, re-asserted in one line:")

# escape curve: phi(T) -> infinity while F(phi(T)) -> finite wall point
T_, y0, c3 = sp.symbols("T y0 c3")
curve = {x: T_, y: y0, z: (2 * T_ - 3 * T_**2 * y0 - c3) / T_**3}
Flim = tuple(sp.limit(sp.simplify(Fi.subs(curve)), T_, sp.oo) for Fi in F)
Flim_expected = (y0**2 * (1 - c3 * y0), y0 * (4 - 3 * c3 * y0), c3)
assert all(sp.expand(u - v) == 0 for u, v in zip(Flim, Flim_expected))
Flim = Flim_expected
assert sp.expand(p.subs(dict(zip(SRC, Flim)))) == 0
print(f"    phi(T) = (T, y0, (2T - 3T^2 y0 - c3)/T^3):  |phi| -> oo,")
print(f"    F(phi(T)) -> {sp.sstr(Flim)}  -- a FINITE point, and p = 0 there.")
print("  Converse link: {p = 0} is exactly the Jelonek non-properness set")
print("  (escape only in x, leading coefficient of the x-eliminant = p;")
print("  scripts/branch_locus.py -- not re-derived here).")

# ===========================================================================
print("\n=== 3. EXACT: Mathai-Quillen partition function, closed form, "
      "finiteness ===")
# chamber-rule ingredient (as in scripts/measure_anomaly.py):
assert sp.expand(4 * q**3 + 27 * p * r**2 - 4 * D0**2) == 0
# finiteness fibration ingredient (as in scripts/damped_partition.py):
coeffs = [sp.Poly(Fi, z).all_coeffs() for Fi in F]
assert all(len(co) == 2 for co in coeffs), "F must be linear in z"
Avec = sp.Matrix([co[0] for co in coeffs])
Bvec = sp.Matrix([co[1] for co in coeffs])
alpha = sp.expand((Avec.T * Avec)[0])
assert sp.expand(alpha - ((1 + x*y)**6 + 9*x**2*(1 + x*y)**4 + x**6)) == 0
print("  Definition:  Z_sigma(J) = (2 pi sigma^2)^{-3/2}")
print("               * int_{R^3} det DF(phi) e^{-|F(phi)-J|^2/(2 sigma^2)} dphi.")
print("  det DF = -2 = -|det DF| comes out of the integral, so Z_sigma =")
print("  -2 (2 pi s^2)^{-3/2} Z^damped_{hbar = s^2}(J), and the EXACT pushforward")
print("  F_* d^3phi = (N/2) d^3J (measure_anomaly.py, DAMPED_PARTITION.md §1.1)")
print("  gives the closed form")
print()
print("      Z_sigma(J) = -E[ N(J + sigma xi) ] = -1 - 2 P[ p(J + sigma xi) < 0 ],")
print("      xi ~ N(0, 1_3).")
print()
print("  CONSEQUENCES (exact):")
print("   (i)   -3 < Z_sigma(J) < -1  for ALL J and ALL sigma > 0: the MQ")
print("         integral is FINITE unconditionally -- wall and cusp included")
print("         (the escaping tube has cross-section |A|^{-1/2} ~ |x|^{-3},")
print("         integrable; same mechanism as DAMPED_PARTITION.md §1.1).")
print("   (ii)  off the wall, |Z_sigma(J) - deg(F, J)| = 2 P[p flips sign]")
print("         <= 2 (1 + sqrt(2/pi) d/sigma) e^{-d^2/(2 sigma^2)},")
print("         d = dist(J, {p=0}):  Z_sigma -> deg(F, J) as sigma -> 0.")
print("   (iii) ON the wall the limit is the two-sided mean -2, NOT -N;")
print("         at the cusp orbit (N = 0) the limit is -1: ALL of it is a")
print("         boundary contribution from solutions at infinity.")
print("  The sigma -> 0 localization is therefore valid PER CHAMBER and")
print("  fails exactly on the non-properness set -- the jump of the index")
print("  and the non-properness of F are the same fact.")

# ===========================================================================
print("\n=== 4. EXACT: SUSY structure -- what exists and what does not ===")
asym = sp.expand(DF[0, 2] - DF[2, 0])
assert asym != 0 and asym.subs({x: 0, y: 0, z: 0}) == -1
print(f"  (a) DF - DF^T != 0:  (DF)_13 - (DF)_31 = {sp.sstr(asym)}")
print("      (nonzero even at phi = 0).  So NO superpotential W with")
print("      F = grad W exists, even locally: the Parisi-Sourlas action")
print("      |grad W - J|^2/2 + psibar Hess W psi is NOT available.")

print("\n  (b) The Mathai-Quillen completion exists for arbitrary F.")
print("      Fields: phi in R^3, auxiliary B in R^3, fermions psi, chibar.")
print("      BRST:  d(phi_i) = psi_i, d(psi_i) = 0, d(chibar_i) = B_i,")
print("      d(B_i) = 0;  S = d[ chibar.( i(F - J) + (s^2/2) B ) ].")
print("      Verified in an explicit 6-generator Grassmann algebra:")


# ---- minimal Grassmann (Berezin) algebra --------------------------------
# generators: psi_1..3 -> indices 0,1,2;  chibar_1..3 -> indices 3,4,5.
# element = {tuple(sorted generator indices): sympy coefficient}

def gmul(u, v):
    out = {}
    for ku, cu in u.items():
        for kv, cv in v.items():
            if set(ku) & set(kv):
                continue
            merged = ku + kv
            sign = 1
            for i in range(len(merged)):
                for j in range(i + 1, len(merged)):
                    if merged[i] > merged[j]:
                        sign = -sign
            key = tuple(sorted(merged))
            out[key] = out.get(key, 0) + sign * cu * cv
    return {k: sp.expand(cc) for k, cc in out.items() if sp.expand(cc) != 0}


def gadd(u, v):
    out = dict(u)
    for k, cc in v.items():
        out[k] = sp.expand(out.get(k, 0) + cc)
    return {k: cc for k, cc in out.items() if cc != 0}


def gscale(u, s):
    return {k: sp.expand(s * cc) for k, cc in u.items() if sp.expand(s * cc) != 0}


def gexp(u):
    """exp of an even element (terminates: 6 generators)."""
    res = {(): sp.Integer(1)}
    term = {(): sp.Integer(1)}
    for n in range(1, 4):
        term = gscale(gmul(term, u), sp.Rational(1, n))
        if not term:
            break
        res = gadd(res, term)
    return res


Bsym = sp.symbols("B1 B2 B3")
Jsym = sp.symbols("J1 J2 J3")
sig = sp.Symbol("sigma", positive=True)
PSIg = [{(i,): sp.Integer(1)} for i in range(3)]
CHIg = [{(3 + i,): sp.Integer(1)} for i in range(3)]
DELTA_RULES = {i: None for i in range(3)}            # d(psi_i) = 0
DELTA_RULES.update({3 + i: Bsym[i] for i in range(3)})  # d(chibar_i) = B_i


def delta(u):
    """The odd derivation: d(phi_i) = psi_i on coefficients, plus the
    generator rules; Leibniz with Grassmann signs."""
    out = {}
    for key, coeff in u.items():
        # d acting on the coefficient: sum_i (d coeff/d phi_i) psi_i * theta
        for i, v in enumerate(PHI):
            dc = sp.diff(coeff, v)
            if dc != 0:
                out = gadd(out, gmul({(i,): dc}, {key: sp.Integer(1)}))
        # d acting on the generators (all odd, images even)
        for pos, g in enumerate(key):
            img = DELTA_RULES[g]
            if img is None:
                continue
            newkey = key[:pos] + key[pos + 1:]
            sgn = (-1) ** pos
            out = gadd(out, {newkey: sgn * coeff * img})
    return out


# gauge fermion and action
Psi_gf = {}
for i in range(3):
    Psi_gf = gadd(Psi_gf, gscale(
        CHIg[i], sp.I * (F[i] - Jsym[i]) + sig**2 / 2 * Bsym[i]))
S_impl = delta(Psi_gf)

S_explicit = {(): sp.expand(
    sum(sp.I * Bsym[i] * (F[i] - Jsym[i]) for i in range(3))
    + sig**2 / 2 * sum(Bi**2 for Bi in Bsym))}
for i in range(3):
    for j in range(3):
        S_explicit = gadd(S_explicit,
                          gscale(gmul(CHIg[i], PSIg[j]), -sp.I * DF[i, j]))
assert gadd(S_impl, gscale(S_explicit, -1)) == {}
print("        S = delta(Psi) = i B.(F-J) + (s^2/2)|B|^2 - i chibar.DF psi  [assert]")
assert delta(S_impl) == {}
assert delta(delta(Psi_gf)) == {}
test_elt = {(0, 4): x**2 * y * z, (1,): x * y, (): z**3}
assert delta(delta(test_elt)) == {}
print("        delta(S) = 0  and  delta^2 = 0 (on Psi and on a generic test")
print("        element)  [assert]  -- the action is BRST-exact, delta-closed.")

# Berezin integral of the fermion weight, generic matrix
M = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"M{i}{j}"))
bil = {}
for i in range(3):
    for j in range(3):
        bil = gadd(bil, gscale(gmul(CHIg[i], PSIg[j]), sp.I * M[i, j]))
top = gexp(bil).get((0, 1, 2, 3, 4, 5), sp.Integer(0))
assert sp.expand(top - sp.I**3 * M.det()) == 0
print("        Berezin:  int dpsi dchibar e^{i chibar.M psi} = -i det M")
print("        (coefficient of psi1 psi2 psi3 chibar1 chibar2 chibar3;")
print("        generic 3x3 M)  [assert]")

# Gaussian B-integral, one component
Bk, dk = sp.symbols("Bk dk", real=True)
IB = sp.integrate(sp.exp(-sp.I * Bk * dk - sig**2 * Bk**2 / 2),
                  (Bk, -sp.oo, sp.oo))
assert sp.simplify(IB - sp.sqrt(2 * sp.pi) / sig
                   * sp.exp(-dk**2 / (2 * sig**2))) == 0
print("        int dB e^{-i B d - s^2 B^2/2} = (sqrt(2 pi)/s) e^{-d^2/2s^2}  [assert]")

# normalization: N * (-i) * (2 pi / s^2)^{3/2} = (2 pi s^2)^{-3/2}
Nmq = sp.I * (2 * sp.pi) ** (-3)
assert sp.simplify(Nmq * (-sp.I) * (2 * sp.pi / sig**2) ** sp.Rational(3, 2)
                   - (2 * sp.pi * sig**2) ** sp.Rational(-3, 2)) == 0
print("        normalization N = i (2 pi)^{-3}:  N int dphi dB dpsi dchibar")
print("        e^{-S} = (2 pi s^2)^{-3/2} int det DF e^{-|F-J|^2/2s^2} dphi")
print("        = Z_sigma(J)  [assert]  -- the MQ model reduces EXACTLY to")
print("        the bosonic Z_sigma of section 3; its 'index' is deg(F, J).")

# ===========================================================================
print("\n=== 5. NUMERICAL EVIDENCE: Z_sigma converges to the index ===")
# semi-analytic evaluator for P[p(J + sigma xi) < 0]: the a-marginal is
# exact (p is quadratic in a), Gauss-Legendre panels in (b, c)
# -- adapted from scripts/damped_partition.py.

_erfc = np.frompyfunc(math.erfc, 1, 1)


def Phi_cdf(t):
    if np.isscalar(t):
        return 0.5 * math.erfc(-t / math.sqrt(2))
    return 0.5 * _erfc(-np.asarray(t) / math.sqrt(2)).astype(float)


def gl_nodes(lo, hi, npan, ngl):
    xg, wg = np.polynomial.legendre.leggauss(ngl)
    ed = np.linspace(lo, hi, npan + 1)
    mid = 0.5 * (ed[:-1] + ed[1:])
    half = 0.5 * np.diff(ed)
    return ((mid[:, None] + half[:, None] * xg[None, :]).ravel(),
            (half[:, None] * wg[None, :]).ravel())


def prob_p_neg(J, s, npan=32, ngl=8, R=8.5):
    t, wq = gl_nodes(-R, R, npan, ngl)
    Bn = J[1] + s * t
    Cn = J[2] + s * t
    W = wq * np.exp(-t**2 / 2) / math.sqrt(2 * math.pi)
    Bg, Cg = np.meshgrid(Bn, Cn, indexing="ij")
    A2 = 27.0 * Cg**2
    A1 = 16.0 - 18.0 * Bg * Cg
    A0 = Bg**3 * Cg - Bg**2
    Disc = A1 * A1 - 4 * A2 * A0
    with np.errstate(all="ignore"):
        sq = np.sqrt(np.maximum(Disc, 0.0))
        qq = -0.5 * (A1 + np.sign(A1) * sq)
        r1 = qq / A2
        r2 = np.where(qq != 0, A0 / qq, 0.0)
        lo = np.minimum(r1, r2)
        hi = np.maximum(r1, r2)
        lin = -A0 / A1
    Pq = np.where(Disc > 0,
                  Phi_cdf((hi - J[0]) / s) - Phi_cdf((lo - J[0]) / s), 0.0)
    Pl = np.where(A1 > 0, Phi_cdf((lin - J[0]) / s),
                  1.0 - Phi_cdf((lin - J[0]) / s))
    P2 = np.where(A2 > 1e-280, Pq, Pl)
    return float(np.einsum("i,j,ij->", W, W, P2))


def Z_mq(J, s, **kw):
    """Closed form: Z_sigma(J) = -1 - 2 P[p(J + sigma xi) < 0]."""
    return -1.0 - 2.0 * prob_p_neg(J, s, **kw)


fA = sp.lambdify((x, y), list(Avec), "numpy")
fB = sp.lambdify((x, y), list(Bvec), "numpy")


def Z_mq_direct(J, s, S=6.0, nx=2401, Ylim=12.0, ny=2401):
    """Direct phi-space quadrature of the MQ integral (no pushforward):
    det DF = -2 out front, exact Gaussian z-integral, sinh-mapped x grid
    (escape tails), trapezoid in y -- adapted from damped_partition.py."""
    hb = s * s
    u = np.linspace(-S, S, nx)
    Xn = np.sinh(u)
    wX = np.cosh(u) * (u[1] - u[0])
    wX[0] *= 0.5
    wX[-1] *= 0.5
    Y = np.linspace(-Ylim, Ylim, ny)
    dY = Y[1] - Y[0]
    total = 0.0
    for i0 in range(0, nx, 256):
        Xb = Xn[i0:i0 + 256][:, None]
        wb = wX[i0:i0 + 256]
        A1v, A2v, A3v = fA(Xb, Y[None, :])
        B1v, B2v, B3v = fB(Xb, Y[None, :])
        A1v, A2v, A3v, B1v, B2v, B3v = np.broadcast_arrays(
            A1v, A2v, A3v, B1v, B2v, B3v)
        al = A1v * A1v + A2v * A2v + A3v * A3v
        e1, e2, e3 = B1v - J[0], B2v - J[1], B3v - J[2]
        cx = A2v * e3 - A3v * e2
        cy = A3v * e1 - A1v * e3
        cz = A1v * e2 - A2v * e1
        m = (cx * cx + cy * cy + cz * cz) / al
        total += float(np.einsum("i,ij->", wb,
                                 np.exp(-m / (2 * hb)) / np.sqrt(al)))
    h = total * dY * math.sqrt(2 * math.pi * hb) / (2 * math.pi * hb) ** 1.5
    return -2.0 * h


print("  5.1 Direct phi-space quadrature of the MQ integral vs closed form")
print("      (validates the integral itself, independent of the pushforward):")
J3a = (-0.25, 0.0, 0.0)
J1a = (1.0, 0.0, 0.0)
Jw = (0.0, 1.0, 1.0)          # generic wall point: p = 0, q = 1, D0 = -1
Jc = (4 / 27, 4 / 3, 1.0)     # cusp orbit: empty fiber, N = 0
for J, s in [(J3a, math.sqrt(0.05)), (J1a, math.sqrt(0.05)),
             (Jw, math.sqrt(0.05))]:
    zd = Z_mq_direct(J, s)
    ze = Z_mq(J, s, npan=48, ngl=12)
    rel = zd / ze - 1
    print(f"      J = {J}, sigma = {s:.4f}:  direct = {zd:.6f},  "
          f"closed form = {ze:.6f},  rel.err = {rel:+.1e}")
    assert abs(rel) < 5e-4, rel
print("      All within 5e-4 (grid resolution).")

print("\n  5.2 Convergence table:  Z_sigma(J)  ->  deg(F, J)")
sig_grid = (0.5, 0.3, 0.2, 0.1, 0.05, 0.03)
chamber_pts = [("(-1/4,0,0)", J3a, -3), ("(0,2,0)", (0.0, 2.0, 0.0), -3),
               ("(1,0,0)", J1a, -1), ("(2,1,1)", (2.0, 1.0, 1.0), -1)]
print(f"      {'J':>12s}  " + " ".join(f"s={s:<5g}" for s in sig_grid)
      + "   deg")
for name, J, dg in chamber_pts:
    row = [Z_mq(J, s) for s in sig_grid]
    print(f"      {name:>12s}  " + " ".join(f"{v:7.4f}" for v in row)
          + f"   {dg}")
    assert abs(row[-1] - dg) < 1e-3, (name, row[-1])
print("      At sigma = 0.03 all four match deg(F, J) to < 1e-3.")

print("\n  5.3 Decay rate: fitted d^2 in |Z + N| ~ e^{-d^2/(2 sigma^2)} vs")
print("      the measured distance to the wall {p = 0} (grid refinement):")


def wall_dist(J, R=3.0, n=201, levels=4):
    """min |J' - J| over {p = 0}: grid over (b, c), a-roots of the
    quadratic p(a; b, c) = 0, iteratively refined."""
    cb, cc = J[1], J[2]
    best = math.inf
    for _ in range(levels):
        bs = np.linspace(cb - R, cb + R, n)
        cs = np.linspace(cc - R, cc + R, n)
        Bg, Cg = np.meshgrid(bs, cs, indexing="ij")
        A2 = 27.0 * Cg**2
        A1 = 16.0 - 18.0 * Bg * Cg
        A0 = Bg**3 * Cg - Bg**2
        with np.errstate(all="ignore"):
            Disc = A1 * A1 - 4 * A2 * A0
            sq = np.sqrt(np.maximum(Disc, 0.0))
            d2best = np.full(Bg.shape, np.inf)
            for sgn in (1.0, -1.0):
                aa = np.where(np.abs(A2) > 1e-12, (-A1 + sgn * sq) / (2 * A2),
                              np.where(sgn > 0, -A0 / A1, np.nan))
                ok = np.isfinite(aa) & ((Disc >= 0) | (np.abs(A2) <= 1e-12))
                d2 = (aa - J[0])**2 + (Bg - J[1])**2 + (Cg - J[2])**2
                d2best = np.minimum(d2best, np.where(ok, d2, np.inf))
        k = int(np.argmin(d2best))
        i, j = divmod(k, n)
        best = min(best, math.sqrt(float(d2best[i, j])))
        cb, cc = bs[i], cs[j]
        R *= 2.5 / n * 4
    return best


for name, J, dg in [("(1,0,0)", J1a, -1), ("(-1/4,0,0)", J3a, -3)]:
    d_grid = wall_dist(J)
    ss = np.geomspace(d_grid / 2.0, d_grid / 5.5, 8)
    errs = np.array([abs(Z_mq(J, s, npan=48, ngl=12) - dg) for s in ss])
    slope = np.polyfit(1.0 / (2.0 * ss[-4:] ** 2), np.log(errs[-4:]), 1)[0]
    d_fit = math.sqrt(-slope)
    print(f"      J = {name}:  d_grid = {d_grid:.4f},  d_fit = {d_fit:.4f}"
          f"   (ratio {d_fit/d_grid:.3f})")
    assert abs(d_fit / d_grid - 1) < 0.15
print("      Fitted Gaussian-decay distances match dist(J, wall) to < 15%")
print("      (polynomial prefactor causes the residual drift; the EXACT")
print("      two-sided bound is in section 3(ii)).")

print("\n  5.4 Wall, vacuum, cusp: the sigma -> 0 limit off the chambers")
for name, J, lim in [("generic wall (0,1,1)", Jw, -2.0),
                     ("vacuum J = 0 (ON the wall)", (0.0, 0.0, 0.0), -2.0),
                     ("cusp (4/27, 4/3, 1), N = 0", Jc, -1.0)]:
    row = [Z_mq(J, s, npan=48, ngl=12) for s in (0.1, 0.01, 1e-3)]
    print(f"      {name:>28s}:  " + "  ".join(f"{v:.5f}" for v in row)
          + f"   -> {lim}")
    assert abs(row[-1] - lim) < 0.02, (name, row[-1])
print("      Wall: -2 (mean of -1 and -3), NOT -N(J_wall) = -1: the escaped")
print("      pair contributes -1 from infinity.  Vacuum: same (J = 0 is ON")
print("      the wall).  Cusp: -1 with an EMPTY fiber -- the whole value is")
print("      a boundary term.  (Rates: c1 sqrt(hbar) and kappa hbar^{1/4},")
print("      hbar = sigma^2 -- DAMPED_PARTITION.md §1.3, not re-derived.)")

print("\n  5.5 Approaching the wall: crossover at sigma ~ dist(J, wall)")
gradp = sp.Matrix([sp.diff(p, v) for v in SRC])
gp_w = np.array([float(g.subs(dict(zip(SRC, Jw)))) for g in gradp])
n_w = gp_w / np.linalg.norm(gp_w)     # normal, points into p > 0
eps = 1e-3
Jeps = tuple(np.array(Jw) - eps * n_w)   # N = 3 side, dist ~ eps
print(f"      J = J_wall - eps*n, eps = {eps:g} (deg = -3 there); flat-wall")
print("      model: Z ~ -1 - 2 Phi(eps/sigma)  (curvature correction O(sigma)):")
for ratio in (5.0, 1.0, 0.2):
    s = ratio * eps
    zv = Z_mq(Jeps, s, npan=48, ngl=12)
    pred = -1.0 - 2.0 * Phi_cdf(eps / s)
    print(f"        sigma = {ratio:3g} eps:  Z = {zv:.4f}   flat-wall pred"
          f" {pred:.4f}")
    assert abs(zv - pred) < 0.01
print("      => the index is resolved only once sigma << dist(J, wall):")
print("      the localization scale must beat the distance to the wall")
print("      (uniformity boundary hbar* ~ dist^2, DAMPED_PARTITION.md §2.3).")

# ===========================================================================
print("\n=== 6. EXACT: contrast with the complex count (no wall over C) ===")
disc = sp.discriminant(cubic, X)
assert sp.expand(disc - (-4 * D0**2 * p)) == 0
print("  disc_X(p X^3 + q X + r) = -4 D0^2 p; monic discriminant -4 D0^2/p^3:")
print("  sign = -sign(p)  =>  3 real roots iff p < 0, 1 iff p > 0 (chamber")
print("  rule, scripts/measure_anomaly.py).  Over C the fiber count is 3 for")
print("  EVERY J off {p = 0} (asserted at all four sample points in section")
print("  1, both chambers); on the wall it degenerates by ESCAPE, not")
print("  collision (det DF = -2 is etale):")
fib_wall = exact_fiber((sp.Rational(2, 27), 1, 1))
fib_cusp = exact_fiber((sp.Rational(4, 27), sp.Rational(4, 3), 1))
assert len(fib_wall) == 1 and len(fib_cusp) == 0
print("    wall point (2/27, 1, 1):  1 preimage (pair escaped);")
print("    cusp point (4/27, 4/3, 1):  0 preimages (all three escaped).")
print("  The complex 'index' (unsigned count 3) sees no wall; the REAL")
print("  signed count jumps.  The wall is a real-locus phenomenon: exactly")
print("  where the discriminant of the real cubic changes sign.")

# ===========================================================================
print("\n=== 7. Synthesis ===")
print("""  The 0D Witten index of the counterexample -- the SUSY-localized
  partition function of its Mathai-Quillen model -- is the Brouwer degree
  deg(F, J) = -N(J): it equals -1 on {p > 0} and -3 on {p < 0}, and JUMPS
  across the Jelonek set {p = 0}.  For a proper map this index would be a
  single integer; its non-invariance here is an exact certificate of
  non-properness, and conversely the non-properness locus is exactly where
  the sigma -> 0 localization acquires boundary terms (wall value -2, cusp
  value -1 over an empty fiber).  No superpotential exists (DF != DF^T);
  the MQ completion is the honest SUSY structure, its action is BRST-exact
  and its 'index' is deg(F, J).  Interpretation (flagged): two classical
  vacua escape through infinity as J crosses the wall -- the 0D caricature
  of Witten-index jumping via vacua at infinity and of wall-crossing.
  All of this is a 0D statement about one map; no claim about D >= 1.""")

print(f"\nALL CHECKS PASSED   (runtime {time.time() - T0:.1f} s)")
