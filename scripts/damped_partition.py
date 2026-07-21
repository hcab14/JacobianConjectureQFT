"""The damped 0D partition function of the counterexample and the wall's
signature in its semiclassical limit (item B5 of docs/OPEN_QUESTIONS.md,
posed in docs/QFT_IMPLICATIONS.md section 5.3):

    Z_hbar(J) = integral over R^3 of exp(-|F(phi) - J|^2 / (2 hbar)) d^3 phi .

Main results (write-up: docs/DAMPED_PARTITION.md):

1. EXACT closed form and well-definedness.  |det DF| = 2 everywhere, so the
   real pushforward of Lebesgue measure is (N(J')/2) d^3 J' (chamber rule,
   scripts/measure_anomaly.py).  Hence for EVERY J in R^3 and EVERY hbar > 0

       Z_hbar(J) = (2 pi hbar)^{3/2} * h(J, hbar),
       h(J, hbar) = E[ N(J + sqrt(hbar) xi) ] / 2
                  = 1/2 + P[ p(J + sqrt(hbar) xi) < 0 ],   xi ~ N(0, 1_3),

   and therefore
       (1/2) (2 pi hbar)^{3/2}  <  Z_hbar(J)  <  (3/2) (2 pi hbar)^{3/2}.
   Z is finite and BOUNDED for all J -- including on the Jelonek set {p=0}
   and on the empty-fiber cusp orbit.  Tube volumes grow polynomially and
   uniformly in J:  (2 pi/3) t^3 <= V_J(t) = vol{|F-J| <= t} <= 2 pi t^3.
   The damped partition function does NOT diverge on the non-properness set;
   its entire wall signature sits in the bounded prefactor h(J, hbar).

2. SEMICLASSICS PER CHAMBER: h(J, hbar) -> N(J)/2 as hbar -> 0 for J off the
   wall, at the Gaussian rate exp(-dist(J, wall)^2 / (2 hbar)).  The leading
   semiclassical coefficient is the chamber function.

3. WALL AND CUSP VALUES (Laplace fails AT the wall, by escaping mass, not by
   Morse degeneration -- det DF = -2 forbids finite coalescence):
     generic wall point:  h = 1 + c1 sqrt(hbar) + O(hbar),
        c1 = -div(grad p/|grad p|) / (2 sqrt(2 pi))   (wall mean curvature);
        note h -> 1 != N(J_wall)/2 = 1/2: the escaped pair carries mass 1/2.
     perturbative vacuum J = 0 (ON the wall): h -> 1 != 1/2 = N(0)/2.
     cusp orbit (empty fiber, N = 0):  h = 1/2 + kappa hbar^{1/4} + ...,
        anomalous quarter-power from the A2 horn; kappa computed exactly
        from the cusp normal form and the source-fluctuation covariance.

4. UNIFORMITY BOUNDARY hbar*(eps) (where the semiclassical value is half-lost,
   h = 1.25, at distance eps from the wall inside the N=3 chamber):
     wall-normal approach:    hbar* = eps^2 / Phi^{-1}(3/4)^2       (gamma = 2),
     medial approach to cusp: hbar* = 27 eps^3/(826 Phi^{-1}(7/8)^2) (gamma = 3),
     wall-tangent approach:   hbar* ~ eps^4  (= dist^2 again; dist ~ eps^2).
   Unifying statement: hbar* ~ dist(J, wall)^2; the cusp exponent 3 = 2*(3/2)
   is the A2 horn-width exponent.

Exact claims carry asserts; numerics carry convergence tables.
Runtime: a few minutes.  Usage:  .venv/bin/python scripts/damped_partition.py
"""

import math
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import sympy as sp

from jcqft.core import D0, F, PHI, SRC, X, p, q, r
from jcqft.fibers import F_num, coef_num, yz_num

T0 = time.time()
x, y, z = PHI
a, b, c = SRC

# ===========================================================================
print("=== 1. EXACT: closed form, boundedness, tube-volume growth ===")

DF = sp.Matrix([[sp.diff(Fi, v) for v in PHI] for Fi in F])
assert sp.simplify(DF.det()) == -2
print("  det DF = -2 (so |det DF| = 2; each sheet contributes 1/2)")

# F is linear in z: F = A(x,y) z + B(x,y).
coeffs = [sp.Poly(Fi, z).all_coeffs() for Fi in F]
assert all(len(co) == 2 for co in coeffs), "F must be linear in z"
Avec = sp.Matrix([co[0] for co in coeffs])
Bvec = sp.Matrix([co[1] for co in coeffs])
alpha = sp.expand((Avec.T * Avec)[0])
# alpha = |A|^2 is a sum of even powers with no common real zero:
assert sp.expand(alpha - ((1 + x*y)**6 + 9*x**2*(1 + x*y)**4 + x**6)) == 0
print("  F = A(x,y) z + B(x,y),  A = ((1+xy)^3, 3x(1+xy)^2, -x^3)")
print("  |A|^2 = (1+xy)^6 + 9x^2(1+xy)^4 + x^6 > 0 everywhere (x=0 => first")
print("  term = 1): the z-integral of Z is an honest 1D Gaussian for all")
print("  (x, y), giving the exact fibration")
print("    Z = sqrt(2 pi hbar) * int dx dy |A|^{-1/2} exp(-m(x,y;J)/(2 hbar)),")
print("    m = |A x d|^2 / |A|^2,  d = B(x,y) - J   (Lagrange identity).")
d1, d2, d3 = sp.symbols("d1 d2 d3")
dv = sp.Matrix([d1, d2, d3])
cross = Avec.cross(dv)
assert sp.expand(alpha*(dv.T*dv)[0] - ((Avec.T*dv)[0])**2
                 - (cross.T*cross)[0]) == 0
print("  Lagrange identity verified symbolically.")

# Chamber-rule ingredient (as in scripts/measure_anomaly.py):
assert sp.expand(4*q**3 + 27*p*r**2 - 4*D0**2) == 0
print("  4q^3 + 27 p r^2 = 4 D0^2  =>  N(J) = 3 iff p < 0, else 1 (a.e.).")
print()
print("  PROPOSITION (exact).  The change-of-variables formula for the")
print("  non-injective local diffeomorphism F (|det DF| = 2, real preimage")
print("  count N) pushes Lebesgue measure to (N(J')/2) d^3J'.  Applying it")
print("  to the Gaussian g(J') = exp(-|J'-J|^2/(2 hbar)):")
print()
print("    Z_hbar(J) = (2 pi hbar)^{3/2} * h(J,hbar),")
print("    h(J,hbar) = E[N(J + sqrt(hbar) xi)]/2 = 1/2 + P[p(J+sqrt(hbar) xi)<0].")
print()
print("  Since 1 <= N <= 3 a.e.:  (2 pi hbar)^{3/2}/2 < Z < 3(2 pi hbar)^{3/2}/2")
print("  for ALL J (Jelonek set and cusp orbit included) and ALL hbar > 0.")
print("  Tube volumes: V_J(t) = int_{B(J,t)} (N/2) d^3J' is squeezed between")
print("  (2 pi/3) t^3 and 2 pi t^3 -- polynomial growth, uniform in J.")
print("  NO divergence detects the Jelonek set; the signature is in h.")

# ===========================================================================
# Numerical machinery for h(J, hbar): the a-integral is done exactly
# (p is quadratic in a), the (b,c)-integral by panel Gauss-Legendre.

_erfc = np.frompyfunc(math.erfc, 1, 1)


def Phi_cdf(t):
    """Standard normal CDF, vectorized."""
    return 0.5 * _erfc(-np.asarray(t) / math.sqrt(2)).astype(float)


def gl_nodes(lo, hi, npan, ngl):
    xg, wg = np.polynomial.legendre.leggauss(ngl)
    ed = np.linspace(lo, hi, npan + 1)
    mid = 0.5 * (ed[:-1] + ed[1:])
    half = 0.5 * np.diff(ed)
    return ((mid[:, None] + half[:, None] * xg[None, :]).ravel(),
            (half[:, None] * wg[None, :]).ravel())


def prob_p_neg(J, hb, npan=32, ngl=8, R=8.5):
    """P[p(J + sqrt(hb) xi) < 0], xi ~ N(0,1_3): exact in a, GL in (b, c)."""
    s = math.sqrt(hb)
    t, wq = gl_nodes(-R, R, npan, ngl)
    Bn = J[1] + s * t
    Cn = J[2] + s * t
    W = wq * np.exp(-t**2 / 2) / math.sqrt(2 * math.pi)
    Bg, Cg = np.meshgrid(Bn, Cn, indexing="ij")
    # p = A2*a^2 + A1*a + A0 with
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
    Pq = np.where(Disc > 0, Phi_cdf((hi - J[0]) / s) - Phi_cdf((lo - J[0]) / s),
                  0.0)
    Pl = np.where(A1 > 0, Phi_cdf((lin - J[0]) / s),
                  1.0 - Phi_cdf((lin - J[0]) / s))
    P2 = np.where(A2 > 1e-280, Pq, Pl)
    return float(np.einsum("i,j,ij->", W, W, P2))


def h_val(J, hb, **kw):
    return 0.5 + prob_p_neg(J, hb, **kw)


p_np = sp.lambdify(SRC, p, "numpy")
fA = sp.lambdify((x, y), list(Avec), "numpy")
fB = sp.lambdify((x, y), list(Bvec), "numpy")


def z_direct(J, hb, S=6.0, nx=2401, Ylim=12.0, ny=2401):
    """Direct quadrature of Z/(2 pi hbar)^{3/2}: exact Gaussian z-integral,
    sinh-mapped trapezoid in x (captures the escape tails), trapezoid in y."""
    s = np.linspace(-S, S, nx)
    Xn = np.sinh(s)
    wX = np.cosh(s) * (s[1] - s[0])
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
        al = A1v*A1v + A2v*A2v + A3v*A3v
        e1, e2, e3 = B1v - J[0], B2v - J[1], B3v - J[2]
        cx = A2v*e3 - A3v*e2
        cy = A3v*e1 - A1v*e3
        cz = A1v*e2 - A2v*e1
        m = (cx*cx + cy*cy + cz*cz) / al
        total += float(np.einsum("i,ij->", wb, np.exp(-m/(2*hb))/np.sqrt(al)))
    return total * dY * math.sqrt(2 * math.pi * hb) / (2 * math.pi * hb)**1.5


# ===========================================================================
print("\n=== 2. Independent numerical check of the closed form ===")
print("  Direct phi-space quadrature (no change of variables) vs the exact h:")
J3a = (-0.25, 0.0, 0.0)      # triple point image, N = 3
J1a = (1.0, 0.0, 0.0)        # N = 1
Jw = (0.0, 1.0, 1.0)         # generic wall point: p = 0, q = 1, D0 = -1
Jc = (4/27, 4/3, 1.0)        # cusp orbit: empty fiber, N = 0
for J, hb in [(J3a, 0.05), (J1a, 0.05), (Jw, 0.05), (Jc, 0.05), (J3a, 0.2)]:
    zd = z_direct(J, hb)
    he = h_val(J, hb, npan=48, ngl=12)
    rel = zd / he - 1
    print(f"    J = {J},  hbar = {hb}:  direct = {zd:.6f},  "
          f"h = {he:.6f},  rel.err = {rel:+.1e}")
    assert abs(rel) < 5e-4, rel
print("  All within 5e-4 (grid resolution).")

print("\n  Cutoff convergence at the WALL point (finiteness is not an")
print("  artifact of truncation; escaping-tube cross-section |A|^{-1/2} ~ ")
print("  |x|^{-3} makes the tail integrable):")
prev = None
for S in (4.0, 5.0, 6.0, 7.0):
    zd = z_direct(Jw, 0.1, S=S)
    tag = "" if prev is None else f"   delta = {zd - prev:+.2e}"
    print(f"    |x| <= sinh({S}) = {math.sinh(S):9.1f}:  "
          f"Z/(2 pi hbar)^{{3/2}} = {zd:.7f}{tag}")
    if prev is not None:
        last_delta = abs(zd - prev)
    prev = zd
assert last_delta < 1e-3 * prev
print("  Stable under cutoff growth: Z_hbar(J_wall) is finite.")

# ===========================================================================
print("\n=== 3. Semiclassics per chamber: prefactor -> N(J)/2 ===")
chamber_pts = [("(-1/4,0,0)  N=3", J3a, 3), ("(0,2,0)     N=3", (0.0, 2.0, 0.0), 3),
               ("(1,0,0)     N=1", J1a, 1), ("(2,1,1)     N=1", (2.0, 1.0, 1.0), 1)]
# verify N by direct root count (as in scripts/measure_anomaly.py)
for name, J, N in chamber_pts:
    pc, qc, rc = coef_num(*J)
    roots = np.roots([pc, 0.0, qc, rc])
    n_direct = 0
    for rt in roots:
        if abs(rt.imag) < 1e-9:
            yr, zr = yz_num(rt.real, *J)
            assert np.linalg.norm(np.array(F_num(rt.real, yr, zr))
                                  - np.array(J)) < 1e-8
            n_direct += 1
    assert n_direct == N and (3 if pc < 0 else 1) == N
print("  Direct preimage counts match the chamber rule at all 4 test points.")
hb_grid = (0.3, 0.1, 0.03, 0.01, 3e-3, 1e-3)
hdr = "  " + " ".join(f"{hb:>9.0e}" for hb in hb_grid)
print("\n  Convergence table:  h(J, hbar) = Z/(2 pi hbar)^{3/2}   ->   N(J)/2")
print(f"  {'J':>16s}" + hdr + "    limit")
for name, J, N in chamber_pts:
    row = [h_val(J, hb) for hb in hb_grid]
    print(f"  {name:>16s}  " + " ".join(f"{v:9.6f}" for v in row)
          + f"    {N/2:.1f}")
    assert abs(row[-1] - N/2) < 1e-4, (name, row[-1])
print("  At hbar = 1e-3 all four agree with N(J)/2 to < 1e-4; the rate is")
print("  the Gaussian tail exp(-dist(J,wall)^2/(2 hbar)) -- e.g. dist = 1/4")
print("  for (-1/4,0,0) gives exp(-31) at hbar = 1e-3.")
print("  => The leading semiclassical normalization IS the chamber function:")
print("     it jumps 1/2 <-> 3/2 across the Jelonek set, a jump invisible to")
print("     the hbar-expansion around any single vacuum (all corrections to")
print("     Laplace around one nondegenerate minimum are powers of hbar with")
print("     smooth coefficients; the jump is a boundary-of-field-space effect).")

# ===========================================================================
print("\n=== 4. On the wall: h -> 1, with an exact sqrt(hbar) coefficient ===")
# escaping-pair location at distance eps from the wall (normal direction):
gradp = sp.Matrix([sp.diff(p, v) for v in SRC])
gp_w = np.array([float(g.subs(dict(zip(SRC, Jw)))) for g in gradp])
n_w = gp_w / np.linalg.norm(gp_w)           # points into p > 0 (N=1 side)
for eps in (1e-3, 1e-5):
    J = np.array(Jw) - eps * n_w            # inside N=3 chamber
    pc, qc, rc = coef_num(*J)
    roots = sorted(np.roots([pc, 0.0, qc, rc]), key=lambda t: -abs(t))[:2]
    for sgn, rt in zip((1, -1), sorted(roots, key=lambda t: -t.real)):
        pred = sgn * math.sqrt(-qc / pc) + rc / (2 * qc)
        # next correction is O(eps) relative
        assert abs(rt.imag) < 1e-9 and abs(rt.real / pred - 1) < 30 * eps
print("  Escaping pair at J = J_wall - eps*n:  x = +-sqrt(-q/p) + r/(2q)")
print("  verified (to the O(eps) next order) at eps = 1e-3, 1e-5;  |x| ~")
print("  (q/|grad p|)^{1/2} eps^{-1/2}: SQUARE-ROOT escape.  No finite minima")
print("  merge (det DF = -2 forbids finite Morse degeneration);")
print("  Laplace 'fails' at the wall only through this escaping mass.")
divn = sum(sp.diff(gradp[i] / sp.sqrt((gradp.T*gradp)[0]), v)
           for i, v in enumerate(SRC))
divn_w = float(divn.subs(dict(zip(SRC, Jw))))
c1_ana = -divn_w / (2 * math.sqrt(2 * math.pi))
print(f"\n  Curvature coefficient: h(J_wall,hbar) = 1 + c1 sqrt(hbar) + O(hbar),")
print(f"  c1 = -div(grad p/|grad p|)/(2 sqrt(2 pi)) = {c1_ana:+.6f} at Jw={Jw}")
for hb in (1e-4, 1e-5, 1e-6):
    hm = h_val(Jw, hb, npan=48, ngl=12)
    c1_m = (hm - 1) / math.sqrt(hb)
    print(f"    hbar = {hb:.0e}:  (h-1)/sqrt(hbar) = {c1_m:+.6f}"
          f"   (residual h-1-c1*sqrt(hbar) = {hm-1-c1_ana*math.sqrt(hb):+.1e})")
assert abs(c1_m / c1_ana - 1) < 5e-3
print("  Measured coefficient matches the analytic one to < 0.5%.")
print("  NOTE h -> 1 = (3+1)/2/2 * 2 -- the mean of the chamber values -- and")
print("  NOT N(J_wall)/2 = 1/2: the escaped pair still carries Gaussian mass")
print("  1/2.  Same at the vacuum J = 0 (which lies ON the wall):")
row = [h_val((0.0, 0.0, 0.0), hb) for hb in (0.1, 0.01, 1e-3, 1e-4, 1e-5)]
print("    h(0, hbar) at hbar = 0.1 ... 1e-5:  "
      + "  ".join(f"{v:.5f}" for v in row))
assert abs(row[-1] - 1) < 5e-3
print("  => Z_hbar(0)/(2 pi hbar)^{3/2} -> 1 = TWICE the perturbative-saddle")
print("     value 1/2: the damped partition function at the perturbative")
print("     vacuum sees the two sheets at infinity, which no order of the")
print("     tree expansion sees.")

# ===========================================================================
print("\n=== 5. At the cusp orbit (empty fiber): anomalous hbar^{1/4} ===")
u_, w_, du_, dw_ = sp.symbols("u w du dw")
P2 = 27*u_**2 + 16*u_ - 18*u_*w_ + w_**3 - w_**2
assert sp.expand(c**2 * p - P2.subs({u_: a*c**2, w_: b*c})) == 0
P2s = sp.expand(P2.subs({u_: sp.Rational(4, 27) + du_,
                         w_: sp.Rational(4, 3) + dw_}))
quad = sum(t for t in P2s.as_ordered_terms()
           if sp.Poly(t, du_, dw_).total_degree() == 2)
assert sp.expand(quad - 3*(3*du_ - dw_)**2) == 0
assert sp.expand(P2s.subs(dw_, 3*du_) - 27*du_**3) == 0
print("  A2 normal form at the cusp (4/27, 4/3) of the reduced wall P2:")
print("    P2 = 3 (3 du - dw)^2 + 27 du^3 + higher order,")
print("  so the N=3 horn at depth s = -dw (medial line dw = 3 du) has")
print("  half-width |v| < (s/3)^{3/2} in v = du - dw/3: exponent 3/2.")
# Exact medial-line factorization (cuspidal-tangent stratum, c = 1):
eps_s = sp.Symbol("varepsilon", positive=True)
med = {a: (9*(sp.Rational(4, 3) - eps_s) - 8) / 27,
       b: sp.Rational(4, 3) - eps_s, c: 1}
assert sp.expand(p.subs(med) + eps_s**3) == 0
assert sp.expand(q.subs(med) - 3*eps_s) == 0 and r.subs(med) == -2
assert sp.expand((p.subs(med)*X**3 + q.subs(med)*X + r.subs(med))
                 + (eps_s*X - 1)**2 * (eps_s*X + 2)) == 0
print("  Medial line (c=1, w = 4/3 - eps):  p = -eps^3, q = 3 eps, r = -2,")
print("    cubic = -(eps X - 1)^2 (eps X + 2):  fiber x = 1/eps (double),")
print("    -2/eps -- ALL THREE preimages escape like eps^{-1} (empty fiber")
print("    in the limit), vs the eps^{-1/2} pair at a generic wall point.")
# Source-fluctuation covariance in cusp coordinates (linearized at Jc):
#   du = da + (8/27) dc,  dw = db + (4/3) dc,  v = du - dw/3.
Vw = sp.Rational(1) + sp.Rational(4, 3)**2
Vv = sp.Rational(1) + sp.Rational(1, 3)**2 + sp.Rational(4, 27)**2
Cvw = -sp.Rational(1, 3) + sp.Rational(-4, 27)*sp.Rational(4, 3)
assert (Vw, Vv, Cvw) == (sp.Rational(25, 9), sp.Rational(826, 729),
                         sp.Rational(-43, 81))
sig2c = sp.nsimplify(Vv - Cvw**2 / Vw)
assert sig2c == sp.Rational(18801, 18225)
Vwf, kf, s2f = float(Vw), float(Cvw / Vw), float(sig2c)
tq = np.linspace(1e-6, 40.0, 400001)
kappa_ana = (2/(3*math.sqrt(3))) * float(np.trapezoid(
    tq**1.5
    * np.exp(-tq**2/(2*Vwf)) / math.sqrt(2*math.pi*Vwf)
    * np.exp(-(kf*tq)**2/(2*s2f)) / math.sqrt(2*math.pi*s2f), tq))
print(f"\n  h(J_cusp, hbar) = 1/2 + kappa hbar^{{1/4}} + ...,   analytic")
print(f"  kappa = (2/(3 sqrt 3)) int t^{{3/2}} phi_{{25/9}}(t) "
      f"phi_{{{s2f:.4f}}}(43t/225) dt = {kappa_ana:.6f}")
for hb in (1e-5, 1e-6, 1e-7, 1e-8):
    km = (h_val(Jc, hb, npan=64, ngl=12) - 0.5) / hb**0.25
    print(f"    hbar = {hb:.0e}:  (h - 1/2)/hbar^{{1/4}} = {km:.5f}")
assert abs(km / kappa_ana - 1) < 0.01
print("  Measured amplitude matches to < 1%.  The prefactor limit at the")
print("  cusp is 1/2, NOT N/2 = 0: the Gaussian mass of a source with NO")
print("  classical solution comes entirely from the tube at infinity.")

# ===========================================================================
print("\n=== 6. Uniformity boundary hbar*(eps): measured exponents gamma ===")


def hstar(J, target=1.25, lo=1e-16, hi=1.0, iters=40):
    if not (h_val(J, lo) > target > h_val(J, hi)):
        return None
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        if h_val(J, mid) > target:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


z75 = float(mp.sqrt(2) * mp.erfinv(mp.mpf(1) / 2))   # Phi^{-1}(3/4)
print(f"  Definition: hbar*(J) solves h(J, hbar*) = 1.25 (half the 3/2 - 1")
print(f"  gap lost).  Flat-wall prediction for normal approach: h = 1/2 +")
print(f"  Phi(eps/sqrt(hbar)) => hbar* = eps^2/Phi^{{-1}}(3/4)^2 "
      f"= {1/z75**2:.5f} eps^2.")
print("\n  (a) normal approach to the generic wall point Jw:")
wall_scan = []
for k in range(3, 9):
    eps = 10**(-k/2)
    J = tuple(np.array(Jw) - eps * n_w)
    assert p_np(*J) < 0
    hs = hstar(J)
    wall_scan.append((eps, hs))
    print(f"    eps = {eps:.3e}:  hbar* = {hs:.4e}   hbar*/eps^2 = "
          f"{hs/eps**2:.5f}")
le = np.log([e for e, _ in wall_scan])
lh = np.log([h for _, h in wall_scan])
gamma_wall = np.polyfit(le[-3:], lh[-3:], 1)[0]
print(f"    fit exponent gamma_wall = {gamma_wall:.4f}  (predicted 2)")
assert abs(gamma_wall - 2) < 0.02
assert abs(wall_scan[-1][1] / wall_scan[-1][0]**2 * z75**2 - 1) < 5e-3
print(f"    constant matches 1/Phi^{{-1}}(3/4)^2 = {1/z75**2:.5f} to < 0.5%.")

print("\n  (b) medial (cuspidal-tangent) approach to the cusp, c = 1:")
cusp_scan = []
for k in range(3, 11):
    eps = 10**(-k/3)
    wv = 4/3 - eps
    J = ((9*wv - 8)/27, wv, 1.0)
    assert p_np(*J) < 0
    hs = hstar(J)
    cusp_scan.append((eps, hs))
    print(f"    eps = {eps:.3e}:  hbar* = {hs:.4e}   hbar*/eps^3 = "
          f"{hs/eps**3:.5f}")
le = np.log([e for e, _ in cusp_scan])
lh = np.log([h for _, h in cusp_scan])
gamma_cusp = np.polyfit(le[-4:], lh[-4:], 1)[0]
C_cusp = cusp_scan[-1][1] / cusp_scan[-1][0]**3
z875 = float(mp.sqrt(2) * mp.erfinv(mp.mpf(3) / 4))   # Phi^{-1}(7/8)
C_pred = (1 / 27) / (float(Vv) * z875**2)
print(f"    fit exponent gamma_cusp = {gamma_cusp:.4f}  (horn prediction 3)")
print(f"    constant: measured C = {C_cusp:.5f}; asymptotic model")
print(f"    C = (horn half-width)^2/eps^3 / (Var(v) Phi^{{-1}}(7/8)^2)")
print(f"      = (1/27) / ((826/729) Phi^{{-1}}(7/8)^2) = {C_pred:.5f}")
print("    [transverse fluctuation v = da - db/3 - (4/27) dc must fit in the")
print("     horn half-width eps^{3/2}/sqrt(27); depth fluctuation ~ sqrt(hbar*)")
print("     ~ eps^{3/2} << eps is subleading, so the 1D crossing rule is exact.]")
assert abs(gamma_cusp - 3) < 0.02
assert abs(C_cusp / C_pred - 1) < 0.01

print("\n  (c) tangential approach to Jw (dist ~ eps^2): expect hbar* ~ eps^4")
t1 = np.array([1.0, 0.0, 0.0]) - n_w * n_w[0]
t1 /= np.linalg.norm(t1)
tan_scan = []
for eps in (0.1, 0.06, 0.03):
    J = tuple(np.array(Jw) + eps * t1)
    if p_np(*J) >= 0:
        continue
    hs = hstar(J)
    tan_scan.append((eps, hs))
    print(f"    eps = {eps:.2f}:  hbar* = {hs:.3e}   hbar*/eps^4 = "
          f"{hs/eps**4:.4f}")
slope_tan = np.polyfit(np.log([e for e, _ in tan_scan]),
                       np.log([h for _, h in tan_scan]), 1)[0]
print(f"    fit exponent = {slope_tan:.2f} (~4: consistent with the unified")
print("    rule hbar* ~ dist(J, wall)^2 in every regime; gamma_cusp = 3 is")
print("    2 x (3/2), the A2 horn-width exponent).")
assert 3.6 < slope_tan < 4.4

# ===========================================================================
print("\n=== 7. Figure ===")
os.makedirs("outputs", exist_ok=True)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.8))
hbs = np.logspace(-8, 0, 25)
curves = [("$J=(-1/4,0,0)$, $N=3$", J3a, "tab:blue"),
          ("$J=(1,0,0)$, $N=1$", J1a, "tab:orange"),
          ("$J=0$ (vacuum, on wall)", (0.0, 0.0, 0.0), "tab:green"),
          ("$J_w=(0,1,1)$ (generic wall)", Jw, "tab:red"),
          ("$J_c=(4/27,4/3,1)$ (cusp, $N=0$)", Jc, "tab:purple")]
for lab, J, col in curves:
    ax1.semilogx(hbs, [h_val(J, hb) for hb in hbs], color=col, label=lab)
for v in (0.5, 1.0, 1.5):
    ax1.axhline(v, color="gray", lw=0.6, ls=":")
ax1.set_xlabel(r"$\hbar$")
ax1.set_ylabel(r"$Z_\hbar(J)\,/\,(2\pi\hbar)^{3/2}$")
ax1.set_title("Semiclassical prefactor: chambers, wall, cusp")
ax1.legend(fontsize=8, loc="upper left")
ew, hw = zip(*wall_scan)
ec, hc = zip(*cusp_scan)
ax2.loglog(ew, hw, "o-", label=r"wall-normal: $\hbar^*=\varepsilon^2/"
           r"\Phi^{-1}(3/4)^2$ ($\gamma=2$)")
ax2.loglog(ec, hc, "s-", label=r"cusp medial: $\hbar^*=C\,\varepsilon^3$ "
           r"($\gamma=3$)")
ax2.loglog(ew, np.array(ew)**2 / z75**2, "k--", lw=0.8)
ax2.loglog(ec, C_cusp * np.array(ec)**3, "k:", lw=0.8)
ax2.set_xlabel(r"$\varepsilon$ (distance parameter to wall/cusp)")
ax2.set_ylabel(r"$\hbar^*(\varepsilon)$")
ax2.set_title(r"Uniformity boundary of the $\hbar$-expansion")
ax2.legend(fontsize=8)
fig.tight_layout()
fig.savefig("outputs/damped_partition.png", dpi=160)
print("  saved outputs/damped_partition.png")

# ===========================================================================
print("\n=== 8. Synthesis ===")
print("""  The damped partition function Z_hbar(J) is a constructively finite,
  hbar-uniformly bounded, real-analytic function of the source for every
  hbar > 0 -- there is nothing to renormalize and nothing diverges, on or
  off the Jelonek set.  Its semiclassical expansion has leading term
  (2 pi hbar)^{3/2} (N(J)/2): a PIECEWISE-CONSTANT prefactor that jumps
  across the non-properness set {p=0}, equals the two-sided mean 1 on the
  wall (vacuum included), and drops to the anomalous value 1/2 + O(hbar^{1/4})
  on the empty-fiber cusp orbit.  Because det DF = -2 is constant, every
  local (perturbative) datum of the model is chamber-independent -- the
  standard lore 'constant Jacobian => trivial semiclassics' holds for each
  vacuum separately -- yet the semiclassical normalization is a step
  function of J.  The hbar-expansion is uniform only for
  hbar << dist(J, wall)^2 (gamma_wall = 2, constant 1/Phi^{-1}(3/4)^2 for
  the half-loss criterion), degrading to hbar << C eps^3 on approach to the
  cusp (gamma_cusp = 3 = 2 x 3/2, the A2 horn exponent).  All of this is a
  0D statement about one map; no claim is made about D >= 1.""")

print(f"\nALL CHECKS PASSED   (runtime {time.time()-T0:.1f} s)")
