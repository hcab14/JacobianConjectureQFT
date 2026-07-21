"""Wall complement = A2-discriminant complement: the braid-group structure
of the Alpoge-Mathew model, verified exactly.

Companion script of docs/WALL_COMPLEMENT.md; answers item B2 of
docs/OPEN_QUESTIONS.md (identify the correct "amplituhedron-analogue"
object: the wall complement with its S3 local system).

Chain of results (exact/symbolic unless marked numerical):

  1. INVARIANT ELIMINANT (I1).  On the chart c != 0, in the C*-invariants
     u = a c^2, w = b c and the weight-0 root variable xi = X / c, the
     x-eliminant  p X^3 + q X + r  becomes, after dividing by c,
         E(xi; u, w) = P2(u,w) xi^3 + (4 - 3w) xi - 2,
     with P2 = c^2 p = 27 u^2 + 16 u - 18 u w + w^3 - w^2 the plane wall.
  2. TAUTOLOGICAL MAP TO DEPRESSED CUBICS (I2).  The scaling eta = P2 xi
     makes the family monic and depressed:
         eta^3 + Qhat eta + Rhat,   (Qhat, Rhat) = ((4-3w) P2, -2 P2^2),
     a POLYNOMIAL map from the invariant plane to the universal space
     {x^3 + Q x + R} of depressed cubics.
  3. DISCRIMINANT PULLBACK (I3).
         4 Qhat^3 + 27 Rhat^2 = 4 P2^3 (27u - 9w + 8)^2,
     driven by the key lemma  (4-3w)^3 + 27 P2 = (27u - 9w + 8)^2.
     The wall {P2 = 0} maps to the CUSP (Q,R) = (0,0) of the universal
     discriminant; the D0-line contributes the second factor.
  4. AFFINE EQUIVALENCE (I4).  The wall is affinely equivalent to the
     standard A2 discriminant curve (semicubical parabola):
         U = (4 - 3w)/9,   W = (27u - 9w + 8)/27
         ==>  P2 = 27 (W^2 - U^3)     (exact polynomial identity).
     Equivalently (Q, R) = (-3U, 2W) is an affine ISOMORPHISM of the
     invariant plane onto the space of depressed cubics carrying the wall
     exactly onto {4Q^3 + 27R^2 = 0}, with multiplicity ONE:
         4 Q^3 + 27 R^2 = 4 P2.
     By Arnold/Brieskorn/Deligne the discriminant complement is K(B3, 1),
     so  pi_1(C^2 \\ {P2 = 0}) = B3, the braid group on three strands.
  5. NUMERICAL MONODROMY in the invariant plane (mpmath root tracking of
     the invariant eliminant): wall meridians -> transpositions, D0-line
     meridians -> id, and a radius-0.05 loop around the CUSP -> a 3-CYCLE
     of order 3 (the image of sigma1*sigma2, a Coxeter element of W(A2));
     a torus-type loop (the full twist (sigma1 sigma2)^3) -> id.  All of
     it matches the universal family x^3 + Q x + R computed side by side.
  6. LOCAL-SYSTEM BOOKKEEPING: the rank-3 sheet local system decomposes
     as trivial (+) standard 2-dim (the reflection representation of
     S3 = Weyl(A2)); the trivial summand is spanned by the trace
     observables (rational, poles only on the wall — the invariant-plane
     version of scripts/trace_pushforward.py); the Euler characteristic
     of the wall complement is  chi(C^2 \\ {P2=0}) = 1 - 1 = 0.

Run:  .venv/bin/python scripts/wall_braid.py     (~1 minute)
"""

import itertools
import time

import mpmath as mp
import sympy as sp

from jcqft import D0, SRC, X, p, q, r

T0 = time.time()
a, b, c = SRC
u, w, xi, eta = sp.symbols("u w xi eta")
m, m1, m2, s = sp.symbols("m m1 m2 s")

P2 = 27 * u**2 + 16 * u - 18 * u * w + w**3 - w**2
qt = 4 - 3 * w                    # q = 4 - 3bc in the invariants
D0r = 27 * u - 9 * w + 8          # D0 (weight 0 already) in the invariants
CUSP = (sp.Rational(4, 27), sp.Rational(4, 3))

# ---------------------------------------------------------------------------
print("=== 1. Invariant eliminant and the tautological map to depressed"
      " cubics ===")

# weight bookkeeping: under (a,b,c,X) ~ (-2,-1,1,1) the eliminant is
# quasi-homogeneous of weight 1, so dividing by c (weight 1) after X = xi*c
# (xi of weight 0) must land in the weight-0 invariants u = a c^2, w = b c.
WTS = {a: -2, b: -1, c: 1, X: 1}
elim = p * X**3 + q * X + r
wset = {sum(e * WTS[g] for e, g in zip(mono, WTS))
        for mono, _ in sp.Poly(elim, *WTS).terms()}
print(f"  weights of the monomials of p*X^3 + q*X + r: {sorted(wset)}")
assert wset == {1}

# reductions of the building blocks
assert sp.expand(c**2 * p - P2.subs({u: a * c**2, w: b * c})) == 0
assert sp.expand(q - qt.subs(w, b * c)) == 0 and sp.expand(r + 2 * c) == 0
assert sp.expand(D0.subs({a: u / c**2, b: w / c}) - D0r) == 0

# (I1): substitute X = xi*c, divide by c
E_inv = sp.expand(elim.subs(X, xi * c) / c)
E_target = P2.subs({u: a * c**2, w: b * c}) * xi**3 + (4 - 3 * b * c) * xi - 2
assert sp.expand(E_inv - E_target) == 0
E = P2 * xi**3 + qt * xi - 2
print("  (I1) VERIFIED:  (p X^3 + q X + r)|_{X = xi c} / c")
print("       = P2(u,w) xi^3 + (4 - 3w) xi - 2   with u = a c^2, w = b c.")

# its discriminant in xi: the invariant-plane shadow of disc_X = -4 D0^2 p
disc_xi = sp.expand(sp.discriminant(E, xi))
assert sp.expand(disc_xi - (-4 * P2 * D0r**2)) == 0
print("  disc_xi(E) = -4 * P2 * (27u - 9w + 8)^2   (invariant form of")
print("  disc_X = -4 D0^2 p; branch locus {P2=0}, collision line {D0r=0}).")

# (I2): depress by eta = P2*xi (multiply E by P2^2)
Qhat = qt * P2
Rhat = -2 * P2**2
lhs = sp.expand(P2**2 * E)
rhs = sp.expand((eta**3 + Qhat * eta + Rhat).subs(eta, P2 * xi))
assert sp.expand(lhs - rhs) == 0
print("  (I2) VERIFIED:  P2^2 * E = eta^3 + Qhat*eta + Rhat  with")
print("       eta = P2*xi,  (Qhat, Rhat) = ((4-3w) P2, -2 P2^2):")
print("       a POLYNOMIAL map (u,w) -> (Qhat,Rhat) to depressed cubics.")

# (I3): pullback of the universal discriminant
lemma = sp.expand(qt**3 + 27 * P2 - D0r**2)
assert lemma == 0
assert sp.expand(4 * Qhat**3 + 27 * Rhat**2 - 4 * P2**3 * D0r**2) == 0
print("  (I3) VERIFIED:  4 Qhat^3 + 27 Rhat^2 = 4 P2^3 (27u - 9w + 8)^2,")
print("       via the exact lemma  (4-3w)^3 + 27 P2 = (27u - 9w + 8)^2.")
print("  => on the wall {P2=0}: (Qhat, Rhat) = (0,0) EXACTLY -- the wall is")
print("     contracted to the CUSP of the universal discriminant; the")
print("     D0-line supplies the second factor (x-collision = double root).")
# sanity: on the D0-line (off the wall) the image lies ON the universal
# discriminant but away from the origin (Qhat != 0 there generically)
onD0 = {u: (9 * w - 8) / 27}
assert sp.simplify((4 * Qhat**3 + 27 * Rhat**2).subs(onD0)) == 0
assert sp.simplify(Qhat.subs(onD0)) != 0

# ---------------------------------------------------------------------------
print("\n=== 2. Explicit affine equivalence with the A2 discriminant"
      " (I4) ===")
# Ansatz: (U, W) affine-linear in (u, w) with  P2 = K (W^2 - U^3),  K != 0.
K, a1, a2, a3, b1, b2, b3 = sp.symbols("K a1 a2 a3 b1 b2 b3")
Ua = a1 * u + a2 * w + a3
Wa = b1 * u + b2 * w + b3
mism = sp.Poly(sp.expand(P2 - K * (Wa**2 - Ua**3)), u, w)
# the u^3 coefficient of the mismatch is K a1^3, so a1 = 0 is forced:
assert sp.expand(mism.coeff_monomial(u**3) - K * a1**3) == 0
print("  u^3-coefficient of P2 - K(W^2 - U^3) is K*a1^3  =>  a1 = 0 forced.")
# residual freedom is the scaling (U,W) -> (t^2 U, t^3 W), K -> K t^6;
# normalize K = 27 (the rational branch) and solve the remaining system:
eqs = [e.subs({a1: 0, K: 27}) for e in mism.coeffs()]
sols = sp.solve(eqs, [a2, a3, b1, b2, b3], dict=True)
rat_sols = [so for so in sols if all(v.is_rational for v in so.values())]
print(f"  solutions with K = 27: {len(sols)} total, rational: {rat_sols}")
assert len(rat_sols) == 2       # (U, W) and (U, -W)
sol = next(so for so in rat_sols if so[b1] == 1)
Uex = Ua.subs({a1: 0, **sol})
Wex = Wa.subs(sol)
assert sp.expand(P2 - 27 * (Wex**2 - Uex**3)) == 0
print(f"  (I4) VERIFIED:  U = {Uex},  W = {Wex},")
print("                  P2 = 27 (W^2 - U^3)   (exact identity).")
assert sp.expand(9 * Uex - qt) == 0 and sp.expand(27 * Wex - D0r) == 0
print("  Structure:  U = (4 - 3w)/9 = q~/9   and   W = (27u - 9w + 8)/27")
print("  = D0r/27:  the eliminant coefficient q and the collision factor D0")
print("  ARE the universal coordinates, up to scale.  (I3)'s lemma is this")
print("  identity restated:  P2 = (D0r^2 - q~^3)/27 = 27 (W^2 - U^3).")

# the induced affine isomorphism onto the space of depressed cubics
Qe = sp.expand(-3 * Uex)
Re = sp.expand(2 * Wex)
assert sp.expand(4 * Qe**3 + 27 * Re**2 - 4 * P2) == 0
lin = sp.Matrix([[sp.diff(Qe, u), sp.diff(Qe, w)],
                 [sp.diff(Re, u), sp.diff(Re, w)]])
assert lin.det() != 0
print(f"  Induced AFFINE ISOMORPHISM  (u,w) -> (Q,R) = ({Qe}, {Re}):")
print(f"    4 Q^3 + 27 R^2 = 4 P2   (multiplicity ONE; det of linear part"
      f" = {lin.det()} != 0).")
print("  => (C^2, wall) ~ (C^2, A2 discriminant) as affine pairs; the wall")
print("     complement IS the discriminant complement of depressed cubics,")
print("     a K(B3,1) [Arnold 1969, Brieskorn 1971, Deligne 1972]:")
print("     pi_1(C^2 \\ {P2=0}) = B3, the braid group on 3 strands.")
# checkpoints of the affine map
assert Qe.subs({u: CUSP[0], w: CUSP[1]}) == 0
assert Re.subs({u: CUSP[0], w: CUSP[1]}) == 0
assert sp.simplify(Re.subs(onD0)) == 0
print("  cusp (4/27, 4/3) -> (Q,R) = (0,0);  D0-line -> {R = 0} (the")
print("  cuspidal tangent of the universal discriminant);  {w = 4/3} ->")
print("  {Q = 0}.")

# ---------------------------------------------------------------------------
print("\n=== 3. Numerical monodromy in the invariant plane ===")
mp.mp.dps = 30
PERMS3 = list(itertools.permutations(range(3)))


def solve_cubic(cf):
    return [mp.mpc(rr) for rr in mp.polyroots(cf, maxsteps=300,
                                              extraprec=120)]


def coeffs_inv(uv, wv):
    """Invariant eliminant P2 xi^3 + (4-3w) xi - 2 at a numeric point."""
    Pv = 27 * uv**2 + 16 * uv - 18 * uv * wv + wv**3 - wv**2
    return [Pv, mp.mpc(0), 4 - 3 * wv, mp.mpc(-2)], Pv


def transport(cur, cfun, t0, t1, depth=0):
    """Advance the 3 tracked roots from parameter t0 to t1 (bisecting on
    trouble): solve the cubic at t1, match by nearest neighbour."""
    if depth > 45:
        raise RuntimeError("bisection depth exceeded")
    cf, lead = cfun(t1)
    if abs(lead) < mp.mpf("1e-14"):
        raise RuntimeError("path passes through the wall")
    new = solve_cubic(cf)
    sep = min(abs(cur[i] - cur[j]) for i in range(3) for j in range(i + 1, 3))
    best = min(PERMS3,
               key=lambda sg: max(abs(new[sg[i]] - cur[i]) for i in range(3)))
    if max(abs(new[best[i]] - cur[i]) for i in range(3)) > 0.35 * sep:
        tm = (t0 + t1) / 2
        cur = transport(cur, cfun, t0, tm, depth + 1)
        return transport(cur, cfun, tm, t1, depth + 1)
    return [new[best[i]] for i in range(3)]


def match_perm(base, cur):
    perm = []
    for cv in cur:
        order = sorted(range(3), key=lambda jj: abs(base[jj] - cv))
        d0_, d1_ = abs(base[order[0]] - cv), abs(base[order[1]] - cv)
        assert d0_ < mp.mpf("1e-12") and d1_ > 1e6 * d0_ + mp.mpf("1e-6"), \
            f"ambiguous sheet matching: d0={mp.nstr(d0_, 5)}"
        perm.append(order[0])
    assert sorted(perm) == [0, 1, 2]
    return tuple(perm)


def track(svals, s_to_coeffs):
    """Track the 3 roots along the polyline of parameter values `svals`;
    return the permutation relative to the starting fiber."""
    base = solve_cubic(s_to_coeffs(svals[0])[0])
    cur = list(base)
    for k in range(len(svals) - 1):
        s0, s1 = svals[k], svals[k + 1]
        cur = transport(cur, lambda t, s0=s0, s1=s1:
                        s_to_coeffs(s0 + (s1 - s0) * t),
                        mp.mpf(0), mp.mpf(1))
    return match_perm(base, cur)


def seg_dist(pt, z0, z1):
    d = z1 - z0
    t = mp.re((pt - z0) * mp.conj(d)) / abs(d) ** 2
    t = min(mp.mpf(1), max(mp.mpf(0), t))
    return abs(pt - (z0 + t * d))


def circle(center, radius, n, start=None):
    ang0 = mp.arg(start - center) if start is not None else mp.mpf(0)
    return [center + radius * mp.exp(mp.mpc(0, 1) *
                                     (ang0 + 2 * mp.pi * mp.mpf(k) / n))
            for k in range(n + 1)]


def loop(basept, center, radius, others, n=120):
    """Basepoint -> best circle start -> CCW circle -> back."""
    best = None
    for k in range(72):
        cand = center + radius * mp.exp(mp.mpc(0, 1) * 2 * mp.pi *
                                        mp.mpf(k) / 72)
        dmin = min([seg_dist(o, basept, cand) for o in others],
                   default=mp.mpf(1))
        if best is None or dmin > best[0]:
            best = (dmin, cand)
    start = best[1]
    approach = [basept + (start - basept) * mp.mpf(k) / 8 for k in range(9)]
    return approach + circle(center, radius, n, start=start)[1:] \
        + approach[::-1][1:]


def cyc(perm):
    seen, out = set(), []
    for i in range(3):
        if i in seen:
            continue
        cl, j = [i], perm[i]
        seen.add(i)
        while j != i:
            cl.append(j)
            seen.add(j)
            j = perm[j]
        if len(cl) > 1:
            out.append("(" + " ".join(str(k + 1) for k in cl) + ")")
    return "".join(out) or "id"


def order_of(perm):
    n, cur = 1, perm
    while cur != (0, 1, 2):
        cur = tuple(perm[cur[i]] for i in range(3))
        n += 1
    return n


def compose(p1_, p2_):
    """(p1 after p2)[i] = p1[p2[i]]  (loop p2 traversed first)."""
    return tuple(p1_[p2_[i]] for i in range(3))


def generated(gens):
    elems = {(0, 1, 2)} | set(gens)
    while True:
        new = {compose(g, h) for g in elems for h in elems} - elems
        if not new:
            return elems
        elems |= new


# --- (a), (b): a generic complex line in the (u,w)-plane --------------------
lu, lw = sp.Rational(-1, 2) + s, sp.Rational(1, 3) + s * (2 + sp.I)
P2_line = sp.Poly(sp.expand(P2.subs({u: lu, w: lw})), s)
D0_line = sp.Poly(sp.expand(D0r.subs({u: lu, w: lw})), s)
assert P2_line.degree() == 3 and D0_line.degree() == 1
assert sp.degree(sp.gcd(P2_line, P2_line.diff(s)), s) == 0   # square-free
assert sp.degree(sp.gcd(P2_line, D0_line), s) == 0           # disjoint loci
assert P2_line.eval(0) != 0
wall_s = [mp.mpc(complex(rr)) for rr in sp.nroots(P2_line, n=20)]
d0_s = [mp.mpc(complex(rr)) for rr in sp.nroots(D0_line, n=20)]
sing = wall_s + d0_s
print("  generic line (u,w) = (-1/2, 1/3) + s*(1, 2+i):")
print(f"    wall points {{P2=0}}: s = "
      f"{', '.join(mp.nstr(rr, 6) for rr in wall_s)}")
print(f"    D0-line point:       s = {mp.nstr(d0_s[0], 6)}")


def line_coeffs(sv):
    return coeffs_inv(mp.mpc(-0.5) + sv, mp.mpf(1) / 3 + sv * mp.mpc(2, 1))


print("\n  (a) small loops around generic wall points -> transpositions:")
wall_perms = []
for ctr in wall_s:
    others = [o for o in sing if abs(o - ctr) > 1e-9]
    rad = mp.mpf("0.3") * min(abs(o - ctr) for o in others)
    perm = track(loop(mp.mpc(0), ctr, rad, others), line_coeffs)
    wall_perms.append(perm)
    print(f"    s* = {mp.nstr(ctr, 6):<24} perm = {cyc(perm)}")
    assert order_of(perm) == 2, "wall loop is not a transposition!"
assert len(set(wall_perms)) >= 2 and len(generated(wall_perms)) == 6
print("    => every wall meridian is a transposition; together they")
print("       generate the full S3 (group order 6).")

print("\n  (b) small loop around the D0-line point (off the wall) -> id:")
ctr = d0_s[0]
others = [o for o in sing if abs(o - ctr) > 1e-9]
rad = mp.mpf("0.3") * min(abs(o - ctr) for o in others)
perm_d0 = track(loop(mp.mpc(0), ctr, rad, others), line_coeffs)
print(f"    s* = {mp.nstr(ctr, 6):<24} perm = {cyc(perm_d0)}")
assert perm_d0 == (0, 1, 2), "D0 loop is not trivial!"

# --- (c): loop around the cusp ---------------------------------------------
print("\n  (c) radius-0.05 loop around the CUSP (4/27, 4/3):")
e1c, e2c = sp.Rational(1) + 3 * sp.I / 10, sp.Rational(-7, 10) + 11 * sp.I / 10
cusp_line = sp.Poly(sp.expand(P2.subs({u: CUSP[0] + s * e1c,
                                       w: CUSP[1] + s * e2c})), s)
# the line through the cusp meets the wall at s = 0 (multiplicity 2, the
# cusp) and at exactly one far point:
assert cusp_line.coeff_monomial(1) == 0 and cusp_line.coeff_monomial(s) == 0
far = -cusp_line.coeff_monomial(s**2) / cusp_line.coeff_monomial(s**3)
far = mp.mpc(complex(far))
print(f"    line through the cusp, direction (1 + 0.3i, -0.7 + 1.1i);")
print(f"    wall intersections: s = 0 (the cusp, double) and s ="
      f" {mp.nstr(far, 6)}")
assert abs(far) > 1        # the circle |s| = 0.05 encloses only the cusp

E1, E2 = mp.mpc(1, "0.3"), mp.mpc("-0.7", "1.1")


def cusp_coeffs(sv):
    return coeffs_inv(mp.mpf(4) / 27 + sv * E1, mp.mpf(4) / 3 + sv * E2)


circ = circle(mp.mpc(0), mp.mpf("0.05"), 240)
minP2 = min(abs(cusp_coeffs(sv)[1]) for sv in circ)
print(f"    min |P2| on the circle: {mp.nstr(minP2, 3)}  (stays off the"
      f" wall)")
assert minP2 > mp.mpf("1e-3")
perm_cusp = track(circ, cusp_coeffs)
print(f"    MEASURED cusp-loop permutation: {cyc(perm_cusp)},"
      f" order {order_of(perm_cusp)}")
# honest determination first (2026-07-21 run: 3-cycle); pinned as regression
assert order_of(perm_cusp) == 3, \
    "cusp loop is not a 3-cycle -- investigate before trusting the docs!"
print("    => a 3-CYCLE (order 3): the image of sigma1*sigma2 in S3, a")
print("       COXETER ELEMENT of W(A2) (order = Coxeter number h = 3).")
print("       NOT the image of sigma1*sigma2*sigma1 (which would be a")
print("       transposition, order 2).  Local model: near the cusp the")
print("       eliminant degenerates to P2 xi^3 = 2 and P2 winds TWICE")
print("       around a cusp loop (its Hessian part is a perfect square),")
print("       so xi ~ (2/P2)^(1/3) is permuted cyclically.")

# decomposition of the cusp loop after perturbing the line off the cusp:
print("\n      perturbed line through (4/27, 4/3 + 1/100), same direction:")
pw0 = CUSP[1] + sp.Rational(1, 100)
pert_line = sp.Poly(sp.expand(P2.subs({u: CUSP[0] + s * e1c,
                                       w: pw0 + s * e2c})), s)
pert_D0 = sp.Poly(sp.expand(D0r.subs({u: CUSP[0] + s * e1c,
                                      w: pw0 + s * e2c})), s)
proots = sorted((mp.mpc(complex(rr)) for rr in sp.nroots(pert_line, n=20)),
                key=abs)
s1, s2 = proots[0], proots[1]
sD = mp.mpc(complex(sp.nroots(pert_D0, n=20)[0]))
assert abs(s1) < 0.01 and abs(s2) < 0.01 and abs(proots[2]) > 1
assert abs(sD) < 0.01      # the D0 point sits between the two wall points


def pert_coeffs(sv):
    return coeffs_inv(mp.mpf(4) / 27 + sv * E1,
                      mp.mpf(4) / 3 + mp.mpf("0.01") + sv * E2)


base_s = mp.mpc("0.05")
mer_rad = mp.mpf("0.25") * min(abs(s1 - s2), abs(s1 - sD), abs(s2 - sD))
pm1 = track(loop(base_s, s1, mer_rad, [s2, sD], n=100), pert_coeffs)
pm2 = track(loop(base_s, s2, mer_rad, [s1, sD], n=100), pert_coeffs)
pmD = track(loop(base_s, sD, mer_rad, [s1, s2], n=100), pert_coeffs)
big = track(loop(base_s, mp.mpc(0), mp.mpf("0.05"), [], n=240), pert_coeffs)
print(f"      the cusp splits into wall points s1 = {mp.nstr(s1, 5)},")
print(f"      s2 = {mp.nstr(s2, 5)}, with the D0 point sD = {mp.nstr(sD, 5)}"
      f" between them;")
print(f"      meridian(s1) = {cyc(pm1)},  meridian(s2) = {cyc(pm2)},"
      f"  meridian(sD) = {cyc(pmD)},")
print(f"      radius-0.05 loop around all three = {cyc(big)}")
assert order_of(pm1) == 2 and order_of(pm2) == 2 and pm1 != pm2
assert pmD == (0, 1, 2)
assert big in (compose(pm1, pm2), compose(pm2, pm1))
assert order_of(big) == 3
print("      => cusp loop = product of two DISTINCT wall transpositions")
print("         (the D0 meridian in between is invisible), hence a 3-cycle.")

# --- universal-family cross-check ------------------------------------------
print("\n  cross-check in the universal family x^3 + Q x + R:")


def univ_generic(sv):
    return [mp.mpf(1), mp.mpc(0), sv * E1, sv * E2], mp.mpf(1)


perm_u = track(circle(mp.mpc(0), mp.mpf("0.05"), 240), univ_generic)
print(f"    radius-0.05 loop around the discriminant cusp (Q,R) = s*(dQ,dR):"
      f" {cyc(perm_u)}, order {order_of(perm_u)}")
assert order_of(perm_u) == 3

EPS = mp.mpf("0.2")


def univ_torus(tv):
    return [mp.mpf(1), mp.mpc(0), EPS**2 * mp.exp(mp.mpc(0, 4) * mp.pi * tv),
            EPS**3 * mp.exp(mp.mpc(0, 6) * mp.pi * tv)], mp.mpf(1)


perm_t = track([mp.mpf(k) / 300 for k in range(301)], univ_torus)
print(f"    torus loop (Q,R) = (eps^2 e^(2it), eps^3 e^(3it)), t: 0 -> 2pi:"
      f" {cyc(perm_t)}")
assert perm_t == (0, 1, 2)


def our_torus(tv):
    # pull the torus loop back through the inverse of the affine map (I4):
    #   u = W - U + 4/27,   w = 4/3 - 3U;  (U,W) = (2 eps^2 e^(2it),
    #   eps^3 e^(3it)) stays off {W^2 = U^3} (W^2 - U^3 = -7 eps^6 e^(6it)).
    Uv = 2 * EPS**2 * mp.exp(mp.mpc(0, 4) * mp.pi * tv)
    Wv = EPS**3 * mp.exp(mp.mpc(0, 6) * mp.pi * tv)
    return coeffs_inv(Wv - Uv + mp.mpf(4) / 27, mp.mpf(4) / 3 - 3 * Uv)


perm_t2 = track([mp.mpf(k) / 300 for k in range(301)], our_torus)
print(f"    same loop pulled back to the invariant plane: {cyc(perm_t2)}")
assert perm_t2 == (0, 1, 2)
print("    => identical monodromy; the torus loop is the FULL TWIST")
print("       (sigma1 sigma2)^3 = Delta^2, the generator of the center of")
print("       B3 -- it maps to id in S3, consistent with (3-cycle)^3 = id.")

# ---------------------------------------------------------------------------
print("\n=== 4. Local system bookkeeping and Euler characteristic ===")
print("  The sheet local system L (rank 3, monodromy = permutation of the")
print("  three eliminant roots) decomposes under S3 = Weyl(A2) as")
print("      L  =  trivial (+) standard,")
print("  by the character identity  chi_perm = chi_triv + chi_std:")
chi_perm = {"id": 3, "transposition": 1, "3-cycle": 0}
chi_triv = {"id": 1, "transposition": 1, "3-cycle": 1}
chi_std = {"id": 2, "transposition": 0, "3-cycle": -1}
for cl in chi_perm:
    assert chi_perm[cl] == chi_triv[cl] + chi_std[cl]
    print(f"      {cl:<14} {chi_perm[cl]:>2} = {chi_triv[cl]} +"
          f" {chi_std[cl]:>2}")

print("\n  TRIVIAL SUMMAND = trace observables (invariant-plane version of")
print("  scripts/trace_pushforward.py): from E = P2 xi^3 + (4-3w) xi - 2,")
e1s, e2s, e3s = sp.Integer(0), qt / P2, 2 / P2
print(f"    e1 = 0,   e2 = (4-3w)/P2,   e3 = 2/P2,")
# consistency with the 3D traces e2 = q/p, e3 = -r/p under xi = X/c:
assert sp.simplify((q / p / c**2).subs({a: u / c**2, b: w / c}) - e2s) == 0
assert sp.simplify((-r / p / c**3).subs({a: u / c**2, b: w / c}) - e3s) == 0
S = {0: sp.Integer(3), 1: e1s}
S[2] = e1s * S[1] - 2 * e2s
S[3] = e1s * S[2] - e2s * S[1] + 3 * e3s
for k in range(4, 7):
    S[k] = e1s * S[k - 1] - e2s * S[k - 2] + e3s * S[k - 3]
for k in range(2, 7):
    den = sp.fraction(sp.cancel(sp.together(S[k])))[1]
    fac = sp.factor_list(den)[1]
    assert all(sp.simplify(f - P2) == 0 for f, _ in fac)
print("    power sums S_2..S_6 over the sheets: RATIONAL, poles only on")
print("    {P2 = 0}  (verified).  All non-rational sheet data lives in the")
print("    2-dim STANDARD summand -- the reflection representation of")
print("    Weyl(A2), i.e. the vanishing cohomology of the A2 singularity.")

print("\n  EULER CHARACTERISTIC of the wall complement:")
# the wall is irreducible ...
assert len(sp.factor_list(P2)[1]) == 1 and sp.factor_list(P2)[1][0][1] == 1
# ... rational with the known parametrization ...
u_m = sp.Rational(4, 27) - 3 * (m - 3) ** 2 / m**3
w_m = sp.Rational(4, 3) - 3 * (m - 3) ** 2 / m**2
assert sp.simplify(P2.subs({u: u_m, w: w_m})) == 0
# ... and the parametrization is INJECTIVE on P^1 \ {m=0}: two parameters
# collide only where both difference-numerators vanish:
cu_ = sp.cancel(sp.numer(sp.together(u_m.subs(m, m1) - u_m.subs(m, m2)))
                / (m1 - m2))
cw_ = sp.cancel(sp.numer(sp.together(w_m.subs(m, m1) - w_m.subs(m, m2)))
                / (m1 - m2))
res = sp.factor(sp.resultant(sp.Poly(cu_, m2), sp.Poly(cw_, m2), m2))
assert sp.expand(res - 2187 * m1**2 * (m1 - 3) ** 2) == 0
# so collisions require m1 in {0, 3}; at the cusp parameter m1 = 3 the only
# common solution is m2 = 3 (unibranch), m = 0 is the point at infinity:
assert sp.solve([cu_.subs(m1, 3), cw_.subs(m1, 3)], m2, dict=True) == \
    [{m2: 3}]
assert sp.simplify(u_m.subs(m, 3) - CUSP[0]) == 0
print("    - the affine wall is an irreducible rational curve; its")
print("      normalization P^1 -> closure is a BIJECTION (resultant of the")
print("      collision equations = 2187 m1^2 (m1-3)^2: parameters collide")
print("      only at m = 0 [the single place at infinity] and m = 3 [the")
print("      cusp, unibranch: only m2 = 3 maps there]).")
print("    - hence the affine curve is homeomorphic to P^1 minus one point")
print("      = C, so  chi(affine wall) = chi(C) = 1.")
print("    - excision/additivity:  chi(C^2 \\ wall) = chi(C^2) -")
print("      chi(wall) = 1 - 1 = 0.")
print("    - cross-check [Arnold 1969]: H*(B3; Z) = (Z, Z, 0, ...), so the")
print("      K(B3,1) wall complement has chi = 1 - 1 + 0 = 0.  Consistent.")

print(f"\nAll assertions passed.  Total wall time: {time.time() - T0:.1f} s")
