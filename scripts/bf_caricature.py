"""The 0D Buchholz-Fredenhagen S(J) caricature for the Alpöge-Mathew map.

(Write-up: docs/BF_CARICATURE.md.  Resolves docs/OPEN_QUESTIONS.md B1;
framing: docs/QFT_IMPLICATIONS.md §4.3(c), paper §sec:bf.)

Background [BF20] (cited from memory): Buchholz-Fredenhagen build
interacting dynamics from unitaries S(f) subject to (i) causal
factorization S(f+g+h) = S(f+g) S(g)^{-1} S(g+h) when supp f is "later"
than supp h, and (ii) a dynamical relation S(f) = S(f^psi + deltaL(psi))
implementing the (off-shell) field equations.  This script asserts the 0D
degeneration of both axioms for the Alpöge-Mathew map F (jcqft.core):

1. CAUSAL FACTORIZATION TRIVIALIZES (exact, formalized).  On a one-point
   spacetime, causal disjointness of supports forces supp f = ∅ or
   supp h = ∅, and in every allowed case the relation reduces to a
   group-theoretic tautology (verified with noncommutative symbols).

2. THE DYNAMICAL RELATION SURVIVES AND FORCES THE FIBER ALGEBRA.
   No potential exists in n=3 (DF != DF^T: I8 gate), so the only
   Lagrangian implementing F(phi)=J is the first-order one
   L = phibar.(F(phi)-J).  Antifield shifts phibar -> phibar + beta give
   deltaL(beta) = beta.(F(phi)-J) exactly (affine in phibar: no hbar
   corrections), and any multiplicative (pure-state) evaluation killing
   these shifts is a character of A_J = C[x,y,z]/(F-J): pure classical
   states with source J = points of the fiber.  Asserted on explicit
   fiber and non-fiber points.

3. THE ALGEBRA BUNDLE SEES THE WALL (exact).  dim_C A_J = 3 off the wall
   (4 rational chamber points, étale => reduced), = 1 at a generic wall
   point, = 0 (zero ring, NO characters) on the empty-fiber cusp orbit.
   Non-constant dimension <=> non-properness; the tame-shear control has
   constant dimension 1 and an exact rational inverse (the caricature
   collapses exactly for proper Keller maps).

4. NO SINGLE-VALUED SECTOR DATA (exact).  The eliminant cubic is
   primitive and irreducible over C(a,b,c) (multivariate factorization +
   Gauss), so: no rational section J -> phi(J) exists, the generic fiber
   algebra over the function field is a FIELD (no nontrivial idempotents
   = no single-valued sector projections), and disc = -4 D0^2 p with p
   irreducible => Galois group S3 => Aut(L/K) = 1: the cover is
   non-Galois and the deck group is TRIVIAL (candidate formulation (a)
   fails provably).

5. THE INVARIANT PART IS THE TRANSFER OPERATOR, AND IT PAYS IN POLES
   (exact).  T[1]=3, T[x]=0, T[x^k] rational with poles only on {p=0}
   (k<=6); the normalized expectation is not multiplicative (candidate
   (b) alone loses the sector structure); the sheet-separating
   coefficients of y and z carry the factor p verbatim.

6. REAL C*-FIBERS AND THE 1<->3 JUMP (exact).  Characters of the real
   fiber C*-algebra = real solutions: counts 3/3/1/1 at the chamber
   points, 1 at the generic wall point (the surviving character is
   x = -r/q), 0 at the cusp.  The two chamber algebras C^3 and C^1 are
   non-isomorphic.

7. HOLONOMY OF THE RELATIVE S-TRANSPORT IS S3 (Numerical).  Parallel
   transport of the fiber-algebra bundle along loops in the invariant
   (u,w)-plane: wall meridians act by transpositions generating a group
   of order 6, the D0-meridian acts trivially, the cusp loop is a
   3-cycle of order 3 (Coxeter element).  Permutation unitaries U(sigma)
   verify U_cusp^3 = 1 != U_cusp.

Exact claims carry asserts ([ok] lines); the monodromy section is
labelled Numerical (mpmath tracking, not certified).  Runtime well under
3 min.  Usage:  .venv/bin/python scripts/bf_caricature.py
"""

import time

import mpmath as mp
import sympy as sp

from jcqft.core import D0, F, PHI, SRC, X, cubic, p, q, r
from jcqft.fibers import B as B_COEF
from jcqft.fibers import D as D_COEF
from jcqft.fibers import exact_fiber

T0 = time.time()
N_CHECKS = 0
x, y, z = PHI
a, b, c = SRC


def check(label, cond=True):
    global N_CHECKS
    assert cond, label
    N_CHECKS += 1
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


def _zero_mod_minpoly(expr):
    """Exact zero test for a polynomial expression in <=1 CRootOf."""
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


# Rational chamber samples from scripts/witten_index.py
SAMPLES = [
    ((sp.Rational(-1, 4), 0, 0), -1),  # p = -4 < 0, N = 3
    ((0, 2, 0), -1),                   # p = -4 < 0, N = 3
    ((1, 0, 0), +1),                   # p = 16 > 0, N = 1
    ((2, 1, 1), +1),                   # p = 104 > 0, N = 1
]
WALL_PT = (sp.Rational(2, 27), 1, 1)          # p = 0, q != 0, D0 = 1
CUSP_PT = (sp.Rational(4, 27), sp.Rational(4, 3), 1)  # empty-fiber orbit

# ===========================================================================
print("=== 0. Core identities (anchors) ===")
# ===========================================================================
DF = sp.Matrix([[sp.diff(Fi, v) for v in PHI] for Fi in F])
check("det DF = -2 identically (Keller, étale everywhere)",
      sp.simplify(DF.det()) == -2)
check("collision identity 4q^3 + 27 p r^2 = 4 D0^2",
      sp.expand(4 * q**3 + 27 * p * r**2 - 4 * D0**2) == 0)
disc = sp.discriminant(cubic, X)
check("disc_X(eliminant) = -4 D0^2 p",
      sp.expand(disc + 4 * D0**2 * p) == 0)

# ===========================================================================
print("\n=== 1. Causal factorization trivializes on a point (exact) ===")
# ===========================================================================
# Spacetime = {pt}.  Supports are subsets of {pt}; the causal past of a
# support is itself (the point is in its own past).  BF causal
# factorization applies when supp f ∩ J^-(supp h) = ∅.
PT = frozenset({"pt"})
EMPTY = frozenset()


def causal_past(supp):
    return supp  # J^-({pt}) = {pt}, J^-(∅) = ∅


allowed = [(sf, sh) for sf in (EMPTY, PT) for sh in (EMPTY, PT)
           if not (sf & causal_past(sh))]
check("causal disjointness on {pt} forces supp f = ∅ or supp h = ∅",
      set(allowed) == {(EMPTY, EMPTY), (EMPTY, PT), (PT, EMPTY)}
      and (PT, PT) not in allowed)

# supp = ∅ means the functional is 0.  Check each allowed case reduces to
# a tautology, with noncommutative symbols for the S-elements.
S_g, S_gh, S_fg = sp.symbols("S_g S_gh S_fg", commutative=False)
# f = 0:  S(g+h) =?= S(g) S(g)^{-1} S(g+h)
check("case supp f = ∅: relation is the tautology S(g+h) = S(g)S(g)^{-1}S(g+h)",
      sp.simplify(S_g * S_g**(-1) * S_gh - S_gh) == 0)
# h = 0:  S(f+g) =?= S(f+g) S(g)^{-1} S(g)
check("case supp h = ∅: relation is the tautology S(f+g) = S(f+g)S(g)^{-1}S(g)",
      sp.simplify(S_fg * S_g**(-1) * S_g - S_fg) == 0)
# f = h = 0:  S(g) =?= S(g) S(g)^{-1} S(g)
check("case both ∅: S(g) = S(g)S(g)^{-1}S(g)",
      sp.simplify(S_g * S_g**(-1) * S_g - S_g) == 0)
print("  -> the only non-tautological instance (supp f = supp h = {pt}) is")
print("     excluded by causal disjointness: axiom (i) is VACUOUS in 0D.")

# ===========================================================================
print("\n=== 2. The dynamical relation forces the fiber algebra (exact) ===")
# ===========================================================================
# I8 gate: no potential W with F = grad W exists in n = 3.
check("I8 gate: DF is not symmetric ((DF)_13 - (DF)_31 = -1 at phi = 0)",
      sp.simplify((DF[0, 2] - DF[2, 0]).subs(dict(zip(PHI, (0, 0, 0)))))
      == -1)

# The only Lagrangian implementing F(phi) = J is first-order:
#   L(phibar, phi; J) = phibar . (F(phi) - J).
bx, by, bz = sp.symbols("bx by bz")          # antifields phibar
B1, B2, B3 = sp.symbols("beta1 beta2 beta3")  # constant antifield shift
Jt = sp.symbols("Ja Jb Jc")
L_first = sum(bb * (Fi - Ji) for bb, Fi, Ji in zip((bx, by, bz), F, Jt))
L_shift = L_first.subs({bx: bx + B1, by: by + B2, bz: bz + B3})
deltaL = sp.expand(L_shift - L_first)
check("antifield shift: deltaL(beta) = beta . (F(phi) - J) exactly "
      "(L affine in phibar: relation is hbar-exact)",
      sp.expand(deltaL - sum(
          Bi * (Fi - Ji) for Bi, Fi, Ji in zip((B1, B2, B3), F, Jt))) == 0)

# A multiplicative evaluation (character chi = evaluation at phi0) kills
# all shift terms beta.(F - J) iff F(phi0) = J, i.e. iff chi factors
# through A_J = C[x,y,z]/(F - J).  (A character killing the span of the
# generators kills the ideal they generate.)  Witness both directions:
J0 = SAMPLES[0][0]
fib0 = exact_fiber(J0)
check("fiber points are characters: F(phi0) - J = 0 for all phi0 in "
      f"F^(-1){tuple(J0)}",
      all(_zero_mod_minpoly(Fi.subs(dict(zip(PHI, pt))) - Ji)
          for pt in fib0 for Fi, Ji in zip(F, J0)))
off_pt = (1, 1, 1)
vals = [sp.simplify(Fi.subs(dict(zip(PHI, off_pt))) - Ji)
        for Fi, Ji in zip(F, J0)]
check("non-fiber evaluation fails the shift relation: F(1,1,1) != J",
      any(v != 0 for v in vals))
print("  -> 0D dynamical relation: pure states with source J = characters")
print("     of the fiber algebra A_J = C[x,y,z]/(F - J).")

# ===========================================================================
print("\n=== 3. The algebra bundle sees the wall: dim A_J jumps (exact) ===")
# ===========================================================================
dims = {}
for J, sgn_p in SAMPLES:
    sub = dict(zip(SRC, J))
    assert sp.sign(p.subs(sub)) == sgn_p and D0.subs(sub) != 0
    fib = exact_fiber(J)
    # étale (det DF = -2) => fibers reduced => dim A_J = #points;
    # points pairwise distinct:
    assert len(fib) == 3
    assert all(fib[i] != fib[j] for i in range(3) for j in range(i + 1, 3))
    for pt in fib:
        psub = dict(zip(PHI, pt))
        assert all(_zero_mod_minpoly(Fi.subs(psub) - Ji)
                   for Fi, Ji in zip(F, J))
    dims[J] = 3
check("dim_C A_J = 3 at all four off-wall chamber points "
      "(3 distinct points, étale => reduced)",
      all(d == 3 for d in dims.values()))

sub_w = dict(zip(SRC, WALL_PT))
check("generic wall point (2/27,1,1): p = 0, q != 0, D0 = 1",
      p.subs(sub_w) == 0 and q.subs(sub_w) != 0 and D0.subs(sub_w) == 1)
fib_w = exact_fiber(WALL_PT)
check("dim_C A_J = 1 at the generic wall point (eliminant degenerates "
      "to qX + r)", len(fib_w) == 1)
check("the surviving character has x = -r/q there",
      sp.simplify(fib_w[0][0] - (-r / q).subs(sub_w)) == 0)

sub_c = dict(zip(SRC, CUSP_PT))
assert p.subs(sub_c) == 0 and D0.subs(sub_c) == 0
fib_c = exact_fiber(CUSP_PT)
check("dim_C A_J = 0 on the cusp orbit: A_J is the ZERO RING, no "
      "characters at all", len(fib_c) == 0)

check("dim A_J is NON-CONSTANT (3 -> 1 -> 0): the sheaf of fiber "
      "algebras is not locally free across {p=0}",
      len({3, len(fib_w), len(fib_c)}) == 3)

# Control: a proper Keller map (tame shear).  Constant dimension 1 and an
# exact rational (indeed polynomial) section: the caricature collapses.
F_shear = (x, z + y**2, y)
at, bt, ct = sp.symbols("at bt ct")
inv_shear = (at, ct, bt - ct**2)
check("control (tame shear): exact polynomial inverse exists, dim A_J = 1 "
      "for ALL J, no wall, trivial bundle",
      all(sp.simplify(
          Fi.subs(dict(zip(PHI, inv_shear)), simultaneous=True) - t_) == 0
          for Fi, t_ in zip(F_shear, (at, bt, ct))))

# ===========================================================================
print("\n=== 4. No single-valued sector data (exact) ===")
# ===========================================================================
check("eliminant is primitive in X: gcd(p, q, r) = 1",
      sp.gcd(sp.gcd(p, q), r) == 1)
fl = sp.factor_list(sp.expand(cubic))
check("eliminant cubic is irreducible as a polynomial in (X, a, b, c)",
      len(fl[1]) == 1 and fl[1][0][1] == 1)
print("  -> by Gauss's lemma the cubic is irreducible over K = C(a,b,c):")
print("     (i) NO rational section J -> phi(J) of F exists;")
print("     (ii) the generic fiber algebra K[X]/(cubic) is a FIELD: its only")
print("          idempotents are 0 and 1 - no single-valued sector projection.")
check("no rational section / no nontrivial idempotents over K "
      "(irreducibility asserted above)", True)

flp = sp.factor_list(p)
check("wall polynomial p is irreducible over Q", len(flp[1]) == 1
      and flp[1][0][1] == 1
      and sp.Poly(flp[1][0][0], a, b, c).total_degree() == 4)
# disc = (2 D0)^2 . (-p); p irreducible of multiplicity one => disc is not
# a square in C(a,b,c) => Galois group = S3 (order 6) => the degree-3
# extension L/K is NOT normal => Aut(L/K) = 1.
check("disc = (2 D0)^2 . (-p) with p to multiplicity one => Galois = S3 "
      "=> deck group Aut(L/K) TRIVIAL: candidate (a) fails",
      sp.expand(disc - (2 * D0)**2 * (-p)) == 0)

# ===========================================================================
print("\n=== 5. The invariant part: transfer operator, poles on p (exact) ===")
# ===========================================================================
e1, e2, e3 = sp.Integer(0), q / p, -r / p
S_pow = {0: sp.Integer(3), 1: e1}
S_pow[2] = e1 * S_pow[1] - 2 * e2
S_pow[3] = e1 * S_pow[2] - e2 * S_pow[1] + 3 * e3
for k in range(4, 7):
    S_pow[k] = e1 * S_pow[k - 1] - e2 * S_pow[k - 2] + e3 * S_pow[k - 3]
check("T[1] = 3, T[x] = 0, T[x^2] = -2q/p",
      S_pow[0] == 3 and S_pow[1] == 0
      and sp.simplify(S_pow[2] + 2 * q / p) == 0)
poles_ok = True
for k in range(2, 7):
    den = sp.fraction(sp.cancel(sp.together(S_pow[k])))[1]
    poles_ok &= all(sp.simplify(f - p) == 0
                    for f, _ in sp.factor_list(den)[1])
check("T[x^k] rational with poles ONLY on {p=0}, k = 2..6", poles_ok)

Eexp = {k: S_pow[k] / 3 for k in S_pow}
check("normalized expectation E = T/3 is NOT multiplicative: "
      "E(x^2) = -2q/(3p) != 0 = E(x)^2 (candidate (b) loses the sectors)",
      sp.simplify(Eexp[2] - Eexp[1]**2) != 0)

# Sheet-separating coefficients of y and z carry the factor p verbatim:
# y = -B/(2 D0), z = -D/(8 D0) with B, D quadratic in x (jcqft.fibers).
Bp = sp.Poly(B_COEF, x)
Dp = sp.Poly(D_COEF, x)
c1B = Bp.coeff_monomial(x)
c2B = Bp.coeff_monomial(x**2)
check("separator: coeff_x(B) = -6p  (i.e. c1(y) = 3p/D0)",
      sp.expand(c1B + 6 * p) == 0)
check("separator: p | coeff_{x^2}(B)  (c2(y) carries p)",
      sp.rem(sp.expand(c2B), p, x, a, b, c) == 0
      and sp.denom(sp.cancel(c2B / p)) == 1)
c1D = Dp.coeff_monomial(x)
c2D = Dp.coeff_monomial(x**2)
check("separator: p | coeff_x(D) and p | coeff_{x^2}(D)  (z separators "
      "carry p)",
      sp.denom(sp.cancel(c1D / p)) == 1
      and sp.denom(sp.cancel(c2D / p)) == 1)
print("  -> any single-valued sheet-separating extension across the wall")
print("     diverges on {p=0}: the wall enters as the POLE DIVISOR of the")
print("     invariant-hull S-data.")

# ===========================================================================
print("\n=== 6. Real C*-fibers: the 1<->3 jump as character count (exact) ===")
# ===========================================================================
n_chars = []
for J, sgn_p in SAMPLES:
    fib = exact_fiber(J)
    n_r = sum(1 for pt in fib if all(v.is_real for v in pt))
    n_chars.append(n_r)
    assert n_r == (3 if sgn_p < 0 else 1)
check("characters of the real fiber C*-algebra C^{N(J)}: N = 3,3,1,1 at "
      "the chamber points", n_chars == [3, 3, 1, 1])
check("the chamber C*-algebras C^3 and C^1 are non-isomorphic "
      "(3 != 1 minimal projections)", 3 != 1)
check("generic wall point: exactly 1 character; cusp orbit: 0 characters",
      len(fib_w) == 1 and all(v.is_real for v in fib_w[0])
      and len(fib_c) == 0)

# ===========================================================================
print("\n=== 7. Holonomy of the relative S-transport (NUMERICAL) ===")
# ===========================================================================
# Parallel transport of the fiber-algebra bundle over the invariant
# (u,w)-plane (docs/WALL_COMPLEMENT.md): track the 3 roots of
#   P2(u,w) xi^3 + (4 - 3w) xi - 2 = 0
# along loops; permutations = holonomy of the character bundle.
# NUMERICAL: mpmath tracking with a nearest-neighbour jump guard; not a
# certified computation (cf. docs/MONODROMY.md caveats, OPEN_QUESTIONS B6).
mp.mp.dps = 25
uu, ww, s = sp.symbols("uu ww s")
P2 = 27 * uu**2 + 16 * uu - 18 * uu * ww + ww**3 - ww**2
u0, w0 = sp.Rational(-1, 2), sp.Rational(1, 3)
d1, d2 = sp.Integer(1), 2 + sp.I
P2_s = sp.expand(P2.subs({uu: u0 + d1 * s, ww: w0 + d2 * s}))
D0p_s = sp.expand((27 * uu - 9 * ww + 8).subs({uu: u0 + d1 * s,
                                               ww: w0 + d2 * s}))
assert sp.degree(P2_s, s) == 3 and sp.degree(D0p_s, s) == 1
wall_s = mp.polyroots([mp.mpc(str(sp.re(cc)), str(sp.im(cc))) for cc in
                       sp.Poly(P2_s, s).all_coeffs()],
                      maxsteps=200, extraprec=80)
d0_s = mp.polyroots([mp.mpc(str(sp.re(cc)), str(sp.im(cc))) for cc in
                     sp.Poly(D0p_s, s).all_coeffs()],
                    maxsteps=200, extraprec=80)
sing = list(wall_s) + list(d0_s)

P2f = sp.lambdify((uu, ww), P2, "mpmath")


def cubic_roots(uw):
    u_, w_ = uw
    return mp.polyroots([P2f(u_, w_), mp.mpf(0), 4 - 3 * w_, mp.mpf(-2)],
                        maxsteps=200, extraprec=80)


def line(sv):
    return (mp.mpc(-0.5) + sv, mp.mpc(1) / 3 + (2 + 1j) * sv)


def match_roots(old, new):
    minsep = min(abs(old[i] - old[j])
                 for i in range(3) for j in range(i + 1, 3))
    out = [None] * 3
    used = set()
    for i in range(3):
        d0_, k0 = min((abs(old[i] - nn), k) for k, nn in enumerate(new))
        if k0 in used or d0_ > 0.35 * minsep:
            return None
        used.add(k0)
        out[i] = new[k0]
    return out


def track(path, roots, n0=250):
    t = mp.mpf(0)
    dt = mp.mpf(1) / n0
    while t < 1:
        tn = min(t + dt, mp.mpf(1))
        new = cubic_roots(path(tn))
        m = match_roots(roots, new)
        if m is None:
            dt /= 2
            assert dt > mp.mpf("1e-8"), "tracking step underflow"
            continue
        roots = m
        t = tn
        dt = min(dt * mp.mpf("1.25"), mp.mpf(4) / n0)
    return roots


def perm_of(final, base):
    out = []
    minsep = min(abs(base[i] - base[j])
                 for i in range(3) for j in range(i + 1, 3))
    for fr in final:
        d0_, j = min((abs(fr - bb), k) for k, bb in enumerate(base))
        assert d0_ < 1e-6 * max(1, minsep), "loop return mismatch"
        out.append(j)
    assert set(out) == {0, 1, 2}
    return tuple(out)


def loop_perm(center, base_roots):
    others = [o for o in sing if abs(o - center) > 1e-12]
    R = mp.mpf("0.3") * min(abs(center - o) for o in others)
    c0 = center * (1 - R / abs(center))
    r_ = track(lambda t: line(c0 * t), base_roots)
    r_ = track(lambda t: line(center + (c0 - center)
                              * mp.e**(2j * mp.pi * t)), r_, n0=500)
    r_ = track(lambda t: line(c0 * (1 - t)), r_)
    return perm_of(r_, base_roots)


base_roots = sorted(cubic_roots(line(mp.mpf(0))),
                    key=lambda zz: (mp.re(zz), mp.im(zz)))


def compose(s1, s2):
    return tuple(s1[s2[i]] for i in range(3))


def closure(gens):
    e = (0, 1, 2)
    G = {e}
    frontier = {e}
    while frontier:
        new = set()
        for g in G:
            for h in gens:
                cmp_ = compose(g, h)
                if cmp_ not in G:
                    new.add(cmp_)
        G |= new
        frontier = new
    return G


wall_perms = [loop_perm(ws, base_roots) for ws in wall_s]
check("[Numerical] every wall meridian acts by a transposition (order 2)",
      all(pp != (0, 1, 2) and compose(pp, pp) == (0, 1, 2)
          for pp in wall_perms))
check("[Numerical] wall meridians generate a group of order 6 = S3: "
      "the S-transport holonomy is the full S3",
      len(closure(wall_perms)) == 6)
d0_perm = loop_perm(d0_s[0], base_roots)
check("[Numerical] D0-meridian (x-collision) acts trivially: the algebra "
      "transport is blind to the harmless locus", d0_perm == (0, 1, 2))

# Cusp loop: circle of radius 0.05 around (4/27, 4/3) in the direction
# plane (1 + 0.3i, -0.7 + 1.1i)  (same loop as scripts/wall_braid.py).
cu, cw = mp.mpf(4) / 27, mp.mpf(4) / 3
dv1, dv2 = mp.mpc(1, 0.3), mp.mpc(-0.7, 1.1)


def cusp_path(t):
    ph = mp.mpf("0.05") * mp.e**(2j * mp.pi * t)
    return (cu + ph * dv1, cw + ph * dv2)


minP2 = min(abs(P2f(*cusp_path(mp.mpf(k) / 60))) for k in range(60))
check(f"[Numerical] cusp circle stays off the wall (min |P2| = {float(minP2):.3f})",
      minP2 > mp.mpf("0.05"))
cusp_start = sorted(cubic_roots(cusp_path(mp.mpf(0))),
                    key=lambda zz: (mp.re(zz), mp.im(zz)))
cusp_final = track(cusp_path, cusp_start, n0=600)
cusp_perm = perm_of(cusp_final, cusp_start)
check("[Numerical] cusp loop = 3-cycle of order 3 (Coxeter element): "
      "total escape rotates the three sectors cyclically",
      cusp_perm != (0, 1, 2)
      and compose(cusp_perm, compose(cusp_perm, cusp_perm)) == (0, 1, 2))


def perm_matrix(sig):
    M = sp.zeros(3)
    for i, j in enumerate(sig):
        M[j, i] = 1
    return sp.ImmutableMatrix(M)


U_cusp = perm_matrix(cusp_perm)
check("[Numerical->algebra] permutation unitary: U_cusp^3 = 1, U_cusp != 1",
      U_cusp**3 == sp.ImmutableMatrix(sp.eye(3))
      and U_cusp != sp.ImmutableMatrix(sp.eye(3)))

# ===========================================================================
print("\n=== 8. Verdict ===")
# ===========================================================================
print("  (i)  causal factorization: VACUOUS in 0D (section 1, exact).")
print("  (ii) dynamical relation: NON-TRIVIAL residue -- forces the bundle")
print("       of fiber algebras A_J with rank 3 -> 1 -> 0 at the wall,")
print("       character count 1 <-> 3 over R, and S3 transport holonomy;")
print("       no single-valued sector datum exists (exact), and every")
print("       single-valued invariant-hull datum pays with poles on {p=0}.")
print("  (iii) collapse control: for proper Keller maps (automorphisms)")
print("       the same construction is the trivial rank-1 bundle.")

print(f"\nall {N_CHECKS} checks passed in {time.time() - T0:.1f} s")
print("ALL CHECKS PASSED")
