"""Twisted cohomology of the wall complement M = C^2 \\ {P2 = 0}: exact
dimensions, jump loci, and the twisted-period / canonical-form connection.

Companion script of docs/TWISTED_PERIODS.md; answers Q3 of
docs/WALL_COMPLEMENT.md section 6 (the concrete follow-up to resolved
question B2 of docs/OPEN_QUESTIONS.md).

METHOD.  M is a K(B3,1) (docs/WALL_COMPLEMENT.md), so twisted cohomology
of M with any local system rho equals group cohomology H^*(B3; rho).
B3 = <x, y | xyx = yxy> is a 1-relator group whose relator is not a proper
power, so by Lyndon's theorem the presentation 2-complex is aspherical and
the Fox free resolution gives the 3-term cochain complex

    0 -> V --d0--> V^2 --d1--> V -> 0,
    d0(v)     = ((rho(x) - 1) v, (rho(y) - 1) v),
    d1(v1,v2) = rho(dr/dx) v1 + rho(dr/dy) v2,

with r = x y x y^-1 x^-1 y^-1 and Fox derivatives
    dr/dx = 1 + xy - xyx y^-1 x^-1,     dr/dy = x - xyx y^-1 - 1,
which the braid relation simplifies to 1 + xy - y and x - yx - 1.
Note chi = h0 - h1 + h2 = 0 automatically (rank-nullity on this shape of
complex), consistent with chi(M) = 0; the CONTENT is in the individual
dimensions and in the agreement of independent methods.

Chain of results (all exact, all asserted):

  1. RANK-1 KUMMER TWIST P2^s (lambda_t: sigma_i -> t, t = e^{2 pi i s}):
     d1 = Delta(t) * (1, -1) with Delta(t) = t^2 - t + 1 the ALEXANDER
     POLYNOMIAL of the trefoil.  Dimensions (h0, h1, h2):
         t = 1        (s in Z):             (1, 1, 0)  = untwisted H^*(M)
         Delta(t) = 0 (s in 1/6+Z, 5/6+Z):  (0, 1, 1)
         all other t:                       (0, 0, 0).
  2. WANG-SEQUENCE CROSS-CHECK.  P2 is affinely equivalent to the
     quasi-homogeneous 27(W^2 - U^3), so P2: M -> C* is the global Milnor
     fibration of the cusp; fiber F has (b0, b1) = (1, 2), monodromy h of
     order 6 with eigenvalues the primitive 6th roots of unity (Brieskorn
     formula for x^2 + y^3).  The Wang sequence reproduces the Fox
     dimensions for every tested t, and the jump determinants agree
     IDENTICALLY as polynomials: det(t h* - 1 | H^1(F)) = Delta(t),
     det(t - 1 | H^0(F)) = t - 1.
  3. S3 LOCAL SYSTEMS (through B3 ->> S3): standard/reflection rep S:
     (0, 1, 1); sign: (0, 0, 0) (= rank-1 case at t = -1); permutation:
     (1, 2, 1) = trivial (+) standard (additivity verified); regular
     C[S3]: (1, 3, 2) = H^*(P3; C) = Betti numbers of the A2 braid
     arrangement complement (Shapiro's lemma cross-check, matching an
     independent Orlik-Solomon/Moebius computation).
  4. REDUCED BURAU rho_t (s1 -> [[-t,1],[0,1]], s2 -> [[1,0],[t,-t]]):
     jump locus EXACTLY {t^3 = 1} (t = 0 excluded):
         generic t:                (0, 0, 0)
         t = 1:                    (0, 1, 1)   [Burau(1) IS the standard
                                                rep -- matrix equality]
         t = primitive cube root:  (1, 2, 1)   [Burau reducible: trivial
                                                sub, quotient = Kummer
                                                character at -t, a
                                                primitive 6th root].
     The gcd of the maximal minors of d1 is t^3 - 1 = (t-1)(t^2+t+1) and
     t^2+t+1 = Delta(-t): the Alexander polynomial appears through the
     standard Burau sign convention t_Burau = -t_Kummer.  Specializations:
     t = 1 recovers case 3 (standard rep); det o rho_t = lambda_{-t}
     recovers case 1; rho_{-1} is the integral SL(2,Z) representation
     (full twist -> -I, sigma_i unipotent), NOT the reflection rep.
  5. DE RHAM SIDE (twisted rational forms, nabla_s = d + s dlog f, with
     f = W^2 - U^3 = P2/27 in the A2 coordinates of WALL_COMPLEMENT (I4)):
     with the Euler field E = 2U d/dU + 3W d/dW and eta_g = g iota_E(
     dU^dW)/f^m (g = U^a W^b of weight k = 2a + 3b),
         nabla_s eta_g = (k + 5 - 6m + 6s) * g dU^dW / f^m,
     so pole-order reduction fails exactly at s = (6m - 5 - k)/6.  On the
     Milnor-ring basis {1, U} of the cusp this picks out s in 1/6 + Z
     (class [dU^dW/f]) and s in 5/6 + Z (class [U dU^dW/f]); dlog f spans
     the s in Z jump (df/f = nabla_s(1/s) otherwise).  Non-basis
     numerators do NOT create jumps: W dU^dW/f^2 = nabla_s(dU/f)/(2(1-s))
     kills the spurious Euler-chain failure at s = 2/3.  At s = 0 the
     would-be canonical form of docs/POSITIVE_GEOMETRY.md is globally
     EXACT: dU^dW/f = d(-eta), eta = (2U dW - 3W dU)/f -- zero periods,
     the cohomological face of the residueless double pole.

Run:  .venv/bin/python scripts/twisted_cohomology.py     (~10 s)
"""

import time
from itertools import combinations

import sympy as sp

T0 = time.time()

t, s = sp.symbols("t s")
U, W, u, w = sp.symbols("U W u w")
I2 = sp.eye(2)

Z6 = sp.Rational(1, 2) + sp.sqrt(3) * sp.I / 2       # primitive 6th root
Z6b = sp.Rational(1, 2) - sp.sqrt(3) * sp.I / 2      # = Z6^5 = Z6^-1
Z3 = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2      # primitive cube root
Z3b = sp.Rational(-1, 2) - sp.sqrt(3) * sp.I / 2     # = Z3^2
DELTA = t**2 - t + 1                                 # trefoil Alexander


def simp0(x):
    return sp.simplify(sp.expand(x)) == 0


def is_zero_mat(M):
    return all(simp0(e) for e in M)


def exact_rank(M):
    return M.rank(iszerofunc=lambda x: sp.simplify(x) == 0)


def fox_maps(R1, R2):
    """The Fox cochain complex 0 -> V -> V^2 -> V -> 0 for B3 = <x,y|r>,
    r = xyx y^-1 x^-1 y^-1, with rho(x) = R1, rho(y) = R2."""
    n = R1.rows
    Id = sp.eye(n)
    assert is_zero_mat(R1 * R2 * R1 - R2 * R1 * R2), "braid relation fails"
    # Fox derivatives straight from the product rule (with inverses) ...
    A1_full = Id + R1 * R2 - R1 * R2 * R1 * R2.inv() * R1.inv()
    A2_full = R1 - R1 * R2 * R1 * R2.inv() - Id
    # ... and the braid-relation-simplified words xyxy^-1x^-1 = y,
    # xyxy^-1 = yx (a genuine internal consistency check):
    A1 = Id + R1 * R2 - R2
    A2 = R1 - R2 * R1 - Id
    assert is_zero_mat(A1_full - A1) and is_zero_mat(A2_full - A2)
    d0 = (R1 - Id).col_join(R2 - Id)            # V -> V^2   (2n x n)
    d1 = A1.row_join(A2)                        # V^2 -> V   (n x 2n)
    assert is_zero_mat(d1 * d0)                 # Fox fundamental identity
    return d0, d1


def coh_dims(R1, R2):
    """(h0, h1, h2) of H^*(B3; rho), exactly."""
    d0, d1 = fox_maps(R1, R2)
    n = R1.rows
    r0, r1 = exact_rank(d0), exact_rank(d1)
    dims = (n - r0, (2 * n - r1) - r0, n - r1)
    assert dims[0] - dims[1] + dims[2] == 0     # chi = 0 (structural)
    return dims


# ---------------------------------------------------------------------------
print("=== 1. Fox calculus set-up; untwisted sanity check ===")
dims_triv = coh_dims(sp.Matrix([[1]]), sp.Matrix([[1]]))
print(f"  trivial local system: (h0, h1, h2) = {dims_triv}")
assert dims_triv == (1, 1, 0)
print("  = (C, C, 0) = H^*(B3; Z) (x) C [Arnold 1969]; chi = 0 matches")
print("  chi(M) = 0 of docs/WALL_COMPLEMENT.md section 4.  (chi = 0 holds")
print("  for EVERY local system below by rank-nullity on the 3-term")
print("  complex -- the real content is the individual dimensions.)")

# ---------------------------------------------------------------------------
print("\n=== 2. Rank-1 Kummer twist P2^s: lambda_t(sigma_1) ="
      " lambda_t(sigma_2) = t ===")
Rt = sp.Matrix([[t]])
d0k, d1k = fox_maps(Rt, Rt)
d0k = d0k.applyfunc(sp.cancel)
d1k = d1k.applyfunc(sp.cancel)
print(f"  d0 = {d0k.T.tolist()[0]}^T,   d1 = {d1k.tolist()[0]}")
assert d0k == sp.Matrix([[t - 1], [t - 1]])
assert d1k == sp.Matrix([[DELTA, -DELTA]])
print("  d1 = Delta(t) * (1, -1)  with  Delta(t) = t^2 - t + 1: the")
print("  ALEXANDER POLYNOMIAL of the trefoil (M is homotopy equivalent to")
print("  the trefoil complement).  First elementary ideal = (Delta):")
g_ell = sp.gcd(d1k[0], d1k[1])
assert simp0(g_ell - DELTA) or simp0(g_ell + DELTA)
print(f"  gcd of the d1 entries = {sp.expand(g_ell)}   (up to a unit).")


def rank1_dims(tv):
    return coh_dims(sp.Matrix([[tv]]), sp.Matrix([[tv]]))


# exact dimension table; Delta vanishes exactly at the primitive 6th roots
assert simp0(DELTA.subs(t, Z6)) and simp0(DELTA.subs(t, Z6b))
assert sp.factor_list(DELTA)[1][0][0] == sp.cyclotomic_poly(6, t)
CASES1 = [
    (sp.Integer(1), "t = 1            (s = 0, untwisted)", (1, 1, 0)),
    (Z6, "t = zeta_6       (s = 1/6)", (0, 1, 1)),
    (Z6b, "t = zeta_6^5     (s = 5/6)", (0, 1, 1)),
    (sp.Integer(-1), "t = -1           (s = 1/2)", (0, 0, 0)),
    (Z3, "t = zeta_3       (s = 1/3)", (0, 0, 0)),
    (Z3b, "t = zeta_3^2     (s = 2/3)", (0, 0, 0)),
    (sp.Integer(2), "t = 2            (generic)", (0, 0, 0)),
    (sp.Rational(7, 5), "t = 7/5          (generic)", (0, 0, 0)),
    (sp.I, "t = i            (s = 1/4)", (0, 0, 0)),
]
print("\n  dimension table (h0, h1, h2), all exact:")
for tv, label, expected in CASES1:
    dims = rank1_dims(tv)
    print(f"    {label:<38} {dims}")
    assert dims == expected, (label, dims)
print("  => JUMP LOCUS in s (t = e^{2 pi i s}, s in [0,1)):")
print("       s = 0:          (1, 1, 0)   [t = 1: untwisted H^*(M)]")
print("       s = 1/6, 5/6:   (0, 1, 1)   [t = primitive 6th roots =")
print("                                    roots of Delta(t)]")
print("       all other s:    (0, 0, 0).")
print("  Generic Kummer twist KILLS all cohomology: 0 twisted 'master")
print("  integrals' at generic s (= |chi(M)| = 0).")

# ---------------------------------------------------------------------------
print("\n=== 3. Independent cross-check: Milnor fibration + Wang"
      " sequence ===")
# premise: the wall is affinely the quasi-homogeneous curve 27(W^2 - U^3)
P2 = 27 * u**2 + 16 * u - 18 * u * w + w**3 - w**2
Uex = (4 - 3 * w) / 9
Wex = (27 * u - 9 * w + 8) / 27
assert sp.expand(P2 - 27 * (Wex**2 - Uex**3)) == 0
print("  P2 = 27 (W^2 - U^3) in the affine A2 coordinates (I4) of")
print("  docs/WALL_COMPLEMENT.md, so P2: M -> C* is the global Milnor")
print("  fibration of the cusp x^2 + y^3 [Milnor 1968, section 9 for")
print("  weighted-homogeneous polynomials]; M ~ mapping torus of the")
print("  geometric monodromy h of the Milnor fiber F, and the Kummer")
print("  system C_t is pulled back from C* (a wall meridian maps to a")
print("  degree-1 loop).")

# monodromy eigenvalues on H^1(F) via the Brieskorn formula for x^a + y^b:
# { zeta_a^i zeta_b^j : 1 <= i <= a-1, 1 <= j <= b-1 } = {-zeta_3^j}
eigs = [-Z3, -Z3b]
assert all(simp0(DELTA.subs(t, e)) for e in eigs)
assert not simp0(eigs[0] - eigs[1])
print("  Brieskorn eigenvalues of h* on H^1(F) for x^2 + y^3:")
print("    {-zeta_3, -zeta_3^2} = the two primitive 6th roots = the roots")
print("    of Delta(t).  mu = (2-1)(3-1) = 2 = b1(F);  b0(F) = 1.")

# realize h* exactly as the companion matrix of Delta
H1 = sp.Matrix([[0, -1], [1, 1]])
assert H1.charpoly(t).as_expr() == DELTA
assert H1**6 == sp.eye(2) and H1**3 == -sp.eye(2)
assert sp.expand((t * H1 - I2).det() - DELTA) == 0
print("  h* on H^1(F) = companion matrix of Delta: order 6, h*^3 = -1;")
print("  det(t h* - 1 | H^1 F) = Delta(t)  and  det(t - 1 | H^0 F) = t-1:")
print("  IDENTICAL to the Fox jump polynomials d1 ~ Delta, d0 ~ (t-1).")


def wang_dims(tv):
    """H^*(M; C_t) from the Wang sequence of F -> M -> C*:
    0 -> coker(t h* - 1 | H^{k-1} F) -> H^k(M) -> ker(t h* - 1 | H^k F)
    -> 0."""
    a2 = (tv * H1 - I2)
    k1 = 2 - exact_rank(a2)              # dim ker = dim coker on H^1(F)
    k0 = 1 if simp0(tv - 1) else 0       # ker/coker of (t - 1) on H^0(F)
    return (k0, k0 + k1, k1)


print("\n  Wang vs Fox, point by point (must agree -- and do):")
for tv, label, _ in CASES1:
    dw, df_ = wang_dims(tv), rank1_dims(tv)
    print(f"    {label:<38} Wang {dw}  Fox {df_}")
    assert dw == df_, (label, dw, df_)
print("  => the two computations AGREE at every tested t, and their jump")
print("     loci agree identically as polynomials.")

# ---------------------------------------------------------------------------
print("\n=== 4. The S3 local systems: standard, sign, permutation,"
      " regular ===")
r1 = sp.Matrix([[-1, 1], [0, 1]])
r2 = sp.Matrix([[1, 0], [1, -1]])
assert r1**2 == sp.eye(2) and r2**2 == sp.eye(2)
assert r1 * r2 * r1 == r2 * r1 * r2
cox = r1 * r2
assert cox**3 == sp.eye(2) and cox != sp.eye(2)
print("  reflection matrices verified: r1^2 = r2^2 = 1, braid relation,")
print("  (r1 r2)^3 = 1 (Coxeter element of order 3) => rho factors")
print("  through B3 ->> S3.")

dims_std = coh_dims(r1, r2)
print(f"  STANDARD (reflection) local system S:   H^* = {dims_std}")
assert dims_std == (0, 1, 1)
print("    => dim H^1(M; S) = dim H^2(M; S) = 1: ONE twisted class in")
print("       each degree -- the 'one master integral' of the reflection")
print("       channel proposed in WALL_COMPLEMENT.md section 5.")

dims_sgn = coh_dims(sp.Matrix([[-1]]), sp.Matrix([[-1]]))
print(f"  sign local system:                      H^* = {dims_sgn}")
assert dims_sgn == (0, 0, 0) == rank1_dims(sp.Integer(-1))
print("    (= the rank-1 case at t = -1; Delta(-1) = 3 != 0.)")

P1 = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])       # (1 2)
P2m = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])      # (2 3)
dims_perm = coh_dims(P1, P2m)
print(f"  permutation (sheet) local system L:     H^* = {dims_perm}")
assert dims_perm == (1, 2, 1)
add = tuple(x + y for x, y in zip(dims_triv, dims_std))
assert dims_perm == add
print(f"    = triv (+) std componentwise: {dims_triv} + {dims_std} ="
      f" {add}: additivity of twisted cohomology under the decomposition")
print("    L = triv (+) S of WALL_COMPLEMENT.md section 4, VERIFIED.")

# regular representation C[S3]: Shapiro cross-check against the pure
# braid group / A2 arrangement complement
elems = [(0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 1, 0), (1, 2, 0), (2, 0, 1)]


def reg_matrix(g):
    Mr = sp.zeros(6, 6)
    for j, h in enumerate(elems):
        gh = tuple(g[h[i]] for i in range(3))
        Mr[elems.index(gh), j] = 1
    return Mr


dims_reg = coh_dims(reg_matrix((1, 0, 2)), reg_matrix((0, 2, 1)))
print(f"  regular local system C[S3]:             H^* = {dims_reg}")
assert dims_reg == (1, 3, 2)
# Shapiro: H^*(B3; C[S3]) = H^*(P3; C), P3 = pure braid group = pi_1 of
# the A2 arrangement complement Y.  Independent computation of H^*(Y) by
# Orlik-Solomon/Whitney from the intersection lattice (3 concurrent
# transverse lines in the sum-zero plane):
lines = [(1, -1), (2, 1), (1, 2)]        # x1-x2, x1-x3, x2-x3 at x3=-x1-x2
for (a1, b1), (a2, b2) in combinations(lines, 2):
    assert a1 * b2 - a2 * b1 != 0        # pairwise transverse
mu_origin = -(1 + 3 * (-1))              # Moebius fn of the rank-2 flat
betti_Y = (1, 3, abs(mu_origin))
assert betti_Y == (1, 3, 2)
qq = sp.symbols("qq")
assert sp.Poly(sp.expand((1 + qq) * (1 + 2 * qq)), qq).all_coeffs()[::-1] \
    == [1, 3, 2]                          # Y ~ C* x (C minus 2 points)
assert dims_reg == betti_Y
add6 = tuple(a_ + b_ + 2 * c_ for a_, b_, c_ in
             zip(dims_triv, dims_sgn, dims_std))
assert dims_reg == add6
print("    SHAPIRO CROSS-CHECK: H^*(B3; C[S3]) = H^*(P3; C) = H^*(A2")
print("    arrangement complement) = (1, 3, 2) [Moebius/Orlik-Solomon:")
print("    mu(origin flat) = 2; also = coefficients of (1+q)(1+2q), the")
print("    model Y ~ C* x (C \\ {0,1})].  And (1,3,2) = triv + sign +")
print("    2*std componentwise: the full character theory is consistent.")

# ---------------------------------------------------------------------------
print("\n=== 5. The reduced Burau family rho_t ===")
B1 = sp.Matrix([[-t, 1], [0, 1]])
B2 = sp.Matrix([[1, 0], [t, -t]])
d0b, d1b = fox_maps(B1, B2)
d0b = d0b.applyfunc(sp.expand)
d1b = d1b.applyfunc(sp.expand)
print("  braid relation holds for ALL t (asserted symbolically);"
      " t = 0 excluded (rho_t must be invertible: det rho_t(s_i) = -t).")
assert sp.expand(B1.det() + t) == 0 and sp.expand(B2.det() + t) == 0
print(f"  d1 = {d1b.tolist()}")
assert d1b == sp.Matrix([[0, -t, -1, 0], [0, 1, t**2, 0]])

# specializations FIRST (what recovers cases 1 and 2):
assert B1.subs(t, 1) == r1 and B2.subs(t, 1) == r2
print("  SPECIALIZATION t = +1: Burau(1) IS the reflection rep, on the")
print("  nose (matrix equality) -- case 3 sits inside the family at t=1.")
print("  SPECIALIZATION det: det o rho_t = lambda_{-t} (sigma_i -> -t):")
print("  the rank-1 Kummer family of section 2 is recovered up to the")
print("  classical Burau sign convention t_Burau = -t_Kummer.")
C1m = (B1 * B2).subs(t, -1)
assert C1m**3 == -sp.eye(2) and C1m**6 == sp.eye(2)
assert (B1.subs(t, -1) - I2) != sp.zeros(2, 2)
assert (B1.subs(t, -1) - I2)**2 == sp.zeros(2, 2)
print("  SPECIALIZATION t = -1 (honest correction of the guess in the")
print("  task): Burau(-1) is NOT the reflection rep -- it is the integral")
print("  SL(2,Z) representation: sigma_i -> unipotent (infinite order),")
print("  full twist (s1 s2)^3 -> -1 != 1 (the center acts nontrivially,")
print("  so rho_{-1} does NOT factor through S3).  It is the homological")
print("  monodromy of the torus-knot picture, not the sheet monodromy.")


def minors_gcd(M, k):
    g = sp.Integer(0)
    for rows in combinations(range(M.rows), k):
        for cols in combinations(range(M.cols), k):
            g = sp.gcd(g, sp.expand(M.extract(rows, cols).det()))
    return sp.expand(g)


g0 = minors_gcd(d0b, 2)
g1 = minors_gcd(d1b, 2)
print(f"  gcd of 2x2 minors:  d0: {g0},   d1: {g1}")
assert simp0(g0 - (t**2 + t + 1)) or simp0(g0 + (t**2 + t + 1))
assert simp0(g1 - (t**3 - 1)) or simp0(g1 + (t**3 - 1))
assert sp.expand(t**3 - 1 - (t - 1) * (t**2 + t + 1)) == 0
assert sp.expand(DELTA.subs(t, -t) - (t**2 + t + 1)) == 0
# d0 and d1 never vanish identically at any t (constant entries present):
assert 1 in set(d0b) and -1 in set(d1b)
print("  => rank d0 drops (2 -> 1) exactly on {t^2 + t + 1 = 0} and")
print("     rank d1 drops (2 -> 1) exactly on {t^3 = 1}; neither map is")
print("     ever zero.  JUMP LOCUS of H^*(B3; rho_t) = {t : t^3 = 1},")
print("     and  t^2 + t + 1 = Delta(-t):  the trefoil Alexander")
print("     polynomial again, through the Burau sign convention.")


def burau_dims(tv):
    return coh_dims(B1.subs(t, tv), B2.subs(t, tv))


CASES5 = [
    (sp.Integer(1), "t = 1    (the reflection rep)", (0, 1, 1)),
    (Z3, "t = zeta_3", (1, 2, 1)),
    (Z3b, "t = zeta_3^2", (1, 2, 1)),
    (sp.Integer(-1), "t = -1   (SL(2,Z) rep)", (0, 0, 0)),
    (Z6, "t = zeta_6", (0, 0, 0)),
    (sp.Integer(2), "t = 2    (generic)", (0, 0, 0)),
    (sp.Rational(1, 2), "t = 1/2  (generic)", (0, 0, 0)),
    (sp.Rational(7, 5), "t = 7/5  (generic)", (0, 0, 0)),
]
print("\n  dimension table (h0, h1, h2), all exact:")
for tv, label, expected in CASES5:
    dims = burau_dims(tv)
    print(f"    {label:<38} {dims}")
    assert dims == expected, (label, dims)
assert burau_dims(sp.Integer(1)) == dims_std

# structure at the cube roots: Burau becomes reducible
vfix = sp.Matrix([1, 1 + Z3])
assert is_zero_mat(B1.subs(t, Z3) * vfix - vfix)
assert is_zero_mat(B2.subs(t, Z3) * vfix - vfix)
assert simp0(DELTA.subs(t, -Z3))          # -zeta_3 is a primitive 6th root
print("  STRUCTURE at t = zeta_3: rho_t fixes v = (1, 1+zeta_3) --")
print("  reducible: 0 -> C_triv -> rho_t -> C_{-zeta_3} -> 0, with -zeta_3")
print("  a PRIMITIVE 6TH ROOT (a root of Delta).  Long exact sequence:")
print("  (1,1,0) and (0,1,1) stack to (1,2,1) -- exactly the measured")
print("  dimensions (the connecting maps vanish).")

# duality check t <-> 1/t (dual local systems): dimensions must match
for tv in [Z6, Z3, sp.Integer(2), sp.Rational(7, 5)]:
    inv = sp.simplify(1 / tv)
    assert rank1_dims(tv) == rank1_dims(inv)
    assert burau_dims(tv) == burau_dims(inv)
print("  duality check: dims at t and 1/t agree for both families (jump")
print("  loci closed under t -> 1/t), as Poincare-Lefschetz duality")
print("  H^k(M; L) ~ H^{4-k}_c(M; L^dual)* requires.")

# ---------------------------------------------------------------------------
print("\n=== 6. De Rham side: explicit Euler primitives and the two")
print("        spectral classes ===")
f = W**2 - U**3
fU, fW = sp.diff(f, U), sp.diff(f, W)
assert sp.expand(2 * U * fU + 3 * W * fW - 6 * f) == 0     # E(f) = 6f
print("  f = W^2 - U^3 (= P2/27 up to the affine iso, section 3);")
print("  weights (U, W) = (2, 3), Euler field E = 2U d_U + 3W d_W,")
print("  E(f) = 6 f.  Twisted differential nabla_s = d + s dlog f.")


def nabla1(A, B):
    """nabla_s(A dU + B dW) = [...] dU^dW; returns the coefficient."""
    return sp.cancel(sp.diff(B, U) - sp.diff(A, W)
                     + s * (fU * B - fW * A) / f)


def nabla0(g):
    """nabla_s(g) = (.) dU + (.) dW on 0-forms."""
    return (sp.cancel(sp.diff(g, U) + s * g * fU / f),
            sp.cancel(sp.diff(g, W) + s * g * fW / f))


# (a) the H^1 story: dlog f is the untwisted generator, killed by any
#     honest twist
gA, gB = nabla0(sp.Integer(1))
assert sp.cancel(gA - s * fU / f) == 0 and sp.cancel(gB - s * fW / f) == 0
print("  (a) nabla_s(1) = s dlog f: for s not in Z the class [dlog f]")
print("      is nabla_s-exact (df/f = nabla_s(1/s)) -- matching the Fox")
print("      result h1 = 1 at t = 1 and (generically) 0 otherwise.")

# (b) the H^2 story: Euler primitives, exactly
etaA, etaB = -3 * W / f, 2 * U / f          # eta = iota_E(dU^dW)/f
assert sp.cancel(nabla1(etaA, etaB) - (6 * s - 1) / f) == 0
assert sp.cancel(nabla1(U * etaA, U * etaB) - (6 * s + 1) * U / f) == 0
assert sp.cancel(nabla1(etaA, etaB).subs(s, 0) + 1 / f) == 0
print("  (b) eta = iota_E(dU^dW)/f = (2U dW - 3W dU)/f:")
print("        nabla_s(eta)   = (6s - 1) dU^dW / f")
print("        nabla_s(U eta) = (6s + 1) U dU^dW / f")
print("      => [dU^dW/f] is nabla_s-exact for s != 1/6, and [U dU^dW/f]")
print("         for s != -1/6 ~ 5/6; at s = 0 the would-be CANONICAL FORM")
print("         dU^dW/f = d(-eta) is EXACT, with the explicit rational")
print("         log-pole primitive -eta.  Zero periods over every closed")
print("         2-cycle: the global face of the residueless double pole")
print("         of docs/POSITIVE_GEOMETRY.md section 3.")

# general quasi-homogeneous bookkeeping: the full (k, m) grid
for a_ in range(3):
    for b_ in range(3):
        for m_ in range(1, 4):
            g_ = U**a_ * W**b_
            k_ = 2 * a_ + 3 * b_
            lhs = nabla1(g_ * etaA * f / f**m_, g_ * etaB * f / f**m_)
            rhs = (k_ + 5 - 6 * m_ + 6 * s) * g_ / f**m_
            assert sp.cancel(lhs - rhs) == 0
print("      grid identity VERIFIED (a, b <= 2, m <= 3):")
print("        nabla_s(U^a W^b iota_E(dU^dW)/f^m)")
print("          = (2a + 3b + 5 - 6m + 6s) U^a W^b dU^dW / f^m,")
print("      so pole-order reduction fails only at s = (6m - 5 - k)/6,")
print("      i.e. s = (1 - k)/6 mod Z with k = 2a + 3b.")

# the Milnor ring of the cusp selects which numerators are genuine
# (f_U, f_W) = (-3U^2, 2W) is already a Groebner basis: leading monomials
# U^2 and W are coprime
assert sp.reduced(U**2, [fU, fW], U, W)[1] == 0
assert sp.reduced(W, [fU, fW], U, W)[1] == 0
assert sp.reduced(U, [fU, fW], U, W)[1] == U
print("      Milnor ring C[U,W]/(f_U, f_W) has basis {1, U} (mu = 2);")
print("      the basis numerators give the resonant exponents")
print("        g = 1 (k=0): s in 1/6 + Z;   g = U (k=2): s in 5/6 + Z,")
print("      matching e^{2 pi i s} = the primitive 6th roots (the Milnor")
print("      monodromy eigenvalues = the cusp spectrum mod Z) -- the jump")
print("      locus of sections 2-3 exactly.")

# non-basis numerators do NOT create jumps: Jacobian-ideal reduction
assert sp.cancel(nabla1(1 / f, sp.Integer(0))
                 - (1 - s) * 2 * W / f**2) == 0
print("      Non-basis check: nabla_s(dU/f) = 2(1-s) W dU^dW/f^2, so the")
print("      class [W dU^dW/f^2] is exact for ALL s != 1 -- in particular")
print("      at s = 2/3 where its Euler chain fails (coefficient 6s - 4):")
print("      Euler-chain failure at non-basis numerators is NOT a")
print("      cohomology jump.  (Caveat, stated honestly: these are")
print("      explicit consistency computations on the de Rham side, not a")
print("      full independent computation of algebraic de Rham H^*; the")
print("      dictionary twisted-de-Rham = local-system cohomology is the")
print("      Deligne/ESV comparison, cited not re-proved.)")

# ---------------------------------------------------------------------------
print("\n=== 7. Summary ===")
print("  local system                     (h0, h1, h2)      jump locus")
print("  " + "-" * 68)
print("  C_t, generic t                   (0, 0, 0)")
print("  C_t, t = 1        [s in Z]       (1, 1, 0)     (t - 1)")
print("  C_t, Delta(t)=0   [s in ±1/6+Z]  (0, 1, 1)     t^2 - t + 1")
print("  sign = C_{-1}                    (0, 0, 0)")
print("  standard S                       (0, 1, 1)")
print("  permutation = triv + std         (1, 2, 1)")
print("  regular C[S3] (= H^* of P3)      (1, 3, 2)")
print("  Burau_t, generic t               (0, 0, 0)")
print("  Burau_1 = standard               (0, 1, 1)     t^3 - 1")
print("  Burau_{zeta_3^(+-1)}             (1, 2, 1)     t^3 - 1")
print("  chi = 0 in every single case (as it must: chi(M) = 0).")
print("  Fox calculus and the Wang sequence AGREE everywhere tested;")
print("  their jump polynomials agree identically.")

print(f"\nAll assertions passed.  Total wall time: {time.time() - T0:.1f} s")
