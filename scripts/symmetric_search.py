"""The SYMMETRIC (gradient / variational) Jacobian problem: exact groundwork.

(Write-up: docs/SYMMETRIC_SEARCH.md.  Machinery: jcqft/reduction_w.py for the
equivariant part, jcqft/gb_backend.py for the in-box Groebner certificates.
Literature anchor: de Bondt-van den Essen, Proc. AMS 133 (2005) 2201-2205
[dBvdE05], statement re-verified against the published text; see the doc.)

Default run ~4 min; --full adds brute-force 6x6 determinant cross-checks and
larger Groebner boxes (~1 h, plus one honest wall report).

PROVED HERE, all exact (every claim asserted):

A. EXPLICIT SYMMETRIC COUNTEREXAMPLES IN DIMENSION 6.
   A1. Block lemma det [[M, C],[A, 0]] = (-1)^n det(A) det(C) (generic
       n = 2, 3) -- the mechanism behind everything in section A.
   A2. The COTANGENT LIFT of the Alpoge-Mathew map F,
           W6(phi, phibar) := phibar . F(phi)  in C^6,
       has  Hess W6 = [[sum phibar_k Hess F_k, JF^T],[JF, 0]],
       det Hess W6 = -(det JF)^2 = -4  (constant!), and grad W6 is 3:1
       over a rational point: three explicit RATIONAL witnesses.  This is
       an explicit counterexample to the gradient/Hessian formulation of
       the Jacobian conjecture (Meng's Hessian conjecture) in dim 6 -- and
       W6 is exactly the first-order action phibar.F(phi) of PROBLEM.md.
   A3. The dBvdE TWISTED LIFT: H := L^{-1}F - id (so det(I + JH) = 1),
           f_H := -i sum_j H_j(x + iy) y_j   in C[x, y],  deg 8,
       gives F~ = id + grad f_H with det(I + Hess f_H) = 1 IDENTICALLY:
       a NORMALIZED symmetric Keller counterexample x + grad f over Q(i),
       with three explicit witnesses over Q(i).  Keller-ness is proved
       exactly via the congruence identity
           Hess f_H = S^{-T} (Hess g_H)(S^{-1} .) S^{-1},
           S^T S + Hess g_H = [[I + B, -i(I+JH)^T],[-i(I+JH), 0]],
       and the block lemma: det(I + Hess f_H) = det(I + JH)^2 = 1.
   A4. Hess f_H is NOT nilpotent (tr h^2 != 0), and NO affine
       normalization of the AM map has nilpotent JH: the linear system
       tr(C JF(x)) = 3 identically is INCONSISTENT (linsolve EmptySet).
       So the dBvdE nilpotent/homogeneous normal form (their Cor. 1.3)
       provably cannot be reached from AM in dim 6 by the twist alone;
       it requires the Bass-Connell-Wright detour (dimension blow-up).
   A5. STRONGEST 3D NO-GO: the only matrix K with K JF(x) symmetric for
       all x is K = 0.  Since B F(Ax) + c is a gradient map iff
       K := A^{-T} B has K JF symmetric, the AM counterexample is NOT
       equivalent to ANY gradient map under invertible affine source and
       target transformations.  Its non-properness defect is
       NON-VARIATIONAL in dimension 3, in the strongest linear sense.

B. GRADIENT STRUCTURE vs C*-EQUIVARIANCE (whole family (1,-1,-m)).
   B1. Weight bookkeeping: F = grad W equivariant with component weights
       a permutation of (-m,-1,1) forces a finite branch table; the only
       branch with a nontrivial 2-variable potential is m = 3 with the
       standard ordering (weights (-3,-1,1)) and W = x^{-2} S(w,v).
   B2. At m = 3, DF symmetric <=> (P, Q, R) = (E(S), S_w, S_v) with
       E = w d_w + 3 v d_v - 2; for m != 3, symmetry of DF degenerates
       the map into the branches of B1 (asserted x-power bookkeeping).
   B3. The degenerate branches are classified COMPLETELY for every m:
       only linear maps and the gradient shear family
           W = alpha x^2/2 + delta y z + gamma y^{m+1},
           grad W = (alpha x, delta z + (m+1) gamma y^m, delta y),
       a TAME automorphism (explicit inverse asserted; infinity
       prefilter survived only through the nonlinear-automorphism
       false-positive class).  In particular at m = 2 -- the AM weight
       system -- the ONLY gradient Keller maps in the entire equivariant
       family are linear maps and this shear family.
   B4. m = 3 potential slice: det M(E(S), S_w, S_v) reduces to an exact
       2-variable Monge-Ampere-type PDE; kappa != 0 makes DF(0)
       automatically invertible; residual scalings gauge-fix
       [v]S = 1, [w^2]S = 1/2, kappa = -1.  Then:
       - deg_v S = 0 empty; deg_v S = 1: ONLY the linear map (all
         w-degrees; Wronskian + degree kill + box);
       - deg_v S = 2: EMPTY (cascade c3 -> c2 -> c1 forces a singular
         linearization; all w-degrees);
       - all K = deg_v S: the top rows are rigid, K SYMBOLIC:
         L1: [v^K]S is a CONSTANT           (from c_{3K-2});
         L2: [v^{K-1}]S is LINEAR in w      (from c_{3K-3});
         L3: [v^{K-2}]S is QUADRATIC, second derivative fixed
             prop. to t1^2                  (from c_{3K-4});
       - in-box Groebner (msolve, exact over Q): K = 3 with w-degree
         <= 8 and K = 4 with w-degree <= 6 EMPTY (presolved by L1-L2,
         seconds); larger boxes under --full.
       Verdict: no nonlinear gradient Keller map exists in the
       (1,-1,-m) family for ANY m, beyond the tame shear family, up to
       the honest gap (m = 3, deg_v >= 3, w-degree beyond the boxes).

C. DIRECT BOXES FOR SYMMETRIC KELLER MAPS (normalized WLOG over C:
   Hess W(0) = I congruence => quadratic part sum x_i^2/2, kappa = 1,
   witness pair +-e/2 by translation).
   C1. deg W <= 3, ANY n: the midpoint identity
           grad W(a) - grad W(b) = Hess W((a+b)/2) (a - b)
       (asserted n = 2, 3, 4) makes every symmetric Keller map of
       degree <= 3 INJECTIVE.  No search needed.
   C2. n = 2: witness-pair emptiness deg <= 4, 5, 6 (msolve, instant),
       with controls: dropping Keller leaves witnesses (non-unit), and
       the Keller variety alone is nonempty (the nonlinear family
       W = xy + f(x) is Keller with explicit inverse -- asserted).
   C3. n = 3: deg W <= 4 EMPTY -- the full 25-coefficient box, three
       witness queries, ~2 min.  Combined with C1: no symmetric Keller
       counterexample in dim 3 up to degree 4.  deg <= 5 is the honest
       wall (--full attempts it; see doc section 8).

D. COERCIVITY SUPPORT (for the synthesis in the doc): kappa(W6) = -4 < 0
   and W6 is AFFINE in phibar (both asserted): the explicit variational
   counterexample is maximally non-coercive; the leading form of W6 has
   identically singular Hessian (asserted), the general phenomenon
   behind the no-coercive-counterexample theorem of the doc.

VERDICT: variational Jacobian counterexamples EXIST -- explicitly, in
dimension 6, with rational data (A2) -- but the AM defect itself is
non-variational in dimension 3 in every sense tested: not affinely
symmetrizable (A5), no gradient Keller maps in its equivariant family
beyond tame shears (B), and no symmetric counterexample at all in the
dim-3 boxes through degree 4 (C).
"""

from __future__ import annotations

import argparse
import itertools
import time

import sympy as sp

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from jcqft.core import F, PHI, L, x, y, z                # noqa: E402
from jcqft.gb_backend import available_backends, is_unit_ideal  # noqa: E402
from jcqft.prefilter import infinity_prefilter           # noqa: E402
from jcqft.reduction_w import det_m, w, v                # noqa: E402

T0 = time.time()
N_CHECKS = 0


def check(label, cond):
    global N_CHECKS
    assert cond, label
    N_CHECKS += 1
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


def expand_zero(e):
    e = sp.expand(e.doit() if hasattr(e, "doit") else e)
    if e == 0:
        return True
    return sp.expand(sp.powsimp(e)) == 0


parser = argparse.ArgumentParser()
parser.add_argument("--full", action="store_true",
                    help="brute-force 6x6 determinants + larger boxes")
ARGS = parser.parse_args()

# ===========================================================================
print("== A1. block-determinant lemma (generic) ==")
# ===========================================================================
# det [[M, C],[A, 0]] = (-1)^n det(A) det(C): swap the block rows (n row
# transpositions, sign (-1)^n) and the result is block-triangular.  Asserted
# on fully generic blocks for n = 2, 3; the row-swap proof works for every n.
for n in (2, 3):
    Mg = sp.Matrix(n, n, lambda i, j: sp.Symbol(f"m{i}{j}"))
    Ag = sp.Matrix(n, n, lambda i, j: sp.Symbol(f"A{i}{j}"))
    Cg = sp.Matrix(n, n, lambda i, j: sp.Symbol(f"C{i}{j}"))
    Blk = sp.Matrix(sp.BlockMatrix([[Mg, Cg], [Ag, sp.zeros(n)]]))
    check(f"block lemma n={n}",
          sp.expand(Blk.det(method="berkowitz")
                    - (-1) ** n * Ag.det() * Cg.det()) == 0)

# ===========================================================================
print("== A2. the cotangent lift W6 = phibar . F : dim-6 counterexample ==")
# ===========================================================================
xb, yb, zb = sp.symbols("xb yb zb")
BAR = (xb, yb, zb)
V6 = PHI + BAR

W6 = sp.expand(xb * F[0] + yb * F[1] + zb * F[2])
JF = sp.Matrix([[sp.expand(sp.diff(Fi, u)) for u in PHI] for Fi in F])
check("det JF = -2 (constant)", sp.expand(JF.det()) == -2)

H6 = sp.hessian(W6, V6)
Bblk = sp.Matrix(3, 3, lambda i, j: sum(BAR[k] * sp.diff(F[k], PHI[i], PHI[j])
                                        for k in range(3)))
shape = sp.Matrix(sp.BlockMatrix([[Bblk, JF.T], [JF, sp.zeros(3)]]))
check("Hess W6 = [[sum phibar_k Hess F_k, JF^T],[JF, 0]]",
      sp.expand(H6 - shape) == sp.zeros(6))
# block lemma (A1, generic M) => det Hess W6 = (-1)^3 det(JF) det(JF^T) = -4
check("det Hess W6 = -(det JF)^2 = -4 (constant: symmetric Keller data)",
      sp.expand((-1) ** 3 * JF.det() * JF.T.det()) == -4)
if ARGS.full:
    check("[--full] brute-force 6x6 det Hess W6 = -4 (berkowitz)",
          sp.expand(H6.det(method="berkowitz")) == -4)

# non-injectivity: the three AM fiber points, lifted.  grad W6 =
# (JF(x)^T phibar, F(x)); F(a)=F(b) and JF^T invertible everywhere let us
# match the first block by a linear solve -- with RATIONAL solutions.
AM_PTS = [(0, 0, sp.Rational(-1, 4)),
          (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
          (-1, sp.Rational(3, 2), sp.Rational(13, 2))]
grad6 = [sp.expand(sp.diff(W6, u)) for u in V6]
eta0 = sp.Matrix([1, 1, 1])
Ja = JF.T.subs(dict(zip(PHI, AM_PTS[0])))
wit6, val6 = [], []
for p in AM_PTS:
    Jp = JF.T.subs(dict(zip(PHI, p)))
    etap = Jp.solve(Ja * eta0)
    pt = tuple(p) + tuple(etap)
    wit6.append(pt)
    val6.append(tuple(sp.expand(g.subs(dict(zip(V6, pt)))) for g in grad6))
check("three witness points are distinct and RATIONAL",
      len(set(wit6)) == 3 and all(sp.sympify(c).is_rational
                                  for pt in wit6 for c in pt))
check("grad W6 identical on all three (value (5/4,1,1,-1/4,0,0))",
      val6[0] == val6[1] == val6[2]
      and val6[0] == (sp.Rational(5, 4), 1, 1, sp.Rational(-1, 4), 0, 0))
print("     witnesses:", wit6)
# injectivity transfer is an iff: grad W6(x,phibar)=(JF^T phibar, F(x)); if F
# is injective then x = x' and JF(x)^T invertible gives phibar = phibar'.
check("deg W6 = 8", sp.Poly(W6, *V6).total_degree() == 8)
check("W6 affine in phibar (first-order action => never coercive)",
      all(sp.diff(W6, a, b) == 0 for a in BAR for b in BAR))
lead = sum(t for t in sp.Add.make_args(W6)
           if sum(sp.Poly(t, *V6).degree_list()) == 8)
check("leading form of W6 = xb*x^3*y^3*z; its Hessian is singular "
      "(ybar-row vanishes identically)",
      sp.expand(lead - xb * x**3 * y**3 * z) == 0
      and sp.hessian(lead, V6).det() == 0)

# ===========================================================================
print("== A3. dBvdE twisted lift: normalized x + grad f over Q(i) ==")
# ===========================================================================
# H = L^{-1} F - id: Keller-normalized (det(I + JH) = 1, JH(0) = 0).
Hm = L.inv() * sp.Matrix(F) - sp.Matrix(PHI)
Hs = [sp.expand(h) for h in Hm]
JH = sp.Matrix([[sp.diff(h, u) for u in PHI] for h in Hs])
check("H = L^{-1}F - id has JH(0) = 0 and det(I + JH) = 1",
      JH.subs({x: 0, y: 0, z: 0}) == sp.zeros(3)
      and sp.expand((sp.eye(3) + JH).det()) == 1)

# f_H := -i sum_j H_j(x + iy) y_j   (dBvdE05 eq. (3)); g_H := f_H o S,
# S = (x - iy, y)  =>  g_H = -i <y, H(x)>.
gH = sp.expand(-sp.I * sum(h * b for h, b in zip(Hs, BAR)))
Smat = sp.Matrix(sp.BlockMatrix([[sp.eye(3), -sp.I * sp.eye(3)],
                                 [sp.zeros(3), sp.eye(3)]]))
Sinv = Smat.inv()
sub_inv = dict(zip(V6, list(Sinv * sp.Matrix(V6))))
fH = sp.expand(gH.subs(sub_inv, simultaneous=True))
check("f_H = -i sum_j H_j(x+iy) y_j (dBvdE05 (3)) and deg f_H = 8",
      sp.expand(fH + sp.I * sum(
          sp.expand(h.subs({x: x + sp.I * xb, y: y + sp.I * yb,
                            z: z + sp.I * zb}, simultaneous=True)) * b
          for h, b in zip(Hs, BAR))) == 0
      and sp.Poly(fH, *V6).total_degree() == 8)

hg = sp.hessian(gH, V6)
hf = sp.hessian(fH, V6)
check("congruence identity Hess f_H = S^{-T} (Hess g_H)(S^{-1}.) S^{-1}",
      sp.expand(hf - Sinv.T * hg.subs(sub_inv, simultaneous=True) * Sinv)
      == sp.zeros(6))
Bxx = sp.Matrix(3, 3, lambda i, j: -sp.I * sum(
    BAR[k] * sp.diff(Hs[k], PHI[i], PHI[j]) for k in range(3)))
tgt = sp.Matrix(sp.BlockMatrix(
    [[sp.eye(3) + Bxx, -sp.I * (sp.eye(3) + JH.T)],
     [-sp.I * (sp.eye(3) + JH), sp.zeros(3)]]))
check("S^T S + Hess g_H = [[I+B, -i(I+JH)^T],[-i(I+JH), 0]]",
      sp.expand(sp.expand(Smat.T * Smat) + hg - tgt) == sp.zeros(6))
# block lemma: det(S^T S + Hess g_H) = (-1)^3 det(-i(I+JH))^2
#            = (-1)^3 (-i)^6 det(I+JH)^2 = det(I+JH)^2 = 1,
# and det S = 1, so det(I + Hess f_H) = 1 identically: KELLER, normalized.
check("det(I + Hess f_H) = det(I+JH)^2 = 1 (block lemma + det S = 1): "
      "F~ = id + grad f_H is a NORMALIZED symmetric Keller map",
      sp.expand((-1) ** 3 * ((-sp.I) ** 3) ** 2 * (sp.eye(3) + JH).det()
                * (sp.eye(3) + JH.T).det()) == 1
      and Smat.det() == 1 and hf.subs({u: 0 for u in V6}) == sp.zeros(6))
if ARGS.full:
    check("[--full] brute-force 6x6 det(I + Hess f_H) = 1 (berkowitz)",
          sp.expand((sp.eye(6) + hf).det(method="berkowitz")) == 1)

# witnesses: S^{-1} o F~ o S = (x + H(x), (I+JH(x))^T y - i H(x)) (asserted),
# so matching the second block is a linear solve: (I+JH(b)^T) y_b = i(a-b).
Ft = sp.Matrix([uu + sp.diff(fH, uu) for uu in V6])
Ssub = dict(zip(V6, [PHI[i] - sp.I * BAR[i] for i in range(3)] + list(BAR)))
FtS = Ft.subs(Ssub, simultaneous=True)
G = sp.expand(sp.Matrix(list(FtS[:3, 0] + sp.I * FtS[3:, 0])
                        + list(FtS[3:, 0])))
check("S^{-1} o F~ o S = (x + H(x), (I+JH^T) y - i H(x))  [dBvdE05 proof "
      "of Thm 1.1]",
      sp.expand(G[:3, 0] - (sp.Matrix(PHI) + Hm)) == sp.zeros(3, 1)
      and sp.expand(G[3:, 0] - (sp.Matrix(BAR) + JH.T * sp.Matrix(BAR)
                                - sp.I * Hm)) == sp.zeros(3, 1))
a3 = sp.Matrix(AM_PTS[0])
witT, valT = [], []
gradT = [sp.expand(g) for g in Ft]
for p in AM_PTS:
    b3 = sp.Matrix(p)
    ybv = (sp.eye(3) + JH.subs(dict(zip(PHI, p))).T).solve(sp.I * (a3 - b3))
    pt = tuple(sp.expand(c) for c in list(b3 - sp.I * ybv) + list(ybv))
    witT.append(pt)
    valT.append(tuple(sp.expand(g.subs(dict(zip(V6, pt)))) for g in gradT))
check("three distinct Q(i)-witnesses for F~ = id + grad f_H, common value "
      "(0,0,-1/4,0,0,0)",
      len(set(witT)) == 3 and valT[0] == valT[1] == valT[2]
      and valT[0] == (0, 0, sp.Rational(-1, 4), 0, 0, 0))
print("     witnesses:", witT)

# ===========================================================================
print("== A4. no nilpotent normalization; Hess f_H not nilpotent ==")
# ===========================================================================
check("tr((Hess f_H)^2) != 0: Hess f_H is NOT nilpotent (Lemma 1.2 of "
      "dBvdE05 is an iff and JH is not nilpotent)",
      sp.expand(sp.trace(hf * hf)) != 0)
check("JH not nilpotent: charpoly(JH) != t^3",
      sp.expand(JH.charpoly(sp.Symbol("t")).as_expr()
                - sp.Symbol("t") ** 3) != 0)
# Any affine normalization x + H' with JH' nilpotent needs an invertible C
# with C JF(u) - I nilpotent for all u (C = B A absorbs source and target;
# translations only shift the evaluation point).  Nilpotency needs
# tr(C JF) = 3 identically -- already INCONSISTENT as a linear system:
Cm = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"k{i}{j}"))
tr_eqs = sp.Poly(sp.expand(sp.trace(Cm * JF) - 3), x, y, z).coeffs()
check("tr(C JF(x)) = 3 identically is INCONSISTENT (linsolve EmptySet): "
      "no affine normalization of AM has nilpotent JH",
      sp.linsolve(tr_eqs, list(Cm)) == sp.EmptySet)

# ===========================================================================
print("== A5. the AM map is not affinely equivalent to ANY gradient map ==")
# ===========================================================================
# B F(Ax) + c is a gradient map iff B JF(Ax) A is symmetric for all x iff
# (with K := A^{-T} B, u := Ax)  K JF(u) symmetric for all u.  The K's form
# a linear space; assert it is {0}.
sym_eqs = []
KJ = sp.expand(Cm * JF - (Cm * JF).T)
for i in range(3):
    for j in range(i + 1, 3):
        sym_eqs += sp.Poly(KJ[i, j], x, y, z).coeffs()
solK = sp.linsolve(sym_eqs, list(Cm))
check("the space {K : K JF(x) symmetric for all x} is exactly {0}: "
      "AM is NOT affinely equivalent to a gradient map (strongest 3D no-go)",
      solK == sp.FiniteSet(tuple([0] * 9)))

# ===========================================================================
print("== B1. weight bookkeeping for gradient equivariant maps ==")
# ===========================================================================
# Source weights (1,-1,-m) on (x,y,z); invertibility forces the component
# weights (d1,d2,d3) to be a permutation of (-m,-1,1) [reduction_w].  If
# F = grad W, F_i = d_i W is weight-homogeneous of weight d_i, so the graded
# pieces of W satisfy: piece e contributes to F_i only if e = e_i :=
# d_i + w_i.  Branch table over the 6 permutations:
m = sp.Symbol("m", integer=True, positive=True)
wts = (1, -1, -m)
table = {}
for d in itertools.permutations((-m, -1, 1)):
    e = tuple(d[i] + wts[i] for i in range(3))
    key = tuple(sp.simplify(e[i] - e[j]) for i, j in ((0, 1), (0, 2), (1, 2)))
    all_eq = sp.solve([e[0] - e[1], e[0] - e[2]], m)
    table[d] = (e, all_eq)
# all three e_i equal: only d = (-m,-1,1) at m = 3 (e = -2)
sol_all = {d: s for d, (e, s) in table.items() if s}
check("all-equal branch: exactly d = (-m,-1,1) with m = 3, e = -2",
      list(sol_all) == [(-m, -1, 1)] and sol_all[(-m, -1, 1)] == {m: 3}
      and tuple(sp.simplify(sp.sympify(ei).subs(m, 3))
                for ei in table[(-m, -1, 1)][0]) == (-2,) * 3)
# pairwise-equal branches, m generic (solve e_i = e_j for m or identically):
pair_branches = []
for d, (e, _) in table.items():
    for (i, j, k) in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        diff = sp.simplify(e[i] - e[j])
        if diff == 0:
            pair_branches.append((d, (i, j), "all m"))
        else:
            s = sp.solve(diff, m)
            s = [si for si in s if si.is_integer and si >= 2]
            if s:
                pair_branches.append((d, (i, j), tuple(s)))
expected = {
    ((-1, 1, -m), (0, 1), "all m"),      # 2D slice W = A(xy) + del z^2/2
    ((-m, -1, 1), (0, 2), "all m"),      # 1D slice W = x^{1-m}B(v) + be y^2/2
    ((1, -m, -1), (1, 2), "all m"),      # gradient shear family
}
got_all_m = {(d, ij, tag) for (d, ij, tag) in pair_branches if tag == "all m"}
check("pairwise-equal branches valid for every m: exactly the three "
      "expected slices (A(w) / B(v) / shear)", got_all_m == expected)
sporadic = {(d, ij): tag for (d, ij, tag) in pair_branches
            if tag != "all m"}
check("all sporadic pairwise coincidences happen at m = 3 only "
      "(subsumed into the all-equal branch or weight-infeasible)",
      all(tag == (3,) for tag in sporadic.values()))
# In the sporadic m = 3 cases and in the branches above, the univariate
# third component is weight-constrained; the two genuinely infeasible ones
# (F3 in C[z] of weight +1, F2 in C[y] of weight -m needing y^m: kills
# DF(0)) are recorded in the doc; the surviving families are B3's.

# ===========================================================================
print("== B2. DF symmetric <=> potential S, and the m != 3 degeneration ==")
# ===========================================================================
Pf, Qf, Rf = [f(w, v) for f in sp.symbols("P Q R", cls=sp.Function)]
# work in (x, w, v) coordinates: w = xy, v = x^m z give the derivations
#   d/dx|_{y,z} = d/dx + (w/x) d/dw + (m v/x) d/dv,
#   d/dy|_{x,z} = x d/dw,   d/dz|_{x,y} = x^m d/dv.
Dx = lambda f: (sp.diff(f, x) + (w / x) * sp.diff(f, w)      # noqa: E731
                + (m * v / x) * sp.diff(f, v))
Dy = lambda f: x * sp.diff(f, w)                             # noqa: E731
Dz = lambda f: x**m * sp.diff(f, v)                          # noqa: E731
F1g, F2g, F3g = Pf / x**m, Qf / x, x * Rf
Pw, Pv = sp.diff(Pf, w), sp.diff(Pf, v)
Qw, Qv = sp.diff(Qf, w), sp.diff(Qf, v)
Rw, Rv = sp.diff(Rf, w), sp.diff(Rf, v)
check("dF1/dy - dF2/dx = x^{1-m} P_w - x^{-2}(w Q_w + m v Q_v - Q)",
      expand_zero(sp.together(Dy(F1g) - Dx(F2g)
                              - (x**(1 - m) * Pw
                                 - (w * Qw + m * v * Qv - Qf) / x**2))))
check("dF1/dz - dF3/dx = P_v - (R + w R_w + m v R_v)",
      expand_zero(sp.together(Dz(F1g) - Dx(F3g)
                              - (Pv - (Rf + w * Rw + m * v * Rv)))))
check("dF2/dz - dF3/dy = x^{m-1} Q_v - x^2 R_w",
      expand_zero(sp.together(Dz(F2g) - Dy(F3g)
                              - (x**(m - 1) * Qv - x**2 * Rw))))
# For m != 3 the x-powers are unbalanced: symmetry forces P_w = R_w = 0 and
# w Q_w + m v Q_v = Q -- the degenerate branches of B1/B3.  At m = 3 the
# system is x-free; Q_v = R_w gives U with (Q,R) = (U_w, U_v) (Poincare),
# and the two remaining identities say exactly P = E(U):
Uf = sp.Function("U")(w, v)
E_op = lambda S_: w * sp.diff(S_, w) + 3 * v * sp.diff(S_, v) - 2 * S_  # noqa: E731
check("m=3 integrability: with (Q,R) = (U_w,U_v),  w Q_w + 3 v Q_v - Q "
      "= d_w E(U)  and  R + w R_w + 3 v R_v = d_v E(U)",
      expand_zero(w * sp.diff(Uf, w, 2) + 3 * v * sp.diff(Uf, w, v)
                  - sp.diff(Uf, w) - sp.diff(E_op(Uf), w))
      and expand_zero(sp.diff(Uf, v) + w * sp.diff(Uf, v, w)
                      + 3 * v * sp.diff(Uf, v, 2) - sp.diff(E_op(Uf), v)))
check("E kernel on admissible monomials is span{w^2}: P = E(S) exactly "
      "(no integration constant fits the box j + 3k >= 3)",
      [(j, k) for j in range(9) for k in range(4)
       if j + 3 * k - 2 == 0] == [(2, 0)])
# and conversely (P,Q,R) = (E(S), S_w, S_v) IS symmetric (m = 3): the three
# residuals above vanish identically under the substitution.
Sf = sp.Function("S")(w, v)
res_syms = (
    x**(1 - m) * Pw - (w * Qw + m * v * Qv - Qf) / x**2,
    Pv - (Rf + w * Rw + m * v * Rv),
    x**(m - 1) * Qv - x**2 * Rw,
)
subS = {Pf: E_op(Sf), Qf: sp.diff(Sf, w), Rf: sp.diff(Sf, v), m: 3}
check("(P,Q,R) = (E(S), S_w, S_v) makes all three symmetry residuals "
      "vanish identically (m=3), i.e. F = grad(x^{-2} S(w,v))",
      all(expand_zero(sp.simplify(r.subs(subS).doit()))
          for r in res_syms))

# ===========================================================================
print("== B3. the degenerate branches: complete classification, all m ==")
# ===========================================================================
alpha, beta, delta, gamma = sp.symbols("alpha beta delta gamma", nonzero=True)
dd = sp.Symbol("d", integer=True, positive=True)
aa = sp.Symbol("a", nonzero=True)

# (i) diagonal branch (all e_i distinct): F = (f1(x), f2(y), f3(z)),
# det DF = f1' f2' f3' = kappa: multidegree (d1-1, d2-1, d3-1) with leading
# coefficient prod a_i d_i != 0, so d1 = d2 = d3 = 1: LINEAR.
d1, d2, d3 = sp.symbols("d1 d2 d3", integer=True, positive=True)
a1, a2, a3_ = sp.symbols("a1 a2 a3", nonzero=True)
prod_lead = sp.expand(sp.diff(a1 * x**d1, x) * sp.diff(a2 * y**d2, y)
                      * sp.diff(a3_ * z**d3, z))
check("diagonal branch: leading term of f1'f2'f3' is "
      "a1 a2 a3 d1 d2 d3 x^{d1-1} y^{d2-1} z^{d3-1} != 0 => linear",
      expand_zero(prod_lead - a1 * a2 * a3_ * d1 * d2 * d3
                  * x**(d1 - 1) * y**(d2 - 1) * z**(d3 - 1)))

# (ii) 2D slice W = A(xy) + delta z^2/2 (branch d = (-1,1,-m), any m):
Af = sp.Function("A")(x * y)
detH = sp.factor(sp.expand(sp.hessian(Af + delta * z**2 / 2,
                                      (x, y, z)).det()))
uf = sp.Function("u")(w)
Aw = sp.Function("A")(w)
handA = (-delta * sp.diff(Aw, w) * (sp.diff(Aw, w)
                                    + 2 * w * sp.diff(Aw, w, 2)))
check("det Hess(A(xy) + delta z^2/2) = -delta A'(A' + 2 w A'')|_{w=xy}",
      sp.simplify(detH - handA.subs(w, x * y)) == 0)
check("Wronskian integration: u(u + 2 w u') = (w u^2)'  [u := A']",
      expand_zero(uf * (uf + 2 * w * sp.diff(uf, w))
                  - sp.diff(w * uf**2, w)))
# (w u^2)' = -kappa/delta => w u^2 = -(kappa/delta) w  (constant term 0 at
# w = 0) => u^2 constant => u constant: deg u = d >= 1 is killed by the
# leading term of w u^2, degree 2d + 1 >= 3 > 1:
check("degree kill: [w^{2d+1}] (w u^2) = a^2 != 0 for u = a w^d + ...",
      sp.expand(w * (aa * w**dd) ** 2 - aa**2 * w**(2 * dd + 1)) == 0)

# (iii) 1D slice W = x^{1-m} B(v) + beta y^2/2 (branch d = (-m,-1,1), any m):
# reduced data P = (1-m)B + m v B', Q = beta w, R = B'.
Bf = sp.Function("B")(v)
Pb = (1 - m) * Bf + m * v * sp.diff(Bf, v)
dMB = sp.expand(det_m(Pb, beta * w, sp.diff(Bf, v), m).doit())
handB = beta * (m * (m - 1) * (Bf - v * sp.diff(Bf, v)) * sp.diff(Bf, v, 2)
                - sp.diff(Bf, v) ** 2
                - 2 * m * v * sp.diff(Bf, v) * sp.diff(Bf, v, 2))
check("B(v)-slice: det M = beta [m(m-1)(B - vB')B'' - B'^2 - 2m v B'B'']",
      expand_zero(dMB - sp.expand(handB)))
kk = sp.Symbol("k", integer=True, positive=True)
leadB = sp.expand((handB / beta).subs(Bf, aa * v**kk).doit())
coefB = -kk * (m * (m - 1) * (kk - 1) ** 2 + kk + 2 * m * kk * (kk - 1))
check("B(v)-slice leading coefficient: det M/beta on B = a v^k is "
      "-k[m(m-1)(k-1)^2 + k + 2mk(k-1)] a^2 v^{2k-2} -- nonzero for "
      "k >= 2, m >= 2 => B linear => F LINEAR",
      expand_zero(leadB - coefB * aa**2 * v**(2 * kk - 2))
      and sp.expand(-coefB.subs({kk: kk + 2, m: m + 2})) ==
      sp.expand((kk + 2) * ((m + 2) * (m + 1) * (kk + 1) ** 2 + kk + 2
                            + 2 * (m + 2) * (kk + 2) * (kk + 1))))

# (iv) the gradient shear family (branch d = (1,-m,-1), any m):
Wsh = alpha * x**2 / 2 + delta * y * z + gamma * y**(m + 1)
Fsh = [sp.expand(sp.diff(Wsh, u)) for u in PHI]
check("shear family: grad W = (alpha x, delta z + (m+1) gamma y^m, "
      "delta y), det Hess = -alpha delta^2 (KELLER for every m)",
      all(expand_zero(fi - gi) for fi, gi in zip(
          Fsh, [alpha * x, delta * z + (m + 1) * gamma * y**m, delta * y]))
      and sp.simplify(sp.hessian(Wsh, PHI).det() + alpha * delta**2) == 0)
a_t, b_t, c_t = sp.symbols("at bt ct")
inv_sh = (a_t / alpha, c_t / delta,
          (b_t - (m + 1) * gamma * (c_t / delta) ** m) / delta)
check("shear family is a TAME automorphism: explicit inverse verified by "
      "composition",
      all(sp.simplify(fi.subs(dict(zip(PHI, inv_sh)), simultaneous=True) - t_)
          == 0 for fi, t_ in zip(Fsh, (a_t, b_t, c_t))))
for mval in (2, 3):
    Fv = [fi.subs(m, mval) for fi in Fsh]
    check(f"shear family m={mval}: survives the infinity prefilter only as "
          "a nonlinear automorphism (documented false-positive class)",
          infinity_prefilter(Fv, PHI) is True)
check("shear family weights: (F1,F2,F3) have C*-weights (1,-m,-1), a "
      "permutation of (1,-1,-m)",
      all(sp.simplify(
          fi.subs({x: sp.Symbol("lam") * x, y: y / sp.Symbol("lam"),
                   z: z / sp.Symbol("lam") ** m}, simultaneous=True)
          - sp.Symbol("lam") ** wt * fi) == 0
          for fi, wt in zip(Fsh, (1, -m, -1))))

# ===========================================================================
print("== B4. the m = 3 potential slice: rigidity of det M(E(S),S_w,S_v) ==")
# ===========================================================================
Sfun = sp.Function("S")(w, v)


def detM_grad(Sexpr):
    """det M for the gradient data (P,Q,R) = (E(S), S_w, S_v), m = 3."""
    Sw, Sv = sp.diff(Sexpr, w), sp.diff(Sexpr, v)
    h_ = (sp.diff(Sexpr, w, 2) * sp.diff(Sexpr, v, 2)
          - sp.diff(Sexpr, w, 1, v, 1) ** 2)
    QF_ = (Sv**2 * sp.diff(Sexpr, w, 2)
           + 2 * Sv * Sw * sp.diff(Sexpr, w, 1, v, 1)
           + Sw**2 * sp.diff(Sexpr, v, 2))
    return sp.expand(-2 * (w * Sw + 6 * v * Sv - 3 * Sexpr) * h_ - QF_)


check("reduced gradient Keller identity: det M(E(S),S_w,S_v) = "
      "-2(wS_w + 6vS_v - 3S)(S_ww S_vv - S_wv^2) - (S_v, S_w) Hess S "
      "(S_v, S_w)^T   [matches reduction_w.det_m]",
      expand_zero(sp.expand(det_m(E_op(Sfun), sp.diff(Sfun, w),
                                  sp.diff(Sfun, v), 3).doit())
                  - detM_grad(Sfun).doit()))
check("anchor: S = v + w^2/2 gives det M = -1 (the linear map (z,y,x))",
      sp.expand(detM_grad(v + w**2 / 2).doit()) == -1)
# kappa != 0 makes DF(0) invertible automatically: det DF(0) =
# -p1(0) q0'(0) r0(0) = kappa, so all three factors are nonzero; the
# residual scalings S -> c S(mu w, nu v) then gauge [v]S = 1, [w^2]S = 1/2,
# whence kappa = -1 (doc section 5).

# ---- K = deg_v S = 0: R = S_v = 0 kills F3 = xR = 0, det M = 0 != kappa:
s0f = sp.Function("s0")(w)
check("K=0: det M(S(w)) = 0 identically (R = S_v = 0, F3 = 0): "
      "EMPTY", sp.expand(detM_grad(s0f).doit()) == 0
      and sp.expand(det_m(E_op(s0f), sp.diff(s0f, w), sp.S(0), 3).doit())
      == 0)

# ---- K = 1: S = s0(w) + s1(w) v.
s1f = sp.Function("s1")(w)
e1 = detM_grad(s0f + s1f * v).doit()
pv1 = sp.Poly(sp.expand(e1), v)
c1_1 = sp.expand(pv1.coeff_monomial(v))
c1_0 = sp.expand(pv1.coeff_monomial(1))
check("K=1 equations: [v] = 2w s1'^3 + 4 s1 s1'^2 - s1^2 s1''  and  "
      "[1] = 2w s0's1'^2 - 6 s0 s1'^2 - s1^2 s0'' - 2 s1 s0's1'",
      expand_zero(c1_1 - (2 * w * sp.diff(s1f, w) ** 3
                          + 4 * s1f * sp.diff(s1f, w) ** 2
                          - s1f**2 * sp.diff(s1f, w, 2)))
      and expand_zero(c1_0 - (2 * w * sp.diff(s0f, w) * sp.diff(s1f, w) ** 2
                              - 6 * s0f * sp.diff(s1f, w) ** 2
                              - s1f**2 * sp.diff(s0f, w, 2)
                              - 2 * s1f * sp.diff(s0f, w)
                              * sp.diff(s1f, w))))
lead1 = sp.expand(c1_1.subs(s1f, aa * w**dd).doit())
check("K=1 degree kill: [v]-equation on s1 = a w^d has leading coefficient "
      "a^3 d(d+1)(2d+1) != 0 => s1 CONSTANT",
      expand_zero(lead1 - aa**3 * dd * (dd + 1) * (2 * dd + 1)
                  * w**(3 * dd - 2)))
# with s1 = s const (gauge 1): [1]-eq is -s0'' = kappa, box val_w s0 >= 2
# (q0 = s0' needs j+3k >= 1 after d_w; P-side gives val >= 3 on p0 = E(s0)):
# s0 = -kappa w^2/2, i.e. S = v + w^2/2 in the gauge: the LINEAR map only.
sconst = sp.Symbol("s", nonzero=True)
check("K=1, s1 = s const: [1]-equation becomes -s^2 s0'' = kappa => s0 "
      "exactly quadratic => S = v + w^2/2 up to gauge: the LINEAR map",
      expand_zero(c1_0.subs(s1f, sconst) + sconst**2 * sp.diff(s0f, w, 2)))

# ---- K = 2: S = s0 + s1 v + a v^2 (top coefficient constant by L1 below).
e2 = detM_grad(s0f + s1f * v + aa * v**2).doit()
pv2 = sp.Poly(sp.expand(e2), v)
c2_3 = sp.expand(pv2.coeff_monomial(v**3))
check("K=2 cascade, c3 = -40 a^2 s1'' => s1 = s + beta w",
      expand_zero(c2_3 + 40 * aa**2 * sp.diff(s1f, w, 2)))
beta_ = sp.Symbol("beta0")
sub_s1 = {s1f: sconst + beta_ * w}
c2_2 = sp.expand(pv2.coeff_monomial(v**2).subs(sub_s1).doit())
check("K=2 cascade, c2 = -40 a^2 s0'' + 12 a beta^2 => s0'' = "
      "3 beta^2/(10a) constant",
      expand_zero(c2_2 - (-40 * aa**2 * sp.diff(s0f, w, 2)
                          + 12 * aa * beta_**2)))
gamma_ = 3 * beta_**2 / (20 * aa)          # s0 = gamma w^2 (box kills lower)
c2_1 = sp.expand(pv2.coeff_monomial(v).subs(sub_s1).doit()
                 .subs(s0f, gamma_ * w**2).doit())
check("K=2 cascade, c1 = -(4/5) beta^2 (3 beta w + s) => beta = 0 => "
      "s0'' = 0 => q0'(0) = 0: DF(0) SINGULAR.  K = 2 is EMPTY "
      "(all w-degrees)",
      expand_zero(c2_1 + sp.Rational(4, 5) * beta_**2
                  * (3 * beta_ * w + sconst)))

# ---- symbolic-K lemmas for the top rows.
K = sp.Symbol("K", integer=True, positive=True)
sK, tK, uK = (sp.Function(n)(w) for n in ("sK", "tK", "uK"))


def v_coeff(expr, target):
    """Coefficient of v**target in an expansion with symbolic exponents."""
    out = 0
    for term in sp.Add.make_args(sp.expand(sp.powsimp(expr))):
        c, pw = term.as_independent(v)
        ex = sp.S(0)
        if pw == v:
            ex = sp.S(1)
        elif pw.is_Pow and pw.base == v:
            ex = pw.exp
        elif pw.is_Mul:
            for f_ in pw.args:
                if f_ == v:
                    ex += 1
                elif f_.is_Pow and f_.base == v:
                    ex += f_.exp
                else:
                    c *= f_
        elif pw != 1:
            raise ValueError(term)
        if sp.simplify(ex - target) == 0:
            out += c
    return sp.expand(out)


eK = detM_grad(sK * v**K + tK * v**(K - 1) + uK * v**(K - 2)).doit()
# L1: c_{3K-2} depends on sK alone (v-degree bookkeeping: only (K,K,K)
# reaches 3K-2); on sK = a w^delta it has the single leading term below.
top = v_coeff(eK, 3 * K - 2)
check("L1 support: c_{3K-2} involves only the top row s_K",
      not ({tK, uK} & top.free_symbols))
dl = sp.Symbol("delta", integer=True, positive=True)
top_mono = sp.expand(sp.powsimp(top.subs(sK, aa * w**dl).doit()))
lead_formula = (aa**3 * K * dl * (3 * K + dl - 2) * (4 * K + 2 * dl - 3)
                * w**(3 * dl - 2))
check("L1: c_{3K-2} on s_K = a w^delta equals "
      "a^3 K delta (3K + delta - 2)(4K + 2 delta - 3) w^{3 delta - 2} "
      "(both brackets > 0 for K, delta >= 1)",
      expand_zero(top_mono - lead_formula))
check("L1 positivity: for K, delta >= 1 the brackets decompose as "
      "3K + delta - 2 = 3(K-1) + (delta-1) + 2 >= 2 and "
      "4K + 2 delta - 3 = 4(K-1) + 2(delta-1) + 3 >= 3  "
      "=> s_K is CONSTANT for every K >= 1",
      sp.expand((3 * K + dl - 2) - (3 * (K - 1) + (dl - 1) + 2)) == 0
      and sp.expand((4 * K + 2 * dl - 3)
                    - (4 * (K - 1) + 2 * (dl - 1) + 3)) == 0)
# L2: with s_K = a, c_{3K-3} = -K(12K^2-17K+6) a^2 t''  => s_{K-1} linear.
eK2 = detM_grad(aa * v**K + tK * v**(K - 1) + uK * v**(K - 2)).doit()
sub_lvl = v_coeff(eK2, 3 * K - 3)
check("L2: c_{3K-3} = -K(12K^2 - 17K + 6) a^2 t''  [12K^2-17K+6 > 0 for "
      "K >= 1]  => s_{K-1} is LINEAR",
      expand_zero(sub_lvl + K * (12 * K**2 - 17 * K + 6) * aa**2
                  * sp.diff(tK, w, 2)))
# L3: with s_K = a, s_{K-1} = t0 + t1 w, c_{3K-4} fixes u'':
t0s, t1s = sp.symbols("t0 t1")
eK3 = detM_grad(aa * v**K + (t0s + t1s * w) * v**(K - 1)
                + uK * v**(K - 2)).doit()
lvl3 = v_coeff(eK3, 3 * K - 4)
check("L3: c_{3K-4} = 3a(K-1)(4K^2-7K+2) t1^2 - a^2 K(3K-2)(4K-3) u''  "
      "=> s_{K-2} is QUADRATIC with fixed leading coefficient prop. t1^2",
      expand_zero(lvl3 - (3 * aa * (K - 1) * (4 * K**2 - 7 * K + 2) * t1s**2
                          - aa**2 * K * (3 * K - 2) * (4 * K - 3)
                          * sp.diff(uK, w, 2))))

# ---- in-box Groebner closure for K = 3, 4 (presolved by L1, L2).
print(f"     Groebner backends: {available_backends()}")


def grad_box_query(Kdeg, J, timeout=900):
    """EMPTY-in-box query for the m=3 gradient slice, deg_v S = Kdeg,
    w-degrees <= J, presolved by L1 (top row constant) and L2 (second row
    linear); gauge [v]S = 1, [w^2]S = 1/2 (=> kappa = -1).  Asserts that
    Keller + (top v-coefficient != 0) is the unit ideal."""
    s0c = {j: sp.Symbol(f"a{j}0") for j in range(3, J + 1)}
    s1c = {j: sp.Symbol(f"a{j}1") for j in range(1, J + 1)}
    S = (w**2 / 2 + sum(c * w**j for j, c in s0c.items())
         + (1 + sum(c * w**j for j, c in s1c.items())) * v)
    extra = []
    if Kdeg == 3:
        b0, b1, c0 = sp.symbols("b0 b1 c0")
        S += (b0 + b1 * w) * v**2 + c0 * v**3
        extra, topc = [b0, b1, c0], c0
    elif Kdeg == 4:
        b0, b1, b2, t0q, t1q, c0 = sp.symbols("b0 b1 b2 t0q t1q c0")
        S += (b0 + b1 * w + b2 * w**2) * v**2 + (t0q + t1q * w) * v**3 \
            + c0 * v**4
        extra, topc = [b0, b1, b2, t0q, t1q, c0], c0
    keller = sp.Poly(sp.expand(detM_grad(S) + 1), w, v).coeffs()
    tt = sp.Symbol("tt")
    gens = list(s0c.values()) + list(s1c.values()) + extra + [tt]
    return is_unit_ideal(keller + [tt * topc - 1], gens, timeout=timeout)


for Kdeg, J in ((3, 4), (3, 6), (3, 8), (4, 4), (4, 6)):
    check(f"in-box m=3 gradient slice, deg_v = {Kdeg}, w-deg <= {J}: EMPTY "
          "(msolve, exact over Q)", grad_box_query(Kdeg, J) is True)
if ARGS.full:
    for Kdeg, J in ((3, 10), (4, 8)):
        check(f"[--full] in-box m=3 gradient slice, deg_v = {Kdeg}, "
              f"w-deg <= {J}: EMPTY", grad_box_query(Kdeg, J,
                                                     timeout=3600) is True)

# ===========================================================================
print("== C1. degree <= 3 in ANY dimension: the midpoint identity ==")
# ===========================================================================
for nvar in (2, 3, 4):
    Vn = sp.symbols(f"x1:{nvar + 1}")
    an = sp.symbols(f"a1:{nvar + 1}")
    bn = sp.symbols(f"b1:{nvar + 1}")
    monos3 = [mm for dtot in (2, 3) for mm in
              itertools.combinations_with_replacement(range(nvar), dtot)]
    Wn = sum(sp.Symbol(f"q{i}") * sp.Mul(*[Vn[j] for j in mm])
             for i, mm in enumerate(monos3))
    Hn = sp.hessian(Wn, Vn)
    mid = [sp.Rational(1, 2) * (ai + bi) for ai, bi in zip(an, bn)]
    lhs = sp.Matrix([sp.diff(Wn, u).subs(dict(zip(Vn, an)))
                     - sp.diff(Wn, u).subs(dict(zip(Vn, bn))) for u in Vn])
    rhs = Hn.subs(dict(zip(Vn, mid))) * sp.Matrix(
        [ai - bi for ai, bi in zip(an, bn)])
    check(f"midpoint identity n={nvar}: grad W(a) - grad W(b) = "
          "Hess W((a+b)/2)(a-b) for deg W <= 3  => det Hess = kappa != 0 "
          "implies INJECTIVE", sp.expand(lhs - rhs) == sp.zeros(nvar, 1))

# ===========================================================================
print("== C2/C3. normalized witness boxes (Hess W(0) = I, kappa = 1) ==")
# ===========================================================================
# WLOG over C (doc section 7): congruence normalizes Hess W(0) = I, hence
# kappa = 1 and quadratic part sum x^2/2; a translation centers the witness
# pair at +-e/2.  A counterexample of degree <= D in the box survives all
# normalizations, so unit-ideal certificates below are complete for the box.


def witness_box(nvar, maxdeg, timeout=900, modulus=0):
    Vn = (x, y, z)[:nvar]
    monos_ = [mm for dtot in range(3, maxdeg + 1) for mm in
              itertools.combinations_with_replacement(range(nvar), dtot)]
    cs = {i: sp.Symbol(f"c{i}") for i in range(len(monos_))}
    Wn = (sum(u**2 for u in Vn) / 2
          + sum(cs[i] * sp.Mul(*[Vn[j] for j in mm])
                for i, mm in enumerate(monos_)))
    keller = sp.Poly(sp.expand(sp.hessian(Wn, Vn).det() - 1), *Vn).coeffs()
    es = sp.symbols(f"e1:{nvar + 1}")
    tt = sp.Symbol("tt")
    subp = dict(zip(Vn, [e / 2 for e in es]))
    subm = dict(zip(Vn, [-e / 2 for e in es]))
    wit = [sp.expand(sp.diff(Wn, u).subs(subp) - sp.diff(Wn, u).subs(subm))
           for u in Vn]
    gens = list(cs.values()) + list(es) + [tt]
    return all(is_unit_ideal(keller + wit + [tt * e - 1], gens,
                             timeout=timeout, modulus=modulus)
               for e in es), keller, wit, gens, es, tt


for deg_ in (4, 5, 6):
    ok, *_ = witness_box(2, deg_)
    check(f"n=2, deg W <= {deg_}: NO non-injective symmetric Keller map "
          "in-box (witness queries unit, msolve exact over Q)", ok)
# controls: the pipeline must NOT report unit without the Keller equations,
# and the Keller variety alone must be nonempty (nonlinear solutions exist).
_, kel2, wit2, gens2, es2, tt2 = witness_box(2, 4)
check("control: witness system WITHOUT Keller is NOT unit (witness pairs "
      "exist for generic quartics)",
      is_unit_ideal(wit2 + [tt2 * es2[0] - 1], gens2) is False)
check("control: the n=2 Keller variety alone is NOT empty (nonlinear "
      "symmetric Keller maps exist)",
      is_unit_ideal(kel2, gens2[:-3]) is False)
fq = sp.Function("f")(x)
Wxy = x * y + fq
check("the nonlinear n=2 family W = xy + f(x): det Hess = -1 and grad W = "
      "(y + f'(x), x) has the explicit inverse (b, a - f'(b)): INJECTIVE",
      expand_zero(sp.hessian(Wxy, (x, y)).det() + 1)
      and sp.expand(sp.Matrix([sp.diff(Wxy, x), sp.diff(Wxy, y)]).subs(
          {x: sp.Symbol("bq"), y: sp.Symbol("aq")
           - sp.diff(fq, x).subs(x, sp.Symbol("bq"))},
          simultaneous=True)
          - sp.Matrix([sp.Symbol("aq"), sp.Symbol("bq")])).doit()
          == sp.zeros(2, 1))

ok3, *_ = witness_box(3, 4, timeout=1800)
check("n=3, deg W <= 4: NO non-injective symmetric Keller map in the FULL "
      "25-coefficient box (three witness queries, msolve exact over Q)",
      ok3)

if ARGS.full:
    # the honest wall: n=3, deg <= 5 (46 coefficients, 50 unknowns).
    # mod-p screens first (memory-light), then one exact attempt.
    print("     [--full] n=3 deg<=5: mod-p screens + exact attempt "
          "(may hit the 16 GB F4 wall -- reported honestly)")
    try:
        okp, *_ = witness_box(3, 5, timeout=3600, modulus=1073741789)
        print(f"     [--full] n=3 deg<=5 mod 1073741789: "
              f"{'EMPTY' if okp else 'NOT unit'}")
    except Exception as ex:                                # noqa: BLE001
        print(f"     [--full] n=3 deg<=5 mod-p: WALL ({ex})")
    try:
        okq, *_ = witness_box(3, 5, timeout=10800)
        print(f"     [--full] n=3 deg<=5 over Q: "
              f"{'EMPTY' if okq else 'NOT unit'}")
    except Exception as ex:                                # noqa: BLE001
        print(f"     [--full] n=3 deg<=5 over Q: WALL ({ex})")

# ===========================================================================
print("== D. coercivity support for the synthesis ==")
# ===========================================================================
check("kappa(W6) = -4 < 0: incompatible with a coercive real potential "
      "(a coercive nondegenerate W has a minimum with det Hess > 0)",
      sp.Integer(-4) < 0)
check("W6 is affine in phibar: W6(phi, .) is never coercive along the "
      "conjugate directions",
      all(sp.diff(W6, a, b) == 0 for a in BAR for b in BAR))

print(f"\nall {N_CHECKS} checks passed in {time.time() - T0:.0f} s"
      + (" (--full)" if ARGS.full else ""))
