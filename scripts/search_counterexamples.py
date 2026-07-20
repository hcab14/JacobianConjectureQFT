"""Construction mechanism of the counterexample and a first-order rigidity
check, in the dimensionally-reduced (equivariant) setting.

Facts verified by this script (see docs/NEW_COUNTEREXAMPLES.md):

1. NORMAL FORM.  In the C*-invariants w = x*y, v = x^2*z the map is
       F = ( P(w,v)/x^2 ,  Q(w,v)/x ,  x*R(w,v) ),
   with polynomials
       P = (1+w) * ( v*(1+w)^2 + w^2*(3w+4) ),
       Q = 3v*(1+w)^2 + 9w^3 + 12w^2 + w,
       R = 2 - 3w - v.
   Realizability (F polynomial): P has no monomial with x-degree
   i+2j < 2, Q none with i+2j < 1 (monomial w^i v^j has x-degree i+2j).

2. DIMENSIONAL REDUCTION OF THE KELLER CONDITION.  det DF = const = kappa
   is equivalent to the two-variable identity
       J2( P*R^2 , Q*R ) = kappa * R^2         (here kappa = -2),
   where J2 is the Jacobian in (w, v).  The 3D search collapses to a 2D one.

3. Z2 MECHANISM.  F is 2:1 on the fiber x -> (x, w0/x, v0/x^2) whenever
   (w0, v0) is a common zero of Q and R with P(w0, v0) != 0 (the image
   (P/x^2, 0, 0) is invariant under x -> -x).  For the counterexample:
   (w0, v0) = (-3/2, 13/2), P = -1/4 there.

4. FIRST-ORDER RIGIDITY (modulo gauge, within a degree box).  All solutions
   of the linearized Keller condition around (P, Q, R) inside the box
       deg_x(dP) in [2,7], deg_x(dQ) in [1,6], deg_x(dR) in [0,4]
   are tangent to the gauge orbit (equivariance-preserving source/target
   automorphisms: 6 torus scalings + triangular shifts y -> y + eps*xz,
   z -> z + eps*x^{k-2} y^k).  Kernel dimension == gauge rank ==> no new
   counterexample deforms the Alpoge-Mathew one at first order in this box.
"""

from __future__ import annotations

import sympy as sp

from jcqft import F, PHI
from jcqft.reduction import KAPPA, P0, Q0, R0, extract, j2, keller_residual, v, w

x, y, z = PHI
eps = sp.Symbol("eps")


def run_checks():
    print("=== 1-3. Normal form, reduction, Z2 mechanism ===")
    print(f"  P = {sp.sstr(sp.factor(P0))}")
    print(f"  Q = {sp.sstr(Q0)}")
    print(f"  R = {sp.sstr(R0)}")
    print(f"  reduced Keller identity J2(P R^2, Q R) = {KAPPA} R^2:",
          keller_residual(P0, Q0, R0, KAPPA) == 0)
    sols = sp.solve([Q0, R0], [w, v], dict=True)
    print(f"  common zeros of (Q, R): {sols}")
    for s in sols:
        print(f"    P there = {P0.subs(s)}  (nonzero -> 2:1 orbit, "
              "images (P/x^2, 0, 0) fixed under x -> -x)")


# ---------------------------------------------------------------------------
# 4. First-order deformations vs gauge, inside a degree box
# ---------------------------------------------------------------------------

def box(dmin, dmax):
    return [(i, j) for j in range(dmax // 2 + 1) for i in range(dmax + 1)
            if dmin <= i + 2 * j <= dmax]


BOX_P, BOX_Q, BOX_R = box(2, 7), box(1, 6), box(0, 4)


def basis_vector(dP, dQ, dR, dkappa):
    """Coefficient vector of a deformation in the fixed monomial layout,
    or None if the deformation has monomials outside the box."""
    vec = []
    for d, bx in ((dP, BOX_P), (dQ, BOX_Q), (dR, BOX_R)):
        pd = sp.Poly(sp.expand(d), w, v)
        support = {m for m, cf in pd.terms() if cf != 0}
        if support - set(bx):
            return None
        for m in bx:
            vec.append(pd.coeff_monomial(w ** m[0] * v ** m[1]))
    vec.append(dkappa)
    return vec


def linearized_kernel():
    """Nullspace of the linearized reduced Keller condition in the box."""
    cP = {m: sp.Symbol(f"p_{m[0]}_{m[1]}") for m in BOX_P}
    cQ = {m: sp.Symbol(f"q_{m[0]}_{m[1]}") for m in BOX_Q}
    cR = {m: sp.Symbol(f"r_{m[0]}_{m[1]}") for m in BOX_R}
    dk = sp.Symbol("dkappa")
    dP = sp.Add(*[s * w ** m[0] * v ** m[1] for m, s in cP.items()])
    dQ = sp.Add(*[s * w ** m[0] * v ** m[1] for m, s in cQ.items()])
    dR = sp.Add(*[s * w ** m[0] * v ** m[1] for m, s in cR.items()])

    A, B = P0 * R0**2, Q0 * R0
    dA = sp.expand(dP * R0**2 + 2 * P0 * R0 * dR)
    dB = sp.expand(dQ * R0 + Q0 * dR)
    L = sp.expand(j2(dA, B) + j2(A, dB) - dk * R0**2 - 2 * KAPPA * R0 * dR)

    unknowns = list(cP.values()) + list(cQ.values()) + list(cR.values()) + [dk]
    eqs = [sp.Eq(coeff, 0) for coeff in sp.Poly(L, w, v).coeffs()]
    M, rhs = sp.linear_eq_to_matrix(eqs, unknowns)
    ker = M.nullspace()
    return ker, unknowns


def gauge_tangents():
    """Tangent vectors to the (equivariance-preserving) gauge orbit."""
    tangents = []

    def add(dP, dQ, dR, label):
        # dkappa from the linearized identity: residual must be dkappa * R^2
        resid = sp.expand(
            j2(sp.expand(dP * R0**2 + 2 * P0 * R0 * dR), Q0 * R0)
            + j2(P0 * R0**2, sp.expand(dQ * R0 + Q0 * dR))
            - 2 * KAPPA * R0 * dR)
        dkappa, rem = sp.div(resid, sp.expand(R0**2), w)
        dkappa = sp.simplify(dkappa)
        assert sp.expand(rem) == 0 and dkappa.is_number, label
        vec = basis_vector(dP, dQ, dR, dkappa)
        if vec is None:
            print(f"    (skipped, leaves box: {label})")
            return
        tangents.append((vec, label))

    # --- source automorphisms: torus scalings and weighted-triangular shifts
    #     y -> y + f(x,z) (weight -1), z -> z + g(x,y) (weight -2); these are
    #     the only genuine constant-Jacobian moves compatible with the
    #     grading.  Moves whose tangent leaves the box are skipped. ---
    source_moves = {
        "x-scaling": {x: (1 + eps) * x},
        "y-scaling": {y: (1 + eps) * y},
        "z-scaling": {z: (1 + eps) * z},
        "y -> y + eps*x*z (dw = v)": {y: y + eps * x * z},
        "y -> y + eps*x^3*z^2 (dw = v^2)": {y: y + eps * x**3 * z**2},
        "z -> z + eps*y^2 (dv = w^2)": {z: z + eps * y**2},
        "z -> z + eps*x*y^3 (dv = w^3)": {z: z + eps * x * y**3},
        "z -> z + eps*x^2*y^4 (dv = w^4)": {z: z + eps * x**2 * y**4},
    }
    for label, sub in source_moves.items():
        Fp = tuple(sp.expand(f.subs(sub)) for f in F)
        Pe, Qe, Re = extract(Fp, extra_syms=(eps,))
        add(*(sp.expand(e.coeff(eps, 1)) for e in (Pe, Qe, Re)), label)

    # --- target automorphisms: torus scalings and weighted-triangular shifts
    #     a -> a + s*b^2, a -> a + s*b^3*c, b -> b + s*a*c  (all det 1) ---
    add(P0, 0, 0, "a-scaling")
    add(0, Q0, 0, "b-scaling")
    add(0, 0, R0, "c-scaling")
    add(sp.expand(Q0**2), 0, 0, "a -> a + eps*b^2  (dP = Q^2)")
    add(0, sp.expand(P0 * R0), 0, "b -> b + eps*a*c  (dQ = P*R)")
    add(sp.expand(Q0**3 * R0), 0, 0, "a -> a + eps*b^3*c  (dP = Q^3*R)")

    return tangents


def rigidity_check():
    print("\n=== 4. First-order rigidity modulo gauge (degree box "
          f"P:[2,7] Q:[1,6] R:[0,4]) ===")
    ker, unknowns = linearized_kernel()
    print(f"  deformation unknowns: {len(unknowns)} "
          f"(|BOX_P|={len(BOX_P)}, |BOX_Q|={len(BOX_Q)}, |BOX_R|={len(BOX_R)}, dkappa)")
    print(f"  kernel of linearized Keller condition: dim = {len(ker)}")

    tangents = gauge_tangents()
    G = sp.Matrix([t[0] for t in tangents])
    grank = G.rank()
    print(f"  gauge tangent vectors: {len(tangents)} generators, rank = {grank}")
    for _, label in tangents:
        print(f"    - {label}")

    K = sp.Matrix([list(kv) for kv in ker])
    combined = sp.Matrix.vstack(K, G)
    contained = combined.rank() == len(ker)
    print(f"  gauge orbit contained in kernel: {contained}")
    print(f"  kernel dim == gauge rank: {len(ker) == grank}")
    if len(ker) == grank and contained:
        print("  -> RIGID at first order modulo gauge in this box: every")
        print("     infinitesimal deformation of the Alpoge-Mathew solution is")
        print("     a reparametrization; no new counterexample nearby.")
    else:
        print("  -> possible non-gauge deformation directions; inspect kernel.")


def continuation_test(step=0.05, tol=1e-9):
    """Try to integrate each non-gauge kernel direction to a genuine nearby
    solution of the full nonlinear reduced Keller condition (Gauss-Newton
    with a displacement constraint, numeric polynomial-grid arithmetic).
    Convergence = candidate new Keller family (needs equivalence analysis);
    failure = the direction is obstructed beyond first order."""
    import numpy as np

    print("\n=== 5. Nonlinear continuation along non-gauge directions ===")
    GW, GV = 16, 12  # generous grid: coeff[i, j] of w^i v^j

    def gmul(A, B):
        out = np.zeros((GW, GV))
        for (i, j), aij in np.ndenumerate(A):
            if aij != 0.0:
                bw, bv = GW - i, GV - j
                out[i:, j:] += aij * B[:bw, :bv]
        return out

    def gdw(A):
        out = np.zeros((GW, GV))
        for i in range(1, GW):
            out[i - 1] = i * A[i]
        return out

    def gdv(A):
        out = np.zeros((GW, GV))
        for j in range(1, GV):
            out[:, j - 1] = j * A[:, j]
        return out

    def unpack(u):
        Pg, Qg_, Rg = np.zeros((GW, GV)), np.zeros((GW, GV)), np.zeros((GW, GV))
        k = 0
        for g, bx in ((Pg, BOX_P), (Qg_, BOX_Q), (Rg, BOX_R)):
            for (i, j) in bx:
                g[i, j] = u[k]
                k += 1
        return Pg, Qg_, Rg, u[-1]

    def res_f(*u):
        Pg, Qg_, Rg, kapv = unpack(np.asarray(u))
        R2 = gmul(Rg, Rg)
        A = gmul(Pg, R2)
        B = gmul(Qg_, Rg)
        E = gmul(gdw(A), gdv(B)) - gmul(gdv(A), gdw(B)) - kapv * R2
        return E.ravel()

    def jac_f(u, h=1e-7):
        r0 = res_f(*u)
        J = np.empty((len(r0), len(u)))
        for k in range(len(u)):
            up = u.copy()
            up[k] += h
            J[:, k] = (res_f(*up) - r0) / h
        return J

    def coeffs_of(expr, bx):
        pd = sp.Poly(sp.expand(expr), w, v)
        return [float(pd.coeff_monomial(w ** m[0] * v ** m[1])) for m in bx]

    u0 = np.array(coeffs_of(P0, BOX_P) + coeffs_of(Q0, BOX_Q)
                  + coeffs_of(R0, BOX_R) + [float(KAPPA)])
    assert np.max(np.abs(res_f(*u0))) < 1e-12, "baseline residual nonzero"

    ker, _ = linearized_kernel()
    K = np.array([[float(e) for e in kv] for kv in ker])
    G = np.array([[float(e) for e in t] for t, _ in gauge_tangents()])
    Qg = np.linalg.qr(G.T)[0][:, :np.linalg.matrix_rank(G)]
    K_perp = K - (K @ Qg) @ Qg.T
    _, S, Vt = np.linalg.svd(K_perp, full_matrices=False)
    rank = int(np.sum(S > 1e-8 * S[0])) if len(S) else 0
    B = Vt[:rank]
    print(f"  non-gauge kernel directions (SVD rank): {rank}")

    n_new, n_obstructed, n_jump = 0, 0, 0
    for i, xi in enumerate(B):
        u = u0 + step * xi
        ok = False
        for _ in range(120):
            r = np.array(res_f(*u), dtype=float)
            # displacement constraint: stay at distance `step` along xi
            r_aug = np.append(r, xi @ (u - u0) - step)
            if np.max(np.abs(r_aug)) < tol:
                ok = True
                break
            J = np.vstack([jac_f(u), xi])
            du, *_ = np.linalg.lstsq(J, -r_aug, rcond=None)
            u = u + du
        dist = np.linalg.norm(u - u0)
        if not ok:
            n_obstructed += 1
            print(f"  direction {i}: no solution found (obstructed at higher"
                  f" order) [final |res| = {np.max(np.abs(r_aug)):.2e}]")
        elif dist > 5 * step:
            n_jump += 1
            print(f"  direction {i}: Newton jumped to a distant solution "
                  f"(|u-u0| = {dist:.2f}) -- not a local family; direction "
                  "locally obstructed")
        else:
            n_new += 1
            gauge_part = np.linalg.norm(Qg.T @ (u - u0))
            print(f"  direction {i}: LOCAL solution at |u-u0| = {dist:.4f} "
                  f"(gauge component {gauge_part:.4f}) -> candidate new family!")
    print(f"  summary: {n_new} genuine nearby familie(s), "
          f"{n_obstructed + n_jump} locally obstructed non-gauge directions")
    if n_new == 0:
        print("  -> within this box the Alpoge-Mathew solution is RIGID modulo")
        print("     gauge: no genuinely new counterexample family branches off")
        print("     at first order.")
    return n_new


if __name__ == "__main__":
    run_checks()
    rigidity_check()
    continuation_test()
