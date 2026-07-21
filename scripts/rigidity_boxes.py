"""Rigidity of the Alpoge-Mathew solution in larger degree boxes.

(Item A4 of docs/OPEN_QUESTIONS.md = docs/NEW_COUNTEREXAMPLES.md §5 step 2.
Extends scripts/search_counterexamples.py, which analyzed the single box
deg_x P in [2,7], Q in [1,6], R in [0,4]; write-up in
docs/RIGIDITY_AND_PREFILTER.md.)

For each degree box (the base box as a cross-check anchor, then two strictly
larger boxes obtained by raising every max weighted x-degree by 1 and by 2),
this script:

1. builds the linearized 2D-reduced Keller condition
       J2(dA, B) + J2(A, dB) = dkappa R^2 + 2 kappa R dR,
   A = P R^2, B = Q R, around the Alpoge-Mathew point (P0, Q0, R0),
   counts equations vs unknowns and computes the exact kernel (sympy);
2. computes the tangent space of the gauge group (equivariance-preserving
   source/target automorphisms) inside the box.  Tangents are computed
   directly in (P, Q, R) coordinates:
       source torus:   x-scaling  d = w d_w + 2v d_v + (-2, -1, +1),
                       y-scaling  d = w d_w,   z-scaling  d = v d_v,
       source shifts:  y -> y + eps x^{2j-1} z^j  gives  d = v^j d_w,
                       z -> z + eps x^{j-2} y^j   gives  d = w^j d_v,
       target torus:   dP = P | dQ = Q | dR = R,
       target shifts:  a -> a + eps b^{j+2} c^j   gives  dP = Q^{j+2} R^j,
                       b -> b + eps a^i c^{2i-1}  gives  dQ = P^i R^{2i-1};
   every tangent is verified to satisfy the linearized Keller condition
   exactly (assert), and tangents leaving the box are skipped (reported);
3. checks gauge orbit tangent \\subseteq kernel (assert), and extracts the
   non-gauge kernel directions by SVD;
4. tests each non-gauge direction for SECOND-ORDER integrability: with
   res(u0 + t xi) = t^2 q2 + O(t^3) (the linear term vanishes), the
   direction integrates to second order iff q2 is in the image of the
   linearization (numpy lstsq on the exact-support grid residual);
5. pushes every non-gauge direction through NONLINEAR CONTINUATION
   (Gauss-Newton with displacement constraint on the full 2D Keller
   system, finite-difference Jacobians) — the same arbiter used by
   scripts/search_counterexamples.py: convergence to a genuinely nearby
   solution = candidate new family; non-convergence or a jump to a distant
   solution = obstructed direction.

VERDICT per box: "RIGID modulo gauge" iff no non-gauge direction continues
to a local solution.  If any direction did continue, the script would
verify Keller exactly at the continued point and test fiber injectivity
(jcqft.fibers) — that branch is present but expected (and observed) to be
unreached: "still rigid, stronger evidence" is the honest outcome.

Caveats (docs/NEW_COUNTEREXAMPLES.md §3): float continuation gives strong
evidence, not proof; exact certification is OPEN_QUESTIONS.md B4.

Run: .venv/bin/python scripts/rigidity_boxes.py   (~50 s)
"""

from __future__ import annotations

import time

import numpy as np
import sympy as sp

from jcqft.reduction import KAPPA, P0, Q0, R0, j2, v, w

STEP = 0.05          # continuation displacement
TOL = 1e-9           # residual tolerance for a continued solution
JUMP_FACTOR = 5      # |u - u0| > JUMP_FACTOR*STEP  =>  distant solution


def box(dmin, dmax):
    """Monomials w^i v^j with dmin <= i + 2j <= dmax (x-degree of w^i v^j
    is i + 2j; the lower bound encodes realizability of F as a polynomial
    map, cf. scripts/search_counterexamples.py)."""
    return [(i, j) for j in range(dmax // 2 + 1) for i in range(dmax + 1)
            if dmin <= i + 2 * j <= dmax]


# base box of scripts/search_counterexamples.py, then +1 and +2 on every
# max weighted degree
BOXES = [
    ("base  (P:[2,7], Q:[1,6], R:[0,4])", (2, 7), (1, 6), (0, 4)),
    ("box+1 (P:[2,8], Q:[1,7], R:[0,5])", (2, 8), (1, 7), (0, 5)),
    ("box+2 (P:[2,9], Q:[1,8], R:[0,6])", (2, 9), (1, 8), (0, 6)),
]


# ---------------------------------------------------------------------------
# linearized kernel and gauge tangents (exact, sympy)
# ---------------------------------------------------------------------------

def basis_vector(dP, dQ, dR, dkappa, boxes):
    """Coefficient vector in the fixed monomial layout, or None if the
    deformation has support outside the box."""
    BOX_P, BOX_Q, BOX_R = boxes
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


def linearized_system(boxes):
    """Matrix of the linearized reduced Keller condition in the box.
    Returns (M, n_equations, unknowns)."""
    BOX_P, BOX_Q, BOX_R = boxes
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
    eqs = [sp.Eq(cf, 0) for cf in sp.Poly(L, w, v).coeffs()]
    M, _ = sp.linear_eq_to_matrix(eqs, unknowns)
    return M, len(eqs), unknowns


def gauge_generators():
    """Symbolic (label, dP, dQ, dR) tangents of the gauge group; generous
    list — per-box fitting is decided by basis_vector."""
    Pw, Pv = sp.diff(P0, w), sp.diff(P0, v)
    Qw, Qv = sp.diff(Q0, w), sp.diff(Q0, v)
    Rw, Rv = sp.diff(R0, w), sp.diff(R0, v)
    gens = [
        ("x-scaling",
         w * Pw + 2 * v * Pv - 2 * P0,
         w * Qw + 2 * v * Qv - Q0,
         w * Rw + 2 * v * Rv + R0),
        ("y-scaling", w * Pw, w * Qw, w * Rw),
        ("z-scaling", v * Pv, v * Qv, v * Rv),
        ("a-scaling", P0, 0, 0),
        ("b-scaling", 0, Q0, 0),
        ("c-scaling", 0, 0, R0),
    ]
    for jj in range(1, 5):        # y -> y + eps x^{2j-1} z^j : dw = v^j
        gens.append((f"y-shift dw = v^{jj}",
                     v**jj * Pw, v**jj * Qw, v**jj * Rw))
    for jj in range(2, 9):        # z -> z + eps x^{j-2} y^j : dv = w^j
        gens.append((f"z-shift dv = w^{jj}",
                     w**jj * Pv, w**jj * Qv, w**jj * Rv))
    for jj in range(0, 3):        # a -> a + eps b^{j+2} c^j
        gens.append((f"target a += b^{jj + 2} c^{jj}",
                     sp.expand(Q0**(jj + 2) * R0**jj), 0, 0))
    for ii in range(1, 3):        # b -> b + eps a^i c^{2i-1}
        gens.append((f"target b += a^{ii} c^{2 * ii - 1}",
                     0, sp.expand(P0**ii * R0**(2 * ii - 1)), 0))
    return gens


def gauge_tangents(boxes):
    """In-box gauge tangent vectors, each verified to satisfy the
    linearized Keller identity exactly.  Returns (vectors, labels,
    skipped_labels)."""
    tangents, labels, skipped = [], [], []
    for label, dP, dQ, dR in gauge_generators():
        resid = sp.expand(
            j2(sp.expand(dP * R0**2 + 2 * P0 * R0 * dR), Q0 * R0)
            + j2(P0 * R0**2, sp.expand(dQ * R0 + Q0 * dR))
            - 2 * KAPPA * R0 * dR)
        dkappa, rem = sp.div(resid, sp.expand(R0**2), w)
        dkappa = sp.simplify(dkappa)
        assert sp.expand(rem) == 0 and dkappa.is_number, \
            f"gauge generator fails linearized Keller: {label}"
        vec = basis_vector(dP, dQ, dR, dkappa, boxes)
        if vec is None:
            skipped.append(label)
        else:
            tangents.append(vec)
            labels.append(label)
    return tangents, labels, skipped


# ---------------------------------------------------------------------------
# numeric residual of the full nonlinear 2D Keller system, on an exact-
# support coefficient grid (sized so no product is ever truncated)
# ---------------------------------------------------------------------------

class KellerGrid:
    def __init__(self, boxes):
        self.boxes = boxes
        (pmin, pmax), (qmin, qmax), (rmin, rmax) = \
            [(min(i + 2 * j for i, j in bx), max(i + 2 * j for i, j in bx))
             for bx in boxes]
        # max w-exponent of d_w(P R^2) * d_v(Q R) etc.; +2 slack
        self.GW = (pmax + 2 * rmax) + (qmax + rmax) + 2
        self.GV = self.GW // 2 + 2
        self.n = sum(len(bx) for bx in boxes) + 1

    def gmul(self, A, B):
        if np.count_nonzero(A) > np.count_nonzero(B):
            A, B = B, A
        GW, GV = self.GW, self.GV
        out = np.zeros((GW, GV))
        for i, j in np.argwhere(A):
            out[i:, j:] += A[i, j] * B[:GW - i, :GV - j]
        return out

    def gdw(self, A):
        out = np.zeros_like(A)
        for i in range(1, self.GW):
            out[i - 1] = i * A[i]
        return out

    def gdv(self, A):
        out = np.zeros_like(A)
        for j in range(1, self.GV):
            out[:, j - 1] = j * A[:, j]
        return out

    def unpack(self, u):
        grids = []
        k = 0
        for bx in self.boxes:
            g = np.zeros((self.GW, self.GV))
            for i, j in bx:
                g[i, j] = u[k]
                k += 1
            grids.append(g)
        return grids[0], grids[1], grids[2], u[-1]

    def res(self, u):
        Pg, Qg, Rg, kap = self.unpack(np.asarray(u, dtype=float))
        R2 = self.gmul(Rg, Rg)
        A = self.gmul(Pg, R2)
        B = self.gmul(Qg, Rg)
        E = (self.gmul(self.gdw(A), self.gdv(B))
             - self.gmul(self.gdv(A), self.gdw(B)) - kap * R2)
        return E.ravel()

    def jac(self, u, h=1e-7):
        r0 = self.res(u)
        J = np.empty((len(r0), len(u)))
        for k in range(len(u)):
            up = u.copy()
            up[k] += h
            J[:, k] = (self.res(up) - r0) / h
        return J

    def base_point(self):
        def coeffs_of(expr, bx):
            pd = sp.Poly(sp.expand(expr), w, v)
            return [float(pd.coeff_monomial(w ** m[0] * v ** m[1]))
                    for m in bx]
        BOX_P, BOX_Q, BOX_R = self.boxes
        return np.array(coeffs_of(P0, BOX_P) + coeffs_of(Q0, BOX_Q)
                        + coeffs_of(R0, BOX_R) + [float(KAPPA)])


def second_order_obstruction(grid, u0, xi, Lnum):
    """t^2 coefficient q2 of res(u0 + t*xi) (exact: res is a quintic in t,
    recovered by solving a 6-point Vandermonde system), and the relative
    residual of the best solution of  Lnum @ eta = -q2.  Small residual =
    direction integrates to second order; large = obstructed at order 2."""
    ts = np.array([-3.0, -2.0, -1.0, 1.0, 2.0, 3.0])
    V = np.vander(ts, 6, increasing=True)  # columns t^0 .. t^5
    R = np.array([grid.res(u0 + t * xi) for t in ts])
    coeffs = np.linalg.solve(V, R)         # rows: t^0 .. t^5 coefficients
    q2 = coeffs[2]
    nq = np.linalg.norm(q2)
    if nq < 1e-9:
        return q2, 0.0
    eta, *_ = np.linalg.lstsq(Lnum, -q2, rcond=None)
    rel = np.linalg.norm(Lnum @ eta + q2) / nq
    return q2, rel


def continue_direction(grid, u0, xi, max_iter=80, stall=20):
    """Gauss-Newton continuation with displacement constraint
    xi . (u - u0) = STEP.  Returns (status, u, final_res, dist) with
    status in {'local', 'jump', 'obstructed'}."""
    u = u0 + STEP * xi
    best = np.inf
    since_best = 0
    for _ in range(max_iter):
        r = grid.res(u)
        r_aug = np.append(r, xi @ (u - u0) - STEP)
        m = np.max(np.abs(r_aug))
        if m < TOL:
            dist = np.linalg.norm(u - u0)
            return ("jump" if dist > JUMP_FACTOR * STEP else "local",
                    u, m, dist)
        if m < best * 0.9:
            best, since_best = m, 0
        else:
            since_best += 1
            if since_best >= stall:
                break
        J = np.vstack([grid.jac(u), xi])
        du, *_ = np.linalg.lstsq(J, -r_aug, rcond=None)
        u = u + du
    return "obstructed", u, m, np.linalg.norm(u - u0)


# ---------------------------------------------------------------------------
# per-box analysis
# ---------------------------------------------------------------------------

def analyze_box(name, spec_P, spec_Q, spec_R):
    t_start = time.perf_counter()
    boxes = (box(*spec_P), box(*spec_Q), box(*spec_R))
    BOX_P, BOX_Q, BOX_R = boxes
    print(f"\n=== {name} ===")
    print(f"  unknowns: |P|={len(BOX_P)} + |Q|={len(BOX_Q)} + "
          f"|R|={len(BOX_R)} + dkappa = "
          f"{len(BOX_P) + len(BOX_Q) + len(BOX_R) + 1}")

    # 1. linearized kernel (exact)
    M, n_eq, unknowns = linearized_system(boxes)
    rank = M.rank()
    ker = M.nullspace()
    print(f"  linearized Keller system: {n_eq} equations x "
          f"{len(unknowns)} unknowns, rank {rank}, kernel dim {len(ker)}")
    assert rank + len(ker) == len(unknowns)

    # 2. gauge tangents in the box (each satisfies the linearized eq)
    tangents, labels, skipped = gauge_tangents(boxes)
    G = sp.Matrix([t for t in tangents])
    grank = G.rank()
    print(f"  gauge generators fitting the box: {len(tangents)} "
          f"(rank {grank}); skipped (leave box): {len(skipped)}")
    for lb in labels:
        print(f"    + {lb}")
    for lb in skipped:
        print(f"    - skipped: {lb}")

    # 3. containment gauge orbit tangent <= kernel, non-gauge directions.
    # Exact containment check: g in ker(M) <=> M g = 0 (avoids an exact
    # rank computation on the nullspace basis, whose rational entries blow
    # up in the larger boxes).
    for t, lb in zip(tangents, labels):
        assert M * sp.Matrix(t) == sp.zeros(M.rows, 1), \
            f"gauge tangent not in kernel: {lb}"
    print(f"  gauge orbit tangent contained in kernel: True")
    n_non_gauge = len(ker) - grank
    print(f"  non-gauge first-order deformations: {len(ker)} - {grank} = "
          f"{n_non_gauge}")

    Kf = np.array([[float(e) for e in kv] for kv in ker])
    Kf /= np.linalg.norm(Kf, axis=1, keepdims=True)
    Gf = np.array([[float(e) for e in t] for t in tangents])
    # orthonormal basis of the gauge span via SVD (Gf is rank-deficient,
    # so unpivoted QR of Gf.T would not span the column space)
    Ug, sg, _ = np.linalg.svd(Gf.T, full_matrices=False)
    ng = int(np.sum(sg > 1e-10 * sg[0]))
    assert ng == grank, "float gauge rank must match exact rank"
    Qg = Ug[:, :ng]
    K_perp = Kf - (Kf @ Qg) @ Qg.T
    _, S, Vt = np.linalg.svd(K_perp, full_matrices=False)
    svd_rank = int(np.sum(S > 1e-8 * S[0])) if len(S) else 0
    assert svd_rank == n_non_gauge, \
        "SVD rank of gauge-orthogonal kernel must equal dim ker - rank gauge"
    dirs = Vt[:svd_rank]

    # 4-5. second-order test + nonlinear continuation, per direction
    grid = KellerGrid(boxes)
    print(f"  numeric grid {grid.GW} x {grid.GV} "
          f"({grid.GW * grid.GV} residual entries)")
    u0 = grid.base_point()
    assert np.max(np.abs(grid.res(u0))) < 1e-10, "baseline residual nonzero"
    Lnum = grid.jac(u0)
    # sanity: exact kernel vectors annihilate the numeric linearization
    for kv in Kf:
        assert np.linalg.norm(Lnum @ kv) < 1e-4 * max(1, np.linalg.norm(kv))

    n_local = 0
    survivors = []
    for i, xi in enumerate(dirs):
        _, rel2 = second_order_obstruction(grid, u0, xi, Lnum)
        tag2 = ("integrates to 2nd order" if rel2 < 1e-6
                else f"OBSTRUCTED at 2nd order (rel. defect {rel2:.2e})")
        status, u, m, dist = continue_direction(grid, u0, xi)
        if status == "local":
            n_local += 1
            survivors.append((xi, u))
            gauge_part = np.linalg.norm(Qg.T @ (u - u0))
            print(f"    dir {i}: {tag2}; continuation -> LOCAL solution "
                  f"|u-u0|={dist:.4f} (gauge part {gauge_part:.4f}) "
                  "-> candidate new family!")
        elif status == "jump":
            print(f"    dir {i}: {tag2}; continuation jumped to distant "
                  f"solution (|u-u0|={dist:.2f}) -> locally obstructed")
        else:
            print(f"    dir {i}: {tag2}; continuation failed "
                  f"(final |res|={m:.2e}) -> obstructed")

    verdict = "RIGID modulo gauge" if n_local == 0 else \
        f"{n_local} candidate new famil{'y' if n_local == 1 else 'ies'}"
    dt = time.perf_counter() - t_start
    print(f"  VERDICT: {verdict}   [{dt:.1f} s]")
    return {"name": name, "n_eq": n_eq, "n_unk": len(unknowns),
            "rank": rank, "ker": len(ker), "gauge": grank,
            "non_gauge": n_non_gauge, "n_local": n_local,
            "survivors": survivors, "grid": grid, "time": dt}


def pursue_survivor(grid, u0, u):
    """Only reached if a continuation produced a genuinely nearby solution:
    verify the 2D Keller identity exactly at a rationalization of u, and
    test injectivity numerically via the Z2 mechanism / fiber machinery."""
    BOX_P, BOX_Q, BOX_R = grid.boxes
    vals = [sp.nsimplify(val, rational=True, tolerance=1e-8) for val in u]
    k = 0
    polys = []
    for bx in (BOX_P, BOX_Q, BOX_R):
        polys.append(sp.Add(*[vals[k + t] * w ** m[0] * v ** m[1]
                              for t, m in enumerate(bx)]))
        k += len(bx)
    Pn, Qn, Rn = polys
    kap = vals[-1]
    from jcqft.reduction import keller_residual
    resid = keller_residual(Pn, Qn, Rn, kap)
    print(f"    exact Keller residual at rationalized point: "
          f"{'ZERO' if resid == 0 else sp.sstr(resid)[:120]}")
    sols = sp.solve([Qn, Rn], [w, v], dict=True)
    print(f"    common zeros of (Q, R) [2:1 orbits if P != 0]: {sols}")
    for s in sols:
        print(f"      P there = {Pn.subs(s)}")


def main():
    t0 = time.perf_counter()
    print("Rigidity of the Alpoge-Mathew solution in growing degree boxes")
    print("(2D-reduced Keller condition; first-order kernel, gauge quotient,")
    print(" second-order integrability, nonlinear continuation)")

    results = []
    # base box: anchor against scripts/search_counterexamples.py
    # (documented: kernel 15, in-box gauge rank 9 there; our direct gauge
    # generator list may fit MORE tangents in-box, so require kernel == 15
    # and gauge rank >= 9)
    r = analyze_box(*BOXES[0])
    assert r["ker"] == 15, "base box kernel must reproduce documented dim 15"
    assert r["gauge"] >= 9, "base box gauge rank must be >= documented 9"
    results.append(r)

    for spec in BOXES[1:]:
        r = analyze_box(*spec)
        results.append(r)
        for xi, u in r["survivors"]:
            pursue_survivor(r["grid"], r["grid"].base_point(), u)

    # strict growth of the boxes
    assert results[0]["n_unk"] < results[1]["n_unk"] < results[2]["n_unk"]

    print("\n=== summary ===")
    print(f"  {'box':38s} {'eqs':>5s} {'unk':>4s} {'rank':>5s} "
          f"{'ker':>4s} {'gauge':>6s} {'nongauge':>9s} {'families':>9s}"
          f" {'s':>6s}")
    for r in results:
        print(f"  {r['name']:38s} {r['n_eq']:5d} {r['n_unk']:4d} "
              f"{r['rank']:5d} {r['ker']:4d} {r['gauge']:6d} "
              f"{r['non_gauge']:9d} {r['n_local']:9d} {r['time']:6.1f}")
    if all(r["n_local"] == 0 for r in results):
        print("\n  -> RIGID modulo gauge in every box tested: no genuinely")
        print("     new counterexample family branches off the Alpoge-Mathew")
        print("     solution within these enlarged ansatz boxes (numerical")
        print("     continuation evidence; exact certification is")
        print("     OPEN_QUESTIONS.md B4).")
    print(f"\nTotal runtime: {time.perf_counter() - t0:.1f} s")
    print("All asserts passed.")


if __name__ == "__main__":
    main()
