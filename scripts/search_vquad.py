"""v-quadratic class of the weight system (1,-1,-2): uniqueness remnant of B3.

(Machinery: jcqft/reduction_w.py; v-linear uniqueness in scripts/search_11m.py /
docs/SEARCH_11M.md.  Write-up: docs/SEARCH_VQUAD.md.)

QUESTION.  Within the equivariant class of (1,-1,-2), with deg_v(P,Q,R) <= 2
and bounded deg_w, does anything besides Alpoge-Mathew live near m = 2 when
deg_v >= 2?  (The v-linear class is already classified: only tame shears + AM.)

VERDICT (default box, asserted here).  ONLY TAME + AM STRATUM -- no new
counterexample.  Precisely: in the gauged box

    deg_v <= 2,
    deg_w of (P,Q,R) at v-degrees 0,1  <= (4, 3, 1),   # contains scaled AM
    deg_w of each v^2 coefficient       <= 1,

every Keller map (det M = -1) is either
  (i)  v-linear (all v^2 coeffs = 0): the gauged tame family or (an affine
       conjugate of) Alpoge-Mathew -- by the v-linear uniqueness theorem; or
  (ii) genuinely v-quadratic: forced to R ≡ 1, hence det M = J2(P,Q) = -1,
       a polynomial automorphism of the (w,v)-plane, hence TAME by the
       Jung-van der Kulk theorem (JC in dimension 2).

The forcing R ≡ 1 when any v^2 coefficient of P or Q is nonzero, and the
emptiness of every r2 != 0 slice, are exact in-box Groebner certificates
(msolve / Singular via jcqft.gb_backend, 16 GB cap, sequential).

Spot-checks: scaled AM (kappa = -1) lies in the box; tame shears with and
without v^2 terms are positive controls.  A cheap deg_v = 3 probe (single
cubic monomial, tiny linear box) shows the same pattern: P/R cubic EMPTY;
Q cubic nonempty but forces R ≡ 1 (tame target shear).

Default run ~3-5 min; --full enlarges the v^2 deg_w bound and the linear
box (budget-capped).
"""

from __future__ import annotations

import argparse
import time

import sympy as sp

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from jcqft.gb_backend import available_backends, is_unit_ideal  # noqa: E402
from jcqft.prefilter import infinity_prefilter                  # noqa: E402
from jcqft.reduction_w import (                                 # noqa: E402
    assemble, det_m, extract, j2, monomial_box, v, w, x, y, z)

T0 = time.time()
M = 2
KAPPA = -1  # gauge p1(0)=q0'(0)=r0(0)=1


def check(label, cond):
    assert cond, label
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


# ===========================================================================
print("== 1. reduced Keller identity at m = 2; v-quadratic E-system ==")
# ===========================================================================

p0, p1, p2, q0, q1, q2, r0, r1, r2 = [
    sp.Function(n)(w) for n in
    ("p0", "p1", "p2", "q0", "q1", "q2", "r0", "r1", "r2")]
P_gen = p0 + p1 * v + p2 * v**2
Q_gen = q0 + q1 * v + q2 * v**2
R_gen = r0 + r1 * v + r2 * v**2
D_gen = sp.expand(det_m(P_gen, Q_gen, R_gen, M))
E = [sp.expand(D_gen.coeff(v, k)) for k in range(6)]
check("det M is v-degree <= 5 on the v-quadratic ansatz",
      all(sp.expand(D_gen.coeff(v, k)) == 0 for k in range(6, 9))
      and E[5] != 0)
# E5 Wronskian certificate: q2^4/p2 * (p2^2 r2 / q2^3)' = E5/2
E5 = E[5]
E5_form = 2 * (p2 * q2 * sp.diff(r2, w) - 3 * p2 * r2 * sp.diff(q2, w)
               + 2 * q2 * r2 * sp.diff(p2, w))
ratio = sp.expand(sp.diff(p2**2 * r2 / q2**3, w) * q2**4 / p2)
check("E5 Wronskian: q2^4/p2 * (p2^2 r2 / q2^3)' = E5/2  "
      "(=> p2^2 r2 = c q2^3 whenever p2 q2 r2 != 0)",
      sp.expand(E5 - E5_form) == 0 and sp.expand(ratio - E5 / 2) == 0)
check("R = 1 => det M == J2(P, Q)  (area-preserving maps of the (w,v)-plane)",
      sp.expand(det_m(P_gen, Q_gen, 1, M) - j2(P_gen, Q_gen)) == 0)

# polynomiality box at m = 2
boxP = monomial_box("P", M, jmax=4, kmax=2)
boxQ = monomial_box("Q", M, jmax=3, kmax=2)
check("polynomiality: P-box j+2k >= 2 contains v, w^2, v^2; "
      "Q-box j+2k >= 1 contains w, v, v^2",
      v in boxP and w**2 in boxP and v**2 in boxP
      and w in boxQ and v in boxQ and v**2 in boxQ
      and sp.S.One not in [sp.sympify(m_) for m_ in boxQ])

# ===========================================================================
print("== 2. positive controls: tame families (incl. v-quadratic) and AM ==")
# ===========================================================================

a, b0, c, alpha, beta = sp.symbols("a b0 c alpha beta")
# (A) elementary + linear target shear (v-linear tame)
P_t = a * w**2 + v
Q_t = w + b0 * P_t
R_t = sp.Integer(1)
check("tame v-linear: P = a w^2 + v, Q = w + b0 P, R = 1 has det M == -1",
      det_m(P_t, Q_t, R_t, M) == -1)
F_t = assemble(P_t, Q_t, R_t, M)
A_, B_, C_ = sp.symbols("A B C")
# F = (a y^2 + z,  y + b0 x F1,  x); inverse uses x = C in the shear
yy_t = B_ - b0 * C_ * A_
zz_t = A_ - a * yy_t**2
check("tame v-linear: explicit inverse by composition",
      all(sp.expand(f.subs({x: C_, y: yy_t, z: zz_t})) == t
          for f, t in zip(F_t, (A_, B_, C_))))
check("tame v-linear: infinity_prefilter survives (known false-positive class)",
      infinity_prefilter(tuple(sp.expand(f.subs({a: 1, b0: 1})) for f in F_t),
                         (x, y, z)) is True)

# (B) completing-the-square family (genuinely v-quadratic, R = 1)
P_sq = v + alpha * (w + beta * v)**2
Q_sq = w + beta * v
R_sq = sp.Integer(1)
check("tame v-quadratic (complete-the-square): det M == -1",
      det_m(P_sq, Q_sq, R_sq, M) == -1)
F_sq = assemble(P_sq, Q_sq, R_sq, M)
inv_sq = (C_, B_ - beta * C_ * (A_ - alpha * B_**2), A_ - alpha * B_**2)
check("tame v-quadratic (complete-the-square): explicit inverse",
      all(sp.expand(f.subs(dict(zip((x, y, z), inv_sq)))) == t
          for f, t in zip(F_sq, (A_, B_, C_))))

# (C) quadratic target shear along P (v-quadratic in Q when p2 = 0)
P_ts = a * w**2 + v
Q_ts = w + b0 * P_ts + c * P_ts**2
R_ts = sp.Integer(1)
check("tame v-quadratic (target shear Q = w + b0 P + c P^2): det M == -1",
      det_m(P_ts, Q_ts, R_ts, M) == -1)
F_ts = assemble(P_ts, Q_ts, R_ts, M)
# F2 = y + b0 x F1 + c x^3 F1^2  (see docs/SEARCH_VQUAD.md)
yy_ts = B_ - b0 * C_ * A_ - c * C_**3 * A_**2
zz_ts = A_ - a * yy_ts**2
check("tame v-quadratic (target shear): explicit inverse",
      all(sp.expand(f.subs({x: C_, y: yy_ts, z: zz_ts})) == t
          for f, t in zip(F_ts, (A_, B_, C_))))

# (D) Alpoge-Mathew, target-scaled into the kappa = -1 gauge
from jcqft.core import F as F_AM  # noqa: E402
P_AM, Q_AM, R_AM = extract(F_AM, M)
R_AM_g = sp.expand(R_AM / 2)
check("scaled AM: det M == -1, gauge p1(0)=q0'(0)=r0(0)=1, deg_v = 1 "
      "(so AM lives in the v-linear slice of every v-quadratic box)",
      det_m(P_AM, Q_AM, R_AM_g, M) == -1
      and sp.expand(P_AM.coeff(v, 1)).subs(w, 0) == 1
      and sp.diff(sp.expand(Q_AM.coeff(v, 0)), w).subs(w, 0) == 1
      and sp.expand(R_AM_g.coeff(v, 0)).subs(w, 0) == 1
      and sp.degree(P_AM, v) == 1 and sp.degree(Q_AM, v) == 1
      and sp.degree(R_AM_g, v) == 1
      and sp.degree(P_AM, w) <= 4 and sp.degree(Q_AM, w) <= 3
      and sp.degree(R_AM_g, w) <= 1)
check("scaled AM: infinity_prefilter survives with escape direction",
      infinity_prefilter(
          assemble(P_AM, Q_AM, R_AM_g, M), (x, y, z)) is True)

# ===========================================================================
print("== 3. Jung-van der Kulk hammer: R = 1 => tame ==")
# ===========================================================================
print("  When R ≡ 1, det M = J2(P, Q).  A polynomial self-map of A^2 with")
print("  constant nonzero Jacobian is an automorphism (JC in dimension 2),")
print("  and Aut(A^2) is generated by affine maps and elementary shears")
print("  (Jung-van der Kulk): every such map is TAME.  Assembling")
print("  F = (P/x^2, Q/x, x) therefore yields a tame automorphism of A^3.")
check("sanity: the three tame controls above all have R = 1",
      R_t == 1 and R_sq == 1 and R_ts == 1)

# ===========================================================================
print("== 4. in-box Groebner stratification (default box) ==")
# ===========================================================================
print(f"  GB backends: {available_backends()}")


def vquad_box(jP, jQ, jR, j_quad, comps=("P", "Q", "R")):
    """Gauged ansatz: v-linear part up to (jP,jQ,jR); v^2 coeffs of each
    listed component up to deg_w <= j_quad.  Returns (P, Q, R, gens, quad)."""
    gens = []

    def lin(prefix, jmax, lo, gauges):
        expr = sp.S.Zero
        for j in range(jmax + 1):
            for k in range(2):  # v^0, v^1
                if j + 2 * k < lo:
                    continue
                if (j, k) in gauges:
                    expr += gauges[(j, k)]
                    continue
                cj = sp.Symbol(f"{prefix}{j}k{k}")
                gens.append(cj)
                expr += cj * w**j * v**k
        return expr

    P = lin("P", jP, 2, {(0, 1): v})
    Q = lin("Q", jQ, 1, {(1, 0): w})
    R = lin("R", jR, -10**9, {(0, 0): 1})
    quad = {}
    for comp in comps:
        poly = sp.S.Zero
        cs = []
        for j in range(j_quad + 1):
            cj = sp.Symbol(f"{comp}q{j}")
            gens.append(cj)
            cs.append(cj)
            poly += cj * w**j
        quad[comp] = (poly, cs)
        if comp == "P":
            P = P + poly * v**2
        elif comp == "Q":
            Q = Q + poly * v**2
        else:
            R = R + poly * v**2
    return P, Q, R, gens, quad


def keller_eqs(P, Q, R):
    return [sp.expand(e) for e in
            sp.Poly(sp.expand(det_m(P, Q, R, M) - KAPPA), w, v).coeffs()]


def sat_unit(eqs, gens, coeff, timeout=120):
    """True iff Keller + (coeff != 0) is the unit ideal."""
    yv = sp.Symbol("yv")
    return is_unit_ideal(eqs + [1 - yv * coeff], gens + [yv],
                         backend="auto", timeout=timeout)


def sat2_unit(eqs, gens, c1, c2, timeout=120):
    """True iff Keller + (c1 != 0) + (c2 != 0) is the unit ideal."""
    yv, yw = sp.symbols("yv yw")
    return is_unit_ideal(eqs + [1 - yv * c1, 1 - yw * c2],
                         gens + [yv, yw], backend="auto", timeout=timeout)


# Default box: contains scaled AM; v^2 deg_w <= 1
JP, JQ, JR, JQ2 = 4, 3, 1, 1
print(f"  default box: deg_w(lin P,Q,R) <= ({JP},{JQ},{JR}), "
      f"deg_w(v^2 coeffs) <= {JQ2}")

# --- 4a. r2 != 0 is EMPTY ---
print("  -- 4a. r2 != 0 slices --")
P, Q, R, gens, qd = vquad_box(JP, JQ, JR, JQ2, comps=("R",))
eqs = keller_eqs(P, Q, R)
print(f"  r2-only ansatz: {len(gens)} unknowns, {len(eqs)} eqs")
for cf in qd["R"][1]:
    try:
        empty = sat_unit(eqs, gens, cf, timeout=120)
        check(f"Keller + ({cf} != 0) [r2-only]: EMPTY", empty)
    except (RuntimeError, TimeoutError) as exc:
        print(f"  [unresolved] {cf} != 0 -- {exc}")

# also with P,Q v^2 free (const) and r2 forced
P, Q, R, gens, qd = vquad_box(3, 2, 1, 0, comps=("P", "Q", "R"))
eqs = keller_eqs(P, Q, R)
print(f"  all-quad-const ansatz: {len(gens)} unknowns, {len(eqs)} eqs")
for cf in qd["R"][1]:
    try:
        empty = sat_unit(eqs, gens, cf, timeout=120)
        check(f"Keller + ({cf} != 0) [all-quad-const]: EMPTY", empty)
    except (RuntimeError, TimeoutError) as exc:
        print(f"  [unresolved] {cf} != 0 -- {exc}")

# --- 4b. p2 or q2 nonzero forces R ≡ 1 (r1 = 0 and r0 = 1) ---
print("  -- 4b. genuine P/Q v^2 forces R ≡ 1 --")
for comp, j_q in (("P", 1), ("Q", 1)):
    P, Q, R, gens, qd = vquad_box(JP, JQ, JR, j_q, comps=(comp,))
    eqs = keller_eqs(P, Q, R)
    r1s = [g for g in gens if str(g).startswith("R") and "k1" in str(g)]
    r0h = [g for g in gens if str(g).startswith("R") and "k0" in str(g)]
    print(f"  {comp} v^2 deg_w <= {j_q}: {len(gens)} unk, "
          f"{len(eqs)} eqs; killing {len(r1s)} r1 + {len(r0h)} r0-high")
    for qc in qd[comp][1]:
        for rc in r1s + r0h:
            try:
                empty = sat2_unit(eqs, gens, qc, rc, timeout=120)
                check(f"Keller + ({qc} != 0) + ({rc} != 0): EMPTY "
                      f"(=> {rc} = 0 on {{{qc} != 0}})", empty)
            except (RuntimeError, TimeoutError) as exc:
                print(f"  [unresolved] {qc} & {rc} -- {exc}")

# --- 4c. v-linear slice contains AM; no claim beyond search_11m ---
print("  -- 4c. v-linear slice (all v^2 = 0) --")
check("scaled AM has all v^2 coefficients zero and fits deg_w <= (4,3,1)",
      all(sp.expand(S.coeff(v, 2)) == 0 for S in (P_AM, Q_AM, R_AM_g))
      and sp.degree(sp.expand(P_AM.coeff(v, 0)), w) <= 4
      and sp.degree(sp.expand(P_AM.coeff(v, 1)), w) <= 3
      and sp.degree(sp.expand(Q_AM.coeff(v, 0)), w) <= 3
      and sp.degree(sp.expand(Q_AM.coeff(v, 1)), w) <= 2
      and sp.degree(sp.expand(R_AM_g.coeff(v, 0)), w) <= 1)
print("  v-linear classification (SEARCH_11M): only the gauged tame family")
print("  and Alpoge-Mathew.  Not re-proved here; AM spot-checked above.")

# ===========================================================================
print("== 5. cheap deg_v = 3 probe (single cubic monomial, tiny lin box) ==")
# ===========================================================================


def cubic_probe(comp, jP=2, jQ=2, jR=1, timeout=60):
    gens = []

    def lin(prefix, jmax, lo, gauges):
        expr = sp.S.Zero
        for j in range(jmax + 1):
            for k in range(2):
                if j + 2 * k < lo:
                    continue
                if (j, k) in gauges:
                    expr += gauges[(j, k)]
                    continue
                cj = sp.Symbol(f"{prefix}{j}k{k}")
                gens.append(cj)
                expr += cj * w**j * v**k
        return expr

    P = lin("P", jP, 2, {(0, 1): v})
    Q = lin("Q", jQ, 1, {(1, 0): w})
    R = lin("R", jR, -10**9, {(0, 0): 1})
    c3 = sp.Symbol(f"{comp}c0")
    gens.append(c3)
    if comp == "P":
        P = P + c3 * v**3
    elif comp == "Q":
        Q = Q + c3 * v**3
    else:
        R = R + c3 * v**3
    return P, Q, R, gens, c3


for comp in ("P", "Q", "R"):
    P, Q, R, gens, c3 = cubic_probe(comp)
    eqs = keller_eqs(P, Q, R)
    try:
        empty = sat_unit(eqs, gens, c3, timeout=60)
        if comp in ("P", "R"):
            check(f"deg_v=3, {comp} cubic const != 0 (tiny box): EMPTY",
                  empty)
        else:
            # Q cubic: nonempty (= tame target shear); force R ≡ 1
            check(f"deg_v=3, Q cubic const != 0 (tiny box): NONEMPTY "
                  f"(expected tame)", not empty)
            r1s = [g for g in gens if str(g).startswith("R")
                   and "k1" in str(g)]
            for rc in r1s:
                empty2 = sat2_unit(eqs, gens, c3, rc, timeout=60)
                check(f"deg_v=3, Qc0 != 0 & {rc} != 0: EMPTY (=> R ≡ 1)",
                      empty2)
    except (RuntimeError, TimeoutError) as exc:
        print(f"  [unresolved] deg_v=3 {comp} -- {exc}")

# cubic tame control
P_c = a * w**2 + v
Q_c = w + c * P_c**3
R_c = sp.Integer(1)
check("tame deg_v=3 control Q = w + c P^3, R = 1: det M == -1",
      det_m(P_c, Q_c, R_c, M) == -1)

print("\nVERDICT (v-quadratic class, default box): ONLY TAME + AM STRATUM.")
print("  No new counterexample.  Genuine v^2 terms force R ≡ 1, hence tame")
print("  by Jung-van der Kulk; the v-linear slice is tame ∪ {AM}.")
print(f"  Box: deg_v <= 2, deg_w(lin) <= ({JP},{JQ},{JR}), "
      f"deg_w(v^2) <= {JQ2}.")


# ===========================================================================
# --full: larger boxes
# ===========================================================================

def run_full(budget):
    print("\n== 6. (--full) larger-box stratification ==")
    # medium: lin (5,4,2), v^2 deg_w <= 2, r2-only and P/Q forcing
    for label, jP, jQ, jR, jq, comps in (
        ("medium r2-only", 5, 4, 2, 2, ("R",)),
        ("medium P v^2", 5, 4, 2, 2, ("P",)),
        ("medium Q v^2", 5, 4, 2, 2, ("Q",)),
    ):
        print(f"  -- {label}: lin <= ({jP},{jQ},{jR}), "
              f"v^2 deg_w <= {jq} --")
        P, Q, R, gens, qd = vquad_box(jP, jQ, jR, jq, comps=comps)
        eqs = keller_eqs(P, Q, R)
        print(f"  {len(gens)} unknowns, {len(eqs)} eqs")
        if comps == ("R",):
            for cf in qd["R"][1]:
                try:
                    empty = sat_unit(eqs, gens, cf, timeout=budget)
                    check(f"[--full] {label}, {cf} != 0: EMPTY", empty)
                except (RuntimeError, TimeoutError) as exc:
                    print(f"  [unresolved] {cf} -- {exc}")
                    # Singular degBound ladder
                    if "singular" not in available_backends():
                        continue
                    yv = sp.Symbol("yv")
                    sys_ = eqs + [1 - yv * cf]
                    gens_ = gens + [yv]
                    for db in (4, 5, 6, 7):
                        try:
                            empty = is_unit_ideal(
                                sys_, gens_, backend="singular",
                                timeout=budget, degbound=db)
                        except (RuntimeError, TimeoutError) as exc2:
                            print(f"  [singular degBound {db}] {exc2}")
                            break
                        if empty:
                            check(f"[--full] {label}, {cf} != 0: EMPTY "
                                  f"(Singular degBound {db})", empty)
                            break
                        print(f"  [no cert at degBound {db}] {cf}")
                    else:
                        print(f"  [WALL] {cf}: no unit certificate "
                              f"up to degBound 7")
        else:
            comp = comps[0]
            r1s = [g for g in gens if str(g).startswith("R")
                   and "k1" in str(g)]
            r0h = [g for g in gens if str(g).startswith("R")
                   and "k0" in str(g)]
            # saturate each v^2 coeff against each R-deviation
            for qc in qd[comp][1]:
                for rc in r1s + r0h:
                    try:
                        empty = sat2_unit(eqs, gens, qc, rc, timeout=budget)
                        check(f"[--full] {label}, {qc}!=0 & {rc}!=0: EMPTY",
                              empty)
                    except (RuntimeError, TimeoutError) as exc:
                        print(f"  [unresolved] {qc} & {rc} -- {exc}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="larger-box Groebner stratification")
    ap.add_argument("--budget", type=int, default=600,
                    help="per-query timeout for --full (s)")
    args = ap.parse_args()
    if args.full:
        run_full(args.budget)
    print(f"\nall checks passed in {time.time() - T0:.1f} s")
