"""Degree-2 box of the (2,-1,-3) search: block elimination + mechanism
emptiness (docs/SEARCH_213.md par. 5.3).

The full degree-2 box (scripts/search_213.py run_deg2: A, B, D generic of
degree <= 2 in (u1, u2, u3) with 8 monomials each, C in span{u3, u3^2},
E in span{1, u1, u1^2}; 29 unknowns, 57 Keller equations) was UNRESOLVED:
sympy too slow, msolve/Singular past the 16 GB cap on all six mechanism
queries (even mod p at 8 F4 threads).  This script shrinks the queries by
exact structural reductions BEFORE calling any GB engine.

REDUCTION 1 -- global A-block elimination (asserted here).  The Keller
numerator is trilinear (degree <= 1 in each of the blocks A | B,C | D,E).
The 57 equations are TRIANGULAR in the A-block: 8 pivot equations have the
form  c*Ai + (terms in already-pivoted Aj and the other blocks),  c a
nonzero integer in {2,...,7}.  Eliminating the A-block is therefore a
GLOBAL polynomial substitution (no denominators, no strata): the Keller
variety projects bijectively onto the variety of the reduced system
(49 equations, 21 unknowns, degree <= 4).  Valid verbatim over any field
of characteristic 0 or p > 7.

REDUCTION 2 -- the saturation condition is automatic (asserted here).
Every term of the reduced determinant

    Delta = p^2 (2 p at J2(bt,et) + bt J2(p at, et) - 3 et J2(p at, bt))

carries at, bt or et as an UNDIFFERENTIATED factor:  Delta = at*X + bt*Y
+ et*Z identically.  At a witness point where two of the three vanish,
Keller (Delta == 1) forces the third to be NONZERO automatically -- the
Rabinowitsch variable r and its equation 1 - r*(third) can be dropped.
Same on the y = 0 stratum: det DF restricted to the orbit {u1 = u2 = 0,
u3 = 1} lies in the ideal (A(0,0,1), C(1), D(0,0,1)) of the coefficient
ring, so with two witness values zero, det DF == 1 forces the third
nonzero.  (Integer coefficients: also valid mod every prime.)

REDUCTION 3 -- linear witness elimination.  The witness equations that
are linear with unit pivots (bt/et values for 2:1, bt for 3:1, C/D values
on y = 0) eliminate 1-2 more unknowns per query, globally.

Net effect per query: 60 equations / 30 unknowns (with r)  -->  49-50
equations / 19-20 unknowns, degree <= 4.  The six queries then go through
jcqft.gb_backend (capped subprocess, JCQFT_GB_MEM_MB = 16 GB default),
strictly ONE at a time: a cheap mod-p screen first (screening evidence
only), then exact runs over Q.

Usage:
    .venv/bin/python scripts/deg2_213_elim.py                # everything
    .venv/bin/python scripts/deg2_213_elim.py --skip-exact   # mod-p only
    .venv/bin/python scripts/deg2_213_elim.py --modp 0       # exact only
Options: --budget (s/query, default 1800), --backend, --threads,
--modp N (primes per query), --modp-original (also screen the unreduced
29-unknown systems, cross-checking the reduction mod p), --seed,
--only SUBSTR (filter queries by label).
"""

from __future__ import annotations

import argparse
import os
import random
import resource
import sys
import time

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from jcqft.gb_backend import available_backends, is_unit_ideal  # noqa: E402
from jcqft.reduction_213 import (                               # noqa: E402
    assemble, chart, delta_chart, det_df, j2, keller_numerator, p, q,
    u1, u2, u3)

T0 = time.time()
r_ = sp.Symbol("r")
MON2 = [u1, u2, u3, u1**2, u1 * u3, u3**2, u1 * u2, u2 * u3]
SUB_Y0 = {u1: 0, u2: 0, u3: 1}


def check(label, cond):
    assert cond, label
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


def peak_child_mb():
    """Max RSS of any finished child so far (MB) -- running upper bound."""
    return resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss // 1024


# ---------------------------------------------------------------------------
# the degree-2 box and its Keller system (identical to run_deg2)
# ---------------------------------------------------------------------------

def gpoly(name, monos):
    return sp.Add(*[sp.Symbol(f"{name}{i}") * m for i, m in enumerate(monos)])


def build_box():
    A_ = 1 + gpoly("A", MON2)
    B_ = 1 + gpoly("B", MON2)
    C_ = gpoly("C", [u3, u3**2])
    D_ = 1 + gpoly("D", MON2)
    E_ = gpoly("E", [1, u1, u1**2])
    at_, bt_, et_ = chart(A_, B_, C_, D_, E_)
    N_ = keller_numerator(at_, bt_, et_, 1)
    keller = sp.Poly(N_, p, q).coeffs()
    a_syms = sorted(A_.free_symbols - {u1, u2, u3}, key=str)
    unknowns = sorted(N_.free_symbols - {p, q}, key=str)
    return dict(A=A_, B=B_, C=C_, D=D_, E=E_, at=at_, bt=bt_, et=et_,
                keller=keller, a_syms=a_syms, unknowns=unknowns)


def wit_at(f, P0, Q0):
    return sp.expand(sp.numer(sp.together(f)).subs({p: P0, q: Q0}))


def mech_queries_full(box):
    """The six ORIGINAL mechanism-emptiness systems with the Rabinowitsch
    saturation variable r (same construction as scripts/search_213.py)."""
    at_, bt_, et_ = box["at"], box["bt"], box["et"]
    A_, C_, D_ = box["A"], box["C"], box["D"]
    keller = box["keller"]
    return {
        "2:1 @ (p,q)=(1,1)": keller + [wit_at(bt_, 1, 1), wit_at(et_, 1, 1),
                                       1 - r_ * wit_at(at_, 1, 1)],
        "2:1 @ (p,q)=(1,0)": keller + [wit_at(bt_, 1, 0), wit_at(et_, 1, 0),
                                       1 - r_ * wit_at(at_, 1, 0)],
        "3:1 @ (p,q)=(1,1)": keller + [wit_at(at_, 1, 1), wit_at(bt_, 1, 1),
                                       1 - r_ * wit_at(et_, 1, 1)],
        "3:1 @ (p,q)=(1,0)": keller + [wit_at(at_, 1, 0), wit_at(bt_, 1, 0),
                                       1 - r_ * wit_at(et_, 1, 0)],
        "2:1 @ y=0, u3=1": keller + [C_.subs(SUB_Y0), D_.subs(SUB_Y0),
                                     1 - r_ * A_.subs(SUB_Y0)],
        "3:1 @ y=0, u3=1": keller + [A_.subs(SUB_Y0), C_.subs(SUB_Y0),
                                     1 - r_ * D_.subs(SUB_Y0)],
    }


# ---------------------------------------------------------------------------
# REDUCTION 1: global triangular elimination of the A-block
# ---------------------------------------------------------------------------

def a_coeff_dict(e, a_syms):
    """{Ai (or None): coefficient} for an equation of degree <= 1 in the
    A-block; asserts that degree bound."""
    P = sp.Poly(e, *a_syms)
    out = {}
    for mono, cf in zip(P.monoms(), P.coeffs()):
        vs = [a_syms[i] for i, ex in enumerate(mono) if ex]
        assert len(vs) <= 1, "equation not linear in the A-block"
        out[vs[0] if vs else None] = cf
    return out


def eliminate_a_block(keller, a_syms):
    """Solve the 57 Keller equations for the A-block by triangular pivots
    with CONSTANT nonzero coefficients (globally valid -- polynomial
    substitution, no strata).  Returns (sol_a, reduced, pivot_coeffs)."""
    sol_a: dict = {}
    pivot_coeffs = []
    remaining = list(keller)
    while len(sol_a) < len(a_syms):
        found = None
        for idx, e in enumerate(remaining):
            lin = a_coeff_dict(e, a_syms)
            for a, cf in lin.items():
                if a is None or a in sol_a:
                    continue
                if sp.expand(cf).is_number and cf != 0 and all(
                        k in sol_a for k in lin if k not in (None, a)):
                    found = (idx, a, cf)
                    break
            if found:
                break
        assert found, "no triangular constant pivot left for the A-block"
        idx, a, cf = found
        e = remaining.pop(idx)
        rest = sp.expand((e - cf * a).subs(sol_a))
        assert not rest.free_symbols & set(a_syms), \
            "pivot rest still mentions unsolved A-coefficients"
        sol_a[a] = sp.expand(-rest / cf)
        pivot_coeffs.append(sp.Integer(cf))
        # exactness of the step: the pivot equation vanishes identically
        # after the substitution (it FORCES a = sol_a[a] given the earlier
        # pivots, because its a-coefficient is the nonzero constant cf)
        assert sp.expand(e.subs(sol_a)) == 0, "pivot not solved exactly"
    reduced = []
    for e in remaining:
        f = sp.expand(e.subs(sol_a))
        assert not f.free_symbols & set(a_syms)
        if f != 0:
            reduced.append(f)
    return sol_a, reduced, pivot_coeffs


# ---------------------------------------------------------------------------
# REDUCTION 2: the saturation conditions are automatic
# ---------------------------------------------------------------------------

def assert_saturation_free(box):
    """(i) Delta = at*X + bt*Y + et*Z as an identity of differential
    polynomials in UNDETERMINED functions: at a chart witness point where
    two of (at, bt, et) vanish, Delta == 1 forces the third nonzero.
    (ii) det DF on the y=0 orbit {u1=u2=0, u3=1} vanishes identically on
    {wA = wC = wD = 0}, i.e. lies in the ideal (wA, wC, wD) of the linear
    witness values: same conclusion on the y=0 stratum.  Both statements
    have integer coefficients, so they hold over Q and over every GF(p)."""
    at, bt, et = [sp.Function(n)(p, q) for n in ("at", "bt", "et")]
    X = p**2 * 2 * p * j2(bt, et)
    Y = p**2 * j2(p * at, et)
    Z = -3 * p**2 * j2(p * at, bt)
    check("Delta == at*X + bt*Y + et*Z (undetermined functions) -- chart "
          "saturation conditions are automatic given Keller",
          sp.simplify(delta_chart(at, bt, et)
                      - (at * X + bt * Y + et * Z)) == 0)

    lam, x_, y_, z_ = sp.symbols("lam x y z")
    F = assemble(box["A"], box["B"], box["C"], box["D"], box["E"])
    g = sp.expand(sp.cancel(det_df(F).subs({x_: lam**2, y_: 0,
                                            z_: lam**-3})))
    wA = box["A"].subs(SUB_Y0)
    wC = box["C"].subs(SUB_Y0)
    wD = box["D"].subs(SUB_Y0)
    kill = sp.solve([wA, wC, wD],
                    [sp.Symbol("A5"), sp.Symbol("C1"), sp.Symbol("D5")])
    check("det DF on the y=0, u3=1 orbit: lam-free and in the ideal "
          "(A(0,0,1), C(1), D(0,0,1)) -- y=0 saturation conditions are "
          "automatic given Keller",
          not g.has(lam) and sp.expand(g.subs(kill)) == 0)


# ---------------------------------------------------------------------------
# REDUCTION 3: linear witness elimination (unit pivots, global)
# ---------------------------------------------------------------------------

def lin_eliminate(eqs, wits, gens):
    """Eliminate one unknown per LINEAR witness equation (unit pivot,
    exact global substitution).  Nonlinear witnesses are kept as
    equations.  Returns (new_eqs, new_gens, n_eliminated)."""
    sol: dict = {}
    kept = []
    for w in wits:
        w = sp.expand(w.subs(sol))
        P = sp.Poly(w, *gens)
        if sp.total_degree(P) != 1:
            kept.append(w)
            continue
        # pick a variable with coefficient +-1 (unit over Z: valid mod p)
        pivot = None
        for mono, cf in zip(P.monoms(), P.coeffs()):
            vs = [gens[i] for i, ex in enumerate(mono) if ex]
            if vs and cf in (1, -1):
                pivot = (vs[0], cf)
                break
        assert pivot, f"linear witness without unit pivot: {w}"
        v, cf = pivot
        sol[v] = sp.expand(v - w / cf)
        assert sp.expand(w.subs({v: sol[v]})) == 0
        sol = {k: sp.expand(val.subs({v: sol[v]})) for k, val in sol.items()}
    new_eqs = [sp.expand(e.subs(sol)) for e in eqs]
    new_eqs = [e for e in new_eqs if e != 0] + [sp.expand(k.subs(sol))
                                                for k in kept]
    new_gens = [g for g in gens if g not in sol]
    assert not set().union(*[e.free_symbols for e in new_eqs]) & set(sol)
    return new_eqs, new_gens, len(sol)


def mech_queries_reduced(box, sol_a, reduced, bcde):
    """The six queries after all three reductions: saturation-free, on the
    reduced Keller system, with linear witnesses eliminated.  Returns
    {label: (eqs, gens)}."""
    at_, bt_, et_ = box["at"], box["bt"], box["et"]
    A_, C_, D_ = box["A"], box["C"], box["D"]

    def s(f):
        return sp.expand(f.subs(sol_a))

    raw = {
        "2:1 @ (p,q)=(1,1)": [wit_at(bt_, 1, 1), wit_at(et_, 1, 1)],
        "2:1 @ (p,q)=(1,0)": [wit_at(bt_, 1, 0), wit_at(et_, 1, 0)],
        "3:1 @ (p,q)=(1,1)": [wit_at(bt_, 1, 1), s(wit_at(at_, 1, 1))],
        "3:1 @ (p,q)=(1,0)": [wit_at(bt_, 1, 0), s(wit_at(at_, 1, 0))],
        "2:1 @ y=0, u3=1": [C_.subs(SUB_Y0), D_.subs(SUB_Y0)],
        "3:1 @ y=0, u3=1": [C_.subs(SUB_Y0), s(A_.subs(SUB_Y0))],
    }
    out = {}
    for label, wits in raw.items():
        assert not set().union(*[sp.sympify(w).free_symbols for w in wits]) \
            & set(sol_a), f"witness of {label} still mentions the A-block"
        eqs, gens, nel = lin_eliminate(reduced, wits, bcde)
        out[label] = (eqs, gens)
        print(f"    {label}: {len(eqs)} equations, {len(gens)} unknowns "
              f"({nel} eliminated by linear witnesses)")
    return out


# ---------------------------------------------------------------------------
# GB runs (all through the capped subprocess machinery, one at a time)
# ---------------------------------------------------------------------------

def run_query(label, sys_, gens, backend, timeout, threads, modulus=0):
    """One is_unit_ideal call; returns True / False / None (unresolved)."""
    tag = f"mod {modulus}" if modulus else "over Q"
    t1 = time.time()
    try:
        empty = is_unit_ideal(sys_, gens, backend=backend, timeout=timeout,
                              threads=threads, modulus=modulus)
    except TimeoutError:
        print(f"  [timeout] {label} ({tag}, > {timeout} s) -- UNRESOLVED")
        return None
    except RuntimeError as exc:
        print(f"  [failed] {label} ({tag}) -- UNRESOLVED "
              f"({str(exc)[:120]})")
        return None
    dt = time.time() - t1
    print(f"  [{'EMPTY' if empty else 'NONEMPTY?'}] {label} ({tag}, "
          f"{dt:.1f} s, peak child RSS so far <= {peak_child_mb()} MB)")
    return empty


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=int, default=1800,
                    help="per-query timeout for the exact runs (s)")
    ap.add_argument("--backend", default="msolve",
                    help="msolve | singular | sympy | auto | all")
    ap.add_argument("--threads", type=int, default=4,
                    help="msolve F4 threads (memory grows with threads)")
    ap.add_argument("--modp", type=int, default=3,
                    help="number of random ~30-bit primes for the screen "
                         "(0 disables)")
    ap.add_argument("--modp-original", action="store_true",
                    help="also screen the unreduced 29-unknown systems "
                         "mod p (cross-checks the reduction)")
    ap.add_argument("--seed", type=int, default=213,
                    help="seed for the random screening primes")
    ap.add_argument("--skip-exact", action="store_true",
                    help="only the mod-p screen (cheap)")
    ap.add_argument("--only", default="",
                    help="run only queries whose label contains this")
    args = ap.parse_args()

    print("== degree-2 box: block elimination "
          f"(backends: {available_backends()}) ==")
    box = build_box()
    keller, a_syms = box["keller"], box["a_syms"]
    unknowns = box["unknowns"]
    bcde = [s for s in unknowns if s not in a_syms]
    check(f"box built: {len(unknowns)} unknowns "
          f"({len(a_syms)} A-block + {len(bcde)} B,C,D,E), "
          f"{len(keller)} Keller equations", len(unknowns) == 29
          and len(a_syms) == 8 and len(keller) == 57)
    for e in keller:
        a_coeff_dict(e, a_syms)          # asserts linearity in the A-block
    check("every Keller equation has degree <= 1 in the A-block", True)

    sol_a, reduced, pivot_coeffs = eliminate_a_block(keller, a_syms)
    check("A-block eliminated GLOBALLY by triangular constant pivots "
          f"{sorted(map(int, pivot_coeffs))} (valid over Q and mod p > 7)",
          set(sol_a) == set(a_syms)
          and all(2 <= c <= 7 for c in map(abs, pivot_coeffs)))
    degs = [sp.total_degree(sp.Poly(f, *bcde)) for f in reduced]
    check(f"reduced Keller system: {len(reduced)} equations, {len(bcde)} "
          f"unknowns, degree <= {max(degs)}", max(degs) <= 4)
    # independent spot check of the equivalence at random rational points:
    # a := sol_a(b,c,d,e) turns the 8 pivot equations into 0 == 0 and maps
    # the remaining original residuals onto the reduced residuals.
    rng = random.Random(args.seed)
    for _ in range(3):
        pt = {s: sp.Rational(rng.randint(-9, 9), rng.randint(1, 7))
              for s in bcde}
        pt.update({a: sol_a[a].subs(pt) for a in a_syms})
        orig_res = {sp.expand(e.subs(pt)) for e in keller}
        red_res = {sp.expand(f.subs(pt)) for f in reduced}
        assert red_res <= orig_res | {sp.S.Zero} and \
            orig_res <= red_res | {sp.S.Zero}
    check("random-point spot check: original and reduced Keller systems "
          "have the same residuals under a = sol_a(b,c,d,e)", True)

    assert_saturation_free(box)
    check("witness data bt, et contain no A-coefficients",
          not (box["bt"].free_symbols | box["et"].free_symbols)
          & set(a_syms))

    print("  -- final reduced queries (saturation-free) --")
    queries_red = mech_queries_reduced(box, sol_a, reduced, bcde)
    queries_full = mech_queries_full(box)
    gens_full = unknowns + [r_]

    if args.only:
        queries_red = {k: v for k, v in queries_red.items()
                       if args.only in k}
        queries_full = {k: v for k, v in queries_full.items()
                        if args.only in k}

    # ---- mod-p screen (SCREENING EVIDENCE ONLY, not a proof over Q) ----
    screen: dict = {}
    if args.modp > 0:
        primes = []
        while len(primes) < args.modp:
            pr = sp.randprime(2**29, 2**30)
            if pr not in primes:
                primes.append(pr)
        print(f"\n== mod-p screen: primes {primes} (seed {args.seed}; "
              "evidence only, NOT exact over Q) ==")
        for label, (sys_, gens) in queries_red.items():
            votes = []
            for pr in primes:
                v = run_query(f"{label} [reduced]", sys_, gens,
                              args.backend, args.budget, args.threads,
                              modulus=pr)
                if args.modp_original:
                    w = run_query(f"{label} [29-unknown original]",
                                  queries_full[label], gens_full,
                                  args.backend, args.budget, args.threads,
                                  modulus=pr)
                    assert v is None or w is None or v == w, \
                        "original and reduced disagree mod p -- " \
                        "reduction bug"
                votes.append(v)
            screen[label] = votes
            if votes and all(v is True for v in votes):
                print(f"  => {label}: EMPTY mod all {len(primes)} primes "
                      "(strong screening evidence)")

    # ---- exact runs over Q, on the reduced systems, strictly serial ----
    verdict: dict = {}
    if not args.skip_exact:
        print(f"\n== exact runs over Q (reduced systems; backend "
              f"{args.backend!r}, {args.threads} threads, budget "
              f"{args.budget} s/query, ONE capped subprocess at a time) ==")
        for label, (sys_, gens) in queries_red.items():
            verdict[label] = run_query(label, sys_, gens, args.backend,
                                       args.budget, args.threads)
            if verdict[label] is True:
                check(f"deg2 {label}: EMPTY, exact over Q (GB == [1] on "
                      "the reduced query == the original query by the "
                      "asserted global reductions)", True)

    # ---- summary ----
    print("\n== SUMMARY ==")
    exact_ok = [k for k, v in verdict.items() if v is True]
    exact_open = [k for k, v in verdict.items() if v is None]
    nonempty = [k for k, v in verdict.items() if v is False]
    assert not nonempty, f"NONEMPTY mechanism locus (counterexample?!) " \
        f"in {nonempty} -- investigate immediately"
    for k in exact_ok:
        print(f"  EXACT over Q: {k} EMPTY")
    for k in exact_open:
        ev = screen.get(k, [])
        print(f"  OPEN over Q: {k}"
              + (f" (mod-p screen: {sum(v is True for v in ev)}/{len(ev)} "
                 "primes say empty)" if ev else ""))
    if not args.skip_exact and not exact_open and len(exact_ok) == 6:
        print("\nDEG2 VERDICT: both orbifold mechanisms (2:1 and 3:1) are"
              "\nEMPTY in the FULL degree-2 box -- exact over Q, via the"
              "\nglobal A-block elimination and the saturation-free"
              "\nreformulation (29+1 -> 19-20 unknowns per query).")
    print(f"\ndone in {time.time() - T0:.1f} s "
          f"(peak child RSS <= {peak_child_mb()} MB)")


if __name__ == "__main__":
    main()
