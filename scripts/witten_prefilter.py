"""Validation of the Witten-index / infinity prefilter (jcqft.prefilter).

Implements and validates docs/SEARCH_STRATEGIES.md §1.3 (item A3 of
docs/OPEN_QUESTIONS.md).  Facts verified by the asserts below:

1. The Alpoge-Mathew counterexample SURVIVES the prefilter, both unweighted
   and with its C* source weights (1,-1,-2) — as it must, being non-proper.
   The unweighted witness direction at infinity is [1:0:0], exactly the
   escape direction of the 2:1 orbits x -> (x, w0/x, v0/x^2) and the unique
   point at infinity (flex) of the wall cubic (docs/POSITIVE_GEOMETRY.md §1).
   The weighted witness lies on {R = 2 - 3w - v = 0} in the invariants
   w = xy, v = x^2 z — the escape condition of the C*-normal form.

2. Proper maps behave as the theory predicts:
   - linear invertible maps are REJECTED (leading forms = the map itself,
     common zero only at 0);
   - a proper nonlinear non-automorphism (power-type map) is REJECTED;
   - NONLINEAR polynomial automorphisms (elementary/tame) SURVIVE — these
     are the unavoidable false positives: an injective map of degree > 1
     must have degenerate leading forms (else Bezout would force generic
     fiber size prod deg F_i > 1), so no leading-form test can reject them.
     Reported honestly below; the prefilter's job is rejecting the
     Bezout-generic bulk of a search space, not automorphisms.

3. Throughput: each call is milliseconds (benchmarked on 200 random cubic
   maps), so thousands of candidates can be screened per minute.

Run: .venv/bin/python scripts/witten_prefilter.py   (~10 s)
"""

from __future__ import annotations

import random
import time

import sympy as sp

from jcqft import F, PHI, SOURCE_WEIGHTS
from jcqft.prefilter import infinity_prefilter, infinity_witness, leading_part

x, y, z = PHI


def is_keller_automorphism_candidate(G):
    """Constant nonzero Jacobian determinant (holds for all our tame maps)."""
    d = sp.Matrix([[sp.diff(g, v) for v in PHI] for g in G]).det()
    return sp.expand(d).is_number and d != 0


def compose(G, H):
    """G after H, expanded."""
    sub = dict(zip(PHI, H))
    return tuple(sp.expand(g.subs(sub, simultaneous=True)) for g in G)


def timed(Fc, weights=None):
    t0 = time.perf_counter()
    verdict = infinity_prefilter(Fc, PHI, weights)
    return verdict, (time.perf_counter() - t0) * 1000


def main():
    rows = []  # (name, expected, verdict, ms, note)

    # ------------------------------------------------------------------
    # (i) the known counterexample must survive
    # ------------------------------------------------------------------
    print("=== (i) Alpoge-Mathew counterexample ===")
    v_un, ms_un = timed(F)
    assert v_un, "counterexample must survive the unweighted prefilter"
    wit, sign = infinity_witness(F, PHI)
    print(f"  unweighted: SURVIVES ({ms_un:.1f} ms), witness direction {wit}")
    # leading forms x^3y^3z, 3x^3y^2z, -x^3z all vanish along the x-axis:
    lead = [leading_part(f, PHI)[0] for f in F]
    assert all(lf.subs({x: 1, y: 0, z: 0}) == 0 for lf in lead)
    assert wit == (1, 0, 0), "witness should be the [1:0:0] escape direction"
    # this is the direction of the 2:1 escape orbits (x, w0/x, v0/x^2),
    # x -> oo, and the flex at infinity of the wall cubic
    # (docs/POSITIVE_GEOMETRY.md §1).
    rows.append(("Alpoge-Mathew (unweighted)", "survive", v_un, ms_un,
                 f"witness [1:0:0] = escape direction"))

    v_w, ms_w = timed(F, SOURCE_WEIGHTS)
    assert v_w, "counterexample must survive the weighted prefilter"
    witw, signw = infinity_witness(F, PHI, SOURCE_WEIGHTS)
    print(f"  weighted (1,-1,-2): SURVIVES ({ms_w:.1f} ms), "
          f"witness {witw}, scaling lam -> {'oo' if signw > 0 else '0'}")
    # the weighted witness must satisfy the normal-form escape condition
    # R(w, v) = 2 - 3w - v = 0 with w = x*y, v = x^2*z:
    wx, wy, wz = witw
    w0, v0 = wx * wy, wx**2 * wz
    assert sp.simplify(2 - 3 * w0 - v0) == 0, \
        "weighted witness must lie on the escape locus {R = 0}"
    print(f"    -> witness has (w, v) = ({w0}, {v0}) with "
          "R = 2 - 3w - v = 0: the C*-normal-form escape condition")
    rows.append(("Alpoge-Mathew (weights 1,-1,-2)", "survive", v_w, ms_w,
                 "witness on {R=0} escape locus"))

    # ------------------------------------------------------------------
    # (ii) known-proper maps
    # ------------------------------------------------------------------
    print("\n=== (ii) known-proper maps ===")
    lin1 = (z, y, 2 * x)                             # DF(0) of the model
    lin2 = (x + y, y + z, z + x)                     # invertible (det 2)
    elem1 = (x + y**2, y + z**2, z)                  # elementary automorphism
    elem2 = (x + (y + z)**3, y, z)                   # elementary automorphism
    tame1 = compose(elem1, compose(elem2, lin2))     # tame compositions
    tame2 = compose(compose(lin1, elem2), (x, y + x * z, z))
    power = (x**3 + y, y**3 + z, z**3 + x)           # proper, NOT injective

    proper_maps = [
        ("linear L = DF(0): (z, y, 2x)", lin1, "reject", None),
        ("linear (x+y, y+z, z+x)", lin2, "reject", None),
        ("elementary (x+y^2, y+z^2, z)", elem1, "survive (false pos.)", None),
        ("elementary (x+(y+z)^3, y, z)", elem2, "survive (false pos.)", None),
        ("tame comp. elem1 o elem2 o lin2", tame1, "survive (false pos.)", None),
        ("tame comp. lin1 o elem2 o shear", tame2, "survive (false pos.)", None),
        ("proper power map (x^3+y, y^3+z, z^3+x)", power, "reject", None),
    ]
    # all except `power` are automorphisms; all seven are proper.
    for name, G, expected, _ in proper_maps:
        if G is not power:
            assert is_keller_automorphism_candidate(G), name
        verdict, ms = timed(G)
        tag = "SURVIVES" if verdict else "rejected"
        print(f"  {name}: {tag} ({ms:.1f} ms)  [expected: {expected}]")
        if expected == "reject":
            assert not verdict, f"{name}: proper map with nondegenerate " \
                                "leading forms must be rejected"
            note = "correctly rejected (proper)"
        else:
            # nonlinear automorphisms MUST survive (Bezout argument):
            assert verdict, f"{name}: a nonlinear automorphism cannot be " \
                            "rejected by any leading-form test"
            note = "false positive, forced by Bezout"
        rows.append((name, expected, verdict, ms, note))

    # also check the weighted variant on the linear maps (still rejected)
    for name, G in [("linear L, weighted", lin1)]:
        verdict, ms = timed(G, SOURCE_WEIGHTS)
        assert not verdict
        print(f"  {name}: rejected ({ms:.1f} ms)")
        rows.append((name, "reject", verdict, ms, "correctly rejected"))

    # ------------------------------------------------------------------
    # (iii) throughput benchmark: random cubic maps
    # ------------------------------------------------------------------
    print("\n=== (iii) throughput on random cubic maps ===")
    rng = random.Random(20260721)
    monos3 = sorted(sp.itermonomials([x, y, z], 3),
                    key=lambda m: sp.default_sort_key(m))

    def random_map():
        return tuple(sp.Add(*[rng.randint(-5, 5) * m for m in monos3])
                     for _ in range(3))

    n_maps, n_survive = 200, 0
    t0 = time.perf_counter()
    for _ in range(n_maps):
        if infinity_prefilter(random_map(), PHI):
            n_survive += 1
    total = time.perf_counter() - t0
    per_map = total / n_maps * 1000
    print(f"  {n_maps} random cubic maps: {n_survive} survive, "
          f"{n_maps - n_survive} rejected as proper "
          f"({per_map:.1f} ms/map, {total:.1f} s total)")
    assert per_map < 1000, "prefilter must run in well under a second per map"
    # generic dense cubics have nondegenerate leading forms, so the vast
    # majority must be rejected — that is the point of the prefilter:
    assert n_survive < n_maps / 4, \
        "prefilter should reject the Bezout-generic bulk of random maps"

    # ------------------------------------------------------------------
    # summary table
    # ------------------------------------------------------------------
    print("\n=== summary ===")
    print(f"  {'map':46s} {'verdict':10s} {'ms':>7s}  note")
    for name, expected, verdict, ms, note in rows:
        tag = "survive" if verdict else "reject"
        print(f"  {name:46s} {tag:10s} {ms:7.1f}  {note}")
    print(f"  random cubic benchmark: {n_survive}/{n_maps} survive, "
          f"{per_map:.1f} ms/map")
    print("\nAll asserts passed.")


if __name__ == "__main__":
    main()
