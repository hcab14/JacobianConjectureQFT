"""Exact wall structure of the 2-site lattice deformation F_eps at eps = 1/4
on the probed segment T1 -> T3 (docs/CLASSICAL_MAP_INVARIANTS.md SS6.4-6.5;
resolves the "exact discriminant / fold-vs-escape" item of SS7.2).

The map (convention of scripts/lattice_chamber.py):

    F_eps(phi0, phi1) = (F(phi0) + eps (phi1 - phi0),
                         F(phi1) + eps (phi0 - phi1)),   eps = 1/4 exactly,

with F the Alpoge-Mathew map (jcqft.core).  Segment
J(t) = (1-t) T1 + t T3, T1 = ((-1/4,0,0),(0,2,0)), T3 = ((1,0,0),(2,1,1)).

COMPUTED HERE, exact over Q unless labelled otherwise:

1. eps = 0 anchors.  det DF_0 == 4 (fold ideal empty; Keller); the eps=0
   escape locus on the segment is p(J0(t)) p(J1(t)) = 0 with
   p(J0(t)) = 4(5t-1) and p(J1(t)) = 107 t^4 + 42 t^3 - 85 t^2 + 44 t - 4:
   exactly two crossings in (0,1), the product chamber count is 9 -> 3 -> 1.

2. FOLD polynomial f(t) (eps = 1/4).  The 7-equation system
   {F_eps(phi) = J(t), det DF_eps(phi) = 0} in (phi, t) is 0-dimensional of
   degree 516, t is a separating element, and its exact eliminant
   f(t) in Z[t] has degree 516, is squarefree and IRREDUCIBLE over Q
   (multi-prime factor-pattern certificate), with 50 real roots, exactly
   14 of them in (0,1) (msolve exact real-root isolation).

3. Exact real chamber counts.  msolve real-root isolation at rational t
   gives the EXACT count sequence across the 14 fold roots:
       18 16 14 16 14 12 10 8 6 8 6 4 2 4 6,
   every jump is +-2 and brackets exactly one fold root: on this segment
   every real chamber wall is FOLD-type.  This supersedes the
   homotopy-continuation counts of SS6.4 (12 -> 13 -> 14 -> ...), which were
   certified lower bounds only: e.g. at T1 (t=0) the exact fiber is
   14 real / 54 complex-distinct, not 12 / 52, and the "odd jumps"
   12->13->14 were completeness artifacts, not escape walls.

4. ESCAPE polynomial e(t) (eps = 1/4).  Escape values on the segment
   (fiber degree of the specialized ideal < generic 66) are the roots of
       e(t) = t * (13 t - 1) * (3 t + 1) * q_y(t) * q_x(t),
   q_y = 29823777 t^4 + 5199180 t^3 + 713782 t^2 - 246740 t + 12337
   (irreducible, NO real roots),
   q_x = 2841875 t^4 + 125650 t^3 - 1157957 t^2 + 512672 t - 54016
   (irreducible, real roots ~ -0.83010 and 0.15629).
   Each factor is certified over Q by an exact fiber-degree drop
   (66 -> 54 at t=0, 66 -> 53 at t=1/13, 66 -> 65 at t=-1/3, and degree
   260 < 264 = 4 x 66 for the ideal with either irreducible quartic
   adjoined -- control quartic t^4+t+1 gives exactly 264); completeness
   of the factor list is certified mod p for two independent ~30-bit
   primes in the x0/y0 (and, recorded, x1/y1) coordinate directions
   (eliminant leading coefficients; rerun behind --full); the
   z-directions exceeded the compute budget (docs SS6.5).
   e(t) has NO root in the old SS6.4 bracket (27/512, 7/128): the
   HC 12->13 'first crossing' there was a completeness artifact.
   Interior escape values do NOT shift the real chamber count (the count
   dips exactly AT the value: N(1/13) = 13 between plateaus of 16); the
   segment endpoint t = 0 IS a real escape point: N(0) = 14 vs
   N(0+) = 18 - four real solutions at infinity over T1 itself.

5. --full: the mod-p completeness rerun (item 4) and a mod-p prescreen of
   the SS7.2 fold surface W(a,b,c) (J1 frozen at (2,1,1)); the prescreen
   itself dies at the 16 GB gb_backend memory cap (~12 min of F4) and is
   reported as the documented honest wall -- the exact surface stays open.

Usage:  .venv/bin/python scripts/lattice_discriminant.py [--full]
        (default ~6-8 min; --full adds ~20-40 min of mod-p eliminations)

Requires external/msolve (jcqft.gb_backend.MSOLVE); Singular optional
(irreducibility certificate skipped gracefully if absent).  All msolve /
Singular subprocesses run under the 16 GB gb_backend memory cap.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from jcqft.core import F, PHI, SRC, p  # noqa: E402
from jcqft.gb_backend import GB_MEM_MB, MSOLVE, _run  # noqa: E402

T0 = time.time()
N_CHECKS = 0
OUT_DIR = "/tmp/lattice_discriminant"
THREADS = str(min(8, os.cpu_count() or 1))
PRIME1, PRIME2 = 1073741789, 1073741783          # ~30-bit, msolve-safe
FACT_PRIMES = [268435399, 268435367, 268435361, 268435337]  # < 2^29 (Singular)

t = sp.Symbol("t")
V0 = sp.symbols("x0 y0 z0")
V1 = sp.symbols("x1 y1 z1")
VARS = list(V0) + list(V1)
EPS = sp.Rational(1, 4)
T1 = [sp.Rational(-1, 4), 0, 0, 0, 2, 0]
T3 = [1, 0, 0, 2, 1, 1]


def check(label, cond=True):
    global N_CHECKS
    assert cond, label
    N_CHECKS += 1
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


def note(label):
    print(f"  [..] {label}   ({time.time() - T0:.1f} s)")


# --------------------------------------------------------------------------
# msolve plumbing (all runs inherit the gb_backend 16 GB address-space cap)
# --------------------------------------------------------------------------

def _poly_str(f, gens):
    P = sp.Poly(sp.expand(f), *gens, domain="QQ")
    _, P = P.clear_denoms()
    s = sp.StrPrinter({"order": "none"}).doprint(P.as_expr())
    return s.replace("**", "^").replace(" ", "")


def write_ms(polys, gens, char, fn):
    with open(fn, "w") as fh:
        fh.write(",".join(map(str, gens)) + f"\n{char}\n")
        fh.write(",\n".join(_poly_str(f, gens) for f in polys) + "\n")


def msolve(args, timeout=3600):
    return _run([MSOLVE] + args, timeout)


_DYADIC = re.compile(r"(-?\d+) / 2\^(\d+)")


def _parse_out(fn):
    s = open(fn).read().strip().rstrip(":")
    s = _DYADIC.sub(lambda m: f"Fraction({m.group(1)}, 2**{m.group(2)})", s)
    return eval(s, {"Fraction": Fraction})       # noqa: S307 (msolve output)


def real_boxes(polys, gens, tag, extra=()):
    """Exact real-solution boxes of a 0-dim system over Q (msolve default
    mode: rational-univariate representation + certified isolation)."""
    fi, fo = f"{OUT_DIR}/{tag}.ms", f"{OUT_DIR}/{tag}.sol"
    write_ms(polys, gens, 0, fi)
    msolve(["-t", THREADS, *extra, "-f", fi, "-o", fo])
    data = _parse_out(fo)
    assert data[0] == 0, f"{tag}: system not zero-dimensional"
    return data[1][1]


def param_run(polys, gens, tag):
    """msolve -P 1 over Q: exact rational parametrization.  Returns
    (ideal degree, separating vector, eliminant coeffs ascending, real
    boxes)."""
    fi, fo = f"{OUT_DIR}/{tag}.ms", f"{OUT_DIR}/{tag}.sol"
    write_ms(polys, gens, 0, fi)
    msolve(["-t", THREADS, "-P", "1", "-f", fi, "-o", fo])
    data = _parse_out(fo)
    assert data[0] == 0, f"{tag}: system not zero-dimensional"
    char, nvars, deg, names, sepvec, par = data[1][:6]
    assert char == 0 and nvars == len(gens)
    w_deg, w_coeffs = par[1][0]
    assert w_deg == len(w_coeffs) - 1
    boxes = data[2][1] if len(data) > 2 else []
    return deg, sepvec, w_coeffs, boxes


def gb_modp(polys, gens, char, elim, tag, timeout=3600):
    """Reduced GB mod char with the first `elim` variables eliminated
    (block grevlex).  Returns the list of polynomial strings."""
    fi, fo = f"{OUT_DIR}/{tag}.ms", f"{OUT_DIR}/{tag}.gb"
    write_ms(polys, gens, char, fi)
    msolve(["-g", "2", "-e", str(elim), "-t", THREADS, "-f", fi, "-o", fo],
           timeout)
    body = "\n".join(ln for ln in open(fo).read().splitlines()
                     if not ln.startswith("#")).strip()
    assert body.startswith("[") and body.endswith("]:")
    return [s.strip() for s in body[1:-2].split(",")]


# --------------------------------------------------------------------------
# the deformed system
# --------------------------------------------------------------------------

def build_system():
    """Fiber equations F_eps(phi) - J(t) and det DF_eps at eps = 1/4."""
    F0 = [Fi.subs(dict(zip(PHI, V0))) for Fi in F]
    F1 = [Fi.subs(dict(zip(PHI, V1))) for Fi in F]
    Feps = ([F0[i] + EPS * (V1[i] - V0[i]) for i in range(3)]
            + [F1[i] + EPS * (V0[i] - V1[i]) for i in range(3)])
    J = [(1 - t) * T1[i] + t * T3[i] for i in range(6)]
    eqs = [sp.expand(Feps[i] - J[i]) for i in range(6)]
    A = sp.Matrix([[sp.diff(Fi, v) for v in V0] for Fi in F0])
    B = sp.Matrix([[sp.diff(Fi, v) for v in V1] for Fi in F1])
    I3 = sp.eye(3)
    # block identity: det [[A-e,e],[e,B-e]] = det((A-e)(B-e) - e^2), since
    # the off-diagonal blocks are scalar.
    det = sp.expand(((A - EPS * I3) * (B - EPS * I3) - EPS**2 * I3).det())
    DFe = sp.Matrix(
        [[sp.diff(Fi, v) for v in VARS] for Fi in Feps])
    sub = dict(zip(VARS, [1, 2, -1, 3, -2, sp.Rational(1, 2)]))
    assert (DFe.subs(sub).det() - det.subs(sub)) == 0
    return eqs, det, Feps


def n_real_at(eqs, tv, tag):
    """Exact number of real solutions of F_eps(phi) = J(tv)."""
    tv = sp.Rational(tv)
    polys = [f.subs(t, tv) for f in eqs]
    return len(real_boxes(polys, VARS, tag))


def ideal_degree_at(eqs, tv, tag):
    """Exact degree of the specialized (0-dim) fiber ideal, and the exact
    number of real solutions."""
    tv = sp.Rational(tv)
    polys = [f.subs(t, tv) for f in eqs]
    deg, _, _, boxes = param_run(polys, VARS, tag)
    return deg, len(boxes)


# --------------------------------------------------------------------------
# main sections
# --------------------------------------------------------------------------

def section1_eps0_anchors():
    print("=== 1. eps=0 anchors (exact, sympy only) ===")
    F0 = [Fi.subs(dict(zip(PHI, V0))) for Fi in F]
    F1 = [Fi.subs(dict(zip(PHI, V1))) for Fi in F]
    DF0 = sp.Matrix([[sp.diff(Fi, v) for v in VARS] for Fi in F0 + F1])
    check("eps=0: det DF_0 == 4 identically (Keller product; the fold "
          "ideal {det = 0} is EMPTY)", sp.expand(DF0.det()) == 4)

    J0 = [(1 - t) * T1[i] + t * T3[i] for i in range(3)]
    J1 = [(1 - t) * T1[i + 3] + t * T3[i + 3] for i in range(3)]
    p0 = sp.expand(p.subs(dict(zip(SRC, J0))))
    p1 = sp.expand(p.subs(dict(zip(SRC, J1))))
    check("eps=0 escape (site 0): p(J0(t)) = 4(5t-1), single root t = 1/5",
          sp.expand(p0 - 4 * (5 * t - 1)) == 0)
    q_seg = 107 * t**4 + 42 * t**3 - 85 * t**2 + 44 * t - 4
    check("eps=0 escape (site 1): p(J1(t)) = 107t^4+42t^3-85t^2+44t-4",
          sp.expand(p1 - q_seg) == 0)
    P1 = sp.Poly(p1, t)
    check("p(J1) has exactly 1 real root in (0,1) (Sturm), in (29/256,"
          " 30/256)",
          P1.count_roots(0, 1) == 1
          and P1.count_roots(sp.Rational(29, 256), sp.Rational(30, 256)) == 1)
    # product chamber sequence 9 -> 3 -> 1 across the two eps=0 walls
    seq = []
    for tv in (sp.Rational(1, 16), sp.Rational(3, 20), sp.Rational(1, 2)):
        n0 = 3 if p0.subs(t, tv) < 0 else 1
        n1 = 3 if p1.subs(t, tv) < 0 else 1
        seq.append(n0 * n1)
    check(f"eps=0 product chamber sequence on the segment: {seq} = [9, 3, 1]",
          seq == [9, 3, 1])
    return q_seg


def section2_fold(eqs, det, args):
    print("=== 2. fold polynomial f(t) at eps=1/4 (exact over Q, msolve) ===")
    det_0 = det.subs({v: 0 for v in VARS})
    det_pt = det.subs(dict(zip(VARS, [1, 2, -1, 3, -2, 1])))
    check(f"eps=1/4: det DF_eps NOT constant ({det_0} at 0 vs {det_pt}) -- "
          "fold walls possible", det_0 != det_pt)
    Pd = sp.Poly(det, *VARS)
    check(f"det DF_eps: total degree {Pd.total_degree()}, "
          f"{len(Pd.terms())} terms", Pd.total_degree() == 19)

    note("msolve -P 1 over Q on the 7-var fold system (~3-4 min) ...")
    deg, sepvec, w, boxes = param_run(eqs + [det], VARS + [t], "fold_q")
    check(f"fold ideal is 0-dimensional of degree {deg} = 516", deg == 516)
    check("t itself is a separating element (eliminant is exactly f(t))",
          list(sepvec) == [0, 0, 0, 0, 0, 0, 1])
    fpoly = sp.Poly(list(reversed(w)), t)      # msolve: ascending order
    check(f"f(t): degree {fpoly.degree()}, height "
          f"{len(str(max(abs(c) for c in w)))} digits, primitive "
          f"(content {sp.igcd(*[int(c) for c in w])})",
          fpoly.degree() == 516 and sp.igcd(*[int(c) for c in w]) == 1)
    with open(f"{OUT_DIR}/fold_poly_f_t.txt", "w") as fh:
        fh.write("# f(t): exact fold eliminant, ascending coefficients\n")
        fh.write("\n".join(str(c) for c in w) + "\n")
    note(f"full coefficient list written to {OUT_DIR}/fold_poly_f_t.txt")

    # cross-check against an independent mod-p elimination-order GB
    gens = gb_modp(eqs + [det], VARS + [t], PRIME1, 6, "fold_p1")
    g0 = sp.Poly(sp.sympify(gens[0].replace("^", "**")), t, modulus=PRIME1)
    fp = sp.Poly([int(c) % PRIME1 for c in reversed(w)], t, modulus=PRIME1)
    check(f"mod-{PRIME1} elimination GB reproduces f(t) (monic match)",
          g0.monic() == fp.monic())
    dfp = fp.diff(t)
    check("f squarefree over Q (gcd(f, f') = 1 mod p, lc(f) != 0 mod p)",
          int(w[-1]) % PRIME1 != 0 and sp.gcd(fp, dfp).degree() == 0)

    # irreducibility over Q by factor patterns mod several primes
    import shutil
    if shutil.which("Singular"):
        pats = []
        for pr in FACT_PRIMES:
            assert int(w[-1]) % pr != 0
            terms = [f"{int(c) % pr}*t^{i}" for i, c in enumerate(w)
                     if int(c) % pr]
            script = (f"ring R = {pr}, (t), dp;\npoly f = {'+'.join(terms)};\n"
                      "list L = factorize(f);\nfor (int i=1; i<=size(L[1]);"
                      " i++) { string(deg(L[1][i]), \":\", L[2][i]); }\n"
                      "quit;\n")
            fn = f"{OUT_DIR}/wf_{pr}.sing"
            open(fn, "w").write(script)
            out = _run(["Singular", "-q", fn], 600)
            degs = []
            for ln in out.splitlines():
                d_, m_ = (int(u) for u in ln.split(":"))
                if d_ > 0:
                    degs += [d_] * m_
            assert sum(degs) == 516, (pr, degs)
            pats.append(sorted(degs))
        common = None
        for degs in pats:
            S = {0}
            for d_ in degs:
                S |= {s + d_ for s in S}
            common = S if common is None else (common & S)
        check(f"f(t) IRREDUCIBLE over Q (factor patterns mod "
              f"{len(FACT_PRIMES)} primes, e.g. {pats[1]}; common subset "
              "sums = {0, 516})", common == {0, 516})
    else:
        note("Singular not found -- irreducibility certificate skipped")

    # real roots: 50 total, 14 in (0,1)
    troots = sorted([box[6] for box in boxes], key=lambda ab: ab[0])
    check(f"f(t) has exactly {len(troots)} = 50 real roots (msolve exact "
          "isolation); every real root has a real fold witness phi",
          len(troots) == 50)
    in01 = [ab for ab in troots if ab[1] > 0 and ab[0] < 1]
    check(f"exactly {len(in01)} = 14 real fold values in (0,1)",
          len(in01) == 14)
    esc_lo, esc_hi = Fraction(27, 512), Fraction(7, 128)
    check("NO fold root inside the old SS6.4 'first crossing' bracket "
          "(27/512, 7/128) -- that jump was not a fold",
          all(ab[1] <= esc_lo or ab[0] >= esc_hi for ab in in01))
    print("  real fold values in (0,1):")
    for lo, hi in in01:
        print(f"    t ~ {float((lo + hi) / 2):.9f}")
    return fpoly, in01


# 15 sample points strictly interleaving the 14 fold roots; EXPECT[i] is
# the exact real count on the open interval between fold roots i and i+1.
SAMPLES = [Fraction(1, 8192), Fraction(1, 1024), Fraction(1, 512),
           Fraction(15, 256), Fraction(43, 512), Fraction(23, 256),
           Fraction(55, 512), Fraction(15, 128), Fraction(63, 512),
           Fraction(1, 8), Fraction(9, 64), Fraction(1, 5),
           Fraction(2, 5), Fraction(17, 32), Fraction(3, 4)]
EXPECT = [18, 16, 14, 16, 14, 12, 10, 8, 6, 8, 6, 4, 2, 4, 6]


def section3_counts(eqs, in01):
    print("=== 3. exact real chamber counts along the segment (msolve) ===")
    roots = [(Fraction(lo), Fraction(hi)) for lo, hi in in01]
    for i, (lo, hi) in enumerate(roots):
        assert SAMPLES[i] < lo and hi < SAMPLES[i + 1], (i, lo, hi)
    check("the 15 sample points strictly interleave the 14 fold-root "
          "isolating intervals (fold i between samples i, i+1)")
    counts = []
    for tv in SAMPLES:
        n = n_real_at(eqs, tv, f"cnt_{tv.numerator}_{tv.denominator}")
        counts.append(n)
        print(f"    t = {str(tv):>7s} ~ {float(tv):.6f}:  N = {n}")
    check(f"exact count sequence across the 14 folds: {counts} == {EXPECT}",
          counts == EXPECT)
    jumps = [counts[i + 1] - counts[i] for i in range(len(counts) - 1)]
    check(f"every chamber jump is +-2 (pure fold; no odd escape jumps): "
          f"{jumps}", all(abs(j) == 2 for j in jumps))
    # interior escape value t = 1/13 lies between fold roots 3 and 4:
    # 15/256 < 1/13 < 5/64 < (fold root 4) -- same plateau on both sides
    n_q = n_real_at(eqs, Fraction(5, 64), "cnt_5_64")
    check(f"N(5/64) = {n_q} = 16 = N(15/256): plateau across the interior "
          "escape value t = 1/13 (no real chamber wall there)",
          n_q == 16 and counts[3] == 16)
    n_t1 = n_real_at(eqs, 1, "cnt_1_1")
    check(f"endpoint t=1 (T3): N = {n_t1} = 6 (matches the certified HC "
          "value of SS6.4)", n_t1 == 6)
    n_t0 = n_real_at(eqs, 0, "cnt_0_1")
    check(f"endpoint t=0 (T1): N = {n_t0} = 14 EXACT -- the SS6.4 HC count "
          "12 was an (admitted) incomplete lower bound", n_t0 == 14)
    check("t=0 is a REAL escape point: N(0) = 14 < N(0+) = 18 (four real "
          "solutions live at infinity over T1)",
          n_t0 == 14 and counts[0] == 18)
    return counts


# exact quartic factors of the escape polynomial, reconstructed by 2-prime
# CRT from the x0/y0 coordinate-eliminant leading coefficients (mod-p GBs;
# reproduced behind --full) and then certified over Q by the fiber-degree
# drops asserted below.
Q_Y = (29823777 * t**4 + 5199180 * t**3 + 713782 * t**2
       - 246740 * t + 12337)                       # no real roots
Q_X = (2841875 * t**4 + 125650 * t**3 - 1157957 * t**2
       + 512672 * t - 54016)         # real roots ~ -0.83010, 0.15629


def section4_escape(eqs, args):
    print("=== 4. escape polynomial e(t) at eps=1/4 ===")
    # generic control: full fiber degree 66
    dg, ng = ideal_degree_at(eqs, sp.Rational(1, 7), "deg_ctrl")
    check(f"generic control t=1/7: ideal degree {dg} = 66 (no escape), "
          f"N = {ng} = 6", dg == 66 and ng == 6)
    d0, n0 = ideal_degree_at(eqs, 0, "deg_t0")
    check(f"t=0: ideal degree {d0} = 54 < 66 -- EXACT escape certificate "
          f"(12 solutions at infinity over T1); N = {n0} = 14",
          d0 == 54 and n0 == 14)
    d13, n13 = ideal_degree_at(eqs, sp.Rational(1, 13), "deg_t13")
    check(f"t=1/13: ideal degree {d13} = 53 < 66 -- EXACT escape "
          f"certificate; N = {n13} = 13 (odd! 3 real solutions escape "
          "exactly AT t=1/13, count dips 16 -> 13 -> 16)",
          d13 == 53 and n13 == 13)
    dm3, nm3 = ideal_degree_at(eqs, sp.Rational(-1, 3), "deg_tm13")
    check(f"t=-1/3 (outside the segment): ideal degree {dm3} = 65 < 66 -- "
          "escape", dm3 == 65)

    for name, q, drop in (("q_y", Q_Y, 260), ("q_x", Q_X, 260)):
        qq = sp.Poly(q, t)
        check(f"{name} irreducible over Q with {len(sp.real_roots(qq))} "
              "real roots",
              len(sp.factor_list(q)[1]) == 1
              and sp.factor_list(q)[1][0][1] == 1)
        # adjoin the quartic: ideal degree < 4*66 = 264 certifies escape
        # over EVERY root (the drop is Galois-invariant for an
        # irreducible q).
        deg_adj, _, _, _ = param_run(
            [f for f in eqs] + [q], VARS + [t], f"deg_{name}")
        check(f"ideal degree with {name} adjoined: {deg_adj} = {drop} "
              "< 264 = 4 x 66 -- EXACT escape certificate for all roots "
              f"of {name} (one solution at infinity per root)",
              deg_adj == drop)
    # negative control: an irreducible quartic with no escape root must
    # give exactly 4 x 66.
    deg_ctrl, _, _, _ = param_run(
        [f for f in eqs] + [t**4 + t + 1], VARS + [t], "deg_ctrlq")
    check(f"control quartic t^4+t+1 adjoined: ideal degree {deg_ctrl} "
          "= 264 exactly (no spurious drops)", deg_ctrl == 264)

    e_t = sp.expand(t * (13 * t - 1) * (3 * t + 1) * Q_Y * Q_X)
    print("  e(t) = t (13t-1) (3t+1) q_y(t) q_x(t),  degree "
          f"{sp.Poly(e_t, t).degree()}")
    ee = sp.Poly(e_t, t)
    r_lo, r_hi = sp.Rational(27, 512), sp.Rational(7, 128)
    check("e(t) has NO root in the old SS6.4 'first crossing' bracket "
          "(27/512, 7/128) -- with section 3 this refutes the SS6.4 "
          "12->13 escape-jump reading (HC completeness artifact)",
          ee.count_roots(r_lo, r_hi) == 0)
    d_br, n_br = ideal_degree_at(eqs, sp.Rational(109, 2048), "deg_br")
    check(f"spot check inside that bracket (t=109/2048): degree {d_br} "
          f"= 66, N = {n_br} = 14 -- generic, nothing happens there",
          d_br == 66 and n_br == 14)
    in01 = ee.count_roots(sp.Rational(1, 10**9), 1)
    check(f"e(t) real roots in (0,1): {in01} (t = 1/13 and the q_x root "
          "~ 0.156290) plus the endpoint t = 0", in01 == 2)
    note("escape values are chamber-invisible in R except AT the value "
         "(measure-zero dips); the real walls of section 3 are all folds")

    if args.full:
        note("--full: recomputing the x0/y0 coordinate-eliminant lc's mod "
             "2 primes (completeness certificate for e(t); ~25 min) ...")
        for pr in (PRIME1, PRIME2):
            for keep, shape in ((V0[0], t**16 * Q_X),
                                (V0[1], t**12 * (3 * t + 1) * (13 * t - 1)
                                 * Q_Y)):
                others = [v for v in VARS if v != keep]
                gens = gb_modp(eqs, others + [keep, t], pr, 5,
                               f"esc_{keep}_{pr}", timeout=3600)
                plane = [g for g in gens
                         if not (set(re.findall(r"[a-z]\d", g))
                                 & {str(v) for v in others})]
                assert len(plane) == 1, f"lc_{keep} mod {pr}"
                E = sp.sympify(plane[0].replace("^", "**"))
                Ev = sp.Poly(E, keep)
                lc = sp.Poly(Ev.LC(), t, modulus=pr).monic()
                check(f"mod {pr}: eliminant deg_{keep} = 66 (separating) "
                      f"and lc_{keep}(t) has exactly the e(t) factors "
                      f"({sp.factor(shape)})",
                      Ev.degree() == 66
                      and lc == sp.Poly(shape, t, modulus=pr).monic())
    else:
        note("completeness of the e(t) factor list rests on the recorded "
             "x0/y0/x1/y1 mod-p eliminations (2 primes; rerun with "
             "--full); the z-direction eliminations exceeded the compute "
             "budget (documented in docs SS6.5)")
    return e_t


def section5_surface(eqs_abc, det, args):
    print("=== 5. (--full) fold surface W(a,b,c), J1 frozen at (2,1,1) ===")
    if not args.full:
        note("skipped (enable with --full); see docs SS6.5 for the "
             "documented state")
        return
    a, b, c = sp.symbols("a b c")
    note(f"mod-{PRIME1} elimination of the 6 field variables from the "
         "7-equation fold system in (phi, a, b, c) -- prescreen only; "
         f"memory-capped at {GB_MEM_MB} MB, 60 min timeout ...")
    try:
        gens = gb_modp(eqs_abc + [det], VARS + [a, b, c], PRIME1, 6,
                       "surface_p1", timeout=3600)
        pure = [g for g in gens
                if not (set(re.findall(r"[a-z]\d", g)))]
        note(f"GB has {len(gens)} elements, {len(pure)} in (a,b,c) only")
        if pure:
            W = sp.Poly(sp.sympify(pure[0].replace("^", "**")), a, b, c)
            check(f"fold surface W(a,b,c) mod p: degree {W.total_degree()},"
                  f" {len(W.terms())} terms (over-Q lift not attempted)",
                  W.total_degree() > 0)
    except (RuntimeError, TimeoutError) as exc:
        note(f"surface elimination hit the documented wall: {exc}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="mod-p completeness rerun + fold-surface prescreen")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    assert os.path.isfile(MSOLVE) and os.access(MSOLVE, os.X_OK), (
        "external/msolve/msolve required (see README); aborting")

    section1_eps0_anchors()
    note("building det DF_eps (3x3 block identity + 6x6 spot check) ...")
    eqs, det, Feps = build_system()
    fpoly, in01 = section2_fold(eqs, det, args)
    section3_counts(eqs, in01)
    section4_escape(eqs, args)
    a, b, c = sp.symbols("a b c")
    eqs_abc = ([sp.expand(Feps[i] - v) for i, v in enumerate((a, b, c))]
               + [sp.expand(Feps[i + 3] - v)
                  for i, v in enumerate((2, 1, 1))])
    section5_surface(eqs_abc, det, args)

    print(f"\nall {N_CHECKS} checks passed in {time.time() - T0:.1f} s")
    print("VERDICT (eps = 1/4, segment T1 -> T3): the real chamber walls "
          "are 14 FOLD roots of the")
    print("irreducible degree-516 f(t) (counts 18 16 14 16 14 12 10 8 6 8 "
          "6 4 2 4 6, all jumps +-2);")
    print("ESCAPE happens at t=0 (T1: 4 real solutions at infinity) and at "
          "the roots of e(t) =")
    print("t (13t-1) (3t+1) q_y q_x, where the count only dips AT the "
          "value; the SS6.4 odd jumps were")
    print("homotopy-continuation completeness artifacts, corrected here "
          "exactly.")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
