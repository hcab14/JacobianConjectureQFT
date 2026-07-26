"""Q2a probe (docs/CLASSICAL_MAP_INVARIANTS.md §5.2): does the real chamber
function N_eps(J) of the 2-site lattice Alpöge–Mathew deformation

    F_eps(phi)_0 = F(phi_0) + eps (phi_1 - phi_0)
    F_eps(phi)_1 = F(phi_1) + eps (phi_0 - phi_1)

stay NON-CONSTANT for kinetic coupling eps > 0 — i.e. does wall-crossing
survive kinetic mixing?  (Convention matches
scripts/classical_map_invariants_probe.py: F_M^K = F_M + K.phi with
K = eps [[-I, I], [I, -I]].)

This wrapper:
  1. asserts the exact eps=0 product structure against jcqft.fibers
     (complex fiber 3x3 = 9; real fiber = N(J_a) N(J_b) in {9, 3, 1});
  2. exports the 9 exact product solutions per target as high-precision
     start points, plus the eps ladder, to /tmp/lattice_chamber/input.jl;
  3. launches julia scripts/hc_lattice_chamber.jl (HomotopyContinuation.jl:
     master monodromy set over joint (eps, J) parameters + eps-parameter
     homotopy from eps=0 + fresh polyhedral solves, reality CERTIFIED
     against the exact rational systems);
  4. parses and re-asserts the cross-checks, prints the results table and
     the Q2a verdict.

Epistemic status: real/complex counts are CERTIFIED (interval arithmetic)
lower bounds at exact rational (eps, J); completeness of the solution
lists is numerical (all D master paths tracked + monodromy stabilization
+ fresh-solve cross-check), not interval-certified.

Usage:  .venv/bin/python scripts/lattice_chamber.py [--full]
        (default ~6 min; --full adds eps in {1/1000, 1/50, 4})

Julia is resolved from PATH or ~/.juliaup/bin.  One-time setup if
HomotopyContinuation.jl is missing:
        julia -e 'using Pkg; Pkg.add("HomotopyContinuation")'
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from jcqft.core import D0, F, PHI, SRC, p  # noqa: E402
from jcqft.fibers import exact_fiber  # noqa: E402

T0 = time.time()
N_CHECKS = 0
OUT_DIR = "/tmp/lattice_chamber"


def check(label, cond=True):
    global N_CHECKS
    assert cond, label
    N_CHECKS += 1
    print(f"  [ok] {label}   ({time.time() - T0:.1f} s)")


# Targets: pairs of per-site rational chamber points (scripts/witten_index.py).
# name, J_site0, J_site1, expected eps=0 real product count
TARGETS = [
    ("T1", (sp.Rational(-1, 4), 0, 0), (0, 2, 0), 9),   # N=3 x N=3
    ("T2", (sp.Rational(-1, 4), 0, 0), (1, 0, 0), 3),   # N=3 x N=1
    ("T3", (1, 0, 0), (2, 1, 1), 1),                    # N=1 x N=1
]
EPS_DEFAULT = ["1//100", "1//10", "1//4", "1//2", "1//1", "2//1"]
EPS_FULL_EXTRA = ["1//1000", "1//50", "4//1"]


def _fraction(v):
    q = sp.Rational(v)
    return Fraction(int(q.p), int(q.q))


def _jl_rat(v):
    f = _fraction(v)
    return f"{f.numerator}//{f.denominator}"


def _jl_c128(v):
    cv = complex(sp.N(v, 30))
    return f"ComplexF64({cv.real!r}, {cv.imag!r})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="extra eps values 1/1000, 1/50, 4")
    args = ap.parse_args()
    eps_list = sorted(EPS_DEFAULT + (EPS_FULL_EXTRA if args.full else []),
                      key=lambda s: Fraction(s.replace("//", "/")))

    # =======================================================================
    print("=== 1. exact eps=0 product structure (jcqft.fibers) ===")
    site_fibers = {}
    for name, Ja, Jb, n0 in TARGETS:
        n_prod, fibs = 1, []
        for J in (Ja, Jb):
            sub = dict(zip(SRC, J))
            assert p.subs(sub) != 0 and D0.subs(sub) != 0
            fib = exact_fiber(J)
            assert len(fib) == 3
            n_r = sum(1 for pt in fib if all(v.is_real for v in pt))
            assert n_r == (3 if p.subs(sub) < 0 else 1)
            n_prod *= n_r
            fibs.append(fib)
        check(f"{name}: eps=0 complex fiber 3x3 = 9, real product "
              f"{n_prod} = expected {n0}", n_prod == n0)
        site_fibers[name] = fibs

    # honesty check: kinetic mixing leaves the Keller class -- det DF_eps
    # is NOT constant for eps > 0 (so real counts may also change through
    # finite fold bifurcations {det DF_eps = 0}, not only via infinity).
    eps_s = sp.Rational(1, 4)
    v0, v1 = sp.symbols("X0 Y0 Z0"), sp.symbols("X1 Y1 Z1")
    allv = list(v0) + list(v1)
    F0 = [Fi.subs(dict(zip(PHI, v0))) for Fi in F]
    F1 = [Fi.subs(dict(zip(PHI, v1))) for Fi in F]
    Feps = ([F0[i] + eps_s * (v1[i] - v0[i]) for i in range(3)]
            + [F1[i] + eps_s * (v0[i] - v1[i]) for i in range(3)])
    DFe = sp.Matrix([[sp.diff(Fi, v) for v in allv] for Fi in Feps])
    det_0 = DFe.subs({v: 0 for v in allv}).det()
    det_pt = DFe.subs(dict(zip(allv, [1, 2, -1, 3, -2, 1]))).det()
    check(f"eps=1/4: det DF_eps NOT constant ({det_0} at 0 vs {det_pt} "
          "at a point) -- F_eps leaves the Keller class",
          sp.simplify(det_0 - det_pt) != 0)

    # =======================================================================
    print("=== 2. export start data for Julia ===")
    os.makedirs(OUT_DIR, exist_ok=True)
    input_jl = os.path.join(OUT_DIR, "input.jl")
    lines = ["# generated by scripts/lattice_chamber.py -- do not edit",
             "eps_list = Rational{Int}[" + ", ".join(eps_list) + "]",
             "targets = ["]
    for name, Ja, Jb, n0 in TARGETS:
        Jflat = ", ".join(_jl_rat(v) for v in tuple(Ja) + tuple(Jb))
        lines.append(f"  (name = \"{name}\", "
                     f"J = Rational{{Int}}[{Jflat}], N0 = {n0}, starts = [")
        fib_a, fib_b = site_fibers[name]
        for pt_a in fib_a:
            for pt_b in fib_b:
                coords = ", ".join(_jl_c128(v) for v in pt_a + pt_b)
                lines.append(f"    ComplexF64[{coords}],")
        lines.append("  ]),")
    lines.append("]")
    lines.append("bisect = (eps = 1//4, from = \"T1\", to = \"T3\", "
                 "tol = 1//512)")
    with open(input_jl, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    check(f"wrote {input_jl} (3 targets x 9 exact product starts, "
          f"{len(eps_list)} eps values)")

    # =======================================================================
    print("=== 3. run HomotopyContinuation.jl driver ===")
    julia = shutil.which("julia") or os.path.expanduser("~/.juliaup/bin/julia")
    assert os.path.exists(julia), (
        "julia not found; install via juliaup, then "
        "julia -e 'using Pkg; Pkg.add(\"HomotopyContinuation\")'")
    driver = os.path.join(os.path.dirname(__file__), "hc_lattice_chamber.jl")
    proc = subprocess.Popen([julia, driver, input_jl],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True)
    out_lines = []
    for line in proc.stdout:
        line = line.rstrip("\n")
        out_lines.append(line)
        print(f"    | {line}")
    rc = proc.wait()
    check(f"Julia driver exited 0 (all internal cross-checks passed)",
          rc == 0)

    # =======================================================================
    print("=== 4. parse, re-assert, verdict ===")
    results = {}   # (target, eps) -> dict
    master_D = None
    bisects = []
    for line in out_lines:
        f = line.split("\t")
        if f[0] == "MASTER":
            master_D = int(f[1].split("=")[1])
        elif f[0] == "RESULT":
            results[(f[1], f[2])] = dict(
                D=int(f[3]), n_complex=int(f[4]), n_real=int(f[5]),
                n_eps0=int(f[6]), n_eps0_real=int(f[7]), n_fresh=int(f[8]))
        elif f[0] == "BISECT":
            bisects.append(f[1:])
    check(f"master degree over generic complex (eps, J): D = {master_D}",
          master_D is not None and master_D >= 9)
    check(f"all {len(TARGETS)} x {len(eps_list)} (target, eps) results "
          "present",
          all((n, e) in results for n, _, _, _ in TARGETS for e in eps_list))
    for (name, eps), r in sorted(results.items()):
        assert 0 <= r["n_real"] <= r["n_complex"] <= master_D, (name, eps, r)
        assert (r["n_complex"] - r["n_real"]) % 2 == 0, (name, eps, r)
        assert r["n_eps0_real"] <= r["n_real"], (name, eps, r)
    check("per (target, eps): 0 <= n_real <= n_complex <= D, conjugation "
          "parity, eps0-descendant reals <= total reals")

    # the Q2a verdict: is N_eps non-constant across targets at each eps?
    print()
    print("  results table: certified real count  (complex count | "
          "real descending from eps=0)")
    hdr = "    eps       " + "".join(f"{n:>16s}" for n, _, _, _ in TARGETS)
    print(hdr)
    print("    0 (exact) " + "".join(f"{n0:>10d} (9|{n0})"
                                     for _, _, _, n0 in TARGETS))
    nonconst = {}
    for eps in eps_list:
        row, counts = f"    {eps:<10s}", []
        for name, _, _, _ in TARGETS:
            r = results[(name, eps)]
            counts.append(r["n_real"])
            row += f"{r['n_real']:>10d} ({r['n_complex']}|{r['n_eps0_real']})"
        nonconst[eps] = len(set(counts)) > 1
        print(row + ("   non-constant" if nonconst[eps] else "   CONSTANT"))
    check("Q2a: N_eps(J) is NON-CONSTANT at every probed eps > 0 "
          "(wall-crossing survives kinetic mixing)",
          all(nonconst.values()))
    for b in bisects:
        print(f"  bisection at eps={b[0]}: real count jumps {b[5]} -> {b[6]} "
              f"on {b[1]}->{b[2]} segment, t in ({b[3]}, {b[4]})")
    check("bisection bracketed at least one wall crossing in J at fixed "
          "eps > 0", len(bisects) >= 1)

    print(f"\nall {N_CHECKS} checks passed in {time.time() - T0:.1f} s")
    print("Q2a VERDICT (probed eps and J; certified reality, numerical "
          "completeness): YES —")
    print("the real chamber function N_eps(J) remains non-constant for "
          "every probed eps > 0.")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
