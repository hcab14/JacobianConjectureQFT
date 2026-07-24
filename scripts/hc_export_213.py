"""Export the six reduced degree-2 (2,-1,-3) mechanism queries for
numerical solving with HomotopyContinuation.jl (docs/SEARCH_213.md §5.3).

The exact Groebner campaign left all six queries unresolved (F4 memory
wall / std time wall).  This flips the question from "prove empty" to
"find a point": homotopy continuation is memory-light; if it finds a
(certifiable) solution we have a counterexample candidate, if it
robustly finds none that is strong (non-rigorous) evidence of emptiness
and tells us where to spend exact effort.

Writes one Julia file per query to /tmp/hc_213/:  query_<i>.jl defines
`polys` (Vector of expressions) and `vars`, then the driver
scripts/hc_drive_213.jl measures the polyhedral path count first and
only tracks within budget.

Usage:  .venv/bin/python scripts/hc_export_213.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp

from deg2_213_elim import build_box, eliminate_a_block, mech_queries_reduced

OUT = "/tmp/hc_213"


def jl(e):
    """sympy -> Julia source (integer/rational coefficients kept exact)."""
    s = sp.StrPrinter({"order": "none"}).doprint(sp.expand(e))
    return s.replace("**", "^")


def main():
    os.makedirs(OUT, exist_ok=True)
    box = build_box()
    a_syms = box["a_syms"]
    bcde = [s for s in box["unknowns"] if s not in a_syms]
    sol_a, reduced, _ = eliminate_a_block(box["keller"], a_syms)
    queries = mech_queries_reduced(box, sol_a, reduced, bcde)
    index = []
    for i, (label, (sys_, gens)) in enumerate(sorted(queries.items())):
        path = os.path.join(OUT, f"query_{i}.jl")
        with open(path, "w") as f:
            f.write(f"# {label}: {len(sys_)} equations, {len(gens)} "
                    "unknowns (see scripts/hc_export_213.py)\n")
            f.write("@var " + " ".join(str(g) for g in gens) + "\n")
            f.write("vars = [" + ", ".join(str(g) for g in gens) + "]\n")
            f.write("polys = [\n")
            for e in sys_:
                f.write("    " + jl(e) + ",\n")
            f.write("]\n")
            f.write(f'label = "{label}"\n')
        index.append(f"{i}\t{label}\t{len(sys_)}\t{len(gens)}")
        print(f"  wrote {path}  ({label}: {len(sys_)} eqs, "
              f"{len(gens)} unknowns)")
    with open(os.path.join(OUT, "index.tsv"), "w") as f:
        f.write("\n".join(index) + "\n")
    print(f"done -> {OUT}")


if __name__ == "__main__":
    main()
