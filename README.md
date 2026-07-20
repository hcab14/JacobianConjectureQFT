# JacobianConjectureQFT

Exploration of the 2026 Alpöge–Mathew counterexample to the Jacobian
conjecture through the lens of zero-dimensional perturbative quantum field
theory (tree-graph expansions, branch loci, monodromy, non-perturbative
"vacua at infinity"), with an eye toward implications for mathematically
rigorous QFT.

- Problem statement: [`docs/PROBLEM.md`](docs/PROBLEM.md)
- Running progress log with all results: [`docs/PROGRESS.md`](docs/PROGRESS.md)

## Setup

```bash
uv venv .venv && uv pip install --python .venv/bin/python -r requirements.txt
# or: python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

## Scripts

| Script | Purpose | Runtime |
|---|---|---|
| `verify_counterexample.py` | Independent verification: det DF = -2, non-injectivity, generic 3-point fiber | seconds |
| `counterexample.py` | Shared definitions (the map, propagator, cubic eliminant) | library |
| `tree_expansion.py` | Tree-Feynman-graph inversion to order 10; ray series to t^60; radius of convergence vs exact branch points; Newton branch tracking | ~1 min |
| `branch_locus.py` | Exact geometry: eliminants, discriminant classification, escape-to-infinity locus, exact fibers | seconds |
| `monodromy.py` | Numerical monodromy of the 3 sheets (see `docs/MONODROMY.md`) | minutes |
| `measure_anomaly.py` | Field-redefinition measure anomaly: exact chamber rule for the real preimage count, Monte Carlo anomaly factor, rational sum-over-sheets observables | ~6 s |
| `search_counterexamples.py` | Construction mechanism + search for new counterexamples (see `docs/NEW_COUNTEREXAMPLES.md`) | minutes |

Run any script with `.venv/bin/python <script>.py`.
