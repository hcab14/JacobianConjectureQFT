# JacobianConjectureQFT

Exploration of the 2026 Alpöge–Mathew counterexample to the Jacobian
conjecture through the lens of zero-dimensional perturbative quantum field
theory (tree-graph expansions, branch loci, monodromy, non-perturbative
"vacua at infinity"), with an eye toward implications for mathematically
rigorous QFT.

- Problem statement: [`docs/PROBLEM.md`](docs/PROBLEM.md)
- Running progress log with all results: [`docs/PROGRESS.md`](docs/PROGRESS.md)
- Implications for rigorous QFT: [`docs/QFT_IMPLICATIONS.md`](docs/QFT_IMPLICATIONS.md)
- Connection to the amplitudes program: [`docs/AMPLITUDES_CONNECTION.md`](docs/AMPLITUDES_CONNECTION.md)
- Monodromy of the three sheets: [`docs/MONODROMY.md`](docs/MONODROMY.md)
- Search strategies for new counterexamples: [`docs/SEARCH_STRATEGIES.md`](docs/SEARCH_STRATEGIES.md)
- Chamber geometry and the positive-geometry verdict: [`docs/POSITIVE_GEOMETRY.md`](docs/POSITIVE_GEOMETRY.md)
- The wall complement as the braid-group classifying space: [`docs/WALL_COMPLEMENT.md`](docs/WALL_COMPLEMENT.md)
- All open questions, consolidated and ranked by tractability: [`docs/OPEN_QUESTIONS.md`](docs/OPEN_QUESTIONS.md)

## Setup

```bash
uv venv .venv && uv pip install --python .venv/bin/python -e .
# or: python3 -m venv .venv && .venv/bin/pip install -e .
```

This installs the `jcqft` package (editable) together with its dependencies
(sympy, mpmath, numpy).

## Layout

Shared library code lives in the `jcqft/` package; every analysis is a
runnable script in `scripts/`.

| Module | Contents |
|---|---|
| `jcqft/core.py` | The map F, propagator L, x-eliminant cubic (p, q, r), discriminant factor D0, C\*-weights, real chamber rule |
| `jcqft/truncated.py` | Sparse truncated-ring polynomial arithmetic; the tree-graph (Picard) formal inverse |
| `jcqft/fibers.py` | Lex Groebner fiber parametrization, exact fiber solver, lambdified numeric helpers |
| `jcqft/reduction.py` | C\*-equivariant normal form (P, Q, R) and the 2D-reduced Keller condition |

## Scripts

| Script | Purpose | Runtime |
|---|---|---|
| `scripts/verify_counterexample.py` | Independent verification: det DF = -2, non-injectivity, generic 3-point fiber | seconds |
| `scripts/tree_expansion.py` | Tree-Feynman-graph inversion to order 10; ray series to t^60; radius of convergence vs exact branch points; Newton branch tracking | ~1 min |
| `scripts/branch_locus.py` | Exact geometry: eliminants, discriminant classification, escape-to-infinity locus, exact fibers | seconds |
| `scripts/monodromy.py` | Numerical monodromy of the 3 sheets (see `docs/MONODROMY.md`) | minutes |
| `scripts/measure_anomaly.py` | Field-redefinition measure anomaly: exact chamber rule for the real preimage count, Monte Carlo anomaly factor, rational sum-over-sheets observables | ~6 s |
| `scripts/trace_pushforward.py` | Amplitudes-program structures: trace/pushforward rationality, boundary factorization at the wall (see `docs/AMPLITUDES_CONNECTION.md`) | ~1 s |
| `scripts/search_counterexamples.py` | Construction mechanism (2D reduction of the Keller condition) + first-order rigidity and continuation analysis (see `docs/NEW_COUNTEREXAMPLES.md`) | ~40 s |
| `scripts/positive_geometry.py` | C\*-reduction of the wall to a cuspidal plane cubic; exact non-surjectivity locus; canonical-form test: the chamber is NOT a positive geometry (see `docs/POSITIVE_GEOMETRY.md`) | ~2 s |
| `scripts/plot_chamber.py` | Figure of the reduced chamber geometry (`docs/img/chamber_geometry.png`) | ~2 s |
| `scripts/cusp_traces.py` | Exact trace asymptotics at the wall and cusp: floor(k/2) law, cube-root-of-unity cancellation, exact solvability on the cuspidal tangent (see `docs/POSITIVE_GEOMETRY.md` §6) | ~2 s |
| `scripts/pushforward_forms.py` | Pushforward of general forms: poles only on the wall, D0-cancellation between sheets, pole-order law (closes `docs/AMPLITUDES_CONNECTION.md` §2.4 Q2) | ~8 s |
| `scripts/wall_braid.py` | The wall is affinely the A2 discriminant: invariant eliminant, explicit affine isomorphism, cusp-loop = Coxeter element, local-system decomposition (see `docs/WALL_COMPLEMENT.md`) | ~6 s |

Run any script from the repository root with
`.venv/bin/python scripts/<name>.py`.
