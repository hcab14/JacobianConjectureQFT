# JacobianConjectureQFT

Exploration of the 2026 Alpöge–Mathew counterexample to the Jacobian
conjecture through the lens of zero-dimensional perturbative quantum field
theory (tree-graph expansions, branch loci, monodromy, non-perturbative
"vacua at infinity"), with an eye toward implications for mathematically
rigorous QFT — especially **AQFT** and **pAQFT**.

## Paper

[![DOI](https://zenodo.org/badge/1306864803.svg)](https://doi.org/10.5281/zenodo.21569061)

**Preprint (PDF):** [`paper/main.pdf`](paper/main.pdf)
([LaTeX source](paper/main.tex))

*A zero-dimensional calibration of field redefinitions for AQFT and pAQFT:
the Alpöge–Mathew counterexample and uniqueness of its defect*
(Christoph Mayer, July 2026).

**Cite this repository (Zenodo):**
[https://doi.org/10.5281/zenodo.21569062](https://doi.org/10.5281/zenodo.21569062)
(concept DOI → always the latest release). This version:
[https://doi.org/10.5281/zenodo.21569061](https://doi.org/10.5281/zenodo.21569061).

Headline results: exact structure of the observable algebra \(\mathrm{im}\,F^*\);
calibration of the formal deferral in pAQFT; uniqueness of the
non-properness defect in the equivariant family \((1,-1,-m)\) for
\(m\ge 3\). Every exact claim is asserted by a reproducible script below.
AI-assisted exploratory research — see provenance.

Reproduce the uniqueness theorem in ~1 minute:

```bash
.venv/bin/python scripts/search_11m.py
```

## Status and provenance

This is **AI-assisted exploratory research**, produced in interactive
sessions with AI coding agents and **not peer-reviewed**. The counterexample
itself is due to **Alpöge and Mathew** (announced July 19, 2026); this
repository only studies it. Read the claims with the following calibration,
which the documents themselves state case by case:

- **Exact/symbolic results** (most of the repository): every displayed
  identity is asserted by a reproducible sympy script listed below; "verified"
  means the script passes, within the stated scope (often box- or
  chart-limited — the scope lines matter).
- **Numerical results** (monodromy, rigidity continuation, Monte Carlo):
  labelled as evidence, not proof.
- **Interpretation** (QFT/amplitudes readings): clearly flagged opinion;
  see the claims-vs-limitations ledger in
  [`docs/QFT_IMPLICATIONS.md`](docs/QFT_IMPLICATIONS.md).

Corrections and counter-arguments are very welcome — please open an issue.

- Problem statement: [`docs/PROBLEM.md`](docs/PROBLEM.md)
- Running progress log with all results: [`docs/PROGRESS.md`](docs/PROGRESS.md)
- Implications for rigorous QFT: [`docs/QFT_IMPLICATIONS.md`](docs/QFT_IMPLICATIONS.md)
- Connection to the amplitudes program: [`docs/AMPLITUDES_CONNECTION.md`](docs/AMPLITUDES_CONNECTION.md)
- Monodromy of the three sheets: [`docs/MONODROMY.md`](docs/MONODROMY.md)
- Search strategies for new counterexamples: [`docs/SEARCH_STRATEGIES.md`](docs/SEARCH_STRATEGIES.md)
- Chamber geometry and the positive-geometry verdict: [`docs/POSITIVE_GEOMETRY.md`](docs/POSITIVE_GEOMETRY.md)
- The wall complement as the braid-group classifying space: [`docs/WALL_COMPLEMENT.md`](docs/WALL_COMPLEMENT.md)
- Twisted cohomology of the wall complement and the fate of the period pairing: [`docs/TWISTED_PERIODS.md`](docs/TWISTED_PERIODS.md)
- The infinity prefilter and rigidity in larger degree boxes: [`docs/RIGIDITY_AND_PREFILTER.md`](docs/RIGIDITY_AND_PREFILTER.md)
- The damped partition function: finiteness, chamber prefactor, uniformity exponents: [`docs/DAMPED_PARTITION.md`](docs/DAMPED_PARTITION.md)
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
| `jcqft/prefilter.py` | Infinity prefilter: millisecond properness-at-infinity rejection test for search candidates (plain and weighted) |
| `jcqft/reduction_w.py` | The (1,-1,-m) family: invariants, equivariant normal form, exact reduced Keller identity J₂(PRᵐ,QR)=κRᵐ, polynomiality boxes |
| `jcqft/reduction_213.py` | The (2,-1,-3) system: A₁-cone invariant theory, weight modules, chart-based reduced Keller identity, orbifold-mechanism witnesses |

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
| `scripts/twisted_cohomology.py` | Exact twisted cohomology of the wall complement via Fox calculus + Wang-sequence cross-check; jump loci = trefoil Alexander roots; integer-twist exactness of the canonical form (see `docs/TWISTED_PERIODS.md`) | ~3 s |
| `scripts/witten_prefilter.py` | Validates the infinity prefilter: Alpöge–Mathew survives with witness [1:0:0], proper maps rejected in ~1 ms, 200-map benchmark (see `docs/RIGIDITY_AND_PREFILTER.md`) | ~6 s |
| `scripts/rigidity_boxes.py` | Rigidity of the counterexample in strictly larger degree boxes: all non-gauge first-order deformations obstructed at 2nd order, continuation confirms (see `docs/RIGIDITY_AND_PREFILTER.md`) | ~50 s |
| `scripts/damped_partition.py` | Damped partition function: exact closed form and finiteness for all J, prefactor = N(J)/2 per chamber, wall/cusp corrections, uniformity exponents (see `docs/DAMPED_PARTITION.md`) | ~30 s |
| `scripts/missing_observables.py` | Exact structure of im F\*: basis {1,x,x²}, normal form and membership criterion, non-integrality escape certificate, fiber exhaustiveness (see `docs/MISSING_OBSERVABLES.md`) | ~15 s |
| `scripts/reduction_113.py` | Exact reduced Keller identity for all weight systems (1,-1,-m): generic-function proof, compact form J₂(PRᵐ,QR)=κRᵐ, polynomiality boxes (module `jcqft/reduction_w.py`) | ~4 s |
| `scripts/search_213.py` | Weight system (2,-1,-3): exact reduction on the A₁-cone quotient, complete degree-1 Keller classification (all tame automorphisms), Z₂ and Z₃ orbifold mechanisms provably empty (module `jcqft/reduction_213.py`, see `docs/SEARCH_213.md`) | ~40 s (`--full` ~8 min) |
| `scripts/search_113.py` | Weight system (1,-1,-3): complete Keller classification of the v-linear class for ALL w-degrees (Wronskian stratification, all tame automorphisms), 3:1 orbifold mechanism empty, Alpöge–Mathew stratum numerologically obstructed at m=3 (see `docs/SEARCH_113.md`) | ~20 min (`--full` ~1 h) |
| `scripts/search_11m.py` | THE UNIQUENESS THEOREM: for every m ≥ 3, the v-linear class of (1,-1,-m) contains only tame automorphisms — Alpöge–Mathew is the unique member of its equivariant family. 78 assertions, m symbolic where possible, spot checks m = 2,3,4,5,7 (see `docs/SEARCH_11M.md`) | ~1 min (`--full` ~72 min) |
| `scripts/witten_index.py` | The 0D Witten index: Brouwer degree = −N(J), jumps −1 ↔ −3 across the wall (exact non-properness certificate); Mathai–Quillen SUSY completion (no superpotential exists: DF ≠ DFᵀ), closed form and finiteness of Z_σ, localization numerics (see `docs/WITTEN_INDEX.md`) | ~9 s |
| `scripts/symmetric_search.py` | The symmetric (variational) problem: two explicit dimension-6 gradient Jacobian counterexamples from AM (cotangent lift det Hess ≡ −4, dBvdE twisted lift); gradient no-go in dimension 3 (K·DF symmetric ⟹ K=0; (1,-1,-m) gradient family only tame shears); direct boxes n=2 ≤ deg 6, n=3 ≤ deg 4 empty; coercivity obstruction (see `docs/SYMMETRIC_SEARCH.md`) | ~3.5 min (`--full` hours) |
| `scripts/classical_map_invariants_probe.py` | Classical-map invariants as Lagrangian data: ultralocal $F^{\times N}$ tensors Keller/fibers/wall/Witten index; Galerkin + kinetic mixing — non-properness survives on equal-mode slice; packing API on AM vs tame shear (see `docs/CLASSICAL_MAP_INVARIANTS.md` §6) | ~2 s |
| `scripts/lattice_chamber.py` | Q2a answered for the 2-site lattice: the real chamber function N_ε(J) stays non-constant for all probed ε > 0 (certified reality via HomotopyContinuation.jl, driver `scripts/hc_lattice_chamber.jl`); up to 57 of the 66 master solutions return from infinity at ε > 0, escape-type ±1 wall jumps bisected in J (see `docs/CLASSICAL_MAP_INVARIANTS.md` §6.4) | ~3.5 min (`--full` ~4 min; needs julia) |
| `scripts/bf_caricature.py` | The 0D Buchholz–Fredenhagen S(J) caricature: causal factorization provably vacuous on a point; the dynamical relation forces the fiber-algebra bundle — wall = rank jump 3→1→0 + pole divisor, 1↔3 = character count, S₃ = transport holonomy (deck group provably trivial); obstruction dichotomy + proper-map collapse control (see `docs/BF_CARICATURE.md`) | ~4 s |

Run any script from the repository root with
`.venv/bin/python scripts/<name>.py`.
