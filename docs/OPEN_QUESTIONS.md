# Open Questions — Consolidated and Ranked

*(2026-07-21. Every open question, next step, and future-work item from the
other documents, merged and de-duplicated, then sorted by ease and by the
likelihood of fast progress. Inventory compiled by a documentation sweep of
`PROBLEM.md`, `PROGRESS.md`, `QFT_IMPLICATIONS.md`, `AMPLITUDES_CONNECTION.md`,
`MONODROMY.md`, `SEARCH_STRATEGIES.md`, `NEW_COUNTEREXAMPLES.md`,
`POSITIVE_GEOMETRY.md`. Resolved items are listed at the end so nothing is
silently dropped.)*

Ranking key — the two axes the ordering combines:

- **Effort**: expected wall-clock work with the tools already in this
  repository (sympy/numpy machinery, the 2D Keller reduction, the fiber and
  trace modules).
- **P(progress)**: subjective probability that a focused session produces a
  *documentable* result (a theorem-level exact statement, a verified
  computation, or a decisive negative).

---

## Tier A — hours-to-a-day each; high probability of a clean result

**A1. Trace observables along the empty-fiber (cusp) orbit.**
[`POSITIVE_GEOMETRY.md` §5 Q3]
Over the cusp orbit the field $\phi(J)$ has no branch at all, yet the trace
observables $e_2 = q/p$, $e_3 = -r/p$ stay finite (removable behavior along
the cusp direction). Compute their exact limits and the limit of the full
pushforward algebra there; determine what "the sum over no solutions is
finite" means operatorially. *Small, sharp, exact computation with existing
modules; essentially guaranteed to terminate.*

**A2. Close out the pushforward question (AMPLITUDES_CONNECTION §2.4 Q2).**
[`AMPLITUDES_CONNECTION.md` §2.4 Q2; partial result in `POSITIVE_GEOMETRY.md` §4]
The holomorphic/real pushforward dichotomy ($-\tfrac32\,d^3J$ vs
$\tfrac{N}{2}\,d^3J$) already answers the constant-form case; what remains is
the pushforward of *non-constant* forms $g(\phi)\,d^3\phi$ (poles land on
$\{p=0\}$ by trace rationality) and a precise statement of which spurious
boundaries cancel. Mostly bookkeeping on top of `jcqft` trace machinery, then
an "Update" banner like the one Q1 received.

**A3. Witten-index / infinity prefilter for searches.**
[`SEARCH_STRATEGIES.md` §1.3, §3]
A linear-algebra test (degeneration of the leading homogeneous part on the
hyperplane at infinity) that rejects proper maps before any symbolic Keller
work. Half a day of implementation; immediately reusable by every search in
Tier B. *(Code — deferred until we decide to resume script work.)*

**A4. Larger degree boxes in the $(1,-1,-2)$ class.**
[`NEW_COUNTEREXAMPLES.md` §5 step 2]
Re-run the existing rigidity analysis in bigger boxes. Cheap because the 2D
reduction is already implemented; outcome is either "still rigid" (extends
the evidence) or a new family (a major find). *(Code — deferred.)*

---

## Tier B — days-to-weeks; realistic chance of publishable-grade progress

**B1. The Buchholz–Fredenhagen $S(f)$ caricature in 0D.**
[`QFT_IMPLICATIONS.md` §4.3(c)]
Formulate the 0D analogue of the causal-factorization relations for this $F$
and locate where the wall $\{p=0\}$, the $1\leftrightarrow 3$ sheet count,
and the $S_3$ monodromy enter the algebraic data. Pure structure, no heavy
computation; the risk is conceptual (the caricature may trivialize), which is
itself a documentable outcome.

**B2. Canonical form for the wall complement + $S_3$ local system.**
[`POSITIVE_GEOMETRY.md` §5 Q2; `MONODROMY.md`]
The real chamber failed positivity, but $\mathbb{C}^2\setminus\{P_2=0\}$ with
the rank-3 sheet local system is the natural home for log forms
($dP_2/P_2$); computing the twisted (co)homology and its intersection pairing
is standard machinery for a cuspidal cubic complement. Likely the right
"amplituhedron-analogue" object for this model.

**B3. Other $\mathbb{C}^*$ weight systems, especially $\mathbb{Z}_3$ orbits.**
[`NEW_COUNTEREXAMPLES.md` §5 step 1; `SEARCH_STRATEGIES.md` §1.1, §2]
Enumerate weight systems (e.g. $(1,-1,-3)$, $(2,-1,-3)$, 4-field gradings),
reduce each to a finite 2D Keller problem, search for residual-$\mathbb{Z}_3$
orbits (3:1 non-injectivity, possible $\mathbb{Z}_3$ monodromy). Uses the
`jcqft.reduction` template directly. Could produce the second-ever
counterexample; could also legitimately come up empty. *(Code — deferred.)*

**B4. Exact Gröbner certification of local rigidity.**
[`NEW_COUNTEREXAMPLES.md` §5 step 4, §3]
Replace float continuation by exact ideal computation for the in-box 2D
Keller system, upgrading "rigid (numerical evidence)" to a theorem. Risk:
Gröbner blow-up; mitigations (weights, truncation) are understood. *(Code —
deferred.)*

**B5. Damped 0D integrals and $\hbar \to 0$ sheet resolution.**
[`QFT_IMPLICATIONS.md` §5.3 first bullet]
Study $Z_\hbar(J) = \int e^{-|F(\phi)-J|^2/2\hbar} d^3\phi$: convergence,
Laplace asymptotics per chamber, boundary contributions from the sheets at
infinity. Numerics are easy; the interesting (and harder) part is extracting
an exact statement about the wall's signature in the $\hbar$-expansion.

**B6. Certified monodromy.**
[`MONODROMY.md` caveats]
Upgrade the numerical $S_3$ result to a certified one (interval arithmetic /
Smale alpha-theory along the loops, plus a genericity check for the chosen
line). Mechanical but somewhat laborious; the result is already believed, so
the payoff is rigor rather than novelty. *(Code — deferred.)*

---

## Tier C — weeks-to-months, research-level; progress possible but not fast

**C1. Bootstrap on inverse data.**
[`SEARCH_STRATEGIES.md` §1.2, §3; `NEW_COUNTEREXAMPLES.md` §5 step 3]
Prescribe eliminant data (perfect-square collision factor, monic
non-escaping eliminants, uniruled $\{p=0\}$) and try to *realize* a
polynomial map with that data. The realization (inverse elimination) step
has no known algorithm; genuinely hard.

**C2. Which Keller chambers are positive geometries — is the cusp forced?**
[`POSITIVE_GEOMETRY.md` §5 Q1, §3]
Nodal walls would pass, our cuspidal wall fails; deciding whether the Keller
condition in this weight system *forces* the cusp would tie the rigidity
evidence to the vertex collision. Needs either new examples (B3/C1) or a
structural argument; a good medium-term target with real depth.

**C3. Positive-geometry structure on the 3D cone itself.**
[`POSITIVE_GEOMETRY.md` §3, non-claim]
The verdict was for the C\*-reduced plane. A non-equivariant treatment of
the 3D chamber (the cone) has no off-the-shelf theory ([KPR+25] is planar;
[BD25] genus-zero pairs might apply). Theory-building required.

**C4. Higher fiber degree ($d = 4$) counterexamples.**
[`SEARCH_STRATEGIES.md` §2 item 3]
Would realize new monodromy groups ($\mathbb{Z}_4, D_4, A_4, S_4$) and
richer chamber functions than $N \in \{1,3\}$. Search space grows quickly;
depends on B3/C1 tooling.

**C5. Explicit cubic-homogeneous nilpotent (BCW) counterexample.**
[`QFT_IMPLICATIONS.md` §5.5]
Push Alpöge–Mathew through the Bass–Connell–Wright reduction and exhibit the
resulting non-invertible "$\phi^3$-type" model with nilpotent $DH$
explicitly. Conceptually mechanical, but the dimension blow-up of the
reduction makes the bookkeeping serious. High interpretive value for QFT
(BRST-like triangular vertex structure).

**C6. Explicit gradient (Hessian) counterexample — "the biggest prize".**
[`QFT_IMPLICATIONS.md` §5.3; `SEARCH_STRATEGIES.md` §2 item 1]
A non-injective Keller map $F = \nabla S$ would make partition functions,
Lefschetz thimbles, and resurgence directly well-posed. Must exist in some
dimension (de Bondt–van den Essen), but the symmetric reduction changes $n$
and no explicit example is known. Hard; enormous payoff.

**C7. Vacua at infinity in a genuine functional integral ($D \ge 1$).**
[`QFT_IMPLICATIONS.md` §4.1(iii); `AMPLITUDES_CONNECTION.md` §1.3]
Can the boundary-of-field-space mechanism survive when the measure
suppresses infinite field values? Conceptual; needs a $D \ge 1$ model class
to even pose sharply. Long-horizon.

**C8. The $n = 2$ Jacobian conjecture in QFT language.**
[`QFT_IMPLICATIONS.md` §5.4]
Can a 2-component Keller map have nonempty non-properness set? This is the
surviving Jacobian conjecture itself; treat as background motivation, not a
sprint target.

**C9. Programmatic umbrellas.**
Classification of low-degree $n=3$ counterexamples
[`QFT_IMPLICATIONS.md` §5.2] (advanced by A4/B3/B4/C1/C4); classical-map
invariants as Lagrangian data in $D \ge 1$ [`AMPLITUDES_CONNECTION.md` §1.1];
the model as a laboratory for second-type Landau singularities
[`AMPLITUDES_CONNECTION.md` §2.2]. These are directions, not single
questions; they progress when the items above do.

---

## Resolved (for the record)

| Question | Posed in | Resolution |
|---|---|---|
| Monodromy $\mathbb{Z}/3$ or $S_3$? | `QFT_IMPLICATIONS.md` §5.1 | $S_3$ (numerical, two independent lines; Galois group $S_3$ symbolically). `docs/MONODROMY.md`. Certification remains as B6. |
| Is the $N=3$ chamber a positive geometry? | `AMPLITUDES_CONNECTION.md` §2.4 Q1 | **No** — cuspidal-cubic wall, residueless double pole, vertex collision. `docs/POSITIVE_GEOMETRY.md` (2026-07-21). |
| Where does the tree expansion fail? | `PROBLEM.md` Goal 1 | Converges with finite radius to one branch of a degree-3 algebraic function; failure is non-properness, not divergence. `docs/PROGRESS.md`. |
| Circulating AI claims (zero radius, Borel, "explains $D=4$") | `PROBLEM.md`; `QFT_IMPLICATIONS.md` §3 | Corrected as false / misleading / category error. |
