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
✅ **Resolved 2026-07-21** (`scripts/cusp_traces.py`;
`POSITIVE_GEOMETRY.md` §6). The premise of the question was partly wrong:
$e_2, e_3$ *diverge* at the cusp. True exact statements: $S_1 \equiv 0$;
pole order $\lfloor k/2\rfloor$ at generic wall points; at the cusp the
divergence rates drop below naive $2k/3$ exactly when $3 \nmid k$
($\omega$-cancellation of the escaping $\mathbb{Z}_3$-triple); along the
cuspidal tangent the model is exactly solvable,
$27\varepsilon^3X^3 - 9\varepsilon X - 2 = (3\varepsilon X-2)(3\varepsilon X+1)^2$.

**A2. Close out the pushforward question (AMPLITUDES_CONNECTION §2.4 Q2).**
✅ **Resolved 2026-07-21** (`scripts/pushforward_forms.py`;
`AMPLITUDES_CONNECTION.md` §2.4 update banner). $F_*(g\,d^3\phi)$ is
rational with poles only on $\{p=0\}$ for all polynomial $g$; all
collision-locus ($D_0$) singularities cancel between sheets (verified on a
15-observable basket; forced by étaleness off the wall);
$\mathrm{ord}_p T[x^k] = \lfloor k/2\rfloor$; $F_*(d^3\phi) = -\tfrac32 d^3J$
and $F_*(x\,d^3\phi) = 0$ exactly.

**A3. Witten-index / infinity prefilter for searches.**
✅ **Resolved 2026-07-21** (`jcqft/prefilter.py`;
`scripts/witten_prefilter.py`; `docs/RIGIDITY_AND_PREFILTER.md` §1).
Maps whose leading forms have only the trivial common zero are provably
proper, hence rejected in ~1–30 ms before any symbolic Keller work (199/200
random cubic maps rejected). The Alpöge–Mathew map survives with witness
$[1:0:0]$ — its escape direction. Necessary-not-sufficient: nonlinear
automorphisms always survive (Bézout), documented as false positives.
Plain and weighted variants available for every Tier-B search.

**A4. Larger degree boxes in the $(1,-1,-2)$ class.**
✅ **Resolved 2026-07-21** (`scripts/rigidity_boxes.py`;
`docs/RIGIDITY_AND_PREFILTER.md` §2). Two strictly larger boxes (up to 74
equations × 69 unknowns): in every box the kernel of the linearized Keller
condition splits as gauge tangents ⊕ obstructed directions; **all 30
non-gauge first-order deformations across the enlarged boxes are obstructed
at second order** and nonlinear continuation confirms none integrates.
Still rigid modulo gauge — stronger evidence, no new family. (Numerical
evidence; exact certification remains B4.)

---

## Tier B — days-to-weeks; realistic chance of publishable-grade progress

**B1. The Buchholz–Fredenhagen $S(f)$ caricature in 0D.**
✅ **Resolved 2026-07-26 (split verdict)** (`scripts/bf_caricature.py`;
`docs/BF_CARICATURE.md`). The *causal-factorization* half trivializes,
provably: on a one-point spacetime causal disjointness forces
$\mathrm{supp}\,f=\emptyset$ or $\mathrm{supp}\,h=\emptyset$ and every
allowed instance of $S(f{+}g{+}h)=S(f{+}g)S(g)^{-1}S(g{+}h)$ is a
tautology (formalized + asserted). The *dynamical* half survives and is
nontrivial: with the antifield transcription (forced by I8 — no potential
in $n=3$), the relation lands every S-datum in the fiber algebra
$A_J=\mathbb{C}[x,y,z]/(F-J)$, and the invariants occupy precise slots —
the wall $\{p=0\}$ is the rank-jump locus $\dim A_J = 3\to1\to0$ of the
(non-locally-free) fiber-algebra bundle AND the pole divisor of
single-valued sector data (separators carry $p$ verbatim); the
$1\leftrightarrow3$ count is the number of characters of the real
C\*-fiber $\mathbb{C}^{N(J)}$; $S_3$ is the transport holonomy (wall
meridians = transpositions, cusp loop = order-3 Coxeter element), with
provably NO global deck action ($\mathrm{Aut}(L/K)=1$ — the cover is
non-Galois, so the "deck action" formulation fails). Obstruction
dichotomy: any BF-style $J\mapsto S(J)$ resolving sectors is multi-valued
or singular on $\{p=0\}$; anything single-valued and pole-free is blind
to sectors. Collapse control: for proper Keller maps the construction is
the trivial rank-1 bundle. Captures I1–I4, I6; misses I5, I7; I8 enters
as construction input.

**B2. Canonical form for the wall complement + $S_3$ local system.**
✅ **Resolved 2026-07-21** (`scripts/wall_braid.py`;
`docs/WALL_COMPLEMENT.md`). The wall pair is *affinely isomorphic* to the
$A_2$-discriminant pair via $(Q,R) = (w - \tfrac43,\, 2u - \tfrac{2w}{3} +
\tfrac{16}{27})$ with $4Q^3 + 27R^2 = 4P_2$ exactly; hence the wall
complement is a $K(B_3,1)$ (Arnold–Brieskorn–Deligne) and the $S_3$ sheet
monodromy is the canonical $B_3 \twoheadrightarrow S_3$. Cusp-loop monodromy
measured: a **Coxeter element** (3-cycle, order 3), image of
$\sigma_1\sigma_2$. Sheet local system $= \mathrm{triv} \oplus
\mathrm{std}$: trace rationality is the trivial summand, all
multivaluedness the reflection representation. Proposed
amplituhedron-analogue: (complement, standard local system) with twisted
periods; the concrete follow-up (twisted $H^1$ + intersection pairing) is
`WALL_COMPLEMENT.md` §6 Q3.

**B3. Other $\mathbb{C}^*$ weight systems, especially $\mathbb{Z}_3$ orbits.**
[`NEW_COUNTEREXAMPLES.md` §5 step 1; `SEARCH_STRATEGIES.md` §1.1, §2]
Enumerate weight systems (e.g. $(1,-1,-3)$, $(2,-1,-3)$, 4-field gradings),
reduce each to a finite 2D Keller problem, search for residual-$\mathbb{Z}_3$
orbits (3:1 non-injectivity, possible $\mathbb{Z}_3$ monodromy). Could
produce the second-ever counterexample; could also legitimately come up
empty. *(LARGELY RESOLVED 2026-07-22/24: exact reduced Keller identity
for the whole family $(1,-1,-m)$ — $J_2(PR^m, QR) = \kappa R^m$,
`jcqft/reduction_w.py` + `scripts/reduction_113.py`. **$(2,-1,-3)$
RESOLVED at low degree**: complete degree-1 classification, all Keller
maps tame automorphisms, $\mathbb{Z}_2$ AND $\mathbb{Z}_3$ mechanisms
provably empty (`docs/SEARCH_213.md`). **$(1,-1,-3)$ RESOLVED in the
v-linear class, all $w$-degrees**: Wronskian stratification, every
Keller map a tame automorphism, the 3:1 mechanism empty, and the
Alpöge–Mathew stratum shown to be numerologically obstructed at $m=3$
(`docs/SEARCH_113.md`) — the hoped-for $\mathbb{Z}_3$ anomaly class does
not materialize in either sibling. **UNIQUENESS THEOREM 2026-07-25**:
for EVERY $m \ge 3$ the $v$-linear class of $(1,-1,-m)$ contains only
tame automorphisms — for $m \ge 4$ with no box at all in the D strata —
so Alpöge–Mathew is the unique member of its equivariant family within
the $v$-linear class (`docs/SEARCH_11M.md`, same D3′ non-squarefree gap,
box-closed $m \le 5$). Remaining: $v$-degree $\ge 2$ ansätze, 4-field
gradings, the D3′ corner.)*

**B4. Exact Gröbner certification of local rigidity.**
[`NEW_COUNTEREXAMPLES.md` §5 step 4, §3]
Replace float continuation by exact ideal computation for the in-box 2D
Keller system, upgrading "rigid (numerical evidence)" to a theorem. Risk:
Gröbner blow-up; mitigations (weights, truncation) are understood. *(Code —
deferred.)*

**B5. Damped 0D integrals and $\hbar \to 0$ sheet resolution.**
✅ **Resolved 2026-07-21** (`scripts/damped_partition.py`;
`docs/DAMPED_PARTITION.md`). Exact closed form: $Z_\hbar(J) =
(2\pi\hbar)^{3/2}\bigl(\tfrac12 + \mathbb{P}[p(J+\sqrt{\hbar}\,\xi)<0]\bigr)$
— finite and uniformly bounded for ALL $J$, Jelonek set included (the
escaping tube has finite volume; no divergence signature). The wall's
signature is the **piecewise-constant semiclassical prefactor** $N(J)/2$,
jumping across $\{p=0\}$, equal to the two-sided mean $1$ on the wall
(vacuum $J=0$ included: twice the perturbative saddle) and to $\tfrac12 +
\kappa\hbar^{1/4}$ (anomalous exponent, $\kappa$ in closed form) at the
empty-fiber cusp. Uniformity boundary $\hbar^* \sim
\mathrm{dist}(J,\text{wall})^2$: $\gamma_{\text{wall}} = 2$,
$\gamma_{\text{cusp}} = 3 = 2\cdot\tfrac32$ (the $A_2$ horn exponent) —
set by target-space chamber geometry, not by field-space escape rates.
Sheet resolution beyond total mass remains open.

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
✅ **Resolved 2026-07-25** (`scripts/symmetric_search.py`;
`docs/SYMMETRIC_SEARCH.md`). Two explicit variational counterexamples in
dimension 6, both built from Alpöge–Mathew: the cotangent lift
$W_6 = \bar\varphi\cdot F(\varphi)$ ($\det\operatorname{Hess}W_6 \equiv -4$,
three rational witnesses) and the de Bondt–van den Essen twisted lift
$\widetilde F = \mathrm{id} + \nabla f_H$ (normalized Keller, witnesses over
$\mathbb{Q}(i)$). Complement: in dimension 3 the defect is provably
**non-variational** — no matrix $K$ makes $K\,DF$ symmetric, and the
$(1,-1,-m)$ gradient family contains only tame shears (complete for
$m \neq 3$; $m = 3$ box-closed). Coercivity remains provably out of reach
($\kappa \le 0$ forced), so the partition-function payoff is *conditional*:
an action exists, a good (coercive) action does not. **Remaining open:**
minimal dimension of a symmetric counterexample (4, 5, or 6); the $n = 3$
degree-5 box (16 GB F4 wall); the $m = 3$ slice beyond $\deg_v = 4$.

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
[`QFT_IMPLICATIONS.md` §5.2] (advanced by A4/B3/B4/C1/C4); the model as
a laboratory for second-type Landau singularities
[`AMPLITUDES_CONNECTION.md` §2.2]. These are directions, not single
questions; they progress when the items above do.

**C9′. Classical-map invariants as Lagrangian data in $D\ge 1$.**
◐ **Opened / partially resolved 2026-07-26**
(`docs/CLASSICAL_MAP_INVARIANTS.md`;
`scripts/classical_map_invariants_probe.py`). Packaged the 0D dictionary
(I1–I8) as axiomatic Lagrangian data; ranked $D\ge 1$ model classes with
sharp yes/no questions Q1 (QM path measure) and Q2a/b (lattice +
kinetics). Probe result: ultralocal product *tensors* the invariants;
**linear kinetic mixing does not wash out non-properness** (equal-mode
escape + leading-form prefilter). **Q2a answered YES for the probed
range** (`scripts/lattice_chamber.py` + `scripts/hc_lattice_chamber.jl`,
certified reality): $N_\varepsilon$ is non-constant at every probed
$\varepsilon\in[1/1000,4]$, with up to 57 of the 66 master solutions
returning from infinity at $\varepsilon>0$. **Exact walls on the
probed segment** (`scripts/lattice_discriminant.py`,
`CLASSICAL_MAP_INVARIANTS.md` §6.5): at $\varepsilon=1/4$ the fold
eliminant $f(t)$ has degree 516 with 14 real roots in $(0,1)$ = all
chamber walls (counts $18\ldots6$, every jump $\pm2$); escape
$e(t)=t(13t-1)(3t+1)q_y q_x$ only dips the count pointwise
($N(1/13)=13$; $T_1$ itself has 4 real solutions at infinity). The
§6.4 HC ``odd $\pm1$ escape jumps'' were completeness artifacts.
**Q1 answered YES for the finite-mode truncation**
(`scripts/d1_index_modes.py`; `docs/D1_INDEX.md`): the MQ index jump
survives the path measure (exact saddle factorization). Still open:
a uniform $\varepsilon_*>0$ statement, $L>2$, Q2b at $\varepsilon>0$,
signed counts for the non-Keller $F_\varepsilon$, fold hypersurface,
continuum $M\to\infty$ order of limits.

---

## Resolved (for the record)

| Question | Posed in | Resolution |
|---|---|---|
| Monodromy $\mathbb{Z}/3$ or $S_3$? | `QFT_IMPLICATIONS.md` §5.1 | $S_3$ (numerical, two independent lines; Galois group $S_3$ symbolically). `docs/MONODROMY.md`. Certification remains as B6. |
| Is the $N=3$ chamber a positive geometry? | `AMPLITUDES_CONNECTION.md` §2.4 Q1 | **No** — cuspidal-cubic wall, residueless double pole, vertex collision. `docs/POSITIVE_GEOMETRY.md` (2026-07-21). |
| Where does the tree expansion fail? | `PROBLEM.md` Goal 1 | Converges with finite radius to one branch of a degree-3 algebraic function; failure is non-properness, not divergence. `docs/PROGRESS.md`. |
| Circulating AI claims (zero radius, Borel, "explains $D=4$") | `PROBLEM.md`; `QFT_IMPLICATIONS.md` §3 | Corrected as false / misleading / category error. |
