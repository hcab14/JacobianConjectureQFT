# Classical-map invariants as Lagrangian data

*(2026-07-26. Tier-C9 from `docs/OPEN_QUESTIONS.md` / `docs/AMPLITUDES_CONNECTION.md` §1.1.
§§0–2: 0D dictionary (Exact/Numerical/Interpretive). §§3–5: $D\ge 1$ model-class
proposals — **no continuum computation claimed**. §6: runnable probe
`scripts/classical_map_invariants_probe.py`. Merged from parallel investigations.)*

## 0. One-sentence claim

The classical field map $F\colon\mathrm{fields}\to\mathrm{sources}$ of a
Lagrangian carries a finite list of global invariants — fiber degree,
monodromy/Galois data, non-properness divisor, real chamber function,
signed (Witten/Brouwer) index, observable-algebra defect, torus weights,
and a variationality flag — that are invisible to every local or
perturbative diagnostic yet control solution multiplicity, the radius of
perturbation theory, and the validity of field redefinitions; the
Alpöge–Mathew 0D map is the first case where all of them are computed
exactly (or, for geometric monodromy, to high numerical precision) while
$\det DF$ is a nonzero constant.

---

## 1. The 0D dictionary (exact)

**Setting.** $F\colon\mathbb{C}^3\to\mathbb{C}^3$ is the Alpöge–Mathew map
(`jcqft/core.py`), with $\det DF\equiv -2$
(`scripts/verify_counterexample.py`). Sources $J=(a,b,c)$; fields
$\phi=(x,y,z)$. Eliminant cubic
$p\,X^3+q\,X+r=0$ with
$$
p=27a^2c^2-18abc+16a+b^3c-b^2,\quad q=4-3bc,\quad r=-2c,
$$
and $D_0=27ac^2-9bc+8$ (`scripts/branch_locus.py`). Off $\{p=0\}$, $F$ is
a finite étale cover of degree 3.

**House-style tags.** **Exact** = symbolic identity asserted in a script.
**Numerical** = high-precision floating computation with documented
tolerances, not interval-certified. **Interpretive** = QFT/physics reading
of an Exact or Numerical fact; not an independent theorem.

### 1.1 Table of invariants

| Name | Algebraic definition | AM value | Verification | Local/perturbative data cannot see |
|---|---|---|---|---|
| Fiber degree $d$ | Generic $\# F^{-1}(J)$ over $\mathbb{C}$ (equiv. $[\mathbb{C}(x,y,z):\mathbb{C}(a,b,c)]$ via $F^*$) | $d=3$ | **Exact:** `scripts/verify_counterexample.py` (fiber over $(1,2,3)$); `scripts/branch_locus.py` (cubic eliminant, irreducible) | Tree expansion around $J=0$ sees one sheet; $\det DF=\mathrm{const}$ only says locally $1$-to-$1$ |
| Monodromy / Galois | Image of $\pi_1(\mathbb{C}^3\setminus\{p=0\})\to S_d$; Galois group of the $x$-cubic over $\mathbb{Q}(a,b,c)$ | Geometric monodromy $=S_3$; Galois $=S_3$ | **Exact (Galois):** `disc_X=-4D_0^2\,p` not a square (`scripts/branch_locus.py`, `docs/MONODROMY.md` §Symbolic). **Numerical (geometric):** `scripts/monodromy.py` (two lines; all $p$-loops = transpositions; group order 6). Certification deferred: OPEN_QUESTIONS B6 | Any finite jet of the local inverse; the $D_0$-collision locus (ramification of the $x$-projection only) |
| Non-properness divisor (Jelonek set) $S_F$ | $\{J:$ no nbhd has compact preimage closure$\}$; for this Keller map, leading coeff of the escaping eliminant | $S_F=\{p=0\}$ (quartic hypersurface); escape in $x$ only; $J=0\in S_F$ | **Exact:** `scripts/branch_locus.py`; escape curve + non-integrality of $x$ in `scripts/missing_observables.py` | Étale-ness ($\det DF\neq 0$); any properness test based only on $\det DF$ |
| Chamber function $N(J)$ | $\#\{\phi\in\mathbb{R}^3:F(\phi)=J\}$ | $N(J)=3$ iff $p(J)<0$; $N(J)=1$ iff $p(J)>0$; on $\{p=0\}$: $N=1$ generic wall, $N=0$ at empty-fiber cusp orbit | **Exact (off wall):** `scripts/measure_anomaly.py` (chamber rule from monic disc $-4D_0^2/p^3$; 300-target spot check). Wall fibers: **Exact** `scripts/missing_observables.py` §6. Equivalent off-wall form $N=2-\mathrm{sgn}\,p$ as in `AMPLITUDES_CONNECTION.md` §1.1 — see the §2.3 gap | Formal power-series inverse; complex fiber count (constantly $3$ off wall) |
| Witten / Brouwer index | $\deg(F,J)=\sum_{\phi\in F^{-1}(J)}\mathrm{sign}\,\det DF(\phi)$ | $\deg(F,J)=-N(J)\in\{-1,-3\}$ off wall (since $\mathrm{sign}\,\det DF\equiv -1$); jump is a non-properness certificate | **Exact:** `scripts/witten_index.py` (rational fibers per chamber + map-back). MQ mollification closed form same script; quadrature **Numerical** | Complex count; any proper-map degree theory that assumes a single integer |
| Observable-algebra defect | Cokernel of $F^*\colon\mathbb{C}[a,b,c]\hookrightarrow\mathbb{C}[x,y,z]$; field-level module structure of $\mathbb{C}(x,y,z)/\mathrm{im}\,F^*$ | $F^*$ mono not auto; $[L:K]=3$ with basis $\{1,x,x^2\}$; normal form $\mathcal{O}=c_0(F)+c_1(F)x+c_2(F)x^2$; separators carry factor $p$; $x$ not integral over $\mathrm{im}\,F^*$ | **Exact:** `scripts/missing_observables.py` (`docs/MISSING_OBSERVABLES.md`) | Surjectivity tests that only check $\det DF\neq 0$; polynomial membership of coordinates without fiber separation |
| $\mathbb{C}^*$ weights | Equivariance $F(\lambda\cdot\phi)=\lambda\cdot F(\phi)$ for a $\mathbb{C}^*$-grading | Source weights $(1,-1,-2)$ on $(x,y,z)$; target weights $(-2,-1,1)$ on $(a,b,c)$ | **Exact:** `jcqft/reduction.py` (`extract(F)` asserts the normal form $F=(P/x^2,Q/x,x\,R)$ on import; Keller residual $0$); used by `scripts/search_counterexamples.py`. Constants in `jcqft/core.py` | Untowered coefficient lists; any search that ignores grading |
| Variationality | Exists $W$ with $F=\nabla W$ (equivalently $DF=DF^{\mathsf T}$), or after affine frame change | **No** in $n=3$: $\{K:K\,DF\text{ symmetric}\,\forall\phi\}=\{0\}$. (Dim-6 cotangent / dBvdE lifts *are* variational — different map; see `docs/SYMMETRIC_SEARCH.md`) | **Exact:** `scripts/symmetric_search.py`; one-line $DF\neq DF^{\mathsf T}$ also in `scripts/witten_index.py`, `scripts/verify_counterexample.py` session notes | Existence of a formal first-order action $\bar\phi\cdot(F-J)$ (always available; not a potential for $F$) |

### 1.2 Companion facts (same dictionary, not primary rows)

| Fact | Tag | Script / doc |
|---|---|---|
| $\det DF\equiv -2$ (Keller) | **Exact** | `scripts/verify_counterexample.py` |
| Empty-fiber cusp orbit: $ac^2=4/27$, $bc=4/3$ | **Exact** | `scripts/positive_geometry.py` |
| Wall complement $\simeq A_2$-discriminant complement; $\pi_1=B_3\twoheadrightarrow S_3$ | **Exact** (affine iso) + **Numerical** (cusp-loop = 3-cycle) | `scripts/wall_braid.py`, `docs/WALL_COMPLEMENT.md` |
| Trace / pushforward rationality, poles only on $\{p=0\}$ | **Exact** | `scripts/trace_pushforward.py`, `scripts/pushforward_forms.py` |
| Measure anomaly $A(\sigma)=\langle N\rangle\to 2$ as $\sigma\to 0$ | **Numerical** (MC) on top of **Exact** chamber rule | `scripts/measure_anomaly.py` |
| Tree series = Taylor of one algebraic branch; $R\approx 0.302$ on test ray | **Exact** (cubic identity to truncation) + **Numerical** (radius) | `scripts/tree_expansion.py` |

### 1.3 What "local / perturbative" means here

For every row of §1.1, the right-hand column is scoped as follows:

- **Local:** any open in field space on which $F$ is biholomorphic / a
  $C^\infty$-diffeomorphism (guaranteed everywhere by $\det DF=\mathrm{const}\neq 0$).
- **Perturbative:** the rooted-tree formal inverse at $J=0$
  (`scripts/tree_expansion.py`), or any finite jet thereof.
- Neither sees sheets at infinity, wall-crossing of $N$ or $\deg$,
  Galois/$S_3$ permutation of vacua, the cokernel of $F^*$, or the
  failure of a potential to exist.

---

## 2. Axiomatic packaging

Proposal: attach to any polynomial (or polynomial-like) classical field
map $F\colon\mathbb{A}^n_{\mathrm{fields}}\to\mathbb{A}^n_{\mathrm{sources}}$
the following minimal invariant package. Tags mark the epistemic status
*as presently computed for AM*; the definitions themselves are algebraic /
topological.

### 2.1 Minimal list

| # | Invariant | Birational / affine? | Needs $\mathbb{R}$? | Needs orientation / volume? | Needs potential (Hessian)? | Status (AM) |
|---|---|---|---|---|---|---|
| I1 | Generic fiber degree $d$ (field extension degree) | Birational invariant of $F$ as a dominant map; unchanged by affine automorphisms of source/target composed on either side *only after* accounting for degree of those automorphisms ($=1$) | No (complex) | No | No | **Exact** |
| I2 | Galois group of a primitive element (e.g. eliminant) / geometric monodromy in $S_d$ | Birational (covers isomorphic after birational base change preserving the branch divisor) | No for Galois; geometric monodromy is over $\mathbb{C}$ | No | No | Galois **Exact**; geometric monodromy **Numerical** |
| I3 | Non-properness / Jelonek divisor $S_F$ (support + ideally the equation) | Affine-geometric (properness is not birational); transforms covariantly under affine automorphisms | No for the complex hypersurface; real points of $S_F$ need $\mathbb{R}$ | No | No | **Exact** |
| I4 | Real chamber function $N(J)$ (stratification of $\mathbb{R}^n_{\mathrm{sources}}$) | Not birational; affine-covariant over $\mathbb{R}$ | **Yes** | Uses Lebesgue class only through counting preimages (no volume form beyond the standard one for "$\#$") | No | **Exact** off wall |
| I5 | Signed Brouwer / Witten index $\deg(F,J)$ | Topological degree data; affine-covariant with sign from $\det$ of linear parts | **Yes** | **Yes** — $\mathrm{sign}\,\det DF$ needs an orientation (equivalently a volume form up to positive scale) on domain and codomain | No for the bare signed count; a *superpotential* is needed only for the Parisi–Sourlas reading (absent for AM) | **Exact** |
| I6 | Observable defect: whether $F^*$ is an automorphism; if not, module generators / extension degree of $\mathrm{Frac}(\mathrm{im}\,F^*)$ | Birational at field level; polynomial-module version is affine (finite module $\Leftrightarrow$ proper for these Keller maps) | No | No | No | **Exact** |
| I7 | Torus / $\mathbb{C}^*$ weight system (if equivariant); else "$\emptyset$" | Grading data; not a birational invariant of the underlying ungraded map (a different grading may exist after coordinate change) | No | No | No | **Exact** (weights exist and are used) |
| I8 | Variationality flag: $\exists W$ with $F=\nabla W$, optionally after linear/affine frame change; optionally coercivity of $W$ | Frame-dependent unless quantified ("exists $K$ affine s.t. …"); coercivity needs $\mathbb{R}$ | Coercivity: **Yes** | Hessian determinant uses the same volume/orientation as I5 | **Yes** — this *is* the potential datum | **Exact** (AM: no in $n=3$) |

### 2.2 Packaging rules (proposed)

1. **Always compute I1–I3, I6** for any candidate Keller / classical map —
   they are the complex-algebraic core and match `AMPLITUDES_CONNECTION.md`
   §1.1's "new layer of theory data".
2. **Add I4–I5** whenever a real structure is part of the QFT reading
   (measure anomaly, wall-crossing, SUSY localization). I5 is the signed
   refinement of I4 when $\mathrm{sign}\,\det DF$ is constant.
3. **Record I7** when the search/construction uses a grading (orbifold /
   residual-finite-stabilizer mechanism); absence is itself data.
4. **Record I8** before writing $\int e^{-W}$: if false, use
   Mathai–Quillen / first-order $\bar\phi\cdot F$ completions instead of
   Parisi–Sourlas (`docs/WITTEN_INDEX.md`, `docs/SYMMETRIC_SEARCH.md`).
5. **Interpretive layer** (not package data): "vacua at infinity",
   "Gribov without horizon", "second-type Landau singularity laboratory"
   — useful slogans, always downstream of I1–I8.

### 2.3 Honest gaps in the 0D sources

- **Geometric monodromy vs Galois.** Galois $=S_3$ is **Exact**;
  geometric monodromy $=S_3$ is **Numerical** on two lines
  (`docs/MONODROMY.md` caveats; OPEN_QUESTIONS B6). The package should
  keep these as separate slots until certification.
- **Chamber formula on the wall.** `AMPLITUDES_CONNECTION.md` §1.1 writes
  $N(J)=2-\mathrm{sgn}\,p$. Off $\{p=0\}$ this matches the chamber rule.
  On the wall, $\mathrm{sgn}\,0=0$ would give $N=2$, but actual real
  fibers have $N=1$ (generic) or $N=0$ (cusp orbit)
  (`scripts/missing_observables.py`). The value $2$ *does* appear as the
  two-sided mean / MQ wall limit $-\deg_{\mathrm{wall}}=2$
  (`scripts/witten_index.py`, `docs/DAMPED_PARTITION.md`) — do not conflate
  with $N|_{p=0}$.
- **Complex vs real degree.** Over $\mathbb{C}$, the fiber count is $3$
  for every $J\notin S_F$ (no jump). The jumping invariants I4–I5 are
  genuinely real.
- **$\mathbb{C}^*$ equivariance script surface.** The identity is asserted
  by `jcqft/reduction.py`'s `extract(F)` (import-time), not by a
  dedicated named check in `scripts/verify_counterexample.py`. Functionally
  Exact; discoverability is slightly thinner than for $p$ or $N$.
- **Variationality scope.** I8 for *this* $F$ is negative in $n=3$;
  positive examples in $n=6$ are different maps built from $F$
  (`docs/SYMMETRIC_SEARCH.md`). The dictionary entry is about the
  classical map under study, not about existence of some related action.
- **No claim that the list is complete.** Trace-pushforward pole orders,
  positive-geometry failure, and twisted periods
  (`docs/WALL_COMPLEMENT.md`, `docs/TWISTED_PERIODS.md`) are downstream
  geometry of $(S_F,\text{local system})$, not additional primary
  Lagrangian labels.

---

---

## 3. Why the $D\ge 1$ lift is nontrivial (evaporation mechanisms)

The repo is explicit that $D=0$ does **not** imply continuum statements
(`QFT_IMPLICATIONS.md` §4.2; paper outlook; `DAMPED_PARTITION.md` §3).
The mechanisms by which the 0D dictionary can *evaporate* — cease to be
well-posed, become trivial, or become undetectable — are as follows.
None of these has been computed away; each is a reason a fuzzy "$D\ge 1$
analogue" fails to be a research question.

### 3.1 Measure suppression of infinite field values

In $D=0$, escaping preimages sit at $|\phi|\to\infty$ with $F(\phi)$
finite (`MISSING_OBSERVABLES.md` escape curve; Jelonek set $\{p=0\}$).
The damped integral
$Z_\hbar(J)=\int e^{-|F-J|^2/2\hbar}\,d^3\phi$ is nevertheless finite for
**every** $J$, and the wall is seen by a piecewise-constant
semiclassical prefactor $N(J)/2$, not by divergence
(`DAMPED_PARTITION.md`). The signed Mathai–Quillen index likewise stays
finite and jumps (`WITTEN_INDEX.md`).

**Evaporation risk for $D\ge 1$.** Genuine path / field measures
dynamically suppress large field values (kinetic energy, continuum
Sobolev embeddings, lattice gradients). Open question C7 asks whether
the boundary-of-field-space mechanism can survive that suppression. A
model class that cannot even *state* a yes/no version of this is useless
for C7.

### 3.2 Discrete fibers vs moduli of classical solutions

In $D=0$, fibers $F^{-1}(J)$ are finite (degree 3 off the wall). In
$D\ge 1$, the classical field equations typically have *moduli*:
translational zero modes, gauge orbits, soliton moduli, etc. Fiber
"degree" ceases to be an integer count unless the model is cut down
(static sector, gauge-fixed slice, finite volume, discrete lattice).

**Evaporation risk.** Without a compactness / rigidity mechanism, "fiber
degree" becomes a stacky or measure-theoretic object and the chamber
function is no longer a step function of an algebraic wall.

### 3.3 Local diagnostics that falsely pass

The 0D lesson is that every *local* invertibility test can succeed while
the global covering is multi-sheeted:

| Local diagnostic | 0D status for Alpöge–Mathew |
|---|---|
| $\det DF\equiv\mathrm{const}\neq 0$ | passes ($-2$) |
| Formal / convergent inverse near vacuum | passes (tree series) |
| Hessian of $\Gamma$ nondegenerate at the perturbative vacuum | passes (no finite-distance branch point) |
| Faddeev–Popov-type determinant never zero | passes — "no Gribov horizon" |

(`AMPLITUDES_CONNECTION.md` §1.2; `QFT_IMPLICATIONS.md` §6.1 item 2.)

**Evaporation risk.** In $D\ge 1$ one can "certify" invertibility by
checking that the Fréchet derivative $D_\phi F$ is invertible on a
function space (elliptic estimate, FP operator invertible on the
fundamental modular region, etc.) and miss copies whose distinguishing
feature is escape in field-space norm — exactly the no-horizon Gribov
pattern named in the repo.

### 3.4 Effective-action multi-branching without a Maxwell flag

$\Gamma$ is built by inverting $J\mapsto\bar\phi(J)$ and Legendre
transforming. In 0D the inverse is globally 3-to-1 with branches meeting
only through infinity — no Hessian degeneration at finite distance
(`AMPLITUDES_CONNECTION.md` §1.2).

**Evaporation risk.** Continuum constructions that define $\Gamma$
perturbatively, or via a single chosen background, silently pick a
sheet. Without an independent global section-count, multi-branching is
invisible by construction.

### 3.5 Kinetic mixing destroys ultralocal product structure

An ultralocal product of 0D maps lifts the dictionary *exactly* site by
site. Any kinetic coupling (discrete Laplacian, continuum $\lvert\nabla\phi\rvert^2$)
mixes sites and can:

- glue escaping configurations into infinite-action paths,
- restore properness of the spacetime map,
- or convert discrete sheet-counts into continuous moduli.

**Evaporation risk.** The "obvious" lattice lift may be exactly 0D in
disguise ($\varepsilon=0$), or may kill the defect for every $\varepsilon>0$.
Both outcomes are scientifically useful *if posed as a dichotomy*;
neither is automatic.

### 3.6 Variational / coercivity obstruction

The honest single-field action does not exist for Alpöge–Mathew in
dimension 3 ($DF\neq DF^{\mathsf T}$; `PROBLEM.md`). The cotangent lift
$W_6=\bar\varphi\cdot F(\varphi)$ *is* a variational Keller counterexample
in dimension 6, but is necessarily non-coercive ($\kappa=-4$; affine in
$\bar\varphi$; degenerate leading form — `SYMMETRIC_SEARCH.md` §7).

**Evaporation risk.** Promoting $W_6$ by adding a $D\ge 1$ kinetic term
does not repair coercivity along conjugate directions; $\int e^{-S}$
remains ill-posed as an absolutely convergent integral. Stationary-phase
/ BRST-localized formulations can still be sharp (cf. 0D MQ), but
"partition function sees the sheets" must be asked of a *specified*
regularization, not of "$e^{-W_6}$".

### 3.7 UV / continuum category errors (explicit non-goals)

Per `QFT_IMPLICATIONS.md` §4.2 and the paper's "what the model does not
imply": nothing here bears on UV renormalization, Borel summability of
$\phi^4_{2,3}$, or $D=4$ existence/triviality. A valid $D\ge 1$ model
class for classical-map invariants is allowed to be IR / classical /
finite-volume / lattice — it must **not** be advertised as a UV result.

---

## 4. Candidate model classes (ranked)

Ranking axes (as requested):

- **(a)** How closely the 0D dictionary lifts (fiber degree, wall,
  monodromy, chamber function, signed/Witten index).
- **(b)** Whether anything is computable on a laptop (ODE BVP, small
  lattice, algebraic leftovers from 0D).
- **(c)** Honesty: what counts as a positive vs negative result.

Scores are qualitative (`+++` / `++` / `+` / `—`).

| Rank | Class | (a) lift | (b) laptop | (c) sharp ± | Verdict |
|---|---|---|---|---|---|
| 1 | Multi-component QM / static equilibria of a Keller force map | +++ | +++ | +++ | **Accept — Top-1** |
| 2 | Ultralocal lattice + controlled kinetic deformation | +++→++ | +++ | +++ | **Accept — Top-2** |
| 3 | Elliptic PDE classical maps (scalar "Gribov") | ++ | + | ++ | Accept as Tier-2 continuum target |
| 4 | Formal pAQFT / classical Møller invertibility | ++ (algebraic) | ++ (0D already done) | ++ for C9; — for C7 | Accept as algebraic avatar, not for measure |
| — | Yang–Mills Gribov horizon as *primary* lift | + (cousin) | — | + | **Reject as primary** (contrast case) |
| — | Cotangent $W_6$ + $D\ge 1$ kinetic term | + | + | + | **Reject as primary** (coercivity / first-order obstruction) |

### 4.1 Rank 1 — $D=1$ multi-component mechanics (ODE classical map)

**Accept.** Take a polynomial Keller map $F:\mathbb{R}^n\to\mathbb{R}^n$
(Alpöge–Mathew $n=3$, or the gradient $\nabla W_6$ for a variational
variant) as the *force map*. Classical static equations $F(q)=J$ are
literally the 0D fiber problem. Dynamics
$\ddot q = -F(q)+J(t)$ (or SUSY QM when $F=\nabla W$) supply a genuine
$D=1$ measure: Wiener / Feynman–Kac / Witten index.

| Axis | Assessment |
|---|---|
| (a) | Dictionary lifts *verbatim* on the static sector; monodromy of equilibria under slow $J$-cycles is the 0D $S_3$; the new content is whether the path-space measure / index sees escaping equilibria. |
| (b) | Laptop: fiber geometry is already in-repo; SUSY QM Witten index for small $n$ and polynomial $W$ is standard numerically; static chamber plots are free. |
| (c) | Positive: an index or semiclassical path-integral diagnostic that jumps across the 0D wall $\{p=0\}$ after continuum ($D=1$) measure suppression is imposed. Negative: the index / mass is $\sigma$-independent and equals the *proper* degree (no jump), i.e. escaping equilibria are measure-killed. Either outcome is publishable-grade for C7. |

Closest repo contact: `WITTEN_INDEX.md` §4 already names Witten-index
jumping in SUSY QM with non-compact targets as the structural analogy —
this class *is* that analogy made into a model.

### 4.2 Rank 2 — Ultralocal / lattice $0+1$ with per-site nonlinear map

**Accept.** On a finite lattice $\Lambda$ (periodic chain for $0+1$, or
a few sites), set
$$
\bigl(F_\varepsilon(\phi)\bigr)_x \;=\; F(\phi_x)\;-\;\varepsilon\,(\Delta_{\mathrm{disc}}\phi)_x,
$$
with $F$ the Alpöge–Mathew (or $W_6$-gradient) map and
$\varepsilon\ge 0$ a kinetic coupling.

| Axis | Assessment |
|---|---|
| (a) | At $\varepsilon=0$, the 0D dictionary is an exact product: $N_\Lambda(J)=\prod_x N(J_x)$, wall $=\bigcup_x\{p(J_x)=0\}$, monodromy $=S_3^{\lvert\Lambda\rvert}$ on generic chambers. Finite $\varepsilon$ is the first deformation that can evaporate or preserve the defect. |
| (b) | Laptop: $\lvert\Lambda\rvert\le 4$, $n=3$ is a polynomial map $\mathbb{R}^{3\lvert\Lambda\rvert}\to\mathbb{R}^{3\lvert\Lambda\rvert}$; homotopy / Gröbner / Monte-Carlo pushforward as in `measure_anomaly.py` are realistic. |
| (c) | Positive: $\exists\,\varepsilon_0>0$ such that for all $0\le\varepsilon<\varepsilon_0$ the chamber function of $F_\varepsilon$ still jumps (non-properness survives kinetic mixing). Negative: for every $\varepsilon>0$, $F_\varepsilon$ is proper (or $N$ constant) — kinetic term restores properness. The $\varepsilon\to 0$ order-of-limits question vs continuum measure is the C7 avatar on a lattice. |

### 4.3 Rank 3 — Elliptic PDE classical maps

**Accept as secondary continuum target, not Top-2.** Prototype:
$$
F(\phi) \;=\; -\Delta\phi + V'(\phi) \;=\; J
$$
on a compact Riemannian manifold (or torus), possibly multi-component
$\phi:\mathbb{T}^d\to\mathbb{R}^n$ with a Keller-type nonlinearity
replacing $V'$.

| Axis | Assessment |
|---|---|
| (a) | Fiber degree $\leftrightarrow$ finite solution count for generic $J$ (possible when the nonlinearity is proper and coercive — which Keller counterexamples *fail*); non-properness $\leftrightarrow$ sequences $\phi_k$ with $\|\phi_k\|_{H^1}\to\infty$ but $F(\phi_k)\to J_\infty$; Gribov copies = distinct solutions with invertible linearization. |
| (b) | Laptop only for severe truncations (spectral Galerkin with few modes) — essentially Rank 2. Full PDE is research-level analysis. |
| (c) | Positive: an explicit nonlinearity with $\det D_\phi F$ nowhere zero on a Banach space of fields, yet $\ge 2$ solutions for an open set of $J$, with a wall in $J$-space where solutions escape in norm. Negative: every constant-Jacobian (or elliptic-invertible) polynomial nonlinearity is proper on the chosen function space. |

Honest caveat: standard scalar $\phi^4$ is coercive and proper; the
interesting case needs a Keller-type multi-component nonlinearity, i.e.
one must *import* a 0D counterexample into the zeroth-order term. That
is well-posed but no longer "off-the-shelf PDE."

### 4.4 Rank 4 — Formal pAQFT / classical Møller map invertibility

**Accept as the algebraic avatar of C9; reject as a C7 vehicle.**

The classical Møller map of pAQFT inverts the nonlinear field equation
as a formal (or convergent) power series (`QFT_IMPLICATIONS.md`
§4.3(b); paper §pAQFT). The 0D calibration is already done: the formal
inverse can converge and still miss sheets at infinity.

| Axis | Assessment |
|---|---|
| (a) | Lifts the *algebraic* half of the dictionary (observable pullback $F^*$, non-surjectivity, Galois/$S_3$) without spacetime dynamics. Does **not** lift measure suppression, chamber functions of real solutions, or Witten-index jumping under a path measure. |
| (b) | The sharp 0D caricature is exactly B1 (Buchholz–Fredenhagen $S(f)$) plus the already-verified $F^*$ statements — laptop-scale, conceptual risk that the caricature trivializes (itself a documentable outcome). |
| (c) | Positive for C9: a formulation in which "classical-map invariants" appear as axioms/data of an algebraic Lagrangian alongside anomaly coefficients. Negative: the caricature collapses to "check properness," adding no structure. **Does not answer C7.** |

### 4.5 Reject — Yang–Mills Gribov horizon as primary model class

**Reject as the primary lift; retain as contrast.**

The repo's point (`QFT_IMPLICATIONS.md` §6.1 item 2; `PROGRESS.md`) is
**Gribov copies without a Gribov horizon**: FP-type determinant
$\equiv -2$, copies anyway, horizon retreated to infinity. Standard
Yang–Mills Gribov / Singer geometry has a *horizon* where the FP
operator ceases to be invertible — the local diagnostic *fails* at
finite field strength.

| Why reject as primary | Why keep as contact |
|---|---|
| Wrong failure mode relative to Alpöge–Mathew (horizon present vs absent). | Best-known QFT name for "many solutions of a gauge-fixing map." |
| Continuum YM is not laptop-computable; no AM-style eliminant. | Useful *contrast* sentence: AM is the no-horizon cousin; YM is the horizon cousin. |
| Gauge orbits / infinite-dimensional moduli swamp the finite-fiber dictionary. | Motivates asking for copies in models where every local FP test passes. |

### 4.6 Reject — Cotangent $W_6$ promoted by a $D\ge 1$ kinetic term

**Reject as primary; record why the temptation fails.**

`SYMMETRIC_SEARCH.md` constructs $W_6=\bar\varphi\cdot F(\varphi)$ as a
dimension-6 variational Keller counterexample (constant Hessian
determinant $-4$, explicit 3:1 witnesses). The QFT temptation is
$$
S[\bar\varphi,\varphi] \;=\; \int\!\mathrm{d}^Dx\;\Bigl(
  \tfrac12\partial_\mu\bar\varphi\cdot\partial^\mu\varphi
  + W_6(\bar\varphi,\varphi) - \bar\varphi\cdot J\Bigr).
$$

**Why this does not make sense as a primary C7/C9 model:**

1. **Coercivity is structurally impossible.** $W_6$ is affine in
   $\bar\varphi$ and has $\kappa=-4<0$ (`SYMMETRIC_SEARCH.md` §7). A
   kinetic term does not bound the action from below along conjugate
   rays; $\int e^{-S}$ is not an absolutely convergent functional
   integral for any $D$.
2. **It is still a first-order (BF-type) action.** The "potential" is
   the auxiliary-field trick of `PROBLEM.md`; promoting it in $D$
   yields a first-order field theory whose bosonic integral needs the
   same imaginary-contour / BRST / MQ treatment already used in 0D —
   not a new coercive scalar QFT.
3. **What *would* make sense** is the MQ / Parisi–Sourlas localization
   of a $D=1$ or lattice model whose *force map* is $\nabla W_6$
   (Rank 1/2 variational branch). That is already covered; calling it
   "$W_6$ QFT in $D\ge 1$" overclaims.

**Verdict:** reject the naive kinetic promotion; allow $\nabla W_6$ as
the variational force map inside Rank 1–2.

---

## 5. Top-2 formulations (sharp yes/no questions)

### 5.1 Top-1 — Multi-component mechanics with Keller force map

#### Setup

Fix a polynomial Keller map $F:\mathbb{R}^n\to\mathbb{R}^n$ with
$\det DF\equiv\kappa\neq 0$ and nonempty Jelonek set $S_F$ (running
example: Alpöge–Mathew, $n=3$, $\kappa=-2$, $S_F=\{p=0\}$). Optional
variational branch: replace $F$ by $\nabla W$ for a symmetric Keller
potential $W$ (e.g. $W_6$ on $\mathbb{R}^6$).

**Classical map (static sector).**
$$
F_{\mathrm{cl}}:\ \mathbb{R}^n\to\mathbb{R}^n,\qquad
q\;\longmapsto\; F(q).
$$
Sources $J\in\mathbb{R}^n$ label constant external forces. Equilibria
are solutions of $F(q)=J$.

**Dynamics (genuine $D=1$).** Either:

- (Bosonic) $\ddot q = -F(q) + J$, path measure formally
  $\mathcal{D}q\,\exp\bigl(-\int\tfrac12\lvert\dot q\rvert^2+V_J(q)\bigr)$
  when $F=\nabla V_J$ is variational; or
- (SUSY / MQ) the $N=2$ SUSY QM (variational case) or the
  Mathai–Quillen / Nicolai-type deformation with bosonic integral
  $$
  Z_\sigma(J)
  \;=\;
  \mathcal{N}_\sigma
  \int_{\mathrm{path}}
  \det\!\bigl(D_q F\bigr)\,
  \exp\Bigl(-\tfrac{1}{2\sigma^2}\lVert F(q)-J\rVert_{L^2}^2
  -\tfrac12\lVert\dot q\rVert_{L^2}^2\Bigr)
  $$
  (precise function space: e.g. loops of period $\beta$, or maps
  $[0,\beta]\to\mathbb{R}^n$ with fixed BCs — **part of the
  definition**, not an afterthought).

*(Interpretation flag.)* The SUSY/MQ route is the honest lift of
`WITTEN_INDEX.md`; the bosonic potential route is available only on the
variational branch and inherits non-coercivity for $W_6$.

#### Dictionary analogues

| 0D invariant | $D=1$ analogue |
|---|---|
| Fiber degree $d$ | Number of *static* equilibria $\# F^{-1}(J)$ (complex or real); for AM, $d=3$ off $\{p=0\}$. |
| Non-properness set $S_F$ | Same Jelonek set in source space: values $J$ approached by $\lvert q\rvert\to\infty$. |
| Monodromy | Permutation representation of $\pi_1(\mathbb{C}^n\setminus S_F^{\mathbb{C}})$ on the equilibrium sheet cover — for AM, the measured $S_3$. Adiabatic $J(t)$ cycles in the mechanics model realize it as vacuum monodromy. |
| Chamber function $N(J)$ | Real equilibrium count; for AM, $N=3$ ($p<0$) / $N=1$ ($p>0$) off wall (see §2.3: do not use $N=2-\mathrm{sgn}\,p$ on $\{p=0\}$). |
| Witten / signed index | $\deg(F,J)=\sum_{q:F(q)=J}\mathrm{sign}\det DF(q)$ (static); or $\mathrm{Tr}\,(-1)^F e^{-\beta H_J}$ / MQ path integral $Z_\sigma(J)$ (dynamic). In 0D these jump across $S_F$; the $D=1$ question is whether the dynamic index still jumps. |

#### Local diagnostic that can falsely pass

- $\det DF(q)=\kappa\neq 0$ for every finite $q$ (no mechanical
  "Gribov horizon" in configuration space).
- Linearization of the ODE about the perturbative equilibrium is
  nondegenerate.
- Hessian of the effective potential / pointwise Hessian of $\Gamma$
  along the perturbative branch is invertible.
- Formal gradient expansion / tree expansion for $q(J)$ converges in a
  neighborhood of $J=0$.

All of these pass for AM in 0D; none detects $N:1\leftrightarrow 3$.

#### Sharp yes/no — "vacua at infinity survive the measure"

> **Q1 (Top-1).** Fix AM (or $\nabla W_6$) as force map and a
> specified MQ / SUSY path-integral regularization $Z_\sigma(J)$
> (period $\beta$, BC class, and $\sigma>0$ fixed then $\to 0$).
> Let $J_\pm$ be regular values in the two real chambers
> ($N=3$ and $N=1$). Does
> $$
> \lim_{\sigma\to 0} Z_\sigma(J_+)
> \;\neq\;
> \lim_{\sigma\to 0} Z_\sigma(J_-)
> $$
> equal to the static signed-count jump, or are both limits equal to
> one and the same proper-degree value (escaping equilibria
> measure-killed)?

- **Yes (survive):** $\lim Z_\sigma$ reproduces $\deg(F,J_\pm)=-N(J_\pm)$
  (AM: $-3$ vs $-1$), so the $D=1$ measure still sees vacua at infinity.
  This is a *positive* resolution of C7 in this model class.
- **No (killed):** $\lim Z_\sigma$ is chamber-independent. Escaping
  static solutions do not contribute once path kinetic energy is
  present. This is a *negative* resolution of C7 for this class —
  still a theorem-shaped outcome, and it calibrates how $D=0$
  overstates the continuum.

**Status (2026-07-28) — answered for the finite-mode truncation:
SURVIVE.** `scripts/d1_index_modes.py` (37 checks, ~4.5 min; write-up
`docs/D1_INDEX.md`) computes $Z_\sigma(J; M, \beta)$ for periodic
Fourier-truncated paths (modes $\le M$, $n = 3(2M{+}1)$ real
dimensions). Exact anchors: $M=0$ reduces *exactly* to the 0D MQ
integral with $\sigma_{\rm eff} = \sigma/\sqrt\beta$ (symbolic assert),
and at any constant-path saddle the mode determinant factorizes exactly
as $\det DF(q^*)\cdot\prod_{k\le M}\lvert\det(DF(q^*) + i\omega_k
I)\rvert^2$ — a $J$-independent *positive* factor — so the $\sigma\to0$
saddle sum is $\deg(F,J) = -N(J)$ independent of $M$ and $\beta$.
Numerics (MC with error bars, seeds fixed): $M\in\{1,2\}$,
$\beta\in\{0.5,1,2\}$, $\sigma\downarrow0.05$ give $Z\to-N(J)$ per
chamber with jump ratio $\to3$ within $0.5\%$, tracking the 0D crossover
at the wall; the far-mass probes see no new escape channel in the
nonconstant modes. Honest gaps: nonconstant truncated-flow zeros are
excluded *exactly* only on the gradient branch (a 1920-start Newton
probe finds none for AM); mass from infinity in mode space is probed,
not bounded; and the order of limits $M\to\infty$ vs $\sigma\to0$ is
untouched. So the truncated model takes the **"yes (survive)"** branch
of Q1; the continuum Q1 remains open and is **not** claimed.

---

### 5.2 Top-2 — Ultralocal lattice with kinetic deformation

#### Setup

Let $\Lambda=\{1,\ldots,L\}$ be a periodic chain (or a complete graph
on $L$ sites for a mean-field kinetic term). Fields
$\phi=(\phi_x)_{x\in\Lambda}\in(\mathbb{R}^n)^L$. Fix a 0D Keller map
$F:\mathbb{R}^n\to\mathbb{R}^n$ as above and $\varepsilon\ge 0$. Define
the **lattice classical map**
$$
F_\varepsilon:
(\mathbb{R}^n)^L\to(\mathbb{R}^n)^L,
\qquad
\bigl(F_\varepsilon(\phi)\bigr)_x
\;=\;
F(\phi_x)\;-\;\varepsilon\sum_{y\sim x}(\phi_y-\phi_x).
$$
Sources $J=(J_x)\in(\mathbb{R}^n)^L$. Classical solutions:
$F_\varepsilon(\phi)=J$.

**Measure side (optional second layer).** The lattice MQ / damped
integral
$$
Z_{\hbar,\varepsilon}(J)
\;=\;
\int_{(\mathbb{R}^n)^L}
\exp\Bigl(-\tfrac{1}{2\hbar}\lvert F_\varepsilon(\phi)-J\rvert^2\Bigr)
\,d\phi
$$
is the exact finite-dimensional lift of `DAMPED_PARTITION.md`; no
continuum path space required.

#### Dictionary analogues

| 0D invariant | Lattice analogue |
|---|---|
| Fiber degree | $\# F_\varepsilon^{-1}(J)$; at $\varepsilon=0$, product of per-site degrees. |
| Non-properness set | Jelonek set $S_{F_\varepsilon}\subset(\mathbb{R}^n)^L$ (or its complexification). At $\varepsilon=0$, $\bigcup_x \pi_x^{-1}(S_F)$. |
| Monodromy | Galois / geometric monodromy of the cover defined by $F_\varepsilon$ over the complement of $S_{F_\varepsilon}$. |
| Chamber function | $N_\varepsilon(J)=\#\{\phi\in(\mathbb{R}^n)^L: F_\varepsilon(\phi)=J\}$. |
| Witten / signed index | $\deg(F_\varepsilon,J)=\sum \mathrm{sign}\det DF_\varepsilon$; mollified by $Z_{\hbar,\varepsilon}$ exactly as in `WITTEN_INDEX.md`. |

#### Local diagnostic that can falsely pass

- $\det D_\phi F_\varepsilon$ nowhere zero on $(\mathbb{R}^n)^L$ for
  small $\varepsilon$ (perturbation of $\kappa^L$).
- Sitewise Jacobian tests / mean-field Hessian of a putative $\Gamma$.
- Convergence of the coupled tree expansion about $\phi=0$.

#### Sharp yes/no — two nested questions

> **Q2a (classical, laptop-hard).** Does there exist $\varepsilon_*>0$
> such that for all $\varepsilon\in(0,\varepsilon_*)$ the chamber
> function $N_\varepsilon$ is *non-constant* on $(\mathbb{R}^n)^L$
> (equivalently: $F_\varepsilon$ is non-proper)?

- **Yes:** kinetic mixing preserves a non-properness defect for weak
  coupling — the 0D wall thickens into a lattice wall. **Positive**
  for the classical half of C9 in $D\ge 1$.
- **No:** every $\varepsilon>0$ restores properness / constant $N$.
  **Negative** — the ultralocal defect is unstable to arbitrarily weak
  kinetics. Also decisive.

> **Q2b (measure, C7 on the lattice).** Restrict to the ultralocal line
> $\varepsilon=0$ (product of 0D models) *or* to a fixed
> $\varepsilon$ for which Q2a is "yes." For sources $J$ with all sites
> in the $N=3$ chamber vs all in the $N=1$ chamber, does
> $$
> \lim_{\hbar\to 0}
> \frac{Z_{\hbar,\varepsilon}(J)}{(2\pi\hbar)^{nL/2}}
> $$
> equal $N_\varepsilon(J)/\lvert\kappa\rvert^{L}$ (pushforward
> prediction), including the contribution of configurations with some
> $\lvert\phi_x\rvert\to\infty$ as $\hbar\to 0$ near the wall?

- **Yes:** vacua at infinity survive the (lattice) measure — the
  `DAMPED_PARTITION.md` mechanism is stable under product / weak
  coupling. **Positive** for C7 in the only setting where the measure
  is rigorously finite-dimensional.
- **No:** the semiclassical prefactor collapses to a chamber-independent
  constant once $L\ge 2$ or $\varepsilon>0$. **Negative** for C7.

**Computability note / status.** Q2a for $L=1$ is the 0D theorem.
For $L=2$, $n=3$, **linear** kinetic mixing, §6 answers the
*properness half* of Q2a in the affirmative: leading forms are
unchanged, `infinity_prefilter` still fires, and an equal-mode escape
curve still hits a finite wall point
(`scripts/classical_map_invariants_probe.py`). The remaining half —
does $N_\varepsilon$ stay non-constant as a function on
$(\mathbb{R}^3)^2$ off the equal-mode slice? — is now answered **YES
(Numerical, certified reality)** for
$\varepsilon\in\{1/100,1/10,1/4,1/2,1,2\}$ by the homotopy-continuation
computation of §6.4 (`scripts/lattice_chamber.py` +
`scripts/hc_lattice_chamber.jl`); open for all $\varepsilon$
simultaneously and for $L>2$. Q2b at $\varepsilon=0$ reduces to the
product of known 0D closed forms (`DAMPED_PARTITION.md`) — a
**calibration** — and becomes interesting only for $\varepsilon>0$.

---

### 5.3 What would count as progress (honesty checklist)

| Outcome | Documents as | Does **not** license |
|---|---|---|
| Q1 "survive" for AM force + MQ paths | Positive C7 in $D=1$ mechanics | Anything about $D=4$ YM or UV |
| Q1 "killed" | Negative C7 in $D=1$; sharpens when 0D overstates | "No non-perturbative sectors in QFT" |
| Q2a yes for some $L,\varepsilon$ | Classical-map invariants exist past ultralocal | Continuum PDE statement |
| Q2a no for all $\varepsilon>0$ | Kinetics restore properness | That YM has no Gribov copies |
| Q2b yes at $\varepsilon>0$ | Measure anomaly stable to coupling | Absolute convergence of $\int e^{-W_6}$ |
| pAQFT/B1 caricature nontrivial | Progress on C9 algebraic half | Progress on C7 |
| YM horizon analysis | Contrast / motivation | A lift of the no-horizon AM defect |

---

## 6. Computational probe

Runnable certificates: `scripts/classical_map_invariants_probe.py`
(`.venv/bin/python`, 35 asserts, ~2 s; §§6.1–6.3),
`scripts/lattice_chamber.py` + `scripts/hc_lattice_chamber.jl`
(HomotopyContinuation.jl, ~3.5 min; §6.4), and
`scripts/lattice_discriminant.py` (msolve exact Gröbner/eliminants,
~6 min; §6.5).

### 6.1 What it asserts

**A. Ultralocal lattice $F^{\times N}$ (exact, $N=2,3$).**

- $\det D(F^{\times N})=(-2)^N$ (Keller).
- Real fiber cardinality multiplies: $N_{\mathrm{total}}=\prod_i N(J_i)$ over
  products of rational chamber points from `scripts/witten_index.py`.
- Non-properness set is the union of per-site walls
  $\{\prod_i p(J_i)=0\}$; one-site AM escape curves hit the wall with other
  sites at the origin.
- Witten/Brouwer index $=(-1)^N\prod_i N(J_i)$ (equivalently the product of
  per-site degrees $-N(J_i)$).

So the invariants *tensor* under ultralocal product: they remain Lagrangian
data after "adding sites."

**B. Finite-mode Galerkin + kinetic mixing (exact, $M=2$).**

Diagonal mode map $F_M(\phi_0,\phi_1)=(F(\phi_0),F(\phi_1))$ plus
$F_M^K=F_M+K\cdot\phi$ with discrete Laplacian
$K\cdot(\phi_0,\phi_1)=\varepsilon(\phi_1-\phi_0,\,\phi_0-\phi_1)$.

- Local Keller: $\det(DF_M^K)(0)\neq 0$ for small $\varepsilon$ (constant term
  $4$ at $\varepsilon=0$).
- Leading forms of $F_M^K$ equal those of $F_M$; `infinity_prefilter` still
  returns `True` — linear kinetics cannot kill the leading-form degeneration.
  (The prefilter is necessary-not-sufficient; the actual non-properness
  certificate is the next item.)
- On the equal-mode slice $\phi_0=\phi_1$, $K\cdot\phi\equiv 0$, so the AM
  escape curve — a curve in the *full* space $(\mathbb{R}^3)^2$ — still sends
  $\|\phi\|\to\infty$ to a finite image $(J_{\mathrm{wall}},J_{\mathrm{wall}})$:
  **$F_M^K$ is non-proper for every $\varepsilon$** (the witness lives on the
  diagonal; the conclusion is about the full map).
- Side effect: a $K=0$ product preimage $(\phi_0^*,\phi_1^*)$ with
  $\phi_0^*\neq\phi_1^*$ is *not* a preimage under $F_M^K$ — product-fiber
  factorization is broken by mixing.

**C. Packing API.** `classical_invariants(F, sample_Js)` returns

```text
fiber_degree_complex, chamber_N, wall_polynomial_or_None,
brouwer_index, is_gradient_at_0
```

Exercised on AM (degree $3$, wall $p$, index $-N(J)$, not a gradient) and on
the tame shear control (degree $1$, no wall, gradient at $0$).

### 6.2 What it does *not* claim

- No continuum $D\ge 1$ QFT, no renormalization, no functional determinants.
- No statement about Sobolev / distributional classical maps beyond the
  finite-dimensional Galerkin caricature.
- No reconstruction of $S_3$ monodromy for $F_M^K$ (only fiber / wall /
  Keller probes).
- Damped measures (`docs/DAMPED_PARTITION.md`) are a different layer: they
  can stay finite even when the classical map is non-proper; the probe here
  is about the *map*, not the measure.

### 6.3 Kinetic-deformation verdict

**Survives.** Linear kinetic mixing does not restore properness: leading forms
are unchanged, and an exact equal-mode escape curve still hits a finite wall
point — so $F_M^K$ is non-proper for *every* $\varepsilon$, not only in a
slice-restricted sense.  Local Keller near $0$ also survives for small
$\varepsilon$.

If one scores *all* ultralocal invariants equally, the outcome is mixed in a
side-effect sense: product-fiber factorization is washed out by mode mixing.
The primary defect proposed as Lagrangian data — the non-properness wall —
is not washed out.

### 6.4 Chamber function after mixing: Q2a answered for $L=2$

*(2026-07-26. `scripts/lattice_chamber.py` (Python wrapper, exact
$\varepsilon=0$ data + asserts) driving `scripts/hc_lattice_chamber.jl`
(HomotopyContinuation.jl v2.21); ~3.5 min default, `--full` adds
$\varepsilon\in\{1/1000,1/50,4\}$. Status: **Numerical with certified
reality** — see the scope paragraph below.)*

**Setting.** The §5.2 lattice map for $L=2$, $n=3$, in the §6.B
convention $F_\varepsilon = F_M + K\cdot\phi$,
$K=\varepsilon\begin{pmatrix}-I&I\\I&-I\end{pmatrix}$:
$$
F_\varepsilon(\phi_0,\phi_1)
=\bigl(F(\phi_0)+\varepsilon(\phi_1-\phi_0),\;
       F(\phi_1)+\varepsilon(\phi_0-\phi_1)\bigr),
\qquad N_\varepsilon(J)=\#F_\varepsilon^{-1}(J)\cap(\mathbb{R}^3)^2 .
$$
Sample targets built from the `scripts/witten_index.py` rational chamber
points: $T_1=((-1/4,0,0),(0,2,0))$ (both sites $N=3$),
$T_2=((-1/4,0,0),(1,0,0))$ (mixed), $T_3=((1,0,0),(2,1,1))$ (both
$N=1$). At $\varepsilon=0$: complex fiber $3\times3=9$, real fiber
$N(J_a)N(J_b)=9/3/1$ (**Exact**, `jcqft.fibers.exact_fiber`, asserted by
the wrapper).

**Method (three mutually cross-checking engines).** (m) A *master*
solution set over the joint parameter space $(\varepsilon,J)\in\mathbb{C}^7$:
polyhedral solve at a generic complex point plus monodromy stabilization
gives $D=66$ paths (stable under 20 no-progress monodromy loops and
under transport to an independent generic point); each rational
$(\varepsilon,J)$ is then reached by parameter homotopy along **three**
independent routes (direct + via two random complex midpoints, unioned;
deterministic per-point seeds) — this sees solutions invisible from
$\varepsilon=0$. (a) Parameter
homotopy in $\varepsilon$ from the 9 exact product solutions at
$\varepsilon=0$ (complex detour); every finite endpoint must reappear in
(m). (b) A fresh polyhedral solve at each fixed rational
$(\varepsilon,J)$; any solution outside (m) would prove (m) incomplete
(hard failure — never triggered). Reality and distinctness are
**certified** (`HC.certify`, interval arithmetic) against the exact
rational systems; conjugation parity of every fiber is asserted.

**Results** (certified real count; in parentheses: certified distinct
complex count $|$ number of real solutions descending from the
$\varepsilon=0$ product fiber):

| $\varepsilon$ | $T_1$ ($3\times3$) | $T_2$ ($3\times1$) | $T_3$ ($1\times1$) |
|---|---|---|---|
| $0$ (exact) | $9\;(9\,\vert\,9)$ | $3\;(9\,\vert\,3)$ | $1\;(9\,\vert\,1)$ |
| $1/100$ | $20\;(52\,\vert\,9)$ | $12\;(54\,\vert\,3)$ | $8\;(66\,\vert\,1)$ |
| $1/10$ | $16\;(54\,\vert\,4)$ | $4\;(54\,\vert\,1)$ | $10\;(66\,\vert\,1)$ |
| $1/4$ | $12\;(52\,\vert\,3)$ | $4\;(54\,\vert\,1)$ | $6\;(66\,\vert\,1)$ |
| $1/2$ | $14\;(46\,\vert\,4)$ | $8\;(40\,\vert\,1)$ | $6\;(60\,\vert\,1)$ |
| $1$ | $12\;(54\,\vert\,3)$ | $10\;(54\,\vert\,2)$ | $8\;(56\,\vert\,3)$ |
| $2$ | $10\;(50\,\vert\,2)$ | $8\;(44\,\vert\,1)$ | $17\;(61\,\vert\,5)$ |

**Q2a verdict (probed range): YES.** $N_\varepsilon$ is non-constant on
$(\mathbb{R}^3)^2$ at every probed $\varepsilon>0$ — wall-crossing
survives kinetic mixing; the counts even differ more violently than at
$\varepsilon=0$.

**Findings around the verdict.**

1. **$\varepsilon=0$ is a degenerate member of the family.** The joint
   generic degree is $D=66$; at $\varepsilon=0$ only $9$ solutions are
   finite (the product fiber), the rest sit at infinity. For any probed
   $\varepsilon>0$ dozens return from infinity, many of them *real*:
   at $(T_1,\varepsilon=1/100)$ the real count 20 = 9 continuations of
   the product fiber + 11 solutions that came in from infinity.
   Perturbing the kinetic term *on* rather than off is the singular
   direction — the ultralocal theory understates the solution content
   of its own neighborhood.
2. **Complex non-constancy too.** At fixed $\varepsilon$ the certified
   distinct complex count varies with $J$ (e.g. $52/54/66$ at
   $\varepsilon=1/100$) — for a proper étale map it would be constant,
   so (modulo numerical completeness) $F_\varepsilon$ stays non-proper
   over $\mathbb{C}$ at finite $\varepsilon$, consistent with the exact
   equal-mode escape certificate of §6.3. The missing paths diverge in
   the tracker (evidence of escape, not certified).
3. **$F_\varepsilon$ leaves the Keller class.** $\det DF_\varepsilon$ is
   *not* constant for $\varepsilon>0$ (**Exact**, asserted at
   $\varepsilon=1/4$: $7/4$ at the origin vs $23439371/4$ at an integer
   point; at the origin
   $\det DF_\varepsilon(0)=-2(1-2\varepsilon)(4\varepsilon^2-2)$, which
   vanishes at $\varepsilon=1/2$). So real counts can change both by
   escape through infinity *and* by finite fold bifurcations
   $\{\det DF_\varepsilon=0\}$ — the parity of the real count is no
   longer pinned to $N \bmod 2$ per chamber as in 0D, and the signed
   (Witten) count need not equal $\pm N_\varepsilon$. The signed count
   was not computed here.
4. **Wall motion in $\varepsilon$ and $J$.** Counts move with
   $\varepsilon$ at fixed $J$ (e.g. $T_1$: $20\to10$; $T_3$: $8\to17$,
   overtaking $T_1$ at $\varepsilon=2$; with `--full`,
   $\varepsilon=1/1000$ already gives $18/10/7$ and $\varepsilon=4$
   gives $8/6/18$). At fixed $\varepsilon=1/4$, bisection along the
   straight segment $T_1\to T_3$ finds certified real counts
   $12\to13\to14\to8\to4\to2\to6$ — several walls — with the first
   crossing bracketed at $t\in(27/512,\,7/128)$ (counts 12 vs 13). The
   **odd** jumps $12\to13\to14$ cannot be fold bifurcations (folds
   change the real count by $\pm2$): they are escape-type walls, i.e.
   real solutions arriving from infinity at finite $\varepsilon$ —
   the non-properness mechanism caught in the act.
   **⚠ Corrected by the exact computation of §6.5:** on this segment
   the HC counts were incomplete lower bounds at several points (exact
   $t=0$ fiber: $14$ real / $54$ complex-distinct, not $12/52$); the
   *exact* count sequence is $18\,16\,14\,16\,14\,12\,10\,8\,6\,8\,6\,
   4\,2\,4\,6$, every jump is $\pm2$ and fold-type, and the "odd jumps"
   (in particular the $12\to13$ bracket at $t\in(27/512,7/128)$) were
   completeness artifacts, not walls. Escape is real on the segment,
   but pointwise: see §6.5.
5. **Special $\varepsilon$ values.** *(Side observation, not asserted.)*
   Spot checks at generic complex $J$ gave complex counts 66 at
   $\varepsilon\in\{1/100,1/10,1/4\}$ but 60 at $\varepsilon=1/2$ and 56
   at $\varepsilon=1$ — the generic fiber degree itself appears to drop
   at special $\varepsilon$ (note $\det DF_\varepsilon(0)=0$ exactly at
   $\varepsilon=1/2$).

**Honest scope.** Certified: reality and distinctness of every reported
solution, against exact rational input ($\varepsilon$, $J$, and the
integer-coefficient map). Numerical (not certified): *completeness* of
the solution lists — it rests on the monodromy-stabilized master set,
all-paths transport along three independent routes, the fresh polyhedral
cross-solves, and conjugation parity, but no interval certificate of
exhaustiveness exists. **The completeness caveat is not hypothetical:**
the exact Gröbner computation of §6.5 shows the HC counts on the
$T_1\to T_3$ segment at $\varepsilon=1/4$ undercounted at several points
(e.g. $T_1$ itself: exact $14$ real / $54$ complex vs HC $12/52$; $T_3$
and $T_2$ were exact). Probed: $L=2$, $n=3$, six rational
$\varepsilon\in[1/100,2]$ (nine with `--full`), three rational targets
plus one bisection segment. Open: all $\varepsilon$ simultaneously
(a uniform $\varepsilon_*>0$ statement), $L>2$, the continuum limit,
signed/Witten counts for $F_\varepsilon$. The exact wall structure on
the probed segment is now computed — §6.5.

### 6.5 Exact walls on the probed segment: fold polynomial, escape polynomial

*(2026-07-28. `scripts/lattice_discriminant.py` (35 checks, ~6 min
default, `--full` adds the mod-p completeness rerun and the §7.2
fold-surface prescreen). Engines: msolve (`jcqft.gb_backend`, 16 GB
memory cap) for exact Gröbner bases, eliminants, rational
parametrizations and certified real-root isolation over $\mathbb{Q}$;
Singular for mod-$p$ factor patterns. Status: **Exact** unless labelled.
Resolves the "exact discriminant / wall-crossing / fold-vs-escape" item
of §7.2 for the probed segment.)*

**Setting.** $\varepsilon=1/4$ exactly; the segment
$J(t)=(1-t)\,T_1+t\,T_3$, $t\in[0,1]$, in the 6-dimensional source
space, with $T_1, T_3$ as in §6.4. Two exact univariate polynomials
control the walls:

**Fold polynomial $f(t)$.** The fold system
$\{F_\varepsilon(\phi)=J(t),\ \det DF_\varepsilon(\phi)=0\}$ — 7
polynomial equations in $(\phi,t)\in\mathbb{C}^7$, with
$\det DF_\varepsilon$ of total degree 19 (953 terms, computed via the
block identity
$\det\begin{pmatrix}A-\varepsilon&\varepsilon\\\varepsilon&B-\varepsilon\end{pmatrix}
=\det\bigl((A-\varepsilon)(B-\varepsilon)-\varepsilon^2\bigr)$) — is
0-dimensional of degree **516**, and $t$ is a separating element: the
exact eliminant $f(t)\in\mathbb{Z}[t]$ (msolve rational parametrization
over $\mathbb{Q}$, cross-checked against an independent mod-$p$
elimination-order Gröbner basis) has

- degree $516$, primitive, coefficient height $\le 742$ digits (full
  coefficient list regenerated by the script into
  `/tmp/lattice_discriminant/fold_poly_f_t.txt`);
- **squarefree** (gcd$(f,f')=1$ mod a good prime) and **irreducible
  over $\mathbb{Q}$** (factor patterns mod 4 primes, e.g. $31+485$ mod
  $268435367$, have common subset-sums $\{0,516\}$ only);
- exactly $50$ real roots, of which **14 lie in $(0,1)$** — and since
  $t$ is separating with a rational parametrization over $\mathbb{Q}$,
  every real root carries a *real* fold witness $\phi$:

$$t \approx 0.000435,\ 0.001873,\ 0.057820,\ 0.082786,\ 0.087906,\
0.105876,\ 0.110984,$$
$$0.121133,\ 0.124557,\ 0.126780,\ 0.156577,\ 0.300560,\ 0.511438,\
0.560931.$$

**Exact chamber counts.** msolve real-root isolation at 15 rational
sample points strictly interleaving the 14 fold roots (plus the
endpoints) gives the exact real chamber sequence at $\varepsilon=1/4$:

$$N:\quad 18 \to 16 \to 14 \to 16 \to 14 \to 12 \to 10 \to 8 \to 6 \to
8 \to 6 \to 4 \to 2 \to 4 \to 6 ,$$

with $N(1)=6$ at $T_3$ (agreeing with §6.4) and $N(0)=14$ *at* $T_1$.
Every jump is $\pm2$ and brackets exactly one fold root; conversely
every fold root in $(0,1)$ produces a jump. **On this segment every
real chamber wall is fold-type.**

**Escape polynomial $e(t)$.** The escape values (fiber degree of the
specialized ideal drops below the generic $66$) are the roots of
$$e(t) = t\,(13t-1)\,(3t+1)\,q_y(t)\,q_x(t),$$
$$q_y = 29823777\,t^4 + 5199180\,t^3 + 713782\,t^2 - 246740\,t + 12337
\quad(\text{irreducible, no real roots}),$$
$$q_x = 2841875\,t^4 + 125650\,t^3 - 1157957\,t^2 + 512672\,t - 54016
\quad(\text{irreducible, real roots}\approx -0.83010,\ 0.15629).$$
Per-factor certificates are **Exact over $\mathbb{Q}$**: ideal degree
$54<66$ at $t=0$, $53<66$ at $t=1/13$, $65<66$ at $t=-1/3$, and
$260<264=4\times66$ with either irreducible quartic adjoined (control
quartic $t^4+t+1$: exactly $264$). *Completeness* of the factor list is
certified mod $p$ for two independent 30-bit primes via the
$x_0/y_0$-eliminant leading coefficients (`--full` recomputes; the
recorded $x_1/y_1$ eliminants give the same factors); the
$z$-coordinate and generic-linear-form eliminations did not finish
within the compute budget (mod-$p$ F4, killed after $>1$ h / $>40$ min
without output — the honest wall), so a $z$-only escape channel with
additional roots cannot be excluded, though every certified escape
value above already appears in the $x/y$ directions.

**Fold-vs-escape decomposition (the §6.4 sequence corrected).**

- All 14 real chamber walls are **folds** ($\pm2$ jumps at roots of
  $f$).
- **Escape is pointwise on the reals**: at an escape value the count
  *dips at the point without changing the adjacent chambers* — exact
  examples $N(1/13)=13$ (odd!) between plateaus of $16$ ($3$ real
  solutions escape exactly at $t=1/13$, ideal degree $53$), and
  $N(0)=14$ vs $N(0^+)=18$: the endpoint $T_1$ *is* an escape point
  with $12$ of $66$ solutions at infinity, $4$ of them real. This is
  the non-properness mechanism caught exactly — but it is invisible to
  chamber counting off a measure-zero set of $t$.
- The §6.4 HC "odd jumps" $12\to13\to14$ (first crossing
  $t\in(27/512,7/128)$) were **completeness artifacts**: neither
  $f(t)$ nor $e(t)$ has a root in $(27/512,\,7/128)$, and the exact
  ideal degree at $t=109/2048$ inside the bracket is the generic $66$
  with $N=14$. The certified-lower-bound counts simply missed
  solutions there; the first genuine wall is the fold at
  $t\approx0.057820\in(29/512,\,15/256)$.

**$\varepsilon=0$ anchors (exact).** $\det DF_0\equiv 4$ (fold ideal
empty — folds are a strictly-$\varepsilon>0$ phenomenon), and the
$\varepsilon=0$ escape locus on the segment is
$p(J_0(t))\,p(J_1(t))=0$ with $p(J_0(t))=4(5t-1)$ and
$p(J_1(t))=107t^4+42t^3-85t^2+44t-4$ (one real root each in $(0,1)$:
$t=1/5$ and $t\in(29/256,30/256)$), reproducing the product chamber
sequence $9\to3\to1$. Note the contrast: at $\varepsilon=0$ all walls
are escape-type (Keller), at $\varepsilon=1/4$ all *chamber* walls on
this segment are fold-type and escape only dips.

**Open after this computation** (moved from §7.2): the fold
*hypersurface* $W(a,b,c)$ in one source block ($J_1$ frozen at
$(2,1,1)$): the mod-$p$ elimination prescreen is behind `--full`; the
exact surface over $\mathbb{Q}$ was not obtained. A wall-crossing
*formula* (predicting the $\pm$ sign of each fold jump from local data)
also remains open; the present result is the exact wall *list*, not yet
a formula.

---

## 7. Open questions refined

1. **Beyond linear $K$.** Does a nonlocal / integral kinetic operator (still
   finite-mode) create new escape channels or kill the equal-mode wall?  The
   leading-form argument covers any lower-order perturbation of a polynomial
   $F_M$; non-polynomial kinetics are open.

2. **Chamber function after mixing.** *(Answered on the probed segment —
   §§6.4–6.5.)* The real chamber function of $F_M^K$ ($M=2$) is
   non-constant at every probed $\varepsilon>0$ (§6.4), and on the
   $T_1\to T_3$ segment at $\varepsilon=1/4$ the exact
   discriminant/eliminant *is now computed* (§6.5,
   `scripts/lattice_discriminant.py`): fold polynomial $f(t)$ of degree
   516 (irreducible, squarefree, 14 real roots in $(0,1)$ = the 14
   chamber walls, all $\pm2$ fold jumps), escape polynomial
   $e(t)=t(13t-1)(3t+1)q_yq_x$ of degree 11 (escape dips the count
   pointwise, never shifts a chamber; $T_1$ itself is an escape point).
   Still open: the fold *hypersurface* $W(a,b,c)$ over a frozen source
   block (exact over-$\mathbb{Q}$ elimination out of reach; mod-$p$
   prescreen behind `--full`), a wall-crossing *formula* (sign of each
   fold jump from local data), the analogous exact structure at other
   $\varepsilon$ and on other segments, and the $z$-direction
   completeness certificate for $e(t)$ (§6.5 honest scope).

3. **Monodromy lift.** Does the $S_3$ local system of AM extend as an
   $S_3^{\times N}$ (or braid quotient) for $F^{\times N}$, and how does
   kinetic mixing entangle the factors?

4. **Variational lifts.** The cotangent / dBvdE lifts of AM
   (`scripts/symmetric_search.py`) are gradient Keller maps in dimension $6$.
   Packing their classical-map invariants under ultralocal product is a
   natural variational counterpart of §6.A.

5. **Measure vs map.** Reconcile the survival of the wall for the classical
   map (§6.B) with the unconditional finiteness of damped / Mathai–Quillen
   partition functions (`docs/DAMPED_PARTITION.md`, `docs/WITTEN_INDEX.md`):
   the measure can be well-defined while the inverse source–field map remains
   multi-branched and non-proper.

---

## 8. Reproduce

```bash
.venv/bin/python scripts/verify_counterexample.py
.venv/bin/python scripts/branch_locus.py
.venv/bin/python scripts/measure_anomaly.py
.venv/bin/python scripts/missing_observables.py
.venv/bin/python scripts/witten_index.py
.venv/bin/python scripts/symmetric_search.py          # ~3.5 min; I8
.venv/bin/python scripts/classical_map_invariants_probe.py   # ultralocal + kinetic
.venv/bin/python scripts/lattice_chamber.py           # ~3.5 min; §6.4 Q2a (needs julia + HomotopyContinuation.jl)
.venv/bin/python scripts/d1_index_modes.py            # ~4.5 min; §5.1 Q1, finite-mode truncation
.venv/bin/python scripts/lattice_discriminant.py      # ~6 min; §6.5 exact walls (needs external/msolve)
```

