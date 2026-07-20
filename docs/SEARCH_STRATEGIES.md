# Hunting New Counterexamples with QFT/Amplitudes Ideas — and What Finding One Would Mean

*(Strategy document. The mechanism facts referenced here are verified in this
repo: the $\mathbb{C}^*$-equivariance/orbifold mechanism in `PROGRESS.md`, the
discriminant structure $\operatorname{disc}_X = -4D_0^2\,p$ in
`branch_locus.py`/`monodromy.py`, the trace rationality in
`trace_pushforward.py`. The ansatz/deformation search itself is reported in
`docs/NEW_COUNTEREXAMPLES.md`.)*

## 1. Three search strategies imported from physics

### 1.1 The graded-model ansatz (from the orbifold mechanism)

The Alpöge–Mathew map is equivariant under the weighted scaling
$(x,y,z) \to (\lambda x, \lambda^{-1}y, \lambda^{-2}z)$ with target weights
$(-2,-1,1)$, and its non-injectivity is an orbifold phenomenon: the 2:1 curve
is a $\mathbb{C}^*$-orbit whose image lies in the stratum where only
even-weight target coordinates are nonzero, so $\pm\lambda$ become
indistinguishable — a residual $\mathbb{Z}_2$ acting freely on fields but
trivially on sources.

**Strategy.** Enumerate weight systems $(w_1,\dots,w_n) \to (v_1,\dots,v_n)$
for graded scalar models; write the general equivariant polynomial map
(a *finite-dimensional* ansatz per weight system, organized by the scaling
invariants such as $w = xy$); impose $\det DF = \text{const}$ symbolically;
then look for orbits whose images land in strata where the acting
$\mathbb{C}^*$ descends with a nontrivial kernel $\mathbb{Z}_k$
($k \ge 2$). Each hit is automatically $k$-to-1 on that orbit.
In physics terms: classify 0D scalar theories with a dilatation-type symmetry
whose discrete remnant is a "hidden gauge identification" on the observables.

### 1.2 The bootstrap move (geometry first, map second)

The amplituhedron/S-matrix philosophy: prescribe the global object and its
consistency conditions, then reconstruct the local model. Here, prescribe the
*inverse data* — the eliminant
$p(J)\,X^d + \dots$ of a candidate degree-$d$ covering — and impose the
consistency conditions this project proved are forced for Keller maps:

1. **Perfect-square collision factor.** Étaleness ($\det DF = $ const)
   forbids finite ramification, so every discriminant component *not* on the
   escape locus must appear to even order. (For Alpöge–Mathew:
   $\operatorname{disc}_X = -4\,D_0^2\,p$ — the collision factor $D_0^2$ is
   exactly a square.)
2. **Integrality in all but the escaping directions.** Coordinates that do
   not escape must have monic eliminants over $\mathbb{C}[J]$ (here: $y,z$
   monic, only $x$ escapes, only on $\{p=0\}$).
3. **Uniruled non-properness divisor.** By Jelonek's theory, the
   non-properness set of a polynomial map is covered by images of polynomial
   curves; candidate $p$'s must define such hypersurfaces.
4. **Vieta rationality.** Symmetric functions of the sheets must be rational
   with poles only on $\{p=0\}$ (automatic from the eliminant form; a cheap
   filter for candidate observables).

Solve these constraints for the coefficient polynomials $(p, q, r, \dots)$,
then attempt to realize the covering by an actual polynomial map (the hard
step — an inverse problem in elimination theory, but now finitely
parametrized).

### 1.3 The index heuristic (Witten-index logic)

For a real Keller map all preimages carry the same sign of $\det DF$, so the
signed solution count $\pm N(J)$ is a would-be topological invariant. It can
jump only through boundary terms at field-space infinity — the same mechanism
that makes the Witten index jump in supersymmetric quantum mechanics when
states leak to the boundary of field space. **Cheap pre-filter:** scan
candidate maps for degeneration of the leading homogeneous part on the
hyperplane at infinity (a linear-algebra condition), *before* running any
symbolic Keller check; maps whose leading part is nondegenerate at infinity
are proper and cannot be counterexamples.

## 2. What a new, inequivalent counterexample would imply for the toy QFT

Each genuinely new example (not related by composition with linear/tame
automorphisms) is a new point in a **phase diagram of non-properness
defects**, with measurable invariants: fiber degree $d$, monodromy group
$\subseteq S_d$, geometry of the escape divisor, chamber function $N(J)$,
equivariance group. The most informative possible finds, in order of impact:

1. **A gradient counterexample** ($DF$ symmetric, $F = \nabla S$). The big
   prize: it comes with an honest action, making partition functions,
   Lefschetz-thimble decompositions and resurgence questions directly
   well-posed — everything currently blocked by the non-symmetric $DF$ of
   Alpöge–Mathew. Existence in *some* dimension is guaranteed by the de
   Bondt–van den Essen symmetric reduction; no explicit one is known.
2. **Different monodromy.** Alpöge–Mathew realizes the full $S_3$, generated
   by transpositions (square-root escape, simple zeros of $p$ along the
   wall). A $\mathbb{Z}_3$ example would require cube-root escape (higher-order
   degeneration of the leading coefficient) and would change which
   sum-over-vacua observables are single-valued — effectively a different
   *global anomaly class* for these theories.
3. **Higher fiber degree / different chamber combinatorics.** $d = 3$ is
   minimal (degree 1 is Keller's birational theorem; degree-2 covers are
   Galois and étale degree-2 covers of $\mathbb{C}^n$ are trivial). A $d=4$
   example could realize $\mathbb{Z}_4$, $D_4$, $A_4$ or $S_4$, and richer
   real chamber functions than $N \in \{1,3\}$.

One normalization insight worth keeping in mind when comparing examples: the
dramatic "vacuum on the wall" feature of Alpöge–Mathew ($J=0 \in \{p=0\}$,
single vacuum, competitors at infinity) is a *translation artifact* —
composing with a translation by a preimage of the triple point puts the
vacuum inside the $N=3$ chamber, where the theory has three finite competing
perturbative branches connected by monodromy. The invariant content is the
wall and the chamber structure, not the chamber the vacuum happens to sit in.
Crossing the wall as the source is dialed is *vacuum pair creation from
infinity*.

## 3. Status

The first pass of strategy 1.1 (equivariant ansatz + first-order deformation
theory around the known solution) is reported in
`docs/NEW_COUNTEREXAMPLES.md`, including the gauge-transformation analysis
and the local-rigidity evidence within the searched ansatz box.
Strategies 1.2 and 1.3 are not yet implemented.
