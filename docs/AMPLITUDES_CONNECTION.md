# What the Counterexample Really Means for QFT, and the Connection to the Amplitudes Program

*(Companion to `docs/QFT_IMPLICATIONS.md` §6.3, going deeper. All toy-model
statements below are verified in this repository: `scripts/tree_expansion.py`,
`scripts/branch_locus.py`, `scripts/monodromy.py`, `scripts/measure_anomaly.py`,
`scripts/trace_pushforward.py`.)*

## 0. The one-sentence versions

**For QFT:** the classical field map of a theory carries global invariants —
fiber degree, monodromy group, non-properness divisor, chamber function —
that are invisible to perturbation theory yet control its radius of
convergence, the multiplicity of solutions, and the validity of field
redefinitions; the Alpöge–Mathew map is the first exactly solvable example
where all of these are nontrivial while every local diagnostic is trivial.

**For the amplitudes program:** the toy model realizes, as small theorems,
three of its central intuitions — diagrams are a local triangulation of a
global geometric object; the dangerous singularities live at infinity; and
summing over *all* solutions of the defining equations, not one branch,
restores rationality with poles on the physical divisor.

---

## 1. What it really means for QFT

### 1.1 A new layer of "theory data"

Standard classification data of a field theory: field content, symmetries,
couplings, anomalies. The counterexample shows there is a further layer,
attached to the *classical field map* $\phi \mapsto F(\phi)$ (equations of
motion as a map from fields to sources), which is not determined by any of
the above and not visible at any perturbative order:

| Global invariant | Value for Alpöge–Mathew | Physical meaning |
|---|---|---|
| Fiber degree $d$ | $3$ | number of solution sectors ("vacua") per source |
| Monodromy group | $S_3$ (measured) | how source cycles permute vacua |
| Non-properness divisor | $\{p=0\}$, $p = 27a^2c^2{-}18abc{+}16a{+}b^3c{-}b^2$ | where vacua escape to/return from infinity |
| Chamber function $N(J)$ | $2 - \mathrm{sgn} p$ | number of *real* solutions; wall-crossing |
| Field-space symmetry | $\mathbb{C}^*$, weights $(1,-1,-2) \to (-2,-1,1)$ | orbifold mechanism of non-injectivity |

In $D = 0$ these are algorithmically computable (this repo is the
demonstration). The programmatic proposal: treat them as *computable
invariants of a Lagrangian*, on the same footing as anomaly coefficients. For $D \ge 1$ the
analogous objects live on the solution manifold of the field equations —
where they make contact with the Gribov problem and with moduli of classical
solutions.

### 1.2 The effective action is the object that fails

The deepest practical point for QFT. The 1PI effective action $\Gamma$ is
constructed by *inverting the source–field relation* $J \mapsto \bar\phi(J)$
and Legendre transforming. In the 0D model at tree level, $\bar\phi(J)$ *is*
the inverse map we studied. The counterexample therefore says:

> The source–field relation underlying the effective-action construction can
> be globally 3-to-1 — so $\Gamma$ is multi-branched — even when its Jacobian
> is a nonzero constant, i.e. even when every local convexity/invertibility
> diagnostic passes.

Multi-branched Legendre transforms are familiar (Maxwell construction), but
there the local diagnostic (Hessian degeneration) flags the branch point. Here
there is **no finite-distance flag at all**: branches meet only through
infinity in field space. Any construction that certifies the inversion
locally — which is what formal-power-series QFT does by design — silently
selects one branch. This includes pAQFT's formal constructions; they are
correct as formal statements, and the toy model calibrates precisely what
their convergent completion would additionally have to control.

### 1.3 Non-perturbative physics without instantons

The known catalogue of non-perturbative effects — instantons, solitons,
renormalons, Stokes phenomena of asymptotic series — shares one feature:
finite-action configurations or singularities of Borel transforms, at finite
distance in configuration space. The toy model exhibits a mechanism absent
from this catalogue:

- the perturbative series **converges** (no renormalons, nothing to Borel-sum),
- there are no finite-action tunneling configurations connecting the sectors
  ($F$ is étale — sectors never meet at finite distance),
- yet the theory has extra solution sectors, connected to the perturbative
  one by monodromy **through infinity in field space**, and they dominate the
  global structure: wall-crossing in the source produces *vacuum pair
  creation from infinity* ($N: 1 \to 3$ across $\{p = 0\}$).

The measured order parameter is the anomaly factor $A(\sigma) \to 2$ of
`scripts/measure_anomaly.py`: an $O(1)$ effect for sources arbitrarily close to the
vacuum, carried by field configurations at distance $\sim \sigma^{-1/2}$.
Boundary-of-field-space effects of this type are what compactified/projective
formulations (or weighted compactifications respecting the $\mathbb{C}^*$
grading) are designed to capture; the toy model gives them an exactly
solvable home.

---

## 2. The amplituhedron connection, made precise

The amplitudes program of Arkani-Hamed and collaborators is organized around
a few structural discoveries about perturbative QFT. Each has an exact
counterpart in this model — not as analogy, but as the same mathematical
operation in lower dimension.

### 2.1 "Diagrams triangulate a geometry" — here it is a theorem

The amplituhedron [AT14] repackages the BCFW sum of diagrams as one canonical
form of one geometric object; individual diagram terms are cells of a
triangulation, and their spurious poles cancel in the sum. In the toy model
the statement is exact and proved: the rooted-tree Feynman expansion is the
Taylor series, at one point, of one branch of the algebraic covering
$\Sigma = \{(J,\phi) : F(\phi) = J\}$, a smooth affine threefold étale of
degree 3 over source space. Everything the diagrams cannot express —
the other two sheets, the $S_3$ monodromy, the escape divisor $\{p=0\}$ —
is manifest in $\Sigma$. The "geometry first, expansion second" doctrine is
here simply the difference between $\Sigma$ and its Taylor data.

### 2.2 Singularities at infinity — a solvable model of second-type Landau singularities

The singularities of amplitudes are governed by the Landau equations; their
modern algebro-geometric form is the *Landau discriminant* of Mizera–Telen
[MT22]: eliminate the internal variables from a polynomial system and ask for
which external parameters solutions degenerate. That is *literally* what
`scripts/branch_locus.py` does for the field equations: eliminate $(x,y,z)$, obtain
the eliminant cubic, and stratify its discriminant
$-4D_0^2\,p$.

The dictionary is exact:

| Amplitudes | Toy model |
|---|---|
| kinematic invariants | sources $(a,b,c)$ |
| Feynman/loop variables | fields $(x,y,z)$ |
| Landau equations | field equations $F(\phi) = J$ |
| Landau discriminant | $\mathrm{disc}_X = -4D_0^2\,p$ |
| first-type (pinch) singularities | **absent** (étale: $\det DF = -2$) |
| second-type singularities (loop momenta $\to \infty$, $\mathcal{U}_G = 0$) | escape locus $\{p = 0\}$ ($x \to \infty$) |

Second-type singularities — solutions running to infinity in loop-momentum
space — are precisely the ones that graph combinatorics handles worst and
that are "often excluded from standard discriminant-based analyses" [MT22,
FMT24]. The counterexample is a minimal, completely solvable model in which
*only* second-type singularities exist, and they do everything: set the
radius of convergence (measured: $0.3018 \approx |t| = 0.302028$ on the test
ray), carry all the monodromy ($S_3$, all transpositions), and support the
measure anomaly. If one wants a laboratory for taming singularities at
infinity, this is it.

### 2.3 Sum over all solutions — the CHY/pushforward mechanism, exactly

CHY [CHY14] computes amplitudes as sums over *all* $(n-3)!$ solutions of the
scattering equations; only the sum is rational. The canonical-form technology
of positive geometries [ABL17] uses the same operation: *pushforward* along a
finite map, i.e. the algebraic trace. In the toy model
(`scripts/trace_pushforward.py`, exact):

- Any single sheet: degree-3 algebraic, multivalued, $S_3$ monodromy.
- The trace (sum over the three sheets) of every polynomial observable is
  **rational**, with poles confined to the non-properness divisor:
  $$e_1 = 0,\qquad e_2 = \frac{q}{p},\qquad e_3 = -\frac{r}{p},\qquad
    S_k = \sum_i x_i^k \in \mathbb{Q}(a,b,c)\ \text{with denominators}\ p^{\lceil (k-1)/3 \rceil}.$$
  Equivalently: the pushforward of the field-space measure,
  $F_*\bigl(d^3\phi\bigr) = \tfrac{1}{2}\,N(J)\, d^3J$, and of all its moments,
  is rational/piecewise-constant with walls on $\{p=0\}$.
- **Boundary factorization.** Near the wall, the fiber splits into an
  escaping pair $x_\pm \sim \pm\sqrt{-q/p}$ and a finite sheet
  $x_3 \to -r/q$, and the polar part of trace observables factorizes
  Vieta-exactly:
  $$e_3 \;=\; \underbrace{(x_+x_-)}_{q/p\, +\, O(1)} \cdot
              \underbrace{x_3}_{-r/q\, +\, O(p)} .$$
  A pole whose residue factorizes into a "divergent-channel invariant" times
  an "on-wall lower amplitude" — the exact pattern of amplitude factorization
  on physical poles, here driven entirely by infinity in field space.

The moral in both fields is the same and here it is provable: *individual
solutions/branches are complicated and multivalued; the physical object is
the global sum, and its analytic simplicity (rationality, factorization) is a
consequence of the geometry of the full solution variety.* Perturbation
theory's failure in the counterexample is exactly the failure to sum over all
solutions — the tree expansion is a "one-solution CHY formula".

### 2.4 Chambers, walls, and a positive-geometry question

The real source space stratifies into chambers $\{p<0\}$ ($N=3$) and
$\{p>0\}$ ($N=1$) separated by the wall $\{p=0\}$, with the vacuum sitting on
the wall in the standard normalization (movable by translation). Sign-chamber
stratifications with wall-crossing are the raw material positive geometries
are built from (the amplituhedron itself is carved by sign conditions; its
boundaries are the physical singularities). Two well-posed questions come out
of this:

1. Is the chamber $\{p<0\}$ (or each connected component) a *positive
   geometry* in the sense of [ABL17], i.e. does it admit a canonical form
   with logarithmic poles on $\{p=0\}$ and recursive boundary structure? The
   natural candidates are built from $dp/p$ and the trace observables
   $e_2, e_3$.
2. Does the pushforward map $F_*$ send a canonical form of field space
   (trivial, $d^3\phi$, as the canonical form of the whole $\mathbb{R}^3$) to
   a canonical-form-like object on source space whose "spurious" boundaries
   cancel between sheets except on $\{p=0\}$? The measured $N(J)$ jump is the
   scalar shadow of exactly such a cancellation failure.

Neither is answered here; both are now concrete algebra problems on a
degree-3 cover of $\mathbb{C}^3$, rather than metaphors — and *any* Keller
counterexample defines such a stratification, so an answer would classify
"positive geometries of Keller maps".

> **Update (2026-07-21): question 1 is now answered — negatively, with an
> instructive mechanism.** The C\*-symmetry reduces the wall to a *cuspidal
> plane cubic* in the invariants $(u,w)=(ac^2,bc)$; the residue of every
> candidate canonical form has a residueless double pole at the cusp, and
> the chamber is **not a positive geometry** — it fails by a vertex
> collision (node $\to$ cusp), not by genus. The cusp is simultaneously the
> exact non-surjectivity locus of $F$ (empty fiber: all three sheets at
> infinity). See `docs/POSITIVE_GEOMETRY.md` and
> `scripts/positive_geometry.py`.

> **Update (2026-07-21): question 2 is now answered — affirmatively, and
> exactly** (`scripts/pushforward_forms.py`). For every polynomial
> observable $g$, $F_*\bigl(g\,d^3\phi\bigr) = -\tfrac12\,T[g]\,d^3J$ is a
> rational form whose poles lie **only** on the non-properness wall
> $\{p=0\}$: all singularities of the per-sheet values on the collision
> locus $\{D_0=0\}$ are spurious and cancel in the sheet sum (verified on
> a basket of 15 observables; conceptually forced because $F$ is finite
> étale off the wall). Pole orders obey
> $\mathrm{ord}_p\,T[x^k] = \lfloor k/2\rfloor$, wall poles are sourced
> exclusively by $x$-powers $\ge 2$ (only $x$ escapes), and two exact
> extremes hold: $F_*(d^3\phi) = -\tfrac32\,d^3J$ (constant — the wall is
> *invisible* to the holomorphic pushforward, in contrast to the real
> pushforward $\tfrac{N(J)}{2}\,d^3J$, the gap being the measured
> anomaly), and $F_*(x\,d^3\phi) = 0$ identically. The "cancellation
> between sheets except on $\{p=0\}$" asked about is thus a verified exact
> mechanism, not an analogy.

### 2.5 Honest limits of the connection

- The amplituhedron computes loop integrands of a specific interacting
  theory (planar $\mathcal{N}=4$ SYM); our model has no loops, no kinematics,
  no Grassmannian — the connection is at the level of *mechanisms*
  (triangulation vs global object, pushforward rationality, boundary
  stratification), not of objects.
- Nothing here implies positive geometry "explains" the Jacobian conjecture
  or vice versa; the overlap is that both subjects study finite maps,
  discriminants, and forms with prescribed boundary behavior — and that the
  counterexample supplies the simplest nontrivial testing ground for those
  tools outside kinematic space.
- To keep the question falsifiable, here is what would settle it either way.
  *Positive outcome:* an explicit canonical form for the chamber $\{p<0\}$
  (log poles on $\{p=0\}$, recursive boundary structure in the sense of
  [ABL17]) whose residues on the wall reproduce the boundary factorization
  verified in `scripts/trace_pushforward.py`, together with a pushforward
  identity under $F_*$ accounting for the measured $N(J)$ jump. *Negative
  outcome:* a proof that the chamber violates the positive-geometry axioms
  (e.g. its natural forms have non-logarithmic singularities on the wall).
  Both outcomes would be informative; neither is claimed here.
  *(Resolved 2026-07-21: the negative outcome holds, in exactly this form —
  the residue on the wall is a non-logarithmic double pole at the cusp of
  the reduced wall curve. See `docs/POSITIVE_GEOMETRY.md`.)*

## 3. References (addendum-specific)

- [ABL17] N. Arkani-Hamed, Y. Bai, T. Lam, *Positive geometries and
  canonical forms*, JHEP **11** (2017) 039. arXiv:1703.04541.
- [AT14] N. Arkani-Hamed, J. Trnka, *The Amplituhedron*, JHEP **10** (2014)
  030. arXiv:1312.2007.
- [CHY14] F. Cachazo, S. He, E. Y. Yuan, *Scattering of massless particles
  in arbitrary dimensions*, Phys. Rev. Lett. **113** (2014) 171601.
  arXiv:1307.2199.
- [MT22] S. Mizera, S. Telen, *Landau discriminants*, JHEP **08** (2022)
  200. arXiv:2109.08036. (Second-type singularities as solutions escaping to
  infinity; elimination-theoretic computation of singularity loci — the
  method of `scripts/branch_locus.py`.)
- [FMT24] C. Fevola, S. Mizera, S. Telen, *Landau singularities revisited:
  computational algebraic geometry for Feynman integrals*, Phys. Rev. Lett.
  **132** (2024) 101601. arXiv:2311.14669.
