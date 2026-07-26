# The 0D Buchholz–Fredenhagen S(J) Caricature: Where the Wall, the Sheets, and S3 Sit in BF-Style Algebraic Data

*(2026-07-26. Verified by `scripts/bf_caricature.py` (38 checks, ~4 s; every
Exact claim asserted there, the monodromy section labelled Numerical). This
resolves item **B1** of `docs/OPEN_QUESTIONS.md` — formulate the 0D analogue
of the Buchholz–Fredenhagen $S(f)$ relations for the Alpöge–Mathew map and
locate where the wall $\{p=0\}$, the $1\leftrightarrow 3$ real sheet count,
and the $S_3$ monodromy enter the algebraic data. Framing:
`docs/QFT_IMPLICATIONS.md` §4.3(c), paper §"Non-perturbative algebraic
constructions". Inputs: the invariant dictionary I1–I8 of
`docs/CLASSICAL_MAP_INVARIANTS.md`; the $\mathrm{im}\,F^*$ structure theorem
of `docs/MISSING_OBSERVABLES.md`; the $B_3\twoheadrightarrow S_3$ wall
complement of `docs/WALL_COMPLEMENT.md`; `docs/MONODROMY.md`;
`docs/WITTEN_INDEX.md`. Background [BF20] cited from memory — see the
disclaimer in §7.)*

**Summary of verdicts.**

1. **Causal factorization is provably vacuous in 0D** (Exact, formalized in
   the script). On a one-point spacetime, causal disjointness of supports
   forces $\mathrm{supp}\,f=\emptyset$ or $\mathrm{supp}\,h=\emptyset$, and
   in every allowed case the BF relation
   $S(f+g+h)=S(f+g)\,S(g)^{-1}\,S(g+h)$ reduces to a group-theoretic
   tautology. The one instance with content is exactly the one causal
   disjointness excludes. So the *causal* half of [BF20] trivializes — said
   cleanly, as anticipated.
2. **The dynamical relation does *not* trivialize.** Because no potential
   exists in $n=3$ ($DF\neq DF^{\mathsf T}$ — invariant I8), the only
   Lagrangian implementing $F(\phi)=J$ is the first-order one
   $L=\bar\phi\cdot(F(\phi)-J)$. The BF field-shift relation applied to the
   antifield $\bar\phi$ produces $\delta L(\beta)=\beta\cdot(F(\phi)-J)$
   *exactly* ($L$ is affine in $\bar\phi$: no higher corrections), and any
   multiplicative (pure-state) evaluation compatible with it is a character
   of the **fiber algebra** $A_J=\mathbb{C}[x,y,z]/(F-J)$. The surviving
   BF-datum in 0D is the bundle of fiber algebras over source space plus its
   parallel transport ("relative S-elements").
3. **Where the invariants enter** (the B1 question):
   the **wall** $\{p=0\}$ is (i) the jump locus of
   $\dim_{\mathbb{C}}A_J = 3\to 1\to 0$ — the sheaf of fiber algebras fails
   to be locally free exactly there — and (ii) the **pole divisor** of every
   single-valued sheet-separating datum (separator coefficients and
   would-be sector idempotents carry the factor $p$); the
   **$1\leftrightarrow 3$ count** is the number of characters (= minimal
   projections) of the real fiber C\*-algebra $\mathbb{C}^{N(J)}$; **$S_3$**
   is the holonomy of the parallel transport of the bundle — wall meridians
   are transpositions, the empty-fiber cusp loop is the order-3 Coxeter
   element (Numerical), and there is provably **no global deck action**
   (the cover is non-Galois, $\mathrm{Aut}(L/K)=1$ — Exact).
4. **Sharpest theorem-shaped statement (the obstruction dichotomy).** Any
   BF-style assignment $J\mapsto S(J)$ satisfying the 0D dynamical relation
   and resolving the classical sectors is **multi-valued** (transitive $S_3$
   holonomy; no rational section exists — Exact via irreducibility of the
   eliminant) **or singular on $\{p=0\}$** (single-valued data factors
   through the trace subalgebra, whose sheet-separating elements have poles
   exactly on the wall — Exact). Conversely, everything single-valued and
   pole-free is pulled back through $F^*$ and is blind to the sectors.
5. **Collapse control.** For a *proper* Keller map (an automorphism; tame
   shear tested) the identical construction gives the trivial rank-1 bundle
   with a global polynomial section, trivial holonomy, no wall. The
   caricature therefore does not reduce to "check properness by hand": it
   *re-expresses* non-properness as (rank jump, holonomy, pole divisor) of
   the algebraic S-data, and it collapses exactly when the defect is absent.
6. **B1 status: resolved (split verdict).** Causal half: trivializes,
   provably. Dynamical half: non-trivial, with the invariants located as in
   3–4. See §6 for the I1–I8 scorecard (captured: I1, I2, I3, I4, I6;
   missed: I5, I7; I8 enters as construction input).

---

## 0. Setting and the BF axioms in 0D

**The 0D theory.** $F:\mathbb{C}^3\to\mathbb{C}^3$ is the Alpöge–Mathew map
(`jcqft/core.py`), $\det DF\equiv-2$; classical field equation $F(\phi)=J$
with sources $J=(a,b,c)$; $x$-eliminant $p\,X^3+qX+r$ with
$p=27a^2c^2-18abc+b^3c-b^2+16a$, $q=4-3bc$, $r=-2c$;
$\mathrm{disc}_X=-4D_0^2\,p$ (asserted). $\{p=0\}$ is the Jelonek
non-properness set; $N(J)\in\{1,3\}$ off the wall.

**The BF axioms** ([BF20], from memory). Interacting dynamics is generated
by unitaries $S(f)$, $f$ ranging over local interaction functionals,
subject to

- **(i) causal factorization:** $S(f+g+h)=S(f+g)\,S(g)^{-1}\,S(g+h)$
  whenever $\mathrm{supp}\,f$ lies outside the causal past of
  $\mathrm{supp}\,h$;
- **(ii) dynamical relation:** $S(f)=S(f^\psi+\delta L(\psi))$ for field
  shifts $\psi$, where $f^\psi(\phi)=f(\phi+\psi)$ and
  $\delta L(\psi)=L(\phi+\psi)-L(\phi)$ — the off-shell implementation of
  the field equations.

In $D=0$ there is no spacetime: "spacetime" is a single point, supports are
subsets of $\{\mathrm{pt}\}$, and there is no time ordering. The question
posed by B1 is what, if anything, of (i)–(ii) survives, and where the
classical-map invariants land in the surviving data.

**House-style tags.** **Exact** = symbolic identity asserted in
`scripts/bf_caricature.py`. **Numerical** = floating computation with a
jump-guarded tracker, not certified. **Interpretive** = physics reading, not
an independent theorem.

---

## 1. Axiom (i) trivializes on a point (Exact)

On $\{\mathrm{pt}\}$ the causal past of a support is the support itself.
The disjointness condition
$\mathrm{supp}\,f\cap J^-(\mathrm{supp}\,h)=\emptyset$ therefore admits only
the pairs $(\emptyset,\emptyset)$, $(\emptyset,\{\mathrm{pt}\})$,
$(\{\mathrm{pt}\},\emptyset)$ — i.e. $f=0$ or $h=0$ (a functional with empty
support vanishes). Substituting into the relation:

| case | relation becomes | status |
|---|---|---|
| $f=0$ | $S(g+h)=S(g)\,S(g)^{-1}\,S(g+h)$ | tautology |
| $h=0$ | $S(f+g)=S(f+g)\,S(g)^{-1}\,S(g)$ | tautology |
| $f=h=0$ | $S(g)=S(g)\,S(g)^{-1}\,S(g)$ | tautology |
| $f,h$ both supported | **excluded by causal disjointness** | — |

All three tautologies are asserted with noncommutative symbols; the support
enumeration is asserted on the two-element support lattice. **The only
instance of (i) with content is the one 0D causality excludes.** This is
the clean "yes, it trivializes" for the causal half — and it is a statement
about the *degeneration of the axiom*, not about the [BF20] program (§7).

## 2. Choosing the 0D $S(J)$: three candidates

The task requires evaluating at least two formulations. All three from the
B1 brief were examined; they turn out to be facets of one object, with (c)
the correct primary and (a) provably defective.

**(a) Multiplication elements on the solution variety $\Sigma$, dynamical
relation as deck action — fails (Exact).** $\Sigma=\{(J,\phi):F(\phi)=J\}$
is the graph of $F$, isomorphic to $\phi$-space, and $\Sigma\to J$-space is
the degree-3 étale cover off the wall. The natural "dynamical action" would
be the deck group of the cover — but the cover is **non-Galois**: the
eliminant is irreducible with $\mathrm{disc}=-4D_0^2p$, and $p$ is
irreducible of multiplicity one, so the discriminant is not a square in
$\mathbb{C}(a,b,c)$, the Galois group of the closure is $S_3$, the degree-3
extension $L/K$ is not normal, and $\mathrm{Aut}(L/K)=1$ (all algebraic
ingredients asserted; the Galois-theory glue is standard and cited). **The
deck action a formulation (a) would be built on is the trivial group.**
What acts is the *monodromy on the Galois closure* — i.e. the holonomy of
formulation (c), not automorphisms of $\Sigma$ itself.

**(b) $S(J)$ from the transfer operator $T$ — the invariant part only.**
$T[\mathcal O]=\sum_{\text{sheets}}\mathcal O$
(`scripts/trace_pushforward.py`) is single-valued and rational with poles
only on $\{p=0\}$ — asserted here for $T[1]=3$, $T[x]=0$, $T[x^k]$,
$k\le 6$. But the normalized expectation $E=\tfrac13 T$ is **not
multiplicative** ($E(x^2)=-2q/3p\neq 0=E(x)^2$, asserted): $T$ retains no
product structure and no sector resolution. In local-system language
(`docs/WALL_COMPLEMENT.md` §4) it is exactly the **trivial summand** of
$\mathrm{triv}\oplus\mathrm{std}$. As the sole formulation it would make
the caricature *appear* to trivialize; it is kept as the
conditional-expectation layer of (c).

**(c) The bundle of fiber algebras with its parallel transport — chosen.**
For each $J$, the fiber algebra $A_J=\mathbb{C}[x,y,z]/(F-J)$ (functions on
the classical solutions; over $\mathbb{R}$ off the wall, the
finite-dimensional commutative C\*-algebra $\mathbb{C}^{N(J)}$ of functions
on the real solutions). The BF "net of local algebras" — which has no room
to exist in 0D — is replaced by the **bundle $\{A_J\}_J$ over source
space**; the BF-style S-data is the fiberwise unitary data together with
the **relative S-transport** $S_\gamma(J',J):A_J\to A_{J'}$, the parallel
transport of the local system along paths $\gamma$ of sources. §3 derives
this from the dynamical relation; §§4–5 locate the invariants in it.

## 3. The caricature: the dynamical relation forces the fiber algebra

**The I8 gate.** $DF\neq DF^{\mathsf T}$ (asserted; the $(1,3)$–$(3,1)$
defect equals $-1$ at $\phi=0$), so no potential $W$ with $\nabla W=F$
exists in $n=3$ (`docs/SYMMETRIC_SEARCH.md`) and axiom (ii) *cannot even be
transcribed* with a single-field Lagrangian. The honest Lagrangian is the
first-order (cotangent/antifield) one used throughout the repository
(`docs/QFT_IMPLICATIONS.md` §5.3, [Abd03]):

$$
L(\bar\phi,\phi;J)\;=\;\bar\phi\cdot\bigl(F(\phi)-J\bigr).
$$

**The 0D dynamical relation.** In 0D, "field shifts" are constant shifts.
Shifting the antifield, $\bar\phi\mapsto\bar\phi+\beta$, gives

$$
\delta L(\beta)\;=\;\beta\cdot\bigl(F(\phi)-J\bigr)
\qquad\text{exactly (asserted)},
$$

with no higher terms because $L$ is affine in $\bar\phi$ — the 0D avatar of
"the Schwinger–Dyson relation is exact, not $\hbar$-corrected", parallel to
the constant one-loop determinant. Axiom (ii) then says: the S-labels are
insensitive to adding $\beta\cdot(F(\phi)-J)$, for all $\beta\in\mathbb{C}^3$.

**Consequence for states.** A *multiplicative* evaluation of observables (a
pure-state evaluation; a character $\chi$ of $\mathbb{C}[x,y,z]$, i.e. a
point $\phi_0$ of field space) is compatible with the relation iff
$\chi(F_i)-J_i=0$ for $i=1,2,3$, i.e. iff $F(\phi_0)=J$ — and a character
that kills the generators kills the whole ideal $(F-J)$. Asserted in both
directions on explicit fiber and non-fiber points. Hence:

> **Pure classical states with source $J$ $=$ characters of
> $A_J=\mathbb{C}[x,y,z]/(F-J)$ $=$ points of the fiber $F^{-1}(J)$.**

The 0D dynamical relation, transcribed through the only available
Lagrangian, lands every S-datum in the fiber algebra. (Honest scope note:
constant shifts generate only the *linear span* of the $(F-J)_i$; the
passage to the full ideal holds for multiplicative evaluations — pure
states/sectors, which is what the sector statements below are about. Mixed
states on $A_J$ = probability weights on the fiber. See §7.)

**Relative S-elements and unitaries.** Off the wall the fibers assemble
into the rank-3 local system already computed in this repository; define
the relative S-transport $S_\gamma(J',J):A_J\to A_{J'}$ by continuation of
the solutions along $\gamma$ (over $\mathbb{R}$ within a chamber:
permutation unitaries of $\mathbb{C}^{N}$ in the character basis). The
fiberwise unitary group ($3$ phases over each off-wall $J$) is
unconstrained — there is no time ordering in 0D, so nothing fixes phases;
*all* the invariant content of the caricature sits in the transport and in
how the bundle degenerates. That is where B1's three items land.

## 4. Where the wall, the sheet count, and $S_3$ enter

### 4.1 The wall $\{p=0\}$ — rank jump and pole divisor (Exact)

Two independent appearances, both asserted:

1. **Jump locus of the fiber-algebra dimension.** Since $F$ is étale
   everywhere ($\det DF\equiv-2$), fibers are reduced and
   $\dim_{\mathbb{C}}A_J=\#F^{-1}(J)$:
   $\dim A_J=3$ at all four rational chamber points (three *distinct*
   points, exact map-back), $\dim A_J=1$ at the generic wall point
   $(\tfrac{2}{27},1,1)$ (eliminant degenerates to $qX+r$; the surviving
   character has $x=-r/q$), and $\dim A_J=0$ on the empty-fiber cusp orbit
   $(\tfrac{4}{27},\tfrac43,1)$: **$A_J$ is the zero ring — over the cusp
   orbit the caricature has no states at all.** The dimension function
   $3\to1\to0$ is non-constant, i.e. the sheaf of fiber algebras
   ($\mathbb{C}[x,y,z]$ as a module over $\mathbb{C}[a,b,c]$ via $F^*$) is
   **not locally free across $\{p=0\}$** — the module-theoretic face of
   non-properness ($x$ not integral over $\mathrm{im}F^*$,
   `docs/MISSING_OBSERVABLES.md` §2). For a finite (proper) map the
   pushforward is locally free of constant rank; here the escaping pair
   simply *leaves the algebra*, and the rank drop is exactly the escape.
2. **Pole divisor of single-valued sector data.** Anything single-valued
   that separates sheets pays in poles on $\{p=0\}$, with the wall
   polynomial appearing *verbatim*: the sheet-separating coefficients of
   $y$ and $z$ in the $\{1,x,x^2\}$ normal form carry the factor $p$
   (asserted at the level of the Gröbner data: $\mathrm{coeff}_x(B)=-6p$,
   $p\mid\mathrm{coeff}_{x^2}(B)$, $p\mid\mathrm{coeff}_{x}(D)$,
   $p\mid\mathrm{coeff}_{x^2}(D)$), and the transfer denominators are
   powers of $p$ alone ($k\le6$ asserted). Would-be sector idempotents
   (Lagrange separators through the eliminant roots) live over the wall
   complement localized at $p$ — they cannot be extended across the wall.

So the wall enters the algebraic data twice: as the **non-flat locus of
the bundle** and as the **pole divisor of the invariant hull**. It does
*not* enter as a "boundary of a net" (no net exists in 0D) nor as an ideal
of a global algebra of sections in any way that would hide it: both
appearances are forced.

### 4.2 The $1\leftrightarrow3$ count — characters of the real C\*-fiber (Exact)

Over $\mathbb{R}$ off the wall, the fiber C\*-algebra is
$C(F^{-1}(J)_{\mathbb{R}})=\mathbb{C}^{N(J)}$: the real sheet count is the
**number of characters = number of minimal projections = rank of the
identity** of the fiber algebra. Asserted: $N=3,3,1,1$ at the four chamber
points; $N=1$ at the generic wall point; $N=0$ at the cusp. The two chamber
algebras $\mathbb{C}^3$ and $\mathbb{C}^1$ are non-isomorphic, so **the
bundle of real fiber algebras is not locally trivial across the wall** —
the 0D shadow of "the algebraic data must represent the solution-count
jump" from `docs/QFT_IMPLICATIONS.md` §4.3(c). (The complex fiber algebra
has constant dimension 3 off the wall — the $1\leftrightarrow3$ jump is a
real-form phenomenon, exactly as for the Witten index,
`docs/WITTEN_INDEX.md` §6.)

### 4.3 $S_3$ — holonomy, not deck transformations

- **No global deck action (Exact).** $\mathrm{Aut}(L/K)=1$ (§2(a)): $S_3$
  cannot act by automorphisms of the solution variety over $J$-space.
- **Holonomy of the S-transport (Numerical).** Tracking the fiber around
  loops in the invariant $(u,w)$-plane (same tracker style as
  `scripts/wall_braid.py`; jump-guarded, 25 digits): the three wall
  meridians act by transpositions generating a group of order $6=S_3$; the
  $D_0$-meridian (harmless $x$-collision) acts trivially — the algebra
  transport is blind to the non-wall locus, as étaleness demands; the loop
  around the empty-fiber cusp orbit acts as a **3-cycle of order 3** — the
  Coxeter element, image of $\sigma_1\sigma_2$ under
  $B_3\twoheadrightarrow S_3$ (`docs/WALL_COMPLEMENT.md` §3). The
  permutation unitary $U_{\mathrm{cusp}}$ satisfies
  $U_{\mathrm{cusp}}^3=\mathbf 1\neq U_{\mathrm{cusp}}$ (asserted).
- **Operator meaning of the cusp loop** (Interpretive; answers
  `docs/WALL_COMPLEMENT.md` §6 Q5 within the caricature): the source cycle
  around the *total escape* point — over which the theory has no states at
  all — returns the theory with its three vacuum sectors cyclically
  rotated: "total escape = Coxeter rotation of the vacuum set" is realized
  as the monodromy unitary of the relative S-transport.

## 5. The verdict, sharply

**(i) Causal content: trivializes.** §1, Exact and final for $D=0$.

**(ii) Dynamical content: does not trivialize.** The residual relation
forces the fiber-algebra bundle, and the bundle is *obstructed*:

> **Obstruction dichotomy (theorem-shaped; Exact except the labelled
> step).** Let $J\mapsto S(J)$ be any BF-style assignment satisfying the 0D
> dynamical relation, so that its values are data in/on
> $A_J=\mathbb{C}[x,y,z]/(F-J)$. Then:
>
> 1. *(Home already sees the wall.)* $\dim_{\mathbb{C}}A_J$ takes the
>    values $3/1/0$ with jump locus exactly $\{p=0\}$; for a proper Keller
>    map it is identically $1$. Non-properness is visible before any choice
>    of $S$ is made.
> 2. *(No single-valued sector selection.)* If $S(J)$ selects a sector — is
>    character-valued or minimal-idempotent-valued — and is algebraically
>    single-valued in $J$, it does not exist: the eliminant cubic is
>    irreducible over $\mathbb{C}(a,b,c)$ (Exact: primitivity +
>    multivariate irreducibility + Gauss), so the generic fiber algebra is
>    a *field* — its only idempotents are $0$ and $1$ — and $F$ has no
>    rational section. If single-valuedness is weakened to continuity on
>    $\mathbb{C}^3\setminus\{p=0\}$, it still does not exist, by
>    transitivity of the $S_3$ holonomy (Numerical; Galois-side Exact).
> 3. *(Single-valued $\Rightarrow$ poles on the wall.)* Any single-valued
>    $S$-datum factors through the holonomy invariants — the trace/transfer
>    subalgebra — and its sheet-separating components then carry the factor
>    $p$ in their denominators: poles exactly on $\{p=0\}$ (Exact).
>
> **Hence every BF-style $J\mapsto S(J)$ that resolves the classical
> sectors is multi-valued across loops around $\{p=0\}$ or singular on
> $\{p=0\}$; and everything single-valued and pole-free is a pullback
> through $F^*$, blind to the sectors.** This is precisely the
> $\mathrm{triv}\oplus\mathrm{std}$ dichotomy of the sheet local system
> (`docs/WALL_COMPLEMENT.md` §4), now stated as a constraint on S-data.

**(iii) The caricature does not collapse to "check properness by hand".**
The collapse control (Exact): for the tame shear
$F_{\mathrm{shear}}=(x,\,z+y^2,\,y)$ — a proper Keller map — the identical
construction yields the constant rank-1 bundle with a global polynomial
section (the explicit inverse), trivial holonomy, and no pole divisor. So
the caricature's data is non-trivial *iff* the map is non-proper: it does
not merely presuppose a properness check, it **is** one, packaged as
(rank function, holonomy, pole divisor) of the algebraic S-data. What it
does not do is *reconstruct* BF dynamics — in 0D there is no time ordering
and the fiberwise phases are unconstrained; see the limits in §7.

## 6. Scorecard: I1–I8 and the B1 ledger

Explicit mapping of the outcome onto the invariant dictionary of
`docs/CLASSICAL_MAP_INVARIANTS.md` §2:

| # | Invariant | In the BF-caricature data? | Where |
|---|---|---|---|
| I1 | fiber degree $d=3$ | **captured** | generic rank of the fiber-algebra bundle, $\dim_{\mathbb{C}}A_J=3$ |
| I2 | Galois / geometric monodromy $S_3$ | **captured** | holonomy of the relative S-transport (Numerical); no-rational-section / non-normality of $L/K$ (Exact) |
| I3 | Jelonek divisor $\{p=0\}$, with its equation | **captured** | rank-jump locus of the bundle **and** pole divisor of single-valued sector data — $p$ appears verbatim in separator coefficients and transfer denominators |
| I4 | real chamber function $N(J)$ | **captured** | character count / number of minimal projections of the real fiber C\*-algebra $\mathbb{C}^{N(J)}$ |
| I5 | Witten / Brouwer index $-N(J)$ | **missed** | the commutative fiber algebra is unoriented; $\mathrm{sign}\det DF$ requires the orientation/volume datum of the MQ completion (`docs/WITTEN_INDEX.md`), which is not part of the BF-style data here |
| I6 | observable defect ($\mathrm{im}F^*$, basis $\{1,x,x^2\}$) | **captured — it *is* the caricature** | the bundle is the globalization of the $\mathrm{im}F^*$-module structure; holonomy invariance $=$ membership in (the fraction field of) $\mathrm{im}F^*$ |
| I7 | $\mathbb{C}^*$ weight system | **missed** | nothing in the fiber-algebra data references the grading; equivariance can be *imposed* as extra structure but is not encoded |
| I8 | variationality flag | **input, not output** | the failure of I8 forces the antifield transcription of the dynamical relation (§3); the caricature records I8 at the construction step |

**B1 ledger.** *Resolved 2026-07-26 (split verdict).* The causal half of
the $S(f)$ relations trivializes in 0D, provably and formalizably (§1). The
dynamical half survives and is non-trivial: it forces the bundle of fiber
algebras, in which wall / sheet count / $S_3$ occupy the precise slots of
§4, with the obstruction dichotomy of §5 as the sharpest statement. The
"risk that the caricature trivializes" materialized for exactly one axiom
and one candidate formulation ((a): deck group provably trivial), and both
are documented as such.

## 7. Honest limits

- **This is a 0D caricature, not a statement about [BF20].** The
  Buchholz–Fredenhagen program lives in $D\geq 1$, where causal
  factorization is its engine; nothing here bears on its correctness,
  power, or constructions. What is shown is how *its axioms degenerate* on
  a point and which slots of the surviving algebraic skeleton the
  Alpöge–Mathew invariants occupy. No net, no causality, no Hilbert-space
  dynamics, no vacuum state, no superselection theory exists in 0D.
- **[BF20] is cited from memory** (unitaries $S(f)$; causal factorization;
  dynamical/field-shift relation). The 0D transcription of the dynamical
  relation — antifield shifts of the first-order Lagrangian — is *our*
  principled choice, forced by I8; other transcriptions were evaluated
  (§2) but the list is not exhaustive.
- **Span vs ideal.** Constant antifield shifts generate the linear span of
  the $(F-J)_i$, not the ideal. The step to $A_J$ is exact for
  *multiplicative* evaluations (characters kill an ideal once they kill
  its generators — this covers all sector statements). For general mixed
  states the constraint is weaker; no claim is made there.
- **Unitarity is bookkeeping here.** The fiberwise unitary groups are
  unconstrained (no time ordering in 0D); "S-unitaries" carry content only
  through the permutation part of the transport. Calling the transport
  matrices "unitaries" is Interpretive dressing on an Exact/Numerical
  permutation statement.
- **Geometric monodromy remains Numerical** (jump-guarded tracking at 25
  digits, same caveats as `docs/MONODROMY.md`; certification is
  OPEN_QUESTIONS B6). The Galois-side statements (irreducibility,
  non-normality, trivial deck group) are Exact.
- **Fiber dimensions on the wall** are counted through the exhaustive
  finite-fiber description (`docs/MISSING_OBSERVABLES.md` §4) plus
  étaleness (reduced fibers); the escaped sheets leave *no trace inside*
  $A_J$ — their memory sits in the rank drop and in the $p$-poles of the
  invariant hull, which is exactly the point, but it means $A_{J_{\rm
  wall}}$ alone cannot reconstruct what escaped.
- **One map, $D=0$.** As with every document in this repository: no claim
  about $D\geq1$, about other Keller maps beyond the tested control, or
  about UV physics of any kind.

## 8. Reproduce

```bash
.venv/bin/python scripts/bf_caricature.py   # 38 checks, ~4 s
```

## References

- [BF20] D. Buchholz, K. Fredenhagen, *A C\*-algebraic approach to
  interacting quantum field theories*, Comm. Math. Phys. **377** (2020)
  947–969. arXiv:1902.06062. *(Cited from memory; see §7.)*
- [Abd03] A. Abdesselam, *The Jacobian conjecture as a problem of
  perturbative quantum field theory*, Ann. Henri Poincaré **4** (2003)
  199–215. (First-order / conjugate-field representation used in §3.)
- Repo-internal: `docs/QFT_IMPLICATIONS.md` §4.3, `docs/OPEN_QUESTIONS.md`
  B1, `docs/CLASSICAL_MAP_INVARIANTS.md` (I1–I8),
  `docs/MISSING_OBSERVABLES.md`, `docs/WALL_COMPLEMENT.md`,
  `docs/MONODROMY.md`, `docs/WITTEN_INDEX.md`,
  `scripts/bf_caricature.py` (all verifications for this document).
