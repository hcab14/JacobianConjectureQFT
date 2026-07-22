# Progress log

All results below are produced by scripts in this repository (exact symbolic
computation with SymPy unless stated otherwise) and can be reproduced with
the commands given. Environment: `.venv` (Python 3.12, `requirements.txt`).

## 2026-07-20 — Verification of the counterexample

`scripts/verify_counterexample.py`

- $\det DF \equiv -2$ exactly (symbolic determinant).
- $F(0,0,-\tfrac14) = F(1,-\tfrac32,\tfrac{13}{2}) = F(-1,\tfrac32,\tfrac{13}{2})
  = (-\tfrac14, 0, 0)$: $F$ is not injective. **Jacobian conjecture false for
  $n \ge 3$.**
- Generic fiber (e.g. over $(1,2,3)$) has **3 points**: $F$ is generically
  3-to-1; the associated field extension is cubic.
- Structural observation: with $g = 1+xy$ and $u = g^2 z + y^2(4+3xy)$,
  $$F = \bigl(g\,u,\;\; y + 3x\,u,\;\; x\,(2 - 3xy - x^2 z)\bigr).$$
- $DF$ is **not symmetric**: $F$ is not a gradient field, so there is no
  potential $S$ with $\nabla S = F$; QFT formulations must use the
  source/auxiliary-field form, not a naive action.

## 2026-07-20 — Exact algebra of the inverse (Gröbner elimination)

`scripts/branch_locus.py`

- **x-eliminant.** Every preimage of a target $(a,b,c)$ has $x$-coordinate
  satisfying the cubic
  $$p\,X^3 + q\,X + r = 0,\qquad
    p = 27a^2c^2 - 18abc + 16a + b^3c - b^2,\quad q = 4-3bc,\quad r = -2c.$$
- **Rational fiber parametrization.** Given a root $x$, the lex Gröbner basis
  yields $y$ and $z$ as rational functions of $(x,a,b,c)$ with denominator
  $D_0 = 27ac^2 - 9bc + 8$ (independent of $x$).
- **$y$ and $z$ never escape.** Their minimal cubics over $\mathbb{C}[a,b,c]$
  are (up to constants) *monic* — $y, z$ are integral over the target
  coordinates. Escape to infinity happens **only in $x$**, exactly on the
  hypersurface $\{p = 0\}$ (the Jelonek non-properness set of $F$, in the
  $x$-direction).
- **Discriminant classified.** $\operatorname{disc}_X = -p\,(4q^3 + 27 p r^2)$.
  - $\{p = 0\}$: true branch locus — one sheet at infinity; the fiber drops
    from 3 to $\le 2$ finite points; monodromy can occur around it.
  - $\{4q^3 + 27 p r^2 = 0,\ p \neq 0\}$: two *distinct* fiber points merely
    share an $x$-coordinate (verified at the rational point
    $(a,b,c) = (\tfrac1{27},1,1)$: fiber
    $\{(3,\cdot,\cdot),(3,\cdot,\cdot),(-6,\tfrac56,\tfrac{103}{216})\}$,
    3 distinct points). The covering is unramified there ($\det DF \neq 0$
    forbids merging); the $x$-*projection* ramifies, the covering does not.
- **The origin lies on the branch locus:** $p(0,0,0) = 0$. Over the
  perturbative vacuum $J = 0$ the fiber is a *single* point $(0,0,0)$ — two
  of the three sheets already sit at infinity.
- **Vacua at infinity, exactly.** The curve
  $$\varphi(s) = \bigl(s,\; -\tfrac{3}{2s},\; \tfrac{13}{2s^2}\bigr),\qquad
    F(\varphi(s)) = \bigl(-\tfrac{1}{4s^2},\, 0,\, 0\bigr)$$
  satisfies $\varphi(s) \neq \varphi(-s)$ with equal images: $F$ is exactly
  2-to-1 on this curve. Along the segment $J = (a,0,0)$, $a \to 0^-$, the two
  non-perturbative sheets are $\varphi(\pm 1/(2\sqrt{-a}))$ and recede to
  infinity like $(-a)^{-1/2}$; at $a = -\tfrac14$ they land at
  $(\pm 1, \mp\tfrac32, \tfrac{13}{2})$ — the famous triple point. They are
  invisible to *every* order of perturbation theory around $J = 0$.

## 2026-07-20 — Tree-graph (perturbative) expansion

`scripts/tree_expansion.py` (sparse truncated-ring arithmetic; Picard iteration =
sum over rooted tree Feynman graphs; runtime ~1 min)

- Formal inverse computed to total order 10 in $(a,b,c)$;
  $F(G(J)) = J$ verified to that order. Low orders:
  $x = \tfrac{c}{2} + \tfrac{3bc^2}{8} + \dots$, $y = b - \tfrac{3ac}{2} + \dots$,
  $z = a - 4b^2 + \tfrac{21abc}{2} + \dots$
- **Every order $1..10$ carries nonzero coefficients** — the series never
  terminates; there is no polynomial inverse (consistent with non-injectivity).
- The $x$-series **satisfies the exact cubic** $p x^3 + q x + r = 0$ to order
  10: perturbation theory is computing the Taylor series of a *degree-3
  algebraic function* — the most tractable "non-perturbative completion"
  imaginable.
- **Finite, nonzero radius of convergence** (refuting the "zero radius /
  asymptotic series" claim). On the ray $J = t\,(1,2,3)$, series to $t^{60}$:
  Domb–Sykes extrapolation gives $R \approx 0.3018$ for all three components,
  matching the nearest zero $|t| = 0.302028$ of $p(t)$ — **not** the nearest
  zero $|t| = 0.261013$ of the other discriminant factor.
- Newton path-tracking confirms the mechanism: continuing the perturbative
  branch toward $t = -0.261$ it stays finite (that singularity belongs to
  other sheets' $x$-collision); toward $t = -0.302$ it **escapes to
  infinity**. The radius of convergence of the tree expansion is set by the
  *escape-to-infinity (non-properness) locus*, not by sheet collisions and
  not by combinatorial growth of tree graphs.
- On the special ray $J = t\,(-1,0,0)$ (through the triple point) the series
  *terminates*: the perturbative branch is exactly $\varphi_A = (0,0,-t)$ and
  analytically continues to the preimage $(0,0,-\tfrac14)$ — while the other
  two preimages approach from infinity, unseen.

**QFT reading so far.** In this toy model perturbation theory is convergent
and correct — but only for one local branch. The "non-perturbative effects"
are not exponentially small corrections to the series; they are entire
solution sectors located at infinity in field space, connected to the
perturbative sector only through the boundary. Global invertibility (the
conjecture) fails through non-properness, while every local/perturbative
diagnostic (constant Jacobian, convergent trees) looks perfectly healthy.

## Corrections to circulating AI-generated claims

1. *"The tree expansion has zero radius of convergence / is asymptotic
   because trees grow like $N!$."* **False.** Rooted-tree counts grow only
   exponentially after symmetry factors; we measured a finite radius
   ($\approx 0.302$ on the test ray) matching the exact branch locus.
2. *"Borel resummation is impossible; transcendental cuts."* **Misleading.**
   The inverse is an explicit algebraic function of degree 3; no resummation
   is needed at all. The failure is *globality*, not summability.
3. *"This explains why $D = 4$ QFT is harder than $D \le 3$."* **Category
   error.** The $D=4$ difficulty is UV renormalization / triviality; the
   present phenomenon is tree-level, global, and dimension-0. See
   `docs/QFT_IMPLICATIONS.md`.

## 2026-07-20 — QFT implications document

`docs/QFT_IMPLICATIONS.md` (456 lines; all 16 references verified against the
literature). Key theses: the tree expansion converges to one *local* branch of
an explicit degree-3 algebraic function; the obstruction is Jelonek
non-properness, not divergence; constant-Jacobian field redefinitions require
global-injectivity checks before use as exact changes of variables
(equivalence-theorem caveat); "vacua at infinity" are a non-perturbative
mechanism genuinely distinct from instantons; no bearing on UV
renormalization, $\phi^4_{2,3}$ Borel summability, $\phi^4_4$ triviality, or
Yang–Mills.

## 2026-07-20 — Monodromy of the three sheets

`scripts/monodromy.py`, `docs/MONODROMY.md` (high-precision predictor–corrector
continuation of all 3 sheets, mp.dps = 30, residuals $\le 10^{-22}$;
cross-checked on a second generic line)

- **The geometric monodromy group is the full $S_3$.** On a generic line in
  source space, each loop around a root of $p(t)$ yields a *transposition*
  (the two sheets with $x \sim \pm\sqrt{-q/p}$ swap while the finite sheet
  stays fixed); the four loops gave $(1\,2), (1\,3), (2\,3), (2\,3)$ —
  three distinct transpositions, generating all of $S_3$.
- **Exact identity** (verified symbolically): with $D_0 = 27ac^2 - 9bc + 8$
  (the denominator of the $y,z$ fiber parametrization),
  $$4q^3 + 27 p r^2 = 4\,D_0^2, \qquad\text{so}\qquad
    \operatorname{disc}_X = -4\,D_0^2\,p .$$
  The $x$-collision locus is exactly $\{D_0 = 0\}$, and the discriminant is
  manifestly a non-square (multiplicity-1 factor $p$), so the Galois group of
  the cubic over $\mathbb{Q}(a,b,c)$ is $S_3$ — realized geometrically by the
  numerics.
- **Trivial monodromy confirmed** around $\{D_0 = 0\}$ (identity permutation,
  as forced by $\det DF = -2$: no ramification, only shared $x$-coordinates)
  and around a big loop enclosing all singularities on the line (the four
  transpositions cancel in path order).

Physics reading: dialing the external source around the escape locus
$\{p=0\}$ exchanges the two "vacua at infinity" with one another or with the
perturbative sheet — a source-space analogue of vacuum interchange, mediated
entirely through the boundary of field space. Reproduce:
`.venv/bin/python scripts/monodromy.py` (~80 s; `monodromy.py 7` for the check line).

## 2026-07-20 — Field-redefinition "measure anomaly" (non-properness defect)

`scripts/measure_anomaly.py` (~6 s); write-up in `docs/QFT_IMPLICATIONS.md` §6.

- $F$ is also a *real* local diffeomorphism $\mathbb{R}^3\to\mathbb{R}^3$
  (Jacobian $-2$, three real preimages). **Chamber rule (exact):** the number
  of real preimages is
  $$N(J) = 3 \iff p(a,b,c) < 0, \qquad N(J) = 1 \iff p > 0,$$
  since the monic discriminant of the $x$-cubic is $-4D_0^2/p^3$. Verified
  against direct preimage computation at 300 random targets (0 mismatches).
- Change of variables for non-injective local diffeos carries the factor
  $N(J)$: the naive equivalence-theorem substitution $\phi' = F(\phi)$
  ("constant Jacobian, harmless") is wrong by a **step function**. Measured
  anomaly factor $A(\sigma) = \langle N\rangle$ over Gaussian source
  ensembles ($10^7$ samples): $A(10) = 1.205$, $A(1) = 1.696$,
  $A(0.1) = 2.006$, $A(0.01) = 2.001$, $A(0.001) = 2.000$.
- **$A \to 2$ as $\sigma \to 0$**: since $J=0$ lies on the wall $\{p=0\}$
  (linear part of $p$ is $16a$), the defect is $O(1)$ arbitrarily close to
  the perturbative vacuum — invisible to every perturbative order.
- Upstream practices affected: nonlinear field redefinitions certified only
  by Jacobians; Gribov copies *without* a Gribov horizon (FP-type
  determinant $\equiv -2$, copies anyway — the horizon sits at infinity);
  Nicolai maps (global invertibility is an unchecked hypothesis).
- **Sum over all sheets restores single-valuedness (CHY analogy):** the
  symmetric functions of the three $x$-roots are rational,
  $e_1 = 0,\ e_2 = q/p,\ e_3 = -r/p$, with poles exactly on $\{p=0\}$.
  Monodromy-invariant "sum over all vacua" observables are single-valued;
  the multivaluedness is an artifact of selecting one branch.

## 2026-07-20 — $\mathbb{C}^*$-equivariance: the orbifold mechanism

(verified symbolically; one-liner in the session log)

- The map is **equivariant under a weighted scaling action**: with source
  weights $(1,-1,-2)$ and target weights $(-2,-1,1)$,
  $$F(\lambda x,\ \lambda^{-1}y,\ \lambda^{-2}z)
    = (\lambda^{-2}F_1,\ \lambda^{-1}F_2,\ \lambda F_3).$$
  (Consistent with all interactions being organized by the invariant
  $w = xy$ and the building blocks $g = 1+xy$, $u$.)
- The 2:1 curve $\varphi(s) = (s, -\tfrac{3}{2s}, \tfrac{13}{2s^2})$ is
  precisely the $\mathbb{C}^*$-**orbit** of $(1,-\tfrac32,\tfrac{13}{2})$;
  its image $(-\tfrac{1}{4s^2},0,0)$ lies on the target stratum where only
  the weight-$(-2)$ coordinate is nonzero, so $\lambda \mapsto \lambda^{-2}$
  is 2-to-1 there: $\pm\lambda$ are indistinguishable.
- **Physics reading:** the non-injectivity is an *orbifold phenomenon* — a
  residual $\mathbb{Z}_2 \subset \mathbb{C}^*$ that acts freely on field
  space but trivially on the sources. This suggests a systematic search
  principle for new counterexamples: enumerate weight systems (graded scalar
  models), impose the Keller condition on equivariant maps (a
  finite-dimensional ansatz per weight system), and look for orbits whose
  images land in even-weight strata.

## 2026-07-20 — Trace/pushforward structure and the amplitudes dictionary

`scripts/trace_pushforward.py` (~1 s); full discussion in
`docs/AMPLITUDES_CONNECTION.md`.

- **Trace-map rationality (CHY/pushforward mechanism, exact):** the sum over
  all three sheets of any polynomial observable is rational with poles
  confined to $\{p=0\}$; verified for power sums $S_2..S_6$ (denominators
  $p^{\lceil(k-1)/3\rceil}$). Single sheets: multivalued, $S_3$ monodromy.
  The tree expansion is a "one-solution CHY formula" — its failure is the
  failure to sum over all solutions.
- **Boundary factorization at the wall:** near $p \to 0$ the fiber splits
  into an escaping pair $x_\pm \sim \pm\sqrt{-q/p}$ and a finite sheet
  $x_3 = -r/q + p\,r^3/q^4 + O(p^2)$; the polar part factorizes Vieta-exactly,
  $e_3 = (x_+x_-)\cdot x_3 = [q/p + O(1)]\cdot[-r/q + O(p)]$ — the analogue
  of amplitude factorization on a physical pole.
- **Landau-discriminant dictionary (Mizera–Telen):** eliminating fields from
  the field equations and stratifying the discriminant is methodologically
  identical to computing Landau discriminants of Feynman integrals; the
  model has *no* first-type (pinch) singularities (étale) and *only*
  second-type ones (escape to infinity) — a solvable laboratory for the
  least-understood class of amplitude singularities.
- Effective-action reading: the source–field relation $J \mapsto \bar\phi(J)$
  underlying $\Gamma$ is exactly the object that is globally 3:1 here; the
  1PI construction silently selects a branch with no finite-distance flag.

## 2026-07-20 — Search strategies and defect phase diagram

`docs/SEARCH_STRATEGIES.md`: three physics-derived strategies for finding new
counterexamples (graded/orbifold ansatz; bootstrap on the inverse data with
the forced consistency conditions — perfect-square collision factor,
integrality off the escape directions, Jelonek uniruledness; Witten-index
prefilter at infinity), and the "phase diagram of non-properness defects"
reading of what new examples would mean (gradient example = biggest prize;
$\mathbb{Z}_3$ vs $S_3$ monodromy classes; $d \ge 4$ chamber combinatorics).
Also records the normalization insight: "vacuum on the wall" is a translation
artifact; the invariant content is the wall + chamber structure.

## 2026-07-20 — Construction mechanism and rigidity (search for new examples)

`scripts/search_counterexamples.py` (~40 s); full report in
`docs/NEW_COUNTEREXAMPLES.md`. (Sub-agent exploration consolidated and
re-verified independently after two session timeouts.)

- **Normal form:** in the invariants $w = xy$, $v = x^2z$,
  $F = (P/x^2,\ Q/x,\ xR)$ with explicit $P, Q, R$; realizability =
  $x$-degree conditions on the monomials of $P, Q$.
- **2D reduction of the Keller condition (exact):**
  $\det DF = \kappa \iff \partial(PR^2, QR)/\partial(w,v) = \kappa R^2$.
  The 3D search inside the symmetry class is a two-variable problem.
- **Non-injectivity criterion:** a common zero of $(Q, R)$ with $P \neq 0$
  gives a 2:1 orbit (residual $\mathbb{Z}_2$); reproduces the triple point
  from $(w_0,v_0) = (-3/2, 13/2)$, $P = -1/4$.
- **Rigidity evidence:** in the box $\deg_x P\in[2,7]$, $Q\in[1,6]$,
  $R\in[0,4]$: linearized-Keller kernel dim 15; in-box gauge rank 9 (all
  inside the kernel); all 6 non-gauge directions fail nonlinear continuation
  (obstructed or jump to distant gauge images). **Zero nearby families** —
  the Alpöge–Mathew map is locally rigid modulo gauge in this class.
  (Sub-agent cross-check in a different box: kernel dim = gauge dim = 11.)
- **No new counterexample found**; dead ends recorded (a $\kappa = 0$
  branch collapse). Next: other weight systems (potential $\mathbb{Z}_3$
  monodromy), larger boxes, the bootstrap route, exact Groebner
  certification of rigidity.

## 2026-07-20 — Refactor: `jcqft` package + `scripts/`; AQFT implications fleshed out

Code reorganization (behavior-preserving; all seven scripts re-run and
outputs verified identical — radius ≈ 0.302, monodromy $S_3$, $A(\sigma)\to 2$,
rigidity result unchanged):

- Shared library extracted into the installable package `jcqft/`
  (`core.py` map/eliminant/chamber rule, `truncated.py` tree-graph inverse,
  `fibers.py` Groebner fiber parametrization + numeric helpers,
  `reduction.py` equivariant normal form + 2D Keller reduction). The
  package asserts its own consistency on import (eliminant matches the
  Groebner basis, $A = 2D_0$, reduced Keller identity holds).
- All analyses moved to `scripts/`; duplicated Groebner/parametrization code
  in `branch_locus`/`measure_anomaly` now comes from `jcqft.fibers`.
- `pyproject.toml` added; install with `uv pip install -e .`
  (replaces `requirements.txt`).

Documentation:

- `docs/QFT_IMPLICATIONS.md` §4.3 (new): implications for algebraic QFT,
  layer by layer — the exact 0D statement that $F^*$ is a monomorphism but
  not an automorphism of observable algebras (degree-3 extension, Galois
  group $S_3$), vindicating the algebra-first viewpoint [HK64]; what the
  example does and does not say about pAQFT's formal-series constructions
  [FR12, FR16] (calibration of the deferral, no inconsistency; "formal to
  convergent" is not the missing step — properness is); a well-posed 0D
  exercise for the non-perturbative Buchholz–Fredenhagen approach [BF20];
  and an explicit list of what the 0D model cannot address.
- `docs/QFT_IMPLICATIONS.md` §4.4 (new): **claims ledger** — every verified
  claim, the script that verifies it, and the nearest overstatement it does
  *not* license.

## 2026-07-21 — Chamber geometry: cuspidal cubic wall, exact non-surjectivity locus, positive-geometry verdict

`scripts/positive_geometry.py` (exact symbolic, ~2 s); figure
`scripts/plot_chamber.py`; full write-up in `docs/POSITIVE_GEOMETRY.md`.
This answers question (1) of `docs/AMPLITUDES_CONNECTION.md` §2.4.

- **C\*-reduction of the wall:** $p$ is quasi-homogeneous (weight $-2$)
  under the source weights $(-2,-1,1)$. In the invariants $u = ac^2$,
  $w = bc$: $c^2 p = P_2(u,w) = 27u^2 + 16u - 18uw + w^3 - w^2$ — the
  quartic wall in $\mathbb{C}^3$ is a **plane cubic**.
- **The reduced wall is a cuspidal cubic** (unique $A_2$ singularity at
  $(4/27, 4/3)$, hence rational), with exact parametrization
  $u = 4/27 - 3(m-3)^2/m^3$, $w = 4/3 - 3(m-3)^2/m^2$. Its **cuspidal
  tangent is the $\{D_0=0\}$ line** — the geometric origin of the identity
  $4q^3 + 27pr^2 = 4D_0^2$.
- **Exact non-surjectivity locus:** $F(\mathbb{C}^3)$ misses precisely the
  C\*-orbit $\{ac^2 = 4/27,\ bc = 4/3\}$ (the cusp): there $p = q = 0$,
  $r \neq 0$, the fiber is empty, **all three sheets are at infinity** —
  sources with no classical solution at all. Generic wall points lose 2
  sheets; the cusp loses 3.
- **The $N=3$ chamber is NOT a positive geometry:** the residue of the
  candidate canonical form on the normalized wall is
  $-(\kappa/3)\,dm/(m-3)^2$ — a residueless **double pole at the cusp**,
  not a log form; admitting an adjoint line (KPR+25 polypol framework,
  unbounded chamber with flex at infinity) leads to a linear system with
  only the zero solution. Failure mechanism: **vertex collision**
  (node → cusp degeneration of the would-be interval form), *not* genus —
  nodal-wall chambers would pass. The negative branch of the
  falsifiability criteria in `docs/AMPLITUDES_CONNECTION.md` §2.5 is
  realized in exactly the predicted form.
- **Pushforward dichotomy (exact):** $F_*(d^3\phi) = -\tfrac32\,d^3J$
  (holomorphic — constant, wall invisible) vs $\tfrac{N(J)}{2}\,d^3J$
  (real — jumps $\tfrac12 \leftrightarrow \tfrac32$). The measured anomaly
  $A(\sigma) \to 2$ is exactly the gap between the two pushforwards.
- Cross-check: the ray $J = t(1,2,3)$ maps to $(9t^3, 6t^2)$ in the
  invariant plane; nearest wall crossing $|t^*| = 0.302028$ = the verified
  tree-series radius.

## 2026-07-21 — Lagrangian stated; open questions consolidated and ranked

- **`docs/PROBLEM.md`, new section "The model, stated as a Lagrangian":**
  the honest statement that no single-field action exists ($DF$ non-symmetric
  already at the origin: $\partial F_1/\partial z(0) = 1 \neq 2 =
  \partial F_3/\partial x(0)$), the correct first-order (conjugate-field)
  action $S = \bar\varphi\cdot(F(\varphi) - J)$ à la Abdesselam [Abd03] with
  trivial ghost sector (constant $\det DF = -2$), the propagator
  $\langle\varphi_i\bar\varphi_j\rangle = (L^{-1})_{ij}$, all **13
  interaction vertices** listed explicitly (cubic to octic, one
  $\bar\varphi$-leg each), the $\mathbb{C}^*$-grading extended to
  $\bar\varphi$ (weights $(2,1,-1)$; every term of $S$ has weight 0), and
  the identification of the tree expansion as the one-point function.
  All expansions verified symbolically against `jcqft/core.py`.
- **`docs/OPEN_QUESTIONS.md` (new):** every open question from all documents,
  merged, de-duplicated, and ranked into tiers by effort and likelihood of
  fast progress (Tier A: cusp-orbit trace observables, pushforward closure,
  infinity prefilter, larger boxes; Tier B: Buchholz–Fredenhagen 0D
  caricature, wall-complement local system, other weight systems, exact
  rigidity certification, damped integrals, certified monodromy; Tier C:
  bootstrap realization, chamber classification, 3D cone, $d=4$, BCW
  reduction, gradient counterexample, $n=2$, functional-integral vacua at
  infinity). One stale claims-ledger row corrected
  (`QFT_IMPLICATIONS.md`: chamber positive-geometry question now marked
  resolved for this map).

## 2026-07-21 — A1 and A2 resolved: trace asymptotics at wall/cusp; pushforward of general forms

**A1** (`scripts/cusp_traces.py`, ~2 s; write-up `docs/POSITIVE_GEOMETRY.md`
§6). The original Q3 parenthetical was wrong ($e_2, e_3$ do NOT stay finite
at the cusp). Exact results:

- $S_1 \equiv 0$ (the only finite trace); $S_2 = -2q/p$, $S_3 = -3r/p$, ….
- Generic wall: escaping pair $\pm\sqrt{-q/p}$, pole order
  $\lfloor k/2 \rfloor$ for $S_k$ (odd-power pair cancellation).
- Cusp, generic approach: three sheets escape as asymptotic cube roots of
  unity; divergence rate of $S_k$ drops below naive $2k/3$ **iff**
  $3 \nmid k$ ($\omega$-cancellation); exact rate
  $\max\{a+2b : 2a+3b=k\}$, verified $k \le 9$.
- Cusp, tangent (D0-line) approach: exactly solvable for all $\varepsilon$:
  $27\varepsilon^3X^3\!-\!9\varepsilon X\!-\!2 =
  (3\varepsilon X\!-\!2)(3\varepsilon X\!+\!1)^2$; collided-$x$ pair
  persists into the cusp, fiber stays étale (distinct $(y,z)$).

**A2** (`scripts/pushforward_forms.py`, ~8 s; banner in
`docs/AMPLITUDES_CONNECTION.md` §2.4). Question 2 answered affirmatively:
$F_*(g\,d^3\phi) = -\tfrac12 T[g]\,d^3J$ is rational with poles ONLY on the
wall for every polynomial $g$; all $D_0$-collision singularities of the
per-sheet parametrization cancel in the sheet sum (verified on 15
observables; conceptually forced by étaleness off the wall). Pole-order law
$\mathrm{ord}_p T[x^k] = \lfloor k/2 \rfloor$; wall poles sourced only by
$x$-powers $\ge 2$ (only $x$ escapes). Exact extremes:
$F_*(d^3\phi) = -\tfrac32\,d^3J$ (wall invisible) and
$F_*(x\,d^3\phi) = 0$ identically.

## 2026-07-21 — B2 resolved: the wall complement is the braid-group classifying space

`scripts/wall_braid.py` (~6 s; symbolic identities exact, monodromy
numerical with three cross-checks); full write-up in
`docs/WALL_COMPLEMENT.md`. (Sub-agent work, independently re-run and
integrated.)

- **Invariant eliminant (exact):** $\xi = X/c$ turns the eliminant into
  $P_2(u,w)\,\xi^3 + (4-3w)\,\xi - 2$ — the whole covering is a cubic
  family over the invariant plane.
- **The wall IS the $A_2$ discriminant (exact, affine):**
  $(Q,R) = (w - \tfrac43,\ 2u - \tfrac{2w}{3} + \tfrac{16}{27})$ satisfies
  $4Q^3 + 27R^2 = 4P_2$ — an affine isomorphism of pairs sending cusp to
  origin and the $D_0$-line to the cuspidal tangent $\{R=0\}$. The
  universal coordinates are $q/9$ and $D_0/27$ themselves.
- **Consequence:** $\pi_1(\mathbb{C}^2 \setminus \{P_2=0\}) = B_3$ and the
  complement is a $K(B_3,1)$ (Arnold, Brieskorn, Deligne); the verified
  $S_3$ sheet monodromy is the canonical $B_3 \twoheadrightarrow S_3$.
- **Cusp-loop monodromy measured: a Coxeter element** — 3-cycle of order
  $h = 3$, image of $\sigma_1\sigma_2$ (not of the half-twist); the full
  twist $(\sigma_1\sigma_2)^3$ acts trivially (verified). A loop around the
  empty-fiber orbit cyclically rotates the three vacua.
- **Local system:** rank-3 sheet system $=$ trivial $\oplus$ standard;
  trace rationality $=$ trivial summand; ALL multivaluedness lives in the
  reflection representation of $W(A_2)$. $\chi(\text{complement}) = 0$.
- **Proposal (marked as such):** the amplituhedron-analogue is
  (wall complement, standard local system) $=$ the $K(B_3,1)$ with its
  reflection local system; twisted periods replace canonical-form
  integrals. Concrete follow-up: twisted $H^1$ + intersection pairing
  (`docs/WALL_COMPLEMENT.md` §6).

## 2026-07-21 — Twisted cohomology of the wall complement (WALL_COMPLEMENT §6 Q3)

`scripts/twisted_cohomology.py` (~3 s; exact sympy — Fox calculus on the
1-relator presentation of $B_3$, cross-checked against the Wang sequence of
the global Milnor fibration $P_2 : M \to \mathbb{C}^*$; the two methods
agree everywhere, with identical jump polynomials). Full write-up in
`docs/TWISTED_PERIODS.md`. (Sub-agent work, independently re-run and
integrated.)

- **Complete exact dimension table** for $H^*(B_3;\rho) = H^*(M;\rho)$:
  Kummer twist $P_2^s$ jumps ONLY at $s \in \mathbb{Z}$ ($(1,1,0)$,
  untwisted) and $s \in \pm\tfrac16 + \mathbb{Z}$ ($(0,1,1)$) — the roots
  of the trefoil Alexander polynomial $\Delta(t) = t^2 - t + 1$, which
  appears verbatim as the Fox differential. Standard (reflection) local
  system: $(0,1,1)$ — the proposed amplituhedron-analogue pair carries
  exactly one twisted class per degree. Reduced Burau: jump locus
  $\{t^3 = 1\}$; Burau at $t = +1$ IS the reflection rep (at $t = -1$ it is
  the $SL(2,\mathbb{Z})$ rep — a corrected expectation, flagged honestly).
- **The failed canonical form, explained cohomologically:** at integer
  twist the form is globally twisted-exact,
  $dU{\wedge}dW/f = d(-(2U\,dW - 3W\,dU)/f)$ — every period vanishes;
  $h^2 = 0$. The residueless double pole of `docs/POSITIVE_GEOMETRY.md`
  is the local face of this global exactness. The classes revive exactly
  at $s \equiv \pm\tfrac16$: the singularity spectrum of the cusp.
- **Generic twists carry NOTHING:** $(0,0,0)$, since $|\chi(M)| = 0$ —
  this model has zero generic-twist master integrals; all content is
  resonant. (Contrast with generic hyperplane-arrangement amplitudes,
  where $|\chi|$ counts master integrals.)

## 2026-07-21 — A3 + A4 resolved: infinity prefilter; rigidity survives larger boxes

`scripts/witten_prefilter.py` (~6 s) and `scripts/rigidity_boxes.py`
(~50 s); reusable module `jcqft/prefilter.py`; write-up in
`docs/RIGIDITY_AND_PREFILTER.md`. (Sub-agent work, independently re-run
and integrated.)

- **A3, infinity prefilter:** if the leading forms of a polynomial map
  have only the trivial common zero, the map is provably proper — and a
  proper Keller map is injective, so no counterexample is proper. The
  filter rejects such maps in ~1–30 ms (one grevlex Gröbner basis of the
  leading forms), before any symbolic Keller work; 199/200 random cubic
  maps rejected. The Alpöge–Mathew map survives with witness $[1:0:0]$ at
  infinity — exactly its escape direction. Honest limitation: the test is
  necessary-not-sufficient; nonlinear polynomial automorphisms always
  survive (Bézout argument), so it screens, it does not decide. Plain and
  weighted variants.
- **A4, larger degree boxes:** the first-order/continuation rigidity
  analysis re-run in two strictly larger ansatz boxes in the
  $(1,-1,-2)$ class (57×43 → 65×55 → 74×69 equations × unknowns). In
  every box the linearized-Keller kernel = gauge tangents ⊕ obstructed
  directions; **every non-gauge direction (6, 9, 15 per box) is
  obstructed at second order**, and nonlinear continuation confirms none
  integrates to a nearby family. Verdict: still rigid modulo gauge —
  materially stronger evidence that the counterexample is isolated in
  this weight class. (Numerical; exact Gröbner certification remains
  OPEN_QUESTIONS B4.)

## 2026-07-21 — B5 resolved: the damped partition function, exactly

`scripts/damped_partition.py` (~30 s; exact claims asserted, numerics with
convergence tables; figure `outputs/damped_partition.png`); write-up in
`docs/DAMPED_PARTITION.md`. (Sub-agent work, independently re-run and
integrated.)

- **Exact closed form, finiteness unconditional:** via the verified
  pushforward $F_*(d^3\phi) = \tfrac{N(J)}{2} d^3J$,
  $Z_\hbar(J) = (2\pi\hbar)^{3/2}\bigl(\tfrac12 +
  \mathbb{P}[p(J+\sqrt{\hbar}\,\xi)<0]\bigr)$ with $\xi$ standard normal —
  finite and uniformly bounded for EVERY $J$, Jelonek set and cusp
  included. The escaping tube is unbounded but its cross-section decays
  like $|x|^{-3}$; no divergence detects the wall. (Corrects the guessed
  mechanism: the signature is in the prefactor, not in finiteness.)
- **Semiclassical prefactor = chamber function:** $Z/(2\pi\hbar)^{3/2}
  \to N(J)/2$ per chamber at a pure Gaussian rate (verified to $10^{-4}$);
  on the wall it tends to the two-sided mean $1$ — at the perturbative
  vacuum $J=0$ the damped partition function counts TWICE the perturbative
  saddle (it sees the sheets at infinity); at the empty-fiber cusp it
  tends to $\tfrac12$ with an anomalous $\hbar^{1/4}$ correction whose
  amplitude is computed in closed form from the $A_2$ normal form.
- **Uniformity boundary:** the $\hbar$-expansion around one vacuum is
  uniform only for $\hbar \ll \mathrm{dist}(J, \text{wall})^2$; measured
  $\gamma_{\text{wall}} = 2.001$ and $\gamma_{\text{cusp}} = 2.9996$ with
  closed-form constants matched to $<1\%$. $\gamma_{\text{cusp}} = 3 =
  2\times\tfrac32$ is the $A_2$ horn-width exponent: the boundary is set
  by target-space chamber geometry, NOT by the field-space escape rates.
- **QFT synthesis:** a constructively finite partition function whose
  semiclassical normalization is a step function of the source, jumping
  on the non-properness set, while every local perturbative datum is
  chamber-independent (det $DF$ constant). "Constant Jacobian implies
  trivial semiclassics" holds per vacuum, fails globally. 0D statement
  only; no claim about $D \ge 1$.

## 2026-07-22 — Missing observables made explicit: normal form, generators, non-integrality

`scripts/missing_observables.py` (~15 s, exact; all identities asserted);
write-up `docs/MISSING_OBSERVABLES.md`; `docs/QFT_IMPLICATIONS.md` §4.3(a)
and the claims ledger updated.

- **Field-level structure theorem:** $\mathbb{C}(x,y,z) = K(x)$ with
  $K = \operatorname{Frac}(\operatorname{im}F^*)$; the eliminant cubic is
  irreducible over $K$, so $\{1, x, x^2\}$ is a basis. Every observable
  has a unique normal form $c_0(F) + c_1(F)x + c_2(F)x^2$ (one polynomial
  division); membership in $\operatorname{im}F^*$ iff $c_1 = c_2 = 0$.
  The missing observables are the module generated by $\{x, x^2\}$ — the
  sheet separators. The separator coefficients of $y$ and $z$ carry the
  non-properness polynomial $p$ verbatim.
- **No finite-module version, provably:** $x$ is NOT integral over
  $\operatorname{im}F^*$ — exact escape certificate
  $\phi(T) = (T, y_0, (2T - 3T^2y_0 - c_3)/T^3)$ with $F(\phi(T))$
  converging to a finite wall point; incidentally a rational
  parametrization of the Jelonek set by $(y_0, c_3)$. Non-properness =
  non-finiteness of the ring extension.
- **Cheap membership certificates:** the exact fiber over $F(1,2,3)$ has
  pairwise-distinct $x$-, $y$-, $z$-coordinates, so none of the source
  coordinates is a function of the redefined field (fiber-separation
  beats the too-slow tag-variable Groebner test).
- **Fiber exhaustiveness spot-verified** on generic / image / wall /
  cusp targets (3, 3, 1, 0 points; wall eliminant degenerates to
  $qX + r$).

## 2026-07-22 — Exact reduced Keller identity for ALL weight systems (1,-1,-m)

`scripts/reduction_113.py` (~4 s, exact); reusable module
`jcqft/reduction_w.py` (generalizes `jcqft/reduction.py`, m = 2).
OPEN_QUESTIONS B3, step 1: the reduction that any (1,-1,-3) search needs.

- **Generic-function proof** (undetermined functions $P, Q, R$ of the
  invariants $w = xy$, $v = x^m z$; sympy, m = 1..4): for
  $F = (P/x^m,\ Q/x,\ xR)$,
  $$\det DF = -mP\,J_2(Q,R) + Q\,J_2(P,R) + R\,J_2(P,Q)$$
  — a function of $(w,v)$ alone, with the compact form
  $$R^m \det DF = J_2(P R^m,\ Q R).$$
  **Reduced Keller identity: $\det DF = \kappa \iff J_2(PR^m, QR) =
  \kappa R^m$** ($m=2$ recovers `jcqft.reduction` exactly, verified on
  Alpöge–Mathew).
- **Polynomiality boxes:** $P$ needs $j + mk \ge m$, $Q$ needs
  $j + mk \ge 1$ per monomial $w^jv^k$, $R$ unconstrained; invertible
  $DF(0)$ forces $v \in P$, $w \in Q$, $1 \in R$. Verified exhaustively
  in a box (m = 3) including the rejected complements.
- Random-ansatz consistency check at m = 3: 3D symbolic determinant ==
  reduced 2D form.
