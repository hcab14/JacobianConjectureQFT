# Progress log

All results below are produced by scripts in this repository (exact symbolic
computation with SymPy unless stated otherwise) and can be reproduced with
the commands given. Environment: `.venv` (Python 3.12, `requirements.txt`).

## 2026-07-20 — Verification of the counterexample

`verify_counterexample.py`

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

`branch_locus.py`

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

`tree_expansion.py` (sparse truncated-ring arithmetic; Picard iteration =
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

`monodromy.py`, `docs/MONODROMY.md` (high-precision predictor–corrector
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
`.venv/bin/python monodromy.py` (~80 s; `monodromy.py 7` for the check line).

## 2026-07-20 — Field-redefinition "measure anomaly" (non-properness defect)

`measure_anomaly.py` (~6 s); write-up in `docs/QFT_IMPLICATIONS.md` §6.

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

`trace_pushforward.py` (~1 s); full discussion in
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

`search_counterexamples.py` (~40 s); full report in
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
