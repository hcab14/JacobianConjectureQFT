# The Infinity Prefilter, and Rigidity in Larger Degree Boxes

*(2026-07-21. Resolves items A3 and A4 of `docs/OPEN_QUESTIONS.md`. All
claims verified by asserts in `scripts/witten_prefilter.py` (~6 s) and
`scripts/rigidity_boxes.py` (~50 s); the reusable prefilter lives in
`jcqft.prefilter`. Strategy context: `docs/SEARCH_STRATEGIES.md` §1.3;
prior rigidity evidence: `docs/NEW_COUNTEREXAMPLES.md` §3.)*

**Summary of the two results.**

1. **The Witten-index / infinity prefilter is implemented and validated**
   (`jcqft.prefilter.infinity_prefilter`). It rejects, in milliseconds and
   before any symbolic Keller work, every map whose leading forms have no
   common zero at the hyperplane at infinity — such maps are *provably
   proper*, hence cannot be counterexamples. On a benchmark of 200 random
   cubic maps it rejects 199 (~27 ms/map), so it can screen thousands of
   ansatz candidates per minute. The Alpöge–Mathew map survives, with
   witness direction $[1:0:0]$ at infinity — exactly the escape direction
   of its 2:1 orbits and the flex-at-infinity of its wall cubic. The test
   is necessary-not-sufficient: nonlinear polynomial *automorphisms*
   provably always survive (a Bézout argument, see §1.3), and do.
2. **The Alpöge–Mathew solution stays rigid in strictly larger degree
   boxes.** The first-order/continuation analysis of
   `scripts/search_counterexamples.py` was re-run in two strictly larger
   ansatz boxes (max weighted degree raised by 1 and by 2; up to 69
   coefficient unknowns, 74 equations). In every box: kernel of the
   linearized Keller condition splits as (gauge tangents) $\oplus$
   (obstructed directions); **every non-gauge direction is already
   obstructed at second order**, and nonlinear continuation confirms none
   integrates to a nearby family. **No new deformation direction appeared;
   the expected outcome — "still rigid, stronger evidence" — is the one
   realized.**

---

## 1. The infinity prefilter (A3)

### 1.1 Exact statement

Let $F = (F_1, \dots, F_n) : \mathbb{C}^n \to \mathbb{C}^n$ be polynomial
and let $L_i$ be the leading form (top-degree homogeneous part) of $F_i$.

> **Properness criterion.** If $V(L_1, \dots, L_n) = \{0\}$, then $F$ is
> proper: $|F(\phi)| \to \infty$ whenever $|\phi| \to \infty$.

*Proof sketch:* by homogeneity and compactness of the unit sphere,
$\max_i |L_i(u)| \ge c > 0$ on $|u| = 1$, so along any ray
$\max_i |F_i(r u)| \ge c\, r^{\deg F_i} - O(r^{\deg F_i - 1})$. $\square$

A proper Keller map is a finite étale cover of the simply connected
$\mathbb{C}^n$, hence injective — no counterexample is proper (its Jelonek
non-properness set must be nonempty; cf. the escape locus $\{p=0\}$ of the
Alpöge–Mathew map). Contrapositive, as a filter:

- `infinity_prefilter(F, vars) == False` — the leading forms have only the
  trivial common zero: $F$ is **provably proper, reject** before any
  symbolic work.
- `infinity_prefilter(F, vars) == True` — the leading forms degenerate
  somewhere at infinity: the candidate **survives** to the expensive
  Keller/injectivity analysis.

Implementation: one grevlex Gröbner basis of the leading forms; since they
are homogeneous, $V = \{0\}$ iff the ideal is zero-dimensional (checked via
`sympy`'s zero-dimensionality test on the basis). `infinity_witness`
additionally produces the common projective zero for reporting.

### 1.2 Weighted variant

For a $\mathbb{C}^*$-equivariant search class one passes the source weight
vector (for the Alpöge–Mathew system: $(1,-1,-2)$ on $(x,y,z)$, the
convention of `jcqft.core`/`jcqft.reduction`). The relevant escape curves
are weighted orbits $\phi_j(\lambda) = \lambda^{w_j} c_j$; along such a
curve $F_i = \lambda^{d_i} L^w_i(c) + \dots$ with $L^w_i$ the $w$-leading
part and $d_i$ its weighted degree. The test asks whether the
positive-weighted-degree leading parts share a zero $c$ with some nonzero
positive-weight coordinate (Rabinowitsch/radical-membership per coordinate,
since with mixed-sign weights $V$ need not be a cone); both scalings
$\lambda \to \infty$ and $\lambda \to 0$ are checked. **Scope caveat**: a
weighted `False` only rules out first-order escape along $w$-orbits; the
unweighted call (the $w = (1,\dots,1)$ case) is the one that proves
properness outright.

### 1.3 Validation (all asserted in `scripts/witten_prefilter.py`)

| map | verdict | ms | comment |
|---|---|---|---|
| Alpöge–Mathew, unweighted | **survives** | 7.6 | witness $[1:0:0]$ = escape direction (below) |
| Alpöge–Mathew, weights $(1,-1,-2)$ | **survives** | 5.0 | witness on the escape locus $\{R = 0\}$ |
| linear $L = DF(0) = (z, y, 2x)$ | rejected | 0.9 | correct (proper) |
| linear $(x+y,\, y+z,\, z+x)$ | rejected | 1.2 | correct (proper) |
| elementary $(x + y^2,\, y + z^2,\, z)$ | survives | 0.9 | **false positive** (automorphism) |
| elementary $(x + (y+z)^3,\, y,\, z)$ | survives | 1.7 | **false positive** (automorphism) |
| tame composition #1 (elem∘elem∘linear) | survives | 4.7 | **false positive** (automorphism) |
| tame composition #2 (linear∘elem∘shear) | survives | 2.6 | **false positive** (automorphism) |
| proper power map $(x^3{+}y,\, y^3{+}z,\, z^3{+}x)$ | rejected | 1.2 | correct (proper, non-injective, non-Keller) |
| 200 random dense cubic maps | 199 rejected, 1 survives | 27/map | the Bézout-generic bulk is filtered out |

Two structural facts, both verified:

- **The witness matches the known escape story.** Unweighted leading forms
  of the counterexample are $(x^3y^3z,\ 3x^3y^2z,\ -x^3z)$, all vanishing
  along the $x$-axis: the witness $[1:0:0]$ is precisely the direction of
  the 2:1 escape orbits $x \mapsto (x, w_0/x, v_0/x^2)$, $x \to \infty$,
  and the unique point at infinity (the flex) of the wall cubic
  (`docs/POSITIVE_GEOMETRY.md` §1). The weighted witness satisfies
  $R(w, v) = 2 - 3w - v = 0$ in the invariants $w = xy$, $v = x^2 z$ — the
  normal-form escape condition of `jcqft.reduction`.
- **The false positives are forced, not a bug.** Any injective polynomial
  map with some $\deg F_i > 1$ *must* survive every leading-form test: if
  its leading forms had only the trivial common zero, Bézout would give
  generic fiber cardinality $\prod_i \deg F_i > 1$, contradicting
  injectivity. So no prefilter of this type can separate nonlinear
  automorphisms from counterexamples; its job — done, per the benchmark —
  is to discard the generic bulk of a search space (99.5% of random dense
  cubics here) at negligible cost.

**Honest limitations.** (i) Necessary-not-sufficient: surviving means only
"not provably proper by leading forms"; a survivor can still be proper (the
degeneration at infinity might not be realized by actual escape) or an
automorphism. (ii) The weighted variant is a first-order-in-$\lambda$
screen along torus orbits, not a properness proof. (iii) Timings are for
degree $\le$ 7-ish inputs in 3 variables; Gröbner cost grows with degree
and dimension, though leading-form systems are small.

## 2. Rigidity in larger degree boxes (A4)

### 2.1 Setup

Same machinery as `scripts/search_counterexamples.py` (all in the
2D-reduced picture of `jcqft.reduction`: $F = (P/x^2,\, Q/x,\, xR)$ in
$w = xy$, $v = x^2z$; Keller $\Leftrightarrow J_2(PR^2, QR) = \kappa R^2$).
Ansatz boxes prescribe the allowed monomials $w^i v^j$ by weighted
$x$-degree $i + 2j$; the lower bounds ($\ge 2$ for $P$, $\ge 1$ for $Q$)
encode realizability of $F$ as a polynomial map. The original box is kept
as an anchor, and every max degree is raised by 1 and then by 2:

| box | eqs | unknowns | rank | kernel | gauge rank | non-gauge | new families | time |
|---|---|---|---|---|---|---|---|---|
| base $P{:}[2,7]\ Q{:}[1,6]\ R{:}[0,4]$ | 57 | 43 | 28 | 15 | 9 | 6 | **0** | 9 s |
| +1 $\ P{:}[2,8]\ Q{:}[1,7]\ R{:}[0,5]$ | 65 | 55 | 33 | 22 | 13 | 9 | **0** | 11 s |
| +2 $\ P{:}[2,9]\ Q{:}[1,8]\ R{:}[0,6]$ | 74 | 69 | 40 | 29 | 14 | 15 | **0** | 27 s |

(Unknowns = box monomial coefficients of $(\delta P, \delta Q, \delta R)$
plus $\delta\kappa$; equations = coefficient equations of the linearized
Keller identity; all exact rational linear algebra. The base row reproduces
the documented kernel dimension 15 and gauge rank 9 of
`docs/NEW_COUNTEREXAMPLES.md` §3 — asserted.)

The gauge tangents are computed directly in $(P,Q,R)$-coordinates from the
generators of `docs/NEW_COUNTEREXAMPLES.md` §2 (source/target torus, source
shifts $\delta w = v^j$, $\delta v = w^j$, target shifts
$\delta P = Q^{j+2}R^j$, $\delta Q = P^i R^{2i-1}$); each is verified
exactly against the linearized Keller identity and against membership in
the kernel (asserts). Generators whose tangent leaves the box are skipped
and reported — as the boxes grow, more of them fit (9 → 13 → 14), which is
why the kernel grows without new geometry appearing.

### 2.2 Per-direction analysis and verdicts

For each box, the SVD-orthogonal complement of the gauge span inside the
kernel gives the candidate non-gauge deformation directions (6, 9, 15).
Each was subjected to two independent tests:

1. **Second-order integrability** (new relative to the original script):
   writing $\mathrm{res}(u_0 + t\xi) = t^2 q_2 + O(t^3)$, the direction
   integrates to second order iff $q_2 \in \mathrm{im}\, L$
   ($L$ = linearization at $u_0$; $q_2$ extracted exactly via a 6-point
   Vandermonde solve, since the residual is quintic in $t$). **Result:
   every one of the 30 non-gauge directions across the three boxes is
   obstructed already at second order** (relative defect of the best
   least-squares solution between $1.4\times10^{-3}$ and
   $4\times10^{-1}$, far above the $\sim10^{-6}$ float noise floor).
2. **Nonlinear continuation** (the arbiter also used by the original
   script): Gauss–Newton on the full 2D Keller system with displacement
   constraint $\xi \cdot (u - u_0) = 0.05$. Every direction either fails
   to converge or jumps to a distant solution ($|u - u_0| > 4$, i.e.
   $> 80\times$ the requested step — a far-away gauge image, not a branch
   through the base point).

**Verdict, all three boxes: RIGID modulo gauge.** No direction survived,
so the pursuit branch (exact Keller verification at a continued point +
numerical injectivity via `jcqft.fibers`) was never triggered — the code
path exists in `scripts/rigidity_boxes.py` (`pursue_survivor`) and is
exercised never, honestly. The second-order test strengthens the previous
evidence qualitatively: the obstruction is visible *algebraically at order
2*, not only through the behavior of a float Newton iteration.

### 2.3 What this does and does not establish

- Within the $(1,-1,-2)$ weight class and ansatz boxes up to
  $\deg_x P \le 9$, $\deg_x Q \le 8$, $\deg_x R \le 6$ — strictly
  containing everything searched before — the Alpöge–Mathew solution
  admits **no infinitesimal deformation beyond gauge, and every candidate
  direction is obstructed at second order**. The map continues to behave
  as an isolated, distinguished algebraic object rather than a point on a
  positive-dimensional moduli space.
- This is **numerical evidence, not proof**, in exactly two places: the
  float SVD split of the kernel into gauge + complement (exact: kernel
  dimension, gauge ranks, gauge $\subseteq$ kernel), and the
  second-order/continuation obstruction tests (float lstsq /
  finite-difference Gauss–Newton). Exact certification — Gröbner
  computation of the in-box solution ideal — remains open as
  `docs/OPEN_QUESTIONS.md` B4.
- Rigidity remains relative to the weight system: other gradings
  (`OPEN_QUESTIONS.md` B3) are untouched by this analysis. The prefilter
  of §1 is the intended cheap first pass for those searches.

## 3. Files

- `jcqft/prefilter.py` — `leading_part`, `infinity_prefilter`,
  `infinity_witness` (re-exported as `jcqft.infinity_prefilter`).
- `scripts/witten_prefilter.py` — validation + benchmark of §1 (~6 s,
  all claims asserted).
- `scripts/rigidity_boxes.py` — the three-box analysis of §2 (~50 s,
  all claims asserted; prints the table above).
