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

## Next (delegated to sub-agents, results to be merged here)

- Numerical monodromy of the 3 sheets around $\{p=0\}$ → `monodromy.py`,
  `docs/MONODROMY.md`.
- Reverse-engineering of the construction + search for new inequivalent
  counterexamples → `search_counterexamples.py`, `docs/NEW_COUNTEREXAMPLES.md`.
- Careful assessment of implications for rigorous QFT ($D=4$, pAQFT, field
  redefinitions) → `docs/QFT_IMPLICATIONS.md`.
