# The Symmetric (Variational) Jacobian Problem: Explicit Counterexamples in Dimension 6, Gradient No-Go in Dimension 3

*(A **symmetric Keller map** is $F = \nabla W$ with $\det\mathrm{Hess} W$
a nonzero constant — the "0D QFT with an action": the partition function
$\int e^{-W}$ has genuinely stationary phases, SUSY/Parisi–Sourlas
completions exist, coercivity is meaningful. `docs/PROBLEM.md` and
`docs/QFT_IMPLICATIONS.md` §5.3 record that the Alpöge–Mathew map $F$ is
**not** a gradient ($DF \neq DF^{\mathsf T}$); this document settles what
that means structurally. Every identity and every emptiness claim below is
proved by assertion in `scripts/symmetric_search.py` (72 checks):*

```
.venv/bin/python scripts/symmetric_search.py          # ~3.5 min, everything
.venv/bin/python scripts/symmetric_search.py --full   # + brute-force 6x6
                                                      #   determinants, larger
                                                      #   Groebner boxes, and
                                                      #   the n=3 deg-5 wall
                                                      #   attempt (hours)
```

*Machinery: `jcqft/reduction_w.py` (the $(1,-1,-m)$ equivariant reduction,
proved in `scripts/reduction_113.py`), `jcqft/gb_backend.py` (msolve /
Singular, 16 GB cap). Literature anchor: M. de Bondt, A. van den Essen,
"A reduction of the Jacobian Conjecture to the symmetric case", Proc. AMS
**133** (2005) 2201–2205 [dBvdE05].)*

## 1. Verdicts

1. **Variational Jacobian counterexamples exist, explicitly, in dimension
   6** (§3): two artifacts built from the Alpöge–Mathew map — the
   cotangent lift $W_6 = \bar\varphi \cdot F(\varphi)$ (rational, $\det
   \mathrm{Hess} W_6 \equiv -4$, three explicit rational witnesses)
   and the dBvdE twisted lift $\widetilde F = \mathrm{id} + \nabla f_H$
   (normalized, $\det(I + \mathrm{Hess} f_H) \equiv 1$, witnesses
   over $\mathbb{Q}(i)$). To our knowledge these are the first explicit
   symmetric/variational Jacobian counterexamples written down.
2. **The Alpöge–Mathew defect itself is non-variational in dimension 3, in
   the strongest linear sense** (§4): the space of matrices $K$ with
   $K\,DF(\varphi)$ symmetric for all $\varphi$ is exactly $\{0\}$, so no
   affine source/target change makes $F$ a gradient.
3. **Gradient no-go in the equivariant family** (§5): for every $m \ge 2$,
   the only gradient Keller maps equivariant for
   $\lambda\cdot(x,y,z) = (\lambda x, \lambda^{-1}y, \lambda^{-m}z)$ are
   linear maps and one explicit **tame** shear family — outright for every
   $m \neq 3$; at $m = 3$ up to a precisely-mapped box gap
   ($\deg_v \ge 3$ beyond the Gröbner boxes of §5.5). In particular at
   $m = 2$, the Alpöge–Mathew weight system, the exclusion is **complete**:
   the non-properness defect of the family is non-variational.
4. **Direct boxes** (§6): no non-injective symmetric Keller map exists for
   $\deg W \le 3$ in any dimension (midpoint identity), for $n = 2$ up to
   $\deg W = 6$, or for $n = 3$ up to $\deg W = 4$ (full coefficient box,
   exact over $\mathbb{Q}$). $n = 3$, $\deg W = 5$ is the honest wall.
5. **Coercivity** (§7): a real symmetric Keller counterexample is never
   coercive when $\kappa \le 0$ or $\deg W$ is odd; $\det\mathrm{Hess}
   W = $ const forces the Hessian of the leading form to be singular. The
   explicit $W_6$ fails coercivity in every one of these ways at once.

So the "0D theory with an action" exists — one dimension-doubling away from
the Alpöge–Mathew theory, and necessarily non-coercive; in the original
three fields there is provably no action to be found.

## 2. The de Bondt–van den Essen reduction: exact statement (Task A)

What [dBvdE05] proves (statement reconstructed and **re-proved from first
principles below**; the two bracketed details are NOTE-unverified against
the printed text, which we could not consult — everything actually *used*
is asserted in the script):

> **Theorem (dBvdE05, Thm 1.1).** If the Jacobian Conjecture holds for all
> symmetric Keller maps $x + \nabla f$, $f \in \mathbb{C}[x_1,\dots,x_n]$,
> all $n$, then it holds for all Keller maps, all $n$. The reduction
> **doubles the dimension**: to a Keller map $F = x + H$ on
> $\mathbb{C}^n$ it associates
> $$
> f_H(x, y) \;=\; -\,i \sum_{j=1}^{n} H_j(x + iy)\, y_j
> \;\in\; \mathbb{C}[x_1,\dots,x_n,y_1,\dots,y_n],
> $$
> and $\widetilde F = \mathrm{id} + \nabla f_H$ on $\mathbb{C}^{2n}$ is a
> symmetric Keller map, invertible iff $F$ is. [NOTE-unverified: the
> paper's sign/scaling convention for $f_H$ may differ; ours is fixed by
> the assertions below. NOTE-unverified: dBvdE05 further combine this with
> the Bass–Connell–Wright degree reduction to reach $f$ homogeneous of
> degree 4 with $\mathrm{Hess} f$ nilpotent, at the cost of a further
> dimension blow-up; we use only the twisted lift.]

The proof is three exact identities, all asserted (script §A1–A3) with
$g_H := -i\,\langle y, H(x)\rangle$ and
$S = \begin{pmatrix} I & -iI \\ 0 & I \end{pmatrix}$ (so $f_H = g_H \circ
S^{-1}$, $\det S = 1$):

$$
\mathrm{Hess} f_H = S^{-\mathsf T}\,
\bigl(\mathrm{Hess} g_H\bigr)(S^{-1}\cdot)\, S^{-1},
\qquad
S^{\mathsf T} S + \mathrm{Hess} g_H =
\begin{pmatrix} I + B & -i\,(I + JH)^{\mathsf T} \\ -i\,(I + JH) & 0
\end{pmatrix},
$$

with $B_{jk} = -i \sum_l y_l\, \partial_j\partial_k H_l$, together with the
**block-determinant lemma** $\det\begin{pmatrix} M & C \\ A & 0
\end{pmatrix} = (-1)^n \det A \,\det C$ (asserted on fully generic blocks
for $n = 2, 3$; the $n$ row-block transpositions prove it for every $n$).
Chaining them:

$$
\boxed{\ \det\bigl(I + \mathrm{Hess} f_H\bigr)
 = \det\bigl(I + JH\bigr)^2 \ }
$$

— Keller-ness transfers, squared. Invertibility transfers through the
asserted conjugation formula (script §A3)

$$
S^{-1} \circ \widetilde F \circ S
= \bigl(\,x + H(x),\ \ (I + JH(x))^{\mathsf T}\, y - i\,H(x)\,\bigr),
$$

whose first block is $F$ itself and whose second block is linear in $y$
with everywhere-invertible coefficient matrix: $\widetilde F$ is injective
iff $F$ is. The reduction is **constructive and over $\mathbb{C}$** (it
needs $i$; the real form is not addressed here). **Dimension bookkeeping:
the Alpöge–Mathew counterexample in dimension 3 produces an explicit
symmetric counterexample in dimension 6.**

## 3. The two explicit dimension-6 artifacts (Task A, constructive)

### 3.1 The cotangent lift: $W_6 = \bar\varphi \cdot F(\varphi)$

The first-order action of `docs/PROBLEM.md` — the auxiliary-field trick —
*is* a symmetric counterexample, without any normalization. With
$\varphi = (x,y,z)$, $\bar\varphi = (\bar x, \bar y, \bar z)$ and $F$ the
Alpöge–Mathew map, all asserted:

- $W_6 := \bar x F_1 + \bar y F_2 + \bar z F_3$, $\deg W_6 = 8$;
- $\mathrm{Hess} W_6 = \begin{pmatrix} \sum_k \bar\varphi_k
  \mathrm{Hess} F_k & JF^{\mathsf T} \\ JF & 0 \end{pmatrix}$, so by
  the block lemma
  $\det \mathrm{Hess} W_6 = -(\det JF)^2 = -4$ — **constant**;
- $\nabla W_6 = (JF(\varphi)^{\mathsf T}\bar\varphi,\ F(\varphi))$ is
  **3:1 over a rational point**: since $F(a) = F(b) = F(c)$ on the
  Alpöge–Mathew triple and $JF^{\mathsf T}$ is invertible everywhere, a
  linear solve matches the first block. The three witnesses are rational:
  $$
  \bigl(0, 0, -\tfrac14;\; 1, 1, 1\bigr),\quad
  \bigl(1, -\tfrac32, \tfrac{13}{2};\; -\tfrac{41}{2}, -\tfrac12,
  \tfrac{19}{16}\bigr),\quad
  \bigl(-1, \tfrac32, \tfrac{13}{2};\; -\tfrac{63}{2}, 1,
  -\tfrac{35}{16}\bigr),
  $$
  all mapped by $\nabla W_6$ to $\bigl(\tfrac54, 1, 1, -\tfrac14, 0,
  0\bigr)$ (asserted). Injectivity transfer is an iff, by the same
  block-triangular structure.

This is a counterexample to the gradient/Hessian formulation of the
Jacobian Conjecture in dimension 6 with **rational data of degree 8** —
and its potential is exactly the QFT action $\bar\varphi\cdot F(\varphi)$:
the conjugate fields are the cotangent directions, and the non-properness
of $F$ is inherited by the action's gradient flow.

### 3.2 The normalized twisted lift: $\widetilde F = \mathrm{id} + \nabla f_H$ over $\mathbb{Q}(i)$

For the normalized form ($\kappa = 1$, identity linear part), apply §2 to
$H := L^{-1}F - \mathrm{id}$ where $L = DF(0)$ (so $JH(0) = 0$,
$\det(I + JH) = 1$; asserted). Then, all asserted:

- $f_H = -i\sum_j H_j(x + iy)\,y_j$, $\deg f_H = 8$, coefficients in
  $\mathbb{Q}(i)$, $\mathrm{Hess} f_H(0) = 0$;
- $\det(I + \mathrm{Hess} f_H) = 1$ **identically**: $\widetilde F =
  \mathrm{id} + \nabla f_H$ is a normalized symmetric Keller map on
  $\mathbb{C}^6$;
- three distinct witnesses with common image $(0, 0, -\tfrac14, 0, 0, 0)$,
  obtained from the Alpöge–Mathew triple through the conjugation formula
  (second block solved linearly, then pushed through $S$):
  $$
  \bigl(0,0,-\tfrac14;\,0,0,0\bigr),\quad
  \bigl(-\tfrac{699}{32}, \tfrac{63}{16}, \tfrac{369}{2};\,
        -\tfrac{731}{32}i, \tfrac{87}{16}i, 178\,i\bigr),\quad
  \bigl(\tfrac{699}{32}, -\tfrac{63}{16}, \tfrac{369}{2};\,
        \tfrac{731}{32}i, -\tfrac{87}{16}i, 178\,i\bigr).
  $$

### 3.3 What the twisted lift cannot do from here

[dBvdE05, Lemma 1.2] states $\mathrm{Hess} f_H$ nilpotent iff $JH$
nilpotent. Both fail here, provably (script §A4):

- $\mathrm{tr}\bigl((\mathrm{Hess} f_H)^2\bigr) \neq 0$ and
  $\mathrm{charpoly}(JH) \neq t^3$: **neither is nilpotent**;
- no affine normalization helps: $x + H'$ with $JH'$ nilpotent would need
  an invertible $C$ with $C\,JF(u) - I$ nilpotent for all $u$, hence
  $\mathrm{tr}(C\,JF) = 3$ identically — asserted **inconsistent**
  as a linear system (`linsolve` returns the empty set).

So the dBvdE nilpotent/homogeneous normal form (their Cor. 1.3) provably
cannot be reached from the Alpöge–Mathew map by the twist alone; it
requires the Bass–Connell–Wright stabilization first, with its dimension
blow-up (no explicit dimension computed here).

## 4. Dimension 3: the strongest no-go for the AM map itself (Task A→B)

$B\,F(A\varphi) + c$ is a gradient map iff $B\,JF(A\varphi)\,A$ is
symmetric for all $\varphi$, iff (substituting $u = A\varphi$ and setting
$K := A^{-\mathsf T} B$) $K\,JF(u)$ is symmetric for all $u$. The set of
such $K$ is a linear space; asserted (script §A5):

> **Proposition.** $\{K \in \mathbb{C}^{3\times 3} : K\,JF(u) =
> (K\,JF(u))^{\mathsf T}\ \forall u\} = \{0\}$. Hence the Alpöge–Mathew
> map is not equivalent to **any** gradient map under invertible affine
> source and target transformations.

This strengthens `QFT_IMPLICATIONS.md` §5.3 ($DF \neq DF^{\mathsf T}$) from
"not a gradient as given" to "not a gradient in any affine frame".

## 5. Gradient structure vs $\mathbb{C}^*$-equivariance (Task B)

Throughout: source action $\lambda\cdot(x,y,z) = (\lambda x,
\lambda^{-1}y, \lambda^{-m}z)$, source weights $\mathbf w = (1,-1,-m)$,
and the invertibility constraint that the component weights
$(d_1,d_2,d_3)$ of $F$ are a permutation of $(-m,-1,1)$
(`scripts/reduction_113.py`).

### 5.1 Weight bookkeeping: the branch table

If $F = \nabla W$ then $F_i = \partial_i W$, and the weight-$e$ graded
piece $W_e$ contributes to $F_i$ with weight $e - w_i$. Matching
$d_i = e_i - w_i$ over the six permutations (asserted, script §B1):

| $(d_1,d_2,d_3)$ | $(e_1,e_2,e_3)$ | coincidences | potential |
|---|---|---|---|
| $(-m,-1,1)$ | $(1{-}m,-2,1{-}m)$ | $e_1 = e_3$, all $m$; all equal iff $m=3$ | $W = x^{1-m}B(v) + \beta y^2/2$; at $m{=}3$: $W = x^{-2}S(w,v)$ |
| $(-1,1,-m)$ | $(0,0,-2m)$ | $e_1 = e_2$, all $m$ | $W = A(w) + \delta z^2/2$ |
| $(1,-m,-1)$ | $(2,-m{-}1,-m{-}1)$ | $e_2 = e_3$, all $m$ | $W = \alpha x^2/2 + \delta yz + \gamma y^{m+1}$ (shear) |
| other three | all distinct (sporadic coincidences only at $m=3$, subsumed) | — | diagonal $F = (f_1(x), f_2(y), f_3(z))$ |

($w = xy$, $v = x^m z$, the free invariants.) The **only branch with a
nontrivial two-variable potential is $m = 3$**, $W = x^{-2}S(w, v)$ — a
*rational* potential whose gradient is polynomial exactly on the
admissibility box $j + 3k \ge 3$ for the monomials $w^j v^k$ of
$P = E(S)$ below. The sporadic pairwise coincidences of the last three
permutations happen at $m = 3$ only (asserted) and are subsumed; the two
weight-infeasible univariate slots ($F_3 \in \mathbb{C}[z]$ of weight
$+1$; $F_2 \in \mathbb{C}[y]$ of weight $-m$, which forces
$F_2 \propto y^m$ and kills $DF(0)$) are empty.

### 5.2 Symmetry of $DF$ in reduced coordinates

In coordinates $(x, w, v)$ the derivations are $\partial_y = x\,
\partial_w$, $\partial_z = x^m \partial_v$, $\partial_x|_{y,z} =
\partial_x + (w/x)\partial_w + (mv/x)\partial_v$. For the equivariant
normal form $F = (P/x^m,\ Q/x,\ xR)$, the three symmetry residuals of
$DF$ are asserted to be exactly

$$
\partial_y F_1 - \partial_x F_2 = x^{1-m} P_w - x^{-2}\,(w Q_w + m v Q_v - Q),
\quad
\partial_z F_1 - \partial_x F_3 = P_v - (R + w R_w + m v R_v),
$$
$$
\partial_z F_2 - \partial_y F_3 = x^{m-1} Q_v - x^2 R_w .
$$

For $m \neq 3$ the $x$-powers are unbalanced: symmetry forces $P_w = R_w
= 0$ and $wQ_w + mvQ_v = Q$ — precisely the degenerate branches of §5.1.
At $m = 3$ the system is $x$-free and integrates (Poincaré on the third
residual, then the first two; asserted):

> **Proposition ($m = 3$).** $DF$ symmetric $\iff$
> $(P, Q, R) = (E(S),\, S_w,\, S_v)$ for a single $S(w,v)$, where
> $E = w\,\partial_w + 3v\,\partial_v - 2$, i.e.
> $F = \nabla\bigl(x^{-2} S(w,v)\bigr)$. The kernel of $E$ on admissible
> monomials is $\mathrm{span}\{w^2\}$, inside the gauge slice — no
> integration-constant ambiguity survives the box (asserted).

### 5.3 The degenerate branches: complete classification, every $m$

All asserted (script §B3):

- **diagonal**: $\det DF = f_1'f_2'f_3' = \kappa$ has leading term
  $\prod a_i d_i\, x^{d_1-1}y^{d_2-1}z^{d_3-1} \neq 0$, so all $d_i = 1$:
  **linear**.
- **$A(w)$-slice** ($W = A(xy) + \delta z^2/2$):
  $\det\mathrm{Hess} = -\delta\,A'(A' + 2wA'')\big|_{w=xy}$, and the
  Wronskian $u(u + 2wu') = (wu^2)'$ with $u := A'$ integrates the Keller
  condition to $w u^2 = -(\kappa/\delta)\,w$: $u$ constant (degree kill
  asserted): **linear**.
- **$B(v)$-slice** ($W = x^{1-m}B(v) + \beta y^2/2$, reduced data
  $P = (1{-}m)B + mvB'$, $Q = \beta w$, $R = B'$):
  $\det M = \beta\,[\,m(m{-}1)(B - vB')B'' - B'^2 - 2mvB'B''\,]$
  (asserted against `reduction_w.det_m`); on $B = a v^k + \dots$ the
  leading coefficient is $-k\,[\,m(m{-}1)(k{-}1)^2 + k + 2mk(k{-}1)\,]
  a^2 \neq 0$ for $k \ge 2$, $m \ge 2$ (positivity asserted by shift):
  $B$ **linear**, $F$ **linear**.
- **shear family** ($W = \alpha x^2/2 + \delta yz + \gamma y^{m+1}$):
  $$
  \nabla W = \bigl(\alpha x,\ \delta z + (m{+}1)\gamma y^m,\ \delta y
  \bigr), \qquad \det\mathrm{Hess} W = -\alpha\delta^2,
  $$
  a Keller map for **every** $m$ and every $\gamma$ — and a **tame
  automorphism**: explicit inverse verified by composition; the infinity
  prefilter (`jcqft/prefilter.py`) is survived only through its documented
  nonlinear-automorphism false-positive class (asserted at $m = 2, 3$).
  Component weights $(1, -m, -1)$: a permutation of $(1,-1,-m)$, as
  required.

### 5.4 The $m = 3$ potential slice: a rigid Monge–Ampère problem

Substituting $(P,Q,R) = (E(S), S_w, S_v)$ into the reduced Keller
determinant of `reduction_w` gives the asserted identity

$$
\det M = -2\,\bigl(wS_w + 6vS_v - 3S\bigr)
\bigl(S_{ww}S_{vv} - S_{wv}^2\bigr)
- \bigl(S_v,\ S_w\bigr)\mathrm{Hess}S\,\bigl(S_v,\ S_w
\bigr)^{\mathsf T} \;=\; \kappa,
$$

a two-variable Monge–Ampère-type PDE. $\kappa \neq 0$ makes $DF(0)$
invertible automatically; the residual scalings gauge-fix $[v]S = 1$,
$[w^2]S = \tfrac12$, whence $\kappa = -1$ and the **anchor** $S = v +
w^2/2 \mapsto \det M = -1$ (the linear map $(z, y, x)$; asserted). Writing
$S = \sum_k s_k(w)\,v^k$, $K = \deg_v S$:

- **$K = 0$**: $\det M \equiv 0 \neq \kappa$ ($F_3 = xS_v = 0$): empty.
- **$K = 1$**: the $[v]$-equation $2ws_1'^3 + 4s_1s_1'^2 - s_1^2s_1'' = 0$
  has leading coefficient $a^3 d(d{+}1)(2d{+}1) \neq 0$ on $s_1 = aw^d$:
  $s_1$ constant; then $-s_1^2 s_0'' = \kappa$ forces $s_0$ exactly
  quadratic: **only the linear map**, all $w$-degrees (asserted chain).
- **$K = 2$**: cascade $c_3 = -40a^2s_1'' \Rightarrow s_1 = s + \beta w$;
  $c_2 \Rightarrow s_0'' = \tfrac{3\beta^2}{10a}$ const;
  $c_1 = -\tfrac45\beta^2(3\beta w + s) \Rightarrow \beta = 0 \Rightarrow
  s_0'' = 0 \Rightarrow q_0'(0) = 0$: $DF(0)$ singular. **Empty, all
  $w$-degrees** (asserted cascade).
- **Symbolic $K$** — the top rows are rigid for *every* $K \ge 1$
  (asserted with $K$ a sympy `Symbol`):
  - **L1**: $c_{3K-2}$ involves only $s_K$, and on $s_K = a w^\delta$
    equals $a^3K\delta\,(3K{+}\delta{-}2)(4K{+}2\delta{-}3)\,
    w^{3\delta-2}$, both brackets $> 0$ for $K, \delta \ge 1$:
    $s_K$ **constant**;
  - **L2**: $c_{3K-3} = -K(3K{-}2)(4K{-}3)\,a^2\,s_{K-1}''$:
    $s_{K-1}$ **linear**;
  - **L3**: $c_{3K-4} = 3a(K{-}1)(4K^2{-}7K{+}2)\,t_1^2 -
    a^2K(3K{-}2)(4K{-}3)\,s_{K-2}''$ (with $s_{K-1} = t_0 + t_1w$):
    $s_{K-2}$ **quadratic** with fixed leading coefficient
    $\propto t_1^2$.

### 5.5 In-box Gröbner closure for $K = 3, 4$

Presolving with L1–L2 (top row constant, second row linear) and the gauge,
the remaining Keller system plus the Rabinowitsch certificate
$[v^K]S \neq 0$ is asserted to be the **unit ideal** (msolve, exact over
$\mathbb{Q}$, seconds each):

| box | default run | `--full` |
|---|---|---|
| $K = 3$, $\deg_w \le 4, 6, 8$ | EMPTY | $\deg_w \le 10$ |
| $K = 4$, $\deg_w \le 4, 6$ | EMPTY | $\deg_w \le 8$ |

> **Theorem (gradient no-go, $(1,-1,-m)$ family).** For every integer
> $m \ge 2$, every gradient Keller map equivariant for the weights
> $(1,-1,-m)$ is either linear or a member of the tame shear family
> $W = \alpha x^2/2 + \delta yz + \gamma y^{m+1}$ — outright for every
> $m \neq 3$; for $m = 3$ outright at $\deg_v S \le 2$, in-box at
> $\deg_v S = 3$ ($\deg_w \le 8$) and $4$ ($\deg_w \le 6$), with the top
> three $v$-rows rigid for all $\deg_v$ (L1–L3).
>
> In particular ($m = 2$, complete): **the non-properness defect of the
> Alpöge–Mathew family is non-variational** — no equivariant potential
> exists for any nonlinear non-tame member, in any $w$-degree.

*Gap qualification*: $m = 3$, $\deg_v S = 3$ with $\deg_w > 8$;
$\deg_v S = 4$ with $\deg_w > 6$; $\deg_v S \ge 5$ beyond L1–L3. All other
branches are closed with $m$ symbolic.

## 6. Direct boxes for symmetric Keller maps (Task C)

### 6.1 Normalization (WLOG over $\mathbb{C}$, complete for the boxes)

If $\nabla W(a) = \nabla W(b)$, $a \neq b$, translate so the witness pair
is $\pm e/2$; the substitution $x \to Ax$ acts on the Hessian by
congruence ($\mathrm{Hess} \to A^{\mathsf T}(\mathrm{Hess})A$),
preserves the pair form, and over $\mathbb{C}$ normalizes the (invertible,
symmetric) $\mathrm{Hess} W(0)$ to $I$ — hence quadratic part
$\sum x_i^2/2$ and $\kappa = 1$; the linear part of $W$ drops out of the
witness equations. Degrees are preserved, so **unit-ideal certificates in
the normalized box are complete** for the stated degree.

### 6.2 $\deg W \le 3$, any dimension: no search needed

The **midpoint identity** (asserted on full generic cubics, $n = 2, 3, 4$;
the polarization argument is dimension-free):

$$
\nabla W(a) - \nabla W(b) = \mathrm{Hess} W\!\Bigl(\frac{a+b}{2}
\Bigr)\,(a - b) \qquad (\deg W \le 3),
$$

so $\det\mathrm{Hess} \equiv \kappa \neq 0$ forces injectivity.
Every symmetric Keller map of degree $\le 3$ is injective, in every
dimension. (Contrast: the *non-symmetric* cubic case is the full JC, by
Bass–Connell–Wright.)

### 6.3 $n = 2$: calibration

Witness-pair emptiness asserted for $\deg W \le 4, 5, 6$ (msolve, exact
over $\mathbb{Q}$, seconds), with both controls:

- dropping the Keller equations leaves the witness system **non-unit**
  (witness pairs exist for generic quartics) — the pipeline cannot
  spuriously report emptiness;
- the Keller variety alone is **nonempty**: nonlinear symmetric Keller
  maps exist, e.g. $W = xy + f(x)$ with $\det\mathrm{Hess} = -1$ and
  the explicit inverse $(a, b) \mapsto (b,\ a - f'(b))$ (asserted) — the
  gradient analogue of a triangular automorphism.

### 6.4 $n = 3$: degree 4 closed, degree 5 the wall

- $\deg W \le 4$: the **full 25-coefficient box** (10 cubic + 15 quartic
  coefficients, 29 unknowns with witness vector and Rabinowitsch
  variable), three witness queries, msolve exact over $\mathbb{Q}$,
  ~3 min: **EMPTY**. Combined with §6.2: **no symmetric Keller
  counterexample exists in dimension 3 through degree 4.**
- $\deg W \le 5$: 46 coefficients, 50 unknowns, 222 equations. This is
  the honest wall: even the **mod-$p$ screen** ($p = 1073741789$,
  memory-light coefficients) exceeded the 16 GB msolve F4 cap after
  17 min on this machine (2026-07-26); the exact run over $\mathbb{Q}$
  is strictly harder. Status: **open, certified wall** (`--full`
  re-attempts it and reports honestly).
- $n = 4$, $\deg W = 4$ (45 coefficients): not attempted — strictly
  larger than the $n = 3$, $\deg 5$ wall. Note $\deg W \le 3$ *is*
  closed for $n = 4$ by §6.2.

## 7. Coercivity: why a variational counterexample cannot be a good action (Task D)

For real data, all three mechanisms are relevant; the first two are
elementary and exact, the third is asserted for $W_6$ and proved in
general by one line of degree bookkeeping:

1. **Sign of $\kappa$.** A coercive $C^2$ potential $W$ attains a minimum
   $x_0$; nondegeneracy ($\kappa \neq 0$) makes $\mathrm{Hess}
   W(x_0)$ positive definite, so $\kappa = \det\mathrm{Hess} W > 0$.
   Any symmetric Keller data with $\kappa \le 0$ — e.g. $W_6$, with
   $\kappa = -4$ (asserted) — is **never coercive**.
2. **Parity.** Odd $\deg W$ kills coercivity outright ($W(tx) \to
   -\infty$ along some ray). The cotangent lift is worse: $W_6$ is
   **affine in $\bar\varphi$** (asserted) — a first-order action is never
   coercive along the conjugate directions, for any base map.
3. **The leading form is degenerate.** If $\deg W = D \ge 3$, the entries
   of $\mathrm{Hess} W$ have degree $\le D - 2$ with top-degree
   parts $\mathrm{Hess} W_D$ ($W_D$ the leading form), so the
   degree-$n(D-2)$ part of $\det\mathrm{Hess} W$ is
   $\det\mathrm{Hess} W_D$; constancy of $\det\mathrm{Hess} W$
   forces $\det\mathrm{Hess} W_D \equiv 0$. So the leading form of
   *every* nonquadratic symmetric Keller potential is a degenerate
   (parabolic) form — the growth of $W$ at infinity is never uniformly
   convex. Asserted instance: the leading form of $W_6$ is
   $\bar x\, x^3y^3z$, whose Hessian has an identically vanishing
   $\bar y$-row.

**Synthesis.** Variational Jacobian counterexamples: **exist in dimension
6, explicitly and rationally** (§3); **excluded in dimension 3** in the
controlled equivariant families (§5, complete at $m \neq 3$) and in the
degree boxes (§6, through degree 4); **open in between** (dimensions 4, 5,
and $n = 3$ beyond degree 4). The physics reading: the "0D theory with an
action" exists one cotangent doubling away — the auxiliary-field
first-order action $\bar\varphi \cdot F(\varphi)$ is itself the potential —
but such an action is *necessarily* non-coercive, with degenerate leading
form, so the pathology (stationary phases escaping to infinity, fiber
jumping) persists verbatim in the variational setting. The coercivity
screen of the paper is not circumvented: an action exists, a *good*
(coercive, stable) action does not.

## 8. Honest limitations

1. **Literature statement.** The exact printed hypotheses and conventions
   of [dBvdE05] (sign of $f_H$, their real-case remarks, the homogeneous
   degree-4 corollary) are reconstructed, not checked against the text;
   flagged NOTE-unverified in §2. Everything *used* is proved by
   assertion from first principles.
2. **The B gap** (§5.5): $m = 3$, $\deg_v S = 3$ / $4$ beyond $w$-degree
   8 / 6, and $\deg_v S \ge 5$ (top three rows rigid, remainder open).
3. **The C wall** (§6.4): $n = 3$, $\deg W = 5$ unresolved either way;
   $n = 4$ nonlinear boxes untouched.
4. **Block lemma scope.** Asserted on generic blocks for $n = 2, 3$ only
   (the row-swap proof is general but not formalized for symbolic $n$);
   all uses here are $n = 3$.
5. **Minimality open.** Dimension 6 is what the constructive reduction
   gives from dimension 3; whether symmetric counterexamples exist in
   dimensions 4 or 5 is open (they cannot exist in dimension $\le 2$ if
   JC$_2$ holds, which is itself open; for $\deg W \le 3$ they are
   excluded in every dimension by §6.2).
6. **Coercivity argument scope.** §7.1 assumes real coefficients and
   $C^2$; the complex artifacts are commented on via their real forms
   only qualitatively, except where asserted ($\kappa$, affineness,
   leading form).

## 9. QFT reading

The Alpöge–Mathew theory has no action in its own field space — not even
after affine field redefinitions (§4), and not anywhere in its equivariant
deformation class (§5). But its **phase-space (first-order) formulation is
an action**: $W_6 = \bar\varphi\cdot F(\varphi)$ on the doubled fields, and
that action is an honest variational Jacobian counterexample — gradient,
constant Hessian determinant, non-injective, with rational stationary
data. The doubling is exactly the auxiliary-field trick of
`docs/PROBLEM.md`; what this document adds is that the trick is *forced*
(dimension 3 admits no potential) and *sufficient* (dimension 6 does), and
that any such action is necessarily non-coercive with parabolic leading
form (§7) — the variational formulation buys stationary-phase technology,
never stability. The remaining structural question is the minimal
dimension of a variational counterexample (4, 5, or 6), and whether the
$n = 3$ degree-5 box hides a *symmetric* non-injective map that the
equivariant no-go does not see.
