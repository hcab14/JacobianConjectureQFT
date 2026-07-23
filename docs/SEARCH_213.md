# The Weight System (2, −1, −3): Exact Reduction and Search

*(Machinery: `jcqft/reduction_213.py`. Every identity and every emptiness
claim below is proved by assertion in `scripts/search_213.py`:*

```
.venv/bin/python scripts/search_213.py            # ~40 s, all core claims
.venv/bin/python scripts/search_213.py --full     # + CE2 box no-go, ~8 min
.venv/bin/python scripts/search_213.py --deg2     # experimental, hours
```

*This executes the program of `docs/SEARCH_STRATEGIES.md` §1.1 for the
weight system (2,−1,−3) suggested in `docs/NEW_COUNTEREXAMPLES.md` §5.)*

## 1. Why this weight system is interesting

For $\lambda.(x,y,z) = (\lambda^2 x,\ \lambda^{-1} y,\ \lambda^{-3} z)$ the
orbit space has **two** orbifold strata instead of Alpöge–Mathew's one:

- the target $x$-axis has stabilizer $\mathbb{Z}_2$ ($\lambda^2 = 1$) — the
  same 2:1 mechanism that powers the Alpöge–Mathew counterexample;
- the target $z$-axis has stabilizer $\mathbb{Z}_3$ ($\lambda^3 = 1$) — a
  free orbit mapping into it would be **3:1**, giving cube-root escape and a
  $\mathbb{Z}_3$ monodromy class, the "different global anomaly" scenario of
  `docs/SEARCH_STRATEGIES.md` §2.2.

Unlike the $(1,-1,-m)$ family, the invariant ring is **not free**, so the
reduction of `jcqft/reduction_w.py` does not apply and had to be rebuilt.

## 2. Invariant theory (proved exactly)

Monomials $x^i y^j z^k$ are invariant iff $2i - j - 3k = 0$. The invariant
ring is

$$
R = \mathbb{C}[u_1, u_2, u_3]/(u_2^2 - u_1 u_3),\qquad
u_1 = xy^2,\quad u_2 = x^2yz,\quad u_3 = x^3z^2,
$$

an $A_1$ quadric cone (2-dimensional; the relation is verified and the
quotient presentation is checked through the **unique normal form**
$u_1^a u_2^b u_3^c$, $b \in \{0,1\}$, with $b = j \bmod 2$, $a = (j-b)/2$,
$c = (k-b)/2$ — existence and uniqueness verified exhaustively for all
invariant monomials of degree ≤ 12).

**Weight modules.** The space $M_d$ of weight-$d$ polynomials is an
$R$-module. Verified exhaustively (degree ≤ 12), with syzygy-free direct
splittings:

| weight | module | generators | canonical split |
|---|---|---|---|
| $+2$ | $M_2 = x\,R$ | $x$ | free of rank 1 |
| $-1$ | $M_{-1} = y\,R + xz\,R$ | $y,\ xz$ | $y\,R \oplus xz\,\mathbb{C}[u_3]$ |
| $-3$ | $M_{-3} = z\,R + y^3 R$ | $z,\ y^3$ | $z\,R \oplus y^3\,\mathbb{C}[u_1]$ |

(Syzygies: $u_1\cdot xz = u_2\cdot y$, $u_2\cdot xz = u_3\cdot y$;
$u_1^2\cdot z = u_2\cdot y^3$, $u_1u_2\cdot z = u_3\cdot y^3$. The
canonical splits absorb them: monomials with $j \ge 1$ vs. $j = 0$ in
$M_{-1}$, with $k \ge 1$ vs. $k = 0$ in $M_{-3}$.)

**Component weights (question (a)).** $x, y, z$ carry the pairwise distinct
weights $2, -1, -3$, and a linear monomial of weight $d$ can only occur in a
component of weight exactly $d$. Invertibility of $DF(0)$ requires each
variable to occur linearly somewhere, so with only three components the
component weights are forced to be a **permutation of $(2,-1,-3)$** — no
other assignment admits an invertible linear part — and $DF(0)$ is then
diagonal (each component contains exactly one linear variable). Fixing the
permutation (a target relabeling) gives the normal form of every candidate:

$$
F = \bigl(\,x\,A(u),\ \ y\,B(u) + xz\,C(u_3),\ \ z\,D(u) + y^3 E(u_1)\,\bigr),
\qquad DF(0) = \mathrm{diag}(A(0),\,B(0),\,D(0)).
$$

Every polynomial choice of $(A,B,C,D,E)$ assembles to a polynomial map — no
polynomiality side conditions, unlike $(1,-1,-m)$ where monomials had to
satisfy $j + mk \ge m$.

## 3. The exact reduced Keller identity (question (c))

This is a **chart-based reduction** (dense chart, denominators cleared),
as anticipated: the C* action is free on $\{y \neq 0\}$, and on
$\{x\neq0,\ y\neq 0\}$ the slice parametrization

$$
\Phi(t, p, q) = \Bigl(\frac{p}{t^2},\ t,\ \frac{q\,t^3}{p^2}\Bigr),
\qquad t = y,\ p = u_1,\ q = u_2,
$$

pulls the invariants back to exactly $(u_1, u_2, u_3)\circ\Phi =
(p, q, q^2/p)$ — the cone in its chart $u_3 = u_2^2/u_1$. Writing the chart
data

$$
\tilde a = A\Bigl(p, q, \tfrac{q^2}{p}\Bigr),\qquad
\tilde b = B(\cdot) + \tfrac{q}{p}\,C\Bigl(\tfrac{q^2}{p}\Bigr),\qquad
\tilde e = \tfrac{q}{p^2}\,D(\cdot) + E(p),
$$

one gets $F\circ\Phi = (p\tilde a/t^2,\ t\,\tilde b,\ t^3\tilde e)$ and
$\det D\Phi = -t/p^2$, hence by the chain rule the

$$
\boxed{\ \det DF \;=\; \Delta(p,q)
\;=\; p^2\Bigl(2\,p\tilde a\,J_2(\tilde b,\tilde e)
\;+\; \tilde b\,J_2(p\tilde a,\tilde e)
\;-\; 3\,\tilde e\,J_2(p\tilde a,\tilde b)\Bigr)\ }
$$

with $J_2$ the Jacobian in $(p,q)$: a rational function of the two
invariant chart coordinates alone, with denominator a pure power of $p$.
The coefficients $(2, 1, -3)$ are the $t$-exponents $(-2, 1, 3)$ of
$F\circ\Phi$ in the general pattern
$k_1S_1J_2(S_2,S_3) - k_2S_2J_2(S_1,S_3) + k_3S_3J_2(S_1,S_2)$ — the direct
analogue of `reduction_w`'s $\det M$ (there $k = (-m,-1,1)$; the compact
form $J_2(PR^m, QR) = \kappa R^m$ does not survive here because the
exponents $-k_1/k_2 = 2$, $-k_3/k_2 = -3$ make the corresponding product
combination rational, not polynomial).

**Status of rigor.** The identity is proved in `scripts/search_213.py` with
*undetermined functions* (sympy `Function` objects — an identity of
differential polynomials, not a spot check), and independently re-verified
on the fully symbolic degree-1 ansatz by comparing the raw 3-D
$\det DF$ with $\Delta(u_1, u_2)$, and against the intrinsic normal-form
determinant on the cone (`det_df_normal_form`, computed by toric
rewriting). Since $\Phi$ parametrizes a dense open set,

$$
\det DF \equiv \kappa
\quad\Longleftrightarrow\quad
\operatorname{numer}\bigl(\Delta - \kappa\bigr) \equiv 0 \ \text{in}\
\mathbb{C}[p,q],
$$

which is the finite polynomial system the search solves. The system is
**trilinear**: each equation has degree ≤ 1 in each of the blocks
$(A\,|\,B,C\,|\,D,E)$.

**Gauge fixing.** Target scalings normalize $A(0) = B(0) = D(0) = 1$, hence
$\kappa = \det DF(0) = 1$ (no lost generality: $\kappa \ne 0$ rescales
away). The constant $xz$-shear $y \to y + s\,xz$ normalizes $C(0) = 0$
(box-preserving). $E(0)$ **cannot** be gauged away box-preservingly and is
kept as an unknown. A residual 2-torus
$(u_1,u_2,u_3) \to (\mu_1 u_1, \mu_2 u_2, \mu_2^2\mu_1^{-1}u_3)$ survives
and is used below to normalize witness points.

## 4. Orbifold mechanisms (2:1 **and** 3:1)

A free orbit over the cone point $(p_0, q_0)$, $p_0 \ne 0$, maps

- **2:1** onto the target $x$-axis iff $\tilde b = \tilde e = 0 \neq \tilde a$
  at $(p_0,q_0)$ — image $(p_0\tilde a/t^2, 0, 0)$, preimages
  $\Phi(\pm t_0,\cdot)$;
- **3:1** onto the target $z$-axis iff $\tilde a = \tilde b = 0 \neq \tilde e$
  — image $(0,0,t^3\tilde e)$, preimages $\Phi(\omega^k t_0,\cdot)$,
  $\omega^3 = 1$;

plus the analogous conditions on the $y=0$ stratum (free orbits with
$u_1 = u_2 = 0$, $u_3 \neq 0$): 2:1 iff $C = D = 0 \neq A$ at
$(0,0,u_3)$, 3:1 iff $A = C = 0 \neq D$. The $x=0$ stratum supports no
mechanism ($F_1$ vanishes identically there while $F_2 = y$).
Any Keller solution with such a locus is automatically a counterexample.
This is one algebraic condition system per mechanism — the exact analogue
of "common zero of $(Q,R)$ off $\{P=0\}$" in `docs/NEW_COUNTEREXAMPLES.md`
§1.3.

## 5. Search results

### 5.1 Degree-1 box: complete solution (theorem for the box)

Box: $A, B, D$ of degree ≤ 1 in $(u_1,u_2,u_3)$, $C = c_3u_3$,
$E = e_0 + e_1u_1$ (canonical, syzygy-free), gauge as above; 12 unknowns,
20 trilinear equations from $\Delta \equiv 1$.

The three linear equations force the $A$-block:
$a_1 = -\tfrac{3b_1 + d_1}{2}$, $a_2 = -\tfrac{2(b_2+d_2)}{3}$,
$a_3 = -\tfrac{b_3 + 3d_3}{4}$. After substitution, the Keller variety is
**2-dimensional** and decomposes *exactly* (not just up to closure) as

$$
V(I) \;=\; V_1 \,\cup\, V_2 .
$$

- $V_1$ (elementary $z$-shears): $F = (x,\ y,\ z + y^3(e_0 + e_1u_1))$,
  i.e. $z \mapsto z + y^3 h(xy^2)$ — triangular, det $= 1$.
- $V_2$ ($w$-shears): with $w = xz - \nu y$ (a weight-$(-1)$ invariant of
  the flow) and parameters $(\lambda,\nu)$,
  $$
  F = \bigl(x,\ y + \lambda x w^3,\ z + \lambda\nu w^3\bigr),
  \qquad w\circ F = w .
  $$
  In coefficients: $C = \lambda u_3$, $B = 1 - \lambda\nu^3u_1 +
  3\lambda\nu^2u_2 - 3\lambda\nu u_3$, $D = 1 + 3\lambda\nu^3u_1 -
  3\lambda\nu^2u_2 + \lambda\nu u_3$, $E = -\lambda\nu^4$. Since
  $(x,y,z)\mapsto(x, xz-\nu y, z)$ is triangular, $V_2$ is an elementary
  shear in the coordinates $(x, w, z)$ — **tame**.

**Completeness proof** (all in the script, exact): the six relations
$b_3+3d_3$, $c_3b_2-3d_3^2$, $c_3d_2+3d_3^2$, $c_3^2b_1+d_3^3$,
$c_3^2d_1-3d_3^3$, $c_3^3e_0+d_3^4$ have powers *inside* the ideal $I$
(Gröbner reduction certificates), so they vanish on the whole variety;
$c_3e_1^3 \in I$ splits the variety into $\{c_3 \neq 0\}$ (where the
relations pin every coordinate to $V_2$ with $\lambda = c_3$,
$\nu = d_3/c_3$), $\{e_1 \neq 0\}$ (where $e_1\cdot v^k \in I$ forces all
other coordinates to zero: $V_1$), and $\{c_3 = e_1 = 0\}$ (where
$v^k \in I + (c_3,e_1)$ forces $b_\ast = d_\ast = 0$: $V_1$ again).

Both families have explicit polynomial inverses (verified by composition),
generic fiber cardinality 1 (Gröbner elimination as in `jcqft/fibers.py`),
no 2:1/3:1 witness, and survive the infinity prefilter only as members of
its known false-positive class (nonlinear automorphisms,
`jcqft/prefilter.py`).

Independently of the decomposition, all six mechanism systems
(2:1 and 3:1, at $(p,q) = (1,1)$, $(1,0)$, and on the $y=0$ stratum) are
proved **empty** by Gröbner bases $= [1]$. The residual 2-torus acts freely
on witness points $(p_0, q_0)$, $p_0 \neq 0$ (and scales $u_3$ on the
$y = 0$ stratum), preserving the box and the gauge, so these pointwise
checks exhaust all witnesses.

> **Verdict (degree-1 box).** The weight system (2,−1,−3) *does* support
> nonlinear Keller maps — two 2-parameter families — but **every Keller map
> in the box is a tame automorphism**. No counterexample; both orbifold
> mechanisms are provably empty. This is the exact analogue of the
> rigidity statement of `docs/NEW_COUNTEREXAMPLES.md` §3, but here it is a
> *complete classification*, not a first-order/numerical statement.

### 5.2 CE2 box (`--full`): mechanism no-go

Box enlarged in the syzygy-free directions: $C \in
\mathrm{span}\{u_3, u_3^2\}$, $E \in \mathrm{span}\{1, u_1, u_1^2\}$
($A,B,D$ still degree ≤ 1); 14 unknowns, 34 equations. The six mechanism
emptiness queries are again all **EMPTY** (exact, over $\mathbb{Q}$;
longest single Gröbner basis ≈ 5 min). So no 2:1 or 3:1
orbifold-mechanism counterexample exists in the CE2 box either. (A full
variety decomposition of CE2 was not attempted — see limitations.)

### 5.3 Full degree-2 box: attempted with three engines, unresolved

29 unknowns, 57 Keller equations plus witnesses. Status after a serious
tooling campaign (2026-07-23): **unresolved, with the obstruction now
precisely mapped.**

1. **sympy Buchberger** (`--deg2`): no query finishes (>25 min each, also
   mod $p$).
2. **msolve F4** (v0.10.1, submodule, via `jcqft/gb_backend.py`, all runs
   under a 16 GB address-space cap after an uncapped run froze the
   machine): all six queries exceed the cap over $\mathbb{Q}$ — and,
   decisively, **also mod three independent ~30-bit primes**, even after
   the exact reductions below. The memory blow-up is intrinsic to F4 on
   this system, not a coefficient-growth artifact.
3. **Exact structural reductions** (`scripts/deg2_213_elim.py`, all
   asserted): (i) the 57 Keller equations are *triangular* in the A-block
   with integer pivots in $\{2,\dots,7\}$ — the A-block eliminates by a
   global polynomial substitution, no denominators, valid over any field
   of characteristic $0$ or $p > 7$ (49 equations, 21 unknowns, degree
   $\le 4$); (ii) the Rabinowitsch saturation variable is unnecessary:
   $\Delta = \tilde a X + \tilde b Y + \tilde e Z$ identically, so at a
   witness point the Keller condition itself forces the third chart
   function nonzero; (iii) linear witness equations eliminate 1–2 more
   unknowns per query. Net: 19–20 unknowns per query. **Still past the
   16 GB cap for F4, even mod $p$.**
4. **Singular truncated `std()` ladder** (degBound, one-sided: finding
   $1$ below the bound is an exact certificate; memory-frugal, ~2 GB):
   mod $p$, degBound 4 and 5 complete in under a minute — no certificate;
   degBound 6 and 7 exceed 20-minute budgets. The unit certificate (if
   the ideal is unit) has degree $\ge 6$, where the truncated basis is
   already past this machine's practical time budget.

**Honest verdict.** The six degree-2 mechanism queries are beyond direct
Gröbner methods on 30 GB hardware: F4 hits a memory wall (>16 GB even
mod $p$ at 20 unknowns), Buchberger/std a time wall (certificate degree
$\ge 6$). The box remains **unresolved, not ruled out**. Plausible
routes: much larger memory, further stratification (by vanishing leading
block coefficients), or numerical algebraic geometry
(HomotopyContinuation.jl) for discovery with exact certification of any
find. Reproduce: `.venv/bin/python scripts/deg2_213_elim.py`
(reductions + capped screens; `--skip-exact` for the mod-$p$ screen
only).

## 6. Honest limitations

1. **Box-relative.** Rigidity is proved for the degree-1 box (complete)
   and, for the orbifold mechanisms only, the CE2 box. Higher-degree
   ansätze are open — the degree-2 box already exceeds sympy's Gröbner
   engine (no msolve/Singular available in this environment).
2. **Mechanism no-go ≠ injectivity proof (CE2).** In the CE2 box only the
   *stabilizer-jump* (2:1/3:1 orbit) mechanisms are excluded. An
   equivariant Keller map could in principle be non-injective by
   identifying two *distinct free orbits* (non-injectivity of the induced
   map on the quotient cone); ruling that out needs the full variety
   decomposition, done here only for the degree-1 box. (For Alpöge–Mathew
   the realized mechanism *is* the stabilizer jump, which motivates
   prioritizing it.)
3. **Gauge slice.** $C(0) = 0$ uses the constant $xz$-shear and is
   box-preserving, but the higher shear gauges (e.g. $z \to z + y^3h$)
   are *not* box-preserving; consequently the two families $V_1, V_2$
   intersect the gauge group nontrivially rather than being a transversal.
   This affects bookkeeping, not the classification (membership in
   $V_1 \cup V_2$ and tameness are proved directly).
4. The identity $\det DF = \Delta(p,q)$ is chart-based (dense open chart,
   denominators cleared) — as permitted; the intrinsic normal-form
   determinant on the cone is also implemented (`det_df_normal_form`) and
   checked to agree, but the *search* runs on the chart form.

## 7. QFT reading

The $\mathbb{Z}_3$ stratum made this the most promising place to look for a
new global-anomaly class (3:1 covering, cube-root vacuum escape,
$\mathbb{Z}_3$ monodromy). The answer in the searched boxes is a clean
**no-go**: the Keller (unit-determinant) constraint is compatible with the
$(2,-1,-3)$ grading only through shear-type field redefinitions — the
theory admits no non-properness defect here. The contrast with
$(1,-1,-2)$, where the Alpöge–Mathew defect exists already at low degree,
suggests the defect needs the *free* invariant ring (the $(1,-1,-m)$
family) or larger boxes; the $A_1$ cone singularity of the (2,−1,−3)
quotient appears to rigidify the low-degree Keller moduli.
