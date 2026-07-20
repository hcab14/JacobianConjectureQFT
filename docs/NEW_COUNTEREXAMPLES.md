# Reverse-Engineering the Counterexample and the Search for New Ones

*(Results of the ansatz/deformation search. Reproducible via
`.venv/bin/python search_counterexamples.py` (~40 s). Strategy context in
`docs/SEARCH_STRATEGIES.md`. This writeup consolidates a long sub-agent
exploration — which twice exceeded its session limit — with an independent
re-implementation and verification of its key claims.)*

## 1. The construction mechanism, fully reduced

### 1.1 Normal form in the scaling invariants

The map is equivariant under $(x,y,z) \to (\lambda x, \lambda^{-1}y,
\lambda^{-2}z)$. In the invariants $w = xy$, $v = x^2 z$ it takes the normal
form (verified exactly):

$$
F \;=\; \Bigl(\ \frac{P(w,v)}{x^2},\ \ \frac{Q(w,v)}{x},\ \ x\,R(w,v)\ \Bigr),
$$

$$
P = (1+w)\bigl(v(1+w)^2 + w^2(3w+4)\bigr),\qquad
Q = 3v(1+w)^2 + 9w^3 + 12w^2 + w,\qquad
R = 2 - 3w - v .
$$

*Realizability* (that $F$ be polynomial in $(x,y,z)$) is the condition that
every monomial $w^iv^j$ of $P$ has $x$-degree $i + 2j \ge 2$, of $Q$ has
$i + 2j \ge 1$ (and $R$ is unconstrained).

### 1.2 The Keller condition collapses to two dimensions

With $A = P R^2$, $B = QR$ (these are $x^2F_1\cdot(F_3/x)^2$-type invariant
combinations), one computes (verified exactly):

$$
\det DF = \text{const} = \kappa
\quad\Longleftrightarrow\quad
\frac{\partial(A,B)}{\partial(w,v)} \;=\; \kappa\, R^2 .
$$

For the Alpöge–Mathew data, $J_2(PR^2, QR) = -2R^2$ identically. **The 3D
Keller search within this symmetry class is a 2D problem** — this is the
main structural discovery of the search, and the analogue of how symmetry
reduces field equations in physics.

### 1.3 The non-injectivity criterion is one algebraic condition

$F$ is 2:1 on the orbit $x \mapsto (x,\, w_0/x,\, v_0/x^2)$ **iff**
$(w_0, v_0)$ is a common zero of $Q$ and $R$ with $P(w_0,v_0) \neq 0$: the
image is then $(P/x^2,\,0,\,0)$, invariant under $x \to -x$ (the residual
$\mathbb{Z}_2$). For the counterexample: $(w_0,v_0) = (-\tfrac32,\tfrac{13}{2})$,
$P = -\tfrac14$ there — reproducing the triple point.

So a *recipe* for counterexamples of this class: find $(P,Q,R)$ with (i) the
2D Keller identity, (ii) realizability, (iii) a common zero of $(Q,R)$ off
$\{P = 0\}$.

## 2. Gauge group of the ansatz

Transformations preserving the class and mapping counterexamples to
equivalent ones (compositions with polynomial automorphisms):

- **Torus scalings** of source $(x,y,z)$ and target $(a,b,c)$ (6 parameters,
  one redundancy — rank 5 combined with the ansatz normalization).
- **Weighted-triangular source moves**: $y \to y + f(x,z)$ with $f$ of weight
  $-1$ (e.g. $xz,\ x^3z^2$: $\delta w = v, v^2$), $z \to z + g(x,y)$ of
  weight $-2$ (e.g. $y^2,\ xy^3,\ x^2y^4$: $\delta v = w^2, w^3, w^4$).
- **Weighted-triangular target moves**: $a \to a + s\,b^2$, $a \to a + s\,b^3c$,
  $b \to b + s\,ac$ ($\delta P = Q^2,\ Q^3R$, $\delta Q = PR$).

All are unit-Jacobian (or torus) automorphisms; the sub-agent additionally
checked several would-be "obstructions" ($\delta E \propto v^2$, $v^2w$) and
identified them as gauge directions of this type.

## 3. Rigidity of the Alpöge–Mathew solution (evidence)

First-order deformation theory around $(P,Q,R)$ inside the degree box
$\deg_x P \in [2,7]$, $\deg_x Q \in [1,6]$, $\deg_x R \in [0,4]$
(43 coefficient unknowns incl. $\delta\kappa$), from `search_counterexamples.py`:

- Kernel of the linearized Keller condition: **dimension 15**.
- In-box gauge tangents: **rank 9**, all inside the kernel (4 further gauge
  generators have tangents leaving the box and could not be counted).
- The 6 kernel directions orthogonal to the in-box gauge orbit were pushed
  through **nonlinear continuation** (Gauss–Newton with displacement
  constraint on the full 2D Keller system): every one either fails to
  converge or jumps to a distant solution ($|u - u_0| > 15$ for step $0.05$
  — a far-away gauge image, not a branch through the base point).
- **Conclusion: 0 genuine nearby families.** Within this ansatz box the
  Alpöge–Mathew solution is *locally rigid modulo gauge* — it does not sit
  inside a continuous family of inequivalent counterexamples.

The independent sub-agent run, with a differently-shaped moduli box, found
kernel dimension $=$ gauge dimension $= 11$ exactly, i.e. rigidity already at
the linear level in that box — consistent with the above.

Caveats: numerical continuation (float Gauss–Newton, finite-difference
Jacobians) gives strong evidence, not proof, of higher-order obstruction; and
rigidity is relative to the degree box and to this weight system.

## 4. Dead ends (recorded so they are not retried)

- **$\kappa = 0$ collapse.** A promising-looking solution branch of the 2D
  Keller equations with $P_2 \propto R_1^2$, $R_1 = \alpha + \beta R_0^2$
  survives all algebraic constraints but, upon imposing realizability
  ($P$ of $x$-degree $\ge 2$ at the origin), forces two rows of $DF(0)$ to
  degenerate — i.e. $\det DF \equiv 0$. Not a Keller map; the branch is empty.
- **Naive coefficient scans** (varying the integer coefficients of $P, Q, R$
  without the gauge/deformation framework) reproduce only gauge images of
  the original or non-Keller maps.

## 5. Honest conclusion and next steps

**No genuinely new counterexample was found.** The search instead produced
a strong structural statement: within its symmetry class and search box, the
Alpöge–Mathew map is an isolated (rigid) solution modulo the gauge group —
consistent with it being a distinguished algebraic object rather than a
point of a continuous moduli space.

Next steps, in order of expected value:

1. **Other weight systems** (`docs/SEARCH_STRATEGIES.md` §1.1): the
   $(1,-1,-2)$ grading is one choice; systems like $(1,-1,-3)$,
   $(2,-1,-3)$, or 4-field gradings with $\mathbb{Z}_3$ residual symmetry
   (3:1 orbits, cube-root escape, possible $\mathbb{Z}_3$ monodromy) are
   untouched and each gives a *finite* 2D-reduced problem like §1.2.
2. **Larger degree boxes** in the present weight system (rule out
   higher-degree relatives; the 2D reduction keeps this cheap).
3. **The bootstrap route** (`SEARCH_STRATEGIES.md` §1.2): prescribe the
   degree-3 eliminant data directly and solve the forced consistency
   conditions ($\operatorname{disc} = -(\text{square})\cdot p$, integrality,
   uniruled $\{p=0\}$).
4. **Exact certification**: replace the numerical continuation by Groebner
   computation of the solution ideal of the 2D Keller system in the box,
   to turn "rigid (numerical evidence)" into a theorem.
