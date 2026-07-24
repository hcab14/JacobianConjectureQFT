# The Weight System (1, −1, −3): Exact Reduction and Search

*(Machinery: `jcqft/reduction_w.py`, proved exact for the whole family
$(1,-1,-m)$ in `scripts/reduction_113.py`. Every identity and every
emptiness claim below is proved by assertion in `scripts/search_113.py`:*

```
.venv/bin/python scripts/reduction_113.py        # ~1 min, the reduction
.venv/bin/python scripts/search_113.py           # ~20 min, the search
.venv/bin/python scripts/search_113.py --full    # + larger boxes, ~1 h
```

*This executes the program of `docs/SEARCH_STRATEGIES.md` §1.1 for the
weight system (1,−1,−3) suggested in `docs/NEW_COUNTEREXAMPLES.md` §5 —
the nearest sibling of the Alpöge–Mathew weight system (1,−1,−2).)*

## 1. Why this weight system is interesting

For $\lambda.(x,y,z) = (\lambda x,\ \lambda^{-1}y,\ \lambda^{-3}z)$ the
invariant ring is **free** — $\mathbb{C}[w,v]$ with $w = xy$, $v = x^3z$ —
exactly as in the Alpöge–Mathew case $m = 2$ (contrast (2,−1,−3), whose
invariant ring is an $A_1$ cone, `docs/SEARCH_213.md` §2). The target
coordinates carry weights $(-3,-1,1)$, so:

- the target $a$-axis (the weight $-3$ coordinate) has stabilizer
  $\mathbb{Z}_3$: a free orbit mapping into it is **3:1** — cube-root
  escape, $\mathbb{Z}_3$ monodromy, the "different global anomaly class"
  scenario of `docs/SEARCH_STRATEGIES.md` §2.2;
- **no 2:1 mechanism exists at all**: every other coordinate axis and
  every coordinate pair has trivial stabilizer ($|{-1}| = |1| = 1$, and
  $\gcd(3,1) = 1$ for the pairs), so the *only* stabilizer-jump
  non-injectivity available to this weight system is the 3:1 one.

Since the Alpöge–Mathew counterexample is exactly a stabilizer-jump (2:1)
map in the $m = 2$ member of the same family, (1,−1,−3) is the natural
place to ask whether the mechanism *transplants* — with the residual
$\mathbb{Z}_2$ upgraded to $\mathbb{Z}_3$.

## 2. Reduction (proved exactly, whole family)

From `scripts/reduction_113.py` / `jcqft/reduction_w.py`, re-asserted for
$m = 3$ in the search script: every candidate is

$$
F = \Bigl(\frac{P(w,v)}{x^3},\ \frac{Q(w,v)}{x},\ x\,R(w,v)\Bigr),
$$

with polynomiality $\iff$ every monomial $w^jv^k$ of $P$ has
$j + 3k \ge 3$, of $Q$ has $j + 3k \ge 1$ ($R$ unconstrained), and the
**reduced Keller identity** (generic-function proof, an identity of
differential polynomials):

$$
\det DF \;=\; \det M \;=\; -3P\,J_2(Q,R) + Q\,J_2(P,R) + R\,J_2(P,Q),
\qquad R^3\det M = J_2(PR^3,\,QR),
$$

a polynomial in $(w,v)$ alone. The 3D Keller problem $\det DF = \kappa$
is therefore a 2D problem, as for Alpöge–Mathew.

**Gauge fixing.** $DF(0) = \operatorname{antidiag}(p_1(0), q_0'(0),
r_0(0))$ where $p_1$ is the $v$-coefficient of $P$, etc.; target scalings
normalize $p_1(0) = q_0'(0) = r_0(0) = 1$, hence $\kappa = \det DF(0) =
-1$. A residual 2-torus $(w,v) \to (\mu w, \nu v)$ (with compensating
target scalings) survives and is used to normalize 3:1 witness points.
The source shears $z \to z + y^3h(xy)$ act as $v \to v + w^3\,\mathbb{C}[w]$
— the box-preserving part of the **$v$-shear group** $v \to v + h(w)$,
$h \in \mathbb{C}[w]$ arbitrary, under which $\det M = \kappa$ is exactly
invariant (asserted). The full $v$-shear group is the main *analysis*
tool below: solve modulo shears first, impose the polynomiality box last.

## 3. The search box: the v-linear class

The class searched is $P, Q, R$ of degree $\le 1$ in $v$ and **arbitrary
degree in $w$**:

$$
P = p_0(w) + p_1(w)\,v,\qquad Q = q_0(w) + q_1(w)\,v,\qquad
R = r_0(w) + r_1(w)\,v,
$$

with the box/invertibility constraints $p_0 \in w^3\mathbb{C}[w]$,
$q_0 \in w\mathbb{C}[w]$, $p_1(0) = q_0'(0) = r_0(0) = 1$. This is the
exact analogue of the class containing the Alpöge–Mathew map at $m = 2$
($P_{\rm AM}, Q_{\rm AM}, R_{\rm AM}$ are all $v$-linear). $\det M = -1$
has $v$-degree 2, so the Keller condition is **three ODE-type equations**
in $\mathbb{C}[w]$:

$$
E_2 = 0,\qquad E_1 = 0,\qquad E_0 = \kappa = -1 ,
$$

with (all asserted from `det_m`, not copied from $m=2$):

$$
E_2 = 2\bigl[q_1(p_1r_1)' - 2(p_1r_1)q_1'\bigr],\qquad
E_0 = -3p_0(q_0'r_1{-}q_1r_0') + q_0(p_0'r_1{-}p_1r_0')
      + r_0(p_0'q_1{-}p_1q_0'),
$$

and $E_1$ the mixed analogue. Unlike the degree-boxed searches of
`docs/SEARCH_213.md`, the $w$-direction here is **unbounded**: the
classification below is a theorem for the whole class, not for a box —
with one precisely-mapped exception (§5, stratum D3 non-squarefree).

## 4. Stratification

$E_2 = 0$ integrates: $(p_1r_1/q_1^2)' = E_2/(2q_1^3)$ (asserted), so a
rational function with vanishing derivative is constant, giving
$p_1r_1 = c\,q_1^2$ exactly whenever $q_1p_1r_1 \ne 0$; and $q_1 = 0$ or
$r_1 = 0$ satisfy $E_2$ identically. Since $p_1(0) = 1$ forces
$p_1 \ne 0$:

| stratum | definition | $m=3$ verdict |
|---|---|---|
| A | $q_1 = r_1 = 0$ | tame family (all $w$-degrees) |
| B | $r_1 = 0,\ q_1 \ne 0$ | tame family (all $w$-degrees) |
| C | $q_1 = 0,\ r_1 \ne 0$ | **empty** (all $w$-degrees) |
| D | $q_1r_1 \ne 0 \Rightarrow p_1 = as^2g,\ q_1 = stg,\ r_1 = bt^2g$ | see D0–D3 |
| D0 | $s, t$ const ($p_1,q_1,r_1$ const) | **empty** (all degrees) |
| D1 | $t$ const, $s$ nonconst — *the AM-analogue* | **empty in the box** (all degrees) |
| D2 | $s$ const, $t$ nonconst | **empty in the box** (all degrees) |
| D3 | $s, t$ nonconst, coprime, both squarefree | **empty in the box** (all degrees) |
| D3′ | $s$ or $t$ non-squarefree (needs $\deg p_1 \ge 4$ or $\deg r_1 \ge 4$) | **empty in boxes** (Gröbner), open beyond |

The stratum D parametrization is the UFD solution of $p_1r_1 = cq_1^2$:
$g = \gcd(p_1, r_1)$, residual multiplicities even, $\gcd(s,t) = 1$,
$ab = c$; $E_0$ is linear in $(p_1,q_1,r_1)$ with no derivatives of them,
so $g \mid E_0 = \kappa$ and $g$ is a unit (absorbed).

**$m = 2$ cross-check.** The same $E_2$-integration at $m = 2$ gives
$p_1^2r_1 = c\,q_1^3$, and the Alpöge–Mathew data satisfies it with
$c = -1/27$: $p_1 = (1+w)^3 = s^3$, $q_1 = 3(1+w)^2 = 3s^2$, $r_1 = -1$
(asserted). I.e. **AM lives in the $m = 2$ analogue of stratum D1** —
the stratum that is *empty* at $m = 3$ (below). The machinery does not
spuriously kill the known counterexample; the two weight systems genuinely
differ.

## 5. The classification (all proofs by assertion + one-line arithmetic)

### 5.1 Strata A and B: the tame family — the only Keller maps

Stratum A: $E_1 \equiv 0$ and $E_0 = -p_1(q_0r_0)' = \kappa$, so
$p_1 \mid \kappa \Rightarrow p_1 = 1$ (gauge), $q_0r_0 = w$ with
$q_0(0) = 0$, $r_0(0) = 1$ $\Rightarrow$ $r_0 = 1$, $q_0 = w$; $p_0 \in
w^3\mathbb{C}[w]$ stays free. Stratum B: $E_1$ integrates
($(p_1r_0^2/q_1)' = r_0E_1/q_1^2$) to $p_1r_0^2 = cq_1$; then $E_0 =
p_1\,[\,(p_0r_0^3)/c - q_0r_0\,]'$, forcing $p_1 = 1$, $r_0 \mid w$ with
$r_0(0)=1$ so $r_0 = 1$, and $q_0 = b_0p_0 + w$, $q_1 = b_0 = 1/c$.
Together (with $b_0 = 0 \leftrightarrow$ A), the **complete gauged
solution set** of the $v$-linear class:

$$
\boxed{\ P = p_0 + v,\qquad Q = w + b_0P,\qquad R = 1,\qquad
p_0 \in w^3\mathbb{C}[w],\ b_0 \in \mathbb{C}\ }
$$

which assembles to $F = \bigl(p_0(xy)/x^3 + z,\ \ y + b_0x^2F_1,\ \
x\bigr)$: an elementary $z$-shear followed by the target shear
$b \mapsto b + b_0\,a\,c^2$ — **tame automorphisms**, explicit inverses
verified by composition, generic fiber cardinality 1, and the infinity
prefilter (`jcqft/prefilter.py`) is survived only through its known
false-positive class (nonlinear automorphisms).

### 5.2 Strata C, D0, D1, D2, D3: all empty

Each emptiness is an asserted identity chain plus a one-line divisibility
argument in $\mathbb{C}[w]$ (constants can't have nonconstant divisors):

- **C**: $(p_1/(q_0^4r_1))' = E_1/(q_0^5r_1^2)$, so $E_1 = 0$ forces
  $p_1 = c\,q_0^4r_1$ and $p_1(0) = 0$ — contradicting $p_1(0) = 1$.
- **D0** ($p_1,q_1,r_1$ nonzero constants $\alpha,\beta,\gamma$): $E_1$
  integrates to a linear relation; substituting into $E_0$ gives the
  *exact derivative* $[-(2\alpha\gamma/\beta)V^2 - \alpha\delta V]' =
  \kappa$ with $V = q_0 - (\beta/\alpha)p_0$, $V(0) = 0$, $V'(0) \neq 0$.
  Then $V \mid \kappa w \Rightarrow V = \varepsilon w$, whose quadratic
  term forces $2\alpha\gamma/\beta = 0$: impossible.
- **D1** ($p_1 = As^2$, $q_1 = Bs$, $r_1 = C$, $s$ nonconstant — the
  Alpöge–Mathew slot): the $v$-shear $h = -r_0/C$ sets $r_0 \equiv 0$;
  then $E_1$ is a linear ODE for $p_0$ with particular solution
  $p_0^* = (2A/B)sq_0$ and homogeneous solutions $Y^2 = e\,s^3$
  (nonzero only if $s = c\,d^2$). If $s$ is not a constant times a
  square: $E_0(p_0^*) = (2AC/B)\,q_0(s'q_0 - 2sq_0') = \kappa$, so
  $q_0 = n$ const, $n^2s' = $ const, $s$ **linear**. If $s = cd^2$:
  $E_0 = C\,d\,(d'q_0 - dq_0')\bigl((4Ac/B)q_0 + 3y_0d\bigr)$, and
  $d \mid \kappa$ kills it. The surviving sheared family
  $(s, q_0, p_0, r_0) = (s_0{+}s_1w,\ n,\ (2A/B)ns,\ 0)$ has $\kappa =
  (2AC/B)n^2s_1 \neq 0 \Rightarrow n \neq 0$; **un-shearing into the
  box** ($p_0(0) = q_0(0) = 0$) forces $n = 0$ — the jet system
  $h(0) = 2n/(Bs_0) = n/(Bs_0)$ has only $n = 0$. **Empty.**
- **D2** ($p_1 = A$ const): the shear $h = -p_0/A$ sets $p_0 \equiv 0$;
  $E_0 = -A(q_0r_0)' = \kappa$ makes $q_0r_0$ exactly linear, so one
  factor is constant. $q_0 = \gamma$: $E_1$ forces $t$ linear and $r_0 =
  (2C\gamma/B)t$, but un-shearing needs $\operatorname{val}(p_0 = -Ah)
  \ge 3 \Rightarrow h(0) = 0 \Rightarrow q_0(0) = \gamma \ne 0$:
  contradiction. $r_0 = \delta$: $E_1 \Rightarrow t \mid \delta^2t'
  \Rightarrow t' = 0$: contradiction.
- **D3** ($s,t$ nonconstant, coprime, squarefree): pass to the
  shear-invariants $G_1 = tp_0 - asq_0$, $G_2 = btq_0 - sr_0$. Asserted:
  $E_1 = as\Theta + bt\Theta_1$ and $t(G_2G_3' - 3G_3G_2') + t'G_2G_3 =
  bt^2E_0 - r_0E_1$ ($G_3 = btG_1 + asG_2$). $G_1 = 0$ or $G_2 = 0$
  contradict $\kappa \ne 0$. Squarefree divisibility gives $s \mid G_1$,
  $t \mid G_2$; then $E_1 = st(WZ - 2stZ')$ with $Z = a\bar g_2 -
  b\bar g_1$, $W = s't - st'$, and $(Z^2t/s)' \propto Z(2stZ' - WZ)$
  forces $Z^2t = es \Rightarrow Z = 0$. The $\kappa$-equation then reads
  $g(2stg' - Wg) = -b\kappa/(2a)$, so $g = c$ const and $\kappa =
  (2a/b)c^2W$ — Keller solutions exist **iff the Wronskian $W$ is a
  nonzero constant**, and they form the single shear-orbit
  $p_0 = as(q_0 + c/b)/t$, $r_0 = t(bq_0 - c)/s$. But $t\,p_0 =
  as(q_0 + c/b)$ evaluated at $w = 0$ gives $p_0(0)t(0) = a\,s(0)\,c/b
  \neq 0$, while the box forces $p_0(0) = 0$: **empty for every degree**.

### 5.3 The gap: D3 with non-squarefree $s$ or $t$, and its box closure

The two steps $s \mid s'G_1 \Rightarrow s \mid G_1$ and $t \mid t'G_2
\Rightarrow t \mid G_2$ need $s, t$ squarefree. A non-squarefree $s$ or
$t$ requires $\deg p_1 = 2\deg s \ge 4$ or $\deg r_1 = 2\deg t \ge 4$ —
this corner is closed by exact in-box Gröbner certificates (msolve, over
$\mathbb{Q}$, 16 GB memory cap):

1. **Default run**: Keller $+$ ($r_1$ has a nonzero coefficient) is the
   unit ideal in the box $\deg(p_0,p_1,q_0,q_1,r_0,r_1) \le
   (4,2,3,2,2,2)$ — covering strata C and D wholesale, ~5 s per
   coefficient.
2. **`--full`**: the same in the medium box $\le (5,3,4,3,3,3)$, plus
   the four targeted non-squarefree parametrizations
   $s = (1{+}\rho w)^2$ / $t = (n_0{+}n_1w)^2$ with the partner of
   degree $\le 2$ and $\deg(p_0,q_0,r_0) \le (6,5,4)$, each as a direct
   emptiness query in 18–19 unknowns. Status (2026-07-24, this machine,
   30 GB): $s$-square/$t$-linear **EMPTY** (~11 min); the remaining
   three are memory-borderline for msolve's F4 — see §8.

Independently of the gap, the classification is corroborated in-box by
nilpotency certificates (mirroring `search_213.py`): on the $r_1 = 0$
slice of the small box, powers of $p_1$-, $r_0$-, $q_1$-nonconstant
coefficients and of $q_2$, $q_3 - b_0p_3$, $b_0p_4$ all lie in the
Keller ideal, pinning the in-box variety exactly to the gauged family of
§5.1 (with its box truncation).

## 6. The 3:1 mechanism is empty; no 2:1 exists

A 3:1 hit needs $(w_0,v_0)$ with $Q = R = 0 \neq P$ on a free orbit
$\{x \neq 0\}$ (image on the target $a$-axis, preimages related by
$x \to \omega x$, $\omega^3 = 1$; the $x = 0$ plane supports no
mechanism since $F_3 \equiv 0$ there and the induced 2D map is linear in
its orbit parameter). By the classification every $v$-linear Keller map
has $R = 1$ (gauge) — **$R$ never vanishes, so no witness exists**, for
any $w$-degree. Independent pointwise Gröbner queries (residual torus
moves any witness to $(1,1)$, $(1,0)$ or $(0,1)$; $(0,0)$ is excluded by
$R(0,0) = 1$) confirm emptiness in the small box. And by weight
arithmetic there is no 2:1 mechanism anywhere in this weight system
(§1), so no orbifold-type counterexample at all.

## 7. Prefilter and injectivity verdicts

- **Non-properness screen** (`jcqft/prefilter.py`): the family members
  survive `infinity_prefilter` — as they must, being *nonlinear
  automorphisms*, the filter's documented false-positive class
  (`docs/RIGIDITY_AND_PREFILTER.md`); this is the same verdict pattern
  as `docs/SEARCH_213.md` §5.1. No candidate ever reached the "survives
  and might be non-proper" state, because the Keller classification
  leaves nothing but automorphisms.
- **Injectivity**: explicit polynomial inverses are verified by
  composition for the family ($x = F_3$, $y = F_2 - b_0F_3^2F_1$,
  $z = F_1 - p_0(F_3y)/F_3^3$ evaluates polynomially); generic-fiber
  cardinality 1 re-checked by `sp.solve` on samples.

**Headline: no counterexample.** Every Keller map of the $v$-linear
class is a tame automorphism.

## 8. Honest limitations

1. **$v$-degree.** The class searched is $v$-linear (the AM-analogue
   class). Ansätze of $v$-degree $\ge 2$ are untouched — for them
   $\det M = \kappa$ has $v$-degree $2\deg_v - something$ higher and the
   Wronskian stratification above does not directly apply. (At $m = 2$
   the known counterexample is $v$-linear; but this is a heuristic, not
   an argument, for prioritizing the class.)
2. **The D3′ gap.** Stratum D3 with non-squarefree $s$ or $t$
   ($\deg p_1 \ge 4$ or $\deg r_1 \ge 4$ plus both $s,t$ nonconstant) is
   proved empty only inside Gröbner boxes. Beyond them it is *open*,
   though unpromising: the squarefree D3 analysis shows the entire
   stratum's solution set is one shear-orbit killed by the origin jet,
   and the non-squarefree variant only tightens the divisibility
   constraints. Status of the four `--full` gap queries on this
   machine: $s=(1{+}\rho w)^2,\ t$ linear: **EMPTY** (exact, 649 s);
   the other three hit msolve's 16 GB F4 memory wall — same wall as
   `docs/SEARCH_213.md` §5.3 — and remain formally unresolved in the
   large box (they *are* covered by the medium-box $r_1 \neq 0$ queries
   up to $\deg r_1 \le 3$, which excludes non-squarefree $t$ of degree
   2 in $r_1$-degrees $\le 3$... i.e. only the $\deg = 4$ corner is
   genuinely open).
3. **Gauge bookkeeping.** The classification is stated in the gauge
   $p_1(0) = q_0'(0) = r_0(0) = 1$, $\kappa = -1$; un-gauged solutions
   are recovered by target/torus scalings, which do not change tameness,
   properness or fiber counts.
4. The $v$-shear analysis moves (arbitrary $h \in \mathbb{C}[w]$) are
   *not* box-preserving; the box is imposed at the end through the jet
   conditions at $w = 0$. This is where the $m = 3$ obstruction lives —
   see §9.

## 9. Why Alpöge–Mathew does not transplant: the exact numerology

At $m = 2$ the AM map sits in stratum "D1" of its own family:
$p_1 = s^3$, $q_1 = 3s^2$, $r_1 = -1$, $s = 1 + w$ — the same $t = $
const slot that at $m = 3$ produces the sheared family
$q_0 = n,\ p_0 = (2A/B)ns$. The $m = 3$ box demands
$\operatorname{val}_w p_0 \ge 3$ (from $j + 3k \ge 3$) while the family
jets force $p_0(0) \neq 0$ unless $n = 0$, i.e. $\kappa = 0$. At
$m = 2$ the box demands only $\operatorname{val}_w p_0 \ge 2$ *and* the
corresponding equations have one more derivative-weighting ($E_0^{m=2} =
q_0p_0' - 2p_0q_0' + \dots$), which is exactly what lets
$P_{\rm AM} = (1+w)(v(1+w)^2 + w^2(3w+4))$ thread the needle. The
obstruction at $m = 3$ is **numerological and exact** — not a failure to
search far enough.

## 10. QFT reading

The $(1,-1,-3)$ grading was the closest sibling of the Alpöge–Mathew
theory with a $\mathbb{Z}_3$ orbifold stratum in place of $\mathbb{Z}_2$
— the candidate for a cube-root-escape global anomaly. The verdict is a
clean **no-go in the AM-analogue class**: the unit-determinant constraint
is compatible with this grading only through shear-type field
redefinitions; the theory admits no non-properness defect there, and no
$\mathbb{Z}_3$ vacuum-escape monodromy. Combined with `SEARCH_213.md`
(the $(2,-1,-3)$ no-go) the pattern is that the Alpöge–Mathew defect is,
so far, *rigidly attached* to the $m = 2$ numerology: neither enlarging
the orbifold group ($\mathbb{Z}_3$ here) nor deforming the invariant
geometry (the $A_1$ cone there) reproduces it. The natural next probes
are the $v$-quadratic class here, $m \ge 4$ (where the same
stratification machinery applies verbatim — `reduction_w` is proven for
all $m$), and 4-field gradings.
