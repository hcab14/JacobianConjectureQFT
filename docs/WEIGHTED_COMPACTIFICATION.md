# Weighted Chart Compactification: Vacua at Infinity as Boundary Points

*(2026-08-05. All exact claims verified in `scripts/weighted_compactification.py`
(~few seconds, exact sympy). Inputs: C\* weights and normal form
(`jcqft/core.py`, `jcqft/reduction.py`); escape curve and non-integrality of
$x$ (`docs/MISSING_OBSERVABLES.md`); C\*-reduced wall in $(u,w)$
(`docs/POSITIVE_GEOMETRY.md`, `docs/WALL_COMPLEMENT.md`); invariants I3/I7
(`docs/CLASSICAL_MAP_INVARIANTS.md`); escape geometry
(`scripts/branch_locus.py`).)*

**Summary.**

1. **Choice of compactification.** Because the Alpöge–Mathew $\mathbb{C}^*$
   weights on fields $(1,-1,-2)$ are *mixed-sign*, ordinary weighted
   projective space $\mathbb{P}(w_1,w_2,w_3)$ (all weights same sign) is the
   wrong object. The minimal structure that turns the known escape curves
   into ordinary boundary points is a **two-chart partial compactification**
   $\overline{X} = U_0 \cup U_\infty$ of field space, with Cartier divisor
   $D_\infty = \{s=0\}$ the "infinity in the positive-weight direction $x$".
   Source space stays affine $\mathbb{A}^3_{(a,b,c)}$: the extended map lands
   on the ordinary Jelonek hypersurface $\{p=0\}$. (The C\*-reduced wall
   already has its own projective closure in $\mathbb{P}^2_{(u:w:1)}$,
   `docs/POSITIVE_GEOMETRY.md`.)
2. **Exact:** the AM escape curve extends to a regular map from a
   $\mathbb{P}^1$-chart into $\overline{X}$, landing on $D_\infty$; the
   classical map extends to a **polynomial** morphism
   $\overline{F}\colon U_\infty\to\mathbb{A}^3$; and
   $\overline{F}(D_\infty)=\{p=0\}$ set-theoretically (Jelonek = image of
   the boundary divisor).
3. **Exact negative:** ordinary $\mathbb{P}^3$ homogenization fails — the
   components of $F$ have total degrees $(7,6,4)$, so there is no morphism
   $\mathbb{P}^3\to\mathbb{P}^3$ by single-variable homogenization; the
   total-degree leading cone mixes $\mathbb{C}^*$-weights and marks fake
   escape directions (e.g. the ray $t\mapsto(t,0,0)$).
4. **Interpretive** (flagged): Gribov-at-infinity, second-type Landau, and
   lattice "solutions from infinity" are ordinary boundary evaluations of
   $\overline{F}$ on $D_\infty$.
5. **Limits:** no continuum $D\ge 1$ claim.

---

## 1. Weights and why not $\mathbb{P}(1,1,2)$

Field $\varphi=(x,y,z)$ and source $J=(a,b,c)$ carry the Alpöge–Mathew
grading (`jcqft/core.py`)

$$
\lambda\cdot(x,y,z)=(\lambda x,\lambda^{-1}y,\lambda^{-2}z),\qquad
\lambda\cdot(a,b,c)=(\lambda^{-2}a,\lambda^{-1}b,\lambda\,c),
$$

with $F(\lambda\cdot\varphi)=\lambda\cdot F(\varphi)$. Equivalently
(`jcqft/reduction.py`): $F=(P(w,v)/x^2,\,Q(w,v)/x,\,x\,R(w,v))$ in the
weight-zero invariants $w=xy$, $v=x^2z$.

Standard weighted projective space requires all weights the same sign.
Here the field weights $(1,-1,-2)$ are mixed: the $\mathbb{C}^*$-action is
a hyperbolic torus action, not a positive GIT linearization on
$\mathbb{A}^3\setminus\{0\}$. The correct coarse geometric object is a
**toric / charted partial compactification along the escaping ray**
(positive weight on $x$), not $\mathbb{P}(1,1,2)$ or $\mathbb{P}^3$.

---

## 2. The escape-chart compactification $\overline{X}$

### 2.1 Charts

- **Affine chart** $U_0=\operatorname{Spec}\mathbb{C}[x,y,z]$ — ordinary
  field space.
- **Escape chart** $U_\infty=\operatorname{Spec}\mathbb{C}[s,y,\gamma]$ with
  transition on $\{x\neq 0\}=\{s\neq 0\}$
  $$
  s=\frac1x,\qquad
  \gamma \;=\; 2x-3x^2 y-x^3 z \;=\; F_3(\varphi),\qquad
  z \;=\; 2s^2-3ys-\gamma s^3.
  $$
  Invertible: $x=1/s$, and $\gamma$ is the third component of $F$ itself.

Glue $U_0\cup U_\infty$ along the transition to get the **escape-chart
compactification** $\overline{X}$. The **boundary divisor** is the Cartier
divisor

$$
D_\infty \;=\; \{s=0\}\subset U_\infty \;\cong\; \mathbb{A}^2_{(y,\gamma)}.
$$

**Weights in the escape chart** (asserted): $(s,y,\gamma)$ carry
$\mathbb{C}^*$-weights $(-1,-1,+1)$. On $D_\infty$ the residual action is
$\lambda\cdot(y,\gamma)=(\lambda^{-1}y,\lambda\gamma)$ with invariant
$y\gamma$.

**Honest scope.** $\overline{X}$ only compactifies escape in the $x$-direction
(the only escaping coordinate: `scripts/branch_locus.py`,
`docs/MISSING_OBSERVABLES.md`). Directions with $|y|\to\infty$ or
$|z|\to\infty$ at bounded $x$ are not added; they are not needed for the
Jelonek set of this map.

### 2.2 Source space

Leave sources as $\mathbb{A}^3_{(a,b,c)}$. No source-boundary is required for
escape: $\overline{F}(D_\infty)$ lands in affine space on $\{p=0\}$.

Compatible reduced compactification (already in the repo): under source
weights $(-2,-1,1)$ the wall descends to the plane cubic $P_2(u,w)=0$ in
invariants $(u,w)=(ac^2,bc)$, whose closure in $\mathbb{P}^2$ is a cuspidal
cubic with flex at infinity (`docs/POSITIVE_GEOMETRY.md`). That is the
toric/projective face of the *source* wall, not a replacement for
$\overline{X}$.

---

## 3. Exact statements

*(Every item below is asserted in `scripts/weighted_compactification.py`.)*

### Theorem A — Escape curve extends to $D_\infty$

For parameters $(y_0,c_3)\in\mathbb{A}^2$, the AM escape curve
(`docs/MISSING_OBSERVABLES.md` §2)

$$
\phi(T)=\Bigl(T,\; y_0,\; \frac{2T-3T^2 y_0-c_3}{T^3}\Bigr)
\qquad(T\neq 0)
$$

is, in escape-chart coordinates, the regular map from the
$\mathbb{A}^1_s$-chart of $\mathbb{P}^1$ (with $s=1/T$)

$$
\widehat\phi(s) \;=\; \bigl(s,\; y_0,\; c_3\bigr)\in U_\infty.
$$

It extends across $s=0$ to the boundary point
$(0,\,y_0,\,c_3)\in D_\infty$. Equivalently: a punctured disk
$0<|s|<\varepsilon$ maps into $U_0\cap U_\infty$, and $s=0$ is an ordinary
point of $D_\infty$.

### Theorem B — $\overline{F}$ is polynomial on $U_\infty$; boundary hits the wall

On $U_\infty$, writing $\varphi=(1/s,\,y,\,2s^2-3ys-\gamma s^3)$,

$$
\begin{aligned}
\overline{F}_1(s,y,\gamma)
  &= -(s+y)\,(\gamma s^2+2\gamma s y+\gamma y^2-2s-y),\\
\overline{F}_2(s,y,\gamma)
  &= -3\gamma s^2-6\gamma s y-3\gamma y^2+6s+4y,\\
\overline{F}_3(s,y,\gamma)
  &= \gamma
\end{aligned}
$$

are **polynomials** in $(s,y,\gamma)$ (no poles on $D_\infty$). Restriction
to the boundary:

$$
\overline{F}\big|_{D_\infty}(y,\gamma)
  \;=\;
  \bigl(\,y^2(1-\gamma y),\;\; y(4-3\gamma y),\;\; \gamma\,\bigr).
$$

This is exactly the $(y_0,c_3)$-parametrization of the Jelonek set from
`docs/MISSING_OBSERVABLES.md` §2, and

$$
p\bigl(\overline{F}(0,y,\gamma)\bigr)\;=\;0
\quad\text{identically in }(y,\gamma).
$$

### Theorem C — Jelonek set is the image of the boundary

The morphism $\overline{F}|_{D_\infty}\colon D_\infty\to\mathbb{A}^3$ has
image equal to the non-properness hypersurface $\{p=0\}$:

- $p\circ\overline{F}|_{D_\infty}\equiv 0$ (Theorem B);
- every point of the standard rational parametrization of $\{p=0\}$ is hit;
- eliminating $(y,\gamma)$ from
  $a=\overline{F}_1|_{s=0}$, $b=\overline{F}_2|_{s=0}$,
  $c=\overline{F}_3|_{s=0}$ recovers the ideal $(p)$ up to the expected
  leading-coefficient factor $c$ of the eliminant (asserted:
  $\operatorname{Resultant}_y=c\cdot p$).

**Compactified-fiber reading (exact + numerical cross-check).** Off the
wall, $F\colon\mathbb{A}^3\to\mathbb{A}^3$ is a degree-$3$ étale cover. At a
generic wall point the affine fiber has one finite point; the escaping pair
approaches **the same** point of $D_\infty$ (both sheets $\to$ one boundary
point with $s\to 0$, $y\to y_*$, $\gamma\to c$). At the empty-fiber cusp
orbit all three sheets sit on $D_\infty$
(`docs/POSITIVE_GEOMETRY.md` §2). Thus $\{p=0\}$ is simultaneously:

- the Jelonek non-properness divisor (I3),
- the image $\overline{F}(D_\infty)$,
- the locus where the compactified cover draws preimages from $D_\infty$.

### Theorem D — Ordinary $\mathbb{P}^3$ homogenization fails / obscures C\*

Concrete negatives (all asserted):

1. **Unequal total degrees.** The components of $F$ have total degrees
   $(7,6,4)$. Single-variable homogenization cannot define a morphism
   $\mathbb{P}^3\to\mathbb{P}^3$.
2. **Leading cone mixes weights.** The total-degree leading forms are
   $(x^3 y^3 z,\; 3x^3 y^2 z,\; -x^3 z)$; their common zero locus is the
   cone $\{x^3 z=0\}$, far larger than the actual escape set. In
   particular the ray $\varphi(t)=(t,0,0)$ lies in that cone, yet
   $F(t,0,0)=(0,0,2t)\to\infty$ — not an escape to a finite source.
3. **C\* grading discarded.** Each $F_i$ is weighted-homogeneous of weight
   equal to the target weight (entire component $=$ weighted leading part),
   but total-degree truncation of $F_3=2x-3x^2 y-x^3 z$ keeps only
   $-x^3 z$ and **throws away** $2x-3x^2 y$, which have the *same*
   $\mathbb{C}^*$-weight $+1$. Ordinary homogenization erases the orbifold
   mechanism (I7) that organizes escape.

The weighted / charted compactification of §2 is exactly the structure that
keeps those same-weight terms and makes escape regular.

---

## 4. Interpretive layer (flagged, not theorems)

*Downstream slogans of I3+I7 (`docs/CLASSICAL_MAP_INVARIANTS.md` §2.2);
not independent claims.*

- **Gribov copies without a horizon.** Finite Gribov-type multiplicity
  ($N=3$ chamber) coexists with $\det DF\equiv -2\neq 0$ everywhere on
  $U_0$. The "missing" copies at the wall are ordinary points of
  $D_\infty$, not zeros of a Faddeev–Popov determinant.
- **Second-type Landau laboratory.** Sheet loss at $\{p=0\}$ is
  pinching against the compactification divisor $D_\infty$, matching the
  amplitudes-program reading of non-properness as a second-type
  singularity (`docs/AMPLITUDES_CONNECTION.md`).
- **Lattice "solutions from infinity".** In the ultralocal / kinetic
  lattice probes (`docs/CLASSICAL_MAP_INVARIANTS.md` §6), real solutions
  that appear or disappear without a fold are the lattice shadow of
  paths hitting $D_\infty$; the escape-chart makes that literal.

---

## 5. Honest limits — what is NOT claimed

- No statement about continuum $D\ge 1$ QFT, AQFT nets, or constructive
  measures. Compactification here is of the **0D classical field map**.
- $\overline{X}$ is a *partial* compactification (only the $x$-escape
  divisor). It is not a complete toric variety, not $\mathbb{P}(1,1,2)$,
  and not claimed to be the unique maximal compactification.
- $\overline{F}\colon\overline{X}\to\mathbb{A}^3$ is not claimed proper or
  finite as a global map of varieties — only that it is a morphism on the
  escape chart and that $\overline{F}(D_\infty)=S_F$.
- Geometric monodromy of the compactified cover is not re-derived here
  (affine $S_3$ monodromy is Exact: `docs/MONODROMY.md`,
  `scripts/certified_monodromy.py`).
- No claim that every Keller counterexample admits an equally simple
  two-chart form — only that AM's C\* weights make this one canonical.

---

## 6. Reproduce

```bash
.venv/bin/python scripts/weighted_compactification.py   # ~few s, all assertions
```

Related: `scripts/missing_observables.py` (escape certificate),
`scripts/branch_locus.py` (Jelonek $=\{p=0\}$),
`scripts/positive_geometry.py` (reduced wall in $\mathbb{P}^2$),
`jcqft/prefilter.py` (weighted vs unweighted infinity tests).
