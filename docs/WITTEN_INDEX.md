# The 0D Witten Index: Signed Count, Mathai–Quillen Localization, and Index Jumping at the Non-Properness Wall

*(2026-07-25. Verified by `scripts/witten_index.py` (~10 s; every displayed
identity is asserted there, numerics labelled). This makes exact the
"Witten index" reading of the counterexample: the SUSY-localized partition
function of the 0D Mathai–Quillen model computes the signed solution count
(Brouwer degree), which jumps across the non-properness wall — index
non-invariance **is** non-properness. Background: chamber rule
`scripts/measure_anomaly.py`; damped partition function
`docs/DAMPED_PARTITION.md`; wall = Jelonek set `scripts/branch_locus.py`;
no action exists `docs/QFT_IMPLICATIONS.md` §5.3. Not to be confused with
the properness **prefilter** informally called "Witten-index prefilter" in
`docs/RIGIDITY_AND_PREFILTER.md` §1 — this note is about the index itself.)*

**Summary.** For the Alpöge–Mathew map $F:\mathbb{R}^3\to\mathbb{R}^3$
(det $DF \equiv -2$), with $N(J)$ the real preimage count and
$p = 27a^2c^2 - 18abc + b^3c - b^2 + 16a$ the wall polynomial:

1. **Exact.** The signed solution count ("0D Witten index") is
   $$\deg(F, J) \;=\; \sum_{\phi\,\in\,F^{-1}(J)} \operatorname{sign}\det DF(\phi)
   \;=\; -N(J) \;=\; \begin{cases} -1 & p(J) > 0,\\ -3 & p(J) < 0,\end{cases}$$
   asserted via exact fibers at rational points per chamber (every
   preimage verified to map back exactly).
2. **Exact.** The jump $-1 \to -3$ is a *certificate of non-properness*:
   every $J$ is a regular value and every fiber is finite, and for a
   proper map the signed count is a single integer (degree of the
   extension $S^3 \to S^3$). This certificate is independent of the
   escape-curve certificate of `docs/MISSING_OBSERVABLES.md` §2
   (re-asserted in one line). Conversely, the wall $\{p=0\}$ *is* the
   Jelonek non-properness set (`scripts/branch_locus.py`).
3. **Exact.** The bosonic form of the Mathai–Quillen partition function,
   $Z_\sigma(J) = (2\pi\sigma^2)^{-3/2}\int \det DF\, e^{-|F-J|^2/2\sigma^2} d^3\phi$,
   has the closed form $Z_\sigma(J) = -\,\mathbb{E}[N(J+\sigma\xi)]$
   ($\xi \sim \mathcal N(0,\mathbf 1_3)$): it is the **Gaussian
   mollification of the index**, finite for every $J$ and $\sigma > 0$
   with $-3 < Z_\sigma < -1$, and $Z_\sigma(J) \to \deg(F,J)$ off the wall
   at the exact rate $e^{-\mathrm{dist}(J,\{p=0\})^2/2\sigma^2}$. On the
   wall the limit is $-2$; at the empty-fiber cusp orbit it is $-1$ —
   pure boundary-at-infinity contributions.
4. **Exact.** No superpotential exists ($DF \neq DF^{\mathsf T}$, one
   line), so there is no Parisi–Sourlas model; the Mathai–Quillen
   completion (nilpotent BRST charge, $Q$-exact action) exists for this
   arbitrary $F$ and is verified generator-by-generator in an explicit
   Grassmann algebra, including the Berezin integral
   $\int e^{i\bar\chi M\psi} = -i\det M$ and the normalization that
   reduces the MQ integral to $Z_\sigma$.
5. **Numerical evidence.** Direct $\phi$-space quadrature of the MQ
   integral matches the closed form to $\lesssim 5\times10^{-5}$;
   convergence tables $Z_\sigma \to -1 / -3$ per chamber; the measured
   Gaussian-decay distance matches $\mathrm{dist}(J, \text{wall})$ to
   ~2%; the near-wall crossover matches the flat-wall profile
   $-1 - 2\Phi(\varepsilon/\sigma)$.
6. **Exact contrast.** Over $\mathbb{C}$ the fiber count is 3 for every
   $J$ off $\{p=0\}$ — no jump; the wall is a real-locus phenomenon
   (discriminant sign of the eliminant cubic).

---

## 0. Setting

$F$ is the Alpöge–Mathew map (`jcqft/core.py`), $\det DF \equiv -2$, so
$F$ is an orientation-reversing local diffeomorphism of $\mathbb{R}^3$.
Every preimage of $J = (a,b,c)$ has $x$-coordinate a root of the eliminant
cubic $p\,X^3 + q\,X + r$ ($q = 4 - 3bc$, $r = -2c$), with $y, z$ rational
in the root off $\{D_0 = 0\}$; the finite-fiber description is exhaustive
(`jcqft/fibers.py`, `docs/MISSING_OBSERVABLES.md` §4). The **chamber
rule** (exact, `scripts/measure_anomaly.py`): $N(J) = 3$ iff $p(J) < 0$,
$N(J) = 1$ iff $p(J) > 0$, from $4q^3 + 27pr^2 = 4D_0^2$.

## 1. Exact part

### 1.1 The signed count is $-N(J)$: chamber-wise Brouwer degree

Since $\det DF = -2$ *identically*, $\operatorname{sign}\det DF(\phi) = -1$
at every point — in particular at every real solution of $F(\phi) = J$ —
and every $J \in \mathbb{R}^3$ is a regular value. Hence the signed count
is defined for every $J$ and equals $-N(J)$: the index is minus the vacuum
count, with $|\deg| = N$ and the sign recording that $F$ reverses
orientation. Asserted at two rational targets per chamber, all with
$D_0 \neq 0$ and exact map-back of every fiber point (reduction modulo the
minimal polynomial of the eliminant root):

| $J$ | $p(J)$ | complex fiber | real fiber $N$ | $\deg(F,J)$ |
|---|---|---|---|---|
| $(-\tfrac14, 0, 0)$ | $-4$ | 3 | 3 | $-3$ |
| $(0, 2, 0)$ | $-4$ | 3 | 3 | $-3$ |
| $(1, 0, 0)$ | $16$ | 3 | 1 | $-1$ |
| $(2, 1, 1)$ | $104$ | 3 | 1 | $-1$ |

(The identity $\deg = -N$ holds for **all** $J$ off the wall by the
chamber rule, since the sign factor is constant; the table is the assert
anchor, not the proof perimeter.)

### 1.2 The index jump is an exact certificate of non-properness

**Degree theory (standard).** A proper $C^1$ map $f:\mathbb{R}^n \to
\mathbb{R}^n$ extends continuously to the one-point compactifications
$S^n \to S^n$ (properness is exactly what sends $\infty \mapsto \infty$),
and at every regular value with finite fiber the signed count equals the
topological degree of the extension — one integer, independent of the
value [OR09, Mil65]. For non-proper $f$ the signed count is still defined
pointwise (as here) but no invariance is claimed — and none holds.

**Contrapositive certificate.** For the AM map every $J$ is regular and
every fiber is finite, yet the signed count takes both values $-1$ and
$-3$. Therefore $F$ is **not proper** — a certificate independent of the
escape-curve certificate of `docs/MISSING_OBSERVABLES.md` §2, which is
nevertheless re-asserted in the script in one line: along
$\phi(T) = (T,\, y_0,\, (2T - 3T^2y_0 - c_3)/T^3)$ the field escapes while
$F(\phi(T)) \to \bigl(y_0^2(1 - c_3y_0),\, y_0(4 - 3c_3y_0),\, c_3\bigr)$,
a finite point with $p = 0$.

**Converse link.** The wall $\{p = 0\}$ is exactly the Jelonek set: escape
happens only in $x$, and the leading coefficient of the $x$-eliminant is
$p$ (`scripts/branch_locus.py` — cross-referenced, not re-derived). So
"where the index can jump" and "where $F$ fails to be proper" are the
same hypersurface, and both statements are exact.

### 1.3 The Mathai–Quillen partition function: closed form and finiteness

Define (bosonic form; the fermionic origin is §2.2)

$$
Z_\sigma(J) \;=\; (2\pi\sigma^2)^{-3/2} \int_{\mathbb{R}^3}
\det DF(\phi)\; e^{-|F(\phi) - J|^2/2\sigma^2}\, d^3\phi .
$$

Because $\det DF = -2 = -|\det DF|$ comes out of the integral, the exact
pushforward $F_*(d^3\phi) = \tfrac{N}{2}\,d^3J$
(`scripts/measure_anomaly.py`, `docs/DAMPED_PARTITION.md` §1.1) gives the
closed form

$$
\boxed{\;
Z_\sigma(J) \;=\; -\,\mathbb{E}\bigl[N(J + \sigma\xi)\bigr]
\;=\; \int_{\mathbb{R}^3} \deg(F, J')\,\varphi_\sigma(J' - J)\, d^3J' ,
\qquad \xi \sim \mathcal N(0, \mathbf 1_3),
\;}
$$

i.e. **the MQ partition function is the Gaussian mollification of the
index**. Consequences, all exact:

- **Finiteness is unconditional**: $-3 < Z_\sigma(J) < -1$ for every $J$
  (wall and cusp included) and every $\sigma > 0$. The worry that the
  integral diverges along escape directions — $|F(\phi) - J|$ does *not*
  tend to $\infty$ along them — is resolved exactly as for the damped
  partition function: the escaping tube is unbounded but its
  cross-section decays like $|A|^{-1/2} \sim |x|^{-3}$
  (`docs/DAMPED_PARTITION.md` §1.1; the fibration ingredients are
  re-asserted in the script). So the Jelonek set is *not* detected by
  divergence here either.
- **Localization per chamber**: for $J$ off the wall, with
  $d = \mathrm{dist}(J, \{p=0\})$,
  $$\bigl|Z_\sigma(J) - \deg(F, J)\bigr| \;=\; 2\,\mathbb{P}[p \text{ flips sign}]
  \;\le\; 2\Bigl(1 + \sqrt{\tfrac{2}{\pi}}\,\tfrac{d}{\sigma}\Bigr)
  e^{-d^2/2\sigma^2} \;\longrightarrow\; 0 ,$$
  so $Z_\sigma(J) \to \deg(F,J)$ as $\sigma \to 0$: the SUSY-localized
  partition function computes the Brouwer degree, chamber by chamber.
- **Boundary contributions exactly on the wall**: at a generic wall point
  (and at the vacuum $J = 0$, which lies on the wall)
  $Z_\sigma \to -2$, the two-sided mean, *not* $-N(J_{\rm wall}) = -1$:
  the escaped pair still carries index $-1$ "from infinity". At the cusp
  orbit ($N = 0$, empty fiber) $Z_\sigma \to -1$: the *entire* limit is a
  boundary term. (These are the $h \to 1$, $h \to \tfrac12$ limits of
  `docs/DAMPED_PARTITION.md` §1.3 times $-2$; corrections
  $c_1\sqrt{\hbar}$, $\kappa\,\hbar^{1/4}$ with $\hbar = \sigma^2$.)
- **Properness = $\sigma$-independence.** If $F$ were proper, the signed
  count would be a constant $\deg$ and the mollification formula would
  give $Z_\sigma(J) = \deg$ *exactly, for every $\sigma$ and $J$* (as it
  does for linear maps). The $\sigma$-dependence of $Z_\sigma$ — the
  crossover profile of §3.4 — is therefore a direct measurement of the
  wall. In BRST language: $\partial_{\sigma^2} Z$ is the integral of a
  $\delta$-exact form, which vanishes only if integration by parts in
  field space produces no boundary terms — precisely properness. (This
  sentence is the formal/interpretive face; the closed form above is the
  exact statement.)

## 2. The SUSY structure that exists — and the one that does not

### 2.1 No superpotential: $DF \neq DF^{\mathsf T}$ (exact)

$$(DF)_{13} - (DF)_{31} = x^3y^3 + 3x^2y^2 + 3x^2z + 9xy - 1 \;\not\equiv\; 0$$

(equal to $-1$ at $\phi = 0$: not symmetric even at the vacuum). So there
is no $W$ with $F = \nabla W$, even locally, and the standard 0D
Parisi–Sourlas completion — action $\tfrac12|\nabla W - J|^2$ with fermion
bilinear $\bar\psi\,(\mathrm{Hess}\,W)\,\psi$ [PS79, PS82] — **does not
exist for this map** (cf. `docs/QFT_IMPLICATIONS.md` §5.3: no action, the
dynamics is not variational).

### 2.2 The Mathai–Quillen completion (exists for arbitrary $F$)

Fields: $\phi \in \mathbb{R}^3$, auxiliary $B \in \mathbb{R}^3$, fermions
$\psi_i, \bar\chi_i$ ($i = 1,2,3$). BRST-type charge $\delta$:

$$
\delta\phi_i = \psi_i, \quad \delta\psi_i = 0, \quad
\delta\bar\chi_i = B_i, \quad \delta B_i = 0
\qquad (\delta^2 = 0 \text{ on everything}).
$$

Gauge fermion and ($Q$-exact) action:

$$
S \;=\; \delta\Bigl[\,\bar\chi\cdot\Bigl(i\bigl(F(\phi) - J\bigr)
+ \tfrac{\sigma^2}{2}B\Bigr)\Bigr]
\;=\; i\,B\cdot(F(\phi) - J) \;+\; \tfrac{\sigma^2}{2}|B|^2
\;-\; i\,\bar\chi\cdot DF(\phi)\,\psi ,
$$

$$
Z^{\rm MQ}_\sigma(J) \;=\; \frac{i}{(2\pi)^3} \int
d^3\phi\, d^3B\, d^3\psi\, d^3\bar\chi\;\, e^{-S} .
$$

All of this is verified in the script with an explicit 6-generator
Grassmann (Berezin) algebra — no formal manipulation is trusted:

- $\delta^2 = 0$ on the gauge fermion and on a generic test element,
  and $\delta S = 0$ (the action is $\delta$-closed because it is
  $\delta$-exact) — **asserted**;
- $S = \delta\Psi$ reproduces the displayed component form — **asserted**
  (note the would-be Yukawa term $\propto \partial_k(DF)_{ij}\,
  \psi_k\psi_j\bar\chi_i$ cancels by symmetry of second derivatives:
  this is $\delta^2 = 0$ at work);
- Berezin integral $\int d^3\psi\, d^3\bar\chi\; e^{\,i\bar\chi M\psi}
  = -i\det M$ for a *generic* symbolic $3\times3$ matrix $M$
  (convention: coefficient of $\psi_1\psi_2\psi_3\bar\chi_1\bar\chi_2
  \bar\chi_3$) — **asserted**;
- Gaussian $B$-integral
  $\int e^{-iBd - \sigma^2B^2/2}\,dB = \tfrac{\sqrt{2\pi}}{\sigma}
  e^{-d^2/2\sigma^2}$ — **asserted** (sympy);
- the normalization $\tfrac{i}{(2\pi)^3}$ makes the composite exactly
  $Z^{\rm MQ}_\sigma = Z_\sigma$ of §1.3 — **asserted**.

Because $\det DF \equiv -2$ is constant, the ghost bilinear integrates to
a *constant*: the fermionic sector of this model is trivial (it
contributes the factor $-2$ and the overall sign of the index), which is
what makes the model exactly solvable — the entire nontrivial content
sits in the bosonic escape geometry.

### 2.3 The index of the MQ model

The "Witten index" of this 0D model — its SUSY-localized partition
function $\lim_{\sigma\to0} Z^{\rm MQ}_\sigma(J)$ — is $\deg(F, J)$: it
equals $-1$ or $-3$ depending on the chamber of $J$, is *not* a deformation
invariant in $J$, and its jump locus is exactly the Jelonek set. This is
the content of §§1.1–1.3 assembled; the MQ formalism [MQ86, Bla93] is what
licenses calling the signed count an "index" of a supersymmetric model in
the first place.

## 3. Numerical part

All items in this section are **numerical evidence** (quadrature at
finitely many points), reproducing `scripts/witten_index.py` §5.

### 3.1 Direct $\phi$-space check of the MQ integral

Direct quadrature of $Z_\sigma$ (exact Gaussian $z$-integral, sinh-mapped
$x$-grid for the escape tails, *no* change of variables — the method
cross-validated in `docs/DAMPED_PARTITION.md` §2.1) against the closed
form, at $\sigma^2 = 0.05$: relative error $-4.2\times10^{-6}$ at
$(-\tfrac14,0,0)$, $-1.1\times10^{-9}$ at $(1,0,0)$, $-5.4\times10^{-5}$
at the wall point $(0,1,1)$. The MQ integral itself converges and matches,
independently of the pushforward identity.

### 3.2 Convergence to the index

| $J$ | $\sigma=0.5$ | $0.3$ | $0.2$ | $0.1$ | $0.05$ | $0.03$ | $\deg$ |
|---|---|---|---|---|---|---|---|
| $(-\tfrac14,0,0)$ | $-2.3249$ | $-2.6036$ | $-2.7932$ | $-2.9878$ | $-3.0000$ | $-3.0000$ | $-3$ |
| $(0,2,0)$ | $-2.1385$ | $-2.5016$ | $-2.7661$ | $-2.9850$ | $-3.0000$ | $-3.0000$ | $-3$ |
| $(1,0,0)$ | $-1.0484$ | $-1.0009$ | $-1.0000$ | $-1.0000$ | $-1.0000$ | $-1.0000$ | $-1$ |
| $(2,1,1)$ | $-1.0001$ | $-1.0000$ | $-1.0000$ | $-1.0000$ | $-1.0000$ | $-1.0000$ | $-1$ |

(asserted to $<10^{-3}$ at $\sigma = 0.03$).

### 3.3 Rate

Fitting $|Z_\sigma - \deg| \sim e^{-d^2/2\sigma^2}$: $d_{\rm fit} = 1.022$
vs. $\mathrm{dist}(J,\text{wall}) = 1.000$ at $(1,0,0)$, and $0.2555$ vs.
$0.2500$ at $(-\tfrac14,0,0)$ (grid-refined distance; both ratios $1.022$,
the drift being the polynomial prefactor). The exact two-sided bound is
§1.3.

### 3.4 Wall, vacuum, cusp, and the crossover

$Z_\sigma$ at $\sigma = 0.1,\ 10^{-2},\ 10^{-3}$:

| $J$ | $0.1$ | $10^{-2}$ | $10^{-3}$ | limit |
|---|---|---|---|---|
| generic wall $(0,1,1)$ | $-1.60285$ | $-2.00597$ | $-2.00060$ | $-2$ |
| vacuum $J=0$ (on wall) | $-2.00499$ | $-2.00050$ | $-2.00005$ | $-2$ |
| cusp $(\tfrac{4}{27},\tfrac43,1)$, $N=0$ | $-1.09135$ | $-1.02526$ | $-1.00788$ | $-1$ |

Approach at $J = J_{\rm wall} - \varepsilon\,n$ ($N=3$ side,
$\varepsilon = 10^{-3}$): $Z = -2.1614,\ -2.6831,\ -3.0000$ at
$\sigma = 5\varepsilon,\ \varepsilon,\ 0.2\varepsilon$, matching the
flat-wall profile $-1 - 2\Phi(\varepsilon/\sigma)$ to $< 0.01$. The index
is resolved only once $\sigma \ll \mathrm{dist}(J, \text{wall})$ — the
same uniformity boundary $\hbar^* \sim \mathrm{dist}^2$ measured in
`docs/DAMPED_PARTITION.md` §2.3.

## 4. Interpretation (flagged as interpretation)

The index jumps because **two classical vacua escape through infinity**
as $J$ crosses the wall: on the $N=3$ side the extra pair sits at
$x \approx \pm\sqrt{-q/p}$, recedes like $\mathrm{dist}^{-1/2}$, and its
index contribution $-2$ is handed from the finite-vacuum sum to the
boundary of field space (wall value $-2$ = mean; cusp value $-1$ entirely
from infinity). This is the 0D caricature of the standard mechanisms:

- **Witten index jumping in SUSY QM with non-compact targets** — the
  index is invariant under deformations that keep the asymptotics tame,
  and can jump when vacua run off to infinity [Wit82]; here "deforming
  the theory" is moving the source $J$, and the tame-asymptotics
  hypothesis fails exactly on the Jelonek set.
- **Wall-crossing** — a piecewise-constant integer invariant of a family,
  jumping on a real-codimension-one wall in parameter space
  [KS08, CV93]; here the invariant is $\deg(F,J) = -N(J)$, the wall is
  $\{p = 0\}$, and the jump $\Delta = 2$ is the escaping pair. No claim
  of a KS-type wall-crossing *formula* is made — in this 0D model the
  jump is simply the count of escaping solutions.

These parallels are structural analogies, not derivations; everything
quantitative in this note is §§1–3.

## 5. Relation to wall-crossing and SUSY QM: references

- Witten index and its non-invariance under changes of asymptotic data:
  [Wit82] (esp. the discussion of systems with non-compact field space
  and states running to infinity); functional-measure viewpoint [CG82].
- Mathai–Quillen formalism: [MQ86]; QFT translation and the finite-
  dimensional "index = degree" toy model: [Bla93] (the AM map is exactly
  such a toy model, with the novelty that non-properness makes the
  degree value-dependent).
- Parisi–Sourlas / stochastic quantization form of the localized
  integral (needs $F = \nabla W$; unavailable here by §2.1):
  [PS79, PS82].
- Degree theory used in §1.2: [OR09] (degree on open sets and proper
  maps), [Mil65] (compact case).
- Wall-crossing as jumping of integer invariants on real walls: [KS08];
  BPS-count jumping in 2d $(2,2)$ theories: [CV93].

## 6. Limitations

- The chamber-wise degree is asserted at two rational points per chamber;
  constancy on each chamber follows from the exact chamber rule
  ($\deg = -N$, `scripts/measure_anomaly.py`), not from sampling. The
  degree-theory input of §1.2 (proper $\Rightarrow$ value-independent
  signed count) is cited, not formalized in sympy.
- The closed form for $Z_\sigma$ rests on the same two exact inputs as
  `docs/DAMPED_PARTITION.md`: the chamber rule and the change-of-variables
  formula for non-injective local diffeomorphisms. The direct quadrature
  check (§3.1) is independent of both but is numerics at grid precision,
  at finitely many $(J, \sigma)$.
- The Grassmann verification of §2.2 establishes the *algebra* (BRST
  exactness, Berezin integrals, normalization) exactly; the *measure
  theory* of the $B$- and fermion-integrals is finite-dimensional and
  elementary, but the interchange of the $B$-Gaussian with the
  $\phi$-integral is justified by absolute convergence
  ($|e^{-iBd}| = 1$, Fubini on a product of Gaussians), not asserted
  symbolically.
- The wall-distance values in §3.3 come from grid refinement (4 levels);
  they are numerically clean ($1.0000$, $0.2500$) but not certified
  exact, and the rate fit carries the usual prefactor drift (~2%).
- The $\sigma$-independence-breaking statement in §1.3 (boundary terms in
  the $Q$-exactness argument) is stated as the standard formal mechanism;
  the exact content is the closed form, and no independent verification
  of the boundary-term identity is attempted.
- Everything is $D = 0$ and about this one map; no claim about $D \ge 1$,
  genuine SUSY QM, or wall-crossing formulas.

## 7. Reproduce

```bash
.venv/bin/python scripts/witten_index.py    # ~10 s, all assertions + tables
```

## References

- [Wit82] E. Witten, *Constraints on supersymmetry breaking*, Nucl. Phys.
  B **202** (1982) 253–316.
- [CG82] S. Cecotti, L. Girardello, *Functional measure, topology and
  dynamical supersymmetry breaking*, Phys. Lett. B **110** (1982) 39–43.
- [MQ86] V. Mathai, D. Quillen, *Superconnections, Thom classes, and
  equivariant differential forms*, Topology **25** (1986) 85–110.
- [Bla93] M. Blau, *The Mathai–Quillen formalism and topological field
  theory*, J. Geom. Phys. **11** (1993) 95–127 (hep-th/9203026).
- [PS79] G. Parisi, N. Sourlas, *Random magnetic fields, supersymmetry,
  and negative dimensions*, Phys. Rev. Lett. **43** (1979) 744–745.
- [PS82] G. Parisi, N. Sourlas, *Supersymmetric field theories and
  stochastic differential equations*, Nucl. Phys. B **206** (1982)
  321–332.
- [KS08] M. Kontsevich, Y. Soibelman, *Stability structures, motivic
  Donaldson–Thomas invariants and cluster transformations*,
  arXiv:0811.2435 (2008).
- [CV93] S. Cecotti, C. Vafa, *On classification of N=2 supersymmetric
  theories*, Comm. Math. Phys. **158** (1993) 569–644.
- [Mil65] J. Milnor, *Topology from the Differentiable Viewpoint*, Univ.
  Press of Virginia (1965).
- [OR09] E. Outerelo, J. M. Ruiz, *Mapping Degree Theory*, Graduate
  Studies in Mathematics **108**, AMS (2009).
- [Jel93] Z. Jelonek, *The set of points at which a polynomial map is not
  proper*, Ann. Polon. Math. **58** (1993) 259–266.
- *(NOTE: bibliographic data best-effort from memory; page ranges for
  [CG82] and [CV93] should be double-checked before external citation.)*
- Repository: `scripts/witten_index.py` (this note),
  `scripts/measure_anomaly.py` (chamber rule),
  `scripts/branch_locus.py` (Jelonek set),
  `scripts/missing_observables.py` (escape certificate),
  `docs/DAMPED_PARTITION.md` (pushforward, finiteness machinery, wall/cusp
  corrections), `docs/QFT_IMPLICATIONS.md` §5.3 (no action exists),
  `scripts/pushforward_forms.py` (all-sheets/complex pushforward
  $F_*(d^3\phi) = -\tfrac32\,d^3J$: the constant, wall-blind complex
  counterpart of the jumping real index).
