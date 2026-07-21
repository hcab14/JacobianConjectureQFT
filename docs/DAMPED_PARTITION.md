# The Damped Partition Function: Finiteness, the Chamber Prefactor, and the Uniformity Exponents of the Semiclassical Limit

*(2026-07-21. Resolves item B5 of `docs/OPEN_QUESTIONS.md`, posed in
`docs/QFT_IMPLICATIONS.md` §5.3 first bullet. All exact claims are
assert-verified and all numerics reproduced by
`scripts/damped_partition.py` (~30 s); figure `outputs/damped_partition.png`.
Background: chamber rule and pushforward `scripts/measure_anomaly.py`,
`docs/POSITIVE_GEOMETRY.md`; escape strata `docs/POSITIVE_GEOMETRY.md` §6.)*

**Summary.** The naive partition function of the counterexample is doubly
ill-posed ($DF$ is not symmetric, so no action exists; the contour/delta
versions are distributional). The *damped* partition function

$$
Z_\hbar(J) \;=\; \int_{\mathbb{R}^3} e^{-|F(\phi)-J|^2/2\hbar}\, d^3\phi
$$

is instead a perfectly behaved constructive object, and we determine it
essentially completely:

1. **Finiteness is unconditional — the Jelonek set is *not* detected by
   divergence.** For every $J \in \mathbb{R}^3$ and every $\hbar > 0$,
   exactly:

$$
Z_\hbar(J) = (2\pi\hbar)^{3/2}\, h(J,\hbar),
\qquad
h(J,\hbar) = \tfrac{1}{2}\,\mathbb{E}\bigl[N(J+\sqrt{\hbar}\,\xi)\bigr]
           = \tfrac12 + \mathbb{P}\bigl[p(J+\sqrt{\hbar}\,\xi)<0\bigr],
$$

   with $\xi \sim \mathcal N(0,\mathbf 1_3)$, hence the uniform two-sided
   bound $\tfrac12 (2\pi\hbar)^{3/2} < Z_\hbar(J) < \tfrac32 (2\pi\hbar)^{3/2}$
   — on the wall and on the empty-fiber cusp orbit included. Tube volumes
   grow like $t^3$ uniformly in $J$.
2. **The leading semiclassical coefficient is the chamber function.**
   Off the wall, $h(J,\hbar) \to N(J)/2 \in \{\tfrac12, \tfrac32\}$ at the
   Gaussian rate $e^{-\mathrm{dist}(J,\{p=0\})^2/2\hbar}$; the jump of the
   prefactor across the Jelonek set is the wall's signature, invisible to
   the $\hbar$-expansion around any single vacuum.
3. **On the wall the limit is the two-sided mean, not the fiber count**:
   $h \to 1$ (with an exact mean-curvature $\sqrt{\hbar}$ correction), and
   at the cusp orbit $h \to \tfrac12$ with an **anomalous $\hbar^{1/4}$**
   correction whose amplitude we compute in closed form. In particular
   $Z_\hbar(0)/(2\pi\hbar)^{3/2} \to 1$: at the perturbative vacuum the
   damped partition function counts *twice* the perturbative saddle — it
   sees the two sheets at infinity.
4. **Uniformity boundary of the semiclassical expansion**:
   $\hbar^*(\varepsilon) \sim \varepsilon^{\gamma}$ with measured
   $\gamma_{\text{wall}} = 2.001$ (exact constant
   $1/\Phi^{-1}(3/4)^2$ matched to $<0.5\%$) and
   $\gamma_{\text{cusp}} = 3.000$ (exact constant
   $27/(826\,\Phi^{-1}(7/8)^2)$ matched to $<1\%$); the unifying rule is
   $\hbar^* \sim \mathrm{dist}(J,\text{wall})^2$, and $\gamma_{\text{cusp}}
   = 2\cdot\tfrac32$ is the $A_2$ horn-width exponent.

Throughout, $N(J)$ is the number of real preimages
($N = 3$ iff $p < 0$, $N = 1$ iff $p > 0$; `jcqft.core.n_real_preimages`),
$p = 27a^2c^2 - 18abc + b^3c - b^2 + 16a$ is the Jelonek/wall polynomial,
and $\Phi$ is the standard normal CDF.

---

## 1. Exact part

### 1.1 Well-definedness: a closed form, not just a bound

The question posed in B5 — does $Z_\hbar$ diverge for $J$ on the
non-properness set, where preimages escape to infinity while $F$ stays
bounded? — has a clean negative answer, and the mechanism is worth stating
precisely because *both* halves of the intuition are correct and they
cancel: the escaping tube $\{|F - \phi| < t\}$ near a wall target is indeed
unbounded in $\phi$-space, but its cross-section shrinks fast enough that
its volume is finite, and in fact exactly accounted for by the target-space
picture.

**Proposition (exact).** $|\det DF| \equiv 2$, so $F:\mathbb{R}^3 \to
\mathbb{R}^3$ is a local diffeomorphism and the co-area/change-of-variables
formula for non-injective local diffeomorphisms gives, for any $g \ge 0$
measurable,

$$
\int_{\mathbb{R}^3} g(F(\phi))\, d^3\phi
= \int_{\mathbb{R}^3} g(J')\, \frac{N(J')}{2}\, d^3J' .
$$

(This is the real pushforward $F_*(d^3\phi) = \tfrac{N}{2}\, d^3J$ of
`scripts/measure_anomaly.py` / `docs/POSITIVE_GEOMETRY.md` §4.) Applying it
to $g = e^{-|\cdot - J|^2/2\hbar}$:

$$
\boxed{\;
Z_\hbar(J) = (2\pi\hbar)^{3/2}
\Bigl(\tfrac12 + \mathbb{P}\bigl[p(J+\sqrt{\hbar}\,\xi) < 0\bigr]\Bigr),
\qquad \xi \sim \mathcal N(0, \mathbf 1_3).
\;}
$$

Since $1 \le N \le 3$ a.e., $Z_\hbar(J)$ is finite for **all** $J$ and all
$\hbar > 0$, bounded between $\tfrac12$ and $\tfrac32$ times the free
Gaussian volume, and (being a Gaussian convolution of a bounded function)
real-analytic in $(J, \hbar)$. Likewise the fiber-tube volumes satisfy

$$
\tfrac{2\pi}{3}\, t^3 \;\le\; V_J(t) = \mathrm{vol}\{\phi : |F(\phi)-J| \le t\}
= \int_{B(J,t)} \tfrac{N}{2}\, d^3J' \;\le\; 2\pi\, t^3
\qquad \text{uniformly in } J .
$$

So: **the damped partition function does not detect the Jelonek set by
divergence** — there is no divergence rate, no cutoff sensitivity, and no
$\hbar$-dependence in the finiteness statement. Its entire sensitivity to
the wall sits in the bounded, analytic prefactor $h(J,\hbar) =
Z_\hbar/(2\pi\hbar)^{3/2}$, whose $\hbar \to 0$ limits are computed below.
(This corrects the natural guess in the task statement; the outcome is
arguably sharper — a *finite* constructive object whose semiclassical data
is discontinuous.)

Two structural remarks, both verified symbolically in the script:

- **Why the escaping tube has finite volume.** $F$ is linear in $z$:
  $F = A(x,y)\,z + B(x,y)$ with
  $A = ((1+xy)^3,\, 3x(1+xy)^2,\, -x^3)$ and
  $|A|^2 = (1+xy)^6 + 9x^2(1+xy)^4 + x^6 > 0$ everywhere. The
  $z$-integral is therefore an exact 1D Gaussian for every $(x,y)$, giving
  the fibration
  $Z = \sqrt{2\pi\hbar}\int |A|^{-1/2} e^{-m(x,y;J)/2\hbar}\, dx\, dy$
  with $m = |A \times (B - J)|^2/|A|^2$ (Lagrange identity). Along the
  escape valley $|A|^{-1/2} \sim |x|^{-3}$: the tube's cross-section decays
  cubically, which is integrability with room to spare — the geometric
  face of the pushforward bound.
- **A Fubini caveat for signed integrands.** The proposition uses $g \ge 0$.
  The same pushforward applied to oscillatory integrands (the contour
  version of the partition function) is exactly what fails to converge;
  the damping is what buys the theorem.

### 1.2 Laplace asymptotics per chamber

For $J$ off the wall the standard Laplace analysis applies to each of the
$N(J)$ global minima of $|F(\phi)-J|^2$ (all nondegenerate: Hessian
$2\,DF^{\mathsf T}DF \succ 0$ with $\det(DF^{\mathsf T}DF) = 4$), giving
per minimum a contribution $(2\pi\hbar)^{3/2}/|\det DF| =
(2\pi\hbar)^{3/2}/2$, hence $Z_\hbar \sim (2\pi\hbar)^{3/2}\, N(J)/2$. The
closed form makes this exact and quantifies the error: since
$h(J,\hbar) - N(J)/2$ is the Gaussian measure of the far chamber,

$$
\Bigl| h(J,\hbar) - \tfrac{N(J)}{2} \Bigr|
\;\le\; e^{-\mathrm{dist}(J,\{p=0\})^2 / 2\hbar} \cdot \text{(polynomial)},
$$

i.e. the semiclassical prefactor is $N(J)/2$ **with no power corrections
at all** — the entire perturbative $\hbar$-series of $\log Z$ beyond the
leading term is chamber-independent (in fact trivial), exactly the
"constant Jacobian $\Rightarrow$ trivial semiclassics" lore — and the only
$J$-dependence of the asymptotics is the exponentially small wall tail
plus the *discontinuous jump of the leading constant*.

### 1.3 On the wall and at the cusp: exact limiting values and corrections

At the wall the fiber count drops ($N = 1$ at generic wall points by the
escape of the $x \sim \pm\sqrt{-q/p}$ pair, $N = 0$ on the cusp orbit where
all three sheets escape like $\varepsilon^{-1}$ — the strata of
`docs/POSITIVE_GEOMETRY.md` §6, re-verified in the script including the
exact medial-line factorization
$pX^3+qX+r = -(\varepsilon X-1)^2(\varepsilon X+2)$ at
$p=-\varepsilon^3,\ q=3\varepsilon,\ r=-2$). But the Gaussian smoothing
averages the *two-sided* chamber structure:

- **Generic wall point** ($p = 0$, $\nabla p \ne 0$): the wall is locally a
  smooth hypersurface, so
  $\mathbb P[p(J + \sqrt\hbar\,\xi) < 0] \to \tfrac12$ and

$$
h(J_{\mathrm{wall}}, \hbar) = 1 + c_1 \sqrt{\hbar} + O(\hbar),
\qquad
c_1 = -\frac{\operatorname{div}\bigl(\nabla p / |\nabla p|\bigr)}
            {2\sqrt{2\pi}}
$$

  (mean curvature of the wall; the sign convention makes $c_1 > 0$ when the
  wall bulges away from the $N=3$ side). At $J_w = (0,1,1)$:
  $c_1 = +0.298590$ analytically; measured $(h-1)/\sqrt\hbar = 0.298590$
  at $\hbar = 10^{-6}$. Note $h \to 1 \neq N(J_w)/2 = \tfrac12$: **the
  escaped pair still carries Gaussian mass $\tfrac12$.** In particular at
  the perturbative vacuum $J = 0$ (which lies on the wall):
  $h(0,\hbar) \to 1$ — the damped partition function assigns the vacuum
  *twice* the perturbative-saddle value $\tfrac12$; the excess is the mass
  of the two sheets at infinity, invisible to every order of the tree
  expansion.
- **Cusp orbit** $J_c$ (empty fiber, $N = 0$): in the C\*-invariants the
  wall has the $A_2$ normal form $P_2 = 27v^2 + 27\,\delta u^3 + \dots$
  ($v = \delta u - \delta w/3$), so the $N=3$ region is a horn of
  half-width $|v| < (s/3)^{3/2}$ at depth $s$. A Gaussian centered at the
  cusp tip gives the **anomalous quarter power**

$$
h(J_c, \hbar) = \tfrac12 + \kappa\, \hbar^{1/4} + o(\hbar^{1/4}),
\qquad
\kappa = \frac{2}{3\sqrt3}\int_0^\infty t^{3/2}\,
\varphi_{25/9}(t)\, \varphi_{18801/18225}\!\bigl(\tfrac{43t}{225}\bigr)\, dt
= 0.124405\,,
$$

  where $\varphi_{V}$ is the centered normal density of variance $V$, and
  the constants are the exact covariances of the source fluctuation in cusp
  coordinates at $J_c = (4/27, 4/3, 1)$: $\operatorname{Var}(\delta w) =
  \tfrac{25}{9}$, $\operatorname{Var}(v) = \tfrac{826}{729}$,
  $\operatorname{Cov}(v, \delta w) = -\tfrac{43}{81}$ (all per unit
  $\hbar$). Measured: $(h - \tfrac12)/\hbar^{1/4} = 0.12444$ at
  $\hbar = 10^{-8}$, matching to $< 1\%$. So the limit value $\tfrac12$
  equals the one-finite-sheet value even though $N(J_c) = 0$ — **the
  entire semiclassical mass over an empty fiber comes from the tube at
  infinity** — and the *rate* $\hbar^{1/4}$ (vs. $\hbar^{1/2}$ generically)
  is the wall-stratum fingerprint, in the same spirit as the trace
  divergence rates of `docs/POSITIVE_GEOMETRY.md` §6.

The stratified dictionary for the prefactor limit:

| stratum | $N(J)$ | $\lim_{\hbar\to0} h$ | correction |
|---|---|---|---|
| chamber $p>0$ | 1 | $1/2$ | $e^{-d^2/2\hbar}$ |
| chamber $p<0$ | 3 | $3/2$ | $e^{-d^2/2\hbar}$ |
| generic wall | 1 | $1$ | $c_1\sqrt{\hbar}$, $c_1$ = mean curvature |
| cusp orbit | 0 | $1/2$ | $\kappa\,\hbar^{1/4}$, $\kappa$ exact above |

$h$ is everywhere the *Gaussian-average* of $N/2$, so its limit is the
Lebesgue-density of $N/2$ at $J$ — chamber value in the interior, two-sided
mean on the wall, horn-weighted mean at the cusp. The fiber count itself is
recovered only as the *non-tangential interior limit*.

### 1.4 How Laplace fails at the wall (escape, not degeneration)

$\det DF = -2$ forbids finite critical-point coalescence: the Hessian of
$|F - J|^2$ at every minimum is uniformly nondegenerate for every $J$. What
breaks as $J \to$ wall from the $N=3$ side is the *domain* of localization:
two minima recede as $x = \pm\sqrt{-q/p} + r/2q + O(p^{1/2})$, i.e.
$|\phi| \sim (q/|\nabla p|)^{1/2}\, \varepsilon^{-1/2}$ at wall distance
$\varepsilon$ (verified numerically to the stated order); near the cusp all
three minima recede like $\varepsilon^{-1}$ along the medial line. Each
receding minimum keeps its full Laplace weight $\tfrac12 (2\pi\hbar)^{3/2}$
until its Gaussian neighborhood collides with the wall-image constraint —
whence a *uniformity* breakdown rather than a coefficient change.

## 2. Numerical part

### 2.1 Method and cross-validation

$h(J,\hbar)$ is computed semi-analytically: $p$ is quadratic in $a$, so the
$a$-marginal probability is exact (normal CDF between the two real roots),
and the $(b,c)$-integral is panel Gauss–Legendre ($32\times 8$ panels/nodes
over $\pm 8.5\sigma$; doubling the resolution changes results at the
$10^{-8}$ level). This was validated against direct $\phi$-space quadrature
of $Z_\hbar$ (exact Gaussian $z$-integral, sinh-mapped grid in $x$ to
capture the escape tails, no change of variables): agreement to
$\lesssim 5\times 10^{-5}$ relative (grid resolution) at chamber, wall, and
cusp points, and stability of the wall-point integral under cutoff growth
$|x| \le 27 \to 548$ at the $10^{-5}$ level — an independent numerical
confirmation of finiteness that does not assume the pushforward identity.

### 2.2 Convergence table: prefactor $\to N(J)/2$

$h(J,\hbar) = Z_\hbar/(2\pi\hbar)^{3/2}$ at two points per chamber
(chamber membership re-verified by direct root count):

| $J$ | $N$ | $\hbar=0.3$ | $0.1$ | $0.03$ | $0.01$ | $0.003$ | $0.001$ | limit |
|---|---|---|---|---|---|---|---|---|
| $(-\tfrac14,0,0)$ | 3 | 1.1308 | 1.2892 | 1.4270 | 1.4939 | 1.50000 | 1.50000 | $3/2$ |
| $(0,2,0)$ | 3 | 1.0411 | 1.2314 | 1.4188 | 1.4925 | 1.49999 | 1.50000 | $3/2$ |
| $(1,0,0)$ | 1 | 0.5357 | 0.5008 | 0.50000 | 0.50000 | 0.50000 | 0.50000 | $1/2$ |
| $(2,1,1)$ | 1 | 0.5001 | 0.50000 | 0.50000 | 0.50000 | 0.50000 | 0.50000 | $1/2$ |

(asserted to $<10^{-4}$ at $\hbar = 10^{-3}$; the approach rate matches the
Gaussian wall-distance tail of §1.2.)

### 2.3 Uniformity exponents $\gamma$

Definition: $\hbar^*(\varepsilon)$ solves $h(J(\varepsilon), \hbar^*) =
1.25$ — the semiclassical value half-lost, $J(\varepsilon)$ in the $N=3$
chamber at distance parameter $\varepsilon$ from the wall/cusp. Bisection
in $\log\hbar$; $h$ is monotone through the crossing.

- **Wall-normal approach** to $J_w = (0,1,1)$
  ($J = J_w - \varepsilon\, \nabla p/|\nabla p|$), $\varepsilon$ from
  $10^{-1.5}$ down to $10^{-4}$:
  $\hbar^*/\varepsilon^2 \to 2.1988$, fit exponent
  $\boxed{\gamma_{\text{wall}} = 2.001}$, and the constant matches the
  exact flat-wall prediction $1/\Phi^{-1}(3/4)^2 = 2.19811$ to $< 0.5\%$
  ($h \approx \tfrac12 + \Phi(\varepsilon/\sqrt{\hbar})$, curvature
  subleading). Equivalently: the $\hbar$-expansion at a chamber point is
  trustworthy iff $\hbar \ll \mathrm{dist}(J, \text{wall})^2$ — which is
  also the naive criterion "the escaping minima at
  $|\phi| \sim \varepsilon^{-1/2}$ sit many Gaussian widths apart",
  self-consistently.
- **Medial (cuspidal-tangent) approach to the cusp**
  ($c = 1$, $w = \tfrac43 - \varepsilon$, $u$ on the $D_0$-line, so
  $p = -\varepsilon^3$ exactly), $\varepsilon$ from $10^{-1}$ down to
  $10^{-10/3}$: $\hbar^*/\varepsilon^3 \to 0.024700$, fit exponent
  $\boxed{\gamma_{\text{cusp}} = 3.000}$, with the constant matching the
  exact horn model

$$
\hbar^* = \frac{\varepsilon^3/27}{\operatorname{Var}(v)\,\Phi^{-1}(7/8)^2}
        = \frac{27\,\varepsilon^3}{826\,\Phi^{-1}(7/8)^2}
        = 0.024697\,\varepsilon^3
$$

  to $< 1\%$ (asserted): the transverse source fluctuation
  $v = \delta a - \tfrac{\delta b}{3} - \tfrac{4\,\delta c}{27}$, of
  standard deviation $\sqrt{826\hbar/729}$, must fit inside the horn
  half-width $\varepsilon^{3/2}/\sqrt{27}$; the depth fluctuation
  $\sim \sqrt{\hbar^*} \sim \varepsilon^{3/2} \ll \varepsilon$ is
  subleading, so the 1D crossing rule is asymptotically exact.
- **Wall-tangent approach** to $J_w$ (distance to wall
  $\sim \varepsilon^2$): fit exponent $4.03 \approx 4$ — consistent with
  the unified rule.

**Unified statement:** $\hbar^* \asymp \mathrm{dist}(J, \{p=0\})^2$ in
every regime tested; the cusp exponent $\gamma_{\text{cusp}} = 3 =
2 \times \tfrac32$ is the square of the $A_2$ horn half-width exponent,
i.e. of the *transverse* distance $\varepsilon^{3/2}/\sqrt{27}$ (the
distance to the wall from a medial point is $\Theta(\varepsilon^{3/2})$
near the cusp, not $\Theta(\varepsilon)$). Note the escape rates
($\varepsilon^{-1/2}$ pair at the wall, $\varepsilon^{-1}$ triple at the
cusp) do **not** enter $\gamma$ — the uniformity boundary is a
target-space (chamber-geometry) quantity, not a field-space one; the
task's a-priori guess that the escape rate sets $\gamma$ is refuted by the
exact constants.

## 3. Synthesis for QFT

In this 0D model there exists a *constructively finite* partition function
— finite for **every** source and every $\hbar$, uniformly bounded by free
Gaussians, real-analytic, no cutoffs, no counterterms — whose semiclassical
expansion nevertheless has a **piecewise-constant leading prefactor**
$N(J)/2$ jumping on the Jelonek set of the classical field map, taking the
anomalous boundary values $1$ (wall) and $\tfrac12$ (empty-fiber cusp,
entirely from field configurations near infinity), with corrections
$c_1\sqrt\hbar$ and $\kappa\,\hbar^{1/4}$ that fingerprint the wall
stratum. Because $\det DF$ is constant, *every local semiclassical datum is
chamber-independent* — around any single vacuum the loop expansion is exact
at one loop and sees nothing — so the jump is carried by no term of the
$\hbar$-series: it is a redistribution of $O(1)$ Gaussian mass between
finite vacua and the boundary of field space, with uniformity boundary
$\hbar \ll \mathrm{dist}(J, \text{wall})^2$ ($\gamma_{\text{wall}} = 2$,
$\gamma_{\text{cusp}} = 3$). The standard lore "constant Jacobian
$\Rightarrow$ trivial semiclassics" is thus true pointwise and false
uniformly: the non-properness of the classical map converts trivial local
semiclassics into a discontinuous global semiclassical normalization. In
Glimm–Jaffe-style constructive terms [GJ87], the model shows that a
well-defined, finite functional integral can have semiclassical data that
no perturbative or single-sector analysis reconstructs — here with the
obstruction located exactly (the wall $\{p=0\}$, its cusp, and the
exponents above). Per the claims-ledger discipline of
`docs/QFT_IMPLICATIONS.md` §4.4: this is a statement about one map in
$D = 0$; it implies nothing about UV renormalization, $D \ge 1$ measures,
or genuine functional integrals, where infinite-field configurations are
suppressed dynamically.

## 4. Limitations

- The closed form for $Z_\hbar$ rests on the chamber rule $N = 3$ iff
  $p<0$ (exact, `scripts/measure_anomaly.py`) and the change-of-variables
  formula for locally-finite non-injective local diffeomorphisms; the
  independent $\phi$-space quadrature confirms it numerically but only at
  grid precision ($\sim 5\times10^{-5}$ relative) and at finitely many
  points.
- The $\sqrt\hbar$ wall coefficient and the $\hbar^{1/4}$ cusp amplitude
  are asymptotic statements verified numerically to $<0.5\%$/$<1\%$ at the
  smallest accessible $\hbar$ ($10^{-6}$/$10^{-8}$); the error terms
  ($O(\hbar)$, $o(\hbar^{1/4})$) are observed, not proved.
- $\gamma$ values are measured over 3–7 decades of $\varepsilon$ with fit
  exponents within $10^{-2}$ of the integers 2, 3 and constants matching
  the exact predictions; they are not certified proofs. The tangential
  exponent 4 is tested over one decade only ($\varepsilon \le 0.1$ is
  forced: the tangent line exits the chamber).
- The half-loss threshold $h = 1.25$ in the definition of $\hbar^*$ is a
  convention; any fixed level in $(1, \tfrac32)$ changes the constants
  ($\Phi^{-1}$ values) but not the exponents.
- Nothing here resolves the individual sheets at infinity beyond their
  total mass; a refinement (e.g. clustering the Gaussian mass by
  field-space region as $\hbar \to 0$) would separate the receding pair
  from the finite sheet and is left open.

## 5. References

- [GJ87] J. Glimm, A. Jaffe, *Quantum Physics: A Functional Integral Point
  of View*, 2nd ed., Springer (1987). (Constructive framing: finiteness
  and bounds first, asymptotics second.)
- [Jel93] Z. Jelonek, *The set of points at which a polynomial map is not
  proper*, Ann. Polon. Math. **58** (1993) 259–266. (The wall $\{p=0\}$ as
  the non-properness set.)
- Repository: `scripts/measure_anomaly.py` (chamber rule, pushforward),
  `docs/POSITIVE_GEOMETRY.md` (cusp geometry, escape strata §6),
  `docs/MONODROMY.md` ($S_3$ sheet structure),
  `docs/QFT_IMPLICATIONS.md` §5.3 (the question), `docs/OPEN_QUESTIONS.md`
  B5.
