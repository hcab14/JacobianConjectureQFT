# The D=1 Mathai–Quillen Index in a Finite-Mode Truncation: the Jump Survives the Path Measure

*(2026-07-28. Answers Q1 of `docs/CLASSICAL_MAP_INVARIANTS.md` §5.1 in its
smallest honest setting — a finite Fourier-mode (Galerkin) truncation of
the $D=1$ SUSY-QM Mathai–Quillen model for the Alpöge–Mathew force map.
All exact claims assert-verified and all numerics reproduced by
`scripts/d1_index_modes.py` (37 checks, ~4.5 min; `--full` ~25 min).
Background: 0D index `docs/WITTEN_INDEX.md` / `scripts/witten_index.py`;
damped measures `docs/DAMPED_PARTITION.md`; the Q1 posing
`docs/CLASSICAL_MAP_INVARIANTS.md` §5.1.)*

**Summary.** For periodic paths $q:[0,\beta]\to\mathbb{R}^3$ truncated to
Fourier modes $\le M$, the Mathai–Quillen partition function of the
first-order flow $\dot q + F(q) - J = 0$ (AM force map $F$,
$\det DF\equiv-2$) is a finite-dimensional integral in
$n = 3(2M{+}1)$ mode coordinates. We determine it:

1. **Exact.** At $M=0$ the model reduces *exactly* to the 0D MQ integral:
   $Z_\sigma(J; 0, \beta) = Z^{0\mathrm{D}}_{\sigma/\sqrt{\beta}}(J)$
   (symbolic assert; the kinetic term contributes nothing on constant
   paths and the orthonormal mode normalization produces precisely the
   effective coupling $\sigma_{\rm eff} = \sigma/\sqrt\beta$).
2. **Exact (saddle factorization — the theorem-shaped core).** At a
   constant-path zero $q^*$ ($F(q^*)=J$) the linearized map is
   block-diagonal over modes with $k$-block
   $\bigl[\begin{smallmatrix}A & \omega_k I\\ -\omega_k I &
   A\end{smallmatrix}\bigr]$, $A = DF(q^*)$,
   $\omega_k = 2\pi k/\beta$, and
   $$\det\begin{pmatrix}A & \omega I\\ -\omega I & A\end{pmatrix}
   = \det(A^2 + \omega^2 I) = \bigl|\det(A + i\omega I)\bigr|^2 \;>\; 0$$
   (generic $3\times3$ $A$, real-$A$ positivity, both symbolic). Hence
   the mode fluctuation determinant is a **$J$-independent positive
   factor**: $\operatorname{sign}\det DG(u^*) = \operatorname{sign}\det
   DF(q^*)$, and the $\sigma\to0$ constant-saddle contribution is
   $\sum_{q^*}\operatorname{sign}\det DF(q^*) = \deg(F,J) = -N(J)$,
   **independent of $M$ and $\beta$** — provided
   $\det(A^2+\omega_k^2 I)\neq 0$ (proviso P3, checked) and no
   nonconstant zeros or boundary mass intervene (provisos P1–P2, probed).
3. **Exact (gradient rigidity) + the AM loophole.** At *any* zero of the
   truncated flow map, $\lVert\dot q\rVert^2 = -\langle\dot q,
   F(q)\rangle$. For a gradient force $F=\nabla W$ the right side is a
   loop integral of $dW$ and vanishes (asserted symbolically for the tame
   shear at $M=1$, generic coefficients): all truncated zeros are constant
   equilibria, every $M,\beta$. For AM ($DF\neq DF^{\mathsf T}$) the cross
   term is nonzero on an explicit witness loop ($75\pi/16$, exact) — the
   rigidity proof does **not** apply, so nonconstant zeros are excluded
   only numerically (multi-start Newton, none found).
4. **Numerical (MC with error bars).** On the grid $\sigma\in\{0.5, 0.25,
   0.1, 0.05\}$, $\beta\in\{0.5,1,2\}$, $M\in\{1,2\}$, at the four
   rational chamber points of `scripts/witten_index.py`:
   $Z_\sigma(J; M,\beta)\to -N(J)$ per chamber and the jump ratio
   $Z(J_-)/Z(J_+)\to 3$ (at $\sigma=0.05$: $2.994$–$3.006\pm0.015$ for
   every probed $(M,\beta)$). The values track the 0D closed form at
   $\sigma_{\rm eff}$ to a few percent — the higher modes are
   near-spectators. Proper-map controls (linear $L$, tame shear) give
   $Z\equiv-1$, fixing the normalization.
5. **Numerical (wall approach).** Along $J=(a,0,0)$, $a\to0^-$, the MC
   values match the 0D crossover $Z^{0\mathrm{D}}_{\sigma_{\rm eff}}$
   within error bars while the escaping pair recedes
   ($x^*=\pm(−4a)^{-1/2}\dots$ up to $8$); the $|w|$-mass beyond the
   equilibria in the escape coordinate stays $\lesssim 3\%$ and the mass
   in genuinely nonconstant directions stays at MC-noise level: no new
   $D=1$ escape channel opens in the truncation.
6. **Verdict for Q1 (truncated model): the jump SURVIVES.** The
   finite-mode path measure does not kill the vacua at infinity. No
   continuum ($M\to\infty$) claim is made.

---

## 0. The model (all conventions explicit)

Periodic paths truncated to $M$ Fourier harmonics,
$$
q(\tau) = q_0 + \sum_{k=1}^{M}\Bigl[a_k\cos\tfrac{2\pi k\tau}{\beta}
 + b_k\sin\tfrac{2\pi k\tau}{\beta}\Bigr],\qquad
q_0, a_k, b_k \in\mathbb{R}^3,
$$
written in the $L^2([0,\beta])$-**orthonormal** scalar basis
$e_0=\beta^{-1/2}$, $e_k^c = (2/\beta)^{1/2}\cos\omega_k\tau$,
$e_k^s=(2/\beta)^{1/2}\sin\omega_k\tau$, $\omega_k = 2\pi k/\beta$, with
coordinates $u\in\mathbb{R}^n$, $n = 3(2M{+}1)$.

**Truncated Nicolai / flow map** (the $D=1$ MQ localization datum for
$\delta q = \dot q + F(q) - J$; sign convention $\dot q + F - J$, *not*
$J - F$):
$$
G_{J,\beta,M}:\mathbb{R}^n\to\mathbb{R}^n,\qquad
G(u)_{i\alpha} = \bigl\langle e_\alpha,\ (\dot q + F(q) -
J)_i\bigr\rangle_{L^2},
$$
i.e. the orthonormal coefficients of $P_M[\dot q + F(q) - J]$, $P_M$ the
orthogonal projection onto modes $\le M$ (needed to make the map square).

**Partition function** (finite-dimensional MQ integral; $\det$, not
$|\det|$):
$$
\boxed{\;Z_\sigma(J; M, \beta) = (2\pi\sigma^2)^{-n/2}
\int_{\mathbb{R}^n} \det DG(u)\,
e^{-|G(u)|^2/2\sigma^2}\, d^n u.\;}
$$
The constant $c = (2\pi\sigma^2)^{-n/2}$ is *fixed by the proper-map
controls*: for the linear map $L\phi=(z,y,2x)$, $G$ is affine with
constant $\det DG = \det L\prod_k\det(L^2+\omega_k^2 I) < 0$, and for the
tame shear $\nabla W = (x, z+y^2, y)$ (a polynomial automorphism), the MC
value is $-1$ within error bars at every probed $(M,\beta,\sigma)$
(§2.1) — the truncated MQ index of a proper map is the constant $-1$.

Mode projections are **exact** in the numerics: all integrands are trig
polynomials of degree $\le 8M$ ($\deg F = 7$), and the rectangle rule on
$N_\tau = 64$ points integrates trig polynomials of degree $<64$ exactly
(asserted against $N_\tau=512$ to $3\times10^{-14}$).

**Structural frame (exact, cited co-area input).** Pushing $u\mapsto G(u)$
forward, $Z_\sigma(J;M,\beta) = \mathbb{E}_\xi\bigl[\deg_G(\sigma\xi)\bigr]$
with $\xi\sim\mathcal N(0,\mathbf 1_n)$ and $\deg_G(v) =
\sum_{u\in G^{-1}(v)}\operatorname{sign}\det DG(u)$ — the Gaussian
mollification of the degree function of $G$, exactly as in 0D
(`docs/WITTEN_INDEX.md` §1.3). Absolute convergence: generic fibers of the
dominant polynomial map $G$ are finite and bounded by the Bézout number,
so $\int \#G^{-1}(v)\,\varphi_\sigma(v)\,dv < \infty$. (Co-area formula
and Bézout bound cited, not formalized; the direct MC estimates are
finite and stable, consistent with this.)

## 1. Exact part (all asserted symbolically in the script)

### 1.1 A1 — $M=0$ reduces exactly to the 0D MQ integral

On constant paths $q\equiv q_0$ (mode-0 coordinate $u_0 = \sqrt\beta\,
q_0$): $\dot q = 0$, the residual is constant, and
$$
G(u_0) = \sqrt\beta\,\bigl(F(u_0/\sqrt\beta) - J\bigr),\qquad
DG(u_0) = DF(q_0),\qquad |G|^2 = \beta\,|F(q_0)-J|^2,
$$
(chain-rule $\sqrt\beta$ factors cancel in $DG$ — asserted), and
$(2\pi\sigma^2)^{-3/2}\beta^{3/2} = (2\pi\sigma_{\rm eff}^2)^{-3/2}$,
$\beta/\sigma^2 = 1/\sigma_{\rm eff}^2$ with $\sigma_{\rm eff} =
\sigma/\sqrt\beta$. Hence
$$
Z_\sigma(J; 0, \beta) \;=\; Z^{0\mathrm{D}}_{\sigma/\sqrt\beta}(J)
\;=\; -\,\mathbb{E}\bigl[N(J + \sigma_{\rm eff}\,\xi)\bigr]
$$
*exactly*, and the $M=0$ baseline reproduces the 0D jump: at
$\sigma_{\rm eff}=0.03$ the closed form gives $-3.0000/-3.0000/-1.0000/
-1.0000$ at the four chamber points (asserted to $<10^{-3}$). Larger
$\beta$ = smaller effective coupling: the kinetic term *helps* the
localization, it does not fight it.

### 1.2 A2 — saddle factorization (the theorem-shaped core)

**Statement (exact for the truncated model).** Let $q^*$ be a
constant-path zero of $G$, i.e. $F(q^*)=J$ (on constant paths the
truncated flow equation reduces to the 0D fiber equation — mode 0 is the
only nonzero component of the residual). Then $DG(u^*)$ is block-diagonal
over modes: the $k=0$ block is $A = DF(q^*)$ and the $k\ge1$ block in the
$(\cos,\sin)$ basis is $\bigl[\begin{smallmatrix}A & \omega_k I\\
-\omega_k I & A\end{smallmatrix}\bigr]$. For a *generic* $3\times3$
matrix $A$ (nine symbols) and $\omega>0$:
$$
\det\begin{pmatrix}A & \omega I\\ -\omega I & A\end{pmatrix}
= \det(A^2+\omega^2 I)
= \det(A+i\omega I)\det(A-i\omega I)
\;\overset{A\ \text{real}}{=}\; (\mathrm{Re})^2 + (\mathrm{Im})^2 \ \ge\ 0,
$$
all three equalities asserted symbolically (`A2(i)–(iii)`), with equality
to $0$ iff $\pm i\omega\in\operatorname{spec}A$. Consequently
$$
\det DG(u^*) = \det DF(q^*)\cdot\!\!\prod_{k=1}^{M}\bigl|\det(DF(q^*) +
i\omega_k I)\bigr|^2,\qquad
\operatorname{sign}\det DG(u^*) = \operatorname{sign}\det DF(q^*),
$$
provided $\det(DF(q^*)^2+\omega_k^2 I)\neq0$ for all $k\le M$ (**proviso
P3**; the numeric factorization is re-asserted against the assembled
$DG$ to $3\times10^{-14}$ at all four chamber points, $\beta\in\{0.5,1,
2\}$, $M\in\{1,2\}$, and the spectral gap $\min_k|\operatorname{spec}
DF(q^*)\mp i\omega_k| = 3.14$ is asserted $>10^{-6}$).

**Corollary.** The $\sigma\to0$ limit of $Z_\sigma(J;M,\beta)$ restricted
to the constant-path saddles is
$$
\sum_{q^*:F(q^*)=J}\operatorname{sign}\det DF(q^*) = \deg(F,J) = -N(J),
$$
independent of $M$ and $\beta$: **the mode fluctuation determinant is a
$J$-independent positive factor** that the normalization removes. This is
the exact mechanism by which the truncation preserves the jump. The
statement becomes the full $\sigma\to0$ limit *provided* (P1) no
nonconstant zeros of $G$ contribute and (P2) no mass sneaks in from
infinity in mode space — the honest gaps, stated next and probed in §2.

### 1.3 A3 — gradient rigidity, and why AM escapes it

At any zero of $G$: $0 = \langle\dot q, P_M\delta q\rangle = \langle\dot
q, \delta q\rangle$ ($P_M$ self-adjoint, $\dot q\in\operatorname{ran}
P_M$), and $\langle\dot q, J\rangle = 0$ on periodic loops, so
$$
\lVert\dot q\rVert_{L^2}^2 \;=\; -\,\bigl\langle\dot q, F(q)\bigr\rangle .
$$
- **Gradient force $F=\nabla W$:** $\langle\dot q,\nabla W(q)\rangle =
  \oint dW = 0$ for *every* truncated loop (asserted symbolically for the
  tame shear, $M=1$, all nine coefficients generic), so $\dot q\equiv0$:
  **the truncated zeros are exactly the constant equilibria, for every
  $M,\beta$** — proviso P1 is a theorem on the variational branch.
- **AM force:** on the witness loop $q(s) = (\cos2\pi s,\ \sin2\pi s,\
  \tfrac12)$ the loop integral equals $75\pi/16\neq0$ (exact) — the
  rigidity proof does not apply to the non-gradient AM map (**proviso
  P1** stays open exactly; probed numerically in §2.4: 8 combos
  $\times$ 240 Newton starts, every converged zero is a constant path at
  a known equilibrium).

## 2. Numerical part (MC importance sampling)

**Method.** Deterministic seeds per run (reported), error bars = stderr
of the importance weights. Proposal: mixture of an exact "tube"
component — sample $J' = J + \sigma_{\rm eff}\xi$, pick one of the
$N(J')$ real preimages uniformly (eliminant cubic via batched companion
matrices, rational $y,z$ parametrization, Newton polish), lift to a
constant path, add higher-mode Gaussians with the *exact* block
covariance $\sigma^2(B_k^{\mathsf T}B_k)^{-1}$ — plus a $3\times$
defensive wide tube and a broad isotropic Gaussian (probes mass neither
tube sees, e.g. nonconstant zeros). All component densities are in
closed form, so the estimator is unbiased at every $\sigma$. At $M=0$
the tube proposal makes the estimator essentially exact; the end-to-end
check against the 0D closed form passes at $\sigma=0.25$ to
$<0.02$.

### 2.1 Proper-map controls ($Z \equiv -1$; fixes the normalization)

| control | $(M,\beta,\sigma)$ | $Z$ (MC) |
|---|---|---|
| linear $L$ | $(1,1,0.25)$ / $(2,0.5,0.25)$ / $(1,2,0.5)$ | $-1.000\pm0.009$, $-1.000\pm0.009$, $-1.002\pm0.009$ |
| shear $\nabla W$ | $(1,1,0.25)$ / $(2,0.5,0.25)$ / $(1,2,0.5)$ | $-1.000\pm0.009$, $-1.000\pm0.010$, $-0.998\pm0.009$ |

Also exact: $\det DG$ constant for the affine $G$ of $L$ (asserted
numerically to $10^{-9}$).

### 2.2 The main tables: $Z_\sigma(J; M, \beta)$, MC $\pm$ stderr, $[\,0\mathrm{D}$ value at $\sigma_{\rm eff}]$

$M=1$ ($n=9$), $S = 24{,}000$ samples per cell ($\times2$ at
$\sigma\ge0.25$), default seeds:

**$\beta = 0.5$:**

| $J$ | $\sigma=0.5$ | $\sigma=0.25$ | $\sigma=0.1$ | $\sigma=0.05$ |
|---|---|---|---|---|
| $(-\tfrac14,0,0)$ | $-1.908\pm0.164\ [-2.066]$ | $-2.485\pm0.035\ [-2.524]$ | $-2.923\pm0.014\ [-2.924]$ | $-3.001\pm0.014\ [-3.000]$ |
| $(0,2,0)$ | $-1.893\pm0.033\ [-1.941]$ | $-2.464\pm0.112\ [-2.380]$ | $-2.924\pm0.014\ [-2.915]$ | $-3.005\pm0.014\ [-2.999]$ |
| $(1,0,0)$ | $-0.964\pm0.161\ [-1.153]$ | $-1.002\pm0.005\ [-1.005]$ | $-1.002\pm0.005\ [-1.000]$ | $-1.000\pm0.005\ [-1.000]$ |
| $(2,1,1)$ | $-1.086\pm0.081\ [-1.003]$ | $-1.903\pm0.951\ [-1.000]$ | $-0.919\pm0.024\ [-1.000]$ | $-0.998\pm0.005\ [-1.000]$ |

**$\beta = 1$:**

| $J$ | $\sigma=0.5$ | $\sigma=0.25$ | $\sigma=0.1$ | $\sigma=0.05$ |
|---|---|---|---|---|
| $(-\tfrac14,0,0)$ | $-2.000\pm0.145\ [-2.325]$ | $-2.665\pm0.033\ [-2.690]$ | $-2.979\pm0.016\ [-2.988]$ | $-2.997\pm0.014\ [-3.000]$ |
| $(0,2,0)$ | $-2.403\pm0.318\ [-2.138]$ | $-2.573\pm0.019\ [-2.630]$ | $-2.984\pm0.015\ [-2.985]$ | $-3.002\pm0.014\ [-3.000]$ |
| $(1,0,0)$ | $-2.622\pm1.529\ [-1.048]$ | $-0.991\pm0.004\ [-1.000]$ | $-1.001\pm0.005\ [-1.000]$ | $-1.001\pm0.005\ [-1.000]$ |
| $(2,1,1)$ | $-1.111\pm0.085\ [-1.000]$ | $-1.002\pm0.004\ [-1.000]$ | $-1.001\pm0.005\ [-1.000]$ | $-0.999\pm0.005\ [-1.000]$ |

**$\beta = 2$:**

| $J$ | $\sigma=0.5$ | $\sigma=0.25$ | $\sigma=0.1$ | $\sigma=0.05$ |
|---|---|---|---|---|
| $(-\tfrac14,0,0)$ | $-2.357\pm0.173\ [-2.524]$ | $-2.647\pm0.039\ [-2.846]$ | $-2.996\pm0.023\ [-3.000]$ | $-3.001\pm0.016\ [-3.000]$ |
| $(0,2,0)$ | $-2.013\pm0.114\ [-2.380]$ | $-3.009\pm0.169\ [-2.828]$ | $-2.979\pm0.028\ [-2.999]$ | $-2.989\pm0.016\ [-3.000]$ |
| $(1,0,0)$ | $-0.925\pm0.026\ [-1.005]$ | $-0.993\pm0.058\ [-1.000]$ | $-1.002\pm0.005\ [-1.000]$ | $-1.000\pm0.005\ [-1.000]$ |
| $(2,1,1)$ | $-1.011\pm0.019\ [-1.000]$ | $-1.002\pm0.005\ [-1.000]$ | $-1.000\pm0.005\ [-1.000]$ | $-1.000\pm0.005\ [-1.000]$ |

$M=2$ ($n=15$), $\beta=1$ (all three $\beta$ with `--full`):

| $J$ | $\sigma=0.5$ | $\sigma=0.25$ | $\sigma=0.1$ | $\sigma=0.05$ |
|---|---|---|---|---|
| $(-\tfrac14,0,0)$ | $-2.304\pm0.373\ [-2.325]$ | $-2.575\pm0.032\ [-2.690]$ | $-2.981\pm0.018\ [-2.988]$ | $-3.005\pm0.015\ [-3.000]$ |
| $(0,2,0)$ | $-1.813\pm0.093\ [-2.138]$ | $-2.500\pm0.027\ [-2.630]$ | $-2.980\pm0.018\ [-2.985]$ | $-3.003\pm0.015\ [-3.000]$ |
| $(1,0,0)$ | $-0.973\pm0.021\ [-1.048]$ | $-1.022\pm0.031\ [-1.000]$ | $-1.000\pm0.005\ [-1.000]$ | $-1.001\pm0.005\ [-1.000]$ |
| $(2,1,1)$ | $-1.035\pm0.163\ [-1.000]$ | $-1.000\pm0.191\ [-1.000]$ | $-0.969\pm0.007\ [-1.000]$ | $-1.002\pm0.005\ [-1.000]$ |

**Verdict metric** $R = \bar Z(N{=}3\text{ pts})/\bar Z(N{=}1\text{
pts})$ (jump survives $\iff R\to3$):

| $(M,\beta)$ | $\sigma=0.5$ | $\sigma=0.25$ | $\sigma=0.1$ | $\sigma=0.05$ |
|---|---|---|---|---|
| $(1, 0.5)$ | $1.85\pm0.18$ | $1.70\pm0.56$ | $3.05\pm0.04$ | $3.006\pm0.014$ |
| $(1, 1)$ | $1.18\pm0.49$ | $2.63\pm0.02$ | $2.98\pm0.01$ | $2.998\pm0.014$ |
| $(1, 2)$ | $2.26\pm0.11$ | $2.84\pm0.12$ | $2.98\pm0.02$ | $2.994\pm0.015$ |
| $(2, 1)$ | $2.05\pm0.25$ | $2.51\pm0.24$ | $3.03\pm0.02$ | $2.999\pm0.015$ |

Readings, in order of importance:

1. **$\sigma\to0$: $Z\to-N(J)$ and $R\to3$ for every probed $(M,\beta)$**
   (asserted at $\sigma=0.05$, $\beta=1$, $M=1,2$: all eight cells within
   $\max(5\,\mathrm{se}, 0.05)$ of $-N$, ratio within $5\%$ of 3). The
   jump survives the truncation uniformly in the probed $M$ and $\beta$.
2. **The higher modes are near-spectators**: the values track the *0D*
   closed form at $\sigma_{\rm eff}=\sigma/\sqrt\beta$ (bracketed) to a
   few percent throughout — the content of the saddle factorization plus
   small nonlinear mode couplings. The tracking is *not* an identity for
   $M\ge1$: genuine deviations appear at moderate $\sigma$ (largest:
   $-2.647\pm0.039$ vs $[-2.846]$ at $\beta=2$, $\sigma=0.25$) and die
   off as $\sigma\to0$.
3. **A genuine finite-$\sigma$, $M$-dependent transient** (double-seeded,
   reproducible): at $J=(2,1,1)$, $M=2$, $\beta=1$, $\sigma=0.1$ the value
   dips to $-0.962\pm0.010$ / $-0.961\pm0.016$ (two seeds, $S=60{,}000$)
   and returns to $-1.002\pm0.003$ by $\sigma=0.07$. $DF(q^*)$ there has
   the strongly rotating spectrum $-3.80\pm8.45i$ (comparable to
   $\omega_{1,2}$), so higher-mode nonlinearities are maximal — a
   transient of the localization, not a shift of the limit.
4. **Heavy tails at $\sigma\ge0.25$**: a few cells (e.g. $(2,1,1)$ at
   $(M{=}1,\beta{=}0.5,\sigma{=}0.25)$: $-1.90\pm0.95$; with `--full`,
   $(0,2,0)$ at $(M{=}2,\beta{=}0.5,\sigma{=}0.5)$: $-5.07\pm3.27$) are
   dominated by rare large weights; the quoted stderr is then itself
   noisy. These cells carry no verdict weight — the claim lives at
   $\sigma\le0.1$.

**`--full` run** ($S=60{,}000$, $\sigma$ grid extended to $0.025$, all
three $\beta$ at $M=2$; 37 checks, 1340 s): every $(M,\beta)$ combo has
$Z = -N(J)\pm0.01$ and ratio $R = 2.998$–$3.006\pm0.009$ at
$\sigma=0.025$ — six independent $(M,\beta)$ confirmations of the limit.
The $(2,1,1)$, $M=2$ transient of reading 3 is $\beta$-dependent and
larger at $\beta=0.5$ ($-0.686\pm0.015$ at $\sigma=0.25$,
$-0.906\pm0.006$ at $\sigma=0.1$, $-0.998\pm0.003$ at $\sigma=0.05$):
higher-mode nonlinear corrections are real and $M$-, $\beta$-dependent
at finite $\sigma$, and vanish in the limit. Wall approach at
$\sigma=0.025$ tracks the 0D crossover to $\sim0.01$ with the $x$-far
mass fully localized at $a=-1/4$ ($10^{-22}$) and growing to
$3.4\times10^{-2}$ at $a=-1/64$; $h$-far stays $0$. Newton hunt at 480
starts per combo (3840 total): all converged zeros constant.

### 2.3 Wall approach and the far-mass probe (proviso P2)

$J = (a,0,0)$, $p = 16a\to0^-$; the escaping pair sits at $x^* =
\pm\tfrac12(-a)^{-1/2}$. $M=1$, $\beta=1$, $S=24{,}000$:

| $a$ | $x^*$ | $\sigma$ | $Z$ (MC) | $Z^{0\mathrm{D}}_{\sigma_{\rm eff}}$ | $x$-far $\lvert w\rvert$-frac | $h$-far $\lvert w\rvert$-frac |
|---|---|---|---|---|---|---|
| $-1/4$ | 1 | 0.1 / 0.05 | $-2.996\pm0.016$ / $-3.004\pm0.014$ | $-2.988$ / $-3.000$ | $4\cdot10^{-3}$ / $2\cdot10^{-6}$ | 0 / 0 |
| $-1/16$ | 2 | 0.1 / 0.05 | $-2.480\pm0.014$ / $-2.786\pm0.014$ | $-2.472$ / $-2.790$ | $2.9\cdot10^{-2}$ / $2.8\cdot10^{-2}$ | 0 / 0 |
| $-1/64$ | 4 | 0.1 / 0.05 | $-2.137\pm0.013$ / $-2.258\pm0.013$ | $-2.129$ / $-2.248$ | $1.0\cdot10^{-2}$ / $1.9\cdot10^{-2}$ | 0 / 0 |
| $-1/256$ | 8 | 0.1 / 0.05 | $-2.038\pm0.013$ / $-2.060\pm0.013$ | $-2.036$ / $-2.065$ | $4\cdot10^{-3}$ / $6\cdot10^{-3}$ | 0 / 0 |

($x$-far = $|w|$-mass with $|x(q_0)| > 2x^*+1$, i.e. beyond the known
equilibria in the escape coordinate — escape is exactly in $x$,
`scripts/branch_locus.py`; $h$-far = $|w|$-mass with nonconstant-mode
norm $>10\sigma$.) The truncated $D=1$ model reproduces the 0D crossover
quantitatively; the index degrades toward $-2$ near the wall **exactly as
in 0D**, i.e. by chamber-boundary mass, and the mass fraction living in
*nonconstant* mode directions is zero at MC resolution: **the truncation
opens no new escape channel, and it closes none** — vacua at infinity
keep their index contribution until $\sigma\gg\mathrm{dist}(J,\text{wall})_{\rm eff}$.

### 2.4 Nonconstant-zero hunt (proviso P1) and spectral gap (P3)

Damped multi-start Newton on $G$ (start scales $0.5/1.5/4.0$, 240 starts
per combo; $M=1$: $\beta\in\{0.5,1,2\}$; $M=2$: $\beta=1$; two chamber
points each): 24–72 starts converge per combo ($\lVert G\rVert<10^{-10}$)
and **every** converged zero is a constant path at a known equilibrium
(higher-mode norm $<10^{-7}$, matched to the exact fiber to $10^{-6}$).
No nonconstant zero of the truncated flow was found — a probe, not a
proof (A3). Spectral proviso P3: $\min|\operatorname{spec}DF(q^*)\mp
i\omega_k| = 3.14$ over all equilibria, $\beta$ in the grid, $k\le2$.

## 3. Q1 verdict (truncated model)

> **The index jump survives.** In the finite-mode truncation of the $D=1$
> MQ model, $\lim_{\sigma\to0}Z_\sigma(J;M,\beta) = \deg(F,J) = -N(J)$
> at all probed chamber points, for $M\in\{1,2\}$ and
> $\beta\in\{0.5,1,2\}$, with the exact saddle factorization (§1.2)
> identifying the mechanism: the path measure's mode determinant is a
> $J$-independent positive factor. $\lim_{\sigma\to0}Z(J_+)\ne
> \lim_{\sigma\to0}Z(J_-)$, ratio $3$ — the "yes (survive)" branch of Q1,
> **for the truncation**. The kinetic term only renormalizes the coupling
> ($\sigma\mapsto\sigma/\sqrt\beta$ at $M=0$, spectator determinants at
> $M\ge1$); it does not suppress the equilibria escaping to infinity.

In the §5.3 honesty checklist of `docs/CLASSICAL_MAP_INVARIANTS.md`, this
documents "positive C7 in $D=1$ mechanics" *at finite mode truncation*;
it licenses nothing about the continuum, $D\ge2$, or UV.

## 4. Honest limits

- **Finite modes $\neq$ continuum.** Everything here is a statement about
  the family of finite-dimensional integrals indexed by $(M,\beta)$. The
  $\sigma\to0$ limit is taken at *fixed* $M$; the order of limits
  $M\to\infty$ vs $\sigma\to0$ is exactly where a continuum measure could
  still kill the jump ("uniformly as $M$ grows" is probed only at
  $M\le2$, where no $M$-trend is visible in the limit values — but the
  finite-$\sigma$ transient of §2.2(3) *is* $M$-dependent, so uniformity
  in $\sigma$ is not automatic). No continuum claim.
- **Proviso P1 (nonconstant zeros)** is a theorem only on the gradient
  branch (§1.3). For AM it rests on a 1920-start Newton probe at 8
  parameter combos. A nonconstant zero elsewhere in $(J,\beta,M)$-space
  would add its own $\operatorname{sign}\det DG$ to the limit.
- **Proviso P2 (mass from infinity)** is probed by the wide-tube/broad
  mixture components and the far-mass diagnostics, not bounded
  analytically. The 0D experience (`docs/WITTEN_INDEX.md` §1.3) is that
  this mass is precisely the wall crossover; the truncated model matches
  that pattern at MC resolution.
- **MC error bars** are stderr of importance weights with deterministic
  seeds; at $\sigma\ge0.25$ some cells have heavy-tailed weights and the
  stderr itself is noisy (flagged in §2.2). Verdict-carrying cells
  ($\sigma\le0.1$) have relative errors $\lesssim0.5\%$ and pass
  5-sigma asserts.
- **Exactness perimeter.** Symbolic asserts cover: the $M=0$ reduction,
  the block-determinant factorization for generic $3\times3$ $A$, the
  gradient-rigidity loop integral at $M=1$ (generic coefficients), and
  the AM witness. The co-area/Bézout finiteness frame and the standard
  Gaussian-localization limit at nondegenerate zeros are cited, not
  formalized. $N_\tau$-exactness of the mode projection is asserted
  numerically ($64$ vs $512$ points, $3\times10^{-14}$).
- The wall-approach comparison asserts agreement with the 0D crossover to
  $5\,\mathrm{se}+0.1$ — a tracking statement, not an exact identity;
  genuine $M\ge1$ deviations at moderate $\sigma$ are expected and seen
  (§2.2(3)).

## 5. Reproduce

```bash
.venv/bin/python scripts/d1_index_modes.py          # ~4.5 min, 37 checks
.venv/bin/python scripts/d1_index_modes.py --full   # ~25 min: sigma=0.025, all beta for M=2, S=60k
```

## References

- [Wit82] E. Witten, *Constraints on supersymmetry breaking*, Nucl. Phys.
  B **202** (1982) 253–316 (index jumping via vacua at infinity).
- [MQ86] V. Mathai, D. Quillen, *Superconnections, Thom classes, and
  equivariant differential forms*, Topology **25** (1986) 85–110.
- [Bla93] M. Blau, *The Mathai–Quillen formalism and topological field
  theory*, J. Geom. Phys. **11** (1993) 95–127.
- [Nic80] H. Nicolai, *On a new characterization of scalar supersymmetric
  theories*, Phys. Lett. B **89** (1980) 341 (Nicolai map; bibliographic
  data best-effort from memory).
- Repository: `scripts/d1_index_modes.py` (this note),
  `scripts/witten_index.py` / `docs/WITTEN_INDEX.md` (0D index, closed
  form), `docs/DAMPED_PARTITION.md` (measure layer),
  `docs/CLASSICAL_MAP_INVARIANTS.md` §5.1 (Q1), `jcqft/fibers.py` (exact
  fibers / eliminant parametrization).
