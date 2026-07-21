# Twisted Cohomology of the Wall Complement: Exact Dimensions, Jump Loci, and the Fate of the Period Pairing

*(2026-07-21. All exact claims verified in `scripts/twisted_cohomology.py`
(exact sympy — rationals and cyclotomics, symbolic `t` where possible; runs
in ~3 s); this answers item **Q3** of `docs/WALL_COMPLEMENT.md` §6, the
concrete follow-up to resolved question B2 of `docs/OPEN_QUESTIONS.md`.
Repo-internal inputs: the $K(B_3,1)$ structure and affine $A_2$-equivalence
of `docs/WALL_COMPLEMENT.md`, the failed chamber canonical form of
`docs/POSITIVE_GEOMETRY.md`, the $S_3$ monodromy of `docs/MONODROMY.md`.)*

**Summary of the main results.** Throughout, $M = \mathbb{C}^2 \setminus
\{P_2 = 0\}$ is the wall complement, a $K(B_3,1)$, so twisted cohomology of
$M$ with any local system $\rho$ *is* group cohomology $H^*(B_3;\rho)$,
computed exactly by Fox calculus on the 1-relator presentation
$B_3 = \langle \sigma_1,\sigma_2 \mid \sigma_1\sigma_2\sigma_1 =
\sigma_2\sigma_1\sigma_2\rangle$.

1. **Complete dimension table** (exact; $(h^0,h^1,h^2)$, and
   $\chi = h^0 - h^1 + h^2 = 0$ in every case, as $\chi(M)=0$ forces):

   | local system | $(h^0,h^1,h^2)$ | jump polynomial |
   |---|---|---|
   | Kummer $\mathbb{C}_t$ ($P_2^s$, $t = e^{2\pi i s}$), generic $t$ | $(0,0,0)$ | — |
   | $\mathbb{C}_t$, $t = 1$ (i.e. $s \in \mathbb{Z}$: untwisted) | $(1,1,0)$ | $t-1$ |
   | $\mathbb{C}_t$, $\Delta(t)=0$ (i.e. $s \in \pm\tfrac16 + \mathbb{Z}$) | $(0,1,1)$ | $\Delta(t)=t^2-t+1$ |
   | sign $= \mathbb{C}_{-1}$ ($s = \tfrac12$) | $(0,0,0)$ | — |
   | **standard (reflection) $\mathcal{S}$** | $(0,1,1)$ | — |
   | permutation (sheet) $\mathcal{L}$ | $(1,2,1)$ | — |
   | regular $\mathbb{C}[S_3]$ | $(1,3,2)$ | — |
   | reduced Burau $\rho_t$, generic $t$ | $(0,0,0)$ | — |
   | Burau, $t = 1$ ($=$ standard rep, on the nose) | $(0,1,1)$ | $t^3-1$ |
   | Burau, $t$ a primitive cube root | $(1,2,1)$ | $t^3-1$ |

   The Kummer jump locus in $s \in [0,1)$ is exactly
   $\{0\} \cup \{\tfrac16, \tfrac56\}$: the integers (untwisted) and the
   arguments of the primitive 6th roots of unity — the roots of the
   **Alexander polynomial of the trefoil** $\Delta(t) = t^2 - t + 1$, which
   appears verbatim as the Fox differential $d^1 = \Delta(t)\,(1,-1)$.
2. **Two independent methods agree.** Fox calculus and the **Wang sequence**
   of the global Milnor fibration $P_2 : M \to \mathbb{C}^*$ of the cusp
   ($\mu = 2$, monodromy of order 6 with the primitive 6th roots as
   eigenvalues, by the Brieskorn formula for $x^2 + y^3$) give the *same*
   dimensions at every tested $t$, and the same jump polynomials
   *identically*: $\det(t\,h^* - 1 \mid H^1 F) = \Delta(t)$,
   $\det(t - 1 \mid H^0 F) = t - 1$.
3. **The proposed "amplituhedron-analogue" pair carries exactly one class
   per degree**: $\dim H^1(M;\mathcal{S}) = \dim H^2(M;\mathcal{S}) = 1$ —
   the reflection channel of `docs/WALL_COMPLEMENT.md` §5 has a
   one-dimensional space of twisted "master integrals". Additivity
   $\mathcal{L} = \mathrm{triv} \oplus \mathcal{S}$ is verified on
   cohomology: $(1,2,1) = (1,1,0) + (0,1,1)$. Shapiro cross-check:
   $H^*(B_3;\mathbb{C}[S_3]) = (1,3,2) = H^*(\text{pure braid space})$,
   matching an independent Orlik–Solomon/Möbius computation.
4. **Item 5 answered (canonical-form connection):** at $s \in \mathbb{Z}$
   (integer twist, $t=1$) the period pairing in degree 2 degenerates
   *completely*: the would-be canonical form is globally twisted-exact,
   $\dfrac{dU \wedge dW}{f} = d\left(-\dfrac{2U\,dW - 3W\,dU}{f}\right)$
   ($f = W^2 - U^3 = P_2/27$ in the $A_2$ coordinates) — every period
   vanishes; $h^2 = 0$. This is the global cohomological face of the
   residueless double pole of `docs/POSITIVE_GEOMETRY.md` §3. The two
   $s$-values where $[dU{\wedge}dW/f]$ and $[U\,dU{\wedge}dW/f]$ become
   honest nonzero classes are    $s \equiv \tfrac16$ and $s \equiv \tfrac56$ — the singularity
   spectrum of the cusp modulo $\mathbb{Z}$
   ($e^{2\pi is}$ = the Milnor monodromy eigenvalues), realized on the
   Milnor-ring basis $\{1, U\}$. At generic $s$ the intersection theory is
   *vacuously* perfect: there are **zero** twisted classes
   ($|\chi(M)| = 0$), so this model has no generic-twist master integrals
   at all; everything interesting is resonant.

**Honest deviation from the task's expectations:** reduced Burau at
$t = -1$ is **not** the reflection representation — it is the integral
$SL(2,\mathbb{Z})$ representation ($\sigma_i \mapsto$ unipotent, infinite
order; the full twist $\mapsto -I \neq I$, so it does not factor through
$S_3$). The reflection representation is Burau at $t = +1$ (exact matrix
equality). Accordingly the Burau jump locus is $\{t^3 = 1\}$, not
$\{\Delta(t) = 0\}$; the Alexander polynomial enters as
$t^2 + t + 1 = \Delta(-t)$, through the Burau sign convention
$t_{\text{Burau}} = -t_{\text{Kummer}}$ (visible in
$\det \rho_t(\sigma_i) = -t$).

---

## 0. Method: Fox calculus on the 1-relator group, and why it computes $H^*(M)$

$M$ is a $K(B_3,1)$ (`docs/WALL_COMPLEMENT.md`; Arnold–Brieskorn–Deligne),
so $H^*(M;\rho) = H^*(B_3;\rho)$ for every local system $\rho$ — this is
the definition-level reduction, no comparison theorem needed. $B_3$ is a
2-generator 1-relator group with relator
$r = xyx\,y^{-1}x^{-1}y^{-1}$ not a proper power, so by Lyndon's theorem
the presentation 2-complex is aspherical and the free (Fox) resolution has
length 2. For $\rho: B_3 \to GL(V)$ the cochain complex is

$$
0 \to V \xrightarrow{\;d^0\;} V^2 \xrightarrow{\;d^1\;} V \to 0,
\qquad
d^0 v = \bigl((\rho(x)-1)v,\ (\rho(y)-1)v\bigr),
\qquad
d^1(v_1,v_2) = \rho\!\left(\tfrac{\partial r}{\partial x}\right) v_1
             + \rho\!\left(\tfrac{\partial r}{\partial y}\right) v_2 ,
$$

with Fox derivatives (computed from the product rule, then simplified with
the braid relation — the script verifies both forms agree):

$$
\frac{\partial r}{\partial x} = 1 + xy - xyxy^{-1}x^{-1} = 1 + xy - y,
\qquad
\frac{\partial r}{\partial y} = x - xyxy^{-1} - 1 = x - yx - 1 .
$$

The fundamental identity $d^1 d^0 = \rho(r) - 1 = 0$ is asserted for every
representation used. Note that on this shape of complex
$\chi = h^0 - h^1 + h^2 = \dim V - 2\dim V + \dim V = 0$ *identically* (by
rank–nullity, independent of $\rho$) — consistent with
$\chi(M) = 0$ (`docs/WALL_COMPLEMENT.md` §4) but carrying no information;
the content of this document is the individual dimensions. Untwisted sanity
check: the trivial system gives $(1,1,0) = H^*(B_3;\mathbb{Z}) \otimes
\mathbb{C}$ [Arn69].

All ranks are computed exactly (sympy, exact rationals and cyclotomics;
symbolic $t$ wherever a statement is claimed for all $t$).

## 1. The rank-1 Kummer twist $P_2^s$: the Alexander polynomial governs the jumps

The multiplicative twist by $P_2^s$ is the character
$\lambda_t : B_3 \to \mathbb{C}^*$, $\sigma_1, \sigma_2 \mapsto t$
($t = e^{2\pi i s}$: both generators are wall meridians and the twist has
meridian monodromy $t$; $\lambda_t$ is the composite of abelianization
$B_3 \twoheadrightarrow \mathbb{Z}$ with $1 \mapsto t$). Fox calculus gives,
symbolically in $t$:

$$
d^0 = \begin{pmatrix} t-1 \\ t-1 \end{pmatrix},
\qquad
d^1 = \Delta(t)\,\bigl(1,\ -1\bigr),
\qquad
\Delta(t) = t^2 - t + 1 = \Phi_6(t).
$$

$\Delta$ is the Alexander polynomial of the trefoil — as it must be: $M$ is
homotopy equivalent to the trefoil complement ($B_3$ is the trefoil group),
and $d^1$ presents the Alexander module. Exact dimension table (each row an
assert in the script, §2):

| $t$ | $s \in [0,1)$ | $(h^0, h^1, h^2)$ |
|---|---|---|
| $1$ | $0$ | $(1, 1, 0)$ — untwisted $H^*(M)$ |
| $\zeta_6,\ \zeta_6^5$ | $\tfrac16,\ \tfrac56$ | $(0, 1, 1)$ |
| $-1$ | $\tfrac12$ | $(0,0,0)$ |
| $\zeta_3,\ \zeta_3^2$ | $\tfrac13,\ \tfrac23$ | $(0,0,0)$ |
| $2,\ \tfrac75,\ i$ | (generic) | $(0,0,0)$ |

**Jump locus in $s$: exactly $s \in \mathbb{Z}$ (dimensions $(1,1,0)$) and
$s \in \pm\tfrac16 + \mathbb{Z}$ (dimensions $(0,1,1)$, i.e. one class in
$H^1$ and one in $H^2$); at every other $s$ the twisted cohomology
vanishes identically.** A generic Kummer twist kills everything — the
generic number of twisted classes equals $|\chi(M)| = 0$.

## 2. Independent cross-check: the Milnor fibration and the Wang sequence

Because $P_2 = 27(W^2 - U^3)$ in the affine coordinates (I4) of
`docs/WALL_COMPLEMENT.md` (exact identity, re-verified), $P_2$ is affinely
a quasi-homogeneous polynomial, and $P_2 : M \to \mathbb{C}^*$ is the
**global Milnor fibration** of the cusp $x^2 + y^3$ [Mil68, §9]. $M$ is
the mapping torus of the geometric monodromy $h$ of the Milnor fiber $F$,
with $(b_0, b_1)(F) = (1, \mu) = (1, 2)$ and $h^*$ of order 6 on
$H^1(F)$ with eigenvalues $\{-\zeta_3, -\zeta_3^2\}$ = the primitive 6th
roots of unity (Brieskorn/Pham formula for $x^a + y^b$:
eigenvalues $\zeta_a^i \zeta_b^j$, $1\le i\le a-1$, $1\le j \le b-1$). The
Kummer system $\mathbb{C}_t$ is pulled back from $\mathbb{C}^*$ (wall
meridians map to degree-1 loops), so the Wang sequence gives

$$
0 \to \operatorname{coker}\bigl(t\,h^* - 1 \mid H^{k-1}(F)\bigr)
  \to H^k(M;\mathbb{C}_t)
  \to \ker\bigl(t\,h^* - 1 \mid H^{k}(F)\bigr) \to 0 .
$$

Realizing $h^*$ on $H^1(F)$ exactly as the companion matrix of $\Delta$
(order 6, $h^{*3} = -1$; verified), the script computes the Wang
dimensions at all nine test values of $t$ and finds **agreement with Fox
calculus at every point**; moreover the jump determinants agree
*identically as polynomials in $t$*:

$$
\det\bigl(t\,h^* - 1 \mid H^1 F\bigr) = t^2 - t + 1 = \Delta(t),
\qquad
\det\bigl(t - 1 \mid H^0 F\bigr) = t - 1 .
$$

So the two methods do not just agree numerically — they agree at the level
of the polynomials whose vanishing controls the jumps. (This is the
requested second, independent computation; it is topological, uses no group
presentation, and imports only the Milnor fibration and the Brieskorn
eigenvalues.)

## 3. The $S_3$ local systems: the reflection channel has exactly one class per degree

With the exact reflection matrices $r_1 = \begin{pmatrix} -1&1\\0&1
\end{pmatrix}$, $r_2 = \begin{pmatrix} 1&0\\1&-1 \end{pmatrix}$ (braid
relation, $r_i^2 = 1$, $(r_1r_2)^3 = 1$ all asserted — the representation
factors through $B_3 \twoheadrightarrow S_3$):

$$
H^*(M; \mathcal{S}) = (0,\ 1,\ 1) .
$$

**The standard local system — the proposed amplituhedron-analogue carrier
of `docs/WALL_COMPLEMENT.md` §5 — supports exactly one twisted class in
$H^1$ and one in $H^2$.** Supporting computations, each exact:

- **Sign system** ($\sigma_i \mapsto -1$): $(0,0,0)$ — coincides with the
  Kummer case $t = -1$, as it must ($\Delta(-1) = 3 \neq 0$).
- **Permutation (sheet) system $\mathcal{L}$** (rank 3): $(1,2,1)$, and
  additivity under $\mathcal{L} = \mathrm{triv} \oplus \mathcal{S}$ holds
  componentwise: $(1,2,1) = (1,1,0) + (0,1,1)$. The invariant part
  $(1,1,0)$ is the trace-observable sector; everything else is the
  reflection channel.
- **Regular system $\mathbb{C}[S_3]$** (rank 6): $(1,3,2)$. Shapiro's
  lemma says this must equal $H^*$ of the index-6 subgroup $P_3$ (the pure
  braid group) with trivial coefficients, i.e. the Betti numbers of the
  $A_2$ braid arrangement complement; an independent
  Orlik–Solomon/Möbius-function computation (three concurrent transverse
  lines: $\mu(\text{origin flat}) = 2$; Poincaré polynomial
  $(1+q)(1+2q)$) gives $(1,3,2)$. Match. Character bookkeeping
  $\mathbb{C}[S_3] = \mathrm{triv} \oplus \mathrm{sign} \oplus 2\,
  \mathcal{S}$ also matches: $(1,1,0) + (0,0,0) + 2(0,1,1) = (1,3,2)$.

## 4. The reduced Burau family: jumps on $t^3 = 1$, and which specializations recover what

For $\rho_t(\sigma_1) = \begin{pmatrix} -t&1\\0&1 \end{pmatrix}$,
$\rho_t(\sigma_2) = \begin{pmatrix} 1&0\\t&-t \end{pmatrix}$ (reduced
Burau of $B_3$ [Bir74, KT08]; braid relation asserted symbolically for all
$t$; $t = 0$ excluded since $\det\rho_t(\sigma_i) = -t$):

$$
d^1 = \begin{pmatrix} 0 & -t & -1 & 0 \\ 0 & 1 & t^2 & 0 \end{pmatrix},
\qquad
\gcd\bigl(2{\times}2 \text{ minors of } d^1\bigr) = t^3 - 1,
\qquad
\gcd\bigl(2{\times}2 \text{ minors of } d^0\bigr) = t^2 + t + 1 .
$$

Neither differential ever vanishes identically, so the **jump locus is
exactly $\{t : t^3 = 1\}$**, with exact dimensions

| $t$ | $(h^0,h^1,h^2)$ | structure |
|---|---|---|
| generic (incl. $-1$, $\zeta_6$, $2$, $\tfrac12$, $\tfrac75$) | $(0,0,0)$ | irreducible, non-resonant |
| $1$ | $(0,1,1)$ | $\rho_1 = $ reflection rep (matrix equality) |
| $\zeta_3^{\pm 1}$ | $(1,2,1)$ | reducible: $0 \to \mathbb{C}_{\mathrm{triv}} \to \rho_t \to \mathbb{C}_{-t} \to 0$ |

At $t = \zeta_3$ the fixed vector is $v = (1,\ 1+\zeta_3)$ (asserted), and
the quotient character is $\lambda_{-\zeta_3}$ with $-\zeta_3$ a
**primitive 6th root** — a root of $\Delta$. The long exact sequence stacks
$(1,1,0)$ (trivial sub) and $(0,1,1)$ (resonant Kummer quotient) to
$(1,2,1)$, exactly the measured dimensions.

Specializations, precisely:

- **$t = 1$ recovers case 3 (standard rep)** — not just isomorphic:
  the matrices are equal.
- **$\det \circ \rho_t = \lambda_{-t}$ recovers case 1** up to the Burau
  sign convention $t_{\text{Burau}} = -t_{\text{Kummer}}$; this is also why
  the Alexander polynomial appears in the Burau jump data as
  $t^2 + t + 1 = \Delta(-t)$, a factor of $t^3 - 1$.
- **$t = -1$ is *not* the reflection rep** (correcting the guess in the
  task prompt): $\rho_{-1}$ is the integral $SL(2,\mathbb{Z})$
  representation — $\sigma_i \mapsto$ unipotent matrices of infinite
  order, full twist $(\sigma_1\sigma_2)^3 \mapsto -I \neq I$, so the
  center acts nontrivially and $\rho_{-1}$ does not factor through $S_3$.
  It is the homological monodromy representation of the torus-knot
  picture; its twisted cohomology vanishes, $(0,0,0)$.

Duality sanity check (asserted for both families): dimensions at $t$ and
$t^{-1}$ (dual local system) coincide, and the jump loci are closed under
$t \mapsto t^{-1}$, as Poincaré–Lefschetz duality
$H^k(M;\mathcal{L}) \cong H^{4-k}_c(M;\mathcal{L}^\vee)^*$ requires.

## 5. The de Rham side and the failed canonical form: where the period pairing degenerates

Work in the $A_2$ coordinates $f = W^2 - U^3$ ($= P_2/27$ under (I4);
weights $(U,W) = (2,3)$, Euler field $E = 2U\partial_U + 3W\partial_W$,
$E(f) = 6f$), with the twisted differential $\nabla_s = d + s\,d\!\log f$
on rational forms with poles on the wall. All of the following are exact
one-line identities, asserted in the script (§6):

$$
\nabla_s(1) = s\,\frac{df}{f},
\qquad
\nabla_s\,\eta = (6s-1)\,\frac{dU\wedge dW}{f},
\qquad
\nabla_s(U\eta) = (6s+1)\,\frac{U\,dU\wedge dW}{f},
\qquad
\eta = \frac{\iota_E(dU\wedge dW)}{f} = \frac{2U\,dW - 3W\,dU}{f},
$$

and in general, for $g = U^aW^b$ of weight $k = 2a+3b$ (verified on the
grid $a,b \le 2$, $m \le 3$):

$$
\nabla_s\!\left(\frac{g\,\iota_E(dU\wedge dW)}{f^m}\right)
= \bigl(k + 5 - 6m + 6s\bigr)\,\frac{g\,dU\wedge dW}{f^m} .
$$

Consequences, and the precise answer to item 5:

1. **At $s \in \mathbb{Z}$ ($t=1$) the degree-2 period pairing degenerates
   completely — verified, not analogy.** Setting $s = 0$: the would-be
   canonical form of the chamber is globally exact with an explicit
   rational (log-pole) primitive,
   $\dfrac{dU\wedge dW}{f} = d(-\eta)$, so *every* period of it over
   *every* closed 2-cycle vanishes, and indeed $h^2(t{=}1) = 0$. This is
   the global, cohomological restatement of the residueless double pole
   that killed the chamber canonical form
   (`docs/POSITIVE_GEOMETRY.md` §3): there the residue of
   $du\wedge dw/P_2$ on the normalization was $-\tfrac{\kappa}{3}
   \tfrac{dm}{(m-3)^2}$, a double pole with zero residue; here the same
   form is exhibited as a twisted coboundary. What survives at $t = 1$ is
   only $h^1 = 1$, spanned by $[d\!\log f]$ — which is $\nabla_s$-exact
   ($= \nabla_s(1/s)$) the moment $s$ leaves $\mathbb{Z}$.
2. **The classes come back exactly at the spectral exponents.** Pole-order
   reduction by the Euler chain fails only when the coefficient
   $k + 5 - 6m + 6s$ vanishes; on the Milnor-ring basis $\{1, U\}$ of the
   cusp ($\mathbb{C}[U,W]/(f_U, f_W)$, $\mu = 2$; asserted) this happens
   at $s \in \tfrac16 + \mathbb{Z}$ (class $[dU{\wedge}dW/f]$) and
   $s \in \tfrac56 + \mathbb{Z}$ (class $[U\,dU{\wedge}dW/f]$) — matching
   $e^{2\pi i s} = $ primitive 6th roots, i.e. the topological jump locus
   of §§1–2. These resonant $s$-values agree modulo $\mathbb{Z}$ with the
   **singularity spectrum of the cusp** $x^2 + y^3$, which is
   $\{\tfrac56, \tfrac76\}$ (Steenbrink; for $x^a + y^b$ the spectral
   numbers are $\tfrac{i}{a} + \tfrac{j}{b}$) — equivalently,
   $e^{2\pi i s}$ runs over the eigenvalues of the Milnor monodromy.
   Non-basis numerators create no jumps — e.g. $W\,dU{\wedge}dW/f^2$ has Euler coefficient $6s-4$
   (failing at $s = \tfrac23$) but is nonetheless exact for all
   $s \neq 1$ via the non-Euler primitive $\nabla_s(dU/f) =
   2(1-s)\,W\,dU{\wedge}dW/f^2$ (asserted); and indeed §1 found *no* jump
   at $s = \tfrac23$.
3. **At generic $s$ the intersection theory is vacuously perfect.**
   $\dim H^1 = \dim H^2$ holds at every $t$ (structural $\chi = 0$;
   asserted case by case), and for generic $s$ both are **zero**: the
   intersection matrix is the empty matrix, "nondegenerate" in the trivial
   sense, and the model has **no generic-twist master integrals at all**.
   This is the honest form of the Aomoto–Kita genericity mechanism here
   [AK11; MM19]: in that theory the generic number of master integrals is
   $|\chi|$, and our $\chi = 0$. The interesting, nonempty pairings are
   the *resonant* ones — $(0,1,1)$ at $s \equiv \pm\tfrac16$ and for the
   reflection system $\mathcal{S}$ — each a $1\times 1$ pairing between a
   line in $H^1$ (or $H^2$) and a line of twisted cycles.

**Proven vs proposed in this section, explicitly.** Proven (exact,
asserted): the exactness identities, the Euler-chain grid formula, the
Milnor-ring basis, the location of the Euler-chain failures, the vanishing
of all periods of the canonical form at $s = 0$, and the match of the
resonant $s$-values with the topological jump locus. Cited, not re-proved:
that twisted algebraic de Rham cohomology of $(M, \nabla_s)$ computes
$H^*(M;\mathbb{C}_t)$ (Deligne's comparison; Esnault–Schechtman–Viehweg
non-resonance theory). Not computed: actual intersection *numbers* of
regularized twisted cycles (the matrices here are $0\times 0$ or
$1\times 1$, so nondegeneracy is forced by dimension count and duality
whenever the entry is nonzero; exhibiting the entry — e.g. by a
Mastrolia–Mizera-style residue computation at the resonant $s$ — is left
as new question 2 below).

## 6. What this buys the amplitudes analogy

- **The master-integral count of the model is now a theorem-level
  statement.** In the intersection-theory formulation of amplitudes
  [MM19], the number of master integrals is $\dim H^n_{\text{twisted}}$,
  generically $|\chi|$. For this counterexample-QFT the count is **0 at
  generic twist and 1 in each of the three special channels** (untwisted
  $H^1$; resonant Kummer $s \equiv \pm\tfrac16$; reflection system
  $\mathcal{S}$). The entire "amplitude content" of the wall complement is
  resonance phenomena — there is no generic-position sector.
- **The failed canonical form and the twisted-period proposal are now one
  statement.** The chamber form failed by a residueless double pole
  (`docs/POSITIVE_GEOMETRY.md`); we now see this as: the canonical class
  is twisted-exact at integer $s$, and becomes a *bona fide* period
  precisely when the Kummer exponent hits the spectral exponents
  $s \equiv \pm\tfrac16 \pmod{\mathbb{Z}}$ of the cusp. "Turn on the regulator $P_2^s$" is
  the exact mechanism that resurrects the dead canonical form — the
  standard dimensional-regularization-like move of twisted period theory,
  realized verbatim in this model.
- **The reflection channel is as small as it could be without dying:**
  one class in $H^1$, one in $H^2$. If a "wall amplitude" exists in this
  model, it is a single number (per degree) — a sharp, falsifiable target
  for any future pairing computation.

## 7. New questions

1. **Compute the resonant periods themselves.** At $s = \tfrac16$: the
   pairing of $[dU{\wedge}dW/f^{\,1-s}]$-type classes with the (unique up
   to scale) twisted 2-cycle — e.g. via the Milnor-fiber realization,
   where the period should reduce to Beta-function values (the cusp is
   Brieskorn–Pham, its periods are classical). A one-afternoon computation
   that would turn the $1\times 1$ "intersection matrix is nonzero" claim
   into an exact number.
2. **Intersection numbers à la Mastrolia–Mizera.** Implement the
   1-dimensional intersection-number recursion on the fibration
   $P_2: M \to \mathbb{C}^*$ (fiber $= 3$-punctured plane after root
   normalization?) and verify nondegeneracy of the resonant pairings
   directly.
3. **The reflection channel vs the resonant Kummer channel.** Both have
   $(0,1,1)$. Is there a preferred isomorphism (e.g. is $\mathcal{S}$
   cohomologically equivalent to a direct image from a Kummer system on a
   double cover, branched over the wall)? The $A_2$ picture suggests:
   $\mathcal{S} \otimes \mathbb{C}(\pm\tfrac16$-twist$)$ vs the rank-2
   pieces of pushforwards from the 3-sheeted cover.
4. **Degree-$d$ scaling.** For a hypothetical degree-4 counterexample
   (`docs/OPEN_QUESTIONS.md` C4), the same machinery predicts: jump locus
   = roots of the $A_3$-discriminant Alexander-type polynomials, spectrum
   $\{\tfrac{k+5}{6}\}$ replaced by the $A_3$ spectrum, and generic count
   still $|\chi|$. Concrete falsifiable predictions the moment such an
   example exists.
5. **The $h^2$ class at resonance as a "wall amplitude".** Does the
   resonant period of §7 Q1 have a QFT meaning in the 0D model — e.g. as
   the $s$-regularized version of the (divergent) chamber volume
   $\int_{\text{chamber}} P_2^s\,du\,dw$, whose analytic continuation in
   $s$ should have poles exactly at the spectral exponents? (This is the
   classical Mellin/Bernstein picture; the Bernstein–Sato roots of
   $W^2 - U^3$ are $-\tfrac56, -1, -\tfrac76$, matching the spectrum —
   *stated from the literature, not verified in the repo*.)

## 8. References

- [AK11] K. Aomoto, M. Kita, *Theory of Hypergeometric Functions*,
  Springer Monographs in Mathematics, 2011. (Twisted (co)homology,
  regularized cycles, intersection theory; the genericity mechanism of
  §5.3.)
- [MM19] P. Mastrolia, S. Mizera, *Feynman integrals and intersection
  theory*, JHEP **02** (2019) 139. arXiv:1810.03818. (Master-integral
  counts as twisted Betti numbers; intersection-number formulation of
  amplitudes.)
- [Mil68] J. Milnor, *Singular Points of Complex Hypersurfaces*, Ann. of
  Math. Studies **61**, Princeton, 1968. (Milnor fibration — §9 for the
  global fibration of (quasi-)homogeneous polynomials; Wang sequence;
  $\mu = (a-1)(b-1)$ and monodromy eigenvalues for $x^a + y^b$ via the
  Brieskorn–Pham description.)
- [Bir74] J. Birman, *Braids, Links, and Mapping Class Groups*, Ann. of
  Math. Studies **82**, Princeton, 1974. (Burau representation;
  $B_3$ structure.)
- [KT08] C. Kassel, V. Turaev, *Braid Groups*, GTM **247**, Springer,
  2008. (Reduced Burau conventions; reducibility at roots of unity;
  Burau(1) = reflection rep of the symmetric group.)
- [Arn69] V. I. Arnold, *The cohomology ring of the colored braid group*,
  Mat. Zametki **5** (1969). ($H^*(B_3) = (\mathbb{Z},\mathbb{Z},0)$;
  pure braid Poincaré polynomial $(1+q)(1+2q)$ used in the Shapiro
  cross-check.)
- Lyndon's theorem (aspherical presentation of 1-relator groups with
  non-power relator): R. Lyndon, *Cohomology theory of groups with a
  single defining relation*, Ann. of Math. **52** (1950) 650–665.
- Deligne comparison / non-resonance: P. Deligne, *Équations
  différentielles à points singuliers réguliers*, LNM 163, 1970;
  H. Esnault, V. Schechtman, E. Viehweg, *Cohomology of local systems on
  the complement of hyperplanes*, Invent. Math. **109** (1992) 557–561.
  (Cited for the de Rham dictionary of §5; not re-proved.)
- Repo-internal: `docs/WALL_COMPLEMENT.md` ($K(B_3,1)$, affine
  $A_2$-equivalence (I4), $\chi = 0$, the §6 Q3 prompt),
  `docs/POSITIVE_GEOMETRY.md` (failed chamber form, residueless double
  pole), `docs/MONODROMY.md` ($S_3$ sheet monodromy),
  `scripts/twisted_cohomology.py` (all verifications for this document).
