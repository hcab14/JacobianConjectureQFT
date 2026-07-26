# The Alpöge–Mathew Counterexample as a Zero-Dimensional Field Theory: What It Does and Does Not Teach Us About Rigorous QFT

**Status.** Working notes accompanying the computations in this repository
(`scripts/verify_counterexample.py`, `scripts/tree_expansion.py`, `scripts/branch_locus.py`). All
mathematical claims about the specific map below have been verified
symbolically or numerically in these scripts; statements about the literature
carry references checked against the published record.

## Abstract

The polynomial map $F:\mathbb{C}^3\to\mathbb{C}^3$ announced by L. Alpöge on
July 19, 2026 (found with the LLM Claude Fable 5, in response to a question of
Akhil Mathew) has constant Jacobian determinant $-2$ but is not injective,
disproving the Jacobian conjecture in dimension $n\ge 3$. Read as a
zero-dimensional scalar field theory with three field components, $F(\phi)=J$
is the classical field equation with external source $J$, and the standard
tree-graph (perturbative) inversion computes a formal inverse. We show, by
exact computation in this repository, that this perturbation series is a
*convergent* series with finite nonzero radius of convergence — the Taylor
expansion of an explicit degree-3 algebraic function — and that what fails is
not summability but *globality*: two of the three solution branches sit at
infinity in field space over the perturbative vacuum $J=0$, recede there along
an explicit curve, and are invisible to every order of perturbation theory.
The obstruction is precisely non-properness of $F$ in the sense of Jelonek,
not divergence of the expansion. We correct several errors in a
widely-circulated AI-generated commentary on the result, and we give an honest
assessment of what this toy-model phenomenon does and does not imply for
rigorous quantum field theory in four dimensions, in the spirit of
constructive QFT and perturbative algebraic QFT.

---

## 1. Setup: the map, the dictionary, the counterexample

### 1.1 The counterexample

The map is

$$
F(x,y,z) \;=\; \bigl(\,(1+xy)^3 z + y^2(1+xy)(4+3xy),\;\;
y + 3x(1+xy)^2 z + 3xy^2(4+3xy),\;\;
2x - 3x^2y - x^3 z\,\bigr).
$$

Direct computation (`scripts/verify_counterexample.py`) gives
$\det DF \equiv -2$, so $F$ is a Keller map (constant nonzero Jacobian). Yet

$$
F\bigl(0,0,-\tfrac14\bigr)
= F\bigl(1,-\tfrac32,\tfrac{13}{2}\bigr)
= F\bigl(-1,\tfrac32,\tfrac{13}{2}\bigr)
= \bigl(-\tfrac14,0,0\bigr),
$$

so $F$ is not injective and a fortiori has no polynomial (or any) inverse.
Adjoining identity coordinates gives counterexamples in every dimension
$n\ge 3$; the case $n=2$ remains open. The example was announced on X by
[@\_\_alpoge\_\_](https://x.com/__alpoge__) on July 19, 2026, and a Lean 4
formalization exists ([deancureton/jacobian](https://github.com/deancureton/jacobian),
where a rescaled variant with Jacobian determinant $1$ is also verified). Note
that all three colliding points are real: the restriction $F:\mathbb{R}^3\to\mathbb{R}^3$
is a real polynomial local diffeomorphism with *constant* Jacobian that fails
to be injective — a stronger failure mode than Pinchuk's 1994 counterexample
to the strong real Jacobian conjecture in the plane, where the Jacobian is
nonvanishing but not constant [Pin94].

### 1.2 The zero-dimensional QFT dictionary

The QFT reading of the formal-inverse approach to the Jacobian conjecture goes
back, in precise form, to Abdesselam [Abd03]; the underlying rooted-tree
expansion of the formal inverse is due to Bass–Connell–Wright [BCW82] and
Wright [Wri87, Wri89]. For the present map the dictionary is:

- **Fields.** $\phi = (x,y,z)$, a 3-component scalar in $0$ spacetime
  dimensions (no position dependence; integrals over spacetime are absent).
- **Field equations.** $F(\phi) = J$, with external source
  $J = (a,b,c)$. One must be careful with the word "action" here: $DF$ is
  *not* symmetric (e.g. $\partial F_1/\partial z = (1+xy)^3$ while
  $\partial F_3/\partial x = 2 - 6xy - 3x^2z$), so **there is no potential
  $S(\phi)$ with $\nabla S = F$**. The equations of motion exist; a
  single-scalar-potential Lagrangian for them does not. See §5.3.
- **Propagator.** The linearization at the origin,
  $DF(0) = L$ with $L(x,y,z) = (z,\, y,\, 2x)$, is the inverse propagator;
  $L^{-1}(a,b,c) = (c/2,\, b,\, a)$.
- **Vertices.** The nonlinear (degree $\ge 2$) terms of $F$ are the
  interaction vertices.
- **Perturbation theory.** Solving $F(\phi)=J$ by Picard iteration
  $\phi \mapsto \phi + L^{-1}\bigl(J - F(\phi)\bigr)$
  generates exactly the sum over rooted tree Feynman graphs: each
  substitution of $\phi$ into a vertex grafts subtrees, each internal line
  carries $L^{-1}$, each leaf carries a source insertion $J$.
- **Loops.** Constant Jacobian $\Leftrightarrow$ the one-loop determinant
  $\log\det DF(\phi)$ is field-independent $\Leftrightarrow$ all loop
  corrections are constants. The Jacobian conjecture was thus the statement:
  *if the loop corrections of a polynomial theory are trivial, the tree-level
  theory is globally, polynomially solvable.* That statement is now false
  for $n \ge 3$.

`scripts/tree_expansion.py` implements this iteration in a truncated polynomial ring
and confirms $F(G(J)) = J$ to total order $10$ in $(a,b,c)$ and to order
$t^{60}$ along rays $J = t\,v$.

---

## 2. What we proved in the toy model: the physics reading

### 2.1 The tree series converges — to an algebraic function

The formal inverse $G(J)$ never terminates: every order $1\le d\le 10$ of the
multivariate series carries nonzero coefficients (`scripts/tree_expansion.py`, check
2), consistent with the nonexistence of a polynomial inverse. But it is very
far from being a "merely formal" series. Gröbner elimination
(`scripts/branch_locus.py`) shows that the $x$-component of any preimage of
$(a,b,c)$ satisfies the cubic

$$
p\,x^3 + q\,x + r = 0,\qquad
p = 27a^2c^2 - 18abc + b^3c - b^2 + 16a,\quad
q = 4 - 3bc,\quad
r = -2c,
$$

and the tree series satisfies this relation identically to all computed
orders (`scripts/tree_expansion.py`, check 3). **Perturbation theory is computing the
Taylor series of an explicit degree-3 algebraic function.** The corresponding
minimal cubics for $y$ and $z$ over $\mathbb{C}[a,b,c]$ are *monic* — this
will matter below.

### 2.2 The radius of convergence is set by escape to infinity

Along the generic ray $J = t\,(1,2,3)$, the coefficients of the series in $t$
(computed exactly to order $t^{60}$) give a Domb–Sykes extrapolated radius of
convergence $\approx 0.3018$. The discriminant of the cubic factors as

$$
\mathrm{disc}_x = -\,p\,\bigl(4q^3 + 27\,p\,r^2\bigr),
$$

and along this ray the nearest zero of $p(t)$ sits at $|t| = 0.302028$, while
the nearest zero of the second factor sits at $|t| = 0.261$. The measured
radius matches the zero of $p$, *not* the nearer zero of the other factor.
Newton path-tracking of the perturbative branch confirms this directly: the
branch passes through $t = -0.261$ with finite field values but **escapes to
infinity** as $t \to -0.302$. The convergence of the tree expansion is
limited by the escape-to-infinity locus, not by collisions of solution
sheets.

### 2.3 Exact geometry: non-properness, not ramification

Because $\det DF \equiv -2$, the map $F$ is everywhere a local biholomorphism
(étale): distinct solution sheets can never merge at finite points of
$\mathbb{C}^3$. The generic fiber has exactly 3 points. The two discriminant
components have entirely different characters:

- $\{4q^3 + 27pr^2 = 0,\ p \neq 0\}$: two *distinct* fiber points share the
  same $x$-coordinate. The projection to the $x$-axis ramifies; the covering
  itself does not. Harmless.
- $\{p = 0\}$: the cubic drops degree, i.e. **a sheet escapes to infinity**.
  Since the $y$- and $z$-eliminants are monic, escape happens only in the
  $x$-direction. This is exactly the *non-properness set* $S_F$ of $F$ in the
  sense of Jelonek [Jel93]: the set of targets $J$ admitting no neighborhood
  with compact preimage closure. Jelonek proved that for a dominant polynomial
  map $\mathbb{C}^n\to\mathbb{C}^n$ this set is either empty or a uniruled
  hypersurface, computable effectively; here it is the explicit quartic
  hypersurface $\{p = 0\}$. Off $S_F$, the map $F$ restricts to a genuine
  3-sheeted covering of $\mathbb{C}^3\setminus\{p=0\}$.

The right slogan: *for a Keller map the only possible pathology is
non-properness*, and the counterexample realizes it.

### 2.4 Vacua at infinity, made exact

The origin satisfies $p(0,0,0) = 0$: **the perturbative vacuum $J=0$ lies on
the non-properness set.** The fiber over $J=0$ contains a single finite point
$\phi = 0$; the other two sheets of the generic 3-sheeted structure already
sit at infinity over the vacuum. They can be exhibited exactly. On the source
ray $(a,0,0)$ the cubic degenerates to $16a\,X^3 + 4X = 0$, and the fiber
consists of the perturbative point $\phi_A = (0,0,a)$ together with the two
points of the curve

$$
\phi(s) = \Bigl(\pm s,\ \mp\frac{3}{2s},\ \frac{13}{2s^2}\Bigr),
\qquad
F(\phi(s)) = \Bigl(-\frac{1}{4s^2},\,0,\,0\Bigr),
$$

on which $F$ is 2-to-1. As $a\to 0^-$ these sheets recede like
$|x| = 1/(2\sqrt{-a}) \to \infty$; at $a = -\tfrac14$ (i.e. $s = \pm 1$) they
land at $(\pm 1, \mp\tfrac32, \tfrac{13}{2})$ and collide *in the target*
with the image of the perturbative branch — the famous triple point. These
two "vacua at infinity" carry no signature in the coefficients of the tree
expansion at any order: perturbation theory around $J=0$ converges, sums to a
perfectly good local branch, and is structurally blind to two thirds of the
solution set.

**Summary of the physics reading.** In this model, perturbation theory does
not fail by divergence. It fails by *locality in field space*: the sum over
tree graphs reconstructs one sheet of a three-sheeted algebraic covering, and
the remaining sheets are separated from it not by an energy barrier but by
the boundary at infinity of field space, entering the finite region only
across the non-properness hypersurface $\{p = 0\}$.

---

## 3. Corrections to a circulating commentary

A widely-shared AI-generated commentary (produced by Gemini) on the
counterexample and its QFT interpretation contains three errors worth
correcting explicitly, since they invert the actual lessons of the example.

1. **"The tree series has zero radius of convergence / is asymptotic because
   the number of trees grows like $N!$."** False. The number of rooted trees
   relevant here grows only exponentially once symmetry factors are included
   (the $N!$ of labeled trees is cancelled by the $1/N!$ Taylor factors —
   this is the standard Cayley/Catalan bookkeeping, cf. [BCW82, Wri89]), and
   in any case the question is settled empirically in this repository: the
   series along $J = t(1,2,3)$ has radius $\approx 0.302 \neq 0$, matching
   the exact algebraic singularity. The tree expansion is a convergent
   series. Factorial growth and Borel-summable divergence are phenomena of
   *loop* expansions; this model has no field-dependent loops at all.

2. **"Borel resummation is impossible because of transcendental cuts."**
   Misleading twice over. No Borel summation is needed for a convergent
   series; and the non-perturbative completion here is not some intractable
   transcendental object but an explicit degree-3 *algebraic* function, the
   root of $p\,x^3+qx+r$ with polynomial coefficients — about the most
   tractable non-perturbative completion imaginable. Its singularities are
   algebraic branch points along explicitly computable hypersurfaces.

3. **"This explains why $D=4$ QFT is harder than $D=2,3$."** A category
   error. What separates $D=4$ from $D\le 3$ in constructive QFT is
   ultraviolet behavior — super-renormalizability of $\phi^4_2,\phi^4_3$
   versus strict renormalizability (and suspected triviality) in $D=4$
   [GJ87, EMS75, MS77]. The phenomenon exhibited by the counterexample is a
   statement about *global field-space geometry at tree level*, in zero
   dimensions, with no UV divergences anywhere in sight. It is logically
   independent of renormalization. See §4.2.

---

## 4. Implications for rigorous QFT in $D=4$: an honest assessment

The motivation for these notes is the question of what, if anything, the
counterexample means for the program of constructing interacting QFTs
rigorously — constructive QFT in the Glimm–Jaffe tradition [GJ87] and
perturbative algebraic QFT (pAQFT) in the Fredenhagen–Rejzner tradition
[FR12, FR16, Rej16]. We separate what can be said from what cannot.

### 4.1 What can be said

**(i) A sharp lesson on formal-power-series constructions.** The example is a
rigorous, minimal demonstration that local/perturbative data — even with
*trivial loop structure* and even when the perturbation series *converges* —
cannot control the global solution manifold of nonlinear field equations. In
pAQFT, observables of interacting theories are constructed as formal power
series in $\hbar$ and the coupling [FR12, FR16, Rej16]; convergence and
global field-space questions are deliberately deferred, and this is a
feature, not a bug, of the framework. The counterexample calibrates what that
deferral costs in the worst case: there exist polynomial field equations for
which the perturbative branch is analytic and complete in itself, yet the
exact theory has additional solution sectors at infinite field strength that
no reorganization of perturbation theory around the given vacuum can detect.
Any eventual "convergence + globality" upgrade of formal constructions has to
supply an input beyond the perturbative data — here, properness of the
classical field map, which is exactly what fails.

**(ii) Field redefinitions must be checked for global invertibility.**
Polynomial and formal field redefinitions $\phi \mapsto \phi' = F(\phi)$ with
$F(\phi) = L\phi + O(\phi^2)$, $\det L \neq 0$, are used throughout
renormalization theory, and the equivalence theorem for the $S$-matrix rests
on their invertibility, usually established at the level of formal power
series — where it is automatic. The counterexample shows that for $\ge 3$
field components a *polynomial* redefinition with *constant unit Jacobian*
(rescale one component of $F$ by $-\tfrac12$) need not be globally injective:
the formal inverse exists, converges near $\phi = 0$, and is still only a
local section of a 3-sheeted covering. Purely perturbative equivalence
statements are untouched. But any argument that treats such a redefinition as
an exact change of variables in a functional integral or a lattice measure —
"$\det DF = 1$, so the measure is preserved" — is using a false step: the
pushforward under a non-injective $F$ overcounts sheets, and the domain of
integration is not mapped bijectively. Global injectivity is an honest
hypothesis that must be verified, not read off from the Jacobian.

**(iii) A non-perturbative mechanism distinct from instantons.** The standard
catalogue of non-perturbative effects centers on finite-action saddle points
(instantons): configurations at finite distance in field space, weighted by
$e^{-S_{\mathrm{inst}}/\hbar}$, connected to the vacuum by finite-action
interpolation. The mechanism here is different in kind. The extra sectors are
not saddle points of the same action at finite field values; they live on the
boundary at infinity of field space and communicate with the finite region
only across the non-properness hypersurface of the field map — transitions
"through infinity." In the toy model this is not a metaphor: the curve
$\phi(s)$ of §2.4 realizes it exactly. Whether an analogous
boundary-of-field-space mechanism can survive in a genuine functional-integral
setting (where infinite field values are suppressed by the measure) is open,
but the example establishes that at the level of classical solution manifolds
the mechanism exists and is compatible with completely trivial loop
structure.

### 4.2 What cannot be said

Plainly:

- **Nothing about UV renormalization.** The model is zero-dimensional; there
  are no short-distance singularities, no counterterms, no running. The
  phenomenon is invariant under everything renormalization theory cares
  about.
- **Nothing about Borel summability of $\phi^4_2$ or $\phi^4_3$.** The
  divergence-with-Borel-summability of those perturbation series is a
  loop/large-order phenomenon, proved by Eckmann–Magnen–Sénéor for $P(\phi)_2$
  [EMS75] and Magnen–Sénéor for $\phi^4_3$ [MS77]. Our series is convergent;
  the two situations are disjoint, and the counterexample neither threatens
  nor illuminates those theorems.
- **Nothing about the core $D=4$ difficulties.** Suspected triviality of
  $\phi^4_4$ and the Yang–Mills existence/mass-gap problem are ultraviolet
  and measure-theoretic problems. The counterexample does not bear on them,
  and claims to the contrary (§3, item 3) should be resisted.

The honest headline is narrower and still interesting: *constant-Jacobian
polynomial dynamics can hide non-perturbative sectors at infinity, and no
amount of perturbative or even locally-analytic information detects them.*

### 4.3 Algebraic QFT specifically (Haag–Kastler; Fredenhagen–Rejzner; Buchholz–Fredenhagen)

Since the question "what does this mean for AQFT?" is naturally raised, we
spell out the contact points layer by layer, in decreasing order of rigor of
what we can actually assert.

**(a) The algebraic viewpoint is vindicated by an exact 0D statement.**
Haag–Kastler AQFT [HK64] takes the *algebras of observables*, not field
coordinates, as the primary data; fields are regarded as interchangeable
"coordinates" on the theory (Borchers classes). The counterexample turns this
philosophical preference into a theorem-sized fact in the only setting we
control completely. The pullback

$$
F^*:\ \mathbb{C}[a,b,c] \longrightarrow \mathbb{C}[x,y,z],
\qquad
F^*(\mathcal{O}) = \mathcal{O}\circ F,
$$

is an injective algebra homomorphism that is **not surjective**: the
polynomial observables of the "redefined field" $F(\phi)$ form a *proper*
subalgebra of the observables of $\phi$, of index measured by the degree-3
field extension $\mathbb{C}(x,y,z)\,/\,\mathbb{C}(a,b,c)$ with Galois closure
group $S_3$ (this is exactly the verified eliminant cubic and monodromy;
`scripts/branch_locus.py`, `scripts/monodromy.py`). The observables missing
from the image are precisely those that separate the three sheets — and this
is now fully explicit: $\mathbb{C}(x,y,z)$ is free with basis $\{1, x, x^2\}$
over the fraction field of $\operatorname{im}F^*$, every observable has a
unique normal form $c_0(F) + c_1(F)x + c_2(F)x^2$ (membership iff
$c_1 = c_2 = 0$), and no *finite-module* statement is possible because $x$ is
not integral over $\operatorname{im}F^*$ — an exact escape-curve certificate,
equivalent to non-properness (`scripts/missing_observables.py`,
`docs/MISSING_OBSERVABLES.md`). So: *a
polynomial field redefinition with invertible propagator and constant unit
Jacobian can be a monomorphism, but not an automorphism, of the observable
algebra.* Any framework whose objects are algebras and whose morphisms must
be checked for surjectivity is structurally protected against the mistake
this example punishes; any framework that treats "invertible field
redefinition" as certified by the Jacobian is not. This is, we believe, the
cleanest AQFT-flavored lesson of the counterexample, and it is exact.

**(b) pAQFT (Fredenhagen–Rejzner): the deferral is calibrated, not
contradicted.** In pAQFT [FR12, FR16, Rej16] interacting observables are
constructed from free ones by Møller-type maps built as formal power series
in the coupling and $\hbar$; classical inversions of nonlinear field
equations enter through the (classical) Møller map, and invertibility holds
automatically at the formal-series level. Nothing in the counterexample
exhibits an inconsistency in this: the constructions are internally coherent
as formal deformation quantization, and their formal character is a stated
feature of the framework. What the counterexample adds is a *calibration of
the worst case* of that deferral, sharper than previously available: the
formal inverse here is not merely consistent but **convergent**, analytic on
a neighborhood of the vacuum, and satisfies every identity perturbation
theory can formulate — and still describes only one of three solution
sectors, the others being invisible at every order because they sit at
infinite field strength. So "upgrade formal to convergent" is *not* the
missing step between pAQFT and a non-perturbative construction; the missing
input is global (properness of the classical dynamics), and it must come
from outside perturbation theory. Relatedly, statements of *perturbative
agreement* (independence of the split of the action into free and
interacting parts, [HW05, DHP17]) are field-redefinition-adjacent moves
established as formal-series identities; the example is a reminder of
exactly which global questions such identities do not decide.

**(c) Non-perturbative algebraic constructions: a minimal test case, not a
threat.** The Buchholz–Fredenhagen $C^*$-algebraic approach [BF20] builds
interacting dynamics from unitaries $S(f)$ subject to causal factorization
relations, bypassing formal series altogether. Nothing in a 0D polynomial
map bears on the correctness of that program. The counterexample instead
supplies the minimal instance of a phenomenon any non-perturbative framework
must *represent* somewhere in its data: classical dynamics whose solution
count jumps (here $1 \leftrightarrow 3$) as an external source crosses an
algebraic wall, with the extra solutions entering from infinity in field
space and permuted by an $S_3$ monodromy under cycles of the source. A
well-posed exercise — open, and we claim nothing about its outcome — is to
formulate the 0D caricature of the $S(f)$ relations for this $F$ and see how
the wall $\{p=0\}$ and the sheet structure are encoded.

**(d) What the 0D model cannot say about AQFT.** There is no spacetime, no
net of local algebras, no causality, no Hilbert space, no states, no
superselection theory here. Analogies between "extra solution sheets" and
"sectors" are structural motivation, not theorems. In particular the
existence results of constructive AQFT in $D=2,3$ and the open problems in
$D=4$ are entirely untouched, per §4.2.

### 4.4 Claims ledger

Because results of this kind are easily inflated in transmission (§3 is the
proof), we keep a ledger: each claim, where it is verified, and the nearest
overstatement that it does **not** license.

| Verified claim | Where | Does **not** imply |
|---|---|---|
| $\det DF \equiv -2$, yet $F$ is not injective (3 preimages of $(-\tfrac14,0,0)$) | `scripts/verify_counterexample.py` | anything about existence/triviality of interacting QFT in $D \ge 1$ |
| The tree expansion converges with finite nonzero radius (≈ 0.302 along $t(1,2,3)$), set by the exact branch locus | `scripts/tree_expansion.py` | "zero radius of convergence"; any Borel-summability statement; any contact with large-order/loop divergences |
| Escape locus is $\{p=0\}$; only the $x$-sheet escapes; the vacuum $J=0$ lies on the wall | `scripts/branch_locus.py` | a tunneling/instanton interpretation — there is no action here ($DF$ is not symmetric, the dynamics is not variational) |
| Geometric monodromy of the 3 sheets is the full $S_3$ | `scripts/monodromy.py`, `docs/MONODROMY.md` | an anomaly in the cohomological sense; no symmetry is broken and nothing is loop-generated |
| Real change-of-variables defect: $N(J) = 3$ iff $p<0$, measured $A(\sigma)\to 2$ as $\sigma\to 0$ | `scripts/measure_anomaly.py` | failure of the perturbative equivalence theorem — formal-series equivalence statements are untouched |
| Sheet-summed observables are rational with poles only on $\{p=0\}$; residues factorize at the wall | `scripts/trace_pushforward.py` | existence of a positive geometry / canonical form for the chambers — for *this* map that question is now settled negatively (`docs/POSITIVE_GEOMETRY.md`); the classification across Keller maps remains open |
| $F^*$ is a monomorphism, not an automorphism, of polynomial observable algebras; extension degree 3, Galois group $S_3$ | `scripts/branch_locus.py` + `scripts/monodromy.py` | any statement about automorphisms of local nets or Borchers classes in $D\ge 1$ |
| Missing observables made explicit: basis $\{1,x,x^2\}$ over $\operatorname{im}F^*$, unique normal form $c_0(F)+c_1(F)x+c_2(F)x^2$, separator coefficients carry the factor $p$; $x$ not integral over $\operatorname{im}F^*$ (exact escape certificate) | `scripts/missing_observables.py`, `docs/MISSING_OBSERVABLES.md` | a finite-module (polynomial-level) structure theorem — non-properness forbids it, provably |
| Alpöge–Mathew is rigid modulo gauge to first order within the equivariant degree box | `scripts/search_counterexamples.py` | uniqueness of counterexamples — the check is first-order, box-limited, and within one equivariance class |
| The wall $\{p=0\}$ reduces, by C\*-invariance, to a cuspidal plane cubic; the cuspidal tangent is $\{D_0=0\}$ | `scripts/positive_geometry.py` | any statement about walls of other (hypothetical) counterexamples |
| $F(\mathbb{C}^3)$ misses exactly one C\*-orbit ($ac^2=\tfrac{4}{27}$, $bc=\tfrac43$): sources with empty fiber | `scripts/positive_geometry.py` | a physical "instability" — there is no dynamics or energy here, only an exact solvability statement |
| The $N=3$ chamber admits **no** canonical form: it is not a positive geometry (residueless double pole at the cusp; adjoint-line system has only the zero solution) | `scripts/positive_geometry.py`, `docs/POSITIVE_GEOMETRY.md` | irrelevance of positive geometry — nodal-wall chambers would pass; the failure mode is itself structured |

---

## 5. Open questions and a small research program

### 5.1 Monodromy of the three sheets

Off $\{p=0\}$, $F$ is a 3-sheeted covering of
$\mathbb{C}^3\setminus\{p=0\}$ with connected total space, so its monodromy
group is a transitive subgroup of $S_3$ — either $\mathbb{Z}/3$ or $S_3$.
Which one? Loops around which strata of $\{p=0\}$ exchange the perturbative
sheet with the sheets at infinity? A concrete computation (track the three
roots of the cubic around loops in the $(a,b,c)$-space) would make the
"tunneling through infinity" picture of §4.1(iii) quantitative and would
determine the full analytic continuation of the tree expansion.

### 5.2 Classification of counterexamples

Is the Alpöge–Mathew map isolated, or the first point of a family? Natural
invariants suggested by the QFT reading: the generic fiber degree ($3$ here —
the minimum possible for a non-injective Keller map?), the geometry of the
non-properness hypersurface $S_F$ (degree, singularities, position relative
to $F(\mathbb{C}^n)$), and whether the vacuum $J=0$ can be moved off $S_F$ by
composition with affine maps. A classification for low degree and $n = 3$
now looks like a finite, attackable problem.

### 5.3 The 0D functional integral — with the gradient subtlety stated

One would like to ask: does the partition function
"$Z(J) = \int d^3\phi\, e^{-S(\phi) + J\cdot\phi}$" see all three sheets?
**The question as stated is ill-posed for this map**: as noted in §1.2, $DF$
is not symmetric, so there is *no* potential $S$ with $\nabla S = F$; the
field equations are not variational in the naive sense. The well-posed
formulations use auxiliary fields, e.g. the conjugate-field (source)
representation that Abdesselam introduced for exactly this class of problems
[Abd03], where the formal inverse appears as a correlation function of a
model with action $\bar\phi\cdot(F(\phi) - J)$ plus a (here trivial, since
$\det DF$ is constant) fermionic representation of the Jacobian determinant;
the associated tree/forest combinatorics is the constructive-QFT toolkit of
Abdesselam–Rivasseau [AR95]. Two sharp versions of the question:

- The distributional integral
  $\int_{\mathbb{R}^3} d^3\phi\;\delta^3\bigl(F(\phi)-J\bigr)$ equals
  $N(J)/2$, where $N(J)$ is the number of finite real preimages. It *jumps*
  across the real points of $\{p=0\}$ — the crudest possible "observable"
  already sees the sheet loss. What smooth damped versions, e.g.
  $Z_\hbar(J) = \int d^3\phi\, e^{-|F(\phi)-J|^2/2\hbar}$, converge, and do
  their $\hbar\to 0$ asymptotics resolve the individual sheets, including
  the ones at infinity?
  *(✅ Resolved 2026-07-21, `docs/DAMPED_PARTITION.md`: $Z_\hbar$ is finite
  for ALL $J$ — exactly $(2\pi\hbar)^{3/2}(\tfrac12 +
  \mathbb{P}[p(J+\sqrt{\hbar}\xi)<0])$ — the Jelonek set is detected not by
  divergence but by a piecewise-constant semiclassical prefactor $N(J)/2$
  jumping across the wall, with uniformity boundary
  $\hbar \ll \mathrm{dist}(J,\text{wall})^2$. Sheet resolution beyond total
  mass remains open.)*
- Is there a *gradient* counterexample, i.e. a Keller map $F = \nabla S$
  with symmetric $DF$ that is non-injective? ✅ **Resolved 2026-07-25**
  (`scripts/symmetric_search.py`; `docs/SYMMETRIC_SEARCH.md`): **yes,
  explicitly, in dimension 6** — the cotangent lift
  $W_6 = \bar\varphi\cdot F(\varphi)$ (the first-order action of §1.2
  itself!) has $\det\operatorname{Hess}W_6 \equiv -4$ and $\nabla W_6$ is
  3:1 with rational witnesses; a normalized dBvdE twisted lift over
  $\mathbb{Q}(i)$ is also written down. Complement: in dimension 3 the AM
  map is not a gradient in ANY affine frame ($K\,DF$ symmetric $\implies
  K = 0$), and the $(1,-1,-m)$ gradient family contains only tame shears.
  The measure-theoretic payoff is *conditional*, though: coercivity is
  provably unreachable ($\kappa \le 0$ is forced, and $W_6$ is affine in
  $\bar\varphi$), so the honest action exists but $\int e^{-W}$ is never
  absolutely convergent. See `docs/SYMMETRIC_SEARCH.md` §7.

### 5.4 The $n=2$ case

The two-variable Jacobian conjecture is untouched. In the QFT language: does
a 2-component 0D scalar model with constant Jacobian always have a proper
classical field map? Equivalently, can the non-properness set of a Keller map
in dimension 2 be nonempty? The structure theory of $S_F$ in dimension 2
(Jelonek [Jel93] and subsequent work) makes this a focused question about
curves of non-properness, and any insight would now carry unusual weight in
both directions.

### 5.5 Bass–Connell–Wright reduction in QFT terms

The degree-reduction theorem [BCW82] says it suffices to consider maps
$F = \mathrm{id} + H$ with $H$ cubic homogeneous and $DH$ nilpotent, at the
price of increasing $n$; de Bondt–van den Essen sharpen this to $DH$
nilpotent *and symmetric* [dBvdE05]. In QFT terms: every Keller theory is
stably equivalent to a theory with only cubic vertices and nilpotent
(one-loop-exact, in fact zero) mass mixing — a "$\phi^3$-type" normal form —
and even to a gradient such theory. The counterexample therefore guarantees,
by contraposition, explicit non-invertible cubic-homogeneous nilpotent models
in some higher dimension. Exhibiting one, and understanding what its
nilpotent vertex structure means physically (a BRST-like triangular
structure that nevertheless fails to integrate to a global field
redefinition), seems to us the most concrete next step in the QFT direction.

---

## 6. Addendum: anomaly-adjacent structures and a measured field-redefinition defect

*(Added after the monodromy computation; quantitative results from
`scripts/measure_anomaly.py`.)*

### 6.1 A field-redefinition "measure anomaly", measured

$F$ restricts to a real map $\mathbb{R}^3\to\mathbb{R}^3$ with constant
Jacobian $-2$ — a local diffeomorphism everywhere, and the three preimages of
$(-\tfrac14,0,0)$ are all real. For a non-injective local diffeomorphism the
change-of-variables formula carries a multiplicity:

$$
\int f(F(\phi))\,\lvert\det DF\rvert\, d^3\phi
\;=\; \int f(J)\, N(J)\, d^3J,
\qquad N(J) = \#\{\text{real preimages of } J\}.
$$

For this map the multiplicity has a closed form (**chamber rule**, exact):
the monic discriminant of the $x$-eliminant is $-4D_0^2/p^3$ with
$D_0 = 27ac^2-9bc+8$, hence

$$
N(J) \;=\; \begin{cases} 3, & p(a,b,c) < 0,\\ 1, & p(a,b,c) > 0.\end{cases}
$$

The naive equivalence-theorem manipulation ("the Jacobian is a harmless
constant, so substitute $\phi' = F(\phi)$") is therefore wrong by the factor
$N(J)$, a *step function* in source space. Monte Carlo over Gaussian source
ensembles of width $\sigma$ gives the measured anomaly factor
$A(\sigma) = \langle N\rangle = 1 + 2\,P[p<0]$:

| $\sigma$ | $10$ | $1$ | $0.1$ | $0.01$ | $0.001$ |
|---|---|---|---|---|---|
| $A(\sigma)$ | $1.205$ | $1.696$ | $2.006$ | $2.001$ | $2.000$ |

(statistical errors $\pm 3\cdot 10^{-4}$, $10^7$ samples per width). The
crucial feature: **$A \to 2$ as $\sigma \to 0$.** Because the vacuum $J=0$
lies *on* the wall $\{p=0\}$ (the linear part of $p$ is $16a$), the defect is
$O(1)$ for sources concentrated arbitrarily close to the perturbative vacuum;
it fades only for wide ensembles dominated by the positive quartic term
$27a^2c^2$ in $p$. No order of perturbation theory around $J = 0$ detects any
of this.

Three standard practices sit directly upstream of this failure mode:

1. **Nonlinear field redefinitions / equivalence theorems.** Invertibility is
   in practice certified through the Jacobian. The counterexample shows that
   for $\ge 3$ field components a polynomial redefinition with *constant*
   Jacobian can be non-invertible, and §6.1 quantifies the resulting measure
   defect. Global injectivity is an independent hypothesis that must be
   checked by global means (e.g. properness/degree arguments), not by the
   Jacobian.
2. **Gribov copies without a Gribov horizon.** Gauge fixing assumes a slice
   map whose Faddeev–Popov determinant certifies local uniqueness; Gribov's
   copies [Gri78, Sin78] are traditionally associated with the *horizon*,
   where that determinant vanishes. Here is a map whose FP-type determinant
   is the constant $-2$ — there is no horizon anywhere in finite field space
   — and which nevertheless has copies: the horizon has retreated to infinity.
   Absence of a horizon does not exclude copies.
3. **Nicolai maps.** In supersymmetric models the Nicolai map [Nic80] is a
   nonlinear field transformation whose Jacobian determinant cancels the
   fermion determinant; its global invertibility is generally taken for
   granted. A non-injective Nicolai-type map with the properties above would
   make the "free-field representation" silently miss solution sectors.

### 6.2 Is this a new class of anomalies?

Not in the strict (cohomological) sense: no symmetry is broken, nothing is
generated by loops or by regularization, and there is no Wess–Zumino-type
consistency structure. What it shares with anomalies is the defining pattern
— *a manipulation valid by every local criterion fails through the global
structure of configuration space/measure* — and its closest relatives in the
modern literature are global rather than perturbative anomalies:

- **Anomalies in parameter space** [CFLS20]: dialing external couplings
  around a closed cycle implements a nontrivial transformation. Our verified
  $S_3$ monodromy (see `docs/MONODROMY.md`) is exactly this pattern at tree
  level: a closed loop of the external source returns the same theory with
  the solution branches permuted — vacuum monodromy as in Seiberg–Witten
  theory [SW94], but mediated entirely through infinity in field space, with
  no massless states and no degeneration at finite distance.
- We propose *non-properness defect* as the accurate name: a global
  consistency condition (properness of the classical field map, equivalently
  $N(J)$ locally constant $= \deg F$) that is strictly stronger than all
  local conditions (nonvanishing Jacobian) and that standard arguments
  implicitly assume.

### 6.3 Relation to the amplitudes program (positive geometry, CHY)

The findings resonate with — and give a sharp solvable caricature of — three
themes of the modern S-matrix program:

1. **Diagrams are a local triangulation of a global object.** The
   amplituhedron program [AT14, ABL17] replaced sums of Feynman diagrams by a
   single geometric object whose triangulations reproduce the diagrammatic
   expansion. Here the analogous statement is a theorem: the rooted-tree
   expansion is the Taylor triangulation, around one point, of one branch of
   the global covering $\{F(\phi)=J\}$; the sheets, the monodromy and the
   escape locus are properties of the geometry that no collection of diagrams
   sees.
2. **Singularities at infinity.** The radius of convergence of the tree
   expansion is set by $\{p=0\}$, where a solution *escapes to infinity in
   field space* — not by any finite-distance degeneration. This is the
   field-space analogue of *second-type Landau singularities*, which arise
   from pinches at infinity in loop-momentum space and are notoriously the
   ones not visible from the graph combinatorics. The chamber structure cut
   by the discriminant into $N=1$ and $N=3$ regions, with wall-crossing at
   $\{p=0\}$, is precisely the kind of stratified real geometry the positive
   geometry program axiomatizes (boundaries $\leftrightarrow$ singularities).
3. **Summing over all solutions restores simplicity (CHY analogy).** In the
   CHY formalism [CHY14] amplitudes are sums over *all* solutions of the
   scattering equations, and only that sum is a rational function. Identically
   here: the single perturbative branch is multivalued with $S_3$ monodromy,
   but the elementary symmetric functions over all three sheets are
   *rational* in the sources,
   $$e_1 = 0,\qquad e_2 = \frac{q}{p},\qquad e_3 = -\frac{r}{p},$$
   with poles exactly on the non-properness locus $\{p=0\}$. "Sum over all
   vacua" observables are single-valued; perturbation theory's pathology is
   an artifact of selecting one solution of the field equations.

Whether a genuine positive-geometry structure (a canonical form with
logarithmic singularities on the chamber walls) underlies the $N(J)$
stratification of Keller maps is an open — and now well-posed — question.

---

## 7. References

- [Abd03] A. Abdesselam, *The Jacobian conjecture as a problem of
  perturbative quantum field theory*, Ann. Henri Poincaré **4** (2003)
  199–215. arXiv:math/0208173. DOI:10.1007/s00023-003-0127-7.
- [AR95] A. Abdesselam, V. Rivasseau, *Trees, forests and jungles: a
  botanical garden for cluster expansions*, in *Constructive Physics*,
  Lecture Notes in Physics **446**, Springer (1995) 7–36.
  arXiv:hep-th/9409094.
- [BCW82] H. Bass, E. H. Connell, D. Wright, *The Jacobian conjecture:
  reduction of degree and formal expansion of the inverse*, Bull. Amer.
  Math. Soc. (N.S.) **7** (1982) 287–330.
  DOI:10.1090/S0273-0979-1982-15032-7.
- [BF20] D. Buchholz, K. Fredenhagen, *A C\*-algebraic approach to
  interacting quantum field theories*, Comm. Math. Phys. **377** (2020)
  947–969. arXiv:1902.06062.
- [dBvdE05] M. de Bondt, A. van den Essen, *A reduction of the Jacobian
  conjecture to the symmetric case*, Proc. Amer. Math. Soc. **133** (2005)
  2201–2205. DOI:10.1090/S0002-9939-05-07570-2.
- [DHP17] N. Drago, T.-P. Hack, N. Pinamonti, *The generalised principle of
  perturbative agreement and the thermal mass*, Ann. Henri Poincaré **18**
  (2017) 807–868. arXiv:1502.02705.
- [EMS75] J.-P. Eckmann, J. Magnen, R. Sénéor, *Decay properties and Borel
  summability for the Schwinger functions in $P(\phi)_2$ theories*, Comm.
  Math. Phys. **39** (1975) 251–271. DOI:10.1007/BF01705374.
- [FR12] K. Fredenhagen, K. Rejzner, *Batalin–Vilkovisky formalism in
  perturbative algebraic quantum field theory*, Comm. Math. Phys. **317**
  (2012) 697–725. arXiv:1110.5232.
- [FR16] K. Fredenhagen, K. Rejzner, *Perturbative algebraic quantum field
  theory*, in *Mathematical Aspects of Quantum Field Theories*, Springer
  (2016). arXiv:1208.1428.
- [GJ87] J. Glimm, A. Jaffe, *Quantum Physics: A Functional Integral Point
  of View*, 2nd ed., Springer, New York (1987).
- [HK64] R. Haag, D. Kastler, *An algebraic approach to quantum field
  theory*, J. Math. Phys. **5** (1964) 848–861. DOI:10.1063/1.1704187.
- [HW05] S. Hollands, R. M. Wald, *Conservation of the stress tensor in
  perturbative interacting quantum field theory in curved spacetimes*, Rev.
  Math. Phys. **17** (2005) 227–311. arXiv:gr-qc/0404074.
- [Jel93] Z. Jelonek, *The set of points at which a polynomial map is not
  proper*, Ann. Polon. Math. **58** (1993) 259–266.
  DOI:10.4064/ap-58-3-259-266.
- [Kel39] O.-H. Keller, *Ganze Cremona-Transformationen*, Monatsh. Math.
  Phys. **47** (1939) 299–306. DOI:10.1007/BF01695502.
- [MS77] J. Magnen, R. Sénéor, *Phase space cell expansion and Borel
  summability for the Euclidean $\phi^4_3$ theory*, Comm. Math. Phys.
  **56** (1977) 237–276. DOI:10.1007/BF01614211.
- [Pin94] S. Pinchuk, *A counterexample to the strong real Jacobian
  conjecture*, Math. Z. **217** (1994) 1–4. DOI:10.1007/BF02571929.
- [Rej16] K. Rejzner, *Perturbative Algebraic Quantum Field Theory: An
  Introduction for Mathematicians*, Mathematical Physics Studies, Springer
  (2016).
- [vdE00] A. van den Essen, *Polynomial Automorphisms and the Jacobian
  Conjecture*, Progress in Mathematics **190**, Birkhäuser, Basel (2000).
- [Wri87] D. Wright, *Formal inverse expansion and the Jacobian conjecture*,
  J. Pure Appl. Algebra **48** (1987) 199–219.
- [Wri89] D. Wright, *The tree formulas for reversion of power series*,
  J. Pure Appl. Algebra **57** (1989) 191–211.
- [AM26] L. Alpöge, announcement of the counterexample (with the question
  posed by A. Mathew and the map produced by Claude Fable 5), X post by
  @\_\_alpoge\_\_, July 19, 2026; see also the entry "Jacobian conjecture"
  on Wikipedia and the Lean 4 formalization at
  [deancureton/jacobian](https://github.com/deancureton/jacobian).

References for the addendum (§6):

- [ABL17] N. Arkani-Hamed, Y. Bai, T. Lam, *Positive geometries and
  canonical forms*, JHEP **11** (2017) 039. arXiv:1703.04541.
- [AT14] N. Arkani-Hamed, J. Trnka, *The Amplituhedron*, JHEP **10** (2014)
  030. arXiv:1312.2007.
- [CFLS20] C. Córdova, D. S. Freed, H. T. Lam, N. Seiberg, *Anomalies in the
  space of coupling constants and their dynamical applications I*, SciPost
  Phys. **8** (2020) 001. arXiv:1905.09315.
- [CHY14] F. Cachazo, S. He, E. Y. Yuan, *Scattering of massless particles
  in arbitrary dimensions*, Phys. Rev. Lett. **113** (2014) 171601.
  arXiv:1307.2199.
- [Gri78] V. N. Gribov, *Quantization of non-Abelian gauge theories*, Nucl.
  Phys. B **139** (1978) 1–19.
- [Nic80] H. Nicolai, *On a new characterization of scalar supersymmetric
  theories*, Phys. Lett. B **89** (1980) 341–346.
- [Sin78] I. M. Singer, *Some remarks on the Gribov ambiguity*, Comm. Math.
  Phys. **60** (1978) 7–12.
- [SW94] N. Seiberg, E. Witten, *Electric-magnetic duality, monopole
  condensation, and confinement in N=2 supersymmetric Yang-Mills theory*,
  Nucl. Phys. B **426** (1994) 19–52. arXiv:hep-th/9407087.
