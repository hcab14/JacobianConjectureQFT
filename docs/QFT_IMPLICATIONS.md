# The Alpöge–Mathew Counterexample as a Zero-Dimensional Field Theory: What It Does and Does Not Teach Us About Rigorous QFT

**Status.** Working notes accompanying the computations in this repository
(`verify_counterexample.py`, `tree_expansion.py`, `branch_locus.py`). All
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

Direct computation (`verify_counterexample.py`) gives
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

`tree_expansion.py` implements this iteration in a truncated polynomial ring
and confirms $F(G(J)) = J$ to total order $10$ in $(a,b,c)$ and to order
$t^{60}$ along rays $J = t\,v$.

---

## 2. What we proved in the toy model: the physics reading

### 2.1 The tree series converges — to an algebraic function

The formal inverse $G(J)$ never terminates: every order $1\le d\le 10$ of the
multivariate series carries nonzero coefficients (`tree_expansion.py`, check
2), consistent with the nonexistence of a polynomial inverse. But it is very
far from being a "merely formal" series. Gröbner elimination
(`branch_locus.py`) shows that the $x$-component of any preimage of
$(a,b,c)$ satisfies the cubic

$$
p\,x^3 + q\,x + r = 0,\qquad
p = 27a^2c^2 - 18abc + b^3c - b^2 + 16a,\quad
q = 4 - 3bc,\quad
r = -2c,
$$

and the tree series satisfies this relation identically to all computed
orders (`tree_expansion.py`, check 3). **Perturbation theory is computing the
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
- Is there a *gradient* counterexample, i.e. a Keller map $F = \nabla S$
  with symmetric $DF$ that is non-injective? By de Bondt–van den Essen
  [dBvdE05], the Jacobian conjecture is equivalent to its symmetric
  (Hessian) case, so counterexamples of gradient type must exist in some
  dimension — but the reduction changes $n$, and no explicit one is known.
  A gradient example would make the measure-theoretic questions direct,
  with an honest action $S$.

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

## 6. References

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
- [dBvdE05] M. de Bondt, A. van den Essen, *A reduction of the Jacobian
  conjecture to the symmetric case*, Proc. Amer. Math. Soc. **133** (2005)
  2201–2205. DOI:10.1090/S0002-9939-05-07570-2.
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
