# Problem statement

## Background

On July 19, 2026, Levent Alpöge (with Akhil Mathew, assisted by the LLM Claude
Fable 5) announced a counterexample to the **Jacobian conjecture** (Keller,
1939): the polynomial map $F : \mathbb{C}^3 \to \mathbb{C}^3$,

$$
F(x,y,z) = \bigl(\,(1+xy)^3 z + y^2(1+xy)(4+3xy),\;
                y + 3x(1+xy)^2 z + 3xy^2(4+3xy),\;
                2x - 3x^2y - x^3 z \,\bigr),
$$

has **constant Jacobian determinant $-2$** yet is **not injective**: the three
points $(0,0,-\tfrac14)$, $(1,-\tfrac32,\tfrac{13}{2})$,
$(-1,\tfrac32,\tfrac{13}{2})$ all map to $(-\tfrac14,0,0)$. Hence the
conjecture is false in dimension $n \ge 3$ (the case $n = 2$ remains open).
The example was verified by the community and formalized in Lean 4
(`deancureton/jacobian`). All claims above are independently re-verified in
this repository (`scripts/verify_counterexample.py`).

## Why this matters for QFT

The Jacobian conjecture has a well-known reformulation in zero-dimensional
perturbative quantum field theory (Abdesselam; Abdesselam–Rivasseau; building
on the tree formulas of Bass–Connell–Wright and Wright): the equations
$F(\varphi) = J$ are the classical field equations of a 3-component scalar
model in $D = 0$ with external sources $J$, the linearization $DF(0)$ is the
inverse propagator, the nonlinear monomials of $F$ are interaction vertices,
and the perturbative inversion $\varphi(J)$ is exactly the **sum over rooted
tree Feynman graphs**. A constant Jacobian determinant means the loop
corrections are field-independent constants — the theory looks "free" at loop
level.

A genuine counterexample therefore provides, for the first time, an *exactly
solvable* toy model in which perturbation theory is well defined at every
order, loop corrections are trivial, and yet the theory fails globally.

## The model, stated as a Lagrangian

*(Reference labels as in `docs/QFT_IMPLICATIONS.md` §7. All expansions below
are verified symbolically against `jcqft/core.py`.)*

### No single-field action exists

One would like to write $F(\varphi) = J$ as the stationarity condition of an
action $S(\varphi)$. That is impossible: $F$ is the gradient of a potential
iff $DF$ is symmetric, and here it is not — already at the origin,

$$
\frac{\partial F_1}{\partial z}\Big|_{\varphi=0} = 1
\;\neq\;
2 = \frac{\partial F_3}{\partial x}\Big|_{\varphi=0},
$$

since $\partial F_1/\partial z = (1+xy)^3$ while
$\partial F_3/\partial x = 2 - 6xy - 3x^2 z$. Even the linearization
$L = DF(0)$ is non-symmetric, so no correction by nonlinear terms can help;
the antisymmetric part at the origin is
$\tfrac12(L - L^{\mathsf T}) = \tfrac12\bigl(\begin{smallmatrix}0&0&-1\\0&0&0\\1&0&0\end{smallmatrix}\bigr)$.
So the naive "Lagrangian of the counterexample" is ill-posed; see
`docs/QFT_IMPLICATIONS.md` §5.3 for what this blocks (honest partition
functions, Lefschetz thimbles) and for the gradient-counterexample question.

### The correct (first-order) formulation

The standard remedy, used by Abdesselam for exactly this class of problems
[Abd03], is a conjugate (auxiliary) field
$\bar\varphi = (\bar x, \bar y, \bar z)$ and the first-order action

$$
S(\bar\varphi, \varphi) \;=\; \bar\varphi \cdot \bigl(F(\varphi) - J\bigr),
\qquad
Z(J) \;=\; \int d^3\bar\varphi\, d^3\varphi\; e^{-S(\bar\varphi,\varphi)},
$$

where the $\bar\varphi$-integral (along the imaginary axis, in the usual
formal sense) localizes $Z(J)$ on $\delta^3\bigl(F(\varphi) - J\bigr)$. In
general one must also include a fermionic pair $(\bar\psi, \psi)$ with action
$\bar\psi \cdot DF(\varphi)\, \psi$ to produce the Jacobian factor
$\det DF(\varphi)$; here that ghost sector is trivial, since
$\det DF \equiv -2$ means the fermions contribute only the constant $-2$.
The tree/forest combinatorics needed to control such expansions
constructively is that of Abdesselam–Rivasseau [AR95].

### Free part, propagator, and interaction vertices

Split $S = S_{\mathrm{free}} + S_{\mathrm{int}}$ with

$$
S_{\mathrm{free}} = \bar\varphi \cdot L\,\varphi - \bar\varphi\cdot J,
\qquad
L = \begin{pmatrix} 0 & 0 & 1\\ 0 & 1 & 0\\ 2 & 0 & 0\end{pmatrix},
\qquad
L^{-1} = \begin{pmatrix} 0 & 0 & \tfrac12\\ 0 & 1 & 0\\ 1 & 0 & 0\end{pmatrix},
$$

so the only nonzero propagators are
$\langle \varphi_i\, \bar\varphi_j \rangle = (L^{-1})_{ij}$, i.e.
$\langle x\bar z\rangle = \tfrac12$, $\langle y\bar y\rangle = 1$,
$\langle z\bar x\rangle = 1$. The interaction is
$S_{\mathrm{int}} = \bar\varphi\cdot V(\varphi)$ with
$V(\varphi) = F(\varphi) - L\varphi$; expanding into monomials, the complete
list of vertices is

$$
S_{\mathrm{int}}
= \bar x\,\bigl(4y^2 + 3xyz + 7xy^3 + 3x^2y^2z + 3x^2y^4 + x^3y^3z\bigr)
$$
$$
\;+\; \bar y\,\bigl(3xz + 12xy^2 + 6x^2yz + 9x^2y^3 + 3x^3y^2z\bigr)
\;-\; \bar z\,\bigl(3x^2y + x^3z\bigr).
$$

Thirteen vertices in all, each with exactly one $\bar\varphi$-leg, ranging
from cubic ($\bar x y^2$, $\bar y xz$, $\bar z x^2 y$) to octic
($\bar x\, x^3 y^3 z$).

### $\mathbb{C}^*$-grading

The model is equivariant under the $\mathbb{C}^*$-action recorded in
`jcqft/core.py`: assigning weights $(1, -1, -2)$ to the fields $(x, y, z)$,
the three components of $F$ scale with weights $(-2, -1, 1)$, so the sources
$(a, b, c)$ carry weights $(-2, -1, 1)$ and the conjugate fields
$(\bar x, \bar y, \bar z)$ must carry the opposite weights $(2, 1, -1)$.
Every term of $S$ — each of the three $\bar\varphi J$ terms, each entry of
$\bar\varphi L \varphi$, and each of the thirteen vertices — then has weight
exactly $0$; the interaction preserves the grading vertex by vertex.

### What the tree expansion computes

In this language the one-point function
$\langle \varphi \rangle = G(J)$ is the formal inverse of $F$: expanding
$e^{-S_{\mathrm{int}}}$ and Wick-contracting with the propagator $L^{-1}$
generates precisely the sum over rooted trees of [BCW82, Wri87] — root leg
$L^{-1}$, internal lines $L^{-1}$, branchings given by the thirteen vertices
above, leaves given by source insertions $J$ — and there are no loop
corrections beyond the constant $\det DF = -2$, which cancels between
numerator and normalization. Everything established in `PROGRESS.md` about
the tree series (convergence, algebraicity, blindness to the two sheets at
infinity) is a statement about the correlation functions of this first-order
model.

## Goals of this project

1. **Diagnose the failure precisely in the QFT language.** Where exactly does
   the tree-graph expansion stop representing the physics? (Answered so far:
   it converges on a finite domain to a *local branch* of an algebraic
   function; the obstruction is *non-properness* — field configurations
   escaping to infinity — not divergence of the series. See `PROGRESS.md`.)
2. **Global, non-perturbative structure of the toy model.** The three fiber
   points ("sheets", loosely: competing saddle/vacuum configurations) and the
   monodromy that permutes them; the exact locus in source space where sheets
   escape to infinity; "vacua at infinity" as a non-perturbative mechanism
   distinct from instantons.
3. **New, independent counterexamples.** Reverse-engineer the construction
   mechanism of the Alpöge–Mathew map and search for genuinely inequivalent
   Keller-but-non-injective maps (not related by composition with linear or
   tame automorphisms).
4. **Implications for mathematically rigorous QFT**, in the spirit of
   constructive QFT (Glimm–Jaffe) and perturbative algebraic QFT
   (Fredenhagen–Rejzner): what does a tree-level, $D=0$ global failure teach
   us about formal-series constructions, field redefinitions, and the
   prospects for rigorous interacting QFT in $D = 4$? (With an honest
   delineation of what it does *not* teach: UV renormalization and the
   $D = 4$ constructive difficulties are a separate phenomenon.)

## Non-goals / cautions

- No claim that this resolves or directly attacks $\phi^4_4$ triviality,
  Yang–Mills existence, or any $D \ge 2$ constructive problem.
- Popular AI-generated summaries of the result circulate with errors (e.g.
  "the tree series has zero radius of convergence", "Borel resummation is
  impossible"). This project computes the true statements; corrections are
  recorded in `PROGRESS.md` and `docs/QFT_IMPLICATIONS.md`.
