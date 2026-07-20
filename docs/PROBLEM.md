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
