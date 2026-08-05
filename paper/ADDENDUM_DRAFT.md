# Addendum draft: paste-ready LaTeX subsections

*(2026-08-05. **Do not treat this file as part of the compiled paper.**
These are draft subsections intended for manual pasting into
`paper/main.tex` after editorial review. Every continuum / gauge / BV
reading below carries an explicit house disclaimer. No fake citations:
only anchors already used in `main.tex` or the repository docs.
Companion interpretive packaging: `docs/SPECULATIVE_IMPLICATIONS.md`.)*

**How to use.** Each block below is a self-contained `\subsection{...}`
(or `\paragraph`) in the paper's voice. Suggested homes:

| Draft subsection | Suggested paste location |
|---|---|
| Møller invertibility after kinetics | after `\ref{sec:d1index}` / end of lattice+$D=1$ block |
| BV antifield reading | after variational / $W_6$ material, or after `\ref{sec:bfcar}` |
| Fiber-algebra bundle as 0D AQFT datum | expand or follow `\ref{sec:bfcar}` / `\ref{sec:bf}` |

**House disclaimers (include when pasting).** These paragraphs calibrate
structural analogies; they are not theorems about $D\ge 1$ nets, continuum
path integrals, or Yang–Mills. UV renormalization, Borel summability, and
$D=4$ existence are out of scope (cf.\ Section~\ref{sec:qftread}).

---

## Ready to paste (1): Møller invertibility after kinetics

```latex
\subsection{M{\o}ller invertibility after kinetics
(lattice and finite-mode $D=1$)}
\label{sec:moller-kinetics}

% HOUSE DISCLAIMER: structural calibration inside verified truncations;
% no continuum (M\to\infty / a\to0) claim.

In pAQFT the classical M{\o}ller map inverts the nonlinear field
equation as a (typically formal) power series
(Section~\ref{sec:paqft}). The zero-dimensional calibration is already
sharp: that inverse may converge and still miss sheets at infinity.
The first honest question off $D=0$ is whether \emph{kinetics} ---
lattice gradients or a $D=1$ path measure --- restore global
invertibility of the classical map, or kill the escaping equilibria
that drive the index jump of Section~\ref{sec:witten}. Within the
truncations computed in this repository, the answer is negative on
both counts: kinetics do not restore a single-sheeted M{\o}ller map.

\paragraph{Lattice.}
On an $L$-site chain the ultralocal product $F^{\times L}$ tensors the
classical-map invariants of Section~\ref{sec:dictionary} exactly
(\texttt{scripts/classical\_map\_invariants\_probe.py}). The kinetic
deformation
\[
\bigl(F_\varepsilon(\phi)\bigr)_x
\;=\;
F(\phi_x)\;-\;\varepsilon\,(\Delta_{\mathrm{disc}}\phi)_x
\]
was probed for $L=2$: the real chamber function $N_\varepsilon$ remains
non-constant for every tested $\varepsilon>0$
(\texttt{scripts/lattice\_chamber.py}, homotopy-continuation
certification of reality), so $F_\varepsilon$ stays non-injective after
mixing. At $\varepsilon=1/4$ on the rational segment joining the
one-site chamber points $T_1\to T_3$, the discriminant computation of
\texttt{scripts/lattice\_discriminant.py} gives an exact fold
eliminant of degree $516$ whose $14$ real roots in $(0,1)$ account for
\emph{all} chamber walls on the segment (every jump $\pm 2$); escape
contributes only pointwise dips of the count. Thus weak kinetics
preserve multi-valuedness while exchanging the wall mechanism
(escape-type at $\varepsilon=0$ versus fold-type at
$\varepsilon=1/4$). A formal or convergent M{\o}ller series about the
ultralocal vacuum therefore continues, after kinetics, to describe at
most one local sheet of a still multi-sheeted classical map.

\paragraph{Finite-mode $D=1$.}
The Mathai--Quillen index of the first-order flow
$\dot q+F(q)-J$ on Fourier-truncated periodic paths
(Section~\ref{sec:d1index}) jumps with the zero-dimensional degree:
at constant-path equilibria the mode fluctuation determinant
factorizes as a $J$-independent positive spectator
$\prod_k\lvert\det(DF+i\omega_k I)\rvert^2$, so
\[
\mathrm{sign}\det DG(u^*)
\;=\;
\mathrm{sign}\det DF(q^*)
\]
and the $\sigma\to0$ saddle sum remains $\deg(F,J)=-N(J)$, independent
of mode cutoff $M$ and inverse temperature $\beta$ (symbolic asserts in
\texttt{scripts/d1\_index\_modes.py}). Monte-Carlo integration confirms
the chamber ratio $3{:}1$ for $M\le 2$. Path-space kinetics, in this
truncation, do not suppress the vacua at infinity and do not
single-value the classical inverse.

\paragraph{Reading for pAQFT.}
Upgrading a formal M{\o}ller map to a convergent one is still not the
missing step toward a globally single-valued classical inverse: after
the first kinetic deformations one still needs an independent
properness (or degree) input. We emphasize the scope: no statement is
made about the continuum limits $M\to\infty$ or lattice spacing
$a\to0$, nor about interacting theories in $D\ge 2$.
```

---

## Ready to paste (2): BV antifield reading of the first-order action

```latex
\subsection{BV antifield reading of the first-order action}
\label{sec:bv-antifield}

% HOUSE DISCLAIMER: finite-dimensional BV vocabulary for Exact
% identities already in Sections on the gradient no-go, W_6, and the
% BF caricature. Not a theorem about continuum BV--BRST.

The Alp\"oge--Mathew map admits no potential in three fields
(Proposition~\ref{prop:nogradient}): there is no $W$ with
$\nabla W=F$. The only Lagrangian implementing the classical equation
$F(\phi)=J$ is therefore the first-order (cotangent) density
\[
L\bigl(\bar\phi,\phi;J\bigr)
\;=\;
\bar\phi\cdot\bigl(F(\phi)-J\bigr),
\]
already used by Abdesselam in the QFT reading of the Jacobian
conjecture \cite{Abd03} and forced here by the variationality flag I8
of Section~\ref{sec:dictionary}. In Batalin--Vilkovisky language
$\bar\phi$ is the antifield of $\phi$, and $L$ is the classical
master-action density for the equation $F-J=0$. Two Exact faces of the
same object appear elsewhere in the paper:

\begin{itemize}
\item \emph{Variational counterexample in dimension $6$.}
  The cotangent lift $W_6=\bar\varphi\cdot F(\varphi)$ is an honest
  gradient Keller map,
  $\det\mathrm{Hess}W_6\equiv-4$, with explicit rational
  $3{:}1$ witnesses (Section~\ref{sec:witten},
  \texttt{docs/SYMMETRIC\_SEARCH.md}). The doubling that BV requires
  for non-gradient data is thus also the constructive route to a
  symmetric Jacobian counterexample. Coercivity is nevertheless
  unreachable ($\kappa=-4$; $W_6$ is affine in $\bar\varphi$), so
  $\int e^{-W_6}$ is never absolutely convergent --- stationary-phase
  and BRST-localized formulations remain the well-posed ones.
\item \emph{Dynamical relation in the $0$D BF caricature.}
  Antifield shifts give
  $\delta L(\beta)=\beta\cdot(F(\phi)-J)$ with no higher corrections
  ($L$ is affine in $\bar\phi$), and multiplicative evaluations
  compatible with the relation are precisely the characters of the
  fiber algebra $A_J=\C[x,y,z]/(F-J)$
  (Section~\ref{sec:bfcar}).
\end{itemize}

\paragraph{Master equation versus global invertibility.}
The constant Jacobian $\det DF\equiv-2$ makes the one-loop
(Faddeev--Popov / Berezin) determinant field-independent, and the
classical master equation for $L$ is the Koszul consistency of
$F-J=0$ --- both are local. Neither forces the antifield equations to
select a unique sheet: the fiber of $F$ is generically three points,
the geometric monodromy is $S_3$, and the non-properness wall
$\{p=0\}$ is the jump locus of $\dim A_J$. In slogans adapted to BV:
\emph{the master equation can hold with constant FP determinant while
the antifield fiber bundle still carries nontrivial holonomy.} Global
invertibility of the classical BV complex is an independent axiom,
equivalent here to properness of $F$ (or constancy of the Brouwer
degree). The Mathai--Quillen BRST charge of Section~\ref{sec:witten}
localizes to the signed count $-N(J)$ and \emph{sees} the wall; it
does not remove it.

We claim nothing about the existence of a continuum BV action with
these properties, nor about anomalies of the BV complex in $D\ge 1$.
The finite-dimensional lesson is only that antifield consistency and
global single-valuedness of the classical inverse are separate
inputs.
```

---

## Ready to paste (3): Fiber-algebra bundle as 0D AQFT datum

```latex
\subsection{The fiber-algebra bundle as a $0$D AQFT datum}
\label{sec:fiber-bundle-aqft}

% HOUSE DISCLAIMER: 0D caricature of Buchholz--Fredenhagen data
% (Section~\ref{sec:bfcar}); not a construction of a Haag--Kastler net
% and not a theorem about \cite{BF20}.

Haag--Kastler AQFT takes algebras of observables, rather than field
coordinates, as primary data (Section~\ref{sec:aqft}). In $D=0$ there
is no net and causal factorization trivializes
(Section~\ref{sec:bfcar}). What survives of a Buchholz--Fredenhagen-style
assignment $J\mapsto S(J)$ is nevertheless an algebraic object with a
precise dictionary to the classical-map invariants of
Section~\ref{sec:dictionary}: the \emph{bundle of fiber algebras}
\[
J\;\longmapsto\;
A_J\;=\;\C[x,y,z]\big/(F-J)
\]
over source space, together with the parallel transport of its local
system along paths of sources (``relative $S$-elements'').

\paragraph{Where the invariants sit.}
Off the non-properness wall $\{p=0\}$ one has $\dim_\C A_J=3$; on a
generic wall point the dimension drops to $1$; over the empty-fiber
cusp orbit $A_J$ is the zero ring. Thus $\{p=0\}$ is the non-flat
locus of the bundle --- the module-theoretic face of non-properness
already visible in the failure of $x$ to be integral over
$\mathrm{im}F^*$ (Section~\ref{sec:obs}). Over $\R$, off the
wall, $A_J\cong\C^{N(J)}$ as a finite-dimensional commutative
$C^*$-algebra, so the chamber function $N(J)$ is the character count.
The $S_3$ monodromy is the holonomy of the transport (wall meridians
act by transpositions; the cusp loop by the Coxeter $3$-cycle); there
is provably no global deck action, because the degree-$3$ extension is
non-normal. Single-valued sheet-separating data necessarily introduce
poles along $\{p=0\}$ (separator coefficients carry the factor $p$).

\paragraph{Obstruction dichotomy.}
Any assignment $J\mapsto S(J)$ that satisfies the $0$D dynamical
relation and resolves classical sectors is therefore multi-valued
across loops around the wall, or singular on the wall; everything
single-valued and pole-free pulls back through $F^*$ and is blind to
the sectors. For a proper Keller control (a tame automorphism) the
same construction collapses to the trivial rank-$1$ bundle --- the
caricature re-expresses non-properness as
$(\text{rank jump},\,\text{holonomy},\,\text{pole divisor})$ rather
than assuming it.

\paragraph{Partition functions as single-valued shadows.}
The damped and Mathai--Quillen integrals of
Sections~\ref{sec:witten}--\ref{sec:d1index} supply canonical
\emph{numerical} sections built from the character theory of $A_J$
(unsigned count $N(J)$, signed count $-N(J)$), finite for every source
and blind to algebraic sector choice. In this sense the physical
generating function $Z(J)$ is the single-valued shadow of the
fiber-algebra bundle --- an AQFT-flavoured packaging of the same Exact
defect that the observable pullback $F^*$ detects algebraically.

This remains a $0$D caricature: no spacetime net, no causality, no
Hilbert-space dynamics. Its value is dictionary-shaped. It isolates
the algebraic slots that any non-perturbative framework must fill if
it is to represent classical dynamics whose solution count jumps
across an algebraic wall with extras entering from infinity.
```

---

## Optional glue paragraph (Interpretive; use only with §0 banner)

If the addendum opens a short interpretive block, the following
one-paragraph banner matches `docs/SPECULATIVE_IMPLICATIONS.md` §0:

```latex
\paragraph{Epistemic status.}
The three subsections below reorganize Exact and numerically supported
identities already verified in the companion repository
(\texttt{docs/D1\_INDEX.md}, \texttt{docs/BF\_CARICATURE.md},
\texttt{docs/SYMMETRIC\_SEARCH.md},
\texttt{docs/CLASSICAL\_MAP\_INVARIANTS.md}). Continuum limits,
Yang--Mills Gribov geometry, and $D=4$ constructive existence are
explicit non-goals; where a sentence extrapolates beyond a named
truncation it is interpretive calibration in the sense of
Section~\ref{sec:qftread}, not an additional theorem.
```

---

## Checklist before pasting into `main.tex`

- [ ] Confirm `\label{prop:nogradient}`, `\ref{sec:paqft}`,
      `\ref{sec:dictionary}`, `\ref{sec:bfcar}`, `\ref{sec:obs}`,
      `\ref{sec:witten}`, `\ref{sec:d1index}`, `\ref{sec:aqft}`,
      `\ref{sec:qftread}` still match the compiled paper.
- [ ] `\cite{Abd03}` and `\cite{BF20}` already in the bibliography.
- [ ] Decide whether lattice fold-degree $516$ and HC chamber results
      are already summarized earlier (avoid duplication; cross-ref
      instead of repeating tables).
- [ ] Keep the three house-disclaimer comments (or promote them to
      visible footnotes).
- [ ] Do **not** add citations not already vetted in
      `docs/QFT_IMPLICATIONS.md` §7 / the paper `.bib`.
