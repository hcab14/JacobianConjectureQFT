# Speculative (pA)QFT / Gauge / BV / Amplitudes Implications

*(2026-08-05. Packaging document. **Every claim in this file is
Interpretive** unless it explicitly cites an Exact or Numerical result
already verified elsewhere in the repository. No continuum $D\ge 1$
theorem is asserted here. Companion Exact sources:
`docs/QFT_IMPLICATIONS.md`, `docs/AMPLITUDES_CONNECTION.md`,
`docs/BF_CARICATURE.md`, `docs/CLASSICAL_MAP_INVARIANTS.md` §§3–6,
`docs/D1_INDEX.md`, `docs/WITTEN_INDEX.md`, `docs/SYMMETRIC_SEARCH.md`,
`docs/TWISTED_PERIODS.md`. Paste-ready paper paragraphs:
`paper/ADDENDUM_DRAFT.md`.)*

---

## 0. Banner: epistemic tags

| Tag | Meaning in this document |
|---|---|
| **Exact** | Symbolic identity asserted in a named script; cited only as input. |
| **Numerical** | High-precision floating evidence with documented tolerances; not interval-certified. |
| **Interpretive** | Physics / (pA)QFT / amplitudes *reading* of Exact or Numerical facts. **Default tag for every paragraph below.** |

**Rule.** If a sentence does not name a script or an Exact/Numerical
document section, treat it as Interpretive. Overstatements already
corrected in `docs/QFT_IMPLICATIONS.md` §3 (zero radius, Borel,
"$D=4$ is harder") are not repeated here as claims.

**What this file is for.** The repository already separates Exact
computations from QFT readings (`QFT_IMPLICATIONS.md`,
`AMPLITUDES_CONNECTION.md`, paper §`sec:qftread`). This note collects
the *further* speculative layer — gauge theory, BV–BRST, lattice→QFT
limits, amplitudes beyond the verified pushforward, effective-action
geometry, anomalies in coupling space — and lists concrete theorems that
would promote each item to Exact.

---

## 1. Gauge theory: Gribov without horizon; background-field as Keller

**Exact inputs.** $\det DF\equiv -2$
(`scripts/verify_counterexample.py`); three real preimages of
$(-1/4,0,0)$; chamber rule $N(J)=3$ iff $p<0$
(`scripts/measure_anomaly.py`); no finite-distance Jacobian zero.

**Interpretive reading.**

1. **Gribov copies without a Gribov horizon.** Gauge fixing assumes a
   slice whose Faddeev–Popov determinant certifies *local* uniqueness.
   Classically, Gribov copies [Gri78, Sin78] are associated with the
   *horizon*, where that determinant vanishes. Alpöge–Mathew is the
   complementary pathology named in `QFT_IMPLICATIONS.md` §6.1: the
   FP-type determinant is the constant $-2$ — there is no horizon in
   finite field space — yet copies exist (the triple point). The horizon
   has retreated to infinity (Jelonek set $\{p=0\}$). Absence of a
   horizon does **not** exclude copies.
2. **Contrast, not a YM theorem.** Continuum Yang–Mills Gribov geometry
   has a horizon at finite field strength — the wrong failure mode
   relative to AM (`CLASSICAL_MAP_INVARIANTS.md` §4.5: reject YM as
   *primary* lift; retain as named contrast). Nothing here proves
   existence or absence of copies in continuum YM.
3. **Background-field gauge as a Keller problem.** In the
   background-field method one solves a nonlinear map
   $F_{\mathrm{bg}}(q)=J$ for the quantum fluctuation $q$ at fixed
   background. Local invertibility of $D_q F_{\mathrm{bg}}$ is the
   usual FP / elliptic diagnostic. The AM lesson (Interpretive): that
   diagnostic can pass everywhere while the global fiber remains
   multi-sheeted, with extra solutions entering only through infinity.
   Treating "invertible background-field operator $\Rightarrow$ unique
   quantum field" is the same false step as reading global invertibility
   off $\det DF=\mathrm{const}$.

**Strongest speculative claim of this section.** There exists a
structurally distinct class of gauge-fixing pathologies — *copies with
nowhere-vanishing FP determinant* — of which AM is the finite-dimensional
prototype; continuum models should be scanned for this class separately
from the horizon class.

---

## 2. BV–BRST: antifield pairing, master equation, BRST at infinity

**Exact inputs.** No potential in $n=3$: $\{K:K\,DF\text{ symmetric}\}=\{0\}$
(`scripts/symmetric_search.py`; `SYMMETRIC_SEARCH.md`). Cotangent lift
$W_6=\bar\varphi\cdot F(\varphi)$ has $\det\mathrm{Hess}W_6\equiv-4$
and is 3:1 (`SYMMETRIC_SEARCH.md` §3.1). First-order Lagrangian
$L=\bar\phi\cdot(F(\phi)-J)$ is the only honest implementation of
$F(\phi)=J$ in three fields (`BF_CARICATURE.md` §3; [Abd03]). Coercivity
unreachable for $W_6$ ($\kappa=-4$; affine in $\bar\varphi$;
`SYMMETRIC_SEARCH.md` §7). MQ / nilpotent BRST completion exists without
a superpotential (`WITTEN_INDEX.md`).

**Interpretive reading.**

1. **Antifield pairing is forced, not optional.** In BV language the
   conjugate field $\bar\phi$ is the antifield of $\phi$. The pairing
   $L=\bar\phi\cdot(F-J)$ is the 0D classical master-action density for
   the equation $F-J=0$. Because AM is non-variational (I8), *every*
   Lagrangian/BV transcription must pass through this doubling — there
   is no single-field master action for $F$ itself.
2. **Classical master equation vs global invertibility.** The classical
   master equation $(S,S)=0$ and the constant Jacobian (trivial
   one-loop/FP determinant) control *local* consistency of the gauge /
   equation complex. They do **not** imply that the antifield equations
   select a unique sheet. Exact parallel: $\det DF=\mathrm{const}$ and
   the tree inverse both pass, yet the fiber is 3-sheeted.
3. **BRST at infinity.** The Mathai–Quillen model supplies a nilpotent
   BRST charge $Q$ with $Q$-exact action even though no Parisi–Sourlas
   superpotential exists (`WITTEN_INDEX.md`). Localization computes the
   signed index $-N(J)$, which jumps across $\{p=0\}$. Interpretive
   slogan: *BRST localization sees the wall; it does not restore
   properness.* The escaping sheets contribute through the boundary of
   field space (on-wall limits $-2$ / cusp $-1$ in the 0D closed form),
   not through finite-distance fixed points of $Q$.
4. **Coercivity barrier for the honest action.** Promoting $W_6$ (or
   $L$) by kinetic terms does not repair $\kappa\le 0$ along conjugate
   directions (`CLASSICAL_MAP_INVARIANTS.md` §3.6, §4.6). Stationary-phase
   / BRST-localized formulations remain well-posed; absolute convergence
   of $\int e^{-W_6}$ does not.

**Strongest speculative claim.** In BV for non-gradient Keller data, the
master equation can hold with constant FP determinant while the
antifield fiber bundle still carries nontrivial holonomy and a
non-properness wall — global invertibility is an independent axiom on
the classical BV complex.

---

## 3. Lattice → QFT double limit; transfer-matrix bridge to D1

**Exact / Numerical inputs.**

- Ultralocal product: invariants tensor sitewise
  (`CLASSICAL_MAP_INVARIANTS.md` §6; `scripts/classical_map_invariants_probe.py`).
- Kinetic deformation $F_\varepsilon(\phi)_x=F(\phi_x)-\varepsilon(\Delta\phi)_x$:
  for $L=2$, $N_\varepsilon$ stays non-constant for probed
  $\varepsilon>0$ (**Numerical**, HC-certified;
  `scripts/lattice_chamber.py`); at $\varepsilon=1/4$ on the $T_1\to T_3$
  segment, all chamber walls are *fold*-type with exact fold polynomial
  of degree 516 (**Exact**, `scripts/lattice_discriminant.py`) —
  kinetic mixing *exchanges* the wall mechanism (escape → fold).
- Finite-mode $D=1$ MQ index: jump $-1\leftrightarrow-3$ survives for
  $M\le 2$ (**Exact** saddle factorization + **Numerical** MC;
  `docs/D1_INDEX.md`). $M=0$ reduces exactly to 0D MQ with
  $\sigma_{\mathrm{eff}}=\sigma/\sqrt\beta$.

**Interpretive reading — the double limit.**

```
ultralocal (ε=0)  →  kinetic lattice (ε>0, L finite)  →  continuum (L→∞, a→0)
         │                        │                            │
    Exact 0D⊗L              Exact/Num walls               OPEN
    product dictionary      (fold vs escape)              (no claim)
```

1. **Ultralocal → kinetic.** The 0D defect does not evaporate at the
   first kinetic deformation: chamber non-constancy survives, but the
   *mechanism* of walls changes (escape-dominated → fold-dominated on
   the probed segment). Interpretive: kinetics can convert
   "vacua at infinity" into ordinary finite-distance folds while
   preserving multi-valuedness — a different physical face of the same
   non-injectivity.
2. **Kinetic → continuum.** No continuum claim. The honest bridge is
   the finite-mode $D=1$ model of `D1_INDEX.md`: path-space kinetics
   enter only as positive spectator determinants at constant saddles,
   so the index jump is $M$-independent inside the truncation.
3. **Transfer-matrix bridge (Interpretive).** The $D=1$ periodic MQ
   model is the continuum avatar of a transfer-matrix / quantum-mechanics
   quantization of the Keller force map. The Exact $M=0$ reduction says
   the transfer-matrix ground-sector localization *contains* the 0D
   index; the saddle factorization says higher Matsubara modes cannot
   flip $\mathrm{sign}\det DG$ at equilibria. Promoting this to a
   continuum theorem (no nonconstant zeros; dominated convergence
   $M\to\infty$) would be the first Exact $D=1$ statement beyond
   truncation — listed in §8.

**Strongest speculative claim.** The non-properness defect of a Keller
force map is stable under the first two rungs of the lattice→QM ladder
(ultralocal product, weak kinetics / finite Fourier modes); the open
question is only whether the continuum rung kills it.

---

## 4. Amplitudes: Jelonek = second-type Landau; twisted periods; Stokes; $Z(J)$ as section

**Exact / Numerical inputs (cite, do not redo).**

- Landau dictionary: first-type singularities absent (étale);
  $\{p=0\}$ = second-type / escape (`AMPLITUDES_CONNECTION.md` §2.2;
  method of `scripts/branch_locus.py`; cf. [MT22]).
- Trace / pushforward rationality; poles only on $\{p=0\}$; boundary
  factorization (`scripts/trace_pushforward.py`,
  `scripts/pushforward_forms.py`).
- Chamber is **not** a positive geometry (cuspidal double pole;
  `docs/POSITIVE_GEOMETRY.md`).
- Twisted cohomology of the wall complement: dimensions, Alexander /
  Burau jump loci, integer-twist exactness of the would-be canonical
  form (`docs/TWISTED_PERIODS.md` — **cite only**).
- Geometric monodromy $S_3$ (**Numerical**, `scripts/monodromy.py`);
  Galois $S_3$ (**Exact**).

**Interpretive reading.**

1. **Jelonek = second-type Landau.** The non-properness divisor is the
   field-space avatar of Landau singularities at infinity: solutions
   escape in the eliminated variables while external parameters
   (sources) stay finite. AM is a laboratory where *only* second-type
   singularities exist and they control radius, monodromy, and measure
   anomaly.
2. **Twisted periods.** The one-dimensional spaces
   $H^{1,2}(M;\mathcal{S})$ of `TWISTED_PERIODS.md` are the natural
   homes for "master integrals" of wall-crossing data in the reflection
   channel. Interpretive: physical single-valued observables land in
   the trivial summand; sheet-separating data are twisted periods with
   poles / jumps on $\{p=0\}$. Integer twist makes the naive canonical
   form twisted-exact — the cohomological face of the positive-geometry
   failure.
3. **Stokes through infinity.** Analytic continuation of the
   perturbative branch around $\{p=0\}$ permutes sheets (Stokes /
   monodromy data = $S_3$). Because sheets never meet at finite
   distance, every nontrivial Stokes jump is mediated by infinity —
   "Stokes without Borel" in a convergent series.
4. **$Z(J)$ as a section of the fiber-algebra bundle.** Combining
   `BF_CARICATURE.md` with `DAMPED_PARTITION.md` / `WITTEN_INDEX.md`:
   the partition function (damped or MQ) is a single-valued
   *numerical* section built from the character theory of $A_J$
   (sheet count / signed count), while any *algebraic* sector-resolving
   section is multi-valued or polar on the wall. Interpretive: $Z(J)$
   is the canonical single-valued shadow of the fiber-algebra bundle —
   the AQFT datum of §4 of `BF_CARICATURE.md` evaluated in a state.

**Strongest speculative claim.** The physical generating function $Z(J)$
is best viewed as a flat/unitary section of a bundle of fiber algebras
(or of its character sheaf), with Stokes data the $S_3$ holonomy and
singular locus the Jelonek divisor — not as a single holomorphic
function of the sources with only finite-distance branch points.

---

## 5. Effective action: Maxwell-type vs Jelonek-type branch points of $\Gamma$

**Exact inputs.** Tree inverse = one algebraic branch; radius set by
$\{p=0\}$ escape, not by $D_0$-collisions
(`scripts/tree_expansion.py`, `scripts/branch_locus.py`).
$\det DF\equiv\mathrm{const}$ ⇒ no finite-distance Hessian
degeneration of a potential (and no potential exists in $n=3$).

**Interpretive reading.**

| Type | Local flag | Global meeting of branches | Present in AM? |
|---|---|---|---|
| **Maxwell-type** | Hessian / $\Gamma''$ degenerates at finite field | branches meet at finite $\bar\phi$ | **No** |
| **Jelonek-type** | none (Jacobian constant, étale) | branches meet only through infinity | **Yes** |

The 1PI effective action $\Gamma$ is defined by inverting
$J\mapsto\bar\phi(J)$ and Legendre transforming
(`AMPLITUDES_CONNECTION.md` §1.2). Formal and convergent constructions
silently select the perturbative sheet. Interpretive consequence for
pAQFT / background-field effective actions: *passing every local
convexity and invertibility test does not certify single-valuedness of
$\Gamma$.* Multi-branched Legendre transforms are familiar from Maxwell
construction; AM exhibits a second, horizonless branching mode.

**Strongest speculative claim.** Effective-action constructions need a
dichotomy in their global hypotheses: Maxwell-type (detectable by
$\Gamma''$) vs Jelonek-type (detectable only by properness / degree /
escape diagnostics).

---

## 6. Anomalies in the space of couplings ($S_3$ holonomy)

**Exact / Numerical inputs.** Geometric monodromy $=S_3$ (Numerical);
Galois $=S_3$ (Exact); wall complement $\simeq A_2$ discriminant
complement with $\pi_1=B_3\twoheadrightarrow S_3$
(`docs/WALL_COMPLEMENT.md`); measure anomaly $A(\sigma)\to 2$ as
$\sigma\to 0$ (Numerical on Exact chamber rule).

**Interpretive reading.**

Not a cohomological anomaly (no broken symmetry, not loop-generated —
`QFT_IMPLICATIONS.md` §6.2). Pattern match to *anomalies in the space of
coupling constants* [CFLS20]: transporting the external source
(= coupling) around a closed cycle returns the same local theory with
vacua permuted. Here the transport is entirely through infinity; there
are no massless states and no finite-distance degeneration.

In the BF caricature this is literally the holonomy of relative
$S$-transport on the fiber-algebra bundle (`BF_CARICATURE.md` §4.3):
wall meridians = transpositions; empty-fiber cusp loop = Coxeter
3-cycle. Interpretive name already in the repo: **non-properness
defect** / parameter-space monodromy.

**Strongest speculative claim.** $S_3$ holonomy of Keller counterexamples
is a tree-level, algebraically exact instance of an anomaly in coupling
space, with the Jelonek divisor as its defect locus.

---

## 7. Explicit non-goals

The following are **out of scope** for every Interpretive claim in this
file (restating `QFT_IMPLICATIONS.md` §4.2 and
`CLASSICAL_MAP_INVARIANTS.md` §3.7):

| Non-goal | Why |
|---|---|
| UV renormalization / counterterms / running | Model is $D=0$ or finite-mode / finite lattice; no short-distance singularities. |
| Borel summability of $\phi^4_2$, $\phi^4_3$ | AM tree series *converges*; loop/large-order divergence is a disjoint phenomenon. |
| Existence, triviality, or mass gap of interacting $D=4$ QFT | Ultraviolet and measure-theoretic; logically independent of global field-space geometry at tree level. |
| Continuum $M\to\infty$ or $a\to 0$ limits | Explicitly open in `D1_INDEX.md` and lattice §§. |
| Continuum Yang–Mills Gribov problem | Contrast only (§1). |
| Positive geometry of the $N=3$ chamber | Already settled **negatively** (`POSITIVE_GEOMETRY.md`); not re-opened here. |

---

## 8. Concrete next theorems (Interpretive → Exact)

Each item names a *theorem-shaped* statement that would promote a
section above from Interpretive packaging to an Exact (or certified
Numerical) repository result.

| # | Proposed theorem | Would promote | Suggested attack |
|---|---|---|---|
| T1 | **No nonconstant zeros** of the truncated AM flow $G_{J,\beta,M}$ for all $M$, or a sharp bound on their contribution to $Z_\sigma$ | §3 (D1 bridge); closes the AM loophole in `D1_INDEX.md` | Algebraic degree bounds on $G$; resultant criteria; or interval-Newton certification beyond multi-start |
| T2 | **Continuum $D=1$ index jump:** $\lim_{M\to\infty}Z_\sigma(J;M,\beta)=\deg(F,J)$ off the wall for small $\sigma$ | §3 strongest claim | Dominated convergence for mode tails + T1 |
| T3 | **Lattice wall hypersurface** $W_\varepsilon(a,b,c)=0$ for frozen second site, with irreducible factors separating fold vs escape | §3 kinetic mechanism | Extend `lattice_discriminant.py` beyond the segment (needs more GB / modular methods) |
| T4 | **Q2b at $\varepsilon>0$:** semiclassical prefactor of $Z_{\hbar,\varepsilon}$ equals $N_\varepsilon/|\kappa|^L$ including escaping mass | §3 measure half of C7 | Lattice MQ MC with escape-coordinate mass diagnostics (as in `D1_INDEX.md` §5) |
| T5 | **Certified geometric monodromy** $=S_3$ | §4 Stokes; §6 anomaly | ✅ Done: OPEN_QUESTIONS B6 / `scripts/certified_monodromy.py` (Exact Puiseux + Gal; not interval tracking) |
| T6 | **Period pairing theorem:** a basis of $H^1(M;\mathcal{S})$ whose integrals recover wall-crossing of $T[x^k]$ or separator residues | §4 twisted periods | Build on `TWISTED_PERIODS.md` + `trace_pushforward.py` boundary factorization |
| T7 | **BV master identity** for $L=\bar\phi\cdot(F-J)$ in the finite-dimensional BV complex: $(S,S)=0$, Koszul homology $=A_J$, and holonomy obstruction as in `BF_CARICATURE.md` §5 restated in antibracket language | §2 | Homological algebra on the Koszul complex of $(F-J)$; mostly reorganization of Exact facts |
| T8 | **Effective-action dichotomy lemma** (0D): for polynomial Keller $F$, either $\Gamma$ is single-valued polynomial (automorphism case) or every global branch of $\Gamma$ has a Jelonek-type singularity on $S_F$ with no Maxwell point | §5 | Follows from Jelonek + étaleness once $\Gamma$ is defined sheetwise via Legendre transform on each branch |
| T9 | **Background-field Keller criterion:** a finite-dimensional background-field map with $\det D_q F_{\mathrm{bg}}=\mathrm{const}\neq 0$ is globally invertible iff it is proper | §1 | Special case of Hadamard–Lévy / Hadamard global inversion for local diffeomorphisms — cite and specialize |
| T10 | **Minimal variational dimension:** no gradient Keller counterexample in dimension $<6$ (or an explicit example in 4 or 5) | §2 coercivity / $W_6$ | Extend `symmetric_search.py` boxes; dBvdE dimension accounting |

**Minimal promotion path for a paper addendum.** T7 + T8 are largely
reorganization of existing Exact results (high value / low risk). T1–T2
are the scientific bottleneck for the $D=1$ story. T5 is the
certification bottleneck for every Stokes / anomaly paragraph.

---

## 9. Pointers

| Topic | Primary Exact doc | Paper / draft |
|---|---|---|
| Gauge / Gribov contrast | `QFT_IMPLICATIONS.md` §6.1 | `ADDENDUM_DRAFT.md` (context only) |
| BV / antifield / $W_6$ | `SYMMETRIC_SEARCH.md`, `BF_CARICATURE.md` §3 | `ADDENDUM_DRAFT.md` § BV |
| Lattice + D1 | `CLASSICAL_MAP_INVARIANTS.md` §§5–6, `D1_INDEX.md` | `ADDENDUM_DRAFT.md` § Møller after kinetics |
| Amplitudes / Landau | `AMPLITUDES_CONNECTION.md`, `TWISTED_PERIODS.md` | paper amplitudes remarks; cite twisted periods |
| Fiber-algebra AQFT datum | `BF_CARICATURE.md` | `ADDENDUM_DRAFT.md` § fiber-algebra bundle |
| Non-goals ledger | `QFT_IMPLICATIONS.md` §4.2–4.4 | paper §`sec:qftread` |

---

## References (external; already used in-repo)

No new citations invented. Standard anchors already in
`QFT_IMPLICATIONS.md` / `AMPLITUDES_CONNECTION.md` / `BF_CARICATURE.md`:
[Abd03], [BF20], [CFLS20], [CHY14], [FR12, FR16, Rej16], [Gri78],
[MT22], [Nic80], [Sin78], [SW94].
