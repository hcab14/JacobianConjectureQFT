# The Wall Complement Is the Braid-Group Classifying Space: the A2-Discriminant Structure of the Counterexample

*(2026-07-21. All exact claims verified in `scripts/wall_braid.py` (symbolic
identities + numerical labelled monodromy; runs in ~6 s); this answers item
**B2** of `docs/OPEN_QUESTIONS.md` — identify the correct
"amplituhedron-analogue" object after the failure of the real-chamber
positive geometry (`docs/POSITIVE_GEOMETRY.md`). Repo-internal inputs: the
C\*-reduction and cusp geometry of `docs/POSITIVE_GEOMETRY.md`, the $S_3$
monodromy of `docs/MONODROMY.md` (group-level statement Exact as of B6 /
`scripts/certified_monodromy.py`), the trace rationality of
`scripts/trace_pushforward.py`.)*

**Summary of the main results.**

1. **The model is affinely THE universal depressed cubic.** In the
   C\*-invariants $(u,w) = (ac^2, bc)$ the eliminant becomes
   $P_2\,\xi^3 + (4-3w)\,\xi - 2$, and the explicit affine isomorphism
   $(Q,R) = \bigl(w - \tfrac43,\; 2u - \tfrac{2w}{3} + \tfrac{16}{27}\bigr)$
   carries the invariant plane onto the space of depressed cubics
   $x^3 + Qx + R$ with $4Q^3 + 27R^2 = 4P_2$ **exactly, with multiplicity
   one**: the wall $\{P_2 = 0\}$ *is* the $A_2$ discriminant
   $\{4Q^3+27R^2=0\}$, the cusp goes to the origin, and the $D_0$-line to
   the cuspidal tangent $\{R=0\}$.
2. **Therefore $\pi_1(\mathbb{C}^2 \setminus \{P_2=0\}) = B_3$**, the braid
   group on three strands, and the wall complement is a $K(B_3,1)$
   [Arn69, Bri71, Del72]. The $S_3$ sheet monodromy is the image of the
   canonical surjection $B_3 \twoheadrightarrow S_3$ (**Exact**, B6).
3. **The cusp-loop monodromy is a 3-cycle of order 3** (Exact via local
   model; Numerical path labels in §3): the image of the braid
   $\sigma_1\sigma_2$, a Coxeter element of $W(A_2) = S_3$ — *not* a
   transposition. The full twist $(\sigma_1\sigma_2)^3 = \Delta^2$ (a
   torus-type loop) acts trivially, as it must.
4. **The rank-3 sheet local system splits as trivial $\oplus$ standard.**
   The trivial summand is exactly the trace observables (rational, poles
   only on the wall); *all* non-rational sheet data lives in the 2-dim
   standard local system, the reflection representation of $W(A_2)$. The
   wall complement has Euler characteristic $0$.
5. **Proposal** (clearly marked as such in §5): the pair (wall complement,
   standard local system) — equivalently the $K(B_3,1)$ structure with its
   reflection representation — is the correct replacement for the failed
   real-chamber positive geometry; twisted periods over it play the role of
   canonical-form integrals.

---

## 1. The invariant eliminant and the tautological map to depressed cubics

All statements in this section are exact polynomial identities, verified in
`scripts/wall_braid.py` §1.

Under the C\*-action with source weights $(a,b,c) \sim (-2,-1,1)$ and root
weight $X \sim 1$, the eliminant $pX^3 + qX + r$ is quasi-homogeneous of
weight $1$ (all seven monomials; verified). Substituting the weight-0 root
variable $\xi = X/c$ and dividing by $c$ therefore lands in the invariants
$u = ac^2$, $w = bc$:

$$
\textbf{(I1)}\qquad
\frac{pX^3 + qX + r}{c}\Big|_{X = \xi c}
\;=\; P_2(u,w)\,\xi^3 + (4 - 3w)\,\xi - 2,
$$

with $P_2 = c^2 p = 27u^2 + 16u - 18uw + w^3 - w^2$ the plane wall of
`docs/POSITIVE_GEOMETRY.md`, and $q = 4-3bc = 4-3w$, $r = -2c$. The
three-dimensional covering data reduces *entirely* to this one-parameter
cubic family over the $(u,w)$-plane; its discriminant in $\xi$ is
$-4\,P_2\,(27u-9w+8)^2$, the invariant form of
$\mathrm{disc}_X = -4\,D_0^2\,p$ from `docs/MONODROMY.md`.

Depressing by $\eta = P_2\,\xi$ (multiply (I1) by $P_2^2$) gives a *monic*
depressed cubic whose coefficients are **polynomial** in $(u,w)$:

$$
\textbf{(I2)}\qquad
\eta^3 + \hat{Q}\,\eta + \hat{R},
\qquad (\hat{Q}, \hat{R}) = \bigl((4-3w)\,P_2,\; -2\,P_2^2\bigr),
$$

a tautological map from the invariant plane to the universal space
$\{x^3 + Qx + R\}$ of depressed cubics. Its pullback of the universal
discriminant $4Q^3 + 27R^2$ factors exactly:

$$
\textbf{(I3)}\qquad
4\hat{Q}^3 + 27\hat{R}^2 \;=\; 4\,P_2^3\,(27u - 9w + 8)^2,
\qquad\text{driven by the lemma}\quad
(4-3w)^3 + 27P_2 = (27u - 9w + 8)^2 .
$$

Interpretation. On the wall $\{P_2 = 0\}$ we get $\hat Q = \hat R = 0$
*exactly*: the tautological map **contracts the whole wall to the cusp**
$(Q,R) = (0,0)$ of the universal discriminant (the scaling $\eta = P_2\xi$
degenerates there — this is the algebraic face of "two roots escape, and
after rescaling all three collide at $0$"). The $D_0$-line accounts for the
second factor of (I3): on $\{27u-9w+8 = 0\}$, off the wall, the image lies
on the universal discriminant away from the origin (a genuine double root
of the monic cubic — the x-collision of `docs/MONODROMY.md`, which the
covering itself does not feel). So (I2) is natural but crushes the wall; the
faithful identification is the affine one of §2.

## 2. The explicit affine equivalence with the A2 discriminant

The reduced wall and the $A_2$ discriminant curve
$\{W^2 = U^3\}$ are both cuspidal cubics whose unique flex tangent is the
line at infinity (our curve meets infinity only at $[1:0:0]$ with contact
3 — `docs/POSITIVE_GEOMETRY.md` §1), so an *affine* equivalence can exist.
It does, and it is rational and explicit (`scripts/wall_braid.py` §2; the
ansatz $P_2 = K(W^2 - U^3)$ with affine-linear $U, W$ forces $U$ free of
$u$, and with the normalization $K = 27$ has exactly two rational
solutions, $(U, W)$ and $(U, -W)$):

$$
\textbf{(I4)}\qquad
U = \frac{4 - 3w}{9}, \qquad W = \frac{27u - 9w + 8}{27}
\qquad\Longrightarrow\qquad
P_2 \;=\; 27\,\bigl(W^2 - U^3\bigr)
\quad\text{(exact identity).}
$$

The universal coordinates are thus **the eliminant coefficient $q$ and the
collision factor $D_0$ themselves**, up to scale: $U = q/9$, $W = D_0'/27$
(primes denoting the plane reductions). The lemma in (I3) *is* (I4)
restated. Equivalently, in the coordinates of the universal family:

$$
(Q, R) \;=\; (-3U,\; 2W) \;=\;
\Bigl(w - \tfrac43,\;\; 2u - \tfrac{2w}{3} + \tfrac{16}{27}\Bigr),
\qquad
4Q^3 + 27R^2 \;=\; 4\,P_2 .
$$

This is an affine **isomorphism** $\mathbb{C}^2_{(u,w)} \to
\mathbb{C}^2_{(Q,R)}$ (linear part of determinant $-2 \neq 0$) pulling the
universal discriminant back to $4 P_2$ — multiplicity **one**, so it maps
the pair (plane, wall) isomorphically onto the pair (plane, $A_2$
discriminant). Checkpoints (all exact): cusp $(4/27, 4/3) \mapsto (0,0)$;
$D_0$-line $\mapsto \{R = 0\}$, the cuspidal tangent of the universal
discriminant; $\{w = 4/3\} \mapsto \{Q = 0\}$.

Two remarks on what (I4) is and is not:

- It is an isomorphism of *pairs of affine varieties*, not of cubic
  families: the fiberwise cubic of (I1)/(I2) over $(u,w)$ is **not** the
  universal cubic over the image point (their discriminants differ by the
  factor $P_2^2 D_0'^2$). The two structures are complementary: (I2) is the
  tautological family map (wall $\to$ cusp), (I4) the topological
  identification (wall $\xrightarrow{\sim}$ discriminant).
- Nothing here required a projective transformation: the equivalence is
  affine because the flex-tangent-at-infinity configurations match. This is
  why the *affine* complement inherits the full braid-group structure.

**Consequence.** The wall complement is homeomorphic (indeed affinely
isomorphic) to the complement of the $A_2$ discriminant — the space of
depressed cubics with three distinct roots, i.e. the space of regular
orbits of the reflection group $W(A_2) = S_3$. By Arnold [Arn69], Brieskorn
[Bri71], and Deligne [Del72] this space is a $K(\pi,1)$ with

$$
\pi_1\bigl(\mathbb{C}^2 \setminus \{P_2 = 0\}\bigr)
\;=\; B_3
\;=\; \langle \sigma_1, \sigma_2 \mid
\sigma_1\sigma_2\sigma_1 = \sigma_2\sigma_1\sigma_2 \rangle ,
$$

the braid group on three strands (equivalently the trefoil-knot group;
Brieskorn's theorem covers all finite complex reflection groups, Deligne
all simplicial arrangements — for $A_2$ the statement goes back to the
braid-group case of Fox–Neuwirth/Arnold). The sheet monodromy of the
covering, verified to be the full $S_3$ in `docs/MONODROMY.md` and
re-verified in the invariant plane here, is then necessarily the image of
the canonical surjection $B_3 \twoheadrightarrow S_3$,
$\sigma_i \mapsto (i\;\, i{+}1)$: meridians of the wall map to
transpositions, and our model realizes the universal picture with no
deformation at all.

For contrast (and to forestall a common confusion): the fundamental group
of the **projective** complement of a cuspidal cubic is
$\pi_1(\mathbb{P}^2 \setminus C) = \mathbb{Z}/2 * \mathbb{Z}/3$
[Zar29] — the quotient of $B_3$ by its center
$\langle (\sigma_1\sigma_2)^3 \rangle$ (i.e. $\mathrm{PSL}_2(\mathbb{Z})$),
the center being generated by the class of a loop around the line at
infinity. The affine statement, $B_3$ itself, is the one relevant here: the
invariant plane is honestly affine, and the full twist acts nontrivially in
it (as a braid, though trivially in $S_3$; see §3).

## 3. Numerical monodromy in the invariant plane

`scripts/wall_braid.py` §3 tracks the three roots of the invariant
eliminant $P_2\,\xi^3 + (4-3w)\,\xi - 2 = 0$ (mpmath, 30 digits;
solve-and-match-by-nearest-neighbour continuation with adaptive bisection
and a 0.35-of-separation jump guard, plus an unambiguity test at the final
matching) around three kinds of loops. Results:

| loop | locus | permutation | order |
|---|---|---|---|
| generic line, meridian of wall pt 1 | $P_2 = 0$ | (2 3) | 2 |
| generic line, meridian of wall pt 2 | $P_2 = 0$ | (1 2) | 2 |
| generic line, meridian of wall pt 3 | $P_2 = 0$ | (1 3) | 2 |
| generic line, meridian of $D_0$-line pt | $D_0' = 0$, $P_2 \neq 0$ | id | 1 |
| **radius-0.05 loop around the cusp** $(4/27, 4/3)$ | — | **(1 2 3)** | **3** |
| torus loop $(Q,R) = (\varepsilon^2 e^{2it}, \varepsilon^3 e^{3it})$ pulled back | — | id | 1 |

(Generic line $(u,w) = (-1/2, 1/3) + s\,(1, 2+i)$, checked square-free and
disjoint from the $D_0$ point; wall meridians generate the full $S_3$,
order 6. The cusp loop is a circle of radius $0.05$ in the line through the
cusp with direction $(1 + 0.3i,\, -0.7 + 1.1i)$, which meets the wall only
at the cusp ($s=0$, double) and at a far point $|s| \approx 18.6$;
$\min |P_2| \approx 0.10$ on the circle.)

**The cusp-loop monodromy is a 3-cycle, of order 3.** Honest reporting
against the prior guesses: it is *not* the image of
$\sigma_1\sigma_2\sigma_1$ (the half-twist $\Delta$, whose image would be a
transposition of order 2); it is the image of $\sigma_1\sigma_2$, a
**Coxeter element** of $W(A_2)$, whose order is the Coxeter number
$h = 3$. Three independent confirmations, all in the script:

1. *Local model.* Near the cusp the eliminant degenerates to
   $P_2\,\xi^3 = 2$ (since $q \to 0$ there too — the cusp is the
   total-escape point, `docs/POSITIVE_GEOMETRY.md` §2), and $P_2$ winds
   **twice** around a small cusp loop (its quadratic part is the perfect
   square $3(3\delta u - \delta w)^2$); but $\xi \sim (2/P_2)^{1/3}$ has
   three determinations that a single winding of $P_2$ already permutes
   cyclically. (A double winding of $P_2$ gives the *square* of a 3-cycle —
   still a 3-cycle.)
2. *Splitting the cusp.* On the perturbed line through
   $(4/27,\, 4/3 + 10^{-2})$ the cusp splits into two nearby wall points
   with the $D_0$ point *between* them; their meridians are the distinct
   transpositions $(1\,2)$ and $(1\,3)$, the $D_0$ meridian is id, and the
   radius-0.05 loop around all three is $(1\,2\,3)$ = the product of the
   two transpositions. In $B_3$: the cusp loop is (conjugate to)
   $\sigma_1\sigma_2$.
3. *Universal family side-by-side.* The same tracker on $x^3 + Qx + R$
   around the cusp of $\{4Q^3 + 27R^2 = 0\}$ gives a 3-cycle of order 3 as
   well; and the torus-type loop
   $(Q,R) = (\varepsilon^2 e^{2it},\, \varepsilon^3 e^{3it})$, $t: 0 \to
   2\pi$ — the **full twist** $(\sigma_1\sigma_2)^3 = \Delta^2$, generator
   of the center of $B_3$ — gives the identity in both the universal family
   and its pullback to the invariant plane, consistent with
   $(3\text{-cycle})^3 = \mathrm{id}$.

Caveat on the *tables above*: labelled permutations are numerical
(nearest-neighbour matching at 30 digits).  The *group* they generate
and the Coxeter image of the cusp are now **Exact**
(`scripts/certified_monodromy.py`, resolving OPEN_QUESTIONS B6; see
`docs/MONODROMY.md`). The identities (I1)–(I4) and everything in §§1–2
and 4 remain exact symbolic algebra.

## 4. The local system, its decomposition, and what it explains

The sheet local system $\mathcal{L}$ on the wall complement (rank 3, fiber
= the three eliminant roots, monodromy = the $S_3$ permutation
representation composed with $B_3 \twoheadrightarrow S_3$) decomposes as

$$
\mathcal{L} \;=\; \underline{\mathbb{C}} \,\oplus\, \mathcal{S},
\qquad
\chi_{\mathrm{perm}} = \chi_{\mathrm{triv}} + \chi_{\mathrm{std}}
\quad (3 = 1 + 2,\; 1 = 1 + 0,\; 0 = 1 - 1),
$$

with $\mathcal{S}$ the rank-2 **standard local system** — the reflection
representation of $W(A_2)$, i.e. the vanishing cohomology of the $A_2$
singularity. This bookkeeping is exactly aligned with the two behaviours
already verified in this repository:

- **Trivial summand = trace observables.** The invariants of the fiber,
  $e_1 = 0$, $e_2 = (4-3w)/P_2$, $e_3 = 2/P_2$ and all power sums
  $S_k = \sum_i \xi_i^k$, are single-valued **rational** functions with
  poles only on the wall (verified for $k \le 6$; the invariant-plane
  version of `scripts/trace_pushforward.py`, with which it is consistent
  under $\xi = X/c$). The trivial summand is precisely the part of the
  sheet data that the CHY-type "sum over all solutions" mechanism sees.
- **Standard summand = everything else.** Any function of the sheets that
  is *not* symmetric (a single vacuum branch, a difference of branches,
  the perturbative sheet itself) transforms in $\mathcal{S}$ and is
  genuinely multivalued with the full $S_3$ monodromy. There is no room
  for anything in between: $\mathrm{triv} \oplus \mathrm{std}$ is the
  complete decomposition, so *trace rationality and $S_3$ multivaluedness
  are the only two behaviours the model admits*, and both are now
  identified representation-theoretically.

**Euler characteristic of the complement** (exact, script §4): the affine
wall curve is irreducible and rational; its normalization
$\mathbb{P}^1_m \to \overline{\{P_2=0\}}$ (the parametrization
$u(m), w(m)$ of `docs/POSITIVE_GEOMETRY.md`) is a *bijection* — the
resultant of the two parameter-collision equations is
$2187\,m_1^2(m_1-3)^2$, so identifications could only occur at $m = 0$ (the
unique place at infinity) and $m = 3$ (the cusp, where the only solution is
$m_2 = 3$: unibranch). Hence the affine curve is homeomorphic to
$\mathbb{P}^1$ minus one point $= \mathbb{C}$, and by additivity of the
(compactly-supported = topological, for complex algebraic varieties) Euler
characteristic:

$$
\chi\bigl(\mathbb{C}^2 \setminus \{P_2 = 0\}\bigr)
= \chi(\mathbb{C}^2) - \chi(\text{affine wall})
= 1 - 1 = 0 .
$$

Cross-check: $H^*(B_3;\mathbb{Z}) = (\mathbb{Z}, \mathbb{Z}, 0, \dots)$
[Arn69], so the $K(B_3,1)$ has $\chi = 1 - 1 = 0$. Consistent. (One
consequence recorded for future use: for *any* rank-$r$ local system on the
complement the twisted Euler characteristic is $r \cdot \chi = 0$, so
twisted $H^1$ and $H^2$ always have equal dimensions — the expected home of
a perfect intersection pairing.)

## 5. Interpretation for the amplitudes program — what is proven, what is proposed

**Proven here (exact unless noted):**

- The identities (I1)–(I3): the invariant eliminant, the tautological
  polynomial map to depressed cubics, the discriminant pullback
  $4\hat Q^3 + 27\hat R^2 = 4P_2^3(27u-9w+8)^2$, wall $\mapsto$ cusp.
- (I4): the explicit affine isomorphism of pairs
  $(\mathbb{C}^2, \text{wall}) \cong (\mathbb{C}^2, A_2\text{-discriminant})$,
  hence $\pi_1(\text{wall complement}) = B_3$ and the $K(\pi,1)$ property
  (by [Arn69, Bri71, Del72] — imported theorems, not re-proved here).
- The sheet monodromy factors through the canonical
  $B_3 \twoheadrightarrow S_3$; wall meridians $\mapsto$ transpositions,
  $D_0$-meridians $\mapsto$ id, cusp loop $\mapsto$ 3-cycle (Coxeter
  element), full twist $\mapsto$ id.  **Exact** as of B6
  (`scripts/certified_monodromy.py`; local Puiseux + irreducibility +
  cusp leading model); the labelled-permutation tables in §3 remain a
  Numerical cross-check.
- The decomposition $\mathcal{L} = \underline{\mathbb{C}} \oplus
  \mathcal{S}$, trace rationality as the trivial summand, and
  $\chi(\text{complement}) = 0$.

**Proposal (not proven, offered as the B2 answer):** the correct
"amplituhedron-analogue" object for this model is the pair

$$
\Bigl(\;\mathbb{C}^2 \setminus \{P_2 = 0\}\;,\;\; \mathcal{S}\;\Bigr)
\;\;=\;\; \text{(the } K(B_3,1)\text{, its reflection local system)},
$$

*not* the real chamber. The grounds: the real chamber provably carries no
canonical form (`docs/POSITIVE_GEOMETRY.md` — vertex collision at the
cusp), whereas the complement-with-local-system is exactly the structure
that survived every computation: log forms $dP_2/P_2$ live on it, the
$S_3$-monodromy acts on it, trace observables are its invariant sector, and
its topology is *universal* — the classifying space of $B_3$, independent
of the particular counterexample realization. In this reading the analogue
of a canonical-form integral is a **twisted period**: a pairing between
twisted cycles and $\mathcal{S}$-valued (or Kummer-twisted
$P_2^{s}$-weighted) forms on the complement, in the established framework
of twisted (co)homology and intersection numbers for hypergeometric-type
arrangements [AK11; cf. the intersection-number formulation of amplitudes,
MM19-style]. The concrete next computation this proposal pins down:
$\dim H^1(\mathbb{C}^2 \setminus \{P_2=0\}, \mathcal{S})$ and its
intersection pairing (finite-dimensional linear algebra over the cuspidal
cubic complement — genuinely computable with the tools in this repo,
cf. `docs/OPEN_QUESTIONS.md` B2). What is *not* claimed: any statement
about physical amplitudes, kinematic spaces, or the amplituhedron proper —
the connection remains at the level of mechanisms
(`docs/AMPLITUDES_CONNECTION.md` §2.5), now with the correct geometric
carrier identified on the toy-model side.

## 6. New questions this opens

1. **Universality over Keller maps.** For any counterexample with fiber
   degree $d$, elimination gives a degree-$d$ eliminant; after the
   C\*-reduction (when a torus action exists), does the wall always map to
   the space of degree-$d$ polynomials pulling back the $A_{d-1}$
   discriminant — and is the wall pair always *affinely equivalent* to (a
   linear section of) the $A_{d-1}$ discriminant, making the wall
   complement a $K(\pi,1)$ for (a subgroup of) $B_d$? Here this held in the
   strongest possible form (equality, multiplicity one). A degree-4
   example (`docs/OPEN_QUESTIONS.md` C4) would be the first real test:
   for $d \geq 4$ the discriminant hypersurface of
   $x^d + Q_2 x^{d-2} + \dots + Q_d$ lives in $\mathbb{C}^{d-1}$, so a 2D
   invariant plane could only meet it in a *section*, and generic sections
   of $K(\pi,1)$ discriminants are $K(\pi,1)$ only in favorable cases
   (Zariski–Lefschetz gives $\pi_1$ surjectivity; asphericity is extra).
2. **Is the braid-group structure forced by non-properness?** The
   $A_2$-equivalence used only: cuspidal cubic + flex tangent at infinity.
   Our rigidity evidence (`docs/NEW_COUNTEREXAMPLES.md`) constrains the
   weight system; does the Keller condition itself force the cusp (and
   hence $B_3$), or could a nodal-wall counterexample exist (whose
   complement would have a different, non-braid $\pi_1$)? Same fork as
   `docs/OPEN_QUESTIONS.md` C2, now with a topological reformulation.
3. **Compute the twisted cohomology.** $H^*(\text{complement},
   \mathcal{S})$ and the $P_2^s$-twisted versions, with intersection
   pairing; compare against the residueless double pole that killed the
   chamber form — does the obstruction reappear as a degenerate
   intersection matrix at $s \in \mathbb{Z}$, and is it invertible for
   generic $s$ (the usual genericity mechanism of twisted period
   integrals)?
   *(✅ Resolved 2026-07-21: `scripts/twisted_cohomology.py`,
   `docs/TWISTED_PERIODS.md`. $H^*(M;\mathcal{S}) = (0,1,1)$; the Kummer
   twist jumps exactly at $s \in \mathbb{Z}$ and $s \in \pm\frac16 +
   \mathbb{Z}$ (trefoil Alexander polynomial / cusp spectrum); at integer
   $s$ the canonical form is globally twisted-exact — the cohomological
   face of the residueless double pole; generic twists have NO classes at
   all ($|\chi| = 0$), so everything interesting is resonant.)*
4. **Lift to the 3D cone.** The identification lives in the invariant
   plane. The 3D wall is the C\*-cone over it; is
   $\mathbb{C}^3 \setminus \{p = 0\}$ a $K(\pi,1)$ as well (e.g. a
   $\mathbb{C}^*$-bundle-like extension of $B_3$ up to the codimension-2
   locus $\{c = 0\}$), and what group plays the role of $B_3$ there? The
   chart decomposition $\{c \neq 0\} \cong \mathbb{C}^* \times
   \mathbb{C}^2_{(u,w)}$ suggests $\mathbb{Z} \times B_3$ on the chart, but
   gluing across $\{c = 0\}$ needs care.
5. **Operator meaning of the Coxeter element.** The cusp loop — around the
   *empty-fiber* orbit, where the classical theory has no solution at all —
   acts as a Coxeter element of order $h = 3$ cyclically rotating the three
   vacua. Is there a 0D-QFT statement ("total escape = Coxeter rotation of
   the vacuum set") that survives in the $S(f)$-caricature program
   (`docs/OPEN_QUESTIONS.md` B1)?

## 7. References

- [Arn69] V. I. Arnold, *The cohomology ring of the colored braid group*,
  Mat. Zametki **5** (1969) 227–231; English transl. Math. Notes **5**
  (1969) 138–140. (Cohomology of braid/pure-braid groups;
  $H^*(B_3) = (\mathbb{Z},\mathbb{Z},0,\dots)$.)
- [Bri71] E. Brieskorn, *Die Fundamentalgruppe des Raumes der regulären
  Orbits einer endlichen komplexen Spiegelungsgruppe*, Invent. Math. **12**
  (1971) 57–61. ($\pi_1$ of regular-orbit spaces of finite complex
  reflection groups = generalized braid groups; for $W(A_2)$: $B_3$.)
- [Del72] P. Deligne, *Les immeubles des groupes de tresses généralisés*,
  Invent. Math. **17** (1972) 273–302. (The $K(\pi,1)$ property for
  complements of simplicial arrangements, covering all finite Coxeter
  types.)
- [Zar29] O. Zariski, *On the problem of existence of algebraic functions
  of two variables possessing a given branch curve*, Amer. J. Math. **51**
  (1929) 305–328. ($\pi_1(\mathbb{P}^2 \setminus \text{cuspidal cubic}) =
  \mathbb{Z}/2 * \mathbb{Z}/3$ — the projective statement, a proper
  quotient of the affine $B_3$.)
- [AK11] K. Aomoto, M. Kita, *Theory of Hypergeometric Functions*,
  Springer Monographs in Mathematics, 2011. (Twisted (co)homology and
  intersection theory for arrangement complements — the framework proposed
  in §5.)
- [MM19] P. Mastrolia, S. Mizera, *Feynman integrals and intersection
  theory*, JHEP **02** (2019) 139. arXiv:1810.03818. (Amplitudes as
  intersection numbers of twisted periods — the amplitudes-side face of
  the same framework.)
- Repo-internal: `docs/POSITIVE_GEOMETRY.md` (C\*-reduction, cusp
  geometry, failure of the chamber form), `docs/MONODROMY.md` ($S_3$
  monodromy in 3D), `docs/AMPLITUDES_CONNECTION.md` (mechanism-level
  dictionary), `scripts/trace_pushforward.py` (trace rationality),
  `scripts/wall_braid.py` (all verifications for this document).
