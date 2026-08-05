# Monodromy of the Alpöge–Mathew covering

**Result (Exact): the geometric monodromy group of the degree-3 covering
is the full symmetric group $S_3$.**  Wall meridians act as transpositions,
$D_0$-meridians act as the identity, and (via the affine $A_2$
identification of `docs/WALL_COMPLEMENT.md`) the monodromy representation
is the canonical surjection $B_3 \twoheadrightarrow S_3$, with the cusp
loop mapping to a Coxeter element (3-cycle).

Certification script: [`scripts/certified_monodromy.py`](../scripts/certified_monodromy.py)
(~1 s, 38 `[ok]` asserts).  The older high-precision tracker
[`scripts/monodromy.py`](../scripts/monodromy.py) remains as a Numerical
cross-check of labelled permutations on a Lefschetz line:

```bash
.venv/bin/python scripts/certified_monodromy.py   # Exact proof, ~1 s
.venv/bin/python scripts/monodromy.py             # Numerical tracker, ~80 s
.venv/bin/python scripts/monodromy.py 7           # independent cross-check line
```

## Setup

The map `F : C^3 -> C^3` (see `jcqft/core.py`) has `det DF = -2`
identically, so `F` is a local biholomorphism everywhere and fiber points can
never merge at finite points.  The generic fiber has 3 points.  Every
preimage of a target `(a, b, c)` has x-coordinate solving the eliminant cubic

```
p*X^3 + q*X + r = 0,   p = 27a²c² − 18abc + 16a + b³c − b²,
                       q = 4 − 3bc,   r = −2c,
```

and, given a root `x`, the remaining coordinates are rational functions
`y = −B/A`, `z = −D/C` from the lex Gröbner basis, with
`A = C/4 = 2·(27ac² − 9bc + 8)` independent of `x`.

The discriminant of the cubic factors as

```
disc_X = −p · (4q³ + 27pr²)  =  −4 · (27ac² − 9bc + 8)² · p ,
```

i.e. the x-collision factor `4q³ + 27pr²` is exactly `4·D0²` with
`D0 = 27ac² − 9bc + 8` — the *square* of the parametrization denominator.
So the singular locus of the covering restricted to a line has two kinds of
points:

* **`p = 0` — true branch locus.**  The x-eliminant drops degree; sheets
  escape to infinity in the x-direction only (the y- and z-eliminants are
  monic cubics, so y, z stay bounded).  Since `x ~ ±sqrt(−q/p)` for the
  escaping roots, a loop around a simple zero of `p` swaps the two
  escaping sheets.
* **`D0 = 0` (`p ≠ 0`) — x-collision locus.**  Two *distinct* fiber points
  share the same x-coordinate; the projection to x ramifies but the covering
  does not (`det DF = −2` forbids merging).  Trivial monodromy.

## Exact certification (B6 resolved)

Script: `scripts/certified_monodromy.py`.  The argument is algebraic; no
interval tracking is required.

1. **Galois group $S_3$ (Exact).**  The eliminant is irreducible over
   $K=\mathbb{Q}(a,b,c)$, and $\operatorname{disc}_X=-4\,D_0^2\,p$ with $p$
   irreducible of multiplicity one, hence not a square in $K$.  Classical
   cubic criterion: $\mathrm{Gal}(\text{Galois closure}/K)=S_3$.
2. **Local wall monodromy is a transposition (Exact).**  At the smooth
   wall point $(a,b,c)=(0,1,1)$ ($p=0$, $q=1$, $D_0=-1$, $\nabla p\neq0$),
   the transverse line $J(t)=(t,1,1)$ has $p(t)=t(27t-2)$.  The finite
   sheet $X_0=-r/q=2$ continues holomorphically by the implicit-function
   theorem ($\partial(\text{cubic})/\partial X|_{t=0}=q\neq0$).  The escaping
   sheets admit the Puiseux ansatz $X=Z/\sqrt{t}$; the cleared equation
   $(27s^2-2)Z^3+Z-2s=0$ ($s=\sqrt{t}$) is odd in $(Z,s)$, so a positive
   meridian $t\mapsto e^{2\pi i}t$ (i.e. $s\mapsto-s$) swaps the two
   escaping branches.  The same Puiseux equation governs the invariant
   eliminant on the cut $(u,w)=(t,1)$.
3. **Geometric monodromy $=S_3$ (Exact).**  Irreducibility $\Rightarrow$
   $\mathrm{Mon}$ transitive in $S_3$; the only transitive subgroups are
   $A_3$ and $S_3$; a transposition is odd $\Rightarrow\mathrm{Mon}=S_3$.
   Combined with $\mathrm{Mon}\subseteq\mathrm{Gal}=S_3$, equality is
   forced both ways.  The default Lefschetz line of `monodromy.py`
   ($J_0=(-5/6,2/5,0)$, $v=(3/2,-7/8,-3/2)$) is verified square-free /
   transverse, so the same conclusion holds after restriction to a line.
4. **Canonical $B_3\twoheadrightarrow S_3$ (Exact).**  Re-verify the
   affine isomorphism (I4) of `docs/WALL_COMPLEMENT.md`:
   $\pi_1(\mathbb{C}^2\setminus\{P_2=0\})=B_3$.  Wall meridians map to
   transpositions, so the representation is (conjugate to) the canonical
   surjection.  At the cusp, the quadratic part of $P_2$ is the perfect
   square $3(3\delta u-\delta w)^2$ and $q=0$, so the leading model is
   $\xi\sim(2/P_2)^{1/3}$ with $P_2$ winding twice on a small cusp loop:
   monodromy is a 3-cycle (Coxeter element of $W(A_2)$).
5. **$D_0$-meridians $= \mathrm{id}$ (Exact).**  Forced by $\det DF=-2$;
   confirmed at the rational x-collision point $(1/27,1,1)$ (fiber: three
   distinct points in $\mathbb{C}^3$).

### What remains Numerical

Individual *labelled* permutations on a chosen homotopy basis of loops
(approach-segment conjugacy), as tabulated by `scripts/monodromy.py` and
`scripts/wall_braid.py` §3.  Those labels depend on path choices; the
generated group and the Coxeter image of the cusp do not, and are Exact
above.

## Numerical cross-check (historical tracker)

Method of `scripts/monodromy.py` (unchanged):

1. **Generic line.**  Restrict to `J(t) = J0 + t·v` with small random
   rational `J0, v` (seeded; conditions checked symbolically: `p`, `D0`
   nonzero at `t = 0`, `p(J(t))` of full degree 4 and square-free, `D0(J(t))`
   of degree 3 and square-free, no common roots).  Singular parameters `t`
   are the roots of these two univariate polynomials (mpmath `polyroots` on
   exact rational coefficients).
2. **Base fiber.**  At `t = 0`: roots of the cubic in `x` (mpmath
   `polyroots`), then `y, z` from the rational parametrization, then a few
   Newton steps on the full system `F(φ) = J0`.  Working precision
   `mp.dps = 30`.
3. **Continuation.**  All three sheets are transported simultaneously along
   each path by tangent-predictor (`DF·φ' = v`) / Newton-corrector steps on
   `F(φ) = J(t)`.  Newton is uniformly well-conditioned because
   `det DF = −2` everywhere.  Steps bisect adaptively whenever Newton fails,
   the corrector moves a sheet more than 20 % of the current inter-sheet
   separation, or two tracked sheets approach within `1e−8` (sheet-jump
   guard).  The residual `|F(φ) − J(t)|` is checked at every accepted step
   and stayed below `1e−22` in all runs (required: `1e−20`).
4. **Loops.**  For each singular `t*`: basepoint → circle start → full
   counterclockwise circle of radius `0.3 ×` (distance to the nearest other
   singularity) → back.  The approach segment is chosen to stay as far as
   possible from the other singularities.  The transported fiber is matched
   to the base fiber by nearest neighbour with a `1e−10` tolerance and a
   1000× separation-ratio test — matching was unambiguous in every loop
   (typical return distances `~1e−20`).
5. **Big loop.**  One circle `|t| = 2·max|t*|` enclosing *all*
   singularities, with the radial approach segment steered away from the
   real roots of `p(t)`.

### Numerical results

Default line `J0 = (−5/6, 2/5, 0)`, `v = (3/2, −7/8, −3/2)`: `p(J(t))` has 4
simple roots, `D0(J(t))` has 3.  Base fiber sheets (labelled 1–3 by real part
of x): `x ≈ −0.5445, 0, +0.5445`; minimal sheet separation 5.7.

| loop around | locus | permutation |
|---|---|---|
| `t* ≈ −0.3442` | `p = 0` | **(1 2)** |
| `t* ≈ +0.5559` | `p = 0` | **(1 3)** |
| `t* ≈ 0.5794 − 0.4201i` | `p = 0` | **(2 3)** |
| `t* ≈ 0.5794 + 0.4201i` | `p = 0` | **(2 3)** |
| `t* ≈ −0.2733` | `D0 = 0` | id |
| `t* ≈ 0.4792 − 0.3026i` | `D0 = 0` | id |
| big loop `|t| ≈ 1.431` | all | id |

Generated group: **order 6 = S3**.  Max tracking residual over all loops:
`1e−22`; total wall time ≈ 80 s.

An independent line (seed 7: `J0 = (1/3, 3, −7/9)`, `v = (−1, 9, 7/4)`) gives
the same picture: transpositions (1 3), (1 2), (2 3), (2 3) around the four
`p`-roots, id around the `D0`-roots, id for the big loop, group S3.

Observations (now consequences of the Exact certification):

* **Every simple branch point is a transposition** — local model
  `x ~ ±sqrt(−q/p)`.
* **The x-collision locus is harmless**: identity permutation, as forced
  by `det DF = −2`.
* **Big-loop permutation = id** on both test lines: the product of the four
  branch-point transpositions (in path order, suitably conjugated) cancels,
  i.e. the covering restricted to these lines is unramified over `t = ∞`.

## Note for paper authors

The geometric monodromy group may now be stated as a theorem (not a
numerical observation): $\mathrm{Mon}=S_3=\mathrm{Gal}$, via local Puiseux
at a smooth wall point plus irreducibility of the eliminant; equivalently,
the sheet local system realizes the canonical $B_3\twoheadrightarrow S_3$
on the wall complement ($K(B_3,1)$ by the affine $A_2$ isomorphism of
`docs/WALL_COMPLEMENT.md`).  Cite `scripts/certified_monodromy.py`.  The
tables of labelled permutations on specific loops remain Numerical
illustrations.
