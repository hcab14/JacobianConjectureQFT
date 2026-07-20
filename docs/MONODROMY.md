# Monodromy of the Alpöge–Mathew covering

**Result (numerical): the monodromy group of the degree-3 covering is the full
symmetric group S3.**  Every simple loop around the branch locus `{p = 0}`
acts as a transposition, loops around the x-collision locus
`{4q^3 + 27pr^2 = 0}` act trivially, and different branch points give
*different* transpositions, so together they generate all of S3.

Script: [`scripts/monodromy.py`](../scripts/monodromy.py).  Reproduce with

```bash
.venv/bin/python scripts/monodromy.py          # default line (seed 20260720), ~80 s
.venv/bin/python scripts/monodromy.py 7        # independent cross-check line
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
  escaping roots, a loop around a simple zero of `p` can swap the two
  escaping sheets.
* **`D0 = 0` (`p ≠ 0`) — x-collision locus.**  Two *distinct* fiber points
  share the same x-coordinate; the projection to x ramifies but the covering
  does not (`det DF = −2` forbids merging).  Trivial monodromy expected.

## Method

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

## Numerical results

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

Observations:

* **Every branch point is a simple transposition** — exactly the behaviour
  predicted by the local model `x ~ ±sqrt(−q/p)`: the two sheets that escape
  to infinity swap; the third (finite) sheet is fixed.  Note that *two*
  x-roots blow up as `p → 0` along a generic path, but generically only one
  fiber point is lost over `{p = 0}` itself (the two escaping x-branches
  merge at infinity like a square root).
* **The x-collision locus is confirmed harmless**: identity permutation,
  as forced by `det DF = −2` (fiber points there are distinct in C^3; only
  their x-coordinates collide).
* **Big-loop permutation = id** on both test lines: the product of the four
  branch-point transpositions (in path order, suitably conjugated) cancels,
  i.e. the covering restricted to these lines is unramified over `t = ∞`.
  Consistency check: an even product of four transpositions can be id, and
  the parity matches automatically.

## Symbolic cross-check

`disc_X = −4·D0²·p` is *not* a square in `C(a,b,c)` (the factor `p` occurs
with multiplicity 1), so the Galois group of the x-cubic over `Q(a,b,c)` is
S3.  The geometric monodromy group is a priori only a subgroup of it; the
numerics above show it actually *equals* S3.  In particular the covering is
connected of degree 3 with no intermediate subcover: no sheet is globally
distinguishable, even though over the perturbative vacuum `J = 0` (which lies
*on* `{p = 0}`) only the "perturbative" sheet is visible.

## Caveats

* This is a **numerical** computation: permutations are read off by
  nearest-neighbour matching (tolerance `1e−10`) of fibers tracked at 30
  significant digits with residuals `≤ 1e−22`.  The margin between residual
  and matching tolerance is ~12 orders of magnitude, and the sheet-jump
  guard bounds each corrector step by 20 % of the sheet separation, but it is
  not a certified (interval-arithmetic) proof.
* Monodromy is computed on **one complex line** (plus one cross-check line).
  For a generic line this captures the full monodromy group of the covering
  over the complement of the discriminant (Zariski/Lefschetz-type argument),
  and S3 is in any case maximal for a degree-3 covering; but the code does
  not verify the genericity of the line beyond the explicit square-freeness
  and transversality checks listed above.
* Individual permutation *labels* depend on the homotopy class of the chosen
  approach segments (changing them conjugates the permutation); the generated
  group does not.
