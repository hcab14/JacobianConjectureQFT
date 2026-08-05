# The Family (1, −1, −m), m ≥ 3: Uniform Classification of the v-Linear Class

*(Machinery: `jcqft/reduction_w.py`, proved exact for the whole family
$(1,-1,-m)$ in `scripts/reduction_113.py`. The $m = 3$ case, whose gauge,
stratification and divisibility chains this document generalizes, is
`docs/SEARCH_113.md` / `scripts/search_113.py`. Every identity and every
emptiness claim below is proved by assertion in `scripts/search_11m.py`,
with $m$ **symbolic** wherever the statement is uniform and $k$ symbolic
($m = 2k{+}1$ / $m = 2k$) where the UFD bookkeeping branches on parity:*

```
.venv/bin/python scripts/search_11m.py           # ~2 min, the whole proof
.venv/bin/python scripts/search_11m.py --full    # + larger boxes and the
                                                 #   targeted m = 5 gap
                                                 #   queries, ~70 min

```

*Concrete spot-checks run at $m = 2$ — the Alpöge–Mathew member, which
**must and does survive** — and at $m = 3$ (reproducing `search_113`),
$4$, $5$, $7$.)*

## 1. The theorem

> **Theorem (v-linear class, all $m \ge 3$).** For every integer
> $m \ge 3$, every Keller map in the $v$-linear class of the weight
> system $(1,-1,-m)$ — i.e. every polynomial $F = (P/x^m,\ Q/x,\ xR)$
> with $P, Q, R \in \mathbb{C}[w, v]$ of degree $\le 1$ in $v = x^mz$,
> arbitrary degree in $w = xy$, and $\det DF$ a nonzero constant — is,
> up to the gauge of §3, the **tame automorphism**
>
> $$
> P = p_0 + v,\qquad Q = w + b_0P,\qquad R = 1,\qquad
> p_0 \in w^m\,\mathbb{C}[w],\ b_0 \in \mathbb{C}.
> $$
>
> Consequently the $v$-linear class of $(1,-1,-m)$ contains **no
> counterexample and no $m{:}1$ orbifold covering for any $m \ge 3$**,
> and the Alpöge–Mathew counterexample ($m = 2$) is the **unique**
> member of its equivariant family within the $v$-linear class.
>
> *Gap qualification*: exactly as at $m = 3$, the structural proof does
> not reach stratum D3 with non-squarefree $s$ or $t$ (§7); that corner
> requires $\deg p_1 \ge m{+}1$ or $\deg r_1 \ge m{+}1$ ($m$ odd), resp.
> $\ge 2(m{+}1)$ ($m$ even), and is closed by exact in-box Gröbner
> certificates for $m \le 5$.

The theorem's content is the **contrast with $m = 2$**: the proof
machinery, run at $m = 2$, does *not* kill the Alpöge–Mathew map — it
pinpoints the two exact places where $m = 2$ escapes (§6).

## 2. Reduction and the $v$-linear system (all $m$ symbolic)

From `reduction_w` (asserted with $m$ a sympy `Symbol` on undetermined
functions — identities of differential polynomials, not spot checks):

$$
\det DF = \det M = -m\,P\,J_2(Q,R) + Q\,J_2(P,R) + R\,J_2(P,Q),\qquad
R^m \det M = J_2(PR^m,\,QR),
$$

a function of $(w,v)$ alone; note $\det M$ is **linear in $m$**, which is
what lets the whole derivation run with symbolic $m$. Polynomiality:
monomials $w^jv^k$ of $P$ need $j + mk \ge m$, of $Q$ need $j + mk \ge 1$
($R$ unconstrained). $DF(0) = \mathrm{antidiag}(p_1(0), q_0'(0),
r_0(0))$; target scalings fix the **gauge** $p_1(0) = q_0'(0) = r_0(0) =
1$, hence $\kappa = \det DF = -1$. The $v$-shear group $v \to v + h(w)$
leaves $\det M$ invariant (asserted, $m$ symbolic) and is the main
analysis tool; the polynomiality box is imposed *last*, through jets at
$w = 0$.

With $P = p_0 + p_1v$, $Q = q_0 + q_1v$, $R = r_0 + r_1v$, the Keller
condition is $E_2 = E_1 = 0$, $E_0 = \kappa$, with (asserted, $m$
symbolic)

$$
E_2 = 2p_1'q_1r_1 - (m{+}1)p_1q_1'r_1 + (m{-}1)p_1q_1r_1',\qquad
E_0 = -m p_0(q_0'r_1 - q_1r_0') + q_0(p_0'r_1 - p_1r_0') + r_0(p_0'q_1 - p_1q_0'),
$$

and $E_1$ the mixed analogue.

## 3. The general integrated constraint

$E_2$ is a pure log-derivative relation: dividing by $p_1q_1r_1$,

$$
2\,\frac{p_1'}{p_1} \;-\; (m{+}1)\,\frac{q_1'}{q_1} \;+\; (m{-}1)\,\frac{r_1'}{r_1} \;=\; 0 ,
$$

so the exponent vector is $(\alpha,\beta,\gamma) = (2,\,-(m{+}1),\,m{-}1)$.
The asserted Wronskian certificate, $m$ symbolic,

$$
\Bigl(\frac{p_1^2\,r_1^{m-1}}{q_1^{m+1}}\Bigr)'
= \frac{p_1\,r_1^{m-2}}{q_1^{m+2}}\,E_2
$$

integrates it in the UFD $\mathbb{C}[w]$:

$$
\boxed{\ p_1^2\,r_1^{m-1} \;=\; c\,q_1^{m+1},\qquad c \neq 0\ }
\qquad\text{whenever } q_1p_1r_1 \neq 0 .
$$

**Anchors** (both asserted): $m = 2$ gives $p_1^2r_1 = cq_1^3$, satisfied
by Alpöge–Mathew with $c = -1/27$ ($p_1 = (1{+}w)^3$, $q_1 = 3(1{+}w)^2$,
$r_1 = -1$); $m = 3$ gives $(p_1r_1)^2 = c(q_1^2)^2$, the square of
`search_113`'s $p_1r_1 = cq_1^2$.

## 4. Stratification and the parity-dependent UFD parametrization

Since $p_1(0) = 1$ forces $p_1 \ne 0$, the strata are A ($q_1 = r_1 =
0$), B ($r_1 = 0, q_1 \ne 0$), C ($q_1 = 0, r_1 \ne 0$), D ($q_1r_1 \ne
0$). In D: $E_0$ is linear in $(p_1, q_1, r_1)$ with no derivatives of
them (asserted, $m$ symbolic), so $\gcd(p_1,q_1,r_1) \mid E_0 = \kappa$
is a constant — absorb it. Then per irreducible $\pi$ with multiplicities
$(A,B,C)$ in $(p_1,q_1,r_1)$: $2A + (m{-}1)C = (m{+}1)B$ and
$\min(A,B,C) = 0$; $B = 0$ forces $A = C = 0$, so $\mathrm{rad}q_1
= \mathrm{rad}(p_1r_1)$ and $\gcd(p_1,r_1) = 1$. The minimal
solutions of $2A = (m{+}1)B$ resp. $(m{-}1)C = (m{+}1)B$ depend on the
parity of $m$ (exponent arithmetic asserted with $k$ symbolic):

$$
m = 2k{+}1:\quad p_1 = a\,s^{k+1},\quad q_1 = s\,t^{k},\quad r_1 = b\,t^{k+1};
\qquad
m = 2k:\quad p_1 = a\,s^{2k+1},\quad q_1 = s^2t^{2k-1},\quad r_1 = b\,t^{2k+1},
$$

with $\gcd(s,t) = 1$ (one normalization constant absorbed into $q_1$).
At $m = 3$ this is exactly `search_113`'s $(as^2, st, bt^2)$; at $m = 2$
it is $(as^3, s^2t, bt^3)$, housing Alpöge–Mathew in the $t = $ const
slot. Substrata: **D0** ($s,t$ const), **D1** ($t$ const, $s$ nonconst —
*the Alpöge–Mathew slot*), **D2** ($s$ const), **D3** (both nonconst).

| stratum | verdict for $m \ge 3$ | uniform in $m$? |
|---|---|---|
| A, B | the tame family of §1 (all $w$-degrees) | yes ($m$ symbolic) |
| C | **empty** | yes ($m$ symbolic) |
| D0 | **empty** | yes ($m$ symbolic) |
| D1 | **empty** | $m$ ≥ 4: outright, $k$ symbolic per parity; $m = 3$: box jet |
| D2 | **empty** | yes ($k$ symbolic per parity; one case uses the box) |
| D3, $s,t$ squarefree | **empty** | $m$ ≥ 4: outright; $m = 3$: box |
| D3′, $s$ or $t$ non-squarefree | **empty in boxes** ($m \le 5$), open beyond | the gap, §7 |

## 5. The classification chains (every identity asserted)

**A and B — the tame family, $m$ symbolic.** A: $E_1 \equiv 0$, $E_0 =
-p_1(q_0r_0)' = \kappa \Rightarrow p_1 \mid \kappa \Rightarrow p_1 = 1$,
$q_0r_0 = w \Rightarrow r_0 = 1$, $q_0 = w$. B: the certificate
$(p_1r_0^{m-1}/q_1)' = r_0^{m-2}E_1/q_1^2$ integrates to $p_1r_0^{m-1} =
cq_1$; then $E_0 = p_1\,[\,(p_0r_0^m)/c - q_0r_0\,]'$, forcing $p_1 = 1$,
$r_0 \mid w$ with $r_0(0) = 1 \Rightarrow r_0 = 1$, $q_0 = b_0p_0 + w$,
$q_1 = b_0$. The family assembles to $F = (p_0(xy)/x^m + z,\ y +
b_0x^{m-1}F_1,\ x)$: an elementary shear composed with a target shear
$b \mapsto b + b_0ac^{m-1}$ — **tame**, explicit inverse verified by
composition, generic fiber 1 point, prefilter survived only through its
known false-positive class (all asserted at $m = 3, 4, 5, 7$).

**C, $m$ symbolic.** $(p_1/(q_0^{m+1}r_1))' = E_1/(q_0^{m+2}r_1^2)$, so
$E_1 = 0$ forces $p_1 = cq_0^{m+1}r_1$, whence $p_1(0) = 0$ ($q_0(0)=0$):
contradiction. Empty for every $m \ge 2$.

**D0, $m$ symbolic.** $E_1 = [2\beta\gamma p_0 - (m{+}1)\alpha\gamma q_0 +
(m{-}1)\alpha\beta r_0]'$; solving and substituting into $E_0$ gives the
exact derivative $[-\tfrac{m+1}{m-1}\tfrac{\alpha\gamma}{\beta}V^2 -
\alpha\delta V]' = \kappa$ with $V = q_0 - (\beta/\alpha)p_0$, $V(0)=0$,
$V'(0)=1$. Then $V \mid \kappa w \Rightarrow V = w$, and the $w^2$
coefficient forces $\tfrac{m+1}{m-1}\tfrac{\alpha\gamma}{\beta} = 0$ —
impossible for every $m \ge 2$.

**D1 — the Alpöge–Mathew slot ($t$ const, $u := s$ nonconst, $r_1 = C$).**
The shear $h = -r_0/C$ sets $r_0 = 0$. $E_1$ is a linear ODE for $p_0$
with particular solution $p_0^* = \tfrac{m+1}{2}\tfrac{A}{B}u^{e-f}q_0$
and homogeneous solutions $Y^2 = e\,q_1^m$ (certificates asserted, $k$
symbolic, both parities):

- *$m$ odd* ($p_1 = Au^{k+1}$, $q_1 = Bu$): $Y \neq 0$ needs $u = cd^2$.
  With $Y = 0$: $E_0 = \tfrac{m+1}{2}k\,C\tfrac{A}{B}\,
  \mathbf{u^{(m-3)/2}}\,q_0(u'q_0 - 2uq_0')$; for $m \ge 5$ the bold
  factor forces $u \mid \kappa$: **empty outright**. At $m = 3$ the factor
  is trivial; then $q_0 = n$ const, $u$ linear, and the un-shear jet
  system $p_0(0) = q_0(0) = 0$ forces $n = 0$, $\kappa = 0$: **empty in
  the box** (`search_113`'s argument, re-asserted). With $u = cd^2$:
  $E_0 = C\,\mathbf{d^{m-2}}(d'q_0 - dq_0')(2k(k{+}1)\tfrac{A}{B}c^kq_0 +
  m y_0 d)$, and $d \mid \kappa$ kills it for every odd $m \ge 3$.
- *$m$ even* ($p_1 = Au^{2k+1}$, $q_1 = Bu^2$): the homogeneous direction
  $Y = \tilde{y}\,u^m$ is **always** polynomial, and the complete
  $\kappa$-equation reads
  $$
  E_0 \;=\; C\,\mathbf{u^{m-2}}\,(u'q_0 - uq_0')
  \Bigl[\tfrac{(m-1)(m+1)}{2}\tfrac{A}{B}\,q_0 + m\,\tilde{y}\,u\Bigr]
  \;=\; \kappa .
  $$
  For every even $m \ge 4$ the factor $u^{m-2}$ forces $u \mid \kappa$:
  **empty outright**. At $m = 2$ the factor is trivial — see §6.

**D2 ($s$ const, $p_1 = A$).** The shear $h = -p_0/A$ sets $p_0 = 0$;
then $E_0 = -A(q_0r_0)'$ for *arbitrary* $q_1, r_1$ and symbolic $m$
(asserted), so $q_0r_0$ is exactly linear: one factor a nonzero constant,
the other exactly linear. Case $r_0 = \delta$ const: the $E_1$ closed
form gives $t \mid \delta t'$ ($t^2 \mid$ for $m$ even) $\Rightarrow
\delta = 0$ or $t$ const: empty. Case $q_0 = \gamma$ const: for $m$ odd
the complete $E_1$-solution is $r_0 = \rho_0t + Y$ with $\rho_0 =
\tfrac{(k+1)\gamma C}{kB} \neq 0$ fixed and $Y^2 = et$ (so $\deg Y <
\deg t$; $Y \neq 0$ only for $t = ce^2$, and $r_0 = \rho_0t + y_0e$ is
asserted to solve $E_1$ there); hence $\deg r_0 = \deg t$, and $r_0$
exactly linear forces $t$ **linear** (then $Y = 0$). Un-shearing needs
$\mathrm{val}(p_0 = -Ah) \ge m \Rightarrow h(0) = 0 \Rightarrow
q_0(0) = \gamma \ne 0$, violating the box: empty (uses only
$\mathrm{val}p_0 \ge 1$, so uniform). For $m$
even the complete solution $r_0 = \sigma t^2 + \rho t$ with $\sigma =
\tfrac{(m+1)\gamma C}{(m-1)B} \neq 0$ has degree $2\deg t \ge 2$: never
exactly linear, empty without the box.

**D3 ($s, t$ nonconstant, coprime, squarefree).** Shear-invariants
($\Lambda_1 = q_1p_0 - p_1q_0$ etc., stripped of their forced $s,t$
factors):

$$
m\ \text{odd}:\ G_1 = t^kp_0 - as^kq_0,\quad G_2 = btq_0 - sr_0;
\qquad
m\ \text{even}:\ G_1 = t^{m-1}p_0 - as^{m-1}q_0,\quad G_2 = bt^2q_0 - s^2r_0 .
$$

Asserted, $k$ symbolic: $E_1 = aks^kt^{k-1}\Theta + bt^k\Theta_1$ ($m$
odd) resp. $E_1 = (m{-}1)as^{2k}t^{2k-2}\tilde\Theta + 2bst^{2k}C_0$ ($m$
even), with $\Theta$-type combinations of $G_2$ and $\Theta_1/C_0$-type
combinations of $G_1$. From $E_1 = 0$:

1. *$t$-side (one step)*: $t \mid t'G_2 \Rightarrow t \mid G_2$
   ($t$ squarefree).
2. *$s$-side (iteration)*: the operator identity $B_j(sh) = sB_{j+1}(h)$
   with $B_j(h) \equiv -(m{-}2j)s'th \pmod{s}$ ($m$ odd; $(m{-}j)$ for
   $m$ even) forces $G_1 = s^kg$ resp. $G_1 = s^{m-1}g$ — the iteration
   coefficients $m{-}2j$ (odd) never vanish, and $m, m{-}1, \dots, 2$
   (even) never vanish, **for every $m \ge 3$**. (At $m = 2$ the even
   chain also completes; the divergence is *not* here.)
3. *Integration*: $m$ odd: $E_1 = s^kt^k(WZ - 2stZ')$ with $W = s't -
   st'$, $Z = akg_2 - bg_k$, and $(Z^2t/s)' \propto Z(2stZ' - WZ)$ gives
   $Z^2t = es \Rightarrow Z = 0$. $m$ even: $E_1 = s^mt^{2k}(s'Z - sZ')$
   with $Z = (m{-}1)ag_2 - 2btg$, so $(Z/s)' = 0$, $Z = \varepsilon s$.
4. *The $\kappa$-equation* (asserted, $k$ symbolic; $q_0$ drops out — it
   is the shear direction):
   $$
   m\ \text{odd}:\ E_0 = -\tfrac{ak(k+1)}{b}\,\mathbf{s^{k-1}}\,
   g\,(2stg' - Wg);\qquad
   m\ \text{even}:\ E_0 = -\tfrac{\mathbf{s^{m-2}}\,
   \bigl((2m{+}2)btg + \varepsilon s\bigr)(stg' - Wg)}{(m{-}1)a} .
   $$
   For $m \ge 5$ odd the factor $s^{k-1}$, and for $m \ge 4$ even the
   factor $s^{m-2}$, force $s \mid \kappa$: **empty outright, no box
   needed**. At $m = 3$ ($k = 1$) the equation reduces to $g(2stg' - Wg)
   = $ const, forcing $g = c$, $\kappa = (2a/b)c^2W$, $W$ const — a
   single shear-orbit killed by the box: $tp_0 = as(q_0 + c/b)$ gives
   $p_0(0)t(0) \ne 0$ while the box demands $p_0(0) = 0$ (re-asserted
   from `search_113`). At $m = 2$ the factor $s^{m-2}$ is trivial — the
   same escape hatch as D1.

## 6. The exact $m = 2$ divergence point

The Alpöge–Mathew map lives in stratum **D1 with $m$ even**: $u = 1{+}w$,
$p_1 = u^3$, $q_1 = 3u^2$, $r_1 = -1$. All of the following is asserted
end-to-end in the script:

- sheared by $h = -r_0/C = 2 - 3w$: $p_0^{\rm sh} = (w{+}1)(w{+}2)$,
  $q_0^{\rm sh} = 6 + 4w$;
- $p_0^{\rm sh} - p_0^* = -(1{+}w)^2 = \tilde{y}u^m$ with $\tilde{y} =
  -1$: the AM map **uses the even-$m$ homogeneous direction**;
- the $\kappa$-equation of §5-D1 evaluates to $E_0 = C\cdot u^{0}\cdot
  (u'q_0 - uq_0')\,[\tfrac{3}{2}\tfrac{A}{B}q_0 + 2\tilde{y}u] =
  (-1)\cdot 2\cdot(3 + 2w - 2 - 2w) = -2 = \kappa_{\rm AM}$ exactly.

So $m = 2$ escapes through **two** simultaneous degeneracies of the same
formula that kills every $m \ge 4$:

1. the killing factor $u^{m-2}$ (resp. $s^{m-2}$, $s^{k-1}$, $d^{m-2}$,
   $u^{(m-3)/2}$ in the sibling substrata) is **trivial precisely at
   $m = 2$** (and, in its weaker odd form, at $m = 3$);
2. the homogeneous $E_1$-direction $Y = \tilde{y}u^m$ exists only for
   $m$ **even**, and AM needs $\tilde{y} = -1 \neq 0$ (with $\tilde y =
   0$ the remaining equation forces $q_0$ const and the box jet kills it
   even at $m = 2$).

At $m = 3$ neither killing factor is active either, and the residue is
disposed of by the **box valuation** $\mathrm{val}_wp_0 \ge m$ (the
jet systems of §5) — this is `search_113`'s "numerological obstruction",
now exhibited as the boundary case of the uniform mechanism: for
$m \ge 4$ the polynomial factor kills D1/D3 before the box is even
consulted, at $m = 3$ the box does it, at $m = 2$ nothing does — and
Alpöge–Mathew exists.

## 7. The gap: D3 with non-squarefree $s$ or $t$

Steps 1–2 of §5-D3 use squarefreeness of $t$ (once) and of $s$ (at every
iteration step). A non-squarefree $s$ or $t$ needs $\deg s \ge 2$ or
$\deg t \ge 2$, i.e.

$$
m\ \text{odd}:\ \deg p_1 \ge m{+}1\ \text{or}\ \deg r_1 \ge m{+}1;\qquad
m\ \text{even}:\ \deg p_1 \ge 2(m{+}1)\ \text{or}\ \deg r_1 \ge 2(m{+}1),
$$

exactly the D3′ gap of `docs/SEARCH_113.md` §5.3 (threshold $4$ at
$m = 3$), with thresholds that **grow with $m$**. In-box closure
(msolve, exact over $\mathbb{Q}$, 16 GB cap):

1. **Default run**: Keller $+$ ($r_1$ has a nonzero coefficient) is the
   unit ideal in the boxes $\deg(p_0,p_1,q_0,q_1,r_0,r_1) \le
   (m{+}1,2,3,2,2,2)$ for $m = 4, 5$ — covering strata C and D wholesale
   in-box; and Rabinowitsch certificates on the $r_1 = 0$ slice pin the
   in-box variety exactly to the gauged family of §1 (with its box
   truncation $b_0p_j = 0$).
2. **`--full`**: the same in medium boxes $\le (m{+}2,3,4,3,3,3)$
   (20 unknowns, 32–34 equations), plus the two targeted $m = 5$
   non-squarefree parametrizations $s = (1{+}\rho w)^2$ /
   $t = (e_0{+}e_1w)^2$ with the partner linear and
   $\deg(p_0,q_0,r_0) \le (7,5,4)$. Status (2026-07-24, this machine,
   30 GB, total run 72 min): medium boxes $m = 4$ **EMPTY** (four
   queries, 6–13 min each) and $m = 5$ **EMPTY** (four queries, 2–10
   min); targeted $m = 5$ gap, $s$-square/$t$-linear: **EMPTY** (exact,
   150 s); targeted $t$-square/$s$-linear: msolve exceeds the 16 GB F4
   cap and the Singular degBound ladder (4–7) finds no unit certificate
   — formally **unresolved**, the same wall as the $m = 3$ analogues
   (`docs/SEARCH_113.md` §8). At $m = 4$ the gap thresholds
   ($\deg p_1 \ge 10$ or $\deg r_1 \ge 10$) are beyond any feasible
   targeted box, but for the same reason the gap only opens at
   $w$-degrees $\ge 10$ there.
3. Note the small boxes ($\deg p_1, \deg r_1 \le 2 < m{+}1$) contain
   **no** non-squarefree configuration at all for $m \ge 4$ — in-box the
   classification is complete *without* any gap; the gap is a statement
   about degrees $\ge m{+}1$ only.

## 8. The $m{:}1$ orbifold mechanism is empty; no other exists

Target coordinates carry weights $(-m,-1,1)$. The stabilizer of a point
on the $a$-axis is $\mathbb{Z}_m$; every other coordinate axis and every
pair has trivial stabilizer ($\gcd(m,1) = 1$), and for **composite $m$**
there is no intermediate $\mathbb{Z}_d$ stratum — coordinate stabilizers
are $\mathbb{Z}_m$ or trivial, so any $d{:}1$ sub-escape would still need
a witness on the $a$-axis: $Q = R = 0 \neq P$ on a free orbit. By the
classification every $v$-linear Keller map has $R = 1$ (gauge): **$R$
never vanishes, no witness exists, for any $m \ge 3$** (up to the
box-closed D3′ gap). Independent pointwise Gröbner queries at the
torus-normalized witness points $(1,1)$, $(1,0)$, $(0,1)$ confirm
emptiness in-box at $m = 4, 5$.

## 9. Honest limitations

1. **$v$-degree.** The class is $v$-linear (the AM-analogue class).
   At $m = 2$ the $v$-quadratic class is closed in-box — only tame + AM
   (`docs/SEARCH_VQUAD.md`); $v$-degree $\ge 3$ and $v$-quadratic at
   $m \ge 3$ remain open.
2. **The D3′ gap** (§7): structural proof reaches only squarefree $s,t$;
   beyond the boxes the non-squarefree corner is open — though for
   $m \ge 4$ the $\kappa$-equation kill needs *no box at all* on the
   squarefree part, and the gap thresholds grow linearly in $m$. Of the
   two targeted $m = 5$ queries one is closed exactly, the other is
   msolve-memory-unresolved, exactly like its $m = 3$ siblings.
3. **Gauge bookkeeping.** Everything is stated in the gauge $p_1(0) =
   q_0'(0) = r_0(0) = 1$, $\kappa = -1$; un-gauged solutions are
   recovered by target/torus scalings (tameness, properness, fiber
   counts unchanged). The $v$-shear moves are not box-preserving; the
   box enters only through jets at $w = 0$ — and is *needed* only at
   $m = 3$ (D1, D2-odd, D3) and in D2-odd for general $m$.
4. **Parity bookkeeping is exact but branch-wise.** All identities are
   asserted with $k$ symbolic within each parity branch ($m = 2k{+}1$,
   $m = 2k$); no statement is extrapolated across the branch point.

## 10. QFT reading

The $(1,-1,-m)$ gradings are the complete family of free-invariant-ring
siblings of the Alpöge–Mathew theory, with the $\mathbb{Z}_2$ orbifold
stratum replaced by $\mathbb{Z}_m$ — the candidate $m$-th-root-escape
global anomalies. The verdict is a uniform **no-go for all $m \ge 3$**:
the unit-determinant constraint is compatible with the grading only
through shear-type field redefinitions; no non-properness defect, no
$\mathbb{Z}_m$ vacuum-escape monodromy, for any $m$. Combined with
`SEARCH_213.md` (the $A_1$-cone sibling) the pattern sharpens into a
statement about the *defect itself*: within this entire equivariant
landscape the Alpöge–Mathew anomaly is not one instance of a family but
an **isolated numerological coincidence of $m = 2$** — the unique weight
where the killing factor $u^{m-2}$ degenerates *and* the even-parity
homogeneous direction exists. The natural remaining probes are the
$v$-quadratic class at $m = 2$ (closed in-box: `SEARCH_VQUAD.md`),
higher $v$-degree, and 4-field gradings.
