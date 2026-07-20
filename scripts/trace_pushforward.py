"""Amplitudes-program structures in the counterexample, verified exactly.

Two statements connecting the toy model to the S-matrix/positive-geometry
program (see docs/AMPLITUDES_CONNECTION.md):

1. TRACE MAP / PUSHFORWARD RATIONALITY.  Single-sheet observables are
   multivalued algebraic functions with S3 monodromy; but the pushforward
   (sum over all three sheets, i.e. the algebraic trace) of any polynomial
   observable is a RATIONAL function of the sources with poles only on the
   non-properness divisor {p = 0}.  This is the toy-model version of the
   CHY sum over all scattering-equation solutions and of the canonical-form
   pushforward of Arkani-Hamed-Bai-Lam.  Verified here for the power sums
   S_k = sum_i x_i^k, k = 2..6, via Newton's identities.

2. BOUNDARY FACTORIZATION AT THE WALL.  Near {p = 0} the fiber splits into
   an escaping pair x_+- ~ +-sqrt(-q/p) and a finite sheet x_3 -> -r/q.
   The polar part of the trace observables factorizes, Vieta-exactly:
       e_3 = (x_+ x_-) * x_3  =  [q/p + O(1)] * [-r/q + O(p)],
   i.e. [divergent pair invariant] x [on-wall finite-sheet value] --
   the analogue of amplitude factorization on a pole into lower "amplitudes".
   Verified here by perturbative expansion of the roots in p.
"""

import sympy as sp

from jcqft import SRC, p, q, r

a, b, c = SRC

print("=== 1. Trace-map rationality: S_k = sum over sheets of x^k ===")
e1, e2, e3 = sp.Integer(0), q / p, -r / p
S = {0: sp.Integer(3), 1: e1}
S[2] = e1 * S[1] - 2 * e2
S[3] = e1 * S[2] - e2 * S[1] + 3 * e3
for k in range(4, 7):
    S[k] = e1 * S[k - 1] - e2 * S[k - 2] + e3 * S[k - 3]
for k in range(2, 7):
    num, den = sp.fraction(sp.cancel(sp.together(S[k])))
    fac = sp.factor_list(den)[1]
    dens = " * ".join(f"({sp.sstr(f)})^{m}" for f, m in fac)
    only_p = all(sp.simplify(f - p) == 0 for f, _ in fac)
    print(f"  S_{k}: denominator = {dens}   [poles only on p = 0: {only_p}]")
    assert only_p
print("  -> every sum-over-sheets observable is rational; poles confined to")
print("     the non-properness divisor.  Single sheets: multivalued (S3).")

print("\n=== 2. Boundary factorization near the wall p -> 0 ===")
eps = sp.Symbol("eps")  # stand-in for p as the wall-approach parameter
ORDER = 3
x3 = -r / q
for _ in range(ORDER + 1):  # iterate x = -(r + eps x^3)/q
    x3 = sp.expand(-(r + eps * x3**3) / q) + sp.O(eps**ORDER, eps)
    x3 = x3.removeO()
x3_series = sp.expand(x3 + sp.O(eps**2, eps)).removeO()
print(f"  finite sheet:    x_3(p->0) = {sp.sstr(sp.simplify(x3_series))} + O(p^2)")
expected = -r / q + eps * r**3 / q**4
assert sp.simplify(x3_series - expected) == 0

pair = sp.cancel((-r / eps) / x3_series)  # Vieta: x+ x- = e3_total / x3
pair_lead = sp.simplify(sp.limit(pair * eps, eps, 0))
print(f"  escaping pair:   x_+ x_-  =  ({sp.sstr(pair_lead)})/p + O(1)")
assert sp.simplify(pair_lead - q) == 0
print("  exact Vieta:     e_3 = (x_+ x_-) * x_3 = -r/p")
print("  polar factorization:  [q/p] * [-r/q]  =  -r/p  ✓")
print("  -> the pole of the trace observable factorizes into the escaping-pair")
print("     invariant times the on-wall value of the surviving sheet, the")
print("     analogue of amplitude factorization on a physical pole.")
