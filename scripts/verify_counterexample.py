"""Independent verification of the Alpöge–Mathew counterexample to the
Jacobian conjecture (announced July 19, 2026).

Checks:
  1. The Jacobian determinant of F is identically the constant -2.
  2. Three distinct points of C^3 map to the same target (-1/4, 0, 0),
     so F is not injective, hence not invertible.
  3. (Diagnostic) fiber over a generic rational target point, to probe the
     generic degree of the map.
"""

import sympy as sp

from jcqft import F, PHI
from jcqft.fibers import exact_fiber

print("=== 1. Jacobian determinant ===")
detJ = sp.expand(sp.Matrix(F).jacobian(PHI).det())
print(f"det DF = {detJ}")
assert detJ == -2, "Jacobian determinant is NOT constant -2!"

print("\n=== 2. Non-injectivity: three preimages of (-1/4, 0, 0) ===")
points = [
    (0, 0, sp.Rational(-1, 4)),
    (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
    (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
]
for pt in points:
    image = tuple(sp.simplify(f.subs(dict(zip(PHI, pt)))) for f in F)
    print(f"F{pt} = {image}")
    assert image == (sp.Rational(-1, 4), 0, 0)

print("\n=== 3. Fiber over a generic target (a,b,c) = (1, 2, 3) ===")
fiber = exact_fiber((1, 2, 3))
print(f"Number of solutions over C: {len(fiber)}")
for pt in fiber:
    print("  ", tuple(sp.N(v, 8) for v in pt))

print("\nAll checks passed: F is a Keller map (det DF = -2) that is not injective.")
print("The Jacobian conjecture is FALSE in dimension >= 3.")
