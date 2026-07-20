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

x, y, z = sp.symbols("x y z")

F = sp.Matrix(
    [
        (1 + x * y) ** 3 * z + y**2 * (1 + x * y) * (4 + 3 * x * y),
        y + 3 * x * (1 + x * y) ** 2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    ]
)

print("=== 1. Jacobian determinant ===")
J = F.jacobian([x, y, z])
detJ = sp.expand(J.det())
print(f"det DF = {detJ}")
assert detJ == -2, "Jacobian determinant is NOT constant -2!"

print("\n=== 2. Non-injectivity: three preimages of (-1/4, 0, 0) ===")
points = [
    (0, 0, sp.Rational(-1, 4)),
    (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
    (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
]
for p in points:
    image = tuple(sp.simplify(f.subs(dict(zip((x, y, z), p)))) for f in F)
    print(f"F{p} = {image}")
    assert image == (sp.Rational(-1, 4), 0, 0)

print("\n=== 3. Fiber over a generic target (a,b,c) = (1, 2, 3) ===")
eqs = [F[0] - 1, F[1] - 2, F[2] - 3]
gb = sp.groebner(eqs, x, y, z, order="lex")
sols = sp.solve(eqs, [x, y, z], dict=True)
print(f"Number of solutions over C: {len(sols)}")
for s in sols:
    approx = {k: sp.nsimplify(v, rational=False) for k, v in s.items()}
    print("  ", {k: sp.N(v, 8) for k, v in s.items()})

print("\nAll checks passed: F is a Keller map (det DF = -2) that is not injective.")
print("The Jacobian conjecture is FALSE in dimension >= 3.")
