"""Shared definitions for the Alpöge–Mathew counterexample to the Jacobian
conjecture (announced July 19, 2026) and its 0D-QFT reformulation.

The map F : C^3 -> C^3 has constant Jacobian determinant -2 but is not
injective.  In the 0D field-theory picture the components of F are the
classical equations of motion dS/dphi_i = J_i for fields phi = (x, y, z)
coupled to external sources J = (a, b, c); inverting F perturbatively is
the sum over rooted tree Feynman graphs.
"""

import sympy as sp

x, y, z = sp.symbols("x y z")
a, b, c = sp.symbols("a b c")
PHI = (x, y, z)
SRC = (a, b, c)

F = (
    (1 + x * y) ** 3 * z + y**2 * (1 + x * y) * (4 + 3 * x * y),
    y + 3 * x * (1 + x * y) ** 2 * z + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
)

# Linearization at the origin ("free propagator"):
#   DF(0) = L with L(x, y, z) = (z, y, 2x),  so  L^{-1}(a, b, c) = (c/2, b, a).
L = sp.Matrix([[0, 0, 1], [0, 1, 0], [2, 0, 0]])


def L_inv(v):
    """Apply L^{-1} to a 3-tuple (works for sympy expressions)."""
    return (v[2] / 2, v[1], v[0])


# x-eliminant (computed by Groebner elimination of y, z from F(phi) = (a,b,c)):
# every preimage (x, y, z) of a target (a, b, c) has x-coordinate satisfying
#     p(a,b,c) * X^3 + q(a,b,c) * X + r(a,b,c) = 0.
X = sp.Symbol("X")
p = 27 * a**2 * c**2 - 18 * a * b * c + b**3 * c - b**2 + 16 * a
q = 4 - 3 * b * c
r = -2 * c
cubic = p * X**3 + q * X + r
