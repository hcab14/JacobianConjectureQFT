"""Core definitions: the Alpöge–Mathew counterexample to the Jacobian
conjecture (announced July 19, 2026) and its 0D-QFT reading.

The map F : C^3 -> C^3 has constant Jacobian determinant -2 but is not
injective.  In the field-theory picture the components of F are the classical
equations of motion for fields phi = (x, y, z) coupled to external sources
J = (a, b, c); DF(0) = L is the inverse propagator and the nonlinear terms of
F are the interaction vertices.
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


# C*-equivariance ("orbifold mechanism"):
#   F(l*x, y/l, z/l^2) = (F1/l^2, F2/l, l*F3).
SOURCE_WEIGHTS = (1, -1, -2)
TARGET_WEIGHTS = (-2, -1, 1)

# x-eliminant (Groebner elimination of y, z from F(phi) = (a,b,c)):
# every preimage of a target (a, b, c) has x-coordinate satisfying
#     p(a,b,c) * X^3 + q(a,b,c) * X + r(a,b,c) = 0.
X = sp.Symbol("X")
p = 27 * a**2 * c**2 - 18 * a * b * c + b**3 * c - b**2 + 16 * a
q = 4 - 3 * b * c
r = -2 * c
cubic = p * X**3 + q * X + r

# Discriminant structure (verified in scripts/branch_locus.py, monodromy.py):
#   disc_X(cubic) = -4 * D0^2 * p,   D0 = denominator of the y,z fiber
# parametrization.  {p = 0} is the non-properness (escape) locus; {D0 = 0}
# is a harmless x-collision locus (the covering is etale everywhere).
D0 = 27 * a * c**2 - 9 * b * c + 8


def n_real_preimages(pa: float) -> int:
    """Chamber rule for the REAL map R^3 -> R^3: the number of real
    preimages of a target with p-value `pa` is 3 iff p < 0, else 1
    (monic discriminant of the cubic is -4*D0^2/p^3)."""
    return 3 if pa < 0 else 1
