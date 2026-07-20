"""jcqft — the Alpöge–Mathew counterexample as a 0D quantum field theory.

Submodules:
    core       the map, propagator, eliminant cubic, chamber rule
    truncated  sparse truncated-ring arithmetic and the tree-graph inverse
    fibers     exact/numeric fiber computation via Groebner parametrization
    reduction  C*-equivariant normal form and the 2D-reduced Keller condition
"""

from jcqft.core import (
    D0,
    F,
    L,
    L_inv,
    PHI,
    SOURCE_WEIGHTS,
    SRC,
    TARGET_WEIGHTS,
    X,
    cubic,
    n_real_preimages,
    p,
    q,
    r,
)

__all__ = [
    "D0", "F", "L", "L_inv", "PHI", "SOURCE_WEIGHTS", "SRC",
    "TARGET_WEIGHTS", "X", "cubic", "n_real_preimages", "p", "q", "r",
]
