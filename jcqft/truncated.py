"""Sparse truncated polynomial arithmetic and the perturbative (tree-graph)
inversion of the counterexample map.

Polynomials are dicts {exponent tuple: Rational coefficient}; all operations
truncate at a fixed total degree, so computations stay in the truncated ring.
The Picard iteration

    phi^(0)   = L^{-1} J
    phi^(k+1) = phi^(k) + L^{-1} (J - F(phi^(k)))

is exactly the sum over rooted tree Feynman graphs with up to N leaves.
"""

import sympy as sp

from jcqft.core import F, PHI


def from_expr(expr, gens):
    poly = sp.Poly(sp.expand(expr), *gens)
    return {m: v for m, v in poly.as_dict().items() if v}


def to_expr(d, gens):
    return sp.Add(
        *[v * sp.Mul(*[g**e for g, e in zip(gens, m)]) for m, v in d.items()]
    )


def trunc(d, N):
    return {m: v for m, v in d.items() if sum(m) <= N}


def padd(f, g):
    out = dict(f)
    for m, v in g.items():
        s = out.get(m, 0) + v
        if s:
            out[m] = s
        elif m in out:
            del out[m]
    return out


def pscale(f, s):
    if not s:
        return {}
    return {m: v * s for m, v in f.items()}


def pmul(f, g, N):
    out = {}
    for m1, v1 in f.items():
        d1 = sum(m1)
        for m2, v2 in g.items():
            if d1 + sum(m2) > N:
                continue
            m = tuple(i + j for i, j in zip(m1, m2))
            s = out.get(m, 0) + v1 * v2
            if s:
                out[m] = s
            elif m in out:
                del out[m]
    return out


F_DICTS = [from_expr(f, PHI) for f in F]
# highest power of each field appearing in any component of F
MAX_POW = [max(m[i] for fd in F_DICTS for m in fd) for i in range(3)]


def compose_F(phi, N, nsrc):
    """Evaluate F at phi = (phi_x, phi_y, phi_z), each a truncated series."""
    one = {(0,) * nsrc: sp.Integer(1)}
    pows = []
    for i in range(3):
        lst = [one]
        for _ in range(MAX_POW[i]):
            lst.append(pmul(lst[-1], phi[i], N))
        pows.append(lst)
    out = []
    for fd in F_DICTS:
        acc = {}
        for (ex, ey, ez), coeff in fd.items():
            term = pscale(
                pmul(pmul(pows[0][ex], pows[1][ey], N), pows[2][ez], N), coeff
            )
            acc = padd(acc, term)
        out.append(acc)
    return tuple(out)


def formal_inverse(J, N, nsrc, max_iter=None):
    """Truncated formal inverse: phi with F(phi) = J mod degree N+1."""
    phi = tuple(trunc(d, N) for d in (pscale(J[2], sp.Rational(1, 2)), J[1], J[0]))
    max_iter = max_iter or N + 5
    for _ in range(max_iter):
        err = tuple(
            padd(Ji, pscale(Fi, -1)) for Ji, Fi in zip(J, compose_F(phi, N, nsrc))
        )
        if all(not e for e in err):
            return phi
        corr = (pscale(err[2], sp.Rational(1, 2)), err[1], err[0])
        phi = tuple(trunc(padd(ph, co), N) for ph, co in zip(phi, corr))
    raise RuntimeError("Picard iteration did not reach the truncated fixed point")
