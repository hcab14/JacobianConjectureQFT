"""Witten-index / infinity prefilter for counterexample searches.

(Implements docs/SEARCH_STRATEGIES.md §1.3 = docs/OPEN_QUESTIONS.md A3;
validation in scripts/witten_prefilter.py, write-up in
docs/RIGIDITY_AND_PREFILTER.md.)

MATHEMATICAL BASIS.  A counterexample to the Jacobian conjecture must be
non-proper (nonempty Jelonek non-properness set: a proper Keller map is a
finite etale covering of the simply connected C^n, hence injective).
Properness has a cheap sufficient criterion at infinity.  Let L_i be the
leading form of F_i (top-degree homogeneous part).  If the L_i have no
common zero other than the origin, then by homogeneity and compactness of
the unit sphere, along every direction u at infinity some |F_i(r u)| grows
like r^{deg F_i}: F is PROPER, |F| -> infinity whenever |phi| -> infinity,
and F cannot be a counterexample.

Contrapositive (the prefilter): F can be non-proper only if its leading
forms degenerate at the hyperplane at infinity, i.e. share a common
projective zero.  Testing this is one small Groebner/linear-algebra
computation on the leading forms — no symbolic Keller work — so thousands
of candidates can be screened before anything expensive runs.

    infinity_prefilter(F, vars) == False  =>  F provably proper: REJECT.
    infinity_prefilter(F, vars) == True   =>  candidate survives (the test
                                              is necessary, NOT sufficient,
                                              for non-properness).

KNOWN FALSE-POSITIVE CLASS (unavoidable): every polynomial automorphism of
degree > 1 survives.  (If its leading forms had only the trivial common
zero, Bezout would give generic fiber cardinality prod_i deg F_i > 1,
contradicting injectivity.)  The filter therefore rejects linear maps and
the Bezout-generic bulk of a search space; it cannot separate nonlinear
automorphisms from counterexamples — that is the Keller search's job.

WEIGHTED VARIANT.  For a C*-equivariant search class, pass the source
weight vector (entries may be negative, e.g. the Alpoge-Mathew weights
(1,-1,-2)).  The relevant escape curves are then weighted orbits
phi_j(lam) = lam^{w_j} c_j with lam -> infinity, which escape iff c_j != 0
for some j of positive weight.  Along such a curve, F_i(phi(lam)) =
lam^{d_i} L_i(c) + lower order, where L_i is the w-leading part of F_i and
d_i its weighted degree; F_i stays bounded to leading order iff d_i <= 0 or
L_i(c) = 0.  The weighted test therefore asks: do the leading parts of the
POSITIVE-weighted-degree components share a zero c having some nonzero
positive-weight coordinate?  (Both scaling directions lam -> infinity and
lam -> 0 are checked, the latter via the negated weight vector.)  A False
verdict rules out first-order escape along w-orbits only — for a full
properness proof use the unweighted call, which is the w = (1,...,1) case
and covers all escape directions.  Note that on maps that are exactly
equivariant (each component weighted-homogeneous, as produced by
jcqft.reduction), the w-leading part is the whole component.
"""

from __future__ import annotations

import sympy as sp


def leading_part(f, variables, weights=None):
    """Top (weighted-)degree part of the polynomial f, and its degree.

    weights: per-variable integer weights (default: all 1, the standard
    total degree).  The weighted degree of a monomial prod x_i^{e_i} is
    sum_i weights[i] * e_i.  Returns (leading_part, weighted_degree);
    (0, None) for f == 0.
    """
    f = sp.expand(f)
    if f == 0:
        return sp.S.Zero, None
    if weights is None:
        weights = (1,) * len(variables)
    terms = sp.Poly(f, *variables).terms()
    wdeg = lambda mono: sum(wi * ei for wi, ei in zip(weights, mono))
    top = max(wdeg(m) for m, _ in terms)
    lead = sp.Add(*[
        coeff * sp.Mul(*[xv**e for xv, e in zip(variables, m)])
        for m, coeff in terms if wdeg(m) == top
    ])
    return lead, top


def _survives(F_components, variables, weights):
    """One scaling direction (lam -> infinity for the given weights): is
    there a common zero c of the positive-degree leading parts with some
    positive-weight coordinate of c nonzero?"""
    escape_vars = [xv for xv, wi in zip(variables, weights) if wi > 0]
    if not escape_vars:
        return False  # no coordinate escapes under this scaling
    forms = []
    for f in F_components:
        lead, deg = leading_part(f, variables, weights)
        if deg is not None and deg > 0:
            forms.append(lead)
    if not forms:
        return True  # nothing blows up: every escape direction survives
    gb = sp.groebner(forms, *variables, order="grevlex")
    if list(gb.exprs) == [1]:
        return False  # empty zero set
    if all(wi == weights[0] for wi in weights):
        # equal weights: forms homogeneous, V a cone through 0, and every
        # coordinate is an escape coordinate, so a nontrivial common zero
        # exists iff dim V >= 1
        return not gb.is_zero_dimensional
    # mixed weights: Rabinowitsch test "exists c in V with c_j != 0" for
    # each escape coordinate x_j:  1 not in I + (1 - t*x_j)
    t = sp.Dummy("t")
    for xj in escape_vars:
        gb_j = sp.groebner(list(gb.exprs) + [1 - t * xj], *variables, t,
                           order="grevlex")
        if list(gb_j.exprs) != [1]:
            return True
    return False


def infinity_prefilter(F_components, variables, weights=None):
    """Witten-index prefilter.  True = candidate survives (its leading
    forms degenerate at infinity); False = reject.

    Unweighted (weights=None): False PROVES F proper, hence not a
    counterexample.  Weighted: False rules out first-order escape along
    orbits of the given C* weight vector (both lam -> infinity and
    lam -> 0); see the module docstring for the exact scope.
    """
    if weights is None:
        return _survives(F_components, variables, (1,) * len(variables))
    return (_survives(F_components, variables, tuple(weights))
            or _survives(F_components, variables,
                         tuple(-wi for wi in weights)))


def infinity_witness(F_components, variables, weights=None):
    """A witness for a surviving candidate: a point c with some nonzero
    positive-weight coordinate at which all positive-degree leading parts
    vanish — the direction at infinity along which properness can fail.
    Returns (witness_tuple, scaling_sign) with scaling_sign = +1 for
    lam -> infinity, -1 for the negated weights, or None if the candidate
    is rejected.  Chart-by-chart sp.solve; for validation and reporting
    (the fast screening path is infinity_prefilter)."""
    if weights is None:
        weights = (1,) * len(variables)
    for sign in (+1, -1):
        wts = tuple(sign * wi for wi in weights)
        escape_vars = [xv for xv, wi in zip(variables, wts) if wi > 0]
        forms = []
        for f in F_components:
            lead, deg = leading_part(f, variables, wts)
            if deg is not None and deg > 0:
                forms.append(sp.expand(lead))
        for chart in escape_vars:
            others = [xv for xv in variables if xv != chart]
            eqs = [f.subs(chart, 1) for f in forms]
            try:
                sols = sp.solve(eqs, others, dict=True)
            except Exception:
                continue
            if not forms:
                sols = [dict()]  # unconstrained: any point on the chart
            for s in sols:
                pt = {chart: sp.Integer(1)}
                for xv in others:
                    val = sp.sympify(s.get(xv, sp.Integer(0)))
                    if val.free_symbols:
                        val = val.subs({fs: 0 for fs in val.free_symbols})
                    pt[xv] = val
                if all(sp.expand(f.subs(pt)) == 0 for f in forms):
                    return tuple(pt[xv] for xv in variables), sign
    return None
