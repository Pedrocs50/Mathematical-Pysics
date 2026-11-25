"""Cálculo variacional simbólico (Tópico 2): equação de Euler-Lagrange com SymPy."""
from __future__ import annotations

import sympy as sp


def euler_lagrange(L: sp.Expr, y: sp.Function, x: sp.Symbol, simplify: bool = True) -> sp.Expr:
    """Lado esquerdo da equação de Euler-Lagrange:  d/dx (dL/dy') - dL/dy.

    `y` deve ser uma função simbólica de x (ex.: y = sp.Function("y")(x)).
    A equação de movimento é  euler_lagrange(L, y, x) = 0.
    """
    yp = sp.diff(y, x)
    expr = sp.diff(sp.diff(L, yp), x) - sp.diff(L, y)
    return sp.simplify(expr) if simplify else expr


def action(L: sp.Expr, x: sp.Symbol, x0, x1) -> sp.Integral:
    """Ação S[y] = integral de L(x, y, y') dx entre x0 e x1 (não calculada)."""
    return sp.Integral(L, (x, x0, x1))


def euler_lagrange_with_rayleigh(L: sp.Expr, q: sp.Function, t: sp.Symbol, R: sp.Expr) -> sp.Expr:
    """Euler-Lagrange com dissipação de Rayleigh:  d/dt(dL/dq') - dL/dq + dR/dq' = 0."""
    qp = sp.diff(q, t)
    return sp.simplify(sp.diff(sp.diff(L, qp), t) - sp.diff(L, q) + sp.diff(R, qp))
