"""Tensor de condutividade 2D: parte simbólica (SymPy) e numérica (NumPy).

Convenção do slide do Hands-On 1:  sigma' = R(theta) · sigma · R(theta)^T.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp

# Símbolos reutilizados em todo o projeto
theta = sp.Symbol("theta", real=True)
sigma_x = sp.Symbol("sigma_x", positive=True)
sigma_y = sp.Symbol("sigma_y", positive=True)


def sigma_symbolic(sx: sp.Expr = sigma_x, sy: sp.Expr = sigma_y) -> sp.Matrix:
    """Tensor diagonal nos eixos de coordenadas: diag(sigma_x, sigma_y)."""
    return sp.diag(sx, sy)


def rotation_matrix(t: sp.Expr = theta) -> sp.Matrix:
    """Matriz de rotação 2D R(t)."""
    return sp.Matrix([[sp.cos(t), -sp.sin(t)], [sp.sin(t), sp.cos(t)]])


def rotated_sigma(sigma: sp.Matrix, t: sp.Expr = theta, convention: str = "slide") -> sp.Matrix:
    """Componentes do tensor depois da rotação, já simplificadas.

    convention="slide": sigma' = R sigma R^T  (a do slide: tensor girado de theta).
    convention="axes" : sigma' = R^T sigma R  (componentes nos eixos x', y' girados de theta).
    As duas diferem apenas pelo sinal de theta.
    """
    R = rotation_matrix(t)
    if convention == "slide":
        m = R * sigma * R.T
    elif convention == "axes":
        m = R.T * sigma * R
    else:
        raise ValueError("convention deve ser 'slide' ou 'axes'")
    return m.applyfunc(sp.simplify)


def rotated_sigma_double_angle(sx: sp.Expr = sigma_x, sy: sp.Expr = sigma_y, t: sp.Expr = theta) -> sp.Matrix:
    """Mesma matriz sigma' escrita com ângulo duplo (forma mais reveladora)."""
    media = (sx + sy) / 2
    dif = (sx - sy) / 2
    return sp.Matrix(
        [
            [media + dif * sp.cos(2 * t), dif * sp.sin(2 * t)],
            [dif * sp.sin(2 * t), media - dif * sp.cos(2 * t)],
        ]
    )


def invariants_difference(sigma: sp.Matrix, sigma_p: sp.Matrix) -> tuple[sp.Expr, sp.Expr]:
    """Retorna (tr(sigma') - tr(sigma), det(sigma') - det(sigma)); ambos devem ser 0."""
    d_tr = sp.simplify(sigma_p.trace() - sigma.trace())
    d_det = sp.simplify(sigma_p.det() - sigma.det())
    return d_tr, d_det


@lru_cache(maxsize=None)
def _sigma_prime_function():
    """Converte a expressão simbólica de sigma' numa função NumPy (SymPy -> NumPy)."""
    expr = rotated_sigma(sigma_symbolic())
    return sp.lambdify((theta, sigma_x, sigma_y), expr, "numpy")


def sigma_numeric(sx_val: float, sy_val: float, theta_deg: float) -> np.ndarray:
    """sigma'(theta) como array NumPy 2x2, com theta em graus."""
    f = _sigma_prime_function()
    return np.array(f(np.deg2rad(theta_deg), sx_val, sy_val), dtype=float)


def current_density(sigma_num: np.ndarray, E) -> np.ndarray:
    """J = sigma' E."""
    return sigma_num @ np.asarray(E, dtype=float)


def angle_between_deg(u, v) -> float:
    """Ângulo com sinal (graus) de u para v, em (-180, 180]."""
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    cross = u[0] * v[1] - u[1] * v[0]
    dot = float(u @ v)
    return float(np.degrees(np.arctan2(cross, dot)))
