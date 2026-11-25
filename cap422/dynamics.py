"""Dinâmica numérica (Tópico 2): poços de potencial e pêndulo com scipy.integrate.solve_ivp."""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk

RTOL, ATOL = 1e-10, 1e-12


def integrate(rhs, s0, t_end: float, n: int = 4000, events=None):
    """Integra ds/dt = rhs(t, s) de 0 a t_end com tolerância apertada (DOP853)."""
    t_eval = np.linspace(0.0, t_end, n)
    return solve_ivp(rhs, (0.0, t_end), s0, t_eval=t_eval, method="DOP853",
                     rtol=RTOL, atol=ATOL, events=events)


# ---------- Hands-On 2: partícula em poço de potencial, com atrito de Rayleigh ----------
def well_rhs(dV, m: float = 1.0, gamma: float = 0.0):
    """Sistema de 1ª ordem para  m x'' + V'(x) + gamma x' = 0  (estado = [x, v])."""
    def rhs(t, s):
        x, v = s
        return [v, -(dV(x) + gamma * v) / m]
    return rhs


def mech_energy(x, v, V, m: float = 1.0):
    return 0.5 * m * np.asarray(v) ** 2 + V(np.asarray(x))


# ---------- Hands-On 3: pêndulo simples ----------
def pendulum_rhs(g: float = 9.81, l: float = 1.0, linear: bool = False, damping: float = 0.0):
    """theta'' = -(g/l) sin(theta)  (ou  -(g/l) theta  se linear) - damping * theta'."""
    w2 = g / l

    def rhs(t, s):
        th, om = s
        return [om, -w2 * (th if linear else np.sin(th)) - damping * om]
    return rhs


def pendulum_energy(th, om, m: float = 1.0, g: float = 9.81, l: float = 1.0):
    """E = T + V, com V = 0 na posição mais baixa."""
    return 0.5 * m * l**2 * np.asarray(om) ** 2 + m * g * l * (1.0 - np.cos(th))


def period_exact(theta0: float, g: float = 9.81, l: float = 1.0) -> float:
    """Período exato do pêndulo não linear: T = 4 K(m)/omega0, m = sin^2(theta0/2)."""
    return 4.0 * ellipk(np.sin(theta0 / 2.0) ** 2) * np.sqrt(l / g)


def period_small_angle(g: float = 9.81, l: float = 1.0) -> float:
    return 2.0 * np.pi * np.sqrt(l / g)


def period_measured(rhs, theta0: float, t_end: float) -> float:
    """Período medido numericamente: tempo entre duas passagens consecutivas
    por theta = 0 no mesmo sentido (descendo)."""
    def cross(t, s):
        return s[0]
    cross.direction = -1.0
    sol = integrate(rhs, [theta0, 0.0], t_end, n=10, events=cross)
    te = sol.t_events[0]
    if len(te) < 2:
        raise RuntimeError("t_end curto demais para medir o período")
    return float(te[1] - te[0])
