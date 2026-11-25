import numpy as np
import pytest
import sympy as sp
from sympy.calculus.euler import euler_equations

from cap422.dynamics import (integrate, mech_energy, pendulum_energy, pendulum_rhs,
                             period_exact, period_measured, period_small_angle, well_rhs)
from cap422.variational import euler_lagrange, euler_lagrange_with_rayleigh

x = sp.Symbol("x", real=True)
y = sp.Function("y")(x)
yp = sp.diff(y, x)


def test_euler_lagrange_oscilador_harmonico():
    m, k = sp.symbols("m k", positive=True)
    L = sp.Rational(1, 2) * m * yp**2 - sp.Rational(1, 2) * k * y**2
    assert sp.simplify(euler_lagrange(L, y, x) - (m * sp.diff(y, x, 2) + k * y)) == 0


@pytest.mark.parametrize("L", [
    sp.Rational(1, 2) * yp**2 + y**2,
    yp**2 * y + sp.sin(y),
    sp.sqrt(1 + yp**2),  # comprimento de arco: reta
])
def test_euler_lagrange_confere_com_sympy(L):
    ref = euler_equations(L, y, x)[0].lhs  # sinal oposto ao nosso
    assert sp.simplify(euler_lagrange(L, y, x) + ref) == 0


def test_rayleigh_adiciona_termo_de_atrito():
    t = sp.Symbol("t")
    q = sp.Function("q")(t)
    g = sp.Symbol("gamma", positive=True)
    L = sp.Rational(1, 2) * sp.diff(q, t) ** 2 - sp.Rational(1, 2) * q**2
    R = sp.Rational(1, 2) * g * sp.diff(q, t) ** 2
    eq = euler_lagrange_with_rayleigh(L, q, t, R)
    assert sp.simplify(eq - (sp.diff(q, t, 2) + q + g * sp.diff(q, t))) == 0


def test_energia_conservada_sem_atrito_e_decresce_com_atrito():
    V = lambda x: 0.25 * (x**2 - 1) ** 2
    dV = lambda x: x * (x**2 - 1)
    sol = integrate(well_rhs(dV, 1.0, 0.0), [1.6, 0.0], 30.0)
    E = mech_energy(sol.y[0], sol.y[1], V)
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-8
    sol = integrate(well_rhs(dV, 1.0, 0.2), [1.6, 0.0], 30.0)
    E = mech_energy(sol.y[0], sol.y[1], V)
    assert np.all(np.diff(E) < 1e-9)


@pytest.mark.parametrize("ang", [5, 30, 90])
def test_periodo_do_pendulo_bate_com_integral_eliptica(ang):
    th0 = np.radians(ang)
    medido = period_measured(pendulum_rhs(), th0, 10.0)
    assert medido == pytest.approx(period_exact(th0), rel=1e-6)


def test_pequeno_angulo_e_energia_do_pendulo():
    th0 = np.radians(5)
    assert period_exact(th0) == pytest.approx(period_small_angle(), rel=1e-3)
    sol = integrate(pendulum_rhs(), [np.radians(90), 0.0], 10.0)
    E = pendulum_energy(sol.y[0], sol.y[1])
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-8
