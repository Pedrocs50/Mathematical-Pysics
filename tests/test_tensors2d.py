import numpy as np
import sympy as sp

from cap422 import tensors2d as t


def test_rotation_is_orthogonal():
    R = t.rotation_matrix()
    assert sp.simplify(R * R.T - sp.eye(2)) == sp.zeros(2)
    assert sp.simplify(R.det()) == 1


def test_invariants_symbolic():
    s = t.sigma_symbolic()
    sp_ = t.rotated_sigma(s)
    d_tr, d_det = t.invariants_difference(s, sp_)
    assert d_tr == 0 and d_det == 0


def test_double_angle_form_matches():
    a = t.rotated_sigma(t.sigma_symbolic())
    b = t.rotated_sigma_double_angle()
    assert sp.simplify(a - b) == sp.zeros(2)


def test_numeric_values():
    assert np.allclose(t.sigma_numeric(2, 1, 0), [[2, 0], [0, 1]])
    s30 = t.sigma_numeric(2, 1, 30)
    assert np.allclose(s30, [[1.75, np.sqrt(3) / 4], [np.sqrt(3) / 4, 1.25]])
    assert np.allclose(t.sigma_numeric(2, 1, 45), [[1.5, 0.5], [0.5, 1.5]])


def test_numeric_invariants():
    for th in (0, 30, 45, 77):
        s = t.sigma_numeric(2, 1, th)
        assert np.isclose(np.trace(s), 3)
        assert np.isclose(np.linalg.det(s), 2)


def test_J_for_E_11():
    E = (1, 1)
    assert np.allclose(t.current_density(t.sigma_numeric(2, 1, 0), E), [2, 1])
    # a 45 graus, E=(1,1) é autovetor (autovalor 2): J paralelo a E
    J45 = t.current_density(t.sigma_numeric(2, 1, 45), E)
    assert np.allclose(J45, [2, 2])
    assert np.isclose(t.angle_between_deg(E, J45), 0)


def test_isotropic_is_rotation_invariant():
    for th in (0, 13, 30, 45, 90):
        assert np.allclose(t.sigma_numeric(1.5, 1.5, th), 1.5 * np.eye(2))


def test_axes_convention_is_minus_theta():
    s = t.sigma_symbolic()
    a = t.rotated_sigma(s, convention="axes")
    b = t.rotated_sigma(s, -t.theta, convention="slide")
    assert sp.simplify(a - b) == sp.zeros(2)
