"""Regression tests for the checks THERMPY's README labels 'Validated'.

Each test pins a number the docs already quote, so a change that silently
breaks a validated result fails CI instead of shipping.
"""
import numpy as np
import pytest

import conduction as cd
import environment as env

Q_ABS = 98.0   # W/m2, worst-beta deployed-edge-on panel at 800 km (radiator_fin_submodel.md)


def test_view_factor_integrator_matches_sin2rho():
    # README: integrator vs sin^2(rho) at theta = 0 to 8e-6
    assert env.validate_view_factor() < 1e-5


def test_eclipse_fraction_matches_phase1():
    # README / altitude_selection_30deg.md: 500 km, beta = 0 -> 37.75%, 35.7 min
    c = env.Case(500, 0.0)
    f = env.eclipse_fraction(c)
    assert f == pytest.approx(0.3775, abs=5e-4)
    assert f * c.period_s / 60 == pytest.approx(35.7, abs=0.2)


@pytest.mark.parametrize("L", [0.050, 0.100, 0.200])
def test_fin_collocation_matches_newton_fd(L):
    # radiator_fin_submodel.md section 3.3: -0.015% at 50 mm to -0.051% at 200 mm
    p = cd.Panel(half_length=L, q_abs=Q_ABS)
    eta_bvp, _ = cd.fin_efficiency(p)
    eta_fd, _, ok = cd.fd_fin_newton(p)
    assert ok, "Newton FD did not converge"
    assert abs(eta_bvp - eta_fd) / eta_fd < 6e-4


needs_2d = pytest.mark.skipif(not hasattr(cd, "fd_fin_2d_end"),
                              reason="fd_fin_2d_end not yet merged into conduction.py")


@needs_2d
def test_2d_with_pipe_to_edge_is_exactly_1d():
    # No overhang -> the 2-D problem has a 1-D solution; per-y net must be flat
    # and match the 1-D result to discretization.
    p = cd.Panel(half_length=0.088, t_face=2.5e-4, q_abs=Q_ABS)
    r = cd.fd_fin_2d_end(p, overhang_m=0.0)
    assert np.ptp(r["net_per_y"]) / r["net_per_y"].mean() < 1e-10
    assert r["net_per_area"] == pytest.approx(cd.net_per_area_1d(p), rel=5e-4)


@needs_2d
def test_2d_end_effect_grid_converged():
    # 88 mm overhang: equivalent lost pipe length ~38 mm, stable under refinement
    p = cd.Panel(half_length=0.088, t_face=2.5e-4, q_abs=Q_ABS)
    coarse = cd.fd_fin_2d_end(p, overhang_m=0.088, h=0.002)["lost_length_m"]
    fine = cd.fd_fin_2d_end(p, overhang_m=0.088, h=0.001)["lost_length_m"]
    assert fine == pytest.approx(0.0384, abs=1e-3)
    assert abs(coarse - fine) / fine < 0.01
