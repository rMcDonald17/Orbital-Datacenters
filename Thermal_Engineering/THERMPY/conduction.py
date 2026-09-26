"""
conduction.py -- Fourier sub-models (Part C Layer 3).

Analytical conduction models that produce the lumped parameters Layer 2 needs,
each checked against a discretized solve. The radiator fin is the important one:
it sets the effective radiating temperature (area goes as T^-4) and, through the
fin mass optimization, it replaces sizing.MassModel.radiator_kg_per_m2 -- the
highest-leverage unsourced assumption in the study.

Panel model: aluminium honeycomb sandwich, two facesheets, heat pipes embedded
at spacing 2L, radiating from both faces. Conduction is carried by the
facesheets only; the core is mass without conductance.

Per-slice energy balance on a fin of half-length L:

    2 k t_f T'' = 2 eps sigma (T^4 - T_sp^4)      ->   k t_f T'' = eps sig (T^4 - T_sp^4)

with T(0) = T_root and T'(L) = 0 (adiabatic tip, symmetry between heat pipes).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_bvp
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

SIGMA = 5.670374419e-8
T_SPACE = 3.0

MATERIALS = {                       # k W/m/K, rho kg/m3
    'Al 6061-T6':   dict(k=167.0, rho=2700.0),
    'Al 2024-T3':   dict(k=121.0, rho=2780.0),
    'Al 1100':      dict(k=222.0, rho=2710.0),
    'K1100 composite': dict(k=500.0, rho=1800.0),   # pitch fibre, in-plane
}


@dataclass
class Panel:
    """Honeycomb sandwich radiator panel."""
    material: str = 'Al 6061-T6'
    t_face: float = 0.0003          # m, per facesheet
    half_length: float = 0.075      # m, heat-pipe half-spacing
    core_height: float = 0.015      # m
    core_density: float = 50.0      # kg/m3
    heatpipe_kg_per_m: float = 0.35 # axially grooved Al/ammonia
    misc_kg_per_m2: float = 0.5     # coating, adhesive, doublers
    eps: float = 0.85
    t_root: float = 318.0
    q_abs: float = 0.0              # mean absorbed flux over the radiating faces, W/m2

    @classmethod
    def from_environment(cls, alt_km, beta_deg=None, config='deployed-edge-on', **kw):
        """
        Build a Panel whose q_abs comes from environment.panel_absorbed_flux
        rather than a hardcoded number, so the two stay consistent when the
        view-factor iteration lands. Pass coupling=Coupling(...) through kw to
        include inter-surface terms once geometry is defined.
        """
        import environment as env
        panel_kw = {k: kw.pop(k) for k in list(kw)
                    if k in cls.__dataclass_fields__ and k != 'q_abs'}
        t_root = panel_kw.get('t_root', cls.t_root)
        if beta_deg is None:                      # worst beta for this config
            betas = np.linspace(0.0, env.beta_max(), 28)
            pf = max((env.panel_absorbed_flux(alt_km, b, config, t_root, **kw)
                      for b in betas), key=lambda r: r.mean)
        else:
            pf = env.panel_absorbed_flux(alt_km, beta_deg, config, t_root, **kw)
        return cls(q_abs=pf.mean, **panel_kw)

    @property
    def t_equilibrium(self):
        """Fin asymptote: net removal reaches zero here. Bounds useful fin length."""
        if self.q_abs <= 0:
            return T_SPACE
        return (self.q_abs / (self.eps * SIGMA) + T_SPACE ** 4) ** 0.25

    @property
    def k(self): return MATERIALS[self.material]['k']

    @property
    def rho(self): return MATERIALS[self.material]['rho']

    @property
    def areal_mass(self):
        """kg per m^2 of PHYSICAL panel."""
        faces = 2 * self.rho * self.t_face
        core = self.core_density * self.core_height
        pipes = self.heatpipe_kg_per_m / (2 * self.half_length)
        return faces + core + pipes + self.misc_kg_per_m2


# ------------------------------------------------------------ linearized fin

def fin_efficiency_linear(p: Panel, t_bar=None):
    """
    Standard fin solution with a radiative film coefficient.

        h_rad = 4 eps sigma Tbar^3,  m = sqrt(h_rad / (k t_f)),  eta = tanh(mL)/(mL)

    Exact only for small dT along the fin -- which a radiator is not.
    """
    t_bar = t_bar or p.t_root
    h_rad = 4 * p.eps * SIGMA * t_bar ** 3
    m = np.sqrt(h_rad / (p.k * p.t_face))
    mL = m * p.half_length
    return float(np.tanh(mL) / mL), float(mL)


# ------------------------------------------------------------ nonlinear fin

def fin_profile(p: Panel, n=200, tol=1e-6):
    """
    Solve the nonlinear BVP with the absorbed flux as a distributed source:

        k t_f T'' = eps sigma (T^4 - T_sp^4) - q_abs

    A uniform initial guess at T_root exceeds the mesh node limit without
    converging -- the solver refines the mesh rather than moving the profile.
    The linearized solution is used as the starting iterate instead.
    """
    c = p.eps * SIGMA / (p.k * p.t_face)
    s = p.q_abs / (p.k * p.t_face)
    L = p.half_length

    rhs = lambda x, y: np.vstack([y[1], c * (y[0] ** 4 - T_SPACE ** 4) - s])
    bc = lambda ya, yb: np.array([ya[0] - p.t_root, yb[1]])

    t_eq = p.t_equilibrium
    m = np.sqrt(4 * p.eps * SIGMA * p.t_root ** 3 / (p.k * p.t_face))
    x = np.linspace(0, L, n)
    g = t_eq + (p.t_root - t_eq) * np.cosh(m * (L - x)) / np.cosh(m * L)
    dg = -(p.t_root - t_eq) * m * np.sinh(m * (L - x)) / np.cosh(m * L)

    sol = solve_bvp(rhs, bc, x, np.vstack([g, dg]), tol=tol, max_nodes=100000)
    if not sol.success:
        raise RuntimeError(f"fin BVP failed at L={L*1e3:.0f} mm: {sol.message}")
    xx = np.linspace(0, L, 400)
    return xx, sol.sol(xx)[0]


def fin_efficiency(p: Panel):
    """
    NET fin efficiency: net heat actually removed, over the net an isothermal
    fin at T_root would remove. q_abs appears in BOTH numerator and denominator
    -- that is why the source term makes efficiency slightly worse rather than
    better (see docs/radiator_fin_submodel.md §4).

    Returns (eta_net, tip temperature).
    """
    x, T = fin_profile(p)
    actual = np.trapezoid(p.eps * SIGMA * (T ** 4 - T_SPACE ** 4) - p.q_abs, x)
    ideal = (p.eps * SIGMA * (p.t_root ** 4 - T_SPACE ** 4) - p.q_abs) * p.half_length
    return float(actual / ideal), float(T[-1])


def fd_fin_newton(p: Panel, n=401):
    """
    Independent finite-difference solve of the same BVP, Newton with an analytic
    Jacobian. Exists so a wrong answer from collocation does not go unnoticed.

    Do NOT substitute Jacobi iteration here: it converges as n^2 and returns an
    unconverged profile that looks like a physical result rather than raising.
    """
    from scipy.optimize import root
    L = p.half_length
    x = np.linspace(0, L, n); dx = x[1] - x[0]
    c = p.eps * SIGMA / (p.k * p.t_face)
    s = p.q_abs / (p.k * p.t_face)

    def resid(T):
        R = np.empty_like(T)
        R[0] = T[0] - p.t_root
        R[1:-1] = (T[2:] - 2 * T[1:-1] + T[:-2]) / dx ** 2 - c * (T[1:-1] ** 4 - T_SPACE ** 4) + s
        R[-1] = (T[-1] - T[-2]) / dx
        return R

    def jac(T):
        J = np.zeros((len(T), len(T)))
        J[0, 0] = 1.0
        idx = np.arange(1, len(T) - 1)
        J[idx, idx - 1] = 1 / dx ** 2
        J[idx, idx + 1] = 1 / dx ** 2
        J[idx, idx] = -2 / dx ** 2 - 4 * c * T[idx] ** 3
        J[-1, -1] = 1 / dx; J[-1, -2] = -1 / dx
        return J

    sol = root(resid, np.full(n, p.t_root * 0.95), jac=jac, method='hybr', tol=1e-12)
    T = sol.x
    actual = np.trapezoid(p.eps * SIGMA * (T ** 4 - T_SPACE ** 4) - p.q_abs, x)
    ideal = (p.eps * SIGMA * (p.t_root ** 4 - T_SPACE ** 4) - p.q_abs) * L
    return float(actual / ideal), float(T[-1]), bool(sol.success)


def heat_per_panel_area(p: Panel, n_faces=2):
    """Net W per m^2 of PHYSICAL panel, fin efficiency applied."""
    eta, _ = fin_efficiency(p)
    gross = p.eps * SIGMA * (p.t_root ** 4 - T_SPACE ** 4)
    return n_faces * eta * (gross - p.q_abs), eta


def specific_mass(p: Panel, n_faces=2):
    """kg per kW rejected -- the quantity the fin geometry optimizes."""
    q, eta = heat_per_panel_area(p, n_faces)
    return p.areal_mass / q * 1e3, q, eta


def _neumann_d2(n, h):
    """1-D second-difference operator, both ends adiabatic via mirrored ghost node."""
    main = np.full(n, -2.0)
    upper = np.ones(n - 1); lower = np.ones(n - 1)
    upper[0] = 2.0          # ghost T[-1] = T[1]
    lower[-1] = 2.0         # ghost T[n] = T[n-2]
    return sp.diags([lower, main, upper], [-1, 0, 1]) / h ** 2
 
 
def fd_fin_2d_end(p: Panel, overhang_m, run_m=0.25, h=0.002, max_iter=50):
    """
    2-D fin patch with a heat pipe that ENDS before the panel edge.
 
    x in [0, L] across the fin (pipe at x = 0), y in [0, run_m + overhang_m]
    along the pipe. T = T_root on x = 0 for y <= run_m only. Everything else
    adiabatic: x = 0 beyond the pipe end (symmetry line), x = L (between pipes),
    y = 0 (far from the end; run_m = 0.25 m is ~3 conduction lengths and the
    result is unchanged at 0.5 m), y = Ly (panel edge). Second-order ghost-node
    Neumann on all adiabatic edges. Newton with a sparse analytic Jacobian;
    raises rather than returning an unconverged field.
 
    Returns dict:
        Q              net W rejected by the patch, per face
        net_per_y      W/m per face at each y station
        net_per_area   Q / patch area, W/m2 per face
        lost_length_m  end deficit as equivalent pipe length:
                       (q_line * Ly - Q) / q_line, q_line from the piped run
    With overhang_m = 0 the solution is exactly 1-D and lost_length_m = 0.
    """
    L, Ly = p.half_length, run_m + overhang_m
    nx, ny = int(round(L / h)) + 1, int(round(Ly / h)) + 1
    x, y = np.linspace(0, L, nx), np.linspace(0, Ly, ny)
    c = p.eps * SIGMA / (p.k * p.t_face)
    s = p.q_abs / (p.k * p.t_face)
 
    A = (sp.kron(_neumann_d2(nx, x[1] - x[0]), sp.identity(ny))
         + sp.kron(sp.identity(nx), _neumann_d2(ny, y[1] - y[0]))).tocsr()
    pipe = np.zeros((nx, ny), bool)
    pipe[0, y <= run_m + 1e-12] = True
    d = pipe.ravel()
    keep = sp.diags((~d).astype(float))
 
    T = np.full(nx * ny, p.t_root - 15.0); T[d] = p.t_root
    for _ in range(max_iter):
        R = A @ T - c * (T ** 4 - T_SPACE ** 4) + s
        R[d] = T[d] - p.t_root
        J = keep @ (A - sp.diags(4 * c * T ** 3)) + sp.diags(d.astype(float))
        dT = spsolve(J.tocsc(), R)
        T -= dT
        if np.max(np.abs(dT)) < 1e-9:
            break
    else:
        raise RuntimeError(f"fd_fin_2d_end: Newton did not converge in {max_iter} iterations")
 
    T = T.reshape(nx, ny)
    net_per_y = np.trapezoid(p.eps * SIGMA * (T ** 4 - T_SPACE ** 4) - p.q_abs, x, axis=0)
    Q = float(np.trapezoid(net_per_y, y))
    piped = y <= run_m - 3 * L                     # well upstream of the end
    q_line = float(net_per_y[piped].mean()) if piped.any() else float(net_per_y[0])
    return dict(Q=Q, net_per_y=net_per_y, y=y, T=T,
                net_per_area=Q / (L * Ly),
                lost_length_m=(q_line * Ly - Q) / q_line)


def net_per_area_1d(p: Panel):
    """Net W/m^2 of one radiating face, fin efficiency applied. 1-D reference."""
    eta, _ = fin_efficiency(p)
    return eta * (p.eps * SIGMA * (p.t_root ** 4 - T_SPACE ** 4) - p.q_abs)
