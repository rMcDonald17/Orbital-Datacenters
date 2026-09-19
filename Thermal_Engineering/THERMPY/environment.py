"""
environment.py -- orbital thermal environment for the ODRS datacenter study.

View factors, orbit geometry, eclipse, and absorbed environmental flux for an
arbitrary surface on a nadir-locked bus. Physics only: no mass models, no
sizing. Everything here is validated against a closed form or against the
Phase 1 SPENVIS results.

    from environment import Case, net_rejection, FACES

    c = Case(alt_km=800, beta_deg=0.0)
    q, ecl = absorbed_flux(c, 'zenith')
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# ------------------------------------------------------------------ constants

R_EARTH = 6371.0            # km, mean
MU = 398600.4418            # km^3/s^2
SIGMA = 5.670374419e-8      # W m^-2 K^-4
OBLIQUITY = 23.44           # deg

# Gilmore, Spacecraft Thermal Control Handbook -- LEO values
ENV = {
    'cold':    dict(q_sol=1322.0, q_ir=218.0, albedo=0.25),
    'nominal': dict(q_sol=1361.0, q_ir=237.0, albedo=0.30),
    'hot':     dict(q_sol=1414.0, q_ir=258.0, albedo=0.35),
}

# thermo-optical properties. OSR / silvered teflon radiator.
OSR_BOL = dict(alpha_s=0.09, eps_ir=0.85)
OSR_EOL = dict(alpha_s=0.18, eps_ir=0.80)

ALTS_B = [500, 700, 1000, 1200, 1500, 2000]     # Phase 1 B-case matrix
INC_B = 30.0


# ------------------------------------------------------------------ geometry

def earth_angular_radius(alt_km):
    """rho = arcsin(Re/(Re+h)), radians."""
    return np.arcsin(R_EARTH / (R_EARTH + np.asarray(alt_km, float)))


def period(alt_km):
    """Orbital period of a circular orbit, seconds."""
    return 2 * np.pi * np.sqrt((R_EARTH + np.asarray(alt_km, float)) ** 3 / MU)


def beta_max(inc_deg=INC_B):
    """Hard ceiling on |beta| for any RAAN and epoch: i + obliquity."""
    return inc_deg + OBLIQUITY


def beta_critical(alt_km):
    """|beta| above which a circular orbit never enters the cylindrical umbra, deg."""
    return np.degrees(earth_angular_radius(alt_km))


def view_factor(theta_deg, alt_km, n=250):
    """
    View factor from a flat plate to the Earth sphere, by integration over the
    Earth's angular cone.

    theta_deg : angle between the plate normal and the nadir vector.
                0 = nadir-facing, 90 = edge-on, 180 = zenith-facing.

    The sin^2(rho) closed form is only valid at theta=0; every tilted surface
    needs this integral. Validated against sin^2(rho) to 0.001% (see
    validate_view_factor).
    """
    rho = earth_angular_radius(alt_km)
    th = np.radians(theta_deg)
    normal = np.array([np.sin(th), 0.0, np.cos(th)])

    polar, azim = np.meshgrid(np.linspace(0, rho, n),
                              np.linspace(0, 2 * np.pi, 2 * n), indexing='ij')
    direction = np.stack([np.sin(polar) * np.cos(azim),
                          np.sin(polar) * np.sin(azim),
                          np.cos(polar)])
    cos_a = np.einsum('i...,i->...', direction, normal)
    integrand = np.where(cos_a > 0, cos_a * np.sin(polar), 0.0)
    return np.trapezoid(np.trapezoid(integrand, azim[0], axis=1), polar[:, 0]) / np.pi


def validate_view_factor(alts=ALTS_B):
    """Integrator vs sin^2(rho) at theta=0. Returns max relative error."""
    err = []
    for h in alts:
        exact = np.sin(earth_angular_radius(h)) ** 2
        err.append(abs(view_factor(0.0, h) - exact) / exact)
    return max(err)


# ------------------------------------------------------------------ surfaces

# label -> (angle from nadir, normal builder(r_hat, v_hat))
FACES = {
    'nadir':        (0.0,   lambda r, v: -r),
    'zenith':       (180.0, lambda r, v:  r),
    'orbit-normal': (90.0,  lambda r, v:  np.stack([np.zeros_like(r[0]),
                                                    np.zeros_like(r[0]),
                                                    np.ones_like(r[0])])),
    'ram':          (90.0,  lambda r, v:  v),
    'wake':         (90.0,  lambda r, v: -v),
}


@dataclass
class Case:
    """One orbit condition. Attitude is nadir-locked (see limitations)."""
    alt_km: float
    beta_deg: float = 0.0
    env: str = 'nominal'
    optics: dict = field(default_factory=lambda: dict(OSR_BOL))
    n_steps: int = 2000

    @property
    def q_sol(self): return ENV[self.env]['q_sol']

    @property
    def q_ir(self): return ENV[self.env]['q_ir']

    @property
    def albedo(self): return ENV[self.env]['albedo']

    @property
    def period_s(self): return float(period(self.alt_km))


def trajectory(case: Case):
    """
    One orbit in the orbit frame. Sun phase fixed at 0 (worst alignment).

    Returns u, eclipse mask, r_hat, v_hat, sun vector, cos(solar zenith at
    sub-satellite point).
    """
    u = np.linspace(0, 2 * np.pi, case.n_steps, endpoint=False)
    b = np.radians(case.beta_deg)
    sun = np.array([np.cos(b), 0.0, np.sin(b)])

    r_hat = np.stack([np.cos(u), np.sin(u), np.zeros_like(u)])
    v_hat = np.stack([-np.sin(u), np.cos(u), np.zeros_like(u)])

    cos_zenith = np.einsum('i...,i->...', r_hat, sun)
    d_perp = (R_EARTH + case.alt_km) * np.sqrt(np.clip(1 - cos_zenith ** 2, 0, None))
    eclipse = (cos_zenith < 0) & (d_perp < R_EARTH)
    return u, eclipse, r_hat, v_hat, sun, cos_zenith


def eclipse_fraction(case: Case):
    """Fraction of the revolution spent in umbra."""
    return float(trajectory(case)[1].mean())


def mean_eclipse_fraction(alt_km, n_beta=28, **kw):
    """Eclipse fraction averaged over the annual beta sweep."""
    betas = np.linspace(0.0, beta_max(), n_beta)
    return float(np.mean([eclipse_fraction(Case(alt_km, b, **kw)) for b in betas]))


# ------------------------------------------------------------------ flux

def absorbed_flux(case: Case, face: str):
    """
    Absorbed environmental flux on one surface over the orbit, W/m^2.

        q_abs = alpha_s*q_sol*max(0,cos_theta_sun)     direct solar
              + alpha_s*a*q_sol*F*K                    albedo
              + eps*q_ir*F                             Earth IR

    Both solar terms vanish in eclipse. Returns (q_abs array, eclipse mask).
    """
    theta_nadir, build_normal = FACES[face]
    u, eclipse, r_hat, v_hat, sun, cos_zenith = trajectory(case)

    F = view_factor(theta_nadir, case.alt_km)
    normal = build_normal(r_hat, v_hat)
    cos_sun = np.einsum('i...,i->...', normal, sun)

    a_s, eps = case.optics['alpha_s'], case.optics['eps_ir']
    lit = ~eclipse
    K = np.clip(cos_zenith, 0, None) * lit

    solar = a_s * case.q_sol * np.clip(cos_sun, 0, None) * lit
    albedo = a_s * case.albedo * case.q_sol * F * K
    earth_ir = eps * case.q_ir * F
    return solar + albedo + earth_ir, eclipse


def gross_emission(t_rad_k, eps_ir=OSR_BOL['eps_ir']):
    """eps*sigma*T^4, W/m^2."""
    return eps_ir * SIGMA * np.asarray(t_rad_k, float) ** 4


def net_rejection(case: Case, face: str, t_rad_k=318.0):
    """Worst-instant net rejection over the orbit, W/m^2."""
    q_abs, _ = absorbed_flux(case, face)
    return float(gross_emission(t_rad_k, case.optics['eps_ir']) - q_abs.max())


def net_rejection_worst_beta(alt_km, face, t_rad_k=318.0, n_beta=28, **kw):
    """
    Worst net rejection over BOTH orbit position and the full beta sweep.

    A 30 deg shell's node regresses 3.3-6.6 deg/day, so beta cycles through its
    whole range every 47-84 days. A fixed radiator must survive the worst beta,
    not a chosen one.
    """
    betas = np.linspace(0.0, beta_max(), n_beta)
    return min(net_rejection(Case(alt_km, b, **kw), face, t_rad_k) for b in betas)


def best_orientation(alt_km, t_rad_k=318.0, **kw):
    """(face, net_rejection) for the best fixed orientation at worst beta."""
    scores = {f: net_rejection_worst_beta(alt_km, f, t_rad_k, **kw) for f in FACES}
    best = max(scores, key=scores.get)
    return best, scores[best], scores
