"""
orbital.py -- drag, orbit lifetime, station-keeping and disposal for the ODRS
datacenter study.

Drag decay and disposal are the same physics with the sign of the altitude
change reversed, and they share the density model, so they live together.

SCREENING GRADE. Exponential atmosphere with a static solar-activity
multiplier, not NRLMSISE-00 or JB2008, and no phased solar cycle. The shape of
the results -- a cliff below ~600 km, a shallow propulsion minimum near 800 km,
flatness across 700-1200 km -- is robust. Individual year counts are not.

    from orbital import lifetime_years, station_keeping_dv, disposal_dv
"""
from __future__ import annotations

import numpy as np

from environment import R_EARTH, MU

# Vallado exponential atmosphere, mean solar activity: (h0 km, rho0 kg/m3, H km)
_ATMOS = [(450, 1.585e-12, 62.20), (500, 6.967e-13, 71.835),
          (600, 1.454e-13, 88.667), (700, 3.614e-14, 124.64),
          (800, 1.170e-14, 181.05), (900, 5.245e-15, 268.00),
          (1000, 3.019e-15, 268.00)]

# crude multiplier on density by solar activity. Worst case for DECAY is 'max';
# worst case for DISPOSAL is 'min' (thin air, slow decay).
SOLAR = {'min': 0.25, 'mean': 1.0, 'max': 6.0}

CD = 2.2
SEC_PER_YEAR = 3.15576e7
G0 = 9.80665


def density(alt_km, solar='mean'):
    """Exponential atmosphere density, kg/m^3."""
    h0, rho0, H = _ATMOS[0]
    for a, r, s in _ATMOS:
        if alt_km >= a:
            h0, rho0, H = a, r, s
    return rho0 * np.exp(-(alt_km - h0) / H) * SOLAR[solar]


def area_to_mass(area_m2, mass_kg):
    """A/m in m^2/kg. Starlink V2 Mini ~0.03-0.04; ISS ~0.006."""
    return area_m2 / mass_kg


def ballistic_coefficient(area_m2, mass_kg, cd=CD):
    """B = m/(Cd*A), kg/m^2. Higher is better."""
    return mass_kg / (cd * area_m2)


def lifetime_years(alt_km, aom, solar='mean', cd=CD, h_reentry=200.0,
                   max_years=25.0, dt=86400.0):
    """
    Years to decay from a circular orbit to h_reentry. Returns None if the
    vehicle survives past max_years.

        da/dt = -rho * Cd * (A/m) * sqrt(mu*a)
    """
    a = (R_EARTH + alt_km) * 1e3
    re_m, mu_m = R_EARTH * 1e3, MU * 1e9
    t = 0.0
    while t < max_years * SEC_PER_YEAR:
        h = (a - re_m) / 1e3
        if h < h_reentry:
            return t / SEC_PER_YEAR
        a -= density(h, solar) * cd * aom * np.sqrt(mu_m * a) * dt
        t += dt
    return None


def min_altitude_for_lifetime(aom, years=5.0, solar='max', lo=450.0, hi=1400.0):
    """Lowest altitude meeting a passive-lifetime requirement, km."""
    from scipy.optimize import brentq
    f = lambda h: (lifetime_years(h, aom, solar, max_years=years * 8) or years * 8) - years
    try:
        return float(brentq(f, lo, hi, xtol=5.0))
    except ValueError:
        return float('nan')


def station_keeping_dv(alt_km, aom, years=5.0, solar='max', cd=CD):
    """Drag make-up delta-V over the mission, m/s."""
    a = (R_EARTH + alt_km) * 1e3
    v = np.sqrt(MU * 1e9 / a)
    return density(alt_km, solar) * cd * aom * v ** 2 * years * SEC_PER_YEAR


def _hohmann_dv(h1, h2):
    """Two-burn transfer between circular orbits, m/s."""
    r1, r2 = R_EARTH + h1, R_EARTH + h2
    at = 0.5 * (r1 + r2)
    dv1 = abs(np.sqrt(MU / r1) - np.sqrt(MU * (2 / r1 - 1 / at)))
    dv2 = abs(np.sqrt(MU / r2) - np.sqrt(MU * (2 / r2 - 1 / at)))
    return (dv1 + dv2) * 1e3


def _perigee_lowering_dv(h1, h_perigee):
    """Single burn lowering perigee for controlled reentry, m/s."""
    r1, rp = R_EARTH + h1, R_EARTH + h_perigee
    at = 0.5 * (r1 + rp)
    return (np.sqrt(MU / r1) - np.sqrt(MU * (2 / r1 - 1 / at))) * 1e3


def disposal_dv(alt_km, mode='drag_assist', h_target=550.0, h_perigee=60.0):
    """
    Post-mission disposal delta-V, m/s.

    FCC 22-74 (47 CFR 25.283) requires disposal within 5 years of mission end.
    The same order caps reentry casualty risk at 1 in 10,000, which a multi-
    tonne vehicle may not meet by uncontrolled reentry -- in which case
    mode='controlled' applies and chemical propulsion is needed, since electric
    thrust cannot target an impact footprint.

    drag_assist : lower to h_target and let drag finish the job
    controlled  : lower perigee to h_perigee for a targeted reentry
    """
    if mode == 'drag_assist':
        return _hohmann_dv(alt_km, h_target)
    if mode == 'controlled':
        return _perigee_lowering_dv(alt_km, h_perigee)
    raise ValueError(f"unknown disposal mode: {mode}")


def propellant_mass(dry_mass_kg, dv_ms, isp_s):
    """Rocket equation. isp ~1500 s electric, ~220 s monoprop chemical."""
    return dry_mass_kg * (1 - np.exp(-dv_ms / (isp_s * G0)))


def propulsion_budget(alt_km, aom, dry_mass_kg, years=5.0, isp_s=1500.0,
                      disposal_mode='drag_assist'):
    """Combined station-keeping + disposal budget. Returns a dict."""
    sk = station_keeping_dv(alt_km, aom, years)
    dp = disposal_dv(alt_km, disposal_mode)
    total = sk + dp
    m = propellant_mass(dry_mass_kg, total, isp_s)
    return dict(alt_km=alt_km, sk_dv=sk, disposal_dv=dp, total_dv=total,
                propellant_kg=m, frac_dry=m / dry_mass_kg)
