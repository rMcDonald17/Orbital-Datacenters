"""
sizing.py -- radiator, array and battery sizing and mass models.

Deliberately separated from environment.py. That module computes W/m^2 and is
validated against closed forms; this one converts W/m^2 into m^2 and kg using
mass models that are currently UNSOURCED ASSUMPTIONS. The two have very
different confidence levels and should be replaceable independently.

Every default in MassModel is an assumption pending a sourced value. The
radiator areal mass dominates both the thermal figure of merit and the A/m that
drives orbit lifetime, so it deserves the most scrutiny.

    from sizing import MassModel, size_vehicle
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from environment import Case, mean_eclipse_fraction, eclipse_fraction, period


@dataclass
class MassModel:
    """All ASSUMPTION grade. Replace with sourced values before publication."""
    radiator_kg_per_m2: float = 10.0      # per m^2 of PHYSICAL panel (5 optimistic, 10-15 w/ fluid+structure)
    array_w_per_kg: float = 150.0         # array level, flexible blanket 100-200
    cell_wh_per_kg: float = 150.0         # cell level
    depth_of_discharge: float = 0.70
    rack_kg_per_kw: float = 11.3          # from the NVL72 datapoint, 1.36 t / 120 kW
    array_w_per_m2: float = 312.0         # ~30% cells, packing 0.9, EOL 0.85
    structure_fraction: float = 0.25      # of subsystem subtotal
    eta_charge: float = 0.85
    eta_pmad: float = 0.90
    area_margin: float = 1.25             # on radiator area


def bus_power_kw(payload_kw, optical_kw=1.0, avionics_kw=0.8, pump_frac=0.015):
    """Total heat to reject. Essentially all electrical input becomes heat."""
    return payload_kw * (1 + pump_frac) + optical_kw + avionics_kw


def array_factor(alt_km, mm: MassModel, n_beta=28, **kw):
    """
    Phi = 1 + f_ecl / ((1 - f_ecl) * eta_charge)

    Ratio of array power to load power: the array carries the load AND
    recharges the battery during the sunlit fraction. Uses the beta-averaged
    eclipse fraction, not the beta=0 value.
    """
    f = mean_eclipse_fraction(alt_km, n_beta, **kw)
    return 1 + f / ((1 - f) * mm.eta_charge), f


def radiator_area(q_total_kw, net_rejection_wm2, mm: MassModel):
    """Radiating area (m^2) including margin, and physical panel area."""
    a_rad = q_total_kw * 1e3 / net_rejection_wm2 * mm.area_margin
    return a_rad, a_rad / 2.0          # double-sided deployed panel


def battery_mass(q_total_kw, alt_km, mm: MassModel, beta_deg=0.0):
    """Sized to the worst-case eclipse pass (beta=0), not the average."""
    f = eclipse_fraction(Case(alt_km, beta_deg))
    energy_kwh = q_total_kw * f * float(period(alt_km)) / 3600.0
    return energy_kwh * 1e3 / mm.cell_wh_per_kg / mm.depth_of_discharge, energy_kwh


def size_vehicle(payload_kw, alt_km, net_rejection_wm2, mm: MassModel = None, **kw):
    """
    Full vehicle sizing at one altitude and power. Returns a dict.

    net_rejection_wm2 comes from environment.best_orientation() -- passed in
    rather than computed here so the physics and the mass models stay decoupled.
    """
    mm = mm or MassModel()
    q = bus_power_kw(payload_kw)

    a_rad, a_phys = radiator_area(q, net_rejection_wm2, mm)
    m_rad = a_phys * mm.radiator_kg_per_m2

    phi, f_bar = array_factor(alt_km, mm, **kw)
    p_array = q * phi / mm.eta_pmad
    m_arr = p_array * 1e3 / mm.array_w_per_kg
    a_arr = p_array * 1e3 / mm.array_w_per_m2

    m_bat, e_kwh = battery_mass(q, alt_km, mm)
    m_rack = payload_kw * mm.rack_kg_per_kw

    thermal_subtotal = m_rad + m_arr + m_bat
    subtotal = thermal_subtotal + m_rack
    total = subtotal * (1 + mm.structure_fraction)

    return dict(
        payload_kw=payload_kw, alt_km=alt_km, q_total_kw=q,
        net_rejection=net_rejection_wm2, phi=phi, f_ecl_mean=f_bar,
        a_rad_m2=a_rad, a_panel_m2=a_phys, a_array_m2=a_arr,
        m_radiator=m_rad, m_array=m_arr, m_battery=m_bat, m_rack=m_rack,
        battery_kwh=e_kwh,
        thermal_kg=thermal_subtotal, thermal_kg_per_kw=thermal_subtotal / payload_kw,
        total_kg=total, kg_per_kw=total / payload_kw,
        cycles_per_year=365 * 86400 / float(period(alt_km)),
    )


def drag_area(v: dict, bus_m2=12.0, attitude='tumble'):
    """
    Projected area for drag, from a size_vehicle() result.

    broadside : all surfaces normal to velocity (worst)
    tumble    : half of each flat plate (screening average)
    edge      : 15% of plate area (best, requires active attitude control)
    """
    plates = v['a_array_m2'] + v['a_panel_m2']
    factor = {'broadside': 1.0, 'tumble': 0.5, 'edge': 0.15}[attitude]
    return plates * factor + bus_m2
