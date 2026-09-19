"""
run_part_a.py -- Part A: altitude trade for the 30 deg shell.

Reproduces every table in docs/altitude_selection_30deg.md.
Run top to bottom, or cell by cell.

Physics lives in environment.py / orbital.py / sizing.py. This is just the driver.
"""
# %% ---- setup
import numpy as np

import environment as env
import orbital as orb
import sizing as sz

PAYLOAD_KW = 40.0
T_RAD = 318.0           # held fixed for the trade; selected in Part D
SELECTED_ALT = 800.0    # not a Phase 1 B-case -- see write-up limitations

rule = lambda t: print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)

# %% ---- validation
rule("VALIDATION")
print(f"view factor integrator vs sin^2(rho) at theta=0: "
      f"max rel err {env.validate_view_factor():.2e}")

f500 = env.eclipse_fraction(env.Case(500, 0.0))
print(f"eclipse at 500 km, beta=0: {f500*100:.2f}%  "
      f"({f500*env.Case(500).period_s/60:.1f} min)   [Phase 1: 37.8%, 35.8 min]")

# %% ---- view factors
rule("VIEW FACTORS TO EARTH (nadir-locked)")
print(f"{'h (km)':>7s}{'nadir':>10s}{'zenith':>10s}{'edge-on':>10s}")
for h in env.ALTS_B:
    print(f"{h:7d}{env.view_factor(0, h):10.4f}"
          f"{env.view_factor(180, h):10.4f}{env.view_factor(90, h):10.4f}")

# %% ---- eclipse
rule("ECLIPSE")
print(f"{'h':>6s}{'period':>9s}{'f(b=0)':>9s}{'dur min':>9s}{'f(bmax)':>9s}")
for h in env.ALTS_B:
    T = env.Case(h).period_s / 60
    f0 = env.eclipse_fraction(env.Case(h, 0.0))
    fm = env.eclipse_fraction(env.Case(h, env.beta_max()))
    print(f"{h:6d}{T:9.1f}{f0:9.3f}{f0*T:9.1f}{fm:9.3f}")

# %% ---- radiator orientation
# A 30 deg node regresses 3.3-6.6 deg/day, so beta cycles through its full range
# every 47-84 days. A fixed radiator must survive the worst beta, not a chosen one.
rule(f"WORST-BETA NET REJECTION (W/m2) at Trad={T_RAD:.0f} K")
print(f"{'h (km)':>7s}" + "".join(f"{f:>14s}" for f in env.FACES))

nets = {}
for h in env.ALTS_B:
    best, val, scores = env.best_orientation(h, T_RAD)
    nets[h] = val
    print(f"{h:7d}" + "".join(f"{scores[f]:14.0f}" for f in env.FACES)
          + f"   <- {best}")

# %% ---- thermal figure of merit
rule(f"THERMAL FIGURE OF MERIT ({PAYLOAD_KW:.0f} kW) -- assumption-grade mass models")
print(f"{'h':>6s}{'Phi':>7s}{'A_rad':>8s}{'m_rad':>8s}{'m_arr':>8s}"
      f"{'m_bat':>8s}{'therm':>8s}{'kg/kW':>8s}{'cyc/yr':>8s}")
for h in env.ALTS_B:
    v = sz.size_vehicle(PAYLOAD_KW, h, nets[h])
    print(f"{h:6d}{v['phi']:7.2f}{v['a_rad_m2']:8.0f}{v['m_radiator']:8.0f}"
          f"{v['m_array']:8.0f}{v['m_battery']:8.0f}{v['thermal_kg']:8.0f}"
          f"{v['thermal_kg_per_kw']:8.1f}{v['cycles_per_year']:8.0f}")

# %% ---- the selected vehicle, sized self-consistently
net_sel = env.best_orientation(SELECTED_ALT, T_RAD)[1]
veh = sz.size_vehicle(PAYLOAD_KW, SELECTED_ALT, net_sel)
aom = orb.area_to_mass(sz.drag_area(veh), veh['total_kg'])

rule(f"SELECTED VEHICLE -- {SELECTED_ALT:.0f} km, {PAYLOAD_KW:.0f} kW")
print(f"  mass {veh['total_kg']:.0f} kg | radiating {veh['a_rad_m2']:.0f} m2 "
      f"| array {veh['a_array_m2']:.0f} m2 | battery {veh['battery_kwh']:.0f} kWh")
for att in ('broadside', 'tumble', 'edge'):
    A = sz.drag_area(veh, attitude=att)
    print(f"  {att:10s} A={A:6.1f} m2  A/m={orb.area_to_mass(A, veh['total_kg']):.4f}"
          f"  B={orb.ballistic_coefficient(A, veh['total_kg']):6.1f} kg/m2")
print("  reference: Starlink V2 Mini ~0.03-0.04, ISS ~0.006 m2/kg")

# %% ---- passive lifetime
rule(f"PASSIVE LIFETIME (yr to 200 km), A/m={aom:.4f}")
print(f"{'h':>6s}{'solar min':>12s}{'mean':>12s}{'max':>12s}")
for h in [500, 600, 700, 800, 1000, 1200]:
    row = f"{h:6d}"
    for s in ('min', 'mean', 'max'):
        y = orb.lifetime_years(h, aom, s)
        row += f"{'>25' if y is None else f'{y:.2f}':>12s}"
    print(row)
print(f"\n  min altitude for 5-yr passive life, solar max: "
      f"{orb.min_altitude_for_lifetime(aom):.0f} km")

# %% ---- propulsion budget
rule("PROPULSION BUDGET, 5 yr (Isp 1500 s, drag-assisted disposal)")
print(f"{'h':>6s}{'SK dV':>9s}{'disp dV':>10s}{'total':>9s}{'prop kg':>10s}{'% dry':>8s}")
for h in [700, 800, 1000, 1200]:
    b = orb.propulsion_budget(h, aom, veh['total_kg'])
    print(f"{h:6d}{b['sk_dv']:9.0f}{b['disposal_dv']:10.0f}{b['total_dv']:9.0f}"
          f"{b['propellant_kg']:10.0f}{b['frac_dry']*100:8.1f}")

# %% ---- disposal mode sensitivity
# FCC 22-74 caps reentry casualty risk at 1 in 10,000. If a multi-tonne vehicle
# cannot meet that by uncontrolled reentry, controlled reentry is required --
# chemical, and the altitude gradient reverses to favour lower altitudes.
rule("IF CONTROLLED REENTRY IS REQUIRED (chemical, Isp 220 s)")
print(f"{'h':>6s}{'disp dV':>10s}{'prop kg':>10s}{'% dry':>8s}")
for h in [700, 800, 1000, 1200]:
    dv = orb.disposal_dv(h, 'controlled')
    m = orb.propellant_mass(veh['total_kg'], dv, 220.0)
    print(f"{h:6d}{dv:10.0f}{m:10.0f}{m/veh['total_kg']*100:8.1f}")

# %%
