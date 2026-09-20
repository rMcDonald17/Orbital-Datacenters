# THERMPY

Thermal and orbital analysis for the 30° shell of the orbital datacenter study.
Companion to `SPENVIS/SPENPY/`, same conventions: modules named for what they do,
no work at import time, one driver per deliverable.

## Modules

| Module | Contents | Confidence |
|---|---|---|
| `environment.py` | View factors, orbit geometry, eclipse, absorbed flux, net rejection | **Validated** — integrator vs `sin²ρ` to 8e-6; eclipse vs Phase 1 SPENVIS |
| `orbital.py` | Atmosphere, decay, station-keeping, disposal ΔV, propellant | **Screening** — exponential atmosphere, static solar multiplier |
| `sizing.py` | Radiator/array/battery area and mass; vehicle rollup | **Assumption** — every mass model unsourced |
| `conduction.py` | Fin ODE, efficiency, panel mass optimisation | **Validated** — collocation vs Newton FD to 0.05% |
| `run_altitude_selection_30deg.py` | Driver → `docs/altitude_selection_30deg.md` | — |
| `run_radiator_fin_submodel.py` | Driver → `docs/radiator_fin_submodel.md` | — |

Drivers are named for the document they generate, one each. They are flat scripts
with `# %%` cells, not importable modules — nothing else depends on them.

The split between `environment.py` and `sizing.py` is deliberate. The first computes
W/m² and is validated; the second converts W/m² into m² and kg using assumptions.
Mass models can be replaced without touching validated physics.

## Use

```bash
python run_altitude_selection_30deg.py     # altitude trade, drag, disposal
python run_radiator_fin_submodel.py       # fin ODE, validation, mass optimisation
```

```python
import environment as env, sizing as sz, orbital as orb

face, net, all_faces = env.best_orientation(800)      # -> ('zenith', 370.4, {...})
v = sz.size_vehicle(payload_kw=40, alt_km=800, net_rejection_wm2=net)
aom = orb.area_to_mass(sz.drag_area(v), v['total_kg'])
orb.propulsion_budget(800, aom, v['total_kg'])
```

## Assumptions to replace before publication

Collected in `sizing.MassModel`. In order of leverage:

1. `radiator_kg_per_m2 = 10.0` — largest single mass item and a driver of A/m
2. `array_w_per_kg = 150.0` — range 100–200 in the literature
3. `array_w_per_m2 = 312.0` — sets drag area, currently a guess
4. `rack_kg_per_kw = 11.3` — extrapolated from one NVL72 datapoint

`orbital.SOLAR` multipliers are the crudest thing in the package. Replace with
NRLMSISE-00 or JB2008 with a phased solar cycle before any published lifetime claim.
