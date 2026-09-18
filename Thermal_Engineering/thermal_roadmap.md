# Orbital Datacenter Feasibility Study — Thermal Aside
## Thermal Management vs. Altitude and Power for the 30° Shell

**Author:** rMcDonald17 · **Status:** Aside to Phase 1 (thermal environment, altitude trade, node model)
**Reference scenario:** SpaceX FCC filing, Jan 30, 2026 — orbital data center constellation, LEO, ~30° inclination shells (demand-peak capacity)
**Scope restriction:** 30° inclination only. Sun-synchronous shells are out of scope for this aside.

---

## 1. Purpose

Phase 1 answers the radiation cost of altitude. This aside answers the *thermal* cost of altitude for the same filed range, then builds an analytical thermal model of a single satellite at the altitude that trade selects.

**Design questions this aside answers:**

1. Within the filed 500–2,000 km range, which altitude minimizes the thermal burden on a 30° shell, accounting for Earth IR, albedo, direct solar, and eclipse?
2. At that altitude, what radiator and solar array area, orientation, and mass does a 40 kW compute payload require, and how does that scale to 600 kW?
3. What are the bounding hot and cold cases, and does the vehicle survive both?

**Why 30° only.** A 30° orbit has $|\beta| \le i + 23.44^\circ = 53.4^\circ$, and the no-eclipse threshold

$$
\beta^* = \arcsin\!\left( \frac{R_\oplus}{R_\oplus + h} \right)
$$

does not drop below $53.4^\circ$ until $h \approx 1{,}565$ km. Every altitude in the filed range from 500 km to ~1,565 km therefore eclipses year-round, and above that only briefly at peak $\beta$. The sunlight argument that favors high altitude for the sun-synchronous shells does not apply here, so the altitude trade for this shell is a clean contest between radiation (prefers low) and thermal environment (to be determined by Part A).

---

## 2. Fixed Inputs

### 2.1 Reference scenario (from Phase 1, unchanged)

| Parameter | Value | Source basis |
|---|---|---|
| Altitude range | 500–2,000 km | FCC filing |
| Inclination | ~30° | FCC filing |
| Bus heritage | Starlink V3-class flat-panel, scaled | Musk statements / filing context |
| Launch vehicle | Starship, ~100 t to LEO, ~9 m fairing | Public Starship specs |
| Reference compute payload | H100 SXM-class rack, ~40 kW | Analyst choice (consistent with Phase 1) |
| Power sweep $\dot{Q}_{\text{payload}}$ | 40 / 120 / 190 / 370 / 600 kW | H100 rack → GB200 → Vera Rubin → Rubin CPX → Rubin Ultra rack-class power (public roadmap figures) |
| Replacement cadence | ~5 years | Analyst assumption |

### 2.2 Environmental constants (Gilmore, *Spacecraft Thermal Control Handbook*, LEO values)

| Quantity | Cold case | Nominal | Hot case |
|---|---|---|---|
| Solar constant $q_{\text{sol}}$ (W/m²) | 1,322 (aphelion) | 1,361 | 1,414 (perihelion) |
| Earth IR $q_{\text{IR}}$, top of atmosphere (W/m²) | 218 | 237 | 258 |
| Albedo $a$ | 0.25 | 0.30 | 0.35 |
| Space sink temperature $T_{\text{space}}$ (K) | 3 | 3 | 3 |

### 2.3 Thermo-optical properties

| Surface | $\alpha_s$ (BOL) | $\alpha_s$ (EOL) | $\varepsilon_{\text{IR}}$ | Notes |
|---|---|---|---|---|
| Radiator (OSR / silvered Teflon) | 0.08–0.10 | 0.15–0.20 | 0.80–0.85 | EOL degradation from UV + contamination; EOL used in hot case |
| Solar array (front, cell side) | ~0.90 | ~0.92 | ~0.85 | Absorbed minus electrical output is the heat load |
| Solar array (back, substrate) | — | — | 0.80–0.85 | Radiates to space; sets array temperature |
| Bus MLI (external) | 0.3–0.5 | — | $\varepsilon^* \approx 0.01\text{–}0.03$ | Effective emittance, not surface emittance |

### 2.4 Payload thermal constraints (public datasheets)

| Component | Allowable range | Driver |
|---|---|---|
| H100 SXM junction | $T_j \le 85\,^\circ\mathrm{C}$ | Hot case: sets coolant outlet ceiling |
| HBM | $T \le 95\,^\circ\mathrm{C}$ (typ.) | Confirm against datasheet |
| Li-ion battery, operating | $0$ to $+40\,^\circ\mathrm{C}$ | **Cold case driver** — narrowest allowable range on the vehicle |
| Li-ion battery, charge | $+5$ to $+35\,^\circ\mathrm{C}$ | Constrains post-eclipse recharge |
| Coolant (ammonia) | $T_{\text{freeze}} = -77\,^\circ\mathrm{C}$ | Cold-case check on deployed lines |
| Optical terminal | vendor TBD | Placeholder $\pm 20\,^\circ\mathrm{C}$, stability-limited |

---

## 3. Part A — Altitude Trade for Thermal Management

**Objective:** produce a thermal figure of merit vs. altitude for the 30° shell and select one altitude for Parts B–E.

### 3.1 Altitude cases

Use the Phase 1 B-case altitudes so thermal results align with the radiation database: **500 / 700 / 1,000 / 1,200 / 1,500 / 2,000 km**. Add 600 km if the trade shows a gradient worth resolving between B1 and B2.

### 3.2 Orbit geometry per altitude

For each altitude compute, over the annual $\beta$ sweep ($-53.4^\circ$ to $+53.4^\circ$):

- Orbital period, Earth angular radius $\rho = \arcsin\!\left( R_\oplus / (R_\oplus + h) \right)$
- Eclipse fraction $f_{\text{ecl}}$ and duration vs. $\beta$ (cylindrical shadow; penumbra noted but neglected)
- Annual sunlit fraction
- Worst-case ($\beta = 0$) and best-case ($|\beta| = 53.4^\circ$) eclipse durations

*Scoping note:* eclipse **duration** at $\beta = 0$ is nearly flat (~35 min) from 500 to 1,200 km because period lengthens as the shadow fraction shrinks. Expect battery energy to be weakly altitude-dependent; expect the trade to be driven by Earth IR and albedo view factors instead.

### 3.3 Environmental flux per altitude

For a set of reference surface orientations (nadir, zenith, orbit-normal, ram, wake), compute per altitude and per orbit position:

- View factor to Earth, $F(\theta, \rho)$, by numerical integration over the Earth disk (do not use the $\sin^2\!\rho$ shortcut except for $\theta = 0$ validation)
- Earth IR incident: $q_{\text{IR}} F$
- Albedo incident: $a\, q_{\text{sol}} F\, K(\zeta)$, with $K$ the solar-zenith correction at the sub-satellite point
- Direct solar: $q_{\text{sol}} \max(0, \cos\theta_{\text{sun}})$, zero in eclipse

### 3.4 Figure of merit

Two metrics, reported separately, then combined:

1. **Net radiator rejection, worst instant** (W/m²): $\varepsilon\sigma T^4$ minus absorbed environmental flux, evaluated at the orbit position that minimizes it, for the best available radiator orientation at that altitude. Higher is better.
2. **Array energy-balance factor**: required array power divided by load power,

$$
\Phi = 1 + \frac{f_{\text{ecl}}}{\left( 1 - f_{\text{ecl}} \right) \eta_{\text{charge}}}
$$

   Lower is better.

Combined: **thermal-driven specific mass (kg/kW)** = radiator mass + array mass + battery mass per kW of compute, using the mass models in §7. This is the number that goes on the altitude trade chart alongside the Phase 1 TID curve.

### 3.5 Deliverable T1 — Altitude trade

- Table: period, eclipse, $F_{\text{earth}}$ by orientation, $q_{\text{IR}}$, $q_{\text{alb}}$, net rejection, energy-balance factor, thermal kg/kW, for all six altitudes at $\beta = 0$ and $|\beta| = 53.4^\circ$
- Chart: thermal kg/kW vs. altitude, overlaid on the Phase 1 TID/year curve for the 30° shell (same x-axis) — *this is the aside's money plot*
- One-paragraph selection rationale

**Decision gate:** altitude selection is reviewed and approved before Part B begins.

---

## 4. Part B — Bounding Cases at the Selected Altitude

Define once, apply to every node-model run.

| | **Hot case** | **Cold case** |
|---|---|---|
| Solar constant $q_{\text{sol}}$ | 1,414 W/m² | 1,322 W/m² |
| Earth IR $q_{\text{IR}}$ | 258 W/m² | 218 W/m² |
| Albedo $a$ | 0.35 | 0.25 |
| $\beta$ | Value maximizing solar on radiator + array (from Part A sweep) | $0^\circ$ (longest eclipse) |
| Optical properties | EOL ($\alpha_s$ degraded, $\varepsilon$ at low bound) | BOL |
| Payload state | Full power, all racks | **Safe mode / payload off** through eclipse |
| Bus loads | Max | Min (survival) |
| Array | EOL, degraded per Phase 1 DDD | BOL, cold, high-voltage check |
| Sizing question answered | Radiator area | Survival heater power, battery temperature, coolant freeze margin |

The cold case is not "eclipse with payload on" — the payload's own dissipation keeps the loop warm in that case. The bounding cold case is a payload outage or safe-mode entry during eclipse, which is exactly when a fault is most likely to occur and when heater power competes with a discharging battery.

Add a **nominal case** ($\beta = 0$, nominal constants, BOL, full power) for the parametric sweep and for validation against the algebraic sizing.

---

## 5. Part C — Analytical Thermal Model

Built up in four layers, each validating the next. Terminology follows Thermal Desktop / SINDA conventions so the model maps directly onto that tool if the study later moves there.

### 5.1 Layer 1 — Steady-state energy balance (algebraic)

Whole-vehicle balance, one equation per surface class:

$$
\dot{Q}_{\text{payload}} + \dot{Q}_{\text{bus}}
+ \sum_i \alpha_{s,i} A_i \left( q_{\text{sol}} \cos\theta_i + a\, q_{\text{sol}} F_i K \right)
+ \sum_i \varepsilon_i A_i q_{\text{IR}} F_i
= \sum_i \varepsilon_i A_i \sigma T_i^4
$$

Solve for radiator area at fixed $T_{\text{rad}}$. This layer alone produces the first radiator-area-vs-power curve and is the cross-check for everything below. Reproduce the Phase 1 scoping results ($40$ kW gives $A_{\text{rad}} \approx 117$ m² at $T_{\text{rad}} = 318$ K, deployed panel, 1,200 km) as the regression test, then re-run at the Part A altitude.

### 5.2 Layer 2 — Lumped-parameter node network

Thermal-Desktop-style network. Each node has capacitance $C_i = m_i c_{p,i}$; conductors link nodes; environmental loads apply per node from the Part A flux model.

**Node list (40 kW baseline):**

| Node | Type | Contents | Notes |
|---|---|---|---|
| N1 | Diffusion | GPU die + package + cold plate (lumped) | Lumping avoids stiffness; Biot check justifies it |
| N2 | Diffusion | Coolant loop inventory | Single-phase pumped ammonia; $\Delta T \approx 10$ K design |
| N3 | Diffusion | Radiator root / manifold | Linear conductor to N2 (convective, $hA$) |
| N4 | Diffusion | Radiator panel (fin) | Either one node with fin efficiency, or resolved 1-D (Layer 3) |
| N5 | Diffusion | Bus structure + MLI-wrapped electronics | Low-power, MLI-isolated |
| N6 | Diffusion | Solar array | Radiatively coupled only; own energy balance |
| N7 | Diffusion | Battery | **Narrowest allowable — tracked explicitly** |
| N8 | Boundary | Deep space | 3 K, infinite capacitance |
| — | Heat load | Earth IR, albedo, solar | Time-varying per node from Part A |

**Conductor types:**

- Linear: $G_{ij} = kA/L$ (conduction) or $hA$ (convective, loop-to-wall). Units W/K.
- Radiation: $R_{ij} = \sigma \varepsilon_i A_i F_{ij}$ — "radks" in SINDA; space is the dominant sink for N4, N6. **Note the convention:** $\sigma$ is absorbed into $R_{ij}$, so radiation conductors carry units of $\mathrm{W/K^4}$ while linear conductors carry $\mathrm{W/K}$.
- Contact: $G = h_c A$ at bolted/bonded interfaces (cold plate to package)

**Governing ODE per node:**

$$
C_i \frac{dT_i}{dt} = \dot{Q}_i(t)
+ \sum_j G_{ij} \left( T_j - T_i \right)
+ \sum_j R_{ij} \left( T_j^4 - T_i^4 \right)
$$

where the node source term separates constant internal dissipation from the time-varying environment:

$$
\dot{Q}_i(t) = \dot{q}_{\text{int},i}
+ \alpha_{s,i} A_i \left[ q_{\text{sol}} \cos\theta_i(t) + a\, q_{\text{sol}} F_i(t) K(t) \right]
+ \varepsilon_i A_i q_{\text{IR}} F_i(t)
$$

with $\cos\theta_i$ clamped at zero and both solar terms vanishing in eclipse. Setting $dT_i/dt = 0$ and summing over $i$ recovers the Layer 1 balance of §5.1.

Integrate over $\ge 3$ orbits to a periodic limit cycle; report the last orbit. Use a stiff-capable integrator (Radau / BDF) if N1 is split into die and plate.

**Deliverable T2:** node network diagram + conductor table + capacitance table, as a standalone `nodes.md` that a Thermal Desktop user could rebuild from.

### 5.3 Layer 3 — Fourier conduction sub-models

Analytical 1-D solutions used to derive the lumped parameters in Layer 2, and to check them:

- **Radiator fin:** rectangular fin radiating to space, root at $T_{\text{root}}$. The governing equation

$$
 k t \frac{d^2 T}{dx^2} = \varepsilon \sigma \left( T^4 - T_{\text{space}}^4 \right)
$$

  has no closed form, so use (a) the linearized fin solution with radiative film coefficient $h_{\text{rad}} = 4\varepsilon\sigma \bar{T}^3$ for the efficiency estimate, and (b) a shooting or finite-difference solution of the nonlinear form as the check. Output: fin efficiency $\eta_{\text{fin}}$ vs. fin length and thickness, feeding N4 and the radiator mass model.
- **Cold plate spreading:** point/area heat source on a finite plate — spreading resistance per Lee et al. or Yovanovich. Output: die-to-coolant resistance $R_{\text{sp}}$, giving conductor $G_{1,2}$.
- **Heat pipe / loop conductance:** effective conductance of the transport element between root and fin, giving conductor $G_{3,4}$.
- **MLI:** effective emittance $\varepsilon^*$ model for the bus wrap, giving radiation conductor $R_{5,8}$.

**Deliverable T3:** `conduction_models.md` with each derivation, its assumptions, and the resulting conductor/efficiency value used in Layer 2.

### 5.4 Layer 4 — Discretized cross-check

Finite-difference (or finite-element, if tooling is available) solution of the radiator fin and cold plate, compared against the Layer 3 analytical results. Purpose: confirm that the analytical fin efficiency and spreading resistance are within tolerance of a discretized solve before trusting them in the sweep. Document the discrepancy and its cause, as with the Phase 1 transport-vs-SHIELDOSE validation.

*Author to confirm the intended scope of "FTE analysis" — interpreted here as this discretized finite-element/finite-difference thermal cross-check. Adjust if a different method was meant.*

### 5.5 Power sweep

Run Layers 1–2 at 40 / 120 / 190 / 370 / 600 kW, hot / cold / nominal. Layer 3 parameters are held fixed unless fin geometry changes with scale. Report:

- Radiator area, physical panel dimensions, mass
- Array area and mass
- Battery capacity and mass
- Node temperatures over the limit-cycle orbit, hot and cold
- Margin against every allowable in §2.4

---

## 6. Part D — Radiator Sizing and Orientation

### 6.1 Orientation trade

Evaluate at the Part A altitude, $\beta = 0$ and $\beta = 53.4^\circ$, hot case:

| Configuration | Sun exposure | Earth view | Deployment |
|---|---|---|---|
| Body-mounted, zenith face | Direct at noon | $F \approx 0$ | None |
| Body-mounted, nadir face | None | $F = \sin^2\!\rho$ (max) | None |
| Body-mounted, ram/wake faces | Once per orbit | Partial | None |
| Deployed, edge-on in orbit plane | Grazing at $\beta = 0$; rises with $\beta$ | $F(\theta = 90^\circ) \approx 0.17$ at 1,200 km — *not zero* | Single or multi-fold |
| Deployed, sun-tracking edge-on | Grazing always | Varies | Adds articulation |

Size to the **worst instant** in the orbit, not the orbit average. Report both so the penalty of the conservative choice is visible.

### 6.2 Radiator temperature selection

Sweep $T_{\text{rad}} = 290 / 300 / 310 / 318 / 333 / 350$ K. Radiator area scales as

$$
A_{\text{rad}} = \frac{\dot{Q}_{\text{total}}}{\varepsilon \sigma T_{\text{rad}}^4 - q_{\text{abs}}}
$$

so report area and mass at each point. Tie $T_{\text{rad}}$ to the coolant chain: $T_{\text{rad}} \approx T_{\text{coolant,in}} - \Delta T_{\text{root}}$, with $T_{\text{coolant,out}}$ bounded by the junction limit. State the recommended operating point and the leakage/reliability trade it implies for the silicon.

### 6.3 Geometry

Double-sided deployed panels: $A_{\text{rad}} = 2 A_{\text{phys}}$. Check physical dimensions against the ~9 m fairing envelope; note where multi-fold deployment becomes necessary in the power sweep.

### 6.4 Mass model

Deployed radiator assembly: areal mass $\mu_{\text{rad}}$ in kg per m² of **physical panel**, range 5 (panel only, optimistic) to 10–15 (with fluid, manifolds, deployment structure). Source flight data (ISS PVR/HRS radiators, commercial deployables) and record the basis. The radiator is expected to be the largest single mass item; this number deserves the most scrutiny in the aside.

---

## 7. Part E — Solar Array Sizing and Pointing

### 7.1 Energy balance

$$
P_{\text{array}} = \frac{\dot{Q}_{\text{payload}} + \dot{Q}_{\text{bus}}}{\eta_{\text{PMAD}}}
\left[ 1 + \frac{f_{\text{ecl}}}{\left( 1 - f_{\text{ecl}} \right) \eta_{\text{charge}}} \right]
$$

Size at EOL using the Phase 1 DDD result for cell degradation (cross-link to `environment/` D1 database and the MC-SCREAM / EQFLUX runs already in the SPENVIS workflow).

### 7.2 Pointing

At $\beta = 0$ the sun lies in the orbit plane, so a fixed array on a nadir-locked bus sees $\cos$-varying incidence around the orbit. Options:

- Single-axis tracking about orbit-normal (removes the in-plane cos loss; standard for LEO)
- Second axis or seasonal offset for the $\beta$ sweep ($\pm 53.4^\circ$) — or accept the $\cos\beta$ loss and oversize
- Fixed, oversized array (simplest, heaviest)

Report array area and mass for each; select one for the sweep.

### 7.3 Array thermal node

N6 energy balance, absorbed solar minus electrical output, radiating from both faces:

$$
\left[ \alpha_s - \eta_{\text{cell}}(T) \right] q_{\text{sol}} \cos\theta
= \left( \varepsilon_f + \varepsilon_b \right) \sigma T^4 + \dot{q}_{\text{cond,bus}}
$$

Include the cell efficiency temperature coefficient $d\eta_{\text{cell}}/dT$ so array output and temperature are solved together. Cold case: post-eclipse open-circuit voltage at minimum temperature.

### 7.4 Mass model

Specific power (W/kg) at array level, flexible-blanket vs. rigid, with source. Range 100–200 W/kg; record basis.

---

## 8. Deliverables

| ID | Deliverable | Part |
|---|---|---|
| T1 | Altitude trade table + thermal kg/kW vs. altitude chart overlaid on TID curve | A |
| T2 | Node network diagram, conductor and capacitance tables (`nodes.md`) | C |
| T3 | Fourier sub-model derivations (`conduction_models.md`) | C |
| T4 | Discretized cross-check vs. analytical, with discrepancy discussion | C |
| T5 | Hot/cold/nominal case definitions (`cases.md`) | B |
| T6 | Radiator orientation and temperature trade | D |
| T7 | Array sizing and pointing trade | E |
| T8 | Power sweep results: area, mass, temperatures, margins, 40–600 kW | C–E |
| T9 | Report section (~8–12 pages): assumptions, methods, results, open risks | All |

---

## 9. Schedule (target: 5 weeks part-time)

| Week | Work |
|---|---|
| 1 | Part A: orbit geometry, view-factor integrator, flux model for all six altitudes; T1 draft; **altitude decision gate** |
| 2 | Part B case definitions; Layer 1 balance; Layer 2 node network built and validated against Layer 1 at 40 kW |
| 3 | Layer 3 Fourier sub-models; Layer 4 cross-check; radiator orientation and $T_{\text{rad}}$ trades |
| 4 | Array sizing and pointing; power sweep 40–600 kW, hot/cold/nominal |
| 5 | T9 report section; repo cleanup; write-up post |

---

## 10. Repo Additions

```
orbital-dc-radiation/
├── thermal/
│   ├── README.md                 # this plan
│   ├── altitude_trade/           # Part A: flux model, view factors, T1
│   ├── cases.md                  # Part B: hot/cold/nominal definitions
│   ├── nodes.md                  # T2: node network and conductor tables
│   ├── conduction_models.md      # T3: Fourier derivations
│   ├── model/                    # Layer 1–4 code (Python: numpy/scipy/matplotlib)
│   └── sweep/                    # T8: power sweep outputs, per-case settings.md
└── report/
    └── section_thermal.md        # T9
```

Every run records its inputs in a `settings.md` beside the outputs, as in Phase 1.

---

## 11. Risks and Guardrails

- **Altitude decision before modeling.** Do not build the node model until T1 is reviewed and an altitude is selected. The scoping work in this conversation used 1,200 km by default; that is not the selected altitude.
- **Worst-instant vs. average.** Sizing to the orbit-average rejection is a common error that under-sizes the radiator; always report worst-instant and state which was used.
- **View factors.** The $\sin^2\!\rho$ shortcut is only valid for a nadir-facing plate. Any tilted or edge-on surface needs the integrated view factor. Validate the integrator against $\sin^2\!\rho$ at $\theta = 0$ before use.
- **Cold case honesty.** The bounding cold case is payload-off in eclipse, not payload-on. Do not let the payload's own dissipation mask the survival-heater and battery-temperature problem.
- **Mass model sourcing.** Radiator $\mu_{\text{rad}}$ (kg/m²) and array specific power (W/kg) dominate the answer. Each carries a cited basis or is labeled an assumption.
- **Scope.** No CFD, no two-phase loop modeling, no detailed deployment mechanism design, no economics. Park in `phase2_backlog.md`.
- **Public-data discipline.** As in Phase 1 — every reference number traces to the filing, a datasheet, or a handbook.

---

## 12. Open Items for the Author

1. Confirm the meaning and intended scope of "FTE analysis" (§5.4).
2. Confirm altitude set for Part A: the six Phase 1 B-case altitudes, or add 600 km from the start?
3. Confirm the cold case definition: payload-off in eclipse (proposed) vs. payload-on at minimum load.
4. Coolant working fluid: ammonia (proposed) vs. water-glycol vs. two-phase.
5. Whether the battery mass and array mass models belong in this aside or are deferred to Phase 2 with only the thermal-side sizing kept here.
