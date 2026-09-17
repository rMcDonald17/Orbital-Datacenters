1. Purpose

Phase 1 answers the radiation cost of altitude. This aside answers the thermal cost of altitude for the same filed range, then builds an analytical thermal model of a single satellite at the altitude that trade selects.

Design questions this aside answers:

Within the filed 500–2,000 km range, which altitude minimizes the thermal burden on a 30° shell, accounting for Earth IR, albedo, direct solar, and eclipse?
At that altitude, what radiator and solar array area, orientation, and mass does a 40 kW compute payload require, and how does that scale to 600 kW?
What are the bounding hot and cold cases, and does the vehicle survive both?

Why 30° only. A 30° orbit has |β| ≤ i + 23.44° = 53.4°, and the no-eclipse threshold β* = arcsin(Rₑ/(Rₑ+h)) does not drop below 53.4° until h ≈ 1,565 km. Every altitude in the filed range from 500 km to ~1,565 km therefore eclipses year-round, and above that only briefly at peak β. The sunlight argument that favors high altitude for the sun-synchronous shells does not apply here, so the altitude trade for this shell is a clean contest between radiation (prefers low) and thermal environment (to be determined by Part A).

2. Fixed Inputs
2.1 Reference scenario (from Phase 1, unchanged)
Parameter	Value	Source basis
Altitude range	500–2,000 km	FCC filing
Inclination	~30°	FCC filing
Bus heritage	Starlink V3-class flat-panel, scaled	Musk statements / filing context
Launch vehicle	Starship, ~100 t to LEO, ~9 m fairing	Public Starship specs
Reference compute payload	H100 SXM-class rack, ~40 kW	Analyst choice (consistent with Phase 1)
Power sweep	40 / 120 / 190 / 370 / 600 kW	H100 rack → GB200 → Vera Rubin → Rubin CPX → Rubin Ultra rack-class power (public roadmap figures)
Replacement cadence	~5 years	Analyst assumption
2.2 Environmental constants (Gilmore, Spacecraft Thermal Control Handbook, LEO values)
Quantity	Cold case	Nominal	Hot case
Solar constant (W/m²)	1,322 (aphelion)	1,361	1,414 (perihelion)
Earth IR, top of atmosphere (W/m²)	218	237	258
Albedo	0.25	0.30	0.35
Space sink temperature (K)	3	3	3
2.3 Thermo-optical properties
Surface	α_s (BOL)	α_s (EOL)	ε_IR	Notes
Radiator (OSR / silvered Teflon)	0.08–0.10	0.15–0.20	0.80–0.85	EOL degradation from UV + contamination; EOL used in hot case
Solar array (front, cell side)	~0.90	~0.92	~0.85	Absorbed minus electrical output is the heat load
Solar array (back, substrate)	—	—	0.80–0.85	Radiates to space; sets array temperature
Bus MLI (external)	0.3–0.5	—	ε* ≈ 0.01–0.03	Effective emittance, not surface emittance
2.4 Payload thermal constraints (public datasheets)
Component	Allowable range	Driver
H100 SXM junction	≤ 85 °C	Hot case: sets coolant outlet ceiling
HBM	≤ 95 °C (typ.)	Confirm against datasheet
Li-ion battery, operating	0 to +40 °C	Cold case driver — narrowest allowable range on the vehicle
Li-ion battery, charge	+5 to +35 °C	Constrains post-eclipse recharge
Coolant (ammonia)	> −77 °C freeze	Cold-case check on deployed lines
Optical terminal	vendor TBD	Placeholder ±20 °C, stability-limited
3. Part A — Altitude Trade for Thermal Management

Objective: produce a thermal figure of merit vs. altitude for the 30° shell and select one altitude for Parts B–E.

3.1 Altitude cases

Use the Phase 1 B-case altitudes so thermal results align with the radiation database: 500 / 700 / 1,000 / 1,200 / 1,500 / 2,000 km. Add 600 km if the trade shows a gradient worth resolving between B1 and B2.

3.2 Orbit geometry per altitude

For each altitude compute, over the annual β sweep (−53.4° to +53.4°):

Orbital period, Earth angular radius ρ = arcsin(Rₑ/(Rₑ+h))
Eclipse fraction and duration vs. β (cylindrical shadow; penumbra noted but neglected)
Annual sunlit fraction
Worst-case (β = 0) and best-case (|β| = 53.4°) eclipse durations

Scoping note: eclipse duration at β = 0 is nearly flat (~35 min) from 500 to 1,200 km because period lengthens as the shadow fraction shrinks. Expect battery energy to be weakly altitude-dependent; expect the trade to be driven by Earth IR and albedo view factors instead.

3.3 Environmental flux per altitude

For a set of reference surface orientations (nadir, zenith, orbit-normal, ram, wake), compute per altitude and per orbit position:

View factor to Earth, F(θ, ρ), by numerical integration over the Earth disk (do not use the sin²ρ shortcut except for θ = 0 validation)
Earth IR incident: q_IR · F
Albedo incident: a · q_sol · F · K(ζ), with K the solar-zenith correction at the sub-satellite point
Direct solar: q_sol · max(0, cos θ_sun), zero in eclipse
3.4 Figure of merit

Two metrics, reported separately, then combined:

Net radiator rejection, worst instant (W/m²): εσT⁴ − absorbed environmental flux, evaluated at the orbit position that minimizes it, for the best available radiator orientation at that altitude. Higher is better.
Array energy-balance factor: required array power ÷ load power, = 1 + (f_ecl / (1 − f_ecl)) / η_charge. Lower is better.

Combined: thermal-driven specific mass (kg/kW) = radiator mass + array mass + battery mass per kW of compute, using the mass models in §7. This is the number that goes on the altitude trade chart alongside the Phase 1 TID curve.

3.5 Deliverable T1 — Altitude trade
Table: period, eclipse, F_earth by orientation, q_IR, q_alb, net rejection, energy-balance factor, thermal kg/kW, for all six altitudes at β = 0 and |β| = 53.4°
Chart: thermal kg/kW vs. altitude, overlaid on the Phase 1 TID/year curve for the 30° shell (same x-axis) — this is the aside's money plot
One-paragraph selection rationale

Decision gate: altitude selection is reviewed and approved before Part B begins.

4. Part B — Bounding Cases at the Selected Altitude

Define once, apply to every node-model run.

	Hot case	Cold case
Solar constant	1,414 W/m²	1,322 W/m²
Earth IR	258 W/m²	218 W/m²
Albedo	0.35	0.25
β	Value maximizing solar on radiator + array (from Part A sweep)	0° (longest eclipse)
Optical properties	EOL (α_s degraded, ε at low bound)	BOL
Payload state	Full power, all racks	Safe mode / payload off through eclipse
Bus loads	Max	Min (survival)
Array	EOL, degraded per Phase 1 DDD	BOL, cold, high-voltage check
Sizing question answered	Radiator area	Survival heater power, battery temperature, coolant freeze margin

The cold case is not "eclipse with payload on" — the payload's own dissipation keeps the loop warm in that case. The bounding cold case is a payload outage or safe-mode entry during eclipse, which is exactly when a fault is most likely to occur and when heater power competes with a discharging battery.

Add a nominal case (β = 0, nominal constants, BOL, full power) for the parametric sweep and for validation against the algebraic sizing.

5. Part C — Analytical Thermal Model

Built up in four layers, each validating the next. Terminology follows Thermal Desktop / SINDA conventions so the model maps directly onto that tool if the study later moves there.

5.1 Layer 1 — Steady-state energy balance (algebraic)

Whole-vehicle balance, one equation per surface class:

Q_payload + Q_bus + Σ α_s A_i (q_sol cos θ_i + q_alb F_i) + Σ ε_i A_i q_IR F_i = Σ ε_i A_i σ T_i⁴

Solve for radiator area at fixed T_rad. This layer alone produces the first radiator-area-vs-power curve and is the cross-check for everything below. Reproduce the Phase 1 scoping results (40 kW → ~117 m² radiating surface at 318 K, deployed panel, 1,200 km) as the regression test, then re-run at the Part A altitude.

5.2 Layer 2 — Lumped-parameter node network

Thermal-Desktop-style network. Each node has capacitance C_i = m_i c_p,i; conductors link nodes; environmental loads apply per node from the Part A flux model.

Node list (40 kW baseline):

Node	Type	Contents	Notes
N1	Diffusion	GPU die + package + cold plate (lumped)	Lumping avoids stiffness; Biot check justifies it
N2	Diffusion	Coolant loop inventory	Single-phase pumped ammonia; ΔT ≈ 10 K design
N3	Diffusion	Radiator root / manifold	Linear conductor to N2 (convective, h·A)
N4	Diffusion	Radiator panel (fin)	Either one node with fin efficiency, or resolved 1-D (Layer 3)
N5	Diffusion	Bus structure + MLI-wrapped electronics	Low-power, MLI-isolated
N6	Diffusion	Solar array	Radiatively coupled only; own energy balance
N7	Diffusion	Battery	Narrowest allowable — tracked explicitly
N8	Boundary	Deep space	3 K, infinite capacitance
—	Heat load	Earth IR, albedo, solar	Time-varying per node from Part A

Conductor types:

Linear: G_ij = k A / L (conduction) or h A (convective, loop-to-wall)
Radiation: R_ij = ε_i A_i F_ij σ — "radks" in SINDA; space is the dominant sink for N4, N6
Contact: G = h_c A at bolted/bonded interfaces (cold plate to package)

Governing ODE per node:

C_i dT_i/dt = Q_i(t) + Σ_j G_ij (T_j − T_i) + Σ_j R_ij (T_j⁴ − T_i⁴)

Integrate over ≥ 3 orbits to a periodic limit cycle; report the last orbit. Use a stiff-capable integrator (Radau / BDF) if N1 is split into die and plate.

Deliverable T2: node network diagram + conductor table + capacitance table, as a standalone nodes.md that a Thermal Desktop user could rebuild from.

5.3 Layer 3 — Fourier conduction sub-models

Analytical 1-D solutions used to derive the lumped parameters in Layer 2, and to check them:

Radiator fin: rectangular fin radiating to space, root at T_root. Governing equation k t d²T/dx² = εσ(T⁴ − T_space⁴) has no closed form for T⁴, so use (a) the linearized fin solution with h_rad = 4εσT̄³ for the efficiency estimate, and (b) a shooting or finite-difference solution of the nonlinear form as the check. Output: fin efficiency η_fin vs. fin length and thickness → feeds N4 and the radiator mass model.
Cold plate spreading: point/area heat source on a finite plate — spreading resistance per Lee et al. or Yovanovich. Output: die-to-coolant resistance → conductor N1–N2.
Heat pipe / loop conductance: effective conductance of the transport element between root and fin → conductor N3–N4.
MLI: effective emittance ε* model for the bus wrap → radiation conductor N5–N8.

Deliverable T3: conduction_models.md with each derivation, its assumptions, and the resulting conductor/efficiency value used in Layer 2.

5.4 Layer 4 — Discretized cross-check

Finite-difference (or finite-element, if tooling is available) solution of the radiator fin and cold plate, compared against the Layer 3 analytical results. Purpose: confirm that the analytical fin efficiency and spreading resistance are within tolerance of a discretized solve before trusting them in the sweep. Document the discrepancy and its cause, as with the Phase 1 transport-vs-SHIELDOSE validation.

Author to confirm the intended scope of "FTE analysis" — interpreted here as this discretized finite-element/finite-difference thermal cross-check. Adjust if a different method was meant.

5.5 Power sweep

Run Layers 1–2 at 40 / 120 / 190 / 370 / 600 kW, hot / cold / nominal. Layer 3 parameters are held fixed unless fin geometry changes with scale. Report:

Radiator area, physical panel dimensions, mass
Array area and mass
Battery capacity and mass
Node temperatures over the limit-cycle orbit, hot and cold
Margin against every allowable in §2.4
6. Part D — Radiator Sizing and Orientation
6.1 Orientation trade

Evaluate at the Part A altitude, β = 0 and β = 53.4°, hot case:

Configuration	Sun exposure	Earth view	Deployment
Body-mounted, zenith face	Direct at noon	F ≈ 0	None
Body-mounted, nadir face	None	F = sin²ρ (max)	None
Body-mounted, ram/wake faces	Once per orbit	Partial	None
Deployed, edge-on in orbit plane	Grazing at β = 0; rises with β	F(θ = 90°) ≈ 0.17 at 1,200 km — not zero	Single or multi-fold
Deployed, sun-tracking edge-on	Grazing always	Varies	Adds articulation

Size to the worst instant in the orbit, not the orbit average. Report both so the penalty of the conservative choice is visible.

6.2 Radiator temperature selection

Sweep T_rad = 290 / 300 / 310 / 318 / 333 / 350 K. Radiator area scales as 1/(εσT⁴ − q_abs); report the area and mass at each. Tie T_rad to the coolant chain: T_rad ≈ T_coolant,in − ΔT_root, with T_coolant,out bounded by the junction limit. State the recommended operating point and the leakage/reliability trade it implies for the silicon.

6.3 Geometry

Double-sided deployed panels: radiating area = 2 × physical area. Check physical dimensions against the ~9 m fairing envelope; note where multi-fold deployment becomes necessary in the power sweep.

6.4 Mass model

Deployed radiator assembly: areal mass in kg per m² of physical panel, range 5 (panel only, optimistic) to 10–15 (with fluid, manifolds, deployment structure). Source flight data (ISS PVR/HRS radiators, commercial deployables) and record the basis. The radiator is expected to be the largest single mass item; this number deserves the most scrutiny in the aside.

7. Part E — Solar Array Sizing and Pointing
7.1 Energy balance
P_array,required = (Q_payload + Q_bus) × [1 + f_ecl / ((1 − f_ecl) η_charge)] / η_PMAD

Size at EOL using the Phase 1 DDD result for cell degradation (cross-link to environment/ D1 database and the MC-SCREAM / EQFLUX runs already in the SPENVIS workflow).

7.2 Pointing

At β = 0 the sun lies in the orbit plane, so a fixed array on a nadir-locked bus sees cos-varying incidence around the orbit. Options:

Single-axis tracking about orbit-normal (removes the in-plane cos loss; standard for LEO)
Second axis or seasonal offset for the β sweep (±53.4°) — or accept cos β loss and oversize
Fixed, oversized array (simplest, heaviest)

Report array area and mass for each; select one for the sweep.

7.3 Array thermal node

N6 energy balance: absorbed solar α_s q_sol cos θ − electrical output η_cell(T) q_sol cos θ = ε_f σ T⁴ + ε_b σ T⁴ + small conduction to bus. Include the cell efficiency temperature coefficient so array output and temperature are solved together. Cold case: post-eclipse open-circuit voltage at minimum temperature.

7.4 Mass model

Specific power (W/kg) at array level, flexible-blanket vs. rigid, with source. Range 100–200 W/kg; record basis.

8. Deliverables
ID	Deliverable	Part
T1	Altitude trade table + thermal kg/kW vs. altitude chart overlaid on TID curve	A
T2	Node network diagram, conductor and capacitance tables (nodes.md)	C
T3	Fourier sub-model derivations (conduction_models.md)	C
T4	Discretized cross-check vs. analytical, with discrepancy discussion	C
T5	Hot/cold/nominal case definitions (cases.md)	B
T6	Radiator orientation and temperature trade	D
T7	Array sizing and pointing trade	E
T8	Power sweep results: area, mass, temperatures, margins, 40–600 kW	C–E
T9	Report section (~8–12 pages): assumptions, methods, results, open risks	All
