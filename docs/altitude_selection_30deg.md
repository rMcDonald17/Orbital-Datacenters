# Altitude Selection for the 30° Thermal Study

**Status:** Thermal aside, Part A — deliverable T1
**Scope:** 30° inclination shells (Phase 1 B-cases) only
**Decision:** 800 km, carried forward to Parts B–E
**Code:** [`flux.py`](flux.py), [`fom.py`](fom.py), [`drag.py`](drag.py), [`disposal.py`](disposal.py)

---

## 1. The question, and why the answer is not thermal

Part A was scoped to find the altitude within the filed 500–2,000 km range that minimises
the thermal burden on a 30° shell, accounting for Earth IR, albedo, direct solar and
eclipse.

It does not have a thermal answer. The thermal figure of merit varies by **7.6% across a
4× altitude range**, which is smaller than the uncertainty in any of the mass models
feeding it. Thermal environment does not discriminate between altitudes for this orbit
family.

What does discriminate is propulsion. Atmospheric drag imposes a hard floor near 600–700
km, and the FCC post-mission disposal requirement imposes a soft ceiling. The altitude
selection is therefore made on a drag-and-disposal trade with a thermal tiebreak, which
is not what this part was expected to produce.

---

## 2. Method

### 2.1 View factor model

Surface-to-Earth view factors are computed by numerical integration of

$$F = \frac{1}{\pi} \int_{\Omega_\oplus \cap \, \hat{n} \cdot \hat{v} > 0} \cos\alpha \, d\Omega$$

over the Earth's angular cone of half-angle $\rho = \arcsin\left(R_\oplus/(R_\oplus+h)\right)$,
where $\alpha$ is the angle from the surface normal. The $\sin^2\!\rho$ shortcut is used
only as a validation check at $\theta = 0$, not as the working model, because every
candidate radiator orientation except nadir is tilted or edge-on.

**Validation.** The integrator reproduces $\sin^2\!\rho$ at $\theta = 0$ to within 0.001%
at all six altitudes.

| $h$ (km) | $F_{\text{int}}$ | $\sin^2\!\rho$ | error |
|---:|---:|---:|---:|
| 500 | 0.85975 | 0.85976 | 0.001% |
| 1,000 | 0.74707 | 0.74707 | 0.001% |
| 2,000 | 0.57924 | 0.57924 | 0.000% |

The result that matters is at $\theta = 90^\circ$: a panel lying edge-on in the orbit
plane still sees $F = 0.118$ to $0.267$ depending on altitude, **not zero**. Earth's
angular radius spills past the panel's own horizon.

| $h$ (km) | nadir | zenith | edge-on ($\theta = 90^\circ$) |
|---:|---:|---:|---:|
| 500 | 0.8597 | 0 | 0.2673 |
| 700 | 0.8118 | 0 | 0.2328 |
| 1,000 | 0.7471 | 0 | 0.1939 |
| 1,200 | 0.7081 | 0 | 0.1736 |
| 1,500 | 0.6552 | 0 | 0.1489 |
| 2,000 | 0.5792 | 0 | 0.1182 |

### 2.2 Orbit and flux model

Circular orbit, $i = 30^\circ$, nadir-locked attitude. Cylindrical umbra. Absorbed flux
per surface:

$$q_{\text{abs}} = \alpha_s q_{\text{sol}} \max(0, \cos\theta_{\text{sun}}) + \alpha_s a\, q_{\text{sol}} F K + \varepsilon q_{\text{IR}} F$$

with the three terms being direct solar, albedo, and Earth IR. Both solar terms are set to
zero in eclipse.

| Symbol | Meaning | Value used |
|---|---|---|
| $q_{\text{abs}}$ | Absorbed environmental flux on the surface | computed, W/m² |
| $\alpha_s$ | Solar absorptance of the surface | 0.09 (OSR) |
| $\varepsilon$ | IR emittance of the surface | 0.85 (OSR) |
| $q_{\text{sol}}$ | Solar constant | 1,361 W/m² (nominal) |
| $q_{\text{IR}}$ | Earth outgoing longwave flux, top of atmosphere | 237 W/m² (nominal) |
| $a$ | Earth albedo | 0.30 (nominal) |
| $F$ | View factor from the surface to Earth (§2.1) | 0 to 0.86 |
| $\theta_{\text{sun}}$ | Angle between surface normal and the Sun direction | varies over orbit |
| $\zeta$ | Solar zenith angle at the sub-satellite point | varies over orbit |
| $K$ | Solar-zenith correction on albedo, $K = \max(0, \cos\zeta)$ | 0 to 1 |
| $\sigma$ | Stefan–Boltzmann constant | $5.6704 \times 10^{-8}$ W m⁻² K⁻⁴ |
| $T_{\text{rad}}$ | Radiator effective radiating temperature | 318 K (held fixed) |
| $\beta$ | Angle between the Sun direction and the orbit plane | $-53.44^\circ$ to $+53.44^\circ$ |
| $\theta$ | Angle between a surface normal and the nadir vector | 0° nadir, 180° zenith |
| $\rho$ | Earth angular radius, $\arcsin\left(R_\oplus/(R_\oplus+h)\right)$ | 49.6° to 68.0° |
| $h$ | Orbit altitude above the equatorial radius | 500–2,000 km |
| $R_\oplus$ | Earth equatorial radius | 6,371 km |
| $i$ | Orbit inclination | 30° |
| $\epsilon$ | Obliquity of the ecliptic | 23.44° |

Note the two distinct uses of $\varepsilon$ and $\epsilon$: the former is surface IR
emittance, the latter obliquity. Nominal environmental values are used throughout Part A;
the hot and cold bounding values are applied in Part B.

Net rejection per unit radiator area, the quantity ranked in §3, is then

$$q_{\text{net}} = \varepsilon \sigma T_{\text{rad}}^4 - q_{\text{abs}}$$

evaluated at the orbit position and $\beta$ that minimise it.

**Radiator temperature is held fixed at $T_{\text{rad}} = 318$ K for the trade.** Selecting
$T_{\text{rad}}$ is a Part D question, and using a Part D output as a Part A input would be
circular. Holding it fixed is defensible because $\varepsilon\sigma T^4$ is
altitude-independent while $q_{\text{abs}}$ is not, so the altitude *ranking* is
independent of the value chosen; only the absolute rejection level moves.

### 2.3 Eclipse

Reproduces the Phase 1 eclipse geometry results to within rounding (500 km: 37.75% /
35.7 min here against 37.8% / 35.8 min in `docs/eclipse_geometry.md`), confirming the two
models agree. Eclipse geometry was not recomputed as new work.

$f_{\text{ecl}}$ is the **eclipse fraction** — the fraction of each revolution spent in
the Earth's umbra, given for a circular orbit by

$$f_{\text{ecl}} = \frac{1}{\pi} \arccos\!\left[ \frac{\sqrt{h^2 + 2 R_\oplus h}}{(R_\oplus + h) \cos\beta} \right]$$

and zero when the argument reaches unity. $\beta_{\max} = i + \epsilon = 53.44^\circ$ is
the hard ceiling on the sun/orbit-plane angle for a 30° orbit, with $\epsilon = 23.44^\circ$
the obliquity. The two columns therefore bracket the annual range: $\beta = 0$ is the
deepest shadow pass, $\beta_{\max}$ the shallowest.

| $h$ (km) | Period (min) | $f_{\text{ecl}}$ at $\beta = 0$ | Duration (min) | $f_{\text{ecl}}$ at $\beta_{\max}$ |
|---:|---:|---:|---:|---:|
| 500 | 94.5 | 0.378 | 35.7 | 0.284 |
| 700 | 98.6 | 0.358 | 35.3 | 0.241 |
| 1,000 | 105.0 | 0.333 | 34.9 | 0.181 |
| 1,200 | 109.3 | 0.319 | 34.8 | 0.139 |
| 1,500 | 115.8 | 0.301 | 34.8 | 0.054 |
| 2,000 | 127.0 | 0.276 | 35.0 | 0.000 |

---

## 3. Radiator orientation: zenith wins at every altitude

A 30° shell's node regresses 3.3–6.6°/day, cycling $\beta$ through its full range every
47–84 days. A fixed radiator must therefore survive the **worst** $\beta$, not a chosen
one. Evaluating each orientation at its own worst $\beta$ and worst orbit position:

Worst-instant net rejection, W/m², at $T_{\text{rad}} = 318$ K:

| $h$ (km) | nadir | zenith | orbit-normal | ram / wake |
|---:|---:|---:|---:|---:|
| 500 | 274 | **370** | 335 | 316 |
| 700 | 276 | **370** | 343 | 323 |
| 1,000 | 281 | **370** | 351 | 331 |
| 1,200 | 284 | **370** | 356 | 335 |
| 1,500 | 289 | **370** | 361 | 340 |
| 2,000 | 297 | **370** | 368 | 347 |

At $\beta = 0$ the deployed orbit-normal panel appears best — 452 W/m² at 1,200 km,
because the sun grazes it and its Earth view is small. But its solar load scales as
$\sin\beta$, and by $\beta = 53.4^\circ$ it has degraded below zenith at every altitude.

**The zenith face is $\beta$-invariant at 370 W/m² everywhere.** It sees no Earth
($F = 0$ exactly for a nadir-locked bus) and its peak solar load is bounded at
$\alpha_s q_{\text{sol}} = 122$ W/m² regardless of geometry. It is the only orientation
whose sizing does not move with the nodal cycle, and that robustness — not its peak
performance — is why it wins.

**Consequence: radiator area is altitude-invariant.** Since zenith wins at 370 W/m² at
every altitude, $A_{\text{rad}}$ is identical across the range. The largest single thermal
mass item has no altitude dependence at all.

---

## 4. Thermal figure of merit

40 kW payload, $\dot{Q}_{\text{total}} = 42.4$ kW, zenith radiator, +25% area margin.
Mass models are **assumption-grade and unsourced**: radiator 10 kg/m² of physical panel,
array 150 W/kg, cells 150 Wh/kg at 70% DoD.

$\Phi$ is the **array energy-balance factor** — the ratio of array power required to load
power, accounting for the array having to carry the load *and* recharge the battery during
the sunlit fraction of each orbit:

$$\Phi = 1 + \frac{\bar{f}_{\text{ecl}}}{\left( 1 - \bar{f}_{\text{ecl}} \right) \eta_{\text{charge}}}$$

with $\bar{f}_{\text{ecl}}$ the orbit-averaged eclipse fraction over the annual $\beta$
sweep and $\eta_{\text{charge}} = 0.85$. Lower is better; $\Phi = 1$ would mean a
continuously sunlit orbit needing no storage.

| $h$ (km) | $\Phi$ | $A_{\text{rad}}$ (m²) | $m_{\text{rad}}$ | $m_{\text{arr}}$ | $m_{\text{bat}}$ | Total (kg) | kg/kW | cycles/yr |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 1.64 | 143 | 715 | 515 | 240 | 1,470 | 36.8 | 5,564 |
| 700 | 1.57 | 143 | 715 | 493 | 237 | 1,445 | 36.1 | 5,329 |
| 1,000 | 1.49 | 143 | 715 | 467 | 235 | 1,417 | 35.4 | 5,007 |
| 1,200 | 1.44 | 143 | 715 | 453 | 234 | 1,403 | 35.1 | 4,810 |
| 1,500 | 1.38 | 143 | 715 | 435 | 234 | 1,384 | 34.6 | 4,538 |
| 2,000 | 1.30 | 143 | 715 | 410 | 236 | 1,361 | 34.0 | 4,137 |

Three structural findings:

**Radiator mass is altitude-invariant** — §3, and it is 51% of the thermal total.

**Battery mass is altitude-invariant.** Eclipse *duration* holds at 34.8–35.7 min across
the whole range because the falling eclipse fraction is cancelled by the lengthening
orbital period. Capacity scales with duration × power, so it does not move. What improves
with altitude is cycle count — 5,564/yr at 500 km against 4,137/yr at 2,000 km, a 26%
reduction against the 5-year cadence. That is a battery *life* argument, not a mass one,
and it is the only genuine thermal-side reason to prefer altitude.

**Only the array varies.** Orbit-averaged $f_{\text{ecl}}$ falls from 0.352 to 0.205, so
$\Phi$ improves from 1.64 to 1.30 and array mass drops 20%. This is the entire 7.6%
gradient.

---

## 5. Drag: the constraint that actually bounds the bottom

Thermal and power hardware is area, and area is drag. The 40 kW vehicle carries ~247 m²
of array and ~72 m² of physical radiator panel against ~2,100 kg:

| Attitude | $A$ (m²) | $A/m$ (m²/kg) | $B = m/(C_D A)$ (kg/m²) |
|---|---:|---:|---:|
| Broadside (worst) | 330 | 0.157 | 2.9 |
| Tumble-average | 171 | 0.082 | 5.6 |
| Edge-on (best) | 60 | 0.029 | 16.0 |

For reference, Starlink V2 Mini is roughly 0.03–0.04 m²/kg and ISS about 0.006. **This
vehicle is two to three times draggier than Starlink per unit mass**, and Starlink
station-keeps continuously to hold 480 km.

Passive orbit lifetime, tumble-average $A/m$, years to 200 km:

| $h$ (km) | Solar min | Solar mean | Solar max |
|---:|---:|---:|---:|
| 500 | 1.19 | 0.30 | 0.05 |
| 600 | 5.38 | 1.35 | 0.23 |
| 700 | 22.4 | 5.60 | 0.94 |
| 800 | >25 | 19.7 | 3.29 |
| 1,000 | >25 | >25 | 22.1 |

**500 km is excluded outright.** At mean solar activity the vehicle deorbits in 0.30
years; at solar max, 18 days. Minimum altitude for 5-year passive lifetime is 594 km at
solar min, 690 km at mean, and **833 km at solar max**. A five-year mission spans a
substantial fraction of the solar cycle, so it must be sized against solar max.

Station-keeping drag make-up, worst case, Isp 1500 s:

| $h$ (km) | $\Delta V$/yr (max) | 5-yr $\Delta V$ | Propellant |
|---:|---:|---:|---:|
| 500 | 1,373 | 6,865 | 783 kg |
| 600 | 282 | 1,412 | 192 kg |
| 700 | 69 | 346 | 49 kg |
| 800 | 22 | 110 | 16 kg |
| 1,000 | 5.5 | 28 | 4 kg |

The knee is sharp: a factor of 16 in propellant between 500 and 700 km.

---

## 6. Disposal: the constraint that bounds the top

FCC 22-74 (adopted 29 Sep 2022, effective 29 Sep 2024, codified at 47 CFR §25.283)
requires LEO spacecraft to be disposed of within five years of mission end. It applies to
applications filed after the effective date, which includes the January 2026 reference
filing. The same order requires casualty risk from reentering debris below 1 in 10,000.

Worst case for *disposal* is solar **minimum** — thin atmosphere, slow decay — the inverse
of the station-keeping case. Decay time at solar min sets the disposal target:

| Drop to | Decay (solar min) |
|---:|---:|
| 500 km | 1.19 yr |
| 550 km | 2.58 yr |
| 600 km | 5.38 yr |
| 650 km | 11.6 yr |

550 km is the disposal target with margin. Disposal cost from each operating altitude:

| $h$ (km) | Drag-assist to 550 km | Controlled, perigee 60 km | EP prop (Isp 1500) | Chem prop (Isp 220) |
|---:|---:|---:|---:|---:|
| 700 | 81 m/s | 180 m/s | 12 kg | 168 kg |
| 800 | 133 | 206 | 19 kg | 191 kg |
| 1,000 | 235 | 255 | 33 kg | 234 kg |
| 1,200 | 333 | 302 | 47 kg | 274 kg |
| 1,500 | 472 | 368 | 66 kg | 329 kg |

---

## 7. The selection: 800 km

Combining station-keeping and disposal over a 5-year mission, electric propulsion:

| $h$ (km) | SK $\Delta V$ | Disposal $\Delta V$ | Total | Propellant | % of dry mass |
|---:|---:|---:|---:|---:|---:|
| 700 | 346 | 81 | 427 | 60 kg | 2.9% |
| **800** | **110** | **133** | **244** | **35 kg** | **1.6%** |
| 1,000 | 28 | 235 | 263 | 37 kg | 1.8% |
| 1,200 | 13 | 333 | 346 | 49 kg | 2.3% |

**800 km is the minimum of the propulsion trade.** Station-keeping falls steeply with
altitude while disposal rises roughly linearly; the two cross near 800 km. The minimum is
shallow — everything from 700 to 1,200 km sits between 1.6% and 2.9% of dry mass — but it
is a real minimum, and it is the only quantitative discriminator available, since thermal
is flat to 7.6% and radiation at these altitudes is far inside budget.

Supporting arguments:

- **Radiation is not binding.** Interpolating log-linearly between the Phase 1 B-case
  results at 700 km (0.76 krad(Si)/yr) and 1,000 km (3.02) gives roughly **1.2 krad/yr at
  800 km behind 5 mm Al**, or ~6 krad over five years against the assumed 100 krad COTS
  budget. Every altitude in 700–1,200 km is comfortable; radiation does not force the
  choice within the band.
- **Thermal is indifferent.** Interpolated FOM at 800 km is ~35.9 kg/kW, against 36.1 at
  700 and 35.4 at 1,000. The spread is 2%.
- **Drag margin is adequate.** At 800 km the vehicle survives 3.3 years passively even at
  solar max, so a station-keeping outage is recoverable rather than immediately fatal —
  unlike at 500–600 km, where loss of propulsion means loss of vehicle within months.

### What 800 km costs

It is **below the 833 km solar-max passive-decay floor**, so station-keeping is mandatory
rather than optional. The vehicle cannot be flown as a free-flyer. At 110 m/s over five
years this is a small propellant line, but it is a reliability requirement on the
propulsion system for the full mission duration.

It is also **not a Phase 1 case altitude.** The B-case matrix runs 500 / 700 / 1,000 /
1,200 / 1,500 / 2,000 km, so 800 km has no radiation counterpart and its dose figure above
is an interpolation, not a SPENVIS result. Either a B-case run should be added at 800 km,
or the thermal results should be reported alongside the 700 and 1,000 km radiation cases
as a bracket.

---

## 8. Limitations

1. **Single-surface radiation model.** Every orientation in §3 was evaluated as an
   isolated flat plate with a full hemisphere of space view minus Earth. No bus
   self-blockage, no array-to-radiator exchange, no inter-panel view factors. A solar
   array at 60–80 °C with back-face $\varepsilon \approx 0.8$ occupying 10% of a
   radiator's field of view contributes roughly 58 W/m² of absorbed flux — larger than the
   entire albedo term computed here. **370 W/m² is therefore an optimistic bound**, and
   realistic values with a defined vehicle geometry may fall to 300–330 W/m². All areas
   and masses in §4 scale accordingly.
2. **Nadir-locked attitude assumed, not justified.** Optical intersatellite links and
   Ka-band TT&C do not obviously require nadir pointing. Attitude is a free design
   variable that simultaneously sets array pointing, radiator pointing, and what each
   sees. It was fixed by assumption here and a surface was optimised within it.
3. **The zenith conclusion is the most fragile result.** It rests on $F_{\text{earth}} = 0$
   exactly, which is true for an isolated plate and false once a bus sits underneath it.
4. **Screening-grade atmosphere.** Exponential model with a static solar-activity
   multiplier, not NRLMSISE-00 or JB2008, and no phased solar cycle. The shape — cliff
   below 600 km, shallow minimum, flatness across 700–1,200 — is robust. The 700-versus-800
   km distinction is well inside the error bar.
5. **$A/m$ is derived from the unsourced mass models of §4** and a guessed array specific
   area. Orbit lifetime scales linearly with it; a 25% error moves the drag floor by
   50–100 km.
6. **Disposal mode unresolved.** Whether drag-assisted disposal satisfies the 1-in-10,000
   casualty requirement for a 2.1 t (40 kW) to 31 t (600 kW) vehicle has not been
   assessed. If controlled reentry is required, chemical propellant runs 191 kg at 800 km
   — 9% of dry mass — and the altitude gradient reverses to favour lower altitudes
   monotonically. **This could change the selection.**
7. **Anchored to 40 kW.** Area and mass both scale roughly linearly with power, so $A/m$
   should hold to 600 kW, but the radiator-to-array area ratio does shift with power and
   this has not been checked.

---

## 9. Future work — dedicated altitude study

The altitude selection above is adequate to proceed with Parts B–E but is **not a
defensible final answer**. It is a screening result built on an optimistic thermal model,
a screening atmosphere, and unsourced mass models, and it was driven by a propulsion trade
that this study was not scoped to perform.

A proper altitude study is tracked as a separate work item and should cover:

- **Real atmosphere model.** NRLMSISE-00 or JB2008 with a phased solar cycle over the
  mission window, rather than a static multiplier. Monte Carlo over launch epoch, since
  where a 5-year mission falls in the cycle changes drag by an order of magnitude.
- **Vehicle-level radiation model.** Inter-surface view factors, array-to-radiator
  coupling, bus self-blockage, and attitude as a free variable rather than an assumption.
  This is the largest single source of error in §3 and §4.
- **Demise and casualty analysis** to resolve the disposal mode. This is a binary that
  swings the altitude recommendation.
- **Sourced mass models** for radiator areal mass, array specific power, and cell energy
  density. Radiator kg/m² dominates both the thermal FOM and $A/m$.
- **Conjunction environment.** SpaceX lowered ~4,400 Starlink satellites from 550 to 480
  km during 2026 explicitly to reduce collision risk and shorten disposal time. The
  700–1,200 km band is a historically congested region. Object density versus altitude is
  a selection criterion this study has not considered at all.
- **Coupled optimisation at 600 kW**, where the vehicle is ~31 t with roughly 15× the area
  and the drag and disposal numbers may not scale as assumed.
- **Extension to the sun-synchronous shells**, where the eclipse-free condition above
  ~1,200 km changes the trade entirely and the drag floor is far less binding.

Until that work is done, **800 km should be read as a working assumption for the thermal
model, not as a result of this study.**
