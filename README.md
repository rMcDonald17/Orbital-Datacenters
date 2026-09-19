# Orbital-Datacenters
Where you can actually put a datacenter in LEO, and what it costs you.

SpaceX's January 2026 FCC filing proposes an orbital compute constellation spanning
500–2,000 km, in sun-synchronous shells for continuously-sunlit baseline compute and
~30° shells for demand-peak capacity. This is an independent study of what that altitude
range costs in radiation, and whether the "over 99% sunlit" premise holds for the
baseline shells it was claimed for. The 30° shells are assessed against their own filed
role — demand-peak capacity — not against a continuity standard the filing never set for
them. Every input traces to public filings, published models, or datasheets.

Two constraints bound the answer from opposite ends. **Eclipse geometry rules out the
bottom of the range**: below ~1,200 km even a dawn–dusk sun-synchronous orbit eclipses,
and the 30° shells eclipse on essentially every revolution at every altitude. **Dose
scaling makes the top expensive**: annual TID rises 56× across the range for the SSO
shells and 269× for the 30° shells. What survives both is a narrower band than the one
filed.

A third constraint scales differently from either. Heat rejection rises linearly with
payload power against a fixed ceiling of a few hundred W/m² of blackbody radiating
capacity, and unlike shield mass it has no weak-lever escape. That work is under way in
`Thermal_Engineering/`, and its first result was not a thermal one: a compute satellite
is mostly radiator and solar array, which makes it two to three times draggier per unit
mass than Starlink, moving the altitude floor for the 30° shells from 500 km to ~800.

## The two orbit families

Twelve cases: six altitudes (500, 700, 1,000, 1,200, 1,500, 2,000 km) at two
inclinations, all circular.

| | **A-cases** | **B-cases** |
|---|---|---|
| Inclination | ~97–105°, sun-synchronous | 30° |
| Local time of ascending node | 6h (dawn–dusk) | n/a |
| Filing role | continuously-sunlit baseline compute | demand-peak capacity |
| Ground coverage | global, polar | ±30° latitude — the populated, high-demand band |
| Sees the poles | yes, every orbit | never |
| South Atlantic Anomaly | clips the edge on some orbits | crosses it on most orbits |
| Solar particle events | fully exposed over the poles | **immune** — cutoff excludes them entirely |
| Eclipse | none above ~1,400 km | every orbit, year-round — 64–76% sunlit |

The two differ in almost every way that matters for radiation, and not in the direction
intuition suggests. **A-cases take more solar proton dose; B-cases take more trapped
proton dose.** Which one is worse depends on altitude — see below.

## Eclipse geometry

The reference scenario's sun-synchronous shells are premised on being **over 99%
sunlit**. That premise is altitude-dependent, and it fails across the lower half of
the filed range.

A dawn–dusk SSO avoids eclipse when the sun/orbit-plane angle β exceeds
`asin(R⊕/(R⊕+h))`. Climbing lowers that threshold — but sun-synchronicity ties
inclination to altitude (97.4° at 500 km, 104.9° at 2,000 km), and higher inclination
lowers β. **The two effects oppose each other.**

For a dawn–dusk SSO the penalty has an exact closed form: `β_min = 156.56° − i`. Every
degree of inclination costs one degree of winter margin, so the 7.49° inclination rise
across the filed range costs 7.49° of β. Altitude still wins, because the threshold
falls 18.44° over the same span — but only by a factor of 2.5, and the returns diminish
steeply with height.

![SSO eclipse geometry](figures/beta_eclipse_sweep_sso.png)

| Altitude | Eclipse days/yr | Max eclipse | Annual sunlit |
|---:|---:|---:|---:|
| 500 km | 100 | 23.9% | 94.7% |
| 700 km | 80 | 19.0% | 96.7% |
| 1,000 km | 54 | 12.4% | 98.6% |
| 1,200 km | 34 | 7.9% | 99.4% |
| 1,500 km | 0 | — | 100% |
| 2,000 km | 0 | — | 100% |

Eclipse is a **winter phenomenon**, concentrated around the December solstice — every
shell is fully sunlit from roughly day 45 to day 300. A 500 km shell therefore needs
storage sized for a 22.7-minute worst-case pass, but exercises it ~1,520 times a year
rather than every revolution.

### The 30° shells never clear the shadow

Their RAAN is not frozen — it regresses 6.6°/day at 500 km, cycling the orbit plane
through every sun angle in ~47 days rather than following the seasons. And β at 30°
inclination cannot exceed `i + ε = 53.4°`, a hard geometric ceiling, while the
eclipse-free threshold runs from 68.0° down to 49.6°.

![30 degree eclipse geometry](figures/beta_eclipse_sweep_lowinc.png)

| Altitude | Eclipse-free days/yr | Max eclipse | Worst pass | Annual sunlit |
|---:|---:|---:|---:|---:|
| 500 km | 0 | 37.8% | 35.8 min | 64.0% |
| 1,000 km | 0 | 33.2% | 34.9 min | 69.2% |
| 1,500 km | 0 | 30.0% | 34.8 min | 73.2% |
| 2,000 km | 11 | 27.5% | 35.0 min | 76.4% |

**64–76% sunlit, eclipsing on essentially every revolution.** Worst-case pass stays at
~35 minutes at every altitude — the eclipse fraction falls but the orbit period grows,
and the two cancel. Climbing buys the 30° shells nothing on storage, and they cycle
their batteries ~5,000 times a year against ~1,520 for a 500 km SSO.

The >99%-sunlit premise is achievable only in dawn–dusk sun-synchronous orbit.
The model is validated against the SPENVIS orbit generator across all twelve cases:
inclination agrees to within 0.003°, and β to 0.014° for the sun-synchronous cases.
The larger residuals for the 30° cases are traced to J2 nodal regression, which the
analytic model neglects by construction.

Full derivation, limitations, and validation: **[docs/eclipse_geometry.md](docs/eclipse_geometry.md)**

## Radiation

Total ionising dose behind 5 mm Al, at the centre of an aluminium sphere, in
krad(Si)/yr:

| Altitude | SSO (A) | 30° (B) |
|---:|---:|---:|
| 500 km | 0.38 | 0.18 |
| 700 km | 0.70 | 0.76 |
| 1,000 km | 1.79 | 3.02 |
| 1,200 km | 3.31 | 6.57 |
| 1,500 km | 7.82 | 17.49 |
| 2,000 km | 21.51 | 47.81 |

![TID vs altitude](figures/tid_vs_altitude.png)

Dose rises **56× across the filed range for the SSO shells and 269× for the 30° shells**.
There is no sharp knee at the inner-belt onset; the climb is close to log-linear, roughly
doubling every 300 km above 1,000 km.

**The inclination trade reverses at ~600 km.** Below it the 30° shells are cleaner; above
it they take roughly twice the dose. The reason is that on trapped protons alone the 30°
orbit is worse at *every* altitude — 175 vs 97 rad/yr even at 500 km — because it crosses
the South Atlantic Anomaly on most revolutions while a polar orbit only clips its edge.
The SSO total is higher at 500 km solely because of a solar proton contribution the 30°
shells do not receive at all.

**The 30° shells receive exactly zero solar proton dose at any altitude.** At 30°
inclination the orbit never reaches invariant latitudes where geomagnetic cutoff rigidity
falls low enough for solar protons to arrive — SAPPHIRE's attenuation factor is
identically zero at every energy up to its 1 GeV ceiling. Any requirement for the compute
payload to power-safe during a solar particle event applies to the sun-synchronous shells
only.

Assuming a 100 krad COTS tolerance over a five-year replacement cadence — both analyst
assumptions — **radiation alone excludes no part of the filed range, provided shield
depth scales with altitude**: 2,000 km SSO needs ~10 mm Al where 500 km needs under 2 mm.
But shield mass is a weak lever up there. Going 5→10 mm at 2,000 km reduces annual dose by
only 36%, because the environment hardens with altitude as the electron component gives
way to a proton spectrum that thick aluminium cannot stop.

Full method, decomposition, and limitations:
**[docs/radiation_altitude_trade.md](docs/radiation_altitude_trade.md)**

## Drag and disposal

Radiation and eclipse say where a satellite survives and when it has power. Neither says
whether it stays in orbit.

A 40 kW compute satellite needs roughly 143 m² of radiating surface and 232 m² of solar
array against about 2.4 t — an area-to-mass ratio near 0.07 m²/kg, against 0.03–0.04
for
Starlink V2 Mini and 0.006 for the ISS. **Thermal and power hardware is area, and area is
drag.**

Passive orbit lifetime for the 30° shells — years to decay to 200 km with no
station-keeping, at three levels of solar activity:

| Altitude | Solar min (yr) | Solar mean (yr) | Solar max (yr) |
|---:|---:|---:|---:|
| 500 km | 1.39 | 0.35 | 0.07 |
| 600 km | 6.31 | 1.59 | 0.27 |
| 700 km | >25 | 6.57 | 1.11 |
| 800 km | >25 | 23.1 | 3.86 |
| 1,000 km | >25 | >25 | >25 |

**500 km is excluded outright** — the vehicle deorbits in months at mean solar activity
and in weeks at solar max, and holding it there would cost 774 kg of propellant over five
years, a third of its dry mass. The five-year passive floor sits at 820 km at solar max.

Above that, FCC 22-74 bounds the other end: LEO spacecraft must be disposed of within five
years of mission end, with reentry casualty risk below 1 in 10,000. Station-keeping cost
falls steeply with altitude while disposal cost rises, and the two cross near 800 km at
about 1.5% of dry mass. The minimum is shallow — everything from 700 to 1,200 km sits
between 1.5% and 2.5%.

**These are propulsion constraints, not thermal ones**, and they turn out to bound the
30° band more tightly than either radiation or heat rejection. The thermal environment
itself varies by only 7.6% across the whole filed range and does not discriminate between
altitudes at all.

Full method, validation, and limitations:
**[docs/altitude_selection_30deg.md](docs/altitude_selection_30deg.md)**

## Where that leaves the trade

The two constraints bound the useful band from opposite directions. Eclipse geometry
rules out the bottom of the filed range; dose scaling and shield mass make the top
expensive. Roughly **1,200–1,500 km, sun-synchronous** satisfies both — 99.4–100%
sunlit, 3.3–7.8 krad(Si)/yr at 5 mm, comfortably inside a 100 krad five-year budget.

The filed 2,000 km ceiling buys 0.6 percentage points of additional sunlight for 2.75×
the annual dose relative to 1,500 km.

The 30° shells are a different question, not a worse answer to the same one. They cover
the ±30° population band that polar shells reach only briefly, and the filing scopes them
as demand-peak capacity rather than continuous baseline — an intermittent role that a
64–76% sunlit duty cycle may suit rather than preclude. What they carry is a storage
penalty: ~5,000 charge cycles a year against ~1,520 for a 500 km SSO, and a ~35-minute
worst-case pass that climbing does not shorten.

That mass question now has a first answer. At 40 kW the storage penalty is ~236 kg of
cells, about 10% of vehicle mass — real but not disqualifying. What does constrain the
30° shells is drag: the area their radiators and arrays require rules out the bottom half
of the filed range, leaving roughly **700–1,000 km** once disposal cost is counted at the
top. That band is set by propulsion, not by the eclipse penalty the shells were doubted
for.

## What's next

**Thermal management.** The results here say where a satellite survives and when it has
power, not whether it can shed the heat its payload makes. Rack power is the open variable
and it is moving fast: ~40 kW for an H100-class rack, ~120 kW for GB200 NVL72, and
NVIDIA's public roadmap runs to ~600 kW with Rubin Ultra. Radiator area scales linearly
with every watt of it against a fixed radiating ceiling, so the thermal answer may bind
before dose does. `Thermal_Engineering/` takes the 30° shells first — they eclipse every
revolution, making them the harder and more informative thermal case. The altitude
trade is done
(above); what remains is a lumped-node model at 800 km with bounding hot and cold cases,
Fourier sub-models for the radiator fin and cold plate, and radiator and array sizing
across 40–600 kW. The open question there is where radiator area stops fitting a Starship
fairing: at 600 kW the vehicle is ~34 t with a 29 m radiator panel, one rack per launch.

**IRENE/AE9-AP9.** AP-8 and AE-8 are legacy models — epoch-limited, built largely from
1960s–70s data, and known to be conservative for electrons. More importantly they
represent the solar cycle with two static maps, so the MIN/MAX pair in the results above
is a spread between two model snapshots rather than a confidence interval. AE9/AP9 is a
statistical model with genuine percentile bands, and it is Distribution A, so its spectra
can be published directly where SPENVIS output cannot. Re-running the twelve cases
through IRENE closes both gaps at once: it quantifies how much conservatism AP-8/AE-8
carry, and it produces an environment database this repo can actually ship.

**Geant4 transport.** Every dose here is computed at the centre of a solid aluminium
sphere — a screening geometry, not a spacecraft. The next step is rebuilding the slab
problem in Geant4 against the same source spectra, verifying against SHIELDOSE-2, then
moving to a representative flat-panel bus with electronics behind realistic mass
distribution.

**A dedicated altitude study.** The drag and disposal results above are screening grade —
an exponential atmosphere with a static solar-activity multiplier, and mass models that
are still assumptions. A real answer needs NRLMSISE-00 or JB2008 with a phased solar
cycle, a demise analysis to settle whether controlled reentry is required (it would move
the selection to the bottom of the band), and the conjunction environment, which this
study has not considered at all — notably that SpaceX lowered ~4,400 Starlink satellites
from 550 to 480 km during 2026 to reduce collision risk, while the filed compute
constellation sits entirely above that band.

**Single-event effects.** Total dose says nothing about upset rates, which for a compute
payload may bind well before TID does. >30 MeV integral proton flux is already tabulated
as the precursor metric; turning it into an upset rate needs a device cross-section
survey.

## Acknowledgements

**AI tooling.** Parts of this study were developed with Anthropic's Claude — primarily
Claude Opus 5, with some work on Claude Fable 5. The models were used for derivation
checking, code structuring and refactoring, literature and regulatory lookup, and
drafting prose from results. Every physical model, numerical result and conclusion in this
repository was specified, reviewed and verified by the author; the view-factor integrator
is validated against a closed form, the eclipse model against the SPENVIS orbit generator,
and all environmental constants and regulatory citations are traced to the primary sources
listed inline. Errors are the author's.

<!-- Usage split not tracked per session; if you want a percentage here, log it going
     forward rather than estimating retrospectively. -->

Radiation environment data generated using SPENVIS (www.spenvis.oma.be), an ESA operational software system maintained by the Royal Belgian Institute for Space Aeronomy (BIRA-IASB).
