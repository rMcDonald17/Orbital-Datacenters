# Orbital-Datacenters
Where you can actually put a datacenter in LEO, and what it costs you.

SpaceX's January 2026 FCC filing proposes an orbital compute constellation spanning 500–2,000 km, in sun-synchronous shells for continuously-sunlit baseline compute and ~30° shells for demand-peak capacity. This is an independent study of what that altitude range costs in radiation, and whether the "over 99% sunlit" premise holds across it. Every input traces to public filings, published models, or datasheets.

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
| Solar particle events | fully exposed over the poles | **completely immune** |
| Eclipse | none above ~1,400 km | eclipses year-round — no eclipse-free season |

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

![Eclipse geometry vs altitude](figures/beta_eclipse_sweep.png)

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

**The 30° shells are not analysed here.** Their RAAN is not frozen — it regresses
6.6°/day at 500 km, cycling the orbit plane through the full range of sun angles every
~47 days rather than following the seasons. Because β never reaches the eclipse-free
threshold at 30° inclination, those shells eclipse on essentially every revolution
year-round, and their sunlit fraction is a time-average over the nodal cycle rather than
a seasonal curve. That calculation is deferred; the >99%-sunlit premise in the filing
applies to the sun-synchronous shells regardless.

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

## Where that leaves the trade

The two constraints bound the useful band from opposite directions. Eclipse geometry
rules out the bottom of the filed range for the sun-synchronous shells; dose scaling and
shield mass make the top expensive. Roughly **1,200–1,500 km** satisfies both — 99.4–100%
sunlit, 3.3–7.8 krad(Si)/yr at 5 mm, comfortably inside a 100 krad five-year budget.

The filed 2,000 km ceiling buys 0.6 percentage points of additional sunlight for 2.75× the
annual dose relative to 1,500 km.
