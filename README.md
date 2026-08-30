# Orbital-Datacenters
Where you can actually put a datacenter in LEO, and what it costs you.

SpaceX's January 2026 FCC filing proposes an orbital compute constellation spanning 500–2,000 km, in sun-synchronous shells for continuously-sunlit baseline compute and ~30° shells for demand-peak capacity. This is an independent study of what that altitude range costs in radiation, and whether the "over 99% sunlit" premise holds across it. Every input traces to public filings, published models, or datasheets.

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

The model is validated against the SPENVIS orbit generator across all twelve cases:
inclination agrees to within 0.003°, and β to 0.014° for the sun-synchronous cases.
The larger residuals for the 30° cases are traced to J2 nodal regression, which the
analytic model neglects by construction.

Full derivation, limitations, and validation: **[docs/eclipse_geometry.md](docs/eclipse_geometry.md)**
