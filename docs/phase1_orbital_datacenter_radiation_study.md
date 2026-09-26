# Orbital Datacenter Feasibility Study — Phase 1 Work Plan
## Radiation Environment vs. Altitude for a SpaceX-Class Compute Constellation

**Author:** rMcDonald17 · **Status:** Phase 1 (environment definition & altitude trade)
**Reference scenario:** SpaceX FCC filing, Jan 30, 2026 (public record) — orbital data center constellation, LEO

---

## 1. Reference Scenario (fixed inputs, from public filings/statements)

These parameters anchor the study to what SpaceX has actually proposed. Cite the FCC filing and press coverage in the report; do not represent any of this as insider information.

| Parameter | Value | Source basis |
|---|---|---|
| Altitude range | 500–2,000 km | FCC filing |
| Orbit types | Sun-synchronous (>99% sunlit, baseline compute) + ~30° inclination (demand-peak capacity) | FCC filing |
| Bus heritage | Starlink V3-class flat-panel satellite, scaled up | Musk statements / filing context |
| Launch vehicle | Starship, ~100 t to LEO per flight | Public Starship specs |
| Comms | Optical intersatellite links (~1 Tbps class) relaying via Starlink; Ka-band TT&C | FCC filing |
| Constellation layout | Clusters ~50 km apart, serving different latency profiles | FCC filing coverage |
| Replacement cadence assumption | ~5 years (Starlink-like) | Analyst assumption — flag as such |
| Reference compute payload | H100 SXM-class accelerator + HBM (public datasheet; on-orbit precedent exists) | Analyst choice |

**Design question Phase 1 answers:** *What is the radiation cost of each altitude in the filed 500–2,000 km range, and does any part of that range make COTS GPU compute untenable?*

---

## 2. Phase 1 Scope

**In scope**
- Trapped proton/electron environment across the altitude sweep, both inclinations
- GCR and solar particle event (SPE) design environments
- Dose-depth curves behind representative shielding (slab model first, box model second)
- TID/year and displacement damage dose vs. altitude → implied electronics lifetime vs. 5-yr replacement cadence
- SEE-relevant flux metrics vs. altitude (>30 MeV proton flux, LET spectra behind shielding)
- South Atlantic Anomaly exposure fraction vs. altitude/inclination

**Out of scope until Phase 2**
- Device-level SEU rate predictions (needs cross-section data survey — Phase 2)
- Thermal/power/mass coupling
- Hardware testing

---

## 3. Orbit Case Matrix

Sweep altitude at two inclinations. The 1,200–2,000 km band clips the inner proton belt — expect orders-of-magnitude gradients there, so sample it more densely.

| Case | Altitude (km) | Inclination | Why |
|---|---|---|---|
| A1 | 500 | ~97.4° SSO | Starlink-like baseline, benign end |
| A2 | 700 | SSO | Common SSO operating band |
| A3 | 1,000 | SSO | Belt fringe begins |
| A4 | 1,200 | SSO | Inner belt onset — expect knee in dose curve |
| A5 | 1,500 | SSO | Deep in proton belt |
| A6 | 2,000 | SSO | Filed upper limit — worst case |
| B1–B6 | same six altitudes | 30° | Filed low-inclination shells; strong SAA sensitivity |

Orbit propagation: circular orbits, orbit-averaged fluxes over a full solar cycle position (run solar max and solar min variants; AP9/AE9 handles this via epoch settings).

---

## 4. Environment Models & Tools (all free)

| Model | Purpose | Access |
|---|---|---|
| **IRENE AP9/AE9** | Trapped protons/electrons (the modern standard; supersedes AP8/AE8) | Free from AFRL/VDL — request via the IRENE distribution site |
| **SPENVIS** | Web front-end: orbit generator, AP8/AE8 cross-check, SHIELDOSE-2 quick dose, NIEL/DDD | Free ESA account |
| **CREME96 / CREME-MC** | GCR spectra, LET spectra, SPE worst-day/worst-week models | Free registration (Vanderbilt) |
| **OMERE** | Alternative all-in-one (environment + dose + SEE) — good cross-check | Free from TRAD |
| **ESP/PSYCHIC** (in SPENVIS) | Statistical SPE fluence design case | via SPENVIS |

**Deliverable D1 — Environment database:** orbit-averaged differential + integral spectra (trapped p+/e-, GCR, design SPE) for all 12 cases, exported as CSV, committed to the repo with the exact model settings recorded in a `settings.md` per run. Reproducibility is part of the portfolio value.

---

## 5. Transport & Dose Analysis

**Step 1 — Slab screening (fast, week 1–2 of transport work).**
SHIELDOSE-2 (via SPENVIS) dose-depth curves for planar aluminum, 0.5–20 mm, all 12 cases. This gives the headline chart cheaply and validates the later transport runs.

**Step 2 — Full transport cross-check (your differentiator).**
Rebuild the slab problem in Geant4 (or MCNP/OpenMC if license/familiarity favors it) using the D1 spectra as sources:
- Aluminum slab, tally dose in a silicon layer at depths matching Step 1 → verify against SHIELDOSE-2 (expect agreement within ~tens of %; document discrepancies and why)
- Then a representative geometry: flat-panel bus, 2 mm Al face sheets, electronics bay behind 2 / 5 / 10 mm Al-equivalent, silicon dose tallies
- Add one graded-Z variant (Al/Ta/Al) at fixed areal mass to quantify whether graded shielding buys anything against trapped electrons at the high-altitude cases

**Step 3 — Displacement damage.**
NIEL-weighted proton fluence → DDD vs. altitude (SPENVIS provides NIEL curves). Feeds the Phase 2 solar array degradation model — note the cross-link.

**Deliverable D2 — Transport repo:** Geant4/MCNP input decks + run scripts + validation notebook comparing transport vs. SHIELDOSE-2.

---

## 6. Phase 1 Headline Products

1. **"Radiation cost of altitude" chart:** TID (krad-Si/year) vs. altitude, curves for 2/5/10 mm Al, both inclinations, solar max/min bands. *This is the money plot.*
2. **Lifetime map:** altitude bands where a 100 krad-tolerant COTS assumption survives the 5-year replacement cadence vs. bands where TID, not obsolescence, ends the satellite.
3. **SEE precursor metrics:** >30 MeV trapped proton flux and behind-shield LET spectra vs. altitude (Phase 2 input), plus SAA dwell fraction for the 30° shells.
4. **SPE ride-through note:** dose and flux transient for the design SPE; first-order statement on whether compute must power-safe.
5. **Report Section 1 (~10–15 pages):** written like a Phase-0 mission study section — assumptions table, methods, results, open risks.

---

## 7. Schedule (target: 6 weeks part-time)

| Week | Work |
|---|---|
| 1 | Tool access (IRENE request, SPENVIS + CREME accounts); repo skeleton; lock case matrix |
| 2 | Run all 12 environment cases; build D1 database; sanity-check vs. published Starlink-altitude values |
| 3 | SHIELDOSE-2 sweep; draft the headline chart from screening data |
| 4 | Geant4/MCNP slab validation vs. SHIELDOSE-2 |
| 5 | Box geometry + graded-Z variant; DDD calc |
| 6 | Write Report Section 1; publish repo + a short write-up post |

---

## 8. Repo Skeleton (github.com/rMcDonald17)

```
orbital-dc-radiation/
├── README.md              # project framing + headline chart up front
├── environment/           # D1: spectra CSVs + per-run settings.md
├── transport/
│   ├── slab_validation/   # Geant4/MCNP decks vs SHIELDOSE-2
│   └── bus_model/         # box geometry, graded-Z variant
├── notebooks/             # analysis + plots (Python: numpy/matplotlib)
├── report/                # Section 1 markdown/LaTeX
└── docs/assumptions.md    # the reference-scenario table, with sources
```

---

## 9. Risks & Guardrails

- **IRENE access lag:** the AP9/AE9 request can take days–weeks. Start week 1; use SPENVIS AP8/AE8 to prototype meanwhile and note the model swap.
- **Public-data discipline:** every reference-scenario number traces to the FCC filing, press coverage, or datasheets. Frame throughout as an *independent feasibility study using the filing as reference scenario*.
- **Scope creep:** no device SEU rates, no thermal, no economics in Phase 1. Park ideas in a `phase2_backlog.md`.
- **Validation honesty:** where transport disagrees with SHIELDOSE-2, investigate and write it up — the discrepancy discussion is a feature, not a bug.

---

## 10. Phase 2 Preview (for the backlog)

Device cross-section survey (HBM/DRAM/logic SEE literature) → SEU rates per satellite per orbit → constellation-scale silent-data-corruption estimate; solar array EOL power from Phase 1 DDD; radiator/power/mass coupling; TID bench test of a Jetson-class board against Phase 1 dose predictions.
