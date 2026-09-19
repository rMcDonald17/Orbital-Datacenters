# Radiator Fin Conduction Sub-Model

**Status:** Thermal aside, Part C Layer 3 — deliverable T3
**Scope:** 30° shell, 800 km, zenith/deployed radiator panels
**Code:** [`conduction.py`](../THERMPY/conduction.py)

---

> **This sub-model requires one iteration once the panel and bus geometry are defined.**
> Every absorbed-flux figure below treats each radiating surface as an isolated flat
> plate seeing space and Earth only. No inter-surface view factors are included: no
> array-to-radiator radiative coupling, no bus blockage of the radiator's view to space,
> no panel-to-panel exchange. The array coupling alone is worth between 0% and +161% of
> the radiator's entire environmental load depending on boom layout (§7), which is larger
> than any other correction in this document. **The fin formulation and the fin
> optimisation are unaffected — only the value of $q_{\text{abs}}$ fed into them changes.**
> Re-run this sub-model when the configuration is fixed.

---

## 1. What this sub-model produces

Three things the rest of the study needs:

1. **Fin efficiency $\eta_{\text{net}}$** — the factor by which a real radiator underperforms
   an isothermal plate at root temperature. Part A implicitly assumed $\eta = 1$.
2. **Radiator areal mass $\mu_{\text{rad}}$** — derived from a fin geometry optimisation,
   replacing the assumed `radiator_kg_per_m2 = 10.0` in `sizing.MassModel`.
3. **The conductor value** feeding node N4 in the Layer 2 network (Part C).

---

## 2. Governing equation

Energy balance on a differential slice of the sandwich panel, unit width, thickness $dx$
along the fin. Two facesheets of thickness $t_f$ carry all conduction; both outer surfaces
radiate and both absorb environmental flux.

$$\frac{d}{dx}\left(2 k t_f \frac{dT}{dx}\right) dx = \left[\varepsilon\sigma\left(T^4 - T_{sp}^4\right) - q_A\right] dx + \left[\varepsilon\sigma\left(T^4 - T_{sp}^4\right) - q_B\right] dx$$

where $q_A$ and $q_B$ are the absorbed environmental fluxes on the two faces. Dividing by 2:

$$k t_f \frac{d^2T}{dx^2} = \varepsilon\sigma\left(T^4 - T_{sp}^4\right) - q_{\text{abs}}, \qquad q_{\text{abs}} \equiv \frac{q_A + q_B}{2}$$

**The two faces enter only through their mean.** A panel with an illuminated face and a
shaded face behaves identically to one with both faces at the average — the facesheets are
isothermal through thickness (§5) so the two loads simply sum into one slice. This is why
the asymmetric edge-on configuration is treated with a single $q_{\text{abs}}$.

![Differential energy balance on a fin slice](figures/fin_energy_balance.svg)

*(a) The fin spans from one embedded heat pipe to the next. The solved domain runs from the
pipe contact at $x = 0$ to the symmetry plane at $x = L$; everything beyond is the mirror
image. The Dirichlet condition fixes temperature at the pipe contact, the Neumann condition
sets zero gradient at the symmetry plane, and the sketched profile $T(x)$ shows why — the
minimum sits at the plane, so the gradient there vanishes by construction.
(b) A slice of width $dx$ cut from that domain, with each term of the balance mapped to the
flux it represents. Blue: conduction along the facesheets. Coral: radiation leaving both
outer surfaces. Amber: environmental flux absorbed on both. Dividing through by $2\,dx$
gives the working form.*

This is Fourier conduction with a nonlinear distributed sink. It is the classic fin
equation with $h(T - T_\infty)$ replaced by $\varepsilon\sigma(T^4 - T_{sp}^4)$, and that
replacement is what removes the closed-form solution.

### Boundary conditions

$$T(0) = T_{\text{root}}, \qquad \left.\frac{dT}{dx}\right|_{x=L} = 0$$

$L$ is the **half**-spacing between heat pipes, so $x = L$ is a symmetry plane midway
between two pipes. The adiabatic tip is therefore exact, not an approximation — conditions
are mirror-identical either side, so no heat can cross. The domain is one half-span,
solved once and tiled across the panel.

The root condition assumes the heat pipe is isothermal along its length. See §8.

### Fin efficiency

Defined conventionally — net heat actually removed, over net removed by an isothermal fin
at root temperature:

$$\eta_{\text{net}} = \frac{\int_0^L \left[\varepsilon\sigma\left(T^4 - T_{sp}^4\right) - q_{\text{abs}}\right] dx}{\left[\varepsilon\sigma\left(T_{\text{root}}^4 - T_{sp}^4\right) - q_{\text{abs}}\right] L}$$

Note the denominator carries $q_{\text{abs}}$ too. This matters — see §4.

---

## 3. Solution methods and validation

### 3.1 Linearised

Writing $\varepsilon\sigma(T^4 - T_{sp}^4) \approx h_{\text{rad}}(T - T_{sp})$ with
$h_{\text{rad}} = 4\varepsilon\sigma\bar{T}^3$ recovers the standard fin solution:

$$\eta_{\text{lin}} = \frac{\tanh(mL)}{mL}, \qquad m = \sqrt{\frac{h_{\text{rad}}}{k t_f}}$$

$m$ is the **fin parameter**, with units of m$^{-1}$: $h_{\text{rad}}$ is W m$^{-2}$ K$^{-1}$
and $k t_f$ is W K$^{-1}$, so their ratio is m$^{-2}$. Its reciprocal $1/m$ is the
**conduction length** — the distance heat travels through the facesheets before radiation
drains it. The product $mL$ is therefore **dimensionless**: the fin half-length measured in
conduction lengths. It is the single parameter governing fin behaviour, and every result in
§4 and §9 collapses onto it. For the Al 6061 optimum, $1/m = 82$ mm and $L = 88$ mm, so
$mL = 1.07$.

### 3.2 Nonlinear

Solved as a two-point **boundary value problem** (BVP) — an ODE with conditions imposed at
both ends of the domain (here $T$ at the root and $dT/dx$ at the tip) rather than both at
one end, so it cannot be marched forward and must be solved across the whole domain at once.
Collocation is used (`scipy.integrate.solve_bvp`). A uniform initial guess
at $T_{\text{root}}$ exceeds the mesh node limit without converging, because it places the
whole domain far from the solution and the solver refines the mesh rather than moving the
profile. The linearised solution of §3.1,
$T(x) = T_{eq} + (T_{\text{root}} - T_{eq})\cosh\!\left[m(L-x)\right]/\cosh(mL)$,
is used as the starting iterate instead. This is a change of *initial guess*, not of method
— §3.3 below is a separate solver written independently to check this one.

### 3.3 Finite-difference cross-check

An independent second solver, so that a wrong answer from collocation would not go
unnoticed: the same BVP discretised with second-order central differences and solved by
Newton's method with an analytic Jacobian. Newton is used here because the first attempt at
this check — Jacobi iteration — failed (see the note below), not because of anything in
§3.2. Agreement with the collocation result:

| $L$ (mm) | $\eta_{\text{bvp}}$ | $\eta_{\text{fd}}$ | Difference |
|---:|---:|---:|---:|
| 50 | 0.911521 | 0.911659 | −0.015% |
| 75 | 0.829701 | 0.829922 | −0.027% |
| 100 | 0.746292 | 0.746561 | −0.036% |
| 150 | 0.602160 | 0.602440 | −0.047% |
| 200 | 0.494432 | 0.494678 | −0.050% |

Grid convergence at $L = 100$ mm is monotonic toward the collocation value — 0.74849 at
$n = 51$, 0.74643 at $n = 801$, against 0.74629 from `solve_bvp` — which is the expected
first-order behaviour and confirms both methods solve the same problem.

*Note: the first version of this cross-check used Jacobi iteration and produced nonsense —
profiles that barely cooled below root temperature, disagreeing with §3.2 by up to 50%.
Jacobi convergence scales as $n^2$; at $n = 4001$ it needs roughly 16 million sweeps and was
given 200,000, so it returned an unconverged profile that looked like a physical result. The
Newton solve above is the valid comparison. Recorded because an unconverged iterative solve
fails silently — it returns numbers, not an error.*

---

## 4. When the linearisation fails

| $L$ (mm) | $mL$ | $\eta_{\text{lin}}$ | $\eta_{\text{nl}}$ | Error | Tip (°C) |
|---:|---:|---:|---:|---:|---:|
| 50 | 0.43 | 0.9424 | 0.9438 | −0.15% | 38.0 |
| 75 | 0.83 | 0.8184 | 0.8297 | −1.36% | 23.0 |
| 100 | 1.11 | 0.7236 | 0.7463 | −3.04% | 11.0 |
| 150 | 1.67 | 0.5582 | 0.6022 | −7.31% | −12.6 |
| 200 | 2.73 | 0.3639 | 0.4222 | −13.8% | −49.3 |
| 300 | 4.09 | 0.2445 | 0.2973 | −17.8% | −82.7 |

The error is **always negative** — linearising always under-predicts efficiency — and it
tracks $mL$ cleanly. The reason it is conservative rather than random: $h_{\text{rad}}$ is
evaluated at $\bar{T} = T_{\text{root}}$, over-estimating radiative loss everywhere
downstream where the fin is colder.

**Working rule: linearise for $mL < 1$ (error under 3%); solve the BVP beyond it.**

### The source term does not help efficiency

Adding $q_{\text{abs}}$ makes the fin substantially warmer — tip rises from 11.6 °C to
20.1 °C at $L = 90$ mm, and by 42 K at $L = 300$ mm — yet $\eta_{\text{net}}$ gets slightly
**worse**:

| $L$ (mm) | $\eta$, no source | $\eta_{\text{net}}$, $q_{\text{abs}} = 122$ | Change |
|---:|---:|---:|---:|
| 50 | 0.8967 | 0.8956 | −0.1% |
| 90 | 0.7509 | 0.7453 | −0.7% |
| 150 | 0.5678 | 0.5543 | −2.4% |
| 300 | 0.3287 | 0.3084 | −6.2% |

Efficiency is a ratio, and the source shrinks the denominator faster than it helps the
numerator. The reference quantity becomes
$\varepsilon\sigma T_{\text{root}}^4 - q_{\text{abs}} = 493 - 122 = 371$ W/m² instead of
the full 493, so the same absolute droop is a larger *fractional* loss.

### What the source term does reveal

A hard ceiling on useful fin length that does not exist without it. With
$q_{\text{abs}} = 122$ W/m², radiative equilibrium sits at **−48.8 °C** — a fin approaching
that temperature removes *zero* net heat, so beyond some length the panel is mass that does
nothing. Without the source the fin contributes however long it grows, merely with
diminishing returns. This ceiling is what actually bounds fin length on a sunlit radiator.

---

## 5. Assumptions

| Assumption | Justification | Confidence |
|---|---|---|
| 1-D through thickness | $\text{Bi} = h_{\text{rad}} t_f / k = 9\times10^{-6}$ | Safe by 5 orders of magnitude |
| 1-D along the pipe | 66:1 fin aspect ratio; 2-D solve gives +0.25% (§6) | Strong |
| Core carries no heat | Honeycomb in-plane conductivity negligible vs facesheets | Conservative |
| Isothermal root | Heat pipes near-isothermal by construction | See §8 |
| Steady state | Fin $\tau \approx 8$ min vs 109 min orbit | Holds for hot case; **fails for the cold-case transient** |
| Uniform $q_{\text{abs}}$ | No shadowing modelled | See §7 |
| Both faces at mean flux | Exact for the 1-D slice balance (§2) | Exact |

---

## 6. Dimensionality: 1-D is sufficient

Heat pipes run the full panel length, so each fin is 88 mm across and ~5.8 m long — a
**66:1 aspect ratio**. There is no gradient to resolve along the pipe: it is isothermal by
construction, the environment is uniform along it, and the far edges are adiabatic. Only
the strip within one conduction length (82 mm) of a pipe *end* is genuinely 2-D, which is
**2.8% of panel area**.

A 2-D finite-difference solve of a panel strip against the 1-D result:

| Strip length $L_y$ | 2-D net (W/m²) | 1-D net (W/m²) | Difference |
|---:|---:|---:|---:|
| 88 mm | 302.2 | 301.5 | +0.25% |
| 250 mm | 302.2 | 301.5 | +0.25% |
| 1,000 mm | 302.2 | 301.5 | +0.25% |
| 2,900 mm | 302.2 | 301.5 | +0.25% |

The residual is constant across a 33× range of strip lengths — that is a discretisation
floor, not an edge effect. **2-D buys nothing here.**

Where 2-D *would* matter, in descending order of magnitude:

1. **Partial shadowing** by the bus or opposite panel. Shadow boundaries do not align with
   the fin direction. The fix is a 2-D flux map feeding a set of 1-D fins with local
   $q_{\text{abs}}$, not a 2-D fin.
2. **Heat pipe axial gradient** (§8) — ~15× larger than the geometric 2-D effect.
3. **Manifold and header regions** — genuinely 2-D, small area, properly a
   spreading-resistance problem.

---

## 7. Inter-surface view factors — the missing term

Fins within a panel are **coplanar**, so their mutual view factor is exactly zero. Same for
two deployed panels lying in the same plane. Nothing is missing there.

What is missing is coupling between different structures. A solar array's back face runs at
~33 °C and emits 423 W/m². Its coupling to a radiator depends entirely on boom layout:

| Array–radiator separation | $F$ | Absorbed | vs. 93 W/m² baseline |
|---:|---:|---:|---:|
| 3 m | 0.415 | 149 W/m² | +161% |
| 6 m | 0.200 | 72 W/m² | +77% |
| 10 m | 0.093 | 34 W/m² | +36% |
| 20 m | 0.027 | 10 W/m² | +10% |
| Coplanar | 0 | 0 | — |

Sun-tracking pushes toward the bad case: at high $\beta$ the array normal points along orbit
normal, and so do the radiator normals, so the two are parallel at exactly the $\beta$ where
the radiator is already working hardest. Parallel plates either face each other or sit
coplanar, and which one occurs is purely a matter of boom placement.

**Design rule falling out of this: keep the array coplanar with the radiators, laterally
offset, not stacked axially in front of them.** It costs nothing structurally and is worth
up to 77% of the radiator's environmental budget.

Bus blockage removes solid angle from the radiator's view of space — a loss, not a gain —
but the bus is small relative to 6 m panels and mostly occludes grazing directions that
contribute little. Not yet quantified.

---

## 8. Heat pipe assumption

$T_{\text{root}}$ is taken constant along the pipe. Real axially grooved heat pipes hold a
few kelvin end to end, and considerably more near their capillary or entrainment limits. A
3 K droop over 5.8 m changes local rejection by roughly 3.8% — **fifteen times the 2-D
geometric effect of §6**, and larger than every numerical error in this document.

No heat pipe transport model exists in this study. Capillary limit, entrainment limit,
sonic limit and boiling limit are all unchecked against the transported power. This is the
largest unmodelled effect inside the radiator itself, as distinct from §7 which is the
largest unmodelled effect outside it.

---

## 9. Panel model and mass optimisation

Aluminium honeycomb sandwich: two facesheets, core, embedded heat pipes at spacing $2L$,
radiating from both faces.

$$\mu_{\text{rad}} = 2\rho t_f + \rho_{\text{core}} h_{\text{core}} + \frac{\mu_{\text{hp}}}{2L} + \mu_{\text{misc}}$$

Defaults: core 15 mm at 50 kg/m³, heat pipes 0.35 kg/m (axially grooved Al/ammonia),
misc 0.5 kg/m² (coating, adhesive, doublers).

Objective: minimise kg per kW rejected, where heat per m² of **physical** panel is

$$\dot{q}_{\text{panel}} = 2\,\eta_{\text{net}}\left[\varepsilon\sigma\left(T_{\text{root}}^4 - T_{sp}^4\right) - q_{\text{abs}}\right]$$

Optimised over $(t_f, L)$ at $T_{\text{root}} = 318$ K, $q_{\text{abs}} = 122$ W/m²:

| Material | $k$ (W/m·K) | $t_f$ (mm) | $L$ (mm) | $\eta_{\text{net}}$ | kg/m² | W/m² panel | kg/kW |
|---|---:|---:|---:|---:|---:|---:|---:|
| Al 6061-T6 | 167 | 0.250 | 87.8 | 0.754 | 4.59 | 559 | 8.22 |
| Al 2024-T3 | 121 | 0.250 | 78.4 | 0.738 | 4.87 | 547 | 8.91 |
| Al 1100 | 222 | 0.250 | 96.5 | 0.769 | 4.42 | 570 | 7.75 |
| K1100 composite | 500 | 0.250 | 135.0 | 0.790 | **3.45** | 586 | **5.88** |

### Reading the optimum

**Facesheet thickness hits the 0.25 mm floor in every case.** Thinner is always lighter; the
binding constraint is handling and manufacturability, not physics. The optimum sits on a
constraint boundary, so sourcing thinner facesheets continues to reduce mass.

**The optimum sits at $mL \approx 1.07$.** For Al 6061, $k t_f = 0.042$ W/K and
$h_{\text{rad}} = 6.20$ W/m²K, giving a conduction length $1/m = 82$ mm against an optimal
$L$ of 88 mm. The design statement: **make the fin about as long as its conduction length,
and no longer.** Shorter and you carry a heat pipe every few centimetres; longer and the tip
approaches the −48.8 °C equilibrium where it stops removing heat.

This is also, not coincidentally, exactly where the linearisation starts to degrade (§4) —
both are set by the same conduction-versus-radiation balance.

**The optimum is shallow in $L$** (6.6 to 7.0 kg/kW across 80–125 mm) while $\eta$ varies
from 0.74 to 0.58 over the same span. There is real freedom to trade heat-pipe count
against fin performance at nearly constant mass.

**K1100 pitch-fibre composite wins by 25%** on mass, from higher in-plane conductivity at
lower density, and permits 54% wider pipe spacing — fewer heat pipes, fewer joints, fewer
leak paths. Not costed or qualified here.

---

## 10. Impact on Part A

Part A used 371 W/m² net rejection at $T_{\text{rad}} = 318$ K with **no fin efficiency
applied** — an implicit $\eta = 1$. With $\eta_{\text{net}} = 0.754$ the figure is 280 W/m².

- **Radiator area understated by 33%.**
- Part A treated 45 °C as the *effective* radiating temperature. With a real fin, a 45 °C
  root gives an effective temperature of **22.8 °C**.
- Achieving 45 °C effective would need a **68.5 °C root**, pushing coolant outlet toward
  80 °C and breaking the 85 °C junction limit.

Working against this, the derived panel areal mass is **4.59 kg/m²** rather than the assumed
10 — but that assumption covered the *assembly*. This sub-model computes panel plus heat
pipes only, excluding deployment mechanism, hinges, manifolds, pump, accumulator and fluid
inventory. The assumption is therefore **decomposed, not replaced**: 4.6 kg/m² derived,
~5.4 kg/m² still unsourced.

Net effect on vehicle mass is not obvious without re-running, and it feeds back into $A/m$,
drag, and the altitude selection. **Not yet propagated** — see the banner at the top; this
should be done in one pass together with the face-counting correction, not two.

---

## 11. Open items

1. **Re-run with real view factors** once panel and bus geometry are defined (§7). Largest
   single uncertainty.
2. **Heat pipe transport model** — axial gradient and operating limits (§8).
3. **Source thinner facesheets** — the optimum is on a manufacturing constraint, not a
   physical one (§9).
4. **Non-panel radiator mass** — deployment, manifolds, pump, fluid. ~5.4 kg/m² unsourced.
5. **Cold-case transient** — steady state fails when the payload is off in eclipse. The fin
   will need 2–3 nodes across its span in Layer 2, since a single lumped node would report a
   temperature that exists nowhere on the panel.
6. **K1100 cost and qualification** if the 25% mass saving is to be claimed.
