"""
run_radiator_fin_submodel.py

Prints the tables that appear in docs/radiator_fin_submodel.md (thermal aside,
Part C Layer 3, deliverable T3). One driver per document -- the filename is the
document stem. Output is copied into the doc by hand; the doc is plain markdown.
Run top to bottom, or cell by cell.

q_abs comes from environment.panel_absorbed_flux, not a hardcoded number, so a
change of configuration or altitude propagates. Set COUPLING below when panel
and bus geometry are defined -- that is the iteration the doc's banner calls for.
"""
# %% ---- setup
import numpy as np
from scipy.optimize import minimize

import conduction as cd
import environment as env

ALT = 800.0                     # km, selected in docs/altitude_selection_30deg.md
CONFIG = 'deployed-edge-on'     # see the configuration table below
T_ROOT = 318.0                  # K, radiator root; selected in Part D
TF_FLOOR = 2.5e-4               # m, manufacturable facesheet minimum
L_BOUNDS = (0.02, 0.45)         # m, fin half-length search range

# THE ITERATION HOOK. Zero until panel and bus geometry are defined.
# coplanar -> 0.0 | 10 m separation -> 0.093 | 6 m -> 0.200 | 3 m -> 0.415
COUPLING = env.Coupling()

rule = lambda t: print("\n" + "=" * 80 + f"\n{t}\n" + "=" * 80)

def panel(**kw):
    kw.setdefault('t_root', T_ROOT)
    return cd.Panel.from_environment(ALT, config=CONFIG, coupling=COUPLING, **kw)

# %% ---- environment: where q_abs comes from
rule(f"PANEL CONFIGURATION at {ALT:.0f} km   (doc table, section 7)")
best, net, scores = env.best_panel_config(ALT, T_ROOT, coupling=COUPLING)
print(f"{'configuration':>22s}{'faces':>7s}{'per-face q':>16s}{'mean':>8s}{'W/m2 panel':>13s}")
worst = {}
for c, v in sorted(scores.items(), key=lambda t: -t[1]):
    worst[c] = max((env.panel_absorbed_flux(ALT, b, c, T_ROOT, coupling=COUPLING)
                    for b in np.linspace(0, env.beta_max(), 28)), key=lambda z: z.mean)
    faces = ', '.join(f'{q:.0f}' for q in worst[c].per_face)
    print(f"{c:>22s}{worst[c].n_faces:7d}{faces:>16s}{worst[c].mean:8.0f}{v:13.0f}")
print(f"\n  best: {best}   (this driver uses {CONFIG})")

pf = worst[CONFIG]
print(f"  mean fed to the fin: {pf.mean:.0f} W/m2")
if COUPLING.incident() > 0:
    print(f"  includes {env.OSR_BOL['eps_ir']*COUPLING.incident():.0f} W/m2 of inter-surface coupling")

# %% ---- validation: two independent solvers
rule("VALIDATION -- collocation vs finite-difference Newton")
print(f"{'L (mm)':>8s}{'eta_bvp':>11s}{'eta_fd':>11s}{'diff %':>10s}{'tip C':>9s}")
for L in [0.05, 0.075, 0.10, 0.15, 0.20]:
    p = panel(half_length=L)
    eb, tb = cd.fin_efficiency(p)
    ef, _, ok = cd.fd_fin_newton(p)
    print(f"{L*1e3:8.0f}{eb:11.6f}{ef:11.6f}{(eb-ef)/ef*100:+10.4f}{tb-273.15:9.2f}")

print("\ngrid convergence at L = 100 mm:")
p = panel(half_length=0.10)
for n in [51, 101, 201, 401, 801]:
    print(f"  n={n:4d}  eta={cd.fd_fin_newton(p, n)[0]:.8f}")
print(f"  bvp    eta={cd.fin_efficiency(p)[0]:.8f}")

# %% ---- where the linearization fails
# Fixed t_f so mL is the only variable. The collapse check below justifies
# quoting the rule in terms of mL at all.
rule("LINEARISATION ERROR vs mL   (t_f fixed)")
print(f"{'L (mm)':>8s}{'mL':>8s}{'eta_lin':>10s}{'eta_nl':>9s}{'err %':>9s}{'tip C':>9s}")
L_SWEEP = [0.02, 0.04, 0.06, 0.088, 0.12, 0.18, 0.25, 0.35]
for L in L_SWEEP:
    p = panel(half_length=L, t_face=TF_FLOOR)
    el, mL = cd.fin_efficiency_linear(p)
    en, tip = cd.fin_efficiency(p)
    print(f"{L*1e3:8.0f}{mL:8.3f}{el:10.4f}{en:9.4f}{(el-en)/en*100:+9.2f}{tip-273.15:9.1f}")
print(f"\n  t_f = {TF_FLOOR*1e3:.2f} mm throughout")
print("  rule: linearize for mL < 1 (error under 3%); solve the BVP beyond it")

# %% ---- collapse check: is mL really the only group?
# Non-dimensionalising with xi = x/L, theta = T/T_root gives
#     theta'' = (mL^2/4)[theta^4 - (Tsp/Tr)^4 - q_abs/(eps sig Tr^4)]
# so the solution depends on mL and on q_abs/(eps sig Tr^4) alone. The second
# group is fixed here, so every row at a given mL must agree exactly.
rule("COLLAPSE CHECK -- same mL reached by different (L, t_f)")
_p = panel()
g2 = pf.mean / (_p.eps * cd.SIGMA * T_ROOT ** 4)
print(f"  second group q_abs/(eps sig T_root^4) = {g2:.4f}, fixed\n")
print(f"{'target mL':>10s}{'t_f mm':>9s}{'L mm':>8s}{'mL':>8s}{'eta_lin':>10s}{'eta_nl':>9s}{'err %':>9s}{'tip C':>8s}")
for target in [0.60, 1.07, 2.00]:
    for tf in [2.0e-4, 3.0e-4, 5.0e-4, 1.0e-3]:
        L = target / np.sqrt(4 * _p.eps * cd.SIGMA * T_ROOT ** 3 / (_p.k * tf))
        p = panel(half_length=L, t_face=tf)
        el, mL = cd.fin_efficiency_linear(p)
        en, tip = cd.fin_efficiency(p)
        print(f"{target:10.2f}{tf*1e3:9.2f}{L*1e3:8.1f}{mL:8.3f}{el:10.4f}"
              f"{en:9.4f}{(el-en)/en*100:+9.2f}{tip-273.15:8.1f}")
    print()
print("  identical to 4 dp at each mL -> the collapse is exact, and the mL rule")
print("  is a statement about the governing equation, not an empirical fit")

# %% ---- what the source term does
rule("SOURCE TERM -- warmer fin, slightly worse NET efficiency")
print(f"{'L (mm)':>8s}{'mL':>8s}{'eta q=0':>10s}{'tip C':>8s}{'eta_net':>10s}{'tip C':>8s}{'change':>9s}")
for L in L_SWEEP:
    p1 = panel(half_length=L, t_face=TF_FLOOR)
    p0 = cd.Panel(half_length=L, t_face=TF_FLOOR, t_root=T_ROOT, q_abs=0.0)
    _, mL = cd.fin_efficiency_linear(p1)
    e0, t0 = cd.fin_efficiency(p0)
    e1, t1 = cd.fin_efficiency(p1)
    print(f"{L*1e3:8.0f}{mL:8.3f}{e0:10.4f}{t0-273.15:8.1f}{e1:10.4f}"
          f"{t1-273.15:8.1f}{(e1-e0)/e0*100:+8.1f}%")
print(f"\n  radiative equilibrium: {panel().t_equilibrium-273.15:.1f} C")
print("  a fin reaching this temperature removes zero net heat -- the length ceiling")

# Why does a warmer fin give WORSE net efficiency? Write eta = 1 - D/(G_root - q)
# with D the deficit lost to droop. The source cuts D, but it also shrinks the
# reference from G_root to G_root - q. Efficiency improves only if the deficit
# cut exceeds q/G_root. It is a near-tie, and it loses narrowly.
_g = 0.85 * cd.SIGMA * (T_ROOT ** 4 - cd.T_SPACE ** 4)
print(f"\n  DEFICIT DECOMPOSITION   (break-even = q/G_root = {pf.mean/_g*100:.2f}%)")
print(f"{'L mm':>6s}{'mL':>7s}{'D q=0':>9s}{'D q=98':>9s}{'cut':>9s}{'verdict':>9s}")
for L in [0.02, 0.06, 0.088, 0.18, 0.35]:
    p1 = panel(half_length=L, t_face=TF_FLOOR)
    p0 = cd.Panel(half_length=L, t_face=TF_FLOOR, t_root=T_ROOT, q_abs=0.0)
    e0, _ = cd.fin_efficiency(p0); e1, _ = cd.fin_efficiency(p1)
    _, mL = cd.fin_efficiency_linear(p1)
    d0 = _g - e0 * _g
    d1 = _g - (e1 * (_g - pf.mean) + pf.mean)
    cut = 1 - d1 / d0
    print(f"{L*1e3:6.0f}{mL:7.3f}{d0:9.1f}{d1:9.1f}{cut*100:8.2f}%"
          f"{('  drops' if cut < pf.mean/_g else '  rises'):>9s}")
print("  T^4 means a kelvin gained at the cold tip buys less emission than one")
print("  near the root, so the deficit removed is worth less than proportional")

# %% ---- dimensionality: does the fin need 2-D?   (doc table, section 6)
# Heat pipes run the full panel length, so each fin is ~88 mm across and metres
# long. Only the strip within one conduction length of a pipe END is really 2-D.
rule("2-D vs 1-D   (doc table, section 6)")
_p2 = panel(half_length=0.088, t_face=TF_FLOOR)
q1d = cd.net_per_area_1d(_p2)
Lc = 1.0 / np.sqrt(4 * _p2.eps * cd.SIGMA * T_ROOT ** 3 / (_p2.k * _p2.t_face))
print(f"{'span (mm)':>11s}{'2-D net':>10s}{'1-D net':>10s}{'diff':>9s}")
for span in [0.088, 0.25, 1.0, 2.9]:
    q2d, _ = cd.fd_fin_2d(_p2, span)
    print(f"{span*1e3:11.0f}{q2d:10.1f}{q1d:10.1f}{(q2d-q1d)/q1d*100:+8.3f}%")
print(f"\n  conduction length 1/m = {Lc*1e3:.0f} mm")
print(f"  fin aspect ratio at a 5.8 m panel: {5.8/_p2.half_length:.0f}:1")
print(f"  area within one conduction length of a pipe end: {2*Lc/5.8*100:.1f}%")
print("  constant across a 33x range of span -> a discretization floor, not an")
print("  edge effect. 2-D buys nothing here.")

# %% ---- mass optimization
rule("FIN MASS OPTIMISATION -- minimise kg/kW rejected")
nf = len(env.PANEL_CONFIGS[CONFIG])
print(f"{'material':>18s}{'t_f mm':>9s}{'L mm':>8s}{'mL':>7s}{'eta':>7s}{'kg/m2':>8s}{'W/m2':>8s}{'kg/kW':>8s}")
opt = {}
for mat in cd.MATERIALS:
    f = lambda v: cd.specific_mass(panel(material=mat, t_face=v[0], half_length=v[1]), nf)[0]
    r = minimize(f, [3e-4, 0.09], bounds=[(TF_FLOOR, 2e-3), L_BOUNDS], method="L-BFGS-B")
    p = panel(material=mat, t_face=r.x[0], half_length=r.x[1])
    sm, q, eta = cd.specific_mass(p, nf)
    _, mL = cd.fin_efficiency_linear(p)
    opt[mat] = p
    print(f"{mat:>18s}{p.t_face*1e3:9.3f}{p.half_length*1e3:8.1f}{mL:7.2f}"
          f"{eta:7.3f}{p.areal_mass:8.2f}{q:8.0f}{sm:8.2f}")
print(f"  t_f floored at {TF_FLOOR*1e3:.2f} mm -- a manufacturing constraint, not a physical optimum")
print("  note mL ~ 1 at every optimum: fin length tracks the conduction length 1/m")

# %% ---- where the areal mass actually comes from   (doc table, section 9)
best_p = opt['Al 6061-T6']
rule("AREAL MASS PROVENANCE -- how much of it is derived")
_rho = cd.MATERIALS[best_p.material]['rho']
_parts = [("facesheets 2*rho*t_f", 2*_rho*best_p.t_face, "DERIVED (t_f optimized, floored)"),
          ("core",                 best_p.core_density*best_p.core_height, "ASSUMED 15 mm @ 50 kg/m3"),
          ("heat pipes / 2L",      best_p.heatpipe_kg_per_m/(2*best_p.half_length), "ASSUMED 0.35 kg/m"),
          ("misc",                 best_p.misc_kg_per_m2, "ASSUMED coating/adhesive")]
_tot = best_p.areal_mass
for lab, v, src in _parts:
    print(f"  {lab:22s}{v:6.2f} kg/m2 {v/_tot*100:6.1f}%   {src}")
print(f"  {'TOTAL':22s}{_tot:6.2f} kg/m2   "
      f"assumed fraction {sum(v for _,v,t in _parts if t.startswith('ASSUMED'))/_tot*100:.0f}%")

# The heat-pipe rate is not merely additive: it sets the optimum.
print("\nSENSITIVITY TO THE HEAT-PIPE MASS RATE (largest assumed term)")
print(f"{'kg/m':>7s}{'opt L mm':>10s}{'eta':>8s}{'kg/m2':>8s}{'kg/kW':>8s}")
for hp in [0.15, 0.25, 0.35, 0.50, 0.75]:
    f = lambda v: cd.specific_mass(panel(material=best_p.material, t_face=v[0],
                                         half_length=v[1], heatpipe_kg_per_m=hp), nf)[0]
    r = minimize(f, [3e-4, 0.09], bounds=[(TF_FLOOR, 2e-3), L_BOUNDS], method="L-BFGS-B")
    pp = panel(material=best_p.material, t_face=r.x[0], half_length=r.x[1],
               heatpipe_kg_per_m=hp)
    sm, _, eta = cd.specific_mass(pp, nf)
    print(f"{hp:7.2f}{pp.half_length*1e3:10.1f}{eta:8.3f}{pp.areal_mass:8.2f}{sm:8.2f}")
print("  eta runs 0.84 to 0.69 and specific mass doubles across this range")

# %% ---- sensitivity about the optimum
rule(f"SENSITIVITY -- {best_p.material}, one parameter at a time")
print(f"{'L (mm)':>9s}{'eta':>8s}{'kg/m2':>9s}{'kg/kW':>9s}")
for L in [0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.20, 0.30]:
    p = panel(material=best_p.material, t_face=best_p.t_face, half_length=L)
    sm, _, eta = cd.specific_mass(p, nf)
    print(f"{L*1e3:9.0f}{eta:8.3f}{p.areal_mass:9.2f}{sm:9.2f}")
print(f"\n{'t_f (mm)':>9s}{'eta':>8s}{'kg/m2':>9s}{'kg/kW':>9s}")
for tf in [2.5e-4, 3e-4, 4e-4, 6e-4, 8e-4, 1.2e-3]:
    p = panel(material=best_p.material, t_face=tf, half_length=best_p.half_length)
    sm, _, eta = cd.specific_mass(p, nf)
    print(f"{tf*1e3:9.2f}{eta:8.3f}{p.areal_mass:9.2f}{sm:9.2f}")

# %% ---- coupling sensitivity: what the iteration will cost   (doc table, section 7)
# F is computed from the boom layout, not assumed. Panels are taken as 6 x 6 m
# directly opposed; a lateral offset or a coplanar arrangement gives F = 0.
PANEL_SIDE = 6.0
SEPARATIONS = [None, 20.0, 10.0, 6.0, 3.0]      # None = coplanar

t_arr = env.array_temperature()
rule("COUPLING SENSITIVITY -- the pending view-factor iteration")
print(f"  array back face at {t_arr-273.15:.1f} C, emitting "
      f"{0.85*cd.SIGMA*t_arr**4:.0f} W/m2 (env.array_temperature)")
print(f"  F from env.view_factor_parallel for {PANEL_SIDE:.0f} x {PANEL_SIDE:.0f} m "
      f"opposed plates\n")
print(f"{'separation':>12s}{'F_array':>9s}{'q_abs':>8s}{'eta':>8s}"
      f"{'W/m2 panel':>13s}{'kg/kW':>9s}{'penalty':>10s}")
base = None
for sep in SEPARATIONS:
    f = 0.0 if sep is None else env.view_factor_parallel(PANEL_SIDE, PANEL_SIDE, sep)
    p = cd.Panel.from_environment(ALT, config=CONFIG,
                                  coupling=env.Coupling(f_array=f, t_array=t_arr),
                                  material=best_p.material, t_face=best_p.t_face,
                                  half_length=best_p.half_length, t_root=T_ROOT)
    sm, q, eta = cd.specific_mass(p, nf)
    base = base or sm
    lab = 'coplanar' if sep is None else f'{sep:.0f} m'
    pen = '--' if sep is None else f'+{(sm/base-1)*100:.0f}%'
    print(f"{lab:>12s}{f:9.3f}{p.q_abs:8.0f}{eta:8.3f}{q:13.0f}{sm:9.2f}{pen:>10s}")
print("\n  eta barely moves: the coupling hurts by LOADING the panel, not by")
print("  degrading the fin, so the section 9 optimum survives the iteration")

# %% ---- what this does to Part A
rule("IMPACT ON PART A")
eta, _ = cd.fin_efficiency(best_p)
gross = best_p.eps * cd.SIGMA * (T_ROOT ** 4 - cd.T_SPACE ** 4)
print(f"  Part A net rejection (implicit eta=1):  {gross-best_p.q_abs:5.0f} W/m2 radiating")
print(f"  with fin efficiency eta={eta:.3f}:         {(gross-best_p.q_abs)*eta:5.0f} W/m2 radiating")
print(f"  radiator area understated by            {(1/eta-1)*100:5.0f}%")
print(f"  effective radiating temperature:        {T_ROOT*eta**0.25-273.15:5.1f} C")
print(f"  root needed for {T_ROOT-273.15:.0f} C effective:        {T_ROOT/eta**0.25-273.15:5.1f} C")
print(f"\n  derived panel areal mass: {best_p.areal_mass:.2f} kg/m2")
print(f"  sizing.MassModel assumes: 10.00 kg/m2  (assembly, not panel)")
print(f"  -> {10.0-best_p.areal_mass:.2f} kg/m2 of deployment, manifolds, pump, fluid: unsourced")

# %%

# %%
