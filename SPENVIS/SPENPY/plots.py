"""
ODRS / Datacenter Study -- Phase 1 figures.

    from plots import load_long, dose_table, plot_tid_vs_altitude, plot_dose_decomposition

    long = load_long("environment/odrs_long.parquet")
    dose = dose_table(long)
    plot_tid_vs_altitude(dose)
    plot_dose_decomposition(dose, cases=["A1", "A6"])

Run directly to regenerate every figure into figures/.
"""
from __future__ import annotations

import datetime as _dt
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import LogLocator, NullFormatter

FIGDIR = Path("figures")

COTS_KRAD = 100.0          # assumed COTS TID tolerance
LIFE_YEARS = 5.0           # Starlink-like replacement cadence
BUDGET = COTS_KRAD / LIFE_YEARS

DEPTHS_MM = [2.0, 5.0, 10.0]

SET_STYLE = {
    "A": ("SSO (~97.4 deg)", "#1f4e79"),
    "B": ("30 deg", "#b03a2e"),
}

# SP_worst deliberately absent: its Total is a mixed quantity
# (annual trapped dose + single-event solar dose), dropped upstream.
CYCLE = {
    "SHIELDOSE_MinMin_SP_total": "min",
    "SHIELDOSE_MaxMax_SP_total": "max",
}

COMPONENTS = [
    ("Trapped Protons", "#1f4e79", "-"),
    ("Electrons", "#b03a2e", "-"),
    ("Bremsstrahlung", "#7d6608", "-."),
    ("Solar Protons", "#117a65", "--"),
]


# ------------------------------------------------------------------ loading

def load_long(path):
    path = str(path)
    return pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path)


def dose_table(long):
    """All SHIELDOSE-2 rows, with solar-cycle state and krad/yr attached."""
    d = long[(long.quantity == "ionizing_dose")
             & (long.step_tag.isin(CYCLE))].copy()
    if d.empty:
        raise SystemExit("no ionizing_dose rows matched -- check the long table")
    d["cycle"] = d.step_tag.map(CYCLE)
    d["krad_yr"] = d.value / 1000.0
    return d


def _log_minor(ax):
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(LogLocator(base=10, subs="auto", numticks=20))
        axis.set_minor_formatter(NullFormatter())
    ax.grid(True, which="major", alpha=0.30)
    ax.grid(True, which="minor", alpha=0.12)


def _save(fig, stem):
    FIGDIR.mkdir(exist_ok=True)
    png = FIGDIR / f"{stem}.png"
    fig.savefig(png, dpi=200)
    fig.savefig(FIGDIR / f"{stem}.pdf")
    print("wrote", png)
    return fig


# ------------------------------------------------------- fig 1: TID vs alt

def plot_tid_vs_altitude(dose, depths=DEPTHS_MM, stem="tid_vs_altitude"):
    sel = dose[(dose.series == "Total") & (dose.depth_mm.isin(depths))]
    fig, axes = plt.subplots(1, len(depths), figsize=(13, 4.6),
                             sharey=True, constrained_layout=True)

    for ax, depth in zip(np.atleast_1d(axes), depths):
        sub = sel[sel.depth_mm == depth]
        for s, (label, colour) in SET_STYLE.items():
            w = (sub[sub.set == s]
                 .pivot_table(index="alt_km", columns="cycle",
                              values="krad_yr", aggfunc="first").sort_index())
            if w.empty:
                continue
            # MIN and MAX cross with depth (electron- vs proton-dominated):
            # an envelope, not a bound.
            if {"min", "max"} <= set(w.columns):
                ax.fill_between(w.index, w["min"], w["max"],
                                color=colour, alpha=0.18, lw=0)
            for cyc, ls in (("min", "-"), ("max", "--")):
                if cyc in w:
                    ax.plot(w.index, w[cyc], ls, color=colour, lw=1.8,
                            marker="o", ms=4, label=f"{label}, solar {cyc}")

        ax.axhline(BUDGET, color="0.35", lw=1.2, ls=":")
        ax.annotate(f"{COTS_KRAD:.0f} krad / {LIFE_YEARS:.0f} yr",
                    xy=(0.02, BUDGET), xycoords=("axes fraction", "data"),
                    va="bottom", fontsize=8, color="0.35")
        ax.set_yscale("log")
        ax.set_title(f"{depth:g} mm Al", fontsize=11)
        ax.set_xlabel("Altitude (km)")
        ax.grid(True, which="major", alpha=0.30)
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs="auto", numticks=20))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.grid(True, which="minor", alpha=0.12)

    np.atleast_1d(axes)[0].set_ylabel("TID (krad(Si)/yr)")
    np.atleast_1d(axes)[0].legend(fontsize=8, framealpha=0.9)
    fig.suptitle("Radiation cost of altitude: TID at centre of Al sphere, "
                 "AP-8/AE-8 + SAPPHIRE annual fluence", fontsize=12)
    return _save(fig, stem)


def plot_ddd_vs_altitude(long, depths=(2.0, 5.0, 10.0), stem="ddd_vs_altitude"):
    """Displacement damage dose vs altitude. AP-8 MIN only; protons only."""
    d = long[(long.quantity == "niel_dose_MeV_per_g")
             & (long.depth_mm.isin(depths))]
    if d.empty:
        raise SystemExit("no NIEL rows matched")

    fig, axes = plt.subplots(1, len(depths), figsize=(13, 4.6),
                             sharey=True, constrained_layout=True)
    for ax, depth in zip(np.atleast_1d(axes), depths):
        sub = d[d.depth_mm == depth]
        for s, (label, colour) in SET_STYLE.items():
            for series, ls, lw in (("Total", "-", 2.0),
                                   ("Trapped Protons", "--", 1.3),
                                   ("Solar Protons", ":", 1.3)):
                w = (sub[(sub.set == s) & (sub.series == series)]
                     .sort_values("alt_km"))
                if w.empty or not (w.value > 0).any():
                    continue
                ax.plot(w.alt_km, w.value.where(w.value > 0), ls, color=colour,
                        lw=lw, marker="o" if series == "Total" else None, ms=4,
                        label=f"{label}, {series.lower()}")
        pos = sub.value[sub.value > 0]
        if len(pos):
            ax.set_ylim(pos.min() / 3, pos.max() * 3)
        ax.set_yscale("log")
        ax.set_title(f"{depth:g} mm Al", fontsize=11)
        ax.set_xlabel("Altitude (km)")
        ax.grid(True, which="major", alpha=0.30)
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs="auto", numticks=20))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.grid(True, which="minor", alpha=0.12)

    np.atleast_1d(axes)[0].set_ylabel("DDD (MeV g$^{-1}$ yr$^{-1}$)")
    np.atleast_1d(axes)[0].legend(fontsize=7, framealpha=0.9)
    fig.suptitle("Displacement damage dose vs altitude "
                 "(NIEL, AP-8 MIN, protons only, centre of Al sphere)", fontsize=12)
    return _save(fig, stem)

# --------------------------------------------- fig 2: dose decomposition

def _crossover(x, a, b):
    """Depth where curves a and b swap rank, by log-log interpolation."""
    a = np.where(np.asarray(a, float) > 0, a, np.nan)
    b = np.where(np.asarray(b, float) > 0, b, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.log(a) - np.log(b)
    s = np.sign(d)
    idx = np.where(np.diff(s) != 0)[0]
    if not len(idx) or not np.all(np.isfinite(d[idx[0]:idx[0] + 2])):
        return None
    i = idx[0]
    lx = np.log(x[i:i + 2])
    return float(np.exp(lx[0] - d[i] * (lx[1] - lx[0]) / (d[i + 1] - d[i])))


def plot_dose_decomposition(dose, cases=("A1", "A6"), cycle="min",
                            stem="dose_decomposition"):
    """Dose-depth curves split by contributing species, one panel per case."""
    sel = dose[dose.cycle == cycle]
    cases = [c for c in cases if c in set(sel.case)]
    if not cases:
        raise SystemExit(f"none of the requested cases present for solar {cycle}")

    fig, axes = plt.subplots(1, len(cases),
                             figsize=(max(6.8, 5.2 * len(cases)), 4.8),
                             sharey=True, constrained_layout=True)

    for ax, case in zip(np.atleast_1d(axes), cases):
        sub = sel[sel.case == case]
        w = sub.pivot_table(index="depth_mm", columns="series",
                            values="krad_yr", aggfunc="first").sort_index()
        x = w.index.values

        if "Total" in w:
            ax.plot(x, w["Total"], color="0.15", lw=2.4, label="Total", zorder=5)
        for name, colour, ls in COMPONENTS:
            if name in w and w[name].gt(0).any():
                ax.plot(x, w[name].where(w[name] > 0), ls, color=colour,
                        lw=1.7, label=name)

        # the electron/proton swap is what makes the MIN/MAX envelope invert
        if {"Electrons", "Trapped Protons"} <= set(w.columns):
            xc = _crossover(x, w["Electrons"].values, w["Trapped Protons"].values)
            if xc:
                ax.axvline(xc, color="0.55", lw=1.0, ls=":")
                ax.annotate(f"e$^-$ / p$^+$ crossover\n{xc:.1f} mm",
                            xy=(xc, 0.06), xycoords=("data", "axes fraction"),
                            fontsize=8, color="0.4", ha="left",
                            xytext=(4, 0), textcoords="offset points")

        # AE-8 electron dose collapses ~10 decades by 20 mm; without a floor
        # the axis is all empty space. Show four decades below the Total.
        if "Total" in w:
            ax.set_ylim(w["Total"].min() / 1e4, w["Total"].max() * 2.5)

        alt = sub.alt_km.iloc[0]
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"{case} -- {alt:g} km, {'SSO' if case[0]=='A' else '30 deg'}",
                     fontsize=11)
        ax.set_xlabel("Al shield thickness (mm)")
        _log_minor(ax)

    np.atleast_1d(axes)[0].set_ylabel("Dose (krad(Si)/yr)")
    np.atleast_1d(axes)[0].legend(fontsize=8, framealpha=0.9, loc="upper right")
    fig.suptitle(f"Dose-depth decomposition by species (solar {cycle})",
                 fontsize=12)
    fig.supxlabel("AP-8/AE-8 + SAPPHIRE annual fluence, centre of Al sphere, "
                  "dose in Si", fontsize=8, color="0.4")
    return _save(fig, f"{stem}_solar{cycle}")


# ------------------------------------------- fig 3: beta angle / eclipse

MU_EARTH = 398600.4418      # km^3/s^2
R_EARTH = 6378.137          # km, equatorial
J2 = 1.08262668e-3
SSO_NODAL_RATE = 2 * np.pi / 365.2422 / 86400.0     # rad/s

ALTS_KM = [500, 700, 1000, 1200, 1500, 2000]


def sso_inclination(alt_km):
    """Inclination whose J2 nodal precession matches the mean solar rate.

    Reproduces SPENVIS's 97.4 deg at 500 km. Inclination RISES with altitude,
    which lowers beta and partly cancels the geometric benefit of going up.
    """
    a = R_EARTH + np.asarray(alt_km, float)
    n = np.sqrt(MU_EARTH / a ** 3)
    return np.degrees(np.arccos(-SSO_NODAL_RATE / (1.5 * J2 * (R_EARTH / a) ** 2 * n)))


def _days_from_j2000(doy, year):
    """Days from J2000.0 (2000-01-01 12:00 UT) to 00:00 UT on the given day."""
    epoch = _dt.datetime(2000, 1, 1, 12)
    base = (_dt.datetime(year, 1, 1) - epoch).total_seconds() / 86400.0
    return base + (np.asarray(doy, float) - 1.0)


def solar_declination(doy, year=2026):
    """Apparent solar declination, low-precision series (~0.2 deg)."""
    n = _days_from_j2000(doy, year)
    L = np.radians((280.460 + 0.9856474 * n) % 360)
    g = np.radians((357.528 + 0.9856003 * n) % 360)
    lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
    eps = np.radians(23.439 - 3.56e-7 * n)
    return np.degrees(np.arcsin(np.sin(eps) * np.sin(lam)))


def beta_angle(inc_deg, dec_deg, ltan_hr=6.0):
    """Sun/orbit-plane angle. For a frozen-LTAN SSO, RAAN - RA_sun is constant."""
    dR = np.radians((ltan_hr - 12.0) * 15.0)
    i, d = np.radians(inc_deg), np.radians(dec_deg)
    return np.degrees(np.arcsin(np.cos(d) * np.sin(i) * np.sin(dR)
                                + np.sin(d) * np.cos(i)))


def beta_critical(alt_km):
    """|beta| above which a circular orbit never enters the cylindrical umbra."""
    return np.degrees(np.arcsin(R_EARTH / (R_EARTH + np.asarray(alt_km, float))))


def eclipse_fraction(alt_km, beta_deg):
    """Fraction of each orbit in shadow (cylindrical umbra, circular orbit)."""
    h = np.asarray(alt_km, float)
    b = np.radians(np.abs(np.asarray(beta_deg, float)))
    arg = np.sqrt(h ** 2 + 2 * R_EARTH * h) / ((R_EARTH + h) * np.cos(b))
    return np.where(arg >= 1.0, 0.0, np.arccos(np.clip(arg, -1, 1)) / np.pi)

def solar_ra_dec(doy, year=2026):
    """Apparent solar RA and declination, low-precision series (~0.2 deg)."""
    n = _days_from_j2000(doy, year)
    L = np.radians((280.460 + 0.9856474 * n) % 360)
    g = np.radians((357.528 + 0.9856003 * n) % 360)
    lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
    eps = np.radians(23.439 - 3.56e-7 * n)
    return (np.degrees(np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))) % 360,
            np.degrees(np.arcsin(np.sin(eps) * np.sin(lam))))


def nodal_rate(alt_km, inc_deg):
    """Secular J2 nodal regression, deg/day. Negative for prograde orbits."""
    a = R_EARTH + np.asarray(alt_km, float)
    n = np.sqrt(MU_EARTH / a ** 3)
    return np.degrees(-1.5 * J2 * (R_EARTH / a) ** 2 * n
                      * np.cos(np.radians(inc_deg))) * 86400


def _beta_general(inc_deg, raan_deg, ra_deg, dec_deg):
    """Beta from the full expression, for orbits with no frozen sun angle."""
    dR = np.radians(raan_deg - ra_deg)
    i, d = np.radians(inc_deg), np.radians(dec_deg)
    return np.degrees(np.arcsin(np.cos(d) * np.sin(i) * np.sin(dR)
                                + np.sin(d) * np.cos(i)))

def _beta_family(fam, alts, t, ra, dec, ltan_hr, inc_b, raan_b, axb, axe):
    """Draw one orbit family onto a (beta, eclipse-fraction) axis pair."""
    cmap = plt.cm.viridis(np.linspace(0.05, 0.9, len(alts)))
    rows, free = [], []
    for alt, c in zip(alts, cmap):
        if fam == "A":
            inc = float(sso_inclination(alt))
            beta = np.abs(beta_angle(inc, dec, ltan_hr))
            lab = f"{alt} km (i={inc:.1f}\u00b0)"
        else:
            inc = inc_b
            raan = (raan_b + nodal_rate(alt, inc) * (t - 1)) % 360
            beta = np.abs(_beta_general(inc, raan, ra, dec))
            lab = f"{alt} km"
        bc = float(beta_critical(alt))
        frac = eclipse_fraction(alt, beta)
        lw = 1.6 if fam == "A" else 0.9

        axb.plot(t, beta, color=c, lw=lw, label=lab)
        axb.axhline(bc, color=c, lw=1.0, ls="--", alpha=0.8)
        # a flat zero line carries no information and drags the y-limits
        if frac.max() > 0:
            axe.plot(t, 100 * frac, color=c, lw=lw, label=f"{alt} km")
        else:
            free.append(alt)

        T = 2 * np.pi * np.sqrt((R_EARTH + alt) ** 3 / MU_EARTH) / 60
        rows.append(dict(family=fam, alt_km=alt, inc_deg=round(inc, 2),
                         beta_min=round(beta.min(), 2),
                         beta_max=round(beta.max(), 2),
                         beta_crit=round(bc, 2),
                         eclipse_free_days=round(float((frac == 0).sum()
                                                       * (t[1] - t[0])), 1),
                         max_eclipse_pct=round(100 * frac.max(), 2),
                         max_eclipse_min=round(frac.max() * T, 1),
                         period_min=round(T, 1),
                         annual_sunlit_pct=round(100 * (1 - frac.mean()), 2)))

    if free:
        axe.annotate("eclipse-free all year:\n" + ", ".join(f"{a} km" for a in free),
                     xy=(0.5, 0.55), xycoords="axes fraction",
                     ha="center", fontsize=10, color="0.35")
    axb.set_title("Sun/orbit-plane angle vs eclipse-free threshold", fontsize=11)
    axe.set_title("Eclipse fraction per orbit", fontsize=11)
    axb.set_ylabel("|beta| (deg)")
    axe.set_ylabel("Orbit in eclipse (%)")
    axb.annotate("dashed = threshold for that altitude", xy=(0.02, 0.03),
                 xycoords="axes fraction", fontsize=8, color="0.4")
    for a_ in (axb, axe):
        a_.set_xlabel("Day of year")
        a_.grid(alpha=0.3)
        a_.legend(fontsize=7.5, ncol=1, framealpha=0.9)
    return pd.DataFrame(rows)


def plot_beta_sweep_sso(alts=ALTS_KM, ltan_hr=6.0, year=2026, step=0.25,
                        stem="beta_eclipse_sweep_sso"):
    """Dawn-dusk sun-synchronous: beta follows the season."""
    t = np.arange(0, 365, step) + 1
    ra, dec = solar_ra_dec(t, year)
    fig, (axb, axe) = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    df = _beta_family("A", alts, t, ra, dec, ltan_hr, None, None, axb, axe)
    fig.suptitle(f"Eclipse geometry vs altitude, {year}: "
                 f"dawn-dusk SSO (LTAN {ltan_hr:g}h)", fontsize=12)
    _save(fig, stem)
    return df


def plot_beta_sweep_low_inc(alts=ALTS_KM, inc_b=30.0, raan_b=0.0, year=2026,
                            step=0.25, stem="beta_eclipse_sweep_lowinc"):
    """General low-inclination orbit: beta follows the nodal cycle."""
    t = np.arange(0, 365, step) + 1
    ra, dec = solar_ra_dec(t, year)
    fig, (axb, axe) = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    df = _beta_family("B", alts, t, ra, dec, None, inc_b, raan_b, axb, axe)
    fig.suptitle(f"Eclipse geometry vs altitude, {year}: "
                 f"{inc_b:g}\u00b0 inclination, RAAN {raan_b:g}\u00b0 at epoch",
                 fontsize=12)
    _save(fig, stem)
    return df


def plot_beta_sweep(**kw):
    """Both families. Returns the combined summary table."""
    sso = {k: v for k, v in kw.items() if k in ("alts", "ltan_hr", "year", "step")}
    low = {k: v for k, v in kw.items()
           if k in ("alts", "inc_b", "raan_b", "year", "step")}
    return pd.concat([plot_beta_sweep_sso(**sso), plot_beta_sweep_low_inc(**low)],
                     ignore_index=True)

def solar_ra_dec(doy, year=2026):
    """Apparent solar RA and declination, low-precision series (~0.2 deg)."""
    n = _days_from_j2000(doy, year)
    L = np.radians((280.460 + 0.9856474 * n) % 360)
    g = np.radians((357.528 + 0.9856003 * n) % 360)
    lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
    eps = np.radians(23.439 - 3.56e-7 * n)
    ra = np.degrees(np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))) % 360
    dec = np.degrees(np.arcsin(np.sin(eps) * np.sin(lam)))
    return ra, dec


def _read_att_beta(path, root):
    """Measured beta from an att file's orbit-frame Sun vector."""
    from spenvis_io import read_spenvis          # adjust to your layout
    from spenvis_index import parse_path
    b = read_spenvis(path)[0]
    m, d = b.meta, b.data
    sy = d["SatSun_Y"].values
    doy = (_dt.date(int(m["ORB_YEA"]), int(m["ORB_MON"]), int(m["ORB_DAY"]))
           - _dt.date(int(m["ORB_YEA"]), 1, 1)).days + 1
    return dict(**{k: parse_path(path, root)[k] for k in ("case", "set", "alt_km")},
                orb_typ=m["ORB_TYP"], inc_spenvis=m["ORB_INC"], raan=m["ORB_RAA"],
                year=int(m["ORB_YEA"]), doy=doy,
                beta_meas=np.degrees(np.arcsin(abs(sy.mean()))),
                sy_spread=float(sy.std()), n_pts=len(d))


def collect_validation(root, ltan_hr=6.0):
    root = Path(root); rows = []
    for p in sorted(root.rglob("spenvis_att.txt")):
        r = _read_att_beta(p, root)
        ra, dec = solar_ra_dec(r["doy"], r["year"])
        if r["orb_typ"] == "HEL":
            r["inc_model"] = float(sso_inclination(r["alt_km"]))
            r["beta_model"] = abs(float(beta_angle(r["inc_model"], dec, ltan_hr)))
        else:   # general orbit: use its own RAAN against the true solar RA
            r["inc_model"] = r["inc_spenvis"]
            dR = np.radians(r["raan"] - ra)
            i, dd = np.radians(r["inc_spenvis"]), np.radians(dec)
            r["beta_model"] = abs(np.degrees(np.arcsin(
                np.cos(dd) * np.sin(i) * np.sin(dR) + np.sin(dd) * np.cos(i))))
        r["d_inc"] = r["inc_spenvis"] - r["inc_model"]
        r["d_beta"] = r["beta_meas"] - r["beta_model"]
        rows.append(r)
    return pd.DataFrame(rows).sort_values(["set", "alt_km"])


def plot_beta_validation(root, ltan_hr=6.0, stem="beta_model_validation"):
    df = collect_validation(root, ltan_hr)
    a, b = df[df.set == "A"], df[df.set == "B"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)

    h = np.linspace(400, 2100, 200)
    ax1.plot(h, sso_inclination(h), "-", color="0.4", lw=1.6, label="analytic J2 model")
    ax1.plot(a.alt_km, a.inc_spenvis, "o", color="#1f4e79", ms=8,
             label="SPENVIS orbit generator")
    for _, r in a.iterrows():
        ax1.annotate(f"{r.d_inc:+.3f}\u00b0", (r.alt_km, r.inc_spenvis),
                     textcoords="offset points", xytext=(0, -16),
                     ha="center", fontsize=7.5, color="0.4")
    ax1.set_xlabel("Altitude (km)"); ax1.set_ylabel("Inclination (deg)")
    ax1.set_title("Sun-synchronous inclination", fontsize=11)
    ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

    for sub, colour, lab in ((a, "#1f4e79", "SSO (A-cases)"),
                             (b, "#b03a2e", "30 deg (B-cases)")):
        if sub.empty:
            continue
        ax2.plot(sub.alt_km, sub.beta_model, "-", color=colour, lw=1.4,
                 alpha=0.6, label=f"{lab}, model")
        ax2.plot(sub.alt_km, sub.beta_meas, "o", color=colour, ms=8,
                 label=f"{lab}, SPENVIS SatSun")
        dy = -16 if lab.startswith("SSO") else 12
        for _, r in sub.iterrows():
            ax2.annotate(f"{r.d_beta:+.3f}\u00b0", (r.alt_km, r.beta_meas),
                         textcoords="offset points", xytext=(0, dy),
                         ha="center", fontsize=7.5, color="0.4",
                         va="bottom" if dy > 0 else "top")
    ax2.set_xlabel("Altitude (km)"); ax2.set_ylabel("|beta| (deg)")
    ax2.set_title("Sun/orbit-plane angle, 1 Jan 2026", fontsize=11)
    ax2.legend(fontsize=8); ax2.grid(alpha=0.3)

    fig.suptitle("Analytic model vs SPENVIS orbit output (labels = model residual)",
                 fontsize=12)
    _save(fig, stem)
    return df