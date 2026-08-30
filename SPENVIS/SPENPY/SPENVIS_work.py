

# %% ---- long table: load from disk if it exists, rebuild only when asked
from pathlib import Path
import pandas as pd

ROOT = Path(r"C:\Projects\Orbital_Data\SPENVIS")
OUT = Path("environment"); OUT.mkdir(exist_ok=True)
REBUILD = False          # flip to True after touching spenvis_io / adapters

def find_table(stem="odrs_long", out=OUT):
    for ext in ("parquet", "csv"):
        p = out / f"{stem}.{ext}"
        if p.exists():
            return p
    return None

def save_table(df, stem):
    try:
        p = OUT / f"{stem}.parquet"; df.to_parquet(p); return p
    except ImportError:
        p = OUT / f"{stem}.csv"; df.to_csv(p, index=False); return p

main = find_table()
if REBUILD or main is None:
    from spenvis_adapters import build_long
    long, wide = build_long(ROOT)
    main = save_table(long, "odrs_long")
    wide.to_csv(OUT / "_unmelted_blocks.csv", index=False)
    for q, g in long.groupby("quantity"):
        g.to_csv(OUT / f"{q}.csv", index=False)
    print(f"rebuilt -> {main} ({len(long):,} rows, {len(wide)} unmelted)")
else:
    from plots import load_long
    long = load_long(main)
    wide = pd.read_csv(OUT / "_unmelted_blocks.csv")
    newest = max((p.stat().st_mtime for p in ROOT.rglob("spenvis_*.txt")), default=0)
    stale = newest > main.stat().st_mtime
    print(f"reusing {main} ({len(long):,} rows)"
          + ("   [STALE: SPENVIS files are newer]" if stale else ""))

# %%

from plots import load_long, dose_table, plot_tid_vs_altitude, plot_dose_decomposition

dose = dose_table(load_long("environment/odrs_long.parquet"))
plot_tid_vs_altitude(dose)

# %%
plot_dose_decomposition(dose, cases=["A1", "A6", "B1", "B6"], cycle="min")
plot_dose_decomposition(dose, cases=["A1", "A6", "B1", "B6"], cycle="max")

#  %%

from plots import plot_beta_sweep, plot_beta_validation

summary = plot_beta_sweep(alts=[500, 700, 1000, 1200, 1500, 2000],
                          ltan_hr=6.0, inc_b=30.0, raan_b=0.0, year=2026)
print(summary.to_string(index=False))

# independent check against SPENVIS orbit output
val = plot_beta_validation(r"C:\Projects\Orbital_Data\SPENVIS")
# %%
