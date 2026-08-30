
from spenvis_io import read_spenvis
# inventory.py — put next to spenvis_io.py, run from that directory
import pandas as pd
from pathlib import Path

ROOT = Path(r"C:\Projects\Orbital_Data\SPENVIS")
rows = []
for p in sorted(ROOT.rglob("spenvis_*.txt")):
    rel = str(p.relative_to(ROOT))
    try:
        for b in read_spenvis(p):
            rows.append({"rel": rel, "suffix": p.stem.replace("spenvis_", ""),
                         "block": b.block_index, "mod_abb": b.model,
                         "plt_hdr": b.meta.get("PLT_HDR"),
                         "trp_mod": b.meta.get("TRP_MOD"),
                         "species": b.meta.get("SPECIES"),
                         "run_time": b.run_time, "n_rows": len(b.data),
                         "columns": ";".join(b.columns),
                         "col_units": ";".join(b.col_units),
                         "plt_typ": b.meta.get("PLT_TYP"),
                         "meta_keys": ";".join(sorted(b.meta)), "error": ""},
                         )
    except Exception as e:
        rows.append({"rel": rel, "suffix": p.stem.replace("spenvis_", ""), "error": repr(e)})
pd.DataFrame(rows).to_csv("inventory.csv", index=False)
print(len(rows), "blocks;", sum(1 for r in rows if r.get("error")), "errors")