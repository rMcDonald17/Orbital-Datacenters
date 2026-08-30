
from spenvis_io import read_spenvis

"""spenvis_index.py"""
import re
from pathlib import Path
import pandas as pd

ALTS  = [500, 700, 1000, 1200, 1500, 2000]
_SET  = re.compile(r"^([AB])_cases$", re.I)
_ALT  = re.compile(r"^(\d+)\s*km$", re.I)
_STEP = re.compile(r"^(\d+)[_\s-]+(.*)$")

def parse_path(p: Path, root: Path):
    parts = list(p.relative_to(root).parts)
    out = {"rel_path": str(p.relative_to(root)), "file": p.name}
    ai = next((i for i, s in enumerate(parts) if _ALT.match(s)), None)
    if ai is None:
        raise ValueError(f"no <NNN>km directory in {out['rel_path']}")
    alt = int(_ALT.match(parts[ai]).group(1))
    si = next((s for s in parts[:ai] if _SET.match(s)), None)
    if si is None:
        raise ValueError(f"no <X>_cases directory in {out['rel_path']}")
    letter = _SET.match(si).group(1).upper()
    m = _STEP.match(parts[ai + 1]) if len(parts) > ai + 1 else None
    if not m:
        raise ValueError(f"no <N>_<tag> step directory in {out['rel_path']}")
    out.update(case=f"{letter}{ALTS.index(alt)+1}" if alt in ALTS else f"{letter}?{alt}",
               set=letter, alt_km=alt, step=int(m.group(1)), step_tag=m.group(2),
               subpath="/".join(parts[ai + 2:-1]))
    return out

def check_step_order(df):
    bad = []
    for case, g in df.groupby("case"):
        t = g.groupby("step").run_time.min().sort_index()
        if not t.is_monotonic_increasing:
            off = [f"  step {s}: {v}" for s, v in t.items()]
            bad.append(f"{case}: run times out of step order\n" + "\n".join(off))
    return bad

def build_manifest(root, pattern="spenvis_*.txt"):
    root = Path(root); rows, errs = [], []
    for p in sorted(root.rglob(pattern)):
        try:
            base = parse_path(p, root)
        except ValueError as e:
            errs.append(str(e)); continue
        for b in read_spenvis(p):
            m = b.meta
            rows.append({**base, "block": b.block_index, "mod_abb": b.model,
                         "trp_mod": m.get("TRP_MOD"), "species": m.get("SPECIES"),
                         "plt_hdr": m.get("PLT_HDR"), "orb_typ": m.get("ORB_TYP"),
                         "orb_apo": m.get("ORB_APO"), "orb_per": m.get("ORB_PER"),
                         "orb_inc": m.get("ORB_INC"), "orb_raa": m.get("ORB_RAA"),
                         "mis_dur": m.get("MIS_DUR"), "prj": m.get("PRJ_DEF"),
                         "version": b.version, "run_time": b.run_time,
                         "n_rows": len(b.data), "columns": ";".join(b.columns)})
    return pd.DataFrame(rows), errs