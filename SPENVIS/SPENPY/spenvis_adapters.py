


"""
Layer 2 for the ODRS / Datacenter Study SPENVIS pipeline.

Turns parsed Blocks (spenvis_io.read_spenvis) into one long-format table:

    case, set, alt_km, step, step_tag, quantity, model, species,
    index_name, index_value, index_units,
    series, depth_mm, energy_MeV, ion,
    value, units, run_time, source_file, block, flags

Design note: SPENVIS column *groups* already carry the structure we need
(name, unit, ncols, label-list), so there is one generic melter rather than
one adapter per file type. Per-suffix knowledge is confined to LABELS below.
"""
from __future__ import annotations

from spenvis_io import read_spenvis
from spenvis_index import parse_path

import re
import warnings
import pandas as pd

# ---------------------------------------------------------------- constants

# Group names that act as the independent variable rather than a measurement.
INDEX_NAMES = {"Energy", "Thick", "B", "L", "MJD", "Layer"}

# Elemental symbols H..U as they appear in GCR / SEP ion tables.
_ELEMENTS = set(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb "
    "Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re "
    "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U".split()
)

# suffix -> (quantity, note).  Refined per-block by _refine_quantity().
LABELS = {
    "att": ("attitude_sun_vector", ""),
    "sao": ("trajectory", ""),
    "tri": ("trapped_orbit_avg_spectrum", ""),
    "spp": ("trapped_proton_flux_BL", ""),
    "spe": ("trapped_electron_flux_BL", ""),
    "sef": ("solar_proton_fluence", ""),
    "seo": ("solar_proton_attenuation", ""),
    "sefflare": ("solar_proton_peak_flux", ""),
    "seoflare": ("solar_proton_peak_attenuation", ""),
    "gcf": ("gcr_flux", ""),
    "s2o": ("ionizing_dose", ""),
    "nio": ("niel", ""),
    "efo": ("eqflux", "phase2"),
    "sco": ("mc_scream", "phase2"),
}

# NIEL blocks share column names and are separated only by units.
_NIEL_BY_UNITS = {
    "MeV g^-1": "niel_dose_MeV_per_g",     # <- Phase 1 DDD deliverable
    "cm^-2": "niel_damage_equiv_fluence",
    "": "niel_relative_damage",            # DDD x damage factor -> dimensionless
}

_DEPTH = re.compile(r"^([0-9.eE+-]+)\s*(mm|micron|um|g/cm2|g cm\^-2)$")
_NUM = re.compile(r"^[0-9.eE+-]+$")

def _scalar(v):
    """Metadata values may be lists (SPECIES is 92 ion names in GCR files).
    A list belongs in the per-row `ion` column, not broadcast as a scalar."""
    return None if isinstance(v, list) else v

# ---------------------------------------------------------------- labelling

def _split_label(group: str, label: str, unit: str):
    """Decompose a melted column name into (series, depth_mm, energy_MeV, ion)."""
    series = label[len(group) + 1:] if label.startswith(group + "_") else label
    depth = energy = None
    ion = None

    m = _DEPTH.match(series)
    if m:
        v = float(m.group(1))
        u = m.group(2)
        depth = v if u == "mm" else v / 1000.0 if u in ("micron", "um") else None
        return series, depth, None, None

    if series in _ELEMENTS:
        return series, None, None, series

    if _NUM.match(series):
        v = float(series)
        # Thickness-indexed groups label by depth; flux/attenuation label by energy.
        if group.lower().startswith("dose") and "mm" in unit:
            return series, v, None, None
        return series, None, v, None

    return series, None, None, None


def _refine_quantity(block, quantity: str) -> str:
    """Disambiguate blocks that share a suffix but mean different things."""
    hdr = (block.meta.get("PLT_HDR") or "").lower()
    units = set(block.col_units)

    if quantity == "niel":
        if "spectrum" in hdr:
            return "niel_shielded_spectrum"
        for u, name in _NIEL_BY_UNITS.items():
            if u in units:
                return name
        return "niel_unknown"

    if quantity.startswith("trapped_orbit_avg"):
        mod = (block.meta.get("TRP_MOD") or "").split("-")[0].lower()
        return f"trapped_{'proton' if mod == 'ap' else 'electron'}_orbit_avg_spectrum"

    if quantity in ("solar_proton_peak_flux", "gcr_flux") and any(
        c.endswith("_He") for c in block.columns
    ):
        return quantity.replace("proton_", "") + "_ion_table"

    return quantity


# ---------------------------------------------------------------- melting

def to_long(block, prov: dict) -> pd.DataFrame:
    """Melt one Block into long format. Returns empty frame if not meltable."""
    if block.data is None or block.data.empty:
        return pd.DataFrame()

    idx_groups, val_groups = [], []
    pos = 0
    for name, unit, ncols, _desc, _ref in block.groups:
        cols = block.columns[pos:pos + ncols]
        units = block.col_units[pos:pos + ncols]
        (idx_groups if name in INDEX_NAMES and ncols == 1 else val_groups).append(
            (name, unit, cols, units)
        )
        pos += ncols

    if not idx_groups or not val_groups:
        return pd.DataFrame()  # trajectory / attenuation tables stay wide

    id_cols = [c for _n, _u, cols, _us in idx_groups for c in cols]
    suffix = prov.get("suffix", "")
    quantity = _refine_quantity(block, LABELS.get(suffix, (suffix, ""))[0])

    # Step 11 mixes annual trapped dose with a single-event solar dose:
    # the per-species columns are fine, the Total is not a physical quantity.
    mixed = "worst" in str(prov.get("step_tag", "")).lower()

    out = []
    for gname, _gunit, cols, units in val_groups:
        for col, unit in zip(cols, units):
            series, depth, energy, ion = _split_label(gname, col, unit)
            if mixed and series.lower() == "total":
                continue
            df = block.data[id_cols].copy()
            df["group"] = gname
            df["series"] = series
            df["depth_mm"] = depth
            df["energy_MeV"] = energy
            df["ion"] = ion
            df["value"] = block.data[col].values
            df["units"] = unit
            out.append(df)

    if not out:
        return pd.DataFrame()

    long = pd.concat(out, ignore_index=True)
    idx_name = idx_groups[0][0]
    long = long.rename(columns={idx_name: "index_value"})
    long["index_name"] = idx_name
    long["index_units"] = idx_groups[0][1]

    # Energy/Thick index doubles as the physical coordinate when not already set.
    if idx_name == "Energy":
        long["energy_MeV"] = long.energy_MeV.fillna(long.index_value)
    elif idx_name == "Thick":
        long["depth_mm"] = long.depth_mm.fillna(long.index_value)

    long["quantity"] = quantity
    long["model"] = _scalar(block.meta.get("TRP_MOD")) or _scalar(block.meta.get("PLT_HDR"))
    long["species"] = _scalar(block.meta.get("SPECIES"))
    long["run_time"] = block.run_time
    long["block"] = block.block_index
    long["flags"] = "total_dropped_mixed" if mixed else ""
    for k, v in prov.items():
        long[k] = v
    return long


def build_long(root, pattern="spenvis_*.txt", skip=("efo", "sco")):
    """Walk the tree and return (long_table, skipped_wide_blocks)."""
    from pathlib import Path

    root = Path(root)
    frames, wide = [], []
    for p in sorted(root.rglob(pattern)):
        try:
            prov = parse_path(p, root)
        except ValueError as e:
            warnings.warn(str(e))
            continue
        prov["suffix"] = p.stem.replace("spenvis_", "")
        if prov["suffix"] in skip:
            continue
        for b in read_spenvis(p):
            df = to_long(b, prov)
            if df.empty:
                wide.append({**prov, "block": b.block_index,
                             "columns": ";".join(b.columns)})
            else:
                frames.append(df)

    cols = ["case", "set", "alt_km", "step", "step_tag", "subpath", "quantity",
            "model", "species", "index_name", "index_value", "index_units", "group",
            "series", "depth_mm", "energy_MeV", "ion", "value", "units",
            "run_time", "rel_path", "block", "flags"]
    long = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    long = long[[c for c in cols if c in long.columns]]
    return long, pd.DataFrame(wide)