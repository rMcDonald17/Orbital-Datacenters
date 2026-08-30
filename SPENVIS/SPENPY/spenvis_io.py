

"""spenvis_io.py — format-aware, quantity-agnostic reader."""
import csv, re
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any
import pandas as pd

from datetime import datetime
_VER = re.compile(r"^'SPENVIS\s+(\S+)\s+-\s+(.+?)'\s*$")
_IDL = [(re.compile(r"!u(.*?)!n"), r"^\1"), (re.compile(r"!d(.*?)!n"), r"_\1")]
_KEY = re.compile(r"^[A-Z][A-Z0-9_]{2,7}$")

def _fields(line):
    return next(csv.reader([line], quotechar="'", skipinitialspace=True))

def _clean(s):
    s = str(s)
    for pat, rep in _IDL:
        s = pat.sub(rep, s)
    return s.strip()

def _num(tok):
    try:
        f = float(tok)
    except ValueError:
        return tok.strip()
    return int(f) if f.is_integer() and "." not in tok and "E" not in tok.upper() else f



@dataclass
class Block:
    version: str = ""
    run_time: Any = None
    meta: dict = field(default_factory=dict)
    units: dict = field(default_factory=dict)
    annotations: list = field(default_factory=list)
    columns: list = field(default_factory=list)
    col_units: list = field(default_factory=list)
    groups: list = field(default_factory=list)
    data: Any = None
    source_file: str = ""
    block_index: int = 0
    n_rows_declared: int = -1
    more_blocks: bool = False

    @property
    def model(self):
        return self.meta.get("MOD_ABB")

def _parse_meta_line(line, meta, units):
    f = _fields(line)
    key = f[0].strip()
    if not _KEY.match(key):
        return False
    n = int(float(f[1]))
    if n < 0:                                   # negative -> |n| strings
        vals = [_clean(v) for v in f[2:2 - n]]
        meta[key] = vals[0] if n == -1 else vals
    else:                                       # positive -> n numbers + unit
        vals = [_num(v) for v in f[2:2 + n]]
        meta[key] = vals[0] if n == 1 else vals
        if len(f) > 2 + n and _clean(f[2 + n]):
            units[key] = _clean(f[2 + n])
    return True

def _expand(groups, meta):
    cols, cu = [], []
    for name, unit, ncols, desc, ref in groups:
        if ncols == 1:
            labels = [name]
        else:
            lab = None
            if ref and isinstance(meta.get(ref), list) and len(meta[ref]) == ncols:
                lab = meta[ref]
            else:                               # fall back: any string list of right length
                for v in meta.values():
                    if (isinstance(v, list) and len(v) == ncols
                            and all(isinstance(x, str) for x in v)):
                        lab = v
                        break
            labels = ([f"{name}_{l}" for l in lab] if lab
                      else [f"{name}_{i+1}" for i in range(ncols)])
        cols += labels
        cu += [unit] * ncols
    return cols, cu

def _read_block(lines, i, path, bidx):
    h = _fields(lines[i])
    # was:  n_head, _nb, _nm, _na, n_desc, n_col, n_row = (int(float(x)) for x in h[1:8])
    n_head, _nb, _nm, _na, n_desc, n_col, n_row, more = (int(float(x)) for x in h[1:9])
    head = lines[i + 1: i + n_head]
    b = Block(source_file=str(path), block_index=bidx)

    desc_lines = head[-n_desc:] if n_desc else []
    for ln in head[:len(head) - n_desc]:
        if not ln.strip():
            continue
        mv = _VER.match(ln.strip())
        if mv:
            b.version = mv.group(1)
            b.run_time = datetime.strptime(mv.group(2), "%d-%b-%Y %H:%M:%S")
            continue
        if (ln.lstrip().startswith("'PS Annotation'")
                or not _parse_meta_line(ln, b.meta, b.units)):
            b.annotations.append(ln.strip())

    for ln in desc_lines:
        f = _fields(ln)
        b.groups.append((_clean(f[0]), _clean(f[1]), int(float(f[2])),
                         _clean(f[3]) if len(f) > 3 else "",
                         _clean(f[4]) if len(f) > 4 else None))
    b.columns, b.col_units = _expand(b.groups, b.meta)

    rows, j = [], i + n_head
    while j < len(lines):
        s = lines[j].strip()
        if not s or s.startswith("'"):
            break
        rows.append([float(x) for x in s.split(",")])
        j += 1
        if n_row > 0 and len(rows) == n_row:
            break
    b.data = pd.DataFrame(rows, columns=b.columns if len(b.columns) == n_col else None)   # <-- restore
    b.n_rows_declared, b.more_blocks = n_row, bool(more)
    if n_row > 0 and len(rows) != n_row:
        raise ValueError(f"{path} block {bidx}: declared {n_row} rows, read {len(rows)}")
    return b, j

def read_spenvis(path):
    """Return a list of Blocks; one file may hold several."""
    lines = Path(path).read_text(encoding="utf-8-sig", errors="replace").splitlines()
    blocks, i, bidx = [], 0, 0
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("'*'"):
            blk, i = _read_block(lines, i, path, bidx)
            blocks.append(blk)
            bidx += 1
        elif s.startswith("'End of File'"):
            break
        else:
            i += 1
    return blocks

