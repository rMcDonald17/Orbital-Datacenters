"""Parser tests for spenvis_io.read_spenvis.

Real SPENVIS output is not redistributable, so the fixture is a hand-built
minimal file in SPENVIS format: one block, two meta keys, two column groups,
three data rows. It exercises the header layout, meta parsing (string and
numeric-with-unit), IDL markup cleaning in units, and the declared-row check.
"""
from pathlib import Path

import pytest

from spenvis_io import read_spenvis

FIX = Path(__file__).parent / "fixtures" / "spenvis_min.txt"
HEADER = "'*', 6, 0, 2, 0, 2, 2, 3, 0"


def test_reads_one_block():
    (b,) = read_spenvis(FIX)
    assert b.version == "4.6.10"
    assert b.run_time.year == 2026
    assert b.model == "AP-8"
    assert b.meta["ORB_ALT"] == 500.0
    assert b.units["ORB_ALT"] == "km"
    assert b.columns == ["Energy", "Flux"]
    assert b.col_units == ["MeV", "cm^-2s^-1"]      # IDL !u...!n markup cleaned
    assert b.data.shape == (3, 2)
    assert b.data["Flux"].iloc[0] == pytest.approx(1.0e5)


def test_row_count_mismatch_raises(tmp_path):
    bad = FIX.read_text().replace(HEADER, HEADER.replace(", 3, 0", ", 4, 0"))
    p = tmp_path / "spenvis_bad.txt"
    p.write_text(bad)
    with pytest.raises(ValueError, match="declared 4 rows, read 3"):
        read_spenvis(p)
