"""scripts/standard_placement.py: blastn rows in reference coordinates, only for reliably placed standards. Synthetic sequences; skipped without blastn."""
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import standard_placement as sp  # noqa: E402
from mfe_common import reverse_complement  # noqa: E402

pytestmark = pytest.mark.skipif(sp.find_blastn() is None, reason="blastn not available")


def make(tmp_path, n=1500, seed=4):
    rng = random.Random(seed)
    ref = "".join(rng.choice("ACGT") for _ in range(n))
    (tmp_path / "ref.fasta").write_text(">REF\n" + ref + "\n")
    return ref


def test_forward_and_reverse_standards_land_on_the_right_columns(tmp_path):
    ref = make(tmp_path)
    fwd, rev = ref[200:500], reverse_complement(ref[900:1200])
    rows, reversed_ids = sp.blast_rows(sp.find_blastn(), tmp_path / "ref.fasta", {"F1": fwd, "R1": rev})
    assert rows["F1"][200] != "-" and rows["F1"][199] == "-" and rows["F1"][499] != "-" and rows["F1"][500] == "-"
    assert rows["F1"][200:500] == ref[200:500] and "F1" not in reversed_ids
    assert reversed_ids == ["R1"] and rows["R1"][900:1200] == ref[900:1200]            # written in reference orientation, flagged as reversed


def test_a_standard_that_barely_matches_is_not_placed(tmp_path):
    ref = make(tmp_path)
    junk = "".join(random.Random(9).choice("ACGT") for _ in range(300))
    part = ref[100:160] + junk[:240]                                                 # only 20% of it belongs to the reference
    rows, _ = sp.blast_rows(sp.find_blastn(), tmp_path / "ref.fasta", {"J": junk, "P": part})
    assert "J" not in rows and "P" not in rows                                       # the caller keeps the mafft row for these
