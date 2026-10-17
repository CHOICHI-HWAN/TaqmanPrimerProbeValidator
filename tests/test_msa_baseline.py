"""msa.json carries what the UI needs to re-draw the dots against either baseline (oligo or consensus)."""
import json
from pathlib import Path

import pytest

from mfe_common import IUPAC

ROOT = Path(__file__).resolve().parent.parent
FILES = sorted((ROOT / "output" / "primers").glob("P*/msa/msa.json"))
pytestmark = pytest.mark.skipif(not FILES, reason="no msa.json computed yet")


@pytest.mark.parametrize("path", FILES, ids=lambda p: p.parent.parent.name)
def test_rows_carry_exact_bases_and_the_oligo_baseline(path):
    msa = json.loads(path.read_text())
    if not msa["rows"]:
        pytest.skip("empty alignment")
    n, ref = msa["columns"], msa["ref"]
    if "seq" not in msa["rows"][0]:
        pytest.skip("msa.json predates the baseline switch (re-run msa_alignment)")
    assert len(ref) == len(msa["consensus"]) == n
    for r in msa["rows"]:
        assert len(r["seq"]) == len(r["dots"]) == n
        for i, (d, b) in enumerate(zip(r["dots"], r["seq"])):
            if d == ".":                                   # a dot means: the base is covered by the baseline (IUPAC-aware)
                assert b == "-" == ref[i] or b.upper() in IUPAC.get(ref[i].upper(), ref[i].upper()), (r["genome_id"], i)
            else:
                assert d == b
    for key in ("forward", "probe", "reverse"):             # inside an oligo span the baseline IS the oligo
        rr = msa.get("ref_rows", {}).get(key)
        if rr:
            a, b = rr["col"]
            assert all(ref[i] == rr["aligned"][i] for i in range(a, b) if rr["aligned"][i] != "-")
