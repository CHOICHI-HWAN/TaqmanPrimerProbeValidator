"""The Go probe engine must reproduce scripts/probe_designer.py (own candidates) exactly. Synthetic alignment only -- no real sequences."""
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import engine_check  # noqa: E402

pytestmark = pytest.mark.skipif(engine_check.ensure_built() is None, reason="Go engine not available (install Go to run this test)")


def synthetic_msa(seed, n_rows=80, n_cols=150):
    rng = random.Random(seed)
    cons = "".join(rng.choice("ACGT") for _ in range(n_cols))
    rows = []
    for _ in range(n_rows):
        dots = ["."] * n_cols
        for c in range(n_cols):
            if rng.random() < 0.04:
                dots[c] = rng.choice([b for b in "ACGT" if b != cons[c]] + ["-"])
        rows.append({"genome_id": f"g{len(rows)}", "dots": "".join(dots)})
    return {"consensus": cons, "f_col": [0, 19], "r_col": [n_cols - 20, n_cols - 1], "rows": rows}


@pytest.mark.parametrize("seed", [1, 2, 3])
@pytest.mark.parametrize("max_iupac", [0, 2, 4])
def test_go_engine_matches_python(seed, max_iupac):
    res = engine_check.compare(synthetic_msa(seed), top=15, max_iupac=max_iupac, tm_above=(-30.0, 30.0))
    assert res["n_candidates"] > 100
    bad = [c for c in res["checks"] if not c["ok"]]
    assert not bad, bad[:3]
