"""The Go probe engine must reproduce scripts/probe_designer.py (own candidates) exactly. Synthetic alignment only -- no real sequences."""
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import engine_check  # noqa: E402
engine_check.mc = __import__("mfe_common")

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


def _scan_case(seed):
    rng = random.Random(seed)
    oligo = "".join(rng.choice("ACGT") for _ in range(22))
    genomes = {}
    for i in range(30):
        g = [rng.choice("ACGT") for _ in range(900)]
        for site, strand in ((100 + i, "+"), (400 + 2 * i, "-")):
            seq = oligo if strand == "+" else engine_check.mc.reverse_complement(oligo)
            for k, c in enumerate(seq):
                g[site + k] = c if rng.random() > 0.05 else rng.choice("ACGTN")
        genomes[f"g{i}"] = "".join(g)
    return oligo, genomes


@pytest.mark.skipif(engine_check.ensure_built("scan") is None, reason="Go scan engine not available")
@pytest.mark.parametrize("seed", [1, 2])
def test_go_scan_matches_python(seed):
    mc = engine_check.mc
    oligo, genomes = _scan_case(seed)
    other = mc.reverse_complement(oligo[:20])
    for max_mm in (0, 2, 4):
        req = {"mode": "fr", "f": oligo[:20], "r": other, "three_prime_n": 1, "max_mm": max_mm, "max_size": 2000,
               "genomes": [{"id": k, "seq": v} for k, v in genomes.items()]}
        res, _ = engine_check._scan_go(req)
        go = {x["id"]: x.get("hit") for x in res}
        py = {k: mc.scan_fr(oligo[:20], other, v, 1, max_mm) for k, v in genomes.items()}
        assert any(py.values())
        assert {k: py[k] for k in py} == go
        deg = oligo[:6] + "R" + oligo[7:14] + "Y" + oligo[15:]
        res, _ = engine_check._scan_go({"mode": "probe", "probe": deg, "max_mm": max_mm + 2, "genomes": [{"id": k, "seq": v} for k, v in genomes.items()]})
        assert {x["id"]: x.get("hit") for x in res} == {k: mc.scan_probe(deg, v, None, None, max_mm + 2) for k, v in genomes.items()}
