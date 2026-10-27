"""Coverage-first redesign (scripts/redesign_coverage.py): synthetic alignments only."""
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import redesign_coverage as rc  # noqa: E402


def alignment(seed=1, n=40, length=420, rate=0.03):
    """A conserved reference with a few variable sites; every genome differs from it at `rate` of the sites."""
    rng = random.Random(seed)
    ref = "".join(rng.choice("ACGT") for _ in range(length))
    rows = []
    for _ in range(n):
        r = list(ref)
        for i in range(length):
            if rng.random() < rate:
                r[i] = rng.choice([b for b in "ACGT" if b != ref[i]])
        rows.append("".join(r))
    return [f"g{i}" for i in range(n)], rows, ref


def test_design_covers_most_genomes_and_keeps_geometry():
    ids, rows, ref = alignment()
    res = rc.design_from_alignment(ids, rows, ambig=2, probe_ambig=2, top=2, primer_tm=(40.0, 80.0), probe_above=(-30.0, 30.0), max_dtm=30.0)
    assert res["status"] == "ok" and res["designs"]
    d = res["designs"][0]
    assert d["est_cov"] >= 0.6 * len(ids)
    assert 70 <= d["amplicon_len"] <= 200
    assert rc.n_ambig(d["F"]) <= 2 and rc.n_ambig(d["R"]) <= 2 and rc.n_ambig(d["P"]) <= 2
    assert d["P"][0] != "G" or d["P_strand"] == "-"            # no 5' G on the probe as designed (checked on every resolution)


def test_estimate_equals_a_direct_recount():
    """est_cov must be exactly what the pipeline's rule gives for the designed oligos on the same alignment."""
    ids, rows, ref = alignment(seed=3, n=30)
    res = rc.design_from_alignment(ids, rows, top=1, primer_tm=(40.0, 80.0), probe_above=(-30.0, 30.0), max_dtm=30.0)
    d = res["designs"][0]
    from mfe_common import IUPAC
    cons = rc.Aligned(ids, rows).cons

    def locate(seq, rc_strand):
        import mfe_common as mc
        q = mc.reverse_complement(seq) if rc_strand else seq
        for i in range(len(cons) - len(q) + 1):
            if all(c in IUPAC.get(x, "ACGT") for c, x in zip(cons[i:i + len(q)], q)) and "-" not in cons[i:i + len(q)]:
                return i, q
        return None

    f, r, p = locate(d["F"], False), locate(d["R"], True), locate(d["P"], d["P_strand"] == "-")
    assert f and r and p

    def ok(row, start, q, kind):
        seg = row[start:start + len(q)]
        bad = [j for j, (b, c) in enumerate(zip(seg, q)) if b not in IUPAC.get(c, "ACGT")]
        if kind == "P":
            return not bad
        end = {len(q) - 1, len(q) - 2} if kind == "F" else {0, 1}
        return len(bad) <= 2 and not (set(bad) & end)

    n = sum(ok(rw, f[0], f[1], "F") and ok(rw, r[0], r[1], "R") and ok(rw, p[0], p[1], "P") for rw in rows)
    assert n == d["est_cov"]


def test_second_design_adds_genomes_or_is_dropped():
    ids, rows, ref = alignment(seed=5, n=60, rate=0.06)
    res = rc.design_from_alignment(ids, rows, top=3, primer_tm=(40.0, 80.0), probe_above=(-30.0, 30.0), max_dtm=30.0)
    assert all(d["est_added"] > 0 for d in res["designs"])


def test_no_candidates_gives_a_clear_message():
    ids, rows, ref = alignment(n=10, length=60)           # too short for a 70 bp amplicon
    res = rc.design_from_alignment(ids, rows, top=1)
    assert res["status"] == "none" and res["message"]
