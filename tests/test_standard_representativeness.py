"""Standard-material representativeness: synthetic reference + genomes with known answers."""
import random

from standard_representativeness import compute, grade_for, main_block

PARAMS = {"identity": 0.9, "min_coverage": 0.8, "good": 90, "critical": 50, "site_tol_primer": 1,
          "site_tol_probe": 0, "identity_sensitivity": [0.95, 0.9], "site_sensitivity": [(0, 0), (1, 0)]}


def _world(n_same=6, n_diff=4):
    rnd = random.Random(1)
    ref = "".join(rnd.choice("ACGT") for _ in range(200))
    std = ref[20:100]                       # standard = ref[21..100]
    samples, bases = [], {}
    genomes = []
    for i in range(n_same + n_diff):
        g = list(ref)
        if i >= n_same:                     # 'divergent' genomes: every 4th base inside the span changed
            for p in range(20, 100, 4):
                g[p] = "A" if g[p] != "A" else "C"
        genomes.append("".join(g))
        samples.append(f"G{i}")
    for pos in range(1, 201):
        bases[pos] = [g[pos - 1] for g in genomes]
    aligned = {"STD": "-" * 20 + std + "-" * 100}
    lineage = {s: ("LA" if i < n_same else "LB") for i, s in enumerate(samples)}
    windows = [{"start": p, "end": p, "pi": 0.3 if 20 < p <= 100 else 0.05, "no_sites": 1} for p in range(1, 201)]
    return ref, std, aligned, samples, bases, lineage, windows


def _binding(start_end):
    return [{"role": r, "found": "1", "start": str(a), "end": str(b)} for r, (a, b) in start_end.items()]


def test_grade_thresholds():
    assert [grade_for(p, 90, 50) for p in (95, 90, 89.9, 50, 49.9, None)] == ["good", "good", "review", "review", "critical", "na"]


def test_main_block_drops_stray_fragments():
    cols = {p: "A" for p in range(100, 140)}
    cols.update({5: "C", 400: "G"})
    assert set(main_block(cols, 40)) == set(range(100, 140))


def test_layer1_union_and_lineages():
    ref, std, aligned, samples, bases, lineage, windows = _world()
    r = compute("s", {"STD": std}, "REF", ref, aligned, samples, bases, lineage, windows, {}, PARAMS)
    u = r["union"]
    assert (u["assessable"], u["matched"], u["pct"], u["grade"]) == (10, 6, 60.0, "review")
    assert u["by_lineage"] == {"LA": {"n": 6, "matched": 6}, "LB": {"n": 4, "matched": 0}}
    assert {x["genome_id"] for x in u["unmatched_detail"]} == {"G6", "G7", "G8", "G9"}
    assert r["standards"]["STD"]["ref_start"] == 21 and r["standards"]["STD"]["location"] in ("variable", "moderate")


def test_layer2_sites_and_not_applicable():
    ref, std, aligned, samples, bases, lineage, windows = _world()
    inside = _binding({"F": (24, 40), "P": (45, 60), "R": (80, 96)})
    outside = _binding({"F": (130, 145), "P": (150, 160), "R": (170, 185)})
    missing = [{"role": "F", "found": "0", "start": "", "end": ""}]
    r = compute("s", {"STD": std}, "REF", ref, aligned, samples, bases, lineage, windows,
                {"P1": inside, "P2": outside, "P3": missing}, PARAMS)
    p1 = r["sets"]["P1"]
    assert p1["applicable"] and p1["covering_standards"] == ["STD"]
    assert p1["sites"]["F"]["location"] in ("variable", "moderate")
    assert p1["pct"] is not None and p1["grade"] in ("good", "review", "critical")
    assert r["sets"]["P2"] == {"applicable": False, "reason": "sites_outside_standard", "covering_standards": ["STD"][:0],
                               "sites": r["sets"]["P2"]["sites"], "grade": "na"}
    assert r["sets"]["P3"]["reason"] == "no_binding_on_reference"
