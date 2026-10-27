"""Coverage curve over the tolerated number of mismatches (scope_view.View.mismatch_curve). Synthetic tables only."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
sys.path.insert(0, str(ROOT / "scripts"))
import scope_view  # noqa: E402

COLS = "genome_id\tf_mismatch\tp_mismatch\tr_mismatch\tf_3p_ok\tr_3p_ok\tamplicon_start\tamplicon_end\tprobe_start\tprobe_end\n"


def make(tmp):
    (tmp / "data" / "raw").mkdir(parents=True)
    (tmp / "data" / "filtered").mkdir(parents=True)
    ids = [f"g{i}" for i in range(6)]
    (tmp / "data" / "raw" / "s_metadata.tsv").write_text("accession\thost\n" + "".join(f"{i}\tHomo sapiens\n" for i in ids))
    (tmp / "data" / "filtered" / "s_filtered_metadata.tsv").write_text("accession\tlineage\n" + "".join(f"{i}\tLII\n" for i in ids))
    rows = [("g0", 0, 0, 0), ("g1", 1, 0, 0), ("g2", 2, 0, 1), ("g3", 0, 1, 0), ("g4", 3, 0, 0), ("g5", 1, 2, 2)]
    t = tmp / "output" / "primers" / "X1" / "alignment"
    t.mkdir(parents=True)
    t.joinpath("mismatch_table.tsv").write_text(COLS + "".join(f"{g}\t{f}\t{p}\t{r}\t1\t1\t100\t200\t120\t140\n" for g, f, p, r in rows))


def test_curve_counts(tmp_path):
    make(tmp_path)
    v = scope_view.View(tmp_path, {"accession": "accession", "host": "host"}, 2, 0)
    c = v.mismatch_curve("X1", "s", kmax=3, scope="all")
    assert c["total"] == 6
    # probe fixed at 0: a genome needs F and R <= k (and a perfect probe): g0 k0; g1 k1; g2 k2 (r=1, f=2); g4 k3; g3 and g5 never (probe mismatch)
    assert c["fixed"] == [1, 2, 3, 4]
    # F, R and probe all tolerated up to k: g3 joins at k1 (probe 1), g5 at k2 (max 2)
    assert c["same"] == [1, 3, 5, 6]
    assert all(b >= a for a, b in zip(c["fixed"], c["fixed"][1:])) and all(b >= a for a, b in zip(c["same"], c["same"][1:]))      # a curve never falls
