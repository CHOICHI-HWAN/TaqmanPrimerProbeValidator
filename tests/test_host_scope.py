"""Host classes / scopes and the shared per-genome pass/fail rule (scripts/host_scope.py)."""
import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import host_scope as H  # noqa: E402
from mfe_common import probe_inside_amplicon  # noqa: E402


@pytest.mark.parametrize("host,cls", [
    ("Homo sapiens", "human"), ("human", "human"), ("Human ", "human"), ("patient", "human"), ("Homo sapiens sapiens", "human"),
    ("Mastomys natalensis", "animal"), ("Pteropus lylei", "animal"), ("Rodentia", "animal"), ("Sus scrofa domesticus", "animal"),
    ("non-human primate", "animal"), ("Non-Human Primate", "animal"), ("nonhuman primate", "animal"),
    ("", "unknown"), (None, "unknown"), ("unknown", "unknown"), ("NA", "unknown"), ("Not collected", "unknown"), ("-", "unknown"),
])
def test_classify(host, cls):
    assert H.classify(host) == cls


def test_lab_flag_and_taxid():
    assert H.classify("Homo sapiens", lab="true") == "lab"          # a lab host is excluded even if the text says human
    assert H.classify("Homo sapiens", lab="") == "human"
    assert H.classify("something odd", taxid=9606) == "human"       # NCBI host tax id wins over free text
    assert H.classify("", taxid="9606") == "human"
    assert H.classify("Pteropus", taxid=9402) == "animal"


def test_scopes():
    for c in ("human", "animal", "unknown"):
        assert H.in_scope(c, "all") and H.in_scope(c, c)
    assert not H.in_scope("animal", "human") and not H.in_scope("unknown", "human")
    assert not any(H.in_scope("lab", s) for s in H.SCOPES)         # lab hosts are in no scope


def test_primary_scope_falls_back_to_all_without_humans():
    assert H.primary_scope(["a", "b"], {"a": "human", "b": "animal"}) == "human"
    assert H.primary_scope(["a", "b"], {"a": "animal", "b": "unknown"}) == "all"
    assert H.primary_scope(["a"], {}) == "all"


def test_load_host_classes(tmp_path):
    f = tmp_path / "m.tsv"
    f.write_text("accession\thost\tis_lab_host\nA\tHomo sapiens\t\nB\tMastomys natalensis\t\nC\t\t\nD\tHomo sapiens\ttrue\n", encoding="utf-8")
    assert H.load_host_classes(f, {"accession": "accession"}) == {"A": "human", "B": "animal", "C": "unknown", "D": "lab"}
    assert H.load_host_classes(tmp_path / "missing.tsv", {}) == {}
    assert H.load_host_classes(f, {"accession": "accession", "host": "no_such_column"}) == {"A": "unknown", "B": "unknown", "C": "unknown", "D": "lab"}


def _row(f, p, r, f3="1", r3="1", **kw):
    return {"f_mismatch": f, "p_mismatch": p, "r_mismatch": r, "f_3p_ok": f3, "r_3p_ok": r3, **kw}


def test_pass_fail_rule():
    ok = _row("0", "0", "1")
    assert H.genome_flags(ok, 1, 0, probe_inside_amplicon)["primer"] and H.genome_flags(ok, 1, 0, probe_inside_amplicon)["probe"]
    assert not H.primer_ok(_row("0", "0", "2"), 1) and H.primer_ok(_row("0", "0", "2"), 2)
    assert not H.primer_ok(_row("0", "0", "0", r3="0"), 2)                              # a 3'-end mismatch fails whatever the count
    assert H.primer_ok(_row("0", "0", "0", f3="", r3=""), 1)                           # blank = not checked
    assert not H.primer_ok(_row("", "0", "0"), 1)                                       # no F match
    assert not H.probe_ok(_row("0", "1", "0"), 0, probe_inside_amplicon)
    outside = _row("0", "0", "0", amplicon_start="10", amplicon_end="100", probe_start="200", probe_end="220")
    assert not H.probe_ok(outside, 0, probe_inside_amplicon)                            # a probe outside the amplicon gives no signal


def test_typed_rows_behave_like_tsv_strings():
    """The web app passes rows whose numbers/booleans are already typed."""
    s = _row("0", "0", "1")
    t = _row(0, 0, 1, f3=True, r3=True)
    assert H.genome_flags(s, 1, 0, probe_inside_amplicon) == H.genome_flags(t, 1, 0, probe_inside_amplicon)
    assert not H.primer_ok(_row(0, 0, 1, f3=False), 1)


def test_summarize_counts_and_percentages():
    rows = {"a": _row("0", "0", "0"), "b": _row("1", "0", "1"), "c": _row("0", "2", "0"), "d": _row("3", "0", "0")}
    s = H.summarize(["a", "b", "c", "d"], rows, 1, 0, probe_inside_amplicon)
    assert s["total_genomes"] == 4 and s["validated"] == 2 and s["probe_dropout"] == 1 and s["primer_matched"] == 3
    assert s["validated_pct"] == 50.0 and s["primer_coverage_mm0_pct"] == 50.0 and s["primer_coverage_mm1_pct"] == 75.0
    assert s["_validated_ids"] == {"a", "b"} and s["_dropout_ids"] == {"c"} and s["_failed_ids"] == {"d"}
    assert H.summarize([], rows, 1, 0, probe_inside_amplicon)["validated_pct"] == 0.0


def test_pipeline_summary_matches_a_recount_of_the_real_results():
    """summary.tsv is the PRIMARY host scope; counting the same genomes with this module reproduces it for every set."""
    summary = ROOT / "output" / "summary" / "summary.tsv"
    if not summary.exists():
        pytest.skip("pipeline results not computed yet")
    import re
    cfg = (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    mm = int(re.search(r"^primer_mismatch_pass:\s*(\d+)", cfg, re.M).group(1))
    pm = int(re.search(r"^probe_mismatch:\s*(\d+)", cfg, re.M).group(1))
    for row in csv.DictReader(open(summary, encoding="utf-8"), delimiter="\t"):
        if "host_scope" not in row:
            pytest.skip("summary.tsv predates the host scope column (re-run the pipeline)")
        seg = row["segment"].lower()
        ids = [ln[1:].split()[0] for ln in open(ROOT / f"data/filtered/{seg}_filtered.fasta", encoding="utf-8") if ln.startswith(">")]
        raw = ROOT / f"data/raw/{seg}_metadata.tsv"
        classes = H.load_host_classes(raw, {"accession": "accession"})
        scope = H.primary_scope(ids, classes)
        assert scope == row["host_scope"], row["set_id"]
        rows = {r["genome_id"]: r for r in csv.DictReader(open(ROOT / f"output/primers/{row['set_id']}/alignment/mismatch_table.tsv", encoding="utf-8"), delimiter="\t")}
        s = H.summarize(H.ids_in_scope(sorted(set(ids)), classes, scope), rows, mm, pm, probe_inside_amplicon)
        for k in ("total_genomes", "primer_matched", "probe_dropout", "validated"):
            assert s[k] == int(row[k]), (row["set_id"], k)
        assert s["validated_pct"] == float(row["validated_pct"]), row["set_id"]
