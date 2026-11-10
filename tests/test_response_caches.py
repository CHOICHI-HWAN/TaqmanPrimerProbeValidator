"""The API caches (TSV parse, per-set numbers) must follow the files: edit a file -> the next read sees it."""
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
sys.path.insert(0, str(ROOT / "scripts"))
import server  # noqa: E402


def test_read_tsv_follows_the_file_and_hands_out_copies(tmp_path):
    f = tmp_path / "t.tsv"
    f.write_text("a\tb\n1\t2\n")
    rows = server.read_tsv(f)
    rows[0]["a"] = "changed"                         # a caller editing its rows must not poison the cache
    assert server.read_tsv(f)[0]["a"] == "1"
    time.sleep(0.01)
    f.write_text("a\tb\n9\t2\n3\t4\n")
    os.utime(f, None)
    assert [r["a"] for r in server.read_tsv(f)] == ["9", "3"]
    assert server.read_tsv(tmp_path / "missing.tsv") == []


def test_load_params_carries_the_scan_engine():
    assert server.load_params()["scan_engine"] in ("python", "go")


def test_genome_meta_has_when_where_lineage_host():
    seg = server.seg_keys()[0]
    acc = next(iter(server.scope_ctx().ids(seg)))
    d = server.load_genome_meta(seg, acc)
    assert d and d["accession"] == acc and d["analysed"]
    assert {"date", "country", "host"} <= set(d["raw"]) and d["host_class"] in ("human", "animal", "unknown", "lab")
    assert server.load_genome_meta(seg, "NO_SUCH_ACCESSION") is None


def test_redesign_pair_ambiguity_filter():
    designs = [{"name": "a", "F": "ACGTRY", "R": "ACGTAC"},       # 2 + 0
               {"name": "b", "F": "ACGTRY", "R": "ACGNAC"},       # 2 + 1
               {"name": "c", "F": "ACGTAC", "R": "ACGTAC"}]       # 0
    res = {"status": "ok", "message": "m", "designs": designs}
    assert server.filter_pair_ambig(res, 0) is res                 # 0 = no limit
    out = server.filter_pair_ambig(res, 2)
    assert [d["name"] for d in out["designs"]] == ["a", "c"] and "1개" in out["message"]
    assert server.filter_pair_ambig(res, 5) is res                 # nothing dropped: untouched
    assert server.n_ambiguous("acgtn-R") == 3
