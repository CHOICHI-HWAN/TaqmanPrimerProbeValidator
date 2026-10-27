"""Genome map data (webapp/genome_map.py): genes, standards placed by alignment, binding of the sets. Synthetic files only."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
import genome_map as gm  # noqa: E402


def make(tmp):
    (tmp / "config").mkdir(parents=True)
    (tmp / "config" / "refgenome_genes.tsv").write_text(
        "segment\taccession\tgenome_len\tgene\tstart\tend\tstrand\n"
        "s\tREF1\t3000\tN\t100\t1800\t+\ns\tREF1\t3000\tGPC\t1900\t2900\t-\nl\tREF2\t7000\tL\t400\t6900\t-\n")
    (tmp / "output" / "standards" / "s").mkdir(parents=True)
    (tmp / "output" / "standards" / "s" / "representativeness.json").write_text(json.dumps({
        "standards": {"STD1": {"mapped": True, "ref_start": 101, "ref_end": 900, "length": 800, "mapped_fraction": 1.0, "reversed": False, "identity_to_reference": 0.99},
                      "STD2": {"mapped": False, "length": 800}}, "per_standard": {"STD1": {"pct": 12.5}}}))
    b = tmp / "output" / "primers" / "A1" / "refgenome_binding"
    b.mkdir(parents=True)
    b.joinpath("binding.tsv").write_text("role\toligo_id\tseq\tfound\tchrom\tstart\tend\tstrand\tmismatches\n"
                                         "F\tA1_F\tAAA\t1\tREF1\t2000\t2019\t+\t0\nP\tA1_P\tCCC\t1\tREF1\t2040\t2060\t-\t1\nR\tA1_R\tGGG\t0\t\t\t\t\t\n")


def test_build(tmp_path):
    make(tmp_path)
    d = gm.build(tmp_path, "s", "REF1", [{"set_id": "A1", "name": "Set A"}, {"set_id": "A2", "name": "No result"}],
                 [{"group_id": "G1", "label": "Group", "members": ["A1"]}])
    assert d["reference"] == "REF1" and d["length"] == 3000
    assert [g["gene"] for g in d["genes"]] == ["N", "GPC"] and d["genes"][1]["strand"] == "-"
    std = {s["name"]: s for s in d["standards"]}
    assert std["STD1"]["mapped"] and (std["STD1"]["start"], std["STD1"]["end"]) == (101, 900) and std["STD1"]["matched_pct"] == 12.5
    assert not std["STD2"]["mapped"] and std["STD2"]["start"] is None             # an unaligned standard is listed, not placed
    a1 = d["sets"][0]
    assert a1["hits"]["F"]["found"] and a1["hits"]["P"]["strand"] == "-" and not a1["hits"]["R"]["found"]
    assert a1["amplicon"] == [2000, 2060]                                         # spans the oligos that were found
    assert d["sets"][1]["hits"] == {} and d["sets"][1]["amplicon"] is None         # a set without a binding file is still listed
    assert d["groups"][0]["members"] == ["A1"]


def test_other_segment_only_gets_its_genes(tmp_path):
    make(tmp_path)
    d = gm.build(tmp_path, "l", "REF2", [], [])
    assert d["length"] == 7000 and [g["gene"] for g in d["genes"]] == ["L"] and d["standards"] == []


def test_merge_hsps_prefers_the_best_strand_and_flags_weak_placements():
    rows = [  # qseqid sstart send sstrand qstart qend pident length evalue bitscore
        ["A", "2603", "3402", "minus", "1", "800", "99.75", "800", "0", "1434"],        # the real home: all 800 nt on the reverse strand
        ["A", "1", "32", "plus", "2", "34", "81.8", "33", "0.002", "30"],                # a short stray hit on the other strand
        ["B", "100", "300", "plus", "1", "200", "70.0", "200", "1e-5", "80"],            # weak and partial
    ]
    got = gm.merge_hsps(rows, {"A": 800, "B": 800})
    assert (got["A"]["start"], got["A"]["end"], got["A"]["reversed"]) == (2603, 3402, True)
    assert got["A"]["reliable"] and got["A"]["mapped_fraction"] == 1.0
    assert not got["B"]["reliable"] and got["B"]["mapped_fraction"] == 0.25            # low identity and only a quarter of it aligned


def test_a_scattered_mafft_placement_is_not_trusted_when_blastn_is_missing(tmp_path, monkeypatch):
    make(tmp_path)
    monkeypatch.setattr(gm, "find_blastn", lambda root: None)
    gm._BLAST_CACHE.clear()
    rep = tmp_path / "output" / "standards" / "s" / "representativeness.json"
    rep.write_text(json.dumps({"standards": {"X": {"mapped": True, "ref_start": 100, "ref_end": 1000, "length": 800, "mapped_fraction": 0.5, "identity_to_reference": 0.51}}}))
    d = gm.build(tmp_path, "s", "REF1", [], [])
    x = d["standards"][0]
    assert x["method"] == "mafft" and x["reliable"] is False
