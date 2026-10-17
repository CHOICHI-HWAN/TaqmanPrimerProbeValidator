"""Invariants of the annotated metadata (QC verdict + clade) -- run after the data stage."""
import csv
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
ANN = ROOT / "data" / "filtered" / "genome_metadata_annotated.tsv"
pytestmark = pytest.mark.skipif(not ANN.exists(), reason="Nipah data stage not run")


def _rows():
    return list(csv.DictReader(open(ANN, encoding="utf-8"), delimiter="\t"))


def test_every_record_has_a_verdict_and_every_fail_a_reason():
    rows = _rows()
    assert {r["qc_status"] for r in rows} == {"QC PASS", "QC FAIL"}
    assert all(r["qc_reason"] for r in rows if r["qc_status"] == "QC FAIL")
    assert all(not r["qc_reason"] for r in rows if r["qc_status"] == "QC PASS")


def test_partial_sequences_are_qc_fail_and_passes_match_the_population():
    rows = _rows()
    assert all(r["qc_status"] == "QC FAIL" for r in rows if int(r["length_seq"]) < 17000)
    kept = {r["accession"] for r in csv.DictReader(open(ROOT / "data" / "filtered" / "genome_filtered_metadata.tsv"), delimiter="\t")}
    assert {r["accession"] for r in rows if r["qc_status"] == "QC PASS"} == kept


def test_clades_come_from_the_local_nextstrain_build():
    rows = [r for r in _rows() if r["qc_status"] == "QC PASS"]
    assert {r["clade"] for r in rows} <= {"I", "Ia", "Ib", "II", "unclassified"}
    assert sum(r["clade"] != "unclassified" for r in rows) / len(rows) > 0.9
