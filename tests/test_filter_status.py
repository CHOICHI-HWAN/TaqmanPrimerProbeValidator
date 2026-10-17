"""data/filtered/<seg>_filter_status.tsv must describe exactly what filter.py kept."""
import csv
from pathlib import Path

import pytest

from filter_status import classify

ROOT = Path(__file__).resolve().parent.parent


def test_classify_applies_filter_rules_in_order():
    recs = [("A1", "ACGT" * 800), ("A2", "ACGN" * 800), ("A3", "ACG"), ("A1 dup", "ACGT" * 800),
            ("A5", "ACGT" * 800), ("A6|x", "ACGT" * 800)]
    out = classify(recs, {"A1", "A2", "A3", "A5"}, 3000)
    assert [(a, s) for a, _, s in out] == [
        ("A1", "passed"), ("A2", "dropped_ambiguousN"), ("A3", "dropped_short"),
        ("A1", "dropped_dup_accession"), ("A5", "passed"), ("A6", "dropped_no_metadata")]


@pytest.mark.skipif(not (ROOT / "data" / "filtered" / "s_filter_status.tsv").exists(), reason="pipeline not run")
@pytest.mark.parametrize("seg", ["s", "l"])
def test_status_file_matches_filtered_population(seg):
    status = list(csv.DictReader(open(ROOT / "data" / "filtered" / f"{seg}_filter_status.tsv"), delimiter="\t"))
    kept = {r["accession"] for r in csv.DictReader(open(ROOT / "data" / "filtered" / f"{seg}_filtered_metadata.tsv"), delimiter="\t")}
    assert {r["accession"] for r in status if r["status"] == "passed"} == kept
    report = (ROOT / "data" / "filtered" / f"{seg}_filter_report.txt").read_text()
    for key, label in (("dropped_short", "dropped_short"), ("dropped_ambiguousN", "dropped_ambiguousN")):
        n = int(report.split(f"{label}")[1].split("=")[1].split()[0]) if label in report else None
        assert n == sum(1 for r in status if r["status"] == key)
