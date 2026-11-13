"""One organism, several reference genomes: the downloaded population is split into analysis units by the lineage column (scripts/split_raw_by_lineage.py)."""
import csv
import runpy
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _run(tmp_path, values):
    meta = tmp_path / "all.tsv"
    meta.write_text("accession\tcountry_lineage\n" + "".join(f"{a}\t{c}\n" for a, c in
                    [("A1", "Ia"), ("A2", "Ib"), ("A3", "Ib/IIb"), ("A4", "IIb"), ("A5", "unclassified"), ("A6", "IIa")]), encoding="utf-8")
    fa = tmp_path / "all.fasta"
    fa.write_text("".join(f">{a}\nACGT\n" for a in ("A1", "A2", "A3", "A4", "A5", "A6")), encoding="utf-8")
    snake = types.SimpleNamespace(input=types.SimpleNamespace(fasta=str(fa), meta=str(meta)),
                                  output=types.SimpleNamespace(fasta=str(tmp_path / "o.fasta"), meta=str(tmp_path / "o.tsv")),
                                  params=types.SimpleNamespace(values=values, col="country_lineage"))
    runpy.run_path(str(ROOT / "scripts" / "split_raw_by_lineage.py"), init_globals={"snakemake": snake})
    rows = list(csv.DictReader(open(tmp_path / "o.tsv", encoding="utf-8"), delimiter="\t"))
    heads = [ln[1:].strip() for ln in (tmp_path / "o.fasta").read_text(encoding="utf-8").splitlines() if ln.startswith(">")]
    return rows, heads


def test_unit_keeps_its_values_and_reads_a_mixed_call_as_its_first_part(tmp_path):
    rows, heads = _run(tmp_path, ["Ib"])
    assert heads == ["A2", "A3"] and [r["country_lineage"] for r in rows] == ["Ib", "Ib"]


def test_records_no_unit_lists_are_left_out(tmp_path):
    rows, heads = _run(tmp_path, ["II", "IIa", "IIb"])
    assert heads == ["A4", "A6"] and "unclassified" not in [r["country_lineage"] for r in rows]
