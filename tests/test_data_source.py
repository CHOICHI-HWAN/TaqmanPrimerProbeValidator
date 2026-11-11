"""Every genome carries where it came from (webapp/data_source.py): the record label and the profile-level summary."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
from data_source import info, record_source  # noqa: E402


def _p(kind, **extra):
    return {"source": {"type": kind, "base_url": "https://example.org/x", **extra}, "ui": {"web_url": "https://web.example.org"}}


def test_record_label_from_column_or_accession():
    assert record_source({"accession": "KY1", "data_source": "Pathoplexus · GenBank"}, _p("pathoplexus", organism="m")) == "Pathoplexus · GenBank"
    assert record_source({"accession": "NC_001608"}, _p("pathoplexus", organism="m")) == "Pathoplexus · RefSeq"
    assert record_source({"accession": "AB123456"}, _p("nextstrain")) == "Nextstrain · GenBank"
    assert record_source({"accession": "PP_000ABC"}, _p("pathoplexus", organism="m")).startswith("Pathoplexus ·")
    assert record_source({"accession": "U23417"}, _p("ncbi", taxon=186539)) == "NCBI Datasets · GenBank"


def test_profile_summary_counts_records_by_source(tmp_path):
    raw = tmp_path / "data" / "raw"
    raw.mkdir(parents=True)
    (raw / "genome_metadata.tsv").write_text("accession\tdate_updated\tdate_released\nAB1\t2025-01-02\t2024-05-01\nNC_1\t2026-03-04\t\nAB2\t\t2026-01-01\n", encoding="utf-8")
    def read(p):
        lines = p.read_text(encoding="utf-8").splitlines()
        head = lines[0].split("\t")
        return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]
    d = info(_p("pathoplexus", organism="mpox", lineage_field="clade"), tmp_path, read, ["genome"])
    assert d["name"] == "Pathoplexus" and "mpox" in d["detail"] and d["url"].endswith("/mpox/search")
    assert d["latest_update"] == "2026-03-04" and d["latest_release"] == "2026-01-01"
    assert d["records"] == 3 and d["by_record"] == {"Pathoplexus · GenBank": 2, "Pathoplexus · RefSeq": 1} and d["fetched_at"]
    assert info(_p("ncbi", taxon=186541), tmp_path, read, [])["records"] == 0
    assert info(_p("nextstrain"), tmp_path, read, [])["name"] == "Nextstrain"
