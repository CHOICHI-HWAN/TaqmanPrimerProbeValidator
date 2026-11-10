"""The reference-standard manifest (config/standards.tsv) and its sync to the legacy FASTA files."""
from standards_common import (
    ensure_standards_manifest, read_manifest, sync_standards, validate_rows, write_manifest,
)

SEQ = "ACGTRYACGT" * 6


def _row(**kw):
    base = {"name": "NP1", "segment": "S", "sequence": SEQ, "source": "KDCA 2024", "remark": "", "proceed": "True", "add_date": ""}
    base.update(kw)
    return base


def test_validation_rules():
    assert validate_rows([_row(), {c: "" for c in ("name", "segment", "sequence")}]) == []
    errs = validate_rows([_row(), _row(),                                # duplicate name
                          _row(name="bad name"), _row(name="X", segment="Q"),
                          _row(name="Y", sequence="ACGT!!" + SEQ), _row(name="Z", sequence="ACGT"),
                          _row(name="W", sequence="")])
    fields = {(e["row"], e["field"]) for e in errs}
    assert {(1, "name"), (2, "name"), (3, "segment"), (4, "sequence"), (5, "sequence"), (6, "sequence")} <= fields


def test_active_standard_required_for_segments_with_primers():
    rows = [_row(), _row(name="L1", segment="L", proceed="False")]
    errs = validate_rows(rows, required_segments=("S", "L"))
    assert [e for e in errs if e["row"] == -1 and "L" in e["message"]]
    assert not [e for e in errs if e["row"] == -1 and "segment S" in e["message"]]


def test_write_stamps_date_and_sync_generates_pipeline_files(tmp_path):
    m = tmp_path / "standards.tsv"
    out = write_manifest(m, [_row(sequence=" acgt\nacgt " + SEQ), _row(name="L1", segment="L", proceed="False"),
                             _row(name="GC1", source="")])
    assert out[0]["add_date"] and out[0]["sequence"].startswith("ACGTACGT")
    changed = sync_standards(m, tmp_path)
    assert sorted(changed) == ["kdca_standard_s.fasta", "kdca_standard_s_metadata.tsv"]   # inactive L row writes nothing
    assert (tmp_path / "kdca_standard_s.fasta").read_text().splitlines()[0] == ">NP1"
    meta = (tmp_path / "kdca_standard_s_metadata.tsv").read_bytes().decode().split("\r\n")
    assert meta[0] == "accession\tlineage\tdate\tcountry\tregion" and meta[1].startswith("NP1\tNP1\tKDCA_2024\t")
    assert meta[2].startswith("GC1\tGC1\tSTANDARD\t")
    assert sync_standards(m, tmp_path) == []          # unchanged -> nothing rewritten (no Snakemake re-run)


def test_legacy_files_are_migrated_once(tmp_path):
    (tmp_path / "kdca_standard_s.fasta").write_text(f">A1\n{SEQ}\n>A2\n{SEQ}\n")
    (tmp_path / "kdca_standard_s_metadata.tsv").write_bytes(b"accession\tlineage\tdate\tcountry\tregion\r\nA1\tA1\tOLD_SRC\t-\t-\r\n")
    m = tmp_path / "standards.tsv"
    rows = ensure_standards_manifest(m, tmp_path)
    assert [(r["name"], r["segment"], r["source"]) for r in rows] == [("A1", "S", "OLD_SRC"), ("A2", "S", "KDCA_2024_STANDARD")]
    before = (tmp_path / "kdca_standard_s.fasta").read_bytes()
    sync_standards(m, tmp_path)
    assert (tmp_path / "kdca_standard_s.fasta").read_bytes() == before
    (tmp_path / "kdca_standard_s.fasta").write_text("changed")
    assert len(ensure_standards_manifest(m, tmp_path)) == 2 and read_manifest(m)[0]["name"] == "A1"   # existing manifest wins


def test_a_segment_whose_standards_were_all_removed_loses_its_files(tmp_path):
    m = tmp_path / "standards.tsv"
    write_manifest(m, [_row(), _row(name="L1", segment="L")])
    sync_standards(m, tmp_path)
    (tmp_path / "kdca_standard_l.fasta.fai").write_text("index")          # mfeprimer's index beside it
    write_manifest(m, [_row(), _row(name="L1", segment="L", proceed="False")])
    changed = sync_standards(m, tmp_path)
    assert {"kdca_standard_l.fasta", "kdca_standard_l.fasta.fai", "kdca_standard_l_metadata.tsv"} <= set(changed)
    assert not (tmp_path / "kdca_standard_l.fasta").exists() and (tmp_path / "kdca_standard_s.fasta").exists()
    write_manifest(m, [])                                                  # no standards at all
    sync_standards(m, tmp_path)
    assert not list(tmp_path.glob("kdca_standard_*"))
