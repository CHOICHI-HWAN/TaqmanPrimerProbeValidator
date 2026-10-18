"""Literature primer lists (webapp/literature.py): sources are mandatory, 'verified' needs a note, nothing is trusted by default."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "webapp"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import literature as L  # noqa: E402
from seqgen import synth  # noqa: E402


def row(**kw):
    r = {"name": "set one", "segment": "S", "forward_primer": synth(20, 1), "reverse_primer": synth(21, 2), "probe": synth(24, 3),
         "citation": "Author A, et al. Journal. 2020", "doi": "", "pmid": "", "verified": "no", "verified_note": ""}
    r.update(kw)
    return r


def test_new_entries_are_unverified_and_get_ids():
    rows = L.normalise([row(verified="maybe"), row(name="two")], "lasv")
    assert [r["lit_id"] for r in rows] == ["LIT-LASV-001", "LIT-LASV-002"]
    assert all(r["verified"] == "no" for r in rows)


def test_sequences_are_normalised():
    r = L.normalise([row(forward_primer=" acg-tu \n")], "lasv")[0]
    assert r["forward_primer"] == "ACGTT"


def test_a_source_is_required():
    errs = L.validate(L.normalise([row(citation="")], "lasv"), "lasv")
    assert any(e["field"] == "citation" for e in errs)
    assert not L.validate(L.normalise([row(citation="", doi="10.1128/JCM.00000-00")], "lasv"), "lasv")


def test_verified_needs_a_note():
    errs = L.validate(L.normalise([row(verified="yes")], "lasv"), "lasv")
    assert [e["field"] for e in errs] == ["verified_note"]
    assert not L.validate(L.normalise([row(verified="yes", verified_note="Table 2")], "lasv"), "lasv")


def test_segment_depends_on_the_virus():
    assert L.validate(L.normalise([row(segment="GENOME")], "lasv"), "lasv")
    assert not L.validate(L.normalise([row(segment="GENOME")], "nipah"), "nipah")
    assert L.validate(L.normalise([row(segment="S")], "nipah"), "nipah")


def test_bad_characters_duplicates_and_formats():
    rows = L.normalise([row(forward_primer="ACGTXZ"), row(), row(doi="nope", pmid="12a")], "lasv")
    fields = {(e["row"], e["field"]) for e in L.validate(rows, "lasv")}
    assert (0, "forward_primer") in fields and (2, "name") in fields and (2, "doi") in fields and (2, "pmid") in fields


def test_save_and_load_roundtrip(tmp_path):
    rows = L.normalise([row(verified="yes", verified_note="Table 2", doi="10.1/abc")], "lasv")
    L.save(tmp_path, "lasv", rows)
    assert L.load(tmp_path, "lasv") == rows
    assert L.load(tmp_path, "nipah") == []                       # no file = empty list, never an error


def test_manifest_row_states_the_verification_status():
    r = L.normalise([row()], "lasv")[0]
    m = L.manifest_row(r)
    assert "대조 안 됨" in m["remark"] and m["source"].startswith("문헌:") and m["proceed"] == "True"
    r["verified"] = "yes"
    assert "원문 대조 완료" in L.manifest_row(r)["remark"]


def test_same_sequences_ignores_name_and_case():
    r = L.normalise([row()], "lasv")[0]
    assert L.same_sequences([{"forward_primer": r["forward_primer"].lower(), "reverse_primer": r["reverse_primer"], "probe": r["probe"]}], r)
    assert not L.same_sequences([{"forward_primer": r["forward_primer"], "reverse_primer": r["reverse_primer"], "probe": synth(24, 9)}], r)


def test_incomplete_rows_cannot_be_analysed():
    assert not L.complete(L.normalise([row(probe="")], "lasv")[0])
