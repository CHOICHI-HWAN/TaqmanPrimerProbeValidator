import io

import pytest

from seqgen import synth
from xlsx_lite import read_xlsx, write_xlsx


def test_roundtrip_keeps_text_exactly_including_unicode_and_blanks():
    rows = [["primer_name", "forward_primer", "remark", "proceed"],
            ["demo-set", synth(21, 7), "한글 비고 & <tag> \"quote\"", "True"],
            ["", "", "", ""],
            ["B", synth(10, 8), "", "False"]]
    name, back = read_xlsx(write_xlsx([("Primers", rows, [20, 30, 30, 10])]))
    assert name == "Primers" and back[0] == rows[0] and back[1] == rows[1] and back[3] == rows[3]


def test_a_bad_file_is_rejected_with_a_clear_message():
    with pytest.raises(ValueError, match="xlsx"):
        read_xlsx(b"not a zip")
    with pytest.raises(ValueError):
        read_xlsx(b"PK\x03\x04junk")


def test_files_written_by_excel_like_tools_are_read_too():
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "S"
    ws.append(["name", "n", "flag"])
    ws.append(["CTL1", 12, True])                    # number + boolean cells (stored as <v>, shared strings for text)
    ws.append(["CTL2", 3.5, False])
    buf = io.BytesIO()
    wb.save(buf)
    _name, rows = read_xlsx(buf.getvalue())
    assert rows == [["name", "n", "flag"], ["CTL1", "12", "TRUE"], ["CTL2", "3.5", "FALSE"]]


def test_our_files_open_in_a_standard_library_of_record():
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.load_workbook(io.BytesIO(write_xlsx([("A", [["h1", "h2"], ["x", "y"]], [10, 10]), ("안내", [["c"], ["d"]], [10])])))
    assert wb.sheetnames == ["A", "안내"] and [c.value for c in wb["A"][2]] == ["x", "y"]


def _server():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "webapp"))
    import server
    return server


def test_upload_parser_matches_headers_loosely_and_cleans_values():
    server = _server()
    fwd = synth(10, 9)
    spaced = fwd[:5] + " " + fwd[5:]                 # typed with a space, lower case: cleaned to the bare upper-case string
    rows = [["이름", "Seg", "Forward", "Reverse", "Probe", "분석진행", "bogus"],
            ["A", "", spaced.lower(), "ACGT", "ACGTACGT", "FALSE", "x"], ["", "", "", "", "", "", ""], ["B", "s", "AC", "GT", "AC", "TRUE", ""]]
    data = server.xlsx_parse("primers", write_xlsx([("S", rows, [])]))
    assert [r["primer_name"] for r in data["rows"]] == ["A", "B"]
    assert data["rows"][0]["forward_primer"] == fwd and data["rows"][0]["proceed"] == "False" and data["rows"][1]["segment"] == "S"
    assert any("bogus" in w for w in data["warnings"]) and any("빈 행" in w for w in data["warnings"])
    with pytest.raises(ValueError, match="헤더"):
        server.xlsx_parse("standards", write_xlsx([("S", [["a", "b"], ["c", "d"]], [])]))


def test_export_then_parse_gives_the_same_table(tmp_path, monkeypatch):
    server = _server()
    root = tmp_path
    (root / "config").mkdir()
    (root / "config" / "standards.tsv").write_text("name\tsegment\tsequence\tsource\tremark\tproceed\tadd_date\nSTD1\tS\t" + "ACGT" * 10 + "\tLab\tmemo\tTrue\t2026-10-15\n")
    monkeypatch.setattr(server, "ROOT", root)
    monkeypatch.setattr(server, "seg_keys", lambda: ["s"])
    back = server.xlsx_parse("standards", server.xlsx_export("standards"))["rows"]
    assert back == [{"name": "STD1", "segment": "S", "sequence": "ACGT" * 10, "source": "Lab", "remark": "memo", "proceed": "True", "add_date": "2026-10-15"}]


def test_validate_reports_the_same_problems_as_saving_without_saving(tmp_path, monkeypatch):
    server = _server()
    monkeypatch.setattr(server, "seg_keys", lambda: ["genome"])
    rows = [{"primer_name": "A", "forward_primer": "ACGTACGTAC", "reverse_primer": "ACGTACGTAC", "probe": "ACGTACGTAC"},
            {"primer_name": "A", "forward_primer": "ACGT1CGTAC", "reverse_primer": "ACGTACGTAC", "probe": "ACGTACGTAC"}]
    errors = server.xlsx_validate("primers", rows)
    assert {(e["row"], e["field"]) for e in errors} >= {(1, "forward_primer")} and all(isinstance(e["row"], int) for e in errors)
    assert server.xlsx_validate("standards", [{"name": "S1", "sequence": "ACGT" * 10}]) == []
    assert any(e["field"] == "name" for e in server.xlsx_validate("standards", [{"name": "bad name", "sequence": "ACGT" * 10}]))
