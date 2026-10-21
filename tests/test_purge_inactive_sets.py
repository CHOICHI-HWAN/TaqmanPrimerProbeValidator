"""Saving the primer table with "분석진행" off deletes that set's results at once (webapp/server.py purge_set_results)."""
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "webapp"))
_spec = importlib.util.spec_from_file_location("dashboard_server_purge_test", ROOT / "webapp" / "server.py")
server = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(server)


def _make(tmp_path):
    for sid in ("P001", "P002"):
        d = tmp_path / "output" / "primers" / sid / "final"
        d.mkdir(parents=True)
        (d / "summary_row.tsv").write_text("x", encoding="utf-8")
    (tmp_path / "output" / "probe_design" / "P002").mkdir(parents=True)
    (tmp_path / "output" / "summary").mkdir(parents=True)
    (tmp_path / "output" / "summary" / "summary.tsv").write_text("set_id\tv\nP001\t1\nP002\t2\n", encoding="utf-8")
    (tmp_path / "output" / "summary" / "gene_summary.tsv").write_text("set_id\tg\nP002\tNP\n", encoding="utf-8")
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "final_results.json").write_text(json.dumps({"P002": {"verdict": "hold", "note": "", "updated": "t"},
                                                                          "P001": {"verdict": "adopt", "note": "", "updated": "t"}}), encoding="utf-8")


def test_purge_removes_only_the_switched_off_set(tmp_path, monkeypatch):
    _make(tmp_path)
    monkeypatch.setattr(server, "ROOT", tmp_path)
    report = server.purge_set_results(["P002"])
    assert report == [{"id": "P002", "dirs": 2, "rows": 2}]
    assert not (tmp_path / "output" / "primers" / "P002").exists() and not (tmp_path / "output" / "probe_design" / "P002").exists()
    assert (tmp_path / "output" / "primers" / "P001" / "final" / "summary_row.tsv").exists()
    assert (tmp_path / "output" / "summary" / "summary.tsv").read_text(encoding="utf-8") == "set_id\tv\nP001\t1\n"
    assert "P002" not in (tmp_path / "config" / "final_results.json").read_text(encoding="utf-8")
    assert server.purge_set_results(["P002"]) == []                                  # nothing left: nothing reported


def test_purge_ignores_ids_that_are_not_plain_ids(tmp_path, monkeypatch):
    _make(tmp_path)
    monkeypatch.setattr(server, "ROOT", tmp_path)
    assert server.purge_set_results(["../config", "P0 01", "", None]) == []
    assert (tmp_path / "config" / "final_results.json").exists()
