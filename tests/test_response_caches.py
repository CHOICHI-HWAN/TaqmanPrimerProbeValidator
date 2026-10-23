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
