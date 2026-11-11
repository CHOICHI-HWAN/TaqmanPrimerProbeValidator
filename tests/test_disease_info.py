"""The shipped virus profile carries the text of the 감염병 정보 / 감염 경과 pages (profiles/<id>.json "disease" and "course")."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _profile():
    cfg = (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    name = re.search(r'^profile:\s*"([^"]+)"', cfg, re.M).group(1)
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def test_disease_block_is_complete_and_links_are_https():
    d = _profile().get("disease")
    if not d:                                   # RVFV: the info pages are not written yet
        return
    assert d["title"] and d["summary"] and d["note"]
    assert len(d["facts"]) >= 4 and all(f["k"] and f["v"] for f in d["facts"])
    assert d["sections"] and all(s["h"] and s["p"] for s in d["sections"])
    assert d["refs"] and all(r["url"].startswith("https://") and r["label"] for r in d["refs"])


def test_course_numbers_are_consistent_when_present():
    c = _profile().get("course")
    if not c:                                   # viruses without a described human course (e.g. Reston, Tai Forest)
        return
    lo, hi = c["incubation"]
    assert 0 < lo <= c["onset"] <= hi
    assert c["detect"][0] <= c["peak"][0] <= c["peak"][1] <= c["detect"][1]
    assert 0 <= c["igm"] <= c["igg"]
    assert c["specimens"] and all(s["name"] and s["note"] for s in c["specimens"]) and c["note"]
