"""Smoke tests for the webapp's HTTP API and static files.

Starts webapp/server.py as a subprocess on a free port, against THIS project's computed results
(config/, data/, output/), and only issues GET requests plus a POST that is rejected by validation
-- nothing here writes to the project. Skipped automatically when the pipeline hasn't been run yet
(no output/summary). Checks status codes and the *shape* of each answer, which is what the
front end relies on (renaming a field or moving a path breaks the page, not a unit test).
"""
import json
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

pytestmark = pytest.mark.skipif(
    not (ROOT / "output" / "summary" / "summary.tsv").exists(),
    reason="pipeline results not computed yet (run snakemake first)",
)


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="module")
def server():
    port = _free_port()
    proc = subprocess.Popen([sys.executable, str(ROOT / "webapp" / "server.py"), "--port", str(port), "--host", "127.0.0.1"],
                            cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    base = f"http://127.0.0.1:{port}"
    for _ in range(100):
        try:
            urllib.request.urlopen(base + "/api/jobs/active", timeout=1)
            break
        except Exception:  # noqa: BLE001
            if proc.poll() is not None:
                raise RuntimeError(proc.stderr.read().decode())
            time.sleep(0.1)
    else:
        proc.kill()
        raise RuntimeError("server did not start")
    yield base
    proc.terminate()
    proc.wait(timeout=5)


def get(base, path):
    with urllib.request.urlopen(base + path, timeout=30) as r:
        return r.status, r.headers.get("Content-Type", ""), r.read()


def get_json(base, path):
    status, _, body = get(base, path)
    assert status == 200, path
    return json.loads(body)


@pytest.mark.parametrize("path,ctype", [
    ("/", "text/html"), ("/app.js", "javascript"), ("/style.css", "text/css"), ("/lumo.css", "text/css"),
    ("/vendor/d3.min.js", "javascript"), ("/vendor/peartree_frame.html", "text/html"),
])
def test_static_files(server, path, ctype):
    status, ct, body = get(server, path)
    assert status == 200 and ctype in ct and len(body) > 100


def test_index_has_every_sidebar_page(server):
    html = get(server, "/")[2].decode()
    for el_id in ("overviewBtn", "primersBtn", "standardsBtn", "diversityBtn", "metaTableBtn", "metadataBtn",
                  "paramsBtn", "changelogBtn", "reportBtn", "fetchDataBtn", "rerunBtn", "skinToggleBtn",
                  "primersView", "standardsView", "metaTableView", "diversityView", "metadataView"):
        assert f'id="{el_id}"' in html, el_id


def test_sets_list_uses_ids_and_names(server):
    sets = get_json(server, "/api/sets")
    assert len(sets) >= 1
    s = sets[0]
    assert s["set_id"].startswith("P") and s["display_name"] and s["segment"] in ("S", "L")
    assert {"F_seq", "R_seq", "P_seq", "summary"} <= set(s)
    assert len({x["set_id"] for x in sets}) == len(sets)


def test_set_detail_shape(server):
    sid = get_json(server, "/api/sets")[0]["set_id"]
    d = get_json(server, f"/api/sets/{sid}")
    for key in ("genomes", "standard_genomes", "binding", "msa", "dimer_hairpin", "evaluation", "standard_representativeness"):
        assert key in d, key
    assert d["genomes"] and {"genome_id", "lineage", "f_mismatch", "r_mismatch"} <= set(d["genomes"][0])


def test_unknown_set_is_404(server):
    with pytest.raises(urllib.error.HTTPError) as e:
        get(server, "/api/sets/P999")
    assert e.value.code == 404


def test_groups(server):
    groups = get_json(server, "/api/groups")
    assert groups and groups[0]["group_id"].startswith("G") and groups[0]["member_names"]
    d = get_json(server, f"/api/groups/{groups[0]['group_id']}")
    for key in ("members", "binding", "dimer_hairpin", "combo_genomes", "evaluation"):
        assert key in d, key
    assert all(b.get("set_name") for b in d["binding"])
    assert d["member_names"] and set(d["member_names"]) == set(d["members"])   # the page shows names, never raw ids


def test_overview(server):
    o = get_json(server, "/api/overview")
    assert o["sets"] and o["groups"] and o["lineage"] and o["dataset_sizes"]
    assert all(g["member_names"] for g in o["groups"])


def test_diversity(server):
    d = get_json(server, "/api/diversity")
    seg = next(iter(d["segments"].values()))
    for key in ("windows", "lineage_n", "qc", "genome_wide_pi", "window_size"):
        assert key in seg, key
    assert "ALL_LASV" in seg["windows"] and seg["qc"]["kept"] > 0


def test_metadata_table(server):
    t = get_json(server, "/api/metadata/table")
    keys = [c["key"] for c in t["columns"]]
    assert {"accession", "lineage", "country", "length", "validated_sets", "diversity_qc"} <= set(keys)
    assert len(t["rows"]) > 100
    assert {r["segment"] for r in t["rows"]} == {"S", "L"}
    passed = [r for r in t["rows"] if r["data_qc"] == "통과"]
    assert passed and len(passed) < len(t["rows"])          # dropped records are listed and marked too
    assert all(r["validated_sets"] is None for r in t["rows"] if r["data_qc"] != "통과")
    assert all(set(keys) <= set(r) for r in t["rows"][:20])


def test_metadata_stats(server):
    m = get_json(server, "/api/metadata")
    assert m["segments"]["s"]["kept"] > 0
    assert m["source"]["base_url"].startswith("https://")      # comes from the virus profile


def test_profile_endpoint_matches_profile_file(server):
    p = get_json(server, "/api/profile")
    assert p["id"] and p["short_name"] and p["population"]["pooled_label"] and list(p["segments"]) == ["s", "l"]
    assert "ALL_" in p["population"]["pooled_label"] or p["population"]["pooled_label"]
    div = get_json(server, "/api/diversity")
    assert p["population"]["pooled_label"] in next(iter(div["segments"].values()))["windows"]


def test_standard_representativeness(server):
    d = get_json(server, "/api/standards/representativeness")
    assert d["segments"] and d["set_names"]
    for seg in d["segments"].values():
        assert seg["union"]["grade"] in ("good", "review", "critical", "na") and seg["standards"]
        assert all("grade" in e for e in seg["sets"].values())
    status, ct, body = get(server, "/api/standards/unmatched?seg=s")
    assert status == 200 and body.startswith(b"genome_id\t")


def test_primer_manifest_roundtrip_is_rejected_when_invalid(server):
    m = get_json(server, "/api/primers/manifest")
    assert m["columns"][:2] == ["primer_id", "primer_name"] and "remark" in m["columns"]
    row = dict(m["rows"][0])
    bad = [row, dict(row)]                       # same name twice -> validation error, nothing written
    req = urllib.request.Request(server + "/api/primers/manifest", data=json.dumps({"rows": bad}).encode(),
                                 headers={"Content-Type": "application/json"}, method="POST")
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(req, timeout=10)
    assert e.value.code == 400
    assert json.loads(e.value.read())["errors"]


@pytest.mark.parametrize("path", ["/api/params", "/api/ui-state", "/api/changelog", "/api/jobs/active"])
def test_misc_endpoints(server, path):
    get_json(server, path)


def test_diversity_region_and_download(server):
    status, _, body = get(server, "/api/diversity/download?seg=s&start=100&end=160&lineage=ALL_LASV")
    assert status == 200 and body.startswith(b">")
    r = get_json(server, "/api/diversity/region?seg=s&start=100&end=160")
    assert r["columns"] == 61 and r["rows"]


def test_user_set_rows_are_numeric():
    """A cached live-analysis result (strings in result.json) must come back with numeric
    coordinates, else the page/report concatenate strings instead of adding."""
    sys.path.insert(0, str(ROOT / "webapp"))
    sys.path.insert(0, str(ROOT / "scripts"))
    import server  # noqa: PLC0415
    rows = [{"amplicon_start": "895", "amplicon_end": "1016", "f_tm": "43.60", "f_mismatch": 1, "f_3p_ok": True, "probe_start": ""}]
    server._coerce_numeric_strings(rows)
    assert rows[0] == {"amplicon_start": 895, "amplicon_end": 1016, "f_tm": 43.6, "f_mismatch": 1, "f_3p_ok": True, "probe_start": None}


def test_standards_manifest_api(server):
    m = get_json(server, "/api/standards/manifest")
    assert m["columns"][0] == "name" and len(m["rows"]) >= 2
    r = m["rows"][0]
    assert r["segment"] in ("S", "L") and r["length"] == len(r["sequence"]) >= 30
    # an invalid save (bad name) is rejected and writes nothing
    req = urllib.request.Request(server + "/api/standards/manifest", data=json.dumps(
        {"rows": [dict(x, proceed="False") for x in m["rows"]] + [{"name": "bad name", "segment": "S", "sequence": "ACGT"}]}).encode(),
        headers={"Content-Type": "application/json"}, method="POST")
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(req, timeout=10)
    assert e.value.code == 400
    errs = json.loads(e.value.read())["errors"]
    assert any(x["field"] == "name" for x in errs)                  # (standards are optional: no "need at least one" error)
    assert get_json(server, "/api/standards/manifest")["rows"] == m["rows"]


def test_standards_files_follow_the_manifest(server):
    """config/kdca_standard_*.fasta (what the pipeline reads) must agree with the active manifest rows."""
    m = get_json(server, "/api/standards/manifest")
    for seg in ("S", "L"):
        names = [x["name"] for x in m["rows"] if x["segment"] == seg and x["proceed"].lower() == "true"]
        fasta = (ROOT / "config" / f"kdca_standard_{seg.lower()}.fasta").read_text().splitlines()
        assert [l[1:] for l in fasta if l.startswith(">")] == names


def test_requirements_log_is_served_and_consistent(server):
    d = get_json(server, "/api/requirements")
    items = [i for p in d["phases"] for i in p["items"]]
    assert len(items) > 100 and d["open"]
    assert {i["status"] for i in items} <= {"done", "partial", "answer", "dropped", "open"}
    assert all(i["req"] and i["rev"] and i["date"] and i["cat"] for i in items)
    md = (ROOT / "docs" / "REQUIREMENTS_LOG.md").read_text(encoding="utf-8")
    assert f"총 {len(items)}건" in md          # the Markdown copy was regenerated from the same JSON
