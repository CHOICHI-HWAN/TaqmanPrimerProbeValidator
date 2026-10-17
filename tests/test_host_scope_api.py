"""The web API restricts every number and table to the requested host scope (?host=human|animal|unknown|all)."""
import csv
import json
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCOPES = ("human", "animal", "unknown", "all")

pytestmark = pytest.mark.skipif(
    not (ROOT / "output" / "summary" / "summary.tsv").exists(), reason="pipeline results not computed yet (run snakemake first)")


@pytest.fixture(scope="module")
def server():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
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
    yield base
    proc.terminate()
    proc.wait(timeout=5)


def get(base, path):
    with urllib.request.urlopen(base + path, timeout=60) as r:
        return json.loads(r.read())


def summary_rows():
    with open(ROOT / "output" / "summary" / "summary.tsv", newline="", encoding="utf-8") as fh:
        return {r["set_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_default_scope_equals_the_pipeline_summary(server):
    pipe = summary_rows()
    for s in get(server, "/api/sets"):
        if s["set_id"] in pipe:
            p = pipe[s["set_id"]]
            assert s["summary"]["total_genomes"] == int(p["total_genomes"]) and s["summary"]["validated"] == int(p["validated"])
            assert s["summary"]["validated_pct"] == float(p["validated_pct"])


def test_scopes_add_up_and_lab_hosts_are_in_none(server):
    by = {sc: {s["set_id"]: s["summary"] for s in get(server, f"/api/sets?host={sc}")} for sc in SCOPES}
    for sid in summary_rows():
        for k in ("total_genomes", "validated", "primer_matched", "probe_dropout"):
            assert by["human"][sid][k] + by["animal"][sid][k] + by["unknown"][sid][k] == by["all"][sid][k], (sid, k)
    info = get(server, "/api/host-scope")["segments"]
    for seg, i in info.items():
        assert i["scope_counts"]["all"] == i["counts"]["human"] + i["counts"]["animal"] + i["counts"]["unknown"]    # lab hosts are not counted


def test_unknown_scope_value_means_primary(server):
    a = get(server, "/api/sets?host=bogus")
    b = get(server, "/api/sets")
    assert [s["summary"] for s in a] == [s["summary"] for s in b]


def test_set_detail_lists_only_the_scope_and_tags_host_class(server):
    sid = sorted(summary_rows())[0]
    for sc in SCOPES:
        d = get(server, f"/api/sets/{sid}?host={sc}")
        classes = {g["host_class"] for g in d["genomes"]}
        assert "lab" not in classes
        if sc != "all":
            assert classes <= {sc}
        assert len(d["genomes"]) == d["summary"]["total_genomes"]
        assert len(d["population_ids"]) == d["summary"]["total_genomes"]
        if d.get("msa"):
            assert all(r["host_class"] in ("human", "animal", "unknown") for r in d["msa"]["rows"])
            assert len(d["msa"]["rows"]) + len(d["msa"].get("excluded") or []) <= d["summary"]["total_genomes"]   # only genomes with a match are aligned
        assert sum(d["lineage_totals"].values()) == d["summary"]["total_genomes"]


def test_overview_and_groups_follow_the_scope(server):
    for sc in SCOPES:
        o = get(server, f"/api/overview?host={sc}")
        info = o["host_scope"]
        for seg, n in o["dataset_sizes"].items():
            assert n == info[seg]["scope_counts"][sc]
        for g in o["groups"]:
            assert g["total_genomes"] == info[g["segment"]]["scope_counts"][sc]
            assert g["validated"] <= g["validated_with_cross"] <= g["total_genomes"]
        # the per-lineage rows of a set add up to its validated count
        by_set = {}
        for r in o["lineage"]:
            by_set[r["set_id"]] = by_set.get(r["set_id"], 0) + int(r["validated_count"])
        for s in o["sets"]:
            if s["set_id"] in by_set and s.get("builtin", True) and s["set_id"].startswith("P"):
                assert by_set[s["set_id"]] == s["summary"]["validated"], (sc, s["set_id"])


def test_group_detail_counts_inside_the_scope(server):
    groups = get(server, "/api/groups")
    if not groups:
        pytest.skip("no groups")
    gid = groups[0]["group_id"]
    for sc in SCOPES:
        d = get(server, f"/api/groups/{gid}?host={sc}")
        keep = set(d["population_ids"])
        assert d["summary"]["total_genomes"] == len(keep)
        assert set(d["validated_ids"]) <= keep and set(d["validated_ids_with_cross"]) <= keep
        assert all(r["genome_id"] in keep for r in d["combo_genomes"])
        assert sum(d["lineage_totals"].values()) == len(keep)


def test_metadata_table_has_the_host_class(server):
    t = get(server, "/api/metadata/table")
    assert "host_class" in [c["key"] for c in t["columns"]]
    classes = {r["host_class"] for r in t["rows"]}
    assert "사람" in classes and classes <= {"사람", "동물", "미상", "실험 숙주 (제외)", ""}
