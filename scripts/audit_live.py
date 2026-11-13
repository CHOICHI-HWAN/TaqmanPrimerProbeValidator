#!/usr/bin/env python3
"""Does the web app's live analysis ("새 프라이머 세트 분석") give exactly what the Snakemake pipeline gave for the same set?

Re-analyses chosen primer sets with scripts/live_analyze.py (real mfeprimer, in a temp folder) and compares F / R / probe
mismatches and the validated count with output/primers/<ID>/alignment/mismatch_table.tsv and output/summary/summary.tsv.
    python3 scripts/audit_live.py            # the first two active sets
    python3 scripts/audit_live.py P001 P007
Exit status 1 on any difference.
"""
import csv
import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import live_analyze as L  # noqa: E402


def main():
    wanted = sys.argv[1:]
    rows = [r for r in csv.DictReader(open(ROOT / "config" / "primers.tsv", encoding="utf-8"), delimiter="\t")
            if r.get("proceed", "True").lower() == "true"]
    chosen = [r for r in rows if r["primer_id"] in wanted] if wanted else rows[:2]
    summary = {r["set_id"]: r for r in csv.DictReader(open(ROOT / "output" / "summary" / "summary.tsv", encoding="utf-8"), delimiter="\t")}
    bad = 0
    for r in chosen:
        pid = r["primer_id"]
        table = ROOT / "output" / "primers" / pid / "alignment" / "mismatch_table.tsv"
        if not table.exists():
            print(f"{pid}: no pipeline result yet -- skipped")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            res = L.analyze_set(name=f"audit-{pid}", segment=r["segment"], f_id="F", f_seq=r["forward_primer"], r_id="R", r_seq=r["reverse_primer"],
                                p_id="P", p_seq=r["probe"], root=ROOT, out_dir=Path(tmp), three_prime_n=2, cpu=4)
        live = {g["genome_id"]: g for g in res["genomes"]}
        n = diff = 0
        for x in csv.DictReader(open(table, encoding="utf-8"), delimiter="\t"):
            g = live.get(x["genome_id"])
            if g is None or x["f_mismatch"] == "":
                continue
            n += 1
            for k in ("f_mismatch", "r_mismatch", "p_mismatch"):
                a, b = g[k], x[k]
                if (a is None and b != "") or (a is not None and str(a) != b):
                    diff += 1
        # Same population on both sides: the live analysis counts EVERY genome, so compare with the pipeline's all-genome id list;
        # summary.tsv counts the primary host scope, so compare that against the live rows restricted to the same scope.
        ids_file = ROOT / "output" / "primers" / pid / "final" / "validated_ids.txt"
        v_pipe = len([ln for ln in ids_file.read_text(encoding="utf-8").splitlines() if ln.strip()]) if ids_file.exists() else None
        v_live = res["summary"].get("validated")
        scope_note = ""
        scope_ok = True
        if pid in summary and summary[pid].get("host_scope"):
            import host_scope as H  # noqa: PLC0415
            from mfe_common import probe_inside_amplicon  # noqa: PLC0415
            cfg = (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
            mm = int(re.search(r"^primer_mismatch_pass:\s*(\d+)", cfg, re.M).group(1))
            pm = int(re.search(r"^probe_mismatch:\s*(\d+)", cfg, re.M).group(1))
            seg = summary[pid]["segment"].lower()
            classes = H.load_host_classes(ROOT / "data" / "raw" / f"{seg}_metadata.tsv", {"accession": "accession"})
            keep = H.ids_in_scope(sorted(live), classes, summary[pid]["host_scope"])
            sc = H.summarize(keep, live, mm, pm, probe_inside_amplicon)
            scope_ok = sc["validated"] == int(summary[pid]["validated"]) and sc["total_genomes"] == int(summary[pid]["total_genomes"])
            scope_note = f"  [{summary[pid]['host_scope']} scope: live={sc['validated']}/{sc['total_genomes']} summary.tsv={summary[pid]['validated']}/{summary[pid]['total_genomes']}]"
        ok = diff == 0 and v_live == v_pipe and scope_ok
        bad += not ok
        print(f"{pid:6s} genomes={n:5d}  differing F/R/P values={diff}  validated live={v_live} pipeline={v_pipe}{scope_note}  {'OK' if ok else 'DIFFERENT'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
