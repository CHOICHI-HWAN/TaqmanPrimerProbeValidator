#!/usr/bin/env python3
"""One row per primer set: how often do mfeprimer and the whole-genome scan agree on the primer-pair verdict?"""
import csv
from pathlib import Path


def main():
    out = []
    for path in sorted(snakemake.input.comparisons):
        set_id = Path(path).parent.parent.name
        rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
        n = len(rows)
        both_pass = sum(1 for r in rows if r["verdict_mfe"] == "1" and r["verdict_own"] == "1")
        both_fail = sum(1 for r in rows if r["verdict_mfe"] == "0" and r["verdict_own"] == "0")
        mfe_only = sum(1 for r in rows if r["verdict_mfe"] == "1" and r["verdict_own"] == "0")
        own_only = sum(1 for r in rows if r["verdict_mfe"] == "0" and r["verdict_own"] == "1")
        scan_used = sum(1 for r in rows if r["source"] == "genome-scan")
        pm_only = sum(1 for r in rows if r["probe_verdict_mfe"] == "1" and r["probe_verdict_own"] == "0")
        po_only = sum(1 for r in rows if r["probe_verdict_mfe"] == "0" and r["probe_verdict_own"] == "1")
        p_scan = sum(1 for r in rows if r["probe_source"] == "genome-scan")
        out.append([set_id, n, both_pass, both_fail, mfe_only, own_only, scan_used,
                    round(100 * (both_pass + both_fail) / n, 2) if n else "", pm_only, po_only, p_scan])
    with open(snakemake.output.table, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["set_id", "genomes", "both_pass", "both_fail", "mfeprimer_only_pass", "scan_only_pass",
                    "pairs_taken_from_scan", "agreement_pct", "probe_mfeprimer_only_perfect",
                    "probe_scan_only_perfect", "probes_taken_from_scan"])
        w.writerows(out)


main()
