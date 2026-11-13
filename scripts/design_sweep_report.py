#!/usr/bin/env python3
"""Summary of output/design_sweep/<segment>/results.jsonl: the best single sets and the best small panels, with the parameters that gave them.

    python3 scripts/design_sweep_report.py s [--target 90]
Coverage here is the design engine's estimate on the alignment (the pipeline's rule on the aligned genomes); the real number comes from registering the set
and running the pipeline.
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load(seg):
    out = []
    for f in sorted((ROOT / "output" / "design_sweep" / seg).glob("results*.jsonl")):
        out += [json.loads(l) for l in f.read_text().splitlines() if l.strip()]
    return out


def key(d):
    return (d["F"], d["P"], d["R"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("segment")
    ap.add_argument("--target", type=float, default=90.0)
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--mm", type=int, default=None, help="only the runs with this primer mismatch allowance")
    a = ap.parse_args()
    runs = load(a.segment)
    if a.mm is not None:
        runs = [r for r in runs if r["params"]["pass_mm"] == a.mm]
    n = runs[0]["n"]
    print(f"{a.segment.upper()}: {len(runs)} runs, {n} genomes")
    singles = {}
    for r in runs:
        for d in r["designs"]:
            k = key(d)
            if k not in singles or d["est_cov"] > singles[k][0]["est_cov"]:
                singles[k] = (d, r["params"])
    best = sorted(singles.values(), key=lambda x: -x[0]["est_cov"])
    print(f"\nbest single sets ({len(singles)} distinct designs; target {a.target}%):")
    for d, p in best[:a.top]:
        print(f"  {d['est_cov']}/{n} = {d['est_cov_pct']}%  amp {d['amplicon_len']}bp  F {d['F']} | P {d['P']} | R {d['R']}  (ambig {p['ambig']}/{p['probe_ambig']} mm{p['pass_mm']} amp{p['amp']})")
    ok1 = [x for x in best if x[0]["est_cov_pct"] >= a.target]
    print(f"  single sets at or above {a.target}%: {len(ok1)}")
    print("\nbest panels (greedy order within one run):")
    for size in (2, 3):
        panels = []
        for r in runs:
            if len(r["designs"]) >= size:
                cov = sum(d["est_added"] for d in r["designs"][:size])
                panels.append((cov, r))
        panels.sort(key=lambda x: -x[0])
        if panels:
            cov, r = panels[0]
            print(f"  {size} sets: {cov}/{n} = {round(100.0 * cov / n, 1)}%  (ambig {r['params']['ambig']}/{r['params']['probe_ambig']} mm{r['params']['pass_mm']})")
            for d in r["designs"][:size]:
                print(f"     {d['name']} +{d['est_added']}  {d['amplicon_len']}bp  F {d['F']} | P {d['P']} | R {d['R']}")


if __name__ == "__main__":
    main()
