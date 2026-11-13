#!/usr/bin/env python3
"""Picks a diverse short list of candidates from output/design_sweep/<segment>/results.jsonl and (optionally) registers them as primer sets.

    python3 scripts/design_sweep_select.py s --mm 1 --top 12                 # print the short list
    python3 scripts/design_sweep_select.py s --mm 1 --top 12 --register      # add the sets to config/primers.tsv (source AUTO-DESIGN)

Ranking: estimated coverage first, then fewer IUPAC codes (a degenerate oligo is a mixture, so fewer is better), an amplicon length that suits qPCR (70-200 bp),
and balanced F/R Tm. One candidate per distinct probe, so the list covers different places of the genome.
"""
import argparse
import csv
import datetime
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
IUPAC_N = set("RYSWKMBDHVN")
MIN_PROBE_ABOVE = 2.0          # the probe Tm must be higher than both primers (required rule)


def prefix():
    """Short name of this project's virus (from its profile), used in the names of the registered sets."""
    import re
    m = re.search(r'^profile:\s*"([^"]+)"', (ROOT / "config" / "config.yaml").read_text(encoding="utf-8"), re.M)
    return json.loads((ROOT / m.group(1)).read_text(encoding="utf-8"))["short_name"]


def degeneracy(seq):
    return sum(1 for c in seq if c in IUPAC_N)


def load(seg, mm):
    runs = []
    for f in sorted((ROOT / "output" / "design_sweep" / seg).glob("results*.jsonl")):
        runs += [json.loads(l) for l in f.read_text().splitlines() if l.strip()]
    return [r for r in runs if mm is None or r["params"]["pass_mm"] == mm]


def candidates(seg, mm):
    best = {}
    for r in load(seg, mm):
        for d in r["designs"]:
            k = (d["F"], d["P"], d["R"])
            if k not in best or d["est_cov"] > best[k]["est_cov"]:
                best[k] = dict(d, n=r["n"], params=r["params"])
    return list(best.values())


def score(d):
    deg = degeneracy(d["F"]) + degeneracy(d["R"]) + degeneracy(d["P"])
    amp_pen = 0.0 if 70 <= d["amplicon_len"] <= 200 else min(5.0, (d["amplicon_len"] - 200) / 20.0 if d["amplicon_len"] > 200 else 2.0)
    dtm = abs(d["F_tm"] - d["R_tm"])
    return d["est_cov_pct"] - 0.35 * deg - amp_pen - 0.3 * max(0.0, dtm - 2.0)


def baseline(seg):
    """Best validated coverage (%) of the existing (not auto-designed) sets of the segment, from output/summary/summary.tsv; 0 when unknown."""
    f = ROOT / "output" / "summary" / "summary.tsv"
    if not f.exists():
        return 0.0
    best = 0.0
    for r in csv.DictReader(f.open(encoding="utf-8"), delimiter="\t"):
        if r["segment"].lower() == seg.lower() and r["source"] != "AUTO-DESIGN" and r.get("validated_pct"):
            best = max(best, float(r["validated_pct"]))
    return best


def prune(seg):
    """Removes registered AUTO-DESIGN sets whose validated coverage is not higher than the best existing set of the segment."""
    f = ROOT / "output" / "summary" / "summary.tsv"
    base = baseline(seg)
    val = {}
    for r in csv.DictReader(f.open(encoding="utf-8"), delimiter="\t"):
        if r["segment"].lower() == seg.lower() and r["source"] == "AUTO-DESIGN" and r.get("validated_pct"):
            val[r["set_id"]] = max(val.get(r["set_id"], 0.0), float(r["validated_pct"]))
    drop = {k for k, v in val.items() if v <= base}
    manifest = ROOT / "config" / "primers.tsv"
    rows = list(csv.DictReader(manifest.open(encoding="utf-8"), delimiter="\t"))
    with manifest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(r for r in rows if r["primer_id"] not in drop)
    return sorted(drop), base


def shortlist(seg, mm, top, min_tm=0.0):
    base = baseline(seg)
    ranked = sorted((d for d in candidates(seg, mm) if min(d["F_tm"], d["R_tm"]) >= min_tm and d["P_tm"] >= max(d["F_tm"], d["R_tm"]) + MIN_PROBE_ABOVE), key=lambda d: -score(d))
    seen, out = set(), []
    for d in ranked:
        if d["P"] in seen or d["est_cov_pct"] < 80.0 or d["est_cov_pct"] <= base:
            continue
        seen.add(d["P"])
        out.append(d)
        if len(out) >= top:
            break
    return out


def register(seg, items):
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    manifest = ROOT / "config" / "primers.tsv"
    reg_file = ROOT / "config" / "id_registry.json"
    reg = json.loads(reg_file.read_text()) if reg_file.exists() else {"primer_next": 1}
    rows = list(csv.DictReader(manifest.open(encoding="utf-8"), delimiter="\t"))
    fields = list(rows[0].keys()) if rows else ["primer_id", "primer_name", "segment", "forward_primer", "reverse_primer", "probe", "fluorophore", "source", "remark", "proceed", "add_date"]
    have = {(r["forward_primer"], r["reverse_primer"], r["probe"]) for r in rows}
    n = 0
    taken = sum(1 for r in rows if r["primer_name"].startswith(f"{prefix()}-{seg.upper()}-AD"))
    for d in items:
        if (d["F"], d["R"], d["P"]) in have:
            continue
        taken += 1
        i = taken
        pid = f"P{reg['primer_next']:03d}"
        reg["primer_next"] += 1
        rows.append({"primer_id": pid, "primer_name": f"{prefix()}-{seg.upper()}-AD{i:02d} ★", "segment": seg.upper(), "forward_primer": d["F"], "reverse_primer": d["R"], "probe": d["P"],
                     "fluorophore": "FAM", "source": "AUTO-DESIGN", "remark": f"coverage-first design, est {d['est_cov_pct']}% (ambig {d['params']['ambig']}/{d['params']['probe_ambig']} mm{d['params']['pass_mm']}), Tm F/P/R {d['F_tm']}/{d['P_tm']}/{d['R_tm']}",
                     "proceed": "True", "add_date": datetime.date.today().isoformat()})
        n += 1
    with manifest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    reg_file.write_text(json.dumps(reg, indent=1))
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("segment")
    ap.add_argument("--mm", type=int, default=1)
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--min-tm", type=float, default=0.0, help="both primers at least this Tm (the engine's estimate)")
    ap.add_argument("--register", action="store_true")
    ap.add_argument("--prune", action="store_true", help="remove registered auto-designed sets not better than the best existing set (needs a finished pipeline run)")
    a = ap.parse_args()
    if a.prune:
        dropped, base = prune(a.segment)
        print(f"existing best {base}%: removed {len(dropped)} auto-designed sets {dropped}")
        return
    items = shortlist(a.segment, a.mm, a.top, a.min_tm)
    for i, d in enumerate(items, 1):
        print(f"{i:2d}. {d['est_cov']}/{d['n']} = {d['est_cov_pct']}%  amp {d['amplicon_len']}bp  degeneracy {degeneracy(d['F']) + degeneracy(d['R']) + degeneracy(d['P'])}  Tm F/P/R {d['F_tm']}/{d['P_tm']}/{d['R_tm']}")
        print(f"      F {d['F']}  P {d['P']}  R {d['R']}")
    if a.register:
        print(f"registered {register(a.segment, items)} sets in config/primers.tsv")


if __name__ == "__main__":
    main()
