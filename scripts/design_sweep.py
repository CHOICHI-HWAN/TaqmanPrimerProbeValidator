#!/usr/bin/env python3
"""Coverage-first primer/probe design sweep for one segment: many parameter combinations, each giving a greedy panel of designs.

    python3 scripts/design_sweep.py s --scope all --jobs 32            # S segment, every non-lab genome, 32 parallel runs
    python3 scripts/design_sweep.py s --calibrate                      # one default run, to see the time

For every combination of (IUPAC codes allowed in a primer, in the probe, mismatches a primer may have, amplicon length range, Tm window) the design engine
(scripts/redesign_coverage.py) picks F / R / probe so that the genomes passing ALL THREE under the pipeline's rule are as many as possible, then adds
designs for the genomes the earlier ones miss. Results go to output/design_sweep/<segment>/results.jsonl (one line per run) and aligned.fa.
The genomes are the population of the segment (data/filtered), minus laboratory-host genomes; --scope human|animal|unknown|all selects the host classes.
"""
import argparse
import itertools
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import host_scope as hs  # noqa: E402
import redesign_coverage as rc  # noqa: E402
from mfe_common import find_mafft, read_fasta_seqs  # noqa: E402


def scope_ids(seg, scope):
    """Accessions of the segment's population that belong to the host scope (laboratory hosts are never in any scope)."""
    keep = set(read_fasta_seqs(ROOT / "data" / "filtered" / f"{seg}_filtered.fasta"))
    import csv
    out = []
    with open(ROOT / "data" / "raw" / f"{seg}_metadata.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["accession"] in keep:
                cls = hs.classify(r.get("host"), r.get("is_lab_host"), r.get("taxid"))
                if hs.in_scope(cls, scope):
                    out.append(r["accession"])
    return sorted(out)


def sample_ids(seg, ids, n):
    """At most n genomes, spread evenly over the lineages (round-robin over a hash order, so the sample is reproducible)."""
    import csv
    import hashlib
    lin = {}
    with open(ROOT / "data" / "filtered" / f"{seg}_filtered_metadata.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lin[r["accession"]] = r.get("lineage") or "NA"
    groups = {}
    for i in ids:
        groups.setdefault(lin.get(i, "NA"), []).append(i)
    for g in groups.values():
        g.sort(key=lambda x: hashlib.sha1(x.encode()).hexdigest())
    out = []
    while len(out) < n and any(groups.values()):
        for k in sorted(groups):
            if groups[k] and len(out) < n:
                out.append(groups[k].pop())
    return sorted(out)


def align(seg, ids, threads, sample=0):
    """MAFFT of the genomes. With a sample (large genomes such as a poxvirus): the sample is aligned to the reference genome in parallel chunks
    (mafft --keeplength --add, the reference's own coordinates), so the design works on a few hundred genomes and is validated on all of them later."""
    out = ROOT / "output" / "design_sweep" / seg
    out.mkdir(parents=True, exist_ok=True)
    if sample and len(ids) > sample:
        ids = sample_ids(seg, ids, sample)
    aln = out / "aligned.fa"
    if aln.exists() and aln.stat().st_size > 0 and (out / "ids.json").exists() and json.loads((out / "ids.json").read_text()) == ids:
        return aln, ids
    seqs = read_fasta_seqs(ROOT / "data" / "filtered" / f"{seg}_filtered.fasta")
    (out / "members.fa").write_text("".join(f">{k}\n{seqs[k]}\n" for k in ids), encoding="utf-8")
    if not sample:
        with open(aln, "w", encoding="utf-8") as fh:
            subprocess.run([find_mafft(), "--auto", "--thread", str(threads), str(out / "members.fa")], stdout=fh, stderr=subprocess.DEVNULL, check=True)
    else:
        ref = ROOT / "data" / "refgenome" / f"{seg}_ref.fasta"
        chunks = [ids[i::max(1, threads)] for i in range(max(1, threads))]
        procs = []
        for ci, ch in enumerate(c for c in chunks if c):
            fa = out / f"chunk{ci}.fa"
            fa.write_text("".join(f">{k}\n{seqs[k]}\n" for k in ch), encoding="utf-8")
            res = open(out / f"chunk{ci}.aln", "w", encoding="utf-8")
            procs.append((subprocess.Popen([find_mafft(), "--quiet", "--keeplength", "--thread", "1", "--add", str(fa), str(ref)], stdout=res, stderr=subprocess.DEVNULL), res, out / f"chunk{ci}.aln"))
        rows = {}
        for pr, res, path in procs:
            pr.wait()
            res.close()
            if pr.returncode:
                raise SystemExit(f"mafft failed on {path.name}")
            cur = None
            for ln in path.read_text(encoding="utf-8").splitlines():
                if ln.startswith(">"):
                    cur = ln[1:].split()[0]
                    rows[cur] = []
                elif cur:
                    rows[cur].append(ln.strip())
        with open(aln, "w", encoding="utf-8") as fh:
            for k in ids:
                fh.write(f">{k}\n{''.join(rows[k])}\n")
    (out / "ids.json").write_text(json.dumps(ids))
    return aln, ids


def run_one(seg, aln, params, top):
    ids, rows = rc.read_aligned(aln)
    t = time.time()
    p = dict(params)
    lp, lq = p.pop("len_primer_max", None), p.pop("len_probe_max", None)      # longer oligos reach a higher Tm in AT-rich conserved regions
    if lp:
        rc.LEN_PRIMER = range(18, int(lp) + 1)
    if lq:
        rc.LEN_PROBE = range(20, int(lq) + 1)
    p["amp"] = tuple(p["amp"])
    p["primer_tm"] = tuple(p["primer_tm"])
    if "probe_above" in p:
        p["probe_above"] = tuple(p["probe_above"])
    res = rc.design_from_alignment(ids, rows, top=top, **p)
    return {"segment": seg, "params": params, "n": res["n"], "status": res["status"], "designs": res["designs"], "seconds": round(time.time() - t, 1)}


def combos(seg, grid_name="wide"):
    if grid_name == "hightm":       # TaqMan-like melting temperatures: primers 56-68 C, probe 6-16 C above them (pipeline rule: primers mm<=1)
        grid = {"ambig": [2, 3, 4, 5, 6], "probe_ambig": [2, 3, 4], "pass_mm": [1], "tp_n": [2], "amp": [(70, 200), (70, 250)],
                "primer_tm": [(56.0, 66.0), (58.0, 68.0)], "probe_above": [(6.0, 14.0), (8.0, 16.0)]}
    elif grid_name == "longoligo":
        grid = {"ambig": [3, 4, 5, 6], "probe_ambig": [2, 3, 4], "pass_mm": [1], "tp_n": [2], "amp": [(70, 200), (70, 250)], "len_primer_max": [30], "len_probe_max": [36],
                "primer_tm": [(54.0, 68.0), (56.0, 68.0)], "probe_above": [(3.0, 16.0), (6.0, 16.0)]}
    elif grid_name == "attm":       # AT-rich genome (poxvirus): primers reach only ~46-60 C even at 25-28 nt, the probe is longer and 2-16 C above them
        grid = {"ambig": [1, 2, 3, 4], "probe_ambig": [1, 2, 3], "pass_mm": [1], "tp_n": [2], "amp": [(70, 200), (70, 300)], "len_primer_max": [28], "len_probe_max": [34],
                "primer_tm": [(46.0, 62.0), (48.0, 64.0)], "probe_above": [(2.0, 14.0), (3.0, 16.0)]}
    elif grid_name == "midtm":
        grid = {"ambig": [3, 4, 5, 6], "probe_ambig": [2, 3, 4], "pass_mm": [1], "tp_n": [2], "amp": [(70, 200), (70, 250)],
                "primer_tm": [(54.0, 66.0), (55.0, 67.0)], "probe_above": [(3.0, 14.0), (5.0, 14.0)]}
    else:
        grid = {"ambig": [1, 2, 3, 4, 5], "probe_ambig": [1, 2, 3, 4], "pass_mm": [1, 2], "tp_n": [2], "amp": [(70, 200), (70, 300)], "primer_tm": [(52.0, 66.0), (50.0, 68.0)]}
    keys = list(grid)
    for vals in itertools.product(*(grid[k] for k in keys)):
        yield dict(zip(keys, vals))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("segment")
    ap.add_argument("--scope", default="all")
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--sample", type=int, default=0, help="design on at most this many genomes (large genomes): aligned to the reference in parallel chunks")
    ap.add_argument("--grid", default="wide", choices=("wide", "hightm", "midtm", "longoligo", "attm"))
    a = ap.parse_args()
    ids = scope_ids(a.segment, a.scope)
    print(f"{a.segment}: {len(ids)} genomes in scope '{a.scope}'", flush=True)
    aln, ids = align(a.segment, ids, threads=min(32, a.jobs * 2), sample=a.sample)
    print(f"{a.segment}: aligned {len(ids)} genomes", flush=True)
    out = ROOT / "output" / "design_sweep" / a.segment / ("calibrate.jsonl" if a.calibrate else ("results.jsonl" if a.grid == "wide" else f"results_{a.grid}.jsonl"))
    todo = [dict(ambig=2, probe_ambig=2, pass_mm=1, tp_n=2, amp=(70, 200), primer_tm=(52.0, 66.0))] if a.calibrate else list(combos(a.segment, a.grid))
    done = set()
    if out.exists() and not a.calibrate:
        done = {json.dumps(json.loads(l)["params"], sort_keys=True) for l in out.read_text().splitlines() if l.strip()}
    todo = [p for p in todo if json.dumps(json.loads(json.dumps(p)), sort_keys=True) not in done]
    print(f"{len(todo)} runs to do ({len(done)} already done), {a.jobs} in parallel", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex, open(out, "a", encoding="utf-8") as fh:
        futs = {ex.submit(run_one, a.segment, str(aln), p, a.top): p for p in todo}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                r = f.result()
            except Exception as exc:  # noqa: BLE001
                print("run failed:", futs[f], exc, flush=True)
                continue
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            n = max(r["n"], 1)
            single = max((d["est_cov_pct"] for d in r["designs"]), default=0)
            cum = round(100.0 * sum(d["est_added"] for d in r["designs"]) / n, 1)
            print(f"[{i}/{len(todo)}] {r['seconds']}s  amb{r['params']['ambig']}/{r['params']['probe_ambig']} mm{r['params']['pass_mm']}  best single {single}%  panel of {len(r['designs'])}: {cum}%  (elapsed {time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
