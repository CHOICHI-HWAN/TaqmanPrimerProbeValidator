"""Controlled black-box experiments against `mfeprimer search` (single-oligo search, used for the probe and for
reference-genome binding).

One run = one random oligo + many genome records, each holding one designed site (mismatches at chosen positions,
optional insertion/deletion, either strand, chosen flanks). The hits mfeprimer reports are joined back to the design.
usage: synth_search.py <experiment> <outdir> [--seed N] [--runs N] [--mm 8]
  experiment: mm (0..10 random substitutions) | indel | flank | iupac
"""
import argparse
import csv
import os
import random
import subprocess
from pathlib import Path

MFE = os.environ.get("MFEPRIMER_BIN",
                     os.environ.get("MFEPRIMER", "mfeprimer"))
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n, gc=0.5):
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n))


def design(exp, p, rng):
    n = len(p)
    s = list(p)
    label = exp
    if exp == "mm":
        k = rng.randint(0, 10)
        for i in rng.sample(range(n), k):
            s[i] = rng.choice([b for b in "ACGT" if b != p[i]])
        label = f"mm:{k}"
    elif exp == "indel":
        i = rng.randint(1, n - 2)
        if rng.random() < 0.5:
            s.insert(i, rng.choice("ACGT"))
            label = f"ins:{i}"
        else:
            del s[i]
            label = f"del:{i}"
        for j in rng.sample(range(len(s)), rng.randint(0, 2)):
            s[j] = rng.choice([b for b in "ACGT" if b != s[j]])
    elif exp == "flank":
        for i in rng.sample(range(n), rng.randint(0, 3)):
            s[i] = rng.choice([b for b in "ACGT" if b != p[i]])
    elif exp == "iupac":
        for i in rng.sample(range(n), rng.randint(1, 2)):
            s[i] = rng.choice("NRYKMSWBDHV")
        label = "iupac"
    return label, "".join(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment")
    ap.add_argument("outdir")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--recs", type=int, default=300)
    ap.add_argument("--mm", type=int, default=8)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = Path(a.outdir)
    joined = []
    for run in range(a.runs):
        d = out / f"run_{run}"
        d.mkdir(parents=True, exist_ok=True)
        p = rand_seq(rng, rng.randint(int(os.environ.get("QMIN", 18)), int(os.environ.get("QMAX", 30))), rng.uniform(0.35, 0.65))
        rows = []
        with open(d / "db.fa", "w") as db:
            for i in range(a.recs):
                label, site = design(a.experiment, p, rng)
                l5, r3 = rand_seq(rng, 40), rand_seq(rng, 40)   # plus-strand flanks of the site as the oligo reads it
                strand = rng.choice("+-")
                block = l5 + site + r3
                g = block if strand == "+" else rc(block)
                db.write(f">rec{i}\n{g}\n")
                start = 40 if strand == "+" else 40
                rows.append({"rec": f"rec{i}", "label": label, "oligo": p, "site": site, "strand": strand,
                             "l5": l5[-1], "r3": r3[0], "glen": len(g), "start0": start, "sitelen": len(site)})
        (d / "q.fa").write_text(f">Q\n{p}\n")
        subprocess.run([MFE, "index", "-i", "db.fa"], cwd=d, check=True, capture_output=True)
        if (d / "out.txt").exists():
            (d / "out.txt").unlink()
        subprocess.run([MFE, "search", "db.fa", "q.fa", "out.txt", "-m", str(a.mm)], cwd=d, check=True,
                       capture_output=True)
        hits = {}
        for l in open(d / "out.txt"):
            if l.startswith("#") or not l.strip():
                continue
            c = l.rstrip("\n").split("\t")
            hits.setdefault(c[0], []).append(c)
        for r in rows:
            hs = hits.get(r["rec"], [])
            best = min(hs, key=lambda c: int(c[6])) if hs else None
            joined.append({**r, "run": run, "n_hits": len(hs), "found": int(bool(hs)),
                           **({"h_start": best[1], "h_end": best[2], "h_name": best[3], "h_size": best[4],
                               "h_strand": best[5], "h_mm": best[6], "h_dg": best[7], "h_tm": best[8],
                               "h_q": best[9], "h_hit": best[10]} if best else {})})
    keys = list(joined[0].keys())
    for r in joined:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(out / "joined.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(joined)
    print(f"{a.experiment}: {len(joined)} records, {sum(r['found'] for r in joined)} found -> {out / 'joined.tsv'}")


if __name__ == "__main__":
    main()
