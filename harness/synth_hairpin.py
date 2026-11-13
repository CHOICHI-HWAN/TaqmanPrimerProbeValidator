"""Controlled black-box experiments against `mfeprimer hairpin` (JSON output: full-precision Tm/dG).

Each oligo = 5'flank + stem + loop + rc(stem) + 3'flank, flank columns facing each other made non-complementary
(or, with --mode mm, stems carrying chosen internal mismatches). mfeprimer reports at most one hairpin per oligo.
usage: synth_hairpin.py <outdir> [--n 400] [--seed 1] [--mode perfect|mm|random] [extra mfeprimer flags after --]
"""
import os
import argparse
import json
import random
import subprocess
from pathlib import Path

MFE = os.environ.get("MFEPRIMER", "mfeprimer")
C = {"A": "T", "T": "A", "G": "C", "C": "G"}


def rc(s):
    return "".join(C[b] for b in reversed(s))


def rnd(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def noncomp(rng, b):
    return rng.choice([x for x in "ACGT" if x != C[b]])


def design(rng, mode):
    if mode == "free":         # unstructured random oligos: weak hairpins, many equal-score alternatives
        return rnd(rng, rng.randint(18, 36))
    P = rng.randint(4, 12)
    L = rng.randint(3, 14)
    a, b = rng.randint(0, 6), rng.randint(0, 6)
    if mode == "dg":           # no flanks: the stem runs from end to end
        a = b = 0
        P = rng.randint(5, 12)
    if mode == "dgL":          # long loops
        a = b = 0
        P = rng.randint(8, 14)
        L = rng.randint(10, 34)
        loop = rnd(rng, L)
    stem = rnd(rng, P)
    loop = rnd(rng, L)
    right = list(rc(stem))
    if mode == "dgL":
        loop = "".join(rng.choice("AC") for _ in range(L))   # cannot pair inside the loop
    if mode == "mm":           # mismatches inside the stem (on the 3' arm)
        for k in rng.sample(range(P), rng.randint(1, min(3, P - 2))):
            right[P - 1 - k] = noncomp(rng, stem[k])
    left5 = rnd(rng, a)
    r3 = []
    for k in range(b):        # column k outside the stem: left5[a-1-k] faces r3[k]
        r3.append(noncomp(rng, left5[a - 1 - k]) if k < a and mode != "random" else rng.choice("ACGT"))
    return left5 + stem + loop + "".join(right) + "".join(r3)


def run(seqs, d, extra):
    d.mkdir(parents=True, exist_ok=True)
    (d / "in.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
    subprocess.run([MFE, "hairpin", "-i", "in.fa", "-o", "out", "-j", "-s", "1", "-d", "100", "-c", "1", *extra],
                   cwd=d, check=True, capture_output=True)
    js = json.load(open(d / "out.json"))
    got = {h["Seq"]["ID"]: h for h in js["HairpinList"] or []}
    rows = []
    for i, s in enumerate(seqs):
        h = got.get(f"s{i}")
        rows.append({"id": f"s{i}", "seq": s, **({"score": h["Score"], "tm": h["Tm"], "dg": h["Dg"], "struct": h["Aseq"]}
                                                if h else {})})
    with open(d / "rows.json", "w") as fh:
        json.dump(rows, fh)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mode", default="perfect")
    ap.add_argument("--extra", default="")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    rows = run([design(rng, a.mode) for _ in range(a.n)], Path(a.outdir), a.extra.split())
    print(f"{a.outdir}: {sum('score' in r for r in rows)}/{len(rows)} hairpins")


if __name__ == "__main__":
    main()
