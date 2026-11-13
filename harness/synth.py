"""Controlled black-box experiments against mfeprimer's primer search.

One run = one random F/R pair + many synthetic genome records. Every record carries one amplicon whose F binding site
is a chosen variant (mismatches at chosen positions, chosen flanking bases); the R site is always perfect. Running
mfeprimer on that database and joining its .spec.tsv back to the design tells exactly how each variant was scored
(found or not, fpTm, fpDg).

usage: python3 synth.py <experiment> <outdir> [--seed N] [--pairs N] [--mm 2]
  experiment: flank | mm1 | mm2 | mm3p | random
Writes <outdir>/run_<i>/{db.fa,pr.fa,design.tsv,out.txt*} and <outdir>/joined.tsv
"""
import argparse
import csv
import itertools
import os
import random
import subprocess
from pathlib import Path

MFE = os.environ.get("MFEPRIMER_BIN",
                     os.environ.get("MFEPRIMER", "mfeprimer"))
COMP = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n, gc=0.5):
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n))


def mutate(seq, positions, rng, subs=None):
    """Substitute seq at the given 0-based positions (subs: explicit replacement bases, else random different base)."""
    s = list(seq)
    for j, p in enumerate(positions):
        s[p] = subs[j] if subs else rng.choice([b for b in "ACGT" if b != s[p]])
    return "".join(s)


def variants(exp, f, rng, mm):
    """Yield (label, site, left_flank, right_flank) for the F site. Positions are 0-based from the primer's 5' end."""
    n = len(f)
    if exp == "flank":
        for l5, r3 in itertools.product("ACGT", repeat=2):
            for rep in range(3):
                yield f"flank:{l5}{r3}", f, l5, r3
    elif exp == "mm1":            # one mismatch at every position, every substitution
        for p in range(n):
            for b in "ACGT":
                if b != f[p]:
                    yield f"mm1:{p}:{b}", mutate(f, [p], rng, [b]), None, None
    elif exp == "mm2":            # two mismatches, random positions anywhere
        for _ in range(200):
            ps = sorted(rng.sample(range(n), 2))
            yield f"mm2:{ps[0]},{ps[1]}", mutate(f, ps, rng), None, None
    elif exp == "mm3p":           # mismatches only in the 3' 9-nt window: 1..mm+1 of them
        for k in range(1, mm + 2):
            for _ in range(60):
                ps = sorted(rng.sample(range(n - 9, n), k))
                yield f"mm3p:{k}:" + ",".join(map(str, ps)), mutate(f, ps, rng), None, None
    elif exp == "tandem":         # 2-3 adjacent mismatches outside the 3' window, every substitution pattern sampled
        for _ in range(240):
            k = rng.choice([2, 2, 3])
            p = rng.randint(1, n - 9 - k)
            yield f"tandem:{k}:{p}", mutate(f, list(range(p, p + k)), rng), None, None
    elif exp == "t5":             # 5'-terminal mismatches (position 0, and 0+1) with every flanking base
        for b in "ACGT":
            if b != f[0]:
                for l5 in "ACGT":
                    yield f"t5:0:{b}{l5}", mutate(f, [0], rng, [b]), l5, None
                    yield f"t5:01:{b}{l5}", mutate(mutate(f, [0], rng, [b]), [1], rng), l5, None
                    yield f"t5:1:{b}{l5}", mutate(f, [1], rng), l5, None
    elif exp == "iupac":          # genome ambiguity codes (IUPAC_CODES, default N) in the site, its flank and the window
        codes = os.environ.get("IUPAC_CODES", "N")
        for _ in range(300):
            c = rng.choice(codes)
            where = rng.choice(["out", "out", "out", "out2", "win", "l5", "first"])
            site = f
            l5 = None
            if where == "out":
                site = mutate(f, [rng.randint(1, n - 10)], rng, [c])
            elif where == "out2":
                ps = sorted(rng.sample(range(1, n - 9), 2))
                site = mutate(f, ps, rng, [c, rng.choice(codes)])
            elif where == "win":
                site = mutate(f, [rng.randint(n - 9, n - 1)], rng, [c])
            elif where == "first":
                site = mutate(f, [0], rng, [c])
            else:
                l5 = c
            yield f"iupac:{where}:{c}", site, l5, None
    elif exp == "iupac":          # genome ambiguity codes (IUPAC_CODES, default N) in the site, its flank and the window
        codes = os.environ.get("IUPAC_CODES", "N")
        for _ in range(300):
            c = rng.choice(codes)
            where = rng.choice(["out", "out", "out", "out2", "win", "l5", "first"])
            site = f
            l5 = None
            if where == "out":
                site = mutate(f, [rng.randint(1, n - 10)], rng, [c])
            elif where == "out2":
                ps = sorted(rng.sample(range(1, n - 9), 2))
                site = mutate(f, ps, rng, [c, rng.choice(codes)])
            elif where == "win":
                site = mutate(f, [rng.randint(n - 9, n - 1)], rng, [c])
            elif where == "first":
                site = mutate(f, [0], rng, [c])
            else:
                l5 = c
            yield f"iupac:{where}:{c}", site, l5, None
    elif exp == "tandemN":        # runs of 2-4 consecutive mismatches outside the window, one of them an N
        for _ in range(300):
            k = rng.choice([2, 3, 3, 4])
            p = rng.randint(1, n - 9 - k)
            site = mutate(f, list(range(p, p + k)), rng)
            q = rng.randint(p, p + k - 1)
            site = site[:q] + "N" + site[q + 1:]
            yield f"tandemN:{k}:{p}:{q}", site, None, None
    elif exp == "fill":           # 5'-terminal stacks with N (xN, Nx, NN), s0 = N with an N flank, N beyond the 3' end
        for _ in range(300):
            pat = rng.choice(["xN", "Nx", "NN", "N0", "N0", "r3"])
            site = list(f)
            l5 = rng.choice("ACGTN")
            r3 = None
            if pat == "r3":
                r3 = "N"
                site = list(mutate(f, [rng.randint(1, n - 10)], rng)) if rng.random() < .5 else site
            elif pat == "N0":
                site[0] = "N"
                l5 = "N"
            else:
                for i, c in enumerate(pat):
                    site[i] = "N" if c == "N" else rng.choice([b for b in "ACGT" if b != f[i]])
            yield f"fill:{pat}", "".join(site), l5, r3
    elif exp in ("edge", "random"):  # 0..6 mismatches anywhere (edge: the site touches the genome end)
        for _ in range(300):
            k = rng.randint(0, 6)
            ps = sorted(rng.sample(range(n), k))
            yield f"rand:{k}:" + ",".join(map(str, ps)), mutate(f, ps, rng), None, None
    else:
        raise SystemExit(f"unknown experiment {exp}")


def one_run(exp, rdir, rng, mm, side="F"):
    rdir.mkdir(parents=True, exist_ok=True)
    f = rand_seq(rng, rng.randint(18, 26), rng.uniform(0.35, 0.65))
    r = rand_seq(rng, rng.randint(18, 26), rng.uniform(0.35, 0.65))
    mid = rand_seq(rng, rng.randint(40, 120))
    r_flank5, r_flank3 = rand_seq(rng, 30), rand_seq(rng, 30)
    rows = []
    with open(rdir / "db.fa", "w") as db:
        test = f if side == "F" else r          # the primer whose site is varied
        for i, (label, site, l5, r3) in enumerate(variants(exp, test, rng, mm)):
            l5 = l5 or rng.choice("ACGT")
            r3 = r3 or rng.choice("ACGT")
            rec = f"rec{i}"
            edge = exp == "edge"
            if side == "F":
                left = "" if edge else rand_seq(rng, 29) + l5
                genome = left + site + r3 + mid + rc(r) + r_flank3
                start = len(left) + 1
                if edge:
                    l5 = "-"
            else:                               # R binds the plus strand: on the minus strand it reads l5 + site + r3
                left = rand_seq(rng, 30) + f + mid + rc(r3)
                genome = left + rc(site) + ("" if edge else rc(l5) + r_flank3)
                start = len(left) + 1
                if edge:
                    l5 = "-"
            db.write(f">{rec}\n{genome}\n")
            mism = [p for p in range(len(test)) if site[p] != test[p]]
            rows.append({"rec": rec, "label": label, "side": side, "f": test, "r": r if side == "F" else f, "site": site,
                         "l5": l5, "r3": r3, "fp_start": start, "n_mm": len(mism), "mm_pos": ",".join(map(str, mism))})
    with open(rdir / "pr.fa", "w") as pr:
        pr.write(f">X_F\n{f}\n>X_R\n{r}\n")
    with open(rdir / "design.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    subprocess.run([MFE, "index", "-i", "db.fa"], cwd=rdir, check=True, capture_output=True)
    for p in rdir.glob("out.txt*"):
        p.unlink()
    subprocess.run([MFE, "-d", "db.fa", "-i", "pr.fa", "-o", "out.txt", "--fasta", "--virus", "--misMatch", str(mm),
                    "--cpu", "4"], cwd=rdir, check=True, capture_output=True)
    return rows


def read_spec(path):
    out = {}
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    head = ["name", "chrom", "ampStart", "ampEnd", "ampGC", "ampSize", "fpName", "fpStart", "fpEnd", "fpSeq", "fpTm",
            "fpGC", "fpDg", "rpName", "rpEnd", "rpStart", "rpSeq", "rpTm", "rpGC", "rpDg", "note"]
    for l in lines:
        v = dict(zip(head, l.rstrip("\n").split("\t")))
        out.setdefault(v["chrom"], []).append(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment")
    ap.add_argument("outdir")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=3)
    ap.add_argument("--mm", type=int, default=2)
    ap.add_argument("--side", choices=["F", "R"], default="F")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = Path(a.outdir)
    joined = []
    for i in range(a.pairs):
        rdir = out / f"run_{i}"
        rows = one_run(a.experiment, rdir, rng, a.mm, a.side)
        spec = read_spec(rdir / "out.txt.spec.tsv")
        for row in rows:
            hits = [h for h in spec.get(row["rec"], []) if h["fpName"].startswith("X_F") and h["rpName"].startswith("X_R")
                    and (int(h["fpStart"]) if a.side == "F" else min(int(h["rpStart"]), int(h["rpEnd"]))) == row["fp_start"]]
            h = hits[0] if hits else {}
            t, o = ("fp", "rp") if a.side == "F" else ("rp", "fp")   # columns of the tested / the other primer
            joined.append({**row, "run": i, "found": int(bool(hits)), "n_hits": len(spec.get(row["rec"], [])),
                           "fpName": h.get(t + "Name", ""), "fpSeq": (h.get("fpSeq", "") if a.side == "F" else rc(h.get("rpSeq", ""))),  # rpSeq is printed plus-strand
                           "fpTm": h.get(t + "Tm", ""),
                           "fpDg": h.get(t + "Dg", ""), "rpTm": h.get(o + "Tm", ""), "rpDg": h.get(o + "Dg", "")})
    with open(out / "joined.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(joined[0]), delimiter="\t")
        w.writeheader()
        w.writerows(joined)
    found = sum(r["found"] for r in joined)
    print(f"{a.experiment}: {len(joined)} records, {found} found as the designed F/R amplicon -> {out / 'joined.tsv'}")


if __name__ == "__main__":
    main()
