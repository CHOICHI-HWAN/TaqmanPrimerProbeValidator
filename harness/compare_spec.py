"""Compare a primer-engine run with mfeprimer's outputs for the same db / primers.

usage: compare_spec.py <mfeprimer out.txt> <engine out.txt> [--show N]
Rows are matched on (chrom, fpName, fpStart, rpName, rpStart); the Amp_ number is ignored (mfeprimer's order is not
deterministic). Prints counts of missing / extra rows by primer-pair kind and the field differences of matched rows.
"""
import argparse
import collections
import re

COLS = ["name", "chrom", "ampStart", "ampEnd", "ampGC", "ampSize", "fpName", "fpStart", "fpEnd", "fpSeq", "fpTm",
        "fpGC", "fpDg", "rpName", "rpEnd", "rpStart", "rpSeq", "rpTm", "rpGC", "rpDg", "note"]
NUM = {"fpTm", "fpDg", "rpTm", "rpDg", "ampGC", "fpGC", "rpGC"}


def rows(path):
    out = {}
    for l in open(path + ".spec.tsv"):
        if l.startswith("#"):
            continue
        r = dict(zip(COLS, l.rstrip("\n").split("\t")))
        out[(r["chrom"], r["fpName"], r["fpStart"], r["rpName"], r["rpStart"])] = r
    return out


def kind(r):
    f = re.sub(r"\..*", "", r["fpName"])[-1]
    rr = re.sub(r"\..*", "", r["rpName"])[-1]
    return f + rr


def failed(path):
    try:
        return [l.split()[0] for l in open(path + ".virus.Failed.txt") if l.strip()]
    except FileNotFoundError:
        return []


def coverage(path):
    m = re.search(r"Coverage .*", open(path).read())
    return m.group(0) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("--show", type=int, default=3)
    a = ap.parse_args()
    A, B = rows(a.ref), rows(a.new)
    miss = [k for k in A if k not in B]
    extra = [k for k in B if k not in A]
    both = [k for k in A if k in B]
    diffs = collections.Counter()
    worst = collections.defaultdict(float)
    for k in both:
        for c in COLS[1:]:
            x, y = A[k][c], B[k][c]
            if x == y:
                continue
            if c in NUM:
                d = abs(float(x) - float(y))
                worst[c] = max(worst[c], d)
                if d > 0.0101:
                    diffs[c + ">0.01"] += 1
                else:
                    diffs[c + "±0.01"] += 1
            else:
                diffs[c] += 1
    ck = lambda r: (r["chrom"], r["ampStart"], r["ampEnd"], r["fpSeq"], r["rpSeq"])   # the product itself, names aside
    CA, CB = {ck(r) for r in A.values()}, {ck(r) for r in B.values()}
    print(f"products: ref {len(CA)}  new {len(CB)}  missing {len(CA - CB)}  extra {len(CB - CA)}")
    fa, fb = failed(a.ref), failed(a.new)
    print(f"rows: ref {len(A)}  new {len(B)}  matched {len(both)}  missing {len(miss)}  extra {len(extra)}")
    print("  missing by kind:", dict(collections.Counter(kind(A[k]) for k in miss)))
    print("  extra by kind:  ", dict(collections.Counter(kind(B[k]) for k in extra)))
    print("  field diffs:", dict(diffs), " worst:", {k: round(v, 3) for k, v in worst.items()})
    print(f"failed: same set {set(fa) == set(fb)}  same order {fa == fb}  ref {len(fa)} new {len(fb)}")
    print("coverage ref:", coverage(a.ref))
    print("coverage new:", coverage(a.new))
    with open(a.new + ".missing.tsv", "w") as fh:
        for k in miss:
            fh.write("\t".join(A[k][c] for c in COLS) + "\n")
    with open(a.new + ".extra.tsv", "w") as fh:
        for k in extra:
            fh.write("\t".join(B[k][c] for c in COLS) + "\n")
    for k in miss[:a.show]:
        print("MISSING", "\t".join(A[k][c] for c in COLS[1:9] + ["fpTm", "rpName", "rpEnd", "rpStart", "rpTm"]))
    for k in extra[:a.show]:
        print("EXTRA  ", "\t".join(B[k][c] for c in COLS[1:9] + ["fpTm", "rpName", "rpEnd", "rpStart", "rpTm"]))


if __name__ == "__main__":
    main()
