"""Compare `pe search` output with `mfeprimer search` output (rows keyed by chrom, start, end, query, strand).
usage: compare_search.py <mfeprimer out> <engine out> [--quiet]"""
import sys


def rows(path):
    out = {}
    for l in open(path):
        if l.startswith("#") or not l.strip():
            continue
        c = l.rstrip("\n").split("\t")
        out[(c[0], c[1], c[2], c[3], c[5])] = c
    return out


def main():
    A, B = rows(sys.argv[1]), rows(sys.argv[2])
    miss = [k for k in A if k not in B]
    extra = [k for k in B if k not in A]
    diff = {"mm": 0, "seq": 0, "tm>0.01": 0, "dg>0.01": 0}
    worst = 0.0
    for k in A.keys() & B.keys():
        a, b = A[k], B[k]
        diff["mm"] += a[6] != b[6]
        diff["seq"] += (a[9], a[10]) != (b[9], b[10])
        e = abs(float(a[8]) - float(b[8]))
        worst = max(worst, e)
        diff["tm>0.01"] += e > 0.0101
        diff["dg>0.01"] += abs(float(a[7]) - float(b[7])) > 0.0101
    print(f"hits ref {len(A)} new {len(B)} missing {len(miss)} extra {len(extra)} | diffs {diff} worst Tm {worst:.3f}")
    if "--quiet" not in sys.argv:
        for k in miss[:3]:
            print("MISSING", "\t".join(A[k]))
        for k in extra[:3]:
            print("EXTRA  ", "\t".join(B[k]))


if __name__ == "__main__":
    main()
