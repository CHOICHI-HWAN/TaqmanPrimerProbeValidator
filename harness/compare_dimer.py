"""Compare `pe dimer` with `mfeprimer dimer` reports: entries keyed by (A, B, alignment text), Score and dG; plus the
pipeline's own per-role-pair summary (worst dG / its score / dG <= -9 count, run_dimer_hairpin.py).
usage: compare_dimer.py <mfeprimer report> <engine report>"""
import collections
import sys

import dh_parse as P


def load(path):
    out = collections.Counter()
    vals = {}
    for d in P.dimers(open(path).read()):
        key = (d["a_name"], d["b_name"], d["shift"], d["marks"].rstrip())
        out[key] += 1
        vals[key] = (d["score"], d["dg"])
    return out, vals


def summary(path):
    by = {}
    for d in P.dimers(open(path).read()):
        pair = "-".join(sorted((d["a_name"].split(".")[0], d["b_name"].split(".")[0])))
        s = by.setdefault(pair, [None, None, 0])
        if d["dg"] <= -9.0:
            s[2] += 1
        if s[0] is None or d["dg"] < s[0]:
            s[0], s[1] = d["dg"], d["score"]
    return by


def main():
    A, va = load(sys.argv[1])
    B, vb = load(sys.argv[2])
    miss = sum((A - B).values())
    extra = sum((B - A).values())
    dgd = sum(1 for k in A.keys() & B.keys() if abs(va[k][1] - vb[k][1]) > 0.0101 or va[k][0] != vb[k][0])
    sa, sb = summary(sys.argv[1]), summary(sys.argv[2])
    same = sa.keys() == sb.keys() and all(abs(sa[k][0] - sb[k][0]) < 0.0101 and sa[k][1] == sb[k][1] and sa[k][2] == sb[k][2]
                                          for k in sa)
    print(f"entries ref {sum(A.values())} new {sum(B.values())} missing {miss} extra {extra} value-diff {dgd} | "
          f"pipeline summary same {same}")


if __name__ == "__main__":
    main()
