"""Compare `pe hairpin` with `mfeprimer hairpin` reports: the oligo table, every entry (name, Score, Tm, dG, structure,
sequence, in order) and the pipeline's per-role summary (worst dG / its Tm / score, dG <= -9 count;
run_dimer_hairpin.py). usage: compare_hairpin.py <mfeprimer report> <engine report>"""
import re
import sys

import dh_parse as P

TABLE = re.compile(r"^(\S+)\s+([A-Za-z]+)\s+(\d+)\s+([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)$", re.M)


def entries(text):
    return [(h["name"], h["score"], h["tm"], h["dg"], h["struct"], h["seq"]) for h in P.hairpins(text)]


def summary(es):
    by = {}
    for name, score, tm, dg, *_ in es:
        s = by.setdefault(name.split(".")[0], [None, None, None, 0])
        if dg <= -9.0:
            s[3] += 1
        if s[0] is None or dg < s[0]:
            s[0], s[1], s[2] = dg, tm, score
    return by


def main():
    a, b = open(sys.argv[1]).read(), open(sys.argv[2]).read()
    ta, tb = TABLE.findall(a.split("Hairpin List")[0]), TABLE.findall(b.split("Hairpin List")[0])
    ea, eb = entries(a), entries(b)
    n = re.search(r"Hairpin List \((\d+)\)", a).group(1), re.search(r"Hairpin List \((\d+)\)", b).group(1)
    same_entries = len(ea) == len(eb) and all(x[0] == y[0] and x[1] == y[1] and abs(x[2] - y[2]) < 0.0101 and
                                              abs(x[3] - y[3]) < 0.0101 and x[4:] == y[4:] for x, y in zip(ea, eb))
    exact = ea == eb
    same_set = sorted(ea) == sorted(eb)
    same_sum = summary(ea) == summary(eb)
    print(f"{sys.argv[1]}\ttable={'same' if ta == tb else 'DIFF'}\tcount={n[0]}/{n[1]}\t"
          f"entries={'exact' if exact else ('same' if same_entries else ('order' if same_set else 'DIFF'))}\tsummary={'same' if same_sum else 'DIFF'}")
    if not same_entries and not same_set:
        for x, y in zip(ea + [None] * len(eb), eb + [None] * len(ea)):
            if x != y:
                print("   mfe", x, "\n   pe ", y)
                break


if __name__ == "__main__":
    main()
