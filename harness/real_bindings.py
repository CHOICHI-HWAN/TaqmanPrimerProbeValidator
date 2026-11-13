"""Turn mfeprimer's real .spec.tsv rows into binding records (joined.tsv layout) for fitting / validation.

Every row gives two bindings (left primer on the minus strand, right primer on the plus strand) with the genome
context read from the database it was run on. Identical bindings are kept once.
usage: real_bindings.py <out.tsv> <spec.tsv> [<spec.tsv> ...]   (the db is taken from the sibling .mfe.log)
Output stays under work/ (contains primer sequences).
"""
import csv
import re
import sys
from pathlib import Path

COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_db(path):
    seqs, cur, buf = {}, None, []
    for l in open(path):
        if l.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf).upper()
            cur, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
    if cur:
        seqs[cur] = "".join(buf).upper()
    return seqs


def main():
    out = sys.argv[1]
    seen = set()
    dbs = {}
    rows = []
    for spec in sys.argv[2:]:
        log = open(spec[: -len(".spec.tsv")] + ".mfe.log").read()
        db = Path(re.search(r"Dir: (\S+)", log).group(1)) / re.search(r" -d (\S+)", log).group(1)
        if db not in dbs:
            dbs.clear()
            dbs[db] = read_db(db)
        g = dbs[db]
        for l in open(spec):
            if l.startswith("#"):
                continue
            r = l.rstrip("\n").split("\t")
            seq = g[r[1]]
            n = len(seq)
            fs, fe = int(r[7]) - 1, int(r[8])
            f = r[9]
            site = seq[fs:fe]
            l5 = seq[fs - 1] if fs > 0 else "-"
            r3 = seq[fe] if fe < n else "-"
            b1 = ("F", f, site, l5, r3, r[10], r[12])
            re_, rs = int(r[14]) - 1, int(r[15])
            rp = rc(r[16])
            rsite = rc(seq[re_:rs])
            rl5 = rc(seq[rs]) if rs < n else "-"
            rr3 = rc(seq[re_ - 1]) if re_ > 0 else "-"
            b2 = ("R", rp, rsite, rl5, rr3, r[17], r[19])
            for b in (b1, b2):
                key = b[:5]
                if key in seen:
                    continue
                seen.add(key)
                rows.append({"side": b[0], "found": "1", "fpSeq": b[1], "site": b[2], "l5": b[3], "r3": b[4],
                             "fpTm": b[5], "fpDg": b[6], "f": b[1]})
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["side", "found", "f", "fpSeq", "site", "l5", "r3", "fpTm", "fpDg"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} distinct bindings -> {out}")


if __name__ == "__main__":
    main()
