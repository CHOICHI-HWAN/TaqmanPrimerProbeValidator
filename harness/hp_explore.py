"""Exploration of mfeprimer hairpin records: decode structure, list features."""
import glob, sys
import dh_parse

def decode(st, seq):
    op = [i for i, c in enumerate(st) if c == "/"]
    cl = [i for i, c in enumerate(st) if c == "\\"]
    pairs = list(zip(op, reversed(cl)))
    return pairs

def recs(pat="../work/dh/hh*.txt"):
    out = []
    for f in sorted(glob.glob(pat)):
        for h in dh_parse.hairpins(open(f).read()):
            h["file"] = f
            h["pairs"] = decode(h["struct"], h["seq"])
            out.append(h)
    return out

if __name__ == "__main__":
    for h in recs()[:40]:
        p = h["pairs"]
        loop = p[-1][1] - p[-1][0] - 1
        print(h["score"], h["tm"], h["dg"], len(p), loop, h["struct"], h["seq"], sep="\t")
