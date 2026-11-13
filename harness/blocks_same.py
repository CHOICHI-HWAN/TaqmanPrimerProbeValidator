"""True when two dimer/hairpin reports hold the same entries (each entry's whole text block, numbering aside) and
differ at most in the order of entries with equal dG. usage: blocks_same.py <a> <b>"""
import re
import sys


def blocks(path):
    t = open(path).read()
    t = t[t.index(" List ("):t.index("\nParameters")]
    out = []
    for b in re.split(r"\n\n\n(?=(?:Dimer|Hairpin) \d+:)", t.split("\n", 2)[2].strip("\n")):
        out.append(re.sub(r"^(Dimer|Hairpin) \d+:", r"\1:", b.strip("\n")))
    return out


def dgs(bs):
    return [re.search(r"Delta G = (-?[\d.]+)", b).group(1) if "Delta G" in b else "" for b in bs]


a, b = blocks(sys.argv[1]), blocks(sys.argv[2])
print("order-only" if sorted(a) == sorted(b) and dgs(a) == dgs(b) else ("same" if a == b else "DIFF"))
