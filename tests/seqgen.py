"""Synthetic nucleotide strings for tests.

Nothing in this repository may be a real primer, probe, standard or genome sequence (they are confidential and stay out of git), so
tests build the strings they need here: deterministic for a given (length, seed), no homopolymer runs of 3+ so a scan over them is
unambiguous. `fill()` swaps "@@S<k>:<n>@@" placeholders inside captured report text for such strings.
"""
import random
import re


def synth(n, seed):
    rng = random.Random(seed)
    out = []
    while len(out) < n:
        c = rng.choice("ACGT")
        if len(out) >= 2 and out[-1] == out[-2] == c:
            continue
        out.append(c)
    return "".join(out)


def comp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))


def fill(text):
    return re.sub(r"@@S(\d+):(\d+)@@", lambda m: synth(int(m.group(2)), 100 + int(m.group(1))), text)
