"""Full-precision binding records from `mfeprimer -j` reruns of the synth.py experiments (work/exp_*/run_*/oj.txt.json).

Each amplicon gives an F and an R binding (cut > 0: mfeprimer printed a shortened duplex, the primer's 5' part
cut off after a run of mismatches): primer (5'->3'), site in primer sense, the template base beyond the
primer's 5' end (l5) and beyond its 3' end (r3) in primer sense ('-' at the genome end), Tm and dG.
usage: bind_json.py <out.json> <run dir> [...]
"""
import json
import sys
from pathlib import Path

COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fa(p):
    seqs, cur, buf = {}, None, []
    for l in open(p):
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
    seen, out, bad = set(), [], 0
    for d in sys.argv[2:]:
        d = Path(d)
        js = d / "oj.txt.json"
        if not js.exists():
            continue
        g = read_fa(d / "db.fa")
        for a in json.load(open(js)).get("AmpList") or []:
            seq = g[a["Chr"]]
            n = len(seq)
            F, R = a["F"], a["R"]
            L = len(F["Sseq"])
            s = next((x for x in range(F["Start"] - 2, F["End"] + 3) if seq[x:x + L] == F["Sseq"].upper()), None)
            if s is None:
                bad += 1
                continue
            e = s + L - 1
            s -= len(F["Seq"]["Seq"]) - L                 # a shortened duplex keeps the primer's 3' part
            recs = [("F", F["Seq"]["Seq"].upper(), seq[s:e + 1], seq[s - 1] if s > 0 else "-", seq[e + 1] if e + 1 < n else "-",
                     F["Tm"], F["Dg"], len(F["Seq"]["Seq"]) - len(F["Qseq"]))]
            e2, s2 = R["End"], R["Start"]         # plus-strand site is seq[e2:s2] or seq[e2+1:s2+1]
            L = len(R["Sseq"])
            lo = next((x for x in range(min(e2, s2) - L, max(e2, s2) + 2) if seq[x:x + L] == R["Sseq"].upper()), None)
            if lo is None:
                bad += 1
                continue
            hi = lo + L + len(R["Seq"]["Seq"]) - L
            recs.append(("R", rc(R["Seq"]["Seq"].upper()), rc(seq[lo:hi]), rc(seq[hi]) if hi < n else "-",
                         rc(seq[lo - 1]) if lo > 0 else "-", R["Tm"], R["Dg"], len(R["Seq"]["Seq"]) - len(R["Qseq"])))
            for r in recs:
                if r[:5] not in seen:
                    seen.add(r[:5])
                    out.append(dict(zip(("side", "f", "s", "l5", "r3", "tm", "dg", "cut"), r)))
    json.dump(out, open(sys.argv[1], "w"))
    print(f"{len(out)} distinct bindings, {bad} amplicons not located -> {sys.argv[1]}")


if __name__ == "__main__":
    main()
