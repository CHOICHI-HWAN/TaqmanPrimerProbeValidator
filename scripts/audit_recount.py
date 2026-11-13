#!/usr/bin/env python3
"""Independent audit of the per-genome mismatch tables.

Re-counts F / R / probe mismatches for every primer set by brute force (every window, both strands, IUPAC-aware,
numpy) and compares them with output/primers/<ID>/alignment/mismatch_table.tsv -- the numbers the dashboard and the
final verdict are built from. Only an UNDER-count by the pipeline (table < brute force: impossible for a real site)
is an error. A higher table value is normal: the pipeline reports the best site that forms a valid amplicon with the
partner primer, which can be worse than the best site anywhere in the genome (reported as "pair-limited"). Run from the project folder:

    python3 scripts/audit_recount.py            # every active set
    python3 scripts/audit_recount.py P001 P003  # chosen sets

Exit status 1 if any count differs. Needs numpy (it is in the conda env). A genome base that is itself an IUPAC
ambiguity code (R, Y, N ...) counts as a mismatch here, i.e. the audit is strict.
"""
import collections
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from mfe_common import IUPAC, reverse_complement as rc  # noqa: E402

LUT = np.zeros(256, dtype=np.uint8)                 # A,C,G,T -> 1..4 ; everything else 0 (never matches)
for _i, _c in enumerate("ACGT"):
    LUT[ord(_c)] = _i + 1


def read_fasta(path):
    seqs, key = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line.startswith(">"):
            key = line[1:].split()[0]
            seqs[key] = []
        elif key:
            seqs[key].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def best_mismatch(oligo, codes, both_strands=True):
    """Fewest mismatches of `oligo` against any window of the genome (codes = LUT-encoded)."""
    best = 99
    for o in (oligo, rc(oligo)) if both_strands else (oligo,):
        n = len(codes) - len(o) + 1
        if n <= 0:
            continue
        mm = np.zeros(n, dtype=np.int16)
        for j, ch in enumerate(o):
            allowed = np.zeros(5, dtype=bool)
            for b in IUPAC[ch]:
                allowed["ACGT".index(b) + 1] = True
            mm += ~allowed[codes[j:j + n]]
        best = min(best, int(mm.min()))
    return best


def main():
    wanted = set(sys.argv[1:])
    seg_file = {}
    for p in (ROOT / "data" / "filtered").glob("*_filtered.fasta"):
        seg_file[p.name.split("_")[0].upper()] = p
    total_diff = 0
    for r in csv.DictReader(open(ROOT / "config" / "primers.tsv", encoding="utf-8"), delimiter="\t"):
        pid = r["primer_id"]
        if r.get("proceed", "True").lower() != "true" or (wanted and pid not in wanted):
            continue
        table = ROOT / "output" / "primers" / pid / "alignment" / "mismatch_table.tsv"
        fasta = seg_file.get(r["segment"].upper())
        if not table.exists() or fasta is None:
            continue
        genomes = read_fasta(fasta)
        rows = {x["genome_id"]: x for x in csv.DictReader(open(table, encoding="utf-8"), delimiter="\t")}
        diff, limited, n = collections.Counter(), collections.Counter(), 0
        for gid, g in genomes.items():
            x = rows.get(gid)
            if not x or "" in (x["f_mismatch"], x["r_mismatch"], x["p_mismatch"]):
                continue
            n += 1
            codes = LUT[np.frombuffer(g.encode(), dtype=np.uint8)]
            truth = {"F": best_mismatch(r["forward_primer"].upper(), codes),
                     "R": best_mismatch(r["reverse_primer"].upper(), codes),
                     "P": best_mismatch(r["probe"].upper(), codes)}
            for k, col in (("F", "f_mismatch"), ("R", "r_mismatch"), ("P", "p_mismatch")):
                if truth[k] > int(x[col]):
                    diff[k] += 1                      # the table claims fewer mismatches than the genome allows
                elif truth[k] < int(x[col]):
                    limited[k] += 1
        total_diff += sum(diff.values())
        print(f"{pid:6s} {r['primer_name'][:28]:28s} genomes={n:5d}  under-counted F={diff['F']} R={diff['R']} P={diff['P']}"
              f"  (pair-limited F={limited['F']} R={limited['R']} P={limited['P']})  {'OK' if not diff else 'ERROR'}")
    return 1 if total_diff else 0


if __name__ == "__main__":
    sys.exit(main())
