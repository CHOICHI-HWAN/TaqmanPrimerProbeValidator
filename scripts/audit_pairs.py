"""Audit of the F/R pair search: for every genome, is there a valid amplicon (each primer <= 1 mismatch, intact 3' end,
<= 2000 bp, facing each other) anywhere in the genome by brute force -- and does the pipeline agree?

Catches sites mfeprimer's 3' k-mer seed misses. Run from the project folder: python3 scripts/audit_pairs.py
"""
import csv, sys, glob, collections
from pathlib import Path
import numpy as np
root = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(root / "scripts"))
import audit_recount as A
from mfe_common import IUPAC, reverse_complement as rc
MAXSIZE = 2000
def mm_array(o, codes):
    n = len(codes) - len(o) + 1
    if n <= 0: return np.zeros(0, dtype=np.int16)
    mm = np.zeros(n, dtype=np.int16)
    for j, ch in enumerate(o):
        al = np.zeros(5, dtype=bool)
        for b in IUPAC[ch]: al["ACGT".index(b) + 1] = True
        mm += ~al[codes[j:j + n]]
    return mm
def sites(o, codes, thr, n3):
    """(start, end_exclusive, strand, mm) of windows with mm<=thr and a perfect 3' end (last n3 bases of the oligo)"""
    out = []
    L = len(o)
    for strand, q in (("+", o), ("-", rc(o))):
        mm = mm_array(q, codes)
        idx = np.nonzero(mm <= thr)[0]
        for i in idx:
            # 3' end of the oligo: last n3 chars of q on '+', first n3 chars on '-'
            seg = q[-n3:] if strand == "+" else q[:n3]
            tgt = codes[i + L - n3:i + L] if strand == "+" else codes[i:i + n3]
            ok = all(("ACGT"[t - 1] if t else "?") in IUPAC[c] for c, t in zip(seg, tgt))
            if ok: out.append((int(i), int(i) + L, strand, int(mm[i])))
    return out
def pair_ok(fs, rs):
    for f in fs:
        for r in rs:
            if f[2] == "+" and r[2] == "-" and r[1] > f[0] and r[1] - f[0] <= MAXSIZE and r[0] >= f[0]: return True
            if r[2] == "+" and f[2] == "-" and f[1] > r[0] and f[1] - r[0] <= MAXSIZE and f[0] >= r[0]: return True
    return False
seen = set(); tot = collections.Counter()
for r in csv.DictReader(open(root / "config/primers.tsv"), delimiter="\t"):
    if r.get("proceed", "True").lower() != "true": continue
    pid, seg = r["primer_id"], r["segment"].lower()
    key = (r["forward_primer"], r["reverse_primer"], seg)
    if key in seen: continue
    seen.add(key)
    genomes = A.read_fasta(root / f"data/filtered/{seg}_filtered.fasta")
    tab = {x["genome_id"]: x for x in csv.DictReader(open(root / f"output/primers/{pid}/alignment/mismatch_table.tsv"), delimiter="\t")}
    spec = glob.glob(str(root / f"output/primers/{pid}/primer_qc/*mm2*.spec.tsv"))
    found2 = set()
    if spec:
        for line in open(spec[0]):
            c = line.split("\t")
            if not line.startswith("#") and len(c) > 1: found2.add(c[1])
    F, R = r["forward_primer"].upper(), r["reverse_primer"].upper()
    truth = pipe = both = missed = 0; nospec_but_true = 0; ex = []
    for gid, g in genomes.items():
        codes = A.LUT[np.frombuffer(g.encode(), dtype=np.uint8)]
        t = pair_ok(sites(F, codes, 1, 2), sites(R, codes, 1, 2))
        x = tab.get(gid)
        p = bool(x and x["f_mismatch"] != "" and int(x["f_mismatch"]) <= 1 and int(x["r_mismatch"]) <= 1 and x["f_3p_ok"] in ("", "1") and x["r_3p_ok"] in ("", "1"))
        truth += t; pipe += p; both += (t and p)
        if t and not p:
            missed += 1; nospec_but_true += (gid not in found2)
            if len(ex) < 2: ex.append(gid)
    print(f"{pid:5s} {seg} brute-force valid pair (<=1 mm each, 3' ok, <=2000bp): {truth:4d} | pipeline: {pipe:4d} | missed by pipeline: {missed:3d} (of which mfeprimer returned NO amplicon at all: {nospec_but_true}) | pipeline-only: {pipe-both}", ex)
