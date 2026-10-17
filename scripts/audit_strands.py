"""Strand / reverse-complement audit: for every genuine primer pair (F, R <= 2 mismatches, probe <= 8), F, R and the
probe are re-counted against the STORED amplicon sequence on both strands and must reproduce the table's mismatch
values; F and R must sit at the amplicon ends on opposite strands. Run from the project folder."""
import csv, sys, collections
from pathlib import Path
root = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(root / "scripts"))
from mfe_common import _oligo_sites
seen = set(); T = collections.Counter()
for r in csv.DictReader(open(root / "config/primers.tsv"), delimiter="\t"):
    if r.get("proceed", "True").lower() != "true": continue
    pid = r["primer_id"]; key = (r["forward_primer"], r["reverse_primer"], r["probe"], r["segment"])
    if key in seen: continue
    seen.add(key)
    F, R, P = r["forward_primer"].upper(), r["reverse_primer"].upper(), r["probe"].upper()
    bad = collections.Counter(); n = 0; strands = collections.Counter()
    for x in csv.DictReader(open(root / f"output/primers/{pid}/alignment/mismatch_table.tsv"), delimiter="\t"):
        if not x["amplicon_seq"] or x["f_mismatch"] == "" or x["p_mismatch"] == "": continue
        if int(x["f_mismatch"]) > 2 or int(x["r_mismatch"]) > 2 or int(x["p_mismatch"]) > 8: continue      # real pairs only
        amp = x["amplicon_seq"].upper(); n += 1
        # F, R and probe must each bind SOMEWHERE inside the stored amplicon with exactly the reported mismatches
        for name, o, col in (("F", F, "f_mismatch"), ("R", R, "r_mismatch"), ("P", P, "p_mismatch")):
            sites = _oligo_sites(o, amp, 8)
            best = min((s[3] for s in sites), default=99)
            if best != int(x[col]): bad[name] += 1
        # the primers must face each other: F and R sites on opposite strands, ends of the amplicon
        fs = [s for s in _oligo_sites(F, amp, 2) if s[0] == 0 or s[1] == len(amp)]
        rs = [s for s in _oligo_sites(R, amp, 2) if s[0] == 0 or s[1] == len(amp)]
        ok = any(f[2] != r_[2] for f in fs for r_ in rs)
        if not ok: bad["geometry"] += 1
        ps = sorted(_oligo_sites(P, amp, 8), key=lambda s: s[3])
        if ps: strands[ps[0][2]] += 1
    T["sets"] += 1; T["rows"] += n; T["bad"] += sum(bad.values())
    print(f"{pid:5s} amplicons checked={n:4d}  disagreements F/R/P/geometry = {bad['F']}/{bad['R']}/{bad['P']}/{bad['geometry']}   probe strand (+/-) = {strands['+']}/{strands['-']}")
print("TOTAL", dict(T))
