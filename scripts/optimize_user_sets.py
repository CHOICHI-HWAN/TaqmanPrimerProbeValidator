#!/usr/bin/env python3
"""Can a user's primer/probe set cover more genomes with a small change? Tries, for every non-auto-designed set of a segment:
probe only (any probe inside the amplicon), primers only (F and/or R shifted a few bases, IUPAC codes added), and both together.

    python3 scripts/optimize_user_sets.py s --mm 2            # needs output/design_sweep/<seg>/aligned.fa (design_sweep.py makes it)

Rules for a variant: primer mismatch <= mm with the 3' end intact, probe perfect (IUPAC-aware), probe inside the amplicon, amplicon <= --amp-max,
primer Tm within --dtm of the original, probe Tm at least 2 C above both primers, at most --max-amb IUPAC codes added to a primer.
Result: output/design_sweep/<seg>/OPTIMIZE_<seg>.tsv (contains sequences: it stays in the git-ignored output folder).
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import redesign_coverage as rc  # noqa: E402
from mfe_common import IUPAC, reverse_complement  # noqa: E402


def cfg_int(key, default):
    m = re.search(rf"^{key}:\s*(\d+)", (ROOT / "config" / "config.yaml").read_text(encoding="utf-8"), re.M)
    return int(m.group(1)) if m else default


def tab_of(oligo):
    t = np.zeros((len(oligo), 256), dtype=bool)
    for j, ch in enumerate(oligo):
        for b in IUPAC.get(ch, "ACGT"):
            t[j, ord(b)] = True
    return t


def locate(al, oligo):
    """(start, genomes) of the window of the alignment that most genomes match with <= 3 mismatches (IUPAC-aware), None if no genome does.
    Works on the genomes themselves, not on the consensus, so oligos of diverse viruses are found too."""
    k = len(oligo)
    t = tab_of(oligo)
    n = al.L - k + 1
    if n <= 0:
        return None
    best = (0, -1)
    for c0 in range(0, n, 20000):
        c1 = min(n, c0 + 20000)
        mism = np.zeros((al.G, c1 - c0), dtype=np.int8)
        for j in range(k):
            mism += ~t[j, al.A[:, c0 + j:c1 + j]]
        cnt = (mism <= 3).sum(axis=0)
        i = int(cnt.argmax())
        if cnt[i] > best[1]:
            best = (c0 + i, int(cnt[i]))
    return best if best[1] > 0 else None


def pass_bits(al, oligo, s, kind, mm, tp_n):
    k = len(oligo)
    bd = ~tab_of(oligo)[np.arange(k)[None, :], al.A[:, s:s + k]]
    nb = bd.sum(axis=1)
    if kind == "P":
        return rc._bits(nb == 0)
    end = bd[:, k - tp_n:] if kind == "F" else bd[:, :tp_n]
    return rc._bits((nb <= mm) & ~end.any(axis=1))


def tm_of(seq):
    st = rc._stats(seq)
    return st[1] if st else 0.0


def mk(al, kind, s, k, amb, mm, tp_n):
    oligo, bits, _ = al.evaluate(s, k, kind, amb, mm, tp_n)
    seq = oligo if kind == "F" else reverse_complement(oligo)
    return {"s": s, "e": s + k, "seq": seq, "bits": bits, "cnt": bits.bit_count(), "tm": tm_of(seq), "amb": rc.n_ambig(seq), "kind": kind}


def primer_variants(al, orig, kind, shift, max_amb, mm, tp_n, dtm, fixed_end):
    """Windows near the original primer; the original itself is always in the list."""
    out = {(orig["s"], orig["e"], orig["seq"]): dict(orig, changed=False)}
    for k in rc.LEN_PRIMER:
        for s in range(orig["s"] - shift, orig["s"] + shift + 1):
            if s < 0 or s + k > al.L or not (al.gapcs[s + k] - al.gapcs[s] == 0):
                continue
            if (kind == "F" and s + k > fixed_end) or (kind == "R" and s < fixed_end):
                continue
            for amb in range(0, max_amb + 1):
                v = mk(al, kind, s, k, amb, mm, tp_n)
                if abs(v["tm"] - orig["tm"]) > dtm or v["amb"] > max_amb + orig["amb"]:
                    continue
                st = rc._stats(v["seq"])
                if not st or any(pd_run(x) >= 5 for x in st[0]):
                    continue
                out[(v["s"], v["e"], v["seq"])] = dict(v, changed=True)
    top = sorted(out.values(), key=lambda x: (-x["cnt"], x["amb"]))[:8]
    return top if any(not x["changed"] for x in top) else top + [out[(orig["s"], orig["e"], orig["seq"])]]


def pd_run(x):
    return rc.pd.longest_run(x)


def probe_variants(al, lo, hi, mm, tp_n, orig):
    pool = {}
    for k in rc.LEN_PROBE:
        starts = [int(x) for x in al.windows(k) if lo <= x and x + k <= hi]
        if not starts:
            continue
        starts = np.array(starts)
        ub = al.upper_bound(k, starts, 3)
        for i in np.argsort(-ub)[:60]:
            for amb in range(0, 4):
                v = mk(al, "P", int(starts[i]), k, amb, mm, tp_n)
                for seq in (v["seq"], reverse_complement(v["seq"])):
                    st = rc._stats(seq)
                    if not st or any(x[0] == "G" for x in st[0]) or any(pd_run(x) >= 4 for x in st[0]):
                        continue
                    pool[(v["s"], v["e"], seq)] = dict(v, seq=seq, tm=st[1], changed=True)
    if orig:
        pool[(orig["s"], orig["e"], orig["seq"])] = dict(orig, changed=False)
    top = sorted(pool.values(), key=lambda x: (-x["cnt"], x["amb"]))[:60]
    return top if not orig or any(not x["changed"] for x in top) else top + [pool[(orig["s"], orig["e"], orig["seq"])]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("segment")
    ap.add_argument("--mm", type=int, default=cfg_int("primer_mismatch_pass", 1))
    ap.add_argument("--shift", type=int, default=8)
    ap.add_argument("--max-amb", type=int, default=2)
    ap.add_argument("--amp-max", type=int, default=250)
    ap.add_argument("--dtm", type=float, default=3.0)
    a = ap.parse_args()
    tp_n = cfg_int("three_prime_consecutive_match", 2)
    d = ROOT / "output" / "design_sweep" / a.segment
    ids, rows = rc.read_aligned(d / "aligned.fa")
    al = rc.Aligned(ids, rows)
    manifest = list(csv.DictReader((ROOT / "config" / "primers.tsv").open(encoding="utf-8"), delimiter="\t"))
    sets = [r for r in manifest if r["segment"].lower() == a.segment.lower() and r["source"] != "AUTO-DESIGN"]
    out = []
    print(f"{a.segment}: {al.G} genomes, {len(sets)} user sets, mm {a.mm}", flush=True)
    for r in sets:
        fw, rv, pr = r["forward_primer"].upper(), r["reverse_primer"].upper(), r["probe"].upper()
        # the F/R columns of a manifest are sometimes swapped: find each primer on either strand, the one in genome sense is the forward primer
        loc = {}
        for name, seq in (("a", fw), ("b", rv)):
            plus, minus = locate(al, seq), locate(al, reverse_complement(seq))
            if plus and (not minus or plus[1] >= minus[1]):
                loc[name] = ("+", plus)
            elif minus:
                loc[name] = ("-", minus)
        if len(loc) < 2 or loc["a"][0] == loc["b"][0]:
            loc = {}
        if loc and loc["a"][0] == "-":
            fw, rv = rv, fw
            loc = {"a": loc["b"], "b": loc["a"]}
        lf, lr = (loc["a"][1], loc["b"][1]) if loc else (None, None)
        lp = [(locate(al, pr), pr), (locate(al, reverse_complement(pr)), reverse_complement(pr))]
        lp = [(x, p) for x, p in lp if x]
        if not (lf and lr and lp):
            print(f"  {r['primer_id']} {r['primer_name']}: not located in the alignment (F {bool(lf)} R {bool(lr)} P {bool(lp)})", flush=True)
            continue
        (ps, _), pw = max(lp, key=lambda x: x[0][1])
        F = {"s": lf[0], "e": lf[0] + len(fw), "seq": fw, "bits": pass_bits(al, fw, lf[0], "F", a.mm, tp_n), "tm": tm_of(fw), "amb": rc.n_ambig(fw), "kind": "F", "cnt": 0}
        R = {"s": lr[0], "e": lr[0] + len(rv), "seq": rv, "bits": pass_bits(al, reverse_complement(rv), lr[0], "R", a.mm, tp_n), "tm": tm_of(rv), "amb": rc.n_ambig(rv), "kind": "R", "cnt": 0}
        P = {"s": ps, "e": ps + len(pw), "seq": pr if pw == pr else pr, "bits": pass_bits(al, pw, ps, "P", a.mm, tp_n), "tm": tm_of(pr), "amb": rc.n_ambig(pr), "kind": "P", "cnt": 0}
        for x in (F, R, P):
            x["cnt"] = x["bits"].bit_count()
        base = (F["bits"] & R["bits"] & P["bits"]).bit_count()
        Fs = primer_variants(al, F, "F", a.shift, a.max_amb, a.mm, tp_n, a.dtm, P["s"])
        Rs = primer_variants(al, R, "R", a.shift, a.max_amb, a.mm, tp_n, a.dtm, P["e"])
        Ps = probe_variants(al, F["s"] - a.shift, R["e"] + a.shift, a.mm, tp_n, P)
        best = {}                                                    # change class -> best triple
        for f in Fs:
            for rr in Rs:
                for p in Ps:
                    if p["s"] < f["e"] or p["e"] > rr["s"] or rr["e"] - f["s"] > a.amp_max:
                        continue
                    if p["tm"] < max(f["tm"], rr["tm"]) + 2.0 and (f["changed"] or rr["changed"] or p["changed"]):
                        continue
                    cnt = (f["bits"] & rr["bits"] & p["bits"]).bit_count()
                    cls = ("F" if f["changed"] else "") + ("R" if rr["changed"] else "") + ("P" if p["changed"] else "")
                    key = (cnt, -(f["amb"] + rr["amb"] + p["amb"]))
                    if cls and (cls not in best or key > best[cls][0]):
                        best[cls] = (key, f, rr, p)
        print(f"  {r['primer_id']} {r['primer_name']} [{r['source'][:14]}]: now {base}/{al.G} = {100 * base / al.G:.1f}%", flush=True)
        for cls in ("P", "F", "R", "FR", "FP", "RP", "FRP"):
            if cls in best and best[cls][0][0] > base:
                (cnt, _), f, rr, p = best[cls]
                print(f"      change {cls:3s}: {cnt}/{al.G} = {100 * cnt / al.G:.1f}%  (+{cnt - base})", flush=True)
                out.append({"set_id": r["primer_id"], "name": r["primer_name"], "change": cls, "now_pct": round(100 * base / al.G, 2), "new_pct": round(100 * cnt / al.G, 2),
                            "F": f["seq"], "P": p["seq"], "R": rr["seq"], "F_tm": round(f["tm"], 1), "P_tm": round(p["tm"], 1), "R_tm": round(rr["tm"], 1), "amplicon": rr["e"] - f["s"]})
    with open(d / f"OPTIMIZE_{a.segment}.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["set_id", "name", "change", "now_pct", "new_pct", "F", "P", "R", "F_tm", "P_tm", "R_tm", "amplicon"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} improvements -> {d / ('OPTIMIZE_' + a.segment + '.tsv')}")


if __name__ == "__main__":
    main()
