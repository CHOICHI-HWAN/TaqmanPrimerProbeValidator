#!/usr/bin/env python3
"""Coverage-first F / Probe / R design from a lineage's alignment.

varVAMP and PrimerHunter pick oligos from a consensus or a melting-temperature model, which is not the rule this pipeline validates
with: primers may differ in at most `primer_mismatch_pass` positions with the last bases at their 3' end intact, the probe must match
perfectly (IUPAC-aware). Here that rule is the objective: for every candidate window of the lineage alignment the genomes it would pass
are computed (a bit per genome), F, R and probe are combined so that the genomes passing ALL THREE are as many as possible, and further
designs are chosen for the genomes the earlier ones miss (greedy set cover). Every design is validated afterwards with the real
pipeline validation, exactly like the other engines' designs.

    python3 scripts/redesign_coverage.py output/redesign/s/human/LIII/aligned.fa [--top 3 --ambig 2 --probe-ambig 2]
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import probe_designer as pd  # noqa: E402
from mfe_common import IUPAC, find_mafft, reverse_complement  # noqa: E402

LEN_PRIMER = range(18, 26)
LEN_PROBE = range(20, 31)
POOL_BOUND = 4000          # windows per oligo type that get the exact (IUPAC) evaluation, best upper bound first
POOL_KEEP = 1200           # evaluated windows per type kept for the combination step


def read_aligned(path):
    ids, rows = [], []
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if ln.startswith(">"):
            ids.append(ln[1:].split()[0])
            rows.append([])
        elif ln and rows:
            rows[-1].append(ln.upper())
    return ids, ["".join(r) for r in rows]


def n_ambig(seq):
    return sum(c not in "ACGT" for c in seq)


def _bits(ok):
    return int.from_bytes(np.packbits(ok, bitorder="little").tobytes(), "little")


class Aligned:
    def __init__(self, ids, rows, std_rows=None):
        self.ids, self.G, self.L = ids, len(rows), len(rows[0])
        self.S = np.frombuffer("".join(std_rows).encode(), dtype=np.uint8).reshape(len(std_rows), self.L) if std_rows else None
        self.A = np.frombuffer("".join(rows).encode(), dtype=np.uint8).reshape(self.G, self.L)
        cons = []
        for c in range(self.L):
            col = self.A[:, c]
            counts = {b: int((col == ord(b)).sum()) for b in "ACGT"}
            top = max(counts, key=counts.get)
            cons.append(top if counts[top] >= 0.5 * self.G else "-")       # a column most genomes have a gap/N in is not designed over
        self.cons = "".join(cons)
        self.cons_u8 = np.frombuffer(self.cons.encode(), dtype=np.uint8)
        self.cs = np.concatenate([np.zeros((self.G, 1), dtype=np.int32), np.cumsum(self.A != self.cons_u8, axis=1, dtype=np.int32)], axis=1)
        self.gapcs = np.concatenate([[0], np.cumsum([c == "-" for c in self.cons])])

    def windows(self, k):
        """Start columns whose k-window is entirely conserved-column (no '-' in the consensus)."""
        n = self.L - k + 1
        return np.nonzero((self.gapcs[k:k + n] - self.gapcs[:n]) == 0)[0] if n > 0 else np.array([], dtype=int)

    def upper_bound(self, k, starts, thr):
        """Genomes with <= thr mismatches against the consensus window, for many windows at once."""
        mm = self.cs[:, starts + k] - self.cs[:, starts]
        return (mm <= thr).sum(axis=0)

    def evaluate(self, s, k, kind, amb, pass_mm, tp_n):
        """(oligo, bits) for the window [s, s+k): the consensus with up to `amb` greedy IUPAC codes, and the genomes it passes under the
        pipeline's rule (kind F / R / P). R's 3' end is the window's LEFT end (the oligo is the reverse complement of the window)."""
        window = self.cons[s:s + k]
        sub = self.A[:, s:s + k]
        diff = sub != self.cons_u8[s:s + k]
        row_mm = [{int(j): chr(sub[g, j]) for j in np.nonzero(diff[g])[0]} if diff[g].any() else {} for g in range(self.G)]
        oligo, _ = pd.degenerate_window(window, row_mm, amb)
        tab = np.zeros((k, 256), dtype=bool)
        for j, ch in enumerate(oligo):
            for b in IUPAC.get(ch, "ACGT"):
                tab[j, ord(b)] = True
        def passes(block):
            bad = ~tab[np.arange(k)[None, :], block]
            nbad = bad.sum(axis=1)
            if kind == "P":
                return nbad == 0
            end = slice(k - tp_n, k) if kind == "F" else slice(0, tp_n)
            return (nbad <= pass_mm) & ~bad[:, end].any(axis=1)
        sbits = _bits(passes(self.S[:, s:s + k])) if self.S is not None else 0
        return oligo, _bits(passes(sub)), sbits


def _stats(seq):
    res = pd.resolutions(seq, 64)
    if not res:
        return None
    return res, sum(pd.tm_nn(x) for x in res) / len(res), sum(pd.gc_percent(x) for x in res) / len(res)


def candidates(al, kind, amb, pass_mm, tp_n, primer_tm, gc, say):
    thr = (amb if kind == "P" else pass_mm + amb)
    pool = []
    for k in (LEN_PROBE if kind == "P" else LEN_PRIMER):
        starts = al.windows(k)
        if not len(starts):
            continue
        ub = al.upper_bound(k, starts, thr)
        pool.extend(zip(ub.tolist(), starts.tolist(), [k] * len(starts)))
    pool.sort(key=lambda x: (-x[0], x[1]))
    out = []
    for ub, s, k in pool[:POOL_BOUND]:
        oligo, bits, sbits = al.evaluate(s, k, kind, amb, pass_mm, tp_n)
        cnt = bits.bit_count()
        strands = [("+", oligo), ("-", reverse_complement(oligo))] if kind == "P" else [("+", oligo if kind == "F" else reverse_complement(oligo))]
        for strand, seq in strands:
            st = _stats(seq)
            if not st:
                continue
            res, tm, gcp = st
            if kind == "P":
                if any(x[0] == "G" for x in res) or any(pd.longest_run(x) >= 4 for x in res):
                    continue
            else:
                if not (primer_tm[0] <= tm <= primer_tm[1]) or not (gc[0] <= gcp <= gc[1]) or any(pd.longest_run(x) >= 5 for x in res):
                    continue
            out.append({"cnt": cnt, "s": s, "e": s + k, "seq": seq, "strand": strand, "tm": tm, "gc": gcp, "bits": bits, "sbits": sbits, "amb": n_ambig(seq)})
    out.sort(key=lambda x: (-x["cnt"], x["amb"], x["s"]))
    say(f"{kind} 후보 {len(out[:POOL_KEEP])}개 (최대 {out[0]['cnt'] if out else 0}개 유전체 통과)")
    return out[:POOL_KEEP]


def best_triples(F, R, P, covered, amp, probe_above, max_dtm, want=12, require_std=False):
    """Triples (F, R, P) ranked by the NEW genomes they pass (those not in `covered`), then by the genomes they pass in all. The lists are
    sorted by genomes passed (best first), so a branch whose best possible count cannot beat the current `want`-th best is cut off."""
    import heapq
    heap, seq = [], 0                                      # min-heap of (new, total, -amb, tiebreak, f, r, p)
    floor = 0
    for p in P:
        if p["cnt"] <= floor:
            break
        fs = [f for f in F if f["e"] <= p["s"] and p["e"] - f["s"] <= amp[1]]
        rs = [r for r in R if r["s"] >= p["e"] and r["e"] - p["s"] <= amp[1]]
        for f in fs:
            if f["cnt"] <= floor:
                break
            fp = f["bits"] & p["bits"] & ~covered
            if fp.bit_count() <= floor:
                continue
            for r in rs:
                if r["cnt"] <= floor:
                    break
                if not (amp[0] <= r["e"] - f["s"] <= amp[1]) or abs(f["tm"] - r["tm"]) > max_dtm:
                    continue
                hi = max(f["tm"], r["tm"])
                if not (hi + probe_above[0] <= p["tm"] <= hi + probe_above[1]):
                    continue
                new = (fp & r["bits"]).bit_count()
                if new <= floor:
                    continue
                if require_std and not (f["sbits"] & r["sbits"] & p["sbits"]):      # at least one standard material passes the whole set
                    continue
                tot = (f["bits"] & p["bits"] & r["bits"]).bit_count()
                seq += 1
                item = (new, tot, -(f["amb"] + r["amb"] + p["amb"]), -abs(f["tm"] - r["tm"]), -seq, f, r, p)
                if len(heap) < want:
                    heapq.heappush(heap, item)
                else:
                    heapq.heappushpop(heap, item)
                if len(heap) >= want:
                    floor = heap[0][0] - 1 if heap[0][0] > 0 else 0       # equal counts may still win on the tie-breakers
    out = sorted(heap, reverse=True)
    return [(x[0], x[1], x[5], x[6], x[7]) for x in out]


def design_from_alignment(ids, rows, ambig=2, probe_ambig=2, top=3, pass_mm=2, tp_n=2, amp=(70, 200), primer_tm=(52.0, 66.0),
                          gc=(30.0, 75.0), probe_above=(2.0, 14.0), max_dtm=5.0, progress=None, std_rows=None, require_std=False):
    """-> {"status": "ok"|"none", "n", "designs": [...], "message"}; each design: name, F, P, R, F_tm .. R_gc, amplicon_len,
    est_cov (genomes of the alignment it passes under the pipeline's rule), est_added (new ones beyond the earlier designs)."""
    def say(msg):
        if progress:
            progress(msg)
    al = Aligned(ids, rows, std_rows)
    require_std = bool(require_std and std_rows)
    t0 = time.time()
    F = candidates(al, "F", ambig, pass_mm, tp_n, primer_tm, gc, say)
    R = candidates(al, "R", ambig, pass_mm, tp_n, primer_tm, gc, say)
    P = candidates(al, "P", probe_ambig, pass_mm, tp_n, primer_tm, gc, say)
    if not (F and R and P):
        return {"status": "none", "n": al.G, "designs": [], "message": "조건(Tm·GC·모호 염기)을 만족하는 F/R/Probe 후보가 부족해 설계를 내지 못했습니다."}
    designs, covered = [], 0
    for i in range(top):
        say(f"조합 탐색 {i + 1}/{top}...")
        cand = best_triples(F, R, P, covered, amp, probe_above, max_dtm, require_std=require_std)
        # the next design must differ from the earlier ones in at least one oligo (no repeats of the same amplicon)
        cand = [c for c in cand if not any(abs(c[2]["s"] - d["_f"]) < 6 and abs(c[3]["e"] - d["_r"]) < 6 for d in designs)]
        if not cand:
            break
        new, total, f, r, p = cand[0]
        covered |= f["bits"] & r["bits"] & p["bits"]
        designs.append({"name": f"COV_{i}", "F": f["seq"], "P": p["seq"], "R": r["seq"],
                        "F_tm": round(f["tm"], 1), "P_tm": round(p["tm"], 1), "R_tm": round(r["tm"], 1),
                        "F_gc": round(f["gc"], 1), "P_gc": round(p["gc"], 1), "R_gc": round(r["gc"], 1),
                        "amplicon_len": int(sum(1 for c in al.cons[f["s"]:r["e"]] if c != "-")),
                        "std_est": (bool(f["sbits"] & r["sbits"] & p["sbits"]) if std_rows else None),
                        "est_cov": total, "est_cov_pct": round(100.0 * total / al.G, 1), "est_added": new, "P_strand": p["strand"],
                        "_f": f["s"], "_r": r["e"]})
    for d in designs:
        d.pop("_f"), d.pop("_r")
    if not designs:
        why = "표준물질이 통과하면서 " if require_std else ""
        return {"status": "none", "n": al.G, "designs": [], "message": f"{why}F·R·Probe가 함께 맞는 조합을 찾지 못했습니다(증폭산물 길이·Tm 조건을 확인하세요)."}
    say(f"설계 {len(designs)}개 ({time.time() - t0:.0f}초)")
    return {"status": "ok", "n": al.G, "designs": designs, "message": "", "threshold": None}


def design_lineage(members, out_dir, ambig=2, probe_ambig=2, top=3, threads=4, progress=None, standards=None, require_std=False, **kw):
    """members = {accession: sequence}: MAFFT alignment, then design_from_alignment(). `standards` = {id: sequence} of the reference
    standard material: aligned onto the same columns (mafft --keeplength --add) and judged with the same rule; `require_std` demands that
    at least one of them passes the designed set."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n = len(members)
    if n < 3:
        return {"status": "error", "n": n, "designs": [], "message": f"게놈이 {n}개뿐이라 정렬·설계를 할 수 없습니다 (3개 이상 필요)."}
    fa, aln = out_dir / "members.fa", out_dir / "aligned.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in members.items()), encoding="utf-8")
    if progress:
        progress(f"MAFFT 정렬 중 ({n}개)...")
    with open(aln, "w", encoding="utf-8") as fh:
        subprocess.run([find_mafft(), "--auto", "--thread", str(threads), str(fa)], stdout=fh, stderr=subprocess.DEVNULL, check=True)
    ids, rows = read_aligned(aln)
    std_rows = None
    if standards:
        sfa, saln = out_dir / "standards.fa", out_dir / "aligned_with_standards.fa"
        sfa.write_text("".join(f">{k}\n{v}\n" for k, v in standards.items()), encoding="utf-8")
        with open(saln, "w", encoding="utf-8") as fh:
            subprocess.run([find_mafft(), "--keeplength", "--thread", str(threads), "--add", str(sfa), str(aln)], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        all_ids, all_rows = read_aligned(saln)
        by = dict(zip(all_ids, all_rows))
        std_rows = [by[k] for k in standards if k in by]
    return design_from_alignment(ids, rows, ambig=ambig, probe_ambig=probe_ambig, top=top, progress=progress, std_rows=std_rows,
                                 require_std=require_std, **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("aligned_fasta")
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--ambig", type=int, default=2)
    ap.add_argument("--probe-ambig", type=int, default=2)
    a = ap.parse_args()
    ids, rows = read_aligned(a.aligned_fasta)
    res = design_from_alignment(ids, rows, ambig=a.ambig, probe_ambig=a.probe_ambig, top=a.top, progress=print)
    for d in res["designs"]:
        print(f"{d['name']}  {d['est_cov']}/{res['n']} ({d['est_cov_pct']}%) +{d['est_added']}  F {d['F']}  P {d['P']}  R {d['R']}  {d['amplicon_len']}bp")
    if res["status"] != "ok":
        print(res["message"])


if __name__ == "__main__":
    main()
