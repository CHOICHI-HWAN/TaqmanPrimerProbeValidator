#!/usr/bin/env python3
"""Automatic TaqMan probe design + test for one assay set (F and R stay as registered).

1. Take the set's amplicon alignment (output/primers/<id>/msa/msa.json): the consensus between the F and R windows.
2. Slide 18-30 nt windows over it, on both strands, and keep the ones that pass the TaqMan rules
   (length, GC, Tm, no 5' G, no long homopolymer / dinucleotide runs, no overlap with F/R).
3. Score every candidate on the alignment itself (genomes matching with 0 / 1 mismatches), keep the best ones.
4. Re-test those with the pipeline's own MFEprimer probe search against the whole filtered population
   (config probe_mismatch, default 0) -- the same test a registered probe goes through -- and write a ranked table.

    python3 scripts/probe_designer.py P025 [--top 10] [--min-len 18 --max-len 30] [--no-mfe]

Output: output/probe_design/<set_id>/candidates.tsv (+ candidates.fa, mfeprimer reports).
Candidate sequences are primer-grade data: that folder is git-ignored like the rest of output/.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from mfe_common import (  # noqa: E402
    IUPAC, ensure_index, find_mfeprimer, read_fasta_ids, reverse_complement, run_probe_search, write_oligo_fasta,
)

# SantaLucia (1998) unified nearest-neighbour parameters: dH kcal/mol, dS cal/mol/K
_NN = {"AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3), "CA": (-8.5, -22.7),
       "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4), "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
       "GA": (-8.2, -22.2), "TC": (-8.2, -22.2), "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9),
       "CC": (-8.0, -19.9)}


def tm_nn(seq: str, oligo_nm: float = 250.0, na_mm: float = 50.0) -> float:
    """Nearest-neighbour Tm (deg C) of a plain ACGT oligo, 250 nM oligo / 50 mM Na+ (non-self-complementary)."""
    seq = seq.upper()
    dh = sum(_NN[seq[i:i + 2]][0] for i in range(len(seq) - 1))
    ds = sum(_NN[seq[i:i + 2]][1] for i in range(len(seq) - 1))
    for end in (seq[0], seq[-1]):                       # terminal initiation
        dh, ds = (dh + 2.3, ds + 4.1) if end in "GC" else (dh + 0.1, ds - 2.8)
    ds += 0.368 * (len(seq) - 1) * math.log(na_mm / 1000.0)
    return dh * 1000.0 / (ds + 1.987 * math.log(oligo_nm * 1e-9 / 4.0)) - 273.15


PROBE_TM_ABOVE = (3.0, 10.0)       # probe Tm window relative to the HIGHER primer Tm (deg C)
PROBE_TM_TARGET = 6.0              # among equally good candidates prefer this much above the higher primer


def gc_percent(seq: str) -> float:
    return 100.0 * sum(c in "GC" for c in seq) / len(seq)


def longest_run(seq: str) -> int:
    best = cur = 1
    for a, b in zip(seq, seq[1:]):
        cur = cur + 1 if a == b else 1
        best = max(best, cur)
    return best


def resolutions(seq: str, limit: int = 0) -> list:
    """Every plain-ACGT oligo an IUPAC oligo stands for (2 per two-base code, so 2**k for k degenerate positions).
    With `limit`, an oligo that would stand for more than `limit` of them gives [] without expanding them (a window full of N
    would otherwise build 4**n strings)."""
    if limit:
        n = 1
        for c in seq.upper():
            n *= max(1, len(IUPAC.get(c, "")))
            if n > limit:
                return []
    outs = [""]
    for c in seq.upper():
        bases = IUPAC.get(c, "")
        outs = [o + b for o in outs for b in bases] if bases else []
    return outs


_CODE_OF = {frozenset(v): k for k, v in IUPAC.items() if len(v) == 2}     # {A,G} -> R, {C,T} -> Y, ...


def rule_problems(seq: str, tm_range, gc_range) -> list:
    """Why an oligo is not a usable TaqMan probe (empty list = usable). IUPAC positions are judged on every
    resolution: Tm and GC on their mean, the 5' G and run rules on the worst one."""
    res = resolutions(seq, 64)
    if not res or len(res) > 64:
        return ["non-ACGT"]
    why = []
    tm, gc = sum(tm_nn(r) for r in res) / len(res), sum(gc_percent(r) for r in res) / len(res)
    if not tm_range[0] <= tm <= tm_range[1]:
        why.append(f"Tm {tm:.1f}")
    if not gc_range[0] <= gc <= gc_range[1]:
        why.append(f"GC {gc:.0f}%")
    if any(r[0] == "G" for r in res):
        why.append("5' G (quenches the reporter)")
    if any(longest_run(r) >= 4 for r in res):
        why.append("4+ identical bases in a row")
    return why


def degenerate_window(window: str, row_mm: list, max_iupac: int):
    """Greedily turn up to `max_iupac` positions of the consensus window into two-base IUPAC codes, each time at the
    (position, base) that the most genomes still need. row_mm[i] = {pos: base} of genome i's mismatches against the window.
    -> (oligo, [uncovered mismatch count per genome])"""
    seq = list(window)
    allowed = [{c} for c in window]
    left = [dict(m) for m in row_mm]
    for k in range(max_iupac, 0, -1):
        score = {}
        for m in left:
            if 0 < len(m) <= k:                           # a genome that could still be fully covered with the codes left
                for pos, base in m.items():
                    if base in "ACGT" and base != window[pos] and len(allowed[pos]) == 1 and window[pos] in "ACGT":
                        score[(pos, base)] = score.get((pos, base), 0) + 1
        if not score:
            break
        (pos, base), _ = max(score.items(), key=lambda kv: (kv[1], -kv[0][0]))
        allowed[pos].add(base)
        seq[pos] = _CODE_OF[frozenset(allowed[pos])]
        for m in left:
            if m.get(pos) == base:
                del m[pos]
    return "".join(seq), [len(m) for m in left]


def load_msa(set_id: str) -> dict:
    path = ROOT / "output" / "primers" / set_id / "msa" / "msa.json"
    if not path.exists():
        sys.exit(f"[error] {path} not found -- run the set's analysis first.")
    return json.loads(path.read_text(encoding="utf-8"))


def design(set_id: str, min_len=18, max_len=30, tm_range=None, gc_range=(30.0, 80.0), top=10, max_iupac=2, tm_above=None, method="both", ids=None, msa=None):
    msa = msa or load_msa(set_id)
    cons = msa["consensus"].upper()
    f_col, r_col = msa.get("f_col"), msa.get("r_col")
    if not f_col or not r_col:
        sys.exit("[error] the alignment has no F/R window for this set.")
    # A probe should melt above BOTH primers (so it is bound before they extend): the reference is the higher of the two.
    # Their Tm is taken from the consensus under their windows (plain ACGT, so IUPAC primers work) with the same model
    # as the candidates.
    primer_windows = [w for w in (cons[c[0]:c[1] + 1].replace("-", "") for c in (f_col, r_col)) if w]   # a consensus gap inside a primer window is skipped
    if not primer_windows:
        sys.exit("[error] the consensus has no bases under the F/R primer windows.")
    primer_tm = max(tm_nn(w) for w in primer_windows)
    above = tm_above or PROBE_TM_ABOVE
    if tm_range is None:
        tm_range = (primer_tm + above[0], primer_tm + above[1])
    lo = min(f_col[1], r_col[1]) + 1                      # first column after the left primer
    hi = max(f_col[0], r_col[0]) - 1                      # last column before the right primer
    rows = msa["rows"]
    if ids is not None:                                    # host scope: only the genomes of the active scope count
        rows = [r for r in rows if r["genome_id"] in ids]
    n_rows = len(rows)
    out, rejected = [], 0
    use_own, use_p3 = method in ("own", "both"), method in ("primer3", "both")
    for length in (range(min_len, max_len + 1) if use_own else ()):
        for s in range(lo, hi - length + 2):
            e = s + length                                  # python slice end
            window = cons[s:e]
            # every genome's mismatches against the consensus window ('.' = same, a letter / '-' = different)
            row_mm = [{i - s: ch for i, ch in enumerate(r["dots"][s:e], s) if ch != "."} for r in rows]
            window, left = degenerate_window(window, row_mm, max_iupac)
            n_iupac = sum(c not in "ACGT" for c in window)
            cov0 = sum(m == 0 for m in left)
            cov1 = sum(m <= 1 for m in left)
            for strand in ("+", "-"):
                oligo = window if strand == "+" else reverse_complement(window)
                why = rule_problems(oligo, tm_range, gc_range)
                if why:
                    rejected += 1
                    continue
                res = resolutions(oligo)
                out.append({
                    "seq": oligo, "strand": strand, "col_start": s + 1, "col_end": e, "length": length,
                    "n_iupac": n_iupac, "tm": round(sum(tm_nn(r) for r in res) / len(res), 1),
                    "gc": round(sum(gc_percent(r) for r in res) / len(res), 1),
                    "msa_cov0": cov0, "msa_cov1": cov1, "msa_n": n_rows, "source": "own",
                })
    quota_own = top if method == "own" else (top - top // 2 if method == "both" else 0)
    quota_p3 = top if method == "primer3" else (top // 2 if method == "both" else 0)
    picked = _pick(out, quota_own, primer_tm)
    p3_info = {}
    if use_p3:
        p3_out, p3_info = design_primer3(cons, rows, lo, hi, f_col, r_col, min_len, max_len, gc_range, above, max_iupac,
                                         quota_p3, n_rows)
        picked += _pick(p3_out, quota_p3, primer_tm)
    merged = {}
    for c in picked:                                        # the same oligo found by both methods is listed once
        if c["seq"] in merged:
            merged[c["seq"]]["source"] = "own+primer3"
            for k in ("p3_penalty", "p3_hairpin"):
                if k in c:
                    merged[c["seq"]][k] = c[k]
        else:
            merged[c["seq"]] = c
    picked = sorted(merged.values(), key=lambda c: (-c["msa_cov0"], -c["msa_cov1"], c["n_iupac"], c["length"]))
    for c in picked:
        c["tm_p3"] = tm_primer3(c["seq"])
    return {"primer_tm": primer_tm, "tm_range": tm_range, "window": (lo + 1, hi + 1), "n_total": len(out), "n_rejected": rejected,
            "candidates": picked, "n_rows": n_rows, "existing": msa.get("ref_rows", {}).get("probe", {}), "method": method,
            "primer3": p3_info}


def _pick(cands, quota, primer_tm):
    """Best first (most genomes matching perfectly, then with 1 mismatch, fewer IUPAC, Tm near the target, shorter), without
    near-duplicates (same strand, start within 3 columns of a better candidate)."""
    cands = sorted(cands, key=lambda c: (-c["msa_cov0"], -c["msa_cov1"], c["n_iupac"],
                                         abs(c["tm"] - (primer_tm + PROBE_TM_TARGET)), c["length"]))
    picked = []
    for c in cands:
        if len(picked) >= quota:
            break
        if all(not (p["strand"] == c["strand"] and abs(p["col_start"] - c["col_start"]) <= 3) for p in picked):
            picked.append(c)
    return picked


def tm_primer3(seq: str):
    """Primer3's Tm of an oligo (mean over IUPAC resolutions), or None when primer3-py is not installed."""
    try:
        import primer3  # noqa: PLC0415
    except ImportError:
        return None
    res = resolutions(seq)
    if not res:
        return None
    return round(sum(primer3.calc_tm(r, mv_conc=50, dv_conc=0, dntp_conc=0, dna_conc=250) for r in res) / len(res), 1)


def design_primer3(cons, rows, lo, hi, f_col, r_col, min_len, max_len, gc_range, above, max_iupac, quota, n_rows):
    """Candidates picked by Primer3 itself (internal-oligo design on the amplicon consensus between F and R: its own Tm,
    GC, poly-X, hairpin and self-complementarity filters and penalty). Primer3 only knows plain ACGT, so the IUPAC codes are
    added afterwards exactly as for the own method, and both strands are offered. Tm window = same offsets above the higher
    primer, measured with Primer3's own Tm model.
    -> (candidates, info) ; info carries Primer3's own explanation of what it rejected."""
    try:
        import primer3  # noqa: PLC0415
    except ImportError:
        return [], {"error": "primer3-py is not installed (pip install primer3-py)"}
    conc = dict(mv_conc=50, dv_conc=0, dntp_conc=0, dna_conc=250)
    ref = max(primer3.calc_tm(cons[c[0]:c[1] + 1].replace("-", ""), **conc) for c in (f_col, r_col))
    cons_n = cons.replace("-", "N")                          # a consensus gap must not reach Primer3 (it only accepts ACGTN)
    seq_args = {"SEQUENCE_ID": "probe", "SEQUENCE_TEMPLATE": cons_n,
                "SEQUENCE_INTERNAL_EXCLUDED_REGION": [[0, lo], [hi + 1, len(cons) - hi - 1]]}
    glob = {"PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 0,
            "PRIMER_PICK_INTERNAL_OLIGO": 1, "PRIMER_INTERNAL_MIN_SIZE": min_len, "PRIMER_INTERNAL_MAX_SIZE": max_len,
            "PRIMER_INTERNAL_OPT_SIZE": min(max_len, max(min_len, 22)),
            "PRIMER_INTERNAL_MIN_TM": ref + above[0], "PRIMER_INTERNAL_MAX_TM": ref + above[1],
            "PRIMER_INTERNAL_OPT_TM": ref + PROBE_TM_TARGET,
            "PRIMER_INTERNAL_MIN_GC": gc_range[0], "PRIMER_INTERNAL_MAX_GC": gc_range[1], "PRIMER_INTERNAL_MAX_POLY_X": 3,
            "PRIMER_NUM_RETURN": max(10, min(50, quota * 4)), "PRIMER_INTERNAL_DNA_CONC": 250,
            "PRIMER_INTERNAL_SALT_MONOVALENT": 50, "PRIMER_INTERNAL_SALT_DIVALENT": 0, "PRIMER_INTERNAL_DNTP_CONC": 0}
    res = primer3.design_primers(seq_args, glob)
    info = {"ref_tm": round(ref, 1), "tm_range": [round(ref + above[0], 1), round(ref + above[1], 1)],
            "explain": res.get("PRIMER_INTERNAL_EXPLAIN", ""), "returned": res.get("PRIMER_INTERNAL_NUM_RETURNED", 0)}
    out = []
    for i in range(info["returned"]):
        start, length = res[f"PRIMER_INTERNAL_{i}"]
        s, e = start, start + length
        row_mm = [{j - s: ch for j, ch in enumerate(r["dots"][s:e], s) if ch != "."} for r in rows]
        if "-" in cons[s:e]:
            continue
        window, left = degenerate_window(cons[s:e], row_mm, max_iupac)
        cov0, cov1 = sum(m == 0 for m in left), sum(m <= 1 for m in left)
        for strand in ("+", "-"):
            oligo = window if strand == "+" else reverse_complement(window)
            rs = resolutions(oligo)
            if not rs or any(r[0] == "G" for r in rs):      # a 5' G quenches the reporter -- Primer3 has no such rule
                continue
            out.append({
                "seq": oligo, "strand": strand, "col_start": s + 1, "col_end": e, "length": length,
                "n_iupac": sum(c not in "ACGT" for c in window),
                "tm": round(res[f"PRIMER_INTERNAL_{i}_TM"], 1), "gc": round(res[f"PRIMER_INTERNAL_{i}_GC_PERCENT"], 1),
                "msa_cov0": cov0, "msa_cov1": cov1, "msa_n": n_rows, "source": "primer3",
                "p3_penalty": round(res[f"PRIMER_INTERNAL_{i}_PENALTY"], 2),
                "p3_hairpin": round(res.get(f"PRIMER_INTERNAL_{i}_HAIRPIN_TH", 0.0), 1),
            })
    return out, info


def mfe_verify(set_id: str, segment: str, candidates: list, outdir: Path, probe_mismatch: int, ids=None):
    """The pipeline's own probe test for each candidate: MFEprimer `search` over the whole filtered population."""
    db = ROOT / "data" / "filtered" / f"{segment.lower()}_filtered.fasta"
    total = len(read_fasta_ids(db)) if ids is None else len(ids)
    mfe = find_mfeprimer("", ROOT)
    ensure_index(mfe, db)
    for i, c in enumerate(candidates, 1):
        fa, txt = outdir / f"cand{i:02d}.fa", outdir / f"cand{i:02d}_mfe.txt"
        write_oligo_fasta(fa, [(f"{set_id}_CAND{i:02d}", c["seq"])])
        hits = run_probe_search(mfe, db, fa, txt, probe_mismatch)
        if ids is not None:
            hits = hits & set(ids)
        c["mfe_hits"], c["population"] = len(hits), total
        c["mfe_cov_pct"] = round(100.0 * len(hits) / total, 1) if total else 0.0


CAND_COLS = ["rank", "seq", "source", "strand", "length", "n_iupac", "tm", "tm_p3", "p3_penalty", "p3_hairpin", "gc", "col_start", "col_end", "msa_cov0", "msa_cov1", "msa_n",
             "mfe_hits", "population", "mfe_cov_pct"]


def run_design(set_id: str, top=10, min_len=18, max_len=30, tm_range=None, gc_range=(30.0, 80.0), probe_mismatch=0,
               use_mfe=True, max_iupac=2, tm_above=None, method="both", ids=None, scope=None) -> dict:
    """Design + test + save for one set; returns the summary dict that is also written to result.json."""
    import time  # noqa: PLC0415
    manifest = ROOT / "input_primers" / f"{set_id}.tsv"
    if not manifest.exists():
        raise FileNotFoundError(f"{set_id}: input_primers/{set_id}.tsv not found")
    with open(manifest, encoding="utf-8") as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    res = design(set_id, min_len, max_len, tm_range, gc_range, top, max_iupac, tm_above, method, ids)
    outdir = ROOT / "output" / "probe_design" / set_id
    outdir.mkdir(parents=True, exist_ok=True)
    if use_mfe and res["candidates"]:
        mfe_verify(set_id, row["segment"], res["candidates"], outdir, probe_mismatch, ids)
        res["candidates"].sort(key=lambda c: (-c["mfe_hits"], -c["msa_cov1"], abs(c["tm"] - (res["primer_tm"] + PROBE_TM_TARGET))))
    for i, c in enumerate(res["candidates"], 1):
        c["rank"] = i
    with open(outdir / "candidates.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CAND_COLS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(res["candidates"])
    msa_rows = [r for r in load_msa(set_id)["rows"] if ids is None or r["genome_id"] in ids]
    probe = row["probe"].upper()
    summary = {
        "set_id": set_id, "generated": time.strftime("%Y-%m-%d %H:%M:%S"), "host_scope": scope,
        "population_n": None if ids is None else len(ids),
        "primer_tm": round(res["primer_tm"], 1), "tm_range": [round(x, 1) for x in res["tm_range"]],
        "gc_range": list(gc_range), "len_range": [min_len, max_len], "probe_mismatch": probe_mismatch, "max_iupac": max_iupac, "method": method, "primer3": res.get("primer3"),
        "tm_above": [round(res["tm_range"][0] - res["primer_tm"], 1), round(res["tm_range"][1] - res["primer_tm"], 1)],
        "col_window": list(res["window"]), "n_pass_rules": res["n_total"], "n_rejected": res["n_rejected"],
        "msa_n": res["n_rows"], "mfe_tested": bool(use_mfe),
        "registered": {"seq": probe,
                       "tm": round(tm_nn(probe), 1) if set(probe) <= set("ACGT") else None,
                       "msa_cov0": sum(str(r.get("p_mismatch")) == "0" for r in msa_rows)},
        "candidates": res["candidates"],
    }
    (outdir / "result.json").write_text(json.dumps(summary, ensure_ascii=False), encoding="utf-8")
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("set_id")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--min-len", type=int, default=18)
    ap.add_argument("--max-len", type=int, default=30)
    ap.add_argument("--tm", type=float, nargs=2, default=None, metavar=("MIN", "MAX"),
                    help="absolute probe Tm window (overrides --tm-above)")
    ap.add_argument("--tm-above", type=float, nargs=2, default=None, metavar=("MIN", "MAX"),
                    help="probe Tm window relative to the HIGHER primer Tm (default: 3 10)")
    ap.add_argument("--gc", type=float, nargs=2, default=(30.0, 80.0), metavar=("MIN", "MAX"))
    ap.add_argument("--max-iupac", type=int, default=2, help="max IUPAC (degenerate) positions per probe (default 2; 0 = plain consensus)")
    ap.add_argument("--method", choices=("own", "primer3", "both"), default="both",
                    help="candidate source: this module's own rules, Primer3's internal-oligo design, or both (default)")
    ap.add_argument("--probe-mismatch", type=int, default=0)
    ap.add_argument("--no-mfe", action="store_true", help="skip the MFEprimer population test (alignment score only)")
    a = ap.parse_args()
    r = run_design(a.set_id, a.top, a.min_len, a.max_len, tuple(a.tm) if a.tm else None, tuple(a.gc),
                   a.probe_mismatch, not a.no_mfe, a.max_iupac,
                   tuple(a.tm_above) if a.tm_above else None, a.method)
    print(f"{a.set_id}: columns {r['col_window'][0]}-{r['col_window'][1]}, {r['n_pass_rules']} windows pass "
          f"({r['n_rejected']} rejected), higher primer Tm {r['primer_tm']}, probe Tm {r['tm_range'][0]}-{r['tm_range'][1]}")
    print(f"registered probe: {r['registered']['seq']}  Tm {r['registered']['tm']}  MSA 0mm {r['registered']['msa_cov0']}/{r['msa_n']}")
    for c in r["candidates"]:
        mfe = f"  MFE {c['mfe_hits']}/{c['population']} ({c['mfe_cov_pct']}%)" if "mfe_hits" in c else ""
        print(f"{c['rank']:>2}. [{c['source']}] {c['seq']:<32} {c['strand']} len{c['length']} IUPAC{c['n_iupac']} Tm{c['tm']} GC{c['gc']:.0f}%  "
              f"MSA 0mm {c['msa_cov0']}/{c['msa_n']}, <=1mm {c['msa_cov1']}{mfe}")


if __name__ == "__main__":
    main()
