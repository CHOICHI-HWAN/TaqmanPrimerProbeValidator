#!/usr/bin/env python3
"""Marks the top 3 sets of every segment of this project with the medals at the end of the set name (config/primers.tsv), from the finished results.

    python3 scripts/rank_top3.py                # rank, write warnings into the remarks
    python3 scripts/rank_top3.py --dry-run
    python3 scripts/rank_top3.py --max-amplicon 0 --no-rules      # the plain ranking (coverage, probe Tm) only

Ranking per segment (host scope human, else the first scope the summary has): validated coverage (output/summary/summary.tsv, whole percent) first,
sets within the same whole percent are ordered by the probe Tm (higher first). A set can rank only if it validates >= 50% of the genomes, its probe Tm is above
both primers (a MGB probe gets +10 C, a rough value), its amplicon on the reference is <= 200 bp (--max-amplicon) and it passes the hard rule of
oligo_rules.py (probe 5' end is not G with a FAM-like reporter). The other rules are only written as a warning into the remark of the
tool-made sets (AUTO-DESIGN / -OPT-); a user's own sets get only the amplicon warning appended and keep their sequences, names and remarks.
Old medals are removed first, so the tool can be run again after every analysis.
"""
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import oligo_rules  # noqa: E402
import probe_designer as pd  # noqa: E402

MEDALS = ["\U0001F947", "\U0001F948", "\U0001F949"]
MGB_OFFSET = 10.0
MIN_PCT = 50.0                                    # a set that validates fewer genomes than this gets no medal, however few sets there are
MEDAL_RE = re.compile("\\s*(?:\U0001F947|\U0001F948|\U0001F949)")              # old medals anywhere in the name (also one copied into an -OPT- name)


def tm(seq):
    res = pd.resolutions(seq, 64)
    return sum(pd.tm_nn(x) for x in res) / len(res) if res else 0.0


def amplicon_len(proj, pid):
    """F start .. R end on the reference (output/primers/<id>/refgenome_binding/binding.tsv); None when unknown."""
    bf = proj / "output" / "primers" / pid / "refgenome_binding" / "binding.tsv"
    try:
        b = {x["role"]: x for x in csv.DictReader(bf.open(), delimiter="\t")}
        pos = [int(b[k][e]) for k in "FR" for e in ("start", "end")]
        return max(pos) - min(pos) + 1
    except (OSError, KeyError, ValueError):
        return None


def main():
    dry = "--dry-run" in sys.argv
    max_amp = int(sys.argv[sys.argv.index("--max-amplicon") + 1]) if "--max-amplicon" in sys.argv else 200   # 0 = no amplicon limit
    rules = "--no-rules" not in sys.argv
    proj = HERE.parent
    for f in [proj / "config" / "primers.tsv"]:
        summ = proj / "output" / "summary" / "summary.tsv"
        rows = list(csv.DictReader(f.open(encoding="utf-8"), delimiter="\t"))
        if not rows or not summ.exists():
            continue
        for r in rows:
            r["primer_name"] = MEDAL_RE.sub("", r["primer_name"])
        by_id = {r["primer_id"]: r for r in rows}
        sm = list(csv.DictReader(summ.open(encoding="utf-8"), delimiter="\t"))
        scopes = [s["host_scope"] for s in sm]
        scope = "human" if "human" in scopes else (scopes[0] if scopes else "")
        medals = {}
        for seg in sorted({s["segment"] for s in sm}):
            cand = []
            for s in sm:
                r = by_id.get(s["set_id"])
                if not r or s["segment"] != seg or s["host_scope"] != scope or not s.get("validated_pct"):
                    continue
                if rules and oligo_rules.check(r)[0]:
                    continue
                f_tm, r_tm, p_tm = tm(r["forward_primer"]), tm(r["reverse_primer"]), tm(r["probe"])
                if "MGB" in r.get("fluorophore", "").upper():
                    p_tm += MGB_OFFSET                          # a minor-groove binder raises the probe Tm (rough, same value as chichipes --probe-chem mgb)
                if p_tm <= max(f_tm, r_tm):
                    continue
                if max_amp:
                    n = amplicon_len(proj, r["primer_id"])
                    if n and n > max_amp:
                        continue
                pct = float(s["validated_pct"])
                if pct >= MIN_PCT:
                    cand.append((-round(pct), -p_tm, r["primer_id"], pct, p_tm))
            cand.sort()
            seen, uniq = set(), []
            for c in cand:                                    # the same oligo trio under several names counts once
                r = by_id[c[2]]
                key = (r["forward_primer"], r["reverse_primer"], r["probe"])
                if key not in seen:
                    seen.add(key)
                    uniq.append(c)
            for i, c in enumerate(uniq[:3]):
                medals[c[2]] = MEDALS[i]
                print(f"{seg:7} {MEDALS[i]} {c[2]} {by_id[c[2]]['primer_name'][:40]:40} {c[3]:6.2f}%  probe Tm {c[4]:.1f}")
        if rules or max_amp:
            for r in rows:
                own = r["source"] != "AUTO-DESIGN" and "-OPT-" not in r["primer_name"]   # a user's own set keeps its remark as typed
                w = []
                if rules and not own and r["segment"] and r["proceed"].strip().lower() == "true":
                    h, w = oligo_rules.check(r, tm)
                    w = [oligo_rules.HARD_TEXT[c] + " (메달 제외)" for c in h] + w
                if max_amp and r["segment"] and r["proceed"].strip().lower() == "true":
                    n = amplicon_len(proj, r["primer_id"])
                    if n and n > max_amp:
                        w.append(f"amplicon {n} bp (상한 {max_amp} bp 초과, 메달 제외)")
                r["remark"] = oligo_rules.with_warnings(r.get("remark", ""), w)
        for pid, m in medals.items():
            by_id[pid]["primer_name"] += f" {m}"
        if not dry:
            with f.open("w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)


if __name__ == "__main__":
    main()
