"""Probe / primer design rules checked on a registered set (used by tools/rank_top3.py --rules).

HARD rules: a set that breaks one cannot get a medal.
  H1 probe 5' end is G while the reporter is FAM-like (G next to FAM quenches it)
WARN rules: only a warning text is written into the set's remark (no effect on the ranking).
  W1 probe has more G than C
  W2 run of >= 4 identical bases in an oligo (G runs in the probe are the worst)
  W3 more than 3 G/C among the last 5 bases of a primer (GC clamp)
  W4 forward/reverse primer Tm differ by more than 2 C
  W5 primer GC outside 40-60 % / probe GC outside 30-70 %
  W6 IUPAC code among the 3' last 2..5 bases of a primer
  W0 the primer's 3'-terminal base is an IUPAC code (the 3' end decides extension; kept a warning because the best
     LASV L sets reach their coverage with it, so a hard rule would remove them)
  W7 an oligo expands to more than 48 sequences (many ambiguity codes dilute each variant)
"""
import re

MARK = "⚠ 설계 규칙 경고:"
AMB = {"R": 2, "Y": 2, "S": 2, "W": 2, "K": 2, "M": 2, "B": 3, "D": 3, "H": 3, "V": 3, "N": 4}
FAM_LIKE = ("FAM", "SYBR")
HARD_TEXT = {"H1": "프로브 5' 말단이 G (FAM 소광)"}


def gc(s):
    s = s.upper()
    return 100.0 * sum(c in "GCS" for c in s) / len(s) if s else 0.0


def multiplicity(s):
    n = 1
    for c in s.upper():
        n *= AMB.get(c, 1)
    return n


def has_run(s, n=4):
    return re.search(r"([ACGT])\1{%d,}" % (n - 1), s.upper()) is not None


def check(row, tm=None):
    """row: a primers.tsv row (dict). tm: optional function seq -> Tm. Returns (hard codes, warning texts)."""
    f, r, p = (row["forward_primer"].upper(), row["reverse_primer"].upper(), row["probe"].upper())
    hard, warn = [], []
    if p and p[0] == "G" and any(k in row.get("fluorophore", "").upper() for k in FAM_LIKE):
        hard.append("H1")
    for name, pr in (("F", f), ("R", r)):
        if pr and pr[-1] in AMB:
            warn.append(f"{name} 프라이머 3' 말단 염기가 혼합 염기 (3' 말단은 신장을 좌우함, 실험 확인 권장)")
    if p.count("G") > p.count("C"):
        warn.append("프로브 G가 C보다 많음 (형광 소광 우려)")
    for name, o in (("F", f), ("R", r), ("P", p)):
        if has_run(o, 4):
            warn.append(f"{name}에 같은 염기 4개 이상 연속")
    for name, o in (("F", f), ("R", r)):
        if sum(c in "GCS" for c in o[-5:]) > 3:
            warn.append(f"{name} 3' 말단 5염기 중 G/C {sum(c in 'GCS' for c in o[-5:])}개 (GC clamp 과다)")
        if any(c in AMB for c in o[-5:-1]):
            warn.append(f"{name} 3' 말단 근처(2~5번째)에 혼합 염기")
        if not 40 <= gc(o) <= 60:
            warn.append(f"{name} GC {gc(o):.0f}% (40~60 밖)")
    if p and not 30 <= gc(p) <= 70:
        warn.append(f"프로브 GC {gc(p):.0f}% (30~70 밖)")
    if tm:
        d = abs(tm(f) - tm(r))
        if d > 2.0:
            warn.append(f"F/R Tm 차이 {d:.1f}℃ (2℃ 초과)")
    for name, o in (("F", f), ("R", r), ("P", p)):
        if multiplicity(o) > 48:
            warn.append(f"{name} 조합 수 {multiplicity(o)} (48 초과)")
    return hard, warn


def strip_marker(remark):
    return re.sub(r"\s*\|?\s*" + re.escape(MARK) + r".*?(?=\s\|\s|$)", "", remark or "").strip()


def with_warnings(remark, warn):
    base = strip_marker(remark)
    if not warn:
        return base
    return (base + " | " if base else "") + MARK + " " + "; ".join(warn)
