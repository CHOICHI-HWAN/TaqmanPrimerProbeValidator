"""Parse `mfeprimer dimer` / `hairpin` reports into records (for the black-box analysis of their rules).

Dimer entry:
    Dimer 1: A x B
    <blank>
        Score: 6, Delta G = -3.61 kcal/mol
    <blank>
    5' AAAA...A 3'          <- A, 5'->3'
           ::::::           <- marks under A's columns (':' pair, '.' mismatch)
       3' BBBB...B 5'       <- B written 3'->5' (reversed), shifted
Returns per dimer: names, score, dg, a (seq), b (seq, 5'->3'), shift (column of B's first written base minus A's
first base; B is written reversed), marks (string aligned to A's columns, '' padded).
"""
import re

DIMER_RE = re.compile(r"Dimer (\d+): (\S+) x (\S+)\s*\n\s*Score: (\d+), Delta G = (-?[\d.]+) kcal/mol\s*\n[ \t]*\n"
                      r"( *5' .*? 3')\n(.*)\n( *3' .*? 5')")
HAIRPIN_RE = re.compile(r"Hairpin (\d+): (\S+)\s*\n\s*Score: (\d+), Tm = (-?[\d.]+) .C, Delta G = (-?[\d.]+) kcal/mol"
                        r"\s*\n[ \t]*\n[ \t]*(\S+)\n[ \t]*(\S+)")


def dimers(text):
    out = []
    for m in DIMER_RE.finditer(text):
        _, a_name, b_name, score, dg, l1, l2, l3 = m.groups()
        a1 = l1.lstrip()
        alead = len(l1) - len(a1)      # column of "5'"
        a = a1[3:-3]
        b3 = l3.lstrip()
        blead = len(l3) - len(b3)
        brev = b3[3:-3]
        a0, b0 = alead + 3, blead + 3  # columns of the first base of A and of reversed B
        shift = b0 - a0                # B-rev column j-shift faces A column j
        marks = l2[a0:] if len(l2) > a0 else ""
        out.append({"a_name": a_name, "b_name": b_name, "score": int(score), "dg": float(dg), "a": a,
                    "b": brev[::-1], "brev": brev, "shift": shift, "marks": marks})
    return out


def hairpins(text):
    out = []
    for m in HAIRPIN_RE.finditer(text):
        _, name, score, tm, dg, struct, seq = m.groups()
        out.append({"name": name, "score": int(score), "tm": float(tm), "dg": float(dg), "struct": struct, "seq": seq})
    return out
