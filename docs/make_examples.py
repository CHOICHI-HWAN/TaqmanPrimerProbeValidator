#!/usr/bin/env python3
"""Draw the two "Exact mismatch counting" examples of the README from the pipeline's own results.

    python3 docs/make_examples.py            # writes docs/media/example_hidden_mismatch.svg and example_missed_site.svg
    python3 docs/make_examples.py --list     # only print which sets / genomes were picked

The pictures contain NO nucleotides: each primer position is a square (grey = matches the genome, red = does not), so the real
binding result can be shown in the public repository without publishing a primer or genome sequence. Set names and vendors
are not printed either (sets are called "set A"); only the public genome accession is shown.

  1. hidden mismatch: mfeprimer's report text is the primer itself, so only mismatches inside its 3' seed show up in it; the
     recount against the genome bases finds the mismatch further toward the 5' end.
  2. missed site: mfeprimer found no amplicon at all in this genome (its 3' k-mer seed does not hit), the independent
     whole-genome scan finds a primer pair that passes.
Needs the results of a finished run (data/filtered, output/primers).
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import mfe_common as M  # noqa: E402

SEED = 9            # mfeprimer's default 3' seed (--misEnd): mismatches inside it are visible in its report text
THREE_P = 2         # config three_prime_consecutive_match


def read_tsv(p):
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def vector(oligo, genome, pos):
    """Per-position match list of an oligo against the genome site `pos` (1-based inclusive, strand '+'/'-'), read 5'->3'
    along the OLIGO. True = matches."""
    site = genome[pos["start"] - 1:pos["end"]].upper()
    q = oligo.upper() if pos["strand"] == "+" else M.reverse_complement(oligo.upper())
    ok = [t in M.IUPAC.get(c, "ACGT") for c, t in zip(q, site)]
    return ok if pos["strand"] == "+" else ok[::-1]        # '-' hits are compared on the top strand: flip back to 5'->3' of the oligo


def load(set_id, prim):
    seg = prim[set_id]["segment"].lower()
    seqs = M.read_fasta_seqs(ROOT / f"data/filtered/{seg}_filtered.fasta")
    return seqs, prim[set_id]["forward_primer"].upper(), prim[set_id]["reverse_primer"].upper(), prim[set_id]["probe"].upper()


def pick(prim):
    hidden = missed = None
    for set_id in sorted(prim):
        spec = ROOT / f"output/primers/{set_id}/primer_qc/{set_id}_mm2.txt.spec.tsv"
        cmp_ = ROOT / f"output/primers/{set_id}/alignment/engine_comparison.tsv"
        if not (spec.exists() and cmp_.exists()):
            continue
        seqs, F, R, P = load(set_id, prim)
        if True:
            rep = M.parse_spec_best(spec, f"{set_id}_F", F, f"{set_id}_R", R, THREE_P)
            rec = M.parse_spec_best(spec, f"{set_id}_F", F, f"{set_id}_R", R, THREE_P, genome_seqs=seqs)
            cand = []
            for g in sorted(rep):
                if g in rec and rep[g]["f_mismatch"] <= 1 and rep[g]["r_mismatch"] <= 1 and rec[g]["f_3p_ok"] and rec[g]["r_3p_ok"]:
                    d = (rec[g]["f_mismatch"] - rep[g]["f_mismatch"]) + (rec[g]["r_mismatch"] - rep[g]["r_mismatch"])
                    if d == 1 and max(rec[g]["f_mismatch"], rec[g]["r_mismatch"]) == 2:
                        cand.append((g, rep[g], rec[g]))
            if cand and (hidden is None or len(cand) > hidden[1]):
                hidden = (set_id, len(cand), cand[0])
        if missed is None:
            rows = [r for r in read_tsv(cmp_) if r["agree"] == "0" and r["source"] == "genome-scan" and not r["mfe_f"]
                    and r["own_f"] and int(r["own_f"]) <= 1 and int(r["own_r"]) <= 1 and r["own_3p_ok"] == "1"]
            if rows:
                missed = (set_id, len(rows), rows[0])
    return hidden, missed


# ---------------------------------------------------------------- drawing
GREY, RED, INK, SOFT = "#c9d1d9", "#d1242f", "#1f2328", "#59636e"
BLUE, GREEN, ORANGE = "#0969da", "#1a7f37", "#bc4c00"


def squares(x, y, ok, size=13, gap=3, mark3=True):
    out = []
    for i, good in enumerate(ok):
        out.append(f'<rect x="{x + i * (size + gap)}" y="{y}" width="{size}" height="{size}" rx="2" fill="{GREY if good else RED}"/>')
    return "".join(out), len(ok) * (size + gap) - gap


def text(x, y, s, size=13, weight="400", fill=INK, anchor="start"):
    s = s.replace("&", "&amp;").replace("<", "&lt;")
    return f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{s}</text>'


def primer_row(x, y, label, ok, note, seed=True):
    sq, w = squares(x + 150, y, ok)
    n = len(ok)
    s = [text(x, y + 11, label, 13, "600"), sq]
    # 5' / 3' labels and the seed bracket
    s.append(text(x + 150, y - 5, "5′", 11, "400", SOFT))
    s.append(text(x + 150 + w, y - 5, "3′", 11, "400", SOFT, "end"))
    if seed:
        sx = x + 150 + (n - SEED) * 16
        s.append(f'<line x1="{sx}" y1="{y + 21}" x2="{x + 150 + w}" y2="{y + 21}" stroke="{SOFT}" stroke-width="1.5"/>')
        s.append(text(sx, y + 34, f"3′ seed ({SEED} nt): the only part mfeprimer's text reflects", 11, "400", SOFT))
    s.append(text(x + 150 + w + 16, y + 11, note, 13, "400", INK))
    return "".join(s)


def amplicon(x, y, width, amp, probe, flen, rlen, fpos, rpos):
    a0, a1 = amp
    span = max(a1 - a0, 1)
    sx = lambda p: x + (p - a0) / span * width   # noqa: E731
    s = [f'<line x1="{x}" y1="{y + 8}" x2="{x + width}" y2="{y + 8}" stroke="{SOFT}" stroke-width="2"/>']
    for (p0, p1, col, lab) in ((fpos[0], fpos[1], BLUE, "F"), (probe[0], probe[1], GREEN, "Probe"), (rpos[0], rpos[1], ORANGE, "R")):
        x0, x1 = sx(p0), sx(p1)
        s.append(f'<rect x="{x0:.1f}" y="{y}" width="{max(x1 - x0, 4):.1f}" height="16" rx="2" fill="{col}"/>')
        s.append(text(f"{(x0 + x1) / 2:.1f}", y + 32, lab, 11, "400", col, "middle"))
    flipped = "  ·  this genome is deposited in the opposite orientation, so R is drawn on the left" if fpos[0] > rpos[0] else ""
    s.append(text(x, y - 6, f"amplicon {a0:,}–{a1:,}  ({a1 - a0 + 1} bp){flipped}", 11, "400", SOFT))
    return "".join(s)


def card(title, sub, body, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 820 {h}" width="820" height="{h}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">'
            f'<rect x="0.5" y="0.5" width="819" height="{h - 1}" rx="8" fill="#ffffff" stroke="#d0d7de"/>'
            + text(20, 28, title, 15, "700") + text(20, 48, sub, 12, "400", SOFT) + body + "</svg>\n")


def verdict(x, y, who, lines, good, color):
    return (f'<rect x="{x}" y="{y}" width="360" height="{34 + 18 * len(lines)}" rx="6" fill="none" stroke="{color}" stroke-width="1.5"/>'
            + text(x + 12, y + 20, who, 13, "700", color)
            + "".join(text(x + 12, y + 40 + 18 * i, l, 12, "400", INK) for i, l in enumerate(lines)))


def example_hidden(set_id, genome, rep, rec, prim):
    seqs, F, R, P = load(set_id, prim)
    g = seqs[genome]
    # the oligo whose count differs
    role = "F" if rec["f_mismatch"] > rep["f_mismatch"] else "R"
    oligo = F if role == "F" else R
    pos = rec["f_pos"] if role == "F" else rec["r_pos"]
    real = vector(oligo, g, pos)
    seen = [ok if i >= len(real) - SEED else True for i, ok in enumerate(real)]     # what the report text can show
    mm_pos = [i + 1 for i, ok in enumerate(real) if not ok]
    row = "forward primer" if role == "F" else "reverse primer"
    amp = (int(rec["amp_start"]), int(rec["amp_end"]))
    pr = [r for r in read_tsv(ROOT / f"output/primers/{set_id}/alignment/mismatch_table.tsv") if r["genome_id"] == genome][0]
    body = amplicon(40, 100, 740, amp, (int(pr["probe_start"]), int(pr["probe_end"])), len(F), len(R),
                    (rec["f_pos"]["start"], rec["f_pos"]["end"]), (rec["r_pos"]["start"], rec["r_pos"]["end"]))
    body += primer_row(20, 165, f"as mfeprimer reports it", seen, f"{row}: {rep['f_mismatch' if role == 'F' else 'r_mismatch']} mismatch", seed=False)
    body += primer_row(20, 215, "on the genome bases", real, f"{row}: {len(mm_pos)} mismatches")
    body += verdict(20, 275, "mfeprimer only", [f"forward {rep['f_mismatch']} / reverse {rep['r_mismatch']} mismatches", "→ passes the ≤ 1 rule"], True, "#cf222e")
    body += verdict(430, 275, "this tool (recount)", [f"forward {rec['f_mismatch']} / reverse {rec['r_mismatch']} mismatches", "→ NOT validated (needs ≤ 1 each)"], False, "#1a7f37")
    return card(f"1 · A mismatch that mfeprimer's report does not show — genome {genome}",
                "Real result of set A against a public genome. Grey = base matches the genome, red = mismatch. No sequences are drawn.", body, 365)


def example_missed(set_id, genome, row, prim):
    seqs, F, R, P = load(set_id, prim)
    g = seqs[genome]
    own = M.scan_fr(F, R, g, three_prime_n=THREE_P)
    fv, rv = vector(F, g, own["f_pos"]), vector(R, g, own["r_pos"])
    pr = [r for r in read_tsv(ROOT / f"output/primers/{set_id}/alignment/mismatch_table.tsv") if r["genome_id"] == genome][0]
    amp = (int(own["amp_start"]), int(own["amp_end"]))
    body = amplicon(40, 100, 740, amp, (int(pr["probe_start"]), int(pr["probe_end"])), len(F), len(R),
                    (own["f_pos"]["start"], own["f_pos"]["end"]), (own["r_pos"]["start"], own["r_pos"]["end"]))
    body += primer_row(20, 165, "forward primer", fv, f"{sum(1 for o in fv if not o)} mismatch, 3′ end intact", seed=False)
    body += primer_row(20, 200, "reverse primer", rv, f"{sum(1 for o in rv if not o)} mismatch, 3′ end intact", seed=False)
    body += verdict(20, 250, "mfeprimer only", ["no amplicon found in this genome", "→ counted as NOT detected"], False, "#cf222e")
    body += verdict(430, 250, "this tool (genome scan)", [f"forward {own['f_mismatch']} / reverse {own['r_mismatch']} mismatches, 3′ ends intact",
                                                         "probe exact inside the amplicon → validated"], True, "#1a7f37")
    return card(f"2 · A binding site that mfeprimer misses — genome {genome}",
                "Real result of set B. mfeprimer finds sites through a 3′ k-mer seed; the whole-genome scan checks every window.", body, 345)


def main():
    prim = {r["primer_id"]: r for r in read_tsv(ROOT / "config/primers.tsv")}
    hidden, missed = pick(prim)
    print("hidden mismatch:", hidden and (hidden[0], hidden[1], hidden[2][0]))
    print("missed site    :", missed and (missed[0], missed[1], missed[2]["genome_id"]))
    if "--list" in sys.argv:
        return
    out = ROOT / "docs" / "media"
    if hidden:
        set_id, n, (genome, rep, rec) = hidden
        (out / "example_hidden_mismatch.svg").write_text(example_hidden(set_id, genome, rep, rec, prim), encoding="utf-8")
        print(f"  {n} genomes of that set behave like this")
    if missed:
        set_id, n, row = missed
        (out / "example_missed_site.svg").write_text(example_missed(set_id, row["genome_id"], row, prim), encoding="utf-8")
        print(f"  {n} genomes of that set behave like this")


if __name__ == "__main__":
    main()
