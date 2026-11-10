#!/usr/bin/env python3
"""Primer/probe secondary-structure QC: self-dimers, cross-dimers (F-R, F-P, P-R) and
hairpins, using mfeprimer's own `dimer`/`hairpin` subcommands (the same binary already used
for population matching elsewhere in this pipeline).

Why this exists: every other check in this pipeline asks "does this oligo bind the intended
target" -- mismatch count, Tm, dG against the population/reference genome. None of that
catches the classic wet-lab failure mode where the oligos bind EACH OTHER (or fold back on
themselves) instead of the target, which can kill amplification efficiency even when every
population-match metric looks perfect. mfeprimer expands IUPAC degenerate bases into every
concrete variant before checking (a degenerate primer really is synthesized as that pool of
variants), so a single set can produce many pairwise combinations; this script summarizes
that down to "worst case per role-pair" rather than dumping every combination.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import find_mfeprimer  # noqa: E402

# Rule-of-thumb flagging threshold (kcal/mol) -- a dimer/hairpin more stable (more negative)
# than this is generally considered a real risk to amplification efficiency by common primer
# QC tools (IDT OligoAnalyzer, Primer3Plus use similar cutoffs); like the target-binding dG
# threshold already used elsewhere in this app, it's a reference line, not a hard pass/fail.
DG_FLAG_THRESHOLD = -9.0

ROLE_LABELS = {"F": "Forward", "P": "Probe", "R": "Reverse"}


def run_mfeprimer(mfeprimer_bin, subcmd, fasta_path, out_path, extra_args=None):
    cmd = [mfeprimer_bin, subcmd, "-i", str(fasta_path), "-o", str(out_path)] + (extra_args or [])
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] mfeprimer {subcmd} failed: {result.stderr}")
    return out_path.read_text(encoding="utf-8")


def _role(oligo_id):
    return oligo_id.split(".")[0]


def parse_dimer_report(text):
    m = re.search(r"Dimer List \((\d+)\)", text)
    if not m or int(m.group(1)) == 0:
        return []
    body_m = re.search(r"Dimer List \(\d+\)\n\n(.*?)\n+Parameters", text, re.DOTALL)
    if not body_m:
        return []
    entries = re.split(r"\n\n(?=Dimer \d+:)", body_m.group(1).strip())
    out = []
    for entry in entries:
        em = re.match(
            r"Dimer \d+: (\S+) x (\S+)\s*\n\s*Score: (\d+), Delta G = (-?[\d.]+) kcal/mol\s*(.*)",
            entry, re.DOTALL,
        )
        if not em:
            continue
        id1, id2, score, dg, align = em.groups()
        out.append({
            "role1": _role(id1), "role2": _role(id2),
            "score": int(score), "dg": float(dg), "align": align.strip(),
        })
    return out


def parse_hairpin_report(text):
    m = re.search(r"Hairpin List \((\d+)\)", text)
    if not m or int(m.group(1)) == 0:
        return []
    body_m = re.search(r"Hairpin List \(\d+\)\n\n(.*?)\n+Parameters", text, re.DOTALL)
    if not body_m:
        return []
    entries = re.split(r"\n\n(?=Hairpin \d+:)", body_m.group(1).strip())
    out = []
    for entry in entries:
        em = re.match(
            r"Hairpin \d+: (\S+)\s*\n\s*Score: (\d+), Tm = (-?[\d.]+) .C, Delta G = (-?[\d.]+) kcal/mol\s*(.*)",
            entry, re.DOTALL,
        )
        if not em:
            continue
        oid, score, tm, dg, align = em.groups()
        out.append({
            "role": _role(oid), "score": int(score),
            "tm": float(tm), "dg": float(dg), "align": align.strip(),
        })
    return out


def summarize_dimers(dimers):
    """Worst (most negative dG) entry per unordered role-pair, plus a flagged count."""
    by_pair = {}
    for d in dimers:
        pair = "-".join(sorted((d["role1"], d["role2"])))
        slot = by_pair.setdefault(pair, {"worst": None, "flagged": 0})
        if d["dg"] <= DG_FLAG_THRESHOLD:
            slot["flagged"] += 1
        if slot["worst"] is None or d["dg"] < slot["worst"]["dg"]:
            slot["worst"] = d
    return by_pair


def summarize_hairpins(hairpins):
    by_role = {}
    for h in hairpins:
        slot = by_role.setdefault(h["role"], {"worst": None, "flagged": 0})
        if h["dg"] <= DG_FLAG_THRESHOLD:
            slot["flagged"] += 1
        if slot["worst"] is None or h["dg"] < slot["worst"]["dg"]:
            slot["worst"] = h
    return by_role


def main():
    f_seq = str(getattr(snakemake.params, "f_seq", "") or "")
    p_seq = str(getattr(snakemake.params, "p_seq", "") or "")
    r_seq = str(getattr(snakemake.params, "r_seq", "") or "")
    pipeline_root = Path(snakemake.scriptdir).parent
    mfeprimer_bin = find_mfeprimer(snakemake.config.get("mfeprimer_bin", ""), pipeline_root)
    cpu = int(snakemake.config.get("mfeprimer_cpu", 0) or 0) or None

    out_json = Path(snakemake.output.json)
    out_dir = out_json.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    records = [(role, seq) for role, seq in (("F", f_seq), ("P", p_seq), ("R", r_seq)) if seq]
    if not records:
        out_json.write_text(json.dumps({"dimer": {}, "hairpin": {}, "dg_flag_threshold": DG_FLAG_THRESHOLD}), encoding="utf-8")
        return

    fasta_path = out_dir / "oligos.fasta"
    fasta_path.write_text("".join(f">{role}\n{seq}\n" for role, seq in records), encoding="utf-8")

    dimer_txt = out_dir / "dimer.txt"
    hairpin_txt = out_dir / "hairpin.txt"
    # "hairpin" has no --cpu flag at all; "dimer" does (default 2).
    dimer_extra = ["--cpu", str(cpu)] if cpu else None
    dimer_report = run_mfeprimer(mfeprimer_bin, "dimer", fasta_path, dimer_txt, extra_args=dimer_extra)
    hairpin_report = run_mfeprimer(mfeprimer_bin, "hairpin", fasta_path, hairpin_txt)

    dimers = parse_dimer_report(dimer_report)
    hairpins = parse_hairpin_report(hairpin_report)

    payload = {
        "dimer": {
            pair: {
                "worst_dg": v["worst"]["dg"], "worst_score": v["worst"]["score"],
                "align": v["worst"]["align"], "flagged_count": v["flagged"],
            }
            for pair, v in summarize_dimers(dimers).items()
        },
        "hairpin": {
            role: {
                "worst_dg": v["worst"]["dg"], "worst_tm": v["worst"]["tm"],
                "worst_score": v["worst"]["score"], "align": v["worst"]["align"],
                "flagged_count": v["flagged"],
            }
            for role, v in summarize_hairpins(hairpins).items()
        },
        "dg_flag_threshold": DG_FLAG_THRESHOLD,
    }
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"[dimer_hairpin] {len(dimers)} dimer combinations, {len(hairpins)} hairpin combinations "
          f"-> {out_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
