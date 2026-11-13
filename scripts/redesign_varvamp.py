#!/usr/bin/env python3
"""Primer/probe redesign from the whole population with varVAMP (qPCR mode), one lineage (clade) at a time.

For each lineage: take its genomes (optionally only those of a host scope), align them with MAFFT, run
`varvamp qpcr` (degenerate F/R primers + an internal TaqMan probe that is one ambiguity stricter than the primers),
lowering the consensus threshold step by step until varVAMP finds a scheme. The designs are only *candidates*: the server
then runs every one of them through the pipeline's own validation (live_analyze.analyze_set).

    python3 scripts/redesign_varvamp.py S LII [--threshold 0.9 --min-threshold 0.6 --ambig 2 --probe-ambig 1 --top 3]

Everything is written under output/redesign/<segment>/<lineage>/ (git-ignored like the rest of output/).
"""
import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from mfe_common import find_mafft, read_fasta_seqs  # noqa: E402

_ANSI = re.compile(r"\x1b\[[0-9;]*m")


def find_varvamp() -> str:
    """varvamp from the same environment as the running Python (pip install varvamp), else PATH."""
    cand = Path(sys.executable).parent / "varvamp"
    if cand.exists():
        return str(cand)
    found = shutil.which("varvamp")
    if not found:
        sys.exit("[error] varvamp not found (pip install varvamp, see environment.yml).")
    return found


def _clean(text: str) -> str:
    return _ANSI.sub("", text.replace("\r", "\n"))


def _parse(out: Path) -> list:
    """qpcr_primers.tsv (+ qpcr_design.tsv for amplicon penalty/deltaG) -> [{name, F, P, R, ...}] in varVAMP's order."""
    designs = {}
    with open(out / "qpcr_primers.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d = designs.setdefault(r["qpcr_scheme"], {"name": r["qpcr_scheme"]})
            key = {"LEFT": "F", "PROBE": "P", "RIGHT": "R"}[r["oligo_type"]]
            d[key] = r["seq"].upper()
            d[key + "_tm"] = float(r["temp_best"])
            d[key + "_gc"] = float(r["gc_best"])
            d[key + "_pen"] = float(r["penalty"])
            d[key + "_start"], d[key + "_stop"] = int(r["start"]), int(r["stop"])
    design_tsv = out / "qpcr_design.tsv"
    if design_tsv.exists():
        with open(design_tsv, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                d = designs.get(r["qpcr_scheme"])
                if d:
                    d["penalty"], d["deltaG"], d["amplicon_len"] = float(r["penalty"]), float(r["deltaG"]), int(r["length"])
    return [d for d in designs.values() if all(k in d for k in ("F", "P", "R"))]


def design_lineage(members: dict, out_dir: Path, threshold=0.9, min_threshold=0.6, step=0.1, ambig=2, probe_ambig=1,
                   top=3, threads=4, progress=None) -> dict:
    """members = {accession: sequence}. -> {"status": "ok"|"none"|"error", "threshold", "n", "designs": [...], "message"}"""
    def say(msg):
        if progress:
            progress(msg)

    out_dir.mkdir(parents=True, exist_ok=True)
    n = len(members)
    if n < 3:
        return {"status": "error", "n": n, "designs": [], "message": f"게놈이 {n}개뿐이라 정렬·설계를 할 수 없습니다 (3개 이상 필요)."}
    fa, aln = out_dir / "members.fa", out_dir / "aligned.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in members.items()), encoding="utf-8")
    say(f"MAFFT 정렬 중 ({n}개)...")
    with open(aln, "w", encoding="utf-8") as fh:
        subprocess.run([find_mafft(), "--auto", "--thread", str(threads), str(fa)], stdout=fh, stderr=subprocess.DEVNULL, check=True)
    vv, t, message = find_varvamp(), threshold, ""
    while t >= min_threshold - 1e-9:
        say(f"varVAMP 설계 중 (consensus 임계값 {t:.1f})...")
        run_dir = out_dir / "varvamp"
        if run_dir.exists():
            shutil.rmtree(run_dir)
        proc = subprocess.run([vv, "qpcr", "-t", f"{t:.2f}", "-a", str(ambig), "-pa", str(probe_ambig), "-th", str(threads),
                               "-n", "50", str(aln), str(run_dir)], capture_output=True, text=True)
        if (run_dir / "qpcr_primers.tsv").exists():
            designs = _parse(run_dir)[:top]
            if designs:
                return {"status": "ok", "threshold": round(t, 2), "n": n, "designs": designs, "message": ""}
        err = [ln.strip() for ln in _clean(proc.stdout + proc.stderr).splitlines() if "ERROR" in ln]
        message = err[-1].replace("ERROR:", "").strip() if err else "varVAMP가 설계를 내지 못했습니다."
        t -= step
    return {"status": "none", "n": n, "designs": [], "message": f"임계값 {threshold:.1f}~{min_threshold:.1f} 어디서도 설계가 없습니다: {message}"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("segment")
    ap.add_argument("lineage", help="lineage name as in the metadata, or ALL")
    ap.add_argument("--threshold", type=float, default=0.9)
    ap.add_argument("--min-threshold", type=float, default=0.6)
    ap.add_argument("--ambig", type=int, default=2)
    ap.add_argument("--probe-ambig", type=int, default=1)
    ap.add_argument("--top", type=int, default=3)
    a = ap.parse_args()
    sys.path.insert(0, str(ROOT / "webapp"))
    import json  # noqa: PLC0415
    import scope_view  # noqa: PLC0415
    seg = a.segment.lower()
    prof = json.loads((ROOT / "profiles" / "lassa.json").read_text(encoding="utf-8"))
    v = scope_view.View(ROOT, prof.get("metadata") or {}, 1, 0)
    ids, lin = v.ids(seg, "human"), v.lineage_of(seg)
    seqs = read_fasta_seqs(ROOT / "data" / "filtered" / f"{seg}_filtered.fasta")
    members = {k: seqs[k] for k in ids if k in seqs and (a.lineage == "ALL" or lin.get(k) == a.lineage)}
    res = design_lineage(members, ROOT / "output" / "redesign" / seg / a.lineage, a.threshold, a.min_threshold,
                         ambig=a.ambig, probe_ambig=a.probe_ambig, top=a.top, progress=print)
    print(res["status"], res.get("threshold"), res["message"])
    for d in res["designs"]:
        print(d["name"], d["F"], d["P"], d["R"], "amplicon", d.get("amplicon_len"), "penalty", d.get("penalty"))


if __name__ == "__main__":
    main()
