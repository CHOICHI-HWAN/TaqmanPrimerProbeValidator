#!/usr/bin/env python3
"""Primer-pair design from a lineage's genomes with PrimerHunter (Duitama et al., NAR 2009; GPL, built under tools/primerhunter).

PrimerHunter picks forward/reverse primers whose predicted melting temperature against the target genomes is high enough
(mismatch-tolerant nearest-neighbour model), pairs them inside a product-length window and, when no pair covers all targets, picks a
small set of pairs that together cover them (greedy set cover). It designs primers only -- the probe is added by the caller.

PrimerHunter is slow (every candidate is aligned against every target), so this module
  * runs it on a small, maximally diverse sample of the lineage (`pick_representatives`; the full lineage is validated afterwards), and
  * splits the source sequence into position windows and runs one PrimerHunter process per window in parallel.

    python3 scripts/redesign_primerhunter.py output/redesign/s/LIV_V/aligned.fa --n 30 --threads 12
"""
import argparse
import hashlib
import math
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PH_BIN = ROOT / "tools" / "primerhunter" / "PrimerHunter" / "primerhunter"


def find_primerhunter() -> str:
    if not PH_BIN.exists():
        sys.exit(f"[error] PrimerHunter is not built: {PH_BIN} (see tools/primerhunter/README.md)")
    return str(PH_BIN)


def read_fasta(path: Path) -> dict:
    out, key = {}, None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            out[key] = []
        elif key:
            out[key].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def pick_representatives(aligned: dict, n: int) -> list:
    """n accessions that cover the lineage's diversity: the medoid first (it is PrimerHunter's coordinate system), then repeatedly the
    genome farthest (Hamming distance on the alignment, gaps/N ignored) from everything already chosen."""
    import numpy as np  # noqa: PLC0415
    ids = list(aligned)
    if len(ids) <= n:
        return ids
    arr = np.frombuffer("".join(aligned[i].upper() for i in ids).encode(), dtype=np.uint8).reshape(len(ids), -1)
    valid = (arr != ord("-")) & (arr != ord("N"))

    def dist_to(k):                                       # distance of every genome to genome k
        both = valid & valid[k]
        return ((arr != arr[k]) & both).sum(1) / np.maximum(both.sum(1), 1)

    rng = np.random.default_rng(0)
    sample = rng.choice(len(ids), size=min(len(ids), 60), replace=False)
    first = int(min(sample, key=lambda k: dist_to(k)[sample].sum()))
    chosen, mind = [first], dist_to(first)
    while len(chosen) < n:
        k = int(mind.argmax())
        chosen.append(k)
        mind = np.minimum(mind, dist_to(k))
    return [ids[k] for k in chosen]


def parse_pairs(text: str) -> list:
    """PrimerHunter pair lines: 'Fwd seq pos minTm avgTm maxTm nonTm cov% nonCov% - Rev seq pos ... - pairCov% pairNonCov%'."""
    pairs = []
    for line in text.splitlines():
        t = line.split()
        if len(t) != 22 or not t[0].startswith("Forward_") or t[10][:8] != "Reverse_":
            continue
        try:
            fpos, rpos = int(t[2]), int(t[12])
            pairs.append({
                "F": t[1].upper(), "R": t[11].upper(), "F_pos": fpos, "R_pos": rpos, "amplicon_len": abs(rpos - fpos) + 1,
                "F_tm": float(t[4]), "R_tm": float(t[14]), "F_cov": float(t[7].rstrip("%")), "R_cov": float(t[17].rstrip("%")),
                "pair_cov": float(t[20].rstrip("%"))})
        except ValueError:
            continue
    return pairs


def design_pairs(targets: dict, out_dir: Path, threads=12, min_prod=70, max_prod=200, min_len=18, max_len=25,
                 min_tm=50.0, max_tm=65.0, min_cov=80.0, gc=(30, 70), timeout=3600, progress=None, top=3) -> dict:
    """targets = {accession: sequence}, the first one is the coordinate source. -> {"status", "pairs": [...], "message"}"""
    def say(msg):
        if progress:
            progress(msg)

    out_dir.mkdir(parents=True, exist_ok=True)
    ph = find_primerhunter()
    tf = out_dir / "targets.fa"
    tf.write_text("".join(f">{k}\n{v}\n" for k, v in targets.items()), encoding="utf-8")
    first = next(iter(targets.values()))
    size = max(60, math.ceil(len(first) / max(1, threads)))
    procs = []
    # an earlier run with the same targets and settings is reused (PrimerHunter takes minutes per lineage)
    key = f"{threads}|{min_prod}|{max_prod}|{min_len}|{max_len}|{min_tm}|{max_tm}|{min_cov}|{gc}|" + hashlib.sha1(
        "".join(f"{k}{v}" for k, v in targets.items()).encode()).hexdigest()
    stamp = out_dir / "run_key.txt"
    cached = sorted(out_dir.glob("chunk_*.txt"))
    if (cached and stamp.exists() and stamp.read_text() == key
            and all("Statistics" in c.read_text(encoding="utf-8", errors="replace") for c in cached)):
        say("이전 PrimerHunter 결과를 재사용합니다")
        return _finish([(None, c) for c in cached], targets, top)
    stamp.write_text(key, encoding="utf-8")
    for old in cached:
        old.unlink()
    for i, b in enumerate(range(0, len(first), size)):
        outp = out_dir / f"chunk_{i:02d}.txt"
        cmd = [ph, "-tf", str(tf), "-minCoverageTargets", f"{min_cov:g}", "-minProdLength", str(min_prod),
               "-maxProdLength", str(max_prod), "-minPrimerLength", str(min_len), "-maxPrimerLength", str(max_len),
               "-minTempTargets", f"{min_tm:g}", "-maxTempTargets", f"{max_tm:g}", "-saltConc", "0.05",
               "-primersConc", "0.00000025", "-minGCContent", str(gc[0]), "-maxGCContent", str(gc[1]),
               "-beginPosForward", str(b), "-endPosForward", str(b + size),
               "-beginPosReverse", str(b + min_prod), "-endPosReverse", str(b + size + max_prod)]
        procs.append((subprocess.Popen(cmd, stdout=open(outp, "w"), stderr=subprocess.DEVNULL), outp))
    t0, n_total = time.time(), len(procs)
    while True:
        running = [p for p, _ in procs if p.poll() is None]
        say(f"PrimerHunter 실행 중: 구간 {n_total - len(running)}/{n_total} 완료 ({int(time.time() - t0)}초)")
        if not running:
            break
        if time.time() - t0 > timeout:
            for p in running:
                p.kill()
            say("PrimerHunter 시간 초과 — 끝난 구간의 결과만 사용합니다")
            break
        time.sleep(5)
    return _finish(procs, targets, top)


def _finish(procs, targets, top):
    pairs = []
    for _, outp in procs:
        pairs += parse_pairs(outp.read_text(encoding="utf-8", errors="replace"))
    if not pairs:
        return {"status": "none", "pairs": [], "n_targets": len(targets),
                "message": "PrimerHunter가 조건(Tm·길이·커버리지)을 만족하는 프라이머 쌍을 찾지 못했습니다."}
    # best first: most targets amplified by the pair, then the two Tm's closest to each other and to ~58 degC
    pairs.sort(key=lambda p: (-p["pair_cov"], abs(p["F_tm"] - p["R_tm"]) + abs((p["F_tm"] + p["R_tm"]) / 2 - 58.0)))
    picked = []
    for p in pairs:
        if all(abs(p["F_pos"] - q["F_pos"]) > 15 for q in picked):          # distinct locations, not shifted copies of one pair
            picked.append(p)
        if len(picked) >= top:
            break
    for i, p in enumerate(picked):
        p["name"] = f"PH_{i}"
    return {"status": "ok", "pairs": picked, "n_targets": len(targets), "n_pairs_total": len(pairs), "message": ""}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("aligned_fasta")
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--out", default=str(ROOT / "output" / "redesign" / "_ph_test"))
    a = ap.parse_args()
    aligned = read_fasta(Path(a.aligned_fasta))
    reps = pick_representatives(aligned, a.n)
    raw = {k: aligned[k].replace("-", "") for k in reps}
    res = design_pairs(raw, Path(a.out), threads=a.threads, progress=print)
    print(res["status"], res.get("n_pairs_total"), res["message"])
    for p in res["pairs"]:
        print(p)


if __name__ == "__main__":
    main()
