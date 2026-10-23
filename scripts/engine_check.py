#!/usr/bin/env python3
"""Does the Go probe engine (go/probedesign) give exactly what scripts/probe_designer.py gives?

The Python result is the answer key. Both engines get the same alignment and parameters (rule-based "own" candidates only --
the Primer3 and MFEprimer steps stay in Python); the candidate lists must agree on sequence, strand, position, IUPAC count,
Tm/GC (1 decimal), coverage, source and order, and so must the header numbers. Timing of both is reported.

    python3 scripts/engine_check.py P001 [--top 10] [--max-iupac 2]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import probe_designer as pd  # noqa: E402

GO_DIR = ROOT / "go" / "probedesign"
GO_BIN = GO_DIR / "probedesign"
FIELDS = ("seq", "strand", "col_start", "col_end", "length", "n_iupac", "tm", "gc", "msa_cov0", "msa_cov1", "msa_n", "source")


def find_go():
    return shutil.which("go") or next((str(p) for p in (Path.home() / ".local/go/bin/go", Path("/usr/local/go/bin/go")) if p.exists()), None)


def ensure_built(name="probedesign"):
    """Path of the compiled Go engine `go/<name>` (built on demand when Go is installed), or None."""
    gdir = ROOT / "go" / name
    exe = gdir / name if name != "scan" else gdir / "scan"
    src = max((p.stat().st_mtime for p in gdir.glob("*.go")), default=0)
    if exe.exists() and exe.stat().st_mtime >= src:
        return str(exe)
    go = find_go()
    if not go or not gdir.exists():
        return str(exe) if exe.exists() else None
    r = subprocess.run([go, "build", "-o", str(exe), "."], cwd=gdir, capture_output=True, text=True, timeout=300)
    return str(exe) if r.returncode == 0 else None


def run_go(msa, params):
    exe = ensure_built()
    if not exe:
        raise RuntimeError("Go engine not available (install Go or build go/probedesign)")
    req = {"consensus": msa["consensus"], "f_col": msa["f_col"], "r_col": msa["r_col"], "dots": [r["dots"] for r in msa["rows"]],
           "min_len": params["min_len"], "max_len": params["max_len"], "tm_range": params["tm_range"], "gc_range": list(params["gc_range"]),
           "top": params["top"], "max_iupac": params["max_iupac"], "tm_above": params["tm_above"]}
    t = time.perf_counter()
    r = subprocess.run([exe], input=json.dumps(req), capture_output=True, text=True, timeout=600)
    dt = time.perf_counter() - t
    if r.returncode:
        raise RuntimeError("Go engine failed: " + (r.stderr.strip() or f"rc={r.returncode}"))
    return json.loads(r.stdout), dt


def compare(msa, top=10, min_len=18, max_len=30, max_iupac=2, gc_range=(30.0, 80.0), tm_range=None, tm_above=None, label=""):
    params = dict(top=top, min_len=min_len, max_len=max_len, max_iupac=max_iupac, gc_range=gc_range, tm_range=tm_range, tm_above=tm_above)
    t = time.perf_counter()
    py = pd.design("", min_len=min_len, max_len=max_len, tm_range=tm_range, gc_range=gc_range, top=top, max_iupac=max_iupac,
                   tm_above=tm_above, method="own", msa=msa)
    py_sec = time.perf_counter() - t
    go, go_sec = run_go(msa, params)
    checks = []

    def chk(name, a, b):
        checks.append({"name": name, "ok": a == b, "python": a, "go": b})

    chk("primer Tm", round(py["primer_tm"], 9), round(go["primer_tm"], 9))
    chk("probe Tm window", [round(x, 9) for x in py["tm_range"]], [round(x, 9) for x in go["tm_range"]])
    chk("search window", list(py["window"]), list(go["window"]))
    chk("candidates tested / rejected", [py["n_total"], py["n_rejected"]], [go["n_total"], go["n_rejected"]])
    chk("genomes", py["n_rows"], go["n_rows"])
    pc = [{k: c[k] for k in FIELDS} for c in py["candidates"]]
    gc = [{k: c[k] for k in FIELDS} for c in go["candidates"]]
    chk("number of picked candidates", len(pc), len(gc))
    for i, (a, b) in enumerate(zip(pc, gc), 1):
        diff = [k for k in FIELDS if a[k] != b[k]]
        checks.append({"name": f"candidate {i}", "ok": not diff, "python": a if diff else a["seq"], "go": b if diff else b["seq"], "diff": diff})
    return {"label": label, "ok": all(c["ok"] for c in checks), "checks": checks, "python_sec": round(py_sec, 3), "go_sec": round(go_sec, 3),
            "n_rows": py["n_rows"], "n_candidates": py["n_total"] + py["n_rejected"], "params": {k: params[k] for k in ("top", "min_len", "max_len", "max_iupac")}}


# ---- whole-genome scan engine (go/scan vs mfe_common._oligo_sites / scan_fr / scan_probe) ----
def _scan_go(req):
    exe = ensure_built("scan")
    if not exe:
        raise RuntimeError("Go scan engine not available (install Go or build go/scan)")
    t = time.perf_counter()
    r = subprocess.run([exe], input=json.dumps(req), capture_output=True, text=True, timeout=1800)
    dt = time.perf_counter() - t
    if r.returncode:
        raise RuntimeError("Go scan engine failed: " + (r.stderr.strip() or f"rc={r.returncode}"))
    return json.loads(r.stdout)["results"], dt


def _load_set_genomes(msa, limit=None):
    """(genome id -> sequence) of the population this set was analysed on: the filtered fasta that holds its alignment rows."""
    import mfe_common as mc
    ids = [r["genome_id"] for r in msa["rows"]]
    for seg in ("s", "l"):
        f = ROOT / "data" / "filtered" / f"{seg}_filtered.fasta"
        if f.exists():
            seqs = mc.read_fasta_seqs(f)
            if ids and ids[0] in seqs:
                return {i: seqs[i] for i in (ids[:limit] if limit else ids) if i in seqs}
    raise RuntimeError("population fasta of this set not found under data/filtered")


def _scan_oligos(msa):
    import mfe_common as mc
    cons = msa["consensus"]
    seg = lambda c: cons[c[0]:c[1] + 1].replace("-", "")
    f, r = seg(msa["f_col"]), mc.reverse_complement(seg(msa["r_col"]))
    p = seg(msa["p_col"]) if msa.get("p_col") else cons[len(cons) // 2 - 12:len(cons) // 2 + 12].replace("-", "")
    return f, r, p


def compare_scan(msa, three_prime_n=1, max_mm_fr=2, max_mm_probe=8, limit=None, label=""):
    """Same genomes, same oligos through mfe_common.scan_fr / scan_probe and through go/scan; every genome must agree."""
    import mfe_common as mc
    genomes = _load_set_genomes(msa, limit)
    f, r, p = _scan_oligos(msa)
    # a degenerate copy of the probe exercises the IUPAC path too
    p_deg = p[:5] + "R" + p[6:12] + "Y" + p[13:] if len(p) > 14 else p
    checks, py_total, go_total = [], 0.0, 0.0

    def one(name, py_fn, req, key):
        nonlocal py_total, go_total
        t = time.perf_counter()
        py = {gid: py_fn(gid, seq) for gid, seq in genomes.items()}
        py_total += time.perf_counter() - t
        res, dt = _scan_go(dict(req, genomes=[{"id": g, "seq": s, **({"amp_start": None, "amp_end": None})} for g, s in genomes.items()]))
        go_total += dt
        go = {x["id"]: x.get("hit") for x in res}
        bad = [g for g in genomes if json.loads(json.dumps(py[g])) != go.get(g)]
        hits = sum(1 for v in py.values() if v)
        checks.append({"name": name, "ok": not bad, "python": f"{hits}/{len(genomes)} genomes with a hit", "go": ("identical" if not bad else f"{len(bad)} differ, e.g. {bad[0]}: py={py[bad[0]]} go={go.get(bad[0])}")})

    one("F/R pair search (scan_fr)", lambda g, s: mc.scan_fr(f, r, s, three_prime_n, max_mm_fr),
        {"mode": "fr", "f": f, "r": r, "three_prime_n": three_prime_n, "max_mm": max_mm_fr, "max_size": 2000}, "fr")
    one("probe search (scan_probe)", lambda g, s: mc.scan_probe(p, s, None, None, max_mm_probe),
        {"mode": "probe", "probe": p, "max_mm": max_mm_probe}, "probe")
    one("probe with IUPAC codes", lambda g, s: mc.scan_probe(p_deg, s, None, None, max_mm_probe),
        {"mode": "probe", "probe": p_deg, "max_mm": max_mm_probe}, "probe")
    one("F/R search, 0 mismatches", lambda g, s: mc.scan_fr(f, r, s, 0, 0),
        {"mode": "fr", "f": f, "r": r, "three_prime_n": 0, "max_mm": 0, "max_size": 2000}, "fr")
    return {"label": label, "ok": all(c["ok"] for c in checks), "checks": checks, "python_sec": round(py_total, 3), "go_sec": round(go_total, 3),
            "n_rows": len(genomes), "n_candidates": len(genomes) * len(checks), "kind": "scan"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("set_id")
    ap.add_argument("--scan", action="store_true", help="compare the whole-genome scan engine instead of probe design")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--max-iupac", type=int, default=2)
    a = ap.parse_args()
    res = compare_scan(pd.load_msa(a.set_id), label=a.set_id) if a.scan else compare(pd.load_msa(a.set_id), top=a.top, max_iupac=a.max_iupac, label=a.set_id)
    for c in res["checks"]:
        print("PASS" if c["ok"] else "FAIL", c["name"], "" if c["ok"] else f"python={c['python']} go={c['go']}")
    print(f"python {res['python_sec']}s  go {res['go_sec']}s  -> {'IDENTICAL' if res['ok'] else 'DIFFERENT'}")
    sys.exit(0 if res["ok"] else 1)


if __name__ == "__main__":
    main()
