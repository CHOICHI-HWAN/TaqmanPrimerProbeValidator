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


def ensure_built():
    """Path of the compiled engine (built on demand when Go is installed), or None."""
    src = max(p.stat().st_mtime for p in GO_DIR.glob("*.go")) if GO_DIR.exists() else 0
    if GO_BIN.exists() and GO_BIN.stat().st_mtime >= src:
        return str(GO_BIN)
    go = find_go()
    if not go or not GO_DIR.exists():
        return str(GO_BIN) if GO_BIN.exists() else None
    r = subprocess.run([go, "build", "-o", str(GO_BIN), "."], cwd=GO_DIR, capture_output=True, text=True, timeout=300)
    return str(GO_BIN) if r.returncode == 0 else None


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("set_id")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--max-iupac", type=int, default=2)
    a = ap.parse_args()
    res = compare(pd.load_msa(a.set_id), top=a.top, max_iupac=a.max_iupac, label=a.set_id)
    for c in res["checks"]:
        print("PASS" if c["ok"] else "FAIL", c["name"], "" if c["ok"] else f"python={c['python']} go={c['go']}")
    print(f"python {res['python_sec']}s  go {res['go_sec']}s  -> {'IDENTICAL' if res['ok'] else 'DIFFERENT'}")
    sys.exit(0 if res["ok"] else 1)


if __name__ == "__main__":
    main()
