#!/usr/bin/env python3
"""Interactive local web dashboard for the LASV primer/probe validation pipeline.

Serves the 14 built-in assay sets (reading the already-computed Snakemake pipeline
outputs) and lets the user submit new F/R/Probe sequences for live on-demand analysis
against the same pre-built population / standard-material / reference-genome databases.

Run:
    conda activate nextstrain   # needs only python3 stdlib + pyyaml (for config.yaml)
    python3 webapp/server.py [--port 8765]
"""
import argparse
import base64
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import traceback
import urllib.request
import uuid
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

def _detect_root() -> Path:
    """Project root (holds config/, data/, output/) -- normally the current working
    directory (this is how the webapp has always been invoked: `cd pipeline_snakemake &&
    python3 webapp/server.py`). When installed as a package and launched via the
    `taqmanvalidator-webapp` console script from a scaffolded project directory, cwd is
    still the right answer for the same reason. Falls back to the source tree's own
    location only if cwd doesn't look like a project (e.g. invoked from an unrelated
    directory by mistake), so the error the user sees is "no such file" from a sane
    default rather than a confusing one from an arbitrary cwd."""
    cwd = Path.cwd()
    if (cwd / "config" / "config.yaml").exists():
        return cwd
    return Path(__file__).resolve().parent.parent


ROOT = _detect_root()  # user's project directory: config/, data/, output/, logs/
AUSPICE_CACHE_DIR = ROOT / "data" / "auspice"


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # process exists, we just don't own it
    except OSError:
        return False
    return True


def fetch_auspice(seg: str) -> bytes:
    """Nextstrain's public auspice v2 tree JSON for this segment (data.nextstrain.org has
    no CORS headers, so the browser can't fetch it directly -- proxy + cache it here)."""
    AUSPICE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = AUSPICE_CACHE_DIR / f"{seg}.json"
    if cache_path.exists():
        return cache_path.read_bytes()
    url = f"https://data.nextstrain.org/lassa_{seg}.json"
    req = urllib.request.Request(url, headers={"Accept-Encoding": "gzip"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        raw = resp.read()
        if resp.headers.get("Content-Encoding") == "gzip":
            raw = gzip.decompress(raw)
    cache_path.write_bytes(raw)
    return raw


def _fetch_state_path():
    return ROOT / "config" / "fetch_state.json"


def load_fetch_state():
    path = _fetch_state_path()
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    except (json.JSONDecodeError, OSError):
        return {}


def save_fetch_state(state):
    path = _fetch_state_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")


_PROFILE_CACHE = {}


def _profile_path():
    from virus_profile import find_profile_path  # noqa: PLC0415 -- scripts/ is added to sys.path above
    return find_profile_path(ROOT)


def get_profile():
    """The validated virus profile (reloaded when the file changes)."""
    from virus_profile import load_profile  # noqa: PLC0415 -- scripts/ is added to sys.path above
    path = _profile_path()
    stamp = (str(path), path.stat().st_mtime_ns)
    if _PROFILE_CACHE.get("stamp") != stamp:
        _PROFILE_CACHE.update(stamp=stamp, profile=load_profile(path))
    return _PROFILE_CACHE["profile"]


def pooled_label():
    return get_profile()["population"]["pooled_label"]


def scope_ctx():
    """The host-scope view (human / animal / unknown / all) of the current request: see scope_view.py. The pass/fail rule and the
    mismatch thresholds are the pipeline's own (scripts/host_scope.py, config.yaml)."""
    text = _config_yaml_path().read_text(encoding="utf-8")
    vals = {}
    for key in ("primer_mismatch_pass", "probe_mismatch"):
        m = re.search(rf"^{key}:\s*(\d+)", text, re.MULTILINE)
        vals[key] = int(m.group(1)) if m else PARAM_DEFAULTS[key]
    return scope_view.View(ROOT, get_profile().get("metadata") or {}, vals["primer_mismatch_pass"], vals["probe_mismatch"])


def host_scope_info():
    """Per segment: the active scope, the primary scope and the genome counts of each scope (what the toggle shows)."""
    v = scope_ctx()
    return {seg.upper(): v.info(seg) for seg in seg_keys()}


def seg_keys():
    """Segment keys of the active virus profile (lowercase, e.g. ["s", "l"])."""
    return list(get_profile()["segments"])


def probe_nextstrain_files():
    """HEAD-checks every raw file `fetch_raw` would download, without downloading it --
    Nextstrain's S3 bucket publishes each object's content sha256 as the custom
    `x-amz-meta-sha256sum` header (falls back to the S3 ETag, also content-derived for
    these small non-multipart uploads, if that's ever missing), so this settles "did the
    remote dataset actually change since our last fetch" for ~1KB of HEAD response instead
    of pulling the ~400KB/~40KB .zst files themselves. Raises on a network/HTTP error --
    callers decide how to handle "couldn't check" (fail open: treat as changed, never block
    a real update on a flaky check)."""
    prof = get_profile()
    base_url = prof["source"]["base_url"].rstrip("/")
    out = {}
    for seg in prof["segments"]:
        out[seg] = {}
        for kind, fname in (("sequences", "sequences.fasta.zst"), ("metadata", "metadata.tsv.zst")):
            url = f"{base_url}/{seg}/{fname}"
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=30) as resp:
                h = resp.headers
                out[seg][kind] = {
                    "sha256": h.get("x-amz-meta-sha256sum") or "",
                    "etag": (h.get("ETag") or "").strip('"'),
                    "record_count": h.get("x-amz-meta-recordcount") or "",
                    "last_modified": h.get("Last-Modified") or "",
                }
    return out


def _fetch_fingerprint(entry):
    """One comparable value for a probed file -- the source's own content sha256 when
    present (exact), else the S3 ETag (still content-derived for these small
    non-multipart uploads)."""
    return (entry or {}).get("sha256") or (entry or {}).get("etag") or ""


def save_fetch_baseline_from_remote():
    """Records what the remote just reported as the new freshness baseline -- called right
    after a successful fetch_raw run, so the next 데이터 업데이트 click compares against what
    was actually just downloaded, not a stale earlier check. Best-effort: a probe failure
    here must never fail the pipeline job that just succeeded."""
    try:
        save_fetch_state({"checked_at": datetime.now().isoformat(timespec="seconds"),
                           "remote": probe_nextstrain_files()})
    except Exception:  # noqa: BLE001
        traceback.print_exc()


STATIC_DIR = Path(__file__).resolve().parent / "static"  # this package's own bundled assets
sys.path.insert(0, str(ROOT / "scripts"))
import scope_view  # noqa: E402 -- host-scope view of the results (webapp/scope_view.py)
import literature  # noqa: E402 -- literature primer lists (webapp/literature.py)

from live_analyze import analyze_set  # noqa: E402

JOBS = {}
_FRESHNESS_CACHE = {}
JOBS_LOCK = threading.Lock()


def find_snakemake() -> str:
    from mfe_common import env_tool  # noqa: PLC0415 -- scripts/ is added to sys.path above
    configured = env_tool(ROOT, "snakemake")            # the conda env that setup.sh recorded (config.yaml: env_prefix)
    if configured:
        return configured
    found = shutil.which("snakemake")
    if found:
        return found
    # lasv_pipeline is this project's own dedicated env (snakemake + mafft, nothing borrowed
    # from an unrelated tool's env); the rest are only a fallback for older checkouts.
    for cand in ("lasv_pipeline", "hybracter_0.11.2", "squirrel", "comparem2", "verkko"):
        p = Path(f"/BI/conda/envs/{cand}/bin/snakemake")
        if p.exists():
            return str(p)
    sys.exit("[error] snakemake not found. Run ./setup.sh (creates the conda env), or put snakemake on PATH.")


SNAKEMAKE_BIN = None  # resolved lazily on first use (avoids a slow glob at import time)


def run_pipeline_job(job_id, args, label, on_success=None):
    # This runs in a daemon thread with no caller to propagate an exception to, and
    # _handle_pipeline refuses to start a *new* job while any job's status is still
    # "running" -- so if anything here raises before the final status write (snakemake
    # binary missing, Popen failing, disk error, ...) the job is stuck "running" forever
    # and blocks every future 데이터 업데이트/분석 다시하기 click. Every exit path below must
    # reach the final JOBS[job_id]["status"] write.
    log_path = ROOT / "logs" / f"webapp_{job_id}.log"
    try:
        global SNAKEMAKE_BIN
        if SNAKEMAKE_BIN is None:
            SNAKEMAKE_BIN = find_snakemake()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = [SNAKEMAKE_BIN, "--cores", "16", "-p", "--rerun-incomplete", "--latency-wait", "10"] + args
        with JOBS_LOCK:
            JOBS[job_id]["cmd"] = " ".join(cmd)
        with open(log_path, "w", encoding="utf-8") as logf:
            child_env = dict(os.environ, PATH=str(Path(SNAKEMAKE_BIN).parent) + os.pathsep + os.environ.get("PATH", ""))
            proc = subprocess.Popen(cmd, cwd=str(ROOT), stdout=logf, stderr=subprocess.STDOUT, text=True, env=child_env)
            with JOBS_LOCK:
                JOBS[job_id]["pid"] = proc.pid
                JOBS[job_id]["log_path"] = str(log_path)
            proc.wait()
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "done" if proc.returncode == 0 else "error"
            JOBS[job_id]["returncode"] = proc.returncode
            JOBS[job_id]["progress"] = tail_progress(log_path)
        if proc.returncode == 0 and on_success is not None:
            on_success()
    except Exception as exc:  # noqa: BLE001
        traceback.print_exc()
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error"
            JOBS[job_id]["returncode"] = None
            JOBS[job_id]["error"] = str(exc)
            JOBS[job_id]["progress"] = f"실행 실패: {exc}"
            JOBS[job_id]["log_path"] = str(log_path)


def run_report_job(job_id, base_url):
    # generate_report.py only needs stdlib + a Chrome/Chromium binary on PATH (its own
    # find_chrome() tries several common names, not just "google-chrome" -- no
    # snakemake/mfeprimer/pyyaml), so it runs under this same server's own Python
    # (sys.executable) rather than needing find_snakemake()'s conda-env lookup.
    log_path = ROOT / "logs" / f"webapp_report_{job_id}.log"
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        script = ROOT / "scripts" / "generate_report.py"
        out_pdf = ROOT / "output" / "report" / "report.pdf"
        cmd = [sys.executable, str(script), "--base-url", base_url, "--out", str(out_pdf)]
        with JOBS_LOCK:
            JOBS[job_id]["cmd"] = " ".join(cmd)
        with open(log_path, "w", encoding="utf-8") as logf:
            proc = subprocess.Popen(cmd, cwd=str(ROOT), stdout=logf, stderr=subprocess.STDOUT, text=True)
            with JOBS_LOCK:
                JOBS[job_id]["pid"] = proc.pid
                JOBS[job_id]["log_path"] = str(log_path)
            proc.wait()
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "done" if proc.returncode == 0 and out_pdf.exists() else "error"
            JOBS[job_id]["returncode"] = proc.returncode
            JOBS[job_id]["progress"] = tail_progress(log_path)
            if JOBS[job_id]["status"] == "done":
                JOBS[job_id]["download_url"] = "/api/report/download"
    except Exception as exc:  # noqa: BLE001
        traceback.print_exc()
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error"
            JOBS[job_id]["returncode"] = None
            JOBS[job_id]["error"] = str(exc)
            JOBS[job_id]["progress"] = f"실행 실패: {exc}"
            JOBS[job_id]["log_path"] = str(log_path)


_PD_SEQ_RE = re.compile(r"^[ACGTRYSWKMBDHVN]{10,60}$")
_PD_ID_RE = re.compile(r"^[A-Za-z0-9_\-]+$")


def _pd_slug(seq):
    import hashlib  # noqa: PLC0415
    return hashlib.sha1(seq.encode()).hexdigest()[:10]


def pd_scope(set_id, host):
    """(scope, ids) for a probe-design request: the same host scope the dashboard toggle uses (human by default; anything
    unknown means the primary scope), and the accessions of that scope. Background jobs get it passed in because they
    outlive the request."""
    with open(ROOT / "input_primers" / f"{set_id}.tsv", encoding="utf-8") as fh:
        seg = next(csv.DictReader(fh, delimiter="\t"))["segment"].lower()
    v = scope_ctx()
    scope = host if host in scope_view.H.SCOPES else v.primary(seg)
    return scope, v.ids(seg, scope)


def validate_oligos(segment, f_seq, r_seq, p_seq, out_dir, name, progress=None, scope=None, ids=None):
    """Runs the real validation (live_analyze.analyze_set: the same mfeprimer searches, MSA, KDCA standard check and
    dimer/hairpin check as a registered set) for one F/R/Probe triple, counting only the genomes of the host scope.
    result.json (everything, for the MSA modal) and summary.json (the numbers) go to out_dir. -> the summary dict"""
    params = load_params()
    result = analyze_set(
        name=name, segment=segment, f_id=f"{name}_F", f_seq=f_seq, r_id=f"{name}_R", r_seq=r_seq, p_id=f"{name}_P", p_seq=p_seq,
        root=ROOT, out_dir=out_dir, primer_mismatch_pass=params["primer_mismatch_pass"],
        probe_mismatch=params["probe_mismatch"], three_prime_n=params["three_prime_consecutive_match"],
        cpu=params["mfeprimer_cpu"], progress=progress)
    pass_mm, probe_mm = params["primer_mismatch_pass"], params["probe_mismatch"]
    genomes = result.get("genomes") or []
    if ids is not None:                                    # host scope: count only the genomes of the active scope
        genomes = [g for g in genomes if g.get("genome_id") in ids]
        scoped = scope_ctx().summary_from_rows(segment.lower(), {g["genome_id"]: g for g in genomes}, scope)
        result["summary"] = {**(result.get("summary") or {}), **scoped}

    def num(v):
        try:
            return int(v)
        except (TypeError, ValueError):
            return None
    lineages, validated_ids = {}, []
    for g in genomes:
        f, p, r = num(g.get("f_mismatch")), num(g.get("p_mismatch")), num(g.get("r_mismatch"))
        ok = (f is not None and f <= pass_mm and r is not None and r <= pass_mm and p is not None and p <= probe_mm
              and g.get("f_3p_ok") is not False and g.get("r_3p_ok") is not False)
        rec = lineages.setdefault(g.get("lineage") or "NA", [0, 0])
        rec[0] += 1
        rec[1] += 1 if ok else 0
        if ok:
            validated_ids.append(g["genome_id"])
    summary = {"validated_ids": validated_ids,
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"), "host_scope": scope, "summary": result.get("summary"),
        "standard": result.get("standard"), "dimer_hairpin": result.get("dimer_hairpin"),
        "probe_mismatch": probe_mm, "primer_mismatch_pass": pass_mm,
        "lineages": [{"lineage": k, "total": v[0], "validated": v[1]} for k, v in sorted(lineages.items())],
    }
    (out_dir / "result.json").write_text(json.dumps(result), encoding="utf-8")
    (out_dir / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    return summary


def probe_validate_one(set_id, seq, progress=None, scope=None, ids=None):
    """The set's registered F/R plus one designed probe. The result lives ONLY under output/probe_design/<set>/validation/<slug>/
    -- it is not a set, so it never reaches the overview, summaries or reports."""
    with open(ROOT / "input_primers" / f"{set_id}.tsv", encoding="utf-8") as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    out_dir = ROOT / "output" / "probe_design" / set_id / "validation" / _pd_slug(seq)
    summary = validate_oligos(row["segment"], row["forward_primer"], row["reverse_primer"], seq, out_dir, f"{set_id}_PD",
                              progress, scope, ids)
    summary.update({"set_id": set_id, "seq": seq, "slug": _pd_slug(seq)})
    (out_dir / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    return summary


def run_probe_validate_job(job_id, set_id, seq, scope=None, ids=None):
    def progress(msg):
        with JOBS_LOCK:
            JOBS[job_id]["progress"] = msg

    try:
        probe_validate_one(set_id, seq, progress, scope, ids)
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "done"
            JOBS[job_id]["progress"] = "완료"
    except Exception as exc:  # noqa: BLE001
        traceback.print_exc()
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error"
            JOBS[job_id]["error"] = str(exc)
            JOBS[job_id]["progress"] = f"실행 실패: {exc}"


def _pd_rank_key(v):
    sm, std = v.get("summary") or {}, v.get("standard") or {}
    return (-(sm.get("validated") or 0), 0 if std.get("ok") else 1, -(sm.get("probe_coverage_mm0_pct") or 0),
            sm.get("probe_dropout") or 0, v.get("n_iupac", 0))


def run_probe_best_job(job_id, set_id, opts, scope=None, ids=None):
    """Design candidates, run the real validation for EVERY one of them, and rank the results: most validated genomes first,
    then standard-material match, Probe coverage, fewer dropouts, fewer IUPAC positions. Everything stays in the probe-design card."""
    import probe_designer  # noqa: PLC0415 -- scripts/ is on sys.path

    def progress(msg):
        with JOBS_LOCK:
            JOBS[job_id]["progress"] = msg

    try:
        progress("후보 설계 중...")
        design = probe_designer.run_design(set_id, ids=ids, scope=scope, **opts)
        cands = design["candidates"]
        ranking = []
        for i, c in enumerate(cands, 1):
            prefix = f"{i}/{len(cands)} 후보 검증 중 ({c['seq']}): "
            v = probe_validate_one(set_id, c["seq"], lambda m, prefix=prefix: progress(prefix + m), scope, ids)
            v["n_iupac"], v["tm"], v["source"] = c.get("n_iupac", 0), c.get("tm"), c.get("source", "own")
            ranking.append(v)
        ranking.sort(key=_pd_rank_key)
        best = {
            "set_id": set_id, "generated": time.strftime("%Y-%m-%d %H:%M:%S"), "host_scope": scope, "n_tested": len(ranking),
            "options": {k: opts[k] for k in ("top", "max_iupac", "probe_mismatch", "method") if k in opts},
            "tm_above": design.get("tm_above"),
            "ranking": [{"seq": v["seq"], "tm": v.get("tm"), "source": v.get("source", "own"), "n_iupac": v.get("n_iupac", 0),
                         "validated": (v.get("summary") or {}).get("validated"),
                         "validated_pct": (v.get("summary") or {}).get("validated_pct"),
                         "probe_pct": (v.get("summary") or {}).get("probe_coverage_mm0_pct"),
                         "dropout": (v.get("summary") or {}).get("probe_dropout"),
                         "std_ok": bool((v.get("standard") or {}).get("ok"))} for v in ranking],
        }
        (ROOT / "output" / "probe_design" / set_id / "best.json").write_text(json.dumps(best), encoding="utf-8")
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "done"
            JOBS[job_id]["progress"] = "완료"
    except Exception as exc:  # noqa: BLE001
        traceback.print_exc()
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error"
            JOBS[job_id]["error"] = str(exc)
            JOBS[job_id]["progress"] = f"실행 실패: {exc}"


def _ph_designs(rph, members, out_dir, opts, seg, ids, scope, progress):
    """PrimerHunter F/R pairs for one lineage (on a diverse sample of it), each completed with a designed probe."""
    out_dir.mkdir(parents=True, exist_ok=True)
    n_sample = opts.get("ph_sample", 25)
    if len(members) > n_sample:
        progress(f"MAFFT 정렬로 대표 게놈 {n_sample}개 고르는 중 ({len(members)}개 중)...")
        from mfe_common import find_mafft  # noqa: PLC0415
        fa, aln = out_dir / "members.fa", out_dir / "aligned.fa"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in members.items()), encoding="utf-8")
        with open(aln, "w", encoding="utf-8") as fh:
            subprocess.run([find_mafft(), "--auto", "--thread", "4", str(fa)], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        reps = rph.pick_representatives(rph.read_fasta(aln), n_sample)
        targets = {k: members[k] for k in reps}
    else:
        targets = dict(members)
    res = rph.design_pairs(targets, out_dir / "primerhunter", threads=opts.get("ph_threads", 12), top=opts["top"],
                           min_cov=opts.get("ph_min_cov", 80.0), min_tm=opts.get("ph_min_tm", 56.0), max_tm=opts.get("ph_min_tm", 56.0) + 14, progress=progress)
    designs = []
    for i, pr in enumerate(res["pairs"], 1):
        progress(f"프라이머 쌍 {i}/{len(res['pairs'])}: Probe 설계 중...")
        probe = probe_for_pair(seg, pr["F"], pr["R"], ids, scope, out_dir / "probe_work" / pr["name"])
        if not probe:
            continue
        designs.append({**pr, "P": probe["seq"], "P_tm": probe["tm"], "P_source": probe["source"], "P_iupac": probe["n_iupac"], "P_note": probe["note"]})
    status = "ok" if designs else ("none" if res["status"] == "ok" else res["status"])
    msg = res["message"] or ("" if designs else "프라이머 쌍은 나왔지만 이 쌍들 사이에서 조건에 맞는 Probe를 찾지 못했습니다.")
    return {"status": status, "designs": designs, "message": msg, "sample_n": len(targets), "threshold": None}


def redesign_scope(seg, host):
    v = scope_ctx()
    scope = host if host in scope_view.H.SCOPES else v.primary(seg)
    return scope, v.ids(seg, scope)


def _rd_dir(seg, scope, engine="varvamp"):
    base = ROOT / "output" / "redesign" / seg / scope
    return base / "primerhunter" if engine == "primerhunter" else base


def probe_for_pair(seg, f_seq, r_seq, ids, scope, work_dir):
    """PrimerHunter designs primers only. For a F/R pair: run the pipeline once with a placeholder probe to get the amplicon alignment,
    then let the probe designer (own rules + Primer3, IUPAC allowed) pick the probe that matches the most genomes of the scope.
    -> {"seq", "tm", "source", "n_iupac", "msa_cov0"} or None"""
    import probe_designer  # noqa: PLC0415 -- scripts/ is on sys.path
    params = load_params()
    res = analyze_set(
        name="ph_probe", segment=seg, f_id="F", f_seq=f_seq, r_id="R", r_seq=r_seq, p_id="P", p_seq=f_seq, root=ROOT,
        out_dir=work_dir, primer_mismatch_pass=params["primer_mismatch_pass"], probe_mismatch=params["probe_mismatch"],
        three_prime_n=params["three_prime_consecutive_match"], cpu=params["mfeprimer_cpu"])
    msa = res.get("msa")
    if not msa or not msa.get("rows"):
        return None
    # PrimerHunter's primers are short-amplicon, relatively high-Tm: ask for a probe above the primers first, then accept one a bit
    # below them (flagged) rather than giving no probe at all.
    for above, note in (((0.0, 12.0), ""), ((-3.0, 15.0), "Probe Tm이 프라이머와 비슷하거나 약간 낮음 (조건을 완화해 설계)")):
        try:
            d = probe_designer.design("", top=6, max_iupac=2, method="both", ids=ids, msa=msa, tm_above=above)
        except SystemExit:
            return None
        if d["candidates"]:
            c = d["candidates"][0]
            return {"seq": c["seq"], "tm": c.get("tm"), "source": c.get("source"), "n_iupac": c.get("n_iupac", 0),
                    "msa_cov0": c.get("msa_cov0"), "note": note, "primer_tm": round(d["primer_tm"], 1)}
    return None


def _set_cover(ids_by_design, targets):
    """Greedy set cover over the designs of one lineage: which designs together pass (real validation) the most genomes of
    `targets`, in the order they add genomes. -> {"steps": [{"name", "added", "covered", "pct"}], "covered", "total", "pct"}"""
    left, steps, covered = set(targets), [], set()
    pool = {k: v & targets for k, v in ids_by_design.items()}
    while pool:
        name, best = max(pool.items(), key=lambda kv: len(kv[1] & left))
        gain = best & left
        if not gain:
            break
        covered |= gain
        left -= gain
        steps.append({"name": name, "added": len(gain), "covered": len(covered),
                      "pct": round(100.0 * len(covered) / len(targets), 1) if targets else 0.0})
        del pool[name]
    total = len(targets)
    return {"steps": steps, "covered": len(covered), "total": total, "pct": round(100.0 * len(covered) / total, 1) if total else 0.0}


def run_redesign_job(job_id, seg, scope, ids, opts):
    """varVAMP redesign per lineage (and optionally for all genomes together), every design validated with the real pipeline
    validation. Results: output/redesign/<seg>/<scope>/run.json (+ one folder per lineage/design) -- shown only in the redesign page."""
    import redesign_varvamp  # noqa: PLC0415 -- scripts/ is on sys.path
    import redesign_primerhunter  # noqa: PLC0415
    from mfe_common import read_fasta_seqs  # noqa: PLC0415
    engine = opts.get("engine", "varvamp")

    def progress(msg):
        with JOBS_LOCK:
            JOBS[job_id]["progress"] = msg

    try:
        v = scope_ctx()
        lin = v.lineage_of(seg)
        seqs = read_fasta_seqs(ROOT / "data" / "filtered" / f"{seg}_filtered.fasta")
        counts = {}
        for i in ids:
            counts[lin.get(i, "NA")] = counts.get(lin.get(i, "NA"), 0) + 1
        wanted = opts.get("lineages") or sorted(counts)
        groups = [(L, sorted(i for i in ids if lin.get(i, "NA") == L)) for L in wanted if L in counts]       # sorted: reproducible sampling
        if opts.get("include_all"):
            groups.append(("ALL", sorted(ids)))
        base = _rd_dir(seg, scope, engine)
        base.mkdir(parents=True, exist_ok=True)
        run = {"engine": engine, "segment": seg.upper(), "host_scope": scope, "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
               "population_n": len(ids), "params": {k: opts[k] for k in ("threshold", "min_threshold", "ambig", "probe_ambig", "top", "min_n", "ph_sample", "ph_min_cov", "ph_min_tm") if k in opts},
               "lineages": [], "done": False}
        (base / "run.json").write_text(json.dumps(run), encoding="utf-8")
        for gi, (L, members) in enumerate(groups, 1):
            tag = f"[{gi}/{len(groups)}] {L} ({len(members)}개) "
            entry = {"lineage": L, "n": len(members), "designs": []}
            if len(members) < opts["min_n"]:
                entry.update(status="skipped", message=f"게놈이 {len(members)}개라 최소 {opts['min_n']}개 기준에 못 미쳐 건너뜀")
            else:
                ids_by_design = {}
                if engine == "primerhunter":
                    res = _ph_designs(redesign_primerhunter, {k: seqs[k] for k in members if k in seqs}, base / L, opts, seg, ids, scope,
                                      lambda m, tag=tag: progress(tag + m))
                else:
                    res = redesign_varvamp.design_lineage(
                        {k: seqs[k] for k in members if k in seqs}, base / L, opts["threshold"], opts["min_threshold"],
                        ambig=opts["ambig"], probe_ambig=opts["probe_ambig"], top=opts["top"],
                        threads=4, progress=lambda m, tag=tag: progress(tag + m))
                entry.update(status=res["status"], threshold=res.get("threshold"), message=res["message"],
                             sample_n=res.get("sample_n"))
                for di, d in enumerate(res["designs"], 1):
                    out = base / L / "designs" / d["name"]
                    val = validate_oligos(seg, d["F"], d["R"], d["P"], out, f"RD_{seg}_{L}_{d['name']}",
                                          lambda m, tag=tag, di=di, n=len(res["designs"]): progress(f"{tag}설계 {di}/{n} 검증: {m}"),
                                          scope, ids)
                    ids_by_design[d["name"]] = set(val.get("validated_ids") or [])
                    per = {x["lineage"]: x for x in val["lineages"]}
                    own = per.get(L) if L != "ALL" else None
                    sm = val.get("summary") or {}
                    entry["designs"].append({**d, "validated": sm.get("validated"), "validated_pct": sm.get("validated_pct"),
                                             "probe_pct": sm.get("probe_coverage_mm0_pct"), "dropout": sm.get("probe_dropout"),
                                             "std_ok": bool((val.get("standard") or {}).get("ok")),
                                             "own_validated": own["validated"] if own else None, "own_total": own["total"] if own else None,
                                             "lineages": val["lineages"], "dir": f"{L}/designs/{d['name']}"})
                entry["cover"] = _set_cover(ids_by_design, set(members))
                entry["designs"].sort(key=lambda x: (-(x.get("own_validated") if x.get("own_validated") is not None else x.get("validated") or 0),
                                                     -(x.get("validated") or 0)))
            run["lineages"].append(entry)
            (base / "run.json").write_text(json.dumps(run), encoding="utf-8")
        run["done"] = True
        (base / "run.json").write_text(json.dumps(run), encoding="utf-8")
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "done"
            JOBS[job_id]["progress"] = "완료"
    except Exception as exc:  # noqa: BLE001
        traceback.print_exc()
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error"
            JOBS[job_id]["error"] = str(exc)
            JOBS[job_id]["progress"] = f"실행 실패: {exc}"


def tail_progress(log_path: Path) -> str:
    try:
        text = log_path.read_text(encoding="utf-8", errors="replace")
    except FileNotFoundError:
        return ""
    matches = re.findall(r"\d+ of \d+ steps \(\d+%\) done", text)
    if matches:
        return matches[-1]
    for line in reversed(text.splitlines()):
        if line.strip():
            return line.strip()[:200]
    return ""


def tail_log_lines(log_path: Path, n: int = 400) -> str:
    """Raw tail of the snakemake log for the progress modal -- showing the actual log
    (rather than just the parsed 'N of M steps' line) is what the user asked for after
    finding the old compact one-line status too hard to follow."""
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except FileNotFoundError:
        return ""
    return "\n".join(lines[-n:])


def read_tsv(path: Path):
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_manifest():
    """Only proceed=True rows, legacy-shaped (set_id/segment/source/F_id/F_seq/R_id/R_seq/
    P_id/P_seq/fluorophore/add_date) -- the same filtered view the Snakefile itself builds
    its whole rule-wildcard universe from, so a row someone unchecked in the "프라이머
    입력/관리" spreadsheet disappears from the dashboard exactly like it does from the
    pipeline (it has no results/ to show anyway, since nothing computed them)."""
    from mfe_common import load_primers_manifest  # noqa: PLC0415 -- see load_builtin_sets for why
    return {r["set_id"]: r for r in load_primers_manifest(ROOT / "config" / "primers.tsv")}


# Analysis parameters editable from the webapp's settings page. Stored as plain lines in
# config/config.yaml (the same file Snakemake itself reads via `configfile:`), so changing
# a value here and clicking 분석 다시하기 actually re-runs `reconcile.py` with the new
# threshold -- not just a display-only setting. A tiny regex read/write is used instead of
# a YAML library (PyYAML isn't installed in this server's Python env) since config.yaml's
# top-level scalars are simple "key: value" lines; this only ever touches known keys and
# preserves every comment and the nested `segments:` block untouched.
PARAM_KEYS = {
    "primer_mismatch_pass": int,
    "probe_mismatch": int,
    "target_coverage_pct": float,
    "three_prime_consecutive_match": int,
    "mfeprimer_cpu": int,
}
PARAM_DEFAULTS = {
    "primer_mismatch_pass": 1, "probe_mismatch": 0,
    "target_coverage_pct": 95.0, "three_prime_consecutive_match": 2,
    "mfeprimer_cpu": 4,
}


def PROFILE_ID_FOR_FILES():                 # noqa: N802 -- file-name prefix of the downloads (the virus id)
    try:
        return get_profile()["id"]
    except Exception:  # noqa: BLE001
        return "primer_validation"


XLSX_ALIASES = {
    "primers": {
        "primer_id": ("id", "primerid"), "primer_name": ("name", "이름", "프라이머이름", "primername", "세트이름"),
        "segment": ("seg", "분절", "세그먼트"), "forward_primer": ("forward", "fwd", "f", "forwardprimer", "정방향"),
        "reverse_primer": ("reverse", "rev", "r", "reverseprimer", "역방향"), "probe": ("프로브",),
        "fluorophore": ("형광물질", "fluor", "dye", "형광"), "source": ("출처", "vendor", "제조사", "회사"),
        "remark": ("비고", "note", "notes", "메모"), "proceed": ("분석진행", "run", "analyze", "analyse"),
        "add_date": ("추가일", "date", "등록일"),
    },
    "standards": {
        "name": ("id", "이름", "표준물질", "표준물질이름"), "segment": ("seg", "분절", "세그먼트"),
        "sequence": ("seq", "서열", "염기서열"), "source": ("출처", "vendor", "제조사"),
        "remark": ("비고", "note", "notes", "메모"), "proceed": ("분석진행", "run", "analyze", "analyse"),
        "add_date": ("추가일", "date", "등록일"),
    },
}
XLSX_HELP = {
    "primers": [("primer_id", "자동 부여되는 영구 ID(P001…). 새 행은 비워 두세요. 기존 ID를 그대로 두면 그 프라이머를 수정합니다."),
                ("primer_name", "표시 이름(중복 불가)"), ("segment", "분절(단일 게놈 바이러스는 비워도 됨)"),
                ("forward_primer / reverse_primer / probe", "5'→3' 서열. IUPAC 코드(A C G T R Y W S K M B D H V N)만 허용"),
                ("fluorophore / source / remark", "형광물질 / 출처(제조사) / 비고"),
                ("proceed", "True 또는 False. False인 행은 분석에서 제외"), ("add_date", "저장할 때 자동 기록(수정 불필요)")],
    "standards": [("name", "표준물질 이름 = 결과·실험 기록의 ID(영문·숫자·_ . -, 중복 불가)"), ("segment", "분절(단일 게놈 바이러스는 비워도 됨)"),
                  ("sequence", "서열(IUPAC, 공백·줄바꿈은 자동 제거)"), ("source / remark", "출처 / 비고"),
                  ("proceed", "True 또는 False. False인 행은 분석에서 제외"), ("add_date", "저장할 때 자동 기록(수정 불필요)")],
}


def xlsx_export(kind):
    """The current primer / standard-material table as an .xlsx (sheet 1 = the data in the exact column names the
    upload reads back, sheet 2 = what each column means)."""
    from xlsx_lite import write_xlsx  # noqa: PLC0415
    if kind == "primers":
        from mfe_common import PRIMER_MANIFEST_COLUMNS as cols, load_primers_manifest_raw  # noqa: PLC0415
        rows = load_primers_manifest_raw(ROOT / "config" / "primers.tsv")
        title, widths = "Primers", [10, 24, 10, 34, 34, 40, 12, 14, 24, 10, 12]
    else:
        from standards_common import STANDARD_COLUMNS as cols, ensure_standards_manifest  # noqa: PLC0415
        rows = ensure_standards_manifest(ROOT / "config" / "standards.tsv", ROOT / "config", seg_keys())
        title, widths = "Standards", [20, 10, 60, 14, 24, 10, 12]
    data = [list(cols)] + [[(r.get(c) or "") for c in cols] for r in rows]
    help_rows = [["column", "설명"]] + [[a, b] for a, b in XLSX_HELP[kind]]
    return write_xlsx([(title, data, widths[:len(cols)]), ("안내", help_rows, [34, 90])])


def xlsx_validate(kind, rows):
    """Same checks as saving the table (names unique, valid IUPAC sequences, segment, ...), without saving anything:
    -> [{row, field, message}] with row = index into `rows`."""
    if kind == "primers":
        from mfe_common import PRIMER_MANIFEST_COLUMNS as cols, validate_primer_manifest_rows  # noqa: PLC0415
        rows = [{c: str(r.get(c, "")) for c in cols} for r in rows]
        if len(seg_keys()) == 1:
            for r in rows:
                r["segment"] = seg_keys()[0].upper()
        return validate_primer_manifest_rows(rows, segments=tuple(k.upper() for k in seg_keys()))
    from standards_common import STANDARD_COLUMNS as cols, validate_rows  # noqa: PLC0415
    rows = [{c: str(r.get(c, "")) for c in cols} for r in rows]
    for r in rows:
        r["name"] = r["name"].strip()
        r["segment"] = seg_keys()[0].upper() if len(seg_keys()) == 1 else r["segment"].strip().upper()
    return validate_rows(rows, segments=tuple(k.upper() for k in seg_keys()))


def xlsx_parse(kind, data):
    """Uploaded .xlsx -> ({rows: [{column: text}], warnings: [...]}). Header names are matched loosely (case, spaces,
    underscores; Korean and common English aliases). Nothing is saved: the page loads the rows into its table for review."""
    import re as _re  # noqa: PLC0415
    from xlsx_lite import read_xlsx  # noqa: PLC0415
    if kind == "primers":
        from mfe_common import PRIMER_MANIFEST_COLUMNS as cols  # noqa: PLC0415
    else:
        from standards_common import STANDARD_COLUMNS as cols  # noqa: PLC0415
    _name, table = read_xlsx(data)
    norm = lambda s: _re.sub(r"[\s_\-()/]+", "", str(s or "").lower())      # noqa: E731
    lookup = {norm(c): c for c in cols}
    for col, aliases in XLSX_ALIASES[kind].items():
        for a in aliases:
            lookup.setdefault(norm(a), col)
    header_i = next((i for i, r in enumerate(table) if sum(1 for c in r if norm(c) in lookup) >= 2), None)
    if header_i is None:
        raise ValueError("헤더 행을 찾지 못했습니다. 첫 줄에 열 이름(" + ", ".join(cols) + ")을 넣어 주세요 — 먼저 '엑셀 다운로드'로 양식을 받으면 편합니다.")
    colmap = {j: lookup[norm(h)] for j, h in enumerate(table[header_i]) if norm(h) in lookup}
    unknown = [str(h) for j, h in enumerate(table[header_i]) if str(h).strip() and j not in colmap]
    warnings = [f"인식하지 못해 무시한 열: {', '.join(unknown)}"] if unknown else []
    seq_cols = {"forward_primer", "reverse_primer", "probe", "sequence"}
    out, blank = [], 0
    for r in table[header_i + 1:]:
        rec = {c: "" for c in cols}
        for j, col in colmap.items():
            v = (r[j] if j < len(r) else "").strip()
            if col in seq_cols:
                v = _re.sub(r"\s+", "", v).upper()
            elif col == "segment":
                v = v.upper()
            elif col == "proceed":
                v = "False" if v.lower() in ("false", "0", "no", "n", "x", "아니오") else "True"
            rec[col] = v
        if not any(v for k, v in rec.items() if k not in ("proceed", "segment")):
            blank += 1
            continue
        out.append(rec)
    if blank:
        warnings.append(f"빈 행 {blank}개는 건너뛰었습니다.")
    return {"rows": out, "warnings": warnings}


def _config_yaml_path():
    return ROOT / "config" / "config.yaml"


def load_params():
    from mfe_common import check_mfeprimer, find_mfeprimer, read_config_str  # noqa: PLC0415
    text = _config_yaml_path().read_text(encoding="utf-8")
    out = {}
    for key, caster in PARAM_KEYS.items():
        m = re.search(rf"^{re.escape(key)}:\s*([^\s#]+)", text, re.MULTILINE)
        out[key] = caster(m.group(1)) if m else PARAM_DEFAULTS[key]
    # mfeprimer: the configured path (may be empty = auto-discover) and what actually resolves right now
    out["mfeprimer_bin"] = read_config_str(ROOT, "mfeprimer_bin")
    try:
        used = find_mfeprimer("", ROOT)
    except SystemExit:
        used = ""
    out["mfeprimer_status"] = check_mfeprimer(used)
    from mfe_common import find_mafft  # noqa: PLC0415

    def resolved(finder):
        try:
            return finder()
        except SystemExit:                                # the finders exit when nothing is found
            return ""
    out["tools"] = {"env_prefix": read_config_str(ROOT, "env_prefix"),
                    "snakemake": resolved(find_snakemake), "mafft": resolved(find_mafft)}
    return out


def _set_config_str(text, key, value):
    """Set a top-level `key: "value"` line (quoted, so paths with spaces survive), appending it if absent."""
    line = f'{key}: "{value}"'
    pattern = rf"^{re.escape(key)}:.*$"
    if re.search(pattern, text, re.MULTILINE):
        return re.sub(pattern, lambda _m: line, text, count=1, flags=re.MULTILINE)
    return text.rstrip("\n") + f"\n{line}\n"


def save_params(values):
    from mfe_common import check_mfeprimer  # noqa: PLC0415
    path = _config_yaml_path()
    text = path.read_text(encoding="utf-8")
    if "mfeprimer_bin" in values:
        mfe = str(values["mfeprimer_bin"] or "").strip().strip('"')
        if mfe:                                          # empty = auto-discover; anything else must really be mfeprimer
            chk = check_mfeprimer(mfe)
            if not chk["ok"]:
                raise ValueError(f"mfeprimer 경로: {chk['error']}")
        text = _set_config_str(text, "mfeprimer_bin", mfe)
    for key, caster in PARAM_KEYS.items():
        if key not in values:
            continue
        value = caster(values[key])
        line = f"{key}: {value}"
        pattern = rf"^{re.escape(key)}:.*$"
        if re.search(pattern, text, re.MULTILINE):
            text = re.sub(pattern, line, text, count=1, flags=re.MULTILINE)
        else:
            text = text.rstrip("\n") + f"\n{line}\n"
    path.write_text(text, encoding="utf-8")


def _ui_state_path():
    return ROOT / "config" / "ui_state.json"


UI_STATE_DEFAULTS = {"notes": {}, "set_order": [], "group_order": [], "note_images": {}, "verifications": {},
                     "citations": {}}


def load_ui_state():
    path = _ui_state_path()
    try:
        data = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    except (json.JSONDecodeError, OSError):
        data = {}
    for key, default in UI_STATE_DEFAULTS.items():
        data.setdefault(key, dict(default) if isinstance(default, dict) else list(default))
    return data


def save_ui_state(state):
    path = _ui_state_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")


EVALUATION_LABEL = {"ok": "통과", "warn": "보류", "bad": "실패"}


def _set_parents_path():
    return ROOT / "config" / "set_parents.json"


def load_set_parents():
    """Manual parent assignments {child_set_id: parent_set_id or ""} (an empty string = "no parent", overriding the name-based guess)."""
    try:
        data = json.loads(_set_parents_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def save_set_parent(child, parent):
    data = load_set_parents()
    data[child] = parent
    tmp = _set_parents_path().with_suffix(".tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
    tmp.replace(_set_parents_path())


def annotate_set_relations(sets):
    """Adds the hierarchy facts the sidebar/overview use, without touching any pipeline data:
    parent_id (+ parent_origin "name"|"manual"): a variant is a set named "<parent name>-<n>" (or set by hand), one level deep only;
    duplicate_of / duplicate_ids: sets whose F, R and probe are identical to another set's (only the lowest id has duplicate_of = None)."""
    by_id = {x["set_id"]: x for x in sets}
    names = {(x.get("display_name") or x["set_id"]): x["set_id"] for x in sets}
    manual = load_set_parents()
    for x in sets:
        pid, origin = None, ""
        if x["set_id"] in manual:
            pid, origin = (manual[x["set_id"]] or None), ("manual" if manual[x["set_id"]] else "")
        else:
            m = re.match(r"^(.+)-(\d+)$", x.get("display_name") or "")
            if m and m.group(1) in names and names[m.group(1)] != x["set_id"]:
                pid, origin = names[m.group(1)], "name"
        if pid and (pid not in by_id or by_id[pid].get("segment") != x.get("segment") or pid == x["set_id"]):
            pid, origin = None, ""
        x["parent_id"], x["parent_origin"] = pid, origin
    for x in sets:                                              # one level only: a parent that is itself a variant -> its root
        seen = {x["set_id"]}
        while x["parent_id"] and by_id[x["parent_id"]].get("parent_id") and x["parent_id"] not in seen:
            seen.add(x["parent_id"])
            x["parent_id"] = by_id[x["parent_id"]]["parent_id"]
    def num(i):
        m = re.search(r"\d+", i)
        return (int(m.group()) if m else 10**9, i)
    same = {}
    for x in sets:
        key = tuple((x.get(k) or "").strip().upper() for k in ("F_seq", "R_seq", "P_seq"))
        if all(key):
            same.setdefault(key, []).append(x["set_id"])
    for ids in same.values():
        ids.sort(key=num)
        for sid in ids:
            by_id[sid]["duplicate_of"] = ids[0] if sid != ids[0] else None
            by_id[sid]["duplicate_ids"] = [i for i in ids if i != sid]
    return sets


def load_builtin_sets():
    """Merge all the precomputed pipeline outputs into the same shape /api/sets returns
    for a live-analyzed user set, so the frontend can treat both uniformly."""
    from generate_report import build_assessment  # noqa: PLC0415 -- avoid import at module
    # load time: generate_report.py is a "script" first and a module second, so keep this
    # scoped to where it's actually used rather than risk import-order surprises elsewhere.

    manifest = load_manifest()
    summary = {r["set_id"]: r for r in read_tsv(ROOT / "output" / "summary" / "summary.tsv")}
    std_summary = {r["set_id"]: r for r in read_tsv(ROOT / "output" / "summary" / "standard_summary.tsv")}
    std_matches = {}
    for r in read_tsv(ROOT / "output" / "summary" / "standard_matches.tsv"):
        if int(r["validated_count"]) > 0:
            std_matches.setdefault(r["set_id"], []).append(r["lineage"])
    gene = {r["set_id"]: r for r in read_tsv(ROOT / "output" / "summary" / "gene_summary.tsv")}
    ui_state = load_ui_state()
    notes = ui_state["notes"]
    note_images = ui_state.get("note_images", {})
    verifications = ui_state.get("verifications", {})
    citations = ui_state.get("citations", {})

    sets = []
    for set_id, row in manifest.items():
        s = summary.get(set_id, {})
        std = std_summary.get(set_id, {})
        g = gene.get(set_id, {})
        entry = {
            "set_id": set_id, "segment": row["segment"], "source": row["source"],
            "F_id": row["F_id"], "F_seq": row["F_seq"], "R_id": row["R_id"], "R_seq": row["R_seq"],
            "P_id": row["P_id"], "P_seq": row["P_seq"],
            "note": notes.get(set_id, ""),
            "note_images": note_images.get(set_id, []),
            "verifications": verifications.get(set_id, []),
            "display_name": row.get("name") or set_id,
            "fluorophore": row.get("fluorophore", ""),
            "citation": citations.get(set_id),
            "builtin": True,
            "summary": {
                "total_genomes": int(s.get("total_genomes", 0) or 0),
                "primer_coverage_mm0_pct": float(s.get("primer_coverage_mm0_pct", 0) or 0),
                "primer_coverage_mm1_pct": float(s.get("primer_coverage_mm1_pct", 0) or 0),
                "primer_coverage_mm2_pct": float(s.get("primer_coverage_mm2_pct", 0) or 0),
                "probe_coverage_mm0_pct": float(s.get("probe_coverage_mm0_pct", 0) or 0),
                "primer_matched": int(s.get("primer_matched", 0) or 0),
                "probe_dropout": int(s.get("probe_dropout", 0) or 0),
                "validated": int(s.get("validated", 0) or 0),
                "validated_pct": float(s.get("validated_pct", 0) or 0),
            },
            "standard": {
                "available": bool(std),               # False when this set's segment has no standard material (optional)
                "ok": int(std.get("validated", 0) or 0) > 0,
                "matched": sorted(std_matches.get(set_id, [])),
                "perfect_primer": float(std.get("primer_coverage_mm0_pct", 0) or 0) > 0,
            },
            "gene": {
                "target_gene": g.get("target_gene", "NA"),
                "amplicon_start": g.get("amplicon_start", ""),
                "amplicon_end": g.get("amplicon_end", ""),
                "amplicon_size": g.get("amplicon_size", ""),
                "chrom": g.get("chrom", ""),
            },
        }
        # numbers of the active host scope (human by default), counted with the pipeline's own rule
        seg_key = row["segment"].lower()
        entry["summary"].update(scope_ctx().set_summary(set_id, seg_key))
        entry["host_scope"] = scope_ctx().info(seg_key)
        # binding/dimer_hairpin aren't loaded here (kept lean for the list view -- see
        # load_set_detail for the full per-set payload), so this reuses the same
        # build_assessment() heuristic on a partial set_data; it just skips those two
        # checks (both default to "nothing flagged") and is refined with the full picture
        # once a set is opened (load_set_detail recomputes it with binding/dimer_hairpin).
        _, cls, _ = build_assessment(entry)
        entry["evaluation"] = {"cls": cls, "label": EVALUATION_LABEL[cls]}
        sets.append(entry)

    for user_dir in sorted((ROOT / "output" / "user_sets").glob("*")):
        result_path = user_dir / "result.json"
        if result_path.exists():
            try:
                data = json.loads(result_path.read_text(encoding="utf-8"))
                data["builtin"] = False
                _coerce_numeric_strings(data.get("genomes"))
                seg_u = (data.get("segment") or "").lower()
                if data.get("genomes") and seg_u in seg_keys():     # recount for the active host scope
                    data["summary"] = {**(data.get("summary") or {}), **scope_ctx().summary_from_rows(seg_u, {g["genome_id"]: g for g in data["genomes"]})}
                    data["host_scope"] = scope_ctx().info(seg_u)
                _, cls, _ = build_assessment(data)  # full binding/dimer_hairpin already in result.json
                data["evaluation"] = {"cls": cls, "label": EVALUATION_LABEL[cls]}
                sets.append(_slim(data))
            except Exception:
                continue
    return annotate_set_relations(sets)


def _slim(full):
    """Strip the (potentially large) per-genome list for the list view."""
    return {k: v for k, v in full.items() if k != "genomes"}


def parse_filter_report_kept(path: Path):
    if not path.exists():
        return 0
    m = re.search(r"kept=(\d+)", path.read_text(encoding="utf-8"))
    return int(m.group(1)) if m else 0


def load_population_ids(segment: str):
    """Every accession of the curated population for this segment IN THE ACTIVE HOST SCOPE -- used to tell apart 'this set failed to
    detect it' from 'not in our curated DB at all' when overlaying validation status onto the full nextstrain tree (which includes
    genomes we filtered out for length/N content, and now also the genomes of the other host scopes)."""
    v = scope_ctx()
    seg_key = segment.lower()
    keep = v.ids(seg_key)
    return [a for a in v.pop_ids(seg_key) if a in keep]


def load_lineage_totals():
    """Lineage counts of the population per segment in the active host scope (independent of any assay set) -- the denominator for
    "what % of this lineage does this set cover"."""
    v = scope_ctx()
    return {seg.upper(): v.lineage_totals(seg) for seg in seg_keys()}


def _raw_groups():
    from mfe_common import load_groups_raw  # noqa: PLC0415 -- see load_builtin_sets for why
    return load_groups_raw(ROOT / "config" / "set_groups.tsv")


def load_groups_config():
    """Groups as the dashboard/pipeline see them: each group's member list reduced to the
    primers that are currently active (proceed=True), and a group left with none omitted --
    the same rule sync_input_primers() applies to input_primers/groups.tsv. Edits (create /
    rename / delete) work on _raw_groups() instead so an unchecked member isn't lost from the
    file just because it is temporarily inactive."""
    active = load_manifest()
    out = []
    for r in _raw_groups():
        members = [m for m in (r.get("member_ids") or "").split(";") if m in active]
        if members:
            out.append(dict(r, member_ids=";".join(members)))
    return out


def literature_edition():
    """'lasv' / 'nipah' = the virus this web app analyses (decided by its segments), or None."""
    segs = {k.upper() for k in seg_keys()}
    for key, info in literature.VIRUSES.items():
        if segs == set(info["segments"]):
            return key
    return None


def literature_payload():
    """Both viruses' literature lists (read-only data for the menu) plus what can be added to THIS edition's analysis."""
    from mfe_common import load_primers_manifest_raw  # noqa: PLC0415
    manifest = load_primers_manifest_raw(ROOT / "config" / "primers.tsv")
    edition = literature_edition()
    out = {"edition": edition, "columns": literature.COLUMNS, "viruses": {}}
    for key, info in literature.VIRUSES.items():
        if key != edition:                      # each dashboard lists the literature of its own virus only
            continue
        rows = literature.load(ROOT, key)
        for r in rows:
            r["complete"] = literature.complete(r)
            r["in_analysis"] = bool(literature.complete(r) and literature.same_sequences(manifest, r))
            r["can_add"] = bool(key == edition and r["complete"] and not r["in_analysis"])
        out["viruses"][key] = {"label": info["label"], "segments": list(info["segments"]), "rows": rows,
                               "file": str(literature.file_path(ROOT, key).relative_to(ROOT))}
    return out


def sync_pipeline_inputs():
    """Refresh input_primers/ (what the pipeline actually reads) from the spreadsheet and the
    group definitions -- called after every edit so the folder is never stale, even before
    Snakemake's own start-up sync would have fixed it."""
    from mfe_common import sync_input_primers  # noqa: PLC0415
    sync_input_primers(ROOT / "config" / "primers.tsv", ROOT / "config" / "set_groups.tsv", ROOT / "input_primers")


def _entity_exists(entity_id):
    """True for a builtin set_id OR a group_id -- notes/note-images/wet-lab verifications
    are stored in ui_state.json as flat {id: ...} dicts with no concept of "set" vs "group",
    so the same save/upload handlers serve both /api/sets/<id>/... and /api/groups/<id>/...
    once existence is checked against both id spaces instead of load_manifest() alone."""
    if entity_id in load_manifest():
        return True
    return any(r["group_id"] == entity_id for r in load_groups_config())


def load_group_summary():
    """group_summary.tsv recounted for the active host scope: the population is the scope's genomes, validated = the union of the
    members' validated genomes inside it (the cross-binding extra comes on top in _with_cross)."""
    rows = {r["group_id"]: r for r in read_tsv(ROOT / "output" / "summary" / "group_summary.tsv")}
    view = scope_ctx()
    for gid, r in rows.items():
        seg = (r.get("segment") or "").lower()
        if seg not in seg_keys():
            continue
        keep = view.ids(seg)
        n = len(_read_validated_ids(gid) & keep)
        r["total_genomes"], r["validated"] = len(keep), n
        r["validated_pct"] = round(100 * n / len(keep), 2) if keep else 0.0
        r["host_scope"] = view.scope(seg)
    return rows


def _group_segment(group_id):
    for row in load_groups_config():
        if row["group_id"] == group_id:
            return row["segment"].lower()
    return None


def group_cross_ids(group_id):
    """Genomes passed by some cross-binding pairing of the group's oligos (set 1's F + set 2's R + a probe of either set:
    F/R mismatch <= 1 with intact 3' ends, probe perfect) -- invisible to the plain union of each set's own result."""
    path = ROOT / "output" / "groups" / group_id / "combo_mismatch_table.tsv"
    ids = set()
    params = load_params()                                     # the same pass rule as the set results (primer_mismatch_pass, probe_mismatch)
    pass_mm, probe_mm = params["primer_mismatch_pass"], params["probe_mismatch"]
    if path.exists():
        for r in read_tsv(path):
            try:
                f, rr, p = int(r["f_mismatch"]), int(r["r_mismatch"]), int(r["p_mismatch"])
            except (KeyError, TypeError, ValueError):
                continue
            if f <= pass_mm and rr <= pass_mm and p <= probe_mm and r.get("f_3p_ok", "") != "0" and r.get("r_3p_ok", "") != "0":
                ids.add(r["genome_id"])
    return ids


def _read_validated_ids(group_id):
    """A group's union: the genomes any member set validates, recounted from the member sets' own mismatch tables (the same numbers
    as the set cards). The group's validated_ids.txt is only a fallback: it is a snapshot of the last group run and goes stale whenever
    the sets are recomputed (e.g. after the mismatch parameter was changed) without re-running the group."""
    members = next((r["member_ids"].split(";") for r in load_groups_config() if r["group_id"] == group_id), [])
    seg = _group_segment(group_id)
    if members and seg in seg_keys():
        view = scope_ctx()
        try:
            return set().union(*[view.set_validated_ids(m, seg, "all") for m in members if (ROOT / "output" / "primers" / m / "alignment" / "mismatch_table.tsv").exists()])
        except (OSError, KeyError):
            pass
    path = ROOT / "output" / "groups" / group_id / "validated_ids.txt"
    if not path.exists():
        return set()
    return {ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()}


def _with_cross(summary, group_id):
    """summary dict + the cross-inclusive count/percentage next to the union's, counted inside the active host scope."""
    seg = _group_segment(group_id)
    ids = _read_validated_ids(group_id) | group_cross_ids(group_id)
    if seg in seg_keys():
        keep = scope_ctx().ids(seg)
        ids &= keep
        total = len(keep)
    else:
        total = int(summary.get("total_genomes", 0) or 0)
    n = len(ids)
    return {**summary, "validated_with_cross": n, "validated_with_cross_pct": round(100 * n / total, 2) if total else 0.0}


def load_group_list():
    summaries = load_group_summary()
    sets_by_id = {s["set_id"]: s for s in load_builtin_sets()}
    out = []
    for row in load_groups_config():
        gid = row["group_id"]
        s = summaries.get(gid, {})
        out.append({
            "group_id": gid, "label": row["label"], "segment": row["segment"],
            "members": row["member_ids"].split(";"),
            "member_names": {m: (sets_by_id[m]["display_name"] if m in sets_by_id else m) for m in row["member_ids"].split(";")},
            "member_count": len(row["member_ids"].split(";")),
            "summary": _with_cross({
                "total_genomes": int(s.get("total_genomes", 0) or 0),
                "validated": int(s.get("validated", 0) or 0),
                "validated_pct": float(s.get("validated_pct", 0) or 0),
            }, gid),
        })
    return out


def _load_group_combo_table(path):
    """Shared row-shape for combo_mismatch_table.tsv / standard_combo_table.tsv /
    refgenome_combo_table.tsv -- all written by compute_group_combo_table() with the same
    columns, just against different databases (population / KDCA standard / reference
    genome)."""
    if not path.exists():
        return []
    rows = read_tsv(path)
    # Rows carry permanent ids (P007-F x P008-R); add human-readable counterparts so the UI can
    # show primer names without ever having to know the id scheme.
    names = {sid: (r.get("name") or sid) for sid, r in load_manifest().items()}
    pretty = lambda text: re.sub(r"\bP\d{3,}\b", lambda m: names.get(m.group(0), m.group(0)), text or "")  # noqa: E731
    for g in rows:
        g["combo_label"] = pretty(g.get("combo_id", ""))
        for k in ("f_set_id", "r_set_id", "p_set_id"):
            g[k.replace("_id", "_name")] = pretty(g.get(k, ""))
    int_keys = ("f_mismatch", "r_mismatch", "p_mismatch", "amplicon_start", "amplicon_end",
                "probe_start", "probe_end")
    float_keys = ("f_tm", "f_dg", "r_tm", "r_dg", "p_tm", "p_dg")
    for g in rows:
        g["is_same_set"] = g.get("is_same_set") == "1"
        for k in ("f_3p_ok", "r_3p_ok"):                 # "1" / "0" / "" -> True / False / None (verdictFor reads === false)
            g[k] = {"1": True, "0": False}.get(g.get(k, ""), None)
        for k in int_keys:
            if g.get(k) == "":
                g[k] = None
            elif g.get(k) is not None:
                try:
                    g[k] = int(g[k])
                except ValueError:
                    pass
        for k in float_keys:
            if g.get(k) == "":
                g[k] = None
            elif g.get(k) is not None:
                try:
                    g[k] = float(g[k])
                except ValueError:
                    pass
    return rows


def load_standard_material_ids(segment):
    """Every standard-material (KDCA positive-control) id for a segment, straight from the
    fasta index -- the same universe of materials every set's own standard_genomes table is
    built against, so it's a valid "표준물질 선택" list for a group even though a group has
    no single combined standard_genomes table of its own."""
    fai_path = ROOT / "config" / f"kdca_standard_{segment.lower()}.fasta.fai"
    if not fai_path.exists():
        return []
    return [line.split("\t")[0] for line in fai_path.read_text(encoding="utf-8").splitlines() if line.strip()]


def load_group_detail(group_id):
    cfg = {r["group_id"]: r for r in load_groups_config()}
    row = cfg.get(group_id)
    if row is None:
        return None
    seg = row["segment"]
    members = row["member_ids"].split(";")
    summary = load_group_summary().get(group_id, {})

    validated_ids = _read_validated_ids(group_id)           # the members' union, recounted (see there)
    view = scope_ctx()
    keep_ids = view.ids(seg.lower())                      # the genomes of the active host scope
    validated_ids &= keep_ids

    meta_path = ROOT / "data" / "filtered" / f"{seg.lower()}_filtered_metadata.tsv"
    meta = {r["accession"]: r for r in read_tsv(meta_path)}
    lineage_counts = {}
    for gid in validated_ids:
        lin = meta.get(gid, {}).get("lineage", "NA")
        lineage_counts[lin] = lineage_counts.get(lin, 0) + 1

    sets_by_id = {s["set_id"]: s for s in load_builtin_sets()}
    member_details = [dict(sets_by_id[m]) for m in members if m in sets_by_id]
    for m in member_details:
        binding = read_tsv(ROOT / "output" / "primers" / m["set_id"] / "refgenome_binding" / "binding.tsv")
        for b in binding:
            b["found"] = b.get("found") == "1"
        m["binding"] = binding
    # One flat binding list across every member set (each row still tagged with its own
    # set_id) -- lets the group page reuse the exact same pi-overlay/fine-tuning-suggestion
    # code the individual set page uses (bindingOverlaysFromSet(), compute_shift_suggestions())
    # unmodified, just fed every member's F/P/R sites at once instead of one set's.
    group_binding = [dict(b, set_id=m["set_id"], set_name=m.get("display_name") or m["set_id"]) for m in member_details for b in m["binding"]]

    group_dir = ROOT / "output" / "groups" / group_id
    combo_genomes = _load_group_combo_table(group_dir / "combo_mismatch_table.tsv")
    standard_combo = _load_group_combo_table(group_dir / "standard_combo_table.tsv")
    refgenome_combo = _load_group_combo_table(group_dir / "refgenome_combo_table.tsv")

    # A genome also counts for the group when ANY pairing of its oligos passes (set 1's F + set 2's R + a probe of either
    # set: F/R mismatch <= 1 with intact 3' ends, probe perfect) -- the plain union of each set's own result cannot see that.
    cross_ids = {r["genome_id"] for r in combo_genomes
                 if r.get("f_mismatch") not in (None, "") and r.get("r_mismatch") not in (None, "") and r.get("p_mismatch") not in (None, "")
                 and r["f_mismatch"] <= 1 and r["r_mismatch"] <= 1 and r["p_mismatch"] == 0
                 and r.get("f_3p_ok") is not False and r.get("r_3p_ok") is not False}
    cross_ids &= keep_ids
    combo_genomes = view.tag_and_filter(seg.lower(), combo_genomes)      # only the genomes of the active host scope
    validated_with_cross = validated_ids | cross_ids
    summary_with_cross = _with_cross({"total_genomes": int(summary.get("total_genomes", 0) or 0)}, group_id)
    lineage_counts_with_cross = {}
    for gid in validated_with_cross:
        lin = meta.get(gid, {}).get("lineage", "NA")
        lineage_counts_with_cross[lin] = lineage_counts_with_cross.get(lin, 0) + 1

    size_filter_path = group_dir / "size_filter.json"
    size_filter = json.loads(size_filter_path.read_text(encoding="utf-8")) if size_filter_path.exists() else None

    ui_state = load_ui_state()
    group_summary = {
        "total_genomes": int(summary.get("total_genomes", 0) or 0),
        "validated": int(summary.get("validated", 0) or 0),
        "validated_pct": float(summary.get("validated_pct", 0) or 0),
        "validated_with_cross": summary_with_cross["validated_with_cross"],
        "validated_with_cross_pct": summary_with_cross["validated_with_cross_pct"],
    }
    dimer_hairpin = compute_group_dimer_hairpin(group_id, member_details)

    return {
        "group_id": group_id, "label": row["label"], "segment": seg, "members": members,
        "member_names": {m["set_id"]: (m.get("display_name") or m["set_id"]) for m in member_details},
        "summary": group_summary,
        "lineage_counts": lineage_counts,
        "lineage_totals": load_lineage_totals().get(seg, {}),
        "member_details": member_details,
        "validated_ids": sorted(validated_ids),
        "validated_ids_with_cross": sorted(validated_with_cross),
        "lineage_counts_with_cross": lineage_counts_with_cross,
        "population_ids": load_population_ids(seg),
        "host_scope": view.info(seg.lower()),
        "combo_genomes": combo_genomes,
        "size_filter": size_filter,
        "standard_combo": standard_combo,
        "refgenome_combo": refgenome_combo,
        "dimer_hairpin": dimer_hairpin,
        "evaluation": build_group_evaluation(group_summary, member_details, dimer_hairpin, combo_genomes),
        "binding": group_binding,
        "shift_suggestions": compute_shift_suggestions(seg, group_binding),
        "standard_materials": load_standard_material_ids(seg),
        "note": ui_state["notes"].get(group_id, ""),
        "note_images": ui_state.get("note_images", {}).get(group_id, []),
        "verifications": ui_state.get("verifications", {}).get(group_id, []),
    }


def compute_group_dimer_hairpin(group_id, member_details):
    """Self/cross-dimer & hairpin risk across every F/P/R oligo from every member set in
    this group -- the same mfeprimer dimer/hairpin check run_dimer_hairpin.py runs for one
    set, just on the combined oligo pool so set A's Probe accidentally binding set B's
    Primer (not just its own) is caught too. Computed on demand (two small subprocess
    calls, no population DB involved) rather than as a cached Snakemake rule -- same
    on-demand pattern live_analyze.py already uses for an ad-hoc single set."""
    from mfe_common import find_mfeprimer  # noqa: PLC0415 -- see load_builtin_sets for why
    from run_dimer_hairpin import (  # noqa: PLC0415
        run_mfeprimer, parse_dimer_report, parse_hairpin_report,
        summarize_dimers, summarize_hairpins, DG_FLAG_THRESHOLD,
    )
    records = []
    for s in member_details:
        for role, seq in (("F", s.get("F_seq")), ("P", s.get("P_seq")), ("R", s.get("R_seq"))):
            if seq:
                records.append((f"{s['set_id']}:{role}", seq))
    if len(records) < 2:
        return {"dimer": {}, "hairpin": {}, "dg_flag_threshold": DG_FLAG_THRESHOLD}

    mfeprimer_bin = find_mfeprimer("", ROOT)
    cpu = load_params()["mfeprimer_cpu"]
    work_dir = ROOT / "output" / "groups" / group_id / "dimer_hairpin"
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = work_dir / "oligos.fasta"
    fasta_path.write_text("".join(f">{oid}\n{seq}\n" for oid, seq in records), encoding="utf-8")
    dimer_txt = run_mfeprimer(mfeprimer_bin, "dimer", fasta_path, work_dir / "dimer.txt",
                               extra_args=["--cpu", str(cpu)] if cpu else None)
    hairpin_txt = run_mfeprimer(mfeprimer_bin, "hairpin", fasta_path, work_dir / "hairpin.txt")
    return {
        "dimer": {
            pair: {"worst_dg": v["worst"]["dg"], "worst_score": v["worst"]["score"],
                   "align": v["worst"]["align"], "flagged_count": v["flagged"]}
            for pair, v in summarize_dimers(parse_dimer_report(dimer_txt)).items()
        },
        "hairpin": {
            role: {"worst_dg": v["worst"]["dg"], "worst_tm": v["worst"]["tm"],
                   "worst_score": v["worst"]["score"], "align": v["worst"]["align"],
                   "flagged_count": v["flagged"]}
            for role, v in summarize_hairpins(parse_hairpin_report(hairpin_txt)).items()
        },
        "dg_flag_threshold": DG_FLAG_THRESHOLD,
    }


def build_group_evaluation(group_summary, member_details, dimer_hairpin, combo_genomes):
    """Reuses generate_report.build_assessment()'s exact scoring/wording for a group,
    feeding it a group-shaped set_data: population coverage from the union's own
    validated_pct, "standard material match" as true if ANY member set already matches
    (consistent with how a validated genome only needs one member set to pass to count
    toward the union), the group-wide dimer_hairpin computed above, and combo_genomes (the
    population-wide cross-binding search from group_cross_binding.py) so unintended
    cross-set combos actually count against the verdict instead of only being shown as a
    raw table. Per-set-only checks that don't generalize to a union (probe dropout, single
    amplicon size, single-genome binding-not-found) are left out by passing empty/absent
    inputs for them -- build_assessment() already skips a check gracefully when its inputs
    are missing, so this doesn't need its own parallel heuristic."""
    from generate_report import build_assessment  # noqa: PLC0415 -- see load_builtin_sets
    std_ok = any((m.get("standard") or {}).get("ok") for m in member_details)
    matched = sorted({lin for m in member_details if (m.get("standard") or {}).get("ok")
                       for lin in ((m.get("standard") or {}).get("matched") or [])})
    group_set_data = {
        "summary": {"validated_pct": group_summary["validated_pct"]},
        "standard": {"ok": std_ok, "matched": matched},
        "gene": {},
        "binding": [],
        "dimer_hairpin": dimer_hairpin,
        "cross_binding": combo_genomes,
    }
    _, cls, _ = build_assessment(group_set_data)
    return {"cls": cls, "label": EVALUATION_LABEL[cls]}


def append_group(label, segment, member_ids):
    """New group with the next free permanent id (G001, ...). Returns that id."""
    from mfe_common import write_groups  # noqa: PLC0415
    path = ROOT / "config" / "set_groups.tsv"
    rows = _raw_groups()
    rows.append({"group_id": "", "label": label, "segment": segment, "member_ids": ";".join(member_ids)})
    saved = write_groups(path, rows)
    sync_pipeline_inputs()
    return saved[-1]["group_id"]


def remove_group(group_id):
    """Delete a group: drop its row from set_groups.tsv and its precomputed results.
    Returns False if the group didn't exist."""
    from mfe_common import write_groups  # noqa: PLC0415
    rows = _raw_groups()
    remaining = [r for r in rows if r["group_id"] != group_id]
    if len(remaining) == len(rows):
        return False
    write_groups(ROOT / "config" / "set_groups.tsv", remaining)
    group_dir = ROOT / "output" / "groups" / group_id
    if group_dir.exists():
        shutil.rmtree(group_dir)
    sync_pipeline_inputs()
    return True


def rename_group(group_id, new_label):
    """Change a group's display label in place -- group_id (and every file path derived
    from it under output/groups/<group_id>/) is untouched, so this never needs to
    move or recompute anything. Returns False if the group didn't exist."""
    from mfe_common import write_groups  # noqa: PLC0415
    rows = _raw_groups()
    row = next((r for r in rows if r["group_id"] == group_id), None)
    if row is None:
        return False
    row["label"] = new_label
    write_groups(ROOT / "config" / "set_groups.tsv", rows)
    sync_pipeline_inputs()
    return True


def load_group_lineage_rows():
    """Same {set_id, lineage, validated_count} row shape as summary_by_lineage.tsv (which
    only covers individual sets) but for each group's own union of validated genomes --
    computed on the fly from each group's validated_ids.txt + segment metadata (mirrors
    load_group_detail()'s per-group lineage_counts, just for every group at once so the
    overview's merged coverage+lineage table can include groups too, not just sets)."""
    meta_by_seg = {}
    rows = []
    for row in load_groups_config():
        gid, seg = row["group_id"], row["segment"]
        validated_ids = _read_validated_ids(gid)
        if not validated_ids and not (ROOT / "output" / "groups" / gid / "validated_ids.txt").exists():
            continue
        validated_ids |= group_cross_ids(gid)             # the overview heatmap shows the cross-inclusive coverage
        if seg.lower() in seg_keys():
            validated_ids &= scope_ctx().ids(seg.lower())         # ... inside the active host scope
        if seg not in meta_by_seg:
            meta_path = ROOT / "data" / "filtered" / f"{seg.lower()}_filtered_metadata.tsv"
            meta_by_seg[seg] = {r["accession"]: r for r in read_tsv(meta_path)}
        meta = meta_by_seg[seg]
        counts = {}
        for accid in validated_ids:
            lin = meta.get(accid, {}).get("lineage", "NA") or "NA"
            counts[lin] = counts.get(lin, 0) + 1
        for lin, cnt in sorted(counts.items()):
            rows.append({"set_id": gid, "lineage": lin, "validated_count": str(cnt)})
    return rows


def load_set_lineage_rows():
    """{set_id, segment, lineage, validated_count} for every set: its validated genomes inside the active host scope per lineage
    (summary_by_lineage.tsv holds the same for the primary scope only)."""
    view = scope_ctx()
    rows = []
    for sid, r in load_manifest().items():
        seg = r["segment"].lower()
        p = ROOT / "output" / "primers" / sid / "final" / "validated_ids.txt"
        if seg not in seg_keys() or not p.exists():
            continue
        ids = {ln.strip() for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()} & view.ids(seg)
        for lin, cnt in sorted(view.lineage_counts(seg, ids).items(), key=lambda kv: -kv[1]):
            rows.append({"set_id": sid, "segment": r["segment"].upper(), "lineage": lin, "validated_count": str(cnt)})
    return rows


FINAL_VERDICTS = ("adopt", "hold", "reject")


def _final_results_path():
    return ROOT / "config" / "final_results.json"


_PURGE_SUMMARIES = ("summary.tsv", "summary_by_lineage.tsv", "gene_summary.tsv", "engine_concordance.tsv", "standard_summary.tsv",
                    "standard_matches.tsv")


def purge_set_results(set_ids):
    """Deletes what the pipeline and the dashboard computed for sets whose "분석진행" is off: output/primers/<id>/, the probe-design folder and
    the set's rows in the summary tables (the Parameters-independent results can be recomputed by switching the set back on and re-running).
    The primer sequences stay in config/primers.tsv. -> [{"id", "dirs", "rows"}] for the sets that actually had something to delete."""
    report = []
    ids = {i for i in set_ids if re.fullmatch(r"[A-Za-z0-9_\-]+", i or "")}
    if not ids:
        return report
    removed = {i: {"dirs": 0, "rows": 0} for i in ids}
    for i in ids:
        for d in (ROOT / "output" / "primers" / i, ROOT / "output" / "probe_design" / i):
            if d.is_dir():
                shutil.rmtree(d)
                removed[i]["dirs"] += 1
    for name in _PURGE_SUMMARIES:
        f = ROOT / "output" / "summary" / name
        if not f.exists():
            continue
        lines = f.read_text(encoding="utf-8").split("\n")
        head = lines[0].split("\t") if lines and lines[0] else []
        if "set_id" not in head:
            continue
        col = head.index("set_id")
        keep, changed = [], False
        for ln in lines:
            parts = ln.split("\t")
            if ln and ln is not lines[0] and len(parts) > col and parts[col] in ids:
                removed[parts[col]]["rows"] += 1
                changed = True
                continue
            keep.append(ln)
        if changed:
            tmp = f.with_suffix(".tmp")
            tmp.write_text("\n".join(keep), encoding="utf-8")
            tmp.replace(f)
    fr = load_final_results()                                   # the user's own marks and parent links of a removed set go too
    for i in ids:
        if i in fr:
            save_final_result(i, "", "")
    for i, v in removed.items():
        if v["dirs"] or v["rows"]:
            report.append({"id": i, **v})
    return report


def load_final_results():
    """The user's own final call per primer set / group: {id: {"verdict": adopt|hold|reject, "note": str, "updated": "YYYY-MM-DD HH:MM"}}.
    Purely a label the people using the dashboard set; nothing in the pipeline reads it."""
    try:
        data = json.loads(_final_results_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def save_final_result(item_id, verdict, note):
    data = load_final_results()
    if verdict in FINAL_VERDICTS or (note or "").strip():
        data[item_id] = {"verdict": verdict if verdict in FINAL_VERDICTS else "", "note": (note or "").strip()[:500],
                         "updated": time.strftime("%Y-%m-%d %H:%M")}
    else:
        data.pop(item_id, None)                                   # nothing set any more: the entry disappears
    tmp = _final_results_path().with_suffix(".tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
    tmp.replace(_final_results_path())
    return data.get(item_id)


def load_pipelines():
    """The dashboards listed on the front-door page (entry/pipelines.json, two levels above a project folder), for the LASV / NiV switch in
    the sidebar: [{id, name, port, current}]. Empty when this project is used on its own."""
    reg = ROOT.parent.parent / "entry" / "pipelines.json"
    try:
        items = json.loads(reg.read_text(encoding="utf-8")).get("pipelines", [])
    except (OSError, ValueError):
        return []
    here = f"{ROOT.parent.name}/{ROOT.name}"
    return [{"id": p["id"], "name": p.get("short") or p["name"], "port": p["port"], "current": p.get("path") == here} for p in items]


def _venn_saved_path():
    return ROOT / "config" / "venn_saved.json"


def load_venn_saved():
    """Named selections of the Venn page ([{id, name, seg, keys, lineage, created}]), shared by everyone who uses this dashboard."""
    try:
        data = json.loads(_venn_saved_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    return data if isinstance(data, list) else []


def _write_venn_saved(items):
    tmp = _venn_saved_path().with_suffix(".tmp")
    tmp.write_text(json.dumps(items, ensure_ascii=False, indent=1), encoding="utf-8")
    tmp.replace(_venn_saved_path())


def load_venn(seg):
    """Everything the coverage-Venn page needs for one segment, in the active host scope: the population (accession, lineage, country)
    and, per primer set and per group, the genomes it validates (a group counts its cross-binding genomes too, like the overview)."""
    v = scope_ctx()
    keep = v.ids(seg)
    lin = v.lineage_of(seg)
    pop = {r["accession"]: r for r in v.population(seg)}
    genomes = [[i, lin.get(i, "NA"), (pop.get(i) or {}).get("country") or ""] for i in sorted(keep)]
    items = []
    for st in load_builtin_sets():
        if str(st.get("segment", "")).lower() != seg:
            continue
        try:                                  # recounted like the set card (same mismatch table, same parameters, same host scope)
            ids = v.set_validated_ids(st["set_id"], seg)
        except (OSError, KeyError):
            ids = set()
        items.append({"key": st["set_id"], "kind": "set", "name": st.get("display_name") or st["set_id"], "ids": sorted(ids & keep)})
    for g in load_groups_config():
        if str(g.get("segment", "")).lower() != seg:
            continue
        ids = (_read_validated_ids(g["group_id"]) | group_cross_ids(g["group_id"])) & keep
        items.append({"key": g["group_id"], "kind": "group", "name": g.get("label") or g["group_id"], "ids": sorted(ids)})
    return {"segment": seg.upper(), "host_scope": v.info(seg), "genomes": genomes, "items": items}


def load_overview():
    sets = load_builtin_sets()
    groups = list(load_group_summary().values())
    names = {s["set_id"]: s["display_name"] for s in sets}
    for g in groups:
        g["member_names"] = [names.get(m, m) for m in (g.get("members") or "").split(";") if m]
        g.update(_with_cross({"total_genomes": g.get("total_genomes", 0)}, g["group_id"]))
        g["validated_with_cross"], g["validated_with_cross_pct"] = g.pop("validated_with_cross"), g.pop("validated_with_cross_pct")
    lineage_rows = load_set_lineage_rows() + load_group_lineage_rows()
    view = scope_ctx()
    dataset_sizes = {k.upper(): len(view.ids(k)) for k in seg_keys()}      # genomes in the active host scope
    return {"sets": sets, "groups": groups, "lineage": lineage_rows, "dataset_sizes": dataset_sizes,
            "lineage_totals": load_lineage_totals(), "host_scope": host_scope_info()}


_DIVERSITY_FLIP_CACHE = {}


def _detect_diversity_orientation_flip(seg, genome_len):
    """Whether run_pixy's population VCF ended up with its reference row reverse-complemented
    relative to the raw reference fasta (the coordinate system refgenome_genes.tsv and every
    set's refgenome_binding/*/binding.tsv use). build_diversity_vcf.py's normalize_orientation()
    is specifically meant to prevent this (mafft --adjustdirection can flip ANY row, including
    the reference, to whatever shared orientation it lands on), but empirically the S segment's
    most recent run has it flipped anyway: the VCF's own REF column, reassembled into a
    sequence, is the exact reverse complement of data/refgenome/s_ref.fasta (100% identity
    either way -- confirmed by direct comparison, not a fuzzy call). Detected here rather than
    assumed correct, by sampling REF alleles straight out of the (already bgzipped, so a plain
    gzip reader works fine sequentially) VCF and comparing them against the raw reference in
    both orientations -- so this keeps working (as a no-op) whenever a future pipeline re-run
    fixes the root cause, and would also catch the L segment if it ever happened there too."""
    vcf_path = ROOT / "output" / "diversity" / seg / "allsites.vcf.gz"
    ref_fasta_path = ROOT / "data" / "refgenome" / f"{seg}_ref.fasta"
    if not vcf_path.exists() or not ref_fasta_path.exists():
        return False
    ref_seq = "".join(
        line.strip() for line in ref_fasta_path.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    ).upper()
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    same, flipped = 0, 0
    with gzip.open(vcf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 5)
            pos, ref_allele = int(fields[1]), fields[3]
            if len(ref_allele) != 1 or ref_allele not in complement:
                continue
            if 1 <= pos <= len(ref_seq) and ref_seq[pos - 1] == ref_allele:
                same += 1
            mirrored = genome_len - pos + 1
            if 1 <= mirrored <= len(ref_seq) and ref_seq[mirrored - 1] == complement[ref_allele]:
                flipped += 1
            if same + flipped >= 300:
                break
    return flipped > same


def _diversity_is_flipped(seg, genome_len):
    if seg not in _DIVERSITY_FLIP_CACHE:
        _DIVERSITY_FLIP_CACHE[seg] = _detect_diversity_orientation_flip(seg, genome_len)
    return _DIVERSITY_FLIP_CACHE[seg]


def _mirror_diversity_windows(windows_by_lineage, genome_len):
    """Flips every lineage's window list from the VCF's (wrong-way) coordinate system back
    into the reference's own -- start/end swap-and-mirror (genome_len - end + 1, genome_len -
    start + 1) since the interval itself reverses direction, re-sorted so consumers can still
    assume ascending position order."""
    out = {}
    for lin, windows in windows_by_lineage.items():
        mirrored = [dict(w, start=genome_len - w["end"] + 1, end=genome_len - w["start"] + 1) for w in windows]
        mirrored.sort(key=lambda w: w["start"])
        out[lin] = mirrored
    return out


def _load_standard_rep_raw(seg):
    path = ROOT / "output" / "standards" / seg / "representativeness.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def load_standard_representativeness():
    """Both segments' representativeness.json (written by the standard_representativeness rule),
    trimmed for the page: the long per-genome / per-set unmatched lists stay server-side and are
    served as downloads instead."""
    names = {sid: (r.get("name") or sid) for sid, r in load_manifest().items()}
    out = {"segments": {}, "set_names": names}
    for seg in seg_keys():
        d = _load_standard_rep_raw(seg)
        if d is None:
            continue
        d = dict(d)
        d.pop("genome_lineage", None)
        d["union"] = dict(d["union"])
        d["union"]["unmatched_detail"] = d["union"].get("unmatched_detail", [])[:50]
        sets = {}
        for sid, e in d["sets"].items():
            e = dict(e)
            if "union" in e:
                e["union"] = {k: v for k, v in e["union"].items() if k != "unmatched"}
            sets[sid] = e
        d["sets"] = sets
        out["segments"][seg] = d
    return out


def standard_rep_for_set(set_id, segment):
    """Compact per-set verdict for the set detail page (None until the rule has run)."""
    d = _load_standard_rep_raw(segment.lower())
    if d is None:
        return None
    e = dict(d["sets"].get(set_id) or {"applicable": False, "reason": "unknown", "grade": "na"})
    if "union" in e:
        e["union"] = {k: v for k, v in e["union"].items() if k != "unmatched"}
    return {"segment_union": {k: d["union"][k] for k in ("assessable", "matched", "pct", "grade", "sensitivity", "by_lineage")},
            "params": d["params"], "standards": d["standards"], "set": e}


def standard_rep_unmatched_tsv(seg, set_id=None):
    d = _load_standard_rep_raw(seg)
    if d is None:
        return None
    lin = d.get("genome_lineage", {})
    if set_id:
        e = d["sets"].get(set_id) or {}
        ids = (e.get("union") or {}).get("unmatched", [])
        return "genome_id\tlineage\n" + "".join(f"{g}\t{lin.get(g, 'NA')}\n" for g in ids)
    rows = d["union"].get("unmatched_detail", [])
    return "genome_id\tlineage\tbest_standard\tbest_identity\n" + "".join(
        f"{r['genome_id']}\t{r['lineage']}\t{r['best_standard']}\t{r['best_identity']}\n" for r in rows)


def load_diversity():
    """Per-lineage nucleotide diversity (pixy pi), windowed across the S/L segment --
    written by the run_pixy Snakemake rule. Segments not yet computed (pipeline hasn't been
    rerun since this feature was added) are simply omitted rather than erroring."""
    out = {}
    for seg in seg_keys():
        path = ROOT / "output" / "diversity" / seg / "pixy_pi.json"
        if path.exists():
            out[seg] = json.loads(path.read_text(encoding="utf-8"))
    genes = [r for r in read_tsv(ROOT / "config" / "refgenome_genes.tsv")]
    genes_by_seg = {}
    genome_len_by_seg = {}
    for r in genes:
        seg = r["segment"].lower()
        genes_by_seg.setdefault(seg, []).append(
            {"gene": r["gene"], "start": int(r["start"]), "end": int(r["end"]), "strand": r["strand"]})
        genome_len_by_seg[seg] = int(r["genome_len"])
    for seg, seg_data in out.items():
        genome_len = genome_len_by_seg.get(seg)
        if genome_len and "windows" in seg_data and _diversity_is_flipped(seg, genome_len):
            seg_data["windows"] = _mirror_diversity_windows(seg_data["windows"], genome_len)
    lineage_totals = load_lineage_totals()
    lineage_totals_by_seg = {"s": lineage_totals.get("S", {}), "l": lineage_totals.get("L", {})}
    return {"segments": out, "genes": genes_by_seg, "genome_len": genome_len_by_seg,
            "lineage_totals": lineage_totals_by_seg}


# ---------------------------------------------------------------------------------------
# Fine-tuning shift suggestion: for each F/P/R binding site, search a small +-40bp window
# along the population-pooled ALL_LASV pi track (now correctly oriented -- see
# _detect_diversity_orientation_flip above) for a same-length spot that's meaningfully less
# diverse, and reconstruct what oligo sequence would actually sit there (IUPAC-coded at
# positions where >=3% of the population disagrees). Lives in server.py (not the frontend)
# because building the candidate sequence needs to read the population VCF directly --
# something only this process has cheap filesystem access to do.
SHIFT_SEARCH_RANGE = 40
SHIFT_MIN_ABS_PI = 0.02
SHIFT_MIN_IMPROVEMENT_RATIO = 0.5
SHIFT_IUPAC_MIN_FRAC = 0.03  # only call an alt allele if >=3% of the population carries it

_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
_IUPAC_FROM_BASES = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H", frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}
_IUPAC_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


def _revcomp_iupac(seq):
    return "".join(_IUPAC_COMPLEMENT[c] for c in reversed(seq))


def _true_oligo_range(binding_row):
    """binding.tsv's start is 0-based (mfeprimer's own convention, distinct from the 1-based
    GenBank numbering refgenome_genes.tsv/pixy use) -- true 1-based inclusive range is
    [start+1, start+len(seq)]. Confirmed empirically: scanning the whole (orientation-
    corrected) reference for a known "+"-strand primer's exact sequence finds a perfect
    match only at this offset, one off from treating `start` as already 1-based."""
    lo = int(binding_row["start"]) + 1
    hi = lo + len(binding_row["seq"]) - 1
    return lo, hi


def _avg_pi_for_range(windows, lo, hi):
    vals = [w["pi"] for w in windows if w["pi"] is not None and w["end"] >= lo and w["start"] <= hi]
    return sum(vals) / len(vals) if vals else None


def _query_vcf_allele_counts(seg, vlo, vhi):
    """Per-VCF-position {base: count} across every sample's called genotype, for POS in
    [vlo, vhi] -- reads the bgzipped VCF sequentially with plain gzip (no bcftools/tabix
    dependency; bgzip is a valid concatenated-gzip stream) since this is only ever called for
    a small (<=~120bp) window."""
    vcf_path = ROOT / "output" / "diversity" / seg / "allsites.vcf.gz"
    out = {}
    with gzip.open(vcf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            pos = int(fields[1])
            if pos < vlo:
                continue
            if pos > vhi:
                break
            ref_allele = fields[3]
            alt_alleles = fields[4].split(",") if fields[4] != "." else []
            alleles = [ref_allele] + alt_alleles
            counts = {}
            for gt in fields[9:]:
                gt = gt.strip()
                if gt in (".", "./."):
                    continue
                try:
                    idx = int(gt.split("/")[0].split("|")[0])
                except ValueError:
                    continue
                base = alleles[idx] if idx < len(alleles) else None
                if base in "ACGT":
                    counts[base] = counts.get(base, 0) + 1
            out[pos] = counts
    return out


def _true_orientation_allele_counts(seg, genome_len, lo, hi):
    """Per-true-reference-position (ascending lo..hi) allele counts, translated out of the
    VCF's own coordinate/orientation system when that segment's pixy run flipped it."""
    if not _diversity_is_flipped(seg, genome_len):
        return _query_vcf_allele_counts(seg, lo, hi)
    vlo, vhi = genome_len - hi + 1, genome_len - lo + 1
    vcf_counts = _query_vcf_allele_counts(seg, vlo, vhi)
    out = {}
    for t in range(lo, hi + 1):
        v = genome_len - t + 1
        raw = vcf_counts.get(v, {})
        out[t] = {_COMPLEMENT[b]: c for b, c in raw.items()}
    return out


def _query_vcf_sample_bases(seg, vlo, vhi):
    """Per-sample base at each VCF POS in [vlo, vhi] -- {sample_id: {pos: base}}, "N" for a
    missing/no-call genotype. Same sequential-gzip-scan approach as _query_vcf_allele_counts,
    just keeping each sample's own call instead of pooling them into counts."""
    vcf_path = ROOT / "output" / "diversity" / seg / "allsites.vcf.gz"
    samples = []
    out = {}
    with gzip.open(vcf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                out = {s: {} for s in samples}
                continue
            fields = line.rstrip("\n").split("\t")
            pos = int(fields[1])
            if pos < vlo:
                continue
            if pos > vhi:
                break
            ref_allele = fields[3]
            alt_alleles = fields[4].split(",") if fields[4] != "." else []
            alleles = [ref_allele] + alt_alleles
            for sample, gt in zip(samples, fields[9:]):
                gt = gt.strip()
                base = "N"
                if gt not in (".", "./."):
                    try:
                        idx = int(gt.split("/")[0].split("|")[0])
                        cand = alleles[idx] if idx < len(alleles) else None
                        if cand in "ACGT":
                            base = cand
                    except ValueError:
                        pass
                out[sample][pos] = base
    return out


def _true_orientation_sample_bases(seg, genome_len, lo, hi):
    """Per-sample base at each true-reference position lo..hi (ascending), translated out of
    the VCF's own coordinate/orientation system when that segment's pixy run flipped it --
    same logic as _true_orientation_allele_counts but keeping sample identity."""
    if not _diversity_is_flipped(seg, genome_len):
        return _query_vcf_sample_bases(seg, lo, hi)
    vlo, vhi = genome_len - hi + 1, genome_len - lo + 1
    vcf_bases = _query_vcf_sample_bases(seg, vlo, vhi)
    out = {}
    for sample, by_pos in vcf_bases.items():
        out[sample] = {}
        for t in range(lo, hi + 1):
            v = genome_len - t + 1
            raw = by_pos.get(v, "N")
            out[sample][t] = _COMPLEMENT.get(raw, "N")
    return out




def build_diversity_download_fasta(seg, lo, hi, lineage=None):
    """FASTA of every population genome's actual called sequence over reference range
    [lo, hi] (1-based inclusive, same coordinate system as the diversity chart's x-axis),
    optionally restricted to one lineage. Missing/no-call positions come back as 'N' -- this
    reads calls straight out of the pixy allsites VCF (already one column per genome, one row
    per reference position) rather than re-deriving alignment coordinates from raw genome
    FASTA files."""
    genes = {r["segment"].lower(): int(r["genome_len"]) for r in read_tsv(ROOT / "config" / "refgenome_genes.tsv")}
    genome_len = genes.get(seg)
    if not genome_len:
        return None
    lo = max(1, lo)
    hi = min(genome_len, hi)
    if lo > hi:
        return None
    meta_path = ROOT / "data" / "filtered" / f"{seg}_filtered_metadata.tsv"
    meta = {r["accession"]: r for r in read_tsv(meta_path)} if meta_path.exists() else {}
    bases_by_sample = _true_orientation_sample_bases(seg, genome_len, lo, hi)
    lines = []
    for sample in sorted(bases_by_sample):
        sample_lineage = (meta.get(sample, {}).get("lineage") or "NA")
        if lineage and lineage != pooled_label() and sample_lineage != lineage:
            continue
        by_pos = bases_by_sample[sample]
        seq = "".join(by_pos.get(p, "N") for p in range(lo, hi + 1))
        lines.append(f">{sample} lineage={sample_lineage} {seg.upper()}:{lo}-{hi}")
        lines.append(seq)
    return "\n".join(lines) + ("\n" if lines else "")


REGION_ALIGNMENT_MAX_WIDTH = 2000
REGION_ALIGNMENT_MAX_ROWS = 800


def build_region_alignment(seg, lo, hi, lineage=None):
    """Same per-genome, per-reference-position sequence reconstruction as
    build_diversity_download_fasta(), but shaped for the on-page alignment viewer (a
    detail-page-style MSA card) instead of a FASTA download: a population-majority consensus
    row plus per-genome rows already dot-diffed against it, exactly like build_msa.py's
    per-set amplicon MSA (scripts/msa_common.py's dotify()) so the same .msa-row/.msa-seq
    frontend rendering can be reused unmodified. Capped at REGION_ALIGNMENT_MAX_WIDTH columns
    -- unlike a FASTA download (which a browser or script just saves), a region this wide
    would make an unusable wall of text in the page itself."""
    genes = {r["segment"].lower(): int(r["genome_len"]) for r in read_tsv(ROOT / "config" / "refgenome_genes.tsv")}
    genome_len = genes.get(seg)
    if not genome_len:
        return {"error": "unknown segment"}
    lo = max(1, lo)
    hi = min(genome_len, hi)
    if lo > hi:
        return {"error": "invalid range"}
    if hi - lo + 1 > REGION_ALIGNMENT_MAX_WIDTH:
        return {"error": f"range too wide (max {REGION_ALIGNMENT_MAX_WIDTH}bp, requested {hi - lo + 1}bp)"}
    meta_path = ROOT / "data" / "filtered" / f"{seg}_filtered_metadata.tsv"
    meta = {r["accession"]: r for r in read_tsv(meta_path)} if meta_path.exists() else {}
    bases_by_sample = _true_orientation_sample_bases(seg, genome_len, lo, hi)

    samples = []
    for sample in sorted(bases_by_sample):
        sample_lineage = meta.get(sample, {}).get("lineage") or "NA"
        if lineage and lineage != pooled_label() and sample_lineage != lineage:
            continue
        samples.append(sample)

    positions = list(range(lo, hi + 1))
    seqs = {s: [bases_by_sample[s].get(p, "N") for p in positions] for s in samples}

    consensus_chars = []
    for i in range(len(positions)):
        counts = {}
        for s in samples:
            b = seqs[s][i]
            if b in "ACGT":
                counts[b] = counts.get(b, 0) + 1
        consensus_chars.append(max(counts.items(), key=lambda kv: kv[1])[0] if counts else "N")
    consensus = "".join(consensus_chars)

    rows = []
    for s in samples:
        seq = seqs[s]
        dots = "".join("." if seq[i] == consensus[i] else seq[i] for i in range(len(positions)))
        rows.append({
            "genome_id": s,
            "lineage": meta.get(s, {}).get("lineage") or "NA",
            "country": meta.get(s, {}).get("country", ""),
            "dots": dots,
        })
    rows.sort(key=lambda r: (r["lineage"], r["genome_id"]))
    total_count = len(rows)
    truncated = total_count > REGION_ALIGNMENT_MAX_ROWS
    if truncated:
        rows = rows[:REGION_ALIGNMENT_MAX_ROWS]

    return {
        "seg": seg, "start": lo, "end": hi, "columns": len(positions),
        "consensus": consensus, "rows": rows, "total_count": total_count,
        "truncated": truncated,
    }


def _iupac_consensus(counts_by_pos, lo, hi, min_frac):
    seq = []
    for pos in range(lo, hi + 1):
        counts = counts_by_pos.get(pos, {})
        total = sum(counts.values())
        if not total:
            seq.append("N")
            continue
        bases = frozenset(b for b, c in counts.items() if c / total >= min_frac)
        if not bases:
            bases = frozenset([max(counts, key=counts.get)])
        seq.append(_IUPAC_FROM_BASES.get(bases, "N"))
    return "".join(seq)


def _candidate_sequence_for_window(seg, genome_len, lo, hi, strand):
    counts = _true_orientation_allele_counts(seg, genome_len, lo, hi)
    forward = _iupac_consensus(counts, lo, hi, SHIFT_IUPAC_MIN_FRAC)
    return _revcomp_iupac(forward) if strand == "-" else forward


def compute_shift_suggestions(segment, binding):
    """For every found F/P/R binding site, the best nearby (+-40bp, same length) position by
    population-pooled pi, if it's a meaningful (<=50%) improvement over a non-trivial
    (>=0.02) current value -- each with the actual candidate sequence at that new position."""
    seg = segment.lower()
    diversity = load_diversity()
    genome_len = diversity["genome_len"].get(seg)
    seg_data = diversity["segments"].get(seg)
    if not genome_len or not seg_data:
        return []
    windows = seg_data.get("windows", {}).get(pooled_label())
    if not windows:
        return []

    suggestions = []
    for b in binding:
        if not b.get("found"):
            continue
        lo, hi = _true_oligo_range(b)
        current_pi = _avg_pi_for_range(windows, lo, hi)
        if current_pi is None or current_pi < SHIFT_MIN_ABS_PI:
            continue
        best = None
        for shift in range(-SHIFT_SEARCH_RANGE, SHIFT_SEARCH_RANGE + 1):
            if shift == 0:
                continue
            new_lo, new_hi = lo + shift, hi + shift
            if new_lo < 1 or new_hi > genome_len:
                continue
            pi = _avg_pi_for_range(windows, new_lo, new_hi)
            if pi is None:
                continue
            if best is None or pi < best[1]:
                best = (shift, pi)
        if best is None or best[1] > current_pi * SHIFT_MIN_IMPROVEMENT_RATIO:
            continue
        shift, new_pi = best
        new_lo, new_hi = lo + shift, hi + shift
        candidate_seq = _candidate_sequence_for_window(seg, genome_len, new_lo, new_hi, b["strand"])
        suggestions.append({
            "role": b["role"], "strand": b["strand"], "shift": shift,
            "current_pi": current_pi, "new_pi": new_pi, "candidate_seq": candidate_seq,
        })
    return suggestions


COUNTRY_NAME_ALIASES = {
    # Only aliases actually needed to join this dataset's country strings (from Nextstrain
    # metadata) against world-atlas/Natural-Earth country names used by the vendored
    # countries-110m.json topojson -- not a general-purpose country-name normalizer.
    "USA": "United States of America",
    "Cote d'Ivoire": "Côte d'Ivoire",
    "Ivory Coast": "Côte d'Ivoire",
    "DRC": "Dem. Rep. Congo",
    "Democratic Republic of the Congo": "Dem. Rep. Congo",
}
UNMAPPED_COUNTRY_LABELS = {"Unknown/Missing", "Missing", "", "NA"}


def _parse_filter_report_full(path: Path):
    if not path.exists():
        return {}
    text = path.read_text(encoding="utf-8")
    out = {}
    for key in ("kept", "dropped_short", "dropped_ambiguousN", "dropped_dup_accession", "dropped_no_metadata"):
        m = re.search(rf"{key}(?:\([^)]*\))?=(\d+)", text)
        out[key] = int(m.group(1)) if m else 0
    m = re.search(r"min_len=(\d+)", text)
    out["min_len"] = int(m.group(1)) if m else None
    return out


METADATA_TABLE_COLUMNS = [
    # key, label, type ("text" | "num")
    ("segment", "Seg", "text"), ("accession", "Accession", "text"), ("data_qc", "데이터 QC", "text"),
    ("strain", "Strain", "text"),
    ("isolate", "Isolate", "text"), ("lineage", "Lineage", "text"), ("sub_lineage", "Sub-lineage", "text"),
    ("date", "Date", "text"), ("year", "Year", "num"), ("country", "Country", "text"),
    ("region", "Region", "text"), ("division", "Division", "text"), ("location", "Location", "text"),
    ("host", "Host", "text"), ("host_class", "숙주 분류", "text"), ("length", "Length (bp)", "num"),
    ("validated_sets", "검증 통과 세트 수", "num"), ("validated_pct", "검증 통과율 (%)", "num"),
    ("diversity_qc", "Diversity QC", "text"), ("identity_to_ref", "Identity vs ref", "num"),
    ("coverage", "Coverage vs ref", "num"), ("date_released", "Released", "text"),
    ("authors", "Authors", "text"), ("institution", "Institution", "text"),
]


HOST_CLASS_LABEL = {"human": "사람", "animal": "동물", "unknown": "미상", "lab": "실험 숙주 (제외)"}

DATA_QC_LABEL = {
    "passed": "통과", "dropped_short": "제외: 길이 미달", "dropped_ambiguousN": "제외: 모호염기(N)",
    "dropped_dup_accession": "제외: 중복", "dropped_no_metadata": "제외: 메타데이터 없음",
}


def load_metadata_table():
    """One row per genome of the analysed population (data/filtered), enriched with the raw
    Nextstrain metadata columns, how many of the segment's primer sets validate the genome, and
    the diversity-QC verdict -- the data behind the sortable/filterable 메타데이터 표 page."""
    num = lambda v: None if v in (None, "") else (float(v) if "." in str(v) else int(v))  # noqa: E731
    rows = []
    for seg in seg_keys():
        raw = {r["accession"]: r for r in read_tsv(ROOT / "data" / "raw" / f"{seg}_metadata.tsv")}
        qc_path = ROOT / "output" / "diversity" / seg / "qc.json"
        qc = json.loads(qc_path.read_text(encoding="utf-8")) if qc_path.exists() else {}
        per_genome = qc.get("per_genome", {})
        low_cov = set(qc.get("excluded_low_coverage", []))
        low_id = set(qc.get("excluded_low_identity", []))
        seg_sets = [sid for sid, r in load_manifest().items() if r["segment"].lower() == seg]
        counts = {}
        for sid in seg_sets:
            try:                                              # recounted like the set cards
                for gid_ in scope_ctx().set_validated_ids(sid, seg, "all"):
                    counts[gid_] = counts.get(gid_, 0) + 1
            except (OSError, KeyError):
                pass
        filtered = {f["accession"]: f for f in read_tsv(ROOT / "data" / "filtered" / f"{seg}_filtered_metadata.tsv")}
        host_cls = scope_ctx().classes(seg)
        curated = {c["accession"]: c for c in read_tsv(ROOT / "data" / "renamed" / f"{seg}_metadata_curated.tsv")}
        # Every raw record with the outcome of the data QC filter (data/filtered/<seg>_filter_status.tsv);
        # without that file (older results) only the genomes that passed are listed.
        status_rows = read_tsv(ROOT / "data" / "filtered" / f"{seg}_filter_status.tsv") or [
            {"accession": a, "length": "", "status": "passed"} for a in filtered]
        for st in status_rows:
            acc = st["accession"]
            r = raw.get(acc, {})
            f = filtered.get(acc) or curated.get(acc) or {}
            passed = st["status"] == "passed"
            date = f.get("date") or r.get("date") or ""
            year = date[:4] if date[:4].isdigit() else ""
            g = per_genome.get(acc)
            if not passed:
                dq = "해당 없음 (데이터 QC 제외)"
            elif g:
                dq = "사용"
            elif acc in low_cov:
                dq = "제외 (coverage 미달)"
            elif acc in low_id:
                dq = "제외 (identity 미달)"
            else:
                dq = "해당 없음"
            vs = counts.get(acc, 0) if (passed and seg_sets) else None
            rows.append({
                "segment": seg.upper(), "accession": acc, "data_qc": DATA_QC_LABEL.get(st["status"], st["status"]),
                "strain": r.get("strain", ""), "isolate": r.get("isolate", ""),
                "lineage": f.get("lineage") or r.get(f"{seg}_lineage") or r.get("gpc_lineage") or "NA", "sub_lineage": r.get("sub_lineage", ""),
                "date": date, "year": num(year), "country": f.get("country") or r.get("country", ""),
                "region": f.get("region") or r.get("region", ""), "division": r.get("division", ""),
                "location": r.get("location", ""), "host": r.get("host", ""),
                "host_class": HOST_CLASS_LABEL.get(host_cls.get(acc, "unknown"), ""),
                "length": num(r.get("length")) if r.get("length") else num(st.get("length")),
                "validated_sets": vs,
                "validated_pct": round(100.0 * vs / len(seg_sets), 1) if (passed and seg_sets) else None,
                "diversity_qc": dq, "identity_to_ref": g["identity"] if g else None,
                "coverage": g["coverage"] if g else None, "date_released": r.get("date_released", ""),
                "authors": r.get("abbr_authors") or r.get("authors", ""), "institution": r.get("institution", ""),
            })
    return {"columns": [{"key": k, "label": l, "type": t} for k, l, t in METADATA_TABLE_COLUMNS], "rows": rows}


def load_metadata_stats():
    """Basic provenance/coverage statistics for the underlying Nextstrain population dataset
    (where it came from, how many genomes before/after QC filtering) plus per-country,
    per-year genome counts for the world-map + trend-line view -- everything here is derived
    straight from data/raw + data/filtered, independent of any primer set."""
    base_url = get_profile()["source"]["base_url"]

    segments = {}
    by_country = {}
    by_country_year = {}
    by_lineage_year = {}
    by_lineage = {k: {} for k in seg_keys()}
    for seg in seg_keys():
        report = _parse_filter_report_full(ROOT / "data" / "filtered" / f"{seg}_filter_report.txt")
        raw_fasta = ROOT / "data" / "raw" / f"{seg}_sequences.fasta"
        raw_total = 0
        if raw_fasta.exists():
            with open(raw_fasta, encoding="utf-8") as fh:
                raw_total = sum(1 for line in fh if line.startswith(">"))
        raw_meta_path = ROOT / "data" / "raw" / f"{seg}_metadata.tsv"
        fetched_at = None
        if raw_meta_path.exists():
            fetched_at = datetime.fromtimestamp(raw_meta_path.stat().st_mtime).isoformat(timespec="seconds")
        segments[seg] = {
            "min_len": report.get("min_len"),
            "raw_total": raw_total,
            "kept": report.get("kept", 0),
            "dropped_short": report.get("dropped_short", 0),
            "dropped_ambiguousN": report.get("dropped_ambiguousN", 0),
            "dropped_dup_accession": report.get("dropped_dup_accession", 0),
            "dropped_no_metadata": report.get("dropped_no_metadata", 0),
            "fetched_at": fetched_at,
        }

        for r in read_tsv(ROOT / "data" / "filtered" / f"{seg}_filtered_metadata.tsv"):
            country = (r.get("country") or "").strip() or "NA"
            lineage = (r.get("lineage") or "NA") or "NA"
            by_lineage[seg][lineage] = by_lineage[seg].get(lineage, 0) + 1
            slot = by_country.setdefault(country, {"s": 0, "l": 0})
            slot[seg] += 1
            date = (r.get("date") or "").strip()
            year = date[:4] if len(date) >= 4 and date[:4].isdigit() else None
            if year:
                by_country_year.setdefault(country, {}).setdefault(year, {"s": 0, "l": 0})[seg] += 1
                by_lineage_year.setdefault(lineage, {}).setdefault(year, {"s": 0, "l": 0})[seg] += 1

    by_country_rows = []
    for country, counts in by_country.items():
        mapped_name = COUNTRY_NAME_ALIASES.get(country, country)
        by_country_rows.append({
            "country": country,
            "map_name": None if country in UNMAPPED_COUNTRY_LABELS else mapped_name,
            "s": counts["s"], "l": counts["l"], "total": counts["s"] + counts["l"],
        })
    by_country_rows.sort(key=lambda r: -r["total"])

    # Same per-year s/l counts as by_country_year, just summed across every country -- the
    # population-wide "genomes reported per year" trend line (map hover already shows this
    # per-country; this is the page-level aggregate).
    by_year = {}
    for country_years in by_country_year.values():
        for year, counts in country_years.items():
            slot = by_year.setdefault(year, {"s": 0, "l": 0})
            slot["s"] += counts["s"]
            slot["l"] += counts["l"]

    return {
        "source": {"base_url": base_url},
        "segments": segments,
        "by_country": by_country_rows,
        "by_country_year": by_country_year,
        "by_year": by_year,
        "by_lineage": by_lineage,
        "by_lineage_year": by_lineage_year,
    }


MISMATCH_TABLE_INT_KEYS = ("f_mismatch", "p_mismatch", "r_mismatch", "amplicon_start", "amplicon_end",
                            "probe_start", "probe_end")
MISMATCH_TABLE_FLOAT_KEYS = ("f_tm", "f_dg", "p_tm", "p_dg", "r_tm", "r_dg")
MISMATCH_TABLE_BOOL_KEYS = ("f_3p_ok", "r_3p_ok")


def _coerce_numeric_strings(rows):
    """A user-added set's cached result.json (live_analyze.py) stores coordinates and Tm/dG as
    strings, where a built-in set's rows go through _coerce_mismatch_table_numbers() -- so the
    page did `amplicon_start + length` as string concatenation and generate_report.py raised
    TypeError. Converts the numeric columns in place (bool/flag columns are left alone)."""
    for g in rows or []:
        for k in MISMATCH_TABLE_INT_KEYS + MISMATCH_TABLE_FLOAT_KEYS:
            v = g.get(k)
            if isinstance(v, str):
                try:
                    g[k] = None if v == "" else (int(v) if k in MISMATCH_TABLE_INT_KEYS else float(v))
                except ValueError:
                    pass


def _coerce_mismatch_table_numbers(rows):
    for g in rows:
        for k in MISMATCH_TABLE_INT_KEYS:
            if g.get(k) == "":
                g[k] = None
            elif g.get(k) is not None:
                try:
                    g[k] = int(g[k])
                except ValueError:
                    pass
        for k in MISMATCH_TABLE_FLOAT_KEYS:
            if g.get(k) == "":
                g[k] = None
            elif g.get(k) is not None:
                try:
                    g[k] = float(g[k])
                except ValueError:
                    pass
        for k in MISMATCH_TABLE_BOOL_KEYS:
            # blank = not checked (three_prime_consecutive_match was 0, or this row predates
            # the column) -- kept as None rather than coerced to True/False so callers can
            # tell "passed the 3' check" apart from "3' check wasn't run".
            g[k] = None if g.get(k) in ("", None) else g.get(k) == "1"
    return rows


def load_set_detail(set_id):
    manifest = load_manifest()
    if set_id in manifest:
        sets = {s["set_id"]: s for s in load_builtin_sets()}
        base = sets.get(set_id)
        if base is None:
            return None
        genomes = _coerce_mismatch_table_numbers(
            read_tsv(ROOT / "output" / "primers" / set_id / "alignment" / "mismatch_table.tsv"))
        standard_path = ROOT / "output" / "primers" / set_id / "standard_alignment" / "mismatch_table.tsv"
        standard_genomes = _coerce_mismatch_table_numbers(read_tsv(standard_path)) if standard_path.exists() else []
        binding = read_tsv(ROOT / "output" / "primers" / set_id / "refgenome_binding" / "binding.tsv")
        for b in binding:
            b["found"] = b.get("found") == "1"
        msa_path = ROOT / "output" / "primers" / set_id / "msa" / "msa.json"
        msa = json.loads(msa_path.read_text(encoding="utf-8")) if msa_path.exists() else None
        dh_path = ROOT / "output" / "primers" / set_id / "dimer_hairpin" / "dimer_hairpin.json"
        dimer_hairpin = json.loads(dh_path.read_text(encoding="utf-8")) if dh_path.exists() else None
        base = dict(base)
        seg_k = base["segment"].lower()
        view = scope_ctx()
        base["genomes"] = view.tag_and_filter(seg_k, genomes)             # only the genomes of the active host scope
        if msa:
            msa["rows"] = view.tag_and_filter(seg_k, msa.get("rows"))
            msa["excluded"] = view.tag_and_filter(seg_k, msa.get("excluded"))
        base["host_scope"] = view.info(seg_k)
        base["standard_genomes"] = standard_genomes
        base["binding"] = binding
        base["msa"] = msa
        base["dimer_hairpin"] = dimer_hairpin
        base["standard_representativeness"] = standard_rep_for_set(set_id, base["segment"])
        base["shift_suggestions"] = compute_shift_suggestions(base["segment"], binding)
        # load_builtin_sets() already set a lightweight "evaluation" (no binding/dimer_hairpin
        # available there) -- redo it now that the full set_data is in hand, so the detail
        # view's evaluation always reflects every build_assessment() check.
        from generate_report import build_assessment  # noqa: PLC0415 -- see load_builtin_sets
        _, cls, _ = build_assessment(base)
        base["evaluation"] = {"cls": cls, "label": EVALUATION_LABEL[cls]}
        return base
    user_path = ROOT / "output" / "user_sets" / set_id / "result.json"
    if user_path.exists():
        result = json.loads(user_path.read_text(encoding="utf-8"))
        _coerce_numeric_strings(result.get("genomes"))
        _coerce_numeric_strings(result.get("standard_genomes"))
        seg_u = (result.get("segment") or "").lower()
        if result.get("genomes") and seg_u in seg_keys():
            view = scope_ctx()
            result["summary"] = {**(result.get("summary") or {}), **view.summary_from_rows(seg_u, {g["genome_id"]: g for g in result["genomes"]})}
            result["genomes"] = view.tag_and_filter(seg_u, result["genomes"])
            result["host_scope"] = view.info(seg_u)
        return result
    return None


def slugify(name: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    slug = slug.strip("_") or "user_set"
    return f"USER_{slug}"[:80]


def run_analysis_job(job_id, payload):
    def progress(msg):
        with JOBS_LOCK:
            JOBS[job_id]["progress"] = msg

    try:
        name = payload.get("name") or "user_set"
        set_id = slugify(name)
        out_dir = ROOT / "output" / "user_sets" / set_id
        params = load_params()
        result = analyze_set(
            name=set_id, segment=payload["segment"],
            f_id=f"{set_id}_F", f_seq=payload["f_seq"],
            r_id=f"{set_id}_R", r_seq=payload["r_seq"],
            p_id=f"{set_id}_P", p_seq=payload["p_seq"],
            root=ROOT, out_dir=out_dir,
            primer_mismatch_pass=params["primer_mismatch_pass"], probe_mismatch=params["probe_mismatch"],
            three_prime_n=params["three_prime_consecutive_match"], cpu=params["mfeprimer_cpu"],
            progress=progress,
        )
        result["display_name"] = name
        (out_dir / "result.json").write_text(json.dumps(result), encoding="utf-8")
        if payload.get("save_to_manifest"):
            _append_to_manifest(result)
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "done"
            JOBS[job_id]["result"] = _slim(result)
            JOBS[job_id]["set_id"] = set_id
    except Exception as exc:  # noqa: BLE001
        traceback.print_exc()
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error"
            JOBS[job_id]["error"] = str(exc)


def _append_to_manifest(result):
    """Adds a live-analyzed set to config/primers.tsv (new permanent id, proceed=True, add_date
    auto-stamped), skipping it if a primer with that name already exists."""
    from mfe_common import load_primers_manifest_raw, write_primers_manifest  # noqa: PLC0415
    manifest_path = ROOT / "config" / "primers.tsv"
    rows = load_primers_manifest_raw(manifest_path)
    name = result.get("display_name") or result["set_id"]
    if any(r.get("primer_name") == name for r in rows):
        return
    rows.append({
        "primer_id": "", "primer_name": name, "segment": str(result["segment"]).upper(),
        "forward_primer": result["F_seq"], "reverse_primer": result["R_seq"], "probe": result["P_seq"],
        "fluorophore": result.get("fluorophore", ""), "source": "user", "proceed": "True", "add_date": "",
    })
    write_primers_manifest(manifest_path, rows)
    sync_pipeline_inputs()


class Handler(BaseHTTPRequestHandler):
    server_version = "LASVDashboard/1.0"

    def log_message(self, fmt, *args):  # noqa: A003
        sys.stderr.write(f"[webapp] {self.address_string()} {fmt % args}\n")

    def _read_json_body(self):
        length = int(self.headers.get("Content-Length", 0))
        return json.loads(self.rfile.read(length) or b"{}")

    def _send_json(self, obj, status=200):
        body = json.dumps(obj).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_file(self, path: Path, content_type: str, download_name: str = None):
        if not path.exists():
            self.send_error(404)
            return
        body = path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        # revalidate on every load: without this browsers keep an old app.js/style.css after an update (new buttons "missing")
        self.send_header("Cache-Control", "no-cache")
        if download_name:
            self.send_header("Content-Disposition", f'attachment; filename="{download_name}"')
        self.end_headers()
        self.wfile.write(body)

    def _send_bytes(self, body: bytes, content_type: str, download_name: str = None):
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        if download_name:
            self.send_header("Content-Disposition", f'attachment; filename="{download_name}"')
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):  # noqa: N802
        self._guarded(self._dispatch_get)

    def do_POST(self):  # noqa: N802
        self._guarded(self._dispatch_post)

    def _guarded(self, dispatch):
        """Any unhandled exception inside a request handler otherwise kills this thread
        silently -- BaseHTTPRequestHandler has no built-in try/except around do_GET/do_POST,
        so the client just sees the TCP connection drop with no response at all (e.g.
        Python's urllib surfaces this as http.client.RemoteDisconnected), and the actual
        traceback only ever reached stderr, invisible unless someone was watching the
        terminal at that exact moment. Catches everything, logs the full traceback to a
        file that survives (not just stderr), and still returns a real HTTP response so the
        caller gets an actual error to report/debug instead of a mysterious disconnect."""
        try:
            dispatch()
        except (Exception, SystemExit) as exc:  # noqa: BLE001
            # SystemExit too: helper modules shared with the Snakemake scripts report fatal
            # problems with sys.exit("[error] ...") (e.g. mfeprimer not found), and that is NOT an
            # Exception -- uncaught it silently ends this handler thread and drops the connection.
            tb = traceback.format_exc()
            sys.stderr.write(tb)
            try:
                log_path = ROOT / "logs" / "webapp_server_errors.log"
                log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(log_path, "a", encoding="utf-8") as fh:
                    fh.write(f"\n--- {datetime.now().isoformat(timespec='seconds')} {self.command} {self.path} ---\n{tb}")
            except Exception:  # noqa: BLE001
                pass  # logging the error must never itself crash the handler
            try:
                self._send_json({"error": f"internal server error: {exc}"}, 500)
            except Exception:  # noqa: BLE001
                pass  # headers/body may already be partially sent -- nothing more we can do

    def _dispatch_get(self):
        parsed = urlparse(self.path)
        path = parsed.path
        qs = parse_qs(parsed.query)
        scope_view.set_request_scope((qs.get("host") or [None])[0])
        if path == "/" or path == "/index.html":
            return self._send_file(STATIC_DIR / "index.html", "text/html; charset=utf-8")
        if path == "/app.js":
            # Serves the esbuild-minified bundle when present (see webapp/package.json's
            # "build" script) and falls back to the raw source otherwise -- so index.html
            # never needs to know which one exists, and the app still works unminified on a
            # machine without Node/esbuild installed.
            minified = STATIC_DIR / "app.min.js"
            src = minified if minified.exists() else STATIC_DIR / "app.js"
            return self._send_file(src, "application/javascript; charset=utf-8")
        if path == "/style.css":
            return self._send_file(STATIC_DIR / "style.css", "text/css; charset=utf-8")
        if path == "/publish.js":                       # "Firebase 업로드" button (see firebase/ in the project folder above)
            return self._send_file(STATIC_DIR / "publish.js", "application/javascript; charset=utf-8")
        if path == "/api/firebase/config":              # absent file -> no button; contains only the public web config
            cfg_file = ROOT / "config" / "firebase.json"
            if not cfg_file.exists():
                return self._send_json({"error": "not configured"}, 404)
            return self._send_json(json.loads(cfg_file.read_text(encoding="utf-8")))
        if path == "/lumo.css":
            return self._send_file(STATIC_DIR / "lumo.css", "text/css; charset=utf-8")

        if path == "/api/report/download":
            return self._send_file(ROOT / "output" / "report" / "report.pdf", "application/pdf",
                                    download_name=f"{get_profile()['short_name']}_primer_validation_report.pdf")

        m = re.match(r"^/vendor/([A-Za-z0-9_.-]+)$", path)
        if m:
            fname = m.group(1)
            if fname.endswith(".js"):
                content_type = "application/javascript; charset=utf-8"
            elif fname.endswith(".html"):
                content_type = "text/html; charset=utf-8"
            elif fname.endswith(".json"):
                content_type = "application/json; charset=utf-8"
            else:
                content_type = "application/octet-stream"
            return self._send_file(STATIC_DIR / "vendor" / fname, content_type)

        if path == "/api/sets":
            return self._send_json(load_builtin_sets())

        if path == "/api/host-scope":
            return self._send_json({"scopes": ["human", "animal", "unknown", "all"], "segments": host_scope_info()})
        if path == "/api/overview":
            return self._send_json(load_overview())
        if path == "/api/final-results":
            return self._send_json(load_final_results())
        if path == "/api/venn-saved":
            return self._send_json(load_venn_saved())
        if path == "/api/pipelines":
            return self._send_json(load_pipelines())
        m = re.match(r"^/api/venn/(\w+)$", path)
        if m:
            seg = m.group(1).lower()
            if seg not in seg_keys():
                return self._send_json({"error": "bad segment"}, 400)
            return self._send_json(load_venn(seg))

        if path == "/api/requirements":
            f = ROOT / "docs" / "requirements_log.json"
            return self._send_json(json.loads(f.read_text(encoding="utf-8")) if f.exists() else {"phases": [], "open": []})

        if path == "/api/profile":
            try:
                version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
            except OSError:
                version = ""
            return self._send_json(dict(get_profile(), version=version))        # the footer shows the VERSION file, not a typed-in number

        if path == "/api/diversity":
            return self._send_json(load_diversity())

        if path == "/api/standards/manifest":
            from standards_common import STANDARD_COLUMNS, ensure_standards_manifest  # noqa: PLC0415
            rows = ensure_standards_manifest(ROOT / "config" / "standards.tsv", ROOT / "config", seg_keys())
            for r in rows:
                r["length"] = len((r.get("sequence") or "").strip())
            return self._send_json({"columns": STANDARD_COLUMNS, "rows": rows})

        if path == "/api/standards/representativeness":
            return self._send_json(load_standard_representativeness())

        if path == "/api/standards/unmatched":
            seg = (qs.get("seg", [""])[0] or "").lower()
            set_id = qs.get("set", [None])[0]
            if seg not in seg_keys():
                return self._send_json({"error": "seg must be s or l"}, status=400)
            text = standard_rep_unmatched_tsv(seg, set_id)
            if text is None:
                return self._send_json({"error": "not computed yet"}, status=404)
            name = f"standard_unmatched_{set_id or seg.upper()}.tsv"
            return self._send_bytes(text.encode("utf-8"), "text/tab-separated-values; charset=utf-8", download_name=name)

        if path == "/api/diversity/download":
            seg = (qs.get("seg", [""])[0] or "").lower()
            lineage = qs.get("lineage", [None])[0]
            if seg not in seg_keys():
                return self._send_json({"error": "seg must be s or l"}, status=400)
            try:
                lo = int(qs.get("start", [""])[0])
                hi = int(qs.get("end", [""])[0])
            except (TypeError, ValueError):
                return self._send_json({"error": "start/end must be integers"}, status=400)
            fasta = build_diversity_download_fasta(seg, lo, hi, lineage)
            if fasta is None:
                return self._send_json({"error": "invalid range"}, status=400)
            name = f"{get_profile()['short_name']}_{seg.upper()}_{lo}-{hi}{'_' + lineage if lineage and lineage != pooled_label() else ''}.fasta"
            return self._send_bytes(fasta.encode("utf-8"), "text/plain; charset=utf-8", download_name=name)

        if path == "/api/diversity/region":
            seg = (qs.get("seg", [""])[0] or "").lower()
            lineage = qs.get("lineage", [None])[0]
            if seg not in seg_keys():
                return self._send_json({"error": "seg must be s or l"}, status=400)
            try:
                lo = int(qs.get("start", [""])[0])
                hi = int(qs.get("end", [""])[0])
            except (TypeError, ValueError):
                return self._send_json({"error": "start/end must be integers"}, status=400)
            result = build_region_alignment(seg, lo, hi, lineage)
            status = 400 if "error" in result else 200
            return self._send_json(result, status=status)

        if path == "/api/metadata":
            return self._send_json(load_metadata_stats())

        if path == "/api/metadata/table":
            return self._send_json(load_metadata_table())

        if path == "/api/pipeline/freshness":
            # The page asks once on every load (lazily, after it has rendered), so remember the
            # answer for a few minutes -- several people opening the dashboard shouldn't each
            # fire four HEAD/probe requests at Nextstrain. A finished 데이터 업데이트 rewrites
            # the baseline (new checked_at), which changes the cache key and so invalidates it.
            baseline_key = (load_fetch_state() or {}).get("checked_at")
            cached = _FRESHNESS_CACHE.get("entry")
            if cached and cached["key"] == baseline_key and time.time() - cached["t"] < 300:
                return self._send_json(cached["payload"])
            try:
                remote = probe_nextstrain_files()
            except Exception as exc:  # noqa: BLE001
                return self._send_json({"error": f"Nextstrain 서버 확인 실패: {exc}"}, 502)
            stored = (load_fetch_state() or {}).get("remote") or {}
            seg_label = {"s": "S segment", "l": "L segment"}
            kind_label = {"sequences": "서열", "metadata": "메타데이터"}
            changed_files, unchanged_files = [], []
            for seg in seg_keys():
                for kind in ("sequences", "metadata"):
                    cur = remote.get(seg, {}).get(kind, {})
                    prev = stored.get(seg, {}).get(kind, {})
                    label = f"{seg_label[seg]} {kind_label[kind]}"
                    entry = {"label": label, "record_count": cur.get("record_count")}
                    if not prev:
                        # No prior baseline (first check, or state file predates this
                        # feature) -- can't say "unchanged" with nothing to compare
                        # against, so fail open and let the update proceed normally.
                        changed_files.append(dict(entry, reason="no_baseline"))
                    elif _fetch_fingerprint(cur) != _fetch_fingerprint(prev):
                        changed_files.append(dict(entry, reason="content_changed",
                                                   prev_record_count=prev.get("record_count")))
                    else:
                        unchanged_files.append(entry)
            payload = {
                "changed": len(changed_files) > 0,
                "changed_files": changed_files,
                "unchanged_files": unchanged_files,
                "checked_at": (load_fetch_state() or {}).get("checked_at"),
            }
            _FRESHNESS_CACHE["entry"] = {"key": baseline_key, "t": time.time(), "payload": payload}
            return self._send_json(payload)

        if path == "/api/groups":
            return self._send_json(load_group_list())

        if path == "/api/ui-state":
            return self._send_json(load_ui_state())

        if path in ("/api/primers/xlsx", "/api/standards/xlsx"):
            kind = "primers" if "primers" in path else "standards"
            body = xlsx_export(kind)
            self.send_response(200)
            self.send_header("Content-Type", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            self.send_header("Content-Disposition", f'attachment; filename="{PROFILE_ID_FOR_FILES()}_{kind}.xlsx"')
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return None

        if path == "/api/params":
            return self._send_json(load_params())

        if path == "/api/params/check-mfeprimer":
            from mfe_common import check_mfeprimer  # noqa: PLC0415
            q = parse_qs(urlparse(self.path).query)
            return self._send_json(check_mfeprimer((q.get("path") or [""])[0].strip()))

        if path == "/api/changelog":
            path_obj = ROOT / "CHANGELOG.md"
            text = path_obj.read_text(encoding="utf-8") if path_obj.exists() else "# Changelog\n\n(아직 없음)"
            return self._send_json({"markdown": text})

        if path == "/api/literature":
            return self._send_json(literature_payload())
        if path == "/api/primers/manifest":
            from mfe_common import load_primers_manifest_raw, PRIMER_MANIFEST_COLUMNS  # noqa: PLC0415
            rows = load_primers_manifest_raw(ROOT / "config" / "primers.tsv")
            return self._send_json({"columns": PRIMER_MANIFEST_COLUMNS, "rows": rows})

        m = re.match(r"^/api/groups/([^/]+)$", path)
        if m:
            detail = load_group_detail(m.group(1))
            if detail is None:
                return self._send_json({"error": "not found"}, 404)
            return self._send_json(detail)

        m = re.match(r"^/api/auspice/(s|l)$", path)
        if m:
            try:
                raw = fetch_auspice(m.group(1))
            except Exception as exc:  # noqa: BLE001
                return self._send_json({"error": f"auspice fetch failed: {exc}"}, 502)
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(raw)))
            self.end_headers()
            self.wfile.write(raw)
            return None

        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/note/images/([A-Za-z0-9]+)$", path)
        if m:
            return self._handle_get_note_image(m.group(1), m.group(2))

        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/verify/([A-Za-z0-9]+)/images/([A-Za-z0-9]+)$", path)
        if m:
            return self._handle_get_verification_image(m.group(1), m.group(2), m.group(3))

        m = re.match(r"^/api/redesign/(\w+)/run$", path)
        if m:
            seg = m.group(1).lower()
            if seg not in seg_keys():
                return self._send_json({"error": "bad segment"}, 400)
            scope, _ids = redesign_scope(seg, (qs.get("host") or [None])[0])
            engine = "primerhunter" if (qs.get("engine") or [""])[0] == "primerhunter" else "varvamp"
            f = _rd_dir(seg, scope, engine) / "run.json"
            out = json.loads(f.read_text(encoding="utf-8")) if f.exists() else {"generated": None}
            out["host_scope"] = out.get("host_scope") or scope
            counts = {}
            lin = scope_ctx().lineage_of(seg)
            for i in _ids:
                counts[lin.get(i, "NA")] = counts.get(lin.get(i, "NA"), 0) + 1
            out["available"] = dict(sorted(counts.items()))        # lineage -> genomes in the active host scope
            out["active_scope"] = scope
            return self._send_json(out)

        m = re.match(r"^/api/redesign/(\w+)/msa$", path)
        if m:
            seg = m.group(1).lower()
            rel = (qs.get("dir") or [""])[0]
            if seg not in seg_keys() or not re.fullmatch(r"[A-Za-z0-9_\-]+/designs/[A-Za-z0-9_\-]+", rel):
                return self._send_json({"error": "bad request"}, 400)
            scope, _ids = redesign_scope(seg, (qs.get("host") or [None])[0])
            engine = "primerhunter" if (qs.get("engine") or [""])[0] == "primerhunter" else "varvamp"
            f = _rd_dir(seg, scope, engine) / rel / "result.json"
            if not f.exists():
                return self._send_json({"error": "결과가 없습니다."}, 404)
            data = json.loads(f.read_text(encoding="utf-8"))
            v = scope_ctx()
            msa = dict(data.get("msa") or {})
            msa["rows"] = v.tag_and_filter(seg, msa.get("rows") or [])
            genomes = v.tag_and_filter(seg, data.get("genomes") or [])
            _coerce_numeric_strings(msa["rows"])
            _coerce_numeric_strings(genomes)
            return self._send_json({"set_id": rel, "segment": data.get("segment"), "msa": msa, "genomes": genomes,
                                    "gene": data.get("gene"), "summary": data.get("summary"), "host_scope": v.info(seg)})

        m = re.match(r"^/api/sets/([^/]+)/probe-design/msa$", path)
        if m:
            seq = (qs.get("seq") or [""])[0].strip().upper()
            f = ROOT / "output" / "probe_design" / m.group(1) / "validation" / (_pd_slug(seq) if _PD_SEQ_RE.match(seq) else "-") / "result.json"
            if not _PD_ID_RE.fullmatch(m.group(1)) or not f.exists():
                return self._send_json({"error": "이 Probe의 검증 결과가 없습니다. 먼저 \"검증\"을 실행하세요."}, 404)
            data = json.loads(f.read_text(encoding="utf-8"))
            seg = (data.get("segment") or "").lower()
            v = scope_ctx()
            msa = dict(data.get("msa") or {})
            msa["rows"] = v.tag_and_filter(seg, msa.get("rows") or [])         # only the genomes of the active host scope
            genomes = v.tag_and_filter(seg, data.get("genomes") or [])
            _coerce_numeric_strings(msa["rows"])
            _coerce_numeric_strings(genomes)
            return self._send_json({"set_id": m.group(1), "segment": data.get("segment"), "msa": msa, "genomes": genomes,
                                    "gene": data.get("gene"), "summary": data.get("summary"), "P_seq": seq,
                                    "host_scope": v.info(seg)})

        m = re.match(r"^/api/sets/([^/]+)/probe-design/best$", path)
        if m:
            if not _PD_ID_RE.fullmatch(m.group(1)):
                return self._send_json({"error": "bad id"}, 400)
            f = ROOT / "output" / "probe_design" / m.group(1) / "best.json"
            return self._send_json(json.loads(f.read_text(encoding="utf-8")) if f.exists() else {"generated": None})

        m = re.match(r"^/api/sets/([^/]+)/probe-design/validations$", path)
        if m:
            if not _PD_ID_RE.fullmatch(m.group(1)):
                return self._send_json({"error": "bad id"}, 400)
            items = []
            for f in (ROOT / "output" / "probe_design" / m.group(1) / "validation").glob("*/summary.json"):
                try:
                    items.append(json.loads(f.read_text(encoding="utf-8")))
                except ValueError:
                    pass
            items.sort(key=lambda x: x.get("generated", ""), reverse=True)
            return self._send_json(items)

        m = re.match(r"^/api/sets/([^/]+)/probe-design$", path)
        if m:
            f = ROOT / "output" / "probe_design" / m.group(1) / "result.json"
            if not re.fullmatch(r"[A-Za-z0-9_\-]+", m.group(1)):
                return self._send_json({"error": "bad id"}, 400)
            return self._send_json(json.loads(f.read_text(encoding="utf-8")) if f.exists() else {"set_id": m.group(1), "generated": None})

        m = re.match(r"^/api/sets/([^/]+)$", path)
        if m:
            detail = load_set_detail(m.group(1))
            if detail is None:
                return self._send_json({"error": "not found"}, 404)
            detail = dict(detail)
            detail["lineage_totals"] = load_lineage_totals().get(detail.get("segment"), {})
            detail["population_ids"] = load_population_ids(detail.get("segment", ""))
            return self._send_json(detail)

        if path == "/api/jobs/active":
            # Polled by every connected browser (not just the one that clicked a button) so
            # anyone looking at the dashboard from a different machine can see that a job is
            # already running elsewhere -- this pipeline/webapp has no concept of "users",
            # so this is the only way a second person avoids stepping on an in-progress run.
            with JOBS_LOCK:
                active = [dict(j, job_id=jid) for jid, j in JOBS.items() if j.get("status") == "running"]
            for j in active:
                j["elapsed_sec"] = round(time.time() - j["started"]) if j.get("started") else None
                j.pop("log_tail", None)
            return self._send_json(active)

        m = re.match(r"^/api/jobs/([^/]+)$", path)
        if m:
            with JOBS_LOCK:
                job = JOBS.get(m.group(1))
                if job is not None:
                    job = dict(job)  # shallow copy so we can safely read outside the lock
            if job is None:
                return self._send_json({"error": "not found"}, 404)
            if job.get("kind") == "pipeline" and job.get("log_path"):
                job["log_tail"] = tail_log_lines(Path(job["log_path"]))
                if job.get("status") == "running":
                    job["progress"] = tail_progress(Path(job["log_path"]))
                    job["elapsed_sec"] = round(time.time() - job["started"])
            return self._send_json(job)

        self.send_error(404)

    def _dispatch_post(self):
        if self.path == "/api/analyze":
            return self._handle_analyze()
        if self.path == "/api/pipeline/fetch-data":
            return self._handle_pipeline(
                label="데이터 업데이트 (nextstrain fetch + mfeprimer indexing)",
                args=["--forcerun", "fetch_raw", "--",
                      "data/filtered/s_filtered.fasta.fai", "data/filtered/l_filtered.fasta.fai"],
                on_success=save_fetch_baseline_from_remote,
            )
        m = re.match(r"^/api/xlsx/validate\?kind=(primers|standards)$", self.path)
        if m:
            try:
                payload = self._read_json_body()
            except json.JSONDecodeError:
                return self._send_json({"error": "invalid JSON"}, 400)
            rows = payload.get("rows")
            if not isinstance(rows, list):
                return self._send_json({"error": "rows는 배열이어야 합니다"}, 400)
            return self._send_json({"errors": xlsx_validate(m.group(1), rows)})
        m = re.match(r"^/api/xlsx/parse\?kind=(primers|standards)$", self.path)
        if m:
            length = int(self.headers.get("Content-Length") or 0)
            if length <= 0 or length > 8 * 1024 * 1024:
                return self._send_json({"error": "파일이 비었거나 너무 큽니다 (최대 8MB)."}, 400)
            try:
                return self._send_json(xlsx_parse(m.group(1), self.rfile.read(length)))
            except ValueError as exc:
                return self._send_json({"error": str(exc)}, 400)
        m = re.match(r"^/api/groups/([^/?]+)/analyze(?:\?(.*))?$", self.path)
        if m:
            return self._handle_group_analyze(m.group(1), "force=1" in (m.group(2) or ""))
        if self.path == "/api/redesign":
            return self._handle_redesign()
        m = re.match(r"^/api/sets/([^/?]+)/probe-design/best$", self.path)
        if m:
            return self._handle_probe_best(m.group(1))
        m = re.match(r"^/api/sets/([^/?]+)/probe-design/validate$", self.path)
        if m:
            return self._handle_probe_validate(m.group(1))
        m = re.match(r"^/api/sets/([^/?]+)/probe-design$", self.path)
        if m:
            return self._handle_probe_design(m.group(1))
        m = re.match(r"^/api/sets/([^/?]+)/analyze(?:\?(.*))?$", self.path)
        if m:
            return self._handle_set_analyze(m.group(1), "force=1" in (m.group(2) or ""))
        if self.path == "/api/pipeline/rerun-analysis":
            return self._handle_pipeline(label="전체 분석 다시 실행", args=[])
        if self.path == "/api/report/generate":
            return self._handle_generate_report()
        if self.path == "/api/groups":
            return self._handle_create_group()
        if self.path == "/api/venn-saved":
            try:
                body = self._read_json_body()
            except json.JSONDecodeError:
                return self._send_json({"error": "invalid JSON"}, 400)
            name = str(body.get("name", "")).strip()[:80]
            keys = [k for k in (body.get("keys") or []) if isinstance(k, str) and re.fullmatch(r"[A-Za-z0-9_\-]+", k)][:12]
            seg = str(body.get("seg", "")).lower()
            if not name or not keys or seg not in seg_keys():
                return self._send_json({"error": "이름과 한 개 이상의 세트·그룹이 필요합니다."}, 400)
            items = [x for x in load_venn_saved() if not (x.get("name") == name and x.get("seg") == seg)]    # same name = replace
            item = {"id": uuid.uuid4().hex[:8], "name": name, "seg": seg, "keys": keys, "lineage": str(body.get("lineage", "ALL"))[:40],
                    "created": time.strftime("%Y-%m-%d %H:%M")}
            items.append(item)
            _write_venn_saved(items)
            return self._send_json(item)
        m = re.match(r"^/api/venn-saved/([0-9a-f]{8})/delete$", self.path)
        if m:
            _write_venn_saved([x for x in load_venn_saved() if x.get("id") != m.group(1)])
            return self._send_json({"ok": True})
        m = re.match(r"^/api/set-parent/([A-Za-z0-9_\-]+)$", self.path)
        if m:
            try:
                payload = self._read_json_body()
            except json.JSONDecodeError:
                return self._send_json({"error": "invalid JSON"}, 400)
            parent = str(payload.get("parent", "")).strip()
            sets = {x["set_id"]: x for x in load_builtin_sets()}
            if m.group(1) not in sets or (parent and (parent not in sets or parent == m.group(1) or sets[parent].get("segment") != sets[m.group(1)].get("segment"))):
                return self._send_json({"error": "부모 세트가 올바르지 않습니다 (같은 segment의 다른 세트여야 합니다)."}, 400)
            if parent and sets[parent].get("parent_id") == m.group(1):
                return self._send_json({"error": "서로를 부모로 지정할 수 없습니다."}, 400)
            save_set_parent(m.group(1), parent)
            return self._send_json({"ok": True})
        m = re.match(r"^/api/final-results/([A-Za-z0-9_\-]+)$", self.path)
        if m:
            try:
                payload = self._read_json_body()
            except json.JSONDecodeError:
                return self._send_json({"error": "invalid JSON"}, 400)
            return self._send_json({"id": m.group(1), "result": save_final_result(m.group(1), str(payload.get("verdict", "")), str(payload.get("note", "")))})
        if self.path == "/api/ui-state/order":
            return self._handle_save_order()
        if self.path == "/api/params":
            try:
                payload = self._read_json_body()
            except json.JSONDecodeError:
                return self._send_json({"error": "invalid JSON"}, 400)
            try:
                save_params(payload)
            except (TypeError, ValueError) as exc:
                return self._send_json({"error": f"invalid parameter value: {exc}"}, 400)
            return self._send_json(load_params())
        if self.path == "/api/standards/manifest":
            return self._handle_save_standards()
        if self.path in ("/api/literature/save", "/api/literature/add"):
            try:
                payload = self._read_json_body()
            except json.JSONDecodeError:
                return self._send_json({"error": "invalid JSON"}, 400)
            virus = literature_edition()
            if virus is None or payload.get("virus", virus) != virus:
                return self._send_json({"error": "이 대시보드의 바이러스 문헌만 다룰 수 있습니다"}, 400)
            if self.path.endswith("/save"):
                rows = payload.get("rows")
                if not isinstance(rows, list):
                    return self._send_json({"error": "rows는 배열이어야 합니다"}, 400)
                rows = literature.normalise(rows, virus)
                errors = literature.validate(rows, virus)
                if errors:
                    return self._send_json({"errors": errors}, 400)
                literature.save(ROOT, virus, rows)
                return self._send_json({"ok": True, "payload": literature_payload()})
            # add one entry to this edition's analysis list (config/primers.tsv)
            from mfe_common import (  # noqa: PLC0415
                PRIMER_MANIFEST_COLUMNS, load_primers_manifest_raw, validate_primer_manifest_rows, write_primers_manifest,
            )
            row = next((r for r in literature.load(ROOT, virus) if r["lit_id"] == payload.get("lit_id")), None)
            if row is None:
                return self._send_json({"error": "해당 항목을 찾을 수 없습니다"}, 404)
            if not literature.complete(row):
                return self._send_json({"error": "Forward/Reverse/Probe 서열이 모두 있어야 분석에 추가할 수 있습니다"}, 400)
            manifest = load_primers_manifest_raw(ROOT / "config" / "primers.tsv")
            if literature.same_sequences(manifest, row):
                return self._send_json({"error": "같은 서열의 세트가 이미 분석 목록에 있습니다"}, 400)
            new = literature.manifest_row(row)
            names = {m.get("primer_name") for m in manifest}
            base, n = new["primer_name"], 2
            while new["primer_name"] in names:
                new["primer_name"], n = f"{base} ({n})", n + 1
            rows = [{c: str(m.get(c, "")) for c in PRIMER_MANIFEST_COLUMNS} for m in manifest]
            rows.append({c: str(new.get(c, "")) for c in PRIMER_MANIFEST_COLUMNS})
            errors = validate_primer_manifest_rows(rows, segments=tuple(k.upper() for k in seg_keys()))
            if errors:
                return self._send_json({"errors": errors}, 400)
            saved = write_primers_manifest(ROOT / "config" / "primers.tsv", rows)
            sync_pipeline_inputs()
            return self._send_json({"ok": True, "name": new["primer_name"], "rows": len(saved),
                                    "note": "분석 목록에 추가했습니다. 사이드바의 '분석 다시하기'를 눌러야 결과가 계산됩니다."})
        if self.path == "/api/primers/manifest":
            from mfe_common import (  # noqa: PLC0415 -- see load_builtin_sets for why
                PRIMER_MANIFEST_COLUMNS, validate_primer_manifest_rows, write_primers_manifest,
            )
            try:
                payload = self._read_json_body()
            except json.JSONDecodeError:
                return self._send_json({"error": "invalid JSON"}, 400)
            rows = payload.get("rows")
            if not isinstance(rows, list):
                return self._send_json({"error": "rows는 배열이어야 합니다"}, 400)
            rows = [{col: str(r.get(col, "")) for col in PRIMER_MANIFEST_COLUMNS} for r in rows]
            errors = validate_primer_manifest_rows(rows, segments=tuple(k.upper() for k in seg_keys()))
            if errors:
                return self._send_json({"errors": errors}, 400)
            # The spreadsheet UI always keeps one trailing blank row ready to type a new
            # entry into -- validate_primer_manifest_rows() above already treats a fully
            # empty row as "not an error", so drop it here rather than persisting a
            # near-empty row to the file on every single save.
            rows = [r for r in rows if r["primer_name"].strip() or r["forward_primer"].strip()
                    or r["reverse_primer"].strip() or r["probe"].strip()]
            saved = write_primers_manifest(ROOT / "config" / "primers.tsv", rows)
            sync_pipeline_inputs()
            # a set whose "분석진행" is off is gone from the dashboard at once: its results are deleted right here
            off = [r["primer_id"] for r in saved if r.get("primer_id") and str(r.get("proceed", "")).strip().lower() not in ("true", "1", "yes", "y")]
            purged = purge_set_results(off)
            return self._send_json({"ok": True, "rows": saved, "purged": purged})
        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/note$", self.path)
        if m:
            return self._handle_save_note(m.group(1))
        m = re.match(r"^/api/sets/([^/]+)/rename$", self.path)
        if m:
            return self._handle_rename_set(m.group(1))
        m = re.match(r"^/api/groups/([^/]+)/rename$", self.path)
        if m:
            return self._handle_rename_group(m.group(1))
        m = re.match(r"^/api/sets/([^/]+)/fluorophore$", self.path)
        if m:
            return self._handle_set_fluorophore(m.group(1))
        m = re.match(r"^/api/sets/([^/]+)/citation$", self.path)
        if m:
            return self._handle_set_citation(m.group(1))
        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/note/images$", self.path)
        if m:
            return self._handle_upload_note_image(m.group(1))
        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/verify$", self.path)
        if m:
            return self._handle_create_verification(m.group(1))
        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/verify/([A-Za-z0-9]+)/images$", self.path)
        if m:
            return self._handle_upload_verification_image(m.group(1), m.group(2))
        return self.send_error(404)

    def _handle_save_note(self, set_id):
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        if not _entity_exists(set_id):
            return self._send_json({"error": "not found"}, 404)
        state = load_ui_state()
        note = str(payload.get("note", ""))
        if note:
            state["notes"][set_id] = note
        else:
            state["notes"].pop(set_id, None)
        save_ui_state(state)
        return self._send_json({"ok": True})

    NOTE_IMAGE_TYPES = {
        "image/png": "png", "image/jpeg": "jpg", "image/jpg": "jpg",
        "image/gif": "gif", "image/webp": "webp",
    }
    NOTE_IMAGE_MAX_BYTES = 8 * 1024 * 1024  # 8MB per image -- generous for a gel photo

    def _handle_upload_note_image(self, set_id):
        if not _entity_exists(set_id):
            return self._send_json({"error": "not found"}, 404)
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        content_type = str(payload.get("content_type", "")).lower()
        ext = self.NOTE_IMAGE_TYPES.get(content_type)
        if not ext:
            return self._send_json({"error": "지원하지 않는 이미지 형식입니다 (png/jpg/gif/webp만 가능)"}, 400)
        data_b64 = payload.get("data_base64", "")
        try:
            raw = base64.b64decode(data_b64, validate=True)
        except (ValueError, TypeError):
            return self._send_json({"error": "이미지 데이터를 읽을 수 없습니다"}, 400)
        if not raw or len(raw) > self.NOTE_IMAGE_MAX_BYTES:
            return self._send_json({"error": "이미지 크기가 너무 크거나 비어 있습니다 (8MB 이하만 가능)"}, 400)

        image_id = uuid.uuid4().hex[:16]
        img_dir = ROOT / "config" / "note_images" / set_id
        img_dir.mkdir(parents=True, exist_ok=True)
        (img_dir / f"{image_id}.{ext}").write_bytes(raw)

        state = load_ui_state()
        state["note_images"].setdefault(set_id, [])
        entry = {
            "id": image_id, "ext": ext,
            "caption": str(payload.get("caption", ""))[:200],
            "uploaded": time.strftime("%Y-%m-%d %H:%M"),
        }
        state["note_images"][set_id].append(entry)
        save_ui_state(state)
        return self._send_json({"ok": True, "image": entry})

    def _handle_get_note_image(self, set_id, image_id):
        state = load_ui_state()
        entries = state.get("note_images", {}).get(set_id, [])
        entry = next((e for e in entries if e["id"] == image_id), None)
        if entry is None:
            return self._send_json({"error": "not found"}, 404)
        img_path = ROOT / "config" / "note_images" / set_id / f"{image_id}.{entry['ext']}"
        content_type = {v: k for k, v in self.NOTE_IMAGE_TYPES.items()}.get(entry["ext"], "application/octet-stream")
        return self._send_file(img_path, content_type)

    def _handle_delete_note_image(self, set_id, image_id):
        state = load_ui_state()
        entries = state.get("note_images", {}).get(set_id, [])
        entry = next((e for e in entries if e["id"] == image_id), None)
        if entry is None:
            return self._send_json({"error": "not found"}, 404)
        img_path = ROOT / "config" / "note_images" / set_id / f"{image_id}.{entry['ext']}"
        img_path.unlink(missing_ok=True)
        state["note_images"][set_id] = [e for e in entries if e["id"] != image_id]
        save_ui_state(state)
        return self._send_json({"ok": True})

    def _handle_rename_set(self, set_id):
        # The primer name is just a label -- the permanent id (P001...) is what results, notes,
        # wet-lab records and groups hang off, and the name isn't part of input_primers/<id>.tsv
        # -- so renaming here edits the spreadsheet's primer_name only and recomputes nothing.
        from mfe_common import load_primers_manifest_raw, write_primers_manifest  # noqa: PLC0415
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        manifest_path = ROOT / "config" / "primers.tsv"
        rows = load_primers_manifest_raw(manifest_path)
        row = next((r for r in rows if r.get("primer_id") == set_id), None)
        if row is None:
            return self._send_json({"error": "not found"}, 404)
        name = str(payload.get("name", "")).strip()[:120]
        if not name:
            return self._send_json({"error": "이름은 비울 수 없습니다"}, 400)
        if any(r.get("primer_name") == name and r is not row for r in rows):
            return self._send_json({"error": f"'{name}' 이름이 이미 다른 프라이머에 사용 중입니다"}, 400)
        row["primer_name"] = name
        write_primers_manifest(manifest_path, rows)
        return self._send_json({"ok": True, "display_name": name})

    def _handle_set_fluorophore(self, set_id):
        # Writes straight into config/primers.tsv's own fluorophore column now (the
        # manifest is the single source of truth since the "프라이머 입력/관리" spreadsheet
        # page can also edit this same field) rather than the old separate ui_state.json
        # side-table, so the two editing paths can never disagree with each other.
        from mfe_common import load_primers_manifest_raw, write_primers_manifest  # noqa: PLC0415
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        manifest_path = ROOT / "config" / "primers.tsv"
        rows = load_primers_manifest_raw(manifest_path)
        row = next((r for r in rows if r.get("primer_id") == set_id), None)
        if row is None:
            return self._send_json({"error": "not found"}, 404)
        value = str(payload.get("fluorophore", "")).strip()[:40]
        row["fluorophore"] = value
        write_primers_manifest(manifest_path, rows)
        return self._send_json({"ok": True, "fluorophore": value})

    def _handle_set_citation(self, set_id):
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        if set_id not in load_manifest():
            return self._send_json({"error": "not found"}, 404)
        text = str(payload.get("text", "")).strip()[:600]
        url = str(payload.get("url", "")).strip()[:300]
        state = load_ui_state()
        if text or url:
            state["citations"][set_id] = {"text": text, "url": url}
        else:
            state["citations"].pop(set_id, None)
        save_ui_state(state)
        return self._send_json({"ok": True, "citation": state["citations"].get(set_id)})

    def _handle_rename_group(self, group_id):
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        name = str(payload.get("name", "")).strip()[:120]
        if not name:
            return self._send_json({"error": "이름을 입력하세요"}, 400)
        if not rename_group(group_id, name):
            return self._send_json({"error": "not found"}, 404)
        return self._send_json({"ok": True, "label": name})

    def _handle_create_verification(self, set_id):
        if not _entity_exists(set_id):
            return self._send_json({"error": "not found"}, 404)
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        material = str(payload.get("standard_material", "")).strip()
        if not material:
            return self._send_json({"error": "표준물질을 선택하세요"}, 400)
        detected_raw = payload.get("detected")
        detected = None if detected_raw in (None, "") else bool(detected_raw)
        entry = {
            "id": uuid.uuid4().hex[:16],
            "created_at": time.strftime("%Y-%m-%d %H:%M"),
            "standard_material": material,
            "detected": detected,
            "ct_value": str(payload.get("ct_value", ""))[:40],
            "fluorophore": str(payload.get("fluorophore", "")).strip()[:40],
            "note": str(payload.get("note", ""))[:2000],
            "images": [],
        }
        state = load_ui_state()
        state["verifications"].setdefault(set_id, [])
        state["verifications"][set_id].append(entry)
        save_ui_state(state)
        return self._send_json({"ok": True, "verification": entry})

    def _handle_delete_verification(self, set_id, vid):
        state = load_ui_state()
        entries = state.get("verifications", {}).get(set_id, [])
        entry = next((e for e in entries if e["id"] == vid), None)
        if entry is None:
            return self._send_json({"error": "not found"}, 404)
        img_dir = ROOT / "config" / "verification_images" / set_id / vid
        if img_dir.exists():
            shutil.rmtree(img_dir)
        state["verifications"][set_id] = [e for e in entries if e["id"] != vid]
        save_ui_state(state)
        return self._send_json({"ok": True})

    def _handle_upload_verification_image(self, set_id, vid):
        state = load_ui_state()
        entries = state.get("verifications", {}).get(set_id, [])
        entry = next((e for e in entries if e["id"] == vid), None)
        if entry is None:
            return self._send_json({"error": "not found"}, 404)
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        content_type = str(payload.get("content_type", "")).lower()
        ext = self.NOTE_IMAGE_TYPES.get(content_type)
        if not ext:
            return self._send_json({"error": "지원하지 않는 이미지 형식입니다 (png/jpg/gif/webp만 가능)"}, 400)
        data_b64 = payload.get("data_base64", "")
        try:
            raw = base64.b64decode(data_b64, validate=True)
        except (ValueError, TypeError):
            return self._send_json({"error": "이미지 데이터를 읽을 수 없습니다"}, 400)
        if not raw or len(raw) > self.NOTE_IMAGE_MAX_BYTES:
            return self._send_json({"error": "이미지 크기가 너무 크거나 비어 있습니다 (8MB 이하만 가능)"}, 400)

        image_id = uuid.uuid4().hex[:16]
        img_dir = ROOT / "config" / "verification_images" / set_id / vid
        img_dir.mkdir(parents=True, exist_ok=True)
        (img_dir / f"{image_id}.{ext}").write_bytes(raw)

        image_entry = {
            "id": image_id, "ext": ext,
            "caption": str(payload.get("caption", ""))[:200],
            "uploaded": time.strftime("%Y-%m-%d %H:%M"),
        }
        entry.setdefault("images", []).append(image_entry)
        save_ui_state(state)
        return self._send_json({"ok": True, "image": image_entry})

    def _handle_get_verification_image(self, set_id, vid, image_id):
        state = load_ui_state()
        entries = state.get("verifications", {}).get(set_id, [])
        entry = next((e for e in entries if e["id"] == vid), None)
        img = next((i for i in (entry or {}).get("images", [])if i["id"] == image_id), None)
        if img is None:
            return self._send_json({"error": "not found"}, 404)
        img_path = ROOT / "config" / "verification_images" / set_id / vid / f"{image_id}.{img['ext']}"
        content_type = {v: k for k, v in self.NOTE_IMAGE_TYPES.items()}.get(img["ext"], "application/octet-stream")
        return self._send_file(img_path, content_type)

    def _handle_delete_verification_image(self, set_id, vid, image_id):
        state = load_ui_state()
        entries = state.get("verifications", {}).get(set_id, [])
        entry = next((e for e in entries if e["id"] == vid), None)
        if entry is None:
            return self._send_json({"error": "not found"}, 404)
        img = next((i for i in entry.get("images", []) if i["id"] == image_id), None)
        if img is None:
            return self._send_json({"error": "not found"}, 404)
        img_path = ROOT / "config" / "verification_images" / set_id / vid / f"{image_id}.{img['ext']}"
        img_path.unlink(missing_ok=True)
        entry["images"] = [i for i in entry["images"] if i["id"] != image_id]
        save_ui_state(state)
        return self._send_json({"ok": True})

    def _handle_save_order(self):
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        kind = payload.get("kind")
        order = payload.get("order")
        if kind not in ("set", "group") or not isinstance(order, list):
            return self._send_json({"error": "kind('set'|'group')와 order(list)가 필요합니다."}, 400)
        state = load_ui_state()
        state["set_order" if kind == "set" else "group_order"] = [str(x) for x in order]
        save_ui_state(state)
        return self._send_json({"ok": True})

    def do_DELETE(self):  # noqa: N802
        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/note/images/([A-Za-z0-9]+)$", self.path)
        if m:
            return self._handle_delete_note_image(m.group(1), m.group(2))
        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/verify/([A-Za-z0-9]+)/images/([A-Za-z0-9]+)$", self.path)
        if m:
            return self._handle_delete_verification_image(m.group(1), m.group(2), m.group(3))
        m = re.match(r"^/api/(?:sets|groups)/([^/]+)/verify/([A-Za-z0-9]+)$", self.path)
        if m:
            return self._handle_delete_verification(m.group(1), m.group(2))
        m = re.match(r"^/api/groups/([^/]+)$", self.path)
        if m:
            ok = remove_group(m.group(1))
            if not ok:
                return self._send_json({"error": "not found"}, 404)
            return self._send_json({"ok": True})
        return self.send_error(404)

    def _handle_save_standards(self):
        """Saves the 표준물질 spreadsheet (config/standards.tsv) and regenerates the FASTA/metadata files the
        pipeline reads. Refuses to leave a segment that has active primer sets without any active standard."""
        from standards_common import (  # noqa: PLC0415
            STANDARD_COLUMNS, sync_standards, validate_rows, write_manifest,
        )
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)
        rows = payload.get("rows")
        if not isinstance(rows, list):
            return self._send_json({"error": "rows는 배열이어야 합니다"}, 400)
        rows = [{c: str(r.get(c, "")) for c in STANDARD_COLUMNS} for r in rows]
        for r in rows:
            r["name"] = r["name"].strip()
            r["segment"] = r["segment"].strip().upper()
        # standard material is optional: a segment may have none (its sets then just have no standard results)
        errors = validate_rows(rows, segments=tuple(k.upper() for k in seg_keys()))
        if errors:
            return self._send_json({"errors": errors}, 400)
        rows = [r for r in rows if r["name"] or r["sequence"].strip()]      # drop the spare blank row
        saved = write_manifest(ROOT / "config" / "standards.tsv", rows)
        changed = sync_standards(ROOT / "config" / "standards.tsv", ROOT / "config")
        for r in saved:
            r["length"] = len(r["sequence"])
        return self._send_json({"ok": True, "rows": saved, "files_changed": changed})

    def _handle_create_group(self):
        length = int(self.headers.get("Content-Length", 0))
        try:
            payload = json.loads(self.rfile.read(length) or b"{}")
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)

        label = (payload.get("label") or "").strip()
        member_ids = payload.get("member_set_ids") or []
        if not label or len(member_ids) < 2:
            return self._send_json({"error": "label과 2개 이상의 member_set_ids가 필요합니다."}, 400)

        manifest = load_manifest()
        segs = {manifest[m]["segment"] for m in member_ids if m in manifest}
        missing = [m for m in member_ids if m not in manifest]
        if missing:
            return self._send_json({"error": f"존재하지 않는 세트: {', '.join(missing)}"}, 400)
        if len(segs) != 1:
            return self._send_json({"error": "묶으려는 세트는 모두 같은 segment(S 또는 L)여야 합니다."}, 400)
        segment = next(iter(segs))

        group_id = append_group(label, segment, member_ids)

        with JOBS_LOCK:
            already_running = any(
                j.get("kind") == "pipeline" and j.get("status") == "running" for j in JOBS.values()
            )
        if already_running:
            return self._send_json({"group_id": group_id, "warning": "그룹은 만들어졌지만 다른 파이프라인 작업이 실행 중이라 결과 계산은 그 작업이 끝난 뒤 '분석 다시하기'로 실행해야 합니다."})

        job_id = uuid.uuid4().hex[:12]
        with JOBS_LOCK:
            JOBS[job_id] = {"status": "running", "kind": "pipeline", "label": f"그룹 '{label}' 계산",
                             "progress": "시작...", "started": time.time()}
        thread = threading.Thread(
            target=run_pipeline_job, args=(job_id, [
                "output/summary/group_summary.tsv", "output/dashboard/index.html",
                f"output/groups/{group_id}/combo_mismatch_table.tsv",
                f"output/groups/{group_id}/standard_combo_table.tsv",
                f"output/groups/{group_id}/refgenome_combo_table.tsv",
            ], f"그룹 '{label}' 계산"),
            daemon=True)
        thread.start()
        return self._send_json({"group_id": group_id, "job_id": job_id})

    def _handle_analyze(self):
        length = int(self.headers.get("Content-Length", 0))
        try:
            payload = json.loads(self.rfile.read(length) or b"{}")
        except json.JSONDecodeError:
            return self._send_json({"error": "invalid JSON"}, 400)

        required = ["name", "segment", "f_seq", "r_seq", "p_seq"]
        missing = [k for k in required if not payload.get(k)]
        if missing:
            return self._send_json({"error": f"missing fields: {', '.join(missing)}"}, 400)

        job_id = uuid.uuid4().hex[:12]
        with JOBS_LOCK:
            JOBS[job_id] = {"status": "running", "kind": "analyze", "progress": "시작...",
                             "label": f"새 세트 분석: {payload.get('name', '')}", "started": time.time()}
        thread = threading.Thread(target=run_analysis_job, args=(job_id, payload), daemon=True)
        thread.start()
        return self._send_json({"job_id": job_id})

    def _handle_pipeline(self, label, args, on_success=None):
        with JOBS_LOCK:
            for jid, j in JOBS.items():
                if j.get("kind") != "pipeline" or j.get("status") != "running":
                    continue
                pid = j.get("pid")
                # A "running" job whose process is actually gone (server restarted mid-job
                # without the process surviving, or it was killed outside our own exception
                # handling) would otherwise block every future run forever -- self-heal it.
                if pid is not None and not _pid_alive(pid):
                    j["status"] = "error"
                    j["error"] = "프로세스가 더 이상 존재하지 않습니다 (서버 재시작 등으로 중단된 것으로 보입니다)."
                    continue
                return self._send_json({"error": "다른 파이프라인 작업이 이미 실행 중입니다. 완료 후 다시 시도하세요."}, 409)
        job_id = uuid.uuid4().hex[:12]
        with JOBS_LOCK:
            JOBS[job_id] = {"status": "running", "kind": "pipeline", "label": label,
                             "progress": "시작...", "started": time.time()}
        thread = threading.Thread(target=run_pipeline_job, args=(job_id, args, label),
                                   kwargs={"on_success": on_success}, daemon=True)
        thread.start()
        return self._send_json({"job_id": job_id})

    def _handle_redesign(self):
        length = int(self.headers.get("Content-Length", 0))
        try:
            o = json.loads(self.rfile.read(length) or b"{}")
        except ValueError:
            o = {}
        seg = str(o.get("segment", "")).lower()
        if seg not in seg_keys():
            return self._send_json({"error": "segment가 올바르지 않습니다."}, 400)
        with JOBS_LOCK:
            if any(j.get("kind") == "redesign" and j.get("status") == "running" for j in JOBS.values()):
                return self._send_json({"error": "다른 재설계가 이미 실행 중입니다."}, 409)

        def num(key, default, lo, hi, cast=float):
            try:
                return min(hi, max(lo, cast(o.get(key, default))))
            except (TypeError, ValueError):
                return default
        thr = num("threshold", 0.9, 0.3, 1.0)
        opts = {"threshold": thr, "min_threshold": min(thr, num("min_threshold", 0.6, 0.3, 1.0)),
                "ambig": num("ambig", 2, 0, 8, int), "probe_ambig": num("probe_ambig", 1, 0, 6, int),
                "top": num("top", 3, 1, 8, int), "min_n": num("min_n", 5, 3, 1000, int),
                "include_all": bool(o.get("include_all")), "engine": "primerhunter" if o.get("engine") == "primerhunter" else "varvamp",
                "ph_sample": num("ph_sample", 25, 5, 100, int), "ph_threads": num("ph_threads", 12, 1, 32, int),
                "ph_min_cov": num("ph_min_cov", 80.0, 10.0, 100.0), "ph_min_tm": num("ph_min_tm", 56.0, 40.0, 70.0),
                "lineages": [x for x in (o.get("lineages") or []) if isinstance(x, str) and re.fullmatch(r"[A-Za-z0-9_\-]+", x)]}
        scope, ids = redesign_scope(seg, o.get("host"))
        job_id = uuid.uuid4().hex[:12]
        with JOBS_LOCK:
            JOBS[job_id] = {"status": "running", "kind": "redesign", "label": f"{seg.upper()} 리니지별 프라이머·프로브 재설계 ({'PrimerHunter' if opts['engine'] == 'primerhunter' else 'varVAMP'})",
                            "progress": "시작...", "started": time.time()}
        threading.Thread(target=run_redesign_job, args=(job_id, seg, scope, ids, opts), daemon=True).start()
        return self._send_json({"job_id": job_id})

    def _handle_probe_best(self, set_id):
        length = int(self.headers.get("Content-Length", 0))
        try:
            o = json.loads(self.rfile.read(length) or b"{}")
        except ValueError:
            o = {}
        if not _PD_ID_RE.fullmatch(set_id) or not (ROOT / "input_primers" / f"{set_id}.tsv").exists():
            return self._send_json({"error": "세트를 찾을 수 없습니다."}, 404)
        with JOBS_LOCK:
            if any(j.get("kind") == "probe_best" and j.get("status") == "running" and j.get("set_id") == set_id for j in JOBS.values()):
                return self._send_json({"error": "이 세트의 최적 Probe 탐색이 이미 실행 중입니다."}, 409)

        def num(key, default, lo, hi, cast=float):
            try:
                return min(hi, max(lo, cast(o.get(key, default))))
            except (TypeError, ValueError):
                return default
        lo = num("tm_above_min", 3.0, -5.0, 30.0)
        opts = {"top": num("top", 20, 1, 30, int), "max_iupac": num("max_iupac", 2, 0, 4, int),
                "probe_mismatch": num("probe_mismatch", 0, 0, 3, int),
                "method": o.get("method") if o.get("method") in ("own", "primer3", "both") else "both",
                "tm_above": (lo, max(lo, num("tm_above_max", 10.0, -5.0, 30.0)))}
        job_id = uuid.uuid4().hex[:12]
        with JOBS_LOCK:
            JOBS[job_id] = {"status": "running", "kind": "probe_best", "set_id": set_id,
                            "label": f"{set_id} 최적 Probe 탐색 (후보 {opts['top']}개 전체 검증)",
                            "progress": "시작...", "started": time.time()}
        scope, ids = pd_scope(set_id, o.get("host"))
        threading.Thread(target=run_probe_best_job, args=(job_id, set_id, opts, scope, ids), daemon=True).start()
        return self._send_json({"job_id": job_id})

    def _handle_probe_validate(self, set_id):
        length = int(self.headers.get("Content-Length", 0))
        try:
            body = json.loads(self.rfile.read(length) or b"{}")
            seq, host = str(body.get("seq", "")).strip().upper(), body.get("host")
        except ValueError:
            seq, host = "", None
        if not _PD_ID_RE.fullmatch(set_id) or not (ROOT / "input_primers" / f"{set_id}.tsv").exists():
            return self._send_json({"error": "세트를 찾을 수 없습니다."}, 404)
        if not _PD_SEQ_RE.match(seq):
            return self._send_json({"error": "Probe 서열은 10-60nt의 IUPAC 문자(ACGTRYSWKMBDHVN)여야 합니다."}, 400)
        job_id = uuid.uuid4().hex[:12]
        with JOBS_LOCK:
            JOBS[job_id] = {"status": "running", "kind": "probe_validate", "label": f"{set_id} 설계 Probe 검증 ({seq})",
                            "progress": "시작...", "started": time.time()}
        scope, ids = pd_scope(set_id, host)
        threading.Thread(target=run_probe_validate_job, args=(job_id, set_id, seq, scope, ids), daemon=True).start()
        return self._send_json({"job_id": job_id})

    def _handle_probe_design(self, set_id):
        """Design TaqMan probe candidates for this set's F/R (scripts/probe_designer.py) and test them against the population.
        Takes about a second, so it runs inside the request instead of as a background job."""
        import probe_designer  # noqa: PLC0415 -- scripts/ is on sys.path
        length = int(self.headers.get("Content-Length", 0))
        try:
            opts = json.loads(self.rfile.read(length) or b"{}")
        except ValueError:
            opts = {}
        if not re.fullmatch(r"[A-Za-z0-9_\-]+", set_id):
            return self._send_json({"error": "bad id"}, 400)
        def num(key, default, lo, hi, cast=float):
            try:
                return min(hi, max(lo, cast(opts.get(key, default))))
            except (TypeError, ValueError):
                return default
        try:
            scope, ids = pd_scope(set_id, opts.get("host"))
            result = probe_designer.run_design(
                set_id, ids=ids, scope=scope, top=num("top", 10, 1, 30, int), min_len=num("min_len", 18, 15, 40, int),
                max_len=num("max_len", 30, 15, 40, int), probe_mismatch=num("probe_mismatch", 0, 0, 3, int),
                max_iupac=num("max_iupac", 2, 0, 4, int), method=(opts.get("method") if opts.get("method") in ("own", "primer3", "both") else "both"),
                tm_above=(num("tm_above_min", 3.0, -5.0, 30.0), max(num("tm_above_min", 3.0, -5.0, 30.0), num("tm_above_max", 10.0, -5.0, 30.0))))
        except (FileNotFoundError, SystemExit) as exc:
            return self._send_json({"error": str(exc) or "분석 결과(MSA)가 없습니다. 먼저 이 세트를 분석하세요."}, 400)
        return self._send_json(result)

    def _handle_set_analyze(self, set_id, force):
        """Analyse (or force-re-analyse) ONE primer set: only its own jobs plus the pipeline-wide summary tables that
        depend on it. Snakemake skips everything that is already up to date, so a new set costs just that set."""
        from mfe_common import load_primers_manifest  # noqa: PLC0415
        sync_pipeline_inputs()
        row = next((r for r in load_primers_manifest(ROOT / "config" / "primers.tsv") if r["set_id"] == set_id), None)
        if row is None:
            return self._send_json({"error": f"프라이머 세트 {set_id}를 찾을 수 없습니다 (입력 페이지에 저장된 세트만 개별 분석할 수 있습니다)."}, 404)
        seg = row["segment"].lower()
        out = f"output/primers/{set_id}"
        has_std = (ROOT / "config" / f"kdca_standard_{seg}.fasta").exists()
        targets = [f"{out}/final/summary_row.tsv", f"{out}/final/lineage_rows.tsv", f"{out}/alignment/mismatch_table.tsv",
                   f"{out}/msa/msa.json", f"{out}/refgenome_binding/binding.tsv", f"{out}/dimer_hairpin/dimer_hairpin.json"]
        if has_std:
            targets += [f"{out}/standard_alignment/mismatch_table.tsv", f"{out}/standard_final/summary_row.tsv",
                        f"{out}/standard_final/lineage_rows.tsv"]
        targets += ["output/summary/summary.tsv", "output/summary/summary_by_lineage.tsv", "output/summary/gene_summary.tsv",
                    "output/summary/engine_concordance.tsv", "output/summary/standard_summary.tsv",
                    "output/summary/standard_matches.tsv"]
        args = list(targets)
        if force:
            m = re.search(r"^primer_mismatch_levels:\s*\[([^\]]*)\]", _config_yaml_path().read_text(encoding="utf-8"), re.M)
            top = max([int(x) for x in re.findall(r"\d+", m.group(1))] or [2]) if m else 2
            forced = [f"{out}/primer_qc/{set_id}_mm{top}.txt"] + ([f"{out}/standard_qc/{set_id}_mm{top}.txt"] if has_std else [])
            args = ["--forcerun", *forced, "--", *targets]
        label = f"세트 {set_id} " + ("재분석 (강제)" if force else "분석")
        return self._handle_pipeline(label=label, args=args)

    def _handle_group_analyze(self, group_id, force):
        """Analyse (or force-re-analyse) ONE primer group: its union, cross-binding, standard and reference-genome
        combination tables and the group summary. Member sets are brought up to date first if they are outdated."""
        sync_pipeline_inputs()
        group = next((g for g in load_group_list() if g["group_id"] == group_id), None)
        if group is None:
            return self._send_json({"error": f"그룹 {group_id}를 찾을 수 없습니다."}, 404)
        seg = str(group.get("segment", "")).lower()
        out = f"output/groups/{group_id}"
        has_std = (ROOT / "config" / f"kdca_standard_{seg}.fasta").exists()
        combos = [f"{out}/combo_mismatch_table.tsv", f"{out}/refgenome_combo_table.tsv"] + ([f"{out}/standard_combo_table.tsv"] if has_std else [])
        targets = [f"{out}/validated_ids.txt", f"{out}/summary_row.tsv", *combos, "output/summary/group_summary.tsv"]
        args = list(targets)
        if force:
            args = ["--forcerun", f"{out}/validated_ids.txt", *combos, "--", *targets]
        return self._handle_pipeline(label=f"그룹 {group_id} " + ("재분석 (강제)" if force else "분석"), args=args)

    def _handle_generate_report(self):
        with JOBS_LOCK:
            for jid, j in JOBS.items():
                if j.get("kind") != "report" or j.get("status") != "running":
                    continue
                pid = j.get("pid")
                if pid is not None and not _pid_alive(pid):
                    j["status"] = "error"
                    j["error"] = "프로세스가 더 이상 존재하지 않습니다."
                    continue
                return self._send_json({"error": "이미 리포트를 생성하는 중입니다. 완료 후 다시 시도하세요."}, 409)
        job_id = uuid.uuid4().hex[:12]
        # Build the base URL from the request's own Host header, so this hits the exact
        # same host:port the browser is already using (no hardcoded port to drift out of
        # sync with --port on the command line).
        base_url = f"http://{self.headers.get('Host', 'localhost:8765')}"
        with JOBS_LOCK:
            JOBS[job_id] = {"status": "running", "kind": "report", "label": "종합 리포트 PDF 생성",
                             "progress": "시작...", "started": time.time()}
        thread = threading.Thread(target=run_report_job, args=(job_id, base_url), daemon=True)
        thread.start()
        return self._send_json({"job_id": job_id})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8765)
    # 127.0.0.1 by default: the dashboard holds confidential sequences and has write endpoints but no login, so it must not be
    # reachable from the network unless asked for (--host 0.0.0.0, e.g. to open it from another PC on a trusted LAN).
    ap.add_argument("--host", default="127.0.0.1", help="address to listen on (default 127.0.0.1 = this computer only)")
    args = ap.parse_args()
    (ROOT / "output" / "user_sets").mkdir(parents=True, exist_ok=True)
    httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"[webapp] serving on http://{args.host}:{args.port}  (root={ROOT})")
    httpd.serve_forever()


if __name__ == "__main__":
    main()
