#!/usr/bin/env python3
"""Front door for the per-pathogen primer/probe validation pipelines (LASV, Nipah, ...).

Reads pipelines.json, shows every pipeline as a card with live status (is its web app up?) and a few
numbers read straight from that project's own result files (so they show while the app is stopped),
and can start / stop a pipeline's web app. Only apps started from here can be stopped from here; an
app someone started by hand is shown as "running (external)" and left alone. Pure standard library.

Run:
    python3 server.py [--port 8760]
"""
import argparse
import json
import os
import re
import signal
import socket
import subprocess
import sys
import threading
import time
import urllib.parse
import urllib.request
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent
STATIC_DIR = HERE / "static"
LOG_DIR = HERE / "logs"
BASE = HERE.parent                     # the folder that holds LASV/, NIPAH/, entry/
STARTED = {}                           # pipeline id -> subprocess.Popen started from this page

sys.path.insert(0, str(BASE / "tools"))
import backup_restore as br  # noqa: E402  -- ALL editions' data backup/restore (per-edition dialog: webapp/static/backup.js)

JOBS = {}                              # job_id -> {"status": running|done|error, "progress": str, "kind": "backup"|"restore", ...}
JOBS_LOCK = threading.Lock()


def _job_log(job_id):
    def log(msg):
        with JOBS_LOCK:
            JOBS[job_id]["progress"] = str(msg)[-400:]
    return log


def run_full_backup_job(job_id, data, output, fast):
    log = _job_log(job_id)
    try:
        res = br.run_backup(data=data, output=output, fast=fast, log=log)
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error" if res["failed"] else "done"
            n = len(res["editions"])
            JOBS[job_id]["progress"] = f"{n}개 에디션 완료" + (f", {len(res['failed'])}개 실패: {', '.join(x[0] for x in res['failed'])}" if res["failed"] else "")
    except (Exception, SystemExit) as exc:  # noqa: BLE001
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error"
            JOBS[job_id]["progress"] = f"백업 실패: {exc}"


def run_full_restore_job(job_id, folder, shared):
    log = _job_log(job_id)
    try:
        res = br.run_restore(folder, shared=shared, log=log)
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error" if res["skipped"] else "done"
            JOBS[job_id]["progress"] = f"{len(res['restored'])}개 복원" + (f", {len(res['skipped'])}개 건너뜀: {', '.join(x[0] for x in res['skipped'])}" if res["skipped"] else "")
    except (Exception, SystemExit) as exc:  # noqa: BLE001
        with JOBS_LOCK:
            JOBS[job_id]["status"] = "error"
            JOBS[job_id]["progress"] = f"복원 실패: {exc}"


def _any_job_running():
    with JOBS_LOCK:
        return any(j.get("status") == "running" for j in JOBS.values())


def registry():
    return json.loads((HERE / "pipelines.json").read_text(encoding="utf-8"))


def find_pipeline(pid):
    return next((p for p in registry()["pipelines"] if p["id"] == pid), None)


def port_open(port, timeout=0.4):
    try:
        with socket.create_connection(("127.0.0.1", port), timeout=timeout):
            return True
    except OSError:
        return False


def app_profile(port):
    """The running app's own /api/profile (None if it is not a validation web app)."""
    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/api/profile", timeout=1.5) as r:
            return json.loads(r.read())
    except Exception:  # noqa: BLE001
        return None


def _read(path):
    try:
        return Path(path).read_text(encoding="utf-8")
    except OSError:
        return ""


def project_stats(root: Path):
    """Numbers from the project's result files; any missing file just leaves that number out."""
    out = {"exists": root.is_dir()}
    if not out["exists"]:
        return out
    ver = _read(root / "VERSION").strip()
    if ver:
        out["version"] = ver
    cfg = _read(root / "config" / "config.yaml")
    m = re.search(r'^\s*profile:\s*"?([^"\s#]+)"?', cfg, re.M)
    prof = {}
    try:
        prof = json.loads(_read(root / (m.group(1) if m else "profiles/lassa.json")))
    except ValueError:
        pass
    out["lineage_label"] = (prof.get("lineage") or {}).get("label", "Lineage")
    out["taxonomy"] = prof.get("taxonomy")                      # ICTV ranks from the project profile (shown on the card)
    out["segments"] = [s.get("label", k) for k, s in (prof.get("segments") or {}).items()]
    summary = root / "output" / "summary" / "summary.tsv"
    rows = []
    text = _read(summary)
    if text:
        lines = [l.split("\t") for l in text.splitlines() if l.strip()]
        head = lines[0]
        rows = [dict(zip(head, r)) for r in lines[1:]]
        out["analysed_at"] = time.strftime("%Y-%m-%d %H:%M", time.localtime(summary.stat().st_mtime))
    out["sets"] = len({r.get("set_id") for r in rows})
    best = None
    for r in rows:
        try:
            v = float(r.get("validated_pct", ""))
        except ValueError:
            continue
        if best is None or v > best:
            best = v
    if best is not None:
        out["best_validated_pct"] = round(best, 1)
    kept, lineages = 0, set()
    for rep in sorted((root / "data" / "filtered").glob("*_filter_report.txt")):
        t = _read(rep)
        km = re.search(r"kept=(\d+)", t)
        kept += int(km.group(1)) if km else 0
        sec = t.split("lineage counts:")[-1].split("country counts:")[0] if "lineage counts:" in t else ""
        lineages.update(l.split("\t")[0].strip() for l in sec.splitlines() if "\t" in l)
    if kept:
        out["genomes"] = kept
    if lineages:
        out["lineages"] = len(lineages - {"unclassified", "NA", "Unknown"})
    return out


def pipeline_status(p):
    root = BASE / p["path"]
    up = port_open(p["port"])
    prof = app_profile(p["port"]) if up else None
    proc = STARTED.get(p["id"])
    managed = bool(proc and proc.poll() is None)
    if proc and proc.poll() is not None:
        STARTED.pop(p["id"], None)
    return {
        "id": p["id"], "name": p["name"], "sci": p["sci"], "desc": p["desc"], "port": p["port"], "path": p["path"],
        "state": ("running" if prof else "port_busy") if up else "stopped",
        "managed": managed, "app_virus": (prof or {}).get("short_name"),
        "stats": project_stats(root),
    }


def start_pipeline(p):
    root = BASE / p["path"]
    if not (root / "webapp" / "server.py").exists():
        return 404, {"error": f"{p['path']}/webapp/server.py 가 없습니다."}
    if port_open(p["port"]):
        return 409, {"error": f"포트 {p['port']}이 이미 사용 중입니다."}
    LOG_DIR.mkdir(exist_ok=True)
    log = open(LOG_DIR / f"{p['id']}.log", "ab")
    proc = subprocess.Popen([sys.executable, "webapp/server.py", "--port", str(p["port"])], cwd=root,
                            stdout=log, stderr=log, stdin=subprocess.DEVNULL, start_new_session=True)
    STARTED[p["id"]] = proc
    for _ in range(50):                                   # up to ~10 s for the port to open
        if proc.poll() is not None:
            STARTED.pop(p["id"], None)
            return 500, {"error": f"서버가 바로 종료되었습니다. 로그: entry/logs/{p['id']}.log"}
        if port_open(p["port"]):
            return 200, {"ok": True}
        time.sleep(0.2)
    return 504, {"error": "서버가 10초 안에 열리지 않았습니다.", "log": f"entry/logs/{p['id']}.log"}


def stop_pipeline(p):
    proc = STARTED.get(p["id"])
    if not proc or proc.poll() is not None:
        return 409, {"error": "이 페이지에서 시작한 서버가 아니라 여기서 중지할 수 없습니다."}
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()
    STARTED.pop(p["id"], None)
    return 200, {"ok": True}


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):  # noqa: A003
        line = fmt % args
        if "GET /api/pipelines " in line or "favicon.ico" in line:     # the page polls every 5 s; keep the log readable
            return
        print(f"[entry] {self.address_string()} {line}")

    def _send(self, code, body: bytes, ctype):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj, ensure_ascii=False).encode(), "application/json; charset=utf-8")

    def _read_json_body(self):
        n = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(n) if n else b"{}"
        return json.loads(raw or b"{}")

    def do_GET(self):  # noqa: N802
        path = self.path.split("?")[0]
        if path in ("/", "/index.html"):
            return self._send(200, (STATIC_DIR / "index.html").read_bytes(), "text/html; charset=utf-8")
        if path == "/api/pipelines":
            reg = registry()
            return self._json({"title": reg.get("title", ""), "pipelines": [pipeline_status(p) for p in reg["pipelines"]],
                               "planned": reg.get("planned", [])})
        if path == "/api/backup/sets":
            # Every pathogen panel already has its OWN backup/restore dialog (webapp/static/backup.js) for just that edition; this is the
            # all-editions version: one backup SET (tools/backup_restore.py) across every pipeline plus the shared files (Firebase allow
            # lists, this page's own pipeline list, the dashboards' systemd units).
            return self._json({"dir": str(br.bd.default_out()), "editions": sorted(br.editions()), "sets": br.list_sets()})
        if path == "/api/jobs/active":
            with JOBS_LOCK:
                active = [dict(j, job_id=jid) for jid, j in JOBS.items() if j.get("status") == "running"]
            return self._json(active)
        m = re.fullmatch(r"/api/jobs/([^/]+)", path)
        if m:
            with JOBS_LOCK:
                j = JOBS.get(m.group(1))
            return self._json(j, 200) if j else self._json({"error": "unknown job"}, 404)
        self.send_error(404)

    def do_POST(self):  # noqa: N802
        path = self.path.split("?")[0]
        if path == "/api/backup/create":
            try:
                body = self._read_json_body()
            except json.JSONDecodeError:
                return self._json({"error": "invalid JSON"}, 400)
            if _any_job_running():
                return self._json({"error": "다른 전체 백업·복원 작업이 실행 중입니다. 끝난 뒤 다시 시도하세요."}, 409)
            job_id = uuid.uuid4().hex[:12]
            with JOBS_LOCK:
                JOBS[job_id] = {"status": "running", "kind": "backup", "label": "전체 데이터 백업", "progress": "시작..."}
            threading.Thread(target=run_full_backup_job, args=(job_id, bool(body.get("data", True)), bool(body.get("output", True)),
                                                                bool(body.get("fast", False))), daemon=True).start()
            return self._json({"job_id": job_id})
        if path == "/api/backup/restore":
            try:
                body = self._read_json_body()
            except json.JSONDecodeError:
                return self._json({"error": "invalid JSON"}, 400)
            folder_name = str(body.get("folder", ""))
            match = next((s for s in br.list_sets() if s["folder"] == folder_name), None)
            if match is None or not re.fullmatch(r"[0-9]{8}-[0-9]{6}", folder_name):
                return self._json({"error": "백업 세트를 찾을 수 없습니다."}, 404)
            if body.get("confirm") is not True:
                return self._json({"error": "복원은 확인(confirm)이 필요합니다."}, 400)
            if _any_job_running():
                return self._json({"error": "다른 전체 백업·복원 작업이 실행 중입니다. 끝난 뒤 다시 시도하세요."}, 409)
            job_id = uuid.uuid4().hex[:12]
            with JOBS_LOCK:
                JOBS[job_id] = {"status": "running", "kind": "restore", "label": f"전체 데이터 복원 ({folder_name})", "progress": "시작..."}
            threading.Thread(target=run_full_restore_job, args=(job_id, match["path"], bool(body.get("shared", False))), daemon=True).start()
            return self._json({"job_id": job_id})
        m = re.fullmatch(r"/api/pipelines/([a-z0-9_-]+)/(start|stop)", path)
        p = find_pipeline(m.group(1)) if m else None
        if not p:
            return self._json({"error": "unknown pipeline"}, 404)
        code, body = (start_pipeline if m.group(2) == "start" else stop_pipeline)(p)
        self._json(body, code)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8760)
    # 127.0.0.1 by default: the dashboard holds confidential sequences and has write endpoints but no login, so it must not be
    # reachable from the network unless asked for (--host 0.0.0.0, e.g. to open it from another PC on a trusted LAN).
    ap.add_argument("--host", default="127.0.0.1", help="address to listen on (default 127.0.0.1 = this computer only)")
    args = ap.parse_args()
    try:
        httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    except OSError as e:
        sys.exit(f"[entry] 포트 {args.port}을 열 수 없습니다 ({e.strerror}). 이미 실행 중이면 http://localhost:{args.port} 를 여세요, "
                 f"아니면 다른 포트로: ./start.sh <port>")
    print(f"[entry] serving on http://{args.host}:{args.port}")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
