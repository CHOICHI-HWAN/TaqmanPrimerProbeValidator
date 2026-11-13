#!/usr/bin/env python3
"""Backs up and restores ALL data of one edition (this project folder) as a single .tar.xz plus a small .json manifest.

    python3 scripts/backup_data.py backup  [--out DIR] [--no-data] [--no-output] [--threads N] [--level 6] [--fast]
    python3 scripts/backup_data.py list    [--out DIR]
    python3 scripts/backup_data.py verify  FILE.tar.xz|FILE.tar.zst
    python3 scripts/backup_data.py restore FILE.tar.xz|FILE.tar.zst [--dry-run] [--port PORT]

What is in the archive (relative to the project folder): config/ (primer sets, standards, groups, notes, ids, settings), input_primers/, profiles/, VERSION,
data/ (downloaded genomes, references, filtered sets), output/ (all results) and .snakemake/metadata (so a restored folder is not re-run needlessly).
Code is NOT included: it lives in git. `--no-data` leaves out data/ (it can be downloaded again), `--no-output` leaves out output/.

The backup goes to --out, else $PRIMER_BACKUP_DIR, else ~/primer_validation_backups: NEVER inside a git repository, because it holds primer sequences
and vendor names (confidential). xz is run with all cores (-T0): smallest file, ~10 min for 80 GB. --fast uses zstd -10 instead (.tar.zst):
about 10x faster but the file is ~1.7x bigger (measured on EBOV: xz 155 MB in 11 s, zstd -10 258 MB in 1.2 s). Restore reads both by the file extension.

Restore never deletes anything: it extracts to a temporary folder next to the project, checks the sha256 from the manifest, then renames every
existing top-level folder it replaces to <project>/.pre_restore_<time>/ (an instant rename on the same disk). Delete that folder yourself when you are happy.
Restoring while a dashboard job is running is refused when --port is given and /api/jobs/active is not empty.
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOP = ["config", "input_primers", "profiles", "VERSION", "data", "output", ".snakemake/metadata"]
SKIP_SUFFIX = (".tmp",)


def default_out():
    return Path(os.environ.get("PRIMER_BACKUP_DIR") or Path.home() / "primer_validation_backups")


def edition_name(root=ROOT):
    return root.name.replace("pipeline_snakemake_v0.19.0", "lasv").replace("pipeline_", "")


def human(n):
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or u == "TB":
            return f"{n:.1f} {u}" if u != "B" else f"{n} B"
        n /= 1024


def tree_size(p):
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file() and not f.is_symlink())


def sha256(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def inside_git(path):
    p = path.resolve()
    return any((q / ".git").exists() for q in [p, *p.parents])


def backup(root=ROOT, out=None, data=True, output=True, threads=0, level=6, fast=False, log=print):
    """Creates <out>/<edition>_<time>.tar.xz (or .tar.zst with fast=True) (+ .json). Returns the manifest dict."""
    out = Path(out) if out else default_out()
    if inside_git(out):
        raise SystemExit(f"[error] {out} is inside a git repository; a backup holds primer sequences and must stay outside git")
    out.mkdir(parents=True, exist_ok=True)
    names = [n for n in TOP if (root / n).exists() and (data or n != "data") and (output or n != "output")]
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    ext = "tar.zst" if fast else "tar.xz"
    target = out / f"{edition_name(root)}_{stamp}{'-fast' if fast else ''}.{ext}"          # one manifest (<name>.json) per archive
    sizes = {n: tree_size(root / n) for n in names}
    log(f"backup {edition_name(root)}: {', '.join(f'{n} {human(s)}' for n, s in sizes.items())}")
    started = time.time()
    tar = subprocess.Popen(["tar", "-cf", "-", "-C", str(root), *[f"--exclude=*{s}" for s in SKIP_SUFFIX], *names], stdout=subprocess.PIPE)
    tmp = target.with_name(target.name + ".part")
    with open(tmp, "wb") as fh:
        comp = ["zstd", f"-T{threads}", "-10", "-c"] if fast else ["xz", f"-T{threads}", f"-{level}", "-c"]
        xz = subprocess.Popen(comp, stdin=tar.stdout, stdout=fh)
        tar.stdout.close()
        rc_xz = xz.wait()
    rc_tar = tar.wait()
    if rc_tar not in (0, 1) or rc_xz != 0:              # tar exits 1 when a file changed while it was read (a running dashboard): the archive is still usable
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"[error] backup failed (tar {rc_tar}, xz {rc_xz})")
    tmp.rename(target)
    manifest = {"edition": edition_name(root), "created": datetime.datetime.now().isoformat(timespec="seconds"),
                "version": (root / "VERSION").read_text().strip() if (root / "VERSION").exists() else "",
                "includes": names, "compression": "zstd -10" if fast else f"xz -{level}", "archive": target.name, "uncompressed_bytes": sum(sizes.values()), "archive_bytes": target.stat().st_size,
                "sha256": sha256(target), "seconds": round(time.time() - started, 1), "tar_note": "files changed during the backup" if rc_tar == 1 else ""}
    target.with_suffix("").with_suffix(".json").write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
    log(f"-> {target} ({human(manifest['archive_bytes'])} from {human(manifest['uncompressed_bytes'])}, {manifest['seconds']} s)")
    return {**manifest, "file": str(target)}


def manifest_of(archive):
    m = Path(archive).with_suffix("").with_suffix(".json")
    return json.loads(m.read_text(encoding="utf-8")) if m.exists() else None


def verify(archive, log=print):
    archive = Path(archive)
    m = manifest_of(archive)
    ok = True
    if m and m.get("sha256"):
        good = sha256(archive) == m["sha256"]
        log(f"sha256 {'ok' if good else 'MISMATCH'}")
        ok = good
    else:
        log("no manifest: sha256 not checked")
    test = subprocess.run(["zstd", "-tq", str(archive)] if archive.name.endswith(".zst") else ["xz", "-t", str(archive)], capture_output=True, text=True)
    kind = "zstd" if archive.name.endswith(".zst") else "xz"
    log(f"{kind} integrity ok" if test.returncode == 0 else f"{kind} integrity FAILED: {test.stderr.strip()}")
    return ok and test.returncode == 0


def safe_members(archive):
    """Names in the archive; refuses absolute paths, '..' and links that point outside."""
    names = []
    if str(archive).endswith(".zst"):                         # python's tarfile cannot read zstd (3.12): list it with tar
        listing = subprocess.run(["tar", "--zstd", "-tvf", str(archive)], capture_output=True, text=True, check=True).stdout.splitlines()
        for line in listing:
            name, link = line.split(None, 5)[5], ""
            for sep in (" -> ", " link to "):
                if sep in name:
                    name, link = name.split(sep, 1)
            p = Path(name)
            if p.is_absolute() or ".." in p.parts or (link and (Path(link).is_absolute() or ".." in Path(link).parts)):
                raise SystemExit(f"[error] unsafe path in archive: {line}")
            names.append(name)
        return names
    with tarfile.open(archive, "r:xz") as tf:
        for m in tf:
            p = Path(m.name)
            if p.is_absolute() or ".." in p.parts or ((m.issym() or m.islnk()) and (Path(m.linkname).is_absolute() or ".." in Path(m.linkname).parts)):
                raise SystemExit(f"[error] unsafe path in archive: {m.name}")
            names.append(m.name)
    return names


def restore(archive, root=ROOT, dry_run=False, port=None, log=print):
    archive = Path(archive)
    if port:
        try:
            busy = urllib.request.urlopen(f"http://127.0.0.1:{port}/api/jobs/active", timeout=3).read().decode().strip()
        except Exception:
            busy = "[]"                                      # dashboard not running: nothing to disturb
        if busy not in ("[]", ""):
            raise SystemExit(f"[error] a dashboard job is running on :{port}; restore after it finished ({busy[:80]})")
    if not verify(archive, log):
        raise SystemExit("[error] the archive is damaged or does not match its manifest; nothing was changed")
    names = safe_members(archive)
    tops = sorted({n.split("/")[0] if not n.startswith(".snakemake/") else ".snakemake/metadata" for n in names})
    m = manifest_of(archive)
    log(f"restore {archive.name}: replaces {', '.join(tops)} in {root}" + (f" (edition {m['edition']}, created {m['created']})" if m else ""))
    if m and m["edition"] != edition_name(root):
        raise SystemExit(f"[error] this archive is edition {m['edition']}, the folder is {edition_name(root)}")
    if m:
        need = m["uncompressed_bytes"]
        free = shutil.disk_usage(root).free
        if free < need * 1.1:
            raise SystemExit(f"[error] not enough free space: need {human(need)}, free {human(free)}")
    if dry_run:
        log("dry run: nothing changed")
        return {"replaced": tops}
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S-%f")          # microseconds: two restores in one second must not share a folder
    tmp, old = root / f".restore_tmp_{stamp}", root / f".pre_restore_{stamp}"
    tmp.mkdir()
    try:
        subprocess.run(["tar", "--zstd" if archive.name.endswith(".zst") else "-J", "-xf", str(archive), "-C", str(tmp)], check=True)
        old.mkdir()
        for t in tops:
            src, dst = tmp / t, root / t
            if not src.exists():
                continue
            if dst.exists():
                (old / t).parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(dst), str(old / t))
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    log(f"restored. The previous versions are in {old} (delete it when you do not need it)")
    return {"replaced": tops, "previous": str(old)}


def listing(out=None):
    out = Path(out) if out else default_out()
    rows = []
    for j in sorted(list(out.glob(f"{edition_name()}_*.json")) + list(out.glob(f"*/{edition_name()}_*.json"))):
        m = json.loads(j.read_text(encoding="utf-8"))
        m["file"] = str(j.with_name(m.get("archive") or (j.stem + ".tar.xz")))
        m["exists"] = Path(m["file"]).exists()
        rows.append(m)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("backup")
    b.add_argument("--out")
    b.add_argument("--no-data", action="store_true")
    b.add_argument("--no-output", action="store_true")
    b.add_argument("--threads", type=int, default=0, help="xz threads (0 = all cores)")
    b.add_argument("--level", type=int, default=6)
    b.add_argument("--fast", action="store_true", help="zstd -10: ~10x faster, ~1.7x bigger than xz")
    for name in ("verify", "restore"):
        s = sub.add_parser(name)
        s.add_argument("file")
        if name == "restore":
            s.add_argument("--dry-run", action="store_true")
            s.add_argument("--port", type=int)
    sub.add_parser("list").add_argument("--out")
    a = ap.parse_args()
    if a.cmd == "backup":
        backup(out=a.out, data=not a.no_data, output=not a.no_output, threads=a.threads, level=a.level, fast=a.fast)
    elif a.cmd == "verify":
        sys.exit(0 if verify(a.file) else 1)
    elif a.cmd == "restore":
        restore(a.file, dry_run=a.dry_run, port=a.port)
    else:
        for m in listing(a.out):
            print(f"{m['created']}  {m['version']:8} {human(m['archive_bytes']):>10}  {','.join(m['includes'])}  {m['file']}")


if __name__ == "__main__":
    main()
