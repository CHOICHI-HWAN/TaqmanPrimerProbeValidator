"""Reference-standard material (표준물질) definitions: the spreadsheet `config/standards.tsv`.

The standards used to be two hand-placed FASTA files (config/kdca_standard_<seg>.fasta plus a
metadata TSV). They are now a manifest people edit in the web app, exactly like the primer
manifest; `sync_standards()` writes the very same two files from it, so the Snakemake rules, the
live "new primer set" analysis and the tests keep reading the paths they always read. A file is
rewritten only when its bytes would change, so saving an unchanged manifest re-runs nothing.

A standard's *name* is its record id in every result (standard_genomes, the wet-lab verification
records, the representativeness report), so it must be unique, filename-safe and is treated as a
key: renaming one orphans results/verifications that mention the old name.
"""
import csv
import datetime
import io
import re
from pathlib import Path

STANDARD_COLUMNS = ["name", "segment", "sequence", "source", "remark", "proceed", "add_date"]
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")
VALID_BASES = set("ACGTURYSWKMBDHVN")
MIN_LEN = 30
LEGACY_SOURCE = "KDCA_2024_STANDARD"


def _truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes", "y")


def read_manifest(path):
    path = Path(path)
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def _write_if_changed(path, text):
    path = Path(path)
    if path.exists() and path.read_bytes() == text.encode("utf-8"):
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(text.encode("utf-8"))
    return True


def write_manifest(path, rows):
    """Stamps add_date on rows that have a name but no date yet, then writes the TSV (only if changed)."""
    today = datetime.date.today().isoformat()
    out = [dict(r) for r in rows]
    for r in out:
        if (r.get("name") or "").strip() and not (r.get("add_date") or "").strip():
            r["add_date"] = today
        r["sequence"] = re.sub(r"\s+", "", r.get("sequence") or "").upper()
    buf = io.StringIO(newline="")
    w = csv.DictWriter(buf, fieldnames=STANDARD_COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for r in out:
        w.writerow({c: r.get(c, "") for c in STANDARD_COLUMNS})
    _write_if_changed(path, buf.getvalue())
    return out


def active_rows(rows):
    return [r for r in rows if (r.get("name") or "").strip() and (r.get("sequence") or "").strip() and _truthy(r.get("proceed", "True"))]


def validate_rows(rows, segments=("S", "L"), required_segments=()):
    """[{row, field, message}] -- empty list = OK. Fully blank rows are ignored (the page keeps a spare
    one). `required_segments`: segments that have active primer sets and so need >= 1 active standard."""
    errors, seen = [], {}
    for i, r in enumerate(rows):
        name = (r.get("name") or "").strip()
        seq = re.sub(r"\s+", "", r.get("sequence") or "").upper()
        if not (name or seq or (r.get("source") or "").strip()):
            continue
        if not name:
            errors.append({"row": i, "field": "name", "message": "이름이 비어 있습니다"})
        elif not NAME_RE.match(name):
            errors.append({"row": i, "field": "name", "message": "이름에는 영문·숫자·_ . - 만 쓸 수 있습니다 (결과 파일의 ID로 쓰입니다)"})
        elif name.lower() in seen:
            errors.append({"row": i, "field": "name", "message": f"이름이 중복됩니다 (행 {seen[name.lower()] + 1})"})
        else:
            seen[name.lower()] = i
        if (r.get("segment") or "").strip().upper() not in segments:
            errors.append({"row": i, "field": "segment", "message": f"segment는 {' 또는 '.join(segments)} 이어야 합니다"})
        if not seq:
            errors.append({"row": i, "field": "sequence", "message": "서열이 비어 있습니다"})
        else:
            bad = sorted({c for c in seq if c not in VALID_BASES})
            if bad:
                errors.append({"row": i, "field": "sequence", "message": f"올바르지 않은 문자: {', '.join(bad)} (IUPAC 염기만 허용)"})
            elif len(seq) < MIN_LEN:
                errors.append({"row": i, "field": "sequence", "message": f"서열이 너무 짧습니다 ({len(seq)}bp, 최소 {MIN_LEN}bp)"})
        if str(r.get("proceed", "True")).strip().lower() not in ("true", "false", "1", "0", "yes", "no", "y", "n", ""):
            errors.append({"row": i, "field": "proceed", "message": "분석진행은 True 또는 False"})
    active = active_rows(rows)
    have = {(r.get("segment") or "").strip().upper() for r in active}
    for seg in required_segments:
        if seg.upper() not in have:
            errors.append({"row": -1, "field": "segment", "message": f"segment {seg.upper()}에는 분석 대상 프라이머 세트가 있어서 표준물질이 최소 1개(분석진행=True) 필요합니다"})
    return errors


def sync_standards(manifest_path, config_dir):
    """Writes config/kdca_standard_<seg>.fasta and _metadata.tsv from the active manifest rows. A segment with no
    active standard has its files (and the mfeprimer index beside them) removed, so the pipeline treats standard
    material as absent for it. Returns the files that changed."""
    config_dir = Path(config_dir)
    by_seg = {}
    for r in active_rows(read_manifest(manifest_path)):
        by_seg.setdefault(r["segment"].strip().lower(), []).append(r)
    changed = []
    for seg, rows in by_seg.items():
        fasta = "".join(f">{r['name'].strip()}\n{r['sequence'].strip()}\n" for r in rows)
        buf = io.StringIO(newline="")
        w = csv.writer(buf, delimiter="\t")       # default \r\n line ends, as the original hand-made file had
        w.writerow(["accession", "lineage", "date", "country", "region"])
        for r in rows:
            src = re.sub(r"\s+", "_", (r.get("source") or "").strip()) or "STANDARD"
            w.writerow([r["name"].strip(), r["name"].strip(), src, "-", "-"])
        if _write_if_changed(config_dir / f"kdca_standard_{seg}.fasta", fasta):
            changed.append(f"kdca_standard_{seg}.fasta")
        if _write_if_changed(config_dir / f"kdca_standard_{seg}_metadata.tsv", buf.getvalue()):
            changed.append(f"kdca_standard_{seg}_metadata.tsv")
    for fasta in sorted(config_dir.glob("kdca_standard_*.fasta")):          # segments whose standards were all removed
        seg = fasta.name[len("kdca_standard_"):-len(".fasta")]
        if seg in by_seg:
            continue
        for stale in config_dir.glob(f"kdca_standard_{seg}.fasta*"):
            stale.unlink()
            changed.append(stale.name)
        meta = config_dir / f"kdca_standard_{seg}_metadata.tsv"
        if meta.exists():
            meta.unlink()
            changed.append(meta.name)
    return changed


def ensure_standards_manifest(manifest_path, config_dir, segments=("s", "l")):
    """First run after the manifest was introduced: builds it from the legacy FASTA/metadata files
    (one row per record, source taken from the old metadata's date column). No-op if it exists."""
    manifest_path, config_dir = Path(manifest_path), Path(config_dir)
    if manifest_path.exists():
        return read_manifest(manifest_path)
    rows = []
    for seg in segments:
        fa = config_dir / f"kdca_standard_{seg}.fasta"
        if not fa.exists():
            continue
        meta = {}
        mp = config_dir / f"kdca_standard_{seg}_metadata.tsv"
        if mp.exists():
            with open(mp, newline="", encoding="utf-8") as fh:
                meta = {m["accession"]: m for m in csv.DictReader(fh, delimiter="\t")}
        name, chunks = None, []
        def flush():
            if name:
                rows.append({"name": name, "segment": seg.upper(), "sequence": "".join(chunks).upper(),
                             "source": (meta.get(name) or {}).get("date") or LEGACY_SOURCE, "remark": "", "proceed": "True", "add_date": ""})
        for line in fa.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                flush(); name, chunks = line[1:].split()[0], []
            elif line.strip():
                chunks.append(line.strip())
        flush()
    return write_manifest(manifest_path, rows)
