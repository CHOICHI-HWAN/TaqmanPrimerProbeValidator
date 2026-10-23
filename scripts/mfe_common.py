"""Shared helpers for invoking the mfeprimer static binary from Snakemake script: blocks."""
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

COVERAGE_RE = re.compile(
    r"Coverage\s+([\d.]+)%,\s*total records in database:\s*(\d+),\s*amplicon number:\s*(\d+),\s*missed number:\s*(\d+)"
)


def read_config_str(pipeline_root, key: str) -> str:
    """A top-level `key: "value"` line of <project>/config/config.yaml (regex: callers such as the web server have no
    YAML library). Empty string when the key or the file is missing."""
    cfg = Path(pipeline_root) / "config" / "config.yaml"
    if cfg.exists():
        m = re.search(rf'^{re.escape(key)}:\s*"?([^"\n#]*?)"?\s*(?:#.*)?$', cfg.read_text(encoding="utf-8"), re.MULTILINE)
        if m:
            return m.group(1).strip()
    return ""


def env_tool(pipeline_root, name: str) -> str:
    """<env_prefix>/bin/<name> for the conda env that setup.sh recorded in config.yaml (`env_prefix`), else ''."""
    prefix = read_config_str(pipeline_root, "env_prefix")
    if prefix:
        p = Path(prefix) / "bin" / name
        if p.exists():
            return str(p)
    return ""


def find_mafft() -> str:
    configured = env_tool(Path(__file__).resolve().parent.parent, "mafft")   # the env setup.sh built wins over PATH
    if configured:
        return configured
    found = shutil.which("mafft")
    if found:
        return found
    dedicated = Path("/BI/conda/envs/lasv_pipeline/bin/mafft")
    if dedicated.exists():
        return str(dedicated)
    for cand in sorted(Path("/BI/conda/envs").glob("*/bin/mafft")):
        return str(cand)
    sys.exit("[error] mafft binary not found. Run ./setup.sh, or put mafft on PATH.")


def _configured_mfeprimer(pipeline_root: Path) -> str:
    """`mfeprimer_bin:` from <project>/config/config.yaml, so the web server and live analysis find the same
    binary the Snakemake jobs get from config."""
    return read_config_str(pipeline_root, "mfeprimer_bin")


def check_mfeprimer(path: str, timeout: float = 8.0) -> dict:
    """Does `path` point at a runnable mfeprimer? -> {ok, path, version, error}. Never raises."""
    out = {"ok": False, "path": path, "version": "", "error": ""}
    p = Path(path).expanduser() if path else None
    if not p or not p.is_file():
        out["error"] = "파일을 찾을 수 없습니다" if path else "경로가 비어 있습니다"
        return out
    if not os.access(p, os.X_OK):
        out["error"] = "실행 권한이 없습니다 (chmod +x 필요)"
        return out
    try:
        r = subprocess.run([str(p), "version"], capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        out["error"] = f"실행할 수 없습니다: {exc}"
        return out
    text = (r.stdout + r.stderr).strip()
    if r.returncode != 0 or "mfeprimer" not in text.lower():
        out["error"] = "mfeprimer 실행 파일이 아닌 것 같습니다 (`version` 출력이 이상함)"
        return out
    out.update(ok=True, version=text.splitlines()[0].strip())
    return out


def find_mfeprimer(configured: str, pipeline_root: Path) -> str:
    configured = configured or _configured_mfeprimer(pipeline_root)
    if configured:
        p = Path(configured)
        if p.exists():
            return str(p)
    env = os.environ.get("MFEPRIMER_BIN")
    if env and Path(env).exists():
        return env
    candidates = sorted([c for pat in ("linux", "darwin") for c in pipeline_root.parent.glob(f"*/mfeprimer/mfeprimer-*/mfeprimer-*-{pat}-*")])
    if candidates:
        return str(candidates[0])
    import shutil
    found = shutil.which("mfeprimer")
    if found:
        return found
    sys.exit("[error] mfeprimer binary not found. Run ./setup.sh, or set mfeprimer_bin on the web app's parameter page / $MFEPRIMER_BIN.")


def ensure_index(mfeprimer_bin: str, db_fasta: Path, cpu: int = None):
    fai = Path(str(db_fasta) + ".fai")
    if fai.exists():
        return
    cmd = [mfeprimer_bin, "index", "-i", str(db_fasta)]
    if cpu:
        cmd += ["--cpu", str(cpu)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] mfeprimer index failed for {db_fasta}:\n{result.stderr}")


def write_oligo_fasta(path: Path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for oid, seq in records:
            fh.write(f">{oid}\n{seq}\n")


def run_primer_search(mfeprimer_bin: str, db_fasta: Path, primer_fa: Path, out_txt: Path, mismatch: int,
                       cpu: int = None):
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    # mfeprimer refuses to overwrite an existing output file -- clean up any sibling
    # files from a previous run (.spec.tsv, .virus.Failed.txt, etc.) even if the main
    # report itself is absent, since a prior run can leave siblings without the report.
    for p in out_txt.parent.glob(out_txt.name + "*"):
        p.unlink()
    cmd = [
        mfeprimer_bin, "-d", str(db_fasta), "-i", str(primer_fa), "-o", str(out_txt),
        "--fasta", "--virus", "--misMatch", str(mismatch),
    ]
    # mfeprimer's own --virus/spec search defaults to just 4 CPUs regardless of what the
    # machine has -- this is the single most expensive step in the whole pipeline (a
    # degenerate F/R pair search can otherwise take many minutes even on a large machine),
    # so an explicit --cpu is worth passing whenever the caller has one configured.
    if cpu:
        cmd += ["--cpu", str(cpu)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] mfeprimer primer search failed ({primer_fa}, mm={mismatch}):\n{result.stderr}")
    return out_txt


def parse_primer_coverage(out_txt: Path):
    text = out_txt.read_text(encoding="utf-8", errors="replace")
    m = COVERAGE_RE.search(text)
    failed_file = Path(str(out_txt) + ".virus.Failed.txt")
    failed_ids = set()
    if failed_file.exists():
        for line in failed_file.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                failed_ids.add(line.split()[0])
    if not m:
        return {"coverage_pct": 0.0, "total": 0, "missed": 0, "failed_ids": failed_ids}
    coverage_pct, total, _amp, missed = m.groups()
    return {"coverage_pct": float(coverage_pct), "total": int(total), "missed": int(missed), "failed_ids": failed_ids}


def run_probe_search(mfeprimer_bin: str, db_fasta: Path, probe_fa: Path, out_txt: Path, mismatch: int = 0):
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    if out_txt.exists():
        out_txt.unlink()
    cmd = [mfeprimer_bin, "search", str(db_fasta), str(probe_fa), str(out_txt), "-m", str(mismatch)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] mfeprimer probe search failed ({probe_fa}):\n{result.stderr}")
    matched_ids = set()
    with open(out_txt, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            matched_ids.add(line.split("\t")[0])
    return matched_ids


def read_fasta_ids(path: Path):
    ids = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def read_fasta_seqs(path: Path) -> dict:
    """id -> full sequence. Used to slice amplicon sequences directly from the genome by
    coordinate (spec.tsv's ampStart/ampEnd are always complete, unlike mfeprimer's own
    --fasta amplicon output, which silently truncates after --max-amp-count candidates --
    see build_alignment_table.py)."""
    seqs = {}
    header = None
    chunks = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            seqs[header] = "".join(chunks)
    return seqs


# ---------------------------------------------------------------------------------------
# Primer manifest (config/primers.tsv) -- the spreadsheet the user edits -- and the
# `input_primers/` snapshot the pipeline actually reads.
#
# * Every primer gets a short, permanent internal id (P001, P002, ...; groups G001, ...) that is
#   what the whole pipeline, the output/ folders, notes, wet-lab records and group membership
#   are keyed by. The human-facing `primer_name` is just a label: it can be long, and renaming
#   it never orphans anything or triggers a recompute. Ids come from a counter in
#   config/id_registry.json that only ever goes up, so a deleted primer's id is never reused.
# * The pipeline never reads config/primers.tsv directly. sync_input_primers() copies each
#   ACTIVE (proceed=True) primer into input_primers/<id>.tsv containing only what affects the
#   computation (segment + the three sequences) and rewrites a file only when its content really
#   changed -- so Snakemake re-runs exactly the primers whose sequences changed, not everything
#   whenever any spreadsheet cell is touched.
# ---------------------------------------------------------------------------------------
PRIMER_MANIFEST_COLUMNS = [
    "primer_id", "primer_name", "segment", "forward_primer", "reverse_primer", "probe",
    "fluorophore", "source", "remark", "proceed", "add_date",
]
INPUT_PRIMER_COLUMNS = ["primer_id", "segment", "forward_primer", "reverse_primer", "probe"]
GROUP_COLUMNS = ["group_id", "label", "segment", "member_ids"]
ID_PREFIX = {"primer": "P", "group": "G"}


def _row_proceeds(row: dict) -> bool:
    return str(row.get("proceed", "")).strip().lower() in ("true", "1", "yes", "y")


def _read_tsv_rows(path: Path):
    if not Path(path).exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv_if_changed(path: Path, columns, rows) -> bool:
    """Writes a TSV only if the bytes would differ from what's already there (so an
    unchanged file keeps its mtime and Snakemake doesn't treat it as new input)."""
    import io
    buf = io.StringIO(newline="")
    writer = csv.DictWriter(buf, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for r in rows:
        writer.writerow({c: r.get(c, "") for c in columns})
    text = buf.getvalue()
    path = Path(path)
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return True


def _registry_path(manifest_path: Path) -> Path:
    return Path(manifest_path).parent / "id_registry.json"


def _id_number(value: str, prefix: str):
    value = (value or "").strip()
    if value.startswith(prefix) and value[len(prefix):].isdigit():
        return int(value[len(prefix):])
    return None


def allocate_ids(manifest_path: Path, kind: str, existing_ids, count: int = 1):
    """`count` fresh ids of `kind` ("primer"/"group"). The counter lives in
    config/id_registry.json and is also bumped past any id already in use (e.g. one typed in
    by hand), so an id is never handed out twice -- not even after its primer was deleted."""
    import json
    reg_path = _registry_path(manifest_path)
    reg = json.loads(reg_path.read_text(encoding="utf-8")) if reg_path.exists() else {}
    prefix = ID_PREFIX[kind]
    nxt = int(reg.get(f"{kind}_next", 1))
    used = [n for n in (_id_number(i, prefix) for i in existing_ids) if n is not None]
    nxt = max([nxt] + [n + 1 for n in used])
    ids = [f"{prefix}{nxt + i:03d}" for i in range(count)]
    reg[f"{kind}_next"] = nxt + count
    reg_path.write_text(json.dumps(reg, indent=1), encoding="utf-8")
    return ids


def load_primers_manifest_raw(path: Path):
    """Every row exactly as stored (including proceed=False ones) -- for the spreadsheet
    page, which needs to show and let someone re-enable a not-yet-ready candidate, not just
    the rows the pipeline currently acts on."""
    return _read_tsv_rows(path)


def write_primers_manifest(path: Path, rows):
    """Writes the spreadsheet back out. A row that has any content but no primer_id yet gets
    the next free id, and one with a name but no add_date gets today's date -- existing ids
    and dates are never touched, so this is safe to call on every save."""
    import datetime
    path = Path(path)
    today = datetime.date.today().isoformat()
    out_rows = [dict(r) for r in rows]
    need = [r for r in out_rows if not (r.get("primer_id") or "").strip()
            and any((r.get(k) or "").strip() for k in ("primer_name", "forward_primer", "reverse_primer", "probe"))]
    if need:
        new_ids = allocate_ids(path, "primer", [r.get("primer_id", "") for r in out_rows], len(need))
        for r, pid in zip(need, new_ids):
            r["primer_id"] = pid
    for r in out_rows:
        if (r.get("primer_name") or "").strip() and not (r.get("add_date") or "").strip():
            r["add_date"] = today
    _write_tsv_if_changed(path, PRIMER_MANIFEST_COLUMNS, out_rows)
    return out_rows


def ensure_primer_ids(path: Path):
    """Gives every row of a hand-edited manifest that has no id one (rewrites the file only if
    something was missing). Called when the pipeline starts, so a primer added straight into
    config/primers.tsv works exactly like one added through the webapp."""
    path = Path(path)
    rows = _read_tsv_rows(path)
    if any(not (r.get("primer_id") or "").strip() and any((r.get(k) or "").strip() for k in ("primer_name", "forward_primer"))
           for r in rows):
        return write_primers_manifest(path, rows)
    return rows


def load_groups_raw(path: Path):
    return _read_tsv_rows(path)


def write_groups(path: Path, rows):
    rows = [dict(r) for r in rows]
    need = [r for r in rows if not (r.get("group_id") or "").strip()]
    if need:
        new_ids = allocate_ids(Path(path), "group", [r.get("group_id", "") for r in rows], len(need))
        for r, gid in zip(need, new_ids):
            r["group_id"] = gid
    _write_tsv_if_changed(Path(path), GROUP_COLUMNS, rows)
    return rows


def sync_input_primers(manifest_path: Path, groups_path: Path, input_dir: Path):
    """Brings input_primers/ in line with the spreadsheet: one <id>.tsv per active primer (only
    segment + sequences), stale files for primers that are no longer active removed, and
    groups.tsv holding each group's still-active members only (a group with none left is
    dropped). Returns the list of active primer ids."""
    input_dir = Path(input_dir)
    input_dir.mkdir(parents=True, exist_ok=True)
    active = [r for r in _read_tsv_rows(manifest_path)
              if _row_proceeds(r) and (r.get("primer_id") or "").strip() and (r.get("primer_name") or "").strip()]
    for r in active:
        r = dict(r)
        r["segment"] = (r.get("segment") or "").strip().upper()
        _write_tsv_if_changed(input_dir / f"{r['primer_id'].strip()}.tsv", INPUT_PRIMER_COLUMNS, [r])
    active_ids = {r["primer_id"].strip() for r in active}
    for f in input_dir.glob("P*.tsv"):
        if f.stem not in active_ids:
            f.unlink()
    groups = []
    for g in load_groups_raw(groups_path):
        members = [m for m in (g.get("member_ids") or "").split(";") if m in active_ids]
        if members:
            groups.append({"group_id": g["group_id"], "label": g.get("label", ""),
                           "segment": (g.get("segment") or "").upper(), "member_ids": ";".join(members)})
    _write_tsv_if_changed(input_dir / "groups.tsv", GROUP_COLUMNS, groups)
    return sorted(active_ids)


def validate_primer_manifest_rows(rows, segments=("S", "L")):
    """Per-row validation for the spreadsheet's save endpoint -- returns a list of
    {row, field, message} problems (empty = all rows OK). Blank rows (no name and no
    sequences at all) are silently skipped rather than flagged, since the spreadsheet UI
    always keeps one trailing blank row to type a new entry into."""
    errors = []
    seen_names = {}
    seen_ids = {}
    valid_bases = set("ACGTUWSMKRYBDHVN")
    for i, r in enumerate(rows):
        name = (r.get("primer_name") or "").strip()
        f_seq = (r.get("forward_primer") or "").strip()
        r_seq = (r.get("reverse_primer") or "").strip()
        p_seq = (r.get("probe") or "").strip()
        if not name and not f_seq and not r_seq and not p_seq:
            continue  # a genuinely empty row -- not an error, just not a real entry yet
        pid = (r.get("primer_id") or "").strip()
        if pid:
            if _id_number(pid, "P") is None:
                errors.append({"row": i, "field": "primer_id", "message": "primer_id 형식이 올바르지 않습니다 (P001 형식)"})
            elif pid in seen_ids:
                errors.append({"row": i, "field": "primer_id", "message": f"primer_id가 중복됩니다 (행 {seen_ids[pid] + 1}과 동일)"})
            else:
                seen_ids[pid] = i
        if not name:
            errors.append({"row": i, "field": "primer_name", "message": "프라이머 이름이 비어 있습니다"})
        elif name in seen_names:
            errors.append({"row": i, "field": "primer_name", "message": f"이름이 중복됩니다 (행 {seen_names[name] + 1}과 동일)"})
        else:
            seen_names[name] = i
        seg = (r.get("segment") or "").strip().upper()
        if seg not in segments:
            errors.append({"row": i, "field": "segment", "message": f"segment는 {' 또는 '.join(segments)} 이어야 합니다"})
        for field, seq, label in (("forward_primer", f_seq, "Forward"), ("reverse_primer", r_seq, "Reverse"), ("probe", p_seq, "Probe")):
            if not seq:
                errors.append({"row": i, "field": field, "message": f"{label} 서열이 비어 있습니다"})
            elif set(seq.upper()) - valid_bases:
                bad = sorted(set(seq.upper()) - valid_bases)
                errors.append({"row": i, "field": field, "message": f"{label} 서열에 올바르지 않은 문자가 있습니다: {', '.join(bad)}"})
        proceed = (r.get("proceed") or "").strip().lower()
        if proceed not in ("true", "false", "1", "0", "yes", "no", "y", "n", ""):
            errors.append({"row": i, "field": "proceed", "message": "분석진행여부는 True/False여야 합니다"})
    return errors


def load_primers_manifest(path_or_paths):
    """Reads the manifest (config/primers.tsv) or one-or-more per-primer input_primers/<id>.tsv
    files and returns one dict per ACTIVE primer shaped the way the rest of the pipeline
    already expects: set_id (= the permanent primer id) / segment / source plus
    F_id/F_seq/R_id/R_seq/P_id/P_seq, and `name` (the human-facing label, empty for the
    per-primer files, which deliberately don't carry it). A row whose `proceed` column exists
    and isn't true is dropped; a row without an id can't be addressed and is skipped."""
    paths = [path_or_paths] if isinstance(path_or_paths, (str, Path)) else list(path_or_paths)
    out = []
    for path in paths:
        for r in _read_tsv_rows(Path(path)):
            if "proceed" in r and not _row_proceeds(r):
                continue
            pid = (r.get("primer_id") or "").strip()
            if not pid:
                continue
            out.append({
                "set_id": pid,
                "name": (r.get("primer_name") or "").strip(),
                "segment": (r.get("segment") or "").strip().upper(),
                "source": (r.get("source") or "").strip(),
                "F_id": f"{pid}_F", "F_seq": (r.get("forward_primer") or "").strip(),
                "R_id": f"{pid}_R", "R_seq": (r.get("reverse_primer") or "").strip(),
                "P_id": f"{pid}_P", "P_seq": (r.get("probe") or "").strip(),
                "fluorophore": (r.get("fluorophore") or "").strip(),
                "add_date": (r.get("add_date") or "").strip(),
            })
    return out


def write_ids(path: Path, ids):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for i in sorted(ids):
            fh.write(i + "\n")


def read_ids(path: Path):
    with open(path, encoding="utf-8") as fh:
        return {line.strip() for line in fh if line.strip()}


IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "D": "H", "H": "D", "B": "V", "V": "B",
}


def reverse_complement(seq: str) -> str:
    return "".join(_COMPLEMENT.get(c, "N") for c in reversed(seq.upper()))


def iupac_mismatch_count(query: str, target: str) -> int:
    """Count positions where `target` (plain ACGT, e.g. a genome's matched site) doesn't
    satisfy the ambiguity code in `query` (the original IUPAC primer/probe sequence)."""
    query, target = query.upper(), target.upper()
    n = min(len(query), len(target))
    mism = sum(1 for q, t in zip(query[:n], target[:n]) if t not in IUPAC.get(q, "ACGT"))
    mism += abs(len(query) - len(target))
    return mism


def direct_best_match(query: str, target_seq: str, both_strands: bool = True):
    """Best-matching window for `query` (an IUPAC oligo) directly within `target_seq` (plain
    ACGT), by brute-force scan rather than mfeprimer -- used as a fallback for genomes where
    mfeprimer's own search reports no hit at all despite a real, low-mismatch site existing.
    mfeprimer first finds candidates via a short seed (k-mer) match and only THEN applies its
    mismatch filter, so a genuinely low-mismatch site can be missed entirely if the true
    mismatch happens to fall inside whatever region it seeds on -- confirmed concretely (see
    CHANGELOG) for both a probe (mfeprimer `search -m 8` finding nothing for a genome with an
    actual 1-mismatch site) and, more rarely, a primer (`--virus --misMatch 2` finding no
    amplicon at all for a genome with an actual 0-mismatch site). Returns
    (mismatch, start_0based, strand) with strand "+" if `query` matched directly or "-" if
    its reverse complement matched (only tried when both_strands=True), or None if
    `target_seq` is shorter than `query`. No Tm/deltaG here -- those need mfeprimer's own
    thermodynamic calculation, which this intentionally doesn't try to replicate."""
    target_seq = target_seq.upper()
    candidates = [("+", query)]
    if both_strands:
        candidates.append(("-", reverse_complement(query)))
    best = None
    for strand, seq in candidates:
        n = len(seq)
        if len(target_seq) < n:
            continue
        for i in range(len(target_seq) - n + 1):
            mm = iupac_mismatch_count(seq, target_seq[i:i + n])
            if best is None or mm < best[0]:
                best = (mm, i, strand)
    return best


# F/R fallback tuning: matches the widest F/R search this pipeline already runs
# (--virus --misMatch 2), so the fallback never accepts a genome the primary search itself
# would have rejected on mismatch count alone -- it only rescues a genuinely-missed seed.
# Amplicon size bounds are calibrated from this pipeline's own observed amplicons (~24-173bp
# across all 25 sets' reference-genome binding), padded well past that on both sides so a
# real, slightly unusual amplicon isn't rejected, while still ruling out two independently
# well-matching but unrelated F/R hits landing far apart by chance.
FR_FALLBACK_MAX_MISMATCH = 2
FR_FALLBACK_MIN_AMP_SIZE = 15
FR_FALLBACK_MAX_AMP_SIZE = 2000    # same ceiling as mfeprimer's own --maxSize default


def direct_best_fr_match(f_seq: str, r_seq: str, genome_seq: str, three_prime_n: int = 0):
    """Fallback for a genome where mfeprimer's --virus geometric search found no amplicon at
    all, despite a real, low-mismatch F+R pair existing on that genome -- the same
    seed-finding-miss class direct_best_match() handles for a single oligo (see its
    docstring), but for the primer PAIR. --virus mode requires both F and R to independently
    seed before it will even attempt pairing them, so if either one's true mismatch happens
    to fall inside mfeprimer's own seed region, --virus reports no amplicon at all for that
    genome -- confirmed concretely for one genome, MH215287, in set
    VendorD_GPC_PREVIOUS1_LASV1_P2 (1/659 population genomes for that set).

    Independently brute-force-scans the WHOLE genome for each of F and R (both strands --
    cheap here since these genomes are only a few kb, unlike the population-wide search
    mfeprimer's seed index exists to keep tractable), then accepts the result only if the two
    hits form a geometrically valid amplicon: one on the top strand, the other on the bottom,
    correctly ordered (F upstream of R in top-strand coordinates), within a sane amplicon size
    range, and each within FR_FALLBACK_MAX_MISMATCH. Returns a dict shaped like
    parse_spec_best()'s per-chrom entries (amp_start/amp_end 1-based inclusive, matching
    build_alignment_table.py's convention), or None if no valid pair is found."""
    f_hit = direct_best_match(f_seq, genome_seq)
    r_hit = direct_best_match(r_seq, genome_seq)
    if not f_hit or not r_hit:
        return None
    f_mm, f_i, f_strand = f_hit
    r_mm, r_i, r_strand = r_hit
    if f_mm > FR_FALLBACK_MAX_MISMATCH or r_mm > FR_FALLBACK_MAX_MISMATCH:
        return None
    if f_strand == "+" and r_strand == "-":
        fp_i, fp_len = f_i, len(f_seq)
        rp_i, rp_len = r_i, len(r_seq)
    elif f_strand == "-" and r_strand == "+":
        fp_i, fp_len = r_i, len(r_seq)
        rp_i, rp_len = f_i, len(f_seq)
    else:
        return None  # both hit the same strand -- not a valid PCR orientation
    if fp_i >= rp_i:
        return None  # the "forward" side must sit upstream of the "reverse" side
    amp_start_0, amp_end_0 = fp_i, rp_i + rp_len  # 0-based, half-open
    amp_size = amp_end_0 - amp_start_0
    if not (FR_FALLBACK_MIN_AMP_SIZE <= amp_size <= FR_FALLBACK_MAX_AMP_SIZE):
        return None
    f_query = f_seq if f_strand == "+" else reverse_complement(f_seq)
    f_target = genome_seq[f_i:f_i + len(f_seq)].upper()
    r_query = r_seq if r_strand == "+" else reverse_complement(r_seq)
    r_target = genome_seq[r_i:r_i + len(r_seq)].upper()
    return {
        "f_mismatch": f_mm, "r_mismatch": r_mm, "total": f_mm + r_mm,
        "f_3p_ok": terminal_match_ok(f_query, f_target, f_strand == "+", three_prime_n),
        "r_3p_ok": terminal_match_ok(r_query, r_target, r_strand == "+", three_prime_n),
        "amp_start": str(amp_start_0 + 1), "amp_end": str(amp_end_0),
        # No Tm/deltaG here -- direct_best_match() doesn't run mfeprimer's own thermodynamic
        # calculation, same limitation noted on the probe fallback that already ships.
        "f_tm": "", "f_dg": "", "r_tm": "", "r_dg": "",
    }


def terminal_match_ok(query: str, target: str, is_fp: bool, n: int) -> bool:
    """True if the primer's own 3' end (the last n bases in the primer's own 5'->3' reading
    direction, NOT necessarily the last n characters of `query`) has zero mismatches against
    the matched genome sequence -- a mismatch right at the 3' end is disproportionately
    disruptive to polymerase extension (the base the polymerase is about to extend from),
    so even a primer that passes on total mismatch count can still fail here.

    `query`/`target` are exactly what iupac_mismatch_count() compares (see parse_spec_best):
    for an "fp" (direct-strand) hit, query is the oligo's own sequence unchanged, so its 3'
    end is query's own last character. For an "rp" (reverse-complemented) hit, query is
    reverse_complement(oligo_seq) -- reversing a string puts the oligo's own 3'-end base
    (complemented) at index 0, so the 3' end there is query's FIRST character instead."""
    n = min(n, len(query), len(target))
    if n <= 0:
        return True
    q, t = (query[-n:], target[-n:]) if is_fp else (query[:n], target[:n])
    return iupac_mismatch_count(q, t) == 0


# ---- optional Go engine for the whole-genome scan (config `scan_engine: python|go`, env SCAN_ENGINE overrides) ------------------
# The scans are results-identical in both engines (go/scan is compared with the Python code on every set, see scripts/engine_check.py).
# "go" works through a cache: prefetch_sites() runs go/scan once for all oligos x genomes of a step, and _oligo_sites() then answers from
# it. Anything not prefetched -- or any failure of the Go binary -- simply runs the Python code, so "go" can never change a result.
_ROOT = Path(__file__).resolve().parent.parent
_SITE_CACHE = {}
SCAN_ENGINES = ("python", "go")


def scan_engine(pipeline_root=None) -> str:
    name = (os.environ.get("SCAN_ENGINE") or read_config_str(pipeline_root or _ROOT, "scan_engine") or "python").strip().lower()
    return name if name in SCAN_ENGINES else "python"


def _go_scan_bin(pipeline_root=None):
    root = Path(pipeline_root or _ROOT)
    gdir = root / "go" / "scan"
    exe = gdir / "scan"
    src = max((p.stat().st_mtime for p in gdir.glob("*.go")), default=0)
    if exe.exists() and exe.stat().st_mtime >= src:
        return str(exe)
    go = shutil.which("go") or next((str(p) for p in (Path.home() / ".local/go/bin/go", Path("/usr/local/go/bin/go")) if p.exists()), None)
    if go and gdir.exists():
        r = subprocess.run([go, "build", "-o", str(exe), "."], cwd=gdir, capture_output=True, text=True, timeout=300)
        if r.returncode == 0:
            return str(exe)
    return str(exe) if exe.exists() else None


def clear_site_cache():
    _SITE_CACHE.clear()


def prefetch_sites(oligos, genome_seqs, pipeline_root=None) -> bool:
    """With the Go engine selected: find the sites of every (oligo, max_mm) in every genome in one parallel pass and keep them for
    _oligo_sites(). `oligos` = [(oligo_seq, max_mm)], `genome_seqs` = sequences. Returns True when the cache was filled."""
    import json as _json
    if scan_engine(pipeline_root) != "go":
        return False
    uniq_oligos = list(dict.fromkeys((o.upper(), int(mm)) for o, mm in oligos if o))
    genomes = list(dict.fromkeys(g.upper() for g in genome_seqs if g))
    if not uniq_oligos or not genomes:
        return False
    exe = _go_scan_bin(pipeline_root)
    if not exe:
        print("[warn] scan_engine=go but go/scan is not built and Go is not installed -- using the Python scan", file=sys.stderr)
        return False
    req = {"mode": "multi", "oligos": [{"seq": o, "max_mm": mm} for o, mm in uniq_oligos],
           "genomes": [{"id": str(i), "seq": g} for i, g in enumerate(genomes)]}
    try:
        r = subprocess.run([exe], input=_json.dumps(req), capture_output=True, text=True, timeout=3600)
        if r.returncode:
            raise RuntimeError(r.stderr.strip() or f"rc={r.returncode}")
        results = _json.loads(r.stdout)["results"]
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError, KeyError) as exc:
        print(f"[warn] go/scan failed ({exc}) -- using the Python scan", file=sys.stderr)
        return False
    for gi, res in enumerate(results):
        for k, (o, mm) in enumerate(uniq_oligos):
            _SITE_CACHE[(o, genomes[gi], mm)] = [(s["start"], s["end"], s["strand"], s["mm"]) for s in res["multi"][k]]
    return True


def _oligo_sites(oligo: str, genome: str, max_mm: int):
    """Every window of `genome` where `oligo` (IUPAC) has <= max_mm mismatches, on both strands:
    [(start0, end0_exclusive, strand, mismatches)]. Pure Python with an early exit per window."""
    if _SITE_CACHE:
        hit = _SITE_CACHE.get((oligo, genome, max_mm))
        if hit is not None:
            return list(hit)
    sites = []
    n = len(genome)
    for strand, q in (("+", oligo), ("-", reverse_complement(oligo))):
        sets = [IUPAC.get(c, "ACGT") for c in q]
        L = len(q)
        for i in range(n - L + 1):
            m = 0
            for j in range(L):
                if genome[i + j] not in sets[j]:
                    m += 1
                    if m > max_mm:
                        break
            else:
                sites.append((i, i + L, strand, m))
    return sites


def scan_fr(f_seq: str, r_seq: str, genome_seq: str, three_prime_n: int = 0, max_mm: int = 2,
            max_size: int = 2000):
    """Independent whole-genome F/R pair search (no mfeprimer): every site of each primer with <= max_mm
    mismatches on either strand, every pairing that forms a valid amplicon (opposite strands facing each other,
    <= max_size bp), best pair by (3' ends intact, total mismatches, shorter amplicon). Same dict shape as
    parse_spec_best() entries (without Tm/dG), or None."""
    g = genome_seq.upper()
    fs, rs = _oligo_sites(f_seq.upper(), g, max_mm), _oligo_sites(r_seq.upper(), g, max_mm)
    best = None
    for f in fs:
        for r in rs:
            if f[2] == "+" and r[2] == "-" and f[0] <= r[0] and r[1] > f[0] and r[1] - f[0] <= max_size:
                left, right = f, r
            elif r[2] == "+" and f[2] == "-" and r[0] <= f[0] and f[1] > r[0] and f[1] - r[0] <= max_size:
                left, right = r, f
            else:
                continue
            f_q = f_seq.upper() if f[2] == "+" else reverse_complement(f_seq.upper())
            r_q = r_seq.upper() if r[2] == "+" else reverse_complement(r_seq.upper())
            f_ok = terminal_match_ok(f_q, g[f[0]:f[1]], f[2] == "+", three_prime_n)
            r_ok = terminal_match_ok(r_q, g[r[0]:r[1]], r[2] == "+", three_prime_n)
            key = (0 if (f_ok and r_ok) else 1, f[3] + r[3], right[1] - left[0])
            if best is None or key < best[0]:
                best = (key, {
                    "f_mismatch": f[3], "r_mismatch": r[3], "total": f[3] + r[3], "f_3p_ok": f_ok, "r_3p_ok": r_ok,
                    "amp_start": str(left[0] + 1), "amp_end": str(right[1]),
                    "f_tm": "", "f_dg": "", "r_tm": "", "r_dg": "",
                    "f_pos": {"start": f[0] + 1, "end": f[1], "strand": f[2]},
                    "r_pos": {"start": r[0] + 1, "end": r[1], "strand": r[2]}})
    return best[1] if best else None


def scan_probe(probe_seq: str, genome_seq: str, amp_start: int = None, amp_end: int = None, max_mm: int = 8):
    """Independent whole-genome probe search (no mfeprimer): the best window of the genome (both strands, IUPAC-aware)
    that lies inside the amplicon [amp_start, amp_end] (1-based, inclusive) when one is given, else anywhere.
    -> {"mismatch", "start", "end" (1-based inclusive), "tm", "dg"} or None."""
    best = None
    for s0, e0, _strand, mm in _oligo_sites(probe_seq.upper(), genome_seq.upper(), max_mm):
        if amp_start is not None and not (amp_start <= s0 + 1 and e0 <= amp_end):
            continue
        if best is None or mm < best["mismatch"]:
            best = {"mismatch": mm, "start": str(s0 + 1), "end": str(e0), "tm": "", "dg": ""}
    return best


def choose_probe(mfe_hit, own_hit):
    """Better (fewer mismatches) of mfeprimer's best probe hit and the genome scan -> (hit, source); a tie keeps
    mfeprimer's (it carries Tm/dG)."""
    if mfe_hit is None and own_hit is None:
        return None, ""
    if mfe_hit is None:
        return own_hit, "genome-scan"
    if own_hit is not None and own_hit["mismatch"] < mfe_hit["mismatch"]:
        return own_hit, "genome-scan"
    return mfe_hit, "mfeprimer"


def _fr_key(x):
    return (0 if (x.get("f_3p_ok") is not False and x.get("r_3p_ok") is not False) else 1, x["total"])


def choose_fr(mfe_fr, own_fr):
    """Combine the two engines' best pair for one genome -> (pair, source). The genome scan wins only when it is
    strictly better (mfeprimer's seed missed a site); on a tie mfeprimer's entry is kept because it carries Tm/dG."""
    if mfe_fr is None and own_fr is None:
        return None, ""
    if mfe_fr is None:
        return own_fr, "genome-scan"
    if own_fr is not None and _fr_key(own_fr) < _fr_key(mfe_fr):
        return own_fr, "genome-scan"
    return mfe_fr, "mfeprimer"


def probe_inside_amplicon(row: dict) -> bool:
    """A TaqMan probe only counts when it sits between the primers. Rows lacking coordinates are not
    judged here (True), so older tables keep working."""
    try:
        a, b = int(row["amplicon_start"]), int(row["amplicon_end"])
        ps, pe = int(row["probe_start"]), int(row["probe_end"])
    except (KeyError, TypeError, ValueError):
        return True
    return a <= min(ps, pe) and max(ps, pe) <= b


def parse_spec_best(spec_path: Path, f_id: str, f_seq: str, r_id: str, r_seq: str, three_prime_n: int = 0,
                    genome_seqs: dict = None):
    """chrom -> {f_mismatch, r_mismatch, f_3p_ok, r_3p_ok, amp_start, amp_end, f_tm, f_dg,
    r_tm, r_dg}, keeping the best (lowest total mismatch) row per chrom.

    mfeprimer always reports "fp" as whichever of the pair sits left (lower coordinate,
    matches the top strand directly) and "rp" as whichever sits right (matches the top
    strand only after reverse-complementing) -- this is a geometric amplicon convention,
    not necessarily the same as which oligo we called F vs R, so match by name rather
    than assuming fpSeq<->F_seq. Shared by the Snakemake pipeline (build_alignment_table.py)
    and the interactive webapp's live analysis (live_analyze.py) so both compute mismatch
    counts and 3'-end status identically -- these used to be two separately-maintained
    copies of this same function, which is exactly how they drifted apart before.

    genome_seqs (chrom -> full sequence): when given, every mismatch count and 3'-end check is made against the
    GENOME bases at the reported coordinates. mfeprimer's own fpSeq/rpSeq column is the primer (a degenerate
    primer expanded into plain-ACGT variants) and only reflects mismatches inside the 3' k-mer it seeds on
    (--misStart/--misEnd, 9 bases by default): a mismatch further toward the 5' end is not in that text, so
    counting against it reports 0 where the genome really differs. Always pass genome_seqs when the sequences
    are at hand; without it the (under-counting) reported text is used, as before."""
    best = {}
    if not spec_path.exists():
        # mfeprimer only writes this sibling file when it finds >=1 amplicon -- against a
        # tiny db (e.g. the 2-4 record KDCA standard-material set) a primer pair can
        # legitimately find none at all, which is itself a valid "no match" result.
        return best
    with open(spec_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            # name chrom ampStart ampEnd ampGC ampSize fpName fpStart fpEnd fpSeq fpTm fpGC fpDg rpName rpEnd rpStart rpSeq ...
            chrom = cols[1]
            amp_start, amp_end = cols[2], cols[3]
            fp_name, fp_start, fp_end, fp_seq = cols[6], cols[7], cols[8], cols[9]
            rp_name, rp_end, rp_start, rp_seq = cols[13], cols[14], cols[15], cols[16]

            # fp = whichever oligo sits left (matches top strand directly, strand "+");
            # rp = whichever sits right (matches only after reverse-complementing, strand
            # "-") -- rp's own coordinate pair comes out of mfeprimer reversed (rpEnd before
            # rpStart, with rpStart > rpEnd) since it's numbered along the oligo's own 5'->3'
            # direction rather than top-strand order; min/max below undoes that for a plain
            # ascending (start, end) range matching fp's own convention.
            # Figure out which of our F/R that is by name, then store under the right label.
            fp_tm, fp_dg = cols[10], cols[12]
            rp_tm, rp_dg = cols[17], cols[19]

            mismatches = {}
            terminal_ok = {}
            thermo = {}
            positions = {}
            for name, seq, is_fp, tm, dg, pos_a, pos_b in (
                (fp_name, fp_seq, True, fp_tm, fp_dg, fp_start, fp_end),
                (rp_name, rp_seq, False, rp_tm, rp_dg, rp_end, rp_start),
            ):
                if name.startswith(f_id):
                    role, oligo_seq = "F", f_seq
                elif name.startswith(r_id):
                    role, oligo_seq = "R", r_seq
                else:
                    mismatches = None
                    break
                query = oligo_seq if is_fp else reverse_complement(oligo_seq)
                target = seq
                if genome_seqs is not None and chrom in genome_seqs:
                    lo, hi = min(int(pos_a), int(pos_b)), max(int(pos_a), int(pos_b))
                    target = genome_seqs[chrom][lo - 1:hi].upper()      # 1-based inclusive, top strand
                mismatches[role] = iupac_mismatch_count(query, target)
                terminal_ok[role] = terminal_match_ok(query, target, is_fp, three_prime_n)
                thermo[role] = {"tm": tm, "dg": dg}
                positions[role] = {"start": min(int(pos_a), int(pos_b)), "end": max(int(pos_a), int(pos_b)),
                                    "strand": "+" if is_fp else "-"}
            if not mismatches or "F" not in mismatches or "R" not in mismatches:
                continue  # unexpected primer name, skip defensively

            f_mm, r_mm = mismatches["F"], mismatches["R"]
            total = f_mm + r_mm
            prev = best.get(chrom)
            if prev is None or total < prev["total"]:
                best[chrom] = {"f_mismatch": f_mm, "r_mismatch": r_mm, "total": total,
                                "f_3p_ok": terminal_ok["F"], "r_3p_ok": terminal_ok["R"],
                                "amp_start": amp_start, "amp_end": amp_end,
                                "f_tm": thermo["F"]["tm"], "f_dg": thermo["F"]["dg"],
                                "r_tm": thermo["R"]["tm"], "r_dg": thermo["R"]["dg"],
                                "f_pos": positions["F"], "r_pos": positions["R"]}
    return best
