"""Host of each genome -> host class, and the host SCOPES the results can be viewed in.

classes   human | animal | unknown | lab
            lab = a laboratory / cell-culture host (metadata flag): such genomes are excluded from every scope
scopes    human | animal | unknown | all      (all = human + animal + unknown, never lab)

The primary scope is "human" (diagnostic assays are meant for patient samples); a virus whose metadata has no human genome at all
falls back to "all". The pass/fail logic of one genome (primer + probe rules) lives here too, so the pipeline (reconcile.py) and
the web app count with the SAME code when they restrict a result to a scope.

This file is shared verbatim by the LASV and Nipah pipelines (and copied into ADHOC).
"""
import csv
import re
from pathlib import Path

SCOPES = ("human", "animal", "unknown", "all")
CLASSES = ("human", "animal", "unknown", "lab")
SCOPE_LABEL_KO = {"human": "사람", "animal": "동물", "unknown": "미상", "all": "전체", "lab": "실험 숙주"}
NCBI_HUMAN_TAXID = 9606

_UNKNOWN = {"", "unknown", "na", "n/a", "nan", "none", "null", "missing", "not applicable", "not collected", "not available",
            "not provided", "unspecified", "-", "?", "unknown/missing"}
_HUMAN_EXACT = {"homo sapiens", "human", "humans", "patient", "patients", "homo sapiens sapiens", "man", "woman", "person"}
_HUMAN_WORD = re.compile(r"(?<![a-z-])(human|homo sapiens|patient)(?![a-z])")


def _truthy(v):
    return str(v).strip().lower() in ("1", "true", "t", "yes", "y")


def classify(host, lab=False, taxid=None):
    """Host class of one genome from its metadata. `taxid` (NCBI host tax id) wins over the free text when given."""
    if _truthy(lab):
        return "lab"
    if taxid not in (None, ""):
        try:
            if int(taxid) == NCBI_HUMAN_TAXID:
                return "human"
        except (TypeError, ValueError):
            pass
    h = re.sub(r"\s+", " ", str(host or "").strip().lower())
    if h in _UNKNOWN:
        return "unknown"
    if "non-human" in h or "nonhuman" in h or "non human" in h:
        return "animal"
    if h in _HUMAN_EXACT or h.startswith("homo sapiens") or _HUMAN_WORD.search(h):
        return "human"
    return "animal"


def in_scope(cls, scope):
    """Is a genome of host class `cls` part of `scope`? Lab hosts are never in any scope."""
    if cls == "lab" or cls is None:
        return False
    return True if scope == "all" else cls == scope


def load_host_classes(path, columns):
    """{accession: class} from a source metadata TSV. `columns` is the profile's metadata mapping (accession column, optional
    'host' and 'lab_host' columns; defaults 'host' / 'is_lab_host'). Missing file or columns -> every genome is 'unknown'."""
    path = Path(path)
    if not path.exists():
        return {}
    acc_col = (columns or {}).get("accession", "accession")
    host_col = (columns or {}).get("host", "host")
    lab_col = (columns or {}).get("lab_host", "is_lab_host")
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get(acc_col) or "").strip()
            if acc:
                out[acc] = classify(row.get(host_col), row.get(lab_col))
    return out


def host_counts(ids, classes):
    """{class: n} over `ids` (accessions without a metadata row count as unknown)."""
    out = {c: 0 for c in CLASSES}
    for i in ids:
        out[classes.get(i, "unknown")] += 1
    return out


def primary_scope(ids, classes):
    """'human' when at least one genome is human, else 'all' (e.g. a virus whose metadata carries no host)."""
    return "human" if any(classes.get(i) == "human" for i in ids) else "all"


def ids_in_scope(ids, classes, scope):
    return [i for i in ids if in_scope(classes.get(i, "unknown"), scope)]


# ---- pass/fail of one genome (moved here from reconcile.py, unchanged in meaning) ------------------------------------------------
def _num(v):
    if v is None or v == "":
        return None
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _three_prime_ok(v):
    """blank/None = not checked (do not fail); '1'/1/True = ok; anything else = a mismatch at the primer's own 3' end."""
    if v is None or v == "":
        return True
    return v in ("1", 1, True) or str(v).strip().lower() in ("1", "true")


def primer_ok(row, mm):
    """F and R both found with at most `mm` mismatches and no mismatch at either primer's 3' end."""
    f, r = _num(row.get("f_mismatch")), _num(row.get("r_mismatch"))
    if f is None or r is None or f > mm or r > mm:
        return False
    return _three_prime_ok(row.get("f_3p_ok")) and _three_prime_ok(row.get("r_3p_ok"))


def probe_ok(row, mm, inside=None):
    """The probe found with at most `mm` mismatches AND lying between the primers (a probe binding elsewhere gives no signal)."""
    p = _num(row.get("p_mismatch"))
    if p is None or p > mm:
        return False
    if inside is None:
        from mfe_common import probe_inside_amplicon  # noqa: PLC0415 -- only needed when the caller did not pass its own check
        inside = probe_inside_amplicon
    return inside(row)


def genome_flags(row, primer_mm_pass, probe_mm, inside=None):
    return {
        "mm0": primer_ok(row, 0), "mm1": primer_ok(row, 1), "mm2": primer_ok(row, 2),
        "primer": primer_ok(row, primer_mm_pass), "probe": probe_ok(row, probe_mm, inside),
    }


def summarize(ids, rows_by_id, primer_mm_pass, probe_mm, inside=None):
    """Coverage numbers for the genomes `ids` (same fields and rounding as reconcile.py's summary row), plus the id sets."""
    flags = {i: genome_flags(rows_by_id.get(i, {}), primer_mm_pass, probe_mm, inside) for i in ids}
    n = len(flags)
    pct = lambda k: round(100 * k / n, 2) if n else 0.0  # noqa: E731
    primer = {i for i, f in flags.items() if f["primer"]}
    probe = {i for i, f in flags.items() if f["probe"]}
    validated = primer & probe
    dropout = primer - probe
    return {
        "total_genomes": n,
        "primer_coverage_mm0_pct": pct(sum(f["mm0"] for f in flags.values())),
        "primer_coverage_mm1_pct": pct(sum(f["mm1"] for f in flags.values())),
        "primer_coverage_mm2_pct": pct(sum(f["mm2"] for f in flags.values())),
        "probe_coverage_mm0_pct": pct(len(probe)),
        "primer_matched": len(primer), "probe_dropout": len(dropout), "validated": len(validated), "validated_pct": pct(len(validated)),
        "_validated_ids": validated, "_dropout_ids": dropout, "_primer_ids": primer, "_failed_ids": set(flags) - primer,
    }
