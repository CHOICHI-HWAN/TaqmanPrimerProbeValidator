"""Species/virus profile: everything the pipeline needs to know about *which virus* it is analysing.

The profile (profiles/<id>.json) is the one place that says where the data comes from, how its
metadata columns map onto the pipeline's canonical table (accession, lineage, date, country, region),
what the segments are (length floor, lineage column, reference genome), and the species-specific
quality thresholds. Rules and scripts receive these values as parameters instead of hard-coding them,
so another virus is a new profile file rather than a code change. Only the standard library is used
(JSON) so the web server -- which runs on a bare python3 -- can read the same file.
"""
import json
from pathlib import Path

PROFILE_VERSION = 1
CANONICAL_METADATA = ("accession", "lineage_fallback", "date", "country", "region")
LINEAGE_MODES = ("metadata_column", "auspice_clades")
QC_DEFAULTS = {"max_n_fraction": 0.0, "exclude_list": None}
DIVERSITY_DEFAULTS = {"window_size": 1, "min_coverage": 0.95, "min_identity": 0.65, "min_lineage_n": 10}
SOURCE_TYPES = ("nextstrain",)


class ProfileError(ValueError):
    pass


def _need(cond, msg):
    if not cond:
        raise ProfileError(msg)


def validate(p):
    _need(isinstance(p, dict), "profile must be a JSON object")
    _need(p.get("profile_version") == PROFILE_VERSION, f"profile_version must be {PROFILE_VERSION}")
    for k in ("id", "name", "short_name"):
        _need(isinstance(p.get(k), str) and p[k].strip(), f"'{k}' is required")
    _need(p["id"].replace("_", "").replace("-", "").isalnum(), "'id' may only contain letters, digits, - and _")
    src = p.get("source") or {}
    _need(src.get("type") in SOURCE_TYPES, f"source.type must be one of {SOURCE_TYPES}")
    _need(isinstance(src.get("base_url"), str) and src["base_url"].startswith("http"), "source.base_url must be a URL")
    src["segment_subdirs"] = bool(src.get("segment_subdirs", True))
    for f in ("sequences_file", "metadata_file"):
        src.setdefault(f, "sequences.fasta.zst" if f == "sequences_file" else "metadata.tsv.zst")
    p["source"] = src
    segs = p.get("segments")
    _need(isinstance(segs, dict) and segs, "'segments' needs at least one segment")
    for key, s in segs.items():
        _need(key.isalnum() and key == key.lower(), f"segment key '{key}' must be lowercase letters/digits")
        _need(isinstance(s.get("min_len"), int) and s["min_len"] > 0, f"segments.{key}.min_len must be a positive integer")
        for f in ("lineage_col", "reference_accession", "label"):
            _need(isinstance(s.get(f), str) and s[f].strip(), f"segments.{key}.{f} is required")
    lin = {"label": "Lineage", "mode": "metadata_column", **(p.get("lineage") or {})}
    _need(lin["mode"] in LINEAGE_MODES, f"lineage.mode must be one of {LINEAGE_MODES}")
    if lin["mode"] == "auspice_clades":
        # clades come from a locally built Auspice JSON (nextstrain build); {seg} is replaced by the segment key
        _need(isinstance(lin.get("auspice"), str) and lin["auspice"].strip(), "lineage.auspice (path to the Auspice JSON) is required for lineage.mode=auspice_clades")
        lin.setdefault("attr", "clade_membership")
        lin.setdefault("missing", "unclassified")
    p["lineage"] = lin
    qc = {**QC_DEFAULTS, **(p.get("qc") or {})}
    _need(0 <= qc["max_n_fraction"] < 1, "qc.max_n_fraction must be within 0..1")
    p["qc"] = qc
    md = p.get("metadata") or {}
    for f in CANONICAL_METADATA:
        _need(isinstance(md.get(f), str) and md[f].strip(), f"metadata.{f} (source column name) is required")
    pop = p.get("population") or {}
    for f in ("pooled_label", "chrom_prefix"):
        _need(isinstance(pop.get(f), str) and pop[f].strip(), f"population.{f} is required")
    div = {**DIVERSITY_DEFAULTS, **(p.get("diversity") or {})}
    _need(div["window_size"] >= 1 and div["min_lineage_n"] >= 1, "diversity.window_size / min_lineage_n must be >= 1")
    _need(0 <= div["min_coverage"] <= 1 and 0 <= div["min_identity"] <= 1, "diversity.min_coverage / min_identity must be within 0..1")
    p["diversity"] = div
    # wording shown in the web app / reports (optional; every key has a neutral default)
    ui = {"name_ko": p["short_name"], "web_url": p["source"]["base_url"], "lineage_ko": "리니지",
          "lineage_examples": "", "lineage_def": "", "orientation_note": "", **(p.get("ui") or {})}
    p["ui"] = ui
    return p


def load_profile(path):
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise ProfileError(f"profile not found: {path}") from None
    except json.JSONDecodeError as e:
        raise ProfileError(f"{path.name} is not valid JSON: {e}") from None
    return validate(data)


def find_profile_path(root):
    """Path of the profile selected by `profile:` in <root>/config/config.yaml (read with a regex so the
    web server and helper scripts do not need a YAML library). Falls back to profiles/lassa.json."""
    import re
    root = Path(root)
    cfg = root / "config" / "config.yaml"
    text = cfg.read_text(encoding="utf-8") if cfg.exists() else ""
    m = re.search(r'^profile:\s*"?([^"\n#]+?)"?\s*(?:#.*)?$', text, re.MULTILINE)
    return root / (m.group(1).strip() if m else "profiles/lassa.json")


def load_for_root(root):
    return load_profile(find_profile_path(root))


def source_urls(p, seg):
    """(sequences_url, metadata_url) for one segment, honouring source.segment_subdirs."""
    s = p["source"]
    prefix = s["base_url"].rstrip("/") + (f"/{seg}" if s["segment_subdirs"] else "")
    return f"{prefix}/{s['sequences_file']}", f"{prefix}/{s['metadata_file']}"


def segment_keys(p):
    return list(p["segments"])


def chrom_name(p, seg):
    """Contig name used in the all-sites VCF (e.g. LASV_S)."""
    return f"{p['population']['chrom_prefix']}_{seg.upper()}"
