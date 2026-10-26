"""DRAFT link between the segments of a segmented virus (Lassa: S and L), built from the results -- the analysis pipeline is not touched.

Each segment is analysed on its own (its own population, its own primer sets). Here the records of the SAME virus are tied together by
strain name, so one "genome" = one virus with a record per segment, and the per-segment results can be read together: which viruses an
S set and an L set catch, which ones the pair together catches, and which ones nothing catches.

The link is a DRAFT: a plain, exact (case- and spacing-insensitive) match of the metadata `strain` column. It is checked against the
records' own data (country, lineage) and every doubtful case is listed instead of guessed:
    linked      exactly one record in every segment's analysed population
    ambiguous   several records of one strain in some segment (all of them are kept; any validated record counts for the virus)
    <seg>_only  the strain has a record in this segment only (it cannot be judged by the other segment's primer sets)
Sequences are NOT joined: the segments are separate RNA molecules, a concatenation would let primers pair across the junction.
"""
import csv
import re
from pathlib import Path

_CACHE = {}


def _stamp(p):
    p = Path(p)
    return p.stat().st_mtime_ns if p.exists() else None


def norm(strain):
    """Key of a strain name: case and runs of spaces do not matter, nothing else is guessed."""
    return re.sub(r"\s+", " ", str(strain or "").strip()).casefold()


def _read_raw(path):
    p = Path(path)
    if not p.exists():
        return {}
    with open(p, newline="", encoding="utf-8") as fh:
        return {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t") if r.get("accession")}


def build(root, segs, ids_by_seg, lineage_by_seg):
    """segs: ordered segment keys; ids_by_seg[seg]: accessions of the analysed population (in the host scope);
    lineage_by_seg[seg]: {accession: lineage}. -> {"segments", "viruses": [...], "counts": {...}, "conflicts": [...]}"""
    root = Path(root)
    raw = {s: _read_raw(root / "data" / "raw" / f"{s}_metadata.tsv") for s in segs}
    key = (str(root), tuple(segs), tuple(_stamp(root / "data" / "raw" / f"{s}_metadata.tsv") for s in segs),
           tuple((s, hash(frozenset(ids_by_seg[s])), hash(tuple(sorted(lineage_by_seg[s].items())))) for s in segs))
    if key in _CACHE:
        return _CACHE[key]
    by_name = {}                                           # normalised strain -> {"name": first spelling, "acc": {seg: [accessions]}}
    for s in segs:
        for acc in sorted(ids_by_seg[s]):
            r = raw[s].get(acc)
            if not r or not norm(r.get("strain")):
                continue
            v = by_name.setdefault(norm(r["strain"]), {"name": r["strain"].strip(), "acc": {x: [] for x in segs}})
            v["acc"][s].append(acc)
    viruses, conflicts = [], []
    for k, v in sorted(by_name.items()):
        have = [s for s in segs if v["acc"][s]]
        if len(have) == len(segs):
            status = "ambiguous" if any(len(v["acc"][s]) > 1 for s in segs) else "linked"
        else:
            status = f"{have[0]}_only"
        first = {s: (v["acc"][s][0] if v["acc"][s] else None) for s in segs}
        lineages = {s: lineage_by_seg[s].get(first[s], "NA") for s in segs if first[s]}
        countries = {s: (raw[s].get(first[s], {}).get("country") or "").strip() for s in segs if first[s]}
        lin = next((x for x in lineages.values() if x and x not in ("NA", "unclassified")), next(iter(lineages.values()), "NA"))
        country = next((c for c in countries.values() if c), "")
        evidence = {"country_agrees": (len({c for c in countries.values() if c}) <= 1), "lineage_agrees": (len({x for x in lineages.values() if x not in ("NA", "unclassified")}) <= 1)}
        if status in ("linked", "ambiguous") and not (evidence["country_agrees"] and evidence["lineage_agrees"]):
            conflicts.append({"virus_id": k, "name": v["name"], "countries": countries, "lineages": lineages})
        viruses.append({"virus_id": k, "name": v["name"], "acc": v["acc"], "status": status, "lineage": lin, "country": country, "evidence": evidence})
    counts = {"viruses": len(viruses)}
    for s in ("linked", "ambiguous"):
        counts[s] = sum(1 for x in viruses if x["status"] == s)
    for s in segs:
        counts[f"{s}_only"] = sum(1 for x in viruses if x["status"] == f"{s}_only")
    counts["records"] = {s: len(ids_by_seg[s]) for s in segs}
    counts["conflicts"] = len(conflicts)
    out = {"segments": list(segs), "viruses": viruses, "counts": counts, "conflicts": conflicts}
    if len(_CACHE) > 20:
        _CACHE.clear()
    _CACHE[key] = out
    return out


def to_tsv(link):
    segs = link["segments"]
    head = ["virus_id", "strain", "status", "lineage", "country"] + [f"{s}_accessions" for s in segs] + ["country_agrees", "lineage_agrees"]
    lines = ["\t".join(head)]
    for v in link["viruses"]:
        lines.append("\t".join([v["virus_id"], v["name"], v["status"], v["lineage"], v["country"]] + [";".join(v["acc"][s]) for s in segs]
                               + [str(v["evidence"]["country_agrees"]).lower(), str(v["evidence"]["lineage_agrees"]).lower()]))
    return "\n".join(lines) + "\n"


def universe(link):
    """The viruses a combined (all-segment) judgement can be made on: a record in every segment (linked + ambiguous)."""
    return [v for v in link["viruses"] if v["status"] in ("linked", "ambiguous")]


def map_ids(link, seg, validated_ids):
    """Virus ids (of the combined universe) for which at least one record of segment `seg` is in `validated_ids`."""
    vs = set(validated_ids)
    return sorted(v["virus_id"] for v in universe(link) if any(a in vs for a in v["acc"][seg]))
