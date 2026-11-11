"""Where the analysed genomes come from -- shown with every genome table and on the overview (shared verbatim by every edition).

Two levels:
  profile level  info(profile, root, read_tsv)  -> {"name", "detail", "url", "upstream", "fetched_at", "records", "by_record"}
  record level   record_source(row, profile)    -> "Pathoplexus . GenBank", "Nextstrain . RefSeq", "NCBI Datasets . GenBank" ...
The record label uses the `data_source` column the fetch scripts write (Pathoplexus / NCBI Datasets: the upstream database of that record); records
fetched before that column existed fall back to what the accession and the profile say.
"""
from datetime import datetime
from pathlib import Path

NAMES = {"nextstrain": "Nextstrain", "pathoplexus": "Pathoplexus", "ncbi": "NCBI Datasets"}


def _upstream(acc):
    a = str(acc or "")
    if a.startswith(("NC_", "NZ_", "NG_", "NM_")):
        return "RefSeq"
    if a.startswith("PP_"):
        return ""                                               # a Pathoplexus-only submission (no INSDC accession)
    return "GenBank"


def record_source(row, profile):
    """Source label of one genome record (row = the raw metadata row)."""
    v = (row.get("data_source") or "").strip()
    if v:
        return v
    name = NAMES.get(profile["source"]["type"], profile["source"]["type"])
    up = _upstream(row.get("accession"))
    return f"{name} · {up}" if up else f"{name} · 자체 등록"


def info(profile, root, read_tsv, seg_keys):
    s = profile["source"]
    kind = s["type"]
    name = NAMES.get(kind, kind)
    if kind == "pathoplexus":
        detail = f"organism {s['organism']}" + (f" (lineage = {s['lineage_field']})" if s.get("lineage_field") else "")
        url = f"https://pathoplexus.org/{s['organism']}/search"
        upstream = "INSDC (GenBank / ENA / DDBJ) + Pathoplexus submissions"
    elif kind == "ncbi":
        taxa = s["taxon"] if isinstance(s["taxon"], list) else [s["taxon"]]
        detail = "taxon " + ", ".join(str(t) for t in taxa) + (f" (collected since {s['since']})" if s.get("since") else "")
        url = f"https://www.ncbi.nlm.nih.gov/datasets/taxonomy/{taxa[0]}/"
        upstream = "NCBI GenBank / RefSeq"
    else:
        detail = s["base_url"].replace("https://", "")
        url = profile.get("ui", {}).get("web_url") or s["base_url"]
        upstream = "NCBI GenBank / RefSeq (via Nextstrain)"
    fetched, records, by, upd, rel = None, 0, {}, "", ""
    for seg in seg_keys:
        p = Path(root) / "data" / "raw" / f"{seg}_metadata.tsv"
        if not p.exists():
            continue
        t = datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds")
        fetched = max(fetched, t) if fetched else t
        for r in read_tsv(p):
            lab = record_source(r, profile)
            by[lab] = by.get(lab, 0) + 1
            records += 1
            upd = max(upd, (r.get("date_updated") or "")[:10])           # the source's own dates: when a record was last updated / first released
            rel = max(rel, (r.get("date_released") or "")[:10])
    return {"name": name, "detail": detail, "url": url, "upstream": upstream, "fetched_at": fetched, "latest_update": upd or None, "latest_release": rel or None, "records": records,
            "by_record": dict(sorted(by.items(), key=lambda kv: -kv[1]))}
