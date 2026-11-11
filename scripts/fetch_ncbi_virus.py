#!/usr/bin/env python3
"""Fetch every genome of one virus taxon from the NCBI Datasets API (public, no key) and write them in the shape the pipeline's other fetch steps
produce: a FASTA (header = accession without version) and a metadata TSV (accession, date, country, region, host, is_lab_host, ...).

usage: fetch_ncbi_virus.py <datasets_api_url> <taxon_id>[/<SEGMENT>] <out.fasta> <out_metadata.tsv>      e.g. https://api.ncbi.nlm.nih.gov/datasets/v2 186539 ...
Several taxa at once with filters (influenza A: H1N1 + H3N2 + H5N1, one NCBI segment number, collected since a year, host filter per taxon):
  <api> "114727,119210,102793/4?since=2015&host=114727:9606,119210:9606&sub=H3&group=subtype" out.fasta out.tsv
(`sub` keeps the records whose subtype has that HA or NA part, `group=subtype` puts the subtype into the lineage column instead of the country).
A segmented virus is fetched one segment at a time (`11588/S`): the record's `segment` field decides; a record without it is placed by its length when that
is close to a full segment (SEG_LEN), so unlabeled complete segments are not lost.
The extra column `country_lineage` (= country) lets a virus without a clade scheme group its genomes by the country they were collected in.
"""
import csv
import fcntl
import gzip
import hashlib
import io
import json
import os
import re
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

AFRICA = {"Democratic Republic of the Congo", "Uganda", "Sierra Leone", "Guinea", "Liberia", "Gabon", "Republic of the Congo", "Sudan", "South Sudan",
          "Kenya", "Nigeria", "Mali", "Senegal", "Ivory Coast", "Cote d'Ivoire", "Cameroon", "Central African Republic", "Tanzania", "Angola", "Ethiopia", "Rwanda"}


def _get(url, tries=3):
    last = None
    for _ in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "primer-validation-pipeline"}), timeout=600) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
    raise SystemExit(f"fetch failed: {url[:120]}: {last}")


SEG_LEN = {"11588": {"S": (1600, 1800), "M": (3750, 4050), "L": (6250, 6550)}}      # Rift Valley fever virus: full-length windows of an unlabeled record


def main(base, taxon, out_fa, out_meta):
    base = base.rstrip("/")
    taxon, _, seg = taxon.partition("/")
    reports, token = [], ""
    while True:
        d = json.loads(_get(f"{base}/virus/taxon/{taxon}/dataset_report?page_size=1000" + (f"&page_token={token}" if token else "")))
        reports += d.get("reports", [])
        token = d.get("next_page_token") or ""
        if not token:
            break
    z = zipfile.ZipFile(io.BytesIO(_get(f"{base}/virus/taxon/{taxon}/genome/download?include_sequence=GENOME")))
    fasta = z.read("ncbi_dataset/data/genomic.fna").decode("utf-8")
    seqs, name = {}, None
    for ln in fasta.splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0].split(".")[0]
            seqs[name] = []
        elif name:
            seqs[name].append(ln.strip())
    cols = ["accession", "strain", "date", "date_released", "date_updated", "region", "country", "division", "location", "length", "host", "taxid",
            "is_lab_host", "completeness", "authors", "institution", "submitter_country", "sample_type", "country_lineage", "url"]
    n = 0
    with open(out_fa, "w", encoding="utf-8") as fa, open(out_meta, "w", encoding="utf-8", newline="") as fm:
        w = csv.writer(fm, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in reports:
            acc = r["accession"].split(".")[0]
            if acc not in seqs:
                continue
            if seg:
                lo, hi = SEG_LEN.get(taxon, {}).get(seg, (0, 0))
                have = (r.get("segment") or "").strip().upper()
                if have != seg and not (not have and lo <= int(r.get("length") or 0) <= hi):
                    continue
            fa.write(f">{acc}\n" + "\n".join(seqs[acc]) + "\n")
            n += 1
            iso, loc, host = r.get("isolate") or {}, r.get("location") or {}, r.get("host") or {}
            geo = (loc.get("geographic_location") or "").strip()
            country, _, rest = geo.partition(":")
            country = country.strip()
            if country.lower() in ("missing", "not collected", "not available", "na", "unknown"):
                country = ""
            sub = r.get("submitter") or {}
            w.writerow([acc, iso.get("name") or "", iso.get("collection_date") or "", (r.get("release_date") or "")[:10], (r.get("update_date") or "")[:10],
                        "Africa" if country in AFRICA else "", country, rest.strip(), "", r.get("length") or "", host.get("organism_name") or "", host.get("tax_id") or "",
                        "false", r.get("completeness") or "", "; ".join(sub.get("names") or []), sub.get("affiliation") or "", sub.get("country") or "", "",
                        country or "unknown", f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"])
    print(f"taxon {taxon}: {len(reports)} reports, {n} sequences written")


# ---- several taxa with filters (influenza A) ----------------------------------------------------------------------------------------------
REPORT_KEEP = ("accession", "isolate", "host", "location", "length", "completeness", "segment", "release_date", "update_date", "source_database", "submitter")


def _post_zip(url, body, tries=3):
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(url, data=json.dumps(body).encode(), headers={"Content-Type": "application/json", "Accept": "application/zip",
                                                                                      "User-Agent": "primer-validation-pipeline"})
            with urllib.request.urlopen(req, timeout=600) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(3)
    raise SystemExit(f"fetch failed: {url[:100]}: {last}")


def _subtype(r):
    name = (r.get("virus") or {}).get("organism_name") or ""
    m = re.search(r"\((H\d{1,2}N\d{1,2})\)\)?\s*$", name)
    if m:
        return m.group(1)
    for n in (r.get("virus") or {}).get("lineage") or []:
        m = re.match(r"^(H\d{1,2}N\d{1,2}) subtype$", n.get("name") or "")
        if m:
            return m.group(1)
    return "unclassified"


def _year(r):
    d = (r.get("isolate") or {}).get("collection_date") or ""
    if d[:4].isdigit():
        return int(d[:4])
    m = re.search(r"/(\d{4})(?:\(|$)", (r.get("isolate") or {}).get("name") or "")
    if m:
        return int(m.group(1))
    rd = (r.get("release_date") or "")[:4]
    return int(rd) if rd.isdigit() else 0


def _cached_reports(base, taxa, hostmap, since, cache_dir):
    """All dataset reports of the taxa (host filter per taxon, released since the year), cached for 12 h so the per-segment jobs share one download."""
    key = hashlib.sha1(json.dumps([base, taxa, sorted(hostmap.items()), since]).encode()).hexdigest()[:12]
    cache = Path(cache_dir) / f".ncbi_reports_{key}.json.gz"
    lock = Path(cache_dir) / f".ncbi_reports_{key}.lock"
    with open(lock, "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)                                   # one job downloads, the others wait and reuse
        if cache.exists() and time.time() - cache.stat().st_mtime < 12 * 3600:
            return json.loads(gzip.open(cache, "rt", encoding="utf-8").read())
        out = []
        for t in taxa:
            token, flt = "", ""
            if hostmap.get(t):
                flt += f"&filter.host={hostmap[t]}"
            if since:
                flt += f"&filter.released_since={since}-01-01T00:00:00.000Z"
            while True:
                d = json.loads(_get(f"{base}/virus/taxon/{t}/dataset_report?page_size=1000{flt}" + (f"&page_token={token}" if token else "")))
                for r in d.get("reports", []):
                    x = {k: r[k] for k in REPORT_KEEP if k in r}
                    x["_subtype"], x["_year"] = _subtype(r), _year(r)
                    out.append(x)
                token = d.get("next_page_token") or ""
                if not token:
                    break
        tmp = cache.with_suffix(".tmp")
        with gzip.open(tmp, "wt", encoding="utf-8") as fh:
            fh.write(json.dumps(out))
        os.replace(tmp, cache)
        return out


# ---- influenza A: subtype of records that carry none, and a capped, strain-consistent sample ---------------------------------------------------
HA_REFS = {"H1": "NC_026433.1", "H3": "NC_007366.1", "H5": "NC_007362.1"}              # reference HA of the subtypes that are kept
HA_SUBTYPE = {"H1": "H1N1", "H3": "H3N2", "H5": "H5N1"}                               # human/avian seasonal-type assumption: HA decides the subtype name
KEEP_SUBTYPES = ("H1N1", "H3N2", "H5N1")


def _fasta_of(base, accessions):
    out = {}
    for i in range(0, len(accessions), 2000):
        z = zipfile.ZipFile(io.BytesIO(_post_zip(f"{base}/virus/genome/download", {"accessions": accessions[i:i + 2000], "include_sequence": ["GENOME"]})))
        name = None
        for ln in z.read("ncbi_dataset/data/genomic.fna").decode("utf-8").splitlines():
            if ln.startswith(">"):
                name = ln[1:].split()[0]
                out[name] = []
            elif name:
                out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def _blastn():
    import shutil
    import subprocess  # noqa: F401
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from standard_placement import find_blastn
        exe = find_blastn(Path(__file__).resolve().parent.parent)
    except Exception:  # noqa: BLE001
        exe = None
    return exe or shutil.which("blastn")


def _classify_ha(base, ha_records, tmp_dir):
    """{accession (versioned): 'H1'|'H3'|'H5'} by megablast of the HA sequences against the three reference HAs."""
    import subprocess
    exe = _blastn()
    if not exe:
        raise SystemExit("blastn not found (needed to name the subtype of records that carry none): set $BLASTN_BIN or install blast")
    refs = _fasta_of(base, list(HA_REFS.values()))
    sub = Path(tmp_dir) / ".flu_ha_refs.fasta"
    name_to_h = {v.split(".")[0]: k for k, v in HA_REFS.items()}
    with open(sub, "w") as fh:
        for acc, seq in refs.items():
            fh.write(f">{acc.split('.')[0]}\n{seq}\n")
    calls = {}
    accs = [r["accession"] for r in ha_records]
    for i in range(0, len(accs), 20000):
        seqs = _fasta_of(base, accs[i:i + 20000])
        q = Path(tmp_dir) / ".flu_ha_query.fasta"
        with open(q, "w") as fh:
            for a, sq in seqs.items():
                fh.write(f">{a}\n{sq}\n")
        res = subprocess.run([exe, "-query", str(q), "-subject", str(sub), "-task", "megablast", "-outfmt", "6 qseqid sseqid pident length bitscore"],
                             capture_output=True, text=True)
        best = {}
        for ln in res.stdout.splitlines():
            a, b, pid, ln_, bits = ln.split("\t")
            if float(pid) >= 75 and int(ln_) >= 300 and (a not in best or float(bits) > best[a][0]):
                best[a] = (float(bits), name_to_h[b])
        calls.update({a: v[1] for a, v in best.items()})
        print(f"  HA classified {min(i + 20000, len(accs))}/{len(accs)}", flush=True)
    return calls


def _flu_selection(base, taxa, hostmap, since, cap, cache_dir):
    """Selection of records for all segments, built once (locked, cached 12 h): reports of the taxa, subtype named from the strain's HA where the record has none,
    then at most `cap` strains per (subtype, collection year) -- chosen by a hash of the strain name, so every segment of a chosen strain is in."""
    key = hashlib.sha1(json.dumps(["sel", base, taxa, sorted(hostmap.items()), since, cap]).encode()).hexdigest()[:12]
    cache = Path(cache_dir) / f".ncbi_select_{key}.json.gz"
    lock = Path(cache_dir) / f".ncbi_select_{key}.lock"
    with open(lock, "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if cache.exists() and time.time() - cache.stat().st_mtime < 12 * 3600:
            return json.loads(gzip.open(cache, "rt", encoding="utf-8").read())
        reports = _cached_reports(base, taxa, hostmap, since, cache_dir)
        recs, seen = [], set()
        for r in reports:
            if r["accession"] in seen or not (r.get("segment") or "").strip() or r["_year"] < since:
                continue
            seen.add(r["accession"])
            r["_strain"] = ((r.get("isolate") or {}).get("name") or r["accession"]).strip()
            recs.append(r)
        unlabeled_ha = [r for r in recs if r["_subtype"] == "unclassified" and r["segment"] == "4"]
        calls = _classify_ha(base, unlabeled_ha, cache_dir)
        strain_sub = {}
        for r in unlabeled_ha:
            h = calls.get(r["accession"])
            if h:
                strain_sub.setdefault(r["_strain"], HA_SUBTYPE[h])
        for r in recs:
            if r["_subtype"] == "unclassified":
                r["_subtype"] = strain_sub.get(r["_strain"], "unclassified")
        recs = [r for r in recs if r["_subtype"] in KEEP_SUBTYPES]
        groups = {}
        for r in recs:
            groups.setdefault((r["_subtype"], r["_year"]), set()).add(r["_strain"])
        chosen = set()
        for strains in groups.values():
            chosen |= set(sorted(strains, key=lambda x: hashlib.sha1(x.encode()).hexdigest())[:cap])
        out = [r for r in recs if r["_strain"] in chosen]
        print(f"  influenza selection: {len(reports)} reports -> {len(recs)} typed records -> {len(out)} selected ({len(chosen)} strains)", flush=True)
        tmp = cache.with_suffix(".tmp")
        with gzip.open(tmp, "wt", encoding="utf-8") as fh:
            fh.write(json.dumps(out))
        os.replace(tmp, cache)
        return out


def main_multi(base, spec, out_fa, out_meta):
    base = base.rstrip("/")
    head, _, qs = spec.partition("?")
    taxa_s, _, seg = head.partition("/")
    taxa = [t for t in taxa_s.split(",") if t]
    par = dict(p.split("=", 1) for p in qs.split("&") if "=" in p)
    since = int(par.get("since") or 0)
    hostmap = dict(x.split(":", 1) for x in (par.get("host") or "").split(",") if ":" in x)
    sub, group = par.get("sub") or "", par.get("group") or "country"
    rows, seen = [], set()
    source = _flu_selection(base, taxa, hostmap, since, int(par.get("cap") or 800), Path(out_fa).parent) if par.get("classify") else _cached_reports(base, taxa, hostmap, since, Path(out_fa).parent)
    for r in source:
        if r["accession"] in seen or (seg and (r.get("segment") or "") != seg) or (since and r["_year"] < since):
            continue
        st = r["_subtype"]
        if sub and sub not in (re.match(r"(H\d+)", st) or [None, ""])[1:] + (re.search(r"(N\d+)$", st) or [None, ""])[1:]:
            continue
        seen.add(r["accession"])
        rows.append(r)
    seqs = {}
    for i in range(0, len(rows), 2000):
        z = zipfile.ZipFile(io.BytesIO(_post_zip(f"{base}/virus/genome/download", {"accessions": [r["accession"] for r in rows[i:i + 2000]], "include_sequence": ["GENOME"]})))
        name = None
        for ln in z.read("ncbi_dataset/data/genomic.fna").decode("utf-8").splitlines():
            if ln.startswith(">"):
                name = ln[1:].split()[0].split(".")[0]
                seqs[name] = []
            elif name:
                seqs[name].append(ln.strip())
    cols = ["accession", "strain", "date", "date_released", "date_updated", "region", "country", "division", "location", "length", "host", "taxid",
            "is_lab_host", "completeness", "authors", "institution", "submitter_country", "sample_type", "country_lineage", "url", "subtype"]
    n = 0
    with open(out_fa, "w", encoding="utf-8") as fa, open(out_meta, "w", encoding="utf-8", newline="") as fm:
        w = csv.writer(fm, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            acc = r["accession"].split(".")[0]
            if acc not in seqs:
                continue
            fa.write(f">{acc}\n" + "\n".join(seqs[acc]) + "\n")
            n += 1
            iso, host, sub_ = r.get("isolate") or {}, r.get("host") or {}, r.get("submitter") or {}
            geo = ((r.get("location") or {}).get("geographic_location") or "").strip()
            country, _, rest = geo.partition(":")
            country = country.strip()
            if country.lower() in ("missing", "not collected", "not available", "na", "unknown"):
                country = ""
            date = iso.get("collection_date") or ""
            w.writerow([acc, iso.get("name") or "", date, (r.get("release_date") or "")[:10], (r.get("update_date") or "")[:10], "Africa" if country in AFRICA else "",
                        country, rest.strip(), "", r.get("length") or "", host.get("organism_name") or "", host.get("tax_id") or "", "false", r.get("completeness") or "",
                        "; ".join(sub_.get("names") or []), sub_.get("affiliation") or "", sub_.get("country") or "", "",
                        (r["_subtype"] if group == "subtype" else (country or "unknown")), f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}", r["_subtype"]])
    print(f"taxa {','.join(taxa)} segment {seg or '-'}: {len(rows)} records selected, {n} sequences written")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    if "," in sys.argv[2] or "?" in sys.argv[2]:
        main_multi(*sys.argv[1:])
    else:
        main(*sys.argv[1:])
