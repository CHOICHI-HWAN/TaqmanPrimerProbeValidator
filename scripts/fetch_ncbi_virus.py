#!/usr/bin/env python3
"""Fetch every genome of one virus taxon from the NCBI Datasets API (public, no key) and write them in the shape the pipeline's other fetch steps
produce: a FASTA (header = accession without version) and a metadata TSV (accession, date, country, region, host, is_lab_host, ...).

usage: fetch_ncbi_virus.py <datasets_api_url> <taxon_id> <out.fasta> <out_metadata.tsv>      e.g. https://api.ncbi.nlm.nih.gov/datasets/v2 186539 ...
The extra column `country_lineage` (= country) lets a virus without a clade scheme group its genomes by the country they were collected in.
"""
import csv
import io
import json
import sys
import urllib.request
import zipfile

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


def main(base, taxon, out_fa, out_meta):
    base = base.rstrip("/")
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


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
