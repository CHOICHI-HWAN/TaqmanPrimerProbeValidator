#!/usr/bin/env python3
"""Fetch one organism's genomes from Pathoplexus (LAPIS, open data) and write them in the shape the pipeline's Nextstrain fetch produces:
a FASTA with one record per genome (header = accession) and a metadata TSV with the columns the virus profile maps (accession, date,
country, region, host, is_lab_host, ...). Only open-access, latest, non-revoked versions are requested (what the API serves without a login).

usage: fetch_pathoplexus.py <lapis_url> <organism>[@<lineage column>][/<SEGMENT>] <out.fasta> <out_metadata.tsv>      e.g. https://lapis.pathoplexus.org ebola-sudan ...
A segmented organism (cchf) is fetched one segment at a time: `cchf/S` takes the S sequences and the per-segment columns (insdcAccessionBase_S, length_S,
lineage_S ...); the lineage_S column (S-segment genotype) becomes the `country_lineage` column of every segment of that record is not possible
(reassortment), so for segmented organisms `country_lineage` holds the genotype of the S segment when the record has one, else "unclassified".
The accession is the INSDC one (e.g. NC_006432, KY...) when the record has it, otherwise the Pathoplexus accession (PP_...).
The extra column `country_lineage` (= country) lets a virus without a clade scheme group its genomes by outbreak country.
"""
import csv
import io
import sys
import urllib.parse
import urllib.request

FIELDS = ["accession", "accessionVersion", "insdcAccessionBase", "insdcAccessionFull", "sampleCollectionDate", "geoLocCountry", "geoLocAdmin1",
          "geoLocAdmin2", "hostNameScientific", "hostNameCommon", "hostTaxonId", "isLabHost", "length", "authors", "authorAffiliations",
          "ncbiSubmitterCountry", "releasedDate", "ncbiUpdateDate", "sampleType", "completeness", "displayName"]
PAGE, PAGE_ABOVE = 1000, 2500          # sequences per request / row count from which the sequences are fetched in pages
BASE_FILTER = {"versionStatus": "LATEST_VERSION", "isRevocation": "false"}
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


SEG_FIELDS = ("insdcAccessionBase", "insdcAccessionFull", "length", "completeness", "ncbiUpdateDate")


def main(base, organism, out_fa, out_meta):
    organism, _, seg = organism.partition("/")
    organism, _, lf = organism.partition("@")                  # organism@clade: that metadata column becomes the lineage column (mpox)
    sfx = f"_{seg}" if seg else ""
    q = urllib.parse.urlencode(BASE_FILTER)
    root = f"{base.rstrip('/')}/{organism}/sample"
    fields = [f + sfx if f in SEG_FIELDS else f for f in FIELDS] + (["lineage_S"] if seg else []) + ([lf] if lf else [])
    det = _get(f"{root}/details?{q}&dataFormat=tsv&fields={','.join(fields)}").decode("utf-8")
    rows = [dict(r) for r in csv.DictReader(io.StringIO(det), delimiter="\t")]
    for r in rows:                                             # one canonical column set whatever the organism's segment layout
        for f in SEG_FIELDS:
            r[f] = r.get(f + sfx, "")
    if seg:
        rows = [r for r in rows if (r.get("insdcAccessionBase") or "").strip() or int(float(r.get("length") or 0)) > 0]
        rows = [r for r in rows if int(float(r.get("length") or 0)) > 0]
    ids, seen = {}, set()
    for r in rows:
        a = (r.get("insdcAccessionBase") or "").strip() or r["accession"]
        if a in seen:
            a = r["accession"]
        seen.add(a)
        ids[r["accession"]] = a
        r["_id"] = a
    seq_url = f"{root}/unalignedNucleotideSequences{'/' + seg if seg else ''}?{q}"
    n = 0
    with open(out_fa, "w", encoding="utf-8") as fh:
        # a big organism (mpox: ~15 000 genomes of ~197 kb) is fetched page by page in a fixed order, so no single answer is gigabytes
        # (one request for everything stalled and timed out); a small one stays a single request
        page = PAGE if len(rows) > PAGE_ABOVE else None
        offset = 0
        while True:
            fa = _get(seq_url + (f"&orderBy=accessionVersion&limit={page}&offset={offset}" if page else "")).decode("utf-8")
            got = 0
            keep = False
            for ln in fa.splitlines():
                if ln.startswith(">"):
                    got += 1
                    name = ids.get(ln[1:].strip().split()[0].split("|")[0].rsplit(".", 1)[0])
                    keep = name is not None
                    if keep:
                        fh.write(f">{name}\n")
                        n += 1
                elif keep:
                    fh.write(ln + "\n")
            if not page or got < page:
                break
            offset += page
            print(f"  {organism}: {offset} sequences read ...", flush=True)
    cols = ["accession", "strain", "date", "date_released", "date_updated", "region", "country", "division", "location", "length", "host", "taxid",
            "is_lab_host", "completeness", "authors", "institution", "submitter_country", "sample_type", "country_lineage", "url"]
    with open(out_meta, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            a = r["_id"]
            c = (r.get("geoLocCountry") or "").strip()
            lab = (r.get("isLabHost") or "").strip().lower()
            w.writerow([a, r.get("displayName") or "", r.get("sampleCollectionDate") or "", r.get("releasedDate") or "", r.get("ncbiUpdateDate") or "",
                        "Africa" if c in AFRICA else "", c, r.get("geoLocAdmin1") or "", r.get("geoLocAdmin2") or "", r.get("length") or "",
                        (r.get("hostNameScientific") or r.get("hostNameCommon") or "").strip(), r.get("hostTaxonId") or "",
                        "true" if lab in ("true", "1", "yes") else "false", r.get("completeness") or "", r.get("authors") or "",
                        r.get("authorAffiliations") or "", r.get("ncbiSubmitterCountry") or "", r.get("sampleType") or "", ((r.get(lf) or "unclassified") if lf else (r.get("lineage_S") or "unclassified") if seg else (c or "unknown")),
                        f"https://pathoplexus.org/seq/{r['accessionVersion']}"])
    print(f"{organism}: {len(rows)} genomes in the table, {n} sequences written")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
