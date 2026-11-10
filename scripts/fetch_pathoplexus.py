#!/usr/bin/env python3
"""Fetch one organism's genomes from Pathoplexus (LAPIS, open data) and write them in the shape the pipeline's Nextstrain fetch produces:
a FASTA with one record per genome (header = accession) and a metadata TSV with the columns the virus profile maps (accession, date,
country, region, host, is_lab_host, ...). Only open-access, latest, non-revoked versions are requested (what the API serves without a login).

usage: fetch_pathoplexus.py <lapis_url> <organism> <out.fasta> <out_metadata.tsv>      e.g. https://lapis.pathoplexus.org ebola-sudan ...
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


def main(base, organism, out_fa, out_meta):
    q = urllib.parse.urlencode(BASE_FILTER)
    root = f"{base.rstrip('/')}/{organism}/sample"
    det = _get(f"{root}/details?{q}&dataFormat=tsv&fields={','.join(FIELDS)}").decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(det), delimiter="\t"))
    ids, seen = {}, set()
    for r in rows:
        a = (r.get("insdcAccessionBase") or "").strip() or r["accession"]
        if a in seen:
            a = r["accession"]
        seen.add(a)
        ids[r["accessionVersion"]] = a
    fa = _get(f"{root}/unalignedNucleotideSequences?{q}").decode("utf-8")
    n = 0
    with open(out_fa, "w", encoding="utf-8") as fh:
        keep, name = False, None
        for ln in fa.splitlines():
            if ln.startswith(">"):
                name = ids.get(ln[1:].strip().split()[0])
                keep = name is not None
                if keep:
                    fh.write(f">{name}\n")
                    n += 1
            elif keep:
                fh.write(ln + "\n")
    cols = ["accession", "strain", "date", "date_released", "date_updated", "region", "country", "division", "location", "length", "host", "taxid",
            "is_lab_host", "completeness", "authors", "institution", "submitter_country", "sample_type", "country_lineage", "url"]
    with open(out_meta, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            a = ids[r["accessionVersion"]]
            c = (r.get("geoLocCountry") or "").strip()
            lab = (r.get("isLabHost") or "").strip().lower()
            w.writerow([a, r.get("displayName") or "", r.get("sampleCollectionDate") or "", r.get("releasedDate") or "", r.get("ncbiUpdateDate") or "",
                        "Africa" if c in AFRICA else "", c, r.get("geoLocAdmin1") or "", r.get("geoLocAdmin2") or "", r.get("length") or "",
                        (r.get("hostNameScientific") or r.get("hostNameCommon") or "").strip(), r.get("hostTaxonId") or "",
                        "true" if lab in ("true", "1", "yes") else "false", r.get("completeness") or "", r.get("authors") or "",
                        r.get("authorAffiliations") or "", r.get("ncbiSubmitterCountry") or "", r.get("sampleType") or "", c or "unknown",
                        f"https://pathoplexus.org/seq/{r['accessionVersion']}"])
    print(f"{organism}: {len(rows)} genomes in the table, {n} sequences written")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
