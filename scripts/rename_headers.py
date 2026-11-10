#!/usr/bin/env python3
"""Step 1b: build a human-readable, metadata-annotated copy of the raw Nextstrain fasta.

Header format: "LINEAGE | ACCESSION | DATE | COUNTRY | REGION" (matches the convention
used in the previous manual analysis, kept here purely for human review / dashboard use).

Also emits a per-accession metadata sidecar TSV (accession, lineage, date, country, region)
that the filter step joins against -- mfeprimer itself is run later against a *separate*
bare-accession fasta (see 02_filter.py) because a "LINEAGE | ACC | ..." header would make
mfeprimer key every hit on the (non-unique) lineage prefix instead of the accession.
"""
import csv
from pathlib import Path


def read_fasta(path: Path):
    header, chunks = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def load_metadata(meta_path: Path, lineage_col: str, columns: dict, lineage_override: dict = None, lineage_missing: str = "NA") -> dict:
    """columns: the profile's mapping from the source metadata's column names to the canonical fields."""
    mapping = {}
    with open(meta_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row.get(columns["accession"], "").strip()
            if not acc:
                continue
            if lineage_override is not None:
                lineage = lineage_override.get(acc) or lineage_missing      # clade from the Auspice build
            else:
                lineage = (row.get(lineage_col) or row.get(columns["lineage_fallback"]) or "NA").strip() or "NA"
            date = (row.get(columns["date"]) or "XXXX-XX-XX").strip() or "XXXX-XX-XX"
            country = (row.get(columns["country"]) or "").strip() or "Unknown/Missing"
            region = (row.get(columns["region"]) or "").strip() or "Unknown/Missing"
            mapping[acc] = {"lineage": lineage, "date": date, "country": country, "region": region}
    return mapping


def main():
    raw_fasta = Path(snakemake.input.fasta)
    raw_meta = Path(snakemake.input.meta)
    lineage_col = snakemake.params.lineage_col
    out_renamed = Path(snakemake.output.renamed_fasta)
    out_meta = Path(snakemake.output.metadata_curated)

    override = None
    table = snakemake.input.get("lineage_table") if hasattr(snakemake.input, "get") else None
    if table:
        with open(table, newline="", encoding="utf-8") as fh:
            override = {r["accession"]: r["lineage"] for r in csv.DictReader(fh, delimiter="\t")}
    mapping = load_metadata(raw_meta, lineage_col, dict(snakemake.params.columns), override, snakemake.params.lineage_missing)

    n_written, n_no_meta = 0, 0
    with open(out_renamed, "w", encoding="utf-8") as out_fa, \
         open(out_meta, "w", newline="", encoding="utf-8") as out_tsv:
        writer = csv.writer(out_tsv, delimiter="\t")
        writer.writerow(["accession", "lineage", "date", "country", "region"])
        for acc, seq in read_fasta(raw_fasta):
            acc_clean = acc.split()[0].split("|")[0].strip()
            m = mapping.get(acc_clean)
            if m is None:
                n_no_meta += 1
                continue
            header = f"{m['lineage']} | {acc_clean} | {m['date']} | {m['country']} | {m['region']}"
            out_fa.write(f">{header}\n{seq.upper()}\n")
            writer.writerow([acc_clean, m["lineage"], m["date"], m["country"], m["region"]])
            n_written += 1

    with open(snakemake.log[0], "w", encoding="utf-8") as logf:
        logf.write(f"segment={snakemake.wildcards.seg} written={n_written} no_metadata={n_no_meta}\n")


if __name__ == "__main__":
    main()
