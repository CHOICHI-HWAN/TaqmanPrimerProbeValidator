#!/usr/bin/env python3
"""Splits the one downloaded population (data/raw/_all_*) into the analysis units whose profile segment lists `lineage_in`: a record goes to the unit whose
list contains the value of its lineage column. A mixed call such as "Ib/IIb" is read as its first part ("Ib"). Records with a value no unit lists are left out.

Used when one organism needs several reference genomes (the units are separate analysis segments, each with its own reference_accession).
"""
import csv
from pathlib import Path


def first_part(v):
    return (v or "").split("/")[0].strip()


def main():
    allowed = set(snakemake.params.values)
    col = snakemake.params.col
    keep = set()
    with open(snakemake.input.meta, newline="", encoding="utf-8") as fin, open(snakemake.output.meta, "w", newline="", encoding="utf-8") as fout:
        r = csv.DictReader(fin, delimiter="\t")
        w = csv.DictWriter(fout, fieldnames=r.fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in r:
            v = first_part(row.get(col))
            if v in allowed:
                row[col] = v
                keep.add(row["accession"])
                w.writerow(row)
    n = 0
    with open(snakemake.input.fasta, encoding="utf-8") as fin, open(snakemake.output.fasta, "w", encoding="utf-8") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                write = line[1:].split()[0] in keep
                n += write
            if write:
                fout.write(line)
    print(f"{Path(snakemake.output.fasta).name}: {n} sequences of {len(keep)} records ({', '.join(sorted(allowed))})")


main()
