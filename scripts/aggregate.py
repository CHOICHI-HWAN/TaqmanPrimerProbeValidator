#!/usr/bin/env python3
"""Step 5 (final): concatenate every set's per-set summary/lineage rows into the two
pipeline-wide result tables."""
import csv
from pathlib import Path


def concat(paths, out_path):
    rows = []
    header = None
    for p in paths:
        with open(p, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            h = next(reader)
            header = header or h
            rows.extend(list(reader))
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    concat(snakemake.input.summaries, Path(snakemake.output.summary))
    concat(snakemake.input.lineages, Path(snakemake.output.summary_by_lineage))


if __name__ == "__main__":
    main()
