#!/usr/bin/env python3
"""Concatenate every combined-panel group's summary row into one table."""
import csv

GROUP_SUMMARY_HEADER = ["group_id", "label", "segment", "members", "total_genomes", "validated", "validated_pct"]


def main():
    rows = []
    header = None   # no groups yet -> still write the (empty) table with its header
    for p in snakemake.input.summaries:
        with open(p, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            h = next(reader)
            header = header or h
            rows.extend(list(reader))
    header = header or GROUP_SUMMARY_HEADER
    with open(snakemake.output.summary, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
