#!/usr/bin/env python3
"""config/refgenome_genes.tsv (segment, accession, genome_len, gene, start, end, strand) from a GFF3
gene map: python3 scripts/genes_from_gff.py <genemap.gff> <segment key> <out.tsv>.
Only CDS features with a `gene=` attribute are used (1-based inclusive coordinates, as in the GFF)."""
import csv
import re
import sys


def main(gff, seg, out):
    genome_len, acc, rows = None, None, []
    for line in open(gff, encoding="utf-8"):
        if line.startswith("##sequence-region"):
            _, acc, _s, genome_len = line.split()[:4]
            acc = acc.split(".")[0]      # the reference is extracted from the raw data under its unversioned accession
        elif not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[2] == "CDS":
                m = re.search(r"gene=([^;]+)", f[8])
                if m:
                    rows.append((seg, acc, genome_len, m.group(1), f[3], f[4], f[6]))
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["segment", "accession", "genome_len", "gene", "start", "end", "strand"])
        w.writerows(rows)
    print(f"{len(rows)} genes -> {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
