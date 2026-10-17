#!/usr/bin/env python3
"""Summarize, per set: which gene it targets and where the amplicon falls on the
reference genome (start..end, size) -- answers "전체 게놈부위에서 어느 부위에 바인딩하는지".
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import load_primers_manifest  # noqa: E402


def load_genes(path):
    genes = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genes.setdefault(row["segment"], []).append(row)
    return genes


def genes_overlapping(genes_for_seg, start, end):
    hits = []
    for g in genes_for_seg:
        gs, ge = int(g["start"]), int(g["end"])
        if start <= ge and end >= gs:
            hits.append(g["gene"])
    return hits


def main():
    genes = load_genes(snakemake.input.genes)
    manifest = {r["set_id"]: r for r in load_primers_manifest(snakemake.input.manifest)}

    with open(snakemake.output.summary, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["set_id", "segment", "target_gene", "amplicon_start", "amplicon_end",
                          "amplicon_size", "chrom"])
        for binding_path in snakemake.input.bindings:
            set_id = Path(binding_path).parent.parent.name  # output/primers/<id>/refgenome_binding/binding.tsv
            seg = manifest[set_id]["segment"].lower()
            starts, ends, chrom = [], [], ""
            with open(binding_path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if row["found"] == "1":
                        starts.append(int(row["start"]))
                        ends.append(int(row["end"]))
                        chrom = row["chrom"]
            if not starts:
                writer.writerow([set_id, seg.upper(), "NA", "", "", "", ""])
                continue
            amp_start, amp_end = min(starts), max(ends)
            gene_hits = genes_overlapping(genes.get(seg, []), amp_start, amp_end)
            target_gene = "/".join(dict.fromkeys(gene_hits)) if gene_hits else "NA"
            writer.writerow([set_id, seg.upper(), target_gene, amp_start, amp_end,
                              amp_end - amp_start, chrom])


if __name__ == "__main__":
    main()
