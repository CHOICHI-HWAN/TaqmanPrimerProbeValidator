#!/usr/bin/env python3
"""Build a genuine multiple-sequence alignment of every matched genome's amplicon for one
assay set, in the style of the reference Nipah henipavirus screenshot: all genomes shown
as rows at once (not one at a time), a consensus row, and the Forward/Probe/Reverse
regions marked as fixed column bands across every row.

Uses mafft (fast, <1s for ~700 short 100-300bp sequences) on the per-genome amplicons
already extracted by the primer_match step. The actual alignment logic lives in
scripts/msa_common.py so the interactive webapp's live analysis (live_analyze.py) can
build the exact same view for a user-submitted set.
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import load_primers_manifest  # noqa: E402
from msa_common import build_msa  # noqa: E402


def main():
    set_id = snakemake.wildcards.set_id
    manifest = load_primers_manifest(Path(snakemake.input.manifest))
    row = next(r for r in manifest if r["set_id"] == set_id)

    genome_rows = []
    with open(snakemake.input.mismatch_table, newline="", encoding="utf-8") as fh:
        for g in csv.DictReader(fh, delimiter="\t"):
            genome_rows.append(g)

    out_dir = Path(snakemake.output.msa).parent
    result = build_msa(genome_rows, row["F_seq"], row["P_seq"], row["R_seq"], out_dir)
    Path(snakemake.output.msa).write_text(json.dumps(result), encoding="utf-8")


if __name__ == "__main__":
    main()
