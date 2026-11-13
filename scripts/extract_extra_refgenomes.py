#!/usr/bin/env python3
"""Writes the extra reference genomes of a unit (profile segments.<unit>.extra_references, a list of accessions) into one multi-record FASTA.
The file is empty when the unit has none."""
from pathlib import Path

from mfe_common import read_fasta_seqs


def main():
    wanted = list(snakemake.params.accessions)
    out = Path(snakemake.output.ref)
    seqs = {k.split()[0].split("|")[0].strip(): v for k, v in read_fasta_seqs(Path(snakemake.input.fasta)).items()}
    missing = [a for a in wanted if a not in seqs]
    if missing:
        raise SystemExit(f"[error] extra reference accession(s) not found in {snakemake.input.fasta}: {', '.join(missing)}")
    with open(out, "w", encoding="utf-8") as fh:
        for a in wanted:
            fh.write(f">{a}\n{seqs[a].upper()}\n")


if __name__ == "__main__":
    main()
