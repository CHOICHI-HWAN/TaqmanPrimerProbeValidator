#!/usr/bin/env python3
"""Pull one named accession's full sequence out of the raw Nextstrain fasta, to use as the
single canonical reference genome for the per-oligo binding-site maps (mirrors the
'참조유전체비교' AY628203/AY628202 slides)."""
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


def main():
    fasta = Path(snakemake.input.fasta)
    accession = snakemake.params.accession
    out = Path(snakemake.output.ref)

    for header, seq in read_fasta(fasta):
        acc_clean = header.split()[0].split("|")[0].strip()
        if acc_clean == accession:
            with open(out, "w", encoding="utf-8") as fh:
                fh.write(f">{accession}\n{seq.upper()}\n")
            return
    raise SystemExit(f"[error] reference accession {accession!r} not found in {fasta}")


if __name__ == "__main__":
    main()
