#!/usr/bin/env python3
"""Step 3: mfeprimer F/R primer-pair search (--virus mode) for one assay set at one
mismatch level. Writes the mfeprimer report + a parsed one-line coverage summary."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import (  # noqa: E402
    find_mfeprimer, ensure_index, write_oligo_fasta, run_primer_search,
    parse_primer_coverage, load_primers_manifest,
)


def main():
    set_id = snakemake.wildcards.set_id
    mismatch = int(snakemake.wildcards.mismatch)
    db_fasta = Path(snakemake.input.db)
    manifest = load_primers_manifest(Path(snakemake.input.manifest))
    row = next(r for r in manifest if r["set_id"] == set_id)

    pipeline_root = Path(snakemake.scriptdir).parent
    mfeprimer_bin = find_mfeprimer(snakemake.config.get("mfeprimer_bin", ""), pipeline_root)
    cpu = int(snakemake.config.get("mfeprimer_cpu", 0) or 0) or None
    ensure_index(mfeprimer_bin, db_fasta, cpu=cpu)

    fr_fa = Path(snakemake.output.report).with_name(f"{set_id}_FR.fa")
    write_oligo_fasta(fr_fa, [(row["F_id"], row["F_seq"]), (row["R_id"], row["R_seq"])])

    out_txt = Path(snakemake.output.report)
    run_primer_search(mfeprimer_bin, db_fasta, fr_fa, out_txt, mismatch, cpu=cpu)
    stats = parse_primer_coverage(out_txt)

    with open(snakemake.output.coverage, "w", encoding="utf-8") as fh:
        fh.write("set_id\tmismatch\tcoverage_pct\ttotal\tmissed\n")
        fh.write(f"{set_id}\t{mismatch}\t{stats['coverage_pct']}\t{stats['total']}\t{stats['missed']}\n")

    # mfeprimer skips writing .amp.fasta/.spec.tsv/.virus.Failed.txt entirely when it
    # finds zero amplicons (e.g. a primer set with ~0% coverage at mismatch=0) -- touch
    # empty placeholders so Snakemake's declared outputs always exist and downstream
    # rules (which already treat "no rows" as "no matches") don't need special-casing.
    # (Only rules that actually declare these extra outputs, e.g. the population
    # primer_match rule, need this; standard_primer_match doesn't declare them.)
    for out_name in ("amp", "spec", "failed"):
        out_path = getattr(snakemake.output, out_name, None)
        if out_path is None:
            continue
        p = Path(out_path)
        if not p.exists():
            p.touch()


if __name__ == "__main__":
    main()
