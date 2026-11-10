#!/usr/bin/env python3
"""Step 5: reconcile primer-pair matches with perfect-match probe matches for one set.

validated  = (genomes where F/R matched at config['primer_mismatch_pass'])
             AND (genomes where the probe matched at config['probe_mismatch'], default 0)
probe_dropout = primer-matched genomes that the probe fails to bind perfectly --
                excluded from "validated" because a mismatched TaqMan probe can fail to
                generate signal even when the flanking primers amplify fine.

Derives everything from the already-computed per-genome mismatch_table.tsv (built by
build_alignment_table.py at mismatch<=2 for F/R and <=8 for the probe, with the exact
per-genome mismatch count independently recomputed from the actual matched sequence)
rather than running a second, separate mfeprimer search directly at
--misMatch=primer_mismatch_pass. Those two used to disagree: mfeprimer's own seed-finding
step happens BEFORE its mismatch filtering, so a direct search run at a low --misMatch can
fail to find a site at all (report it as a hard miss) in cases the wider mismatch_table.tsv
search -- and an independent Python-side recount of the real mismatch value -- correctly
resolves as a genuine 0/1-mismatch match. Confirmed concretely for
VendorD_GPC_PREVIOUS1_LASV1_P2/AY628201: a real 0-mismatch F/P/R match that a direct
--misMatch 1 search nonetheless reported as failed. Deriving "validated" from the same
mismatch_table.tsv the per-genome detail table and MSA already use makes reconcile.py's
counts agree with what the dashboard's own detail views show, instead of silently
undercounting."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from host_scope import ids_in_scope, load_host_classes, primary_scope, summarize  # noqa: E402
from mfe_common import probe_inside_amplicon, read_fasta_ids, write_ids  # noqa: E402


def main():
    set_id = snakemake.wildcards.set_id
    segment = snakemake.params.segment
    source = snakemake.params.source
    primer_mismatch_pass = int(snakemake.params.primer_mismatch_pass)
    probe_mismatch = int(snakemake.params.probe_mismatch)

    db_fasta = Path(snakemake.input.db)
    all_ids = sorted(set(read_fasta_ids(db_fasta)))

    rows_by_id = {}
    with open(snakemake.input.mismatch_table, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows_by_id[r["genome_id"]] = r

    # The id lists cover EVERY genome of the population (the tree overlay and the group unions need them all); the pass/fail rule
    # itself is host_scope.genome_flags, shared with the web app.
    full = summarize(all_ids, rows_by_id, primer_mismatch_pass, probe_mismatch, probe_inside_amplicon)
    write_ids(Path(snakemake.output.validated), full["_validated_ids"])
    write_ids(Path(snakemake.output.probe_dropout), full["_dropout_ids"])
    write_ids(Path(snakemake.output.primer_failed), full["_failed_ids"])

    # The summary / lineage tables are for the PRIMARY host scope (human genomes; "all" when the metadata has no human genome).
    # Laboratory hosts are in no scope. The web app can show the other scopes (animal / unknown / all) with the same code.
    # (standard_reconcile reuses this script on the standard-material sequences, which have no host: no host input -> everything counts)
    if "host_meta" in snakemake.input.keys():
        classes = load_host_classes(snakemake.input.host_meta, snakemake.params.host_columns)
    else:
        classes = {}
    scope = primary_scope(all_ids, classes)
    sub = summarize(ids_in_scope(all_ids, classes, scope), rows_by_id, primer_mismatch_pass, probe_mismatch, probe_inside_amplicon)

    with open(snakemake.output.summary_row, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "set_id", "segment", "source", "total_genomes",
            "primer_coverage_mm0_pct", "primer_coverage_mm1_pct", "primer_coverage_mm2_pct",
            "probe_coverage_mm0_pct", "primer_mismatch_used",
            "primer_matched", "probe_dropout", "validated", "validated_pct", "host_scope",
        ])
        writer.writerow([
            set_id, segment, source, sub["total_genomes"],
            sub["primer_coverage_mm0_pct"], sub["primer_coverage_mm1_pct"], sub["primer_coverage_mm2_pct"],
            sub["probe_coverage_mm0_pct"], primer_mismatch_pass,
            sub["primer_matched"], sub["probe_dropout"], sub["validated"], sub["validated_pct"], scope,
        ])

    # lineage breakdown of the validated set (primary scope)
    meta = {}
    with open(snakemake.input.metadata, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["accession"]] = row
    lineage_counts = {}
    for acc in sub["_validated_ids"]:
        lineage = meta.get(acc, {}).get("lineage", "NA")
        lineage_counts[lineage] = lineage_counts.get(lineage, 0) + 1
    with open(snakemake.output.lineage_rows, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["set_id", "segment", "lineage", "validated_count"])
        for lineage, cnt in sorted(lineage_counts.items(), key=lambda x: -x[1]):
            writer.writerow([set_id, segment, lineage, cnt])


if __name__ == "__main__":
    main()
