#!/usr/bin/env python3
"""Union the validated-genome sets of several assay sets that would be used together as a
panel (e.g. "use S1_P1 OR S1_P2 OR S2_P1 OR S2_P2"), matching the '<<ALL' / '<<P1+P2'
combined columns in the reference report."""
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from host_scope import ids_in_scope, load_host_classes, primary_scope  # noqa: E402
from mfe_common import read_fasta_ids, read_ids, write_ids  # noqa: E402


def main():
    group_id = snakemake.wildcards.group_id
    label = snakemake.params.label
    segment = snakemake.params.segment
    members = snakemake.params.members

    all_ids = set(read_fasta_ids(Path(snakemake.input.db)))
    union_ids = set()
    for member_path in snakemake.input.members:
        union_ids |= read_ids(Path(member_path))

    write_ids(Path(snakemake.output.validated), union_ids)

    # The validated_ids.txt above stays the union over EVERY genome (the tree overlay and the cross-binding step need them all);
    # the summary row counts the PRIMARY host scope (human genomes; "all" without host data), like the member sets' summary rows.
    classes = load_host_classes(snakemake.input.host_meta, snakemake.params.host_columns)
    scope = primary_scope(all_ids, classes)
    scoped = set(ids_in_scope(all_ids, classes, scope))
    total = len(scoped)
    validated = len(union_ids & scoped)
    with open(snakemake.output.summary_row, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["group_id", "label", "segment", "members", "total_genomes", "validated", "validated_pct", "host_scope"])
        writer.writerow([
            group_id, label, segment, ";".join(members), total, validated,
            round(100 * validated / total, 2) if total else 0.0, scope,
        ])

if __name__ == "__main__":
    main()
