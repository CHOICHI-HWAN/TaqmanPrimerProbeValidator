#!/usr/bin/env python3
"""Same pooled cross-set combo search as group_cross_binding.py, but against the tiny
KDCA standard-material db and the single reference genome instead of the population db --
so a group's "표준물질 매치" and "참조유전체 결합위치" views show the same cross-set-aware
picture (e.g. Set2's probe landing inside Set1's amplicon) that the population-level
cross-binding analysis already does, instead of only ever checking each member set against
itself independently."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import find_mfeprimer, load_primers_manifest  # noqa: E402
from group_cross_binding import compute_group_combo_table  # noqa: E402


def main():
    manifest = load_primers_manifest(snakemake.input.manifest)
    member_ids = snakemake.params.members
    rows_by_set = {r["set_id"]: r for r in manifest if r["set_id"] in member_ids}

    out_dir = Path(snakemake.output.table).parent
    pipeline_root = Path(snakemake.scriptdir).parent
    mfeprimer_bin = find_mfeprimer(snakemake.config.get("mfeprimer_bin", ""), pipeline_root)

    metadata_path = Path(snakemake.input.metadata) if hasattr(snakemake.input, "metadata") else None

    compute_group_combo_table(
        db_fasta=Path(snakemake.input.db),
        metadata_path=metadata_path,
        rows_by_set=rows_by_set,
        mfeprimer_bin=mfeprimer_bin,
        out_dir=out_dir,
        output_table=Path(snakemake.output.table),
        gene_summary_path=None,
        apply_size_filter=False,
        tag=snakemake.params.tag,
    )


if __name__ == "__main__":
    main()
