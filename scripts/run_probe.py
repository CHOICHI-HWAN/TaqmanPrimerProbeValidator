#!/usr/bin/env python3
"""Step 4: mfeprimer single-oligo `search` for one assay set's probe, restricted to
perfect matches (mismatch=0 by default via config['probe_mismatch']). A TaqMan probe is
far more mismatch-sensitive than its flanking primers, so this is checked separately
from -- and at a stricter tolerance than -- the primer-pair search in step 3."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import (  # noqa: E402
    find_mfeprimer, ensure_index, write_oligo_fasta, run_probe_search, load_primers_manifest, write_ids,
)


def main():
    set_id = snakemake.wildcards.set_id
    mismatch = int(snakemake.params.probe_mismatch)
    db_fasta = Path(snakemake.input.db)
    manifest = load_primers_manifest(Path(snakemake.input.manifest))
    row = next(r for r in manifest if r["set_id"] == set_id)

    pipeline_root = Path(snakemake.scriptdir).parent
    mfeprimer_bin = find_mfeprimer(snakemake.config.get("mfeprimer_bin", ""), pipeline_root)
    ensure_index(mfeprimer_bin, db_fasta)

    probe_fa = Path(snakemake.output.report).with_name(f"{set_id}_probe.fa")
    write_oligo_fasta(probe_fa, [(row["P_id"], row["P_seq"])])

    out_txt = Path(snakemake.output.report)
    matched_ids = run_probe_search(mfeprimer_bin, db_fasta, probe_fa, out_txt, mismatch)
    write_ids(Path(snakemake.output.matched_ids), matched_ids)


if __name__ == "__main__":
    main()
