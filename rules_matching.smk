# Primer/probe matching against the population and the final reconciliation (validated = primer AND probe).
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 3. primer match ---------------------------------------------------------
rule primer_match:
    input:
        db=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered.fasta",
        db_fai=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered.fasta.fai",
        manifest=primer_input,
    output:
        report="output/primers/{set_id}/primer_qc/{set_id}_mm{mismatch}.txt",
        coverage="output/primers/{set_id}/primer_qc/{set_id}_mm{mismatch}_coverage.tsv",
        amp="output/primers/{set_id}/primer_qc/{set_id}_mm{mismatch}.txt.amp.fasta",
        spec="output/primers/{set_id}/primer_qc/{set_id}_mm{mismatch}.txt.spec.tsv",
        failed="output/primers/{set_id}/primer_qc/{set_id}_mm{mismatch}.txt.virus.Failed.txt",
    threads: config["mfeprimer_cpu"]
    log:
        "logs/primer_match_{set_id}_mm{mismatch}.log",
    script:
        "scripts/run_primer.py"


# ---- 4. probe match (perfect match only) -------------------------------------
rule probe_match:
    input:
        db=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered.fasta",
        db_fai=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered.fasta.fai",
        manifest=primer_input,
    output:
        report="output/primers/{set_id}/probe_qc/{set_id}_probe_mm0.txt",
        matched_ids="output/primers/{set_id}/probe_qc/{set_id}_probe_matched_ids.txt",
    params:
        probe_mismatch=config["probe_mismatch"],
    log:
        "logs/probe_match_{set_id}.log",
    script:
        "scripts/run_probe.py"


# ---- 5. reconcile + aggregate -------------------------------------------------
rule reconcile:
    input:
        db=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered.fasta",
        metadata=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered_metadata.tsv",
        mismatch_table="output/primers/{set_id}/alignment/mismatch_table.tsv",
        host_meta=lambda wc: f"data/raw/{SET_SEG[wc.set_id]}_metadata.tsv",     # host / lab-host columns -> host scope
    output:
        validated="output/primers/{set_id}/final/validated_ids.txt",
        probe_dropout="output/primers/{set_id}/final/probe_dropout_ids.txt",
        primer_failed="output/primers/{set_id}/final/primer_failed_ids.txt",
        summary_row="output/primers/{set_id}/final/summary_row.tsv",
        lineage_rows="output/primers/{set_id}/final/lineage_rows.tsv",
    params:
        segment=lambda wc: SET_SEG[wc.set_id].upper(),
        source=lambda wc: SET_SOURCE[wc.set_id],
        primer_mismatch_pass=PRIMER_MISMATCH_PASS,
        probe_mismatch=config["probe_mismatch"],
        host_columns=PROFILE["metadata"],
    log:
        "logs/reconcile_{set_id}.log",
    script:
        "scripts/reconcile.py"


rule aggregate:
    input:
        summaries=expand("output/primers/{set_id}/final/summary_row.tsv", set_id=SET_IDS),
        lineages=expand("output/primers/{set_id}/final/lineage_rows.tsv", set_id=SET_IDS),
    output:
        summary="output/summary/summary.tsv",
        summary_by_lineage="output/summary/summary_by_lineage.tsv",
    log:
        "logs/aggregate.log",
    script:
        "scripts/aggregate.py"
