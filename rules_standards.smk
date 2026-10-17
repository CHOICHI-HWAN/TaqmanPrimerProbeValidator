# Reference-standard-material (표준물질) matching and how representative the standards are of the population.
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 5b. 표준물질(KDCA reference-standard-material) match -----------------------
# Exactly the same primer_match/probe_match/reconcile logic as steps 3-5, just pointed
# at the tiny (4- or 2-record) standard-material db instead of the population db.
rule standards_index:
    input:
        lambda wc: STANDARD_FASTA[wc.seg],
    output:
        fai="config/kdca_standard_{seg}.fasta.fai",
        json="config/kdca_standard_{seg}.fasta.json",
        log_="config/kdca_standard_{seg}.fasta.log",
        primerqc="config/kdca_standard_{seg}.fasta.primerqc",
        primerqc_fai="config/kdca_standard_{seg}.fasta.primerqc.fai",
    threads: config["mfeprimer_cpu"]
    log:
        "logs/standards_index_{seg}.log",
    script:
        "scripts/build_index.py"


rule standard_primer_match:
    input:
        db=lambda wc: STANDARD_FASTA[SET_SEG[wc.set_id]],
        db_fai=lambda wc: STANDARD_FASTA[SET_SEG[wc.set_id]] + ".fai",
        manifest=primer_input,
    output:
        report="output/primers/{set_id}/standard_qc/{set_id}_mm{mismatch}.txt",
        coverage="output/primers/{set_id}/standard_qc/{set_id}_mm{mismatch}_coverage.tsv",
    threads: config["mfeprimer_cpu"]
    log:
        "logs/standard_primer_match_{set_id}_mm{mismatch}.log",
    script:
        "scripts/run_primer.py"


rule standard_probe_match:
    input:
        db=lambda wc: STANDARD_FASTA[SET_SEG[wc.set_id]],
        db_fai=lambda wc: STANDARD_FASTA[SET_SEG[wc.set_id]] + ".fai",
        manifest=primer_input,
    output:
        report="output/primers/{set_id}/standard_qc/{set_id}_probe_mm0.txt",
        matched_ids="output/primers/{set_id}/standard_qc/{set_id}_probe_matched_ids.txt",
    params:
        probe_mismatch=config["probe_mismatch"],
    log:
        "logs/standard_probe_match_{set_id}.log",
    script:
        "scripts/run_probe.py"


rule standard_reconcile:
    input:
        db=lambda wc: STANDARD_FASTA[SET_SEG[wc.set_id]],
        metadata=lambda wc: STANDARD_META[SET_SEG[wc.set_id]],
        mismatch_table="output/primers/{set_id}/standard_alignment/mismatch_table.tsv",
    output:
        validated="output/primers/{set_id}/standard_final/validated_ids.txt",
        probe_dropout="output/primers/{set_id}/standard_final/probe_dropout_ids.txt",
        primer_failed="output/primers/{set_id}/standard_final/primer_failed_ids.txt",
        summary_row="output/primers/{set_id}/standard_final/summary_row.tsv",
        lineage_rows="output/primers/{set_id}/standard_final/lineage_rows.tsv",
    params:
        segment=lambda wc: SET_SEG[wc.set_id].upper(),
        source=lambda wc: SET_SOURCE[wc.set_id],
        primer_mismatch_pass=PRIMER_MISMATCH_PASS,
        probe_mismatch=config["probe_mismatch"],
    log:
        "logs/standard_reconcile_{set_id}.log",
    script:
        "scripts/reconcile.py"


if SET_IDS_STD:
    rule standard_aggregate:
        input:
            summaries=expand("output/primers/{set_id}/standard_final/summary_row.tsv", set_id=SET_IDS_STD),
            lineages=expand("output/primers/{set_id}/standard_final/lineage_rows.tsv", set_id=SET_IDS_STD),
        output:
            summary="output/summary/standard_summary.tsv",
            summary_by_lineage="output/summary/standard_matches.tsv",
        log:
            "logs/standard_aggregate.log",
        script:
            "scripts/aggregate.py"
else:
    # no standard material for any segment: the two tables still exist (header only) so everything downstream
    # (dashboard, web app) works unchanged
    rule standard_aggregate:
        output:
            summary="output/summary/standard_summary.tsv",
            summary_by_lineage="output/summary/standard_matches.tsv",
        run:
            with open(output.summary, "w", encoding="utf-8") as fh:
                fh.write("set_id\tsegment\tsource\ttotal_genomes\tprimer_coverage_mm0_pct\tprimer_coverage_mm1_pct\t"
                         "primer_coverage_mm2_pct\tprobe_coverage_mm0_pct\tprimer_mismatch_used\tprimer_matched\t"
                         "probe_dropout\tvalidated\tvalidated_pct\n")
            with open(output.summary_by_lineage, "w", encoding="utf-8") as fh:
                fh.write("set_id\tsegment\tlineage\tvalidated_count\n")


# ---- 5b-2. per-material F/Probe/R mismatch detail for the standard-material db, same
#            shape/script as the population alignment_table but with every material kept
#            (include_all) instead of dropped when it fails to bind entirely ------------
rule standard_alignment_table:
    input:
        db=lambda wc: STANDARD_FASTA[SET_SEG[wc.set_id]],
        db_fai=lambda wc: STANDARD_FASTA[SET_SEG[wc.set_id]] + ".fai",
        # mfeprimer only writes the .spec.tsv/.amp.fasta siblings when it finds at least
        # one amplicon -- against the tiny 2-4 record standard db a primer set can
        # legitimately find NONE, so those siblings can't be a hard Snakemake input (the
        # DAG would refuse to build). The always-present main report is the real
        # dependency; the script derives the (maybe-missing) sibling paths from it.
        report=lambda wc: f"output/primers/{wc.set_id}/standard_qc/{wc.set_id}_mm2.txt",
        metadata=lambda wc: STANDARD_META[SET_SEG[wc.set_id]],
        manifest=primer_input,
    output:
        table="output/primers/{set_id}/standard_alignment/mismatch_table.tsv",
        comparison="output/primers/{set_id}/standard_alignment/engine_comparison.tsv",
    params:
        include_all=True,
        three_prime_n=config["three_prime_consecutive_match"],
    log:
        "logs/standard_alignment_table_{set_id}.log",
    script:
        "scripts/build_alignment_table.py"


# ---- 5h. how representative is the standard material of the real population? ------
#          (standard-wide match, per-set assay-site match, diversity of the location)
rule standard_representativeness:
    input:
        standards=lambda wc: STANDARD_FASTA[wc.seg],
        ref_fasta="data/refgenome/{seg}_ref.fasta",
        vcf="output/diversity/{seg}/allsites.vcf.gz",
        populations="output/diversity/{seg}/populations.txt",
        pixy="output/diversity/{seg}/pixy_pi.json",
        bindings=lambda wc: [f"output/primers/{s}/refgenome_binding/binding.tsv" for s in SET_IDS if SET_SEG[s] == wc.seg],
    output:
        json="output/standards/{seg}/representativeness.json",
    params:
        identity=config["standard_rep_identity"],
        min_coverage=config["standard_rep_min_coverage"],
        site_tol_primer=config["standard_rep_site_tol_primer"],
        site_tol_probe=config["standard_rep_site_tol_probe"],
        good=config["standard_rep_good"],
        critical=config["standard_rep_critical"],
    log:
        "logs/standard_representativeness_{seg}.log",
    script:
        "scripts/standard_representativeness.py"
