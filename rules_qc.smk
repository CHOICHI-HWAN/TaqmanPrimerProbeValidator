# Per-set QC views: secondary structure, per-genome mismatch table, multiple-sequence alignment.
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 4c. primer/probe secondary-structure QC (self-dimer, cross-dimer, hairpin) --------
# Independent of population matching -- catches oligos that bind EACH OTHER or fold back on
# themselves, which can kill amplification efficiency even when every population-match
# metric (mismatch/Tm/dG against the target) looks perfect.
rule dimer_hairpin_check:
    input:
        manifest=primer_input,
    output:
        json="output/primers/{set_id}/dimer_hairpin/dimer_hairpin.json",
    params:
        f_seq=lambda wc: SET_ROW[wc.set_id]["F_seq"],
        p_seq=lambda wc: SET_ROW[wc.set_id]["P_seq"],
        r_seq=lambda wc: SET_ROW[wc.set_id]["R_seq"],
    threads: config["mfeprimer_cpu"]
    log:
        "logs/dimer_hairpin_{set_id}.log",
    script:
        "scripts/run_dimer_hairpin.py"


# ---- 5e. per-genome F/Probe/R mismatch table (for the interactive alignment/dot view) --
rule alignment_table:
    input:
        db=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered.fasta",
        db_fai=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered.fasta.fai",
        spec=lambda wc: f"output/primers/{wc.set_id}/primer_qc/{wc.set_id}_mm2.txt.spec.tsv",
        metadata=lambda wc: f"data/filtered/{SET_SEG[wc.set_id]}_filtered_metadata.tsv",
        manifest=primer_input,
    output:
        table="output/primers/{set_id}/alignment/mismatch_table.tsv",
        comparison="output/primers/{set_id}/alignment/engine_comparison.tsv",
    params:
        three_prime_n=config["three_prime_consecutive_match"],
    log:
        "logs/alignment_table_{set_id}.log",
    script:
        "scripts/build_alignment_table.py"


# ---- 5f. real multiple-sequence-alignment view (mafft) of every matched genome's
#          amplicon, for the "see all genomes' F/Probe/R binding at once" view --------
rule engine_concordance:
    """mfeprimer vs the whole-genome scan: per set, how many genomes get the same primer-pair verdict from both."""
    input:
        comparisons=expand("output/primers/{set_id}/alignment/engine_comparison.tsv", set_id=SET_IDS),
    output:
        table="output/summary/engine_concordance.tsv",
    log:
        "logs/engine_concordance.log",
    script:
        "scripts/engine_concordance.py"


rule msa_alignment:
    input:
        mismatch_table="output/primers/{set_id}/alignment/mismatch_table.tsv",
        manifest=primer_input,
    output:
        msa="output/primers/{set_id}/msa/msa.json",
    log:
        "logs/msa_alignment_{set_id}.log",
    script:
        "scripts/build_msa.py"
