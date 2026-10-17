# Genome-wide nucleotide diversity (pixy) by lineage.
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 5g. per-lineage genome-wide nucleotide diversity (pixy), scanned across the whole
#          S/L segment -- a Nextstrain-entropy-panel-style "diversity by position" curve
#          per lineage, not tied to any one primer set. See scripts/build_diversity_vcf.py
#          for why this needs its own whole-population mafft alignment (with orientation
#          correction) rather than reusing the per-set amplicon MSA machinery. ------------
rule build_diversity_vcf:
    input:
        population="data/filtered/{seg}_filtered.fasta",
        metadata="data/filtered/{seg}_filtered_metadata.tsv",
        ref_fasta="data/refgenome/{seg}_ref.fasta",
    output:
        vcf="output/diversity/{seg}/allsites.vcf.gz",
        populations="output/diversity/{seg}/populations.txt",
        qc="output/diversity/{seg}/qc.json",
    params:
        ref_accession=lambda wc: REFGENOME_ACCESSION[wc.seg],
        chrom_prefix=PROFILE["population"]["chrom_prefix"],
        min_coverage=PROFILE["diversity"]["min_coverage"],
        min_identity=PROFILE["diversity"]["min_identity"],
    threads: 8
    log:
        "logs/build_diversity_vcf_{seg}.log",
    script:
        "scripts/build_diversity_vcf.py"


rule run_pixy:
    input:
        vcf="output/diversity/{seg}/allsites.vcf.gz",
        populations="output/diversity/{seg}/populations.txt",
        qc="output/diversity/{seg}/qc.json",
    output:
        json="output/diversity/{seg}/pixy_pi.json",
    params:
        pooled_label=PROFILE["population"]["pooled_label"],
        window_size=PROFILE["diversity"]["window_size"],
        min_lineage_n=PROFILE["diversity"]["min_lineage_n"],
    log:
        "logs/run_pixy_{seg}.log",
    script:
        "scripts/run_pixy.py"
