# Data acquisition and the population filter: fetch -> rename headers -> length/N filter -> mfeprimer index.
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 1. fetch --------------------------------------------------------------
rule fetch_raw:
    output:
        fasta="data/raw/{seg}_sequences.fasta",
        meta="data/raw/{seg}_metadata.tsv",
    params:
        base_url=PROFILE["source"]["base_url"],
    log:
        "logs/fetch_raw_{seg}.log",
    shell:
        "bash scripts/fetch_one.sh {wildcards.seg} {params.base_url} {output.fasta} {output.meta} "
        "> {log} 2>&1"


rule rename_headers:
    input:
        fasta="data/raw/{seg}_sequences.fasta",
        meta="data/raw/{seg}_metadata.tsv",
    output:
        renamed_fasta="data/renamed/{seg}_sequences.renamed.fasta",
        metadata_curated="data/renamed/{seg}_metadata_curated.tsv",
    params:
        lineage_col=lambda wc: PROFILE["segments"][wc.seg]["lineage_col"],
        columns=PROFILE["metadata"],
    log:
        "logs/rename_headers_{seg}.log",
    script:
        "scripts/rename_headers.py"


# ---- 2. filter --------------------------------------------------------------
rule filter:
    input:
        fasta="data/raw/{seg}_sequences.fasta",
        metadata_curated="data/renamed/{seg}_metadata_curated.tsv",
    output:
        filtered_fasta="data/filtered/{seg}_filtered.fasta",
        filtered_metadata="data/filtered/{seg}_filtered_metadata.tsv",
        report="data/filtered/{seg}_filter_report.txt",
        orientation="data/filtered/{seg}_orientation.tsv",
    params:
        min_len=lambda wc: PROFILE["segments"][wc.seg]["min_len"],
        reference_accession=lambda wc: REFGENOME_ACCESSION[wc.seg],
    log:
        "logs/filter_{seg}.log",
    script:
        "scripts/filter.py"


# ---- 2a. per-genome record of the filter above (which raw genomes passed the data QC and, for
#          the rest, why not) -- feeds the "데이터 QC" column of the metadata table. A separate rule
#          (rather than another output of `filter`) so adding it doesn't re-trigger filtering and
#          everything downstream of it.
rule filter_status:
    input:
        fasta="data/raw/{seg}_sequences.fasta",
        metadata_curated="data/renamed/{seg}_metadata_curated.tsv",
        filtered_metadata="data/filtered/{seg}_filtered_metadata.tsv",
    output:
        status="data/filtered/{seg}_filter_status.tsv",
    params:
        min_len=lambda wc: PROFILE["segments"][wc.seg]["min_len"],
    log:
        "logs/filter_status_{seg}.log",
    script:
        "scripts/filter_status.py"


# ---- 2b. mfeprimer index (built once per segment, ahead of the many parallel
#          primer_match/probe_match jobs, so they don't race to build it themselves) ---
rule mfeprimer_index:
    input:
        "data/filtered/{seg}_filtered.fasta",
    output:
        fai="data/filtered/{seg}_filtered.fasta.fai",
        json="data/filtered/{seg}_filtered.fasta.json",
        log_="data/filtered/{seg}_filtered.fasta.log",
        primerqc="data/filtered/{seg}_filtered.fasta.primerqc",
        primerqc_fai="data/filtered/{seg}_filtered.fasta.primerqc.fai",
    threads: config["mfeprimer_cpu"]
    log:
        "logs/mfeprimer_index_{seg}.log",
    script:
        "scripts/build_index.py"
