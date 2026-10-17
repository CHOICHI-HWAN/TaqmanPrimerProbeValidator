# Data acquisition and the population filter: fetch -> rename headers -> length/N filter -> mfeprimer index.
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 1. fetch --------------------------------------------------------------
rule fetch_raw:
    output:
        fasta="data/raw/{seg}_sequences.fasta",
        meta="data/raw/{seg}_metadata.tsv",
    params:
        seq_url=lambda wc: source_urls(PROFILE, wc.seg)[0],
        meta_url=lambda wc: source_urls(PROFILE, wc.seg)[1],
    log:
        "logs/fetch_raw_{seg}.log",
    shell:
        "bash scripts/fetch_one.sh {params.seq_url} {params.meta_url} {output.fasta} {output.meta} "
        "> {log} 2>&1"


# Clade/lineage from a locally built Auspice JSON (profile lineage.mode = "auspice_clades"): one
# accession -> clade table per segment, merged into the curated metadata by rename_headers.
LINEAGE_FROM_AUSPICE = PROFILE["lineage"]["mode"] == "auspice_clades"

if LINEAGE_FROM_AUSPICE:
    rule lineage_table:
        input:
            auspice=lambda wc: PROFILE["lineage"]["auspice"].format(seg=wc.seg),
        output:
            table="data/lineage/{seg}_lineage.tsv",
        params:
            attr=PROFILE["lineage"]["attr"],
        log:
            "logs/lineage_table_{seg}.log",
        script:
            "scripts/lineage_from_auspice.py"


def rename_inputs(wc):
    d = {"fasta": f"data/raw/{wc.seg}_sequences.fasta", "meta": f"data/raw/{wc.seg}_metadata.tsv"}
    if LINEAGE_FROM_AUSPICE:
        d["lineage_table"] = f"data/lineage/{wc.seg}_lineage.tsv"
    return d


rule rename_headers:
    input:
        unpack(rename_inputs),
    output:
        renamed_fasta="data/renamed/{seg}_sequences.renamed.fasta",
        metadata_curated="data/renamed/{seg}_metadata_curated.tsv",
    params:
        lineage_col=lambda wc: PROFILE["segments"][wc.seg]["lineage_col"],
        columns=PROFILE["metadata"],
        lineage_missing=PROFILE["lineage"].get("missing", "NA"),
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
        max_n_fraction=PROFILE["qc"]["max_n_fraction"],
        exclude_list=PROFILE["qc"]["exclude_list"] or "",
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
        max_n_fraction=PROFILE["qc"]["max_n_fraction"],
        exclude_list=PROFILE["qc"]["exclude_list"] or "",
    log:
        "logs/filter_status_{seg}.log",
    script:
        "scripts/filter_status.py"


# ---- 2a-2. the full raw metadata with the clade/lineage and an explicit QC verdict per record
#            (QC PASS / QC FAIL + reason) -- the "annotated metadata" people read and filter.
rule annotate_metadata:
    input:
        meta="data/raw/{seg}_metadata.tsv",
        curated="data/renamed/{seg}_metadata_curated.tsv",
        status="data/filtered/{seg}_filter_status.tsv",
    output:
        annotated="data/filtered/{seg}_metadata_annotated.tsv",
    params:
        min_len=lambda wc: PROFILE["segments"][wc.seg]["min_len"],
        lineage_label=PROFILE["lineage"]["label"].lower(),
        accession_col=PROFILE["metadata"]["accession"],
    log:
        "logs/annotate_metadata_{seg}.log",
    script:
        "scripts/annotate_metadata.py"


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
