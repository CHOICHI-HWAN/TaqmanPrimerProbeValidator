# Reference genome: extraction, indexing, per-oligo binding sites and gene summary.
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 5d. reference-genome (AY628203/AY628202) per-oligo binding-site maps --------
rule extract_refgenome:
    input:
        fasta="data/raw/{seg}_sequences.fasta",
    output:
        ref="data/refgenome/{seg}_ref.fasta",
    params:
        accession=lambda wc: REFGENOME_ACCESSION[wc.seg],
    log:
        "logs/extract_refgenome_{seg}.log",
    script:
        "scripts/extract_refgenome.py"


rule refgenome_index:
    input:
        "data/refgenome/{seg}_ref.fasta",
    output:
        fai="data/refgenome/{seg}_ref.fasta.fai",
        json="data/refgenome/{seg}_ref.fasta.json",
        log_="data/refgenome/{seg}_ref.fasta.log",
        primerqc="data/refgenome/{seg}_ref.fasta.primerqc",
        primerqc_fai="data/refgenome/{seg}_ref.fasta.primerqc.fai",
    threads: config["mfeprimer_cpu"]
    log:
        "logs/refgenome_index_{seg}.log",
    script:
        "scripts/build_index.py"


rule refgenome_binding:
    input:
        ref=lambda wc: f"data/refgenome/{SET_SEG[wc.set_id]}_ref.fasta",
        ref_fai=lambda wc: f"data/refgenome/{SET_SEG[wc.set_id]}_ref.fasta.fai",
        manifest=primer_input,
        # Already-validated population-wide F/R/Probe results for this same set -- the
        # reference genome is itself just a regular population member, so if it already
        # has a result here (the common case), reuse it instead of re-deriving mismatch
        # counts/positions from a separate, narrower single-genome search. See
        # run_refgenome_binding.py's module docstring for why this is preferred over the
        # independent mfeprimer search this rule also still runs as a fallback.
        spec=lambda wc: f"output/primers/{wc.set_id}/primer_qc/{wc.set_id}_mm2.txt.spec.tsv",
        alignment=lambda wc: f"output/primers/{wc.set_id}/alignment/mismatch_table.tsv",
    output:
        binding="output/primers/{set_id}/refgenome_binding/binding.tsv",
    params:
        max_mismatch_search=6,
    log:
        "logs/refgenome_binding_{set_id}.log",
    script:
        "scripts/run_refgenome_binding.py"


rule gene_summary:
    input:
        bindings=expand("output/primers/{set_id}/refgenome_binding/binding.tsv", set_id=SET_IDS),
        genes="config/refgenome_genes.tsv",
        manifest=expand("input_primers/{set_id}.tsv", set_id=SET_IDS),
    output:
        summary="output/summary/gene_summary.tsv",
    log:
        "logs/gene_summary.log",
    script:
        "scripts/build_gene_summary.py"
