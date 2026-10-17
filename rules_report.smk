# Static dashboard.
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 6. dashboard --------------------------------------------------------------
rule dashboard:
    input:
        summary="output/summary/summary.tsv",
        summary_by_lineage="output/summary/summary_by_lineage.tsv",
        filter_reports=expand("data/filtered/{seg}_filter_report.txt", seg=SEGMENTS),
        standard_summary="output/summary/standard_summary.tsv",
        standard_matches="output/summary/standard_matches.tsv",
        group_summary="output/summary/group_summary.tsv",
        refgenome_bindings=expand("output/primers/{set_id}/refgenome_binding/binding.tsv", set_id=SET_IDS),
        refgenomes=expand("data/refgenome/{seg}_ref.fasta", seg=SEGMENTS),
    output:
        html="output/dashboard/index.html",
    params:
        primer_mismatch_pass=PRIMER_MISMATCH_PASS,
        probe_mismatch=config["probe_mismatch"],
        virus_short=PROFILE["short_name"],
        source_url=PROFILE["ui"]["web_url"],
        seg_keys=SEGMENTS,
        seg_labels={seg: s["label"] for seg, s in PROFILE["segments"].items()},
        set_seg=SET_SEG,
        set_name=SET_NAME,
        refgenome_accession=REFGENOME_ACCESSION,
    log:
        "logs/dashboard.log",
    script:
        "scripts/dashboard.py"
