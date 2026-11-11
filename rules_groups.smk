# Primer-set groups: union coverage and cross-set (F of one set x R of another) binding.
# (Included from the Snakefile; shares its globals. Script paths are relative to the project root (this file sits next to the Snakefile on purpose: Snakemake resolves a script relative to the file defining the rule, and changing that path text would make every existing result look outdated).)

# ---- 5c. combined-panel ("<<ALL" / "<<P1+P2") union coverage ---------------------
rule group_union:
    input:
        db=lambda wc: f"data/filtered/{GROUP_SEG[wc.group_id]}_filtered.fasta",
        members=lambda wc: [f"output/primers/{sid}/final/validated_ids.txt" for sid in GROUP_MEMBERS[wc.group_id]],
        host_meta=lambda wc: f"data/raw/{GROUP_SEG[wc.group_id]}_metadata.tsv",
    output:
        validated="output/groups/{group_id}/validated_ids.txt",
        summary_row="output/groups/{group_id}/summary_row.tsv",
    params:
        label=lambda wc: GROUP_LABEL[wc.group_id],
        segment=lambda wc: GROUP_SEG[wc.group_id].upper(),
        members=lambda wc: GROUP_MEMBERS[wc.group_id],
        host_columns=PROFILE["metadata"],
    log:
        "logs/group_union_{group_id}.log",
    script:
        "scripts/group_union.py"


rule aggregate_groups:
    input:
        summaries=expand("output/groups/{group_id}/summary_row.tsv", group_id=GROUP_IDS),
    output:
        summary="output/summary/group_summary.tsv",
    log:
        "logs/aggregate_groups.log",
    script:
        "scripts/aggregate_groups.py"


# ---- 5c-2. cross-set binding detection: pools every member set's F/R/Probe into ONE
#            mfeprimer search so unintended cross-set pairings (Set1-F x Set2-R, or
#            Set2's probe landing inside Set1's amplicon) are actually found, instead of
#            assuming each set only ever binds with itself -----------------------------
rule group_cross_binding:
    input:
        db=lambda wc: f"data/filtered/{GROUP_SEG[wc.group_id]}_filtered.fasta",
        db_fai=lambda wc: f"data/filtered/{GROUP_SEG[wc.group_id]}_filtered.fasta.fai",
        metadata=lambda wc: f"data/filtered/{GROUP_SEG[wc.group_id]}_filtered_metadata.tsv",
        manifest=group_inputs,
        gene_summary="output/summary/gene_summary.tsv",
    output:
        table="output/groups/{group_id}/combo_mismatch_table.tsv",
    params:
        members=lambda wc: GROUP_MEMBERS[wc.group_id],
    log:
        "logs/group_cross_binding_{group_id}.log",
    script:
        "scripts/group_cross_binding.py"


# ---- 5c-3/4. same cross-set-aware pooled search as group_cross_binding, but against the
#              tiny standard-material / single reference-genome databases instead of the
#              population db -- so a group's KDCA match and reference binding views are
#              cross-set-aware too, not just each member set checked independently -------
rule group_standard_binding:
    input:
        db=lambda wc: STANDARD_FASTA[GROUP_SEG[wc.group_id]],
        db_fai=lambda wc: STANDARD_FASTA[GROUP_SEG[wc.group_id]] + ".fai",
        metadata=lambda wc: STANDARD_META[GROUP_SEG[wc.group_id]],
        manifest=group_inputs,
    output:
        table="output/groups/{group_id}/standard_combo_table.tsv",
    params:
        members=lambda wc: GROUP_MEMBERS[wc.group_id],
        tag="standard_FR",
    log:
        "logs/group_standard_binding_{group_id}.log",
    script:
        "scripts/group_standard_refgenome_binding.py"


rule group_refgenome_binding:
    input:
        db=lambda wc: f"data/refgenome/{GROUP_SEG[wc.group_id]}_ref.fasta",
        db_fai=lambda wc: f"data/refgenome/{GROUP_SEG[wc.group_id]}_ref.fasta.fai",
        manifest=group_inputs,
    output:
        table="output/groups/{group_id}/refgenome_combo_table.tsv",
    params:
        members=lambda wc: GROUP_MEMBERS[wc.group_id],
        tag="refgenome_FR",
    log:
        "logs/group_refgenome_binding_{group_id}.log",
    script:
        "scripts/group_standard_refgenome_binding.py"
