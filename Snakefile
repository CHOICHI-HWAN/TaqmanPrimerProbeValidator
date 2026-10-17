"""LASV TaqMan primer/probe validation pipeline (v2, Snakemake).

    1. fetch      : download S/L sequences+metadata from Nextstrain, and build a
                     human-readable, lineage-annotated copy of the sequences
    2. filter     : keep only records with no ambiguous 'N' bases and length above the
                     per-segment threshold (S>=3000bp, L>=6000bp)
    3. primer_match: mfeprimer F/R primer-pair search (--virus mode), mismatch 0/1/2
    4. probe_match : mfeprimer single-oligo `search` for the probe, mismatch=0 only
    5. reconcile/aggregate: validated = (primer matched) AND (probe perfect match)
    6. dashboard  : charts + tables -> output/dashboard/index.html

Run from this directory:  snakemake --cores 4
"""
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from mfe_common import (  # noqa: E402
    ensure_primer_ids, load_primers_manifest, sync_input_primers,
)
from standards_common import ensure_standards_manifest, sync_standards  # noqa: E402
from virus_profile import load_profile  # noqa: E402

configfile: "config/config.yaml"

# --- virus profile ----------------------------------------------------------------------
# Everything species-specific (data source, segments + length floors + lineage columns + reference
# accessions, metadata column mapping, diversity QC thresholds) lives in profiles/<id>.json; the
# rules below read it from here instead of hard-coding the virus. See profiles/README.md.
PROFILE = load_profile(Path(workflow.basedir) / config["profile"])

# --- primer input ---------------------------------------------------------------------
# config/primers.tsv is the spreadsheet people edit (name, sequences, fluorophore, proceed...).
# The pipeline never reads it directly: every active primer is copied into input_primers/<id>.tsv
# (segment + the three sequences only) and every rule reads that. Each primer has a short
# permanent id (P001, ...) which is what set_id means everywhere below and what its results folder
# (output/primers/<id>/) is named; the human-facing name is only a label, so renaming a primer
# changes nothing here. A file is rewritten only when its content changed, so editing one primer
# re-runs just that primer.
MANIFEST_PATH = Path(workflow.basedir) / config["primers_manifest"]
GROUPS_SRC = Path(workflow.basedir) / config["groups_manifest"]
INPUT_DIR = Path(workflow.basedir) / "input_primers"
ensure_primer_ids(MANIFEST_PATH)                      # ids for rows added by hand-editing the TSV
sync_input_primers(MANIFEST_PATH, GROUPS_SRC, INPUT_DIR)
PRIMERS = load_primers_manifest(MANIFEST_PATH)        # active (proceed=True) primers only


def primer_input(wc):
    return f"input_primers/{wc.set_id}.tsv"


SET_IDS = [r["set_id"] for r in PRIMERS]
SET_SEG = {r["set_id"]: r["segment"].lower() for r in PRIMERS}
SET_NAME = {r["set_id"]: r["name"] or r["set_id"] for r in PRIMERS}  # display only -- never part of a path
SET_SOURCE = {r["set_id"]: r["source"] for r in PRIMERS}
SET_ROW = {r["set_id"]: r for r in PRIMERS}
SEGMENTS = list(PROFILE["segments"])
MISMATCH_LEVELS = config["primer_mismatch_levels"]
PRIMER_MISMATCH_PASS = config["primer_mismatch_pass"]

# KDCA 2024 reference-standard-material sequences (표준물질): a handful of short (800bp)
# gene-region snippets the physical wet-lab positive controls are derived from. A set is
# only reported "O" if it can be validated (primer mm<=pass AND probe mm=0) against at
# least one of these -- reuses the exact same primer_match/probe_match/reconcile rules,
# just pointed at this tiny db instead of the population db.
# The standards are defined in the spreadsheet config/standards.tsv (side menu "표준물질 입력 / 관리");
# the two files below are generated from it at start-up (rewritten only when their content changes).
ensure_standards_manifest(Path(workflow.basedir) / config["standards_manifest"], Path(workflow.basedir) / "config", SEGMENTS)
sync_standards(Path(workflow.basedir) / config["standards_manifest"], Path(workflow.basedir) / "config")
STANDARD_FASTA = {seg: f"config/kdca_standard_{seg}.fasta" for seg in SEGMENTS}
STANDARD_META = {seg: f"config/kdca_standard_{seg}_metadata.tsv" for seg in SEGMENTS}

# single canonical reference genome per segment, used for the per-oligo binding-site maps
REFGENOME_ACCESSION = {seg: s["reference_accession"] for seg, s in PROFILE["segments"].items()}

# Groups come from the input_primers/ snapshot too: member ids of primers that are still active
# (a member that was unchecked/deleted is already gone, and a group left with none is dropped).
with open(INPUT_DIR / "groups.tsv", newline="", encoding="utf-8") as _fh:
    GROUPS = list(csv.DictReader(_fh, delimiter="\t"))
GROUP_IDS = [g["group_id"] for g in GROUPS]
GROUP_MEMBERS = {g["group_id"]: g["member_ids"].split(";") for g in GROUPS}
GROUP_SEG = {g["group_id"]: g["segment"].lower() for g in GROUPS}
GROUP_LABEL = {g["group_id"]: g["label"] for g in GROUPS}


def group_inputs(wc):
    return [f"input_primers/{m}.tsv" for m in GROUP_MEMBERS[wc.group_id]]


wildcard_constraints:
    seg="|".join(SEGMENTS),
    set_id="|".join(SET_IDS),
    mismatch="|".join(str(m) for m in MISMATCH_LEVELS),
    group_id="|".join(GROUP_IDS),


# Standard material is optional: sets / groups whose segment has none simply have no standard results.
STD_SEGS = [seg for seg in SEGMENTS if (Path(workflow.basedir) / STANDARD_FASTA[seg]).exists()]
SET_IDS_STD = [s for s in SET_IDS if SET_SEG[s] in STD_SEGS]
GROUP_IDS_STD = [g for g in GROUP_IDS if GROUP_SEG[g] in STD_SEGS]


rule all:
    input:
        "output/dashboard/index.html",
        "output/summary/summary.tsv",
        "output/summary/summary_by_lineage.tsv",
        "output/summary/standard_summary.tsv",
        "output/summary/standard_matches.tsv",
        "output/summary/group_summary.tsv",
        "output/summary/gene_summary.tsv",
        "output/summary/engine_concordance.tsv",
        expand("output/primers/{set_id}/alignment/mismatch_table.tsv", set_id=SET_IDS),
        expand("output/primers/{set_id}/msa/msa.json", set_id=SET_IDS),
        expand("output/groups/{group_id}/combo_mismatch_table.tsv", group_id=GROUP_IDS),
        expand("output/groups/{group_id}/standard_combo_table.tsv", group_id=GROUP_IDS_STD),
        expand("output/groups/{group_id}/refgenome_combo_table.tsv", group_id=GROUP_IDS),
        expand("output/primers/{set_id}/standard_alignment/mismatch_table.tsv", set_id=SET_IDS_STD),
        expand("output/diversity/{seg}/pixy_pi.json", seg=SEGMENTS),
        expand("output/standards/{seg}/representativeness.json", seg=STD_SEGS),
        expand("data/filtered/{seg}_filter_status.tsv", seg=SEGMENTS),
        expand("output/primers/{set_id}/dimer_hairpin/dimer_hairpin.json", set_id=SET_IDS),


# ---- rule modules (one file per stage: rules_*.smk) ------------------------------------
include: "rules_data.smk"
include: "rules_matching.smk"
include: "rules_standards.smk"
include: "rules_groups.smk"
include: "rules_reference.smk"
include: "rules_qc.smk"
include: "rules_diversity.smk"
include: "rules_report.smk"
