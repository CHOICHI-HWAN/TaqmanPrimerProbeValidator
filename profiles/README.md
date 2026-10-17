# Virus profiles

A profile is a JSON file that tells the pipeline which virus it is analysing. `config/config.yaml`
selects one with `profile: "profiles/<id>.json"`. Rules and scripts get their species-specific values
from it instead of hard-coding them (`scripts/virus_profile.py` loads and validates it).

| key | meaning |
|---|---|
| `id`, `name`, `short_name` | identifier (letters/digits/-/_), full name, short label used in titles (e.g. LASV) |
| `source.type`, `source.base_url` | where the data comes from. `nextstrain` = `<base_url>/<segment>/sequences.fasta.zst` and `metadata.tsv.zst` |
| `segments.<key>` | one entry per genome segment (single-segment viruses have one). `label` (display), `min_len` (records shorter than this are dropped), `lineage_col` (metadata column holding the lineage/clade), `reference_accession` (the reference genome for binding maps and the diversity VCF) |
| `metadata.*` | names of the source metadata columns that map to the pipeline's canonical table (`accession`, `date`, `country`, `region`, and `lineage_fallback` used when `lineage_col` is empty). Optional: `host` (default column `host`) and `lab_host` (default `is_lab_host`) feed the host scope (human / animal / unknown; lab hosts are excluded everywhere) |
| `population.pooled_label` | name of the all-genomes-pooled population in the diversity output (e.g. `ALL_LASV`) |
| `population.chrom_prefix` | prefix of the contig name in the all-sites VCF (`LASV` → `LASV_S`) |
| `diversity.window_size` | pixy window (bp); 1 = per-site |
| `diversity.min_coverage` / `min_identity` | a genome needs a real base at ≥ this fraction of reference positions, and ≥ this fraction of those equal to the reference, to enter the diversity analysis |
| `diversity.min_lineage_n` | lineages with fewer genomes are left out of the per-lineage curves (they still count in the pooled curve) |

Species-independent run settings (mismatch limits, CPU cores, representativeness thresholds…) stay in
`config/config.yaml`. The profile contains no primer or standard-material sequences, so it can be shared.
