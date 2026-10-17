# TaqMan Primer/Probe Validator — LASV edition

![License](https://img.shields.io/badge/license-MIT-blue) ![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.27-brightgreen) ![Platform](https://img.shields.io/badge/platform-linux--64-lightgrey)

[한국어 README](README.ko.md)  ·  Nipah edition: branch [`nipah`](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator/tree/nipah)

This tool measures **how well a TaqMan primer/probe set detects the known genomes of a virus** (Lassa virus, LASV), and presents the result in an interactive dashboard.
Use it when you design a new set or check an existing kit: *"what fraction of the known genomes does this set detect, and which lineages does it miss?"*

The LASV edition handles the two genome segments (S and L); genomes and metadata come from the public Nextstrain dataset. A Nipah edition lives on the `nipah` branch.

Detailed documentation (pipeline steps, settings, dashboard features, directory layout): [docs/PIPELINE.md](docs/PIPELINE.md)

## Features

1) **Exact mismatch counting.** mfeprimer's reported primer text hides mismatches toward the 5′ end and sites its seed step misses. Every site is therefore recounted on the **genome bases**, and mfeprimer and an independent genome scan are both run (the better result is kept). A probe only counts when it lies inside the amplicon.

2) **A simple, explicit pass/fail rule.** *Validated* = forward and reverse primer each ≤ 1 mismatch with an intact 3′ end **and** an exact probe match. Genomes that bind but only just miss the rule are shown separately (“caution”).

3) **Host-aware coverage.** Coverage is counted on **human-derived genomes by default**; switch to animal, unknown or all with one click. Laboratory-host genomes are excluded everywhere.

4) **Resolved by lineage (LI–LVII).** Coverage per lineage (LI–LVII), strand-normalised multiple alignment of the amplicon with F/Probe/R positions, reference coordinates and a nucleotide-diversity (π) track.

5) **Panels of several sets.** Group sets, see the union coverage, and detect **cross-binding** (set 1’s forward primer with set 2’s reverse primer).

6) **Standard-material check** including partial matches (red = probe or primer 3′-end mismatch, yellow = tolerated primer mismatch), self-/cross-dimer and hairpin screening, gene-conservation view, Excel import/export.

7) **Self-testing.** `./test.sh` runs synthetic genomes with a known answer through the real mfeprimer; `./test.sh --full` also audits the real results.

## Schematic

![Schematic](docs/schematic/pipeline_schematic.svg)

## Interactive dashboard

![Overview](docs/media/overview.png)

![Alignment](docs/media/alignment.png)

<details><summary>Animated tour (host scope · set detail · standard materials)</summary>

![Host scope](docs/media/host-scope.gif)

![Set detail](docs/media/set-detail.gif)

![Standard materials](docs/media/standard-match.gif)

</details>

*(Vendor names, primer/probe/standard sequences and standard-material names are blurred in these screenshots with `docs/demo_blur.js`.)*

## Installation

### Current version

Code: **v0.19.0**   ·   Tested on Ubuntu 24.04, x86_64

### Requirements

| | |
|---|---|
| conda / mamba / micromamba | creates the environment (snakemake, mafft, pixy, htslib …) |
| **mfeprimer 4.2.4** | **not included and not on conda — install it yourself (steps 3–4 below)** |
| internet | genome download; the optional tree viewer |
| Chrome (optional) | only for the PDF report export |

**1)  Get the code**

```bash
git clone https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git
cd TaqmanPrimerProbeValidator
```

**2)  Create the environment** (asks for the mfeprimer path in step 4; safe to re-run)

```bash
./setup.sh                       # = mamba env create -f environment.yml, then records the paths in config/config.yaml
conda activate primer_validation # only needed to run snakemake yourself; ./start.sh does not need it
```

**3)  Download mfeprimer** — a single static binary (MFEprimer v4.2.4, linux-amd64) from the official site <https://www.mfeprimer.com/>.
It is not redistributed here because its licence does not say that redistribution is allowed.

```bash
chmod +x /path/to/mfeprimer-4.2.4-linux-amd64
/path/to/mfeprimer-4.2.4-linux-amd64 version     # prints the version if it works
```

**4)  Tell the pipeline where mfeprimer is** (any one of these):

- answer the prompt of `./setup.sh`, or run `./setup.sh --mfeprimer /path/to/mfeprimer-4.2.4-linux-amd64 --yes`
- set `mfeprimer_bin: "/path/to/mfeprimer-4.2.4-linux-amd64"` in `config/config.yaml`
- type it on the dashboard’s **Analysis Parameters** page (it is executed once to check it)
- export `MFEPRIMER_BIN=/path/to/mfeprimer-4.2.4-linux-amd64`
- leave it empty to auto-discover a binary under `../*/mfeprimer/` or on `PATH`

**5)  Check the installation** (about 10 seconds; every line should be ✔)

```bash
./test.sh
```

**6)  Optional: tree viewer** (used by the “Nextstrain tree” card; MIT-licensed, so it is not bundled)

```bash
curl -L -o webapp/static/vendor/peartree.bundle.min.js \
  https://github.com/artic-network/peartree/releases/latest/download/peartree.bundle.min.js
```

# Running the tool

*Required input:* one or more **primer/probe sets** (and, optionally, **standard-material** sequences). The genomes are downloaded for you.

**Start the dashboard**

```bash
./start.sh            # http://localhost:8765   (127.0.0.1 only — there is no login)
```

**Enter your primer/probe sets** on the dashboard page **📝 Primer input / manage** (“📝 프라이머 입력 / 관리”; Excel download/upload supported), or edit `config/primers.tsv` (tab-separated):

| column | meaning |
|---|---|
| `primer_name` | display name (a permanent id `P001, P002 …` is assigned on save) |
| `segment` | `S` or `L` |
| `forward_primer`, `reverse_primer`, `probe` | 5′→3′ sequence, IUPAC degenerate bases allowed |
| `fluorophore`, `source`, `remark` | dye, origin (vendor / paper), free note |
| `proceed` | `True` = include in the analysis |

**Download the genomes and analyse everything** — use the buttons **Update Data** and **Re-run Analysis** in the sidebar, or on the command line:

```bash
snakemake --cores 4
```

**Re-run a single set**

```bash
snakemake --cores 4 output/primers/P001/final/validated_ids.txt
```

**Analyse one new set quickly** (no pipeline run): sidebar **+ New Primer Set Analysis**.

**Choose the host scope** — sidebar “Host scope” (human · animal · unknown · all). The same choice is available to scripts:

```bash
curl 'http://localhost:8765/api/sets?host=animal'
```

**Open the dashboard from another computer on a trusted network**

```bash
python3 webapp/server.py --host 0.0.0.0 --port 8765
```

## Output

| file | content |
|---|---|
| `output/summary/summary.tsv` | coverage and validated % per set (primary host scope, column `host_scope`) |
| `output/summary/summary_by_lineage.tsv` | validated genomes per lineage (LI–LVII) |
| `output/summary/group_summary.tsv` | union coverage of each panel |
| `output/primers/<ID>/alignment/mismatch_table.tsv` | per genome: F / probe / R mismatches, 3′-end status, positions |
| `output/primers/<ID>/final/validated_ids.txt` | accessions that pass |
| `output/primers/<ID>/msa/msa.json` | amplicon alignment shown in the dashboard |

## What is not in this repository

Primer, probe, standard-material and genome (FASTA) **sequences and analysis results are never committed** — they live in `config/`, `data/` and `output/`, which are git-ignored,
and `tests/test_no_sequences_in_repo.py` fails if a sequence slips in. Tests generate their nucleotide strings (`tests/seqgen.py`). Vendor names in this public copy are neutral (VendorA …).

# Citation

There is no publication for this tool yet. Please cite the software as below, **and** the tools it builds on.

**This software** — on GitHub, **“Cite this repository”** (right sidebar, built from [CITATION.cff](CITATION.cff)) gives APA and BibTeX; or use:

```bibtex
@software{choi_taqman_primer_probe_validator_2026,
  author  = {Choi, Chihwan},
  title   = {TaqMan Primer/Probe Validator},
  year    = {2026},
  version = {0.19.0},
  url     = {https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator}
}
```

**Tools it builds on:** **MFEprimer-3.0** (Wang K, et al. *Nucleic Acids Res.* 2019, doi:10.1093/nar/gkz351), MAFFT (Katoh & Standley 2013), pixy (Korunes & Samuk 2021), Snakemake (Mölder et al. 2021), Nextstrain (Hadfield et al. 2018).

Third-party licences: [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) · [CHANGELOG.md](CHANGELOG.md) · MIT licence, © 2026 Chihwan Choi.
