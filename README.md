# chichipes

![License](https://img.shields.io/badge/license-GPL--3.0--or--later-blue) ![Go](https://img.shields.io/badge/go-%E2%89%A51.22-00ADD8) ![Dependencies](https://img.shields.io/badge/dependencies-none-brightgreen)

[한국어](README.ko.md)

**chichipes** is a fast, single-binary primer / probe search engine that is a drop-in replacement for
[MFEprimer 4.2.4](https://www.mfeprimer.com) in the commands a TaqMan validation pipeline needs: primer-pair
(amplicon) search with virus coverage, probe / oligo search, dimer and hairpin checks. Same command line, same output
files, same numbers — usually tens to hundreds of times faster.

It is the search engine of the [TaqMan Primer/Probe Validator](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator)
pipeline, published on its own branch for anyone who only wants the engine.

- **No dependencies.** Pure Go standard library; one static executable for Linux, macOS and Windows.
- **No index.** Databases are read directly; `chichipes index` only writes the `.fai`/`.json` files tools expect.
- **Clean-room.** Built from published thermodynamics and black-box comparison of MFEprimer's outputs; no MFEprimer
  code was used or decompiled. See [docs/FINDINGS.md](docs/FINDINGS.md).

## Install

Download the archive for your platform from the
[Releases page](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator/releases) (tags `chichipes-v*`),
unpack it and put `chichipes` on your `PATH`:

```bash
tar xzf chichipes-v0.3.0-linux-amd64.tar.gz
sudo install chichipes-v0.3.0-linux-amd64/chichipes /usr/local/bin/
chichipes version
```

Or build from source (Go 1.22 or newer):

```bash
git clone -b chichipes https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git chichipes
cd chichipes
go build -o chichipes ./cmd/chichipes
go test ./...
```

## Usage

The options are MFEprimer's, so existing scripts only need the binary path changed.

| Command | MFEprimer equivalent | What it does |
|---|---|---|
| `chichipes -d db.fa -i FR.fa -o out.txt --fasta --virus --misMatch 2 [--cpu N]` | `mfeprimer` (spec, virus mode) | Every amplicon a primer set makes in every genome, plus coverage and the genomes it misses |
| `chichipes search db.fa query.fa out.txt -m 8` | `mfeprimer search` | Where an oligo (probe) binds, with edit-distance mismatches, Tm and ΔG |
| `chichipes dimer -i oligos.fa -o dimer.txt` | `mfeprimer dimer` | Primer–primer / primer–probe dimers |
| `chichipes hairpin -i oligos.fa -o hairpin.txt` | `mfeprimer hairpin` | Hairpins |
| `chichipes index -i db.fa` | `mfeprimer index` | Writes `db.fa.fai` / `db.fa.json` (the search itself needs no index) |
| `chichipes reads -i oligos.fa -r reads.fastq.gz -o out` | — | Which sequencing reads carry each primer / probe (below) |
| `chichipes version` | `mfeprimer version` | Version |

Degenerate (IUPAC) primers are expanded, and up to `--misMatch` mismatches in the 3'-end window are tolerated, exactly
as MFEprimer does. Buffer options (`--mono`, `--diva`, `--dntp`, `--oligo`) and cut-offs (`-t`, `-T`, `-s`, `-S`)
are accepted with MFEprimer's defaults.

### Output files

`out.txt` (report with the `Coverage …` line), `out.txt.spec.tsv` (one row per amplicon), `out.txt.amp.fasta`,
`out.txt.virus.Failed.txt` — the same names and formats as MFEprimer — plus `out.txt.borderline.tsv`
(bindings within ±0.02 °C of the Tm cut-off, for reference). Details: [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md) (Korean).

### Using it in place of MFEprimer

Point the tool's MFEprimer path at `chichipes`. In the TaqMan Primer/Probe Validator pipeline that is one line in
`config/config.yaml`:

```yaml
mfeprimer_bin: "/usr/local/bin/chichipes"
```

## How close to MFEprimer is it?

Checked against MFEprimer 4.2.4's own outputs on 182 real primer sets from 17 virus panels:

| Function | Result |
|---|---|
| Primer-pair search | 2,295,176 distinct amplicons: none missing, none extra; F/R Tm and ΔG identical |
| Probe search (`-m 8`, 182 runs) and reference search (`-m 6`, 812 runs) | identical rows |
| Dimer (182 reports) | identical entries and values (equal-ΔG entries may be listed in another order, as MFEprimer itself does with `--cpu` > 1) |
| Hairpin (182 reports) | byte-identical apart from the header and footer |
| Full-precision check (`-j` JSON) | 135,851 bindings, 15,864 dimers, 29,509 hairpins: ΔG and Tm exact |

Where chichipes **deliberately differs** (it is more correct):

- **Coverage beyond 10,000 amplicons.** MFEprimer computes coverage and the Failed list from its first 10,000 products
  only, so highly degenerate sets can be reported as missing genomes they actually amplify. chichipes counts all of them.
- **Duplicate oligos in pooled input.** When every expansion of an oligo already occurs in an earlier oligo, MFEprimer
  skips that oligo's mismatch variants and loses amplicons; chichipes keeps every variant, and reports every primer
  name a sequence belongs to.

Known small differences: a probe Tm can differ by 0.01 °C (rounding of a half-way value); the `ampGC` column
differs by a few G/C in some rows; hairpin ties with equal score are resolved differently in about 1 of 7,000 random
oligos. Non-default buffer conditions and non-default hairpin `--maxl`/`--minl` are not yet verified.

### Speed

| Job | MFEprimer 4.2.4 | chichipes |
|---|---|---|
| One degenerate primer pair against a 663-genome segment | 18.2 s | 0.1 s |
| Whole pipeline run of an Ebola-virus panel (2,900 genomes, 29 steps) | 14 min 58 s | 16 s |

## Inosine (I)

Oligos may contain deoxyinosine, written `I` (not an IUPAC code, but common in degenerate primers). chichipes treats
it as what it is — **one** base that pairs with every base, unevenly — rather than as `N`:

| | `N` | `I` |
|---|---|---|
| What is synthesised | a mix of 4 oligos | one oligo |
| Expanded into variants | ×4 | no |
| Counts as a mismatch | no | no (against A, C, G or T) |
| Tm / ΔG at that position | each variant pairs perfectly | I·C > I·A > I·T ≈ I·G, all weaker than a Watson–Crick pair |

The inosine stacks are the nearest-neighbour parameters of Watkins & SantaLucia (2005), *Nucleic Acids Res.*
33:6258–6267 (Table 2, I·C, I·A, I·T, I·G, I·I). A stack the paper has no value for (an inosine next to a mismatch or
an `N`, or at a dangling end) adds nothing beyond its salt term, the convention already used for `N`. Inosine in the
3′ window is still found (one index key per facing base). For dimer and hairpin detection I·C and I·A count as pairs.
MFEprimer does not accept `I`, so these results have no MFEprimer counterpart; oligos without `I` give byte-identical
results to before.

## Install the dashboards (UI)

`setup.sh` turns chichipes into the full TaqMan Primer/Probe Validator: it installs the pathogen dashboards (panels)
from this repository, points them at chichipes, downloads their genome data and starts them.

```bash
git clone -b chichipes https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git chichipes
cd chichipes
./setup.sh                          # asks which panels; then open http://localhost:8760
./setup.sh --panels lasv,ebov --yes # no questions
```

| Your computer | Command | What happens |
|---|---|---|
| Ubuntu / Debian / other Linux | `./setup.sh` | mamba environment (Miniforge installed if missing) |
| macOS (Apple Silicon or Intel) | `./setup.sh` | mamba environment; `--backend brew` uses Homebrew instead |
| Windows | `.\setup.ps1` in PowerShell | Docker Desktop (default); `.\setup.ps1 -Wsl` runs `setup.sh` inside WSL2 Ubuntu instead |
| Any, with Docker | `./setup.sh --backend docker` | everything in one container; panels and data stay in the install folder |

Snakemake, MAFFT and htslib have no Windows builds, so on Windows the panels run in WSL2 or Docker (both need
virtualisation enabled in the BIOS and administrator rights the first time).

Main options: `--panels LIST|all|everything` (`all` leaves out the very large mpox and influenza panels),
`--dir` (default `~/chichipes-ui`), `--port-base` (front page port, default 8760), `--no-data`,
`--primers FILE` (primer sets exported with the pipeline's `tools/export_primer_sets.py`), `--update`. The install
folder gets `start.sh`, `stop.sh` and `fetch_data.sh`; data download runs in the background
(`logs/data_status.tsv`). The public panels contain no primer sets: add yours on each panel's primer-set page or
with `--primers`.

## Checking sequencing reads (FASTQ)

`chichipes reads` uses FASTQ files the way the primer search uses a genome FASTA: it reports which reads carry each
primer / probe. Single-end, paired-end and long reads; FASTQ or FASTA, gzip or not.

```bash
chichipes reads -i oligos.fa -r sample.fastq.gz -o sample                      # single-end / long reads
chichipes reads -i oligos.fa -1 S_R1.fastq.gz -2 S_R2.fastq.gz -o sample        # paired-end (a pair counts once)
chichipes reads -i oligos.fa -r ont.fastq.gz -o ont -m 3 --indel                # error-prone long reads (ONT, PacBio CLR)
```

Name the oligos `SET_F`, `SET_R`, `SET_P` to also get per-set counts.

| File | Contents |
|---|---|
| `<prefix>.summary.tsv` | per oligo: reads (or read pairs) carrying it, %, how many with 0, 1, 2… mismatches, strand |
| `<prefix>.sets.tsv` | per set: reads with F, R, P, with F and R, with F, R and P |
| `<prefix>.positions.tsv` | per oligo position: the bases the reads have there (A/C/G/T/other/deleted) and % matching — shows which position a variant hits |
| `<prefix>.hits.tsv` | with `--hits`: one line per read and oligo (position, strand, mismatches, read bases) |

`-m` is mismatches per oligo (substitutions; with `--indel`, edits). The whole oligo must lie inside the read. IUPAC
codes and inosine match every base they stand for; a read `N` is a mismatch. Each oligo is cut into m+1 or more
disjoint exact seeds (at least one is untouched when there are at most m errors), and every seed hit is verified
over the whole oligo, so nothing within `-m` is missed: on simulated Illumina single-end (20,000 reads), paired-end
(10,000 pairs) and ONT-like long reads (400 reads, 4 % errors, `--indel`) the reads found were exactly those a
brute-force check finds. About 1 million 150-bp gzipped reads in 2 s on 16 threads.

## Thermodynamics

SantaLucia & Hicks (2004) nearest-neighbour, mismatch, dangling-end and loop parameters, salt correction after
von Ahsen et al. (2001), Tm = 1000·ΔH / (ΔS + R ln(Ct/4)) − 273.15. The exact rules (3'-window lookup, mismatch-run
truncation, end terms, hairpin search) are written up in [docs/FINDINGS.md](docs/FINDINGS.md).

## Repository layout

| Path | Contents |
|---|---|
| `cmd/chichipes/` | command line |
| `internal/thermo/` | thermodynamics: perfect duplex, binding to a target, dimer and hairpin energies |
| `internal/spec/` | primer-pair, probe, dimer and hairpin search and their output files |
| `internal/primers/` | IUPAC expansion, 3'-window mismatch variants and their names |
| `harness/` | the black-box experiments and comparison scripts used to derive and check the rules (need `MFEPRIMER`, `EDITIONS`) |
| `docs/` | `FINDINGS.md` (rules, English), `OUTPUT_FORMAT.md` (output formats, Korean) |
| `scripts/release.sh` | builds the release archives for every platform |

## Releasing

```bash
git tag chichipes-v0.3.1 && git push origin chichipes-v0.3.1
```

Pushing a `chichipes-v*` tag runs `.github/workflows/release.yml`, which tests, builds Linux / macOS / Windows archives with
`scripts/release.sh` and attaches them with `SHA256SUMS.txt` to a GitHub Release. `scripts/release.sh v0.3.1` builds
the same archives locally into `dist/`.

## License and citation

GPL-3.0-or-later ([LICENSE](LICENSE)); third-party notices in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
chichipes is not affiliated with MFEprimer. If you publish results produced in MFEprimer-compatible mode, please also
cite MFEprimer: Wang K. et al., *Nucleic Acids Res.* 2019;47(W1):W610–W613, https://doi.org/10.1093/nar/gkz351.
