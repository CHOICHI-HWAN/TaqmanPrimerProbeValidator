# chichipes

![License](https://img.shields.io/badge/license-GPL--3.0--or--later-blue) ![Rust](https://img.shields.io/badge/rust-stable-B7410E)

[한국어](README.ko.md)

**chichipes** is a fast, single-binary primer / probe search engine that is a drop-in replacement for
[MFEprimer 4.2.4](https://www.mfeprimer.com) in the commands a TaqMan validation pipeline needs: primer-pair
(amplicon) search with virus coverage, probe / oligo search, dimer and hairpin checks. Same command line, same output
files, same numbers — usually tens to hundreds of times faster.

It is the search engine of the [TaqMan Primer/Probe Validator](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator)
pipeline, published on its own branch for anyone who only wants the engine.

- **One executable.** Written in Rust (v0.4 and earlier: Go); a single file for Linux, macOS and Windows, nothing to
  install alongside. The Linux build runs on glibc 2.17 or newer (CentOS 7 / Ubuntu 14.04 onwards).
- **No index.** Databases are read directly; `chichipes index` only writes the `.fai`/`.json` files tools expect.
- **Clean-room.** Built from published thermodynamics and black-box comparison of MFEprimer's outputs; no MFEprimer
  code was used or decompiled. See [docs/FINDINGS.md](docs/FINDINGS.md).

## Install

Download the archive for your platform from the
[Releases page](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator/releases) (tags `chichipes-v*`),
unpack it and put `chichipes` on your `PATH`:

```bash
tar xzf chichipes-v0.5.0-linux-amd64.tar.gz
sudo install chichipes-v0.5.0-linux-amd64/chichipes /usr/local/bin/
chichipes version
```

Or build from source ([Rust](https://rustup.rs) stable and a C linker, `cc`):

```bash
git clone -b chichipes https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git chichipes
cd chichipes
cargo build --release          # -> target/release/chichipes
cargo test --release           # the MFEprimer reference tests run when work/ holds the reference runs
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

The Rust engine (v0.5.0) gives byte-identical output to the Go engine (v0.4.0) and is faster again: the 182 real
primer searches take 70 s instead of 326 s on one CPU (20 s instead of 42 s on 24), the largest set (228,000
amplicons) 0.29 s instead of 0.84 s with 57 MB instead of 140 MB, and 4,000 long reads with `--indel` 0.5 s
instead of 13.5 s.

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

## Designing primers and probes

`chichipes design` picks F / R / probe from a lineage alignment (aligned FASTA, e.g. from MAFFT) so that the genomes
that PASS are as many as possible. A genome passes when each primer has at most `--mm` mismatches with its last
`--tp` bases at the 3' end intact, and the probe matches perfectly. IUPAC codes count as matches for every base they
stand for.

```bash
chichipes design cover -a aligned.fa --top 3                       # coverage-first panel
chichipes design sweep -a aligned.fa --grid mega -o results.jsonl  # thousands of parameter combinations
chichipes design shift -a aligned.fa --sets primers.tsv --shift 50 # improve existing sets
chichipes design tile  -a aligned.fa --step 10                     # coverage map along the genome
chichipes design labels --sets primers.tsv --groups set_groups.tsv # reporter / quencher check
```

- **cover**: for every window of the alignment, the consensus gets up to `--ambig` two-base IUPAC codes, placed
  where they add the most passing genomes. The genomes each candidate passes are kept as a bit set. F, R and probe
  are combined to pass the most genomes, and further sets are chosen for the genomes the earlier ones miss (greedy
  set cover).
  - This is a port of the panels' `scripts/redesign_coverage.py`. On LASV S and L, all 240 runs of the `hightm` grid
    gave identical designs in every field (sequences, coverage, Tm, GC).
  - About 16–21× faster per run: 1.4 s instead of 22.7 s (S), 0.8 s instead of 16.6 s (L), at 2 threads per run.
  - The Rust engine gives the same output files as the earlier Go one (cover, shift, tile, exclusion, labels on LASV) and
    is about 1.5× faster.
- **sweep**: runs the grids of `design_sweep.py` (`wide`, `hightm`, `midtm`, `longoligo`, `mega`, `mega2`,
  `mega3`). It writes the same `results.jsonl` lines and skips runs already in the file.
- **shift** (brute force): moves each oligo of a set up to `--shift` bases either way, in `--step` steps. It tries every
  length of the length ranges and 0 to N IUPAC codes, and scores every variant against every genome. It keeps the
  `--keep` best of each oligo, combines all of them, and reports the best set for each change class (F, R, P, FR,
  FP, RP, FRP).
  - Each primer stays within 3 °C of its original Tm, and the probe at least `--probe-above` above the primers.
  - On the 25 user sets of LASV S, ±8 bases was never worse than `optimize_user_sets.py` (3.2 s instead of 45 s).
    ±50 bases took 4.9 s and raised several sets from about 0 % to 55–68 %.
- **tile** (brute force): evaluates every conserved window of every length with 0 to N IUPAC codes. For each region
  `[x, x + amp_max)`, with x stepping by `--step`, it finds the best set. `tile.tsv` is the resulting coverage map;
  `tile.panel.tsv` holds the `--top` sets that together cover the most genomes. On the whole LASV S alignment
  (655 genomes) this took 1.4 s.
- **Non-target exclusion** (`--exclude other.fa`): other viruses, the host, and so on.
  - With `cover`, a set whose primers would amplify a non-target sequence with the probe inside is skipped, and the
    next best set is taken.
  - Every reported set gets an `off_target` column: products, products with probe, and oligo hits.
  - A primer binds when it has at most `--mm` mismatches with its 3' bases intact; the search uses the read
    search's exact seeds.
  - Checked against 129 New World arenavirus S genomes (no hits) and against LASV itself (every product found).
- **Check of each set**: chichipes Tm (SantaLucia & Hicks 2004, spec-search buffer), the most stable hairpin, and
  dimers. Only dimers at a primer's 3' end count as extendable, because the probe's 3' end carries the quencher.
  Each set also gets a varVAMP-like penalty (lower is better) and flags.

### Fluorophore and quencher

`--probe-chem` sets the probe format:

| Format | Example | Probe length | Probe Tm offset |
|---|---|---|---|
| `linear` (default) | FAM / BHQ1, FAM / TAMRA | ≤ 30 nt (quenching needs the dye and quencher close) | 0 |
| `zen` | FAM / ZEN / IBFQ (internal quencher) | ≤ 40 nt | 0 |
| `mgb` | FAM / MGB-NFQ | 13–22 nt | +10 °C (rough; set `--probe-tm-offset` from the vendor's Tm of a probe you use) |

Every format also follows these rules:

- No G at the probe's 5' end, because it quenches the reporter.
- No run of four or more identical bases.
- The probe goes on the strand with more C than G; coverage does not change, because a perfect-match probe passes on
  either strand.

`design labels` reads the `fluorophore` column. It flags a quencher that does not cover its reporter (for example
FAM / BHQ2), and a group whose sets share a detection channel if they would run in one tube. Vendor spellings such as
`6FAM / MGBNFQ`, `SUN` and `CY5` are understood.

## Engines of the validation pipeline

`chichipes scan` and `chichipes probedesign` read one JSON request on stdin and write the answer on stdout, exactly
as the pipeline's former Go programs `go/scan` and `go/probedesign` did (byte-identical output on 992 requests from
17 editions).

```bash
echo '{"mode":"fr","f":"ACGT...","r":"TTGC...","three_prime_n":2,"max_mm":2,"max_size":2000,"genomes":[{"id":"x","seq":"..."}]}' | chichipes scan
```

Modes: `fr` (best F/R pair per genome), `probe` (best probe window, optionally inside an amplicon), `sites` and
`multi` (every window with at most `max_mm` mismatches, for one or several oligos). The pipeline calls `scan`
through `scripts/mfe_common.py` when its `mfeprimer_bin` points at a chichipes that has the command, and falls back
to `go/scan` and then to the Python scan, which all give the same result.

## Thermodynamics

SantaLucia & Hicks (2004) nearest-neighbour, mismatch, dangling-end and loop parameters, salt correction after
von Ahsen et al. (2001), Tm = 1000·ΔH / (ΔS + R ln(Ct/4)) − 273.15. The exact rules (3'-window lookup, mismatch-run
truncation, end terms, hairpin search) are written up in [docs/FINDINGS.md](docs/FINDINGS.md).

## Repository layout

| Path | Contents |
|---|---|
| `src/main.rs`, `src/flags.rs` | command line (MFEprimer's / Go `flag` syntax) |
| `src/thermo.rs`, `data/hairpin_tables.json` | thermodynamics: perfect duplex, binding to a target, dimer and hairpin energies |
| `src/spec.rs`, `src/search.rs`, `src/dimer.rs`, `src/hairpin.rs` | primer-pair, probe, dimer and hairpin search and their output files |
| `src/primers.rs` | IUPAC expansion, 3'-window mismatch variants and their names |
| `src/reads.rs`, `src/seqio.rs` | sequencing-read check; FASTA / FASTQ (gzip) reading |
| `src/pipe/` | `scan` and `probedesign`: the validation pipeline's whole-genome scan and rule-based probe design (JSON on stdin and stdout; formerly the pipeline's Go programs `go/scan` and `go/probedesign`) |
| `src/design/` | primer / probe design: coverage-first, sweep, shift, tile, exclusion, probe formats, labels |
| `src/golog.rs` | Go's `math.Log`, so every printed value matches the earlier releases bit for bit |
| `harness/` | the black-box experiments and comparison scripts used to derive and check the rules (need `MFEPRIMER`, `EDITIONS`) |
| `docs/` | `FINDINGS.md` (rules, English), `OUTPUT_FORMAT.md` (output formats, Korean) |
| `scripts/release.sh` | packs one built binary into a release archive |

## Releasing

```bash
git tag chichipes-v0.5.1 && git push origin chichipes-v0.5.1
```

`.github/workflows/release.yml` tests and builds on every push to the `chichipes` branch (Linux amd64/arm64 with
`cargo zigbuild` against glibc 2.17, macOS amd64/arm64, Windows amd64; archives kept as workflow artifacts). A
`chichipes-v*` tag also attaches the archives with `SHA256SUMS.txt` to a GitHub Release. Locally:
`cargo build --release && scripts/release.sh v0.5.1 linux amd64 target/release/chichipes` writes the archive to `dist/`.

## License and citation

GPL-3.0-or-later ([LICENSE](LICENSE)); third-party notices in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
chichipes is not affiliated with MFEprimer. If you publish results produced in MFEprimer-compatible mode, please also
cite MFEprimer: Wang K. et al., *Nucleic Acids Res.* 2019;47(W1):W610–W613, https://doi.org/10.1093/nar/gkz351.
