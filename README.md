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
tar xzf chichipes-v0.1.0-linux-amd64.tar.gz
sudo install chichipes-v0.1.0-linux-amd64/chichipes /usr/local/bin/
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
git tag chichipes-v0.1.1 && git push origin chichipes-v0.1.1
```

Pushing a `chichipes-v*` tag runs `.github/workflows/release.yml`, which tests, builds Linux / macOS / Windows archives with
`scripts/release.sh` and attaches them with `SHA256SUMS.txt` to a GitHub Release. `scripts/release.sh v0.1.1` builds
the same archives locally into `dist/`.

## License and citation

GPL-3.0-or-later ([LICENSE](LICENSE)); third-party notices in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
chichipes is not affiliated with MFEprimer. If you publish results produced in MFEprimer-compatible mode, please also
cite MFEprimer: Wang K. et al., *Nucleic Acids Res.* 2019;47(W1):W610–W613, https://doi.org/10.1093/nar/gkz351.
