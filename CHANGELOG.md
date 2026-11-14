# chichipes v0.6.0 (2026-11-14)

- **The engine is now written in Rust.** Same commands, options and output files; every output is byte-identical to
  v0.4.0/v0.5.0 on 182 real primer sets, 182 probe and 812 reference searches, dimer/hairpin reports (real and 81,000
  synthetic oligos, non-default options), inosine oligos, 83 databases (`index`) and simulated reads (`reads`).
  The MFEprimer reference tests (full-precision binding, dimer and hairpin energies) moved along: `cargo test`.
- **Faster:** the 182 primer searches 326 s → 70 s on one CPU; the largest set 0.84 s → 0.29 s and 140 → 57 MB;
  4,000 long reads with `--indel` 13.5 s → 0.5 s (read batches now also close at 256 kb, so long-read files use every
  thread).
- `reads`: hits are taken in oligo order, so when both mates of a pair carry an oligo with equally many mismatches the
  first mate's hit is always the one counted (Go picked either), and `.hits.tsv` lines come in input order.
- Linux release builds need glibc 2.17 or newer. Building from source needs Rust (`cargo build --release`).
- **`chichipes design` in Rust as well:** the same output files as the Go v0.5.0 engine (cover, shift, tile, exclusion,
  labels on LASV) and about 1.5× faster; sweep results equal the Python engine's in every field on 240 LASV runs.
- `design`: Tm, GC and coverage % are rounded like Python's `round()` (ties to even), so sweep results match the
  Python engine field for field.

# chichipes v0.5.0 (2026-11-14)

- **`chichipes design`: primer and probe design.**
  - `cover`: coverage-first design, ported from the panels' Python engine. All 240 LASV `hightm` runs gave identical
    designs, about 10–25× faster per run.
  - `sweep`: the `design_sweep.py` grids, with compatible `results.jsonl`.
  - `shift`: exhaustive ±N bp / length / IUPAC variants of existing sets.
  - `tile`: a genome-wide sliding-window coverage map plus a greedy panel.
  - `labels`: reporter / quencher and channel checks.
- `--exclude`: non-target product check, and veto of cross-reactive sets in `cover`.
- `--probe-chem linear|zen|mgb`: probe length limits and Tm offset.
- Each set is checked with chichipes Tm, hairpin and 3'-extendable dimers.

# chichipes v0.4.0 (2026-11-14)

- **`setup.sh`: one-command dashboard install.** Installs the chosen pathogen panels and the front page from this
  repository, points them at chichipes, downloads their data in the background and starts them; `start.sh`,
  `stop.sh`, `fetch_data.sh` in the install folder. Backends: mamba (default; Miniforge installed when missing),
  Homebrew, Docker. Per-OS scripts in `install/` (Ubuntu/Debian/WSL, macOS, Windows), `setup.ps1` for PowerShell
  (WSL2 or Docker Desktop). Release archives now include the installer.

# chichipes v0.3.0 (2026-11-14)

- **`chichipes reads`**: checks sequencing reads (FASTQ or FASTA, gzip or not; single-end, paired-end with `-1/-2`,
  long reads, `--indel` for ONT-type errors) for the primers and probes — reads carrying each oligo and how many
  mismatches, per-set counts (F, R, P, F+R, F+R+P), and the bases the reads have at every oligo position. Exact seeds +
  full verification: matches a brute-force check on simulated Illumina SE/PE and ONT reads.

# chichipes v0.2.0 (2026-11-14)

- **Inosine (`I`) in oligos.** One base that pairs with A, C, G and T, never expanded into variants and never a
  mismatch; its binding, probe, dimer and hairpin energies use the deoxyinosine nearest-neighbour parameters of
  Watkins & SantaLucia (2005) Nucleic Acids Res. 33:6258 (I·C > I·A > I·T ≈ I·G). An inosine in the 3' window is still
  found. Oligos without `I` give byte-identical results to v0.1.0.

# chichipes v0.1.0 (2026-11-14)

First release: a clean-room, MFEprimer-4.2.4-compatible primer / probe search engine in one dependency-free binary.

- Primer-pair search with virus coverage (`chichipes -d … -i … -o … --virus --misMatch N`), probe / oligo search
  (`chichipes search`), `dimer`, `hairpin`, and `index` (writes `.fai` / `.json`), with MFEprimer's options and output files.
- Matches MFEprimer 4.2.4 on 182 real primer sets from 17 virus panels (2,295,176 amplicons: none missing or extra).
- Counts every amplicon for coverage (MFEprimer stops at 10,000) and keeps the mismatch variants of duplicate oligos
  in pooled input (MFEprimer drops them).
- `--version` / `version`; the version is stamped into release builds.

Downloads: `chichipes-v0.1.0-<os>-<arch>.tar.gz` (Linux, macOS) or `.zip` (Windows); verify with `SHA256SUMS.txt`.
