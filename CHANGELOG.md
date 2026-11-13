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
