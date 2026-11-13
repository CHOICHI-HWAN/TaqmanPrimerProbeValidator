# Third-party notices

chichipes contains no third-party code; it uses only the Go standard library.

**Thermodynamic parameters.** The nearest-neighbour, mismatch, dangling-end, hairpin-loop and tri/tetraloop values
are the published parameters of SantaLucia J. Jr. & Hicks D. (2004), *The thermodynamics of DNA structural motifs*,
Annu. Rev. Biophys. Biomol. Struct. 33:415–440. They were read from the parameter tables distributed with primer3
(https://github.com/primer3-org/primer3) and converted into `internal/thermo/hairpin_tables.json`; a few entries were
replaced by values measured against MFEprimer's output (see `docs/FINDINGS.md` §7–8).

**Deoxyinosine parameters.** The inosine stacks in `internal/thermo/inosine.go` are Table 2 of Watkins N.E. Jr. &
SantaLucia J. Jr. (2005), *Nearest-neighbor thermodynamics of deoxyinosine pairs in DNA duplexes*, Nucleic Acids Res.
33(19):6258–6267, https://doi.org/10.1093/nar/gki918.

**MFEprimer.** chichipes reproduces the command line and output files of MFEprimer 4.2.4
(Wang K. et al., *MFEprimer-3.0: quality control for PCR primers*, Nucleic Acids Res. 2019;47(W1):W610–W613,
https://doi.org/10.1093/nar/gkz351). It is a clean-room implementation: the rules were found only from published
papers, the program's `--help` text and comparisons of its outputs; the MFEprimer binary was never decompiled or
disassembled, and no MFEprimer code is included. MFEprimer is not affiliated with this project. Please cite MFEprimer
as well when results produced in MFEprimer-compatible mode are published.
