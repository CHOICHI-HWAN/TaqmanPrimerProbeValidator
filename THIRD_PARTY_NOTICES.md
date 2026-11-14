# Third-party notices

chichipes links two crates from crates.io, both MIT / Apache-2.0: `flate2` (with `miniz_oxide`; reading gzip
files) and `chrono` (local time in the dimer and hairpin report headers).

**Thermodynamic parameters.** The nearest-neighbour, mismatch, dangling-end, hairpin-loop and tri/tetraloop values
are the published parameters of SantaLucia J. Jr. & Hicks D. (2004), *The thermodynamics of DNA structural motifs*,
Annu. Rev. Biophys. Biomol. Struct. 33:415–440. They were read from the parameter tables distributed with primer3
(https://github.com/primer3-org/primer3) and converted into `data/hairpin_tables.json`; a few entries were
replaced by values measured against MFEprimer's output (see `docs/FINDINGS.md` §7–8).

**Deoxyinosine parameters.** The inosine stacks in `src/thermo.rs` are Table 2 of Watkins N.E. Jr. &
SantaLucia J. Jr. (2005), *Nearest-neighbor thermodynamics of deoxyinosine pairs in DNA duplexes*, Nucleic Acids Res.
33(19):6258–6267, https://doi.org/10.1093/nar/gki918.

**MFEprimer.** chichipes reproduces the command line and output files of MFEprimer 4.2.4
(Wang K. et al., *MFEprimer-3.0: quality control for PCR primers*, Nucleic Acids Res. 2019;47(W1):W610–W613,
https://doi.org/10.1093/nar/gkz351). It is a clean-room implementation: the rules were found only from published
papers, the program's `--help` text and comparisons of its outputs; the MFEprimer binary was never decompiled or
disassembled, and no MFEprimer code is included. MFEprimer is not affiliated with this project. Please cite MFEprimer
as well when results produced in MFEprimer-compatible mode are published.

**Go `math.Log`.** `src/golog.rs` is a translation of Go's `src/math/log.go`, so that every printed Tm and dG matches
the earlier Go releases bit for bit. It is distributed under Go's licence:

```
Copyright 2009 The Go Authors.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.
   * Neither the name of Google LLC nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

The algorithm and constants of that file come from FreeBSD's `e_log.c`:

```
Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.

Developed at SunPro, a Sun Microsystems, Inc. business.
Permission to use, copy, modify, and distribute this
software is freely granted, provided that this notice
is preserved.
```
