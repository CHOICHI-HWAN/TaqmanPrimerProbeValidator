# PrimerHunter (third-party, not part of this repository's code)

PrimerHunter 1.0.2 -- Duitama J, Kumar DM, Hemphill E, Khan M, Mandoiu II, Nelson CE. *PrimerHunter: a primer design tool for
PCR-based virus subtype identification.* Nucleic Acids Res 37(8):2483-92 (2009). GNU GPL.

Source: https://dna.engr.uconn.edu/software/PrimerHunter/PrimerHunter-1.0.2.tar.gz (141 KB, downloaded 2026-10-20,
sha256 15dce8b4075cf48158bbb444b76576c848d8e8a1cc13aea144c5d6a7ec2a2a3b).

Build (needs g++ and make; the 2009 sources compile on current gcc with warnings):

    tar -xzf PrimerHunter-1.0.2.tar.gz && cd PrimerHunter && make

`scripts/redesign_primerhunter.py` looks for `tools/primerhunter/PrimerHunter/primerhunter`.
