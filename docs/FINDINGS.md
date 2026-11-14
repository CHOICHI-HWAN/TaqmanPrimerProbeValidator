# Black-box findings about MFEprimer 4.2.4

Clean-room rule: everything here comes from published papers, the program's own `--help` text, and comparing its
outputs on inputs we generated. The binary is never decompiled or disassembled.

Reference binary: `mfeprimer-4.2.4-linux-amd64` (downloaded from the official release page).

## 1. Perfect-match thermodynamics (`mfeprimer thermo`) — reproduced

Data: 3,000 random oligos (12–40 nt, GC 0–100 %) × 5 buffer conditions (`work/r_*.thermo.tsv`).

| item | model |
|---|---|
| nearest-neighbour table | SantaLucia & Hicks 2004 (AA/TT dH −7.6, dS −21.3; rest same as 1998) |
| initiation | dH +0.2, dS −5.7, plus terminal A·T penalty dH +2.2, dS +6.9 per A/T end |
| self-complementary | dS −1.4, and Tm uses ln(Ct) instead of ln(Ct/4) |
| salt | Na_eq = mono + 120·√(diva − dNTP) (von Ahsen 2001), applied as dS += 0.368·(n−1)·ln(Na_eq [M]) |
| dG | dH − 310.15·dS/1000 (37 °C) |
| Tm | 1000·dH / (dS + R·ln(Ct/4)) − 273.15, R = 1.987, Ct = oligo nM · 1e-9 |

Error vs reference: Tm max 0.01 °C, dG max 0.00 kcal/mol (rounding level) on every condition.
Script: `harness/thermo_fit.py`.

## 2. Primer search in `--virus --misMatch m` mode — observed

- Each input oligo is expanded first: degenerate (IUPAC) bases into all concrete sequences (`P001_F.1`, `.2`, …),
  then every concrete sequence into variants with up to `m` substitutions inside the 3'-end window
  (`--misStart 1` … `--misEnd 9`, i.e. the last k = 9 bases): `P001_F.1.1` … For m = 2 on one sequence:
  9·3 single + C(9,2)·9 double = 351 variants. Each variant is looked up exactly by its 3'-end k-mer.
- Outside the 3' window mismatches are not enumerated; they are only limited by the binding Tm cut-off
  (`-t/--tm`, default 30 °C).
- The binding dG / Tm reported per amplicon in `.spec.tsv` differs from the oligo's own `thermo` value even for a
  perfect site (−23.12 vs −22.84 for a 23-mer): the target context (dangling ends, mismatches) enters the binding
  energy. → Section 3 (in progress).

## 3. Binding thermodynamics against a target — reproduced (exact, SantaLucia tables)
Full-precision data: `mfeprimer -j` reruns of the 856 synth.py experiment runs (`harness/bind_json.py` ->
`work/bind_all.json`, 135,851 distinct F and R bindings incl. mismatches, N/ambiguity codes, genome edges and cut
duplexes). Model `harness/bind_tables.py` / `thermo::binding` (src/thermo.rs), no fitted values:
- initiation 0.2 / -5.7; one stack per adjacent column pair: perfect = NN (section 1), mismatch = the mismatch-stack
  table of §7/§8, any N (every ambiguity code is N) = 0 but still one salt unit.
- **Cut**: a mismatch run (N counts) touching the 5' end, or the 3'-most run of >= 3 mismatches, is dropped up to its
  last column, which stays as the duplex's first column; that 5' end is then blunt (no dangle, no terminal-AT term).
  mfeprimer prints the shortened duplex (Qseq/Sseq) in that case.
- Ends: a dangling template base (half a salt unit; an N flank gives only the half unit); at the genome end a
  terminal-AT term (2.2 / 6.9) for A/T. **Forward primer**: 5' dangle read with the template base uncomplemented
  (`dangle half 0 [comp f0][f0][l5]`), 3' dangle taken from the primer's own last base (`half 1 [f3][comp f3][f3]`),
  so r3 never matters. **Reverse primer**: read on the plus strand (rc(primer) over the plus-strand site), so each
  mismatch is scored as the complementary one (A·G -> T·C); dangles `half 1 [comp f0][f0][comp l5]` and
  `half 0 [f3][comp f3][comp r3]`.
- salt 0.368 ln[Na+] per unit; Tm = 1000 dH / (dS + R ln(Ct/4)) - 273.15 with R = 1.9872.
135,850 of 135,851 bindings exact in dG and Tm (the other sits in a poly-G run where the harness cannot place the
site unambiguously). The earlier fitted model (`harness/fit_binding.py`, `harness/params.json`, ±0.005-0.01) is kept
in the harness for the record only.

## 4. Detection rule — reproduced (F primer side)

A designed amplicon is reported iff
(a) the 3'-end 9-nt window has <= misMatch mismatches, and
(b) the binding Tm (section 3, using the variant that matches the window) is >= `--tm` (30 °C).
Hold-out + two training sets: 26,600 records, 0 misclassified (found and not-found).

## 4b. R primer side — reproduced
- `.spec.tsv` prints rpSeq on the plus strand (reverse complement of the oligo).
- Mismatch stacks are read on the other strand: key = rc(primer dinuc) / rc(site dinuc), sharing the F table.
- End terms differ from F: 5' end keyed by (r0, s0, l5) as for F but a separate table, and the 3' end keyed by
  (last base, the genome base beyond the 3' end) — on F that 3' flank has no effect.
- Truncation (>= 3 consecutive mismatches) and detection rule as for F.
Joint F+R fit (39,505 records, 628 keys): every set at rounding level. Held-out R random (7,069 found): Tm max 0.020,
dG max 0.010; detection over R hold-out + tandem: 0 misclassified.

## 4c. Genome ambiguity codes — reproduced
- A site whose 3' window holds any non-A/C/G/T letter is never found (the window lookup is exact).
- Elsewhere in the site or as a flanking base, every IUPAC ambiguity code (R, Y, K, M, S, W, B, D, H, V) is scored
  exactly like N: N gets its own stack / end parameters (fitted), and counts as a mismatch for the truncation rule.
  Check: 15,418 ambiguity-code bindings scored with the N parameters, Tm max 0.011, dG max 0.008.

## 4d. Search set and amplicon rules (from the 182 real reports and spec tables)
- Expansion and variant naming: §OUTPUT_FORMAT; an expansion whose sequence equals an earlier variant is not searched
  again (its sites are reported under the earlier name).
- Two primers of a product may overlap by at most k-1 = 8 bases.
- Sites touching the genome ends are searched; the missing flank is its own parameter value.
- **Coverage (`--virus`) uses only the first `-M` (10,000) products** after mfeprimer's internal sort/dedup, and the
  report's "amplicon number" is then 10000. Genomes whose products all fall beyond that cut are listed as Failed
  although `.spec.tsv` has products for them. In the current pipeline data this under-reports coverage in 22 of 182
  sets (e.g. MPXV P002 98.09 % vs 99.95 % actual: 149 genomes wrongly failed; LASV P048 94.57 % vs 97.89 %).
  The engine counts every product (the correct coverage); the exact mfeprimer cut is not reproduced.
- Beyond 10,000 products, mfeprimer also stops de-duplicating: the same product appears again under the alternative
  names of identical sequences (CCHFV P032: 37,562 rows, 27,164 distinct products). Comparison is therefore done on
  distinct products (chrom, start, end, fpSeq, rpSeq).

## 4e. Real-data check (2026-11-14, table model)
Full re-run of all 182 primer searches of 17 panels against mfeprimer's stored outputs:
| measure | result |
|---|---|
| distinct products | 2,295,176 in mfeprimer; 0 missing, 0 extra (fitted model: 32 missing, 7,546 extra) |
| fpTm / fpDg / rpTm / rpDg of matched rows | identical (no difference even at the printed 2 decimals) |
| Failed-genome list | identical in 160 of 182; the other 22 are exactly the sets with more than 10,000 products, where mfeprimer builds the list from the first 10,000 only (section 4d) |
| ampGC | differs by up to a few GC counts in ~58 % of rows (pre-existing, see below) |

Pooled multi-set runs (the pipeline's `group_cross_binding.py` puts every member set's F/R into one run; ids
`SET::F`): BDBV 8 sets (16 oligos) and LASV 6 sets (12 oligos) — the sets of products (genome, positions, primer
sequences; all F/R/F-F/R-R combinations across sets) and their Tm/dG are identical (186,465 and 370,344 products).
Only the labels differ: when several names give the same concrete sequence (two sets with the same primer, or two
degenerate expansions/variants giving one sequence) the engine reports every name, while mfeprimer reports one of them
per product and, past its 10,000-product sort, sometimes another (so one sequence shows up under two names).
Decision (2026-11-14, user): the engine keeps reporting every name — the group check needs to see which sets pair.

LASV groups as the pipeline defines them (`config/set_groups.tsv` of the user's setup; members by set name, A/B sets
only; `harness/run_groups.sh`, `harness/compare_group_labels.py`): G004, G006, G007, G008 (4 sets), G009, G010, G011
— products, Tm/dG and the (F set, R set) pairs per product identical (11,193 to 182,038 products each).
G012 (three sets of one vendor, L segment): mfeprimer misses 2,904 products the engine finds (all with the third set's R primer). Cause,
an mfeprimer defect in pooled input: when **every expansion of an oligo already occurs as a sequence of an earlier
oligo** (identical primer, or variants that coincide), mfeprimer lists the expansions but generates **none of that
oligo's mismatch variants** (P009 R: 8 instead of 1,412, of which 1,092 are new sequences), so products that need those
variants are lost. Same pattern in the BDBV and LASV pools (an identical P004 primer: 1 / 16 entries instead of 352 /
1,420; no products lost there because the variants also coincide). The engine keeps every variant (not reproduced).

ampGC (open): mfeprimer's amplicon in the JSON output is the genome region (sites as in the genome, so mismatched
primer bases do not count) and its GC% counts S as G/C; on synthetic data that rule matches 31,909 of 31,909
amplicons, but on the real tables it matches only ~57 % of rows (the engine's primer+middle+primer rule ~69 %).

## 6. `mfeprimer search` (probe, reference-genome binding) — reproduced
Experiments: `harness/synth_search.py` (substitutions, indels, ambiguity codes, long queries; both strands).
- A hit needs the query's **5'-end and 3'-end k-mers (k = 9) both exact** on the genome (either strand; an ambiguity
  code inside either k-mer means no hit). Every pair of such anchors whose span is within ±m of the query length is
  aligned to the query; the **MisMatch column is the edit distance (Levenshtein)** and the hit is kept when <= `-m`.
- Span length = query length: Tm and dG = the perfect duplex of the query with each mismatched stack replaced by the
  mismatch-stack table as the primer search (read as written, N = 0; no end/flank terms, no truncation of runs),
  even when the edit distance is smaller than the Hamming count. Span length differs (gapped): Tm = dG = 0.00.
- Ambiguity codes in the hit count as mismatches and are scored as N; on the minus strand they are complemented with
  the IUPAC table (M->K, B->V, ...).
- Naming as the primer search (NAME, or NAME.k for degenerate queries); coordinates 0-based half-open.
Validation: synthetic 1,657 hits; all 182 pipeline probe searches (-m 8, 65,962 hits); 812 reference-genome
searches (-m 6, 589 hits) — every row identical in position, strand, mismatch, Tm, dG and sequences.

## 7. `mfeprimer dimer` — reproduced (exact)
Experiments: 12 × 300 random oligos (47,689 dimers), `harness/dh_parse.py`; full-precision dG/Tm from `dimer -j`
(`work/dh/dj*.json`, 15,864 dimers).
- Every pair of expanded oligos (later one printed first as A; each oligo also with itself) is tried at every shift
  of B written 3'->5' under A. Per shift the candidate is the whole span from the first to the last complementary
  column; it is a dimer when mismatches in the span <= `-m` (2), Score = pairs - mismatches >= `-s` (5) and dG <= `-d`
  (0). One shift = at most one entry. (135,450 oligo pairs: identical dimer sets.)
- dG / Tm use the same tables as the hairpin (§8), no fitted values: initiation 0.2 / -5.7; one stack per adjacent
  column pair in the span (perfect: the binding model's NN; one or two mismatched columns: primer3's SantaLucia
  mismatch-stack table with mfeprimer's few differing entries, two mismatches in a row mostly 0); each end: a terminal
  mismatch stack when both strands overhang, a dangling end (half a salt unit) when only one does, nothing when blunt;
  no terminal-AT term; salt 0.368 ln[Na+] per stack; symmetry -1.4 (and Tm at Ct instead of Ct/4) when the span reads
  the same on both strands **and** the overhanging neighbours mirror (always so for an oligo with itself).
  Tm = 1000 dH / (dS + R ln(Ct/4)) - 273.15 with R = 1.9872. 15,864 dimers: dG max error 4e-14, Tm exact.
  (The earlier fitted end terms, `harness/dimer_ends.json`, `harness/dimer_model.py`, were an approximation of the dangles/terminal
  mismatches.)
- List sorted by dG ascending. Order inside equal dG is not reproducible: with `--cpu` > 1 mfeprimer's own order
  changes between runs (the pipeline uses `--cpu 4`); it decides which equal-dG picture the pipeline shows.
- Report text: oligo table `%-31s%-36s%8d%8.2f%8.2f%10.2f` with the table Tm at R = 1.9872 (the `thermo` command uses
  1.987); `Dimer List (0)` is followed by `No dimer found.`.
- Validation: all 12 synthetic sets and all 182 pipeline dimer reports: same entries and values (no missing, no
  extra); 116 reports byte-identical up to the parameter block, the other 66 differ only in the order of equal-dG
  entries; the pipeline's per-role-pair summary identical in all 182.

## 8. `mfeprimer hairpin` — reproduced (exact energies; structure search bar rare ties)
Experiments: `harness/synth_hairpin.py` (stems with chosen loops/flanks/mismatches, random oligos; `-j` gives full
precision), ~74,000 oligos over default and non-default options; `harness/hp_model.py`, `harness/hp_energy.py`.
- One hairpin (or none) per expanded oligo. Structure line: `/` 5' arm, `\` 3' arm, `-` unpaired; all pairs lie on
  one diagonal i + j = K (a "gap" is a mismatched column, so gaps are symmetric).
- Search, per diagonal K: walk the columns (c, K-c) from the outer end inwards while at least `--minl` bases stay
  unpaired. Mismatched columns are gaps; when a segment collects more than `--maxgap` gaps it is emitted at that
  column with its end gaps dropped and scored by its pairs alone, and a new segment starts after that column (this
  is why `n` mismatched outer columns cost `n mod 3`). The segment that reaches the loop keeps its leading gaps,
  drops loop-side pieces (gaps plus the helix after them) that do not raise its score, and scores pairs - gaps.
- Choice: candidates in the order they are found (column, then diagonal); the first one with the highest score.
  The cut-offs (`-s` 5, `-d` 0: dG <= d, `-t` 0: Tm >= t) are applied to that one only, so a weaker alternative is
  never reported instead. `--mindh` has no visible effect; `--maxl` matters only when the loop would exceed it (only
  the default 12 is verified; with 8 or minl 5 about 2-3 % of oligos differ).
- Energy (no fitted values): initiation 0.2 / -5.7; NN / mismatch stacks over the paired columns as in §7; outer end
  as a dimer end (terminal mismatch / dangle / nothing); at the loop end one more stack read as 5'-X X-3' / 3'-Y b-5'
  (X·Y the closing pair, b the loop base 5' of Y — the loop base after X is never read); salt per stack (half per
  dangle); dG adds the SantaLucia hairpin loop term (3:3.5 4:3.5 5:3.3 6:4.0 7:4.2 8:4.3 9:4.5 10:4.6; longer loops
  take the next listed size up, 12:5.0 14:5.1 16:5.3 18:5.5 20:5.7 25:6.1 30:6.3, none from 30 on) and the
  tri/tetraloop bonus (primer3's table, bonus rounded to 0.1, GTTTTC +0.5). Loop and bonus do not enter Tm.
  Tm = 1000 dH / (dS + R ln(Ct/4)) - 273.15, R = 1.9872. 29,509 hairpins: dG and Tm exact.
- Table entries where mfeprimer differs from primer3's files (measured): TA/AA dH 4.7, CT/CA and AC/TC dH 0,
  TT/TG and GT/TT 5.8/16.3, TA/TG and GT/AT -0.1/-1.7; all other tandem mismatches 0.
- List sorted by dG ascending; `Hairpin List (0)` is followed by `No hairpins found.`.
- Known residual: among candidates with the same highest score the found-first rule picks another one than
  mfeprimer in ~1 of 7,000 random oligos (10 of 71,000 with default options); not explained by position, dG or Tm.
- Validation: all 182 pipeline hairpin reports (14 hairpins) byte-identical up to the parameter block; default-
  option synthetic runs: 18,500 oligos, 2 differ (the tie above).

## 5. Still to check
- Binding terms under non-default buffer conditions (the pipeline uses the defaults).
- ampGC rule on real data (§4e).
- Reverse primer at the genome end beyond its 3' end (terminal-AT term assumed by symmetry, never observed).
- The hairpin equal-score tie-break (§8).
