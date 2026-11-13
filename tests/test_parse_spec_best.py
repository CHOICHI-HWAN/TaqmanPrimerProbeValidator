"""Mismatch counts must come from the GENOME bases, not from mfeprimer's reported primer text.

mfeprimer lists the primer (degenerate ones expanded into ACGT variants) and only reflects mismatches inside the
3' k-mer it seeds on; a mismatch further toward the 5' end is not in that text, so recounting against it gives 0.
"""
from mfe_common import parse_spec_best, reverse_complement
from seqgen import synth

F = synth(21, 1)
R = synth(26, 2)
MID = synth(39, 3)


def write_spec(path, chrom, f_start, f_end, r_start, r_end):
    row = ["Amp_1", chrom, str(f_start), str(r_end), "49.0", str(r_end - f_start + 1),
           "T_F", str(f_start), str(f_end), F, "62.0", "52.0", "-23.0",
           "T_R", str(r_end), str(r_start), reverse_complement(R), "58.0", "38.0", "-23.0", ""]
    path.write_text("#header\n" + "\t".join(row) + "\n")


def genome_with(r_site):
    g = "GGGG" + F + MID + r_site + "CCCC"
    f_start = 5
    f_end = f_start + len(F) - 1
    r_start = f_end + len(MID) + 1
    return g, f_start, f_end, r_start, r_start + len(R) - 1


def test_a_mismatch_outside_the_3p_kmer_is_counted_from_the_genome(tmp_path):
    site = reverse_complement(R)
    site = site[:2] + ("G" if site[2] != "G" else "A") + site[3:]          # 1 substitution near the primer's 5' end
    g, fs, fe, rs, re_ = genome_with(site)
    spec = tmp_path / "x.spec.tsv"
    write_spec(spec, "chr", fs, fe, rs, re_)
    old = parse_spec_best(spec, "T_F", F, "T_R", R)["chr"]                  # reported text only: blind to it
    new = parse_spec_best(spec, "T_F", F, "T_R", R, genome_seqs={"chr": g})["chr"]
    assert (old["f_mismatch"], old["r_mismatch"]) == (0, 0)
    assert (new["f_mismatch"], new["r_mismatch"]) == (0, 1)


def test_a_perfect_site_and_degenerate_primers_stay_zero(tmp_path):
    g, fs, fe, rs, re_ = genome_with(reverse_complement(R))
    spec = tmp_path / "x.spec.tsv"
    write_spec(spec, "chr", fs, fe, rs, re_)
    best = parse_spec_best(spec, "T_F", F, "T_R", R, genome_seqs={"chr": g})["chr"]
    assert (best["f_mismatch"], best["r_mismatch"]) == (0, 0) and best["f_3p_ok"] is not False
    # an IUPAC code in the primer covers either genome base
    deg = F[:5] + "Y" + F[6:]
    spec2 = tmp_path / "y.spec.tsv"
    write_spec(spec2, "chr", fs, fe, rs, re_)
    best2 = parse_spec_best(spec2, "T_F", deg, "T_R", R, genome_seqs={"chr": g})["chr"]
    assert best2["f_mismatch"] == (0 if F[5] in "CT" else 1)


def test_probe_must_lie_between_the_primers():
    from mfe_common import probe_inside_amplicon
    row = {"amplicon_start": "100", "amplicon_end": "220", "probe_start": "140", "probe_end": "170"}
    assert probe_inside_amplicon(row)
    assert not probe_inside_amplicon(dict(row, probe_start="300", probe_end="330"))     # binds elsewhere
    assert not probe_inside_amplicon(dict(row, probe_start="90", probe_end="120"))      # sticks out of the amplicon
    assert probe_inside_amplicon({"amplicon_start": "", "probe_start": ""})              # no coordinates -> not judged


def test_scan_fr_finds_the_pair_and_choose_fr_prefers_the_better_engine():
    from mfe_common import choose_fr, scan_fr
    g = "TTTT" + F + "A" * 60 + reverse_complement(R) + "GGGG"
    own = scan_fr(F, R, g, three_prime_n=2)
    assert own["total"] == 0 and own["amp_start"] == "5" and own["amp_end"] == str(4 + len(F) + 60 + len(R))
    worse = {"f_mismatch": 5, "r_mismatch": 4, "total": 9, "f_3p_ok": True, "r_3p_ok": True}   # a poor pair mfeprimer listed
    assert choose_fr(worse, own) == (own, "genome-scan")
    tie = dict(own, f_tm="62.0")
    assert choose_fr(tie, own) == (tie, "mfeprimer")                                            # tie: keep mfeprimer (has Tm/dG)
    assert choose_fr(None, own) == (own, "genome-scan") and choose_fr(None, None) == (None, "")
    assert scan_fr(F, R, "ACGT" * 30) is None


def test_scan_fr_requires_facing_primers_within_the_size_limit():
    from mfe_common import scan_fr
    same_strand = "TTTT" + F + "A" * 60 + R + "GGGG"                  # R written on the same strand as F: no amplicon
    assert scan_fr(F, R, same_strand) is None
    far = "TTTT" + F + "A" * 2100 + reverse_complement(R)
    assert scan_fr(F, R, far) is None and scan_fr(F, R, far, max_size=3000) is not None


def test_scan_probe_prefers_a_site_inside_the_amplicon_and_choose_probe_takes_the_better():
    from mfe_common import choose_probe, scan_probe
    probe = synth(21, 4)
    probe = probe[:10] + ("A" if probe[10] == "C" else probe[10]) + probe[11:]      # the mismatch test below substitutes a "C" there
    other = "T" * 30 + probe + "T" * 30                      # a perfect site at 31-51, outside the amplicon below
    g = other + "G" * 40 + probe[:10] + "C" + probe[11:] + "G" * 20       # 1-mismatch site further right
    inside = scan_probe(probe, g, amp_start=100, amp_end=200)
    assert inside["mismatch"] == 1 and 100 <= int(inside["start"]) and int(inside["end"]) <= 200
    assert scan_probe(probe, g)["mismatch"] == 0               # anywhere: the perfect one
    assert scan_probe(probe, g, amp_start=1, amp_end=20) is None
    mfe = {"mismatch": 2, "start": "100", "end": "120", "tm": "60", "dg": "-20"}
    assert choose_probe(mfe, inside)[1] == "genome-scan" and choose_probe(inside, dict(inside))[1] == "mfeprimer"
    assert choose_probe(None, None) == (None, "")


def test_scan_combos_finds_cross_set_pairs_and_counts_5prime_mismatches():
    from group_cross_binding import scan_combos
    F2, R2 = synth(22, 5), synth(22, 6)
    oligos = {"A::F": {"set_id": "A", "seq": F}, "A::R": {"set_id": "A", "seq": R},
              "B::F": {"set_id": "B", "seq": F2}, "B::R": {"set_id": "B", "seq": R2}}
    site = reverse_complement(R2)                                  # (R2 as it sits on the genome's top strand)
    site_r2 = site[:2] + next(c for c in "ACGT" if c != site[2]) + site[3:]     # 1 mismatch near B::R's 5' end
    g = "TTTT" + F + "A" * 50 + site_r2 + "GGGG"                   # A's F with B's R: a cross-set amplicon
    combos = scan_combos(oligos, g, three_prime_n=2)
    cross = combos[("A::F", "B::R")]
    assert (cross["f_mismatch"], cross["r_mismatch"]) == (0, 1) and cross["f_set_id"] == "A" and cross["r_set_id"] == "B"
    assert ("A::F", "A::R") not in combos and ("B::F", "B::R") not in combos     # those sites are not in this genome
