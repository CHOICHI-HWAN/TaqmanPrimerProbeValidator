"""Regression tests for scripts/mfe_common.py.

terminal_match_ok's is_fp direction-flip is the exact kind of bug this pipeline hit for real
(the OQ919513 F/R-orientation issue): silently reading the WRONG end of the primer as its 3'
end would still run without erroring, just quietly produce wrong terminal_ok values -- these
tests pin the direction-dependent behavior down so a future refactor can't reintroduce it
without a test failing.
"""
import textwrap

from mfe_common import (
    iupac_mismatch_count,
    parse_spec_best,
    reverse_complement,
    terminal_match_ok,
)


def test_reverse_complement_basic():
    assert reverse_complement("AAGG") == "CCTT"
    assert reverse_complement("acgt") == "ACGT"  # case-insensitive


def test_reverse_complement_iupac_codes():
    # each ambiguity code's own complement, order reversed
    assert reverse_complement("RYSWKM") == "KMWSRY"


def test_iupac_mismatch_count_exact_match():
    assert iupac_mismatch_count("ACGT", "ACGT") == 0


def test_iupac_mismatch_count_degenerate_match_is_not_a_mismatch():
    # R = purine (A or G); a genome base of A or G against primer R must not count
    assert iupac_mismatch_count("RCGT", "ACGT") == 0
    assert iupac_mismatch_count("RCGT", "GCGT") == 0
    # but a pyrimidine (C/T) at that position is a real mismatch
    assert iupac_mismatch_count("RCGT", "CCGT") == 1


def test_iupac_mismatch_count_length_difference_penalized():
    # 3 chars compared (all match) + abs(4-3)=1 length penalty
    assert iupac_mismatch_count("ACGT", "ACG") == 1


def test_terminal_match_ok_fp_checks_last_n_chars():
    # is_fp=True: query's OWN last n characters are its 3' end.
    query = "AAAAACCCCC"
    target_ok = "TTTTTCCCCC"  # differs everywhere except the last 3 -> 3' end intact
    target_bad = "AAAAACCCGG"  # 3' end (last 3) mismatches
    assert terminal_match_ok(query, target_ok, is_fp=True, n=3) is True
    assert terminal_match_ok(query, target_bad, is_fp=True, n=3) is False


def test_terminal_match_ok_rp_checks_first_n_chars():
    # is_fp=False: query has already been reverse_complement()'d by the caller, which puts
    # the oligo's own 3'-end base at index 0 -- so the 3' end here is the FIRST n chars, not
    # the last. This is the exact direction that was wrong before the OQ919513 fix.
    query = "CCCCCAAAAA"
    target_ok = "CCCGGAAAAA"  # first 3 chars (the 3' end) match; rest of the string differs
    target_bad = "GGCCCAAAAA"  # first 3 chars mismatch
    assert terminal_match_ok(query, target_ok, is_fp=False, n=3) is True
    assert terminal_match_ok(query, target_bad, is_fp=False, n=3) is False


def test_terminal_match_ok_n_zero_always_true():
    assert terminal_match_ok("AAAA", "TTTT", is_fp=True, n=0) is True
    assert terminal_match_ok("AAAA", "TTTT", is_fp=False, n=0) is True


def test_terminal_match_ok_n_larger_than_length_is_clamped():
    # n=10 against 4-char strings should compare the whole string, not crash/over-slice
    assert terminal_match_ok("ACGT", "ACGT", is_fp=True, n=10) is True
    assert terminal_match_ok("ACGT", "ACGA", is_fp=True, n=10) is False


def _write_spec(tmp_path, rows_text):
    p = tmp_path / "spec.tsv"
    p.write_text(textwrap.dedent(rows_text), encoding="utf-8")
    return p


SPEC_HEADER = (
    "#1-based coordinate\n"
    "#name\tchrom\tampStart\tampEnd\tampGC\tampSize\tfpName\tfpStart\tfpEnd\tfpSeq\t"
    "fpTm\tfpGC\tfpDg\trpName\trpEnd\trpStart\trpSeq\trpTm\trpGC\trpDg\tnote\n"
)


def test_parse_spec_best_picks_lowest_total_mismatch_row(tmp_path):
    # Two candidate amplicons on the same chrom for the same F/R pair; the second (perfect
    # match) has a lower total mismatch and must win even though it appears later.
    f_id, f_seq = "F1", "ACGTACGTAC"
    r_id, r_seq = "R1", "TTTTGGGGCC"
    r_seq_matched_on_genome = reverse_complement(r_seq)  # what mfeprimer reports as rpSeq
    rows = SPEC_HEADER + (
        f"Amp_1\tCHR1\t1\t100\t50\t100\t{f_id}.1\t1\t10\tACGTACGTAA\t60\t50\t-20\t"
        f"{r_id}.1\t90\t100\t{r_seq_matched_on_genome}\t60\t50\t-20\tnote\n"
        f"Amp_2\tCHR1\t1\t100\t50\t100\t{f_id}.2\t1\t10\t{f_seq}\t60\t50\t-20\t"
        f"{r_id}.2\t90\t100\t{r_seq_matched_on_genome}\t60\t50\t-20\tnote\n"
    )
    spec_path = _write_spec(tmp_path, rows)
    best = parse_spec_best(spec_path, f_id, f_seq, r_id, r_seq)
    assert "CHR1" in best
    assert best["CHR1"]["f_mismatch"] == 0
    assert best["CHR1"]["r_mismatch"] == 0
    assert best["CHR1"]["total"] == 0


def test_parse_spec_best_missing_file_returns_empty(tmp_path):
    assert parse_spec_best(tmp_path / "does_not_exist.tsv", "F1", "ACGT", "R1", "ACGT") == {}


def test_parse_spec_best_three_prime_check_can_fail_despite_low_total_mismatch(tmp_path):
    # A single internal mismatch (not at the 3' end) should still pass the 3' check; a
    # mismatch placed at F's own 3' end must fail f_3p_ok even though total mismatch is low.
    f_id, f_seq = "F1", "AAAAAAAAAT"  # 3' end is the final "T"
    r_id, r_seq = "R1", "CCCCCCCCCC"
    genome_f_seq_internal_mm = "AAAAATAAAT"  # mismatch in the middle, 3' end "T" intact
    genome_f_seq_terminal_mm = "AAAAAAAAAG"  # mismatch AT the 3' end
    r_matched = reverse_complement(r_seq)

    rows_internal = SPEC_HEADER + (
        f"Amp_1\tCHR1\t1\t100\t50\t100\t{f_id}.1\t1\t10\t{genome_f_seq_internal_mm}\t60\t50\t-20\t"
        f"{r_id}.1\t90\t100\t{r_matched}\t60\t50\t-20\tnote\n"
    )
    best_internal = parse_spec_best(_write_spec(tmp_path, rows_internal), f_id, f_seq, r_id, r_seq, three_prime_n=2)
    assert best_internal["CHR1"]["f_mismatch"] == 1
    assert best_internal["CHR1"]["f_3p_ok"] is True

    rows_terminal = SPEC_HEADER + (
        f"Amp_1\tCHR1\t1\t100\t50\t100\t{f_id}.1\t1\t10\t{genome_f_seq_terminal_mm}\t60\t50\t-20\t"
        f"{r_id}.1\t90\t100\t{r_matched}\t60\t50\t-20\tnote\n"
    )
    spec_path2 = tmp_path / "spec2.tsv"
    spec_path2.write_text(rows_terminal, encoding="utf-8")
    best_terminal = parse_spec_best(spec_path2, f_id, f_seq, r_id, r_seq, three_prime_n=2)
    assert best_terminal["CHR1"]["f_mismatch"] == 1
    assert best_terminal["CHR1"]["f_3p_ok"] is False
