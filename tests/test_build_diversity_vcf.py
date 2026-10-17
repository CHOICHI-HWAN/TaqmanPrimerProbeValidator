"""Regression tests for scripts/build_diversity_vcf.py.

`mafft --adjustdirection --add` tags every ADDED sequence it had to reverse-complement with an
"_R_" header prefix and returns that row ALREADY flipped into the reference's coordinates, so
resolve_add_flips must only strip the tag. It once reverse-complemented the row a second time,
which un-aligned it (a spot-checked reversed genome fell from 77.7% to 26.5% identity to the
reference) and silently fed random sequence into pixy for ~94% of S and 43 L genomes.
"""
from build_diversity_vcf import build_vcf_rows, resolve_add_flips, sample_qc


def test_resolve_add_flips_only_strips_the_tag_and_never_flips_the_row_again():
    aligned = {"REF": "ACGTACGTAC", "_R_S1": "ACGTTCGTAC", "S2": "GGGGACGTAC"}
    out, reversed_ids = resolve_add_flips(aligned)
    assert set(out) == {"REF", "S1", "S2"}
    assert out["S1"] == "ACGTTCGTAC"      # exactly as mafft aligned it -- NOT reverse-complemented
    assert out["S2"] == "GGGGACGTAC"
    assert out["REF"] == "ACGTACGTAC"
    assert reversed_ids == ["S1"]


def test_sample_qc_coverage_and_identity():
    aligned = {
        "REF": "ACGTACGTAC",
        "S1": "ACGTACGTAC",   # perfect
        "S2": "ACG-ACNTAC",   # gap + N -> 8/10 covered, all covered bases match
        "S3": "TGCATGCATG",   # fully covered but everything differs
    }
    qc = sample_qc(aligned, ["S1", "S2", "S3"], "REF")
    assert qc["S1"] == (1.0, 1.0)
    assert qc["S2"] == (0.8, 1.0)
    assert qc["S3"][0] == 1.0 and qc["S3"][1] < 0.2


def test_build_vcf_rows_skips_reference_gap_columns():
    # column 2 (0-indexed) has a gap in the reference -- a pure insertion in S1 relative to
    # the reference -- and must not consume a reference-coordinate position.
    aligned = {
        "REF": "AC-GT",
        "S1": "ACAGT",
    }
    rows = list(build_vcf_rows(aligned, ["REF", "S1"], "REF"))
    positions = [r[0] for r in rows]
    assert positions == [1, 2, 3, 4]  # 4 real reference bases, not 5 alignment columns


def test_build_vcf_rows_ref_allele_and_alt_alleles():
    aligned = {
        "REF": "ACGT",
        "S1": "ACGT",   # matches ref everywhere
        "S2": "ACTT",   # differs at position 3 (G->T)
    }
    rows = list(build_vcf_rows(aligned, ["REF", "S1", "S2"], "REF"))
    pos3 = next(r for r in rows if r[0] == 3)
    ref_pos, ref_base, alt_alleles, genotypes = pos3
    assert ref_base == "G"
    assert alt_alleles == ["T"]
    assert genotypes["REF"] == 0
    assert genotypes["S1"] == 0
    assert genotypes["S2"] == 1


def test_build_vcf_rows_missing_genotype_for_gap_or_ambiguous():
    aligned = {
        "REF": "ACGT",
        "S1": "AC-T",   # gap at position 3
        "S2": "ACNT",   # ambiguous base at position 3
    }
    rows = list(build_vcf_rows(aligned, ["REF", "S1", "S2"], "REF"))
    pos3 = next(r for r in rows if r[0] == 3)
    _, _, _, genotypes = pos3
    assert genotypes["S1"] is None
    assert genotypes["S2"] is None
