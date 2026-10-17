"""Regression tests for scripts/run_dimer_hairpin.py's text-report parser.

mfeprimer's `dimer`/`hairpin` subcommands emit a fixed-width human-readable report, not JSON
(the -j/--json flag was tried and found to be a no-op on the installed build) -- these tests
pin the parser to real captured output so a future mfeprimer version with slightly different
spacing/wording doesn't silently start dropping every entry (as actually happened once during
development: the hairpin report's whitespace between the Score line and the alignment block
isn't a clean blank line, which an earlier version of this parser's regex missed entirely).
"""
from seqgen import fill
from run_dimer_hairpin import (
    parse_dimer_report,
    parse_hairpin_report,
    summarize_dimers,
    summarize_hairpins,
)

# Trimmed but structurally faithful captures of real mfeprimer output (see scripts/
# run_dimer_hairpin.py's module docstring for why the format is parsed this way).
DIMER_REPORT = fill("""MFEprimer-3.0 Dimer Reports (2026-09-18 06:43:41)

Primer ID                      Sequence (5'-->3')                    Length     GC      Tm      Dg
                                                                      (bp)      (%)    (°C)  (kcal/mol)

F                               @@S0:20@@                      21   50.00   60.00    -20.00
R                               @@S1:20@@                      21   50.00   60.00    -20.00



Dimer List (2)


Dimer 1: R x F

    Score: 20, Delta G = -21.38 kcal/mol

5' @@S0:20@@ 3'
   ::::::::::::::::::::
3' @@S2:20@@ 5'


Dimer 2: R x R

    Score: 6, Delta G = -5.70 kcal/mol

     5' @@S1:20@@ 3'
        ::::::
3' @@S2:20@@ 5'



Parameters

                   Primer file: test.fasta
""")

HAIRPIN_REPORT = fill("""MFEprimer-3.0 Hairpin Reports (2026-09-18 06:44:54)

Primer ID                      Sequence (5'-->3')                    Length     GC      Tm      Dg
                                                                      (bp)      (%)    (°C)  (kcal/mol)

H1                              @@S3:19@@                       19   63.16   66.45    -25.37



Hairpin List (1)


Hairpin 1: H1

    Score: 7, Tm = 32.15 °C, Delta G = -6.86 kcal/mol

    ///////-----\\\\\\\\\\\\\\
    @@S3:19@@



Parameters

                   Primer file: hp_test.fasta
""")

HAIRPIN_REPORT_NONE = fill("""MFEprimer-3.0 Hairpin Reports (2026-09-18 06:43:41)

Primer ID                      Sequence (5'-->3')                    Length     GC      Tm      Dg
                                                                      (bp)      (%)    (°C)  (kcal/mol)

F                               @@S0:20@@                      21   50.00   60.00    -20.00



Hairpin List (0)


No hairpins found.


Parameters

                   Primer file: test.fasta
""")


def test_parse_dimer_report_extracts_role_score_and_dg():
    dimers = parse_dimer_report(DIMER_REPORT)
    assert len(dimers) == 2
    d1 = dimers[0]
    assert d1["role1"] == "R" and d1["role2"] == "F"
    assert d1["score"] == 20
    assert d1["dg"] == -21.38
    assert fill("5' @@S0:20@@ 3'") in d1["align"]


def test_parse_dimer_report_empty_list():
    text = DIMER_REPORT.replace("Dimer List (2)", "Dimer List (0)").split("Dimer 1:")[0] + "No dimers found.\n"
    assert parse_dimer_report(text) == []


def test_parse_hairpin_report_handles_non_blank_whitespace_before_alignment():
    # regression: the gap between the Score line and the alignment block is
    # "\n    \n    " (newline + indent), not a clean "\n\n" -- a naive \n\n-based split used
    # to silently produce zero parsed hairpins here.
    hairpins = parse_hairpin_report(HAIRPIN_REPORT)
    assert len(hairpins) == 1
    h = hairpins[0]
    assert h["role"] == "H1"
    assert h["score"] == 7
    assert h["tm"] == 32.15
    assert h["dg"] == -6.86


def test_parse_hairpin_report_none_found():
    assert parse_hairpin_report(HAIRPIN_REPORT_NONE) == []


def test_role_extracted_from_degenerate_expansion_suffix():
    # mfeprimer names each expanded IUPAC variant "ROLE.N" (e.g. "F.3"); the role is
    # whatever precedes the first dot, and must still work with no dot at all (no
    # degenerate bases -> mfeprimer leaves the id bare, e.g. plain "F").
    text = DIMER_REPORT.replace("R x F", "R.3 x F.12")
    dimers = parse_dimer_report(text)
    assert dimers[0]["role1"] == "R"
    assert dimers[0]["role2"] == "F"


def test_summarize_dimers_picks_worst_per_unordered_role_pair():
    dimers = parse_dimer_report(DIMER_REPORT)  # R-F at -21.38, R-R at -5.70
    summary = summarize_dimers(dimers)
    assert set(summary.keys()) == {"F-R", "R-R"}
    assert summary["F-R"]["worst"]["dg"] == -21.38
    assert summary["R-R"]["worst"]["dg"] == -5.70


def test_summarize_dimers_flags_entries_at_or_below_threshold():
    dimers = parse_dimer_report(DIMER_REPORT)
    summary = summarize_dimers(dimers)
    # DG_FLAG_THRESHOLD is -9.0: R-F (-21.38) flagged, R-R (-5.70) not
    assert summary["F-R"]["flagged"] == 1
    assert summary["R-R"]["flagged"] == 0


def test_summarize_hairpins_worst_and_flagged():
    hairpins = parse_hairpin_report(HAIRPIN_REPORT)
    summary = summarize_hairpins(hairpins)
    assert summary["H1"]["worst"]["dg"] == -6.86
    assert summary["H1"]["flagged"] == 0  # -6.86 is above (less negative than) -9.0
