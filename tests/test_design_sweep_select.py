"""Ranking helpers of the design sweep (scripts/design_sweep_select.py): synthetic designs only, no real sequences."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import design_sweep_select as sel  # noqa: E402


def _d(cov, deg_chars="", amp=120, ftm=55.0, rtm=55.0):
    seq = "A" * 18 + deg_chars
    return {"F": seq, "P": "C" * 20, "R": "G" * 18, "est_cov_pct": cov, "amplicon_len": amp, "F_tm": ftm, "R_tm": rtm}


def test_degeneracy_counts_iupac_codes_only():
    assert sel.degeneracy("ACGT") == 0
    assert sel.degeneracy("ARCYGN") == 3


def test_coverage_wins_but_degeneracy_and_long_amplicons_cost():
    plain, degenerate, long_amp = _d(95.0), _d(95.0, "RYRYRYRY"), _d(95.0, amp=320)
    assert sel.score(plain) > sel.score(degenerate)
    assert sel.score(plain) > sel.score(long_amp)
    assert sel.score(_d(96.0, "RYRYRYRY")) < sel.score(_d(95.5))       # a little coverage does not pay for a lot of degeneracy


def test_unbalanced_primer_tm_is_penalised():
    assert sel.score(_d(95.0, ftm=55.0, rtm=62.0)) < sel.score(_d(95.0, ftm=55.0, rtm=56.0))
