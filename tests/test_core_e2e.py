"""End-to-end check of the core analysis on synthetic genomes whose true answer is known by construction.

Runs the real mfeprimer (skipped if it is not configured), the parsing + genome recount, the whole-genome scan, the
engine choice, the probe search and the orientation step, and compares the verdict for every genome with the truth.
"""
import random
import subprocess
from pathlib import Path

import pytest

from mfe_common import (check_mfeprimer, choose_fr, choose_probe, find_mfeprimer, parse_spec_best,
                        probe_inside_amplicon, reverse_complement, scan_fr, scan_probe)
from orientation import kmer_set, orient

ROOT = Path(__file__).resolve().parent.parent
try:
    MFE = find_mfeprimer("", ROOT)
except SystemExit:
    MFE = ""
pytestmark = pytest.mark.skipif(not MFE or not check_mfeprimer(MFE)["ok"], reason="mfeprimer not available")

rnd = random.Random(11)
dna = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
F, P, R = dna(22), dna(26), dna(24)                    # primers/probe (R is written 5'->3' as ordered)
LEFT, GAP1, GAP2, RIGHT = dna(300), dna(40), dna(40), dna(300)


def sub(seq, pos, new=None):
    """substitute one base at pos (0-based) with a different one"""
    other = new or next(c for c in "ACGT" if c != seq[pos])
    return seq[:pos] + other + seq[pos + 1:]


def genome(f=F, p=P, r=R, p_inside=True):
    r_site = reverse_complement(r)
    probe_part = p if p_inside else dna(len(p))
    g = LEFT + f + GAP1 + probe_part + GAP2 + r_site + RIGHT
    return g if p_inside else g + dna(100) + p + dna(100)      # probe placed outside the amplicon instead


GENOMES = {
    "perfect": genome(),
    "r_5prime_mm": genome(r=sub(R, 2)),                         # 1 mismatch near R's 5' end (invisible to mfeprimer's k-mer)
    "f_two_mm": genome(f=sub(sub(F, 3), 7)),                    # 2 mismatches in F, far from its 3' end
    "f_3prime_mm": genome(f=sub(F, len(F) - 1)),                # F's last (3') base differs
    "probe_mm": genome(p=sub(P, 10)),                           # probe has 1 mismatch
    "probe_outside": genome(p_inside=False),                    # probe only binds outside the amplicon
    "unrelated": dna(1100),
}
# truth: (f mismatches, r mismatches, probe mismatches or None = no probe site inside the amplicon, validated)
TRUTH = {
    "perfect": (0, 0, 0, True), "r_5prime_mm": (0, 1, 0, True), "f_two_mm": (2, 0, 0, False),
    "f_3prime_mm": (1, 0, 0, False), "probe_mm": (0, 0, 1, False), "probe_outside": (0, 0, None, False),     # no probe site inside the amplicon at all
}


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    d = tmp_path_factory.mktemp("e2e")
    fa = d / "db.fasta"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in GENOMES.items()))
    (d / "fr.fa").write_text(f">T_F\n{F}\n>T_R\n{R}\n")
    subprocess.run([MFE, "index", "-i", str(fa)], check=True, capture_output=True)
    out = d / "mm2.txt"
    subprocess.run([MFE, "-d", str(fa), "-i", str(d / "fr.fa"), "-o", str(out), "--fasta", "--virus", "--misMatch", "2",
                    "--cpu", "2"], check=True, capture_output=True)
    mfe = parse_spec_best(Path(str(out) + ".spec.tsv"), "T_F", F, "T_R", R, three_prime_n=2, genome_seqs=GENOMES)
    res = {}
    for gid, g in GENOMES.items():
        fr, source = choose_fr(mfe.get(gid), scan_fr(F, R, g, three_prime_n=2))
        if fr is None:
            res[gid] = None
            continue
        a, b = int(fr["amp_start"]), int(fr["amp_end"])
        probe, _ = choose_probe(None, scan_probe(P, g, a, b, max_mm=8))
        validated = (fr["f_mismatch"] <= 1 and fr["r_mismatch"] <= 1 and fr["f_3p_ok"] and fr["r_3p_ok"]
                     and probe is not None and probe["mismatch"] == 0)
        res[gid] = (fr["f_mismatch"], fr["r_mismatch"], probe["mismatch"] if probe else None, bool(validated))
    return res, mfe


@pytest.mark.parametrize("gid", list(TRUTH))
def test_mismatch_counts_and_verdict_match_the_truth(results, gid):
    res, _ = results
    assert res[gid] == TRUTH[gid], (gid, res[gid])


def test_unrelated_genome_has_no_pair(results):
    assert results[0]["unrelated"] is None


def test_mfeprimer_alone_is_blind_to_the_5prime_mismatch_which_the_genome_recount_catches(results):
    _, mfe = results
    # parse_spec_best was given the genome, so it counts the real mismatch; mfeprimer's own reported text would say 0
    assert mfe["r_5prime_mm"]["r_mismatch"] == 1


def test_probe_outside_the_amplicon_does_not_count():
    g = GENOMES["probe_outside"]
    fr = scan_fr(F, R, g, three_prime_n=2)
    row = {"amplicon_start": fr["amp_start"], "amplicon_end": fr["amp_end"]}
    anywhere = scan_probe(P, g)                                   # exists in the genome ...
    assert anywhere["mismatch"] == 0
    row.update(probe_start=anywhere["start"], probe_end=anywhere["end"])
    assert not probe_inside_amplicon(row)                         # ... but outside the amplicon
    assert scan_probe(P, g, int(fr["amp_start"]), int(fr["amp_end"])) is None or \
        scan_probe(P, g, int(fr["amp_start"]), int(fr["amp_end"]))["mismatch"] > 0


def test_a_reverse_complemented_genome_is_brought_into_the_reference_orientation_and_gives_the_same_answer():
    ref_kmers = kmer_set(GENOMES["perfect"])
    flipped_input = reverse_complement(GENOMES["r_5prime_mm"])
    seq, status, *_ = orient(flipped_input, ref_kmers)
    assert status == "flipped" and seq == GENOMES["r_5prime_mm"]
    fr = scan_fr(F, R, seq, three_prime_n=2)
    assert (fr["f_mismatch"], fr["r_mismatch"]) == (0, 1)
    fr2 = scan_fr(F, R, flipped_input, three_prime_n=2)           # the engines are strand-agnostic: same counts unflipped
    assert (fr2["f_mismatch"], fr2["r_mismatch"]) == (0, 1)
