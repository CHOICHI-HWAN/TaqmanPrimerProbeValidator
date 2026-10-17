import random

from mfe_common import reverse_complement
from orientation import kmer_set, orient


def genome(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def test_reversed_genomes_are_flipped_into_the_reference_orientation():
    ref = genome(3000, 1)
    ref_k = kmer_set(ref)
    divergent = list(ref)                                   # a ~8% divergent relative
    rnd = random.Random(2)
    for i in range(0, len(divergent), 13):
        divergent[i] = rnd.choice("ACGT")
    divergent = "".join(divergent)
    seq, status, f, r = orient(divergent, ref_k)
    assert status == "ok" and seq == divergent and f > r
    seq, status, f, r = orient(reverse_complement(divergent), ref_k)
    assert status == "flipped" and seq == divergent and r > f             # back to exactly the original strand
    seq, status, *_ = orient(genome(3000, 99), ref_k)                     # unrelated sequence: no call, left unchanged
    assert status == "undetermined"
