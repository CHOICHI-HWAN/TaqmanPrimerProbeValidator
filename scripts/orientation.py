"""Bring every genome into the reference genome's orientation before any analysis.

Public databases hold some genomes on the opposite strand. Matching works on either strand, but coordinates, amplicon
sequences, the MSA and gene positions are only comparable when all genomes read the same way, so each genome is
reverse-complemented when it shares clearly more k-mers with the reference's reverse complement than with the reference.
"""
from mfe_common import reverse_complement

K = 11
MIN_SUPPORT = 20          # shared k-mers needed before an orientation call is trusted
SAMPLE_ABOVE = 60000      # a genome longer than this (poxviruses, ~197 kb) is judged on its central stretch only: hundreds of thousands of k-mers per
SAMPLE_LEN = 20000        # genome take hours for thousands of genomes, and 20 kb of unique sequence is far more than the call needs


def kmer_set(seq: str, k: int = K) -> set:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def orientation_scores(seq: str, ref_kmers: set, k: int = K):
    """(shared k-mers as given, shared k-mers after reverse-complementing) -- both measured against the reference."""
    if len(seq) > SAMPLE_ABOVE:                       # (the ends of such genomes are inverted repeats, so the middle is the informative part)
        mid = len(seq) // 2
        seq = seq[mid - SAMPLE_LEN // 2: mid + SAMPLE_LEN // 2]
    forward = len(kmer_set(seq, k) & ref_kmers)
    reverse = len(kmer_set(reverse_complement(seq), k) & ref_kmers)
    return forward, reverse


def orient(seq: str, ref_kmers: set, k: int = K):
    """-> (sequence in reference orientation, status, forward_score, reverse_score);
    status: "ok" (kept), "flipped" (reverse-complemented) or "undetermined" (too little shared sequence; kept as is)."""
    f, r = orientation_scores(seq, ref_kmers, k)
    if max(f, r) < MIN_SUPPORT:
        return seq, "undetermined", f, r
    if r > f:
        return reverse_complement(seq), "flipped", f, r
    return seq, "ok", f, r
