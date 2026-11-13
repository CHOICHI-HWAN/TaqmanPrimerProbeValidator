"""Per product (genome, start, end, primer sequences): the set pairs (F's set, R's set) each tool reports, as the
pipeline's group_cross_binding.py reads them (set = id before '::'). usage: compare_group_labels.py <mfe out> <pe out>"""
import collections
import sys

import compare_spec as CS


def load(p):
    d = collections.defaultdict(set)
    for r in CS.rows(p).values():
        k = (r["chrom"], r["ampStart"], r["ampEnd"], r["fpSeq"], r["rpSeq"])
        d[k].add((r["fpName"].split("::")[0], r["rpName"].split("::")[0]))
    return d


A, B = load(sys.argv[1]), load(sys.argv[2])
same = sum(1 for k in A if A[k] == B.get(k))
sup = sum(1 for k in A if k in B and A[k] < B[k])
other = sum(1 for k in A if k in B and not A[k] <= B[k])
genomes = lambda d: {k[0] for k in d}
pairs = lambda d: collections.Counter(p for v in d.values() for p in v)
print(f"products {len(A)}/{len(B)}  same set pairs {same}  engine adds pairs {sup}  mfeprimer pair missing in engine {other}"
      f"  genomes {len(genomes(A))}/{len(genomes(B))}")
print("  set pairs mfe:", dict(sorted(pairs(A).items())))
print("  set pairs pe: ", dict(sorted(pairs(B).items())))
