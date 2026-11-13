"""Totals over a compare_spec.py summary file (work/cmp/summary.txt). usage: agg_spec.py <summary>"""
import collections
import re
import sys

t = open(sys.argv[1]).read()
tot, fd, worst = collections.Counter(), collections.Counter(), {}
for m in re.finditer(r"products: ref (\d+)\s+new (\d+)\s+missing (\d+)\s+extra (\d+)", t):
    a, b, mi, ex = map(int, m.groups())
    tot["products"] += a
    tot["missing"] += mi
    tot["extra"] += ex
    tot["sets with a difference"] += bool(mi or ex)
for m in re.finditer(r"field diffs: (\{.*?\})  worst: (\{.*?\})", t):
    for k, v in eval(m.group(1)).items():
        fd[k] += v
    for k, v in eval(m.group(2)).items():
        worst[k] = max(worst.get(k, 0), v)
print(dict(tot), dict(fd), worst, f"failed-genome set same {len(re.findall('failed: same set True', t))}/{t.count('### ')}")
