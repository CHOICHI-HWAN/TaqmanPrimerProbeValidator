#!/usr/bin/env python3
"""Clade/lineage per accession from a locally built Auspice JSON (`nextstrain build`, augur clades).

Writes data/lineage/<seg>_lineage.tsv (accession, lineage). Tips are matched on their `accession`
node attribute (falling back to the node name). Genomes that are not tips of the tree are simply
absent here; rename_headers labels them with the profile's `lineage.missing` value.
"""
import csv
import json
from pathlib import Path


def tips(node, out):
    kids = node.get("children")
    if kids:
        for k in kids:
            tips(k, out)
    else:
        out.append(node)


def attr_value(attrs, name):
    v = attrs.get(name)
    return v.get("value") if isinstance(v, dict) else v


def main():
    auspice = json.loads(Path(snakemake.input.auspice).read_text(encoding="utf-8"))
    attr = snakemake.params.attr
    t = []
    tips(auspice["tree"], t)
    rows, seen = [], set()
    for n in t:
        attrs = n.get("node_attrs", {})
        acc = (attr_value(attrs, "accession") or n.get("name") or "").strip()
        val = attr_value(attrs, attr)
        if acc and val and acc not in seen:
            seen.add(acc)
            rows.append((acc, val))
    out = Path(snakemake.output.table)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "lineage"])
        w.writerows(rows)
    print(f"[lineage] {len(rows)} tips with {attr} from {snakemake.input.auspice}")


if __name__ == "__main__":
    main()
