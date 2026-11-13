#!/usr/bin/env python3
"""Per-genome outcome of the data QC filter (scripts/filter.py): one row per raw record with
status `passed` or the reason it was dropped, so the web UI can mark which genomes made it into the
analysed population. Applies filter.py's rules in the same order (length, ambiguous N, duplicate
accession, missing metadata) and checks the result against the filter's own output -- if the two
ever disagree, this fails loudly instead of showing wrong marks.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from filter import read_fasta  # noqa: E402

STATUSES = ("passed", "dropped_short", "dropped_excluded", "dropped_ambiguousN", "dropped_dup_accession", "dropped_no_metadata")


def classify(records, curated_accessions, min_len, max_n_fraction=0.0, excluded=frozenset()):
    """records: iterable of (accession, sequence) -> list of (accession, length, status)."""
    seen, out = set(), []
    for acc, seq in records:
        acc_clean = acc.split()[0].split("|")[0].strip()
        seq_u = seq.upper()
        if len(seq_u) < min_len:
            status = "dropped_short"
        elif acc_clean in excluded:
            status = "dropped_excluded"
        elif seq_u.count("N") > max_n_fraction * len(seq_u):
            status = "dropped_ambiguousN"
        elif acc_clean in seen:
            status = "dropped_dup_accession"
        elif acc_clean not in curated_accessions:
            status = "dropped_no_metadata"
        else:
            seen.add(acc_clean)
            status = "passed"
        out.append((acc_clean, len(seq_u), status))
    return out


def main():
    min_len = int(snakemake.params.min_len)
    with open(snakemake.input.metadata_curated, newline="", encoding="utf-8") as fh:
        curated = {r["accession"] for r in csv.DictReader(fh, delimiter="\t")}
    excluded = set()
    if snakemake.params.exclude_list:
        for line in Path(snakemake.params.exclude_list).read_text(encoding="utf-8").splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                excluded.add(line)
    rows = classify(read_fasta(Path(snakemake.input.fasta)), curated, min_len, float(snakemake.params.max_n_fraction), excluded)
    with open(snakemake.input.filtered_metadata, newline="", encoding="utf-8") as fh:
        kept_by_filter = {r["accession"] for r in csv.DictReader(fh, delimiter="\t")}
    passed = {acc for acc, _, st in rows if st == "passed"}
    if passed != kept_by_filter:
        sys.exit(f"[error] filter_status disagrees with filter.py: {len(passed ^ kept_by_filter)} accessions differ")
    with open(snakemake.output.status, "w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "length", "status"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
