#!/usr/bin/env python3
"""data/filtered/<seg>_metadata_annotated.tsv: every raw record with its clade/lineage and an explicit
QC verdict. Columns = the raw metadata columns + <lineage label> + qc_status + qc_reason (+ length_seq,
the real sequence length, which is what the length filter looks at). Records that only exist as
sequences (no metadata row) are listed too, flagged QC FAIL.
"""
import csv

REASONS = {
    "dropped_short": "길이 미달 (부분 서열, < {min_len}bp)",
    "dropped_ambiguousN": "모호염기(N) 과다",
    "dropped_dup_accession": "중복 accession",
    "dropped_no_metadata": "메타데이터 없음",
    "dropped_excluded": "제외 목록(exclude list)",
}


def main():
    with open(snakemake.input.meta, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        raw_cols = rd.fieldnames or []
        raw = {r[snakemake.params.accession_col]: r for r in rd}
    with open(snakemake.input.curated, newline="", encoding="utf-8") as fh:
        lin = {r["accession"]: r["lineage"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(snakemake.input.status, newline="", encoding="utf-8") as fh:
        status = list(csv.DictReader(fh, delimiter="\t"))
    lab = snakemake.params.lineage_label
    cols = list(raw_cols) + [lab, "qc_status", "qc_reason", "length_seq"]
    acc_col = snakemake.params.accession_col
    with open(snakemake.output.annotated, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for st in status:
            acc = st["accession"]
            row = dict(raw.get(acc, {acc_col: acc}))
            row[lab] = lin.get(acc, "")
            ok = st["status"] == "passed"
            row["qc_status"] = "QC PASS" if ok else "QC FAIL"
            row["qc_reason"] = "" if ok else REASONS.get(st["status"], st["status"]).format(min_len=snakemake.params.min_len)
            row["length_seq"] = st["length"]
            w.writerow(row)


if __name__ == "__main__":
    main()
