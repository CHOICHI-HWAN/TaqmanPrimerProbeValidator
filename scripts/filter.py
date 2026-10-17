#!/usr/bin/env python3
"""Step 2: keep only records with no ambiguous 'N' bases and length >= segment threshold
(S >= 3000bp, L >= 6000bp), deduplicated by accession. Output header = bare accession
(required so mfeprimer's genome/"Chrom" IDs stay unique in later steps).
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orientation import kmer_set, orient  # noqa: E402


def read_fasta(path: Path):
    header, chunks = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def main():
    raw_fasta = Path(snakemake.input.fasta)
    metadata_curated = Path(snakemake.input.metadata_curated)
    min_len = int(snakemake.params.min_len)
    seg = snakemake.wildcards.seg

    # every kept genome is brought into the reference genome's orientation (see orientation.py)
    ref_acc = str(snakemake.params.reference_accession)
    ref_seq = next((s.upper() for h, s in read_fasta(raw_fasta) if h.split()[0].split("|")[0].strip() == ref_acc), None)
    ref_kmers = kmer_set(ref_seq) if ref_seq else None
    orient_rows, orient_counts = [], {"ok": 0, "flipped": 0, "undetermined": 0}

    out_fasta = Path(snakemake.output.filtered_fasta)
    out_meta = Path(snakemake.output.filtered_metadata)
    out_report = Path(snakemake.output.report)

    meta = {}
    with open(metadata_curated, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["accession"]] = row

    seen = set()
    kept = dropped_short = dropped_n = dropped_dup = dropped_nometa = 0
    lineage_counts, country_counts = {}, {}

    with open(out_fasta, "w", encoding="utf-8") as out_fa, \
         open(out_meta, "w", newline="", encoding="utf-8") as out_tsv:
        writer = csv.writer(out_tsv, delimiter="\t")
        writer.writerow(["accession", "lineage", "date", "country", "region"])
        for acc, seq in read_fasta(raw_fasta):
            acc_clean = acc.split()[0].split("|")[0].strip()
            seq_u = seq.upper()
            if len(seq_u) < min_len:
                dropped_short += 1
                continue
            if "N" in seq_u:
                dropped_n += 1
                continue
            if acc_clean in seen:
                dropped_dup += 1
                continue
            m = meta.get(acc_clean)
            if m is None:
                dropped_nometa += 1
                continue
            seen.add(acc_clean)
            if ref_kmers is None:                  # reference not in the data: the first kept genome defines the orientation
                ref_kmers = kmer_set(seq_u)
            seq_u, status, f_score, r_score = orient(seq_u, ref_kmers)
            orient_counts[status] += 1
            orient_rows.append([acc_clean, status, f_score, r_score])
            out_fa.write(f">{acc_clean}\n{seq_u}\n")
            writer.writerow([acc_clean, m["lineage"], m["date"], m["country"], m["region"]])
            kept += 1
            lineage_counts[m["lineage"]] = lineage_counts.get(m["lineage"], 0) + 1
            country_counts[m["country"]] = country_counts.get(m["country"], 0) + 1

    with open(snakemake.output.orientation, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "orientation", "kmers_as_given", "kmers_reverse_complemented"])
        w.writerows(orient_rows)

    with open(out_report, "w", encoding="utf-8") as rep:
        rep.write(f"segment={seg} min_len={min_len}\n")
        rep.write(f"orientation: reversed_to_reference={orient_counts['flipped']} undetermined={orient_counts['undetermined']}\n")
        rep.write(f"kept={kept} dropped_short(<{min_len}bp)={dropped_short} dropped_ambiguousN={dropped_n} "
                  f"dropped_dup_accession={dropped_dup} dropped_no_metadata={dropped_nometa}\n\n")
        rep.write("lineage counts:\n")
        for k, v in sorted(lineage_counts.items(), key=lambda x: -x[1]):
            rep.write(f"  {k}\t{v}\n")
        rep.write("\ncountry counts:\n")
        for k, v in sorted(country_counts.items(), key=lambda x: -x[1]):
            rep.write(f"  {k}\t{v}\n")


if __name__ == "__main__":
    main()
