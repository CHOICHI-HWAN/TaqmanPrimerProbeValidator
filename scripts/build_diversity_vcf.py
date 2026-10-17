#!/usr/bin/env python3
"""Build an all-sites (variant + invariant) VCF from a whole-population multiple-sequence
alignment, for pixy (https://github.com/ksamuk/pixy) to compute nucleotide diversity (pi)
per lineage across the S/L segment.

Genome-wide population diversity needs every population genome co-aligned to the reference
once (not just within short primer amplicons like the rest of this pipeline) -- and,
critically, needs every genome's own 5'->3' orientation reconciled first: LASV S/L records
are ambisense and, exactly like the individual amplicon-orientation bug found and fixed
earlier in this pipeline (see msa_common.py / highlight_amplicon()'s "F/R orientation"
handling), a real fraction of GenBank submissions for otherwise-comparable genomes are
deposited in the opposite strand convention from the rest.

An EARLIER version of this script aligned all population genomes together in one joint
mafft MSA (`mafft --adjustdirection --auto`) and read VCF sites straight off that alignment's
columns. That produced badly inflated missingness: spot-checked against the population's own
raw FASTA (zero ambiguous/N bases in the source genomes), the joint-MSA-derived VCF showed
~35-47% "missing" genotype calls for nearly every non-reference-like genome -- a joint
progressive alignment of 662 divergent, multi-lineage sequences fragments badly for anything
that isn't close to the guide tree's early branches, gapping out real, alignable sequence
rather than genuinely failing to align it. Re-aligning the SAME genomes one-at-a-time against
just the reference (pairwise, or via mafft's own `--keeplength --add` reference-anchored
mode) dropped that to ~1-12% for the same test genomes -- confirmed empirically during
development. So this script now builds one reference-anchored alignment: every population
genome is added, in a single `--add` call, to a fixed one-sequence profile containing only
the reference (`--keeplength` guarantees the output alignment stays exactly reference-length,
so every other genome's own insertions relative to the reference are dropped rather than
shifting anyone's column numbering -- exactly what a reference-coordinate VCF needs anyway).
Each genome's alignment quality is therefore independent of every other genome's, unlike the
old joint MSA where one divergent sequence could distort the whole guide tree. `mafft`'s own
`--adjustdirection` flag -- previously avoided in this pipeline for ~20nt primer alignments
(too little signal for its heuristic to work reliably there) -- IS reliable at whole-genome
length, and combined with `--add` only ever flips the ADDED sequence (the reference profile
itself, being the fixed base, is never touched), so there is no separate orientation
reconciliation step to run afterward.

pixy needs a real VCF (bgzipped + tabix-indexed) that records EVERY site, not just variant
ones -- that's the entire point of pixy over naive per-SNP diversity calculators, since
correctly normalizing pi requires knowing the true number of comparable sites per window.
Rather than a read-aligner + variant-caller chain (minimap2/bcftools), which -- tested during
development -- struggles badly with LASV's cross-lineage divergence (many genomes simply
fail to seed an alignment against a single reference at usable identity), this script builds
the VCF directly from the reference-anchored mafft alignment: each column becomes one VCF
record, with the reference row as the REF allele and every other distinct base at that column
as an ALT allele; gap/N/IUPAC-ambiguous bases become missing genotypes rather than forcing a
called allele.
"""
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import find_mafft, reverse_complement  # noqa: E402

IUPAC_VALID = set("ACGT")


def find_tool(name, dedicated_path=None):
    found = shutil.which(name)
    if found:
        return found
    if dedicated_path and Path(dedicated_path).exists():
        return dedicated_path
    for cand in sorted(Path("/BI/conda/envs").glob(f"*/bin/{name}")):
        return str(cand)
    sys.exit(f"[error] {name} binary not found on PATH or under /BI/conda/envs/*/bin/{name}")


def read_fasta(path):
    seqs = {}
    header = None
    chunks = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs[header] = "".join(chunks)
    return seqs


def run_mafft_reference_anchored(mafft_bin, ref_fasta: Path, add_fasta: Path, out_fasta: Path, threads=8):
    """Aligns every sequence in add_fasta to the FIXED single-sequence profile ref_fasta,
    one at a time (mafft's own --add machinery), rather than building one joint MSA across
    the whole population -- see the module docstring for why the joint approach badly
    fragments divergent genomes. --keeplength forces the output alignment width to stay
    exactly the reference's own length: any insertion a genome has relative to the reference
    is dropped rather than adding a column (which would shift everyone else's numbering) --
    exactly right for a reference-coordinate VCF, which has no use for a position that isn't
    in the reference anyway. --adjustdirection can only flip sequences being added, never the
    fixed reference profile."""
    cmd = [mafft_bin, "--keeplength", "--adjustdirection", "--add", str(add_fasta),
           "--quiet", "--thread", str(threads), str(ref_fasta)]
    with open(out_fasta, "w", encoding="utf-8") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] mafft failed: {result.stderr}")


def resolve_add_flips(aligned: dict):
    """Returns (rows, reversed_ids): `mafft --adjustdirection` tags every ADDED sequence it had
    to reverse-complement with an "_R_" header prefix, and the aligned row it hands back is
    ALREADY the flipped sequence laid out in the reference's own coordinates -- so all that's
    left to do here is strip the tag. (An earlier version reverse-complemented the row a second
    time "to restore the original orientation", which un-aligned it: a spot-checked reversed
    genome went from 77.7% identity to the reference back down to 26.5% -- random -- and every
    reversed genome, i.e. ~94% of S and 43 of 366 L, silently fed garbage into pixy.)"""
    out = {}
    reversed_ids = []
    for header, seq in aligned.items():
        if header.startswith("_R_"):
            out[header[3:]] = seq
            reversed_ids.append(header[3:])
        else:
            out[header] = seq
    return out, reversed_ids


def sample_qc(aligned: dict, sample_ids, ref_id: str):
    """Per-genome (coverage, identity) against the reference row of a reference-anchored
    alignment: coverage = fraction of reference positions where this genome has a real A/C/G/T
    call; identity = fraction of those called positions that equal the reference base. Low
    identity means a misaligned/foreign sequence, low coverage a truncated one -- both drag
    the per-window pi estimate toward noise."""
    ref_row = aligned[ref_id].upper()
    ref_cols = [i for i, c in enumerate(ref_row) if c in IUPAC_VALID]
    out = {}
    for sid in sample_ids:
        row = aligned[sid].upper()
        called = same = 0
        for i in ref_cols:
            c = row[i] if i < len(row) else "-"
            if c in IUPAC_VALID:
                called += 1
                same += (c == ref_row[i])
        out[sid] = (called / len(ref_cols) if ref_cols else 0.0, same / called if called else 0.0)
    return out


def build_vcf_rows(aligned: dict, sample_order: list, ref_id: str):
    """Yields (ref_pos, ref_allele, alt_alleles, genotypes) per alignment column that has a
    real (non-gap, non-ambiguous) base in the reference row. --keeplength guarantees the
    reference row itself is gap-free and exactly genome-length, so ref_pos normally just
    tracks col+1 -- the gap check stays only as a defensive fallback."""
    n_cols = len(aligned[ref_id])
    ref_row = aligned[ref_id]
    ref_pos = 0
    for col in range(n_cols):
        ref_base = ref_row[col].upper() if col < len(ref_row) else "-"
        if ref_base == "-":
            continue
        ref_pos += 1
        if ref_base not in IUPAC_VALID:
            continue
        alt_alleles = []
        genotypes = {}
        for sid in sample_order:
            row = aligned[sid]
            c = row[col].upper() if col < len(row) else "-"
            if c not in IUPAC_VALID:
                genotypes[sid] = None
                continue
            if c == ref_base:
                genotypes[sid] = 0
            else:
                if c not in alt_alleles:
                    alt_alleles.append(c)
                genotypes[sid] = 1 + alt_alleles.index(c)
        yield ref_pos, ref_base, alt_alleles, genotypes


def write_vcf(rows, sample_order, chrom, chrom_len, out_path: Path):
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={chrom},length={chrom_len}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_order) + "\n")
        for pos, ref_allele, alt_alleles, genotypes in rows:
            alt_str = ",".join(alt_alleles) if alt_alleles else "."
            gts = ["." if genotypes.get(sid) is None else str(genotypes[sid]) for sid in sample_order]
            fh.write(f"{chrom}\t{pos}\t.\t{ref_allele}\t{alt_str}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")


def bgzip_and_index(vcf_path: Path, bgzip_bin, tabix_bin):
    subprocess.run([bgzip_bin, "-f", str(vcf_path)], check=True)
    gz_path = Path(str(vcf_path) + ".gz")
    subprocess.run([tabix_bin, "-f", "-p", "vcf", str(gz_path)], check=True)
    return gz_path


def write_populations_file(sample_order, meta_by_id, out_path: Path):
    with open(out_path, "w", encoding="utf-8") as fh:
        for sid in sample_order:
            lineage = meta_by_id.get(sid, {}).get("lineage") or "NA"
            fh.write(f"{sid}\t{lineage}\n")


def main():
    import csv

    seg = snakemake.wildcards.seg
    pop_fasta = Path(snakemake.input.population)
    metadata_path = Path(snakemake.input.metadata)
    ref_fasta_path = Path(snakemake.input.ref_fasta)
    ref_accession = snakemake.params.ref_accession
    threads = int(snakemake.threads) if hasattr(snakemake, "threads") and snakemake.threads else 8

    out_vcf_gz = Path(snakemake.output.vcf)
    out_pops = Path(snakemake.output.populations)
    out_vcf_gz.parent.mkdir(parents=True, exist_ok=True)

    mafft_bin = find_mafft()
    bgzip_bin = find_tool("bgzip")
    tabix_bin = find_tool("tabix")

    # The reference accession is itself a normal member of the population fasta -- excluded
    # here from the "add" side since it's already the fixed base profile (ref_fasta_path);
    # aligning it against itself too would just be a wasted (0-gap) --add call.
    pop_seqs = read_fasta(pop_fasta)
    add_fasta = out_vcf_gz.parent / f"{seg}_add.fasta"
    with open(add_fasta, "w", encoding="utf-8") as fh:
        for gid, seq in pop_seqs.items():
            if gid == ref_accession:
                continue
            fh.write(f">{gid}\n{seq}\n")

    aligned_fasta = out_vcf_gz.parent / f"{seg}_aligned.fasta"
    print(f"[diversity] aligning {len(pop_seqs) - 1} genomes to reference {ref_accession} "
          f"with mafft --keeplength --add ...", file=sys.stderr)
    run_mafft_reference_anchored(mafft_bin, ref_fasta_path, add_fasta, aligned_fasta, threads=threads)

    aligned = read_fasta(aligned_fasta)
    aligned, reversed_ids = resolve_add_flips(aligned)
    min_cov = float(snakemake.params.min_coverage)
    min_ident = float(snakemake.params.min_identity)
    qc = sample_qc(aligned, [g for g in aligned if g != ref_accession], ref_accession)
    excluded = {"low_coverage": [], "low_identity": []}
    kept = [ref_accession]
    for sid, (cov, ident) in sorted(qc.items()):
        if cov < min_cov:
            excluded["low_coverage"].append(sid)
        elif ident < min_ident:
            excluded["low_identity"].append(sid)
        else:
            kept.append(sid)
    aligned = {g: aligned[g] for g in kept}
    sample_order = sorted(aligned.keys())
    print(f"[diversity] {len(reversed_ids)} genomes were reverse-complemented by mafft; QC kept "
          f"{len(sample_order)}/{len(qc) + 1} (coverage>={min_cov}, identity>={min_ident}; dropped "
          f"{len(excluded['low_coverage'])} low-coverage, {len(excluded['low_identity'])} low-identity); "
          f"alignment width {len(aligned[ref_accession])}", file=sys.stderr)
    Path(snakemake.output.qc).write_text(json.dumps({
        "input_genomes": len(qc) + 1, "kept": len(sample_order), "reversed_fixed": len(reversed_ids),
        "min_coverage": min_cov, "min_identity": min_ident,
        "excluded_low_coverage": excluded["low_coverage"], "excluded_low_identity": excluded["low_identity"],
        "per_genome": {g: {"coverage": round(c, 4), "identity": round(i, 4)} for g, (c, i) in qc.items()},
    }, indent=1), encoding="utf-8")

    chrom = f"{snakemake.params.chrom_prefix}_{seg.upper()}"
    ref_len = len(aligned[ref_accession].replace("-", ""))
    rows = list(build_vcf_rows(aligned, sample_order, ref_accession))
    print(f"[diversity] {len(rows)} reference-coordinate sites", file=sys.stderr)

    raw_vcf = out_vcf_gz.with_suffix("")  # e.g. .../allsites.vcf (bgzip appends .gz itself)
    write_vcf(rows, sample_order, chrom, ref_len, raw_vcf)
    bgzip_and_index(raw_vcf, bgzip_bin, tabix_bin)

    meta_by_id = {}
    with open(metadata_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta_by_id[r["accession"]] = r
    write_populations_file(sample_order, meta_by_id, out_pops)

    aligned_fasta.unlink(missing_ok=True)
    add_fasta.unlink(missing_ok=True)
    print(f"[diversity] wrote {out_vcf_gz} and {out_pops}", file=sys.stderr)


if __name__ == "__main__":
    main()
