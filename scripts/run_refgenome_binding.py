#!/usr/bin/env python3
"""For one assay set, locate F/R/Probe binding sites individually on the single reference
genome (AY628203 for S, AY628202 for L) -- mirrors the '참조유전체비교' primer-binding-map
slides.

The reference genome is itself just a regular member of the population database, so its F/R
mismatch/position is usually already sitting, already validated, in this same set's own
population-wide results (primer_qc's --virus spec.tsv for F/R, the alignment table's probe
search for Probe) -- reusing those is both more accurate and cheaper than re-deriving
mismatch counts from a separate, narrower single-genome search: a plain fixed-window
substring comparison (what this script's own independent `mfeprimer search` and its
direct_best_match() fallback both are) can only count substitutions, not indels, and mfeprimer's
--virus geometric search has already been seen to find a real site that a same-primer,
same-genome `search -m <wide>` call entirely misses (confirmed concretely for
VendorB_S_NG_1/VendorB_S_SL_2's F/R on AY628203: population mismatch 1 each, independent search -m 6
finding nothing at all). Only when the reference genome ISN'T already in that population-wide
result (mismatch too high for the mismatch_levels this pipeline otherwise searches, e.g. >2 for
F/R) does this script fall back to its own independent, wider-tolerance search, and only when
even that fails does it fall back further to the direct brute-force scan."""
import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import (  # noqa: E402
    find_mfeprimer, write_oligo_fasta, load_primers_manifest, direct_best_match, read_fasta_seqs,
    parse_spec_best, iupac_mismatch_count, reverse_complement,
)


def search_best_hit(mfeprimer_bin, ref_fasta, oligo_fa, out_txt, max_mismatch):
    if out_txt.exists():
        out_txt.unlink()
    cmd = [mfeprimer_bin, "search", str(ref_fasta), str(oligo_fa), str(out_txt), "-m", str(max_mismatch)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] mfeprimer search failed ({oligo_fa}):\n{result.stderr}")
    best = None
    with open(out_txt, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, _qname, _size, strand, mismatch = line.rstrip("\n").split("\t")[:7]
            mismatch = int(mismatch)
            if best is None or mismatch < best["mismatches"]:
                best = {"chrom": chrom, "start": int(start), "end": int(end), "strand": strand, "mismatches": mismatch}
    return best


def load_probe_row(alignment_path: Path, ref_chrom: str, probe_seq: str, ref_genome_seq: str):
    """The reference genome's own row from this set's already-computed population alignment
    table (which already has its own direct_best_match() fallback baked in for the probe --
    see build_alignment_table.py), if present and the probe actually matched there. The
    table itself doesn't track strand (the geometric amplicon search resolves F/R roles
    without needing to expose it), so it's re-derived here the cheap way -- the position and
    mismatch count are already known-correct, so whichever orientation reproduces that same
    mismatch count against the genome is the real one."""
    if not alignment_path.exists():
        return None
    with open(alignment_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["genome_id"] != ref_chrom:
                continue
            if not r.get("probe_start") or not r.get("p_mismatch"):
                return None  # row exists but probe wasn't found for this genome either
            start, end = int(r["probe_start"]), int(r["probe_end"])
            window = ref_genome_seq[start - 1:end]  # 1-based inclusive -> 0-based slice
            mm_plus = iupac_mismatch_count(probe_seq, window)
            mm_minus = iupac_mismatch_count(probe_seq, reverse_complement(window))
            strand = "+" if mm_plus <= mm_minus else "-"
            return {"start": start, "end": end, "mismatches": int(r["p_mismatch"]), "strand": strand}
    return None


def main():
    set_id = snakemake.wildcards.set_id
    ref_fasta = Path(snakemake.input.ref)
    max_mismatch = int(snakemake.params.max_mismatch_search)
    manifest = load_primers_manifest(Path(snakemake.input.manifest))
    row = next(r for r in manifest if r["set_id"] == set_id)

    pipeline_root = Path(snakemake.scriptdir).parent
    mfeprimer_bin = find_mfeprimer(snakemake.config.get("mfeprimer_bin", ""), pipeline_root)

    out_dir = Path(snakemake.output.binding).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_seqs = read_fasta_seqs(ref_fasta)  # only one record, but keyed by id for the fallbacks below
    ref_chrom = next(iter(ref_seqs)) if ref_seqs else None

    fr_validated = parse_spec_best(Path(snakemake.input.spec), row["F_id"], row["F_seq"],
                                    row["R_id"], row["R_seq"], genome_seqs=ref_seqs).get(ref_chrom) if ref_chrom else None
    probe_validated = (load_probe_row(Path(snakemake.input.alignment), ref_chrom, row["P_seq"], ref_seqs[ref_chrom])
                        if ref_chrom else None)

    oligos = [("F", row["F_id"], row["F_seq"]), ("P", row["P_id"], row["P_seq"]), ("R", row["R_id"], row["R_seq"])]

    with open(snakemake.output.binding, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["role", "oligo_id", "seq", "found", "chrom", "start", "end", "strand", "mismatches"])
        for role, oid, seq in oligos:
            best = None
            if role in ("F", "R") and fr_validated is not None:
                pos = fr_validated[f"{role.lower()}_pos"]
                best = {"chrom": ref_chrom, "start": pos["start"], "end": pos["end"],
                        "strand": pos["strand"], "mismatches": fr_validated[f"{role.lower()}_mismatch"]}
            elif role == "P" and probe_validated is not None:
                best = {"chrom": ref_chrom, "start": probe_validated["start"], "end": probe_validated["end"],
                        "strand": probe_validated["strand"], "mismatches": probe_validated["mismatches"]}

            if best is None:
                # Not in this set's own already-validated population results (mismatch too
                # high for what that search covers, or genuinely never found there either) --
                # fall back to this script's own independent, wider-tolerance search.
                oligo_fa = out_dir / f"{oid}.fa"
                write_oligo_fasta(oligo_fa, [(oid, seq)])
                out_txt = out_dir / f"{oid}_search.txt"
                best = search_best_hit(mfeprimer_bin, ref_fasta, oligo_fa, out_txt, max_mismatch)
                if best is None and ref_chrom is not None:
                    # Same seed-finding-miss class direct_best_match() handles elsewhere in
                    # this pipeline -- last resort when even the wide independent search above
                    # finds nothing at all.
                    fallback = direct_best_match(seq, ref_seqs[ref_chrom])
                    if fallback and fallback[0] <= max_mismatch:
                        mm, i, strand = fallback
                        best = {"chrom": ref_chrom, "start": i, "end": i + len(seq), "strand": strand, "mismatches": mm}

            if best is None:
                writer.writerow([role, oid, seq, "0", "", "", "", "", ""])
            else:
                writer.writerow([role, oid, seq, "1", best["chrom"], best["start"], best["end"],
                                  best["strand"], best["mismatches"]])


if __name__ == "__main__":
    main()
