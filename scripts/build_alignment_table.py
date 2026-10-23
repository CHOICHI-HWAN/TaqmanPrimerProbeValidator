#!/usr/bin/env python3
"""Build a per-genome F/Probe/R mismatch table for one assay set against the full
population database -- the data backing the interactive "how much does this set mismatch
across every genome" view (dot-plot / problems list), modeled on the reference Nipah
henipavirus example screenshot.

Reuses the primer-pair spec.tsv + amplicon fasta already produced by primer_match (mm=2)
for F/R, and runs one additional wide-mismatch probe `search` for P, rather than
recomputing anything from scratch.
"""
import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import (  # noqa: E402
    find_mfeprimer, write_oligo_fasta, load_primers_manifest, parse_spec_best,
    read_fasta_ids, read_fasta_seqs, prefetch_sites, clear_site_cache, direct_best_match, direct_best_fr_match, scan_fr, choose_fr, scan_probe, choose_probe,
)

# A probe missing from mfeprimer's own `search -m 8` results despite the amplicon (from F/R)
# existing means its seed-finding step failed, not that the site is truly absent -- direct_
# best_match() re-checks by brute force within the amplicon itself. Only bother falling back
# when the direct hit is this tight; a genuinely absent probe can still align somewhere in a
# ~100bp window at 10+ mismatches by chance, which isn't a meaningful "it actually binds".
PROBE_FALLBACK_MAX_MISMATCH = 8

PROBE_WIDE_MISMATCH = 8


def probe_wide_search(mfeprimer_bin, db_fasta, probe_id, probe_seq, out_txt):
    probe_fa = out_txt.with_suffix(".fa")
    write_oligo_fasta(probe_fa, [(probe_id, probe_seq)])
    if out_txt.exists():
        out_txt.unlink()
    cmd = [mfeprimer_bin, "search", str(db_fasta), str(probe_fa), str(out_txt), "-m", str(PROBE_WIDE_MISMATCH)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] probe wide search failed: {result.stderr}")
    best = {}
    with open(out_txt, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, start, end, mismatch = cols[0], cols[1], cols[2], int(cols[6])
            dg, tm = cols[7], cols[8]
            # mfeprimer `search`'s own header says "Start(0-based)" (half-open [start,end)),
            # but amp_start/amp_end (from the F/R --virus spec.tsv) are 1-based inclusive --
            # convert here so probe_start/probe_end land in the SAME coordinate system,
            # otherwise probe_start ends up shifted by one base against the amplicon it's
            # supposed to be a sub-span of (start+1 makes it 1-based inclusive; end is
            # already numerically correct since a half-open end equals the 1-based
            # inclusive end).
            start = str(int(start) + 1)
            if chrom not in best or mismatch < best[chrom]["mismatch"]:
                best[chrom] = {"mismatch": mismatch, "start": start, "end": end, "tm": tm, "dg": dg}
    return best


COMPARISON_COLUMNS = ["genome_id", "mfe_f", "mfe_r", "mfe_3p_ok", "mfe_amplicon", "own_f", "own_r", "own_3p_ok",
                      "own_amplicon", "source", "verdict_mfe", "verdict_own", "agree",
                      "mfe_probe", "own_probe", "probe_source", "probe_verdict_mfe", "probe_verdict_own", "probe_agree"]


def _verdict(x, level=1):
    """Primer-pair pass at `level` mismatches each with intact 3' ends (what the final verdict uses)."""
    return bool(x and x["f_mismatch"] <= level and x["r_mismatch"] <= level
                and x.get("f_3p_ok") is not False and x.get("r_3p_ok") is not False)


def engine_comparison_row(gid, mfe_fr, own_fr, source):
    def cols(x):
        if not x:
            return ["", "", "", ""]
        return [x["f_mismatch"], x["r_mismatch"], int(x.get("f_3p_ok") is not False and x.get("r_3p_ok") is not False),
                f"{x['amp_start']}-{x['amp_end']}"]
    vm, vo = _verdict(mfe_fr), _verdict(own_fr)
    return [gid] + cols(mfe_fr) + cols(own_fr) + [source, int(vm), int(vo), int(vm == vo)]


def probe_comparison_cols(mfe_hit, own_hit, source):
    """mfe_hit: mfeprimer's best probe hit that lies in the amplicon (or None); own_hit: the genome scan's."""
    pm = None if mfe_hit is None else mfe_hit["mismatch"]
    po = None if own_hit is None else own_hit["mismatch"]
    vm, vo = int(pm == 0), int(po == 0)
    return ["" if pm is None else pm, "" if po is None else po, source, vm, vo, int(vm == vo)]


def main():
    set_id = snakemake.wildcards.set_id
    manifest = load_primers_manifest(Path(snakemake.input.manifest))
    row = next(r for r in manifest if r["set_id"] == set_id)
    db_fasta = Path(snakemake.input.db)
    if hasattr(snakemake.input, "spec"):
        spec_path = Path(snakemake.input.spec)
    else:
        # standard_alignment_table passes the always-present main report instead (its
        # .spec.tsv sibling is only conditionally written by mfeprimer, so it can't be a
        # hard Snakemake input -- see that rule's comment).
        report_path = Path(snakemake.input.report)
        spec_path = Path(str(report_path) + ".spec.tsv")
    metadata_path = Path(snakemake.input.metadata)

    pipeline_root = Path(snakemake.scriptdir).parent
    mfeprimer_bin = find_mfeprimer(snakemake.config.get("mfeprimer_bin", ""), pipeline_root)
    # For the tiny (2-4 record) KDCA standard-material db, a material that fails to bind
    # under ANY tolerance we searched is itself the finding -- so standard_alignment_table
    # passes include_all=True to keep every record (blank mismatch/Tm/dG columns for a
    # total non-match) instead of silently omitting it the way the population table does.
    include_all = bool(getattr(snakemake.params, "include_all", False))
    three_prime_n = int(getattr(snakemake.params, "three_prime_n", 0) or 0)

    all_ids = sorted(read_fasta_ids(db_fasta))

    db_seqs = read_fasta_seqs(db_fasta)
    fr_best = parse_spec_best(spec_path, row["F_id"], row["F_seq"], row["R_id"], row["R_seq"],
                               three_prime_n=three_prime_n, genome_seqs=db_seqs)
    # amp_path/parse_amp_fasta_all is NOT used for amplicon_seq: mfeprimer's --fasta output
    # silently truncates after --max-amp-count (default 10000) amplicon candidates, and a
    # highly degenerate primer set (many IUPAC-expanded variants) can produce far more
    # candidate rows than that in spec.tsv alone -- e.g. VendorC_L_1 had ~197k spec.tsv rows
    # but only the first 10000 ever got written to .amp.fasta, silently dropping the
    # amplicon sequence for every genome processed after that cutoff (confirmed: 48/364
    # genomes for VendorC_L_1, all with otherwise-valid f/p/r_mismatch data -- the missing
    # amplicon_seq made them vanish from the MSA without being counted in `excluded`,
    # while still showing up normally in the per-genome mismatch table, i.e. exactly the
    # "these two views don't agree" report). spec.tsv's ampStart/ampEnd are never
    # truncated, so slicing the amplicon directly out of the population genome fasta by
    # coordinate sidesteps mfeprimer's cap entirely, is exact, and needs no separate file.
    probe_out = Path(snakemake.output.table).with_name(f"{set_id}_probe_wide_search.txt")
    comparison_rows = []
    p_best = probe_wide_search(mfeprimer_bin, db_fasta, row["P_id"], row["P_seq"], probe_out)

    meta = {}
    with open(metadata_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["accession"]] = r

    with open(snakemake.output.table, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["genome_id", "lineage", "date", "country", "region",
                          "f_mismatch", "p_mismatch", "r_mismatch", "f_3p_ok", "r_3p_ok",
                          "amplicon_start", "amplicon_end", "probe_start", "probe_end", "amplicon_seq",
                          "f_tm", "f_dg", "p_tm", "p_dg", "r_tm", "r_dg", "fr_source"])
        prefetch_sites([(row["F_seq"], 2), (row["R_seq"], 2), (row["P_seq"], PROBE_FALLBACK_MAX_MISMATCH)], db_seqs.values())   # config scan_engine: go
        for gid in all_ids:
            # Run BOTH engines on every genome and keep the better pair: mfeprimer finds sites through a 3' k-mer
            # seed and can miss a primer's true best site; scan_fr() checks every window of the genome.
            mfe_fr = fr_best.get(gid)
            own_fr = scan_fr(row["F_seq"], row["R_seq"], db_seqs.get(gid, ""), three_prime_n=three_prime_n)
            fr, fr_source = choose_fr(mfe_fr, own_fr)
            comparison_rows.append(engine_comparison_row(gid, mfe_fr, own_fr, fr_source))
            f_mm = fr["f_mismatch"] if fr else None
            r_mm = fr["r_mismatch"] if fr else None
            f_3p_ok = fr["f_3p_ok"] if fr else None
            r_3p_ok = fr["r_3p_ok"] if fr else None
            amp_start = fr["amp_start"] if fr else ""
            amp_end = fr["amp_end"] if fr else ""
            f_tm = fr["f_tm"] if fr else ""
            f_dg = fr["f_dg"] if fr else ""
            r_tm = fr["r_tm"] if fr else ""
            r_dg = fr["r_dg"] if fr else ""
            p_hit = p_best.get(gid)
            p_mm = p_hit["mismatch"] if p_hit else None
            p_start = p_hit["start"] if p_hit else ""
            p_end = p_hit["end"] if p_hit else ""
            p_tm = p_hit["tm"] if p_hit else ""
            p_dg = p_hit["dg"] if p_hit else ""
            amp_seq = ""
            if fr:
                full_seq = db_seqs.get(gid, "")
                a_start, a_end = int(amp_start), int(amp_end)
                if full_seq and 1 <= a_start <= a_end <= len(full_seq):
                    amp_seq = full_seq[a_start - 1:a_end]
            # Probe, both engines: mfeprimer's best hit and a whole-genome scan restricted to the amplicon. A probe that
            # binds outside the amplicon does not count, so mfeprimer's best-anywhere hit is only used when it lies inside.
            full_seq = db_seqs.get(gid, "")
            a_lo, a_hi = (int(amp_start), int(amp_end)) if fr else (None, None)
            def inside(h):
                return h is not None and (a_lo is None or (a_lo <= int(h["start"]) and int(h["end"]) <= a_hi))
            mfe_probe = p_hit if inside(p_hit) else None
            own_probe = scan_probe(row["P_seq"], full_seq, a_lo, a_hi, max_mm=PROBE_FALLBACK_MAX_MISMATCH) if full_seq else None
            chosen, probe_source = choose_probe(mfe_probe, own_probe)
            if chosen is None and p_hit is not None:
                chosen, probe_source = p_hit, "mfeprimer"          # only outside the amplicon: shown, but reconcile rejects it
            comparison_rows[-1] += probe_comparison_cols(mfe_probe, own_probe, probe_source)
            if chosen is not None:
                p_mm, p_start, p_end, p_tm, p_dg = chosen["mismatch"], chosen["start"], chosen["end"], chosen["tm"], chosen["dg"]
            if fr is None and chosen is None and not include_all:
                continue  # not found under any tolerance we searched -- omit rather than pad with fake data
            m = meta.get(gid, {})
            writer.writerow([
                gid, m.get("lineage", "NA"), m.get("date", ""), m.get("country", ""), m.get("region", ""),
                f_mm if f_mm is not None else "", p_mm if p_mm is not None else "",
                r_mm if r_mm is not None else "",
                "1" if f_3p_ok else ("" if f_3p_ok is None else "0"),
                "1" if r_3p_ok else ("" if r_3p_ok is None else "0"),
                amp_start, amp_end, p_start, p_end, amp_seq,
                f_tm, f_dg, p_tm, p_dg, r_tm, r_dg, fr_source,
            ])
    with open(snakemake.output.comparison, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(COMPARISON_COLUMNS)
        writer.writerows(comparison_rows)


if __name__ == "__main__":
    main()
