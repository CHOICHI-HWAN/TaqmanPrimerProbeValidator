#!/usr/bin/env python3
"""Detect ALL geometrically valid primer/probe binding combinations across every member
set of a group, not just each set's own independently-computed F<->R<->P match.

A simple union of each set's own validated genomes silently assumes Set1-F only ever
pairs with Set1-R, and that Set1-P only ever sits inside Set1's own amplicon. When two
sets' binding sites are close together that assumption can be wrong -- Set1-F may pair
with Set2-R (or vice versa), and Set2's probe may land inside Set1's amplicon. This
script pools every member set's F/R oligos into ONE mfeprimer --virus search (which finds
every geometrically compatible pair regardless of which set each oligo came from) and
every member set's probe into one wide `search`, then reports, per genome, every amplicon
found -- labeled with exactly which set's F paired with which set's R, and which set's
probe (if any) falls inside it.
"""
import csv
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import (  # noqa: E402
    find_mfeprimer, write_oligo_fasta, load_primers_manifest, iupac_mismatch_count,
    reverse_complement, run_primer_search, ensure_index, read_fasta_seqs,
    _oligo_sites, prefetch_sites, terminal_match_ok, scan_probe, choose_probe,
)
from build_alignment_table import probe_wide_search  # noqa: E402

ID_SEP = "::"


def mismatch_for(oligo_seq: str, matched_seq: str, is_fp: bool) -> int:
    query = oligo_seq if is_fp else reverse_complement(oligo_seq)
    return iupac_mismatch_count(query, matched_seq)


def resolve_oligo(name: str, oligo_ids_by_len):
    """mfeprimer suffixes a `.{variant}` index onto whatever id we gave it (degenerate
    primers expand into many concrete variants) -- match by longest known id that's a
    prefix of the reported name, so e.g. 'SET::F.3' resolves to 'SET::F'."""
    for oid in oligo_ids_by_len:
        if name == oid or name.startswith(oid + "."):
            return oid
    return None


def scan_combos(oligos: dict, genome: str, three_prime_n: int, max_mm: int = 2, max_size: int = 2000):
    """Independent whole-genome search of EVERY pairing of the group's F/R oligos (no mfeprimer): each oligo's sites
    (both strands, <= max_mm mismatches) are found once, then every facing pair within max_size is a combo. The best
    pair per (left oligo, right oligo) is kept -- 3' ends intact first, then fewer mismatches, then the shorter
    amplicon. Same dict shape as parse_combo_spec()'s combos (Tm/dG empty: only mfeprimer computes them)."""
    g = genome.upper()
    sites = {oid: _oligo_sites(info["seq"].upper(), g, max_mm) for oid, info in oligos.items()}
    best = {}
    for a, a_sites in sites.items():
        for b, b_sites in sites.items():
            if a == b:
                continue
            for left in a_sites:
                if left[2] != "+":
                    continue
                for right in b_sites:
                    if right[2] != "-" or right[0] < left[0] or right[1] <= left[0] or right[1] - left[0] > max_size:
                        continue
                    a_q, b_q = oligos[a]["seq"].upper(), reverse_complement(oligos[b]["seq"].upper())
                    a_ok = terminal_match_ok(a_q, g[left[0]:left[1]], True, three_prime_n)
                    b_ok = terminal_match_ok(b_q, g[right[0]:right[1]], False, three_prime_n)
                    key = (0 if (a_ok and b_ok) else 1, left[3] + right[3], right[1] - left[0])
                    if (a, b) not in best or key < best[(a, b)][0]:
                        best[(a, b)] = (key, {
                            "total": left[3] + right[3],
                            "f_oligo_id": a, "f_set_id": oligos[a]["set_id"], "f_mismatch": left[3], "f_tm": "", "f_dg": "", "f_3p_ok": a_ok,
                            "r_oligo_id": b, "r_set_id": oligos[b]["set_id"], "r_mismatch": right[3], "r_tm": "", "r_dg": "", "r_3p_ok": b_ok,
                            "amp_start": str(left[0] + 1), "amp_end": str(right[1]), "source": "genome-scan"})
    return {k: v[1] for k, v in best.items()}


def _combo_key(c):
    return (0 if (c.get("f_3p_ok") is not False and c.get("r_3p_ok") is not False) else 1, c["total"])


def parse_combo_spec(spec_path: Path, oligos: dict, genome_seqs: dict = None, three_prime_n: int = 0):
    """chrom -> list of combo dicts, one per distinct (f_oligo_id, r_oligo_id) pair
    mfeprimer reports for that genome (a genome can have more than one geometrically
    valid pairing, e.g. its own set's pair AND an unintended cross-set pair).

    A degenerate oligo expands into many concrete primer variants, and mfeprimer reports
    each variant's own match separately -- so the SAME (chrom, f_oligo_id, r_oligo_id)
    combo can appear dozens of times, differing only in which variant (and thus Tm/dG)
    happened to match. Keeping every row exploded a single 4-set group's table to ~300k
    rows for ~660 genomes; keying by (chrom, f_oligo_id, r_oligo_id) and keeping only the
    lowest-total-mismatch row per key collapses that back to one row per genuinely
    distinct combo per genome, mirroring how the per-set alignment table already keeps
    only the best row per chrom."""
    oligo_ids_by_len = sorted(oligos.keys(), key=len, reverse=True)
    best = {}
    if not spec_path.exists():
        return {}
    with open(spec_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            chrom = cols[1]
            amp_start, amp_end = cols[2], cols[3]
            fp_name, fp_seq = cols[6], cols[9]
            fp_tm, fp_dg = cols[10], cols[12]
            rp_name, rp_seq = cols[13], cols[16]
            rp_tm, rp_dg = cols[17], cols[19]
            if genome_seqs is not None and chrom in genome_seqs:
                # mfeprimer's fpSeq/rpSeq is the primer (and only reflects mismatches in its 3' k-mer): count against the
                # GENOME bases at the reported coordinates instead (1-based inclusive, top strand)
                fp_seq = genome_seqs[chrom][int(cols[7]) - 1:int(cols[8])].upper()
                lo, hi = sorted((int(cols[14]), int(cols[15])))
                rp_seq = genome_seqs[chrom][lo - 1:hi].upper()

            f_oid = resolve_oligo(fp_name, oligo_ids_by_len)
            r_oid = resolve_oligo(rp_name, oligo_ids_by_len)
            if f_oid is None or r_oid is None:
                continue
            f_info, r_info = oligos[f_oid], oligos[r_oid]
            f_mm = mismatch_for(f_info["seq"], fp_seq, True)
            r_mm = mismatch_for(r_info["seq"], rp_seq, False)
            f_ok = terminal_match_ok(f_info["seq"].upper(), fp_seq.upper(), True, three_prime_n)
            r_ok = terminal_match_ok(reverse_complement(r_info["seq"].upper()), rp_seq.upper(), False, three_prime_n)
            key = (chrom, f_oid, r_oid)
            total = f_mm + r_mm
            cand = {
                "total": total,
                "f_oligo_id": f_oid, "f_set_id": f_info["set_id"], "f_mismatch": f_mm, "f_3p_ok": f_ok,
                "f_tm": fp_tm, "f_dg": fp_dg,
                "r_oligo_id": r_oid, "r_set_id": r_info["set_id"], "r_mismatch": r_mm, "r_3p_ok": r_ok,
                "r_tm": rp_tm, "r_dg": rp_dg,
                "amp_start": amp_start, "amp_end": amp_end, "source": "mfeprimer",
            }
            prev = best.get(key)
            if prev is not None and _combo_key(prev) <= _combo_key(cand):
                continue
            best[key] = cand

    out = {}
    for (chrom, _f_oid, _r_oid), combo in best.items():
        out.setdefault(chrom, []).append(combo)
    return out


def compute_group_combo_table(*, db_fasta: Path, metadata_path: Path, rows_by_set: dict,
                               mfeprimer_bin: str, out_dir: Path, output_table: Path,
                               gene_summary_path: Path = None, apply_size_filter: bool = True,
                               tag: str = "combined_FR", three_prime_n: int = 2):
    """Pool every member set's F/R/Probe, search against `db_fasta`, and write the
    per-genome combo table to `output_table`. Used for the population db (group_cross_binding
    rule, with amplicon-size filtering) and reused as-is for the tiny standard-material and
    single-reference-genome databases (no size filtering needed there -- any match on a
    2-4 record or 1-record db is rare and worth surfacing directly, and there's no
    population of same-set matches to derive a size baseline from anyway)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    ensure_index(mfeprimer_bin, db_fasta)

    fr_oligos = {}
    fr_records = []
    for set_id, row in rows_by_set.items():
        f_full = f"{set_id}{ID_SEP}F"
        r_full = f"{set_id}{ID_SEP}R"
        fr_oligos[f_full] = {"set_id": set_id, "role": "F", "seq": row["F_seq"]}
        fr_oligos[r_full] = {"set_id": set_id, "role": "R", "seq": row["R_seq"]}
        fr_records.append((f_full, row["F_seq"]))
        fr_records.append((r_full, row["R_seq"]))

    fr_fa = out_dir / f"{tag}.fa"
    write_oligo_fasta(fr_fa, fr_records)
    fr_out = out_dir / f"{tag}_mm2.txt"
    run_primer_search(mfeprimer_bin, db_fasta, fr_fa, fr_out, mismatch=2)
    spec_path = Path(str(fr_out) + ".spec.tsv")

    # amplicon_seq is sliced directly from the db by (spec.tsv-reported) coordinate rather
    # than looked up in mfeprimer's own --fasta output, which silently truncates after
    # --max-amp-count (default 10000) amplicons -- see build_alignment_table.py's main()
    # for the same fix applied to the per-set table.
    db_seqs = read_fasta_seqs(db_fasta)
    combos_by_chrom = parse_combo_spec(spec_path, fr_oligos, genome_seqs=db_seqs, three_prime_n=three_prime_n)
    # Both engines on every genome (as for the single sets): mfeprimer finds sites through a 3' k-mer seed and can
    # miss one, so also scan the whole genome for every pairing and keep the better combo per (left, right) oligo.
    merged = {}
    for chrom, combos in combos_by_chrom.items():
        for c in combos:
            merged[(chrom, c["f_oligo_id"], c["r_oligo_id"])] = c
    prefetch_sites([(info["seq"], 2) for info in fr_oligos.values()] + [(rows_by_set[s]["P_seq"], 8) for s in rows_by_set], db_seqs.values())   # config scan_engine: go
    for chrom, seq in db_seqs.items():
        for (a, b), c in scan_combos(fr_oligos, seq, three_prime_n).items():
            prev = merged.get((chrom, a, b))
            if prev is None or _combo_key(c) < _combo_key(prev):
                merged[(chrom, a, b)] = c
    combos_by_chrom = {}
    for (chrom, _a, _b), c in merged.items():
        combos_by_chrom.setdefault(chrom, []).append(c)

    size_filter_info = None
    if apply_size_filter:
        # -- discover-then-filter: the first (only) mfeprimer run above used its default,
        # generous 2000bp size ceiling, so it also reports geometrically-valid-but-absurd
        # pairings (e.g. two primers from different genes, 1-2kb apart) that could never
        # actually amplify. There's no single "too big" size that applies to every group --
        # different assays (different genes/segments) have genuinely different real
        # amplicon sizes -- so instead of guessing a fixed cutoff, derive it from this
        # group's OWN correctly-designed pairs. A cross-set combo is only kept if its size
        # is within a generous margin of that. This is a second, distinct filtering pass;
        # re-running mfeprimer itself with a smaller -S/--maxSize would find nothing this
        # doesn't already find, since a stricter size cap can only be a subset of the first
        # (already maximally permissive) run.
        #
        # Prefer each member set's reference-genome-verified amplicon_size (gene_summary.tsv,
        # from a single curated reference genome -- exact, not a statistic) over deriving a
        # size from this group's own population same-set matches. Population-derived sizes
        # looked reliable at first but turned out not to be: a degenerate primer can match
        # multiple, genuinely different loci across the population, and those get mixed into
        # the same "same-set" bucket, inflating the very median meant to define "normal"
        # (seen e.g. for VendorA_L_ALL: population median 855bp vs a true reference-verified
        # ~110bp). Population same-set sizes are now only a fallback for a set
        # gene_summary.tsv has no reference-genome match for at all.
        ref_sizes_by_set = {}
        if gene_summary_path and gene_summary_path.exists():
            with open(gene_summary_path, newline="", encoding="utf-8") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r["set_id"] in rows_by_set and r.get("amplicon_size", "NA") not in ("", "NA"):
                        try:
                            ref_sizes_by_set[r["set_id"]] = int(r["amplicon_size"])
                        except ValueError:
                            pass

        own_sizes_by_set = {}
        for combos in combos_by_chrom.values():
            for combo in combos:
                if combo["f_set_id"] == combo["r_set_id"]:
                    size = int(combo["amp_end"]) - int(combo["amp_start"]) + 1
                    own_sizes_by_set.setdefault(combo["f_set_id"], []).append(size)

        size_sources = {}  # set_id -> (size, source)
        for set_id in rows_by_set:
            if set_id in ref_sizes_by_set:
                size_sources[set_id] = (ref_sizes_by_set[set_id], "refgenome")
            elif set_id in own_sizes_by_set:
                size_sources[set_id] = (statistics.median(own_sizes_by_set[set_id]), "population_median")

        if size_sources:
            anchor_size = max(v[0] for v in size_sources.values())
            ideal_max_size = max(300, round(anchor_size * 2))
        else:
            anchor_size = None
            ideal_max_size = 500  # no size info for any member at all -- generous fallback

        all_combo_count = sum(len(v) for v in combos_by_chrom.values())
        for chrom in list(combos_by_chrom.keys()):
            combos_by_chrom[chrom] = [
                c for c in combos_by_chrom[chrom]
                if (int(c["amp_end"]) - int(c["amp_start"]) + 1) <= ideal_max_size
            ]
            if not combos_by_chrom[chrom]:
                del combos_by_chrom[chrom]
        kept_combo_count = sum(len(v) for v in combos_by_chrom.values())

        size_filter_info = {
            "per_set_sizes": {sid: {"size": size, "source": src} for sid, (size, src) in size_sources.items()},
            "anchor_size": anchor_size,
            "ideal_max_size": ideal_max_size,
            "combos_before_size_filter": all_combo_count,
            "combos_after_size_filter": kept_combo_count,
        }

    # -- every member set's probe, searched independently (each probe is its own oligo,
    # unaffected by which F/R pair it ends up co-located with) --------------------------
    probe_hits_by_set = {}
    for set_id, row in rows_by_set.items():
        probe_out = out_dir / f"{tag}_{set_id}_probe_wide.txt"
        probe_hits_by_set[set_id] = probe_wide_search(
            mfeprimer_bin, db_fasta, f"{set_id}{ID_SEP}P", row["P_seq"], probe_out)

    meta = {}
    if metadata_path and metadata_path.exists():
        with open(metadata_path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                meta[r["accession"]] = r

    with open(output_table, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "genome_id", "lineage", "date", "country", "region",
            "combo_id", "is_same_set",
            "f_set_id", "f_oligo_id", "f_mismatch", "f_tm", "f_dg",
            "r_set_id", "r_oligo_id", "r_mismatch", "r_tm", "r_dg",
            "amplicon_start", "amplicon_end", "amplicon_seq",
            "p_set_id", "p_oligo_id", "p_mismatch", "p_tm", "p_dg", "probe_start", "probe_end",
            "f_3p_ok", "r_3p_ok", "source",
        ])
        for chrom in sorted(combos_by_chrom):
            m = meta.get(chrom, {})
            for combo in combos_by_chrom[chrom]:
                amp_start, amp_end = combo["amp_start"], combo["amp_end"]
                amp_seq = ""
                full_seq = db_seqs.get(chrom, "")
                a_start, a_end = int(amp_start), int(amp_end)
                if full_seq and 1 <= a_start <= a_end <= len(full_seq):
                    amp_seq = full_seq[a_start - 1:a_end]

                # Each member probe, both engines: mfeprimer's best hit (only usable when it lies inside THIS amplicon)
                # and a whole-genome scan restricted to the amplicon; the best probe over all members wins.
                best_probe = None
                for set_id, hits in probe_hits_by_set.items():
                    hit = hits.get(chrom)
                    if hit and not (int(hit["start"]) >= a_start and int(hit["end"]) <= a_end):
                        hit = None
                    own = scan_probe(rows_by_set[set_id]["P_seq"], full_seq, a_start, a_end, max_mm=8) if full_seq else None
                    chosen, _src = choose_probe(hit, own)
                    if chosen is not None and (best_probe is None or chosen["mismatch"] < best_probe["mismatch"]):
                        best_probe = {"set_id": set_id, "oligo_id": f"{set_id}{ID_SEP}P", **chosen}

                # The oligo that ends up playing the geometric "left/fp" role isn't always
                # that set's own F -- its own R can end up there instead (fp/rp is purely
                # about which end sits left vs right on the top strand, not our F/R label).
                # f_oligo_id/r_oligo_id already carry the REAL resolved role (e.g.
                # "SET::R"), so read it from there rather than hardcoding "-F"/"-R" -- doing
                # that previously made two genuinely different pairings (e.g. a set's own
                # F+R normally, vs the same two oligos in the opposite/reversed roles at a
                # different genomic span) render with an identical, misleading label.
                f_role = combo["f_oligo_id"].split(ID_SEP)[-1]
                r_role = combo["r_oligo_id"].split(ID_SEP)[-1]
                writer.writerow([
                    chrom, m.get("lineage", "NA"), m.get("date", ""), m.get("country", ""), m.get("region", ""),
                    f"{combo['f_set_id']}-{f_role} x {combo['r_set_id']}-{r_role}",
                    "1" if combo["f_set_id"] == combo["r_set_id"] else "0",
                    combo["f_set_id"], combo["f_oligo_id"], combo["f_mismatch"], combo["f_tm"], combo["f_dg"],
                    combo["r_set_id"], combo["r_oligo_id"], combo["r_mismatch"], combo["r_tm"], combo["r_dg"],
                    amp_start, amp_end, amp_seq,
                    best_probe["set_id"] if best_probe else "",
                    best_probe["oligo_id"] if best_probe else "",
                    best_probe["mismatch"] if best_probe else "",
                    best_probe["tm"] if best_probe else "",
                    best_probe["dg"] if best_probe else "",
                    best_probe["start"] if best_probe else "",
                    best_probe["end"] if best_probe else "",
                    "" if combo.get("f_3p_ok") is None else int(bool(combo["f_3p_ok"])),
                    "" if combo.get("r_3p_ok") is None else int(bool(combo["r_3p_ok"])),
                    combo.get("source", ""),
                ])

    return size_filter_info


def main():
    group_id = snakemake.wildcards.group_id
    manifest = load_primers_manifest(snakemake.input.manifest)
    member_ids = snakemake.params.members
    rows_by_set = {r["set_id"]: r for r in manifest if r["set_id"] in member_ids}

    out_dir = Path(snakemake.output.table).parent
    pipeline_root = Path(snakemake.scriptdir).parent
    mfeprimer_bin = find_mfeprimer(snakemake.config.get("mfeprimer_bin", ""), pipeline_root)

    size_filter_info = compute_group_combo_table(
        db_fasta=Path(snakemake.input.db),
        metadata_path=Path(snakemake.input.metadata),
        rows_by_set=rows_by_set,
        mfeprimer_bin=mfeprimer_bin,
        out_dir=out_dir,
        output_table=Path(snakemake.output.table),
        gene_summary_path=Path(snakemake.input.gene_summary),
        apply_size_filter=True,
        tag="combined_FR",
        three_prime_n=int(snakemake.config.get("three_prime_consecutive_match", 2)),
    )
    (out_dir / "size_filter.json").write_text(json.dumps(size_filter_info, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
