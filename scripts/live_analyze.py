#!/usr/bin/env python3
"""Run the full validation pipeline for one ad-hoc (user-submitted) F/R/Probe set against
the already-built population / standard-material / reference-genome databases, without
going through Snakemake. Used by the interactive webapp so a user can type in a new
primer set and see results in the same shape as the 14 built-in sets.

Everything here reuses the exact same mfeprimer invocations and comparison logic as the
Snakemake pipeline (scripts/run_primer.py, run_probe.py, reconcile.py,
run_refgenome_binding.py, build_gene_summary.py, build_alignment_table.py) -- just called
as plain functions instead of as separate Snakemake `script:` steps, since the whole
point is to run interactively without a DAG.
"""
import csv
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import (  # noqa: E402
    find_mfeprimer, ensure_index, write_oligo_fasta, run_primer_search, parse_spec_best,
    read_fasta_ids, read_fasta_seqs, direct_best_match, direct_best_fr_match, scan_fr, choose_fr, scan_probe, choose_probe,
)
from msa_common import build_msa  # noqa: E402
from run_dimer_hairpin import (  # noqa: E402
    run_mfeprimer, parse_dimer_report, parse_hairpin_report,
    summarize_dimers, summarize_hairpins, DG_FLAG_THRESHOLD,
)

VALID_SEQ_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)
GENE_TABLE_CACHE = {}
# See build_alignment_table.py's identical constant/comment: mfeprimer's probe `search` can
# fail to find a genuinely low-mismatch site (seed-finding miss, not a real absence) -- this
# caps how far direct_best_match()'s brute-force fallback is trusted before treating a hit
# as coincidental rather than real.
PROBE_FALLBACK_MAX_MISMATCH = 8


def _probe_hit_with_fallback(p_hit, p_seq, amp_seq, amp_start, genome_seq=None, amp_end=None):
    """The probe hit for one genome, the same way the Snakemake pipeline decides it: mfeprimer's best hit (only usable when
    it lies inside the F/R amplicon) and a whole-genome scan restricted to the amplicon, the better one wins. Without an
    amplicon (no F/R pair) mfeprimer's hit is returned as is."""
    if amp_start is None or not genome_seq or amp_end is None:
        return p_hit
    inside = p_hit is not None and amp_start <= int(p_hit["start"]) and int(p_hit["end"]) <= amp_end
    own = scan_probe(p_seq, genome_seq, amp_start, amp_end, max_mm=PROBE_FALLBACK_MAX_MISMATCH)
    chosen, _src = choose_probe(p_hit if inside else None, own)
    if chosen is None:
        return p_hit                                  # only outside the amplicon: shown, but it does not count as bound
    return {"mismatch": chosen["mismatch"], "start": chosen["start"], "end": chosen["end"],
            "tm": chosen.get("tm") or None, "dg": chosen.get("dg") or None}


def validate_oligo(seq: str, label: str):
    seq = (seq or "").strip().upper()
    if not seq:
        raise ValueError(f"{label} sequence is empty")
    if not VALID_SEQ_RE.match(seq):
        raise ValueError(f"{label} sequence contains characters other than IUPAC nucleotide codes (ACGTRYSWKMBDHVN): {seq}")
    if len(seq) < 10 or len(seq) > 60:
        raise ValueError(f"{label} sequence length looks unusual ({len(seq)}bp) -- expected roughly 15-40bp")
    return seq


def load_genes(genes_tsv: Path):
    if genes_tsv in GENE_TABLE_CACHE:
        return GENE_TABLE_CACHE[genes_tsv]
    genes = {}
    with open(genes_tsv, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genes.setdefault(row["segment"], []).append(row)
    GENE_TABLE_CACHE[genes_tsv] = genes
    return genes


def genes_overlapping(genes_for_seg, start, end):
    hits = []
    for g in genes_for_seg:
        gs, ge = int(g["start"]), int(g["end"])
        if start <= ge and end >= gs:
            hits.append(g["gene"])
    return list(dict.fromkeys(hits))


## parse_spec_best now lives in mfe_common.py, shared with build_alignment_table.py so the
## Snakemake pipeline and this live/interactive path can never again independently drift
## (see mfe_common.parse_spec_best's docstring -- that's exactly how these two copies
## disagreed before).


def probe_search_best(mfeprimer_bin, db_fasta, probe_id, probe_seq, out_txt, mismatch):
    probe_fa = out_txt.with_suffix(".fa")
    write_oligo_fasta(probe_fa, [(probe_id, probe_seq)])
    if out_txt.exists():
        out_txt.unlink()
    cmd = [mfeprimer_bin, "search", str(db_fasta), str(probe_fa), str(out_txt), "-m", str(mismatch)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mfeprimer search failed: {result.stderr}")
    best = {}
    with open(out_txt, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, start, end, mm = cols[0], cols[1], cols[2], int(cols[6])
            dg, tm = cols[7], cols[8]
            if chrom not in best or mm < best[chrom]["mismatch"]:
                best[chrom] = {"mismatch": mm, "start": start, "end": end, "tm": tm, "dg": dg}
    return best


def refgenome_binding(mfeprimer_bin, ref_fasta, oligos, out_dir, max_mismatch=6):
    """oligos: list of (role, oligo_id, seq). Returns list of hit dicts."""
    ensure_index(mfeprimer_bin, ref_fasta)
    hits = []
    for role, oid, seq in oligos:
        oligo_fa = out_dir / f"ref_{role}.fa"
        write_oligo_fasta(oligo_fa, [(oid, seq)])
        out_txt = out_dir / f"ref_{role}_search.txt"
        if out_txt.exists():
            out_txt.unlink()
        cmd = [mfeprimer_bin, "search", str(ref_fasta), str(oligo_fa), str(out_txt), "-m", str(max_mismatch)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"refgenome search failed: {result.stderr}")
        best = None
        with open(out_txt, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                chrom, start, end, _q, _sz, strand, mm = line.rstrip("\n").split("\t")[:7]
                mm = int(mm)
                if best is None or mm < best["mismatches"]:
                    best = {"chrom": chrom, "start": int(start), "end": int(end), "strand": strand, "mismatches": mm}
        if best is None:
            hits.append({"role": role, "oligo_id": oid, "found": False})
        else:
            hits.append({"role": role, "oligo_id": oid, "found": True, **best})
    return hits


def analyze_set(*, name, segment, f_id, f_seq, r_id, r_seq, p_id, p_seq, root: Path, out_dir: Path,
                 primer_mismatch_pass=1, probe_mismatch=0, three_prime_n=0, refgenome_accession=None,
                 cpu=None, progress=None):
    """Returns a dict: {meta, summary, standard, gene, binding, genomes: [...]}."""
    def report(msg):
        if progress:
            progress(msg)

    seg = segment.lower()
    from virus_profile import load_for_root  # noqa: PLC0415
    allowed = list(load_for_root(root)["segments"])
    if seg not in allowed:
        raise ValueError(f"segment must be one of {allowed}")
    f_seq = validate_oligo(f_seq, "Forward primer")
    r_seq = validate_oligo(r_seq, "Reverse primer")
    p_seq = validate_oligo(p_seq, "Probe")

    out_dir.mkdir(parents=True, exist_ok=True)
    pipeline_root = root / "scripts"
    mfeprimer_bin = find_mfeprimer("", pipeline_root.parent)

    db_fasta = root / "data" / "filtered" / f"{seg}_filtered.fasta"
    db_meta_path = root / "data" / "filtered" / f"{seg}_filtered_metadata.tsv"
    standard_fasta = root / "config" / f"kdca_standard_{seg}.fasta"
    standard_meta_path = root / "config" / f"kdca_standard_{seg}_metadata.tsv"
    ref_fasta = root / "data" / "refgenome" / f"{seg}_ref.fasta"
    genes_tsv = root / "config" / "refgenome_genes.tsv"

    all_ids = set(read_fasta_ids(db_fasta))
    total = len(all_ids)

    # Everything below (summary coverage %, validated/probe_dropout ids, the KDCA standard
    # check) is derived from ONE wide search (F/R at mismatch<=2, probe at mismatch<=8)
    # with the exact per-genome mismatch count independently recomputed in Python --
    # exactly the same data the per-genome mismatch table/MSA use -- rather than running
    # separate direct searches at each specific mismatch level. Those used to disagree:
    # mfeprimer's own seed-finding step happens BEFORE its mismatch filtering, so a direct
    # search run at a low --misMatch can fail to find a site at all (report it as a hard
    # miss) in cases the wider search -- and the independent Python-side recount -- convict
    # as a genuine 0/1-mismatch match. See scripts/reconcile.py for the population-pipeline
    # side of the same fix, confirmed concretely against real data.
    report("F/R primer 매치 실행 중 (mismatch<=2, 정확한 개수는 재계산)...")
    fr_fa = out_dir / "FR.fa"
    write_oligo_fasta(fr_fa, [(f_id, f_seq), (r_id, r_seq)])
    ensure_index(mfeprimer_bin, db_fasta, cpu=cpu)
    fr_mm2_txt = out_dir / "primer_mm2.txt"
    run_primer_search(mfeprimer_bin, db_fasta, fr_fa, fr_mm2_txt, 2, cpu=cpu)
    db_seqs = read_fasta_seqs(db_fasta)
    fr_best = parse_spec_best(Path(str(fr_mm2_txt) + ".spec.tsv"), f_id, f_seq, r_id, r_seq,
                               three_prime_n=three_prime_n, genome_seqs=db_seqs)

    report("probe 매치 확인 중 (mismatch<=8, 정확한 개수는 재계산)...")
    probe_fa = out_dir / "probe.fa"
    write_oligo_fasta(probe_fa, [(p_id, p_seq)])
    probe_wide_out = out_dir / "probe_wide.txt"
    p_best = probe_search_best(mfeprimer_bin, db_fasta, p_id, p_seq, probe_wide_out, mismatch=8)

    db_meta = {}
    if db_meta_path.exists():
        with open(db_meta_path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                db_meta[r["accession"]] = r

    genomes = []
    for gid in sorted(all_ids):
        fr, _src = choose_fr(fr_best.get(gid), scan_fr(f_seq, r_seq, db_seqs.get(gid, ""), three_prime_n=three_prime_n))
        p_hit = p_best.get(gid)
        amp_seq = ""
        amp_start_int = None
        if fr:
            full_seq = db_seqs.get(gid, "")
            amp_start_int, a_end = int(fr["amp_start"]), int(fr["amp_end"])
            if full_seq and 1 <= amp_start_int <= a_end <= len(full_seq):
                amp_seq = full_seq[amp_start_int - 1:a_end]
        p_hit = _probe_hit_with_fallback(p_hit, p_seq, amp_seq, amp_start_int, db_seqs.get(gid, ""), int(fr["amp_end"]) if fr else None)
        if fr is None and p_hit is None:
            continue
        m = db_meta.get(gid, {})
        genomes.append({
            "genome_id": gid, "lineage": m.get("lineage", "NA"), "date": m.get("date", ""),
            "country": m.get("country", ""), "region": m.get("region", ""),
            "f_mismatch": fr["f_mismatch"] if fr else None,
            "p_mismatch": p_hit["mismatch"] if p_hit else None,
            "r_mismatch": fr["r_mismatch"] if fr else None,
            "f_3p_ok": fr["f_3p_ok"] if fr else None,
            "r_3p_ok": fr["r_3p_ok"] if fr else None,
            "amplicon_start": fr["amp_start"] if fr else None,
            "amplicon_end": fr["amp_end"] if fr else None,
            "probe_start": p_hit["start"] if p_hit else None,
            "probe_end": p_hit["end"] if p_hit else None,
            "amplicon_seq": amp_seq,
            "f_tm": fr["f_tm"] if fr else None, "f_dg": fr["f_dg"] if fr else None,
            "p_tm": p_hit["tm"] if p_hit else None, "p_dg": p_hit["dg"] if p_hit else None,
            "r_tm": fr["r_tm"] if fr else None, "r_dg": fr["r_dg"] if fr else None,
        })

    def _fr_ok(g, mm):
        return (g["f_mismatch"] is not None and g["r_mismatch"] is not None
                and g["f_mismatch"] <= mm and g["r_mismatch"] <= mm
                and g.get("f_3p_ok", True) and g.get("r_3p_ok", True))

    def _probe_ok(g, mm):
        return g["p_mismatch"] is not None and g["p_mismatch"] <= mm

    primer_matched_by_mm = {mm: {g["genome_id"] for g in genomes if _fr_ok(g, mm)} for mm in (0, 1, 2)}
    primer_matched_ids = primer_matched_by_mm[primer_mismatch_pass]
    probe_matched_ids = {g["genome_id"] for g in genomes if _probe_ok(g, probe_mismatch)}

    validated_ids = primer_matched_ids & probe_matched_ids
    probe_dropout_ids = primer_matched_ids - probe_matched_ids

    summary = {
        "total_genomes": total,
        "primer_coverage_mm0_pct": round(100 * len(primer_matched_by_mm[0]) / total, 2) if total else 0.0,
        "primer_coverage_mm1_pct": round(100 * len(primer_matched_by_mm[1]) / total, 2) if total else 0.0,
        "primer_coverage_mm2_pct": round(100 * len(primer_matched_by_mm[2]) / total, 2) if total else 0.0,
        "probe_coverage_mm0_pct": round(100 * len(probe_matched_ids) / total, 2) if total else 0.0,
        "primer_mismatch_used": primer_mismatch_pass,
        "primer_matched": len(primer_matched_ids),
        "probe_dropout": len(probe_dropout_ids),
        "validated": len(validated_ids),
        "validated_pct": round(100 * len(validated_ids) / total, 2) if total else 0.0,
    }

    report("표준물질(KDCA) 매치 확인 중...")
    standard = {"ok": False, "matched": []}
    standard_genomes = []
    if standard_fasta.exists():
        ensure_index(mfeprimer_bin, standard_fasta, cpu=cpu)
        std_ids = set(read_fasta_ids(standard_fasta))
        std_meta = {}
        if standard_meta_path.exists():
            with open(standard_meta_path, newline="", encoding="utf-8") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    std_meta[r["accession"]] = r

        # Same wide-search-then-recount approach as the population db above (see the
        # comment there) -- avoids the same direct-narrow-search seed-finding blind spot
        # for the "O/X 매치" pass/fail determination.
        std_fr2_txt = out_dir / "standard_primer_mm2.txt"
        run_primer_search(mfeprimer_bin, standard_fasta, fr_fa, std_fr2_txt, 2, cpu=cpu)
        std_db_seqs = read_fasta_seqs(standard_fasta)
        std_fr_best = parse_spec_best(Path(str(std_fr2_txt) + ".spec.tsv"), f_id, f_seq, r_id, r_seq,
                                       three_prime_n=three_prime_n, genome_seqs=std_db_seqs)
        std_probe_wide_out = out_dir / "standard_probe_wide.txt"
        std_p_best = probe_search_best(mfeprimer_bin, standard_fasta, p_id, p_seq, std_probe_wide_out, mismatch=8)
        for gid in sorted(std_ids):
            fr, _src = choose_fr(std_fr_best.get(gid), scan_fr(f_seq, r_seq, std_db_seqs.get(gid, ""), three_prime_n=three_prime_n))
            p_hit = std_p_best.get(gid)
            m = std_meta.get(gid, {})
            amp_seq = ""
            amp_start_int = None
            if fr:
                full_seq = std_db_seqs.get(gid, "")
                amp_start_int, a_end = int(fr["amp_start"]), int(fr["amp_end"])
                if full_seq and 1 <= amp_start_int <= a_end <= len(full_seq):
                    amp_seq = full_seq[amp_start_int - 1:a_end]
            p_hit = _probe_hit_with_fallback(p_hit, p_seq, amp_seq, amp_start_int, std_db_seqs.get(gid, ""), int(fr["amp_end"]) if fr else None)
            standard_genomes.append({
                "genome_id": gid, "lineage": m.get("lineage", gid), "date": m.get("date", ""),
                "country": m.get("country", ""), "region": m.get("region", ""),
                "f_mismatch": fr["f_mismatch"] if fr else None,
                "p_mismatch": p_hit["mismatch"] if p_hit else None,
                "r_mismatch": fr["r_mismatch"] if fr else None,
                "f_3p_ok": fr["f_3p_ok"] if fr else None,
                "r_3p_ok": fr["r_3p_ok"] if fr else None,
                "amplicon_start": fr["amp_start"] if fr else None,
                "amplicon_end": fr["amp_end"] if fr else None,
                "probe_start": p_hit["start"] if p_hit else None,
                "probe_end": p_hit["end"] if p_hit else None,
                "amplicon_seq": amp_seq,
                "f_tm": fr["f_tm"] if fr else None, "f_dg": fr["f_dg"] if fr else None,
                "p_tm": p_hit["tm"] if p_hit else None, "p_dg": p_hit["dg"] if p_hit else None,
                "r_tm": fr["r_tm"] if fr else None, "r_dg": fr["r_dg"] if fr else None,
            })
        std_validated = {g["genome_id"] for g in standard_genomes if _fr_ok(g, primer_mismatch_pass) and _probe_ok(g, probe_mismatch)}
        std_perfect_primer = any(_fr_ok(g, 0) for g in standard_genomes)
        standard = {
            "ok": len(std_validated) > 0,
            "matched": sorted(std_validated),
            "perfect_primer": std_perfect_primer,
        }

    report("참조유전체 결합위치 조회 중...")
    binding = []
    gene_info = {"target_gene": "NA", "amplicon_start": "", "amplicon_end": "", "amplicon_size": "", "chrom": ""}
    if ref_fasta.exists():
        binding = refgenome_binding(mfeprimer_bin, ref_fasta, [("F", f_id, f_seq), ("P", p_id, p_seq), ("R", r_id, r_seq)], out_dir)
        found = [h for h in binding if h["found"]]
        if found:
            starts = [h["start"] for h in found]
            ends = [h["end"] for h in found]
            genes = load_genes(genes_tsv)
            gene_hits = genes_overlapping(genes.get(seg, []), min(starts), max(ends))
            gene_info = {
                "target_gene": "/".join(gene_hits) if gene_hits else "NA",
                "amplicon_start": min(starts), "amplicon_end": max(ends),
                "amplicon_size": max(ends) - min(starts),
                "chrom": found[0]["chrom"],
            }

    report("프라이머/프로브 이차구조(dimer/hairpin) 검사 중...")
    dh_dir = out_dir / "dimer_hairpin"
    dh_dir.mkdir(parents=True, exist_ok=True)
    dh_fasta = dh_dir / "oligos.fasta"
    dh_records = [(role, s) for role, s in (("F", f_seq), ("P", p_seq), ("R", r_seq)) if s]
    dh_fasta.write_text("".join(f">{role}\n{s}\n" for role, s in dh_records), encoding="utf-8")
    dimer_report_txt = run_mfeprimer(mfeprimer_bin, "dimer", dh_fasta, dh_dir / "dimer.txt",
                                      extra_args=["--cpu", str(cpu)] if cpu else None)
    hairpin_report_txt = run_mfeprimer(mfeprimer_bin, "hairpin", dh_fasta, dh_dir / "hairpin.txt")
    dimer_hairpin = {
        "dimer": {
            pair: {"worst_dg": v["worst"]["dg"], "worst_score": v["worst"]["score"],
                   "align": v["worst"]["align"], "flagged_count": v["flagged"]}
            for pair, v in summarize_dimers(parse_dimer_report(dimer_report_txt)).items()
        },
        "hairpin": {
            role: {"worst_dg": v["worst"]["dg"], "worst_tm": v["worst"]["tm"],
                   "worst_score": v["worst"]["score"], "align": v["worst"]["align"],
                   "flagged_count": v["flagged"]}
            for role, v in summarize_hairpins(parse_hairpin_report(hairpin_report_txt)).items()
        },
        "dg_flag_threshold": DG_FLAG_THRESHOLD,
    }

    report("전체 게놈 정렬(MSA) 생성 중 (mafft)...")
    msa = build_msa(genomes, f_seq, p_seq, r_seq, out_dir / "msa")

    return {
        "set_id": name, "segment": seg.upper(), "source": "user",
        "F_id": f_id, "F_seq": f_seq, "R_id": r_id, "R_seq": r_seq, "P_id": p_id, "P_seq": p_seq,
        "summary": summary, "standard": standard, "gene": gene_info,
        "binding": binding, "genomes": genomes, "msa": msa,
        "standard_genomes": standard_genomes, "dimer_hairpin": dimer_hairpin,
    }
