#!/usr/bin/env python3
"""Run pixy (https://github.com/ksamuk/pixy) --stats pi on the all-sites VCF built by
build_diversity_vcf.py, then parse its windowed per-lineage output into one compact JSON
the webapp can render directly as a position-vs-diversity chart (mirroring Nextstrain's own
per-position "entropy" panel, but using pixy's unbiased pi estimator per lineage instead of
raw per-site Shannon entropy across the whole population).

Runs pixy TWICE: once with the real per-lineage populations file (the existing per-lineage
curves), and once with every sample pooled into a single synthetic population ("ALL_LASV") --
this second run is NOT the same as averaging the per-lineage curves (which only counts
within-lineage pairs); pooling captures inter-lineage divergence too, giving a genuine
"across the whole filtered population" baseline to plot alongside the per-lineage lines."""
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

ALL_LABEL = None   # name of the pooled all-genomes population; set from the profile in main()


def find_pixy():
    found = shutil.which("pixy")
    if found:
        return [found]
    # pixy 2.x ships as a plain Python package with a console-script entry point; if that
    # entry point isn't on PATH (e.g. installed into a different interpreter than the one
    # running this script), fall back to invoking its CLI module directly.
    for cand in sorted(Path("/BI/conda/envs").glob("*/bin/pixy")):
        return [str(cand)]
    return [sys.executable, "-m", "pixy.__main__"]


def run_pixy_pi(pixy_cmd_base, vcf_path, pops_path, window_size, out_dir, prefix):
    pixy_cmd = pixy_cmd_base + [
        "--stats", "pi",
        "--vcf", str(vcf_path),
        "--populations", str(pops_path),
        "--window_size", str(window_size),
        "--output_folder", str(out_dir),
        "--output_prefix", prefix,
        # pixy defaults to biallelic-SNPs-and-invariant-sites-only, silently dropping any
        # site with 2+ ALT alleles from the pi calculation (no_sites=0 for that site/window)
        # -- confirmed by direct testing: in a population spanning 7 distinct LASV lineages,
        # a genuinely diverse site very often has more than one lineage carrying its own
        # distinct substitution, so multiallelic sites are the NORM here, not a rare edge
        # case, and dropping them was silently discarding most of the real diversity signal
        # (74% of sites excluded in one spot-checked 165bp window). pixy's own pi estimator
        # (pixy/calc.py) is explicitly correct at multiallelic sites -- this flag only
        # controls whether they're included, not whether the math handles them.
        "--include_multiallelic_snps",
    ]
    print(f"[pixy] running: {' '.join(pixy_cmd)}", file=sys.stderr)
    result = subprocess.run(pixy_cmd, capture_output=True, text=True)
    print(result.stdout, file=sys.stderr)
    print(result.stderr, file=sys.stderr)
    if result.returncode != 0:
        sys.exit(f"[error] pixy failed (returncode={result.returncode})")
    return out_dir / f"{prefix}_pi.txt"


def parse_pi_txt(pi_path):
    """Returns (by_pop windows dict, genome-wide summary dict) from one pixy *_pi.txt."""
    by_pop = {}
    totals = {}  # pop -> [sum of avg_pi*no_sites, sum of no_sites]
    with open(pi_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pop = row["pop"]
            no_sites = int(row["no_sites"]) if row["no_sites"] not in ("", "NA") else 0
            avg_pi = row["avg_pi"]
            avg_pi = float(avg_pi) if avg_pi not in ("", "NA") else None
            by_pop.setdefault(pop, []).append({
                "start": int(row["window_pos_1"]), "end": int(row["window_pos_2"]),
                "pi": avg_pi, "no_sites": no_sites,
            })
            if avg_pi is not None:
                t = totals.setdefault(pop, [0.0, 0])
                t[0] += avg_pi * no_sites
                t[1] += no_sites
    summary = {pop: (round(sums[0] / sums[1], 5) if sums[1] else None) for pop, sums in totals.items()}
    return by_pop, summary


def write_all_populations_file(pops_path, out_path):
    """Every sample from the real (per-lineage) populations file, remapped to one shared
    synthetic population label -- so pixy computes a single pooled-population pi curve."""
    with open(pops_path, encoding="utf-8") as fh, open(out_path, "w", encoding="utf-8") as out_fh:
        for line in fh:
            sample = line.rstrip("\n").split("\t")[0]
            if sample:
                out_fh.write(f"{sample}\t{ALL_LABEL}\n")


def split_populations(pops_path, min_n, per_lineage_path):
    """Writes the per-lineage populations file pixy gets for the per-lineage curves, leaving out
    "NA" and any lineage with fewer than `min_n` genomes -- pi from a handful of sequences is
    mostly sampling noise (n=2 can only ever see one pair), so those lineages are dropped from
    the per-lineage curves (their genomes still count in the pooled ALL_LASV curve, which uses
    the untouched populations file). Returns (kept {lineage: n}, hidden {lineage: n})."""
    rows = []
    with open(pops_path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0]:
                rows.append((parts[0], parts[1]))
    counts = {}
    for _, lin in rows:
        counts[lin] = counts.get(lin, 0) + 1
    kept = {l: n for l, n in counts.items() if l != "NA" and n >= min_n}
    hidden = {l: n for l, n in counts.items() if l not in kept}
    with open(per_lineage_path, "w", encoding="utf-8") as out:
        for sid, lin in rows:
            if lin in kept:
                out.write(f"{sid}\t{lin}\n")
    return kept, hidden, len(rows)


def main():
    vcf_path = Path(snakemake.input.vcf)
    pops_path = Path(snakemake.input.populations)
    window_size = int(snakemake.params.window_size)
    global ALL_LABEL
    ALL_LABEL = snakemake.params.pooled_label
    min_n = int(snakemake.params.min_lineage_n)
    out_json = Path(snakemake.output.json)
    out_dir = out_json.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    seg = snakemake.wildcards.seg
    pixy_cmd_base = find_pixy()

    per_lineage_pops = out_dir / f"populations_lineage_{seg}.txt"
    kept, hidden, n_total = split_populations(pops_path, min_n, per_lineage_pops)
    lineage_pi_path = run_pixy_pi(pixy_cmd_base, vcf_path, per_lineage_pops, window_size, out_dir, f"pixy_{seg}")
    by_lineage, summary = parse_pi_txt(lineage_pi_path)
    per_lineage_pops.unlink(missing_ok=True)

    all_pops_path = out_dir / f"populations_all_{seg}.txt"
    write_all_populations_file(pops_path, all_pops_path)
    all_pi_path = run_pixy_pi(pixy_cmd_base, vcf_path, all_pops_path, window_size, out_dir, f"pixy_{seg}_all")
    by_all, summary_all = parse_pi_txt(all_pi_path)
    by_lineage.update(by_all)  # adds the single ALL_LASV key
    summary.update(summary_all)
    all_pops_path.unlink(missing_ok=True)

    qc = json.loads(Path(snakemake.input.qc).read_text(encoding="utf-8"))
    qc.pop("per_genome", None)  # the per-genome table stays in qc.json; the chart only needs the summary
    payload = {
        "segment": seg,
        "window_size": window_size,
        "all_label": ALL_LABEL,
        "windows": by_lineage,
        "genome_wide_pi": summary,
        "lineage_n": {**kept, ALL_LABEL: n_total},
        "lineages_hidden_small_n": hidden,
        "min_lineage_n": min_n,
        "qc": qc,
    }
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"[pixy] wrote {out_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
