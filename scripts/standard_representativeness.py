#!/usr/bin/env python3
"""How well does the KDCA reference-standard material (표준물질) represent the real LASV
population? -- per segment, written to output/standards/<seg>/representativeness.json.

The wet-lab positive control is one (or a few) short gene snippets. If most real genomes
differ from it in the region the assay looks at, a "positive" on the control says little about
whether the assay will detect field viruses. This script quantifies that, in three parts:

1. Standard-wide match (layer 1). Each standard is aligned to the reference genome (mafft
   --keeplength --add, same reference-anchored mode the diversity VCF uses), which gives it
   reference coordinates. Every QC-passed population genome's base at those positions is then
   read from the pixy all-sites VCF (output/diversity/<seg>/allsites.vcf.gz -- so "trusted
   genomes" means exactly the set that passed the diversity QC, reversed genomes already
   corrected) and compared with the standard: a genome "matches" when identity over the
   comparable positions >= standard_rep_identity. A genome that doesn't cover enough of the
   standard's span (< standard_rep_min_coverage) can't be judged and is left out of the
   denominator (reported separately).

2. Assay-site match per primer set (layer 2). If the standard fully covers the set's F/R/(P)
   binding sites on the reference, the genome's bases at exactly those sites are compared with the
   standard's: a genome matches when F and R differ by <= primer_mismatch_pass bases and the
   probe by <= probe_mismatch -- i.e. the assay would "see" this genome like it sees the
   standard. (Tolerances: standard_rep_site_tol_primer / _probe -- primer sites of a
   divergent virus almost never match a single sequence exactly, so the pipeline's own
   mm<=1 / probe mm=0 pass criteria would be far too strict for a genome-vs-genome comparison.) A set whose sites aren't inside any standard is reported as not applicable.

3. Diversity of the location. Mean pixy pi (ALL_LASV curve) over the standard's span and over
   each oligo site, ranked against same-width windows elsewhere on the segment (percentile
   <= 25 -> conserved, >= 75 -> variable), plus how often the standard carries the population's
   majority base at the oligo sites.

Grades use the union over all standards that cover the region: >= standard_rep_good (90%)
-> good, >= standard_rep_critical (50%) -> review, below -> critical.
"""
import csv
import gzip
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_diversity_vcf import read_fasta, resolve_add_flips  # noqa: E402
from mfe_common import find_mafft  # noqa: E402

BASES = set("ACGT")
MIN_SITE_CALLED = 0.9   # a genome must have a real base at >= 90% of the oligo-site positions to be judged


def grade_for(pct, good, critical):
    if pct is None:
        return "na"
    return "good" if pct >= good else ("review" if pct >= critical else "critical")


def align_standards_to_reference(mafft_bin, ref_fasta, standards, threads=2):
    """{standard_id: aligned row in reference coordinates (len == ref length, '-' outside)}"""
    with tempfile.TemporaryDirectory() as tmp:
        add = Path(tmp) / "standards.fasta"
        add.write_text("".join(f">{k}\n{v}\n" for k, v in standards.items()), encoding="utf-8")
        cmd = [mafft_bin, "--keeplength", "--adjustdirection", "--add", str(add), "--quiet",
               "--thread", str(threads), str(ref_fasta)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            sys.exit(f"[error] mafft failed: {res.stderr}")
    rows = {}
    header, chunks = None, []
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            if header is not None:
                rows[header] = "".join(chunks)
            header, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if header is not None:
        rows[header] = "".join(chunks)
    rows, reversed_ids = resolve_add_flips(rows)
    return rows, reversed_ids


def load_vcf_bases(vcf_path, positions):
    """samples, {pos: [base-or-None per sample]} for the requested reference positions only."""
    samples, bases = [], {}
    with gzip.open(vcf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            pos = int(f[1])
            if pos not in positions:
                continue
            alleles = [f[3]] + ([] if f[4] == "." else f[4].split(","))
            bases[pos] = [None if g == "." else alleles[int(g)] for g in f[9:]]
    return samples, bases


def main_block(cols, std_len, slack=1.25):
    """mafft sometimes scatters a few stray bases of a short, divergent standard across the
    reference (a standard's real home is one contiguous stretch). Keeps only the densest
    window of ~std_len*slack reference positions and drops everything outside it."""
    if not cols:
        return {}
    positions = sorted(cols)
    width = int(std_len * slack)
    best, best_i, j = 0, 0, 0
    for i, p in enumerate(positions):
        while positions[j] < p - width + 1:
            j += 1
        if i - j + 1 > best:
            best, best_i = i - j + 1, j
    keep = positions[best_i:best_i + best]
    return {p: cols[p] for p in keep}


class PiTrack:
    """Prefix sums over the ALL_LASV pi curve (weighted by pixy's no_sites)."""

    def __init__(self, windows, genome_len):
        num = [0.0] * (genome_len + 2)
        den = [0.0] * (genome_len + 2)
        for w in windows:
            if w.get("pi") is None or not w.get("no_sites"):
                continue
            span = max(1, w["end"] - w["start"] + 1)
            for p in range(max(1, w["start"]), min(genome_len, w["end"]) + 1):
                num[p] += w["pi"] * w["no_sites"] / span
                den[p] += w["no_sites"] / span
        self.n = genome_len
        self.cn = [0.0]
        self.cd = [0.0]
        for p in range(1, genome_len + 1):
            self.cn.append(self.cn[-1] + num[p])
            self.cd.append(self.cd[-1] + den[p])
        self._dist = {}

    def mean(self, lo, hi):
        lo, hi = max(1, lo), min(self.n, hi)
        d = self.cd[hi] - self.cd[lo - 1]
        return (self.cn[hi] - self.cn[lo - 1]) / d if d > 0 else None

    def _distribution(self, width):
        if width not in self._dist:
            step = max(1, width // 4)
            vals = []
            for lo in range(1, self.n - width + 2, step):
                m = self.mean(lo, lo + width - 1)
                if m is not None:
                    vals.append(m)
            self._dist[width] = sorted(vals)
        return self._dist[width]

    def describe(self, lo, hi):
        """{mean_pi, percentile (0-100), location: conserved|moderate|variable} for [lo, hi]."""
        m = self.mean(lo, hi)
        if m is None:
            return {"mean_pi": None, "percentile": None, "location": "unknown"}
        dist = self._distribution(hi - lo + 1)
        pct = 100.0 * sum(1 for v in dist if v < m) / len(dist) if dist else 50.0
        loc = "conserved" if pct <= 25 else ("variable" if pct >= 75 else "moderate")
        return {"mean_pi": round(m, 5), "percentile": round(pct, 1), "location": loc}


def summarize(match_by_genome, lineage_of, good, critical):
    """match_by_genome: {genome: True/False/None(None = can't be judged)}."""
    judged = {g: v for g, v in match_by_genome.items() if v is not None}
    matched = sum(1 for v in judged.values() if v)
    by_lineage = {}
    for g, v in judged.items():
        b = by_lineage.setdefault(lineage_of.get(g, "NA"), {"n": 0, "matched": 0})
        b["n"] += 1
        b["matched"] += 1 if v else 0
    pct = round(100.0 * matched / len(judged), 2) if judged else None
    return {
        "assessable": len(judged), "not_assessable": len(match_by_genome) - len(judged),
        "matched": matched, "pct": pct, "grade": grade_for(pct, good, critical),
        "by_lineage": by_lineage,
        "unmatched": sorted(g for g, v in judged.items() if not v),
    }


def compute(seg, standards, ref_id, ref_seq, aligned_rows, samples, bases, lineage_of, pi_windows,
            bindings_by_set, params, reversed_ids=None):
    good, critical = params["good"], params["critical"]
    min_ident, min_cov = params["identity"], params["min_coverage"]
    tol_primer, tol_probe = params["site_tol_primer"], params["site_tol_probe"]
    genome_len = len(ref_seq)
    pi = PiTrack(pi_windows, genome_len)
    n_s = len(samples)

    # --- standards on reference coordinates
    std_info, std_cols = {}, {}   # std_cols[id] = {pos: base} for called positions
    for sid in standards:
        row = aligned_rows.get(sid, "").upper()
        raw = {i + 1: c for i, c in enumerate(row) if c in BASES}
        cols = main_block(raw, len(standards[sid]))
        std_cols[sid] = cols
        if not cols:
            std_info[sid] = {"mapped": False}
            continue
        lo, hi = min(cols), max(cols)
        ref_same = sum(1 for p, c in cols.items() if ref_seq[p - 1] == c)
        std_info[sid] = {
            "mapped": True, "ref_start": lo, "ref_end": hi, "mapped_bases": len(cols),
            "length": len(standards[sid]), "mapped_fraction": round(len(cols) / len(standards[sid]), 3),
            "reversed": sid in (reversed_ids or []),
            "identity_to_reference": round(ref_same / len(cols), 4),
            **{k: v for k, v in pi.describe(lo, hi).items()},
        }

    # --- layer 1: per standard, per genome
    layer1_ident = {sid: {} for sid in standards}   # sid -> genome -> (coverage, identity)
    for sid, cols in std_cols.items():
        if not cols:
            continue
        for gi in range(n_s):
            called = same = 0
            for p, c in cols.items():
                col = bases.get(p)
                b = col[gi] if col else None
                if b is not None:
                    called += 1
                    same += (b == c)
            layer1_ident[sid][samples[gi]] = (called / len(cols), same / called if called else 0.0)

    def l1_match(sid, g):
        cov, ident = layer1_ident[sid][g]
        return None if cov < min_cov else ident >= min_ident

    per_standard = {}
    for sid in standards:
        if not std_info[sid]["mapped"]:
            per_standard[sid] = {"assessable": 0, "matched": 0, "pct": None, "grade": "na"}
            continue
        per_standard[sid] = summarize({g: l1_match(sid, g) for g in samples}, lineage_of, good, critical)
        idents = sorted(i for (c, i) in layer1_ident[sid].values() if c >= min_cov)
        per_standard[sid]["identity_median"] = round(idents[len(idents) // 2], 4) if idents else None
        per_standard[sid].pop("unmatched", None)
    union_match, best = {}, {}
    for g in samples:
        vals = [l1_match(s, g) for s in standards if std_info[s]["mapped"]]
        judged = [v for v in vals if v is not None]
        union_match[g] = None if not judged else any(judged)
        cands = [(layer1_ident[s][g][1], s) for s in standards
                 if std_info[s]["mapped"] and layer1_ident[s][g][0] >= min_cov]
        best[g] = max(cands) if cands else None
    union = summarize(union_match, lineage_of, good, critical)
    union["sensitivity"] = []
    for thr in params["identity_sensitivity"]:
        mu = {}
        for g in samples:
            vals = [(layer1_ident[s][g][1] >= thr) for s in standards
                    if std_info[s]["mapped"] and layer1_ident[s][g][0] >= min_cov]
            mu[g] = any(vals) if vals else None
        su = summarize(mu, lineage_of, good, critical)
        union["sensitivity"].append({"identity": thr, "pct": su["pct"], "grade": su["grade"]})
    union["unmatched_detail"] = [
        {"genome_id": g, "lineage": lineage_of.get(g, "NA"),
         "best_standard": best[g][1] if best[g] else None,
         "best_identity": round(best[g][0], 4) if best[g] else None}
        for g in union.pop("unmatched")
    ]

    # --- layer 2: per primer set
    majority = {}
    def maj(p):
        if p not in majority:
            col = bases.get(p) or []
            cnt = Counter(b for b in col if b is not None)
            majority[p] = cnt.most_common(1)[0][0] if cnt else None
        return majority[p]

    sets_out = {}
    for set_id, hits in bindings_by_set.items():
        found = {h["role"]: h for h in hits if h.get("found") == "1" and h.get("start") and h.get("end")}
        if "F" not in found or "R" not in found:
            sets_out[set_id] = {"applicable": False, "reason": "no_binding_on_reference", "grade": "na"}
            continue
        sites = {r: (int(h["start"]), int(h["end"])) for r, h in found.items()}
        site_pi = {r: dict(pi.describe(lo, hi), start=lo, end=hi) for r, (lo, hi) in sites.items()}
        covering = [sid for sid in standards if std_info[sid]["mapped"] and all(
            p in std_cols[sid] for lo, hi in sites.values() for p in range(lo, hi + 1))]
        entry = {"applicable": bool(covering), "sites": site_pi, "covering_standards": covering}
        if not covering:
            entry.update(reason="sites_outside_standard", grade="na")
            sets_out[set_id] = entry
            continue
        per_std, match_union = {}, {}
        site_diffs = {sid: {} for sid in covering}
        for sid in covering:
            ms = {}
            for gi, g in enumerate(samples):
                tot = called = 0
                diffs = {r: 0 for r in sites}
                for r, (lo, hi) in sites.items():
                    for p in range(lo, hi + 1):
                        tot += 1
                        col = bases.get(p)
                        b = col[gi] if col else None
                        if b is not None:
                            called += 1
                            diffs[r] += (b != std_cols[sid][p])
                if called / tot < MIN_SITE_CALLED:
                    ms[g] = None
                else:
                    ms[g] = (diffs["F"] <= tol_primer and diffs["R"] <= tol_primer and diffs.get("P", 0) <= tol_probe)
                    site_diffs[sid][g] = diffs
            s = summarize(ms, lineage_of, good, critical)
            s.pop("unmatched", None)
            s.pop("by_lineage", None)
            s["majority_agreement"] = {
                r: round(sum(1 for p in range(lo, hi + 1) if maj(p) == std_cols[sid][p]) / (hi - lo + 1), 3)
                for r, (lo, hi) in sites.items()}
            per_std[sid] = s
            for g, v in ms.items():
                if v is None:
                    match_union.setdefault(g, None)
                else:
                    match_union[g] = bool(match_union.get(g)) or v
        u = summarize(match_union, lineage_of, good, critical)
        # how the verdict moves with the tolerance (so a viewer sees it isn't a knife-edge call)
        sens = []
        for tp, tq in params["site_sensitivity"]:
            mu = {}
            for g in samples:
                vals = [(d[g]["F"] <= tp and d[g]["R"] <= tp and d[g].get("P", 0) <= tq)
                        for d in site_diffs.values() if g in d]
                mu[g] = any(vals) if vals else None
            su = summarize(mu, lineage_of, good, critical)
            sens.append({"primer": tp, "probe": tq, "pct": su["pct"], "grade": su["grade"]})
        entry.update(per_standard=per_std, union=u, grade=u["grade"], pct=u["pct"], sensitivity=sens)
        sets_out[set_id] = entry

    return {
        "segment": seg, "reference": ref_id, "genome_count": n_s,
        "params": {k: params[k] for k in ("identity", "min_coverage", "good", "critical", "site_tol_primer", "site_tol_probe")},
        "standards": std_info, "per_standard": per_standard, "union": union,
        "genome_lineage": {g: lineage_of.get(g, "NA") for g in samples},
        "sets": sets_out,
    }


def read_bindings(paths):
    out = {}
    for p in paths:
        with open(p, newline="", encoding="utf-8") as fh:
            out[Path(p).parent.parent.name] = list(csv.DictReader(fh, delimiter="\t"))
    return out


def main():
    seg = snakemake.wildcards.seg
    standards = read_fasta(snakemake.input.standards)
    ref = read_fasta(snakemake.input.ref_fasta)
    ref_id, ref_seq = next(iter(ref.items()))
    ref_seq = ref_seq.upper()
    rows, reversed_ids = align_standards_to_reference(find_mafft(), snakemake.input.ref_fasta, standards)
    lineage_of = {}
    with open(snakemake.input.populations, encoding="utf-8") as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) >= 2:
                lineage_of[a[0]] = a[1]
    positions = {p for row in rows.values() for p, c in enumerate(row.upper(), 1) if c in BASES}
    samples, bases = load_vcf_bases(snakemake.input.vcf, positions)
    pixy = json.loads(Path(snakemake.input.pixy).read_text(encoding="utf-8"))
    p = snakemake.params
    params = {"identity": float(p.identity), "min_coverage": float(p.min_coverage),
              "good": float(p.good), "critical": float(p.critical),
              "site_tol_primer": int(p.site_tol_primer), "site_tol_probe": int(p.site_tol_probe),
              "identity_sensitivity": [0.95, 0.90, 0.85, 0.80],
              "site_sensitivity": [(1, 0), (2, 1), (3, 2), (4, 3)]}
    # rows holds the standards only (mafft output also has the reference row -- harmless)
    std_rows = {k: v for k, v in rows.items() if k in standards}
    result = compute(seg, standards, ref_id, ref_seq, std_rows, samples, bases, lineage_of,
                     pixy["windows"].get(pixy.get("all_label", "ALL_LASV"), []),
                     read_bindings(snakemake.input.bindings), params, reversed_ids)
    out = Path(snakemake.output.json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
    print(f"[standard_rep] {seg}: union {result['union']['matched']}/{result['union']['assessable']} "
          f"({result['union']['pct']}%) grade={result['union']['grade']}", file=sys.stderr)


if __name__ == "__main__":
    main()
