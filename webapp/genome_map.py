"""Data of the genome map page: one reference genome with its genes, the standard material placed on it by alignment, and where
the primer sets bind. Results that already exist are read, the standard material is placed with a quick blastn; the analysis itself is untouched.

* genes       config/refgenome_genes.tsv                 (segment, accession, genome_len, gene, start, end, strand)
* standards   placed with blastn (local alignment, both strands) of config/kdca_standard_<seg>.fasta against the reference. The pipeline's own
              mafft placement (output/standards/<seg>/representativeness.json) is kept as the fallback and for comparison: mafft --keeplength
              can scatter a divergent standard over the wrong stretch (GC2 of Lassa S: 99.8% identical to 2603-3402, mafft put it at 1423-2400)
* sets        output/primers/<set>/refgenome_binding/binding.tsv (role, found, start, end, strand, mismatches)
"""
import csv
import json
import subprocess
import sys
from pathlib import Path


def _read_tsv(path):
    p = Path(path)
    if not p.exists():
        return []
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def find_blastn(root):
    """blastn (see scripts/standard_placement.py, which the pipeline uses too); None when there is none."""
    sys.path.insert(0, str(Path(root) / "scripts"))
    from standard_placement import find_blastn as _find  # noqa: PLC0415
    return _find(root)


def merge_hsps(rows, qlens, min_identity=80.0, min_cover=0.8):
    """blastn rows (qseqid, sstart, send, sstrand, qstart, qend, pident, length, evalue, bitscore) -> {query: placement}. Per query the best-scoring
    strand is taken, the HSPs on it that lie on the same diagonal band as the best one are merged: reference span, share of the query that
    aligned and the identity over the aligned bases. `reliable` = identity and aligned share both reach the thresholds."""
    by_q = {}
    for r in rows:
        q = r[0]
        by_q.setdefault(q, []).append({"q": q, "s0": int(r[1]), "s1": int(r[2]), "strand": r[3], "q0": int(r[4]), "q1": int(r[5]), "pid": float(r[6]),
                                       "len": int(r[7]), "bits": float(r[9])})
    out = {}
    for q, hs in by_q.items():
        best = max(hs, key=lambda h: h["bits"])
        diag = lambda h: (min(h["s0"], h["s1"]) - h["q0"]) if h["strand"] == "plus" else (min(h["s0"], h["s1"]) + h["q0"])   # noqa: E731
        keep = [h for h in hs if h["strand"] == best["strand"] and abs(diag(h) - diag(best)) <= 60]
        lo = min(min(h["s0"], h["s1"]) for h in keep)
        hi = max(max(h["s0"], h["s1"]) for h in keep)
        covered = sorted({i for h in keep for i in range(h["q0"], h["q1"] + 1)})
        qlen = qlens.get(q) or max(h["q1"] for h in keep)
        ident = sum(h["pid"] * h["len"] for h in keep) / max(1, sum(h["len"] for h in keep))
        share = len(covered) / qlen if qlen else 0.0
        out[q] = {"start": lo, "end": hi, "reversed": best["strand"] == "minus", "identity": round(ident / 100.0, 4), "mapped_fraction": round(share, 3),
                  "reliable": ident >= min_identity and share >= min_cover}
    return out


def blast_standards(root, seg):
    """{standard: placement} from blastn of the registered standard material against the reference, or {} when blastn / the files are missing."""
    root = Path(root)
    std, ref = root / "config" / f"kdca_standard_{seg}.fasta", root / "data" / "refgenome" / f"{seg}_ref.fasta"
    exe = find_blastn(root)
    if not (exe and std.exists() and ref.exists()):
        return {}
    qlens, name = {}, None
    for ln in std.read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            qlens[name] = 0
        elif name:
            qlens[name] += len(ln.strip())
    try:
        res = subprocess.run([exe, "-query", str(std), "-subject", str(ref), "-task", "blastn", "-outfmt", "6 qseqid sstart send sstrand qstart qend pident length evalue bitscore"],
                             capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return {}
    if res.returncode:
        return {}
    return merge_hsps([ln.split("\t") for ln in res.stdout.splitlines() if ln.strip()], qlens)


_BLAST_CACHE = {}


def build(root, seg, reference_accession, sets, groups, active_standards=None):
    """sets: [{"set_id", "name", ...}] of this segment; groups: [{"group_id", "label", "members": [set ids]}] of this segment.
    active_standards: names of the standard materials switched on in the standards manifest (분석진행); the others are left out (None = all).
    -> {"segment", "reference", "length", "genes", "standards", "sets", "groups"}"""
    root = Path(root)
    genes = [r for r in _read_tsv(root / "config" / "refgenome_genes.tsv") if (r.get("segment") or "").lower() == seg.lower()]
    length = max([_int(r.get("genome_len")) or 0 for r in genes] + [0])
    ref = (genes[0].get("accession") if genes else "") or reference_accession or ""
    out_genes = [{"gene": r["gene"], "start": _int(r["start"]), "end": _int(r["end"]), "strand": r.get("strand") or "+"}
                 for r in genes if _int(r.get("start")) and _int(r.get("end"))]
    standards = []
    rep = root / "output" / "standards" / seg.lower() / "representativeness.json"
    data = {}
    if rep.exists():
        try:
            data = json.loads(rep.read_text(encoding="utf-8"))
        except ValueError:
            data = {}
    per = data.get("per_standard") or {}
    stamp = tuple(p.stat().st_mtime_ns if p.exists() else None for p in (root / "config" / f"kdca_standard_{seg}.fasta", root / "data" / "refgenome" / f"{seg}_ref.fasta"))
    key = (str(root), seg, stamp)
    if key not in _BLAST_CACHE:
        _BLAST_CACHE.clear()
        _BLAST_CACHE[key] = blast_standards(root, seg.lower())
    blast = _BLAST_CACHE[key]
    names = [n for n in dict.fromkeys(list((data.get("standards") or {}).keys()) + list(blast)) if active_standards is None or n in active_standards]
    for name in names:
        s = (data.get("standards") or {}).get(name) or {}
        piped = {"start": _int(s.get("ref_start")), "end": _int(s.get("ref_end")), "identity": s.get("identity_to_reference"), "mapped_fraction": s.get("mapped_fraction"),
                 "reversed": bool(s.get("reversed"))} if s.get("mapped") and _int(s.get("ref_start")) else None
        b = blast.get(name)
        if b:
            placed = dict(b, method="blastn")
        elif piped:
            ok = (piped["identity"] or 0) >= 0.8 and (piped["mapped_fraction"] or 0) >= 0.8          # same bar as blastn's: a scattered mafft placement is not trusted
            placed = dict(piped, method="mafft", reliable=ok)
        else:
            placed = None
        entry = {"name": name, "mapped": bool(placed), "start": placed["start"] if placed else None, "end": placed["end"] if placed else None,
                 "length": _int(s.get("length")), "mapped_fraction": placed["mapped_fraction"] if placed else None,
                 "identity": placed["identity"] if placed else None, "reversed": bool(placed and placed["reversed"]), "location": s.get("location"),
                 "matched_pct": (per.get(name) or {}).get("pct"), "method": placed["method"] if placed else None, "reliable": bool(placed and placed["reliable"]),
                 "pipeline_start": piped["start"] if piped else None, "pipeline_end": piped["end"] if piped else None}
        # the pipeline's mafft placement disagrees with the one shown (more than 100 bp off): say so
        entry["pipeline_differs"] = bool(placed and piped and placed["method"] == "blastn" and (abs(piped["start"] - placed["start"]) > 100 or abs(piped["end"] - placed["end"]) > 100))
        standards.append(entry)
    out_sets = []
    for st in sets:
        hits = {}
        for r in _read_tsv(root / "output" / "primers" / st["set_id"] / "refgenome_binding" / "binding.tsv"):
            role = r.get("role")
            if role in ("F", "P", "R"):
                found = str(r.get("found")) in ("1", "true", "True")
                hits[role] = {"found": found, "start": _int(r.get("start")) if found else None, "end": _int(r.get("end")) if found else None,
                              "strand": r.get("strand") or "", "mismatches": _int(r.get("mismatches"))}
        located = [h for h in hits.values() if h["found"] and h["start"] and h["end"]]
        amp = [min(min(h["start"], h["end"]) for h in located), max(max(h["start"], h["end"]) for h in located)] if located else None
        out_sets.append({"set_id": st["set_id"], "name": st.get("name") or st["set_id"], "hits": hits, "amplicon": amp})
    return {"segment": seg.upper(), "reference": ref, "length": length, "genes": out_genes, "standards": standards,
            "sets": out_sets, "groups": [{"group_id": g["group_id"], "label": g.get("label") or g["group_id"], "members": g.get("members") or []} for g in groups]}
