"""Where a standard material sits on the reference genome, decided with blastn (local alignment, both strands).

The reference-anchored mafft alignment (`mafft --keeplength --add`) is fine for a standard that is close to the reference, but it can scatter a
divergent standard over the wrong stretch of the genome: Lassa S standard GC2 (800 nt) is 99.8% identical to reference positions 2603-3402 on the
opposite strand, and mafft placed 405 of its bases at 1423-2400 with 51% identity. blastn finds the one real home.

`blast_rows()` returns each standard as a row in REFERENCE coordinates (the same shape as the mafft rows, '-' outside), but only for the standards
blastn places RELIABLY (identity and aligned share both >= 80%); the caller keeps the mafft row for every other standard, so a standard that only
partly matches (e.g. a short conserved stretch) is not shrunk.
"""
import collections
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def find_blastn(root=None):
    """blastn from $BLASTN_BIN, the PATH, the pipeline's conda env or any conda env next to the running Python; None when there is none."""
    cands = [os.environ.get("BLASTN_BIN"), shutil.which("blastn")]
    try:
        from mfe_common import read_config_str  # noqa: PLC0415
        prefix = read_config_str(root or Path(__file__).resolve().parent.parent, "env_prefix")
        if prefix:
            cands.append(str(Path(prefix) / "bin" / "blastn"))
    except Exception:  # noqa: BLE001
        pass
    pref = Path(sys.prefix)                      # a conda base (<base>/envs/*) or one of its envs (<base>/envs/<env>)
    base = pref.parent.parent if pref.parent.name == "envs" else pref
    for pat in (pref / "bin" / "blastn", base / "envs" / "*" / "bin" / "blastn", base / "bin" / "blastn"):
        cands += sorted(glob.glob(str(pat)))
    return next((c for c in cands if c and Path(c).exists()), None)


def blast_rows(blastn_bin, ref_fasta, standards, min_identity=80.0, min_cover=0.8):
    """standards: {id: sequence}. -> ({id: row in reference coordinates}, [ids placed on the opposite strand]) for the reliably placed ones."""
    with __import__("tempfile").TemporaryDirectory() as tmp:
        q = Path(tmp) / "standards.fasta"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in standards.items()), encoding="utf-8")
        res = subprocess.run([blastn_bin, "-query", str(q), "-subject", str(ref_fasta), "-task", "blastn",
                              "-outfmt", "6 qseqid sstart send sstrand qstart qend pident length bitscore qseq sseq"], capture_output=True, text=True, timeout=300)
    if res.returncode:
        return {}, []
    ref_len = sum(len(ln.strip()) for ln in Path(ref_fasta).read_text(encoding="utf-8").splitlines() if not ln.startswith(">"))
    hsps = collections.defaultdict(list)
    for ln in res.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 11:
            continue
        hsps[f[0]].append({"s0": int(f[1]), "strand": f[3], "q0": int(f[4]), "q1": int(f[5]), "pid": float(f[6]), "len": int(f[7]), "bits": float(f[8]), "qseq": f[9], "sseq": f[10]})
    rows, reversed_ids = {}, []
    for sid, hs in hsps.items():
        best = max(hs, key=lambda h: h["bits"])
        diag = lambda h: (h["s0"] - h["q0"]) if h["strand"] == "plus" else (h["s0"] + h["q0"])        # noqa: E731
        keep = [h for h in hs if h["strand"] == best["strand"] and abs(diag(h) - diag(best)) <= 60]
        covered = {i for h in keep for i in range(h["q0"], h["q1"] + 1)}
        qlen = len(standards.get(sid, ""))
        ident = sum(h["pid"] * h["len"] for h in keep) / max(1, sum(h["len"] for h in keep))
        if not qlen or ident < min_identity or len(covered) / qlen < min_cover:
            continue                                                                                     # not a reliable placement: the caller keeps the mafft row
        row = ["-"] * ref_len
        for h in keep:
            pos = h["s0"]
            for qc, sc in zip(h["qseq"].upper(), h["sseq"].upper()):
                if sc != "-":                                                                            # a column of the reference
                    if qc in "ACGT" and 1 <= pos <= ref_len:
                        row[pos - 1] = _COMP[qc] if h["strand"] == "minus" else qc
                    pos += -1 if h["strand"] == "minus" else 1
        rows[sid] = "".join(row)
        if best["strand"] == "minus":
            reversed_ids.append(sid)
    return rows, reversed_ids
