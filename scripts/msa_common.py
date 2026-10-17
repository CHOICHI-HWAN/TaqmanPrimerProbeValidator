"""Shared multi-sequence-alignment builder, used by both the Snakemake `msa_alignment`
rule (scripts/build_msa.py) and the interactive webapp's live analysis
(scripts/live_analyze.py) so built-in and user-submitted sets get identical treatment.

The Forward/Probe/Reverse reference rows are built by feeding the oligo sequences into
the *same* mafft run as the genome amplicons (as extra short records) so mafft places
real, correctly gapped copies of them into the true alignment coordinates itself, rather
than us guessing their position from a representative genome.

F/R orientation is NOT auto-detected by mafft (`--adjustdirection`) -- that was tried
first and found unreliable for these short (~20nt), often highly-degenerate oligos: with
hundreds of much longer amplicons in the same run, mafft's fast distance-estimation step
can anchor a 20nt query to the wrong end of the alignment (observed for VendorD_S_PREV_2,
where the reverse primer's un-flipped form got placed on top of the forward primer's own
site). Instead each genome's own amplicon_seq is normalized (see
_normalize_amplicon_orientation) BEFORE mafft ever sees it, so the forward oligo always
matches it as given (no RC needed) and the reverse oligo always needs
reverse-complementing to match its 3' end -- true by mfeprimer's own amplicon-extraction
convention for most genomes, but not guaranteed for one deposited under the opposite
strand convention (confirmed for OQ919513/VendorD_GPC_PREVIOUS1_LASV1_P2), which is why this
is actively checked per genome rather than assumed. Only the probe's orientation is
genuinely ambiguous a priori, so both orientations are aligned and the one that actually
fits the real-genome consensus (closer to its own length, fewer mismatches) is kept."""
import statistics
import subprocess
import sys
from pathlib import Path

from mfe_common import IUPAC, find_mafft, iupac_mismatch_count, reverse_complement

REF_F_ID = "__REF_FORWARD__"
REF_P_FWD_ID = "__REF_PROBE_FWD__"
REF_P_RC_ID = "__REF_PROBE_RC__"
REF_R_ID = "__REF_REVERSE__"


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


def _normalize_amplicon_orientation(amp_seq: str, f_seq: str, r_seq: str) -> str:
    """mfeprimer's amplicon extraction assumes the genome record is deposited on the same
    strand as the reference used to design F/R, so F always sits at the amplicon's own
    start and R (reverse-complemented) at its end -- true for most genomes, but some
    accessions are deposited under the opposite strand convention, which flips this
    (confirmed for OQ919513/VendorD_GPC_PREVIOUS1_LASV1_P2: R matched the amplicon's own first
    bases directly and RC(F) matched its last bases -- exactly backwards). Feeding such a
    genome's un-normalized amplicon into the same mafft run as everyone else's
    correctly-oriented ones misaligns it (mafft has no way to know it should be
    reverse-complemented first -- see the module docstring on why mafft's own
    --adjustdirection isn't used for this). Self-detect and reverse-complement to the
    orientation that actually fits, exactly the same comparison highlightAmplicon() in
    app.js makes for the per-genome detail view, so both agree."""
    f_len, r_len = len(f_seq), len(r_seq)
    n = len(amp_seq)
    if not amp_seq or not f_len or not r_len or n < max(f_len, r_len):
        return amp_seq
    r_rc = reverse_complement(r_seq)
    direct_mm = (iupac_mismatch_count(f_seq, amp_seq[:f_len])
                 + iupac_mismatch_count(r_rc, amp_seq[max(0, n - r_len):]))
    flipped_seq = reverse_complement(amp_seq)
    flipped_mm = (iupac_mismatch_count(f_seq, flipped_seq[:f_len])
                  + iupac_mismatch_count(r_rc, flipped_seq[max(0, n - r_len):]))
    return flipped_seq if flipped_mm < direct_mm else amp_seq


def span_of(aligned_seq: str):
    """[start, end) of the non-gap run in an otherwise all-gap aligned row (used for the
    short synthetic F/P/R reference sequences, which mafft places as one contiguous
    ungapped run inside a sea of leading/trailing gap columns)."""
    start = None
    end = None
    for i, ch in enumerate(aligned_seq):
        if ch != "-":
            if start is None:
                start = i
            end = i + 1
    if start is None:
        return None
    return [start, end]


def _score_probe_orientation(aligned_seq, orig_seq, consensus):
    """Lower is better. Primary signal: how far the oligo's real aligned span is from its
    own length (a correctly-oriented short oligo aligns as one tight, ungapped block;
    a wrong-orientation one tends to get fragmented or padded with internal gaps).
    Secondary signal: IUPAC-aware mismatches against the real-genome consensus."""
    span = span_of(aligned_seq)
    if span is None:
        return None
    a, b = span
    width = b - a
    frag_penalty = abs(width - len(orig_seq))
    mism = iupac_mismatch_count(orig_seq, consensus[a:b].replace("-", ""))
    return (frag_penalty, mism)


def _run_mafft_once(mafft_bin, records, f_seq, p_seq, r_seq, out_dir: Path, tag: str):
    """records: list of (genome_id, amplicon_seq). Returns {header: aligned_seq} including
    the synthetic F/P/R reference records."""
    r_seq_rc = reverse_complement(r_seq) if r_seq else ""
    dedup_fa = out_dir / f"dedup_{tag}.fa"
    with open(dedup_fa, "w", encoding="utf-8") as fh:
        for gid, seq in records:
            fh.write(f">{gid}\n{seq}\n")
        fh.write(f">{REF_F_ID}\n{f_seq}\n")
        if p_seq:
            fh.write(f">{REF_P_FWD_ID}\n{p_seq}\n")
            fh.write(f">{REF_P_RC_ID}\n{reverse_complement(p_seq)}\n")
        fh.write(f">{REF_R_ID}\n{r_seq_rc}\n")

    aligned_fa = out_dir / f"aligned_{tag}.fa"
    with open(aligned_fa, "w", encoding="utf-8") as out_f:
        result = subprocess.run([mafft_bin, "--auto", "--quiet", str(dedup_fa)],
                                 stdout=out_f, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] mafft failed: {result.stderr}")

    aligned = {hdr: seq.upper() for hdr, seq in read_fasta(aligned_fa)}
    dedup_fa.unlink(missing_ok=True)
    aligned_fa.unlink(missing_ok=True)
    return aligned


def _find_disruptive_genomes(aligned: dict, ref_ids: set, min_private_cols=12, min_private_frac=0.08):
    """A handful of genuinely divergent/mis-assembled amplicons can force mafft to open
    huge gap columns that everyone else pads through, stretching the whole alignment into
    an unreadable wall of gaps. Detect the specific sequences responsible: for every
    column, track whether exactly one sequence has real (non-gap) content there ("private"
    to it) -- a sequence single-handedly "owning" many such columns is the one inserting
    content nobody else shares, i.e. the thing dragging the alignment out of shape.
    Returns {genome_id: private_column_count} for genomes over the threshold."""
    genome_items = [(gid, seq) for gid, seq in aligned.items() if gid not in ref_ids]
    if not genome_items:
        return {}
    n_cols = len(genome_items[0][1])
    owner = [None] * n_cols  # None=no one yet, gid=sole owner so far, False=multiple owners
    for gid, seq in genome_items:
        for col in range(n_cols):
            if seq[col] == "-":
                continue
            if owner[col] is None:
                owner[col] = gid
            elif owner[col] is not False and owner[col] != gid:
                owner[col] = False

    private_counts = {}
    for col_owner in owner:
        if col_owner and col_owner is not False:
            private_counts[col_owner] = private_counts.get(col_owner, 0) + 1

    disruptive = {}
    seq_by_id = dict(genome_items)
    for gid, cnt in private_counts.items():
        own_len = len(seq_by_id[gid].replace("-", ""))
        if cnt >= min_private_cols and (own_len == 0 or cnt >= min_private_frac * own_len):
            disruptive[gid] = cnt
    return disruptive


MAX_DISRUPTIVE_PASSES = 4  # bounded refinement loop, see build_msa


def build_msa(genome_rows, f_seq, p_seq, r_seq, out_dir: Path):
    """genome_rows: list of dicts with genome_id, lineage, date, country,
    f_mismatch, p_mismatch, r_mismatch, amplicon_seq (strings, as produced by
    build_alignment_table.py / live_analyze.py). Returns the msa.json-shaped dict.

    Iterative mafft refinement: each pass finds genomes single-handedly stretching the
    alignment with content no one else shares (see _find_disruptive_genomes), pulls them
    out, and re-aligns the rest. Removing the worst offender can reveal a second, milder
    one that wasn't detectable while the first one dominated the alignment (many mildly
    divergent sequences can collectively bloat it as badly as one bad outlier), so this
    repeats -- bounded by MAX_DISRUPTIVE_PASSES -- until a pass finds nothing left to
    remove. Every pulled-out genome is reported in `excluded` (with a reason) instead of
    silently dropped, so the UI can still list them."""
    all_rows = [g for g in genome_rows if g.get("amplicon_seq")]
    out_dir.mkdir(parents=True, exist_ok=True)

    if not all_rows:
        return {"columns": 0, "f_col": None, "p_col": None, "r_col": None, "consensus": "", "rows": [],
                "excluded": [], "ref_rows": {}}

    median_len = statistics.median(len(g["amplicon_seq"]) for g in all_rows)
    tolerance = max(20, median_len * 0.2)
    rows = [g for g in all_rows if abs(len(g["amplicon_seq"]) - median_len) <= tolerance]
    excluded = [{"genome_id": g["genome_id"], "lineage": g.get("lineage", "NA"), "reason": "length_outlier"}
                for g in all_rows if g not in rows]

    mafft_bin = find_mafft()
    ref_ids = {REF_F_ID, REF_P_FWD_ID, REF_P_RC_ID, REF_R_ID}

    aligned = {}
    for i in range(MAX_DISRUPTIVE_PASSES):
        if not rows:
            break
        records = [(g["genome_id"], _normalize_amplicon_orientation(g["amplicon_seq"], f_seq, r_seq)) for g in rows]
        aligned = _run_mafft_once(mafft_bin, records, f_seq, p_seq, r_seq, out_dir, f"pass{i}")
        disruptive = _find_disruptive_genomes(aligned, ref_ids)
        if not disruptive:
            break
        by_id = {g["genome_id"]: g for g in rows}
        for gid in disruptive:
            excluded.append({"genome_id": gid, "lineage": by_id[gid].get("lineage", "NA"),
                              "reason": "disrupts_alignment"})
        rows = [g for g in rows if g["genome_id"] not in disruptive]

    if not rows or not aligned:
        return {"columns": 0, "f_col": None, "p_col": None, "r_col": None, "consensus": "", "rows": [],
                "excluded": excluded, "ref_rows": {}}

    f_aligned = aligned.pop(REF_F_ID, None)
    p_fwd_aligned = aligned.pop(REF_P_FWD_ID, None)
    p_rc_aligned = aligned.pop(REF_P_RC_ID, None)
    r_aligned = aligned.pop(REF_R_ID, None)

    n_cols = len(next(iter(aligned.values())))
    seqs = list(aligned.values())
    consensus_chars = []
    for col in range(n_cols):
        counts = {}
        for s in seqs:
            c = s[col]
            counts[c] = counts.get(c, 0) + 1
        consensus_chars.append(max(counts.items(), key=lambda kv: kv[1])[0])
    consensus = "".join(consensus_chars)

    # F always matches as given, R was fed in already reverse-complemented -- both are
    # "+"/"-" strand labels describing that fixed, known relationship to the oligo as the
    # user wrote it (not auto-detected), so the UI can still show forward/reverse arrows.
    ref_rows = {
        "forward": {"strand": "+", "aligned": f_aligned, "col": span_of(f_aligned)} if f_aligned else None,
        "reverse": {"strand": "-", "aligned": r_aligned, "col": span_of(r_aligned)} if r_aligned else None,
    }
    if p_fwd_aligned or p_rc_aligned:
        candidates = []
        if p_fwd_aligned:
            score = _score_probe_orientation(p_fwd_aligned, p_seq, consensus)
            if score is not None:
                candidates.append((score, "+", p_fwd_aligned))
        if p_rc_aligned:
            score = _score_probe_orientation(p_rc_aligned, reverse_complement(p_seq), consensus)
            if score is not None:
                candidates.append((score, "-", p_rc_aligned))
        if candidates:
            candidates.sort(key=lambda c: c[0])
            _, best_strand, best_aligned = candidates[0]
            ref_rows["probe"] = {"strand": best_strand, "aligned": best_aligned, "col": span_of(best_aligned)}
        else:
            ref_rows["probe"] = None
    else:
        ref_rows["probe"] = None

    cols = {label: (r["col"] if r else None) for label, r in ref_rows.items()}

    # Per-genome dot rows compared plainly against the population consensus, so at any
    # column an oligo's own IUPAC ambiguity code covers (e.g. Y = C or T), a genome whose
    # real base satisfies the primer/probe but merely differs from the plain
    # population-majority consensus (e.g. genome=T where consensus=C) rendered as a
    # spurious "mismatch" letter instead of ".", even though that genome is a perfectly
    # fine match for the actual oligo. Build a combined per-column reference that's the
    # oligo's own (possibly ambiguous) aligned character within its span, falling back to
    # plain consensus elsewhere, and compare against THAT (IUPAC-aware) instead.
    ref_chars = list(consensus)
    for label in ("forward", "probe", "reverse"):
        rr = ref_rows.get(label)
        if not rr or not rr["col"]:
            continue
        a, b = rr["col"]
        for i in range(a, min(b, n_cols)):
            ch = rr["aligned"][i]
            if ch != "-":
                ref_chars[i] = ch
    ref_seq = "".join(ref_chars)

    def _dot_matches(genome_char, ref_char):
        if genome_char == "-" or ref_char == "-":
            return genome_char == ref_char
        return genome_char.upper() in IUPAC.get(ref_char.upper(), ref_char.upper())

    def dotify(seq):
        return "".join("." if _dot_matches(seq[i], ref_seq[i]) else seq[i] for i in range(n_cols))

    meta = {g["genome_id"]: g for g in rows}
    out_rows = []
    for gid, aseq in aligned.items():
        g = meta.get(gid, {})
        out_rows.append({
            "genome_id": gid, "lineage": g.get("lineage", "NA"), "date": g.get("date", ""),
            "country": g.get("country", ""),
            "f_mismatch": g.get("f_mismatch"), "p_mismatch": g.get("p_mismatch"), "r_mismatch": g.get("r_mismatch"),
            "dots": dotify(aseq),
            "seq": aseq,        # exact aligned bases: lets the UI re-dot against any baseline
        })
    out_rows.sort(key=lambda r: (r["date"] or "", r["genome_id"]))

    return {
        "columns": n_cols,
        "f_col": cols.get("forward"), "p_col": cols.get("probe"), "r_col": cols.get("reverse"),
        "consensus": consensus, "ref": ref_seq, "rows": out_rows, "excluded": excluded,
        "ref_rows": {k: v for k, v in ref_rows.items() if v is not None},
    }
