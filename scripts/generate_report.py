#!/usr/bin/env python3
"""Build a single comprehensive PDF report covering every primer/probe set -- one section
per set (coverage, standard-material match, reference-genome binding, lineage
distribution, an automated "overall assessment" narrative, and the user's own memo/notes
including any attached photos, e.g. real RT-PCR gel images).

Pulls data from the already-running local webapp (so it reuses exactly the same
already-computed, already-fixed numbers the dashboard itself shows -- no separate
computation path to drift out of sync with) and renders a print-optimized HTML file,
then drives headless Chrome's DevTools Protocol directly (Page.printToPDF) to turn that
into an actual PDF -- not the `--print-to-pdf` CLI flag, whose header/footer-suppression
flag (--print-to-pdf-no-header) turned out not to work on this Chrome build; the CDP
`displayHeaderFooter` param does.
Confidential: written under output/ (gitignored), never committed.

Usage:
    python3 scripts/generate_report.py [--base-url http://localhost:8765] [--out output/report/report.pdf]
"""
import argparse
import base64
import html
import json
import os
import shutil
import socket
import struct
import subprocess
import tempfile
import sys
import time
import urllib.request
from pathlib import Path

def _detect_root() -> Path:
    """Same detection as webapp/server.py: prefer cwd (the project directory holding
    config/output) so the default --out path is correct whether this script is run from
    a source checkout or as an installed console script from a scaffolded project dir."""
    cwd = Path.cwd()
    if (cwd / "config" / "config.yaml").exists():
        return cwd
    return Path(__file__).resolve().parent.parent


ROOT = _detect_root()


def fetch_json(base_url, path):
    with urllib.request.urlopen(f"{base_url}{path}", timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch_bytes(base_url, path):
    with urllib.request.urlopen(f"{base_url}{path}", timeout=60) as resp:
        return resp.read(), resp.headers.get("Content-Type", "application/octet-stream")


def esc(s):
    return html.escape(str(s if s is not None else ""))

def scope_note(row):
    """' · 숙주 기준: 사람' -- which host scope the summary row counts (the pipeline's primary scope; see scripts/host_scope.py)."""
    label = {"human": "사람", "animal": "동물", "unknown": "미상", "all": "전체"}
    sc = (row or {}).get("host_scope")
    return f" · 숙주 기준: {label.get(sc, sc)}" if sc else ""


def mask_seq(seq):
    """Hide the literal primer/probe sequence in the report (this PDF tends to travel
    further/more loosely than the live dashboard, which stays behind localhost) -- shown as
    one bullet per base so the reader can still see its length, never the actual bases."""
    seq = seq or ""
    return "&bull;" * len(seq) if seq else "&mdash;"


def render_citation(citation):
    if not citation or not citation.get("text"):
        return ""
    link = f" <a href='{esc(citation['url'])}'>{esc(citation['url'])}</a>" if citation.get("url") else ""
    return f"""<div class="citation-box">
      <b>출처 (Reference)</b><br>{esc(citation['text'])}{link}
    </div>"""


def fmt_pct(v):
    return f"{v:.1f}%" if isinstance(v, (int, float)) else "—"


def is_validated_genome(g):
    return (
        g.get("f_mismatch") is not None and g["f_mismatch"] <= 1
        and g.get("p_mismatch") == 0
        and g.get("r_mismatch") is not None and g["r_mismatch"] <= 1
        and g.get("f_3p_ok") is not False
        and g.get("r_3p_ok") is not False
    )


def build_assessment(set_data):
    """Python port of the dashboard's own narrateSetAssessment() (app.js) -- same
    heuristics, same thresholds, so the PDF's verdict matches what the live page shows."""
    s = set_data.get("summary") or {}
    std = set_data.get("standard") or {}
    gene = set_data.get("gene") or {}
    binding = set_data.get("binding") or []
    points = []
    score = 0

    vpct = s.get("validated_pct") or 0
    if vpct >= 95:
        points.append((True, f"Population 커버리지 {fmt_pct(vpct)} — 현재 알려진 유전적 다양성 대부분을 커버합니다."))
        score += 2
    elif vpct >= 80:
        points.append((True, f"Population 커버리지 {fmt_pct(vpct)} — 양호한 편입니다."))
        score += 1
    elif vpct >= 50:
        points.append((None, f"Population 커버리지 {fmt_pct(vpct)} — 절반 안팎에 그쳐, 상당수 lineage/변이체를 놓칠 수 있습니다."))
        score -= 1
    else:
        points.append((False, f"Population 커버리지 {fmt_pct(vpct)} — 낮은 편이라 실제 시료의 상당수를 검출하지 못할 위험이 있습니다."))
        score -= 2

    dropout = s.get("probe_dropout") or 0
    primer_matched = s.get("primer_matched") or 0
    dropout_frac = (dropout / primer_matched) if primer_matched else 0
    if dropout_frac > 0.3:
        points.append((False, f"Primer는 매치되지만 Probe가 실패하는 게놈이 {dropout}개(약 {round(dropout_frac*100)}%)로 많습니다 — Probe 서열 재검토를 권장합니다."))
        score -= 1
    elif dropout_frac > 0.1:
        points.append((None, f"Probe dropout {dropout}개(약 {round(dropout_frac*100)}%) — 일부 변이체에서 Probe 결합이 실패합니다."))
    elif primer_matched:
        points.append((True, f"Probe dropout {dropout}개로 적어, Primer와 Probe 설계가 서로 잘 맞습니다."))
        score += 1

    if std.get("available") is False:
        points.append((None, "표준물질이 등록되어 있지 않아 실험실 양성대조군과의 비교는 평가에서 제외했습니다(표준물질은 선택 사항)."))
    elif std.get("ok"):
        points.append((True, f"표준물질(KDCA) 매치 O ({', '.join(std.get('matched') or [])}) — 실제 wet-lab 양성대조군에서도 동작할 것으로 기대됩니다."))
        score += 1
    else:
        points.append((False, "표준물질(KDCA) 매치 X — 실험실 양성대조군 검증에 실패할 가능성이 있습니다."))
        score -= 2

    not_found = [b for b in binding if not b.get("found")]
    if not_found:
        roles = "/".join(b.get("role", "?") for b in not_found)
        points.append((False, f"참조유전체에서 {roles} 결합위치를 찾지 못했습니다 — 서열 설계 자체를 재검토할 필요가 있습니다."))
        score -= 2
    elif gene.get("amplicon_size"):
        try:
            size = int(gene["amplicon_size"])
            if size > 300:
                points.append((None, f"Amplicon 크기 {size}bp — TaqMan qPCR 권장 범위(약 60~150bp)보다 큰 편입니다."))
            else:
                points.append((True, f"Amplicon 크기 {size}bp — qPCR에 적합한 범위입니다."))
                score += 1
        except (TypeError, ValueError):
            pass

    dh = set_data.get("dimer_hairpin")
    if dh:
        threshold = dh.get("dg_flag_threshold", -9)
        flagged_dimers = [(pair, v) for pair, v in (dh.get("dimer") or {}).items() if v.get("worst_dg", 0) <= threshold]
        flagged_hairpins = [(role, v) for role, v in (dh.get("hairpin") or {}).items() if v.get("worst_dg", 0) <= threshold]
        if flagged_dimers or flagged_hairpins:
            parts = [f"{pair}(&#916;G {v['worst_dg']:.1f})" for pair, v in flagged_dimers]
            parts += [f"{role} hairpin(&#916;G {v['worst_dg']:.1f})" for role, v in flagged_hairpins]
            points.append((False, f"이차구조 위험: {', '.join(parts)} — 자가/교차 결합이나 헤어핀 형성으로 증폭 효율이 떨어질 수 있습니다."))
            score -= 1
        else:
            points.append((True, "프라이머/프로브 간 뚜렷한 self-dimer·cross-dimer·hairpin 위험은 발견되지 않았습니다."))
            score += 1

    # Only a group's set_data carries this -- an individual set has no "other member sets"
    # to cross-check against, so this key is simply absent there and the check is skipped,
    # same pattern as gene/binding above.
    cross_binding = set_data.get("cross_binding")
    if cross_binding is not None:
        cross_genomes = {g.get("genome_id") for g in cross_binding
                          if not g.get("is_same_set") and g.get("genome_id")}
        if cross_genomes:
            points.append((False, f"그룹 구성 세트 간 의도치 않은 교차 결합이 population {len(cross_genomes)}개 게놈에서 발견됐습니다 "
                                   "— 비특이적 증폭(off-target) 가능성이 있으니 조합을 재검토하세요."))
            score -= 2
        else:
            points.append((True, "그룹 구성 세트 간 의도치 않은 교차 결합은 population 전체에서 발견되지 않았습니다."))
            score += 1

    if score >= 3:
        verdict, cls = "✅ 전반적으로 설계가 양호합니다.", "ok"
    elif score >= 0:
        verdict, cls = "ℹ️ 대체로 사용 가능하나 일부 보완이 필요합니다.", "warn"
    else:
        verdict, cls = "⚠️ 여러 문제가 발견되어 재설계 또는 보완을 권장합니다.", "bad"
    return verdict, cls, points


def lineage_distribution(genomes, lineage_totals=None):
    """[(lineage, validated_count, population_total, pct), ...] -- includes every lineage
    that appears in lineage_totals (the segment's full population breakdown) even if this
    set validated zero genomes there, so "0 out of 40" is visible rather than the lineage
    silently not appearing at all."""
    lineage_totals = lineage_totals or {}
    counts = {}
    for g in genomes:
        if is_validated_genome(g):
            counts[g.get("lineage", "NA")] = counts.get(g.get("lineage", "NA"), 0) + 1
    lineages = set(counts) | set(lineage_totals)
    rows = []
    for lin in lineages:
        c = counts.get(lin, 0)
        tot = lineage_totals.get(lin)
        pct = (100 * c / tot) if tot else None
        rows.append((lin, c, tot, pct))
    rows.sort(key=lambda r: -r[1])
    return rows


LOW_N_THRESHOLD = 20  # same rule-of-thumb cutoff used by the dashboard (webapp/static/app.js)


def render_lineage_table(lineage_rows):
    if not lineage_rows:
        return "<p class='muted'>없음</p>"
    body = "".join(
        f"<tr><td>{esc(lin)}{' <span class=\"warn\" style=\"font-size:8.5px;\">(n&lt;20)</span>' if tot is not None and tot < LOW_N_THRESHOLD else ''}</td>"
        f"<td>{c} / {esc(tot) if tot is not None else '?'}</td>"
        f"<td>{f'{pct:.0f}%' if pct is not None else '—'}</td></tr>"
        for lin, c, tot, pct in lineage_rows
    )
    return f"""<table class="mini-table">
      <thead><tr><th>Lineage</th><th>검출 / 전체 게놈 수</th><th>검출율</th></tr></thead>
      <tbody>{body}</tbody></table>
      <div class="muted" style="font-size:8.5px;margin-top:3px;">(n&lt;20) = population 내 이 리니지 표본 수가 적어 검출율 수치의 신뢰도가 낮음</div>"""


# ---- amplicon binding visualization (Python port of app.js's highlightAmplicon()) ----
# The dashboard's per-genome/per-standard-material detail view isn't just a mismatch/Tm/DG
# table -- it draws the actual F/Probe/R binding aligned character-for-character against the
# amplicon sequence (colored spans, 5'/3' direction arrows, a ":"/X match row). The PDF report
# originally only carried the table form; this reproduces that same alignment view so the
# report shows exactly what the dashboard shows, not a text summary of it.

IUPAC_RC = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
IUPAC_SET = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _rc(seq):
    return "".join(IUPAC_RC.get(c, "N") for c in reversed(seq.upper()))


def _iupac_matches(oligo_char, other_char):
    if oligo_char == "-" or other_char == "-":
        return True
    allowed = IUPAC_SET.get(oligo_char.upper(), oligo_char.upper())
    return other_char.upper() in allowed


def _count_iupac_mismatches(query, target):
    return sum(1 for a, b in zip(query, target) if not _iupac_matches(a, b))


def _mm_class(v):
    if v is None:
        return "mm-na"
    n = int(v)
    if n == 0:
        return "mm-0"
    if n == 1:
        return "mm-1"
    return "mm-2plus"


def _mm_text(v):
    return "—" if v is None else str(v)


def _strand_label(s):
    return {"+": "정방향", "-": "역방향"}.get(s, "")


def _thermo_spans(tm, dg):
    tm_html = f"<span class='amp-thermo mono'>Tm {round(float(tm), 1)}&deg;C</span>" if tm is not None else ""
    dg_html = f"<span class='amp-thermo mono'>&Delta;G {round(float(dg), 2)} kcal/mol</span>" if dg is not None else ""
    return tm_html + dg_html


# A standard-material fragment (~800bp) is far wider than the ~20nt F/P/R spans sitting
# inside it, so the literal alignment is mostly "-" filler between spans -- printing that
# filler literally either overflows the page width or (if wrapped by the browser) wraps
# each of the 4 stacked rows at a different character position since they don't wrap
# identically (long runs of "-"/" " vs. <mark>-wrapped sequence letters don't break at the
# same point), which is what made the block "안 보이는" -- rows drifting out of column
# alignment. Fixed two ways: (1) collapse any gap longer than GAP_KEEP_THRESHOLD into a
# short "...(Nbp 생략)..." placeholder (that's the actual, useful information -- "there's
# nothing here" -- not hundreds of literal dashes), which alone gets nearly everything onto
# one line; (2) as a safety net for whatever's left, chunk into fixed-width groups computed
# in Python (not CSS wrapping) so all 4 rows in a group always break at the exact same index.
GAP_KEEP_THRESHOLD = 8
AMP_ALIGN_CHUNK_WIDTH = 120


def _compact_amp_gaps(role, dir_chars, primer_chars, colon_chars, seq_chars):
    n = len(role)
    out_role, out_dir, out_primer, out_colon, out_seq = [], [], [], [], []
    i = 0
    while i < n:
        j = i
        while j < n and role[j] == role[i]:
            j += 1
        run_len = j - i
        if role[i] is None and run_len > GAP_KEEP_THRESHOLD:
            # ASCII-only: a Korean-glyph placeholder here (used to say "...(Nbp 생략)...")
            # rendered wider per-character than the plain spaces the other 3 rows get for
            # this same span (the monospace font has no Hangul glyphs, so the browser
            # substitutes a CJK fallback font whose characters are wider), which visibly
            # drifted every row out of column alignment after the first such gap -- see the
            # legend line above for what "(Nbp)" means instead of spelling it out here.
            label = f"...({run_len}bp)..."
            out_role.extend([None] * len(label))
            out_dir.extend([" "] * len(label))
            out_primer.extend([" "] * len(label))
            out_colon.extend([" "] * len(label))
            out_seq.extend(list(label))
        else:
            out_role.extend(role[i:j])
            out_dir.extend(dir_chars[i:j])
            out_primer.extend(primer_chars[i:j])
            out_colon.extend(colon_chars[i:j])
            out_seq.extend(seq_chars[i:j])
        i = j
    return out_role, out_dir, out_primer, out_colon, out_seq


def _wrap_role_range(chars, role_slice):
    out, i, m = [], 0, len(chars)
    while i < m:
        j = i
        while j < m and role_slice[j] == role_slice[i]:
            j += 1
        chunk = esc("".join(chars[i:j]))
        out.append(f'<mark class="{role_slice[i]}">{chunk}</mark>' if role_slice[i] else chunk)
        i = j
    return "".join(out)


def _render_amp_align(role, dir_chars, primer_chars, colon_chars, seq_chars):
    role, dir_chars, primer_chars, colon_chars, seq_chars = _compact_amp_gaps(
        role, dir_chars, primer_chars, colon_chars, seq_chars)
    n = len(role)
    groups = []
    for start in range(0, n, AMP_ALIGN_CHUNK_WIDTH):
        end = min(n, start + AMP_ALIGN_CHUNK_WIDTH)
        role_slice = role[start:end]
        groups.append(f"""
        <div class="amp-align-group">
          <div class="amp-align-row amp-align-dir">{esc("".join(dir_chars[start:end]))}</div>
          <div class="amp-align-row">{_wrap_role_range(primer_chars[start:end], role_slice)}</div>
          <div class="amp-align-row amp-align-colon">{esc("".join(colon_chars[start:end]))}</div>
          <div class="amp-align-row">{_wrap_role_range(seq_chars[start:end], role_slice)}</div>
        </div>""")
    return "".join(groups)


def highlight_amplicon(g, set_data, labels=None):
    labels = labels or {"f": "FORWARD", "p": "PROBE", "r": "REVERSE"}
    verdict_ok = is_validated_genome(g)
    verdict_text = "✅ 검증 통과 — F/R mismatch ≤1, Probe 완전일치" if verdict_ok else "⚠️ 검증 실패"
    seq = g.get("amplicon_seq") or ""

    if not seq:
        def no_amp_row(role_class, label, mm, tm, dg):
            return f"""<div class="amp-role-row">
              <span class="amp-role-tag {role_class}">{esc(label)}</span>
              <span class="msa-mm-tag {_mm_class(mm)}">mismatch {esc(_mm_text(mm))}</span>
              {_thermo_spans(tm, dg)}
            </div>"""
        return f"""
          <div class="amp-verdict {'ok' if verdict_ok else 'bad'}">{verdict_text}</div>
          <div class="amp-role-summary">
            {no_amp_row("f", labels["f"], g.get("f_mismatch"), g.get("f_tm"), g.get("f_dg"))}
            {no_amp_row("p", labels["p"], g.get("p_mismatch"), g.get("p_tm"), g.get("p_dg"))}
            {no_amp_row("r", labels["r"], g.get("r_mismatch"), g.get("r_tm"), g.get("r_dg"))}
          </div>
          <p class="muted">F/R이 짝을 이루지 못해 amplicon 서열을 얻지 못했습니다.</p>"""

    n = len(seq)
    raw_seq = seq
    amp_start = g.get("amplicon_start")
    amp_end = g.get("amplicon_end")
    p_start_abs = g.get("probe_start")
    p_seq = set_data.get("P_seq") or ""
    f_seq = set_data.get("F_seq") or ""
    r_seq = set_data.get("R_seq") or ""
    f_len, r_len = len(f_seq), len(r_seq)

    # mfeprimer's amplicon-extraction convention assumes the genome record is deposited on
    # the same strand as the reference used to design F/R, so F always sits at the
    # amplicon's own start and R (reverse-complemented) at its end -- true for most
    # genomes, but not guaranteed for one deposited under the opposite strand convention
    # (confirmed for OQ919513/VendorD_GPC_PREVIOUS1_LASV1_P2). Self-detect which orientation
    # actually fits (same comparison the dashboard's highlightAmplicon() makes) so both
    # agree, rather than assuming the convention always holds.
    flipped = False
    if f_len and r_len and n:
        r_rc = _rc(r_seq)
        direct_mm = _count_iupac_mismatches(f_seq, raw_seq[:f_len]) + _count_iupac_mismatches(r_rc, raw_seq[max(0, n - r_len):])
        flipped_check = _rc(raw_seq)
        flipped_mm = _count_iupac_mismatches(f_seq, flipped_check[:f_len]) + _count_iupac_mismatches(r_rc, flipped_check[max(0, n - r_len):])
        flipped = flipped_mm < direct_mm
    seq = _rc(raw_seq) if flipped else raw_seq

    # F/R coordinate LABELS always describe true genome positions -- when the display got
    # flipped, F is actually the oligo sitting at the genome's amp_end side, and vice versa.
    f_coord_start = (amp_end - f_len + 1) if (flipped and amp_end is not None) else amp_start
    f_coord_end = amp_end if flipped else ((amp_start + f_len - 1) if amp_start is not None else None)
    r_coord_start = amp_start if flipped else ((amp_end - r_len + 1) if amp_end is not None else None)
    r_coord_end = (amp_start + r_len - 1) if (flipped and amp_start is not None) else amp_end

    # The probe's own absolute genome coordinate (used for its label) is unaffected by the
    # display flip -- only its LOCAL character offset within the now-possibly-reoriented
    # `seq` string needs to account for it, mirrored around the amplicon span (reversing an
    # n-length string maps old index i -> n-1-i, so a [a,b) span becomes [n-b, n-a)).
    p_local_start = p_local_end = None
    if p_start_abs is not None and amp_start is not None and amp_end is not None:
        if flipped:
            p_local_end = amp_end - p_start_abs + 1
            p_local_start = p_local_end - len(p_seq)
        else:
            p_local_start = p_start_abs - amp_start
            p_local_end = p_local_start + len(p_seq)

    r_start = max(0, n - r_len)
    role = [None] * n
    for k in range(min(f_len, n)):
        role[k] = "f"
    for k in range(r_start, n):
        role[k] = "r"

    p_strand = "+"
    has_probe = p_local_start is not None and p_local_end is not None and 0 <= p_local_start < p_local_end <= n
    if has_probe:
        p_hit = seq[p_local_start:p_local_end]
        if _count_iupac_mismatches(p_seq, _rc(p_hit)) < _count_iupac_mismatches(p_seq, p_hit):
            p_strand = "-"
        for k in range(p_local_start, p_local_end):
            role[k] = "p"

    r_display = _rc(r_seq)
    p_display = _rc(p_seq) if p_strand == "-" else p_seq
    primer_chars = ["-"] * n
    colon_chars = [" "] * n
    dir_chars = [" "] * n
    for k in range(n):
        if role[k] == "f" and k < len(f_seq):
            primer_chars[k] = f_seq[k]
        elif role[k] == "r" and (k - r_start) < len(r_display):
            primer_chars[k] = r_display[k - r_start]
        elif role[k] == "p" and (k - p_local_start) < len(p_display):
            primer_chars[k] = p_display[k - p_local_start]
        if role[k]:
            colon_chars[k] = ":" if _iupac_matches(primer_chars[k], seq[k]) else "X"

    def write_direction(start, end, strand):
        if end - start < 4:
            return
        arrow = "<" if strand == "-" else ">"
        for k in range(start, end):
            dir_chars[k] = arrow
        left = "3'" if strand == "-" else "5'"
        right = "5'" if strand == "-" else "3'"
        dir_chars[start], dir_chars[start + 1] = left[0], left[1]
        dir_chars[end - 2], dir_chars[end - 1] = right[0], right[1]

    write_direction(0, min(f_len, n), "+")
    write_direction(r_start, n, "-")
    if has_probe:
        write_direction(p_local_start, p_local_end, p_strand)

    def role_row(role_class, label, strand, mm, start, end, tm, dg):
        strand_html = f" ({_strand_label(strand)})" if strand else ""
        pos_html = f"<span class='mono'>{start}-{end}</span>" if start is not None else "<span class='muted'>매치 안됨</span>"
        return (f"""<span class="amp-role-tag {role_class}">{esc(label)}{strand_html}</span>
          {pos_html}
          <span class="msa-mm-tag {_mm_class(mm)}">mismatch {esc(_mm_text(mm))}</span>
          {_thermo_spans(tm, dg)}""")

    p_end_pos = (p_start_abs + len(p_seq) - 1) if (has_probe and p_start_abs is not None) else None

    align_html = _render_amp_align(role, dir_chars, primer_chars, colon_chars, list(seq))

    return f"""
      <div class="amp-verdict {'ok' if verdict_ok else 'bad'}">{verdict_text}</div>
      <div class="amp-role-summary">
        <div class="amp-role-row">{role_row("f", labels["f"], "-" if flipped else "+", g.get("f_mismatch"), f_coord_start, f_coord_end, g.get("f_tm"), g.get("f_dg"))}</div>
        <div class="amp-role-row">{role_row("p", labels["p"], p_strand if has_probe else None, g.get("p_mismatch"), p_start_abs if has_probe else None, p_end_pos, g.get("p_tm"), g.get("p_dg"))}</div>
        <div class="amp-role-row">{role_row("r", labels["r"], "+" if flipped else "-", g.get("r_mismatch"), r_coord_start, r_coord_end, g.get("r_tm"), g.get("r_dg"))}</div>
      </div>
      <div class="amp-align">{align_html}</div>"""


def render_standard_detail(standard_genomes, set_data):
    """Only the standard materials that actually bind (validated) -- a print report is
    read as "what works", not a diagnostic surface, so materials this set/group can't
    bind at all (e.g. GC2, NP1/NP2, L1/L2 when only GC1 matches) are omitted entirely
    rather than listed as failures (the interactive dashboard already covers that).
    Renders the same F/Probe/R binding alignment visualization as the dashboard's
    renderStandardDetail()/highlightAmplicon(), not just a mismatch/Tm/DG table."""
    matched = [g for g in standard_genomes if is_validated_genome(g)]
    if not matched:
        return "<p class='muted'>매치되는 표준물질 없음</p>"
    legend = (
        "<p class='amp-legend'>"
        "<span class='amp-role-tag f'>F</span>정방향/역방향은 이 표준물질 서열 기준 결합 방향&nbsp;&middot;&nbsp;"
        "<b>:</b>=염기 매치, <b>X</b>=불일치, <b>-</b>=다른 프라이머 결합 구간&nbsp;&middot;&nbsp;"
        "<b>...(Nbp)...</b>=프라이머 사이 먼 구간(N염기, 정보 없음)을 줄여 표시"
        "</p>"
    )
    blocks = [
        f"""<div class="standard-material-block">
          <div class="standard-material-head">
            <span class="mono">{esc(g.get('genome_id'))}</span>
          </div>
          {highlight_amplicon(g, set_data)}
        </div>"""
        for g in matched
    ]
    return legend + "".join(blocks)


def render_combo_detail(combo_rows, sets_by_id, label_suffix=""):
    """Same idea as render_standard_detail but for group cross-binding combo rows (F/R/Probe
    possibly from three DIFFERENT member sets) -- validated combos only, one representative
    genome per distinct combo (a report page isn't the place for hundreds of per-genome rows
    -- the interactive dashboard has those). Renders the same F/Probe/R binding alignment
    visualization as render_standard_detail()/the dashboard's renderComboGenomeTable(),
    labeling each span with its OWN origin set since a cross-set combo's F and R can come
    from different designs."""
    matched = [g for g in combo_rows if is_validated_genome(g)]
    if not matched:
        return f"<p class='muted'>매치되는 조합 없음{label_suffix}</p>"
    by_combo = {}
    for g in matched:
        by_combo.setdefault(g.get("combo_id", ""), g)

    def seq_for_role(set_id, role):
        s = sets_by_id.get(set_id) or {}
        return s.get("R_seq" if role == "R" else "F_seq") or ""

    blocks = []
    for combo_id, g in sorted(by_combo.items()):
        same = g.get("is_same_set")
        tag = "" if same else " <span class='pill pill-bad'>&#9888; 교차 결합(다른 세트)</span>"
        # f_oligo_id/r_oligo_id are "<set_id>::F" or "<set_id>::R" -- the oligo that ends up
        # playing the geometric left/right role isn't always that set's own F/R (see
        # group_cross_binding.py's ID_SEP comment), so which sequence to display must be
        # read from the resolved role, not assumed from f_set_id/r_set_id alone.
        f_role = (g.get("f_oligo_id") or "").split("::")[-1] or "F"
        r_role = (g.get("r_oligo_id") or "").split("::")[-1] or "R"
        combo_set_data = {
            "F_seq": seq_for_role(g.get("f_set_id"), f_role),
            "R_seq": seq_for_role(g.get("r_set_id"), r_role),
            "P_seq": (sets_by_id.get(g.get("p_set_id")) or {}).get("P_seq", "") if g.get("p_set_id") else "",
        }
        labels = {
            "f": f"FORWARD ({g.get('f_set_name') or g.get('f_set_id')})",
            "p": f"PROBE ({g.get('p_set_name') or g.get('p_set_id')})" if g.get("p_set_id") else "PROBE",
            "r": f"REVERSE ({g.get('r_set_name') or g.get('r_set_id')})",
        }
        blocks.append(f"""<div class="standard-material-block">
          <div class="standard-material-head">
            <span class="mono">{esc(g.get("combo_label") or combo_id)}</span>{tag}
          </div>
          {highlight_amplicon(g, combo_set_data, labels)}
        </div>""")
    return "".join(blocks)


def render_note_images(base_url, set_id, note_images):
    if not note_images:
        return ""
    figs = []
    for img in note_images:
        try:
            raw, content_type = fetch_bytes(base_url, f"/api/sets/{set_id}/note/images/{img['id']}")
        except Exception:  # noqa: BLE001
            continue
        b64 = base64.b64encode(raw).decode("ascii")
        caption = esc(img.get("caption", ""))
        figs.append(f'<figure><img src="data:{content_type};base64,{b64}"><figcaption>{caption}</figcaption></figure>')
    if not figs:
        return ""
    return f'<div class="note-images">{"".join(figs)}</div>'


def render_verification_images(base_url, set_id, vid, images):
    if not images:
        return ""
    figs = []
    for img in images:
        try:
            raw, content_type = fetch_bytes(base_url, f"/api/sets/{set_id}/verify/{vid}/images/{img['id']}")
        except Exception:  # noqa: BLE001
            continue
        b64 = base64.b64encode(raw).decode("ascii")
        caption = esc(img.get("caption", ""))
        figs.append(f'<figure><img src="data:{content_type};base64,{b64}"><figcaption>{caption}</figcaption></figure>')
    if not figs:
        return ""
    return f'<div class="note-images">{"".join(figs)}</div>'


def render_verifications(base_url, set_id, verifications, standard_genomes):
    """Wet-lab (real RT-PCR) verification records -- standard material used, whether the
    experiment actually detected it, and how that compares to the in-silico prediction for
    that same material (is_validated_genome(), now including the 3'-end check), plus any
    attached result images (amplification curve, gel photo, etc.)."""
    if not verifications:
        return ""
    std_by_id = {g.get("genome_id"): g for g in (standard_genomes or [])}
    blocks = []
    for v in sorted(verifications, key=lambda x: x.get("created_at", ""), reverse=True):
        material = v.get("standard_material", "")
        g = std_by_id.get(material)
        in_silico_ok = is_validated_genome(g) if g else None
        wet_ok = v.get("detected")
        wet_label = "검출(양성)" if wet_ok is True else "미검출(음성)" if wet_ok is False else "판정 보류"
        wet_cls = "pill-ok" if wet_ok is True else "pill-bad" if wet_ok is False else ""
        in_silico_label = "O 매치" if in_silico_ok else ("X 불일치" if in_silico_ok is False else "데이터없음")
        in_silico_cls = "pill-ok" if in_silico_ok else "pill-bad" if in_silico_ok is False else ""
        agree_html = ""
        if in_silico_ok is not None and wet_ok is not None:
            agree_html = ("<span class='pill pill-ok'>&#10003; 인실리코-실험 일치</span>" if in_silico_ok == wet_ok
                           else "<span class='pill pill-bad'>&#9888; 인실리코-실험 불일치</span>")
        meta_bits = [esc(v.get("created_at", ""))]
        if v.get("ct_value"):
            meta_bits.append(f"Ct {esc(v['ct_value'])}")
        images_html = render_verification_images(base_url, set_id, v.get("id", ""), v.get("images") or [])
        blocks.append(f"""<div class="standard-material-block">
          <div class="standard-material-head">
            <span class="mono">{esc(material)}</span>
            {f"<span class='pill'>{esc(v.get('fluorophore'))}</span>" if v.get('fluorophore') else ""}
            <span class="pill {wet_cls}">실험: {wet_label}</span>
            <span class="pill {in_silico_cls}">인실리코: {in_silico_label}</span>
            {agree_html}
            <span class="muted" style="margin-left:auto;">{' &middot; '.join(meta_bits)}</span>
          </div>
          {f"<p style='margin:4px 0 0;'>{esc(v.get('note'))}</p>" if v.get('note') else ""}
          {images_html}
        </div>""")
    return "".join(blocks)


def _pi_for_range(diversity_data, segment, start, end):
    """Population-wide (all-lineage) average pixy pi across a reference-coordinate span --
    the report's print-friendly stand-in for the dashboard's interactive per-lineage overlay
    chart: one number per binding site instead of a chart, but same underlying data."""
    if not diversity_data or start in (None, "") or end in (None, ""):
        return None
    seg_data = (diversity_data.get("segments") or {}).get((segment or "").lower())
    if not seg_data:
        return None
    windows = seg_data.get("windows") or {}
    window_size = seg_data.get("window_size") or 1
    lineages = list(windows.keys())
    if not lineages:
        return None
    try:
        lo, hi = sorted((int(start), int(end)))
    except (TypeError, ValueError):
        return None
    idx_lo = max(0, (lo - 1) // window_size)
    idx_hi = (hi - 1) // window_size
    vals = []
    for lin in lineages:
        wins = windows.get(lin) or []
        for idx in range(idx_lo, min(idx_hi, len(wins) - 1) + 1):
            pi = wins[idx].get("pi")
            if pi is not None:
                vals.append(pi)
    return sum(vals) / len(vals) if vals else None


def _pi_class(pi):
    if pi is None:
        return ""
    if pi < 0.05:
        return "ok"
    if pi < 0.15:
        return "warn"
    return "bad"


def render_binding(binding, segment=None, diversity_data=None):
    if not binding:
        return "<p class='muted'>참조유전체 결합 정보 없음</p>"
    rows = []
    for b in binding:
        pi = _pi_for_range(diversity_data, segment, b.get("start"), b.get("end")) if diversity_data else None
        pi_cell = f"<td class='{_pi_class(pi)}'>&pi;&nbsp;{pi:.3f}</td>" if pi is not None else "<td class='muted'>&pi; &mdash;</td>"
        if not b.get("found"):
            rows.append(f"<tr><td>{esc(b.get('role'))}</td><td colspan='3' class='bad'>NO MATCH</td></tr>")
        else:
            rows.append(
                f"<tr><td>{esc(b.get('role'))}</td><td>{esc(b.get('chrom'))}:{esc(b.get('start'))}-{esc(b.get('end'))}"
                f" ({esc(b.get('strand'))})</td><td>mismatch {esc(b.get('mismatches'))}</td>{pi_cell}</tr>"
            )
    note = ""
    if diversity_data is not None:
        note = ("<div class='muted' style='font-size:8.5px;margin-top:4px;'>&pi; = population 전체 리니지 평균 "
                "nucleotide diversity(pixy, window=1bp) &mdash; 낮을수록(초록) 보존, 높을수록(빨강) 다양</div>")
    return f"<table class='mini-table'><tbody>{''.join(rows)}</tbody></table>{note}"


def render_set_section(base_url, set_data, index, total, diversity_data=None):
    sid = set_data["set_id"]
    display_name = set_data.get("display_name") or sid
    s = set_data.get("summary") or {}
    std = set_data.get("standard") or {}
    gene = set_data.get("gene") or {}
    genomes = set_data.get("genomes") or []
    verdict, verdict_cls, points = build_assessment(set_data)
    lineage_rows = lineage_distribution(genomes, set_data.get("lineage_totals") or {})
    note = set_data.get("note") or ""
    note_images_html = render_note_images(base_url, sid, set_data.get("note_images") or [])
    verifications_html = render_verifications(base_url, sid, set_data.get("verifications") or [],
                                               set_data.get("standard_genomes") or [])

    points_html = "".join(
        f"<li class='{'ok' if ok is True else 'bad' if ok is False else 'warn'}'>"
        f"{'✅' if ok is True else '⚠️' if ok is False else 'ℹ️'} {esc(text)}</li>"
        for ok, text in points
    )

    return f"""
    <section class="report-page">
      <div class="page-header">
        <span class="page-index">{index} / {total}</span>
        <h1>{esc(display_name)}</h1>
        <span class="pill pill-{esc(set_data.get('segment', '').lower())}">{esc(set_data.get('segment'))}</span>
        <span class="pill">{esc(set_data.get('source'))}</span>
      </div>
      <div class="detail-sub">
        {f"ID: <span class='mono'>{esc(sid)}</span> &middot; " if display_name != sid else ""}Target gene: <b>{esc(gene.get('target_gene', 'NA'))}</b>
        {f" &middot; amplicon {esc(gene.get('chrom',''))}:{esc(gene.get('amplicon_start'))}..{esc(gene.get('amplicon_end'))} ({esc(gene.get('amplicon_size'))}bp)" if gene.get('amplicon_start') else ""}
      </div>

      <div class="assessment-box {esc(verdict_cls)}">
        <div class="assessment-verdict">{esc(verdict)}</div>
        <ul class="assessment-points">{points_html}</ul>
      </div>

      <div class="grid2">
        <div class="box">
          <h3>프라이머 / 프로브 서열</h3>
          <div class="seq-row"><span class="tag">Forward</span><span class="mono seq-redacted">{mask_seq(set_data.get('F_seq'))} <span class="muted">({len(set_data.get('F_seq') or '')}nt)</span></span></div>
          <div class="seq-row"><span class="tag">Probe</span><span class="mono seq-redacted">{mask_seq(set_data.get('P_seq'))} <span class="muted">({len(set_data.get('P_seq') or '')}nt)</span></span></div>
          <div class="seq-row"><span class="tag">Reverse</span><span class="mono seq-redacted">{mask_seq(set_data.get('R_seq'))} <span class="muted">({len(set_data.get('R_seq') or '')}nt)</span></span></div>
          <div class="note" style="margin-top:4px;font-size:9.5px;color:#889;">서열은 대외비 처리를 위해 리포트에서 가렸습니다 — 대시보드에서 확인하세요.</div>
          {render_citation(set_data.get('citation'))}
        </div>
        <div class="box">
          <h3>Population 커버리지</h3>
          <table class="mini-table">
            <tr><td>Primer mismatch=0</td><td>{fmt_pct(s.get('primer_coverage_mm0_pct'))}</td></tr>
            <tr><td>Primer mismatch&le;1</td><td>{fmt_pct(s.get('primer_coverage_mm1_pct'))}</td></tr>
            <tr><td>Primer mismatch&le;2</td><td>{fmt_pct(s.get('primer_coverage_mm2_pct'))}</td></tr>
            <tr><td>Probe 완전일치</td><td>{fmt_pct(s.get('probe_coverage_mm0_pct'))}</td></tr>
            <tr><td>Probe dropout</td><td>{esc(s.get('probe_dropout'))}</td></tr>
            <tr class="highlight"><td><b>최종 Validated</b></td><td><b>{esc(s.get('validated'))} / {esc(s.get('total_genomes'))} ({fmt_pct(s.get('validated_pct'))})</b>{esc(scope_note(s))}</td></tr>
          </table>
        </div>
        <div class="box">
          <h3>표준물질(KDCA) 매치</h3>
          <div class="big-ox {'ok' if std.get('ok') else 'bad'}">{'O' if std.get('ok') else 'X'}</div>
          <div class="muted">{esc(', '.join(std.get('matched') or []) or '매치 없음')}</div>
        </div>
        <div class="box">
          <h3>참조유전체 결합위치</h3>
          {render_binding(set_data.get('binding'), set_data.get('segment'), diversity_data)}
        </div>
      </div>

      <div class="box">
        <h3>검증 통과 게놈의 리니지 분포</h3>
        {render_lineage_table(lineage_rows)}
      </div>

      <div class="box">
        <h3>표준물질(KDCA) 상세 결합 정보 (매치되는 것만)</h3>
        {render_standard_detail(set_data.get('standard_genomes') or [], set_data)}
      </div>

      {f'''<div class="box">
        <h3>실험(Wet-lab) 검증 기록</h3>
        <p class="muted" style="margin:0 0 8px;">표준물질(양성대조물질)로 실제 프라이머를 테스트한 결과와 인실리코 예측 비교입니다.</p>
        {verifications_html}
      </div>''' if verifications_html else ""}

      {f'<div class="box"><h3>메모 / 노트</h3><p class="note-text">{esc(note)}</p>{note_images_html}</div>' if note or note_images_html else ""}
    </section>
    """


def render_group_member_coverage(member_details):
    if not member_details:
        return "<p class='muted'>세트 정보 없음</p>"
    rows = []
    for s in member_details:
        sm = s.get("summary") or {}
        rows.append(
            f"<tr><td class='mono'>{esc(s.get('display_name') or s.get('set_id'))}</td>"
            f"<td>{fmt_pct(sm.get('primer_coverage_mm1_pct'))}</td>"
            f"<td>{fmt_pct(sm.get('probe_coverage_mm0_pct'))}</td>"
            f"<td><b>{fmt_pct(sm.get('validated_pct'))}</b></td></tr>"
        )
    return f"""<table class="mini-table">
      <thead><tr><th>프라이머</th><th>Primer mm&le;1</th><th>Probe 완전일치</th><th>Validated</th></tr></thead>
      <tbody>{''.join(rows)}</tbody></table>"""


def render_combo_summary_chips(combo_rows):
    """Distinct combo_ids found, with genome counts and same/cross-set flag -- a compact
    overview (not a per-genome table, which the interactive dashboard already provides in
    full) of every geometric F/R pairing pooled search actually found across the whole
    population, cross-set combos included."""
    if not combo_rows:
        return "<p class='muted'>발견된 결합 조합 없음</p>"
    by_combo = {}
    for g in combo_rows:
        key = g.get("combo_id", "")
        by_combo.setdefault(key, {"count": 0, "same": g.get("is_same_set"), "label": g.get("combo_label") or key})
        by_combo[key]["count"] += 1
    rows = []
    for combo_id, info in sorted(by_combo.items(), key=lambda kv: kv[0]):
        cls = "" if info["same"] else "bad"
        tag = "정상 (같은 세트)" if info["same"] else "⚠ 교차 결합 (다른 세트)"
        rows.append(f"<tr><td class='mono {cls}'>{esc(info["label"])}</td><td>{tag}</td><td>{info['count']}개 게놈</td></tr>")
    return f"<table class='mini-table small'><tbody>{''.join(rows)}</tbody></table>"


def render_group_section(group_data, index, total):
    gid = group_data["group_id"]
    members = group_data.get("members") or []
    lineage_counts = group_data.get("lineage_counts") or {}
    lineage_totals = group_data.get("lineage_totals") or {}
    lineages = set(lineage_counts) | set(lineage_totals)
    lineage_rows = []
    for lin in lineages:
        c = lineage_counts.get(lin, 0)
        tot = lineage_totals.get(lin)
        pct = (100 * c / tot) if tot else None
        lineage_rows.append((lin, c, tot, pct))
    lineage_rows.sort(key=lambda r: -r[1])

    summary = group_data.get("summary") or {}
    sets_by_id = {s["set_id"]: s for s in (group_data.get("member_details") or [])}

    member_details = group_data.get("member_details") or []
    std_ok = any((m.get("standard") or {}).get("ok") for m in member_details)
    matched = sorted({lin for m in member_details if (m.get("standard") or {}).get("ok")
                       for lin in ((m.get("standard") or {}).get("matched") or [])})
    group_set_data = {
        "summary": {"validated_pct": summary.get("validated_pct")},
        "standard": {"ok": std_ok, "matched": matched},
        "gene": {},
        "binding": [],
        "dimer_hairpin": group_data.get("dimer_hairpin"),
        "cross_binding": group_data.get("combo_genomes") or [],
    }
    verdict, verdict_cls, points = build_assessment(group_set_data)
    points_html = "".join(
        f"<li class='{'ok' if ok is True else 'bad' if ok is False else 'warn'}'>"
        f"{'✅' if ok is True else '⚠️' if ok is False else 'ℹ️'} {esc(text)}</li>"
        for ok, text in points
    )

    return f"""
    <section class="report-page">
      <div class="page-header">
        <span class="page-index">그룹 {index} / {total}</span>
        <h1>&#9633; {esc(group_data.get('label'))}</h1>
        <span class="pill pill-{esc(group_data.get('segment', '').lower())}">{esc(group_data.get('segment'))}</span>
      </div>
      <div class="detail-sub">구성 세트: {esc(', '.join((group_data.get('member_names') or {}).get(m, m) for m in members))}</div>

      <div class="assessment-box {verdict_cls}">
        <div class="assessment-verdict">합집합 커버리지: {esc(summary.get('validated'))} / {esc(summary.get('total_genomes'))}{esc(scope_note(summary))}
          ({fmt_pct(summary.get('validated_pct'))}) &mdash; {esc(verdict)}</div>
        <div class="muted">이 그룹의 세트 중 하나라도 통과(validated)하면 그 게놈은 합집합에 포함됩니다.</div>
        <ul class="assessment-points">{points_html}</ul>
      </div>

      <div class="grid2">
        <div class="box">
          <h3>세트별 Population 커버리지</h3>
          {render_group_member_coverage(group_data.get('member_details'))}
        </div>
        <div class="box">
          <h3>검증 통과 게놈의 리니지 분포</h3>
          {render_lineage_table(lineage_rows)}
        </div>
      </div>

      <div class="box">
        <h3>표준물질(KDCA) 매치 (교차 결합 포함, 매치되는 것만)</h3>
        {render_combo_detail(group_data.get('standard_combo') or [], sets_by_id)}
      </div>

      <div class="box">
        <h3>참조유전체 결합위치 (교차 결합 포함, 매치되는 것만)</h3>
        {render_combo_detail(group_data.get('refgenome_combo') or [], sets_by_id)}
      </div>

      <div class="box">
        <h3>Population 교차 결합 분석 요약</h3>
        <div class="muted" style="margin-bottom:6px;">그룹의 모든 F/R/Probe를 하나로 모아 population 전체에서 실제로 발견된
          결합 조합입니다(세트 간 의도치 않은 교차 결합 포함). 게놈별 전체 상세는 웹 대시보드에서 확인할 수 있습니다.</div>
        {render_combo_summary_chips(group_data.get('combo_genomes') or [])}
      </div>
    </section>
    """


WORKFLOW_STAGES = [
    ("1", "데이터 수집 (fetch)",
     "Nextstrain(nextstrain.org/lassa)에서 LASV S segment / L segment의 최신 서열과 메타데이터"
     "(lineage, 국가, 날짜 등)를 내려받습니다."),
    ("2", "필터링 (filter)",
     "모호 염기(N)가 하나라도 있거나 길이가 기준(S&ge;3000bp, L&ge;6000bp) 미만인 게놈은 제외해, "
     "mfeprimer 검색에 쓰이는 population DB를 만듭니다."),
    ("3~4", "F/R · 프로브 결합 탐색 (primer_match / probe_match / alignment_table)",
     "각 세트의 Forward/Reverse 프라이머 쌍은 두 엔진으로 모든 게놈에서 찾습니다: mfeprimer --virus(mismatch&le;2)와 "
     "게놈 전체 직접 탐색(모든 위치&middot;양쪽 가닥). Probe도 mfeprimer search와 직접 탐색을 함께 씁니다. 두 결과를 비교해 더 나은 "
     "쪽을 쓰고 어느 쪽인지 기록합니다. mismatch는 mfeprimer가 보고한 서열이 아니라 그 위치의 실제 게놈 염기와 프라이머를 "
     "비교해 세며(IUPAC 축퇴염기 인식), F/R 쌍은 서로 마주 보는 유효한 amplicon(최대 2000 bp)과 프라이머별 3&prime; 끝 일치를 "
     "확인하고, Probe는 amplicon 안에 있어야 합니다. mismatch 0/1/2 단계별 커버리지는 이 게놈별 표에서 계산합니다."),
    ("5", "판정 (reconcile / aggregate)",
     "validated = (primer mismatch &le; 설정값) AND (probe mismatch = 설정값) 조건을 만족하는 게놈만 "
     "최종 통과로 집계합니다. mfeprimer는 후보 결합 위치를 찾는 데 쓰고, mismatch 개수는 그 위치의 실제 게놈 염기와 "
     "프라이머를 직접 비교해 다시 셉니다(mfeprimer의 --misMatch는 3&prime; 끝 9염기만 반영해 5&prime; 쪽 불일치를 놓치기 때문). "
     "이 재계산은 게놈 전체를 무차별 대입으로 다시 센 독립 계산과 모든 세트에서 일치함을 확인했습니다. "
     "프로브는 두 프라이머 사이에 있어야 통과로 셉니다."),
    ("6", "부가 분석 (같은 탐색 로직 재사용)",
     "판정과 별개로 검증에 도움이 되는 추가 정보를 계산합니다 &mdash; 표준물질 매치(선택 사항: 등록하지 않으면 건너뜀), "
     "참조유전체 결합위치(대표 유전체 1개에서 F/Probe/R 결합 좌표), 게놈별 F/Probe/R 미스매치&middot;Tm&middot;&Delta;G, "
     "전체 게놈 정렬(mafft MSA), 그룹 교차 결합 분석(그룹의 모든 올리고를 서로 짝지어 두 엔진으로 탐색; 그룹 커버리지는 "
     "세트별 합집합과 교차 결합 포함 두 값으로 표시), 유전자 보존영역(자리별 &pi;에서 계산)."),
    ("7", "리니지별 유전체 다양성 (pixy &pi;)",
     "개별 세트의 mfeprimer 매치와는 별개로, population 게놈 전체를 mafft --adjustdirection으로 한 번에 "
     "정렬한 뒤(LASV S/L 게놈은 GenBank 등록 방향이 서로 달라 이 단계가 필요합니다) 참조유전체 좌표계로 "
     "정규화하고, 참조 행을 기준으로 all-sites VCF를 직접 만들어 pixy로 리니지별 nucleotide diversity"
     "(&pi;)를 1bp(염기 단위) 슬라이딩 윈도우로 계산합니다. amplicon 크기(보통 100&ndash;170bp)에 맞춰 "
     "가장 촘촘한 해상도로 계산해, 결합 위치가 보존 영역인지 정확히 어느 염기에서 변이가 있는지 세트별로 "
     "바로 확인할 수 있습니다."),
    ("8", "출력 (대시보드 / 리포트)",
     "위 모든 결과를 인터랙티브 웹 대시보드에서 세트별&middot;그룹별로 드릴다운하며 볼 수 있고(세트 상세 "
     "페이지에는 &pi; 오버레이 차트와 MSA 위 &pi; 히스토그램도 포함), 국가별 보고 현황을 보여주는 지도 "
     "페이지도 별도로 제공됩니다. 이 PDF처럼 정리된 리포트로도 내보낼 수 있습니다."),
]


def render_workflow_page():
    stages_html = "".join(f"""
      <div class="wf-stage">
        <div class="wf-num">{esc(n)}</div>
        <div class="wf-body"><div class="wf-title">{esc(title)}</div><div class="wf-desc">{desc}</div></div>
      </div>
      {'<div class="wf-arrow">&darr;</div>' if i < len(WORKFLOW_STAGES) - 1 else ''}"""
        for i, (n, title, desc) in enumerate(WORKFLOW_STAGES))
    return f"""
    <section class="report-page">
      <h1 style="font-size:18px;">분석 워크플로</h1>
      <p class="detail-sub">이 리포트가 어떤 순서로, 무엇을 근거로 계산됐는지 정리한 페이지입니다.</p>
      <div class="box" style="margin-bottom:12px;">
        <p style="margin:0;">LASV(라사바이러스) TaqMan primer/probe 세트가 Nextstrain의 실제 population
        유전체 데이터에 대해 얼마나 잘 매치되는지를 자동으로 검증하는 Snakemake 파이프라인의 결과입니다.
        설계한 프라이머/프로브가 실제 유전적 다양성(lineage별 변이 등)을 얼마나 커버하는지, 어떤 게놈에서
        mismatch가 발생하는지, 표준물질과는 맞는지, 여러 세트를 함께 쓸 때 의도치 않게 서로 결합하지는
        않는지를 확인합니다.</p>
      </div>
      <div class="wf-flow">{stages_html}</div>
      <div class="mini-table-note">사용 프로그램: Snakemake(워크플로 엔진) &middot; mfeprimer v4.2.4(결합 위치/mismatch/Tm/&Delta;G 계산 엔진)
      &middot; mafft(다중 서열 정렬) &middot; Nextstrain(공개 서열 및 계통 메타데이터 출처)</div>
    </section>
    """


def render_summary_page(sets_list, gene_by_set, groups_list):
    def std_pill(s):
        ok = (s.get("standard") or {}).get("ok")
        return "<span class='pill pill-ok'>O</span>" if ok else "<span class='pill pill-bad'>X</span>"

    set_rows = "".join(f"""
      <tr>
        <td>{esc(s.get('display_name') or s['set_id'])}</td>
        <td><span class="pill pill-{esc(s['segment'].lower())}">{esc(s['segment'])}</span></td>
        <td>{esc(s.get('source'))}</td>
        <td>{esc(gene_by_set.get(s['set_id'], '—'))}</td>
        <td>{fmt_pct((s.get('summary') or {}).get('validated_pct'))}</td>
        <td>{std_pill(s)}</td>
      </tr>""" for s in sets_list)

    group_rows = "".join(f"""
      <tr>
        <td>{esc(g.get('label'))}</td>
        <td><span class="pill pill-{esc(g['segment'].lower())}">{esc(g['segment'])}</span></td>
        <td>{g.get('member_count')}</td>
        <td>{fmt_pct((g.get('summary') or {}).get('validated_pct'))}</td>
      </tr>""" for g in groups_list)

    return f"""
    <section class="report-page">
      <h1 style="font-size:18px;">전체 요약</h1>
      <p class="detail-sub">{len(sets_list)}개 세트 &middot; {len(groups_list)}개 그룹의 핵심 지표를 한눈에 볼 수 있는
      요약표입니다. 상세 결과는 다음 페이지부터 세트/그룹별로 이어집니다.</p>
      <div class="box" style="margin-bottom:10px;">
        <h3>세트별 요약 (이름순 &mdash; 다음 페이지부터 이 순서대로 이어집니다)</h3>
        <table class="mini-table">
          <thead><tr><th>세트</th><th>segment</th><th>source</th><th>target gene</th>
          <th>Validated</th><th>표준물질</th></tr></thead>
          <tbody>{set_rows}</tbody>
        </table>
      </div>
      {f'''<div class="box">
        <h3>그룹(세트 조합 "ALL") 요약</h3>
        <table class="mini-table">
          <thead><tr><th>그룹</th><th>segment</th><th>구성 세트 수</th><th>Validated</th></tr></thead>
          <tbody>{group_rows}</tbody>
        </table>
      </div>''' if groups_list else ""}
    </section>
    """


def render_cover(overview_data, sets_count, generated_at, groups_count=0):
    ds = overview_data.get("dataset_sizes", {})
    return f"""
    <section class="report-page cover-page">
      <div class="confidential-banner">🔒 대외비 (CONFIDENTIAL) — 외부 공유 금지</div>
      <h1 class="cover-title">{VIRUS_SHORT} Primer/Probe Validation<br>종합 분석 리포트</h1>
      <div class="cover-meta">생성 일시: {esc(generated_at)}</div>
      <div class="cover-stats">
        <div class="stat"><div class="v">{esc(ds.get('S', '—'))}</div><div class="l">S segment 게놈</div></div>
        <div class="stat"><div class="v">{esc(ds.get('L', '—'))}</div><div class="l">L segment 게놈</div></div>
        <div class="stat"><div class="v">{sets_count}</div><div class="l">등록된 세트</div></div>
        <div class="stat"><div class="v">{groups_count}</div><div class="l">등록된 그룹</div></div>
      </div>
      <p class="cover-note">이 리포트의 프라이머/프로브 서열 및 분석 결과는 대외비 정보입니다.
      각 세트별 상세 결과는 다음 페이지부터 이어지며, 마지막에 그룹(여러 세트 조합) 결과가
      이어집니다.</p>
    </section>
    """


REPORT_CSS = """
@page { size: A4; margin: 16mm 14mm; }
* { box-sizing: border-box; }
body { font-family: -apple-system, "Noto Sans KR", "Malgun Gothic", sans-serif; color: #1c2430; font-size: 11px; line-height: 1.35; }
.mono { font-family: "Consolas", "Menlo", monospace; word-break: break-all; }
.report-page { page-break-after: always; }
.report-page:last-child { page-break-after: auto; }
.page-header { display: flex; align-items: center; gap: 10px; margin-bottom: 4px; }
.page-header h1 { font-size: 20px; margin: 0; }
.page-index { margin-left: auto; color: #888; font-size: 10px; }
.pill { display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 10px; font-weight: 700; background: #eef1f5; color: #445; }
.pill-s { background: #e3ecf7; color: #2f5fa8; }
.pill-l { background: #fbe8e7; color: #b23f3e; }
.detail-sub { color: #667; font-size: 11px; margin-bottom: 6px; }
.assessment-box { border-left: 4px solid #ccc; border-radius: 6px; padding: 6px 10px; margin-bottom: 8px; background: #fafbfc; }
.assessment-box.ok { border-left-color: #4a9c4a; }
.assessment-box.warn { border-left-color: #d99a3d; }
.assessment-box.bad { border-left-color: #c0504d; }
.assessment-verdict { font-weight: 800; font-size: 12.5px; margin-bottom: 6px; }
.assessment-points { margin: 0; padding-left: 18px; }
.assessment-points li { margin-bottom: 3px; }
.assessment-points li.ok { color: #2c7a2c; }
.assessment-points li.bad { color: #b23f3e; }
.grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; margin-bottom: 8px; }
.box { border: 1px solid #dde3ea; border-radius: 8px; padding: 8px 10px; margin-bottom: 8px; break-inside: avoid; }
.grid2 .box { margin-bottom: 0; }
.box h3 { font-size: 12px; margin: 0 0 6px; }
.seq-row { display: flex; gap: 8px; margin-bottom: 3px; font-size: 10px; }
.seq-row .tag { flex-shrink: 0; width: 60px; font-weight: 700; color: #556; }
.mini-table { width: 100%; border-collapse: collapse; font-size: 10.5px; }
.mini-table td, .mini-table th { padding: 2px 6px; border-bottom: 1px solid #eee; text-align: left; }
.mini-table tr.highlight td { background: #eef5ff; font-weight: 700; }
.mini-table td.bad, .mini-table td.mono.bad { color: #b23f3e; font-weight: 700; }
.mini-table td.warn { color: #b3792e; font-weight: 700; }
.mini-table td.ok { color: #3d7d34; font-weight: 700; }
.mini-table.small { font-size: 9.5px; }
.mini-table.small td, .mini-table.small th { padding: 2px 5px; }
.amp-legend { font-size: 9px; color: #889; margin: 0 0 8px; }
.amp-legend .amp-role-tag { font-size: 8px; padding: 0 4px; margin-right: 2px; }
.standard-material-block { border: 1px solid #dde3ea; border-radius: 8px; padding: 6px 10px; margin-bottom: 6px; break-inside: avoid; }
.standard-material-head { display: flex; flex-wrap: wrap; align-items: center; gap: 5px; font-weight: 700; font-size: 11px; margin-bottom: 6px; }
.mm-0 { background: rgba(84,162,74,.16); color: #3d7d34; }
.mm-1 { background: rgba(242,166,90,.20); color: #b5731a; }
.mm-2plus { background: rgba(228,87,86,.18); color: #b23f3e; }
.mm-na { color: #889; }
.amp-verdict { font-size: 11px; font-weight: 700; padding: 5px 8px; border-radius: 6px; margin-bottom: 6px; }
.amp-verdict.ok { background: rgba(84,162,74,.14); color: #3d7d34; }
.amp-verdict.bad { background: rgba(228,87,86,.14); color: #b23f3e; }
.amp-role-summary { display: flex; flex-direction: column; gap: 3px; margin-bottom: 6px; }
.amp-role-row { display: flex; align-items: center; flex-wrap: wrap; gap: 4px 8px; font-size: 9.5px; }
.amp-role-tag { font-weight: 700; padding: 1px 5px; border-radius: 4px; background: #eef1f5; }
.amp-role-tag.f { color: #4C78A8; }
.amp-role-tag.p { color: #b5731a; }
.amp-role-tag.r { color: #3d7d34; }
.amp-thermo { color: #889; font-size: 9px; }
.msa-mm-tag { font-size: 8.5px; padding: 0 3px; border-radius: 3px; font-weight: 700; }
.amp-align { font-family: "Consolas", "Menlo", monospace; font-size: 7.5px; letter-spacing: 0; line-height: 1.4; }
.amp-align-group { margin-bottom: 4px; }
.amp-align-group:last-child { margin-bottom: 0; }
.amp-align-row { white-space: pre; }
.amp-align-colon { color: #889; }
.amp-align-dir { color: #4C78A8; font-weight: 700; }
.amp-align mark.f { background: rgba(76,120,168,.35); color: inherit; }
.amp-align mark.r { background: rgba(84,162,74,.35); color: inherit; }
.amp-align mark.p { background: rgba(242,166,90,.45); color: inherit; }
.big-ox { font-size: 28px; font-weight: 800; }
.big-ox.ok { color: #2c7a2c; }
.big-ox.bad { color: #b23f3e; }
.muted { color: #889; font-size: 10.5px; }
.note-text { white-space: pre-wrap; font-size: 11px; background: #fafbfc; border-radius: 6px; padding: 8px 10px; }
.note-images { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 8px; }
.note-images figure { margin: 0; width: 140px; }
.note-images img { width: 100%; border-radius: 6px; border: 1px solid #dde3ea; }
.note-images figcaption { font-size: 9px; color: #889; margin-top: 2px; }
.confidential-banner { display: inline-block; background: #fbe8e7; color: #b23f3e; font-weight: 700; padding: 6px 14px; border-radius: 8px; margin-bottom: 40px; }
.cover-page { display: flex; flex-direction: column; align-items: center; justify-content: center; text-align: center; height: 240mm; }
.cover-title { font-size: 26px; line-height: 1.4; }
.cover-meta { color: #778; margin: 10px 0 30px; }
.cover-stats { display: flex; gap: 30px; margin-bottom: 30px; }
.cover-stats .stat .v { font-size: 26px; font-weight: 800; }
.cover-stats .stat .l { font-size: 11px; color: #778; }
.cover-note { max-width: 420px; color: #667; font-size: 11px; }
.wf-flow { display: flex; flex-direction: column; align-items: stretch; }
.wf-stage { display: flex; gap: 10px; border: 1px solid #dde3ea; border-radius: 8px; padding: 8px 10px; break-inside: avoid; }
.wf-num { flex-shrink: 0; width: 26px; height: 26px; border-radius: 50%; background: #eef1f5; color: #2f5fa8;
  font-weight: 800; font-size: 10.5px; display: flex; align-items: center; justify-content: center; }
.wf-title { font-weight: 700; font-size: 12px; margin-bottom: 3px; }
.wf-desc { font-size: 10.5px; color: #556; line-height: 1.5; }
.wf-arrow { text-align: center; color: #aab; font-size: 13px; margin: 2px 0; }
.mini-table-note { font-size: 10px; color: #889; margin-top: 10px; }
.seq-redacted { color: #aab; letter-spacing: 1px; }
.citation-box { margin-top: 8px; padding: 6px 9px; border: 1px dashed #dde3ea; border-radius: 6px; font-size: 9px; color: #778; }
.citation-box b { color: #445; }
.citation-box a { color: #2f5fa8; word-break: break-all; }
"""


def _ws_handshake(sock, host, port, path):
    key = base64.b64encode(b"lasvreportlasvrpt").decode()
    req = (
        f"GET {path} HTTP/1.1\r\nHost: {host}:{port}\r\nUpgrade: websocket\r\n"
        f"Connection: Upgrade\r\nSec-WebSocket-Key: {key}\r\nSec-WebSocket-Version: 13\r\n\r\n"
    )
    sock.sendall(req.encode())
    resp = b""
    while b"\r\n\r\n" not in resp:
        resp += sock.recv(4096)
    if b" 101 " not in resp.split(b"\r\n", 1)[0]:
        raise RuntimeError(f"CDP websocket handshake failed: {resp[:200]!r}")


def _ws_send(sock, payload: bytes):
    header = bytearray([0x81])
    length = len(payload)
    mask = b"\x00\x01\x02\x03"
    if length <= 125:
        header.append(0x80 | length)
    elif length <= 65535:
        header.append(0x80 | 126)
        header += struct.pack(">H", length)
    else:
        header.append(0x80 | 127)
        header += struct.pack(">Q", length)
    header += mask
    masked = bytes(b ^ mask[i % 4] for i, b in enumerate(payload))
    sock.sendall(bytes(header) + masked)


def _ws_recv(sock):
    def recv_exact(n):
        buf = b""
        while len(buf) < n:
            chunk = sock.recv(n - len(buf))
            if not chunk:
                raise ConnectionError("CDP websocket closed unexpectedly")
            buf += chunk
        return buf

    while True:
        b1, b2 = recv_exact(2)
        opcode = b1 & 0x0F
        masked = b2 & 0x80
        length = b2 & 0x7F
        if length == 126:
            length = struct.unpack(">H", recv_exact(2))[0]
        elif length == 127:
            length = struct.unpack(">Q", recv_exact(8))[0]
        mask_key = recv_exact(4) if masked else None
        payload = recv_exact(length)
        if mask_key:
            payload = bytes(b ^ mask_key[i % 4] for i, b in enumerate(payload))
        if opcode == 0x9:  # ping frame -- ignore and keep reading
            continue
        if opcode == 0x1:
            return payload


CHROME_CANDIDATES = (
    "google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
    "chromium-browser-stable", "/snap/bin/chromium",
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome", "/Applications/Chromium.app/Contents/MacOS/Chromium",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
)


def find_chrome(configured: str) -> str:
    """Resolves the headless-Chrome-compatible binary used for PDF export -- tries the
    explicit --chrome value first (if it's not just the argparse default), then a handful
    of common names/paths across distros, since "google-chrome" specifically isn't
    installed everywhere a Chromium-based browser is (e.g. Debian/Ubuntu's own chromium
    package ships as `chromium`/`chromium-browser` instead). Raises a clear, actionable
    message rather than letting a bare FileNotFoundError from subprocess.Popen surface."""
    import shutil
    if configured and configured != "google-chrome" and shutil.which(configured):
        return configured
    for name in CHROME_CANDIDATES:
        found = shutil.which(name) or (name if name.startswith("/") and Path(name).exists() else None)
        if found:
            return found
    sys.exit(
        "[report] PDF 생성에 필요한 Chrome/Chromium을 찾을 수 없습니다. "
        "google-chrome 또는 chromium을 설치한 뒤 다시 시도하거나(예: sudo apt install chromium-browser), "
        "--chrome <경로>로 직접 지정하세요. (대시보드의 다른 기능들에는 영향 없음 -- PDF 내보내기만 이 바이너리가 필요합니다.)"
    )


def print_to_pdf_via_cdp(chrome_bin, html_path: Path, out_pdf: Path, port=9333):
    """Render html_path to out_pdf using headless Chrome's DevTools Protocol directly
    (Page.printToPDF, displayHeaderFooter=False) -- the `--print-to-pdf-no-header` CLI
    flag doesn't suppress the header/footer on this Chrome build, but the CDP parameter
    does."""
    profile_dir = Path(tempfile.gettempdir()) / f"lasv-report-chrome-profile-{os.getpid()}"
    proc = subprocess.Popen([
        chrome_bin, "--headless", "--disable-gpu", "--no-sandbox",
        f"--remote-debugging-port={port}", "--remote-allow-origins=*",
        f"--user-data-dir={profile_dir}",
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(100):
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{port}/json/version", timeout=1)
                break
            except Exception:  # noqa: BLE001
                time.sleep(0.2)
        else:
            raise RuntimeError("headless Chrome DevTools endpoint never came up")

        target_url = html_path.resolve().as_uri()
        req = urllib.request.Request(f"http://127.0.0.1:{port}/json/new?{target_url}", method="PUT")
        with urllib.request.urlopen(req, timeout=10) as r:
            target = json.loads(r.read())
        ws_url = target["webSocketDebuggerUrl"]
        hostport, path = ws_url.split("://", 1)[1].split("/", 1)
        host, wport = hostport.split(":")
        path = "/" + path

        sock = socket.create_connection((host, int(wport)), timeout=15)
        try:
            _ws_handshake(sock, host, int(wport), path)
            _ws_send(sock, json.dumps({"id": 1, "method": "Page.enable"}).encode())
            _ws_recv(sock)
            time.sleep(1.0)  # let the already-navigated page (via /json/new?url) settle
            _ws_send(sock, json.dumps({
                "id": 2, "method": "Page.printToPDF",
                "params": {"printBackground": True, "displayHeaderFooter": False, "preferCSSPageSize": True},
            }).encode())
            while True:
                msg = json.loads(_ws_recv(sock))
                if msg.get("id") == 2:
                    break
            if "error" in msg:
                raise RuntimeError(f"Page.printToPDF failed: {msg['error']}")
            out_pdf.write_bytes(base64.b64decode(msg["result"]["data"]))
        finally:
            sock.close()
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:  # noqa: BLE001
            proc.kill()
        shutil.rmtree(profile_dir, ignore_errors=True)


VIRUS_SHORT = "LASV"   # replaced in main() by the profile's short_name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-url", default="http://localhost:8765")
    ap.add_argument("--out", default=str(ROOT / "output" / "report" / "report.pdf"))
    ap.add_argument("--chrome", default="google-chrome")
    args = ap.parse_args()
    global VIRUS_SHORT
    try:
        VIRUS_SHORT = fetch_json(args.base_url, "/api/profile")["short_name"]
    except Exception:  # noqa: BLE001 -- an older server without /api/profile: keep the default
        pass

    chrome_bin = find_chrome(args.chrome)  # fail fast -- before spending time fetching/rendering 25+ sets

    base_url = args.base_url.rstrip("/")
    print(f"[report] fetching data from {base_url} ...", file=sys.stderr)
    overview = fetch_json(base_url, "/api/overview")
    sets_list = fetch_json(base_url, "/api/sets")
    sets_list = sorted(sets_list, key=lambda s: s["set_id"])
    groups_list = fetch_json(base_url, "/api/groups")
    groups_list = sorted(groups_list, key=lambda g: g["group_id"])
    try:
        diversity_data = fetch_json(base_url, "/api/diversity")
    except Exception as e:  # noqa: BLE001
        print(f"[report] warning: /api/diversity unavailable ({e}), skipping pi column", file=sys.stderr)
        diversity_data = None

    sections = [render_cover(overview, len(sets_list), time.strftime("%Y-%m-%d %H:%M"), len(groups_list)),
                render_workflow_page()]
    summary_page_index = len(sections)
    sections.append("")  # placeholder -- filled in once gene_by_set is known, after the set loop below

    gene_by_set = {}
    for i, s in enumerate(sets_list, 1):
        print(f"[report] set {i}/{len(sets_list)} {s['set_id']}", file=sys.stderr)
        detail = fetch_json(base_url, f"/api/sets/{s['set_id']}")
        gene_by_set[s["set_id"]] = (detail.get("gene") or {}).get("target_gene") or "—"
        sections.append(render_set_section(base_url, detail, i, len(sets_list), diversity_data))
    sections[summary_page_index] = render_summary_page(sets_list, gene_by_set, groups_list)

    for i, g in enumerate(groups_list, 1):
        print(f"[report] group {i}/{len(groups_list)} {g['group_id']}", file=sys.stderr)
        detail = fetch_json(base_url, f"/api/groups/{g['group_id']}")
        sections.append(render_group_section(detail, i, len(groups_list)))

    html_doc = f"""<!doctype html><html lang="ko"><head><meta charset="utf-8">
<title>{VIRUS_SHORT} Primer/Probe Validation Report</title><style>{REPORT_CSS}</style></head>
<body>{''.join(sections)}</body></html>"""

    out_pdf = Path(args.out)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    out_html = out_pdf.with_suffix(".html")
    out_html.write_text(html_doc, encoding="utf-8")
    print(f"[report] wrote {out_html}", file=sys.stderr)

    print_to_pdf_via_cdp(chrome_bin, out_html, out_pdf)
    if not out_pdf.exists():
        sys.exit("[report] chrome print-to-pdf (CDP) failed: no output written")
    print(f"[report] wrote {out_pdf} ({out_pdf.stat().st_size} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
