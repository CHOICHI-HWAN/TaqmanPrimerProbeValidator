#!/usr/bin/env python3
"""The summary report as an editable Hancom HWPX document (the PDF's text content: assessment, sequences, coverage, binding, lineage, curve values).

Same data and same wording as scripts/generate_report.py (it is imported, so the verdicts cannot drift); the layout is a plain document
(headings, bullet lists, tables), not a picture of the PDF -- the coverage curve becomes a value table, the amplicon alignment drawings are left out.

    python3 scripts/report_hwpx.py [--base-url http://localhost:8765] [--out output/report/report.hwpx] [--unmasked]

Primer/probe sequences are masked (black blocks) unless --unmasked is given, exactly like the PDF. Needs the `python-hwpx` package
(pip install python-hwpx, Apache-2.0); without it the script exits with a message and the PDF report is unaffected.
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generate_report as gr  # noqa: E402


def cell(v):
    """One table cell: no pipes / line breaks (they would break the Markdown table)."""
    s = "—" if v is None or v == "" else str(v)
    return s.replace("|", "/").replace("\r", " ").replace("\n", " ").strip()


def table(header, rows):
    if not rows:
        return ""
    out = ["| " + " | ".join(cell(h) for h in header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(cell(c) for c in r) + " |" for r in rows]
    return "\n".join(out) + "\n"


def bullets(points):
    mark = {True: "[O]", False: "[!]", None: "[i]"}
    return "".join(f"{mark.get(ok, '[i]')} {cell(text)}\n\n" for ok, text in points)   # plain paragraphs: real HWPX list styles are not needed here


def seq_text(seq):
    seq = seq or ""
    if not seq:
        return "—"
    return "█" * len(seq) if not gr.UNMASK else seq


def lineage_rows(rows):
    return [(lin, f"{c} / {tot if tot is not None else '?'}", f"{pct:.0f}%" if pct is not None else "—",
             f"n<{gr.LOW_N_THRESHOLD}" if tot is not None and tot < gr.LOW_N_THRESHOLD else "") for lin, c, tot, pct in rows]


def curve_rows(base_url, set_data):
    seg = str(set_data.get("segment") or "").lower()
    try:
        if seg not in gr._CURVES:
            gr._CURVES[seg] = gr.fetch_json(base_url, f"/api/mismatch-curve/{seg}")
        c = gr._CURVES[seg]
    except Exception:  # noqa: BLE001 -- older server: no curve
        return None
    s = next((x for x in c.get("series") or [] if x.get("set_id") == set_data.get("set_id")), None)
    ks, vals, total = c.get("x") or [], (s or {}).get("fixed") or [], c.get("total") or 0
    if not s or not ks or not total:
        return None
    pp = c.get("primer_pass")
    return [(f"k = {k}" + (" (현재 기준)" if k == pp else ""), f"{vals[i] if i < len(vals) else 0} / {total}",
             f"{100.0 * (vals[i] if i < len(vals) else 0) / total:.1f}%") for i, k in enumerate(ks)]


def set_md(base_url, d, i, n, diversity):
    sid = d["set_id"]
    s, std, gene = d.get("summary") or {}, d.get("standard") or {}, d.get("gene") or {}
    verdict, _cls, points = gr.build_assessment(d)
    name = d.get("display_name") or sid
    md = [f"### {i} / {n}. {cell(name)} ({cell(d.get('segment'))}, {cell(d.get('source'))})\n"]
    amp = (f" · amplicon {gene.get('chrom', '')}:{gene.get('amplicon_start')}..{gene.get('amplicon_end')} ({gene.get('amplicon_size')}bp)"
           if gene.get("amplicon_start") else "")
    md.append(f"ID {sid} · Target gene: {cell(gene.get('target_gene', 'NA'))}{amp}\n")
    md.append(f"**{cell(verdict)}**\n")
    md.append(bullets(points))
    seq_rows = [("Forward", seq_text(d.get("F_seq")), len(d.get("F_seq") or "")),
                ("Probe", seq_text(d.get("P_seq")), len(d.get("P_seq") or "")),
                ("Reverse", seq_text(d.get("R_seq")), len(d.get("R_seq") or ""))]
    md.append("#### 프라이머 / 프로브 서열\n")
    md.append(table(["구분", "서열 (5'→3')", "길이(nt)"], seq_rows))
    if not gr.UNMASK:
        md.append("서열은 리포트에서 가렸습니다 — 대시보드에서 확인하세요.\n")
    cit = (d.get("citation") or {}).get("text")
    if cit:
        md.append(f"출처: {cell(cit)}\n")
    md.append("#### Population 커버리지\n")
    md.append(table(["항목", "값"], [
        ("Primer mismatch=0", gr.fmt_pct(s.get("primer_coverage_mm0_pct"))),
        ("Primer mismatch≤1", gr.fmt_pct(s.get("primer_coverage_mm1_pct"))),
        ("Primer mismatch≤2", gr.fmt_pct(s.get("primer_coverage_mm2_pct"))),
        ("Probe 완전일치", gr.fmt_pct(s.get("probe_coverage_mm0_pct"))),
        ("Probe dropout", s.get("probe_dropout")),
        ("최종 Validated", f"{s.get('validated')} / {s.get('total_genomes')} ({gr.fmt_pct(s.get('validated_pct'))}){gr.scope_note(s)}")]))
    cr = curve_rows(base_url, d)
    if cr:
        md.append("커버리지 곡선 값 (허용 mismatch 수 k → 통과 유전체):\n")
        md.append(table(["허용 mismatch", "통과 / 전체", "%"], cr))
    md.append("#### 표준물질(KDCA) 매치\n")
    md.append(f"{'O' if std.get('ok') else 'X'} — {cell(', '.join(std.get('matched') or []) or '매치 없음')}\n")
    md.append("#### 참조유전체 결합위치\n")
    brows = []
    for b in d.get("binding") or []:
        pi = gr._pi_for_range(diversity, d.get("segment"), b.get("start"), b.get("end")) if diversity else None
        if not b.get("found"):
            brows.append((b.get("role"), "NO MATCH", "", ""))
        else:
            brows.append((b.get("role"), f"{b.get('chrom')}:{b.get('start')}-{b.get('end')} ({b.get('strand')})",
                          f"mismatch {b.get('mismatches')}", f"π {pi:.3f}" if pi is not None else "π —"))
    md.append(table(["구분", "위치", "mismatch", "π (다양성)"], brows) if brows else "참조유전체 결합 정보 없음\n")
    md.append("#### 검증 통과 게놈의 리니지 분포\n")
    lr = gr.lineage_distribution(d.get("genomes") or [], d.get("lineage_totals") or {})
    md.append(table(["Lineage", "검출 / 전체 게놈 수", "검출율", "비고"], lineage_rows(lr)) if lr else "없음\n")
    if d.get("note"):
        md.append("#### 메모 / 노트\n")
        md.append(cell(d["note"]) + "\n")
    return "\n".join(md)


def group_md(g, i, n):
    summary = g.get("summary") or {}
    details = g.get("member_details") or []
    std_ok = any((m.get("standard") or {}).get("ok") for m in details)
    matched = sorted({x for m in details if (m.get("standard") or {}).get("ok") for x in ((m.get("standard") or {}).get("matched") or [])})
    verdict, _c, points = gr.build_assessment({"summary": {"validated_pct": summary.get("validated_pct")}, "standard": {"ok": std_ok, "matched": matched},
                                              "gene": {}, "binding": [], "dimer_hairpin": g.get("dimer_hairpin"), "cross_binding": g.get("combo_genomes") or []})
    counts, totals = g.get("lineage_counts") or {}, g.get("lineage_totals") or {}
    rows = []
    for lin in set(counts) | set(totals):
        c, tot = counts.get(lin, 0), totals.get(lin)
        rows.append((lin, c, tot, (100 * c / tot) if tot else None))
    rows.sort(key=lambda r: -r[1])
    names = ", ".join((g.get("member_names") or {}).get(m, m) for m in (g.get("members") or []))
    md = [f"### 그룹 {i} / {n}. {cell(g.get('label'))} ({cell(g.get('segment'))})\n", f"구성 세트: {cell(names)}\n",
          f"**합집합 커버리지: {summary.get('validated')} / {summary.get('total_genomes')}{gr.scope_note(summary)} ({gr.fmt_pct(summary.get('validated_pct'))}) — {cell(verdict)}**\n",
          "이 그룹의 세트 중 하나라도 통과(validated)하면 그 게놈은 합집합에 포함됩니다.\n", bullets(points)]
    if rows:
        md.append(table(["Lineage", "검출 / 전체 게놈 수", "검출율", "비고"], lineage_rows(rows)))
    return "\n".join(md)


def build_markdown(base_url, unmasked):
    gr.UNMASK = unmasked
    try:
        gr.PROFILE_INFO.update(gr.fetch_json(base_url, "/api/profile"))
        gr.VIRUS_SHORT = gr.PROFILE_INFO["short_name"]
    except Exception:  # noqa: BLE001
        pass
    overview = gr.fetch_json(base_url, "/api/overview")
    sets_list = sorted(gr.fetch_json(base_url, "/api/sets"), key=lambda s: s["set_id"])
    groups_list = sorted(gr.fetch_json(base_url, "/api/groups"), key=lambda g: g["group_id"])
    try:
        diversity = gr.fetch_json(base_url, "/api/diversity")
    except Exception:  # noqa: BLE001
        diversity = None
    ds = overview.get("dataset_sizes", {})
    md = [f"# {gr.VIRUS_SHORT} Primer/Probe Validation 종합 분석 리포트\n",
          f"생성 일시: {time.strftime('%Y-%m-%d %H:%M')}\n",
          "ⓘ in silico 예측 — 실험으로 확인 필요 (결과는 서열 기반 컴퓨터 예측이며 실제 결합·진단 성능을 확정하지 않습니다). "
          "외부에 공유하기 전에 서열과 결과의 공개 가능 여부를 확인하세요.\n",
          table(["게놈 수 (segment)", "값"], [(f"{k} segment 게놈", v) for k, v in ds.items()] + [("등록된 세트", len(sets_list)), ("등록된 그룹", len(groups_list))])]
    details, gene_by_set = [], {}
    for s in sets_list:
        print(f"[hwpx] set {s['set_id']}", file=sys.stderr)
        d = gr.fetch_json(base_url, f"/api/sets/{s['set_id']}")
        details.append(d)
        gene_by_set[s["set_id"]] = (d.get("gene") or {}).get("target_gene") or "—"
    md.append("## 전체 요약\n")
    md.append(table(["세트", "segment", "source", "target gene", "Validated", "표준물질"],
                    [(s.get("display_name") or s["set_id"], s.get("segment"), s.get("source"), gene_by_set.get(s["set_id"]),
                      gr.fmt_pct((s.get("summary") or {}).get("validated_pct")), "O" if (s.get("standard") or {}).get("ok") else "X") for s in sets_list]))
    if groups_list:
        md.append(table(["그룹", "segment", "구성 세트 수", "Validated"],
                        [(g.get("label"), g.get("segment"), g.get("member_count"), gr.fmt_pct((g.get("summary") or {}).get("validated_pct"))) for g in groups_list]))
    md.append("## 세트별 상세\n")
    for i, d in enumerate(details, 1):
        md.append(set_md(base_url, d, i, len(details), diversity))
    if groups_list:
        md.append("## 그룹(세트 조합) 결과\n")
        for i, g in enumerate(groups_list, 1):
            md.append(group_md(gr.fetch_json(base_url, f"/api/groups/{g['group_id']}"), i, len(groups_list)))
    return "\n".join(md)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-url", default="http://localhost:8765")
    ap.add_argument("--out", default=str(gr.ROOT / "output" / "report" / "report.hwpx"))
    ap.add_argument("--unmasked", action="store_true", help="write the real primer/probe sequences (default: black blocks)")
    args = ap.parse_args()
    try:
        from hwpx.tools.markdown_import import from_markdown
    except ImportError:
        sys.exit("[hwpx] the python-hwpx package is not installed (pip install python-hwpx) -- no HWPX written; the PDF report is not affected")
    md = build_markdown(args.base_url.rstrip("/"), args.unmasked)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    from_markdown(md).save_to_path(str(out))
    print(f"[hwpx] wrote {out} ({out.stat().st_size} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
