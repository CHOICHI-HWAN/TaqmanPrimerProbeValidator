#!/usr/bin/env python3
"""Step 6: render the final HTML dashboard (charts + tables) from the pipeline outputs."""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from svg_charts import grouped_hbar_chart, simple_bar_chart, heatmap_table  # noqa: E402
from genome_track import genome_track_svg  # noqa: E402

SEG_COLOR = {"S": "#4C78A8", "L": "#E45756"}
PRIMER_COLOR = "#B0B7C1"
PROBE_COLOR = "#F2A65A"
VALIDATED_COLOR = "#4C78A8"


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_filter_report(path):
    text = Path(path).read_text(encoding="utf-8")
    m = re.search(r"kept=(\d+)", text)
    return int(m.group(1)) if m else 0


def read_single_fasta_seq(path):
    seq_lines = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                seq_lines.append(line.strip())
    return "".join(seq_lines)


def main():
    summary = read_tsv(snakemake.input.summary)
    lineage_rows = read_tsv(snakemake.input.summary_by_lineage)
    kept_by_seg = {
        "S": parse_filter_report(snakemake.input.filter_report_s),
        "L": parse_filter_report(snakemake.input.filter_report_l),
    }
    primer_mismatch_pass = snakemake.params.primer_mismatch_pass
    probe_mismatch = snakemake.params.probe_mismatch

    names = dict(snakemake.params.set_name)
    virus_short = snakemake.params.virus_short
    ref_accessions = "/".join(snakemake.params.refgenome_accession.values())
    summary.sort(key=lambda r: -float(r["validated_pct"]))
    categories = [f"{names.get(r['set_id'], r['set_id'])} ({r['segment']})" for r in summary]
    primer_label = f"Primer match (mm≤{primer_mismatch_pass})"
    probe_label = f"Probe perfect match (mm={probe_mismatch})"
    series = {
        primer_label: [float(r[f"primer_coverage_mm{primer_mismatch_pass}_pct"]) for r in summary],
        probe_label: [float(r["probe_coverage_mm0_pct"]) for r in summary],
        "Validated (final)": [float(r["validated_pct"]) for r in summary],
    }
    series_colors = {primer_label: PRIMER_COLOR, probe_label: PROBE_COLOR, "Validated (final)": VALIDATED_COLOR}
    chart_svg = grouped_hbar_chart(
        categories, series, series_colors,
        row_height=58, title="Set-by-set coverage: primer vs. probe vs. final validated",
    )

    dataset_svg = simple_bar_chart(
        ["S segment", "L segment"], [kept_by_seg["S"], kept_by_seg["L"]],
        colors=[SEG_COLOR["S"], SEG_COLOR["L"]], unit="", width=360, height=220,
        title="Curated database size (post N/length filter)",
    )

    # lineage heatmap: rows = set_id sorted as in summary, cols = union of lineages (by total desc)
    lineage_totals = {}
    per_set_lineage = {}
    for row in lineage_rows:
        per_set_lineage.setdefault(row["set_id"], {})[row["lineage"]] = int(row["validated_count"])
        lineage_totals[row["lineage"]] = lineage_totals.get(row["lineage"], 0) + int(row["validated_count"])
    lineages_sorted = [l for l, _ in sorted(lineage_totals.items(), key=lambda x: -x[1])]
    heat_rows = [f"{names.get(r['set_id'], r['set_id'])} ({r['segment']})" for r in summary]
    heat_matrix = [[per_set_lineage.get(r["set_id"], {}).get(l, 0) for l in lineages_sorted] for r in summary]
    heatmap_html = heatmap_table(heat_rows, lineages_sorted, heat_matrix)

    def fmt_pct(v):
        return f"{float(v):.1f}%"

    table_rows_html = []
    for r in summary:
        table_rows_html.append(
            "<tr>"
            f"<td class='mono'>{names.get(r['set_id'], r['set_id'])}</td>"
            f"<td><span class='pill pill-{r['segment']}'>{r['segment']}</span></td>"
            f"<td>{r['source']}</td>"
            f"<td>{r['total_genomes']}</td>"
            f"<td>{fmt_pct(r['primer_coverage_mm0_pct'])}</td>"
            f"<td>{fmt_pct(r['primer_coverage_mm1_pct'])}</td>"
            f"<td>{fmt_pct(r['primer_coverage_mm2_pct'])}</td>"
            f"<td>{fmt_pct(r['probe_coverage_mm0_pct'])}</td>"
            f"<td>{r['probe_dropout']}</td>"
            f"<td class='validated'>{fmt_pct(r['validated_pct'])}</td>"
            "</tr>"
        )

    # ---- combined-panel ("<<ALL" / "<<P1+P2") union rows ----
    groups = read_tsv(snakemake.input.group_summary)
    group_rows_html = []
    for g in groups:
        group_rows_html.append(
            "<tr>"
            f"<td class='mono'>{g['label']}</td>"
            f"<td><span class='pill pill-{g['segment']}'>{g['segment']}</span></td>"
            f"<td class='mono' style='font-size:11px;'>{g['members'].replace(';', ', ')}</td>"
            f"<td>{g['total_genomes']}</td>"
            f"<td class='validated'>{g['validated']} ({fmt_pct(g['validated_pct'])})</td>"
            "</tr>"
        )

    # ---- 표준물질(KDCA standard-material) O/X match table ----
    standard_summary = {r["set_id"]: r for r in read_tsv(snakemake.input.standard_summary)}
    standard_matches = read_tsv(snakemake.input.standard_matches)
    matched_standards_by_set = {}
    for row in standard_matches:
        if int(row["validated_count"]) > 0:
            matched_standards_by_set.setdefault(row["set_id"], []).append(row["lineage"])

    standard_rows_html = []
    for r in summary:  # keep same validated_pct-sorted order
        sid = r["set_id"]
        srow = standard_summary.get(sid)
        ok = srow is not None and int(srow["validated"]) > 0
        matched = ", ".join(sorted(matched_standards_by_set.get(sid, [])))
        if ok:
            perfect_primer = srow and float(srow["primer_coverage_mm0_pct"]) > 0
            qualifier = "완전일치" if perfect_primer else f"primer mismatch&le;{primer_mismatch_pass}"
            status_html = f"<span class='ox ox-o'>O</span> {matched} <span class='note' style='display:inline'>({qualifier})</span>"
        else:
            status_html = "<span class='ox ox-x'>X</span>"
        standard_rows_html.append(
            "<tr>"
            f"<td class='mono'>{names.get(sid, sid)}</td>"
            f"<td><span class='pill pill-{r['segment']}'>{r['segment']}</span></td>"
            f"<td>{status_html}</td>"
            "</tr>"
        )

    # ---- reference-genome (AY628203/AY628202) per-oligo binding-site maps ----
    refgenome_seq = {
        "s": read_single_fasta_seq(snakemake.input.refgenome_s),
        "l": read_single_fasta_seq(snakemake.input.refgenome_l),
    }
    refgenome_accession = snakemake.params.refgenome_accession
    set_seg = snakemake.params.set_seg
    binding_svgs_html = []
    for path in snakemake.input.refgenome_bindings:
        set_id = Path(path).parent.parent.name  # output/primers/<id>/refgenome_binding/binding.tsv
        seg = set_seg[set_id]
        chrom_label = refgenome_accession[seg]
        rows = read_tsv(path)
        hits = []
        for row in rows:
            found = row["found"] == "1"
            hits.append({
                "role": row["role"], "oligo_id": row["oligo_id"], "found": found,
                "chrom_label": chrom_label,
                "start": int(row["start"]) if found else None,
                "end": int(row["end"]) if found else None,
                "strand": row["strand"] if found else None,
                "mismatches": int(row["mismatches"]) if found else None,
            })
        svg = genome_track_svg(snakemake.params.set_name.get(set_id, set_id), refgenome_seq[seg], hits)
        binding_svgs_html.append(
            f"<div class='track-card'>{svg}</div>"
        )

    html = f"""<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<title>{virus_short} Primer/Probe Validation Dashboard</title>
<style>
  :root {{
    --bg: #ffffff; --fg: #1a1d23; --muted: #6b7280; --border: #e5e7eb;
    --card: #f8fafc; --accent: #4C78A8;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{ --bg: #14161a; --fg: #e7e9ee; --muted: #9aa3af; --border: #2a2e35; --card: #1b1e24; }}
  }}
  * {{ box-sizing: border-box; }}
  body {{ background: var(--bg); color: var(--fg); font-family: -apple-system, "Segoe UI", Pretendard, sans-serif;
          margin: 0; padding: 24px 20px 60px; line-height: 1.5; }}
  h1 {{ font-size: 22px; margin: 0 0 4px; }}
  .subtitle {{ color: var(--muted); font-size: 13px; margin-bottom: 24px; }}
  .grid {{ display: grid; grid-template-columns: 2fr 1fr; gap: 16px; margin-bottom: 24px; }}
  @media (max-width: 900px) {{ .grid {{ grid-template-columns: 1fr; }} }}
  .card {{ background: var(--card); border: 1px solid var(--border); border-radius: 12px; padding: 16px; overflow-x: auto; }}
  .card h2 {{ font-size: 15px; margin: 0 0 12px; }}
  svg {{ max-width: 100%; height: auto; color: var(--fg); }}
  table {{ border-collapse: collapse; width: 100%; font-size: 12.5px; }}
  th, td {{ padding: 6px 10px; border-bottom: 1px solid var(--border); text-align: right; white-space: nowrap; }}
  th {{ color: var(--muted); font-weight: 600; text-align: right; }}
  td.mono, th:first-child, td:first-child {{ text-align: left; font-family: ui-monospace, monospace; }}
  td.validated {{ font-weight: 700; color: var(--accent); }}
  .pill {{ display:inline-block; padding: 1px 8px; border-radius: 999px; font-size: 11px; font-weight: 600; }}
  .pill-S {{ background: rgba(76,120,168,.15); color: #4C78A8; }}
  .pill-L {{ background: rgba(228,87,86,.15); color: #E45756; }}
  .heatmap td, .heatmap th {{ text-align: center; padding: 5px 8px; }}
  .note {{ color: var(--muted); font-size: 12px; margin-top: 8px; }}
  footer {{ color: var(--muted); font-size: 11.5px; margin-top: 30px; border-top: 1px solid var(--border); padding-top: 12px; }}
  .ox {{ display:inline-block; width:18px; height:18px; line-height:18px; text-align:center;
         border-radius:5px; font-weight:800; font-size:12px; margin-right:6px; }}
  .ox-o {{ background: rgba(84,162,74,.18); color: #54A24B; }}
  .ox-x {{ background: rgba(228,87,86,.18); color: #E45756; }}
  .track-grid {{ display: grid; grid-template-columns: 1fr; gap: 12px; }}
  .track-card {{ background: var(--card); border: 1px solid var(--border); border-radius: 10px;
                 padding: 10px 14px; overflow-x: auto; }}
  section.card + section.card, .card + .card {{ margin-top: 24px; }}
</style>
</head>
<body>
  <h1>{virus_short} Primer / Probe Validation Dashboard</h1>
  <div class="subtitle">
    Nextstrain(lassa/s, lassa/l) 최신 데이터셋 기준 &middot;
    S segment: {kept_by_seg['S']}개 게놈 &middot; L segment: {kept_by_seg['L']}개 게놈 &middot;
    판정 기준: primer mismatch&le;{primer_mismatch_pass}, probe mismatch={probe_mismatch}(완전일치)만 인정
  </div>

  <div class="grid">
    <div class="card">
      <h2>세트별 커버리지 (primer만 vs probe 완전일치 vs 최종 통과)</h2>
      {chart_svg}
      <div class="note">Validated = primer 매치 게놈 ∩ probe 완전일치(mismatch=0) 게놈. probe에 mismatch가 있으면 primer가 매치해도 최종 통과에서 제외됩니다.</div>
    </div>
    <div class="card">
      <h2>큐레이션된 데이터셋 규모</h2>
      {dataset_svg}
      <div class="note">N 미포함 + 길이 필터(S&ge;3000bp, L&ge;6000bp) 통과 게놈 수 (아래 커버리지는 사람 유래 게놈 기준, 숙주 정보가 없으면 전체)</div>
    </div>
  </div>

  <div class="card" style="margin-bottom:24px;">
    <h2>세트별 상세 요약 (validated_pct 내림차순)</h2>
    <table>
      <thead><tr>
        <th>primer</th><th>segment</th><th>source</th><th>total</th>
        <th>primer mm=0</th><th>primer mm=1</th><th>primer mm=2</th>
        <th>probe mm=0</th><th>probe dropout</th><th>validated</th>
      </tr></thead>
      <tbody>
        {''.join(table_rows_html)}
      </tbody>
    </table>
  </div>

  <div class="card" style="margin-bottom:24px;">
    <h2>세트 조합("ALL") 커버리지 &mdash; 같은 벤더의 여러 세트를 함께 쓸 때 합집합</h2>
    <table>
      <thead><tr><th>조합</th><th>segment</th><th>구성 세트</th><th>total</th><th>validated (합집합)</th></tr></thead>
      <tbody>{''.join(group_rows_html)}</tbody>
    </table>
  </div>

  <div class="card" style="margin-bottom:24px;">
    <h2>표준물질(KDCA 2024) 매치여부</h2>
    <table>
      <thead><tr><th>primer</th><th>segment</th><th>매치 (일치한 표준물질)</th></tr></thead>
      <tbody>{''.join(standard_rows_html)}</tbody>
    </table>
    <div class="note">O = primer(mismatch&le;{primer_mismatch_pass}) AND probe(mismatch=0) 둘 다 표준물질 서열(800bp 구간) 안에서 발견됨.
    표준물질은 유전체 일부 구간(NP1/NP2/GC1/GC2: N/GPC 유전자, 800bp, L1/L2: L 유전자, 800bp)만 담고 있어서,
    세트의 타겟 위치가 그 구간 밖이면 X로 표시됩니다.</div>
  </div>

  <div class="card">
    <h2>최종 통과 게놈의 리니지 분포 (세트 &times; 리니지)</h2>
    {heatmap_html}
  </div>

  <div class="card" style="margin-top:24px;">
    <h2>참조유전체({ref_accessions}) 결합위치 맵</h2>
    <div class="note" style="margin-top:0;margin-bottom:10px;">각 세트의 F/Probe/R 프라이머가 참조유전체에서 어디에 결합하는지, 몇 개의 미스매치로 결합하는지 보여줍니다. 결합 위치를 찾지 못하면 빨간 점선 박스로 표시됩니다.</div>
    <div class="track-grid">
      {''.join(binding_svgs_html)}
    </div>
  </div>

  <footer>
    Generated by the {virus_short} primer/probe validation Snakemake pipeline (v2). Source data: nextstrain.org/lassa/s, /l.
    mfeprimer(--virus mode) for primer-pair coverage; mfeprimer search for probe perfect-match coverage.
  </footer>
</body>
</html>
"""
    Path(snakemake.output.html).write_text(html, encoding="utf-8")


if __name__ == "__main__":
    main()
