#!/usr/bin/env python3
"""Regenerates docs/REQUIREMENTS_LOG.md from docs/requirements_log.json (the data the web app's
"요구사항 기록" page also reads). Edit the JSON, then run this."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LABEL = {"done": "완료", "partial": "일부 미검증", "answer": "답변·제안", "dropped": "철회·불필요", "open": "진행 중"}


def main():
    d = json.loads((ROOT / "docs" / "requirements_log.json").read_text(encoding="utf-8"))
    items = [i for p in d["phases"] for i in p["items"]]
    counts = {k: sum(1 for i in items if i["status"] == k) for k in LABEL}
    esc = lambda s: (s or "").replace("|", "\\|")  # noqa: E731
    md = [f"# {d['title']} (프로젝트 시작 ~ 현재)\n",
          f"갱신일 {d['updated']} · 총 {len(items)}건 — " + " · ".join(f"{LABEL[k]} {counts[k]}" for k in LABEL if counts[k]) + "\n",
          "> 출처: 프로젝트 시작(2026-09-14)부터의 대화 기록에서 사용자가 요구한 사항을 시간순으로 정리했습니다. 상태는 당시 세션의 처리 결과와 CHANGELOG, 점검 결과 기준이며 오래된 항목을 하나하나 다시 실행해 확인한 것은 아닙니다. 요구 사항 칸의 괄호 번호는 사용자가 직접 붙인 번호입니다.\n"]
    n = 0
    for p in d["phases"]:
        md += ["", f"## {p['title']}", "", "| # | 날짜 | 구분 | 요구 사항 | 검토 · 수정 내용 | 상태 | 비고 |", "|---|---|---|---|---|---|---|"]
        for i in p["items"]:
            n += 1
            md.append(f"| {n} | {i['date']} | {i['cat']} | {esc(i['req'])} | {esc(i['rev'])} | {LABEL[i['status']]} | {esc(i.get('note'))} |")
    md += ["", "## 남은 확인 사항", "", "| 항목 | 내용 |", "|---|---|"] + [f"| {o['item']} | {o['note']} |" for o in d["open"]]
    (ROOT / "docs" / "REQUIREMENTS_LOG.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"{len(items)} items -> docs/REQUIREMENTS_LOG.md")


if __name__ == "__main__":
    main()
