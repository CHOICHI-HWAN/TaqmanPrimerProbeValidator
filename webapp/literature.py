"""Primer/probe sets collected from the literature (the dashboard menu "문헌 프라이머" > LASV / Nipah).

The sequences live ONLY in config/literature/<virus>.tsv of the user's own computer (config/* is git-ignored); this module and the
menu ship with no data. Every entry must name its source (citation / DOI / PMID), and "verified" means: the sequences were compared
with the primary source, and `verified_note` says where. A new entry is "no" until someone has done that -- nothing here is
trusted just because it is in the list.

This file is shared verbatim by the LASV and Nipah web apps.
"""
import csv
import os
import re
from pathlib import Path

COLUMNS = ["lit_id", "name", "segment", "forward_primer", "reverse_primer", "probe", "probe_dyes", "target",
           "citation", "doi", "pmid", "verified", "verified_note", "remark"]
VIRUSES = {"lasv": {"label": "LASV", "segments": ("S", "L"), "prefix": "LIT-LASV-"},
           "nipah": {"label": "Nipah", "segments": ("GENOME",), "prefix": "LIT-NIPAH-"}}
SEQ_FIELDS = (("forward_primer", "Forward"), ("reverse_primer", "Reverse"), ("probe", "Probe"))
_BASES = set("ACGTUWSMKRYBDHVN")


def file_path(root, virus):
    return Path(root) / "config" / "literature" / f"{virus}.tsv"


def _clean(v):
    return re.sub(r"\s+", " ", str(v or "")).strip()


def _seq(v):
    return re.sub(r"[\s\-]", "", str(v or "")).upper().replace("U", "T")


def load(root, virus):
    """Rows of config/literature/<virus>.tsv (every column present, strings). Missing file = empty list."""
    p = file_path(root, virus)
    if virus not in VIRUSES or not p.exists():
        return []
    with open(p, newline="", encoding="utf-8") as fh:
        return [{c: (r.get(c) or "").strip() for c in COLUMNS} for r in csv.DictReader(fh, delimiter="\t")]


def next_id(rows, virus):
    pre = VIRUSES[virus]["prefix"]
    nums = [int(m.group(1)) for r in rows for m in [re.fullmatch(re.escape(pre) + r"(\d+)", r.get("lit_id", ""))] if m]
    return f"{pre}{(max(nums) + 1 if nums else 1):03d}"


def normalise(rows, virus):
    """Trimmed, upper-cased sequences, ids assigned to new rows, 'verified' reduced to yes/no. Fully empty rows are dropped."""
    out, work = [], [dict(r) for r in rows or []]
    for r in work:
        n = {c: _clean(r.get(c)) for c in COLUMNS}
        for f, _ in SEQ_FIELDS:
            n[f] = _seq(r.get(f))
        n["segment"] = n["segment"].upper()
        n["verified"] = "yes" if n["verified"].lower() in ("yes", "y", "true", "1") else "no"
        if not any(n[c] for c in ("name", "forward_primer", "reverse_primer", "probe", "citation")):
            continue
        out.append(n)
    for n in out:
        if not n["lit_id"]:
            n["lit_id"] = next_id(out, virus)
    return out


def validate(rows, virus):
    """[{row, field, message}] -- empty list = fine. Rows should already be normalised."""
    segs = VIRUSES[virus]["segments"]
    errors, names, ids = [], {}, {}

    def err(i, field, msg):
        errors.append({"row": i, "field": field, "message": msg})

    for i, r in enumerate(rows):
        if r["lit_id"] in ids:
            err(i, "lit_id", f"id가 중복됩니다 (행 {ids[r['lit_id']] + 1}과 동일)")
        ids[r["lit_id"]] = i
        if not r["name"]:
            err(i, "name", "이름이 비어 있습니다")
        elif r["name"] in names:
            err(i, "name", f"이름이 중복됩니다 (행 {names[r['name']] + 1}과 동일)")
        else:
            names[r["name"]] = i
        if r["segment"] not in segs:
            err(i, "segment", f"segment는 {' 또는 '.join(segs)} 이어야 합니다")
        for f, label in SEQ_FIELDS:
            bad = sorted(set(r[f]) - _BASES)
            if bad:
                err(i, f, f"{label} 서열에 올바르지 않은 문자가 있습니다: {', '.join(bad)}")
        if not any(r[f] for f, _ in SEQ_FIELDS):
            err(i, "forward_primer", "서열이 하나도 없습니다")
        if not (r["citation"] or r["doi"] or r["pmid"]):
            err(i, "citation", "출처(논문, DOI 또는 PMID)가 필요합니다 — 출처 없는 서열은 목록에 둘 수 없습니다")
        if r["doi"] and not re.fullmatch(r"10\.\d{4,9}/\S+", r["doi"]):
            err(i, "doi", "DOI 형식이 올바르지 않습니다 (10.xxxx/...)")
        if r["pmid"] and not r["pmid"].isdigit():
            err(i, "pmid", "PMID는 숫자여야 합니다")
        if r["verified"] == "yes" and not r["verified_note"]:
            err(i, "verified_note", "'원문 대조 완료'로 표시하려면 어디서 확인했는지(예: 논문 Table 2, 본문 p.3) 적어야 합니다")
    return errors


def save(root, virus, rows):
    p = file_path(root, virus)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLUMNS})
    os.replace(tmp, p)


def complete(row):
    """Can it be analysed at all? (the pipeline needs F, R and a probe)"""
    return all(row.get(f) for f, _ in SEQ_FIELDS)


def source_label(row):
    cite = row["citation"] or (f"doi:{row['doi']}" if row["doi"] else f"PMID {row['pmid']}")
    return f"문헌: {cite}"


def manifest_row(row):
    """A row for config/primers.tsv (the analysis list). The remark says whether the sequences were checked against the paper."""
    check = "원문 대조 완료" if row["verified"] == "yes" else "원문 대조 안 됨 — 서열 확인 필요"
    bits = [check, row["lit_id"]]
    if row["doi"]:
        bits.append("doi:" + row["doi"])
    if row["pmid"]:
        bits.append("PMID " + row["pmid"])
    if row["remark"]:
        bits.append(row["remark"])
    return {"primer_name": row["name"], "segment": row["segment"], "forward_primer": row["forward_primer"],
            "reverse_primer": row["reverse_primer"], "probe": row["probe"], "fluorophore": row["probe_dyes"],
            "source": source_label(row), "remark": " · ".join(bits), "proceed": "True"}


def same_sequences(manifest_rows, row):
    """The analysis list already holds this exact F/R/probe combination (any name)."""
    key = tuple(_seq(row.get(f)) for f, _ in SEQ_FIELDS)
    return any(tuple(_seq(m.get(f)) for f, _ in SEQ_FIELDS) == key for m in manifest_rows)
