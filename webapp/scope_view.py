"""Host-scope view of the results for the web API (scripts/host_scope.py defines the classes, the scopes and the counting).

A request may carry ?host=human|animal|unknown|all. Without it (or with anything else) the PRIMARY scope of each segment is used:
human genomes, or "all" when the metadata has no human genome. Laboratory hosts are in no scope. Everything below restricts what the
pipeline already computed to the genomes of the scope and recounts with the SAME code the pipeline uses (host_scope.summarize), so the
numbers agree with summary.tsv for the primary scope.

This file is shared verbatim by the LASV and Nipah web apps.
"""
import csv
import threading
from pathlib import Path

import host_scope as H
from mfe_common import probe_inside_amplicon

_REQ = threading.local()
_CACHE = {}


def set_request_scope(value):
    """Called once per request. Anything that is not one of host_scope.SCOPES means 'the primary scope'."""
    _REQ.scope = value if value in H.SCOPES else None


def requested_scope():
    return getattr(_REQ, "scope", None)


_SUMMARIES = {}


def _stamp(p):
    p = Path(p)
    return p.stat().st_mtime_ns if p.exists() else None


def _cached(path, extra, loader):
    """loader(path) result, re-read only when the file changes."""
    p = Path(path)
    stamp = p.stat().st_mtime_ns if p.exists() else None
    key = (str(p), repr(extra))
    hit = _CACHE.get(key)
    if hit and hit[0] == stamp:
        return hit[1]
    value = loader(p)
    _CACHE[key] = (stamp, value)
    return value


def _read_tsv(p):
    if not p.exists():
        return []
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class View:
    def __init__(self, root, columns, primer_mismatch_pass, probe_mismatch):
        self.root = Path(root)
        self.columns = dict(columns or {})
        self.mm = int(primer_mismatch_pass)
        self.pm = int(probe_mismatch)

    # ---- the population and the host class of each genome -----------------------------------------------------------------
    def classes(self, seg):
        return _cached(self.root / "data" / "raw" / f"{seg.lower()}_metadata.tsv", self.columns,
                       lambda p: H.load_host_classes(p, self.columns))

    def population(self, seg):
        """Rows of the filtered metadata (the population every set was matched against)."""
        return _cached(self.root / "data" / "filtered" / f"{seg.lower()}_filtered_metadata.tsv", None, _read_tsv)

    def pop_ids(self, seg):
        return [r["accession"] for r in self.population(seg)]

    def primary(self, seg):
        return H.primary_scope(self.pop_ids(seg), self.classes(seg))

    def scope(self, seg):
        return requested_scope() or self.primary(seg)

    def ids(self, seg, scope=None):
        return set(H.ids_in_scope(self.pop_ids(seg), self.classes(seg), scope or self.scope(seg)))

    def host_class(self, seg, accession):
        return self.classes(seg).get(accession, "unknown")

    def info(self, seg):
        """What the toggle needs: the active scope, the primary one and how many genomes each scope holds."""
        counts = H.host_counts(self.pop_ids(seg), self.classes(seg))
        n = {"human": counts["human"], "animal": counts["animal"], "unknown": counts["unknown"],
             "all": counts["human"] + counts["animal"] + counts["unknown"]}
        return {"scope": self.scope(seg), "primary": self.primary(seg), "counts": counts, "scope_counts": n,
                "excluded_lab": counts["lab"]}

    # ---- numbers --------------------------------------------------------------------------------------------------------------
    def summary_from_rows(self, seg, rows_by_id, scope=None):
        ids = sorted(self.ids(seg, scope))
        s = H.summarize(ids, rows_by_id, self.mm, self.pm, probe_inside_amplicon)
        return {k: v for k, v in s.items() if not k.startswith("_")}

    def _set_numbers(self, set_id, seg, scope):
        """H.summarize() of a set's mismatch table inside `scope`, kept until the table, the population/host files or the pass rule change
        (the overview, the Venn page and every group card ask for the same numbers again and again)."""
        path = self.root / "output" / "primers" / set_id / "alignment" / "mismatch_table.tsv"
        scope = scope or self.scope(seg)
        sp = lambda *a: self.root.joinpath(*a)  # noqa: E731
        key = (str(path), _stamp(path), seg.lower(), scope, self.mm, self.pm, repr(self.columns),
               _stamp(sp("data", "filtered", f"{seg.lower()}_filtered_metadata.tsv")), _stamp(sp("data", "raw", f"{seg.lower()}_metadata.tsv")))
        hit = _SUMMARIES.get(key)
        if hit is None:
            rows = _cached(path, None, lambda p: {r["genome_id"]: r for r in _read_tsv(p)})
            hit = H.summarize(sorted(self.ids(seg, scope)), rows, self.mm, self.pm, probe_inside_amplicon)
            if len(_SUMMARIES) > 600:
                _SUMMARIES.clear()
            _SUMMARIES[key] = hit
        return hit

    def set_summary(self, set_id, seg, scope=None):
        return {k: v for k, v in self._set_numbers(set_id, seg, scope).items() if not k.startswith("_")}

    def set_validated_ids(self, set_id, seg, scope=None):
        """The genomes a set validates, recounted from its mismatch table with the CURRENT parameters inside `scope` -- exactly the rule
        behind the set card's number (set_summary), so a figure built from these ids can never disagree with the card."""
        return set(self._set_numbers(set_id, seg, scope)["_validated_ids"])

    def mismatch_curve(self, set_id, seg, kmax=4, scope=None):
        """How many genomes of the scope a set validates when `k` mismatches are tolerated, k = 0..kmax, in two readings: `fixed` = F and R may differ
        in k positions (3' ends intact) while the probe keeps the pipeline's own setting, `same` = F, R and the probe all may differ in k positions.
        -> {"total", "fixed": [...], "same": [...], "by_lineage": {lineage: {"total", "fixed", "same"}}} (counts)"""
        path = self.root / "output" / "primers" / set_id / "alignment" / "mismatch_table.tsv"
        scope = scope or self.scope(seg)
        key = ("curve", str(path), _stamp(path), seg.lower(), scope, kmax, self.pm, repr(self.columns),
               _stamp(self.root / "data" / "filtered" / f"{seg.lower()}_filtered_metadata.tsv"), _stamp(self.root / "data" / "raw" / f"{seg.lower()}_metadata.tsv"))
        hit = _SUMMARIES.get(key)
        if hit is None:
            rows = _cached(path, None, lambda p: {r["genome_id"]: r for r in _read_tsv(p)})
            ids = sorted(self.ids(seg, scope))
            ok = lambda subset, a, b: H.summarize(subset, rows, a, b, probe_inside_amplicon)["validated"]  # noqa: E731
            curve = lambda subset: {"total": len(subset), "fixed": [ok(subset, k, self.pm) for k in range(kmax + 1)], "same": [ok(subset, k, k) for k in range(kmax + 1)]}  # noqa: E731
            lin = self.lineage_of(seg)
            by = {}
            for L in sorted({lin.get(i, "NA") for i in ids}):
                sub = [i for i in ids if lin.get(i, "NA") == L]
                if len(sub) >= 3:                                  # a lineage of one or two genomes makes no curve
                    by[L] = curve(sub)
            hit = dict(curve(ids), by_lineage=by)
            if len(_SUMMARIES) > 600:
                _SUMMARIES.clear()
            _SUMMARIES[key] = hit
        return hit

    def lineage_of(self, seg):
        return _cached(self.root / "data" / "filtered" / f"{seg.lower()}_filtered_metadata.tsv", "lineage",
                       lambda p: {r["accession"]: (r.get("lineage") or "NA") for r in _read_tsv(p)})

    def lineage_counts(self, seg, ids):
        lin = self.lineage_of(seg)
        out = {}
        for i in ids:
            k = lin.get(i, "NA")
            out[k] = out.get(k, 0) + 1
        return out

    def lineage_totals(self, seg, scope=None):
        return self.lineage_counts(seg, self.ids(seg, scope))

    # ---- restricting per-genome payloads ----------------------------------------------------------------------------------------
    def tag_and_filter(self, seg, rows, key="genome_id"):
        """Rows of the scope only, each with its host_class (the toggle then shows exactly what the numbers count)."""
        ids, cls = self.ids(seg), self.classes(seg)
        out = []
        for r in rows or []:
            if r.get(key) in ids:
                r["host_class"] = cls.get(r[key], "unknown")
                out.append(r)
        return out

    def tag_all(self, seg, rows, key="genome_id"):
        """Every row keeps its place, lab hosts included, with host_class (for tables that list the whole population)."""
        cls = self.classes(seg)
        for r in rows or []:
            r["host_class"] = cls.get(r.get(key), "unknown")
        return rows
