"""Regression tests for the spreadsheet-format primer manifest (config/primers.tsv), the
permanent id registry, and the input_primers/ snapshot the pipeline actually reads."""
import json

from mfe_common import (
    PRIMER_MANIFEST_COLUMNS, allocate_ids, load_primers_manifest, load_primers_manifest_raw,
    sync_input_primers, validate_primer_manifest_rows, write_groups, write_primers_manifest,
)


def _row(**kw):
    base = {"primer_id": "", "primer_name": "SET1", "segment": "S", "forward_primer": "ACGTRY",
            "reverse_primer": "TTGGCC", "probe": "GGCCAA", "fluorophore": "FAM", "source": "VendorD",
            "proceed": "True", "add_date": "2026-01-01"}
    base.update(kw)
    return base


def test_ids_are_assigned_once_and_never_reused(tmp_path):
    p = tmp_path / "primers.tsv"
    rows = write_primers_manifest(p, [_row(), _row(primer_name="SET2")])
    assert [r["primer_id"] for r in rows] == ["P001", "P002"]
    # delete P002, add a new one: the new row must not take P002 back
    rows = write_primers_manifest(p, [rows[0], _row(primer_name="SET3")])
    assert [r["primer_id"] for r in rows] == ["P001", "P003"]
    assert json.loads((tmp_path / "id_registry.json").read_text())["primer_next"] == 4


def test_rename_keeps_id_and_hand_typed_ids_bump_the_counter(tmp_path):
    p = tmp_path / "primers.tsv"
    rows = write_primers_manifest(p, [_row(primer_id="P010")])
    rows[0]["primer_name"] = "RENAMED"
    rows = write_primers_manifest(p, rows + [_row(primer_name="NEW")])
    assert [r["primer_id"] for r in rows] == ["P010", "P011"]
    assert allocate_ids(p, "primer", [], 1) == ["P012"]


def test_proceed_false_rows_are_dropped_and_ids_derived(tmp_path):
    p = tmp_path / "primers.tsv"
    write_primers_manifest(p, [_row(), _row(primer_name="SET2", proceed="False")])
    rows = load_primers_manifest(p)
    assert [r["set_id"] for r in rows] == ["P001"]
    r = rows[0]
    assert r["name"] == "SET1"
    assert (r["F_id"], r["R_id"], r["P_id"]) == ("P001_F", "P001_R", "P001_P")
    assert (r["F_seq"], r["R_seq"], r["P_seq"]) == ("ACGTRY", "TTGGCC", "GGCCAA")
    assert len(load_primers_manifest_raw(p)) == 2  # the spreadsheet still sees the unchecked row


def test_add_date_is_stamped_only_when_missing(tmp_path):
    p = tmp_path / "primers.tsv"
    out = write_primers_manifest(p, [_row(add_date=""), _row(primer_name="OLD", add_date="2020-05-05")])
    assert out[0]["add_date"] != ""
    assert out[1]["add_date"] == "2020-05-05"


def test_sync_input_primers_only_touches_changed_primers(tmp_path):
    manifest, groups, inp = tmp_path / "primers.tsv", tmp_path / "set_groups.tsv", tmp_path / "input_primers"
    write_primers_manifest(manifest, [_row(), _row(primer_name="SET2", forward_primer="AAAAAA"),
                                      _row(primer_name="OFF", proceed="False")])
    write_groups(groups, [{"group_id": "", "label": "ALL", "segment": "S", "member_ids": "P001;P002;P003"}])
    assert sync_input_primers(manifest, groups, inp) == ["P001", "P002"]
    assert sorted(f.name for f in inp.glob("*.tsv")) == ["P001.tsv", "P002.tsv", "groups.tsv"]
    # the snapshot carries sequences only -- no name, so a pure rename can't trigger recomputation
    assert "SET1" not in (inp / "P001.tsv").read_text()
    # inactive member is dropped from the group snapshot
    assert (inp / "groups.tsv").read_text().splitlines()[1] == "G001\tALL\tS\tP001;P002"

    mtimes = {f.name: f.stat().st_mtime_ns for f in inp.glob("*.tsv")}
    rows = load_primers_manifest_raw(manifest)
    rows[0]["primer_name"] = "RENAMED"
    rows[0]["fluorophore"] = "HEX"
    rows[1]["reverse_primer"] = "CCCCCC"
    write_primers_manifest(manifest, rows)
    sync_input_primers(manifest, groups, inp)
    after = {f.name: f.stat().st_mtime_ns for f in inp.glob("*.tsv")}
    assert after["P001.tsv"] == mtimes["P001.tsv"]   # rename/fluorophore: untouched
    assert after["P002.tsv"] != mtimes["P002.tsv"]   # sequence edit: rewritten


def test_group_without_active_members_is_dropped(tmp_path):
    manifest, groups, inp = tmp_path / "primers.tsv", tmp_path / "set_groups.tsv", tmp_path / "input_primers"
    write_primers_manifest(manifest, [_row(proceed="False")])
    write_groups(groups, [{"group_id": "", "label": "ALL", "segment": "S", "member_ids": "P001"}])
    sync_input_primers(manifest, groups, inp)
    assert (inp / "groups.tsv").read_text().strip().count("\n") == 0  # header only


def test_blank_row_is_not_an_error_but_bad_rows_are():
    assert validate_primer_manifest_rows([_row(), {c: "" for c in PRIMER_MANIFEST_COLUMNS}]) == []
    errs = validate_primer_manifest_rows([
        _row(), _row(),                       # duplicate name
        _row(primer_name="X", segment="Q"),   # bad segment
        _row(primer_name="Y", forward_primer="ACGT!"),  # bad character
        _row(primer_name="Z", probe=""),      # missing probe
        _row(primer_name="W", proceed="maybe"),
    ])
    fields = {(e["row"], e["field"]) for e in errs}
    assert {(1, "primer_name"), (2, "segment"), (3, "forward_primer"), (4, "probe"), (5, "proceed")} <= fields


def test_remark_is_kept_in_manifest_but_not_in_the_input_snapshot(tmp_path):
    manifest, groups, inp = tmp_path / "primers.tsv", tmp_path / "set_groups.tsv", tmp_path / "input_primers"
    write_primers_manifest(manifest, [_row(remark="VendorD oligo #4471")])
    assert load_primers_manifest_raw(manifest)[0]["remark"] == "VendorD oligo #4471"
    sync_input_primers(manifest, groups, inp)
    assert "4471" not in (inp / "P001.tsv").read_text()  # editing a remark must never trigger recomputation
    # a manifest written before the column existed still loads
    old = tmp_path / "old.tsv"
    old.write_text("primer_id\tprimer_name\tsegment\tforward_primer\treverse_primer\tprobe\tproceed\nP001\tA\tS\tAC\tGT\tAA\tTrue\n")
    assert load_primers_manifest(old)[0]["set_id"] == "P001"
