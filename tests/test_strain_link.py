"""Draft S/L link by strain name (webapp/strain_link.py). Synthetic metadata only."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
import strain_link as sl  # noqa: E402

COLS = ["accession", "strain", "country"]


def write(root, seg, rows):
    p = Path(root) / "data" / "raw"
    p.mkdir(parents=True, exist_ok=True)
    (p / f"{seg}_metadata.tsv").write_text("\n".join(["\t".join(COLS)] + ["\t".join(r) for r in rows]) + "\n")


def make(tmp_path):
    write(tmp_path, "s", [("S1", "Josiah", "Sierra Leone"), ("S2", "NIG08", "Nigeria"), ("S3", "dup", "Mali"), ("S4", "dup", "Mali"),
                          ("S5", "onlyS", "Guinea"), ("S6", "Clash", "Liberia"), ("S7", "padded  name", "Togo")])
    write(tmp_path, "l", [("L1", "JOSIAH", "Sierra Leone"), ("L2", "NIG08", "Nigeria"), ("L3", "dup", "Mali"), ("L6", "Clash", "Nigeria"),
                          ("L8", "onlyL", "Benin"), ("L7", "padded name", "Togo")])
    ids = {"s": {"S1", "S2", "S3", "S4", "S5", "S6", "S7"}, "l": {"L1", "L2", "L3", "L6", "L8", "L7"}}
    lin = {"s": {"S1": "LIV_V", "S2": "LII", "S3": "LII", "S4": "LII", "S5": "LII", "S6": "LIII", "S7": "LVII"},
           "l": {"L1": "LIV_V", "L2": "LII", "L3": "LII", "L6": "LII", "L8": "LIII", "L7": "LVII"}}
    return sl.build(tmp_path, ["s", "l"], ids, lin)


def test_statuses_and_counts(tmp_path):
    link = make(tmp_path)
    by = {v["virus_id"]: v for v in link["viruses"]}
    assert by["josiah"]["status"] == "linked"                # case does not matter
    assert by["padded name"]["status"] == "linked"           # runs of spaces do not matter
    assert by["dup"]["status"] == "ambiguous" and by["dup"]["acc"]["s"] == ["S3", "S4"]
    assert by["onlys"]["status"] == "s_only" and by["onlyl"]["status"] == "l_only"
    c = link["counts"]
    assert (c["linked"], c["ambiguous"], c["s_only"], c["l_only"]) == (4, 1, 1, 1)
    assert c["records"] == {"s": 7, "l": 6}


def test_conflicts_are_listed_not_hidden(tmp_path):
    link = make(tmp_path)
    assert [x["virus_id"] for x in link["conflicts"]] == ["clash"]       # different country AND lineage
    clash = next(v for v in link["viruses"] if v["virus_id"] == "clash")
    assert clash["evidence"] == {"country_agrees": False, "lineage_agrees": False}


def test_universe_and_id_mapping(tmp_path):
    link = make(tmp_path)
    uni = {v["virus_id"] for v in sl.universe(link)}
    assert uni == {"josiah", "nig08", "dup", "clash", "padded name"}     # a record in every segment; the single-segment ones are not judged
    assert sl.map_ids(link, "s", {"S1", "S4", "S5"}) == ["dup", "josiah"]   # S5 is S-only: outside the universe; any record of a virus counts
    assert sl.map_ids(link, "l", {"L2"}) == ["nig08"]


def test_tsv_has_one_row_per_virus(tmp_path):
    link = make(tmp_path)
    lines = sl.to_tsv(link).strip().split("\n")
    assert lines[0].split("\t")[:3] == ["virus_id", "strain", "status"] and len(lines) == len(link["viruses"]) + 1
