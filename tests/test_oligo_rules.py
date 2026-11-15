"""oligo_rules.py: the hard rule (probe 5' G with a FAM reporter) and the remark warnings. Synthetic short oligos only."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import oligo_rules as o  # noqa: E402


def row(f="ACGTACGTAC", r="TGCATGCATG", p="CACGTACGTA", fl="FAM"):
    return {"forward_primer": f, "reverse_primer": r, "probe": p, "fluorophore": fl}


def test_probe_five_prime_g_with_fam_is_hard():
    assert o.check(row(p="GACGTACGTA"))[0] == ["H1"]
    assert o.check(row(p="GACGTACGTA", fl="CY5"))[0] == []
    assert o.check(row())[0] == []


def test_warnings():
    assert any("G가 C보다" in w for w in o.check(row(p="CGGTACGTAG"))[1])
    assert any("같은 염기 4개" in w for w in o.check(row(f="ACGGGGTACA"))[1])
    assert any("혼합 염기" in w for w in o.check(row(f="ACGTACGTAR"))[1])


def test_remark_marker_is_replaced_not_stacked():
    a = o.with_warnings("note", ["x"])
    b = o.with_warnings(a, ["y"])
    assert b.startswith("note") and "x" not in b and "y" in b
    assert o.with_warnings(b, []) == "note"
