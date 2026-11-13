"""The subtype part filter of the influenza NCBI fetch (HA part H3, NA part N2 ...)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from fetch_ncbi_virus import subtype_has  # noqa: E402


def test_ha_and_na_parts_of_a_subtype():
    assert subtype_has("H3N2", "H3") and subtype_has("H3N2", "N2")
    assert subtype_has("H5N1", "N1") and not subtype_has("H1N1", "N2")
    assert not subtype_has("H3N2", "H1") and not subtype_has("unclassified", "H1")
