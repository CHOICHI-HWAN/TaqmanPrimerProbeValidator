"""The MSA card of a set page must have its export / snapshot / copy buttons the moment the page is drawn -- not only after the
first filter change re-renders the header (the SeaLion FASTA button was missing on first view)."""
import re
from pathlib import Path

JS = (Path(__file__).resolve().parent.parent / "webapp" / "static" / "app.js").read_text(encoding="utf-8")
BUTTONS = ("msaExportBtn", "msaSnapBtn", "msaCopyImgBtn", "msaCopyTxtBtn", "msaFullscreenBtn")


def test_first_render_header_has_all_msa_buttons():
    start = JS.index('<div class="card msa-card" style="margin-bottom:20px;">')
    header = JS[start:JS.index("renderMSA(setData)", start)]
    missing = [b for b in BUTTONS if f'id="{b}"' not in header]
    assert not missing, f"missing in the first-render MSA header: {missing}"


def test_rerendered_header_has_the_same_buttons():
    start = JS.index("const MSA_HEADER")
    header = JS[start:JS.index("renderMSA(setData)", start)]
    assert all(f'id="{b}"' in header for b in BUTTONS)
