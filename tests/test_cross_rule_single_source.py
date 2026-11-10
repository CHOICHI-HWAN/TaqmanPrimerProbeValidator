"""The cross-binding pass rule must come from ONE place (group_cross_ids, current primer/probe limits).

Group detail (lineage card), summary, overview and the Venn page used to disagree on a project whose primer_mismatch_pass is not 1,
because the group detail carried its own hard-coded "F/R <= 1, probe == 0" copy.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _body(src, name):
    m = re.search(rf"^def {name}\(.*?(?=^def |\Z)", src, re.S | re.M)
    assert m, name
    return m.group(0)


def test_group_detail_uses_shared_cross_rule():
    src = (ROOT / "webapp" / "server.py").read_text(encoding="utf-8")
    body = _body(src, "load_group_detail")
    assert "group_cross_ids(group_id)" in body
    assert not re.search(r"f_mismatch\"\]\s*<=\s*\d", body), "hard-coded mismatch limit in load_group_detail"


def test_ui_pass_text_follows_params():
    js = (ROOT / "webapp" / "static" / "app.js").read_text(encoding="utf-8")
    assert "mismatch&le;${PASS.mm}" in js
