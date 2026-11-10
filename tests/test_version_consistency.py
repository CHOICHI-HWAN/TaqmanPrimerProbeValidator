"""The version is written in several files by hand: they must all say what VERSION says (the NiV dashboard once showed v0.21.0 next to 0.22.x)."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _find(path, pattern):
    return re.findall(pattern, (ROOT / path).read_text(encoding="utf-8"), re.M)


def test_all_version_strings_agree():
    v = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    assert _find("CITATION.cff", r'^version:\s*"([^"]+)"') == [v]
    assert _find("webapp/package.json", r'"version":\s*"([^"]+)"') == [v]
    for readme in ("README.md", "README.ko.md"):
        assert _find(readme, r"\*\*v(\d+\.\d+\.\d+)") == [v], readme
        assert _find(readme, r"version = \{([^}]+)\}") == [v], readme
    assert _find("webapp/static/index.html", r'id="footerVersion">([^<]+)<') == [v]
    assert _find("CHANGELOG.md", r"^## \[(\d+\.\d+\.\d+)\]")[0] == v
