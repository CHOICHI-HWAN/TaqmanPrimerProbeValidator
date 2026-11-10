"""CITATION.cff (GitHub's "Cite this repository") must be valid and agree with the project."""
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent


def load():
    return yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))


def test_required_fields():
    c = load()
    for key in ("cff-version", "message", "title", "authors", "version", "date-released", "license", "repository-code"):
        assert c.get(key), key
    assert str(c["cff-version"]).startswith("1.2")
    assert all(a.get("family-names") and a.get("given-names") for a in c["authors"])


def test_version_matches_the_project():
    assert str(load()["version"]) == (ROOT / "VERSION").read_text(encoding="utf-8").strip()


def test_the_readmes_show_the_same_reference():
    c = load()
    for name in ("README.md", "README.ko.md"):
        text = (ROOT / name).read_text(encoding="utf-8")
        assert "@software" in text and c["repository-code"] in text and f"version = {{{c['version']}}}" in text, name
        assert "CITATION.cff" in text
