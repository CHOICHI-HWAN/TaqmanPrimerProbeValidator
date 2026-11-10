"""Every local image / file the READMEs point to must exist (a broken picture on GitHub is easy to miss)."""
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.mark.parametrize("name", ["README.md", "README.ko.md", "docs/PIPELINE.md"])
def test_local_targets_exist(name):
    path = ROOT / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    text = path.read_text(encoding="utf-8")
    targets = set(re.findall(r"\]\((?!https?:|#|mailto:)([^)\s#]+)", text)) | set(re.findall(r'src="(?!https?:)([^"]+)"', text))
    base = path.parent
    missing = [t for t in targets if not (base / t).exists()]
    assert not missing, f"{name}: missing {missing}"
