"""Runs every tests/js/*.test.mjs with node (skipped where node is not installed)."""
import shutil
import subprocess
from pathlib import Path

import pytest

JS_DIR = Path(__file__).resolve().parent / "js"


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_msa_pi_bars_use_reference_positions():
    r = subprocess.run(["node", "--test", *map(str, sorted(JS_DIR.glob("*.test.mjs")))], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stdout[-1500:] + r.stderr[-500:]
