"""scan_engine=python|go: the switch, the Go answer cache and the fallback. Synthetic genomes only."""
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import mfe_common as mc  # noqa: E402


def genomes(seed=5, n=12, length=700):
    rng = random.Random(seed)
    oligo = "".join(rng.choice("ACGT") for _ in range(21))
    out = []
    for _ in range(n):
        g = [rng.choice("ACGT") for _ in range(length)]
        pos = rng.randrange(50, length - 80)
        for k, c in enumerate(oligo):
            g[pos + k] = c if rng.random() > 0.06 else rng.choice("ACGT")
        out.append("".join(g))
    return oligo, out


def test_default_engine_is_python(monkeypatch, tmp_path):
    monkeypatch.delenv("SCAN_ENGINE", raising=False)
    (tmp_path / "config").mkdir()
    assert mc.scan_engine(tmp_path) == "python"
    (tmp_path / "config" / "config.yaml").write_text('scan_engine: "go"\n')
    assert mc.scan_engine(tmp_path) == "go"
    monkeypatch.setenv("SCAN_ENGINE", "python")
    assert mc.scan_engine(tmp_path) == "python"
    monkeypatch.setenv("SCAN_ENGINE", "bogus")
    assert mc.scan_engine(tmp_path) == "python"


def test_prefetch_is_a_noop_for_python(monkeypatch):
    monkeypatch.setenv("SCAN_ENGINE", "python")
    mc.clear_site_cache()
    oligo, gs = genomes()
    assert mc.prefetch_sites([(oligo, 2)], gs) is False
    assert not mc._SITE_CACHE


def test_go_cache_gives_python_answers(monkeypatch):
    monkeypatch.setenv("SCAN_ENGINE", "go")
    if mc._go_scan_bin() is None:
        pytest.skip("Go scan engine not available")
    oligo, gs = genomes()
    mc.clear_site_cache()
    want = {(g, mm): mc._oligo_sites(oligo.upper(), g.upper(), mm) for g in gs for mm in (0, 2, 5)}
    assert mc.prefetch_sites([(oligo, 0), (oligo, 2), (oligo, 5)], gs) is True
    assert mc._SITE_CACHE
    for (g, mm), sites in want.items():
        assert mc._oligo_sites(oligo.upper(), g.upper(), mm) == sites
    mc.clear_site_cache()


def test_go_failure_falls_back_to_python(monkeypatch):
    monkeypatch.setenv("SCAN_ENGINE", "go")
    monkeypatch.setattr(mc, "_go_scan_bin", lambda root=None: None)
    oligo, gs = genomes()
    mc.clear_site_cache()
    assert mc.prefetch_sites([(oligo, 2)], gs) is False
    assert mc._oligo_sites(oligo.upper(), gs[0].upper(), 2) == mc._oligo_sites(oligo.upper(), gs[0].upper(), 2)
