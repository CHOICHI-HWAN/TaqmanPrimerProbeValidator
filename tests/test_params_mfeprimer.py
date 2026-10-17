"""mfeprimer path + env_prefix in config.yaml: reading, validating on save, and the tool lookup that uses them."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))

from mfe_common import check_mfeprimer, env_tool, find_mfeprimer, read_config_str  # noqa: E402

CONFIG = 'primer_mismatch_pass: 1\nenv_prefix: ""\nmfeprimer_bin: ""   # path\nmfeprimer_cpu: 4\n'


def fake_mfeprimer(tmp_path, body="echo 'mfeprimer v9.9.9'", mode=0o755):
    p = tmp_path / "mfeprimer-fake"
    p.write_text(f"#!/bin/sh\n{body}\n")
    p.chmod(mode)
    return p


def project(tmp_path, text=CONFIG):
    (tmp_path / "config").mkdir(exist_ok=True)
    (tmp_path / "config" / "config.yaml").write_text(text)
    return tmp_path


def test_check_mfeprimer_accepts_a_runnable_binary_and_explains_failures(tmp_path):
    ok = check_mfeprimer(str(fake_mfeprimer(tmp_path)))
    assert ok["ok"] and ok["version"] == "mfeprimer v9.9.9"
    assert "찾을 수 없" in check_mfeprimer(str(tmp_path / "nope"))["error"]
    assert "비어" in check_mfeprimer("")["error"]
    assert "실행 권한" in check_mfeprimer(str(fake_mfeprimer(tmp_path, mode=0o644)))["error"]
    assert "아닌" in check_mfeprimer(str(fake_mfeprimer(tmp_path, body="echo hello")))["error"]
    assert not check_mfeprimer(str(fake_mfeprimer(tmp_path, body="exit 3")))["ok"]


def test_read_config_str_handles_quotes_comments_and_missing_keys(tmp_path):
    root = project(tmp_path, 'a: "/with space/x"  # why\nb: plain\nc: ""\n')
    assert read_config_str(root, "a") == "/with space/x"
    assert read_config_str(root, "b") == "plain" and read_config_str(root, "c") == "" and read_config_str(root, "zz") == ""


def test_env_tool_and_find_mfeprimer_follow_the_config(tmp_path):
    (tmp_path / "env" / "bin").mkdir(parents=True)
    (tmp_path / "env" / "bin" / "snakemake").write_text("")
    mfe = fake_mfeprimer(tmp_path)
    proj = tmp_path / "proj"
    proj.mkdir()
    root = project(proj, f'env_prefix: "{tmp_path / "env"}"\nmfeprimer_bin: "{mfe}"\n')
    assert env_tool(root, "snakemake") == str(tmp_path / "env" / "bin" / "snakemake") and env_tool(root, "mafft") == ""
    assert find_mfeprimer("", root) == str(mfe)


def test_save_params_validates_the_mfeprimer_path_and_keeps_the_rest(tmp_path, monkeypatch):
    import server
    root = project(tmp_path)
    monkeypatch.setattr(server, "ROOT", root)
    with pytest.raises(ValueError, match="mfeprimer 경로"):
        server.save_params({"mfeprimer_bin": str(tmp_path / "missing")})
    assert (root / "config" / "config.yaml").read_text() == CONFIG                 # a rejected value writes nothing
    good = fake_mfeprimer(tmp_path)
    server.save_params({"mfeprimer_bin": str(good), "primer_mismatch_pass": 2})
    text = (root / "config" / "config.yaml").read_text()
    assert f'mfeprimer_bin: "{good}"' in text and "primer_mismatch_pass: 2" in text and "env_prefix" in text
    server.save_params({"mfeprimer_bin": ""})                                      # empty = auto-discover again
    assert 'mfeprimer_bin: ""' in (root / "config" / "config.yaml").read_text()
    p = server.load_params()
    assert p["mfeprimer_bin"] == "" and set(p["tools"]) == {"env_prefix", "snakemake", "mafft"}
