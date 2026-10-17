"""The dashboard holds confidential sequences and has write endpoints without a login: it must listen on this computer only
unless --host is given explicitly."""
import re
from pathlib import Path

SRC = (Path(__file__).resolve().parent.parent / "webapp" / "server.py").read_text(encoding="utf-8")


def test_default_listen_address_is_loopback():
    m = re.search(r'add_argument\("--host",\s*default="([^"]+)"', SRC)
    assert m and m.group(1) == "127.0.0.1"
