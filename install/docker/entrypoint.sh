#!/usr/bin/env bash
# Container start: install / update the panels in /work with this image's environment, start them, then stay up.
set -e
mkdir -p "$HOME"
git config --global --add safe.directory '*' 2>/dev/null || true
bash /opt/chichipes/setup.sh --backend container --dir /work --yes "$@"
echo "chichipes UI running (container). Logs: /work/logs"
exec sleep infinity
